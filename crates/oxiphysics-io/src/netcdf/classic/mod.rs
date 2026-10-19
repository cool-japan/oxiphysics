// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! NetCDF classic binary format (`CDF-1`, `CDF-2` and `CDF-5`).
//!
//! This module reads and writes the on-disk encoding defined by the
//! [NetCDF classic format specification], which is what AMBER, cpptraj, VMD,
//! MDAnalysis, `ncdump` and every netCDF-3 library expect:
//!
//! * [`CdfVersion::Classic`] — magic `CDF\x01`, 32-bit `begin` offsets;
//! * [`CdfVersion::Offset64`] — magic `CDF\x02`, 64-bit `begin` offsets (the
//!   variant required by the AMBER trajectory convention);
//! * [`CdfVersion::Data64`] — magic `CDF\x05` (CDF-5 / PnetCDF), 64-bit counts
//!   and offsets plus the unsigned and 64-bit integer types.
//!
//! The lossless data model is [`CdfDataset`]: typed dimensions, attributes
//! and variables.  [`write_cdf_dataset`] / [`read_cdf_dataset`] serialise it.
//! [`write_cdf`] / [`read_cdf`] bridge the simpler, `f64`-valued
//! [`NetCdfFile`](super::NetCdfFile) model.
//!
//! [NetCDF classic format specification]: https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html

mod bridge;
mod reader;
#[cfg(test)]
mod tests;
mod writer;

pub use bridge::{read_cdf, read_cdf_file, write_cdf, write_cdf_file};
pub use reader::read_cdf_dataset;
pub use writer::write_cdf_dataset;

use crate::error::{Error, Result};

/// Header tag introducing the dimension list.
pub(crate) const NC_DIMENSION: u32 = 0x0A;
/// Header tag introducing a variable list.
pub(crate) const NC_VARIABLE: u32 = 0x0B;
/// Header tag introducing an attribute list.
pub(crate) const NC_ATTRIBUTE: u32 = 0x0C;

/// On-disk variant of the classic format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CdfVersion {
    /// Classic format (`CDF\x01`): 32-bit offsets, files below 2 GiB of
    /// fixed-size data.
    Classic,
    /// 64-bit offset format (`CDF\x02`).
    Offset64,
    /// 64-bit data format, CDF-5 (`CDF\x05`): 64-bit counts and offsets and
    /// the extended integer types.
    Data64,
}

impl CdfVersion {
    /// The fourth magic byte (`1`, `2` or `5`).
    pub fn version_byte(self) -> u8 {
        match self {
            CdfVersion::Classic => 1,
            CdfVersion::Offset64 => 2,
            CdfVersion::Data64 => 5,
        }
    }

    /// Variant identified by the fourth magic byte.
    pub fn from_version_byte(byte: u8) -> Option<Self> {
        match byte {
            1 => Some(CdfVersion::Classic),
            2 => Some(CdfVersion::Offset64),
            5 => Some(CdfVersion::Data64),
            _ => None,
        }
    }

    /// The four magic bytes that start a file of this variant.
    pub fn magic(self) -> [u8; 4] {
        [b'C', b'D', b'F', self.version_byte()]
    }

    /// Whether variables and attributes of type `nc_type` can be stored.
    pub fn supports(self, nc_type: NcType) -> bool {
        !nc_type.is_extended() || self == CdfVersion::Data64
    }

    /// Width in bytes of the `begin` offset field.
    pub(crate) fn offset_width(self) -> usize {
        match self {
            CdfVersion::Classic => 4,
            CdfVersion::Offset64 | CdfVersion::Data64 => 8,
        }
    }

    /// Width in bytes of the `NON_NEG` count fields (lengths, `nelems`,
    /// `numrecs`, dimension ids, `vsize`).
    pub(crate) fn count_width(self) -> usize {
        match self {
            CdfVersion::Classic | CdfVersion::Offset64 => 4,
            CdfVersion::Data64 => 8,
        }
    }

    /// Largest `begin` offset representable (`OFFSET` is a signed integer).
    pub(crate) fn max_offset(self) -> u64 {
        match self {
            CdfVersion::Classic => i32::MAX as u64,
            CdfVersion::Offset64 | CdfVersion::Data64 => i64::MAX as u64,
        }
    }

    /// Largest count written to a `NON_NEG` field.  CDF-1/CDF-2 fields are
    /// read as unsigned 32-bit values by netCDF-C; the all-ones pattern is
    /// reserved (streaming `numrecs`, oversized `vsize`).
    pub(crate) fn max_count(self) -> u64 {
        match self {
            CdfVersion::Classic | CdfVersion::Offset64 => u64::from(u32::MAX - 1),
            CdfVersion::Data64 => i64::MAX as u64,
        }
    }
}

/// External data type of a NetCDF variable or attribute (`nc_type`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum NcType {
    /// `NC_BYTE`: 8-bit signed integer.
    Byte,
    /// `NC_CHAR`: 8-bit character (text).
    Char,
    /// `NC_SHORT`: 16-bit signed integer.
    Short,
    /// `NC_INT`: 32-bit signed integer.
    Int,
    /// `NC_FLOAT`: 32-bit IEEE 754 floating point.
    Float,
    /// `NC_DOUBLE`: 64-bit IEEE 754 floating point.
    #[default]
    Double,
    /// `NC_UBYTE`: 8-bit unsigned integer (CDF-5 only).
    UByte,
    /// `NC_USHORT`: 16-bit unsigned integer (CDF-5 only).
    UShort,
    /// `NC_UINT`: 32-bit unsigned integer (CDF-5 only).
    UInt,
    /// `NC_INT64`: 64-bit signed integer (CDF-5 only).
    Int64,
    /// `NC_UINT64`: 64-bit unsigned integer (CDF-5 only).
    UInt64,
}

impl NcType {
    /// Every type, in `nc_type` code order.
    pub const ALL: [NcType; 11] = [
        NcType::Byte,
        NcType::Char,
        NcType::Short,
        NcType::Int,
        NcType::Float,
        NcType::Double,
        NcType::UByte,
        NcType::UShort,
        NcType::UInt,
        NcType::Int64,
        NcType::UInt64,
    ];

    /// The `nc_type` code stored in the file header.
    pub fn code(self) -> u32 {
        match self {
            NcType::Byte => 1,
            NcType::Char => 2,
            NcType::Short => 3,
            NcType::Int => 4,
            NcType::Float => 5,
            NcType::Double => 6,
            NcType::UByte => 7,
            NcType::UShort => 8,
            NcType::UInt => 9,
            NcType::Int64 => 10,
            NcType::UInt64 => 11,
        }
    }

    /// The type identified by a header `nc_type` code.
    pub fn from_code(code: u32) -> Option<Self> {
        NcType::ALL.iter().copied().find(|t| t.code() == code)
    }

    /// Size in bytes of one external value.
    pub fn size(self) -> usize {
        match self {
            NcType::Byte | NcType::Char | NcType::UByte => 1,
            NcType::Short | NcType::UShort => 2,
            NcType::Int | NcType::Float | NcType::UInt => 4,
            NcType::Double | NcType::Int64 | NcType::UInt64 => 8,
        }
    }

    /// The CDL keyword (`ncdump` spelling) of this type.
    pub fn cdl_name(self) -> &'static str {
        match self {
            NcType::Byte => "byte",
            NcType::Char => "char",
            NcType::Short => "short",
            NcType::Int => "int",
            NcType::Float => "float",
            NcType::Double => "double",
            NcType::UByte => "ubyte",
            NcType::UShort => "ushort",
            NcType::UInt => "uint",
            NcType::Int64 => "int64",
            NcType::UInt64 => "uint64",
        }
    }

    /// The type named by a CDL keyword (`long` and `real` are accepted as the
    /// legacy synonyms of `int` and `float`).
    pub fn from_cdl_name(name: &str) -> Option<Self> {
        match name {
            "long" => Some(NcType::Int),
            "real" => Some(NcType::Float),
            other => NcType::ALL.iter().copied().find(|t| t.cdl_name() == other),
        }
    }

    /// Whether the type exists only in the CDF-5 format.
    pub fn is_extended(self) -> bool {
        matches!(
            self,
            NcType::UByte | NcType::UShort | NcType::UInt | NcType::Int64 | NcType::UInt64
        )
    }
}

/// A typed array of values, as stored in a variable or an attribute.
#[derive(Debug, Clone, PartialEq)]
pub enum NcValues {
    /// `NC_BYTE` values.
    Byte(Vec<i8>),
    /// `NC_CHAR` values (raw bytes; text is usually UTF-8 or ASCII).
    Char(Vec<u8>),
    /// `NC_SHORT` values.
    Short(Vec<i16>),
    /// `NC_INT` values.
    Int(Vec<i32>),
    /// `NC_FLOAT` values.
    Float(Vec<f32>),
    /// `NC_DOUBLE` values.
    Double(Vec<f64>),
    /// `NC_UBYTE` values.
    UByte(Vec<u8>),
    /// `NC_USHORT` values.
    UShort(Vec<u16>),
    /// `NC_UINT` values.
    UInt(Vec<u32>),
    /// `NC_INT64` values.
    Int64(Vec<i64>),
    /// `NC_UINT64` values.
    UInt64(Vec<u64>),
}

impl NcValues {
    /// `NC_CHAR` values holding the bytes of `text`.
    pub fn text(text: &str) -> Self {
        NcValues::Char(text.as_bytes().to_vec())
    }

    /// The external type of the values.
    pub fn nc_type(&self) -> NcType {
        match self {
            NcValues::Byte(_) => NcType::Byte,
            NcValues::Char(_) => NcType::Char,
            NcValues::Short(_) => NcType::Short,
            NcValues::Int(_) => NcType::Int,
            NcValues::Float(_) => NcType::Float,
            NcValues::Double(_) => NcType::Double,
            NcValues::UByte(_) => NcType::UByte,
            NcValues::UShort(_) => NcType::UShort,
            NcValues::UInt(_) => NcType::UInt,
            NcValues::Int64(_) => NcType::Int64,
            NcValues::UInt64(_) => NcType::UInt64,
        }
    }

    /// Number of values.
    pub fn len(&self) -> usize {
        match self {
            NcValues::Byte(v) => v.len(),
            NcValues::Char(v) | NcValues::UByte(v) => v.len(),
            NcValues::Short(v) => v.len(),
            NcValues::Int(v) => v.len(),
            NcValues::Float(v) => v.len(),
            NcValues::Double(v) => v.len(),
            NcValues::UShort(v) => v.len(),
            NcValues::UInt(v) => v.len(),
            NcValues::Int64(v) => v.len(),
            NcValues::UInt64(v) => v.len(),
        }
    }

    /// Whether there are no values.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// An empty array of type `nc_type`.
    pub fn empty(nc_type: NcType) -> Self {
        match nc_type {
            NcType::Byte => NcValues::Byte(Vec::new()),
            NcType::Char => NcValues::Char(Vec::new()),
            NcType::Short => NcValues::Short(Vec::new()),
            NcType::Int => NcValues::Int(Vec::new()),
            NcType::Float => NcValues::Float(Vec::new()),
            NcType::Double => NcValues::Double(Vec::new()),
            NcType::UByte => NcValues::UByte(Vec::new()),
            NcType::UShort => NcValues::UShort(Vec::new()),
            NcType::UInt => NcValues::UInt(Vec::new()),
            NcType::Int64 => NcValues::Int64(Vec::new()),
            NcType::UInt64 => NcValues::UInt64(Vec::new()),
        }
    }

    /// The values as text, when they are `NC_CHAR`.
    ///
    /// Trailing NUL padding is removed and invalid UTF-8 is replaced by
    /// `U+FFFD`.
    pub fn as_text(&self) -> Option<String> {
        match self {
            NcValues::Char(bytes) => {
                let end = bytes.iter().rposition(|&b| b != 0).map_or(0, |p| p + 1);
                Some(String::from_utf8_lossy(&bytes[..end]).into_owned())
            }
            _ => None,
        }
    }

    /// Every value widened to `f64`.
    ///
    /// The conversion is exact for every type except `NC_INT64` / `NC_UINT64`
    /// magnitudes above 2^53, which round to the nearest `f64`.  `NC_CHAR`
    /// values become their byte codes.
    pub fn to_f64_vec(&self) -> Vec<f64> {
        match self {
            NcValues::Byte(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::Char(v) | NcValues::UByte(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::Short(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::Int(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::Float(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::Double(v) => v.clone(),
            NcValues::UShort(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::UInt(v) => v.iter().map(|&x| f64::from(x)).collect(),
            NcValues::Int64(v) => v.iter().map(|&x| x as f64).collect(),
            NcValues::UInt64(v) => v.iter().map(|&x| x as f64).collect(),
        }
    }

    /// Convert `f64` values to type `nc_type`.
    ///
    /// `NC_FLOAT` rounds to the nearest `f32` (values beyond the `f32` range
    /// are rejected, infinities and NaN are kept).  Integer and `NC_CHAR`
    /// targets require every value to be a finite integer inside the target
    /// range; anything else is an error rather than a silent truncation.
    pub fn from_f64(nc_type: NcType, values: &[f64]) -> Result<Self> {
        fn ints<T>(
            values: &[f64],
            lo: f64,
            hi: f64,
            ty: NcType,
            f: fn(f64) -> T,
        ) -> Result<Vec<T>> {
            values
                .iter()
                .enumerate()
                .map(|(i, &v)| {
                    if v.is_finite() && v.fract() == 0.0 && v >= lo && v <= hi {
                        Ok(f(v))
                    } else {
                        Err(Error::General(format!(
                            "value {v} at index {i} is not representable as NetCDF {}",
                            ty.cdl_name()
                        )))
                    }
                })
                .collect()
        }
        Ok(match nc_type {
            NcType::Byte => NcValues::Byte(ints(values, -128.0, 127.0, nc_type, |v| v as i8)?),
            NcType::Char => NcValues::Char(ints(values, 0.0, 255.0, nc_type, |v| v as u8)?),
            NcType::Short => NcValues::Short(ints(
                values,
                f64::from(i16::MIN),
                f64::from(i16::MAX),
                nc_type,
                |v| v as i16,
            )?),
            NcType::Int => NcValues::Int(ints(
                values,
                f64::from(i32::MIN),
                f64::from(i32::MAX),
                nc_type,
                |v| v as i32,
            )?),
            NcType::Float => NcValues::Float(
                values
                    .iter()
                    .enumerate()
                    .map(|(i, &v)| {
                        if v.is_finite() && v.abs() > f64::from(f32::MAX) {
                            Err(Error::General(format!(
                                "value {v} at index {i} exceeds the NetCDF float range"
                            )))
                        } else {
                            Ok(v as f32)
                        }
                    })
                    .collect::<Result<Vec<f32>>>()?,
            ),
            NcType::Double => NcValues::Double(values.to_vec()),
            NcType::UByte => NcValues::UByte(ints(values, 0.0, 255.0, nc_type, |v| v as u8)?),
            NcType::UShort => {
                NcValues::UShort(ints(values, 0.0, f64::from(u16::MAX), nc_type, |v| {
                    v as u16
                })?)
            }
            NcType::UInt => NcValues::UInt(ints(values, 0.0, f64::from(u32::MAX), nc_type, |v| {
                v as u32
            })?),
            // 2^63 and 2^64 are exactly representable; the bounds are
            // exclusive at the top so that the casts cannot saturate.
            NcType::Int64 => NcValues::Int64(ints(
                values,
                -9_223_372_036_854_775_808.0,
                9_223_372_036_854_774_784.0,
                nc_type,
                |v| v as i64,
            )?),
            NcType::UInt64 => NcValues::UInt64(ints(
                values,
                0.0,
                18_446_744_073_709_549_568.0,
                nc_type,
                |v| v as u64,
            )?),
        })
    }

    /// Render the values the way `ncdump` lists them (`", "`-separated
    /// numbers, or the text for `NC_CHAR`).
    pub fn to_display_string(&self) -> String {
        fn join<T: std::fmt::Display>(v: &[T]) -> String {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        }
        match self {
            NcValues::Char(_) => self.as_text().unwrap_or_default(),
            NcValues::Byte(v) => join(v),
            NcValues::Short(v) => join(v),
            NcValues::Int(v) => join(v),
            NcValues::Float(v) => join(v),
            NcValues::Double(v) => join(v),
            NcValues::UByte(v) => join(v),
            NcValues::UShort(v) => join(v),
            NcValues::UInt(v) => join(v),
            NcValues::Int64(v) => join(v),
            NcValues::UInt64(v) => join(v),
        }
    }
}

/// A dimension of a [`CdfDataset`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CdfDimension {
    /// Dimension name.
    pub name: String,
    /// Length.  For the unlimited (record) dimension this is the number of
    /// records currently stored (`numrecs`).
    pub len: usize,
    /// Whether this is the unlimited (record) dimension.
    pub unlimited: bool,
}

/// A named, typed attribute.
#[derive(Debug, Clone, PartialEq)]
pub struct CdfAttribute {
    /// Attribute name.
    pub name: String,
    /// Attribute values.
    pub values: NcValues,
}

impl CdfAttribute {
    /// A typed attribute.
    pub fn new(name: &str, values: NcValues) -> Self {
        CdfAttribute {
            name: name.to_string(),
            values,
        }
    }

    /// An `NC_CHAR` (text) attribute.
    pub fn text(name: &str, text: &str) -> Self {
        CdfAttribute::new(name, NcValues::text(text))
    }
}

/// A variable of a [`CdfDataset`].
#[derive(Debug, Clone, PartialEq)]
pub struct CdfVariable {
    /// Variable name.
    pub name: String,
    /// Indices into [`CdfDataset::dimensions`], slowest-varying first.
    pub dim_ids: Vec<usize>,
    /// Variable attributes.
    pub attributes: Vec<CdfAttribute>,
    /// Row-major values (all records for a record variable).
    pub values: NcValues,
}

impl CdfVariable {
    /// The attribute called `name`.
    pub fn attribute(&self, name: &str) -> Option<&CdfAttribute> {
        self.attributes.iter().find(|a| a.name == name)
    }

    /// The text of the `NC_CHAR` attribute called `name`.
    pub fn text_attribute(&self, name: &str) -> Option<String> {
        self.attribute(name).and_then(|a| a.values.as_text())
    }

    /// Add (or replace) an attribute.
    pub fn set_attribute(&mut self, attribute: CdfAttribute) {
        set_attribute_in(&mut self.attributes, attribute);
    }

    /// External type of the variable.
    pub fn nc_type(&self) -> NcType {
        self.values.nc_type()
    }
}

fn set_attribute_in(list: &mut Vec<CdfAttribute>, attribute: CdfAttribute) {
    if let Some(existing) = list.iter_mut().find(|a| a.name == attribute.name) {
        *existing = attribute;
    } else {
        list.push(attribute);
    }
}

/// A lossless in-memory NetCDF classic dataset.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct CdfDataset {
    /// Dimensions in definition order (their index is the dimension id).
    pub dimensions: Vec<CdfDimension>,
    /// Global attributes.
    pub attributes: Vec<CdfAttribute>,
    /// Variables in definition order.
    pub variables: Vec<CdfVariable>,
}

impl CdfDataset {
    /// An empty dataset.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a fixed-length dimension and return its id.
    ///
    /// # Errors
    ///
    /// Fails when the name is invalid or already used, or `len` is zero (a
    /// zero length denotes the unlimited dimension in the classic format).
    pub fn add_dimension(&mut self, name: &str, len: usize) -> Result<usize> {
        if len == 0 {
            return Err(Error::General(format!(
                "fixed dimension '{name}' must have a non-zero length"
            )));
        }
        self.push_dimension(name, len, false)
    }

    /// Add the unlimited (record) dimension holding `num_records` records and
    /// return its id.
    ///
    /// # Errors
    ///
    /// Fails when the name is invalid or already used, or the dataset
    /// already has an unlimited dimension.
    pub fn add_unlimited_dimension(&mut self, name: &str, num_records: usize) -> Result<usize> {
        if self.record_dimension().is_some() {
            return Err(Error::General(
                "a classic NetCDF dataset has at most one unlimited dimension".to_string(),
            ));
        }
        self.push_dimension(name, num_records, true)
    }

    fn push_dimension(&mut self, name: &str, len: usize, unlimited: bool) -> Result<usize> {
        validate_name(name)?;
        if self.dimension_id(name).is_some() {
            return Err(Error::General(format!("duplicate dimension '{name}'")));
        }
        self.dimensions.push(CdfDimension {
            name: name.to_string(),
            len,
            unlimited,
        });
        Ok(self.dimensions.len() - 1)
    }

    /// Add (or replace) a global attribute.
    pub fn set_attribute(&mut self, attribute: CdfAttribute) {
        set_attribute_in(&mut self.attributes, attribute);
    }

    /// Add (or replace) a global text attribute.
    pub fn set_text_attribute(&mut self, name: &str, text: &str) {
        self.set_attribute(CdfAttribute::text(name, text));
    }

    /// Add a variable over the named dimensions and return its index.
    ///
    /// # Errors
    ///
    /// Fails when the name is invalid or used, a dimension is unknown, the
    /// unlimited dimension is not the first one, or the number of values does
    /// not match the shape.
    pub fn add_variable(&mut self, name: &str, dims: &[&str], values: NcValues) -> Result<usize> {
        validate_name(name)?;
        if self.variable(name).is_some() {
            return Err(Error::General(format!("duplicate variable '{name}'")));
        }
        let dim_ids = dims
            .iter()
            .map(|d| {
                self.dimension_id(d).ok_or_else(|| {
                    Error::General(format!("variable '{name}' uses unknown dimension '{d}'"))
                })
            })
            .collect::<Result<Vec<usize>>>()?;
        let var = CdfVariable {
            name: name.to_string(),
            dim_ids,
            attributes: Vec::new(),
            values,
        };
        self.check_variable(&var)?;
        self.variables.push(var);
        Ok(self.variables.len() - 1)
    }

    /// Id of the dimension called `name`.
    pub fn dimension_id(&self, name: &str) -> Option<usize> {
        self.dimensions.iter().position(|d| d.name == name)
    }

    /// The dimension called `name`.
    pub fn dimension(&self, name: &str) -> Option<&CdfDimension> {
        self.dimensions.iter().find(|d| d.name == name)
    }

    /// Id of the unlimited dimension, if any.
    pub fn record_dimension(&self) -> Option<usize> {
        self.dimensions.iter().position(|d| d.unlimited)
    }

    /// Number of records (length of the unlimited dimension, 0 without one).
    pub fn num_records(&self) -> usize {
        self.record_dimension()
            .and_then(|id| self.dimensions.get(id))
            .map_or(0, |d| d.len)
    }

    /// The global attribute called `name`.
    pub fn attribute(&self, name: &str) -> Option<&CdfAttribute> {
        self.attributes.iter().find(|a| a.name == name)
    }

    /// The text of the global `NC_CHAR` attribute called `name`.
    pub fn text_attribute(&self, name: &str) -> Option<String> {
        self.attribute(name).and_then(|a| a.values.as_text())
    }

    /// The variable called `name`.
    pub fn variable(&self, name: &str) -> Option<&CdfVariable> {
        self.variables.iter().find(|v| v.name == name)
    }

    /// The variable called `name`, mutably.
    pub fn variable_mut(&mut self, name: &str) -> Option<&mut CdfVariable> {
        self.variables.iter_mut().find(|v| v.name == name)
    }

    /// Whether `var` is a record variable (its first dimension is unlimited).
    pub fn is_record_variable(&self, var: &CdfVariable) -> bool {
        var.dim_ids
            .first()
            .and_then(|&id| self.dimensions.get(id))
            .is_some_and(|d| d.unlimited)
    }

    /// Dimension lengths of `var` (the record dimension contributes the
    /// current number of records).
    pub fn variable_shape(&self, var: &CdfVariable) -> Vec<usize> {
        var.dim_ids
            .iter()
            .map(|&id| self.dimensions.get(id).map_or(0, |d| d.len))
            .collect()
    }

    /// Check every structural invariant the classic format relies on.
    ///
    /// # Errors
    ///
    /// Reports the first violation: invalid or duplicate names, more than one
    /// unlimited dimension, a zero-length fixed dimension, an unknown
    /// dimension id, the unlimited dimension in a non-leading position, or a
    /// value count that does not match the variable shape.
    pub fn validate(&self) -> Result<()> {
        let mut unlimited = 0usize;
        for (i, dim) in self.dimensions.iter().enumerate() {
            validate_name(&dim.name)?;
            if self.dimensions[..i].iter().any(|d| d.name == dim.name) {
                return Err(Error::General(format!(
                    "duplicate dimension '{}'",
                    dim.name
                )));
            }
            if dim.unlimited {
                unlimited += 1;
            } else if dim.len == 0 {
                return Err(Error::General(format!(
                    "fixed dimension '{}' must have a non-zero length",
                    dim.name
                )));
            }
        }
        if unlimited > 1 {
            return Err(Error::General(
                "a classic NetCDF dataset has at most one unlimited dimension".to_string(),
            ));
        }
        validate_attributes(&self.attributes, "global attributes")?;
        for (i, var) in self.variables.iter().enumerate() {
            validate_name(&var.name)?;
            if self.variables[..i].iter().any(|v| v.name == var.name) {
                return Err(Error::General(format!("duplicate variable '{}'", var.name)));
            }
            self.check_variable(var)?;
            validate_attributes(&var.attributes, &var.name)?;
        }
        Ok(())
    }

    fn check_variable(&self, var: &CdfVariable) -> Result<()> {
        let mut expected = 1usize;
        for (pos, &id) in var.dim_ids.iter().enumerate() {
            let dim = self.dimensions.get(id).ok_or_else(|| {
                Error::General(format!(
                    "variable '{}' refers to missing dimension id {id}",
                    var.name
                ))
            })?;
            if dim.unlimited && pos != 0 {
                return Err(Error::General(format!(
                    "variable '{}': the unlimited dimension '{}' must be its first dimension",
                    var.name, dim.name
                )));
            }
            expected = expected
                .checked_mul(dim.len)
                .ok_or_else(|| Error::General(format!("variable '{}' is too large", var.name)))?;
        }
        if var.values.len() != expected {
            return Err(Error::General(format!(
                "variable '{}' holds {} values but its shape {:?} needs {expected}",
                var.name,
                var.values.len(),
                self.variable_shape(var)
            )));
        }
        Ok(())
    }
}

fn validate_attributes(list: &[CdfAttribute], owner: &str) -> Result<()> {
    for (i, attr) in list.iter().enumerate() {
        validate_name(&attr.name)?;
        if list[..i].iter().any(|a| a.name == attr.name) {
            return Err(Error::General(format!(
                "duplicate attribute '{}' in {owner}",
                attr.name
            )));
        }
    }
    Ok(())
}

/// Check a dimension, variable or attribute name against the NetCDF rules:
/// non-empty, starting with an ASCII letter or digit, `_` or a non-ASCII
/// character, containing no `/` or control characters, and without trailing
/// whitespace.
pub fn validate_name(name: &str) -> Result<()> {
    let first_ok = name
        .chars()
        .next()
        .is_some_and(|c| c.is_ascii_alphanumeric() || c == '_' || !c.is_ascii());
    let body_ok = name.chars().all(|c| c != '/' && !c.is_control());
    let tail_ok = !name.ends_with(|c: char| c.is_whitespace());
    if first_ok && body_ok && tail_ok {
        Ok(())
    } else {
        Err(Error::General(format!("invalid NetCDF name '{name}'")))
    }
}
