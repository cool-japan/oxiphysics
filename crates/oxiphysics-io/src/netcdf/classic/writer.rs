// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Encoder for the NetCDF classic binary format.

use super::{
    CdfAttribute, CdfDataset, CdfVariable, CdfVersion, NC_ATTRIBUTE, NC_DIMENSION, NC_VARIABLE,
    NcValues,
};
use crate::error::{Error, Result};

/// `vsize` written for a variable too large for a 32-bit `vsize` field.
const VSIZE_OVERFLOW: u64 = 0xFFFF_FFFF;
/// Largest `vsize` a non-final variable may have in CDF-1 / CDF-2.
const VSIZE_LIMIT: u64 = 0xFFFF_FFFC;

/// Round `n` up to a multiple of four.
pub(crate) fn round_up4(n: u64) -> Option<u64> {
    n.checked_add(3).map(|v| v & !3)
}

/// Per-variable layout information shared by the writer and the reader.
#[derive(Debug, Clone, Copy)]
pub(crate) struct VarLayout {
    /// Unpadded bytes of the variable (of one record for a record variable).
    pub(crate) bytes: u64,
    /// `bytes` rounded up to a multiple of four.
    pub(crate) vsize: u64,
    /// Whether the variable is a record variable.
    pub(crate) is_record: bool,
}

/// Compute the layout of every variable from the dimension lengths.
///
/// `dim_lens` holds each dimension length (the record dimension's entry is
/// ignored) and `record_dim` the id of the unlimited dimension.
pub(crate) fn variable_layout(
    dim_lens: &[u64],
    record_dim: Option<usize>,
    dim_ids: &[usize],
    type_size: usize,
) -> Option<VarLayout> {
    let is_record = dim_ids.first().is_some_and(|&id| Some(id) == record_dim);
    let skip = usize::from(is_record);
    let mut count: u64 = 1;
    for &id in dim_ids.iter().skip(skip) {
        count = count.checked_mul(*dim_lens.get(id)?)?;
    }
    let bytes = count.checked_mul(type_size as u64)?;
    Some(VarLayout {
        bytes,
        vsize: round_up4(bytes)?,
        is_record,
    })
}

/// Record stride: the sum of the record variables' `vsize`, except that a
/// lone record variable is stored without padding between records.
pub(crate) fn record_size(layouts: &[VarLayout]) -> Option<u64> {
    let records: Vec<&VarLayout> = layouts.iter().filter(|l| l.is_record).collect();
    match records.as_slice() {
        [only] => Some(only.bytes),
        many => many
            .iter()
            .try_fold(0u64, |acc, l| acc.checked_add(l.vsize)),
    }
}

struct HeaderWriter {
    buf: Vec<u8>,
    version: CdfVersion,
}

impl HeaderWriter {
    fn put_u32(&mut self, v: u32) {
        self.buf.extend_from_slice(&v.to_be_bytes());
    }

    fn put_non_neg(&mut self, v: u64, what: &str) -> Result<()> {
        if v > self.version.max_count() {
            return Err(Error::General(format!(
                "{what} = {v} does not fit the {:?} NetCDF format",
                self.version
            )));
        }
        match self.version.count_width() {
            4 => self.put_u32(v as u32),
            _ => self.buf.extend_from_slice(&v.to_be_bytes()),
        }
        Ok(())
    }

    fn put_padding(&mut self) {
        while !self.buf.len().is_multiple_of(4) {
            self.buf.push(0);
        }
    }

    fn put_name(&mut self, name: &str) -> Result<()> {
        self.put_non_neg(name.len() as u64, "name length")?;
        self.buf.extend_from_slice(name.as_bytes());
        self.put_padding();
        Ok(())
    }

    fn put_list_header(&mut self, tag: u32, count: usize, what: &str) -> Result<()> {
        // An empty list is ABSENT: a zero tag followed by a zero count.
        self.put_u32(if count == 0 { 0 } else { tag });
        self.put_non_neg(count as u64, what)
    }

    fn put_attributes(&mut self, attrs: &[CdfAttribute]) -> Result<()> {
        self.put_list_header(NC_ATTRIBUTE, attrs.len(), "attribute count")?;
        for attr in attrs {
            self.put_name(&attr.name)?;
            self.put_u32(attr.values.nc_type().code());
            self.put_non_neg(attr.values.len() as u64, "attribute length")?;
            encode_values(&attr.values, &mut self.buf);
            self.put_padding();
        }
        Ok(())
    }

    /// Reserve a `begin` field and return its position.
    fn reserve_offset(&mut self) -> usize {
        let pos = self.buf.len();
        self.buf
            .extend(std::iter::repeat_n(0u8, self.version.offset_width()));
        pos
    }
}

/// Append the big-endian external representation of `values` to `out`.
pub(crate) fn encode_values(values: &NcValues, out: &mut Vec<u8>) {
    encode_range(values, 0, values.len(), out);
}

/// Append values `start..end` of `values` in big-endian order.
fn encode_range(values: &NcValues, start: usize, end: usize, out: &mut Vec<u8>) {
    macro_rules! put {
        ($v:expr) => {
            for x in &$v[start..end] {
                out.extend_from_slice(&x.to_be_bytes());
            }
        };
    }
    match values {
        NcValues::Byte(v) => put!(v),
        NcValues::Char(v) | NcValues::UByte(v) => out.extend_from_slice(&v[start..end]),
        NcValues::Short(v) => put!(v),
        NcValues::Int(v) => put!(v),
        NcValues::Float(v) => put!(v),
        NcValues::Double(v) => put!(v),
        NcValues::UShort(v) => put!(v),
        NcValues::UInt(v) => put!(v),
        NcValues::Int64(v) => put!(v),
        NcValues::UInt64(v) => put!(v),
    }
}

fn pad_to(out: &mut Vec<u8>, start: usize, len: u64) {
    let target = start + len as usize;
    if out.len() < target {
        out.resize(target, 0);
    }
}

fn check_types(ds: &CdfDataset, version: CdfVersion) -> Result<()> {
    let unsupported = |owner: &str, ty: super::NcType| {
        Error::General(format!(
            "{owner} has type {} which the {version:?} NetCDF format cannot store (use CDF-5)",
            ty.cdl_name()
        ))
    };
    for attr in &ds.attributes {
        if !version.supports(attr.values.nc_type()) {
            return Err(unsupported(
                &format!("global attribute '{}'", attr.name),
                attr.values.nc_type(),
            ));
        }
    }
    for var in &ds.variables {
        if !version.supports(var.nc_type()) {
            return Err(unsupported(
                &format!("variable '{}'", var.name),
                var.nc_type(),
            ));
        }
        for attr in &var.attributes {
            if !version.supports(attr.values.nc_type()) {
                return Err(unsupported(
                    &format!("attribute '{}:{}'", var.name, attr.name),
                    attr.values.nc_type(),
                ));
            }
        }
    }
    Ok(())
}

/// Encode `vsize`, applying the classic large-variable rule: in CDF-1/CDF-2 a
/// variable larger than 4 GiB − 4 is only allowed as the last fixed-size
/// variable (when there are no record variables) or the last record
/// variable, and its `vsize` field is written as `2^32 − 1`.
fn vsize_field(
    version: CdfVersion,
    var: &CdfVariable,
    layout: &VarLayout,
    may_be_large: bool,
) -> Result<u64> {
    if version == CdfVersion::Data64 || layout.vsize <= VSIZE_LIMIT {
        return Ok(layout.vsize);
    }
    if may_be_large {
        Ok(VSIZE_OVERFLOW)
    } else {
        Err(Error::General(format!(
            "variable '{}' needs {} bytes; in the {version:?} format only the last \
             variable may exceed 4 GiB (use CDF-5)",
            var.name, layout.vsize
        )))
    }
}

/// Serialise `ds` in the NetCDF classic binary format `version`.
///
/// Fixed-size variables are stored contiguously in definition order after
/// the header; record variables follow, interleaved record by record.  All
/// values are big-endian and every variable (or per-record slab) is padded
/// to a four-byte boundary, except the slabs of a lone record variable.
///
/// # Errors
///
/// Fails when [`CdfDataset::validate`] fails, a type needs CDF-5, or a size
/// or offset does not fit the chosen format.
pub fn write_cdf_dataset(ds: &CdfDataset, version: CdfVersion) -> Result<Vec<u8>> {
    ds.validate()?;
    check_types(ds, version)?;
    let record_dim = ds.record_dimension();
    let num_records = ds.num_records() as u64;
    let dim_lens: Vec<u64> = ds.dimensions.iter().map(|d| d.len as u64).collect();
    let layouts = ds
        .variables
        .iter()
        .map(|v| {
            variable_layout(&dim_lens, record_dim, &v.dim_ids, v.nc_type().size())
                .ok_or_else(|| Error::General(format!("variable '{}' is too large", v.name)))
        })
        .collect::<Result<Vec<VarLayout>>>()?;
    let n_record_vars = layouts.iter().filter(|l| l.is_record).count();
    let last_fixed = layouts.iter().rposition(|l| !l.is_record);
    let last_record = layouts.iter().rposition(|l| l.is_record);

    let mut hw = HeaderWriter {
        buf: Vec::new(),
        version,
    };
    hw.buf.extend_from_slice(&version.magic());
    hw.put_non_neg(num_records, "number of records")?;
    hw.put_list_header(NC_DIMENSION, ds.dimensions.len(), "dimension count")?;
    for dim in &ds.dimensions {
        hw.put_name(&dim.name)?;
        hw.put_non_neg(
            if dim.unlimited { 0 } else { dim.len as u64 },
            "dimension length",
        )?;
    }
    hw.put_attributes(&ds.attributes)?;
    hw.put_list_header(NC_VARIABLE, ds.variables.len(), "variable count")?;
    let mut begin_slots = Vec::with_capacity(ds.variables.len());
    for (i, (var, layout)) in ds.variables.iter().zip(&layouts).enumerate() {
        hw.put_name(&var.name)?;
        hw.put_non_neg(var.dim_ids.len() as u64, "dimension count")?;
        for &id in &var.dim_ids {
            hw.put_non_neg(id as u64, "dimension id")?;
        }
        hw.put_attributes(&var.attributes)?;
        hw.put_u32(var.nc_type().code());
        let may_be_large = if layout.is_record {
            Some(i) == last_record
        } else {
            n_record_vars == 0 && Some(i) == last_fixed
        };
        let vsize = vsize_field(version, var, layout, may_be_large)?;
        // vsize is a NON_NEG field but VSIZE_OVERFLOW is the documented
        // all-ones marker, so write it directly.
        if vsize == VSIZE_OVERFLOW && version != CdfVersion::Data64 {
            hw.put_u32(u32::MAX);
        } else {
            hw.put_non_neg(vsize, "variable size")?;
        }
        begin_slots.push(hw.reserve_offset());
    }

    let mut out = hw.buf;
    let mut begins = vec![0u64; ds.variables.len()];
    // Fixed-size variables, contiguous in definition order.
    for (i, (var, layout)) in ds.variables.iter().zip(&layouts).enumerate() {
        if layout.is_record {
            continue;
        }
        let start = out.len();
        begins[i] = start as u64;
        encode_values(&var.values, &mut out);
        pad_to(&mut out, start, layout.vsize);
    }
    // Record variables: begin is the offset of the slab in record 0.
    let record_start = out.len() as u64;
    let recsize =
        record_size(&layouts).ok_or_else(|| Error::General("record size overflows".to_string()))?;
    let mut offset_in_record = 0u64;
    for (i, layout) in layouts.iter().enumerate() {
        if layout.is_record {
            begins[i] = record_start
                .checked_add(offset_in_record)
                .ok_or_else(|| Error::General("record offset overflows".to_string()))?;
            offset_in_record = offset_in_record
                .checked_add(layout.vsize)
                .ok_or_else(|| Error::General("record offset overflows".to_string()))?;
        }
    }
    let record_vars: Vec<(&CdfVariable, &VarLayout)> = ds
        .variables
        .iter()
        .zip(&layouts)
        .filter(|(_, l)| l.is_record)
        .collect();
    for rec in 0..num_records {
        let rec_begin = out.len();
        for (var, layout) in &record_vars {
            let per_record = (layout.bytes / var.nc_type().size() as u64) as usize;
            let first = rec as usize * per_record;
            let slab_start = out.len();
            encode_range(&var.values, first, first + per_record, &mut out);
            if n_record_vars > 1 {
                pad_to(&mut out, slab_start, layout.vsize);
            }
        }
        debug_assert_eq!((out.len() - rec_begin) as u64, recsize);
    }

    for (slot, &begin) in begin_slots.iter().zip(&begins) {
        if begin > version.max_offset() {
            return Err(Error::General(format!(
                "data offset {begin} does not fit the {version:?} NetCDF format \
                 (use CdfVersion::Offset64)"
            )));
        }
        match version.offset_width() {
            4 => out[*slot..*slot + 4].copy_from_slice(&(begin as u32).to_be_bytes()),
            _ => out[*slot..*slot + 8].copy_from_slice(&begin.to_be_bytes()),
        }
    }
    Ok(out)
}
