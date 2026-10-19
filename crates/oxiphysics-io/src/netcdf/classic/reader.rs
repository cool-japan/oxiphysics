// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Decoder for the NetCDF classic binary format.

use super::writer::{VarLayout, record_size, round_up4, variable_layout};
use super::{
    CdfAttribute, CdfDataset, CdfDimension, CdfVariable, CdfVersion, NC_ATTRIBUTE, NC_DIMENSION,
    NC_VARIABLE, NcType, NcValues,
};
use crate::error::{Error, Result};

fn parse_err(msg: impl Into<String>) -> Error {
    Error::Parse(format!("NetCDF classic: {}", msg.into()))
}

struct Cursor<'a> {
    bytes: &'a [u8],
    pos: usize,
    version: CdfVersion,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize, what: &str) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(n)
            .filter(|&e| e <= self.bytes.len())
            .ok_or_else(|| parse_err(format!("truncated header while reading {what}")))?;
        let slice = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(slice)
    }

    fn remaining(&self) -> usize {
        self.bytes.len().saturating_sub(self.pos)
    }

    fn u32(&mut self, what: &str) -> Result<u32> {
        let b = self.take(4, what)?;
        Ok(u32::from_be_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn u64(&mut self, what: &str) -> Result<u64> {
        let b = self.take(8, what)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(b);
        Ok(u64::from_be_bytes(a))
    }

    /// A `NON_NEG` count (4 bytes, or 8 in CDF-5), returned raw.
    fn non_neg_raw(&mut self, what: &str) -> Result<u64> {
        match self.version.count_width() {
            4 => self.u32(what).map(u64::from),
            _ => {
                let v = self.u64(what)?;
                if v > i64::MAX as u64 {
                    return Err(parse_err(format!("negative {what}")));
                }
                Ok(v)
            }
        }
    }

    fn non_neg(&mut self, what: &str) -> Result<usize> {
        let v = self.non_neg_raw(what)?;
        usize::try_from(v).map_err(|_| parse_err(format!("{what} {v} is too large")))
    }

    /// A count of header elements, each of which needs at least `min_bytes`.
    fn element_count(&mut self, what: &str, min_bytes: usize) -> Result<usize> {
        let n = self.non_neg(what)?;
        if n.checked_mul(min_bytes)
            .is_none_or(|need| need > self.remaining())
        {
            return Err(parse_err(format!(
                "{what} {n} exceeds the remaining header bytes"
            )));
        }
        Ok(n)
    }

    fn offset(&mut self, what: &str) -> Result<u64> {
        match self.version.offset_width() {
            4 => {
                let v = self.u32(what)?;
                if v > i32::MAX as u32 {
                    return Err(parse_err(format!("negative {what}")));
                }
                Ok(u64::from(v))
            }
            _ => {
                let v = self.u64(what)?;
                if v > i64::MAX as u64 {
                    return Err(parse_err(format!("negative {what}")));
                }
                Ok(v)
            }
        }
    }

    fn skip_padding(&mut self, len: usize, what: &str) -> Result<()> {
        let padded =
            round_up4(len as u64).ok_or_else(|| parse_err(format!("{what} too long")))? as usize;
        self.take(padded - len, what).map(|_| ())
    }

    fn name(&mut self, what: &str) -> Result<String> {
        let len = self.element_count(what, 1)?;
        let raw = self.take(len, what)?;
        let name = std::str::from_utf8(raw)
            .map_err(|_| parse_err(format!("{what} is not valid UTF-8")))?
            .to_string();
        self.skip_padding(len, what)?;
        Ok(name)
    }

    /// Read a list header: returns the element count, checking the tag.
    fn list_header(&mut self, tag: u32, what: &str, min_bytes: usize) -> Result<usize> {
        let found = self.u32(what)?;
        if found == 0 {
            let n = self.non_neg(what)?;
            if n != 0 {
                return Err(parse_err(format!("ABSENT {what} with non-zero count {n}")));
            }
            return Ok(0);
        }
        if found != tag {
            return Err(parse_err(format!(
                "expected {what} tag {tag:#x}, found {found:#x}"
            )));
        }
        self.element_count(what, min_bytes)
    }

    fn nc_type(&mut self, what: &str) -> Result<NcType> {
        let code = self.u32(what)?;
        let ty = NcType::from_code(code)
            .ok_or_else(|| parse_err(format!("unknown nc_type {code} for {what}")))?;
        if !self.version.supports(ty) {
            return Err(parse_err(format!(
                "nc_type {} is not allowed in {:?} files",
                ty.cdl_name(),
                self.version
            )));
        }
        Ok(ty)
    }

    fn attributes(&mut self, owner: &str) -> Result<Vec<CdfAttribute>> {
        let what = format!("attribute list of {owner}");
        // Smallest attribute: name (count + 0 bytes) + type + count.
        let min = 2 * self.version.count_width() + 4;
        let n = self.list_header(NC_ATTRIBUTE, &what, min)?;
        let mut attrs = Vec::with_capacity(n);
        for _ in 0..n {
            let name = self.name("attribute name")?;
            let ty = self.nc_type("attribute type")?;
            let count = self.non_neg("attribute length")?;
            let nbytes = count
                .checked_mul(ty.size())
                .ok_or_else(|| parse_err(format!("attribute '{name}' is too large")))?;
            let raw = self.take(nbytes, "attribute values")?;
            let values = decode_values(ty, raw);
            self.skip_padding(nbytes, "attribute padding")?;
            attrs.push(CdfAttribute { name, values });
        }
        Ok(attrs)
    }
}

/// Decode big-endian external values of type `ty` (`raw.len()` must be a
/// multiple of the type size).
pub(crate) fn decode_values(ty: NcType, raw: &[u8]) -> NcValues {
    macro_rules! get {
        ($variant:ident, $t:ty, $n:expr) => {
            NcValues::$variant(
                raw.chunks_exact($n)
                    .map(|c| {
                        let mut a = [0u8; $n];
                        a.copy_from_slice(c);
                        <$t>::from_be_bytes(a)
                    })
                    .collect(),
            )
        };
    }
    match ty {
        NcType::Byte => NcValues::Byte(raw.iter().map(|&b| b as i8).collect()),
        NcType::Char => NcValues::Char(raw.to_vec()),
        NcType::UByte => NcValues::UByte(raw.to_vec()),
        NcType::Short => get!(Short, i16, 2),
        NcType::UShort => get!(UShort, u16, 2),
        NcType::Int => get!(Int, i32, 4),
        NcType::UInt => get!(UInt, u32, 4),
        NcType::Float => get!(Float, f32, 4),
        NcType::Double => get!(Double, f64, 8),
        NcType::Int64 => get!(Int64, i64, 8),
        NcType::UInt64 => get!(UInt64, u64, 8),
    }
}

/// Append decoded values to an existing array of the same type.
fn extend_values(dst: &mut NcValues, src: NcValues) {
    match (dst, src) {
        (NcValues::Byte(a), NcValues::Byte(b)) => a.extend(b),
        (NcValues::Char(a), NcValues::Char(b)) => a.extend(b),
        (NcValues::Short(a), NcValues::Short(b)) => a.extend(b),
        (NcValues::Int(a), NcValues::Int(b)) => a.extend(b),
        (NcValues::Float(a), NcValues::Float(b)) => a.extend(b),
        (NcValues::Double(a), NcValues::Double(b)) => a.extend(b),
        (NcValues::UByte(a), NcValues::UByte(b)) => a.extend(b),
        (NcValues::UShort(a), NcValues::UShort(b)) => a.extend(b),
        (NcValues::UInt(a), NcValues::UInt(b)) => a.extend(b),
        (NcValues::Int64(a), NcValues::Int64(b)) => a.extend(b),
        (NcValues::UInt64(a), NcValues::UInt64(b)) => a.extend(b),
        // Callers always pass matching types (both come from the same nc_type).
        _ => {}
    }
}

struct VarHeader {
    name: String,
    dim_ids: Vec<usize>,
    attributes: Vec<CdfAttribute>,
    nc_type: NcType,
    begin: u64,
}

fn data_slice<'a>(bytes: &'a [u8], begin: u64, len: u64, name: &str) -> Result<&'a [u8]> {
    let end = begin
        .checked_add(len)
        .filter(|&e| e <= bytes.len() as u64)
        .ok_or_else(|| parse_err(format!("data of variable '{name}' lies beyond end of file")))?;
    Ok(&bytes[begin as usize..end as usize])
}

/// Parse a NetCDF classic (`CDF-1`), 64-bit offset (`CDF-2`) or 64-bit data
/// (`CDF-5`) file.
///
/// A streaming header (`numrecs` = all ones) is resolved from the file
/// length.  Every count, offset and size is bounds-checked, so malformed or
/// truncated input yields an error instead of a panic or a huge allocation.
///
/// # Errors
///
/// [`Error::Parse`] for a bad magic number, an unknown type or tag, a
/// truncated header or data section, inconsistent dimensions, or data that
/// lies outside the file.
pub fn read_cdf_dataset(bytes: &[u8]) -> Result<CdfDataset> {
    let magic = bytes
        .get(..4)
        .ok_or_else(|| parse_err("file shorter than the 4-byte magic number"))?;
    if &magic[..3] != b"CDF" {
        return Err(parse_err("missing 'CDF' magic number"));
    }
    let version = CdfVersion::from_version_byte(magic[3])
        .ok_or_else(|| parse_err(format!("unsupported format version byte {}", magic[3])))?;
    let mut cur = Cursor {
        bytes,
        pos: 4,
        version,
    };
    let numrecs_raw = cur.non_neg_raw_or_streaming()?;

    // Dimensions.
    let cw = version.count_width();
    let n_dims = cur.list_header(NC_DIMENSION, "dimension list", 2 * cw)?;
    let mut dimensions = Vec::with_capacity(n_dims);
    for _ in 0..n_dims {
        let name = cur.name("dimension name")?;
        let len = cur.non_neg("dimension length")?;
        dimensions.push(CdfDimension {
            name,
            len,
            unlimited: len == 0,
        });
    }
    let record_dim = dimensions.iter().position(|d| d.unlimited);
    if dimensions.iter().filter(|d| d.unlimited).count() > 1 {
        return Err(parse_err("more than one unlimited dimension"));
    }

    let attributes = cur.attributes("the dataset")?;

    // Variables.
    let min_var = 3 * cw + 8 + 4 + cw + version.offset_width();
    let n_vars = cur.list_header(NC_VARIABLE, "variable list", min_var)?;
    let mut headers = Vec::with_capacity(n_vars);
    for _ in 0..n_vars {
        let name = cur.name("variable name")?;
        let n = cur.element_count("variable rank", cw)?;
        let mut dim_ids = Vec::with_capacity(n);
        for _ in 0..n {
            let id = cur.non_neg("dimension id")?;
            if id >= dimensions.len() {
                return Err(parse_err(format!(
                    "variable '{name}' refers to missing dimension id {id}"
                )));
            }
            dim_ids.push(id);
        }
        let attributes = cur.attributes(&name)?;
        let nc_type = cur.nc_type("variable type")?;
        // vsize is redundant (and saturated for huge variables); sizes are
        // recomputed from the shape below.
        cur.non_neg_raw("variable size")?;
        let begin = cur.offset("variable begin")?;
        headers.push(VarHeader {
            name,
            dim_ids,
            attributes,
            nc_type,
            begin,
        });
    }
    let header_end = cur.pos as u64;

    let dim_lens: Vec<u64> = dimensions.iter().map(|d| d.len as u64).collect();
    let layouts = headers
        .iter()
        .map(|h| {
            if let Some(pos) = h
                .dim_ids
                .iter()
                .skip(1)
                .position(|&id| Some(id) == record_dim)
            {
                return Err(parse_err(format!(
                    "variable '{}' uses the unlimited dimension at position {}",
                    h.name,
                    pos + 1
                )));
            }
            variable_layout(&dim_lens, record_dim, &h.dim_ids, h.nc_type.size())
                .ok_or_else(|| parse_err(format!("variable '{}' is too large", h.name)))
        })
        .collect::<Result<Vec<VarLayout>>>()?;
    let recsize = record_size(&layouts).ok_or_else(|| parse_err("record size overflows"))?;
    let first_record_begin = headers
        .iter()
        .zip(&layouts)
        .filter(|(_, l)| l.is_record)
        .map(|(h, _)| h.begin)
        .min();

    let num_records: u64 = match numrecs_raw {
        Some(n) => n,
        None => match first_record_begin {
            Some(begin) if recsize > 0 => (bytes.len() as u64).saturating_sub(begin) / recsize,
            _ => 0,
        },
    };
    // Without an unlimited dimension `numrecs` describes nothing and is
    // ignored, as netCDF-C does.
    if let Some(id) = record_dim {
        dimensions[id].len = usize::try_from(num_records)
            .map_err(|_| parse_err("number of records is too large"))?;
    }

    let mut variables = Vec::with_capacity(headers.len());
    for (h, layout) in headers.into_iter().zip(layouts) {
        if h.begin < header_end {
            return Err(parse_err(format!(
                "data of variable '{}' overlaps the header",
                h.name
            )));
        }
        let values = if layout.is_record {
            let mut values = NcValues::empty(h.nc_type);
            if num_records > 0 {
                let last_start = (num_records - 1)
                    .checked_mul(recsize)
                    .and_then(|o| o.checked_add(h.begin))
                    .ok_or_else(|| parse_err("record offset overflows"))?;
                data_slice(bytes, last_start, layout.bytes, &h.name)?;
                for rec in 0..num_records {
                    let start = h.begin + rec * recsize;
                    let raw = data_slice(bytes, start, layout.bytes, &h.name)?;
                    extend_values(&mut values, decode_values(h.nc_type, raw));
                }
            }
            values
        } else {
            decode_values(
                h.nc_type,
                data_slice(bytes, h.begin, layout.bytes, &h.name)?,
            )
        };
        variables.push(CdfVariable {
            name: h.name,
            dim_ids: h.dim_ids,
            attributes: h.attributes,
            values,
        });
    }

    let ds = CdfDataset {
        dimensions,
        attributes,
        variables,
    };
    ds.validate()
        .map_err(|e| parse_err(format!("invalid dataset: {e}")))?;
    Ok(ds)
}

impl Cursor<'_> {
    /// `numrecs`: a count, or `None` for the streaming marker (all ones).
    fn non_neg_raw_or_streaming(&mut self) -> Result<Option<u64>> {
        match self.version.count_width() {
            4 => {
                let v = self.u32("number of records")?;
                Ok((v != u32::MAX).then_some(u64::from(v)))
            }
            _ => {
                let v = self.u64("number of records")?;
                if v == u64::MAX {
                    Ok(None)
                } else if v > i64::MAX as u64 {
                    Err(parse_err("negative number of records"))
                } else {
                    Ok(Some(v))
                }
            }
        }
    }
}
