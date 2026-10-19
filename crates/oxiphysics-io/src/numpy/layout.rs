// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Element layout of `.npy` arrays: dtype descriptors, byte order and memory
//! order.

/// Byte order of the elements of an `.npy` array.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ByteOrder {
    /// Least significant byte first (`<`).
    Little,
    /// Most significant byte first (`>`).
    Big,
}

impl ByteOrder {
    /// Byte order of the machine this code runs on (`=` in a NumPy dtype).
    pub(super) const NATIVE: ByteOrder = if cfg!(target_endian = "little") {
        ByteOrder::Little
    } else {
        ByteOrder::Big
    };
}

/// A parsed simple (non-structured) dtype descriptor such as `<f8`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct SimpleDescr {
    /// Byte order of the elements.  Irrelevant for one-byte types.
    pub(super) byte_order: ByteOrder,
    /// NumPy kind character: `f` float, `i` signed, `u` unsigned, `b` boolean.
    pub(super) kind: char,
    /// Element size in bytes.
    pub(super) size: usize,
}

/// Parse a simple dtype descriptor.
///
/// The optional first character selects the byte order: `<` little endian,
/// `>` big endian, `=` (or no marker at all) the native order and `|` "not
/// applicable", which NumPy only writes for one-byte types.  It is followed by
/// the kind character and the element size in bytes (`f8`, `i4`, `u1`, `b1`).
/// A bare `?` denotes a boolean.
///
/// # Errors
///
/// Returns a description of the problem for structured or otherwise
/// unsupported descriptors, a zero element size, or `|` on a multi-byte type.
pub(super) fn parse_simple_descr(descr: &str) -> Result<SimpleDescr, String> {
    let (marker, body) = match descr.chars().next() {
        Some(c @ ('<' | '>' | '=' | '|')) => (Some(c), &descr[c.len_utf8()..]),
        _ => (None, descr),
    };
    let mut chars = body.chars();
    let kind = chars
        .next()
        .ok_or_else(|| format!("dtype '{descr}' has no type character"))?;
    let digits = chars.as_str();
    let size = if kind == '?' && digits.is_empty() {
        1
    } else if !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit()) {
        digits
            .parse::<usize>()
            .map_err(|e| format!("bad element size in dtype '{descr}': {e}"))?
    } else {
        return Err(format!("unsupported dtype '{descr}'"));
    };
    if size == 0 {
        return Err(format!("dtype '{descr}' has a zero element size"));
    }
    let byte_order = match marker {
        Some('<') => ByteOrder::Little,
        Some('>') => ByteOrder::Big,
        Some('|') if size == 1 => ByteOrder::Little,
        Some('|') => {
            return Err(format!(
                "dtype '{descr}': the '|' marker is only valid for one-byte types"
            ));
        }
        _ => ByteOrder::NATIVE,
    };
    Ok(SimpleDescr {
        byte_order,
        kind,
        size,
    })
}

/// Number of elements of an array of the given shape, or an error when the
/// product does not fit in `usize`.
pub(super) fn checked_element_count(shape: &[usize]) -> Result<usize, String> {
    shape
        .iter()
        .try_fold(1usize, |count, &dim| count.checked_mul(dim))
        .ok_or_else(|| format!("array shape {shape:?} has more elements than can be addressed"))
}

/// Reorder `data`, stored in Fortran (column-major, first axis fastest) order
/// for `shape`, into C (row-major, last axis fastest) order.
///
/// # Errors
///
/// Returns an error when `data.len()` differs from the number of elements of
/// `shape`.
pub(super) fn fortran_to_c_order<T: Copy>(data: &[T], shape: &[usize]) -> Result<Vec<T>, String> {
    let total = checked_element_count(shape)?;
    if data.len() != total {
        return Err(format!(
            "shape {shape:?} needs {total} elements but {} were supplied",
            data.len()
        ));
    }
    if shape.len() < 2 {
        return Ok(data.to_vec());
    }
    // Offset step of each axis in the Fortran layout.
    let mut strides = Vec::with_capacity(shape.len());
    let mut stride = 1usize;
    for &dim in shape {
        strides.push(stride);
        stride = stride.saturating_mul(dim);
    }
    let mut index = vec![0usize; shape.len()];
    let mut offset = 0usize;
    let mut out = Vec::with_capacity(total);
    for _ in 0..total {
        let value = data
            .get(offset)
            .copied()
            .ok_or_else(|| format!("offset {offset} outside the array data"))?;
        out.push(value);
        // Advance the multi-index in C order (last axis fastest) and keep the
        // Fortran offset in step with it.
        for axis in (0..shape.len()).rev() {
            index[axis] += 1;
            offset += strides[axis];
            if index[axis] < shape[axis] {
                break;
            }
            offset -= strides[axis] * shape[axis];
            index[axis] = 0;
        }
    }
    Ok(out)
}
