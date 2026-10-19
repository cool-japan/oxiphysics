//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::layout::{ByteOrder, checked_element_count, fortran_to_c_order, parse_simple_descr};
use super::types::NpyDtype;

pub(super) const NPY_MAGIC: &[u8; 6] = b"\x93NUMPY";
pub(super) const NPY_MAJOR: u8 = 1;
pub(super) const NPY_MINOR: u8 = 0;
/// Validate that a shape and data length are consistent.
pub fn validate_shape(shape: &[usize], data_len: usize) -> Result<(), String> {
    let expected: usize = shape.iter().product();
    if expected != data_len {
        Err(format!(
            "shape {shape:?} requires {expected} elements but got {data_len}"
        ))
    } else {
        Ok(())
    }
}
/// Compute the flat index from multi-dimensional indices (row-major).
pub fn flat_index(indices: &[usize], shape: &[usize]) -> Result<usize, String> {
    if indices.len() != shape.len() {
        return Err(format!(
            "index dimensionality {} != shape dimensionality {}",
            indices.len(),
            shape.len()
        ));
    }
    let mut idx = 0usize;
    let mut stride = 1usize;
    for i in (0..shape.len()).rev() {
        if indices[i] >= shape[i] {
            return Err(format!(
                "index {} out of range for axis {} with size {}",
                indices[i], i, shape[i]
            ));
        }
        idx += indices[i] * stride;
        stride *= shape[i];
    }
    Ok(idx)
}
/// Compute multi-dimensional indices from a flat index (row-major).
pub fn unravel_index(flat: usize, shape: &[usize]) -> Result<Vec<usize>, String> {
    let total: usize = shape.iter().product();
    if flat >= total {
        return Err(format!("flat index {flat} out of range for total {total}"));
    }
    let mut indices = vec![0usize; shape.len()];
    let mut remaining = flat;
    for i in (0..shape.len()).rev() {
        indices[i] = remaining % shape[i];
        remaining /= shape[i];
    }
    Ok(indices)
}
/// Build the NPY v1.0 header string for the given dtype and shape.
pub(super) fn build_npy_header(dtype_str: &str, shape: &[usize]) -> Vec<u8> {
    let shape_str = if shape.is_empty() {
        "()".to_string()
    } else if shape.len() == 1 {
        format!("({},)", shape[0])
    } else {
        let inner: Vec<String> = shape.iter().map(|d| d.to_string()).collect();
        format!("({})", inner.join(", "))
    };
    let dict = format!(
        "{{'descr': '{}', 'fortran_order': False, 'shape': {}, }}",
        dtype_str, shape_str
    );
    let mut header_bytes = dict.into_bytes();
    header_bytes.push(b'\n');
    while (10 + header_bytes.len()) % 64 != 0 {
        let last = header_bytes.len() - 1;
        header_bytes.insert(last, b' ');
    }
    header_bytes
}
/// Assemble a complete `.npy` v1.0 byte sequence for a `f64` array.
pub fn write_npy_f64(shape: &[usize], data: &[f64]) -> Vec<u8> {
    let header_bytes = build_npy_header("<f8", shape);
    let header_len = header_bytes.len() as u16;
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(NPY_MAGIC);
    out.push(NPY_MAJOR);
    out.push(NPY_MINOR);
    out.extend_from_slice(&header_len.to_le_bytes());
    out.extend_from_slice(&header_bytes);
    for &v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}
/// Assemble a complete `.npy` v1.0 byte sequence for an `f32` array.
pub fn write_npy_f32(shape: &[usize], data: &[f32]) -> Vec<u8> {
    let header_bytes = build_npy_header("<f4", shape);
    let header_len = header_bytes.len() as u16;
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(NPY_MAGIC);
    out.push(NPY_MAJOR);
    out.push(NPY_MINOR);
    out.extend_from_slice(&header_len.to_le_bytes());
    out.extend_from_slice(&header_bytes);
    for &v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}
/// Assemble a complete `.npy` v1.0 byte sequence for an `i32` array.
pub fn write_npy_i32(shape: &[usize], data: &[i32]) -> Vec<u8> {
    let header_bytes = build_npy_header("<i4", shape);
    let header_len = header_bytes.len() as u16;
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(NPY_MAGIC);
    out.push(NPY_MAJOR);
    out.push(NPY_MINOR);
    out.extend_from_slice(&header_len.to_le_bytes());
    out.extend_from_slice(&header_bytes);
    for &v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}
/// Assemble a complete `.npy` v1.0 byte sequence for an `i64` array.
pub fn write_npy_i64(shape: &[usize], data: &[i64]) -> Vec<u8> {
    let header_bytes = build_npy_header("<i8", shape);
    let header_len = header_bytes.len() as u16;
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(NPY_MAGIC);
    out.push(NPY_MAJOR);
    out.push(NPY_MINOR);
    out.extend_from_slice(&header_len.to_le_bytes());
    out.extend_from_slice(&header_bytes);
    for &v in data {
        out.extend_from_slice(&v.to_le_bytes());
    }
    out
}
/// Parsed header of an `.npy` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct NpyHeader {
    /// The `descr` entry, for example `<f8`.
    pub(super) descr: String,
    /// Whether the data are stored in Fortran (column-major) order.
    pub(super) fortran_order: bool,
    /// Logical array shape.
    pub(super) shape: Vec<usize>,
    /// Byte offset of the first data element.
    pub(super) data_start: usize,
}
/// Parse the header of an `.npy` file (format versions 1.0, 2.0 and 3.0).
pub(super) fn parse_npy_header(bytes: &[u8]) -> Result<NpyHeader, String> {
    let magic = bytes.get(..6).ok_or("npy data too short")?;
    if magic != NPY_MAGIC {
        return Err(format!("bad npy magic: {magic:?}"));
    }
    let (major, minor) = match (bytes.get(6), bytes.get(7)) {
        (Some(&major), Some(&minor)) => (major, minor),
        _ => return Err("npy data too short".to_string()),
    };
    // Versions 2.0 and 3.0 widen the header length field from 16 to 32 bits;
    // version 3.0 additionally encodes the header as UTF-8 (Latin-1 before).
    let length_field_bytes = match (major, minor) {
        (1, 0) => 2,
        (2, 0) | (3, 0) => 4,
        _ => return Err(format!("unsupported npy version: {major}.{minor}")),
    };
    let header_start = 8 + length_field_bytes;
    let length_field = bytes.get(8..header_start).ok_or("npy data too short")?;
    let header_len = match length_field {
        [a, b] => usize::from(u16::from_le_bytes([*a, *b])),
        [a, b, c, d] => usize::try_from(u32::from_le_bytes([*a, *b, *c, *d]))
            .map_err(|_| "npy header length does not fit in memory".to_string())?,
        _ => return Err("npy data too short".to_string()),
    };
    let data_start = header_start
        .checked_add(header_len)
        .ok_or("npy header length overflows")?;
    let raw_header = bytes
        .get(header_start..data_start)
        .ok_or("npy header truncated")?;
    let header_text: String = if major >= 3 {
        std::str::from_utf8(raw_header)
            .map_err(|e| format!("npy header not utf-8: {e}"))?
            .to_string()
    } else {
        raw_header.iter().map(|&b| char::from(b)).collect()
    };
    let header_str = header_text.trim();
    let descr = extract_dict_value(header_str, "descr")?;
    let fortran_order = match extract_dict_value(header_str, "fortran_order")?.as_str() {
        "True" => true,
        "False" => false,
        other => return Err(format!("bad fortran_order value '{other}' in npy header")),
    };
    let shape = parse_shape_tuple(&extract_dict_value(header_str, "shape")?)?;
    Ok(NpyHeader {
        descr,
        fortran_order,
        shape,
        data_start,
    })
}
/// Extract a value from a Python-style dict string for a given key.
pub(super) fn extract_dict_value(header: &str, key: &str) -> Result<String, String> {
    let search = format!("'{key}'");
    let pos = header
        .find(&search)
        .ok_or_else(|| format!("key '{key}' not found in npy header"))?;
    let rest = header[pos + search.len()..].trim_start();
    let rest = rest
        .strip_prefix(':')
        .ok_or("missing ':' after key")?
        .trim_start();
    if let Some(inner) = rest.strip_prefix('\'') {
        let end = inner.find('\'').ok_or("unterminated string value")?;
        Ok(inner[..end].to_string())
    } else if rest.starts_with(['(', '[']) {
        Ok(rest[..bracketed_len(rest)?].to_string())
    } else {
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        Ok(rest[..end].trim().to_string())
    }
}
/// Byte length of the bracketed group `(...)` or `[...]` that starts at the
/// beginning of `s`, nested groups included.
fn bracketed_len(s: &str) -> Result<usize, String> {
    let mut depth = 0usize;
    for (i, c) in s.char_indices() {
        match c {
            '(' | '[' => depth += 1,
            ')' | ']' => {
                depth = depth
                    .checked_sub(1)
                    .ok_or("unbalanced brackets in npy header")?;
                if depth == 0 {
                    return Ok(i + c.len_utf8());
                }
            }
            _ => {}
        }
    }
    Err("unterminated bracketed value in npy header".to_string())
}
/// Parse a Python tuple string like `"(3, 2)"` or `"(6,)"` or `"()"` into `Vec`usize`.
pub(super) fn parse_shape_tuple(s: &str) -> Result<Vec<usize>, String> {
    let inner = s.trim();
    let inner = inner
        .strip_prefix('(')
        .ok_or("shape missing '('")?
        .strip_suffix(')')
        .ok_or("shape missing ')'")?;
    if inner.trim().is_empty() {
        return Ok(vec![]);
    }
    let mut dims = Vec::new();
    for part in inner.split(',') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let d: usize = part
            .parse()
            .map_err(|e| format!("bad shape dimension '{part}': {e}"))?;
        dims.push(d);
    }
    Ok(dims)
}
/// Decode the elements of an `.npy` array whose elements are `N` bytes wide.
///
/// `expected` is the canonical dtype string used in error messages and `kind`
/// the NumPy kind character the caller accepts.  Both byte orders are
/// supported.  Data stored in Fortran order (`fortran_order: True`, as written
/// by `numpy.save` for transposed and other column-major arrays) are returned
/// in C (row-major) order, so the returned vector always matches the returned
/// shape in row-major layout.
fn read_npy_elements<T: Copy, const N: usize>(
    bytes: &[u8],
    expected: &str,
    kind: char,
    decode_le: fn([u8; N]) -> T,
    decode_be: fn([u8; N]) -> T,
) -> Result<(Vec<usize>, Vec<T>), String> {
    let header = parse_npy_header(bytes)?;
    let wrong_dtype = || format!("expected dtype '{expected}', got '{}'", header.descr);
    let descr = parse_simple_descr(&header.descr).map_err(|_| wrong_dtype())?;
    if descr.kind != kind || descr.size != N {
        return Err(wrong_dtype());
    }
    let n_elems = checked_element_count(&header.shape)?;
    let payload_len = n_elems.checked_mul(N).ok_or("npy payload size overflows")?;
    let end = header
        .data_start
        .checked_add(payload_len)
        .ok_or("npy payload size overflows")?;
    let payload = bytes
        .get(header.data_start..end)
        .ok_or_else(|| format!("data truncated: expected {end} bytes, got {}", bytes.len()))?;
    let decode = match descr.byte_order {
        ByteOrder::Little => decode_le,
        ByteOrder::Big => decode_be,
    };
    // `payload` holds exactly `n_elems * N` bytes, so this allocation is
    // bounded by the size of the input.
    let mut data = Vec::with_capacity(n_elems);
    for chunk in payload.chunks_exact(N) {
        let raw = <[u8; N]>::try_from(chunk)
            .map_err(|_| "internal error: element chunk has the wrong width".to_string())?;
        data.push(decode(raw));
    }
    if header.fortran_order {
        data = fortran_to_c_order(&data, &header.shape)?;
    }
    Ok((header.shape, data))
}
/// Parse a `.npy` byte buffer and return `(shape, f64 data)`.
///
/// Accepts `<f8`, `>f8` and `=f8` arrays in C or Fortran order; the data are
/// returned in C order.
pub fn read_npy_f64(bytes: &[u8]) -> Result<(Vec<usize>, Vec<f64>), String> {
    read_npy_elements::<f64, 8>(bytes, "<f8", 'f', f64::from_le_bytes, f64::from_be_bytes)
}
/// Parse a `.npy` byte buffer and return `(shape, f32 data)`.
///
/// Accepts `<f4`, `>f4` and `=f4` arrays in C or Fortran order; the data are
/// returned in C order.
pub fn read_npy_f32(bytes: &[u8]) -> Result<(Vec<usize>, Vec<f32>), String> {
    read_npy_elements::<f32, 4>(bytes, "<f4", 'f', f32::from_le_bytes, f32::from_be_bytes)
}
/// Parse a `.npy` byte buffer and return `(shape, i32 data)`.
///
/// Accepts `<i4`, `>i4` and `=i4` arrays in C or Fortran order; the data are
/// returned in C order.
pub fn read_npy_i32(bytes: &[u8]) -> Result<(Vec<usize>, Vec<i32>), String> {
    read_npy_elements::<i32, 4>(bytes, "<i4", 'i', i32::from_le_bytes, i32::from_be_bytes)
}
/// Parse a `.npy` byte buffer and return `(shape, i64 data)`.
///
/// Accepts `<i8`, `>i8` and `=i8` arrays in C or Fortran order; the data are
/// returned in C order.
pub fn read_npy_i64(bytes: &[u8]) -> Result<(Vec<usize>, Vec<i64>), String> {
    read_npy_elements::<i64, 8>(bytes, "<i8", 'i', i64::from_le_bytes, i64::from_be_bytes)
}
/// Auto-detect dtype from header and return the NpyDtype.
///
/// The byte order of the stored elements is not part of the result; the typed
/// readers handle both byte orders.
pub fn detect_npy_dtype(bytes: &[u8]) -> Result<NpyDtype, String> {
    let header = parse_npy_header(bytes)?;
    NpyDtype::from_numpy_str(&header.descr)
}
/// Auto-detect and return the shape from a .npy byte buffer.
pub fn read_npy_shape(bytes: &[u8]) -> Result<Vec<usize>, String> {
    Ok(parse_npy_header(bytes)?.shape)
}
pub(super) fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32, String> {
    let raw = pos
        .checked_add(4)
        .and_then(|end| data.get(*pos..end))
        .and_then(|slice| <[u8; 4]>::try_from(slice).ok())
        .ok_or_else(|| format!("unexpected EOF reading u32 at offset {pos}"))?;
    *pos += 4;
    Ok(u32::from_le_bytes(raw))
}
/// Compute the arithmetic mean of a slice.
///
/// Returns `None` if the slice is empty.
pub fn slice_mean(data: &[f64]) -> Option<f64> {
    if data.is_empty() {
        return None;
    }
    Some(data.iter().sum::<f64>() / data.len() as f64)
}
/// Compute the variance of a slice (population variance, ddof=0).
pub fn slice_var(data: &[f64]) -> Option<f64> {
    let mean = slice_mean(data)?;
    let var = data.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / data.len() as f64;
    Some(var)
}
/// Compute the standard deviation (population, ddof=0).
pub fn slice_std(data: &[f64]) -> Option<f64> {
    Some(slice_var(data)?.sqrt())
}
/// Compute min, max, and their flat indices.
pub fn slice_min_max(data: &[f64]) -> Option<(f64, usize, f64, usize)> {
    if data.is_empty() {
        return None;
    }
    let mut min_val = data[0];
    let mut max_val = data[0];
    let mut min_idx = 0;
    let mut max_idx = 0;
    for (i, &v) in data.iter().enumerate() {
        if v < min_val {
            min_val = v;
            min_idx = i;
        }
        if v > max_val {
            max_val = v;
            max_idx = i;
        }
    }
    Some((min_val, min_idx, max_val, max_idx))
}
/// Compute the p-th percentile of a slice using linear interpolation.
///
/// `p` must be in `\[0, 100\]`.
pub fn slice_percentile(data: &[f64], p: f64) -> std::result::Result<f64, String> {
    if data.is_empty() {
        return Err("slice_percentile: empty slice".to_string());
    }
    if !(0.0..=100.0).contains(&p) {
        return Err(format!("percentile p={p} not in [0,100]"));
    }
    let mut sorted = data.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = sorted.len();
    let idx = p / 100.0 * (n - 1) as f64;
    let lo = idx.floor() as usize;
    let hi = idx.ceil() as usize;
    if lo == hi {
        return Ok(sorted[lo]);
    }
    let frac = idx - lo as f64;
    Ok(sorted[lo] * (1.0 - frac) + sorted[hi] * frac)
}
/// Clip values to `\[lo, hi\]`.
pub fn slice_clip(data: &[f64], lo: f64, hi: f64) -> Vec<f64> {
    data.iter().map(|&v| v.clamp(lo, hi)).collect()
}
/// Element-wise sum of two equal-length slices.
pub fn slice_add(a: &[f64], b: &[f64]) -> std::result::Result<Vec<f64>, String> {
    if a.len() != b.len() {
        return Err(format!(
            "slice_add: length mismatch {} vs {}",
            a.len(),
            b.len()
        ));
    }
    Ok(a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect())
}
/// Element-wise product of two equal-length slices.
pub fn slice_mul(a: &[f64], b: &[f64]) -> std::result::Result<Vec<f64>, String> {
    if a.len() != b.len() {
        return Err(format!(
            "slice_mul: length mismatch {} vs {}",
            a.len(),
            b.len()
        ));
    }
    Ok(a.iter().zip(b.iter()).map(|(&x, &y)| x * y).collect())
}
/// Dot product of two equal-length slices.
pub fn slice_dot(a: &[f64], b: &[f64]) -> std::result::Result<f64, String> {
    Ok(slice_mul(a, b)?.iter().sum())
}
/// Generate `n` equally-spaced values from `start` to `stop` (inclusive).
///
/// Equivalent to `numpy.linspace(start, stop, num=n)`.
pub fn linspace(start: f64, stop: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (stop - start) * i as f64 / (n - 1) as f64)
        .collect()
}
/// Generate values from `start` to `stop` (exclusive) with step `step`.
///
/// Equivalent to `numpy.arange(start, stop, step)`.
pub fn arange(start: f64, stop: f64, step: f64) -> std::result::Result<Vec<f64>, String> {
    if step == 0.0 {
        return Err("arange: step cannot be zero".to_string());
    }
    if (stop - start) / step < 0.0 {
        return Ok(Vec::new());
    }
    let n = ((stop - start) / step).ceil() as usize;
    Ok((0..n).map(|i| start + i as f64 * step).collect())
}
/// Generate `n` log-spaced values from `10^start` to `10^stop`.
///
/// Equivalent to `numpy.logspace(start, stop, num=n)`.
pub fn logspace(start: f64, stop: f64, n: usize) -> Vec<f64> {
    linspace(start, stop, n)
        .into_iter()
        .map(|v| 10.0_f64.powf(v))
        .collect()
}
/// Transpose a 2-D row-major matrix stored as a flat `Vec`f64`.
///
/// `shape` must be `[nrows, ncols]`. Returns `(transposed_data, [ncols, nrows])`.
pub fn transpose_2d(
    data: &[f64],
    shape: &[usize],
) -> std::result::Result<(Vec<f64>, Vec<usize>), String> {
    if shape.len() != 2 {
        return Err(format!(
            "transpose_2d requires 2-D shape, got {}D",
            shape.len()
        ));
    }
    let nrows = shape[0];
    let ncols = shape[1];
    if data.len() != nrows * ncols {
        return Err(format!(
            "transpose_2d: data length {} != {}*{}",
            data.len(),
            nrows,
            ncols
        ));
    }
    let mut out = vec![0.0_f64; nrows * ncols];
    for r in 0..nrows {
        for c in 0..ncols {
            out[c * nrows + r] = data[r * ncols + c];
        }
    }
    Ok((out, vec![ncols, nrows]))
}
/// Compute the matrix product C = A * B where A is (m×k) and B is (k×n),
/// both stored as flat row-major `f64` slices.
pub fn matmul(
    a: &[f64],
    a_shape: &[usize],
    b: &[f64],
    b_shape: &[usize],
) -> std::result::Result<(Vec<f64>, Vec<usize>), String> {
    if a_shape.len() != 2 || b_shape.len() != 2 {
        return Err("matmul: both inputs must be 2-D".to_string());
    }
    let (m, k_a) = (a_shape[0], a_shape[1]);
    let (k_b, n) = (b_shape[0], b_shape[1]);
    if k_a != k_b {
        return Err(format!(
            "matmul: inner dimensions mismatch ({k_a} vs {k_b})"
        ));
    }
    if a.len() != m * k_a || b.len() != k_b * n {
        return Err("matmul: data length does not match shape".to_string());
    }
    let mut c = vec![0.0_f64; m * n];
    for i in 0..m {
        for j in 0..n {
            let mut s = 0.0_f64;
            for kk in 0..k_a {
                s += a[i * k_a + kk] * b[kk * n + j];
            }
            c[i * n + j] = s;
        }
    }
    Ok((c, vec![m, n]))
}
