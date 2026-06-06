//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

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
/// Parse the NPY v1.0 header and return `(dtype_str, shape, data_offset)`.
pub(super) fn parse_npy_header(bytes: &[u8]) -> Result<(String, Vec<usize>, usize), String> {
    if bytes.len() < 10 {
        return Err("npy data too short".to_string());
    }
    if &bytes[0..6] != NPY_MAGIC {
        return Err(format!("bad npy magic: {:?}", &bytes[0..6]));
    }
    let major = bytes[6];
    let minor = bytes[7];
    if major != 1 || minor != 0 {
        return Err(format!("unsupported npy version: {major}.{minor}"));
    }
    let header_len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
    let data_start = 10 + header_len;
    if bytes.len() < data_start {
        return Err("npy header truncated".to_string());
    }
    let header_str = std::str::from_utf8(&bytes[10..data_start])
        .map_err(|e| format!("npy header not utf-8: {e}"))?
        .trim();
    let dtype_str = extract_dict_value(header_str, "descr")?;
    let shape_str = extract_dict_value(header_str, "shape")?;
    let shape = parse_shape_tuple(&shape_str)?;
    Ok((dtype_str, shape, data_start))
}
/// Extract a value from a Python-style dict string for a given key.
pub(super) fn extract_dict_value(header: &str, key: &str) -> Result<String, String> {
    let search = format!("'{key}'");
    let pos = header
        .find(&search)
        .ok_or_else(|| format!("key '{key}' not found in npy header"))?;
    let rest = &header[pos + search.len()..];
    let rest = rest.trim_start();
    let rest = rest
        .strip_prefix(':')
        .ok_or("missing ':' after key")?
        .trim_start();
    if rest.starts_with('\'') {
        let inner = rest.strip_prefix('\'').expect("prefix should be present");
        let end = inner.find('\'').ok_or("unterminated string value")?;
        Ok(inner[..end].to_string())
    } else if rest.starts_with('(') {
        let end = rest.find(')').ok_or("unterminated tuple value")? + 1;
        Ok(rest[..end].to_string())
    } else {
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        Ok(rest[..end].trim().to_string())
    }
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
/// Parse a `.npy` byte buffer and return `(shape, f64 data)`.
pub fn read_npy_f64(bytes: &[u8]) -> Result<(Vec<usize>, Vec<f64>), String> {
    let (dtype_str, shape, data_start) = parse_npy_header(bytes)?;
    if dtype_str != "<f8" {
        return Err(format!("expected dtype '<f8', got '{dtype_str}'"));
    }
    let n_elems: usize = shape.iter().product();
    let expected_bytes = data_start + n_elems * 8;
    if bytes.len() < expected_bytes {
        return Err(format!(
            "data truncated: expected {expected_bytes} bytes, got {}",
            bytes.len()
        ));
    }
    let mut data = Vec::with_capacity(n_elems);
    let mut pos = data_start;
    for _ in 0..n_elems {
        let v = f64::from_le_bytes(
            bytes[pos..pos + 8]
                .try_into()
                .expect("slice length must match"),
        );
        pos += 8;
        data.push(v);
    }
    Ok((shape, data))
}
/// Parse a `.npy` byte buffer and return `(shape, f32 data)`.
pub fn read_npy_f32(bytes: &[u8]) -> Result<(Vec<usize>, Vec<f32>), String> {
    let (dtype_str, shape, data_start) = parse_npy_header(bytes)?;
    if dtype_str != "<f4" {
        return Err(format!("expected dtype '<f4', got '{dtype_str}'"));
    }
    let n_elems: usize = shape.iter().product();
    let expected_bytes = data_start + n_elems * 4;
    if bytes.len() < expected_bytes {
        return Err(format!(
            "data truncated: expected {expected_bytes} bytes, got {}",
            bytes.len()
        ));
    }
    let mut data = Vec::with_capacity(n_elems);
    let mut pos = data_start;
    for _ in 0..n_elems {
        let v = f32::from_le_bytes(
            bytes[pos..pos + 4]
                .try_into()
                .expect("slice length must match"),
        );
        pos += 4;
        data.push(v);
    }
    Ok((shape, data))
}
/// Parse a `.npy` byte buffer and return `(shape, i32 data)`.
pub fn read_npy_i32(bytes: &[u8]) -> Result<(Vec<usize>, Vec<i32>), String> {
    let (dtype_str, shape, data_start) = parse_npy_header(bytes)?;
    if dtype_str != "<i4" {
        return Err(format!("expected dtype '<i4', got '{dtype_str}'"));
    }
    let n_elems: usize = shape.iter().product();
    let expected_bytes = data_start + n_elems * 4;
    if bytes.len() < expected_bytes {
        return Err(format!(
            "data truncated: expected {expected_bytes} bytes, got {}",
            bytes.len()
        ));
    }
    let mut data = Vec::with_capacity(n_elems);
    let mut pos = data_start;
    for _ in 0..n_elems {
        let v = i32::from_le_bytes(
            bytes[pos..pos + 4]
                .try_into()
                .expect("slice length must match"),
        );
        pos += 4;
        data.push(v);
    }
    Ok((shape, data))
}
/// Parse a `.npy` byte buffer and return `(shape, i64 data)`.
pub fn read_npy_i64(bytes: &[u8]) -> Result<(Vec<usize>, Vec<i64>), String> {
    let (dtype_str, shape, data_start) = parse_npy_header(bytes)?;
    if dtype_str != "<i8" {
        return Err(format!("expected dtype '<i8', got '{dtype_str}'"));
    }
    let n_elems: usize = shape.iter().product();
    let expected_bytes = data_start + n_elems * 8;
    if bytes.len() < expected_bytes {
        return Err(format!(
            "data truncated: expected {expected_bytes} bytes, got {}",
            bytes.len()
        ));
    }
    let mut data = Vec::with_capacity(n_elems);
    let mut pos = data_start;
    for _ in 0..n_elems {
        let v = i64::from_le_bytes(
            bytes[pos..pos + 8]
                .try_into()
                .expect("slice length must match"),
        );
        pos += 8;
        data.push(v);
    }
    Ok((shape, data))
}
/// Auto-detect dtype from header and return the NpyDtype.
pub fn detect_npy_dtype(bytes: &[u8]) -> Result<NpyDtype, String> {
    let (dtype_str, _, _) = parse_npy_header(bytes)?;
    NpyDtype::from_numpy_str(&dtype_str)
}
/// Auto-detect and return the shape from a .npy byte buffer.
pub fn read_npy_shape(bytes: &[u8]) -> Result<Vec<usize>, String> {
    let (_, shape, _) = parse_npy_header(bytes)?;
    Ok(shape)
}
pub(super) fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32, String> {
    if *pos + 4 > data.len() {
        return Err(format!("unexpected EOF reading u32 at offset {pos}"));
    }
    let v = u32::from_le_bytes(
        data[*pos..*pos + 4]
            .try_into()
            .expect("slice length must match"),
    );
    *pos += 4;
    Ok(v)
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
