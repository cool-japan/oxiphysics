//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::type_complexity)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{NPY_MAGIC, NPY_MAJOR, NPY_MINOR};
#[allow(unused_imports)]
use super::functions_2::*;

/// A non-owning view into a contiguous slice of `f64` data with shape metadata.
///
/// Useful for reading a row, column, or arbitrary slab from a multi-dimensional
/// array without copying data.
#[allow(dead_code)]
pub struct NpySlice<'a> {
    /// Underlying data slice.
    pub data: &'a [f64],
    /// Shape of this view.
    pub shape: Vec<usize>,
}
#[allow(dead_code)]
impl<'a> NpySlice<'a> {
    /// Create a new view.
    pub fn new(data: &'a [f64], shape: Vec<usize>) -> std::result::Result<Self, String> {
        let expected: usize = shape.iter().product();
        if expected != data.len() {
            return Err(format!(
                "NpySlice: data length {} != shape product {}",
                data.len(),
                expected
            ));
        }
        Ok(NpySlice { data, shape })
    }
    /// Number of dimensions.
    pub fn ndim(&self) -> usize {
        self.shape.len()
    }
    /// Total number of elements.
    pub fn numel(&self) -> usize {
        self.shape.iter().product()
    }
    /// Extract a single row (axis 0) from a 2-D slice.
    pub fn row(&self, row_idx: usize) -> std::result::Result<&[f64], String> {
        if self.shape.len() != 2 {
            return Err(format!(
                "row() requires 2-D slice, got {}D",
                self.shape.len()
            ));
        }
        let ncols = self.shape[1];
        if row_idx >= self.shape[0] {
            return Err(format!(
                "row {} out of bounds (shape[0]={})",
                row_idx, self.shape[0]
            ));
        }
        Ok(&self.data[row_idx * ncols..(row_idx + 1) * ncols])
    }
    /// Get element at multi-dimensional index (row-major).
    pub fn get(&self, indices: &[usize]) -> std::result::Result<f64, String> {
        let flat = flat_index(indices, &self.shape)?;
        Ok(self.data[flat])
    }
}
/// A masked NumPy-style array: elements where `mask[i]` is `true` are considered
/// invalid/missing (following NumPy `ma` conventions).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NpyMaskedArray {
    /// Underlying data.
    pub data: Vec<f64>,
    /// Per-element mask. `true` = masked (invalid).
    pub mask: Vec<bool>,
    /// Fill value used when accessing masked elements.
    pub fill_value: f64,
    /// Shape of the array.
    pub shape: Vec<usize>,
}
#[allow(dead_code)]
impl NpyMaskedArray {
    /// Create a masked array from data and mask.
    pub fn new(
        data: Vec<f64>,
        mask: Vec<bool>,
        shape: Vec<usize>,
        fill_value: f64,
    ) -> std::result::Result<Self, String> {
        let n: usize = shape.iter().product();
        if data.len() != n {
            return Err(format!("data length {} != shape product {}", data.len(), n));
        }
        if mask.len() != n {
            return Err(format!("mask length {} != shape product {}", mask.len(), n));
        }
        Ok(Self {
            data,
            mask,
            fill_value,
            shape,
        })
    }
    /// Create with all elements unmasked.
    pub fn from_data(data: Vec<f64>, shape: Vec<usize>) -> std::result::Result<Self, String> {
        let n: usize = shape.iter().product();
        if data.len() != n {
            return Err(format!("data length {} != shape product {}", data.len(), n));
        }
        let mask = vec![false; n];
        Ok(Self {
            data,
            mask,
            fill_value: 1e20,
            shape,
        })
    }
    /// Get element, returning `fill_value` if masked.
    pub fn get_filled(&self, idx: usize) -> f64 {
        if self.mask[idx] {
            self.fill_value
        } else {
            self.data[idx]
        }
    }
    /// Number of valid (unmasked) elements.
    pub fn count_valid(&self) -> usize {
        self.mask.iter().filter(|&&m| !m).count()
    }
    /// Mean of valid (unmasked) elements. Returns `None` if all masked.
    pub fn mean_valid(&self) -> Option<f64> {
        let (sum, count) = self
            .data
            .iter()
            .zip(self.mask.iter())
            .filter(|&(_, &m)| !m)
            .fold((0.0_f64, 0_usize), |(s, c), (&v, _)| (s + v, c + 1));
        if count == 0 {
            None
        } else {
            Some(sum / count as f64)
        }
    }
    /// Fill masked values with `fill_value` and return a plain `Vec`f64`.
    pub fn filled(&self) -> Vec<f64> {
        self.data
            .iter()
            .zip(self.mask.iter())
            .map(|(&v, &m)| if m { self.fill_value } else { v })
            .collect()
    }
    /// Apply a threshold mask: mask elements where `|data\[i\]| > threshold`.
    pub fn mask_greater_than(&mut self, threshold: f64) {
        for (m, &v) in self.mask.iter_mut().zip(self.data.iter()) {
            if v.abs() > threshold {
                *m = true;
            }
        }
    }
    /// Unmask all elements.
    pub fn unmask_all(&mut self) {
        self.mask.iter_mut().for_each(|m| *m = false);
    }
}
/// A simple structured / record array: multiple named columns each stored
/// as a flat `Vec`f64` (all fields promoted to f64 for simplicity).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NpyRecordArray {
    /// Field definitions.
    pub fields: Vec<NpyField>,
    /// Flat column data, one `Vec`f64` per field (length = n_records * field.count).
    pub columns: Vec<Vec<f64>>,
    /// Number of records (rows).
    pub n_records: usize,
}
#[allow(dead_code)]
impl NpyRecordArray {
    /// Create an empty record array with given field schema.
    pub fn new(fields: Vec<NpyField>) -> Self {
        let columns = vec![Vec::new(); fields.len()];
        Self {
            fields,
            columns,
            n_records: 0,
        }
    }
    /// Push one record; `values` must have one entry per field (scalars) or
    /// `field.count` entries for vector fields.
    pub fn push_record(&mut self, values: &[f64]) -> std::result::Result<(), String> {
        let total: usize = self.fields.iter().map(|f| f.count).sum();
        if values.len() != total {
            return Err(format!(
                "push_record: expected {total} values, got {}",
                values.len()
            ));
        }
        let mut offset = 0;
        for (col, field) in self.columns.iter_mut().zip(self.fields.iter()) {
            col.extend_from_slice(&values[offset..offset + field.count]);
            offset += field.count;
        }
        self.n_records += 1;
        Ok(())
    }
    /// Get a column by field name.
    pub fn column(&self, name: &str) -> Option<&[f64]> {
        self.fields
            .iter()
            .position(|f| f.name == name)
            .map(|i| self.columns[i].as_slice())
    }
    /// Get a single scalar value: `(record_idx, field_name)`.
    pub fn get_scalar(&self, record: usize, name: &str) -> std::result::Result<f64, String> {
        let fi = self
            .fields
            .iter()
            .position(|f| f.name == name)
            .ok_or_else(|| format!("field '{name}' not found"))?;
        let field = &self.fields[fi];
        if field.count != 1 {
            return Err(format!(
                "field '{name}' is not scalar (count={})",
                field.count
            ));
        }
        if record >= self.n_records {
            return Err(format!(
                "record {record} out of range (n_records={})",
                self.n_records
            ));
        }
        Ok(self.columns[fi][record])
    }
}
/// Element type stored in a `.npy` file.
#[derive(Debug, Clone, PartialEq)]
pub enum NpyDtype {
    /// 64-bit IEEE 754 float (little-endian `<f8`).
    Float64,
    /// 32-bit IEEE 754 float (little-endian `<f4`).
    Float32,
    /// 32-bit signed integer (little-endian `<i4`).
    Int32,
    /// 64-bit signed integer (little-endian `<i8`).
    Int64,
    /// Boolean (`?`).
    Bool,
    /// Unsigned 8-bit integer (`|u1`).
    Uint8,
}
impl NpyDtype {
    /// Returns the NumPy dtype string, e.g. `"<f8"` for `Float64`.
    pub fn numpy_str(&self) -> &str {
        match self {
            NpyDtype::Float64 => "<f8",
            NpyDtype::Float32 => "<f4",
            NpyDtype::Int32 => "<i4",
            NpyDtype::Int64 => "<i8",
            NpyDtype::Bool => "?",
            NpyDtype::Uint8 => "|u1",
        }
    }
    /// Number of bytes per element.
    pub fn element_size(&self) -> usize {
        match self {
            NpyDtype::Float64 => 8,
            NpyDtype::Float32 => 4,
            NpyDtype::Int32 => 4,
            NpyDtype::Int64 => 8,
            NpyDtype::Bool => 1,
            NpyDtype::Uint8 => 1,
        }
    }
    /// Parse a NumPy dtype string into an `NpyDtype`.
    pub fn from_numpy_str(s: &str) -> Result<Self, String> {
        match s {
            "<f8" => Ok(NpyDtype::Float64),
            "<f4" => Ok(NpyDtype::Float32),
            "<i4" => Ok(NpyDtype::Int32),
            "<i8" => Ok(NpyDtype::Int64),
            "?" => Ok(NpyDtype::Bool),
            "|u1" => Ok(NpyDtype::Uint8),
            _ => Err(format!("unsupported dtype: '{s}'")),
        }
    }
}
/// An in-memory NumPy array.
///
/// Only one of `data_f64`, `data_f32`, or `data_i32` is populated,
/// depending on `dtype`.
#[derive(Debug, Clone)]
pub struct NpyArray {
    /// Element type.
    pub dtype: NpyDtype,
    /// Shape (row-major), e.g. `\[3, 2\]` for a 3×2 matrix.
    pub shape: Vec<usize>,
    /// Float64 payload (populated when `dtype == Float64`).
    pub data_f64: Vec<f64>,
    /// Float32 payload (populated when `dtype == Float32`).
    pub data_f32: Vec<f32>,
    /// Int32 payload (populated when `dtype == Int32`).
    pub data_i32: Vec<i32>,
}
impl NpyArray {
    /// Total number of elements (product of shape dimensions).
    pub fn numel(&self) -> usize {
        self.shape.iter().product()
    }
    /// Number of dimensions (length of shape).
    pub fn ndim(&self) -> usize {
        self.shape.len()
    }
    /// Validate that the data length matches the shape.
    pub fn validate(&self) -> Result<(), String> {
        let expected = self.numel();
        let actual = match self.dtype {
            NpyDtype::Float64 => self.data_f64.len(),
            NpyDtype::Float32 => self.data_f32.len(),
            NpyDtype::Int32 => self.data_i32.len(),
            _ => expected,
        };
        if actual != expected {
            Err(format!(
                "shape {:?} expects {} elements, but data has {}",
                self.shape, expected, actual
            ))
        } else {
            Ok(())
        }
    }
    /// Create an NpyArray from f64 data and shape.
    pub fn from_f64(shape: Vec<usize>, data: Vec<f64>) -> Self {
        Self {
            dtype: NpyDtype::Float64,
            shape,
            data_f64: data,
            data_f32: Vec::new(),
            data_i32: Vec::new(),
        }
    }
    /// Create an NpyArray from f32 data and shape.
    pub fn from_f32(shape: Vec<usize>, data: Vec<f32>) -> Self {
        Self {
            dtype: NpyDtype::Float32,
            shape,
            data_f64: Vec::new(),
            data_f32: data,
            data_i32: Vec::new(),
        }
    }
    /// Create an NpyArray from i32 data and shape.
    pub fn from_i32(shape: Vec<usize>, data: Vec<i32>) -> Self {
        Self {
            dtype: NpyDtype::Int32,
            shape,
            data_f64: Vec::new(),
            data_f32: Vec::new(),
            data_i32: data,
        }
    }
    /// Reshape the array (does not change data, just shape metadata).
    pub fn reshape(&mut self, new_shape: Vec<usize>) -> Result<(), String> {
        let old_numel = self.numel();
        let new_numel: usize = new_shape.iter().product();
        if old_numel != new_numel {
            return Err(format!(
                "cannot reshape: old numel={old_numel}, new numel={new_numel}"
            ));
        }
        self.shape = new_shape;
        Ok(())
    }
}
impl NpyArray {
    /// Serialize a structured array with named fields to NPY v1.0 bytes.
    ///
    /// `fields` is a slice of `(name, dtype_str)` pairs, e.g.
    /// `&\[("x", "<f8"), ("y", "<f8")\]`.  The data must be already packed
    /// in record order (all fields of record 0, then record 1, …).
    ///
    /// Returns the NPY bytes suitable for writing to a `.npy` file.
    #[allow(dead_code)]
    pub fn save_structured(
        fields: &[(&str, &str)],
        n_records: usize,
        data_bytes: &[u8],
    ) -> std::result::Result<Vec<u8>, String> {
        if fields.is_empty() {
            return Err("save_structured: field list is empty".into());
        }
        let dtype_parts: Vec<String> = fields
            .iter()
            .map(|(name, dt)| format!("('{}', '{}')", name, dt))
            .collect();
        let dtype_str = format!("[{}]", dtype_parts.join(", "));
        let header_dict = format!(
            "{{'descr': {}, 'fortran_order': False, 'shape': ({},), }}",
            dtype_str, n_records
        );
        let raw_len = header_dict.len() + 1;
        let pad_to = raw_len.div_ceil(64) * 64;
        let padding = pad_to - raw_len;
        let mut header_bytes = header_dict.into_bytes();
        header_bytes.extend(std::iter::repeat_n(b' ', padding));
        header_bytes.push(b'\n');
        let header_len = header_bytes.len() as u16;
        let mut out = Vec::new();
        out.extend_from_slice(NPY_MAGIC);
        out.push(NPY_MAJOR);
        out.push(NPY_MINOR);
        out.extend_from_slice(&header_len.to_le_bytes());
        out.extend_from_slice(&header_bytes);
        out.extend_from_slice(data_bytes);
        Ok(out)
    }
}
/// A higher-level NPZ archive that stores typed `NpyArray` objects by name.
///
/// Unlike `NpzWriter`, this type stores fully parsed `NpyArray` values and
/// supports reading them back without specifying the dtype at call time.
#[allow(dead_code)]
#[derive(Debug, Clone, Default)]
pub struct NpzArchive {
    /// Named arrays stored in the archive.
    pub arrays: Vec<(String, NpyArray)>,
}
#[allow(dead_code)]
impl NpzArchive {
    /// Create an empty archive.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add an `NpyArray` under the given name.
    pub fn insert(&mut self, name: &str, array: NpyArray) {
        self.arrays.push((name.to_string(), array));
    }
    /// Retrieve a reference to the array with the given name.
    pub fn get(&self, name: &str) -> Option<&NpyArray> {
        self.arrays.iter().find(|(n, _)| n == name).map(|(_, a)| a)
    }
    /// List all array names.
    pub fn names(&self) -> Vec<&str> {
        self.arrays.iter().map(|(n, _)| n.as_str()).collect()
    }
    /// Remove an array by name; returns `true` if it was present.
    pub fn remove(&mut self, name: &str) -> bool {
        let before = self.arrays.len();
        self.arrays.retain(|(n, _)| n != name);
        self.arrays.len() < before
    }
    /// Number of arrays stored.
    pub fn len(&self) -> usize {
        self.arrays.len()
    }
    /// Whether the archive is empty.
    pub fn is_empty(&self) -> bool {
        self.arrays.is_empty()
    }
    /// Serialize to the NPZ container format (compatible with `NpzWriter`).
    pub fn to_bytes(&self) -> std::result::Result<Vec<u8>, String> {
        let mut writer = NpzWriter::new();
        for (name, array) in &self.arrays {
            match array.dtype {
                NpyDtype::Float64 => {
                    writer.add_array_f64(name, &array.shape, &array.data_f64);
                }
                NpyDtype::Float32 => {
                    writer.add_array_f32(name, &array.shape, &array.data_f32);
                }
                NpyDtype::Int32 => {
                    writer.add_array_i32(name, &array.shape, &array.data_i32);
                }
                _ => {
                    return Err(format!(
                        "NpzArchive::to_bytes: unsupported dtype {:?}",
                        array.dtype
                    ));
                }
            }
        }
        Ok(writer.to_bytes())
    }
    /// Deserialize from NPZ container bytes.
    pub fn from_bytes(data: &[u8]) -> std::result::Result<Self, String> {
        let writer = NpzWriter::from_bytes(data)?;
        let mut archive = NpzArchive::new();
        for (name, npy_bytes) in &writer.files {
            let dtype = detect_npy_dtype(npy_bytes)?;
            let array = match dtype {
                NpyDtype::Float64 => {
                    let (shape, data_f64) = read_npy_f64(npy_bytes)?;
                    NpyArray::from_f64(shape, data_f64)
                }
                NpyDtype::Float32 => {
                    let (shape, data_f32) = read_npy_f32(npy_bytes)?;
                    NpyArray::from_f32(shape, data_f32)
                }
                NpyDtype::Int32 => {
                    let (shape, data_i32) = read_npy_i32(npy_bytes)?;
                    NpyArray::from_i32(shape, data_i32)
                }
                other => {
                    return Err(format!(
                        "NpzArchive::from_bytes: unsupported dtype {:?} in '{name}'",
                        other
                    ));
                }
            };
            archive.insert(name, array);
        }
        Ok(archive)
    }
}
impl NpzArchive {
    /// Add a pre-built `NpyArray` under `name`, replacing any existing entry
    /// with that name.
    #[allow(dead_code)]
    pub fn add_array(&mut self, name: &str, array: NpyArray) {
        self.arrays.retain(|(n, _)| n.as_str() != name);
        self.arrays.push((name.to_string(), array));
    }
    /// Load all arrays from raw NPZ bytes and return them as a new archive.
    ///
    /// This is an alias for [`NpzArchive::from_bytes`] with a more descriptive
    /// name to match the "load_all" specification.
    #[allow(dead_code)]
    pub fn load_all(data: &[u8]) -> std::result::Result<Self, String> {
        Self::from_bytes(data)
    }
    /// Return an iterator over `(name, &NpyArray)` pairs.
    #[allow(dead_code)]
    pub fn iter(&self) -> impl Iterator<Item = (&str, &NpyArray)> {
        self.arrays.iter().map(|(n, a)| (n.as_str(), a))
    }
    /// Merge another archive into `self`, overwriting duplicate names.
    #[allow(dead_code)]
    pub fn merge(&mut self, other: NpzArchive) {
        for (name, array) in other.arrays {
            self.add_array(&name, array);
        }
    }
    /// Return the total number of elements across all stored arrays.
    #[allow(dead_code)]
    pub fn total_elements(&self) -> usize {
        self.arrays.iter().map(|(_, a)| a.numel()).sum()
    }
}
/// Writer for the simplified `.npz` container format.
///
/// Each array is stored as a raw `.npy` blob.  The container layout is:
/// ```text
/// [count: u32 LE]
/// foreach:
///   [name_len: u32 LE][name UTF-8 bytes]
///   [npy_len: u32 LE][npy bytes]
/// ```
#[derive(Debug, Clone)]
pub struct NpzWriter {
    /// Stored `(name, npy_bytes)` pairs.
    pub files: Vec<(String, Vec<u8>)>,
}
impl NpzWriter {
    /// Create an empty [`NpzWriter`].
    pub fn new() -> Self {
        NpzWriter { files: Vec::new() }
    }
    /// Append a `f64` array under `name`.
    pub fn add_array_f64(&mut self, name: &str, shape: &[usize], data: &[f64]) {
        let npy = write_npy_f64(shape, data);
        self.files.push((name.to_string(), npy));
    }
    /// Append an `f32` array under `name`.
    pub fn add_array_f32(&mut self, name: &str, shape: &[usize], data: &[f32]) {
        let npy = write_npy_f32(shape, data);
        self.files.push((name.to_string(), npy));
    }
    /// Append an `i32` array under `name`.
    pub fn add_array_i32(&mut self, name: &str, shape: &[usize], data: &[i32]) {
        let npy = write_npy_i32(shape, data);
        self.files.push((name.to_string(), npy));
    }
    /// Append an `i64` array under `name`.
    pub fn add_array_i64(&mut self, name: &str, shape: &[usize], data: &[i64]) {
        let npy = write_npy_i64(shape, data);
        self.files.push((name.to_string(), npy));
    }
    /// Number of arrays stored.
    pub fn len(&self) -> usize {
        self.files.len()
    }
    /// Whether the archive is empty.
    pub fn is_empty(&self) -> bool {
        self.files.is_empty()
    }
    /// List all array names.
    pub fn names(&self) -> Vec<&str> {
        self.files.iter().map(|(n, _)| n.as_str()).collect()
    }
    /// Check if an array with the given name exists.
    pub fn contains(&self, name: &str) -> bool {
        self.files.iter().any(|(n, _)| n == name)
    }
    /// Remove an array by name. Returns true if found and removed.
    pub fn remove(&mut self, name: &str) -> bool {
        let before = self.files.len();
        self.files.retain(|(n, _)| n != name);
        self.files.len() < before
    }
    /// Serialize all stored arrays to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(&(self.files.len() as u32).to_le_bytes());
        for (name, npy) in &self.files {
            let name_bytes = name.as_bytes();
            out.extend_from_slice(&(name_bytes.len() as u32).to_le_bytes());
            out.extend_from_slice(name_bytes);
            out.extend_from_slice(&(npy.len() as u32).to_le_bytes());
            out.extend_from_slice(npy);
        }
        out
    }
    /// Deserialize from bytes produced by [`NpzWriter::to_bytes`].
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        let mut pos = 0usize;
        let count = read_u32(data, &mut pos)? as usize;
        let mut files = Vec::with_capacity(count);
        for _ in 0..count {
            let name_len = read_u32(data, &mut pos)? as usize;
            if pos + name_len > data.len() {
                return Err("name out of bounds".to_string());
            }
            let name = std::str::from_utf8(&data[pos..pos + name_len])
                .map_err(|e| format!("invalid UTF-8 in name: {e}"))?
                .to_string();
            pos += name_len;
            let npy_len = read_u32(data, &mut pos)? as usize;
            if pos + npy_len > data.len() {
                return Err("npy payload out of bounds".to_string());
            }
            let npy = data[pos..pos + npy_len].to_vec();
            pos += npy_len;
            files.push((name, npy));
        }
        Ok(NpzWriter { files })
    }
    /// Retrieve a `f64` array by name, returning `(shape, data)`.
    pub fn get_f64(&self, name: &str) -> Option<Result<(Vec<usize>, Vec<f64>), String>> {
        self.files
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, npy)| read_npy_f64(npy))
    }
    /// Retrieve an `f32` array by name, returning `(shape, data)`.
    pub fn get_f32(&self, name: &str) -> Option<Result<(Vec<usize>, Vec<f32>), String>> {
        self.files
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, npy)| read_npy_f32(npy))
    }
    /// Retrieve an `i32` array by name, returning `(shape, data)`.
    pub fn get_i32(&self, name: &str) -> Option<Result<(Vec<usize>, Vec<i32>), String>> {
        self.files
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, npy)| read_npy_i32(npy))
    }
    /// Retrieve an `i64` array by name, returning `(shape, data)`.
    pub fn get_i64(&self, name: &str) -> Option<Result<(Vec<usize>, Vec<i64>), String>> {
        self.files
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, npy)| read_npy_i64(npy))
    }
}
/// A field definition in a structured / record array.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NpyField {
    /// Field name (column name).
    pub name: String,
    /// Data type for this field.
    pub dtype: NpyDtype,
    /// Number of elements per record (1 for scalar, >1 for vector fields).
    pub count: usize,
}
#[allow(dead_code)]
impl NpyField {
    /// Create a scalar field definition.
    pub fn scalar(name: &str, dtype: NpyDtype) -> Self {
        Self {
            name: name.to_string(),
            dtype,
            count: 1,
        }
    }
    /// Create a vector field definition.
    pub fn vector(name: &str, dtype: NpyDtype, count: usize) -> Self {
        Self {
            name: name.to_string(),
            dtype,
            count,
        }
    }
    /// Bytes per record entry for this field.
    pub fn byte_size(&self) -> usize {
        self.dtype.element_size() * self.count
    }
}
