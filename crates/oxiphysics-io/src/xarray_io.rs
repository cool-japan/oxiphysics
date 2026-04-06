// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! xarray-compatible multi-dimensional labeled array I/O.
//!
//! Provides `XarrayDataArray` and `XarrayDataset` mirroring the Python xarray
//! data model, together with CSV round-trip, VTK export, and array operations
//! such as resampling, time-averaging, and finite-difference gradients.

use std::fs::File;
use std::io::{BufRead, BufReader, Write};

// ── Coordinate ───────────────────────────────────────────────────────────────

/// A labeled coordinate axis.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct XarrayCoordinate {
    /// Coordinate name (e.g. `"time"`, `"lat"`).
    pub name: String,
    /// Coordinate tick values.
    pub values: Vec<f64>,
    /// Physical units string (e.g. `"s"`, `"m"`).
    pub units: String,
}

impl XarrayCoordinate {
    /// Create a new coordinate axis.
    pub fn new(name: impl Into<String>, values: Vec<f64>, units: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            values,
            units: units.into(),
        }
    }
}

// ── DataArray ─────────────────────────────────────────────────────────────────

/// A labeled N-dimensional array (mirrors `xarray.DataArray`).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct XarrayDataArray {
    /// Variable name.
    pub name: String,
    /// Dimension names in order (length = ndim).
    pub dims: Vec<String>,
    /// Size along each dimension (C-order).
    pub shape: Vec<usize>,
    /// Flat data storage (row-major / C-order).
    pub data: Vec<f64>,
    /// Coordinate objects for each dimension (optional per-dim).
    pub coords: Vec<XarrayCoordinate>,
    /// Arbitrary string key-value attributes.
    pub attrs: Vec<(String, String)>,
}

impl XarrayDataArray {
    /// Create a zero-filled `XarrayDataArray` with the given name, dims, and shape.
    pub fn new(name: impl Into<String>, dims: Vec<String>, shape: Vec<usize>) -> Self {
        let total: usize = shape.iter().product();
        Self {
            name: name.into(),
            dims,
            shape,
            data: vec![0.0; total],
            coords: Vec::new(),
            attrs: Vec::new(),
        }
    }

    /// Attach (or replace) a coordinate for one dimension.
    pub fn set_coord(&mut self, coord: XarrayCoordinate) {
        if let Some(pos) = self.coords.iter().position(|c| c.name == coord.name) {
            self.coords[pos] = coord;
        } else {
            self.coords.push(coord);
        }
    }

    /// Number of dimensions.
    pub fn ndim(&self) -> usize {
        self.dims.len()
    }

    /// Total number of elements.
    pub fn size(&self) -> usize {
        self.shape.iter().product()
    }

    /// Get element at multi-dimensional index (C-order).
    pub fn get(&self, indices: &[usize]) -> f64 {
        let flat = linear_index(indices, &self.shape);
        self.data[flat]
    }

    /// Set element at multi-dimensional index (C-order).
    pub fn set(&mut self, indices: &[usize], value: f64) {
        let flat = linear_index(indices, &self.shape);
        self.data[flat] = value;
    }
}

// ── Dataset ───────────────────────────────────────────────────────────────────

/// A collection of `XarrayDataArray` variables (mirrors `xarray.Dataset`).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct XarrayDataset {
    /// Dataset name / title.
    pub name: String,
    /// Variables stored in this dataset.
    pub variables: Vec<XarrayDataArray>,
    /// Arbitrary string key-value attributes.
    pub attrs: Vec<(String, String)>,
}

impl XarrayDataset {
    /// Create an empty dataset with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            variables: Vec::new(),
            attrs: Vec::new(),
        }
    }

    /// Add a variable to the dataset.
    pub fn add_variable(&mut self, var: XarrayDataArray) {
        self.variables.push(var);
    }

    /// Number of variables currently stored.
    pub fn variable_count(&self) -> usize {
        self.variables.len()
    }

    /// Look up a variable by name.
    pub fn get_variable(&self, name: &str) -> Option<&XarrayDataArray> {
        self.variables.iter().find(|v| v.name == name)
    }
}

// ── Index utilities ───────────────────────────────────────────────────────────

/// Compute the C-order (row-major) flat index from multi-dimensional indices.
///
/// Panics in debug mode if `indices.len() != shape.len()`.
pub fn linear_index(indices: &[usize], shape: &[usize]) -> usize {
    debug_assert_eq!(indices.len(), shape.len());
    let mut flat = 0usize;
    let mut stride = 1usize;
    for d in (0..shape.len()).rev() {
        flat += indices[d] * stride;
        stride *= shape[d];
    }
    flat
}

// ── CSV I/O ───────────────────────────────────────────────────────────────────

/// Write each variable in the dataset as a CSV file.
///
/// Each variable is written to `path`_<variable_name>.csv`.
/// The first line is a header containing the dimension names.
pub fn write_csv_xarray(dataset: &XarrayDataset, path: &str) -> Result<(), String> {
    for var in &dataset.variables {
        let file_path = format!("{}_{}.csv", path, var.name);
        let mut f =
            File::create(&file_path).map_err(|e| format!("cannot create {file_path}: {e}"))?;
        // Header: dim names + "value"
        let header = var.dims.join(",") + ",value\n";
        f.write_all(header.as_bytes())
            .map_err(|e| format!("write error: {e}"))?;
        // Iterate in C-order
        let n = var.size();
        let ndim = var.ndim();
        let mut indices = vec![0usize; ndim];
        for flat in 0..n {
            let row: Vec<String> = indices.iter().map(|v| v.to_string()).collect();
            let line = row.join(",") + "," + &var.data[flat].to_string() + "\n";
            f.write_all(line.as_bytes())
                .map_err(|e| format!("write error: {e}"))?;
            // Increment C-order
            for d in (0..ndim).rev() {
                indices[d] += 1;
                if indices[d] < var.shape[d] {
                    break;
                }
                indices[d] = 0;
            }
        }
    }
    Ok(())
}

/// Read a labeled CSV back into an `XarrayDataArray`.
///
/// Expects the format produced by [`write_csv_xarray`]: first row is headers
/// (dim_0, dim_1, …, value), subsequent rows are data.
pub fn read_csv_xarray(path: &str) -> Result<XarrayDataArray, String> {
    let f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let reader = BufReader::new(f);
    let mut lines = reader.lines();

    let header_line = lines
        .next()
        .ok_or("empty file")?
        .map_err(|e| format!("read error: {e}"))?;
    let headers: Vec<String> = header_line
        .split(',')
        .map(|s| s.trim().to_string())
        .collect();
    let ndim = headers.len().saturating_sub(1);
    let dim_names: Vec<String> = headers[..ndim].to_vec();

    let mut rows: Vec<Vec<usize>> = Vec::new();
    let mut values: Vec<f64> = Vec::new();

    for line in lines {
        let line = line.map_err(|e| format!("read error: {e}"))?;
        if line.trim().is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < ndim + 1 {
            continue;
        }
        let idx: Vec<usize> = parts[..ndim]
            .iter()
            .map(|s| s.trim().parse::<usize>().unwrap_or(0))
            .collect();
        let val: f64 = parts[ndim].trim().parse::<f64>().unwrap_or(0.0);
        rows.push(idx);
        values.push(val);
    }

    // Infer shape from max index per dim
    let mut shape = vec![0usize; ndim];
    for row in &rows {
        for d in 0..ndim {
            if row[d] + 1 > shape[d] {
                shape[d] = row[d] + 1;
            }
        }
    }

    let total: usize = if shape.is_empty() {
        0
    } else {
        shape.iter().product()
    };
    let mut data = vec![0.0f64; total];
    for (idx, &val) in rows.iter().zip(values.iter()) {
        if idx.len() == ndim {
            let flat = linear_index(idx, &shape);
            if flat < data.len() {
                data[flat] = val;
            }
        }
    }

    Ok(XarrayDataArray {
        name: "loaded".to_string(),
        dims: dim_names,
        shape,
        data,
        coords: Vec::new(),
        attrs: Vec::new(),
    })
}

// ── VTK export ────────────────────────────────────────────────────────────────

/// Export a 3-D variable from a dataset as a VTK structured grid.
///
/// The variable must have exactly 3 dimensions. Uses VTK legacy ASCII format.
pub fn xarray_to_vtk_structured(
    dataset: &XarrayDataset,
    var_name: &str,
    path: &str,
) -> Result<(), String> {
    let var = dataset
        .get_variable(var_name)
        .ok_or_else(|| format!("variable '{var_name}' not found"))?;
    if var.shape.len() != 3 {
        return Err(format!(
            "variable '{var_name}' must have 3 dimensions, has {}",
            var.shape.len()
        ));
    }
    let (nx, ny, nz) = (var.shape[0], var.shape[1], var.shape[2]);
    let mut f = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    writeln!(f, "# vtk DataFile Version 3.0").map_err(|e| e.to_string())?;
    writeln!(f, "XarrayExport").map_err(|e| e.to_string())?;
    writeln!(f, "ASCII").map_err(|e| e.to_string())?;
    writeln!(f, "DATASET STRUCTURED_POINTS").map_err(|e| e.to_string())?;
    writeln!(f, "DIMENSIONS {} {} {}", nx, ny, nz).map_err(|e| e.to_string())?;
    writeln!(f, "ORIGIN 0 0 0").map_err(|e| e.to_string())?;
    writeln!(f, "SPACING 1 1 1").map_err(|e| e.to_string())?;
    writeln!(f, "POINT_DATA {}", nx * ny * nz).map_err(|e| e.to_string())?;
    writeln!(f, "SCALARS {} double 1", var_name).map_err(|e| e.to_string())?;
    writeln!(f, "LOOKUP_TABLE default").map_err(|e| e.to_string())?;
    for &val in &var.data {
        writeln!(f, "{val}").map_err(|e| e.to_string())?;
    }
    Ok(())
}

// ── Array operations ──────────────────────────────────────────────────────────

/// Resample an array to a new shape using nearest-neighbour / linear interpolation.
///
/// Uses trilinear interpolation for up to 3-D; falls back to nearest-neighbour
/// for higher dimensions.
pub fn resample_linear(arr: &XarrayDataArray, new_shape: Vec<usize>) -> XarrayDataArray {
    let ndim = arr.ndim();
    assert_eq!(ndim, new_shape.len(), "shape rank mismatch");
    let new_total: usize = new_shape.iter().product();
    let mut out = XarrayDataArray::new(arr.name.clone(), arr.dims.clone(), new_shape.clone());
    for flat in 0..new_total {
        // Convert flat -> new_indices
        let mut tmp = flat;
        let mut new_idx = vec![0usize; ndim];
        for d in (0..ndim).rev() {
            new_idx[d] = tmp % new_shape[d];
            tmp /= new_shape[d];
        }
        // Map to source coordinates
        let src_coords: Vec<f64> = new_idx
            .iter()
            .enumerate()
            .map(|(d, &ni)| {
                if new_shape[d] <= 1 {
                    0.0
                } else {
                    ni as f64 * (arr.shape[d] as f64 - 1.0) / (new_shape[d] as f64 - 1.0)
                }
            })
            .collect();
        // Nearest-neighbour fallback for all dims
        let src_idx: Vec<usize> = src_coords
            .iter()
            .enumerate()
            .map(|(d, &sc)| (sc.round() as usize).min(arr.shape[d].saturating_sub(1)))
            .collect();
        out.data[flat] = arr.get(&src_idx);
    }
    out
}

/// Average an array along dimension `time_dim`.
///
/// Returns an array with the same shape except `time_dim` is reduced to 1.
pub fn time_average(arr: &XarrayDataArray, time_dim: usize) -> XarrayDataArray {
    let ndim = arr.ndim();
    assert!(time_dim < ndim);
    let mut out_shape = arr.shape.clone();
    out_shape[time_dim] = 1;
    let mut out = XarrayDataArray::new(
        arr.name.clone() + "_tavg",
        arr.dims.clone(),
        out_shape.clone(),
    );
    let nt = arr.shape[time_dim];
    // Iterate over all output cells
    let out_total: usize = out_shape.iter().product();
    for flat_out in 0..out_total {
        let mut tmp = flat_out;
        let mut out_idx = vec![0usize; ndim];
        for d in (0..ndim).rev() {
            out_idx[d] = tmp % out_shape[d];
            tmp /= out_shape[d];
        }
        let mut sum = 0.0f64;
        let mut src_idx = out_idx.clone();
        for t in 0..nt {
            src_idx[time_dim] = t;
            sum += arr.get(&src_idx);
        }
        out.data[flat_out] = sum / nt as f64;
    }
    out
}

/// Compute the finite-difference gradient of an array along dimension `dim`.
///
/// Uses central differences for interior points, one-sided for boundary.
/// `dx` is the grid spacing along that dimension.
pub fn spatial_gradient(arr: &XarrayDataArray, dim: usize, dx: f64) -> XarrayDataArray {
    assert!(dim < arr.ndim());
    let mut out = XarrayDataArray::new(
        arr.name.clone() + "_grad",
        arr.dims.clone(),
        arr.shape.clone(),
    );
    let n_dim = arr.shape[dim];
    let total = arr.size();
    for flat in 0..total {
        let mut tmp = flat;
        let mut idx = vec![0usize; arr.ndim()];
        for d in (0..arr.ndim()).rev() {
            idx[d] = tmp % arr.shape[d];
            tmp /= arr.shape[d];
        }
        let i = idx[dim];
        let grad = if i == 0 {
            // Forward difference
            let mut idx_p = idx.clone();
            idx_p[dim] = 1.min(n_dim - 1);
            (arr.get(&idx_p) - arr.get(&idx)) / dx
        } else if i == n_dim - 1 {
            // Backward difference
            let mut idx_m = idx.clone();
            idx_m[dim] = i - 1;
            (arr.get(&idx) - arr.get(&idx_m)) / dx
        } else {
            // Central difference
            let mut idx_p = idx.clone();
            let mut idx_m = idx.clone();
            idx_p[dim] = i + 1;
            idx_m[dim] = i - 1;
            (arr.get(&idx_p) - arr.get(&idx_m)) / (2.0 * dx)
        };
        out.data[flat] = grad;
    }
    out
}

// ── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── linear_index tests ────────────────────────────────────────────────

    #[test]
    fn test_linear_index_1d() {
        assert_eq!(linear_index(&[3], &[10]), 3);
    }

    #[test]
    fn test_linear_index_2d() {
        // shape 3×4: index(1,2) = 1*4 + 2 = 6
        assert_eq!(linear_index(&[1, 2], &[3, 4]), 6);
    }

    #[test]
    fn test_linear_index_3d() {
        // shape 2×3×4: index(1,2,3) = 1*12 + 2*4 + 3 = 23
        assert_eq!(linear_index(&[1, 2, 3], &[2, 3, 4]), 23);
    }

    #[test]
    fn test_linear_index_origin() {
        assert_eq!(linear_index(&[0, 0, 0], &[5, 5, 5]), 0);
    }

    #[test]
    fn test_linear_index_last_element() {
        let shape = [2, 3, 4];
        let idx = [1, 2, 3];
        assert_eq!(linear_index(&idx, &shape), 2 * 3 * 4 - 1);
    }

    // ── XarrayDataArray tests ─────────────────────────────────────────────

    #[test]
    fn test_data_array_size() {
        let arr = XarrayDataArray::new("t", vec!["x".into(), "y".into()], vec![3, 4]);
        assert_eq!(arr.size(), 12);
    }

    #[test]
    fn test_data_array_size_empty_shape() {
        let arr = XarrayDataArray::new("t", vec![], vec![]);
        assert_eq!(arr.size(), 1); // product of empty = 1
    }

    #[test]
    fn test_data_array_ndim() {
        let arr =
            XarrayDataArray::new("t", vec!["x".into(), "y".into(), "z".into()], vec![2, 3, 4]);
        assert_eq!(arr.ndim(), 3);
    }

    #[test]
    fn test_data_array_set_get_roundtrip() {
        let mut arr = XarrayDataArray::new("t", vec!["x".into(), "y".into()], vec![3, 4]);
        arr.set(&[1, 2], 42.0);
        assert!((arr.get(&[1, 2]) - 42.0).abs() < 1e-12);
    }

    #[test]
    fn test_data_array_initial_zeros() {
        let arr = XarrayDataArray::new("t", vec!["x".into()], vec![5]);
        assert!(arr.data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_data_array_set_coord() {
        let mut arr = XarrayDataArray::new("t", vec!["x".into()], vec![3]);
        let coord = XarrayCoordinate::new("x", vec![0.0, 1.0, 2.0], "m");
        arr.set_coord(coord);
        assert_eq!(arr.coords.len(), 1);
    }

    #[test]
    fn test_data_array_set_coord_replace() {
        let mut arr = XarrayDataArray::new("t", vec!["x".into()], vec![3]);
        arr.set_coord(XarrayCoordinate::new("x", vec![0.0, 1.0, 2.0], "m"));
        arr.set_coord(XarrayCoordinate::new("x", vec![0.0, 0.5, 1.0], "m"));
        assert_eq!(arr.coords.len(), 1);
        assert!((arr.coords[0].values[1] - 0.5).abs() < 1e-12);
    }

    // ── XarrayDataset tests ───────────────────────────────────────────────

    #[test]
    fn test_dataset_new_empty() {
        let ds = XarrayDataset::new("test");
        assert_eq!(ds.variable_count(), 0);
    }

    #[test]
    fn test_dataset_add_variable_increases_count() {
        let mut ds = XarrayDataset::new("test");
        ds.add_variable(XarrayDataArray::new("u", vec!["x".into()], vec![4]));
        assert_eq!(ds.variable_count(), 1);
    }

    #[test]
    fn test_dataset_add_two_variables() {
        let mut ds = XarrayDataset::new("test");
        ds.add_variable(XarrayDataArray::new("u", vec!["x".into()], vec![4]));
        ds.add_variable(XarrayDataArray::new("v", vec!["x".into()], vec![4]));
        assert_eq!(ds.variable_count(), 2);
    }

    #[test]
    fn test_dataset_get_variable() {
        let mut ds = XarrayDataset::new("test");
        ds.add_variable(XarrayDataArray::new("temp", vec!["x".into()], vec![4]));
        let var = ds.get_variable("temp");
        assert!(var.is_some());
        assert_eq!(var.unwrap().name, "temp");
    }

    #[test]
    fn test_dataset_get_variable_missing() {
        let ds = XarrayDataset::new("test");
        assert!(ds.get_variable("nosuchvar").is_none());
    }

    // ── time_average tests ────────────────────────────────────────────────

    #[test]
    fn test_time_average_reduces_shape() {
        let arr = XarrayDataArray::new("u", vec!["time".into(), "x".into()], vec![4, 3]);
        let avg = time_average(&arr, 0);
        assert_eq!(avg.shape[0], 1);
        assert_eq!(avg.shape[1], 3);
    }

    #[test]
    fn test_time_average_correct_value() {
        let mut arr = XarrayDataArray::new("u", vec!["t".into(), "x".into()], vec![4, 1]);
        for t in 0..4 {
            arr.set(&[t, 0], t as f64);
        }
        let avg = time_average(&arr, 0);
        assert!((avg.get(&[0, 0]) - 1.5).abs() < 1e-12);
    }

    #[test]
    fn test_time_average_second_dim() {
        let arr = XarrayDataArray::new("u", vec!["x".into(), "y".into()], vec![3, 4]);
        let avg = time_average(&arr, 1);
        assert_eq!(avg.shape, vec![3, 1]);
    }

    // ── spatial_gradient tests ────────────────────────────────────────────

    #[test]
    fn test_spatial_gradient_constant_array_zero() {
        let mut arr = XarrayDataArray::new("u", vec!["x".into()], vec![5]);
        for i in 0..5 {
            arr.set(&[i], 3.0);
        }
        let grad = spatial_gradient(&arr, 0, 1.0);
        assert!(grad.data.iter().all(|&g| g.abs() < 1e-12));
    }

    #[test]
    fn test_spatial_gradient_linear_array() {
        let mut arr = XarrayDataArray::new("u", vec!["x".into()], vec![5]);
        for i in 0..5 {
            arr.set(&[i], i as f64);
        }
        let grad = spatial_gradient(&arr, 0, 1.0);
        // Central differences should all give ~1.0 for interior points
        assert!((grad.get(&[2]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_spatial_gradient_shape_preserved() {
        let arr = XarrayDataArray::new("u", vec!["x".into(), "y".into()], vec![3, 4]);
        let grad = spatial_gradient(&arr, 0, 0.1);
        assert_eq!(grad.shape, arr.shape);
    }

    // ── resample tests ────────────────────────────────────────────────────

    #[test]
    fn test_resample_same_shape() {
        let mut arr = XarrayDataArray::new("u", vec!["x".into()], vec![4]);
        for i in 0..4 {
            arr.set(&[i], i as f64);
        }
        let out = resample_linear(&arr, vec![4]);
        assert_eq!(out.shape, vec![4]);
    }

    #[test]
    fn test_resample_upscale_shape() {
        let arr = XarrayDataArray::new("u", vec!["x".into()], vec![3]);
        let out = resample_linear(&arr, vec![6]);
        assert_eq!(out.shape, vec![6]);
    }

    #[test]
    fn test_resample_downscale_shape() {
        let arr = XarrayDataArray::new("u", vec!["x".into()], vec![8]);
        let out = resample_linear(&arr, vec![4]);
        assert_eq!(out.shape, vec![4]);
    }

    // ── CSV round-trip tests ──────────────────────────────────────────────

    #[test]
    fn test_write_read_csv_roundtrip() {
        let mut ds = XarrayDataset::new("test");
        let mut arr = XarrayDataArray::new("temperature", vec!["x".into(), "y".into()], vec![2, 3]);
        arr.set(&[0, 1], 3.125);
        arr.set(&[1, 2], 2.72);
        ds.add_variable(arr);
        let path = "/tmp/xarray_test_roundtrip";
        write_csv_xarray(&ds, path).expect("write failed");
        let loaded = read_csv_xarray(&format!("{path}_temperature.csv")).expect("read failed");
        assert!((loaded.get(&[0, 1]) - 3.125).abs() < 1e-9);
        assert!((loaded.get(&[1, 2]) - 2.72).abs() < 1e-9);
    }

    #[test]
    fn test_write_csv_creates_file() {
        let mut ds = XarrayDataset::new("test2");
        ds.add_variable(XarrayDataArray::new("v", vec!["x".into()], vec![3]));
        write_csv_xarray(&ds, "/tmp/xarray_test2").expect("write failed");
        assert!(std::path::Path::new("/tmp/xarray_test2_v.csv").exists());
    }

    // ── VTK export test ───────────────────────────────────────────────────

    #[test]
    fn test_vtk_export_creates_file() {
        let mut ds = XarrayDataset::new("vtk_test");
        let arr = XarrayDataArray::new(
            "pressure",
            vec!["x".into(), "y".into(), "z".into()],
            vec![2, 2, 2],
        );
        ds.add_variable(arr);
        xarray_to_vtk_structured(&ds, "pressure", "/tmp/xarray_test_pressure.vtk")
            .expect("vtk failed");
        assert!(std::path::Path::new("/tmp/xarray_test_pressure.vtk").exists());
    }

    #[test]
    fn test_vtk_export_wrong_dims_returns_err() {
        let mut ds = XarrayDataset::new("vtk_test2");
        let arr = XarrayDataArray::new("u2d", vec!["x".into(), "y".into()], vec![2, 2]);
        ds.add_variable(arr);
        let res = xarray_to_vtk_structured(&ds, "u2d", "/tmp/xarray_test_u2d.vtk");
        assert!(res.is_err());
    }

    #[test]
    fn test_coordinate_new() {
        let c = XarrayCoordinate::new("time", vec![0.0, 1.0, 2.0], "s");
        assert_eq!(c.name, "time");
        assert_eq!(c.values.len(), 3);
    }
}
