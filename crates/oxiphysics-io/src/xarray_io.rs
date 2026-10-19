// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! xarray-compatible multi-dimensional labeled array I/O.
//!
//! Provides `XarrayDataArray` and `XarrayDataset` mirroring the Python xarray
//! data model, together with CSV round-trip, VTK export, and array operations
//! such as resampling, time-averaging, and finite-difference gradients.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

// ── Coordinate ───────────────────────────────────────────────────────────────

/// A labeled coordinate axis.
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

/// Geometry of one axis of an exported VTK grid.
enum VtkAxis {
    /// Evenly spaced points `origin + i * spacing`. Also the default index
    /// axis (`origin = 0`, `spacing = 1`) of a dimension without a coordinate.
    Uniform { origin: f64, spacing: f64 },
    /// Strictly monotonic points that are not evenly spaced.
    Explicit(Vec<f64>),
}

impl VtkAxis {
    /// Coordinate of each of the `n` points along the axis.
    fn points(&self, n: usize) -> Vec<f64> {
        match self {
            VtkAxis::Uniform { origin, spacing } => {
                (0..n).map(|i| origin + spacing * i as f64).collect()
            }
            VtkAxis::Explicit(values) => values.clone(),
        }
    }
}

/// Resolve the grid axis for dimension `dim` of `var`.
///
/// The coordinate belonging to the dimension is the one whose name equals the
/// dimension name. Without a coordinate the axis is the index axis
/// (`0, 1, 2, ...`). A coordinate must have exactly one value per point of its
/// dimension and be strictly monotonic (increasing or decreasing).
fn vtk_axis(var: &XarrayDataArray, dim: usize) -> Result<VtkAxis, String> {
    let index_axis = VtkAxis::Uniform {
        origin: 0.0,
        spacing: 1.0,
    };
    let n = var.shape[dim];
    let Some(dim_name) = var.dims.get(dim) else {
        return Ok(index_axis);
    };
    let Some(coord) = var.coords.iter().find(|c| &c.name == dim_name) else {
        return Ok(index_axis);
    };
    let values = &coord.values;
    if values.len() != n {
        return Err(format!(
            "coordinate '{}' has {} values but dimension '{dim_name}' has size {n}",
            coord.name,
            values.len()
        ));
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(format!(
            "coordinate '{}' contains non-finite values",
            coord.name
        ));
    }
    let (Some(&first), Some(&last)) = (values.first(), values.last()) else {
        return Ok(index_axis);
    };
    if n == 1 {
        // A single point has no spacing: keep its position, use unit spacing.
        return Ok(VtkAxis::Uniform {
            origin: first,
            spacing: 1.0,
        });
    }
    let step = (last - first) / (n - 1) as f64;
    let increasing = step > 0.0;
    let monotonic = step != 0.0
        && values
            .windows(2)
            .all(|w| if increasing { w[1] > w[0] } else { w[1] < w[0] });
    if !monotonic {
        return Err(format!(
            "coordinate '{}' is not strictly monotonic and cannot describe a VTK grid axis",
            coord.name
        ));
    }
    let tolerance = 1e-9 * (last - first).abs().max(first.abs()).max(last.abs());
    let uniform = values
        .iter()
        .enumerate()
        .all(|(i, &v)| (v - (first + step * i as f64)).abs() <= tolerance);
    if uniform {
        Ok(VtkAxis::Uniform {
            origin: first,
            spacing: step,
        })
    } else {
        Ok(VtkAxis::Explicit(values.clone()))
    }
}

/// A VTK legacy data name is a single token: whitespace is replaced by `_`.
fn vtk_data_name(name: &str) -> String {
    let token: String = name
        .chars()
        .map(|c| if c.is_whitespace() { '_' } else { c })
        .collect();
    if token.is_empty() {
        "data".to_string()
    } else {
        token
    }
}

/// Write the legacy-VTK document for `var`.
fn write_vtk_document<W: Write>(
    w: &mut W,
    var: &XarrayDataArray,
    axes: &[VtkAxis; 3],
    data_name: &str,
) -> std::io::Result<()> {
    let (nx, ny, nz) = (var.shape[0], var.shape[1], var.shape[2]);
    writeln!(w, "# vtk DataFile Version 3.0")?;
    writeln!(w, "XarrayExport")?;
    writeln!(w, "ASCII")?;
    match axes {
        [
            VtkAxis::Uniform {
                origin: ox,
                spacing: sx,
            },
            VtkAxis::Uniform {
                origin: oy,
                spacing: sy,
            },
            VtkAxis::Uniform {
                origin: oz,
                spacing: sz,
            },
        ] => {
            writeln!(w, "DATASET STRUCTURED_POINTS")?;
            writeln!(w, "DIMENSIONS {nx} {ny} {nz}")?;
            writeln!(w, "ORIGIN {ox} {oy} {oz}")?;
            writeln!(w, "SPACING {sx} {sy} {sz}")?;
        }
        _ => {
            writeln!(w, "DATASET RECTILINEAR_GRID")?;
            writeln!(w, "DIMENSIONS {nx} {ny} {nz}")?;
            for (label, axis, n) in [
                ("X", &axes[0], nx),
                ("Y", &axes[1], ny),
                ("Z", &axes[2], nz),
            ] {
                writeln!(w, "{label}_COORDINATES {n} double")?;
                let line: Vec<String> = axis.points(n).iter().map(|v| v.to_string()).collect();
                writeln!(w, "{}", line.join(" "))?;
            }
        }
    }
    writeln!(w, "POINT_DATA {}", nx * ny * nz)?;
    writeln!(w, "SCALARS {data_name} double 1")?;
    writeln!(w, "LOOKUP_TABLE default")?;
    // The array is C-order (the last dimension varies fastest) but VTK
    // requires the first DIMENSIONS entry (x) to vary fastest, so the values
    // are emitted transposed: x innermost, z outermost.
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                writeln!(w, "{}", var.data[linear_index(&[i, j, k], &var.shape)])?;
            }
        }
    }
    Ok(())
}

/// Export a 3-D variable from a dataset as a VTK legacy ASCII grid.
///
/// The variable must have exactly three non-empty dimensions and exactly
/// `shape[0] * shape[1] * shape[2]` values, otherwise an `Err` is returned.
/// Dimension 0 becomes the VTK x axis, dimension 1 y and dimension 2 z.
///
/// The data are stored C-order (the last dimension varies fastest) while VTK
/// expects x, the first entry of `DIMENSIONS`, to vary fastest, so the values
/// are written transposed.
///
/// The grid geometry comes from the coordinates whose name equals the
/// dimension name (a dimension without one uses the index `0, 1, 2, ...`):
///
/// * if every axis is evenly spaced the output is a `STRUCTURED_POINTS`
///   dataset with `ORIGIN` and `SPACING` taken from the coordinates;
/// * if any coordinate is not evenly spaced the output is a
///   `RECTILINEAR_GRID` dataset with explicit `X/Y/Z_COORDINATES`.
///
/// A coordinate must have one value per point of its dimension and be
/// strictly monotonic, else an `Err` is returned. Whitespace in the variable
/// name is replaced by `_` because VTK data names are single tokens.
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
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(format!(
            "variable '{var_name}' has an empty dimension (shape {:?})",
            var.shape
        ));
    }
    let n_points = nx
        .checked_mul(ny)
        .and_then(|n| n.checked_mul(nz))
        .ok_or_else(|| format!("variable '{var_name}': shape {:?} overflows", var.shape))?;
    if var.data.len() != n_points {
        return Err(format!(
            "variable '{var_name}' holds {} values but shape {:?} requires {n_points}",
            var.data.len(),
            var.shape
        ));
    }
    let axes = [vtk_axis(var, 0)?, vtk_axis(var, 1)?, vtk_axis(var, 2)?];

    let file = File::create(path).map_err(|e| format!("cannot create {path}: {e}"))?;
    let mut w = BufWriter::new(file);
    write_vtk_document(&mut w, var, &axes, &vtk_data_name(var_name))
        .map_err(|e| format!("cannot write {path}: {e}"))?;
    // A BufWriter drops flush errors silently: surface them.
    w.flush().map_err(|e| format!("cannot write {path}: {e}"))
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
        let tmpdir = std::env::temp_dir();
        let path = tmpdir
            .join("xarray_test_roundtrip")
            .to_str()
            .unwrap_or("")
            .to_string();
        write_csv_xarray(&ds, &path).expect("write failed");
        let loaded = read_csv_xarray(&format!("{path}_temperature.csv")).expect("read failed");
        assert!((loaded.get(&[0, 1]) - 3.125).abs() < 1e-9);
        assert!((loaded.get(&[1, 2]) - 2.72).abs() < 1e-9);
    }

    #[test]
    fn test_write_csv_creates_file() {
        let mut ds = XarrayDataset::new("test2");
        ds.add_variable(XarrayDataArray::new("v", vec!["x".into()], vec![3]));
        let base = std::env::temp_dir()
            .join("xarray_test2")
            .to_str()
            .unwrap_or("")
            .to_string();
        write_csv_xarray(&ds, &base).expect("write failed");
        let expected = std::env::temp_dir().join("xarray_test2_v.csv");
        assert!(expected.exists());
    }

    // ── VTK export tests ──────────────────────────────────────────────────

    /// Scratch path unique to this process and test.
    fn vtk_path(tag: &str) -> std::path::PathBuf {
        std::env::temp_dir().join(format!("xarray_vtk_{tag}_{}.vtk", std::process::id()))
    }

    fn dims_xyz() -> Vec<String> {
        vec!["x".into(), "y".into(), "z".into()]
    }

    /// A 3-D variable `name` of the given shape whose value at flat C-order
    /// index `n` is `n`.
    fn counting_variable(name: &str, shape: [usize; 3]) -> XarrayDataArray {
        let mut var = XarrayDataArray::new(name, dims_xyz(), shape.to_vec());
        for (n, v) in var.data.iter_mut().enumerate() {
            *v = n as f64;
        }
        var
    }

    /// Export `var` and return the written text; the scratch file is removed.
    fn export_text(var: XarrayDataArray, tag: &str) -> Result<String, String> {
        let name = var.name.clone();
        let mut ds = XarrayDataset::new("vtk_test");
        ds.add_variable(var);
        let path = vtk_path(tag);
        let result = xarray_to_vtk_structured(&ds, &name, path.to_str().unwrap_or(""));
        let text = std::fs::read_to_string(&path).unwrap_or_default();
        let _ = std::fs::remove_file(&path);
        result.map(|()| text)
    }

    /// The scalar values listed after `LOOKUP_TABLE default`.
    fn vtk_values(text: &str) -> Vec<f64> {
        text.lines()
            .skip_while(|l| *l != "LOOKUP_TABLE default")
            .skip(1)
            .map(|l| l.parse::<f64>().unwrap())
            .collect()
    }

    #[test]
    fn test_vtk_export_creates_file() {
        let mut ds = XarrayDataset::new("vtk_test");
        let arr = XarrayDataArray::new(
            "pressure",
            vec!["x".into(), "y".into(), "z".into()],
            vec![2, 2, 2],
        );
        ds.add_variable(arr);
        let path = std::env::temp_dir().join("xarray_test_pressure.vtk");
        xarray_to_vtk_structured(&ds, "pressure", path.to_str().unwrap_or("")).expect("vtk failed");
        assert!(path.exists());
    }

    #[test]
    fn test_vtk_export_wrong_dims_returns_err() {
        let mut ds = XarrayDataset::new("vtk_test2");
        let arr = XarrayDataArray::new("u2d", vec!["x".into(), "y".into()], vec![2, 2]);
        ds.add_variable(arr);
        let path = std::env::temp_dir().join("xarray_test_u2d.vtk");
        let res = xarray_to_vtk_structured(&ds, "u2d", path.to_str().unwrap_or(""));
        assert!(res.is_err());
    }

    #[test]
    fn test_vtk_export_missing_variable_returns_err() {
        let ds = XarrayDataset::new("vtk_empty");
        let path = vtk_path("missing");
        let err = xarray_to_vtk_structured(&ds, "nope", path.to_str().unwrap_or("")).unwrap_err();
        assert!(err.contains("'nope' not found"), "{err}");
        assert!(!path.exists());
    }

    #[test]
    fn test_vtk_export_writes_x_fastest() {
        // shape [nx=2, ny=1, nz=3], C-order data 0..6: flat index of (i, 0, k)
        // is 3 * i + k. VTK order (x fastest) visits (0,0,0) (1,0,0) (0,0,1)
        // (1,0,1) (0,0,2) (1,0,2) = flat 0, 3, 1, 4, 2, 5.
        let text = export_text(counting_variable("f", [2, 1, 3]), "xfastest").unwrap();
        assert!(text.contains("DIMENSIONS 2 1 3\n"), "{text}");
        assert!(text.contains("POINT_DATA 6\n"), "{text}");
        assert_eq!(vtk_values(&text), vec![0.0, 3.0, 1.0, 4.0, 2.0, 5.0]);
    }

    #[test]
    fn test_vtk_export_non_cubic_array_is_not_scrambled() {
        let (nx, ny, nz) = (2usize, 3usize, 4usize);
        let var = counting_variable("f", [nx, ny, nz]);
        let shape = var.shape.clone();
        let text = export_text(var, "noncubic").unwrap();
        let values = vtk_values(&text);
        assert_eq!(values.len(), nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let vtk_pos = i + nx * (j + ny * k);
                    let expected = linear_index(&[i, j, k], &shape) as f64;
                    assert!(
                        (values[vtk_pos] - expected).abs() < 1e-12,
                        "point ({i},{j},{k}): got {} expected {expected}",
                        values[vtk_pos]
                    );
                }
            }
        }
    }

    #[test]
    fn test_vtk_export_default_geometry_is_the_index_grid() {
        let text = export_text(counting_variable("f", [2, 2, 2]), "defaultgeom").unwrap();
        assert!(text.contains("DATASET STRUCTURED_POINTS\n"), "{text}");
        assert!(text.contains("ORIGIN 0 0 0\n"), "{text}");
        assert!(text.contains("SPACING 1 1 1\n"), "{text}");
        assert!(text.contains("SCALARS f double 1\n"), "{text}");
    }

    #[test]
    fn test_vtk_export_uniform_coordinates_set_origin_and_spacing() {
        let mut var = counting_variable("f", [2, 1, 3]);
        var.set_coord(XarrayCoordinate::new("x", vec![10.0, 12.0], "m"));
        var.set_coord(XarrayCoordinate::new("y", vec![5.0], "m"));
        var.set_coord(XarrayCoordinate::new("z", vec![0.0, 0.5, 1.0], "m"));
        let text = export_text(var, "uniformcoords").unwrap();
        assert!(text.contains("DATASET STRUCTURED_POINTS\n"), "{text}");
        // A single-point axis keeps its position and gets unit spacing.
        assert!(text.contains("ORIGIN 10 5 0\n"), "{text}");
        assert!(text.contains("SPACING 2 1 0.5\n"), "{text}");
    }

    #[test]
    fn test_vtk_export_descending_uniform_coordinate_has_negative_spacing() {
        let mut var = counting_variable("f", [1, 1, 3]);
        var.set_coord(XarrayCoordinate::new("z", vec![1.0, 0.5, 0.0], "m"));
        let text = export_text(var, "descending").unwrap();
        assert!(text.contains("ORIGIN 0 0 1\n"), "{text}");
        assert!(text.contains("SPACING 1 1 -0.5\n"), "{text}");
    }

    #[test]
    fn test_vtk_export_non_uniform_coordinates_write_a_rectilinear_grid() {
        let mut var = counting_variable("f", [2, 1, 3]);
        var.set_coord(XarrayCoordinate::new("z", vec![0.0, 1.0, 3.0], "m"));
        let text = export_text(var, "rectilinear").unwrap();
        assert!(text.contains("DATASET RECTILINEAR_GRID\n"), "{text}");
        assert!(!text.contains("SPACING"), "{text}");
        // x and y have no coordinate: index axes.
        assert!(text.contains("X_COORDINATES 2 double\n0 1\n"), "{text}");
        assert!(text.contains("Y_COORDINATES 1 double\n0\n"), "{text}");
        assert!(text.contains("Z_COORDINATES 3 double\n0 1 3\n"), "{text}");
        // The data are still written x fastest.
        assert_eq!(vtk_values(&text), vec![0.0, 3.0, 1.0, 4.0, 2.0, 5.0]);
    }

    #[test]
    fn test_vtk_export_ignores_coordinates_not_named_after_a_dimension() {
        let mut var = counting_variable("f", [2, 2, 2]);
        var.set_coord(XarrayCoordinate::new(
            "longitude",
            vec![100.0, 200.0],
            "deg",
        ));
        let text = export_text(var, "unmatched").unwrap();
        assert!(text.contains("ORIGIN 0 0 0\n"), "{text}");
    }

    #[test]
    fn test_vtk_export_rejects_data_length_mismatch() {
        let mut var = counting_variable("f", [2, 2, 2]);
        var.data.pop();
        let err = export_text(var, "badlen").unwrap_err();
        assert!(err.contains("holds 7 values"), "{err}");
        assert!(err.contains("requires 8"), "{err}");
    }

    #[test]
    fn test_vtk_export_rejects_empty_dimension() {
        let err = export_text(counting_variable("f", [2, 0, 2]), "empty").unwrap_err();
        assert!(err.contains("empty dimension"), "{err}");
    }

    #[test]
    fn test_vtk_export_rejects_coordinate_of_wrong_length() {
        let mut var = counting_variable("f", [2, 2, 2]);
        var.set_coord(XarrayCoordinate::new("x", vec![0.0, 1.0, 2.0], "m"));
        let err = export_text(var, "coordlen").unwrap_err();
        assert!(err.contains("coordinate 'x' has 3 values"), "{err}");
        assert!(err.contains("size 2"), "{err}");
    }

    #[test]
    fn test_vtk_export_rejects_non_monotonic_or_non_finite_coordinates() {
        let mut zigzag = counting_variable("f", [1, 1, 3]);
        zigzag.set_coord(XarrayCoordinate::new("z", vec![0.0, 2.0, 1.0], "m"));
        let err = export_text(zigzag, "zigzag").unwrap_err();
        assert!(err.contains("not strictly monotonic"), "{err}");

        let mut repeated = counting_variable("f", [1, 1, 3]);
        repeated.set_coord(XarrayCoordinate::new("z", vec![0.0, 0.0, 0.0], "m"));
        assert!(export_text(repeated, "repeated").is_err());

        let mut nan = counting_variable("f", [1, 1, 2]);
        nan.set_coord(XarrayCoordinate::new("z", vec![0.0, f64::NAN], "m"));
        let err = export_text(nan, "nan").unwrap_err();
        assert!(err.contains("non-finite"), "{err}");
    }

    #[test]
    fn test_vtk_export_data_name_is_a_single_token() {
        let text =
            export_text(counting_variable("air temperature", [1, 1, 2]), "dataname").unwrap();
        assert!(
            text.contains("SCALARS air_temperature double 1\n"),
            "{text}"
        );
    }

    #[test]
    fn test_vtk_export_unwritable_path_returns_err() {
        let mut ds = XarrayDataset::new("vtk_test");
        ds.add_variable(counting_variable("f", [1, 1, 1]));
        let path = std::env::temp_dir()
            .join(format!("xarray_no_such_dir_{}", std::process::id()))
            .join("out.vtk");
        let err = xarray_to_vtk_structured(&ds, "f", path.to_str().unwrap_or("")).unwrap_err();
        assert!(err.contains("cannot create"), "{err}");
    }

    #[test]
    fn test_coordinate_new() {
        let c = XarrayCoordinate::new("time", vec![0.0, 1.0, 2.0], "s");
        assert_eq!(c.name, "time");
        assert_eq!(c.values.len(), 3);
    }
}
