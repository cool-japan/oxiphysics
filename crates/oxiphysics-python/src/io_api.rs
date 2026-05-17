// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! I/O API for Python interop.
//!
//! Provides lightweight in-memory representations of common simulation I/O
//! formats (VTK, CSV, XYZ, LAMMPS dump, HDF5-style, trajectory). All types
//! use plain `f64`, `Vec`f64`, and `String` — no nalgebra — for easy FFI.

use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// PyVtkWriter
// ---------------------------------------------------------------------------

/// In-memory VTK legacy-format writer (ASCII and binary).
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyVtkWriter {
    /// Output filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Named point data arrays.
    pub point_data: Vec<(String, Vec<f64>)>,
    /// Named cell data arrays.
    pub cell_data: Vec<(String, Vec<f64>)>,
    /// Point positions `[x, y, z]` for binary serialisation.
    pub positions: Vec<[f64; 3]>,
}

#[pymethods]
impl PyVtkWriter {
    /// Create a new VTK writer targeting the given file.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            point_data: Vec::new(),
            cell_data: Vec::new(),
            positions: Vec::new(),
        }
    }

    /// Add a point position `[x, y, z]`.
    pub fn add_position(&mut self, x: f64, y: f64, z: f64) {
        self.positions.push([x, y, z]);
    }

    /// Attach a named point-data array.
    pub fn add_point_data(&mut self, name: String, data: Vec<f64>) {
        self.point_data.push((name, data));
    }

    /// Attach a named cell-data array.
    pub fn add_cell_data(&mut self, name: String, data: Vec<f64>) {
        self.cell_data.push((name, data));
    }

    /// Serialise to a minimal VTK ASCII string (stub — does not write disk I/O).
    pub fn write_ascii(&self) -> String {
        let mut out = format!(
            "# vtk DataFile Version 3.0\nOxiPhysics output\nASCII\nDATASET UNSTRUCTURED_GRID\nfile={}\n",
            self.filename
        );
        for (name, data) in &self.point_data {
            out.push_str(&format!("POINT_DATA {} len={}\n", name, data.len()));
        }
        for (name, data) in &self.cell_data {
            out.push_str(&format!("CELL_DATA {} len={}\n", name, data.len()));
        }
        out
    }

    /// Serialise to VTK legacy binary format bytes.
    ///
    /// Produces a valid VTK legacy binary file with `UNSTRUCTURED_GRID` points.
    /// Point coordinates are written as big-endian `f32` values as required by
    /// the VTK legacy binary format specification.
    pub fn to_vtk_binary(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"# vtk DataFile Version 3.0\n");
        buf.extend_from_slice(b"OxiPhysics data\n");
        buf.extend_from_slice(b"BINARY\n");
        buf.extend_from_slice(b"DATASET UNSTRUCTURED_GRID\n");

        let n = self.positions.len();
        buf.extend_from_slice(format!("POINTS {} float\n", n).as_bytes());
        for pos in &self.positions {
            buf.extend_from_slice(&(pos[0] as f32).to_be_bytes());
            buf.extend_from_slice(&(pos[1] as f32).to_be_bytes());
            buf.extend_from_slice(&(pos[2] as f32).to_be_bytes());
        }

        // Write degenerate CELLS and CELL_TYPES sections so the file is valid.
        // Each point is a vertex cell (VTK type 1).
        if n > 0 {
            buf.extend_from_slice(format!("CELLS {} {}\n", n, n * 2).as_bytes());
            for _ in 0..n {
                buf.extend_from_slice(&1_i32.to_be_bytes()); // n_verts = 1
                buf.extend_from_slice(&0_i32.to_be_bytes()); // point index placeholder (0)
            }
            buf.extend_from_slice(format!("CELL_TYPES {}\n", n).as_bytes());
            for _ in 0..n {
                buf.extend_from_slice(&1_i32.to_be_bytes()); // VTK_VERTEX = 1
            }
        }

        // Named point scalar fields.
        if !self.point_data.is_empty() {
            let n_pts = self.positions.len();
            buf.extend_from_slice(format!("POINT_DATA {}\n", n_pts).as_bytes());
            for (name, data) in &self.point_data {
                buf.extend_from_slice(
                    format!("SCALARS {} float 1\nLOOKUP_TABLE default\n", name).as_bytes(),
                );
                for &v in data {
                    buf.extend_from_slice(&(v as f32).to_be_bytes());
                }
            }
        }

        // Named cell scalar fields.
        if !self.cell_data.is_empty() {
            buf.extend_from_slice(format!("CELL_DATA {}\n", n).as_bytes());
            for (name, data) in &self.cell_data {
                buf.extend_from_slice(
                    format!("SCALARS {} float 1\nLOOKUP_TABLE default\n", name).as_bytes(),
                );
                for &v in data {
                    buf.extend_from_slice(&(v as f32).to_be_bytes());
                }
            }
        }

        buf
    }

    /// Write the binary VTK data to `path`.
    pub fn write_to_file(&self, path: &str) -> PyResult<()> {
        std::fs::write(path, self.to_vtk_binary())
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Returns the byte length that a binary VTK file would have (kept for
    /// backward compatibility; prefer `to_vtk_binary().len()` for accuracy).
    pub fn write_binary(&self) -> usize {
        let base = self.filename.len() + 64;
        let pd: usize = self.point_data.iter().map(|(_, v)| v.len() * 8).sum();
        let cd: usize = self.cell_data.iter().map(|(_, v)| v.len() * 8).sum();
        base + pd + cd
    }

    /// Number of point-data arrays attached.
    pub fn n_point_arrays(&self) -> usize {
        self.point_data.len()
    }

    /// Number of cell-data arrays attached.
    pub fn n_cell_arrays(&self) -> usize {
        self.cell_data.len()
    }
}

impl Default for PyVtkWriter {
    fn default() -> Self {
        Self::new("output.vtk".to_string())
    }
}

// ---------------------------------------------------------------------------
// PyCsvReader
// ---------------------------------------------------------------------------

/// In-memory CSV reader.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCsvReader {
    /// Source filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Column header names (populated from the first row if `has_header = true`).
    #[pyo3(get, set)]
    pub headers: Vec<String>,
    /// Row-major data storage (rows × columns).
    #[pyo3(get, set)]
    pub rows: Vec<Vec<f64>>,
}

#[pymethods]
impl PyCsvReader {
    /// Create a new CSV reader backed by an in-memory dataset.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            headers: Vec::new(),
            rows: Vec::new(),
        }
    }

    /// Load in-memory data directly (useful in tests without touching disk).
    pub fn load_data(&mut self, headers: Vec<String>, rows: Vec<Vec<f64>>) {
        self.headers = headers;
        self.rows = rows;
    }

    /// Read a single column by index, returning a `Vec<f64>`.
    pub fn read_column(&self, col: usize) -> Vec<f64> {
        self.rows
            .iter()
            .filter_map(|r| r.get(col).copied())
            .collect()
    }

    /// Return all data as a flat `Vec<f64>` in row-major order.
    pub fn read_all_f64(&self) -> Vec<f64> {
        self.rows.iter().flat_map(|r| r.iter().copied()).collect()
    }

    /// Return the column header names.
    pub fn header_names(&self) -> Vec<String> {
        self.headers.clone()
    }

    /// Number of data rows.
    pub fn n_rows(&self) -> usize {
        self.rows.len()
    }

    /// Number of columns (inferred from the first row, or 0 if empty).
    pub fn n_cols(&self) -> usize {
        self.rows.first().map_or(0, |r| r.len())
    }
}

impl Default for PyCsvReader {
    fn default() -> Self {
        Self::new("input.csv".to_string())
    }
}

// ---------------------------------------------------------------------------
// PyCsvWriter
// ---------------------------------------------------------------------------

/// In-memory CSV writer.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCsvWriter {
    /// Output filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Buffered rows waiting to be flushed.
    #[pyo3(get, set)]
    pub buffer: Vec<Vec<f64>>,
}

#[pymethods]
impl PyCsvWriter {
    /// Create a new CSV writer targeting the given file.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            buffer: Vec::new(),
        }
    }

    /// Append a row of values to the internal buffer.
    pub fn write_row(&mut self, data: Vec<f64>) {
        self.buffer.push(data);
    }

    /// Build a CSV string from the current buffer without clearing it.
    pub fn to_csv_string(&self) -> String {
        self.buffer
            .iter()
            .map(|row| {
                row.iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>()
            .join("\n")
    }

    /// Consume the buffer and return it as a CSV string.
    pub fn flush(&mut self) -> String {
        let csv = self.to_csv_string();
        self.buffer.clear();
        csv
    }

    /// Write the buffered CSV data to `path` without clearing the buffer.
    pub fn write_to_file(&self, path: &str) -> PyResult<()> {
        std::fs::write(path, self.to_csv_string())
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Number of rows currently buffered.
    pub fn buffered_rows(&self) -> usize {
        self.buffer.len()
    }
}

impl Default for PyCsvWriter {
    fn default() -> Self {
        Self::new("output.csv".to_string())
    }
}

// ---------------------------------------------------------------------------
// PyXyzReader
// ---------------------------------------------------------------------------

/// In-memory reader for the XYZ molecular dynamics format.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyXyzReader {
    /// Source filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Atom species strings (e.g., `"C"`, `"H"`, `"O"`).
    #[pyo3(get, set)]
    pub species: Vec<String>,
    /// Flat positions: `[x0, y0, z0, x1, y1, z1, …]`.
    #[pyo3(get, set)]
    pub pos_flat: Vec<f64>,
}

#[pymethods]
impl PyXyzReader {
    /// Create a new XYZ reader.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            species: Vec::new(),
            pos_flat: Vec::new(),
        }
    }

    /// Load in-memory data (useful in tests without touching disk).
    pub fn load_data(&mut self, species: Vec<String>, pos_flat: Vec<f64>) {
        self.species = species;
        self.pos_flat = pos_flat;
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.species.len()
    }

    /// Position array `[x, y, z, …]` for all atoms.
    pub fn positions(&self) -> Vec<f64> {
        self.pos_flat.clone()
    }

    /// Species array.
    pub fn get_species(&self) -> Vec<String> {
        self.species.clone()
    }

    /// Position `[x, y, z]` of atom `i`.
    pub fn position_of(&self, i: usize) -> Option<Vec<f64>> {
        let base = i * 3;
        if base + 2 < self.pos_flat.len() {
            Some(vec![
                self.pos_flat[base],
                self.pos_flat[base + 1],
                self.pos_flat[base + 2],
            ])
        } else {
            None
        }
    }
}

impl Default for PyXyzReader {
    fn default() -> Self {
        Self::new("atoms.xyz".to_string())
    }
}

// ---------------------------------------------------------------------------
// PyXyzWriter
// ---------------------------------------------------------------------------

/// In-memory writer for the XYZ molecular dynamics format.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyXyzWriter {
    /// Output filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Accumulated frames as raw XYZ strings.
    #[pyo3(get, set)]
    pub frames: Vec<String>,
}

#[pymethods]
impl PyXyzWriter {
    /// Create a new XYZ writer.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            frames: Vec::new(),
        }
    }

    /// Append a frame to the internal buffer.
    ///
    /// `positions` is flat `[x0,y0,z0, x1,y1,z1, …]`; `species` has one entry
    /// per atom; `comment` is written to the second header line.
    pub fn write_frame(&mut self, positions: Vec<f64>, species: Vec<String>, comment: String) {
        let n = species.len();
        let mut frame = format!("{}\n{}\n", n, comment);
        for (i, sp) in species.iter().enumerate() {
            let base = i * 3;
            let (x, y, z) = if base + 2 < positions.len() {
                (positions[base], positions[base + 1], positions[base + 2])
            } else {
                (0.0, 0.0, 0.0)
            };
            frame.push_str(&format!("{} {} {} {}\n", sp, x, y, z));
        }
        self.frames.push(frame);
    }

    /// Number of frames written so far.
    pub fn n_frames(&self) -> usize {
        self.frames.len()
    }

    /// Return all frames concatenated as a single string.
    pub fn as_string(&self) -> String {
        self.frames.concat()
    }
}

impl Default for PyXyzWriter {
    fn default() -> Self {
        Self::new("output.xyz".to_string())
    }
}

// ---------------------------------------------------------------------------
// PyLammpsReader
// ---------------------------------------------------------------------------

/// In-memory reader for LAMMPS dump files.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyLammpsReader {
    /// Source filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Flat atom data: `[id, type, x, y, z, …]` per atom (5 fields each).
    #[pyo3(get, set)]
    pub atom_data: Vec<Vec<f64>>,
    /// Simulation box bounds `[[xlo, xhi], [ylo, yhi], [zlo, zhi]]`.
    pub box_bounds: [[f64; 2]; 3],
}

#[pymethods]
impl PyLammpsReader {
    /// Create a new LAMMPS reader.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            atom_data: Vec::new(),
            box_bounds: [[0.0, 1.0]; 3],
        }
    }

    /// Load in-memory atom data and box bounds.
    ///
    /// `box_bounds_flat` is a flat `[xlo, xhi, ylo, yhi, zlo, zhi]` list.
    pub fn load_data(&mut self, atom_data: Vec<Vec<f64>>, box_bounds_flat: Vec<f64>) {
        self.atom_data = atom_data;
        if box_bounds_flat.len() >= 6 {
            self.box_bounds = [
                [box_bounds_flat[0], box_bounds_flat[1]],
                [box_bounds_flat[2], box_bounds_flat[3]],
                [box_bounds_flat[4], box_bounds_flat[5]],
            ];
        }
    }

    /// Return all atom records.
    pub fn read_atoms(&self) -> Vec<Vec<f64>> {
        self.atom_data.clone()
    }

    /// Number of atoms in the last read frame.
    pub fn n_atoms(&self) -> usize {
        self.atom_data.len()
    }

    /// Box bounds as flat `[xlo, xhi, ylo, yhi, zlo, zhi]`.
    pub fn box_bounds_flat(&self) -> Vec<f64> {
        vec![
            self.box_bounds[0][0],
            self.box_bounds[0][1],
            self.box_bounds[1][0],
            self.box_bounds[1][1],
            self.box_bounds[2][0],
            self.box_bounds[2][1],
        ]
    }

    /// Box side lengths `[Lx, Ly, Lz]`.
    pub fn box_lengths(&self) -> Vec<f64> {
        vec![
            self.box_bounds[0][1] - self.box_bounds[0][0],
            self.box_bounds[1][1] - self.box_bounds[1][0],
            self.box_bounds[2][1] - self.box_bounds[2][0],
        ]
    }
}

impl Default for PyLammpsReader {
    fn default() -> Self {
        Self::new("dump.lammpstrj".to_string())
    }
}

// ---------------------------------------------------------------------------
// PyHdf5Writer
// ---------------------------------------------------------------------------

/// In-memory HDF5-style writer (no actual HDF5 dependency).
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyHdf5Writer {
    /// Output filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Named datasets.
    pub datasets: Vec<(String, Vec<f64>)>,
    /// Named scalar attributes.
    pub attributes: Vec<(String, f64)>,
}

#[pymethods]
impl PyHdf5Writer {
    /// Create a new HDF5 writer.
    #[new]
    pub fn new(filename: String) -> Self {
        Self {
            filename,
            datasets: Vec::new(),
            attributes: Vec::new(),
        }
    }

    /// Write a named dataset.
    pub fn write_dataset(&mut self, name: String, data: Vec<f64>) {
        self.datasets.push((name, data));
    }

    /// Write a named scalar attribute.
    pub fn write_attribute(&mut self, name: String, value: f64) {
        self.attributes.push((name, value));
    }

    /// Number of datasets stored.
    pub fn n_datasets(&self) -> usize {
        self.datasets.len()
    }

    /// Number of attributes stored.
    pub fn n_attributes(&self) -> usize {
        self.attributes.len()
    }

    /// Retrieve a dataset by name.
    pub fn get_dataset(&self, name: &str) -> Option<Vec<f64>> {
        self.datasets
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, d)| d.clone())
    }

    /// Retrieve a scalar attribute by name.
    pub fn get_attribute(&self, name: &str) -> Option<f64> {
        self.attributes
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| *v)
    }
}

impl Default for PyHdf5Writer {
    fn default() -> Self {
        Self::new("output.h5".to_string())
    }
}

// ---------------------------------------------------------------------------
// TrajectoryFrame (Rust-internal only, not exposed as pyclass)
// ---------------------------------------------------------------------------

/// A single trajectory frame for structured multi-format output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryFrame {
    /// Simulation step number.
    pub step: usize,
    /// Simulation time.
    pub time: f64,
    /// Atom positions as `[x, y, z]` arrays.
    pub positions: Vec<[f64; 3]>,
    /// Element symbols, one per atom (e.g., `"C"`, `"H"`).
    pub element_symbols: Vec<String>,
    /// Simulation box side lengths `[Lx, Ly, Lz]`, if periodic.
    pub box_lengths: Option<[f64; 3]>,
}

// ---------------------------------------------------------------------------
// PyTrajectoryWriter
// ---------------------------------------------------------------------------

/// Multi-format trajectory writer.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTrajectoryWriter {
    /// Output filename.
    #[pyo3(get, set)]
    pub filename: String,
    /// Format string (e.g., `"xyz"`, `"lammps"`, `"vtk"`).
    #[pyo3(get, set)]
    pub format: String,
    /// Whether the writer has been closed.
    #[pyo3(get)]
    pub closed: bool,
    /// Number of frames written.
    #[pyo3(get)]
    pub frame_count: usize,
    /// Buffered frame strings (legacy raw-string API).
    pub frame_buffer: Vec<String>,
    /// Structured frames for the typed format-dispatch API.
    pub typed_frames: Vec<TrajectoryFrame>,
}

#[pymethods]
impl PyTrajectoryWriter {
    /// Create a new trajectory writer.
    #[new]
    pub fn new(filename: String, format: String) -> Self {
        Self {
            filename,
            format,
            closed: false,
            frame_count: 0,
            frame_buffer: Vec::new(),
            typed_frames: Vec::new(),
        }
    }

    /// Add a structured frame for format-dispatch serialisation.
    ///
    /// `positions_flat` is `[x0,y0,z0, x1,y1,z1, …]`; `element_symbols` has
    /// one entry per atom; `box_lengths_flat` is optional `[Lx,Ly,Lz]`.
    pub fn add_frame_typed(
        &mut self,
        step: usize,
        time: f64,
        positions_flat: Vec<f64>,
        element_symbols: Vec<String>,
        box_lengths_flat: Option<Vec<f64>>,
    ) {
        let positions: Vec<[f64; 3]> = positions_flat
            .chunks(3)
            .map(|c| {
                let x = c.first().copied().unwrap_or(0.0);
                let y = c.get(1).copied().unwrap_or(0.0);
                let z = c.get(2).copied().unwrap_or(0.0);
                [x, y, z]
            })
            .collect();
        let box_lengths = box_lengths_flat.and_then(|v| {
            if v.len() >= 3 {
                Some([v[0], v[1], v[2]])
            } else {
                None
            }
        });
        self.typed_frames.push(TrajectoryFrame {
            step,
            time,
            positions,
            element_symbols,
            box_lengths,
        });
    }

    /// Serialise all typed frames to a string using the configured format.
    pub fn write_to_string(&self) -> PyResult<String> {
        match self.format.to_lowercase().as_str() {
            "xyz" => Ok(self.write_xyz()),
            "lammps" | "lammps-dump" => Ok(self.write_lammps_dump()),
            "vtk" => Ok(self.write_vtk_ascii()),
            fmt => Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Unknown trajectory format: {}",
                fmt
            ))),
        }
    }

    /// Write all typed frames to `path` using the configured format.
    pub fn write_typed_to_file(&self, path: &str) -> PyResult<()> {
        let content = self.write_to_string()?;
        std::fs::write(path, content)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    /// Write a trajectory frame (legacy raw-string API).
    ///
    /// `positions` and `velocities` are flat `[x,y,z, …]` arrays; `step` is
    /// the integer simulation step number.
    pub fn write_frame(&mut self, positions: Vec<f64>, velocities: Vec<f64>, step: u64) {
        if self.closed {
            return;
        }
        let frame = format!(
            "FRAME step={} n_pos={} n_vel={} fmt={}\n",
            step,
            positions.len(),
            velocities.len(),
            self.format
        );
        self.frame_buffer.push(frame);
        self.frame_count += 1;
    }

    /// Close the writer (no further frames can be added).
    pub fn close(&mut self) {
        self.closed = true;
    }

    /// Return true if the writer is closed.
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// Total number of frames written.
    pub fn n_frames(&self) -> usize {
        self.frame_count
    }

    /// Return buffered content as a string.
    pub fn as_string(&self) -> String {
        self.frame_buffer.concat()
    }
}

impl PyTrajectoryWriter {
    /// Add a structured [`TrajectoryFrame`] for format-dispatch serialisation.
    pub fn add_frame(&mut self, frame: TrajectoryFrame) {
        self.typed_frames.push(frame);
    }

    /// Serialise typed frames to XYZ format.
    fn write_xyz(&self) -> String {
        let mut out = String::new();
        for frame in &self.typed_frames {
            out.push_str(&format!("{}\n", frame.positions.len()));
            out.push_str(&format!("frame {} t={:.6}\n", frame.step, frame.time));
            for (sym, pos) in frame.element_symbols.iter().zip(frame.positions.iter()) {
                out.push_str(&format!(
                    "{} {:.6} {:.6} {:.6}\n",
                    sym, pos[0], pos[1], pos[2]
                ));
            }
        }
        out
    }

    /// Serialise typed frames to LAMMPS dump format.
    fn write_lammps_dump(&self) -> String {
        let mut out = String::new();
        for frame in &self.typed_frames {
            out.push_str("ITEM: TIMESTEP\n");
            out.push_str(&format!("{}\n", frame.step));
            out.push_str("ITEM: NUMBER OF ATOMS\n");
            out.push_str(&format!("{}\n", frame.positions.len()));
            if let Some(box_l) = frame.box_lengths {
                out.push_str("ITEM: BOX BOUNDS pp pp pp\n");
                out.push_str(&format!(
                    "0 {:.6}\n0 {:.6}\n0 {:.6}\n",
                    box_l[0], box_l[1], box_l[2]
                ));
            }
            out.push_str("ITEM: ATOMS id type x y z\n");
            for (i, pos) in frame.positions.iter().enumerate() {
                out.push_str(&format!(
                    "{} 1 {:.6} {:.6} {:.6}\n",
                    i + 1,
                    pos[0],
                    pos[1],
                    pos[2]
                ));
            }
        }
        out
    }

    /// Serialise typed frames to VTK ASCII PolyData format.
    fn write_vtk_ascii(&self) -> String {
        let mut out = String::new();
        for (fi, frame) in self.typed_frames.iter().enumerate() {
            let n = frame.positions.len();
            out.push_str(&format!(
                "# vtk DataFile Version 3.0\nframe {} t={:.6}\nASCII\nDATASET POLYDATA\n",
                fi, frame.time
            ));
            out.push_str(&format!("POINTS {} float\n", n));
            for pos in &frame.positions {
                out.push_str(&format!("{:.6} {:.6} {:.6}\n", pos[0], pos[1], pos[2]));
            }
            // Vertex cells so the file is valid.
            out.push_str(&format!("VERTICES {} {}\n", n, n * 2));
            for i in 0..n {
                out.push_str(&format!("1 {}\n", i));
            }
        }
        out
    }
}

impl Default for PyTrajectoryWriter {
    fn default() -> Self {
        Self::new("trajectory.xyz".to_string(), "xyz".to_string())
    }
}

// ---------------------------------------------------------------------------
// Registration helper
// ---------------------------------------------------------------------------

/// Register all `io` classes into a Python sub-module.
///
/// Called from the top-level `#[pymodule]` in `lib.rs`.
pub fn register_io_module(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    use pyo3::types::PyModuleMethods;
    let child = PyModule::new(parent.py(), "io")?;
    child.add_class::<PyVtkWriter>()?;
    child.add_class::<PyCsvReader>()?;
    child.add_class::<PyCsvWriter>()?;
    child.add_class::<PyXyzReader>()?;
    child.add_class::<PyXyzWriter>()?;
    child.add_class::<PyLammpsReader>()?;
    child.add_class::<PyHdf5Writer>()?;
    child.add_class::<PyTrajectoryWriter>()?;
    parent.add_submodule(&child)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- PyVtkWriter ---

    #[test]
    fn test_vtk_new() {
        let w = PyVtkWriter::new("out.vtk".to_string());
        assert_eq!(w.filename, "out.vtk");
        assert_eq!(w.n_point_arrays(), 0);
    }

    #[test]
    fn test_vtk_add_point_data() {
        let mut w = PyVtkWriter::default();
        w.add_point_data("pressure".to_string(), vec![1.0, 2.0, 3.0]);
        assert_eq!(w.n_point_arrays(), 1);
    }

    #[test]
    fn test_vtk_add_cell_data() {
        let mut w = PyVtkWriter::default();
        w.add_cell_data("stress".to_string(), vec![10.0, 20.0]);
        assert_eq!(w.n_cell_arrays(), 1);
    }

    #[test]
    fn test_vtk_write_ascii_contains_header() {
        let w = PyVtkWriter::new("test.vtk".to_string());
        let s = w.write_ascii();
        assert!(s.contains("vtk DataFile"));
    }

    #[test]
    fn test_vtk_write_ascii_contains_point_data_name() {
        let mut w = PyVtkWriter::new("test.vtk".to_string());
        w.add_point_data("velocity".to_string(), vec![1.0, 2.0]);
        let s = w.write_ascii();
        assert!(s.contains("velocity"));
    }

    #[test]
    fn test_vtk_write_binary_size_grows() {
        let mut w = PyVtkWriter::new("test.vtk".to_string());
        let s0 = w.write_binary();
        w.add_point_data("p".to_string(), vec![1.0; 100]);
        let s1 = w.write_binary();
        assert!(s1 > s0);
    }

    #[test]
    fn test_vtk_default() {
        let w = PyVtkWriter::default();
        assert!(w.filename.ends_with(".vtk"));
    }

    // --- PyCsvReader ---

    #[test]
    fn test_csv_reader_new() {
        let r = PyCsvReader::new("data.csv".to_string());
        assert_eq!(r.filename, "data.csv");
        assert_eq!(r.n_rows(), 0);
    }

    #[test]
    fn test_csv_reader_load_and_read_column() {
        let mut r = PyCsvReader::default();
        r.load_data(
            vec!["x".to_string(), "y".to_string()],
            vec![vec![1.0, 2.0], vec![3.0, 4.0]],
        );
        let col0 = r.read_column(0);
        assert_eq!(col0, vec![1.0, 3.0]);
    }

    #[test]
    fn test_csv_reader_read_all_f64() {
        let mut r = PyCsvReader::default();
        r.load_data(vec![], vec![vec![1.0, 2.0], vec![3.0, 4.0]]);
        let all = r.read_all_f64();
        assert_eq!(all.len(), 4);
    }

    #[test]
    fn test_csv_reader_header_names() {
        let mut r = PyCsvReader::default();
        r.load_data(vec!["a".to_string(), "b".to_string()], vec![]);
        assert_eq!(r.header_names().len(), 2);
    }

    #[test]
    fn test_csv_reader_n_cols() {
        let mut r = PyCsvReader::default();
        r.load_data(vec![], vec![vec![1.0, 2.0, 3.0]]);
        assert_eq!(r.n_cols(), 3);
    }

    #[test]
    fn test_csv_reader_empty_n_cols_zero() {
        let r = PyCsvReader::default();
        assert_eq!(r.n_cols(), 0);
    }

    // --- PyCsvWriter ---

    #[test]
    fn test_csv_writer_new() {
        let w = PyCsvWriter::new("out.csv".to_string());
        assert_eq!(w.filename, "out.csv");
        assert_eq!(w.buffered_rows(), 0);
    }

    #[test]
    fn test_csv_writer_write_row() {
        let mut w = PyCsvWriter::default();
        w.write_row(vec![1.0, 2.0, 3.0]);
        assert_eq!(w.buffered_rows(), 1);
    }

    #[test]
    fn test_csv_writer_flush_clears_buffer() {
        let mut w = PyCsvWriter::default();
        w.write_row(vec![1.0]);
        w.flush();
        assert_eq!(w.buffered_rows(), 0);
    }

    #[test]
    fn test_csv_writer_flush_returns_csv() {
        let mut w = PyCsvWriter::default();
        w.write_row(vec![1.0, 2.0]);
        let s = w.flush();
        assert!(s.contains("1") && s.contains("2"));
    }

    #[test]
    fn test_csv_writer_default() {
        let w = PyCsvWriter::default();
        assert!(w.filename.ends_with(".csv"));
    }

    // --- PyXyzReader ---

    #[test]
    fn test_xyz_reader_new() {
        let r = PyXyzReader::new("mol.xyz".to_string());
        assert_eq!(r.filename, "mol.xyz");
        assert_eq!(r.n_atoms(), 0);
    }

    #[test]
    fn test_xyz_reader_load_and_n_atoms() {
        let mut r = PyXyzReader::default();
        r.load_data(
            vec!["C".to_string(), "H".to_string()],
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        );
        assert_eq!(r.n_atoms(), 2);
    }

    #[test]
    fn test_xyz_reader_positions() {
        let mut r = PyXyzReader::default();
        r.load_data(vec!["O".to_string()], vec![1.0, 2.0, 3.0]);
        assert_eq!(r.positions(), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_xyz_reader_species() {
        let mut r = PyXyzReader::default();
        r.load_data(vec!["N".to_string()], vec![0.0, 0.0, 0.0]);
        assert_eq!(r.get_species()[0], "N");
    }

    #[test]
    fn test_xyz_reader_position_of() {
        let mut r = PyXyzReader::default();
        r.load_data(
            vec!["C".to_string(), "H".to_string()],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        );
        let p = r.position_of(0).expect("position exists");
        assert_eq!(p, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_xyz_reader_default() {
        let r = PyXyzReader::default();
        assert!(r.filename.ends_with(".xyz"));
    }

    // --- PyXyzWriter ---

    #[test]
    fn test_xyz_writer_new() {
        let w = PyXyzWriter::new("out.xyz".to_string());
        assert_eq!(w.filename, "out.xyz");
        assert_eq!(w.n_frames(), 0);
    }

    #[test]
    fn test_xyz_writer_write_frame_increments_count() {
        let mut w = PyXyzWriter::default();
        w.write_frame(
            vec![0.0, 0.0, 0.0],
            vec!["C".to_string()],
            "frame 0".to_string(),
        );
        assert_eq!(w.n_frames(), 1);
    }

    #[test]
    fn test_xyz_writer_as_string_contains_n_atoms() {
        let mut w = PyXyzWriter::default();
        w.write_frame(
            vec![0.0, 0.0, 0.0],
            vec!["C".to_string()],
            "test".to_string(),
        );
        assert!(w.as_string().contains('1'));
    }

    #[test]
    fn test_xyz_writer_multiple_frames() {
        let mut w = PyXyzWriter::default();
        for _ in 0..5 {
            w.write_frame(vec![0.0, 0.0, 0.0], vec!["H".to_string()], String::new());
        }
        assert_eq!(w.n_frames(), 5);
    }

    #[test]
    fn test_xyz_writer_default() {
        let w = PyXyzWriter::default();
        assert!(w.filename.ends_with(".xyz"));
    }

    // --- PyLammpsReader ---

    #[test]
    fn test_lammps_reader_new() {
        let r = PyLammpsReader::new("dump.lammps".to_string());
        assert_eq!(r.n_atoms(), 0);
    }

    #[test]
    fn test_lammps_reader_load_and_n_atoms() {
        let mut r = PyLammpsReader::default();
        r.load_data(
            vec![vec![1.0, 1.0, 0.0, 0.5, 0.5]],
            vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0],
        );
        assert_eq!(r.n_atoms(), 1);
    }

    #[test]
    fn test_lammps_reader_box_bounds() {
        let mut r = PyLammpsReader::default();
        r.load_data(vec![], vec![-5.0, 5.0, -5.0, 5.0, -5.0, 5.0]);
        let b = r.box_bounds;
        assert_eq!(b[0], [-5.0, 5.0]);
    }

    #[test]
    fn test_lammps_reader_box_lengths() {
        let mut r = PyLammpsReader::default();
        r.load_data(vec![], vec![0.0, 10.0, 0.0, 20.0, 0.0, 30.0]);
        let l = r.box_lengths();
        assert_eq!(l, vec![10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_lammps_reader_read_atoms() {
        let mut r = PyLammpsReader::default();
        let atom = vec![1.0, 1.0, 0.1, 0.2, 0.3];
        r.load_data(vec![atom.clone()], vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0]);
        assert_eq!(r.read_atoms()[0], atom);
    }

    #[test]
    fn test_lammps_reader_default() {
        let r = PyLammpsReader::default();
        assert!(!r.filename.is_empty());
    }

    // --- PyHdf5Writer ---

    #[test]
    fn test_hdf5_writer_new() {
        let w = PyHdf5Writer::new("out.h5".to_string());
        assert_eq!(w.filename, "out.h5");
        assert_eq!(w.n_datasets(), 0);
    }

    #[test]
    fn test_hdf5_writer_write_dataset() {
        let mut w = PyHdf5Writer::default();
        w.write_dataset("pressure".to_string(), vec![1.0, 2.0, 3.0]);
        assert_eq!(w.n_datasets(), 1);
    }

    #[test]
    fn test_hdf5_writer_write_attribute() {
        let mut w = PyHdf5Writer::default();
        w.write_attribute("timestep".to_string(), 0.001);
        assert_eq!(w.n_attributes(), 1);
    }

    #[test]
    fn test_hdf5_writer_get_dataset() {
        let mut w = PyHdf5Writer::default();
        w.write_dataset("vel".to_string(), vec![1.0, 2.0]);
        let d = w.get_dataset("vel").expect("dataset exists");
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn test_hdf5_writer_get_attribute() {
        let mut w = PyHdf5Writer::default();
        w.write_attribute("dt".to_string(), 1e-4);
        let v = w.get_attribute("dt").expect("attribute exists");
        assert!((v - 1e-4).abs() < 1e-12);
    }

    #[test]
    fn test_hdf5_writer_missing_dataset_none() {
        let w = PyHdf5Writer::default();
        assert!(w.get_dataset("missing").is_none());
    }

    #[test]
    fn test_hdf5_writer_default() {
        let w = PyHdf5Writer::default();
        assert!(w.filename.ends_with(".h5"));
    }

    // --- PyTrajectoryWriter ---

    #[test]
    fn test_trajectory_writer_new() {
        let w = PyTrajectoryWriter::new("traj.xyz".to_string(), "xyz".to_string());
        assert_eq!(w.format, "xyz");
        assert_eq!(w.n_frames(), 0);
    }

    #[test]
    fn test_trajectory_writer_write_frame() {
        let mut w = PyTrajectoryWriter::default();
        w.write_frame(vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0], 0);
        assert_eq!(w.n_frames(), 1);
    }

    #[test]
    fn test_trajectory_writer_close() {
        let mut w = PyTrajectoryWriter::default();
        w.close();
        assert!(w.is_closed());
    }

    #[test]
    fn test_trajectory_writer_no_write_after_close() {
        let mut w = PyTrajectoryWriter::default();
        w.close();
        w.write_frame(vec![0.0], vec![], 1);
        assert_eq!(w.n_frames(), 0);
    }

    #[test]
    fn test_trajectory_writer_as_string_contains_step() {
        let mut w = PyTrajectoryWriter::default();
        w.write_frame(vec![1.0, 2.0, 3.0], vec![0.1, 0.2, 0.3], 42);
        assert!(w.as_string().contains("42"));
    }

    #[test]
    fn test_trajectory_writer_multiple_frames() {
        let mut w = PyTrajectoryWriter::new("t.lammps".to_string(), "lammps".to_string());
        for i in 0..10_u64 {
            w.write_frame(vec![0.0], vec![0.0], i);
        }
        assert_eq!(w.n_frames(), 10);
    }

    #[test]
    fn test_trajectory_writer_default() {
        let w = PyTrajectoryWriter::default();
        assert!(!w.format.is_empty());
    }

    // --- register_io_module is tested at the lib.rs integration level ---

    // --- N1: PyVtkWriter binary serialiser ---

    #[test]
    fn test_vtk_binary_header_starts_correctly() {
        let w = PyVtkWriter::new("test.vtk".to_string());
        let bytes = w.to_vtk_binary();
        let header = std::str::from_utf8(&bytes[..26]).expect("valid utf8 header");
        assert!(header.starts_with("# vtk DataFile"));
    }

    #[test]
    fn test_vtk_binary_length_gt_100() {
        let mut w = PyVtkWriter::new("test.vtk".to_string());
        // Add a position so the binary section exists; combined size must exceed 100.
        w.add_position(1.0, 2.0, 3.0);
        let bytes = w.to_vtk_binary();
        assert!(
            bytes.len() > 100,
            "expected > 100 bytes, got {}",
            bytes.len()
        );
    }

    #[test]
    fn test_vtk_binary_be_f32_roundtrip() {
        let mut w = PyVtkWriter::new("test.vtk".to_string());
        let original = [1.5_f64, -2.25_f64, 0.125_f64];
        w.add_position(original[0], original[1], original[2]);
        let bytes = w.to_vtk_binary();

        // The POINTS header ends with a newline; scan for the binary section
        // by locating "float\n" and reading the 12 bytes that follow.
        let marker = b"float\n";
        let start = bytes
            .windows(marker.len())
            .position(|w| w == marker)
            .expect("POINTS float header")
            + marker.len();
        let x = f32::from_be_bytes(bytes[start..start + 4].try_into().expect("4 bytes"));
        let y = f32::from_be_bytes(bytes[start + 4..start + 8].try_into().expect("4 bytes"));
        let z = f32::from_be_bytes(bytes[start + 8..start + 12].try_into().expect("4 bytes"));

        assert!((x as f64 - original[0]).abs() < 1e-6_f64);
        assert!((y as f64 - original[1]).abs() < 1e-6_f64);
        assert!((z as f64 - original[2]).abs() < 1e-6_f64);
    }

    #[test]
    fn test_vtk_write_to_file_creates_file() {
        let mut w = PyVtkWriter::new("test.vtk".to_string());
        w.add_position(1.0, 2.0, 3.0);
        let tmp = std::env::temp_dir().join("test_vtk_write.vtk");
        let path = tmp.to_str().expect("valid path");
        w.write_to_file(path).expect("write succeeds");
        let on_disk = std::fs::read(path).expect("file readable");
        let expected = w.to_vtk_binary();
        assert_eq!(on_disk, expected);
        let _ = std::fs::remove_file(path);
    }

    // --- N2: PyCsvWriter write_to_file ---

    #[test]
    fn test_csv_write_to_file_creates_readable_file() {
        let mut w = PyCsvWriter::new("out.csv".to_string());
        w.write_row(vec![1.0, 2.0, 3.0]);
        w.write_row(vec![4.0, 5.0, 6.0]);
        let tmp = std::env::temp_dir().join("test_csv_write.csv");
        let path = tmp.to_str().expect("valid path");
        w.write_to_file(path).expect("write succeeds");
        let on_disk = std::fs::read_to_string(path).expect("file readable");
        assert!(on_disk.contains("1") && on_disk.contains("6"));
        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn test_csv_write_to_file_matches_to_csv_string() {
        let mut w = PyCsvWriter::new("out.csv".to_string());
        w.write_row(vec![7.0, 8.0]);
        let tmp = std::env::temp_dir().join("test_csv_match.csv");
        let path = tmp.to_str().expect("valid path");
        w.write_to_file(path).expect("write succeeds");
        let on_disk = std::fs::read_to_string(path).expect("file readable");
        assert_eq!(on_disk, w.to_csv_string());
        let _ = std::fs::remove_file(path);
    }

    // --- N3: PyTrajectoryWriter format dispatch ---

    fn make_two_frame_trajectory(format: &str) -> PyTrajectoryWriter {
        let mut w = PyTrajectoryWriter::new("traj".to_string(), format.to_string());
        let atoms = vec![[0.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 0.5, 0.0]];
        let syms = vec!["C".to_string(), "H".to_string(), "O".to_string()];
        w.add_frame(TrajectoryFrame {
            step: 0,
            time: 0.0,
            positions: atoms.clone(),
            element_symbols: syms.clone(),
            box_lengths: Some([10.0, 10.0, 10.0]),
        });
        w.add_frame(TrajectoryFrame {
            step: 1,
            time: 0.001,
            positions: atoms,
            element_symbols: syms,
            box_lengths: None,
        });
        w
    }

    #[test]
    fn test_trajectory_xyz_atom_count_headers() {
        let w = make_two_frame_trajectory("xyz");
        let out = w.write_to_string().expect("xyz succeeds");
        // Each frame starts with the atom count (3) on its own line.
        let count = out.lines().filter(|l| *l == "3").count();
        assert_eq!(count, 2, "expected 2 atom-count lines, got: {}", count);
    }

    #[test]
    fn test_trajectory_lammps_timestep_per_frame() {
        let w = make_two_frame_trajectory("lammps");
        let out = w.write_to_string().expect("lammps succeeds");
        let count = out.matches("ITEM: TIMESTEP").count();
        assert_eq!(count, 2, "expected 2 ITEM: TIMESTEP, got {}", count);
    }

    #[test]
    fn test_trajectory_unknown_format_returns_err() {
        let w = PyTrajectoryWriter::new("traj".to_string(), "netcdf".to_string());
        let result = w.write_to_string();
        assert!(result.is_err());
    }

    #[test]
    fn test_trajectory_write_typed_to_file() {
        let w = make_two_frame_trajectory("xyz");
        let tmp = std::env::temp_dir().join("test_traj_typed.xyz");
        let path = tmp.to_str().expect("valid path");
        w.write_typed_to_file(path).expect("write succeeds");
        let on_disk = std::fs::read_to_string(path).expect("file readable");
        let expected = w.write_to_string().expect("write_to_string ok");
        assert_eq!(on_disk, expected);
        let _ = std::fs::remove_file(path);
    }
}
