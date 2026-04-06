// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! I/O API for Python interop.
//!
//! Provides lightweight in-memory representations of common simulation I/O
//! formats (VTK, CSV, XYZ, LAMMPS dump, HDF5-style, trajectory). All types
//! use plain `f64`, `Vec`f64`, and `String` — no nalgebra — for easy FFI.

#![allow(missing_docs)]
#![allow(dead_code)]

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// PyVtkWriter
// ---------------------------------------------------------------------------

/// In-memory VTK legacy-format writer (ASCII and binary stubs).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyVtkWriter {
    /// Output filename.
    pub filename: String,
    /// Named point data arrays.
    pub point_data: Vec<(String, Vec<f64>)>,
    /// Named cell data arrays.
    pub cell_data: Vec<(String, Vec<f64>)>,
}

impl PyVtkWriter {
    /// Create a new VTK writer targeting the given file.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            point_data: Vec::new(),
            cell_data: Vec::new(),
        }
    }

    /// Attach a named point-data array.
    pub fn add_point_data(&mut self, name: impl Into<String>, data: Vec<f64>) {
        self.point_data.push((name.into(), data));
    }

    /// Attach a named cell-data array.
    pub fn add_cell_data(&mut self, name: impl Into<String>, data: Vec<f64>) {
        self.cell_data.push((name.into(), data));
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

    /// Stub that returns the byte length that a binary VTK file would have.
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
        Self::new("output.vtk")
    }
}

// ---------------------------------------------------------------------------
// PyCsvReader
// ---------------------------------------------------------------------------

/// In-memory CSV reader.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCsvReader {
    /// Source filename.
    pub filename: String,
    /// Column header names (populated from the first row if `has_header = true`).
    pub headers: Vec<String>,
    /// Row-major data storage (rows × columns).
    pub rows: Vec<Vec<f64>>,
}

impl PyCsvReader {
    /// Create a new CSV reader backed by an in-memory dataset.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            headers: Vec::new(),
            rows: Vec::new(),
        }
    }

    /// Load in-memory data directly (useful in tests without touching disk).
    pub fn load_data(&mut self, headers: Vec<String>, rows: Vec<Vec<f64>>) {
        self.headers = headers;
        self.rows = rows;
    }

    /// Read a single column by index, returning a `Vec`f64`.
    pub fn read_column(&self, col: usize) -> Vec<f64> {
        self.rows
            .iter()
            .filter_map(|r| r.get(col).copied())
            .collect()
    }

    /// Return all data as a flat `Vec`f64` in row-major order.
    pub fn read_all_f64(&self) -> Vec<f64> {
        self.rows.iter().flat_map(|r| r.iter().copied()).collect()
    }

    /// Return the column header names.
    pub fn header_names(&self) -> &[String] {
        &self.headers
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
        Self::new("input.csv")
    }
}

// ---------------------------------------------------------------------------
// PyCsvWriter
// ---------------------------------------------------------------------------

/// In-memory CSV writer.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyCsvWriter {
    /// Output filename.
    pub filename: String,
    /// Buffered rows waiting to be flushed.
    pub buffer: Vec<Vec<f64>>,
}

impl PyCsvWriter {
    /// Create a new CSV writer targeting the given file.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            buffer: Vec::new(),
        }
    }

    /// Append a row of values to the internal buffer.
    pub fn write_row(&mut self, data: Vec<f64>) {
        self.buffer.push(data);
    }

    /// Consume the buffer and return it as a CSV string (stub — no disk I/O).
    pub fn flush(&mut self) -> String {
        let csv = self
            .buffer
            .iter()
            .map(|row| {
                row.iter()
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>()
            .join("\n");
        self.buffer.clear();
        csv
    }

    /// Number of rows currently buffered.
    pub fn buffered_rows(&self) -> usize {
        self.buffer.len()
    }
}

impl Default for PyCsvWriter {
    fn default() -> Self {
        Self::new("output.csv")
    }
}

// ---------------------------------------------------------------------------
// PyXyzReader
// ---------------------------------------------------------------------------

/// In-memory reader for the XYZ molecular dynamics format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyXyzReader {
    /// Source filename.
    pub filename: String,
    /// Atom species strings (e.g., `"C"`, `"H"`, `"O"`).
    pub species: Vec<String>,
    /// Flat positions: `\[x0, y0, z0, x1, y1, z1, …\]`.
    pub pos_flat: Vec<f64>,
}

impl PyXyzReader {
    /// Create a new XYZ reader.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
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

    /// Position array `\[x, y, z, …\]` for all atoms.
    pub fn positions(&self) -> &[f64] {
        &self.pos_flat
    }

    /// Species array.
    pub fn species(&self) -> &[String] {
        &self.species
    }

    /// Position `\[x, y, z\]` of atom `i`.
    pub fn position_of(&self, i: usize) -> Option<[f64; 3]> {
        let base = i * 3;
        if base + 2 < self.pos_flat.len() {
            Some([
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
        Self::new("atoms.xyz")
    }
}

// ---------------------------------------------------------------------------
// PyXyzWriter
// ---------------------------------------------------------------------------

/// In-memory writer for the XYZ molecular dynamics format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyXyzWriter {
    /// Output filename.
    pub filename: String,
    /// Accumulated frames as raw XYZ strings.
    pub frames: Vec<String>,
}

impl PyXyzWriter {
    /// Create a new XYZ writer.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            frames: Vec::new(),
        }
    }

    /// Append a frame to the internal buffer.
    ///
    /// `positions` is flat `\[x0,y0,z0, x1,y1,z1, …\]`; `species` has one entry
    /// per atom; `comment` is written to the second header line.
    pub fn write_frame(
        &mut self,
        positions: &[f64],
        species: &[String],
        comment: impl Into<String>,
    ) {
        let n = species.len();
        let mut frame = format!("{}\n{}\n", n, comment.into());
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
        Self::new("output.xyz")
    }
}

// ---------------------------------------------------------------------------
// PyLammpsReader
// ---------------------------------------------------------------------------

/// In-memory reader for LAMMPS dump files.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyLammpsReader {
    /// Source filename.
    pub filename: String,
    /// Flat atom data: `\[id, type, x, y, z, …\]` per atom (5 fields each).
    pub atom_data: Vec<Vec<f64>>,
    /// Simulation box bounds `\[\[xlo, xhi\\], \[ylo, yhi\], \[zlo, zhi\]]`.
    pub box_bounds: [[f64; 2]; 3],
}

impl PyLammpsReader {
    /// Create a new LAMMPS reader.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            atom_data: Vec::new(),
            box_bounds: [[0.0, 1.0]; 3],
        }
    }

    /// Load in-memory atom data and box bounds.
    pub fn load_data(&mut self, atom_data: Vec<Vec<f64>>, box_bounds: [[f64; 2]; 3]) {
        self.atom_data = atom_data;
        self.box_bounds = box_bounds;
    }

    /// Return a reference to all atom records.
    pub fn read_atoms(&self) -> &[Vec<f64>] {
        &self.atom_data
    }

    /// Number of atoms in the last read frame.
    pub fn n_atoms(&self) -> usize {
        self.atom_data.len()
    }

    /// Box bounds `\[\[xlo,xhi\\],\[ylo,yhi\],\[zlo,zhi\]]`.
    pub fn box_bounds(&self) -> [[f64; 2]; 3] {
        self.box_bounds
    }

    /// Box side lengths `\[Lx, Ly, Lz\]`.
    pub fn box_lengths(&self) -> [f64; 3] {
        [
            self.box_bounds[0][1] - self.box_bounds[0][0],
            self.box_bounds[1][1] - self.box_bounds[1][0],
            self.box_bounds[2][1] - self.box_bounds[2][0],
        ]
    }
}

impl Default for PyLammpsReader {
    fn default() -> Self {
        Self::new("dump.lammpstrj")
    }
}

// ---------------------------------------------------------------------------
// PyHdf5Writer
// ---------------------------------------------------------------------------

/// In-memory HDF5-style writer (no actual HDF5 dependency).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyHdf5Writer {
    /// Output filename.
    pub filename: String,
    /// Named datasets.
    pub datasets: Vec<(String, Vec<f64>)>,
    /// Named scalar attributes.
    pub attributes: Vec<(String, f64)>,
}

impl PyHdf5Writer {
    /// Create a new HDF5 writer.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            datasets: Vec::new(),
            attributes: Vec::new(),
        }
    }

    /// Write a named dataset.
    pub fn write_dataset(&mut self, name: impl Into<String>, data: Vec<f64>) {
        self.datasets.push((name.into(), data));
    }

    /// Write a named scalar attribute.
    pub fn write_attribute(&mut self, name: impl Into<String>, value: f64) {
        self.attributes.push((name.into(), value));
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
    pub fn get_dataset(&self, name: &str) -> Option<&Vec<f64>> {
        self.datasets
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, d)| d)
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
        Self::new("output.h5")
    }
}

// ---------------------------------------------------------------------------
// PyTrajectoryWriter
// ---------------------------------------------------------------------------

/// Multi-format trajectory writer (in-memory stub).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PyTrajectoryWriter {
    /// Output filename.
    pub filename: String,
    /// Format string (e.g., `"xyz"`, `"lammps"`, `"vtk"`).
    pub format: String,
    /// Whether the writer has been closed.
    pub closed: bool,
    /// Number of frames written.
    pub frame_count: usize,
    /// Buffered frame strings.
    pub frame_buffer: Vec<String>,
}

impl PyTrajectoryWriter {
    /// Create a new trajectory writer.
    pub fn new(filename: impl Into<String>, format: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
            format: format.into(),
            closed: false,
            frame_count: 0,
            frame_buffer: Vec::new(),
        }
    }

    /// Write a trajectory frame.
    ///
    /// `positions` and `velocities` are flat `\[x,y,z, …\]` arrays; `step` is
    /// the integer simulation step number.
    pub fn write_frame(&mut self, positions: &[f64], velocities: &[f64], step: u64) {
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

impl Default for PyTrajectoryWriter {
    fn default() -> Self {
        Self::new("trajectory.xyz", "xyz")
    }
}

// ---------------------------------------------------------------------------
// Registration helper
// ---------------------------------------------------------------------------

/// Register all I/O classes into a Python sub-module named `"io"`.
///
/// This is a no-op placeholder that documents the intended PyO3 registration
/// point. When PyO3 is enabled as a dependency the body should call
/// `m.add_class::`PyVtkWriter`()` etc.
pub fn register_io_module(_m: &str) {
    // Placeholder: actual PyO3 registration would happen here.
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
        let w = PyVtkWriter::new("out.vtk");
        assert_eq!(w.filename, "out.vtk");
        assert_eq!(w.n_point_arrays(), 0);
    }

    #[test]
    fn test_vtk_add_point_data() {
        let mut w = PyVtkWriter::default();
        w.add_point_data("pressure", vec![1.0, 2.0, 3.0]);
        assert_eq!(w.n_point_arrays(), 1);
    }

    #[test]
    fn test_vtk_add_cell_data() {
        let mut w = PyVtkWriter::default();
        w.add_cell_data("stress", vec![10.0, 20.0]);
        assert_eq!(w.n_cell_arrays(), 1);
    }

    #[test]
    fn test_vtk_write_ascii_contains_header() {
        let w = PyVtkWriter::new("test.vtk");
        let s = w.write_ascii();
        assert!(s.contains("vtk DataFile"));
    }

    #[test]
    fn test_vtk_write_ascii_contains_point_data_name() {
        let mut w = PyVtkWriter::new("test.vtk");
        w.add_point_data("velocity", vec![1.0, 2.0]);
        let s = w.write_ascii();
        assert!(s.contains("velocity"));
    }

    #[test]
    fn test_vtk_write_binary_size_grows() {
        let mut w = PyVtkWriter::new("test.vtk");
        let s0 = w.write_binary();
        w.add_point_data("p", vec![1.0; 100]);
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
        let r = PyCsvReader::new("data.csv");
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
        let w = PyCsvWriter::new("out.csv");
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
        let r = PyXyzReader::new("mol.xyz");
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
        assert_eq!(r.positions(), &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_xyz_reader_species() {
        let mut r = PyXyzReader::default();
        r.load_data(vec!["N".to_string()], vec![0.0, 0.0, 0.0]);
        assert_eq!(r.species()[0], "N");
    }

    #[test]
    fn test_xyz_reader_position_of() {
        let mut r = PyXyzReader::default();
        r.load_data(
            vec!["C".to_string(), "H".to_string()],
            vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        );
        let p = r.position_of(0).unwrap();
        assert_eq!(p, [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_xyz_reader_default() {
        let r = PyXyzReader::default();
        assert!(r.filename.ends_with(".xyz"));
    }

    // --- PyXyzWriter ---

    #[test]
    fn test_xyz_writer_new() {
        let w = PyXyzWriter::new("out.xyz");
        assert_eq!(w.filename, "out.xyz");
        assert_eq!(w.n_frames(), 0);
    }

    #[test]
    fn test_xyz_writer_write_frame_increments_count() {
        let mut w = PyXyzWriter::default();
        w.write_frame(&[0.0, 0.0, 0.0], &["C".to_string()], "frame 0");
        assert_eq!(w.n_frames(), 1);
    }

    #[test]
    fn test_xyz_writer_as_string_contains_n_atoms() {
        let mut w = PyXyzWriter::default();
        w.write_frame(&[0.0, 0.0, 0.0], &["C".to_string()], "test");
        assert!(w.as_string().contains('1'));
    }

    #[test]
    fn test_xyz_writer_multiple_frames() {
        let mut w = PyXyzWriter::default();
        for _ in 0..5 {
            w.write_frame(&[0.0, 0.0, 0.0], &["H".to_string()], "");
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
        let r = PyLammpsReader::new("dump.lammps");
        assert_eq!(r.n_atoms(), 0);
    }

    #[test]
    fn test_lammps_reader_load_and_n_atoms() {
        let mut r = PyLammpsReader::default();
        r.load_data(vec![vec![1.0, 1.0, 0.0, 0.5, 0.5]], [[0.0, 1.0]; 3]);
        assert_eq!(r.n_atoms(), 1);
    }

    #[test]
    fn test_lammps_reader_box_bounds() {
        let mut r = PyLammpsReader::default();
        r.load_data(vec![], [[-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0]]);
        let b = r.box_bounds();
        assert_eq!(b[0], [-5.0, 5.0]);
    }

    #[test]
    fn test_lammps_reader_box_lengths() {
        let mut r = PyLammpsReader::default();
        r.load_data(vec![], [[0.0, 10.0], [0.0, 20.0], [0.0, 30.0]]);
        let l = r.box_lengths();
        assert_eq!(l, [10.0, 20.0, 30.0]);
    }

    #[test]
    fn test_lammps_reader_read_atoms() {
        let mut r = PyLammpsReader::default();
        let atom = vec![1.0, 1.0, 0.1, 0.2, 0.3];
        r.load_data(vec![atom.clone()], [[0.0, 1.0]; 3]);
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
        let w = PyHdf5Writer::new("out.h5");
        assert_eq!(w.filename, "out.h5");
        assert_eq!(w.n_datasets(), 0);
    }

    #[test]
    fn test_hdf5_writer_write_dataset() {
        let mut w = PyHdf5Writer::default();
        w.write_dataset("pressure", vec![1.0, 2.0, 3.0]);
        assert_eq!(w.n_datasets(), 1);
    }

    #[test]
    fn test_hdf5_writer_write_attribute() {
        let mut w = PyHdf5Writer::default();
        w.write_attribute("timestep", 0.001);
        assert_eq!(w.n_attributes(), 1);
    }

    #[test]
    fn test_hdf5_writer_get_dataset() {
        let mut w = PyHdf5Writer::default();
        w.write_dataset("vel", vec![1.0, 2.0]);
        let d = w.get_dataset("vel").unwrap();
        assert_eq!(d.len(), 2);
    }

    #[test]
    fn test_hdf5_writer_get_attribute() {
        let mut w = PyHdf5Writer::default();
        w.write_attribute("dt", 1e-4);
        let v = w.get_attribute("dt").unwrap();
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
        let w = PyTrajectoryWriter::new("traj.xyz", "xyz");
        assert_eq!(w.format, "xyz");
        assert_eq!(w.n_frames(), 0);
    }

    #[test]
    fn test_trajectory_writer_write_frame() {
        let mut w = PyTrajectoryWriter::default();
        w.write_frame(&[0.0, 0.0, 0.0], &[0.0, 0.0, 0.0], 0);
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
        w.write_frame(&[0.0], &[], 1);
        assert_eq!(w.n_frames(), 0);
    }

    #[test]
    fn test_trajectory_writer_as_string_contains_step() {
        let mut w = PyTrajectoryWriter::default();
        w.write_frame(&[1.0, 2.0, 3.0], &[0.1, 0.2, 0.3], 42);
        assert!(w.as_string().contains("42"));
    }

    #[test]
    fn test_trajectory_writer_multiple_frames() {
        let mut w = PyTrajectoryWriter::new("t.lammps", "lammps");
        for i in 0..10_u64 {
            w.write_frame(&[0.0], &[0.0], i);
        }
        assert_eq!(w.n_frames(), 10);
    }

    #[test]
    fn test_trajectory_writer_default() {
        let w = PyTrajectoryWriter::default();
        assert!(!w.format.is_empty());
    }

    // --- register_io_module ---

    #[test]
    fn test_register_io_module_no_panic() {
        register_io_module("io");
    }
}
