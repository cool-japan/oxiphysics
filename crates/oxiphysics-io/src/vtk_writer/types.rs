//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::io::Write;
/// Writer for the VTK XML Unstructured Grid format (`.vtu`).
pub struct XmlVtuWriter;
impl XmlVtuWriter {
    /// Serialise a complete unstructured grid to a VTU XML string.
    pub fn write(grid: &VtkUnstructuredGrid) -> String {
        let mut s = String::new();
        s.push_str(&Self::write_xml_header());
        s.push_str(&format!(
            "  <UnstructuredGrid>\n    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            grid.n_points(),
            grid.n_cells()
        ));
        s.push_str("      <Points>\n");
        s.push_str(&Self::write_points_xml(&grid.points));
        s.push_str("      </Points>\n");
        s.push_str("      <Cells>\n");
        s.push_str(&Self::write_cells_xml(&grid.cells, &grid.cell_types));
        s.push_str("      </Cells>\n");
        if !grid.point_data.is_empty() {
            s.push_str("      <PointData>\n");
            for arr in &grid.point_data {
                match arr {
                    VtkDataArrayW::Scalars { name, values } => {
                        let flat: Vec<f64> = values.clone();
                        s.push_str(&Self::write_data_array_xml(name, 1, &flat));
                    }
                    VtkDataArrayW::Vectors { name, values } => {
                        let flat: Vec<f64> =
                            values.iter().flat_map(|v| v.iter().copied()).collect();
                        s.push_str(&Self::write_data_array_xml(name, 3, &flat));
                    }
                    VtkDataArrayW::Tensors { name, values } => {
                        let flat: Vec<f64> = values
                            .iter()
                            .flat_map(|t| t.iter().flat_map(|r| r.iter().copied()))
                            .collect();
                        s.push_str(&Self::write_data_array_xml(name, 9, &flat));
                    }
                }
            }
            s.push_str("      </PointData>\n");
        }
        if !grid.cell_data.is_empty() {
            s.push_str("      <CellData>\n");
            for arr in &grid.cell_data {
                match arr {
                    VtkDataArrayW::Scalars { name, values } => {
                        let flat: Vec<f64> = values.clone();
                        s.push_str(&Self::write_data_array_xml(name, 1, &flat));
                    }
                    VtkDataArrayW::Vectors { name, values } => {
                        let flat: Vec<f64> =
                            values.iter().flat_map(|v| v.iter().copied()).collect();
                        s.push_str(&Self::write_data_array_xml(name, 3, &flat));
                    }
                    VtkDataArrayW::Tensors { name, values } => {
                        let flat: Vec<f64> = values
                            .iter()
                            .flat_map(|t| t.iter().flat_map(|r| r.iter().copied()))
                            .collect();
                        s.push_str(&Self::write_data_array_xml(name, 9, &flat));
                    }
                }
            }
            s.push_str("      </CellData>\n");
        }
        s.push_str("    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n");
        s
    }
    /// Generate the VTU XML file header.
    pub fn write_xml_header() -> String {
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            .to_owned()
    }
    /// Serialise the Points DataArray element.
    pub fn write_points_xml(points: &[[f64; 3]]) -> String {
        let mut s = String::from(
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        );
        for p in points {
            s.push_str(&format!("          {} {} {}\n", p[0], p[1], p[2]));
        }
        s.push_str("        </DataArray>\n");
        s
    }
    /// Serialise the Cells section (connectivity, offsets, types).
    pub fn write_cells_xml(cells: &[Vec<usize>], cell_types: &[VtkCellTypeW]) -> String {
        let mut s = String::new();
        s.push_str(
            "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ",
        );
        let mut first = true;
        for conn in cells {
            for &idx in conn {
                if !first {
                    s.push(' ');
                }
                s.push_str(&idx.to_string());
                first = false;
            }
        }
        s.push_str("\n        </DataArray>\n");
        s.push_str(
            "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ",
        );
        let mut offset: usize = 0;
        for (i, conn) in cells.iter().enumerate() {
            if i > 0 {
                s.push(' ');
            }
            offset += conn.len();
            s.push_str(&offset.to_string());
        }
        s.push_str("\n        </DataArray>\n");
        s.push_str(
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ",
        );
        for (i, ct) in cell_types.iter().enumerate() {
            if i > 0 {
                s.push(' ');
            }
            s.push_str(&(*ct as u8).to_string());
        }
        s.push_str("\n        </DataArray>\n");
        s
    }
    /// Serialise a generic DataArray XML element.
    pub fn write_data_array_xml(name: &str, n_components: usize, values: &[f64]) -> String {
        let mut s = format!(
            "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"{}\" format=\"ascii\">\n          ",
            name, n_components
        );
        for (i, v) in values.iter().enumerate() {
            if i > 0 {
                s.push(' ');
            }
            s.push_str(&v.to_string());
        }
        s.push_str("\n        </DataArray>\n");
        s
    }
    /// Encode a byte slice as a Base64 string (no line wrapping).
    pub fn encode_base64(data: &[u8]) -> String {
        const TABLE: &[u8; 64] =
            b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
        for chunk in data.chunks(3) {
            let b0 = chunk[0] as u32;
            let b1 = if chunk.len() > 1 { chunk[1] as u32 } else { 0 };
            let b2 = if chunk.len() > 2 { chunk[2] as u32 } else { 0 };
            let triple = (b0 << 16) | (b1 << 8) | b2;
            out.push(TABLE[((triple >> 18) & 0x3F) as usize] as char);
            out.push(TABLE[((triple >> 12) & 0x3F) as usize] as char);
            if chunk.len() > 1 {
                out.push(TABLE[((triple >> 6) & 0x3F) as usize] as char);
            } else {
                out.push('=');
            }
            if chunk.len() > 2 {
                out.push(TABLE[(triple & 0x3F) as usize] as char);
            } else {
                out.push('=');
            }
        }
        out
    }
}
/// Writer for the VTK legacy ASCII format (`.vtk`).
pub struct LegacyVtkWriter;
impl LegacyVtkWriter {
    /// Serialise a complete unstructured grid to a legacy VTK ASCII string.
    pub fn write(grid: &VtkUnstructuredGrid) -> String {
        let mut s = String::new();
        s.push_str(&Self::write_header());
        s.push_str(&Self::write_points(&grid.points));
        s.push_str(&Self::write_cells(&grid.cells, &grid.cell_types));
        if !grid.point_data.is_empty() {
            s.push_str(&Self::write_point_data(&grid.point_data));
        }
        if !grid.cell_data.is_empty() {
            s.push_str(&Self::write_cell_data(&grid.cell_data));
        }
        s
    }
    /// Generate the standard VTK legacy file header.
    pub fn write_header() -> String {
        "# vtk DataFile Version 3.0\nOxiPhysics unstructured grid\nASCII\nDATASET UNSTRUCTURED_GRID\n"
            .to_owned()
    }
    /// Serialise the POINTS section.
    pub fn write_points(points: &[[f64; 3]]) -> String {
        let mut s = format!("POINTS {} double\n", points.len());
        for p in points {
            s.push_str(&format!("{} {} {}\n", p[0], p[1], p[2]));
        }
        s
    }
    /// Serialise the CELLS and CELL_TYPES sections.
    pub fn write_cells(cells: &[Vec<usize>], cell_types: &[VtkCellTypeW]) -> String {
        let total_entries: usize = cells.iter().map(|c| c.len() + 1).sum();
        let mut s = format!("CELLS {} {}\n", cells.len(), total_entries);
        for conn in cells {
            s.push_str(&conn.len().to_string());
            for &idx in conn {
                s.push(' ');
                s.push_str(&idx.to_string());
            }
            s.push('\n');
        }
        s.push_str(&format!("CELL_TYPES {}\n", cell_types.len()));
        for ct in cell_types {
            s.push_str(&format!("{}\n", *ct as u8));
        }
        s
    }
    /// Serialise the POINT_DATA section.
    pub fn write_point_data(arrays: &[VtkDataArrayW]) -> String {
        if arrays.is_empty() {
            return String::new();
        }
        let n = arrays[0].len();
        let mut s = format!("POINT_DATA {}\n", n);
        for arr in arrays {
            s.push_str(&Self::write_data_array(arr));
        }
        s
    }
    /// Serialise the CELL_DATA section.
    pub fn write_cell_data(arrays: &[VtkDataArrayW]) -> String {
        if arrays.is_empty() {
            return String::new();
        }
        let n = arrays[0].len();
        let mut s = format!("CELL_DATA {}\n", n);
        for arr in arrays {
            s.push_str(&Self::write_data_array(arr));
        }
        s
    }
    /// Serialise a single data array in legacy VTK format.
    pub fn write_data_array(arr: &VtkDataArrayW) -> String {
        let mut s = String::new();
        match arr {
            VtkDataArrayW::Scalars { name, values } => {
                s.push_str(&format!(
                    "SCALARS {} double 1\nLOOKUP_TABLE default\n",
                    name
                ));
                for v in values {
                    s.push_str(&format!("{}\n", v));
                }
            }
            VtkDataArrayW::Vectors { name, values } => {
                s.push_str(&format!("VECTORS {} double\n", name));
                for v in values {
                    s.push_str(&format!("{} {} {}\n", v[0], v[1], v[2]));
                }
            }
            VtkDataArrayW::Tensors { name, values } => {
                s.push_str(&format!("TENSORS {} double\n", name));
                for t in values {
                    for row in t {
                        s.push_str(&format!("{} {} {}\n", row[0], row[1], row[2]));
                    }
                    s.push('\n');
                }
            }
        }
        s
    }
}
/// Time-step writer: manages sequential VTU output with a PVD index.
pub struct TimeStepWriter {
    /// Base directory for output files.
    pub output_dir: String,
    /// Base name for each VTU file (without extension).
    pub base_name: String,
    /// Accumulated PVD entries.
    pub entries: Vec<(f64, String)>,
}
impl TimeStepWriter {
    /// Create a new time-step writer.
    pub fn new(output_dir: impl Into<String>, base_name: impl Into<String>) -> Self {
        Self {
            output_dir: output_dir.into(),
            base_name: base_name.into(),
            entries: Vec::new(),
        }
    }
    /// Get the VTU filename for step `i`.
    pub fn vtu_filename(&self, step: usize) -> String {
        format!("{}_{:06}.vtu", self.base_name, step)
    }
    /// Full path to the VTU file for step `i`.
    pub fn vtu_path(&self, step: usize) -> String {
        format!("{}/{}", self.output_dir, self.vtu_filename(step))
    }
    /// Register a time step (without writing).
    pub fn register_step(&mut self, time: f64, step: usize) {
        self.entries.push((time, self.vtu_filename(step)));
    }
    /// Write the PVD collection file.
    pub fn write_pvd(&self, pvd_name: &str) -> std::io::Result<()> {
        use std::io::Write;
        let path = format!("{}/{}", self.output_dir, pvd_name);
        let file = std::fs::File::create(&path)?;
        let mut w = std::io::BufWriter::new(file);
        writeln!(w, r#"<?xml version="1.0"?>"#)?;
        writeln!(
            w,
            r#"<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">"#
        )?;
        writeln!(w, r#"  <Collection>"#)?;
        for (t, fname) in &self.entries {
            writeln!(
                w,
                r#"    <DataSet timestep="{t:.6e}" part="0" file="{fname}"/>"#
            )?;
        }
        writeln!(w, r#"  </Collection>"#)?;
        writeln!(w, r#"</VTKFile>"#)?;
        w.flush()?;
        Ok(())
    }
    /// Number of registered steps.
    pub fn n_steps(&self) -> usize {
        self.entries.len()
    }
}
/// A VTK POLYDATA dataset for surface meshes (lines and polygons).
pub struct VtkPolyData {
    /// 3-D coordinates of every point.
    pub points: Vec<[f64; 3]>,
    /// Line connectivity lists.
    pub lines: Vec<Vec<usize>>,
    /// Polygon connectivity lists.
    pub polygons: Vec<Vec<usize>>,
    /// Per-point data arrays.
    pub point_data: Vec<VtkDataArrayW>,
}
impl VtkPolyData {
    /// Create an empty poly data set.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            lines: Vec::new(),
            polygons: Vec::new(),
            point_data: Vec::new(),
        }
    }
    /// Serialise to VTK legacy ASCII POLYDATA format.
    pub fn write_legacy(&self) -> String {
        let mut s = String::new();
        s.push_str("# vtk DataFile Version 3.0\nOxiPhysics polydata\nASCII\nDATASET POLYDATA\n");
        s.push_str(&format!("POINTS {} double\n", self.points.len()));
        for p in &self.points {
            s.push_str(&format!("{} {} {}\n", p[0], p[1], p[2]));
        }
        if !self.lines.is_empty() {
            let total: usize = self.lines.iter().map(|l| l.len() + 1).sum();
            s.push_str(&format!("LINES {} {}\n", self.lines.len(), total));
            for line in &self.lines {
                s.push_str(&line.len().to_string());
                for &idx in line {
                    s.push(' ');
                    s.push_str(&idx.to_string());
                }
                s.push('\n');
            }
        }
        if !self.polygons.is_empty() {
            let total: usize = self.polygons.iter().map(|p| p.len() + 1).sum();
            s.push_str(&format!("POLYGONS {} {}\n", self.polygons.len(), total));
            for poly in &self.polygons {
                s.push_str(&poly.len().to_string());
                for &idx in poly {
                    s.push(' ');
                    s.push_str(&idx.to_string());
                }
                s.push('\n');
            }
        }
        if !self.point_data.is_empty() {
            s.push_str(&format!("POINT_DATA {}\n", self.points.len()));
            for arr in &self.point_data {
                s.push_str(&LegacyVtkWriter::write_data_array(arr));
            }
        }
        s
    }
}
/// An unstructured grid for VTK output.
///
/// Build up the grid with [`add_point`](VtkUnstructuredGrid::add_point) and
/// [`add_cell`](VtkUnstructuredGrid::add_cell), attach field data with the
/// `add_point_*` / `add_cell_*` helpers, then serialise with
/// [`LegacyVtkWriter::write`] or [`XmlVtuWriter::write`].
pub struct VtkUnstructuredGrid {
    /// 3-D coordinates of every point.
    pub points: Vec<[f64; 3]>,
    /// Cell connectivity (variable-length point-index lists).
    pub cells: Vec<Vec<usize>>,
    /// VTK cell type for each cell (parallel to `cells`).
    pub cell_types: Vec<VtkCellTypeW>,
    /// Per-vertex data arrays.
    pub point_data: Vec<VtkDataArrayW>,
    /// Per-cell data arrays.
    pub cell_data: Vec<VtkDataArrayW>,
}
impl VtkUnstructuredGrid {
    /// Create an empty unstructured grid.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            cells: Vec::new(),
            cell_types: Vec::new(),
            point_data: Vec::new(),
            cell_data: Vec::new(),
        }
    }
    /// Append a point and return its zero-based index.
    pub fn add_point(&mut self, p: [f64; 3]) -> usize {
        let idx = self.points.len();
        self.points.push(p);
        idx
    }
    /// Append a cell defined by `connectivity` (point indices) and a `cell_type`.
    pub fn add_cell(&mut self, connectivity: Vec<usize>, cell_type: VtkCellTypeW) {
        self.cells.push(connectivity);
        self.cell_types.push(cell_type);
    }
    /// Attach a per-point scalar field.
    pub fn add_point_scalars(&mut self, name: &str, values: Vec<f64>) {
        self.point_data.push(VtkDataArrayW::Scalars {
            name: name.to_owned(),
            values,
        });
    }
    /// Attach a per-point vector field.
    pub fn add_point_vectors(&mut self, name: &str, values: Vec<[f64; 3]>) {
        self.point_data.push(VtkDataArrayW::Vectors {
            name: name.to_owned(),
            values,
        });
    }
    /// Attach a per-cell scalar field.
    pub fn add_cell_scalars(&mut self, name: &str, values: Vec<f64>) {
        self.cell_data.push(VtkDataArrayW::Scalars {
            name: name.to_owned(),
            values,
        });
    }
    /// Returns the number of points in the grid.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
    /// Returns the number of cells in the grid.
    pub fn n_cells(&self) -> usize {
        self.cells.len()
    }
}
/// Partition-information for one MPI rank's data in a parallel VTK output.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VtkPartition {
    /// Rank index.
    pub rank: usize,
    /// File name for this piece.
    pub filename: String,
    /// Global point offset.
    pub point_offset: usize,
    /// Number of points in this piece.
    pub n_points: usize,
    /// Number of cells in this piece.
    pub n_cells: usize,
}
/// VTK cell type identifiers used by the vtk_writer module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VtkCellTypeW {
    /// VTK_VERTEX (1)
    Vertex = 1,
    /// VTK_LINE (3)
    Line = 3,
    /// VTK_TRIANGLE (5)
    Triangle = 5,
    /// VTK_QUAD (9)
    Quad = 9,
    /// VTK_TETRA (10)
    Tetrahedron = 10,
    /// VTK_HEXAHEDRON (12)
    Hexahedron = 12,
    /// VTK_WEDGE (13)
    Wedge = 13,
    /// VTK_PYRAMID (14)
    Pyramid = 14,
}
impl VtkCellTypeW {
    /// Returns the number of points per cell for this cell type.
    pub fn n_points(&self) -> usize {
        match self {
            Self::Vertex => 1,
            Self::Line => 2,
            Self::Triangle => 3,
            Self::Quad => 4,
            Self::Tetrahedron => 4,
            Self::Hexahedron => 8,
            Self::Wedge => 6,
            Self::Pyramid => 5,
        }
    }
}
/// Improved ASCII writer with configurable precision.
pub struct FormattedAsciiWriter {
    /// Number of decimal places for floating-point output.
    pub precision: usize,
}
impl FormattedAsciiWriter {
    /// Create with default precision (6 decimal places).
    pub fn new() -> Self {
        Self { precision: 6 }
    }
    /// Create with specified precision.
    pub fn with_precision(precision: usize) -> Self {
        Self { precision }
    }
    /// Format a single f64 value.
    pub fn format_f64(&self, v: f64) -> String {
        format!("{:.prec$}", v, prec = self.precision)
    }
    /// Format a 3D point.
    pub fn format_point(&self, p: [f64; 3]) -> String {
        format!(
            "{:.prec$} {:.prec$} {:.prec$}",
            p[0],
            p[1],
            p[2],
            prec = self.precision
        )
    }
    /// Write POINTS section with specified precision.
    pub fn write_points(&self, points: &[[f64; 3]]) -> String {
        let mut s = format!("POINTS {} double\n", points.len());
        for p in points {
            s.push_str(&self.format_point(*p));
            s.push('\n');
        }
        s
    }
    /// Write a scalar data section with specified precision.
    pub fn write_scalars(&self, name: &str, values: &[f64]) -> String {
        let mut s = format!("SCALARS {} double 1\nLOOKUP_TABLE default\n", name);
        for v in values {
            s.push_str(&self.format_f64(*v));
            s.push('\n');
        }
        s
    }
    /// Write a full grid using the legacy ASCII format with configurable precision.
    pub fn write_grid(&self, grid: &VtkUnstructuredGrid) -> String {
        let mut s = LegacyVtkWriter::write_header();
        s.push_str(&self.write_points(&grid.points));
        s.push_str(&LegacyVtkWriter::write_cells(&grid.cells, &grid.cell_types));
        if !grid.point_data.is_empty() {
            let n = grid.point_data[0].len();
            s.push_str(&format!("POINT_DATA {}\n", n));
            for arr in &grid.point_data {
                if let VtkDataArrayW::Scalars { name, values } = arr {
                    s.push_str(&self.write_scalars(name, values));
                } else {
                    s.push_str(&LegacyVtkWriter::write_data_array(arr));
                }
            }
        }
        s
    }
}
/// Validation result for a VTU grid.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VtuValidationResult {
    /// Whether the file structure is valid.
    pub is_valid: bool,
    /// List of issues found.
    pub issues: Vec<String>,
}
impl VtuValidationResult {
    pub(super) fn ok() -> Self {
        Self {
            is_valid: true,
            issues: Vec::new(),
        }
    }
    pub(super) fn add_issue(&mut self, msg: impl Into<String>) {
        self.is_valid = false;
        self.issues.push(msg.into());
    }
}
/// Provides simple run-length encoding for f64 data (for illustration).
pub struct VtkCompression;
impl VtkCompression {
    /// Run-length encode a slice of `f64` values.
    ///
    /// Consecutive equal values are stored as `(value, count)` pairs.
    pub fn rle_encode_f64(data: &[f64]) -> Vec<(f64, usize)> {
        if data.is_empty() {
            return vec![];
        }
        let mut runs: Vec<(f64, usize)> = Vec::new();
        let mut cur = data[0];
        let mut count = 1usize;
        for &v in &data[1..] {
            if (v - cur).abs() < 1e-15 {
                count += 1;
            } else {
                runs.push((cur, count));
                cur = v;
                count = 1;
            }
        }
        runs.push((cur, count));
        runs
    }
    /// Decode RLE-encoded f64 data back to a flat vector.
    pub fn rle_decode_f64(runs: &[(f64, usize)]) -> Vec<f64> {
        let total: usize = runs.iter().map(|(_, c)| c).sum();
        let mut out = Vec::with_capacity(total);
        for &(v, c) in runs {
            for _ in 0..c {
                out.push(v);
            }
        }
        out
    }
    /// Compute compression ratio: `original_len / encoded_pairs`.
    pub fn compression_ratio(original_len: usize, n_runs: usize) -> f64 {
        if n_runs == 0 {
            return 0.0;
        }
        original_len as f64 / n_runs as f64
    }
    /// Losslessly encode a u8 slice using a simple byte-pair delta encoding.
    ///
    /// Returns `(delta_bytes, first_value)`. Decode with `delta_decode_u8`.
    pub fn delta_encode_u8(data: &[u8]) -> (Vec<i8>, u8) {
        if data.is_empty() {
            return (vec![], 0);
        }
        let first = data[0];
        let mut deltas: Vec<i8> = Vec::with_capacity(data.len() - 1);
        for w in data.windows(2) {
            deltas.push(w[1].wrapping_sub(w[0]) as i8);
        }
        (deltas, first)
    }
    /// Decode delta-encoded bytes.
    pub fn delta_decode_u8(deltas: &[i8], first: u8) -> Vec<u8> {
        let mut out = Vec::with_capacity(deltas.len() + 1);
        out.push(first);
        let mut cur = first;
        for &d in deltas {
            cur = cur.wrapping_add(d as u8);
            out.push(cur);
        }
        out
    }
}
/// Exports particle data as VTK VERTEX cells.
pub struct ParticleVtkExporter;
impl ParticleVtkExporter {
    /// Export a particle set (each particle as a VTK VERTEX cell).
    ///
    /// Optional velocity vectors and a scalar field can be attached.
    pub fn export_particles(
        positions: &[[f64; 3]],
        velocities: Option<&[[f64; 3]]>,
        scalars: Option<(&str, &[f64])>,
    ) -> String {
        let mut grid = VtkUnstructuredGrid::new();
        for &p in positions {
            let idx = grid.add_point(p);
            grid.add_cell(vec![idx], VtkCellTypeW::Vertex);
        }
        if let Some(vels) = velocities {
            grid.add_point_vectors("velocity", vels.to_vec());
        }
        if let Some((name, vals)) = scalars {
            grid.add_point_scalars(name, vals.to_vec());
        }
        LegacyVtkWriter::write(&grid)
    }
    /// Export SPH particles with density and pressure fields.
    pub fn export_sph_particles(
        positions: &[[f64; 3]],
        densities: &[f64],
        pressures: &[f64],
    ) -> String {
        let mut grid = VtkUnstructuredGrid::new();
        for &p in positions {
            let idx = grid.add_point(p);
            grid.add_cell(vec![idx], VtkCellTypeW::Vertex);
        }
        grid.add_point_scalars("density", densities.to_vec());
        grid.add_point_scalars("pressure", pressures.to_vec());
        LegacyVtkWriter::write(&grid)
    }
}
/// A named data array for VTK point or cell data.
#[derive(Debug, Clone)]
pub enum VtkDataArrayW {
    /// Scalar (1-component) field.
    Scalars {
        /// Field name.
        name: String,
        /// One value per point or cell.
        values: Vec<f64>,
    },
    /// 3-component vector field.
    Vectors {
        /// Field name.
        name: String,
        /// One `[x, y, z]` tuple per point or cell.
        values: Vec<[f64; 3]>,
    },
    /// Symmetric 3x3 tensor field.
    Tensors {
        /// Field name.
        name: String,
        /// One 3x3 tensor per point or cell.
        values: Vec<[[f64; 3]; 3]>,
    },
}
impl VtkDataArrayW {
    /// Returns the name of this data array.
    pub fn name(&self) -> &str {
        match self {
            Self::Scalars { name, .. } => name,
            Self::Vectors { name, .. } => name,
            Self::Tensors { name, .. } => name,
        }
    }
    /// Returns the number of tuples in the array.
    pub fn len(&self) -> usize {
        match self {
            Self::Scalars { values, .. } => values.len(),
            Self::Vectors { values, .. } => values.len(),
            Self::Tensors { values, .. } => values.len(),
        }
    }
    /// Returns `true` if the array contains no tuples.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
/// Write a VTK XML UnstructuredGrid (.vtu) file (ASCII or raw binary inline).
///
/// Produces valid VTU XML that ParaView and VisIt can read.
pub struct VtuXmlWriter;
impl VtuXmlWriter {
    /// Write points and cells to a VTU XML file.
    ///
    /// `point_arrays`: list of (name, values) where values.len() == n_points.
    /// `cell_arrays`: list of (name, values) where values.len() == n_cells.
    #[allow(clippy::too_many_arguments)]
    pub fn write(
        path: &str,
        points: &[[f64; 3]],
        cells: &[Vec<usize>],
        cell_types: &[VtkCellTypeW],
        point_arrays: &[(&str, &[f64])],
        cell_arrays: &[(&str, &[f64])],
    ) -> std::io::Result<()> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        let n_pts = points.len();
        let n_cells = cells.len();
        let offsets: Vec<usize> = {
            let mut off = Vec::with_capacity(n_cells);
            let mut running = 0_usize;
            for c in cells {
                running += c.len();
                off.push(running);
            }
            off
        };
        writeln!(w, r#"<?xml version="1.0"?>"#)?;
        writeln!(
            w,
            r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#
        )?;
        writeln!(w, r#"  <UnstructuredGrid>"#)?;
        writeln!(
            w,
            r#"    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">"#
        )?;
        if !point_arrays.is_empty() {
            writeln!(w, r#"      <PointData>"#)?;
            for &(name, vals) in point_arrays {
                writeln!(
                    w,
                    r#"        <DataArray type="Float64" Name="{name}" format="ascii">"#
                )?;
                let row: Vec<String> = vals.iter().map(|v| format!("{v:.10e}")).collect();
                writeln!(w, "          {}", row.join(" "))?;
                writeln!(w, r#"        </DataArray>"#)?;
            }
            writeln!(w, r#"      </PointData>"#)?;
        }
        if !cell_arrays.is_empty() {
            writeln!(w, r#"      <CellData>"#)?;
            for &(name, vals) in cell_arrays {
                writeln!(
                    w,
                    r#"        <DataArray type="Float64" Name="{name}" format="ascii">"#
                )?;
                let row: Vec<String> = vals.iter().map(|v| format!("{v:.10e}")).collect();
                writeln!(w, "          {}", row.join(" "))?;
                writeln!(w, r#"        </DataArray>"#)?;
            }
            writeln!(w, r#"      </CellData>"#)?;
        }
        writeln!(w, r#"      <Points>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
        )?;
        for p in points {
            writeln!(w, "          {:.10e} {:.10e} {:.10e}", p[0], p[1], p[2])?;
        }
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(w, r#"      </Points>"#)?;
        writeln!(w, r#"      <Cells>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="Int64" Name="connectivity" format="ascii">"#
        )?;
        for c in cells {
            let row: Vec<String> = c.iter().map(|i| i.to_string()).collect();
            writeln!(w, "          {}", row.join(" "))?;
        }
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="Int64" Name="offsets" format="ascii">"#
        )?;
        let off_row: Vec<String> = offsets.iter().map(|o| o.to_string()).collect();
        writeln!(w, "          {}", off_row.join(" "))?;
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="UInt8" Name="types" format="ascii">"#
        )?;
        let types_row: Vec<String> = cell_types.iter().map(|t| (*t as u8).to_string()).collect();
        writeln!(w, "          {}", types_row.join(" "))?;
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(w, r#"      </Cells>"#)?;
        writeln!(w, r#"    </Piece>"#)?;
        writeln!(w, r#"  </UnstructuredGrid>"#)?;
        writeln!(w, r#"</VTKFile>"#)?;
        w.flush()?;
        Ok(())
    }
}
/// Utilities for writing numeric data in little-endian binary format.
pub struct VtkBinaryWriter;
impl VtkBinaryWriter {
    /// Encode a slice of `f64` values as little-endian bytes.
    pub fn encode_f64_le(values: &[f64]) -> Vec<u8> {
        let mut out = Vec::with_capacity(values.len() * 8);
        for &v in values {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }
    /// Encode a slice of `i64` values as little-endian bytes.
    pub fn encode_i64_le(values: &[i64]) -> Vec<u8> {
        let mut out = Vec::with_capacity(values.len() * 8);
        for &v in values {
            out.extend_from_slice(&v.to_le_bytes());
        }
        out
    }
    /// Encode a slice of `u8` values (pass-through).
    pub fn encode_u8(values: &[u8]) -> Vec<u8> {
        values.to_vec()
    }
    /// Base64-encode a byte slice (same algorithm as `XmlVtuWriter::encode_base64`).
    pub fn base64_encode(data: &[u8]) -> String {
        XmlVtuWriter::encode_base64(data)
    }
    /// Build a binary-format DataArray XML element for float64 data.
    ///
    /// The data is base64-encoded with a 4-byte little-endian length prefix.
    pub fn binary_data_array_xml(name: &str, n_components: usize, values: &[f64]) -> String {
        let raw = Self::encode_f64_le(values);
        let len_bytes = (raw.len() as u32).to_le_bytes();
        let mut payload = len_bytes.to_vec();
        payload.extend_from_slice(&raw);
        let encoded = Self::base64_encode(&payload);
        format!(
            "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"{}\" format=\"binary\">\n          {}\n        </DataArray>\n",
            name, n_components, encoded
        )
    }
}
/// Stateless collection of higher-level VTK writing helpers.
#[allow(dead_code)]
pub struct VtkWriter;
#[allow(dead_code)]
impl VtkWriter {
    /// Write an unstructured mesh with an attached vector field to a VTU file.
    ///
    /// # Arguments
    /// * `path`       – output file path (`.vtu`)
    /// * `points`     – node coordinates
    /// * `cells`      – connectivity (each element is a list of 0-based point indices)
    /// * `cell_types` – VTK cell type for each cell (parallel to `cells`)
    /// * `field_name` – name of the vector field (e.g. `"velocity"`)
    /// * `vectors`    – one `[vx, vy, vz]` per *node* (point data)
    pub fn write_unstructured_vector_field(
        path: &str,
        points: &[[f64; 3]],
        cells: &[Vec<usize>],
        cell_types: &[VtkCellTypeW],
        field_name: &str,
        vectors: &[[f64; 3]],
    ) -> std::io::Result<()> {
        let file = std::fs::File::create(path)?;
        let mut w = std::io::BufWriter::new(file);
        let n_pts = points.len();
        let n_cells = cells.len();
        let offsets: Vec<usize> = {
            let mut running = 0_usize;
            cells
                .iter()
                .map(|c| {
                    running += c.len();
                    running
                })
                .collect()
        };
        writeln!(w, r#"<?xml version="1.0"?>"#)?;
        writeln!(
            w,
            r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#
        )?;
        writeln!(w, r#"  <UnstructuredGrid>"#)?;
        writeln!(
            w,
            r#"    <Piece NumberOfPoints="{n_pts}" NumberOfCells="{n_cells}">"#
        )?;
        if !vectors.is_empty() {
            writeln!(w, r#"      <PointData>"#)?;
            writeln!(
                w,
                r#"        <DataArray type="Float64" Name="{field_name}" NumberOfComponents="3" format="ascii">"#
            )?;
            for v in vectors {
                writeln!(w, "          {:.10e} {:.10e} {:.10e}", v[0], v[1], v[2])?;
            }
            writeln!(w, r#"        </DataArray>"#)?;
            writeln!(w, r#"      </PointData>"#)?;
        }
        writeln!(w, r#"      <Points>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
        )?;
        for p in points {
            writeln!(w, "          {:.10e} {:.10e} {:.10e}", p[0], p[1], p[2])?;
        }
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(w, r#"      </Points>"#)?;
        writeln!(w, r#"      <Cells>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="Int64" Name="connectivity" format="ascii">"#
        )?;
        for c in cells {
            let row: Vec<String> = c.iter().map(|i| i.to_string()).collect();
            writeln!(w, "          {}", row.join(" "))?;
        }
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="Int64" Name="offsets" format="ascii">"#
        )?;
        let off_row: Vec<String> = offsets.iter().map(|o| o.to_string()).collect();
        writeln!(w, "          {}", off_row.join(" "))?;
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(
            w,
            r#"        <DataArray type="UInt8" Name="types" format="ascii">"#
        )?;
        let types_row: Vec<String> = cell_types.iter().map(|t| (*t as u8).to_string()).collect();
        writeln!(w, "          {}", types_row.join(" "))?;
        writeln!(w, r#"        </DataArray>"#)?;
        writeln!(w, r#"      </Cells>"#)?;
        writeln!(w, r#"    </Piece>"#)?;
        writeln!(w, r#"  </UnstructuredGrid>"#)?;
        writeln!(w, r#"</VTKFile>"#)?;
        Ok(())
    }
    /// Write a set of streamlines as VTK PolyData in legacy ASCII format.
    ///
    /// Each streamline is a polyline: a sequence of 3-D points.  The output
    /// is a valid legacy VTK file containing one `POLYDATA` dataset.
    ///
    /// # Arguments
    /// * `path`        – output file path (`.vtk`)
    /// * `streamlines` – each inner `Vec` is one polyline (≥ 2 points each)
    pub fn write_streamlines(path: &str, streamlines: &[Vec<[f64; 3]>]) -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        let n_points: usize = streamlines.iter().map(|s| s.len()).sum();
        let n_lines: usize = streamlines.len();
        let total_entries: usize = streamlines.iter().map(|s| s.len() + 1).sum();
        writeln!(f, "# vtk DataFile Version 3.0")?;
        writeln!(f, "Streamlines")?;
        writeln!(f, "ASCII")?;
        writeln!(f, "DATASET POLYDATA")?;
        writeln!(f, "POINTS {} double", n_points)?;
        for sl in streamlines {
            for p in sl {
                writeln!(f, "{} {} {}", p[0], p[1], p[2])?;
            }
        }
        writeln!(f, "LINES {} {}", n_lines, total_entries)?;
        let mut offset = 0usize;
        for sl in streamlines {
            let n = sl.len();
            let indices: Vec<String> = (offset..offset + n).map(|i| i.to_string()).collect();
            writeln!(f, "{} {}", n, indices.join(" "))?;
            offset += n;
        }
        Ok(())
    }
    /// Write a rectilinear grid to a legacy VTK ASCII file.
    ///
    /// A rectilinear grid is defined by three 1-D coordinate arrays
    /// `x_coords`, `y_coords`, `z_coords`.  The grid dimensions are
    /// `nx × ny × nz` where `nx = x_coords.len()` etc.
    ///
    /// Optionally attach a scalar point-data field `scalars` of length
    /// `nx * ny * nz` (pass an empty slice to omit).
    #[allow(clippy::too_many_arguments)]
    pub fn write_rectilinear_grid(
        path: &str,
        x_coords: &[f64],
        y_coords: &[f64],
        z_coords: &[f64],
        field_name: Option<&str>,
        scalars: &[f64],
    ) -> std::io::Result<()> {
        let mut f = std::fs::File::create(path)?;
        let nx = x_coords.len();
        let ny = y_coords.len();
        let nz = z_coords.len();
        let n_pts = nx * ny * nz;
        writeln!(f, "# vtk DataFile Version 3.0")?;
        writeln!(f, "Rectilinear Grid")?;
        writeln!(f, "ASCII")?;
        writeln!(f, "DATASET RECTILINEAR_GRID")?;
        writeln!(f, "DIMENSIONS {} {} {}", nx, ny, nz)?;
        writeln!(f, "X_COORDINATES {} double", nx)?;
        for &x in x_coords {
            write!(f, "{} ", x)?;
        }
        writeln!(f)?;
        writeln!(f, "Y_COORDINATES {} double", ny)?;
        for &y in y_coords {
            write!(f, "{} ", y)?;
        }
        writeln!(f)?;
        writeln!(f, "Z_COORDINATES {} double", nz)?;
        for &z in z_coords {
            write!(f, "{} ", z)?;
        }
        writeln!(f)?;
        if !scalars.is_empty() && scalars.len() == n_pts {
            let name = field_name.unwrap_or("scalar_field");
            writeln!(f, "POINT_DATA {}", n_pts)?;
            writeln!(f, "SCALARS {} double 1", name)?;
            writeln!(f, "LOOKUP_TABLE default")?;
            for &v in scalars {
                writeln!(f, "{}", v)?;
            }
        }
        Ok(())
    }
}
/// A piece descriptor for a parallel VTU file.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PvtuPiece {
    /// Relative path to the piece VTU file.
    pub filename: String,
}
impl PvtuPiece {
    /// Create a new piece entry.
    pub fn new(filename: impl Into<String>) -> Self {
        Self {
            filename: filename.into(),
        }
    }
}
/// Writer for parallel VTK datasets (PVTU format).
pub struct VtkParallelWriter;
impl VtkParallelWriter {
    /// Generate a PVTU XML file that references individual piece files.
    ///
    /// # Arguments
    /// * `partitions` – information about each piece.
    /// * `point_arrays` – list of `(name, n_components)` for point data arrays.
    /// * `cell_arrays` – list of `(name, n_components)` for cell data arrays.
    pub fn write_pvtu(
        partitions: &[VtkPartition],
        point_arrays: &[(&str, usize)],
        cell_arrays: &[(&str, usize)],
    ) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <PUnstructuredGrid GhostLevel=\"0\">\n");
        s.push_str("    <PPoints>\n");
        s.push_str("      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
        s.push_str("    </PPoints>\n");
        if !point_arrays.is_empty() {
            s.push_str("    <PPointData>\n");
            for (name, ncomp) in point_arrays {
                s.push_str(&format!(
                    "      <PDataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"{}\"/>\n",
                    name, ncomp
                ));
            }
            s.push_str("    </PPointData>\n");
        }
        if !cell_arrays.is_empty() {
            s.push_str("    <PCellData>\n");
            for (name, ncomp) in cell_arrays {
                s.push_str(&format!(
                    "      <PDataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"{}\"/>\n",
                    name, ncomp
                ));
            }
            s.push_str("    </PCellData>\n");
        }
        for part in partitions {
            s.push_str(&format!("    <Piece Source=\"{}\"/>\n", part.filename));
        }
        s.push_str("  </PUnstructuredGrid>\n</VTKFile>\n");
        s
    }
    /// Compute a balanced partition of `n_points` across `n_ranks` ranks.
    ///
    /// Returns a list of `(offset, count)` pairs.
    pub fn partition_points(n_points: usize, n_ranks: usize) -> Vec<(usize, usize)> {
        if n_ranks == 0 || n_points == 0 {
            return vec![];
        }
        let base = n_points / n_ranks;
        let extra = n_points % n_ranks;
        let mut partitions = Vec::with_capacity(n_ranks);
        let mut offset = 0;
        for r in 0..n_ranks {
            let count = base + if r < extra { 1 } else { 0 };
            partitions.push((offset, count));
            offset += count;
        }
        partitions
    }
}
