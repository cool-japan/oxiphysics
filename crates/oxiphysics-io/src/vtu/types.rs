//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::{Error, Result};
use oxiphysics_core::math::Vec3;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use super::functions::{VTK_HEX, VTK_QUAD, VTK_TET, VTK_TRIANGLE};

/// Serialize multiple `VtuWriter` pieces into a single VTU XML string
/// (each piece is a separate `Piece` inside `<UnstructuredGrid>`).
pub struct VtuMultiPiece;
impl VtuMultiPiece {
    /// Write multiple writers as pieces in a single VTU XML string.
    pub fn write(pieces: &[&VtuWriter]) -> String {
        let mut s = String::from("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <UnstructuredGrid>\n");
        for w in pieces {
            let xml = w.to_xml();
            if let Some(start) = xml.find("<Piece")
                && let Some(end) = xml.rfind("</Piece>")
            {
                s.push_str("    ");
                s.push_str(&xml[start..end + 8]);
                s.push('\n');
            }
        }
        s.push_str("  </UnstructuredGrid>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Total point count across all pieces.
    pub fn total_points(pieces: &[&VtuWriter]) -> usize {
        pieces.iter().map(|w| w.num_points()).sum()
    }
    /// Total cell count across all pieces.
    pub fn total_cells(pieces: &[&VtuWriter]) -> usize {
        pieces.iter().map(|w| w.num_cells()).sum()
    }
}
/// A simple metadata tag for annotating VTU files.
#[derive(Debug, Clone)]
pub struct VtuAnnotation {
    /// Key name.
    pub key: String,
    /// Value string.
    pub value: String,
}
pub(super) struct PointVector {
    pub(super) name: String,
    pub(super) vectors: Vec<[f64; 3]>,
}
/// Legacy static VTU writer (kept for backward compatibility with existing callers).
///
/// New code should use the builder-style [`VtuWriter`] struct instead.
pub struct VtuWriterLegacy;
impl VtuWriterLegacy {
    /// Write an unstructured grid to a VTU (VTK XML) file.
    ///
    /// Uses ASCII encoding. Cells are assumed to be tetrahedra (type 10).
    pub fn write(
        path: &str,
        positions: &[Vec3],
        cells: &[[usize; 4]],
        point_data: &[(&str, &[f64])],
        cell_data: &[(&str, &[f64])],
    ) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        let npoints = positions.len();
        let ncells = cells.len();
        writeln!(w, "<?xml version=\"1.0\"?>")?;
        writeln!(
            w,
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">"
        )?;
        writeln!(w, "  <UnstructuredGrid>")?;
        writeln!(
            w,
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">",
            npoints, ncells
        )?;
        if !point_data.is_empty() {
            writeln!(w, "      <PointData>")?;
            for (name, vals) in point_data {
                writeln!(
                    w,
                    "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
                    name
                )?;
                write!(w, "          ")?;
                for (i, v) in vals.iter().enumerate() {
                    if i > 0 {
                        write!(w, " ")?;
                    }
                    write!(w, "{}", v)?;
                }
                writeln!(w)?;
                writeln!(w, "        </DataArray>")?;
            }
            writeln!(w, "      </PointData>")?;
        }
        if !cell_data.is_empty() {
            writeln!(w, "      <CellData>")?;
            for (name, vals) in cell_data {
                writeln!(
                    w,
                    "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">",
                    name
                )?;
                write!(w, "          ")?;
                for (i, v) in vals.iter().enumerate() {
                    if i > 0 {
                        write!(w, " ")?;
                    }
                    write!(w, "{}", v)?;
                }
                writeln!(w)?;
                writeln!(w, "        </DataArray>")?;
            }
            writeln!(w, "      </CellData>")?;
        }
        writeln!(w, "      <Points>")?;
        writeln!(
            w,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"
        )?;
        for p in positions {
            writeln!(w, "          {} {} {}", p.x, p.y, p.z)?;
        }
        writeln!(w, "        </DataArray>")?;
        writeln!(w, "      </Points>")?;
        writeln!(w, "      <Cells>")?;
        writeln!(
            w,
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
        )?;
        for c in cells {
            writeln!(w, "          {} {} {} {}", c[0], c[1], c[2], c[3])?;
        }
        writeln!(w, "        </DataArray>")?;
        writeln!(
            w,
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
        )?;
        write!(w, "          ")?;
        for i in 0..ncells {
            if i > 0 {
                write!(w, " ")?;
            }
            write!(w, "{}", (i + 1) * 4)?;
        }
        writeln!(w)?;
        writeln!(w, "        </DataArray>")?;
        writeln!(
            w,
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">"
        )?;
        write!(w, "          ")?;
        for i in 0..ncells {
            if i > 0 {
                write!(w, " ")?;
            }
            write!(w, "10")?;
        }
        writeln!(w)?;
        writeln!(w, "        </DataArray>")?;
        writeln!(w, "      </Cells>")?;
        writeln!(w, "    </Piece>")?;
        writeln!(w, "  </UnstructuredGrid>")?;
        writeln!(w, "</VTKFile>")?;
        w.flush()?;
        Ok(())
    }
}
pub(super) struct CellVector {
    pub(super) name: String,
    pub(super) vectors: Vec<[f64; 3]>,
}
/// Builder-style writer for VTK XML unstructured grid (.vtu) files.
///
/// Add points, cells, and optional point data; then call [`VtuWriter::to_xml`] or
/// [`VtuWriter::write_to_file`].
pub struct VtuWriter {
    pub(super) points: Vec<[f64; 3]>,
    pub(super) cells: Vec<Vec<usize>>,
    pub(super) cell_types: Vec<u8>,
    pub(super) point_scalars: Vec<PointScalar>,
    pub(super) point_vectors: Vec<PointVector>,
    pub(super) cell_scalars: Vec<CellScalar>,
    pub(super) cell_vectors: Vec<CellVector>,
}
impl VtuWriter {
    /// Create an empty VTU writer.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            cells: Vec::new(),
            cell_types: Vec::new(),
            point_scalars: Vec::new(),
            point_vectors: Vec::new(),
            cell_scalars: Vec::new(),
            cell_vectors: Vec::new(),
        }
    }
    /// Add a point and return its index.
    pub fn add_point(&mut self, p: [f64; 3]) -> usize {
        let idx = self.points.len();
        self.points.push(p);
        idx
    }
    /// Add multiple points at once and return the index of the first added point.
    pub fn add_points(&mut self, pts: &[[f64; 3]]) -> usize {
        let first = self.points.len();
        self.points.extend_from_slice(pts);
        first
    }
    /// Return the number of points currently stored.
    pub fn num_points(&self) -> usize {
        self.points.len()
    }
    /// Return the number of cells currently stored.
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }
    /// Add a cell defined by its connectivity (0-based point indices) and VTK cell type.
    pub fn add_cell(&mut self, connectivity: Vec<usize>, cell_type_id: u8) {
        self.cells.push(connectivity);
        self.cell_types.push(cell_type_id);
    }
    /// Add a triangle cell from three point indices.
    pub fn add_triangle(&mut self, i0: usize, i1: usize, i2: usize) {
        self.add_cell(vec![i0, i1, i2], VTK_TRIANGLE);
    }
    /// Add a quad cell from four point indices.
    pub fn add_quad(&mut self, i0: usize, i1: usize, i2: usize, i3: usize) {
        self.add_cell(vec![i0, i1, i2, i3], VTK_QUAD);
    }
    /// Add a tetrahedron cell from four point indices.
    pub fn add_tet(&mut self, i0: usize, i1: usize, i2: usize, i3: usize) {
        self.add_cell(vec![i0, i1, i2, i3], VTK_TET);
    }
    /// Add a hexahedron cell from eight point indices.
    pub fn add_hex(
        &mut self,
        i0: usize,
        i1: usize,
        i2: usize,
        i3: usize,
        i4: usize,
        i5: usize,
        i6: usize,
        i7: usize,
    ) {
        self.add_cell(vec![i0, i1, i2, i3, i4, i5, i6, i7], VTK_HEX);
    }
    /// Attach a named scalar field on points.
    pub fn add_point_data_scalar(&mut self, name: &str, values: &[f64]) {
        self.point_scalars.push(PointScalar {
            name: name.to_string(),
            values: values.to_vec(),
        });
    }
    /// Attach a named vector field on points.
    pub fn add_point_data_vector(&mut self, name: &str, vectors: &[[f64; 3]]) {
        self.point_vectors.push(PointVector {
            name: name.to_string(),
            vectors: vectors.to_vec(),
        });
    }
    /// Write symmetric tensor cell data as a VTK XML snippet.
    ///
    /// Each cell supplies 6 components (Voigt notation): xx, yy, zz, xy, xz, yz.
    /// Returns an XML ``DataArray` element string ready to embed in ``CellData`.
    pub fn write_cell_data_tensor(name: &str, tensors: &[[f64; 6]]) -> String {
        let mut s = String::new();
        s.push_str(
            &format!(
                "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"6\" format=\"ascii\">\n          ",
                name
            ),
        );
        let vals: Vec<String> = tensors
            .iter()
            .flat_map(|t| t.iter())
            .map(|v| format!("{}", v))
            .collect();
        s.push_str(&vals.join(" "));
        s.push_str("\n        </DataArray>\n");
        s
    }
    /// Generate a parallel VTU master file (`.pvtu`) referencing the given piece files.
    ///
    /// `piece_files` is a slice of relative paths to the individual `.vtu` pieces.
    /// Returns the pvtu XML string.
    pub fn write_pvtu_parallel(
        piece_files: &[&str],
        point_data_names: &[&str],
        cell_data_names: &[&str],
    ) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <PUnstructuredGrid GhostLevel=\"0\">\n");
        s.push_str("    <PPoints>\n");
        s.push_str(
            "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>\n",
        );
        s.push_str("    </PPoints>\n");
        if !point_data_names.is_empty() {
            s.push_str("    <PPointData>\n");
            for name in point_data_names {
                s.push_str(&format!(
                    "      <PDataArray type=\"Float64\" Name=\"{}\" format=\"ascii\"/>\n",
                    name
                ));
            }
            s.push_str("    </PPointData>\n");
        }
        if !cell_data_names.is_empty() {
            s.push_str("    <PCellData>\n");
            for name in cell_data_names {
                s.push_str(&format!(
                    "      <PDataArray type=\"Float64\" Name=\"{}\" format=\"ascii\"/>\n",
                    name
                ));
            }
            s.push_str("    </PCellData>\n");
        }
        for file in piece_files {
            s.push_str(&format!("    <Piece Source=\"{}\"/>\n", file));
        }
        s.push_str("  </PUnstructuredGrid>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Attach a named scalar field on cells.
    pub fn add_cell_data_scalar(&mut self, name: &str, values: &[f64]) {
        self.cell_scalars.push(CellScalar {
            name: name.to_string(),
            values: values.to_vec(),
        });
    }
    /// Attach a named vector field on cells.
    pub fn add_cell_data_vector(&mut self, name: &str, vectors: &[[f64; 3]]) {
        self.cell_vectors.push(CellVector {
            name: name.to_string(),
            vectors: vectors.to_vec(),
        });
    }
    /// Compute the bounding box of all points.
    /// Returns `(min, max)` or `None` if there are no points.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.points.is_empty() {
            return None;
        }
        let mut min = self.points[0];
        let mut max = self.points[0];
        for p in &self.points[1..] {
            for d in 0..3 {
                if p[d] < min[d] {
                    min[d] = p[d];
                }
                if p[d] > max[d] {
                    max[d] = p[d];
                }
            }
        }
        Some((min, max))
    }
    /// Generate the complete VTU XML as a `String`.
    pub fn to_xml(&self) -> String {
        let npoints = self.points.len();
        let ncells = self.cells.len();
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <UnstructuredGrid>\n");
        s.push_str(&format!(
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            npoints, ncells
        ));
        let has_pdata = !self.point_scalars.is_empty() || !self.point_vectors.is_empty();
        if has_pdata {
            s.push_str("      <PointData>\n");
            for ps in &self.point_scalars {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">\n          ",
                        ps.name
                    ),
                );
                let vals: Vec<String> = ps.values.iter().map(|v| format!("{}", v)).collect();
                s.push_str(&vals.join(" "));
                s.push_str("\n        </DataArray>\n");
            }
            for pv in &self.point_vectors {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"3\" format=\"ascii\">\n          ",
                        pv.name
                    ),
                );
                let vals: Vec<String> = pv
                    .vectors
                    .iter()
                    .flat_map(|v| v.iter())
                    .map(|x| format!("{}", x))
                    .collect();
                s.push_str(&vals.join(" "));
                s.push_str("\n        </DataArray>\n");
            }
            s.push_str("      </PointData>\n");
        }
        let has_cdata = !self.cell_scalars.is_empty() || !self.cell_vectors.is_empty();
        if has_cdata {
            s.push_str("      <CellData>\n");
            for cs in &self.cell_scalars {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">\n          ",
                        cs.name
                    ),
                );
                let vals: Vec<String> = cs.values.iter().map(|v| format!("{}", v)).collect();
                s.push_str(&vals.join(" "));
                s.push_str("\n        </DataArray>\n");
            }
            for cv in &self.cell_vectors {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"3\" format=\"ascii\">\n          ",
                        cv.name
                    ),
                );
                let vals: Vec<String> = cv
                    .vectors
                    .iter()
                    .flat_map(|v| v.iter())
                    .map(|x| format!("{}", x))
                    .collect();
                s.push_str(&vals.join(" "));
                s.push_str("\n        </DataArray>\n");
            }
            s.push_str("      </CellData>\n");
        }
        s.push_str("      <Points>\n");
        s.push_str(
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        );
        for p in &self.points {
            s.push_str(&format!("          {} {} {}\n", p[0], p[1], p[2]));
        }
        s.push_str("        </DataArray>\n");
        s.push_str("      </Points>\n");
        s.push_str("      <Cells>\n");
        s.push_str(
            "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n          ",
        );
        let conn: Vec<String> = self
            .cells
            .iter()
            .flat_map(|c| c.iter())
            .map(|i| format!("{}", i))
            .collect();
        s.push_str(&conn.join(" "));
        s.push_str("\n        </DataArray>\n");
        s.push_str(
            "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n          ",
        );
        let mut offset = 0usize;
        let offsets: Vec<String> = self
            .cells
            .iter()
            .map(|c| {
                offset += c.len();
                format!("{}", offset)
            })
            .collect();
        s.push_str(&offsets.join(" "));
        s.push_str("\n        </DataArray>\n");
        s.push_str(
            "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ",
        );
        let type_strs: Vec<String> = self.cell_types.iter().map(|t| format!("{}", t)).collect();
        s.push_str(&type_strs.join(" "));
        s.push_str("\n        </DataArray>\n");
        s.push_str("      </Cells>\n");
        s.push_str("    </Piece>\n");
        s.push_str("  </UnstructuredGrid>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Write the VTU XML to a file.
    pub fn write_to_file(&self, path: &str) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        write!(w, "{}", self.to_xml())?;
        w.flush()?;
        Ok(())
    }
    /// Validate the writer state: check that cell connectivity references
    /// only existing point indices and that data arrays have correct length.
    pub fn validate(&self) -> Result<()> {
        let np = self.points.len();
        for (ci, cell) in self.cells.iter().enumerate() {
            for &idx in cell {
                if idx >= np {
                    return Err(Error::Parse(format!(
                        "cell {} references point index {} but only {} points exist",
                        ci, idx, np
                    )));
                }
            }
        }
        for ps in &self.point_scalars {
            if ps.values.len() != np {
                return Err(Error::Parse(format!(
                    "point scalar '{}' has {} values but {} points exist",
                    ps.name,
                    ps.values.len(),
                    np
                )));
            }
        }
        for pv in &self.point_vectors {
            if pv.vectors.len() != np {
                return Err(Error::Parse(format!(
                    "point vector '{}' has {} vectors but {} points exist",
                    pv.name,
                    pv.vectors.len(),
                    np
                )));
            }
        }
        let nc = self.cells.len();
        for cs in &self.cell_scalars {
            if cs.values.len() != nc {
                return Err(Error::Parse(format!(
                    "cell scalar '{}' has {} values but {} cells exist",
                    cs.name,
                    cs.values.len(),
                    nc
                )));
            }
        }
        for cv in &self.cell_vectors {
            if cv.vectors.len() != nc {
                return Err(Error::Parse(format!(
                    "cell vector '{}' has {} vectors but {} cells exist",
                    cv.name,
                    cv.vectors.len(),
                    nc
                )));
            }
        }
        Ok(())
    }
    /// Return the point positions slice.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }
    /// Return the cell types slice.
    pub fn cell_types(&self) -> &[u8] {
        &self.cell_types
    }
    /// Return the per-cell connectivity (each entry is a Vec of point indices).
    pub fn cells_per_cell(&self) -> &[Vec<usize>] {
        &self.cells
    }
    /// Return a flat connectivity array (concatenation of all cell index lists).
    pub fn cells_flat(&self) -> Vec<usize> {
        self.cells.iter().flat_map(|c| c.iter().cloned()).collect()
    }
}
/// Arithmetic operations on VTU scalar fields.
pub struct VtuFieldOps;
impl VtuFieldOps {
    /// Add two fields element-wise, returning the result.
    pub fn add(a: &[f64], b: &[f64]) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
    }
    /// Subtract field `b` from `a` element-wise.
    pub fn sub(a: &[f64], b: &[f64]) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
    }
    /// Multiply two fields element-wise.
    pub fn mul(a: &[f64], b: &[f64]) -> Vec<f64> {
        a.iter().zip(b.iter()).map(|(x, y)| x * y).collect()
    }
    /// Scale a field by a constant factor.
    pub fn scale(a: &[f64], s: f64) -> Vec<f64> {
        a.iter().map(|x| x * s).collect()
    }
    /// Compute the L2 norm (magnitude) of a vector field stored as interleaved [x, y, z, x, y, z, ...].
    pub fn vector_magnitude(v: &[f64]) -> Vec<f64> {
        v.chunks(3)
            .map(|c| {
                if c.len() == 3 {
                    (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
                } else {
                    0.0
                }
            })
            .collect()
    }
    /// Normalize a vector field so each vector has unit length.
    pub fn normalize_vectors(v: &[f64]) -> Vec<f64> {
        let mut out = Vec::with_capacity(v.len());
        for c in v.chunks(3) {
            if c.len() == 3 {
                let mag = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
                if mag > 1e-30 {
                    out.push(c[0] / mag);
                    out.push(c[1] / mag);
                    out.push(c[2] / mag);
                } else {
                    out.push(0.0);
                    out.push(0.0);
                    out.push(1.0);
                }
            }
        }
        out
    }
    /// Clamp field values to [lo, hi].
    pub fn clamp(a: &[f64], lo: f64, hi: f64) -> Vec<f64> {
        a.iter().map(|&x| x.max(lo).min(hi)).collect()
    }
    /// Compute the dot product of two vector fields (both stored as interleaved XYZ).
    pub fn dot_product(a: &[f64], b: &[f64]) -> Vec<f64> {
        a.chunks(3)
            .zip(b.chunks(3))
            .map(|(u, v)| {
                if u.len() == 3 && v.len() == 3 {
                    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
                } else {
                    0.0
                }
            })
            .collect()
    }
    /// Compute the cross product of two vector fields.
    pub fn cross_product(a: &[f64], b: &[f64]) -> Vec<f64> {
        let mut out = Vec::with_capacity(a.len());
        for (u, v) in a.chunks(3).zip(b.chunks(3)) {
            if u.len() == 3 && v.len() == 3 {
                out.push(u[1] * v[2] - u[2] * v[1]);
                out.push(u[2] * v[0] - u[0] * v[2]);
                out.push(u[0] * v[1] - u[1] * v[0]);
            }
        }
        out
    }
    /// Compute min, max, mean of a scalar field.
    pub fn stats(a: &[f64]) -> (f64, f64, f64) {
        if a.is_empty() {
            return (0.0, 0.0, 0.0);
        }
        let mut min_v = a[0];
        let mut max_v = a[0];
        let mut sum = 0.0_f64;
        for &v in a {
            if v < min_v {
                min_v = v;
            }
            if v > max_v {
                max_v = v;
            }
            sum += v;
        }
        (min_v, max_v, sum / a.len() as f64)
    }
    /// Threshold a field: return 1.0 where value > threshold, 0.0 otherwise.
    pub fn threshold(a: &[f64], threshold: f64) -> Vec<f64> {
        a.iter()
            .map(|&v| if v > threshold { 1.0 } else { 0.0 })
            .collect()
    }
}
pub(super) struct PointScalar {
    pub(super) name: String,
    pub(super) values: Vec<f64>,
}
pub(super) struct CellScalar {
    pub(super) name: String,
    pub(super) values: Vec<f64>,
}
/// Simple reader that extracts point coordinates from a VTU XML string.
pub struct VtuReader {
    pub(super) points: Vec<[f64; 3]>,
    pub(super) num_cells: usize,
    pub(super) cell_types: Vec<u8>,
    pub(super) point_scalar_names: Vec<String>,
    pub(super) cell_scalar_names: Vec<String>,
}
impl VtuReader {
    /// Parse a VTU XML string and extract point coordinates.
    pub fn from_xml(data: &str) -> Result<Self> {
        let mut points = Vec::new();
        let mut in_points_array = false;
        let mut cell_types = Vec::new();
        let mut num_cells = 0usize;
        let mut point_scalar_names = Vec::new();
        let mut cell_scalar_names = Vec::new();
        let mut in_cell_types = false;
        let mut in_connectivity = false;
        let mut in_point_data = false;
        let mut in_cell_data = false;
        for line in data.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with("<PointData") {
                in_point_data = true;
                continue;
            }
            if trimmed.starts_with("</PointData>") {
                in_point_data = false;
                continue;
            }
            if trimmed.starts_with("<CellData") {
                in_cell_data = true;
                continue;
            }
            if trimmed.starts_with("</CellData>") {
                in_cell_data = false;
                continue;
            }
            if (in_point_data || in_cell_data)
                && trimmed.contains("Name=\"")
                && let Some(start) = trimmed.find("Name=\"")
            {
                let rest = &trimmed[start + 6..];
                if let Some(end) = rest.find('"') {
                    let name = rest[..end].to_string();
                    if in_point_data {
                        point_scalar_names.push(name);
                    } else {
                        cell_scalar_names.push(name);
                    }
                }
            }
            if trimmed.contains("NumberOfCells=")
                && let Some(start) = trimmed.find("NumberOfCells=\"")
            {
                let rest = &trimmed[start + 15..];
                if let Some(end) = rest.find('"') {
                    num_cells = rest[..end].parse::<usize>().unwrap_or(0);
                }
            }
            if trimmed.contains("NumberOfComponents=\"3\"")
                && trimmed.contains("Float64")
                && !in_point_data
                && !in_cell_data
            {
                in_points_array = true;
                continue;
            }
            if in_points_array && trimmed.starts_with("</DataArray>") {
                in_points_array = false;
                continue;
            }
            if in_points_array && !trimmed.starts_with('<') && !trimmed.is_empty() {
                let nums: Vec<&str> = trimmed.split_whitespace().collect();
                if nums.len() >= 3 {
                    let mut i = 0;
                    while i + 2 < nums.len() {
                        let x = nums[i]
                            .parse::<f64>()
                            .map_err(|e| Error::Parse(e.to_string()))?;
                        let y = nums[i + 1]
                            .parse::<f64>()
                            .map_err(|e| Error::Parse(e.to_string()))?;
                        let z = nums[i + 2]
                            .parse::<f64>()
                            .map_err(|e| Error::Parse(e.to_string()))?;
                        points.push([x, y, z]);
                        i += 3;
                    }
                }
            }
            if trimmed.contains("Name=\"types\"") {
                in_cell_types = true;
                continue;
            }
            if in_cell_types && trimmed.starts_with("</DataArray>") {
                in_cell_types = false;
                continue;
            }
            if in_cell_types && !trimmed.starts_with('<') && !trimmed.is_empty() {
                for tok in trimmed.split_whitespace() {
                    if let Ok(t) = tok.parse::<u8>() {
                        cell_types.push(t);
                    }
                }
            }
            if trimmed.contains("Name=\"connectivity\"") {
                in_connectivity = true;
                continue;
            }
            if in_connectivity && trimmed.starts_with("</DataArray>") {
                in_connectivity = false;
                continue;
            }
            // connectivity data is tracked via in_connectivity state machine but not stored
        }
        Ok(Self {
            points,
            num_cells,
            cell_types,
            point_scalar_names,
            cell_scalar_names,
        })
    }
    /// Return parsed point coordinates.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }
    /// Return the number of cells reported.
    pub fn num_cells(&self) -> usize {
        self.num_cells
    }
    /// Return the cell types parsed from the file.
    pub fn cell_types(&self) -> &[u8] {
        &self.cell_types
    }
    /// Return point data array names found in the file.
    pub fn point_data_names(&self) -> &[String] {
        &self.point_scalar_names
    }
    /// Return cell data array names found in the file.
    pub fn cell_data_names(&self) -> &[String] {
        &self.cell_scalar_names
    }
    /// Read cell data arrays from VTU XML string by name.
    ///
    /// Returns a `Vec`f64` of the scalar values for the requested array,
    /// or an error if the array is not found in the ``CellData` section.
    pub fn read_cell_data(data: &str, array_name: &str) -> std::result::Result<Vec<f64>, String> {
        let mut in_cell_data = false;
        let mut target_found = false;
        let mut collecting = false;
        let mut values = Vec::new();
        for line in data.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with("<CellData") {
                in_cell_data = true;
                continue;
            }
            if trimmed.starts_with("</CellData>") {
                in_cell_data = false;
                target_found = false;
                collecting = false;
                continue;
            }
            if in_cell_data && trimmed.contains("Name=\"") {
                if let Some(start) = trimmed.find("Name=\"") {
                    let rest = &trimmed[start + 6..];
                    if let Some(end) = rest.find('"') {
                        let name = &rest[..end];
                        target_found = name == array_name;
                    }
                }
                if target_found {
                    collecting = true;
                }
                continue;
            }
            if collecting && trimmed.starts_with("</DataArray>") {
                collecting = false;
                continue;
            }
            if collecting && !trimmed.starts_with('<') && !trimmed.is_empty() {
                for tok in trimmed.split_whitespace() {
                    if let Ok(v) = tok.parse::<f64>() {
                        values.push(v);
                    }
                }
            }
        }
        if values.is_empty() && !data.contains(&format!("Name=\"{}\"", array_name)) {
            Err(format!("cell data array '{}' not found", array_name))
        } else {
            Ok(values)
        }
    }
    /// Return the bounding box of the parsed points.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.points.is_empty() {
            return None;
        }
        let mut min = self.points[0];
        let mut max = self.points[0];
        for p in &self.points[1..] {
            for d in 0..3 {
                if p[d] < min[d] {
                    min[d] = p[d];
                }
                if p[d] > max[d] {
                    max[d] = p[d];
                }
            }
        }
        Some((min, max))
    }
}
/// Writer for Parallel VTU (.pvtu) header files.
///
/// A PVTU file references a set of VTU piece files for parallel I/O.
pub struct PvtuWriter;
impl PvtuWriter {
    /// Generate a PVTU XML string referencing the given piece files.
    ///
    /// `point_data_names` lists the scalar point data arrays present in each piece.
    /// `cell_data_names` lists the scalar cell data arrays present in each piece.
    pub fn write(
        piece_files: &[&str],
        point_data_names: &[&str],
        cell_data_names: &[&str],
    ) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <PUnstructuredGrid>\n");
        if !point_data_names.is_empty() {
            s.push_str("    <PPointData>\n");
            for name in point_data_names {
                s.push_str(&format!(
                    "      <PDataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\"/>\n",
                    name
                ));
            }
            s.push_str("    </PPointData>\n");
        }
        if !cell_data_names.is_empty() {
            s.push_str("    <PCellData>\n");
            for name in cell_data_names {
                s.push_str(&format!(
                    "      <PDataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\"/>\n",
                    name
                ));
            }
            s.push_str("    </PCellData>\n");
        }
        s.push_str("    <PPoints>\n");
        s.push_str("      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
        s.push_str("    </PPoints>\n");
        for path in piece_files {
            s.push_str(&format!("    <Piece Source=\"{}\"/>\n", path));
        }
        s.push_str("  </PUnstructuredGrid>\n");
        s.push_str("</VTKFile>\n");
        s
    }
}
/// Convenience struct for writing a time series of VTU files and a PVD collection.
pub struct VtuTimeSeries {
    pub(super) prefix: String,
    pub(super) directory: String,
    pub(super) timesteps: Vec<(f64, String)>,
}
impl VtuTimeSeries {
    /// Create a new time series writer.
    ///
    /// `directory` is the output directory path.
    /// `prefix` is the file name prefix (e.g. "result" -> "result_0000.vtu").
    pub fn new(directory: &str, prefix: &str) -> Self {
        Self {
            prefix: prefix.to_string(),
            directory: directory.to_string(),
            timesteps: Vec::new(),
        }
    }
    /// Write a VTU file for the given time value and return the generated filename.
    pub fn write_frame(&mut self, time: f64, writer: &VtuWriter) -> Result<String> {
        let frame_idx = self.timesteps.len();
        let filename = format!("{}_{:04}.vtu", self.prefix, frame_idx);
        let filepath = format!("{}/{}", self.directory, filename);
        writer.write_to_file(&filepath)?;
        self.timesteps.push((time, filename));
        Ok(filepath)
    }
    /// Generate and write the PVD collection file.
    pub fn write_pvd(&self) -> Result<String> {
        let pvd_path = format!("{}/{}.pvd", self.directory, self.prefix);
        let entries: Vec<(f64, &str)> = self
            .timesteps
            .iter()
            .map(|(t, f)| (*t, f.as_str()))
            .collect();
        let pvd_xml = PvdWriter::write(&entries);
        let file = File::create(Path::new(&pvd_path))?;
        let mut w = BufWriter::new(file);
        write!(w, "{}", pvd_xml)?;
        w.flush()?;
        Ok(pvd_path)
    }
    /// Return the number of frames written so far.
    pub fn num_frames(&self) -> usize {
        self.timesteps.len()
    }
    /// Return the time values recorded so far.
    pub fn times(&self) -> Vec<f64> {
        self.timesteps.iter().map(|(t, _)| *t).collect()
    }
}
/// Quality metric for a single VTU cell.
#[derive(Debug, Clone)]
pub struct VtuCellQuality {
    /// Cell index (0-based).
    pub cell_index: usize,
    /// Cell type.
    pub cell_type: u8,
    /// Aspect ratio (min/max edge length, 1.0 = perfect).
    pub aspect_ratio: f64,
    /// Minimum edge length.
    pub min_edge: f64,
    /// Maximum edge length.
    pub max_edge: f64,
    /// Scaled Jacobian (1.0 = perfect).
    pub scaled_jacobian: f64,
}
/// Writer for ParaView Data (.pvd) collection files.
///
/// A PVD file references a list of VTU files at different timesteps,
/// enabling animation playback in ParaView.
pub struct PvdWriter;
impl PvdWriter {
    /// Generate a PVD XML string referencing the given timesteps.
    ///
    /// `timesteps` is a slice of `(time_value, vtu_file_path)`.
    pub fn write(timesteps: &[(f64, &str)]) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<VTKFile type=\"Collection\" version=\"0.1\">\n");
        s.push_str("  <Collection>\n");
        for (t, path) in timesteps {
            s.push_str(&format!(
                "    <DataSet timestep=\"{}\" part=\"0\" file=\"{}\"/>\n",
                t, path
            ));
        }
        s.push_str("  </Collection>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Generate a PVD XML string with part indices for parallel VTU files.
    ///
    /// Each timestep can have multiple parts (parallel pieces).
    pub fn write_parallel(timesteps: &[(f64, Vec<&str>)]) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<VTKFile type=\"Collection\" version=\"0.1\">\n");
        s.push_str("  <Collection>\n");
        for (t, paths) in timesteps {
            for (part, path) in paths.iter().enumerate() {
                s.push_str(&format!(
                    "    <DataSet timestep=\"{}\" part=\"{}\" file=\"{}\"/>\n",
                    t, part, path
                ));
            }
        }
        s.push_str("  </Collection>\n");
        s.push_str("</VTKFile>\n");
        s
    }
}
