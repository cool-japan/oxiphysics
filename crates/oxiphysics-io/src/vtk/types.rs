//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::Result;
use oxiphysics_core::math::Vec3;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// VTK XML Unstructured Grid.
///
/// Build up a grid with [`VtuGrid::add_point`] / [`VtuGrid::add_cell`], attach
/// field data via the `add_point_*` / `add_cell_*` helpers, then call
/// [`VtuGrid::to_vtu_string`] to obtain a ready-to-write VTU XML string.
pub struct VtuGrid {
    /// 3-D coordinates of every point.
    pub points: Vec<[f64; 3]>,
    /// Connectivity list: each entry is the ordered point indices of one cell.
    pub cells: Vec<Vec<usize>>,
    /// VTK cell type for each cell (parallel to [`VtuGrid::cells`]).
    pub cell_types: Vec<VtkCellType>,
    /// Per-point data arrays.
    pub point_data: Vec<VtkDataArray>,
    /// Per-cell data arrays.
    pub cell_data: Vec<VtkDataArray>,
}
impl VtuGrid {
    /// Create an empty grid.
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
    pub fn add_cell(&mut self, connectivity: Vec<usize>, cell_type: VtkCellType) {
        self.cells.push(connectivity);
        self.cell_types.push(cell_type);
    }
    /// Attach a per-point scalar field.
    pub fn add_point_scalar(&mut self, name: &str, values: Vec<f64>) {
        self.point_data.push(VtkDataArray::Scalar {
            name: name.to_owned(),
            values,
        });
    }
    /// Attach a per-point 3-component vector field.
    pub fn add_point_vector(&mut self, name: &str, values: Vec<[f64; 3]>) {
        self.point_data.push(VtkDataArray::Vector3 {
            name: name.to_owned(),
            values,
        });
    }
    /// Attach a per-cell scalar field.
    pub fn add_cell_scalar(&mut self, name: &str, values: Vec<f64>) {
        self.cell_data.push(VtkDataArray::Scalar {
            name: name.to_owned(),
            values,
        });
    }
    /// Number of points in the grid.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
    /// Number of cells in the grid.
    pub fn n_cells(&self) -> usize {
        self.cells.len()
    }
    /// Serialise the grid to a VTU XML string.
    ///
    /// The output is ASCII-encoded and compatible with ParaView / VisIt.
    pub fn to_vtu_string(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <UnstructuredGrid>\n");
        s.push_str(&format!(
            "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            self.n_points(),
            self.n_cells()
        ));
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
        s.push_str("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
        s.push_str("          ");
        let mut first = true;
        for conn in &self.cells {
            for &idx in conn {
                if !first {
                    s.push(' ');
                }
                s.push_str(&idx.to_string());
                first = false;
            }
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");
        s.push_str("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        s.push_str("          ");
        let mut offset: usize = 0;
        for (i, conn) in self.cells.iter().enumerate() {
            if i > 0 {
                s.push(' ');
            }
            offset += conn.len();
            s.push_str(&offset.to_string());
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");
        s.push_str("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
        s.push_str("          ");
        for (i, ct) in self.cell_types.iter().enumerate() {
            if i > 0 {
                s.push(' ');
            }
            s.push_str(&(*ct as u8).to_string());
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");
        s.push_str("      </Cells>\n");
        if !self.point_data.is_empty() {
            s.push_str("      <PointData>\n");
            for arr in &self.point_data {
                s.push_str(&Self::data_array_xml(arr));
            }
            s.push_str("      </PointData>\n");
        }
        if !self.cell_data.is_empty() {
            s.push_str("      <CellData>\n");
            for arr in &self.cell_data {
                s.push_str(&Self::data_array_xml(arr));
            }
            s.push_str("      </CellData>\n");
        }
        s.push_str("    </Piece>\n");
        s.push_str("  </UnstructuredGrid>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Render a [`VtkDataArray`] as a ``DataArray` XML element.
    fn data_array_xml(arr: &VtkDataArray) -> String {
        let mut s = String::new();
        match arr {
            VtkDataArray::Scalar { name, values } => {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">\n",
                        name
                    ),
                );
                s.push_str("          ");
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        s.push(' ');
                    }
                    s.push_str(&v.to_string());
                }
                s.push('\n');
                s.push_str("        </DataArray>\n");
            }
            VtkDataArray::Vector3 { name, values } => {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"3\" format=\"ascii\">\n",
                        name
                    ),
                );
                for v in values {
                    s.push_str(&format!("          {} {} {}\n", v[0], v[1], v[2]));
                }
                s.push_str("        </DataArray>\n");
            }
            VtkDataArray::Integer { name, values } => {
                s.push_str(
                    &format!(
                        "        <DataArray type=\"Int64\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">\n",
                        name
                    ),
                );
                s.push_str("          ");
                for (i, v) in values.iter().enumerate() {
                    if i > 0 {
                        s.push(' ');
                    }
                    s.push_str(&v.to_string());
                }
                s.push('\n');
                s.push_str("        </DataArray>\n");
            }
        }
        s
    }
    /// Write a PVD collection file that references a series of VTU snapshots.
    ///
    /// Returns the PVD XML as a `String`. The caller is responsible for writing
    /// it to disk and ensuring that the referenced VTU files are co-located.
    ///
    /// # Arguments
    /// * `base_name` – e.g. `"simulation"` (currently used as the collection title).
    /// * `time_steps` – slice of `(time, vtu_filename)` pairs.
    pub fn write_pvd_collection(base_name: &str, time_steps: &[(f64, String)]) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            &format!(
                "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n  <!-- {} -->\n",
                base_name
            ),
        );
        s.push_str("  <Collection>\n");
        for (time, filename) in time_steps {
            s.push_str(&format!(
                "    <DataSet timestep=\"{}\" group=\"\" part=\"0\" file=\"{}\"/>\n",
                time, filename
            ));
        }
        s.push_str("  </Collection>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Create a point-cloud grid: each position becomes a `Vertex` cell.
    pub fn from_points(positions: &[[f64; 3]]) -> Self {
        let mut grid = Self::new();
        for &p in positions {
            let idx = grid.add_point(p);
            grid.add_cell(vec![idx], VtkCellType::Vertex);
        }
        grid
    }
    /// Create a grid from a tetrahedral mesh.
    ///
    /// `nodes` are the 3-D coordinates; `elements` contains 4-node connectivity
    /// arrays, one per tetrahedron.
    pub fn from_tet_mesh(nodes: &[[f64; 3]], elements: &[[usize; 4]]) -> Self {
        let mut grid = Self::new();
        for &n in nodes {
            grid.add_point(n);
        }
        for &[a, b, c, d] in elements {
            grid.add_cell(vec![a, b, c, d], VtkCellType::Tetra);
        }
        grid
    }
}
/// A key-value pair in a VTK field data section.
#[derive(Debug, Clone)]
pub struct VtkFieldRecord {
    /// Field name.
    pub name: String,
    /// Scalar values.
    pub values: Vec<f64>,
}
impl VtkFieldRecord {
    /// Create a field record.
    pub fn new(name: impl Into<String>, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            values,
        }
    }
    /// Number of values.
    pub fn len(&self) -> usize {
        self.values.len()
    }
    /// Returns true if the record has no values.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}
/// VTK cell types for use with [`VtuGrid`].
#[derive(Debug, Clone, Copy)]
pub enum VtkCellType {
    /// VTK_VERTEX (1)
    Vertex = 1,
    /// VTK_LINE (3)
    Line = 3,
    /// VTK_TRIANGLE (5)
    Triangle = 5,
    /// VTK_QUAD (9)
    Quad = 9,
    /// VTK_TETRA (10)
    Tetra = 10,
    /// VTK_HEXAHEDRON (12)
    Hexahedron = 12,
    /// VTK_WEDGE (13)
    Wedge = 13,
    /// VTK_PYRAMID (14)
    Pyramid = 14,
}
/// A parsed VTK legacy ASCII dataset.
#[derive(Debug, Clone)]
pub struct VtkLegacyData {
    /// Dataset title line.
    pub title: String,
    /// Dataset type (e.g. "POLYDATA", "UNSTRUCTURED_GRID").
    pub dataset_type: String,
    /// Point coordinates (x, y, z).
    pub points: Vec<[f64; 3]>,
    /// Named scalar point data arrays.
    pub point_scalars: Vec<(String, Vec<f64>)>,
}
impl VtkLegacyData {
    /// Create an empty VTK legacy data container.
    pub fn empty() -> Self {
        Self {
            title: String::new(),
            dataset_type: String::new(),
            points: Vec::new(),
            point_scalars: Vec::new(),
        }
    }
}
/// A collection of VTK field data records.
///
/// Field data is global (not attached to points or cells) and is useful for
/// simulation metadata such as time, timestep index, and solver parameters.
#[derive(Debug, Clone, Default)]
pub struct VtkFieldData {
    /// The field records in this collection.
    pub records: Vec<VtkFieldRecord>,
}
impl VtkFieldData {
    /// Create an empty field data collection.
    pub fn new() -> Self {
        Self {
            records: Vec::new(),
        }
    }
    /// Add a named field.
    pub fn add(&mut self, name: impl Into<String>, values: Vec<f64>) {
        self.records.push(VtkFieldRecord::new(name, values));
    }
    /// Add a single scalar value as a 1-element field.
    pub fn add_scalar(&mut self, name: impl Into<String>, value: f64) {
        self.add(name, vec![value]);
    }
    /// Look up a field by name.
    pub fn get(&self, name: &str) -> Option<&VtkFieldRecord> {
        self.records.iter().find(|r| r.name == name)
    }
    /// Serialize as a VTK ASCII `FIELD` section.
    pub fn to_vtk_field_string(&self) -> String {
        if self.records.is_empty() {
            return String::new();
        }
        let mut s = String::new();
        s.push_str(&format!("FIELD FieldData {}\n", self.records.len()));
        for rec in &self.records {
            s.push_str(&format!("{} {} 1 float\n", rec.name, rec.values.len()));
            let vals: Vec<String> = rec.values.iter().map(|v| v.to_string()).collect();
            s.push_str(&vals.join(" "));
            s.push('\n');
        }
        s
    }
}
/// Writer for legacy VTK (.vtk) files.
pub struct VtkWriter;
impl VtkWriter {
    /// Write a point cloud to a legacy VTK file.
    pub fn write_points(path: &str, positions: &[Vec3]) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "OxiPhysics point cloud")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET POLYDATA")?;
        writeln!(w, "POINTS {} float", positions.len())?;
        for p in positions {
            writeln!(w, "{} {} {}", p.x, p.y, p.z)?;
        }
        w.flush()?;
        Ok(())
    }
    /// Write an unstructured grid with tetrahedral cells to a legacy VTK file.
    ///
    /// Optionally includes scalar and vector point data.
    pub fn write_unstructured_grid(
        path: &str,
        positions: &[Vec3],
        cells: &[[usize; 4]],
        scalars: Option<(&str, &[f64])>,
        vectors: Option<(&str, &[Vec3])>,
    ) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "OxiPhysics unstructured grid")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET UNSTRUCTURED_GRID")?;
        writeln!(w, "POINTS {} float", positions.len())?;
        for p in positions {
            writeln!(w, "{} {} {}", p.x, p.y, p.z)?;
        }
        let ncells = cells.len();
        let cell_size = ncells * 5;
        writeln!(w, "CELLS {} {}", ncells, cell_size)?;
        for c in cells {
            writeln!(w, "4 {} {} {} {}", c[0], c[1], c[2], c[3])?;
        }
        writeln!(w, "CELL_TYPES {}", ncells)?;
        for _ in 0..ncells {
            writeln!(w, "10")?;
        }
        let has_data = scalars.is_some() || vectors.is_some();
        if has_data {
            writeln!(w, "POINT_DATA {}", positions.len())?;
        }
        if let Some((name, vals)) = scalars {
            writeln!(w, "SCALARS {} float 1", name)?;
            writeln!(w, "LOOKUP_TABLE default")?;
            for v in vals {
                writeln!(w, "{}", v)?;
            }
        }
        if let Some((name, vecs)) = vectors {
            writeln!(w, "VECTORS {} float", name)?;
            for v in vecs {
                writeln!(w, "{} {} {}", v.x, v.y, v.z)?;
            }
        }
        w.flush()?;
        Ok(())
    }
    /// Write polygon (triangle) surface data to a legacy VTK file.
    pub fn write_polydata(path: &str, positions: &[Vec3], triangles: &[[usize; 3]]) -> Result<()> {
        let file = File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, "# vtk DataFile Version 3.0")?;
        writeln!(w, "OxiPhysics polydata")?;
        writeln!(w, "ASCII")?;
        writeln!(w, "DATASET POLYDATA")?;
        writeln!(w, "POINTS {} float", positions.len())?;
        for p in positions {
            writeln!(w, "{} {} {}", p.x, p.y, p.z)?;
        }
        let ntri = triangles.len();
        writeln!(w, "POLYGONS {} {}", ntri, ntri * 4)?;
        for t in triangles {
            writeln!(w, "3 {} {} {}", t[0], t[1], t[2])?;
        }
        w.flush()?;
        Ok(())
    }
}
/// A time series of `VtuGrid` snapshots.
pub struct VtkTimeSeries {
    /// The simulation time of each snapshot.
    pub times: Vec<f64>,
    /// The grid snapshot at each time.
    pub grids: Vec<VtuGrid>,
    /// Base file name (without extension) for output files.
    pub base_name: String,
}
impl VtkTimeSeries {
    /// Create an empty time series.
    pub fn new(base_name: &str) -> Self {
        Self {
            times: Vec::new(),
            grids: Vec::new(),
            base_name: base_name.to_owned(),
        }
    }
    /// Append a snapshot at time `t`.
    pub fn push(&mut self, t: f64, grid: VtuGrid) {
        self.times.push(t);
        self.grids.push(grid);
    }
    /// Number of snapshots.
    pub fn n_steps(&self) -> usize {
        self.times.len()
    }
    /// Generate a PVD collection file referencing all snapshots.
    ///
    /// File names are `{base_name}_{step:05}.vtu`.
    pub fn to_pvd_string(&self) -> String {
        let entries: Vec<(f64, String)> = self
            .times
            .iter()
            .enumerate()
            .map(|(i, &t)| (t, format!("{}_{:05}.vtu", self.base_name, i)))
            .collect();
        VtuGrid::write_pvd_collection(&self.base_name, &entries)
    }
    /// Get the VTU XML string for snapshot `i`.
    pub fn vtu_string(&self, i: usize) -> Option<String> {
        self.grids.get(i).map(|g| g.to_vtu_string())
    }
    /// Return estimated total data size in bytes (rough: 30 chars per point).
    pub fn estimated_size_bytes(&self) -> usize {
        self.grids.iter().map(|g| g.n_points() * 30).sum()
    }
}
/// A single time step entry in a PVD collection.
#[derive(Debug, Clone)]
pub struct PvdEntry {
    /// Simulation time.
    pub time: f64,
    /// Path to the VTU file for this time step.
    pub filename: String,
    /// Part name (optional, default "0").
    pub part: u32,
    /// Group name (optional).
    pub group: String,
}
impl PvdEntry {
    /// Create a new PVD entry.
    pub fn new(time: f64, filename: impl Into<String>) -> Self {
        Self {
            time,
            filename: filename.into(),
            part: 0,
            group: String::new(),
        }
    }
}
/// A ParaView-compatible multi-block dataset.
///
/// Wraps several [`VtuGrid`] blocks, each with a name.  Serializes to a
/// `<vtkMultiBlockDataSet>` XML document.
pub struct VtkMultiBlock {
    /// The named blocks comprising this dataset.
    pub blocks: Vec<VtkBlock>,
    /// Optional dataset title.
    pub title: String,
}
impl VtkMultiBlock {
    /// Create an empty multi-block dataset.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            blocks: Vec::new(),
            title: title.into(),
        }
    }
    /// Add a block.
    pub fn add_block(&mut self, name: impl Into<String>, grid: VtuGrid) {
        self.blocks.push(VtkBlock::new(name, grid));
    }
    /// Number of blocks.
    pub fn n_blocks(&self) -> usize {
        self.blocks.len()
    }
    /// Total number of points across all blocks.
    pub fn total_points(&self) -> usize {
        self.blocks.iter().map(|b| b.grid.n_points()).sum()
    }
    /// Total number of cells across all blocks.
    pub fn total_cells(&self) -> usize {
        self.blocks.iter().map(|b| b.grid.n_cells()).sum()
    }
    /// Serialize to a VTK multi-block VTM XML string.
    ///
    /// Each block references an external VTU file `{block_name}.vtu`.
    pub fn to_vtm_string(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"vtkMultiBlockDataSet\" version=\"1.0\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str("  <vtkMultiBlockDataSet>\n");
        for (i, block) in self.blocks.iter().enumerate() {
            s.push_str(&format!(
                "    <DataSet index=\"{}\" name=\"{}\" file=\"{}.vtu\"/>\n",
                i, block.name, block.name
            ));
        }
        s.push_str("  </vtkMultiBlockDataSet>\n");
        s.push_str("</VTKFile>\n");
        s
    }
    /// Write all VTU block files and a VTM index to the given directory.
    ///
    /// Returns the list of written file paths.
    pub fn write_to_dir(&self, dir: &str) -> crate::Result<Vec<String>> {
        use std::io::Write;
        let mut written = Vec::new();
        for block in &self.blocks {
            let path = format!("{}/{}.vtu", dir, block.name);
            let xml = block.grid.to_vtu_string();
            let mut f = std::fs::File::create(&path)?;
            f.write_all(xml.as_bytes())?;
            written.push(path);
        }
        let vtm_path = format!("{}/{}.vtm", dir, self.title);
        let vtm = self.to_vtm_string();
        let mut f = std::fs::File::create(&vtm_path)?;
        f.write_all(vtm.as_bytes())?;
        written.push(vtm_path);
        Ok(written)
    }
}
/// A VTK polydata dataset storing points, lines, and polygons.
pub struct VtkPolyDataGrid {
    /// 3-D point coordinates.
    pub points: Vec<[f64; 3]>,
    /// Line connectivity.
    pub lines: Vec<[usize; 2]>,
    /// Triangle connectivity.
    pub triangles: Vec<[usize; 3]>,
    /// Per-point data arrays.
    pub point_data: Vec<VtkDataArray>,
    /// Per-cell (poly) data arrays.
    pub cell_data: Vec<VtkDataArray>,
}
impl VtkPolyDataGrid {
    /// Create an empty polydata grid.
    pub fn new() -> Self {
        Self {
            points: Vec::new(),
            lines: Vec::new(),
            triangles: Vec::new(),
            point_data: Vec::new(),
            cell_data: Vec::new(),
        }
    }
    /// Add a point and return its index.
    pub fn add_point(&mut self, p: [f64; 3]) -> usize {
        let idx = self.points.len();
        self.points.push(p);
        idx
    }
    /// Add a triangle.
    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize) {
        self.triangles.push([a, b, c]);
    }
    /// Add a line.
    pub fn add_line(&mut self, a: usize, b: usize) {
        self.lines.push([a, b]);
    }
    /// Number of points.
    pub fn n_points(&self) -> usize {
        self.points.len()
    }
    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }
    /// Compute normals for each triangle (cross product, normalised).
    pub fn compute_triangle_normals(&self) -> Vec<[f64; 3]> {
        self.triangles
            .iter()
            .map(|&[a, b, c]| {
                let pa = self.points[a];
                let pb = self.points[b];
                let pc = self.points[c];
                let ab = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                let ac = [pc[0] - pa[0], pc[1] - pa[1], pc[2] - pa[2]];
                let n = [
                    ab[1] * ac[2] - ab[2] * ac[1],
                    ab[2] * ac[0] - ab[0] * ac[2],
                    ab[0] * ac[1] - ab[1] * ac[0],
                ];
                let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
                if len < 1e-12 {
                    [0.0, 0.0, 1.0]
                } else {
                    [n[0] / len, n[1] / len, n[2] / len]
                }
            })
            .collect()
    }
    /// Convert to a VTU unstructured grid (triangles become Triangle cells).
    pub fn to_vtu_grid(&self) -> VtuGrid {
        let mut g = VtuGrid::new();
        for &p in &self.points {
            g.add_point(p);
        }
        for &[a, b, c] in &self.triangles {
            g.add_cell(vec![a, b, c], VtkCellType::Triangle);
        }
        for &[a, b] in &self.lines {
            g.add_cell(vec![a, b], VtkCellType::Line);
        }
        g
    }
}
/// A scalar, vector, or integer data array attached to points or cells.
#[derive(Debug, Clone)]
pub enum VtkDataArray {
    /// Named scalar (1-component) field.
    Scalar {
        /// Field name.
        name: String,
        /// One value per point or cell.
        values: Vec<f64>,
    },
    /// Named 3-component vector field.
    Vector3 {
        /// Field name.
        name: String,
        /// One `\[x, y, z\]` tuple per point or cell.
        values: Vec<[f64; 3]>,
    },
    /// Named integer field.
    Integer {
        /// Field name.
        name: String,
        /// One value per point or cell.
        values: Vec<i64>,
    },
}
impl VtkDataArray {
    /// Returns the name of the data array.
    pub fn name(&self) -> &str {
        match self {
            Self::Scalar { name, .. } => name,
            Self::Vector3 { name, .. } => name,
            Self::Integer { name, .. } => name,
        }
    }
    /// Returns the number of components per tuple (1 for scalar/integer, 3 for vector).
    pub fn n_components(&self) -> usize {
        match self {
            Self::Scalar { .. } => 1,
            Self::Vector3 { .. } => 3,
            Self::Integer { .. } => 1,
        }
    }
    /// Returns the number of tuples in the array.
    pub fn len(&self) -> usize {
        match self {
            Self::Scalar { values, .. } => values.len(),
            Self::Vector3 { values, .. } => values.len(),
            Self::Integer { values, .. } => values.len(),
        }
    }
    /// Returns `true` if the array contains no tuples.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
/// A named block in a multi-block dataset.
pub struct VtkBlock {
    /// Block name.
    pub name: String,
    /// The grid for this block.
    pub grid: VtuGrid,
}
impl VtkBlock {
    /// Create a named block from a grid.
    pub fn new(name: impl Into<String>, grid: VtuGrid) -> Self {
        Self {
            name: name.into(),
            grid,
        }
    }
}
/// A VTK rectilinear grid (`.vtr`): axes are independently sampled.
pub struct VtkRectilinearGrid {
    /// Coordinate values along the X axis.
    pub x_coords: Vec<f64>,
    /// Coordinate values along the Y axis.
    pub y_coords: Vec<f64>,
    /// Coordinate values along the Z axis.
    pub z_coords: Vec<f64>,
    /// Per-point data arrays (ordered i,j,k like structured grid).
    pub point_data: Vec<VtkDataArray>,
}
impl VtkRectilinearGrid {
    /// Create from independent coordinate arrays.
    pub fn new(x_coords: Vec<f64>, y_coords: Vec<f64>, z_coords: Vec<f64>) -> Self {
        Self {
            x_coords,
            y_coords,
            z_coords,
            point_data: Vec::new(),
        }
    }
    /// Total number of grid points.
    pub fn n_points(&self) -> usize {
        self.x_coords.len() * self.y_coords.len() * self.z_coords.len()
    }
    /// Grid dimensions `\[ni, nj, nk\]`.
    pub fn dims(&self) -> [usize; 3] {
        [
            self.x_coords.len(),
            self.y_coords.len(),
            self.z_coords.len(),
        ]
    }
    /// Append a per-point scalar field.
    pub fn add_point_scalar(&mut self, name: &str, values: Vec<f64>) {
        self.point_data.push(VtkDataArray::Scalar {
            name: name.to_owned(),
            values,
        });
    }
    /// Serialise to VTK XML `.vtr` string.
    pub fn to_vtr_string(&self) -> String {
        let [ni, nj, nk] = self.dims();
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str(&format!(
            "  <RectilinearGrid WholeExtent=\"0 {} 0 {} 0 {}\">\n",
            ni.saturating_sub(1),
            nj.saturating_sub(1),
            nk.saturating_sub(1)
        ));
        s.push_str(&format!(
            "    <Piece Extent=\"0 {} 0 {} 0 {}\">\n",
            ni.saturating_sub(1),
            nj.saturating_sub(1),
            nk.saturating_sub(1)
        ));
        s.push_str("      <Coordinates>\n");
        for (label, coords) in [
            ("x", &self.x_coords),
            ("y", &self.y_coords),
            ("z", &self.z_coords),
        ] {
            s.push_str(&format!(
                "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">\n          ",
                label
            ));
            for (i, v) in coords.iter().enumerate() {
                if i > 0 {
                    s.push(' ');
                }
                s.push_str(&v.to_string());
            }
            s.push_str("\n        </DataArray>\n");
        }
        s.push_str("      </Coordinates>\n");
        if !self.point_data.is_empty() {
            s.push_str("      <PointData>\n");
            for arr in &self.point_data {
                if let VtkDataArray::Scalar { name, values } = arr {
                    s.push_str(
                        &format!(
                            "        <DataArray type=\"Float64\" Name=\"{}\" format=\"ascii\">\n          ",
                            name
                        ),
                    );
                    for (i, v) in values.iter().enumerate() {
                        if i > 0 {
                            s.push(' ');
                        }
                        s.push_str(&v.to_string());
                    }
                    s.push_str("\n        </DataArray>\n");
                }
            }
            s.push_str("      </PointData>\n");
        }
        s.push_str("    </Piece>\n  </RectilinearGrid>\n</VTKFile>\n");
        s
    }
}
/// A VTK structured grid (`.vts`) with explicit point coordinates.
///
/// Dimensions are `(ni, nj, nk)` in index space.  Points are ordered
/// with i varying fastest, then j, then k.
pub struct VtkStructuredGrid {
    /// Dimensions `\[ni, nj, nk\]`.
    pub dims: [usize; 3],
    /// 3-D coordinates of every point (ordered i,j,k).
    pub points: Vec<[f64; 3]>,
    /// Per-point data arrays.
    pub point_data: Vec<VtkDataArray>,
}
impl VtkStructuredGrid {
    /// Create an empty structured grid with the given dimensions.
    pub fn new(ni: usize, nj: usize, nk: usize) -> Self {
        Self {
            dims: [ni, nj, nk],
            points: Vec::with_capacity(ni * nj * nk),
            point_data: Vec::new(),
        }
    }
    /// Number of points (ni × nj × nk).
    pub fn n_points(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }
    /// Append a per-point scalar field.
    pub fn add_point_scalar(&mut self, name: &str, values: Vec<f64>) {
        self.point_data.push(VtkDataArray::Scalar {
            name: name.to_owned(),
            values,
        });
    }
    /// Serialise to VTK XML `.vts` format string.
    pub fn to_vts_string(&self) -> String {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str(
            "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
        );
        s.push_str(&format!(
            "  <StructuredGrid WholeExtent=\"0 {} 0 {} 0 {}\">\n",
            self.dims[0].saturating_sub(1),
            self.dims[1].saturating_sub(1),
            self.dims[2].saturating_sub(1)
        ));
        s.push_str(&format!(
            "    <Piece Extent=\"0 {} 0 {} 0 {}\">\n",
            self.dims[0].saturating_sub(1),
            self.dims[1].saturating_sub(1),
            self.dims[2].saturating_sub(1)
        ));
        s.push_str("      <Points>\n");
        s.push_str(
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        );
        for p in &self.points {
            s.push_str(&format!("          {} {} {}\n", p[0], p[1], p[2]));
        }
        s.push_str("        </DataArray>\n      </Points>\n");
        if !self.point_data.is_empty() {
            s.push_str("      <PointData>\n");
            for arr in &self.point_data {
                if let VtkDataArray::Scalar { name, values } = arr {
                    s.push_str(
                        &format!(
                            "        <DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">\n          ",
                            name
                        ),
                    );
                    for (i, v) in values.iter().enumerate() {
                        if i > 0 {
                            s.push(' ');
                        }
                        s.push_str(&v.to_string());
                    }
                    s.push_str("\n        </DataArray>\n");
                }
            }
            s.push_str("      </PointData>\n");
        }
        s.push_str("    </Piece>\n  </StructuredGrid>\n</VTKFile>\n");
        s
    }
    /// Build a uniform Cartesian grid covering `\[x0,x1\] x \[y0,y1\] x \[z0,z1\]`.
    pub fn uniform(
        x0: f64,
        x1: f64,
        ni: usize,
        y0: f64,
        y1: f64,
        nj: usize,
        z0: f64,
        z1: f64,
        nk: usize,
    ) -> Self {
        let mut g = Self::new(ni, nj, nk);
        let dx = if ni > 1 {
            (x1 - x0) / (ni - 1) as f64
        } else {
            0.0
        };
        let dy = if nj > 1 {
            (y1 - y0) / (nj - 1) as f64
        } else {
            0.0
        };
        let dz = if nk > 1 {
            (z1 - z0) / (nk - 1) as f64
        } else {
            0.0
        };
        for k in 0..nk {
            for j in 0..nj {
                for i in 0..ni {
                    g.points
                        .push([x0 + i as f64 * dx, y0 + j as f64 * dy, z0 + k as f64 * dz]);
                }
            }
        }
        g
    }
}
