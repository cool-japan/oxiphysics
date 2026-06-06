//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::fs::File;
use std::io::{self, Read, Write};

use super::functions::*;
use super::functions::{CHECKPOINT_MAGIC, CHECKPOINT_VERSION};

/// A collection of AMR cells at one refinement level.
#[derive(Debug, Clone)]
pub struct AmrGrid {
    /// Refinement level this grid represents.
    pub level_meta: AmrLevel,
    /// All cells in this grid.
    pub cells: Vec<AmrCell>,
    /// Names of the data fields stored in each cell.
    pub field_names: Vec<String>,
}
impl AmrGrid {
    /// Create a new empty grid.
    pub fn new(level_meta: AmrLevel, field_names: Vec<String>) -> Self {
        Self {
            level_meta,
            cells: Vec::new(),
            field_names,
        }
    }
    /// Add a cell to the grid.
    pub fn add_cell(&mut self, cell: AmrCell) {
        self.cells.push(cell);
    }
    /// Axis-aligned bounding box of all cells: `[i_min,j_min,k_min, i_max,j_max,k_max]`.
    pub fn integer_bounding_box(&self) -> Option<[i64; 6]> {
        if self.cells.is_empty() {
            return None;
        }
        let mut mn = self.cells[0].coords;
        let mut mx = self.cells[0].coords;
        for c in &self.cells {
            for d in 0..3 {
                if c.coords[d] < mn[d] {
                    mn[d] = c.coords[d];
                }
                if c.coords[d] > mx[d] {
                    mx[d] = c.coords[d];
                }
            }
        }
        Some([mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]])
    }
    /// Physical bounding box `[x_min,y_min,z_min, x_max,y_max,z_max]`.
    pub fn physical_bounding_box(&self) -> Option<[f64; 6]> {
        let ibb = self.integer_bounding_box()?;
        let cs = self.level_meta.cell_size;
        let o = &self.level_meta.domain_bounds;
        Some([
            o[0] + ibb[0] as f64 * cs,
            o[1] + ibb[1] as f64 * cs,
            o[2] + ibb[2] as f64 * cs,
            o[0] + (ibb[3] + 1) as f64 * cs,
            o[1] + (ibb[4] + 1) as f64 * cs,
            o[2] + (ibb[5] + 1) as f64 * cs,
        ])
    }
    /// Number of cells in this grid.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }
}
/// Coarsening operations: merge sibling blocks and project field data to the
/// coarser level.
pub struct AmrCoarsening;
impl AmrCoarsening {
    /// Coarsen fine-level cells into a new coarse grid.
    ///
    /// For each group of fine cells whose coordinates map to the same coarse
    /// cell (via integer division by `ratio`), the average field value is
    /// stored in the coarse cell.
    pub fn coarsen_grid(fine: &AmrGrid, coarse_level_meta: AmrLevel, ratio: i64) -> AmrGrid {
        let mut coarse = AmrGrid::new(coarse_level_meta, fine.field_names.clone());
        let num_fields = fine.field_names.len();
        let mut sums: std::collections::HashMap<[i64; 3], (Vec<f64>, usize)> =
            std::collections::HashMap::new();
        for cell in &fine.cells {
            let ck = [
                cell.coords[0] / ratio,
                cell.coords[1] / ratio,
                cell.coords[2] / ratio,
            ];
            let entry = sums
                .entry(ck)
                .or_insert_with(|| (vec![0.0_f64; num_fields], 0));
            for f in 0..num_fields {
                entry.0[f] += cell.data.get(f).copied().unwrap_or(0.0_f64);
            }
            entry.1 += 1;
        }
        let level_val = fine.level_meta.level.saturating_sub(1);
        for (ck, (sum, cnt)) in sums {
            let avg: Vec<f64> = sum.iter().map(|&s| s / cnt as f64).collect();
            coarse.add_cell(AmrCell::new(level_val, ck, avg));
        }
        coarse
    }
    /// Remove cells from a hierarchy level that are covered at finer levels.
    ///
    /// A coarse cell at `coarse_level` is "covered" when all children at
    /// `coarse_level + 1` exist in `fine_grid`.
    pub fn remove_covered_cells(coarse_grid: &mut AmrGrid, fine_grid: &AmrGrid, ratio: i64) {
        let fine_coords: std::collections::HashSet<[i64; 3]> =
            fine_grid.cells.iter().map(|c| c.coords).collect();
        coarse_grid.cells.retain(|cc| {
            let mut all_covered = true;
            for dz in 0..ratio {
                for dy in 0..ratio {
                    for dx in 0..ratio {
                        let child = [
                            cc.coords[0] * ratio + dx,
                            cc.coords[1] * ratio + dy,
                            cc.coords[2] * ratio + dz,
                        ];
                        if !fine_coords.contains(&child) {
                            all_covered = false;
                        }
                    }
                }
            }
            !all_covered
        });
    }
}
/// Export an AMR mesh to VTK UnstructuredGrid format (.vtu).
pub struct AmrMeshExport;
impl AmrMeshExport {
    /// Write an AMR hierarchy as a VTK UnstructuredGrid XML file.
    ///
    /// Each AMR cell becomes a VTK hexahedral cell (type 12).
    pub fn write_vtu(path: &str, hierarchy: &AmrHierarchy, field_idx: usize) -> io::Result<()> {
        let mut all_pts: Vec<[f64; 3]> = Vec::new();
        let mut connectivity: Vec<usize> = Vec::new();
        let mut offsets: Vec<usize> = Vec::new();
        let mut types: Vec<u8> = Vec::new();
        let mut field_vals: Vec<f64> = Vec::new();
        let mut offset = 0usize;
        for grid in &hierarchy.levels {
            let cs = grid.level_meta.cell_size;
            let o = &grid.level_meta.domain_bounds;
            for cell in &grid.cells {
                let x0 = o[0] + cell.coords[0] as f64 * cs;
                let y0 = o[1] + cell.coords[1] as f64 * cs;
                let z0 = o[2] + cell.coords[2] as f64 * cs;
                let x1 = x0 + cs;
                let y1 = y0 + cs;
                let z1 = z0 + cs;
                let base = all_pts.len();
                all_pts.push([x0, y0, z0]);
                all_pts.push([x1, y0, z0]);
                all_pts.push([x1, y1, z0]);
                all_pts.push([x0, y1, z0]);
                all_pts.push([x0, y0, z1]);
                all_pts.push([x1, y0, z1]);
                all_pts.push([x1, y1, z1]);
                all_pts.push([x0, y1, z1]);
                for p in base..base + 8 {
                    connectivity.push(p);
                }
                offset += 8;
                offsets.push(offset);
                types.push(12);
                field_vals.push(cell.data.get(field_idx).copied().unwrap_or(0.0_f64));
            }
        }
        let num_pts = all_pts.len();
        let num_cells = offsets.len();
        let mut f = File::create(path)?;
        writeln!(f, r#"<?xml version="1.0"?>"#)?;
        writeln!(
            f,
            r#"<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">"#
        )?;
        writeln!(f, r#"  <UnstructuredGrid>"#)?;
        writeln!(
            f,
            r#"    <Piece NumberOfPoints="{num_pts}" NumberOfCells="{num_cells}">"#
        )?;
        writeln!(f, r#"      <Points>"#)?;
        writeln!(
            f,
            r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
        )?;
        for pt in &all_pts {
            writeln!(f, "          {:.10} {:.10} {:.10}", pt[0], pt[1], pt[2])?;
        }
        writeln!(f, r#"        </DataArray>"#)?;
        writeln!(f, r#"      </Points>"#)?;
        writeln!(f, r#"      <Cells>"#)?;
        writeln!(
            f,
            r#"        <DataArray type="Int64" Name="connectivity" format="ascii">"#
        )?;
        let conn_str: Vec<String> = connectivity.iter().map(|&c| c.to_string()).collect();
        writeln!(f, "          {}", conn_str.join(" "))?;
        writeln!(f, r#"        </DataArray>"#)?;
        writeln!(
            f,
            r#"        <DataArray type="Int64" Name="offsets" format="ascii">"#
        )?;
        let off_str: Vec<String> = offsets.iter().map(|&o| o.to_string()).collect();
        writeln!(f, "          {}", off_str.join(" "))?;
        writeln!(f, r#"        </DataArray>"#)?;
        writeln!(
            f,
            r#"        <DataArray type="UInt8" Name="types" format="ascii">"#
        )?;
        let typ_str: Vec<String> = types.iter().map(|&t| t.to_string()).collect();
        writeln!(f, "          {}", typ_str.join(" "))?;
        writeln!(f, r#"        </DataArray>"#)?;
        writeln!(f, r#"      </Cells>"#)?;
        writeln!(f, r#"      <CellData>"#)?;
        writeln!(
            f,
            r#"        <DataArray type="Float64" Name="field_{field_idx}" format="ascii">"#
        )?;
        let fv_str: Vec<String> = field_vals.iter().map(|&v| format!("{v:.10}")).collect();
        writeln!(f, "          {}", fv_str.join(" "))?;
        writeln!(f, r#"        </DataArray>"#)?;
        writeln!(f, r#"      </CellData>"#)?;
        writeln!(f, r#"    </Piece>"#)?;
        writeln!(f, r#"  </UnstructuredGrid>"#)?;
        writeln!(f, r#"</VTKFile>"#)?;
        Ok(())
    }
}
/// A single AMR cell.
#[derive(Debug, Clone)]
pub struct AmrCell {
    /// Refinement level this cell belongs to.
    pub level: usize,
    /// Integer grid coordinates `(i, j, k)` within the level grid.
    pub coords: [i64; 3],
    /// Field data stored at this cell (one value per field).
    pub data: Vec<f64>,
}
impl AmrCell {
    /// Create a new cell.
    pub fn new(level: usize, coords: [i64; 3], data: Vec<f64>) -> Self {
        Self {
            level,
            coords,
            data,
        }
    }
    /// Compute the physical centre of this cell given level metadata.
    pub fn physical_center(&self, lvl: &AmrLevel) -> [f64; 3] {
        let cs = lvl.cell_size;
        let ox = lvl.domain_bounds[0];
        let oy = lvl.domain_bounds[1];
        let oz = lvl.domain_bounds[2];
        [
            ox + (self.coords[0] as f64 + 0.5) * cs,
            oy + (self.coords[1] as f64 + 0.5) * cs,
            oz + (self.coords[2] as f64 + 0.5) * cs,
        ]
    }
}
/// Metadata for one refinement level in an AMR hierarchy.
#[derive(Debug, Clone)]
pub struct AmrLevel {
    /// Zero-based level index (0 = coarsest).
    pub level: usize,
    /// Cell size at this refinement level (uniform isotropic).
    pub cell_size: f64,
    /// Domain bounds `[x_min, y_min, z_min, x_max, y_max, z_max]`.
    pub domain_bounds: [f64; 6],
}
impl AmrLevel {
    /// Create a new refinement level.
    pub fn new(level: usize, cell_size: f64, domain_bounds: [f64; 6]) -> Self {
        Self {
            level,
            cell_size,
            domain_bounds,
        }
    }
    /// Refinement ratio relative to the coarsest level (level 0): `2^level`.
    pub fn refinement_ratio(&self) -> usize {
        1 << self.level
    }
}
/// Binary checkpoint/restart for an AMR hierarchy.
///
/// Format (little-endian binary):
/// ```text
/// [magic: u32][version: u32][time: f64][cycle: u64]
/// [num_levels: u32]
/// For each level:
///   [level: u32][cell_size: f64][domain_bounds: 6×f64]
///   [num_field_names: u32]
///   For each field name:
///     [len: u32][bytes: len×u8]
///   [num_cells: u32]
///   For each cell:
///     [coords: 3×i64][num_data: u32][data: num_data×f64]
/// ```
pub struct AmrCheckpoint;
impl AmrCheckpoint {
    /// Write a checkpoint file.
    pub fn write(path: &str, hierarchy: &AmrHierarchy) -> io::Result<()> {
        let mut f = File::create(path)?;
        f.write_all(&CHECKPOINT_MAGIC.to_le_bytes())?;
        f.write_all(&CHECKPOINT_VERSION.to_le_bytes())?;
        f.write_all(&hierarchy.time.to_le_bytes())?;
        f.write_all(&hierarchy.cycle.to_le_bytes())?;
        f.write_all(&(hierarchy.levels.len() as u32).to_le_bytes())?;
        for grid in &hierarchy.levels {
            let lm = &grid.level_meta;
            f.write_all(&(lm.level as u32).to_le_bytes())?;
            f.write_all(&lm.cell_size.to_le_bytes())?;
            for &b in &lm.domain_bounds {
                f.write_all(&b.to_le_bytes())?;
            }
            f.write_all(&(grid.field_names.len() as u32).to_le_bytes())?;
            for name in &grid.field_names {
                let bytes = name.as_bytes();
                f.write_all(&(bytes.len() as u32).to_le_bytes())?;
                f.write_all(bytes)?;
            }
            f.write_all(&(grid.cells.len() as u32).to_le_bytes())?;
            for cell in &grid.cells {
                for &c in &cell.coords {
                    f.write_all(&c.to_le_bytes())?;
                }
                f.write_all(&(cell.data.len() as u32).to_le_bytes())?;
                for &d in &cell.data {
                    f.write_all(&d.to_le_bytes())?;
                }
            }
        }
        Ok(())
    }
    /// Read a checkpoint file written by `AmrCheckpoint::write`.
    pub fn read(path: &str) -> io::Result<AmrHierarchy> {
        let mut raw = Vec::new();
        File::open(path)?.read_to_end(&mut raw)?;
        let mut cursor = 0usize;
        macro_rules! read_u32 {
            () => {{
                let v =
                    u32::from_le_bytes(raw[cursor..cursor + 4].try_into().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "short read u32")
                    })?);
                cursor += 4;
                v
            }};
        }
        macro_rules! read_u64 {
            () => {{
                let v =
                    u64::from_le_bytes(raw[cursor..cursor + 8].try_into().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "short read u64")
                    })?);
                cursor += 8;
                v
            }};
        }
        macro_rules! read_i64 {
            () => {{
                let v =
                    i64::from_le_bytes(raw[cursor..cursor + 8].try_into().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "short read i64")
                    })?);
                cursor += 8;
                v
            }};
        }
        macro_rules! read_f64 {
            () => {{
                let v =
                    f64::from_le_bytes(raw[cursor..cursor + 8].try_into().map_err(|_| {
                        io::Error::new(io::ErrorKind::InvalidData, "short read f64")
                    })?);
                cursor += 8;
                v
            }};
        }
        let magic = read_u32!();
        if magic != CHECKPOINT_MAGIC {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "bad magic"));
        }
        let _version = read_u32!();
        let time = read_f64!();
        let cycle = read_u64!();
        let num_levels = read_u32!() as usize;
        let mut hierarchy = AmrHierarchy::new();
        hierarchy.time = time;
        hierarchy.cycle = cycle;
        for _ in 0..num_levels {
            let level = read_u32!() as usize;
            let cell_size = read_f64!();
            let mut domain_bounds = [0.0f64; 6];
            for b in &mut domain_bounds {
                *b = read_f64!();
            }
            let lm = AmrLevel::new(level, cell_size, domain_bounds);
            let num_fields = read_u32!() as usize;
            let mut field_names = Vec::with_capacity(num_fields);
            for _ in 0..num_fields {
                let len = read_u32!() as usize;
                let bytes = raw[cursor..cursor + len].to_vec();
                cursor += len;
                field_names.push(
                    String::from_utf8(bytes)
                        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "bad utf8"))?,
                );
            }
            let mut grid = AmrGrid::new(lm, field_names);
            let num_cells = read_u32!() as usize;
            for _ in 0..num_cells {
                let i = read_i64!();
                let j = read_i64!();
                let k = read_i64!();
                let nd = read_u32!() as usize;
                let mut data = Vec::with_capacity(nd);
                for _ in 0..nd {
                    data.push(read_f64!());
                }
                grid.add_cell(AmrCell::new(level, [i, j, k], data));
            }
            hierarchy.add_level(grid);
        }
        Ok(hierarchy)
    }
}
/// A multi-level AMR hierarchy.
#[derive(Debug, Clone)]
pub struct AmrHierarchy {
    /// Grids ordered from coarsest (index 0) to finest.
    pub levels: Vec<AmrGrid>,
    /// Simulation time (if applicable).
    pub time: f64,
    /// Cycle number (iteration count).
    pub cycle: u64,
}
impl AmrHierarchy {
    /// Create a new empty hierarchy.
    pub fn new() -> Self {
        Self {
            levels: Vec::new(),
            time: 0.0,
            cycle: 0,
        }
    }
    /// Add a grid level (must be added in order from coarsest to finest).
    pub fn add_level(&mut self, grid: AmrGrid) {
        self.levels.push(grid);
    }
    /// Number of refinement levels.
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }
    /// Total number of cells across all levels.
    pub fn total_cells(&self) -> usize {
        self.levels.iter().map(|g| g.cell_count()).sum()
    }
    /// Coarsen a fine-level field value to the parent level.
    ///
    /// Given a fine-level cell index and field index, returns the
    /// average of all child cells that map to the same coarse cell.
    pub fn coarsen_field(&self, fine_level: usize, field_idx: usize) -> Vec<f64> {
        if fine_level == 0 || fine_level >= self.levels.len() {
            return Vec::new();
        }
        let ratio = 2i64;
        let coarse = &self.levels[fine_level - 1];
        let fine = &self.levels[fine_level];
        let mut sums: std::collections::HashMap<[i64; 3], (f64, usize)> =
            std::collections::HashMap::new();
        for fc in &fine.cells {
            let ck = [
                fc.coords[0] / ratio,
                fc.coords[1] / ratio,
                fc.coords[2] / ratio,
            ];
            let val = fc.data.get(field_idx).copied().unwrap_or(0.0);
            let entry = sums.entry(ck).or_insert((0.0, 0));
            entry.0 += val;
            entry.1 += 1;
        }
        coarse
            .cells
            .iter()
            .map(|cc| {
                if let Some((sum, cnt)) = sums.get(&cc.coords) {
                    sum / *cnt as f64
                } else {
                    cc.data.get(field_idx).copied().unwrap_or(0.0)
                }
            })
            .collect()
    }
    /// Refine (interpolate) a coarse-level field to the next finer level.
    ///
    /// Uses nearest-neighbour injection from parent cell.
    pub fn refine_field(&self, coarse_level: usize, field_idx: usize) -> Vec<f64> {
        if coarse_level + 1 >= self.levels.len() {
            return Vec::new();
        }
        let ratio = 2i64;
        let coarse = &self.levels[coarse_level];
        let fine = &self.levels[coarse_level + 1];
        let lookup: std::collections::HashMap<[i64; 3], f64> = coarse
            .cells
            .iter()
            .map(|c| (c.coords, c.data.get(field_idx).copied().unwrap_or(0.0)))
            .collect();
        fine.cells
            .iter()
            .map(|fc| {
                let ck = [
                    fc.coords[0] / ratio,
                    fc.coords[1] / ratio,
                    fc.coords[2] / ratio,
                ];
                lookup.get(&ck).copied().unwrap_or(0.0)
            })
            .collect()
    }
}
/// Writer for AMR data in VTK XML format.
pub struct AmrWriter;
impl AmrWriter {
    /// Write an AMR hierarchy as a VTK HierarchicalBoxDataSet XML file.
    pub fn write_vtk(path: &str, hierarchy: &AmrHierarchy) -> io::Result<()> {
        let mut f = File::create(path)?;
        writeln!(f, r#"<?xml version="1.0"?>"#)?;
        writeln!(
            f,
            r#"<VTKFile type="vtkHierarchicalBoxDataSet" version="1.1" byte_order="LittleEndian">"#
        )?;
        writeln!(f, r#"  <vtkHierarchicalBoxDataSet>"#)?;
        for (li, grid) in hierarchy.levels.iter().enumerate() {
            writeln!(f, r#"    <Block level="{li}">"#)?;
            for cell in &grid.cells {
                let cs = grid.level_meta.cell_size;
                let o = &grid.level_meta.domain_bounds;
                let x0 = o[0] + cell.coords[0] as f64 * cs;
                let y0 = o[1] + cell.coords[1] as f64 * cs;
                let z0 = o[2] + cell.coords[2] as f64 * cs;
                let x1 = x0 + cs;
                let y1 = y0 + cs;
                let z1 = z0 + cs;
                writeln!(
                    f,
                    r#"      <DataSet index="0" amr_box="{} {} {} {} {} {}">"#,
                    cell.coords[0],
                    cell.coords[1],
                    cell.coords[2],
                    cell.coords[0],
                    cell.coords[1],
                    cell.coords[2]
                )?;
                writeln!(
                    f,
                    r#"        <!-- bounds: {x0:.6} {y0:.6} {z0:.6} {x1:.6} {y1:.6} {z1:.6} -->"#
                )?;
                writeln!(f, r#"        <PointData>"#)?;
                for (fi, fname) in grid.field_names.iter().enumerate() {
                    let val = cell.data.get(fi).copied().unwrap_or(0.0);
                    writeln!(
                        f,
                        r#"          <DataArray Name="{fname}" NumberOfTuples="1" format="ascii">{val:.10}</DataArray>"#
                    )?;
                }
                writeln!(f, r#"        </PointData>"#)?;
                writeln!(f, r#"      </DataSet>"#)?;
            }
            writeln!(f, r#"    </Block>"#)?;
        }
        writeln!(f, r#"  </vtkHierarchicalBoxDataSet>"#)?;
        writeln!(f, r#"</VTKFile>"#)?;
        Ok(())
    }
}
/// A block of field data on an AMR mesh patch.
///
/// Each block covers a rectangular region of cells at a single refinement level.
#[derive(Debug, Clone)]
pub struct AmrBlock {
    /// Refinement level of this block.
    pub level: usize,
    /// Lower-left integer corner `(i, j, k)`.
    pub lo: [i64; 3],
    /// Upper-right integer corner (inclusive) `(i, j, k)`.
    pub hi: [i64; 3],
    /// Number of cells in each direction: `[nx, ny, nz]`.
    pub dims: [usize; 3],
    /// Whether this block should be refined next.
    pub refinement_flag: bool,
    /// Field data stored in flat row-major order (z-major).
    ///
    /// Length = `nx * ny * nz * num_fields`.
    pub data: Vec<f64>,
    /// Number of fields stored per cell.
    pub num_fields: usize,
}
impl AmrBlock {
    /// Create a new zero-filled block.
    pub fn new(level: usize, lo: [i64; 3], hi: [i64; 3], num_fields: usize) -> Self {
        let dims = [
            (hi[0] - lo[0] + 1) as usize,
            (hi[1] - lo[1] + 1) as usize,
            (hi[2] - lo[2] + 1) as usize,
        ];
        let total = dims[0] * dims[1] * dims[2] * num_fields;
        Self {
            level,
            lo,
            hi,
            dims,
            refinement_flag: false,
            data: vec![0.0_f64; total],
            num_fields,
        }
    }
    /// Total number of cells.
    pub fn cell_count(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }
    /// Flat index of cell `(i, j, k)` for field `f`.
    pub fn index(&self, i: usize, j: usize, k: usize, f: usize) -> usize {
        ((k * self.dims[1] + j) * self.dims[0] + i) * self.num_fields + f
    }
    /// Read field `f` at cell `(i, j, k)`.
    pub fn get(&self, i: usize, j: usize, k: usize, f: usize) -> f64 {
        self.data[self.index(i, j, k, f)]
    }
    /// Write field `f` at cell `(i, j, k)`.
    pub fn set(&mut self, i: usize, j: usize, k: usize, f: usize, val: f64) {
        let idx = self.index(i, j, k, f);
        self.data[idx] = val;
    }
    /// Fill all cells of field `f` with `val`.
    pub fn fill_field(&mut self, f: usize, val: f64) {
        for k in 0..self.dims[2] {
            for j in 0..self.dims[1] {
                for i in 0..self.dims[0] {
                    self.set(i, j, k, f, val);
                }
            }
        }
    }
    /// Maximum value of field `f` across all cells.
    pub fn field_max(&self, f: usize) -> f64 {
        let mut m = f64::NEG_INFINITY;
        for k in 0..self.dims[2] {
            for j in 0..self.dims[1] {
                for i in 0..self.dims[0] {
                    let v = self.get(i, j, k, f);
                    if v > m {
                        m = v;
                    }
                }
            }
        }
        m
    }
    /// Minimum value of field `f` across all cells.
    pub fn field_min(&self, f: usize) -> f64 {
        let mut m = f64::INFINITY;
        for k in 0..self.dims[2] {
            for j in 0..self.dims[1] {
                for i in 0..self.dims[0] {
                    let v = self.get(i, j, k, f);
                    if v < m {
                        m = v;
                    }
                }
            }
        }
        m
    }
}
/// Method used to flag cells for refinement.
#[derive(Debug, Clone, PartialEq)]
pub enum RefinementMethod {
    /// Flag cells where the gradient of the specified field exceeds `threshold`.
    GradientBased {
        /// Index of the field to differentiate.
        field_idx: usize,
        /// Gradient magnitude threshold.
        threshold: f64,
    },
    /// Flag cells where the curvature (Laplacian) exceeds `threshold`.
    CurvatureBased {
        /// Index of the field.
        field_idx: usize,
        /// Laplacian threshold.
        threshold: f64,
    },
    /// Flag cells by a user-supplied predicate (cell data must exceed `value`).
    UserDefined {
        /// Index of the field.
        field_idx: usize,
        /// Threshold value.
        value: f64,
    },
}
/// Statistics about an AMR hierarchy.
#[derive(Debug, Clone)]
pub struct AmrStats {
    /// Number of cells per level.
    pub cells_per_level: Vec<usize>,
    /// Total number of cells.
    pub total_cells: usize,
    /// Estimated memory usage in bytes (f64 per field per cell).
    pub memory_bytes: usize,
}
impl AmrStats {
    /// Compute statistics for a hierarchy.
    pub fn compute(hierarchy: &AmrHierarchy) -> Self {
        let cells_per_level: Vec<usize> = hierarchy.levels.iter().map(|g| g.cell_count()).collect();
        let total_cells = cells_per_level.iter().sum();
        let fields = hierarchy
            .levels
            .first()
            .and_then(|g| g.field_names.len().into())
            .unwrap_or(0);
        let memory_bytes = total_cells * fields * std::mem::size_of::<f64>();
        Self {
            cells_per_level,
            total_cells,
            memory_bytes,
        }
    }
}
/// Reader for AMR data in the simple VTK XML format written by `AmrWriter`.
pub struct AmrReader;
impl AmrReader {
    /// Read an AMR hierarchy previously written by `AmrWriter::write_vtk`.
    ///
    /// This is a lightweight text parser for the specific format we write;
    /// it does not aim to support arbitrary VTK XML files.
    pub fn read_vtk(path: &str) -> io::Result<AmrHierarchy> {
        let mut raw = String::new();
        File::open(path)?.read_to_string(&mut raw)?;
        let mut hierarchy = AmrHierarchy::new();
        let mut current_level: Option<usize> = None;
        let mut current_coords: Option<[i64; 3]> = None;
        let mut current_data: Vec<f64> = Vec::new();
        let mut current_fields: Vec<String> = Vec::new();
        let mut grid_map: std::collections::BTreeMap<usize, AmrGrid> =
            std::collections::BTreeMap::new();
        for line in raw.lines() {
            let trimmed = line.trim();
            if trimmed.starts_with(r#"<Block level=""#) {
                if let Some(l) = parse_attr(trimmed, "level") {
                    current_level = l.parse::<usize>().ok();
                    if let Some(lv) = current_level {
                        grid_map.entry(lv).or_insert_with(|| {
                            AmrGrid::new(
                                AmrLevel::new(
                                    lv,
                                    1.0 / (1 << lv) as f64,
                                    [0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                                ),
                                Vec::new(),
                            )
                        });
                    }
                }
            } else if trimmed.starts_with("<DataSet") {
                if let Some(ab) = parse_attr(trimmed, "amr_box") {
                    let nums: Vec<i64> = ab
                        .split_whitespace()
                        .filter_map(|s| s.parse().ok())
                        .collect();
                    if nums.len() >= 3 {
                        current_coords = Some([nums[0], nums[1], nums[2]]);
                        current_data.clear();
                        current_fields.clear();
                    }
                }
            } else if trimmed.starts_with("<DataArray") {
                if let (Some(name), Some(val_str)) =
                    (parse_attr(trimmed, "Name"), parse_inline_value(trimmed))
                {
                    current_fields.push(name.to_string());
                    current_data.push(val_str.trim().parse::<f64>().unwrap_or(0.0));
                }
            } else if trimmed == "</DataSet>"
                && let (Some(lv), Some(coords)) = (current_level, current_coords.take())
                && let Some(grid) = grid_map.get_mut(&lv)
            {
                if grid.field_names.is_empty() {
                    grid.field_names = current_fields.clone();
                }
                grid.cells
                    .push(AmrCell::new(lv, coords, current_data.clone()));
            }
        }
        for (_, grid) in grid_map {
            hierarchy.add_level(grid);
        }
        Ok(hierarchy)
    }
}
/// A node in the `AmrtTree` oct-tree.
#[derive(Debug, Clone)]
pub struct AmrtNode {
    /// Refinement level (0 = coarsest).
    pub level: usize,
    /// Morton code at this node's position.
    pub morton: u64,
    /// Integer grid coordinates `(i, j, k)` at this level.
    pub coords: [u32; 3],
    /// Indices of up to 8 children (oct-tree).  `None` = no child.
    pub children: [Option<usize>; 8],
    /// Index of the parent node (`None` for the root).
    pub parent: Option<usize>,
    /// True if this node is a leaf.
    pub is_leaf: bool,
    /// 6-neighbour indices (±X, ±Y, ±Z). `None` if no same-level neighbour.
    pub neighbors: [Option<usize>; 6],
}
impl AmrtNode {
    /// Create a new leaf node.
    pub fn leaf(level: usize, coords: [u32; 3], parent: Option<usize>) -> Self {
        let [i, j, k] = coords;
        Self {
            level,
            morton: morton_encode(i, j, k),
            coords,
            children: [None; 8],
            parent,
            is_leaf: true,
            neighbors: [None; 6],
        }
    }
    /// Octant index for a child at relative position `(dx, dy, dz)` where each
    /// component is 0 or 1.
    pub fn octant_index(dx: usize, dy: usize, dz: usize) -> usize {
        dx | (dy << 1) | (dz << 2)
    }
}
/// Interpolation utilities for AMR fields.
pub struct AmrFieldInterp;
impl AmrFieldInterp {
    /// Linearly interpolate a field value at a physical point `(x, y, z)`.
    ///
    /// Searches for the finest level cell containing the point and returns
    /// its field value at `field_idx`.  Falls back to coarser levels if the
    /// point is not covered at the finest level.
    pub fn interpolate_at(
        hierarchy: &AmrHierarchy,
        point: [f64; 3],
        field_idx: usize,
    ) -> Option<f64> {
        for grid in hierarchy.levels.iter().rev() {
            let cs = grid.level_meta.cell_size;
            let o = &grid.level_meta.domain_bounds;
            let i = ((point[0] - o[0]) / cs).floor() as i64;
            let j = ((point[1] - o[1]) / cs).floor() as i64;
            let k = ((point[2] - o[2]) / cs).floor() as i64;
            for cell in &grid.cells {
                if cell.coords == [i, j, k] {
                    return cell.data.get(field_idx).copied();
                }
            }
        }
        None
    }
    /// Bilinear interpolation within a grid level at a physical point.
    ///
    /// Uses the four surrounding cells in the XY plane (Z ignored).
    /// Returns `None` if fewer than one corner cell is found.
    pub fn bilinear_xy(grid: &AmrGrid, x: f64, y: f64, field_idx: usize) -> Option<f64> {
        let cs = grid.level_meta.cell_size;
        let ox = grid.level_meta.domain_bounds[0];
        let oy = grid.level_meta.domain_bounds[1];
        let fi = (x - ox) / cs - 0.5;
        let fj = (y - oy) / cs - 0.5;
        let i0 = fi.floor() as i64;
        let j0 = fj.floor() as i64;
        let tx = fi - fi.floor();
        let ty = fj - fj.floor();
        let lookup: std::collections::HashMap<[i64; 3], f64> = grid
            .cells
            .iter()
            .map(|c| (c.coords, c.data.get(field_idx).copied().unwrap_or(0.0)))
            .collect();
        let k0 = 0i64;
        let v00 = lookup.get(&[i0, j0, k0]).copied()?;
        let v10 = lookup.get(&[i0 + 1, j0, k0]).copied().unwrap_or(v00);
        let v01 = lookup.get(&[i0, j0 + 1, k0]).copied().unwrap_or(v00);
        let v11 = lookup.get(&[i0 + 1, j0 + 1, k0]).copied().unwrap_or(v00);
        let v = v00 * (1.0 - tx) * (1.0 - ty)
            + v10 * tx * (1.0 - ty)
            + v01 * (1.0 - tx) * ty
            + v11 * tx * ty;
        Some(v)
    }
}
/// Applies refinement criteria to an `AmrHierarchy`.
///
/// Flags cells for refinement or coarsening based on the chosen method.
pub struct AmrRefinementCriteria {
    /// The method used to flag cells.
    pub method: RefinementMethod,
}
impl AmrRefinementCriteria {
    /// Create refinement criteria with the given method.
    pub fn new(method: RefinementMethod) -> Self {
        Self { method }
    }
    /// Flag cells in `grid` for refinement.
    ///
    /// Returns a `Vec`bool` parallel to `grid.cells` where `true` means the
    /// cell should be refined.
    pub fn flag_refinement(&self, grid: &AmrGrid) -> Vec<bool> {
        match &self.method {
            RefinementMethod::GradientBased {
                field_idx,
                threshold,
            } => self.flag_gradient(grid, *field_idx, *threshold),
            RefinementMethod::CurvatureBased {
                field_idx,
                threshold,
            } => self.flag_curvature(grid, *field_idx, *threshold),
            RefinementMethod::UserDefined { field_idx, value } => grid
                .cells
                .iter()
                .map(|c| c.data.get(*field_idx).copied().unwrap_or(0.0_f64) > *value)
                .collect(),
        }
    }
    fn flag_gradient(&self, grid: &AmrGrid, field_idx: usize, threshold: f64) -> Vec<bool> {
        let lookup: std::collections::HashMap<[i64; 3], f64> = grid
            .cells
            .iter()
            .map(|c| (c.coords, c.data.get(field_idx).copied().unwrap_or(0.0_f64)))
            .collect();
        grid.cells
            .iter()
            .map(|cell| {
                let v = cell.data.get(field_idx).copied().unwrap_or(0.0_f64);
                let mut grad_sq = 0.0_f64;
                for d in 0..3i64 {
                    let mut plus_c = cell.coords;
                    let mut minus_c = cell.coords;
                    plus_c[d as usize] += 1;
                    minus_c[d as usize] -= 1;
                    let vp = lookup.get(&plus_c).copied().unwrap_or(v);
                    let vm = lookup.get(&minus_c).copied().unwrap_or(v);
                    let g = (vp - vm) * 0.5_f64;
                    grad_sq += g * g;
                }
                grad_sq.sqrt() > threshold
            })
            .collect()
    }
    fn flag_curvature(&self, grid: &AmrGrid, field_idx: usize, threshold: f64) -> Vec<bool> {
        let lookup: std::collections::HashMap<[i64; 3], f64> = grid
            .cells
            .iter()
            .map(|c| (c.coords, c.data.get(field_idx).copied().unwrap_or(0.0_f64)))
            .collect();
        grid.cells
            .iter()
            .map(|cell| {
                let v = cell.data.get(field_idx).copied().unwrap_or(0.0_f64);
                let mut laplacian = 0.0_f64;
                for d in 0..3i64 {
                    let mut plus_c = cell.coords;
                    let mut minus_c = cell.coords;
                    plus_c[d as usize] += 1;
                    minus_c[d as usize] -= 1;
                    let vp = lookup.get(&plus_c).copied().unwrap_or(v);
                    let vm = lookup.get(&minus_c).copied().unwrap_or(v);
                    laplacian += vp + vm - 2.0_f64 * v;
                }
                laplacian.abs() > threshold
            })
            .collect()
    }
}
/// Serialise/deserialise an AMR hierarchy to a simple JSON-like text format.
///
/// Format:
/// ```text
/// {"time":`f64`,"cycle":`u64`,"levels":[
///   {"level":`usize`,"cell_size":`f64`,"cells":[
///     {"coords":[i,j,k],"data":[...]},...
///   ]},...
/// ]}
/// ```
pub struct AmrSerializer;
impl AmrSerializer {
    /// Serialise an `AmrHierarchy` to a JSON string.
    pub fn to_json(hierarchy: &AmrHierarchy) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "{{\"time\":{:.15e},\"cycle\":{},\"levels\":[",
            hierarchy.time, hierarchy.cycle
        ));
        for (li, grid) in hierarchy.levels.iter().enumerate() {
            if li > 0 {
                out.push(',');
            }
            out.push_str(&format!(
                "{{\"level\":{},\"cell_size\":{:.15e},\"domain\":[{},{},{},{},{},{}],\"fields\":[",
                grid.level_meta.level,
                grid.level_meta.cell_size,
                grid.level_meta.domain_bounds[0],
                grid.level_meta.domain_bounds[1],
                grid.level_meta.domain_bounds[2],
                grid.level_meta.domain_bounds[3],
                grid.level_meta.domain_bounds[4],
                grid.level_meta.domain_bounds[5],
            ));
            for (fi, name) in grid.field_names.iter().enumerate() {
                if fi > 0 {
                    out.push(',');
                }
                out.push('"');
                out.push_str(name);
                out.push('"');
            }
            out.push_str("],\"cells\":[");
            for (ci, cell) in grid.cells.iter().enumerate() {
                if ci > 0 {
                    out.push(',');
                }
                out.push_str(&format!(
                    "{{\"coords\":[{},{},{}],\"data\":[",
                    cell.coords[0], cell.coords[1], cell.coords[2]
                ));
                for (di, &d) in cell.data.iter().enumerate() {
                    if di > 0 {
                        out.push(',');
                    }
                    out.push_str(&format!("{:.15e}", d));
                }
                out.push_str("]}");
            }
            out.push_str("]}}");
        }
        out.push_str("]}");
        out
    }
    /// Deserialise an `AmrHierarchy` from JSON produced by `to_json`.
    ///
    /// This is a minimal hand-written parser for the specific format we emit.
    /// It does not handle arbitrary JSON.
    pub fn from_json(s: &str) -> Option<AmrHierarchy> {
        let time = extract_f64_after(s, "\"time\":")?;
        let cycle = extract_u64_after(s, "\"cycle\":")?;
        let mut hierarchy = AmrHierarchy::new();
        hierarchy.time = time;
        hierarchy.cycle = cycle;
        let mut search_from = 0;
        loop {
            let level_pos = s[search_from..].find("\"level\":")?;
            let abs_pos = search_from + level_pos;
            let level_val = extract_usize_at(&s[abs_pos + 8..])?;
            let cs_pos = s[abs_pos..].find("\"cell_size\":")?;
            let cell_size = extract_f64_after(&s[abs_pos + cs_pos..], "\"cell_size\":")?;
            let dom_pos = s[abs_pos..].find("\"domain\":")?;
            let dom_start = abs_pos + dom_pos + 9;
            let domain = extract_f64_array_6(&s[dom_start..])?;
            let fields_pos = s[abs_pos..].find("\"fields\":")?;
            let fields_start = abs_pos + fields_pos + 9;
            let field_names = extract_string_array(&s[fields_start..]);
            let lm = AmrLevel::new(level_val, cell_size, domain);
            let mut grid = AmrGrid::new(lm, field_names);
            let cells_pos = s[abs_pos..].find("\"cells\":")?;
            let cells_start = abs_pos + cells_pos + 8;
            let arr_start = cells_start + s[cells_start..].find('[')?;
            let arr_end = find_matching_bracket(&s[arr_start..])? + arr_start;
            let cells_str = &s[arr_start + 1..arr_end];
            let mut cell_from = 0;
            while let Some(obj_start) = cells_str[cell_from..].find('{') {
                let abs_obj = cell_from + obj_start;
                let obj_end_rel = find_matching_brace(&cells_str[abs_obj..])?;
                let obj_str = &cells_str[abs_obj..abs_obj + obj_end_rel + 1];
                let coords = extract_i64_array_3(obj_str)?;
                let data = extract_f64_array_vec(obj_str);
                grid.add_cell(AmrCell::new(level_val, coords, data));
                cell_from = abs_obj + obj_end_rel + 1;
            }
            hierarchy.add_level(grid);
            search_from = abs_pos + 8;
            if s[search_from..].find("\"level\":").is_none() {
                break;
            }
        }
        Some(hierarchy)
    }
}
/// Oct-tree AMR data structure.
///
/// Manages a spatially recursive oct-tree where each node can be refined into
/// 8 children.  Nodes are stored in a flat `Vec` for cache-friendly traversal.
#[derive(Debug, Default)]
pub struct AmrtTree {
    /// Flat storage of all nodes.
    pub nodes: Vec<AmrtNode>,
    /// Maximum allowed refinement level.
    pub max_level: usize,
    /// Domain bounds `\[x_min, y_min, z_min, x_max, y_max, z_max\]`.
    pub domain: [f64; 6],
}
impl AmrtTree {
    /// Create a new empty oct-tree.
    pub fn new(domain: [f64; 6], max_level: usize) -> Self {
        Self {
            nodes: Vec::new(),
            max_level,
            domain,
        }
    }
    /// Initialise the tree with a single root node at level 0.
    pub fn init_root(&mut self) {
        self.nodes.clear();
        self.nodes.push(AmrtNode::leaf(0, [0, 0, 0], None));
    }
    /// Refine the leaf at `node_idx` into 8 children.
    ///
    /// Returns `false` if the node is already at `max_level` or is not a leaf.
    pub fn refine(&mut self, node_idx: usize) -> bool {
        if node_idx >= self.nodes.len() {
            return false;
        }
        if !self.nodes[node_idx].is_leaf {
            return false;
        }
        let parent_level = self.nodes[node_idx].level;
        if parent_level >= self.max_level {
            return false;
        }
        let [pi, pj, pk] = self.nodes[node_idx].coords;
        let child_level = parent_level + 1;
        let mut child_indices = [usize::MAX; 8];
        for dz in 0..2usize {
            for dy in 0..2usize {
                for dx in 0..2usize {
                    let ci = pi * 2 + dx as u32;
                    let cj = pj * 2 + dy as u32;
                    let ck = pk * 2 + dz as u32;
                    let child_node = AmrtNode::leaf(child_level, [ci, cj, ck], Some(node_idx));
                    let child_idx = self.nodes.len();
                    child_indices[AmrtNode::octant_index(dx, dy, dz)] = child_idx;
                    self.nodes.push(child_node);
                }
            }
        }
        let mut children_opts = [None; 8];
        for (k, &ci) in child_indices.iter().enumerate() {
            children_opts[k] = Some(ci);
        }
        self.nodes[node_idx].children = children_opts;
        self.nodes[node_idx].is_leaf = false;
        true
    }
    /// Coarsen a set of 8 sibling leaves back to their parent.
    ///
    /// Returns `false` if any sibling is not a leaf.
    pub fn coarsen(&mut self, parent_idx: usize) -> bool {
        if parent_idx >= self.nodes.len() || self.nodes[parent_idx].is_leaf {
            return false;
        }
        let children: Vec<usize> = self.nodes[parent_idx]
            .children
            .iter()
            .filter_map(|&c| c)
            .collect();
        for &ci in &children {
            if !self.nodes[ci].is_leaf {
                return false;
            }
        }
        for &ci in &children {
            self.nodes[ci].is_leaf = false;
        }
        self.nodes[parent_idx].children = [None; 8];
        self.nodes[parent_idx].is_leaf = true;
        true
    }
    /// Count active (non-deleted) leaf nodes.
    pub fn leaf_count(&self) -> usize {
        self.nodes.iter().filter(|n| n.is_leaf).count()
    }
    /// Depth of the tree (maximum level among all active leaves).
    pub fn max_active_level(&self) -> usize {
        self.nodes
            .iter()
            .filter(|n| n.is_leaf)
            .map(|n| n.level)
            .max()
            .unwrap_or(0)
    }
    /// Physical cell size at a given refinement level.
    pub fn cell_size_at(&self, level: usize) -> [f64; 3] {
        let nx = (self.domain[3] - self.domain[0]) / (1u64 << level) as f64;
        let ny = (self.domain[4] - self.domain[1]) / (1u64 << level) as f64;
        let nz = (self.domain[5] - self.domain[2]) / (1u64 << level) as f64;
        [nx, ny, nz]
    }
    /// Physical centre of the node at `node_idx`.
    pub fn node_centre(&self, node_idx: usize) -> [f64; 3] {
        let n = &self.nodes[node_idx];
        let cs = self.cell_size_at(n.level);
        [
            self.domain[0] + (n.coords[0] as f64 + 0.5) * cs[0],
            self.domain[1] + (n.coords[1] as f64 + 0.5) * cs[1],
            self.domain[2] + (n.coords[2] as f64 + 0.5) * cs[2],
        ]
    }
}
