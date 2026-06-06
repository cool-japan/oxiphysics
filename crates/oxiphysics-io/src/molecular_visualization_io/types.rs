//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::functions::*;
use super::functions::{Lattice, Result, Vec3};

/// Coordinate mode used in a POSCAR file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PoscarCoordMode {
    /// Direct (fractional) coordinates.
    Direct,
    /// Cartesian coordinates.
    Cartesian,
}
/// Parsed contents of a PSF file.
#[derive(Debug, Clone)]
pub struct PsfFile {
    /// Title lines found after the `PSF` header.
    pub title: Vec<String>,
    /// All atom records.
    pub atoms: Vec<PsfAtom>,
    /// Covalent bonds.
    pub bonds: Vec<PsfBond>,
    /// Valence angles.
    pub angles: Vec<PsfAngle>,
    /// Whether this is an XPLOR-style PSF (extended atom types).
    pub xplor: bool,
}
impl PsfFile {
    /// Create an empty [`PsfFile`].
    pub fn empty() -> Self {
        Self {
            title: Vec::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
            xplor: false,
        }
    }
    /// Return the total charge as the sum of all atomic partial charges.
    pub fn total_charge(&self) -> f64 {
        self.atoms.iter().map(|a| a.charge).sum()
    }
    /// Return a map from `segname` to a list of atom indices (0-based).
    pub fn segments(&self) -> HashMap<String, Vec<usize>> {
        let mut map: HashMap<String, Vec<usize>> = HashMap::new();
        for (i, atom) in self.atoms.iter().enumerate() {
            map.entry(atom.segname.clone()).or_default().push(i);
        }
        map
    }
    /// Parse a PSF file from its full text content.
    ///
    /// Supports both standard and XPLOR (EXT) PSF dialects.
    pub fn parse(content: &str) -> Result<Self> {
        let mut psf = Self::empty();
        let mut lines = content.lines().peekable();
        let header =
            lines
                .find(|l| !l.trim().is_empty())
                .ok_or_else(|| MolVizError::ParseError {
                    message: "empty file".into(),
                })?;
        if !header.trim_start().starts_with("PSF") {
            return Err(MolVizError::ParseError {
                message: "not a PSF file".into(),
            });
        }
        psf.xplor = header.contains("EXT") || header.contains("XPLOR");
        let all_lines: Vec<&str> = lines.collect();
        let mut idx = 0usize;
        while idx < all_lines.len() {
            let line = all_lines[idx].trim();
            idx += 1;
            if line.is_empty() {
                continue;
            }
            if line.contains("!NTITLE") {
                let n: usize = line
                    .split_whitespace()
                    .next()
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
                for _ in 0..n {
                    if idx < all_lines.len() {
                        psf.title.push(all_lines[idx].to_string());
                        idx += 1;
                    }
                }
            } else if line.contains("!NATOM") {
                let n: usize = line
                    .split_whitespace()
                    .next()
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
                for _ in 0..n {
                    if idx >= all_lines.len() {
                        break;
                    }
                    let aline = all_lines[idx];
                    idx += 1;
                    let atom = parse_psf_atom_line(aline, psf.xplor)?;
                    psf.atoms.push(atom);
                }
            } else if line.contains("!NBOND") {
                let n_bonds: usize = line
                    .split_whitespace()
                    .next()
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
                let mut collected = 0usize;
                while collected < n_bonds && idx < all_lines.len() {
                    let bline = all_lines[idx];
                    idx += 1;
                    let tokens: Vec<&str> = bline.split_whitespace().collect();
                    let mut ti = 0;
                    while ti + 1 < tokens.len() && collected < n_bonds {
                        let i: u32 = tokens[ti].parse().unwrap_or(0);
                        let j: u32 = tokens[ti + 1].parse().unwrap_or(0);
                        if i > 0 && j > 0 {
                            psf.bonds.push(PsfBond { i, j });
                            collected += 1;
                        }
                        ti += 2;
                    }
                }
            } else if line.contains("!NTHETA") {
                let n_angles: usize = line
                    .split_whitespace()
                    .next()
                    .unwrap_or("0")
                    .parse()
                    .unwrap_or(0);
                let mut collected = 0usize;
                while collected < n_angles && idx < all_lines.len() {
                    let aline = all_lines[idx];
                    idx += 1;
                    let tokens: Vec<&str> = aline.split_whitespace().collect();
                    let mut ti = 0;
                    while ti + 2 < tokens.len() && collected < n_angles {
                        let i: u32 = tokens[ti].parse().unwrap_or(0);
                        let j: u32 = tokens[ti + 1].parse().unwrap_or(0);
                        let k: u32 = tokens[ti + 2].parse().unwrap_or(0);
                        if i > 0 && j > 0 && k > 0 {
                            psf.angles.push(PsfAngle { i, j, k });
                            collected += 1;
                        }
                        ti += 3;
                    }
                }
            }
        }
        Ok(psf)
    }
}
/// A single atom entry from a CIF `_atom_site_*` loop.
#[derive(Debug, Clone)]
pub struct CifAtomSite {
    /// Atom site label (e.g. `"Fe1"`).
    pub label: String,
    /// Element symbol.
    pub element: String,
    /// Fractional x coordinate.
    pub fx: f64,
    /// Fractional y coordinate.
    pub fy: f64,
    /// Fractional z coordinate.
    pub fz: f64,
    /// Site occupancy (0–1).
    pub occupancy: f64,
}
/// Simplified reader for VASP CHGCAR files.
///
/// The CHGCAR format is a POSCAR block followed by a charge density grid.
/// The charge density is stored as ρ × Volume (in electrons).
#[derive(Debug, Clone)]
pub struct ChgcarFile {
    /// Structure part of the CHGCAR.
    pub poscar: PoscarFile,
    /// Grid dimensions nx, ny, nz.
    pub grid_dims: [usize; 3],
    /// Charge density data (ρ × Ω), flattened, fastest index = x.
    pub charge: Vec<f64>,
}
impl ChgcarFile {
    /// Parse a CHGCAR file from text.
    pub fn parse(content: &str) -> Result<Self> {
        let lines: Vec<&str> = content.lines().collect();
        let poscar_lines = {
            let mut n_atoms = 0usize;
            let mut has_selective = false;
            for (i, line) in lines.iter().enumerate() {
                let t = line.trim();
                if i == 6 && t.chars().next().is_some_and(|c| c.is_alphabetic()) {
                    continue;
                }
                if i == 6 || i == 7 {
                    let counts: Vec<usize> = t
                        .split_whitespace()
                        .filter_map(|s| s.parse().ok())
                        .collect();
                    if !counts.is_empty() {
                        n_atoms = counts.iter().sum();
                    }
                }
                if (i == 7 || i == 8) && (t.starts_with('S') || t.starts_with('s')) {
                    has_selective = true;
                }
                if i > 8 {
                    break;
                }
            }
            let coord_offset = if has_selective { 9 } else { 8 };
            coord_offset + n_atoms
        };
        let mut search = poscar_lines;
        while search < lines.len() && !lines[search].trim().is_empty() {
            search += 1;
        }
        let poscar_end = search;
        let poscar_text: String = lines[..poscar_end].join("\n");
        let poscar = PoscarFile::parse(&poscar_text)?;
        let mut grid_idx = poscar_end + 1;
        while grid_idx < lines.len() && lines[grid_idx].trim().is_empty() {
            grid_idx += 1;
        }
        let grid_dims = if grid_idx < lines.len() {
            let gdline = lines[grid_idx].trim();
            let gd: Vec<usize> = gdline
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            [
                gd.first().copied().unwrap_or(1),
                gd.get(1).copied().unwrap_or(1),
                gd.get(2).copied().unwrap_or(1),
            ]
        } else {
            [1, 1, 1]
        };
        let expected = grid_dims[0] * grid_dims[1] * grid_dims[2];
        let mut charge = Vec::with_capacity(expected);
        let start = grid_idx + 1;
        for line in lines[start..].iter() {
            let t = line.trim();
            if t.is_empty() {
                break;
            }
            if t.starts_with("augmentation") || t.starts_with("AUGMENTATION") {
                break;
            }
            for tok in t.split_whitespace() {
                if let Ok(v) = tok.parse::<f64>() {
                    charge.push(v);
                }
            }
        }
        charge.resize(expected, 0.0);
        Ok(ChgcarFile {
            poscar,
            grid_dims,
            charge,
        })
    }
    /// Return the charge density at grid point (ix, iy, iz) divided by the cell volume.
    /// Gives electron density in e/Å³.
    pub fn density(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        let [_nx, ny, nz] = self.grid_dims;
        let idx = ix * ny * nz + iy * nz + iz;
        if idx < self.charge.len() {
            let vol = lattice_volume(self.poscar.lattice);
            if vol > 1e-12 {
                self.charge[idx] / vol
            } else {
                0.0
            }
        } else {
            0.0
        }
    }
}
/// Errors that may arise during molecular I/O.
#[derive(Debug)]
pub enum MolVizError {
    /// A line in the input could not be parsed.
    ParseError {
        /// Human-readable description of what failed.
        message: String,
    },
    /// An expected section was not found in the file.
    MissingSection {
        /// Name of the missing section keyword.
        section: String,
    },
    /// Grid dimensions are inconsistent.
    DimensionMismatch {
        /// Description.
        message: String,
    },
    /// Unknown or unsupported file format.
    UnknownFormat {
        /// Extension or magic string that was detected.
        detected: String,
    },
    /// Generic I/O error (wraps a message for `no_std` compatibility).
    IoError {
        /// Error message.
        message: String,
    },
    /// Index out of bounds in a trajectory or grid.
    OutOfBounds {
        /// Description.
        message: String,
    },
}
/// A trajectory assembled from multiple frames.
#[derive(Debug, Clone)]
pub struct Trajectory {
    /// Number of atoms (constant across frames).
    pub n_atoms: usize,
    /// All frames in order.
    pub frames: Vec<TrajectoryFrame>,
}
impl Trajectory {
    /// Construct an empty trajectory.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            n_atoms,
            frames: Vec::new(),
        }
    }
    /// Append a frame. Returns an error if the frame has the wrong atom count.
    pub fn push_frame(&mut self, frame: TrajectoryFrame) -> Result<()> {
        if frame.positions.len() != self.n_atoms {
            return Err(MolVizError::DimensionMismatch {
                message: format!(
                    "expected {} atoms, got {}",
                    self.n_atoms,
                    frame.positions.len()
                ),
            });
        }
        self.frames.push(frame);
        Ok(())
    }
    /// Return a slice of frames \[start, end).
    pub fn slice(&self, start: usize, end: usize) -> Result<Trajectory> {
        if start > end || end > self.frames.len() {
            return Err(MolVizError::OutOfBounds {
                message: format!(
                    "slice [{start},{end}) out of range for {} frames",
                    self.frames.len()
                ),
            });
        }
        let frames: Vec<TrajectoryFrame> = self.frames[start..end]
            .iter()
            .enumerate()
            .map(|(i, f)| TrajectoryFrame {
                index: i,
                time: f.time,
                positions: f.positions.clone(),
                velocities: f.velocities.clone(),
                cell: f.cell,
            })
            .collect();
        Ok(Trajectory {
            n_atoms: self.n_atoms,
            frames,
        })
    }
    /// Concatenate two trajectories (must have same `n_atoms`).
    pub fn concat(&self, other: &Trajectory) -> Result<Trajectory> {
        if self.n_atoms != other.n_atoms {
            return Err(MolVizError::DimensionMismatch {
                message: format!("n_atoms mismatch: {} vs {}", self.n_atoms, other.n_atoms),
            });
        }
        let offset = self.frames.len();
        let mut frames = self.frames.clone();
        for (i, f) in other.frames.iter().enumerate() {
            frames.push(TrajectoryFrame {
                index: offset + i,
                time: f.time,
                positions: f.positions.clone(),
                velocities: f.velocities.clone(),
                cell: f.cell,
            });
        }
        Ok(Trajectory {
            n_atoms: self.n_atoms,
            frames,
        })
    }
    /// Compute the mean-square displacement for atom `atom_idx` relative to frame 0.
    pub fn msd_atom(&self, atom_idx: usize) -> Result<Vec<f64>> {
        if atom_idx >= self.n_atoms {
            return Err(MolVizError::OutOfBounds {
                message: format!("atom index {atom_idx} >= n_atoms {}", self.n_atoms),
            });
        }
        if self.frames.is_empty() {
            return Ok(Vec::new());
        }
        let ref_pos = self.frames[0].positions[atom_idx];
        let msd = self
            .frames
            .iter()
            .map(|f| {
                let p = f.positions[atom_idx];
                let dx = p[0] - ref_pos[0];
                let dy = p[1] - ref_pos[1];
                let dz = p[2] - ref_pos[2];
                dx * dx + dy * dy + dz * dz
            })
            .collect();
        Ok(msd)
    }
    /// Sub-sample every `stride`-th frame.
    pub fn stride(&self, stride: usize) -> Trajectory {
        let stride = stride.max(1);
        let frames: Vec<TrajectoryFrame> = self
            .frames
            .iter()
            .step_by(stride)
            .enumerate()
            .map(|(i, f)| TrajectoryFrame {
                index: i,
                time: f.time,
                positions: f.positions.clone(),
                velocities: f.velocities.clone(),
                cell: f.cell,
            })
            .collect();
        Trajectory {
            n_atoms: self.n_atoms,
            frames,
        }
    }
    /// Number of frames.
    pub fn len(&self) -> usize {
        self.frames.len()
    }
    /// Return `true` if there are no frames.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}
/// An electrostatics volumetric grid in DX format.
#[derive(Debug, Clone)]
pub struct DxGrid {
    /// Number of grid points in x, y, z.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Number of grid points in z.
    pub nz: usize,
    /// Origin of the grid (angstroms).
    pub origin: Vec3,
    /// Voxel spacing along x (delta x).
    pub delta_x: Vec3,
    /// Voxel spacing along y.
    pub delta_y: Vec3,
    /// Voxel spacing along z.
    pub delta_z: Vec3,
    /// Flattened potential data (row-major: z varies fastest).
    pub data: Vec<f64>,
    /// Name of the quantity stored (e.g. `"electrostatic potential"`).
    pub quantity: String,
    /// Units string as found in the DX file comments.
    pub units: String,
}
impl DxGrid {
    /// Create an empty DX grid.
    pub fn empty(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            nx,
            ny,
            nz,
            origin: [0.0; 3],
            delta_x: [1.0, 0.0, 0.0],
            delta_y: [0.0, 1.0, 0.0],
            delta_z: [0.0, 0.0, 1.0],
            data: vec![0.0; nx * ny * nz],
            quantity: String::from("potential"),
            units: String::new(),
        }
    }
    /// Access a grid value by (ix, iy, iz) indices.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.data[ix * self.ny * self.nz + iy * self.nz + iz]
    }
    /// Set a grid value by (ix, iy, iz) indices.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, value: f64) {
        let idx = ix * self.ny * self.nz + iy * self.nz + iz;
        self.data[idx] = value;
    }
    /// Return the minimum and maximum potential in the grid.
    pub fn range(&self) -> (f64, f64) {
        let min = self.data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = self.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        (min, max)
    }
    /// Cartesian coordinate for grid point (ix, iy, iz).
    pub fn coord(&self, ix: usize, iy: usize, iz: usize) -> Vec3 {
        let i = ix as f64;
        let j = iy as f64;
        let k = iz as f64;
        [
            self.origin[0] + i * self.delta_x[0] + j * self.delta_y[0] + k * self.delta_z[0],
            self.origin[1] + i * self.delta_x[1] + j * self.delta_y[1] + k * self.delta_z[1],
            self.origin[2] + i * self.delta_x[2] + j * self.delta_y[2] + k * self.delta_z[2],
        ]
    }
    /// Parse a DX file from its full text content.
    pub fn parse(content: &str) -> Result<Self> {
        let mut nx = 0usize;
        let mut ny = 0usize;
        let mut nz = 0usize;
        let mut origin = [0.0f64; 3];
        let mut delta_x = [1.0f64, 0.0, 0.0];
        let mut delta_y = [0.0f64, 1.0, 0.0];
        let mut delta_z = [0.0f64, 0.0, 1.0];
        let mut data: Vec<f64> = Vec::new();
        let mut quantity = String::from("potential");
        let mut units = String::new();
        let mut delta_count = 0usize;
        let mut in_data = false;
        for raw_line in content.lines() {
            let line = raw_line.trim();
            if line.starts_with('#') {
                if line.contains("units") {
                    units = line.to_string();
                }
                continue;
            }
            if line.starts_with("object 1") && line.contains("gridpositions") {
                let toks: Vec<&str> = line.split_whitespace().collect();
                if toks.len() >= 8 {
                    nx = toks[5].parse().unwrap_or(0);
                    ny = toks[6].parse().unwrap_or(0);
                    nz = toks[7].parse().unwrap_or(0);
                }
                continue;
            }
            if line.starts_with("origin") {
                let toks: Vec<&str> = line.split_whitespace().collect();
                if toks.len() >= 4 {
                    origin[0] = toks[1].parse().unwrap_or(0.0);
                    origin[1] = toks[2].parse().unwrap_or(0.0);
                    origin[2] = toks[3].parse().unwrap_or(0.0);
                }
                continue;
            }
            if line.starts_with("delta") {
                let toks: Vec<&str> = line.split_whitespace().collect();
                if toks.len() >= 4 {
                    let dx = toks[1].parse().unwrap_or(0.0);
                    let dy = toks[2].parse().unwrap_or(0.0);
                    let dz = toks[3].parse().unwrap_or(0.0);
                    match delta_count {
                        0 => delta_x = [dx, dy, dz],
                        1 => delta_y = [dx, dy, dz],
                        _ => delta_z = [dx, dy, dz],
                    }
                    delta_count += 1;
                }
                continue;
            }
            if line.starts_with("object 2") {
                if line.contains("type double") || line.contains("type float") {
                    quantity = String::from("potential");
                }
                in_data = true;
                continue;
            }
            if line.starts_with("attribute") || line.starts_with("object 3") {
                in_data = false;
                continue;
            }
            if in_data {
                for token in line.split_whitespace() {
                    if let Ok(v) = token.parse::<f64>() {
                        data.push(v);
                    }
                }
            }
        }
        let expected = nx * ny * nz;
        if data.len() != expected && expected > 0 {
            data.resize(expected, 0.0);
        }
        Ok(DxGrid {
            nx,
            ny,
            nz,
            origin,
            delta_x,
            delta_y,
            delta_z,
            data,
            quantity,
            units,
        })
    }
    /// Serialize this grid back to DX text format.
    pub fn to_dx_string(&self) -> String {
        let mut out = String::new();
        out.push_str(&format!(
            "object 1 class gridpositions counts {} {} {}\n",
            self.nx, self.ny, self.nz
        ));
        out.push_str(&format!(
            "origin {:.6} {:.6} {:.6}\n",
            self.origin[0], self.origin[1], self.origin[2]
        ));
        out.push_str(&format!(
            "delta {:.6} {:.6} {:.6}\n",
            self.delta_x[0], self.delta_x[1], self.delta_x[2]
        ));
        out.push_str(&format!(
            "delta {:.6} {:.6} {:.6}\n",
            self.delta_y[0], self.delta_y[1], self.delta_y[2]
        ));
        out.push_str(&format!(
            "delta {:.6} {:.6} {:.6}\n",
            self.delta_z[0], self.delta_z[1], self.delta_z[2]
        ));
        let n = self.nx * self.ny * self.nz;
        out.push_str(&format!(
            "object 2 class array type double rank 0 items {n} data follows\n"
        ));
        for (i, v) in self.data.iter().enumerate() {
            out.push_str(&format!("{:.6}", v));
            if (i + 1) % 3 == 0 {
                out.push('\n');
            } else {
                out.push(' ');
            }
        }
        if !n.is_multiple_of(3) {
            out.push('\n');
        }
        out.push_str(
            "attribute \"dep\" string \"positions\"\nobject 3 class field\ncomponent \"positions\" value 1\ncomponent \"data\" value 2\n",
        );
        out
    }
}
/// Parsed XSF file with optional lattice and atoms.
#[derive(Debug, Clone)]
pub struct XsfFile {
    /// Periodicity class.
    pub periodicity: XsfPeriodicity,
    /// Primitive lattice vectors (rows) in angstroms. `None` for molecules.
    pub lattice: Option<Lattice>,
    /// Conventional cell vectors (may equal `lattice`). `None` for molecules.
    pub conventional: Option<Lattice>,
    /// List of atoms with element, position, and optional force.
    pub atoms: Vec<XsfAtom>,
}
impl XsfFile {
    /// Parse an XSF file from text content.
    pub fn parse(content: &str) -> Result<Self> {
        let mut periodicity = XsfPeriodicity::Molecule;
        let mut lattice: Option<Lattice> = None;
        let mut conventional: Option<Lattice> = None;
        let mut atoms: Vec<XsfAtom> = Vec::new();
        let lines: Vec<&str> = content.lines().collect();
        let mut idx = 0usize;
        while idx < lines.len() {
            let line = lines[idx].trim();
            idx += 1;
            match line {
                "CRYSTAL" => periodicity = XsfPeriodicity::Crystal,
                "SLAB" => periodicity = XsfPeriodicity::Slab,
                "POLYMER" => periodicity = XsfPeriodicity::Polymer,
                "MOLECULE" => periodicity = XsfPeriodicity::Molecule,
                "PRIMVEC" => {
                    let mat = read_3x3_matrix(&lines, &mut idx)?;
                    lattice = Some(mat);
                }
                "CONVVEC" => {
                    let mat = read_3x3_matrix(&lines, &mut idx)?;
                    conventional = Some(mat);
                }
                k if k.starts_with("PRIMCOORD") || k.starts_with("ATOMS") => {
                    if idx >= lines.len() {
                        break;
                    }
                    let count_line = lines[idx].trim();
                    idx += 1;
                    let nat: usize = count_line
                        .split_whitespace()
                        .next()
                        .and_then(|s| s.parse().ok())
                        .unwrap_or(0);
                    for _ in 0..nat {
                        if idx >= lines.len() {
                            break;
                        }
                        let aline = lines[idx].trim();
                        idx += 1;
                        if let Some(atom) = parse_xsf_atom(aline) {
                            atoms.push(atom);
                        }
                    }
                }
                _ => {}
            }
        }
        Ok(XsfFile {
            periodicity,
            lattice,
            conventional,
            atoms,
        })
    }
}
/// Parsed POSCAR or CONTCAR file.
#[derive(Debug, Clone)]
pub struct PoscarFile {
    /// Comment line at the top of the file.
    pub comment: String,
    /// Universal scaling factor.
    pub scale: f64,
    /// Lattice vectors a, b, c (rows, in angstroms after applying scale).
    pub lattice: Lattice,
    /// Species names (VASP 5+ format, or empty if absent).
    pub species: Vec<String>,
    /// Number of atoms per species.
    pub counts: Vec<usize>,
    /// Whether selective dynamics are enabled.
    pub selective_dynamics: bool,
    /// Coordinate mode.
    pub coord_mode: PoscarCoordMode,
    /// Atom positions (fractional or Cartesian depending on `coord_mode`).
    pub positions: Vec<Vec3>,
    /// Selective dynamics flags per atom (true = free to move).
    pub selective: Vec<[bool; 3]>,
}
impl PoscarFile {
    /// Parse a POSCAR/CONTCAR file from text.
    pub fn parse(content: &str) -> Result<Self> {
        let mut lines = content.lines();
        let comment = lines.next().unwrap_or("").to_string();
        let scale: f64 = lines.next().unwrap_or("1.0").trim().parse().unwrap_or(1.0);
        let mut lattice = [[0.0f64; 3]; 3];
        for row in lattice.iter_mut() {
            let l = lines.next().unwrap_or("");
            let t: Vec<f64> = l
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect();
            if t.len() >= 3 {
                row[0] = t[0] * scale;
                row[1] = t[1] * scale;
                row[2] = t[2] * scale;
            }
        }
        let line6 = lines.next().unwrap_or("").to_string();
        let has_species = line6
            .trim_start()
            .chars()
            .next()
            .is_some_and(|c| c.is_alphabetic());
        let (species, count_line_owned): (Vec<String>, String) = if has_species {
            let sp = line6.split_whitespace().map(|s| s.to_string()).collect();
            let cl = lines.next().unwrap_or("").to_string();
            (sp, cl)
        } else {
            (Vec::new(), line6)
        };
        let counts: Vec<usize> = count_line_owned
            .split_whitespace()
            .filter_map(|s| s.parse().ok())
            .collect();
        let n_total: usize = counts.iter().sum();
        let mut selective_dynamics = false;
        let line8 = lines.next().unwrap_or("").to_string();
        let coord_line_owned: String =
            if line8.trim().starts_with('S') || line8.trim().starts_with('s') {
                selective_dynamics = true;
                lines.next().unwrap_or("").to_string()
            } else {
                line8
            };
        let coord_mode = if coord_line_owned.trim().starts_with('D')
            || coord_line_owned.trim().starts_with('d')
        {
            PoscarCoordMode::Direct
        } else {
            PoscarCoordMode::Cartesian
        };
        let mut positions = Vec::with_capacity(n_total);
        let mut selective = Vec::with_capacity(n_total);
        for _ in 0..n_total {
            let aline = lines.next().unwrap_or("");
            let toks: Vec<&str> = aline.split_whitespace().collect();
            let x: f64 = toks.first().and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let y: f64 = toks.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let z: f64 = toks.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            positions.push([x, y, z]);
            if selective_dynamics && toks.len() >= 6 {
                let fx = toks[3] == "T";
                let fy = toks[4] == "T";
                let fz = toks[5] == "T";
                selective.push([fx, fy, fz]);
            } else {
                selective.push([true, true, true]);
            }
        }
        Ok(PoscarFile {
            comment,
            scale,
            lattice,
            species,
            counts,
            selective_dynamics,
            coord_mode,
            positions,
            selective,
        })
    }
    /// Convert all positions to Cartesian coordinates (angstroms).
    /// If already Cartesian, returns a clone of `positions`.
    pub fn cartesian_positions(&self) -> Vec<Vec3> {
        match self.coord_mode {
            PoscarCoordMode::Cartesian => self.positions.clone(),
            PoscarCoordMode::Direct => self
                .positions
                .iter()
                .map(|frac| {
                    let a = self.lattice[0];
                    let b = self.lattice[1];
                    let c = self.lattice[2];
                    [
                        frac[0] * a[0] + frac[1] * b[0] + frac[2] * c[0],
                        frac[0] * a[1] + frac[1] * b[1] + frac[2] * c[1],
                        frac[0] * a[2] + frac[1] * b[2] + frac[2] * c[2],
                    ]
                })
                .collect(),
        }
    }
    /// Total number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.counts.iter().sum()
    }
}
/// A section found in a Molden file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MoldenSection {
    /// `[Atoms]` section (geometry).
    Atoms,
    /// `[GTO]` section (Gaussian-type orbitals).
    Gto,
    /// `[MO]` section (molecular orbitals).
    Mo,
    /// `[FREQ]` section (vibrational frequencies).
    Freq,
    /// `[FR-COORD]` section (coordinates at frequency).
    FrCoord,
    /// Unknown or unparsed section.
    Other(String),
}
/// Parsed contents of a Molden file.
#[derive(Debug, Clone)]
pub struct MoldenFile {
    /// Atoms from the `[Atoms]` section.
    pub atoms: Vec<Atom>,
    /// Coordinate units: `"Angs"` or `"AU"`.
    pub coord_units: String,
    /// Molecular orbitals.
    pub mos: Vec<MoldenMo>,
    /// Vibrational frequencies (cm⁻¹).
    pub frequencies: Vec<f64>,
}
impl MoldenFile {
    /// Parse a Molden file from text.
    pub fn parse(content: &str) -> Result<Self> {
        let mut atoms: Vec<Atom> = Vec::new();
        let mut coord_units = String::from("Angs");
        let mut mos: Vec<MoldenMo> = Vec::new();
        let mut frequencies: Vec<f64> = Vec::new();
        let lines: Vec<&str> = content.lines().collect();
        let mut idx = 0usize;
        while idx < lines.len() {
            let line = lines[idx].trim();
            idx += 1;
            if line.starts_with("[Atoms]") {
                if line.contains("AU") || line.contains("au") {
                    coord_units = String::from("AU");
                } else {
                    coord_units = String::from("Angs");
                }
                while idx < lines.len() {
                    let aline = lines[idx].trim();
                    if aline.starts_with('[') {
                        break;
                    }
                    idx += 1;
                    if aline.is_empty() || aline.starts_with('#') {
                        continue;
                    }
                    let toks: Vec<&str> = aline.split_whitespace().collect();
                    if toks.len() < 5 {
                        continue;
                    }
                    let element = toks[0].to_string();
                    let x: f64 = toks[2].parse().unwrap_or(0.0);
                    let y: f64 = toks[3].parse().unwrap_or(0.0);
                    let z: f64 = toks[4].parse().unwrap_or(0.0);
                    atoms.push(Atom::new(element, [x, y, z]));
                }
            } else if line.starts_with("[MO]") {
                while idx < lines.len() {
                    let mline = lines[idx].trim();
                    if mline.starts_with('[') {
                        break;
                    }
                    idx += 1;
                    if mline.is_empty() {
                        continue;
                    }
                    if mline.starts_with("Sym=") || mline.starts_with("Ene=") {
                        let mut mo = MoldenMo {
                            symmetry: String::new(),
                            energy: 0.0,
                            spin: String::from("Alpha"),
                            occupation: 0.0,
                            coefficients: Vec::new(),
                        };
                        let kv_line = mline;
                        parse_molden_kv(kv_line, &mut mo);
                        while idx < lines.len() {
                            let bline = lines[idx].trim();
                            if bline.starts_with('[')
                                || bline.starts_with("Sym=")
                                || bline.starts_with("Ene=")
                            {
                                break;
                            }
                            idx += 1;
                            if bline.is_empty() {
                                continue;
                            }
                            parse_molden_kv(bline, &mut mo);
                            let ct: Vec<&str> = bline.split_whitespace().collect();
                            if ct.len() == 2
                                && let (Ok(i), Ok(v)) =
                                    (ct[0].parse::<usize>(), ct[1].parse::<f64>())
                            {
                                mo.coefficients.push((i, v));
                            }
                        }
                        mos.push(mo);
                    }
                }
            } else if line.starts_with("[FREQ]") {
                while idx < lines.len() {
                    let fline = lines[idx].trim();
                    if fline.starts_with('[') {
                        break;
                    }
                    idx += 1;
                    if let Ok(f) = fline.parse::<f64>() {
                        frequencies.push(f);
                    }
                }
            }
        }
        Ok(MoldenFile {
            atoms,
            coord_units,
            mos,
            frequencies,
        })
    }
}
/// File format identifiers for auto-detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MolFormat {
    /// PSF (Protein Structure File).
    Psf,
    /// DX electrostatics grid.
    Dx,
    /// Gaussian CUBE file.
    Cube,
    /// XCrysDen XSF file.
    Xsf,
    /// VASP POSCAR / CONTCAR.
    Poscar,
    /// Crystallographic Information File.
    Cif,
    /// Molden orbital / frequency file.
    Molden,
    /// VASP CHGCAR charge density.
    Chgcar,
    /// Unknown format.
    Unknown,
}
/// A single snapshot in a multi-frame trajectory.
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Frame index (0-based).
    pub index: usize,
    /// Simulation time (ps or fs, context-dependent).
    pub time: f64,
    /// Atom positions.
    pub positions: Vec<Vec3>,
    /// Optional velocities.
    pub velocities: Option<Vec<Vec3>>,
    /// Optional box vectors (3 × 3).
    pub cell: Option<Lattice>,
}
/// A simplified CIF data block.
#[derive(Debug, Clone)]
pub struct CifBlock {
    /// Block name (the string after `data_`).
    pub name: String,
    /// Key-value pairs from `_tag value` lines.
    pub tags: HashMap<String, String>,
    /// Atom site loop: each entry is (label, element, x, y, z, occupancy).
    pub atom_sites: Vec<CifAtomSite>,
}
/// A covalent angle triplet (atom indices, 1-based).
#[derive(Debug, Clone, Copy)]
pub struct PsfAngle {
    /// First atom.
    pub i: u32,
    /// Central atom.
    pub j: u32,
    /// Third atom.
    pub k: u32,
}
/// XSF periodicity descriptor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XsfPeriodicity {
    /// Isolated molecule (no periodic boundary conditions).
    Molecule,
    /// Polymer (1-D PBC).
    Polymer,
    /// Slab (2-D PBC).
    Slab,
    /// Crystal (3-D PBC).
    Crystal,
}
/// An atom entry in an XSF file.
#[derive(Debug, Clone)]
pub struct XsfAtom {
    /// Atomic number.
    pub atomic_number: i32,
    /// Position in angstroms.
    pub position: Vec3,
    /// Force vector (eV/Å) if present.
    pub force: Option<Vec3>,
}
/// Parsed contents of a Gaussian CUBE file.
#[derive(Debug, Clone)]
pub struct CubeFile {
    /// First comment line.
    pub comment1: String,
    /// Second comment line (often type of quantity).
    pub comment2: String,
    /// Number of atoms in the molecule.
    pub n_atoms: usize,
    /// Grid origin in Bohr.
    pub origin: Vec3,
    /// Number of voxels along x.
    pub nx: usize,
    /// Number of voxels along y.
    pub ny: usize,
    /// Number of voxels along z.
    pub nz: usize,
    /// Step vector along x (Bohr).
    pub dx: Vec3,
    /// Step vector along y (Bohr).
    pub dy: Vec3,
    /// Step vector along z (Bohr).
    pub dz: Vec3,
    /// Atoms: (atomic_number, charge, x, y, z) in Bohr.
    pub atoms: Vec<(i32, f64, Vec3)>,
    /// Volumetric data (flattened, fastest index = z).
    pub data: Vec<f64>,
}
impl CubeFile {
    /// Parse a CUBE file from text.
    pub fn parse(content: &str) -> Result<Self> {
        let mut lines = content.lines();
        let comment1 = lines.next().unwrap_or("").to_string();
        let comment2 = lines.next().unwrap_or("").to_string();
        let line3 = lines.next().ok_or_else(|| MolVizError::MissingSection {
            section: "natoms/origin".into(),
        })?;
        let toks3: Vec<&str> = line3.split_whitespace().collect();
        if toks3.len() < 4 {
            return Err(MolVizError::ParseError {
                message: "CUBE line 3 too short".into(),
            });
        }
        let n_atoms: usize = toks3[0].parse::<i32>().unwrap_or(0).unsigned_abs() as usize;
        let origin = [
            toks3[1].parse().unwrap_or(0.0),
            toks3[2].parse().unwrap_or(0.0),
            toks3[3].parse().unwrap_or(0.0),
        ];
        let parse_grid_line = |l: &str| -> (usize, Vec3) {
            let t: Vec<&str> = l.split_whitespace().collect();
            let n = t.first().and_then(|s| s.parse().ok()).unwrap_or(0usize);
            let v = [
                t.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0),
                t.get(2).and_then(|s| s.parse().ok()).unwrap_or(0.0),
                t.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.0),
            ];
            (n, v)
        };
        let (nx, dx) = parse_grid_line(lines.next().unwrap_or(""));
        let (ny, dy) = parse_grid_line(lines.next().unwrap_or(""));
        let (nz, dz) = parse_grid_line(lines.next().unwrap_or(""));
        let mut atoms = Vec::with_capacity(n_atoms);
        for _ in 0..n_atoms {
            let aline = lines.next().unwrap_or("");
            let at: Vec<&str> = aline.split_whitespace().collect();
            if at.len() < 5 {
                atoms.push((0, 0.0, [0.0; 3]));
                continue;
            }
            let atomic_num: i32 = at[0].parse().unwrap_or(0);
            let chg: f64 = at[1].parse().unwrap_or(0.0);
            let pos = [
                at[2].parse().unwrap_or(0.0),
                at[3].parse().unwrap_or(0.0),
                at[4].parse().unwrap_or(0.0),
            ];
            atoms.push((atomic_num, chg, pos));
        }
        let expected = nx * ny * nz;
        let mut data = Vec::with_capacity(expected);
        for line in lines {
            for tok in line.split_whitespace() {
                if let Ok(v) = tok.parse::<f64>() {
                    data.push(v);
                }
            }
        }
        data.resize(expected, 0.0);
        Ok(CubeFile {
            comment1,
            comment2,
            n_atoms,
            origin,
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            atoms,
            data,
        })
    }
    /// Electron density integrated over the entire grid (in Bohr³ × value units).
    pub fn integrate(&self) -> f64 {
        let vol = voxel_volume(self.dx, self.dy, self.dz);
        self.data.iter().sum::<f64>() * vol
    }
    /// Return the value at grid indices (ix, iy, iz).
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            self.data[ix * self.ny * self.nz + iy * self.nz + iz]
        } else {
            0.0
        }
    }
    /// Cartesian position (Bohr) for grid point (ix, iy, iz).
    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> Vec3 {
        let i = ix as f64;
        let j = iy as f64;
        let k = iz as f64;
        [
            self.origin[0] + i * self.dx[0] + j * self.dy[0] + k * self.dz[0],
            self.origin[1] + i * self.dx[1] + j * self.dy[1] + k * self.dz[1],
            self.origin[2] + i * self.dx[2] + j * self.dy[2] + k * self.dz[2],
        ]
    }
}
/// A single atom record extracted from a PSF file.
#[derive(Debug, Clone)]
pub struct PsfAtom {
    /// Sequential 1-based atom index as in the PSF file.
    pub serial: u32,
    /// Segment name (up to 4 characters).
    pub segname: String,
    /// Residue sequence number.
    pub resid: i32,
    /// Residue name (e.g. `"ALA"`).
    pub resname: String,
    /// Atom name within the residue (e.g. `"CA"`).
    pub atomname: String,
    /// CHARMM atom type (e.g. `"CT1"`).
    pub atom_type: String,
    /// Partial charge in units of electron charge.
    pub charge: f64,
    /// Mass in atomic mass units.
    pub mass: f64,
}
/// A molecular orbital entry from a Molden `[MO]` section.
#[derive(Debug, Clone)]
pub struct MoldenMo {
    /// Symmetry label (e.g. `"A1"`).
    pub symmetry: String,
    /// MO energy in Hartree.
    pub energy: f64,
    /// Spin: `"Alpha"` or `"Beta"`.
    pub spin: String,
    /// Occupation number.
    pub occupation: f64,
    /// Coefficients indexed by basis function (1-based → value).
    pub coefficients: Vec<(usize, f64)>,
}
/// A covalent bond between two atom indices (1-based, as in PSF).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PsfBond {
    /// First atom serial number.
    pub i: u32,
    /// Second atom serial number.
    pub j: u32,
}
/// A single atom with element symbol, position, and optional partial charge.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Element symbol (e.g. `"C"`, `"N"`, `"O"`).
    pub element: String,
    /// Cartesian coordinates \[x, y, z\].
    pub position: Vec3,
    /// Optional partial charge (e or esu depending on context).
    pub charge: Option<f64>,
    /// Optional atom name (e.g. `"CA"` for alpha-carbon in PDB/PSF notation).
    pub name: Option<String>,
}
impl Atom {
    /// Construct a minimal [`Atom`] from element and position.
    pub fn new(element: impl Into<String>, position: Vec3) -> Self {
        Self {
            element: element.into(),
            position,
            charge: None,
            name: None,
        }
    }
    /// Distance in angstroms to another atom.
    pub fn distance_to(&self, other: &Atom) -> f64 {
        let dx = self.position[0] - other.position[0];
        let dy = self.position[1] - other.position[1];
        let dz = self.position[2] - other.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}
