//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;

/// A sphere primitive representing one atom in ball-and-stick mode.
#[derive(Debug, Clone)]
pub struct AtomSphere {
    /// World-space centre of the sphere.
    pub center: [f64; 3],
    /// Sphere radius (scaled covalent radius).
    pub radius: f64,
    /// Display colour.
    pub color: MolColor,
    /// Source atom index.
    pub atom_index: usize,
}
/// A single (φ, ψ) backbone dihedral angle pair for one residue.
#[derive(Debug, Clone, Copy)]
pub struct DihedralPoint {
    /// φ angle in degrees (−180 to 180).
    pub phi: f64,
    /// ψ angle in degrees (−180 to 180).
    pub psi: f64,
    /// Residue sequence number.
    pub residue_seq: i32,
    /// Secondary structure at this residue.
    pub ss: SecondaryStructure,
}
/// Rendering data for a Ramachandran (φ/ψ) scatter plot.
#[derive(Debug, Clone)]
pub struct RamachandranData {
    /// All (φ, ψ) points in the protein.
    pub points: Vec<DihedralPoint>,
    /// 2-D density histogram for background colouring.
    ///
    /// `density_map[i * bins + j]` = count in bin (i_phi, i_psi).
    pub density_map: Vec<u32>,
    /// Number of bins per axis.
    pub bins: usize,
}
impl RamachandranData {
    /// Construct from a list of dihedral points.
    pub fn new(points: Vec<DihedralPoint>, bins: usize) -> Self {
        let b = bins.max(2);
        let mut density_map = vec![0u32; b * b];
        for p in &points {
            let i_phi = ((p.phi + 180.0) / 360.0 * b as f64).floor() as usize;
            let i_psi = ((p.psi + 180.0) / 360.0 * b as f64).floor() as usize;
            let i_phi = i_phi.min(b - 1);
            let i_psi = i_psi.min(b - 1);
            density_map[i_phi * b + i_psi] += 1;
        }
        Self {
            points,
            density_map,
            bins: b,
        }
    }
    /// Fraction of residues in the most-favoured region (α-helix and β-sheet).
    ///
    /// Uses empirical angle ranges for favoured regions.
    pub fn favoured_fraction(&self) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        let favoured = self
            .points
            .iter()
            .filter(|p| is_favoured(p.phi, p.psi))
            .count();
        favoured as f64 / self.points.len() as f64
    }
    /// Number of outlier points (not in favoured or allowed regions).
    pub fn outlier_count(&self) -> usize {
        self.points
            .iter()
            .filter(|p| !is_allowed(p.phi, p.psi))
            .count()
    }
    /// Separate points by secondary structure assignment.
    pub fn by_ss(&self, ss: SecondaryStructure) -> Vec<&DihedralPoint> {
        self.points.iter().filter(|p| p.ss == ss).collect()
    }
}
/// A covalent bond between two atoms.
#[derive(Debug, Clone)]
pub struct Bond {
    /// Index of the first atom.
    pub atom_a: usize,
    /// Index of the second atom.
    pub atom_b: usize,
    /// Bond order.
    pub order: BondOrder,
}
impl Bond {
    /// Construct a single bond between two atom indices.
    pub fn single(atom_a: usize, atom_b: usize) -> Self {
        Self {
            atom_a,
            atom_b,
            order: BondOrder::Single,
        }
    }
}
/// A sphere representing one atom in the CPK space-filling model.
#[derive(Debug, Clone)]
pub struct CpkSphere {
    /// World-space centre.
    pub center: [f64; 3],
    /// Van der Waals radius.
    pub radius: f64,
    /// CPK colour.
    pub color: MolColor,
    /// Source atom index.
    pub atom_index: usize,
}
/// Bond order classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BondOrder {
    /// Single bond.
    Single,
    /// Double bond.
    Double,
    /// Triple bond.
    Triple,
    /// Aromatic bond.
    Aromatic,
}
/// The seven crystal systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrystalSystem {
    /// Cubic (a=b=c, α=β=γ=90°).
    Cubic,
    /// Tetragonal (a=b≠c, α=β=γ=90°).
    Tetragonal,
    /// Orthorhombic (a≠b≠c, α=β=γ=90°).
    Orthorhombic,
    /// Hexagonal (a=b≠c, α=β=90°, γ=120°).
    Hexagonal,
    /// Trigonal (a=b=c, α=β=γ≠90°).
    Trigonal,
    /// Monoclinic (a≠b≠c, α=γ=90°, β≠90°).
    Monoclinic,
    /// Triclinic (a≠b≠c, α≠β≠γ).
    Triclinic,
}
/// A single triangle in a 3-D isosurface mesh.
#[derive(Debug, Clone, Copy)]
pub struct IsoTriangle {
    /// Vertex positions (Å).
    pub vertices: [[f64; 3]; 3],
    /// Per-vertex normals (unit vectors).
    pub normals: [[f64; 3]; 3],
    /// Isovalue level this triangle belongs to.
    pub level: f64,
}
impl IsoTriangle {
    /// Compute the face normal (cross product of two edges, normalised).
    pub fn face_normal(&self) -> [f64; 3] {
        let e0 = sub3(self.vertices[1], self.vertices[0]);
        let e1 = sub3(self.vertices[2], self.vertices[0]);
        normalize3(cross3(e0, e1))
    }
    /// Centroid of the triangle.
    pub fn centroid(&self) -> [f64; 3] {
        scale3(
            add3(add3(self.vertices[0], self.vertices[1]), self.vertices[2]),
            1.0 / 3.0,
        )
    }
}
/// Radial distribution function `g(r)` plot data.
#[derive(Debug, Clone)]
pub struct RadialDistribution {
    /// Atom type pair label (e.g. `"O-H"`).
    pub pair_label: String,
    /// Bin centres in Ångström.
    pub r: Vec<f64>,
    /// `g(r)` values at each bin.
    pub g_r: Vec<f64>,
    /// Bin width in Ångström.
    pub dr: f64,
    /// Number density of the system (atoms/Å³) used for normalisation.
    pub number_density: f64,
}
impl RadialDistribution {
    /// Compute the RDF from a single trajectory frame.
    ///
    /// - `positions` — all atom positions.
    /// - `indices_a`, `indices_b` — indices of the two atom types.
    /// - `r_max` — maximum radius (Å).
    /// - `n_bins` — number of histogram bins.
    /// - `number_density` — bulk number density (Å⁻³).
    /// - `pair_label` — display label.
    pub fn compute(
        positions: &[[f64; 3]],
        indices_a: &[usize],
        indices_b: &[usize],
        r_max: f64,
        n_bins: usize,
        number_density: f64,
        pair_label: impl Into<String>,
    ) -> Self {
        let dr = r_max / n_bins as f64;
        let mut hist = vec![0u64; n_bins];
        let n_a = indices_a.len();
        let n_b = indices_b.len();
        for &ia in indices_a {
            for &ib in indices_b {
                if ia == ib {
                    continue;
                }
                if let (Some(pa), Some(pb)) = (positions.get(ia), positions.get(ib)) {
                    let d = len3(sub3(*pa, *pb));
                    if d < r_max {
                        let bin = (d / dr) as usize;
                        if bin < n_bins {
                            hist[bin] += 1;
                        }
                    }
                }
            }
        }
        let n_pairs = (n_a * n_b) as f64;
        let r: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * dr).collect();
        let g_r: Vec<f64> = hist
            .iter()
            .enumerate()
            .map(|(i, &count)| {
                let r_lo = i as f64 * dr;
                let r_hi = r_lo + dr;
                let shell_vol = (4.0 / 3.0) * PI * (r_hi.powi(3) - r_lo.powi(3));
                let ideal = number_density * shell_vol * n_pairs;
                if ideal < 1e-30 {
                    0.0
                } else {
                    count as f64 / ideal
                }
            })
            .collect();
        Self {
            pair_label: pair_label.into(),
            r,
            g_r,
            dr,
            number_density,
        }
    }
    /// Index of the first peak in `g(r)` (largest value in the RDF).
    pub fn first_peak_index(&self) -> Option<usize> {
        self.g_r
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
    }
    /// `r` value at the first peak.
    pub fn first_peak_r(&self) -> Option<f64> {
        self.first_peak_index().and_then(|i| self.r.get(i).copied())
    }
}
/// Symmetry type of a molecular orbital.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrbitalSymmetry {
    /// Sigma (σ) orbital.
    Sigma,
    /// Pi (π) orbital.
    Pi,
    /// Delta (δ) orbital.
    Delta,
    /// Non-bonding (n) orbital.
    NonBonding,
}
/// CPK (space-filling) scene.
#[derive(Debug, Clone)]
pub struct CpkScene {
    /// CPK spheres — one per atom, using van der Waals radii.
    pub spheres: Vec<CpkSphere>,
}
impl CpkScene {
    /// Build a CPK scene from a list of atoms.
    ///
    /// An optional `vdw_scale` factor (default 1.0) uniformly scales all radii.
    pub fn build(atoms: &[Atom], vdw_scale: f64) -> Self {
        let spheres = atoms
            .iter()
            .map(|a| CpkSphere {
                center: a.position,
                radius: a.element.vdw_radius() * vdw_scale,
                color: a.element.cpk_color(),
                atom_index: a.index,
            })
            .collect();
        Self { spheres }
    }
    /// Return the total solvent-accessible volume approximation (sum of sphere volumes).
    pub fn total_volume(&self) -> f64 {
        self.spheres
            .iter()
            .map(|s| (4.0 / 3.0) * PI * s.radius.powi(3))
            .sum()
    }
}
/// A single residue control point for the ribbon spline.
#[derive(Debug, Clone)]
pub struct RibbonControlPoint {
    /// Position of the α-carbon (Cα) in Ångström.
    pub ca_position: [f64; 3],
    /// Direction vector of the ribbon at this residue (roughly in-plane normal).
    pub ribbon_normal: [f64; 3],
    /// Secondary structure assignment.
    pub ss: SecondaryStructure,
    /// Residue sequence number.
    pub residue_seq: i32,
}
/// Combined scene containing all molecular rendering primitives.
#[derive(Debug, Clone)]
pub struct MolecularScene {
    /// Ball-and-stick representation (optional).
    pub ball_stick: Option<BallStickScene>,
    /// CPK space-filling representation (optional).
    pub cpk: Option<CpkScene>,
    /// Ribbon diagram (optional, for proteins).
    pub ribbon: Option<RibbonDiagram>,
    /// Electron density isosurface triangles (optional).
    pub density_iso: Vec<IsoTriangle>,
    /// Molecular orbital isosurface triangles: positive phase.
    pub mo_pos: Vec<IsoTriangle>,
    /// Molecular orbital isosurface triangles: negative phase.
    pub mo_neg: Vec<IsoTriangle>,
    /// Crystal cell wireframe edges (optional).
    pub cell_edges: Vec<CellEdge>,
    /// Trajectory metadata (used by external playback engine).
    pub trajectory: Option<TrajectoryAnimation>,
}
impl MolecularScene {
    /// Construct an empty molecular scene.
    pub fn empty() -> Self {
        Self {
            ball_stick: None,
            cpk: None,
            ribbon: None,
            density_iso: Vec::new(),
            mo_pos: Vec::new(),
            mo_neg: Vec::new(),
            cell_edges: Vec::new(),
            trajectory: None,
        }
    }
    /// Total number of rendering primitives in the scene.
    pub fn primitive_count(&self) -> usize {
        let bs = self.ball_stick.as_ref().map_or(0, |s| s.primitive_count());
        let cpk = self.cpk.as_ref().map_or(0, |s| s.spheres.len());
        let rib = self.ribbon.as_ref().map_or(0, |r| r.segments.len());
        bs + cpk
            + rib
            + self.density_iso.len()
            + self.mo_pos.len()
            + self.mo_neg.len()
            + self.cell_edges.len()
    }
}
/// Builder for spatial / chemical atom queries.
#[derive(Debug, Clone)]
pub struct AtomSelector {
    /// Maximum distance from `center` (Å).
    pub distance_cutoff: Option<f64>,
    /// Centre for distance query.
    pub center: Option<[f64; 3]>,
    /// Element whitelist.
    pub elements: Vec<Element>,
    /// Residue index range `(lo, hi)` inclusive.
    pub residue_range: Option<(usize, usize)>,
}
impl AtomSelector {
    /// Construct an empty (select-all) selector.
    pub fn new() -> Self {
        Self {
            distance_cutoff: None,
            center: None,
            elements: Vec::new(),
            residue_range: None,
        }
    }
    /// Restrict to atoms within `radius` Å of `center`.
    pub fn within(mut self, center: [f64; 3], radius: f64) -> Self {
        self.center = Some(center);
        self.distance_cutoff = Some(radius);
        self
    }
    /// Restrict to specific elements.
    pub fn by_elements(mut self, elems: &[Element]) -> Self {
        self.elements = elems.to_vec();
        self
    }
    /// Restrict to residue indices in `[lo, hi]`.
    pub fn by_residues(mut self, lo: usize, hi: usize) -> Self {
        self.residue_range = Some((lo, hi));
        self
    }
    /// Apply the selector to a list of atoms, returning matching indices.
    pub fn select(&self, atoms: &[Atom]) -> Vec<usize> {
        atoms
            .iter()
            .enumerate()
            .filter_map(|(i, a)| {
                if !self.elements.is_empty() && !self.elements.contains(&a.element) {
                    return None;
                }
                if let (Some(ctr), Some(cutoff)) = (self.center, self.distance_cutoff)
                    && len3(sub3(a.position, ctr)) > cutoff
                {
                    return None;
                }
                if let Some((lo, hi)) = self.residue_range {
                    match a.residue_index {
                        Some(ri) if ri >= lo && ri <= hi => {}
                        _ => return None,
                    }
                }
                Some(i)
            })
            .collect()
    }
}
/// A molecular dynamics trajectory: ordered sequence of frames.
#[derive(Debug, Clone)]
pub struct TrajectoryAnimation {
    /// Number of atoms in the system (fixed across frames).
    pub n_atoms: usize,
    /// Trajectory frames in chronological order.
    pub frames: Vec<TrajectoryFrame>,
    /// Frame index currently selected for display.
    pub current_frame: usize,
    /// Playback speed in frames per second.
    pub playback_fps: f64,
    /// Whether the animation loops.
    pub looping: bool,
}
impl TrajectoryAnimation {
    /// Construct a new empty trajectory.
    pub fn new(n_atoms: usize, playback_fps: f64) -> Self {
        Self {
            n_atoms,
            frames: Vec::new(),
            current_frame: 0,
            playback_fps,
            looping: true,
        }
    }
    /// Append a frame to the trajectory.
    pub fn push_frame(&mut self, frame: TrajectoryFrame) {
        self.frames.push(frame);
    }
    /// Simulation duration in picoseconds (last frame time − first frame time).
    pub fn duration_ps(&self) -> f64 {
        match (self.frames.first(), self.frames.last()) {
            (Some(f), Some(l)) => l.time_ps - f.time_ps,
            _ => 0.0,
        }
    }
    /// Advance to the next frame, wrapping if `looping` is enabled.
    ///
    /// Returns `true` if the animation wrapped around.
    pub fn advance(&mut self) -> bool {
        if self.frames.is_empty() {
            return false;
        }
        self.current_frame += 1;
        if self.current_frame >= self.frames.len() {
            if self.looping {
                self.current_frame = 0;
                return true;
            } else {
                self.current_frame = self.frames.len() - 1;
            }
        }
        false
    }
    /// Return the current frame, or `None` if there are no frames.
    pub fn current(&self) -> Option<&TrajectoryFrame> {
        self.frames.get(self.current_frame)
    }
    /// Mean total energy over all frames that have it available.
    pub fn mean_total_energy(&self) -> Option<f64> {
        let energies: Vec<f64> = self
            .frames
            .iter()
            .filter_map(|f| f.total_energy())
            .collect();
        if energies.is_empty() {
            None
        } else {
            Some(energies.iter().sum::<f64>() / energies.len() as f64)
        }
    }
}
/// Crystal unit cell visualization data.
#[derive(Debug, Clone)]
pub struct CrystalCell {
    /// Lattice parameters.
    pub params: LatticeParameters,
    /// Origin of the unit cell (Å, Cartesian).
    pub origin: [f64; 3],
    /// Atoms in the asymmetric unit (fractional coordinates).
    pub asymmetric_atoms: Vec<AtomFractional>,
    /// All symmetry-expanded atoms in the unit cell (Cartesian Å).
    pub unit_cell_atoms: Vec<Atom>,
}
impl CrystalCell {
    /// Construct a crystal cell from lattice parameters and fractional-coordinate atoms.
    pub fn new(
        params: LatticeParameters,
        origin: [f64; 3],
        asymmetric_atoms: Vec<AtomFractional>,
    ) -> Self {
        let (va, vb, vc) = params.lattice_vectors();
        let unit_cell_atoms: Vec<Atom> = asymmetric_atoms
            .iter()
            .enumerate()
            .map(|(i, af)| {
                let pos = [
                    origin[0] + af.fa * va[0] + af.fb * vb[0] + af.fc * vc[0],
                    origin[1] + af.fa * va[1] + af.fb * vb[1] + af.fc * vc[1],
                    origin[2] + af.fa * va[2] + af.fb * vb[2] + af.fc * vc[2],
                ];
                Atom::new(i, pos, af.element)
            })
            .collect();
        Self {
            params,
            origin,
            asymmetric_atoms,
            unit_cell_atoms,
        }
    }
    /// Generate the 12 wireframe edges of the unit cell.
    ///
    /// Uses the CPK colour of each axis vector for visual differentiation.
    pub fn wireframe_edges(&self) -> Vec<CellEdge> {
        let (va, vb, vc) = self.params.lattice_vectors();
        let o = self.origin;
        let corners = [
            o,
            add3(o, va),
            add3(o, vb),
            add3(o, vc),
            add3(add3(o, va), vb),
            add3(add3(o, va), vc),
            add3(add3(o, vb), vc),
            add3(add3(add3(o, va), vb), vc),
        ];
        let red = MolColor::new(0.85, 0.25, 0.25, 1.0);
        let green = MolColor::new(0.25, 0.75, 0.30, 1.0);
        let blue = MolColor::new(0.25, 0.40, 0.90, 1.0);
        vec![
            CellEdge {
                start: corners[0],
                end: corners[1],
                color: red,
            },
            CellEdge {
                start: corners[2],
                end: corners[4],
                color: red,
            },
            CellEdge {
                start: corners[3],
                end: corners[5],
                color: red,
            },
            CellEdge {
                start: corners[6],
                end: corners[7],
                color: red,
            },
            CellEdge {
                start: corners[0],
                end: corners[2],
                color: green,
            },
            CellEdge {
                start: corners[1],
                end: corners[4],
                color: green,
            },
            CellEdge {
                start: corners[3],
                end: corners[6],
                color: green,
            },
            CellEdge {
                start: corners[5],
                end: corners[7],
                color: green,
            },
            CellEdge {
                start: corners[0],
                end: corners[3],
                color: blue,
            },
            CellEdge {
                start: corners[1],
                end: corners[5],
                color: blue,
            },
            CellEdge {
                start: corners[2],
                end: corners[6],
                color: blue,
            },
            CellEdge {
                start: corners[4],
                end: corners[7],
                color: blue,
            },
        ]
    }
    /// Generate a supercell by tiling `(na, nb, nc)` times along each axis.
    pub fn supercell_atoms(&self, na: usize, nb: usize, nc: usize) -> Vec<Atom> {
        let (va, vb, vc) = self.params.lattice_vectors();
        let mut atoms = Vec::new();
        let mut idx = 0;
        for ia in 0..na {
            for ib in 0..nb {
                for ic in 0..nc {
                    let shift = [
                        ia as f64 * va[0] + ib as f64 * vb[0] + ic as f64 * vc[0],
                        ia as f64 * va[1] + ib as f64 * vb[1] + ic as f64 * vc[1],
                        ia as f64 * va[2] + ib as f64 * vb[2] + ic as f64 * vc[2],
                    ];
                    for a in &self.unit_cell_atoms {
                        atoms.push(Atom {
                            index: idx,
                            position: add3(a.position, shift),
                            element: a.element,
                            partial_charge: a.partial_charge,
                            name: a.name.clone(),
                            residue_index: a.residue_index,
                        });
                        idx += 1;
                    }
                }
            }
        }
        atoms
    }
}
/// 3-D regular grid of electron density values (e/Å³).
#[derive(Debug, Clone)]
pub struct ElectronDensityGrid {
    /// Number of grid points along X.
    pub nx: usize,
    /// Number of grid points along Y.
    pub ny: usize,
    /// Number of grid points along Z.
    pub nz: usize,
    /// Grid spacing in Ångström (uniform in all directions).
    pub spacing: f64,
    /// Origin of the grid (Å).
    pub origin: [f64; 3],
    /// Flat array of density values in row-major (x, y, z) order.
    pub values: Vec<f64>,
}
impl ElectronDensityGrid {
    /// Construct a new zero-filled grid.
    pub fn new(nx: usize, ny: usize, nz: usize, spacing: f64, origin: [f64; 3]) -> Self {
        Self {
            nx,
            ny,
            nz,
            spacing,
            origin,
            values: vec![0.0; nx * ny * nz],
        }
    }
    /// Flat index for grid point `(ix, iy, iz)`.
    pub fn index(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.ny * self.nz + iy * self.nz + iz
    }
    /// World-space position of grid point `(ix, iy, iz)`.
    pub fn position(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.spacing,
            self.origin[1] + iy as f64 * self.spacing,
            self.origin[2] + iz as f64 * self.spacing,
        ]
    }
    /// Get the density at grid point `(ix, iy, iz)`.
    ///
    /// Returns `0.0` if any index is out of range.
    pub fn get(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            self.values[self.index(ix, iy, iz)]
        } else {
            0.0
        }
    }
    /// Set the density at grid point `(ix, iy, iz)`.
    pub fn set(&mut self, ix: usize, iy: usize, iz: usize, val: f64) {
        if ix < self.nx && iy < self.ny && iz < self.nz {
            let idx = self.index(ix, iy, iz);
            self.values[idx] = val;
        }
    }
    /// Maximum density value in the grid.
    pub fn max_value(&self) -> f64 {
        self.values
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Minimum density value in the grid.
    pub fn min_value(&self) -> f64 {
        self.values.iter().cloned().fold(f64::INFINITY, f64::min)
    }
    /// Extract a flat slice at a fixed Z index as a `(nx × ny)` Vec.
    pub fn slice_z(&self, iz: usize) -> Vec<f64> {
        let mut out = Vec::with_capacity(self.nx * self.ny);
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                out.push(self.get(ix, iy, iz));
            }
        }
        out
    }
    /// Tri-linearly interpolate density at a world-space position.
    pub fn interpolate(&self, pos: [f64; 3]) -> f64 {
        let fx = (pos[0] - self.origin[0]) / self.spacing;
        let fy = (pos[1] - self.origin[1]) / self.spacing;
        let fz = (pos[2] - self.origin[2]) / self.spacing;
        let ix = fx.floor() as isize;
        let iy = fy.floor() as isize;
        let iz = fz.floor() as isize;
        let tx = fx - fx.floor();
        let ty = fy - fy.floor();
        let tz = fz - fz.floor();
        let get = |i: isize, j: isize, k: isize| -> f64 {
            if i < 0 || j < 0 || k < 0 {
                return 0.0;
            }
            self.get(i as usize, j as usize, k as usize)
        };
        let c000 = get(ix, iy, iz);
        let c100 = get(ix + 1, iy, iz);
        let c010 = get(ix, iy + 1, iz);
        let c110 = get(ix + 1, iy + 1, iz);
        let c001 = get(ix, iy, iz + 1);
        let c101 = get(ix + 1, iy, iz + 1);
        let c011 = get(ix, iy + 1, iz + 1);
        let c111 = get(ix + 1, iy + 1, iz + 1);
        let c00 = lerp(c000, c100, tx);
        let c01 = lerp(c001, c101, tx);
        let c10 = lerp(c010, c110, tx);
        let c11 = lerp(c011, c111, tx);
        let c0 = lerp(c00, c10, ty);
        let c1 = lerp(c01, c11, ty);
        lerp(c0, c1, tz)
    }
    /// Generate a list of isosurface seed points above a given density threshold.
    ///
    /// Returns grid positions where the density exceeds `isovalue`.
    pub fn isosurface_seeds(&self, isovalue: f64) -> Vec<[f64; 3]> {
        let mut seeds = Vec::new();
        for ix in 0..self.nx {
            for iy in 0..self.ny {
                for iz in 0..self.nz {
                    if self.get(ix, iy, iz) >= isovalue {
                        seeds.push(self.position(ix, iy, iz));
                    }
                }
            }
        }
        seeds
    }
}
/// A cylinder primitive representing one bond in ball-and-stick mode.
///
/// The cylinder spans from `start` to `end` with the given radius and colour.
#[derive(Debug, Clone)]
pub struct BondCylinder {
    /// Start point (near atom A).
    pub start: [f64; 3],
    /// End point (near atom B).
    pub end: [f64; 3],
    /// Cylinder radius in Ångström.
    pub radius: f64,
    /// Colour of this half of the bond (half-bond colouring by element).
    pub color: MolColor,
}
/// A single 2-D contour line segment from marching-squares.
#[derive(Debug, Clone, Copy)]
pub struct ContourSegment {
    /// First endpoint (x, y) in data-space.
    pub p0: [f64; 2],
    /// Second endpoint (x, y) in data-space.
    pub p1: [f64; 2],
    /// Contour level.
    pub level: f64,
}
/// A 2-D potential energy surface (PES) on a regular grid.
#[derive(Debug, Clone)]
pub struct EnergyLandscape {
    /// Axis label for the first collective variable (CV1).
    pub cv1_label: String,
    /// Axis label for the second collective variable (CV2).
    pub cv2_label: String,
    /// CV1 axis values (uniformly spaced).
    pub cv1: Vec<f64>,
    /// CV2 axis values (uniformly spaced).
    pub cv2: Vec<f64>,
    /// Energy values in kJ/mol, row-major `[i_cv1 * n_cv2 + i_cv2]`.
    pub energy: Vec<f64>,
    /// Number of grid points along CV1.
    pub n_cv1: usize,
    /// Number of grid points along CV2.
    pub n_cv2: usize,
}
impl EnergyLandscape {
    /// Build from sampled CV data using a 2-D histogram / free energy estimator.
    ///
    /// - `cv1_data`, `cv2_data` — paired CV samples.
    /// - `cv1_range` / `cv2_range` — `(min, max)` of each axis.
    /// - `n_cv1`, `n_cv2` — number of bins per axis.
    /// - `temperature` — temperature in K (used for k_BT units).
    /// - Labels for each axis.
    pub fn from_samples(
        cv1_data: &[f64],
        cv2_data: &[f64],
        cv1_range: (f64, f64),
        cv2_range: (f64, f64),
        n_cv1: usize,
        n_cv2: usize,
        temperature: f64,
        cv1_label: impl Into<String>,
        cv2_label: impl Into<String>,
    ) -> Self {
        let kbt = 8.314e-3 * temperature;
        let n1 = n_cv1.max(2);
        let n2 = n_cv2.max(2);
        let d1 = (cv1_range.1 - cv1_range.0) / n1 as f64;
        let d2 = (cv2_range.1 - cv2_range.0) / n2 as f64;
        let mut counts = vec![0u64; n1 * n2];
        for (&v1, &v2) in cv1_data.iter().zip(cv2_data.iter()) {
            let i1 = ((v1 - cv1_range.0) / d1).floor() as isize;
            let i2 = ((v2 - cv2_range.0) / d2).floor() as isize;
            if i1 >= 0 && i1 < n1 as isize && i2 >= 0 && i2 < n2 as isize {
                counts[i1 as usize * n2 + i2 as usize] += 1;
            }
        }
        let total: u64 = counts.iter().sum();
        let total_f = total as f64;
        let energy: Vec<f64> = counts
            .iter()
            .map(|&c| {
                if c == 0 || total == 0 {
                    f64::INFINITY
                } else {
                    -kbt * ((c as f64) / total_f).ln()
                }
            })
            .collect();
        let e_min = energy
            .iter()
            .cloned()
            .filter(|e| e.is_finite())
            .fold(f64::INFINITY, f64::min);
        let energy: Vec<f64> = energy
            .iter()
            .map(|&e| if e.is_finite() { e - e_min } else { e })
            .collect();
        let cv1: Vec<f64> = (0..n1)
            .map(|i| cv1_range.0 + (i as f64 + 0.5) * d1)
            .collect();
        let cv2: Vec<f64> = (0..n2)
            .map(|i| cv2_range.0 + (i as f64 + 0.5) * d2)
            .collect();
        Self {
            cv1_label: cv1_label.into(),
            cv2_label: cv2_label.into(),
            cv1,
            cv2,
            energy,
            n_cv1: n1,
            n_cv2: n2,
        }
    }
    /// Energy at grid point `(i1, i2)`.
    pub fn get(&self, i1: usize, i2: usize) -> f64 {
        if i1 < self.n_cv1 && i2 < self.n_cv2 {
            self.energy[i1 * self.n_cv2 + i2]
        } else {
            f64::INFINITY
        }
    }
    /// Minimum finite energy and its grid indices.
    pub fn minimum(&self) -> Option<(f64, usize, usize)> {
        self.energy
            .iter()
            .enumerate()
            .filter(|&(_, &e)| e.is_finite())
            .min_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(idx, &e)| (e, idx / self.n_cv2, idx % self.n_cv2))
    }
    /// Generate contour levels at evenly-spaced energies between 0 and `e_max`.
    pub fn contour_levels(&self, n_levels: usize, e_max: f64) -> Vec<f64> {
        let n = n_levels.max(2);
        (0..n).map(|i| e_max * i as f64 / (n - 1) as f64).collect()
    }
}
/// RGBA colour with `f32` components in `[0, 1]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MolColor {
    /// Red component.
    pub r: f32,
    /// Green component.
    pub g: f32,
    /// Blue component.
    pub b: f32,
    /// Alpha component.
    pub a: f32,
}
impl MolColor {
    /// Construct from individual RGBA components.
    pub fn new(r: f32, g: f32, b: f32, a: f32) -> Self {
        Self { r, g, b, a }
    }
    /// Opaque white.
    pub fn white() -> Self {
        Self::new(1.0, 1.0, 1.0, 1.0)
    }
    /// Opaque black.
    pub fn black() -> Self {
        Self::new(0.0, 0.0, 0.0, 1.0)
    }
    /// Linearly interpolate between two colours.
    pub fn blend(a: Self, b: Self, t: f32) -> Self {
        let t = t.clamp(0.0, 1.0);
        Self::new(
            a.r + (b.r - a.r) * t,
            a.g + (b.g - a.g) * t,
            a.b + (b.b - a.b) * t,
            a.a + (b.a - a.a) * t,
        )
    }
}
/// Complete ball-and-stick scene ready for rendering.
#[derive(Debug, Clone)]
pub struct BallStickScene {
    /// Atom spheres.
    pub spheres: Vec<AtomSphere>,
    /// Bond cylinders (each bond produces two half-cylinders).
    pub cylinders: Vec<BondCylinder>,
}
impl BallStickScene {
    /// Build a ball-and-stick scene from atoms and bonds.
    ///
    /// - `atom_scale` — scaling factor applied to covalent radii (typically 0.3–0.5).
    /// - `bond_radius` — cylinder radius in Ångström (typically 0.10–0.15).
    pub fn build(atoms: &[Atom], bonds: &[Bond], atom_scale: f64, bond_radius: f64) -> Self {
        let spheres: Vec<AtomSphere> = atoms
            .iter()
            .map(|a| AtomSphere {
                center: a.position,
                radius: a.element.covalent_radius() * atom_scale,
                color: a.element.cpk_color(),
                atom_index: a.index,
            })
            .collect();
        let mut cylinders = Vec::with_capacity(bonds.len() * 2);
        for bond in bonds {
            let Some(a) = atoms.iter().find(|x| x.index == bond.atom_a) else {
                continue;
            };
            let Some(b) = atoms.iter().find(|x| x.index == bond.atom_b) else {
                continue;
            };
            let mid = scale3(add3(a.position, b.position), 0.5);
            cylinders.push(BondCylinder {
                start: a.position,
                end: mid,
                radius: bond_radius,
                color: a.element.cpk_color(),
            });
            cylinders.push(BondCylinder {
                start: mid,
                end: b.position,
                radius: bond_radius,
                color: b.element.cpk_color(),
            });
        }
        Self { spheres, cylinders }
    }
    /// Total number of primitives (spheres + cylinders).
    pub fn primitive_count(&self) -> usize {
        self.spheres.len() + self.cylinders.len()
    }
    /// Compute the axis-aligned bounding box `(min, max)` of all atom positions.
    pub fn aabb(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.spheres.is_empty() {
            return None;
        }
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for s in &self.spheres {
            for k in 0..3 {
                if s.center[k] - s.radius < mn[k] {
                    mn[k] = s.center[k] - s.radius;
                }
                if s.center[k] + s.radius > mx[k] {
                    mx[k] = s.center[k] + s.radius;
                }
            }
        }
        Some((mn, mx))
    }
}
/// A single atom with position, element, and optional partial charge.
#[derive(Debug, Clone)]
pub struct Atom {
    /// Atom index within its parent molecule.
    pub index: usize,
    /// 3-D position in Ångström.
    pub position: [f64; 3],
    /// Element type.
    pub element: Element,
    /// Optional partial charge (e).
    pub partial_charge: Option<f64>,
    /// Optional atom name (e.g. `"CA"` for α-carbon in PDB files).
    pub name: String,
    /// Optional residue index for protein atoms.
    pub residue_index: Option<usize>,
}
impl Atom {
    /// Construct a new atom.
    pub fn new(index: usize, position: [f64; 3], element: Element) -> Self {
        Self {
            index,
            position,
            element,
            partial_charge: None,
            name: String::new(),
            residue_index: None,
        }
    }
}
/// Secondary structure class for a residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SecondaryStructure {
    /// α-helix.
    Helix,
    /// β-strand.
    Strand,
    /// Loop / coil.
    Loop,
    /// Turn.
    Turn,
}
impl SecondaryStructure {
    /// Default ribbon colour for this secondary structure type.
    pub fn ribbon_color(self) -> MolColor {
        match self {
            SecondaryStructure::Helix => MolColor::new(0.88, 0.30, 0.30, 1.0),
            SecondaryStructure::Strand => MolColor::new(0.25, 0.55, 0.90, 1.0),
            SecondaryStructure::Loop => MolColor::new(0.60, 0.85, 0.50, 1.0),
            SecondaryStructure::Turn => MolColor::new(0.95, 0.80, 0.30, 1.0),
        }
    }
}
/// A ribbon segment connecting two consecutive control points.
#[derive(Debug, Clone)]
pub struct RibbonSegment {
    /// Start position.
    pub start: [f64; 3],
    /// End position.
    pub end: [f64; 3],
    /// Ribbon width at start.
    pub width_start: f64,
    /// Ribbon width at end.
    pub width_end: f64,
    /// Colour at start.
    pub color_start: MolColor,
    /// Colour at end.
    pub color_end: MolColor,
    /// Normal direction for lighting.
    pub normal: [f64; 3],
}
/// Periodic-table element identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Element {
    /// Hydrogen (H).
    H,
    /// Carbon (C).
    C,
    /// Nitrogen (N).
    N,
    /// Oxygen (O).
    O,
    /// Sulfur (S).
    S,
    /// Phosphorus (P).
    P,
    /// Fluorine (F).
    F,
    /// Chlorine (Cl).
    Cl,
    /// Bromine (Br).
    Br,
    /// Iodine (I).
    I,
    /// Iron (Fe).
    Fe,
    /// Calcium (Ca).
    Ca,
    /// Zinc (Zn).
    Zn,
    /// Magnesium (Mg).
    Mg,
    /// Unknown/generic element.
    Unknown,
}
impl Element {
    /// Covalent radius in Ångström for ball-and-stick rendering.
    pub fn covalent_radius(self) -> f64 {
        match self {
            Element::H => 0.31,
            Element::C => 0.77,
            Element::N => 0.75,
            Element::O => 0.73,
            Element::S => 1.02,
            Element::P => 1.06,
            Element::F => 0.71,
            Element::Cl => 0.99,
            Element::Br => 1.14,
            Element::I => 1.33,
            Element::Fe => 1.32,
            Element::Ca => 1.76,
            Element::Zn => 1.22,
            Element::Mg => 1.41,
            Element::Unknown => 0.80,
        }
    }
    /// Van der Waals radius in Ångström for CPK space-filling rendering.
    pub fn vdw_radius(self) -> f64 {
        match self {
            Element::H => 1.20,
            Element::C => 1.70,
            Element::N => 1.55,
            Element::O => 1.52,
            Element::S => 1.80,
            Element::P => 1.80,
            Element::F => 1.47,
            Element::Cl => 1.75,
            Element::Br => 1.85,
            Element::I => 1.98,
            Element::Fe => 2.05,
            Element::Ca => 2.31,
            Element::Zn => 2.10,
            Element::Mg => 1.73,
            Element::Unknown => 1.70,
        }
    }
    /// Standard CPK colour for this element.
    pub fn cpk_color(self) -> MolColor {
        match self {
            Element::H => MolColor::new(1.00, 1.00, 1.00, 1.0),
            Element::C => MolColor::new(0.33, 0.33, 0.33, 1.0),
            Element::N => MolColor::new(0.13, 0.47, 1.00, 1.0),
            Element::O => MolColor::new(0.91, 0.15, 0.15, 1.0),
            Element::S => MolColor::new(1.00, 0.85, 0.00, 1.0),
            Element::P => MolColor::new(1.00, 0.50, 0.00, 1.0),
            Element::F => MolColor::new(0.56, 0.82, 0.32, 1.0),
            Element::Cl => MolColor::new(0.12, 0.94, 0.12, 1.0),
            Element::Br => MolColor::new(0.65, 0.16, 0.16, 1.0),
            Element::I => MolColor::new(0.58, 0.00, 0.58, 1.0),
            Element::Fe => MolColor::new(0.88, 0.40, 0.20, 1.0),
            Element::Ca => MolColor::new(0.24, 1.00, 0.00, 1.0),
            Element::Zn => MolColor::new(0.49, 0.50, 0.69, 1.0),
            Element::Mg => MolColor::new(0.54, 1.00, 0.00, 1.0),
            Element::Unknown => MolColor::new(0.80, 0.80, 0.80, 1.0),
        }
    }
    /// Parse an element symbol string (case-insensitive).
    pub fn from_symbol(sym: &str) -> Self {
        match sym.trim() {
            "H" | "h" => Element::H,
            "C" | "c" => Element::C,
            "N" | "n" => Element::N,
            "O" | "o" => Element::O,
            "S" | "s" => Element::S,
            "P" | "p" => Element::P,
            "F" | "f" => Element::F,
            "Cl" | "CL" | "cl" => Element::Cl,
            "Br" | "BR" | "br" => Element::Br,
            "I" | "i" => Element::I,
            "Fe" | "FE" | "fe" => Element::Fe,
            "Ca" | "CA" | "ca" => Element::Ca,
            "Zn" | "ZN" | "zn" => Element::Zn,
            "Mg" | "MG" | "mg" => Element::Mg,
            _ => Element::Unknown,
        }
    }
}
/// Complete ribbon diagram for a protein chain.
#[derive(Debug, Clone)]
pub struct RibbonDiagram {
    /// Ribbon segments (one per consecutive residue pair after smoothing).
    pub segments: Vec<RibbonSegment>,
    /// Number of residues used to build the diagram.
    pub residue_count: usize,
}
impl RibbonDiagram {
    /// Build a ribbon diagram from control points.
    ///
    /// - `smooth_steps` — number of sub-steps for Catmull-Rom interpolation per segment.
    /// - `helix_width`, `strand_width`, `loop_width` — ribbon widths in Å for each SS type.
    pub fn build(
        control_points: &[RibbonControlPoint],
        smooth_steps: usize,
        helix_width: f64,
        strand_width: f64,
        loop_width: f64,
    ) -> Self {
        if control_points.len() < 2 {
            return Self {
                segments: vec![],
                residue_count: control_points.len(),
            };
        }
        let n = control_points.len();
        let steps = smooth_steps.max(1);
        let mut segments = Vec::with_capacity(n * steps);
        for i in 0..n.saturating_sub(1) {
            let p0 = if i == 0 {
                &control_points[0]
            } else {
                &control_points[i - 1]
            };
            let p1 = &control_points[i];
            let p2 = &control_points[i + 1];
            let p3 = if i + 2 < n {
                &control_points[i + 2]
            } else {
                &control_points[n - 1]
            };
            for step in 0..steps {
                let t0 = step as f64 / steps as f64;
                let t1 = (step + 1) as f64 / steps as f64;
                let pos0 = catmull_rom(
                    p0.ca_position,
                    p1.ca_position,
                    p2.ca_position,
                    p3.ca_position,
                    t0,
                );
                let pos1 = catmull_rom(
                    p0.ca_position,
                    p1.ca_position,
                    p2.ca_position,
                    p3.ca_position,
                    t1,
                );
                let norm0 = normalize3(catmull_rom(
                    p0.ribbon_normal,
                    p1.ribbon_normal,
                    p2.ribbon_normal,
                    p3.ribbon_normal,
                    t0,
                ));
                let width_start = lerp(
                    ss_width(p1.ss, helix_width, strand_width, loop_width),
                    ss_width(p2.ss, helix_width, strand_width, loop_width),
                    t0,
                );
                let width_end = lerp(
                    ss_width(p1.ss, helix_width, strand_width, loop_width),
                    ss_width(p2.ss, helix_width, strand_width, loop_width),
                    t1,
                );
                let col_start =
                    blend_mol_color(p1.ss.ribbon_color(), p2.ss.ribbon_color(), t0 as f32);
                let col_end =
                    blend_mol_color(p1.ss.ribbon_color(), p2.ss.ribbon_color(), t1 as f32);
                segments.push(RibbonSegment {
                    start: pos0,
                    end: pos1,
                    width_start,
                    width_end,
                    color_start: col_start,
                    color_end: col_end,
                    normal: norm0,
                });
            }
        }
        Self {
            segments,
            residue_count: n,
        }
    }
}
/// Crystal lattice parameters (lengths in Å, angles in degrees).
#[derive(Debug, Clone)]
pub struct LatticeParameters {
    /// Lattice parameter *a* (Å).
    pub a: f64,
    /// Lattice parameter *b* (Å).
    pub b: f64,
    /// Lattice parameter *c* (Å).
    pub c: f64,
    /// Angle α (°) between **b** and **c**.
    pub alpha: f64,
    /// Angle β (°) between **a** and **c**.
    pub beta: f64,
    /// Angle γ (°) between **a** and **b**.
    pub gamma: f64,
    /// Crystal system.
    pub system: CrystalSystem,
}
impl LatticeParameters {
    /// Construct cubic lattice parameters (a=b=c, all angles 90°).
    pub fn cubic(a: f64) -> Self {
        Self {
            a,
            b: a,
            c: a,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
            system: CrystalSystem::Cubic,
        }
    }
    /// Construct orthorhombic lattice parameters.
    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Self {
        Self {
            a,
            b,
            c,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
            system: CrystalSystem::Orthorhombic,
        }
    }
    /// Unit-cell volume in Å³.
    pub fn volume(&self) -> f64 {
        let al = self.alpha.to_radians();
        let be = self.beta.to_radians();
        let ga = self.gamma.to_radians();
        let ca = al.cos();
        let cb = be.cos();
        let cg = ga.cos();
        self.a * self.b * self.c * (1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg).sqrt()
    }
    /// Compute the three lattice vectors **a**, **b**, **c** in Cartesian coordinates.
    pub fn lattice_vectors(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let al = self.alpha.to_radians();
        let be = self.beta.to_radians();
        let ga = self.gamma.to_radians();
        let ax = self.a;
        let bx = self.b * ga.cos();
        let by = self.b * ga.sin();
        let cx = self.c * be.cos();
        let cy = self.c * (al.cos() - be.cos() * ga.cos()) / ga.sin();
        let cz = (self.c * self.c - cx * cx - cy * cy).max(0.0).sqrt();
        ([ax, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz])
    }
}
/// A complete molecule with atoms and bonds.
#[derive(Debug, Clone)]
pub struct Molecule {
    /// Molecule name or formula.
    pub name: String,
    /// All atoms.
    pub atoms: Vec<Atom>,
    /// All bonds.
    pub bonds: Vec<Bond>,
}
impl Molecule {
    /// Construct an empty molecule.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            atoms: Vec::new(),
            bonds: Vec::new(),
        }
    }
    /// Add an atom and return its assigned index.
    pub fn add_atom(&mut self, position: [f64; 3], element: Element) -> usize {
        let idx = self.atoms.len();
        self.atoms.push(Atom::new(idx, position, element));
        idx
    }
    /// Add a bond between two atom indices.
    pub fn add_bond(&mut self, a: usize, b: usize, order: BondOrder) {
        self.bonds.push(Bond {
            atom_a: a,
            atom_b: b,
            order,
        });
    }
    /// Auto-detect bonds using covalent radii with the default tolerance of 1.15.
    pub fn detect_bonds(&mut self) {
        self.bonds = detect_bonds(&self.atoms, 1.15);
    }
    /// Build a ball-and-stick scene for this molecule.
    pub fn ball_stick(&self, atom_scale: f64, bond_radius: f64) -> BallStickScene {
        BallStickScene::build(&self.atoms, &self.bonds, atom_scale, bond_radius)
    }
    /// Build a CPK scene for this molecule.
    pub fn cpk(&self, vdw_scale: f64) -> CpkScene {
        CpkScene::build(&self.atoms, vdw_scale)
    }
    /// Molecular formula as a sorted atom-count string (e.g. `"C6H6"`).
    pub fn formula(&self) -> String {
        let mut counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
        for a in &self.atoms {
            let sym = element_symbol(a.element);
            *counts.entry(sym.to_string()).or_insert(0) += 1;
        }
        let mut keys: Vec<&String> = counts.keys().collect();
        keys.sort();
        keys.iter()
            .map(|k| format!("{}{}", k, counts[*k]))
            .collect::<Vec<_>>()
            .join("")
    }
}
/// A molecular orbital (MO) stored as a real-valued wavefunction on a 3-D grid.
#[derive(Debug, Clone)]
pub struct MolecularOrbital {
    /// Human-readable label (e.g. `"HOMO"`, `"LUMO+1"`).
    pub label: String,
    /// Orbital energy in electron-volts (eV).
    pub energy_ev: f64,
    /// Orbital occupation (0.0–2.0 for closed-shell).
    pub occupation: f64,
    /// Symmetry label.
    pub symmetry: OrbitalSymmetry,
    /// Real-valued wavefunction grid (same layout as [`ElectronDensityGrid`]).
    pub wavefunction: ElectronDensityGrid,
}
impl MolecularOrbital {
    /// Construct a new MO with a zero-filled wavefunction grid.
    pub fn new(
        label: impl Into<String>,
        energy_ev: f64,
        occupation: f64,
        symmetry: OrbitalSymmetry,
        nx: usize,
        ny: usize,
        nz: usize,
        spacing: f64,
        origin: [f64; 3],
    ) -> Self {
        Self {
            label: label.into(),
            energy_ev,
            occupation,
            symmetry,
            wavefunction: ElectronDensityGrid::new(nx, ny, nz, spacing, origin),
        }
    }
    /// Compute electron density (ψ²) at the given grid point.
    pub fn density_at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        let psi = self.wavefunction.get(ix, iy, iz);
        psi * psi * self.occupation
    }
    /// Generate positive- and negative-phase isosurface seeds for a given |ψ| threshold.
    ///
    /// Returns `(positive_seeds, negative_seeds)`.
    pub fn phase_seeds(&self, threshold: f64) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let mut pos = Vec::new();
        let mut neg = Vec::new();
        for ix in 0..self.wavefunction.nx {
            for iy in 0..self.wavefunction.ny {
                for iz in 0..self.wavefunction.nz {
                    let v = self.wavefunction.get(ix, iy, iz);
                    let pt = self.wavefunction.position(ix, iy, iz);
                    if v >= threshold {
                        pos.push(pt);
                    } else if v <= -threshold {
                        neg.push(pt);
                    }
                }
            }
        }
        (pos, neg)
    }
}
/// A single frame in a molecular dynamics trajectory.
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Simulation time in picoseconds.
    pub time_ps: f64,
    /// Atomic positions (Å) indexed by atom.
    pub positions: Vec<[f64; 3]>,
    /// Optional per-atom velocities (Å/ps).
    pub velocities: Option<Vec<[f64; 3]>>,
    /// Optional total potential energy (kJ/mol).
    pub potential_energy: Option<f64>,
    /// Optional total kinetic energy (kJ/mol).
    pub kinetic_energy: Option<f64>,
    /// Optional instantaneous temperature (K).
    pub temperature: Option<f64>,
}
impl TrajectoryFrame {
    /// Compute the centre of geometry (unweighted mean position).
    pub fn center_of_geometry(&self) -> [f64; 3] {
        if self.positions.is_empty() {
            return [0.0; 3];
        }
        let n = self.positions.len() as f64;
        let sum = self
            .positions
            .iter()
            .fold([0.0f64; 3], |acc, p| add3(acc, *p));
        scale3(sum, 1.0 / n)
    }
    /// Compute the RMSD relative to a reference set of positions.
    pub fn rmsd(&self, reference: &[[f64; 3]]) -> f64 {
        if self.positions.len() != reference.len() || self.positions.is_empty() {
            return 0.0;
        }
        let n = self.positions.len() as f64;
        let sum_sq: f64 = self
            .positions
            .iter()
            .zip(reference.iter())
            .map(|(a, b)| {
                let d = sub3(*a, *b);
                dot3(d, d)
            })
            .sum();
        (sum_sq / n).sqrt()
    }
    /// Total mechanical energy (potential + kinetic), if both are available.
    pub fn total_energy(&self) -> Option<f64> {
        match (self.potential_energy, self.kinetic_energy) {
            (Some(pe), Some(ke)) => Some(pe + ke),
            _ => None,
        }
    }
}
/// Colour scheme selector for molecular rendering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColorScheme {
    /// Standard CPK element colours.
    Cpk,
    /// Residue index rainbow.
    ResidueRainbow,
    /// Secondary structure (helix/strand/loop).
    SecondaryStructure,
    /// B-factor / temperature factor.
    BFactor,
    /// Uniform single colour.
    Uniform,
}
/// A line edge of the crystal unit cell for wireframe display.
#[derive(Debug, Clone)]
pub struct CellEdge {
    /// Start vertex (Å, Cartesian).
    pub start: [f64; 3],
    /// End vertex (Å, Cartesian).
    pub end: [f64; 3],
    /// Display colour.
    pub color: MolColor,
}
/// An atom described by fractional coordinates in the unit cell.
#[derive(Debug, Clone)]
pub struct AtomFractional {
    /// Fractional coordinate along **a**.
    pub fa: f64,
    /// Fractional coordinate along **b**.
    pub fb: f64,
    /// Fractional coordinate along **c**.
    pub fc: f64,
    /// Element type.
    pub element: Element,
    /// Occupancy (0–1).
    pub occupancy: f64,
}
