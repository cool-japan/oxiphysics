//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{
    KB, cg_angle_energy, cg_bond_energy, interpolate_table, martini_lj_energy, martini_lj_force,
};

/// Iterative Boltzmann Inversion (IBI) potential optimiser.
///
/// Differs from [`IbiPotential`] in that it stores the target RDF and the
/// current tabulated potential separately as plain vectors, and uses
/// temperature (K) directly rather than pre-computed kBT.
#[derive(Debug, Clone)]
pub struct IterativeBoltzmannInversion {
    /// Target radial distribution function (same length as `r_bins`).
    pub target_rdf: Vec<f64>,
    /// Current CG pair potential (kJ/mol).
    pub current_potential: Vec<f64>,
    /// Bin centres for the pair distance (nm).
    pub r_bins: Vec<f64>,
}
impl IterativeBoltzmannInversion {
    /// Create from a target RDF and bin centres.
    ///
    /// The initial potential is U_0(r) = -k_B T ln(g_target(r)).
    pub fn from_target(target_rdf: Vec<f64>, r_bins: Vec<f64>, temperature: f64) -> Self {
        let kbt = KB * temperature / 1.602_176_634e-22;
        let n_a = 6.022_140_76e23_f64;
        let kbt_kj = KB * temperature * n_a / 1000.0;
        let _ = kbt;
        let current_potential: Vec<f64> = target_rdf
            .iter()
            .map(|&g| {
                if g > 1e-10 {
                    -kbt_kj * g.ln()
                } else {
                    kbt_kj * 20.0
                }
            })
            .collect();
        Self {
            target_rdf,
            current_potential,
            r_bins,
        }
    }
    /// Perform one IBI update step.
    ///
    /// U_{n+1}(r) = U_n(r) + k_B T ln(g_current(r) / g_target(r))
    pub fn update_potential(&mut self, current_rdf: &[f64], temperature: f64) {
        let n_a = 6.022_140_76e23_f64;
        let kbt_kj = KB * temperature * n_a / 1000.0;
        for (i, (u, (&g_curr, &g_targ))) in self
            .current_potential
            .iter_mut()
            .zip(current_rdf.iter().zip(self.target_rdf.iter()))
            .enumerate()
        {
            let _ = i;
            if g_curr > 1e-10 && g_targ > 1e-10 {
                *u += kbt_kj * (g_curr / g_targ).ln();
            }
        }
    }
}
/// Simplified MARTINI bead interaction types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MartiniType {
    /// Apolar C1 bead.
    C1,
    /// Apolar C2 bead.
    C2,
    /// Non-polar Na bead.
    Na,
    /// Charged Qa bead (negative).
    Qa,
    /// Charged Qd bead (positive).
    Qd,
    /// Polar P1 bead.
    P1,
    /// Polar P2 bead.
    P2,
}
/// Radial distribution function (RDF) calculation for CG beads.
///
/// Computes g(r) from a set of bead positions in a cubic box.
pub struct RdfCalculator {
    /// Maximum distance for RDF calculation (nm).
    pub r_max: f64,
    /// Number of bins.
    pub n_bins: usize,
    /// Bin width (nm).
    pub dr: f64,
    /// Accumulated histogram counts.
    pub(super) histogram: Vec<f64>,
    /// Number of frames accumulated.
    pub(super) n_frames: usize,
    /// Number of particles in each frame (assumed constant).
    pub(super) n_particles: usize,
    /// Box side length (nm).
    pub(super) box_length: f64,
}
impl RdfCalculator {
    /// Create a new RDF calculator.
    pub fn new(r_max: f64, n_bins: usize, n_particles: usize, box_length: f64) -> Self {
        let dr = r_max / n_bins as f64;
        Self {
            r_max,
            n_bins,
            dr,
            histogram: vec![0.0; n_bins],
            n_frames: 0,
            n_particles,
            box_length,
        }
    }
    /// Accumulate pair distances from a single frame of positions.
    ///
    /// Applies minimum image convention for periodic boundaries.
    pub fn accumulate(&mut self, positions: &[[f64; 3]]) {
        let n = positions.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let mut dr2 = 0.0;
                for (pi, qi) in positions[i].iter().zip(positions[j].iter()) {
                    let mut dx = qi - pi;
                    dx -= self.box_length * (dx / self.box_length).round();
                    dr2 += dx * dx;
                }
                let r = dr2.sqrt();
                if r < self.r_max {
                    let bin = (r / self.dr) as usize;
                    if bin < self.n_bins {
                        self.histogram[bin] += 2.0;
                    }
                }
            }
        }
        self.n_frames += 1;
    }
    /// Compute the normalized RDF g(r).
    ///
    /// g(r) = histogram(r) / (N_frames * N * 4πr²Δr * ρ)
    pub fn compute_rdf(&self) -> (Vec<f64>, Vec<f64>) {
        let volume = self.box_length * self.box_length * self.box_length;
        let density = self.n_particles as f64 / volume;
        let n = self.n_particles as f64;
        let mut r_values = Vec::with_capacity(self.n_bins);
        let mut g_values = Vec::with_capacity(self.n_bins);
        for i in 0..self.n_bins {
            let r = (i as f64 + 0.5) * self.dr;
            r_values.push(r);
            let r_inner = i as f64 * self.dr;
            let r_outer = (i as f64 + 1.0) * self.dr;
            let shell_volume =
                (4.0 / 3.0) * std::f64::consts::PI * (r_outer.powi(3) - r_inner.powi(3));
            let ideal_count = n * density * shell_volume;
            let g = if self.n_frames > 0 && ideal_count > 0.0 {
                self.histogram[i] / (self.n_frames as f64 * ideal_count)
            } else {
                0.0
            };
            g_values.push(g);
        }
        (r_values, g_values)
    }
    /// Reset accumulated data.
    pub fn reset(&mut self) {
        self.histogram.fill(0.0);
        self.n_frames = 0;
    }
}
/// A harmonic bond between two coarse-grained beads.
#[derive(Debug, Clone)]
pub struct CgBond {
    /// Index of first bead.
    pub bead_i: usize,
    /// Index of second bead.
    pub bead_j: usize,
    /// Equilibrium bond length (nm).
    pub r0: f64,
    /// Spring constant (kJ/mol/nm^2).
    pub k: f64,
}
/// A coarse-grained bead with position, mass, type, and charge.
#[derive(Debug, Clone)]
pub struct CgBead {
    /// Position in 3D space (nm).
    pub position: [f64; 3],
    /// Mass of the bead (amu).
    pub mass: f64,
    /// Numeric bead type identifier.
    pub bead_type: u8,
    /// Partial charge on the bead (elementary charges).
    pub charge: f64,
}
/// Region classification for adaptive resolution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResolutionRegion {
    /// Fully atomistic region.
    Atomistic,
    /// Hybrid transition region.
    Hybrid,
    /// Fully coarse-grained region.
    CoarseGrained,
}
/// Force matching (multiscale coarse-graining) result for a pair interaction.
///
/// Given atomistic reference forces and CG coordinates, determines the
/// optimal CG pair force function by least-squares fitting.
pub struct ForceMatchResult {
    /// Bin edges for the pair distance histogram (nm).
    pub bin_edges: Vec<f64>,
    /// Mean force in each distance bin (kJ/mol/nm).
    pub mean_forces: Vec<f64>,
    /// Number of samples in each bin.
    pub counts: Vec<usize>,
}
/// LJ interaction parameters between two MARTINI bead types.
#[derive(Debug, Clone, Copy)]
pub struct MartiniLjParams {
    /// LJ ε (kJ mol^-1).
    pub epsilon: f64,
    /// LJ σ (nm).
    pub sigma: f64,
}
impl MartiniLjParams {
    /// Look up MARTINI 2.x ε and σ between two bead types.
    ///
    /// Returns standard interaction level based on bead hydrophilicity.
    pub fn lookup(bead_a: MartiniBead, bead_b: MartiniBead) -> Self {
        use MartiniBead::*;
        let (a, b) = if (bead_a as u8) <= (bead_b as u8) {
            (bead_a, bead_b)
        } else {
            (bead_b, bead_a)
        };
        let level = match (a, b) {
            (P4, P4) => 5,
            (Nda, P4) | (P4, Nda) => 4,
            (C1, P4) | (P4, C1) => 1,
            (C1, C1) => 5,
            (C3, C3) => 4,
            (Qa, Qd) | (Qd, Qa) => 5,
            (C1, Qd) | (Qd, C1) | (C1, Qa) | (Qa, C1) => 1,
            _ => 3,
        };
        let (epsilon, sigma) = match level {
            5 => (5.6, 0.47),
            4 => (5.0, 0.47),
            3 => (4.5, 0.47),
            2 => (4.0, 0.47),
            1 => (3.1, 0.47),
            _ => (3.5, 0.47),
        };
        MartiniLjParams { epsilon, sigma }
    }
}
/// MARTINI-inspired coarse-graining scheme: N heavy atoms → 1 bead.
#[derive(Debug, Clone)]
pub struct MartiniMappingScheme {
    /// Number of heavy atoms per bead.
    pub n_heavy: usize,
    /// Mass of the resulting CG bead (amu).
    pub bead_mass: f64,
}
impl MartiniMappingScheme {
    /// Create a new `MartiniMappingScheme`.
    pub fn new(n_heavy: usize, bead_mass: f64) -> Self {
        Self { n_heavy, bead_mass }
    }
    /// Map atomistic positions to a single CG bead position using COM.
    ///
    /// Returns the center-of-mass position weighted by `masses`.
    pub fn map_atomistic(&self, positions: &[[f64; 3]], masses: &[f64]) -> [f64; 3] {
        let mut com = [0.0f64; 3];
        let mut total_m = 0.0f64;
        for (pos, &m) in positions.iter().zip(masses.iter()) {
            for d in 0..3 {
                com[d] += m * pos[d];
            }
            total_m += m;
        }
        if total_m > 0.0 {
            for v in &mut com {
                *v /= total_m;
            }
        }
        com
    }
}
/// A harmonic angle defined by three coarse-grained beads.
#[derive(Debug, Clone)]
pub struct CgAngle {
    /// Index of first bead.
    pub i: usize,
    /// Index of central bead.
    pub j: usize,
    /// Index of third bead.
    pub k: usize,
    /// Equilibrium angle (radians).
    pub theta0: f64,
    /// Spring constant (kJ/mol/rad^2).
    pub k_theta: f64,
}
/// Mapping from all-atom representation to coarse-grained beads.
///
/// Each bead corresponds to a group of atoms; positions and velocities
/// are computed as center-of-mass quantities.
pub struct CgMapping {
    /// Indices of atoms belonging to each bead.
    pub atom_indices: Vec<Vec<usize>>,
    /// Total mass of each bead (sum of constituent atom masses).
    pub bead_masses: Vec<f64>,
}
impl CgMapping {
    /// Create an empty mapping.
    pub fn new() -> Self {
        CgMapping {
            atom_indices: Vec::new(),
            bead_masses: Vec::new(),
        }
    }
    /// Add a bead defined by the given atom indices and total mass.
    pub fn add_bead(&mut self, atom_indices: Vec<usize>, total_mass: f64) {
        self.atom_indices.push(atom_indices);
        self.bead_masses.push(total_mass);
    }
    /// Compute center-of-mass positions for all beads.
    ///
    /// # Arguments
    /// * `atom_positions` - slice of per-atom positions.
    /// * `atom_masses`    - slice of per-atom masses.
    pub fn map_positions(&self, atom_positions: &[[f64; 3]], atom_masses: &[f64]) -> Vec<[f64; 3]> {
        self.atom_indices
            .iter()
            .map(|indices| {
                let mut com = [0.0f64; 3];
                let mut total_mass = 0.0f64;
                for &idx in indices {
                    let m = atom_masses[idx];
                    total_mass += m;
                    for d in 0..3 {
                        com[d] += m * atom_positions[idx][d];
                    }
                }
                if total_mass > 0.0 {
                    for v in &mut com {
                        *v /= total_mass;
                    }
                }
                com
            })
            .collect()
    }
    /// Compute center-of-mass velocities for all beads.
    ///
    /// # Arguments
    /// * `atom_velocities` - slice of per-atom velocities.
    /// * `atom_masses`     - slice of per-atom masses.
    pub fn map_velocities(
        &self,
        atom_velocities: &[[f64; 3]],
        atom_masses: &[f64],
    ) -> Vec<[f64; 3]> {
        self.atom_indices
            .iter()
            .map(|indices| {
                let mut com_vel = [0.0f64; 3];
                let mut total_mass = 0.0f64;
                for &idx in indices {
                    let m = atom_masses[idx];
                    total_mass += m;
                    for d in 0..3 {
                        com_vel[d] += m * atom_velocities[idx][d];
                    }
                }
                if total_mass > 0.0 {
                    for v in &mut com_vel {
                        *v /= total_mass;
                    }
                }
                com_vel
            })
            .collect()
    }
    /// Compute center-of-mass forces for all beads.
    ///
    /// The CG force is simply the sum of all atomic forces in each bead group.
    pub fn map_forces(&self, atom_forces: &[[f64; 3]]) -> Vec<[f64; 3]> {
        self.atom_indices
            .iter()
            .map(|indices| {
                let mut com_force = [0.0f64; 3];
                for &idx in indices {
                    for d in 0..3 {
                        com_force[d] += atom_forces[idx][d];
                    }
                }
                com_force
            })
            .collect()
    }
    /// Return the number of coarse-grained beads.
    pub fn n_beads(&self) -> usize {
        self.atom_indices.len()
    }
}
/// A system of coarse-grained beads in a periodic cubic box.
#[derive(Debug, Clone)]
pub struct CgSystem {
    /// All beads in the system.
    pub beads: Vec<CgBead>,
    /// Box side lengths \[Lx, Ly, Lz\] (nm).
    pub box_size: [f64; 3],
}
impl CgSystem {
    /// Create an empty `CgSystem` with the given box dimensions.
    pub fn new(box_size: [f64; 3]) -> Self {
        Self {
            beads: Vec::new(),
            box_size,
        }
    }
    /// Add a `CgBead` to the system.
    pub fn add_bead(&mut self, bead: CgBead) {
        self.beads.push(bead);
    }
    /// Number of beads in the system.
    pub fn n_beads(&self) -> usize {
        self.beads.len()
    }
    /// Total mass of all beads (amu).
    pub fn total_mass(&self) -> f64 {
        self.beads.iter().map(|b| b.mass).sum()
    }
    /// Center of mass of the system.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let mut com = [0.0f64; 3];
        let mut total_m = 0.0f64;
        for bead in &self.beads {
            for (d, c) in com.iter_mut().enumerate() {
                *c += bead.mass * bead.position[d];
            }
            total_m += bead.mass;
        }
        if total_m > 0.0 {
            for v in &mut com {
                *v /= total_m;
            }
        }
        com
    }
}
/// Lennard-Jones non-bonded potential for coarse-grained beads.
#[derive(Debug, Clone)]
pub struct CgNonbonded {
    /// LJ epsilon (kJ/mol).
    pub epsilon: f64,
    /// LJ sigma (nm).
    pub sigma: f64,
    /// Interaction cutoff (nm).
    pub cutoff: f64,
}
impl CgNonbonded {
    /// Create a new `CgNonbonded` interaction.
    pub fn new(epsilon: f64, sigma: f64, cutoff: f64) -> Self {
        Self {
            epsilon,
            sigma,
            cutoff,
        }
    }
    /// LJ energy at distance `r` (returns 0 beyond cutoff).
    pub fn energy(&self, r: f64) -> f64 {
        if r >= self.cutoff || r < 1e-12 {
            return 0.0;
        }
        martini_lj_energy(r, self.sigma, self.epsilon)
    }
    /// Magnitude of the LJ force at distance `r` (positive = repulsive).
    /// Returns 0 beyond cutoff.
    pub fn force(&self, r: f64) -> f64 {
        if r >= self.cutoff || r < 1e-12 {
            return 0.0;
        }
        martini_lj_force(r, self.sigma, self.epsilon)
    }
}
/// MARTINI-like coarse-grained Lennard-Jones potential parameters.
pub struct MartiniPotential {
    /// Sigma (distance) parameter (nm).
    pub sigma: f64,
    /// Epsilon (energy) parameter (kJ/mol).
    pub epsilon: f64,
    /// Bead interaction type.
    pub interaction_type: MartiniType,
}
/// Elastic Network Model built from Cα positions.
///
/// Reference: Tirion (1996), Bahar *et al.* (1997).
#[derive(Debug, Clone)]
pub struct ElasticNetworkModel {
    /// ENM contacts.
    pub contacts: Vec<EnmContact>,
    /// Number of nodes (Cα atoms).
    pub n_nodes: usize,
    /// Cutoff radius for contact formation (nm).
    pub cutoff: f64,
    /// Uniform spring constant γ (kJ mol^-1 nm^-2).
    pub gamma: f64,
}
impl ElasticNetworkModel {
    /// Build an ENM from Cα positions with a given cutoff and spring constant.
    pub fn build(positions: &[[f64; 3]], cutoff: f64, gamma: f64) -> Self {
        let n = positions.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r0 = (dx * dx + dy * dy + dz * dz).sqrt();
                if r0 < cutoff {
                    contacts.push(EnmContact { i, j, r0, gamma });
                }
            }
        }
        ElasticNetworkModel {
            contacts,
            n_nodes: n,
            cutoff,
            gamma,
        }
    }
    /// Compute the elastic energy given current positions.
    ///
    /// E = (γ/2) Σ_{contacts} (r_{ij} - r0_{ij})²
    pub fn energy(&self, positions: &[[f64; 3]]) -> f64 {
        self.contacts
            .iter()
            .map(|c| {
                let dx = positions[c.j][0] - positions[c.i][0];
                let dy = positions[c.j][1] - positions[c.i][1];
                let dz = positions[c.j][2] - positions[c.i][2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                let dr = r - c.r0;
                0.5 * c.gamma * dr * dr
            })
            .sum()
    }
    /// Compute ENM forces on each node.
    pub fn forces(&self, positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut f = vec![[0.0f64; 3]; self.n_nodes];
        for c in &self.contacts {
            let dx = positions[c.j][0] - positions[c.i][0];
            let dy = positions[c.j][1] - positions[c.i][1];
            let dz = positions[c.j][2] - positions[c.i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < 1e-30 {
                continue;
            }
            let dr = r - c.r0;
            let factor = c.gamma * dr / r;
            let fx = factor * dx;
            let fy = factor * dy;
            let fz = factor * dz;
            f[c.j][0] -= fx;
            f[c.j][1] -= fy;
            f[c.j][2] -= fz;
            f[c.i][0] += fx;
            f[c.i][1] += fy;
            f[c.i][2] += fz;
        }
        f
    }
    /// Build the 3N×3N Hessian matrix of the ENM.
    ///
    /// H_{ia,jb} = -γ (r_{ij,a} r_{ij,b}) / r²  for i≠j contacts,
    /// diagonal blocks are sums of off-diagonal elements.
    pub fn hessian(&self, positions: &[[f64; 3]]) -> Vec<f64> {
        let dim = 3 * self.n_nodes;
        let mut h = vec![0.0f64; dim * dim];
        for c in &self.contacts {
            let dx = [
                positions[c.j][0] - positions[c.i][0],
                positions[c.j][1] - positions[c.i][1],
                positions[c.j][2] - positions[c.i][2],
            ];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 < 1e-60 {
                continue;
            }
            for a in 0..3 {
                for b in 0..3 {
                    let val = -c.gamma * dx[a] * dx[b] / r2;
                    let row_ij = 3 * c.i + a;
                    let col_ij = 3 * c.j + b;
                    h[row_ij * dim + col_ij] += val;
                    h[col_ij * dim + row_ij] += val;
                    h[(3 * c.i + a) * dim + (3 * c.i + b)] -= val;
                    h[(3 * c.j + a) * dim + (3 * c.j + b)] -= val;
                }
            }
        }
        h
    }
    /// Mean square fluctuation of each node from ENM theory.
    ///
    /// ⟨Δr_i²⟩ = k_B T * Σ_{k≠0,1,2,3,4,5} (C_ik)² / λ_k
    ///
    /// This is a diagonal approximation using the pseudo-inverse of H.
    /// Here we use a simple power-iteration approximation: we return the
    /// diagonal of H⁺ × (k_B T) scaled by the harmonic spring constant.
    ///
    /// For production, use a full eigensolver; here we provide the contact count
    /// as a proxy for flexibility (more contacts → stiffer → smaller fluctuations).
    pub fn contact_count(&self) -> Vec<usize> {
        let mut counts = vec![0usize; self.n_nodes];
        for c in &self.contacts {
            counts[c.i] += 1;
            counts[c.j] += 1;
        }
        counts
    }
}
/// MARTINI 2.x CG bead type identifiers.
///
/// Reference: Marrink *et al.* (2007); de Jong *et al.* (2013).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MartiniBead {
    /// Polar (water-like).
    P4,
    /// Intermediate polar.
    Nda,
    /// Apolar (hydrophobic).
    C1,
    /// Slightly apolar.
    C3,
    /// Apolar, small ring bead.
    Sc3,
    /// Charged (positive).
    Qd,
    /// Charged (negative).
    Qa,
}
/// Tabulated pair potential for CG interactions.
///
/// Stores U(r) and F(r) on a uniform grid for efficient lookup
/// via linear interpolation.
pub struct TabulatedPotential {
    /// Minimum distance (nm).
    pub r_min: f64,
    /// Maximum distance (nm).
    pub r_max: f64,
    /// Number of table entries.
    pub n_points: usize,
    /// Grid spacing (nm).
    pub dr: f64,
    /// Potential values U(r) (kJ/mol).
    pub energies: Vec<f64>,
    /// Force values F(r) = -dU/dr (kJ/mol/nm).
    pub forces: Vec<f64>,
}
impl TabulatedPotential {
    /// Create a tabulated potential from an analytic potential function.
    ///
    /// The force is computed via central finite differences.
    pub fn from_function<F>(r_min: f64, r_max: f64, n_points: usize, potential_fn: F) -> Self
    where
        F: Fn(f64) -> f64,
    {
        let dr = (r_max - r_min) / (n_points - 1) as f64;
        let mut energies = Vec::with_capacity(n_points);
        let mut forces = Vec::with_capacity(n_points);
        for i in 0..n_points {
            let r = r_min + i as f64 * dr;
            energies.push(potential_fn(r));
        }
        for i in 0..n_points {
            let f = if i == 0 {
                -(energies[1] - energies[0]) / dr
            } else if i == n_points - 1 {
                -(energies[n_points - 1] - energies[n_points - 2]) / dr
            } else {
                -(energies[i + 1] - energies[i - 1]) / (2.0 * dr)
            };
            forces.push(f);
        }
        Self {
            r_min,
            r_max,
            n_points,
            dr,
            energies,
            forces,
        }
    }
    /// Create a tabulated LJ potential.
    pub fn lj(r_min: f64, r_max: f64, n_points: usize, sigma: f64, epsilon: f64) -> Self {
        Self::from_function(r_min, r_max, n_points, |r| {
            martini_lj_energy(r, sigma, epsilon)
        })
    }
    /// Look up energy at distance r via linear interpolation.
    pub fn energy(&self, r: f64) -> f64 {
        self.interpolate(r, &self.energies)
    }
    /// Look up force at distance r via linear interpolation.
    pub fn force(&self, r: f64) -> f64 {
        self.interpolate(r, &self.forces)
    }
    /// Linear interpolation helper.
    fn interpolate(&self, r: f64, table: &[f64]) -> f64 {
        if r <= self.r_min {
            return table[0];
        }
        if r >= self.r_max {
            return *table.last().unwrap_or(&0.0);
        }
        let idx_f = (r - self.r_min) / self.dr;
        let idx = idx_f as usize;
        let frac = idx_f - idx as f64;
        if idx + 1 < table.len() {
            table[idx] * (1.0 - frac) + table[idx + 1] * frac
        } else {
            table[idx]
        }
    }
    /// Shift the potential so that U(r_max) = 0.
    pub fn shift_to_zero_at_cutoff(&mut self) {
        let u_cut = *self.energies.last().unwrap_or(&0.0);
        for e in &mut self.energies {
            *e -= u_cut;
        }
    }
}
/// IBI run accumulator: tracks convergence over iterations.
#[derive(Debug, Clone)]
pub struct IbiRunner {
    /// Internal IBI state.
    pub ibi: IterativeBoltzmannInversion,
    /// History of max |ΔU| per iteration (kJ mol^-1).
    pub convergence_history: Vec<f64>,
}
impl IbiRunner {
    /// Create from a target RDF.
    pub fn new(target_rdf: Vec<f64>, r_bins: Vec<f64>, temperature: f64) -> Self {
        IbiRunner {
            ibi: IterativeBoltzmannInversion::from_target(target_rdf, r_bins, temperature),
            convergence_history: Vec::new(),
        }
    }
    /// Perform one IBI iteration and record the max potential change.
    pub fn iterate(&mut self, current_rdf: &[f64], temperature: f64) {
        let old_potential = self.ibi.current_potential.clone();
        self.ibi.update_potential(current_rdf, temperature);
        let max_delta = old_potential
            .iter()
            .zip(self.ibi.current_potential.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        self.convergence_history.push(max_delta);
    }
    /// Check if the IBI has converged to within `tol` kJ mol^-1.
    pub fn converged(&self, tol: f64) -> bool {
        self.convergence_history
            .last()
            .map(|&d| d < tol)
            .unwrap_or(false)
    }
}
/// Least-squares force matching result for a CG potential table.
///
/// FM minimises Σ_{frames} Σ_{atoms} |F^AA_i - F^CG_i|² over CG potential
/// parameters.  Here we implement the accumulation side: given AA forces and
/// bead positions, build the normal equations A^T A and A^T b.
#[derive(Debug, Clone)]
pub struct FmNormalEquations {
    /// Number of bins in the CG pair potential.
    pub n_bins: usize,
    /// Bin width dr (nm).
    pub dr: f64,
    /// Minimum distance r_min (nm).
    pub r_min: f64,
    /// Gram matrix A^T A (n_bins × n_bins, flattened).
    pub ata: Vec<f64>,
    /// Right-hand side A^T b (length n_bins).
    pub atb: Vec<f64>,
    /// Number of frames accumulated.
    pub n_frames: usize,
}
impl FmNormalEquations {
    /// Create empty normal equations for a pair potential with `n_bins` bins.
    pub fn new(n_bins: usize, r_min: f64, dr: f64) -> Self {
        FmNormalEquations {
            n_bins,
            dr,
            r_min,
            ata: vec![0.0f64; n_bins * n_bins],
            atb: vec![0.0f64; n_bins],
            n_frames: 0,
        }
    }
    /// Accumulate one frame of force-matching data.
    ///
    /// For each pair (i,j) in `bead_pairs`, the piecewise-linear basis function
    /// contribution to the normal equations is added.  `target_forces` is the
    /// projected reference force along the inter-bead axis.
    pub fn accumulate(
        &mut self,
        bead_positions: &[[f64; 3]],
        bead_pairs: &[(usize, usize)],
        target_forces: &[f64],
    ) {
        for (&(i, j), &f_target) in bead_pairs.iter().zip(target_forces.iter()) {
            let dx = bead_positions[j][0] - bead_positions[i][0];
            let dy = bead_positions[j][1] - bead_positions[i][1];
            let dz = bead_positions[j][2] - bead_positions[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < self.r_min || r >= self.r_min + self.n_bins as f64 * self.dr {
                continue;
            }
            let idx_f = (r - self.r_min) / self.dr;
            let k = idx_f as usize;
            let t = idx_f - k as f64;
            if k + 1 >= self.n_bins {
                continue;
            }
            let basis_k = (1.0 - t) / self.dr;
            let basis_kp1 = t / self.dr;
            self.ata[k * self.n_bins + k] += basis_k * basis_k;
            self.ata[k * self.n_bins + (k + 1)] += basis_k * basis_kp1;
            self.ata[(k + 1) * self.n_bins + k] += basis_kp1 * basis_k;
            self.ata[(k + 1) * self.n_bins + (k + 1)] += basis_kp1 * basis_kp1;
            self.atb[k] += basis_k * f_target;
            self.atb[k + 1] += basis_kp1 * f_target;
        }
        self.n_frames += 1;
    }
}
/// Properties for a CG bead type.
#[derive(Debug, Clone)]
pub struct CgBeadProperties {
    /// Bead type.
    pub bead_type: CgBeadType,
    /// Mass (amu).
    pub mass: f64,
    /// Charge (elementary charges).
    pub charge: f64,
    /// LJ sigma parameter (nm).
    pub sigma: f64,
    /// LJ epsilon parameter (kJ/mol).
    pub epsilon: f64,
}
impl CgBeadProperties {
    /// Create a new CG bead type with specified properties.
    pub fn new(bead_type: CgBeadType, mass: f64, charge: f64, sigma: f64, epsilon: f64) -> Self {
        Self {
            bead_type,
            mass,
            charge,
            sigma,
            epsilon,
        }
    }
    /// Create a MARTINI-like backbone bead.
    pub fn martini_backbone() -> Self {
        Self::new(CgBeadType::Backbone, 72.0, 0.0, 0.47, 3.5)
    }
    /// Create a MARTINI-like water bead (4:1 mapping).
    pub fn martini_water() -> Self {
        Self::new(CgBeadType::Water, 72.0, 0.0, 0.47, 3.375)
    }
}
/// Adaptive resolution scheme for multi-scale coupling.
///
/// Provides a smooth transition between atomistic (AA) and coarse-grained (CG)
/// regions using a switching function w(x) ∈ \[0, 1\].
///
/// w = 1 in the AA region, w = 0 in the CG region, and smoothly varies
/// in the hybrid (transition) region.
pub struct AdaptiveResolution {
    /// Center of the AA region (nm).
    pub center: [f64; 3],
    /// Inner radius of the AA region (nm).
    pub r_aa: f64,
    /// Outer radius (start of pure CG region) (nm).
    pub r_cg: f64,
}
impl AdaptiveResolution {
    /// Create a new adaptive resolution scheme.
    ///
    /// The AA region extends from 0 to `r_aa`, the hybrid region from
    /// `r_aa` to `r_cg`, and pure CG beyond `r_cg`.
    pub fn new(center: [f64; 3], r_aa: f64, r_cg: f64) -> Self {
        assert!(r_cg > r_aa, "r_cg must be greater than r_aa");
        Self { center, r_aa, r_cg }
    }
    /// Compute the distance from a point to the center.
    fn distance_from_center(&self, position: &[f64; 3]) -> f64 {
        let dx = position[0] - self.center[0];
        let dy = position[1] - self.center[1];
        let dz = position[2] - self.center[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
    /// Switching function w(r).
    ///
    /// Uses a cosine switching function for smooth interpolation:
    /// w = 1                                   if r <= r_aa
    /// w = 0.5 * (1 + cos(π*(r-r_aa)/(r_cg-r_aa)))  if r_aa < r < r_cg
    /// w = 0                                   if r >= r_cg
    pub fn weight(&self, position: &[f64; 3]) -> f64 {
        let r = self.distance_from_center(position);
        if r <= self.r_aa {
            1.0
        } else if r >= self.r_cg {
            0.0
        } else {
            let s = (r - self.r_aa) / (self.r_cg - self.r_aa);
            0.5 * (1.0 + (std::f64::consts::PI * s).cos())
        }
    }
    /// Blend atomistic and CG forces using the switching function.
    ///
    /// F_hybrid = w * F_AA + (1 - w) * F_CG
    pub fn blend_forces(
        &self,
        position: &[f64; 3],
        force_aa: &[f64; 3],
        force_cg: &[f64; 3],
    ) -> [f64; 3] {
        let w = self.weight(position);
        let mut result = [0.0f64; 3];
        for d in 0..3 {
            result[d] = w * force_aa[d] + (1.0 - w) * force_cg[d];
        }
        result
    }
    /// Classify a position as AA, hybrid, or CG.
    pub fn region(&self, position: &[f64; 3]) -> ResolutionRegion {
        let r = self.distance_from_center(position);
        if r <= self.r_aa {
            ResolutionRegion::Atomistic
        } else if r >= self.r_cg {
            ResolutionRegion::CoarseGrained
        } else {
            ResolutionRegion::Hybrid
        }
    }
}
/// CG pair potential with gradient (force) computed by finite difference.
#[derive(Debug, Clone)]
pub struct CgPairPotential {
    /// Minimum distance r_min (nm).
    pub r_min: f64,
    /// Bin width dr (nm).
    pub dr: f64,
    /// Tabulated potential values U(r) (kJ mol^-1).
    pub u_table: Vec<f64>,
}
impl CgPairPotential {
    /// Create from tabulated values.
    pub fn new(r_min: f64, dr: f64, u_table: Vec<f64>) -> Self {
        CgPairPotential { r_min, dr, u_table }
    }
    /// Evaluate U(r) by linear interpolation.
    pub fn energy(&self, r: f64) -> f64 {
        interpolate_table(r, self.r_min, self.dr, &self.u_table)
    }
    /// Evaluate -dU/dr by central finite difference on the table.
    pub fn force_magnitude(&self, r: f64) -> f64 {
        let u_plus = self.energy(r + self.dr * 0.5);
        let u_minus = self.energy(r - self.dr * 0.5);
        -(u_plus - u_minus) / self.dr
    }
    /// Build a Lennard-Jones pair potential table.
    pub fn lj_table(epsilon: f64, sigma: f64, r_min: f64, r_max: f64, n_bins: usize) -> Self {
        let dr = (r_max - r_min) / n_bins as f64;
        let u_table = (0..n_bins)
            .map(|i| {
                let r = r_min + (i as f64 + 0.5) * dr;
                let sr6 = (sigma / r).powi(6);
                4.0 * epsilon * (sr6 * sr6 - sr6)
            })
            .collect();
        CgPairPotential { r_min, dr, u_table }
    }
}
/// Reverse coarse-graining (back-mapping) configuration.
///
/// Maps CG bead positions back to all-atom (AA) positions using
/// local reference conformations stored per bead.
#[derive(Debug, Clone)]
pub struct BackMapper {
    /// Number of beads.
    pub n_beads: usize,
    /// Reference AA positions in the bead local frame, per bead.
    /// Indexed \[bead\]\[atom_in_bead\]\[xyz\].
    pub reference_frames: Vec<Vec<[f64; 3]>>,
    /// Atom masses per bead, for COM calculation.
    pub atom_masses: Vec<Vec<f64>>,
}
impl BackMapper {
    /// Create a back-mapper from reference frames.
    pub fn new(reference_frames: Vec<Vec<[f64; 3]>>, atom_masses: Vec<Vec<f64>>) -> Self {
        let n = reference_frames.len();
        BackMapper {
            n_beads: n,
            reference_frames,
            atom_masses,
        }
    }
    /// Back-map a CG configuration to AA positions.
    ///
    /// For each bead, places the reference AA atoms around the CG bead centre.
    /// No rotational fitting is performed (identity orientation assumed).
    pub fn back_map(&self, cg_positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
        let mut aa_positions = Vec::new();
        for (bead_idx, bead_pos) in cg_positions.iter().enumerate() {
            for ref_pos in &self.reference_frames[bead_idx] {
                aa_positions.push([
                    bead_pos[0] + ref_pos[0],
                    bead_pos[1] + ref_pos[1],
                    bead_pos[2] + ref_pos[2],
                ]);
            }
        }
        aa_positions
    }
    /// Number of all-atom particles after back-mapping.
    pub fn n_atoms(&self) -> usize {
        self.reference_frames.iter().map(|f| f.len()).sum()
    }
    /// Centre of mass of a single bead's reference frame.
    pub fn bead_com(&self, bead_idx: usize) -> [f64; 3] {
        let frame = &self.reference_frames[bead_idx];
        let masses = &self.atom_masses[bead_idx];
        let total_mass: f64 = masses.iter().sum();
        let mut com = [0.0f64; 3];
        for (pos, &m) in frame.iter().zip(masses.iter()) {
            for d in 0..3 {
                com[d] += m * pos[d];
            }
        }
        if total_mass > 0.0 {
            for v in &mut com {
                *v /= total_mass;
            }
        }
        com
    }
}
/// A coarse-grained bead particle.
pub struct Bead {
    /// Position in 3D space (nm).
    pub position: [f64; 3],
    /// Velocity in 3D space (nm/ps).
    pub velocity: [f64; 3],
    /// Mass of the bead (amu).
    pub mass: f64,
    /// Bead type identifier.
    pub bead_type: u32,
}
/// Enumeration of common coarse-grained bead types for biomolecular simulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum CgBeadType {
    /// Backbone bead (e.g., protein backbone).
    Backbone,
    /// Sidechain bead.
    Sidechain,
    /// Polar headgroup (e.g., lipid).
    PolarHead,
    /// Apolar tail (e.g., lipid tail).
    ApolarTail,
    /// Water bead (typically 4 water molecules).
    Water,
    /// Ion bead.
    Ion,
    /// Custom bead type with an identifier.
    Custom(u32),
}
/// Composed CG force field: bonds + angles + nonbonded tabulated.
#[derive(Debug, Clone)]
pub struct CgForceField {
    /// Bond list and parameters.
    pub bonds: Vec<CgBond>,
    /// Angle list and parameters.
    pub angles: Vec<CgAngle>,
    /// Pair potentials per bead-type pair.
    pub pair_potentials: Vec<(u32, u32, CgPairPotential)>,
}
impl CgForceField {
    /// Create an empty CG force field.
    pub fn new() -> Self {
        CgForceField {
            bonds: Vec::new(),
            angles: Vec::new(),
            pair_potentials: Vec::new(),
        }
    }
    /// Add a bond.
    pub fn add_bond(&mut self, bond: CgBond) {
        self.bonds.push(bond);
    }
    /// Add an angle.
    pub fn add_angle(&mut self, angle: CgAngle) {
        self.angles.push(angle);
    }
    /// Add a pair potential for bead types (ta, tb).
    pub fn add_pair_potential(&mut self, ta: u32, tb: u32, pot: CgPairPotential) {
        self.pair_potentials.push((ta, tb, pot));
    }
    /// Look up pair potential for two bead types (order-insensitive).
    pub fn get_pair_potential(&self, ta: u32, tb: u32) -> Option<&CgPairPotential> {
        self.pair_potentials
            .iter()
            .find(|(a, b, _)| (*a == ta && *b == tb) || (*a == tb && *b == ta))
            .map(|(_, _, pot)| pot)
    }
    /// Total bonded energy for given bead positions.
    pub fn bonded_energy(&self, positions: &[[f64; 3]]) -> f64 {
        cg_bond_energy(positions, &self.bonds) + cg_angle_energy(positions, &self.angles)
    }
    /// Total nonbonded energy for a set of beads with type information.
    pub fn nonbonded_energy(&self, beads: &[CgBead]) -> f64 {
        let n = beads.len();
        let mut e = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(pot) =
                    self.get_pair_potential(beads[i].bead_type as u32, beads[j].bead_type as u32)
                {
                    let dx = beads[j].position[0] - beads[i].position[0];
                    let dy = beads[j].position[1] - beads[i].position[1];
                    let dz = beads[j].position[2] - beads[i].position[2];
                    let r = (dx * dx + dy * dy + dz * dz).sqrt();
                    e += pot.energy(r);
                }
            }
        }
        e
    }
}
/// A cosine dihedral (torsion) potential for four CG beads.
///
/// E = k * (1 + cos(n*φ - φ_0))
pub struct CgDihedral {
    /// Index of first bead.
    pub i: usize,
    /// Index of second bead.
    pub j: usize,
    /// Index of third bead.
    pub k: usize,
    /// Index of fourth bead.
    pub l: usize,
    /// Force constant (kJ/mol).
    pub k_phi: f64,
    /// Multiplicity n.
    pub n: u32,
    /// Phase angle φ_0 (radians).
    pub phi0: f64,
}
/// Accumulator for force matching data collection.
pub struct ForceMatchAccumulator {
    /// Minimum distance (nm).
    pub r_min: f64,
    /// Maximum distance (nm).
    pub r_max: f64,
    /// Number of bins.
    pub n_bins: usize,
    /// Bin width (nm).
    pub dr: f64,
    /// Sum of forces in each bin.
    pub(super) force_sums: Vec<f64>,
    /// Count of samples in each bin.
    pub(super) counts: Vec<usize>,
}
impl ForceMatchAccumulator {
    /// Create a new force matching accumulator.
    pub fn new(r_min: f64, r_max: f64, n_bins: usize) -> Self {
        let dr = (r_max - r_min) / n_bins as f64;
        Self {
            r_min,
            r_max,
            n_bins,
            dr,
            force_sums: vec![0.0; n_bins],
            counts: vec![0; n_bins],
        }
    }
    /// Add a force sample at a given pair distance.
    ///
    /// The radial component of the force (F·r̂) is accumulated.
    pub fn add_sample(&mut self, distance: f64, radial_force: f64) {
        if distance < self.r_min || distance >= self.r_max {
            return;
        }
        let bin = ((distance - self.r_min) / self.dr) as usize;
        if bin < self.n_bins {
            self.force_sums[bin] += radial_force;
            self.counts[bin] += 1;
        }
    }
    /// Compute the mean force in each bin.
    pub fn compute_mean_forces(&self) -> ForceMatchResult {
        let bin_edges: Vec<f64> = (0..=self.n_bins)
            .map(|i| self.r_min + i as f64 * self.dr)
            .collect();
        let mean_forces: Vec<f64> = self
            .force_sums
            .iter()
            .zip(self.counts.iter())
            .map(
                |(&sum, &count)| {
                    if count > 0 { sum / count as f64 } else { 0.0 }
                },
            )
            .collect();
        ForceMatchResult {
            bin_edges,
            mean_forces,
            counts: self.counts.clone(),
        }
    }
    /// Reset all accumulated data.
    pub fn reset(&mut self) {
        self.force_sums.fill(0.0);
        self.counts.fill(0);
    }
}
/// Simple CG velocity-Verlet integrator configuration.
pub struct CgIntegratorConfig {
    /// Time step (ps).
    pub dt: f64,
    /// Number of integration steps.
    pub n_steps: usize,
    /// Thermostat target temperature (K).
    pub temperature: f64,
    /// Thermostat time constant (ps); 0 = no thermostat.
    pub tau_t: f64,
    /// Boltzmann constant in kJ/(mol·K).
    pub kb: f64,
}
impl CgIntegratorConfig {
    /// Create a default config.
    pub fn new(dt: f64, n_steps: usize, temperature: f64) -> Self {
        Self {
            dt,
            n_steps,
            temperature,
            tau_t: 0.0,
            kb: 0.008_314_462_618,
        }
    }
    /// kT in kJ/mol.
    pub fn kt(&self) -> f64 {
        self.kb * self.temperature
    }
}
/// MARTINI bonded parameters between two bonded bead types.
#[derive(Debug, Clone, Copy)]
pub struct MartiniBondParams {
    /// Equilibrium bond length b₀ (nm).
    pub b0: f64,
    /// Spring constant K_b (kJ mol^-1 nm^-2).
    pub kb: f64,
}
/// Gaussian network model (GNM) — 1D version of the ANM/ENM.
///
/// The GNM uses a scalar (not tensorial) connectivity matrix Γ.
/// B-factors are proportional to the diagonal of the pseudo-inverse Γ⁺.
#[derive(Debug, Clone)]
pub struct GaussianNetworkModel {
    /// Kirchhoff (connectivity) matrix Γ, flattened row-major (n×n).
    pub kirchhoff: Vec<f64>,
    /// Number of nodes.
    pub n_nodes: usize,
}
impl GaussianNetworkModel {
    /// Build GNM from Cα positions and cutoff.
    pub fn build(positions: &[[f64; 3]], cutoff: f64) -> Self {
        let n = positions.len();
        let mut k = vec![0.0f64; n * n];
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = positions[j][0] - positions[i][0];
                let dy = positions[j][1] - positions[i][1];
                let dz = positions[j][2] - positions[i][2];
                let r2 = dx * dx + dy * dy + dz * dz;
                if r2 < cutoff * cutoff {
                    k[i * n + j] = -1.0;
                    k[j * n + i] = -1.0;
                    k[i * n + i] += 1.0;
                    k[j * n + j] += 1.0;
                }
            }
        }
        GaussianNetworkModel {
            kirchhoff: k,
            n_nodes: n,
        }
    }
    /// Diagonal of Γ (degree of each node = number of contacts).
    pub fn degree(&self) -> Vec<f64> {
        let n = self.n_nodes;
        (0..n).map(|i| self.kirchhoff[i * n + i]).collect()
    }
    /// Relative B-factor estimate: 1/degree_i (proxy for true Γ⁺ diagonal).
    pub fn relative_b_factors(&self) -> Vec<f64> {
        self.degree()
            .iter()
            .map(|&d| if d > 0.0 { 1.0 / d } else { f64::INFINITY })
            .collect()
    }
}
/// Iterative Boltzmann inversion (IBI) for deriving CG potentials.
///
/// Starting from a Boltzmann-inverted target RDF:
/// U_0(r) = -k_B T ln(g_target(r))
///
/// Iteratively corrects the potential:
/// U_{n+1}(r) = U_n(r) + k_B T ln(g_n(r) / g_target(r))
pub struct IbiPotential {
    /// Bin centers for distance (nm).
    pub r_values: Vec<f64>,
    /// Current potential values (kJ/mol).
    pub potential: Vec<f64>,
    /// Target RDF g_target(r).
    pub target_rdf: Vec<f64>,
    /// Thermal energy k_B T (kJ/mol).
    pub kbt: f64,
}
impl IbiPotential {
    /// Create an initial IBI potential from a target RDF.
    ///
    /// U_0(r) = -k_B T ln(g_target(r))
    ///
    /// # Arguments
    /// * `r_values` - Distance bin centers (nm)
    /// * `target_rdf` - Target radial distribution function
    /// * `kbt` - Thermal energy k_B T (kJ/mol), e.g., 2.494 at 300 K
    pub fn from_rdf(r_values: Vec<f64>, target_rdf: Vec<f64>, kbt: f64) -> Self {
        let potential: Vec<f64> = target_rdf
            .iter()
            .map(|&g| if g > 1e-10 { -kbt * g.ln() } else { kbt * 20.0 })
            .collect();
        Self {
            r_values,
            potential,
            target_rdf,
            kbt,
        }
    }
    /// Perform one IBI update step.
    ///
    /// U_{n+1}(r) = U_n(r) + k_B T ln(g_current(r) / g_target(r))
    ///
    /// # Arguments
    /// * `current_rdf` - RDF measured from the current CG simulation
    /// * `damping` - Damping factor ∈ (0, 1] for stability (1.0 = no damping)
    pub fn update(&mut self, current_rdf: &[f64], damping: f64) {
        let damping = damping.clamp(0.01, 1.0);
        for (i, (u, (&g_curr, &g_targ))) in self
            .potential
            .iter_mut()
            .zip(current_rdf.iter().zip(self.target_rdf.iter()))
            .enumerate()
        {
            let _ = i;
            if g_curr > 1e-10 && g_targ > 1e-10 {
                *u += damping * self.kbt * (g_curr / g_targ).ln();
            }
        }
    }
    /// Compute the potential at a given distance via linear interpolation.
    pub fn interpolate(&self, r: f64) -> f64 {
        if self.r_values.is_empty() {
            return 0.0;
        }
        if r <= self.r_values[0] {
            return self.potential[0];
        }
        let last = self.r_values.len() - 1;
        if r >= self.r_values[last] {
            return self.potential[last];
        }
        for i in 0..last {
            if r >= self.r_values[i] && r < self.r_values[i + 1] {
                let t = (r - self.r_values[i]) / (self.r_values[i + 1] - self.r_values[i]);
                return self.potential[i] * (1.0 - t) + self.potential[i + 1] * t;
            }
        }
        self.potential[last]
    }
    /// Compute the force (negative derivative) via finite differences.
    pub fn force(&self, r: f64) -> f64 {
        let dr = 1e-5;
        let u_plus = self.interpolate(r + dr);
        let u_minus = self.interpolate(r - dr);
        -(u_plus - u_minus) / (2.0 * dr)
    }
    /// Convergence metric: RMS relative difference between current and target RDF.
    pub fn convergence_metric(&self, current_rdf: &[f64]) -> f64 {
        let mut sum_sq = 0.0;
        let mut count = 0;
        for (&g_curr, &g_targ) in current_rdf.iter().zip(self.target_rdf.iter()) {
            if g_targ > 1e-10 {
                let rel_diff = (g_curr - g_targ) / g_targ;
                sum_sq += rel_diff * rel_diff;
                count += 1;
            }
        }
        if count > 0 {
            (sum_sq / count as f64).sqrt()
        } else {
            0.0
        }
    }
}
/// An elastic network model (ENM) contact.
///
/// Two Cα atoms are connected by a spring if their distance is below
/// the cutoff radius.  The spring constant is uniform across all contacts
/// in the standard Tirion ENM.
#[derive(Debug, Clone)]
pub struct EnmContact {
    /// Index of the first Cα atom.
    pub i: usize,
    /// Index of the second Cα atom.
    pub j: usize,
    /// Equilibrium distance r₀ (nm).
    pub r0: f64,
    /// Spring constant γ (kJ mol^-1 nm^-2).
    pub gamma: f64,
}
/// CG velocity-Verlet integrator state.
pub struct CgIntegratorState {
    /// Positions of each bead \[n_beads × 3\].
    pub positions: Vec<[f64; 3]>,
    /// Velocities of each bead \[n_beads × 3\].
    pub velocities: Vec<[f64; 3]>,
    /// Masses of each bead (amu).
    pub masses: Vec<f64>,
    /// Current forces \[n_beads × 3\] (kJ mol⁻¹ nm⁻¹).
    pub forces: Vec<[f64; 3]>,
    /// Current simulation time (ps).
    pub time: f64,
    /// Step counter.
    pub step: usize,
}
impl CgIntegratorState {
    /// Create a new integrator state.
    pub fn new(positions: Vec<[f64; 3]>, velocities: Vec<[f64; 3]>, masses: Vec<f64>) -> Self {
        let n = positions.len();
        Self {
            positions,
            velocities,
            masses,
            forces: vec![[0.0; 3]; n],
            time: 0.0,
            step: 0,
        }
    }
    /// Perform a single velocity-Verlet half-step for velocities (v += 0.5*a*dt).
    pub fn velocity_half_step(&mut self, dt: f64) {
        let n = self.positions.len();
        for i in 0..n {
            let inv_mass = if self.masses[i] > 1e-300 {
                1.0 / self.masses[i]
            } else {
                0.0
            };
            for d in 0..3 {
                self.velocities[i][d] += 0.5 * self.forces[i][d] * inv_mass * dt;
            }
        }
    }
    /// Perform a full position update step (x += v*dt).
    pub fn position_step(&mut self, dt: f64) {
        let n = self.positions.len();
        for i in 0..n {
            for d in 0..3 {
                self.positions[i][d] += self.velocities[i][d] * dt;
            }
        }
    }
    /// Kinetic energy (kJ/mol) assuming velocities in nm/ps, masses in g/mol (≈ amu×kJ factor).
    ///
    /// Uses K = ½ Σ m_i |v_i|², with a unit-conversion factor of 1/1000 (amu nm²/ps² → kJ/mol).
    pub fn kinetic_energy(&self) -> f64 {
        let unit = 1.0 / 1000.0;
        self.masses
            .iter()
            .zip(self.velocities.iter())
            .map(|(&m, v)| 0.5 * m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) * unit)
            .sum()
    }
    /// Instantaneous temperature (K) from equipartition.
    pub fn instantaneous_temperature(&self, n_dof: f64) -> f64 {
        let kb = 0.008_314_462_618_f64;
        if n_dof <= 0.0 || kb <= 0.0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (n_dof * kb)
    }
    /// Rescale velocities to achieve target temperature (velocity rescaling thermostat).
    pub fn rescale_velocities(&mut self, target_temp: f64, n_dof: f64) {
        let current_temp = self.instantaneous_temperature(n_dof);
        if current_temp < 1e-12 {
            return;
        }
        let scale = (target_temp / current_temp).sqrt();
        for v in &mut self.velocities {
            for x in v.iter_mut() {
                *x *= scale;
            }
        }
    }
    /// Apply periodic boundary conditions (minimum image convention).
    ///
    /// # Arguments
    /// * `box_lengths` – simulation box lengths \[lx, ly, lz\] (nm).
    pub fn apply_pbc(&mut self, box_lengths: [f64; 3]) {
        for pos in &mut self.positions {
            for d in 0..3 {
                let l = box_lengths[d];
                if l > 0.0 {
                    pos[d] -= l * (pos[d] / l).floor();
                }
            }
        }
    }
    /// Centre-of-mass position of the system.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let mut com = [0.0_f64; 3];
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-300 {
            return com;
        }
        for (i, pos) in self.positions.iter().enumerate() {
            let m = self.masses[i];
            for d in 0..3 {
                com[d] += m * pos[d];
            }
        }
        for v in &mut com {
            *v /= total_mass;
        }
        com
    }
    /// Centre-of-mass velocity.
    pub fn centre_of_mass_velocity(&self) -> [f64; 3] {
        let mut com_vel = [0.0_f64; 3];
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass < 1e-300 {
            return com_vel;
        }
        for (i, vel) in self.velocities.iter().enumerate() {
            let m = self.masses[i];
            for d in 0..3 {
                com_vel[d] += m * vel[d];
            }
        }
        for v in &mut com_vel {
            *v /= total_mass;
        }
        com_vel
    }
    /// Remove centre-of-mass translation (zero net momentum).
    pub fn remove_com_velocity(&mut self) {
        let com_vel = self.centre_of_mass_velocity();
        for vel in &mut self.velocities {
            for d in 0..3 {
                vel[d] -= com_vel[d];
            }
        }
    }
}
/// MARTINI angle parameters.
#[derive(Debug, Clone, Copy)]
pub struct MartiniAngleParams {
    /// Equilibrium angle θ₀ (degrees).
    pub theta0_deg: f64,
    /// Force constant K_θ (kJ mol^-1 rad^-2).
    pub k_theta: f64,
}
impl MartiniAngleParams {
    /// Bond bending energy for a given angle (radians).
    pub fn energy(&self, theta_rad: f64) -> f64 {
        let theta0 = self.theta0_deg.to_radians();
        let dtheta = theta_rad - theta0;
        0.5 * self.k_theta * dtheta * dtheta
    }
}
