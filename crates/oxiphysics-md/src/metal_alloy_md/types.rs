//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// Stores per-atom initial positions for MSD tracking.
#[derive(Debug, Clone)]
pub struct DiffusionTracker {
    /// Initial (unwrapped) positions.
    pub initial_positions: Vec<[f64; 3]>,
    /// Current unwrapped positions.
    pub unwrapped_positions: Vec<[f64; 3]>,
    /// Previous wrapped positions (for unwrapping).
    pub(super) prev_wrapped: Vec<[f64; 3]>,
    /// Species filter (None = all).
    pub species_filter: Option<usize>,
}
impl DiffusionTracker {
    /// Create a new tracker from system.
    pub fn new(system: &AlloySystem, species_filter: Option<usize>) -> Self {
        let positions: Vec<[f64; 3]> = system.atoms.iter().map(|a| a.position).collect();
        Self {
            initial_positions: positions.clone(),
            unwrapped_positions: positions.clone(),
            prev_wrapped: positions,
            species_filter,
        }
    }
    /// Update unwrapped positions from current system state.
    pub fn update(&mut self, system: &AlloySystem) {
        for (idx, atom) in system.atoms.iter().enumerate() {
            let dr_wrapped = vsub(atom.position, self.prev_wrapped[idx]);
            let dr = system.sim_box.min_image(dr_wrapped);
            self.unwrapped_positions[idx] = vadd(self.unwrapped_positions[idx], dr);
            self.prev_wrapped[idx] = atom.position;
        }
    }
    /// Compute mean square displacement (ų).
    pub fn msd(&self, system: &AlloySystem) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for (idx, atom) in system.atoms.iter().enumerate() {
            if let Some(sf) = self.species_filter
                && atom.species != sf
            {
                continue;
            }
            let dr = vsub(self.unwrapped_positions[idx], self.initial_positions[idx]);
            sum += vdot(dr, dr);
            count += 1;
        }
        if count == 0 {
            return 0.0;
        }
        sum / count as f64
    }
    /// Compute self-diffusion coefficient from MSD (ų/ps).
    /// D = MSD / (6 * t) for 3D.
    pub fn diffusion_coefficient(&self, system: &AlloySystem, time: f64) -> f64 {
        if time <= 0.0 {
            return 0.0;
        }
        self.msd(system) / (6.0 * time)
    }
}
/// CNA signature type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CnaSignature {
    /// Number of common neighbors.
    pub n_common: u32,
    /// Number of bonds among common neighbors.
    pub n_bonds: u32,
    /// Longest chain among those bonds.
    pub longest_chain: u32,
}
/// Alloy system with multiple species.
#[derive(Debug, Clone)]
pub struct AlloySystem {
    /// EAM parameters per species.
    pub species: Vec<EamParams>,
    /// Atoms.
    pub atoms: Vec<MetalAtom>,
    /// Simulation box.
    pub sim_box: SimBox,
    /// Mixing rule for cross-species interactions.
    pub mixing_rule: MixingRule,
    /// Timestep (ps).
    pub dt: f64,
}
impl AlloySystem {
    /// Create a new alloy system.
    pub fn new(species: Vec<EamParams>, sim_box: SimBox, mixing_rule: MixingRule, dt: f64) -> Self {
        Self {
            species,
            atoms: Vec::new(),
            sim_box,
            mixing_rule,
            dt,
        }
    }
    /// Add an atom to the system.
    pub fn add_atom(&mut self, position: [f64; 3], species: usize) {
        self.atoms.push(MetalAtom::new(position, species));
    }
    /// Number of atoms.
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }
    /// Get maximum cutoff across all species.
    pub fn max_cutoff(&self) -> f64 {
        self.species
            .iter()
            .map(|s| s.cutoff)
            .fold(0.0_f64, f64::max)
    }
    /// Build simple neighbor list (O(N^2), brute force).
    /// Returns Vec of (i, j, dr_vec, dist) with i < j.
    pub fn build_neighbor_list(&self) -> Vec<(usize, usize, [f64; 3], f64)> {
        let rc = self.max_cutoff();
        let n = self.atoms.len();
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dr = self
                    .sim_box
                    .min_image(vsub(self.atoms[j].position, self.atoms[i].position));
                let r = vnorm(dr);
                if r < rc && r > 0.1 {
                    pairs.push((i, j, dr, r));
                }
            }
        }
        pairs
    }
    /// Compute forces using EAM.
    pub fn compute_forces(&mut self) {
        let n = self.atoms.len();
        for atom in &mut self.atoms {
            atom.force = [0.0; 3];
            atom.rho_bar = 0.0;
        }
        let pairs = self.build_neighbor_list();
        for &(i, j, _dr, r) in &pairs {
            let si = self.atoms[i].species;
            let sj = self.atoms[j].species;
            self.atoms[i].rho_bar += self.species[sj].electron_density(r);
            self.atoms[j].rho_bar += self.species[si].electron_density(r);
        }
        for atom in &mut self.atoms {
            let s = atom.species;
            atom.df_drho = self.species[s].embedding_deriv(atom.rho_bar);
        }
        let mut force_inc: Vec<[f64; 3]> = vec![[0.0; 3]; n];
        for &(i, j, dr, r) in &pairs {
            let si = self.atoms[i].species;
            let sj = self.atoms[j].species;
            let pair_f = if si == sj {
                self.species[si].pair_force(r)
            } else {
                cross_pair_force(&self.species[si], &self.species[sj], r, self.mixing_rule)
            };
            let drho_ij = self.species[sj].electron_density_deriv(r);
            let drho_ji = self.species[si].electron_density_deriv(r);
            let eam_f = self.atoms[i].df_drho * drho_ij + self.atoms[j].df_drho * drho_ji;
            let f_mag = -(pair_f + eam_f);
            let rhat = vscale(dr, 1.0 / r);
            let fvec = vscale(rhat, f_mag);
            force_inc[i] = vadd(force_inc[i], fvec);
            force_inc[j] = vsub(force_inc[j], fvec);
        }
        for (k, inc) in force_inc.iter().enumerate().take(n) {
            self.atoms[k].force = vadd(self.atoms[k].force, *inc);
        }
    }
    /// Compute total potential energy (eV).
    pub fn potential_energy(&self) -> f64 {
        let pairs = self.build_neighbor_list();
        let mut e_pair = 0.0;
        for &(i, j, _dr, r) in &pairs {
            let si = self.atoms[i].species;
            let sj = self.atoms[j].species;
            if si == sj {
                e_pair += self.species[si].pair_potential(r);
            } else {
                e_pair +=
                    cross_pair_potential(&self.species[si], &self.species[sj], r, self.mixing_rule);
            }
        }
        let mut e_embed = 0.0;
        for atom in &self.atoms {
            e_embed += self.species[atom.species].embedding_energy(atom.rho_bar);
        }
        e_pair + e_embed
    }
    /// Total kinetic energy (eV).
    pub fn kinetic_energy(&self) -> f64 {
        self.atoms
            .iter()
            .map(|a| a.kinetic_energy(self.species[a.species].mass))
            .sum()
    }
    /// Instantaneous temperature (K).
    pub fn temperature(&self) -> f64 {
        let ke = self.kinetic_energy();
        let ndof = (3 * self.atoms.len()) as f64;
        2.0 * ke / (ndof * 8.617333e-5)
    }
    /// Velocity Verlet integration step.
    pub fn step(&mut self) {
        let dt = self.dt;
        let half_dt = 0.5 * dt;
        for atom in &mut self.atoms {
            let mass = 1.0;
            let acc = vscale(atom.force, 1.0 / mass);
            atom.velocity = vadd(atom.velocity, vscale(acc, half_dt));
            atom.position = vadd(atom.position, vscale(atom.velocity, dt));
            atom.position = self.sim_box.wrap(atom.position);
        }
        self.compute_forces();
        for atom in &mut self.atoms {
            let mass = 1.0;
            let acc = vscale(atom.force, 1.0 / mass);
            atom.velocity = vadd(atom.velocity, vscale(acc, half_dt));
        }
    }
    /// Run multiple steps.
    pub fn run(&mut self, nsteps: usize) {
        self.compute_forces();
        for _ in 0..nsteps {
            self.step();
        }
    }
}
/// Result of dislocation core analysis.
#[derive(Debug, Clone)]
pub struct DislocationCore {
    /// Burger's vector (Å).
    pub burgers_vector: [f64; 3],
    /// Core position estimate (Å).
    pub core_center: [f64; 3],
    /// Core width estimate (Å).
    pub core_width: f64,
    /// Displacement field (atom_index, displacement).
    pub displacements: Vec<(usize, [f64; 3])>,
}
/// Periodic boundary conditions box.
#[derive(Debug, Clone)]
pub struct SimBox {
    /// Box lengths along x, y, z (Å).
    pub lengths: [f64; 3],
}
impl SimBox {
    /// Create a new orthorhombic box.
    pub fn new(lx: f64, ly: f64, lz: f64) -> Self {
        Self {
            lengths: [lx, ly, lz],
        }
    }
    /// Cubic box.
    pub fn cubic(l: f64) -> Self {
        Self::new(l, l, l)
    }
    /// Apply minimum image convention.
    pub fn min_image(&self, dr: [f64; 3]) -> [f64; 3] {
        let mut out = dr;
        for (out_k, &len) in out.iter_mut().zip(self.lengths.iter()) {
            while *out_k > 0.5 * len {
                *out_k -= len;
            }
            while *out_k < -0.5 * len {
                *out_k += len;
            }
        }
        out
    }
    /// Wrap position into box.
    pub fn wrap(&self, pos: [f64; 3]) -> [f64; 3] {
        let mut out = pos;
        for (out_k, &len) in out.iter_mut().zip(self.lengths.iter()) {
            while *out_k >= len {
                *out_k -= len;
            }
            while *out_k < 0.0 {
                *out_k += len;
            }
        }
        out
    }
    /// Volume of the box (ų).
    pub fn volume(&self) -> f64 {
        self.lengths[0] * self.lengths[1] * self.lengths[2]
    }
}
/// A metal atom in the alloy simulation.
#[derive(Debug, Clone)]
pub struct MetalAtom {
    /// Position (Å).
    pub position: [f64; 3],
    /// Velocity (Å/ps).
    pub velocity: [f64; 3],
    /// Force (eV/Å).
    pub force: [f64; 3],
    /// Species index (into AlloySystem::species).
    pub species: usize,
    /// Local electron density (accumulated).
    pub rho_bar: f64,
    /// Embedding derivative dF/dρ̄.
    pub df_drho: f64,
}
impl MetalAtom {
    /// Create a new atom of given species at position.
    pub fn new(position: [f64; 3], species: usize) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            species,
            rho_bar: 0.0,
            df_drho: 0.0,
        }
    }
    /// Kinetic energy (eV).  mass in amu, velocity in Å/ps.
    /// KE = 0.5 * m * v^2 * (1 amu * (1e-12 m/s)^2 -> eV conversion).
    pub fn kinetic_energy(&self, mass: f64) -> f64 {
        let v2 = vdot(self.velocity, self.velocity);
        0.5 * mass * v2 * 1.03643e-4
    }
}
/// Mixing rule for cross-species EAM interactions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MixingRule {
    /// Geometric mean for pair parameters.
    Geometric,
    /// Arithmetic mean for pair parameters.
    Arithmetic,
    /// Johnson (1989) universal mixing.
    Johnson,
}
/// Crystal structure type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrystalStructure {
    /// Face-centered cubic.
    FCC,
    /// Body-centered cubic.
    BCC,
    /// Hexagonal close-packed.
    HCP,
}
/// Parameters for a single-element EAM potential.
///
/// Uses the Finnis-Sinclair form:
///   pair potential  φ(r) = A * exp(-α*(r/re - 1)) - B * exp(-β*(r/re - 1))
///   electron density ρ(r) = fe * exp(-β*(r/re - 1))
///   embedding  F(ρ̄) = -F0 * sqrt(ρ̄)   (Finnis-Sinclair)
#[derive(Debug, Clone)]
pub struct EamParams {
    /// Equilibrium nearest-neighbor distance (Å).
    pub re: f64,
    /// Pair repulsion amplitude (eV).
    pub a_coeff: f64,
    /// Pair attraction amplitude (eV).
    pub b_coeff: f64,
    /// Repulsive range parameter.
    pub alpha: f64,
    /// Attractive / density range parameter.
    pub beta: f64,
    /// Electron density scaling.
    pub fe: f64,
    /// Embedding energy scale (eV).
    pub f0: f64,
    /// Cutoff distance (Å).
    pub cutoff: f64,
    /// Atomic mass (amu).
    pub mass: f64,
    /// Lattice constant (Å).
    pub lattice_a: f64,
    /// Element label (e.g. "Cu", "Ni").
    pub label: String,
}
impl EamParams {
    /// Create EAM parameters for copper (Cu).
    pub fn copper() -> Self {
        Self {
            re: 2.556,
            a_coeff: 0.396,
            b_coeff: 0.548,
            alpha: 5.09,
            beta: 3.63,
            fe: 1.554,
            f0: 1.845,
            cutoff: 5.50,
            mass: 63.546,
            lattice_a: 3.615,
            label: "Cu".to_string(),
        }
    }
    /// Create EAM parameters for nickel (Ni).
    pub fn nickel() -> Self {
        Self {
            re: 2.490,
            a_coeff: 0.440,
            b_coeff: 0.633,
            alpha: 5.21,
            beta: 3.73,
            fe: 1.669,
            f0: 2.007,
            cutoff: 5.30,
            mass: 58.693,
            lattice_a: 3.524,
            label: "Ni".to_string(),
        }
    }
    /// Create EAM parameters for aluminum (Al).
    pub fn aluminum() -> Self {
        Self {
            re: 2.864,
            a_coeff: 0.314,
            b_coeff: 0.365,
            alpha: 4.61,
            beta: 2.81,
            fe: 1.070,
            f0: 1.178,
            cutoff: 6.00,
            mass: 26.982,
            lattice_a: 4.050,
            label: "Al".to_string(),
        }
    }
    /// Create EAM parameters for iron (Fe, BCC).
    pub fn iron() -> Self {
        Self {
            re: 2.481,
            a_coeff: 0.392,
            b_coeff: 0.418,
            alpha: 5.12,
            beta: 3.44,
            fe: 1.885,
            f0: 1.692,
            cutoff: 5.30,
            mass: 55.845,
            lattice_a: 2.8665,
            label: "Fe".to_string(),
        }
    }
    /// Pair potential φ(r).
    pub fn pair_potential(&self, r: f64) -> f64 {
        if r >= self.cutoff || r < 0.5 {
            return 0.0;
        }
        let x = r / self.re - 1.0;
        self.a_coeff * (-self.alpha * x).exp() - self.b_coeff * (-self.beta * x).exp()
    }
    /// Derivative of pair potential dφ/dr.
    pub fn pair_force(&self, r: f64) -> f64 {
        if r >= self.cutoff || r < 0.5 {
            return 0.0;
        }
        let x = r / self.re - 1.0;
        let inv_re = 1.0 / self.re;
        -self.a_coeff * self.alpha * inv_re * (-self.alpha * x).exp()
            + self.b_coeff * self.beta * inv_re * (-self.beta * x).exp()
    }
    /// Electron density contribution ρ(r).
    pub fn electron_density(&self, r: f64) -> f64 {
        if r >= self.cutoff || r < 0.5 {
            return 0.0;
        }
        let x = r / self.re - 1.0;
        self.fe * (-self.beta * x).exp()
    }
    /// Derivative dρ/dr.
    pub fn electron_density_deriv(&self, r: f64) -> f64 {
        if r >= self.cutoff || r < 0.5 {
            return 0.0;
        }
        let x = r / self.re - 1.0;
        -self.fe * self.beta / self.re * (-self.beta * x).exp()
    }
    /// Embedding energy F(ρ̄) = -F0 * sqrt(ρ̄).
    pub fn embedding_energy(&self, rho_bar: f64) -> f64 {
        if rho_bar <= 0.0 {
            return 0.0;
        }
        -self.f0 * rho_bar.sqrt()
    }
    /// Derivative dF/dρ̄.
    pub fn embedding_deriv(&self, rho_bar: f64) -> f64 {
        if rho_bar <= 1e-30 {
            return 0.0;
        }
        -self.f0 / (2.0 * rho_bar.sqrt())
    }
}
