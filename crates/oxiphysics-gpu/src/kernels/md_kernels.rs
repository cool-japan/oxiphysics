// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-mock Molecular Dynamics (MD) compute kernels.
//!
//! Mirrors GPU dispatch layout but executes in pure Rust on the CPU.


use std::f64::consts::PI;

// ── MdKernelConfig ────────────────────────────────────────────────────────────

/// Configuration shared across all MD kernels.
#[derive(Debug, Clone)]
pub struct MdKernelConfig {
    /// Lennard-Jones cut-off radius.
    pub cutoff: f64,
    /// Verlet list skin distance (Δ beyond cutoff).
    pub skin_distance: f64,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Simulation box dimensions \[Lx, Ly, Lz\].
    pub box_size: [f64; 3],
    /// Whether to apply periodic boundary conditions.
    pub periodic: bool,
}

impl MdKernelConfig {
    /// Create a default config for a 100-atom LJ gas.
    pub fn new_lj_gas(n_atoms: usize, box_size: [f64; 3]) -> Self {
        Self {
            cutoff: 2.5,
            skin_distance: 0.3,
            n_atoms,
            box_size,
            periodic: true,
        }
    }

    /// Neighbour-list update radius (cutoff + skin).
    pub fn list_cutoff(&self) -> f64 {
        self.cutoff + self.skin_distance
    }

    /// Box volume.
    pub fn volume(&self) -> f64 {
        self.box_size[0] * self.box_size[1] * self.box_size[2]
    }

    /// Number density.
    pub fn number_density(&self) -> f64 {
        self.n_atoms as f64 / self.volume()
    }
}

// ── Helper math ───────────────────────────────────────────────────────────────

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 { a[0]*b[0]+a[1]*b[1]+a[2]*b[2] }
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]-b[0],a[1]-b[1],a[2]-b[2]] }
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0]+b[0],a[1]+b[1],a[2]+b[2]] }
fn scale3(s: f64, v: [f64; 3]) -> [f64; 3] { [s*v[0],s*v[1],s*v[2]] }
fn len3(v: [f64; 3]) -> f64 { dot3(v,v).sqrt() }

/// Apply minimum image convention to a displacement vector.
fn minimum_image(dr: [f64; 3], box_size: [f64; 3]) -> [f64; 3] {
    let mut r = dr;
    for k in 0..3 {
        let l = box_size[k];
        if r[k] > 0.5 * l  { r[k] -= l; }
        if r[k] < -0.5 * l { r[k] += l; }
    }
    r
}

// ── LJ parameters ─────────────────────────────────────────────────────────────

/// Lennard-Jones pair parameters.
#[derive(Debug, Clone, Copy)]
pub struct LjParams {
    /// Energy well depth ε.
    pub epsilon: f64,
    /// Finite distance σ at which pair potential is zero.
    pub sigma: f64,
}

impl LjParams {
    /// Create new LJ parameters.
    pub fn new(epsilon: f64, sigma: f64) -> Self {
        Self { epsilon, sigma }
    }

    /// Standard Argon-like parameters (ε/k_B = 119.8 K, σ = 3.405 Å).
    pub fn argon() -> Self {
        Self { epsilon: 1.0, sigma: 1.0 } // in reduced LJ units
    }
}

// ── LjForceKernel ─────────────────────────────────────────────────────────────

/// Kernel for Lennard-Jones pairwise force computation.
#[derive(Debug, Clone)]
pub struct LjForceKernel {
    /// LJ parameters.
    pub params: LjParams,
    /// Cut-off radius.
    pub cutoff: f64,
    /// Box size for periodic images.
    pub box_size: [f64; 3],
    /// Use periodic boundary conditions.
    pub periodic: bool,
}

impl LjForceKernel {
    /// Create a new LJ force kernel.
    pub fn new(params: LjParams, cutoff: f64, box_size: [f64; 3], periodic: bool) -> Self {
        Self { params, cutoff, box_size, periodic }
    }

    /// LJ potential energy between two particles at distance r.
    pub fn lj_potential(&self, r: f64) -> f64 {
        let s_over_r = self.params.sigma / r;
        let sr6 = s_over_r.powi(6);
        let sr12 = sr6 * sr6;
        4.0 * self.params.epsilon * (sr12 - sr6)
    }

    /// LJ force magnitude df/dr * (-1/r) applied along r-vector.
    pub fn lj_force_scalar(&self, r: f64) -> f64 {
        let s_over_r = self.params.sigma / r;
        let sr6 = s_over_r.powi(6);
        let sr12 = sr6 * sr6;
        24.0 * self.params.epsilon / r * (2.0 * sr12 - sr6)
    }

    /// Compute all LJ forces and total potential energy.
    ///
    /// `positions` — flat list of positions, length n_atoms * 3.
    /// Returns (forces: `Vec<[f64;3]>`, potential_energy: f64).
    pub fn compute_lj_forces(
        &self,
        positions: &[[f64; 3]],
    ) -> (Vec<[f64; 3]>, f64) {
        let n = positions.len();
        let mut forces = vec![[0.0f64; 3]; n];
        let mut energy = 0.0;
        let rc2 = self.cutoff * self.cutoff;

        for i in 0..n {
            for j in i+1..n {
                let mut dr = sub3(positions[i], positions[j]);
                if self.periodic {
                    dr = minimum_image(dr, self.box_size);
                }
                let r2 = dot3(dr, dr);
                if r2 > rc2 || r2 < 1e-20 { continue; }
                let r = r2.sqrt();
                let f_mag = self.lj_force_scalar(r);
                let f = scale3(f_mag / r, dr);
                forces[i] = add3(forces[i], f);
                forces[j] = [forces[j][0]-f[0], forces[j][1]-f[1], forces[j][2]-f[2]];
                energy += self.lj_potential(r);
            }
        }
        (forces, energy)
    }

    /// LJ virial (for pressure calculation): W = sum_ij r_ij * f_ij.
    pub fn compute_virial(&self, positions: &[[f64; 3]]) -> f64 {
        let n = positions.len();
        let mut virial = 0.0;
        let rc2 = self.cutoff * self.cutoff;
        for i in 0..n {
            for j in i+1..n {
                let mut dr = sub3(positions[i], positions[j]);
                if self.periodic { dr = minimum_image(dr, self.box_size); }
                let r2 = dot3(dr, dr);
                if r2 > rc2 || r2 < 1e-20 { continue; }
                let r = r2.sqrt();
                let f_mag = self.lj_force_scalar(r);
                virial += f_mag * r;
            }
        }
        virial / 3.0
    }
}

// ── EwaldKernel ───────────────────────────────────────────────────────────────

/// Kernel for Ewald summation electrostatics.
#[derive(Debug, Clone)]
pub struct EwaldKernel {
    /// Ewald splitting parameter α.
    pub alpha: f64,
    /// Real-space cut-off.
    pub cutoff: f64,
    /// Maximum k-vector in each direction.
    pub kmax: usize,
    /// Box size.
    pub box_size: [f64; 3],
}

impl EwaldKernel {
    /// Create a new Ewald kernel.
    pub fn new(alpha: f64, cutoff: f64, kmax: usize, box_size: [f64; 3]) -> Self {
        Self { alpha, cutoff, kmax, box_size }
    }

    /// Ewald real-space contribution for a pair at distance r.
    pub fn erfc_term(&self, r: f64, qi: f64, qj: f64) -> f64 {
        let k_e = 1.0; // Coulomb constant in reduced units
        k_e * qi * qj * erfc_approx(self.alpha * r) / r
    }

    /// Ewald real-space energy sum.
    pub fn ewald_real_space(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> f64 {
        let n = positions.len();
        let mut energy = 0.0;
        let rc2 = self.cutoff * self.cutoff;
        for i in 0..n {
            for j in i+1..n {
                let mut dr = sub3(positions[i], positions[j]);
                dr = minimum_image(dr, self.box_size);
                let r2 = dot3(dr, dr);
                if r2 > rc2 || r2 < 1e-20 { continue; }
                let r = r2.sqrt();
                energy += self.erfc_term(r, charges[i], charges[j]);
            }
        }
        energy
    }

    /// Ewald reciprocal-space energy sum.
    pub fn ewald_recip_space(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> f64 {
        let vol = self.box_size[0] * self.box_size[1] * self.box_size[2];
        let prefac = 4.0 * PI / vol;
        let mut energy = 0.0;

        let km = self.kmax as i64;
        for kx in -km..=km {
            for ky in -km..=km {
                for kz in -km..=km {
                    if kx == 0 && ky == 0 && kz == 0 { continue; }
                    let k = [
                        2.0 * PI * kx as f64 / self.box_size[0],
                        2.0 * PI * ky as f64 / self.box_size[1],
                        2.0 * PI * kz as f64 / self.box_size[2],
                    ];
                    let k2 = dot3(k, k);
                    if k2 < 1e-20 { continue; }
                    let factor = prefac * (-k2 / (4.0 * self.alpha * self.alpha)).exp() / k2;
                    // Structure factor |S(k)|²
                    let mut s_re = 0.0;
                    let mut s_im = 0.0;
                    for (pos, &q) in positions.iter().zip(charges.iter()) {
                        let kdotr = dot3(k, *pos);
                        s_re += q * kdotr.cos();
                        s_im += q * kdotr.sin();
                    }
                    energy += factor * (s_re * s_re + s_im * s_im);
                }
            }
        }
        0.5 * energy
    }

    /// Ewald self-energy correction.
    pub fn self_energy(&self, charges: &[f64]) -> f64 {
        let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
        -self.alpha / PI.sqrt() * sum_q2
    }
}

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 { return 2.0 - erfc_approx(-x); }
    let p = 0.3275911;
    let t = 1.0 / (1.0 + p * x);
    let poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741
        + t * (-1.453152027 + t * 1.061405429))));
    (-x * x).exp() * poly
}

// ── BondForceKernel ───────────────────────────────────────────────────────────

/// Kernel for covalent bond force computations.
#[derive(Debug, Clone)]
pub struct BondForceKernel;

impl BondForceKernel {
    /// Create a new bond force kernel.
    pub fn new() -> Self {
        Self
    }

    /// Harmonic bond: U = 0.5 * k * (r - r0)². Returns (energy, force magnitude).
    pub fn harmonic_bond(&self, r: f64, r0: f64, k: f64) -> (f64, f64) {
        let dr = r - r0;
        (0.5 * k * dr * dr, -k * dr)
    }

    /// Morse bond: U = D * (1 - exp(-alpha*(r-r0)))². Returns (energy, force).
    pub fn morse_bond(&self, r: f64, r0: f64, d: f64, alpha: f64) -> (f64, f64) {
        let x = 1.0 - (-alpha * (r - r0)).exp();
        let energy = d * x * x;
        let force = -2.0 * d * alpha * x * (-alpha * (r - r0)).exp();
        (energy, force)
    }

    /// FENE (Finitely Extensible Nonlinear Elastic) bond.
    /// U = -0.5 * k * R0² * ln(1 - (r/R0)²).
    /// Returns (energy, force). Returns NaN if r >= R0.
    pub fn fene_bond(&self, r: f64, r0: f64, k: f64, big_r0: f64) -> (f64, f64) {
        let dr = r - r0;
        let fene_arg = 1.0 - (dr / big_r0).powi(2);
        if fene_arg <= 0.0 { return (f64::INFINITY, f64::INFINITY); }
        let energy = -0.5 * k * big_r0 * big_r0 * fene_arg.ln();
        let force = -k * dr / fene_arg;
        (energy, force)
    }

    /// Compute bond force vectors for a list of bonded pairs.
    ///
    /// `pairs` — list of (i, j) atom index pairs.
    /// `r0_k`  — list of (r0, k) for each bond.
    /// Returns forces array and total energy.
    pub fn compute_harmonic_bond_forces(
        &self,
        positions: &[[f64; 3]],
        pairs: &[(usize, usize)],
        r0_k: &[(f64, f64)],
    ) -> (Vec<[f64; 3]>, f64) {
        let mut forces = vec![[0.0f64; 3]; positions.len()];
        let mut energy = 0.0;
        for (&(i, j), &(r0, k)) in pairs.iter().zip(r0_k.iter()) {
            let dr = sub3(positions[i], positions[j]);
            let r = len3(dr);
            if r < 1e-15 { continue; }
            let (e, f_mag) = self.harmonic_bond(r, r0, k);
            energy += e;
            let f = scale3(f_mag / r, dr);
            forces[i] = add3(forces[i], f);
            forces[j] = [forces[j][0]-f[0], forces[j][1]-f[1], forces[j][2]-f[2]];
        }
        (forces, energy)
    }
}

impl Default for BondForceKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── AngleForceKernel ──────────────────────────────────────────────────────────

/// Kernel for three-body angle force computations.
#[derive(Debug, Clone)]
pub struct AngleForceKernel;

impl AngleForceKernel {
    /// Create a new angle force kernel.
    pub fn new() -> Self {
        Self
    }

    /// Harmonic angle: U = 0.5 * k * (θ - θ₀)². Returns (energy, force on angle).
    pub fn harmonic_angle(&self, theta: f64, theta0: f64, k: f64) -> (f64, f64) {
        let dtheta = theta - theta0;
        (0.5 * k * dtheta * dtheta, -k * dtheta)
    }

    /// Urey-Bradley: harmonic + 1-3 distance term.
    /// Returns (energy, force magnitude for 1-3 distance).
    pub fn urey_bradley(&self, theta: f64, theta0: f64, k_theta: f64,
                         r13: f64, r0_ub: f64, k_ub: f64) -> (f64, f64) {
        let (e_angle, _) = self.harmonic_angle(theta, theta0, k_theta);
        let dr = r13 - r0_ub;
        let e_ub = 0.5 * k_ub * dr * dr;
        let f_ub = -k_ub * dr;
        (e_angle + e_ub, f_ub)
    }

    /// Compute angle θ between vectors (i→j) and (k→j).
    pub fn compute_angle(&self, ri: [f64; 3], rj: [f64; 3], rk: [f64; 3]) -> f64 {
        let a = sub3(ri, rj);
        let b = sub3(rk, rj);
        let la = len3(a);
        let lb = len3(b);
        if la < 1e-15 || lb < 1e-15 { return 0.0; }
        let cos_theta = (dot3(a, b) / (la * lb)).clamp(-1.0, 1.0);
        cos_theta.acos()
    }

    /// Compute angle forces for a list of angle triples.
    ///
    /// `triples` — list of (i, j, k) atom indices.
    /// `params`  — list of (theta0, k) per angle.
    pub fn compute_harmonic_angle_forces(
        &self,
        positions: &[[f64; 3]],
        triples: &[(usize, usize, usize)],
        params: &[(f64, f64)],
    ) -> (Vec<[f64; 3]>, f64) {
        let mut forces = vec![[0.0f64; 3]; positions.len()];
        let mut energy = 0.0;
        for (&(i, j, k), &(theta0, kth)) in triples.iter().zip(params.iter()) {
            let theta = self.compute_angle(positions[i], positions[j], positions[k]);
            let (e, tau) = self.harmonic_angle(theta, theta0, kth);
            energy += e;
            // Distribute torque as forces (simplified: equal on i and k)
            let a = sub3(positions[i], positions[j]);
            let b = sub3(positions[k], positions[j]);
            let la = len3(a);
            let lb = len3(b);
            if la < 1e-15 || lb < 1e-15 { continue; }
            let f_i = scale3(tau / la, [a[0]/la, a[1]/la, a[2]/la]);
            let f_k = scale3(tau / lb, [b[0]/lb, b[1]/lb, b[2]/lb]);
            forces[i] = add3(forces[i], f_i);
            forces[k] = add3(forces[k], f_k);
            let fj = [-(f_i[0]+f_k[0]), -(f_i[1]+f_k[1]), -(f_i[2]+f_k[2])];
            forces[j] = add3(forces[j], fj);
        }
        (forces, energy)
    }
}

impl Default for AngleForceKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── DihedralForceKernel ───────────────────────────────────────────────────────

/// Kernel for four-body dihedral (torsion) force computations.
#[derive(Debug, Clone)]
pub struct DihedralForceKernel;

impl DihedralForceKernel {
    /// Create a new dihedral force kernel.
    pub fn new() -> Self {
        Self
    }

    /// Cosine dihedral: U = k * (1 + cos(n*φ - δ)).
    pub fn cosine_dihedral(&self, phi: f64, n: f64, delta: f64, k: f64) -> (f64, f64) {
        let u = k * (1.0 + (n * phi - delta).cos());
        let du = -k * n * (n * phi - delta).sin();
        (u, du)
    }

    /// Ryckaert-Bellemans dihedral: U = sum_n C_n * cos^n(φ).
    /// `coeffs` — polynomial coefficients C_0, C_1, ..., C_5.
    pub fn ryckaert_bellemans(&self, phi: f64, coeffs: &[f64]) -> (f64, f64) {
        let cos_phi = phi.cos();
        let mut u = 0.0;
        let mut du_dphi = 0.0;
        let mut cos_pow = 1.0;
        for (n, &c) in coeffs.iter().enumerate() {
            u += c * cos_pow;
            if n > 0 {
                du_dphi -= c * n as f64 * cos_pow / cos_phi * phi.sin();
            }
            cos_pow *= cos_phi;
        }
        (u, du_dphi)
    }

    /// Compute dihedral angle φ for atoms i-j-k-l.
    pub fn compute_dihedral(
        &self,
        ri: [f64; 3],
        rj: [f64; 3],
        rk: [f64; 3],
        rl: [f64; 3],
    ) -> f64 {
        let b1 = sub3(rj, ri);
        let b2 = sub3(rk, rj);
        let b3 = sub3(rl, rk);
        let n1 = cross3(b1, b2);
        let n2 = cross3(b2, b3);
        let len_n1 = len3(n1);
        let len_n2 = len3(n2);
        if len_n1 < 1e-15 || len_n2 < 1e-15 { return 0.0; }
        let cos_phi = (dot3(n1, n2) / (len_n1 * len_n2)).clamp(-1.0, 1.0);
        let sign = dot3(cross3(n1, n2), b2);
        let phi = cos_phi.acos();
        if sign < 0.0 { -phi } else { phi }
    }
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
}

impl Default for DihedralForceKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── NeighborListKernel ────────────────────────────────────────────────────────

/// A cell list for O(N) neighbour searches.
#[derive(Debug, Clone)]
pub struct CellList {
    /// Number of cells in each direction.
    pub n_cells: [usize; 3],
    /// Cell spacing.
    pub cell_size: [f64; 3],
    /// Cell contents: each cell holds a list of atom indices.
    pub cells: Vec<Vec<usize>>,
}

impl CellList {
    /// Map atom position to cell index.
    pub fn cell_index(&self, pos: [f64; 3], box_size: [f64; 3]) -> [usize; 3] {
        let mut idx = [0usize; 3];
        for k in 0..3 {
            let fi = (pos[k] / box_size[k] * self.n_cells[k] as f64).floor() as i64;
            idx[k] = fi.rem_euclid(self.n_cells[k] as i64) as usize;
        }
        idx
    }
}

/// Kernel for building and querying neighbour lists.
#[derive(Debug, Clone)]
pub struct NeighborListKernel {
    /// Box size.
    pub box_size: [f64; 3],
    /// Use periodic boundaries.
    pub periodic: bool,
}

impl NeighborListKernel {
    /// Create a new neighbour list kernel.
    pub fn new(box_size: [f64; 3], periodic: bool) -> Self {
        Self { box_size, periodic }
    }

    /// Build a cell list for the given positions and cutoff.
    pub fn build_cell_list(&self, positions: &[[f64; 3]], cutoff: f64) -> CellList {
        let n_cells = [
            (self.box_size[0] / cutoff).floor().max(1.0) as usize,
            (self.box_size[1] / cutoff).floor().max(1.0) as usize,
            (self.box_size[2] / cutoff).floor().max(1.0) as usize,
        ];
        let cell_size = [
            self.box_size[0] / n_cells[0] as f64,
            self.box_size[1] / n_cells[1] as f64,
            self.box_size[2] / n_cells[2] as f64,
        ];
        let total_cells = n_cells[0] * n_cells[1] * n_cells[2];
        let mut cells = vec![Vec::new(); total_cells];

        let cl = CellList { n_cells, cell_size, cells: Vec::new() };
        for (i, pos) in positions.iter().enumerate() {
            let idx = cl.cell_index(*pos, self.box_size);
            let flat = idx[0] + n_cells[0] * (idx[1] + n_cells[1] * idx[2]);
            cells[flat].push(i);
        }
        CellList { n_cells, cell_size, cells }
    }

    /// Build a Verlet neighbour list from a cell list.
    ///
    /// Returns `neighbours[i]` = list of j > i within cutoff.
    pub fn build_verlet_list(
        &self,
        positions: &[[f64; 3]],
        cutoff: f64,
        cell_list: &CellList,
    ) -> Vec<Vec<usize>> {
        let n = positions.len();
        let mut neighbours = vec![Vec::new(); n];
        let rc2 = cutoff * cutoff;
        let nc = cell_list.n_cells;

        for i in 0..n {
            let ci = cell_list.cell_index(positions[i], self.box_size);
            // Search neighbouring cells
            for dz in -1i64..=1 {
                for dy in -1i64..=1 {
                    for dx in -1i64..=1 {
                        let nx = (ci[0] as i64 + dx).rem_euclid(nc[0] as i64) as usize;
                        let ny = (ci[1] as i64 + dy).rem_euclid(nc[1] as i64) as usize;
                        let nz = (ci[2] as i64 + dz).rem_euclid(nc[2] as i64) as usize;
                        let flat = nx + nc[0] * (ny + nc[1] * nz);
                        for &j in &cell_list.cells[flat] {
                            if j <= i { continue; }
                            let mut dr = sub3(positions[i], positions[j]);
                            if self.periodic { dr = minimum_image(dr, self.box_size); }
                            if dot3(dr, dr) <= rc2 {
                                neighbours[i].push(j);
                            }
                        }
                    }
                }
            }
        }
        neighbours
    }
}

// ── ThermostatKernel ──────────────────────────────────────────────────────────

/// Kernel for thermostat algorithms.
#[derive(Debug, Clone)]
pub struct ThermostatKernel {
    /// Target temperature.
    pub t_target: f64,
    /// Boltzmann constant (1.0 in reduced units).
    pub k_b: f64,
}

impl ThermostatKernel {
    /// Create a new thermostat kernel.
    pub fn new(t_target: f64, k_b: f64) -> Self {
        Self { t_target, k_b }
    }

    /// Current kinetic temperature from velocities.
    ///
    /// `masses` — per-atom masses.
    /// `velocities` — per-atom velocities.
    /// `n_dof` — number of degrees of freedom.
    pub fn current_temperature(
        &self,
        masses: &[f64],
        velocities: &[[f64; 3]],
        n_dof: usize,
    ) -> f64 {
        let ke: f64 = masses.iter().zip(velocities.iter())
            .map(|(&m, v)| 0.5 * m * dot3(*v, *v))
            .sum();
        2.0 * ke / (n_dof as f64 * self.k_b)
    }

    /// Velocity rescaling thermostat.
    ///
    /// Multiplies all velocities by sqrt(T_target / T_current).
    pub fn velocity_rescale(
        &self,
        velocities: &mut Vec<[f64; 3]>,
        t_current: f64,
    ) {
        if t_current < 1e-15 { return; }
        let scale = (self.t_target / t_current).sqrt();
        for v in velocities.iter_mut() {
            v[0] *= scale;
            v[1] *= scale;
            v[2] *= scale;
        }
    }

    /// Nosé-Hoover chain (NHC) propagation step for a single thermostat variable.
    ///
    /// `xi`     — thermostat position (heat bath variable).
    /// `v_xi`   — thermostat velocity.
    /// `g`      — thermostat "force" = (ke - n_dof/2 * kT).
    /// `q`      — thermostat mass.
    /// `dt`     — time step.
    ///
    /// Returns updated (xi, v_xi, scale_factor_for_velocities).
    pub fn nhc_propagate(
        &self,
        xi: f64,
        v_xi: f64,
        g: f64,
        q: f64,
        dt: f64,
    ) -> (f64, f64, f64) {
        let v_xi_new = v_xi + 0.5 * dt * g / q;
        let scale = (-0.5 * dt * v_xi_new).exp();
        let xi_new = xi + dt * v_xi_new;
        let v_xi_final = v_xi_new + 0.5 * dt * g / q;
        (xi_new, v_xi_final, scale)
    }

    /// Berendsen weak-coupling thermostat scale factor.
    ///
    /// λ = sqrt(1 + dt/tau * (T_target/T_current - 1)).
    pub fn berendsen_scale(&self, t_current: f64, tau: f64, dt: f64) -> f64 {
        if t_current < 1e-15 { return 1.0; }
        (1.0 + dt / tau * (self.t_target / t_current - 1.0)).sqrt().abs()
    }
}

// ── BarostatKernel ────────────────────────────────────────────────────────────

/// Kernel for barostat (pressure coupling) algorithms.
#[derive(Debug, Clone)]
pub struct BarostatKernel {
    /// Target pressure.
    pub p_target: f64,
    /// Isothermal compressibility β_T (1/GPa for water ≈ 4.57e-10).
    pub beta_t: f64,
}

impl BarostatKernel {
    /// Create a new barostat kernel.
    pub fn new(p_target: f64, beta_t: f64) -> Self {
        Self { p_target, beta_t }
    }

    /// Current pressure from the virial theorem.
    ///
    /// P = (N * k_B * T + W) / V
    pub fn current_pressure(
        &self,
        n_atoms: usize,
        k_b: f64,
        temperature: f64,
        virial: f64,
        volume: f64,
    ) -> f64 {
        (n_atoms as f64 * k_b * temperature + virial) / volume
    }

    /// Berendsen isotropic barostat: scale box and positions.
    ///
    /// Returns the scale factor μ.
    pub fn berendsen_scale(
        &self,
        p_current: f64,
        tau_p: f64,
        dt: f64,
    ) -> f64 {
        let mu3 = 1.0 - self.beta_t * dt / tau_p * (self.p_target - p_current);
        mu3.cbrt()
    }

    /// MTTK (Martyna-Tobias-Tuckerman-Klein) barostat velocity update.
    ///
    /// `v_box` — barostat velocity (epsilon_dot).
    /// `g_box` — barostat force = (P_inst - P_ref) * V.
    /// `w_box` — barostat mass.
    /// Returns updated v_box.
    pub fn mttk_barostat_step(
        &self,
        v_box: f64,
        g_box: f64,
        w_box: f64,
        dt: f64,
    ) -> f64 {
        v_box + 0.5 * dt * g_box / w_box
    }

    /// Scale positions and box for an isotropic volume change by μ.
    pub fn scale_positions_and_box(
        &self,
        positions: &mut Vec<[f64; 3]>,
        box_size: &mut [f64; 3],
        mu: f64,
    ) {
        for pos in positions.iter_mut() {
            pos[0] *= mu;
            pos[1] *= mu;
            pos[2] *= mu;
        }
        box_size[0] *= mu;
        box_size[1] *= mu;
        box_size[2] *= mu;
    }
}

// ── PmeKernel ─────────────────────────────────────────────────────────────────

/// Particle Mesh Ewald (PME) kernel for long-range electrostatics.
#[derive(Debug, Clone)]
pub struct PmeKernel {
    /// PME grid dimensions \[nx, ny, nz\].
    pub grid_dims: [usize; 3],
    /// Ewald splitting parameter β.
    pub beta: f64,
    /// B-spline interpolation order.
    pub spline_order: usize,
    /// Box size.
    pub box_size: [f64; 3],
}

impl PmeKernel {
    /// Create a new PME kernel.
    pub fn new(grid_dims: [usize; 3], beta: f64, spline_order: usize, box_size: [f64; 3]) -> Self {
        Self { grid_dims, beta, spline_order, box_size }
    }

    /// Spread charges onto the PME grid using nearest-grid-point assignment.
    ///
    /// `positions` — atom positions.
    /// `charges`   — atom charges.
    /// Returns a flat charge-density grid \[nx * ny * nz\].
    pub fn pme_charge_spreading(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
    ) -> Vec<f64> {
        let [nx, ny, nz] = self.grid_dims;
        let mut grid = vec![0.0f64; nx * ny * nz];
        for (pos, &q) in positions.iter().zip(charges.iter()) {
            let ix = ((pos[0] / self.box_size[0]) * nx as f64).floor() as usize % nx;
            let iy = ((pos[1] / self.box_size[1]) * ny as f64).floor() as usize % ny;
            let iz = ((pos[2] / self.box_size[2]) * nz as f64).floor() as usize % nz;
            grid[iz * ny * nx + iy * nx + ix] += q;
        }
        grid
    }

    /// Compute PME reciprocal energy from a charge grid.
    ///
    /// Uses a simplified Green's function G(k) = exp(-k²/(4β²)) / k².
    pub fn pme_reciprocal_sum(&self, grid: &[f64]) -> f64 {
        let [nx, ny, nz] = self.grid_dims;
        let vol = self.box_size[0] * self.box_size[1] * self.box_size[2];
        let prefac = 4.0 * PI / vol;
        let mut energy = 0.0;

        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let q = grid[iz * ny * nx + iy * nx + ix];
                    if q.abs() < 1e-20 { continue; }
                    let kx = 2.0 * PI * (ix as i64 - if ix > nx/2 { nx as i64 } else { 0 }) as f64
                        / self.box_size[0];
                    let ky = 2.0 * PI * (iy as i64 - if iy > ny/2 { ny as i64 } else { 0 }) as f64
                        / self.box_size[1];
                    let kz = 2.0 * PI * (iz as i64 - if iz > nz/2 { nz as i64 } else { 0 }) as f64
                        / self.box_size[2];
                    let k2 = kx*kx + ky*ky + kz*kz;
                    if k2 < 1e-20 { continue; }
                    let green = prefac * (-k2 / (4.0 * self.beta * self.beta)).exp() / k2;
                    energy += green * q * q;
                }
            }
        }
        0.5 * energy
    }

    /// B-spline basis function of order n at x ∈ \[0, 1\].
    pub fn bspline(&self, x: f64, order: usize) -> f64 {
        if order == 1 {
            if x >= 0.0 && x < 1.0 { 1.0 } else { 0.0 }
        } else {
            let o = order as f64;
            x / (o - 1.0) * self.bspline(x, order - 1)
                + (o - x) / (o - 1.0) * self.bspline(x - 1.0, order - 1)
        }
    }

    /// PME self-energy correction.
    pub fn self_energy(&self, charges: &[f64]) -> f64 {
        let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
        -self.beta / PI.sqrt() * sum_q2
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod md_tests {
    use super::*;

    fn make_positions_2() -> Vec<[f64; 3]> {
        vec![[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]]
    }

    #[test]
    fn test_md_config_volume() {
        let cfg = MdKernelConfig::new_lj_gas(100, [10.0, 10.0, 10.0]);
        assert!((cfg.volume() - 1000.0).abs() < 1e-10);
    }

    #[test]
    fn test_md_config_number_density() {
        let cfg = MdKernelConfig::new_lj_gas(100, [10.0, 10.0, 10.0]);
        assert!((cfg.number_density() - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_md_config_list_cutoff() {
        let cfg = MdKernelConfig::new_lj_gas(100, [10.0, 10.0, 10.0]);
        assert!((cfg.list_cutoff() - 2.8).abs() < 1e-10);
    }

    #[test]
    fn test_lj_potential_minimum() {
        let kern = LjForceKernel::new(
            LjParams::argon(), 3.0, [10.0; 3], false
        );
        // Minimum at r = 2^(1/6) * sigma ≈ 1.122
        let r_min = 2.0f64.powf(1.0 / 6.0);
        let u = kern.lj_potential(r_min);
        assert!((u + 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_lj_forces_conservation() {
        let kern = LjForceKernel::new(
            LjParams::argon(), 3.0, [20.0; 3], false
        );
        let positions = make_positions_2();
        let (forces, _) = kern.compute_lj_forces(&positions);
        // Newton's 3rd law: forces should sum to zero
        for k in 0..3 {
            let sum = forces[0][k] + forces[1][k];
            assert!(sum.abs() < 1e-10, "force not conserved along axis {k}");
        }
    }

    #[test]
    fn test_lj_forces_repulsive_at_close_range() {
        let kern = LjForceKernel::new(LjParams::argon(), 3.0, [20.0; 3], false);
        let positions = vec![[0.0, 0.0, 0.0], [0.8, 0.0, 0.0]];
        let (forces, _) = kern.compute_lj_forces(&positions);
        // Atom 0 pushed in -x
        assert!(forces[0][0] < 0.0);
    }

    #[test]
    fn test_ewald_self_energy() {
        let ew = EwaldKernel::new(0.5, 5.0, 3, [10.0; 3]);
        let charges = vec![1.0, -1.0];
        let se = ew.self_energy(&charges);
        assert!(se < 0.0);
    }

    #[test]
    fn test_ewald_real_space_neutral() {
        let ew = EwaldKernel::new(0.5, 5.0, 2, [10.0; 3]);
        let pos = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let e = ew.ewald_real_space(&pos, &charges);
        assert!(e < 0.0); // opposites attract
    }

    #[test]
    fn test_harmonic_bond_zero_at_r0() {
        let bk = BondForceKernel::new();
        let (e, f) = bk.harmonic_bond(1.5, 1.5, 1000.0);
        assert!(e.abs() < 1e-12);
        assert!(f.abs() < 1e-12);
    }

    #[test]
    fn test_morse_bond_minimum() {
        let bk = BondForceKernel::new();
        let (e, f) = bk.morse_bond(1.0, 1.0, 5.0, 2.0);
        assert!(e.abs() < 1e-12);
        assert!(f.abs() < 1e-12);
    }

    #[test]
    fn test_fene_bond_restoring() {
        let bk = BondForceKernel::new();
        let (e, f) = bk.fene_bond(1.5, 1.0, 10.0, 1.0);
        assert!(e > 0.0);
        assert!(f < 0.0); // restoring toward r0
    }

    #[test]
    fn test_harmonic_angle_zero_at_theta0() {
        let ak = AngleForceKernel::new();
        let (e, f) = ak.harmonic_angle(1.0472, 1.0472, 500.0);
        assert!(e.abs() < 1e-10);
        assert!(f.abs() < 1e-10);
    }

    #[test]
    fn test_compute_angle_linear() {
        let ak = AngleForceKernel::new();
        let theta = ak.compute_angle([1.0,0.0,0.0], [0.0,0.0,0.0], [-1.0,0.0,0.0]);
        assert!((theta - PI).abs() < 1e-10);
    }

    #[test]
    fn test_compute_angle_right() {
        let ak = AngleForceKernel::new();
        let theta = ak.compute_angle([1.0,0.0,0.0], [0.0,0.0,0.0], [0.0,1.0,0.0]);
        assert!((theta - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_cosine_dihedral_energy() {
        let dk = DihedralForceKernel::new();
        let (e, _) = dk.cosine_dihedral(0.0, 1.0, 0.0, 1.0);
        assert!((e - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_ryckaert_bellemans_constant() {
        let dk = DihedralForceKernel::new();
        let (e, du) = dk.ryckaert_bellemans(PI / 2.0, &[3.0, 0.0, 0.0]);
        assert!((e - 3.0).abs() < 1e-10);
        assert!(du.abs() < 1e-10); // derivative at pi/2 is zero for constant
    }

    #[test]
    fn test_build_cell_list() {
        let kern = NeighborListKernel::new([10.0; 3], true);
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64, 0.0, 0.0]).collect();
        let cl = kern.build_cell_list(&positions, 2.5);
        assert!(!cl.cells.is_empty());
    }

    #[test]
    fn test_verlet_list_symmetry() {
        let kern = NeighborListKernel::new([20.0; 3], false);
        let positions = vec![[0.0,0.0,0.0],[2.0,0.0,0.0],[10.0,0.0,0.0]];
        let cl = kern.build_cell_list(&positions, 3.0);
        let nl = kern.build_verlet_list(&positions, 3.0, &cl);
        // Atoms 0 and 1 should be neighbours
        assert!(nl[0].contains(&1), "0 should see 1, nl[0]={:?}", nl[0]);
    }

    #[test]
    fn test_velocity_rescale() {
        let thermo = ThermostatKernel::new(300.0, 1.0);
        let mut velocities = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let t_cur = thermo.current_temperature(&masses, &velocities, 6);
        thermo.velocity_rescale(&mut velocities, t_cur);
        let t_new = thermo.current_temperature(&masses, &velocities, 6);
        assert!((t_new - thermo.t_target).abs() < 1e-6);
    }

    #[test]
    fn test_nhc_propagate() {
        let thermo = ThermostatKernel::new(300.0, 1.0);
        let (xi, v_xi, scale) = thermo.nhc_propagate(0.0, 0.1, 0.5, 10.0, 0.001);
        assert!(xi.is_finite());
        assert!(v_xi.is_finite());
        assert!(scale > 0.0);
    }

    #[test]
    fn test_berendsen_scale_thermostat() {
        let thermo = ThermostatKernel::new(300.0, 1.0);
        let scale = thermo.berendsen_scale(290.0, 0.1, 0.001);
        assert!(scale > 1.0); // need to heat up
    }

    #[test]
    fn test_barostat_berendsen_scale() {
        let baro = BarostatKernel::new(1.0, 4.57e-10);
        let mu = baro.berendsen_scale(1.1, 0.5, 0.001);
        assert!(mu.is_finite());
    }

    #[test]
    fn test_current_pressure() {
        let baro = BarostatKernel::new(1.0, 4.57e-10);
        let p = baro.current_pressure(100, 1.0, 1.0, 50.0, 1000.0);
        assert!((p - 0.15).abs() < 1e-10);
    }

    #[test]
    fn test_pme_charge_spreading() {
        let pme = PmeKernel::new([4, 4, 4], 0.5, 4, [10.0; 3]);
        let pos = vec![[0.0, 0.0, 0.0]];
        let charges = vec![1.0];
        let grid = pme.pme_charge_spreading(&pos, &charges);
        let total: f64 = grid.iter().sum();
        assert!((total - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_pme_reciprocal_sum_positive() {
        let pme = PmeKernel::new([4, 4, 4], 0.5, 4, [10.0; 3]);
        let pos = vec![[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let grid = pme.pme_charge_spreading(&pos, &charges);
        let e = pme.pme_reciprocal_sum(&grid);
        assert!(e.is_finite());
    }

    #[test]
    fn test_pme_self_energy() {
        let pme = PmeKernel::new([8, 8, 8], 0.3, 4, [10.0; 3]);
        let charges = vec![1.0, -1.0];
        let se = pme.self_energy(&charges);
        assert!(se < 0.0);
    }

    #[test]
    fn test_erfc_approx() {
        // erfc(0) = 1
        assert!((erfc_approx(0.0) - 1.0).abs() < 1e-5);
        // erfc(large) ≈ 0
        assert!(erfc_approx(5.0) < 1e-10);
    }

    #[test]
    fn test_minimum_image_wraps() {
        let dr = [6.0, 0.0, 0.0];
        let bx = [10.0, 10.0, 10.0];
        let mdr = minimum_image(dr, bx);
        assert!((mdr[0] + 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_lj_virial() {
        let kern = LjForceKernel::new(LjParams::argon(), 3.0, [20.0; 3], false);
        let positions = vec![[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let virial = kern.compute_virial(&positions);
        assert!(virial.is_finite());
    }
}
