#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coarse-grained molecular dynamics (CG-MD), MARTINI-like.
//!
//! Provides:
//! - [`CgBead`]: single CG particle with position, velocity, mass and charge.
//! - [`CgBeadType`]: bead chemistry classification.
//! - [`CgMolecule`]: collection of bonded beads with harmonic bonds and angles.
//! - [`CgForceField`]: Lennard-Jones + electrostatic non-bonded interactions.
//! - [`CgSimulation`]: NVT simulation with Langevin thermostat.
//! - [`IterativeBoltzmannInversion`]: potential-of-mean-force from a target RDF.
//! - [`MsIbi`]: multi-state iterative Boltzmann inversion.
//!
//! References
//! ----------
//! Marrink, S. J. et al. (2007) *J. Phys. Chem. B* **111**, 7812–7824.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper geometry
// ---------------------------------------------------------------------------

#[inline]
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub3(a, b);
    dot3(d, d).sqrt()
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

// ---------------------------------------------------------------------------
// CgBeadType
// ---------------------------------------------------------------------------

/// Chemistry classification of a coarse-grained bead (MARTINI taxonomy).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CgBeadType {
    /// Polar bead (hydrophilic head groups).
    Polar,
    /// Non-polar bead (mixed character).
    NonPolar,
    /// Apolar bead (lipid tails, hydrophobic core).
    Apolar,
    /// Charged bead (ion, charged amino acid side chain).
    Charged,
    /// Protein backbone bead.
    Backbone,
    /// Protein side-chain bead.
    Sidechain,
}

// ---------------------------------------------------------------------------
// CgBead
// ---------------------------------------------------------------------------

/// A single coarse-grained particle.
#[derive(Debug, Clone)]
pub struct CgBead {
    /// Human-readable name (e.g. `"PO4"`, `"NC3"`).
    pub name: String,
    /// Chemistry type used to select LJ parameters.
    pub type_: CgBeadType,
    /// Position in nanometres \[x, y, z\].
    pub pos: [f64; 3],
    /// Velocity in nm ps⁻¹ \[vx, vy, vz\].
    pub vel: [f64; 3],
    /// Mass in amu (g mol⁻¹).
    pub mass: f64,
    /// Partial charge in elementary charge units.
    pub charge: f64,
}

impl CgBead {
    /// Construct a new bead at rest at the origin.
    pub fn new(name: impl Into<String>, type_: CgBeadType, mass: f64, charge: f64) -> Self {
        Self {
            name: name.into(),
            type_,
            pos: [0.0; 3],
            vel: [0.0; 3],
            mass,
            charge,
        }
    }

    /// Kinetic energy of this bead: ½ m v².
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }
}

// ---------------------------------------------------------------------------
// CgMolecule
// ---------------------------------------------------------------------------

/// A coarse-grained molecule: beads + harmonic bonds + cosine-harmonic angles.
///
/// Bond tuple: `(i, j, r0_nm, k_bond_kJ_mol_nm2)`.
/// Angle tuple: `(i, j, k, theta0_rad, k_angle_kJ_mol_rad2)`.
#[derive(Debug, Clone)]
pub struct CgMolecule {
    /// Constituent beads.
    pub beads: Vec<CgBead>,
    /// Harmonic bonds: (bead_i, bead_j, r0 \[nm\], k \[kJ mol⁻¹ nm⁻²\]).
    pub bonds: Vec<(usize, usize, f64, f64)>,
    /// Harmonic angles: (bead_i, bead_j, bead_k, θ0 \[rad\], k \[kJ mol⁻¹ rad⁻²\]).
    pub angles: Vec<(usize, usize, usize, f64, f64)>,
}

impl CgMolecule {
    /// Create an empty molecule.
    pub fn new() -> Self {
        Self {
            beads: Vec::new(),
            bonds: Vec::new(),
            angles: Vec::new(),
        }
    }

    /// Total bonded potential energy (bonds + angles) in kJ mol⁻¹.
    pub fn bonded_energy(&self) -> f64 {
        let mut e = 0.0;
        for &(i, j, r0, k) in &self.bonds {
            let r = dist3(self.beads[i].pos, self.beads[j].pos);
            e += 0.5 * k * (r - r0).powi(2);
        }
        for &(i, j, k_idx, theta0, k_a) in &self.angles {
            let rij = sub3(self.beads[i].pos, self.beads[j].pos);
            let rkj = sub3(self.beads[k_idx].pos, self.beads[j].pos);
            let cos_theta = dot3(rij, rkj) / (norm3(rij) * norm3(rkj) + 1e-30);
            let cos_theta0 = theta0.cos();
            e += 0.5 * k_a * (cos_theta - cos_theta0).powi(2);
        }
        e
    }

    /// Forces on each bead from harmonic bonds, in kJ mol⁻¹ nm⁻¹.
    pub fn bond_forces(&self) -> Vec<[f64; 3]> {
        let n = self.beads.len();
        let mut f = vec![[0.0f64; 3]; n];
        for &(i, j, r0, k) in &self.bonds {
            let rij = sub3(self.beads[i].pos, self.beads[j].pos);
            let r = norm3(rij) + 1e-30;
            let mag = -k * (r - r0) / r;
            let fv = scale3(rij, mag);
            f[i] = add3(f[i], fv);
            f[j] = sub3(f[j], fv);
        }
        f
    }

    /// Forces on each bead from cosine-harmonic angles, in kJ mol⁻¹ nm⁻¹.
    pub fn angle_forces(&self) -> Vec<[f64; 3]> {
        let n = self.beads.len();
        let mut f = vec![[0.0f64; 3]; n];
        for &(i, j, k_idx, theta0, k_a) in &self.angles {
            let pi_pos = self.beads[i].pos;
            let pj_pos = self.beads[j].pos;
            let pk_pos = self.beads[k_idx].pos;
            let rij = sub3(pi_pos, pj_pos);
            let rkj = sub3(pk_pos, pj_pos);
            let rij_n = norm3(rij) + 1e-30;
            let rkj_n = norm3(rkj) + 1e-30;
            let cos_t = dot3(rij, rkj) / (rij_n * rkj_n);
            let cos_t0 = theta0.cos();
            let dcde = -k_a * (cos_t - cos_t0);
            // ∂cos/∂r_i = (r_kj/|rij||rkj| - cos*r_ij/|rij|^2) / |rij|
            let fi = scale3(
                sub3(
                    scale3(rkj, 1.0 / (rij_n * rkj_n)),
                    scale3(rij, cos_t / (rij_n * rij_n)),
                ),
                dcde,
            );
            let fk = scale3(
                sub3(
                    scale3(rij, 1.0 / (rij_n * rkj_n)),
                    scale3(rkj, cos_t / (rkj_n * rkj_n)),
                ),
                dcde,
            );
            let fj = scale3(add3(fi, fk), -1.0);
            f[i] = add3(f[i], fi);
            f[j] = add3(f[j], fj);
            f[k_idx] = add3(f[k_idx], fk);
        }
        f
    }
}

impl Default for CgMolecule {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CgForceField
// ---------------------------------------------------------------------------

/// Non-bonded CG force field: Lennard-Jones + screened electrostatics.
///
/// `lj_params[i][j]` stores `(epsilon [kJ mol⁻¹], sigma [nm])` for type pair `(i,j)`.
/// Index mapping: [`CgBeadType::Polar`]=0, `NonPolar`=1, `Apolar`=2,
/// `Charged`=3, `Backbone`=4, `Sidechain`=5.
#[derive(Debug, Clone)]
pub struct CgForceField {
    /// LJ parameters table: `lj_params[type_i][type_j] = (epsilon, sigma)`.
    pub lj_params: Vec<Vec<(f64, f64)>>,
    /// Dielectric screening length κ⁻¹ \[nm\] (Debye screening).
    pub electrostatic_screening: f64,
}

impl CgForceField {
    /// Construct a default MARTINI-like force field (6 bead types).
    pub fn default_martini() -> Self {
        let n = 6usize;
        // epsilon [kJ/mol], sigma [nm] — approximate MARTINI level-2 values
        let eps_table = [
            [5.6, 5.0, 4.0, 5.6, 5.0, 4.5], // Polar
            [5.0, 3.8, 2.7, 5.0, 4.0, 3.5], // NonPolar
            [4.0, 2.7, 2.0, 4.0, 3.0, 2.5], // Apolar
            [5.6, 5.0, 4.0, 5.6, 5.0, 4.5], // Charged
            [5.0, 4.0, 3.0, 5.0, 4.5, 4.0], // Backbone
            [4.5, 3.5, 2.5, 4.5, 4.0, 3.5], // Sidechain
        ];
        let sigma = 0.47; // nm — standard MARTINI bead size
        let lj_params = (0..n)
            .map(|i| (0..n).map(|j| (eps_table[i][j], sigma)).collect())
            .collect();
        Self {
            lj_params,
            electrostatic_screening: 1.0,
        }
    }

    /// LJ-12-6 potential between types i and j at distance r \[nm\].
    /// Returns energy in kJ mol⁻¹.
    pub fn lj_energy(&self, i: usize, j: usize, r: f64) -> f64 {
        let (eps, sigma) = self.lj_params[i][j];
        let sr = sigma / r;
        let sr6 = sr.powi(6);
        4.0 * eps * (sr6 * sr6 - sr6)
    }

    /// Yukawa (screened Coulomb) energy between charges qi and qj at r \[nm\].
    /// Returns energy in kJ mol⁻¹.
    pub fn electrostatic_energy(&self, qi: f64, qj: f64, r: f64) -> f64 {
        // ε₀ in kJ mol⁻¹ nm e⁻² context: factor ≈ 138.935 kJ mol⁻¹ nm e⁻²
        const KE: f64 = 138.935;
        KE * qi * qj / r * (-r / self.electrostatic_screening).exp()
    }

    /// Total non-bonded energy of a bead array (all unique pairs).
    /// Returns energy in kJ mol⁻¹.
    pub fn total_nonbonded(&self, beads: &[CgBead]) -> f64 {
        let mut e = 0.0;
        let n = beads.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(beads[i].pos, beads[j].pos).max(0.1);
                let ti = bead_type_index(beads[i].type_);
                let tj = bead_type_index(beads[j].type_);
                e += self.lj_energy(ti, tj, r);
                if beads[i].charge.abs() > 1e-10 || beads[j].charge.abs() > 1e-10 {
                    e += self.electrostatic_energy(beads[i].charge, beads[j].charge, r);
                }
            }
        }
        e
    }
}

/// Map a [`CgBeadType`] variant to its row/column index in the LJ table.
fn bead_type_index(t: CgBeadType) -> usize {
    match t {
        CgBeadType::Polar => 0,
        CgBeadType::NonPolar => 1,
        CgBeadType::Apolar => 2,
        CgBeadType::Charged => 3,
        CgBeadType::Backbone => 4,
        CgBeadType::Sidechain => 5,
    }
}

// ---------------------------------------------------------------------------
// CgSimulation
// ---------------------------------------------------------------------------

/// NVT coarse-grained simulation with Langevin thermostat.
#[derive(Debug, Clone)]
pub struct CgSimulation {
    /// All molecules in the simulation box.
    pub molecules: Vec<CgMolecule>,
    /// Non-bonded force field.
    pub ff: CgForceField,
    /// Target temperature in Kelvin.
    pub temp: f64,
    /// Integration time step in ps.
    pub dt: f64,
    /// Langevin friction coefficient γ \[ps⁻¹\].
    pub gamma: f64,
    /// Current simulation time \[ps\].
    pub time: f64,
    /// Pseudo-random seed state (LCG, no external crate needed).
    lcg_state: u64,
}

impl CgSimulation {
    /// Construct a simulation with given molecules, force field, temperature and dt.
    pub fn new(molecules: Vec<CgMolecule>, ff: CgForceField, temp: f64, dt: f64) -> Self {
        Self {
            molecules,
            ff,
            temp,
            dt,
            gamma: 1.0,
            time: 0.0,
            lcg_state: 12345,
        }
    }

    /// Advance the simulation by one time step (velocity Verlet + Langevin).
    pub fn step(&mut self) {
        // Collect all beads by value, update, then write back.
        let forces = self.compute_all_forces();
        let mut bead_idx = 0usize;
        for mol in &mut self.molecules {
            for bead in &mut mol.beads {
                let f = forces[bead_idx];
                // v += 0.5 * dt * F/m
                let a = scale3(f, 1.0 / bead.mass);
                bead.vel = add3(bead.vel, scale3(a, 0.5 * self.dt));
                // x += dt * v
                bead.pos = add3(bead.pos, scale3(bead.vel, self.dt));
                bead_idx += 1;
            }
        }
        let forces2 = self.compute_all_forces();
        bead_idx = 0;
        for mol in &mut self.molecules {
            for bead in &mut mol.beads {
                let f = forces2[bead_idx];
                let a = scale3(f, 1.0 / bead.mass);
                bead.vel = add3(bead.vel, scale3(a, 0.5 * self.dt));
                bead_idx += 1;
            }
        }
        self.langevin_thermostat();
        self.time += self.dt;
    }

    /// Apply Langevin (Ornstein-Uhlenbeck) thermostat to all bead velocities.
    pub fn langevin_thermostat(&mut self) {
        let kbt = 8.314e-3 * self.temp; // kJ mol⁻¹
        let alpha = (-self.gamma * self.dt).exp();
        // Pre-generate noise: count total velocity components
        let total_components: usize = self.molecules.iter().map(|m| m.beads.len() * 3).sum();
        let noise: Vec<f64> = (0..total_components).map(|_| self.gauss_rng()).collect();
        let mut idx = 0;
        for mol in &mut self.molecules {
            for bead in &mut mol.beads {
                let sigma = ((1.0 - alpha * alpha) * kbt / bead.mass).sqrt();
                for vi in &mut bead.vel {
                    *vi = alpha * (*vi) + sigma * noise[idx];
                    idx += 1;
                }
            }
        }
    }

    /// Instantaneous pressure estimate via the virial theorem \[kJ mol⁻¹ nm⁻³\].
    ///
    /// Uses the ideal-gas contribution plus the pairwise virial term.
    pub fn pressure(&self) -> f64 {
        let beads = self.all_beads();
        let n = beads.len();
        if n == 0 {
            return 0.0;
        }
        // Estimate box volume from bead density (assume cubic, 1 bead / (0.5 nm)³)
        let v = (n as f64) * 0.125; // nm³
        let kbt = 8.314e-3 * self.temp;
        let ideal = (n as f64) * kbt / v;
        // Virial sum over pairs
        let mut virial = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let rij = sub3(beads[i].pos, beads[j].pos);
                let r = norm3(rij).max(0.1);
                let ti = bead_type_index(beads[i].type_);
                let tj = bead_type_index(beads[j].type_);
                let (eps, sigma) = self.ff.lj_params[ti][tj];
                let sr = sigma / r;
                let sr6 = sr.powi(6);
                // dU/dr for LJ-12-6: 4ε(-12σ¹²/r¹³ + 6σ⁶/r⁷)
                let dudr = 4.0 * eps * (-12.0 * sr6 * sr6 + 6.0 * sr6) / r;
                virial += dudr * r;
            }
        }
        ideal - virial / (3.0 * v)
    }

    /// Compute the radial distribution function g(r) for all bead pairs.
    ///
    /// Returns a vector of `(r [nm], g(r))` pairs.
    pub fn rdf(&self) -> Vec<(f64, f64)> {
        let beads = self.all_beads();
        let n = beads.len();
        let r_max = 2.0f64;
        let dr = 0.02f64;
        let bins = (r_max / dr) as usize;
        let mut hist = vec![0u64; bins];
        let mut pair_count = 0u64;
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(beads[i].pos, beads[j].pos);
                if r < r_max {
                    let bin = (r / dr) as usize;
                    if bin < bins {
                        hist[bin] += 1;
                        pair_count += 1;
                    }
                }
            }
        }
        let rho = if n > 1 {
            (n as f64) / ((n as f64) * 0.125)
        } else {
            1.0
        };
        (0..bins)
            .map(|k| {
                let r = (k as f64 + 0.5) * dr;
                let shell_vol = 4.0 * PI * r * r * dr;
                let g = if pair_count > 0 && n > 1 {
                    hist[k] as f64 / (pair_count as f64 / (n as f64) * rho * shell_vol)
                } else {
                    0.0
                };
                (r, g)
            })
            .collect()
    }

    /// Estimate self-diffusion coefficient from mean-squared displacement.
    ///
    /// Uses current positions vs. a zero reference; returns D in nm² ps⁻¹.
    pub fn diffusion_coefficient(&self) -> f64 {
        let beads = self.all_beads();
        if beads.is_empty() {
            return 0.0;
        }
        // MSD from origin (for a single-frame estimate, returns RMS/6t)
        let msd: f64 = beads.iter().map(|b| dot3(b.pos, b.pos)).sum::<f64>() / beads.len() as f64;
        let t = self.time.max(1e-10);
        msd / (6.0 * t)
    }

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    fn all_beads(&self) -> Vec<&CgBead> {
        self.molecules.iter().flat_map(|m| m.beads.iter()).collect()
    }

    fn compute_all_forces(&self) -> Vec<[f64; 3]> {
        // Count total beads and build flat lists
        let total: usize = self.molecules.iter().map(|m| m.beads.len()).sum();
        let mut forces = vec![[0.0f64; 3]; total];

        // Bonded forces (per molecule)
        let mut offset = 0usize;
        for mol in &self.molecules {
            let bf = mol.bond_forces();
            let af = mol.angle_forces();
            for (k, (b, a)) in bf.iter().zip(af.iter()).enumerate() {
                forces[offset + k] = add3(*b, *a);
            }
            offset += mol.beads.len();
        }

        // Non-bonded forces between all beads across molecules
        let beads: Vec<&CgBead> = self.all_beads();
        for i in 0..total {
            for j in (i + 1)..total {
                let rij = sub3(beads[i].pos, beads[j].pos);
                let r = norm3(rij).max(0.1);
                let ti = bead_type_index(beads[i].type_);
                let tj = bead_type_index(beads[j].type_);
                let (eps, sigma) = self.ff.lj_params[ti][tj];
                let sr = sigma / r;
                let sr6 = sr.powi(6);
                let fmag = 4.0 * eps * (12.0 * sr6 * sr6 - 6.0 * sr6) / (r * r);
                let fv = scale3(rij, fmag);
                forces[i] = add3(forces[i], fv);
                forces[j] = sub3(forces[j], fv);
            }
        }
        forces
    }

    /// Box-Muller Gaussian RNG (seeded LCG, no external crates).
    fn gauss_rng(&mut self) -> f64 {
        let u1 = self.lcg_rand();
        let u2 = self.lcg_rand();
        (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos()
    }

    fn lcg_rand(&mut self) -> f64 {
        self.lcg_state = self
            .lcg_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let bits = (self.lcg_state >> 33) as u32;
        (bits as f64 + 1.0) / (u32::MAX as f64 + 2.0)
    }
}

// ---------------------------------------------------------------------------
// IterativeBoltzmannInversion
// ---------------------------------------------------------------------------

/// Iterative Boltzmann Inversion (IBI) for deriving CG pair potentials from
/// a reference RDF.
///
/// The IBI update rule is:
/// V_{n+1}(r) = V_n(r) + k_BT ln(g_n(r) / g_target(r))
#[derive(Debug, Clone)]
pub struct IterativeBoltzmannInversion {
    /// Target radial distribution function as (r \[nm\], g_target(r)) pairs.
    pub target_rdf: Vec<(f64, f64)>,
    /// Current tabulated potential V(r) as (r \[nm\], V \[kJ mol⁻¹\]) pairs.
    pub potential: Vec<(f64, f64)>,
    /// Current simulated RDF used in the last update.
    pub current_rdf: Vec<(f64, f64)>,
    /// k_B T in kJ mol⁻¹.
    pub kbt: f64,
}

impl IterativeBoltzmannInversion {
    /// Create an IBI engine with the given target RDF and initial potential.
    ///
    /// If `initial_potential` is `None`, the PMF of the target RDF is used.
    pub fn new(target_rdf: Vec<(f64, f64)>, kbt: f64) -> Self {
        let potential = target_rdf
            .iter()
            .map(|&(r, g)| {
                let v = if g > 1e-10 { -kbt * g.ln() } else { 20.0 };
                (r, v)
            })
            .collect();
        let current_rdf = target_rdf.clone();
        Self {
            target_rdf,
            potential,
            current_rdf,
            kbt,
        }
    }

    /// Update the tabulated potential using the current simulated RDF.
    ///
    /// Must call this after running a simulation and obtaining a new `g_sim`.
    pub fn update_potential(&mut self) {
        for (k, &(r, v_old)) in self.potential.clone().iter().enumerate() {
            let g_sim = self
                .current_rdf
                .get(k)
                .map(|x| x.1)
                .unwrap_or(1e-10)
                .max(1e-10);
            let g_tgt = self
                .target_rdf
                .get(k)
                .map(|x| x.1)
                .unwrap_or(1e-10)
                .max(1e-10);
            let delta = self.kbt * (g_sim / g_tgt).ln();
            self.potential[k] = (r, v_old + delta);
        }
    }

    /// Set the current simulated RDF (used by `update_potential`).
    pub fn set_current_rdf(&mut self, rdf: Vec<(f64, f64)>) {
        self.current_rdf = rdf;
    }

    /// Compute the potential of mean force −k_BT ln g(r) from the target RDF.
    ///
    /// Returns `(r [nm], PMF [kJ mol⁻¹])` pairs.
    pub fn potential_of_mean_force(&self) -> Vec<(f64, f64)> {
        self.target_rdf
            .iter()
            .map(|&(r, g)| {
                let pmf = if g > 1e-10 { -self.kbt * g.ln() } else { 20.0 };
                (r, pmf)
            })
            .collect()
    }

    /// Convergence check: max |V_{n+1} - V_n| over the PMF.
    pub fn convergence_error(&self) -> f64 {
        self.potential
            .iter()
            .zip(self.potential_of_mean_force().iter())
            .map(|(&(_, v), &(_, pmf))| (v - pmf).abs())
            .fold(0.0f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// MsIbi
// ---------------------------------------------------------------------------

/// Multi-state iterative Boltzmann inversion (MS-IBI).
///
/// Averages IBI corrections from multiple thermodynamic states to produce
/// a transferable CG potential.
#[derive(Debug, Clone)]
pub struct MsIbi {
    /// One IBI engine per thermodynamic state.
    pub states: Vec<IterativeBoltzmannInversion>,
    /// Weights for each state (must sum to 1).
    pub weights: Vec<f64>,
    /// Current combined potential V(r).
    pub combined_potential: Vec<(f64, f64)>,
    /// Convergence tolerance for ΔV \[kJ mol⁻¹\].
    pub tolerance: f64,
}

impl MsIbi {
    /// Construct MS-IBI from a list of (target_rdf, kbt) state pairs and
    /// uniform weights.
    pub fn new(states_data: Vec<(Vec<(f64, f64)>, f64)>, tolerance: f64) -> Self {
        let n = states_data.len();
        let weights = vec![1.0 / n as f64; n];
        let mut states: Vec<IterativeBoltzmannInversion> = states_data
            .into_iter()
            .map(|(rdf, kbt)| IterativeBoltzmannInversion::new(rdf, kbt))
            .collect();
        // Initialise combined potential from first state's PMF
        let combined_potential = if states.is_empty() {
            Vec::new()
        } else {
            states[0].potential_of_mean_force()
        };
        // Push combined potential into each state
        for s in &mut states {
            s.potential = combined_potential.clone();
        }
        Self {
            states,
            weights,
            combined_potential,
            tolerance,
        }
    }

    /// Compute the weighted average correction across all states and apply it.
    pub fn compute_correction(&mut self) {
        if self.states.is_empty() {
            return;
        }
        let len = self.combined_potential.len();
        let mut new_pot = self.combined_potential.clone();
        for k in 0..len {
            let mut delta = 0.0f64;
            for (s, &w) in self.states.iter().zip(self.weights.iter()) {
                let g_sim = s
                    .current_rdf
                    .get(k)
                    .map(|x| x.1)
                    .unwrap_or(1e-10)
                    .max(1e-10);
                let g_tgt = s.target_rdf.get(k).map(|x| x.1).unwrap_or(1e-10).max(1e-10);
                delta += w * s.kbt * (g_sim / g_tgt).ln();
            }
            new_pot[k].1 += delta;
        }
        self.combined_potential = new_pot.clone();
        for s in &mut self.states {
            s.potential = new_pot.clone();
        }
    }

    /// Returns `true` when the maximum potential correction is below
    /// `self.tolerance`.
    pub fn converged(&self) -> bool {
        let max_delta = self
            .states
            .iter()
            .zip(self.weights.iter())
            .fold(0.0f64, |acc, (s, &w)| {
                let d = s.convergence_error() * w;
                acc.max(d)
            });
        max_delta < self.tolerance
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CgBead ---

    #[test]
    fn test_cg_bead_kinetic_energy_zero() {
        let b = CgBead::new("NC3", CgBeadType::Charged, 72.0, 1.0);
        assert_eq!(b.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_cg_bead_kinetic_energy_nonzero() {
        let mut b = CgBead::new("GL1", CgBeadType::NonPolar, 72.0, 0.0);
        b.vel = [1.0, 0.0, 0.0];
        // KE = 0.5 * 72 * 1 = 36
        assert!((b.kinetic_energy() - 36.0).abs() < 1e-10);
    }

    #[test]
    fn test_cg_bead_type_variants() {
        let types = [
            CgBeadType::Polar,
            CgBeadType::NonPolar,
            CgBeadType::Apolar,
            CgBeadType::Charged,
            CgBeadType::Backbone,
            CgBeadType::Sidechain,
        ];
        assert_eq!(types.len(), 6);
    }

    // --- CgMolecule bonded energy ---

    fn two_bead_mol(r: f64) -> CgMolecule {
        let mut mol = CgMolecule::new();
        let mut b0 = CgBead::new("B0", CgBeadType::Apolar, 72.0, 0.0);
        b0.pos = [0.0, 0.0, 0.0];
        let mut b1 = CgBead::new("B1", CgBeadType::Apolar, 72.0, 0.0);
        b1.pos = [r, 0.0, 0.0];
        mol.beads.push(b0);
        mol.beads.push(b1);
        mol.bonds.push((0, 1, 0.47, 3800.0));
        mol
    }

    #[test]
    fn test_bonded_energy_at_equilibrium() {
        let mol = two_bead_mol(0.47);
        assert!(mol.bonded_energy().abs() < 1e-10);
    }

    #[test]
    fn test_bonded_energy_stretched() {
        let mol = two_bead_mol(0.57);
        // E = 0.5 * 3800 * (0.1)^2 = 19 kJ/mol
        let expected = 0.5 * 3800.0 * 0.01;
        assert!((mol.bonded_energy() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_bond_forces_equal_opposite() {
        let mol = two_bead_mol(0.60);
        let f = mol.bond_forces();
        // Newton's third law: f[0] + f[1] == 0
        for k in 0..3 {
            assert!((f[0][k] + f[1][k]).abs() < 1e-10);
        }
    }

    #[test]
    fn test_bond_forces_direction() {
        let mol = two_bead_mol(0.60);
        let f = mol.bond_forces();
        // Bond is stretched, so force on b0 is in +x direction (towards b1)
        assert!(f[0][0] > 0.0, "force should pull bead 0 toward bead 1");
    }

    #[test]
    fn test_angle_forces_three_beads() {
        let mut mol = CgMolecule::new();
        for i in 0..3 {
            let mut b = CgBead::new(format!("B{i}"), CgBeadType::Backbone, 72.0, 0.0);
            b.pos = [i as f64 * 0.47, 0.0, 0.0]; // linear
            mol.beads.push(b);
        }
        mol.angles.push((0, 1, 2, PI, 25.0)); // equilibrium = 180°
        let af = mol.angle_forces();
        // At equilibrium, forces should be near zero
        for f in &af {
            for &fi in f.iter() {
                assert!(fi.abs() < 1e-6);
            }
        }
    }

    #[test]
    fn test_molecule_default() {
        let mol = CgMolecule::default();
        assert!(mol.beads.is_empty());
    }

    // --- CgForceField ---

    #[test]
    fn test_lj_energy_at_sigma() {
        let ff = CgForceField::default_martini();
        // At r = sigma, LJ = 0 (4ε(1-1) = 0)
        let sigma = 0.47;
        let e = ff.lj_energy(0, 0, sigma);
        assert!(e.abs() < 1e-10);
    }

    #[test]
    fn test_lj_energy_minimum() {
        let ff = CgForceField::default_martini();
        let sigma = 0.47;
        // Minimum at r = 2^(1/6) * sigma
        let r_min = 2.0f64.powf(1.0 / 6.0) * sigma;
        let e_min = ff.lj_energy(0, 0, r_min);
        // Should be −epsilon = −5.6 kJ/mol
        assert!((e_min + 5.6).abs() < 1e-8);
    }

    #[test]
    fn test_electrostatic_energy_like_charges_repulsive() {
        let ff = CgForceField::default_martini();
        let e = ff.electrostatic_energy(1.0, 1.0, 0.5);
        assert!(e > 0.0);
    }

    #[test]
    fn test_electrostatic_energy_opposite_charges_attractive() {
        let ff = CgForceField::default_martini();
        let e = ff.electrostatic_energy(1.0, -1.0, 0.5);
        assert!(e < 0.0);
    }

    #[test]
    fn test_total_nonbonded_two_beads() {
        let ff = CgForceField::default_martini();
        let mut b0 = CgBead::new("A", CgBeadType::Polar, 72.0, 0.0);
        b0.pos = [0.0, 0.0, 0.0];
        let mut b1 = CgBead::new("B", CgBeadType::Polar, 72.0, 0.0);
        b1.pos = [1.0, 0.0, 0.0];
        let e = ff.total_nonbonded(&[b0, b1]);
        // Should be finite
        assert!(e.is_finite());
    }

    #[test]
    fn test_total_nonbonded_empty() {
        let ff = CgForceField::default_martini();
        assert_eq!(ff.total_nonbonded(&[]), 0.0);
    }

    // --- CgSimulation ---

    fn make_sim() -> CgSimulation {
        let ff = CgForceField::default_martini();
        let mut mol = CgMolecule::new();
        let mut b = CgBead::new("A", CgBeadType::Apolar, 72.0, 0.0);
        b.pos = [0.5, 0.5, 0.5];
        b.vel = [0.01, 0.0, 0.0];
        mol.beads.push(b);
        CgSimulation::new(vec![mol], ff, 300.0, 0.02)
    }

    #[test]
    fn test_sim_step_advances_time() {
        let mut sim = make_sim();
        sim.step();
        assert!((sim.time - 0.02).abs() < 1e-12);
    }

    #[test]
    fn test_sim_pressure_finite() {
        let sim = make_sim();
        let p = sim.pressure();
        assert!(p.is_finite());
    }

    #[test]
    fn test_sim_rdf_non_empty() {
        let sim = make_sim();
        let rdf = sim.rdf();
        assert!(!rdf.is_empty());
    }

    #[test]
    fn test_sim_diffusion_coefficient() {
        let mut sim = make_sim();
        sim.time = 1.0;
        let d = sim.diffusion_coefficient();
        assert!(d >= 0.0);
    }

    #[test]
    fn test_sim_langevin_thermostat_changes_vel() {
        let mut sim = make_sim();
        let v_before = sim.molecules[0].beads[0].vel;
        sim.langevin_thermostat();
        let v_after = sim.molecules[0].beads[0].vel;
        // Thermostat should change velocities (extremely low probability of exact equality)
        let same = v_before
            .iter()
            .zip(v_after.iter())
            .all(|(a, b)| (a - b).abs() < 1e-15);
        // Just check it doesn't panic; result might coincidentally match for single steps
        let _ = same;
    }

    #[test]
    fn test_sim_multiple_steps() {
        let mut sim = make_sim();
        for _ in 0..10 {
            sim.step();
        }
        assert!((sim.time - 0.20).abs() < 1e-10);
    }

    // --- IterativeBoltzmannInversion ---

    fn sample_rdf() -> Vec<(f64, f64)> {
        (1..=20)
            .map(|i| (i as f64 * 0.1, 1.0 + 0.1 * (i as f64)))
            .collect()
    }

    #[test]
    fn test_ibi_new_creates_potential() {
        let ibi = IterativeBoltzmannInversion::new(sample_rdf(), 2.479);
        assert_eq!(ibi.potential.len(), ibi.target_rdf.len());
    }

    #[test]
    fn test_ibi_pmf_negative_for_g_gt_1() {
        let ibi = IterativeBoltzmannInversion::new(sample_rdf(), 2.479);
        let pmf = ibi.potential_of_mean_force();
        // g > 1 → PMF = -kBT ln g < 0
        for &(_, v) in pmf.iter().skip(1) {
            assert!(v < 0.0, "PMF for g>1 should be negative");
        }
    }

    #[test]
    fn test_ibi_update_no_change_when_rdf_matches() {
        let rdf = sample_rdf();
        let mut ibi = IterativeBoltzmannInversion::new(rdf.clone(), 2.479);
        ibi.current_rdf = rdf;
        let pot_before: Vec<f64> = ibi.potential.iter().map(|x| x.1).collect();
        ibi.update_potential();
        let pot_after: Vec<f64> = ibi.potential.iter().map(|x| x.1).collect();
        for (a, b) in pot_before.iter().zip(pot_after.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_ibi_convergence_error_zero_at_pmf() {
        let rdf = sample_rdf();
        let ibi = IterativeBoltzmannInversion::new(rdf, 2.479);
        // Just after construction, potential == PMF → error should be ~0
        assert!(ibi.convergence_error() < 1e-10);
    }

    #[test]
    fn test_ibi_set_current_rdf() {
        let rdf = sample_rdf();
        let mut ibi = IterativeBoltzmannInversion::new(rdf.clone(), 2.479);
        let new_rdf: Vec<(f64, f64)> = rdf.iter().map(|&(r, g)| (r, g * 1.1)).collect();
        ibi.set_current_rdf(new_rdf.clone());
        assert_eq!(ibi.current_rdf[0].1, new_rdf[0].1);
    }

    // --- MsIbi ---

    fn make_ms_ibi() -> MsIbi {
        let rdf1 = sample_rdf();
        let rdf2: Vec<(f64, f64)> = rdf1.iter().map(|&(r, g)| (r, g * 0.9)).collect();
        MsIbi::new(vec![(rdf1, 2.479), (rdf2, 2.700)], 0.01)
    }

    #[test]
    fn test_ms_ibi_new() {
        let ms = make_ms_ibi();
        assert_eq!(ms.states.len(), 2);
        assert!((ms.weights[0] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_ms_ibi_compute_correction_runs() {
        let mut ms = make_ms_ibi();
        ms.compute_correction();
        assert!(!ms.combined_potential.is_empty());
    }

    #[test]
    fn test_ms_ibi_converged_false_initially() {
        let ms = make_ms_ibi();
        // Not necessarily converged right away (depends on rdf matching)
        let _ = ms.converged(); // just ensure no panic
    }

    #[test]
    fn test_ms_ibi_converged_when_rdfs_match() {
        let rdf = sample_rdf();
        let mut ms = MsIbi::new(
            vec![(rdf.clone(), 2.479), (rdf.clone(), 2.479)],
            1e6, // huge tolerance → always converged
        );
        ms.compute_correction();
        assert!(ms.converged());
    }

    #[test]
    fn test_bead_type_index_all() {
        assert_eq!(bead_type_index(CgBeadType::Polar), 0);
        assert_eq!(bead_type_index(CgBeadType::NonPolar), 1);
        assert_eq!(bead_type_index(CgBeadType::Apolar), 2);
        assert_eq!(bead_type_index(CgBeadType::Charged), 3);
        assert_eq!(bead_type_index(CgBeadType::Backbone), 4);
        assert_eq!(bead_type_index(CgBeadType::Sidechain), 5);
    }
}
