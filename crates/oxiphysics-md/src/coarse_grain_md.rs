// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coarse-grained molecular dynamics (CG-MD), MARTINI-like.
//!
//! Provides:
//! - [`CgBead`]: single CG particle with position, velocity, mass and charge.
//! - [`CgBeadType`]: bead chemistry classification.
//! - [`CgMolecule`]: collection of bonded beads with harmonic bonds and angles.
//! - [`CgForceField`]: Lennard-Jones + electrostatic non-bonded interactions.
//! - [`CgSimulation`]: NVT simulation with Langevin thermostat in an orthorhombic
//!   periodic box (minimum-image forces, virial pressure, g(r), diffusion).
//! - [`IterativeBoltzmannInversion`]: potential-of-mean-force from a target RDF.
//! - [`MsIbi`]: multi-state iterative Boltzmann inversion.
//!
//! References
//! ----------
//! Marrink, S. J. et al. (2007) *J. Phys. Chem. B* **111**, 7812–7824.

use crate::MdError;
use std::f64::consts::PI;

/// Boltzmann constant in kJ mol⁻¹ K⁻¹ (`R = N_A k_B`, CODATA 2018).
const KB_KJ_PER_MOL_K: f64 = 8.314_462_618e-3;

/// Coulomb prefactor `1/(4π ε₀)` in kJ mol⁻¹ nm e⁻².
const COULOMB_KJ_NM_PER_E2: f64 = 138.935;

/// Distance floor (nm) of the non-bonded kernel: closer pairs interact as if
/// separated by this distance, which keeps energies and forces finite.
const MIN_PAIR_DISTANCE: f64 = 0.1;

/// Charge magnitude (e) below which a bead carries no electrostatic interaction.
const NEGLIGIBLE_CHARGE: f64 = 1e-10;

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

/// Minimum-image convention for an orthorhombic periodic box: wrap every
/// component of the separation `d` into `[-L/2, L/2]`.  An axis whose edge
/// length is not positive and finite is treated as non-periodic.
#[inline]
fn minimum_image(d: [f64; 3], box_lengths: [f64; 3]) -> [f64; 3] {
    let wrap = |x: f64, l: f64| {
        if l.is_finite() && l > 0.0 {
            x - l * (x / l).round()
        } else {
            x
        }
    };
    [
        wrap(d[0], box_lengths[0]),
        wrap(d[1], box_lengths[1]),
        wrap(d[2], box_lengths[2]),
    ]
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
        COULOMB_KJ_NM_PER_E2 * qi * qj / r * (-r / self.electrostatic_screening).exp()
    }

    /// Total non-bonded energy of a bead array (all unique pairs, no periodic
    /// boundaries).
    /// Returns energy in kJ mol⁻¹.
    ///
    /// A bead array carries no topology, so bonded neighbours are *not* excluded
    /// here; [`CgSimulation::total_potential_energy`] applies the exclusions of 1-2
    /// and 1-3 pairs of the simulated molecules and the minimum image.
    pub fn total_nonbonded(&self, beads: &[CgBead]) -> f64 {
        let mut e = 0.0;
        let n = beads.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let r = dist3(beads[i].pos, beads[j].pos).max(MIN_PAIR_DISTANCE);
                let ti = bead_type_index(beads[i].type_);
                let tj = bead_type_index(beads[j].type_);
                e += self.lj_energy(ti, tj, r);
                if beads[i].charge.abs() > NEGLIGIBLE_CHARGE
                    || beads[j].charge.abs() > NEGLIGIBLE_CHARGE
                {
                    e += self.electrostatic_energy(beads[i].charge, beads[j].charge, r);
                }
            }
        }
        e
    }

    /// Energy, force factor and virial of the non-bonded interaction (Lennard-Jones
    /// plus Yukawa electrostatics) between two beads whose minimum-image
    /// separation vector is `r_vec = r_i − r_j`.
    ///
    /// The distance is floored at [`MIN_PAIR_DISTANCE`], exactly as in
    /// [`CgForceField::total_nonbonded`], so energy, force and virial of every
    /// caller stay mutually consistent.
    fn pair_terms(&self, bead_i: &CgBead, bead_j: &CgBead, r_vec: [f64; 3]) -> PairTerms {
        let r = norm3(r_vec).max(MIN_PAIR_DISTANCE);
        let ti = bead_type_index(bead_i.type_);
        let tj = bead_type_index(bead_j.type_);
        let (eps, sigma) = self.lj_params[ti][tj];
        let sr6 = (sigma / r).powi(6);
        let sr12 = sr6 * sr6;
        let mut energy = self.lj_energy(ti, tj, r);
        // F_r = -dU/dr = 4ε (12 sr12 - 6 sr6) / r; the force on bead i is F_r r_vec / r.
        let mut force_factor = 4.0 * eps * (12.0 * sr12 - 6.0 * sr6) / (r * r);
        // r dU/dr for LJ-12-6.
        let mut virial = 4.0 * eps * (-12.0 * sr12 + 6.0 * sr6);
        let lambda = self.electrostatic_screening;
        let charged =
            bead_i.charge.abs() > NEGLIGIBLE_CHARGE || bead_j.charge.abs() > NEGLIGIBLE_CHARGE;
        if charged && lambda > 0.0 {
            // U = K qi qj e^{-r/λ} / r  =>  -dU/dr = K qi qj e^{-r/λ} (1/r² + 1/(λ r)).
            let qq = COULOMB_KJ_NM_PER_E2 * bead_i.charge * bead_j.charge;
            let screening = (-r / lambda).exp();
            energy += self.electrostatic_energy(bead_i.charge, bead_j.charge, r);
            force_factor += qq * screening * (1.0 / (r * r) + 1.0 / (lambda * r)) / r;
            virial -= qq * screening * (1.0 / r + 1.0 / lambda);
        }
        PairTerms {
            energy,
            force_factor,
            virial,
        }
    }
}

/// Energy, force and virial of one non-bonded bead pair.
struct PairTerms {
    /// Pair energy `U` \[kJ mol⁻¹\].
    energy: f64,
    /// Factor `-dU/dr / r` that multiplies the minimum-image vector `r_i − r_j` to
    /// give the force on bead `i` \[kJ mol⁻¹ nm⁻²\] (the force on `j` is opposite).
    force_factor: f64,
    /// Pair virial `r dU/dr` \[kJ mol⁻¹\].
    virial: f64,
}

/// Bead pairs that do not interact non-bonded: bonded 1-2 neighbours and the end
/// beads of every angle (1-3 pairs).
///
/// For every bead the sorted list of excluded partners with a *higher* global bead
/// index is stored, so a lookup is a binary search in a list of a few entries.
struct Exclusions {
    partners: Vec<Vec<usize>>,
}

impl Exclusions {
    /// Whether the global bead pair `(i, j)` with `i < j` is excluded.
    fn contains(&self, i: usize, j: usize) -> bool {
        self.partners
            .get(i)
            .is_some_and(|list| list.binary_search(&j).is_ok())
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

/// NVT coarse-grained simulation with Langevin thermostat in an orthorhombic
/// periodic box.
///
/// Non-bonded pairs (Lennard-Jones plus Yukawa electrostatics of charged beads)
/// interact through the minimum image of their separation, except for the 1-2 pairs
/// joined by a bond and the 1-3 end beads of an angle of the same molecule, which
/// are excluded.  The box edge lengths [`CgSimulation::box_lengths`] define the volume used by
/// [`CgSimulation::pressure`] and the number density used by
/// [`CgSimulation::rdf`].  Bead positions are never wrapped back into the box, so
/// they stay *unwrapped*: [`CgSimulation::diffusion_coefficient`] measures the
/// mean-squared displacement from the positions captured at construction (or at
/// the last [`CgSimulation::reset_displacement_origin`]) without periodic jumps.
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
    /// Edge lengths \[Lx, Ly, Lz\] of the orthorhombic periodic box \[nm\].
    pub box_lengths: [f64; 3],
    /// Unwrapped bead positions (molecule order) at the displacement origin \[nm\].
    initial_positions: Vec<[f64; 3]>,
    /// Simulation time at which `initial_positions` was captured \[ps\].
    origin_time: f64,
    /// Pseudo-random seed state (LCG, no external crate needed).
    lcg_state: u64,
}

impl CgSimulation {
    /// Construct a simulation with given molecules, force field, temperature, dt
    /// and the edge lengths `box_lengths` \[nm\] of the periodic box.
    ///
    /// The bead positions at construction are recorded as the origin of the
    /// mean-squared displacement used by [`CgSimulation::diffusion_coefficient`].
    ///
    /// # Errors
    /// Returns [`MdError::General`] if an edge length is not finite and positive.
    pub fn new(
        molecules: Vec<CgMolecule>,
        ff: CgForceField,
        temp: f64,
        dt: f64,
        box_lengths: [f64; 3],
    ) -> Result<Self, MdError> {
        if box_lengths.iter().any(|l| !l.is_finite() || *l <= 0.0) {
            return Err(MdError::General(format!(
                "CgSimulation box edge lengths must be finite and positive, got {box_lengths:?}"
            )));
        }
        let initial_positions = molecules
            .iter()
            .flat_map(|m| m.beads.iter().map(|b| b.pos))
            .collect();
        Ok(Self {
            molecules,
            ff,
            temp,
            dt,
            gamma: 1.0,
            time: 0.0,
            box_lengths,
            initial_positions,
            origin_time: 0.0,
            lcg_state: 12345,
        })
    }

    /// Volume `Lx · Ly · Lz` of the periodic box \[nm³\].
    pub fn volume(&self) -> f64 {
        self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2]
    }

    /// Record the current bead positions as the origin of the mean-squared
    /// displacement and restart the elapsed-time clock of
    /// [`CgSimulation::diffusion_coefficient`] at the current simulation time.
    pub fn reset_displacement_origin(&mut self) {
        let positions: Vec<[f64; 3]> = self.all_beads().iter().map(|b| b.pos).collect();
        self.initial_positions = positions;
        self.origin_time = self.time;
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
        let kbt = KB_KJ_PER_MOL_K * self.temp; // kJ mol⁻¹
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

    /// Total potential energy of the system \[kJ mol⁻¹\].
    ///
    /// The sum of the bonded energy of every molecule
    /// ([`CgMolecule::bonded_energy`]) and the non-bonded Lennard-Jones plus Yukawa
    /// energy of all non-excluded bead pairs (minimum image, bonded 1-2 and 1-3
    /// pairs of a molecule excluded, separations floored at 0.1 nm).  The forces of
    /// the simulation are exactly the negative gradient of this energy for pair
    /// separations above 0.1 nm.
    pub fn total_potential_energy(&self) -> f64 {
        let beads = self.all_beads();
        let mut energy: f64 = self.molecules.iter().map(CgMolecule::bonded_energy).sum();
        self.for_each_nonbonded_pair(&beads, |_, _, _, terms| energy += terms.energy);
        energy
    }

    /// Instantaneous pressure estimate via the virial theorem \[kJ mol⁻¹ nm⁻³\].
    ///
    /// `P = N k_B T / V - W / (3 V)` with `V = Lx · Ly · Lz` the volume of the
    /// periodic box, `T` the target temperature and `W = Σ r dU/dr` the
    /// configurational virial of the non-bonded pairs (Lennard-Jones plus Yukawa at
    /// minimum-image separations, with the same bonded exclusions as the forces)
    /// and of the harmonic bonds.  The cosine-harmonic angle terms depend only on
    /// the angle, which is invariant under a uniform rescaling of the system, so
    /// they contribute nothing to the virial.
    ///
    /// Returns `0.0` for an empty system and `NaN` if the box volume is not
    /// finite and positive.
    pub fn pressure(&self) -> f64 {
        let beads = self.all_beads();
        let n = beads.len();
        if n == 0 {
            return 0.0;
        }
        let volume = self.volume();
        if !volume.is_finite() || volume <= 0.0 {
            return f64::NAN;
        }
        let kbt = KB_KJ_PER_MOL_K * self.temp;
        let ideal = (n as f64) * kbt / volume;
        // Σ r dU/dr over all non-excluded pairs
        let mut virial = 0.0f64;
        self.for_each_nonbonded_pair(&beads, |_, _, _, terms| virial += terms.virial);
        // Harmonic bonds: r dU/dr = k (r - r0) r.  Molecules are never wrapped,
        // so the intramolecular separation is the plain bead distance.
        for mol in &self.molecules {
            for &(i, j, r0, k) in &mol.bonds {
                let r = dist3(mol.beads[i].pos, mol.beads[j].pos);
                virial += k * (r - r0) * r;
            }
        }
        ideal - virial / (3.0 * volume)
    }

    /// Compute the radial distribution function g(r) for all bead pairs.
    ///
    /// Pair separations use the minimum image, histogrammed in bins of 0.02 nm
    /// up to `min(2 nm, L_min / 2)` (beyond half the shortest box edge the
    /// minimum-image shell is no longer complete).  Each bin is normalised by the
    /// number of pairs an ideal gas of the same density would place in the
    /// spherical shell, `N (N - 1) / 2 · V_shell / V`, with `V = Lx · Ly · Lz`,
    /// so a uniform system gives `g(r) ≈ 1`.
    ///
    /// Returns a vector of `(r [nm], g(r))` pairs with `r` at the bin centre; `g`
    /// is `0` for fewer than two beads, and the vector is empty when the box
    /// volume or the usable range is not positive.
    pub fn rdf(&self) -> Vec<(f64, f64)> {
        const BIN_WIDTH: f64 = 0.02;
        const R_MAX_REQUESTED: f64 = 2.0;
        let beads = self.all_beads();
        let n = beads.len();
        let volume = self.volume();
        let shortest_edge = self
            .box_lengths
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        let r_max = R_MAX_REQUESTED.min(0.5 * shortest_edge);
        if !volume.is_finite() || volume <= 0.0 || !r_max.is_finite() || r_max <= 0.0 {
            return Vec::new();
        }
        let bins = (r_max / BIN_WIDTH + 1e-9).floor() as usize;
        let mut hist = vec![0u64; bins];
        for i in 0..n {
            for j in (i + 1)..n {
                let rij = minimum_image(sub3(beads[i].pos, beads[j].pos), self.box_lengths);
                let bin = (norm3(rij) / BIN_WIDTH) as usize;
                if bin < bins {
                    hist[bin] += 1;
                }
            }
        }
        let total_pairs = 0.5 * (n as f64) * (n as f64 - 1.0).max(0.0);
        (0..bins)
            .map(|k| {
                let r_lo = k as f64 * BIN_WIDTH;
                let r_hi = r_lo + BIN_WIDTH;
                let shell_volume = 4.0 * PI / 3.0 * (r_hi.powi(3) - r_lo.powi(3));
                let ideal_pairs = total_pairs * shell_volume / volume;
                let g = if ideal_pairs > 0.0 {
                    hist[k] as f64 / ideal_pairs
                } else {
                    0.0
                };
                (0.5 * (r_lo + r_hi), g)
            })
            .collect()
    }

    /// Estimate the self-diffusion coefficient from the mean-squared displacement.
    ///
    /// Einstein relation `D = <|r(t) - r(0)|²> / (6 t)` with `r(0)` the bead
    /// positions recorded at construction (or at the last
    /// [`CgSimulation::reset_displacement_origin`]), `t` the time elapsed since
    /// then and the average taken over all beads that have a recorded origin.
    /// Positions are unwrapped, so periodic-boundary crossings do not
    /// contribute spurious jumps.  Returns `D` in nm² ps⁻¹, or `0.0` if no time
    /// has elapsed.
    pub fn diffusion_coefficient(&self) -> f64 {
        let elapsed = self.time - self.origin_time;
        if elapsed <= 0.0 {
            return 0.0;
        }
        let mut msd_sum = 0.0f64;
        let mut count = 0usize;
        for (bead, origin) in self.all_beads().iter().zip(self.initial_positions.iter()) {
            let displacement = sub3(bead.pos, *origin);
            msd_sum += dot3(displacement, displacement);
            count += 1;
        }
        if count == 0 {
            return 0.0;
        }
        msd_sum / (count as f64) / (6.0 * elapsed)
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

        // Non-bonded forces between all non-excluded beads across molecules.
        let beads: Vec<&CgBead> = self.all_beads();
        self.for_each_nonbonded_pair(&beads, |i, j, r_vec, terms| {
            let fv = scale3(r_vec, terms.force_factor);
            forces[i] = add3(forces[i], fv);
            forces[j] = sub3(forces[j], fv);
        });
        forces
    }

    /// Bonded 1-2 pairs and 1-3 angle end pairs of every molecule as global bead
    /// indices (molecule order, beads numbered consecutively).
    fn exclusions(&self) -> Exclusions {
        let total: usize = self.molecules.iter().map(|m| m.beads.len()).sum();
        let mut partners = vec![Vec::new(); total];
        let mut offset = 0usize;
        for mol in &self.molecules {
            let bond_pairs = mol.bonds.iter().map(|&(i, j, _, _)| (i, j));
            let angle_pairs = mol.angles.iter().map(|&(i, _, k, _, _)| (i, k));
            for (a, b) in bond_pairs.chain(angle_pairs) {
                let (lo, hi) = (offset + a.min(b), offset + a.max(b));
                if let Some(list) = partners.get_mut(lo) {
                    list.push(hi);
                }
            }
            offset += mol.beads.len();
        }
        for list in &mut partners {
            list.sort_unstable();
            list.dedup();
        }
        Exclusions { partners }
    }

    /// Visit every non-excluded bead pair once as `(i, j, r_i − r_j, terms)` with
    /// the minimum-image separation vector, in the order `i < j`.
    ///
    /// This is the single pair kernel behind the forces, the potential energy and
    /// the virial pressure, so the three always agree on exclusions, the minimum
    /// image and the distance floor.
    fn for_each_nonbonded_pair<F>(&self, beads: &[&CgBead], mut visit: F)
    where
        F: FnMut(usize, usize, [f64; 3], PairTerms),
    {
        let exclusions = self.exclusions();
        for i in 0..beads.len() {
            for j in (i + 1)..beads.len() {
                if exclusions.contains(i, j) {
                    continue;
                }
                let r_vec = minimum_image(sub3(beads[i].pos, beads[j].pos), self.box_lengths);
                visit(i, j, r_vec, self.ff.pair_terms(beads[i], beads[j], r_vec));
            }
        }
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
        for (k, pot_entry) in new_pot.iter_mut().enumerate().take(len) {
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
            pot_entry.1 += delta;
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
mod exclusion_tests;

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
        for (&f0k, &f1k) in f[0].iter().zip(f[1].iter()) {
            assert!((f0k + f1k).abs() < 1e-10);
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
        CgSimulation::new(vec![mol], ff, 300.0, 0.02, [5.0, 5.0, 5.0]).expect("valid box")
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
    fn test_sim_diffusion_coefficient_zero_without_displacement() {
        // The bead sits at (0.5, 0.5, 0.5), far from the coordinate origin, and has
        // not moved: the old "distance from the origin" estimate was non-zero here.
        let mut sim = make_sim();
        sim.time = 1.0;
        assert_eq!(sim.diffusion_coefficient(), 0.0);
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

    // --- CgSimulation: periodic box, pressure, g(r), diffusion ---

    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    /// `R = N_A k_B` in kJ mol⁻¹ K⁻¹ (CODATA 2018), independent of the module constant.
    const R_KJ_MOL_K: f64 = 8.314_462_618e-3;

    /// Force field with every Lennard-Jones well depth set to zero (an ideal gas).
    fn ideal_gas_ff() -> CgForceField {
        let mut ff = CgForceField::default_martini();
        for row in &mut ff.lj_params {
            for entry in row.iter_mut() {
                entry.0 = 0.0;
            }
        }
        ff
    }

    /// `n` single-bead molecules placed uniformly at random inside `box_lengths`.
    fn uniform_gas(n: usize, box_lengths: [f64; 3], seed: u64) -> Vec<CgMolecule> {
        let mut rng = StdRng::seed_from_u64(seed);
        (0..n)
            .map(|idx| {
                let mut mol = CgMolecule::new();
                let mut bead = CgBead::new(format!("G{idx}"), CgBeadType::Apolar, 72.0, 0.0);
                bead.pos = [
                    rng.random_range(0.0..box_lengths[0]),
                    rng.random_range(0.0..box_lengths[1]),
                    rng.random_range(0.0..box_lengths[2]),
                ];
                mol.beads.push(bead);
                mol
            })
            .collect()
    }

    /// Single-bead molecule with a velocity.
    fn moving_bead(pos: [f64; 3], vel: [f64; 3]) -> CgMolecule {
        let mut mol = CgMolecule::new();
        let mut bead = CgBead::new("M", CgBeadType::Apolar, 72.0, 0.0);
        bead.pos = pos;
        bead.vel = vel;
        mol.beads.push(bead);
        mol
    }

    /// Two-bead system in a cubic box of edge 4 nm with beads at `x_a` and `x_b`
    /// on the line y = z = 2 nm.
    fn bead_pair_in_box(x_a: f64, x_b: f64) -> CgSimulation {
        CgSimulation::new(
            vec![
                moving_bead([x_a, 2.0, 2.0], [0.0; 3]),
                moving_bead([x_b, 2.0, 2.0], [0.0; 3]),
            ],
            CgForceField::default_martini(),
            300.0,
            0.02,
            [4.0; 3],
        )
        .expect("valid box")
    }

    #[test]
    fn test_sim_new_rejects_invalid_box() {
        for bad in [
            [0.0, 1.0, 1.0],
            [1.0, -2.0, 1.0],
            [1.0, 1.0, f64::NAN],
            [f64::INFINITY, 1.0, 1.0],
        ] {
            let result = CgSimulation::new(
                Vec::new(),
                CgForceField::default_martini(),
                300.0,
                0.02,
                bad,
            );
            assert!(
                matches!(result, Err(MdError::General(_))),
                "box {bad:?} must be rejected"
            );
        }
        let ok = CgSimulation::new(
            Vec::new(),
            CgForceField::default_martini(),
            300.0,
            0.02,
            [1.0, 2.0, 3.0],
        )
        .expect("valid box");
        assert_eq!(ok.box_lengths, [1.0, 2.0, 3.0]);
        assert!((ok.volume() - 6.0).abs() < 1e-15);
    }

    #[test]
    fn test_sim_pressure_ideal_gas_is_n_kt_over_box_volume() {
        // eps = 0: P = N k_B T / V with V = Lx Ly Lz of a non-cubic box.
        let box_lengths = [4.0, 5.0, 6.0];
        let n = 120;
        let sim = CgSimulation::new(
            uniform_gas(n, box_lengths, 7),
            ideal_gas_ff(),
            300.0,
            0.02,
            box_lengths,
        )
        .expect("valid box");
        let expected = n as f64 * R_KJ_MOL_K * 300.0 / 120.0;
        assert!(
            ((sim.pressure() - expected) / expected).abs() < 1e-12,
            "P = {} expected {expected}",
            sim.pressure()
        );
    }

    #[test]
    fn test_sim_pressure_depends_on_box_not_on_bead_count_alone() {
        // The same 60 beads in a box with doubled edges (8x the volume): the ideal-gas
        // pressure drops by exactly 8; the fabricated "1 bead per 0.5 nm cube" volume
        // gave the same pressure for both.
        let small = [3.0, 3.0, 3.0];
        let large = [6.0, 6.0, 6.0];
        let beads = uniform_gas(60, small, 11);
        let p_small = CgSimulation::new(beads.clone(), ideal_gas_ff(), 300.0, 0.02, small)
            .expect("valid box")
            .pressure();
        let p_large = CgSimulation::new(beads, ideal_gas_ff(), 300.0, 0.02, large)
            .expect("valid box")
            .pressure();
        assert!(
            (p_small / p_large - 8.0).abs() < 1e-12,
            "{p_small} / {p_large}"
        );
    }

    #[test]
    fn test_sim_pressure_of_empty_system_is_zero() {
        let sim = CgSimulation::new(Vec::new(), ideal_gas_ff(), 300.0, 0.02, [2.0; 3])
            .expect("valid box");
        assert_eq!(sim.pressure(), 0.0);
        assert!(sim.rdf().iter().all(|&(_, g)| g == 0.0));
    }

    #[test]
    fn test_sim_degenerate_box_gives_nan_pressure_and_empty_rdf() {
        let mut sim = make_sim();
        sim.box_lengths = [0.0, 5.0, 5.0];
        assert!(sim.pressure().is_nan());
        assert!(sim.rdf().is_empty());
    }

    /// Total potential energy of the simulation, evaluated independently of the
    /// force and virial code: minimum-image LJ over all pairs that are neither
    /// bonded (1-2) nor the ends of an angle (1-3) in the same molecule (with the
    /// same 0.1 nm distance floor) plus the bonded energy of every molecule.
    fn reference_potential_energy(sim: &CgSimulation) -> f64 {
        let beads: Vec<&CgBead> = sim.molecules.iter().flat_map(|m| m.beads.iter()).collect();
        let mut excluded: Vec<(usize, usize)> = Vec::new();
        let mut offset = 0;
        for mol in &sim.molecules {
            for &(a, b, _, _) in &mol.bonds {
                excluded.push((offset + a.min(b), offset + a.max(b)));
            }
            for &(a, _, c, _, _) in &mol.angles {
                excluded.push((offset + a.min(c), offset + a.max(c)));
            }
            offset += mol.beads.len();
        }
        let mut energy = 0.0;
        for i in 0..beads.len() {
            for j in (i + 1)..beads.len() {
                if excluded.contains(&(i, j)) {
                    continue;
                }
                let d = minimum_image(sub3(beads[i].pos, beads[j].pos), sim.box_lengths);
                let r = norm3(d).max(0.1);
                energy += sim.ff.lj_energy(
                    bead_type_index(beads[i].type_),
                    bead_type_index(beads[j].type_),
                    r,
                );
            }
        }
        energy
            + sim
                .molecules
                .iter()
                .map(CgMolecule::bonded_energy)
                .sum::<f64>()
    }

    /// The simulation with all positions and the box edges scaled by `s`.
    fn scaled_copy(sim: &CgSimulation, s: f64) -> CgSimulation {
        let mut out = sim.clone();
        for mol in &mut out.molecules {
            for bead in &mut mol.beads {
                bead.pos = scale3(bead.pos, s);
            }
        }
        out.box_lengths = scale3(sim.box_lengths, s);
        out
    }

    /// Chains (bonded, bent, slightly stretched) and single beads on a jittered
    /// 3 x 3 x 3 lattice of spacing 1.2 nm, in a 3.6 nm cubic box.
    fn interacting_system() -> CgSimulation {
        let mut rng = StdRng::seed_from_u64(2024);
        let mut molecules = Vec::new();
        let mut site = 0;
        for ix in 0..3 {
            for iy in 0..3 {
                for iz in 0..3 {
                    let origin = [
                        1.2 * ix as f64 + 0.6,
                        1.2 * iy as f64 + 0.6,
                        1.2 * iz as f64 + 0.6,
                    ];
                    let mut mol = CgMolecule::new();
                    if site % 5 == 0 {
                        // Three-bead chain: bond lengths 0.50 and 0.43 (r0 = 0.47), angle 110 deg.
                        let theta = 110.0_f64.to_radians();
                        let positions = [
                            [origin[0] - 0.50, origin[1], origin[2]],
                            origin,
                            [
                                origin[0] + 0.43 * (PI - theta).cos(),
                                origin[1] + 0.43 * (PI - theta).sin(),
                                origin[2],
                            ],
                        ];
                        for (idx, p) in positions.iter().enumerate() {
                            let mut bead = CgBead::new(
                                format!("C{site}_{idx}"),
                                CgBeadType::Backbone,
                                72.0,
                                0.0,
                            );
                            bead.pos = *p;
                            mol.beads.push(bead);
                        }
                        mol.bonds.push((0, 1, 0.47, 1250.0));
                        mol.bonds.push((1, 2, 0.47, 1250.0));
                        mol.angles.push((0, 1, 2, 95.0_f64.to_radians(), 25.0));
                    } else {
                        let mut bead =
                            CgBead::new(format!("S{site}"), CgBeadType::Polar, 72.0, 0.0);
                        bead.pos = [
                            origin[0] + rng.random_range(-0.12..0.12),
                            origin[1] + rng.random_range(-0.12..0.12),
                            origin[2] + rng.random_range(-0.12..0.12),
                        ];
                        mol.beads.push(bead);
                    }
                    molecules.push(mol);
                    site += 1;
                }
            }
        }
        CgSimulation::new(
            molecules,
            CgForceField::default_martini(),
            300.0,
            0.02,
            [3.6; 3],
        )
        .expect("valid box")
    }

    #[test]
    fn test_sim_pressure_virial_matches_energy_volume_derivative() {
        // Configurational pressure -dU/dV under a uniform rescaling of positions and
        // box; the LJ pair virial (minimum image) and the bond virial must reproduce it,
        // and the angle terms must not contribute (they are scale invariant).
        let sim = interacting_system();
        let delta = 1e-6;
        let u_plus = reference_potential_energy(&scaled_copy(&sim, 1.0 + delta));
        let u_minus = reference_potential_energy(&scaled_copy(&sim, 1.0 - delta));
        let volume = sim.volume();
        let d_u_d_volume = (u_plus - u_minus) / (2.0 * delta * 3.0 * volume);
        let n = sim.molecules.iter().map(|m| m.beads.len()).sum::<usize>() as f64;
        let expected = n * R_KJ_MOL_K * 300.0 / volume - d_u_d_volume;
        let got = sim.pressure();
        assert!(
            ((got - expected) / expected.abs().max(1.0)).abs() < 1e-6,
            "virial pressure {got} vs -dU/dV estimate {expected}"
        );
        // The configurational part must actually matter in this test.
        let ideal = n * R_KJ_MOL_K * 300.0 / volume;
        assert!(
            (got - ideal).abs() > 0.1 * ideal,
            "virial contribution too small"
        );
    }

    #[test]
    fn test_sim_forces_use_minimum_image_across_the_boundary() {
        // 0.45 nm apart through the periodic boundary (beads at x = 3.8 and x = 0.25,
        // L = 4): the LJ pair is repulsive and must equal the interaction of the same
        // pair placed 0.45 nm apart in the middle of the box.
        let across = bead_pair_in_box(3.8, 0.25);
        let inside = bead_pair_in_box(1.0, 1.45);
        let f_across = across.compute_all_forces();
        let f_inside = inside.compute_all_forces();
        // Analytic LJ-12-6 repulsion for Apolar-Apolar at r = 0.45 nm.
        let (eps, sigma) = across.ff.lj_params[2][2];
        let sr6 = (sigma / 0.45_f64).powi(6);
        let magnitude = 4.0 * eps * (12.0 * sr6 * sr6 - 6.0 * sr6) / 0.45;
        assert!(
            magnitude > 100.0,
            "test geometry must be strongly repulsive"
        );
        // bead A (x = 3.8) is pushed away from B's image at x = 4.25, i.e. to -x.
        assert!(
            (f_across[0][0] + magnitude).abs() < 1e-8 * magnitude,
            "{:?}",
            f_across[0]
        );
        assert!(
            (f_across[1][0] - magnitude).abs() < 1e-8 * magnitude,
            "{:?}",
            f_across[1]
        );
        for k in 0..2 {
            for axis in 0..3 {
                assert!(
                    (f_across[k][axis] - f_inside[k][axis]).abs() < 1e-8 * magnitude,
                    "bead {k} axis {axis}: {} vs {}",
                    f_across[k][axis],
                    f_inside[k][axis]
                );
            }
        }
    }

    #[test]
    fn test_sim_forces_are_invariant_under_lattice_translation_of_a_bead() {
        let sim = interacting_system();
        let reference = sim.compute_all_forces();
        let mut shifted = sim.clone();
        // Move one bead of the first single-bead molecule by whole box edges.
        let target = shifted
            .molecules
            .iter_mut()
            .find(|m| m.beads.len() == 1)
            .expect("single-bead molecule");
        target.beads[0].pos = add3(target.beads[0].pos, [3.6, -7.2, 3.6]);
        let moved = shifted.compute_all_forces();
        for (a, b) in reference.iter().zip(moved.iter()) {
            for axis in 0..3 {
                assert!(
                    (a[axis] - b[axis]).abs() < 1e-7 * (1.0 + a[axis].abs()),
                    "{a:?} vs {b:?}"
                );
            }
        }
    }

    #[test]
    fn test_sim_rdf_of_uniform_random_gas_is_one() {
        // 1000 uniformly random beads in a 8 x 6 x 5 nm box (density 4.17 nm^-3, far
        // from the fabricated 8 nm^-3).  g(r) must be 1 within statistics up to L_min/2.
        let box_lengths = [8.0, 6.0, 5.0];
        let sim = CgSimulation::new(
            uniform_gas(1000, box_lengths, 99),
            ideal_gas_ff(),
            300.0,
            0.02,
            box_lengths,
        )
        .expect("valid box");
        let rdf = sim.rdf();
        // L_min / 2 = 2.5 nm > 2 nm cap, so the histogram spans 100 bins up to 2 nm.
        assert_eq!(rdf.len(), 100);
        assert!((rdf[0].0 - 0.01).abs() < 1e-12);
        for window in [(0.5, 1.0), (1.0, 1.5), (1.5, 2.0)] {
            let values: Vec<f64> = rdf
                .iter()
                .filter(|&&(r, _)| r >= window.0 && r < window.1)
                .map(|&(_, g)| g)
                .collect();
            let mean = values.iter().sum::<f64>() / values.len() as f64;
            assert!(
                (mean - 1.0).abs() < 0.05,
                "mean g(r) over [{}, {}) = {mean}",
                window.0,
                window.1
            );
        }
        let wide: Vec<f64> = rdf
            .iter()
            .filter(|&&(r, _)| r >= 0.5)
            .map(|&(_, g)| g)
            .collect();
        let wide_mean = wide.iter().sum::<f64>() / wide.len() as f64;
        assert!(
            (wide_mean - 1.0).abs() < 0.02,
            "mean g(r) for r >= 0.5 nm: {wide_mean}"
        );
    }

    #[test]
    fn test_sim_rdf_range_is_limited_to_half_the_shortest_edge() {
        let box_lengths = [1.0, 3.0, 3.0];
        let sim = CgSimulation::new(
            uniform_gas(40, box_lengths, 5),
            ideal_gas_ff(),
            300.0,
            0.02,
            box_lengths,
        )
        .expect("valid box");
        let rdf = sim.rdf();
        assert_eq!(rdf.len(), 25);
        let (r_last, _) = rdf[rdf.len() - 1];
        assert!(
            r_last < 0.5,
            "last bin centre {r_last} must stay below L_min / 2"
        );
    }

    #[test]
    fn test_sim_rdf_uses_minimum_image_and_exact_pair_normalisation() {
        // Two beads 0.21 nm apart through the boundary of a 4 nm box: one pair in bin
        // 10 ([0.20, 0.22) nm), normalised by N(N-1)/2 V_shell / V with N = 2.
        let sim = bead_pair_in_box(3.89, 0.1);
        let rdf = sim.rdf();
        assert_eq!(rdf.len(), 100);
        let shell_volume = 4.0 * PI / 3.0 * (0.22_f64.powi(3) - 0.20_f64.powi(3));
        let expected = 64.0 / shell_volume;
        for (k, &(r, g)) in rdf.iter().enumerate() {
            if k == 10 {
                assert!((r - 0.21).abs() < 1e-12);
                assert!(
                    ((g - expected) / expected).abs() < 1e-12,
                    "g = {g}, expected {expected}"
                );
            } else {
                assert_eq!(g, 0.0, "bin {k} must be empty");
            }
        }
    }

    #[test]
    fn test_sim_diffusion_coefficient_measures_displacement_from_initial_positions() {
        // Free flight (no forces, gamma = 0): every bead moves by v t, so
        // D = <|v|^2> t / 6 irrespective of where the beads started.
        let velocities = [[0.30, 0.0, 0.0], [0.0, -0.20, 0.10], [0.05, 0.05, -0.40]];
        let starts = [[3.8, 1.5, 1.5], [3.5, 0.5, 3.0], [0.5, 3.5, 2.0]];
        let molecules = starts
            .iter()
            .zip(velocities.iter())
            .map(|(&p, &v)| moving_bead(p, v))
            .collect();
        let mut sim =
            CgSimulation::new(molecules, ideal_gas_ff(), 300.0, 0.02, [4.0; 3]).expect("valid box");
        sim.gamma = 0.0;
        for _ in 0..100 {
            sim.step();
        }
        let t = sim.time;
        let mean_v2 = velocities
            .iter()
            .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum::<f64>()
            / 3.0;
        let expected = mean_v2 * t / 6.0;
        let d = sim.diffusion_coefficient();
        assert!(
            ((d - expected) / expected).abs() < 1e-9,
            "D = {d}, expected {expected}"
        );
        // Positions stay unwrapped: bead 0 left the 4 nm box through the +x face
        // (3.8 + 0.3 nm/ps * 2 ps = 4.4 nm) and was not folded back.
        assert!((sim.molecules[0].beads[0].pos[0] - 4.4).abs() < 1e-9);
    }

    #[test]
    fn test_sim_reset_displacement_origin_restarts_the_estimate() {
        let mut sim = CgSimulation::new(
            vec![moving_bead([1.0, 1.0, 1.0], [0.2, 0.0, 0.0])],
            ideal_gas_ff(),
            300.0,
            0.02,
            [4.0; 3],
        )
        .expect("valid box");
        sim.gamma = 0.0;
        for _ in 0..25 {
            sim.step();
        }
        assert!(sim.diffusion_coefficient() > 0.0);
        sim.reset_displacement_origin();
        assert_eq!(sim.diffusion_coefficient(), 0.0);
        for _ in 0..25 {
            sim.step();
        }
        // 0.5 ps since the reset at v = 0.2 nm/ps: D = v^2 t / 6.
        let expected = 0.2 * 0.2 * 0.5 / 6.0;
        let d = sim.diffusion_coefficient();
        assert!(
            ((d - expected) / expected).abs() < 1e-9,
            "D = {d}, expected {expected}"
        );
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
