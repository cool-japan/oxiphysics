// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Primitive-model electrolyte MD driver of the electrochemical module.
//!
//! [`ElectrochemicalMdSimulation`] integrates ions in a periodic box with
//! solvent-screened Coulomb forces, a repulsive WCA core and a seeded Langevin
//! thermostat (BAOAB splitting); units are Angstrom, femtosecond, amu and eV.

use super::{
    AVOGADRO, BOLTZMANN, ConstantPotentialMd, ELEMENTARY_CHARGE, EPSILON_0, ElectrochemicalIon,
    GouyChapmanSternModel, IonicSpecies, LithiumIntercalation, MarcusElectronTransfer,
};
use crate::MdError;
use crate::gaussian_noise::standard_normal;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// `e² / (4π ε₀)` in eV·Å: the Coulomb energy of two elementary charges 1 Å apart
/// in vacuum (14.3996 eV·Å).
const COULOMB_EV_ANGSTROM: f64 = ELEMENTARY_CHARGE * 1.0e10 / (4.0 * PI * EPSILON_0);

/// Converts `F / m` in eV Å⁻¹ amu⁻¹ to an acceleration in Å fs⁻² (equally,
/// `k_B T / m` in eV amu⁻¹ to a squared velocity in Å² fs⁻²).
///
/// 1 eV Å⁻¹ amu⁻¹ = e · N_A · 1e13 m s⁻² = 9.648533e-3 Å fs⁻² with
/// 1 amu = 1 g mol⁻¹ / N_A (the molar-mass convention of [`IonicSpecies::mass`]).
const EV_PER_AMU_ANGSTROM_TO_ANGSTROM_PER_FS2: f64 = ELEMENTARY_CHARGE * AVOGADRO * 1.0e-7;

/// Default well depth of the repulsive WCA ion core: `k_B · 298.15 K` in eV.
const DEFAULT_WCA_EPSILON_EV: f64 = BOLTZMANN * 298.15 / ELEMENTARY_CHARGE;

/// Default Langevin friction of [`ElectrochemicalMdSimulation::step`] in fs⁻¹.
const DEFAULT_FRICTION_PER_FS: f64 = 0.01;

/// Edge of the cubic box built by [`ElectrochemicalMdSimulation::new_liion_interface`] in Å.
const LIION_INTERFACE_BOX_ANGSTROM: f64 = 40.0;

/// Seed used by the constructors that take no explicit seed.
const DEFAULT_SEED: u64 = 0x00E1_EC70_5EED;

/// Random trial positions per ion before the initial placement settles for the
/// least-overlapping candidate.
const PLACEMENT_ATTEMPTS: usize = 2000;

// ---------------------------------------------------------------------------
// ElectrochemicalMdSimulation
// ---------------------------------------------------------------------------

/// Inputs the cached ion forces were computed from; any change invalidates the cache.
#[derive(Debug, Clone, PartialEq)]
struct ForceInputs {
    /// `[x, y, z, charge, radius]` of every ion.
    ions: Vec<[f64; 5]>,
    /// Box edges \[Lx, Ly, Lz\] in Angstrom.
    box_size: [f64; 3],
    /// Dielectric constant of the solvent.
    dielectric: f64,
    /// Well depth of the WCA ion core in eV.
    wca_epsilon: f64,
}

/// Interaction of one ion pair (minimum image).
#[derive(Debug, Clone, Copy)]
struct PairInteraction {
    /// Separation `r_i - r_j` wrapped to the minimum image, in Angstrom.
    separation: [f64; 3],
    /// Shifted-force screened Coulomb energy in eV.
    coulomb_energy: f64,
    /// WCA core repulsion energy in eV.
    core_energy: f64,
    /// `F / r` with `F` the radial force in eV/Angstrom (positive = repulsive):
    /// the force on ion `i` is `force_over_r * separation`, on ion `j` the opposite.
    force_over_r: f64,
}

/// High-level driver for electrochemical MD simulations.
///
/// Combines ions, electrode, GCS model, Marcus kinetics, and Li intercalation.
///
/// The ions are a primitive-model electrolyte in an implicit solvent: they move in
/// the periodic box `box_size` under
///
/// * the Coulomb interaction `e² z_i z_j / (4π ε₀ ε r)` screened by the solvent
///   dielectric constant `dielectric`, truncated at half the shortest box edge with
///   the shifted-force scheme (`V(r) = c (1/r + r/r_c² - 2/r_c)`, so that energy
///   and force vanish continuously at the cutoff `r_c`), and
/// * a purely repulsive WCA (Weeks-Chandler-Andersen) core between ions that
///   closes at the contact distance `r_i + r_j` of their ionic radii,
///
/// integrated with the BAOAB Langevin splitting (velocity Verlet for `friction = 0`)
/// and noise from a seeded generator, so that a run is reproducible.
#[derive(Debug, Clone)]
pub struct ElectrochemicalMdSimulation {
    /// Ions in the simulation box.
    pub ions: Vec<ElectrochemicalIon>,
    /// GCS double layer model.
    pub gcs: GouyChapmanSternModel,
    /// CPM electrode.
    pub cpm: ConstantPotentialMd,
    /// Marcus model for electron transfer.
    pub marcus: MarcusElectronTransfer,
    /// Li intercalation model.
    pub intercalation: LithiumIntercalation,
    /// Box size \[Lx, Ly, Lz\] in Angstrom.
    pub box_size: [f64; 3],
    /// Current simulation time in fs.
    pub time: f64,
    /// Integration timestep in fs.
    pub dt: f64,
    /// Temperature in K.
    pub temperature: f64,
    /// Dielectric constant of the solvent.
    pub dielectric: f64,
    /// Langevin friction γ in fs⁻¹ used by [`ElectrochemicalMdSimulation::step`]
    /// (`0` gives microcanonical dynamics).
    pub friction: f64,
    /// Well depth of the repulsive WCA ion core in eV.
    pub wca_epsilon: f64,
    /// Seeded generator behind the thermostat noise.
    rng: SmallRng,
    /// Inputs of the forces currently stored in the ions, if they are up to date.
    forces_valid_for: Option<ForceInputs>,
}

impl ElectrochemicalMdSimulation {
    /// Create a simulation for a Li-ion battery interface.
    ///
    /// Equivalent to [`ElectrochemicalMdSimulation::new_liion_interface_seeded`] with a
    /// fixed default seed, so repeated calls build the same system.
    pub fn new_liion_interface(voltage: f64, temperature: f64, n_ions: usize) -> Self {
        Self::new_liion_interface_seeded(voltage, temperature, n_ions, DEFAULT_SEED)
    }

    /// Create a Li-ion battery interface in a 40 Å cubic box from the generator
    /// seeded with `seed`.
    ///
    /// Li⁺ and PF₆⁻ ions alternate; they are placed at random without core overlap
    /// and receive Maxwell-Boltzmann velocities at `temperature`.  The same seed
    /// reproduces the initial state and the whole trajectory.
    pub fn new_liion_interface_seeded(
        voltage: f64,
        temperature: f64,
        n_ions: usize,
        seed: u64,
    ) -> Self {
        Self::build_liion_interface(
            voltage,
            temperature,
            n_ions,
            LIION_INTERFACE_BOX_ANGSTROM,
            seed,
        )
    }

    /// Like [`ElectrochemicalMdSimulation::new_liion_interface_seeded`] in a cubic box
    /// of edge `box_len` (Å).
    ///
    /// # Errors
    /// Returns [`MdError::General`] if `box_len` is not finite and positive.
    pub fn new_liion_interface_in_box(
        voltage: f64,
        temperature: f64,
        n_ions: usize,
        box_len: f64,
        seed: u64,
    ) -> Result<Self, MdError> {
        if !box_len.is_finite() || box_len <= 0.0 {
            return Err(MdError::General(format!(
                "box edge must be finite and positive, got {box_len}"
            )));
        }
        Ok(Self::build_liion_interface(
            voltage,
            temperature,
            n_ions,
            box_len,
            seed,
        ))
    }

    /// Shared constructor; `box_len` must be finite and positive.
    fn build_liion_interface(
        voltage: f64,
        temperature: f64,
        n_ions: usize,
        box_len: f64,
        seed: u64,
    ) -> Self {
        let box_size = [box_len; 3];
        let mut rng = SmallRng::seed_from_u64(seed);
        let conc = 1.0; // mol/L
        let mut ions: Vec<ElectrochemicalIon> = Vec::with_capacity(n_ions);
        for i in 0..n_ions {
            let species = if i % 2 == 0 {
                IonicSpecies::Li
            } else {
                IonicSpecies::Pf6
            };
            let pos = Self::sample_free_position(&ions, species, box_size, &mut rng);
            let mut ion = ElectrochemicalIon::new(pos, species);
            ion.id = i;
            ions.push(ion);
        }
        // Maxwell-Boltzmann velocities: sigma_v^2 = k_B T / m.
        let kbt_ev = (BOLTZMANN * temperature / ELEMENTARY_CHARGE).max(0.0);
        for ion in &mut ions {
            let sigma_v =
                (kbt_ev / ion.mass_amu() * EV_PER_AMU_ANGSTROM_TO_ANGSTROM_PER_FS2).sqrt();
            for component in &mut ion.vel {
                *component = sigma_v * standard_normal(&mut rng);
            }
        }
        Self {
            ions,
            gcs: GouyChapmanSternModel::aqueous_kcl(conc, voltage, temperature),
            cpm: ConstantPotentialMd::graphene_electrode(50, voltage, temperature),
            marcus: MarcusElectronTransfer::li_couple(voltage, temperature),
            intercalation: LithiumIntercalation::graphite_anode(temperature),
            box_size,
            time: 0.0,
            dt: 1.0,
            temperature,
            dielectric: 78.5,
            friction: DEFAULT_FRICTION_PER_FS,
            wca_epsilon: DEFAULT_WCA_EPSILON_EV,
            rng,
            forces_valid_for: None,
        }
    }

    /// Random position for a new ion of `species` that keeps the contact distance
    /// (sum of ionic radii) to every ion in `placed`; if no trial position is free
    /// within [`PLACEMENT_ATTEMPTS`] attempts, the least-overlapping one is used.
    fn sample_free_position(
        placed: &[ElectrochemicalIon],
        species: IonicSpecies,
        box_size: [f64; 3],
        rng: &mut SmallRng,
    ) -> [f64; 3] {
        let mut best = [0.0; 3];
        let mut best_clearance = f64::NEG_INFINITY;
        for _ in 0..PLACEMENT_ATTEMPTS {
            let candidate = [
                rng.random_range(0.0..box_size[0]),
                rng.random_range(0.0..box_size[1]),
                rng.random_range(0.0..box_size[2]),
            ];
            // Ratio of the separation to the contact distance, smallest over all ions.
            let clearance = placed
                .iter()
                .map(|other| {
                    let d = minimum_image_vector(candidate, other.pos, box_size);
                    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
                        / (species.ionic_radius() + other.radius())
                })
                .fold(f64::INFINITY, f64::min);
            if clearance > best_clearance {
                best_clearance = clearance;
                best = candidate;
            }
            if clearance >= 1.0 {
                break;
            }
        }
        best
    }

    /// Cutoff of the Coulomb interaction in Angstrom: half the shortest box edge, the
    /// largest range for which the minimum image is complete.
    pub fn coulomb_cutoff(&self) -> f64 {
        0.5 * self.box_size.iter().copied().fold(f64::INFINITY, f64::min)
    }

    /// Interaction of the ion pair `(i, j)`, or `None` if they are outside both the
    /// Coulomb cutoff and the WCA core (or coincide).
    fn pair_interaction(&self, i: usize, j: usize, cutoff: f64) -> Option<PairInteraction> {
        let ion_i = &self.ions[i];
        let ion_j = &self.ions[j];
        let separation = minimum_image_vector(ion_i.pos, ion_j.pos, self.box_size);
        let r_sq = separation[0] * separation[0]
            + separation[1] * separation[1]
            + separation[2] * separation[2];
        if r_sq < 1.0e-12 {
            return None;
        }
        let r = r_sq.sqrt();
        let mut coulomb_energy = 0.0;
        let mut core_energy = 0.0;
        let mut force_over_r = 0.0;

        let charge_product = ion_i.charge * ion_j.charge;
        if r < cutoff && charge_product != 0.0 {
            // Shifted-force Coulomb: V = c (1/r + r/r_c² - 2/r_c),
            // F = -dV/dr = c (1/r² - 1/r_c²).
            let c = COULOMB_EV_ANGSTROM * charge_product / self.dielectric;
            let cutoff_sq = cutoff * cutoff;
            coulomb_energy = c * (1.0 / r + r / cutoff_sq - 2.0 / cutoff);
            force_over_r += c * (1.0 / r_sq - 1.0 / cutoff_sq) / r;
        }

        let contact = ion_i.radius() + ion_j.radius();
        if r < contact {
            // WCA: LJ with sigma = contact / 2^(1/6) (minimum at the contact distance),
            // cut at the minimum and shifted by eps.
            let sigma_sq = contact * contact / 2.0_f64.powf(1.0 / 3.0);
            let sr2 = sigma_sq / r_sq;
            let sr6 = sr2 * sr2 * sr2;
            core_energy = 4.0 * self.wca_epsilon * (sr6 * sr6 - sr6) + self.wca_epsilon;
            force_over_r += 24.0 * self.wca_epsilon * (2.0 * sr6 * sr6 - sr6) / r_sq;
        }

        if coulomb_energy == 0.0 && core_energy == 0.0 && force_over_r == 0.0 {
            return None;
        }
        Some(PairInteraction {
            separation,
            coulomb_energy,
            core_energy,
            force_over_r,
        })
    }

    /// Sum of the pair Coulomb and WCA core energies in eV.
    fn interaction_energies(&self) -> (f64, f64) {
        let n = self.ions.len();
        let cutoff = self.coulomb_cutoff();
        let mut coulomb = 0.0;
        let mut core = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(pair) = self.pair_interaction(i, j, cutoff) {
                    coulomb += pair.coulomb_energy;
                    core += pair.core_energy;
                }
            }
        }
        (coulomb, core)
    }

    /// Compute total Coulomb energy in eV.
    ///
    /// Sum over ion pairs of the minimum-image, shifted-force screened Coulomb energy
    /// that [`ElectrochemicalMdSimulation::compute_forces`] differentiates.
    pub fn total_coulomb_energy(&self) -> f64 {
        self.interaction_energies().0
    }

    /// Total potential energy in eV: screened Coulomb plus WCA core repulsion.
    pub fn total_potential_energy_ev(&self) -> f64 {
        let (coulomb, core) = self.interaction_energies();
        coulomb + core
    }

    /// Total energy (kinetic + potential) in eV; conserved by the integrator when
    /// `friction` is zero.
    pub fn total_energy_ev(&self) -> f64 {
        self.total_kinetic_energy_ev() + self.total_potential_energy_ev()
    }

    /// Inputs that determine the stored ion forces.
    fn force_inputs(&self) -> ForceInputs {
        ForceInputs {
            ions: self
                .ions
                .iter()
                .map(|ion| [ion.pos[0], ion.pos[1], ion.pos[2], ion.charge, ion.radius()])
                .collect(),
            box_size: self.box_size,
            dielectric: self.dielectric,
            wca_epsilon: self.wca_epsilon,
        }
    }

    /// Recompute the force on every ion in eV/Angstrom and return the total
    /// potential energy in eV.
    ///
    /// Forces are the exact negative gradient of
    /// [`ElectrochemicalMdSimulation::total_potential_energy_ev`]: the shifted-force
    /// screened Coulomb interaction `F = e² z_i z_j / (4π ε₀ ε) (1/r² - 1/r_c²)`
    /// (14.3996 eV·Å · z_i z_j / (ε r²) at short range) and the WCA core repulsion,
    /// both with minimum-image separations in `box_size`.
    pub fn compute_forces(&mut self) -> f64 {
        let n = self.ions.len();
        let cutoff = self.coulomb_cutoff();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut potential = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let Some(pair) = self.pair_interaction(i, j, cutoff) else {
                    continue;
                };
                potential += pair.coulomb_energy + pair.core_energy;
                let pair_force = pair
                    .separation
                    .map(|component| pair.force_over_r * component);
                for (acc, f) in forces[i].iter_mut().zip(pair_force.iter()) {
                    *acc += f;
                }
                for (acc, f) in forces[j].iter_mut().zip(pair_force.iter()) {
                    *acc -= f;
                }
            }
        }
        for (ion, force) in self.ions.iter_mut().zip(forces.iter()) {
            ion.force = *force;
        }
        self.forces_valid_for = Some(self.force_inputs());
        potential
    }

    /// Apply the Ornstein-Uhlenbeck (Langevin) velocity update for one timestep.
    ///
    /// With friction `gamma` (fs⁻¹) every velocity component becomes
    /// `v <- exp(-gamma dt) v + sqrt((1 - exp(-2 gamma dt)) k_B T / m) xi`,
    /// where `xi` is a standard normal variate from the seeded generator; this is the
    /// exact solution of the friction + white-noise part of the Langevin equation and
    /// leaves the Maxwell-Boltzmann distribution at `temperature` invariant.
    /// A negative or NaN `gamma` leaves the velocities unchanged.
    pub fn apply_langevin(&mut self, gamma: f64) {
        let decay = (-gamma.max(0.0) * self.dt).exp();
        let noise_scale = (1.0 - decay * decay).max(0.0).sqrt();
        let kbt_ev = (BOLTZMANN * self.temperature / ELEMENTARY_CHARGE).max(0.0);
        for ion in &mut self.ions {
            let thermal_speed =
                (kbt_ev / ion.mass_amu() * EV_PER_AMU_ANGSTROM_TO_ANGSTROM_PER_FS2).sqrt();
            for component in &mut ion.vel {
                *component = decay * *component
                    + noise_scale * thermal_speed * standard_normal(&mut self.rng);
            }
        }
    }

    /// Velocity half-kick `v += (dt / 2) F / m` with the stored forces.
    fn half_kick(&mut self) {
        let half_dt = 0.5 * self.dt;
        for ion in &mut self.ions {
            let factor = EV_PER_AMU_ANGSTROM_TO_ANGSTROM_PER_FS2 / ion.mass_amu() * half_dt;
            for axis in 0..3 {
                ion.vel[axis] += ion.force[axis] * factor;
            }
        }
    }

    /// Position half-drift `x += (dt / 2) v`.
    fn half_drift(&mut self) {
        let half_dt = 0.5 * self.dt;
        for ion in &mut self.ions {
            for axis in 0..3 {
                ion.pos[axis] += ion.vel[axis] * half_dt;
            }
        }
    }

    /// Fold every ion back into the periodic box `[0, L)`.
    fn wrap_positions(&mut self) {
        for ion in &mut self.ions {
            for axis in 0..3 {
                let edge = self.box_size[axis];
                let wrapped = ion.pos[axis].rem_euclid(edge);
                // rem_euclid can round up to exactly `edge` for tiny negative inputs.
                ion.pos[axis] = if wrapped >= edge { 0.0 } else { wrapped };
            }
        }
    }

    /// Advance by one timestep with the BAOAB Langevin splitting.
    ///
    /// half kick (B), half drift (A), Ornstein-Uhlenbeck thermostat with
    /// `friction` (O), half drift (A), force update, half kick (B).  For
    /// `friction = 0` this is the velocity-Verlet integrator.  Forces are reused from
    /// the previous step when the configuration has not been modified in between and
    /// recomputed otherwise.
    pub fn step(&mut self) {
        if self.forces_valid_for.as_ref() != Some(&self.force_inputs()) {
            self.compute_forces();
        }
        self.half_kick();
        self.half_drift();
        self.apply_langevin(self.friction);
        self.half_drift();
        self.wrap_positions();
        self.compute_forces();
        self.half_kick();
        self.time += self.dt;
    }

    /// Run for n_steps timesteps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Total kinetic energy in eV.
    pub fn total_kinetic_energy_ev(&self) -> f64 {
        self.ions.iter().map(|i| i.kinetic_energy_ev()).sum()
    }

    /// Temperature estimate from kinetic energy in K.
    pub fn temperature_estimate(&self) -> f64 {
        let ke = self.total_kinetic_energy_ev() * ELEMENTARY_CHARGE; // J
        let n = self.ions.len() as f64;
        if n < 1.0 {
            return 0.0;
        }
        2.0 * ke / (3.0 * n * BOLTZMANN)
    }

    /// Count cation-anion pairs within cutoff in Angstrom.
    pub fn count_ion_pairs(&self, cutoff: f64) -> usize {
        let n = self.ions.len();
        let mut count = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let zi = self.ions[i].species.charge_number();
                let zj = self.ions[j].species.charge_number();
                if zi * zj < 0 && self.ions[i].distance_to(&self.ions[j]) < cutoff {
                    count += 1;
                }
            }
        }
        count
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Separation `a - b` wrapped to the minimum image of an orthorhombic periodic box.
fn minimum_image_vector(a: [f64; 3], b: [f64; 3], box_size: [f64; 3]) -> [f64; 3] {
    let mut d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    for axis in 0..3 {
        d[axis] -= box_size[axis] * (d[axis] / box_size[axis]).round();
    }
    d
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::super::STD_TEMP;
    use super::*;

    // --- ElectrochemicalMdSimulation ---

    #[test]
    fn test_sim_ion_count() {
        let sim = ElectrochemicalMdSimulation::new_liion_interface(1.0, STD_TEMP, 10);
        assert_eq!(sim.ions.len(), 10);
    }

    #[test]
    fn test_sim_time_advances() {
        let mut sim = ElectrochemicalMdSimulation::new_liion_interface(1.0, STD_TEMP, 4);
        sim.run(5);
        assert!(sim.time > 0.0);
    }

    #[test]
    fn test_sim_ions_stay_in_box() {
        let mut sim = ElectrochemicalMdSimulation::new_liion_interface(1.0, STD_TEMP, 8);
        sim.run(10);
        for ion in &sim.ions {
            for k in 0..3 {
                assert!(ion.pos[k] >= 0.0 && ion.pos[k] < sim.box_size[k]);
            }
        }
    }

    // --- ElectrochemicalMdSimulation: interparticle forces and Langevin dynamics ---

    /// Energy/force conversion used by the checks below, written out independently
    /// of the module constants: 1 eV/(Å amu) = 9.648533e-3 Å/fs².
    const ACCEL_FACTOR: f64 = 9.648_533e-3;

    /// An empty simulation in a cubic box of edge `box_len` (Å) at `temperature` whose
    /// ions are supplied by the caller.
    fn simulation_with_ions(
        ions: Vec<ElectrochemicalIon>,
        box_len: f64,
        temperature: f64,
    ) -> ElectrochemicalMdSimulation {
        let mut sim = ElectrochemicalMdSimulation::new_liion_interface_in_box(
            1.0,
            temperature,
            0,
            box_len,
            11,
        )
        .expect("valid box");
        sim.ions = ions;
        sim
    }

    fn ion_at(species: IonicSpecies, pos: [f64; 3]) -> ElectrochemicalIon {
        ElectrochemicalIon::new(pos, species)
    }

    /// Closed-form shifted-force Coulomb force (eV/Å, positive = repulsive) between
    /// charges `qq` (in e²) at distance `r`, for dielectric `eps` and cutoff `rc`.
    fn shifted_force(qq: f64, eps: f64, r: f64, rc: f64) -> f64 {
        14.399_645 * qq / eps * (1.0 / (r * r) - 1.0 / (rc * rc))
    }

    #[test]
    fn test_sim_opposite_ions_attract() {
        // Li+ at x = 10 Å, Cl- at x = 18 Å in a 40 Å box (cutoff 20 Å).
        let ions = vec![
            ion_at(IonicSpecies::Li, [10.0, 20.0, 20.0]),
            ion_at(IonicSpecies::Cl, [18.0, 20.0, 20.0]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.compute_forces();
        let expected = shifted_force(-1.0, 78.5, 8.0, 20.0).abs();
        // The attraction pulls Li+ towards +x (the Cl- side) with |F| = c (1/r² - 1/r_c²).
        assert!(
            (sim.ions[0].force[0] - expected).abs() < 1e-7 * expected,
            "F_Li = {:?}, expected {expected}",
            sim.ions[0].force
        );
        assert!((sim.ions[1].force[0] + expected).abs() < 1e-7 * expected);
        for axis in 1..3 {
            assert!(sim.ions[0].force[axis].abs() < 1e-15);
            assert!(sim.ions[1].force[axis].abs() < 1e-15);
        }
        // Without friction the pair accelerates together and the distance shrinks.
        sim.friction = 0.0;
        sim.ions[0].vel = [0.0; 3];
        sim.ions[1].vel = [0.0; 3];
        sim.run(200);
        let gap = sim.ions[1].pos[0] - sim.ions[0].pos[0];
        assert!(
            gap < 7.95,
            "opposite ions must approach each other: gap = {gap}"
        );
        assert!(gap > 7.5, "unexpectedly large displacement: gap = {gap}");
    }

    #[test]
    fn test_sim_like_ions_repel() {
        let ions = vec![
            ion_at(IonicSpecies::Na, [10.0, 20.0, 20.0]),
            ion_at(IonicSpecies::K, [18.0, 20.0, 20.0]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.compute_forces();
        let expected = shifted_force(1.0, 78.5, 8.0, 20.0);
        assert!(expected > 0.0);
        assert!(
            (sim.ions[0].force[0] + expected).abs() < 1e-7 * expected,
            "{:?}",
            sim.ions[0].force
        );
        assert!((sim.ions[1].force[0] - expected).abs() < 1e-7 * expected);
        // Divalent calcium repels a monovalent cation twice as strongly.
        sim.ions[1] = ion_at(IonicSpecies::Ca, [18.0, 20.0, 20.0]);
        sim.compute_forces();
        assert!((sim.ions[1].force[0] / expected - 2.0).abs() < 1e-7);
    }

    #[test]
    fn test_sim_coulomb_force_vanishes_at_and_beyond_the_cutoff() {
        let rc = 20.0;
        // Just inside the cutoff the shifted force is tiny but positive-definite...
        let inside = vec![
            ion_at(IonicSpecies::Na, [0.0, 0.0, 0.0]),
            ion_at(IonicSpecies::K, [rc - 1e-3, 0.0, 0.0]),
        ];
        let mut sim = simulation_with_ions(inside, 40.0, 300.0);
        sim.compute_forces();
        let f_inside = sim.ions[1].force[0];
        let expected = shifted_force(1.0, 78.5, rc - 1e-3, rc);
        assert!(expected > 0.0 && expected < 1e-6, "closed form {expected}");
        assert!(
            (f_inside - expected).abs() < 1e-3 * expected,
            "force just inside r_c: {f_inside}, expected {expected}"
        );
        // ...and a separation along the box diagonal (34.6 Å) lies outside the cutoff.
        let outside = vec![
            ion_at(IonicSpecies::Na, [0.0, 0.0, 0.0]),
            ion_at(IonicSpecies::K, [20.0, 20.0, 20.0]),
        ];
        let mut sim = simulation_with_ions(outside, 40.0, 300.0);
        let potential = sim.compute_forces();
        assert_eq!(potential, 0.0);
        assert_eq!(sim.ions[0].force, [0.0; 3]);
        assert_eq!(sim.total_coulomb_energy(), 0.0);
    }

    #[test]
    fn test_sim_forces_are_the_negative_gradient_of_the_potential_energy() {
        // Six ions of mixed charge; Li+/PF6- and Na+/Cl- start inside their contact
        // distance (3.3 Å and 2.83 Å) so that the WCA core is active, and the last
        // pair straddles the periodic boundary.
        let ions = vec![
            ion_at(IonicSpecies::Li, [10.0, 10.0, 10.0]),
            ion_at(IonicSpecies::Pf6, [12.4, 11.6, 10.3]),
            ion_at(IonicSpecies::Na, [25.0, 8.0, 14.0]),
            ion_at(IonicSpecies::Cl, [26.9, 9.6, 14.9]),
            ion_at(IonicSpecies::Ca, [39.0, 30.0, 5.0]),
            ion_at(IonicSpecies::Tfsi, [2.5, 31.0, 6.5]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        let potential = sim.compute_forces();
        assert!((potential - sim.total_potential_energy_ev()).abs() < 1e-12);
        assert!(potential != 0.0);
        let analytic: Vec<[f64; 3]> = sim.ions.iter().map(|ion| ion.force).collect();
        let h = 1e-6;
        for (atom, analytic_force) in analytic.iter().enumerate() {
            for (axis, &analytic_component) in analytic_force.iter().enumerate() {
                let origin = sim.ions[atom].pos[axis];
                sim.ions[atom].pos[axis] = origin + h;
                let e_plus = sim.total_potential_energy_ev();
                sim.ions[atom].pos[axis] = origin - h;
                let e_minus = sim.total_potential_energy_ev();
                sim.ions[atom].pos[axis] = origin;
                let numeric = -(e_plus - e_minus) / (2.0 * h);
                assert!(
                    (analytic_component - numeric).abs() < 1e-6 * numeric.abs().max(1e-3),
                    "ion {atom} axis {axis}: analytic {analytic_component} vs numeric {numeric}"
                );
            }
        }
        // Newton's third law.
        for axis in 0..3 {
            let net: f64 = analytic.iter().map(|f| f[axis]).sum();
            assert!(net.abs() < 1e-12, "net force {net} on axis {axis}");
        }
    }

    #[test]
    fn test_sim_wca_core_repels_ions_inside_contact_distance() {
        // Two like-charged Li+ ions 0.6 Å apart (contact 1.52 Å): strong core repulsion
        // on top of the weak Coulomb repulsion.
        let ions = vec![
            ion_at(IonicSpecies::Li, [10.0, 10.0, 10.0]),
            ion_at(IonicSpecies::Li, [10.6, 10.0, 10.0]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.compute_forces();
        let coulomb_only = shifted_force(1.0, 78.5, 0.6, 20.0);
        assert!(sim.ions[1].force[0] > 100.0 * coulomb_only);
        assert!(sim.ions[0].force[0] < 0.0);
        // At the contact distance the core exerts no force, so a Na+ / Cl- pair at
        // 2.83 Å feels the Coulomb attraction alone.
        let ions = vec![
            ion_at(IonicSpecies::Na, [10.0, 10.0, 10.0]),
            ion_at(IonicSpecies::Cl, [12.83, 10.0, 10.0]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.compute_forces();
        let attraction = shifted_force(-1.0, 78.5, 2.83, 20.0).abs();
        assert!((sim.ions[0].force[0] - attraction).abs() < 1e-6 * attraction);
    }

    #[test]
    fn test_sim_nve_conserves_energy_and_momentum() {
        // 16 ions in a 14 Å box (dense enough for many encounters), no thermostat.
        let mut sim =
            ElectrochemicalMdSimulation::new_liion_interface_in_box(1.0, 300.0, 16, 14.0, 5)
                .expect("valid box");
        sim.friction = 0.0;
        sim.dt = 0.5;
        sim.compute_forces();
        let e0 = sim.total_energy_ev();
        let ke0 = sim.total_kinetic_energy_ev();
        let momentum = |sim: &ElectrochemicalMdSimulation| -> [f64; 3] {
            let mut p = [0.0; 3];
            for ion in &sim.ions {
                for (component, velocity) in p.iter_mut().zip(ion.vel.iter()) {
                    *component += ion.mass_amu() * velocity;
                }
            }
            p
        };
        let p0 = momentum(&sim);
        let mut max_drift = 0.0_f64;
        for _ in 0..1000 {
            sim.step();
            max_drift = max_drift.max((sim.total_energy_ev() - e0).abs());
        }
        assert!(
            max_drift < 1e-4 * ke0,
            "energy drift {max_drift} eV against kinetic energy {ke0} eV"
        );
        // The ions interact: the potential energy is not negligible.
        assert!(sim.total_potential_energy_ev().abs() > 1e-3 * ke0);
        let p1 = momentum(&sim);
        let p_scale: f64 = sim.ions.iter().map(|i| i.mass_amu() * i.speed()).sum();
        for axis in 0..3 {
            assert!(
                (p1[axis] - p0[axis]).abs() < 1e-9 * p_scale,
                "momentum drifted on axis {axis}: {} -> {}",
                p0[axis],
                p1[axis]
            );
        }
    }

    #[test]
    fn test_sim_langevin_thermostat_reaches_the_target_temperature() {
        // Start from rest: only the friction + noise can heat the ions to 300 K.
        let mut sim =
            ElectrochemicalMdSimulation::new_liion_interface_in_box(1.0, 300.0, 40, 40.0, 3)
                .expect("valid box");
        for ion in &mut sim.ions {
            ion.vel = [0.0; 3];
        }
        sim.friction = 0.02;
        assert_eq!(sim.temperature_estimate(), 0.0);
        sim.run(2000);
        let mut sum = 0.0;
        let mut samples = 0;
        for block in 0..400 {
            sim.run(20);
            sum += sim.temperature_estimate();
            samples += 1;
            assert!(sim.temperature_estimate() > 0.0, "block {block}");
        }
        let mean_t = sum / f64::from(samples);
        assert!(
            (mean_t / 300.0 - 1.0).abs() < 0.05,
            "mean temperature {mean_t} K, target 300 K"
        );
    }

    #[test]
    fn test_sim_langevin_noise_is_gaussian_with_equipartition_variance() {
        // gamma * dt = 1000: the old velocity is forgotten and each component is a
        // fresh draw with variance k_B T / m (in Å²/fs²), standard normal in shape.
        let n = 4000;
        let ions: Vec<ElectrochemicalIon> = (0..n)
            .map(|_| ion_at(IonicSpecies::Na, [1.0, 2.0, 3.0]))
            .collect();
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.dt = 1.0;
        for ion in &mut sim.ions {
            ion.vel = [5.0, -5.0, 5.0];
        }
        sim.apply_langevin(1.0e3);
        let samples: Vec<f64> = sim.ions.iter().flat_map(|ion| ion.vel).collect();
        let count = samples.len() as f64;
        let mean = samples.iter().sum::<f64>() / count;
        let variance = samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / count;
        let fourth = samples.iter().map(|v| (v - mean).powi(4)).sum::<f64>() / count;
        let kbt_ev = 1.380_649e-23 * 300.0 / 1.602_176_634e-19;
        let expected_variance = kbt_ev / 22.990 * ACCEL_FACTOR;
        assert!(
            (variance / expected_variance - 1.0).abs() < 0.05,
            "variance {variance}, expected k_B T / m = {expected_variance}"
        );
        assert!(
            mean.abs() < 4.0 * (expected_variance / count).sqrt(),
            "mean {mean}"
        );
        // Kurtosis 3 for a Gaussian (1.8 for the former uniform noise).
        let kurtosis = fourth / (variance * variance);
        assert!((kurtosis - 3.0).abs() < 0.3, "kurtosis {kurtosis}");
    }

    #[test]
    fn test_sim_zero_friction_leaves_velocities_untouched() {
        let mut sim = ElectrochemicalMdSimulation::new_liion_interface_seeded(1.0, 300.0, 6, 9);
        let before: Vec<[f64; 3]> = sim.ions.iter().map(|ion| ion.vel).collect();
        sim.apply_langevin(0.0);
        sim.apply_langevin(-1.0);
        sim.apply_langevin(f64::NAN);
        let after: Vec<[f64; 3]> = sim.ions.iter().map(|ion| ion.vel).collect();
        assert_eq!(before, after);
    }

    #[test]
    fn test_sim_acceleration_uses_the_ev_per_amu_angstrom_conversion() {
        // Li+ / Cl- at 8 Å, at rest, friction 0, tiny step: the velocity gained by Li+
        // is F / m_Li in Å/fs² times dt, with the 9.648533e-3 factor of
        // eV Å⁻¹ amu⁻¹ -> Å fs⁻².
        let ions = vec![
            ion_at(IonicSpecies::Li, [10.0, 20.0, 20.0]),
            ion_at(IonicSpecies::Cl, [18.0, 20.0, 20.0]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.friction = 0.0;
        sim.dt = 0.1;
        sim.step();
        let force = shifted_force(-1.0, 78.5, 8.0, 20.0).abs();
        let expected_li = force / 6.941 * ACCEL_FACTOR * 0.1;
        let expected_cl = -force / 35.453 * ACCEL_FACTOR * 0.1;
        assert!(
            (sim.ions[0].vel[0] / expected_li - 1.0).abs() < 1e-3,
            "v_Li = {}, expected {expected_li}",
            sim.ions[0].vel[0]
        );
        assert!((sim.ions[1].vel[0] / expected_cl - 1.0).abs() < 1e-3);
        assert!((sim.time - 0.1).abs() < 1e-15);
    }

    #[test]
    fn test_sim_forces_are_recomputed_after_the_configuration_is_edited() {
        // After a step the forces are cached.  Flipping a charge between steps must not
        // be answered with the stale cache: the kick uses the new attraction/repulsion.
        let ions = vec![
            ion_at(IonicSpecies::Li, [10.0, 20.0, 20.0]),
            ion_at(IonicSpecies::Cl, [18.0, 20.0, 20.0]),
        ];
        let mut sim = simulation_with_ions(ions, 40.0, 300.0);
        sim.friction = 0.0;
        sim.dt = 0.1;
        sim.step();
        sim.ions[0].vel = [0.0; 3];
        sim.ions[1].vel = [0.0; 3];
        sim.ions[1].charge = 1.0; // Cl- -> a cation: the pair now repels
        sim.step();
        assert!(
            sim.ions[0].vel[0] < 0.0,
            "Li+ must be pushed away from the now positive partner: {:?}",
            sim.ions[0].vel
        );
        let force = shifted_force(1.0, 78.5, 8.0, 20.0);
        let expected = -force / 6.941 * ACCEL_FACTOR * 0.1;
        assert!((sim.ions[0].vel[0] / expected - 1.0).abs() < 1e-3);
    }

    #[test]
    fn test_sim_seeded_runs_are_reproducible() {
        let run = |seed: u64| {
            let mut sim =
                ElectrochemicalMdSimulation::new_liion_interface_seeded(1.0, 300.0, 8, seed);
            sim.run(60);
            sim.ions.iter().map(|ion| ion.pos).collect::<Vec<_>>()
        };
        assert_eq!(run(42), run(42));
        assert_ne!(run(42), run(43));
        // The unseeded-looking constructor is deterministic as well.
        let a = ElectrochemicalMdSimulation::new_liion_interface(1.0, 300.0, 6);
        let b = ElectrochemicalMdSimulation::new_liion_interface(1.0, 300.0, 6);
        let positions = |sim: &ElectrochemicalMdSimulation| -> Vec<[f64; 3]> {
            sim.ions.iter().map(|ion| ion.pos).collect()
        };
        assert_eq!(positions(&a), positions(&b));
    }

    #[test]
    fn test_sim_initial_ions_do_not_overlap() {
        let sim = ElectrochemicalMdSimulation::new_liion_interface_seeded(1.0, 300.0, 60, 21);
        assert_eq!(sim.ions.len(), 60);
        for i in 0..sim.ions.len() {
            for j in (i + 1)..sim.ions.len() {
                let d = minimum_image_vector(sim.ions[i].pos, sim.ions[j].pos, sim.box_size);
                let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                let contact = sim.ions[i].radius() + sim.ions[j].radius();
                assert!(
                    r >= contact - 1e-12,
                    "ions {i},{j}: r = {r} < contact {contact}"
                );
            }
        }
        // Ions are labelled by their index and start with Maxwell-Boltzmann velocities.
        assert!(sim.ions.iter().enumerate().all(|(idx, ion)| ion.id == idx));
        assert!(sim.temperature_estimate() > 0.0);
    }

    #[test]
    fn test_ion_mass_amu_is_the_molar_mass_in_grams_per_mole() {
        assert!((ion_at(IonicSpecies::Li, [0.0; 3]).mass_amu() - 6.941).abs() < 1e-9);
        assert!((ion_at(IonicSpecies::Pf6, [0.0; 3]).mass_amu() - 144.96).abs() < 1e-9);
        assert_eq!(ion_at(IonicSpecies::Cl, [0.0; 3]).radius(), 1.81);
    }

    #[test]
    fn test_sim_in_box_constructor_rejects_invalid_edge() {
        for bad in [0.0, -5.0, f64::NAN, f64::INFINITY] {
            let result =
                ElectrochemicalMdSimulation::new_liion_interface_in_box(1.0, 300.0, 4, bad, 1);
            assert!(matches!(result, Err(MdError::General(_))), "edge {bad}");
        }
        let sim = ElectrochemicalMdSimulation::new_liion_interface_in_box(1.0, 300.0, 4, 25.0, 1)
            .expect("valid box");
        assert_eq!(sim.box_size, [25.0; 3]);
        assert!((sim.coulomb_cutoff() - 12.5).abs() < 1e-15);
    }
}
