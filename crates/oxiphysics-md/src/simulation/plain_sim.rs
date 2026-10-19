// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Plain (non-generic) MD simulation and observable tracking.
//!
//! Contains MdSimulationPlain, ObservableRecord, and displacement functions.

use super::advanced::ParrinelloRahmanBarostat;
use super::checkpoint::NptBarostat;
use super::config::{KB_REDUCED, velocity_increment};
use super::forces::{Atom, BarostatType, ThermostatType};
use super::plain_thermostats::{NoseHooverChainState, csvr_kinetic_energy};
use super::stochastic::SimRng;
use super::tracking::KJ_MOL_PER_ANG3_TO_BAR;
use crate::error::{Error, Result};

/// Seed of the random-number stream of a freshly created [`MdSimulationPlain`].
const DEFAULT_SEED: u64 = 0x0DDB_1A5E_5BAD_5EED;

/// Build the error returned for an invalid parameter or configuration.
fn invalid(message: String) -> Error {
    Error::General(message)
}

// ---------------------------------------------------------------------------
// ForceEvaluation — result of one force-routine call
// ---------------------------------------------------------------------------

/// Forces, pair virial and potential energy returned by the force routine of
/// [`MdSimulationPlain::step_with`].
///
/// Only the forces are mandatory.  Pressure coupling needs the scalar pair
/// virial and the Monte-Carlo barostat the potential energy; a step that needs
/// a quantity the routine did not supply fails with an error instead of using
/// a made-up value.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ForceEvaluation {
    /// Force on every atom (kJ mol⁻¹ Å⁻¹), in atom order.
    pub forces: Vec<[f64; 3]>,
    /// Scalar pair virial `W = Σ_{i<j} r_ij · F_ij` (kJ mol⁻¹) with
    /// minimum-image separations `r_ij = r_i − r_j` and `F_ij` the force on atom
    /// `i` due to atom `j`; a repulsive pair contributes `+|F| r`.  Do **not**
    /// pass the atomic sum `Σ_i r_i · F_i` of wrapped coordinates: under periodic
    /// boundary conditions it contains a spurious `L · F_ij` term for every pair
    /// that interacts across a cell face.
    pub virial: Option<f64>,
    /// Potential energy (kJ mol⁻¹) of the configuration.
    pub potential_energy: Option<f64>,
}

impl ForceEvaluation {
    /// An evaluation that supplies only the forces.
    pub fn new(forces: Vec<[f64; 3]>) -> Self {
        Self {
            forces,
            virial: None,
            potential_energy: None,
        }
    }

    /// Attach the scalar pair virial `W` (kJ mol⁻¹), see
    /// [`ForceEvaluation::virial`].
    pub fn with_virial(mut self, virial: f64) -> Self {
        self.virial = Some(virial);
        self
    }

    /// Attach the potential energy (kJ mol⁻¹).
    pub fn with_potential_energy(mut self, potential_energy: f64) -> Self {
        self.potential_energy = Some(potential_energy);
        self
    }
}

// ---------------------------------------------------------------------------
// McVolumeMoves — settings and statistics of the Monte-Carlo barostat
// ---------------------------------------------------------------------------

/// Settings and statistics of the Monte-Carlo volume moves of
/// [`BarostatType::MonteCarlo`].
///
/// A trial move draws `ln(V'/V)` uniformly from `[−ln(1+δ), +ln(1+δ)]` with
/// `δ =` [`McVolumeMoves::max_relative_change`], i.e. the volume grows by at most
/// the factor `1 + δ` and shrinks by at most `1/(1 + δ)`.  Every coordinate is
/// scaled by `(V'/V)^(1/3)` and the move is accepted with the probability
///
/// ```text
/// min(1, exp(−[ΔU + P ΔV] / k_B T + (N + 1) ln(V'/V)))
/// ```
///
/// of the isothermal-isobaric ensemble (McDonald 1972; Frenkel & Smit,
/// *Understanding Molecular Simulation*, §5.4).  The proposal is symmetric in
/// `ln V`, which is why the weight carries `N + 1` rather than `N`.  A proposal
/// whose width is proportional to the current volume but uniform in `V` itself is
/// *not* symmetric (the reverse of a large contraction may be impossible), and
/// without a Hastings correction it samples a volume distribution that is too
/// narrow and shifted to smaller volumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct McVolumeMoves {
    /// Largest relative volume change of a trial move, in `(0, 1)`: the volume
    /// ratio `V'/V` of a move lies in `[1/(1 + δ), 1 + δ]`.
    pub max_relative_change: f64,
    /// A move is attempted every `interval` steps (`0` disables the moves).
    pub interval: u64,
    /// Number of attempted moves so far.
    pub attempted: u64,
    /// Number of accepted moves so far.
    pub accepted: u64,
}

impl Default for McVolumeMoves {
    fn default() -> Self {
        Self {
            max_relative_change: 0.01,
            interval: 25,
            attempted: 0,
            accepted: 0,
        }
    }
}

impl McVolumeMoves {
    /// Fraction of the attempted moves that were accepted (`0.0` before the
    /// first attempt).
    pub fn acceptance_rate(&self) -> f64 {
        if self.attempted == 0 {
            0.0
        } else {
            self.accepted as f64 / self.attempted as f64
        }
    }
}

// ---------------------------------------------------------------------------
// MdSimulationPlain — self-contained simulation on Vec<Atom>
// ---------------------------------------------------------------------------

/// Self-contained MD simulation driver operating on a `Vec<Atom>`.
///
/// Works in the MD unit system of the module (Å, ps, amu, kJ mol⁻¹, K, bar; see
/// [`config`](super::config)).  Provides velocity-Verlet integration with
/// per-step force evaluation ([`MdSimulationPlain::step_with`]), PBC wrapping,
/// temperature, kinetic energy and virial pressure, and the ensemble controls
/// selected by the `thermostat` and `barostat` fields:
///
/// | field value                    | algorithm |
/// |--------------------------------|-----------|
/// | `ThermostatType::Berendsen`    | weak-coupling velocity rescaling (time `tau_t`) |
/// | `ThermostatType::Csvr`         | exact stochastic velocity rescaling of Bussi, Donadio & Parrinello (2007), canonical, time `tau_t` |
/// | `ThermostatType::NoseHoover`   | Nosé-Hoover chain of three thermostats (Martyna *et al.* 1992/1996), period `tau_t` |
/// | `BarostatType::Berendsen`      | isotropic weak-coupling box and coordinate scaling (time `tau_p`) |
/// | `BarostatType::Parrinello`     | Parrinello-Rahman box dynamics with persistent box velocity (isotropic pressure, orthorhombic cell) |
/// | `BarostatType::MonteCarlo`     | Monte-Carlo volume moves, see [`McVolumeMoves`] |
///
/// All three thermostats use `3N` degrees of freedom, the same count as
/// [`MdSimulationPlain::temperature`].  The pressure that the barostats see is
/// `P = (2·KE + W) / (3V)` with the pair virial `W` of the force routine.
///
/// Stochastic components draw from a seeded stream, see
/// [`MdSimulationPlain::set_seed`]; a run is reproducible.
#[derive(Debug, Clone)]
pub struct MdSimulationPlain {
    /// All atoms in the simulation.
    pub atoms: Vec<Atom>,
    /// Orthorhombic box lengths [Lx, Ly, Lz] (Å).
    pub box_lengths: [f64; 3],
    /// Number of completed integration steps.
    pub step_count: u64,
    /// Simulated time (ps).
    pub time: f64,
    /// Optional thermostat for NVT / NPT.
    pub thermostat: Option<ThermostatType>,
    /// Optional barostat for NPT.
    pub barostat: Option<BarostatType>,
    /// Thermostat coupling time (ps): the relaxation time of the Berendsen and
    /// CSVR thermostats and the characteristic time `τ` of the Nosé-Hoover chain
    /// (`Q_1 = N_f k_B T τ²`).
    pub tau_t: f64,
    /// Target temperature for thermostat (K).
    pub target_temp: f64,
    /// Barostat coupling time (ps): the relaxation time of the Berendsen
    /// barostat and the period of the box oscillation of the Parrinello-Rahman
    /// barostat.
    pub tau_p: f64,
    /// Target pressure of the barostat (bar).
    pub target_pressure: f64,
    /// Isothermal compressibility of the system (bar⁻¹) used by the Berendsen
    /// and Parrinello-Rahman barostats (liquid water: `4.5e-5`).
    pub compressibility: f64,
    /// State of the Nosé-Hoover chain (persistent between steps).
    pub nose_hoover: NoseHooverChainState,
    /// Box-edge velocities of the Parrinello-Rahman barostat (Å ps⁻¹).
    pub box_velocity: [f64; 3],
    /// Settings and statistics of the Monte-Carlo barostat.
    pub monte_carlo: McVolumeMoves,
    /// Scalar pair virial (kJ mol⁻¹) of the latest force evaluation, if the
    /// force routine supplied it.
    pub virial: Option<f64>,
    /// Potential energy (kJ mol⁻¹) of the latest force evaluation, if the force
    /// routine supplied it.
    pub potential_energy: Option<f64>,
    /// Do the `force` fields of the atoms belong to the current positions?
    forces_current: bool,
    /// Random-number stream of the stochastic thermostat and barostat.
    rng: SimRng,
}

impl MdSimulationPlain {
    /// Create a new simulation in the NVE ensemble.
    pub fn new(atoms: Vec<Atom>, box_lengths: [f64; 3]) -> Self {
        Self {
            atoms,
            box_lengths,
            step_count: 0,
            time: 0.0,
            thermostat: None,
            barostat: None,
            tau_t: 0.1,
            target_temp: 300.0,
            tau_p: 1.0,
            target_pressure: 1.0,
            compressibility: 4.5e-5,
            nose_hoover: NoseHooverChainState::new(),
            box_velocity: [0.0; 3],
            monte_carlo: McVolumeMoves::default(),
            virial: None,
            potential_energy: None,
            forces_current: false,
            rng: SimRng::new(DEFAULT_SEED),
        }
    }

    /// Restart the random-number stream of the stochastic thermostat (CSVR) and
    /// barostat (Monte-Carlo) from `seed`.
    pub fn set_seed(&mut self, seed: u64) {
        self.rng.reseed(seed);
    }

    /// Builder form of [`MdSimulationPlain::set_seed`].
    pub fn with_seed(mut self, seed: u64) -> Self {
        self.set_seed(seed);
        self
    }

    /// Number of atoms.
    #[inline]
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Number of degrees of freedom `3N` used for the temperature and by the
    /// thermostats.
    #[inline]
    pub fn degrees_of_freedom(&self) -> usize {
        3 * self.atoms.len()
    }

    /// Volume of the orthorhombic box (Å³).
    #[inline]
    pub fn volume(&self) -> f64 {
        self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2]
    }

    // -----------------------------------------------------------------------
    // PBC wrapping
    // -----------------------------------------------------------------------

    /// Wrap all atom positions into the primary box [0, L) using PBC.
    pub fn apply_pbc(&mut self) {
        for atom in &mut self.atoms {
            for a in 0..3 {
                let l = self.box_lengths[a];
                atom.position[a] = atom.position[a].rem_euclid(l);
            }
        }
    }

    // -----------------------------------------------------------------------
    // Energy, temperature and pressure
    // -----------------------------------------------------------------------

    /// Kinetic energy (kJ mol⁻¹): KE = ½ Σ mᵢ |vᵢ|² · 0.01 (masses in amu,
    /// velocities in Å ps⁻¹).
    pub fn kinetic_energy(&self) -> f64 {
        self.atoms.iter().map(|a| a.kinetic_energy()).sum()
    }

    /// Instantaneous temperature (K): T = 2 KE / (3 N k_B).
    pub fn temperature(&self) -> f64 {
        let n = self.n_atoms() as f64;
        if n < 1.0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (3.0 * n * KB_REDUCED)
    }

    /// Instantaneous pressure (bar) using the virial theorem.
    ///
    /// ```text
    /// P = (2·KE + W) / (3 V)
    /// ```
    ///
    /// converted to bar (1 kJ mol⁻¹ Å⁻³ = 16 605.4 bar).  `virial` is the scalar
    /// **pair** virial `W = Σ_{i<j} r_ij · F_ij` (kJ mol⁻¹, minimum-image
    /// separations, see [`ForceEvaluation::virial`]); a non-interacting gas has
    /// `W = 0`.  Returns `0.0` for a degenerate box.
    pub fn pressure_virial(&self, virial: f64) -> f64 {
        let v = self.volume();
        if !(v.is_finite() && v > 1e-30) {
            return 0.0;
        }
        (2.0 * self.kinetic_energy() + virial) / (3.0 * v) * KJ_MOL_PER_ANG3_TO_BAR
    }

    /// Instantaneous pressure (bar) from the virial of the latest force
    /// evaluation, or `None` if the force routine did not supply one.
    pub fn pressure(&self) -> Option<f64> {
        self.virial.map(|w| self.pressure_virial(w))
    }

    /// Total energy KE + PE (kJ mol⁻¹), or `None` if no potential energy is
    /// known.
    pub fn total_energy(&self) -> Option<f64> {
        self.potential_energy.map(|u| u + self.kinetic_energy())
    }

    /// Extended energy `KE + PE + H_bath` (kJ mol⁻¹), or `None` if no potential
    /// energy is known.
    ///
    /// For the Nosé-Hoover chain thermostat `H_bath` is the energy of the chain
    /// ([`NoseHooverChainState::bath_energy`]) and, at constant volume, the sum is
    /// the quantity that the equations of motion conserve (a barostat does work
    /// on the system and adds terms that are not tracked here); for every other
    /// setting `H_bath` is zero and the result equals
    /// [`MdSimulationPlain::total_energy`].
    pub fn extended_energy(&self) -> Option<f64> {
        let bath = match self.thermostat {
            Some(ThermostatType::NoseHoover) => self.nose_hoover.bath_energy(),
            _ => 0.0,
        };
        self.total_energy().map(|e| e + bath)
    }

    // -----------------------------------------------------------------------
    // Velocity-Verlet integration
    // -----------------------------------------------------------------------

    /// Perform one velocity-Verlet step in the constant external forces stored in
    /// the atoms (`atom.force`, zero by default: free particles).
    ///
    /// The forces are *frozen*: they are not recomputed, so the system is a set
    /// of non-interacting particles in a uniform external field, its pair virial
    /// is exactly zero and the potential energy is reported as zero.  The
    /// configured thermostat and barostat act as in
    /// [`MdSimulationPlain::step_with`]; the pressure coupling sees the ideal-gas
    /// pressure `N k_B T / V`.  For an interacting system call
    /// [`MdSimulationPlain::step_with`], which recomputes the forces at the new
    /// positions between the two half kicks.
    ///
    /// # Errors
    ///
    /// Returns an error for a non-positive or non-finite `dt`, a non-positive
    /// mass or box edge, or invalid thermostat or barostat parameters.
    pub fn step(&mut self, dt: f64) -> Result<()> {
        self.validate(dt)?;
        self.virial = Some(0.0);
        self.potential_energy = Some(0.0);
        self.forces_current = true;
        let mut frozen = |atoms: &[Atom], _box_lengths: &[f64; 3]| {
            ForceEvaluation::new(atoms.iter().map(|atom| atom.force).collect())
                .with_virial(0.0)
                .with_potential_energy(0.0)
        };
        let outcome = self.advance(dt, &mut frozen);
        self.forces_current = false;
        outcome
    }

    /// Perform one velocity-Verlet step with forces recomputed by `evaluate`.
    ///
    /// `evaluate(atoms, box_lengths)` returns the [`ForceEvaluation`] of the given
    /// configuration.  One step is
    ///
    /// 1. thermostat half step (Nosé-Hoover chain / CSVR),
    /// 2. half kick `v += F/(2m)·dt`, drift `x += v·dt`,
    /// 3. pressure coupling (Berendsen / Parrinello-Rahman: dilate the box and
    ///    the coordinates by the factors that follow from the pressure of the
    ///    configuration at the start of the step), periodic wrap,
    /// 4. **force evaluation at the new positions**, half kick with the new
    ///    forces,
    /// 5. thermostat half step (Berendsen: one full-step rescaling),
    /// 6. Monte-Carlo volume move every `monte_carlo.interval` steps.
    ///
    /// If the forces of the current positions are not known yet (a new
    /// simulation, or after [`MdSimulationPlain::invalidate_forces`]) they are
    /// evaluated first.
    ///
    /// # Errors
    ///
    /// Returns an error for invalid parameters (see [`MdSimulationPlain::step`]),
    /// for a force routine that returns a wrong number of forces, and when the
    /// configured barostat needs a quantity the routine did not supply (the
    /// virial for Berendsen / Parrinello-Rahman, the potential energy for
    /// Monte-Carlo).  After an error the step is incomplete and the state must
    /// not be used for further stepping.
    pub fn step_with<F>(&mut self, dt: f64, mut evaluate: F) -> Result<()>
    where
        F: FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation,
    {
        self.validate(dt)?;
        if !self.forces_current {
            self.refresh_forces(&mut evaluate)?;
        }
        self.advance(dt, &mut evaluate)
    }

    /// Evaluate the forces at the current positions with `evaluate` and store
    /// them (together with the virial and the potential energy).
    ///
    /// Call it before reading [`MdSimulationPlain::pressure`] on a fresh
    /// simulation, or after changing positions or the box by hand.
    ///
    /// # Errors
    ///
    /// Returns an error if the routine returns a wrong number of forces.
    pub fn refresh_forces<F>(&mut self, mut evaluate: F) -> Result<()>
    where
        F: FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation,
    {
        let evaluation = evaluate(&self.atoms, &self.box_lengths);
        self.store_evaluation(evaluation)
    }

    /// Mark the stored forces as stale, so that the next
    /// [`MdSimulationPlain::step_with`] re-evaluates them first.
    pub fn invalidate_forces(&mut self) {
        self.forces_current = false;
    }

    /// Run `n_steps` frozen-force steps of size `dt`, returning (KE, T) after
    /// each step.
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`MdSimulationPlain::step`].
    pub fn run(&mut self, n_steps: u64, dt: f64) -> Result<Vec<(f64, f64)>> {
        let mut records = Vec::with_capacity(usize::try_from(n_steps).unwrap_or(0));
        for _ in 0..n_steps {
            self.step(dt)?;
            records.push((self.kinetic_energy(), self.temperature()));
        }
        Ok(records)
    }

    /// Run `n_steps` steps of size `dt` with forces from `evaluate`, returning
    /// (KE, T) after each step.
    ///
    /// # Errors
    ///
    /// Propagates the errors of [`MdSimulationPlain::step_with`].
    pub fn run_with<F>(&mut self, n_steps: u64, dt: f64, mut evaluate: F) -> Result<Vec<(f64, f64)>>
    where
        F: FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation,
    {
        let mut records = Vec::with_capacity(usize::try_from(n_steps).unwrap_or(0));
        for _ in 0..n_steps {
            self.step_with(dt, &mut evaluate)?;
            records.push((self.kinetic_energy(), self.temperature()));
        }
        Ok(records)
    }

    // -----------------------------------------------------------------------
    // Step machinery
    // -----------------------------------------------------------------------

    /// Check `dt` and the configuration of the integrator, thermostat and
    /// barostat.
    fn validate(&self, dt: f64) -> Result<()> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(invalid(format!(
                "time step must be positive and finite, got {dt}"
            )));
        }
        if !self.box_lengths.iter().all(|&l| l.is_finite() && l > 0.0) {
            return Err(invalid(format!(
                "box lengths must be positive and finite, got {:?}",
                self.box_lengths
            )));
        }
        if let Some(index) = self
            .atoms
            .iter()
            .position(|atom| !(atom.mass.is_finite() && atom.mass > 0.0))
        {
            return Err(invalid(format!(
                "atom {index} has the non-positive or non-finite mass {}",
                self.atoms[index].mass
            )));
        }
        if let Some(thermostat) = self.thermostat {
            if !(self.tau_t.is_finite() && self.tau_t > 0.0) {
                return Err(invalid(format!(
                    "thermostat coupling time tau_t must be positive and finite, got {}",
                    self.tau_t
                )));
            }
            let needs_positive_temperature = thermostat == ThermostatType::NoseHoover;
            let temperature_ok = self.target_temp.is_finite()
                && (self.target_temp > 0.0
                    || (self.target_temp == 0.0 && !needs_positive_temperature));
            if !temperature_ok {
                return Err(invalid(format!(
                    "target temperature {} K is invalid for the {thermostat:?} thermostat",
                    self.target_temp
                )));
            }
        }
        if let Some(barostat) = self.barostat {
            if !(self.tau_p.is_finite() && self.tau_p > 0.0) {
                return Err(invalid(format!(
                    "barostat coupling time tau_p must be positive and finite, got {}",
                    self.tau_p
                )));
            }
            if !(self.compressibility.is_finite() && self.compressibility > 0.0) {
                return Err(invalid(format!(
                    "compressibility must be positive and finite, got {}",
                    self.compressibility
                )));
            }
            if !self.target_pressure.is_finite() {
                return Err(invalid(format!(
                    "target pressure must be finite, got {}",
                    self.target_pressure
                )));
            }
            if barostat == BarostatType::MonteCarlo {
                let change = self.monte_carlo.max_relative_change;
                if !(change.is_finite() && change > 0.0 && change < 1.0) {
                    return Err(invalid(format!(
                        "Monte-Carlo barostat needs 0 < max_relative_change < 1, got {change}"
                    )));
                }
                if !(self.target_temp.is_finite() && self.target_temp > 0.0) {
                    return Err(invalid(format!(
                        "Monte-Carlo barostat needs a positive target temperature, got {} K",
                        self.target_temp
                    )));
                }
            }
        }
        Ok(())
    }

    /// Store an evaluation: forces into the atoms, virial and potential energy
    /// into the simulation.
    fn store_evaluation(&mut self, evaluation: ForceEvaluation) -> Result<()> {
        if evaluation.forces.len() != self.atoms.len() {
            return Err(invalid(format!(
                "the force routine returned {} forces for {} atoms",
                evaluation.forces.len(),
                self.atoms.len()
            )));
        }
        if self.barostat == Some(BarostatType::MonteCarlo) && evaluation.potential_energy.is_none()
        {
            return Err(invalid(
                "the Monte-Carlo barostat needs the potential energy: \
                 attach it with ForceEvaluation::with_potential_energy"
                    .to_string(),
            ));
        }
        for (atom, force) in self.atoms.iter_mut().zip(evaluation.forces) {
            atom.force = force;
        }
        self.virial = evaluation.virial;
        self.potential_energy = evaluation.potential_energy;
        self.forces_current = true;
        Ok(())
    }

    /// Pressure (bar) of the current configuration for the barostat: the stored
    /// pair virial and the current kinetic energy and volume.
    fn coupling_pressure(&self) -> Result<f64> {
        match self.virial {
            Some(virial) => Ok(self.pressure_virial(virial)),
            None => Err(invalid(
                "pressure coupling needs the pair virial: \
                 attach it with ForceEvaluation::with_virial"
                    .to_string(),
            )),
        }
    }

    /// One velocity-Verlet step for which the forces of the current positions
    /// are known.
    fn advance<F>(&mut self, dt: f64, evaluate: &mut F) -> Result<()>
    where
        F: FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation,
    {
        let half_dt = 0.5 * dt;

        // Pressure of the configuration at the start of the step: virial of the
        // latest evaluation, kinetic energy of the full-step velocities.
        let coupling_pressure = match self.barostat {
            Some(BarostatType::Berendsen) | Some(BarostatType::Parrinello) => {
                Some(self.coupling_pressure()?)
            }
            _ => None,
        };

        self.apply_chain_thermostat(half_dt);
        self.kick(half_dt);
        self.drift(dt);

        if let Some(pressure) = coupling_pressure {
            match self.barostat {
                Some(BarostatType::Berendsen) => self.berendsen_dilate(pressure, dt),
                Some(BarostatType::Parrinello) => self.parrinello_dilate(pressure, dt),
                _ => {}
            }
        }
        self.apply_pbc();

        // New forces at the new positions (and box), then the second half kick.
        let evaluation = evaluate(&self.atoms, &self.box_lengths);
        self.store_evaluation(evaluation)?;
        self.kick(half_dt);

        self.apply_chain_thermostat(half_dt);
        if self.thermostat == Some(ThermostatType::Berendsen) {
            self.apply_berendsen_rescale(self.target_temp, dt);
        }

        if self.barostat == Some(BarostatType::MonteCarlo) {
            let interval = self.monte_carlo.interval;
            if interval > 0 && (self.step_count + 1).is_multiple_of(interval) {
                self.monte_carlo_volume_move(evaluate)?;
            }
        }

        self.step_count += 1;
        self.time += dt;
        Ok(())
    }

    /// Kick every velocity with the stored forces for the duration `dt`
    /// (`Δv = F/m·dt`, with the unit factor `100` of
    /// [`ACCEL_PER_FORCE_OVER_MASS`](super::ACCEL_PER_FORCE_OVER_MASS)).
    fn kick(&mut self, dt: f64) {
        for atom in &mut self.atoms {
            let dv = velocity_increment(&atom.force, atom.mass, dt);
            for (v, d) in atom.velocity.iter_mut().zip(dv.iter()) {
                *v += d;
            }
        }
    }

    /// Advance every position with the current velocity for the duration `dt`.
    fn drift(&mut self, dt: f64) {
        for atom in &mut self.atoms {
            for (x, v) in atom.position.iter_mut().zip(atom.velocity.iter()) {
                *x += dt * v;
            }
        }
    }

    /// Multiply every velocity by `scale`.
    fn scale_velocities(&mut self, scale: f64) {
        for atom in &mut self.atoms {
            for v in atom.velocity.iter_mut() {
                *v *= scale;
            }
        }
    }

    /// Half step of the time-reversible thermostats (Nosé-Hoover chain, CSVR)
    /// for the duration `half_dt`; a no-op for the other settings and for a
    /// system at rest.
    fn apply_chain_thermostat(&mut self, half_dt: f64) {
        let kinetic = self.kinetic_energy();
        if !(kinetic.is_finite() && kinetic > 0.0) {
            return;
        }
        let n_dof = self.degrees_of_freedom() as f64;
        let kt = KB_REDUCED * self.target_temp;
        let scale = match self.thermostat {
            Some(ThermostatType::NoseHoover) => {
                self.nose_hoover
                    .propagate(half_dt, 2.0 * kinetic, n_dof, kt, self.tau_t)
            }
            Some(ThermostatType::Csvr) => {
                let target_kinetic = 0.5 * n_dof * kt;
                let rescaled = csvr_kinetic_energy(
                    kinetic,
                    target_kinetic,
                    n_dof,
                    half_dt,
                    self.tau_t,
                    &mut self.rng,
                );
                (rescaled / kinetic).sqrt()
            }
            _ => return,
        };
        if scale.is_finite() {
            self.scale_velocities(scale);
        }
    }

    /// Apply Berendsen velocity rescaling toward `target_t` (K) with coupling `tau_t`.
    fn apply_berendsen_rescale(&mut self, target_t: f64, dt: f64) {
        let t_current = self.temperature();
        if t_current < 1e-10 {
            return;
        }
        let ratio = target_t / t_current;
        let scale = (1.0 + dt / self.tau_t * (ratio - 1.0)).max(0.0).sqrt();
        self.scale_velocities(scale);
    }

    /// Berendsen barostat: scale the box and every coordinate isotropically.
    fn berendsen_dilate(&mut self, pressure: f64, dt: f64) {
        let barostat = NptBarostat::new(self.target_pressure, self.tau_p, self.compressibility);
        let mu = barostat.scale_factor(pressure, dt);
        for atom in &mut self.atoms {
            for x in atom.position.iter_mut() {
                *x *= mu;
            }
        }
        for edge in self.box_lengths.iter_mut() {
            *edge *= mu;
        }
    }

    /// Parrinello-Rahman barostat: advance the box with its persistent
    /// velocity, dilate the coordinates and scale the peculiar velocities.
    fn parrinello_dilate(&mut self, pressure: f64, dt: f64) {
        let mut barostat = ParrinelloRahmanBarostat {
            target_pressure: self.target_pressure,
            tau_p: self.tau_p,
            compressibility: self.compressibility,
            box_velocity: self.box_velocity,
        };
        let scales = barostat.advance_edges(&mut self.box_lengths, pressure, dt);
        self.box_velocity = barostat.box_velocity;
        for atom in &mut self.atoms {
            for (x, scale) in atom.position.iter_mut().zip(scales.iter()) {
                *x *= scale;
            }
            for (v, scale) in atom.velocity.iter_mut().zip(scales.iter()) {
                *v /= scale;
            }
        }
    }

    /// One Monte-Carlo volume move (see [`McVolumeMoves`]).
    fn monte_carlo_volume_move<F>(&mut self, evaluate: &mut F) -> Result<()>
    where
        F: FnMut(&[Atom], &[f64; 3]) -> ForceEvaluation,
    {
        let potential_old = self.potential_energy.ok_or_else(|| {
            invalid("the Monte-Carlo barostat needs the potential energy".to_string())
        })?;
        let volume_old = self.volume();
        // Symmetric proposal in ln V: ln(V'/V) uniform in [-ln(1+δ), +ln(1+δ)].
        let max_log_ratio = self.monte_carlo.max_relative_change.ln_1p();
        let log_ratio = max_log_ratio * (2.0 * self.rng.uniform() - 1.0);
        let volume_new = volume_old * log_ratio.exp();
        let delta_volume = volume_new - volume_old;
        self.monte_carlo.attempted += 1;
        if !(volume_new.is_finite() && volume_new > 0.0) {
            return Ok(());
        }

        let scale = (log_ratio / 3.0).exp();
        let old_positions: Vec<[f64; 3]> = self.atoms.iter().map(|atom| atom.position).collect();
        let old_box = self.box_lengths;
        for atom in &mut self.atoms {
            for x in atom.position.iter_mut() {
                *x *= scale;
            }
        }
        for edge in self.box_lengths.iter_mut() {
            *edge *= scale;
        }

        let trial = evaluate(&self.atoms, &self.box_lengths);
        let trial_is_valid =
            trial.forces.len() == self.atoms.len() && trial.potential_energy.is_some();
        let potential_new = trial.potential_energy.unwrap_or(f64::NAN);
        let kt = KB_REDUCED * self.target_temp;
        let pressure_energy_density = self.target_pressure / KJ_MOL_PER_ANG3_TO_BAR;
        let exponent = -(potential_new - potential_old + pressure_energy_density * delta_volume)
            / kt
            + (self.n_atoms() as f64 + 1.0) * log_ratio;
        let accept = trial_is_valid && (exponent >= 0.0 || self.rng.uniform() < exponent.exp());

        if accept {
            self.monte_carlo.accepted += 1;
            return self.store_evaluation(trial);
        }

        for (atom, position) in self.atoms.iter_mut().zip(old_positions) {
            atom.position = position;
        }
        self.box_lengths = old_box;
        if trial_is_valid {
            Ok(())
        } else {
            Err(invalid(
                "the force routine returned a wrong number of forces or no potential \
                 energy for the Monte-Carlo trial configuration"
                    .to_string(),
            ))
        }
    }

    /// Compute the virial scalar W = Σᵢ rᵢ · Fᵢ using the stored forces.
    ///
    /// This is the *atomic* virial of the stored coordinates, which equals the
    /// pair virial `Σ_{i<j} r_ij · F_ij` only for an isolated (non-periodic)
    /// system whose forces are pairwise and whose positions are not wrapped.  Under
    /// periodic boundary conditions it contains a spurious `L · F_ij` term for
    /// every pair that interacts across a cell face; supply the pair virial
    /// through [`ForceEvaluation::with_virial`] instead.
    pub fn virial_approx(&self) -> f64 {
        self.atoms
            .iter()
            .map(|a| {
                a.position[0] * a.force[0] + a.position[1] * a.force[1] + a.position[2] * a.force[2]
            })
            .sum()
    }

    /// Remove centre-of-mass velocity.
    pub fn remove_com_velocity(&mut self) {
        let n = self.n_atoms();
        if n == 0 {
            return;
        }
        let mut total_mass = 0.0_f64;
        let mut com_v = [0.0_f64; 3];
        for atom in &self.atoms {
            total_mass += atom.mass;
            for (cv, &v) in com_v.iter_mut().zip(atom.velocity.iter()) {
                *cv += atom.mass * v;
            }
        }
        if total_mass > 1e-20 {
            for cv in com_v.iter_mut() {
                *cv /= total_mass;
            }
            for atom in &mut self.atoms {
                for (v, &cv) in atom.velocity.iter_mut().zip(com_v.iter()) {
                    *v -= cv;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// ObservableRecord — observable history for analysis
// ---------------------------------------------------------------------------

/// Time-series record of MD observables collected during a run.
#[derive(Debug, Clone, Default)]
pub struct ObservableRecord {
    /// Simulation time (ps) at each recorded frame.
    pub time: Vec<f64>,
    /// Kinetic energy (kJ mol⁻¹) per frame.
    pub kinetic_energy: Vec<f64>,
    /// Potential energy (kJ mol⁻¹) per frame.
    pub potential_energy: Vec<f64>,
    /// Total energy (kJ mol⁻¹) per frame.
    pub total_energy: Vec<f64>,
    /// Instantaneous temperature (K) per frame.
    pub temperature: Vec<f64>,
    /// Instantaneous pressure (bar) per frame.
    pub pressure: Vec<f64>,
    /// Mean squared displacement (Å²) per frame.
    pub msd: Vec<f64>,
}

impl ObservableRecord {
    /// Create a new empty record.
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a frame.
    pub fn push(&mut self, time: f64, ke: f64, pe: f64, temp: f64, pressure: f64, msd: f64) {
        self.time.push(time);
        self.kinetic_energy.push(ke);
        self.potential_energy.push(pe);
        self.total_energy.push(ke + pe);
        self.temperature.push(temp);
        self.pressure.push(pressure);
        self.msd.push(msd);
    }

    /// Number of recorded frames.
    pub fn len(&self) -> usize {
        self.time.len()
    }

    /// Is the record empty?
    pub fn is_empty(&self) -> bool {
        self.time.is_empty()
    }

    /// Mean temperature over all frames.
    pub fn mean_temperature(&self) -> f64 {
        if self.temperature.is_empty() {
            return 0.0;
        }
        self.temperature.iter().sum::<f64>() / self.temperature.len() as f64
    }

    /// Mean total energy over all frames.
    pub fn mean_total_energy(&self) -> f64 {
        if self.total_energy.is_empty() {
            return 0.0;
        }
        self.total_energy.iter().sum::<f64>() / self.total_energy.len() as f64
    }

    /// Variance of total energy over all frames.
    ///
    /// Var(E) = `E²` - `E`²
    pub fn energy_variance(&self) -> f64 {
        let n = self.total_energy.len();
        if n < 2 {
            return 0.0;
        }
        let mean = self.mean_total_energy();
        let sum_sq: f64 = self
            .total_energy
            .iter()
            .map(|e| (e - mean) * (e - mean))
            .sum();
        sum_sq / (n as f64 - 1.0)
    }

    /// Maximum absolute energy drift from frame 0.
    pub fn max_energy_drift(&self) -> f64 {
        if self.total_energy.is_empty() {
            return 0.0;
        }
        let e0 = self.total_energy[0];
        self.total_energy
            .iter()
            .map(|&e| (e - e0).abs())
            .fold(0.0_f64, f64::max)
    }

    /// Relative energy drift: max_drift / |E₀|.
    pub fn relative_energy_drift(&self) -> f64 {
        if self.total_energy.is_empty() {
            return 0.0;
        }
        let e0 = self.total_energy[0].abs();
        if e0 < 1e-30 {
            return 0.0;
        }
        self.max_energy_drift() / e0
    }
}

// ---------------------------------------------------------------------------
// MSD — Mean Squared Displacement
// ---------------------------------------------------------------------------

/// Compute the mean squared displacement (Å²) between two sets of positions.
///
/// MSD = (1/N) Σᵢ |rᵢ(t) - rᵢ(0)|²
///
/// Note: this uses unwrapped positions (no PBC folding) for correct MSD.
pub fn compute_msd(positions_now: &[[f64; 3]], positions_ref: &[[f64; 3]]) -> f64 {
    let n = positions_now.len().min(positions_ref.len());
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for i in 0..n {
        for a in 0..3 {
            let d = positions_now[i][a] - positions_ref[i][a];
            sum += d * d;
        }
    }
    sum / n as f64
}

/// Compute per-atom squared displacements (Å²).
///
/// Returns a vector of |rᵢ(t) - rᵢ(0)|² for each atom.
pub fn per_atom_displacement_sq(
    positions_now: &[[f64; 3]],
    positions_ref: &[[f64; 3]],
) -> Vec<f64> {
    let n = positions_now.len().min(positions_ref.len());
    let mut result = Vec::with_capacity(n);
    for i in 0..n {
        let mut d2 = 0.0;
        for a in 0..3 {
            let d = positions_now[i][a] - positions_ref[i][a];
            d2 += d * d;
        }
        result.push(d2);
    }
    result
}
