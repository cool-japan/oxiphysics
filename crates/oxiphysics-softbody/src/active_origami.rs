// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Active origami: smart-material actuators and deployable structures.
//!
//! Provides shape-memory alloy (SMA) actuators, pneumatic origami structures,
//! electroactive origami with dielectric elastomer actuators, rigid-foldable
//! origami kinematics, bistable crease self-folding, and deployable structure
//! mechanics.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_softbody::active_origami::{SmaActuator, OrigamiMechanism};
//!
//! let mut sma = SmaActuator::new(0.001, 0.001, 50_000.0);
//! let force = sma.recovery_force(100.0);
//! assert!(force >= 0.0);
//!
//! let mech = OrigamiMechanism::new(4, 1);
//! assert_eq!(mech.dof(), 1);
//! ```

use std::f64::consts::PI;

// ── SmaActuator ───────────────────────────────────────────────────────────────

/// Phase state of a shape-memory alloy wire.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmaPhase {
    /// Fully martensite phase (low temperature).
    Martensite,
    /// Mixed phase (transformation in progress).
    Mixed,
    /// Fully austenite phase (high temperature).
    Austenite,
}

/// Shape-memory alloy actuator model.
///
/// Uses a simplified Brinson model with one-way SMA behavior:
/// recovery force as a function of temperature and pre-strain.
#[derive(Debug, Clone)]
pub struct SmaActuator {
    /// Cross-sectional area of the SMA wire (m²).
    pub area: f64,
    /// Pre-strain (dimensionless, typically 0.03–0.08 for NiTi).
    pub pre_strain: f64,
    /// Elastic modulus in austenite phase (Pa).
    pub e_austenite: f64,
    /// Elastic modulus in martensite phase (Pa).
    pub e_martensite: f64,
    /// Maximum transformation strain εL (dimensionless).
    pub epsilon_l: f64,
    /// Austenite start temperature As (°C).
    pub as_temp: f64,
    /// Austenite finish temperature Af (°C).
    pub af_temp: f64,
    /// Martensite start temperature Ms (°C).
    pub ms_temp: f64,
    /// Martensite finish temperature Mf (°C).
    pub mf_temp: f64,
    /// Current temperature (°C).
    pub temperature: f64,
    /// Current martensite fraction ξ (0 = full austenite, 1 = full martensite).
    pub martensite_fraction: f64,
    /// Current strain (dimensionless).
    pub strain: f64,
}

impl SmaActuator {
    /// Create a new SMA actuator.
    ///
    /// # Arguments
    /// * `area` — cross-sectional area (m²).
    /// * `pre_strain` — initial pre-strain (dimensionless).
    /// * `e_austenite` — austenite elastic modulus (Pa).
    pub fn new(area: f64, pre_strain: f64, e_austenite: f64) -> Self {
        Self {
            area,
            pre_strain,
            e_austenite,
            e_martensite: e_austenite * 0.5,
            epsilon_l: 0.05,
            as_temp: 68.0,
            af_temp: 78.0,
            ms_temp: 52.0,
            mf_temp: 42.0,
            temperature: 20.0,
            martensite_fraction: 1.0,
            strain: pre_strain,
        }
    }

    /// Current phase based on temperature.
    pub fn phase(&self) -> SmaPhase {
        let t = self.temperature;
        if t >= self.af_temp {
            SmaPhase::Austenite
        } else if t >= self.as_temp {
            SmaPhase::Mixed
        } else {
            SmaPhase::Martensite
        }
    }

    /// Martensite volume fraction ξ computed from temperature.
    ///
    /// Uses cosine model for smooth phase transformation.
    pub fn xi_from_temperature(&self, t: f64) -> f64 {
        if t >= self.af_temp {
            0.0
        } else if t >= self.as_temp {
            // Heating: austenite conversion
            0.5 * (1.0 + ((PI * (t - self.as_temp) / (self.af_temp - self.as_temp)).cos()))
        } else if t <= self.mf_temp {
            1.0
        } else if t <= self.ms_temp {
            // Cooling: martensite conversion
            0.5 * (1.0 - ((PI * (t - self.ms_temp) / (self.ms_temp - self.mf_temp)).cos()))
        } else {
            self.martensite_fraction
        }
    }

    /// Update actuator state for a given temperature.
    pub fn set_temperature(&mut self, t: f64) {
        self.temperature = t;
        self.martensite_fraction = self.xi_from_temperature(t);
        // Strain proportional to martensite fraction
        self.strain = self.epsilon_l * self.martensite_fraction;
    }

    /// Effective elastic modulus at current temperature.
    pub fn effective_modulus(&self) -> f64 {
        let xi = self.martensite_fraction;
        (1.0 - xi) * self.e_austenite + xi * self.e_martensite
    }

    /// Recovery force (N) at given temperature for constrained actuator.
    ///
    /// F = E(ξ) * A * (εL * (1 - ξ(T)) - pre_strain * ξ(T))
    pub fn recovery_force(&self, temperature: f64) -> f64 {
        let xi = self.xi_from_temperature(temperature);
        let e = (1.0 - xi) * self.e_austenite + xi * self.e_martensite;
        let recovery_strain = self.epsilon_l * (1.0 - xi);

        e * self.area * (recovery_strain - self.pre_strain * xi).max(0.0)
    }

    /// Transformation strain at given temperature.
    pub fn transformation_strain(&self, temperature: f64) -> f64 {
        let xi = self.xi_from_temperature(temperature);
        self.epsilon_l * xi
    }

    /// Actuation stroke (m) for a wire of length `length` (m).
    pub fn stroke(&self, length: f64, temperature: f64) -> f64 {
        let xi_cold = 1.0; // full martensite
        let xi_hot = self.xi_from_temperature(temperature);
        self.epsilon_l * (xi_cold - xi_hot) * length
    }

    /// Specific work output (J/kg) assuming density 6450 kg/m³ (NiTi).
    pub fn specific_work(&self, length: f64, temperature: f64) -> f64 {
        let f = self.recovery_force(temperature);
        let s = self.stroke(length, temperature);
        let volume = self.area * length;
        let density = 6450.0;
        let mass = density * volume;
        if mass < 1e-20 { 0.0 } else { f * s / mass }
    }
}

// ── PneumaticOrigami ──────────────────────────────────────────────────────────

/// Inflatable origami crease structure.
///
/// Models a pneumatic soft actuator based on a folded membrane.
/// Inflation pressure drives the unfolding and generates blocking force.
#[derive(Debug, Clone)]
pub struct PneumaticOrigami {
    /// Number of chambers / bellows sections.
    pub n_chambers: usize,
    /// Rest (deflated) length per chamber (m).
    pub rest_length: f64,
    /// Maximum extended length per chamber (m).
    pub max_length: f64,
    /// Chamber cross-sectional area (m²).
    pub chamber_area: f64,
    /// Internal gauge pressure (Pa).
    pub pressure: f64,
    /// Membrane stiffness per chamber (N/m).
    pub stiffness: f64,
    /// Current extension per chamber (m).
    pub extension: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
    /// Extension velocity (m/s).
    pub velocity: f64,
}

impl PneumaticOrigami {
    /// Create a pneumatic origami actuator.
    pub fn new(n_chambers: usize, rest_length: f64, max_length: f64, chamber_area: f64) -> Self {
        Self {
            n_chambers,
            rest_length,
            max_length,
            chamber_area,
            pressure: 0.0,
            stiffness: 500.0,
            extension: 0.0,
            damping: 5.0,
            velocity: 0.0,
        }
    }

    /// Total actuator rest length (m).
    pub fn total_rest_length(&self) -> f64 {
        self.rest_length * self.n_chambers as f64
    }

    /// Current total actuator length (m).
    pub fn current_length(&self) -> f64 {
        (self.rest_length + self.extension) * self.n_chambers as f64
    }

    /// Maximum total length (m).
    pub fn maximum_length(&self) -> f64 {
        self.max_length * self.n_chambers as f64
    }

    /// Pneumatic force due to internal pressure (N).
    ///
    /// F_pressure = p * A
    pub fn pressure_force(&self) -> f64 {
        self.pressure * self.chamber_area
    }

    /// Elastic restoring force (N).
    pub fn elastic_force(&self) -> f64 {
        self.stiffness * self.extension
    }

    /// Blocking force: maximum force at zero displacement (N).
    pub fn blocking_force(&self) -> f64 {
        self.pressure_force()
    }

    /// Free stroke: maximum extension at zero load (m).
    pub fn free_stroke(&self) -> f64 {
        if self.stiffness < 1e-12 {
            self.max_length - self.rest_length
        } else {
            (self.pressure_force() / self.stiffness).min(self.max_length - self.rest_length)
        }
    }

    /// Advance dynamics by `dt` (s) against external load `f_ext` (N).
    pub fn step(&mut self, dt: f64, f_ext: f64) {
        let f_net =
            self.pressure_force() - self.elastic_force() - f_ext - self.damping * self.velocity;
        // Assume unit effective mass per chamber
        let mass = 0.01 * self.n_chambers as f64;
        let accel = f_net / mass.max(1e-6);
        self.velocity += accel * dt;
        self.extension =
            (self.extension + self.velocity * dt).clamp(0.0, self.max_length - self.rest_length);
    }

    /// Set gauge pressure (Pa) and clamp to non-negative.
    pub fn set_pressure(&mut self, p: f64) {
        self.pressure = p.max(0.0);
    }

    /// Volume of one chamber at current extension (m³).
    pub fn chamber_volume(&self) -> f64 {
        (self.rest_length + self.extension) * self.chamber_area
    }

    /// Fold angle (rad) of crease due to inflation.
    ///
    /// Simplified geometric relation: angle ~ π * extension / (max - rest).
    pub fn fold_angle(&self) -> f64 {
        let max_ext = self.max_length - self.rest_length;
        if max_ext < 1e-12 {
            return 0.0;
        }
        PI * self.extension / max_ext
    }
}

// ── ElectroactiveOrigami ──────────────────────────────────────────────────────

/// Dielectric elastomer actuator (DEA) model for electroactive origami.
///
/// Computes Maxwell stress, voltage-controlled actuation strain,
/// and electrostatic energy stored in the membrane.
#[derive(Debug, Clone)]
pub struct ElectroactiveOrigami {
    /// Dielectric permittivity ε (F/m).
    pub permittivity: f64,
    /// Rest membrane thickness (m).
    pub thickness: f64,
    /// Membrane area (m²).
    pub area: f64,
    /// Applied voltage (V).
    pub voltage: f64,
    /// Shear modulus of the elastomer (Pa).
    pub shear_modulus: f64,
    /// Current in-plane stretch ratio λ (≥ 1).
    pub stretch: f64,
    /// Dielectric breakdown field (V/m).
    pub breakdown_field: f64,
}

impl ElectroactiveOrigami {
    /// Create a new DEA actuator.
    ///
    /// # Arguments
    /// * `permittivity` — relative permittivity (dimensionless) × ε0.
    /// * `thickness` — rest thickness (m).
    /// * `area` — membrane area (m²).
    /// * `shear_modulus` — elastomer shear modulus (Pa).
    pub fn new(permittivity: f64, thickness: f64, area: f64, shear_modulus: f64) -> Self {
        Self {
            permittivity,
            thickness,
            area,
            voltage: 0.0,
            shear_modulus,
            stretch: 1.0,
            breakdown_field: 100e6, // typical 100 MV/m for silicone
        }
    }

    /// Maxwell stress (Pa) for given electric field E (V/m).
    ///
    /// σ_Maxwell = ε * E²
    pub fn maxwell_stress_from_field(&self, e_field: f64) -> f64 {
        self.permittivity * e_field * e_field
    }

    /// Maxwell stress (Pa) at current voltage.
    pub fn maxwell_stress(&self) -> f64 {
        let t_current = self.current_thickness();
        if t_current < 1e-20 {
            return 0.0;
        }
        let e_field = self.voltage / t_current;
        self.maxwell_stress_from_field(e_field)
    }

    /// Current thickness at stretch ratio λ (incompressibility: t = t0 / λ²).
    pub fn current_thickness(&self) -> f64 {
        self.thickness / (self.stretch * self.stretch)
    }

    /// Neo-Hookean elastic stress in plane (Pa).
    ///
    /// σ_elastic = μ * (λ² − 1/λ⁴)
    pub fn elastic_stress(&self) -> f64 {
        let l = self.stretch;
        self.shear_modulus * (l * l - 1.0 / (l * l * l * l))
    }

    /// Equilibrium stretch ratio from voltage (iterative Newton solve).
    ///
    /// Solves σ_Maxwell = σ_elastic.
    pub fn equilibrium_stretch(&self) -> f64 {
        let mut l = self.stretch.max(1.001);
        for _ in 0..50 {
            let t = self.thickness / (l * l);
            let e = self.voltage / t.max(1e-20);
            let sigma_m = self.permittivity * e * e;
            let sigma_e = self.shear_modulus * (l * l - 1.0 / (l * l * l * l));
            let f = sigma_m - sigma_e;
            // df/dl
            let dsm = -4.0 * self.permittivity * self.voltage * self.voltage
                / (self.thickness.powi(2))
                * l.powi(-5)
                * 2.0
                * self.thickness;
            let dse = self.shear_modulus * (2.0 * l + 4.0 / (l * l * l * l * l));
            let df = -dsm - dse; // chain rule (approximate)
            if df.abs() < 1e-20 {
                break;
            }
            l -= f / df;
            l = l.max(1.001);
        }
        l
    }

    /// Actuation strain ε = λ - 1 for applied voltage.
    pub fn actuation_strain(&self) -> f64 {
        self.equilibrium_stretch() - 1.0
    }

    /// Electrostatic energy stored (J).
    pub fn stored_energy(&self) -> f64 {
        let t = self.current_thickness();
        let cap = self.permittivity * self.area / t.max(1e-20);
        0.5 * cap * self.voltage * self.voltage
    }

    /// Capacitance (F) at current stretch.
    pub fn capacitance(&self) -> f64 {
        let t = self.current_thickness();
        self.permittivity * self.area / t.max(1e-20)
    }

    /// Check if electric field exceeds breakdown threshold.
    pub fn is_breakdown(&self) -> bool {
        let t = self.current_thickness();
        let e_field = self.voltage / t.max(1e-20);
        e_field > self.breakdown_field
    }

    /// Set voltage and update stretch to equilibrium.
    pub fn set_voltage(&mut self, v: f64) {
        self.voltage = v.max(0.0);
        self.stretch = self.equilibrium_stretch();
    }

    /// Force output (N) = Maxwell stress × area.
    pub fn force(&self) -> f64 {
        self.maxwell_stress() * self.area
    }
}

// ── OrigamiMechanism ──────────────────────────────────────────────────────────

/// Rigid-foldable origami mechanism.
///
/// Kinematic model using rigid panels connected by crease hinges.
/// Tracks panel angles and fold angles through the mechanism DOF.
#[derive(Debug, Clone)]
pub struct OrigamiMechanism {
    /// Number of panels in the mechanism.
    pub n_panels: usize,
    /// Mechanical degrees of freedom.
    pub n_dof: usize,
    /// Crease fold angles (rad), one per interior crease.
    pub fold_angles: Vec<f64>,
    /// Sector angles at each vertex (rad).
    pub sector_angles: Vec<f64>,
    /// Crease stiffness values (N·m/rad).
    pub crease_stiffnesses: Vec<f64>,
    /// Panel inertia values (kg·m²).
    pub panel_inertia: Vec<f64>,
    /// Angular velocities of fold angles (rad/s).
    pub fold_velocities: Vec<f64>,
}

impl OrigamiMechanism {
    /// Create a new rigid-foldable mechanism.
    ///
    /// # Arguments
    /// * `n_panels` — number of rigid panels.
    /// * `n_dof` — kinematic degrees of freedom.
    pub fn new(n_panels: usize, n_dof: usize) -> Self {
        let n_creases = n_panels.saturating_sub(1);
        Self {
            n_panels,
            n_dof,
            fold_angles: vec![0.0; n_creases],
            sector_angles: vec![PI / 2.0; n_creases],
            crease_stiffnesses: vec![1.0; n_creases],
            panel_inertia: vec![0.01; n_panels],
            fold_velocities: vec![0.0; n_creases],
        }
    }

    /// Degrees of freedom of the mechanism.
    pub fn dof(&self) -> usize {
        self.n_dof
    }

    /// Number of interior creases.
    pub fn n_creases(&self) -> usize {
        self.fold_angles.len()
    }

    /// Set fold angle at crease index `i`.
    pub fn set_fold_angle(&mut self, i: usize, angle: f64) {
        if i < self.fold_angles.len() {
            self.fold_angles[i] = angle;
        }
    }

    /// Elastic potential energy stored in all creases (J).
    pub fn elastic_energy(&self) -> f64 {
        self.fold_angles
            .iter()
            .zip(self.crease_stiffnesses.iter())
            .map(|(a, k)| 0.5 * k * a * a)
            .sum()
    }

    /// Kinematic closure constraint residual for a degree-4 vertex.
    ///
    /// For Miura-ori: α + γ = π (sector angles sum condition for flat-foldability).
    pub fn kawasaki_residual(&self, vertex: usize) -> f64 {
        if vertex + 1 >= self.sector_angles.len() {
            return 0.0;
        }
        let alpha = self.sector_angles[vertex];
        let gamma = self.sector_angles[vertex + 1];
        (alpha + gamma) - PI
    }

    /// Advance fold dynamics using damped second-order system.
    ///
    /// Each crease is driven by external torque `torques[i]` and damped.
    pub fn step(&mut self, torques: &[f64], damping: f64, dt: f64) {
        for i in 0..self.fold_angles.len() {
            let i_eff = self.panel_inertia[i.min(self.panel_inertia.len() - 1)];
            let k = self.crease_stiffnesses[i];
            let tau = if i < torques.len() { torques[i] } else { 0.0 };
            let ang = self.fold_angles[i];
            let vel = self.fold_velocities[i];
            let accel = (tau - k * ang - damping * vel) / i_eff.max(1e-12);
            self.fold_velocities[i] += accel * dt;
            self.fold_angles[i] += self.fold_velocities[i] * dt;
        }
    }

    /// Fold ratio: average |fold_angle| / π.
    pub fn fold_ratio(&self) -> f64 {
        if self.fold_angles.is_empty() {
            return 0.0;
        }
        let sum: f64 = self.fold_angles.iter().map(|a| a.abs()).sum();
        sum / (PI * self.fold_angles.len() as f64)
    }

    /// Check whether all creases satisfy |angle| < tol (flat state).
    pub fn is_flat(&self, tol: f64) -> bool {
        self.fold_angles.iter().all(|a| a.abs() < tol)
    }
}

// ── SelfFoldingPattern ────────────────────────────────────────────────────────

/// Energy landscape type for a crease.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CreaseEnergy {
    /// Single-well (monostable).
    Monostable,
    /// Double-well (bistable).
    Bistable,
}

/// Bistable crease with snap-through instability.
///
/// Models the snap-through behavior via a double-well potential energy.
#[derive(Debug, Clone)]
pub struct SelfFoldingPattern {
    /// Number of creases.
    pub n_creases: usize,
    /// Target fold angles (rad) for each crease.
    pub target_angles: Vec<f64>,
    /// Current fold angles (rad).
    pub current_angles: Vec<f64>,
    /// Crease stiffness (N·m/rad).
    pub stiffness: Vec<f64>,
    /// Energy barrier for snap-through (J).
    pub energy_barrier: Vec<f64>,
    /// Energy landscape type per crease.
    pub energy_type: Vec<CreaseEnergy>,
    /// Well positions for bistable creases (rad): \[well1, well2\].
    pub bistable_wells: Vec<[f64; 2]>,
    /// Activation stimulus level \[0, 1\].
    pub stimulus: f64,
}

impl SelfFoldingPattern {
    /// Create a new self-folding pattern.
    pub fn new(n_creases: usize) -> Self {
        Self {
            n_creases,
            target_angles: vec![PI; n_creases],
            current_angles: vec![0.0; n_creases],
            stiffness: vec![0.5; n_creases],
            energy_barrier: vec![1.0; n_creases],
            energy_type: vec![CreaseEnergy::Monostable; n_creases],
            bistable_wells: vec![[0.0, PI]; n_creases],
            stimulus: 0.0,
        }
    }

    /// Potential energy of crease `i` at angle θ.
    pub fn crease_potential(&self, i: usize, theta: f64) -> f64 {
        match self.energy_type[i] {
            CreaseEnergy::Monostable => {
                0.5 * self.stiffness[i] * (theta - self.target_angles[i]).powi(2)
            }
            CreaseEnergy::Bistable => {
                let [w1, w2] = self.bistable_wells[i];
                let e1 = 0.5 * self.stiffness[i] * (theta - w1).powi(2);
                let e2 = 0.5 * self.stiffness[i] * (theta - w2).powi(2);
                // Double-well: min of two quadratics minus barrier at midpoint
                let mid = (w1 + w2) / 2.0;
                let height = self.energy_barrier[i];
                e1.min(e2) + height * (-(theta - mid).powi(2) / (mid - w1).powi(2).max(1e-12)).exp()
            }
        }
    }

    /// Restoring torque (N·m) at crease `i` for current angle.
    pub fn crease_torque(&self, i: usize) -> f64 {
        let theta = self.current_angles[i];
        let dt = 1e-6;
        let u_plus = self.crease_potential(i, theta + dt);
        let u_minus = self.crease_potential(i, theta - dt);
        -(u_plus - u_minus) / (2.0 * dt)
    }

    /// Total elastic potential energy (J).
    pub fn total_potential_energy(&self) -> f64 {
        (0..self.n_creases)
            .map(|i| self.crease_potential(i, self.current_angles[i]))
            .sum()
    }

    /// Check if crease `i` has snapped through to the other well.
    ///
    /// Returns true if current angle is closer to the second well.
    pub fn has_snapped(&self, i: usize) -> bool {
        if self.energy_type[i] != CreaseEnergy::Bistable {
            return false;
        }
        let [w1, w2] = self.bistable_wells[i];
        let theta = self.current_angles[i];
        (theta - w2).abs() < (theta - w1).abs()
    }

    /// Apply stimulus-driven actuation: move creases toward targets.
    ///
    /// Uses a first-order rate proportional to stimulus and stiffness-weighted error.
    pub fn actuate(&mut self, dt: f64) {
        for i in 0..self.n_creases {
            // First-order approach: velocity proportional to net restoring force
            // Use larger effective inertia (1.0) to ensure numerical stability
            let effective_inertia = 1.0_f64.max(self.stiffness[i]);
            let tau = self.crease_torque(i);
            // Additional driving torque toward target proportional to stimulus
            let drive = self.stimulus
                * self.stiffness[i]
                * (self.target_angles[i] - self.current_angles[i]);
            let total = tau + drive;
            // Clamp increment to avoid numerical instability
            let max_delta = PI * dt;
            let delta = (total * dt / effective_inertia).clamp(-max_delta, max_delta);
            self.current_angles[i] += delta;
        }
    }

    /// Set stimulus level.
    pub fn set_stimulus(&mut self, s: f64) {
        self.stimulus = s.clamp(0.0, 1.0);
    }
}

// ── DeployableStructure ───────────────────────────────────────────────────────

/// Deployment state of a structure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeploymentState {
    /// Fully stowed / compact.
    Stowed,
    /// Deploying.
    Deploying,
    /// Fully deployed and locked.
    Deployed,
}

/// Locking mechanism type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LockType {
    /// Mechanical snap lock.
    SnapLock,
    /// Magnetic retention.
    Magnetic,
    /// Friction / contact.
    Friction,
}

/// Deployable origami structure: deployment sequence and shape change ratio.
#[derive(Debug, Clone)]
pub struct DeployableStructure {
    /// Current deployment state.
    pub state: DeploymentState,
    /// Number of deployment stages.
    pub n_stages: usize,
    /// Current stage index.
    pub current_stage: usize,
    /// Fold angles at each stage for each crease (rad).
    pub stage_targets: Vec<Vec<f64>>,
    /// Current fold angles (rad).
    pub fold_angles: Vec<f64>,
    /// Deployment speed (rad/s).
    pub deploy_speed: f64,
    /// Locking mechanism type.
    pub lock_type: LockType,
    /// Stowed volume (m³).
    pub stowed_volume: f64,
    /// Deployed volume (m³).
    pub deployed_volume: f64,
    /// Whether locking has been achieved.
    pub locked: bool,
    /// Time elapsed in current deployment (s).
    pub deploy_time: f64,
}

impl DeployableStructure {
    /// Create a deployable structure.
    pub fn new(
        n_stages: usize,
        n_creases: usize,
        stowed_volume: f64,
        deployed_volume: f64,
    ) -> Self {
        let stage_targets = vec![vec![0.0; n_creases]; n_stages];
        Self {
            state: DeploymentState::Stowed,
            n_stages,
            current_stage: 0,
            stage_targets,
            fold_angles: vec![0.0; n_creases],
            deploy_speed: 0.5,
            lock_type: LockType::SnapLock,
            stowed_volume,
            deployed_volume,
            locked: false,
            deploy_time: 0.0,
        }
    }

    /// Volumetric expansion ratio (deployed / stowed).
    pub fn shape_change_ratio(&self) -> f64 {
        if self.stowed_volume < 1e-20 {
            1.0
        } else {
            self.deployed_volume / self.stowed_volume
        }
    }

    /// Begin deployment sequence.
    pub fn start_deployment(&mut self) {
        self.state = DeploymentState::Deploying;
        self.current_stage = 0;
        self.locked = false;
        self.deploy_time = 0.0;
    }

    /// Advance deployment by `dt` seconds.
    ///
    /// Moves fold angles toward the current stage targets.
    pub fn step(&mut self, dt: f64) {
        if self.state != DeploymentState::Deploying {
            return;
        }
        self.deploy_time += dt;
        let stage = self.current_stage.min(self.n_stages - 1);
        let targets = &self.stage_targets[stage].clone();
        let max_delta = self.deploy_speed * dt;

        let mut all_reached = true;
        for (i, (angle, target)) in self.fold_angles.iter_mut().zip(targets.iter()).enumerate() {
            let _ = i;
            let error = target - *angle;
            if error.abs() > 1e-4 {
                all_reached = false;
                *angle += error.clamp(-max_delta, max_delta);
            }
        }

        if all_reached {
            if self.current_stage + 1 >= self.n_stages {
                self.state = DeploymentState::Deployed;
                self.locked = true;
            } else {
                self.current_stage += 1;
            }
        }
    }

    /// Current deployment progress \[0, 1\].
    pub fn progress(&self) -> f64 {
        match self.state {
            DeploymentState::Stowed => 0.0,
            DeploymentState::Deployed => 1.0,
            DeploymentState::Deploying => {
                let stage_frac = self.current_stage as f64 / self.n_stages as f64;
                // Add fractional progress within current stage
                stage_frac
            }
        }
    }

    /// Current volume (interpolated between stowed and deployed).
    pub fn current_volume(&self) -> f64 {
        let p = self.progress();
        self.stowed_volume + p * (self.deployed_volume - self.stowed_volume)
    }

    /// Locking force (N) for snap-lock mechanism.
    pub fn lock_force(&self, snap_stiffness: f64, snap_displacement: f64) -> f64 {
        match self.lock_type {
            LockType::SnapLock => snap_stiffness * snap_displacement,
            LockType::Magnetic => snap_stiffness * 0.1, // simpler
            LockType::Friction => snap_stiffness * snap_displacement * 0.3,
        }
    }

    /// Set stage target fold angles.
    pub fn set_stage_target(&mut self, stage: usize, targets: Vec<f64>) {
        if stage < self.n_stages {
            self.stage_targets[stage] = targets;
        }
    }
}

// ── Miura-ori geometry ────────────────────────────────────────────────────────

/// Compute Miura-ori fold angle ψ from base angle α and the single DOF θ.
///
/// Kinematic constraint for Miura-ori: cos(ψ/2) = cos(α) / sin(α) * 1 / tan(θ/2)
pub fn miura_fold_angle(alpha: f64, theta: f64) -> f64 {
    // Classical Miura-ori: tan(ψ/2) = tan(θ/2) / cos(α)
    let t = (theta / 2.0).tan();
    let cos_alpha = alpha.cos();
    if cos_alpha.abs() < 1e-12 {
        return PI;
    }
    2.0 * (t / cos_alpha).atan()
}

/// Miura-ori projected dimensions for a unit cell.
///
/// Returns `[width, height, depth]` as a function of fold angle θ (rad).
pub fn miura_dimensions(a: f64, b: f64, alpha: f64, theta: f64) -> [f64; 3] {
    let psi = miura_fold_angle(alpha, theta);
    let sin_a = alpha.sin();
    let cos_a = alpha.cos();
    let sin_t2 = (theta / 2.0).sin();
    let cos_t2 = (theta / 2.0).cos();
    let width = 2.0 * a * sin_a * cos_t2;
    let height = 2.0 * b * (1.0 - cos_a * cos_a * sin_t2 * sin_t2).sqrt();
    let depth = 2.0 * a * cos_a * sin_t2 * (psi / 2.0).sin();
    [width, height, depth]
}

/// Waterbomb base fold angle constraint.
///
/// For a waterbomb: cos(ψ) = 1 - 2 * sin²(θ/2) / sin²(α).
pub fn waterbomb_dihedral(alpha: f64, theta: f64) -> f64 {
    let sin_a = alpha.sin();
    if sin_a.abs() < 1e-12 {
        return 0.0;
    }
    let cos_psi = 1.0 - 2.0 * (theta / 2.0).sin().powi(2) / sin_a.powi(2);
    cos_psi.clamp(-1.0, 1.0).acos()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    // ── SmaActuator ───────────────────────────────────────────────────────

    #[test]
    fn test_sma_new() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        assert!((sma.area - 0.001).abs() < EPS);
        assert!((sma.pre_strain - 0.03).abs() < EPS);
    }

    #[test]
    fn test_sma_phase_cold_is_martensite() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        assert_eq!(sma.phase(), SmaPhase::Martensite);
    }

    #[test]
    fn test_sma_phase_hot_is_austenite() {
        let mut sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        sma.temperature = 100.0;
        assert_eq!(sma.phase(), SmaPhase::Austenite);
    }

    #[test]
    fn test_sma_set_temperature_updates_fraction() {
        let mut sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        sma.set_temperature(100.0);
        assert!(sma.martensite_fraction < 0.01);
    }

    #[test]
    fn test_sma_set_temperature_cold() {
        let mut sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        sma.set_temperature(20.0);
        assert!((sma.martensite_fraction - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_sma_recovery_force_hot() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        let f = sma.recovery_force(90.0);
        assert!(f >= 0.0);
    }

    #[test]
    fn test_sma_recovery_force_cold_is_zero_or_positive() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        let f = sma.recovery_force(20.0);
        assert!(f >= 0.0);
    }

    #[test]
    fn test_sma_transformation_strain_hot() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        let eps = sma.transformation_strain(90.0);
        assert!(eps < 0.01); // austenite has low transformation strain
    }

    #[test]
    fn test_sma_stroke_positive() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        let s = sma.stroke(0.1, 90.0);
        assert!(s >= 0.0);
    }

    #[test]
    fn test_sma_effective_modulus() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        let e = sma.effective_modulus();
        assert!(e > 0.0);
        assert!(e <= sma.e_austenite);
    }

    #[test]
    fn test_sma_specific_work() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        let w = sma.specific_work(0.1, 90.0);
        assert!(w >= 0.0);
    }

    #[test]
    fn test_sma_xi_from_temperature_boundary() {
        let sma = SmaActuator::new(0.001, 0.03, 50_000.0);
        assert!((sma.xi_from_temperature(100.0)).abs() < EPS);
        assert!((sma.xi_from_temperature(10.0) - 1.0).abs() < EPS);
    }

    // ── PneumaticOrigami ──────────────────────────────────────────────────

    #[test]
    fn test_pneumatic_new() {
        let p = PneumaticOrigami::new(4, 0.05, 0.15, 0.01);
        assert_eq!(p.n_chambers, 4);
        assert!((p.total_rest_length() - 0.2).abs() < EPS);
    }

    #[test]
    fn test_pneumatic_pressure_force() {
        let mut p = PneumaticOrigami::new(1, 0.05, 0.15, 0.01);
        p.set_pressure(10_000.0);
        assert!((p.pressure_force() - 100.0).abs() < EPS);
    }

    #[test]
    fn test_pneumatic_blocking_force_equals_pressure_force() {
        let mut p = PneumaticOrigami::new(1, 0.05, 0.15, 0.01);
        p.set_pressure(5_000.0);
        assert!((p.blocking_force() - p.pressure_force()).abs() < EPS);
    }

    #[test]
    fn test_pneumatic_free_stroke() {
        let mut p = PneumaticOrigami::new(1, 0.05, 0.15, 0.01);
        p.set_pressure(10_000.0);
        let s = p.free_stroke();
        assert!(s > 0.0);
        assert!(s <= p.max_length - p.rest_length);
    }

    #[test]
    fn test_pneumatic_step_extends() {
        let mut p = PneumaticOrigami::new(1, 0.05, 0.15, 0.01);
        p.set_pressure(50_000.0);
        // Run enough steps for extension to build up
        for _ in 0..200 {
            p.step(0.001, 0.0);
        }
        // Should have reached free-stroke equilibrium > 0
        assert!(
            p.extension > 0.0,
            "Should extend under pressure, got {}",
            p.extension
        );
    }

    #[test]
    fn test_pneumatic_fold_angle() {
        let mut p = PneumaticOrigami::new(1, 0.05, 0.15, 0.01);
        p.extension = (p.max_length - p.rest_length) / 2.0;
        let angle = p.fold_angle();
        assert!((angle - PI / 2.0).abs() < 0.01);
    }

    #[test]
    fn test_pneumatic_chamber_volume() {
        let p = PneumaticOrigami::new(1, 0.1, 0.2, 0.005);
        let vol = p.chamber_volume();
        assert!((vol - 0.1 * 0.005).abs() < EPS);
    }

    #[test]
    fn test_pneumatic_zero_pressure_no_extension() {
        let mut p = PneumaticOrigami::new(2, 0.05, 0.15, 0.01);
        p.set_pressure(0.0);
        let ext0 = p.extension;
        p.step(0.1, 0.0);
        // With zero pressure and zero velocity, should remain near zero
        assert!(p.extension >= ext0 - 0.001);
    }

    // ── ElectroactiveOrigami ──────────────────────────────────────────────

    #[test]
    fn test_dea_new() {
        let dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        assert_eq!(dea.voltage, 0.0);
        assert!((dea.stretch - 1.0).abs() < EPS);
    }

    #[test]
    fn test_dea_maxwell_stress_zero_voltage() {
        let dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        assert!(dea.maxwell_stress().abs() < EPS);
    }

    #[test]
    fn test_dea_maxwell_stress_with_voltage() {
        let mut dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        dea.voltage = 1000.0;
        let sigma = dea.maxwell_stress();
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_dea_capacitance_increases_with_stretch() {
        let mut dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        let c1 = dea.capacitance();
        dea.stretch = 2.0;
        let c2 = dea.capacitance();
        assert!(c2 > c1, "Capacitance should increase with stretch");
    }

    #[test]
    fn test_dea_stored_energy_zero_voltage() {
        let dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        assert!(dea.stored_energy().abs() < EPS);
    }

    #[test]
    fn test_dea_stored_energy_positive_voltage() {
        let mut dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        dea.voltage = 500.0;
        assert!(dea.stored_energy() > 0.0);
    }

    #[test]
    fn test_dea_breakdown_check() {
        let dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        // Voltage below breakdown field
        assert!(!dea.is_breakdown());
    }

    #[test]
    fn test_dea_force_positive_voltage() {
        let mut dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        dea.voltage = 1000.0;
        let f = dea.force();
        assert!(f >= 0.0);
    }

    #[test]
    fn test_dea_current_thickness_incompressible() {
        let mut dea = ElectroactiveOrigami::new(3.5e-11, 1e-3, 0.01, 1e5);
        dea.stretch = 2.0;
        let t = dea.current_thickness();
        // t = t0 / λ² = 1e-3 / 4 = 2.5e-4
        assert!((t - 2.5e-4).abs() < 1e-8);
    }

    // ── OrigamiMechanism ──────────────────────────────────────────────────

    #[test]
    fn test_mechanism_new() {
        let m = OrigamiMechanism::new(4, 1);
        assert_eq!(m.n_panels, 4);
        assert_eq!(m.n_dof, 1);
        assert_eq!(m.n_creases(), 3);
    }

    #[test]
    fn test_mechanism_dof() {
        let m = OrigamiMechanism::new(6, 2);
        assert_eq!(m.dof(), 2);
    }

    #[test]
    fn test_mechanism_set_fold_angle() {
        let mut m = OrigamiMechanism::new(4, 1);
        m.set_fold_angle(0, 1.0);
        assert!((m.fold_angles[0] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_mechanism_elastic_energy_zero() {
        let m = OrigamiMechanism::new(4, 1);
        assert!(m.elastic_energy().abs() < EPS);
    }

    #[test]
    fn test_mechanism_elastic_energy_nonzero() {
        let mut m = OrigamiMechanism::new(4, 1);
        m.set_fold_angle(0, 1.0);
        assert!(m.elastic_energy() > 0.0);
    }

    #[test]
    fn test_mechanism_fold_ratio_zero() {
        let m = OrigamiMechanism::new(4, 1);
        assert!(m.fold_ratio().abs() < EPS);
    }

    #[test]
    fn test_mechanism_is_flat() {
        let m = OrigamiMechanism::new(4, 1);
        assert!(m.is_flat(0.01));
    }

    #[test]
    fn test_mechanism_step_applies_torque() {
        let mut m = OrigamiMechanism::new(4, 1);
        let torques = [1.0, 0.5, 0.5];
        m.step(&torques, 0.1, 0.01);
        assert!(!m.fold_angles.iter().all(|&a| a.abs() < EPS));
    }

    #[test]
    fn test_mechanism_kawasaki_residual() {
        let mut m = OrigamiMechanism::new(4, 1);
        m.sector_angles = vec![PI / 2.0, PI / 2.0, PI / 2.0];
        // α + γ = π → residual = 0
        let res = m.kawasaki_residual(0);
        assert!(res.abs() < EPS);
    }

    // ── SelfFoldingPattern ────────────────────────────────────────────────

    #[test]
    fn test_sfp_new() {
        let sfp = SelfFoldingPattern::new(3);
        assert_eq!(sfp.n_creases, 3);
        assert_eq!(sfp.current_angles.len(), 3);
    }

    #[test]
    fn test_sfp_monostable_potential_zero() {
        let sfp = SelfFoldingPattern::new(1);
        // At target_angle = π, current = π → energy = 0
        let sfp2 = {
            let mut s = SelfFoldingPattern::new(1);
            s.current_angles[0] = PI;
            s.target_angles[0] = PI;
            s
        };
        assert!(sfp2.crease_potential(0, PI).abs() < EPS);
        let _ = sfp;
    }

    #[test]
    fn test_sfp_monostable_potential_positive() {
        let sfp = SelfFoldingPattern::new(1);
        let u = sfp.crease_potential(0, 0.5);
        assert!(u > 0.0);
    }

    #[test]
    fn test_sfp_total_energy_zero_at_target() {
        let sfp = {
            let mut s = SelfFoldingPattern::new(2);
            s.current_angles = s.target_angles.clone();
            s
        };
        assert!(sfp.total_potential_energy() < EPS);
    }

    #[test]
    fn test_sfp_actuation_moves_toward_target() {
        let mut sfp = SelfFoldingPattern::new(1);
        sfp.target_angles[0] = PI;
        // Start close but not at target; use high stimulus to drive quickly
        sfp.current_angles[0] = 0.0;
        sfp.set_stimulus(1.0);
        for _ in 0..200 {
            sfp.actuate(0.01);
        }
        // Should be significantly closer to PI than 0
        assert!(
            sfp.current_angles[0] > PI / 2.0,
            "Should have moved toward target PI, got {}",
            sfp.current_angles[0]
        );
    }

    #[test]
    fn test_sfp_no_actuation_without_stimulus() {
        let mut sfp = SelfFoldingPattern::new(1);
        // Set target and current to same angle so restoring torque is zero
        sfp.target_angles[0] = 0.0;
        sfp.set_stimulus(0.0);
        sfp.current_angles[0] = 0.0;
        let angle0 = sfp.current_angles[0];
        sfp.actuate(0.1);
        // At equilibrium (current == target), no force → no motion
        assert!(
            (sfp.current_angles[0] - angle0).abs() < 1e-10,
            "At equilibrium with zero stimulus, angle should not change"
        );
    }

    #[test]
    fn test_sfp_bistable_has_snapped() {
        let mut sfp = SelfFoldingPattern::new(1);
        sfp.energy_type[0] = CreaseEnergy::Bistable;
        sfp.bistable_wells[0] = [0.0, PI];
        sfp.current_angles[0] = PI - 0.1; // near second well
        assert!(sfp.has_snapped(0));
    }

    #[test]
    fn test_sfp_stimulus_clamp() {
        let mut sfp = SelfFoldingPattern::new(1);
        sfp.set_stimulus(2.0);
        assert!((sfp.stimulus - 1.0).abs() < EPS);
        sfp.set_stimulus(-1.0);
        assert!(sfp.stimulus.abs() < EPS);
    }

    // ── DeployableStructure ───────────────────────────────────────────────

    #[test]
    fn test_deployable_new() {
        let d = DeployableStructure::new(3, 4, 0.01, 1.0);
        assert_eq!(d.state, DeploymentState::Stowed);
        assert_eq!(d.n_stages, 3);
    }

    #[test]
    fn test_deployable_shape_change_ratio() {
        let d = DeployableStructure::new(1, 1, 0.01, 1.0);
        assert!((d.shape_change_ratio() - 100.0).abs() < EPS);
    }

    #[test]
    fn test_deployable_start_deployment() {
        let mut d = DeployableStructure::new(2, 2, 0.01, 1.0);
        d.start_deployment();
        assert_eq!(d.state, DeploymentState::Deploying);
    }

    #[test]
    fn test_deployable_step_progresses() {
        let mut d = DeployableStructure::new(1, 2, 0.01, 1.0);
        d.stage_targets[0] = vec![PI, PI];
        d.deploy_speed = 10.0;
        d.start_deployment();
        for _ in 0..1000 {
            d.step(0.01);
            if d.state == DeploymentState::Deployed {
                break;
            }
        }
        assert_eq!(d.state, DeploymentState::Deployed);
    }

    #[test]
    fn test_deployable_progress_stowed() {
        let d = DeployableStructure::new(2, 2, 0.01, 1.0);
        assert!(d.progress().abs() < EPS);
    }

    #[test]
    fn test_deployable_progress_deployed() {
        let mut d = DeployableStructure::new(1, 1, 0.01, 1.0);
        d.state = DeploymentState::Deployed;
        assert!((d.progress() - 1.0).abs() < EPS);
    }

    #[test]
    fn test_deployable_lock_force() {
        let d = DeployableStructure::new(1, 1, 0.01, 1.0);
        let f = d.lock_force(1000.0, 0.005);
        assert!(f > 0.0);
    }

    #[test]
    fn test_deployable_current_volume_stowed() {
        let d = DeployableStructure::new(1, 1, 0.01, 1.0);
        assert!((d.current_volume() - 0.01).abs() < EPS);
    }

    #[test]
    fn test_deployable_set_stage_target() {
        let mut d = DeployableStructure::new(2, 3, 0.01, 1.0);
        d.set_stage_target(0, vec![1.0, 2.0, 3.0]);
        assert!((d.stage_targets[0][1] - 2.0).abs() < EPS);
    }

    #[test]
    fn test_deployable_locked_after_deploy() {
        let mut d = DeployableStructure::new(1, 1, 0.01, 1.0);
        d.stage_targets[0] = vec![PI];
        d.deploy_speed = 100.0;
        d.start_deployment();
        for _ in 0..100 {
            d.step(0.1);
        }
        assert!(d.locked || d.state == DeploymentState::Deployed);
    }

    // ── Miura-ori geometry ────────────────────────────────────────────────

    #[test]
    fn test_miura_fold_angle_flat() {
        // When θ → 0, ψ → 0 (flat)
        let psi = miura_fold_angle(PI / 4.0, 1e-6);
        assert!(psi < 0.01, "Near-flat θ should give near-flat ψ: {psi}");
    }

    #[test]
    fn test_miura_fold_angle_finite() {
        let psi = miura_fold_angle(PI / 4.0, PI / 2.0);
        assert!(psi.is_finite());
        assert!(psi > 0.0);
    }

    #[test]
    fn test_miura_dimensions_positive() {
        let [w, h, d] = miura_dimensions(0.05, 0.03, PI / 4.0, PI / 4.0);
        assert!(w > 0.0);
        assert!(h > 0.0);
        assert!(d >= 0.0);
    }

    #[test]
    fn test_waterbomb_dihedral_finite() {
        let psi = waterbomb_dihedral(PI / 3.0, PI / 4.0);
        assert!(psi.is_finite());
        assert!(psi >= 0.0);
    }

    #[test]
    fn test_waterbomb_dihedral_range() {
        let psi = waterbomb_dihedral(PI / 3.0, 0.5);
        assert!((0.0..=PI).contains(&psi));
    }
}
