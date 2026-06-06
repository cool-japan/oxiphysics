// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tendon and ligament soft-body simulation.
//!
//! Provides biomechanical models for tendons and ligaments including:
//!
//! - **Quasi-linear viscoelastic (QLV)** model
//! - **Fung-type** exponential strain energy
//! - **Fiber-reinforced continuum** model
//! - **Crimping / uncrimping** mechanics (toe region)
//! - **Stress-strain curve** with toe, linear, and failure regions
//! - **Strain rate dependence**
//! - **Creep and stress relaxation**
//! - **Fatigue damage** accumulation
//! - **Tendon wrapping** around bones (obstacle set)
//! - **Muscle-tendon unit** (Hill model coupling)
//! - **Enthesis** (bone-tendon insertion) gradient properties
//! - **Pre-tension / initial strain**
//! - **Healing model** (collagen remodeling)

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Constants
// ─────────────────────────────────────────────────────────────────────────────

/// Small tolerance for numerical comparisons.
const EPS: f64 = 1.0e-15;

/// Default tendon elastic modulus (Pa) — roughly 1.5 GPa for mature tendon.
const DEFAULT_TENDON_MODULUS: f64 = 1.5e9;

/// Default tendon cross-sectional area (m^2).
const DEFAULT_CROSS_SECTION: f64 = 50.0e-6;

/// Default tendon rest length (m).
const DEFAULT_REST_LENGTH: f64 = 0.20;

/// Default Fung material constant C (Pa).
const DEFAULT_FUNG_C: f64 = 100.0;

/// Default Fung material constant c1 (dimensionless).
const DEFAULT_FUNG_C1: f64 = 50.0;

/// Default failure strain (dimensionless).
const DEFAULT_FAILURE_STRAIN: f64 = 0.12;

/// Default toe region limit strain (dimensionless).
const DEFAULT_TOE_LIMIT: f64 = 0.03;

/// Default crimp angle (radians).
const DEFAULT_CRIMP_ANGLE: f64 = 0.15;

/// Default density of tendon tissue (kg / m^3).
const DEFAULT_DENSITY: f64 = 1100.0;

/// Default Poisson's ratio for tendon (nearly incompressible).
const DEFAULT_POISSON: f64 = 0.49;

// ─────────────────────────────────────────────────────────────────────────────
// TendonMaterialProperties
// ─────────────────────────────────────────────────────────────────────────────

/// Material properties for a tendon or ligament.
#[derive(Debug, Clone)]
pub struct TendonMaterialProperties {
    /// Elastic modulus of the linear region (Pa).
    pub elastic_modulus: f64,
    /// Cross-sectional area (m^2).
    pub cross_section_area: f64,
    /// Rest (slack) length (m).
    pub rest_length: f64,
    /// Strain at the end of the toe region (dimensionless).
    pub toe_limit: f64,
    /// Ultimate failure strain (dimensionless).
    pub failure_strain: f64,
    /// Ultimate tensile strength (Pa).
    pub ultimate_stress: f64,
    /// Density (kg / m^3).
    pub density: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Pre-tension strain (dimensionless, typically 0.01-0.04).
    pub pretension_strain: f64,
    /// Crimp angle at rest (radians).
    pub crimp_angle: f64,
}

impl Default for TendonMaterialProperties {
    fn default() -> Self {
        Self::new()
    }
}

impl TendonMaterialProperties {
    /// Create new tendon material properties with defaults.
    pub fn new() -> Self {
        Self {
            elastic_modulus: DEFAULT_TENDON_MODULUS,
            cross_section_area: DEFAULT_CROSS_SECTION,
            rest_length: DEFAULT_REST_LENGTH,
            toe_limit: DEFAULT_TOE_LIMIT,
            failure_strain: DEFAULT_FAILURE_STRAIN,
            ultimate_stress: 100.0e6,
            density: DEFAULT_DENSITY,
            poisson_ratio: DEFAULT_POISSON,
            pretension_strain: 0.02,
            crimp_angle: DEFAULT_CRIMP_ANGLE,
        }
    }

    /// Create properties for a human Achilles tendon.
    pub fn achilles() -> Self {
        Self {
            elastic_modulus: 1.2e9,
            cross_section_area: 60.0e-6,
            rest_length: 0.15,
            toe_limit: 0.02,
            failure_strain: 0.08,
            ultimate_stress: 80.0e6,
            density: 1100.0,
            poisson_ratio: 0.49,
            pretension_strain: 0.03,
            crimp_angle: 0.12,
        }
    }

    /// Create properties for a patellar tendon.
    pub fn patellar() -> Self {
        Self {
            elastic_modulus: 1.0e9,
            cross_section_area: 120.0e-6,
            rest_length: 0.05,
            toe_limit: 0.025,
            failure_strain: 0.10,
            ultimate_stress: 65.0e6,
            density: 1100.0,
            poisson_ratio: 0.49,
            pretension_strain: 0.01,
            crimp_angle: 0.10,
        }
    }

    /// Create properties for a ligament (e.g. ACL).
    pub fn ligament_acl() -> Self {
        Self {
            elastic_modulus: 0.3e9,
            cross_section_area: 40.0e-6,
            rest_length: 0.035,
            toe_limit: 0.04,
            failure_strain: 0.15,
            ultimate_stress: 38.0e6,
            density: 1100.0,
            poisson_ratio: 0.45,
            pretension_strain: 0.005,
            crimp_angle: 0.18,
        }
    }

    /// Stiffness k = E A / L.
    pub fn stiffness(&self) -> f64 {
        self.elastic_modulus * self.cross_section_area / self.rest_length
    }

    /// Mass of the tendon.
    pub fn mass(&self) -> f64 {
        self.density * self.cross_section_area * self.rest_length
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CrimpModel
// ─────────────────────────────────────────────────────────────────────────────

/// Models the crimped collagen fiber structure that produces the toe region.
///
/// At rest, collagen fibers have a wavy (crimped) structure. Under small
/// strains the fibers straighten (uncrimp) with low stiffness. Beyond the
/// toe limit the fibers bear load linearly.
#[derive(Debug, Clone)]
pub struct CrimpModel {
    /// Crimp angle at zero strain (radians).
    pub initial_angle: f64,
    /// Strain at which fibers are fully straightened (toe limit).
    pub toe_limit: f64,
    /// Number of crimp waves per unit length.
    pub crimp_frequency: f64,
}

impl CrimpModel {
    /// Create a new crimp model.
    pub fn new(initial_angle: f64, toe_limit: f64) -> Self {
        Self {
            initial_angle,
            toe_limit,
            crimp_frequency: 100.0,
        }
    }

    /// Current crimp angle at a given strain.
    ///
    /// The angle decreases linearly from `initial_angle` to 0 as strain
    /// goes from 0 to `toe_limit`.
    pub fn angle_at_strain(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            self.initial_angle
        } else if strain >= self.toe_limit {
            0.0
        } else {
            self.initial_angle * (1.0 - strain / self.toe_limit)
        }
    }

    /// Fraction of fibers that are fully recruited (straightened).
    ///
    /// Returns a value in \[0, 1\].
    pub fn recruitment_fraction(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            0.0
        } else if strain >= self.toe_limit {
            1.0
        } else {
            strain / self.toe_limit
        }
    }

    /// Effective stiffness multiplier due to crimping.
    ///
    /// In the toe region the stiffness ramps up from ~0 to the full modulus.
    /// Uses a quadratic ramp for a smooth transition.
    pub fn stiffness_multiplier(&self, strain: f64) -> f64 {
        let r = self.recruitment_fraction(strain);
        r * r // quadratic recruitment
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FungStrainEnergy
// ─────────────────────────────────────────────────────────────────────────────

/// Fung-type exponential strain energy density.
///
/// `W = C/2 (exp(c1 E^2) − 1)`
///
/// where E is the Green-Lagrange strain, C and c1 are material constants.
/// This model captures the nonlinear stiffening behaviour of soft tissues.
#[derive(Debug, Clone)]
pub struct FungStrainEnergy {
    /// Material constant C (Pa).
    pub c: f64,
    /// Material constant c1 (dimensionless).
    pub c1: f64,
}

impl FungStrainEnergy {
    /// Create a new Fung model.
    pub fn new(c: f64, c1: f64) -> Self {
        Self { c, c1 }
    }

    /// Create with default tendon parameters.
    pub fn default_tendon() -> Self {
        Self::new(DEFAULT_FUNG_C, DEFAULT_FUNG_C1)
    }

    /// Strain energy density (J / m^3).
    pub fn energy(&self, strain: f64) -> f64 {
        0.5 * self.c * ((self.c1 * strain * strain).exp() - 1.0)
    }

    /// Second Piola-Kirchhoff stress (Pa): `dW/dE = C c1 E exp(c1 E^2)`.
    pub fn stress(&self, strain: f64) -> f64 {
        self.c * self.c1 * strain * (self.c1 * strain * strain).exp()
    }

    /// Tangent modulus: `d²W/dE² = C c1 (1 + 2 c1 E^2) exp(c1 E^2)`.
    pub fn tangent_modulus(&self, strain: f64) -> f64 {
        self.c
            * self.c1
            * (1.0 + 2.0 * self.c1 * strain * strain)
            * (self.c1 * strain * strain).exp()
    }

    /// Complementary energy density (Legendre transform approximation).
    pub fn complementary_energy(&self, stress: f64) -> f64 {
        // Numerical approximation: W* ≈ σε − W
        // Use Newton to find ε from σ, then compute
        let eps = self.inverse_stress(stress);
        stress * eps - self.energy(eps)
    }

    /// Inverse: find strain from stress using Newton-Raphson.
    pub fn inverse_stress(&self, target_stress: f64) -> f64 {
        if target_stress.abs() < EPS {
            return 0.0;
        }
        let mut eps = target_stress / (self.c * self.c1).max(1.0);
        for _ in 0..50 {
            let f = self.stress(eps) - target_stress;
            let df = self.tangent_modulus(eps);
            if df.abs() < EPS {
                break;
            }
            eps -= f / df;
            eps = eps.max(0.0);
        }
        eps
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TendonStressStrainCurve
// ─────────────────────────────────────────────────────────────────────────────

/// Piecewise stress-strain model with three regions:
///
/// 1. **Toe region** (0 ≤ ε ≤ ε_toe): quadratic/exponential ramp
/// 2. **Linear region** (ε_toe ≤ ε ≤ ε_fail): linear elastic
/// 3. **Failure region** (ε > ε_fail): stress drops to zero
#[derive(Debug, Clone)]
pub struct TendonStressStrainCurve {
    /// Material properties.
    pub props: TendonMaterialProperties,
    /// Crimp model for toe region.
    pub crimp: CrimpModel,
    /// Optional damage state (0 = undamaged, 1 = fully failed).
    pub damage: f64,
}

impl TendonStressStrainCurve {
    /// Create from material properties.
    pub fn new(props: TendonMaterialProperties) -> Self {
        let crimp = CrimpModel::new(props.crimp_angle, props.toe_limit);
        Self {
            props,
            crimp,
            damage: 0.0,
        }
    }

    /// Compute stress at a given strain.
    pub fn stress(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            return 0.0; // no compression bearing
        }
        if self.damage >= 1.0 || strain > self.props.failure_strain {
            return 0.0; // failed
        }
        let sigma = if strain <= self.props.toe_limit {
            // Toe region: quadratic
            let mult = self.crimp.stiffness_multiplier(strain);
            self.props.elastic_modulus * mult * strain
        } else {
            // Linear region
            let toe_stress = self.props.elastic_modulus * self.props.toe_limit;
            toe_stress + self.props.elastic_modulus * (strain - self.props.toe_limit)
        };
        sigma * (1.0 - self.damage)
    }

    /// Compute tangent modulus at a given strain.
    pub fn tangent_modulus(&self, strain: f64) -> f64 {
        if strain <= 0.0 || strain > self.props.failure_strain || self.damage >= 1.0 {
            return 0.0;
        }
        if strain <= self.props.toe_limit {
            let r = self.crimp.recruitment_fraction(strain);
            // d(E r^2 ε)/dε where r = ε/ε_toe
            // = E (3ε^2/ε_toe^2) approx
            self.props.elastic_modulus * 3.0 * r * r * (1.0 - self.damage)
        } else {
            self.props.elastic_modulus * (1.0 - self.damage)
        }
    }

    /// Generate stress-strain data points.
    pub fn generate_curve(&self, n_points: usize) -> Vec<(f64, f64)> {
        let n = n_points.max(2);
        let max_strain = self.props.failure_strain * 1.1;
        let mut result = Vec::with_capacity(n);
        for i in 0..n {
            let eps = max_strain * (i as f64) / (n as f64 - 1.0);
            result.push((eps, self.stress(eps)));
        }
        result
    }

    /// Energy stored (integral of stress-strain up to current strain).
    pub fn stored_energy(&self, strain: f64) -> f64 {
        // Numerical integration via trapezoidal rule
        let n = 200;
        let eps_max = strain.clamp(0.0, self.props.failure_strain);
        let de = eps_max / n as f64;
        let mut energy = 0.0;
        for i in 0..n {
            let e0 = de * i as f64;
            let e1 = de * (i + 1) as f64;
            energy += 0.5 * (self.stress(e0) + self.stress(e1)) * de;
        }
        energy
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// QuasiLinearViscoelastic (QLV)
// ─────────────────────────────────────────────────────────────────────────────

/// Quasi-linear viscoelastic (QLV) model for tendons and ligaments.
///
/// The stress at time t is a hereditary integral:
///
/// `σ(t) = ∫₀ᵗ G(t − τ) dσ^e/dτ dτ`
///
/// where `G(t)` is the reduced relaxation function and `σ^e` is the
/// elastic (instantaneous) stress response.
///
/// The relaxation function uses a Prony series:
/// `G(t) = G∞ + Σ gᵢ exp(−t/τᵢ)`
#[derive(Debug, Clone)]
pub struct QLVModel {
    /// Elastic stress function (Fung type).
    pub elastic: FungStrainEnergy,
    /// Prony series coefficients (gᵢ).
    pub prony_g: Vec<f64>,
    /// Prony series time constants (τᵢ, in seconds).
    pub prony_tau: Vec<f64>,
    /// Long-term modulus ratio G∞.
    pub g_infinity: f64,
    /// Internal state: hereditary integral components.
    pub state: Vec<f64>,
    /// Previous elastic stress.
    pub prev_elastic_stress: f64,
}

impl QLVModel {
    /// Create a new QLV model with a 3-term Prony series.
    pub fn new(elastic: FungStrainEnergy) -> Self {
        // Default 3-term Prony series for tendon (from Fung, 1993)
        let prony_g = vec![0.30, 0.20, 0.10];
        let prony_tau = vec![1.0, 10.0, 100.0]; // seconds
        let g_infinity = 1.0 - prony_g.iter().sum::<f64>();
        let n = prony_g.len();
        Self {
            elastic,
            prony_g,
            prony_tau,
            g_infinity,
            state: vec![0.0; n],
            prev_elastic_stress: 0.0,
        }
    }

    /// Create with custom Prony series.
    pub fn with_prony(elastic: FungStrainEnergy, prony_g: Vec<f64>, prony_tau: Vec<f64>) -> Self {
        let g_infinity = 1.0 - prony_g.iter().sum::<f64>();
        let n = prony_g.len();
        Self {
            elastic,
            prony_g,
            prony_tau,
            g_infinity: g_infinity.max(0.0),
            state: vec![0.0; n],
            prev_elastic_stress: 0.0,
        }
    }

    /// Reduced relaxation function G(t).
    pub fn relaxation_function(&self, t: f64) -> f64 {
        let mut g = self.g_infinity;
        for (gi, tau_i) in self.prony_g.iter().zip(self.prony_tau.iter()) {
            g += gi * (-t / tau_i).exp();
        }
        g
    }

    /// Update the QLV state and return the current viscoelastic stress.
    ///
    /// Uses a recursive formula for the hereditary integral.
    pub fn update(&mut self, strain: f64, dt: f64) -> f64 {
        let sigma_e = self.elastic.stress(strain);
        let dsigma = sigma_e - self.prev_elastic_stress;

        // Update Prony series internal variables
        for i in 0..self.prony_g.len() {
            let tau = self.prony_tau[i];
            let decay = (-dt / tau).exp();
            self.state[i] = decay * self.state[i] + self.prony_g[i] * dsigma;
        }

        self.prev_elastic_stress = sigma_e;

        // Total stress
        let mut sigma = self.g_infinity * sigma_e;
        for s in &self.state {
            sigma += s;
        }
        sigma
    }

    /// Reset internal state.
    pub fn reset(&mut self) {
        for s in &mut self.state {
            *s = 0.0;
        }
        self.prev_elastic_stress = 0.0;
    }

    /// Stress relaxation test: apply step strain and return stress vs time.
    pub fn stress_relaxation_test(
        &mut self,
        strain: f64,
        total_time: f64,
        n_steps: usize,
    ) -> Vec<(f64, f64)> {
        self.reset();
        let dt = total_time / n_steps as f64;
        let mut result = Vec::with_capacity(n_steps);
        // Apply step strain at t=0
        let _sigma0 = self.update(strain, dt);
        result.push((0.0, self.prev_elastic_stress));
        for i in 1..n_steps {
            let sigma = self.update(strain, dt);
            result.push((dt * i as f64, sigma));
        }
        result
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CreepModel
// ─────────────────────────────────────────────────────────────────────────────

/// Creep response under constant load.
///
/// Strain increases logarithmically with time:
/// `ε(t) = ε₀ + A ln(1 + t/t₀)`
///
/// where ε₀ is the initial elastic strain, A is the creep coefficient,
/// and t₀ is the characteristic time.
#[derive(Debug, Clone)]
pub struct CreepModel {
    /// Initial elastic strain.
    pub initial_strain: f64,
    /// Creep coefficient A.
    pub creep_coeff: f64,
    /// Characteristic time t₀ (s).
    pub char_time: f64,
    /// Current accumulated creep strain.
    pub creep_strain: f64,
    /// Elapsed time (s).
    pub elapsed_time: f64,
}

impl CreepModel {
    /// Create a new creep model.
    pub fn new(initial_strain: f64, creep_coeff: f64, char_time: f64) -> Self {
        Self {
            initial_strain,
            creep_coeff,
            char_time,
            creep_strain: 0.0,
            elapsed_time: 0.0,
        }
    }

    /// Default creep parameters for tendon.
    pub fn default_tendon(initial_strain: f64) -> Self {
        Self::new(initial_strain, 0.01, 10.0)
    }

    /// Total strain at time t.
    pub fn total_strain(&self, t: f64) -> f64 {
        self.initial_strain + self.creep_coeff * (1.0 + t / self.char_time).ln()
    }

    /// Creep strain rate at time t.
    pub fn creep_rate(&self, t: f64) -> f64 {
        self.creep_coeff / (self.char_time + t)
    }

    /// Advance the creep model by a time step.
    pub fn step(&mut self, dt: f64) {
        self.elapsed_time += dt;
        self.creep_strain = self.creep_coeff * (1.0 + self.elapsed_time / self.char_time).ln();
    }

    /// Current total strain.
    pub fn current_total_strain(&self) -> f64 {
        self.initial_strain + self.creep_strain
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StressRelaxation
// ─────────────────────────────────────────────────────────────────────────────

/// Stress relaxation model for constant-strain holding.
///
/// Stress decays as `σ(t) = σ₀ G(t)` where G(t) is the normalized
/// relaxation function from a Prony series.
#[derive(Debug, Clone)]
pub struct StressRelaxation {
    /// Initial peak stress (Pa).
    pub initial_stress: f64,
    /// Prony series coefficients.
    pub prony_g: Vec<f64>,
    /// Prony series time constants (s).
    pub prony_tau: Vec<f64>,
    /// Equilibrium fraction.
    pub g_eq: f64,
}

impl StressRelaxation {
    /// Create a default stress relaxation model.
    pub fn new(initial_stress: f64) -> Self {
        let prony_g = vec![0.25, 0.15, 0.10];
        let prony_tau = vec![1.0, 10.0, 100.0];
        let g_eq = 1.0 - prony_g.iter().sum::<f64>();
        Self {
            initial_stress,
            prony_g,
            prony_tau,
            g_eq,
        }
    }

    /// Stress at time t.
    pub fn stress_at(&self, t: f64) -> f64 {
        let mut g = self.g_eq;
        for (gi, tau_i) in self.prony_g.iter().zip(self.prony_tau.iter()) {
            g += gi * (-t / tau_i).exp();
        }
        self.initial_stress * g
    }

    /// Half-life: time at which stress drops to 50% of initial.
    pub fn half_life(&self) -> f64 {
        // Newton search for σ(t) = 0.5 σ₀
        let target = 0.5;
        let mut t = 1.0;
        for _ in 0..100 {
            let g = self.stress_at(t) / self.initial_stress;
            let err = g - target;
            if err.abs() < 1.0e-8 {
                break;
            }
            // dg/dt
            let mut dg = 0.0;
            for (gi, tau_i) in self.prony_g.iter().zip(self.prony_tau.iter()) {
                dg += -gi / tau_i * (-t / tau_i).exp();
            }
            if dg.abs() < EPS {
                break;
            }
            t -= err / dg;
            t = t.max(0.001);
        }
        t
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// StrainRateDependence
// ─────────────────────────────────────────────────────────────────────────────

/// Strain rate dependence model for tendon.
///
/// The effective modulus scales with strain rate:
/// `E_eff = E₀ (1 + β ln(1 + ε̇/ε̇₀))`
///
/// where β is the rate sensitivity parameter and ε̇₀ is a reference rate.
#[derive(Debug, Clone)]
pub struct StrainRateDependence {
    /// Base elastic modulus (Pa).
    pub base_modulus: f64,
    /// Rate sensitivity parameter β (dimensionless).
    pub beta: f64,
    /// Reference strain rate (1/s).
    pub ref_rate: f64,
}

impl StrainRateDependence {
    /// Create a new strain rate dependence model.
    pub fn new(base_modulus: f64, beta: f64, ref_rate: f64) -> Self {
        Self {
            base_modulus,
            beta,
            ref_rate,
        }
    }

    /// Default parameters for tendon.
    pub fn default_tendon() -> Self {
        Self::new(DEFAULT_TENDON_MODULUS, 0.05, 0.01)
    }

    /// Effective modulus at a given strain rate.
    pub fn effective_modulus(&self, strain_rate: f64) -> f64 {
        let rate_ratio = strain_rate.abs() / self.ref_rate;
        self.base_modulus * (1.0 + self.beta * (1.0 + rate_ratio).ln())
    }

    /// Rate enhancement factor (dimensionless multiplier).
    pub fn enhancement_factor(&self, strain_rate: f64) -> f64 {
        let rate_ratio = strain_rate.abs() / self.ref_rate;
        1.0 + self.beta * (1.0 + rate_ratio).ln()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FatigueDamage
// ─────────────────────────────────────────────────────────────────────────────

/// Fatigue damage accumulation model using a Miner-type linear rule.
///
/// Damage D accumulates per cycle proportionally to the applied strain
/// amplitude relative to the fatigue limit. When D ≥ 1 the tendon ruptures.
#[derive(Debug, Clone)]
pub struct FatigueDamage {
    /// Current damage level D ∈ \[0, 1\].
    pub damage: f64,
    /// Number of cycles applied.
    pub cycles: u64,
    /// Fatigue limit strain amplitude (below this, no damage accumulates).
    pub fatigue_limit: f64,
    /// Wohler exponent (S-N curve slope).
    pub wohler_exponent: f64,
    /// Reference cycles to failure at reference strain.
    pub ref_cycles: f64,
    /// Reference strain amplitude.
    pub ref_strain: f64,
}

impl Default for FatigueDamage {
    fn default() -> Self {
        Self::new()
    }
}

impl FatigueDamage {
    /// Create a new fatigue damage model with defaults for tendon.
    pub fn new() -> Self {
        Self {
            damage: 0.0,
            cycles: 0,
            fatigue_limit: 0.02,
            wohler_exponent: 5.0,
            ref_cycles: 1_000_000.0,
            ref_strain: 0.04,
        }
    }

    /// Number of cycles to failure at a given strain amplitude.
    pub fn cycles_to_failure(&self, strain_amplitude: f64) -> f64 {
        if strain_amplitude <= self.fatigue_limit {
            return f64::INFINITY;
        }
        let ratio = self.ref_strain / strain_amplitude;
        self.ref_cycles * ratio.powf(self.wohler_exponent)
    }

    /// Damage increment per cycle at a given strain amplitude.
    pub fn damage_per_cycle(&self, strain_amplitude: f64) -> f64 {
        let nf = self.cycles_to_failure(strain_amplitude);
        if nf.is_infinite() { 0.0 } else { 1.0 / nf }
    }

    /// Apply n cycles at the given strain amplitude.
    pub fn apply_cycles(&mut self, n: u64, strain_amplitude: f64) {
        let dd = self.damage_per_cycle(strain_amplitude) * n as f64;
        self.damage = (self.damage + dd).min(1.0);
        self.cycles += n;
    }

    /// Check if the tendon has ruptured.
    pub fn is_ruptured(&self) -> bool {
        self.damage >= 1.0
    }

    /// Remaining life (fraction).
    pub fn remaining_life(&self) -> f64 {
        (1.0 - self.damage).max(0.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// FiberReinforcedModel
// ─────────────────────────────────────────────────────────────────────────────

/// Fiber-reinforced continuum model for tendon.
///
/// Decomposes the strain energy into isotropic (ground substance) and
/// anisotropic (collagen fiber) contributions:
///
/// `W = W_iso(I1, I2) + W_fib(I4)`
///
/// where I4 = a₀ · C · a₀ is the fiber stretch invariant along
/// the preferred direction a₀.
#[derive(Debug, Clone)]
pub struct FiberReinforcedModel {
    /// Ground substance shear modulus (Pa).
    pub mu: f64,
    /// Ground substance bulk modulus (Pa).
    pub kappa: f64,
    /// Fiber direction (unit vector as \[f64; 3\]).
    pub fiber_direction: [f64; 3],
    /// Fiber stiffness parameter k1 (Pa).
    pub k1: f64,
    /// Fiber nonlinearity parameter k2 (dimensionless).
    pub k2: f64,
    /// Fiber dispersion parameter κ ∈ \[0, 1/3\].
    pub fiber_dispersion: f64,
}

impl FiberReinforcedModel {
    /// Create a new fiber-reinforced model.
    pub fn new(mu: f64, kappa: f64, k1: f64, k2: f64) -> Self {
        Self {
            mu,
            kappa,
            fiber_direction: [1.0, 0.0, 0.0],
            k1,
            k2,
            fiber_dispersion: 0.0,
        }
    }

    /// Default tendon parameters.
    pub fn default_tendon() -> Self {
        Self::new(0.3e6, 100.0e6, 50.0e6, 10.0)
    }

    /// Set fiber direction.
    pub fn with_direction(mut self, dir: [f64; 3]) -> Self {
        let norm = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        if norm > EPS {
            self.fiber_direction = [dir[0] / norm, dir[1] / norm, dir[2] / norm];
        }
        self
    }

    /// Isotropic ground substance stress (Neo-Hookean).
    pub fn ground_stress(&self, strain: f64) -> f64 {
        let lambda = 1.0 + strain;
        self.mu * (lambda - 1.0 / (lambda * lambda))
    }

    /// Fiber contribution to stress along the fiber direction.
    pub fn fiber_stress(&self, fiber_stretch: f64) -> f64 {
        let i4 = fiber_stretch * fiber_stretch;
        if i4 <= 1.0 {
            return 0.0; // fibers only resist extension
        }
        let e_term = i4 - 1.0;
        2.0 * self.k1 * e_term * (self.k2 * e_term * e_term).exp()
    }

    /// Total uniaxial stress along fiber direction.
    pub fn total_stress(&self, strain: f64) -> f64 {
        let lambda = 1.0 + strain;
        self.ground_stress(strain) + self.fiber_stress(lambda)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TendonWrapping
// ─────────────────────────────────────────────────────────────────────────────

/// Models tendon wrapping around a cylindrical obstacle (bone).
///
/// The tendon follows a straight-arc-straight path around a cylinder
/// defined by its center, radius, and axis direction.
#[derive(Debug, Clone)]
pub struct TendonWrapping {
    /// Obstacle center position \[x, y, z\].
    pub obstacle_center: [f64; 3],
    /// Obstacle radius (m).
    pub obstacle_radius: f64,
    /// Obstacle axis direction (unit vector).
    pub obstacle_axis: [f64; 3],
    /// Insertion point A \[x, y, z\].
    pub point_a: [f64; 3],
    /// Insertion point B \[x, y, z\].
    pub point_b: [f64; 3],
}

impl TendonWrapping {
    /// Create a new tendon wrapping model.
    pub fn new(
        center: [f64; 3],
        radius: f64,
        axis: [f64; 3],
        point_a: [f64; 3],
        point_b: [f64; 3],
    ) -> Self {
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let axis_n = if norm > EPS {
            [axis[0] / norm, axis[1] / norm, axis[2] / norm]
        } else {
            [0.0, 0.0, 1.0]
        };
        Self {
            obstacle_center: center,
            obstacle_radius: radius,
            obstacle_axis: axis_n,
            point_a,
            point_b,
        }
    }

    /// Compute the straight-line distance between A and B.
    pub fn direct_distance(&self) -> f64 {
        let dx = self.point_b[0] - self.point_a[0];
        let dy = self.point_b[1] - self.point_a[1];
        let dz = self.point_b[2] - self.point_a[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }

    /// Compute the distance from a point to the obstacle center
    /// projected onto the plane perpendicular to the axis.
    fn _projected_distance(&self, point: [f64; 3]) -> f64 {
        let dx = point[0] - self.obstacle_center[0];
        let dy = point[1] - self.obstacle_center[1];
        let dz = point[2] - self.obstacle_center[2];
        // Project onto axis
        let dot =
            dx * self.obstacle_axis[0] + dy * self.obstacle_axis[1] + dz * self.obstacle_axis[2];
        let px = dx - dot * self.obstacle_axis[0];
        let py = dy - dot * self.obstacle_axis[1];
        let pz = dz - dot * self.obstacle_axis[2];
        (px * px + py * py + pz * pz).sqrt()
    }

    /// Check if the straight line AB intersects the obstacle cylinder.
    pub fn needs_wrapping(&self) -> bool {
        // Simplified: check midpoint distance to cylinder
        let mid = [
            0.5 * (self.point_a[0] + self.point_b[0]),
            0.5 * (self.point_a[1] + self.point_b[1]),
            0.5 * (self.point_a[2] + self.point_b[2]),
        ];
        self._projected_distance(mid) < self.obstacle_radius * 1.1
    }

    /// Approximate total path length including wrapping arc.
    ///
    /// Returns `(total_length, wrap_angle)`.
    pub fn wrapped_path_length(&self) -> (f64, f64) {
        let da = self._projected_distance(self.point_a);
        let db = self._projected_distance(self.point_b);
        let r = self.obstacle_radius;

        if da <= r || db <= r {
            // Point is inside cylinder; return direct distance
            return (self.direct_distance(), 0.0);
        }

        // Tangent lengths
        let ta = (da * da - r * r).sqrt();
        let tb = (db * db - r * r).sqrt();

        // Half-angles
        let alpha_a = (r / da).asin();
        let alpha_b = (r / db).asin();

        // Wrap angle (approximate 2D projection)
        let wrap_angle = (PI - alpha_a - alpha_b).max(0.0);

        let arc_len = r * wrap_angle;
        (ta + arc_len + tb, wrap_angle)
    }

    /// Force direction at insertion point A due to tendon tension.
    ///
    /// Returns unit direction vector as \[f64; 3\].
    pub fn force_direction_a(&self) -> [f64; 3] {
        if !self.needs_wrapping() {
            // Direct A → B
            let d = self.direct_distance();
            if d < EPS {
                return [0.0, 0.0, 0.0];
            }
            return [
                (self.point_b[0] - self.point_a[0]) / d,
                (self.point_b[1] - self.point_a[1]) / d,
                (self.point_b[2] - self.point_a[2]) / d,
            ];
        }
        // When wrapping, direction is tangent from A to cylinder
        // Simplified: direction toward cylinder center projected
        let dx = self.obstacle_center[0] - self.point_a[0];
        let dy = self.obstacle_center[1] - self.point_a[1];
        let dz = self.obstacle_center[2] - self.point_a[2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        if d < EPS {
            return [0.0, 0.0, 0.0];
        }
        [dx / d, dy / d, dz / d]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MuscleTendonUnit
// ─────────────────────────────────────────────────────────────────────────────

/// Muscle-tendon unit coupling using a simplified Hill model.
///
/// The muscle-tendon unit consists of a contractile element (CE) in series
/// with a tendon (series elastic element, SEE). Force equilibrium requires
/// that muscle force equals tendon force.
#[derive(Debug, Clone)]
pub struct MuscleTendonUnit {
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fiber length (m).
    pub l_opt: f64,
    /// Tendon slack length (m).
    pub l_slack: f64,
    /// Pennation angle at optimal length (rad).
    pub pennation_angle: f64,
    /// Current activation level \[0, 1\].
    pub activation: f64,
    /// Current muscle fiber length (m).
    pub fiber_length: f64,
    /// Current tendon length (m).
    pub tendon_length: f64,
    /// Tendon stiffness constant (dimensionless).
    pub tendon_stiffness: f64,
    /// Maximum shortening velocity (l_opt / s).
    pub v_max: f64,
    /// Width of active force-length Gaussian.
    pub fl_width: f64,
}

impl MuscleTendonUnit {
    /// Create a new muscle-tendon unit.
    pub fn new(f_max: f64, l_opt: f64, l_slack: f64, pennation: f64) -> Self {
        Self {
            f_max,
            l_opt,
            l_slack,
            pennation_angle: pennation,
            activation: 0.0,
            fiber_length: l_opt,
            tendon_length: l_slack,
            tendon_stiffness: 35.0,
            v_max: 10.0,
            fl_width: 0.56,
        }
    }

    /// Default unit (generic muscle).
    pub fn default_unit() -> Self {
        Self::new(1000.0, 0.10, 0.20, 0.0)
    }

    /// Active force-length relationship (Gaussian).
    pub fn active_force_length(&self, l_norm: f64) -> f64 {
        let x = (l_norm - 1.0) / self.fl_width;
        (-x * x).exp()
    }

    /// Passive force-length relationship (exponential).
    pub fn passive_force_length(&self, l_norm: f64) -> f64 {
        if l_norm <= 1.0 {
            return 0.0;
        }
        let strain = l_norm - 1.0;
        let kpe = 4.0;
        let e0 = 0.6;
        ((-kpe * strain / e0).exp().recip() - 1.0) / ((kpe / e0).exp() - 1.0)
    }

    /// Tendon force-strain relationship.
    pub fn tendon_force(&self, tendon_strain: f64) -> f64 {
        if tendon_strain <= 0.0 {
            return 0.0;
        }
        self.f_max * self.tendon_stiffness * tendon_strain * tendon_strain
    }

    /// Total musculotendon length.
    pub fn total_length(&self) -> f64 {
        self.fiber_length * self.pennation_angle.cos() + self.tendon_length
    }

    /// Muscle force at current state.
    pub fn muscle_force(&self) -> f64 {
        let l_norm = self.fiber_length / self.l_opt;
        let f_active = self.activation * self.active_force_length(l_norm);
        let f_passive = self.passive_force_length(l_norm);
        self.f_max * (f_active + f_passive) * self.pennation_angle.cos()
    }

    /// Set activation level.
    pub fn set_activation(&mut self, a: f64) {
        self.activation = a.clamp(0.0, 1.0);
    }

    /// Simple equilibrium solve: given total length, find fiber/tendon split.
    pub fn solve_equilibrium(&mut self, total_len: f64) {
        // Iterate to find fiber_length such that muscle force = tendon force
        let cos_p = self.pennation_angle.cos();
        for _ in 0..50 {
            let tendon_len = total_len - self.fiber_length * cos_p;
            let tendon_strain = (tendon_len - self.l_slack) / self.l_slack;
            let f_tendon = self.tendon_force(tendon_strain);
            let f_muscle = self.muscle_force();
            let err = f_muscle - f_tendon;
            // Simple bisection-like adjustment
            if err.abs() < 0.01 {
                break;
            }
            let dl = if err > 0.0 { -0.001 } else { 0.001 };
            self.fiber_length += dl * self.l_opt;
            self.fiber_length = self.fiber_length.max(0.5 * self.l_opt);
        }
        self.tendon_length = total_len - self.fiber_length * cos_p;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Enthesis (bone-tendon insertion)
// ─────────────────────────────────────────────────────────────────────────────

/// Enthesis model: gradient properties at the bone-tendon insertion site.
///
/// The enthesis transitions from bone (stiff, mineralized) to tendon
/// (compliant, fibrous) over a graded zone. Material properties vary
/// along a normalized coordinate s ∈ \[0, 1\] where s=0 is bone and
/// s=1 is tendon.
#[derive(Debug, Clone)]
pub struct Enthesis {
    /// Bone modulus (Pa).
    pub bone_modulus: f64,
    /// Tendon modulus (Pa).
    pub tendon_modulus: f64,
    /// Transition zone length (m).
    pub zone_length: f64,
    /// Mineralization gradient exponent.
    pub gradient_exponent: f64,
}

impl Enthesis {
    /// Create a new enthesis model.
    pub fn new(bone_modulus: f64, tendon_modulus: f64, zone_length: f64) -> Self {
        Self {
            bone_modulus,
            tendon_modulus,
            zone_length,
            gradient_exponent: 2.0,
        }
    }

    /// Default enthesis for human tendon.
    pub fn default_human() -> Self {
        Self::new(20.0e9, 1.5e9, 0.002)
    }

    /// Modulus at normalized position s ∈ \[0, 1\].
    ///
    /// Uses a power-law interpolation.
    pub fn modulus_at(&self, s: f64) -> f64 {
        let s = s.clamp(0.0, 1.0);
        let log_e_bone = self.bone_modulus.ln();
        let log_e_tendon = self.tendon_modulus.ln();
        let log_e = log_e_bone + (log_e_tendon - log_e_bone) * s.powf(self.gradient_exponent);
        log_e.exp()
    }

    /// Mineralization fraction at position s.
    pub fn mineralization(&self, s: f64) -> f64 {
        let s = s.clamp(0.0, 1.0);
        1.0 - s.powf(self.gradient_exponent)
    }

    /// Stress concentration factor at the interface.
    ///
    /// Ratio of local stress to far-field tendon stress.
    pub fn stress_concentration(&self, s: f64) -> f64 {
        let e_local = self.modulus_at(s);
        e_local / self.tendon_modulus
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// HealingModel
// ─────────────────────────────────────────────────────────────────────────────

/// Collagen remodeling / healing model for damaged tendon.
///
/// Models the three phases of tendon healing:
/// 1. **Inflammatory** (0-7 days): minimal mechanical contribution
/// 2. **Proliferative** (7-21 days): collagen deposition, stiffness recovery
/// 3. **Remodeling** (21-365+ days): collagen alignment and maturation
///
/// The modulus recovers as: `E(t) = E_healed · h(t)`
/// where h(t) is the healing function.
#[derive(Debug, Clone)]
pub struct HealingModel {
    /// Healthy tendon modulus (Pa).
    pub healthy_modulus: f64,
    /// Inflammatory phase duration (days).
    pub inflammatory_duration: f64,
    /// Proliferative phase duration (days).
    pub proliferative_duration: f64,
    /// Remodeling time constant (days).
    pub remodeling_tau: f64,
    /// Current healing time (days).
    pub healing_time: f64,
    /// Initial damage level at injury (0 = undamaged, 1 = fully damaged).
    pub initial_damage: f64,
    /// Collagen alignment quality \[0, 1\].
    pub alignment_quality: f64,
}

impl HealingModel {
    /// Create a new healing model.
    pub fn new(healthy_modulus: f64, initial_damage: f64) -> Self {
        Self {
            healthy_modulus,
            inflammatory_duration: 7.0,
            proliferative_duration: 21.0,
            remodeling_tau: 90.0,
            healing_time: 0.0,
            initial_damage,
            alignment_quality: 0.0,
        }
    }

    /// Default healing model for tendon.
    pub fn default_tendon() -> Self {
        Self::new(DEFAULT_TENDON_MODULUS, 0.8)
    }

    /// Healing function h(t) ∈ \[0, 1\].
    pub fn healing_fraction(&self) -> f64 {
        let t = self.healing_time;
        if t < self.inflammatory_duration {
            // Inflammatory phase: minimal recovery
            let frac = t / self.inflammatory_duration;
            0.05 * frac // up to 5% during inflammation
        } else if t < self.inflammatory_duration + self.proliferative_duration {
            // Proliferative phase: rapid collagen deposition
            let t_prolif = t - self.inflammatory_duration;
            let frac = t_prolif / self.proliferative_duration;
            0.05 + 0.45 * frac // 5% to 50%
        } else {
            // Remodeling phase: slow maturation
            let t_remodel = t - self.inflammatory_duration - self.proliferative_duration;
            let base = 0.50;
            let remaining = 1.0 - base;
            base + remaining * (1.0 - (-t_remodel / self.remodeling_tau).exp())
        }
    }

    /// Current modulus based on healing state.
    pub fn current_modulus(&self) -> f64 {
        let undamaged_fraction = 1.0 - self.initial_damage;
        let healed_fraction = self.initial_damage * self.healing_fraction();
        self.healthy_modulus * (undamaged_fraction + healed_fraction)
    }

    /// Advance healing by a number of days.
    pub fn advance(&mut self, days: f64) {
        self.healing_time += days;
        // Update alignment quality (improves during remodeling)
        let t_remodel =
            (self.healing_time - self.inflammatory_duration - self.proliferative_duration).max(0.0);
        self.alignment_quality =
            (1.0 - (-t_remodel / (self.remodeling_tau * 2.0)).exp()).clamp(0.0, 1.0);
    }

    /// Phase description at current time.
    pub fn current_phase(&self) -> &'static str {
        if self.healing_time < self.inflammatory_duration {
            "inflammatory"
        } else if self.healing_time < self.inflammatory_duration + self.proliferative_duration {
            "proliferative"
        } else {
            "remodeling"
        }
    }

    /// Estimated time to reach a target modulus fraction (days).
    pub fn time_to_recovery(&self, target_fraction: f64) -> f64 {
        // Binary search
        let mut lo = 0.0_f64;
        let mut hi = 1000.0_f64;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            let saved = self.healing_time;
            let mut temp = self.clone();
            temp.healing_time = mid;
            let frac = temp.healing_fraction();
            let _ = saved; // restore original conceptually
            if frac < target_fraction {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PreTension
// ─────────────────────────────────────────────────────────────────────────────

/// Models pre-tension (initial strain) in a tendon.
///
/// In vivo, tendons are typically under 2-4% pre-strain even at rest,
/// which places them in the toe-to-linear transition zone for rapid
/// force transmission.
#[derive(Debug, Clone)]
pub struct PreTension {
    /// Pre-tension strain (dimensionless).
    pub pretension_strain: f64,
    /// Corresponding pre-tension stress (Pa).
    pub pretension_stress: f64,
    /// Rest length under pre-tension (shorter than slack length).
    pub pretensioned_length: f64,
    /// Slack length (m).
    pub slack_length: f64,
}

impl PreTension {
    /// Create from material properties.
    pub fn new(props: &TendonMaterialProperties) -> Self {
        let pretension_stress = props.elastic_modulus * props.pretension_strain;
        let pretensioned_length = props.rest_length * (1.0 + props.pretension_strain);
        Self {
            pretension_strain: props.pretension_strain,
            pretension_stress,
            pretensioned_length,
            slack_length: props.rest_length,
        }
    }

    /// Effective strain accounting for pre-tension.
    pub fn effective_strain(&self, total_strain: f64) -> f64 {
        total_strain + self.pretension_strain
    }

    /// Force due to pre-tension alone (N).
    pub fn pretension_force(&self, cross_section: f64) -> f64 {
        self.pretension_stress * cross_section
    }

    /// Whether the tendon is slack (below pre-tension).
    pub fn is_slack(&self, current_length: f64) -> bool {
        current_length < self.slack_length
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TendonSimulation (integrated simulation)
// ─────────────────────────────────────────────────────────────────────────────

/// Integrated tendon simulation combining all models.
#[derive(Debug, Clone)]
pub struct TendonSimulation {
    /// Material properties.
    pub props: TendonMaterialProperties,
    /// Stress-strain model.
    pub curve: TendonStressStrainCurve,
    /// QLV model.
    pub qlv: QLVModel,
    /// Fatigue damage.
    pub fatigue: FatigueDamage,
    /// Healing state.
    pub healing: HealingModel,
    /// Current strain.
    pub current_strain: f64,
    /// Current stress (Pa).
    pub current_stress: f64,
    /// Time elapsed (s).
    pub time: f64,
}

impl TendonSimulation {
    /// Create a new integrated tendon simulation.
    pub fn new(props: TendonMaterialProperties) -> Self {
        let curve = TendonStressStrainCurve::new(props.clone());
        let qlv = QLVModel::new(FungStrainEnergy::default_tendon());
        let healing = HealingModel::new(props.elastic_modulus, 0.0);
        Self {
            props,
            curve,
            qlv,
            fatigue: FatigueDamage::new(),
            healing,
            current_strain: 0.0,
            current_stress: 0.0,
            time: 0.0,
        }
    }

    /// Apply a strain increment and advance time.
    pub fn step(&mut self, strain: f64, dt: f64) {
        self.current_strain = strain;
        self.current_stress = self.qlv.update(strain, dt);
        self.curve.damage = self.fatigue.damage;
        self.time += dt;
    }

    /// Get the current force (N).
    pub fn current_force(&self) -> f64 {
        self.current_stress * self.props.cross_section_area
    }

    /// Apply fatigue cycles.
    pub fn apply_fatigue(&mut self, n_cycles: u64, amplitude: f64) {
        self.fatigue.apply_cycles(n_cycles, amplitude);
    }

    /// Check if ruptured.
    pub fn is_ruptured(&self) -> bool {
        self.fatigue.is_ruptured()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1.0e-6;

    // 1. TendonMaterialProperties: stiffness calculation
    #[test]
    fn test_material_stiffness() {
        let props = TendonMaterialProperties::new();
        let k = props.stiffness();
        let expected = props.elastic_modulus * props.cross_section_area / props.rest_length;
        assert!(
            (k - expected).abs() < TOL,
            "Stiffness mismatch: {k} vs {expected}"
        );
    }

    // 2. Mass calculation
    #[test]
    fn test_material_mass() {
        let props = TendonMaterialProperties::new();
        let m = props.mass();
        assert!(m > 0.0, "Mass should be positive");
        let expected = props.density * props.cross_section_area * props.rest_length;
        assert!((m - expected).abs() < TOL);
    }

    // 3. Achilles tendon preset
    #[test]
    fn test_achilles_preset() {
        let props = TendonMaterialProperties::achilles();
        assert!(props.elastic_modulus > 0.0);
        assert!(props.failure_strain > props.toe_limit);
    }

    // 4. CrimpModel: angle at zero strain
    #[test]
    fn test_crimp_zero_strain() {
        let crimp = CrimpModel::new(0.15, 0.03);
        let angle = crimp.angle_at_strain(0.0);
        assert!((angle - 0.15).abs() < TOL, "Angle at zero strain");
    }

    // 5. CrimpModel: angle at toe limit is zero
    #[test]
    fn test_crimp_toe_limit() {
        let crimp = CrimpModel::new(0.15, 0.03);
        let angle = crimp.angle_at_strain(0.03);
        assert!(angle.abs() < TOL, "Angle at toe limit should be 0");
    }

    // 6. CrimpModel: recruitment fraction
    #[test]
    fn test_crimp_recruitment() {
        let crimp = CrimpModel::new(0.15, 0.03);
        assert!(crimp.recruitment_fraction(0.0) < TOL);
        assert!((crimp.recruitment_fraction(0.03) - 1.0).abs() < TOL);
        assert!((crimp.recruitment_fraction(0.015) - 0.5).abs() < TOL);
    }

    // 7. FungStrainEnergy: zero strain gives zero energy
    #[test]
    fn test_fung_zero() {
        let fung = FungStrainEnergy::default_tendon();
        assert!(fung.energy(0.0).abs() < TOL);
        assert!(fung.stress(0.0).abs() < TOL);
    }

    // 8. FungStrainEnergy: positive strain gives positive stress
    #[test]
    fn test_fung_positive() {
        let fung = FungStrainEnergy::default_tendon();
        let sigma = fung.stress(0.05);
        assert!(sigma > 0.0, "Positive strain should give positive stress");
    }

    // 9. FungStrainEnergy: tangent modulus positive
    #[test]
    fn test_fung_tangent() {
        let fung = FungStrainEnergy::default_tendon();
        let e = fung.tangent_modulus(0.02);
        assert!(e > 0.0, "Tangent modulus should be positive");
    }

    // 10. FungStrainEnergy: inverse stress roundtrip
    #[test]
    fn test_fung_inverse() {
        let fung = FungStrainEnergy::default_tendon();
        let strain_in = 0.04;
        let sigma = fung.stress(strain_in);
        let strain_out = fung.inverse_stress(sigma);
        assert!(
            (strain_in - strain_out).abs() < 1.0e-4,
            "Inverse should recover strain: {strain_in} vs {strain_out}"
        );
    }

    // 11. TendonStressStrainCurve: zero stress at zero strain
    #[test]
    fn test_curve_zero() {
        let curve = TendonStressStrainCurve::new(TendonMaterialProperties::new());
        assert!(curve.stress(0.0).abs() < TOL);
    }

    // 12. TendonStressStrainCurve: toe region is softer than linear
    #[test]
    fn test_curve_toe_vs_linear() {
        let curve = TendonStressStrainCurve::new(TendonMaterialProperties::new());
        let s_toe = curve.stress(0.01);
        let s_lin = curve.stress(0.05);
        // Stress per strain should be lower in toe
        let mod_toe = s_toe / 0.01;
        let mod_lin = s_lin / 0.05;
        assert!(mod_toe < mod_lin, "Toe region should be softer");
    }

    // 13. TendonStressStrainCurve: failure returns zero
    #[test]
    fn test_curve_failure() {
        let curve = TendonStressStrainCurve::new(TendonMaterialProperties::new());
        let s = curve.stress(0.15); // beyond failure
        assert!(s.abs() < TOL, "Beyond failure should be zero");
    }

    // 14. TendonStressStrainCurve: stored energy is positive
    #[test]
    fn test_curve_energy() {
        let curve = TendonStressStrainCurve::new(TendonMaterialProperties::new());
        let w = curve.stored_energy(0.05);
        assert!(w > 0.0, "Stored energy should be positive");
    }

    // 15. QLVModel: relaxation function at t=0 is 1.0
    #[test]
    fn test_qlv_g_at_zero() {
        let qlv = QLVModel::new(FungStrainEnergy::default_tendon());
        let g0 = qlv.relaxation_function(0.0);
        assert!((g0 - 1.0).abs() < TOL, "G(0) should be 1.0, got {g0}");
    }

    // 16. QLVModel: relaxation function decreases with time
    #[test]
    fn test_qlv_g_decreases() {
        let qlv = QLVModel::new(FungStrainEnergy::default_tendon());
        let g0 = qlv.relaxation_function(0.0);
        let g10 = qlv.relaxation_function(10.0);
        let g100 = qlv.relaxation_function(100.0);
        assert!(g10 < g0, "G(10) < G(0)");
        assert!(g100 < g10, "G(100) < G(10)");
    }

    // 17. QLVModel: stress relaxation test
    #[test]
    fn test_qlv_relaxation() {
        let mut qlv = QLVModel::new(FungStrainEnergy::default_tendon());
        let result = qlv.stress_relaxation_test(0.03, 100.0, 200);
        assert!(!result.is_empty());
        // Stress should decrease over time
        let first_stress = result[1].1;
        let last_stress = result.last().unwrap().1;
        assert!(
            last_stress < first_stress || first_stress.abs() < TOL,
            "Stress should relax: first={first_stress}, last={last_stress}"
        );
    }

    // 18. CreepModel: strain increases with time
    #[test]
    fn test_creep_increases() {
        let creep = CreepModel::default_tendon(0.03);
        let e0 = creep.total_strain(0.0);
        let e10 = creep.total_strain(10.0);
        let e100 = creep.total_strain(100.0);
        assert!(e10 > e0, "Strain should increase with time");
        assert!(e100 > e10, "Strain should continue increasing");
    }

    // 19. CreepModel: step advances correctly
    #[test]
    fn test_creep_step() {
        let mut creep = CreepModel::default_tendon(0.03);
        creep.step(1.0);
        assert!(creep.elapsed_time > 0.0);
        assert!(creep.creep_strain > 0.0);
    }

    // 20. StressRelaxation: stress at t=0 equals initial
    #[test]
    fn test_relaxation_initial() {
        let sr = StressRelaxation::new(100.0e6);
        let s0 = sr.stress_at(0.0);
        assert!(
            (s0 - 100.0e6).abs() < TOL,
            "Stress at t=0 should be initial"
        );
    }

    // 21. StressRelaxation: stress decreases with time
    #[test]
    fn test_relaxation_decreases() {
        let sr = StressRelaxation::new(100.0e6);
        let s0 = sr.stress_at(0.0);
        let s100 = sr.stress_at(100.0);
        assert!(s100 < s0, "Stress should decrease");
    }

    // 22. StrainRateDependence: higher rate gives higher modulus
    #[test]
    fn test_rate_dependence() {
        let srd = StrainRateDependence::default_tendon();
        let e_slow = srd.effective_modulus(0.001);
        let e_fast = srd.effective_modulus(1.0);
        assert!(e_fast > e_slow, "Fast rate should give higher modulus");
    }

    // 23. StrainRateDependence: enhancement at ref rate
    #[test]
    fn test_rate_enhancement() {
        let srd = StrainRateDependence::default_tendon();
        let f = srd.enhancement_factor(srd.ref_rate);
        assert!(f > 1.0, "Enhancement at ref rate should be > 1");
    }

    // 24. FatigueDamage: no damage below fatigue limit
    #[test]
    fn test_fatigue_below_limit() {
        let mut fd = FatigueDamage::new();
        fd.apply_cycles(1_000_000, 0.01); // below fatigue limit
        assert!(fd.damage < TOL, "No damage below fatigue limit");
    }

    // 25. FatigueDamage: damage accumulates above limit
    #[test]
    fn test_fatigue_above_limit() {
        let mut fd = FatigueDamage::new();
        fd.apply_cycles(100_000, 0.05);
        assert!(fd.damage > 0.0, "Damage should accumulate");
    }

    // 26. FiberReinforcedModel: ground stress at zero is zero
    #[test]
    fn test_fiber_ground_zero() {
        let model = FiberReinforcedModel::default_tendon();
        let s = model.ground_stress(0.0);
        assert!(s.abs() < TOL);
    }

    // 27. FiberReinforcedModel: fiber stress only in tension
    #[test]
    fn test_fiber_tension_only() {
        let model = FiberReinforcedModel::default_tendon();
        let s_comp = model.fiber_stress(0.9); // lambda < 1 → compression
        assert!(s_comp.abs() < TOL, "No fiber stress in compression");
        let s_tens = model.fiber_stress(1.1); // lambda > 1 → tension
        assert!(s_tens > 0.0, "Fiber stress in tension should be positive");
    }

    // 28. TendonWrapping: direct distance
    #[test]
    fn test_wrapping_distance() {
        let tw = TendonWrapping::new(
            [0.0, 0.0, 0.0],
            0.01,
            [0.0, 0.0, 1.0],
            [-0.05, 0.0, 0.0],
            [0.05, 0.0, 0.0],
        );
        assert!((tw.direct_distance() - 0.1).abs() < TOL);
    }

    // 29. MuscleTendonUnit: zero activation → zero active force
    #[test]
    fn test_mtu_zero_activation() {
        let mtu = MuscleTendonUnit::default_unit();
        let f = mtu.muscle_force();
        // At zero activation, only passive force (which is 0 at l_opt)
        assert!(f.abs() < 1.0, "Zero activation should give near-zero force");
    }

    // 30. MuscleTendonUnit: full activation at optimal length
    #[test]
    fn test_mtu_full_activation() {
        let mut mtu = MuscleTendonUnit::default_unit();
        mtu.set_activation(1.0);
        let f = mtu.muscle_force();
        assert!(f > 0.0, "Full activation should give positive force");
    }

    // 31. Enthesis: bone end is stiffer than tendon end
    #[test]
    fn test_enthesis_gradient() {
        let enth = Enthesis::default_human();
        let e_bone = enth.modulus_at(0.0);
        let e_tendon = enth.modulus_at(1.0);
        assert!(e_bone > e_tendon, "Bone end should be stiffer");
    }

    // 32. Enthesis: mineralization at bone end is 1
    #[test]
    fn test_enthesis_mineralization() {
        let enth = Enthesis::default_human();
        let m0 = enth.mineralization(0.0);
        let m1 = enth.mineralization(1.0);
        assert!((m0 - 1.0).abs() < TOL, "Bone end fully mineralized");
        assert!(m1.abs() < TOL, "Tendon end unmineralized");
    }

    // 33. HealingModel: initial modulus reflects damage
    #[test]
    fn test_healing_initial() {
        let heal = HealingModel::new(1.5e9, 0.8);
        let e = heal.current_modulus();
        let expected = 1.5e9 * (1.0 - 0.8 + 0.8 * heal.healing_fraction());
        assert!(
            (e - expected).abs() < 1.0,
            "Initial modulus should reflect damage"
        );
    }

    // 34. HealingModel: modulus increases with healing
    #[test]
    fn test_healing_recovery() {
        let mut heal = HealingModel::default_tendon();
        let e0 = heal.current_modulus();
        heal.advance(60.0);
        let e60 = heal.current_modulus();
        assert!(e60 > e0, "Modulus should increase with healing");
    }

    // 35. PreTension: effective strain includes pretension
    #[test]
    fn test_pretension_effective() {
        let props = TendonMaterialProperties::new();
        let pt = PreTension::new(&props);
        let eff = pt.effective_strain(0.01);
        assert!((eff - (0.01 + props.pretension_strain)).abs() < TOL);
    }

    // 36. TendonSimulation: step doesn't panic
    #[test]
    fn test_simulation_step() {
        let mut sim = TendonSimulation::new(TendonMaterialProperties::new());
        for i in 0..100 {
            let strain = 0.0001 * i as f64;
            sim.step(strain, 0.01);
        }
        assert!(sim.current_stress.is_finite());
        assert!(sim.time > 0.0);
    }

    // 37. TendonSimulation: force is positive under tension
    #[test]
    fn test_simulation_force() {
        let mut sim = TendonSimulation::new(TendonMaterialProperties::new());
        sim.step(0.05, 0.01);
        let f = sim.current_force();
        assert!(f > 0.0, "Force under tension should be positive");
    }

    // 38. Curve generation produces correct number of points
    #[test]
    fn test_curve_generation() {
        let curve = TendonStressStrainCurve::new(TendonMaterialProperties::new());
        let data = curve.generate_curve(50);
        assert_eq!(data.len(), 50);
    }

    // 39. Patellar tendon preset
    #[test]
    fn test_patellar_preset() {
        let props = TendonMaterialProperties::patellar();
        assert!(props.rest_length > 0.0);
        assert!(props.failure_strain > 0.0);
    }

    // 40. Ligament ACL preset
    #[test]
    fn test_acl_preset() {
        let props = TendonMaterialProperties::ligament_acl();
        assert!(
            props.elastic_modulus < DEFAULT_TENDON_MODULUS,
            "ACL is softer than tendon"
        );
    }
}
