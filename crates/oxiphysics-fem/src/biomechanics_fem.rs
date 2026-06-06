// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Biomechanics FEM: musculoskeletal modeling, cardiac mechanics.
//!
//! This module provides finite element models for biological tissues including:
//! - Hill muscle model with force-length and force-velocity relationships
//! - Musculoskeletal multi-muscle systems with redundancy resolution
//! - Cardiac mechanics (Holzapfel-Ogden + Guccione active stress)
//! - Bone remodeling via Wolff's law
//! - Ligament constitutive model (toe + linear regions)
//! - Artery pressurization (Laplace law)
//! - Breathing mechanics (lung compliance)
//! - Spinal segment modeling

// ────────────────────────────────────────────────────────────────────────────
// MuscleFiber – Hill three-element muscle model
// ────────────────────────────────────────────────────────────────────────────

/// Hill muscle model constants for a single fiber bundle.
pub struct MuscleFiberParams {
    /// Maximum isometric force \[N\]
    pub f_max: f64,
    /// Optimal fiber length \[m\]
    pub l_opt: f64,
    /// Maximum shortening velocity \[m/s\]
    pub v_max: f64,
    /// Force-velocity shape constant (Hill constant)
    pub k_hill: f64,
    /// Passive stiffness coefficient \[N/m\]
    pub k_passive: f64,
    /// Activation time constant \[s\]
    pub tau_act: f64,
    /// Deactivation time constant \[s\]
    pub tau_deact: f64,
}

impl Default for MuscleFiberParams {
    fn default() -> Self {
        Self {
            f_max: 1000.0,
            l_opt: 0.10,
            v_max: 0.50,
            k_hill: 0.25,
            k_passive: 5000.0,
            tau_act: 0.010,
            tau_deact: 0.040,
        }
    }
}

/// State of a single Hill muscle fiber.
pub struct MuscleFiber {
    /// Current fiber length \[m\]
    pub length: f64,
    /// Current fiber velocity \[m/s\] (positive = lengthening)
    pub velocity: f64,
    /// Neural activation signal u ∈ \[0, 1\]
    pub excitation: f64,
    /// Activation state a ∈ \[0, 1\]
    pub activation: f64,
    /// Model parameters
    pub params: MuscleFiberParams,
}

impl MuscleFiber {
    /// Create a new `MuscleFiber` with given parameters and initial state.
    pub fn new(params: MuscleFiberParams) -> Self {
        let l0 = params.l_opt;
        Self {
            length: l0,
            velocity: 0.0,
            excitation: 0.0,
            activation: 0.0,
            params,
        }
    }

    /// Normalized fiber length λ = l / l_opt.
    pub fn normalized_length(&self) -> f64 {
        self.length / self.params.l_opt
    }

    /// Force-length multiplier fl(λ): Gaussian centred at λ=1.
    ///
    /// fl(λ) = exp(-((λ - 1) / 0.45)²)
    pub fn force_length(&self) -> f64 {
        let lambda = self.normalized_length();
        let width = 0.45_f64;
        (-((lambda - 1.0) / width).powi(2)).exp()
    }

    /// Force-velocity multiplier fv(v).
    ///
    /// For shortening (v < 0): fv = (v_max + v) / (v_max + k*|v|)  \[Hill concentric\]
    /// For lengthening (v ≥ 0): fv capped at 1.8 (eccentric)
    pub fn force_velocity(&self) -> f64 {
        let v = self.velocity;
        let v_max = self.params.v_max;
        let k = self.params.k_hill;
        if v <= 0.0 {
            // shortening or isometric
            let num = v_max - v.abs();
            let den = v_max + k * v.abs();
            if den <= 0.0 {
                return 0.0;
            }
            (num / den).max(0.0)
        } else {
            // lengthening – eccentric force exceeds isometric
            // Simplified: fv = 1 + 0.8 * (v / v_max), capped at 1.8
            let factor = 1.0 + 0.8 * (v / v_max);
            factor.min(1.8)
        }
    }

    /// Passive elastic force \[N\].
    pub fn passive_force(&self) -> f64 {
        let lambda = self.normalized_length();
        if lambda <= 1.0 {
            0.0
        } else {
            self.params.k_passive * (lambda - 1.0).powi(2)
        }
    }

    /// Total muscle force \[N\] = active + passive.
    pub fn total_force(&self) -> f64 {
        let f_active =
            self.params.f_max * self.activation * self.force_length() * self.force_velocity();
        let f_passive = self.passive_force();
        f_active + f_passive
    }

    /// Advance activation dynamics: da/dt = (u - a) / τ  (Zajac model).
    pub fn update_activation(&mut self, dt: f64) {
        let u = self.excitation.clamp(0.0, 1.0);
        let a = self.activation;
        let tau = if u > a {
            self.params.tau_act
        } else {
            self.params.tau_deact
        };
        let da = (u - a) / tau;
        self.activation = (a + da * dt).clamp(0.0, 1.0);
    }
}

// ────────────────────────────────────────────────────────────────────────────
// MusculoskeletalModel – multi-muscle system around a joint
// ────────────────────────────────────────────────────────────────────────────

/// A single muscle entry in the musculoskeletal model.
pub struct MuscleEntry {
    /// Maximum isometric force \[N\]
    pub f_max: f64,
    /// Moment arm \[m\] (positive = flexion)
    pub moment_arm: f64,
    /// Current normalized force f/f_max ∈ \[0, 1\]
    pub normalized_force: f64,
}

/// Multi-muscle redundancy model for a single-axis joint.
pub struct MusculoskeletalModel {
    /// List of muscles
    pub muscles: Vec<MuscleEntry>,
    /// Required joint torque \[N·m\]
    pub required_torque: f64,
    /// Optimization exponent (static optimization criterion)
    pub n_exp: f64,
}

impl MusculoskeletalModel {
    /// Create a new `MusculoskeletalModel`.
    pub fn new(muscles: Vec<MuscleEntry>, required_torque: f64, n_exp: f64) -> Self {
        Self {
            muscles,
            required_torque,
            n_exp,
        }
    }

    /// Static optimization: minimise Σ(f_m/f_max)^n subject to torque balance.
    ///
    /// Uses proportional sharing weighted by moment arms and f_max (simplified).
    pub fn solve_static(&mut self) {
        // Simple proportional distribution: weight = moment_arm * f_max
        let weights: Vec<f64> = self
            .muscles
            .iter()
            .map(|m| (m.moment_arm * m.f_max).max(0.0))
            .collect();
        let w_sum: f64 = weights.iter().sum();
        if w_sum <= 0.0 {
            return;
        }
        for (m, w) in self.muscles.iter_mut().zip(weights.iter()) {
            let f = self.required_torque * w / (w_sum * m.moment_arm.max(1e-9));
            m.normalized_force = (f / m.f_max).clamp(0.0, 1.0);
        }
    }

    /// Compute total joint torque from current muscle forces \[N·m\].
    pub fn joint_torque(&self) -> f64 {
        self.muscles
            .iter()
            .map(|m| m.normalized_force * m.f_max * m.moment_arm)
            .sum()
    }

    /// Objective function value Σ(f_m/f_max)^n.
    pub fn objective(&self) -> f64 {
        self.muscles
            .iter()
            .map(|m| m.normalized_force.powf(self.n_exp))
            .sum()
    }
}

// ────────────────────────────────────────────────────────────────────────────
// CardiacFem – active stress in cardiac tissue
// ────────────────────────────────────────────────────────────────────────────

/// Parameters for Guccione active stress model.
pub struct CardiacParams {
    /// Maximum active stress \[Pa\]
    pub t_max: f64,
    /// Calcium concentration at half-maximum activation \[μM\]
    pub ca_50: f64,
    /// Hill coefficient for calcium activation
    pub n_hill: f64,
    /// Fiber stretch at which maximum force is generated
    pub lambda_opt: f64,
    /// Passive stiffness scale C \[Pa\]
    pub c_passive: f64,
    /// Holzapfel-Ogden fiber stiffness k1 \[Pa\]
    pub k1: f64,
    /// Holzapfel-Ogden fiber exponential k2 \[-\]
    pub k2: f64,
}

impl Default for CardiacParams {
    fn default() -> Self {
        Self {
            t_max: 135_000.0,
            ca_50: 0.715,
            n_hill: 1.9,
            lambda_opt: 1.1,
            c_passive: 940.0,
            k1: 2360.0,
            k2: 10.8,
        }
    }
}

/// Cardiac tissue FEM state at a single integration point.
pub struct CardiacFem {
    /// Fiber stretch ratio λ_f
    pub lambda_fiber: f64,
    /// Cytosolic calcium concentration \[μM\]
    pub calcium: f64,
    /// Active stress along fiber \[Pa\]
    pub active_stress: f64,
    /// Passive stress along fiber \[Pa\]
    pub passive_stress: f64,
    /// Model parameters
    pub params: CardiacParams,
}

impl CardiacFem {
    /// Create a `CardiacFem` at reference state.
    pub fn new(params: CardiacParams) -> Self {
        Self {
            lambda_fiber: 1.0,
            calcium: 0.0,
            active_stress: 0.0,
            passive_stress: 0.0,
            params,
        }
    }

    /// Calcium activation function f_Ca = Ca^n / (Ca^n + Ca50^n).
    pub fn calcium_activation(&self) -> f64 {
        let ca = self.calcium.max(0.0);
        let ca50 = self.params.ca_50;
        let n = self.params.n_hill;
        if ca <= 0.0 {
            return 0.0;
        }
        let can = ca.powf(n);
        let ca50n = ca50.powf(n);
        can / (can + ca50n)
    }

    /// Force-length modifier for cardiac fiber (parabolic).
    pub fn force_length(&self) -> f64 {
        let lam = self.lambda_fiber;
        let lopt = self.params.lambda_opt;
        let val = 1.0 - 4.0 * (lam - lopt).powi(2);
        val.max(0.0)
    }

    /// Update active and passive stresses at current state.
    pub fn update_stress(&mut self) {
        let f_ca = self.calcium_activation();
        let fl = self.force_length();
        self.active_stress = self.params.t_max * f_ca * fl;

        // Holzapfel-Ogden simplified fiber term: σ_f = 2*k1*(λ²-1)*exp(k2*(λ²-1)²)
        let lam2 = self.lambda_fiber * self.lambda_fiber;
        let e4 = lam2 - 1.0;
        self.passive_stress = 2.0 * self.params.k1 * e4 * (self.params.k2 * e4 * e4).exp();
    }

    /// Total Cauchy stress along fiber direction \[Pa\].
    pub fn total_stress(&self) -> f64 {
        self.active_stress + self.passive_stress
    }
}

// ────────────────────────────────────────────────────────────────────────────
// BoneRemodeling – Wolff's law FEM
// ────────────────────────────────────────────────────────────────────────────

/// Parameters for bone remodeling.
pub struct BoneRemodelingParams {
    /// Remodeling rate coefficient B \[kg/(m³·s·Pa)\]
    pub b_rate: f64,
    /// Reference stimulus strain energy density \[Pa\]
    pub sigma_ref: f64,
    /// Dead zone half-width around σ_ref \[-\]
    pub dead_zone: f64,
    /// Minimum density \[kg/m³\]
    pub rho_min: f64,
    /// Maximum density \[kg/m³\]
    pub rho_max: f64,
}

impl Default for BoneRemodelingParams {
    fn default() -> Self {
        Self {
            b_rate: 1e-4,
            sigma_ref: 3500.0,
            dead_zone: 0.075,
            rho_min: 100.0,
            rho_max: 2000.0,
        }
    }
}

/// Single bone voxel subject to Wolff's remodeling law.
pub struct BoneRemodeling {
    /// Current apparent density \[kg/m³\]
    pub density: f64,
    /// Current effective stress stimulus \[Pa\]
    pub stress_stimulus: f64,
    /// Model parameters
    pub params: BoneRemodelingParams,
}

impl BoneRemodeling {
    /// Create a `BoneRemodeling` element at initial density.
    pub fn new(initial_density: f64, params: BoneRemodelingParams) -> Self {
        Self {
            density: initial_density,
            stress_stimulus: 0.0,
            params,
        }
    }

    /// Remodeling signal S = σ_e/σ_ref - 1 (outside dead zone).
    pub fn remodeling_signal(&self) -> f64 {
        let s = self.stress_stimulus / self.params.sigma_ref - 1.0;
        let dz = self.params.dead_zone;
        if s > dz {
            s - dz
        } else if s < -dz {
            s + dz
        } else {
            0.0
        }
    }

    /// Advance density: dρ/dt = B * S.
    pub fn update_density(&mut self, dt: f64) {
        let signal = self.remodeling_signal();
        let drho = self.params.b_rate * signal * dt;
        self.density = (self.density + drho).clamp(self.params.rho_min, self.params.rho_max);
    }

    /// Elastic modulus from density via power law: E = c * ρ^γ \[Pa\].
    ///
    /// Uses Currey relation E = 3.79e-3 * ρ^3.
    pub fn elastic_modulus(&self) -> f64 {
        3.79e-3 * self.density.powi(3)
    }
}

// ────────────────────────────────────────────────────────────────────────────
// LigamentModel – toe and linear region, transversely isotropic
// ────────────────────────────────────────────────────────────────────────────

/// Ligament constitutive model (toe + linear fiber recruitment).
pub struct LigamentModel {
    /// Ligament cross-sectional area \[m²\]
    pub area: f64,
    /// Toe-region limit strain εt \[-\]
    pub epsilon_toe: f64,
    /// Stiffness in linear region \[Pa\]
    pub e_linear: f64,
    /// Toe-region curvature parameter \[-\]
    pub c_toe: f64,
    /// Current strain ε
    pub strain: f64,
}

impl LigamentModel {
    /// Create a new `LigamentModel`.
    pub fn new(area: f64, epsilon_toe: f64, e_linear: f64, c_toe: f64) -> Self {
        Self {
            area,
            epsilon_toe,
            e_linear,
            c_toe,
            strain: 0.0,
        }
    }

    /// Stress at current strain \[Pa\].
    ///
    /// Compression → 0; toe region → exponential; linear region → linear.
    pub fn stress(&self) -> f64 {
        let eps = self.strain;
        if eps <= 0.0 {
            // ligaments carry no compressive load
            return 0.0;
        }
        if eps <= self.epsilon_toe {
            // toe region: exponential recruitment
            let ratio = eps / self.epsilon_toe;
            self.e_linear * self.epsilon_toe * (self.c_toe * ratio).exp() / self.c_toe
        } else {
            // linear region
            let sigma_toe = self.e_linear * self.epsilon_toe * self.c_toe.exp() / self.c_toe;
            sigma_toe + self.e_linear * (eps - self.epsilon_toe)
        }
    }

    /// Force = stress × area \[N\].
    pub fn force(&self) -> f64 {
        self.stress() * self.area
    }
}

// ────────────────────────────────────────────────────────────────────────────
// ArteryFem – pressurized thick-walled cylinder
// ────────────────────────────────────────────────────────────────────────────

/// Thick-walled arterial tube model (Lamé + Laplace approximation).
pub struct ArteryFem {
    /// Inner radius \[m\]
    pub r_inner: f64,
    /// Outer radius \[m\]
    pub r_outer: f64,
    /// Internal blood pressure \[Pa\]
    pub pressure: f64,
    /// Wall elastic modulus \[Pa\]
    pub elastic_modulus: f64,
    /// Poisson ratio \[-\]
    pub poisson: f64,
}

impl ArteryFem {
    /// Create an `ArteryFem` model.
    pub fn new(
        r_inner: f64,
        r_outer: f64,
        pressure: f64,
        elastic_modulus: f64,
        poisson: f64,
    ) -> Self {
        Self {
            r_inner,
            r_outer,
            pressure,
            elastic_modulus,
            poisson,
        }
    }

    /// Wall thickness \[m\].
    pub fn wall_thickness(&self) -> f64 {
        self.r_outer - self.r_inner
    }

    /// Laplace-law hoop stress σ_θ = p*r_inner / h \[Pa\].
    pub fn laplace_hoop_stress(&self) -> f64 {
        let h = self.wall_thickness();
        if h <= 0.0 {
            return 0.0;
        }
        self.pressure * self.r_inner / h
    }

    /// Lamé exact hoop stress at inner radius \[Pa\].
    pub fn lame_hoop_stress_inner(&self) -> f64 {
        let ri2 = self.r_inner * self.r_inner;
        let ro2 = self.r_outer * self.r_outer;
        let denom = ro2 - ri2;
        if denom <= 0.0 {
            return 0.0;
        }
        self.pressure * ri2 * (ri2 + ro2) / (denom * ri2)
    }

    /// Radial displacement at inner wall \[m\] (thick-wall elasticity).
    pub fn inner_radial_displacement(&self) -> f64 {
        let ri = self.r_inner;
        let ro = self.r_outer;
        let ri2 = ri * ri;
        let ro2 = ro * ro;
        let denom = ro2 - ri2;
        if denom <= 0.0 {
            return 0.0;
        }
        self.pressure * ri * (1.0 - self.poisson * ri2 / ro2 + ri2 / ro2)
            / (self.elastic_modulus * denom / ri2)
    }
}

// ────────────────────────────────────────────────────────────────────────────
// BreathingMechanics – lung parenchyma compliance and breathing work
// ────────────────────────────────────────────────────────────────────────────

/// Lung parenchyma breathing mechanics model.
pub struct BreathingMechanics {
    /// Static lung compliance \[m³/Pa\] (~200 mL/cmH2O)
    pub compliance: f64,
    /// Airway resistance \[Pa·s/m³\]
    pub resistance: f64,
    /// Functional residual capacity \[m³\]
    pub frc: f64,
    /// Current lung volume \[m³\]
    pub volume: f64,
    /// Current pleural pressure \[Pa\]
    pub pleural_pressure: f64,
}

impl BreathingMechanics {
    /// Create a `BreathingMechanics` model at FRC.
    pub fn new(compliance: f64, resistance: f64, frc: f64) -> Self {
        Self {
            compliance,
            resistance,
            frc,
            volume: frc,
            pleural_pressure: -500.0, // ~-5 cmH2O
        }
    }

    /// Tidal volume from pressure change ΔP \[m³\].
    pub fn tidal_volume(&self, delta_pressure: f64) -> f64 {
        self.compliance * delta_pressure.abs()
    }

    /// Elastic recoil pressure \[Pa\].
    pub fn recoil_pressure(&self) -> f64 {
        (self.volume - self.frc) / self.compliance
    }

    /// Resistive pressure drop for given flow rate \[Pa\].
    pub fn resistive_pressure(&self, flow_rate: f64) -> f64 {
        self.resistance * flow_rate
    }

    /// Elastic work of breathing for one tidal breath \[J\].
    pub fn elastic_work(&self, tidal_vol: f64) -> f64 {
        0.5 * tidal_vol * tidal_vol / self.compliance
    }

    /// Resistive work of breathing for sinusoidal flow \[J\].
    pub fn resistive_work(&self, tidal_vol: f64, frequency: f64) -> f64 {
        let flow_peak = std::f64::consts::PI * frequency * tidal_vol;
        0.5 * self.resistance * flow_peak * flow_peak / frequency
    }
}

// ────────────────────────────────────────────────────────────────────────────
// SpineSegment – vertebral body + disc + ligaments
// ────────────────────────────────────────────────────────────────────────────

/// Parameters for a single spinal functional unit.
pub struct SpineSegmentParams {
    /// Disc height \[m\]
    pub disc_height: f64,
    /// Disc cross-section area \[m²\]
    pub disc_area: f64,
    /// Nucleus pulposus bulk modulus \[Pa\]
    pub k_nucleus: f64,
    /// Annulus fibrosus shear modulus \[Pa\]
    pub g_annulus: f64,
    /// Ligament stiffness (combined) \[N/m\]
    pub k_ligaments: f64,
    /// Vertebral body compressive stiffness \[N/m\]
    pub k_vertebra: f64,
}

impl Default for SpineSegmentParams {
    fn default() -> Self {
        Self {
            disc_height: 0.010,
            disc_area: 1.5e-3,
            k_nucleus: 2e6,
            g_annulus: 500e3,
            k_ligaments: 2000.0,
            k_vertebra: 80_000.0,
        }
    }
}

/// Lumped-parameter spinal functional unit (one motion segment).
pub struct SpineSegment {
    /// Compressive displacement \[m\]
    pub compression: f64,
    /// Shear displacement \[m\]
    pub shear: f64,
    /// Flexion angle \[rad\]
    pub flexion_angle: f64,
    /// Model parameters
    pub params: SpineSegmentParams,
}

impl SpineSegment {
    /// Create a `SpineSegment` at neutral position.
    pub fn new(params: SpineSegmentParams) -> Self {
        Self {
            compression: 0.0,
            shear: 0.0,
            flexion_angle: 0.0,
            params,
        }
    }

    /// Disc compressive stiffness \[N/m\].
    pub fn disc_compressive_stiffness(&self) -> f64 {
        self.params.k_nucleus * self.params.disc_area / self.params.disc_height
    }

    /// Intradiscal pressure from axial load \[Pa\].
    pub fn intradiscal_pressure(&self, axial_load: f64) -> f64 {
        axial_load / self.params.disc_area
    }

    /// Compressive force from disc \[N\].
    pub fn compressive_force(&self) -> f64 {
        self.disc_compressive_stiffness() * self.compression
    }

    /// Shear stiffness \[N/m\] from annulus.
    pub fn shear_stiffness(&self) -> f64 {
        self.params.g_annulus * self.params.disc_area / self.params.disc_height
    }

    /// Shear force from disc \[N\].
    pub fn shear_force(&self) -> f64 {
        self.shear_stiffness() * self.shear
    }

    /// Flexion bending stiffness \[N·m/rad\].
    pub fn flexion_stiffness(&self) -> f64 {
        let disc_bend =
            self.params.g_annulus * self.params.disc_area * self.params.disc_height / 12.0;
        disc_bend + self.params.k_ligaments * self.params.disc_height * 0.5
    }

    /// Bending moment from flexion \[N·m\].
    pub fn flexion_moment(&self) -> f64 {
        self.flexion_stiffness() * self.flexion_angle
    }

    /// Combined segment axial stiffness (disc + vertebrae in series) \[N/m\].
    pub fn total_axial_stiffness(&self) -> f64 {
        let k_disc = self.disc_compressive_stiffness();
        let k_vert = self.params.k_vertebra;
        1.0 / (1.0 / k_disc + 2.0 / k_vert)
    }
}

// ────────────────────────────────────────────────────────────────────────────
// Tests
// ────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── MuscleFiber tests ──────────────────────────────────────────────────

    #[test]
    fn test_hill_force_length_at_optimal() {
        let fiber = MuscleFiber::new(MuscleFiberParams::default());
        // At optimal length fl should be 1.0
        let fl = fiber.force_length();
        assert!(
            (fl - 1.0).abs() < 1e-9,
            "fl at optimal should be 1.0, got {:.6}",
            fl
        );
    }

    #[test]
    fn test_hill_force_length_off_optimal() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.length = 0.15; // 50% longer than l_opt = 0.10
        let fl = fiber.force_length();
        assert!(fl < 1.0, "fl off-optimal should be < 1.0, got {:.6}", fl);
        assert!(fl >= 0.0, "fl must be non-negative");
    }

    #[test]
    fn test_hill_force_velocity_isometric() {
        let fiber = MuscleFiber::new(MuscleFiberParams::default());
        // v=0 → fv should be 1.0
        let fv = fiber.force_velocity();
        assert!((fv - 1.0).abs() < 1e-9, "fv isometric = {:.6}", fv);
    }

    #[test]
    fn test_hill_force_velocity_max_shortening() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        // v = -v_max → force should be ~0
        fiber.velocity = -fiber.params.v_max;
        let fv = fiber.force_velocity();
        assert!(fv.abs() < 1e-9, "fv at max shortening = {:.6}", fv);
    }

    #[test]
    fn test_hill_force_velocity_lengthening() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.velocity = 0.1; // lengthening
        let fv = fiber.force_velocity();
        assert!(fv > 1.0, "eccentric fv should exceed 1.0, got {:.6}", fv);
    }

    #[test]
    fn test_passive_force_zero_at_optimal() {
        let fiber = MuscleFiber::new(MuscleFiberParams::default());
        assert_eq!(fiber.passive_force(), 0.0);
    }

    #[test]
    fn test_passive_force_positive_beyond_optimal() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.length = 0.12; // beyond l_opt
        assert!(fiber.passive_force() > 0.0);
    }

    #[test]
    fn test_total_force_zero_activation() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.activation = 0.0;
        assert_eq!(fiber.total_force(), 0.0);
    }

    #[test]
    fn test_activation_dynamics_rise() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.excitation = 1.0;
        fiber.activation = 0.0;
        fiber.update_activation(0.01);
        assert!(fiber.activation > 0.0, "activation should rise");
        assert!(fiber.activation <= 1.0);
    }

    #[test]
    fn test_activation_dynamics_decay_to_zero() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.excitation = 0.0;
        fiber.activation = 1.0;
        for _ in 0..500 {
            fiber.update_activation(0.001);
        }
        assert!(
            fiber.activation < 0.01,
            "activation should decay near 0, got {:.6}",
            fiber.activation
        );
    }

    #[test]
    fn test_activation_clamped() {
        let mut fiber = MuscleFiber::new(MuscleFiberParams::default());
        fiber.excitation = 2.0; // out of range
        fiber.activation = 0.0;
        fiber.update_activation(10.0);
        assert!(fiber.activation <= 1.0);
    }

    // ── MusculoskeletalModel tests ─────────────────────────────────────────

    #[test]
    fn test_musculoskeletal_joint_torque() {
        let muscles = vec![
            MuscleEntry {
                f_max: 1000.0,
                moment_arm: 0.05,
                normalized_force: 0.5,
            },
            MuscleEntry {
                f_max: 800.0,
                moment_arm: 0.04,
                normalized_force: 0.6,
            },
        ];
        let model = MusculoskeletalModel::new(muscles, 40.0, 2.0);
        let torque = model.joint_torque();
        // 0.5*1000*0.05 + 0.6*800*0.04 = 25 + 19.2 = 44.2
        assert!((torque - 44.2).abs() < 1e-9);
    }

    #[test]
    fn test_musculoskeletal_static_solve() {
        let muscles = vec![
            MuscleEntry {
                f_max: 1000.0,
                moment_arm: 0.05,
                normalized_force: 0.0,
            },
            MuscleEntry {
                f_max: 500.0,
                moment_arm: 0.03,
                normalized_force: 0.0,
            },
        ];
        let mut model = MusculoskeletalModel::new(muscles, 30.0, 2.0);
        model.solve_static();
        // forces should be distributed (non-negative)
        for m in &model.muscles {
            assert!(m.normalized_force >= 0.0);
            assert!(m.normalized_force <= 1.0);
        }
    }

    #[test]
    fn test_musculoskeletal_objective_positive() {
        let muscles = vec![MuscleEntry {
            f_max: 1000.0,
            moment_arm: 0.05,
            normalized_force: 0.3,
        }];
        let model = MusculoskeletalModel::new(muscles, 15.0, 2.0);
        assert!(model.objective() > 0.0);
    }

    // ── CardiacFem tests ───────────────────────────────────────────────────

    #[test]
    fn test_cardiac_no_calcium_no_active_stress() {
        let mut cardiac = CardiacFem::new(CardiacParams::default());
        cardiac.calcium = 0.0;
        cardiac.update_stress();
        assert_eq!(cardiac.active_stress, 0.0);
    }

    #[test]
    fn test_cardiac_stress_increases_with_calcium() {
        let mut c1 = CardiacFem::new(CardiacParams::default());
        c1.calcium = 0.5;
        c1.update_stress();
        let mut c2 = CardiacFem::new(CardiacParams::default());
        c2.calcium = 2.0;
        c2.update_stress();
        assert!(
            c2.active_stress > c1.active_stress,
            "higher Ca → higher active stress"
        );
    }

    #[test]
    fn test_cardiac_calcium_activation_half_max() {
        let mut cardiac = CardiacFem::new(CardiacParams::default());
        cardiac.calcium = cardiac.params.ca_50;
        let act = cardiac.calcium_activation();
        assert!(
            (act - 0.5).abs() < 0.01,
            "at Ca50 activation ≈ 0.5, got {:.6}",
            act
        );
    }

    #[test]
    fn test_cardiac_passive_stress_at_reference() {
        let mut cardiac = CardiacFem::new(CardiacParams::default());
        cardiac.lambda_fiber = 1.0;
        cardiac.update_stress();
        assert_eq!(cardiac.passive_stress, 0.0);
    }

    #[test]
    fn test_cardiac_passive_stress_positive_stretch() {
        let mut cardiac = CardiacFem::new(CardiacParams::default());
        cardiac.lambda_fiber = 1.2;
        cardiac.update_stress();
        assert!(cardiac.passive_stress > 0.0);
    }

    #[test]
    fn test_cardiac_total_stress_sum() {
        let mut cardiac = CardiacFem::new(CardiacParams::default());
        cardiac.calcium = 1.0;
        cardiac.lambda_fiber = 1.1;
        cardiac.update_stress();
        let expected = cardiac.active_stress + cardiac.passive_stress;
        assert!((cardiac.total_stress() - expected).abs() < 1e-9);
    }

    // ── BoneRemodeling tests ───────────────────────────────────────────────

    #[test]
    fn test_bone_density_increases_high_strain() {
        let mut bone = BoneRemodeling::new(500.0, BoneRemodelingParams::default());
        bone.stress_stimulus = 5000.0; // above sigma_ref
        bone.update_density(1e6);
        assert!(
            bone.density > 500.0,
            "density should increase under high strain"
        );
    }

    #[test]
    fn test_bone_density_decreases_low_strain() {
        let mut bone = BoneRemodeling::new(800.0, BoneRemodelingParams::default());
        bone.stress_stimulus = 100.0; // below sigma_ref
        bone.update_density(1e6);
        assert!(bone.density < 800.0, "density should decrease under disuse");
    }

    #[test]
    fn test_bone_density_bounds() {
        let mut bone = BoneRemodeling::new(100.0, BoneRemodelingParams::default());
        bone.stress_stimulus = 0.0;
        bone.update_density(1e12); // very long time
        assert!(bone.density >= bone.params.rho_min);
        assert!(bone.density <= bone.params.rho_max);
    }

    #[test]
    fn test_bone_no_change_in_dead_zone() {
        let params = BoneRemodelingParams::default();
        let sigma_ref = params.sigma_ref;
        let mut bone = BoneRemodeling::new(500.0, params);
        bone.stress_stimulus = sigma_ref; // exactly at reference → in dead zone
        let rho_before = bone.density;
        bone.update_density(1e6);
        assert!((bone.density - rho_before).abs() < 1e-9);
    }

    #[test]
    fn test_bone_elastic_modulus_positive() {
        let bone = BoneRemodeling::new(500.0, BoneRemodelingParams::default());
        assert!(bone.elastic_modulus() > 0.0);
    }

    // ── LigamentModel tests ────────────────────────────────────────────────

    #[test]
    fn test_ligament_no_force_in_compression() {
        let mut lig = LigamentModel::new(1e-4, 0.03, 100e6, 3.0);
        lig.strain = -0.05;
        assert_eq!(
            lig.force(),
            0.0,
            "ligament should not carry compressive force"
        );
    }

    #[test]
    fn test_ligament_force_positive_in_tension() {
        let mut lig = LigamentModel::new(1e-4, 0.03, 100e6, 3.0);
        lig.strain = 0.05;
        assert!(lig.force() > 0.0);
    }

    #[test]
    fn test_ligament_linear_region_stiffer() {
        let mut lig = LigamentModel::new(1e-4, 0.03, 100e6, 3.0);
        lig.strain = 0.015; // toe region
        let f_toe = lig.stress();
        lig.strain = 0.05; // linear region
        let f_lin = lig.stress();
        assert!(f_lin > f_toe, "linear region should have higher stress");
    }

    #[test]
    fn test_ligament_zero_strain_zero_force() {
        let mut lig = LigamentModel::new(1e-4, 0.03, 100e6, 3.0);
        lig.strain = 0.0;
        assert_eq!(lig.force(), 0.0);
    }

    // ── ArteryFem tests ────────────────────────────────────────────────────

    #[test]
    fn test_artery_laplace_hoop_stress() {
        // r=5mm, h=0.5mm, p=13.3kPa (100 mmHg)
        let artery = ArteryFem::new(0.005, 0.0055, 13_300.0, 500e3, 0.45);
        let sigma = artery.laplace_hoop_stress();
        let expected = 13_300.0 * 0.005 / 0.0005;
        assert!(
            (sigma - expected).abs() < 1.0,
            "Laplace stress {:.6} vs {:.6}",
            sigma,
            expected
        );
    }

    #[test]
    fn test_artery_hoop_stress_scales_with_pressure() {
        let a1 = ArteryFem::new(0.005, 0.0055, 10_000.0, 500e3, 0.45);
        let a2 = ArteryFem::new(0.005, 0.0055, 20_000.0, 500e3, 0.45);
        assert!((a2.laplace_hoop_stress() / a1.laplace_hoop_stress() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn test_artery_wall_thickness() {
        let artery = ArteryFem::new(0.005, 0.0055, 13_300.0, 500e3, 0.45);
        assert!((artery.wall_thickness() - 0.0005).abs() < 1e-12);
    }

    #[test]
    fn test_artery_lame_hoop_stress_larger_than_laplace() {
        // Lamé inner stress > Laplace thin-wall approx for thick walls
        let artery = ArteryFem::new(0.005, 0.010, 13_300.0, 500e3, 0.45);
        let laplace = artery.laplace_hoop_stress();
        let lame = artery.lame_hoop_stress_inner();
        // Lamé includes pressure magnification factor
        assert!(
            lame > 0.0,
            "Lamé stress should be positive, got {:.6}",
            lame
        );
        let _ = laplace; // used for context
    }

    // ── BreathingMechanics tests ───────────────────────────────────────────

    #[test]
    fn test_breathing_tidal_volume() {
        // compliance = 2e-4 m³/Pa, ΔP = 250 Pa → VT = 2e-4 * 250 = 0.05 m³
        let bm = BreathingMechanics::new(2e-4, 200.0, 2.5e-3);
        let vt = bm.tidal_volume(250.0);
        let expected = 2e-4_f64 * 250.0;
        assert!((vt - expected).abs() < 1e-12, "VT = {:.6}", vt);
    }

    #[test]
    fn test_breathing_recoil_at_frc() {
        let bm = BreathingMechanics::new(2e-4, 200.0, 2.5e-3);
        assert_eq!(bm.recoil_pressure(), 0.0);
    }

    #[test]
    fn test_breathing_elastic_work_positive() {
        let bm = BreathingMechanics::new(2e-4, 200.0, 2.5e-3);
        assert!(bm.elastic_work(0.5e-3) > 0.0);
    }

    #[test]
    fn test_breathing_resistive_pressure() {
        let bm = BreathingMechanics::new(2e-4, 200.0, 2.5e-3);
        let dp = bm.resistive_pressure(1e-4);
        assert!((dp - 200.0 * 1e-4).abs() < 1e-12);
    }

    // ── SpineSegment tests ─────────────────────────────────────────────────

    #[test]
    fn test_spine_compressive_stiffness_positive() {
        let seg = SpineSegment::new(SpineSegmentParams::default());
        assert!(seg.disc_compressive_stiffness() > 0.0);
    }

    #[test]
    fn test_spine_intradiscal_pressure() {
        let seg = SpineSegment::new(SpineSegmentParams::default());
        let p = seg.intradiscal_pressure(1000.0);
        let expected = 1000.0 / SpineSegmentParams::default().disc_area;
        assert!((p - expected).abs() < 1.0);
    }

    #[test]
    fn test_spine_compressive_force() {
        let mut seg = SpineSegment::new(SpineSegmentParams::default());
        seg.compression = 1e-4;
        assert!(seg.compressive_force() > 0.0);
    }

    #[test]
    fn test_spine_shear_force() {
        let mut seg = SpineSegment::new(SpineSegmentParams::default());
        seg.shear = 1e-3;
        assert!(seg.shear_force() > 0.0);
    }

    #[test]
    fn test_spine_flexion_moment() {
        let mut seg = SpineSegment::new(SpineSegmentParams::default());
        seg.flexion_angle = 0.1;
        assert!(seg.flexion_moment() > 0.0);
    }

    #[test]
    fn test_spine_total_axial_stiffness_less_than_disc() {
        let seg = SpineSegment::new(SpineSegmentParams::default());
        let k_disc = seg.disc_compressive_stiffness();
        let k_total = seg.total_axial_stiffness();
        assert!(
            k_total < k_disc,
            "total stiffness (series) must be less than disc alone"
        );
    }
}
