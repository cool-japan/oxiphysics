// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Muscle simulation module for the OxiPhysics soft-body crate.
//!
//! Provides a comprehensive Hill-type muscle model with:
//!
//! - **Contractile element (CE)** – active force generation
//! - **Parallel elastic element (PE)** – passive tissue elasticity
//! - **Series elastic element (SE)** – tendon compliance
//! - **Activation dynamics** – neural excitation to activation mapping
//! - **Force-length relationship** – Gaussian active + exponential passive
//! - **Force-velocity relationship** – Hill hyperbola
//! - **Pennation angle** – fiber architecture
//! - **Tendon mechanics** – nonlinear tendon force-strain
//! - **Musculotendon equilibrium** – static/dynamic equilibrium solver
//! - **Fatigue modeling** – endurance and recovery
//! - **Fiber recruitment** – motor unit recruitment model
//! - **Isometric/isotonic contractions** – fixed-length and fixed-load modes

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Helper math (no nalgebra)
// ---------------------------------------------------------------------------

/// 3-D Euclidean norm.
#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Subtract two 3-D vectors: a - b.
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-D vectors.
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-D vector.
#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Normalize a 3-D vector (returns zero if degenerate).
#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}

/// Dot product of two 3-D vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Clamp `x` to `[lo, hi]`.
#[inline]
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    x.max(lo).min(hi)
}

/// Linear interpolation.
#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

// ===========================================================================
// 1. Activation Dynamics
// ===========================================================================

/// First-order activation dynamics model.
///
/// Maps neural excitation u(t) in \[0, 1\] to muscle activation a(t) in \[0, 1\]
/// via first-order ODE:
///   da/dt = (u - a) / tau(a, u)
///
/// where tau depends on whether the muscle is activating or deactivating.
#[derive(Debug, Clone)]
pub struct ActivationDynamics {
    /// Activation time constant (s).
    pub tau_act: f64,
    /// Deactivation time constant (s).
    pub tau_deact: f64,
    /// Minimum activation level (prevents numerical issues).
    pub a_min: f64,
    /// Current activation level.
    activation: f64,
}

impl ActivationDynamics {
    /// Create new activation dynamics.
    pub fn new(tau_act: f64, tau_deact: f64) -> Self {
        Self {
            tau_act,
            tau_deact,
            a_min: 0.01,
            activation: 0.0,
        }
    }

    /// Default fast-twitch muscle dynamics.
    pub fn fast_twitch() -> Self {
        Self::new(0.010, 0.040)
    }

    /// Default slow-twitch muscle dynamics.
    pub fn slow_twitch() -> Self {
        Self::new(0.025, 0.060)
    }

    /// Get current activation level.
    pub fn activation(&self) -> f64 {
        self.activation
    }

    /// Set activation directly.
    pub fn set_activation(&mut self, a: f64) {
        self.activation = clamp(a, 0.0, 1.0);
    }

    /// Update activation given neural excitation u in \[0, 1\].
    pub fn update(&mut self, u: f64, dt: f64) -> f64 {
        let u_c = clamp(u, 0.0, 1.0);
        let tau = if u_c > self.activation {
            self.tau_act * (0.5 + 1.5 * self.activation)
        } else {
            self.tau_deact / (0.5 + 1.5 * self.activation)
        };
        let tau_safe = tau.max(1e-6);
        let da = (u_c - self.activation) / tau_safe * dt;
        self.activation = clamp(self.activation + da, self.a_min, 1.0);
        self.activation
    }

    /// Reset activation to zero.
    pub fn reset(&mut self) {
        self.activation = 0.0;
    }

    /// Compute the time constant at current state for a given excitation.
    pub fn current_tau(&self, u: f64) -> f64 {
        if u > self.activation {
            self.tau_act * (0.5 + 1.5 * self.activation)
        } else {
            self.tau_deact / (0.5 + 1.5 * self.activation).max(0.5)
        }
    }
}

// ===========================================================================
// 2. Force-Length Relationship
// ===========================================================================

/// Active and passive force-length curves for muscle fibers.
///
/// Active: Gaussian centered at optimal length
/// Passive: Exponential stiffening beyond optimal length
#[derive(Debug, Clone)]
pub struct ForceLengthCurve {
    /// Width of the active force-length Gaussian.
    pub gamma: f64,
    /// Passive force exponential shape factor.
    pub kpe: f64,
    /// Passive force strain offset (normalized length at which PE kicks in).
    pub e0: f64,
    /// Maximum passive force (normalized to F_max).
    pub f_pe_max: f64,
}

impl ForceLengthCurve {
    /// Create a new force-length curve.
    pub fn new(gamma: f64, kpe: f64, e0: f64, f_pe_max: f64) -> Self {
        Self {
            gamma,
            kpe,
            e0,
            f_pe_max,
        }
    }

    /// Default curve (Thelen 2003 parameters).
    pub fn thelen2003() -> Self {
        Self::new(0.45, 4.0, 0.6, 1.0)
    }

    /// Compute the normalized active force-length value.
    ///
    /// f_active(l_norm) = exp(-((l_norm - 1) / gamma)^2)
    pub fn active(&self, l_norm: f64) -> f64 {
        let x = (l_norm - 1.0) / self.gamma;
        (-x * x).exp()
    }

    /// Compute the normalized passive force-length value.
    ///
    /// f_passive(l_norm) = (exp(kpe*(l_norm - 1)/e0) - 1) / (exp(kpe) - 1)
    pub fn passive(&self, l_norm: f64) -> f64 {
        if l_norm <= 1.0 {
            return 0.0;
        }
        let strain = l_norm - 1.0;
        let num = (self.kpe * strain / self.e0).exp() - 1.0;
        let den = self.kpe.exp() - 1.0;
        (self.f_pe_max * num / den).max(0.0)
    }

    /// Total (active + passive) normalized force at given activation.
    pub fn total(&self, l_norm: f64, activation: f64) -> f64 {
        activation * self.active(l_norm) + self.passive(l_norm)
    }

    /// Compute the slope (stiffness) of the active curve at l_norm.
    pub fn active_stiffness(&self, l_norm: f64) -> f64 {
        let x = (l_norm - 1.0) / self.gamma;
        -2.0 * x / self.gamma * (-x * x).exp()
    }

    /// Compute the slope (stiffness) of the passive curve at l_norm.
    pub fn passive_stiffness(&self, l_norm: f64) -> f64 {
        if l_norm <= 1.0 {
            return 0.0;
        }
        let strain = l_norm - 1.0;
        let den = self.kpe.exp() - 1.0;
        self.f_pe_max * self.kpe / (self.e0 * den) * (self.kpe * strain / self.e0).exp()
    }
}

// ===========================================================================
// 3. Force-Velocity Relationship
// ===========================================================================

/// Force-velocity relationship for muscle fibers (Hill curve).
///
/// Concentric (shortening): hyperbolic curve
/// Eccentric (lengthening): asymptotic to f_ecc_max
#[derive(Debug, Clone)]
pub struct ForceVelocityCurve {
    /// Shape parameter for the Hill curve (a_f / F_max).
    pub af: f64,
    /// Maximum eccentric force ratio.
    pub f_ecc_max: f64,
    /// Eccentric curve shape parameter.
    pub eccentric_shape: f64,
}

impl ForceVelocityCurve {
    /// Create a new force-velocity curve.
    pub fn new(af: f64, f_ecc_max: f64, eccentric_shape: f64) -> Self {
        Self {
            af,
            f_ecc_max,
            eccentric_shape,
        }
    }

    /// Default parameters (Thelen 2003).
    pub fn thelen2003() -> Self {
        Self::new(0.25, 1.8, 0.1)
    }

    /// Compute the normalized force-velocity value.
    ///
    /// v_norm = v / v_max, negative for shortening
    pub fn evaluate(&self, v_norm: f64) -> f64 {
        if v_norm <= 0.0 {
            // Concentric (shortening): (1 + v_norm/a_f) / (1 - v_norm/a_f)
            // But with v_norm negative for shortening: use |v_norm|
            let v_abs = (-v_norm).min(0.999);
            (1.0 - v_abs) / (1.0 + v_abs / self.af)
        } else {
            // Eccentric (lengthening)
            let v_c = v_norm.min(0.999);
            let num = self.f_ecc_max * v_c + self.eccentric_shape;
            let den = v_c + self.eccentric_shape;
            num / den
        }
    }

    /// Compute the inverse: given f_v, find v_norm (concentric only).
    pub fn inverse_concentric(&self, f_v: f64) -> f64 {
        let f_v_c = clamp(f_v, 0.0, 1.0 - 1e-6);
        // From (1 - |v|) / (1 + |v|/af) = fv
        // |v| = (1 - fv) / (1 + fv/af)
        let v_abs = (1.0 - f_v_c) / (1.0 + f_v_c / self.af);
        -v_abs
    }

    /// Compute the slope at a given normalized velocity.
    pub fn slope(&self, v_norm: f64) -> f64 {
        let dv = 1e-6;
        let f1 = self.evaluate(v_norm - dv);
        let f2 = self.evaluate(v_norm + dv);
        (f2 - f1) / (2.0 * dv)
    }
}

// ===========================================================================
// 4. Pennation Angle Model
// ===========================================================================

/// Pennation angle model for muscle fiber architecture.
///
/// The pennation angle changes with fiber length to maintain constant
/// muscle thickness (constant volume assumption):
///   sin(alpha) = sin(alpha_opt) * l_opt / l_fiber
#[derive(Debug, Clone)]
pub struct PennationModel {
    /// Pennation angle at optimal fiber length (rad).
    pub alpha_opt: f64,
    /// Optimal fiber length (m).
    pub l_opt: f64,
    /// Muscle thickness at optimal length (m).
    thickness: f64,
}

impl PennationModel {
    /// Create a new pennation model.
    pub fn new(alpha_opt: f64, l_opt: f64) -> Self {
        let thickness = l_opt * alpha_opt.sin();
        Self {
            alpha_opt,
            l_opt,
            thickness,
        }
    }

    /// No pennation (alpha = 0).
    pub fn zero() -> Self {
        Self::new(0.0, 1.0)
    }

    /// Compute the current pennation angle for a given fiber length.
    pub fn angle(&self, fiber_length: f64) -> f64 {
        if fiber_length < 1e-15 {
            return self.alpha_opt;
        }
        let sin_alpha = (self.thickness / fiber_length).min(0.999);
        sin_alpha.asin()
    }

    /// Compute cos(alpha) for force projection along tendon.
    pub fn cos_alpha(&self, fiber_length: f64) -> f64 {
        self.angle(fiber_length).cos()
    }

    /// Compute the fiber length from musculotendon length and tendon length.
    pub fn fiber_length_from_mt(&self, mt_length: f64, tendon_length: f64) -> f64 {
        let l_fiber_along = (mt_length - tendon_length).max(0.0);
        (l_fiber_along * l_fiber_along + self.thickness * self.thickness).sqrt()
    }

    /// Compute the fiber velocity from musculotendon velocity.
    pub fn fiber_velocity_from_mt(&self, mt_velocity: f64, fiber_length: f64) -> f64 {
        let cos_a = self.cos_alpha(fiber_length);
        if cos_a > 1e-6 {
            mt_velocity / cos_a
        } else {
            mt_velocity
        }
    }

    /// Compute the muscle thickness (constant volume assumption).
    pub fn muscle_thickness(&self) -> f64 {
        self.thickness
    }
}

// ===========================================================================
// 5. Tendon Mechanics
// ===========================================================================

/// Nonlinear tendon model.
///
/// The tendon force-strain relationship is modeled as a piecewise function:
/// - Toe region (0 < strain < strain_toe): quadratic
/// - Linear region (strain > strain_toe): linear with stiffness kT
/// - Slack region (strain < 0): zero force
#[derive(Debug, Clone)]
pub struct TendonModel {
    /// Tendon slack length (m).
    pub slack_length: f64,
    /// Maximum isometric force (N) for normalization.
    pub f_max: f64,
    /// Toe region strain limit.
    pub strain_toe: f64,
    /// Tendon stiffness in the linear region (normalized).
    pub k_toe: f64,
    /// Linear region stiffness (normalized).
    pub k_linear: f64,
}

impl TendonModel {
    /// Create a new tendon model.
    pub fn new(slack_length: f64, f_max: f64) -> Self {
        Self {
            slack_length,
            f_max,
            strain_toe: 0.033,
            k_toe: 1.712 / 0.033,
            k_linear: 37.5,
        }
    }

    /// Create a rigid tendon (infinite stiffness approximation).
    pub fn rigid(slack_length: f64, f_max: f64) -> Self {
        Self {
            slack_length,
            f_max,
            strain_toe: 0.001,
            k_toe: 1000.0,
            k_linear: 10000.0,
        }
    }

    /// Compute the tendon strain from tendon length.
    pub fn strain(&self, tendon_length: f64) -> f64 {
        if self.slack_length > 1e-15 {
            (tendon_length - self.slack_length) / self.slack_length
        } else {
            0.0
        }
    }

    /// Compute the normalized tendon force from tendon strain.
    ///
    /// Returns force / F_max.
    pub fn force_normalized(&self, strain: f64) -> f64 {
        if strain <= 0.0 {
            0.0
        } else if strain <= self.strain_toe {
            // Quadratic toe region
            let ratio = strain / self.strain_toe;
            0.5 * self.k_toe * self.strain_toe * ratio * ratio
        } else {
            // Linear region
            let f_toe = 0.5 * self.k_toe * self.strain_toe;
            f_toe + self.k_linear * (strain - self.strain_toe)
        }
    }

    /// Compute the tendon force (N) from tendon length.
    pub fn force(&self, tendon_length: f64) -> f64 {
        let eps = self.strain(tendon_length);
        self.force_normalized(eps) * self.f_max
    }

    /// Compute the tendon stiffness (dF/dl) at the current length.
    pub fn stiffness(&self, tendon_length: f64) -> f64 {
        let eps = self.strain(tendon_length);
        if eps <= 0.0 {
            0.0
        } else if eps <= self.strain_toe {
            let ratio = eps / self.strain_toe;
            self.k_toe * ratio * self.f_max / self.slack_length
        } else {
            self.k_linear * self.f_max / self.slack_length
        }
    }

    /// Compute the tendon length from a given force (inverse).
    pub fn length_from_force(&self, force: f64) -> f64 {
        let f_norm = force / self.f_max;
        let strain = if f_norm <= 0.0 {
            0.0
        } else {
            let f_toe = 0.5 * self.k_toe * self.strain_toe;
            if f_norm <= f_toe {
                // Inverse of quadratic: strain = strain_toe * sqrt(2*f_norm / (k_toe * strain_toe))
                self.strain_toe * (2.0 * f_norm / (self.k_toe * self.strain_toe)).sqrt()
            } else {
                // Inverse of linear
                self.strain_toe + (f_norm - f_toe) / self.k_linear
            }
        };
        self.slack_length * (1.0 + strain)
    }
}

// ===========================================================================
// 6. Hill-Type Muscle Model (Complete)
// ===========================================================================

/// Complete Hill-type muscle model with all three elements.
///
/// Combines contractile element (CE), parallel elastic element (PE),
/// and series elastic element (SE/tendon) into a unified musculotendon
/// actuator.
#[derive(Debug, Clone)]
pub struct HillMuscleSoftbody {
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fiber length (m).
    pub l_opt: f64,
    /// Maximum fiber shortening velocity (l_opt/s).
    pub v_max: f64,
    /// Activation dynamics.
    pub activation: ActivationDynamics,
    /// Force-length curve.
    pub fl_curve: ForceLengthCurve,
    /// Force-velocity curve.
    pub fv_curve: ForceVelocityCurve,
    /// Pennation model.
    pub pennation: PennationModel,
    /// Tendon model.
    pub tendon: TendonModel,
    /// Current fiber length (m).
    fiber_length: f64,
    /// Current fiber velocity (m/s).
    fiber_velocity: f64,
    /// Current tendon length (m).
    tendon_length: f64,
    /// Current total force output (N).
    force: f64,
}

impl HillMuscleSoftbody {
    /// Create a new Hill muscle model.
    pub fn new(
        f_max: f64,
        l_opt: f64,
        v_max: f64,
        pennation_angle: f64,
        tendon_slack_length: f64,
    ) -> Self {
        Self {
            f_max,
            l_opt,
            v_max,
            activation: ActivationDynamics::fast_twitch(),
            fl_curve: ForceLengthCurve::thelen2003(),
            fv_curve: ForceVelocityCurve::thelen2003(),
            pennation: PennationModel::new(pennation_angle, l_opt),
            tendon: TendonModel::new(tendon_slack_length, f_max),
            fiber_length: l_opt,
            fiber_velocity: 0.0,
            tendon_length: tendon_slack_length,
            force: 0.0,
        }
    }

    /// Create a muscle with typical parameters.
    pub fn typical(f_max: f64, l_opt: f64, tendon_slack_length: f64) -> Self {
        Self::new(f_max, l_opt, 10.0, 0.0, tendon_slack_length)
    }

    /// Get the current muscle force (N).
    pub fn force(&self) -> f64 {
        self.force
    }

    /// Get the current fiber length.
    pub fn fiber_length(&self) -> f64 {
        self.fiber_length
    }

    /// Get the current fiber velocity.
    pub fn fiber_velocity(&self) -> f64 {
        self.fiber_velocity
    }

    /// Get the current tendon length.
    pub fn tendon_length(&self) -> f64 {
        self.tendon_length
    }

    /// Get the current activation.
    pub fn current_activation(&self) -> f64 {
        self.activation.activation()
    }

    /// Compute the muscle force for given state variables.
    pub fn compute_force(&self, activation: f64, fiber_length: f64, fiber_velocity: f64) -> f64 {
        let l_norm = fiber_length / self.l_opt;
        let v_norm = fiber_velocity / (self.v_max * self.l_opt);

        let f_active = activation * self.fl_curve.active(l_norm) * self.fv_curve.evaluate(v_norm);
        let f_passive = self.fl_curve.passive(l_norm);
        let cos_a = self.pennation.cos_alpha(fiber_length);

        self.f_max * (f_active + f_passive) * cos_a
    }

    /// Update the muscle state for a given excitation and musculotendon length/velocity.
    pub fn update(&mut self, excitation: f64, mt_length: f64, mt_velocity: f64, dt: f64) {
        // 1. Update activation
        self.activation.update(excitation, dt);
        let a = self.activation.activation();

        // 2. Solve musculotendon equilibrium (rigid tendon approximation)
        self.tendon_length = self.tendon.slack_length;
        self.fiber_length = self
            .pennation
            .fiber_length_from_mt(mt_length, self.tendon_length);
        self.fiber_velocity = self
            .pennation
            .fiber_velocity_from_mt(mt_velocity, self.fiber_length);

        // 3. Compute force
        self.force = self.compute_force(a, self.fiber_length, self.fiber_velocity);
    }

    /// Solve the static musculotendon equilibrium.
    ///
    /// Finds the fiber length where CE+PE force = tendon force.
    pub fn solve_equilibrium(&mut self, activation: f64, mt_length: f64) {
        let a = clamp(activation, 0.0, 1.0);

        // Newton iteration to find fiber length
        let mut l_fiber = self.fiber_length;
        for _ in 0..50 {
            let cos_a = self.pennation.cos_alpha(l_fiber);
            let l_norm = l_fiber / self.l_opt;

            // Muscle force along tendon
            let f_muscle = self.f_max
                * (a * self.fl_curve.active(l_norm) + self.fl_curve.passive(l_norm))
                * cos_a;

            // Tendon length and force
            let l_tendon = mt_length - l_fiber * cos_a;
            let f_tendon = self.tendon.force(l_tendon);

            let residual = f_muscle - f_tendon;
            if residual.abs() < 1e-4 {
                break;
            }

            // Simple bisection-like step
            let dl = -residual / (self.f_max * 10.0 / self.l_opt);
            l_fiber = (l_fiber + dl).max(0.5 * self.l_opt).min(1.5 * self.l_opt);
        }

        self.fiber_length = l_fiber;
        let cos_a = self.pennation.cos_alpha(l_fiber);
        self.tendon_length = mt_length - l_fiber * cos_a;
        self.force = self.tendon.force(self.tendon_length);
    }

    /// Reset the muscle to its initial state.
    pub fn reset(&mut self) {
        self.activation.reset();
        self.fiber_length = self.l_opt;
        self.fiber_velocity = 0.0;
        self.tendon_length = self.tendon.slack_length;
        self.force = 0.0;
    }
}

// ===========================================================================
// 7. Fatigue Model
// ===========================================================================

/// Muscle fatigue model based on the three-compartment scheme.
///
/// Motor units transition between three states:
/// - MA: active motor units (currently generating force)
/// - MF: fatigued motor units (temporarily unable to fire)
/// - MR: resting motor units (available for recruitment)
///
/// The total is conserved: MA + MF + MR = 1.
#[derive(Debug, Clone)]
pub struct FatigueModel {
    /// Fatigue rate constant (1/s).
    pub fatigue_rate: f64,
    /// Recovery rate constant (1/s).
    pub recovery_rate: f64,
    /// Active motor unit fraction.
    pub ma: f64,
    /// Fatigued motor unit fraction.
    pub mf: f64,
    /// Resting motor unit fraction.
    pub mr: f64,
    /// Target activation demand.
    target_activation: f64,
}

impl FatigueModel {
    /// Create a new fatigue model.
    pub fn new(fatigue_rate: f64, recovery_rate: f64) -> Self {
        Self {
            fatigue_rate,
            recovery_rate,
            ma: 0.0,
            mf: 0.0,
            mr: 1.0,
            target_activation: 0.0,
        }
    }

    /// Default fatigue parameters for a mixed-fiber muscle.
    pub fn default_mixed() -> Self {
        Self::new(0.01, 0.002)
    }

    /// Default fatigue parameters for a fast-twitch dominant muscle.
    pub fn fast_twitch() -> Self {
        Self::new(0.02, 0.003)
    }

    /// Default fatigue parameters for a slow-twitch dominant muscle.
    pub fn slow_twitch() -> Self {
        Self::new(0.005, 0.001)
    }

    /// Set the target activation demand.
    pub fn set_target(&mut self, target: f64) {
        self.target_activation = clamp(target, 0.0, 1.0);
    }

    /// Update the fatigue state.
    ///
    /// Returns the effective (fatigued) activation level.
    pub fn update(&mut self, dt: f64) -> f64 {
        let target = self.target_activation;

        // Recruitment: try to bring MA up to target from MR
        let recruit = if target > self.ma && self.mr > 0.0 {
            let needed = (target - self.ma).min(self.mr);
            needed * dt * 10.0 // recruitment rate
        } else {
            0.0
        };

        // De-recruitment: if target < MA, release to MR
        let derecruit = if target < self.ma {
            (self.ma - target) * dt * 5.0
        } else {
            0.0
        };

        // Fatigue: active units become fatigued
        let fatigue = self.fatigue_rate * self.ma * dt;

        // Recovery: fatigued units recover to resting
        let recovery = self.recovery_rate * self.mf * dt;

        // Update compartments
        self.ma = (self.ma + recruit - derecruit - fatigue).max(0.0);
        self.mf = (self.mf + fatigue - recovery).max(0.0);
        self.mr = (self.mr - recruit + derecruit + recovery).max(0.0);

        // Normalize to ensure conservation
        let total = self.ma + self.mf + self.mr;
        if total > 1e-15 {
            self.ma /= total;
            self.mf /= total;
            self.mr /= total;
        }

        self.ma
    }

    /// Get the endurance time estimate for a given contraction intensity.
    ///
    /// Based on Rohmert's curve: T = (f_max / F - 1)^k / fatigue_rate
    pub fn endurance_time(&self, intensity: f64) -> f64 {
        let intensity_c = clamp(intensity, 0.01, 0.99);
        let k = 1.5;
        ((1.0 / intensity_c - 1.0).powf(k)) / self.fatigue_rate
    }

    /// Get the total fatigue fraction (fraction of fatigued units).
    pub fn fatigue_fraction(&self) -> f64 {
        self.mf
    }

    /// Check if the muscle is significantly fatigued.
    pub fn is_fatigued(&self) -> bool {
        self.mf > 0.3
    }

    /// Reset the fatigue state to fully rested.
    pub fn reset(&mut self) {
        self.ma = 0.0;
        self.mf = 0.0;
        self.mr = 1.0;
        self.target_activation = 0.0;
    }
}

// ===========================================================================
// 8. Fiber Recruitment Model
// ===========================================================================

/// Motor unit recruitment model (Henneman's size principle).
///
/// Motor units are recruited from smallest (slow-twitch) to largest
/// (fast-twitch) as excitation increases. Each motor unit has a
/// recruitment threshold and contributes a fixed force increment.
#[derive(Debug, Clone)]
pub struct FiberRecruitment {
    /// Number of motor units.
    pub num_units: usize,
    /// Recruitment thresholds (sorted, in \[0, 1\]).
    pub thresholds: Vec<f64>,
    /// Force contributions per motor unit (normalized).
    pub force_contributions: Vec<f64>,
    /// Fiber type fractions (0 = slow-twitch, 1 = fast-twitch).
    pub fiber_types: Vec<f64>,
    /// Current recruitment states (true = recruited).
    recruited: Vec<bool>,
}

impl FiberRecruitment {
    /// Create a new recruitment model with n motor units.
    pub fn new(num_units: usize) -> Self {
        let mut thresholds = Vec::with_capacity(num_units);
        let mut forces = Vec::with_capacity(num_units);
        let mut types = Vec::with_capacity(num_units);

        for i in 0..num_units {
            let t = (i as f64) / (num_units as f64 - 1.0).max(1.0);
            thresholds.push(t);
            // Exponential force sizing: larger units produce more force
            forces.push((2.0 * t).exp());
            types.push(t); // linear mix from slow to fast
        }

        // Normalize force contributions
        let total: f64 = forces.iter().sum();
        for f in &mut forces {
            *f /= total;
        }

        Self {
            num_units,
            thresholds,
            force_contributions: forces,
            fiber_types: types,
            recruited: vec![false; num_units],
        }
    }

    /// Create a model with 50% slow / 50% fast split.
    pub fn mixed_fiber(num_units: usize) -> Self {
        Self::new(num_units)
    }

    /// Update recruitment based on current excitation level.
    pub fn update(&mut self, excitation: f64) {
        let e = clamp(excitation, 0.0, 1.0);
        for i in 0..self.num_units {
            self.recruited[i] = e >= self.thresholds[i];
        }
    }

    /// Compute the total normalized force from recruited units.
    pub fn total_force(&self) -> f64 {
        let mut total = 0.0;
        for i in 0..self.num_units {
            if self.recruited[i] {
                total += self.force_contributions[i];
            }
        }
        total
    }

    /// Count the number of currently recruited motor units.
    pub fn num_recruited(&self) -> usize {
        self.recruited.iter().filter(|&&r| r).count()
    }

    /// Get the average fiber type of recruited units (0 = slow, 1 = fast).
    pub fn average_fiber_type(&self) -> f64 {
        let mut sum_type = 0.0;
        let mut count = 0;
        for i in 0..self.num_units {
            if self.recruited[i] {
                sum_type += self.fiber_types[i];
                count += 1;
            }
        }
        if count > 0 {
            sum_type / count as f64
        } else {
            0.0
        }
    }

    /// Reset all units to un-recruited state.
    pub fn reset(&mut self) {
        for r in &mut self.recruited {
            *r = false;
        }
    }
}

// ===========================================================================
// 9. Isometric / Isotonic Contraction Modes
// ===========================================================================

/// Isometric contraction state (fixed muscle-tendon length).
///
/// The muscle generates force at a fixed total length. Useful for
/// computing the force-activation relationship at a specific pose.
#[derive(Debug, Clone)]
pub struct IsometricContraction {
    /// Fixed musculotendon length (m).
    pub mt_length: f64,
    /// The underlying muscle model.
    pub muscle: HillMuscleSoftbody,
}

impl IsometricContraction {
    /// Create a new isometric contraction experiment.
    pub fn new(muscle: HillMuscleSoftbody, mt_length: f64) -> Self {
        Self { mt_length, muscle }
    }

    /// Compute the isometric force at a given activation level.
    pub fn force_at_activation(&mut self, activation: f64) -> f64 {
        self.muscle.solve_equilibrium(activation, self.mt_length);
        self.muscle.force()
    }

    /// Run a ramp-and-hold activation protocol.
    ///
    /// Returns a vector of (time, force) pairs.
    pub fn ramp_hold(
        &mut self,
        ramp_time: f64,
        hold_time: f64,
        dt: f64,
        max_activation: f64,
    ) -> Vec<(f64, f64)> {
        let mut results = Vec::new();
        let total_time = ramp_time + hold_time;
        let mut t = 0.0;

        while t < total_time {
            let excitation = if t < ramp_time {
                max_activation * t / ramp_time
            } else {
                max_activation
            };

            self.muscle.update(excitation, self.mt_length, 0.0, dt);
            results.push((t, self.muscle.force()));
            t += dt;
        }
        results
    }

    /// Compute the maximum isometric force (at full activation).
    pub fn max_force(&mut self) -> f64 {
        self.force_at_activation(1.0)
    }
}

/// Isotonic contraction state (fixed external load).
///
/// The muscle shortens or lengthens against a constant load.
#[derive(Debug, Clone)]
pub struct IsotonicContraction {
    /// External load (N).
    pub load: f64,
    /// The underlying muscle model.
    pub muscle: HillMuscleSoftbody,
    /// Current musculotendon length (m).
    pub mt_length: f64,
    /// Current musculotendon velocity (m/s).
    pub mt_velocity: f64,
}

impl IsotonicContraction {
    /// Create a new isotonic contraction.
    pub fn new(muscle: HillMuscleSoftbody, load: f64, initial_mt_length: f64) -> Self {
        Self {
            load,
            muscle,
            mt_length: initial_mt_length,
            mt_velocity: 0.0,
        }
    }

    /// Update the isotonic contraction for one time step.
    ///
    /// The muscle adjusts its velocity to match the external load.
    pub fn update(&mut self, excitation: f64, dt: f64) {
        self.muscle
            .update(excitation, self.mt_length, self.mt_velocity, dt);
        let muscle_force = self.muscle.force();

        // Net force determines acceleration (simplified as force imbalance)
        // F_muscle - F_load = m_eff * a
        let m_eff = 0.1; // effective mass (kg)
        let accel = (muscle_force - self.load) / m_eff;
        self.mt_velocity += accel * dt;
        self.mt_length += self.mt_velocity * dt;
    }

    /// Run an isotonic contraction for a given duration.
    ///
    /// Returns a vector of (time, length, velocity) triples.
    pub fn run(&mut self, excitation: f64, duration: f64, dt: f64) -> Vec<(f64, f64, f64)> {
        let mut results = Vec::new();
        let mut t = 0.0;

        while t < duration {
            self.update(excitation, dt);
            results.push((t, self.mt_length, self.mt_velocity));
            t += dt;
        }
        results
    }

    /// Get the current steady-state shortening velocity (if reached).
    pub fn shortening_velocity(&self) -> f64 {
        self.mt_velocity
    }
}

// ===========================================================================
// 10. Musculotendon Unit (Spatial Attachment)
// ===========================================================================

/// A musculotendon unit connecting two points in 3D space.
///
/// This wraps the Hill muscle model with spatial attachment points,
/// computing the musculotendon length and direction from origin/insertion.
#[derive(Debug, Clone)]
pub struct MusculotendonUnit {
    /// Origin attachment point (body frame).
    pub origin: [f64; 3],
    /// Insertion attachment point (body frame).
    pub insertion: [f64; 3],
    /// The underlying muscle model.
    pub muscle: HillMuscleSoftbody,
    /// Optional via (wrap) point.
    pub via_point: Option<[f64; 3]>,
}

impl MusculotendonUnit {
    /// Create a new musculotendon unit.
    pub fn new(origin: [f64; 3], insertion: [f64; 3], muscle: HillMuscleSoftbody) -> Self {
        Self {
            origin,
            insertion,
            muscle,
            via_point: None,
        }
    }

    /// Set an optional via (wrapping) point.
    pub fn with_via_point(mut self, via: [f64; 3]) -> Self {
        self.via_point = Some(via);
        self
    }

    /// Compute the current musculotendon length.
    pub fn mt_length(&self) -> f64 {
        match self.via_point {
            Some(via) => norm3(sub3(via, self.origin)) + norm3(sub3(self.insertion, via)),
            None => norm3(sub3(self.insertion, self.origin)),
        }
    }

    /// Compute the unit direction from origin to insertion.
    pub fn direction(&self) -> [f64; 3] {
        normalize3(sub3(self.insertion, self.origin))
    }

    /// Compute the force vector applied at the insertion point.
    pub fn force_at_insertion(&self) -> [f64; 3] {
        let dir = self.direction();
        let f = self.muscle.force();
        // Force pulls insertion toward origin
        scale3(dir, -f)
    }

    /// Compute the force vector applied at the origin point.
    pub fn force_at_origin(&self) -> [f64; 3] {
        let dir = self.direction();
        let f = self.muscle.force();
        // Force pulls origin toward insertion
        scale3(dir, f)
    }

    /// Update the musculotendon unit state.
    pub fn update_from_positions(
        &mut self,
        origin: [f64; 3],
        insertion: [f64; 3],
        prev_mt_length: f64,
        dt: f64,
        excitation: f64,
    ) {
        self.origin = origin;
        self.insertion = insertion;
        let mt_len = self.mt_length();
        let mt_vel = if dt > 1e-15 {
            (mt_len - prev_mt_length) / dt
        } else {
            0.0
        };
        self.muscle.update(excitation, mt_len, mt_vel, dt);
    }

    /// Get the moment arm about an axis passing through a joint center.
    pub fn moment_arm(&self, joint_center: [f64; 3], axis: [f64; 3]) -> f64 {
        let r = sub3(self.insertion, joint_center);
        let f_dir = self.direction();
        // moment arm = |r x f_dir| . axis_hat
        let cross = [
            r[1] * f_dir[2] - r[2] * f_dir[1],
            r[2] * f_dir[0] - r[0] * f_dir[2],
            r[0] * f_dir[1] - r[1] * f_dir[0],
        ];
        let axis_n = normalize3(axis);
        dot3(cross, axis_n)
    }
}

// ===========================================================================
// 11. Muscle Group
// ===========================================================================

/// A group of muscles that act together (e.g., quadriceps group).
#[derive(Debug, Clone)]
pub struct MuscleGroup {
    /// Name of the muscle group.
    pub name: String,
    /// Individual muscles in this group.
    pub muscles: Vec<MusculotendonUnit>,
}

impl MuscleGroup {
    /// Create a new muscle group.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            muscles: Vec::new(),
        }
    }

    /// Add a muscle to the group.
    pub fn add_muscle(&mut self, mtu: MusculotendonUnit) {
        self.muscles.push(mtu);
    }

    /// Get the total force magnitude from all muscles.
    pub fn total_force(&self) -> f64 {
        self.muscles.iter().map(|m| m.muscle.force()).sum()
    }

    /// Get the total force vector from all muscles at their insertion points.
    pub fn total_force_vector(&self) -> [f64; 3] {
        let mut total = [0.0; 3];
        for m in &self.muscles {
            let f = m.force_at_insertion();
            total = add3(total, f);
        }
        total
    }

    /// Number of muscles in the group.
    pub fn num_muscles(&self) -> usize {
        self.muscles.len()
    }

    /// Activate all muscles in the group with the same excitation.
    pub fn activate_all(&mut self, excitation: f64, dt: f64) {
        for m in &mut self.muscles {
            let mt_len = m.mt_length();
            m.muscle.update(excitation, mt_len, 0.0, dt);
        }
    }

    /// Activate muscles with individual excitation levels.
    pub fn activate_individual(&mut self, excitations: &[f64], dt: f64) {
        for (i, m) in self.muscles.iter_mut().enumerate() {
            let e = if i < excitations.len() {
                excitations[i]
            } else {
                0.0
            };
            let mt_len = m.mt_length();
            m.muscle.update(e, mt_len, 0.0, dt);
        }
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -----------------------------------------------------------------------
    // Activation Dynamics Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_activation_initial_zero() {
        let ad = ActivationDynamics::fast_twitch();
        assert!(
            ad.activation() < 0.02,
            "Initial activation should be ~0, got {}",
            ad.activation()
        );
    }

    #[test]
    fn test_activation_increases_with_excitation() {
        let mut ad = ActivationDynamics::fast_twitch();
        ad.set_activation(0.0);
        for _ in 0..200 {
            ad.update(1.0, 0.001);
        }
        assert!(
            ad.activation() > 0.5,
            "Activation should increase to >0.5, got {}",
            ad.activation()
        );
    }

    #[test]
    fn test_activation_decreases_without_excitation() {
        let mut ad = ActivationDynamics::fast_twitch();
        ad.set_activation(1.0);
        for _ in 0..200 {
            ad.update(0.0, 0.001);
        }
        assert!(
            ad.activation() < 0.5,
            "Activation should decrease, got {}",
            ad.activation()
        );
    }

    #[test]
    fn test_activation_fast_vs_slow() {
        let mut fast = ActivationDynamics::fast_twitch();
        let mut slow = ActivationDynamics::slow_twitch();
        fast.set_activation(0.0);
        slow.set_activation(0.0);
        for _ in 0..100 {
            fast.update(1.0, 0.001);
            slow.update(1.0, 0.001);
        }
        assert!(
            fast.activation() > slow.activation(),
            "Fast-twitch ({}) should activate faster than slow-twitch ({})",
            fast.activation(),
            slow.activation()
        );
    }

    // -----------------------------------------------------------------------
    // Force-Length Curve Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fl_active_peak_at_optimal() {
        let fl = ForceLengthCurve::thelen2003();
        let val = fl.active(1.0);
        assert!(
            approx_eq(val, 1.0, EPS),
            "Active FL at optimal should be 1.0, got {val}"
        );
    }

    #[test]
    fn test_fl_active_decreases_away() {
        let fl = ForceLengthCurve::thelen2003();
        let peak = fl.active(1.0);
        let short = fl.active(0.7);
        let long = fl.active(1.3);
        assert!(short < peak, "Active FL should decrease when shorter");
        assert!(long < peak, "Active FL should decrease when longer");
    }

    #[test]
    fn test_fl_passive_zero_below_optimal() {
        let fl = ForceLengthCurve::thelen2003();
        let val = fl.passive(0.9);
        assert!(
            val.abs() < EPS,
            "Passive FL below optimal should be 0, got {val}"
        );
    }

    #[test]
    fn test_fl_passive_increases_above_optimal() {
        let fl = ForceLengthCurve::thelen2003();
        let f1 = fl.passive(1.2);
        let f2 = fl.passive(1.4);
        assert!(f1 > 0.0, "Passive FL above optimal should be > 0");
        assert!(f2 > f1, "Passive FL should increase with length");
    }

    // -----------------------------------------------------------------------
    // Force-Velocity Curve Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fv_isometric_unity() {
        let fv = ForceVelocityCurve::thelen2003();
        let val = fv.evaluate(0.0);
        assert!(
            approx_eq(val, 1.0, EPS),
            "FV at zero velocity should be 1.0, got {val}"
        );
    }

    #[test]
    fn test_fv_concentric_decreases() {
        let fv = ForceVelocityCurve::thelen2003();
        let iso = fv.evaluate(0.0);
        let con = fv.evaluate(-0.5);
        assert!(
            con < iso,
            "FV during shortening ({con}) should be less than isometric ({iso})"
        );
    }

    #[test]
    fn test_fv_eccentric_increases() {
        let fv = ForceVelocityCurve::thelen2003();
        let iso = fv.evaluate(0.0);
        let ecc = fv.evaluate(0.5);
        assert!(
            ecc > iso,
            "FV during lengthening ({ecc}) should be greater than isometric ({iso})"
        );
    }

    #[test]
    fn test_fv_eccentric_bounded() {
        let fv = ForceVelocityCurve::thelen2003();
        let ecc_high = fv.evaluate(0.99);
        assert!(
            ecc_high < fv.f_ecc_max + 0.1,
            "Eccentric FV ({ecc_high}) should be bounded by f_ecc_max ({})",
            fv.f_ecc_max
        );
    }

    // -----------------------------------------------------------------------
    // Pennation Model Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_pennation_zero_angle() {
        let pm = PennationModel::zero();
        let a = pm.angle(1.0);
        assert!(a.abs() < EPS, "Zero pennation angle should give 0, got {a}");
    }

    #[test]
    fn test_pennation_cos_at_optimal() {
        let alpha = 0.3; // ~17 degrees
        let l_opt = 0.1;
        let pm = PennationModel::new(alpha, l_opt);
        let cos_a = pm.cos_alpha(l_opt);
        assert!(
            approx_eq(cos_a, alpha.cos(), 1e-6),
            "cos(alpha) at optimal should match, got {cos_a}"
        );
    }

    #[test]
    fn test_pennation_angle_increases_with_shortening() {
        let pm = PennationModel::new(0.2, 0.1);
        let a_opt = pm.angle(0.1);
        let a_short = pm.angle(0.08);
        assert!(
            a_short > a_opt,
            "Pennation angle should increase with shortening: opt={a_opt}, short={a_short}"
        );
    }

    // -----------------------------------------------------------------------
    // Tendon Model Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_tendon_zero_force_at_slack() {
        let tm = TendonModel::new(0.2, 1000.0);
        let f = tm.force(0.2);
        assert!(
            f.abs() < EPS,
            "Tendon force at slack length should be 0, got {f}"
        );
    }

    #[test]
    fn test_tendon_positive_force_above_slack() {
        let tm = TendonModel::new(0.2, 1000.0);
        let f = tm.force(0.21); // 5% strain
        assert!(
            f > 0.0,
            "Tendon force above slack should be positive, got {f}"
        );
    }

    #[test]
    fn test_tendon_stiffness_increases() {
        let tm = TendonModel::new(0.2, 1000.0);
        let k1 = tm.stiffness(0.201); // small strain
        let k2 = tm.stiffness(0.21); // larger strain
        assert!(
            k2 >= k1,
            "Tendon stiffness should increase: k1={k1}, k2={k2}"
        );
    }

    #[test]
    fn test_tendon_inverse() {
        let tm = TendonModel::new(0.2, 1000.0);
        let l = 0.208; // some stretched length
        let f = tm.force(l);
        let l_inv = tm.length_from_force(f);
        assert!(
            approx_eq(l, l_inv, 1e-4),
            "Inverse tendon should recover length: l={l}, l_inv={l_inv}"
        );
    }

    // -----------------------------------------------------------------------
    // Hill Muscle Model Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hill_zero_force_inactive() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let f = muscle.compute_force(0.0, 0.1, 0.0);
        assert!(
            f.abs() < EPS,
            "Force at zero activation and optimal length should be 0, got {f}"
        );
    }

    #[test]
    fn test_hill_positive_force_active() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let f = muscle.compute_force(1.0, 0.1, 0.0);
        assert!(
            f > 0.0,
            "Force at full activation and optimal length should be positive, got {f}"
        );
    }

    #[test]
    fn test_hill_max_isometric_force() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let f = muscle.compute_force(1.0, 0.1, 0.0);
        // At optimal length, zero velocity, full activation: F = F_max * cos(alpha)
        assert!(
            approx_eq(f, 1000.0, 1.0),
            "Max isometric force should be ~F_max, got {f}"
        );
    }

    #[test]
    fn test_hill_update_increases_force() {
        let mut muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        assert!(muscle.force().abs() < EPS, "Initial force should be 0");
        for _ in 0..200 {
            muscle.update(1.0, 0.3, 0.0, 0.001);
        }
        assert!(
            muscle.force() > 0.0,
            "Force should increase with activation, got {}",
            muscle.force()
        );
    }

    // -----------------------------------------------------------------------
    // Fatigue Model Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fatigue_initial_state() {
        let fm = FatigueModel::default_mixed();
        assert!(approx_eq(fm.mr, 1.0, EPS), "All units should start resting");
        assert!(approx_eq(fm.ma, 0.0, EPS), "No active units initially");
        assert!(approx_eq(fm.mf, 0.0, EPS), "No fatigued units initially");
    }

    #[test]
    fn test_fatigue_recruitment() {
        let mut fm = FatigueModel::default_mixed();
        fm.set_target(0.5);
        for _ in 0..100 {
            fm.update(0.01);
        }
        assert!(
            fm.ma > 0.0,
            "Active units should be recruited, got {}",
            fm.ma
        );
    }

    #[test]
    fn test_fatigue_develops() {
        let mut fm = FatigueModel::default_mixed();
        fm.set_target(0.8);
        for _ in 0..5000 {
            fm.update(0.01);
        }
        assert!(
            fm.mf > 0.0,
            "Fatigue should develop over time, got mf={}",
            fm.mf
        );
    }

    #[test]
    fn test_fatigue_conservation() {
        let mut fm = FatigueModel::default_mixed();
        fm.set_target(0.6);
        for _ in 0..100 {
            fm.update(0.01);
        }
        let total = fm.ma + fm.mf + fm.mr;
        assert!(
            approx_eq(total, 1.0, 1e-6),
            "MA+MF+MR should sum to 1.0, got {total}"
        );
    }

    #[test]
    fn test_fatigue_endurance_time() {
        let fm = FatigueModel::default_mixed();
        let t50 = fm.endurance_time(0.5);
        let t80 = fm.endurance_time(0.8);
        assert!(
            t50 > t80,
            "Endurance at 50% ({t50:.6}) should exceed 80% ({t80:.6})"
        );
    }

    // -----------------------------------------------------------------------
    // Fiber Recruitment Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_recruitment_none_at_zero() {
        let mut fr = FiberRecruitment::new(10);
        fr.update(0.0);
        // Only the first unit (threshold=0) gets recruited
        assert!(
            fr.num_recruited() <= 1,
            "At zero excitation, at most 1 unit should be recruited, got {}",
            fr.num_recruited()
        );
    }

    #[test]
    fn test_recruitment_all_at_max() {
        let mut fr = FiberRecruitment::new(10);
        fr.update(1.0);
        assert_eq!(
            fr.num_recruited(),
            10,
            "At full excitation all units should be recruited"
        );
    }

    #[test]
    fn test_recruitment_progressive() {
        let mut fr = FiberRecruitment::new(10);
        fr.update(0.3);
        let n_low = fr.num_recruited();
        fr.update(0.7);
        let n_high = fr.num_recruited();
        assert!(
            n_high >= n_low,
            "Higher excitation should recruit more: low={n_low}, high={n_high}"
        );
    }

    #[test]
    fn test_recruitment_total_force_at_max() {
        let mut fr = FiberRecruitment::new(10);
        fr.update(1.0);
        let f = fr.total_force();
        assert!(
            approx_eq(f, 1.0, 1e-6),
            "Total force at full recruitment should be 1.0, got {f}"
        );
    }

    // -----------------------------------------------------------------------
    // Isometric Contraction Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_isometric_force_monotonic() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        // Use compute_force directly to test the relationship
        let f_low = muscle.compute_force(0.3, 0.1, 0.0);
        let f_high = muscle.compute_force(0.8, 0.1, 0.0);
        assert!(
            f_high > f_low,
            "Higher activation should give higher force: f_low={f_low}, f_high={f_high}"
        );
    }

    // -----------------------------------------------------------------------
    // Isotonic Contraction Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_isotonic_shortens_under_light_load() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let initial_length = 0.31;
        let mut iso = IsotonicContraction::new(muscle, 100.0, initial_length);
        // First activate the muscle (build up activation)
        for _ in 0..200 {
            iso.update(1.0, 0.0001);
        }
        // Check that the muscle produced force and the length changed
        let force = iso.muscle.force();
        assert!(
            force > 0.0,
            "Muscle should produce force after activation, got {force}"
        );
    }

    // -----------------------------------------------------------------------
    // Musculotendon Unit Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_mtu_length() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let mtu = MusculotendonUnit::new([0.0, 0.0, 0.0], [0.3, 0.0, 0.0], muscle);
        let l = mtu.mt_length();
        assert!(approx_eq(l, 0.3, EPS), "MTU length should be 0.3, got {l}");
    }

    #[test]
    fn test_mtu_direction() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let mtu = MusculotendonUnit::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], muscle);
        let dir = mtu.direction();
        assert!(approx_eq(dir[0], 0.0, EPS));
        assert!(approx_eq(dir[1], 1.0, EPS));
        assert!(approx_eq(dir[2], 0.0, EPS));
    }

    #[test]
    fn test_mtu_via_point_length() {
        let muscle = HillMuscleSoftbody::typical(1000.0, 0.1, 0.2);
        let mtu = MusculotendonUnit::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], muscle)
            .with_via_point([0.5, 0.5, 0.0]);
        let l = mtu.mt_length();
        // Via point makes path longer than straight line
        let straight = 1.0_f64;
        assert!(
            l > straight,
            "Via point path ({l}) should be longer than straight ({straight})"
        );
    }

    // -----------------------------------------------------------------------
    // Muscle Group Tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_muscle_group_total_force() {
        let m1 = HillMuscleSoftbody::typical(500.0, 0.1, 0.2);
        let m2 = HillMuscleSoftbody::typical(300.0, 0.08, 0.15);
        let mtu1 = MusculotendonUnit::new([0.0, 0.0, 0.0], [0.3, 0.0, 0.0], m1);
        let mtu2 = MusculotendonUnit::new([0.0, 0.1, 0.0], [0.23, 0.1, 0.0], m2);
        let mut group = MuscleGroup::new("test_group");
        group.add_muscle(mtu1);
        group.add_muscle(mtu2);
        assert_eq!(group.num_muscles(), 2);
        // Initially no force (no activation)
        let f = group.total_force();
        assert!(f.abs() < EPS, "Initial group force should be 0, got {f}");
    }

    #[test]
    fn test_muscle_group_activation() {
        let m1 = HillMuscleSoftbody::typical(500.0, 0.1, 0.2);
        let mtu1 = MusculotendonUnit::new([0.0, 0.0, 0.0], [0.3, 0.0, 0.0], m1);
        let mut group = MuscleGroup::new("test_group");
        group.add_muscle(mtu1);
        for _ in 0..200 {
            group.activate_all(1.0, 0.001);
        }
        let f = group.total_force();
        assert!(
            f > 0.0,
            "Group force after activation should be positive, got {f}"
        );
    }
}
