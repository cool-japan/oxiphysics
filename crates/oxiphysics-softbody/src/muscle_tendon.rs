// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Musculotendon dynamics: Hill-type muscle model, activation dynamics,
//! tendon compliance, force-velocity/force-length relations, pennation
//! angle, muscle redundancy solver, and analysis utilities.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// Constants from std::f64 used inline where needed.

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Maximum isometric muscle force scaling factor (dimensionless).
const DEFAULT_F_MAX: f64 = 1.0;
/// Optimal fiber length default (m).
const DEFAULT_L_OPT: f64 = 0.10;
/// Tendon slack length default (m).
const DEFAULT_L_SLACK: f64 = 0.20;
/// Default pennation angle at optimal length (rad).
const DEFAULT_PENNATION: f64 = 0.0;
/// Hill constant a/F0 ratio (typical skeletal muscle).
const DEFAULT_A_REL: f64 = 0.25;
/// Hill constant b/L_opt ratio (typical skeletal muscle).
const DEFAULT_B_REL: f64 = 0.25;
/// Maximum shortening velocity in L_opt/s.
const DEFAULT_V_MAX: f64 = 10.0;
/// Default activation time constant (s).
const DEFAULT_TAU_ACT: f64 = 0.015;
/// Default deactivation time constant (s).
const DEFAULT_TAU_DEACT: f64 = 0.060;
/// Width of the Gaussian active force-length curve.
const DEFAULT_FL_WIDTH: f64 = 0.56;
/// Shape factor for passive force-length curve.
const DEFAULT_KPE: f64 = 4.0;
/// Passive force strain threshold.
const DEFAULT_E0: f64 = 0.6;
/// Tendon strain at maximum isometric force.
const DEFAULT_E_T: f64 = 0.033;
/// Tendon stiffness shape factor.
const DEFAULT_KT: f64 = 35.0;

// ---------------------------------------------------------------------------
// ActivationDynamics
// ---------------------------------------------------------------------------

/// First-order activation dynamics: da/dt = (u - a) / tau(u, a).
///
/// Models the calcium-driven excitation-contraction coupling delay.
/// The time constant switches between activation (tau_act) when
/// excitation exceeds activation and deactivation (tau_deact) otherwise.
#[derive(Debug, Clone)]
pub struct ActivationDynamics {
    /// Current activation level \[0, 1\].
    pub activation: f64,
    /// Neural excitation input \[0, 1\].
    pub excitation: f64,
    /// Activation time constant (s).
    pub tau_act: f64,
    /// Deactivation time constant (s).
    pub tau_deact: f64,
    /// Minimum activation floor.
    pub a_min: f64,
    /// Accumulated calcium-like intermediate state \[0, 1\].
    pub calcium: f64,
    /// Calcium release rate (1/s).
    pub calcium_release_rate: f64,
    /// Calcium uptake rate (1/s).
    pub calcium_uptake_rate: f64,
}

impl ActivationDynamics {
    /// Create new activation dynamics with default time constants.
    pub fn new() -> Self {
        Self {
            activation: 0.0,
            excitation: 0.0,
            tau_act: DEFAULT_TAU_ACT,
            tau_deact: DEFAULT_TAU_DEACT,
            a_min: 0.01,
            calcium: 0.0,
            calcium_release_rate: 200.0,
            calcium_uptake_rate: 30.0,
        }
    }

    /// Create activation dynamics with custom time constants.
    pub fn with_time_constants(tau_act: f64, tau_deact: f64) -> Self {
        Self {
            tau_act,
            tau_deact,
            ..Self::new()
        }
    }

    /// Set excitation input \[0, 1\].
    pub fn set_excitation(&mut self, u: f64) {
        self.excitation = u.clamp(0.0, 1.0);
    }

    /// Effective time constant based on current excitation and activation.
    pub fn effective_tau(&self) -> f64 {
        if self.excitation > self.activation {
            self.tau_act * (0.5 + 1.5 * self.activation)
        } else {
            self.tau_deact / (0.5 + 1.5 * self.activation)
        }
    }

    /// Step the activation ODE forward by dt seconds (explicit Euler).
    pub fn step(&mut self, dt: f64) {
        let tau = self.effective_tau();
        if tau > 1e-15 {
            let da = (self.excitation - self.activation) / tau;
            self.activation += da * dt;
        }
        self.activation = self.activation.clamp(self.a_min, 1.0);
    }

    /// Step with calcium dynamics intermediate.
    pub fn step_calcium(&mut self, dt: f64) {
        // Calcium release/uptake
        let d_calcium = if self.excitation > 0.5 {
            self.calcium_release_rate * (1.0 - self.calcium)
        } else {
            -self.calcium_uptake_rate * self.calcium
        };
        self.calcium += d_calcium * dt;
        self.calcium = self.calcium.clamp(0.0, 1.0);

        // Activation follows calcium
        let tau = self.effective_tau();
        if tau > 1e-15 {
            let da = (self.calcium - self.activation) / tau;
            self.activation += da * dt;
        }
        self.activation = self.activation.clamp(self.a_min, 1.0);
    }

    /// Steady-state activation for a given excitation.
    pub fn steady_state(&self) -> f64 {
        self.excitation.clamp(self.a_min, 1.0)
    }

    /// Reset activation and calcium to zero.
    pub fn reset(&mut self) {
        self.activation = self.a_min;
        self.calcium = 0.0;
        self.excitation = 0.0;
    }
}

impl Default for ActivationDynamics {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ForceVelocityRelation
// ---------------------------------------------------------------------------

/// Hill's hyperbolic force-velocity relationship.
///
/// Concentric: F = F_max * (v_max - v) / (v_max + v / a_rel)
/// Eccentric: extended beyond isometric force with asymptotic limit.
#[derive(Debug, Clone)]
pub struct ForceVelocityRelation {
    /// Maximum shortening velocity (L_opt/s).
    pub v_max: f64,
    /// Hill's a parameter relative to F_max (dimensionless).
    pub a_rel: f64,
    /// Hill's b parameter relative to L_opt (L_opt/s).
    pub b_rel: f64,
    /// Eccentric force enhancement factor (typically 1.4-1.8).
    pub f_eccentric_max: f64,
    /// Eccentric asymptote shape factor.
    pub eccentric_shape: f64,
}

impl ForceVelocityRelation {
    /// Create with default Hill parameters.
    pub fn new() -> Self {
        Self {
            v_max: DEFAULT_V_MAX,
            a_rel: DEFAULT_A_REL,
            b_rel: DEFAULT_B_REL,
            f_eccentric_max: 1.5,
            eccentric_shape: 0.18,
        }
    }

    /// Create with custom parameters.
    pub fn with_params(v_max: f64, a_rel: f64, b_rel: f64) -> Self {
        Self {
            v_max,
            a_rel,
            b_rel,
            ..Self::new()
        }
    }

    /// Normalized force multiplier for a given normalized velocity.
    ///
    /// `v_norm`: fiber velocity normalized by v_max. Negative = shortening, positive = lengthening.
    ///
    /// Returns force multiplier \[0, f_eccentric_max\].
    pub fn evaluate(&self, v_norm: f64) -> f64 {
        if v_norm <= 0.0 {
            // Concentric (shortening): Hill's equation
            let v_bar = (-v_norm).min(0.999);
            let num = 1.0 - v_bar;
            let den = 1.0 + v_bar / self.a_rel;
            if den.abs() < 1e-15 {
                0.0
            } else {
                (num / den).max(0.0)
            }
        } else {
            // Eccentric (lengthening)
            let v_bar = v_norm;
            let f_asym = self.f_eccentric_max;
            let shape = self.eccentric_shape;
            // Smooth asymptotic rise: F = 1 + (f_asym - 1) * v / (v + shape)
            let f = 1.0 + (f_asym - 1.0) * v_bar / (v_bar + shape);
            f.min(f_asym)
        }
    }

    /// Inverse: given force multiplier, compute normalized velocity.
    pub fn inverse(&self, f_norm: f64) -> f64 {
        if f_norm <= 0.0 {
            return -1.0; // max shortening
        }
        if f_norm >= 1.0 {
            // Eccentric region
            if f_norm >= self.f_eccentric_max {
                return 10.0; // saturated
            }
            let f_asym = self.f_eccentric_max;
            let shape = self.eccentric_shape;
            // f = 1 + (f_asym-1)*v/(v+shape) => v = shape*(f-1)/(f_asym-f)
            let denom = f_asym - f_norm;
            if denom < 1e-15 {
                return 10.0;
            }
            shape * (f_norm - 1.0) / denom
        } else {
            // Concentric: f = (1-v)/(1+v/a_rel) => v = (1-f)/(1+f/a_rel)... rearranged
            let denom = 1.0 + f_norm / self.a_rel;
            if denom < 1e-15 {
                return -1.0;
            }
            let v_bar = (1.0 - f_norm) / denom;
            -v_bar
        }
    }

    /// Isometric force multiplier (v=0).
    pub fn isometric(&self) -> f64 {
        self.evaluate(0.0)
    }
}

impl Default for ForceVelocityRelation {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// ForceLengthRelation
// ---------------------------------------------------------------------------

/// Active and passive force-length relationships for muscle fibers.
///
/// Active curve: Gaussian centered at L_opt.
/// Passive curve: exponential rise beyond L_opt.
#[derive(Debug, Clone)]
pub struct ForceLengthRelation {
    /// Width of the Gaussian active force-length curve.
    pub width: f64,
    /// Passive elastic shape factor.
    pub kpe: f64,
    /// Passive strain threshold (fraction of L_opt).
    pub e0: f64,
    /// Active curve skewness factor (0 = symmetric).
    pub skewness: f64,
}

impl ForceLengthRelation {
    /// Create with default parameters.
    pub fn new() -> Self {
        Self {
            width: DEFAULT_FL_WIDTH,
            kpe: DEFAULT_KPE,
            e0: DEFAULT_E0,
            skewness: 0.0,
        }
    }

    /// Create with custom width and passive parameters.
    pub fn with_params(width: f64, kpe: f64, e0: f64) -> Self {
        Self {
            width,
            kpe,
            e0,
            skewness: 0.0,
        }
    }

    /// Active force-length multiplier.
    ///
    /// `l_norm`: fiber length normalized by L_opt.
    /// Returns value in \[0, 1\] with peak at l_norm = 1.
    pub fn active(&self, l_norm: f64) -> f64 {
        let x = (l_norm - 1.0) / self.width;
        let skew_offset = self.skewness * (l_norm - 1.0);
        let val = (-(x + skew_offset) * (x + skew_offset)).exp();
        val.clamp(0.0, 1.0)
    }

    /// Passive force-length multiplier.
    ///
    /// `l_norm`: fiber length normalized by L_opt.
    /// Returns 0 below slack and exponential rise above.
    pub fn passive(&self, l_norm: f64) -> f64 {
        if l_norm <= 1.0 {
            return 0.0;
        }
        let strain = l_norm - 1.0;
        if strain <= 0.0 {
            return 0.0;
        }
        let val = (self.kpe * strain / self.e0).exp() - 1.0;
        let denom = self.kpe.exp() - 1.0;
        if denom.abs() < 1e-15 {
            return 0.0;
        }
        (val / denom).max(0.0)
    }

    /// Total force-length multiplier (active * activation + passive).
    pub fn total(&self, l_norm: f64, activation: f64) -> f64 {
        activation * self.active(l_norm) + self.passive(l_norm)
    }

    /// Length at which active force is maximal (always 1.0 in normalized units).
    pub fn optimal_length_norm(&self) -> f64 {
        1.0
    }

    /// Range of normalized lengths where active force exceeds a threshold.
    pub fn active_range(&self, threshold: f64) -> (f64, f64) {
        // Gaussian: exp(-x^2/w^2) >= threshold => |x| <= w*sqrt(-ln(threshold))
        let threshold = threshold.clamp(1e-10, 0.999);
        let half_width = self.width * (-threshold.ln()).sqrt();
        (1.0 - half_width, 1.0 + half_width)
    }
}

impl Default for ForceLengthRelation {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TendonModel
// ---------------------------------------------------------------------------

/// Nonlinear tendon force-length model.
///
/// The tendon is modeled as a nonlinear elastic element with a slack length.
/// Below the slack length, force is zero. Above it, force rises exponentially.
#[derive(Debug, Clone)]
pub struct TendonModel {
    /// Tendon slack length (m).
    pub l_slack: f64,
    /// Maximum isometric force that the tendon can transmit (N).
    pub f_max: f64,
    /// Tendon strain at maximum isometric force.
    pub e_t: f64,
    /// Shape factor for stiffness curve.
    pub k_t: f64,
    /// Linear stiffness region threshold (strain).
    pub linear_threshold: f64,
}

impl TendonModel {
    /// Create tendon with given slack length and maximum force.
    pub fn new(l_slack: f64, f_max: f64) -> Self {
        Self {
            l_slack,
            f_max,
            e_t: DEFAULT_E_T,
            k_t: DEFAULT_KT,
            linear_threshold: 0.609 * DEFAULT_E_T,
        }
    }

    /// Create with custom stiffness parameters.
    pub fn with_params(l_slack: f64, f_max: f64, e_t: f64, k_t: f64) -> Self {
        Self {
            l_slack,
            f_max,
            e_t,
            k_t,
            linear_threshold: 0.609 * e_t,
        }
    }

    /// Tendon strain from current length.
    pub fn strain(&self, l_tendon: f64) -> f64 {
        if self.l_slack <= 1e-15 {
            return 0.0;
        }
        (l_tendon - self.l_slack) / self.l_slack
    }

    /// Tendon force given current tendon length.
    pub fn force(&self, l_tendon: f64) -> f64 {
        let e = self.strain(l_tendon);
        if e <= 0.0 {
            return 0.0;
        }
        self.force_from_strain(e)
    }

    /// Tendon force from strain directly.
    pub fn force_from_strain(&self, e: f64) -> f64 {
        if e <= 0.0 {
            return 0.0;
        }
        if e <= self.linear_threshold {
            // Toe region: quadratic
            let c1 = self.f_max / (self.e_t * self.e_t);
            let f = c1 * e * e;
            f.max(0.0)
        } else {
            // Linear region
            let c1 = self.f_max / (self.e_t * self.e_t);
            let f_at_threshold = c1 * self.linear_threshold * self.linear_threshold;
            let k_lin = 2.0 * c1 * self.linear_threshold;
            let f = f_at_threshold + k_lin * (e - self.linear_threshold);
            f.max(0.0)
        }
    }

    /// Tendon stiffness (dF/dL) at current length.
    pub fn stiffness(&self, l_tendon: f64) -> f64 {
        if self.l_slack <= 1e-15 {
            return 0.0;
        }
        let e = self.strain(l_tendon);
        if e <= 0.0 {
            return 0.0;
        }
        let c1 = self.f_max / (self.e_t * self.e_t);
        let de_dl = 1.0 / self.l_slack;
        if e <= self.linear_threshold {
            2.0 * c1 * e * de_dl
        } else {
            2.0 * c1 * self.linear_threshold * de_dl
        }
    }

    /// Length at which tendon force equals a given fraction of F_max.
    pub fn length_at_force_fraction(&self, frac: f64) -> f64 {
        let frac = frac.clamp(0.0, 5.0);
        let target = frac * self.f_max;
        let c1 = self.f_max / (self.e_t * self.e_t);
        let f_threshold = c1 * self.linear_threshold * self.linear_threshold;

        let e = if target <= f_threshold {
            (target / c1).sqrt()
        } else {
            let k_lin = 2.0 * c1 * self.linear_threshold;
            self.linear_threshold + (target - f_threshold) / k_lin
        };
        self.l_slack * (1.0 + e)
    }

    /// Check if tendon is slack (no tension).
    pub fn is_slack(&self, l_tendon: f64) -> bool {
        l_tendon <= self.l_slack
    }
}

// ---------------------------------------------------------------------------
// HillMuscleModel
// ---------------------------------------------------------------------------

/// Hill-type three-element muscle model.
///
/// Consists of:
/// - CE: contractile element (activation-dependent, force-velocity, force-length)
/// - PE: parallel elastic element (passive force-length)
/// - SE: series elastic element (tendon)
///
/// The muscle force is: F_muscle = (a * f_active_FL * f_FV + f_passive_FL) * F_max
#[derive(Debug, Clone)]
pub struct HillMuscleModel {
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fiber length (m).
    pub l_opt: f64,
    /// Force-velocity relation.
    pub fv: ForceVelocityRelation,
    /// Force-length relation.
    pub fl: ForceLengthRelation,
    /// Current fiber length (m).
    pub l_fiber: f64,
    /// Current fiber velocity (m/s, negative = shortening).
    pub v_fiber: f64,
    /// Damping coefficient (N·s/m).
    pub damping: f64,
}

impl HillMuscleModel {
    /// Create a new Hill muscle model.
    pub fn new(f_max: f64, l_opt: f64) -> Self {
        Self {
            f_max,
            l_opt,
            fv: ForceVelocityRelation::new(),
            fl: ForceLengthRelation::new(),
            l_fiber: l_opt,
            v_fiber: 0.0,
            damping: 0.001,
        }
    }

    /// Normalized fiber length.
    pub fn l_norm(&self) -> f64 {
        if self.l_opt > 1e-15 {
            self.l_fiber / self.l_opt
        } else {
            1.0
        }
    }

    /// Normalized fiber velocity.
    pub fn v_norm(&self) -> f64 {
        if self.l_opt > 1e-15 {
            self.v_fiber / (self.fv.v_max * self.l_opt)
        } else {
            0.0
        }
    }

    /// Active force component.
    pub fn active_force(&self, activation: f64) -> f64 {
        let fl_mult = self.fl.active(self.l_norm());
        let fv_mult = self.fv.evaluate(self.v_norm());
        activation * fl_mult * fv_mult * self.f_max
    }

    /// Passive force component.
    pub fn passive_force(&self) -> f64 {
        self.fl.passive(self.l_norm()) * self.f_max
    }

    /// Damping force component.
    pub fn damping_force(&self) -> f64 {
        self.damping * self.v_fiber * self.f_max
    }

    /// Total muscle fiber force.
    pub fn total_force(&self, activation: f64) -> f64 {
        self.active_force(activation) + self.passive_force() + self.damping_force()
    }

    /// Set fiber state directly.
    pub fn set_state(&mut self, l_fiber: f64, v_fiber: f64) {
        self.l_fiber = l_fiber.max(0.001 * self.l_opt);
        self.v_fiber = v_fiber;
    }

    /// Compute contractile element (CE) power.
    pub fn ce_power(&self, activation: f64) -> f64 {
        self.active_force(activation) * self.v_fiber
    }

    /// Instantaneous stiffness dF/dL (numerical).
    pub fn stiffness(&self, activation: f64) -> f64 {
        let dl = 1e-6 * self.l_opt;
        let l_save = self.l_fiber;

        let mut m1 = self.clone();
        m1.l_fiber = l_save + dl;
        let f1 = m1.total_force(activation);

        let mut m2 = self.clone();
        m2.l_fiber = l_save - dl;
        let f2 = m2.total_force(activation);

        (f1 - f2) / (2.0 * dl)
    }
}

// ---------------------------------------------------------------------------
// MusculotendonUnit
// ---------------------------------------------------------------------------

/// Musculotendon unit: muscle fibers in series with a tendon.
///
/// Accounts for pennation angle, tendon compliance, and
/// force equilibrium between muscle and tendon.
#[derive(Debug, Clone)]
pub struct MusculotendonUnit {
    /// Hill muscle model.
    pub muscle: HillMuscleModel,
    /// Tendon model.
    pub tendon: TendonModel,
    /// Activation dynamics.
    pub activation: ActivationDynamics,
    /// Pennation angle at optimal fiber length (rad).
    pub pennation_opt: f64,
    /// Total musculotendon length (m).
    pub l_mt: f64,
    /// Name identifier.
    pub name: String,
    /// Maximum number of Newton iterations for equilibrium.
    pub max_iters: usize,
    /// Newton convergence tolerance.
    pub tol: f64,
}

impl MusculotendonUnit {
    /// Create a new musculotendon unit.
    pub fn new(name: &str, f_max: f64, l_opt: f64, l_slack: f64, pennation_opt: f64) -> Self {
        let l_mt = l_opt * pennation_opt.cos() + l_slack;
        Self {
            muscle: HillMuscleModel::new(f_max, l_opt),
            tendon: TendonModel::new(l_slack, f_max),
            activation: ActivationDynamics::new(),
            pennation_opt,
            l_mt,
            name: name.to_string(),
            max_iters: 100,
            tol: 1e-6,
        }
    }

    /// Current pennation angle from fiber length.
    pub fn pennation_angle(&self) -> f64 {
        let sin_pen = self.pennation_opt.sin() * self.muscle.l_opt / self.muscle.l_fiber;
        let sin_pen = sin_pen.clamp(-1.0, 1.0);
        sin_pen.asin()
    }

    /// Fiber force projected along tendon direction.
    pub fn projected_fiber_force(&self) -> f64 {
        let pen = self.pennation_angle();
        self.muscle.total_force(self.activation.activation) * pen.cos()
    }

    /// Tendon length from total MT length and fiber length.
    pub fn tendon_length(&self) -> f64 {
        let pen = self.pennation_angle();
        self.l_mt - self.muscle.l_fiber * pen.cos()
    }

    /// Set the total musculotendon path length.
    pub fn set_length(&mut self, l_mt: f64) {
        self.l_mt = l_mt;
    }

    /// Solve static equilibrium: fiber force * cos(pen) = tendon force.
    ///
    /// Uses Newton-Raphson to find the fiber length that satisfies
    /// force equilibrium at the current activation and MT length.
    pub fn solve_equilibrium(&mut self) {
        let a = self.activation.activation;
        let mut l_fiber = self.muscle.l_fiber;

        for _iter in 0..self.max_iters {
            // Compute pennation
            let sin_pen = (self.pennation_opt.sin() * self.muscle.l_opt / l_fiber).clamp(-1.0, 1.0);
            let pen = sin_pen.asin();
            let cos_pen = pen.cos();

            // Fiber force
            self.muscle.l_fiber = l_fiber;
            let f_fiber = self.muscle.total_force(a);
            let f_fiber_proj = f_fiber * cos_pen;

            // Tendon force
            let l_t = self.l_mt - l_fiber * cos_pen;
            let f_tendon = self.tendon.force(l_t);

            // Residual
            let residual = f_fiber_proj - f_tendon;
            if residual.abs() < self.tol {
                break;
            }

            // Numerical Jacobian
            let dl = 1e-7 * self.muscle.l_opt;
            let l_plus = l_fiber + dl;
            let sin_pen_p =
                (self.pennation_opt.sin() * self.muscle.l_opt / l_plus).clamp(-1.0, 1.0);
            let pen_p = sin_pen_p.asin();
            let cos_pen_p = pen_p.cos();
            self.muscle.l_fiber = l_plus;
            let f_fiber_p = self.muscle.total_force(a) * cos_pen_p;
            let l_t_p = self.l_mt - l_plus * cos_pen_p;
            let f_tendon_p = self.tendon.force(l_t_p);
            let residual_p = f_fiber_p - f_tendon_p;

            let drdl = (residual_p - residual) / dl;
            if drdl.abs() < 1e-20 {
                break;
            }

            let delta = -residual / drdl;
            l_fiber += delta;
            l_fiber = l_fiber.clamp(0.01 * self.muscle.l_opt, 2.0 * self.muscle.l_opt);
        }

        self.muscle.l_fiber = l_fiber;
    }

    /// Step the musculotendon unit forward in time.
    pub fn step(&mut self, dt: f64) {
        self.activation.step(dt);
        self.solve_equilibrium();
    }

    /// Get the current tendon force.
    pub fn tendon_force(&self) -> f64 {
        self.tendon.force(self.tendon_length())
    }

    /// Get the current total output force along the tendon.
    pub fn output_force(&self) -> f64 {
        self.projected_fiber_force()
    }

    /// Normalized fiber length.
    pub fn normalized_fiber_length(&self) -> f64 {
        self.muscle.l_norm()
    }

    /// Check if tendon is slack.
    pub fn is_tendon_slack(&self) -> bool {
        self.tendon.is_slack(self.tendon_length())
    }
}

// ---------------------------------------------------------------------------
// MuscleRedundancy
// ---------------------------------------------------------------------------

/// Static optimization to resolve muscle redundancy.
///
/// Minimizes the sum of squared activations subject to
/// the constraint that the net joint torque equals the desired torque.
///
/// min sum(a_i^2) s.t. sum(a_i * F_max_i * r_i) = tau_desired
#[derive(Debug, Clone)]
pub struct MuscleRedundancy {
    /// Number of muscles.
    pub n_muscles: usize,
    /// Maximum isometric force for each muscle (N).
    pub f_max: Vec<f64>,
    /// Moment arm for each muscle about the joint (m).
    pub moment_arms: Vec<f64>,
    /// Lower bound on activation.
    pub a_min: f64,
    /// Upper bound on activation.
    pub a_max: f64,
    /// Maximum iterations for the optimizer.
    pub max_iters: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

impl MuscleRedundancy {
    /// Create a new muscle redundancy solver.
    pub fn new(f_max: Vec<f64>, moment_arms: Vec<f64>) -> Self {
        let n = f_max.len();
        Self {
            n_muscles: n,
            f_max,
            moment_arms,
            a_min: 0.01,
            a_max: 1.0,
            max_iters: 200,
            tol: 1e-8,
        }
    }

    /// Solve for optimal activations given desired joint torque.
    ///
    /// Uses Lagrange multiplier approach for quadratic cost with linear constraint.
    /// Returns activation vector.
    pub fn solve(&self, tau_desired: f64) -> Vec<f64> {
        // Quadratic optimization: min sum(a_i^2) s.t. sum(a_i * c_i) = tau
        // where c_i = F_max_i * r_i
        let coeffs: Vec<f64> = self
            .f_max
            .iter()
            .zip(self.moment_arms.iter())
            .map(|(f, r)| f * r)
            .collect();

        // Lagrange multiplier: a_i = lambda * c_i
        // sum(lambda * c_i^2) = tau => lambda = tau / sum(c_i^2)
        let sum_c2: f64 = coeffs.iter().map(|c| c * c).sum();
        if sum_c2 < 1e-20 {
            return vec![0.0; self.n_muscles];
        }

        let lambda = tau_desired / sum_c2;
        let mut activations: Vec<f64> = coeffs.iter().map(|c| lambda * c).collect();

        // Project to bounds and iterate
        for _iter in 0..self.max_iters {
            // Clamp activations
            for a in activations.iter_mut() {
                *a = a.clamp(self.a_min, self.a_max);
            }

            // Check constraint satisfaction
            let tau_current: f64 = activations
                .iter()
                .zip(coeffs.iter())
                .map(|(a, c)| a * c)
                .sum();
            let error = tau_desired - tau_current;

            if error.abs() < self.tol {
                break;
            }

            // Find free muscles (not at bounds) and redistribute
            let free_indices: Vec<usize> = activations
                .iter()
                .enumerate()
                .filter(|(_i, a)| **a > self.a_min + 1e-10 && **a < self.a_max - 1e-10)
                .map(|(i, _)| i)
                .collect();

            if free_indices.is_empty() {
                break;
            }

            let free_c2: f64 = free_indices.iter().map(|&i| coeffs[i] * coeffs[i]).sum();
            if free_c2 < 1e-20 {
                break;
            }

            let d_lambda = error / free_c2;
            for &i in &free_indices {
                activations[i] += d_lambda * coeffs[i];
            }
        }

        // Final clamp
        for a in activations.iter_mut() {
            *a = a.clamp(self.a_min, self.a_max);
        }

        activations
    }

    /// Compute the cost (sum of squared activations).
    pub fn cost(&self, activations: &[f64]) -> f64 {
        activations.iter().map(|a| a * a).sum()
    }

    /// Compute the torque produced by given activations.
    pub fn torque(&self, activations: &[f64]) -> f64 {
        activations
            .iter()
            .zip(self.f_max.iter())
            .zip(self.moment_arms.iter())
            .map(|((a, f), r)| a * f * r)
            .sum()
    }
}

// ---------------------------------------------------------------------------
// MusculotendonAnalysis
// ---------------------------------------------------------------------------

/// Analysis utilities for musculotendon systems.
///
/// Provides moment arm computation, joint torque calculation,
/// mechanical power, work, and fatigue index estimation.
#[derive(Debug, Clone)]
pub struct MusculotendonAnalysis {
    /// Collection of musculotendon units.
    pub units: Vec<MusculotendonUnit>,
    /// Moment arms for each unit about the joint (m).
    pub moment_arms: Vec<f64>,
    /// Accumulated work per unit (J).
    pub work: Vec<f64>,
    /// Fatigue time constant (s).
    pub fatigue_tau: f64,
    /// Fatigue state per unit \[0, 1\] (1 = fully fatigued).
    pub fatigue: Vec<f64>,
    /// Recovery rate for fatigue (1/s).
    pub recovery_rate: f64,
}

impl MusculotendonAnalysis {
    /// Create analysis for a set of musculotendon units.
    pub fn new(units: Vec<MusculotendonUnit>, moment_arms: Vec<f64>) -> Self {
        let n = units.len();
        Self {
            units,
            moment_arms,
            work: vec![0.0; n],
            fatigue_tau: 120.0,
            fatigue: vec![0.0; n],
            recovery_rate: 0.002,
        }
    }

    /// Number of muscles in the analysis.
    pub fn num_muscles(&self) -> usize {
        self.units.len()
    }

    /// Compute joint torque as sum of (force * moment_arm).
    pub fn joint_torque(&self) -> f64 {
        self.units
            .iter()
            .zip(self.moment_arms.iter())
            .map(|(u, r)| u.output_force() * r)
            .sum()
    }

    /// Compute power for each muscle (F * v).
    pub fn powers(&self) -> Vec<f64> {
        self.units
            .iter()
            .map(|u| u.output_force() * u.muscle.v_fiber)
            .collect()
    }

    /// Total mechanical power at the joint.
    pub fn joint_power(&self) -> f64 {
        self.units
            .iter()
            .zip(self.moment_arms.iter())
            .map(|(u, r)| u.output_force() * u.muscle.v_fiber * r.signum())
            .sum()
    }

    /// Step all units and accumulate work and fatigue.
    pub fn step(&mut self, dt: f64) {
        for (i, unit) in self.units.iter_mut().enumerate() {
            unit.step(dt);
            let power = unit.output_force() * unit.muscle.v_fiber;
            self.work[i] += power.abs() * dt;

            // Fatigue model: increases with activation, recovers when inactive
            let a = unit.activation.activation;
            let d_fatigue = (a * a) / self.fatigue_tau - self.recovery_rate * (1.0 - a);
            self.fatigue[i] += d_fatigue * dt;
            self.fatigue[i] = self.fatigue[i].clamp(0.0, 1.0);
        }
    }

    /// Fatigue index for each muscle \[0, 1\].
    pub fn fatigue_indices(&self) -> Vec<f64> {
        self.fatigue.clone()
    }

    /// Total accumulated work (J).
    pub fn total_work(&self) -> f64 {
        self.work.iter().sum()
    }

    /// Set excitation for a specific muscle.
    pub fn set_excitation(&mut self, index: usize, u: f64) {
        if index < self.units.len() {
            self.units[index].activation.set_excitation(u);
        }
    }

    /// Set musculotendon length for a specific muscle.
    pub fn set_mt_length(&mut self, index: usize, l_mt: f64) {
        if index < self.units.len() {
            self.units[index].set_length(l_mt);
        }
    }

    /// Get activations for all muscles.
    pub fn activations(&self) -> Vec<f64> {
        self.units.iter().map(|u| u.activation.activation).collect()
    }

    /// Get forces for all muscles.
    pub fn forces(&self) -> Vec<f64> {
        self.units.iter().map(|u| u.output_force()).collect()
    }

    /// Effective moment arm considering pennation.
    pub fn effective_moment_arm(&self, index: usize) -> f64 {
        if index >= self.units.len() {
            return 0.0;
        }
        let pen = self.units[index].pennation_angle();
        self.moment_arms[index] * pen.cos()
    }

    /// Reset accumulated work and fatigue.
    pub fn reset(&mut self) {
        for w in self.work.iter_mut() {
            *w = 0.0;
        }
        for f in self.fatigue.iter_mut() {
            *f = 0.0;
        }
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // 1. Force-velocity at isometric (v=0) should return 1.0
    #[test]
    fn test_fv_isometric() {
        let fv = ForceVelocityRelation::new();
        let f = fv.evaluate(0.0);
        assert!((f - 1.0).abs() < EPS, "Isometric FV should be 1.0, got {f}");
    }

    // 2. Force-velocity: concentric shortening reduces force
    #[test]
    fn test_fv_concentric() {
        let fv = ForceVelocityRelation::new();
        let f = fv.evaluate(-0.5);
        assert!(
            f > 0.0 && f < 1.0,
            "Concentric FV should be in (0,1), got {f}"
        );
    }

    // 3. Force-velocity: eccentric lengthening increases force
    #[test]
    fn test_fv_eccentric() {
        let fv = ForceVelocityRelation::new();
        let f = fv.evaluate(0.3);
        assert!(f > 1.0, "Eccentric FV should be > 1.0, got {f}");
    }

    // 4. Force-velocity: max eccentric force capped
    #[test]
    fn test_fv_eccentric_max() {
        let fv = ForceVelocityRelation::new();
        let f = fv.evaluate(100.0);
        assert!(
            (f - fv.f_eccentric_max).abs() < 0.01,
            "Large eccentric velocity should saturate at {}, got {f}",
            fv.f_eccentric_max
        );
    }

    // 5. Force-length: peak at optimal length (l_norm=1)
    #[test]
    fn test_fl_peak_at_optimal() {
        let fl = ForceLengthRelation::new();
        let f = fl.active(1.0);
        assert!(
            (f - 1.0).abs() < EPS,
            "Active FL at optimal should be 1.0, got {f}"
        );
    }

    // 6. Force-length: active force drops away from optimal
    #[test]
    fn test_fl_drops_away() {
        let fl = ForceLengthRelation::new();
        let f_short = fl.active(0.5);
        let f_long = fl.active(1.5);
        assert!(
            f_short < 1.0,
            "Active FL at 0.5 should be < 1.0, got {f_short}"
        );
        assert!(
            f_long < 1.0,
            "Active FL at 1.5 should be < 1.0, got {f_long}"
        );
    }

    // 7. Force-length: passive force zero at and below optimal
    #[test]
    fn test_fl_passive_zero_below_optimal() {
        let fl = ForceLengthRelation::new();
        assert!(fl.passive(0.8).abs() < EPS, "Passive at 0.8 should be 0");
        assert!(fl.passive(1.0).abs() < EPS, "Passive at 1.0 should be 0");
    }

    // 8. Force-length: passive force increases above optimal
    #[test]
    fn test_fl_passive_increases() {
        let fl = ForceLengthRelation::new();
        let f1 = fl.passive(1.3);
        let f2 = fl.passive(1.5);
        assert!(f1 > 0.0, "Passive at 1.3 should be > 0, got {f1}");
        assert!(f2 > f1, "Passive at 1.5 ({f2}) should be > at 1.3 ({f1})");
    }

    // 9. Tendon: no force below slack length
    #[test]
    fn test_tendon_slack() {
        let t = TendonModel::new(0.20, 1000.0);
        assert!(
            t.force(0.15).abs() < EPS,
            "Tendon force below slack should be 0"
        );
        assert!(
            t.force(0.20).abs() < EPS,
            "Tendon force at slack should be 0"
        );
    }

    // 10. Tendon: force increases above slack
    #[test]
    fn test_tendon_force_above_slack() {
        let t = TendonModel::new(0.20, 1000.0);
        let f1 = t.force(0.202);
        let f2 = t.force(0.205);
        assert!(f1 > 0.0, "Tendon force above slack should be > 0, got {f1}");
        assert!(f2 > f1, "More stretch should give more force: {f2} > {f1}");
    }

    // 11. Tendon stiffness positive above slack
    #[test]
    fn test_tendon_stiffness() {
        let t = TendonModel::new(0.20, 1000.0);
        let k = t.stiffness(0.205);
        assert!(
            k > 0.0,
            "Tendon stiffness above slack should be > 0, got {k}"
        );
    }

    // 12. Activation: steady state matches excitation
    #[test]
    fn test_activation_steady_state() {
        let mut ad = ActivationDynamics::new();
        ad.set_excitation(0.7);
        for _ in 0..10000 {
            ad.step(0.001);
        }
        assert!(
            (ad.activation - 0.7).abs() < 0.01,
            "Steady state activation should be ~0.7, got {}",
            ad.activation
        );
    }

    // 13. Activation: starts at minimum
    #[test]
    fn test_activation_initial() {
        let ad = ActivationDynamics::new();
        assert!(
            (ad.activation - 0.0).abs() < 0.02,
            "Initial activation should be ~0, got {}",
            ad.activation
        );
    }

    // 14. Activation: deactivation slower than activation
    #[test]
    fn test_activation_deactivation_slower() {
        let mut ad1 = ActivationDynamics::new();
        ad1.set_excitation(1.0);
        for _ in 0..500 {
            ad1.step(0.001);
        }
        let a_up = ad1.activation;

        let mut ad2 = ActivationDynamics::new();
        ad2.activation = 1.0;
        ad2.set_excitation(0.0);
        for _ in 0..500 {
            ad2.step(0.001);
        }
        let a_down = ad2.activation;

        // After same time, deactivation should leave more residual
        assert!(
            a_down > (1.0 - a_up),
            "Deactivation should be slower: a_down={a_down}, 1-a_up={}",
            1.0 - a_up
        );
    }

    // 15. Hill muscle: isometric force at optimal length
    #[test]
    fn test_hill_isometric_at_optimal() {
        let mut m = HillMuscleModel::new(1000.0, 0.10);
        m.set_state(0.10, 0.0);
        let f = m.active_force(1.0);
        assert!(
            (f - 1000.0).abs() < 1.0,
            "Isometric force at optimal should be ~F_max, got {f}"
        );
    }

    // 16. Hill muscle: zero activation gives zero active force
    #[test]
    fn test_hill_zero_activation() {
        let m = HillMuscleModel::new(1000.0, 0.10);
        let f = m.active_force(0.0);
        assert!(
            f.abs() < EPS,
            "Zero activation should give zero force, got {f}"
        );
    }

    // 17. Hill muscle: passive force at long length
    #[test]
    fn test_hill_passive_at_long_length() {
        let mut m = HillMuscleModel::new(1000.0, 0.10);
        m.set_state(0.15, 0.0);
        let f = m.passive_force();
        assert!(f > 0.0, "Passive force at 1.5*L_opt should be > 0, got {f}");
    }

    // 18. Musculotendon: equilibrium convergence
    #[test]
    fn test_mt_equilibrium() {
        let mut mt = MusculotendonUnit::new("test", 1000.0, 0.10, 0.20, 0.0);
        mt.activation.activation = 0.5;
        mt.l_mt = 0.30;
        mt.solve_equilibrium();

        let f_fiber = mt.projected_fiber_force();
        let f_tendon = mt.tendon_force();
        let error = (f_fiber - f_tendon).abs();
        assert!(
            error < 1.0,
            "Equilibrium should balance fiber/tendon forces: fiber={f_fiber}, tendon={f_tendon}"
        );
    }

    // 19. Musculotendon: pennation angle at optimal length
    #[test]
    fn test_mt_pennation() {
        let mt = MusculotendonUnit::new("test", 1000.0, 0.10, 0.20, 0.3);
        let pen = mt.pennation_angle();
        assert!(
            (pen - 0.3).abs() < 0.01,
            "Pennation at optimal should be ~0.3 rad, got {pen}"
        );
    }

    // 20. Musculotendon: zero activation gives minimal force
    #[test]
    fn test_mt_zero_activation() {
        let mut mt = MusculotendonUnit::new("test", 1000.0, 0.10, 0.20, 0.0);
        mt.activation.activation = 0.01;
        mt.l_mt = 0.30;
        mt.solve_equilibrium();
        let f = mt.output_force();
        // With minimal activation, force should be small
        assert!(
            f < 100.0,
            "Minimal activation should give small force, got {f}"
        );
    }

    // 21. Muscle redundancy: two muscles producing required torque
    #[test]
    fn test_redundancy_two_muscles() {
        let solver = MuscleRedundancy::new(vec![1000.0, 1000.0], vec![0.05, 0.05]);
        let activations = solver.solve(50.0);
        assert_eq!(activations.len(), 2);
        let tau = solver.torque(&activations);
        assert!((tau - 50.0).abs() < 1.0, "Torque should be ~50, got {tau}");
    }

    // 22. Muscle redundancy: symmetric muscles get equal activations
    #[test]
    fn test_redundancy_symmetric() {
        let solver = MuscleRedundancy::new(vec![500.0, 500.0, 500.0], vec![0.04, 0.04, 0.04]);
        let activations = solver.solve(30.0);
        let a0 = activations[0];
        for (i, &a) in activations.iter().enumerate() {
            assert!(
                (a - a0).abs() < 0.05,
                "Symmetric muscles should have equal activations: a[0]={a0}, a[{i}]={a}"
            );
        }
    }

    // 23. Muscle redundancy: zero torque gives minimum activation
    #[test]
    fn test_redundancy_zero_torque() {
        let solver = MuscleRedundancy::new(vec![1000.0, 1000.0], vec![0.05, 0.05]);
        let activations = solver.solve(0.0);
        let cost = solver.cost(&activations);
        // Activations should be at minimum
        assert!(
            cost < 0.01,
            "Zero torque should give minimal activations, cost={cost}"
        );
    }

    // 24. Analysis: joint torque computation
    #[test]
    fn test_analysis_joint_torque() {
        let u1 = MusculotendonUnit::new("m1", 1000.0, 0.10, 0.20, 0.0);
        let u2 = MusculotendonUnit::new("m2", 800.0, 0.08, 0.15, 0.0);
        let analysis = MusculotendonAnalysis::new(vec![u1, u2], vec![0.05, -0.03]);
        let _tau = analysis.joint_torque();
        // Just verify it computes without panic
        assert!(analysis.num_muscles() == 2);
    }

    // 25. Analysis: fatigue increases with sustained activation
    #[test]
    fn test_analysis_fatigue_increases() {
        let mut u1 = MusculotendonUnit::new("m1", 1000.0, 0.10, 0.20, 0.0);
        u1.activation.activation = 0.8;
        u1.activation.set_excitation(0.8);

        let mut analysis = MusculotendonAnalysis::new(vec![u1], vec![0.05]);
        analysis.set_excitation(0, 0.8);

        let f_before = analysis.fatigue[0];
        for _ in 0..1000 {
            analysis.step(0.01);
        }
        let f_after = analysis.fatigue[0];
        assert!(
            f_after > f_before,
            "Fatigue should increase: before={f_before}, after={f_after}"
        );
    }

    // 26. Force-velocity inverse consistency
    #[test]
    fn test_fv_inverse_consistency() {
        let fv = ForceVelocityRelation::new();
        let v_test = -0.3;
        let f = fv.evaluate(v_test);
        let v_recovered = fv.inverse(f);
        assert!(
            (v_recovered - v_test).abs() < 0.01,
            "FV inverse should recover velocity: expected {v_test}, got {v_recovered}"
        );
    }

    // 27. Force-length active range
    #[test]
    fn test_fl_active_range() {
        let fl = ForceLengthRelation::new();
        let (lo, hi) = fl.active_range(0.5);
        assert!(lo < 1.0, "Lower range should be < 1.0, got {lo}");
        assert!(hi > 1.0, "Upper range should be > 1.0, got {hi}");
        // Check symmetry (default skewness = 0)
        assert!(
            ((1.0 - lo) - (hi - 1.0)).abs() < 0.01,
            "Range should be symmetric: lo={lo}, hi={hi}"
        );
    }

    // 28. Tendon: is_slack check
    #[test]
    fn test_tendon_is_slack() {
        let t = TendonModel::new(0.20, 1000.0);
        assert!(t.is_slack(0.15), "Below slack should be slack");
        assert!(t.is_slack(0.20), "At slack should be slack");
        assert!(!t.is_slack(0.21), "Above slack should not be slack");
    }

    // 29. Activation: calcium dynamics
    #[test]
    fn test_activation_calcium() {
        let mut ad = ActivationDynamics::new();
        ad.set_excitation(1.0);
        for _ in 0..5000 {
            ad.step_calcium(0.001);
        }
        assert!(
            ad.calcium > 0.9,
            "Calcium should rise with excitation, got {}",
            ad.calcium
        );
        assert!(
            ad.activation > 0.8,
            "Activation should follow calcium, got {}",
            ad.activation
        );
    }

    // 30. Musculotendon: step updates activation
    #[test]
    fn test_mt_step_updates_activation() {
        let mut mt = MusculotendonUnit::new("test", 1000.0, 0.10, 0.20, 0.0);
        mt.activation.set_excitation(0.9);
        let a_before = mt.activation.activation;
        for _ in 0..100 {
            mt.step(0.001);
        }
        let a_after = mt.activation.activation;
        assert!(
            a_after > a_before,
            "Activation should increase: {a_before} -> {a_after}"
        );
    }

    // 31. Hill muscle: stiffness is positive at optimal length
    #[test]
    fn test_hill_stiffness_positive() {
        let mut m = HillMuscleModel::new(1000.0, 0.10);
        m.set_state(0.10, 0.0);
        let k = m.stiffness(0.5);
        assert!(k > 0.0, "Muscle stiffness should be positive, got {k}");
    }

    // 32. Tendon: length_at_force_fraction round-trip
    #[test]
    fn test_tendon_length_force_roundtrip() {
        let t = TendonModel::new(0.20, 1000.0);
        let frac = 0.5;
        let l = t.length_at_force_fraction(frac);
        let f = t.force(l);
        let f_expected = frac * t.f_max;
        assert!(
            (f - f_expected).abs() < 5.0,
            "Round-trip: expected ~{f_expected}, got {f}"
        );
    }

    // 33. Analysis: effective moment arm with pennation
    #[test]
    fn test_effective_moment_arm() {
        let u = MusculotendonUnit::new("test", 1000.0, 0.10, 0.20, 0.5);
        let analysis = MusculotendonAnalysis::new(vec![u], vec![0.05]);
        let ema = analysis.effective_moment_arm(0);
        assert!(
            ema < 0.05,
            "Effective moment arm with pennation should be < raw: {ema}"
        );
        assert!(ema > 0.0, "Effective moment arm should be positive: {ema}");
    }

    // 34. Activation: reset clears state
    #[test]
    fn test_activation_reset() {
        let mut ad = ActivationDynamics::new();
        ad.set_excitation(1.0);
        for _ in 0..1000 {
            ad.step(0.001);
        }
        ad.reset();
        assert!(
            ad.activation <= ad.a_min + EPS,
            "Reset should clear activation, got {}",
            ad.activation
        );
    }

    // 35. Redundancy: larger muscle contributes more torque fraction
    #[test]
    fn test_redundancy_larger_muscle_more_torque() {
        let solver = MuscleRedundancy::new(vec![2000.0, 500.0], vec![0.05, 0.05]);
        let activations = solver.solve(50.0);
        let torque0 = activations[0] * 2000.0 * 0.05;
        let torque1 = activations[1] * 500.0 * 0.05;
        assert!(
            torque0 > torque1,
            "Larger muscle should contribute more torque: t0={torque0}, t1={torque1}"
        );
    }
}
