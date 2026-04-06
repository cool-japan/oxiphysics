//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// First-order muscle activation dynamics model.
///
/// Models the calcium transient and cross-bridge kinetics that translate
/// neural drive (0–1) to mechanical activation (0–1).
#[derive(Debug, Clone)]
pub struct ActivationDynamics {
    /// Activation time constant τ_act (s), typically ~10 ms.
    pub tau_act: f64,
    /// Deactivation time constant τ_deact (s), typically ~40 ms.
    pub tau_deact: f64,
    /// Minimum activation level (dimensionless, ~0.01).
    pub activation_min: f64,
    /// Current activation a ∈ \[0,1\].
    pub activation: f64,
    /// Nonlinear activation-neural drive shape factor b (default 0.1).
    pub shape_factor: f64,
}
impl ActivationDynamics {
    /// Default parameters (Zajac 1989 / Thelen 2003).
    pub fn default_params() -> Self {
        Self {
            tau_act: 0.010,
            tau_deact: 0.040,
            activation_min: 0.01,
            activation: 0.0,
            shape_factor: 0.1,
        }
    }
    /// Fast-twitch fiber dynamics.
    pub fn fast_twitch() -> Self {
        Self {
            tau_act: 0.005,
            tau_deact: 0.020,
            activation_min: 0.01,
            activation: 0.0,
            shape_factor: 0.1,
        }
    }
    /// Advance activation by one time step.
    ///
    /// `u` is the neural excitation (0–1).
    /// `dt` is the time step (s).
    pub fn step(&mut self, u: f64, dt: f64) {
        let u_clamped = clamp(u, 0.0, 1.0);
        let tau = if u_clamped >= self.activation {
            self.tau_act * (0.5 + 1.5 * self.activation)
        } else {
            self.tau_deact / (0.5 + 1.5 * self.activation)
        };
        let da = (u_clamped - self.activation) / tau;
        self.activation = clamp(self.activation + da * dt, self.activation_min, 1.0);
    }
    /// Steady-state activation for a given neural excitation.
    pub fn steady_state(&self, u: f64) -> f64 {
        clamp(u, self.activation_min, 1.0)
    }
    /// Rate of change of activation dA/dt (s⁻¹) at current state.
    pub fn da_dt(&self, u: f64) -> f64 {
        let u_clamped = clamp(u, 0.0, 1.0);
        let tau = if u_clamped >= self.activation {
            self.tau_act * (0.5 + 1.5 * self.activation)
        } else {
            self.tau_deact / (0.5 + 1.5 * self.activation)
        };
        (u_clamped - self.activation) / tau
    }
}
/// Motor unit pool with Henneman's size-principle recruitment.
///
/// Motor units are ordered by recruitment threshold (smallest first).
/// As neural drive increases, progressively larger, faster units are recruited.
#[derive(Debug, Clone)]
pub struct MuscleRecruitment {
    /// Ordered list of motor units (ascending recruitment threshold).
    pub motor_units: Vec<MotorUnit>,
    /// Current neural drive (0–1).
    pub neural_drive: f64,
    /// Maximum pool force (sum of all peak forces, N).
    pub total_peak_force: f64,
}
impl MuscleRecruitment {
    /// Create a motor pool with `n` motor units spanning the given max force.
    ///
    /// Applies Henneman's principle: thresholds spaced logarithmically,
    /// twitch forces increase exponentially with threshold.
    pub fn new(n: usize, max_force: f64) -> Self {
        let mut units = Vec::with_capacity(n);
        let log_range = 3.0_f64;
        for i in 0..n {
            let t = i as f64 / (n as f64 - 1.0).max(1.0);
            let threshold = 10.0_f64.powf(-log_range + log_range * t)
                / 10.0_f64.powf(-log_range + log_range * 0.0);
            let threshold_norm = clamp(t, 0.0, 1.0);
            let twitch_force = max_force / (n as f64) * 10.0_f64.powf(3.0 * t);
            let twitch_force_norm = twitch_force / (max_force / (n as f64) * 10.0_f64.powf(3.0));
            let contraction_time = if t > 0.5 { 0.030 } else { 0.080 };
            let is_fast = t > 0.5;
            let _ = threshold;
            units.push(MotorUnit::new(
                threshold_norm,
                twitch_force_norm * max_force,
                contraction_time,
                is_fast,
            ));
        }
        let total_peak = units.iter().map(|u| u.peak_force()).sum();
        Self {
            motor_units: units,
            neural_drive: 0.0,
            total_peak_force: total_peak,
        }
    }
    /// Update firing rates and forces for all units given a neural drive `u` ∈ \[0,1\].
    pub fn update(&mut self, u: f64) {
        self.neural_drive = clamp(u, 0.0, 1.0);
        for mu in &mut self.motor_units {
            if u >= mu.recruitment_threshold {
                let excess =
                    (u - mu.recruitment_threshold) / (1.0 - mu.recruitment_threshold).max(1e-6);
                let ff = clamp(excess, 0.0, 1.0);
                mu.firing_rate = mu.required_firing_rate(0.1 + 0.89 * ff).max(0.0);
                let sigmoid = 1.0 / (1.0 + (-12.0 * (ff - 0.5)).exp());
                mu.force = mu.peak_force() * sigmoid;
            } else {
                mu.firing_rate = 0.0;
                mu.force = 0.0;
            }
        }
    }
    /// Total force produced by the pool (N).
    pub fn total_force(&self) -> f64 {
        self.motor_units.iter().map(|u| u.force).sum()
    }
    /// Number of recruited motor units.
    pub fn recruited_count(&self) -> usize {
        if self.neural_drive <= 0.0 {
            return 0;
        }
        self.motor_units
            .iter()
            .filter(|u| u.recruitment_threshold <= self.neural_drive)
            .count()
    }
    /// Fraction of the maximum force being produced (0–1).
    pub fn force_fraction(&self) -> f64 {
        if self.total_peak_force < 1e-10 {
            return 0.0;
        }
        self.total_force() / self.total_peak_force
    }
}
/// Muscle optimization routines: static (minimise muscle stress²) and
/// EMG-driven (activation tracks EMG envelope).
#[derive(Debug, Clone)]
pub struct MuscleOptimization {
    /// Muscles in the optimization problem.
    pub muscles: Vec<OptMuscle>,
}
impl MuscleOptimization {
    /// Create a new optimization problem with the given muscles.
    pub fn new(muscles: Vec<OptMuscle>) -> Self {
        Self { muscles }
    }
    /// Static optimization: find activations that produce `required_torque`
    /// (N·m) while minimising the sum of squared muscle activations.
    ///
    /// Uses a greedy proportional allocation scaled by maximum torque.
    pub fn static_optimization(&mut self, required_torque: f64) -> Vec<f64> {
        let n = self.muscles.len();
        if n == 0 {
            return Vec::new();
        }
        let total_max: f64 = self.muscles.iter().map(|m| m.max_torque()).sum();
        if total_max < 1e-12 {
            return vec![0.0; n];
        }
        let fraction = (required_torque / total_max).clamp(0.0, 1.0);
        let activations: Vec<f64> = self
            .muscles
            .iter()
            .map(|m| {
                if m.max_torque() > 1e-12 {
                    (fraction * total_max / (m.max_torque() * n as f64)).clamp(0.0, 1.0)
                } else {
                    0.0
                }
            })
            .collect();
        for (m, &a) in self.muscles.iter_mut().zip(activations.iter()) {
            m.activation = a;
        }
        activations
    }
    /// EMG-driven step: advance activation dynamics driven by an EMG signal.
    ///
    /// `emg` ∈ \[0,1\] is the normalised EMG envelope for each muscle.
    /// Returns the current activation values after the step.
    pub fn emg_driven_step(&mut self, dt: f64, emg: &[f64]) -> Vec<f64> {
        let n = self.muscles.len();
        let mut activations = Vec::with_capacity(n);
        for (i, m) in self.muscles.iter_mut().enumerate() {
            let emg_val = if i < emg.len() {
                emg[i].clamp(0.0, 1.0)
            } else {
                0.0
            };
            let tau_act = 0.01;
            let tau_deact = 0.04;
            let tau = if emg_val > m.activation {
                tau_act
            } else {
                tau_deact
            };
            let da = (emg_val - m.activation) / tau * dt;
            m.activation = clamp(m.activation + da, 0.01, 1.0);
            activations.push(m.activation);
        }
        activations
    }
    /// Compute sum of squared activations (objective function for static opt).
    pub fn stress_objective(&self) -> f64 {
        self.muscles
            .iter()
            .map(|m| m.activation * m.activation)
            .sum()
    }
}
/// Hill force–velocity relationship for skeletal muscle.
///
/// Concentric (shortening) phase follows Hill's hyperbolic equation.
/// Eccentric (lengthening) phase uses an extended model (Thelen 2003).
#[derive(Debug, Clone)]
pub struct ForceVelocityRelation {
    /// Maximum shortening velocity (m/s, positive value).
    pub v_max: f64,
    /// Hill's curve constant a/F0 (dimensionless, typically 0.25).
    pub a_over_f0: f64,
    /// Eccentric force maximum factor (typically 1.5).
    pub eccentric_force_max: f64,
    /// Eccentric curve shape parameter (0–1, typically 0.1).
    pub eccentric_shape: f64,
}
impl ForceVelocityRelation {
    /// Default fast-twitch fiber parameters.
    pub fn fast_twitch() -> Self {
        Self {
            v_max: 10.0,
            a_over_f0: 0.25,
            eccentric_force_max: 1.5,
            eccentric_shape: 0.1,
        }
    }
    /// Default slow-twitch fiber parameters.
    pub fn slow_twitch() -> Self {
        Self {
            v_max: 4.0,
            a_over_f0: 0.25,
            eccentric_force_max: 1.4,
            eccentric_shape: 0.08,
        }
    }
    /// Normalised force–velocity multiplier fv(v̄).
    ///
    /// `v_norm` = v / v_max (negative = shortening, positive = lengthening).
    pub fn fv(&self, v_norm: f64) -> f64 {
        if v_norm <= 0.0 {
            let a = self.a_over_f0;
            let vn = (-v_norm).min(1.0);
            (1.0 - vn) / (1.0 + vn / a)
        } else {
            let vn = v_norm.min(1.0);
            let fmax = self.eccentric_force_max;
            let b = self.eccentric_shape;
            fmax - (fmax - 1.0) * (1.0 - vn) / (1.0 + vn / b)
        }
    }
    /// Optimal (maximum power) shortening velocity normalised to v_max.
    ///
    /// Solved from d(F*v)/dv = 0: v_opt/v_max = (√(a+1) − √a) / √(a+1).
    pub fn optimal_velocity_fraction(&self) -> f64 {
        let a = self.a_over_f0;
        let sqrt_term = (a + 1.0).sqrt();
        (sqrt_term - a.sqrt()) / sqrt_term
    }
}
/// Geometric description of a muscle: fiber direction, pennation angle,
/// volume, and attachment points.
///
/// Pennation angle `alpha` (rad) is the angle between the muscle fiber
/// direction and the muscle line of action.
#[derive(Debug, Clone)]
pub struct MuscleGeometry {
    /// Proximal (origin) attachment point (m).
    pub origin: [f64; 3],
    /// Distal (insertion) attachment point (m).
    pub insertion: [f64; 3],
    /// Muscle fiber direction unit vector (may differ from line-of-action due to pennation).
    pub fiber_direction: [f64; 3],
    /// Pennation angle α at optimal fiber length (rad).
    pub pennation_angle: f64,
    /// Optimal fiber length l_opt (m).
    pub optimal_fiber_length: f64,
    /// Physiological cross-sectional area PCSA (m²).
    pub pcsa: f64,
    /// Muscle volume (m³).
    pub volume: f64,
}
impl MuscleGeometry {
    /// Create a muscle geometry from origin/insertion points and PCSA.
    pub fn new(
        origin: [f64; 3],
        insertion: [f64; 3],
        pennation_angle: f64,
        optimal_fiber_length: f64,
        pcsa: f64,
    ) -> Self {
        let dir = normalize3(sub3(insertion, origin));
        let volume = pcsa * optimal_fiber_length;
        Self {
            origin,
            insertion,
            fiber_direction: dir,
            pennation_angle,
            optimal_fiber_length,
            pcsa,
            volume,
        }
    }
    /// Typical soleus muscle geometry.
    pub fn soleus() -> Self {
        Self::new([0.0, -0.20, 0.0], [0.0, -0.45, 0.0], 0.349, 0.030, 25.0e-4)
    }
    /// Typical rectus femoris muscle geometry.
    pub fn rectus_femoris() -> Self {
        Self::new([0.0, 0.10, 0.0], [0.0, -0.30, 0.0], 0.087, 0.084, 16.0e-4)
    }
    /// Musculotendon unit length (origin-to-insertion straight line, m).
    pub fn musculotendon_length(&self) -> f64 {
        norm3(sub3(self.insertion, self.origin))
    }
    /// Fiber length accounting for pennation angle given musculotendon length
    /// and tendon slack length `l_ts` (m).
    pub fn fiber_length(&self, musculotendon_len: f64, tendon_slack_length: f64) -> f64 {
        let l_mt = musculotendon_len;
        let l_ts = tendon_slack_length;
        let cos_alpha = self.pennation_angle.cos();
        if cos_alpha < 1e-10 {
            return self.optimal_fiber_length;
        }
        (l_mt - l_ts).max(1e-6) / cos_alpha
    }
    /// Projected fiber component along the line of action.
    pub fn fiber_force_projection(&self, fiber_force: f64) -> f64 {
        fiber_force * self.pennation_angle.cos()
    }
}
/// Active and passive force–length relationship for skeletal muscle fibers.
///
/// Active curve: bell-shaped Gaussian around the optimal fiber length.
/// Passive curve: exponential toe-region for stretch beyond optimal.
/// Titin: a third (passive) element providing force at long sarcomere lengths.
#[derive(Debug, Clone)]
pub struct ForceLengthRelation {
    /// Optimal fiber length l_opt (m).
    pub optimal_fiber_length: f64,
    /// Width of the active force–length curve (σ in Gaussian, fraction of l_opt).
    pub active_width: f64,
    /// Passive stiffness coefficient k_p (dimensionless).
    pub passive_stiffness: f64,
    /// Passive toe-region strain ε₀ (fraction of l_opt, stretch onset).
    pub passive_toe_strain: f64,
    /// Titin stiffness (N/m).
    pub titin_stiffness: f64,
    /// Titin engagement length (fraction of l_opt, typically 1.4).
    pub titin_engagement: f64,
}
impl ForceLengthRelation {
    /// Typical gastrocnemius parameters.
    pub fn gastrocnemius() -> Self {
        Self {
            optimal_fiber_length: 0.051,
            active_width: 0.56,
            passive_stiffness: 3.0,
            passive_toe_strain: 0.1,
            titin_stiffness: 10.0,
            titin_engagement: 1.4,
        }
    }
    /// Default generic parameters.
    pub fn default_params() -> Self {
        Self {
            optimal_fiber_length: 0.10,
            active_width: 0.56,
            passive_stiffness: 3.0,
            passive_toe_strain: 0.10,
            titin_stiffness: 8.0,
            titin_engagement: 1.4,
        }
    }
    /// Normalised active force fl(l̄).
    ///
    /// `l_norm` = l_fiber / l_opt.
    pub fn fl_active(&self, l_norm: f64) -> f64 {
        let w = self.active_width;
        (-(l_norm - 1.0).powi(2) / (2.0 * w * w)).exp()
    }
    /// Normalised passive force fp(l̄).
    ///
    /// `l_norm` = l_fiber / l_opt (>1 when stretched beyond optimal).
    pub fn fl_passive(&self, l_norm: f64) -> f64 {
        let eps0 = self.passive_toe_strain;
        if l_norm <= 1.0 + eps0 {
            0.0
        } else {
            let e = l_norm - 1.0 - eps0;
            self.passive_stiffness * e * e
        }
    }
    /// Titin force (N) at a given fiber length `l_fiber` (m).
    pub fn titin_force(&self, l_fiber: f64) -> f64 {
        let l_eng = self.titin_engagement * self.optimal_fiber_length;
        if l_fiber <= l_eng {
            0.0
        } else {
            self.titin_stiffness * (l_fiber - l_eng)
        }
    }
    /// Total passive force multiplier at normalised length `l_norm`.
    pub fn fl_total_passive(&self, l_norm: f64) -> f64 {
        let titin_l = l_norm * self.optimal_fiber_length;
        self.fl_passive(l_norm)
            + self.titin_force(titin_l) / (self.titin_stiffness * self.optimal_fiber_length)
    }
}
/// Multi-muscle whole-body model with redundancy resolution via static
/// optimisation (minimise sum of squared activations subject to joint-moment
/// target, a standard approach in musculoskeletal modelling).
#[derive(Debug, Clone)]
pub struct WholeBodyMuscle {
    /// List of muscles.
    pub muscles: Vec<MuscleEntry>,
}
impl WholeBodyMuscle {
    /// Create a new whole-body model.
    pub fn new(muscles: Vec<MuscleEntry>) -> Self {
        Self { muscles }
    }
    /// Total joint moment produced by all muscles (N·m).
    pub fn total_joint_moment(&self) -> f64 {
        self.muscles.iter().map(|m| m.joint_moment()).sum()
    }
    /// Maximum possible joint moment (N·m) when all muscles at full activation.
    pub fn max_joint_moment(&self) -> f64 {
        self.muscles
            .iter()
            .map(|m| m.mtu.muscle.f_max * m.moment_arm)
            .sum()
    }
    /// Static optimisation: compute activations that produce the target moment
    /// while minimising the sum of squared activations (load sharing criterion).
    ///
    /// Uses a simple proportional allocation (each muscle's contribution
    /// is proportional to its capacity F0·r), which is the analytical solution
    /// to the quadratic cost with equal weighting.
    ///
    /// Returns a `Vec`f64` of activations, one per muscle.
    pub fn static_optimisation(&self, target_moment: f64) -> Vec<f64> {
        let n = self.muscles.len();
        if n == 0 {
            return Vec::new();
        }
        let capacities: Vec<f64> = self
            .muscles
            .iter()
            .map(|m| m.mtu.muscle.f_max * m.moment_arm)
            .collect();
        let sum_c2: f64 = capacities.iter().map(|c| c * c).sum();
        if sum_c2 < 1e-20 {
            return vec![0.0; n];
        }
        let lambda = target_moment / sum_c2;
        capacities
            .iter()
            .zip(self.muscles.iter())
            .map(|(c, m)| clamp(c * lambda, 0.0, m.max_activation))
            .collect()
    }
    /// Redundancy index: number of muscles − degrees of freedom (here 1 joint = 1 DOF).
    pub fn redundancy_index(&self) -> usize {
        self.muscles.len().saturating_sub(1)
    }
    /// Step all MTUs forward by `dt` (s) using target activations.
    pub fn step(&mut self, dt: f64, activations: &[f64], l_mts: &[f64]) {
        for ((m, &a), &l) in self
            .muscles
            .iter_mut()
            .zip(activations.iter())
            .zip(l_mts.iter())
        {
            m.mtu.step(dt, l, a);
        }
    }
}
/// Musculotendon unit combining fiber mechanics, tendon compliance,
/// pennation kinematics, and geometric path integration.
#[derive(Debug, Clone)]
pub struct MusculotendonUnit {
    /// Muscle model.
    pub muscle: HillMuscleModel,
    /// Geometry.
    pub geometry: MuscleGeometry,
    /// Musculotendon length (m).
    pub lmt: f64,
    /// Musculotendon velocity (m/s, negative = shortening).
    pub v_mt: f64,
}
impl MusculotendonUnit {
    /// Create from model and geometry with default length.
    pub fn new(muscle: HillMuscleModel, geometry: MuscleGeometry) -> Self {
        let lmt = geometry.musculotendon_length();
        Self {
            muscle,
            geometry,
            lmt,
            v_mt: 0.0,
        }
    }
    /// Typical soleus MTU.
    pub fn soleus() -> Self {
        let muscle = HillMuscleModel::soleus();
        let geometry = MuscleGeometry::soleus();
        Self::new(muscle, geometry)
    }
    /// Tendon strain (dimensionless).
    pub fn tendon_strain(&self) -> f64 {
        let l_ts = self.muscle.tendon_slack_length;
        let cos_a = self.muscle.pennation_angle.cos();
        let l_ce = self.muscle.l_ce_norm * self.muscle.l_opt;
        let l_t = self.lmt - l_ce * cos_a;
        if l_t <= l_ts {
            0.0
        } else {
            (l_t - l_ts) / l_ts
        }
    }
    /// Fiber length (m).
    pub fn fiber_length(&self) -> f64 {
        self.muscle.l_ce_norm * self.muscle.l_opt
    }
    /// Pennation angle at current fiber length (rad) — constant volume.
    ///
    /// sin(α) = l_opt * sin(α₀) / l_ce  (constant muscle width assumption)
    pub fn current_pennation_angle(&self) -> f64 {
        let sin_a0 = self.muscle.pennation_angle.sin();
        let l_ce = self.fiber_length();
        let sin_a = (self.muscle.l_opt * sin_a0 / l_ce).min(1.0);
        sin_a.asin()
    }
    /// Musculotendon force projected along the line of action (N).
    pub fn mtu_force(&self) -> f64 {
        let cos_a = self.current_pennation_angle().cos();
        self.muscle.ce_force_total() * cos_a
    }
    /// Step the MTU forward by `dt` (s) given new length `l_mt` (m) and activation target.
    pub fn step(&mut self, dt: f64, l_mt: f64, a_target: f64) {
        self.v_mt = (l_mt - self.lmt) / dt;
        self.lmt = l_mt;
        self.muscle.step(dt, l_mt, a_target);
    }
}
/// Single motor unit within the motor pool.
///
/// Characterised by a recruitment threshold and maximum twitch force.
#[derive(Debug, Clone)]
pub struct MotorUnit {
    /// Recruitment threshold (normalised neural drive, 0–1).
    pub recruitment_threshold: f64,
    /// Maximum twitch force (N).
    pub twitch_force: f64,
    /// Twitch contraction time (s).
    pub contraction_time: f64,
    /// Fiber type: `true` = fast-twitch (type II), `false` = slow-twitch (type I).
    pub is_fast_twitch: bool,
    /// Current firing rate (Hz).
    pub firing_rate: f64,
    /// Current force output (N).
    pub force: f64,
}
impl MotorUnit {
    /// Construct a new motor unit.
    pub fn new(
        recruitment_threshold: f64,
        twitch_force: f64,
        contraction_time: f64,
        is_fast_twitch: bool,
    ) -> Self {
        Self {
            recruitment_threshold,
            twitch_force,
            contraction_time,
            is_fast_twitch,
            firing_rate: 0.0,
            force: 0.0,
        }
    }
    /// Maximum force at saturation firing rate (N).
    pub fn peak_force(&self) -> f64 {
        self.twitch_force * 7.0
    }
    /// Firing rate needed to produce a given force fraction (0–1).
    ///
    /// Based on a sigmoidal rate–tension curve centred at the characteristic
    /// firing rate (1/(2 * contraction_time)).
    pub fn required_firing_rate(&self, force_fraction: f64) -> f64 {
        let ff = clamp(force_fraction, 0.01, 0.99);
        let f_char = 1.0 / (2.0 * self.contraction_time);
        f_char * (ff / (1.0 - ff)).ln() + f_char
    }
}
/// Entry in the musculoskeletal model — a single muscle with its moment arm.
#[derive(Debug, Clone)]
pub struct MuscleEntry {
    /// Muscle name (label).
    pub name: String,
    /// Hill muscle model.
    pub mtu: MusculotendonUnit,
    /// Moment arm about the joint of interest (m).
    pub moment_arm: f64,
    /// Maximum activation (0–1), can be used for voluntary effort scaling.
    pub max_activation: f64,
}
impl MuscleEntry {
    /// Create a muscle entry.
    pub fn new(name: &str, mtu: MusculotendonUnit, moment_arm: f64, max_activation: f64) -> Self {
        Self {
            name: name.to_string(),
            mtu,
            moment_arm,
            max_activation,
        }
    }
    /// Moment produced by this muscle about the joint (N·m).
    pub fn joint_moment(&self) -> f64 {
        self.mtu.mtu_force() * self.moment_arm
    }
}
/// Configuration of a skeletal joint including angle limits and position.
#[derive(Debug, Clone)]
pub struct JointConfig {
    /// Name/label of the joint (e.g. "knee").
    pub name: String,
    /// Joint centre position in body-fixed frame (m).
    pub position: [f64; 3],
    /// Current joint angle (rad).
    pub angle: f64,
    /// Minimum joint angle (rad).
    pub min_angle: f64,
    /// Maximum joint angle (rad).
    pub max_angle: f64,
}
impl JointConfig {
    /// Create a new joint configuration.
    pub fn new(name: &str, position: [f64; 3], angle: f64, min_angle: f64, max_angle: f64) -> Self {
        Self {
            name: name.to_string(),
            position,
            angle: clamp(angle, min_angle, max_angle),
            min_angle,
            max_angle,
        }
    }
    /// Set the joint angle, clamped to [min_angle, max_angle].
    pub fn set_angle(&mut self, angle: f64) {
        self.angle = clamp(angle, self.min_angle, self.max_angle);
    }
}
/// A muscle attachment specifying origin, insertion, joint it crosses, and moment arm.
#[derive(Debug, Clone)]
pub struct MuscleAttachment {
    /// Muscle name.
    pub name: String,
    /// Musculotendon unit for this attachment.
    pub mtu: MusculotendonUnit,
    /// Index of the joint this muscle crosses in the parent system.
    pub joint_idx: usize,
    /// Moment arm length (m) at nominal configuration.
    pub moment_arm: f64,
    /// Proximal (origin) point in body frame (m).
    pub origin: [f64; 3],
    /// Distal (insertion) point in body frame (m).
    pub insertion: [f64; 3],
}
impl MuscleAttachment {
    /// Create a new muscle attachment.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        name: &str,
        mtu: MusculotendonUnit,
        joint_idx: usize,
        moment_arm: f64,
        origin: [f64; 3],
        insertion: [f64; 3],
    ) -> Self {
        Self {
            name: name.to_string(),
            mtu,
            joint_idx,
            moment_arm: moment_arm.abs(),
            origin,
            insertion,
        }
    }
    /// Muscle line-of-action vector (unit, from origin to insertion).
    pub fn line_of_action(&self) -> [f64; 3] {
        normalize3(sub3(self.insertion, self.origin))
    }
}
/// Musculoskeletal system: collection of joints and muscles with force paths.
#[derive(Debug, Clone)]
pub struct MusculoskeletalSystem {
    /// Joints in the system (indexed from 0).
    pub joints: Vec<JointConfig>,
    /// Muscles crossing joints.
    pub muscles: Vec<MuscleAttachment>,
}
impl MusculoskeletalSystem {
    /// Create an empty musculoskeletal system.
    pub fn new() -> Self {
        Self {
            joints: Vec::new(),
            muscles: Vec::new(),
        }
    }
    /// Add a joint to the system and return its index.
    pub fn add_joint(&mut self, cfg: JointConfig) -> usize {
        let idx = self.joints.len();
        self.joints.push(cfg);
        idx
    }
    /// Add a muscle attachment to the system.
    pub fn add_muscle(&mut self, att: MuscleAttachment) {
        self.muscles.push(att);
    }
    /// Compute joint torques (N·m) from all muscle forces.
    ///
    /// Returns a `Vec`f64` with one torque value per joint, summing
    /// contributions from all muscles crossing that joint.
    pub fn compute_joint_torques(&self) -> Vec<f64> {
        let mut torques = vec![0.0f64; self.joints.len()];
        for att in &self.muscles {
            if att.joint_idx >= self.joints.len() {
                continue;
            }
            let muscle = &att.mtu.muscle;
            let f_active = muscle.ce_force_active();
            let f_passive = muscle.pee_force();
            let f_total = f_active + f_passive;
            torques[att.joint_idx] += f_total * att.moment_arm;
        }
        torques
    }
    /// Advance all muscle MTUs by `dt` seconds with given activation targets.
    pub fn step(&mut self, dt: f64, activations: &[f64], lmt_values: &[f64]) {
        for (i, att) in self.muscles.iter_mut().enumerate() {
            let a = if i < activations.len() {
                activations[i]
            } else {
                0.0
            };
            let lmt = if i < lmt_values.len() {
                lmt_values[i]
            } else {
                att.mtu.lmt
            };
            att.mtu.step(dt, lmt, a);
        }
    }
}
/// Dual-component fatigue model with central nervous system and peripheral
/// (metabolic) fatigue.
///
/// Based on the three-compartment fatigue model (Ma et al. 2009):
/// active, fatigued, and resting motor unit pools.
#[derive(Debug, Clone)]
pub struct FatigueDynamics {
    /// Central fatigue level ∈ \[0, 1\] (failure of motor drive).
    pub central_fatigue: f64,
    /// Peripheral fatigue level ∈ \[0, 1\] (metabolic/neuromuscular junction).
    pub peripheral_fatigue: f64,
    /// Central fatigue accumulation rate constant (1/s).
    pub central_fatigue_rate: f64,
    /// Peripheral fatigue accumulation rate constant (1/s).
    pub peripheral_fatigue_rate: f64,
    /// Central recovery rate constant (1/s).
    pub central_recovery_rate: f64,
    /// Peripheral recovery rate constant (1/s).
    pub peripheral_recovery_rate: f64,
    /// Effort threshold above which fatigue accumulates (0-1).
    pub effort_threshold: f64,
}
impl FatigueDynamics {
    /// Create a fresh (unfatigued) model with typical parameter values.
    pub fn new() -> Self {
        Self {
            central_fatigue: 0.0,
            peripheral_fatigue: 0.0,
            central_fatigue_rate: 0.005,
            peripheral_fatigue_rate: 0.015,
            central_recovery_rate: 0.002,
            peripheral_recovery_rate: 0.004,
            effort_threshold: 0.2,
        }
    }
    /// Advance fatigue/recovery dynamics for time `dt` at `effort` level \[0,1\].
    pub fn step(&mut self, dt: f64, effort: f64) {
        let e = effort.clamp(0.0, 1.0);
        if e > self.effort_threshold {
            let excess = e - self.effort_threshold;
            self.central_fatigue += self.central_fatigue_rate * excess * dt;
            self.peripheral_fatigue += self.peripheral_fatigue_rate * excess * dt;
        } else {
            self.central_fatigue -= self.central_recovery_rate * dt;
            self.peripheral_fatigue -= self.peripheral_recovery_rate * dt;
        }
        self.central_fatigue = clamp(self.central_fatigue, 0.0, 1.0);
        self.peripheral_fatigue = clamp(self.peripheral_fatigue, 0.0, 1.0);
    }
    /// Combined fatigue index ∈ \[0, 1\]: weighted mean of both components.
    ///
    /// Central fatigue (30%) and peripheral fatigue (70%) contribution.
    pub fn fatigue_index(&self) -> f64 {
        0.3 * self.central_fatigue + 0.7 * self.peripheral_fatigue
    }
    /// Effective force scale accounting for fatigue (1.0 = fresh, 0.0 = fully fatigued).
    ///
    /// `voluntary_effort` ∈ \[0,1\] is the intended activation level.
    pub fn effective_force_scale(&self, voluntary_effort: f64) -> f64 {
        let e = voluntary_effort.clamp(0.0, 1.0);
        let fi = self.fatigue_index();
        (e * (1.0 - fi)).max(0.0)
    }
    /// Time to recovery (seconds) assuming rest (zero effort).
    ///
    /// Estimates time until `fatigue_index() < threshold`.
    pub fn time_to_recovery(&self, threshold: f64) -> f64 {
        let fi = self.fatigue_index();
        if fi <= threshold {
            return 0.0;
        }
        let decay_rate = 0.3 * self.central_recovery_rate + 0.7 * self.peripheral_recovery_rate;
        if decay_rate < 1e-12 {
            return f64::INFINITY;
        }
        -(fi - threshold).ln() / decay_rate
    }
}
/// Three-element Hill muscle model: contractile element (CE), series elastic
/// element (SEE), and parallel elastic element (PEE).
///
/// State: activation `a` ∈ \[0,1\] and normalised fiber length `l_ce`.
#[derive(Debug, Clone)]
pub struct HillMuscleModel {
    /// Maximum isometric force F₀ (N).
    pub f_max: f64,
    /// Optimal fiber length l_opt (m).
    pub l_opt: f64,
    /// Tendon slack length l_ts (m).
    pub tendon_slack_length: f64,
    /// Tendon stiffness k_t (N/m).
    pub tendon_stiffness: f64,
    /// Pennation angle α at optimal length (rad).
    pub pennation_angle: f64,
    /// Force–velocity relation.
    pub fv: ForceVelocityRelation,
    /// Force–length relation.
    pub fl: ForceLengthRelation,
    /// Current activation (0–1).
    pub activation: f64,
    /// Current normalised CE length (l_ce / l_opt).
    pub l_ce_norm: f64,
    /// Current CE velocity (normalised: v / v_max).
    pub v_ce_norm: f64,
}
impl HillMuscleModel {
    /// Create a new Hill muscle model with default fast-twitch parameters.
    pub fn new(f_max: f64, l_opt: f64, tendon_slack_length: f64, pennation_angle: f64) -> Self {
        Self {
            f_max,
            l_opt,
            tendon_slack_length,
            tendon_stiffness: f_max / (0.04 * tendon_slack_length),
            pennation_angle,
            fv: ForceVelocityRelation::fast_twitch(),
            fl: ForceLengthRelation::default_params(),
            activation: 0.0,
            l_ce_norm: 1.0,
            v_ce_norm: 0.0,
        }
    }
    /// Typical tibialis anterior.
    pub fn tibialis_anterior() -> Self {
        Self::new(600.0, 0.098, 0.223, 0.105)
    }
    /// Typical soleus.
    pub fn soleus() -> Self {
        Self::new(3550.0, 0.030, 0.268, 0.349)
    }
    /// Active CE force (N) at current state.
    pub fn ce_force_active(&self) -> f64 {
        let fl = self.fl.fl_active(self.l_ce_norm);
        let fv = self.fv.fv(self.v_ce_norm);
        self.f_max * self.activation * fl * fv
    }
    /// Passive (PEE) force (N) at current CE length.
    pub fn pee_force(&self) -> f64 {
        let fp = self.fl.fl_passive(self.l_ce_norm);
        self.f_max * fp
    }
    /// Total CE force (active + passive, N).
    pub fn ce_force_total(&self) -> f64 {
        self.ce_force_active() + self.pee_force()
    }
    /// SEE (tendon) force (N) given musculotendon length l_mt (m).
    pub fn see_force(&self, l_mt: f64) -> f64 {
        let l_ts = self.tendon_slack_length;
        let cos_alpha = self.pennation_angle.cos();
        let l_ce = self.l_ce_norm * self.l_opt;
        let l_tendon = l_mt - l_ce * cos_alpha;
        if l_tendon <= l_ts {
            0.0
        } else {
            self.tendon_stiffness * (l_tendon - l_ts)
        }
    }
    /// Equilibrium: given l_mt, find the normalised v_ce that balances CE and SEE forces.
    ///
    /// Returns normalised velocity (negative = shortening).
    pub fn compute_v_ce(&self, l_mt: f64) -> f64 {
        let f_see = self.see_force(l_mt);
        let f_ce_passive = self.pee_force();
        let fl = self.fl.fl_active(self.l_ce_norm);
        if fl < 1e-10 || self.activation < 1e-10 {
            return 0.0;
        }
        let f_active_target = (f_see - f_ce_passive).max(0.0);
        let f_iso = self.f_max * self.activation * fl;
        if f_iso < 1e-10 {
            return 0.0;
        }
        let fv_target = f_active_target / f_iso;
        let a = self.fv.a_over_f0;
        let fv_clamped = clamp(fv_target, 0.0, 1.0);
        let vn = (1.0 - fv_clamped) / (1.0 + fv_clamped / a);
        -vn
    }
    /// Integrate the muscle state forward by dt (s) given l_mt (m) and
    /// target activation `a_target`.
    pub fn step(&mut self, dt: f64, l_mt: f64, a_target: f64) {
        let tau = if a_target > self.activation {
            0.01
        } else {
            0.04
        };
        self.activation += (a_target - self.activation) * dt / tau;
        self.activation = clamp(self.activation, 0.0, 1.0);
        self.v_ce_norm = self.compute_v_ce(l_mt);
        let dl = self.v_ce_norm * self.fv.v_max * dt / self.l_opt;
        self.l_ce_norm += dl;
        self.l_ce_norm = clamp(self.l_ce_norm, 0.2, 1.8);
    }
}
/// A single muscle entry in the optimization problem.
#[derive(Debug, Clone)]
pub struct OptMuscle {
    /// Muscle name.
    pub name: String,
    /// Musculotendon unit.
    pub mtu: MusculotendonUnit,
    /// Moment arm (m) for this muscle at the reference joint.
    pub moment_arm: f64,
    /// Current activation (state variable).
    pub activation: f64,
}
impl OptMuscle {
    /// Create a new optimization muscle entry.
    pub fn new(name: &str, mtu: MusculotendonUnit, moment_arm: f64) -> Self {
        Self {
            name: name.to_string(),
            mtu,
            moment_arm: moment_arm.abs(),
            activation: 0.0,
        }
    }
    /// Maximum torque this muscle can produce: F_max * moment_arm (N·m).
    pub fn max_torque(&self) -> f64 {
        self.mtu.muscle.f_max * self.moment_arm
    }
}
/// State of a half-sarcomere including cross-bridge occupancy and length.
///
/// Implements a simplified Huxley (1957) two-state cross-bridge model:
/// actin-myosin cross-bridges transition between detached and attached states
/// with rate constants `f` (attachment) and `g` (detachment) that depend on
/// the relative sliding velocity.
#[derive(Debug, Clone)]
pub struct SarcomereModel {
    /// Fraction of cross-bridges in the attached (force-generating) state.
    pub attached_fraction: f64,
    /// Current sarcomere length (m).
    pub sarcomere_length: f64,
    /// Optimal sarcomere length for maximum force production (m).
    pub optimal_sarcomere_length: f64,
    /// Attachment rate constant f1 (1/s).
    pub attach_rate: f64,
    /// Detachment rate constant g1 during shortening (1/s).
    pub detach_rate_shortening: f64,
    /// Detachment rate constant g2 during lengthening (1/s).
    pub detach_rate_lengthening: f64,
    /// Maximum force per unit cross-sectional area at full attachment (N/m²).
    pub max_isometric_stress: f64,
    /// Half-width of the length-dependent active force plateau (normalised).
    pub length_plateau_width: f64,
}
impl SarcomereModel {
    /// Advance the cross-bridge kinetics by `dt` (s).
    ///
    /// `activation` ∈ \[0,1\] scales the attachment rate (neural drive).
    /// `shortening_velocity` is positive for shortening, negative for lengthening (normalised, m/s / l_opt·s⁻¹).
    pub fn step(&mut self, dt: f64, activation: f64, shortening_velocity: f64) {
        let f = self.attach_rate * activation.clamp(0.0, 1.0);
        let g = if shortening_velocity >= 0.0 {
            self.detach_rate_shortening * (1.0 + shortening_velocity.abs())
        } else {
            self.detach_rate_lengthening * (1.0 + shortening_velocity.abs() * 0.5)
        };
        let dn = f * (1.0 - self.attached_fraction) - g * self.attached_fraction;
        self.attached_fraction = clamp(self.attached_fraction + dn * dt, 0.0, 1.0);
    }
    /// Active force from cross-bridge attachment, scaled by length-dependence.
    ///
    /// Returns force in N per unit area (Pa) — multiply by PCSA for absolute force.
    pub fn active_force(&self) -> f64 {
        let lf = self.length_dependent_force(self.sarcomere_length);
        lf * self.attached_fraction * self.max_isometric_stress
    }
    /// Length-dependent scaling factor for active force (FL relation, \[0,1\]).
    ///
    /// Uses a trapezoidal plateau centred on `optimal_sarcomere_length`.
    pub fn length_dependent_force(&self, length: f64) -> f64 {
        let l_norm = length / self.optimal_sarcomere_length;
        let w = self.length_plateau_width;
        if (l_norm - 1.0).abs() <= w {
            1.0
        } else if l_norm < 1.0 - w {
            ((l_norm - 0.5) / (0.5 - w)).clamp(0.0, 1.0)
        } else {
            ((1.8 - l_norm) / (0.8 - w)).clamp(0.0, 1.0)
        }
    }
    /// Passive elastic force from titin/myosin filaments (normalised, >= 0).
    pub fn passive_force(&self) -> f64 {
        let l_norm = self.sarcomere_length / self.optimal_sarcomere_length;
        if l_norm > 1.0 {
            (l_norm - 1.0).powi(2) * 0.5 * self.max_isometric_stress
        } else {
            0.0
        }
    }
}
