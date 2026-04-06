//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;

/// Parameters for a single-DOF angular motor (e.g., on a revolute joint).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AngularMotorData {
    /// Motor control mode.
    pub mode: AngularMotorMode,
    /// Target angle \[rad\] (used in `Position` mode).
    pub target_angle: f64,
    /// Target angular velocity \[rad/s\] (used in `Velocity` mode).
    pub target_velocity: f64,
    /// Maximum torque the motor can apply \[N·m\].
    pub max_torque: f64,
    /// Proportional gain for position mode.
    pub kp: f64,
    /// Derivative gain for position mode (velocity damping).
    pub kd: f64,
    /// Current motor torque (output).
    pub current_torque: f64,
}
impl AngularMotorData {
    /// Create a new angular motor in `Off` mode.
    #[allow(dead_code)]
    pub fn new(max_torque: f64) -> Self {
        AngularMotorData {
            mode: AngularMotorMode::Off,
            target_angle: 0.0,
            target_velocity: 0.0,
            max_torque,
            kp: 100.0,
            kd: 10.0,
            current_torque: 0.0,
        }
    }
    /// Compute the motor torque for the current state.
    ///
    /// * `current_angle`    — current joint angle \[rad\].
    /// * `current_velocity` — current angular velocity \[rad/s\].
    /// * `dt`               — time step \[s\].
    #[allow(dead_code)]
    pub fn compute_torque(&mut self, current_angle: f64, current_velocity: f64, _dt: f64) -> f64 {
        let torque = match self.mode {
            AngularMotorMode::Off => 0.0,
            AngularMotorMode::Velocity => {
                let err = self.target_velocity - current_velocity;
                self.kd * err
            }
            AngularMotorMode::Position => {
                let pos_err = self.target_angle - current_angle;
                let vel_err = 0.0 - current_velocity;
                self.kp * pos_err + self.kd * vel_err
            }
        };
        let clamped = torque.clamp(-self.max_torque, self.max_torque);
        self.current_torque = clamped;
        clamped
    }
    /// Returns `true` if the joint is within `tolerance` \[rad\] of the target.
    #[allow(dead_code)]
    pub fn at_target(&self, current_angle: f64, tolerance: f64) -> bool {
        (current_angle - self.target_angle).abs() <= tolerance
    }
}
/// Supported norms for measuring constraint error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum ErrorNorm {
    /// L1 norm (sum of absolute values).
    L1,
    /// L2 norm (Euclidean / RMS).
    L2,
    /// L∞ norm (maximum absolute value).
    LInfinity,
    /// Weighted L2 norm (each component scaled by a weight).
    WeightedL2,
}
/// A minimal cache for the previous frame's constraint impulses, used for
/// warm-starting and stability analysis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ImpulseCache {
    /// Cached impulses keyed by constraint index.
    pub(super) entries: Vec<(usize, f64)>,
    /// Maximum number of entries to retain.
    pub(super) capacity: usize,
}
impl ImpulseCache {
    /// Create a new cache with the given capacity.
    #[allow(dead_code)]
    pub fn new(capacity: usize) -> Self {
        ImpulseCache {
            entries: Vec::with_capacity(capacity),
            capacity,
        }
    }
    /// Store an impulse for a constraint index, evicting the oldest entry if
    /// at capacity.
    #[allow(dead_code)]
    pub fn store(&mut self, index: usize, impulse: f64) {
        if let Some(pos) = self.entries.iter().position(|(i, _)| *i == index) {
            self.entries[pos].1 = impulse;
        } else {
            if self.entries.len() >= self.capacity {
                self.entries.remove(0);
            }
            self.entries.push((index, impulse));
        }
    }
    /// Retrieve the cached impulse for a constraint index, returning `None` if
    /// not found.
    #[allow(dead_code)]
    pub fn get(&self, index: usize) -> Option<f64> {
        self.entries
            .iter()
            .find(|(i, _)| *i == index)
            .map(|(_, v)| *v)
    }
    /// Clear all cached impulses.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.entries.clear();
    }
    /// Number of entries currently stored.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Returns `true` if the cache is empty.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}
/// XPBD (eXtended Position Based Dynamics) parameters for a soft constraint.
///
/// The XPBD update rule is:
///
/// ```text
/// Δλ = -(C + α̃·λ + β̃·Ċ) / (∇C M⁻¹ ∇Cᵀ + α̃)
/// α̃ = α / dt²
/// β̃ = β / dt
/// ```
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct XpbdParams {
    /// Compliance α \[m/N or rad/(N·m)\].  0 → rigid.
    pub compliance: f64,
    /// Damping β \[s/m or s·rad/(N·m)\].  0 → undamped.
    pub damping: f64,
}
impl XpbdParams {
    /// Create XPBD params with given compliance and damping.
    #[allow(dead_code)]
    pub fn new(compliance: f64, damping: f64) -> Self {
        XpbdParams {
            compliance,
            damping,
        }
    }
    /// Rigid (zero compliance, zero damping).
    #[allow(dead_code)]
    pub fn rigid() -> Self {
        XpbdParams {
            compliance: 0.0,
            damping: 0.0,
        }
    }
    /// Scaled compliance: α̃ = α / dt².
    #[allow(dead_code)]
    pub fn alpha_tilde(&self, dt: f64) -> f64 {
        if dt.abs() < 1e-15 {
            0.0
        } else {
            self.compliance / (dt * dt)
        }
    }
    /// Scaled damping: β̃ = β / dt.
    #[allow(dead_code)]
    pub fn beta_tilde(&self, dt: f64) -> f64 {
        if dt.abs() < 1e-15 {
            0.0
        } else {
            self.damping / dt
        }
    }
    /// Compute the full XPBD Δλ.
    ///
    /// * `c`          — positional constraint residual C(x).
    /// * `c_dot`      — velocity residual Ċ (used for damping term).
    /// * `jmj`        — J M⁻¹ Jᵀ (generalized inverse mass).
    /// * `lambda`     — current accumulated Lagrange multiplier.
    /// * `dt`         — sub-step time.
    #[allow(dead_code)]
    pub fn delta_lambda(&self, c: f64, c_dot: f64, jmj: f64, lambda: f64, dt: f64) -> f64 {
        let alpha_t = self.alpha_tilde(dt);
        let beta_t = self.beta_tilde(dt);
        let denom = jmj + alpha_t;
        if denom.abs() < 1e-15 {
            return 0.0;
        }
        -(c + alpha_t * lambda + beta_t * c_dot) / denom
    }
}
/// A compact, copy-friendly identifier for the type of a constraint.
///
/// Used for serialization, debug output, and solver statistics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[allow(dead_code)]
pub enum ConstraintKind {
    /// Fixed joint (zero relative DOF).
    Fixed,
    /// Ball-and-socket joint (3 translational DOF locked).
    Ball,
    /// Revolute (hinge) joint.
    Revolute,
    /// Prismatic (sliding) joint.
    Prismatic,
    /// Spring constraint.
    Spring,
    /// Contact / collision constraint.
    Contact,
    /// 6-DOF generic joint.
    SixDof,
    /// Gear ratio coupling.
    Gear,
    /// Pulley constraint.
    Pulley,
    /// Rack-and-pinion coupling.
    RackPinion,
    /// Motor-driven DOF.
    Motor,
    /// PBD distance / shape constraint.
    Pbd,
    /// User-defined extension.
    Custom(u32),
}
/// A named group of constraint indices that should be solved together.
///
/// Groups are used to enforce inter-constraint ordering (e.g. solve all
/// structural joints before contact constraints).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConstraintGroup {
    /// Human-readable group name (for debug output).
    pub name: String,
    /// Indices into the solver's constraint array.
    pub indices: Vec<usize>,
    /// Priority for the group as a whole.
    pub priority: ConstraintPriority,
    /// Whether all constraints in the group must converge for the group to pass.
    pub require_all: bool,
}
impl ConstraintGroup {
    /// Create a new empty group.
    #[allow(dead_code)]
    pub fn new(name: impl Into<String>, priority: ConstraintPriority) -> Self {
        ConstraintGroup {
            name: name.into(),
            indices: Vec::new(),
            priority,
            require_all: true,
        }
    }
    /// Add a constraint index to the group.
    #[allow(dead_code)]
    pub fn push(&mut self, index: usize) {
        self.indices.push(index);
    }
    /// Return the number of constraints in this group.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.indices.len()
    }
    /// Returns `true` if the group contains no constraints.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}
/// A validated compliance range `[min, max]` for use in soft constraints.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub struct ComplianceRange {
    /// Minimum compliance value (most rigid end).
    pub min: f64,
    /// Maximum compliance value (most compliant end).
    pub max: f64,
}
impl ComplianceRange {
    /// Create a new compliance range, ensuring `min ≤ max` and both > 0.
    ///
    /// Panics in debug mode if `min > max` or either value is non-positive.
    #[allow(dead_code)]
    pub fn new(min: f64, max: f64) -> Self {
        debug_assert!(min > 0.0, "compliance min must be positive");
        debug_assert!(max >= min, "compliance max must be >= min");
        ComplianceRange { min, max }
    }
    /// Clamp a compliance value to this range.
    #[allow(dead_code)]
    pub fn clamp(&self, value: f64) -> f64 {
        value.clamp(self.min, self.max)
    }
    /// Linearly interpolate within the range: `t=0` → `min`, `t=1` → `max`.
    #[allow(dead_code)]
    pub fn lerp(&self, t: f64) -> f64 {
        let t = t.clamp(0.0, 1.0);
        self.min + t * (self.max - self.min)
    }
    /// Returns the midpoint of the range.
    #[allow(dead_code)]
    pub fn midpoint(&self) -> f64 {
        (self.min + self.max) * 0.5
    }
}
/// A one-sided (or two-sided) distance constraint between two points.
///
/// The constraint is active when the distance between the anchors violates
/// the limit: `distance < lower_limit` or `distance > upper_limit`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct DistanceLimitData {
    /// Anchor A in body-A local space.
    pub anchor_a: [f64; 3],
    /// Anchor B in body-B local space.
    pub anchor_b: [f64; 3],
    /// Minimum allowed distance (0.0 → no lower limit).
    pub lower_limit: f64,
    /// Maximum allowed distance (f64::INFINITY → no upper limit).
    pub upper_limit: f64,
    /// Unilateral multiplier for the lower limit.
    pub lambda_lower: LagrangeMultiplier,
    /// Unilateral multiplier for the upper limit.
    pub lambda_upper: LagrangeMultiplier,
}
impl DistanceLimitData {
    /// Create a distance limit with given bounds.
    #[allow(dead_code)]
    pub fn new(anchor_a: [f64; 3], anchor_b: [f64; 3], lower: f64, upper: f64) -> Self {
        DistanceLimitData {
            anchor_a,
            anchor_b,
            lower_limit: lower,
            upper_limit: upper,
            lambda_lower: LagrangeMultiplier::unilateral(),
            lambda_upper: LagrangeMultiplier::with_bounds(f64::NEG_INFINITY, 0.0),
        }
    }
    /// Returns the violation at the current distance, or 0.0 if within limits.
    ///
    /// Positive = below lower limit, Negative = above upper limit, 0 = within.
    #[allow(dead_code)]
    pub fn violation(&self, distance: f64) -> f64 {
        if distance < self.lower_limit {
            self.lower_limit - distance
        } else if distance > self.upper_limit {
            self.upper_limit - distance
        } else {
            0.0
        }
    }
    /// Returns `true` if the lower limit is currently active.
    #[allow(dead_code)]
    pub fn lower_active(&self, distance: f64) -> bool {
        distance < self.lower_limit
    }
    /// Returns `true` if the upper limit is currently active.
    #[allow(dead_code)]
    pub fn upper_active(&self, distance: f64) -> bool {
        distance > self.upper_limit
    }
}
/// Tracks residuals for a multi-DOF constraint across iterations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConstraintResidualTracker {
    /// Per-DOF residuals from the last iteration.
    pub residuals: Vec<f64>,
    /// Norm type used for convergence checking.
    pub norm: ErrorNorm,
    /// Convergence threshold.
    pub tolerance: f64,
    /// Number of iterations performed since last reset.
    pub iteration_count: usize,
}
impl ConstraintResidualTracker {
    /// Create a new tracker for `dof_count` DOFs.
    #[allow(dead_code)]
    pub fn new(dof_count: usize, norm: ErrorNorm, tolerance: f64) -> Self {
        ConstraintResidualTracker {
            residuals: vec![0.0; dof_count],
            norm,
            tolerance,
            iteration_count: 0,
        }
    }
    /// Update residuals for the current iteration.
    #[allow(dead_code)]
    pub fn update(&mut self, new_residuals: &[f64]) {
        let len = self.residuals.len().min(new_residuals.len());
        self.residuals[..len].copy_from_slice(&new_residuals[..len]);
        self.iteration_count += 1;
    }
    /// Compute the current error norm.
    #[allow(dead_code)]
    pub fn error_norm(&self) -> f64 {
        compute_error_norm(&self.residuals, self.norm)
    }
    /// Returns `true` if the current error norm is below the tolerance.
    #[allow(dead_code)]
    pub fn has_converged(&self) -> bool {
        self.error_norm() <= self.tolerance
    }
    /// Reset the tracker.
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        for r in &mut self.residuals {
            *r = 0.0;
        }
        self.iteration_count = 0;
    }
}
/// Snapshot of per-constraint solver statistics for one physics step.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct ConstraintMetrics {
    /// Number of velocity-solve iterations performed.
    pub velocity_iterations: usize,
    /// Final velocity residual after convergence.
    pub velocity_residual: f64,
    /// Position correction magnitude applied.
    pub position_correction: f64,
    /// Accumulated impulse magnitude over the step.
    pub accumulated_impulse: f64,
    /// Whether the constraint converged within tolerance.
    pub converged: bool,
}
impl ConstraintMetrics {
    /// Create a new zeroed metrics snapshot.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Returns `true` if the residual is below the given tolerance.
    #[allow(dead_code)]
    pub fn is_within_tolerance(&self, tolerance: f64) -> bool {
        self.velocity_residual <= tolerance
    }
    /// Merge another metrics snapshot into this one (accumulate counts/maxima).
    #[allow(dead_code)]
    pub fn merge(&mut self, other: &ConstraintMetrics) {
        self.velocity_iterations += other.velocity_iterations;
        self.velocity_residual = self.velocity_residual.max(other.velocity_residual);
        self.position_correction = self.position_correction.max(other.position_correction);
        self.accumulated_impulse += other.accumulated_impulse;
        self.converged = self.converged && other.converged;
    }
}
/// Priority level used to order constraints within a solve island.
///
/// Higher-priority constraints are solved first, giving them precedence when
/// the solver cannot satisfy all constraints simultaneously.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default)]
#[allow(dead_code)]
pub enum ConstraintPriority {
    /// Lowest priority — may be violated if higher-priority constraints conflict.
    Low = 0,
    /// Normal priority (default).
    #[default]
    Normal = 1,
    /// High priority — enforced before normal/low constraints.
    High = 2,
    /// Critical — must never be violated (e.g. structural joint, collision).
    Critical = 3,
}
/// Lifecycle state of a constraint within the solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub enum ConstraintState {
    /// Not yet initialised; will be prepared on the next step.
    #[default]
    Uninitialised,
    /// Active and participating in the solve.
    Active,
    /// Temporarily disabled; will be skipped by the solver.
    Disabled,
    /// Permanently deactivated (e.g. broken, out of scope).
    Removed,
    /// Sleeping — bodies are at rest; no solve needed.
    Sleeping,
}
/// A spring constraint between two anchor points on two bodies.
///
/// Enforces `|r_a - r_b| ≈ rest_length` using a linear spring force
/// with optional XPBD compliance.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SpringConstraintData {
    /// Anchor point in body-A local space.
    pub anchor_a: [f64; 3],
    /// Anchor point in body-B local space.
    pub anchor_b: [f64; 3],
    /// Rest length of the spring \[m\].
    pub rest_length: f64,
    /// Spring stiffness \[N/m\].
    pub stiffness: f64,
    /// Damping coefficient \[N·s/m\].
    pub damping_coeff: f64,
    /// XPBD compliance (0 = rigid).
    pub xpbd: XpbdParams,
    /// Accumulated Lagrange multiplier.
    pub lambda: LagrangeMultiplier,
    /// Whether the spring is currently active.
    pub active: bool,
}
impl SpringConstraintData {
    /// Create a new spring constraint.
    #[allow(dead_code)]
    pub fn new(
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        rest_length: f64,
        stiffness: f64,
        damping: f64,
    ) -> Self {
        SpringConstraintData {
            anchor_a,
            anchor_b,
            rest_length,
            stiffness,
            damping_coeff: damping,
            xpbd: XpbdParams::new(1.0 / stiffness.max(1e-12), damping / stiffness.max(1e-12)),
            lambda: LagrangeMultiplier::bilateral(),
            active: true,
        }
    }
    /// Current spring force magnitude given a current separation `distance`.
    #[allow(dead_code)]
    pub fn spring_force(&self, distance: f64) -> f64 {
        self.stiffness * (distance - self.rest_length)
    }
    /// Compute XPBD Δλ for this spring given current extension and velocity.
    #[allow(dead_code)]
    pub fn xpbd_delta(&mut self, extension: f64, ext_rate: f64, jmj: f64, dt: f64) -> f64 {
        self.xpbd
            .delta_lambda(extension, ext_rate, jmj, self.lambda.lambda, dt)
    }
}
/// Events emitted by constraints during the solve.
///
/// Solvers can broadcast these to registered listeners for game logic, sound,
/// visual effects, etc.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub enum ConstraintEvent {
    /// The constraint was broken (force exceeded threshold).
    Broken {
        /// Index of the constraint that broke.
        constraint_index: usize,
        /// Force at the time of breaking \[N\].
        force: f64,
    },
    /// A contact constraint began (bodies started touching).
    ContactBegin {
        /// Index of the contact constraint.
        constraint_index: usize,
    },
    /// A contact constraint ended (bodies separated).
    ContactEnd {
        /// Index of the contact constraint.
        constraint_index: usize,
    },
    /// A joint limit was hit.
    LimitHit {
        /// Index of the limit constraint.
        constraint_index: usize,
        /// Whether the lower (`true`) or upper (`false`) limit was hit.
        is_lower: bool,
    },
    /// The solver failed to converge for this constraint.
    ConvergenceWarning {
        /// Index of the failing constraint.
        constraint_index: usize,
        /// Residual at termination.
        residual: f64,
    },
}
/// Hints that a constraint can supply to the solver for optimisation.
///
/// These are advisory; the solver is free to ignore any hint.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct SolverHints {
    /// Suggested maximum number of velocity-solve iterations for this constraint.
    ///
    /// `None` means "use the solver default".
    pub max_velocity_iterations: Option<usize>,
    /// Whether this constraint benefits from position-level stabilization.
    pub needs_position_solve: bool,
    /// Whether the constraint is one-sided (inequality / contact).
    pub is_unilateral: bool,
    /// Priority override (overrides [`PrioritizedConstraint::priority`] if `Some`).
    pub priority_override: Option<ConstraintPriority>,
}
/// Control mode for an angular motor.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub enum AngularMotorMode {
    /// Drive toward a target angular velocity.
    Velocity,
    /// Drive toward a target angle.
    Position,
    /// No motor (passive joint).
    #[default]
    Off,
}
/// Accumulated Lagrange multiplier state for one constraint DOF.
///
/// Stores the current λ and applies clamping to enforce inequality constraints.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct LagrangeMultiplier {
    /// Current accumulated impulse (λ · Δt).
    pub lambda: f64,
    /// Lower bound on λ (use `f64::NEG_INFINITY` for equality constraints).
    pub lower_bound: f64,
    /// Upper bound on λ (use `f64::INFINITY` for equality constraints).
    pub upper_bound: f64,
}
impl LagrangeMultiplier {
    /// Create a new bilateral (equality) Lagrange multiplier.
    #[allow(dead_code)]
    pub fn bilateral() -> Self {
        LagrangeMultiplier {
            lambda: 0.0,
            lower_bound: f64::NEG_INFINITY,
            upper_bound: f64::INFINITY,
        }
    }
    /// Create a new unilateral (contact / one-sided) Lagrange multiplier.
    #[allow(dead_code)]
    pub fn unilateral() -> Self {
        LagrangeMultiplier {
            lambda: 0.0,
            lower_bound: 0.0,
            upper_bound: f64::INFINITY,
        }
    }
    /// Create with explicit bounds.
    #[allow(dead_code)]
    pub fn with_bounds(lower: f64, upper: f64) -> Self {
        LagrangeMultiplier {
            lambda: 0.0,
            lower_bound: lower,
            upper_bound: upper,
        }
    }
    /// Apply a delta-lambda update and return the actual change (after clamping).
    ///
    /// The accumulated `lambda` is updated in-place.
    #[allow(dead_code)]
    pub fn apply_delta(&mut self, delta: f64) -> f64 {
        let new_lambda = (self.lambda + delta).clamp(self.lower_bound, self.upper_bound);
        let actual_delta = new_lambda - self.lambda;
        self.lambda = new_lambda;
        actual_delta
    }
    /// Reset the multiplier to zero (call at the start of each frame).
    #[allow(dead_code)]
    pub fn reset(&mut self) {
        self.lambda = 0.0;
    }
    /// Apply warm-start: seed λ with a fraction of the previous value.
    #[allow(dead_code)]
    pub fn warm_start(&mut self, prev_lambda: f64, factor: f64) {
        self.lambda = (prev_lambda * factor).clamp(self.lower_bound, self.upper_bound);
    }
}
/// A cone limit restricting the angle between a body axis and a reference axis.
///
/// Commonly used for ball-and-socket joints to prevent over-extension.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ConeLimitData {
    /// Reference axis in world space (unit vector).
    pub reference_axis: [f64; 3],
    /// Half-angle of the cone \[rad\].
    pub half_angle: f64,
    /// Lagrange multiplier for the cone limit.
    pub lambda: LagrangeMultiplier,
    /// Whether the limit is currently active.
    pub active: bool,
}
impl ConeLimitData {
    /// Create a new cone limit.
    #[allow(dead_code)]
    pub fn new(reference_axis: [f64; 3], half_angle: f64) -> Self {
        ConeLimitData {
            reference_axis,
            half_angle,
            lambda: LagrangeMultiplier::unilateral(),
            active: false,
        }
    }
    /// Compute the angle between `body_axis` and the reference axis.
    #[allow(dead_code)]
    pub fn current_angle(&self, body_axis: [f64; 3]) -> f64 {
        let ref_len = (self.reference_axis[0] * self.reference_axis[0]
            + self.reference_axis[1] * self.reference_axis[1]
            + self.reference_axis[2] * self.reference_axis[2])
            .sqrt();
        let body_len = (body_axis[0] * body_axis[0]
            + body_axis[1] * body_axis[1]
            + body_axis[2] * body_axis[2])
            .sqrt();
        if ref_len < 1e-12 || body_len < 1e-12 {
            return 0.0;
        }
        let r = [
            self.reference_axis[0] / ref_len,
            self.reference_axis[1] / ref_len,
            self.reference_axis[2] / ref_len,
        ];
        let b = [
            body_axis[0] / body_len,
            body_axis[1] / body_len,
            body_axis[2] / body_len,
        ];
        let cos_angle = (r[0] * b[0] + r[1] * b[1] + r[2] * b[2]).clamp(-1.0, 1.0);
        cos_angle.acos()
    }
    /// Returns `true` if the cone limit is currently violated.
    #[allow(dead_code)]
    pub fn is_violated(&self, body_axis: [f64; 3]) -> bool {
        self.current_angle(body_axis) > self.half_angle
    }
    /// Cone constraint residual: angle - half_angle (positive when violated).
    #[allow(dead_code)]
    pub fn residual(&self, body_axis: [f64; 3]) -> f64 {
        (self.current_angle(body_axis) - self.half_angle).max(0.0)
    }
    /// Update the active flag based on current axis direction.
    #[allow(dead_code)]
    pub fn update_active(&mut self, body_axis: [f64; 3]) {
        self.active = self.is_violated(body_axis);
    }
}
