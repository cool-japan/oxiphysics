//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::pgs_solver::SolverStats;
use crate::traits::Constraint;
use oxiphysics_core::math::Vec3;
use oxiphysics_rigid::RigidBodySet;

/// TGS solver with adaptive substep selection.
///
/// Chooses the number of substeps based on the maximum velocity magnitude
/// and a target CFL-like condition: `max_velocity * sub_dt < max_displacement`.
#[derive(Debug, Clone)]
pub struct AdaptiveTgsSolver {
    /// Base TGS solver parameters.
    pub base: TgsSolver,
    /// Maximum allowed displacement per substep (metres).
    pub max_displacement_per_substep: f64,
    /// Minimum number of substeps.
    pub min_substeps: usize,
    /// Maximum number of substeps.
    pub max_substeps: usize,
}
impl AdaptiveTgsSolver {
    /// Create a new adaptive TGS solver.
    pub fn new(velocity_iterations: usize, gravity: Vec3, max_displacement: f64) -> Self {
        Self {
            base: TgsSolver::new(1, velocity_iterations, gravity),
            max_displacement_per_substep: max_displacement,
            min_substeps: MIN_SUBSTEPS,
            max_substeps: MAX_SUBSTEPS,
        }
    }
    /// Determine optimal substep count given maximum velocity and dt.
    pub fn compute_substeps(&self, max_velocity: f64, dt: f64) -> usize {
        if max_velocity < 1e-12 || dt < 1e-12 {
            return self.min_substeps;
        }
        let needed = (max_velocity * dt / self.max_displacement_per_substep).ceil() as usize;
        needed.clamp(self.min_substeps, self.max_substeps)
    }
    /// Solve with adaptive substepping.
    pub fn step(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        dt: f64,
    ) -> SolverStats {
        let max_vel = compute_max_velocity(bodies);
        let substeps = self.compute_substeps(max_vel, dt);
        let mut solver = self.base.clone();
        solver.substeps = substeps;
        solver.step(constraints, bodies, dt)
    }
}
/// Block TGS constraint: one row of Jacobian coupling two 6-DOF bodies.
///
/// Stores the constraint bias (Baumgarte + restitution + user bias),
/// the accumulated impulse (warm-start), and the effective mass.
#[derive(Debug, Clone)]
pub struct BlockConstraint6 {
    /// Jacobian row.
    pub jacobian: JacobianRow6,
    /// Mass matrix for body A.
    pub mass_a: MassMatrix6,
    /// Mass matrix for body B.
    pub mass_b: MassMatrix6,
    /// Constraint bias (target velocity).
    pub bias: f64,
    /// Lower impulse bound.
    pub lambda_min: f64,
    /// Upper impulse bound.
    pub lambda_max: f64,
    /// Accumulated impulse (warm-start).
    pub lambda: f64,
    /// Precomputed effective mass.
    pub(super) eff_mass: f64,
}
impl BlockConstraint6 {
    /// Create a new block constraint.
    pub fn new(
        jacobian: JacobianRow6,
        mass_a: MassMatrix6,
        mass_b: MassMatrix6,
        bias: f64,
        lambda_min: f64,
        lambda_max: f64,
    ) -> Self {
        let eff_a = mass_a.effective_mass(&jacobian.j_a);
        let eff_b = mass_b.effective_mass(&jacobian.j_b);
        let eff_mass = eff_a + eff_b;
        Self {
            jacobian,
            mass_a,
            mass_b,
            bias,
            lambda_min,
            lambda_max,
            lambda: 0.0,
            eff_mass,
        }
    }
    /// Solve one iteration of the block constraint.
    ///
    /// Returns `(delta_lambda, impulse_a, impulse_b)`.
    pub fn solve_iteration(&mut self, va: &[f64; 6], vb: &[f64; 6]) -> (f64, [f64; 6], [f64; 6]) {
        if self.eff_mass < 1e-30 {
            return (0.0, [0.0; 6], [0.0; 6]);
        }
        let jv = self.jacobian.constraint_velocity(va, vb);
        let delta = -(jv + self.bias) / self.eff_mass;
        let lambda_prev = self.lambda;
        self.lambda = (self.lambda + delta).clamp(self.lambda_min, self.lambda_max);
        let d_lambda = self.lambda - lambda_prev;
        let imp_a = {
            let inv_jt = self.mass_a.apply_inv_j(&self.jacobian.j_a);
            let mut arr = [0.0f64; 6];
            for i in 0..6 {
                arr[i] = inv_jt[i] * d_lambda;
            }
            arr
        };
        let imp_b = {
            let inv_jt = self.mass_b.apply_inv_j(&self.jacobian.j_b);
            let mut arr = [0.0f64; 6];
            for i in 0..6 {
                arr[i] = inv_jt[i] * d_lambda;
            }
            arr
        };
        (d_lambda, imp_a, imp_b)
    }
    /// Reset accumulated impulse.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
}
/// How position errors are corrected at the velocity level.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PositionCorrectionMode {
    /// Classic Baumgarte: bias velocity added to constraint.
    Baumgarte,
    /// Split impulse: separate pseudo-velocity for position correction.
    SplitImpulse,
    /// NGS (Non-linear Gauss-Seidel): direct position projection.
    NonlinearGaussSeidel,
}
/// A 6×6 block mass matrix (diagonal for simplicity: inv_mass, inv_inertia).
#[derive(Debug, Clone, Copy)]
pub struct MassMatrix6 {
    /// Inverse mass (scalar, applied to linear DOFs).
    pub inv_mass: f64,
    /// Inverse inertia (diagonal 3×3 stored as array \[Ixx, Iyy, Izz\]).
    pub inv_inertia: [f64; 3],
}
impl MassMatrix6 {
    /// Create a new 6×6 diagonal mass matrix.
    pub fn new(inv_mass: f64, inv_inertia: [f64; 3]) -> Self {
        Self {
            inv_mass,
            inv_inertia,
        }
    }
    /// Static/kinematic body: all zeros.
    pub fn static_body() -> Self {
        Self {
            inv_mass: 0.0,
            inv_inertia: [0.0; 3],
        }
    }
    /// Compute M⁻¹ · J^T for a single Jacobian row.
    ///
    /// Returns the 6-vector `M⁻¹ · J^T`.
    pub fn apply_inv_j(&self, j: &[f64; 6]) -> [f64; 6] {
        [
            self.inv_mass * j[0],
            self.inv_mass * j[1],
            self.inv_mass * j[2],
            self.inv_inertia[0] * j[3],
            self.inv_inertia[1] * j[4],
            self.inv_inertia[2] * j[5],
        ]
    }
    /// Compute the effective mass (scalar): J · M⁻¹ · J^T.
    pub fn effective_mass(&self, j: &[f64; 6]) -> f64 {
        let inv_jt = self.apply_inv_j(j);
        j.iter().zip(inv_jt.iter()).map(|(a, b)| a * b).sum()
    }
}
/// Strategy for constraint correction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CorrectionStrategy {
    /// Velocity-level only (Baumgarte bias injected into velocity constraint).
    VelocityLevel,
    /// Position-level only (direct position correction, no velocity bias).
    PositionLevel,
    /// Hybrid: velocity correction for small errors, position for large.
    Hybrid {
        /// Penetration threshold below which velocity correction is used.
        threshold: f64,
    },
}
impl CorrectionStrategy {
    /// Determine whether to apply velocity correction for a given penetration.
    pub fn use_velocity_correction(&self, penetration: f64) -> bool {
        match *self {
            CorrectionStrategy::VelocityLevel => true,
            CorrectionStrategy::PositionLevel => false,
            CorrectionStrategy::Hybrid { threshold } => penetration <= threshold,
        }
    }
    /// Determine whether to apply position correction for a given penetration.
    pub fn use_position_correction(&self, penetration: f64) -> bool {
        match *self {
            CorrectionStrategy::VelocityLevel => false,
            CorrectionStrategy::PositionLevel => true,
            CorrectionStrategy::Hybrid { threshold } => penetration > threshold,
        }
    }
}
/// Descriptor for one simulation island (group of mutually interacting bodies).
#[derive(Debug, Clone)]
pub struct IslandDescriptor {
    /// Handle indices of bodies in this island.
    pub body_indices: Vec<usize>,
    /// Constraint indices belonging to this island.
    pub constraint_indices: Vec<usize>,
    /// Estimated FLOP count (used for load balancing).
    pub estimated_cost: usize,
}
impl IslandDescriptor {
    /// Create a new island descriptor.
    pub fn new(body_indices: Vec<usize>, constraint_indices: Vec<usize>) -> Self {
        let cost = body_indices.len() * 6 + constraint_indices.len() * 12;
        Self {
            body_indices,
            constraint_indices,
            estimated_cost: cost,
        }
    }
    /// Number of bodies in the island.
    pub fn body_count(&self) -> usize {
        self.body_indices.len()
    }
    /// Number of constraints in the island.
    pub fn constraint_count(&self) -> usize {
        self.constraint_indices.len()
    }
}
/// A 6-DOF block row of a constraint Jacobian (3 linear + 3 angular).
///
/// Stores `J_a` and `J_b` as 6-element arrays: `[Jlin_x, Jlin_y, Jlin_z,
/// Jang_x, Jang_y, Jang_z]` for body A and body B respectively.
#[derive(Debug, Clone, Copy)]
pub struct JacobianRow6 {
    /// Jacobian for body A (linear + angular, 6 elements).
    pub j_a: [f64; 6],
    /// Jacobian for body B (linear + angular, 6 elements).
    pub j_b: [f64; 6],
}
impl JacobianRow6 {
    /// Create a zero Jacobian row.
    pub fn zero() -> Self {
        Self {
            j_a: [0.0; 6],
            j_b: [0.0; 6],
        }
    }
    /// Create from separate linear and angular parts.
    pub fn from_parts(lin_a: [f64; 3], ang_a: [f64; 3], lin_b: [f64; 3], ang_b: [f64; 3]) -> Self {
        Self {
            j_a: [lin_a[0], lin_a[1], lin_a[2], ang_a[0], ang_a[1], ang_a[2]],
            j_b: [lin_b[0], lin_b[1], lin_b[2], ang_b[0], ang_b[1], ang_b[2]],
        }
    }
    /// Compute J · v for a 6-vector `v = [lin_x, lin_y, lin_z, ang_x, ang_y, ang_z]`.
    pub fn dot_a(&self, v: &[f64; 6]) -> f64 {
        self.j_a.iter().zip(v.iter()).map(|(j, x)| j * x).sum()
    }
    /// Compute J · v for body B's velocity 6-vector.
    pub fn dot_b(&self, v: &[f64; 6]) -> f64 {
        self.j_b.iter().zip(v.iter()).map(|(j, x)| j * x).sum()
    }
    /// Constraint velocity (relative): J_a · v_a + J_b · v_b.
    pub fn constraint_velocity(&self, va: &[f64; 6], vb: &[f64; 6]) -> f64 {
        self.dot_a(va) + self.dot_b(vb)
    }
}
/// Parallel island solver planner: decomposes a set of islands into
/// independent work units that can be executed in parallel.
#[derive(Debug, Clone, Default)]
pub struct ParallelIslandPlanner {
    /// List of islands to solve.
    pub islands: Vec<IslandDescriptor>,
}
impl ParallelIslandPlanner {
    /// Create an empty planner.
    pub fn new() -> Self {
        Self {
            islands: Vec::new(),
        }
    }
    /// Add an island to the planner.
    pub fn add_island(&mut self, island: IslandDescriptor) {
        self.islands.push(island);
    }
    /// Sort islands by estimated cost (descending) for greedy load balancing.
    pub fn sort_by_cost(&mut self) {
        self.islands
            .sort_by_key(|b| std::cmp::Reverse(b.estimated_cost));
    }
    /// Assign islands to `n_workers` threads using a greedy approach.
    ///
    /// Returns a `Vec` of length `n_workers`, where each element is a list
    /// of island indices assigned to that worker.
    pub fn assign_workers(&self, n_workers: usize) -> Vec<Vec<usize>> {
        if n_workers == 0 {
            return vec![];
        }
        let mut assignments: Vec<Vec<usize>> = vec![Vec::new(); n_workers];
        let mut costs = vec![0usize; n_workers];
        for (i, island) in self.islands.iter().enumerate() {
            let worker = costs
                .iter()
                .enumerate()
                .min_by_key(|&(_, &c)| c)
                .map(|(w, _)| w)
                .unwrap_or(0);
            assignments[worker].push(i);
            costs[worker] += island.estimated_cost;
        }
        assignments
    }
    /// Total estimated work across all islands.
    pub fn total_cost(&self) -> usize {
        self.islands.iter().map(|i| i.estimated_cost).sum()
    }
}
/// Result of a drift correction pass.
#[derive(Debug, Clone, Copy)]
pub struct DriftCorrectionResult {
    /// Maximum position error before correction.
    pub max_error_before: f64,
    /// Maximum position error after correction.
    pub max_error_after: f64,
    /// Number of iterations performed.
    pub iterations: usize,
}
impl DriftCorrectionResult {
    /// Error reduction ratio (before/after). Returns 1.0 if before ≈ 0.
    pub fn error_reduction_ratio(&self) -> f64 {
        if self.max_error_before < 1e-20 {
            return 1.0;
        }
        self.max_error_after / self.max_error_before
    }
    /// True if error reduced by at least the given fraction.
    pub fn reduced_by(&self, fraction: f64) -> bool {
        self.error_reduction_ratio() <= 1.0 - fraction
    }
}
/// Diagnostics collected per substep for analysis.
#[derive(Debug, Clone, Copy)]
pub struct SubstepDiagnostics {
    /// Kinetic energy at start of substep.
    pub ke_start: f64,
    /// Kinetic energy at end of substep.
    pub ke_end: f64,
    /// Maximum velocity magnitude among all bodies.
    pub max_velocity: f64,
    /// Maximum penetration among constraints (estimated).
    pub max_penetration: f64,
}
impl SubstepDiagnostics {
    /// Create zeroed diagnostics.
    pub fn zero() -> Self {
        Self {
            ke_start: 0.0,
            ke_end: 0.0,
            max_velocity: 0.0,
            max_penetration: 0.0,
        }
    }
    /// Energy change ratio (ke_end / ke_start).
    pub fn energy_ratio(&self) -> f64 {
        if self.ke_start < 1e-20 {
            return 1.0;
        }
        self.ke_end / self.ke_start
    }
    /// Whether energy increased beyond the threshold.
    pub fn energy_gained(&self) -> bool {
        self.energy_ratio() > ENERGY_GAIN_THRESHOLD
    }
}
/// Baumgarte stabilization variant selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BaumgarteVariant {
    /// Classic Baumgarte: `bias = β * C / dt`.
    Classic,
    /// Scaled by ERP: `bias = erp * C / dt`.
    Erp,
    /// Two-parameter (overdamped): `bias = 2 * ζ * ω_n * C`.
    TwoParam {
        /// Damping ratio ζ.
        zeta: f64,
        /// Natural frequency ω_n (rad/s).
        omega_n: f64,
    },
}
impl BaumgarteVariant {
    /// Compute the bias velocity for a given penetration and time step.
    pub fn compute_bias(&self, penetration: f64, slop: f64, dt: f64) -> f64 {
        if dt < 1e-12 {
            return 0.0;
        }
        let excess = (penetration - slop).max(0.0);
        match *self {
            BaumgarteVariant::Classic => BAUMGARTE_BETA * excess / dt,
            BaumgarteVariant::Erp => {
                let erp = compute_erp(BAUMGARTE_BETA, dt);
                erp * excess / dt
            }
            BaumgarteVariant::TwoParam { zeta, omega_n } => 2.0 * zeta * omega_n * excess,
        }
    }
}
/// Warm-start state for the TGS solver.
///
/// Stores per-constraint accumulated impulses from the previous frame to
/// seed the next frame's solver, improving convergence for persistent contacts.
#[derive(Debug, Clone)]
pub struct TgsWarmStartState {
    /// Accumulated normal impulses per constraint.
    pub normal_impulses: Vec<f64>,
    /// Accumulated friction impulses per constraint (2 tangent directions).
    pub friction_impulses: Vec<[f64; 2]>,
    /// Previous time step (for dt-ratio scaling).
    pub prev_dt: f64,
    /// Number of frames this state has been carried forward.
    pub age: u32,
}
impl TgsWarmStartState {
    /// Create new warm-start state for `n` constraints.
    pub fn new(n: usize) -> Self {
        Self {
            normal_impulses: vec![0.0; n],
            friction_impulses: vec![[0.0; 2]; n],
            prev_dt: 1.0 / 60.0,
            age: 0,
        }
    }
    /// Reset all cached impulses.
    pub fn reset(&mut self) {
        for l in &mut self.normal_impulses {
            *l = 0.0;
        }
        for f in &mut self.friction_impulses {
            *f = [0.0; 2];
        }
        self.age = 0;
    }
    /// Scale all impulses by a factor (e.g., dt_ratio).
    pub fn scale(&mut self, factor: f64) {
        for l in &mut self.normal_impulses {
            *l *= factor;
        }
        for f in &mut self.friction_impulses {
            f[0] *= factor;
            f[1] *= factor;
        }
    }
    /// Resize to accommodate `n` constraints.
    pub fn resize(&mut self, n: usize) {
        self.normal_impulses.resize(n, 0.0);
        self.friction_impulses.resize(n, [0.0; 2]);
    }
    /// Total accumulated normal impulse magnitude.
    pub fn total_normal_impulse(&self) -> f64 {
        self.normal_impulses.iter().map(|l| l.abs()).sum()
    }
    /// Return the number of active constraints (non-zero normal impulse).
    pub fn active_count(&self) -> usize {
        self.normal_impulses
            .iter()
            .filter(|&&l| l.abs() > 1e-12)
            .count()
    }
    /// Age the state by one frame.
    pub fn advance_age(&mut self) {
        self.age += 1;
    }
    /// Check if the warm-start data is stale (older than `max_age` frames).
    pub fn is_stale(&self, max_age: u32) -> bool {
        self.age > max_age
    }
}
/// Temporal Gauss-Seidel solver.
///
/// Subdivides the time step into smaller substeps and solves constraints
/// at each substep.  This provides better stability for stiff constraints
/// and stacking scenarios compared to a single-step PGS solver.
#[derive(Debug, Clone)]
pub struct TgsSolver {
    /// Number of substeps per time step.
    pub substeps: usize,
    /// Number of velocity iterations per substep.
    pub velocity_iterations: usize,
    /// Number of position iterations per substep.
    pub position_iterations: usize,
    /// Gravity vector applied during integration.
    pub gravity: Vec3,
    /// Baumgarte stabilization coefficient (β).
    pub baumgarte_beta: f64,
    /// Enable energy-gain monitoring and damping.
    pub energy_monitoring: bool,
}
impl TgsSolver {
    /// Create a new TGS solver with explicit substep/iteration counts.
    pub fn new(substeps: usize, velocity_iterations: usize, gravity: Vec3) -> Self {
        Self {
            substeps,
            velocity_iterations,
            position_iterations: 1,
            gravity,
            baumgarte_beta: BAUMGARTE_BETA,
            energy_monitoring: false,
        }
    }
    /// Create a solver with full control over all parameters.
    pub fn with_options(
        substeps: usize,
        velocity_iterations: usize,
        position_iterations: usize,
        gravity: Vec3,
        baumgarte_beta: f64,
        energy_monitoring: bool,
    ) -> Self {
        Self {
            substeps,
            velocity_iterations,
            position_iterations,
            gravity,
            baumgarte_beta,
            energy_monitoring,
        }
    }
    /// Solve constraints over one full time step using substepping.
    ///
    /// Each substep:
    /// 1. Integrate forces for `sub_dt`
    /// 2. Prepare + solve velocity constraints
    /// 3. Solve position constraints (Baumgarte / projection)
    /// 4. Optionally monitor energy and damp if gain detected
    /// 5. Integrate positions for `sub_dt`
    ///
    /// Returns [`SolverStats`] aggregated over all substeps.
    pub fn step(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        dt: f64,
    ) -> SolverStats {
        let sub_dt = dt / self.substeps as f64;
        let mut total_iters = 0usize;
        for _ in 0..self.substeps {
            for (_handle, body) in bodies.iter_mut() {
                body.integrate_forces(sub_dt, &self.gravity);
            }
            for c in constraints.iter_mut() {
                c.prepare(bodies, sub_dt);
            }
            self.solve_velocity_constraints(constraints, bodies, sub_dt);
            total_iters += self.velocity_iterations;
            self.solve_position_constraints(constraints, bodies);
            if self.energy_monitoring {
                self.monitor_energy(bodies);
            }
            for (_handle, body) in bodies.iter_mut() {
                body.integrate_velocity(sub_dt);
            }
        }
        SolverStats {
            iterations_used: total_iters,
            residual: 0.0,
            converged: true,
        }
    }
    /// Solve with warm-start from a previous-frame [`TgsWarmStartState`].
    ///
    /// Scales cached impulses by `dt / prev_dt` before the first substep.
    pub fn step_warm(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        dt: f64,
        state: &mut TgsWarmStartState,
    ) -> SolverStats {
        if state.prev_dt > 1e-12 {
            let dt_ratio = dt / state.prev_dt;
            state.scale(dt_ratio);
        }
        state.advance_age();
        let stats = self.step(constraints, bodies, dt);
        state.prev_dt = dt;
        stats
    }
    /// Step with diagnostics collection.
    ///
    /// Returns both solver stats and a vector of per-substep diagnostics.
    pub fn step_with_diagnostics(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        dt: f64,
    ) -> (SolverStats, Vec<SubstepDiagnostics>) {
        let sub_dt = dt / self.substeps as f64;
        let mut total_iters = 0usize;
        let mut diagnostics = Vec::with_capacity(self.substeps);
        for _ in 0..self.substeps {
            let ke_start = compute_kinetic_energy(bodies);
            for (_handle, body) in bodies.iter_mut() {
                body.integrate_forces(sub_dt, &self.gravity);
            }
            for c in constraints.iter_mut() {
                c.prepare(bodies, sub_dt);
            }
            self.solve_velocity_constraints(constraints, bodies, sub_dt);
            total_iters += self.velocity_iterations;
            self.solve_position_constraints(constraints, bodies);
            if self.energy_monitoring {
                self.monitor_energy(bodies);
            }
            for (_handle, body) in bodies.iter_mut() {
                body.integrate_velocity(sub_dt);
            }
            let ke_end = compute_kinetic_energy(bodies);
            let max_vel = compute_max_velocity(bodies);
            diagnostics.push(SubstepDiagnostics {
                ke_start,
                ke_end,
                max_velocity: max_vel,
                max_penetration: 0.0,
            });
        }
        let stats = SolverStats {
            iterations_used: total_iters,
            residual: 0.0,
            converged: true,
        };
        (stats, diagnostics)
    }
    /// Run velocity-constraint resolution for a single substep.
    ///
    /// Iterates `velocity_iterations` times over all constraints, applying
    /// sequential impulses.
    pub fn solve_velocity_constraints(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        sub_dt: f64,
    ) {
        for _ in 0..self.velocity_iterations {
            for c in constraints.iter_mut() {
                c.solve_velocity(bodies, sub_dt);
            }
        }
    }
    /// Run position-constraint resolution for a single substep.
    ///
    /// Uses the Baumgarte stabilization term `β * penetration / dt` that is
    /// baked into each constraint's `solve_position`.
    pub fn solve_position_constraints(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
    ) {
        for _ in 0..self.position_iterations {
            for c in constraints.iter_mut() {
                c.solve_position(bodies, 0.0);
            }
        }
    }
    /// Solve a single island (a group of mutually interacting bodies).
    ///
    /// Identical to `step` but operates only on the provided body handles
    /// and constraints, enabling parallel island solving in the future.
    pub fn solve_island(
        &self,
        island_bodies: &[oxiphysics_core::BodyHandle],
        island_constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        dt: f64,
    ) -> SolverStats {
        let sub_dt = dt / self.substeps as f64;
        let mut total_iters = 0usize;
        for _ in 0..self.substeps {
            for &handle in island_bodies {
                if let Some(body) = bodies.get_mut(handle) {
                    body.integrate_forces(sub_dt, &self.gravity);
                }
            }
            for c in island_constraints.iter_mut() {
                c.prepare(bodies, sub_dt);
            }
            self.solve_velocity_constraints(island_constraints, bodies, sub_dt);
            total_iters += self.velocity_iterations;
            self.solve_position_constraints(island_constraints, bodies);
            if self.energy_monitoring {
                self.monitor_energy_island(island_bodies, bodies);
            }
            for &handle in island_bodies {
                if let Some(body) = bodies.get_mut(handle) {
                    body.integrate_velocity(sub_dt);
                }
            }
        }
        SolverStats {
            iterations_used: total_iters,
            residual: 0.0,
            converged: true,
        }
    }
    /// Compute the Baumgarte bias term for a penetration depth.
    ///
    /// `bias = β * max(0, penetration - slop) / dt`
    pub fn baumgarte_bias(&self, penetration: f64, dt: f64) -> f64 {
        if dt < 1e-12 {
            return 0.0;
        }
        let excess = (penetration - PENETRATION_SLOP).max(0.0);
        self.baumgarte_beta * excess / dt
    }
    /// Compute total kinetic energy for all bodies in the set.
    pub fn total_kinetic_energy(&self, bodies: &RigidBodySet) -> f64 {
        compute_kinetic_energy(bodies)
    }
    /// Compute the maximum velocity magnitude among all dynamic bodies.
    pub fn max_velocity(&self, bodies: &RigidBodySet) -> f64 {
        compute_max_velocity(bodies)
    }
    /// Detect and remove energy gain for all dynamic bodies.
    fn monitor_energy(&self, bodies: &mut RigidBodySet) {
        let handles: Vec<_> = bodies.iter().map(|(h, _)| h).collect();
        for handle in handles {
            clamp_body_energy_gain(bodies, handle, ENERGY_GAIN_THRESHOLD);
        }
    }
    /// Detect and remove energy gain for island bodies only.
    fn monitor_energy_island(
        &self,
        island_bodies: &[oxiphysics_core::BodyHandle],
        bodies: &mut RigidBodySet,
    ) {
        for &handle in island_bodies {
            clamp_body_energy_gain(bodies, handle, ENERGY_GAIN_THRESHOLD);
        }
    }
    /// Alias for `step` for compatibility with call-sites that use `solve`.
    pub fn solve(
        &self,
        constraints: &mut [Box<dyn Constraint>],
        bodies: &mut RigidBodySet,
        dt: f64,
    ) -> SolverStats {
        self.step(constraints, bodies, dt)
    }
}
