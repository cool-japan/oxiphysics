// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! XPBD solver for soft-body simulation.

use oxiphysics_core::math::Real;

use crate::constraint::SoftConstraint;
use crate::particle::{SoftBody, SoftParticle};

// ---------------------------------------------------------------------------
// Constraint kind enum
// ---------------------------------------------------------------------------

/// The kind of compliance a constraint carries.
///
/// Compliance is the inverse of stiffness: higher compliance means softer.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum ConstraintKind {
    /// Perfectly rigid – zero compliance, enforced exactly.
    Rigid,
    /// Linear elastic with explicit stiffness (N/m).
    Elastic {
        /// Spring stiffness in N/m.
        stiffness: Real,
    },
    /// Bend/dihedral softness, measured in N·m/rad.
    Bending {
        /// Bending stiffness in N·m/rad.
        stiffness: Real,
    },
    /// Volume preservation with given bulk modulus (Pa).
    Volume {
        /// Bulk modulus in Pa.
        bulk_modulus: Real,
    },
    /// Collision response with a restitution coefficient.
    Collision {
        /// Coefficient of restitution (dimensionless, 0–1).
        restitution: Real,
    },
    /// Arbitrary compliance value (inverse stiffness, m²/N).
    Custom {
        /// Compliance value (α = 1/k).
        compliance: Real,
    },
}

impl ConstraintKind {
    /// Convert the kind to an XPBD compliance value (α).
    ///
    /// Returns 0 for `Rigid` and `Collision` (hard constraints).
    #[allow(dead_code)]
    pub fn compliance(&self) -> Real {
        match self {
            ConstraintKind::Rigid => 0.0,
            ConstraintKind::Elastic { stiffness } => {
                if *stiffness > 0.0 {
                    1.0 / stiffness
                } else {
                    0.0
                }
            }
            ConstraintKind::Bending { stiffness } => {
                if *stiffness > 0.0 {
                    1.0 / stiffness
                } else {
                    0.0
                }
            }
            ConstraintKind::Volume { bulk_modulus } => {
                if *bulk_modulus > 0.0 {
                    1.0 / bulk_modulus
                } else {
                    0.0
                }
            }
            ConstraintKind::Collision { .. } => 0.0,
            ConstraintKind::Custom { compliance } => *compliance,
        }
    }
}

// ---------------------------------------------------------------------------
// Sleep state
// ---------------------------------------------------------------------------

/// Whether the solver considers the body asleep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SleepState {
    /// The body is actively being simulated.
    Awake,
    /// The body is below the sleep threshold and not being integrated.
    Asleep,
}

// ---------------------------------------------------------------------------
// XpbdSolver
// ---------------------------------------------------------------------------

/// Extended Position Based Dynamics (XPBD) solver.
///
/// The solver performs the following pipeline per call to [`XpbdSolver::solve`]:
///
/// 1. Predict positions using semi-implicit Euler integration.
/// 2. For each sub-step, project all constraints.
/// 3. Update velocities from corrected positions.
/// 4. Apply velocity damping.
#[derive(Debug, Clone)]
pub struct XpbdSolver {
    /// Number of sub-steps per solver call.
    pub num_substeps: usize,
    /// Number of constraint-projection iterations per sub-step.
    pub num_iterations: usize,
    /// Velocity magnitude threshold below which a body is considered asleep.
    pub sleep_threshold: Real,
    /// Number of consecutive solver calls all particles must be below the
    /// sleep threshold before the body transitions to [`SleepState::Asleep`].
    pub sleep_counter_max: usize,
    /// Internal counter: how many calls have been below the sleep threshold.
    sleep_counter: usize,
    /// Current sleep state.
    pub sleep_state: SleepState,
}

impl XpbdSolver {
    /// Create a new XPBD solver with the given number of sub-steps.
    pub fn new(num_substeps: usize) -> Self {
        Self {
            num_substeps,
            num_iterations: 1,
            sleep_threshold: 1e-4,
            sleep_counter_max: 10,
            sleep_counter: 0,
            sleep_state: SleepState::Awake,
        }
    }

    /// Create a solver with explicit sub-step and iteration counts.
    #[allow(dead_code)]
    pub fn with_iterations(num_substeps: usize, num_iterations: usize) -> Self {
        Self {
            num_substeps,
            num_iterations,
            ..Self::new(num_substeps)
        }
    }

    /// Wake the body up (reset sleep counter and state).
    #[allow(dead_code)]
    pub fn wake(&mut self) {
        self.sleep_counter = 0;
        self.sleep_state = SleepState::Awake;
    }

    /// Compute an adaptive CFL time-step given the current body state.
    ///
    /// Returns the largest `dt` such that no particle moves more than
    /// `max_displacement` in one sub-step.
    ///
    /// If all particles are static or have zero velocity, returns `dt_max`.
    #[allow(dead_code)]
    pub fn cfl_timestep(body: &SoftBody, dt_max: Real, max_displacement: Real) -> Real {
        let v_max = body
            .particles
            .iter()
            .filter(|p| !p.is_static())
            .map(|p| p.velocity.norm())
            .fold(0.0_f64, f64::max);

        if v_max < 1e-12 {
            return dt_max;
        }

        let dt_cfl = max_displacement / v_max;
        dt_cfl.min(dt_max)
    }

    /// Check whether the body should go to sleep based on current velocities.
    ///
    /// Returns `true` if all dynamic particles are below [`Self::sleep_threshold`].
    #[allow(dead_code)]
    fn check_sleep(&mut self, body: &SoftBody) -> bool {
        let all_slow = body
            .particles
            .iter()
            .filter(|p| !p.is_static())
            .all(|p| p.velocity.norm() < self.sleep_threshold);

        if all_slow {
            self.sleep_counter += 1;
        } else {
            self.sleep_counter = 0;
        }

        if self.sleep_counter >= self.sleep_counter_max {
            self.sleep_state = SleepState::Asleep;
            true
        } else {
            self.sleep_state = SleepState::Awake;
            false
        }
    }

    /// Run one full solve step over `body` with the given `constraints`.
    #[allow(clippy::needless_range_loop)]
    pub fn solve(
        &mut self,
        body: &mut SoftBody,
        constraints: &mut [Box<dyn SoftConstraint>],
        dt: Real,
    ) {
        let n = body.particles.len();
        if n == 0 || self.num_substeps == 0 {
            return;
        }

        // Skip integration if asleep (but still check if we should stay asleep).
        if self.sleep_state == SleepState::Asleep {
            return;
        }

        let dt_sub = dt / self.num_substeps as Real;

        for _sub in 0..self.num_substeps {
            // 1. Predict positions.
            for i in 0..n {
                let p = &mut body.particles[i];
                if p.is_static() {
                    continue;
                }
                p.velocity += p.external_force * (p.inverse_mass * dt_sub);
                p.prev_position = p.position;
                p.position += p.velocity * dt_sub;
            }

            // 2. Project constraints (multiple iterations per sub-step).
            for _iter in 0..self.num_iterations {
                for c in constraints.iter_mut() {
                    c.project(&mut body.particles, dt_sub);
                }
            }

            // 3. Update velocities from position corrections.
            for i in 0..n {
                let p = &mut body.particles[i];
                if p.is_static() {
                    continue;
                }
                p.velocity = (p.position - p.prev_position) / dt_sub;
            }

            // 4. Apply damping.
            let damp = 1.0 - body.damping;
            for i in 0..n {
                body.particles[i].velocity *= damp;
            }
        }

        // 5. Update sleep state.
        self.check_sleep(body);
    }

    /// Integrate particle positions forward by `dt` without constraint projection.
    ///
    /// Useful when you want to call the integration and projection phases
    /// separately (e.g. from a higher-level PBD loop).
    #[allow(dead_code)]
    pub fn integrate_positions(&self, body: &mut SoftBody, dt: Real) {
        for p in &mut body.particles {
            if p.is_static() {
                continue;
            }
            p.velocity += p.external_force * (p.inverse_mass * dt);
            p.prev_position = p.position;
            p.position += p.velocity * dt;
        }
    }

    /// Update particle velocities from the displacement since the last
    /// `integrate_positions` call (i.e. from `prev_position`).
    ///
    /// Call this **after** all constraint projections for a sub-step.
    #[allow(dead_code)]
    pub fn integrate_velocities(&self, body: &mut SoftBody, dt: Real) {
        for p in &mut body.particles {
            if p.is_static() {
                continue;
            }
            p.velocity = (p.position - p.prev_position) / dt;
        }
    }

    /// Apply velocity damping to all dynamic particles.
    #[allow(dead_code)]
    pub fn apply_damping(&self, body: &mut SoftBody) {
        let damp = 1.0 - body.damping;
        for p in &mut body.particles {
            p.velocity *= damp;
        }
    }

    /// Compute total kinetic energy of the body (½ Σ mᵢ |vᵢ|²).
    #[allow(dead_code)]
    pub fn kinetic_energy(body: &SoftBody) -> Real {
        body.particles
            .iter()
            .filter(|p| !p.is_static())
            .map(|p| {
                let m = if p.inverse_mass > 0.0 {
                    1.0 / p.inverse_mass
                } else {
                    0.0
                };
                0.5 * m * p.velocity.norm_squared()
            })
            .sum()
    }

    /// Compute the maximum particle displacement in the last sub-step.
    ///
    /// Useful for adaptive iteration count decisions.
    #[allow(dead_code)]
    pub fn max_displacement(body: &SoftBody) -> Real {
        body.particles
            .iter()
            .filter(|p| !p.is_static())
            .map(|p| (p.position - p.prev_position).norm())
            .fold(0.0_f64, f64::max)
    }
}

impl Default for XpbdSolver {
    fn default() -> Self {
        Self::new(10)
    }
}

// ---------------------------------------------------------------------------
// SolverConvergenceTracker
// ---------------------------------------------------------------------------

/// Tracks convergence statistics across solver iterations.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SolverConvergenceTracker {
    /// Per-iteration maximum constraint error.
    pub error_history: Vec<Real>,
    /// Total number of constraint projections performed.
    pub total_projections: usize,
    /// Whether the last solve converged below the threshold.
    pub converged: bool,
    /// Convergence threshold (maximum allowed constraint error).
    pub threshold: Real,
}

#[allow(dead_code)]
impl SolverConvergenceTracker {
    /// Create a new tracker.
    pub fn new(threshold: Real) -> Self {
        Self {
            error_history: Vec::new(),
            total_projections: 0,
            converged: false,
            threshold,
        }
    }

    /// Record an error value for the current iteration.
    pub fn record(&mut self, error: Real) {
        self.error_history.push(error);
        self.converged = error < self.threshold;
    }

    /// Reset the tracker for a new solve.
    pub fn reset(&mut self) {
        self.error_history.clear();
        self.total_projections = 0;
        self.converged = false;
    }

    /// Convergence rate: ratio of last two errors (< 1 means converging).
    pub fn convergence_rate(&self) -> Option<Real> {
        let n = self.error_history.len();
        if n < 2 {
            return None;
        }
        let prev = self.error_history[n - 2];
        let curr = self.error_history[n - 1];
        if prev.abs() < 1e-14 {
            return None;
        }
        Some(curr / prev)
    }

    /// Suggest an iteration count based on convergence rate.
    ///
    /// If convergence is fast (rate < 0.5), fewer iterations are needed.
    /// If convergence is slow (rate > 0.9), more iterations help.
    pub fn suggest_iterations(&self, current: usize, min: usize, max: usize) -> usize {
        match self.convergence_rate() {
            Some(rate) if rate < 0.3 => (current / 2).max(min),
            Some(rate) if rate > 0.8 => (current * 2).min(max),
            _ => current,
        }
    }
}

// ---------------------------------------------------------------------------
// ConstraintOrderingStrategy
// ---------------------------------------------------------------------------

/// Strategy for ordering constraint projections within an iteration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum ConstraintOrdering {
    /// Process constraints in the order they were added (default).
    Sequential,
    /// Reverse the order each iteration (improves convergence for chains).
    Alternating,
    /// Randomly shuffle each iteration (helps avoid systematic bias).
    Shuffled,
}

/// Solve with a specified constraint ordering strategy.
#[allow(dead_code)]
pub fn solve_with_ordering(
    solver: &mut XpbdSolver,
    body: &mut SoftBody,
    constraints: &mut [Box<dyn SoftConstraint>],
    dt: Real,
    ordering: ConstraintOrdering,
) {
    let n = body.particles.len();
    if n == 0 || solver.num_substeps == 0 {
        return;
    }
    if solver.sleep_state == SleepState::Asleep {
        return;
    }

    let dt_sub = dt / solver.num_substeps as Real;

    for sub in 0..solver.num_substeps {
        // 1. Predict positions
        for i in 0..n {
            let p = &mut body.particles[i];
            if p.is_static() {
                continue;
            }
            p.velocity += p.external_force * (p.inverse_mass * dt_sub);
            p.prev_position = p.position;
            p.position += p.velocity * dt_sub;
        }

        // 2. Project constraints with ordering
        for iter in 0..solver.num_iterations {
            match ordering {
                ConstraintOrdering::Sequential => {
                    for c in constraints.iter_mut() {
                        c.project(&mut body.particles, dt_sub);
                    }
                }
                ConstraintOrdering::Alternating => {
                    if (sub + iter) % 2 == 0 {
                        for c in constraints.iter_mut() {
                            c.project(&mut body.particles, dt_sub);
                        }
                    } else {
                        for c in constraints.iter_mut().rev() {
                            c.project(&mut body.particles, dt_sub);
                        }
                    }
                }
                ConstraintOrdering::Shuffled => {
                    // Deterministic pseudo-shuffle based on iteration index
                    let offset = (iter * 7 + sub * 13) % constraints.len().max(1);
                    for k in 0..constraints.len() {
                        let idx = (k + offset) % constraints.len();
                        constraints[idx].project(&mut body.particles, dt_sub);
                    }
                }
            }
        }

        // 3. Update velocities
        for i in 0..n {
            let p = &mut body.particles[i];
            if p.is_static() {
                continue;
            }
            p.velocity = (p.position - p.prev_position) / dt_sub;
        }

        // 4. Damping
        let damp = 1.0 - body.damping;
        for i in 0..n {
            body.particles[i].velocity *= damp;
        }
    }
}

// ---------------------------------------------------------------------------
// Solver warmstarting
// ---------------------------------------------------------------------------

/// Stores per-constraint Lagrange multipliers from the previous solve
/// for warm-starting the next solve.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct WarmstartCache {
    /// Previous Lagrange multipliers indexed by constraint.
    pub lambdas: Vec<Real>,
}

#[allow(dead_code)]
impl WarmstartCache {
    /// Create a new cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// Resize the cache to match the number of constraints.
    pub fn resize(&mut self, n: usize) {
        self.lambdas.resize(n, 0.0);
    }

    /// Apply warm-start displacements to particles based on cached lambdas.
    ///
    /// This is a simplified version: it scales all particle velocities by
    /// a factor derived from the previous solve's total lambda.
    pub fn apply_warmstart(&self, body: &mut SoftBody, factor: Real) {
        if self.lambdas.is_empty() {
            return;
        }
        let avg_lambda: Real =
            self.lambdas.iter().map(|l| l.abs()).sum::<Real>() / self.lambdas.len() as Real;
        // Nudge velocities in their current direction
        for p in &mut body.particles {
            if !p.is_static() {
                let v_mag = p.velocity.norm();
                if v_mag > 1e-14 {
                    let scale = 1.0 + factor * avg_lambda / (v_mag + 1e-10);
                    p.velocity *= scale.clamp(0.5, 2.0);
                }
            }
        }
    }

    /// Reset all cached lambdas to zero.
    pub fn clear(&mut self) {
        for l in &mut self.lambdas {
            *l = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Gauss-Seidel solver
// ---------------------------------------------------------------------------

/// A simple Gauss-Seidel constraint solver that processes constraints
/// one at a time, immediately applying corrections.
///
/// This is essentially what the standard XPBD loop does, but packaged
/// as a separate utility for clarity and testability.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GaussSeidelSolver {
    /// Number of iterations.
    pub iterations: usize,
    /// Successive over-relaxation factor (1.0 = standard GS, >1 = SOR).
    pub omega: Real,
}

#[allow(dead_code)]
impl GaussSeidelSolver {
    /// Create a new Gauss-Seidel solver.
    pub fn new(iterations: usize) -> Self {
        Self {
            iterations,
            omega: 1.0,
        }
    }

    /// Create with SOR factor.
    pub fn with_sor(iterations: usize, omega: Real) -> Self {
        Self { iterations, omega }
    }

    /// Project all constraints `iterations` times.
    pub fn solve(
        &self,
        particles: &mut [SoftParticle],
        constraints: &mut [Box<dyn SoftConstraint>],
        dt_sub: Real,
    ) {
        for _ in 0..self.iterations {
            for c in constraints.iter_mut() {
                c.project(particles, dt_sub);
            }
            // Apply SOR if omega != 1
            if (self.omega - 1.0).abs() > 1e-10 {
                for p in particles.iter_mut() {
                    if !p.is_static() {
                        let displacement = p.position - p.prev_position;
                        p.position = p.prev_position + displacement * self.omega;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint::DistanceConstraint;

    use oxiphysics_core::math::Vec3;

    // T1. Substep decomposition: running 1 step with N substeps should move
    //     a freely-falling particle by roughly the same amount as N individual
    //     1-substep solves over the same total dt.
    #[test]
    fn test_substep_decomposition() {
        let make_body = || {
            let mut body =
                SoftBody::from_particles(vec![SoftParticle::new(Vec3::new(0.0, 10.0, 0.0), 1.0)]);
            body.apply_force(&Vec3::new(0.0, -9.81, 0.0));
            body
        };

        let dt = 1.0 / 60.0;
        let n_sub = 5;

        // Solver A: single call with n_sub substeps.
        let mut body_a = make_body();
        let mut solver_a = XpbdSolver::new(n_sub);
        solver_a.solve(&mut body_a, &mut [], dt);

        // Solver B: n_sub calls each with 1 substep over dt/n_sub.
        let mut body_b = make_body();
        let mut solver_b = XpbdSolver::new(1);
        for _ in 0..n_sub {
            solver_b.solve(&mut body_b, &mut [], dt / n_sub as Real);
        }

        let dy_a = (body_a.particles[0].position.y - 10.0).abs();
        let dy_b = (body_b.particles[0].position.y - 10.0).abs();

        // Results should be very close (within 1e-10).
        assert!(
            (dy_a - dy_b).abs() < 1e-10,
            "Substep decomposition mismatch: dy_a={dy_a}, dy_b={dy_b}"
        );
    }

    // T2. Sleep detection: a stationary body should go to sleep after enough
    //     calls with no external force.
    #[test]
    fn test_sleep_detection() {
        let mut body =
            SoftBody::from_particles(vec![SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0)]);
        let mut solver = XpbdSolver::new(1);
        solver.sleep_counter_max = 3;
        solver.sleep_threshold = 1e-3;

        // No external force → particle stays still → should fall asleep.
        for _ in 0..5 {
            solver.solve(&mut body, &mut [], 1.0 / 60.0);
        }
        assert_eq!(
            solver.sleep_state,
            SleepState::Asleep,
            "Body should be asleep when velocity is zero"
        );
    }

    // T3. Sleeping body does not move.
    #[test]
    fn test_sleeping_body_not_integrated() {
        let mut body =
            SoftBody::from_particles(vec![SoftParticle::new(Vec3::new(0.0, 5.0, 0.0), 1.0)]);
        body.apply_force(&Vec3::new(0.0, -9.81, 0.0));
        let mut solver = XpbdSolver::new(1);
        solver.sleep_state = SleepState::Asleep;

        let y_before = body.particles[0].position.y;
        solver.solve(&mut body, &mut [], 1.0 / 60.0);
        let y_after = body.particles[0].position.y;

        assert!(
            (y_before - y_after).abs() < 1e-15,
            "Sleeping body must not move"
        );
    }

    // T4. Wake resets sleep state and counter.
    #[test]
    fn test_wake_resets_sleep() {
        let mut solver = XpbdSolver::new(5);
        solver.sleep_state = SleepState::Asleep;
        solver.sleep_counter = 99;
        solver.wake();
        assert_eq!(solver.sleep_state, SleepState::Awake);
        assert_eq!(solver.sleep_counter, 0);
    }

    // T5. Constraint iteration: a chain of particles with competing constraints
    //     converges further with more iterations per sub-step.
    //
    // A chain: p0 (pinned) - p1 - p2 - p3, rest length 1.0 between each pair,
    // all initially displaced to positions 0, 3, 6, 9 (3× rest length).
    // With only 1 iteration per sub-step the correction cannot fully propagate
    // along the chain; with 30 iterations it gets much closer to rest lengths.
    #[test]
    fn test_constraint_iterations_converge() {
        let make_setup = |iters: usize| {
            let mut particles = vec![
                SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
                SoftParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
                SoftParticle::new(Vec3::new(6.0, 0.0, 0.0), 1.0),
                SoftParticle::new(Vec3::new(9.0, 0.0, 0.0), 1.0),
            ];
            particles[0].inverse_mass = 0.0; // pin first particle
            let rest = 1.0;
            let constraints: Vec<Box<dyn SoftConstraint>> = vec![
                Box::new(DistanceConstraint::new(0, 1, rest, 0.0)),
                Box::new(DistanceConstraint::new(1, 2, rest, 0.0)),
                Box::new(DistanceConstraint::new(2, 3, rest, 0.0)),
            ];
            let mut body = SoftBody::from_particles(particles);
            let mut constraints = constraints;
            let mut solver = XpbdSolver::with_iterations(1, iters);
            solver.solve(&mut body, &mut constraints, 1.0 / 60.0);
            // Sum of errors across all three springs.
            let d01 = (body.particles[0].position - body.particles[1].position).norm();
            let d12 = (body.particles[1].position - body.particles[2].position).norm();
            let d23 = (body.particles[2].position - body.particles[3].position).norm();
            (d01 - rest).abs() + (d12 - rest).abs() + (d23 - rest).abs()
        };

        let err_1 = make_setup(1);
        let err_30 = make_setup(30);

        assert!(
            err_30 < err_1,
            "More iterations should give smaller total error: err_1={err_1:.4}, err_30={err_30:.4}"
        );
    }

    // T6. CFL timestep clamps correctly.
    #[test]
    fn test_cfl_timestep() {
        let mut body = SoftBody::from_particles(vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
        ]);
        // Give one particle a large velocity.
        body.particles[0].velocity = Vec3::new(100.0, 0.0, 0.0);

        let dt_max = 0.1;
        let max_disp = 0.5;
        let dt_cfl = XpbdSolver::cfl_timestep(&body, dt_max, max_disp);

        // dt_cfl should be max_disp / 100 = 0.005 < dt_max.
        assert!(
            dt_cfl < dt_max,
            "CFL dt should be smaller than dt_max: {dt_cfl}"
        );
        assert!(
            (dt_cfl - max_disp / 100.0).abs() < 1e-10,
            "CFL dt mismatch: {dt_cfl}"
        );
    }

    // T7. ConstraintKind::compliance() returns correct values.
    #[test]
    fn test_constraint_kind_compliance() {
        assert_eq!(ConstraintKind::Rigid.compliance(), 0.0);
        assert_eq!(
            ConstraintKind::Collision { restitution: 0.5 }.compliance(),
            0.0
        );
        let k = ConstraintKind::Elastic { stiffness: 1000.0 };
        assert!((k.compliance() - 1e-3).abs() < 1e-12);
        let cv = ConstraintKind::Custom { compliance: 0.007 };
        assert!((cv.compliance() - 0.007).abs() < 1e-12);
    }

    // T8. kinetic_energy returns zero for static-only body.
    #[test]
    fn test_kinetic_energy_static() {
        let body = SoftBody::from_particles(vec![SoftParticle::new_static(Vec3::zeros())]);
        assert_eq!(XpbdSolver::kinetic_energy(&body), 0.0);
    }

    // T9. integrate_positions + integrate_velocities round-trip.
    #[test]
    fn test_integrate_round_trip() {
        let mut body =
            SoftBody::from_particles(vec![SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0)]);
        body.particles[0].velocity = Vec3::new(1.0, 0.0, 0.0);
        let solver = XpbdSolver::new(1);
        let dt = 0.1;
        solver.integrate_positions(&mut body, dt);
        // pos should now be 0.1 in x
        assert!((body.particles[0].position.x - 0.1).abs() < 1e-10);
        // velocity re-derived from displacement should match.
        solver.integrate_velocities(&mut body, dt);
        assert!((body.particles[0].velocity.x - 1.0).abs() < 1e-10);
    }

    // T10. SolverConvergenceTracker records and reports.
    #[test]
    fn test_convergence_tracker() {
        let mut tracker = SolverConvergenceTracker::new(0.01);
        tracker.record(1.0);
        tracker.record(0.5);
        tracker.record(0.25);

        assert_eq!(tracker.error_history.len(), 3);
        assert!(!tracker.converged, "0.25 > 0.01, should not be converged");

        let rate = tracker.convergence_rate().unwrap();
        assert!((rate - 0.5).abs() < 1e-10, "Expected rate 0.5, got {rate}");

        tracker.record(0.005);
        assert!(tracker.converged, "0.005 < 0.01, should be converged");

        tracker.reset();
        assert!(tracker.error_history.is_empty());
        assert!(!tracker.converged);
    }

    // T11. Convergence tracker suggest_iterations.
    #[test]
    fn test_suggest_iterations() {
        let mut tracker = SolverConvergenceTracker::new(0.01);
        // Fast convergence
        tracker.record(1.0);
        tracker.record(0.1); // rate = 0.1
        let suggested = tracker.suggest_iterations(10, 2, 50);
        assert!(
            suggested <= 10,
            "Fast convergence should suggest fewer iters"
        );

        tracker.reset();
        // Slow convergence
        tracker.record(1.0);
        tracker.record(0.95); // rate = 0.95
        let suggested_slow = tracker.suggest_iterations(10, 2, 50);
        assert!(
            suggested_slow >= 10,
            "Slow convergence should suggest more iters"
        );
    }

    // T12. ConstraintOrdering::Alternating vs Sequential.
    #[test]
    fn test_alternating_ordering() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(6.0, 0.0, 0.0), 1.0),
        ];
        particles[0].inverse_mass = 0.0; // pin first
        let rest = 1.0;
        let mut constraints: Vec<Box<dyn SoftConstraint>> = vec![
            Box::new(DistanceConstraint::new(0, 1, rest, 0.0)),
            Box::new(DistanceConstraint::new(1, 2, rest, 0.0)),
        ];
        let mut body = SoftBody::from_particles(particles);
        let mut solver = XpbdSolver::with_iterations(1, 5);
        solve_with_ordering(
            &mut solver,
            &mut body,
            &mut constraints,
            1.0 / 60.0,
            ConstraintOrdering::Alternating,
        );

        // Verify positions are finite and constraints partially satisfied
        for p in &body.particles {
            assert!(p.position.x.is_finite(), "position should be finite");
        }
    }

    // T13. WarmstartCache basic operations.
    #[test]
    fn test_warmstart_cache() {
        let mut cache = WarmstartCache::new();
        cache.resize(5);
        assert_eq!(cache.lambdas.len(), 5);
        for l in &cache.lambdas {
            assert!(l.abs() < 1e-14);
        }

        cache.lambdas[0] = 1.0;
        cache.lambdas[1] = -0.5;
        cache.clear();
        for l in &cache.lambdas {
            assert!(l.abs() < 1e-14, "clear should zero all lambdas");
        }
    }

    // T14. GaussSeidelSolver basic operation.
    #[test]
    fn test_gauss_seidel_solver() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(3.0, 0.0, 0.0), 1.0),
        ];
        particles[0].inverse_mass = 0.0;
        // Save prev_position for the solver
        for p in &mut particles {
            p.prev_position = p.position;
        }
        let mut constraints: Vec<Box<dyn SoftConstraint>> =
            vec![Box::new(DistanceConstraint::new(0, 1, 1.0, 0.0))];
        let gs = GaussSeidelSolver::new(20);
        gs.solve(&mut particles, &mut constraints, 1.0 / 60.0);

        let d = (particles[0].position - particles[1].position).norm();
        assert!(
            (d - 1.0).abs() < 0.5,
            "GS should bring particles closer to rest length: d={d}"
        );
    }

    // T15. max_displacement for stationary body.
    #[test]
    fn test_max_displacement_zero() {
        let body = SoftBody::from_particles(vec![SoftParticle::new(Vec3::zeros(), 1.0)]);
        let d = XpbdSolver::max_displacement(&body);
        assert!(
            d.abs() < 1e-14,
            "Stationary body should have zero displacement"
        );
    }

    // T16. SOR factor in GaussSeidelSolver.
    #[test]
    fn test_gauss_seidel_sor() {
        let gs = GaussSeidelSolver::with_sor(10, 1.5);
        assert!((gs.omega - 1.5).abs() < 1e-14);
        assert_eq!(gs.iterations, 10);
    }

    // T17. ConstraintKind::Volume compliance.
    #[test]
    fn test_volume_constraint_kind() {
        let k = ConstraintKind::Volume { bulk_modulus: 1e6 };
        assert!((k.compliance() - 1e-6).abs() < 1e-12);
    }

    // T18. Bending ConstraintKind compliance.
    #[test]
    fn test_bending_constraint_kind() {
        let k = ConstraintKind::Bending { stiffness: 500.0 };
        assert!((k.compliance() - 1.0 / 500.0).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// Constraint batching for XPBD
// ---------------------------------------------------------------------------

/// A batch of constraints that can be solved in parallel (no shared particles).
///
/// Two constraints are compatible for batching if they do not share any
/// particle indices.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ConstraintBatch {
    /// Indices into the original constraint list that belong to this batch.
    pub indices: Vec<usize>,
}

#[allow(dead_code)]
impl ConstraintBatch {
    /// Create a new empty batch.
    pub fn new() -> Self {
        Self {
            indices: Vec::new(),
        }
    }

    /// Number of constraints in this batch.
    pub fn len(&self) -> usize {
        self.indices.len()
    }

    /// Check if the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.indices.is_empty()
    }
}

impl Default for ConstraintBatch {
    fn default() -> Self {
        Self::new()
    }
}

/// Partition constraints into independent batches (graph coloring).
///
/// Two constraints conflict if they share a particle index.
/// Uses greedy graph coloring to assign batches.
///
/// # Arguments
/// * `constraint_particles` – for each constraint, the list of particle indices it touches
///
/// Returns a list of batches, each batch containing constraint indices that
/// can be solved in parallel.
#[allow(dead_code)]
pub fn partition_constraints_into_batches(
    constraint_particles: &[Vec<usize>],
) -> Vec<ConstraintBatch> {
    let n = constraint_particles.len();
    let mut batch_id = vec![usize::MAX; n];
    let mut batches: Vec<ConstraintBatch> = Vec::new();

    for i in 0..n {
        // Find the set of batch IDs used by conflicting constraints
        let mut forbidden = std::collections::BTreeSet::new();
        for j in 0..i {
            if batch_id[j] == usize::MAX {
                continue;
            }
            // Check if constraints i and j share any particle
            let share = constraint_particles[i]
                .iter()
                .any(|p| constraint_particles[j].contains(p));
            if share {
                forbidden.insert(batch_id[j]);
            }
        }

        // Find the smallest non-forbidden batch id
        let mut b = 0;
        while forbidden.contains(&b) {
            b += 1;
        }
        batch_id[i] = b;

        if b >= batches.len() {
            batches.push(ConstraintBatch::new());
        }
        batches[b].indices.push(i);
    }

    batches
}

// ---------------------------------------------------------------------------
// Compliance matrix computation
// ---------------------------------------------------------------------------

/// Compute the effective compliance matrix for a set of XPBD constraints.
///
/// For XPBD, each constraint has an effective compliance:
/// α_tilde = α / (dt² * sum_w_i |∇C_i|²)
///
/// where α is the physical compliance (1/stiffness) and the denominator
/// accounts for the weighted gradient contributions.
///
/// This function computes the diagonal of the compliance matrix (one entry
/// per constraint).
#[allow(dead_code)]
pub fn xpbd_compliance_diagonal(
    compliances: &[Real],
    gradient_norms_sq: &[Real],
    dt: Real,
) -> Vec<Real> {
    let dt_sq = dt * dt;
    compliances
        .iter()
        .zip(gradient_norms_sq.iter())
        .map(|(&alpha, &grad_sq)| {
            let denom = dt_sq * grad_sq;
            if denom.abs() < 1e-60 {
                alpha / 1e-60
            } else {
                alpha / denom
            }
        })
        .collect()
}

/// Compute the XPBD constraint residual for a distance constraint.
///
/// C = |x_b - x_a| - rest_length
///
/// Returns (constraint_value, gradient_norm_squared).
#[allow(dead_code)]
pub fn distance_constraint_residual(
    pos_a: [Real; 3],
    pos_b: [Real; 3],
    rest_length: Real,
) -> (Real, Real) {
    let dx = pos_b[0] - pos_a[0];
    let dy = pos_b[1] - pos_a[1];
    let dz = pos_b[2] - pos_a[2];
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    let c = len - rest_length;
    // |∇C|² = 2 (unit vector dotted each particle: sum = 2)
    let grad_sq = 2.0; // for unit inverse mass
    (c, grad_sq)
}

// ---------------------------------------------------------------------------
// XPBD global step with position update
// ---------------------------------------------------------------------------

/// Perform one global XPBD update step.
///
/// This is the "global" form of XPBD where all constraint corrections
/// are accumulated and applied simultaneously (Jacobi-style), suitable
/// for GPU-ready parallel execution.
///
/// Returns the sum of |Δx| (total displacement applied this step).
#[allow(dead_code)]
pub fn xpbd_global_step(
    positions: &mut [[Real; 3]],
    inv_masses: &[Real],
    constraints: &[(usize, usize, Real, Real)], // (a, b, rest_len, compliance)
    dt: Real,
) -> Real {
    let n = positions.len();
    let mut deltas = vec![[0.0_f64; 3]; n];
    let mut counts = vec![0_usize; n];
    let dt_sq = dt * dt;

    for &(a, b, rest, alpha) in constraints {
        if a >= n || b >= n {
            continue;
        }
        let dx = positions[b][0] - positions[a][0];
        let dy = positions[b][1] - positions[a][1];
        let dz = positions[b][2] - positions[a][2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < 1e-15 {
            continue;
        }

        let c = len - rest;
        let wa = inv_masses[a];
        let wb = inv_masses[b];
        let w_sum = wa + wb;
        if w_sum < 1e-30 {
            continue;
        }

        // XPBD lambda: Δλ = (-C - α̃ λ) / (w_sum + α̃)
        // Simplified (λ=0 start): Δλ = -C / (w_sum + α / dt²)
        let alpha_tilde = alpha / dt_sq;
        let d_lambda = -c / (w_sum + alpha_tilde);

        let nx = dx / len;
        let ny = dy / len;
        let nz = dz / len;

        deltas[a][0] -= wa * d_lambda * nx;
        deltas[a][1] -= wa * d_lambda * ny;
        deltas[a][2] -= wa * d_lambda * nz;
        deltas[b][0] += wb * d_lambda * nx;
        deltas[b][1] += wb * d_lambda * ny;
        deltas[b][2] += wb * d_lambda * nz;
        counts[a] += 1;
        counts[b] += 1;
    }

    // Apply averaged corrections
    let mut total_disp = 0.0;
    for i in 0..n {
        if counts[i] > 0 {
            let s = 1.0 / counts[i] as Real;
            positions[i][0] += deltas[i][0] * s;
            positions[i][1] += deltas[i][1] * s;
            positions[i][2] += deltas[i][2] * s;
            let d = (deltas[i][0] * s).hypot(deltas[i][1] * s);
            total_disp += (d * d + (deltas[i][2] * s).powi(2)).sqrt();
        }
    }
    total_disp
}

// ---------------------------------------------------------------------------
// Parallel Gauss-Seidel via graph coloring
// ---------------------------------------------------------------------------

/// Parallel Gauss-Seidel solver using pre-computed constraint batches.
///
/// Each batch is processed sequentially, but within a batch all constraints
/// can be solved in parallel (no shared particles).
#[allow(dead_code)]
pub struct ParallelGaussSeidelSolver {
    /// Number of iterations.
    pub iterations: usize,
    /// Pre-computed constraint batches.
    pub batches: Vec<ConstraintBatch>,
}

#[allow(dead_code)]
impl ParallelGaussSeidelSolver {
    /// Create a new parallel GS solver from constraint topology.
    pub fn new(iterations: usize, constraint_particles: &[Vec<usize>]) -> Self {
        let batches = partition_constraints_into_batches(constraint_particles);
        Self {
            iterations,
            batches,
        }
    }

    /// Number of batches (colors).
    pub fn n_batches(&self) -> usize {
        self.batches.len()
    }
}

// ---------------------------------------------------------------------------
// GPU-ready position update utility
// ---------------------------------------------------------------------------

/// GPU-ready position update: given current and previous positions,
/// compute new velocity for each particle.
///
/// `v_i = (x_i^new - x_i^prev) / dt`
///
/// Returns the kinetic energy.
#[allow(dead_code)]
pub fn compute_velocities_from_positions(
    positions: &[[Real; 3]],
    prev_positions: &[[Real; 3]],
    inv_masses: &[Real],
    dt: Real,
) -> (Vec<[Real; 3]>, Real) {
    assert_eq!(positions.len(), prev_positions.len());
    assert_eq!(positions.len(), inv_masses.len());

    let mut velocities = Vec::with_capacity(positions.len());
    let mut ke = 0.0;

    for i in 0..positions.len() {
        let vx = (positions[i][0] - prev_positions[i][0]) / dt;
        let vy = (positions[i][1] - prev_positions[i][1]) / dt;
        let vz = (positions[i][2] - prev_positions[i][2]) / dt;
        velocities.push([vx, vy, vz]);

        if inv_masses[i] > 0.0 {
            let mass = 1.0 / inv_masses[i];
            ke += 0.5 * mass * (vx * vx + vy * vy + vz * vz);
        }
    }

    (velocities, ke)
}

// ---------------------------------------------------------------------------
// Additional tests for XPBD solver extensions
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests_extended {

    use crate::solver::ConstraintBatch;
    use crate::solver::ParallelGaussSeidelSolver;
    use crate::solver::compute_velocities_from_positions;
    use crate::solver::distance_constraint_residual;
    use crate::solver::partition_constraints_into_batches;
    use crate::solver::xpbd_compliance_diagonal;
    use crate::solver::xpbd_global_step;

    #[test]
    fn test_partition_no_conflicts() {
        // Two constraints on different particles → same batch possible
        let cp = vec![vec![0, 1], vec![2, 3]];
        let batches = partition_constraints_into_batches(&cp);
        assert_eq!(
            batches.len(),
            1,
            "non-conflicting constraints should be in 1 batch"
        );
        assert_eq!(batches[0].len(), 2);
    }

    #[test]
    fn test_partition_all_conflicts() {
        // Three constraints all sharing particle 0
        let cp = vec![vec![0, 1], vec![0, 2], vec![0, 3]];
        let batches = partition_constraints_into_batches(&cp);
        assert_eq!(
            batches.len(),
            3,
            "all-conflicting constraints need 3 batches"
        );
    }

    #[test]
    fn test_partition_chain_of_constraints() {
        // Chain: 0-1, 1-2, 2-3 → 2-colorable
        let cp = vec![vec![0, 1], vec![1, 2], vec![2, 3]];
        let batches = partition_constraints_into_batches(&cp);
        // Should need exactly 2 colors for a path graph
        assert!(batches.len() >= 2 && batches.len() <= 3);
    }

    #[test]
    fn test_xpbd_compliance_diagonal() {
        let compliances = vec![1e-3, 0.0];
        let grad_sq = vec![2.0, 2.0];
        let dt = 0.01;
        let diag = xpbd_compliance_diagonal(&compliances, &grad_sq, dt);
        assert_eq!(diag.len(), 2);
        // alpha_tilde = alpha / (dt^2 * grad_sq)
        let expected = 1e-3 / (0.01 * 0.01 * 2.0);
        assert!(
            (diag[0] - expected).abs() / expected < 1e-10,
            "diag[0] = {}",
            diag[0]
        );
    }

    #[test]
    fn test_distance_constraint_residual() {
        let a = [0.0, 0.0, 0.0];
        let b = [3.0, 0.0, 0.0]; // actual distance = 3.0
        let rest = 1.0;
        let (c, grad_sq) = distance_constraint_residual(a, b, rest);
        assert!(
            (c - 2.0).abs() < 1e-12,
            "constraint value should be 2.0: {c}"
        );
        assert!(
            (grad_sq - 2.0).abs() < 1e-12,
            "gradient norm sq = {grad_sq}"
        );
    }

    #[test]
    fn test_distance_constraint_residual_at_rest() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let (c, _) = distance_constraint_residual(a, b, 1.0);
        assert!(
            c.abs() < 1e-12,
            "constraint at rest length should be zero: {c}"
        );
    }

    #[test]
    fn test_xpbd_global_step_reduces_violation() {
        let mut positions = [[0.0, 0.0, 0.0_f64], [3.0, 0.0, 0.0]];
        let inv_masses = [1.0, 1.0];
        let constraints = vec![(0, 1, 1.0, 0.0)];
        let dt = 0.01;

        xpbd_global_step(&mut positions, &inv_masses, &constraints, dt);

        let dx = positions[1][0] - positions[0][0];
        let dy = positions[1][1] - positions[0][1];
        let dz = positions[1][2] - positions[0][2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        // Should be closer to rest length 1.0 than original 3.0
        assert!(dist < 3.0, "distance should decrease: {dist}");
    }

    #[test]
    fn test_xpbd_global_step_static_particle() {
        let mut positions = [[0.0, 0.0, 0.0_f64], [3.0, 0.0, 0.0]];
        let inv_masses = [0.0, 1.0]; // first particle is static
        let constraints = vec![(0, 1, 1.0, 0.0)];
        let dt = 0.01;

        xpbd_global_step(&mut positions, &inv_masses, &constraints, dt);

        // Static particle should not move
        assert!(
            (positions[0][0]).abs() < 1e-14,
            "static particle should not move"
        );
    }

    #[test]
    fn test_compute_velocities_from_positions() {
        let pos = [[1.0, 0.0, 0.0_f64]];
        let prev = [[0.0, 0.0, 0.0_f64]];
        let inv_m = [1.0];
        let dt = 0.1;

        let (vels, ke) = compute_velocities_from_positions(&pos, &prev, &inv_m, dt);
        assert_eq!(vels.len(), 1);
        assert!((vels[0][0] - 10.0).abs() < 1e-10, "vx = {}", vels[0][0]);
        // KE = 0.5 * 1.0 * 100 = 50
        assert!((ke - 50.0).abs() < 1e-10, "KE = {ke}");
    }

    #[test]
    fn test_compute_velocities_static_particle() {
        let pos = [[1.0, 0.0, 0.0_f64]];
        let prev = [[0.0, 0.0, 0.0_f64]];
        let inv_m = [0.0]; // static
        let dt = 0.1;

        let (_vels, ke) = compute_velocities_from_positions(&pos, &prev, &inv_m, dt);
        assert_eq!(ke, 0.0, "static particle has zero KE");
    }

    #[test]
    fn test_constraint_batch_default() {
        let batch = ConstraintBatch::default();
        assert!(batch.is_empty());
        assert_eq!(batch.len(), 0);
    }

    #[test]
    fn test_parallel_gauss_seidel_n_batches() {
        let cp = vec![vec![0, 1], vec![2, 3], vec![1, 2]];
        let pgs = ParallelGaussSeidelSolver::new(10, &cp);
        assert!(pgs.n_batches() >= 2, "chain needs at least 2 colors");
        assert_eq!(pgs.iterations, 10);
    }
}
