//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
use oxiphysics_core::BodyHandle;
use oxiphysics_core::math::{Real, Vec3};
use oxiphysics_rigid::RigidBodySet;

/// Per-frame statistics for the CCD solver.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct CcdStats {
    /// Number of CCD constraints processed.
    pub constraint_count: usize,
    /// Number of CCD impacts that produced a non-zero impulse.
    pub active_impacts: usize,
    /// Total number of sub-steps taken.
    pub substep_count: usize,
    /// Earliest TOI encountered this frame.
    pub earliest_toi: Real,
    /// Total impulse magnitude applied across all CCD impacts.
    pub total_impulse: Real,
}
impl CcdStats {
    /// Create a zeroed stats snapshot.
    #[allow(dead_code)]
    pub fn new() -> Self {
        CcdStats {
            constraint_count: 0,
            active_impacts: 0,
            substep_count: 0,
            earliest_toi: 1.0,
            total_impulse: 0.0,
        }
    }
    /// Record a CCD event.
    #[allow(dead_code)]
    pub fn record_event(&mut self, toi: Real, impulse: Real) {
        self.constraint_count += 1;
        if impulse > 1e-12 {
            self.active_impacts += 1;
            self.total_impulse += impulse;
        }
        if toi < self.earliest_toi {
            self.earliest_toi = toi;
        }
    }
    /// Merge another stats snapshot.
    #[allow(dead_code)]
    pub fn merge(&mut self, other: &CcdStats) {
        self.constraint_count += other.constraint_count;
        self.active_impacts += other.active_impacts;
        self.substep_count += other.substep_count;
        self.total_impulse += other.total_impulse;
        if other.earliest_toi < self.earliest_toi {
            self.earliest_toi = other.earliest_toi;
        }
    }
    /// Returns `true` if any CCD events were active this frame.
    #[allow(dead_code)]
    pub fn had_active_impacts(&self) -> bool {
        self.active_impacts > 0
    }
}
/// A priority queue of pending TOI events, sorted by ascending TOI.
///
/// Used to process CCD events in the correct temporal order within a frame.
#[derive(Debug, Clone, Default)]
#[allow(dead_code)]
pub struct ToiQueue {
    /// All pending TOI pairs (maintained in ascending-TOI order).
    pub(super) pairs: Vec<ToiPair>,
}
impl ToiQueue {
    /// Create an empty queue.
    #[allow(dead_code)]
    pub fn new() -> Self {
        ToiQueue { pairs: Vec::new() }
    }
    /// Insert a new pair and maintain sorted order.
    #[allow(dead_code)]
    pub fn push(&mut self, pair: ToiPair) {
        let pos = self.pairs.partition_point(|p| p.toi <= pair.toi);
        self.pairs.insert(pos, pair);
    }
    /// Pop the earliest unprocessed pair, or `None` if the queue is empty.
    #[allow(dead_code)]
    pub fn pop_earliest(&mut self) -> Option<ToiPair> {
        if let Some(pos) = self.pairs.iter().position(|p| !p.processed) {
            Some(self.pairs.remove(pos))
        } else {
            None
        }
    }
    /// Number of unprocessed pairs remaining.
    #[allow(dead_code)]
    pub fn pending_count(&self) -> usize {
        self.pairs.iter().filter(|p| !p.processed).count()
    }
    /// Total number of pairs (including processed).
    #[allow(dead_code)]
    pub fn total_count(&self) -> usize {
        self.pairs.len()
    }
    /// Clear all pairs.
    #[allow(dead_code)]
    pub fn clear(&mut self) {
        self.pairs.clear();
    }
    /// Mark a pair as processed without removing it.
    #[allow(dead_code)]
    pub fn mark_processed(&mut self, index: usize) {
        if index < self.pairs.len() {
            self.pairs[index].processed = true;
        }
    }
    /// Returns `true` if the queue has no unprocessed pairs.
    #[allow(dead_code)]
    pub fn is_done(&self) -> bool {
        self.pending_count() == 0
    }
}
/// Sub-step configuration for CCD.
///
/// When a CCD event is detected at TOI t*, the full time step dt is split into:
/// 1. Pre-TOI sub-step: `[0, t*] * dt`
/// 2. Impact resolution
/// 3. Post-TOI sub-step: `[t*, 1] * dt`
///
/// Multiple CCD events create additional sub-steps.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SubStepCcd {
    /// Maximum number of CCD sub-steps per frame.
    pub max_substeps: usize,
    /// Minimum sub-step fraction (skip sub-steps smaller than this).
    pub min_substep_fraction: Real,
    /// Velocity solver iterations at each impact.
    pub velocity_iterations: usize,
}
impl SubStepCcd {
    /// Create a new sub-step CCD configuration.
    #[allow(dead_code)]
    pub fn new(max_substeps: usize, velocity_iterations: usize) -> Self {
        SubStepCcd {
            max_substeps,
            min_substep_fraction: 1e-6,
            velocity_iterations,
        }
    }
    /// Decompose a list of TOI events into a sequence of sub-steps.
    ///
    /// Returns a list of `(t_start, t_end)` intervals covering \[0, 1\], with
    /// each CCD event inserted as a boundary.
    #[allow(dead_code)]
    pub fn build_substeps(&self, toi_events: &[Real]) -> Vec<(Real, Real)> {
        let mut times: Vec<Real> = vec![0.0, 1.0];
        for &t in toi_events {
            let t_clamped = t.clamp(0.0, 1.0);
            if !times.contains(&t_clamped) {
                times.push(t_clamped);
            }
        }
        times.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        times.dedup_by(|a, b| (*a - *b).abs() < 1e-14);
        let max_times = self.max_substeps + 1;
        if times.len() > max_times {
            times.truncate(max_times);
            *times.last_mut().expect("collection should not be empty") = 1.0;
        }
        times
            .windows(2)
            .filter(|w| (w[1] - w[0]) >= self.min_substep_fraction)
            .map(|w| (w[0], w[1]))
            .collect()
    }
}
/// Correct body positions to resolve penetration at the time of impact.
///
/// After advancing to the TOI, bodies may have a small residual penetration
/// due to numerical errors.  This function applies a position correction
/// (similar to Baumgarte stabilization) to push them apart.
///
/// Returns the position correction magnitude applied.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ToiPositionCorrector {
    /// Fraction of penetration to correct per step (Baumgarte factor).
    pub baumgarte: Real,
    /// Slop: penetrations smaller than this are ignored.
    pub slop: Real,
    /// Maximum correction per step.
    pub max_correction: Real,
}
impl ToiPositionCorrector {
    /// Create a new corrector with default Baumgarte settings.
    #[allow(dead_code)]
    pub fn new(baumgarte: Real, slop: Real, max_correction: Real) -> Self {
        ToiPositionCorrector {
            baumgarte,
            slop,
            max_correction,
        }
    }
    /// Compute the position correction impulse magnitude for two bodies.
    ///
    /// Uses the Baumgarte method:
    /// `correction = baumgarte * max(0, depth - slop) / dt`
    #[allow(dead_code)]
    pub fn correction_magnitude(&self, depth: Real, inv_mass_sum: Real, dt: Real) -> Real {
        if inv_mass_sum < 1e-15 || dt < 1e-15 {
            return 0.0;
        }
        let correctable = (depth - self.slop).max(0.0);
        let correction = self.baumgarte * correctable / dt;
        (correction / inv_mass_sum).min(self.max_correction)
    }
}
/// Output of the CCD broadphase: a candidate pair that may collide.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct CcdBroadphasePair {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Conservative upper bound on the TOI (from swept AABB overlap).
    pub toi_upper_bound: Real,
}
/// Solver for CCD constraints using sub-stepping to the time of impact.
///
/// # Algorithm (per constraint, in TOI order)
///
/// 1. Advance all bodies from the current sub-time to `toi * dt`.
/// 2. Apply a single velocity-level impulse to prevent penetration.
/// 3. Continue integration for the remaining `(1 - toi) * dt`.
///
/// Processing constraints in ascending TOI order ensures correctness when
/// multiple CCD events occur within the same frame.
#[derive(Debug, Clone)]
pub struct CcdConstraintSolver {
    /// Number of velocity solver iterations at each CCD impact.
    pub velocity_iterations: usize,
    /// Gravity applied during CCD sub-stepping.
    pub gravity: Vec3,
}
impl CcdConstraintSolver {
    /// Create a new CCD solver.
    pub fn new(velocity_iterations: usize, gravity: Vec3) -> Self {
        Self {
            velocity_iterations,
            gravity,
        }
    }
    /// Solve a list of CCD constraints for the given bodies over a full time step.
    ///
    /// Constraints are processed in ascending TOI order to handle multiple
    /// simultaneous CCD events correctly.
    pub fn solve(&self, constraints: &[CcdConstraint], bodies: &mut RigidBodySet, dt: f64) {
        let mut sorted: Vec<&CcdConstraint> = constraints.iter().collect();
        sorted.sort_by(|a, b| {
            a.toi
                .partial_cmp(&b.toi)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let mut time_cursor: f64 = 0.0;
        for constraint in sorted {
            let toi = constraint.toi.clamp(0.0, 1.0);
            let advance_fraction = (toi - time_cursor).max(0.0);
            if advance_fraction > 1e-12 {
                let sub_dt = advance_fraction * dt;
                self.integrate_all(bodies, sub_dt);
            }
            for _ in 0..self.velocity_iterations {
                self.apply_ccd_impulse(constraint, bodies);
            }
            time_cursor = toi;
        }
        let remaining = (1.0 - time_cursor).max(0.0);
        if remaining > 1e-12 {
            self.integrate_all(bodies, remaining * dt);
        }
    }
    /// Integrate forces and velocities for all bodies over `sub_dt`.
    fn integrate_all(&self, bodies: &mut RigidBodySet, sub_dt: f64) {
        for (_handle, body) in bodies.iter_mut() {
            body.integrate_forces(sub_dt, &self.gravity);
            body.integrate_velocity(sub_dt);
        }
    }
    /// Apply a single corrective impulse for one CCD constraint.
    ///
    /// Computes the relative normal velocity at the contact point and
    /// applies an impulse that brings it to zero (or to the restitution
    /// target if the impact velocity exceeds a threshold).
    fn apply_ccd_impulse(&self, c: &CcdConstraint, bodies: &mut RigidBodySet) {
        let (vel_a, ang_vel_a, pos_a, inv_mass_a, inv_inertia_a) = match bodies.get(c.body_a) {
            Some(b) => (
                b.velocity,
                b.angular_velocity,
                b.transform.position,
                b.inverse_mass,
                b.world_inverse_inertia,
            ),
            None => return,
        };
        let (vel_b, ang_vel_b, pos_b, inv_mass_b, inv_inertia_b) = match bodies.get(c.body_b) {
            Some(b) => (
                b.velocity,
                b.angular_velocity,
                b.transform.position,
                b.inverse_mass,
                b.world_inverse_inertia,
            ),
            None => return,
        };
        let r_a = c.contact_point - pos_a;
        let r_b = c.contact_point - pos_b;
        let n = c.contact_normal;
        let v_a_contact = vel_a + ang_vel_a.cross(&r_a);
        let v_b_contact = vel_b + ang_vel_b.cross(&r_b);
        let rel_vel_n = (v_a_contact - v_b_contact).dot(&n);
        if rel_vel_n >= 0.0 {
            return;
        }
        let rn_a = r_a.cross(&n);
        let rn_b = r_b.cross(&n);
        let k = inv_mass_a
            + inv_mass_b
            + rn_a.dot(&(inv_inertia_a * rn_a))
            + rn_b.dot(&(inv_inertia_b * rn_b));
        if k < 1e-12 {
            return;
        }
        let effective_mass = 1.0 / k;
        const RESTITUTION_THRESHOLD: f64 = 0.5;
        let restitution = if rel_vel_n.abs() > RESTITUTION_THRESHOLD {
            c.restitution
        } else {
            0.0
        };
        let lambda_n = effective_mass * (-rel_vel_n * (1.0 + restitution)).max(0.0);
        let normal_impulse = n * lambda_n;
        if let Some(body_a) = bodies.get_mut(c.body_a) {
            body_a.velocity += normal_impulse * body_a.inverse_mass;
            body_a.angular_velocity += body_a.world_inverse_inertia * r_a.cross(&normal_impulse);
        }
        if let Some(body_b) = bodies.get_mut(c.body_b) {
            body_b.velocity -= normal_impulse * body_b.inverse_mass;
            body_b.angular_velocity -= body_b.world_inverse_inertia * r_b.cross(&normal_impulse);
        }
        if c.friction > 1e-12 {
            self.apply_friction_impulse(c, bodies, lambda_n, r_a, r_b);
        }
    }
    /// Apply a friction impulse opposing the tangential relative velocity.
    fn apply_friction_impulse(
        &self,
        c: &CcdConstraint,
        bodies: &mut RigidBodySet,
        lambda_n: f64,
        r_a: Vec3,
        r_b: Vec3,
    ) {
        let (vel_a, ang_vel_a, inv_mass_a, inv_inertia_a) = match bodies.get(c.body_a) {
            Some(b) => (
                b.velocity,
                b.angular_velocity,
                b.inverse_mass,
                b.world_inverse_inertia,
            ),
            None => return,
        };
        let (vel_b, ang_vel_b, inv_mass_b, inv_inertia_b) = match bodies.get(c.body_b) {
            Some(b) => (
                b.velocity,
                b.angular_velocity,
                b.inverse_mass,
                b.world_inverse_inertia,
            ),
            None => return,
        };
        let n = c.contact_normal;
        let v_rel = (vel_a + ang_vel_a.cross(&r_a)) - (vel_b + ang_vel_b.cross(&r_b));
        let v_tangent = v_rel - n * v_rel.dot(&n);
        let tangent_speed = v_tangent.norm();
        if tangent_speed < 1e-12 {
            return;
        }
        let t = v_tangent / tangent_speed;
        let rt_a = r_a.cross(&t);
        let rt_b = r_b.cross(&t);
        let kt = inv_mass_a
            + inv_mass_b
            + rt_a.dot(&(inv_inertia_a * rt_a))
            + rt_b.dot(&(inv_inertia_b * rt_b));
        if kt < 1e-12 {
            return;
        }
        let lambda_t_raw = -tangent_speed / kt;
        let lambda_t = lambda_t_raw.clamp(-c.friction * lambda_n, c.friction * lambda_n);
        let friction_impulse = t * lambda_t;
        if let Some(body_a) = bodies.get_mut(c.body_a) {
            body_a.velocity += friction_impulse * body_a.inverse_mass;
            body_a.angular_velocity += body_a.world_inverse_inertia * r_a.cross(&friction_impulse);
        }
        if let Some(body_b) = bodies.get_mut(c.body_b) {
            body_b.velocity -= friction_impulse * body_b.inverse_mass;
            body_b.angular_velocity -= body_b.world_inverse_inertia * r_b.cross(&friction_impulse);
        }
    }
}
/// A speculative contact constraint that prevents tunneling by predicting
/// the closest approach over the next time step and generating a constraint
/// to prevent it.
///
/// Unlike a standard CCD constraint (which interrupts integration at the TOI),
/// a speculative constraint is added to the solver at the start of the step
/// and ensures the body does not tunnel through the surface.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SpeculativeContactConstraint {
    /// Body handle.
    pub body: BodyHandle,
    /// Contact normal (pointing into the surface).
    pub normal: Vec3,
    /// Predicted signed distance (negative = would penetrate).
    pub predicted_distance: Real,
    /// Contact point in world space.
    pub contact_point: Vec3,
    /// Coefficient of restitution.
    pub restitution: Real,
    /// Coefficient of friction.
    pub friction: Real,
    /// Accumulated normal impulse (for clamping).
    pub accumulated_impulse: Real,
}
impl SpeculativeContactConstraint {
    /// Create a new speculative contact constraint.
    #[allow(dead_code)]
    pub fn new(
        body: BodyHandle,
        normal: Vec3,
        predicted_distance: Real,
        contact_point: Vec3,
        restitution: Real,
        friction: Real,
    ) -> Self {
        SpeculativeContactConstraint {
            body,
            normal,
            predicted_distance,
            contact_point,
            restitution,
            friction,
            accumulated_impulse: 0.0,
        }
    }
    /// Compute the speculative velocity target along the normal.
    ///
    /// For a body that would penetrate by `predicted_distance` over `dt`,
    /// the required velocity correction along the normal is:
    ///
    /// ```text
    /// v_target = -predicted_distance / dt
    /// ```
    #[allow(dead_code)]
    pub fn target_normal_velocity(&self, dt: Real) -> Real {
        if dt.abs() < 1e-15 {
            return 0.0;
        }
        (self.predicted_distance / dt).min(0.0)
    }
    /// Returns `true` if the constraint is active (predicted penetration).
    #[allow(dead_code)]
    pub fn is_active(&self) -> bool {
        self.predicted_distance < 0.0
    }
}
/// An entry in the TOI pair queue.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct ToiPair {
    /// Handle of body A.
    pub body_a: BodyHandle,
    /// Handle of body B.
    pub body_b: BodyHandle,
    /// Time of impact in \[0, 1\].
    pub toi: Real,
    /// Whether this pair has already been processed.
    pub processed: bool,
}
impl ToiPair {
    /// Create a new TOI pair.
    #[allow(dead_code)]
    pub fn new(body_a: BodyHandle, body_b: BodyHandle, toi: Real) -> Self {
        ToiPair {
            body_a,
            body_b,
            toi,
            processed: false,
        }
    }
}
/// A CCD (Continuous Collision Detection) constraint computed from a TOI query.
///
/// Stores all data needed to advance two bodies to their exact contact moment
/// and then apply a corrective impulse before resuming integration.
#[derive(Debug, Clone)]
pub struct CcdConstraint {
    /// Handle of body A (the fast-moving body).
    pub body_a: BodyHandle,
    /// Handle of body B (typically static or slow).
    pub body_b: BodyHandle,
    /// Time of impact in \[0, 1\] (fraction of the full dt).
    pub toi: Real,
    /// Contact normal at the TOI (pointing from B toward A).
    pub contact_normal: Vec3,
    /// Contact point in world space at the TOI.
    pub contact_point: Vec3,
    /// Coefficient of restitution for the CCD impact.
    pub restitution: Real,
    /// Coefficient of friction at the CCD impact.
    pub friction: Real,
}
impl CcdConstraint {
    /// Create a new CCD constraint.
    pub fn new(
        body_a: BodyHandle,
        body_b: BodyHandle,
        toi: Real,
        contact_normal: Vec3,
        contact_point: Vec3,
        restitution: Real,
        friction: Real,
    ) -> Self {
        debug_assert!(
            (0.0..=1.0).contains(&toi),
            "TOI must be in [0,1], got {toi}"
        );
        Self {
            body_a,
            body_b,
            toi,
            contact_normal,
            contact_point,
            restitution,
            friction,
        }
    }
}
