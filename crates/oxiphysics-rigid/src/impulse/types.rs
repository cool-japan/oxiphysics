//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::*;
/// Restitution model for collision response.
#[derive(Debug, Clone, Copy)]
pub enum RestitutionModel {
    /// Newton's restitution: e = -v_rel_after / v_rel_before (along normal).
    Newton {
        /// Coefficient of restitution \[0, 1\].
        coefficient: f64,
    },
    /// Poisson's restitution: ratio of restitution impulse to compression impulse.
    Poisson {
        /// Coefficient of restitution \[0, 1\].
        coefficient: f64,
    },
    /// Velocity-dependent restitution: e decreases at low impact velocities.
    /// Prevents micro-bouncing (Baumgarte-like threshold).
    VelocityDependent {
        /// Base coefficient of restitution.
        coefficient: f64,
        /// Velocity threshold below which restitution is zeroed.
        threshold: f64,
    },
}
impl RestitutionModel {
    /// Compute the effective restitution coefficient for a given approach velocity.
    pub fn effective_restitution(&self, approach_velocity: f64) -> f64 {
        match *self {
            RestitutionModel::Newton { coefficient } => coefficient,
            RestitutionModel::Poisson { coefficient } => coefficient,
            RestitutionModel::VelocityDependent {
                coefficient,
                threshold,
            } => {
                if approach_velocity.abs() < threshold {
                    0.0
                } else {
                    coefficient
                }
            }
        }
    }
}
/// Corrects restitution impulses across multiple sub-steps.
///
/// When substepping is used, the relative velocity at a contact changes
/// between sub-steps due to ongoing collision resolution.  This corrector
/// rescales accumulated impulses to match the current effective restitution.
pub struct SubStepRestitutionCorrector {
    /// Target coefficient of restitution.
    pub restitution: f64,
    /// Velocity threshold: below this impact speed, restitution is zeroed.
    pub velocity_threshold: f64,
    /// Number of sub-steps (used to scale the per-step restitution).
    pub sub_steps: usize,
}
impl SubStepRestitutionCorrector {
    /// Create a new corrector.
    pub fn new(restitution: f64, velocity_threshold: f64, sub_steps: usize) -> Self {
        Self {
            restitution,
            velocity_threshold,
            sub_steps,
        }
    }
    /// Compute the effective per-sub-step restitution coefficient.
    ///
    /// The full-step restitution is distributed across `sub_steps` so that
    /// the combined effect approximates `restitution`.  For high impact
    /// velocities this is `1 − (1 − e)^(1/N)`.
    pub fn per_substep_restitution(&self) -> f64 {
        let n = self.sub_steps.max(1) as f64;
        1.0 - (1.0 - self.restitution.clamp(0.0, 1.0)).powf(1.0 / n)
    }
    /// Compute the impulse magnitude for one sub-step.
    ///
    /// Returns the normal impulse magnitude (or 0 if the contact is
    /// separating or below the velocity threshold).
    pub fn compute_substep_impulse(
        &self,
        a: &RigidBodyState,
        b: &RigidBodyState,
        manifold: &CollisionManifold,
    ) -> f64 {
        let v_rel = relative_velocity_at_contact(a, b, manifold.contact_point);
        let v_rel_n = dot3(v_rel, manifold.normal);
        if v_rel_n <= 0.0 {
            return 0.0;
        }
        if v_rel_n < self.velocity_threshold {
            return 0.0;
        }
        let e = self.per_substep_restitution();
        let denom = effective_mass_along_direction(a, b, manifold.contact_point, manifold.normal);
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        (1.0 + e) * v_rel_n / denom
    }
}
/// Strategy for clamping impulse magnitudes.
#[derive(Debug, Clone, Copy)]
pub enum ImpulseClampStrategy {
    /// No clamping — apply full computed impulse.
    None,
    /// Clamp normal impulse to be non-negative (no pull).
    NonNegativeNormal,
    /// Clamp total impulse magnitude to a maximum value.
    MaxMagnitude {
        /// Maximum impulse magnitude.
        max_impulse: f64,
    },
    /// Accumulated impulse clamping (for sequential solver).
    /// The delta impulse is computed so the accumulated impulse stays non-negative.
    Accumulated,
}
/// High-level impulse solver that bundles the coefficient of restitution,
/// Coulomb friction coefficient, and a velocity threshold for micro-bounce
/// suppression.
#[derive(Debug, Clone)]
pub struct ImpulseSolver {
    /// Nominal coefficient of restitution (0 = inelastic, 1 = elastic).
    pub restitution: f64,
    /// Coulomb friction coefficient μ.
    pub friction: f64,
    /// Impact velocity below which restitution is set to zero.
    pub velocity_threshold: f64,
}
impl ImpulseSolver {
    /// Create a new solver.
    pub fn new(restitution: f64, friction: f64, velocity_threshold: f64) -> Self {
        Self {
            restitution: restitution.clamp(0.0, 1.0),
            friction: friction.max(0.0),
            velocity_threshold: velocity_threshold.max(0.0),
        }
    }
    /// Compute the effective coefficient of restitution for a given closing
    /// speed.
    ///
    /// Returns `self.restitution` when `|closing_speed| >= velocity_threshold`,
    /// and `0.0` otherwise (inelastic at low speeds to suppress micro-bouncing).
    pub fn compute_restitution_coefficient(&self, closing_speed: f64) -> f64 {
        if closing_speed.abs() >= self.velocity_threshold {
            self.restitution
        } else {
            0.0
        }
    }
    /// Apply a Coulomb friction impulse at the contact point.
    ///
    /// Computes the tangential component of the relative velocity at the
    /// contact and applies a friction impulse that opposes sliding, clamped to
    /// the Coulomb cone `μ * |J_normal|`.
    ///
    /// * `a`, `b` — the two colliding bodies (mutated in place).
    /// * `manifold` — contact geometry.
    /// * `normal_impulse` — magnitude of the already-applied normal impulse
    ///   (used to scale the friction limit).
    pub fn apply_coulomb_friction(
        &self,
        a: &mut RigidBodyState,
        b: &mut RigidBodyState,
        manifold: &CollisionManifold,
        normal_impulse: f64,
    ) {
        let n = manifold.normal;
        let candidate = if n[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let t1 = {
            let c = cross3(n, candidate);
            let l = len3(c);
            if l < 1e-30 {
                return;
            }
            scale3(c, 1.0 / l)
        };
        let t2 = cross3(n, t1);
        let v_rel = relative_velocity_at_contact(a, b, manifold.contact_point);
        let v_t1 = dot3(v_rel, t1);
        let v_t2 = dot3(v_rel, t2);
        let max_friction = self.friction * normal_impulse.abs();
        if max_friction < 1e-30 {
            return;
        }
        let m_eff_t1 = effective_mass_along_direction(a, b, manifold.contact_point, t1);
        let m_eff_t2 = effective_mass_along_direction(a, b, manifold.contact_point, t2);
        let j_t1_raw = if m_eff_t1 > 1e-30 {
            -v_t1 / m_eff_t1
        } else {
            0.0
        };
        let j_t2_raw = if m_eff_t2 > 1e-30 {
            -v_t2 / m_eff_t2
        } else {
            0.0
        };
        let mag = (j_t1_raw * j_t1_raw + j_t2_raw * j_t2_raw).sqrt();
        let scale = if mag > max_friction && mag > 1e-30 {
            max_friction / mag
        } else {
            1.0
        };
        let j_t1 = j_t1_raw * scale;
        let j_t2 = j_t2_raw * scale;
        let friction_vec = add3(scale3(t1, j_t1), scale3(t2, j_t2));
        apply_impulse_b(a, manifold.contact_point, friction_vec);
        apply_impulse_a(b, manifold.contact_point, friction_vec);
    }
    /// Propagate an impulse through a serial chain of rigid bodies.
    ///
    /// Models a sequential articulated chain:
    /// `bodies[0]` receives the external impulse directly; each subsequent
    /// body receives the impulse attenuated by the mass ratio of its
    /// predecessor.  Angular components are ignored (translational chain only).
    ///
    /// * `bodies` — mutable slice of body states in chain order.
    /// * `direction` — unit vector defining the impulse direction (normalised
    ///   internally).
    /// * `impulse_magnitude` — initial impulse magnitude (N·s).
    ///
    /// Returns the impulse magnitudes applied to each link.
    pub fn solve_chain_constraint(
        &self,
        bodies: &mut [RigidBodyState],
        direction: [f64; 3],
        impulse_magnitude: f64,
    ) -> Vec<f64> {
        let dir_len = len3(direction);
        let dir = if dir_len > 1e-30 {
            scale3(direction, 1.0 / dir_len)
        } else {
            return vec![0.0; bodies.len()];
        };
        let mut magnitudes = Vec::with_capacity(bodies.len());
        let mut current_impulse = impulse_magnitude;
        for body in bodies.iter_mut() {
            if body.inv_mass < 1e-30 {
                magnitudes.push(current_impulse);
                current_impulse = 0.0;
                continue;
            }
            let dv = scale3(dir, current_impulse * body.inv_mass);
            body.velocity = add3(body.velocity, dv);
            magnitudes.push(current_impulse);
            let mass = 1.0 / body.inv_mass;
            current_impulse *= mass / (mass + mass);
        }
        magnitudes
    }
}
/// Extended rigid-body state that carries a persistent torque accumulator,
/// used for multi-step angular impulse integration.
pub struct AngularImpulseState {
    /// Underlying body state.
    pub body: RigidBodyState,
    /// Accumulated torque (N·m) applied this frame.
    pub torque: [f64; 3],
}
impl AngularImpulseState {
    /// Wrap a [`RigidBodyState`].
    pub fn new(body: RigidBodyState) -> Self {
        Self {
            body,
            torque: [0.0; 3],
        }
    }
    /// Accumulate a torque contribution.
    pub fn apply_torque(&mut self, t: [f64; 3]) {
        self.torque[0] += t[0];
        self.torque[1] += t[1];
        self.torque[2] += t[2];
    }
    /// Integrate accumulated torque into angular velocity, then reset.
    ///
    /// Δω = I⁻¹ · τ · dt
    pub fn integrate_torque(&mut self, dt: f64) {
        let delta_omega = apply_inv_inertia(self.body.inv_inertia_local, self.torque);
        self.body.angular_velocity[0] += delta_omega[0] * dt;
        self.body.angular_velocity[1] += delta_omega[1] * dt;
        self.body.angular_velocity[2] += delta_omega[2] * dt;
        self.torque = [0.0; 3];
    }
    /// Apply an angular impulse directly (bypasses torque accumulator).
    ///
    /// Δω = I⁻¹ · J_angular
    pub fn apply_angular_impulse(&mut self, j_angular: [f64; 3]) {
        let delta_omega = apply_inv_inertia(self.body.inv_inertia_local, j_angular);
        self.body.angular_velocity[0] += delta_omega[0];
        self.body.angular_velocity[1] += delta_omega[1];
        self.body.angular_velocity[2] += delta_omega[2];
    }
    /// Compute the angular kinetic energy: ½ ω · I · ω (diagonal inertia).
    pub fn angular_kinetic_energy(&self) -> f64 {
        let w = self.body.angular_velocity;
        let ii = self.body.inv_inertia_local;
        let ke: f64 = if ii[0] > 0.0 && ii[1] > 0.0 && ii[2] > 0.0 {
            0.5 * (w[0] * w[0] / ii[0] + w[1] * w[1] / ii[1] + w[2] * w[2] / ii[2])
        } else {
            0.0
        };
        ke
    }
}
/// Sequential impulse solver configuration.
pub struct SequentialImpulseSolverConfig {
    /// Number of velocity iterations.
    pub velocity_iterations: usize,
    /// Baumgarte stabilization factor (typically 0.1-0.3).
    pub baumgarte: f64,
    /// Time step for Baumgarte bias calculation.
    pub dt: f64,
    /// Whether to apply warm-starting.
    pub warm_start: bool,
}
/// Accumulated impulse for a single contact pair, cached between frames for
/// warm-starting.
#[derive(Debug, Clone)]
pub struct ImpulseCache {
    /// Key: (body_a index, body_b index)
    pub key: (usize, usize),
    /// Accumulated normal impulse.
    pub lambda_n: f64,
    /// Accumulated tangent-1 impulse.
    pub lambda_t1: f64,
    /// Accumulated tangent-2 impulse.
    pub lambda_t2: f64,
}
impl ImpulseCache {
    /// Create a zeroed cache entry for the given body pair.
    pub fn new(a: usize, b: usize) -> Self {
        Self {
            key: (a, b),
            lambda_n: 0.0,
            lambda_t1: 0.0,
            lambda_t2: 0.0,
        }
    }
    /// Scale all impulses (e.g. when timestep changes between frames).
    pub fn scale(&mut self, factor: f64) {
        self.lambda_n *= factor;
        self.lambda_t1 *= factor;
        self.lambda_t2 *= factor;
    }
}
/// Collision manifold describing the contact geometry.
pub struct CollisionManifold {
    /// Collision normal (unit vector, pointing from B to A).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// World-space contact point.
    pub contact_point: [f64; 3],
}
/// Maintains a pool of [`ImpulseCache`] entries and exposes a warm-started
/// sequential impulse solve.
pub struct WarmStartSolver {
    /// Solver configuration.
    pub config: SequentialImpulseSolverConfig,
    /// Cached impulses from the previous frame.
    pub cache: Vec<ImpulseCache>,
}
impl WarmStartSolver {
    /// Create with the given config.
    pub fn new(config: SequentialImpulseSolverConfig) -> Self {
        Self {
            config,
            cache: Vec::new(),
        }
    }
    /// Look up the cache entry for a body pair.
    fn find_cache(&self, a: usize, b: usize) -> Option<&ImpulseCache> {
        self.cache
            .iter()
            .find(|e| e.key == (a, b) || e.key == (b, a))
    }
    /// Look up the mutable cache entry, creating it if absent.
    fn get_or_insert_cache(&mut self, a: usize, b: usize) -> &mut ImpulseCache {
        let pos = self
            .cache
            .iter()
            .position(|e| e.key == (a, b) || e.key == (b, a));
        match pos {
            Some(i) => &mut self.cache[i],
            None => {
                self.cache.push(ImpulseCache::new(a, b));
                self.cache
                    .last_mut()
                    .expect("collection should not be empty")
            }
        }
    }
    /// Remove cached entries for body pairs that are no longer in contact.
    pub fn prune_cache(&mut self, active_pairs: &[(usize, usize)]) {
        self.cache.retain(|e| {
            active_pairs
                .iter()
                .any(|&(a, b)| e.key == (a, b) || e.key == (b, a))
        });
    }
    /// Run one frame of warm-started sequential impulse solving.
    ///
    /// 1. Pre-compute effective masses and velocity biases.
    /// 2. Apply warm-start impulses from the cache.
    /// 3. Iterate `config.velocity_iterations` times.
    /// 4. Update the cache with the new accumulated impulses.
    pub fn solve(&mut self, bodies: &mut [RigidBodyState], constraints: &mut [ContactConstraint]) {
        for c in constraints.iter_mut() {
            let (a_slice, b_slice) = bodies.split_at(c.body_b);
            let a = &a_slice[c.body_a];
            let b = &b_slice[0];
            c.prepare(a, b, self.config.baumgarte, self.config.dt);
        }
        if self.config.warm_start {
            for c in constraints.iter_mut() {
                if let Some(cached) = self.find_cache(c.body_a, c.body_b) {
                    c.accumulated_normal_impulse = cached.lambda_n;
                    c.accumulated_friction_impulse[0] = cached.lambda_t1;
                    c.accumulated_friction_impulse[1] = cached.lambda_t2;
                }
                let (a_slice, b_slice) = bodies.split_at_mut(c.body_b);
                let a = &mut a_slice[c.body_a];
                let b = &mut b_slice[0];
                c.warm_start(a, b);
            }
        }
        for _ in 0..self.config.velocity_iterations {
            for c in constraints.iter_mut() {
                let (a_slice, b_slice) = bodies.split_at_mut(c.body_b);
                let a = &mut a_slice[c.body_a];
                let b = &mut b_slice[0];
                c.solve_normal(a, b);
                c.solve_friction_t1(a, b);
                c.solve_friction_t2(a, b);
            }
        }
        for c in constraints.iter() {
            let entry = self.get_or_insert_cache(c.body_a, c.body_b);
            entry.lambda_n = c.accumulated_normal_impulse;
            entry.lambda_t1 = c.accumulated_friction_impulse[0];
            entry.lambda_t2 = c.accumulated_friction_impulse[1];
        }
    }
}
/// Minimal rigid-body state needed for impulse resolution.
pub struct RigidBodyState {
    /// World-space position of the centre of mass.
    pub position: [f64; 3],
    /// Linear velocity.
    pub velocity: [f64; 3],
    /// Angular velocity.
    pub angular_velocity: [f64; 3],
    /// Orientation quaternion \[x, y, z, w\].
    pub orientation: [f64; 4],
    /// Inverse mass (0 for static/kinematic bodies).
    pub inv_mass: f64,
    /// Diagonal of the *local* inverse inertia tensor \[Ixx⁻¹, Iyy⁻¹, Izz⁻¹\].
    pub inv_inertia_local: [f64; 3],
}
/// Contact constraint used by the sequential impulse solver.
pub struct ContactConstraint {
    /// Index of body A in the body array.
    pub body_a: usize,
    /// Index of body B in the body array.
    pub body_b: usize,
    /// Contact normal (from B to A).
    pub normal: [f64; 3],
    /// Contact point in world space.
    pub contact_point: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// Restitution coefficient.
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
    /// Accumulated normal impulse (for warm-starting and clamping).
    pub accumulated_normal_impulse: f64,
    /// Accumulated friction impulse (2D in tangent plane: \[t1, t2\]).
    pub accumulated_friction_impulse: [f64; 2],
    /// Effective mass along normal.
    pub effective_mass_normal: f64,
    /// Tangent directions (two orthonormal vectors in the contact plane).
    pub tangent1: [f64; 3],
    /// Effective tangent direction 2.
    pub tangent2: [f64; 3],
    /// Effective mass along tangent1.
    pub effective_mass_t1: f64,
    /// Effective mass along tangent2.
    pub effective_mass_t2: f64,
    /// Target velocity (restitution bias).
    pub velocity_bias: f64,
}
impl ContactConstraint {
    /// Create a new contact constraint.
    pub fn new(
        body_a: usize,
        body_b: usize,
        normal: [f64; 3],
        contact_point: [f64; 3],
        depth: f64,
        restitution: f64,
        friction: f64,
    ) -> Self {
        let candidate = if normal[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let tangent1 = normalize3(cross3(normal, candidate));
        let tangent2 = cross3(normal, tangent1);
        Self {
            body_a,
            body_b,
            normal,
            contact_point,
            depth,
            restitution,
            friction,
            accumulated_normal_impulse: 0.0,
            accumulated_friction_impulse: [0.0, 0.0],
            effective_mass_normal: 0.0,
            tangent1,
            tangent2,
            effective_mass_t1: 0.0,
            effective_mass_t2: 0.0,
            velocity_bias: 0.0,
        }
    }
    /// Pre-compute effective masses and velocity bias. Call once before solving.
    pub fn prepare(&mut self, a: &RigidBodyState, b: &RigidBodyState, baumgarte: f64, dt: f64) {
        self.effective_mass_normal =
            effective_mass_along_direction(a, b, self.contact_point, self.normal);
        self.effective_mass_t1 =
            effective_mass_along_direction(a, b, self.contact_point, self.tangent1);
        self.effective_mass_t2 =
            effective_mass_along_direction(a, b, self.contact_point, self.tangent2);
        let v_rel = relative_velocity_at_contact(a, b, self.contact_point);
        let v_rel_n = dot3(v_rel, self.normal);
        let restitution_bias = if v_rel_n > 1.0 {
            self.restitution * v_rel_n
        } else {
            0.0
        };
        let penetration_bias = if self.depth > 0.0 {
            baumgarte * self.depth / dt
        } else {
            0.0
        };
        self.velocity_bias = restitution_bias + penetration_bias;
    }
    /// Apply warm-starting impulses (from previous frame's accumulated impulses).
    pub fn warm_start(&self, a: &mut RigidBodyState, b: &mut RigidBodyState) {
        let normal_impulse = scale3(self.normal, self.accumulated_normal_impulse);
        let friction_impulse = add3(
            scale3(self.tangent1, self.accumulated_friction_impulse[0]),
            scale3(self.tangent2, self.accumulated_friction_impulse[1]),
        );
        let total = add3(normal_impulse, friction_impulse);
        apply_impulse_b(a, self.contact_point, total);
        apply_impulse_a(b, self.contact_point, total);
    }
    /// Solve one iteration of the normal impulse constraint.
    pub fn solve_normal(&mut self, a: &mut RigidBodyState, b: &mut RigidBodyState) {
        if self.effective_mass_normal.abs() < 1e-30 {
            return;
        }
        let v_rel = relative_velocity_at_contact(a, b, self.contact_point);
        let v_rel_n = dot3(v_rel, self.normal);
        let dj = (v_rel_n + self.velocity_bias) / self.effective_mass_normal;
        let old_accumulated = self.accumulated_normal_impulse;
        self.accumulated_normal_impulse = (old_accumulated + dj).max(0.0);
        let actual_dj = self.accumulated_normal_impulse - old_accumulated;
        if actual_dj.abs() < 1e-30 {
            return;
        }
        let impulse_vec = scale3(self.normal, actual_dj);
        apply_impulse_b(a, self.contact_point, impulse_vec);
        apply_impulse_a(b, self.contact_point, impulse_vec);
    }
    /// Solve one iteration of the friction constraint along tangent1.
    pub fn solve_friction_t1(&mut self, a: &mut RigidBodyState, b: &mut RigidBodyState) {
        if self.effective_mass_t1.abs() < 1e-30 {
            return;
        }
        let v_rel = relative_velocity_at_contact(a, b, self.contact_point);
        let v_along_t1 = dot3(v_rel, self.tangent1);
        let dj = v_along_t1 / self.effective_mass_t1;
        let max_friction = self.friction * self.accumulated_normal_impulse;
        let old = self.accumulated_friction_impulse[0];
        self.accumulated_friction_impulse[0] = (old + dj).clamp(-max_friction, max_friction);
        let actual_dj = self.accumulated_friction_impulse[0] - old;
        if actual_dj.abs() < 1e-30 {
            return;
        }
        let impulse_vec = scale3(self.tangent1, actual_dj);
        apply_impulse_a(a, self.contact_point, impulse_vec);
        apply_impulse_b(b, self.contact_point, impulse_vec);
    }
    /// Solve one iteration of the friction constraint along tangent2.
    pub fn solve_friction_t2(&mut self, a: &mut RigidBodyState, b: &mut RigidBodyState) {
        if self.effective_mass_t2.abs() < 1e-30 {
            return;
        }
        let v_rel = relative_velocity_at_contact(a, b, self.contact_point);
        let v_along_t2 = dot3(v_rel, self.tangent2);
        let dj = v_along_t2 / self.effective_mass_t2;
        let max_friction = self.friction * self.accumulated_normal_impulse;
        let old = self.accumulated_friction_impulse[1];
        self.accumulated_friction_impulse[1] = (old + dj).clamp(-max_friction, max_friction);
        let actual_dj = self.accumulated_friction_impulse[1] - old;
        if actual_dj.abs() < 1e-30 {
            return;
        }
        let impulse_vec = scale3(self.tangent2, actual_dj);
        apply_impulse_a(a, self.contact_point, impulse_vec);
        apply_impulse_b(b, self.contact_point, impulse_vec);
    }
}
/// Applies direct positional correction (pseudo-velocity) to prevent
/// penetration drift.
///
/// Unlike the velocity-level Baumgarte bias embedded in the constraint solver,
/// this corrector moves body positions directly, bypassing velocity
/// integration.  It is typically used as a post-step position-correction pass.
pub struct PositionalCorrector {
    /// Baumgarte factor β ∈ (0, 1].
    pub beta: f64,
    /// Allowed penetration slop (m); penetrations below this are ignored.
    pub slop: f64,
}
impl PositionalCorrector {
    /// Create a corrector with the given beta and slop.
    pub fn new(beta: f64, slop: f64) -> Self {
        Self { beta, slop }
    }
    /// Apply positional correction to bodies `a` and `b` for a single contact.
    ///
    /// Moves each body away from the penetration by a fraction `beta` of the
    /// excess overlap (overlap minus slop), weighted by their inverse masses.
    ///
    /// # Arguments
    /// * `a` — mutable reference to body A.
    /// * `b` — mutable reference to body B.
    /// * `normal` — contact normal pointing from B to A.
    /// * `depth` — penetration depth (positive = overlap).
    pub fn correct(
        &self,
        a: &mut RigidBodyState,
        b: &mut RigidBodyState,
        normal: [f64; 3],
        depth: f64,
    ) {
        let excess = (depth - self.slop).max(0.0);
        if excess < 1e-12 {
            return;
        }
        let denom = a.inv_mass + b.inv_mass;
        if denom < 1e-30 {
            return;
        }
        let correction_mag = self.beta * excess / denom;
        let correction = scale3(normal, correction_mag);
        a.position = add3(a.position, scale3(correction, a.inv_mass));
        b.position = sub3(b.position, scale3(correction, b.inv_mass));
    }
}
