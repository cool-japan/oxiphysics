//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Models how the friction coefficient changes with temperature.
///
/// Uses a simple piecewise-linear model parameterised by a reference mu at
/// room temperature and a thermal gradient.
///
/// `mu(T) = mu_ref + k_thermal * (T - T_ref)`
///
/// Clamped to \[mu_min, mu_max\] to prevent unphysical values.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ThermalFriction {
    /// Friction coefficient at reference temperature.
    pub mu_ref: f64,
    /// Reference temperature (Kelvin or Celsius, consistent with T inputs).
    pub t_ref: f64,
    /// Rate of change of mu per degree: typically negative (mu drops with heat).
    pub k_thermal: f64,
    /// Minimum allowed friction coefficient.
    pub mu_min: f64,
    /// Maximum allowed friction coefficient.
    pub mu_max: f64,
}
#[allow(dead_code)]
impl ThermalFriction {
    /// Create a thermal friction model.
    pub fn new(mu_ref: f64, t_ref: f64, k_thermal: f64, mu_min: f64, mu_max: f64) -> Self {
        Self {
            mu_ref,
            t_ref,
            k_thermal,
            mu_min,
            mu_max,
        }
    }
    /// Evaluate the friction coefficient at temperature `t`.
    pub fn mu_at(&self, t: f64) -> f64 {
        let mu = self.mu_ref + self.k_thermal * (t - self.t_ref);
        mu.clamp(self.mu_min, self.mu_max)
    }
    /// Compute the friction force magnitude at temperature `t` for normal
    /// force `fn_`.
    pub fn friction_force(&self, fn_: f64, t: f64) -> f64 {
        self.mu_at(t) * fn_.abs()
    }
}
/// Velocity-dependent friction model with Stribeck effect.
///
/// μ(v) = μ_k + (μ_s - μ_k) * exp(-(v/v_s)^2) + μ_v * v
///
/// where:
/// - μ_s is the static friction coefficient
/// - μ_k is the kinetic friction coefficient
/// - v_s is the Stribeck velocity
/// - μ_v is the viscous friction coefficient
#[derive(Debug, Clone)]
pub struct VelocityDependentFriction {
    /// Static friction coefficient.
    pub mu_static: f64,
    /// Kinetic friction coefficient.
    pub mu_kinetic: f64,
    /// Stribeck velocity (m/s).
    pub stribeck_velocity: f64,
    /// Viscous friction coefficient (s/m).
    pub mu_viscous: f64,
}
impl VelocityDependentFriction {
    /// Create a new velocity-dependent friction model.
    pub fn new(mu_static: f64, mu_kinetic: f64, stribeck_velocity: f64, mu_viscous: f64) -> Self {
        Self {
            mu_static,
            mu_kinetic,
            stribeck_velocity,
            mu_viscous,
        }
    }
    /// Compute the effective friction coefficient at a given sliding speed.
    pub fn effective_mu(&self, speed: f64) -> f64 {
        let stribeck =
            (self.mu_static - self.mu_kinetic) * (-(speed / self.stribeck_velocity).powi(2)).exp();
        self.mu_kinetic + stribeck + self.mu_viscous * speed
    }
    /// Compute friction force from tangential velocity and normal force.
    pub fn friction_force(&self, tangential_vel: [f64; 3], normal_force: f64) -> [f64; 3] {
        let speed = length(tangential_vel);
        if speed < 1e-12 {
            return [0.0, 0.0, 0.0];
        }
        let mu_eff = self.effective_mu(speed);
        let mag = mu_eff * normal_force.abs();
        let dir = normalize(tangential_vel);
        scale(dir, -mag)
    }
}
/// A 2-DOF static friction constraint that locks the contact tangential velocity.
///
/// This is used when the contact is in the static regime.  The constraint
/// keeps the relative tangential velocity near zero by applying corrective
/// impulses in both tangential directions independently.
#[derive(Debug, Clone)]
pub struct StaticFrictionConstraint {
    /// Contact normal.
    pub normal: [f64; 3],
    /// First tangent direction.
    pub tangent1: [f64; 3],
    /// Second tangent direction.
    pub tangent2: [f64; 3],
    /// Effective mass along tangent1.
    pub eff_mass_t1: f64,
    /// Effective mass along tangent2.
    pub eff_mass_t2: f64,
    /// Accumulated impulse along tangent1 (warm-start).
    pub lambda_t1: f64,
    /// Accumulated impulse along tangent2 (warm-start).
    pub lambda_t2: f64,
    /// Maximum tangential impulse magnitude (μ * λ_n).
    pub max_tangential: f64,
}
impl StaticFrictionConstraint {
    /// Create a new static friction constraint.
    ///
    /// `inv_mass_sum` is `inv_mass_a + inv_mass_b` (linear contribution only).
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        normal: [f64; 3],
        tangent1: [f64; 3],
        tangent2: [f64; 3],
        inv_mass_sum: f64,
        mu: f64,
        normal_impulse: f64,
    ) -> Self {
        let eff_mass = if inv_mass_sum > 1e-15 {
            1.0 / inv_mass_sum
        } else {
            0.0
        };
        Self {
            normal,
            tangent1,
            tangent2,
            eff_mass_t1: eff_mass,
            eff_mass_t2: eff_mass,
            lambda_t1: 0.0,
            lambda_t2: 0.0,
            max_tangential: mu * normal_impulse.abs(),
        }
    }
    /// Solve one PGS iteration for the 2-DOF static friction constraint.
    ///
    /// Returns the delta impulse vector applied to body A (negate for body B).
    pub fn solve_velocity(&mut self, rel_vel: [f64; 3]) -> [f64; 3] {
        let cdot_t1 = dot(rel_vel, self.tangent1);
        let cdot_t2 = dot(rel_vel, self.tangent2);
        let delta1_raw = -self.eff_mass_t1 * cdot_t1;
        let delta2_raw = -self.eff_mass_t2 * cdot_t2;
        let old1 = self.lambda_t1;
        let old2 = self.lambda_t2;
        let new1 = old1 + delta1_raw;
        let new2 = old2 + delta2_raw;
        let mag = (new1 * new1 + new2 * new2).sqrt();
        let (clamped1, clamped2) = if mag > self.max_tangential && mag > 1e-14 {
            let scale_f = self.max_tangential / mag;
            (new1 * scale_f, new2 * scale_f)
        } else {
            (new1, new2)
        };
        self.lambda_t1 = clamped1;
        self.lambda_t2 = clamped2;
        let actual_delta1 = clamped1 - old1;
        let actual_delta2 = clamped2 - old2;
        add(
            scale(self.tangent1, actual_delta1),
            scale(self.tangent2, actual_delta2),
        )
    }
    /// Reset accumulated impulses (e.g. on a new time step without warm-starting).
    pub fn reset(&mut self) {
        self.lambda_t1 = 0.0;
        self.lambda_t2 = 0.0;
    }
    /// Whether the contact is currently at the static friction limit.
    pub fn is_at_limit(&self) -> bool {
        let mag = (self.lambda_t1 * self.lambda_t1 + self.lambda_t2 * self.lambda_t2).sqrt();
        mag >= self.max_tangential - 1e-10
    }
}
/// A cache of friction contact entries for all active contact pairs.
///
/// Provides O(n) lookup by `(body_a, body_b, contact_id)`.
#[derive(Debug, Clone, Default)]
pub struct FrictionContactCache {
    /// All cached entries.
    pub entries: Vec<FrictionCacheEntry>,
    /// Maximum number of entries before eviction.
    pub max_entries: usize,
}
impl FrictionContactCache {
    /// Create a new cache with the given capacity.
    pub fn new(max_entries: usize) -> Self {
        Self {
            entries: Vec::new(),
            max_entries,
        }
    }
    /// Look up a cache entry by (body_a, body_b, contact_id).
    pub fn get(&self, body_a: u32, body_b: u32, contact_id: u64) -> Option<&FrictionCacheEntry> {
        self.entries.iter().find(|e| {
            ((e.body_a == body_a && e.body_b == body_b)
                || (e.body_a == body_b && e.body_b == body_a))
                && e.contact_id == contact_id
        })
    }
    /// Insert or update a cache entry.
    pub fn insert(&mut self, entry: FrictionCacheEntry) {
        let existing = self.entries.iter_mut().find(|e| {
            ((e.body_a == entry.body_a && e.body_b == entry.body_b)
                || (e.body_a == entry.body_b && e.body_b == entry.body_a))
                && e.contact_id == entry.contact_id
        });
        if let Some(e) = existing {
            *e = entry;
        } else {
            if self.entries.len() >= self.max_entries
                && self.max_entries > 0
                && let Some(oldest_idx) = self
                    .entries
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, e)| e.age)
                    .map(|(i, _)| i)
            {
                self.entries.remove(oldest_idx);
            }
            self.entries.push(entry);
        }
    }
    /// Age all entries by one step.  Entries older than `max_age` are removed.
    pub fn age_and_evict(&mut self, max_age: u32) {
        self.entries.retain(|e| e.age < max_age);
        for e in &mut self.entries {
            e.age += 1;
        }
    }
    /// Number of cached entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Whether the cache is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Clear all entries.
    pub fn clear(&mut self) {
        self.entries.clear();
    }
}
/// Coulomb friction model with static and kinetic coefficients.
#[derive(Debug, Clone)]
pub struct CoulombFriction {
    /// Static friction coefficient μ_s.
    pub mu_static: f64,
    /// Kinetic (dynamic) friction coefficient μ_k.
    pub mu_kinetic: f64,
}
impl CoulombFriction {
    /// Create a new `CoulombFriction` with given coefficients.
    pub fn new(mu_static: f64, mu_kinetic: f64) -> Self {
        Self {
            mu_static,
            mu_kinetic,
        }
    }
    /// Returns `true` if `|F_t| > μ_s * F_n` (sliding condition).
    pub fn is_sliding(&self, tangential_force: f64, normal_force: f64) -> bool {
        tangential_force.abs() > self.mu_static * normal_force.abs()
    }
    /// Compute friction force given tangential velocity and normal force magnitude.
    ///
    /// - Static regime: force opposes tangential motion, clamped to μ_s * F_n.
    /// - Kinetic regime: μ_k * F_n in the direction opposing tangential velocity.
    pub fn friction_force(&self, tangential_vel: [f64; 3], normal_force: f64) -> [f64; 3] {
        let fn_abs = normal_force.abs();
        let speed = length(tangential_vel);
        if speed < 1e-10 {
            [0.0, 0.0, 0.0]
        } else {
            let max_static = self.mu_static * fn_abs;
            let kinetic_mag = self.mu_kinetic * fn_abs;
            let dir = normalize(tangential_vel);
            if kinetic_mag <= max_static {
                scale(dir, -kinetic_mag)
            } else {
                scale(dir, -max_static)
            }
        }
    }
    /// Compute friction impulse from relative velocity at contact.
    ///
    /// Removes normal component of `rel_vel`, then applies Coulomb friction
    /// in the tangential plane.
    pub fn friction_impulse(
        &self,
        rel_vel: [f64; 3],
        normal: [f64; 3],
        normal_impulse: f64,
    ) -> [f64; 3] {
        let vn = dot(rel_vel, normal);
        let normal_part = scale(normal, vn);
        let tangential_vel = sub(rel_vel, normal_part);
        let speed = length(tangential_vel);
        if speed < 1e-10 {
            return [0.0, 0.0, 0.0];
        }
        let dir = normalize(tangential_vel);
        let max_friction = self.mu_kinetic * normal_impulse.abs();
        scale(dir, -max_friction)
    }
}
/// Solves friction for multiple contact pairs sharing the same pair of bodies.
///
/// When two bodies have multiple contact points (e.g. a box resting on a
/// plane), the friction budgets of all contacts must be coordinated.  This
/// solver applies a per-contact sequential solve but post-processes the total
/// impulse to respect a global friction budget.
#[derive(Debug, Clone)]
pub struct PairwiseFrictionSolver {
    /// Kinetic friction coefficient.
    pub mu_kinetic: f64,
    /// Total normal impulse across all contacts in this pair.
    pub total_normal_impulse: f64,
    /// Per-contact accumulated tangential impulses.
    pub contact_impulses: Vec<[f64; 3]>,
}
impl PairwiseFrictionSolver {
    /// Create a new pairwise friction solver.
    pub fn new(mu_kinetic: f64) -> Self {
        Self {
            mu_kinetic,
            total_normal_impulse: 0.0,
            contact_impulses: Vec::new(),
        }
    }
    /// Set the number of contacts and reset impulses.
    pub fn set_contact_count(&mut self, n: usize) {
        self.contact_impulses.clear();
        self.contact_impulses.resize(n, [0.0, 0.0, 0.0]);
    }
    /// Update the total normal impulse.
    pub fn set_total_normal_impulse(&mut self, impulse: f64) {
        self.total_normal_impulse = impulse.abs();
    }
    /// Apply friction for one contact point.
    ///
    /// Returns the friction impulse to apply.  Clamped to `μ * Fn / n_contacts`.
    pub fn apply_contact(
        &mut self,
        contact_idx: usize,
        rel_vel: [f64; 3],
        normal: [f64; 3],
        inv_mass_sum: f64,
    ) -> [f64; 3] {
        let n_contacts = self.contact_impulses.len().max(1);
        let per_contact_budget = self.mu_kinetic * self.total_normal_impulse / n_contacts as f64;
        let vn = dot(rel_vel, normal);
        let normal_part = scale(normal, vn);
        let tangential_vel = sub(rel_vel, normal_part);
        let speed = length(tangential_vel);
        if speed < 1e-14 {
            return [0.0, 0.0, 0.0];
        }
        let eff_mass = if inv_mass_sum > 1e-15 {
            1.0 / inv_mass_sum
        } else {
            0.0
        };
        let raw_magnitude = speed * eff_mass;
        let clamped = raw_magnitude.min(per_contact_budget);
        let dir = normalize(tangential_vel);
        let impulse = scale(dir, -clamped);
        if contact_idx < self.contact_impulses.len() {
            self.contact_impulses[contact_idx] = add(self.contact_impulses[contact_idx], impulse);
        }
        impulse
    }
    /// Get the total accumulated tangential impulse across all contacts.
    pub fn total_tangential_impulse(&self) -> [f64; 3] {
        self.contact_impulses
            .iter()
            .fold([0.0, 0.0, 0.0], |acc, &imp| add(acc, imp))
    }
    /// Reset all accumulated impulses.
    pub fn reset(&mut self) {
        for imp in &mut self.contact_impulses {
            *imp = [0.0, 0.0, 0.0];
        }
    }
}
/// Parameters for the Stribeck friction model.
///
/// The Stribeck curve smoothly transitions from static through kinetic to
/// viscous friction as sliding speed increases.
///
/// ```text
/// mu(v) = mu_k + (mu_s - mu_k) * exp(-(v / v_s)^2)   (static peak)
///       + mu_v * v                                      (viscous term)
/// ```
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct StribeckParams {
    /// Kinetic (Coulomb) friction coefficient.
    pub mu_k: f64,
    /// Peak static friction coefficient (mu_s >= mu_k).
    pub mu_s: f64,
    /// Stribeck velocity: speed at which friction drops to ~37% of the
    /// static-kinetic difference.
    pub v_stribeck: f64,
    /// Viscous friction coefficient (force per unit speed).
    pub mu_viscous: f64,
}
#[allow(dead_code)]
impl StribeckParams {
    /// Create Stribeck parameters.
    pub fn new(mu_k: f64, mu_s: f64, v_stribeck: f64, mu_viscous: f64) -> Self {
        Self {
            mu_k,
            mu_s,
            v_stribeck: v_stribeck.max(1e-12),
            mu_viscous,
        }
    }
    /// Evaluate the Stribeck friction coefficient at sliding speed `v`.
    pub fn mu_at(&self, v: f64) -> f64 {
        let v_abs = v.abs();
        let ratio = v_abs / self.v_stribeck;
        let stribeck_term = (self.mu_s - self.mu_k) * (-ratio * ratio).exp();
        self.mu_k + stribeck_term + self.mu_viscous * v_abs
    }
    /// Compute the friction force magnitude for normal force `fn_` at speed `v`.
    pub fn friction_force(&self, fn_: f64, v: f64) -> f64 {
        self.mu_at(v) * fn_.abs()
    }
    /// Sign-correct friction force opposing motion.
    pub fn friction_force_signed(&self, fn_: f64, v: f64) -> f64 {
        if v.abs() < 1e-14 {
            0.0
        } else {
            -v.signum() * self.friction_force(fn_, v)
        }
    }
}
/// Friction warm-starting cache for iterative solvers.
///
/// Stores accumulated tangential impulses from the previous frame
/// and applies a fraction of them at the start of the current frame
/// to accelerate convergence.
#[derive(Debug, Clone)]
pub struct FrictionWarmstart {
    /// Warmstart factor in \[0, 1\]. 0 = no warmstart, 1 = full warmstart.
    pub warmstart_factor: f64,
    /// Previous frame's tangential impulses, indexed by contact pair id.
    pub cached_impulses: Vec<[f64; 3]>,
}
impl FrictionWarmstart {
    /// Create a new warmstart cache.
    pub fn new(warmstart_factor: f64) -> Self {
        Self {
            warmstart_factor: warmstart_factor.clamp(0.0, 1.0),
            cached_impulses: Vec::new(),
        }
    }
    /// Resize the cache for a given number of contacts.
    pub fn resize(&mut self, n_contacts: usize) {
        self.cached_impulses.resize(n_contacts, [0.0, 0.0, 0.0]);
    }
    /// Get the warmstart impulse for a contact.
    pub fn get_impulse(&self, contact_idx: usize) -> [f64; 3] {
        if contact_idx < self.cached_impulses.len() {
            scale(self.cached_impulses[contact_idx], self.warmstart_factor)
        } else {
            [0.0, 0.0, 0.0]
        }
    }
    /// Update the cached impulse for a contact after solving.
    pub fn update(&mut self, contact_idx: usize, impulse: [f64; 3]) {
        if contact_idx < self.cached_impulses.len() {
            self.cached_impulses[contact_idx] = impulse;
        }
    }
    /// Clear all cached impulses.
    pub fn clear(&mut self) {
        for imp in &mut self.cached_impulses {
            *imp = [0.0, 0.0, 0.0];
        }
    }
}
/// Linearized (polygonal) approximation of the Coulomb friction cone.
///
/// Approximates the circular friction cone with an n-sided polygon in the
/// tangent plane for use in LCP or MIP-based contact solvers.
#[derive(Debug, Clone)]
pub struct FrictionConeLinearization {
    /// Number of linearization facets.
    pub n_facets: usize,
    /// Tangent direction vectors in the contact plane.
    pub tangent_dirs: Vec<[f64; 3]>,
}
impl FrictionConeLinearization {
    /// Create a linearized friction cone with `n_facets` facets.
    ///
    /// Generates tangent directions uniformly distributed in the plane
    /// perpendicular to `normal`.
    pub fn new(normal: [f64; 3], n_facets: usize) -> Self {
        assert!(n_facets >= 3, "Need at least 3 facets");
        let n = normalize(normal);
        let t1 = if n[0].abs() < 0.9 {
            normalize(cross(n, [1.0, 0.0, 0.0]))
        } else {
            normalize(cross(n, [0.0, 1.0, 0.0]))
        };
        let t2 = normalize(cross(n, t1));
        let mut dirs = Vec::with_capacity(n_facets);
        for i in 0..n_facets {
            let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n_facets as f64);
            let cos_a = angle.cos();
            let sin_a = angle.sin();
            dirs.push([
                cos_a * t1[0] + sin_a * t2[0],
                cos_a * t1[1] + sin_a * t2[1],
                cos_a * t1[2] + sin_a * t2[2],
            ]);
        }
        Self {
            n_facets,
            tangent_dirs: dirs,
        }
    }
    /// Project a tangential impulse onto the linearized cone.
    ///
    /// Returns the projected impulse and the active facet index.
    pub fn project(&self, tangential_impulse: [f64; 3], max_friction: f64) -> ([f64; 3], usize) {
        let mut best_dot = f64::NEG_INFINITY;
        let mut best_idx = 0;
        for (i, dir) in self.tangent_dirs.iter().enumerate() {
            let d = dot(tangential_impulse, *dir);
            if d > best_dot {
                best_dot = d;
                best_idx = i;
            }
        }
        let mag = length(tangential_impulse);
        if mag <= max_friction {
            (tangential_impulse, best_idx)
        } else {
            let projected = scale(
                self.tangent_dirs[best_idx],
                max_friction.min(best_dot.max(0.0)),
            );
            (projected, best_idx)
        }
    }
    /// Get all facet normals (for constraint generation).
    pub fn facet_normals(&self) -> &[[f64; 3]] {
        &self.tangent_dirs
    }
}
/// Viscous (linear velocity-dependent) friction model.
#[derive(Debug, Clone)]
pub struct ViscousFriction {
    /// Viscous damping coefficient b.
    pub coefficient: f64,
}
impl ViscousFriction {
    /// Create a new `ViscousFriction`.
    pub fn new(coefficient: f64) -> Self {
        Self { coefficient }
    }
    /// Compute viscous friction force: F = -b * v.
    pub fn force(&self, velocity: [f64; 3]) -> [f64; 3] {
        scale(velocity, -self.coefficient)
    }
    /// Compute viscous friction impulse over a time step: J = -b * v * dt.
    pub fn impulse(&self, velocity: [f64; 3], dt: f64) -> [f64; 3] {
        scale(velocity, -self.coefficient * dt)
    }
}
/// Represents a Coulomb friction cone for impulse projection.
#[derive(Debug, Clone)]
pub struct FrictionCone {
    /// Friction coefficient μ.
    pub mu: f64,
    /// Contact normal, unit vector.
    pub normal: [f64; 3],
}
impl FrictionCone {
    /// Create a new `FrictionCone`.
    pub fn new(mu: f64, normal: [f64; 3]) -> Self {
        Self { mu, normal }
    }
    /// Project `impulse` so that its tangential part satisfies |J_t| ≤ μ * J_n.
    pub fn project_impulse(&self, impulse: [f64; 3], normal_impulse: f64) -> [f64; 3] {
        let jn = dot(impulse, self.normal);
        let normal_part = scale(self.normal, jn);
        let tangential = sub(impulse, normal_part);
        let max_tangential = self.mu * normal_impulse.abs();
        let t_mag = length(tangential);
        let clamped_tangential = if t_mag > max_tangential && t_mag > 1e-12 {
            scale(normalize(tangential), max_tangential)
        } else {
            tangential
        };
        add(normal_part, clamped_tangential)
    }
    /// Returns `true` if the tangential part of `impulse` is inside the cone.
    pub fn inside_cone(&self, impulse: [f64; 3], normal_impulse: f64) -> bool {
        let jn = dot(impulse, self.normal);
        let normal_part = scale(self.normal, jn);
        let tangential = sub(impulse, normal_part);
        length(tangential) <= self.mu * normal_impulse.abs() + 1e-10
    }
}
/// Anisotropic (elliptical) friction cone with different μ in each tangential direction.
#[derive(Debug, Clone)]
pub struct AnisotropicFriction {
    /// Friction coefficient along the first tangent direction.
    pub mu_x: f64,
    /// Friction coefficient along the second tangent direction.
    pub mu_y: f64,
}
impl AnisotropicFriction {
    /// Create a new `AnisotropicFriction`.
    pub fn new(mu_x: f64, mu_y: f64) -> Self {
        Self { mu_x, mu_y }
    }
    /// Project `tangential_impulse` inside the elliptical friction cone:
    ///
    /// (J_t1 / μ_x)² + (J_t2 / μ_y)² ≤ J_n²
    pub fn anisotropic_impulse(
        &self,
        tangential_impulse: [f64; 3],
        tangent1: [f64; 3],
        tangent2: [f64; 3],
        normal_impulse: f64,
    ) -> [f64; 3] {
        let jt1 = dot(tangential_impulse, tangent1);
        let jt2 = dot(tangential_impulse, tangent2);
        let jn_abs = normal_impulse.abs();
        let ellipse_val = if self.mu_x > 1e-12 && self.mu_y > 1e-12 {
            (jt1 / self.mu_x).powi(2) + (jt2 / self.mu_y).powi(2)
        } else {
            0.0
        };
        if ellipse_val <= jn_abs * jn_abs + 1e-20 {
            tangential_impulse
        } else {
            let scale_factor = jn_abs / ellipse_val.sqrt();
            let clamped_jt1 = jt1 * scale_factor;
            let clamped_jt2 = jt2 * scale_factor;
            let t1_part = scale(tangent1, clamped_jt1);
            let t2_part = scale(tangent2, clamped_jt2);
            add(t1_part, t2_part)
        }
    }
    /// Compute the effective friction coefficient in a given direction.
    ///
    /// For direction at angle θ from the first tangent axis:
    /// μ_eff = 1 / √((cos θ / μ_x)² + (sin θ / μ_y)²)
    pub fn effective_mu(&self, angle: f64) -> f64 {
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let denom = (cos_a / self.mu_x).powi(2) + (sin_a / self.mu_y).powi(2);
        if denom > 1e-30 {
            1.0 / denom.sqrt()
        } else {
            0.0
        }
    }
}
/// Anisotropic friction coefficients: different static and kinetic friction
/// along the primary tread axis versus the cross-axis (lateral).
///
/// Typical application: a tire rolling along world-X — longitudinal friction
/// (along rolling direction) and lateral friction (sideways) are independent.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct AnisotropicFrictionCoeffs {
    /// Static friction coefficient along the primary axis.
    pub mu_s_primary: f64,
    /// Kinetic friction coefficient along the primary axis.
    pub mu_k_primary: f64,
    /// Static friction coefficient along the secondary (cross) axis.
    pub mu_s_secondary: f64,
    /// Kinetic friction coefficient along the secondary (cross) axis.
    pub mu_k_secondary: f64,
}
#[allow(dead_code)]
impl AnisotropicFrictionCoeffs {
    /// Create anisotropic coefficients.
    pub fn new(
        mu_s_primary: f64,
        mu_k_primary: f64,
        mu_s_secondary: f64,
        mu_k_secondary: f64,
    ) -> Self {
        Self {
            mu_s_primary,
            mu_k_primary,
            mu_s_secondary,
            mu_k_secondary,
        }
    }
    /// Isotropic shorthand: same coefficients on both axes.
    pub fn isotropic(mu_s: f64, mu_k: f64) -> Self {
        Self::new(mu_s, mu_k, mu_s, mu_k)
    }
}
/// A contact pair that carries friction state for warm-starting.
#[derive(Debug, Clone)]
pub struct FrictionContactPair {
    /// Index of body A.
    pub body_a: u32,
    /// Index of body B.
    pub body_b: u32,
    /// Contact normal (from A toward B), unit vector.
    pub contact_normal: [f64; 3],
    /// Contact point in world space.
    pub contact_point: [f64; 3],
    /// Normal force magnitude at this contact.
    pub normal_force: f64,
    /// Accumulated tangential impulse for warm-starting.
    pub tangent_impulse: [f64; 3],
    /// Current friction state.
    pub state: FrictionState,
}
/// State of a friction contact.
#[derive(Debug, Clone, PartialEq)]
pub enum FrictionState {
    /// Contact is not sliding; static friction is active.
    Static,
    /// Contact is sliding; kinetic friction is active.
    Sliding {
        /// Sliding velocity at the contact point.
        velocity: [f64; 3],
    },
    /// Rolling contact.
    Rolling {
        /// Angular velocity of the rolling body.
        angular_velocity: [f64; 3],
    },
}
/// Impulse-based friction solver using a `CoulombFriction` model.
#[derive(Debug, Clone)]
pub struct FrictionSolver {
    /// Underlying Coulomb friction model.
    pub friction: CoulombFriction,
}
impl FrictionSolver {
    /// Create a new `FrictionSolver`.
    pub fn new(friction: CoulombFriction) -> Self {
        Self { friction }
    }
    /// Solve friction for a contact pair and return the friction impulse.
    ///
    /// `rel_vel` is the relative velocity of body A w.r.t. body B at the contact.
    pub fn solve_contact(
        &self,
        pair: &mut FrictionContactPair,
        rel_vel: [f64; 3],
        _inv_mass_a: f64,
        _inv_mass_b: f64,
        _dt: f64,
    ) -> [f64; 3] {
        self.friction
            .friction_impulse(rel_vel, pair.contact_normal, pair.normal_force)
    }
    /// Update accumulated tangential impulse (warm-start), clamped to friction cone.
    pub fn update_warm_start(&self, pair: &mut FrictionContactPair, new_impulse: [f64; 3]) {
        let acc = add(pair.tangent_impulse, new_impulse);
        let max_mag = self.friction.mu_kinetic * pair.normal_force.abs();
        let mag = length(acc);
        pair.tangent_impulse = if mag > max_mag && mag > 1e-12 {
            scale(normalize(acc), max_mag)
        } else {
            acc
        };
    }
}
/// Tracks energy dissipated by friction over simulation time.
///
/// This is useful for diagnostics and stability monitoring.  Large
/// dissipation spikes may indicate constraint instability.
#[derive(Debug, Clone)]
pub struct FrictionDissipation {
    /// Total energy dissipated so far (Joules).
    pub total_dissipation: f64,
    /// Dissipation in the current time step.
    pub current_step_dissipation: f64,
    /// Maximum dissipation in any single step.
    pub peak_dissipation: f64,
    /// Number of steps tracked.
    pub step_count: usize,
}
impl FrictionDissipation {
    /// Create a new dissipation tracker.
    pub fn new() -> Self {
        Self {
            total_dissipation: 0.0,
            current_step_dissipation: 0.0,
            peak_dissipation: 0.0,
            step_count: 0,
        }
    }
    /// Record friction work done in the current step.
    ///
    /// `friction_impulse` is the applied friction impulse vector.
    /// `tangential_vel` is the relative tangential velocity before the impulse.
    /// Work = -F · v (friction opposes motion, so dissipation is positive).
    pub fn record(&mut self, friction_impulse: [f64; 3], tangential_vel: [f64; 3]) {
        let work = -dot(friction_impulse, tangential_vel);
        let dissipation = work.max(0.0);
        self.current_step_dissipation += dissipation;
    }
    /// Advance to the next time step.
    pub fn advance_step(&mut self) {
        self.total_dissipation += self.current_step_dissipation;
        self.peak_dissipation = self.peak_dissipation.max(self.current_step_dissipation);
        self.step_count += 1;
        self.current_step_dissipation = 0.0;
    }
    /// Average dissipation per step.
    pub fn average_dissipation(&self) -> f64 {
        if self.step_count == 0 {
            0.0
        } else {
            self.total_dissipation / self.step_count as f64
        }
    }
    /// Reset all counters.
    pub fn reset(&mut self) {
        *self = FrictionDissipation::new();
    }
}
/// Physical material properties relevant to friction and collision response.
#[derive(Debug, Clone)]
pub struct FrictionMaterial {
    /// Static friction coefficient.
    pub static_friction: f64,
    /// Kinetic friction coefficient.
    pub kinetic_friction: f64,
    /// Rolling resistance coefficient.
    pub rolling_resistance: f64,
    /// Coefficient of restitution (bounciness).
    pub restitution: f64,
}
impl FrictionMaterial {
    /// Create a new `FrictionMaterial`.
    pub fn new(
        static_friction: f64,
        kinetic_friction: f64,
        rolling_resistance: f64,
        restitution: f64,
    ) -> Self {
        Self {
            static_friction,
            kinetic_friction,
            rolling_resistance,
            restitution,
        }
    }
    /// Combine two materials:
    /// - Friction: geometric mean √(μ_a * μ_b)
    /// - Restitution: min(e_a, e_b)
    pub fn combine(a: &FrictionMaterial, b: &FrictionMaterial) -> FrictionMaterial {
        FrictionMaterial {
            static_friction: (a.static_friction * b.static_friction).sqrt(),
            kinetic_friction: (a.kinetic_friction * b.kinetic_friction).sqrt(),
            rolling_resistance: (a.rolling_resistance * b.rolling_resistance).sqrt(),
            restitution: a.restitution.min(b.restitution),
        }
    }
    /// Combine two materials using arithmetic mean for friction.
    pub fn combine_arithmetic(a: &FrictionMaterial, b: &FrictionMaterial) -> FrictionMaterial {
        FrictionMaterial {
            static_friction: (a.static_friction + b.static_friction) * 0.5,
            kinetic_friction: (a.kinetic_friction + b.kinetic_friction) * 0.5,
            rolling_resistance: (a.rolling_resistance + b.rolling_resistance) * 0.5,
            restitution: a.restitution.max(b.restitution),
        }
    }
    /// Combine two materials using maximum for friction (conservative).
    pub fn combine_max(a: &FrictionMaterial, b: &FrictionMaterial) -> FrictionMaterial {
        FrictionMaterial {
            static_friction: a.static_friction.max(b.static_friction),
            kinetic_friction: a.kinetic_friction.max(b.kinetic_friction),
            rolling_resistance: a.rolling_resistance.max(b.rolling_resistance),
            restitution: a.restitution.min(b.restitution),
        }
    }
}
/// Accumulates friction forces/impulses over multiple contacts for a single body.
///
/// Useful for multi-contact scenarios where individual friction impulses
/// must be accumulated and the total clamped to the overall friction budget.
#[derive(Debug, Clone)]
pub struct FrictionForceAccumulator {
    /// Accumulated total tangential impulse.
    pub total_impulse: [f64; 3],
    /// Maximum allowed total tangential impulse magnitude.
    pub max_impulse: f64,
    /// Number of contributions accumulated.
    pub count: usize,
}
impl FrictionForceAccumulator {
    /// Create a new accumulator with a given friction budget.
    pub fn new(max_impulse: f64) -> Self {
        Self {
            total_impulse: [0.0, 0.0, 0.0],
            max_impulse,
            count: 0,
        }
    }
    /// Add a friction impulse contribution.
    pub fn accumulate(&mut self, impulse: [f64; 3]) {
        self.total_impulse = add(self.total_impulse, impulse);
        self.count += 1;
    }
    /// Get the clamped total impulse.
    pub fn clamped_impulse(&self) -> [f64; 3] {
        let mag = length(self.total_impulse);
        if mag > self.max_impulse && mag > 1e-12 {
            scale(normalize(self.total_impulse), self.max_impulse)
        } else {
            self.total_impulse
        }
    }
    /// Reset the accumulator.
    pub fn reset(&mut self) {
        self.total_impulse = [0.0, 0.0, 0.0];
        self.count = 0;
    }
    /// Check if the accumulated impulse exceeds the budget.
    pub fn is_saturated(&self) -> bool {
        length(self.total_impulse) >= self.max_impulse - 1e-10
    }
}
/// A cached entry for a friction contact pair.
///
/// Used to persist friction impulses across frames for warm-starting and
/// to carry extra state such as contact age.
#[derive(Debug, Clone)]
pub struct FrictionCacheEntry {
    /// Body A index.
    pub body_a: u32,
    /// Body B index.
    pub body_b: u32,
    /// Contact ID (e.g. feature pair hash).
    pub contact_id: u64,
    /// Accumulated tangential impulse from the previous frame.
    pub accumulated_impulse: [f64; 3],
    /// Number of consecutive frames this contact has been active.
    pub age: u32,
    /// Whether the contact was in the static friction regime last frame.
    pub was_static: bool,
}
impl FrictionCacheEntry {
    /// Create a new cache entry.
    pub fn new(body_a: u32, body_b: u32, contact_id: u64) -> Self {
        Self {
            body_a,
            body_b,
            contact_id,
            accumulated_impulse: [0.0, 0.0, 0.0],
            age: 0,
            was_static: false,
        }
    }
    /// Apply warm-start impulse scaled by `factor`.
    pub fn warm_start_impulse(&self, factor: f64) -> [f64; 3] {
        scale(self.accumulated_impulse, factor.clamp(0.0, 1.0))
    }
    /// Update the cached impulse at the end of the time step.
    pub fn update(&mut self, new_impulse: [f64; 3], is_static: bool) {
        self.accumulated_impulse = new_impulse;
        self.was_static = is_static;
        self.age += 1;
    }
    /// Decay the accumulated impulse by `factor` (e.g. 0.9 for 10% decay per frame).
    pub fn decay(&mut self, factor: f64) {
        self.accumulated_impulse = scale(self.accumulated_impulse, factor.clamp(0.0, 1.0));
    }
}
/// Evaluates anisotropic friction limits for a contact with a decomposed
/// tangential velocity.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct AnisotropicFrictionSolver {
    /// The axis-dependent friction coefficients.
    pub coeffs: AnisotropicFrictionCoeffs,
    /// Contact normal force magnitude (non-negative).
    pub normal_force: f64,
}
#[allow(dead_code)]
impl AnisotropicFrictionSolver {
    /// Create the solver.
    pub fn new(coeffs: AnisotropicFrictionCoeffs, normal_force: f64) -> Self {
        Self {
            coeffs,
            normal_force,
        }
    }
    /// Clamp the tangential impulse `(j_primary, j_secondary)` to the
    /// anisotropic friction ellipse and return whether slipping occurred.
    ///
    /// The impulse is clamped independently per axis using the kinetic
    /// coefficient (once the static limit is exceeded along either axis).
    pub fn clamp_impulse(&self, j_primary: f64, j_secondary: f64) -> ([f64; 2], bool) {
        let fn_ = self.normal_force;
        let max_p = self.coeffs.mu_s_primary * fn_;
        let max_s = self.coeffs.mu_s_secondary * fn_;
        let ep = if max_p > 1e-12 {
            j_primary / max_p
        } else {
            0.0
        };
        let es = if max_s > 1e-12 {
            j_secondary / max_s
        } else {
            0.0
        };
        if ep * ep + es * es <= 1.0 {
            return ([j_primary, j_secondary], false);
        }
        let kp = self.coeffs.mu_k_primary * fn_;
        let ks = self.coeffs.mu_k_secondary * fn_;
        let clamped = [j_primary.clamp(-kp, kp), j_secondary.clamp(-ks, ks)];
        (clamped, true)
    }
    /// Compute the friction force vector in 3D given the primary and secondary
    /// tangent axes and the decomposed tangential velocities.
    #[allow(clippy::too_many_arguments)]
    pub fn friction_force_3d(
        &self,
        t_primary: [f64; 3],
        t_secondary: [f64; 3],
        v_primary: f64,
        v_secondary: f64,
        dt: f64,
    ) -> [f64; 3] {
        let j_p = -v_primary / dt;
        let j_s = -v_secondary / dt;
        let (j_clamped, _slipping) = self.clamp_impulse(j_p, j_s);
        let fx = t_primary[0] * j_clamped[0] + t_secondary[0] * j_clamped[1];
        let fy = t_primary[1] * j_clamped[0] + t_secondary[1] * j_clamped[1];
        let fz = t_primary[2] * j_clamped[0] + t_secondary[2] * j_clamped[1];
        [fx, fy, fz]
    }
}
/// Models distributed (patch) friction by integrating Coulomb friction over
/// a circular contact area of radius `r`.
///
/// For a uniform pressure distribution over a circle the effective friction
/// force scales as `mu * fn * (2/3)` (standard result) and the spin-friction
/// torque is `mu * fn * (3*pi/16) * r` (for a Hertzian contact).  This struct
/// provides these standard factors.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct PatchFriction {
    /// Friction coefficient mu.
    pub mu: f64,
    /// Contact patch radius.
    pub radius: f64,
}
#[allow(dead_code)]
impl PatchFriction {
    /// Create a patch-friction model.
    pub fn new(mu: f64, radius: f64) -> Self {
        Self {
            mu,
            radius: radius.max(0.0),
        }
    }
    /// Effective translational friction force for normal force `fn_`.
    ///
    /// Uses the `2/3` factor from integrating uniform Coulomb friction.
    pub fn translational_force(&self, fn_: f64) -> f64 {
        (2.0 / 3.0) * self.mu * fn_.abs()
    }
    /// Spin friction torque about the contact normal for normal force `fn_`.
    ///
    /// Standard result for a circular patch with uniform pressure:
    /// `T_spin = (2/3) * mu * fn_ * r`.
    pub fn spin_torque(&self, fn_: f64) -> f64 {
        (2.0 / 3.0) * self.mu * fn_.abs() * self.radius
    }
    /// Compute the 3D friction force vector opposing the sliding velocity
    /// `v_slide` (already projected onto the contact plane).
    pub fn friction_force_3d(&self, v_slide: [f64; 3], fn_: f64) -> [f64; 3] {
        let mag =
            (v_slide[0] * v_slide[0] + v_slide[1] * v_slide[1] + v_slide[2] * v_slide[2]).sqrt();
        if mag < 1e-14 {
            return [0.0, 0.0, 0.0];
        }
        let f_mag = self.translational_force(fn_);
        [
            -v_slide[0] / mag * f_mag,
            -v_slide[1] / mag * f_mag,
            -v_slide[2] / mag * f_mag,
        ]
    }
}
/// Regularizes a friction impulse using a soft constraint approach.
///
/// Instead of a hard Coulomb cone projection, the regularizer adds a compliance
/// term `α` that allows a small amount of "sliding" even inside the static
/// friction cone.  This improves numerical stability.
///
/// The regularized normal equation is:
///   (K + α * I) * Δλ = -C_dot
///
/// where K is the effective mass and α is the compliance.
#[derive(Debug, Clone)]
pub struct FrictionRegularizer {
    /// Coulomb friction coefficient.
    pub mu: f64,
    /// Compliance term α (0 = standard hard constraint).
    pub compliance: f64,
}
impl FrictionRegularizer {
    /// Create a new regularizer.
    pub fn new(mu: f64, compliance: f64) -> Self {
        Self {
            mu,
            compliance: compliance.max(0.0),
        }
    }
    /// Compute a regularized tangential impulse.
    ///
    /// `tangential_vel` is the relative tangential velocity at the contact.
    /// `inv_eff_mass` is the inverse effective mass for the tangential DOFs.
    /// `normal_impulse` is the magnitude of the normal impulse (used for clamping).
    pub fn solve(
        &self,
        tangential_vel: [f64; 3],
        inv_eff_mass: f64,
        normal_impulse: f64,
    ) -> [f64; 3] {
        let speed = length(tangential_vel);
        if speed < 1e-14 {
            return [0.0, 0.0, 0.0];
        }
        let k_reg = inv_eff_mass + self.compliance;
        if k_reg < 1e-20 {
            return [0.0, 0.0, 0.0];
        }
        let raw_magnitude = speed / k_reg;
        let max_friction = self.mu * normal_impulse.abs();
        let clamped = raw_magnitude.min(max_friction);
        let dir = normalize(tangential_vel);
        scale(dir, -clamped)
    }
}
/// Computes the spin (torsional) friction torque about the contact normal.
///
/// This resists spinning/twisting of one body relative to another at the
/// contact.  The torque opposes the angular relative velocity about the
/// contact normal.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct SpinFriction {
    /// Maximum spin-friction torque = `mu * fn_ * r`.
    pub max_torque: f64,
    /// Regularization velocity (Stribeck-like transition from zero to full
    /// torque).  At `|omega_rel| = reg_vel`, torque = `0.63 * max_torque`.
    pub reg_vel: f64,
}
#[allow(dead_code)]
impl SpinFriction {
    /// Create a spin-friction model.
    pub fn new(max_torque: f64, reg_vel: f64) -> Self {
        Self {
            max_torque,
            reg_vel: reg_vel.max(1e-12),
        }
    }
    /// Torque magnitude for angular relative velocity `omega_rel` about the
    /// contact normal.
    pub fn torque_magnitude(&self, omega_rel: f64) -> f64 {
        let s = omega_rel.abs() / self.reg_vel;
        let tanh_s = if s > 20.0 { 1.0 } else { s.tanh() };
        self.max_torque * tanh_s
    }
    /// Signed torque opposing `omega_rel`.
    pub fn torque_signed(&self, omega_rel: f64) -> f64 {
        if omega_rel.abs() < 1e-14 {
            0.0
        } else {
            -omega_rel.signum() * self.torque_magnitude(omega_rel)
        }
    }
}
/// Rolling resistance model.
#[derive(Debug, Clone)]
pub struct RollingResistance {
    /// Rolling resistance coefficient c_rr.
    pub coefficient: f64,
}
impl RollingResistance {
    /// Create a new `RollingResistance`.
    pub fn new(coefficient: f64) -> Self {
        Self { coefficient }
    }
    /// Compute the rolling resistance torque.
    ///
    /// T = -c_rr * F_n * R * ω_hat
    pub fn torque(&self, normal_force: f64, radius: f64, angular_velocity: [f64; 3]) -> [f64; 3] {
        let omega_hat = normalize(angular_velocity);
        scale(omega_hat, -self.coefficient * normal_force.abs() * radius)
    }
}
/// Cached tangent basis for a contact normal.
///
/// Recomputing the tangent basis every frame for every contact is wasteful.
/// This cache stores the basis for re-use when the normal has not changed
/// significantly.
#[derive(Debug, Clone)]
pub struct TangentBasisCache {
    /// The normal for which the basis was computed.
    pub normal: [f64; 3],
    /// First tangent direction (perpendicular to `normal`).
    pub tangent1: [f64; 3],
    /// Second tangent direction (perpendicular to both `normal` and `tangent1`).
    pub tangent2: [f64; 3],
    /// Dot product threshold below which the cache is considered stale.
    pub freshness_threshold: f64,
}
impl TangentBasisCache {
    /// Create a new tangent basis cache for the given normal.
    pub fn new(normal: [f64; 3]) -> Self {
        let n = normalize(normal);
        let t1 = build_tangent(n);
        let t2 = normalize(cross(n, t1));
        Self {
            normal: n,
            tangent1: t1,
            tangent2: t2,
            freshness_threshold: 0.995,
        }
    }
    /// Returns `true` if the stored basis is still valid for `new_normal`.
    pub fn is_fresh(&self, new_normal: [f64; 3]) -> bool {
        let d = dot(self.normal, normalize(new_normal));
        d >= self.freshness_threshold
    }
    /// Update the cache for a new normal (no-op if the basis is still fresh).
    pub fn update(&mut self, new_normal: [f64; 3]) {
        if !self.is_fresh(new_normal) {
            *self = TangentBasisCache::new(new_normal);
        }
    }
    /// Decompose `impulse` into `(jt1, jt2)` scalar components along the tangent
    /// directions.
    pub fn decompose(&self, impulse: [f64; 3]) -> (f64, f64) {
        (dot(impulse, self.tangent1), dot(impulse, self.tangent2))
    }
    /// Reconstruct an impulse from tangential components `(jt1, jt2)`.
    pub fn reconstruct(&self, jt1: f64, jt2: f64) -> [f64; 3] {
        add(scale(self.tangent1, jt1), scale(self.tangent2, jt2))
    }
}
/// Projects a 3D tangential impulse onto the friction cone defined by `mu * fn_`.
///
/// The impulse is expressed in the contact-tangent plane (2 components).
/// If the impulse magnitude exceeds `mu * fn_`, it is scaled down to lie on the
/// cone boundary.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct ConeFrictionProjector {
    /// Friction coefficient (combined static or kinetic).
    pub mu: f64,
}
#[allow(dead_code)]
impl ConeFrictionProjector {
    /// Create the projector.
    pub fn new(mu: f64) -> Self {
        Self { mu }
    }
    /// Project `tangential_impulse` (a 2-vector in the tangent plane) onto the
    /// friction cone for normal force `fn_`.
    ///
    /// Returns the clamped impulse and a flag indicating whether slipping
    /// occurred.
    pub fn project(&self, tangential_impulse: [f64; 2], fn_: f64) -> ([f64; 2], bool) {
        let max_mag = self.mu * fn_.abs();
        let mag_sq = tangential_impulse[0] * tangential_impulse[0]
            + tangential_impulse[1] * tangential_impulse[1];
        if mag_sq <= max_mag * max_mag {
            (tangential_impulse, false)
        } else {
            let mag = mag_sq.sqrt();
            let scale = max_mag / mag;
            (
                [tangential_impulse[0] * scale, tangential_impulse[1] * scale],
                true,
            )
        }
    }
    /// Project a full 3D tangential impulse `t` (already in the tangent plane,
    /// i.e., perpendicular to the normal) onto the friction cone.
    pub fn project_3d(&self, t: [f64; 3], fn_: f64) -> ([f64; 3], bool) {
        let max_mag = self.mu * fn_.abs();
        let mag_sq = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
        if mag_sq <= max_mag * max_mag {
            (t, false)
        } else {
            let mag = mag_sq.sqrt();
            let scale = max_mag / mag;
            ([t[0] * scale, t[1] * scale, t[2] * scale], true)
        }
    }
}
