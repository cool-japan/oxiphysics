//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    CollisionManifold, ContactConstraint, ImpulseClampStrategy, RestitutionModel, RigidBodyState,
    SequentialImpulseSolverConfig,
};

/// Cross product of two 3-vectors.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Dot product of two 3-vectors.
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Scale a 3-vector by a scalar.
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Add two 3-vectors.
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors (a − b).
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Euclidean length of a 3-vector.
pub fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}
/// Normalize a 3-vector. Returns zero vector if length is negligible.
pub(super) fn normalize3(a: [f64; 3]) -> [f64; 3] {
    let l = len3(a);
    if l < 1e-30 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(a, 1.0 / l)
    }
}
pub(super) fn apply_inv_inertia(inv_inertia_local: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [
        inv_inertia_local[0] * v[0],
        inv_inertia_local[1] * v[1],
        inv_inertia_local[2] * v[2],
    ]
}
/// Clamp a normal impulse magnitude according to the strategy.
pub(super) fn clamp_normal_impulse(j: f64, strategy: ImpulseClampStrategy) -> f64 {
    match strategy {
        ImpulseClampStrategy::None => j,
        ImpulseClampStrategy::NonNegativeNormal => j.max(0.0),
        ImpulseClampStrategy::MaxMagnitude { max_impulse } => j.clamp(-max_impulse, max_impulse),
        ImpulseClampStrategy::Accumulated => j,
    }
}
/// Clamp a friction impulse vector to the Coulomb cone.
pub(super) fn clamp_friction_to_cone(friction_impulse: [f64; 3], max_friction: f64) -> [f64; 3] {
    let mag = len3(friction_impulse);
    if mag <= max_friction {
        friction_impulse
    } else {
        scale3(friction_impulse, max_friction / mag)
    }
}
/// Build a polygonal approximation to the friction cone using `n` directions.
///
/// Returns unit tangent directions in the contact plane.
/// `normal` must be a unit vector.
pub(super) fn build_friction_cone_directions(normal: [f64; 3], n: usize) -> Vec<[f64; 3]> {
    let candidate = if normal[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let t1 = normalize3(cross3(normal, candidate));
    let t2 = cross3(normal, t1);
    let mut directions = Vec::with_capacity(n);
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
        let (s, c) = angle.sin_cos();
        directions.push(add3(scale3(t1, c), scale3(t2, s)));
    }
    directions
}
/// Compute the relative velocity at the contact point.
///
/// v_rel = (v_a + ω_a × r_a) − (v_b + ω_b × r_b)
/// where r_x = contact − position_x
pub fn relative_velocity_at_contact(
    a: &RigidBodyState,
    b: &RigidBodyState,
    contact: [f64; 3],
) -> [f64; 3] {
    let ra = sub3(contact, a.position);
    let rb = sub3(contact, b.position);
    let va_contact = add3(a.velocity, cross3(a.angular_velocity, ra));
    let vb_contact = add3(b.velocity, cross3(b.angular_velocity, rb));
    sub3(va_contact, vb_contact)
}
/// Compute the effective mass along a direction at a contact point.
pub(super) fn effective_mass_along_direction(
    a: &RigidBodyState,
    b: &RigidBodyState,
    contact: [f64; 3],
    direction: [f64; 3],
) -> f64 {
    let ra = sub3(contact, a.position);
    let rb = sub3(contact, b.position);
    let ra_cross_d = cross3(ra, direction);
    let rb_cross_d = cross3(rb, direction);
    let ang_a = dot3(
        ra_cross_d,
        apply_inv_inertia(a.inv_inertia_local, ra_cross_d),
    );
    let ang_b = dot3(
        rb_cross_d,
        apply_inv_inertia(b.inv_inertia_local, rb_cross_d),
    );
    a.inv_mass + b.inv_mass + ang_a + ang_b
}
/// Compute the scalar normal impulse magnitude.
///
/// j = −(1 + e) · v_rel·n
///     ────────────────────────────────────────────────────────────────────
///     1/mA + 1/mB + (rA×n)·IA⁻¹·(rA×n) + (rB×n)·IB⁻¹·(rB×n)
pub fn compute_impulse(
    a: &RigidBodyState,
    b: &RigidBodyState,
    manifold: &CollisionManifold,
    restitution: f64,
) -> f64 {
    let n = manifold.normal;
    let contact = manifold.contact_point;
    let v_rel = relative_velocity_at_contact(a, b, contact);
    let v_rel_n = dot3(v_rel, n);
    if v_rel_n <= 0.0 {
        return 0.0;
    }
    let denom = effective_mass_along_direction(a, b, contact, n);
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    (1.0 + restitution) * v_rel_n / denom
}
/// Compute normal impulse with a specific restitution model.
pub fn compute_impulse_with_model(
    a: &RigidBodyState,
    b: &RigidBodyState,
    manifold: &CollisionManifold,
    model: RestitutionModel,
) -> f64 {
    let n = manifold.normal;
    let contact = manifold.contact_point;
    let v_rel = relative_velocity_at_contact(a, b, contact);
    let v_rel_n = dot3(v_rel, n);
    if v_rel_n <= 0.0 {
        return 0.0;
    }
    let e = model.effective_restitution(v_rel_n);
    let denom = effective_mass_along_direction(a, b, contact, n);
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    (1.0 + e) * v_rel_n / denom
}
/// Apply a normal impulse to body A (adds to velocity and angular velocity).
///
/// Δv_A  = +j·n / mA
/// Δω_A  = +IA⁻¹ · (rA × (j·n))
pub fn apply_impulse_a(state: &mut RigidBodyState, contact: [f64; 3], impulse_vec: [f64; 3]) {
    let ra = sub3(contact, state.position);
    state.velocity = add3(state.velocity, scale3(impulse_vec, state.inv_mass));
    let torque_arm = cross3(ra, impulse_vec);
    let delta_omega = apply_inv_inertia(state.inv_inertia_local, torque_arm);
    state.angular_velocity = add3(state.angular_velocity, delta_omega);
}
/// Apply a normal impulse to body B (subtracts — Newton's 3rd law).
///
/// Δv_B  = −j·n / mB
/// Δω_B  = −IB⁻¹ · (rB × (j·n))
pub fn apply_impulse_b(state: &mut RigidBodyState, contact: [f64; 3], impulse_vec: [f64; 3]) {
    let rb = sub3(contact, state.position);
    state.velocity = sub3(state.velocity, scale3(impulse_vec, state.inv_mass));
    let torque_arm = cross3(rb, impulse_vec);
    let delta_omega = apply_inv_inertia(state.inv_inertia_local, torque_arm);
    state.angular_velocity = sub3(state.angular_velocity, delta_omega);
}
/// Compute the tangential (friction) impulse vector using Coulomb's law.
///
/// * Computes the relative velocity at the contact and removes its normal component
///   to obtain the tangential relative velocity.
/// * If tangential speed is negligible, returns zero.
/// * Otherwise, applies a tangential impulse that opposes relative sliding,
///   clamped to μ · |j_n| (Coulomb cone).
pub fn friction_impulse(
    a: &RigidBodyState,
    b: &RigidBodyState,
    manifold: &CollisionManifold,
    j_n: f64,
    friction: f64,
) -> [f64; 3] {
    let n = manifold.normal;
    let contact = manifold.contact_point;
    let v_rel = relative_velocity_at_contact(a, b, contact);
    let v_rel_n_scalar = dot3(v_rel, n);
    let v_tangential = sub3(v_rel, scale3(n, v_rel_n_scalar));
    let tan_speed = len3(v_tangential);
    if tan_speed < 1e-12 {
        return [0.0, 0.0, 0.0];
    }
    let t_hat = scale3(v_tangential, -1.0 / tan_speed);
    let denom = effective_mass_along_direction(a, b, contact, t_hat);
    if denom.abs() < 1e-30 {
        return [0.0, 0.0, 0.0];
    }
    let j_t_desired = tan_speed / denom;
    let j_t = j_t_desired.min(friction * j_n.abs());
    scale3(t_hat, j_t)
}
/// Compute friction impulse using a pyramidal approximation of the friction cone.
///
/// Uses `n_dirs` directions in the contact tangent plane.
pub fn friction_impulse_pyramidal(
    a: &RigidBodyState,
    b: &RigidBodyState,
    manifold: &CollisionManifold,
    j_n: f64,
    friction: f64,
    n_dirs: usize,
) -> [f64; 3] {
    let n = manifold.normal;
    let contact = manifold.contact_point;
    let v_rel = relative_velocity_at_contact(a, b, contact);
    let v_rel_n_scalar = dot3(v_rel, n);
    let v_tangential = sub3(v_rel, scale3(n, v_rel_n_scalar));
    let tan_speed = len3(v_tangential);
    if tan_speed < 1e-12 {
        return [0.0, 0.0, 0.0];
    }
    let dirs = build_friction_cone_directions(n, n_dirs);
    let max_friction_mag = friction * j_n.abs();
    let mut best_impulse = [0.0; 3];
    let mut best_score = f64::NEG_INFINITY;
    for dir in &dirs {
        let v_along = dot3(v_tangential, *dir);
        if v_along.abs() < 1e-15 {
            continue;
        }
        let denom = effective_mass_along_direction(a, b, contact, *dir);
        if denom.abs() < 1e-30 {
            continue;
        }
        let j_dir = (-v_along / denom).clamp(-max_friction_mag, max_friction_mag);
        let candidate = scale3(*dir, j_dir);
        let score = -dot3(candidate, v_tangential);
        if score > best_score {
            best_score = score;
            best_impulse = candidate;
        }
    }
    clamp_friction_to_cone(best_impulse, max_friction_mag)
}
/// Resolve a collision between two rigid bodies.
///
/// 1. Computes the normal impulse magnitude `j`.
/// 2. Applies the normal impulse to both bodies (Newton's 3rd law).
/// 3. Computes a Coulomb friction impulse and applies it too.
pub fn resolve_collision(
    a: &mut RigidBodyState,
    b: &mut RigidBodyState,
    manifold: &CollisionManifold,
    restitution: f64,
    friction: f64,
) {
    let j = compute_impulse(a, b, manifold, restitution);
    if j.abs() < 1e-30 {
        return;
    }
    let n = manifold.normal;
    let contact = manifold.contact_point;
    let impulse_vec = scale3(n, j);
    apply_impulse_b(a, contact, impulse_vec);
    apply_impulse_a(b, contact, impulse_vec);
    if friction > 0.0 {
        let jf = friction_impulse(a, b, manifold, j, friction);
        apply_impulse_a(a, contact, jf);
        apply_impulse_b(b, contact, jf);
    }
}
/// Resolve using a specific restitution model.
pub fn resolve_collision_with_model(
    a: &mut RigidBodyState,
    b: &mut RigidBodyState,
    manifold: &CollisionManifold,
    model: RestitutionModel,
    friction: f64,
) {
    let j = compute_impulse_with_model(a, b, manifold, model);
    if j.abs() < 1e-30 {
        return;
    }
    let n = manifold.normal;
    let contact = manifold.contact_point;
    let impulse_vec = scale3(n, j);
    apply_impulse_b(a, contact, impulse_vec);
    apply_impulse_a(b, contact, impulse_vec);
    if friction > 0.0 {
        let jf = friction_impulse(a, b, manifold, j, friction);
        apply_impulse_a(a, contact, jf);
        apply_impulse_b(b, contact, jf);
    }
}
/// Run the sequential impulse solver over a set of contact constraints.
///
/// `bodies` is a mutable slice of all rigid body states.
/// `constraints` is a mutable slice of contact constraints (with body indices).
pub fn solve_sequential_impulses(
    bodies: &mut [RigidBodyState],
    constraints: &mut [ContactConstraint],
    config: &SequentialImpulseSolverConfig,
) {
    for c in constraints.iter_mut() {
        let (a_slice, b_slice) = bodies.split_at(c.body_b);
        let a = &a_slice[c.body_a];
        let b = &b_slice[0];
        c.prepare(a, b, config.baumgarte, config.dt);
    }
    if config.warm_start {
        for c in constraints.iter() {
            let (a_slice, b_slice) = bodies.split_at_mut(c.body_b);
            let a = &mut a_slice[c.body_a];
            let b = &mut b_slice[0];
            c.warm_start(a, b);
        }
    }
    for _ in 0..config.velocity_iterations {
        for c in constraints.iter_mut() {
            let (a_slice, b_slice) = bodies.split_at_mut(c.body_b);
            let a = &mut a_slice[c.body_a];
            let b = &mut b_slice[0];
            c.solve_normal(a, b);
            c.solve_friction_t1(a, b);
            c.solve_friction_t2(a, b);
        }
    }
}
/// Apply a restitution impulse pair with a velocity threshold guard.
///
/// If the closing speed `|v_rel·n|` is below `threshold`, the restitution
/// coefficient is set to zero (perfectly inelastic) to prevent micro-bouncing.
///
/// Returns the actual restitution used.
pub fn restitution_with_threshold(
    a: &mut RigidBodyState,
    b: &mut RigidBodyState,
    manifold: &CollisionManifold,
    restitution: f64,
    threshold: f64,
) -> f64 {
    let v_rel = relative_velocity_at_contact(a, b, manifold.contact_point);
    let v_rel_n = dot3(v_rel, manifold.normal);
    let effective_e = if v_rel_n.abs() < threshold {
        0.0
    } else {
        restitution
    };
    let j = {
        if v_rel_n <= 0.0 {
            0.0
        } else {
            let denom =
                effective_mass_along_direction(a, b, manifold.contact_point, manifold.normal);
            if denom.abs() < 1e-30 {
                0.0
            } else {
                (1.0 + effective_e) * v_rel_n / denom
            }
        }
    };
    if j > 1e-30 {
        let impulse_vec = scale3(manifold.normal, j);
        apply_impulse_b(a, manifold.contact_point, impulse_vec);
        apply_impulse_a(b, manifold.contact_point, impulse_vec);
    }
    effective_e
}
#[cfg(test)]
mod tests {
    use super::*;
    fn sphere_state(
        pos: [f64; 3],
        vel: [f64; 3],
        inv_mass: f64,
        inv_inertia: f64,
    ) -> RigidBodyState {
        RigidBodyState {
            position: pos,
            velocity: vel,
            angular_velocity: [0.0, 0.0, 0.0],
            orientation: [0.0, 0.0, 0.0, 1.0],
            inv_mass,
            inv_inertia_local: [inv_inertia, inv_inertia, inv_inertia],
        }
    }
    /// Two equal-mass spheres colliding head-on with e=1 (elastic): velocities swap.
    #[test]
    fn test_head_on_elastic_equal_mass() {
        let mut a = sphere_state([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere_state([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        resolve_collision(&mut a, &mut b, &manifold, 1.0, 0.0);
        assert!(
            (a.velocity[0] - (-1.0)).abs() < 1e-10,
            "a.vx = {}",
            a.velocity[0]
        );
        assert!(
            (b.velocity[0] - 1.0).abs() < 1e-10,
            "b.vx = {}",
            b.velocity[0]
        );
    }
    /// Ball hits a static wall (inv_mass=0): normal velocity component flips.
    #[test]
    fn test_static_wall_reflection() {
        let mut a = sphere_state([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.0);
        let mut wall = sphere_state([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.01,
            contact_point: [0.5, 0.0, 0.0],
        };
        resolve_collision(&mut a, &mut wall, &manifold, 1.0, 0.0);
        assert!(
            (a.velocity[0] - (-2.0)).abs() < 1e-10,
            "a.vx = {}",
            a.velocity[0]
        );
        assert!(wall.velocity[0].abs() < 1e-30);
    }
    /// compute_impulse: two same-mass bodies approaching at 1 m/s with e=0.
    #[test]
    fn test_perfectly_inelastic_impulse() {
        let a = sphere_state([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0);
        let b = sphere_state([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let j = compute_impulse(&a, &b, &manifold, 0.0);
        assert!((j - 0.5).abs() < 1e-10, "j = {}", j);
    }
    /// friction_impulse: zero tangential velocity → zero friction impulse.
    #[test]
    fn test_zero_tangential_friction() {
        let a = sphere_state([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0);
        let b = sphere_state([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let jf = friction_impulse(&a, &b, &manifold, 2.0, 0.5);
        assert!(len3(jf) < 1e-12, "friction should be zero but got {:?}", jf);
    }
    /// apply_impulse_a: linear velocity changes by j*n/m.
    #[test]
    fn test_apply_impulse_a_linear_change() {
        let mut state = sphere_state([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 2.0, 0.0);
        let contact = [0.0, 0.0, 0.0];
        let impulse_vec = [3.0, 0.0, 0.0];
        apply_impulse_a(&mut state, contact, impulse_vec);
        assert!(
            (state.velocity[0] - 6.0).abs() < 1e-10,
            "vx = {}",
            state.velocity[0]
        );
        let mut state2 = sphere_state([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        apply_impulse_a(&mut state2, contact, impulse_vec);
        assert!(
            (state2.velocity[0] - 3.0).abs() < 1e-10,
            "vx = {}",
            state2.velocity[0]
        );
    }
    #[test]
    fn test_normalize3() {
        let v = normalize3([3.0, 4.0, 0.0]);
        assert!((len3(v) - 1.0).abs() < 1e-10);
        assert!((v[0] - 0.6).abs() < 1e-10);
        assert!((v[1] - 0.8).abs() < 1e-10);
        let z = normalize3([0.0, 0.0, 0.0]);
        assert!(len3(z) < 1e-20);
    }
    #[test]
    fn test_cross3() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = cross3(x, y);
        assert!((z[2] - 1.0).abs() < 1e-15);
        assert!(z[0].abs() < 1e-15);
        assert!(z[1].abs() < 1e-15);
    }
    #[test]
    fn test_restitution_model_newton() {
        let model = RestitutionModel::Newton { coefficient: 0.5 };
        assert!((model.effective_restitution(10.0) - 0.5).abs() < 1e-15);
        assert!((model.effective_restitution(0.001) - 0.5).abs() < 1e-15);
    }
    #[test]
    fn test_restitution_model_velocity_dependent() {
        let model = RestitutionModel::VelocityDependent {
            coefficient: 0.8,
            threshold: 0.5,
        };
        assert!((model.effective_restitution(1.0) - 0.8).abs() < 1e-15);
        assert!((model.effective_restitution(0.1) - 0.0).abs() < 1e-15);
    }
    #[test]
    fn test_restitution_model_poisson() {
        let model = RestitutionModel::Poisson { coefficient: 0.6 };
        assert!((model.effective_restitution(5.0) - 0.6).abs() < 1e-15);
    }
    #[test]
    fn test_compute_impulse_with_model_newton() {
        let a = sphere_state([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0);
        let b = sphere_state([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let model = RestitutionModel::Newton { coefficient: 1.0 };
        let j = compute_impulse_with_model(&a, &b, &manifold, model);
        let j2 = compute_impulse(&a, &b, &manifold, 1.0);
        assert!((j - j2).abs() < 1e-10, "j={}, j2={}", j, j2);
    }
    #[test]
    fn test_compute_impulse_with_model_velocity_dependent() {
        let a = sphere_state([0.0, 0.0, 0.0], [0.01, 0.0, 0.0], 1.0, 0.0);
        let b = sphere_state([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let model = RestitutionModel::VelocityDependent {
            coefficient: 0.9,
            threshold: 0.5,
        };
        let j = compute_impulse_with_model(&a, &b, &manifold, model);
        assert!((j - 0.005).abs() < 1e-10, "j = {}", j);
    }
    #[test]
    fn test_impulse_clamp_non_negative() {
        assert!(
            (clamp_normal_impulse(5.0, ImpulseClampStrategy::NonNegativeNormal) - 5.0).abs()
                < 1e-15
        );
        assert!(
            (clamp_normal_impulse(-3.0, ImpulseClampStrategy::NonNegativeNormal) - 0.0).abs()
                < 1e-15
        );
    }
    #[test]
    fn test_impulse_clamp_max_magnitude() {
        let strategy = ImpulseClampStrategy::MaxMagnitude { max_impulse: 10.0 };
        assert!((clamp_normal_impulse(5.0, strategy) - 5.0).abs() < 1e-15);
        assert!((clamp_normal_impulse(15.0, strategy) - 10.0).abs() < 1e-15);
        assert!((clamp_normal_impulse(-15.0, strategy) - (-10.0)).abs() < 1e-15);
    }
    #[test]
    fn test_clamp_friction_to_cone() {
        let f = [3.0, 4.0, 0.0];
        let clamped = clamp_friction_to_cone(f, 5.0);
        assert!((len3(clamped) - 5.0).abs() < 1e-10);
        let clamped2 = clamp_friction_to_cone(f, 2.5);
        assert!((len3(clamped2) - 2.5).abs() < 1e-10);
        let ratio = clamped2[0] / clamped2[1];
        assert!((ratio - 0.75).abs() < 1e-10);
    }
    #[test]
    fn test_build_friction_cone_directions() {
        let normal = [0.0, 0.0, 1.0];
        let dirs = build_friction_cone_directions(normal, 4);
        assert_eq!(dirs.len(), 4);
        for d in &dirs {
            assert!((len3(*d) - 1.0).abs() < 1e-10, "not unit: {:?}", d);
        }
        for d in &dirs {
            assert!(
                dot3(*d, normal).abs() < 1e-10,
                "not perpendicular: {:?}, dot={}",
                d,
                dot3(*d, normal)
            );
        }
    }
    #[test]
    fn test_effective_mass_along_direction() {
        let a = sphere_state([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let b = sphere_state([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0);
        let contact = [1.0, 0.0, 0.0];
        let em = effective_mass_along_direction(&a, &b, contact, [1.0, 0.0, 0.0]);
        assert!((em - 2.0).abs() < 1e-10, "effective mass = {}", em);
    }
    #[test]
    fn test_friction_impulse_pyramidal() {
        let a = sphere_state([0.0, 0.0, 0.0], [1.0, 1.0, 0.0], 1.0, 0.0);
        let b = sphere_state([0.0, 0.0, 2.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [0.0, 0.0, 1.0],
            depth: 0.0,
            contact_point: [0.0, 0.0, 1.0],
        };
        let jf = friction_impulse_pyramidal(&a, &b, &manifold, 10.0, 0.5, 8);
        assert!(len3(jf) <= 5.0 + 1e-10, "friction magnitude = {}", len3(jf));
        assert!(jf[2].abs() < 1e-10, "should have no normal component");
    }
    #[test]
    fn test_resolve_collision_with_model() {
        let mut a = sphere_state([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere_state([2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let model = RestitutionModel::Newton { coefficient: 1.0 };
        resolve_collision_with_model(&mut a, &mut b, &manifold, model, 0.0);
        assert!((a.velocity[0] - (-1.0)).abs() < 1e-10);
        assert!((b.velocity[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_contact_constraint_new() {
        let c = ContactConstraint::new(0, 1, [0.0, 1.0, 0.0], [0.0, 0.0, 0.0], 0.01, 0.5, 0.3);
        assert_eq!(c.body_a, 0);
        assert_eq!(c.body_b, 1);
        assert!((c.accumulated_normal_impulse - 0.0).abs() < 1e-15);
        assert!((len3(c.tangent1) - 1.0).abs() < 1e-10);
        assert!((len3(c.tangent2) - 1.0).abs() < 1e-10);
        assert!(dot3(c.tangent1, c.normal).abs() < 1e-10);
        assert!(dot3(c.tangent2, c.normal).abs() < 1e-10);
        assert!(dot3(c.tangent1, c.tangent2).abs() < 1e-10);
    }
    #[test]
    fn test_sequential_impulse_solver_basic() {
        let mut bodies = vec![
            sphere_state([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.0),
            sphere_state([2.0, 0.0, 0.0], [-2.0, 0.0, 0.0], 1.0, 0.0),
        ];
        let mut constraints = vec![ContactConstraint::new(
            0,
            1,
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.01,
            1.0,
            0.0,
        )];
        let config = SequentialImpulseSolverConfig {
            velocity_iterations: 20,
            baumgarte: 0.2,
            dt: 1.0 / 60.0,
            warm_start: false,
        };
        solve_sequential_impulses(&mut bodies, &mut constraints, &config);
        let v_rel_n = bodies[0].velocity[0] - bodies[1].velocity[0];
        assert!(
            v_rel_n <= 0.1,
            "bodies should be separating or at rest, v_rel_n = {}",
            v_rel_n
        );
    }
    #[test]
    fn test_sequential_solver_momentum_conservation() {
        let mut bodies = vec![
            sphere_state([0.0, 0.0, 0.0], [3.0, 1.0, 0.0], 1.0, 0.0),
            sphere_state([2.0, 0.0, 0.0], [-1.0, -1.0, 0.0], 2.0, 0.0),
        ];
        let m0 = 1.0 / bodies[0].inv_mass;
        let m1 = 1.0 / bodies[1].inv_mass;
        let initial_px = m0 * bodies[0].velocity[0] + m1 * bodies[1].velocity[0];
        let initial_py = m0 * bodies[0].velocity[1] + m1 * bodies[1].velocity[1];
        let mut constraints = vec![ContactConstraint::new(
            0,
            1,
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.0,
            0.5,
            0.0,
        )];
        let config = SequentialImpulseSolverConfig {
            velocity_iterations: 20,
            baumgarte: 0.0,
            dt: 1.0 / 60.0,
            warm_start: false,
        };
        solve_sequential_impulses(&mut bodies, &mut constraints, &config);
        let final_px = m0 * bodies[0].velocity[0] + m1 * bodies[1].velocity[0];
        let final_py = m0 * bodies[0].velocity[1] + m1 * bodies[1].velocity[1];
        assert!(
            (final_px - initial_px).abs() < 1e-6,
            "px: initial={}, final={}",
            initial_px,
            final_px
        );
        assert!(
            (final_py - initial_py).abs() < 1e-6,
            "py: initial={}, final={}",
            initial_py,
            final_py
        );
    }
    #[test]
    fn test_contact_constraint_warm_start() {
        let mut a = sphere_state([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let mut b = sphere_state([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let mut c = ContactConstraint::new(0, 1, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0, 0.5, 0.3);
        c.accumulated_normal_impulse = 5.0;
        c.warm_start(&mut a, &mut b);
        assert!(a.velocity[0] < 0.0, "a should move in -x");
        assert!(b.velocity[0] > 0.0, "b should move in +x");
    }
    #[test]
    fn test_relative_velocity_with_angular() {
        let mut a = sphere_state([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 1.0);
        a.angular_velocity = [0.0, 0.0, 1.0];
        let b = sphere_state([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0);
        let contact = [1.0, 0.0, 0.0];
        let v_rel = relative_velocity_at_contact(&a, &b, contact);
        assert!((v_rel[0]).abs() < 1e-10);
        assert!((v_rel[1] - 1.0).abs() < 1e-10);
        assert!((v_rel[2]).abs() < 1e-10);
    }
    #[test]
    fn test_separating_bodies_no_impulse() {
        let a = sphere_state([0.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0);
        let b = sphere_state([2.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0);
        let manifold = CollisionManifold {
            normal: [1.0, 0.0, 0.0],
            depth: 0.0,
            contact_point: [1.0, 0.0, 0.0],
        };
        let j = compute_impulse(&a, &b, &manifold, 1.0);
        assert!(
            j.abs() < 1e-15,
            "separating bodies should have j=0, got {}",
            j
        );
    }
    #[test]
    fn test_friction_with_tangential_velocity() {
        let a = sphere_state([0.0, 0.0, 0.0], [1.0, 2.0, 0.0], 1.0, 0.0);
        let b = sphere_state([0.0, 0.0, 2.0], [0.0, 0.0, 0.0], 0.0, 0.0);
        let manifold = CollisionManifold {
            normal: [0.0, 0.0, 1.0],
            depth: 0.0,
            contact_point: [0.0, 0.0, 1.0],
        };
        let jf = friction_impulse(&a, &b, &manifold, 10.0, 0.5);
        assert!(len3(jf) <= 5.0 + 1e-10);
        assert!(jf[2].abs() < 1e-10, "should have no normal component");
    }
}
