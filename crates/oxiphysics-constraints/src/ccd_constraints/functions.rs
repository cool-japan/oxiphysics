//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::{Real, Vec3};

use super::types::CcdBroadphasePair;

/// Estimate the time of impact for two linearly-swept point-masses.
///
/// Returns `Some(toi)` in \[0, 1\] if the swept segments come within
/// `sum_radius` of each other, or `None` if they do not collide.
///
/// This is a simplified swept-sphere TOI suitable for CCD broadphase.
pub fn estimate_linear_toi(
    pos_a: Vec3,
    vel_a: Vec3,
    pos_b: Vec3,
    vel_b: Vec3,
    sum_radius: f64,
    dt: f64,
) -> Option<f64> {
    let rel_pos = pos_a - pos_b;
    let rel_vel = (vel_a - vel_b) * dt;
    let a = rel_vel.dot(&rel_vel);
    let b = 2.0 * rel_pos.dot(&rel_vel);
    let c = rel_pos.dot(&rel_pos) - sum_radius * sum_radius;
    if c < 0.0 {
        return Some(0.0);
    }
    if a < 1e-12 {
        return None;
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return None;
    }
    let sqrt_disc = discriminant.sqrt();
    let t = (-b - sqrt_disc) / (2.0 * a);
    if (0.0..=1.0).contains(&t) {
        Some(t)
    } else {
        None
    }
}
/// Conservative advancement for sphere-sphere CCD.
///
/// Iteratively advances two spheres toward each other by a safe step
/// that guarantees no tunneling, until either:
/// - The spheres are within `contact_tolerance` of touching
/// - The advancement reaches `t = 1` (no collision)
/// - The maximum number of iterations is reached
///
/// Returns `Some(toi)` if contact is detected, `None` otherwise.
pub fn conservative_advancement_spheres(
    pos_a: Vec3,
    vel_a: Vec3,
    radius_a: f64,
    pos_b: Vec3,
    vel_b: Vec3,
    radius_b: f64,
    dt: f64,
    contact_tolerance: f64,
    max_iterations: usize,
) -> Option<f64> {
    let sum_radius = radius_a + radius_b;
    let mut t = 0.0;
    let mut p_a = pos_a;
    let mut p_b = pos_b;
    for _ in 0..max_iterations {
        let diff = p_a - p_b;
        let dist = diff.norm();
        let gap = dist - sum_radius;
        if gap <= contact_tolerance {
            return Some(t);
        }
        let rel_vel = vel_a - vel_b;
        let closing_speed = if dist > 1e-12 {
            -(rel_vel.dot(&diff) / dist)
        } else {
            rel_vel.norm()
        };
        if closing_speed <= 1e-12 {
            return None;
        }
        let delta_t = gap / closing_speed;
        let delta_t = delta_t.min((1.0 - t) * dt);
        if delta_t < 1e-14 * dt {
            return Some(t);
        }
        t += delta_t / dt;
        if t >= 1.0 {
            return None;
        }
        p_a += vel_a * delta_t;
        p_b += vel_b * delta_t;
    }
    None
}
/// Refine a coarse TOI estimate using bisection.
///
/// Given a coarse TOI bracket \[t_lo, t_hi\] where:
/// - At t_lo, distance > sum_radius (no contact)
/// - At t_hi, distance ≤ sum_radius (contact/penetration)
///
/// Returns a refined TOI within `tolerance`.
pub fn refine_toi_bisection(
    pos_a: Vec3,
    vel_a: Vec3,
    pos_b: Vec3,
    vel_b: Vec3,
    sum_radius: f64,
    dt: f64,
    t_lo: f64,
    t_hi: f64,
    tolerance: f64,
    max_iterations: usize,
) -> f64 {
    let mut lo = t_lo;
    let mut hi = t_hi;
    for _ in 0..max_iterations {
        let mid = (lo + hi) * 0.5;
        let t_sec = mid * dt;
        let pa = pos_a + vel_a * t_sec;
        let pb = pos_b + vel_b * t_sec;
        let dist = (pa - pb).norm();
        if dist > sum_radius {
            lo = mid;
        } else {
            hi = mid;
        }
        if (hi - lo) < tolerance {
            break;
        }
    }
    (lo + hi) * 0.5
}
/// Refine TOI using Newton-Raphson on the distance function.
///
/// f(t) = |p_a(t) - p_b(t)| - R
/// f'(t) = d/dt |p_a - p_b| = (p_rel · v_rel) / |p_rel| * dt
pub fn refine_toi_newton(
    pos_a: Vec3,
    vel_a: Vec3,
    pos_b: Vec3,
    vel_b: Vec3,
    sum_radius: f64,
    dt: f64,
    t_initial: f64,
    tolerance: f64,
    max_iterations: usize,
) -> f64 {
    let rel_vel = (vel_a - vel_b) * dt;
    let mut t = t_initial;
    for _ in 0..max_iterations {
        let p_rel = (pos_a + vel_a * (t * dt)) - (pos_b + vel_b * (t * dt));
        let dist = p_rel.norm();
        let f_val = dist - sum_radius;
        if f_val.abs() < tolerance {
            break;
        }
        if dist < 1e-12 {
            break;
        }
        let f_prime = p_rel.dot(&rel_vel) / dist;
        if f_prime.abs() < 1e-14 {
            break;
        }
        let delta = f_val / f_prime;
        t -= delta;
        t = t.clamp(0.0, 1.0);
    }
    t
}
/// Swept AABB (Axis-Aligned Bounding Box) for CCD broadphase.
///
/// Computes the bounding box that encloses the entire swept path of
/// a sphere-like object over a time step.
pub fn swept_aabb(pos: Vec3, vel: Vec3, radius: f64, dt: f64) -> (Vec3, Vec3) {
    let end_pos = pos + vel * dt;
    let min_x = pos.x.min(end_pos.x) - radius;
    let min_y = pos.y.min(end_pos.y) - radius;
    let min_z = pos.z.min(end_pos.z) - radius;
    let max_x = pos.x.max(end_pos.x) + radius;
    let max_y = pos.y.max(end_pos.y) + radius;
    let max_z = pos.z.max(end_pos.z) + radius;
    (
        Vec3::new(min_x, min_y, min_z),
        Vec3::new(max_x, max_y, max_z),
    )
}
/// Check if two swept AABBs overlap.
pub fn swept_aabb_overlap(min_a: Vec3, max_a: Vec3, min_b: Vec3, max_b: Vec3) -> bool {
    min_a.x <= max_b.x
        && max_a.x >= min_b.x
        && min_a.y <= max_b.y
        && max_a.y >= min_b.y
        && min_a.z <= max_b.z
        && max_a.z >= min_b.z
}
/// Swept sphere vs. plane intersection test.
///
/// Returns `Some(toi)` in \[0, dt\] if the sphere intersects the plane.
/// The plane is defined by a point and normal.
pub fn swept_sphere_plane(
    sphere_pos: Vec3,
    sphere_vel: Vec3,
    sphere_radius: f64,
    plane_point: Vec3,
    plane_normal: Vec3,
    dt: f64,
) -> Option<f64> {
    let n = plane_normal.normalize();
    let d0 = (sphere_pos - plane_point).dot(&n) - sphere_radius;
    let d1 = ((sphere_pos + sphere_vel * dt) - plane_point).dot(&n) - sphere_radius;
    if d0 >= 0.0 && d1 >= 0.0 {
        return None;
    }
    if d0 < 0.0 {
        return Some(0.0);
    }
    let denom = d0 - d1;
    if denom.abs() < 1e-12 {
        return None;
    }
    let t = d0 / denom;
    if (0.0..=1.0).contains(&t) {
        Some(t * dt)
    } else {
        None
    }
}
/// Estimate maximum displacement of a point on a rotating body.
///
/// For a body with angular velocity ω and a point at distance r from
/// the center of mass, the maximum linear speed due to rotation is |ω| * r.
///
/// Returns the maximum displacement over `dt`.
pub fn max_rotation_displacement(angular_velocity: Vec3, max_radius: f64, dt: f64) -> f64 {
    angular_velocity.norm() * max_radius * dt
}
/// Check if CCD is needed for a body based on its velocity and size.
///
/// CCD is needed when the body's maximum displacement per frame exceeds
/// a fraction of its size (CCD activation threshold).
pub fn ccd_needed(
    velocity: Vec3,
    angular_velocity: Vec3,
    max_radius: f64,
    dt: f64,
    activation_threshold: f64,
) -> bool {
    let linear_disp = velocity.norm() * dt;
    let angular_disp = max_rotation_displacement(angular_velocity, max_radius, dt);
    let total_disp = linear_disp + angular_disp;
    total_disp > activation_threshold * max_radius
}
/// Conservative upper bound on the swept distance including rotation.
///
/// For a rigid body with linear velocity v and angular velocity ω,
/// the maximum displacement of any point on the body over dt is:
///
/// d_max = |v| * dt + |ω| * r_max * dt
pub fn conservative_swept_distance(
    velocity: Vec3,
    angular_velocity: Vec3,
    max_radius: f64,
    dt: f64,
) -> f64 {
    velocity.norm() * dt + angular_velocity.norm() * max_radius * dt
}
/// Bilateral advancement: both bodies advance simultaneously.
///
/// Unlike simple conservative advancement (which advances one body at a time),
/// bilateral advancement considers the relative motion of both bodies
/// to compute a tighter safe step.
///
/// Returns `Some(toi)` if contact is detected, `None` otherwise.
pub fn bilateral_advancement_spheres(
    pos_a: Vec3,
    vel_a: Vec3,
    radius_a: f64,
    pos_b: Vec3,
    vel_b: Vec3,
    radius_b: f64,
    dt: f64,
    contact_tolerance: f64,
    max_iterations: usize,
) -> Option<f64> {
    let sum_radius = radius_a + radius_b;
    let rel_pos_init = pos_a - pos_b;
    let rel_vel = vel_a - vel_b;
    if rel_pos_init.norm() <= sum_radius + contact_tolerance {
        return Some(0.0);
    }
    let rel_speed = rel_vel.norm();
    if rel_speed < 1e-12 {
        return None;
    }
    let mut t = 0.0;
    for _ in 0..max_iterations {
        let current_t = t * dt;
        let p_rel = rel_pos_init + rel_vel * current_t;
        let dist = p_rel.norm();
        let gap = dist - sum_radius;
        if gap <= contact_tolerance {
            return Some(t);
        }
        let closing_speed = if dist > 1e-12 {
            -(p_rel.dot(&rel_vel)) / dist
        } else {
            rel_speed
        };
        if closing_speed <= 1e-12 {
            return None;
        }
        let delta_t = (gap / closing_speed) / dt;
        let delta_t = delta_t.max(1e-14);
        t += delta_t;
        if t >= 1.0 {
            return None;
        }
    }
    None
}
/// Parameters for [`bilateral_advancement_with_rotation`].
#[derive(Debug, Clone, Copy)]
pub struct RotationalAdvancementParams {
    /// Position of body A.
    pub pos_a: Vec3,
    /// Linear velocity of body A.
    pub vel_a: Vec3,
    /// Angular velocity of body A.
    pub ang_vel_a: Vec3,
    /// Bounding radius of body A.
    pub radius_a: f64,
    /// Position of body B.
    pub pos_b: Vec3,
    /// Linear velocity of body B.
    pub vel_b: Vec3,
    /// Angular velocity of body B.
    pub ang_vel_b: Vec3,
    /// Bounding radius of body B.
    pub radius_b: f64,
    /// Time step length.
    pub dt: f64,
    /// Contact tolerance (gap threshold for a hit).
    pub contact_tolerance: f64,
    /// Maximum number of advancement iterations.
    pub max_iterations: usize,
}

/// Bilateral advancement with angular velocity consideration.
///
/// Accounts for the fact that rotating bodies can have points moving
/// faster than their center of mass, requiring smaller safe steps.
pub fn bilateral_advancement_with_rotation(p: RotationalAdvancementParams) -> Option<f64> {
    let RotationalAdvancementParams {
        pos_a,
        vel_a,
        ang_vel_a,
        radius_a,
        pos_b,
        vel_b,
        ang_vel_b,
        radius_b,
        dt,
        contact_tolerance,
        max_iterations,
    } = p;
    let sum_radius = radius_a + radius_b;
    let max_speed_a = vel_a.norm() + ang_vel_a.norm() * radius_a;
    let max_speed_b = vel_b.norm() + ang_vel_b.norm() * radius_b;
    let max_rel_speed = max_speed_a + max_speed_b;
    if max_rel_speed < 1e-12 {
        let dist = (pos_a - pos_b).norm();
        return if dist <= sum_radius + contact_tolerance {
            Some(0.0)
        } else {
            None
        };
    }
    let mut t = 0.0;
    for _ in 0..max_iterations {
        let current_t = t * dt;
        let pa = pos_a + vel_a * current_t;
        let pb = pos_b + vel_b * current_t;
        let diff = pa - pb;
        let dist = diff.norm();
        let gap = dist - sum_radius;
        if gap <= contact_tolerance {
            return Some(t);
        }
        let delta_t = (gap / max_rel_speed) / dt;
        let delta_t = delta_t.max(1e-14);
        t += delta_t;
        if t >= 1.0 {
            return None;
        }
    }
    None
}
/// Computes the restitution impulse at the exact TOI for a sphere-sphere
/// collision (or any two point-mass bodies).
///
/// # Algorithm
///
/// 1. Compute relative normal velocity `v_rel_n` at TOI.
/// 2. Apply Newton's law of restitution: `Δv = -(1 + e) * v_rel_n`.
/// 3. Compute impulse `j = Δv / (1/m_a + 1/m_b)`.
///
/// Returns `(j, v_post_a, v_post_b)`.
pub fn toi_restitution_impulse(
    vel_a: Vec3,
    inv_mass_a: Real,
    vel_b: Vec3,
    inv_mass_b: Real,
    normal: Vec3,
    restitution: Real,
) -> (Real, Vec3, Vec3) {
    let rel_vel_n = (vel_a - vel_b).dot(&normal);
    if rel_vel_n >= 0.0 {
        return (0.0, vel_a, vel_b);
    }
    let inv_mass_sum = inv_mass_a + inv_mass_b;
    if inv_mass_sum < 1e-15 {
        return (0.0, vel_a, vel_b);
    }
    let j = -(1.0 + restitution) * rel_vel_n / inv_mass_sum;
    let impulse = normal * j;
    let v_post_a = vel_a + impulse * inv_mass_a;
    let v_post_b = vel_b - impulse * inv_mass_b;
    (j, v_post_a, v_post_b)
}
/// Compute the kinetic energy change due to a restitution impulse.
///
/// Returns `(KE_before, KE_after)`.
pub fn toi_kinetic_energy_change(
    vel_a: Vec3,
    mass_a: Real,
    vel_b: Vec3,
    mass_b: Real,
    vel_post_a: Vec3,
    vel_post_b: Vec3,
) -> (Real, Real) {
    let ke_before = 0.5 * mass_a * vel_a.dot(&vel_a) + 0.5 * mass_b * vel_b.dot(&vel_b);
    let ke_after =
        0.5 * mass_a * vel_post_a.dot(&vel_post_a) + 0.5 * mass_b * vel_post_b.dot(&vel_post_b);
    (ke_before, ke_after)
}
/// Filter CCD broadphase pairs by a TOI threshold.
///
/// Returns only pairs whose upper-bound TOI is ≤ `toi_threshold`.
pub fn filter_ccd_pairs(
    pairs: &[CcdBroadphasePair],
    toi_threshold: Real,
) -> Vec<&CcdBroadphasePair> {
    pairs
        .iter()
        .filter(|p| p.toi_upper_bound <= toi_threshold)
        .collect()
}
/// Sort broadphase pairs by ascending TOI upper bound.
pub fn sort_ccd_pairs_by_toi(pairs: &mut [CcdBroadphasePair]) {
    pairs.sort_by(|a, b| {
        a.toi_upper_bound
            .partial_cmp(&b.toi_upper_bound)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::CcdConstraint;
    use crate::CcdConstraintSolver;
    use oxiphysics_core::BodyHandle;
    use oxiphysics_rigid::RigidBody;
    use oxiphysics_rigid::RigidBodySet;
    fn make_bodies() -> (RigidBodySet, BodyHandle, BodyHandle) {
        let mut bodies = RigidBodySet::new();
        let mut body_a = RigidBody::new(1.0);
        body_a.transform.position = Vec3::new(0.0, 2.0, 0.0);
        body_a.velocity = Vec3::new(0.0, -100.0, 0.0);
        body_a.linear_damping = 0.0;
        body_a.angular_damping = 0.0;
        body_a.gravity_scale = 0.0;
        let ha = bodies.insert(body_a);
        let mut body_b = RigidBody::new_static();
        body_b.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let hb = bodies.insert(body_b);
        (bodies, ha, hb)
    }
    #[test]
    fn test_ccd_constraint_creation() {
        let (bodies, ha, hb) = make_bodies();
        let _ = bodies;
        let c = CcdConstraint::new(
            ha,
            hb,
            0.02,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.3,
        );
        assert!((c.toi - 0.02).abs() < 1e-12);
        assert_eq!(c.body_a, ha);
        assert_eq!(c.body_b, hb);
    }
    #[test]
    fn test_ccd_prevents_tunneling() {
        let (mut bodies, ha, hb) = make_bodies();
        let dt = 1.0 / 60.0;
        let c = CcdConstraint::new(
            ha,
            hb,
            0.02,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
        );
        let solver = CcdConstraintSolver::new(4, Vec3::zeros());
        solver.solve(&[c], &mut bodies, dt);
        let body_a = bodies.get(ha).unwrap();
        assert!(
            body_a.velocity.y >= -1.0,
            "CCD should have stopped/reversed the downward velocity, vy = {}",
            body_a.velocity.y
        );
    }
    #[test]
    fn test_ccd_multiple_constraints_toi_order() {
        let mut bodies = RigidBodySet::new();
        let mut b1 = RigidBody::new(1.0);
        b1.transform.position = Vec3::new(0.0, 5.0, 0.0);
        b1.velocity = Vec3::new(0.0, -50.0, 0.0);
        b1.linear_damping = 0.0;
        b1.angular_damping = 0.0;
        b1.gravity_scale = 0.0;
        let h1 = bodies.insert(b1);
        let mut b2 = RigidBody::new(1.0);
        b2.transform.position = Vec3::new(1.0, 3.0, 0.0);
        b2.velocity = Vec3::new(-50.0, 0.0, 0.0);
        b2.linear_damping = 0.0;
        b2.angular_damping = 0.0;
        b2.gravity_scale = 0.0;
        let h2 = bodies.insert(b2);
        let ground = RigidBody::new_static();
        let hg = bodies.insert(ground);
        let c1 = CcdConstraint::new(
            h1,
            hg,
            0.1,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
        );
        let c2 = CcdConstraint::new(
            h2,
            hg,
            0.05,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
            0.0,
        );
        let solver = CcdConstraintSolver::new(4, Vec3::zeros());
        solver.solve(&[c1, c2], &mut bodies, 1.0 / 60.0);
        assert!(bodies.get(h1).is_some());
        assert!(bodies.get(h2).is_some());
    }
    #[test]
    fn test_ccd_no_impulse_when_separating() {
        let mut bodies = RigidBodySet::new();
        let mut ba = RigidBody::new(1.0);
        ba.transform.position = Vec3::new(0.0, 1.0, 0.0);
        ba.velocity = Vec3::new(0.0, 10.0, 0.0);
        ba.linear_damping = 0.0;
        ba.angular_damping = 0.0;
        ba.gravity_scale = 0.0;
        let ha = bodies.insert(ba);
        let hb = bodies.insert(RigidBody::new_static());
        let c = CcdConstraint::new(
            ha,
            hb,
            0.5,
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            0.0,
            0.0,
        );
        let solver = CcdConstraintSolver::new(4, Vec3::zeros());
        solver.solve(&[c], &mut bodies, 1.0 / 60.0);
        let final_vy = bodies.get(ha).unwrap().velocity.y;
        assert!(
            final_vy > 0.0,
            "Body moving away should not be reversed: final_vy={final_vy}"
        );
    }
    #[test]
    fn test_linear_toi_estimate_collision() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, -10.0, 0.0);
        let pos_b = Vec3::new(0.0, 0.0, 0.0);
        let vel_b = Vec3::zeros();
        let toi = estimate_linear_toi(pos_a, vel_a, pos_b, vel_b, 0.5, 1.0);
        assert!(toi.is_some(), "Should detect collision");
        let t = toi.unwrap();
        assert!((0.0..=1.0).contains(&t), "TOI must be in [0,1]: {t}");
    }
    #[test]
    fn test_linear_toi_estimate_no_collision() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, 10.0, 0.0);
        let pos_b = Vec3::new(0.0, 0.0, 0.0);
        let vel_b = Vec3::zeros();
        let toi = estimate_linear_toi(pos_a, vel_a, pos_b, vel_b, 0.5, 1.0);
        assert!(toi.is_none(), "Moving-away bodies should not collide");
    }
    #[test]
    fn test_linear_toi_already_overlapping() {
        let pos_a = Vec3::new(0.0, 0.2, 0.0);
        let vel_a = Vec3::zeros();
        let pos_b = Vec3::zeros();
        let vel_b = Vec3::zeros();
        let toi = estimate_linear_toi(pos_a, vel_a, pos_b, vel_b, 0.5, 1.0);
        assert_eq!(toi, Some(0.0), "Overlapping bodies → TOI = 0");
    }
    #[test]
    fn test_conservative_advancement_collision() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, -10.0, 0.0);
        let pos_b = Vec3::new(0.0, 0.0, 0.0);
        let vel_b = Vec3::zeros();
        let toi =
            conservative_advancement_spheres(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 1.0, 1e-6, 100);
        assert!(toi.is_some(), "Should detect collision");
        let t = toi.unwrap();
        assert!(t > 0.0 && t < 1.0, "TOI should be in (0,1): {t}");
    }
    #[test]
    fn test_conservative_advancement_no_collision() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, 10.0, 0.0);
        let pos_b = Vec3::new(0.0, 0.0, 0.0);
        let vel_b = Vec3::zeros();
        let toi =
            conservative_advancement_spheres(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 1.0, 1e-6, 100);
        assert!(toi.is_none(), "Should not detect collision");
    }
    #[test]
    fn test_conservative_advancement_already_touching() {
        let pos_a = Vec3::new(0.0, 1.0, 0.0);
        let vel_a = Vec3::new(0.0, -1.0, 0.0);
        let pos_b = Vec3::zeros();
        let vel_b = Vec3::zeros();
        let toi =
            conservative_advancement_spheres(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 1.0, 0.01, 100);
        assert_eq!(toi, Some(0.0), "Already touching → TOI = 0");
    }
    #[test]
    fn test_toi_refinement_bisection() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, -10.0, 0.0);
        let pos_b = Vec3::zeros();
        let vel_b = Vec3::zeros();
        let sum_r = 1.0;
        let t_refined =
            refine_toi_bisection(pos_a, vel_a, pos_b, vel_b, sum_r, 1.0, 0.0, 1.0, 1e-8, 100);
        let pa = pos_a + vel_a * t_refined;
        let dist = (pa - pos_b).norm();
        assert!(
            (dist - sum_r).abs() < 0.01,
            "At refined TOI, dist should ≈ sum_r: dist={dist}, sum_r={sum_r}"
        );
    }
    #[test]
    fn test_toi_refinement_newton() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, -10.0, 0.0);
        let pos_b = Vec3::zeros();
        let vel_b = Vec3::zeros();
        let sum_r = 1.0;
        let t_refined = refine_toi_newton(pos_a, vel_a, pos_b, vel_b, sum_r, 1.0, 0.4, 1e-8, 20);
        let pa = pos_a + vel_a * (t_refined * 1.0);
        let dist = (pa - pos_b).norm();
        assert!(
            (dist - sum_r).abs() < 0.01,
            "Newton refinement: dist={dist}, expected ≈{sum_r}"
        );
    }
    #[test]
    fn test_swept_aabb() {
        let pos = Vec3::new(0.0, 5.0, 0.0);
        let vel = Vec3::new(10.0, -20.0, 0.0);
        let (min, max) = swept_aabb(pos, vel, 0.5, 1.0);
        assert!(min.x <= -0.5);
        assert!(max.x >= 10.5);
        assert!(min.y <= -15.5);
        assert!(max.y >= 5.5);
    }
    #[test]
    fn test_swept_aabb_overlap_test() {
        let (min_a, max_a) = (Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 2.0, 2.0));
        let (min_b, max_b) = (Vec3::new(1.0, 1.0, 1.0), Vec3::new(3.0, 3.0, 3.0));
        assert!(swept_aabb_overlap(min_a, max_a, min_b, max_b));
        let (min_c, max_c) = (Vec3::new(5.0, 5.0, 5.0), Vec3::new(6.0, 6.0, 6.0));
        assert!(!swept_aabb_overlap(min_a, max_a, min_c, max_c));
    }
    #[test]
    fn test_swept_sphere_plane_collision() {
        let pos = Vec3::new(0.0, 5.0, 0.0);
        let vel = Vec3::new(0.0, -10.0, 0.0);
        let plane_pt = Vec3::zeros();
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        let toi = swept_sphere_plane(pos, vel, 0.5, plane_pt, plane_n, 1.0);
        assert!(toi.is_some(), "Should detect sphere-plane collision");
        let t = toi.unwrap();
        assert!(t > 0.0 && t < 1.0, "TOI should be in (0,1): {t}");
    }
    #[test]
    fn test_swept_sphere_plane_no_collision() {
        let pos = Vec3::new(0.0, 5.0, 0.0);
        let vel = Vec3::new(0.0, 10.0, 0.0);
        let plane_pt = Vec3::zeros();
        let plane_n = Vec3::new(0.0, 1.0, 0.0);
        let toi = swept_sphere_plane(pos, vel, 0.5, plane_pt, plane_n, 1.0);
        assert!(toi.is_none());
    }
    #[test]
    fn test_max_rotation_displacement() {
        let omega = Vec3::new(0.0, 10.0, 0.0);
        let disp = max_rotation_displacement(omega, 2.0, 0.01);
        assert!((disp - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_ccd_needed_slow_body() {
        let vel = Vec3::new(0.1, 0.0, 0.0);
        let omega = Vec3::zeros();
        assert!(!ccd_needed(vel, omega, 1.0, 1.0 / 60.0, 0.5));
    }
    #[test]
    fn test_ccd_needed_fast_body() {
        let vel = Vec3::new(100.0, 0.0, 0.0);
        let omega = Vec3::zeros();
        assert!(ccd_needed(vel, omega, 1.0, 1.0 / 60.0, 0.5));
    }
    #[test]
    fn test_ccd_needed_fast_rotation() {
        let vel = Vec3::zeros();
        let omega = Vec3::new(0.0, 100.0, 0.0);
        assert!(ccd_needed(vel, omega, 2.0, 1.0 / 60.0, 0.5));
    }
    #[test]
    fn test_conservative_swept_distance() {
        let vel = Vec3::new(10.0, 0.0, 0.0);
        let omega = Vec3::new(0.0, 5.0, 0.0);
        let d = conservative_swept_distance(vel, omega, 2.0, 0.01);
        assert!((d - 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_bilateral_advancement_collision() {
        let pos_a = Vec3::new(-5.0, 0.0, 0.0);
        let vel_a = Vec3::new(10.0, 0.0, 0.0);
        let pos_b = Vec3::new(5.0, 0.0, 0.0);
        let vel_b = Vec3::new(-10.0, 0.0, 0.0);
        let toi =
            bilateral_advancement_spheres(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 1.0, 1e-4, 100);
        assert!(toi.is_some(), "Should detect head-on collision");
    }
    #[test]
    fn test_bilateral_advancement_no_collision() {
        let pos_a = Vec3::new(-5.0, 0.0, 0.0);
        let vel_a = Vec3::new(-10.0, 0.0, 0.0);
        let pos_b = Vec3::new(5.0, 0.0, 0.0);
        let vel_b = Vec3::new(10.0, 0.0, 0.0);
        let toi =
            bilateral_advancement_spheres(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 1.0, 1e-4, 100);
        assert!(toi.is_none());
    }
    #[test]
    fn test_bilateral_advancement_with_rotation() {
        let pos_a = Vec3::new(0.0, 5.0, 0.0);
        let vel_a = Vec3::new(0.0, -10.0, 0.0);
        let ang_a = Vec3::new(0.0, 0.0, 5.0);
        let pos_b = Vec3::zeros();
        let vel_b = Vec3::zeros();
        let ang_b = Vec3::zeros();
        let toi = bilateral_advancement_with_rotation(RotationalAdvancementParams {
            pos_a,
            vel_a,
            ang_vel_a: ang_a,
            radius_a: 0.5,
            pos_b,
            vel_b,
            ang_vel_b: ang_b,
            radius_b: 0.5,
            dt: 1.0,
            contact_tolerance: 1e-4,
            max_iterations: 200,
        });
        assert!(toi.is_some(), "Should detect collision with rotation");
    }
    #[test]
    fn test_bilateral_already_touching() {
        let pos_a = Vec3::new(0.0, 0.9, 0.0);
        let vel_a = Vec3::zeros();
        let pos_b = Vec3::zeros();
        let vel_b = Vec3::zeros();
        let toi =
            bilateral_advancement_spheres(pos_a, vel_a, 0.5, pos_b, vel_b, 0.5, 1.0, 0.1, 100);
        assert_eq!(toi, Some(0.0));
    }
}
