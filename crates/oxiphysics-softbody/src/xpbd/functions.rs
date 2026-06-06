//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use oxiphysics_core::math::Vec3;

/// Compute the effective compliance (alpha_tilde) for a given compliance and
/// time step: alpha_tilde = alpha / dt^2.
///
/// Smaller dt makes the effective compliance larger (softer), which is why
/// XPBD properly decouples stiffness from substep count.
pub fn effective_compliance(alpha: f64, dt: f64) -> f64 {
    alpha / (dt * dt)
}
/// Convert a spring stiffness k (N/m) to XPBD compliance alpha (m^2/N).
///
/// alpha = 1/k.
pub fn stiffness_to_compliance(k: f64) -> f64 {
    if k.abs() < 1e-30 { 0.0 } else { 1.0 / k }
}
/// Convert XPBD compliance alpha (m^2/N) to spring stiffness k (N/m).
pub fn compliance_to_stiffness(alpha: f64) -> f64 {
    if alpha.abs() < 1e-30 {
        f64::INFINITY
    } else {
        1.0 / alpha
    }
}
/// Apply a global position correction by shifting all particles toward their
/// centre of mass.
///
/// This is useful for removing drift in simulations with only relative
/// constraints.
pub fn remove_com_drift(positions: &mut [Vec3], inv_masses: &[f64]) {
    if positions.is_empty() {
        return;
    }
    let mut total_mass = 0.0;
    let mut com = Vec3::zeros();
    for (pos, &w) in positions.iter().zip(inv_masses.iter()) {
        let m = if w > 0.0 { 1.0 / w } else { 0.0 };
        com += *pos * m;
        total_mass += m;
    }
    if total_mass < 1e-30 {
        return;
    }
    com /= total_mass;
    for pos in positions.iter_mut() {
        *pos -= com;
    }
}
/// Clamp particle velocities to a maximum magnitude.
///
/// This prevents numerical explosion in unstable simulations.
pub fn clamp_velocities(velocities: &mut [Vec3], max_speed: f64) {
    let max2 = max_speed * max_speed;
    for v in velocities.iter_mut() {
        let s2 = v.norm_squared();
        if s2 > max2 {
            let s = s2.sqrt();
            *v *= max_speed / s;
        }
    }
}
/// Apply global velocity damping to all particles.
///
/// `factor` should be in \[0, 1\]; 0 = no damping, 1 = full damping (freeze).
pub fn apply_global_damping(velocities: &mut [Vec3], factor: f64) {
    let keep = (1.0 - factor).clamp(0.0, 1.0);
    for v in velocities.iter_mut() {
        *v *= keep;
    }
}
/// Polar decomposition of a 3×3 matrix A = R * S.
///
/// Returns the orthogonal factor R using iterative method:
/// `R_{n+1} = 0.5 * (R_n + R_n^{-T})`
///
/// Convergence is guaranteed for non-singular matrices.
/// Falls back to Gram-Schmidt if the matrix is degenerate.
pub(super) fn polar_decompose_3x3(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mat3_transpose = |m: [[f64; 3]; 3]| -> [[f64; 3]; 3] {
        let mut t = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] = m[j][i];
            }
        }
        t
    };
    let mat3_mul = |a: [[f64; 3]; 3], b: [[f64; 3]; 3]| -> [[f64; 3]; 3] {
        let mut c = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        c
    };
    let mat3_add_scaled = |a: [[f64; 3]; 3], b: [[f64; 3]; 3], s: f64| -> [[f64; 3]; 3] {
        let mut c = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                c[i][j] = s * a[i][j] + s * b[i][j];
            }
        }
        c
    };
    let mat3_det = |m: [[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let mat3_inv = |m: [[f64; 3]; 3]| -> Option<[[f64; 3]; 3]> {
        let d = mat3_det(m);
        if d.abs() < 1e-20 {
            return None;
        }
        let inv_d = 1.0 / d;
        Some([
            [
                (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_d,
                (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_d,
                (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_d,
            ],
            [
                (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_d,
                (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_d,
                (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_d,
            ],
            [
                (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_d,
                (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_d,
                (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_d,
            ],
        ])
    };
    let mut r = a;
    let det = mat3_det(r);
    if det.abs() < 1e-20 {
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }
    for _ in 0..20 {
        let r_inv = match mat3_inv(r) {
            Some(inv) => inv,
            None => break,
        };
        let r_inv_t = mat3_transpose(r_inv);
        let r_new = mat3_add_scaled(r, r_inv_t, 0.5);
        let mut err = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                let d = r_new[i][j] - r[i][j];
                err += d * d;
            }
        }
        r = r_new;
        if err.sqrt() < 1e-10 {
            break;
        }
    }
    if mat3_det(r) < 0.0 {
        r[0][2] = -r[0][2];
        r[1][2] = -r[1][2];
        r[2][2] = -r[2][2];
    }
    let _ = mat3_mul;
    col_gs_3x3(r)
}
/// Column-Gram-Schmidt orthogonalization of a 3×3 matrix.
///
/// Returns a rotation matrix R (3×3 row-major) whose columns are orthonormal.
pub(super) fn col_gs_3x3(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut c0 = [a[0][0], a[1][0], a[2][0]];
    let mut c1 = [a[0][1], a[1][1], a[2][1]];

    let norm3 = |v: [f64; 3]| -> f64 { (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt() };
    let dot3 = |u: [f64; 3], v: [f64; 3]| -> f64 { u[0] * v[0] + u[1] * v[1] + u[2] * v[2] };
    let scale3 = |v: [f64; 3], s: f64| -> [f64; 3] { [v[0] * s, v[1] * s, v[2] * s] };
    let sub3 = |u: [f64; 3], v: [f64; 3]| -> [f64; 3] { [u[0] - v[0], u[1] - v[1], u[2] - v[2]] };
    let n0 = norm3(c0);
    c0 = if n0 > 1e-12 {
        scale3(c0, 1.0 / n0)
    } else {
        [1.0, 0.0, 0.0]
    };
    c1 = sub3(c1, scale3(c0, dot3(c0, c1)));
    let n1 = norm3(c1);
    c1 = if n1 > 1e-12 {
        scale3(c1, 1.0 / n1)
    } else {
        if c0[0].abs() < 0.9 {
            let t = sub3([1.0, 0.0, 0.0], scale3(c0, dot3(c0, [1.0, 0.0, 0.0])));
            let n = norm3(t).max(1e-12);
            scale3(t, 1.0 / n)
        } else {
            let t = sub3([0.0, 1.0, 0.0], scale3(c0, dot3(c0, [0.0, 1.0, 0.0])));
            let n = norm3(t).max(1e-12);
            scale3(t, 1.0 / n)
        }
    };
    let c2 = [
        c0[1] * c1[2] - c0[2] * c1[1],
        c0[2] * c1[0] - c0[0] * c1[2],
        c0[0] * c1[1] - c0[1] * c1[0],
    ];
    [
        [c0[0], c1[0], c2[0]],
        [c0[1], c1[1], c2[1]],
        [c0[2], c1[2], c2[2]],
    ]
}
/// Apply a rigid-body coupling constraint between a fixed rigid-body anchor
/// and a soft particle.
///
/// Models a pin joint between a rigid body (kinematically driven, infinite
/// mass) and a single XPBD soft particle.  The correction is one-sided:
/// only the soft particle is displaced.
///
/// # Arguments
/// * `rigid_pos`   – world-space position of the rigid anchor point.
/// * `soft_pos`    – mutable position of the soft particle.
/// * `inv_mass`    – inverse mass of the soft particle (0 → pinned, no effect).
/// * `rest_length` – rest distance between the anchor and the particle.
/// * `compliance`  – XPBD compliance α (m²/N).  0 = rigid coupling.
/// * `dt`          – current sub-step time (s).
pub fn apply_rigid_body_coupling(
    rigid_pos: &Vec3,
    soft_pos: &mut Vec3,
    inv_mass: f64,
    rest_length: f64,
    compliance: f64,
    dt: f64,
) {
    if inv_mass < 1e-30 {
        return;
    }
    let diff = *soft_pos - *rigid_pos;
    let dist = diff.norm();
    if dist < 1e-15 {
        return;
    }
    let c = dist - rest_length;
    let alpha_tilde = compliance / (dt * dt);
    let delta_lambda = -c / (inv_mass + alpha_tilde);
    let n = diff / dist;
    *soft_pos += n * (inv_mass * delta_lambda);
}
/// Compute the RMS (root-mean-square) residual of a set of constraint values.
///
/// Given a slice of scalar constraint evaluations C_i, returns
/// `||C||_rms = sqrt( Σ C_i² / n )`.
///
/// An empty slice returns 0.
pub fn compute_constraint_residual(constraint_values: &[f64]) -> f64 {
    let n = constraint_values.len();
    if n == 0 {
        return 0.0;
    }
    let sum_sq: f64 = constraint_values.iter().map(|c| c * c).sum();
    (sum_sq / n as f64).sqrt()
}
/// Compute strain-dependent (adaptive) XPBD compliance.
///
/// Stiff constraints can cause numerical issues under large deformations.
/// This function increases the compliance smoothly once the local strain
/// exceeds a reference threshold, up to a maximum cap.
///
/// Formula (soft ramp beyond the elastic regime):
/// ```text
/// α_adapt = clamp( α_base * (1 + softening * strain²), α_base, α_max )
/// ```
///
/// # Arguments
/// * `base_compliance`  – nominal compliance α at zero strain.
/// * `strain`           – dimensionless local strain ε = (l − l₀) / l₀.
/// * `max_compliance`   – upper bound on the returned compliance.
/// * `softening`        – dimensionless softening factor (≥ 0).  Higher
///   values ramp up compliance more aggressively.
pub fn adaptive_compliance(
    base_compliance: f64,
    strain: f64,
    max_compliance: f64,
    softening: f64,
) -> f64 {
    let adapted = base_compliance * (1.0 + softening * strain * strain);
    adapted.clamp(base_compliance, max_compliance.max(base_compliance))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::xpbd::SelfCollisionPair;
    use crate::xpbd::SubstepConfig;
    use crate::xpbd::XpbdBendingConstraint;
    use crate::xpbd::XpbdBody;
    use crate::xpbd::XpbdDampingConstraint;
    use crate::xpbd::XpbdDistanceConstraint;
    use crate::xpbd::XpbdSelfCollision;
    use crate::xpbd::XpbdShapeMatchingConstraint;
    use crate::xpbd::XpbdSolverStats;
    use crate::xpbd::XpbdVolumeConstraint;
    #[test]
    fn test_xpbd_distance_constraint_evaluate() {
        let positions = vec![Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0)];
        let c = XpbdDistanceConstraint::new(0, 1, 1.0, 0.0);
        let val = c.evaluate(&positions);
        assert!(val > 0.0, "Stretched constraint must have C > 0, got {val}");
    }
    #[test]
    fn test_xpbd_distance_constraint_solve() {
        let mut positions = vec![Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let rest = 1.0;
        let mut c = XpbdDistanceConstraint::new(0, 1, rest, 0.0);
        let before = c.evaluate(&positions).abs();
        c.solve(&mut positions, &inv_masses, 1.0 / 60.0);
        let after = c.evaluate(&positions).abs();
        assert!(
            after < before,
            "Constraint residual should decrease after solve: before={before}, after={after}"
        );
    }
    #[test]
    fn test_xpbd_volume_constraint_evaluate() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 2.0),
        ];
        let rest_volume = 1.0 / 6.0;
        let vc = XpbdVolumeConstraint::new([0, 1, 2, 3], rest_volume, 0.0);
        let val = vc.evaluate(&positions);
        assert!(
            val.abs() > 1e-10,
            "Displaced tet must have nonzero C, got {val}"
        );
    }
    #[test]
    fn test_xpbd_body_gravity_freefall() {
        let mut body = XpbdBody::new(1);
        body.add_particle(Vec3::new(0.0, 10.0, 0.0), 1.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            body.step(dt, gravity, 1, 1);
        }
        let y = body.positions[0].y;
        assert!(
            y < 9.0,
            "Free particle should have fallen under gravity after ~1 s, y={y}"
        );
    }
    #[test]
    fn test_xpbd_body_pinned_constraint() {
        let mut body = XpbdBody::new(2);
        let _pin = body.add_particle(Vec3::zeros(), 0.0);
        let free = body.add_particle(Vec3::new(1.0, 0.0, 0.0), 1.0);
        body.pin_particle(0);
        body.add_distance_constraint(0, 1, 0.0);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 1.0 / 60.0;
        for _ in 0..120 {
            body.step(dt, gravity, 4, 4);
        }
        let dist = body.positions[free].norm();
        assert!(
            (dist - 1.0).abs() < 0.1,
            "Free particle should remain ~1.0 from pin after constraint, dist={dist}"
        );
    }
    #[test]
    fn test_xpbd_compliance_zero_rigid() {
        let mut positions = vec![Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let rest = 1.0;
        let mut c = XpbdDistanceConstraint::new(0, 1, rest, 0.0);
        for _ in 0..20 {
            c.reset_lambda();
            c.solve(&mut positions, &inv_masses, 1.0 / 60.0);
        }
        let residual = c.evaluate(&positions).abs();
        assert!(
            residual < 1e-6,
            "Rigid constraint should be nearly perfectly satisfied, residual={residual}"
        );
    }
    #[test]
    fn test_xpbd_compliance_soft() {
        let rest = 1.0;
        let dt = 1.0 / 60.0;
        let n_iter = 5;
        let mut pos_rigid = vec![Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let mut c_rigid = XpbdDistanceConstraint::new(0, 1, rest, 0.0);
        c_rigid.reset_lambda();
        for _ in 0..n_iter {
            c_rigid.solve(&mut pos_rigid, &inv_masses, dt);
        }
        let residual_rigid = c_rigid.evaluate(&pos_rigid).abs();
        let mut pos_soft = vec![Vec3::zeros(), Vec3::new(3.0, 0.0, 0.0)];
        let mut c_soft = XpbdDistanceConstraint::new(0, 1, rest, 1.0);
        c_soft.reset_lambda();
        for _ in 0..n_iter {
            c_soft.solve(&mut pos_soft, &inv_masses, dt);
        }
        let residual_soft = c_soft.evaluate(&pos_soft).abs();
        assert!(
            residual_soft > residual_rigid,
            "Soft constraint should have larger residual than rigid: soft={residual_soft}, rigid={residual_rigid}"
        );
    }
    #[test]
    fn test_effective_compliance() {
        let alpha = 0.01;
        let dt = 1.0 / 60.0;
        let at = effective_compliance(alpha, dt);
        let expected = alpha / (dt * dt);
        assert!(
            (at - expected).abs() < 1e-10,
            "effective_compliance={at}, expected={expected}"
        );
    }
    #[test]
    fn test_stiffness_compliance_roundtrip() {
        let k = 1000.0;
        let alpha = stiffness_to_compliance(k);
        let k_back = compliance_to_stiffness(alpha);
        assert!((k - k_back).abs() < 1e-6, "k={k}, round-tripped={k_back}");
    }
    #[test]
    fn test_stiffness_zero() {
        assert_eq!(stiffness_to_compliance(0.0), 0.0);
    }
    #[test]
    fn test_compliance_zero_gives_inf() {
        assert!(compliance_to_stiffness(0.0).is_infinite());
    }
    #[test]
    fn test_damping_constraint_reduces_velocity() {
        let mut positions = vec![Vec3::zeros(), Vec3::new(2.0, 0.0, 0.0)];
        let prev_positions = vec![Vec3::new(0.1, 0.0, 0.0), Vec3::new(1.8, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let dt = 1.0 / 60.0;
        let dc = XpbdDampingConstraint::new(0, 1, 0.5);
        let before = (positions[1] - positions[0]).norm();
        dc.apply(&mut positions, &prev_positions, &inv_masses, dt);
        let after = (positions[1] - positions[0]).norm();
        assert!(
            (after - before).abs() > 1e-10,
            "Damping should modify positions"
        );
    }
    #[test]
    fn test_clamp_velocities() {
        let mut vels = vec![Vec3::new(100.0, 0.0, 0.0), Vec3::new(0.0, 5.0, 0.0)];
        clamp_velocities(&mut vels, 10.0);
        assert!(
            (vels[0].norm() - 10.0).abs() < 1e-6,
            "Should be clamped to 10, got {}",
            vels[0].norm()
        );
        assert!(
            (vels[1].norm() - 5.0).abs() < 1e-6,
            "Should remain 5, got {}",
            vels[1].norm()
        );
    }
    #[test]
    fn test_global_damping_zero() {
        let mut vels = vec![Vec3::new(1.0, 2.0, 3.0)];
        apply_global_damping(&mut vels, 0.0);
        assert!((vels[0].x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_global_damping_full() {
        let mut vels = vec![Vec3::new(1.0, 2.0, 3.0)];
        apply_global_damping(&mut vels, 1.0);
        assert!(vels[0].norm() < 1e-10);
    }
    #[test]
    fn test_step_with_stats() {
        let mut body = XpbdBody::new(2);
        body.add_particle(Vec3::zeros(), 0.0);
        body.add_particle(Vec3::new(2.0, 0.0, 0.0), 1.0);
        body.add_distance_constraint(0, 1, 0.0);
        let stats = body.step_with_stats(1.0 / 60.0, Vec3::new(0.0, -9.81, 0.0), 2, 3);
        assert!(stats.total_projections > 0);
        assert_eq!(stats.substeps, 2);
        assert_eq!(stats.iterations_per_substep, 3);
    }
    #[test]
    fn test_particle_and_constraint_count() {
        let mut body = XpbdBody::new(3);
        body.add_particle(Vec3::zeros(), 1.0);
        body.add_particle(Vec3::new(1.0, 0.0, 0.0), 1.0);
        body.add_particle(Vec3::new(2.0, 0.0, 0.0), 1.0);
        body.add_distance_constraint(0, 1, 0.0);
        body.add_distance_constraint(1, 2, 0.0);
        assert_eq!(body.particle_count(), 3);
        assert_eq!(body.constraint_count(), 2);
    }
    #[test]
    fn test_centre_of_mass() {
        let mut body = XpbdBody::new(2);
        body.add_particle(Vec3::zeros(), 1.0);
        body.add_particle(Vec3::new(4.0, 0.0, 0.0), 1.0);
        let com = body.centre_of_mass();
        assert!((com.x - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_max_velocity_stationary() {
        let mut body = XpbdBody::new(2);
        body.add_particle(Vec3::zeros(), 1.0);
        body.add_particle(Vec3::new(1.0, 0.0, 0.0), 1.0);
        assert!(body.max_velocity() < 1e-10);
    }
    #[test]
    fn test_max_velocity_after_step() {
        let mut body = XpbdBody::new(1);
        body.add_particle(Vec3::new(0.0, 10.0, 0.0), 1.0);
        body.step(1.0 / 60.0, Vec3::new(0.0, -9.81, 0.0), 1, 1);
        assert!(body.max_velocity() > 0.0);
    }
    #[test]
    fn test_bending_constraint_rest_angle() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let bc = XpbdBendingConstraint::new(0, 1, 2, 3, &positions, 0.0);
        let c = bc.evaluate(&positions);
        assert!(c.abs() < 1e-10, "At rest, C should be zero, got {c}");
    }
    #[test]
    fn test_volume_constraint_unit_tet() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let vc = XpbdVolumeConstraint::new([0, 1, 2, 3], 1.0 / 6.0, 0.0);
        let vol = vc.current_volume(&positions);
        assert!((vol - 1.0 / 6.0).abs() < 1e-10, "vol={vol}");
    }
    #[test]
    fn test_remove_com_drift() {
        let mut positions = vec![Vec3::new(10.0, 0.0, 0.0), Vec3::new(12.0, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        remove_com_drift(&mut positions, &inv_masses);
        let com = (positions[0] + positions[1]) / 2.0;
        assert!(com.norm() < 1e-10, "COM should be at origin, got {:?}", com);
    }
    #[test]
    fn test_solver_stats_default() {
        let stats = XpbdSolverStats::default();
        assert_eq!(stats.total_projections, 0);
        assert_eq!(stats.max_residual, 0.0);
    }
    #[test]
    fn test_shape_matching_at_rest_no_change() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ];
        let masses = vec![1.0_f64; 3];
        let sm = XpbdShapeMatchingConstraint::new(vec![0, 1, 2], &positions, &masses, 1.0);
        let inv_masses = vec![1.0_f64; 3];
        let mut pos = positions.clone();
        sm.solve(&mut pos, &inv_masses);
        for i in 0..3 {
            let err = (pos[i] - positions[i]).norm();
            assert!(
                err < 1e-6,
                "At rest, shape-matching should not change positions: err={err}"
            );
        }
    }
    #[test]
    fn test_shape_matching_corrects_deformation() {
        let rest_positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
        ];
        let masses = vec![1.0_f64; 3];
        let sm = XpbdShapeMatchingConstraint::new(vec![0, 1, 2], &rest_positions, &masses, 1.0);
        let inv_masses = vec![1.0_f64; 3];
        let mut pos = rest_positions.clone();
        pos[2] = Vec3::new(0.5, 5.0, 0.0);
        let before_err = {
            let diff = pos[2] - rest_positions[2];
            diff.norm()
        };
        sm.solve(&mut pos, &inv_masses);
        let after_err = (pos[2] - rest_positions[2]).norm();
        assert!(
            after_err < before_err,
            "Shape-matching should reduce deformation error: before={before_err}, after={after_err}"
        );
    }
    #[test]
    fn test_shape_matching_zero_stiffness() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
        ];
        let masses = vec![1.0_f64; 3];
        let sm = XpbdShapeMatchingConstraint::new(vec![0, 1, 2], &positions, &masses, 0.0);
        let inv_masses = vec![1.0_f64; 3];
        let mut pos = positions.clone();
        pos[1] = Vec3::new(3.0, 0.0, 0.0);
        let before = pos[1];
        sm.solve(&mut pos, &inv_masses);
        assert!(
            (pos[1] - before).norm() < 1e-10,
            "Zero stiffness should not change positions"
        );
    }
    #[test]
    fn test_shape_matching_pinned_particle() {
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
        ];
        let masses = vec![1.0_f64; 3];
        let sm = XpbdShapeMatchingConstraint::new(vec![0, 1, 2], &positions, &masses, 1.0);
        let inv_masses = vec![0.0_f64, 1.0, 1.0];
        let mut pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(3.0, 0.0, 0.0),
            Vec3::new(0.5, 1.0, 0.0),
        ];
        let pinned_before = pos[0];
        sm.solve(&mut pos, &inv_masses);
        assert!(
            (pos[0] - pinned_before).norm() < 1e-10,
            "Pinned particle should not be moved by shape-matching"
        );
    }
    #[test]
    fn test_self_collision_detect_overlapping() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.1, 0.0, 0.0)];
        let sc = XpbdSelfCollision::new_uniform(2, 0.1, 0.0);
        let pairs = sc.detect_pairs(&positions);
        assert_eq!(pairs.len(), 1, "Should detect one collision");
    }
    #[test]
    fn test_self_collision_detect_separated() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let sc = XpbdSelfCollision::new_uniform(2, 0.1, 0.0);
        let pairs = sc.detect_pairs(&positions);
        assert!(pairs.is_empty(), "No collision for separated particles");
    }
    #[test]
    fn test_self_collision_resolve_pushes_apart() {
        let mut positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.1, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let sc = XpbdSelfCollision::new_uniform(2, 0.1, 0.0);
        let before_dist = (positions[1] - positions[0]).norm();
        sc.resolve(&mut positions, &inv_masses, 1.0 / 60.0);
        let after_dist = (positions[1] - positions[0]).norm();
        assert!(
            after_dist > before_dist,
            "Particles should be pushed apart: before={before_dist}, after={after_dist}"
        );
    }
    #[test]
    fn test_self_collision_pinned_not_moved() {
        let mut positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.05, 0.0, 0.0)];
        let inv_masses = vec![0.0_f64, 1.0];
        let sc = XpbdSelfCollision::new_uniform(2, 0.1, 0.0);
        let before = positions[0];
        sc.resolve(&mut positions, &inv_masses, 1.0 / 60.0);
        assert!(
            (positions[0] - before).norm() < 1e-10,
            "Pinned particle should not move during self-collision"
        );
    }
    #[test]
    fn test_self_collision_velocity_response() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.15, 0.0, 0.0)];
        let mut velocities = vec![Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0)];
        let inv_masses = vec![1.0, 1.0];
        let sc = XpbdSelfCollision::new_uniform(2, 0.1, 0.0);
        let pair = SelfCollisionPair {
            i: 0,
            j: 1,
            contact_dist: 0.2,
        };
        sc.velocity_response(&pair, &positions, &mut velocities, &inv_masses, 0.5);
        let v_rel = velocities[1] - velocities[0];
        let n = Vec3::new(1.0, 0.0, 0.0);
        let vn = v_rel.dot(&n);
        assert!(
            vn > 0.0,
            "After response, particles should be separating: vn={vn}"
        );
    }
    #[test]
    fn test_substep_count_scales_with_velocity() {
        let cfg = SubstepConfig::default_config();
        let n_slow = cfg.substeps_for_velocity(1.0 / 60.0, 1.0, 0.1);
        let n_fast = cfg.substeps_for_velocity(1.0 / 60.0, 100.0, 0.1);
        assert!(
            n_fast > n_slow,
            "More substeps needed at higher velocity: slow={n_slow}, fast={n_fast}"
        );
    }
    #[test]
    fn test_substep_count_min_clamp() {
        let mut cfg = SubstepConfig::default_config();
        cfg.min_substeps = 3;
        let n = cfg.substeps_for_velocity(1.0 / 60.0, 0.0001, 1.0);
        assert!(n >= 3, "Should be at least min_substeps: {n}");
    }
    #[test]
    fn test_substep_count_max_clamp() {
        let mut cfg = SubstepConfig::default_config();
        cfg.max_substeps = 8;
        let n = cfg.substeps_for_velocity(1.0 / 60.0, 1e9, 1e-6);
        assert!(n <= 8, "Should not exceed max_substeps: {n}");
    }
    #[test]
    fn test_substep_residual_ok() {
        let cfg = SubstepConfig::default_config();
        assert!(cfg.residual_ok(1e-5), "1e-5 should be ok");
        assert!(!cfg.residual_ok(1e-3), "1e-3 should not be ok");
    }
    #[test]
    fn test_substep_dt() {
        let cfg = SubstepConfig::default_config();
        let h = cfg.substep_dt(1.0 / 60.0, 4);
        assert!((h - 1.0 / 240.0).abs() < 1e-12, "substep_dt mismatch: {h}");
    }
    #[test]
    fn test_total_projections_formula() {
        let cfg = SubstepConfig {
            min_substeps: 1,
            max_substeps: 32,
            target_residual: 1e-4,
            iterations: 5,
        };
        let tp = cfg.total_projections(4, 10);
        assert_eq!(tp, 200, "total_projections should be 200, got {tp}");
    }
    #[test]
    fn test_apply_rigid_body_coupling_moves_particle() {
        let rigid_pos = Vec3::new(0.0, 0.0, 0.0);
        let mut soft_pos = Vec3::new(2.0, 0.0, 0.0);
        let inv_mass = 1.0;
        let compliance = 0.0;
        let rest = 0.0;
        let dt = 1.0 / 60.0;
        apply_rigid_body_coupling(&rigid_pos, &mut soft_pos, inv_mass, rest, compliance, dt);
        assert!(
            soft_pos.norm() < 2.0,
            "Soft particle should move closer to rigid body: pos={:?}",
            soft_pos
        );
    }
    #[test]
    fn test_apply_rigid_body_coupling_compliance_effect() {
        let rigid_pos = Vec3::new(0.0, 0.0, 0.0);
        let mut soft_pos_rigid = Vec3::new(3.0, 0.0, 0.0);
        let mut soft_pos_soft = Vec3::new(3.0, 0.0, 0.0);
        let dt = 1.0 / 60.0;
        apply_rigid_body_coupling(&rigid_pos, &mut soft_pos_rigid, 1.0, 0.0, 0.0, dt);
        apply_rigid_body_coupling(&rigid_pos, &mut soft_pos_soft, 1.0, 0.0, 1.0, dt);
        assert!(
            soft_pos_rigid.norm() < soft_pos_soft.norm(),
            "Rigid coupling should move particle more than soft: rigid={}, soft={}",
            soft_pos_rigid.norm(),
            soft_pos_soft.norm()
        );
    }
    #[test]
    fn test_compute_constraint_residual_zero_at_rest() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.0, 0.0)];
        let c = XpbdDistanceConstraint::new(0, 1, 1.0, 0.0);
        let residual = compute_constraint_residual(&[c.evaluate(&positions)]);
        assert!(
            residual < 1e-12,
            "Residual should be ~0 for satisfied constraint, got {residual}"
        );
    }
    #[test]
    fn test_compute_constraint_residual_nonzero() {
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let c = XpbdDistanceConstraint::new(0, 1, 1.0, 0.0);
        let residual = compute_constraint_residual(&[c.evaluate(&positions)]);
        assert!(
            residual > 0.5,
            "Residual should reflect violation, got {residual}"
        );
    }
    #[test]
    fn test_compute_constraint_residual_multiple() {
        let vals = vec![0.1, 0.2, 0.3, 0.0];
        let res = compute_constraint_residual(&vals);
        let expected = (0.01_f64 + 0.04 + 0.09).sqrt() / (4.0_f64.sqrt());
        assert!(
            (res - expected).abs() < 1e-10,
            "Residual mismatch: {res} vs {expected}"
        );
    }
    #[test]
    fn test_adaptive_compliance_increases_with_strain() {
        let base = 1e-4;
        let low_strain = adaptive_compliance(base, 0.01, 1.0, 1.0);
        let high_strain = adaptive_compliance(base, 0.5, 1.0, 1.0);
        assert!(
            high_strain > low_strain,
            "Compliance should increase with strain: low={low_strain}, high={high_strain}"
        );
    }
    #[test]
    fn test_adaptive_compliance_zero_strain() {
        let base = 1e-4;
        let result = adaptive_compliance(base, 0.0, 1.0, 1.0);
        assert!(
            (result - base).abs() < 1e-12,
            "At zero strain, compliance should equal base: {result} vs {base}"
        );
    }
    #[test]
    fn test_adaptive_compliance_clamped_at_max() {
        let base = 1e-4;
        let max_compliance = 1e-2;
        let result = adaptive_compliance(base, 100.0, max_compliance, 1.0);
        assert!(
            result <= max_compliance + 1e-15,
            "Compliance must not exceed max: {result} vs {max_compliance}"
        );
    }
    #[test]
    fn test_apply_rigid_body_coupling_rest_length() {
        let rigid_pos = Vec3::new(0.0, 0.0, 0.0);
        let mut soft_pos = Vec3::new(1.5, 0.0, 0.0);
        let dt = 1.0 / 60.0;
        apply_rigid_body_coupling(&rigid_pos, &mut soft_pos, 1.0, 1.5, 0.0, dt);
        assert!(
            (soft_pos.norm() - 1.5).abs() < 1e-8,
            "At rest length, no correction needed: pos={}",
            soft_pos.norm()
        );
    }
    #[test]
    fn test_compute_constraint_residual_empty() {
        let res = compute_constraint_residual(&[]);
        assert_eq!(res, 0.0, "Empty residual list should return 0");
    }
    #[test]
    fn test_col_gs_3x3_orthonormal() {
        let a = [[1.0, 0.5, 0.1], [0.2, 1.0, 0.3], [0.0, 0.1, 1.0]];
        let r = col_gs_3x3(a);
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| r[k][i] * r[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - expected).abs() < 1e-8,
                    "R^T*R[{i}][{j}] should be {expected}, got {dot}"
                );
            }
        }
    }
}
