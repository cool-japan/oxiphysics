//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
/// Dot product of two 3-vectors.
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Length of a 3-vector.
pub(super) fn len3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// Normalize a 3-vector. Returns zero vector if length < eps.
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        [v[0] / l, v[1] / l, v[2] / l]
    }
}
/// Add two 3-vectors.
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two 3-vectors (a - b).
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a 3-vector.
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// 3x3 matrix-vector multiply.
pub(super) fn mat_vec3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// 3x3 matrix multiply.
pub(super) fn mat_mul3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                r[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    r
}
/// 3x3 identity matrix.
pub(super) fn identity3() -> [[f64; 3]; 3] {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}
/// Skew-symmetric matrix from a 3-vector.
pub(super) fn skew3(v: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]]
}
/// Matrix + matrix (3x3).
pub(super) fn mat_add3(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][j] + b[i][j];
        }
    }
    r
}
/// Matrix scale (3x3).
pub(super) fn mat_scale3(m: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    let mut r = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = m[i][j] * s;
        }
    }
    r
}
/// Determinant of a 3x3 matrix.
pub(super) fn det3(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
/// Inverse of a 3x3 matrix. Returns identity if singular.
pub(super) fn inv3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let d = det3(m);
    if d.abs() < 1e-30 {
        return identity3();
    }
    let inv_d = 1.0 / d;
    [
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
    ]
}
/// Outer product of two 3-vectors.
pub(super) fn outer3(a: [f64; 3], b: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}
/// Exponential map from so(3) to SO(3) using Rodrigues' formula.
///
/// exp(theta * hat(a)) = I + sin(theta)*hat(a) + (1-cos(theta))*hat(a)^2
///
/// where hat(a) is the skew-symmetric matrix of axis a, theta = |v|.
pub fn exp_so3(v: [f64; 3]) -> [[f64; 3]; 3] {
    let theta = len3(v);
    if theta < 1e-14 {
        return mat_add3(identity3(), skew3(v));
    }
    let axis = scale3(v, 1.0 / theta);
    let k = skew3(axis);
    let k2 = mat_mul3(k, k);
    let s = theta.sin();
    let c = 1.0 - theta.cos();
    mat_add3(mat_add3(identity3(), mat_scale3(k, s)), mat_scale3(k2, c))
}
/// Logarithmic map from SO(3) to so(3).
///
/// Returns the rotation vector v such that exp(hat(v)) = R.
pub fn log_so3(r: [[f64; 3]; 3]) -> [f64; 3] {
    let trace = r[0][0] + r[1][1] + r[2][2];
    let cos_theta = (trace - 1.0) * 0.5;
    let cos_clamped = cos_theta.clamp(-1.0, 1.0);
    let theta = cos_clamped.acos();
    if theta.abs() < 1e-10 {
        return [
            0.5 * (r[2][1] - r[1][2]),
            0.5 * (r[0][2] - r[2][0]),
            0.5 * (r[1][0] - r[0][1]),
        ];
    }
    let factor = theta / (2.0 * theta.sin());
    [
        factor * (r[2][1] - r[1][2]),
        factor * (r[0][2] - r[2][0]),
        factor * (r[1][0] - r[0][1]),
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::variational_constraints::BackwardErrorAnalysis;
    use crate::variational_constraints::ConstrainedVariationalIntegrator;
    use crate::variational_constraints::DiscreteEulerLagrange;
    use crate::variational_constraints::DiscreteLagrangian;
    use crate::variational_constraints::DiscreteLegendreTransform;
    use crate::variational_constraints::DiscreteNoether;
    use crate::variational_constraints::DiscreteNullSpace;
    use crate::variational_constraints::EnergyMonitor;
    use crate::variational_constraints::FourthOrderVariational;
    use crate::variational_constraints::HolonomicConstraint;
    use crate::variational_constraints::LieGroupIntegrator;
    use crate::variational_constraints::LieGroupState;
    use crate::variational_constraints::MomentumMap;
    use crate::variational_constraints::NBodyVariational;
    use crate::variational_constraints::SymplecticPRK;
    use crate::variational_constraints::SymplecticVerlet;
    use crate::variational_constraints::VariationalCollision;
    use crate::variational_constraints::VariationalMidpoint;
    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < 1e-10);
    }
    #[test]
    fn test_cross3() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_normalize3() {
        let n = normalize3([3.0, 4.0, 0.0]);
        assert!((len3(n) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_skew_antisymmetric() {
        let v = [1.0, 2.0, 3.0];
        let s = skew3(v);
        for (i, row) in s.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val + s[j][i]).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_inv3_identity() {
        let inv = inv3(identity3());
        for (i, row) in inv.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_discrete_lagrangian_free_particle() {
        let dl = DiscreteLagrangian::new(1.0, 0.0, 0.0);
        let q0 = [0.0; 3];
        let q1 = [1.0, 0.0, 0.0];
        let h = 1.0;
        let ld = dl.evaluate(q0, q1, h);
        assert!((ld - 0.5).abs() < 1e-6, "ld={ld}");
    }
    #[test]
    fn test_discrete_lagrangian_derivatives_consistency() {
        let dl = DiscreteLagrangian::new(2.0, 9.81, 1.0);
        let q0 = [1.0, 2.0, 3.0];
        let q1 = [1.1, 2.1, 3.1];
        let h = 0.01;
        let d1 = dl.d1(q0, q1, h);
        let d2 = dl.d2(q0, q1, h);
        for i in 0..3 {
            assert!(d1[i].is_finite(), "d1[{i}] not finite");
            assert!(d2[i].is_finite(), "d2[{i}] not finite");
        }
    }
    #[test]
    fn test_del_free_particle_uniform_motion() {
        let dl = DiscreteLagrangian::new(1.0, 0.0, 0.0);
        let del = DiscreteEulerLagrange::new(dl, 0.01);
        let q0 = [0.0, 0.0, 0.0];
        let q1 = [0.01, 0.0, 0.0];
        let (q2, iters) = del.solve(q0, q1);
        assert!((q2[0] - 0.02).abs() < 1e-6, "q2[0]={}", q2[0]);
        assert!(q2[1].abs() < 1e-6);
        assert!(iters < 20, "iters={iters}");
    }
    #[test]
    fn test_del_residual_at_solution() {
        let dl = DiscreteLagrangian::new(1.0, 0.0, 0.0);
        let del = DiscreteEulerLagrange::new(dl, 0.01);
        let q0 = [0.0; 3];
        let q1 = [0.01, 0.0, 0.0];
        let q2 = [0.02, 0.0, 0.0];
        let res = del.residual(q0, q1, q2);
        let res_norm = len3(res);
        assert!(res_norm < 1e-4, "residual={res_norm}");
    }
    #[test]
    fn test_symplectic_verlet_free_particle() {
        let sv = SymplecticVerlet::new(1.0, 0.01);
        let q0 = [0.0; 3];
        let p0 = [1.0, 0.0, 0.0];
        let zero_force = |_q: [f64; 3]| -> [f64; 3] { [0.0; 3] };
        let (q1, p1) = sv.step(q0, p0, &zero_force);
        assert!((q1[0] - 0.01).abs() < 1e-10);
        assert!((p1[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_symplectic_verlet_energy_conservation() {
        let sv = SymplecticVerlet::new(1.0, 0.001);
        let q0 = [1.0, 0.0, 0.0];
        let p0 = [0.0, 1.0, 0.0];
        let force = |q: [f64; 3]| -> [f64; 3] { [-q[0], -q[1], -q[2]] };
        let potential = |q: [f64; 3]| -> f64 { 0.5 * dot3(q, q) };
        let trajectory = sv.integrate(q0, p0, 1000, &force);
        let e0 = sv.hamiltonian(p0, potential(q0));
        let (qf, pf) = trajectory[1000];
        let ef = sv.hamiltonian(pf, potential(qf));
        let drift = (ef - e0).abs();
        assert!(drift < 1e-3, "energy drift={drift}");
    }
    #[test]
    fn test_holonomic_distance_evaluate() {
        let c = HolonomicConstraint::distance(0, 1, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let val = c.evaluate(&positions);
        assert!((val - 1.0).abs() < 1e-10, "val={val}");
    }
    #[test]
    fn test_holonomic_plane_evaluate() {
        let c = HolonomicConstraint::plane(0, [0.0, 0.0, 1.0], 0.0);
        let positions = vec![[1.0, 2.0, 0.5]];
        let val = c.evaluate(&positions);
        assert!((val - 0.5).abs() < 1e-10, "val={val}");
    }
    #[test]
    fn test_holonomic_sphere_evaluate() {
        let c = HolonomicConstraint::sphere(0, [0.0; 3], 1.0);
        let positions = vec![[1.5, 0.0, 0.0]];
        let val = c.evaluate(&positions);
        assert!((val - 0.5).abs() < 1e-10, "val={val}");
    }
    #[test]
    fn test_holonomic_jacobian_distance() {
        let c = HolonomicConstraint::distance(0, 1, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let jac = c.jacobian(&positions);
        assert_eq!(jac.len(), 2);
        assert!((jac[0].1[0] - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_null_space_projector_identity() {
        let dns = DiscreteNullSpace::new(0.2);
        let proj = dns.null_space_projector([0.0, 0.0, 1.0]);
        let v = mat_vec3(proj, [0.0, 0.0, 1.0]);
        assert!(len3(v) < 1e-10);
        let v2 = mat_vec3(proj, [1.0, 0.0, 0.0]);
        assert!((v2[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_project_distance_constraint() {
        let dns = DiscreteNullSpace::new(0.2);
        let qa = [0.0, 0.0, 0.0];
        let qb = [2.0, 0.0, 0.0];
        let (qa_new, qb_new) = dns.project_distance(qa, qb, 1.0, 1.0, 1.0);
        let dist = len3(sub3(qa_new, qb_new));
        assert!((dist - 1.0).abs() < 1e-8, "dist={dist}");
    }
    #[test]
    fn test_project_plane() {
        let dns = DiscreteNullSpace::new(0.2);
        let q = [1.0, 2.0, 3.0];
        let q_proj = dns.project_plane(q, [0.0, 0.0, 1.0], 0.0);
        assert!(q_proj[2].abs() < 1e-10);
        assert!((q_proj[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_exp_so3_identity() {
        let r = exp_so3([0.0; 3]);
        for (i, row) in r.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((val - expected).abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_exp_log_roundtrip() {
        let v = [0.1, 0.2, 0.3];
        let r = exp_so3(v);
        let v_back = log_so3(r);
        for i in 0..3 {
            assert!(
                (v[i] - v_back[i]).abs() < 1e-6,
                "v[{i}]={}, v_back={}",
                v[i],
                v_back[i]
            );
        }
    }
    #[test]
    fn test_lie_group_free_body_energy() {
        let state = LieGroupState::with_state(identity3(), [1.0, 0.5, 0.2], [2.0, 3.0, 4.0]);
        let integrator = LieGroupIntegrator::new(0.001);
        let e0 = state.kinetic_energy();
        let trajectory = integrator.integrate_free(&state, 100);
        let ef = trajectory[100].kinetic_energy();
        let drift = (ef - e0).abs() / e0;
        assert!(drift < 0.05, "relative energy drift={drift}");
    }
    #[test]
    fn test_lie_group_angular_momentum() {
        let state = LieGroupState::new([1.0, 1.0, 1.0]);
        let l = state.angular_momentum();
        assert!(len3(l) < 1e-10);
    }
    #[test]
    fn test_collision_plane_no_penetration() {
        let vc = VariationalCollision::new(1.0, 0.0);
        let q = [0.0, 0.0, 1.0];
        let p = [0.0, 0.0, -1.0];
        let (qc, pc) = vc.collide_plane(q, p, 1.0, [0.0, 0.0, 1.0], 0.0);
        assert!((qc[2] - 1.0).abs() < 1e-10);
        assert!((pc[2] - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_collision_plane_with_penetration() {
        let vc = VariationalCollision::new(1.0, 0.0);
        let q = [0.0, 0.0, -0.1];
        let p = [0.0, 0.0, -2.0];
        let (qc, pc) = vc.collide_plane(q, p, 1.0, [0.0, 0.0, 1.0], 0.0);
        assert!(qc[2] >= -1e-6, "qc[2]={}", qc[2]);
        assert!(pc[2] > 0.0, "pc[2]={}", pc[2]);
    }
    #[test]
    fn test_collision_particles_no_overlap() {
        let vc = VariationalCollision::new(1.0, 0.0);
        let qa = [0.0; 3];
        let qb = [5.0, 0.0, 0.0];
        let pa = [1.0, 0.0, 0.0];
        let pb = [-1.0, 0.0, 0.0];
        let (qa2, pa2, qb2, pb2) = vc.collide_particles(qa, pa, 1.0, qb, pb, 1.0, 1.0);
        assert!((qa2[0]).abs() < 1e-10);
        assert!((pa2[0] - 1.0).abs() < 1e-10);
        assert!((qb2[0] - 5.0).abs() < 1e-10);
        assert!((pb2[0] + 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_momentum_linear() {
        let mm = MomentumMap::new(2);
        let momenta = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let total = mm.linear_momentum(&momenta);
        assert!(len3(total) < 1e-10);
    }
    #[test]
    fn test_momentum_angular() {
        let mm = MomentumMap::new(1);
        let positions = vec![[1.0, 0.0, 0.0]];
        let momenta = vec![[0.0, 1.0, 0.0]];
        let l = mm.angular_momentum(&positions, &momenta);
        assert!((l[2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_center_of_mass() {
        let mm = MomentumMap::new(2);
        let positions = vec![[0.0; 3], [2.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let com = mm.center_of_mass(&positions, &masses);
        assert!((com[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_linear_momentum_conservation() {
        let mm = MomentumMap::new(2);
        let p_before = vec![[3.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let p_after = vec![[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert!(mm.is_linear_momentum_conserved(&p_before, &p_after, 1e-10));
    }
    #[test]
    fn test_noether_translation_invariance_free_particle() {
        let dl = DiscreteLagrangian::new(1.0, 0.0, 0.0);
        let noether = DiscreteNoether::new(1e-6);
        let q0 = [1.0, 2.0, 3.0];
        let q1 = [1.1, 2.0, 3.0];
        assert!(noether.check_translation_invariance(&dl, q0, q1, 0.01));
    }
    #[test]
    fn test_noether_translation_invariance_spring() {
        let dl = DiscreteLagrangian::new(1.0, 0.0, 1.0);
        let noether = DiscreteNoether::new(1e-6);
        let q0 = [1.0, 0.0, 0.0];
        let q1 = [1.1, 0.0, 0.0];
        assert!(!noether.check_translation_invariance(&dl, q0, q1, 0.01));
    }
    #[test]
    fn test_noether_verify_conservation() {
        let noether = DiscreteNoether::new(1e-8);
        let p1 = [1.0, 2.0, 3.0];
        let p2 = [1.0, 2.0, 3.0];
        assert!(noether.verify_conservation(p1, p2));
    }
    #[test]
    fn test_backward_error_modified_hamiltonian() {
        let bea = BackwardErrorAnalysis::new(0.01, 2);
        let h_mod = bea.modified_hamiltonian_harmonic(1.0, 1.0, 1.0, 1.0);
        assert!((h_mod - 1.005).abs() < 1e-10, "h_mod={h_mod}");
    }
    #[test]
    fn test_backward_error_order_estimate() {
        let bea = BackwardErrorAnalysis::new(0.01, 2);
        let order = bea.estimate_order(4.0, 1.0);
        assert!((order - 2.0).abs() < 1e-10, "order={order}");
    }
    #[test]
    fn test_backward_error_shadow_drift() {
        let bea = BackwardErrorAnalysis::new(0.01, 2);
        let energies = vec![10.0, 10.001, 9.999, 10.002];
        let drift = bea.shadow_hamiltonian_drift(&energies);
        assert!((drift - 0.002).abs() < 1e-10, "drift={drift}");
    }
    #[test]
    fn test_volume_preservation_perfect() {
        let bea = BackwardErrorAnalysis::new(0.01, 2);
        let err = bea.volume_preservation_error(1.0);
        assert!(err < 1e-15);
    }
    #[test]
    fn test_midpoint_free_particle() {
        let vm = VariationalMidpoint::new(1.0, 0.01);
        let q = [0.0; 3];
        let v = [1.0, 0.0, 0.0];
        let zero_force = |_q: [f64; 3]| -> [f64; 3] { [0.0; 3] };
        let (q1, v1) = vm.step(q, v, &zero_force);
        assert!((q1[0] - 0.01).abs() < 1e-10);
        assert!((v1[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_midpoint_discrete_action() {
        let vm = VariationalMidpoint::new(1.0, 0.01);
        let trajectory = vec![[0.0; 3], [0.01, 0.0, 0.0], [0.02, 0.0, 0.0]];
        let zero_potential = |_q: [f64; 3]| -> f64 { 0.0 };
        let action = vm.discrete_action(&trajectory, &zero_potential);
        assert!((action - 0.01).abs() < 1e-8, "action={action}");
    }
    #[test]
    fn test_constrained_integrator_distance() {
        let constraints = vec![HolonomicConstraint::distance(0, 1, 1.0)];
        let cvi = ConstrainedVariationalIntegrator::new(vec![1.0, 1.0], 0.01, constraints);
        let q_km1 = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let q_k = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let forces = vec![[0.0; 3], [0.0; 3]];
        let q_new = cvi.step(&q_km1, &q_k, &forces);
        let dist = len3(sub3(q_new[0], q_new[1]));
        assert!((dist - 1.0).abs() < 1e-6, "dist={dist}");
    }
    #[test]
    fn test_constrained_velocity() {
        let cvi = ConstrainedVariationalIntegrator::new(vec![1.0], 0.01, vec![]);
        let v = cvi.velocity([0.0; 3], [0.02, 0.0, 0.0]);
        assert!((v[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fourth_order_free_particle() {
        let fov = FourthOrderVariational::new(1.0, 0.01);
        let q0 = [0.0; 3];
        let p0 = [1.0, 0.0, 0.0];
        let zero_force = |_q: [f64; 3]| -> [f64; 3] { [0.0; 3] };
        let (q1, p1) = fov.step(q0, p0, &zero_force);
        assert!((q1[0] - 0.01).abs() < 1e-10, "q1[0]={}", q1[0]);
        assert!((p1[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fourth_order_energy_conservation() {
        let fov = FourthOrderVariational::new(1.0, 0.01);
        let q0 = [1.0, 0.0, 0.0];
        let p0 = [0.0, 1.0, 0.0];
        let force = |q: [f64; 3]| -> [f64; 3] { [-q[0], -q[1], -q[2]] };
        let trajectory = fov.integrate(q0, p0, 1000, &force);
        let energy = |q: [f64; 3], p: [f64; 3]| -> f64 { 0.5 * dot3(p, p) + 0.5 * dot3(q, q) };
        let e0 = energy(q0, p0);
        let (qf, pf) = trajectory[1000];
        let ef = energy(qf, pf);
        let drift = (ef - e0).abs();
        assert!(drift < 1e-4, "energy drift={drift}");
    }
    #[test]
    fn test_sprk_euler() {
        let sprk = SymplecticPRK::symplectic_euler(1.0, 0.01);
        let q = [0.0; 3];
        let p = [1.0, 0.0, 0.0];
        let zero_force = |_q: [f64; 3]| -> [f64; 3] { [0.0; 3] };
        let (q1, p1) = sprk.step(q, p, &zero_force);
        assert!((q1[0] - 0.01).abs() < 1e-10);
        assert!((p1[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_energy_monitor_drift() {
        let mut em = EnergyMonitor::new();
        em.record_energy(10.0);
        em.record_energy(10.01);
        em.record_energy(9.99);
        assert!((em.max_energy_drift() - 0.01).abs() < 1e-10);
    }
    #[test]
    fn test_energy_monitor_relative() {
        let mut em = EnergyMonitor::new();
        em.record_energy(100.0);
        em.record_energy(101.0);
        assert!((em.relative_energy_error() - 0.01).abs() < 1e-10);
    }
    #[test]
    fn test_energy_monitor_bounded() {
        let mut em = EnergyMonitor::new();
        em.record_energy(10.0);
        em.record_energy(10.001);
        assert!(em.is_energy_bounded(0.01));
    }
    #[test]
    fn test_legendre_momentum_free_particle() {
        let dl = DiscreteLagrangian::new(1.0, 0.0, 0.0);
        let dlt = DiscreteLegendreTransform::new(dl, 0.01);
        let q0 = [0.0; 3];
        let q1 = [0.01, 0.0, 0.0];
        let p = dlt.left_transform(q0, q1);
        assert!((p[0] - 1.0).abs() < 1e-4, "p[0]={}", p[0]);
    }
    #[test]
    fn test_nbody_momentum_conservation() {
        let nbody = NBodyVariational::new(vec![1.0, 1.0], 1.0, 0.001);
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let velocities = vec![[0.0, 0.1, 0.0], [0.0, -0.1, 0.0]];
        let p0 = nbody.total_momentum(&velocities);
        let (_, v1) = nbody.step(&positions, &velocities);
        let p1 = nbody.total_momentum(&v1);
        let dp = len3(sub3(p0, p1));
        assert!(dp < 1e-10, "momentum change={dp}");
    }
    #[test]
    fn test_nbody_energy() {
        let nbody = NBodyVariational::new(vec![1.0, 1.0], 1.0, 0.01);
        let positions = vec![[0.0; 3], [2.0, 0.0, 0.0]];
        let velocities = vec![[0.0; 3], [0.0; 3]];
        let e = nbody.total_energy(&positions, &velocities);
        assert!(e < 0.0, "energy should be negative, got {e}");
    }
    #[test]
    fn test_nbody_angular_momentum() {
        let nbody = NBodyVariational::new(vec![1.0, 1.0], 1.0, 0.01);
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let velocities = vec![[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]];
        let l = nbody.total_angular_momentum(&positions, &velocities);
        assert!((l[2] - 1.0).abs() < 1e-10, "Lz={}", l[2]);
    }
}
