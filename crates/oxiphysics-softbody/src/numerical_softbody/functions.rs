//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::ImplicitSpring;

#[inline]
pub(super) fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
#[inline]
pub(super) fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}
#[cfg(test)]
#[inline]
pub(super) fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
#[inline]
pub(super) fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n > 1e-15 {
        v3_scale(a, 1.0 / n)
    } else {
        [0.0, 0.0, 0.0]
    }
}
/// Flat vector dot product.
pub(super) fn vec_dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Flat vector norm.
pub(super) fn vec_norm(a: &[f64]) -> f64 {
    vec_dot(a, a).sqrt()
}
pub(super) fn mat3_det(m: [[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}
pub(super) fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut t = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] = m[j][i];
        }
    }
    t
}
pub(super) fn mat3_inv_scale(m: [[f64; 3]; 3], det: f64) -> [[f64; 3]; 3] {
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}
pub(super) fn mat3_mul_v3(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}
/// Solve `H x = b` using matrix-free CG.
///
/// `H` is provided as a function `apply: &[[f64;3\]] -> &[f64] -> Vec`f64`.
///
/// Returns the solution vector.
pub fn matrix_free_cg<F>(
    positions: &[[f64; 3]],
    b: &[f64],
    apply: F,
    max_iter: usize,
    tol: f64,
) -> Vec<f64>
where
    F: Fn(&[[f64; 3]], &[f64]) -> Vec<f64>,
{
    let n = b.len();
    let mut x = vec![0.0f64; n];
    let mut r = b.to_vec();
    let mut p = r.clone();
    let mut rr = vec_dot(&r, &r);
    for _ in 0..max_iter {
        let ap = apply(positions, &p);
        let pap = vec_dot(&p, &ap);
        if pap.abs() < 1e-30 {
            break;
        }
        let alpha = rr / pap;
        for i in 0..n {
            x[i] += alpha * p[i];
        }
        for i in 0..n {
            r[i] -= alpha * ap[i];
        }
        let rr_new = vec_dot(&r, &r);
        if rr_new.sqrt() < tol {
            break;
        }
        let beta = rr_new / rr.max(1e-30);
        rr = rr_new;
        for i in 0..n {
            p[i] = r[i] + beta * p[i];
        }
    }
    x
}
/// Compute total elastic potential energy for a spring-mass system.
///
/// `E = sum_springs 0.5 * k * (|x_i - x_j| - L)^2`
pub fn spring_elastic_energy(positions: &[[f64; 3]], springs: &[ImplicitSpring]) -> f64 {
    springs
        .iter()
        .map(|s| {
            let diff = v3_sub(positions[s.i], positions[s.j]);
            let len = v3_norm(diff);
            let stretch = len - s.rest_length;
            0.5 * s.stiffness * stretch * stretch
        })
        .sum()
}
/// Zero out gradient entries corresponding to static (pinned) particles.
pub fn apply_static_mask(grad: &mut [f64], is_static: &[bool]) {
    for (i, &fixed) in is_static.iter().enumerate() {
        if fixed {
            grad[3 * i] = 0.0;
            grad[3 * i + 1] = 0.0;
            grad[3 * i + 2] = 0.0;
        }
    }
}
pub(super) const _PI: f64 = PI;
#[cfg(test)]
mod tests {
    use super::*;
    use crate::numerical_softbody::BackwardEulerIntegrator;
    use crate::numerical_softbody::CsrMatrix;
    use crate::numerical_softbody::FastShapeMatchingConstraint;
    use crate::numerical_softbody::FilterLineSearch;
    use crate::numerical_softbody::ImplicitParticle;
    use crate::numerical_softbody::IpcBarrierParams;
    use crate::numerical_softbody::IpcContactSolver;
    use crate::numerical_softbody::LineSearch;
    use crate::numerical_softbody::MatrixFreeHessian;
    use crate::numerical_softbody::NewtonRaphsonSolver;
    use crate::numerical_softbody::PcgParams;
    use crate::numerical_softbody::PcgSolver;
    use crate::numerical_softbody::ProjectiveConstraint;
    use crate::numerical_softbody::ProjectiveDynamicsSolver;
    fn make_two_particle_spring(k: f64, rest: f64) -> (Vec<ImplicitParticle>, Vec<ImplicitSpring>) {
        let particles = vec![
            ImplicitParticle::new([0.0, 0.0, 0.0], 1.0),
            ImplicitParticle::new([2.0, 0.0, 0.0], 1.0),
        ];
        let springs = vec![ImplicitSpring::new(0, 1, rest, k, 0.0)];
        (particles, springs)
    }
    #[test]
    fn test_v3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = v3_add(a, b);
        assert_eq!(c, [5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_v3_dot() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert_eq!(v3_dot(a, b), 0.0);
        let c = [3.0, 4.0, 0.0];
        assert!((v3_norm(c) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_v3_cross() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = v3_cross(x, y);
        assert!((z[0]).abs() < 1e-12);
        assert!((z[1]).abs() < 1e-12);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_csr_matvec_identity() {
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];
        let vals = vec![1.0, 1.0, 1.0];
        let m = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let x = vec![1.0, 2.0, 3.0];
        let y = m.matvec(&x);
        assert!((y[0] - 1.0).abs() < 1e-12);
        assert!((y[1] - 2.0).abs() < 1e-12);
        assert!((y[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_csr_diagonal() {
        let rows = vec![0, 0, 1, 1, 2, 2];
        let cols = vec![0, 1, 0, 1, 1, 2];
        let vals = vec![4.0, -1.0, -1.0, 3.0, -1.0, 5.0];
        let m = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let d = m.diagonal();
        assert!((d[0] - 4.0).abs() < 1e-12);
        assert!((d[1] - 3.0).abs() < 1e-12);
        assert!((d[2] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_pcg_2x2() {
        let rows = vec![0, 0, 1, 1];
        let cols = vec![0, 1, 0, 1];
        let vals = vec![4.0, -1.0, -1.0, 3.0];
        let a = CsrMatrix::from_coo(2, 2, &rows, &cols, &vals);
        let b = vec![5.0, -1.0];
        let solver = PcgSolver::new(PcgParams {
            max_iter: 100,
            tolerance: 1e-10,
        });
        let result = solver.solve(&a, &b, None);
        assert!(result.converged, "PCG should converge on a 2x2 SPD system");
        assert!(
            (result.x[0] - 14.0 / 11.0).abs() < 1e-6,
            "x[0] = {}",
            result.x[0]
        );
        assert!(
            (result.x[1] - 1.0 / 11.0).abs() < 1e-6,
            "x[1] = {}",
            result.x[1]
        );
    }
    #[test]
    fn test_pcg_diagonal_3x3() {
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];
        let vals = vec![2.0, 3.0, 4.0];
        let a = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let b = vec![6.0, 9.0, 12.0];
        let solver = PcgSolver::new(PcgParams::default());
        let result = solver.solve(&a, &b, None);
        assert!(result.converged);
        assert!((result.x[0] - 3.0).abs() < 1e-8);
        assert!((result.x[1] - 3.0).abs() < 1e-8);
        assert!((result.x[2] - 3.0).abs() < 1e-8);
    }
    #[test]
    fn test_spring_energy_at_rest() {
        let (particles, springs) = make_two_particle_spring(100.0, 2.0);
        let e = springs[0].potential_energy(&particles);
        assert!(e.abs() < 1e-12, "Energy at rest should be zero, got {e}");
    }
    #[test]
    fn test_spring_force_at_rest() {
        let (particles, springs) = make_two_particle_spring(100.0, 2.0);
        let f = springs[0].force_on_i(&particles);
        assert!(v3_norm(f) < 1e-12, "Force at rest should be zero");
    }
    #[test]
    fn test_spring_force_stretched() {
        let (particles, springs) = make_two_particle_spring(100.0, 1.0);
        let f = springs[0].force_on_i(&particles);
        assert!(
            v3_norm(f) > 0.0,
            "Stretched spring should have nonzero force"
        );
    }
    #[test]
    fn test_backward_euler_gravity() {
        let mut particles = vec![ImplicitParticle::new([0.0, 10.0, 0.0], 1.0)];
        let springs: Vec<ImplicitSpring> = vec![];
        let integrator = BackwardEulerIntegrator::new(10, 1e-6);
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            integrator.step(&mut particles, &springs, gravity, dt);
        }
        assert!(
            particles[0].position[1] < 9.0,
            "Particle should fall under gravity"
        );
    }
    #[test]
    fn test_backward_euler_static_unmoved() {
        let mut particles = vec![ImplicitParticle::new_static([5.0, 5.0, 5.0])];
        let springs: Vec<ImplicitSpring> = vec![];
        let integrator = BackwardEulerIntegrator::new(5, 1e-6);
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..10 {
            integrator.step(&mut particles, &springs, gravity, 1.0 / 60.0);
        }
        let pos = particles[0].position;
        assert!((pos[0] - 5.0).abs() < 1e-12);
        assert!((pos[1] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_line_search_decreasing() {
        let ls = LineSearch::default();
        let f0 = 10.0;
        let grad = vec![1.0, 0.0];
        let alpha = ls.backtrack(f0, &grad, |a| f0 - a * 1.0);
        assert!(alpha > 0.0, "Line search should return positive step");
    }
    #[test]
    fn test_ipc_barrier_zero_outside() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        assert_eq!(p.barrier(2e-3), 0.0, "Barrier should be zero for d > d_hat");
        assert_eq!(p.barrier(-1.0), 0.0, "Barrier should be zero for d < 0");
    }
    #[test]
    fn test_ipc_barrier_positive_inside() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        let b = p.barrier(5e-4);
        assert!(
            b > 0.0,
            "Barrier should be positive for d in (0, d_hat), got {b}"
        );
    }
    #[test]
    fn test_ipc_barrier_grows_near_zero() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        let b1 = p.barrier(8e-4);
        let b2 = p.barrier(1e-4);
        assert!(b2 > b1, "Barrier should grow as d approaches 0");
    }
    #[test]
    fn test_ipc_barrier_gradient_sign() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1e6,
        };
        let g = p.barrier_gradient(5e-4);
        assert!(
            g < 0.0,
            "Barrier gradient should be negative inside active zone, got {g}"
        );
    }
    #[test]
    fn test_ipc_total_energy() {
        let p = IpcBarrierParams {
            d_hat: 1e-3,
            kappa: 1.0,
        };
        let distances = vec![5e-4, 8e-4, 2e-3];
        let total = p.total_barrier_energy(&distances);
        let expected = p.barrier(5e-4) + p.barrier(8e-4) + p.barrier(2e-3);
        assert!((total - expected).abs() < 1e-15);
    }
    #[test]
    fn test_pd_spring_local_step() {
        let c = ProjectiveConstraint::spring(0, 1, 1.0, 1000.0);
        let positions = vec![[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let proj = c.local_step(&positions);
        let dist = v3_norm(v3_sub(proj[0], proj[1]));
        assert!(
            (dist - 1.0).abs() < 1e-10,
            "Local step should give rest length, got {dist}"
        );
    }
    #[test]
    fn test_pd_anchor_local_step() {
        let target = [1.0, 2.0, 3.0];
        let c = ProjectiveConstraint::anchor(0, target, 1000.0);
        let positions = vec![[0.0, 0.0, 0.0]];
        let proj = c.local_step(&positions);
        assert_eq!(proj[0], target);
    }
    #[test]
    fn test_pd_solver_anchor() {
        let mut solver = ProjectiveDynamicsSolver::new(1, vec![1.0], 1.0 / 60.0, 20);
        solver.add_constraint(ProjectiveConstraint::anchor(0, [0.0, 0.0, 0.0], 1e6));
        let mut positions = vec![[5.0, 0.0, 0.0]];
        let mut velocities = vec![[0.0, 0.0, 0.0]];
        for _ in 0..100 {
            solver.step(&mut positions, &mut velocities, [0.0, 0.0, 0.0]);
        }
        assert!(
            v3_norm(positions[0]) < 1.0,
            "Particle should move close to anchor"
        );
    }
    #[test]
    fn test_pd_inertia_target() {
        let solver = ProjectiveDynamicsSolver::new(1, vec![1.0], 0.1, 5);
        let pos = vec![[0.0, 0.0, 0.0]];
        let vel = vec![[1.0, 0.0, 0.0]];
        let y = solver.inertia_target(&pos, &vel, [0.0, -10.0, 0.0]);
        assert!((y[0][0] - 0.1).abs() < 1e-12);
        assert!((y[0][1] - (-0.1)).abs() < 1e-12);
    }
    #[test]
    fn test_fast_shape_matching_rest() {
        let rest = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let masses = vec![1.0, 1.0, 1.0, 1.0];
        let c = FastShapeMatchingConstraint::new(rest.clone(), masses, 1.0);
        let goals = c.goal_positions(&rest);
        for (i, (g, r)) in goals.iter().zip(rest.iter()).enumerate() {
            let diff = v3_norm(v3_sub(*g, *r));
            assert!(
                diff < 1e-4,
                "Goal {i} should equal rest pos for undeformed config, diff={diff}"
            );
        }
    }
    #[test]
    fn test_center_of_mass() {
        let positions = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let com = FastShapeMatchingConstraint::center_of_mass(&positions, &masses);
        assert!(v3_norm(com) < 1e-12, "CoM should be at origin");
    }
    #[test]
    fn test_matrix_free_hessian_diagonal() {
        let masses = vec![2.0, 3.0];
        let springs: Vec<ImplicitSpring> = vec![];
        let hess = MatrixFreeHessian::new(springs, masses, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let v = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let hv = hess.apply(&positions, &v);
        assert!((hv[0] - 2.0).abs() < 1e-12);
        assert!((hv[4] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_spring_energy_levels() {
        let positions_rest = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let positions_stretch = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let springs = vec![ImplicitSpring::new(0, 1, 1.0, 100.0, 0.0)];
        let e_rest = spring_elastic_energy(&positions_rest, &springs);
        let e_stretch = spring_elastic_energy(&positions_stretch, &springs);
        assert!(e_rest < 1e-12, "Energy at rest should be zero");
        assert!(
            (e_stretch - 50.0).abs() < 1e-10,
            "E = 0.5*100*1^2 = 50, got {e_stretch}"
        );
    }
    #[test]
    fn test_ipc_plane_distances() {
        let params = IpcBarrierParams::default();
        let filter = FilterLineSearch::default();
        let solver = IpcContactSolver::new(params, filter);
        let positions = vec![[0.0, 1.0, 0.0], [0.0, -0.5, 0.0]];
        let normal = [0.0, 1.0, 0.0];
        let offset = 0.0;
        let d = solver.plane_distances(&positions, normal, offset);
        assert!((d[0] - 1.0).abs() < 1e-12);
        assert!((d[1] - (-0.5)).abs() < 1e-12);
    }
    #[test]
    fn test_static_mask() {
        let mut grad = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let is_static = vec![true, false];
        apply_static_mask(&mut grad, &is_static);
        assert_eq!(grad[0], 0.0);
        assert_eq!(grad[1], 0.0);
        assert_eq!(grad[2], 0.0);
        assert!((grad[3] - 4.0).abs() < 1e-12);
        assert!((grad[4] - 5.0).abs() < 1e-12);
        assert!((grad[5] - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_newton_raphson_spring() {
        let (mut particles, springs) = make_two_particle_spring(100_000.0, 1.0);
        let inertia: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let pos_init: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let e0 = spring_elastic_energy(&pos_init, &springs);
        let solver = NewtonRaphsonSolver::new(50, 1e-8);
        let _iters = solver.minimise(&mut particles, &springs, &inertia, 0.1);
        let pos_final: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let e1 = spring_elastic_energy(&pos_final, &springs);
        assert!(
            e1 < e0,
            "Minimiser should reduce spring energy: e0={e0:.2}, e1={e1:.2}"
        );
    }
    #[test]
    fn test_csr_duplicate_sum() {
        let rows = vec![0, 0, 1];
        let cols = vec![0, 0, 1];
        let vals = vec![1.0, 2.0, 3.0];
        let m = CsrMatrix::from_coo(2, 2, &rows, &cols, &vals);
        let d = m.diagonal();
        assert!(
            (d[0] - 3.0).abs() < 1e-12,
            "Duplicate entries should sum, got {}",
            d[0]
        );
        assert!((d[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_filter_line_search_ccd() {
        let ls = FilterLineSearch::default();
        let positions = vec![[0.0, 1.0, 0.0]];
        let direction = vec![[0.0, -1.0, 0.0]];
        let distances = vec![1.0];
        let grads = vec![-1.0];
        let alpha = ls.ccd_step_fraction(&positions, &direction, &distances, &grads);
        assert!(
            alpha > 0.0 && alpha <= ls.tau_max,
            "CCD step should be in (0, tau_max]"
        );
    }
    #[test]
    fn test_mat3_det_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let d = mat3_det(id);
        assert!((d - 1.0).abs() < 1e-12, "det(I) = 1, got {d}");
    }
    #[test]
    fn test_v3_normalize() {
        let v = [3.0, 4.0, 0.0];
        let n = v3_normalize(v);
        assert!(
            (v3_norm(n) - 1.0).abs() < 1e-12,
            "Normalized vector should have unit length"
        );
    }
    #[test]
    fn test_jacobi_preconditioner() {
        let rows = vec![0, 1, 2];
        let cols = vec![0, 1, 2];
        let vals = vec![4.0, 8.0, 2.0];
        let m = CsrMatrix::from_coo(3, 3, &rows, &cols, &vals);
        let prec = m.jacobi_preconditioner();
        assert!((prec[0] - 0.25).abs() < 1e-12);
        assert!((prec[1] - 0.125).abs() < 1e-12);
        assert!((prec[2] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_matrix_free_cg_diagonal() {
        let positions: Vec<[f64; 3]> = vec![];
        let b = vec![4.0, 8.0, 16.0];
        let x = matrix_free_cg(
            &positions,
            &b,
            |_pos, v| vec![2.0 * v[0], 4.0 * v[1], 8.0 * v[2]],
            100,
            1e-10,
        );
        assert!((x[0] - 2.0).abs() < 1e-8, "x[0] = {}", x[0]);
        assert!((x[1] - 2.0).abs() < 1e-8, "x[1] = {}", x[1]);
        assert!((x[2] - 2.0).abs() < 1e-8, "x[2] = {}", x[2]);
    }
}
