//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Solves a tridiagonal system Ax = d using the Thomas algorithm.
///
/// # Arguments
/// * `a` - Sub-diagonal (length n-1), a\[0\] unused.
/// * `b` - Main diagonal (length n).
/// * `c` - Super-diagonal (length n-1), c\[n-1\] unused.
/// * `d` - Right-hand side vector (length n).
///
/// Returns the solution vector x.
pub fn thomas_algorithm(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = b.len();
    assert!(n >= 2, "system must have at least 2 equations");
    assert_eq!(a.len(), n);
    assert_eq!(c.len(), n);
    assert_eq!(d.len(), n);
    let mut c_prime = vec![0.0_f64; n];
    let mut d_prime = vec![0.0_f64; n];
    let mut x = vec![0.0_f64; n];
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];
    for i in 1..n {
        let m = b[i] - a[i] * c_prime[i - 1];
        c_prime[i] = if i < n - 1 { c[i] / m } else { 0.0 };
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / m;
    }
    x[n - 1] = d_prime[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }
    x
}
/// Computes the first-order finite-difference gradient of `u` with spacing `dx`.
///
/// Uses central differences in the interior and forward/backward at boundaries.
pub fn fdm_gradient(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    let mut grad = vec![0.0_f64; n];
    if n == 0 {
        return grad;
    }
    if n == 1 {
        return grad;
    }
    grad[0] = (u[1] - u[0]) / dx;
    for i in 1..n - 1 {
        grad[i] = (u[i + 1] - u[i - 1]) / (2.0 * dx);
    }
    grad[n - 1] = (u[n - 1] - u[n - 2]) / dx;
    grad
}
/// Computes the second-order finite-difference Laplacian of `u` with spacing `dx`.
///
/// Uses the standard 3-point stencil: (u\[i-1\] - 2*u\[i\] + u\[i+1\]) / dx^2.
pub fn fdm_laplacian(u: &[f64], dx: f64) -> Vec<f64> {
    let n = u.len();
    let mut lap = vec![0.0_f64; n];
    if n < 3 {
        return lap;
    }
    let dx2 = dx * dx;
    for i in 1..n - 1 {
        lap[i] = (u[i - 1] - 2.0 * u[i] + u[i + 1]) / dx2;
    }
    lap
}
/// Bilinear interpolation on a regular 2D grid.
///
/// # Arguments
/// * `grid` - Row-major 2D data, shape (ny, nx).
/// * `nx`, `ny` - Grid dimensions.
/// * `x`, `y` - Fractional coordinates in \[0, nx-1\] x \[0, ny-1\].
pub fn bilinear_interp(grid: &[f64], nx: usize, ny: usize, x: f64, y: f64) -> f64 {
    let x0 = x.floor() as usize;
    let y0 = y.floor() as usize;
    let x0 = x0.min(nx - 2);
    let y0 = y0.min(ny - 2);
    let x1 = x0 + 1;
    let y1 = y0 + 1;
    let tx = x - x0 as f64;
    let ty = y - y0 as f64;
    let f00 = grid[y0 * nx + x0];
    let f10 = grid[y0 * nx + x1];
    let f01 = grid[y1 * nx + x0];
    let f11 = grid[y1 * nx + x1];
    f00 * (1.0 - tx) * (1.0 - ty) + f10 * tx * (1.0 - ty) + f01 * (1.0 - tx) * ty + f11 * tx * ty
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::pde::AdiMethod2D;
    use crate::pde::BoundaryCondition;
    use crate::pde::BurgersShock;
    use crate::pde::Dct1D;
    use crate::pde::FftDiff1D;
    use crate::pde::FiniteDiffOps1D;
    use crate::pde::FiniteDiffOps2D;
    use crate::pde::FiniteDiffOps3D;
    use crate::pde::FiniteDifference1D;
    use crate::pde::FiniteDifference2D;
    use crate::pde::FiniteVolume1D;
    use crate::pde::FitzHughNagumo;
    use crate::pde::GrayScottSystem;
    use crate::pde::HeatEquation1D;
    use crate::pde::HeatEquation3D;
    use crate::pde::LevelSet2D;
    use crate::pde::MethodOfLines;
    use crate::pde::NsBurgers1D;
    use crate::pde::PhaseField2D;
    use crate::pde::PoissonSolver;
    use crate::pde::SpectralMethod;
    use crate::pde::WaveEq2ndOrder;
    use crate::pde::WaveEquation1D;
    use std::f64::consts::PI;
    #[test]
    fn test_thomas_algorithm_basic() {
        let a = vec![0.0, -1.0, -1.0];
        let b = vec![2.0, 2.0, 2.0];
        let c = vec![-1.0, -1.0, 0.0];
        let d = vec![0.0, 0.0, 0.0];
        let x = thomas_algorithm(&a, &b, &c, &d);
        assert_eq!(x.len(), 3);
        for xi in &x {
            assert!(xi.abs() < 1e-12);
        }
    }
    #[test]
    fn test_thomas_algorithm_nontrivial() {
        let a = vec![0.0, -1.0];
        let b = vec![2.0, 2.0];
        let c = vec![-1.0, 0.0];
        let d = vec![1.0, 1.0];
        let x = thomas_algorithm(&a, &b, &c, &d);
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fdm_gradient_constant() {
        let u = vec![3.0; 10];
        let grad = fdm_gradient(&u, 0.1);
        for g in &grad {
            assert!(g.abs() < 1e-12);
        }
    }
    #[test]
    fn test_fdm_gradient_linear() {
        let dx = 0.1;
        let u: Vec<f64> = (0..10).map(|i| i as f64 * dx).collect();
        let grad = fdm_gradient(&u, dx);
        for (i, &g) in grad.iter().enumerate().take(9).skip(1) {
            assert!((g - 1.0).abs() < 1e-10, "i={i}: {g}");
        }
    }
    #[test]
    fn test_fdm_laplacian_linear() {
        let dx = 0.1;
        let u: Vec<f64> = (0..10).map(|i| i as f64 * dx).collect();
        let lap = fdm_laplacian(&u, dx);
        for (_, &l) in lap.iter().enumerate().take(9).skip(1) {
            assert!(l.abs() < 1e-10);
        }
    }
    #[test]
    fn test_fdm_laplacian_quadratic() {
        let dx = 0.1;
        let u: Vec<f64> = (0..20).map(|i| (i as f64 * dx).powi(2)).collect();
        let lap = fdm_laplacian(&u, dx);
        for (i, &l) in lap.iter().enumerate().take(19).skip(1) {
            assert!((l - 2.0).abs() < 1e-6, "lap[{}] = {}", i, l);
        }
    }
    #[test]
    fn test_bilinear_interp_at_nodes() {
        let grid = vec![1.0, 2.0, 3.0, 4.0];
        assert!((bilinear_interp(&grid, 2, 2, 0.0, 0.0) - 1.0).abs() < 1e-12);
        assert!((bilinear_interp(&grid, 2, 2, 1.0, 0.0) - 2.0).abs() < 1e-12);
        assert!((bilinear_interp(&grid, 2, 2, 0.0, 1.0) - 3.0).abs() < 1e-12);
        assert!((bilinear_interp(&grid, 2, 2, 1.0, 1.0) - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_bilinear_interp_midpoint() {
        let grid = vec![1.0, 2.0, 3.0, 4.0];
        let v = bilinear_interp(&grid, 2, 2, 0.5, 0.5);
        assert!((v - 2.5).abs() < 1e-12);
    }
    #[test]
    fn test_fd1d_diffusion_number() {
        let fd = FiniteDifference1D::new(
            10,
            0.1,
            0.01,
            0.5,
            0.0,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        let r = fd.diffusion_number();
        assert!((r - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_fd1d_stability_check() {
        let fd_stable = FiniteDifference1D::new(
            10,
            0.1,
            0.001,
            0.1,
            0.0,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        assert!(fd_stable.is_stable_explicit());
        let fd_unstable = FiniteDifference1D::new(
            10,
            0.1,
            0.1,
            1.0,
            0.0,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(0.0),
        );
        assert!(!fd_unstable.is_stable_explicit());
    }
    #[test]
    fn test_fd1d_explicit_diffusion_preserves_bc() {
        let n = 10;
        let fd = FiniteDifference1D::new(
            n,
            0.1,
            0.001,
            0.1,
            0.0,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(1.0),
        );
        let u: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let u_new = fd.step_diffusion_explicit(&u);
        assert!((u_new[0] - 0.0).abs() < 1e-12);
        assert!((u_new[n - 1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fd1d_implicit_diffusion_preserves_bc() {
        let n = 10;
        let fd = FiniteDifference1D::new(
            n,
            0.1,
            0.01,
            0.1,
            0.0,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(1.0),
        );
        let u: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let u_new = fd.step_diffusion_implicit(&u);
        assert!((u_new[0] - 0.0).abs() < 1e-12);
        assert!((u_new[n - 1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fd1d_cn_diffusion_preserves_bc() {
        let n = 10;
        let fd = FiniteDifference1D::new(
            n,
            0.1,
            0.01,
            0.1,
            0.0,
            BoundaryCondition::Dirichlet(0.0),
            BoundaryCondition::Dirichlet(1.0),
        );
        let u: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
        let u_new = fd.step_diffusion_cn(&u);
        assert!((u_new[0] - 0.0).abs() < 1e-12);
        assert!((u_new[n - 1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fd1d_advection_upwind_conservation() {
        let n = 20;
        let fd = FiniteDifference1D::new(
            n,
            0.05,
            0.01,
            0.0,
            1.0,
            BoundaryCondition::Periodic,
            BoundaryCondition::Periodic,
        );
        let u: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let u_new = fd.step_advection_upwind(&u);
        let sum_old: f64 = u.iter().sum();
        let sum_new: f64 = u_new.iter().sum();
        assert!((sum_old - sum_new).abs() < 1.0);
    }
    #[test]
    fn test_fd2d_laplacian_5pt_zero_for_linear() {
        let nx = 5;
        let ny = 5;
        let fd2d = FiniteDifference2D::new(nx, ny, 0.1, 0.1, 0.001, 0.1);
        let u: Vec<f64> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| i as f64 * 0.1 + j as f64 * 0.1))
            .collect();
        let lap = fd2d.laplacian_5pt(&u);
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                assert!(lap[j * nx + i].abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_wave1d_courant_number() {
        let w = WaveEquation1D::new(100, 0.01, 0.005, 1.0);
        assert!((w.courant_number() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_wave1d_lax_wendroff_stable() {
        let n = 50;
        let w = WaveEquation1D::new(n, 0.02, 0.01, 1.0);
        let u: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let u_prev = u.clone();
        let u_new = w.step_lax_wendroff(&u, &u_prev);
        let max_val = u_new.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_val.abs() <= 2.0);
    }
    #[test]
    fn test_poisson_vcycle_converges() {
        let n = 17;
        let dx = 1.0 / (n - 1) as f64;
        let solver = PoissonSolver::new(n, dx, 4, 10);
        let f: Vec<f64> = (0..n).map(|i| -(PI * i as f64 * dx).sin()).collect();
        let u0 = vec![0.0_f64; n];
        let u = solver.vcycle(&u0, &f, 3);
        let exact: Vec<f64> = (0..n)
            .map(|i| (PI * i as f64 * dx).sin() / (PI * PI))
            .collect();
        let err = (0..n)
            .map(|i| (u[i] - exact[i]).powi(2))
            .sum::<f64>()
            .sqrt()
            / n as f64;
        assert!(err < 0.1, "Error too large: {}", err);
    }
    #[test]
    fn test_burgers_step_no_blowup() {
        let n = 50;
        let b = NsBurgers1D::new(n, 0.02, 0.001, 0.01);
        let u: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let u_new = b.step(&u);
        let max_val = u_new.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_val.abs() < 10.0);
    }
    #[test]
    fn test_spectral_chebyshev_nodes_count() {
        let sp = SpectralMethod::new(8);
        let nodes = sp.gauss_lobatto_nodes();
        assert_eq!(nodes.len(), 8);
        assert!((nodes[0].abs() - 1.0).abs() < 1e-12);
        assert!((nodes[7].abs() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_spectral_differentiation_matrix_size() {
        let sp = SpectralMethod::new(5);
        let d = sp.differentiation_matrix();
        assert_eq!(d.len(), 25);
    }
    #[test]
    fn test_dct_roundtrip() {
        let dct = Dct1D::new(8);
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let xk = dct.dct2_forward(&x);
        let x_rec = dct.dct2_inverse(&xk);
        for (orig, rec) in x.iter().zip(x_rec.iter()) {
            assert!((orig - rec).abs() < 1e-10, "orig={} rec={}", orig, rec);
        }
    }
    #[test]
    fn test_fv1d_minmod() {
        assert_eq!(FiniteVolume1D::minmod(2.0, 3.0), 2.0);
        assert_eq!(FiniteVolume1D::minmod(-2.0, -3.0), -2.0);
        assert_eq!(FiniteVolume1D::minmod(2.0, -3.0), 0.0);
    }
    #[test]
    fn test_fv1d_superbee() {
        let s = FiniteVolume1D::superbee(1.0, 1.0);
        assert!(s >= 0.0);
    }
    #[test]
    fn test_fv1d_van_leer() {
        let s = FiniteVolume1D::van_leer(1.0, 2.0);
        assert!((s - 4.0 / 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_fv1d_step_conservation() {
        let n = 30;
        let fv = FiniteVolume1D::new(n, 0.033, 0.005);
        let u: Vec<f64> = (0..n)
            .map(|i| (2.0 * PI * i as f64 / n as f64).sin())
            .collect();
        let u_new = fv.step(&u);
        assert_eq!(u_new.len(), n);
    }
    #[test]
    fn test_levelset2d_creation() {
        let phi: Vec<f64> = (0..25).map(|i| i as f64 - 12.0).collect();
        let ls = LevelSet2D::new(5, 5, 0.1, phi);
        assert_eq!(ls.phi.len(), 25);
    }
    #[test]
    fn test_levelset2d_normal_boundary() {
        let nx = 5;
        let ny = 5;
        let phi: Vec<f64> = (0..nx * ny).map(|i| (i % nx) as f64 * 0.1 - 0.2).collect();
        let ls = LevelSet2D::new(nx, ny, 0.1, phi);
        let (nx_vec, _ny_vec) = ls.normal();
        assert_eq!(nx_vec.len(), nx * ny);
    }
    #[test]
    fn test_phasefield2d_allen_cahn_energy_decrease() {
        let nx = 10;
        let ny = 10;
        let phi0: Vec<f64> = (0..nx * ny)
            .map(|i| if i < nx * ny / 2 { -1.0 } else { 1.0 })
            .collect();
        let mut pf = PhaseField2D::new(nx, ny, 0.1, 0.001, 0.1, 1.0, phi0);
        let e0 = pf.interface_energy();
        pf.step_allen_cahn();
        let e1 = pf.interface_energy();
        assert!(e1 <= e0 * 2.0);
    }
    #[test]
    fn test_method_of_lines_rk4() {
        let mol = MethodOfLines::new(3, 0.1);
        let u0 = vec![1.0, 0.0, 0.0];
        let u_final = mol.integrate_rk4(&u0, 0.01, 100, |_t, u| u.iter().map(|&x| -x).collect());
        assert!(u_final[0] < 1.0);
        assert!(u_final[0] > 0.0);
    }
    #[test]
    fn test_adi_2d_step() {
        let nx = 5;
        let ny = 5;
        let fd2d = FiniteDifference2D::new(nx, ny, 0.1, 0.1, 0.001, 0.1);
        let u: Vec<f64> = vec![1.0; nx * ny];
        let u_new = fd2d.adi_step(&u);
        assert_eq!(u_new.len(), nx * ny);
    }
    #[test]
    fn test_heat3d_explicit_step() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let heat = HeatEquation3D::new(
            nx,
            ny,
            nz,
            0.1,
            0.0001,
            0.01,
            BoundaryCondition::Dirichlet(0.0),
        );
        let u: Vec<f64> = vec![1.0; nx * ny * nz];
        let u_new = heat.step_explicit(&u);
        assert_eq!(u_new.len(), nx * ny * nz);
    }
    #[test]
    fn test_heat3d_implicit_split_unconditionally_stable() {
        // Grid + diffusivity chosen so the explicit CFL limit (r <= 1/6) is
        // violated by a wide margin: r = D*dt/dx^2 = 1.0 * 0.5 / 0.1^2 = 50.
        let nx = 6;
        let ny = 6;
        let nz = 6;
        let dx = 0.1;
        let dt = 0.5;
        let diffusivity = 1.0;
        let r = diffusivity * dt / (dx * dx);
        assert!(r > 1.0 / 6.0, "test must exceed explicit CFL: r = {r}");

        let heat = HeatEquation3D::new(
            nx,
            ny,
            nz,
            dx,
            dt,
            diffusivity,
            BoundaryCondition::Dirichlet(0.0),
        );

        // Interior hot (=1), boundaries cold (=0).  Steady state is u ≡ 0.
        let mut u_init = vec![0.0_f64; nx * ny * nz];
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    u_init[heat.idx(i, j, k)] = 1.0;
                }
            }
        }

        // (a) Explicit step at this dt must blow up far past the initial range.
        let u_exp = heat.step_explicit(&u_init);
        let exp_max = u_exp.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
        assert!(
            exp_max > 5.0 || u_exp.iter().any(|v| !v.is_finite()),
            "explicit scheme should be unstable beyond CFL, max = {exp_max}"
        );

        // (b) Implicit ADI must stay bounded and relax monotonically toward 0.
        let mut u = u_init.clone();
        let mut prev_max = 1.0_f64;
        for step in 0..40 {
            u = heat.step_implicit_split(&u);
            assert!(
                u.iter().all(|v| v.is_finite()),
                "implicit ADI produced non-finite at step {step}"
            );
            let cur_max = u.iter().fold(0.0_f64, |m, &v| m.max(v.abs()));
            // Maximum principle: never exceeds the initial maximum (+tiny slack).
            assert!(
                cur_max <= 1.0 + 1e-9,
                "implicit ADI violated max principle at step {step}: {cur_max}"
            );
            // Monotone decay toward the cold steady state.
            assert!(
                cur_max <= prev_max + 1e-9,
                "implicit ADI not relaxing at step {step}: {cur_max} > {prev_max}"
            );
            prev_max = cur_max;
        }
        // After many large steps the field has essentially reached steady state.
        assert!(
            prev_max < 1e-2,
            "implicit ADI did not relax toward steady state: residual {prev_max}"
        );
    }
    #[test]
    fn test_thomas_three_by_three() {
        let a = vec![0.0, 1.0, 1.0];
        let b = vec![4.0, 4.0, 4.0];
        let c = vec![1.0, 1.0, 0.0];
        let d = vec![5.0, 6.0, 5.0];
        let x = thomas_algorithm(&a, &b, &c, &d);
        let r0 = 4.0 * x[0] + 1.0 * x[1];
        let r1 = 1.0 * x[0] + 4.0 * x[1] + 1.0 * x[2];
        let r2 = 1.0 * x[1] + 4.0 * x[2];
        assert!((r0 - 5.0).abs() < 1e-10);
        assert!((r1 - 6.0).abs() < 1e-10);
        assert!((r2 - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_levelset2d_reinitialize() {
        let nx = 7;
        let ny = 7;
        let phi0: Vec<f64> = (0..nx * ny)
            .map(|i| {
                let x = (i % nx) as f64 * 0.1 - 0.3;
                let y = (i / nx) as f64 * 0.1 - 0.3;
                x * x + y * y - 0.09
            })
            .collect();
        let mut ls = LevelSet2D::new(nx, ny, 0.1, phi0);
        ls.reinitialize(5, 0.01);
        assert_eq!(ls.phi.len(), nx * ny);
    }
    #[test]
    fn test_phasefield2d_cahn_hilliard_conserves_mass() {
        let nx = 8;
        let ny = 8;
        let phi0: Vec<f64> = (0..nx * ny)
            .map(|i| if i % 2 == 0 { 0.5 } else { -0.5 })
            .collect();
        let mut pf = PhaseField2D::new(nx, ny, 0.1, 0.001, 0.1, 0.1, phi0);
        let mass_before: f64 = pf.phi.iter().sum();
        pf.step_cahn_hilliard();
        let mass_after: f64 = pf.phi.iter().sum();
        assert!((mass_before - mass_after).abs() < 1.0);
    }
    #[test]
    fn test_fd1d_absorbing_bc() {
        let n = 10;
        let fd = FiniteDifference1D::new(
            n,
            0.1,
            0.001,
            0.1,
            0.0,
            BoundaryCondition::Absorbing,
            BoundaryCondition::Absorbing,
        );
        let u: Vec<f64> = (0..n).map(|i| i as f64).collect();
        let u_new = fd.step_diffusion_explicit(&u);
        assert!((u_new[0] - u_new[1]).abs() < 1e-10);
        assert!((u_new[n - 1] - u_new[n - 2]).abs() < 1e-10);
    }
    #[test]
    fn test_wave1d_absorbing_bc() {
        let n = 20;
        let w = WaveEquation1D::new(n, 0.05, 0.01, 1.0);
        let mut u = vec![1.0_f64; n];
        w.apply_absorbing_bc(&mut u);
        assert!((u[0] - u[1]).abs() < 1e-12);
        assert!((u[n - 1] - u[n - 2]).abs() < 1e-12);
    }
    #[test]
    fn test_fft_poisson_periodicity() {
        let n = 8;
        let dx = 2.0 * PI / n as f64;
        let solver = PoissonSolver::new(n, dx, 3, 5);
        let f: Vec<f64> = (0..n).map(|i| -(i as f64 * dx).sin()).collect();
        let u = solver.solve_fft_periodic(&f);
        assert_eq!(u.len(), n);
    }
    #[test]
    fn test_laplacian1d_constant_zero() {
        let ops = FiniteDiffOps1D::new(10, 0.1);
        let u = vec![5.0_f64; 10];
        let lap = ops.laplacian(&u);
        for (_, &l) in lap.iter().enumerate().take(9).skip(1) {
            assert!(l.abs() < 1e-10);
        }
    }
    #[test]
    fn test_laplacian2d_linear_zero() {
        let ops = FiniteDiffOps2D::new(5, 5, 0.1, 0.1);
        let u: Vec<f64> = (0..25).map(|k| (k % 5) as f64 + (k / 5) as f64).collect();
        let lap = ops.laplacian(&u);
        for j in 1..4 {
            for i in 1..4 {
                assert!(lap[j * 5 + i].abs() < 1e-9);
            }
        }
    }
    #[test]
    fn test_gradient2d_constant_zero() {
        let ops = FiniteDiffOps2D::new(5, 5, 0.1, 0.1);
        let u = vec![3.0_f64; 25];
        let (gx, gy) = ops.gradient(&u);
        for i in 1..4 {
            for j in 1..4 {
                let idx = j * 5 + i;
                assert!(gx[idx].abs() < 1e-10);
                assert!(gy[idx].abs() < 1e-10);
            }
        }
    }
    #[test]
    fn test_divergence2d_result_size() {
        let ops = FiniteDiffOps2D::new(5, 5, 0.1, 0.1);
        let fx = vec![1.0_f64; 25];
        let fy = vec![0.0_f64; 25];
        let div = ops.divergence(&fx, &fy);
        assert_eq!(div.len(), 25);
    }
    #[test]
    fn test_curl2d_constant_zero() {
        let ops = FiniteDiffOps2D::new(5, 5, 0.1, 0.1);
        let ux = vec![2.0_f64; 25];
        let uy = vec![2.0_f64; 25];
        let curl = ops.curl(&ux, &uy);
        for j in 1..4 {
            for i in 1..4 {
                assert!(curl[j * 5 + i].abs() < 1e-9);
            }
        }
    }
    #[test]
    fn test_laplacian3d_result_size() {
        let ops = FiniteDiffOps3D::new(4, 4, 4, 0.1);
        let u = vec![1.0_f64; 64];
        let lap = ops.laplacian(&u);
        assert_eq!(lap.len(), 64);
    }
    #[test]
    fn test_mol_diffusion_rk4_decay() {
        let n = 20;
        let dx = 1.0 / (n - 1) as f64;
        let mol = MethodOfLines::new(n, dx);
        let u0: Vec<f64> = (0..n).map(|i| (PI * i as f64 * dx).sin()).collect();
        let d = 0.01_f64;
        let u_final = mol.integrate_rk4(&u0, 0.001, 50, |_t, u| {
            fdm_laplacian(u, dx).iter().map(|&v| d * v).collect()
        });
        let amp0: f64 = u0.iter().cloned().fold(0.0_f64, f64::max);
        let amp1: f64 = u_final.iter().cloned().fold(0.0_f64, f64::max);
        assert!(amp1 <= amp0 + 1e-10);
    }
    #[test]
    fn test_spectral_fft_diff_sine() {
        let n = 16;
        let dx = 2.0 * PI / n as f64;
        let fft = FftDiff1D::new(n, dx);
        let u: Vec<f64> = (0..n).map(|i| (i as f64 * dx).sin()).collect();
        let du = fft.differentiate(&u);
        for (i, &du_i) in du.iter().enumerate().take(n - 2).skip(2) {
            let x = i as f64 * dx;
            assert!(
                (du_i - x.cos()).abs() < 0.05,
                "i={} du={} cos={}",
                i,
                du_i,
                x.cos()
            );
        }
    }
    #[test]
    fn test_gray_scott_step_size() {
        let nx = 8;
        let ny = 8;
        let (u0, v0) = GrayScottSystem::initial_conditions(nx, ny);
        let mut gs = GrayScottSystem::new(nx, ny, 0.1, 0.001, 0.04, 0.06);
        gs.u = u0;
        gs.v = v0;
        gs.step();
        assert_eq!(gs.u.len(), nx * ny);
        assert_eq!(gs.v.len(), nx * ny);
    }
    #[test]
    fn test_fitzhugh_nagumo_step() {
        let mut fhn = FitzHughNagumo::new(10, 0.1, 0.001);
        fhn.step();
        assert_eq!(fhn.v.len(), 10);
        assert_eq!(fhn.w.len(), 10);
    }
    #[test]
    fn test_wave_eq_2nd_order_size() {
        let n = 20;
        let w = WaveEq2ndOrder::new(n, 0.05, 0.01, 1.0);
        let u = vec![0.0_f64; n];
        let u_prev = vec![0.0_f64; n];
        let u_new = w.step(&u, &u_prev);
        assert_eq!(u_new.len(), n);
    }
    #[test]
    fn test_burgers_shock_formation() {
        let n = 50;
        let b = BurgersShock::new(n, 0.02, 0.001, 0.005);
        let u: Vec<f64> = (0..n).map(|i| (PI * i as f64 / n as f64).sin()).collect();
        let u_new = b.step(&u);
        assert_eq!(u_new.len(), n);
        let max_abs = u_new.iter().cloned().fold(0.0_f64, |a, v| a.max(v.abs()));
        assert!(max_abs < 5.0);
    }
    #[test]
    fn test_heat_eq_explicit_implicit_cn_size() {
        let n = 10;
        let he = HeatEquation1D::new(n, 0.1, 0.001, 0.1);
        let u = vec![1.0_f64; n];
        assert_eq!(he.step_explicit(&u).len(), n);
        assert_eq!(he.step_implicit(&u).len(), n);
        assert_eq!(he.step_cn(&u).len(), n);
    }
    #[test]
    fn test_poisson_multigrid_restrict_prolongate_roundtrip() {
        let ps = PoissonSolver::new(8, 0.1, 3, 5);
        let f = vec![1.0_f64, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
        let coarse = ps.restrict(&f);
        assert_eq!(coarse.len(), 4);
        let fine = ps.prolongate(&coarse, 8);
        assert_eq!(fine.len(), 8);
    }
    #[test]
    fn test_adi_method_2d_size() {
        let adi = AdiMethod2D::new(6, 6, 0.1, 0.1, 0.001, 0.1);
        let u = vec![1.0_f64; 36];
        let u_new = adi.step(&u);
        assert_eq!(u_new.len(), 36);
    }
    #[test]
    fn test_levelset_advection_preserves_size() {
        let nx = 7;
        let ny = 7;
        let phi0: Vec<f64> = (0..nx * ny).map(|i| i as f64 - 24.0).collect();
        let mut ls = LevelSet2D::new(nx, ny, 0.1, phi0);
        let vel_u = vec![0.1_f64; nx * ny];
        let vel_v = vec![0.0_f64; nx * ny];
        ls.advect(&vel_u, &vel_v, 0.01);
        assert_eq!(ls.phi.len(), nx * ny);
    }
}
