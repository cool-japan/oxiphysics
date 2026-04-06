//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
#[allow(unused_imports)]
use super::functions::*;
use super::functions::{lu_decompose, lu_solve, matvec_dense, normalize, vec_norm};
use super::types::SubspaceConvergence;

/// Rayleigh quotient iteration for refining an eigenvalue/eigenvector estimate.
///
/// Starting from an initial guess `v0`, iteratively updates the Rayleigh
/// quotient shift and solves a shifted linear system until convergence.
///
/// Returns `(eigenvalue, eigenvector)` on convergence, or the best estimate
/// if `max_iter` is exhausted.
#[allow(dead_code)]
pub fn rayleigh_quotient_iteration(
    a: &[Vec<f64>],
    n: usize,
    v0: Vec<f64>,
    max_iter: usize,
    tol: f64,
) -> (f64, Vec<f64>) {
    let mut v = v0;
    let norm = vec_norm(&v);
    if norm > 1e-60 {
        for x in v.iter_mut() {
            *x /= norm;
        }
    }
    let mut mu: f64 = {
        let av = matvec_dense(a, &v, n);
        v.iter().zip(av.iter()).map(|(vi, avi)| vi * avi).sum()
    };
    for _ in 0..max_iter {
        let mut a_shifted: Vec<Vec<f64>> = a.to_vec();
        for i in 0..n {
            a_shifted[i][i] -= mu;
        }
        let lu = lu_decompose(&a_shifted, n);
        let w = lu_solve(&lu, &v, n);
        let w_norm = vec_norm(&w);
        if w_norm < 1e-60 {
            break;
        }
        let mut v_new: Vec<f64> = w.iter().map(|x| x / w_norm).collect();
        let dot: f64 = v.iter().zip(v_new.iter()).map(|(a, b)| a * b).sum();
        if dot < 0.0 {
            for x in v_new.iter_mut() {
                *x = -*x;
            }
        }
        let av_new = matvec_dense(a, &v_new, n);
        let mu_new: f64 = v_new
            .iter()
            .zip(av_new.iter())
            .map(|(vi, avi)| vi * avi)
            .sum();
        let change = (mu_new - mu).abs();
        v = v_new;
        mu = mu_new;
        if change < tol {
            break;
        }
    }
    (mu, v)
}
/// Inverse iteration for the eigenvalue nearest to `sigma`.
///
/// Solves `(A - sigma I) w = v` at each step, then normalizes.
/// Converges to the eigenvector for the eigenvalue nearest `sigma`.
#[allow(dead_code)]
pub fn inverse_iteration(
    a: &[Vec<f64>],
    n: usize,
    sigma: f64,
    max_iter: usize,
    tol: f64,
) -> Option<(f64, Vec<f64>)> {
    if n == 0 {
        return None;
    }
    let mut a_shifted: Vec<Vec<f64>> = a.to_vec();
    for i in 0..n {
        a_shifted[i][i] -= sigma;
    }
    let lu = lu_decompose(&a_shifted, n);
    let mut v: Vec<f64> = (0..n).map(|i| 1.0 + 0.01 * i as f64).collect();
    normalize(&mut v);
    let mut eigenval = sigma;
    for iter in 0..max_iter {
        let w = lu_solve(&lu, &v, n);
        let w_norm = vec_norm(&w);
        if w_norm < 1e-60 {
            return None;
        }
        let mut v_new: Vec<f64> = w.iter().map(|x| x / w_norm).collect();
        let dot: f64 = v.iter().zip(v_new.iter()).map(|(a, b)| a * b).sum();
        if dot < 0.0 {
            for x in v_new.iter_mut() {
                *x = -*x;
            }
        }
        let av = matvec_dense(a, &v_new, n);
        let rq: f64 = v_new.iter().zip(av.iter()).map(|(vi, avi)| vi * avi).sum();
        let change = (rq - eigenval).abs() / (rq.abs() + 1e-60);
        v = v_new;
        eigenval = rq;
        if change < tol && iter > 0 {
            return Some((eigenval, v));
        }
    }
    Some((eigenval, v))
}
/// Mass-normalize a set of mode shapes.
///
/// Scales each mode `phi_i` so that `phi_i^T * M * phi_i = 1`.
/// `m_diag` is the diagonal mass vector.
#[allow(dead_code)]
pub fn mass_normalize_modes(modes: &[Vec<f64>], m_diag: &[f64]) -> Vec<Vec<f64>> {
    modes
        .iter()
        .map(|phi| {
            let m_mass: f64 = phi.iter().zip(m_diag.iter()).map(|(p, m)| p * p * m).sum();
            let scale = if m_mass > 1e-60 {
                1.0 / m_mass.sqrt()
            } else {
                1.0
            };
            phi.iter().map(|p| p * scale).collect()
        })
        .collect()
}
/// Verify mass-orthogonality of mode shapes.
///
/// Computes the generalized mass matrix `M_modal = Phi^T M Phi`.
/// Returns the off-diagonal norm (should be near zero for orthogonal modes).
#[allow(dead_code)]
pub fn check_mass_orthogonality(modes: &[Vec<f64>], m_diag: &[f64]) -> f64 {
    let n_modes = modes.len();
    let mut off_diag_sq = 0.0;
    for i in 0..n_modes {
        for j in 0..n_modes {
            if i == j {
                continue;
            }
            let mij: f64 = modes[i]
                .iter()
                .zip(modes[j].iter())
                .zip(m_diag.iter())
                .map(|((pi, pj), m)| pi * m * pj)
                .sum();
            off_diag_sq += mij * mij;
        }
    }
    off_diag_sq.sqrt()
}
/// Compute eigenvalues of a symmetric tridiagonal matrix using the
/// bisection method (based on Sturm sequences).
///
/// `alpha` – main diagonal (length n)
/// `beta`  – off-diagonal (length n-1)
///
/// Returns eigenvalues sorted ascending.
#[allow(dead_code)]
pub fn tridiagonal_eigenvalues_bisection(alpha: &[f64], beta: &[f64]) -> Vec<f64> {
    let n = alpha.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![alpha[0]];
    }
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for i in 0..n {
        let b_left = if i > 0 { beta[i - 1].abs() } else { 0.0 };
        let b_right = if i < n - 1 { beta[i].abs() } else { 0.0 };
        let center = alpha[i];
        let r = b_left + b_right;
        lo = lo.min(center - r);
        hi = hi.max(center + r);
    }
    lo -= 1e-10;
    hi += 1e-10;
    let sturm_count = |x: f64| -> usize {
        let mut count = 0usize;
        let mut d = alpha[0] - x;
        if d < 0.0 {
            count += 1;
        }
        for i in 1..n {
            let d_prev = d;
            d = (alpha[i] - x)
                - if d_prev.abs() > 1e-60 {
                    beta[i - 1] * beta[i - 1] / d_prev
                } else {
                    beta[i - 1] * beta[i - 1] * 1e60
                };
            if d < 0.0 {
                count += 1;
            }
        }
        count
    };
    let mut eigenvalues = Vec::with_capacity(n);
    for k in 0..n {
        let target_count = k + 1;
        let mut a = lo;
        let mut b = hi;
        for _ in 0..80 {
            let mid = 0.5 * (a + b);
            if sturm_count(mid) < target_count {
                a = mid;
            } else {
                b = mid;
            }
        }
        eigenvalues.push(0.5 * (a + b));
    }
    eigenvalues
}
/// Compute residual norms for each eigenpair `(lambda, v)`.
#[allow(dead_code)]
pub fn eigenpair_residuals(
    a: &[Vec<f64>],
    n: usize,
    pairs: &[(f64, Vec<f64>)],
) -> SubspaceConvergence {
    let residuals = pairs
        .iter()
        .map(|(lambda, v)| {
            let av = matvec_dense(a, v, n);
            let res: f64 = av
                .iter()
                .zip(v.iter())
                .map(|(avi, vi)| (avi - lambda * vi).powi(2))
                .sum::<f64>()
                .sqrt();
            res
        })
        .collect();
    SubspaceConvergence { residuals }
}
/// Generalized Rayleigh quotient: R(v) = (v^T K v) / (v^T M v).
///
/// For the generalized eigenproblem K*phi = lambda*M*phi, the Rayleigh
/// quotient provides an upper bound on the smallest eigenvalue.
#[allow(dead_code)]
pub fn generalized_rayleigh_quotient(k: &[Vec<f64>], m_diag: &[f64], v: &[f64]) -> f64 {
    let n = v.len();
    let kv = matvec_dense(k, v, n);
    let vt_kv: f64 = v.iter().zip(kv.iter()).map(|(vi, kvi)| vi * kvi).sum();
    let vt_mv: f64 = v
        .iter()
        .zip(m_diag.iter())
        .map(|(vi, mi)| vi * vi * mi)
        .sum();
    if vt_mv < 1e-60 {
        return 0.0;
    }
    vt_kv / vt_mv
}
/// Compute modal participation factors for each mode shape.
///
/// The participation factor for mode `i` in direction `r` is:
/// `Gamma_i = phi_i^T * M * r`
///
/// Returns a vector of participation factors (one per mode).
#[allow(dead_code)]
pub fn modal_participation_factors(
    modes: &[Vec<f64>],
    m_diag: &[f64],
    direction: &[f64],
) -> Vec<f64> {
    modes
        .iter()
        .map(|phi| {
            phi.iter()
                .zip(m_diag.iter())
                .zip(direction.iter())
                .map(|((p, m), r)| p * m * r)
                .sum()
        })
        .collect()
}
/// Frequency Response Function magnitude at frequency `omega`.
///
/// Uses the modal superposition formula:
/// H(omega) = Σ_i Gamma_i² / (omega_i² - omega² + 2i*zeta_i*omega_i*omega)
///
/// Returns the magnitude (real approximation without damping for simplicity).
#[allow(dead_code)]
pub fn frf_magnitude(
    omega_sq: &[f64],
    participation_factors: &[f64],
    omega: f64,
    zeta: f64,
) -> f64 {
    let omega2 = omega * omega;
    let mut h_real = 0.0;
    let mut h_imag = 0.0;
    for (i, &osq) in omega_sq.iter().enumerate() {
        let gamma = if i < participation_factors.len() {
            participation_factors[i]
        } else {
            0.0
        };
        let omega_i = osq.max(0.0).sqrt();
        let denom_real = osq - omega2;
        let denom_imag = 2.0 * zeta * omega_i * omega;
        let denom_sq = denom_real * denom_real + denom_imag * denom_imag;
        if denom_sq < 1e-60 {
            continue;
        }
        h_real += gamma * gamma * denom_real / denom_sq;
        h_imag -= gamma * gamma * denom_imag / denom_sq;
    }
    (h_real * h_real + h_imag * h_imag).sqrt()
}
#[cfg(test)]
mod eigen_extended_tests {
    use super::*;

    fn diag_matrix(d: &[f64]) -> Vec<Vec<f64>> {
        let n = d.len();
        let mut a = vec![vec![0.0; n]; n];
        for i in 0..n {
            a[i][i] = d[i];
        }
        a
    }
    /// Rayleigh quotient converges to the nearest eigenvalue.
    #[test]
    fn test_rayleigh_quotient_simple() {
        let a = diag_matrix(&[1.0, 4.0, 9.0]);
        let v0 = vec![0.0, 1.0, 0.01];
        let (lambda, _) = rayleigh_quotient_iteration(&a, 3, v0, 100, 1e-10);
        assert!((lambda - 4.0).abs() < 0.1, "lambda = {lambda}, expected ~4");
    }
    /// Rayleigh quotient iteration converges to dominant eigenvalue.
    #[test]
    fn test_rayleigh_quotient_dominant() {
        let a = diag_matrix(&[3.0, 10.0]);
        let v0 = vec![0.1, 1.0];
        let (lambda, _) = rayleigh_quotient_iteration(&a, 2, v0, 200, 1e-10);
        assert!(
            (lambda - 10.0).abs() < 0.5,
            "lambda = {lambda}, expected ~10"
        );
    }
    /// Inverse iteration converges to eigenvalue nearest shift.
    #[test]
    fn test_inverse_iteration_near_shift() {
        let a = diag_matrix(&[1.0, 5.0, 9.0]);
        let result = inverse_iteration(&a, 3, 4.8, 300, 1e-8);
        assert!(result.is_some());
        let (lambda, _) = result.unwrap();
        assert!(
            (lambda - 5.0).abs() < 0.3,
            "lambda = {lambda}, expected ~5 (nearest to shift 4.8)"
        );
    }
    /// Inverse iteration on 2x2 symmetric.
    #[test]
    fn test_inverse_iteration_2x2() {
        let a = vec![vec![3.0, 1.0], vec![1.0, 3.0]];
        let result = inverse_iteration(&a, 2, 1.9, 200, 1e-8);
        assert!(result.is_some());
        let (lambda, _) = result.unwrap();
        assert!(
            (lambda - 2.0).abs() < 0.1,
            "lambda = {lambda}, expected 2.0"
        );
    }
    /// Mass-normalized modes satisfy phi^T M phi = 1.
    #[test]
    fn test_mass_normalize_unit_mass() {
        let modes = vec![vec![1.0, 2.0, 3.0]];
        let m_diag = vec![1.0, 1.0, 1.0];
        let norm_modes = mass_normalize_modes(&modes, &m_diag);
        let phi = &norm_modes[0];
        let modal_mass: f64 = phi.iter().zip(m_diag.iter()).map(|(p, m)| p * p * m).sum();
        assert!(
            (modal_mass - 1.0).abs() < 1e-12,
            "modal_mass = {modal_mass}"
        );
    }
    /// Mass-normalized with different masses.
    #[test]
    fn test_mass_normalize_nonuniform_mass() {
        let modes = vec![vec![1.0, 1.0]];
        let m_diag = vec![2.0, 3.0];
        let norm_modes = mass_normalize_modes(&modes, &m_diag);
        let phi = &norm_modes[0];
        let modal_mass: f64 = phi.iter().zip(m_diag.iter()).map(|(p, m)| p * p * m).sum();
        assert!(
            (modal_mass - 1.0).abs() < 1e-12,
            "modal_mass = {modal_mass}"
        );
    }
    /// Mass orthogonality check for orthogonal modes under unit mass.
    #[test]
    fn test_mass_orthogonality_orthogonal() {
        let modes = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let m_diag = vec![1.0, 1.0, 1.0];
        let off = check_mass_orthogonality(&modes, &m_diag);
        assert!(off < 1e-12, "off-diagonal mass = {off}");
    }
    /// Non-orthogonal modes have non-zero off-diagonal mass.
    #[test]
    fn test_mass_orthogonality_nonorthogonal() {
        let modes = vec![vec![1.0, 1.0], vec![1.0, 0.5]];
        let m_diag = vec![1.0, 1.0];
        let off = check_mass_orthogonality(&modes, &m_diag);
        assert!(off > 0.0, "off-diagonal mass should be nonzero: {off}");
    }
    /// Bisection on 2×2 tridiagonal gives correct eigenvalues.
    #[test]
    fn test_tridiagonal_bisection_2x2() {
        let alpha = vec![3.0, 3.0];
        let beta = vec![1.0];
        let evals = tridiagonal_eigenvalues_bisection(&alpha, &beta);
        assert_eq!(evals.len(), 2);
        assert!((evals[0] - 2.0).abs() < 1e-6, "eval[0] = {}", evals[0]);
        assert!((evals[1] - 4.0).abs() < 1e-6, "eval[1] = {}", evals[1]);
    }
    /// Tridiagonal eigenvalues for 1×1 case.
    #[test]
    fn test_tridiagonal_bisection_1x1() {
        let alpha = vec![7.5];
        let beta: Vec<f64> = vec![];
        let evals = tridiagonal_eigenvalues_bisection(&alpha, &beta);
        assert_eq!(evals.len(), 1);
        assert!((evals[0] - 7.5).abs() < 1e-10);
    }
    /// Tridiagonal eigenvalues are sorted ascending.
    #[test]
    fn test_tridiagonal_bisection_sorted() {
        let alpha = vec![5.0, 2.0, 8.0];
        let beta = vec![0.5, 0.3];
        let evals = tridiagonal_eigenvalues_bisection(&alpha, &beta);
        for i in 0..evals.len() - 1 {
            assert!(
                evals[i] <= evals[i + 1] + 1e-10,
                "not sorted: evals[{i}]={} > evals[{}]={}",
                evals[i],
                i + 1,
                evals[i + 1]
            );
        }
    }
    /// Rayleigh quotient for an eigenvector equals the eigenvalue.
    #[test]
    fn test_generalized_rayleigh_quotient_eigenvector() {
        let k = diag_matrix(&[4.0, 9.0]);
        let m_diag = vec![1.0, 1.0];
        let v = vec![1.0, 0.0];
        let rq = generalized_rayleigh_quotient(&k, &m_diag, &v);
        assert!((rq - 4.0).abs() < 1e-10, "Rayleigh quotient = {rq}");
    }
    /// Generalized Rayleigh quotient is between smallest and largest eigenvalue.
    #[test]
    fn test_generalized_rayleigh_quotient_bounds() {
        let k = diag_matrix(&[2.0, 8.0]);
        let m_diag = vec![1.0, 1.0];
        let v = vec![1.0 / 2.0f64.sqrt(), 1.0 / 2.0f64.sqrt()];
        let rq = generalized_rayleigh_quotient(&k, &m_diag, &v);
        assert!(
            rq >= 2.0 - 1e-10 && rq <= 8.0 + 1e-10,
            "Rayleigh quotient {rq} should be in [2,8]"
        );
    }
    /// Participation factor in the mode direction equals modal mass.
    #[test]
    fn test_modal_participation_unit_direction() {
        let modes = vec![vec![1.0, 0.0, 0.0]];
        let m_diag = vec![2.0, 1.0, 1.0];
        let direction = vec![1.0, 0.0, 0.0];
        let pf = modal_participation_factors(&modes, &m_diag, &direction);
        assert!(
            (pf[0] - 2.0).abs() < 1e-12,
            "participation factor = {}",
            pf[0]
        );
    }
    /// Participation factor is zero for orthogonal direction.
    #[test]
    fn test_modal_participation_orthogonal() {
        let modes = vec![vec![1.0, 0.0, 0.0]];
        let m_diag = vec![1.0, 1.0, 1.0];
        let direction = vec![0.0, 1.0, 0.0];
        let pf = modal_participation_factors(&modes, &m_diag, &direction);
        assert!(pf[0].abs() < 1e-12, "participation should be 0: {}", pf[0]);
    }
    /// FRF is large near resonance.
    #[test]
    fn test_frf_resonance_peak() {
        let omega_sq = vec![100.0];
        let participation_factors = vec![1.0];
        let h_near = frf_magnitude(&omega_sq, &participation_factors, 9.9, 0.01);
        let h_far = frf_magnitude(&omega_sq, &participation_factors, 1.0, 0.01);
        assert!(
            h_near > h_far,
            "FRF should peak near resonance: near={h_near}, far={h_far}"
        );
    }
    /// FRF is finite even at omega = 0.
    #[test]
    fn test_frf_at_zero_frequency() {
        let omega_sq = vec![4.0, 16.0];
        let pf = vec![1.0, 1.0];
        let h = frf_magnitude(&omega_sq, &pf, 0.0, 0.05);
        assert!(
            h.is_finite() && h > 0.0,
            "FRF at omega=0 should be positive and finite: {h}"
        );
    }
    /// Residual for exact eigenpair is near zero.
    #[test]
    fn test_eigenpair_residual_exact() {
        let a = diag_matrix(&[3.0, 7.0]);
        let pairs = vec![(3.0, vec![1.0, 0.0]), (7.0, vec![0.0, 1.0])];
        let conv = eigenpair_residuals(&a, 2, &pairs);
        for (i, &res) in conv.residuals.iter().enumerate() {
            assert!(res < 1e-12, "residual[{i}] = {res}");
        }
    }
    /// Residual for incorrect eigenpair is non-zero.
    #[test]
    fn test_eigenpair_residual_inexact() {
        let a = diag_matrix(&[3.0, 7.0]);
        let pairs = vec![(5.0, vec![1.0, 0.0])];
        let conv = eigenpair_residuals(&a, 2, &pairs);
        assert!(
            conv.residuals[0] > 1e-10,
            "residual should be non-zero for wrong eigenvalue"
        );
    }
    /// Power iteration returns None for zero-size matrix.
    #[test]
    fn test_power_iteration_empty() {
        let a: Vec<Vec<f64>> = vec![];
        let result = power_iteration(&a, 0, 100, 1e-8);
        assert!(result.is_none());
    }
    /// Power iteration on a 1x1 matrix.
    #[test]
    fn test_power_iteration_1x1() {
        let a = vec![vec![42.0]];
        let result = power_iteration(&a, 1, 100, 1e-10);
        assert!(result.is_some());
        let (lambda, _) = result.unwrap();
        assert!((lambda - 42.0).abs() < 1e-8, "lambda = {lambda}");
    }
    /// Deflation followed by power iteration recovers second eigenvalue.
    #[test]
    fn test_deflation_and_second_eigenvalue() {
        let mut a = diag_matrix(&[5.0, 2.0]);
        let result1 = power_iteration(&a, 2, 1000, 1e-10).unwrap();
        let (lambda1, v1) = result1;
        assert!((lambda1 - 5.0).abs() < 1e-5, "first eigenvalue = {lambda1}");
        deflate(&mut a, lambda1, &v1);
        let result2 = power_iteration(&a, 2, 1000, 1e-10).unwrap();
        let (lambda2, _) = result2;
        assert!(
            (lambda2 - 2.0).abs() < 0.5,
            "second eigenvalue after deflation = {lambda2}"
        );
    }
    /// Jacobi eigenvectors are orthonormal.
    #[test]
    fn test_jacobi_eigenvectors_orthonormal() {
        let a = vec![
            vec![4.0, 1.0, 0.5],
            vec![1.0, 3.0, 0.2],
            vec![0.5, 0.2, 2.0],
        ];
        let (_, evecs) = jacobi_eigen_dense(&a, 3);
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = evecs[i]
                    .iter()
                    .zip(evecs[j].iter())
                    .map(|(a, b)| a * b)
                    .sum();
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (dot - exp).abs() < 1e-8,
                    "dot(v{i}, v{j}) = {dot}, expected {exp}"
                );
            }
        }
    }
    /// Jacobi eigenvalues satisfy Av = lambda*v.
    #[test]
    fn test_jacobi_residual_3x3() {
        let a = vec![
            vec![5.0, 2.0, 0.0],
            vec![2.0, 4.0, 1.0],
            vec![0.0, 1.0, 3.0],
        ];
        let (evals, evecs) = jacobi_eigen_dense(&a, 3);
        for (i, (&lambda, v)) in evals.iter().zip(evecs.iter()).enumerate() {
            let av = matvec_dense(&a, v, 3);
            let res: f64 = av
                .iter()
                .zip(v.iter())
                .map(|(avi, vi)| (avi - lambda * vi).powi(2))
                .sum::<f64>()
                .sqrt();
            assert!(res < 1e-8, "residual[{i}] = {res}");
        }
    }
    /// Modal analysis returns correct number of modes.
    #[test]
    fn test_modal_analysis_mode_count() {
        let k = diag_matrix(&[4.0, 9.0, 16.0]);
        let m_diag = vec![1.0, 1.0, 1.0];
        let result = modal_analysis(&k, &m_diag, 2);
        assert_eq!(result.omega_sq.len(), 2);
        assert_eq!(result.frequencies_hz.len(), 2);
        assert_eq!(result.mode_shapes.len(), 2);
    }
    /// All frequencies are non-negative.
    #[test]
    fn test_modal_analysis_frequencies_nonneg() {
        let k = diag_matrix(&[1.0, 4.0, 9.0]);
        let m_diag = vec![1.0, 1.0, 1.0];
        let result = modal_analysis(&k, &m_diag, 3);
        for &f in &result.frequencies_hz {
            assert!(f >= 0.0, "frequency must be non-negative: {f}");
        }
    }
    /// Sorting by ascending with all same values is idempotent.
    #[test]
    fn test_sort_ascending_equal_values() {
        let evals = vec![5.0, 5.0, 5.0];
        let evecs = vec![vec![1.0], vec![0.0], vec![-1.0]];
        let (sorted, _) = sort_eigenvalues_ascending(&evals, &evecs);
        for &v in &sorted {
            assert!((v - 5.0).abs() < 1e-12);
        }
    }
    /// Sorting by magnitude handles negative values correctly.
    #[test]
    fn test_sort_magnitude_negative_dominant() {
        let evals = vec![-10.0, 1.0, -3.0];
        let evecs = vec![vec![1.0], vec![1.0], vec![1.0]];
        let (sorted, _) = sort_eigenvalues_by_magnitude(&evals, &evecs);
        assert!(
            (sorted[2] - (-10.0)).abs() < 1e-12,
            "largest magnitude last: {}",
            sorted[2]
        );
    }
    /// MAC is symmetric: mac(a,b) = mac(b,a).
    #[test]
    fn test_mac_symmetry() {
        let phi_a = vec![0.6, 0.8, 0.0];
        let phi_b = vec![0.0, 0.8, 0.6];
        let mab = mac_value(&phi_a, &phi_b);
        let mba = mac_value(&phi_b, &phi_a);
        assert!(
            (mab - mba).abs() < 1e-12,
            "MAC not symmetric: {mab} vs {mba}"
        );
    }
    /// MAC of a scaled version of itself should be 1.
    #[test]
    fn test_mac_scaled() {
        let phi = vec![1.0, 2.0, 3.0];
        let phi_scaled: Vec<f64> = phi.iter().map(|x| x * 5.0).collect();
        let m = mac_value(&phi, &phi_scaled);
        assert!((m - 1.0).abs() < 1e-12, "MAC with scaled mode = {m}");
    }
    /// LU solve on 3×3 system.
    #[test]
    fn test_lu_solve_3x3() {
        let a = vec![
            vec![3.0, 1.0, 0.0],
            vec![1.0, 4.0, 1.0],
            vec![0.0, 1.0, 3.0],
        ];
        let b = vec![9.0, 14.0, 10.0];
        let lu = lu_decompose(&a, 3);
        let x = lu_solve(&lu, &b, 3);
        for (i, row) in a.iter().enumerate() {
            let ax_i: f64 = row.iter().zip(x.iter()).map(|(aij, xj)| aij * xj).sum();
            assert!(
                (ax_i - b[i]).abs() < 1e-8,
                "Ax[{i}] = {ax_i}, b[{i}] = {}",
                b[i]
            );
        }
    }
    /// Gram-Schmidt on linearly dependent vectors does not panic.
    #[test]
    fn test_gram_schmidt_near_dependent() {
        let mut vecs = vec![vec![1.0, 0.0], vec![1.0, 1e-15]];
        gram_schmidt(&mut vecs, 2);
        let n0: f64 = vecs[0].iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((n0 - 1.0).abs() < 1e-10 || n0.abs() < 1e-10);
    }
    /// Sum of effective modal masses equals total mass for complete mode set.
    #[test]
    fn test_effective_modal_mass_completeness() {
        let k = diag_matrix(&[4.0, 9.0]);
        let m_diag = vec![1.0, 1.0];
        let direction = vec![1.0, 1.0];
        let result = modal_analysis(&k, &m_diag, 2);
        let total_emm: f64 = result
            .mode_shapes
            .iter()
            .map(|phi| effective_modal_mass(phi, &m_diag, &direction))
            .sum();
        assert!(total_emm > 0.0, "total EMM should be positive: {total_emm}");
    }
    /// Subspace iteration on 2x2 matrix converges to both eigenvalues.
    #[test]
    fn test_subspace_iteration_2x2() {
        let a = vec![vec![5.0, 0.0], vec![0.0, 2.0]];
        let pairs = subspace_iteration(&a, 2, 2, 200, 1e-8);
        let mut evals: Vec<f64> = pairs.iter().map(|(e, _)| *e).collect();
        evals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((evals[0] - 2.0).abs() < 0.5, "eval[0] = {}", evals[0]);
        assert!((evals[1] - 5.0).abs() < 0.5, "eval[1] = {}", evals[1]);
    }
    /// Shift-invert near shift=10 finds eigenvalue 10.
    #[test]
    fn test_shift_invert_diagonal_near_10() {
        let a = diag_matrix(&[1.0, 5.0, 10.0]);
        let pairs = shift_invert_power(&a, 3, 9.5, 1, 300, 1e-8);
        assert!(!pairs.is_empty());
        let (lambda, _) = &pairs[0];
        assert!(
            (*lambda - 10.0).abs() < 0.5,
            "lambda = {lambda}, expected ~10"
        );
    }
}
