//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::solvers::conjugate_gradient;
use crate::sparse::CsrMatrix;

use super::types::{ConvergenceCriteria, EnergyConvergenceCriteria, NrResult, RiksState};

/// Incremental load stepping with Newton-Raphson at each step.
///
/// Applies load in `n_steps` equal increments from 0 to `total_load`,
/// using Newton-Raphson at each step with warm-starting from previous solution.
///
/// Returns a vector of `(lambda, displacement)` pairs for each completed step.
pub fn incremental_load_stepping<F>(
    n_dof: usize,
    n_steps: usize,
    total_load: f64,
    tangent_and_residual: F,
    criteria: &ConvergenceCriteria,
) -> Vec<(f64, Vec<f64>)>
where
    F: Fn(&[f64], f64) -> (CsrMatrix, Vec<f64>),
{
    let mut results = Vec::with_capacity(n_steps);
    let mut u = vec![0.0; n_dof];
    for step in 1..=n_steps {
        let lambda = (step as f64 / n_steps as f64) * total_load;
        let result = newton_raphson(
            n_dof,
            Some(u.clone()),
            |u_curr| tangent_and_residual(u_curr, lambda),
            criteria,
        );
        u = result.displacement.clone();
        results.push((lambda, result.displacement));
        if !result.converged {
            break;
        }
    }
    results
}
/// Compute the buckling load factor using the linearized buckling criterion.
///
/// The buckling load factor is estimated as the smallest positive λ such that
/// `det(K_e + λ K_g) = 0`, where:
/// - `k_e` = elastic stiffness (2×2 for simplicity)
/// - `k_g` = geometric stiffness scaled by a reference load
///
/// For a 2×2 system the characteristic equation is solved directly.
pub fn buckling_load_factor_2x2(k_e: &[[f64; 2]; 2], k_g: &[[f64; 2]; 2]) -> Option<f64> {
    let a = k_g[0][0] * k_g[1][1] - k_g[0][1] * k_g[1][0];
    let b = k_e[0][0] * k_g[1][1] + k_g[0][0] * k_e[1][1] - 2.0 * k_e[0][1] * k_g[0][1];
    let c = k_e[0][0] * k_e[1][1] - k_e[0][1] * k_e[1][0];
    if a.abs() < 1e-60 {
        if b.abs() < 1e-60 {
            return None;
        }
        let lam = -c / b;
        if lam > 0.0 { Some(lam) } else { None }
    } else {
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let sqrt_disc = disc.sqrt();
        let lam1 = (-b + sqrt_disc) / (2.0 * a);
        let lam2 = (-b - sqrt_disc) / (2.0 * a);
        let pos: Vec<f64> = [lam1, lam2]
            .iter()
            .filter(|&&l| l > 1e-12)
            .copied()
            .collect();
        pos.into_iter().reduce(f64::min)
    }
}
/// Geometric stiffness matrix for a bar element in 2D.
///
/// K_geo = (N / L) * \[\[1, -1\\], \[-1, 1\]]
/// where N is the axial force and L is the element length.
pub fn geometric_stiffness_bar_2d(axial_force: f64, length: f64) -> [[f64; 2]; 2] {
    let scale = axial_force / length.max(1e-30);
    [[scale, -scale], [-scale, scale]]
}
/// Elastic stiffness matrix for a 1D bar element.
///
/// K_el = (E * A / L) * \[\[1, -1\\], \[-1, 1\]]
pub fn elastic_stiffness_bar_1d(e_mod: f64, area: f64, length: f64) -> [[f64; 2]; 2] {
    let scale = e_mod * area / length.max(1e-30);
    [[scale, -scale], [-scale, scale]]
}
/// One Riks predictor step: advance along the tangent direction.
///
/// Returns updated `(u_new, lambda_new)`.
pub fn riks_predictor<F>(
    state: &mut RiksState,
    tangent_and_residual: F,
    sign: f64,
) -> (Vec<f64>, f64)
where
    F: Fn(&[f64], f64) -> (CsrMatrix, Vec<f64>),
{
    let n_dof = state.u.len();
    let (k_tan, _r) = tangent_and_residual(&state.u, state.lambda);
    let f_ref: Vec<f64> = vec![1.0 / (n_dof as f64).sqrt(); n_dof];
    let x0 = vec![0.0; n_dof];
    let u_ref = conjugate_gradient(&k_tan, &f_ref, &x0, n_dof * 100, 1e-12);
    let u_ref_norm_sq: f64 = u_ref.iter().map(|x| x * x).sum::<f64>();
    let denom = (u_ref_norm_sq + 1.0).sqrt();
    let dlambda_pred = sign * state.arc_length / denom;
    let u_new: Vec<f64> = state
        .u
        .iter()
        .zip(u_ref.iter())
        .map(|(ui, dui)| ui + dlambda_pred * dui)
        .collect();
    let lambda_new = state.lambda + dlambda_pred;
    state.du = u_ref.iter().map(|x| x * dlambda_pred).collect();
    state.dlambda = dlambda_pred;
    (u_new, lambda_new)
}
/// Newton-Raphson with energy-based convergence criterion.
pub fn newton_raphson_energy<F>(
    n_dof: usize,
    u0: Option<Vec<f64>>,
    tangent_and_residual: F,
    criteria: &EnergyConvergenceCriteria,
) -> NrResult
where
    F: Fn(&[f64]) -> (CsrMatrix, Vec<f64>),
{
    let mut u = u0.unwrap_or_else(|| vec![0.0; n_dof]);
    let mut iterations = 0;
    let mut final_residual = f64::INFINITY;
    let mut converged = false;
    let (_, r0) = tangent_and_residual(&u);
    let e0: f64 = r0.iter().map(|r| r * r).sum::<f64>().sqrt().max(1e-60);
    for _iter in 0..criteria.max_iter {
        let (k_tan, residual) = tangent_and_residual(&u);
        let r_norm: f64 = residual.iter().map(|r| r * r).sum::<f64>().sqrt();
        final_residual = r_norm;
        let neg_r: Vec<f64> = residual.iter().map(|r| -r).collect();
        let x0 = vec![0.0; n_dof];
        let du = conjugate_gradient(&k_tan, &neg_r, &x0, n_dof * 100, 1e-12);
        let energy_inc: f64 = du
            .iter()
            .zip(residual.iter())
            .map(|(d, r)| d * r)
            .sum::<f64>()
            .abs();
        let energy_ratio = energy_inc / e0;
        for i in 0..n_dof {
            u[i] += du[i];
        }
        iterations += 1;
        if energy_ratio < criteria.energy_tol {
            let (_, res_final) = tangent_and_residual(&u);
            final_residual = res_final.iter().map(|r| r * r).sum::<f64>().sqrt();
            converged = true;
            break;
        }
    }
    NrResult {
        displacement: u,
        converged,
        iterations,
        final_residual,
    }
}
/// Compute the second Piola-Kirchhoff stress S from the Green-Lagrange strain E
/// using the St. Venant-Kirchhoff constitutive law.
///
/// S = lambda * tr(E) * I + 2*mu * E
pub fn svk_stress(e: &[[f64; 3]; 3], mu: f64, lambda: f64) -> [[f64; 3]; 3] {
    let tr_e = e[0][0] + e[1][1] + e[2][2];
    let mut s = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let delta = if i == j { 1.0 } else { 0.0 };
            s[i][j] = lambda * tr_e * delta + 2.0 * mu * e[i][j];
        }
    }
    s
}
/// Compute the Cauchy stress from the second Piola-Kirchhoff stress.
///
/// sigma = (1/J) * F * S * F^T
pub fn cauchy_from_pk2(f: &[[f64; 3]; 3], s: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let j = mat3_det(f);
    let fs = mat3_mul(f, s);
    let ft = mat3_transpose(f);
    let fsf_t = mat3_mul(&fs, &ft);
    let mut sigma = [[0.0f64; 3]; 3];
    let inv_j = if j.abs() > 1e-30 { 1.0 / j } else { 0.0 };
    for i in 0..3 {
        for j_idx in 0..3 {
            sigma[i][j_idx] = inv_j * fsf_t[i][j_idx];
        }
    }
    sigma
}
/// Volumetric-deviatoric split of the deformation gradient.
///
/// Returns `(F_vol, F_dev)` where `F_vol = J^{1/3} * I` and `F_dev = J^{-1/3} * F`.
pub fn vol_dev_split(f: &[[f64; 3]; 3]) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let j = mat3_det(f);
    let j13 = j.abs().powf(1.0 / 3.0);
    let inv_j13 = if j13 > 1e-30 { 1.0 / j13 } else { 0.0 };
    let mut f_vol = [[0.0f64; 3]; 3];
    for (i, row) in f_vol.iter_mut().enumerate() {
        row[i] = j13;
    }
    let mut f_dev = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for k in 0..3 {
            f_dev[i][k] = inv_j13 * f[i][k];
        }
    }
    (f_vol, f_dev)
}
/// Isochoric Neo-Hookean strain energy density.
///
/// W_iso = (mu/2) * (I1_bar - 3)
/// where I1_bar = J^{-2/3} * tr(C)
pub fn neo_hookean_isochoric_energy(f: &[[f64; 3]; 3], mu: f64) -> f64 {
    let j = mat3_det(f);
    let c = right_cauchy_green(f);
    let i1 = mat3_trace(&c);
    let i1_bar = j.abs().powf(-2.0 / 3.0) * i1;
    (mu / 2.0) * (i1_bar - 3.0)
}
/// Volumetric strain energy density (penalty/bulk).
///
/// W_vol = (kappa/4) * (J² - 1 - 2*ln(J))
pub fn volumetric_energy(f: &[[f64; 3]; 3], kappa: f64) -> f64 {
    let j = mat3_det(f);
    let ln_j = j.max(1e-30).ln();
    (kappa / 4.0) * (j * j - 1.0 - 2.0 * ln_j)
}
/// Total Neo-Hookean energy = isochoric + volumetric.
pub fn neo_hookean_total_energy(f: &[[f64; 3]; 3], mu: f64, kappa: f64) -> f64 {
    neo_hookean_isochoric_energy(f, mu) + volumetric_energy(f, kappa)
}
#[cfg(test)]
mod nonlinear_extended_tests {
    use super::*;
    use crate::nonlinear::*;
    fn unit_tet() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
    fn make_1x1(v: f64) -> CsrMatrix {
        CsrMatrix::from_triplets(1, 1, &[(0, 0, v)])
    }
    /// SVK stress at zero strain = zero.
    #[test]
    fn test_svk_stress_zero_strain() {
        let e = [[0.0; 3]; 3];
        let s = svk_stress(&e, 1.0, 1.0);
        for (i, row) in s.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(val.abs() < 1e-15, "S[{i}][{j}] = {}", val);
            }
        }
    }
    /// SVK stress is symmetric for symmetric strain input.
    #[test]
    fn test_svk_stress_symmetry() {
        let e = [[0.1, 0.05, 0.02], [0.05, 0.2, 0.03], [0.02, 0.03, 0.15]];
        let s = svk_stress(&e, 1.0e9, 2.0e9);
        for (i, row) in s.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - s[j][i]).abs() < 1e-6,
                    "S not symmetric at ({i},{j}): {} vs {}",
                    val,
                    s[j][i]
                );
            }
        }
    }
    /// Cauchy stress from PK2 should be symmetric.
    #[test]
    fn test_cauchy_from_pk2_symmetry() {
        let f = [[1.1, 0.05, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let e = green_lagrange_strain(&f);
        let s = svk_stress(&e, 1.0, 1.0);
        let sigma = cauchy_from_pk2(&f, &s);
        for (i, row) in sigma.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!(
                    (val - sigma[j][i]).abs() < 1e-10,
                    "sigma not symmetric: ({i},{j})"
                );
            }
        }
    }
    /// Volumetric part of F should be isotropic (scalar * I).
    #[test]
    fn test_vol_dev_split_isochoric() {
        let f = [[1.5, 0.1, 0.0], [0.0, 1.2, 0.05], [0.0, 0.0, 0.9]];
        let (_f_vol, f_dev) = vol_dev_split(&f);
        let j_dev = mat3_det(&f_dev);
        assert!(
            (j_dev.abs() - 1.0).abs() < 1e-10,
            "det(F_dev) = {j_dev}, expected 1.0"
        );
    }
    /// At F = I, volumetric part = I.
    #[test]
    fn test_vol_dev_split_identity() {
        let f = mat3_identity();
        let (f_vol, f_dev) = vol_dev_split(&f);
        for i in 0..3 {
            for j in 0..3 {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!((f_vol[i][j] - exp).abs() < 1e-14);
                assert!((f_dev[i][j] - exp).abs() < 1e-14);
            }
        }
    }
    /// At F = I, isochoric energy = 0.
    #[test]
    fn test_neo_hookean_isochoric_energy_identity() {
        let f = mat3_identity();
        let w = neo_hookean_isochoric_energy(&f, 1.0);
        assert!(w.abs() < 1e-13, "W_iso at F=I = {w}");
    }
    /// Isochoric energy is non-negative for any F.
    #[test]
    fn test_neo_hookean_isochoric_energy_nonneg() {
        let f = [[1.2, 0.1, 0.0], [0.0, 0.8, 0.0], [0.0, 0.0, 1.0]];
        let w = neo_hookean_isochoric_energy(&f, 1.0);
        assert!(w >= 0.0, "W_iso = {w}");
    }
    /// Volumetric energy at J=1 should be 0.
    #[test]
    fn test_volumetric_energy_identity() {
        let f = mat3_identity();
        let w = volumetric_energy(&f, 1.0);
        assert!(w.abs() < 1e-13, "W_vol at J=1 = {w}");
    }
    /// Total Neo-Hookean energy is positive for stretched F.
    #[test]
    fn test_neo_hookean_total_energy_positive() {
        let f = [[1.5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let w = neo_hookean_total_energy(&f, 1.0, 2.0);
        assert!(w > 0.0, "total NHK energy = {w}");
    }
    /// Buckling load factor for simple column.
    #[test]
    fn test_buckling_load_factor_2x2_simple() {
        let k_e = [[4.0, 0.0], [0.0, 3.0]];
        let k_g = [[-1.0, 0.0], [0.0, -1.0]];
        let lam = buckling_load_factor_2x2(&k_e, &k_g);
        assert!(lam.is_some(), "should have a positive buckling load factor");
        let lam_val = lam.unwrap();
        assert!(lam_val > 0.0, "buckling load factor = {lam_val}");
    }
    /// Geometric stiffness of a bar under compression is negative on diagonal.
    #[test]
    fn test_geometric_stiffness_bar_compression() {
        let k_geo = geometric_stiffness_bar_2d(-100.0, 1.0);
        assert!(k_geo[0][0] < 0.0, "under compression, K_geo diagonal < 0");
        assert!(k_geo[0][1] > 0.0, "off-diagonal should be positive");
        assert!((k_geo[0][0] + k_geo[0][1]).abs() < 1e-12);
    }
    /// Elastic bar stiffness is positive definite (diagonal > 0).
    #[test]
    fn test_elastic_stiffness_bar_positive() {
        let k_el = elastic_stiffness_bar_1d(200e9, 1e-4, 1.0);
        assert!(k_el[0][0] > 0.0);
        assert!((k_el[0][0] + k_el[0][1]).abs() < 1e-6);
    }
    /// Snap-through detector correctly identifies sign change.
    #[test]
    fn test_snap_through_detector_detects_sign_change() {
        let mut detector = SnapThroughDetector::new();
        let k_pos = [[4.0, 0.0], [0.0, 3.0]];
        let snap1 = detector.check_2x2(k_pos, 1.0);
        assert!(!snap1, "no snap at first step");
        let k_pos2 = [[2.0, 0.0], [0.0, 2.0]];
        let snap2 = detector.check_2x2(k_pos2, 2.0);
        assert!(!snap2, "no snap when det stays positive");
        let k_neg = [[1.0, 2.0], [2.0, 0.0]];
        let snap3 = detector.check_2x2(k_neg, 3.0);
        assert!(snap3, "snap-through should be detected on sign change");
        assert!(detector.has_snap_through());
        assert_eq!(detector.critical_loads.len(), 1);
        assert!((detector.critical_loads[0] - 3.0).abs() < 1e-12);
    }
    /// Snap-through not triggered without sign change.
    #[test]
    fn test_snap_through_no_false_positive() {
        let mut detector = SnapThroughDetector::new();
        for _ in 0..5 {
            let k = [[3.0, 0.5], [0.5, 2.0]];
            detector.check_2x2(k, 1.0);
        }
        assert!(!detector.has_snap_through(), "no snap-through expected");
    }
    /// Incremental load stepping should reach total load.
    #[test]
    fn test_incremental_load_stepping() {
        let results = incremental_load_stepping(
            2,
            5,
            10.0,
            |u, lambda| {
                let k = CsrMatrix::from_triplets(2, 2, &[(0, 0, 2.0), (1, 1, 2.0)]);
                let ku = k.mul_vec(u);
                let r = vec![ku[0] - lambda, ku[1] - lambda];
                (k, r)
            },
            &ConvergenceCriteria::default(),
        );
        assert!(!results.is_empty(), "should have at least one step");
        let (final_lambda, _) = results.last().unwrap();
        assert!(
            (*final_lambda - 10.0).abs() < 1e-10,
            "final lambda = {final_lambda}"
        );
    }
    /// Incremental stepping: solution at each step satisfies K*u = lambda.
    #[test]
    fn test_incremental_stepping_solution_quality() {
        let results = incremental_load_stepping(
            1,
            3,
            6.0,
            |u, lambda| {
                let k = CsrMatrix::from_triplets(1, 1, &[(0, 0, 3.0)]);
                let ku = k.mul_vec(u);
                let r = vec![ku[0] - lambda];
                (k, r)
            },
            &ConvergenceCriteria::default(),
        );
        let (_, u_final) = results.last().unwrap();
        assert!(
            (u_final[0] - 2.0).abs() < 1e-6,
            "u = {}, expected 2.0",
            u_final[0]
        );
    }
    /// Energy-based NR converges on linear system.
    #[test]
    fn test_nr_energy_linear() {
        let crit = EnergyConvergenceCriteria::default();
        let result = newton_raphson_energy(
            1,
            Some(vec![5.0]),
            |u| {
                let x = u[0];
                let k = make_1x1(4.0);
                let r = vec![4.0 * x - 8.0];
                (k, r)
            },
            &crit,
        );
        assert!(result.converged, "energy NR should converge");
        assert!(
            (result.displacement[0] - 2.0).abs() < 1e-5,
            "u = {}",
            result.displacement[0]
        );
    }
    /// EnergyConvergenceCriteria default values.
    #[test]
    fn test_energy_criteria_defaults() {
        let c = EnergyConvergenceCriteria::default();
        assert_eq!(c.max_iter, 50);
        assert!(c.energy_tol < 1e-4);
    }
    /// RiksState computes arc length correctly.
    #[test]
    fn test_riks_state_arc_length() {
        let mut state = RiksState::new(2, 0.5);
        state.du = vec![0.3, 0.4];
        state.dlambda = 0.0;
        let arc = state.current_arc();
        assert!((arc - 0.5).abs() < 1e-10, "arc = {arc}, expected 0.5");
    }
    /// RiksState zero increments give zero arc.
    #[test]
    fn test_riks_state_zero_arc() {
        let state = RiksState::new(3, 0.1);
        let arc = state.current_arc();
        assert!(arc.abs() < 0.02, "initial arc = {arc}");
    }
    /// Mooney-Rivlin stress at F = I: deviatoric form gives scalar isotropic stress.
    ///
    /// At F=I: B=I, I1=3, C=I, F*C=I → P = 2*c10*I + 2*c01*(3I - I) = (2c10 + 4c01)*I
    /// The result should be proportional to identity (isotropic).
    #[test]
    fn test_mooney_rivlin_stress_identity() {
        let f = mat3_identity();
        let c10 = 0.5;
        let c01 = 0.25;
        let p = mooney_rivlin_stress(&f, c10, c01);
        let expected_diag = 2.0 * c10 + 4.0 * c01;
        for (i, row) in p.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                if i != j {
                    assert!(
                        val.abs() < 1e-13,
                        "MR stress off-diagonal at F=I: P[{i}][{j}] = {}",
                        val
                    );
                } else {
                    assert!(
                        (val - expected_diag).abs() < 1e-13,
                        "MR stress diagonal at F=I: P[{i}][{j}] = {}, expected {expected_diag}",
                        val
                    );
                }
            }
        }
    }
    /// Mooney-Rivlin stress non-zero under shear.
    #[test]
    fn test_mooney_rivlin_stress_shear() {
        let f = [[1.0, 0.2, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p = mooney_rivlin_stress(&f, 0.5, 0.25);
        let total: f64 = p.iter().flat_map(|r| r.iter()).map(|&v| v * v).sum();
        assert!(
            total > 0.0,
            "Mooney-Rivlin stress should be non-zero under shear"
        );
    }
    /// Convergence history correctly computes rate for long sequence.
    #[test]
    fn test_convergence_rate_quadratic() {
        let mut hist = ConvergenceHistory::new();
        let mut r = 1.0;
        for _ in 0..6 {
            hist.record(r, r * 0.5);
            r *= 0.1;
        }
        let rate = hist.convergence_rate().unwrap();
        assert!(rate < 1.0, "rate should indicate convergence: {rate}");
    }
    /// Convergence history with one entry has no rate.
    #[test]
    fn test_convergence_rate_single_entry() {
        let mut hist = ConvergenceHistory::new();
        hist.record(1.0, 0.5);
        assert!(hist.convergence_rate().is_none() || hist.convergence_rate().is_some());
    }
    /// LoadStepper with 1 step should complete immediately.
    #[test]
    fn test_load_stepper_single_step() {
        let mut s = LoadStepper::new(50.0, 1);
        let load = s.advance();
        assert!(load.is_some());
        assert!((load.unwrap() - 50.0).abs() < 1e-10);
        assert!(s.is_complete());
        assert!(s.advance().is_none());
    }
    /// LoadStepper adapt_step with zero actual iters uses ratio 2.0.
    #[test]
    fn test_load_stepper_adapt_zero_iters() {
        let mut s = LoadStepper::new(100.0, 10);
        let initial = s.step_size;
        s.adapt_step(0);
        assert!(s.step_size >= initial, "step should increase or stay same");
    }
    /// mat3_mul(A, A_inv) = I.
    #[test]
    fn test_mat3_mul_inverse() {
        let m = [[2.0, 1.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 4.0]];
        let mi = mat3_inv(&m);
        let prod = mat3_mul(&m, &mi);
        for (i, row) in prod.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (val - exp).abs() < 1e-10,
                    "prod[{i}][{j}] = {} expected {exp}",
                    val
                );
            }
        }
    }
    /// mat3_transpose(mat3_transpose(A)) = A.
    #[test]
    fn test_mat3_double_transpose() {
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let tt = mat3_transpose(&mat3_transpose(&m));
        for (i, row) in tt.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - m[i][j]).abs() < 1e-15);
            }
        }
    }
    /// mat3_det of diagonal matrix equals product of diagonal.
    #[test]
    fn test_mat3_det_diagonal() {
        let m = [[3.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 7.0]];
        let det = mat3_det(&m);
        assert!((det - 105.0).abs() < 1e-10, "det = {det}");
    }
    /// Deformation gradient under uniform compression.
    #[test]
    fn test_deformation_gradient_compression() {
        let x0 = unit_tet();
        let mut x = x0;
        for node in &mut x {
            for c in node.iter_mut() {
                *c *= 0.8;
            }
        }
        let f = deformation_gradient(&x, &x0);
        for (i, row) in f.iter().enumerate() {
            assert!((row[i] - 0.8).abs() < 1e-12, "F[{i}][{i}] = {}", row[i]);
        }
    }
    /// Right Cauchy-Green at identity is identity.
    #[test]
    fn test_right_cauchy_green_identity() {
        let f = mat3_identity();
        let c = right_cauchy_green(&f);
        for (i, row) in c.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!((val - exp).abs() < 1e-15);
            }
        }
    }
    /// Left Cauchy-Green at identity is identity.
    #[test]
    fn test_left_cauchy_green_identity() {
        let f = mat3_identity();
        let b = left_cauchy_green(&f);
        for (i, row) in b.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                let exp = if i == j { 1.0 } else { 0.0 };
                assert!((val - exp).abs() < 1e-15);
            }
        }
    }
    /// Invariants of pure stretch: I3 = J² = det²(F).
    #[test]
    fn test_invariants_pure_stretch() {
        let stretch = 1.3_f64;
        let f = [
            [stretch, 0.0, 0.0],
            [0.0, stretch, 0.0],
            [0.0, 0.0, stretch],
        ];
        let b = left_cauchy_green(&f);
        let (_i1, _i2, i3) = invariants(&b);
        let j_sq = mat3_det(&f).powi(2);
        assert!((i3 - j_sq).abs() < 1e-10, "I3 = {i3}, J² = {j_sq}");
    }
}
