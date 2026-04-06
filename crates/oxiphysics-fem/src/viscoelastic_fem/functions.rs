//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#![allow(clippy::items_after_test_module)]
use super::types::{FrequencySweepPoint, ViscoelasticMaterial, WLFShift};

/// Compute reduced time using WLF time-temperature superposition.
///
/// t_reduced = t / 10^(log a_T)
#[allow(dead_code)]
#[allow(non_snake_case)]
pub fn time_temperature_superposition(t: f64, T: f64, wlf: &WLFShift) -> f64 {
    let log_at = wlf.shift_factor(T);
    let at = 10.0_f64.powf(log_at);
    if at.abs() < 1e-30 { t } else { t / at }
}
/// Solve K u = f by Gaussian elimination with partial pivoting.
pub(super) fn solve_tridiagonal_system(k: &[Vec<f64>], f: &[f64]) -> Vec<f64> {
    let n = f.len();
    let mut a: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = k[i].clone();
            row.push(f[i]);
            row
        })
        .collect();
    for col in 0..n {
        let mut max_row = col;
        for row in (col + 1)..n {
            if a[row][col].abs() > a[max_row][col].abs() {
                max_row = row;
            }
        }
        a.swap(col, max_row);
        let pivot = a[col][col];
        if pivot.abs() < 1e-30 {
            continue;
        }
        for row in (col + 1)..n {
            let factor = a[row][col] / pivot;
            for c in col..=n {
                let val = a[col][c];
                a[row][c] -= factor * val;
            }
        }
    }
    let mut u = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = a[i][n];
        for j in (i + 1)..n {
            sum -= a[i][j] * u[j];
        }
        u[i] = if a[i][i].abs() > 1e-30 {
            sum / a[i][i]
        } else {
            0.0
        };
    }
    u
}
/// Numerical creep compliance J(t) via midpoint-rule integration of 1/E(s).
pub fn creep_compliance_integral(material: &ViscoelasticMaterial, t: f64, n_steps: usize) -> f64 {
    if t <= 0.0 || n_steps == 0 {
        return 1.0 / material.instantaneous_modulus();
    }
    let dt = t / n_steps as f64;
    let sum: f64 = (0..n_steps)
        .map(|i| {
            let t_mid = (i as f64 + 0.5) * dt;
            1.0 / material.relaxation_modulus(t_mid)
        })
        .sum();
    sum * dt / t
}
/// Stress relaxation ratio E(t) / E(0).
pub fn stress_relaxation_ratio(material: &ViscoelasticMaterial, t: f64) -> f64 {
    let e0 = material.instantaneous_modulus();
    if e0 == 0.0 {
        return 1.0;
    }
    material.relaxation_modulus(t) / e0
}
/// Compute the relaxation spectrum from Prony series.
///
/// Returns `(tau_i, H_i)` pairs where H(tau) is the relaxation spectrum.
pub fn relaxation_spectrum(material: &ViscoelasticMaterial) -> Vec<(f64, f64)> {
    material
        .prony_terms
        .iter()
        .map(|term| (term.relaxation_time, term.modulus))
        .collect()
}
/// Compute the retardation spectrum from Prony series (approximation).
///
/// Uses the relationship between relaxation and retardation spectra.
pub fn retardation_spectrum(material: &ViscoelasticMaterial) -> Vec<(f64, f64)> {
    let e0 = material.instantaneous_modulus();
    let e_inf = material.e_inf;
    if e0.abs() < 1e-30 || e_inf.abs() < 1e-30 {
        return Vec::new();
    }
    material
        .prony_terms
        .iter()
        .map(|term| {
            let tau_ret = term.relaxation_time * e_inf / e0;
            let d_i = term.modulus / (e0 * e_inf);
            (tau_ret, d_i)
        })
        .collect()
}
/// Compute dynamic complex modulus magnitude |E*(omega)|.
pub fn complex_modulus_magnitude(material: &ViscoelasticMaterial, omega: f64) -> f64 {
    let e_prime = material.storage_modulus(omega);
    let e_double_prime = material.loss_modulus(omega);
    (e_prime * e_prime + e_double_prime * e_double_prime).sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::viscoelastic_fem::*;
    fn polymer() -> ViscoelasticMaterial {
        ViscoelasticMaterial::polymer()
    }
    #[test]
    fn prony_term_fields() {
        let t = PronyTerm {
            modulus: 1.0e9,
            relaxation_time: 5.0,
        };
        assert_eq!(t.modulus, 1.0e9);
        assert_eq!(t.relaxation_time, 5.0);
    }
    #[test]
    fn polymer_has_three_prony_terms() {
        assert_eq!(polymer().prony_terms.len(), 3);
    }
    #[test]
    fn instantaneous_modulus_is_sum() {
        let mat = polymer();
        let expected = mat.e_inf + mat.prony_terms.iter().map(|t| t.modulus).sum::<f64>();
        assert!((mat.instantaneous_modulus() - expected).abs() < 1e-6);
    }
    #[test]
    fn relaxation_modulus_at_zero_equals_instantaneous() {
        let mat = polymer();
        let e0 = mat.instantaneous_modulus();
        let e_t = mat.relaxation_modulus(0.0);
        assert!((e_t - e0).abs() / e0 < 1e-12);
    }
    #[test]
    fn relaxation_modulus_decays_to_e_inf() {
        let mat = polymer();
        let e_large = mat.relaxation_modulus(1e12);
        assert!((e_large - mat.e_inf).abs() / mat.e_inf < 1e-6);
    }
    #[test]
    fn relaxation_modulus_is_monotone_decreasing() {
        let mat = polymer();
        let times = [0.0, 1.0, 10.0, 100.0, 1000.0];
        let vals: Vec<f64> = times.iter().map(|&t| mat.relaxation_modulus(t)).collect();
        for w in vals.windows(2) {
            assert!(w[0] >= w[1], "E(t) should be non-increasing");
        }
    }
    #[test]
    fn relaxation_ratio_at_zero_is_one() {
        let mat = polymer();
        assert!((stress_relaxation_ratio(&mat, 0.0) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn relaxation_ratio_approaches_e_inf_fraction() {
        let mat = polymer();
        let r = stress_relaxation_ratio(&mat, 1e12);
        let expected = mat.e_inf / mat.instantaneous_modulus();
        assert!(
            (r - expected).abs() < 1e-4,
            "ratio should approach E_inf/E_0 at large t"
        );
    }
    #[test]
    fn internal_variable_initialises_to_zero() {
        let iv = InternalVariable::new(3);
        assert!(iv.values.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn internal_variable_n_terms_matches() {
        let iv = InternalVariable::new(5);
        assert_eq!(iv.n_terms, 5);
        assert_eq!(iv.values.len(), 5);
    }
    #[test]
    fn internal_variable_update_nonzero() {
        let mat = polymer();
        let mut iv = InternalVariable::new(mat.prony_terms.len());
        iv.update(0.01, &mat, 1.0);
        assert!(iv.values.iter().any(|&v| v != 0.0));
    }
    #[test]
    fn internal_variable_decays_without_strain() {
        let mat = polymer();
        let mut iv = InternalVariable::new(mat.prony_terms.len());
        iv.update(0.01, &mat, 1.0);
        let loaded: Vec<f64> = iv.values.clone();
        iv.update(0.0, &mat, 1000.0);
        for (old, new) in loaded.iter().zip(iv.values.iter()) {
            assert!(new.abs() <= old.abs() + 1e-20, "variables must decay");
        }
    }
    #[test]
    fn internal_variable_reset() {
        let mat = polymer();
        let mut iv = InternalVariable::new(mat.prony_terms.len());
        iv.update(0.01, &mat, 1.0);
        iv.reset();
        assert!(iv.values.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn element_stiffness_matrix_symmetric() {
        let mat = polymer();
        let elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        let ke = elem.stiffness_matrix(1.0);
        assert!((ke[0][1] - ke[1][0]).abs() < 1e-10);
    }
    #[test]
    fn element_stiffness_matrix_zero_row_sum() {
        let mat = polymer();
        let elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        let ke = elem.stiffness_matrix(1.0);
        for row in &ke {
            let sum: f64 = row.iter().sum();
            assert!(sum.abs() < 1e-6, "row sum must be zero for rigid-body mode");
        }
    }
    #[test]
    fn element_stiffness_positive_diagonal() {
        let mat = polymer();
        let elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        let ke = elem.stiffness_matrix(0.5);
        assert!(ke[0][0] > 0.0);
        assert!(ke[1][1] > 0.0);
    }
    #[test]
    fn element_stress_update_nonzero() {
        let mat = polymer();
        let mut elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        let sigma = elem.stress_update(0.001, 1.0);
        assert!(sigma.abs() > 0.0);
    }
    #[test]
    fn element_stress_relaxes_over_time() {
        let mat = polymer();
        let mut elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        let sigma0 = elem.stress_update(0.01, 0.01);
        let mut sigma = sigma0;
        for _ in 0..50 {
            sigma = elem.stress_update(0.0, 1.0);
        }
        assert!(sigma < sigma0, "stress should relax under constant strain");
    }
    #[test]
    fn algorithmic_modulus_large_dt_approaches_e_inf() {
        let mat = polymer();
        let elem = ViscoelasticElement1D::new(1.0, 1e-4, mat.clone());
        let e_alg = elem.algorithmic_modulus(1e12);
        assert!((e_alg - mat.e_inf).abs() / mat.e_inf < 1e-4);
    }
    #[test]
    fn algorithmic_modulus_small_dt_approaches_instantaneous() {
        let mat = polymer();
        let elem = ViscoelasticElement1D::new(1.0, 1e-4, mat.clone());
        let e_alg = elem.algorithmic_modulus(1e-16);
        let e0 = mat.instantaneous_modulus();
        assert!((e_alg - e0).abs() / e0 < 1e-9);
    }
    #[test]
    fn beam_n_nodes_correct() {
        let beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        assert_eq!(beam.n_nodes, 5);
    }
    #[test]
    fn beam_stiffness_matrix_size() {
        let n = 3usize;
        let beam = ViscoelasticBeam1D::new(n, 1.0, 1e-4, polymer());
        let k = beam.assemble_stiffness(1.0);
        assert_eq!(k.len(), n + 1);
        assert_eq!(k[0].len(), n + 1);
    }
    #[test]
    fn beam_global_stiffness_symmetric() {
        let beam = ViscoelasticBeam1D::new(3, 1.0, 1e-4, polymer());
        let k = beam.assemble_stiffness(1.0);
        let n = k.len();
        for i in 0..n {
            for j in 0..n {
                assert!((k[i][j] - k[j][i]).abs() < 1e-6);
            }
        }
    }
    #[test]
    fn beam_solve_step_returns_correct_length() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        let u = beam.solve_step(&forces, 1.0);
        assert_eq!(u.len(), 5);
    }
    #[test]
    fn beam_left_end_fixed_displacement_near_zero() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        let u = beam.solve_step(&forces, 1.0);
        assert!(u[0].abs() < 1e-6, "left end must be (nearly) fixed");
    }
    #[test]
    fn beam_right_end_has_positive_displacement() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        let u = beam.solve_step(&forces, 1.0);
        assert!(u[4] > 0.0, "tensile force → positive displacement");
    }
    #[test]
    fn creep_compliance_positive() {
        let mat = polymer();
        let j = creep_compliance_integral(&mat, 100.0, 1000);
        assert!(j > 0.0);
    }
    #[test]
    fn creep_compliance_increases_with_time() {
        let mat = polymer();
        let j1 = creep_compliance_integral(&mat, 1.0, 500);
        let j2 = creep_compliance_integral(&mat, 100.0, 500);
        assert!(j2 >= j1, "creep compliance should increase with time");
    }
    #[test]
    fn creep_compliance_at_zero_time() {
        let mat = polymer();
        let j = creep_compliance_integral(&mat, 0.0, 500);
        let expected = 1.0 / mat.instantaneous_modulus();
        assert!((j - expected).abs() / expected < 1e-10);
    }
    #[test]
    fn test_loss_modulus_positive() {
        let mat = polymer();
        let e_loss = mat.loss_modulus(1.0);
        assert!(e_loss > 0.0, "Loss modulus should be positive");
    }
    #[test]
    fn test_storage_modulus_between_bounds() {
        let mat = polymer();
        let e_storage = mat.storage_modulus(1.0);
        assert!(e_storage >= mat.e_inf, "Storage modulus should be >= E_inf");
        assert!(
            e_storage <= mat.instantaneous_modulus(),
            "Storage modulus should be <= E_0"
        );
    }
    #[test]
    fn test_storage_modulus_limits() {
        let mat = polymer();
        let e_low = mat.storage_modulus(1e-20);
        assert!((e_low - mat.e_inf).abs() / mat.e_inf < 1e-6);
        let e_high = mat.storage_modulus(1e20);
        let e0 = mat.instantaneous_modulus();
        assert!((e_high - e0).abs() / e0 < 1e-4);
    }
    #[test]
    fn test_loss_tangent_nonnegative() {
        let mat = polymer();
        for omega in [0.01, 0.1, 1.0, 10.0, 100.0] {
            let tan_d = mat.loss_tangent(omega);
            assert!(
                tan_d >= 0.0,
                "Loss tangent should be non-negative at omega={omega}"
            );
        }
    }
    #[test]
    fn test_complex_modulus_magnitude() {
        let mat = polymer();
        let mag = complex_modulus_magnitude(&mat, 1.0);
        let e_storage = mat.storage_modulus(1.0);
        assert!(mag >= e_storage, "|E*| >= E'");
    }
    #[test]
    fn test_wlf_shift_at_reference() {
        let wlf = WlfParameters::typical_polymer();
        let at = wlf.shift_factor(wlf.t_ref);
        assert!(
            (at - 1.0).abs() < 1e-10,
            "Shift factor at T_ref should be 1.0"
        );
    }
    #[test]
    fn test_wlf_shift_above_tref() {
        let wlf = WlfParameters::typical_polymer();
        let at = wlf.shift_factor(wlf.t_ref + 20.0);
        assert!(
            at < 1.0,
            "Shift factor above T_ref should be < 1 (faster relaxation)"
        );
    }
    #[test]
    fn test_wlf_shift_below_tref() {
        let wlf = WlfParameters::typical_polymer();
        let at = wlf.shift_factor(wlf.t_ref - 20.0);
        assert!(
            at > 1.0,
            "Shift factor below T_ref should be > 1 (slower relaxation)"
        );
    }
    #[test]
    fn test_wlf_reduced_time() {
        let wlf = WlfParameters::typical_polymer();
        let t_real = 10.0;
        let t_reduced = wlf.reduced_time(t_real, wlf.t_ref);
        assert!(
            (t_reduced - t_real).abs() < 1e-10,
            "At T_ref, reduced time = real time"
        );
    }
    #[test]
    fn test_wlf_relaxation_at_temp() {
        let wlf = WlfParameters::typical_polymer();
        let mat = polymer();
        let e_ref = wlf.relaxation_modulus_at_temp(&mat, 10.0, wlf.t_ref);
        let e_direct = mat.relaxation_modulus(10.0);
        assert!((e_ref - e_direct).abs() / e_direct < 1e-10);
    }
    #[test]
    fn test_arrhenius_shift_at_reference() {
        let arr = ArrheniusParameters {
            t_ref: 300.0,
            delta_h_over_r: 5000.0,
        };
        let at = arr.shift_factor(300.0);
        assert!((at - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_rubber_material() {
        let mat = ViscoelasticMaterial::rubber();
        assert_eq!(mat.num_terms(), 4);
        assert!(mat.instantaneous_modulus() > mat.e_inf);
    }
    #[test]
    fn test_normalized_relaxation() {
        let mat = polymer();
        assert!((mat.normalized_relaxation(0.0) - 1.0).abs() < 1e-12);
        let nr = mat.normalized_relaxation(1e12);
        assert!(nr < 1.0);
    }
    #[test]
    fn test_internal_variable_total_stress() {
        let mat = polymer();
        let mut iv = InternalVariable::new(mat.prony_terms.len());
        iv.update(0.01, &mat, 1.0);
        let total = iv.total_stress();
        assert!(total > 0.0);
    }
    #[test]
    fn test_element_strain_energy() {
        let mat = polymer();
        let mut elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        elem.stress_update(0.01, 1.0);
        let se = elem.strain_energy();
        assert!(se > 0.0);
    }
    #[test]
    fn test_element_internal_force_vector() {
        let mat = polymer();
        let mut elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        elem.stress_update(0.01, 1.0);
        let f = elem.internal_force_vector();
        assert!((f[0] + f[1]).abs() < 1e-10);
    }
    #[test]
    fn test_beam_internal_forces() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        beam.solve_step(&forces, 1.0);
        let f_int = beam.assemble_internal_forces();
        assert_eq!(f_int.len(), 5);
    }
    #[test]
    fn test_beam_total_strain_energy() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        beam.solve_step(&forces, 1.0);
        let se = beam.total_strain_energy();
        assert!(se > 0.0);
    }
    #[test]
    fn test_beam_max_stress() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        beam.solve_step(&forces, 1.0);
        let max_s = beam.max_stress();
        assert!(max_s > 0.0);
    }
    #[test]
    fn test_relaxation_spectrum() {
        let mat = polymer();
        let spectrum = relaxation_spectrum(&mat);
        assert_eq!(spectrum.len(), 3);
        for (tau, h) in &spectrum {
            assert!(*tau > 0.0);
            assert!(*h > 0.0);
        }
    }
    #[test]
    fn test_retardation_spectrum() {
        let mat = polymer();
        let spectrum = retardation_spectrum(&mat);
        assert_eq!(spectrum.len(), 3);
        for (tau, d) in &spectrum {
            assert!(*tau > 0.0);
            assert!(*d > 0.0);
        }
    }
    #[test]
    fn test_stress_update_with_tts() {
        let mat = polymer();
        let mut elem = ViscoelasticElement1D::new(1.0, 1e-4, mat);
        let sigma = elem.stress_update_with_tts(0.01, 1.0, 1.0);
        assert!(sigma.abs() > 0.0);
    }
    #[test]
    fn test_solve_with_temperature() {
        let mut beam = ViscoelasticBeam1D::new(4, 1.0, 1e-4, polymer());
        let forces = vec![0.0, 0.0, 0.0, 0.0, 1000.0];
        let wlf = WlfParameters::typical_polymer();
        let u = beam.solve_step_with_temperature(&forces, 1.0, &wlf, 25.0);
        assert_eq!(u.len(), 5);
        assert!(u[4] > 0.0);
    }
    #[test]
    fn test_prony_relaxation_at_t0() {
        let pr = PronyRelaxation {
            E_inf: 0.5e9,
            elements: vec![
                PronyElement {
                    E_i: 1.0e9,
                    tau_i: 1.0,
                },
                PronyElement {
                    E_i: 0.5e9,
                    tau_i: 10.0,
                },
            ],
        };
        let e0 = pr.relaxation_modulus(0.0);
        let expected = 0.5e9 + 1.0e9 + 0.5e9;
        assert!(
            (e0 - expected).abs() / expected < 1e-12,
            "E(0) should equal E_inf + sum(E_i), got {e0}"
        );
    }
    #[test]
    fn test_prony_relaxation_decays() {
        let pr = PronyRelaxation {
            E_inf: 1e6,
            elements: vec![PronyElement {
                E_i: 1e9,
                tau_i: 1.0,
            }],
        };
        let e0 = pr.relaxation_modulus(0.0);
        let e_late = pr.relaxation_modulus(1e10);
        assert!(e_late < e0, "relaxation modulus must decay");
        assert!((e_late - 1e6).abs() / 1e6 < 1e-4);
    }
    #[test]
    fn test_wlf_shift_at_tref_returns_zero() {
        let wlf = WLFShift {
            C1: 17.44,
            C2: 51.6,
            T_ref: 25.0,
        };
        let log_at = wlf.shift_factor(25.0);
        assert!(
            log_at.abs() < 1e-12,
            "log(a_T) at T_ref should be 0, got {log_at}"
        );
    }
    #[test]
    fn test_time_temperature_superposition_at_tref() {
        let wlf = WLFShift {
            C1: 17.44,
            C2: 51.6,
            T_ref: 25.0,
        };
        let t_red = time_temperature_superposition(10.0, 25.0, &wlf);
        assert!(
            (t_red - 10.0).abs() < 1e-10,
            "reduced time at T_ref should equal real time"
        );
    }
    #[test]
    fn test_creep_compliance_monotone() {
        let creep = Creep {
            J0: 1e-10,
            J_inf: 5e-10,
            tau: 100.0,
        };
        let mut prev = creep.compliance(0.0);
        for t in [1.0, 10.0, 100.0, 1000.0] {
            let j = creep.compliance(t);
            assert!(
                j >= prev - 1e-30,
                "compliance must be non-decreasing at t={t}"
            );
            prev = j;
        }
    }
    #[test]
    fn test_creep_compliance_at_zero() {
        let creep = Creep {
            J0: 2e-10,
            J_inf: 5e-10,
            tau: 50.0,
        };
        let j0 = creep.compliance(0.0);
        assert!((j0 - 2e-10).abs() < 1e-25, "J(0) should equal J0, got {j0}");
    }
    #[test]
    fn test_creep_compliance_at_infinity() {
        let creep = Creep {
            J0: 1e-10,
            J_inf: 5e-10,
            tau: 1.0,
        };
        let j_inf = creep.compliance(1e15);
        let expected = 1e-10 + 5e-10;
        assert!((j_inf - expected).abs() / expected < 1e-6);
    }
    #[test]
    fn test_prony_fem_element_update_changes_state() {
        let prony = PronyRelaxation {
            E_inf: 1e9,
            elements: vec![PronyElement {
                E_i: 0.5e9,
                tau_i: 1.0,
            }],
        };
        let mut elem = PronyFemElement::new([0, 1, 2, 3], prony, 0.1);
        let strain_inc = vec![0.001];
        let stress_before = elem.state.stress.clone();
        let stress_after = elem.update_stresses(&strain_inc);
        assert_ne!(
            stress_before, stress_after,
            "stress should change after strain increment"
        );
    }
    #[test]
    fn test_prony_fem_element_internal_vars_updated() {
        let prony = PronyRelaxation {
            E_inf: 1e9,
            elements: vec![
                PronyElement {
                    E_i: 0.5e9,
                    tau_i: 1.0,
                },
                PronyElement {
                    E_i: 0.3e9,
                    tau_i: 10.0,
                },
            ],
        };
        let mut elem = PronyFemElement::new([0, 1, 2, 3], prony, 0.1);
        let strain_inc = vec![0.001];
        elem.update_stresses(&strain_inc);
        let any_nonzero = elem
            .state
            .internal_vars
            .iter()
            .any(|iv| iv.iter().any(|&v| v.abs() > 1e-30));
        assert!(
            any_nonzero,
            "internal variables should be updated after strain increment"
        );
    }
}
/// Evaluate the hereditary (convolution) integral for viscoelastic stress.
///
/// sigma(t) = integral_0^t E(t-s) * d(epsilon)/ds  ds
///
/// Uses the midpoint rule with uniform time steps.
///
/// # Arguments
/// * `material`      - viscoelastic material
/// * `strain_history` - strain at each time step \[epsilon_0, ..., epsilon_n\]
/// * `dt`            - time step (s)
pub fn hereditary_stress_integral(
    material: &ViscoelasticMaterial,
    strain_history: &[f64],
    dt: f64,
) -> f64 {
    let n = strain_history.len();
    if n < 2 {
        return 0.0;
    }
    let mut sigma = 0.0;
    let t_total = (n - 1) as f64 * dt;
    for k in 0..(n - 1) {
        let d_eps = strain_history[k + 1] - strain_history[k];
        let t_k = (k as f64 + 0.5) * dt;
        let t_minus_t_k = t_total - t_k;
        let e_rel = material.relaxation_modulus(t_minus_t_k);
        sigma += e_rel * d_eps;
    }
    sigma
}
/// Perform a frequency sweep over a logarithmic range.
///
/// Returns a Vec of `FrequencySweepPoint` for each frequency.
///
/// # Arguments
/// * `material`   - viscoelastic material
/// * `omega_min`  - minimum angular frequency (rad/s)
/// * `omega_max`  - maximum angular frequency (rad/s)
/// * `n_points`   - number of frequency points
pub fn frequency_sweep(
    material: &ViscoelasticMaterial,
    omega_min: f64,
    omega_max: f64,
    n_points: usize,
) -> Vec<FrequencySweepPoint> {
    if n_points == 0 || omega_min <= 0.0 || omega_max <= omega_min {
        return Vec::new();
    }
    let log_min = omega_min.ln();
    let log_max = omega_max.ln();
    (0..n_points)
        .map(|i| {
            let log_omega =
                log_min + (i as f64 / (n_points - 1).max(1) as f64) * (log_max - log_min);
            let omega = log_omega.exp();
            let e_prime = material.storage_modulus(omega);
            let e_double_prime = material.loss_modulus(omega);
            let tan_delta = material.loss_tangent(omega);
            let magnitude = (e_prime * e_prime + e_double_prime * e_double_prime).sqrt();
            FrequencySweepPoint {
                omega,
                storage: e_prime,
                loss: e_double_prime,
                tan_delta,
                magnitude,
            }
        })
        .collect()
}
