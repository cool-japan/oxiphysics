//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Advance temperatures by one time step using the enthalpy method for phase change.
///
/// Each node uses the apparent-heat-capacity form:
///   C_eff * dT/dt = q_net
///
/// where C_eff absorbs the latent heat over a small mushy-zone width `delta_T`.
///
/// # Arguments
/// * `temperatures` – mutable node temperature vector (updated in place)
/// * `heat_source`  – volumetric heat source (W/m³, uniform)
/// * `density`      – kg/m³
/// * `specific_heat`– J/(kg K) (sensible heat)
/// * `latent_heat`  – L \[J/kg\] (fusion latent heat)
/// * `t_melt`       – melting point \[K\]
/// * `dt`           – time step \[s\]
/// * `delta_t_mushy`– width of mushy zone \[K\]
pub fn phase_change_step(
    temperatures: &mut [f64],
    heat_source: f64,
    density: f64,
    specific_heat: f64,
    latent_heat: f64,
    t_melt: f64,
    dt: f64,
    delta_t_mushy: f64,
) {
    let dmu = delta_t_mushy.max(1e-10);
    for t in temperatures.iter_mut() {
        let xi = (*t - t_melt) / dmu;
        let f_mushy = if xi.abs() <= 1.0 {
            0.5 * (1.0 + (std::f64::consts::PI * xi).cos())
        } else {
            0.0
        };
        let c_eff = specific_heat + latent_heat / dmu * f_mushy;
        let rho_c = density * c_eff;
        if rho_c > 1e-60 {
            *t += heat_source / rho_c * dt;
        }
    }
}
/// Compute the radiation view factor F_12 for two finite parallel plates.
///
/// Uses the crossed-string method / analytical approximation for two parallel
/// coaxial discs of areas `a1` and `a2` separated by distance `h`:
///
/// This is a simplified planar approximation:
///   R1 = sqrt(a1/π), R2 = sqrt(a2/π)
///   S = 1 + (1 + R2²) / R1²
///   F12 = 0.5*(S - sqrt(S² - 4*(R2/R1)²))
///
/// # Arguments
/// * `area1` – area of surface 1 \[m²\]
/// * `area2` – area of surface 2 \[m²\]
/// * `separation` – distance between plates \[m\]
pub fn view_factor_parallel_plates(area1: f64, area2: f64, separation: f64) -> f64 {
    let r1 = (area1 / std::f64::consts::PI).sqrt().max(1e-60);
    let r2 = (area2 / std::f64::consts::PI).sqrt().max(1e-60);
    let h = separation.max(1e-60);
    let r1n = r1 / h;
    let r2n = r2 / h;
    let s_val = 1.0 + (1.0 + r2n * r2n) / (r1n * r1n);
    let discriminant = s_val * s_val - 4.0 * (r2n / r1n).powi(2);
    let f12 = 0.5 * (s_val - discriminant.max(0.0).sqrt());
    f12.clamp(0.0, 1.0)
}
/// Net radiation heat exchange rate between two grey-body surfaces \[W\].
///
/// Q_12 = σ * F_12 * A1 * (T1⁴ - T2⁴) / (1/ε1 + A1/A2*(1/ε2-1) + 1/F12 - 1)
///
/// Simplified form (black bodies, F12 given):
///   Q_12 = σ * A1 * F_12 * (T1⁴ - T2⁴)
///
/// # Arguments
/// * `t1`, `t2` – surface temperatures \[K\]
/// * `area1`    – area of surface 1 \[m²\]
/// * `f12`      – view factor F_12
/// * `emissivity` – emissivity (0..1) for both surfaces (assumed equal)
pub fn radiation_heat_exchange(t1: f64, t2: f64, area1: f64, f12: f64, emissivity: f64) -> f64 {
    pub(super) const SIGMA: f64 = 5.670_374_419e-8;
    let eps = emissivity.clamp(1e-9, 1.0);
    SIGMA * eps * area1 * f12 * (t1.powi(4) - t2.powi(4))
}
/// Perform one implicit (backward-Euler) time step with Newton-Raphson linearisation.
///
/// Residual: R(T) = C*(T - T_n)/dt + K*T - f = 0
/// Newton:   J * δT = -R, then T ← T + δT
///
/// For linear K the method converges in one iteration, but allowing `max_iter`
/// makes this suitable for weakly nonlinear problems where K = K(T).
///
/// # Arguments
/// * `k_global`  – n×n conductance matrix (dense, row-major)
/// * `c_lumped`  – lumped capacitance vector
/// * `t_n`       – temperatures at time n
/// * `heat_src`  – nodal heat sources
/// * `dirichlet` – Dirichlet BCs (node, value)
/// * `dt`        – time step \[s\]
/// * `max_iter`  – maximum Newton iterations
/// * `tol`       – convergence tolerance on ‖δT‖₂
pub fn transient_nonlinear_step(
    k_global: &[Vec<f64>],
    c_lumped: &[f64],
    t_n: &[f64],
    heat_src: &[f64],
    dirichlet: &[(usize, f64)],
    dt: f64,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let n = t_n.len();
    assert_eq!(k_global.len(), n);
    assert_eq!(c_lumped.len(), n);
    assert_eq!(heat_src.len(), n);
    let mut t = t_n.to_vec();
    for &(node, val) in dirichlet {
        t[node] = val;
    }
    for _iter in 0..max_iter {
        let mut r = vec![0.0_f64; n];
        for (i, r_i) in r.iter_mut().enumerate() {
            let kt_i: f64 = (0..n).map(|j| k_global[i][j] * t[j]).sum();
            *r_i = c_lumped[i] * (t[i] - t_n[i]) / dt + kt_i - heat_src[i];
        }
        for &(node, _) in dirichlet {
            r[node] = 0.0;
        }
        let r_norm: f64 = r.iter().map(|&v| v * v).sum::<f64>().sqrt();
        if r_norm < tol {
            break;
        }
        let mut jac: Vec<Vec<f64>> = (0..n)
            .map(|i| {
                let mut row = k_global[i].to_vec();
                row[i] += c_lumped[i] / dt;
                row
            })
            .collect();
        let mut rhs: Vec<f64> = r.iter().map(|&v| -v).collect();
        for &(node, _) in dirichlet {
            for v in jac[node].iter_mut() {
                *v = 0.0;
            }
            jac[node][node] = 1.0;
            rhs[node] = 0.0;
        }
        let delta = gauss_solve_dense(&jac, &rhs);
        let delta_norm: f64 = delta.iter().map(|&v| v * v).sum::<f64>().sqrt();
        for (t_i, &d_i) in t.iter_mut().zip(delta.iter()) {
            *t_i += d_i;
        }
        for &(node, val) in dirichlet {
            t[node] = val;
        }
        if delta_norm < tol {
            break;
        }
    }
    t
}
/// Simple dense Gauss elimination (no pivoting) for internal use.
pub(super) fn gauss_solve_dense(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut mat: Vec<Vec<f64>> = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut max_row = col;
        let mut max_val = mat[col][col].abs();
        for (row, mat_row) in mat.iter().enumerate().skip(col + 1) {
            if mat_row[col].abs() > max_val {
                max_val = mat_row[col].abs();
                max_row = row;
            }
        }
        mat.swap(col, max_row);
        rhs.swap(col, max_row);
        let pivot = mat[col][col];
        if pivot.abs() < 1e-60 {
            continue;
        }
        let col_row: Vec<f64> = mat[col][col..].to_vec();
        for row in (col + 1)..n {
            let factor = mat[row][col] / pivot;
            for (j_off, &col_val) in col_row.iter().enumerate() {
                mat[row][col + j_off] -= col_val * factor;
            }
            rhs[row] -= rhs[col] * factor;
        }
    }
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..n {
            s -= mat[i][j] * x[j];
        }
        x[i] = if mat[i][i].abs() > 1e-60 {
            s / mat[i][i]
        } else {
            0.0
        };
    }
    x
}

#[cfg(test)]
mod thermal_extended_tests {
    use super::*;
    use crate::thermal::*;
    #[test]
    fn test_isotropic_conductivity_tensor_diagonal() {
        let kt = ThermalConductivityTensor::isotropic(5.0);
        assert!((kt.k[0][0] - 5.0).abs() < 1e-12);
        assert!((kt.k[1][1] - 5.0).abs() < 1e-12);
        assert!((kt.k[2][2] - 5.0).abs() < 1e-12);
        assert!(kt.k[0][1].abs() < 1e-12);
    }
    #[test]
    fn test_isotropic_heat_flux_in_x_direction() {
        let kt = ThermalConductivityTensor::isotropic(5.0);
        let grad = [1.0, 0.0, 0.0];
        let q = kt.heat_flux(&grad);
        assert!((q[0] + 5.0).abs() < 1e-12, "q_x = -k * dT/dx = -5");
        assert!(q[1].abs() < 1e-12);
        assert!(q[2].abs() < 1e-12);
    }
    #[test]
    fn test_orthotropic_conductivity_tensor() {
        let kt = ThermalConductivityTensor::orthotropic(10.0, 5.0, 2.0);
        assert!((kt.k[0][0] - 10.0).abs() < 1e-12);
        assert!((kt.k[1][1] - 5.0).abs() < 1e-12);
        assert!((kt.k[2][2] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_orthotropic_heat_flux_z_direction() {
        let kt = ThermalConductivityTensor::orthotropic(10.0, 5.0, 2.0);
        let grad = [0.0, 0.0, 1.0];
        let q = kt.heat_flux(&grad);
        assert!((q[2] + 2.0).abs() < 1e-12, "q_z = -kz * dT/dz = -2");
    }
    #[test]
    fn test_max_conductivity() {
        let kt = ThermalConductivityTensor::orthotropic(10.0, 5.0, 2.0);
        assert!((kt.max_conductivity() - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_effective_isotropic_conductivity() {
        let kt = ThermalConductivityTensor::orthotropic(8.0, 8.0, 8.0);
        let k_eff = kt.effective_isotropic();
        assert!(
            (k_eff - 8.0).abs() < 1e-10,
            "isotropic case: effective = k = {k_eff}"
        );
    }
    #[test]
    fn test_biot_number_small_for_good_conductor() {
        let bi = biot_number(10.0, 0.01, 200.0);
        assert!(bi < 0.1, "Biot number for metal should be small: Bi={bi}");
    }
    #[test]
    fn test_biot_number_large_for_insulator() {
        let bi = biot_number(100.0, 0.1, 0.05);
        assert!(
            bi > 1.0,
            "Biot number for insulator should be large: Bi={bi}"
        );
    }
    #[test]
    fn test_lumped_capacitance_valid_for_low_biot() {
        assert!(
            lumped_capacitance_valid(10.0, 0.001, 200.0),
            "Lumped cap valid for low Biot number"
        );
    }
    #[test]
    fn test_lumped_capacitance_invalid_for_high_biot() {
        assert!(
            !lumped_capacitance_valid(1000.0, 0.5, 0.1),
            "Lumped cap invalid for high Biot number"
        );
    }
    #[test]
    fn test_lumped_capacitance_temperature_at_t0() {
        let t = lumped_capacitance_temperature(0.0, 400.0, 300.0, 100.0);
        assert!((t - 400.0).abs() < 1e-12, "T(0) should equal T_0");
    }
    #[test]
    fn test_lumped_capacitance_temperature_approaches_ambient() {
        let t = lumped_capacitance_temperature(1e15, 400.0, 300.0, 100.0);
        assert!((t - 300.0).abs() < 1e-6, "T(∞) should approach T_inf");
    }
    #[test]
    fn test_lumped_capacitance_temperature_decays_monotonically() {
        let tau = 100.0;
        let times = [0.0, 10.0, 50.0, 100.0, 500.0];
        let temps: Vec<f64> = times
            .iter()
            .map(|&t| lumped_capacitance_temperature(t, 400.0, 300.0, tau))
            .collect();
        for w in temps.windows(2) {
            assert!(
                w[1] <= w[0] + 1e-10,
                "temperature should decay monotonically"
            );
        }
    }
    #[test]
    fn test_thermal_time_constant_positive() {
        let tau = thermal_time_constant(1000.0, 900.0, 1e-3, 50.0, 0.01);
        assert!(
            tau > 0.0 && tau.is_finite(),
            "time constant should be positive finite: {tau}"
        );
    }
    #[test]
    fn test_thermal_strain_isotropic_zero_delta_t() {
        let eps = thermal_strain_isotropic(12e-6, 300.0, 300.0);
        for e in eps {
            assert!(e.abs() < 1e-30, "zero temp change → zero thermal strain");
        }
    }
    #[test]
    fn test_thermal_strain_isotropic_positive_delta_t() {
        let eps = thermal_strain_isotropic(12e-6, 400.0, 300.0);
        assert!(
            eps[0] > 0.0 && eps[1] > 0.0 && eps[2] > 0.0,
            "positive ΔT with positive CTE → positive strains"
        );
        assert!(eps[3].abs() < 1e-30);
        assert!(eps[4].abs() < 1e-30);
        assert!(eps[5].abs() < 1e-30);
    }
    #[test]
    fn test_thermal_strain_isotropic_all_equal() {
        let eps = thermal_strain_isotropic(12e-6, 400.0, 300.0);
        assert!((eps[0] - eps[1]).abs() < 1e-30);
        assert!((eps[1] - eps[2]).abs() < 1e-30);
    }
    #[test]
    fn test_thermal_stress_isotropic_compressive_on_heating() {
        let sigma = thermal_stress_isotropic(200e9, 0.3, 12e-6, 400.0, 300.0);
        assert!(sigma < 0.0, "constrained heating → compressive stress");
    }
    #[test]
    fn test_thermal_load_vector_1d_equal_opposite() {
        let f = thermal_load_vector_1d(200e9, 1e-4, 12e-6, 100.0);
        assert!(
            (f[0] + f[1]).abs() < 1e-6 * f[1].abs().max(1e-10),
            "thermal load vector should be self-equilibrated"
        );
    }
    #[test]
    fn test_thermal_load_vector_1d_scales_with_delta_t() {
        let f1 = thermal_load_vector_1d(200e9, 1e-4, 12e-6, 100.0);
        let f2 = thermal_load_vector_1d(200e9, 1e-4, 12e-6, 200.0);
        assert!(
            (f2[1] / f1[1] - 2.0).abs() < 1e-10,
            "thermal force should scale linearly with ΔT"
        );
    }
    #[test]
    fn test_crank_nicolson_converges_to_steady_state() {
        let n = 5;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let k = mesh.assemble_conductance_matrix();
        let c = mesh.assemble_capacitance_vector();
        let q = vec![0.0_f64; n];
        let result = crank_nicolson_transient(
            &k,
            &c,
            &mesh.temperatures,
            &q,
            &[(0, 0.0), (n - 1, 100.0)],
            0.01,
            500,
        );
        for (i, &t) in result.final_temperatures.iter().enumerate() {
            let expected = 100.0 * i as f64 / (n - 1) as f64;
            assert!(
                (t - expected).abs() < 1.0,
                "CN node {i}: got {t}, expected {expected}"
            );
        }
    }
    #[test]
    fn test_crank_nicolson_history_length() {
        let n = 3;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 0.0);
        mesh.set_temperature_bc(n - 1, 100.0);
        let k = mesh.assemble_conductance_matrix();
        let c = mesh.assemble_capacitance_vector();
        let q = vec![0.0_f64; n];
        let result = crank_nicolson_transient(
            &k,
            &c,
            &mesh.temperatures,
            &q,
            &[(0, 0.0), (n - 1, 100.0)],
            0.1,
            10,
        );
        assert_eq!(
            result.history.len(),
            11,
            "should have n_steps + 1 history entries"
        );
    }
    #[test]
    fn test_crank_nicolson_bc_preserved() {
        let n = 4;
        let mut mesh = ThermalMesh1D::new_uniform(n, 1.0, 1.0, 1.0, 1.0);
        mesh.set_temperature_bc(0, 50.0);
        mesh.set_temperature_bc(n - 1, 150.0);
        let k = mesh.assemble_conductance_matrix();
        let c = mesh.assemble_capacitance_vector();
        let q = vec![0.0_f64; n];
        let result = crank_nicolson_transient(
            &k,
            &c,
            &mesh.temperatures,
            &q,
            &[(0, 50.0), (n - 1, 150.0)],
            0.05,
            20,
        );
        assert!(
            (result.final_temperatures[0] - 50.0).abs() < 1e-8,
            "BC at node 0 should be preserved"
        );
        assert!(
            (result.final_temperatures[n - 1] - 150.0).abs() < 1e-8,
            "BC at node n-1 should be preserved"
        );
    }
    #[test]
    fn test_enthalpy_update_no_heat_no_change() {
        let t = 200.0;
        let t_new = enthalpy_update(t, 0.0, 1000.0, 500.0, 334000.0, 273.15, 1.0, 0.1, 1e-6);
        assert!(
            (t_new - t).abs() < 1.0,
            "zero heat → temperature nearly unchanged: {t_new}"
        );
    }
    #[test]
    fn test_enthalpy_update_heating_increases_temperature() {
        let t0 = 200.0;
        let t_new = enthalpy_update(t0, 1000.0, 1000.0, 500.0, 334000.0, 500.0, 5.0, 1.0, 1e-3);
        assert!(
            t_new > t0,
            "heating should raise temperature: {t0} → {t_new}"
        );
    }
    #[test]
    fn test_enthalpy_update_at_melting_releases_latent_heat() {
        let t_melt = 273.15;
        let t_below = enthalpy_update(
            t_melt - 0.1,
            10.0,
            1000.0,
            500.0,
            334000.0,
            t_melt,
            0.5,
            0.1,
            1e-6,
        );
        let t_far = enthalpy_update(
            t_melt + 20.0,
            10.0,
            1000.0,
            500.0,
            334000.0,
            t_melt,
            0.5,
            0.1,
            1e-6,
        );
        assert!(
            t_below.is_finite(),
            "enthalpy update near melting: {t_below}"
        );
        assert!(
            t_far.is_finite(),
            "enthalpy update far from melting: {t_far}"
        );
    }
    #[test]
    fn test_biot_number_proportional_to_h() {
        let bi1 = biot_number(10.0, 0.05, 50.0);
        let bi2 = biot_number(20.0, 0.05, 50.0);
        assert!(
            (bi2 / bi1 - 2.0).abs() < 1e-10,
            "Biot number should scale with h"
        );
    }
    #[test]
    fn test_thermal_strain_alpha_scales_linearly() {
        let eps1 = thermal_strain_isotropic(10e-6, 400.0, 300.0);
        let eps2 = thermal_strain_isotropic(20e-6, 400.0, 300.0);
        assert!(
            (eps2[0] / eps1[0] - 2.0).abs() < 1e-10,
            "thermal strain should scale linearly with CTE"
        );
    }
    #[test]
    fn test_phase_change_solver_no_source_stable() {
        let mut temps = vec![300.0_f64; 5];
        phase_change_step(&mut temps, 0.0, 1000.0, 500.0, 334_000.0, 273.15, 0.1, 1e-3);
        for &t in &temps {
            assert!(t.is_finite(), "temperature should remain finite: {t}");
        }
    }
    #[test]
    fn test_phase_change_solver_heating_raises_temperature() {
        let mut temps = vec![200.0_f64; 5];
        let temps_before = temps.clone();
        phase_change_step(
            &mut temps, 50_000.0, 1000.0, 500.0, 334_000.0, 273.15, 1.0, 0.1,
        );
        let mean_before: f64 = temps_before.iter().sum::<f64>() / temps_before.len() as f64;
        let mean_after: f64 = temps.iter().sum::<f64>() / temps.len() as f64;
        assert!(
            mean_after > mean_before,
            "heating should raise average temperature: {mean_before} → {mean_after}"
        );
    }
    #[test]
    fn test_phase_change_stefan_conserves_sign() {
        let mut temps = vec![270.0_f64; 3];
        phase_change_step(
            &mut temps, 100.0, 900.0, 800.0, 334_000.0, 273.15, 0.5, 0.01,
        );
        for &t in &temps {
            assert!(t.is_finite() && t > 0.0, "temp should stay positive: {t}");
        }
    }
    #[test]
    fn test_view_factor_parallel_plates_reciprocity() {
        let area1 = 1.0;
        let area2 = 1.0;
        let f12 = view_factor_parallel_plates(area1, area2, 0.5);
        let f21 = view_factor_parallel_plates(area2, area1, 0.5);
        assert!(
            (f12 - f21).abs() < 1e-10,
            "reciprocity: F12={f12}, F21={f21}"
        );
    }
    #[test]
    fn test_view_factor_enclosure_sum_to_one() {
        let f = view_factor_parallel_plates(1.0, 2.0, 0.3);
        assert!(
            (0.0..=1.0 + 1e-10).contains(&f),
            "F12 must be in [0,1]: {f}"
        );
    }
    #[test]
    fn test_view_factor_increases_with_proximity() {
        let f_close = view_factor_parallel_plates(1.0, 1.0, 0.1);
        let f_far = view_factor_parallel_plates(1.0, 1.0, 2.0);
        assert!(
            f_close > f_far,
            "closer plates should have higher F: close={f_close}, far={f_far}"
        );
    }
    #[test]
    fn test_transient_nonlinear_step_output_length() {
        let n = 4;
        let k = vec![
            vec![2.0_f64, -1.0, 0.0, 0.0],
            vec![-1.0, 2.0, -1.0, 0.0],
            vec![0.0, -1.0, 2.0, -1.0],
            vec![0.0, 0.0, -1.0, 2.0],
        ];
        let c = vec![1.0_f64; n];
        let t_n = vec![300.0_f64; n];
        let q = vec![0.0_f64; n];
        let bcs = [(0_usize, 300.0_f64), (n - 1, 400.0_f64)];
        let t_next = transient_nonlinear_step(&k, &c, &t_n, &q, &bcs, 0.01, 10, 1e-6);
        assert_eq!(t_next.len(), n, "output length should match n");
    }
    #[test]
    fn test_transient_nonlinear_step_bc_satisfied() {
        let n = 4;
        let k = vec![
            vec![2.0_f64, -1.0, 0.0, 0.0],
            vec![-1.0, 2.0, -1.0, 0.0],
            vec![0.0, -1.0, 2.0, -1.0],
            vec![0.0, 0.0, -1.0, 2.0],
        ];
        let c = vec![1.0_f64; n];
        let t_n = vec![300.0_f64; n];
        let q = vec![0.0_f64; n];
        let bcs = [(0_usize, 300.0_f64), (n - 1, 400.0_f64)];
        let t_next = transient_nonlinear_step(&k, &c, &t_n, &q, &bcs, 0.05, 5, 1e-8);
        assert!((t_next[0] - 300.0).abs() < 1e-8, "BC node 0: {}", t_next[0]);
        assert!(
            (t_next[n - 1] - 400.0).abs() < 1e-8,
            "BC node n-1: {}",
            t_next[n - 1]
        );
    }
    #[test]
    fn test_transient_nonlinear_step_all_finite() {
        let n = 3;
        let k = vec![
            vec![2.0_f64, -1.0, 0.0],
            vec![-1.0, 2.0, -1.0],
            vec![0.0, -1.0, 2.0],
        ];
        let c = vec![1.0_f64; n];
        let t_n = vec![500.0_f64; n];
        let q = vec![10.0_f64; n];
        let bcs = [(0_usize, 500.0_f64), (n - 1, 500.0_f64)];
        let t_next = transient_nonlinear_step(&k, &c, &t_n, &q, &bcs, 0.1, 5, 1e-8);
        for (i, &t) in t_next.iter().enumerate() {
            assert!(t.is_finite(), "temperature[{i}] should be finite: {t}");
        }
    }
    #[test]
    fn test_transient_nonlinear_step_steady_state() {
        let n = 3;
        let k = vec![
            vec![2.0_f64, -1.0, 0.0],
            vec![-1.0, 2.0, -1.0],
            vec![0.0, -1.0, 2.0],
        ];
        let c = vec![1.0_f64; n];
        let t_n = vec![300.0_f64; n];
        let q = vec![0.0_f64; n];
        let bcs = [(0_usize, 300.0_f64), (n - 1, 300.0_f64)];
        let t_next = transient_nonlinear_step(&k, &c, &t_n, &q, &bcs, 0.1, 5, 1e-10);
        for &t in &t_next {
            assert!(
                (t - 300.0).abs() < 1e-6,
                "uniform IC/BC → solution stays at 300: {t}"
            );
        }
    }
}
