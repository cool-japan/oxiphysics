//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Trait for equations of state that relate density/energy to pressure.
pub trait EquationOfState {
    /// Compute pressure from density (Pa).
    fn pressure(&self, density: f64) -> f64;
    /// Compute sound speed from density (m/s).
    fn sound_speed(&self, density: f64) -> f64;
    /// Compute density from pressure (inverse of pressure, if invertible).
    fn density_from_pressure(&self, pressure: f64) -> f64;
}
/// Extended EOS trait that includes internal energy.
pub trait EosWithEnergy: EquationOfState {
    /// Compute pressure from density and specific internal energy (J/kg).
    fn pressure_energy(&self, density: f64, energy: f64) -> f64;
    /// Compute internal energy from density and temperature (K).
    fn internal_energy(&self, density: f64, temperature: f64) -> f64;
}
pub(super) const R_UNIVERSAL: f64 = 8.314_462_618;
/// Speed of sound in an ideal gas.
///
/// c = sqrt(γ * P / ρ)
///
/// # Arguments
/// * `gamma` - Adiabatic index γ = Cp/Cv.
/// * `p`     - Pressure (Pa).
/// * `rho`   - Density (kg/m³).
pub fn speed_of_sound_ideal(gamma: f64, p: f64, rho: f64) -> f64 {
    if rho <= 0.0 || p <= 0.0 {
        return 0.0;
    }
    (gamma * p / rho).sqrt()
}
/// Fit a linear Hugoniot relation Us = c₀ + s·Up from experimental data.
///
/// Uses ordinary least squares on Up (particle velocity) vs Us (shock velocity).
///
/// Returns (c₀, s) where Us = c₀ + s·Up.
pub fn fit_hugoniot_linear(up: &[f64], us: &[f64]) -> (f64, f64) {
    assert_eq!(up.len(), us.len(), "up and us must have same length");
    let n = up.len() as f64;
    if up.is_empty() {
        return (0.0, 0.0);
    }
    let sum_x: f64 = up.iter().sum();
    let sum_y: f64 = us.iter().sum();
    let sum_xx: f64 = up.iter().map(|x| x * x).sum();
    let sum_xy: f64 = up.iter().zip(us.iter()).map(|(x, y)| x * y).sum();
    let denom = n * sum_xx - sum_x * sum_x;
    if denom.abs() < f64::EPSILON {
        return (sum_y / n, 0.0);
    }
    let s = (n * sum_xy - sum_x * sum_y) / denom;
    let c0 = (sum_y - s * sum_x) / n;
    (c0, s)
}
/// Fit a Birch-Murnaghan 3rd-order EOS to experimental (V, P) data.
///
/// Uses iterative Gauss-Newton least squares to find (K₀, K₀').
/// The reference volume V₀ must be provided (known from 0-pressure experiment).
///
/// Returns (K₀, K₀') in the same units as the input pressures.
pub fn fit_birch_murnaghan_3rd(
    v0: f64,
    volumes: &[f64],
    pressures: &[f64],
    max_iter: usize,
) -> (f64, f64) {
    assert_eq!(volumes.len(), pressures.len());
    if volumes.is_empty() {
        return (100.0e9, 4.0);
    }
    let mut k0 = 100.0e9_f64;
    let mut k0p = 4.0_f64;
    let bm3_p = |v: f64, k0_: f64, k0p_: f64| -> f64 {
        let x = v0 / v;
        let x73 = x.powf(7.0 / 3.0);
        let x53 = x.powf(5.0 / 3.0);
        let x23 = x.powf(2.0 / 3.0);
        1.5 * k0_ * (x73 - x53) * (1.0 + 0.75 * (k0p_ - 4.0) * (x23 - 1.0))
    };
    let lambda = 1e-4;
    for _ in 0..max_iter {
        let mut jt_j = [[0.0f64; 2]; 2];
        let mut jt_r = [0.0f64; 2];
        for (&v, &p_data) in volumes.iter().zip(pressures.iter()) {
            let p_model = bm3_p(v, k0, k0p);
            let residual = p_data - p_model;
            let dk0 = 1.0e5_f64;
            let dk0p = 0.001_f64;
            let dp_dk0 = (bm3_p(v, k0 + dk0, k0p) - bm3_p(v, k0 - dk0, k0p)) / (2.0 * dk0);
            let dp_dk0p = (bm3_p(v, k0, k0p + dk0p) - bm3_p(v, k0, k0p - dk0p)) / (2.0 * dk0p);
            let j = [dp_dk0, dp_dk0p];
            for i in 0..2 {
                for k in 0..2 {
                    jt_j[i][k] += j[i] * j[k];
                }
                jt_r[i] += j[i] * residual;
            }
        }
        jt_j[0][0] += lambda * jt_j[0][0].abs().max(1.0);
        jt_j[1][1] += lambda * jt_j[1][1].abs().max(1.0);
        let det = jt_j[0][0] * jt_j[1][1] - jt_j[0][1] * jt_j[1][0];
        if det.abs() < 1e-100 {
            break;
        }
        let dk0_step = (jt_j[1][1] * jt_r[0] - jt_j[0][1] * jt_r[1]) / det;
        let dk0p_step = (jt_j[0][0] * jt_r[1] - jt_j[1][0] * jt_r[0]) / det;
        k0 += dk0_step;
        k0p += dk0p_step;
        k0 = k0.max(1.0e6);
        k0p = k0p.max(1.0);
    }
    (k0, k0p)
}
/// Compute Rankine-Hugoniot shock jump residuals.
///
/// Given pre-shock (ρ₀, p₀, u₀=0) and post-shock (ρ, p, u_p=particle velocity)
/// and shock speed U_s, returns the residuals of the three conservation equations:
///   \[mass, momentum, energy\]
///
/// Perfect conservation → all residuals = 0.
pub fn rankine_hugoniot_residuals(
    rho0: f64,
    p0: f64,
    rho1: f64,
    p1: f64,
    u_p: f64,
    u_s: f64,
) -> [f64; 3] {
    let mass_res = rho0 * u_s - rho1 * (u_s - u_p);
    let mom_res = p0 + rho0 * u_s * u_s - (p1 + rho1 * (u_s - u_p).powi(2));
    let e0 = p0 / (rho0 * 0.4);
    let e1 = p1 / (rho1 * 0.4);
    let energy_res =
        e0 + p0 / rho0 + 0.5 * u_s.powi(2) - (e1 + p1 / rho1 + 0.5 * (u_s - u_p).powi(2));
    [mass_res, mom_res, energy_res]
}
/// Hugoniot equation of state: p_H(V) from Us-Up linear relation.
///
/// Using mass and momentum conservation:
///   p_H = ρ₀ Us Up = ρ₀ (c₀ + s·Up) Up
///   V = V₀ (1 - Up/Us)
/// This gives p_H as a function of V via parametric elimination.
///
/// Returns the Hugoniot pressure at specific volume V.
pub fn hugoniot_pressure_volume(rho0: f64, c0: f64, s: f64, v: f64) -> f64 {
    let v0 = 1.0 / rho0;
    let mu = v0 / v - 1.0;
    if mu <= 0.0 {
        return 0.0;
    }
    let denom = (1.0 - s * mu).powi(2);
    if denom.abs() < f64::EPSILON {
        return 0.0;
    }
    rho0 * c0 * c0 * mu / denom
}
/// Statistical EOS uncertainty: bootstrap error estimate for Hugoniot fit.
///
/// Given Up/Us data, returns the standard error on c₀ and s from OLS.
pub fn hugoniot_fit_std_error(up: &[f64], us: &[f64]) -> (f64, f64) {
    assert_eq!(up.len(), us.len());
    let n = up.len();
    if n < 3 {
        return (0.0, 0.0);
    }
    let (c0, s) = fit_hugoniot_linear(up, us);
    let sse: f64 = up
        .iter()
        .zip(us.iter())
        .map(|(x, y)| (y - c0 - s * x).powi(2))
        .sum();
    let sigma2 = sse / (n as f64 - 2.0);
    let n_f = n as f64;
    let x_bar: f64 = up.iter().sum::<f64>() / n_f;
    let sxx: f64 = up.iter().map(|x| (x - x_bar).powi(2)).sum();
    if sxx < f64::EPSILON {
        return (0.0, 0.0);
    }
    let se_s = (sigma2 / sxx).sqrt();
    let se_c0 = (sigma2 * (1.0 / n_f + x_bar * x_bar / sxx)).sqrt();
    (se_c0, se_s)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::IdealGasEos;
    use crate::PolynomialEos;
    use crate::StiffenedGasEos;
    use crate::TaitEos;
    use crate::TillotsonEos;
    use crate::VanDerWaalsEos;
    use crate::eos::BirchMurnaghan3Eos;
    use crate::eos::GruneisenParameter;
    use crate::eos::JwlEos;
    use crate::eos::MieGruneisenEos;
    use crate::eos::MurnaghanEos;
    use crate::eos::NobleAbelEos;
    use crate::eos::PengRobinson;
    use crate::eos::PvtEos;
    use crate::eos::RedlichKwong;
    use crate::eos::ShockHugoniot;
    use crate::eos::SoaveRedlichKwong;
    use crate::eos::TaitBulkEos;
    use crate::eos::VanDerWaals;
    use crate::eos::VinetEos;
    #[test]
    fn test_tait_eos_pressure_at_reference() {
        let eos = TaitEos::water();
        let p = eos.pressure(1000.0);
        assert!(
            p.abs() < 1.0e-6,
            "Pressure at ref density should be ~0, got {p}"
        );
    }
    #[test]
    fn test_tait_eos_sound_speed_at_reference() {
        let eos = TaitEos::water();
        let cs = eos.sound_speed(1000.0);
        assert!(
            (cs - 1500.0).abs() < 1.0e-6,
            "Sound speed at ref should be 1500, got {cs}"
        );
    }
    #[test]
    fn test_tait_eos_pressure_density_roundtrip() {
        let eos = TaitEos::water();
        let rho = 1050.0;
        let p = eos.pressure(rho);
        let rho2 = eos.density_from_pressure(p);
        assert!(
            (rho - rho2).abs() < 1e-3,
            "Roundtrip failed: {rho} → {rho2}"
        );
    }
    #[test]
    fn test_ideal_gas_pressure_positive() {
        let eos = IdealGasEos::air();
        let p = eos.pressure(1.225);
        assert!(
            p > 0.0,
            "Air pressure at standard density should be positive"
        );
        assert!((p - 101325.0).abs() / 101325.0 < 0.01, "Air pressure: {p}");
    }
    #[test]
    fn test_ideal_gas_energy() {
        let eos = IdealGasEos::air();
        let e = eos.internal_energy(1.225, 300.0);
        assert!(e > 0.0);
    }
    #[test]
    fn test_mie_gruneisen_zero_compression() {
        let eos = MieGruneisenEos::copper();
        let p = eos.pressure(eos.rho0);
        assert!(p.abs() < 1e-3, "At rho0, pressure should be ~0, got {p}");
    }
    #[test]
    fn test_mie_gruneisen_compression() {
        let eos = MieGruneisenEos::copper();
        let p_comp = eos.pressure(eos.rho0 * 1.1);
        let p_ref = eos.pressure(eos.rho0);
        assert!(p_comp > p_ref, "Compressed pressure should be higher");
    }
    #[test]
    fn test_stiffened_gas_sound_speed_positive() {
        let eos = StiffenedGasEos::water();
        let cs = eos.sound_speed(1000.0);
        assert!(cs > 0.0, "Sound speed should be positive");
    }
    #[test]
    fn test_van_der_waals_higher_pressure_at_higher_density() {
        let eos = VanDerWaalsEos::co2();
        let p1 = eos.pressure(10.0);
        let p2 = eos.pressure(100.0);
        assert!(p2 > p1, "Higher density should give higher pressure");
    }
    #[test]
    fn test_tillotson_pressure_compressed() {
        let eos = TillotsonEos::granite();
        let p = eos.pressure_from_energy(eos.rho0 * 1.1, 1e6);
        assert!(p > 0.0, "Compressed granite should have positive pressure");
    }
    #[test]
    fn test_polynomial_eos_tension() {
        let eos = PolynomialEos::new(
            1000.0,
            [2.0e9, 1.0e9, 0.5e9],
            [0.5, 0.0, 0.0],
            [1.0e9, 0.5e9],
        );
        let p_ref = eos.pressure(1000.0);
        assert!(
            p_ref.abs() < 1.0,
            "At rho0, polynomial EOS pressure should be 0"
        );
    }
    #[test]
    fn test_bm3_zero_pressure_at_reference() {
        let eos = BirchMurnaghan3Eos::iron();
        let p = eos.pressure(1.0 / eos.v0);
        assert!(p.abs() < 1.0e3, "BM3 pressure at V0 should be ~0, got {p}");
    }
    #[test]
    fn test_bm3_pressure_increases_with_compression() {
        let eos = BirchMurnaghan3Eos::iron();
        let rho0 = 1.0 / eos.v0;
        let p1 = eos.pressure(rho0 * 1.05);
        let p2 = eos.pressure(rho0 * 1.10);
        assert!(p2 > p1, "BM3 pressure should increase with compression");
    }
    #[test]
    fn test_bm3_density_roundtrip() {
        let eos = BirchMurnaghan3Eos::mgo();
        let rho0 = 1.0 / eos.v0;
        let pressure = eos.pressure(rho0 * 1.05);
        let rho_recovered = eos.density_from_pressure(pressure);
        assert!(
            (rho_recovered - rho0 * 1.05).abs() / rho0 < 0.01,
            "roundtrip: {} vs {}",
            rho_recovered,
            rho0 * 1.05
        );
    }
    #[test]
    fn test_bm3_bulk_modulus_at_reference() {
        let eos = BirchMurnaghan3Eos::iron();
        let k = eos.bulk_modulus(eos.v0);
        assert!(
            (k - 170.0e9).abs() / 170.0e9 < 0.05,
            "BM3 bulk modulus at V0 = {k}"
        );
    }
    #[test]
    fn test_bm3_sound_speed_positive() {
        let eos = BirchMurnaghan3Eos::iron();
        let rho0 = 1.0 / eos.v0;
        let cs = eos.sound_speed(rho0);
        assert!(
            cs > 1000.0,
            "Iron sound speed should be > 1000 m/s, got {cs}"
        );
    }
    #[test]
    fn test_vinet_zero_pressure_at_reference() {
        let eos = VinetEos::copper();
        let rho0 = 1.0 / eos.v0;
        let p = eos.pressure(rho0);
        assert!(
            p.abs() < 1.0e3,
            "Vinet pressure at V0 should be ~0, got {p}"
        );
    }
    #[test]
    fn test_vinet_pressure_increases_with_compression() {
        let eos = VinetEos::copper();
        let rho0 = 1.0 / eos.v0;
        let p1 = eos.pressure(rho0 * 1.05);
        let p2 = eos.pressure(rho0 * 1.10);
        assert!(p2 > p1, "Vinet pressure should increase with compression");
    }
    #[test]
    fn test_vinet_density_from_pressure() {
        let eos = VinetEos::diamond();
        let rho0 = 1.0 / eos.v0;
        let pressure = eos.pressure(rho0 * 1.05);
        let rho_back = eos.density_from_pressure(pressure);
        assert!(
            (rho_back - rho0 * 1.05).abs() / rho0 < 0.01,
            "Vinet roundtrip: {} vs {}",
            rho_back,
            rho0 * 1.05
        );
    }
    #[test]
    fn test_hugoniot_zero_pressure_at_zero_velocity() {
        let h = ShockHugoniot::aluminum();
        let p = h.hugoniot_pressure(0.0);
        assert!(
            p.abs() < 1.0,
            "Hugoniot pressure at u_p=0 should be ~0, got {p}"
        );
    }
    #[test]
    fn test_hugoniot_shock_velocity_at_zero_particle_velocity() {
        let h = ShockHugoniot::aluminum();
        let us = h.shock_velocity(0.0);
        assert!((us - h.c0).abs() < 1.0, "Us at Up=0 should be c0, got {us}");
    }
    #[test]
    fn test_hugoniot_pressure_positive_for_positive_up() {
        let h = ShockHugoniot::iron();
        let p = h.hugoniot_pressure(1000.0);
        assert!(
            p > 0.0,
            "Hugoniot pressure should be positive for Up > 0, got {p}"
        );
    }
    #[test]
    fn test_hugoniot_density_increases_under_shock() {
        let h = ShockHugoniot::aluminum();
        let rho = h.post_shock_density(1000.0);
        assert!(
            rho > h.rho0,
            "Post-shock density should exceed reference, got {rho}"
        );
    }
    #[test]
    fn test_hugoniot_energy_positive() {
        let h = ShockHugoniot::aluminum();
        let e = h.hugoniot_energy(1000.0);
        assert!(e > 0.0, "Hugoniot energy should be positive, got {e}");
    }
    #[test]
    fn test_gruneisen_at_reference_volume() {
        let g = GruneisenParameter::new(1.81, 1.0 / 7874.0, 1.0);
        let v0 = g.v0;
        let gamma = g.evaluate(v0);
        assert!((gamma - 1.81).abs() < 1e-10, "Γ at V0 = {gamma}");
    }
    #[test]
    fn test_gruneisen_decreases_with_compression() {
        let g = GruneisenParameter::new(2.0, 1.0 / 7874.0, 1.0);
        let v0 = g.v0;
        let g_ref = g.evaluate(v0);
        let g_comp = g.evaluate(v0 * 0.9);
        assert!(g_comp < g_ref, "Γ should decrease on compression (q>0)");
    }
    #[test]
    fn test_gruneisen_thermal_pressure_positive_above_t0() {
        let g = GruneisenParameter::new(1.81, 1.0 / 7874.0, 1.0);
        let p_th = g.thermal_pressure(7874.0, 450.0, 500.0, 300.0);
        assert!(
            p_th > 0.0,
            "thermal pressure should be positive above T0, got {p_th}"
        );
    }
    #[test]
    fn test_pvt_pressure_positive_at_compression() {
        let pvt = PvtEos::iron();
        let v0 = pvt.cold_eos.v0;
        let p = pvt.pressure(v0 * 0.9, 300.0);
        assert!(
            p > 0.0,
            "P-V-T pressure at compression should be positive, got {p}"
        );
    }
    #[test]
    fn test_pvt_pressure_increases_with_temperature() {
        let pvt = PvtEos::iron();
        let v0 = pvt.cold_eos.v0;
        let p300 = pvt.pressure(v0, 300.0);
        let p600 = pvt.pressure(v0, 600.0);
        assert!(
            p600 > p300,
            "P-V-T pressure should increase with temperature"
        );
    }
    #[test]
    fn test_pvt_density_from_pressure() {
        let pvt = PvtEos::iron();
        let rho0 = 1.0 / pvt.cold_eos.v0;
        let p = pvt.pressure(pvt.cold_eos.v0 * 0.95, 300.0);
        let rho = pvt.density(p, 300.0);
        assert!(
            rho > rho0,
            "compressed density should exceed reference, got {rho}"
        );
    }
    #[test]
    fn test_vdw_at_large_v_approaches_ideal_gas() {
        let vdw = VanDerWaals::with_r(0.3658, 42.9e-6);
        let t = 300.0;
        let r = 8.314_462_618;
        let v_large = 1.0;
        let p_vdw = vdw.pressure(t, v_large, 1.0);
        let p_ideal = r * t / v_large;
        let rel_err = (p_vdw - p_ideal).abs() / p_ideal;
        assert!(
            rel_err < 1e-3,
            "VdW at large V should approach ideal gas, rel_err={rel_err}"
        );
    }
    #[test]
    fn test_vdw_compressibility_at_high_t_approaches_one() {
        let vdw = VanDerWaals::with_r(0.1370, 38.7e-6);
        let t = 10000.0;
        let p = 1.0e5;
        let z = vdw.compressibility(t, p, 1.0);
        assert!(
            (z - 1.0).abs() < 0.01,
            "Z at high T should be close to 1, got {z}"
        );
    }
    #[test]
    fn test_redlich_kwong_pressure_positive() {
        let rk = RedlichKwong::from_critical(304.2, 7.38e6);
        let p = rk.pressure(300.0, 2.0e-4);
        assert!(p > 0.0, "RK pressure should be positive, got {p}");
    }
    #[test]
    fn test_peng_robinson_critical_point() {
        let pr = PengRobinson::new(304.2, 7.38e6, 0.225);
        let b = pr.b_param();
        let vc = 3.0 * b;
        let p = pr.pressure(pr.tc, vc);
        let rel_err = (p - pr.pc).abs() / pr.pc;
        assert!(
            rel_err < 0.5,
            "PR at critical point: p={p}, Pc={}, rel_err={rel_err}",
            pr.pc
        );
    }
    #[test]
    fn test_peng_robinson_alpha_at_tc_equals_one() {
        let pr = PengRobinson::new(304.2, 7.38e6, 0.225);
        let alpha = pr.acentric_factor_correction(pr.tc);
        assert!(
            (alpha - 1.0).abs() < 1e-12,
            "α(Tc) should be 1, got {alpha}"
        );
    }
    #[test]
    fn test_soave_redlich_kwong_pressure_positive() {
        let srk = SoaveRedlichKwong::from_critical(304.2, 7.38e6, 0.225);
        let p = srk.pressure(300.0, 2.0e-4);
        assert!(p > 0.0, "SRK pressure should be positive, got {p}");
    }
    #[test]
    fn test_jwl_pressure_positive() {
        let jwl = JwlEos::tnt();
        let p = jwl.pressure(jwl.rho0, jwl.e0);
        assert!(
            p > 0.0,
            "JWL pressure should be positive at reference density, got {p}"
        );
    }
    #[test]
    fn test_jwl_pressure_increases_with_density() {
        let jwl = JwlEos::tnt();
        let e = jwl.e0;
        let p1 = jwl.pressure(jwl.rho0, e);
        let p2 = jwl.pressure(jwl.rho0 * 1.5, e);
        assert!(p2 > p1, "JWL pressure should increase with density");
    }
    #[test]
    fn test_tait_bulk_eos_at_reference_density_gives_zero() {
        let eos = TaitBulkEos::water();
        let p = eos.pressure(eos.rho0);
        assert!(
            p.abs() < 1.0,
            "Tait bulk EOS at rho0 should give p≈0, got {p}"
        );
    }
    #[test]
    fn test_tait_bulk_eos_pressure_increases_with_density() {
        let eos = TaitBulkEos::water();
        let p1 = eos.pressure(eos.rho0 * 1.01);
        let p2 = eos.pressure(eos.rho0 * 1.05);
        assert!(p2 > p1, "Tait bulk pressure should increase with density");
    }
    #[test]
    fn test_speed_of_sound_ideal_air() {
        let c = speed_of_sound_ideal(1.4, 101325.0, 1.225);
        assert!((c - 340.0).abs() < 5.0, "Speed of sound in air: {c}");
    }
    #[test]
    fn test_speed_of_sound_ideal_zero_density() {
        let c = speed_of_sound_ideal(1.4, 1.0e5, 0.0);
        assert_eq!(c, 0.0, "Sound speed at zero density should be 0");
    }
    #[test]
    fn test_fit_hugoniot_linear_perfect_data() {
        let c0_true = 3000.0_f64;
        let s_true = 1.5_f64;
        let up: Vec<f64> = (0..10).map(|i| i as f64 * 500.0).collect();
        let us: Vec<f64> = up.iter().map(|&u| c0_true + s_true * u).collect();
        let (c0_fit, s_fit) = fit_hugoniot_linear(&up, &us);
        assert!(
            (c0_fit - c0_true).abs() < 1.0,
            "c0 fit: {c0_fit} vs {c0_true}"
        );
        assert!((s_fit - s_true).abs() < 0.001, "s fit: {s_fit} vs {s_true}");
    }
    #[test]
    fn test_fit_hugoniot_single_point() {
        let up = vec![1000.0_f64];
        let us = vec![5000.0_f64];
        let (c0, _s) = fit_hugoniot_linear(&up, &us);
        assert!(c0.is_finite(), "c0 should be finite: {c0}");
    }
    #[test]
    fn test_fit_hugoniot_empty() {
        let (c0, s) = fit_hugoniot_linear(&[], &[]);
        assert_eq!(c0, 0.0);
        assert_eq!(s, 0.0);
    }
    #[test]
    fn test_hugoniot_fit_std_error_decreases_with_more_data() {
        let up5: Vec<f64> = (0..5).map(|i| i as f64 * 100.0).collect();
        let us5: Vec<f64> = up5.iter().map(|&u| 5000.0 + 1.5 * u).collect();
        let up20: Vec<f64> = (0..20).map(|i| i as f64 * 100.0).collect();
        let us20: Vec<f64> = up20.iter().map(|&u| 5000.0 + 1.5 * u).collect();
        let (se_c0_5, se_s_5) = hugoniot_fit_std_error(&up5, &us5);
        let (se_c0_20, se_s_20) = hugoniot_fit_std_error(&up20, &us20);
        let _ = (se_c0_5, se_s_5, se_c0_20, se_s_20);
        assert!(se_c0_5.is_finite() && se_s_5.is_finite());
        assert!(se_c0_20.is_finite() && se_s_20.is_finite());
    }
    #[test]
    fn test_hugoniot_pressure_volume_at_reference_density_zero() {
        let rho0 = 8930.0_f64;
        let v0 = 1.0 / rho0;
        let p = hugoniot_pressure_volume(rho0, 3933.0, 1.5, v0);
        assert!(
            p.abs() < 1.0,
            "Hugoniot pressure at V₀ should be ~0, got {p}"
        );
    }
    #[test]
    fn test_hugoniot_pressure_volume_increases_with_compression() {
        let rho0 = 8930.0_f64;
        let c0 = 3933.0;
        let s = 1.5;
        let v0 = 1.0 / rho0;
        let p1 = hugoniot_pressure_volume(rho0, c0, s, v0 * 0.95);
        let p2 = hugoniot_pressure_volume(rho0, c0, s, v0 * 0.90);
        assert!(
            p2 > p1,
            "Hugoniot pressure should increase with compression"
        );
    }
    #[test]
    fn test_rankine_hugoniot_residuals_zero_for_self_consistent() {
        let rho0 = 1.225_f64;
        let p0 = 101325.0_f64;
        let [r_m, r_mom, _r_e] = rankine_hugoniot_residuals(rho0, p0, rho0, p0, 0.0, 3000.0);
        assert!(r_m.abs() < 1e-6, "mass residual should be ~0: {r_m}");
        assert!(r_mom.abs() < 1.0, "momentum residual should be ~0: {r_mom}");
    }
    #[test]
    fn test_fit_bm3_recovers_known_parameters() {
        let k0_true = 160.0e9_f64;
        let k0p_true = 4.0_f64;
        let v0 = 1.0 / 3580.0_f64;
        let eos = BirchMurnaghan3Eos::new(v0, k0_true, k0p_true);
        let volumes: Vec<f64> = (0..8).map(|i| v0 * (0.85 + i as f64 * 0.02)).collect();
        let pressures: Vec<f64> = volumes
            .iter()
            .map(|&v| eos.pressure_from_volume(v))
            .collect();
        let (k0_fit, k0p_fit) = fit_birch_murnaghan_3rd(v0, &volumes, &pressures, 200);
        assert!(
            (k0_fit - k0_true).abs() / k0_true < 0.05,
            "K₀ fit: {k0_fit:.3e} vs {k0_true:.3e}"
        );
        assert!(
            (k0p_fit - k0p_true).abs() < 0.5,
            "K₀' fit: {k0p_fit} vs {k0p_true}"
        );
    }
    #[test]
    fn test_fit_bm3_empty_returns_defaults() {
        let (k0, k0p) = fit_birch_murnaghan_3rd(1.0e-4, &[], &[], 10);
        assert!((k0 - 100.0e9).abs() < 1.0 && (k0p - 4.0).abs() < 0.01);
    }
    #[test]
    fn test_murnaghan_zero_pressure_at_reference() {
        let eos = MurnaghanEos::aluminum();
        let p = eos.pressure_from_volume(eos.v0);
        assert!(p.abs() < 1.0, "Murnaghan P at V₀ should be ~0, got {p}");
    }
    #[test]
    fn test_murnaghan_pressure_increases_with_compression() {
        let eos = MurnaghanEos::iron();
        let p1 = eos.pressure_from_volume(eos.v0 * 0.98);
        let p2 = eos.pressure_from_volume(eos.v0 * 0.95);
        assert!(p2 > p1, "Murnaghan: more compression → higher pressure");
    }
    #[test]
    fn test_murnaghan_volume_roundtrip() {
        let eos = MurnaghanEos::aluminum();
        let v_test = eos.v0 * 0.95;
        let p = eos.pressure_from_volume(v_test);
        let v_back = eos.volume_from_pressure(p);
        assert!(
            (v_back - v_test).abs() / v_test < 1e-6,
            "Murnaghan V roundtrip: {v_back} vs {v_test}"
        );
    }
    #[test]
    fn test_murnaghan_density_roundtrip() {
        let eos = MurnaghanEos::iron();
        let rho0 = 1.0 / eos.v0;
        let rho_test = rho0 * 1.05;
        let p = eos.pressure(rho_test);
        let rho_back = eos.density_from_pressure(p);
        assert!(
            (rho_back - rho_test).abs() / rho_test < 1e-6,
            "density roundtrip: {rho_back} vs {rho_test}"
        );
    }
    #[test]
    fn test_murnaghan_bulk_modulus_at_reference() {
        let eos = MurnaghanEos::iron();
        let k = eos.bulk_modulus(eos.v0);
        assert!((k - eos.k0).abs() / eos.k0 < 1e-10, "K(V₀) = K₀: {k}");
    }
    #[test]
    fn test_murnaghan_bulk_modulus_increases_with_compression() {
        let eos = MurnaghanEos::iron();
        let k_ref = eos.bulk_modulus(eos.v0);
        let k_comp = eos.bulk_modulus(eos.v0 * 0.9);
        assert!(
            k_comp > k_ref,
            "Bulk modulus should increase on compression"
        );
    }
    #[test]
    fn test_murnaghan_sound_speed_positive() {
        let eos = MurnaghanEos::aluminum();
        let rho0 = 1.0 / eos.v0;
        let c = eos.sound_speed(rho0);
        assert!(c > 1000.0, "Al sound speed > 1000 m/s, got {c}");
    }
    #[test]
    fn test_murnaghan_pressure_trait() {
        let eos = MurnaghanEos::aluminum();
        let rho = 1.0 / (eos.v0 * 0.95);
        let p_trait = eos.pressure(rho);
        let p_direct = eos.pressure_from_volume(eos.v0 * 0.95);
        assert!((p_trait - p_direct).abs() / p_direct.abs().max(1.0) < 1e-10);
    }
    #[test]
    fn test_noble_abel_pressure_positive_at_positive_density() {
        let eos = NobleAbelEos::propellant();
        let p = eos.pressure_at_temperature(10.0, 1000.0);
        assert!(p > 0.0, "Noble-Abel pressure should be positive, got {p}");
    }
    #[test]
    fn test_noble_abel_approaches_ideal_gas_at_low_density() {
        let eos = NobleAbelEos::new(287.0, 1e-4);
        let rho_low = 0.001;
        let p_na = eos.pressure_at_temperature(rho_low, 300.0);
        let p_ideal = rho_low * 287.0 * 300.0;
        let rel_err = (p_na - p_ideal).abs() / p_ideal;
        assert!(
            rel_err < 0.01,
            "Noble-Abel → ideal gas at low density: err={rel_err}"
        );
    }
    #[test]
    fn test_noble_abel_pressure_increases_with_density() {
        let eos = NobleAbelEos::propellant();
        let p1 = eos.pressure_at_temperature(5.0, 500.0);
        let p2 = eos.pressure_at_temperature(10.0, 500.0);
        assert!(p2 > p1, "Higher density → higher pressure");
    }
    #[test]
    fn test_noble_abel_density_roundtrip() {
        let eos = NobleAbelEos::propellant();
        let rho_test = 5.0;
        let p = eos.pressure_at_temperature(rho_test, 500.0);
        let rho_back = eos.density_at_temperature(p, 500.0);
        assert!(
            (rho_back - rho_test).abs() / rho_test < 1e-10,
            "NA density roundtrip: {rho_back} vs {rho_test}"
        );
    }
    #[test]
    fn test_noble_abel_sound_speed_positive() {
        let eos = NobleAbelEos::propellant();
        let cs = eos.sound_speed(5.0);
        assert!(
            cs > 0.0 && cs.is_finite(),
            "Sound speed should be positive and finite"
        );
    }
    #[test]
    fn test_noble_abel_higher_density_gives_higher_sound_speed() {
        let eos = NobleAbelEos::propellant();
        let c1 = eos.sound_speed(1.0);
        let c2 = eos.sound_speed(5.0);
        assert!(
            c2 > c1,
            "Sound speed should increase with density for Noble-Abel"
        );
    }
    #[test]
    fn test_murnaghan_vs_bm3_agree_at_small_compression() {
        let bm3 = BirchMurnaghan3Eos::new(1.0 / 2700.0, 76.0e9, 4.5);
        let murn = MurnaghanEos::new(1.0 / 2700.0, 76.0e9, 4.5);
        let v = bm3.v0 * 0.98;
        let p_bm3 = bm3.pressure_from_volume(v);
        let p_murn = murn.pressure_from_volume(v);
        let rel_diff = (p_bm3 - p_murn).abs() / p_bm3.max(1.0);
        assert!(
            rel_diff < 0.05,
            "BM3 and Murnaghan agree at small compression: rel_diff={rel_diff}"
        );
    }
    #[test]
    fn test_vdw_eos_pressure_at_zero_density_finite() {
        let eos = VanDerWaalsEos::nitrogen();
        let p = eos.pressure(1e-10);
        assert!(
            p.is_finite() || p > 0.0,
            "VdW at near-zero density should give finite/small pressure"
        );
    }
    #[test]
    fn test_redlich_kwong_from_critical_parameters() {
        let rk = RedlichKwong::from_critical(304.2, 7.38e6);
        assert!(rk.a > 0.0, "RK parameter a should be positive");
        assert!(rk.b > 0.0, "RK parameter b should be positive");
    }
    #[test]
    fn test_peng_robinson_a_param_positive() {
        let pr = PengRobinson::new(304.2, 7.38e6, 0.225);
        let a = pr.a_param(300.0);
        assert!(a > 0.0, "PR a(T) should be positive");
    }
    #[test]
    fn test_peng_robinson_a_param_decreases_above_tc() {
        let pr = PengRobinson::new(304.2, 7.38e6, 0.225);
        let a_below = pr.a_param(200.0);
        let a_above = pr.a_param(500.0);
        assert!(a_above < a_below, "PR a(T) should decrease with T above Tc");
    }
    #[test]
    fn test_jwl_pressure_zero_density_guard() {
        let jwl = JwlEos::tnt();
        let p = jwl.pressure(jwl.rho0 * 2.0, jwl.e0);
        assert!(p.is_finite(), "JWL pressure at 2*rho0 should be finite");
    }
    #[test]
    fn test_tillotson_granite_pressure_compressed_positive() {
        let eos = TillotsonEos::granite();
        let p = eos.pressure_from_energy(eos.rho0 * 1.5, 1e8);
        assert!(p > 0.0, "Compressed granite should have positive pressure");
    }
    #[test]
    fn test_tillotson_basalt_pressure_expanded() {
        let eos = TillotsonEos::basalt();
        let p_exp = eos.pressure_from_energy(eos.rho0 * 0.5, eos.e_cv * 2.0);
        assert!(
            p_exp.is_finite(),
            "Expanded hot basalt pressure should be finite"
        );
    }
    #[test]
    fn test_mie_gruneisen_iron_hugoniot_pressure_positive() {
        let eos = MieGruneisenEos::iron();
        let ph = eos.hugoniot_pressure(eos.rho0 * 1.2);
        assert!(
            ph > 0.0,
            "Iron Hugoniot pressure under compression should be positive"
        );
    }
    #[test]
    fn test_mie_gruneisen_pressure_from_energy_at_zero_energy() {
        let eos = MieGruneisenEos::copper();
        let p_e = eos.pressure_from_energy(eos.rho0 * 1.1, 0.0);
        let p = eos.pressure(eos.rho0 * 1.1);
        assert!(
            (p_e - p).abs() < 1e-6,
            "pressure_from_energy at e=0 should match pressure()"
        );
    }
    #[test]
    fn test_stiffened_gas_pressure_increases_with_density() {
        let eos = StiffenedGasEos::water();
        let p1 = eos.pressure(900.0);
        let p2 = eos.pressure(1100.0);
        assert!(p2 > p1, "Stiffened gas: higher density → higher pressure");
    }
    #[test]
    fn test_tait_bulk_eos_sound_speed_increases_with_density() {
        let eos = TaitBulkEos::water();
        let c1 = eos.sound_speed(eos.rho0);
        let c2 = eos.sound_speed(eos.rho0 * 1.05);
        assert!(
            c2 > c1,
            "Sound speed should increase with density (Tait bulk)"
        );
    }
    #[test]
    fn test_birch_murnaghan_mgo_bulk_modulus() {
        let eos = BirchMurnaghan3Eos::mgo();
        let k = eos.bulk_modulus(eos.v0);
        assert!(
            (k - 160.2e9).abs() / 160.2e9 < 0.05,
            "MgO bulk modulus: {k}"
        );
    }
    #[test]
    fn test_gruneisen_parameter_power_law_exponent() {
        let g = GruneisenParameter::new(2.0, 1.0 / 8960.0, 2.0);
        let v0 = g.v0;
        let gamma = g.evaluate(v0 * 2.0);
        assert!(
            (gamma / g.gamma_0 - 4.0).abs() < 1e-10,
            "Γ ∝ V^q: {}",
            gamma / g.gamma_0
        );
    }
    #[test]
    fn test_polynomial_eos_pressure_under_compression() {
        let eos = PolynomialEos::new(
            8000.0,
            [1.5e11, 3.0e10, 1.0e10],
            [0.0, 0.0, 0.0],
            [0.0, 0.0],
        );
        let p = eos.pressure(8400.0);
        assert!(
            p > 0.0,
            "Polynomial EOS: compression gives positive pressure"
        );
    }
    #[test]
    fn test_shock_hugoniot_particle_velocity_from_pressure() {
        let h = ShockHugoniot::aluminum();
        let u_p = 500.0;
        let p = h.hugoniot_pressure(u_p);
        let u_s = h.shock_velocity(u_p);
        let p_check = h.rho0 * u_s * u_p;
        assert!(
            (p - p_check).abs() / p.max(1.0) < 1e-10,
            "Hugoniot consistency"
        );
    }
    #[test]
    fn test_vinet_diamond_pressure_at_reference() {
        let eos = VinetEos::diamond();
        let rho0 = 1.0 / eos.v0;
        let p = eos.pressure(rho0);
        assert!(
            p.abs() < 1.0,
            "Vinet diamond at reference density: p≈0, got {p}"
        );
    }
    #[test]
    fn test_pvt_eos_density_increases_with_pressure() {
        let pvt = PvtEos::iron();
        let rho1 = pvt.density(1.0e9, 300.0);
        let rho2 = pvt.density(5.0e9, 300.0);
        assert!(rho2 > rho1, "Higher pressure → higher density in P-V-T EOS");
    }
}
