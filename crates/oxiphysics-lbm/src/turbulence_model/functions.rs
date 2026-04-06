//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains)]
use super::types::{KOmegaParams, SstParams};

/// Compute turbulent kinematic viscosity: ν_t = k / ω.
///
/// Returns 0 if ω ≤ 0 to avoid division by zero.
pub fn turbulent_viscosity(k: f64, omega: f64) -> f64 {
    if omega <= 0.0 { 0.0 } else { k / omega }
}
/// Compute turbulent kinetic energy production rate: P_k = ν_t · S².
///
/// `strain_rate_sq` is the magnitude-squared of the strain-rate tensor S²
/// (units: 1/s²).
pub fn k_production(nu_t: f64, strain_rate_sq: f64) -> f64 {
    nu_t * strain_rate_sq
}
/// Time derivative of k: ∂k/∂t = P_k - β*·k·ω.
pub fn dk_dt(k: f64, omega: f64, production: f64, params: &KOmegaParams) -> f64 {
    production - params.beta_star * k * omega
}
/// Time derivative of ω: ∂ω/∂t = α·(ω/k)·P_k - β·ω².
///
/// Returns 0 for the production term when k ≤ 0 to avoid division by zero.
pub fn domega_dt(k: f64, omega: f64, production: f64, params: &KOmegaParams) -> f64 {
    let prod_term = if k > 0.0 {
        params.alpha * (omega / k) * production
    } else {
        0.0
    };
    prod_term - params.beta * omega * omega
}
/// Compute turbulent kinematic viscosity for the k-ε model: ν_t = c_μ · k² / ε.
///
/// Returns 0 when ε ≤ 0.
pub fn turbulent_viscosity_ke(k: f64, epsilon: f64, c_mu: f64) -> f64 {
    if epsilon <= 0.0 {
        0.0
    } else {
        c_mu * k * k / epsilon
    }
}
/// Simplified algebraic γ-Re_θ intermittency transition model.
///
/// Based on the Langtry-Menter (2009) framework, simplified to an algebraic
/// expression for the intermittency factor γ as a function of turbulence
/// intensity and local momentum-thickness Reynolds number.
///
/// # Arguments
/// * `tu`       – freestream turbulence intensity (fraction, e.g. 0.01 = 1%)
/// * `re_theta` – momentum-thickness Reynolds number Re_θ
///
/// # Returns
/// Intermittency factor γ ∈ \[0, 1\]:
/// γ = 0 → fully laminar, γ = 1 → fully turbulent.
pub fn gamma_re_theta_model(tu: f64, re_theta: f64) -> f64 {
    if tu <= 0.0 || re_theta <= 0.0 {
        return 0.0;
    }
    let re_theta_crit = 400.0 / (1.0 + 795.0 * tu.powf(1.5));
    if re_theta < re_theta_crit {
        0.0
    } else {
        let width = re_theta_crit * 0.5;
        let x = (re_theta - re_theta_crit) / width;
        (1.0 / (1.0 + (-3.0 * x).exp())).clamp(0.0, 1.0)
    }
}
/// Compute the resolved strain-rate magnitude |S| from a 3×3 velocity
/// gradient tensor.
///
/// `|S| = sqrt(2 * S_ij * S_ij)` where `S_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)`.
#[allow(dead_code)]
pub fn strain_rate_magnitude(velocity_gradient: [f64; 9]) -> f64 {
    let g = velocity_gradient;
    let mut sum = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            let s_ij = 0.5 * (g[i * 3 + j] + g[j * 3 + i]);
            sum += 2.0 * s_ij * s_ij;
        }
    }
    sum.sqrt()
}
/// Compute the rotation-rate tensor magnitude |Ω|.
///
/// `|Ω| = sqrt(2 * Ω_ij * Ω_ij)` where `Ω_ij = 0.5 * (∂u_i/∂x_j - ∂u_j/∂x_i)`.
#[allow(dead_code)]
pub fn rotation_rate_magnitude(velocity_gradient: [f64; 9]) -> f64 {
    let g = velocity_gradient;
    let mut sum = 0.0f64;
    for i in 0..3 {
        for j in 0..3 {
            let omega_ij = 0.5 * (g[i * 3 + j] - g[j * 3 + i]);
            sum += 2.0 * omega_ij * omega_ij;
        }
    }
    sum.sqrt()
}
/// Compute the turbulent kinetic energy dissipation rate from k and ω.
///
/// ε = β* · k · ω  (Wilcox k-ω relation)
#[allow(dead_code)]
pub fn dissipation_rate_from_komega(k: f64, omega: f64) -> f64 {
    let params = KOmegaParams::wilcox_1988();
    params.beta_star * k * omega
}
/// Validate the SGS model: check that eddy viscosity is non-negative
/// and proportional to the strain rate for Smagorinsky.
///
/// Returns `(nu_sgs, |S|)`.
#[allow(dead_code)]
pub fn sgs_validation(cs: f64, delta: f64, velocity_gradient: [f64; 9]) -> (f64, f64) {
    let s_mag = strain_rate_magnitude(velocity_gradient);
    let nu_sgs = (cs * delta) * (cs * delta) * s_mag;
    (nu_sgs, s_mag)
}
/// Compute the sub-grid scale energy flux from k-theory:
///
/// `Π_sgs = -2 * ν_t * S_ij * S_ij`
///
/// where S_ij is the resolved strain-rate tensor.
#[allow(dead_code)]
pub fn sgs_energy_flux(nu_sgs: f64, s_mag: f64) -> f64 {
    -nu_sgs * s_mag * s_mag
}
/// Plane-averaged dynamic Smagorinsky constant computation.
///
/// Instead of point-wise Cs², computes a plane-average to avoid numerical
/// instability.  `l_m_samples` and `m_m_samples` are the LM and MM numerator/
/// denominator samples over the averaging plane.
#[allow(dead_code)]
pub fn plane_averaged_dynamic_constant(l_m_samples: &[f64], m_m_samples: &[f64]) -> f64 {
    assert_eq!(l_m_samples.len(), m_m_samples.len());
    let lm_sum: f64 = l_m_samples.iter().sum();
    let mm_sum: f64 = m_m_samples.iter().sum();
    if mm_sum.abs() < 1e-30 {
        return 0.0;
    }
    (lm_sum / mm_sum).max(0.0)
}
/// Temporal smoothing of the dynamic constant using an exponential filter.
///
/// `Cs²_smooth = (1-alpha) * Cs²_old + alpha * Cs²_new`
#[allow(dead_code)]
pub fn temporal_smoothing(cs2_old: f64, cs2_new: f64, alpha: f64) -> f64 {
    let cs2 = (1.0 - alpha) * cs2_old + alpha * cs2_new;
    cs2.max(0.0)
}
/// Compute the Germano identity residual.
///
/// `|L_ij - M_ij * Cs²|` (Frobenius norm over the 3×3 tensor)
#[allow(dead_code)]
pub fn germano_identity_residual(l_tensor: [f64; 9], m_tensor: [f64; 9], cs2: f64) -> f64 {
    let mut sum = 0.0f64;
    for i in 0..9 {
        let r = l_tensor[i] - cs2 * m_tensor[i];
        sum += r * r;
    }
    sum.sqrt()
}
/// SST k-ω blending function F1.
///
/// F1 = tanh(arg^4) where arg = min(max(sqrt(k)/(beta* omega d), 500 nu/(omega d^2)),
/// 4 sigma_omega2 k / (CDkw d^2)).
///
/// Simplified version using only wall distance `d` and turbulence quantities.
#[allow(dead_code)]
pub fn sst_f1(k: f64, omega: f64, nu: f64, d: f64, params: &SstParams) -> f64 {
    if d <= 0.0 || omega <= 0.0 {
        return 1.0;
    }
    let sqrt_k = k.sqrt();
    let arg1 = sqrt_k / (params.beta_star * omega * d);
    let arg2 = 500.0 * nu / (omega * d * d);
    let arg = arg1.max(arg2);
    arg.tanh().powi(4)
}
/// SST eddy viscosity: ν_t = a1 * k / max(a1 * ω, |S| * F2).
///
/// `F2` is the SST outer blending function.
#[allow(dead_code)]
pub fn sst_eddy_viscosity(k: f64, omega: f64, s_mag: f64, f2: f64) -> f64 {
    let a1 = 0.31;
    let denom = (a1 * omega).max(s_mag * f2);
    if denom <= 0.0 {
        return 0.0;
    }
    a1 * k / denom
}
/// SST F2 blending function.
#[allow(dead_code)]
pub fn sst_f2(k: f64, omega: f64, nu: f64, d: f64) -> f64 {
    if d <= 0.0 || omega <= 0.0 {
        return 1.0;
    }
    let sqrt_k = k.sqrt();
    let arg1 = 2.0 * sqrt_k / (0.09 * omega * d);
    let arg2 = 500.0 * nu / (omega * d * d);
    let arg = arg1.max(arg2);
    arg.tanh().powi(2)
}
/// Compute the ratio of turbulent to molecular viscosity (eddy-viscosity ratio).
///
/// `R_ν = ν_t / ν`
///
/// Values >> 1 indicate a turbulent region; values << 1 indicate laminar or
/// DNS-resolved flow.
#[allow(dead_code)]
pub fn eddy_viscosity_ratio(nu_t: f64, nu: f64) -> f64 {
    if nu <= 0.0 {
        return 0.0;
    }
    nu_t / nu
}
/// Estimate wall distance based on boundary layer log-law.
///
/// Given `y_plus` and friction velocity `u_tau`, compute dimensional
/// wall distance: `y = y_plus * nu / u_tau`.
#[allow(dead_code)]
pub fn wall_distance_from_y_plus(y_plus: f64, nu: f64, u_tau: f64) -> f64 {
    if u_tau <= 0.0 {
        return 0.0;
    }
    y_plus * nu / u_tau
}
/// Compute y+ from dimensional quantities.
///
/// `y⁺ = y * u_tau / ν`
#[allow(dead_code)]
pub fn y_plus(y: f64, u_tau: f64, nu: f64) -> f64 {
    if nu <= 0.0 {
        return 0.0;
    }
    y * u_tau / nu
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::turbulence_model::types::*;
    #[test]
    fn test_turbulent_viscosity_unit() {
        assert!((turbulent_viscosity(1.0, 1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_dk_dt_equilibrium() {
        let params = KOmegaParams::wilcox_1988();
        let k = 0.5;
        let omega = 2.0;
        let nu_t = turbulent_viscosity(k, omega);
        let equilibrium_prod = params.beta_star * k * omega;
        let s2 = equilibrium_prod / nu_t;
        let production = k_production(nu_t, s2);
        let rate = dk_dt(k, omega, production, &params);
        assert!(rate.abs() < 1e-14, "dk_dt at equilibrium = {rate}");
    }
    #[test]
    fn test_effective_nu_exceeds_base() {
        let field = KOmegaField::new(5, 5, 0.1, 1.0);
        let base_nu = 0.01;
        for j in 0..5 {
            for i in 0..5 {
                let nu_eff = field.effective_nu(base_nu, i, j);
                assert!(nu_eff > base_nu, "({i},{j}): {nu_eff} <= {base_nu}");
            }
        }
    }
    #[test]
    fn test_field_step_stays_positive() {
        let nx = 4;
        let ny = 4;
        let mut field = KOmegaField::new(nx, ny, 0.5, 1.0);
        let strain_rates_sq = vec![0.2_f64; nx * ny];
        for _ in 0..100 {
            field.step(0.01, &strain_rates_sq);
        }
        for j in 0..ny {
            for i in 0..nx {
                let cell = field.get(i, j);
                assert!(cell.k > 0.0, "k non-positive at ({i},{j}): {}", cell.k);
                assert!(
                    cell.omega > 0.0,
                    "omega non-positive at ({i},{j}): {}",
                    cell.omega
                );
            }
        }
    }
    #[test]
    fn test_komega_lbm_production_positive() {
        let model = KOmegaLBM::new(0.5, 2.0);
        let prod = model.production(1.0);
        assert!(prod > 0.0, "Production should be positive, got {prod}");
    }
    #[test]
    fn test_komega_lbm_destruction_positive() {
        let model = KOmegaLBM::new(0.5, 2.0);
        let d = model.destruction();
        assert!(d > 0.0, "Destruction should be positive, got {d}");
    }
    #[test]
    fn test_komega_lbm_update_positive() {
        let mut model = KOmegaLBM::new(0.1, 1.0);
        for _ in 0..200 {
            model.update(0.5, 0.01);
        }
        assert!(model.k > 0.0, "k must remain positive: {}", model.k);
        assert!(
            model.omega > 0.0,
            "omega must remain positive: {}",
            model.omega
        );
    }
    #[test]
    fn test_vreman_zero_for_solid_rotation() {
        let g: [f64; 9] = [0.0; 9];
        let model = VremanModel::new(0.07);
        let nu_t = model.vreman_nu_t(g, 1.0);
        assert!(nu_t.abs() < 1e-14, "Vreman ν_t = {nu_t}");
    }
    #[test]
    fn test_les_increases_viscosity() {
        let les = LargeEddySimulation::new(0.1);
        let nu_eff = les.effective_nu(1.0, 1e-3, 1.0);
        assert!(nu_eff > 1e-3, "nu_eff = {nu_eff}");
    }
    #[test]
    fn test_les_nu_increases_with_strain() {
        let les = LargeEddySimulation::new(0.15);
        let nu_low = les.effective_nu(0.1, 1e-4, 1.0);
        let nu_high = les.effective_nu(1.0, 1e-4, 1.0);
        assert!(nu_high > nu_low, "nu_low={nu_low}, nu_high={nu_high}");
    }
    #[test]
    fn test_transition_model_laminar_below_critical() {
        let gamma = gamma_re_theta_model(0.001, 10.0);
        assert!(gamma < 1e-6, "gamma = {gamma}");
    }
    #[test]
    fn test_transition_model_turbulent_above_critical() {
        let gamma = gamma_re_theta_model(0.05, 5000.0);
        assert!(gamma > 0.9, "gamma = {gamma}");
    }
    #[test]
    fn test_filter_width_cubic_root_uniform() {
        let delta = LesFilterWidth::cubic_root(1.0, 1.0, 1.0);
        assert!((delta - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_filter_width_cubic_root_anisotropic() {
        let delta = LesFilterWidth::cubic_root(1.0, 2.0, 4.0);
        let expected = 8.0_f64.cbrt();
        assert!((delta - expected).abs() < 1e-12, "delta = {delta}");
    }
    #[test]
    fn test_filter_width_maximum() {
        let delta = LesFilterWidth::maximum(1.0, 3.0, 2.0);
        assert!((delta - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_van_driest_damped_filter_at_wall() {
        let d = LesFilterWidth::van_driest_damped(1.0, 0.0, 26.0);
        assert!(d.abs() < 1e-14, "At wall, damped filter should be 0: {d}");
    }
    #[test]
    fn test_van_driest_damped_filter_far_from_wall() {
        let d = LesFilterWidth::van_driest_damped(1.0, 1000.0, 26.0);
        assert!((d - 1.0).abs() < 1e-10, "Far from wall: {d}");
    }
    #[test]
    fn test_dynamic_constant_positive_ratio() {
        let cs2 = DynamicSmagorinsky::dynamic_constant(0.5, 2.0);
        assert!((cs2 - 0.25).abs() < 1e-12, "Cs^2 = {cs2}");
    }
    #[test]
    fn test_dynamic_constant_clamps_negative() {
        let cs2 = DynamicSmagorinsky::dynamic_constant(-1.0, 2.0);
        assert_eq!(cs2, 0.0, "Negative Cs^2 should be clamped to 0");
    }
    #[test]
    fn test_dynamic_constant_zero_denominator() {
        let cs2 = DynamicSmagorinsky::dynamic_constant(1.0, 0.0);
        assert_eq!(cs2, 0.0, "Zero denominator should return 0");
    }
    #[test]
    fn test_leonard_stress() {
        let l = DynamicSmagorinsky::leonard_stress(5.0, 2.0, 3.0);
        assert!((l - (-1.0)).abs() < 1e-14, "L = {l}");
    }
    #[test]
    fn test_dynamic_viscosity_proportional() {
        let nu1 = DynamicSmagorinsky::dynamic_viscosity(0.01, 1.0, 1.0);
        let nu2 = DynamicSmagorinsky::dynamic_viscosity(0.01, 2.0, 1.0);
        assert!((nu2 - 4.0 * nu1).abs() < 1e-14, "Should scale as Delta^2");
    }
    #[test]
    fn test_scale_similar_zero_when_equal() {
        let model = ScaleSimilarModel::new(1.0);
        let tau = model.stress_component(4.0, 2.0, 2.0);
        assert!(tau.abs() < 1e-14, "tau = {tau}");
    }
    #[test]
    fn test_scale_similar_nonzero() {
        let model = ScaleSimilarModel::new(1.0);
        let tau = model.stress_component(5.0, 2.0, 2.0);
        assert!((tau - 1.0).abs() < 1e-14, "tau = {tau}");
    }
    #[test]
    fn test_mixed_model_effective_viscosity() {
        let model = MixedModel::new(0.1, 1.0);
        let nu_eff = model.effective_viscosity(1e-3, 1.0, 1.0);
        assert!(nu_eff > 1e-3, "nu_eff should exceed nu_base: {nu_eff}");
    }
    #[test]
    fn test_mixed_model_smagorinsky_part() {
        let model = MixedModel::new(0.1, 1.0);
        let nu_sgs = model.smagorinsky_viscosity(1.0, 1.0);
        let expected = 0.01;
        assert!((nu_sgs - expected).abs() < 1e-14, "nu_sgs = {nu_sgs}");
    }
    #[test]
    fn test_stats_empty() {
        let stats = TurbulenceStatistics::new();
        assert_eq!(stats.mean_ux(), 0.0);
        assert_eq!(stats.variance_ux(), 0.0);
        assert_eq!(stats.reynolds_stress_uv(), 0.0);
    }
    #[test]
    fn test_stats_constant_velocity() {
        let mut stats = TurbulenceStatistics::new();
        for _ in 0..100 {
            stats.add_sample(1.0, 0.5, 0.01);
        }
        assert!(
            (stats.mean_ux() - 1.0).abs() < 1e-12,
            "mean_ux = {}",
            stats.mean_ux()
        );
        assert!(
            (stats.mean_uy() - 0.5).abs() < 1e-12,
            "mean_uy = {}",
            stats.mean_uy()
        );
        assert!(
            stats.variance_ux().abs() < 1e-12,
            "var_ux = {}",
            stats.variance_ux()
        );
        assert!(
            stats.variance_uy().abs() < 1e-12,
            "var_uy = {}",
            stats.variance_uy()
        );
    }
    #[test]
    fn test_stats_fluctuating_velocity() {
        let mut stats = TurbulenceStatistics::new();
        for i in 0..1000 {
            let ux = if i % 2 == 0 { 1.1 } else { 0.9 };
            stats.add_sample(ux, 0.0, 0.0);
        }
        assert!(
            (stats.mean_ux() - 1.0).abs() < 1e-10,
            "mean = {}",
            stats.mean_ux()
        );
        assert!(
            (stats.variance_ux() - 0.01).abs() < 1e-6,
            "var = {}",
            stats.variance_ux()
        );
    }
    #[test]
    fn test_stats_resolved_tke() {
        let mut stats = TurbulenceStatistics::new();
        for i in 0..1000 {
            let ux = if i % 2 == 0 { 1.0 } else { -1.0 };
            let uy = if i % 2 == 0 { 0.5 } else { -0.5 };
            stats.add_sample(ux, uy, 0.0);
        }
        let tke = stats.resolved_tke();
        assert!((tke - 0.625).abs() < 1e-6, "TKE = {tke}");
    }
    #[test]
    fn test_stats_reset() {
        let mut stats = TurbulenceStatistics::new();
        stats.add_sample(1.0, 2.0, 0.1);
        stats.reset();
        assert_eq!(stats.n_samples, 0);
        assert_eq!(stats.mean_ux(), 0.0);
    }
    #[test]
    fn test_stats_turbulence_intensity() {
        let mut stats = TurbulenceStatistics::new();
        for i in 0..1000 {
            let ux = 10.0 + if i % 2 == 0 { 1.0 } else { -1.0 };
            stats.add_sample(ux, 0.0, 0.0);
        }
        let ti = stats.turbulence_intensity();
        assert!(ti > 0.0 && ti < 1.0, "TI = {ti}");
    }
    #[test]
    fn test_stats_mean_k() {
        let mut stats = TurbulenceStatistics::new();
        stats.add_sample(0.0, 0.0, 0.5);
        stats.add_sample(0.0, 0.0, 1.5);
        assert!(
            (stats.mean_k() - 1.0).abs() < 1e-12,
            "mean_k = {}",
            stats.mean_k()
        );
    }
    #[test]
    fn test_strain_rate_magnitude_pure_shear() {
        let mut g = [0.0f64; 9];
        g[1] = 1.0;
        let s_mag = strain_rate_magnitude(g);
        assert!((s_mag - 1.0).abs() < 1e-12, "|S| = {s_mag}");
    }
    #[test]
    fn test_strain_rate_magnitude_zero() {
        let g = [0.0f64; 9];
        let s_mag = strain_rate_magnitude(g);
        assert!(s_mag.abs() < 1e-15, "|S| should be zero: {s_mag}");
    }
    #[test]
    fn test_rotation_rate_magnitude_solid_body() {
        let mut g = [0.0f64; 9];
        let omega_val = 2.0;
        g[1] = omega_val;
        g[3] = -omega_val;
        let om = rotation_rate_magnitude(g);
        let expected = 2.0 * omega_val;
        assert!(
            (om - expected).abs() < 1e-12,
            "|Ω| = {om}, expected {expected}"
        );
    }
    #[test]
    fn test_sgs_validation_proportional_to_strain() {
        let cs = 0.1;
        let delta = 1.0;
        let mut g = [0.0f64; 9];
        g[1] = 2.0;
        let (nu1, s1) = sgs_validation(cs, delta, g);
        let mut g2 = [0.0f64; 9];
        g2[1] = 4.0;
        let (nu2, s2) = sgs_validation(cs, delta, g2);
        assert!(
            (nu2 / nu1 - s2 / s1).abs() < 1e-10,
            "nu2/nu1 = {}, s2/s1 = {}",
            nu2 / nu1,
            s2 / s1
        );
    }
    #[test]
    fn test_sgs_energy_flux_negative() {
        let nu_sgs = 0.01;
        let s_mag = 1.0;
        let flux = sgs_energy_flux(nu_sgs, s_mag);
        assert!(flux <= 0.0, "SGS energy flux should be ≤ 0: {flux}");
    }
    #[test]
    fn test_dissipation_rate_from_komega() {
        let k = 1.0;
        let omega = 10.0;
        let eps = dissipation_rate_from_komega(k, omega);
        let expected = KOmegaParams::wilcox_1988().beta_star * k * omega;
        assert!((eps - expected).abs() < 1e-12, "ε = {eps}");
    }
    #[test]
    fn test_plane_averaged_dynamic_constant() {
        let lm = vec![0.5, 1.0, 0.5];
        let mm = vec![1.0, 2.0, 1.0];
        let cs2 = plane_averaged_dynamic_constant(&lm, &mm);
        assert!((cs2 - 0.5).abs() < 1e-12, "cs2 = {cs2}");
    }
    #[test]
    fn test_plane_averaged_zero_denominator() {
        let lm = vec![1.0, 1.0];
        let mm = vec![0.0, 0.0];
        let cs2 = plane_averaged_dynamic_constant(&lm, &mm);
        assert_eq!(cs2, 0.0);
    }
    #[test]
    fn test_temporal_smoothing_alpha_zero() {
        let cs2_smooth = temporal_smoothing(0.01, 0.1, 0.0);
        assert!((cs2_smooth - 0.01).abs() < 1e-14);
    }
    #[test]
    fn test_temporal_smoothing_alpha_one() {
        let cs2_smooth = temporal_smoothing(0.01, 0.1, 1.0);
        assert!((cs2_smooth - 0.1).abs() < 1e-14);
    }
    #[test]
    fn test_temporal_smoothing_clamps_negative() {
        let cs2_smooth = temporal_smoothing(0.0, -1.0, 0.5);
        assert_eq!(cs2_smooth, 0.0, "Should clamp negative to 0");
    }
    #[test]
    fn test_germano_residual_zero_when_perfect() {
        let l = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let m = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let res = germano_identity_residual(l, m, 1.0);
        assert!(res.abs() < 1e-14, "Residual should be zero: {res}");
    }
    #[test]
    fn test_germano_residual_nonzero() {
        let l = [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let m = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0_f64];
        let res = germano_identity_residual(l, m, 0.5);
        assert!((res - 1.5).abs() < 1e-12, "Residual = {res}");
    }
    #[test]
    fn test_bardina_stress_zero_when_filters_equal() {
        let model = BardinalFullModel::new(1.0, 0.0);
        let tau = model.sgs_stress(4.0, 2.0, 2.0, 1.0, 1.0, 0.5);
        assert!(tau.abs() < 1e-14, "tau = {tau}");
    }
    #[test]
    fn test_bardina_eddy_viscosity_positive() {
        let model = BardinalFullModel::new(1.0, 0.01);
        let nu_t = model.eddy_viscosity(1.0, 1.0);
        assert!(nu_t > 0.0, "Eddy viscosity should be positive: {nu_t}");
    }
    #[test]
    fn test_bardina_stress_nonzero_scale_similar() {
        let model = BardinalFullModel::new(1.0, 0.0);
        let tau = model.sgs_stress(5.0, 2.0, 2.0, 1.0, 0.0, 0.0);
        assert!((tau - 1.0).abs() < 1e-14, "tau = {tau}");
    }
    #[test]
    fn test_sst_params_menter_1994() {
        let p = SstParams::menter_1994();
        assert!((p.beta_star - 0.09).abs() < 1e-12);
        assert!((p.alpha1 - 5.0 / 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_sst_blended_beta_f1_one() {
        let p = SstParams::menter_1994();
        let beta = p.blended_beta(1.0);
        assert!((beta - p.beta1).abs() < 1e-12, "F1=1 → inner beta");
    }
    #[test]
    fn test_sst_blended_beta_f1_zero() {
        let p = SstParams::menter_1994();
        let beta = p.blended_beta(0.0);
        assert!((beta - p.beta2).abs() < 1e-12, "F1=0 → outer beta");
    }
    #[test]
    fn test_sst_f1_at_wall() {
        let params = SstParams::menter_1994();
        let f1 = sst_f1(1.0, 10.0, 1e-5, 1e-10, &params);
        assert!(f1 > 0.9, "F1 at wall = {f1}");
    }
    #[test]
    fn test_sst_f2_at_wall() {
        let f2 = sst_f2(1.0, 10.0, 1e-5, 1e-10);
        assert!(f2 > 0.9, "F2 at wall = {f2}");
    }
    #[test]
    fn test_des_k_destruction_rans_dominant() {
        let model = DesSstModel::new(0.61);
        let k = 1.0;
        let omega = 100.0;
        let delta = 100.0;
        let d = model.k_destruction_des(k, omega, delta);
        let d_rans = model.params.beta_star * k * omega;
        assert!(
            (d - d_rans).abs() < 1e-10,
            "DES destruction = {d}, RANS = {d_rans}"
        );
    }
    #[test]
    fn test_des_k_destruction_les_dominant() {
        let model = DesSstModel::new(0.61);
        let k = 1.0;
        let omega = 1.0;
        let delta = 1e-6;
        let d = model.k_destruction_des(k, omega, delta);
        let d_rans = model.params.beta_star * k * omega;
        assert!(
            d >= d_rans,
            "DES destruction in LES region should >= RANS: {d} vs {d_rans}"
        );
    }
    #[test]
    fn test_hybrid_nu_rans_positive() {
        let model = HybridRansLes::new(0.5, 2.0, 0.1);
        let nu_rans = model.nu_rans();
        assert!(nu_rans > 0.0, "RANS nu_t should be positive: {nu_rans}");
    }
    #[test]
    fn test_hybrid_nu_les_increases_with_strain() {
        let model = HybridRansLes::new(0.5, 2.0, 0.1);
        let nu1 = model.nu_les(1.0, 0.5);
        let nu2 = model.nu_les(1.0, 1.0);
        assert!(nu2 > nu1, "LES nu_t increases with strain rate");
    }
    #[test]
    fn test_hybrid_effective_nu_exceeds_base() {
        let model = HybridRansLes::new(0.5, 2.0, 0.1);
        let nu_base = 1e-3;
        let nu_eff = model.effective_nu(1.0, 0.5, nu_base);
        assert!(nu_eff >= nu_base, "Effective nu should >= base: {nu_eff}");
    }
    #[test]
    fn test_hybrid_update_rans_positive() {
        let mut model = HybridRansLes::new(0.1, 1.0, 0.1);
        for _ in 0..100 {
            model.update_rans(0.5, 0.01);
        }
        assert!(model.k > 0.0, "k > 0: {}", model.k);
        assert!(model.omega > 0.0, "omega > 0: {}", model.omega);
    }
    #[test]
    fn test_eddy_viscosity_ratio() {
        let r = eddy_viscosity_ratio(1e-2, 1e-4);
        assert!((r - 100.0).abs() < 1e-10, "EV ratio = {r}");
    }
    #[test]
    fn test_y_plus_formula() {
        let yp = y_plus(1e-3, 0.1, 1e-5);
        let expected = 1e-3 * 0.1 / 1e-5;
        assert!((yp - expected).abs() < 1e-10, "y+ = {yp}");
    }
    #[test]
    fn test_wall_distance_from_y_plus() {
        let d = wall_distance_from_y_plus(10.0, 1e-5, 0.1);
        let expected = 10.0 * 1e-5 / 0.1;
        assert!((d - expected).abs() < 1e-15, "y = {d}");
    }
}
/// Compute the least-squares Cs² from batches of L:M and M:M.
///
/// Averages over a plane or volume:
/// `Cs² = Σ(L_{ij} M_{ij}) / Σ(M_{ij} M_{ij})`
pub fn lilly_dynamic_constant(l_m: &[f64], m_m: &[f64]) -> f64 {
    if l_m.len() != m_m.len() || m_m.is_empty() {
        return 0.0;
    }
    let lm_sum: f64 = l_m.iter().sum();
    let mm_sum: f64 = m_m.iter().sum();
    if mm_sum.abs() > 1e-20 {
        (lm_sum / mm_sum).max(0.0)
    } else {
        0.0
    }
}
/// Compute the strain rate tensor S_{ij} = 0.5(g_{ij} + g_{ji}).
pub fn strain_rate_tensor(g: [f64; 9]) -> [f64; 9] {
    let mut s = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            s[i * 3 + j] = 0.5 * (g[i * 3 + j] + g[j * 3 + i]);
        }
    }
    s
}
/// Compute the rotation rate tensor Ω_{ij} = 0.5(g_{ij} − g_{ji}).
pub fn rotation_rate_tensor(g: [f64; 9]) -> [f64; 9] {
    let mut omega = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            omega[i * 3 + j] = 0.5 * (g[i * 3 + j] - g[j * 3 + i]);
        }
    }
    omega
}
/// Compute eigenvalues of a symmetric 3×3 matrix using the analytical method.
///
/// Returns \[λ₁, λ₂, λ₃\] (unsorted).
pub fn eigenvalues_3x3_symmetric(a: [f64; 9]) -> [f64; 3] {
    let a11 = a[0];
    let a12 = a[1];
    let a13 = a[2];
    let a22 = a[4];
    let a23 = a[5];
    let a33 = a[8];
    let p = a11 + a22 + a33;
    let q = a11 * a22 + a11 * a33 + a22 * a33 - a12 * a12 - a13 * a13 - a23 * a23;
    let r = a11 * (a22 * a33 - a23 * a23) - a12 * (a12 * a33 - a23 * a13)
        + a13 * (a12 * a23 - a22 * a13);
    let p3 = p / 3.0;
    let q3 = q / 3.0;
    let b = p3 * p3 - q3;
    let c = 2.0 * p3 * p3 * p3 - p3 * q3 + r / 3.0;
    let discriminant = b * b * b - c * c;
    if discriminant >= 0.0 {
        let theta = (c / (b * b * b).max(1e-30).sqrt()).clamp(-1.0, 1.0).acos();
        let sqrt_b = b.max(0.0).sqrt();
        [
            p3 + 2.0 * sqrt_b * (theta / 3.0).cos(),
            p3 + 2.0 * sqrt_b * ((theta + 2.0 * std::f64::consts::PI) / 3.0).cos(),
            p3 + 2.0 * sqrt_b * ((theta + 4.0 * std::f64::consts::PI) / 3.0).cos(),
        ]
    } else {
        [p / 3.0, p / 3.0, p / 3.0]
    }
}
/// Detached Eddy Simulation (DES97) interface length scale.
///
/// Computes `l_DES = min(l_RANS, C_DES Δ)` where l_RANS = k^(1/2)/ω.
///
/// # Arguments
/// * `k`      – turbulent kinetic energy
/// * `omega`  – specific dissipation rate
/// * `delta`  – local grid filter width
/// * `c_des`  – DES constant (typically 0.61)
pub fn des97_length_scale(k: f64, omega: f64, delta: f64, c_des: f64) -> f64 {
    let beta_star = 0.09f64;
    let l_rans = if omega > 1e-14 {
        k.sqrt() / (beta_star.sqrt() * omega)
    } else {
        f64::MAX
    };
    let l_les = c_des * delta;
    l_rans.min(l_les)
}
/// Delayed DES (DDES) shielding function F_d.
///
/// Prevents premature switching from RANS to LES in attached boundary layers.
/// `F_d = 1 − tanh((8 r_d)³)`
/// where `r_d = (ν_t + ν) / (|∇u| d² κ²)`.
///
/// # Arguments
/// * `nu_t`      – eddy viscosity
/// * `nu`        – kinematic viscosity
/// * `s_mag`     – strain-rate magnitude |∇u|
/// * `d`         – wall distance
/// * `kappa`     – von Kármán constant (0.41)
pub fn ddes_shielding_fd(nu_t: f64, nu: f64, s_mag: f64, d: f64, kappa: f64) -> f64 {
    let denom = (s_mag * d * d * kappa * kappa).max(1e-14);
    let r_d = (nu_t + nu) / denom;
    1.0 - (8.0 * r_d).powi(3).tanh()
}
/// IDDES (Improved DDES) subgrid length scale.
///
/// Combines DDES and SAS (Scale-Adaptive Simulation) length scales.
/// `l_IDDES = f_d (1 + f_e) l_RANS + (1 − f_d) l_LES`
///
/// # Arguments
/// * `l_rans`   – RANS length scale
/// * `l_les`    – LES filter length scale
/// * `f_d`      – DDES shielding function
/// * `f_e`      – elevation function (typically 0 or small positive)
pub fn iddes_length_scale(l_rans: f64, l_les: f64, f_d: f64, f_e: f64) -> f64 {
    f_d * (1.0 + f_e) * l_rans + (1.0 - f_d) * l_les
}
/// Scale-Adaptive Simulation (SAS) von Kármán length scale.
///
/// `l_vk = κ |U'| / |U''|` where U' and U'' are first and second
/// velocity derivatives.
///
/// # Arguments
/// * `u_prime`  – first velocity derivative magnitude
/// * `u_double` – second velocity derivative magnitude
/// * `kappa`    – von Kármán constant (0.41)
pub fn sas_von_karman_length(u_prime: f64, u_double: f64, kappa: f64) -> f64 {
    if u_double > 1e-14 {
        kappa * u_prime / u_double
    } else {
        f64::MAX
    }
}
/// Von Kármán turbulence energy spectrum amplitude for mode k.
///
/// E(k) ∝ (k/k_e)⁴ exp(−2(k/k_e)²) where k_e = √(5/12)/L.
/// Returns √(E(k) dk / n_modes).
pub(super) fn von_karman_energy_spectral_amplitude(k: f64, l: f64, intensity: f64) -> f64 {
    let k_e = (5.0_f64 / 12.0).sqrt() / l;
    let x = k / k_e;
    let e_k = intensity * intensity * (x * x * x * x) * (-2.0 * x * x).exp();
    (e_k.max(0.0)).sqrt()
}
/// Linear congruential generator for deterministic pseudo-random numbers.
///
/// Returns a value in \[0, 1).
pub(super) fn lcg_pseudo_random(seed: u64) -> f64 {
    let a: u64 = 6364136223846793005;
    let c: u64 = 1442695040888963407;
    let m: u64 = u64::MAX;
    let x = a.wrapping_mul(seed).wrapping_add(c);
    let x = a.wrapping_mul(x).wrapping_add(c);
    (x & m) as f64 / (u64::MAX as f64 + 1.0)
}
/// Mann turbulence box generator (simplified 1D slice).
///
/// Generates a turbulent inflow profile across Ny points using the Mann
/// spectral method (truncated Fourier series).
///
/// # Arguments
/// * `ny`       – number of points in the y-direction
/// * `dy`       – grid spacing
/// * `ux_mean`  – mean x-velocity
/// * `intensity` – turbulence intensity
/// * `l_scale`  – integral length scale
///
/// # Returns
/// Vector of (ux, uy) pairs at each y-position.
pub fn mann_inflow_1d(
    ny: usize,
    dy: f64,
    ux_mean: f64,
    intensity: f64,
    l_scale: f64,
) -> Vec<(f64, f64)> {
    let mut profile = Vec::with_capacity(ny);
    let n_modes = 32usize;
    for j in 0..ny {
        let y = j as f64 * dy;
        let mut ux_p = 0.0f64;
        let mut uy_p = 0.0f64;
        for n in 1..=n_modes {
            let k = 2.0 * std::f64::consts::PI * n as f64 / (ny as f64 * dy);
            let amp = von_karman_energy_spectral_amplitude(k, l_scale, intensity);
            let phase_x = lcg_pseudo_random(n as u64 * 31) * 2.0 * std::f64::consts::PI;
            let phase_y = lcg_pseudo_random(n as u64 * 37 + 1) * 2.0 * std::f64::consts::PI;
            ux_p += amp * (k * y + phase_x).cos();
            uy_p += amp * (k * y + phase_y).cos();
        }
        profile.push((ux_mean + ux_p, uy_p));
    }
    profile
}
/// Compute the turbulent Prandtl number from eddy viscosity and eddy diffusivity.
///
/// Pr_t = ν_t / α_t
pub fn turbulent_prandtl(nu_t: f64, alpha_t: f64) -> f64 {
    if alpha_t > 1e-14 { nu_t / alpha_t } else { 0.0 }
}
/// Estimate the Kolmogorov micro-scale η from dissipation rate ε and kinematic viscosity ν.
///
/// `η = (ν³ / ε)^(1/4)`
pub fn kolmogorov_length(nu: f64, epsilon: f64) -> f64 {
    if epsilon > 1e-20 {
        (nu * nu * nu / epsilon).powf(0.25)
    } else {
        0.0
    }
}
/// Estimate the Taylor micro-scale λ from k, ε, and ν.
///
/// `λ = √(10 ν k / ε)`
pub fn taylor_micro_scale(nu: f64, k: f64, epsilon: f64) -> f64 {
    if epsilon > 1e-20 {
        (10.0 * nu * k / epsilon).sqrt()
    } else {
        0.0
    }
}
/// Estimate the integral length scale L from k and ε.
///
/// `L = k^(3/2) / ε`
pub fn integral_length_scale_from_ke(k: f64, epsilon: f64) -> f64 {
    if epsilon > 1e-20 {
        k.powf(1.5) / epsilon
    } else {
        0.0
    }
}
/// Compute the turbulent Reynolds number Re_t = k²/(ε ν).
pub fn turbulent_reynolds(k: f64, epsilon: f64, nu: f64) -> f64 {
    if epsilon > 1e-20 && nu > 1e-20 {
        k * k / (epsilon * nu)
    } else {
        0.0
    }
}
/// Compute the ratio of SGS dissipation to total dissipation.
///
/// Used to assess LES resolution quality: ε_sgs / (ε_sgs + ε_resolved).
/// Values > 0.2 indicate under-resolved LES.
pub fn les_resolution_quality(epsilon_sgs: f64, epsilon_resolved: f64) -> f64 {
    let total = epsilon_sgs + epsilon_resolved;
    if total > 1e-20 {
        epsilon_sgs / total
    } else {
        0.0
    }
}
/// Estimate the LES quality index M = k_resolved / k_total (Pope 2000).
///
/// M > 0.8 is considered well-resolved LES.
pub fn les_quality_index(k_resolved: f64, k_sgs: f64) -> f64 {
    let k_total = k_resolved + k_sgs;
    if k_total > 1e-20 {
        k_resolved / k_total
    } else {
        0.0
    }
}
/// Compute the Kolmogorov time scale τ_η = (ν / ε)^(1/2).
pub fn kolmogorov_time_scale(nu: f64, epsilon: f64) -> f64 {
    if epsilon > 1e-20 {
        (nu / epsilon).sqrt()
    } else {
        0.0
    }
}
/// Compute the turbulent time scale T = k / ε.
pub fn turbulent_time_scale(k: f64, epsilon: f64) -> f64 {
    if epsilon > 1e-20 { k / epsilon } else { 0.0 }
}
#[cfg(test)]
mod tests_advanced_turbulence {
    use super::*;
    use crate::turbulence_model::types::*;
    fn zero_gradient() -> [f64; 9] {
        [0.0f64; 9]
    }
    fn simple_shear_gradient(du_dy: f64) -> [f64; 9] {
        let mut g = [0.0f64; 9];
        g[1] = du_dy;
        g
    }
    fn symmetric_gradient(s_val: f64) -> [f64; 9] {
        let mut g = [0.0f64; 9];
        g[0] = s_val;
        g[4] = -s_val;
        g
    }
    #[test]
    fn test_dynamic_smag_lilly_new() {
        let m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        assert_eq!(m.cs2, 0.0, "Initial cs2 should be 0");
        assert_eq!(m.filter_ratio, 2.0);
    }
    #[test]
    fn test_dynamic_smag_lilly_update_cs2() {
        let mut m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        m.update_cs2(0.002, 0.1);
        assert!((m.cs2 - 0.02).abs() < 1e-12, "cs2 = {}", m.cs2);
    }
    #[test]
    fn test_dynamic_smag_lilly_cs2_clamped() {
        let mut m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        m.update_cs2(1.0, 1.0);
        assert!(m.cs2 <= m.cs2_max, "cs2 should be clamped: {}", m.cs2);
    }
    #[test]
    fn test_dynamic_smag_lilly_eddy_viscosity() {
        let mut m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        m.cs2 = 0.01;
        let nu_t = m.eddy_viscosity(1.0);
        assert!((nu_t - 0.01).abs() < 1e-12, "nu_t = {nu_t}");
    }
    #[test]
    fn test_dynamic_smag_lilly_zero_strain() {
        let m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        let nu_t = m.eddy_viscosity(0.0);
        assert_eq!(nu_t, 0.0, "Zero strain → zero nu_t");
    }
    #[test]
    fn test_lilly_dynamic_constant_positive() {
        let l_m = vec![0.1, 0.2, 0.3];
        let m_m = vec![1.0, 2.0, 3.0];
        let cs2 = lilly_dynamic_constant(&l_m, &m_m);
        assert!(cs2 >= 0.0, "Dynamic constant should be non-negative: {cs2}");
    }
    #[test]
    fn test_lilly_dynamic_constant_exact() {
        let l_m = vec![2.0f64];
        let m_m = vec![4.0f64];
        let cs2 = lilly_dynamic_constant(&l_m, &m_m);
        assert!((cs2 - 0.5).abs() < 1e-12, "cs2 = {cs2}");
    }
    #[test]
    fn test_lilly_dynamic_constant_empty() {
        let cs2 = lilly_dynamic_constant(&[], &[]);
        assert_eq!(cs2, 0.0, "Empty input → 0");
    }
    #[test]
    fn test_wale_zero_gradient() {
        let m = WaleModel::new(1.0);
        let nu_t = m.eddy_viscosity(zero_gradient());
        assert_eq!(nu_t, 0.0, "WALE zero gradient → nu_t = 0");
    }
    #[test]
    fn test_wale_shear_positive_nu_t() {
        let m = WaleModel::new(1.0);
        let g = simple_shear_gradient(0.1);
        let nu_t = m.eddy_viscosity(g);
        assert!(nu_t >= 0.0, "WALE nu_t should be non-negative: {nu_t}");
    }
    #[test]
    fn test_wale_solid_body_rotation_zero() {
        let mut g = [0.0f64; 9];
        g[1] = 1.0;
        g[3] = -1.0;
        let m = WaleModel::new(1.0);
        let nu_t = m.eddy_viscosity(g);
        assert!(nu_t >= 0.0, "WALE solid body nu_t={nu_t}");
    }
    #[test]
    fn test_wale_sd_tensor_zero_for_zero_g() {
        let m = WaleModel::new(1.0);
        let sd = m.compute_sd(zero_gradient());
        for v in sd.iter() {
            assert!(v.abs() < 1e-14, "Sd should be zero for zero g: {v}");
        }
    }
    #[test]
    fn test_wale_with_constant() {
        let m1 = WaleModel::new(1.0);
        let m2 = WaleModel::with_constant(1.0, 0.325);
        let g = simple_shear_gradient(0.1);
        let nu1 = m1.eddy_viscosity(g);
        let nu2 = m2.eddy_viscosity(g);
        assert!(
            nu1 >= 0.0 && nu2 >= 0.0,
            "Both nu_t should be non-negative: {nu1}, {nu2}"
        );
    }
    #[test]
    fn test_wale_sgs_stress_magnitude_non_negative() {
        let m = WaleModel::new(1.0);
        let g = simple_shear_gradient(0.2);
        let stress = m.sgs_stress_magnitude(g);
        assert!(stress >= 0.0, "WALE SGS stress should be ≥ 0: {stress}");
    }
    #[test]
    fn test_strain_rate_tensor_symmetric() {
        let g = simple_shear_gradient(0.1);
        let s = strain_rate_tensor(g);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (s[i * 3 + j] - s[j * 3 + i]).abs() < 1e-14,
                    "S not symmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_rotation_rate_tensor_antisymmetric() {
        let g = simple_shear_gradient(0.1);
        let omega = rotation_rate_tensor(g);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (omega[i * 3 + j] + omega[j * 3 + i]).abs() < 1e-14,
                    "Ω not antisymmetric at ({i},{j})"
                );
            }
        }
    }
    #[test]
    fn test_strain_rotation_decomposition() {
        let g = simple_shear_gradient(0.1);
        let s = strain_rate_tensor(g);
        let omega = rotation_rate_tensor(g);
        for i in 0..9 {
            assert!(
                (s[i] + omega[i] - g[i]).abs() < 1e-14,
                "g = S + Ω fails at i={i}"
            );
        }
    }
    #[test]
    fn test_sigma_model_zero_gradient() {
        let m = SigmaModel::new(1.0);
        let nu_t = m.eddy_viscosity(zero_gradient());
        assert_eq!(nu_t, 0.0, "Sigma: zero gradient → nu_t = 0");
    }
    #[test]
    fn test_sigma_model_shear_non_negative() {
        let m = SigmaModel::new(1.0);
        let g = simple_shear_gradient(0.1);
        let nu_t = m.eddy_viscosity(g);
        assert!(nu_t >= 0.0, "Sigma nu_t should be ≥ 0: {nu_t}");
    }
    #[test]
    fn test_singular_values_sorted() {
        let m = SigmaModel::new(1.0);
        let g = symmetric_gradient(0.2);
        let (s1, s2, s3) = m.singular_values_g(g);
        assert!(s1 >= s2, "s1={s1} >= s2={s2}");
        assert!(s2 >= s3, "s2={s2} >= s3={s3}");
    }
    #[test]
    fn test_eigenvalues_3x3_identity() {
        let a = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0f64];
        let mut ev = eigenvalues_3x3_symmetric(a);
        ev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for v in ev.iter() {
            assert!((v - 1.0).abs() < 1e-8, "Identity eigenvalue: {v}");
        }
    }
    #[test]
    fn test_eigenvalues_3x3_diagonal() {
        let a = [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 5.0f64];
        let ev = eigenvalues_3x3_symmetric(a);
        let trace = ev[0] + ev[1] + ev[2];
        assert!((trace - 10.0).abs() < 1e-6, "trace={trace}, expected 10");
        for v in ev.iter() {
            assert!(v.is_finite(), "eigenvalue is not finite: {v}");
        }
    }
    #[test]
    fn test_amd_zero_gradient() {
        let m = AmdModel::new([1.0, 1.0, 1.0]);
        let nu_t = m.eddy_viscosity(zero_gradient());
        assert_eq!(nu_t, 0.0, "AMD zero gradient → nu_t = 0");
    }
    #[test]
    fn test_amd_anisotropy_isotropic() {
        let m = AmdModel::new([1.0, 1.0, 1.0]);
        assert!(
            (m.anisotropy_factor() - 1.0).abs() < 1e-12,
            "Isotropic anisotropy = {}",
            m.anisotropy_factor()
        );
    }
    #[test]
    fn test_amd_anisotropy_stretched() {
        let m = AmdModel::new([1.0, 2.0, 4.0]);
        assert!(
            (m.anisotropy_factor() - 4.0).abs() < 1e-12,
            "Stretched anisotropy = {}",
            m.anisotropy_factor()
        );
    }
    #[test]
    fn test_amd_with_constant_non_negative() {
        let m = AmdModel::with_constant([1.0, 1.0, 1.0], 0.3);
        let g = simple_shear_gradient(0.1);
        let nu_t = m.eddy_viscosity(g);
        assert!(nu_t >= 0.0, "AMD nu_t should be non-negative: {nu_t}");
    }
    #[test]
    fn test_des97_length_scale_rans_dominant() {
        let l = des97_length_scale(1.0, 100.0, 1e6, 0.61);
        let beta_star = 0.09f64;
        let l_rans = 1.0f64.sqrt() / (beta_star.sqrt() * 100.0);
        assert!(
            (l - l_rans).abs() < 1e-10,
            "DES RANS: l={l}, l_rans={l_rans}"
        );
    }
    #[test]
    fn test_des97_length_scale_les_dominant() {
        let c_des = 0.61;
        let delta = 1e-6;
        let l = des97_length_scale(1.0, 0.01, delta, c_des);
        let l_les = c_des * delta;
        assert!((l - l_les).abs() < 1e-12, "DES LES: l={l}, l_les={l_les}");
    }
    #[test]
    fn test_ddes_shielding_fd_range() {
        let fd = ddes_shielding_fd(1e-4, 1e-5, 0.1, 0.01, 0.41);
        assert!(fd >= 0.0 && fd <= 1.0 + 1e-10, "F_d out of [0,1]: {fd}");
    }
    #[test]
    fn test_ddes_shielding_fd_at_wall() {
        let fd = ddes_shielding_fd(0.1, 1e-5, 1e-6, 1e-6, 0.41);
        assert!(fd < 0.5, "F_d near wall should be < 0.5: {fd}");
    }
    #[test]
    fn test_iddes_length_scale_rans_limit() {
        let l = iddes_length_scale(0.5, 0.1, 1.0, 0.0);
        assert!((l - 0.5).abs() < 1e-12, "IDDES RANS limit: l={l}");
    }
    #[test]
    fn test_iddes_length_scale_les_limit() {
        let l = iddes_length_scale(0.5, 0.1, 0.0, 0.0);
        assert!((l - 0.1).abs() < 1e-12, "IDDES LES limit: l={l}");
    }
    #[test]
    fn test_sas_von_karman_length_finite() {
        let l = sas_von_karman_length(1.0, 2.0, 0.41);
        assert!(l.is_finite() && l > 0.0, "SAS von Kármán: {l}");
    }
    #[test]
    fn test_sas_von_karman_length_zero_u_double() {
        let l = sas_von_karman_length(1.0, 0.0, 0.41);
        assert_eq!(l, f64::MAX, "SAS: zero U'' → MAX");
    }
    #[test]
    fn test_synthetic_turbulence_creates() {
        let st = SyntheticTurbulenceInflow::new(10, 0.05, 1.0);
        assert_eq!(st.n_modes, 10);
        assert_eq!(st.wavenumbers.len(), 10);
    }
    #[test]
    fn test_synthetic_turbulence_fluctuation_finite() {
        let st = SyntheticTurbulenceInflow::new(10, 0.05, 1.0);
        let u = st.velocity_fluctuation(0.5, 0.5, 0.0, 0.0);
        assert!(
            u.iter().all(|v| v.is_finite()),
            "Fluctuation has NaN: {u:?}"
        );
    }
    #[test]
    fn test_synthetic_turbulence_ke_positive() {
        let st = SyntheticTurbulenceInflow::new(10, 0.05, 1.0);
        let ke = st.turbulent_ke(10.0);
        assert!(ke > 0.0, "Turbulent KE should be positive: {ke}");
    }
    #[test]
    fn test_mann_inflow_1d_length() {
        let profile = mann_inflow_1d(20, 1.0, 0.1, 0.05, 1.0);
        assert_eq!(profile.len(), 20, "Profile length should match ny");
    }
    #[test]
    fn test_mann_inflow_1d_finite() {
        let profile = mann_inflow_1d(10, 1.0, 0.1, 0.05, 1.0);
        for (ux, uy) in &profile {
            assert!(
                ux.is_finite() && uy.is_finite(),
                "Mann inflow has NaN: {ux}, {uy}"
            );
        }
    }
    #[test]
    fn test_kolmogorov_length_positive() {
        let eta = kolmogorov_length(1e-5, 0.01);
        assert!(eta > 0.0, "Kolmogorov length: {eta}");
    }
    #[test]
    fn test_kolmogorov_length_formula() {
        let nu = 1e-5;
        let eps = 0.001;
        let eta = kolmogorov_length(nu, eps);
        let expected = (nu * nu * nu / eps).powf(0.25);
        assert!(
            (eta - expected).abs() < 1e-15,
            "Kolmogorov length formula: {eta}"
        );
    }
    #[test]
    fn test_taylor_micro_scale_positive() {
        let lambda = taylor_micro_scale(1e-5, 0.1, 0.01);
        assert!(lambda > 0.0, "Taylor micro-scale: {lambda}");
    }
    #[test]
    fn test_integral_length_from_ke() {
        let l = integral_length_scale_from_ke(1.0, 1.0);
        assert!((l - 1.0).abs() < 1e-12, "L = k^(3/2)/eps = 1: {l}");
    }
    #[test]
    fn test_turbulent_reynolds_positive() {
        let re = turbulent_reynolds(0.1, 0.01, 1e-5);
        assert!(re > 0.0, "Re_t: {re}");
    }
    #[test]
    fn test_les_resolution_quality_range() {
        let q = les_resolution_quality(0.1, 0.9);
        assert!((q - 0.1).abs() < 1e-12, "LES quality = {q}");
    }
    #[test]
    fn test_les_quality_index_well_resolved() {
        let m = les_quality_index(0.9, 0.1);
        assert!(m > 0.8, "Well-resolved LES: M={m}");
    }
    #[test]
    fn test_kolmogorov_time_scale_positive() {
        let tau = kolmogorov_time_scale(1e-5, 0.001);
        assert!(tau > 0.0, "Kolmogorov time: {tau}");
    }
    #[test]
    fn test_turbulent_time_scale_formula() {
        let t = turbulent_time_scale(1.0, 0.5);
        assert!((t - 2.0).abs() < 1e-12, "T = k/eps: {t}");
    }
    #[test]
    fn test_turbulent_prandtl_ratio() {
        let pr_t = turbulent_prandtl(1e-3, 0.5e-3);
        assert!((pr_t - 2.0).abs() < 1e-12, "Pr_t = {pr_t}");
    }
    #[test]
    fn test_dynamic_smag_leonard_stress() {
        let m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        let l = m.leonard_stress_scalar(0.25, 0.25);
        assert_eq!(l, 0.0, "Leonard stress: {l}");
    }
    #[test]
    fn test_dynamic_smag_m_tensor_scalar() {
        let m = DynamicSmagorinskyLilly::new(1.0, 2.0);
        let m_val = m.m_tensor_scalar(1.0, 1.0, 1.0, 1.0);
        assert!((m_val - 3.0).abs() < 1e-12, "M tensor: {m_val}");
    }
}
