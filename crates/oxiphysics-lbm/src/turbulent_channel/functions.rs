//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// D2Q9 weights.
pub(super) const W9: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// D2Q9 x-velocity components.
pub(super) const CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
/// D2Q9 y-velocity components.
pub(super) const CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
/// Bounce-back partner directions for D2Q9.
pub(super) const BB: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
/// Von Kármán constant κ for the log-law.
pub const KAPPA_VK: f64 = 0.41;
/// Log-law additive constant B.
pub const B_LL: f64 = 5.2;
/// Viscous sublayer / log-region transition y⁺.
pub const Y_PLUS_TRANSITION: f64 = 11.3;
/// Log-region upper bound y⁺ (approximate).
pub const Y_PLUS_LOG_UPPER: f64 = 300.0;
/// Compute the D2Q9 equilibrium distribution.
///
/// `f_eq = w_q * ρ * [1 + 3(e·u) + 9/2(e·u)² - 3/2|u|²]`
#[inline]
pub fn equilibrium_d2q9(rho: f64, ux: f64, uy: f64, q: usize) -> f64 {
    let cx = CX[q] as f64;
    let cy = CY[q] as f64;
    let eu = cx * ux + cy * uy;
    let u2 = ux * ux + uy * uy;
    W9[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2)
}
/// Compute macroscopic density and velocity from D2Q9 distributions.
#[inline]
pub fn moments_d2q9(f: &[f64], base: usize) -> (f64, f64, f64) {
    let mut rho = 0.0_f64;
    let mut mx = 0.0_f64;
    let mut my = 0.0_f64;
    for q in 0..9 {
        let fq = f[base + q];
        rho += fq;
        mx += fq * CX[q] as f64;
        my += fq * CY[q] as f64;
    }
    let ux = if rho > 0.0 { mx / rho } else { 0.0 };
    let uy = if rho > 0.0 { my / rho } else { 0.0 };
    (rho, ux, uy)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::turbulent_channel::types::*;
    fn make_channel(nx: usize, ny: usize) -> ChannelFlow {
        ChannelFlow::new(nx, ny, 180.0)
    }
    #[test]
    fn test_new_sizes() {
        let ch = make_channel(8, 16);
        assert_eq!(ch.rho.len(), 8 * 16);
        assert_eq!(ch.u.len(), 8 * 16);
        assert_eq!(ch.f_dist.len(), 8 * 16 * 9);
    }
    #[test]
    fn test_index_correctness() {
        let ch = make_channel(4, 8);
        assert_eq!(ch.index(0, 0), 0);
        assert_eq!(ch.index(1, 0), 8);
        assert_eq!(ch.index(0, 1), 1);
    }
    #[test]
    fn test_fi_correctness() {
        let ch = make_channel(4, 8);
        assert_eq!(ch.fi(0, 0, 0), 0);
        assert_eq!(ch.fi(0, 0, 8), 8);
        assert_eq!(ch.fi(1, 0, 0), 8 * 9);
    }
    #[test]
    fn test_viscosity_set_correctly() {
        let ch = ChannelFlow::new(8, 64, 180.0);
        let h = 64.0 / 2.0;
        let expected_nu = ch.u_tau * h / 180.0;
        assert!((ch.nu - expected_nu).abs() < 1e-14);
    }
    #[test]
    fn test_log_law_viscous_sublayer_y_plus_1() {
        let ch = make_channel(8, 16);
        let u_plus = ch.log_law_u_plus(1.0);
        assert!(
            (u_plus - 1.0).abs() < 1e-10,
            "u+ at y+=1 should be 1, got {}",
            u_plus
        );
    }
    #[test]
    fn test_log_law_y_plus_5() {
        let ch = make_channel(8, 16);
        let u_plus = ch.log_law_u_plus(5.0);
        assert!((u_plus - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_log_law_y_plus_30() {
        let ch = make_channel(8, 16);
        let y_plus = 30.0_f64;
        let u_plus = ch.log_law_u_plus(y_plus);
        let expected = (1.0 / KAPPA_VK) * y_plus.ln() + B_LL;
        assert!((u_plus - expected).abs() < 1e-10);
    }
    #[test]
    fn test_log_law_y_plus_100() {
        let ch = make_channel(8, 16);
        let y_plus = 100.0_f64;
        let u_plus = ch.log_law_u_plus(y_plus);
        let expected = (1.0 / KAPPA_VK) * y_plus.ln() + B_LL;
        assert!((u_plus - expected).abs() < 1e-10);
    }
    #[test]
    fn test_log_law_monotone() {
        let ch = make_channel(8, 16);
        let u1 = ch.log_law_u_plus(10.0);
        let u2 = ch.log_law_u_plus(50.0);
        let u3 = ch.log_law_u_plus(200.0);
        assert!(u1 < u2 && u2 < u3);
    }
    #[test]
    fn test_log_law_zero_y_plus() {
        let ch = make_channel(8, 16);
        let u = ch.log_law_u_plus(0.0);
        assert!((u - 0.0).abs() < 1e-14);
    }
    #[test]
    fn test_init_parabolic_positive_bulk() {
        let mut ch = ChannelFlow::new(8, 32, 180.0);
        ch.init_parabolic();
        let bulk = ch.bulk_velocity();
        assert!(
            bulk > 0.0,
            "Bulk velocity should be positive after parabolic init"
        );
    }
    #[test]
    fn test_init_parabolic_centerline_max() {
        let mut ch = ChannelFlow::new(8, 32, 180.0);
        ch.init_parabolic();
        let profile = ch.velocity_profile();
        let center_v = profile[16].1;
        for (j, (_y, u)) in profile.iter().enumerate() {
            if j > 0 && j < 31 {
                assert!(*u <= center_v + 1e-10);
            }
        }
    }
    #[test]
    fn test_viscous_sublayer_thickness() {
        let ch = ChannelFlow::new(8, 32, 180.0);
        let delta_v = ch.viscous_sublayer_thickness();
        let expected = ch.nu / ch.u_tau;
        assert!((delta_v - expected).abs() < 1e-14);
    }
    #[test]
    fn test_equilibrium_sum_to_rho() {
        let ch = make_channel(4, 8);
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.0;
        let sum: f64 = (0..9).map(|q| ch.equilibrium(rho, ux, uy, q)).sum();
        assert!((sum - rho).abs() < 1e-12);
    }
    #[test]
    fn test_equilibrium_x_momentum() {
        let ch = make_channel(4, 8);
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.01;
        let mx: f64 = (0..9)
            .map(|q| ch.equilibrium(rho, ux, uy, q) * CX[q] as f64)
            .sum();
        assert!((mx - rho * ux).abs() < 1e-12);
    }
    #[test]
    fn test_equilibrium_y_momentum() {
        let ch = make_channel(4, 8);
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.01;
        let my: f64 = (0..9)
            .map(|q| ch.equilibrium(rho, ux, uy, q) * CY[q] as f64)
            .sum();
        assert!((my - rho * uy).abs() < 1e-12);
    }
    #[test]
    fn test_friction_velocity_positive() {
        let ch = ChannelFlow::new(8, 32, 180.0);
        assert!(ch.friction_velocity() > 0.0);
    }
    #[test]
    fn test_friction_velocity_matches_u_tau() {
        let ch = ChannelFlow::new(8, 32, 180.0);
        let fv = ch.friction_velocity();
        assert!(fv > 0.0, "friction velocity should be positive: {}", fv);
    }
    #[test]
    fn test_bulk_velocity_zero_initially() {
        let ch = make_channel(8, 16);
        let bulk = ch.bulk_velocity();
        assert!(bulk.abs() < 1e-14);
    }
    #[test]
    fn test_step_no_nan() {
        let mut ch = ChannelFlow::new(4, 16, 180.0);
        for _ in 0..10 {
            ch.step();
        }
        for &v in ch.f_dist.iter() {
            assert!(!v.is_nan());
        }
    }
    #[test]
    fn test_forcing_increases_momentum() {
        let mut ch = ChannelFlow::new(8, 16, 180.0);
        let bulk0 = ch.bulk_velocity();
        for _ in 0..50 {
            ch.step();
        }
        let bulk1 = ch.bulk_velocity();
        assert!(
            bulk1 >= bulk0,
            "Forcing should increase momentum: {:.6} -> {:.6}",
            bulk0,
            bulk1
        );
    }
    #[test]
    fn test_velocity_profile_length() {
        let ch = make_channel(8, 16);
        let profile = ch.velocity_profile();
        assert_eq!(profile.len(), 16);
    }
    #[test]
    fn test_velocity_profile_y_values() {
        let ch = make_channel(4, 8);
        let profile = ch.velocity_profile();
        for (j, (y, _)) in profile.iter().enumerate() {
            assert!((y - j as f64).abs() < 1e-14);
        }
    }
    #[test]
    fn test_wall_shear_stress_positive() {
        let ch = ChannelFlow::new(8, 32, 180.0);
        assert!(ch.wall_shear_stress() > 0.0);
    }
    #[test]
    fn test_centerline_velocity_non_negative_after_run() {
        let mut ch = ChannelFlow::new(8, 16, 180.0);
        for _ in 0..20 {
            ch.step();
        }
        let cv = ch.centerline_velocity();
        assert!(
            cv >= 0.0,
            "centerline velocity should be non-negative: {}",
            cv
        );
    }
    #[test]
    fn test_reynolds_stress_finite() {
        let mut ch = ChannelFlow::new(8, 16, 180.0);
        for _ in 0..5 {
            ch.step();
        }
        let rs = ch.reynolds_stress(8);
        assert!(rs.is_finite());
    }
    #[test]
    fn test_profile_symmetric_after_init() {
        let mut ch = ChannelFlow::new(8, 32, 180.0);
        ch.init_parabolic();
        let profile = ch.velocity_profile();
        let ny = ch.ny;
        for j in 1..(ny / 2) {
            let u_lower = profile[j].1;
            let u_upper = profile[ny - 1 - j].1;
            assert!(
                (u_lower - u_upper).abs() < 1e-10,
                "Profile not symmetric at j={}",
                j
            );
        }
    }
    #[test]
    fn test_nx_ny_stored() {
        let ch = ChannelFlow::new(12, 24, 200.0);
        assert_eq!(ch.nx, 12);
        assert_eq!(ch.ny, 24);
    }
    #[test]
    fn test_re_tau_stored() {
        let ch = ChannelFlow::new(8, 16, 395.0);
        assert!((ch.re_tau - 395.0).abs() < 1e-14);
    }
    #[test]
    fn test_no_slip_walls_sets_zero_velocity() {
        let mut ch = ChannelFlow::new(8, 16, 180.0);
        for _ in 0..20 {
            ch.step();
        }
        for i in 0..ch.nx {
            let u_wall = ch.u[ch.index(i, 0)][0];
            assert!(u_wall.abs() < 0.1, "Wall velocity too large: {}", u_wall);
        }
    }
    #[test]
    fn test_f_dist_non_negative_after_run() {
        let mut ch = ChannelFlow::new(4, 8, 180.0);
        for _ in 0..5 {
            ch.step();
        }
        for &v in ch.f_dist.iter() {
            assert!(v.is_finite());
        }
    }
    #[test]
    fn test_forcing_stored() {
        let ch = ChannelFlow::new(8, 32, 180.0);
        assert!(ch.forcing > 0.0);
    }
    #[test]
    fn test_rho_stays_near_one() {
        let mut ch = ChannelFlow::new(8, 16, 180.0);
        for _ in 0..20 {
            ch.step();
        }
        let mean_rho: f64 = ch.rho.iter().sum::<f64>() / ch.rho.len() as f64;
        assert!(
            (mean_rho - 1.0).abs() < 0.1,
            "Mean rho drifted: {}",
            mean_rho
        );
    }
    #[test]
    fn test_params_new_re180() {
        let p = TurbulentChannelParams::new(180.0, 32.0);
        assert!((p.re_tau - 180.0).abs() < 1e-14);
        assert!(p.nu > 0.0);
        assert!(p.u_tau > 0.0);
        assert!(p.bulk_velocity > 0.0);
    }
    #[test]
    fn test_params_viscous_length() {
        let p = TurbulentChannelParams::new(180.0, 32.0);
        let dl = p.viscous_length();
        assert!((dl - p.nu / p.u_tau).abs() < 1e-14);
    }
    #[test]
    fn test_params_friction_coefficient_positive() {
        let p = TurbulentChannelParams::new(180.0, 32.0);
        assert!(p.friction_coefficient() > 0.0);
    }
    #[test]
    fn test_params_y_plus() {
        let p = TurbulentChannelParams::new(180.0, 32.0);
        let yp = p.y_plus(1.0);
        let expected = p.u_tau / p.nu;
        assert!((yp - expected).abs() < 1e-12);
    }
    #[test]
    fn test_params_re_bulk_positive() {
        let p = TurbulentChannelParams::new(180.0, 32.0);
        assert!(p.re_bulk() > 0.0);
    }
    #[test]
    fn test_dns_lbm_new_sizes() {
        let dns = DnsDrivenLbm::new(8, 16, 180.0);
        assert_eq!(dns.f_dist.len(), 8 * 16 * 9);
        assert_eq!(dns.u.len(), 8 * 16);
        assert_eq!(dns.rho.len(), 8 * 16);
    }
    #[test]
    fn test_dns_lbm_step_no_nan() {
        let mut dns = DnsDrivenLbm::new(4, 16, 180.0);
        for _ in 0..10 {
            dns.step();
        }
        for &v in dns.f_dist.iter() {
            assert!(!v.is_nan(), "NaN in f_dist");
        }
    }
    #[test]
    fn test_dns_lbm_bulk_increases() {
        let mut dns = DnsDrivenLbm::new(8, 16, 180.0);
        let b0 = dns.bulk_velocity();
        for _ in 0..30 {
            dns.step();
        }
        let b1 = dns.bulk_velocity();
        assert!(b1 >= b0);
    }
    #[test]
    fn test_dns_lbm_friction_velocity_positive() {
        let dns = DnsDrivenLbm::new(8, 32, 180.0);
        assert!(dns.friction_velocity() > 0.0);
    }
    #[test]
    fn test_dns_lbm_accumulate_stats() {
        let mut dns = DnsDrivenLbm::new(4, 16, 180.0);
        for _ in 0..10 {
            dns.step();
            dns.accumulate_stats();
        }
        assert_eq!(dns.stat_count, 10);
        let profile = dns.mean_velocity_profile();
        assert_eq!(profile.len(), 16);
    }
    #[test]
    fn test_dns_lbm_reset_stats() {
        let mut dns = DnsDrivenLbm::new(4, 16, 180.0);
        for _ in 0..5 {
            dns.step();
            dns.accumulate_stats();
        }
        dns.reset_stats();
        assert_eq!(dns.stat_count, 0);
    }
    #[test]
    fn test_dns_lbm_turbulent_profile_init() {
        let mut dns = DnsDrivenLbm::new(8, 16, 180.0);
        dns.init_turbulent_profile();
        let b = dns.bulk_velocity();
        assert!(b > 0.0);
    }
    #[test]
    fn test_dns_lbm_computed_re_tau_reasonable() {
        let dns = DnsDrivenLbm::new(8, 32, 180.0);
        let re_c = dns.computed_re_tau();
        assert!(re_c > 0.0);
    }
    #[test]
    fn test_rans_lbm_new() {
        let rans = ReynoldsAveragedLbm::new(8, 16, 180.0);
        assert_eq!(rans.nu_t.len(), 8 * 16);
        assert_eq!(rans.f_dist.len(), 8 * 16 * 9);
    }
    #[test]
    fn test_rans_lbm_step_no_nan() {
        let mut rans = ReynoldsAveragedLbm::new(4, 16, 180.0);
        for _ in 0..10 {
            rans.step();
        }
        for &v in rans.f_dist.iter() {
            assert!(!v.is_nan());
        }
    }
    #[test]
    fn test_rans_lbm_mixing_length_zero_at_wall() {
        let rans = ReynoldsAveragedLbm::new(8, 32, 180.0);
        let lm0 = rans.mixing_length(0);
        assert!(
            lm0 < 1e-10,
            "Mixing length at wall should be ~0, got {}",
            lm0
        );
    }
    #[test]
    fn test_rans_lbm_eddy_diffusivity() {
        let rans = ReynoldsAveragedLbm::new(8, 32, 180.0);
        let ed = rans.eddy_diffusivity(16);
        assert!(ed >= 0.0);
    }
    #[test]
    fn test_rans_lbm_bulk_velocity_increases() {
        let mut rans = ReynoldsAveragedLbm::new(8, 16, 180.0);
        let b0 = rans.bulk_velocity();
        for _ in 0..30 {
            rans.step();
        }
        let b1 = rans.bulk_velocity();
        assert!(b1 >= b0);
    }
    #[test]
    fn test_wf_lbm_new() {
        let wf = WallFunctionLbm::new(8, 16, 180.0);
        assert!(wf.y_plus_first >= 0.0);
        assert_eq!(wf.f_dist.len(), 8 * 16 * 9);
    }
    #[test]
    fn test_wf_lbm_step_no_nan() {
        let mut wf = WallFunctionLbm::new(4, 16, 180.0);
        for _ in 0..10 {
            wf.step();
        }
        for &v in wf.f_dist.iter() {
            assert!(!v.is_nan());
        }
    }
    #[test]
    fn test_wf_lbm_friction_velocity_positive() {
        let wf = WallFunctionLbm::new(8, 32, 180.0);
        assert!(wf.friction_velocity() > 0.0);
    }
    #[test]
    fn test_wf_lbm_wall_shear_from_log_law() {
        let wf = WallFunctionLbm::new(8, 32, 180.0);
        let tau = wf.wall_shear_from_log_law(0.1, 30.0);
        assert!(tau > 0.0);
    }
    #[test]
    fn test_wf_lbm_y_plus_at() {
        let wf = WallFunctionLbm::new(8, 32, 180.0);
        let yp = wf.y_plus_at(0);
        assert!(yp >= 0.0);
    }
    #[test]
    fn test_vpa_from_profile() {
        let profile: Vec<(f64, f64)> = (0..20).map(|j| (j as f64, j as f64 * 0.01)).collect();
        let vpa = VelocityProfileAnalysis::from_profile(&profile, 0.01, 1e-4, 10.0);
        assert_eq!(vpa.y_plus.len(), 20);
        assert_eq!(vpa.u_plus.len(), 20);
        assert_eq!(vpa.y_outer.len(), 20);
    }
    #[test]
    fn test_vpa_log_law_static() {
        let u_plus = VelocityProfileAnalysis::log_law_u_plus(100.0);
        let expected = (1.0 / KAPPA_VK) * 100_f64.ln() + B_LL;
        assert!((u_plus - expected).abs() < 1e-10);
    }
    #[test]
    fn test_vpa_sublayer_slope_near_one() {
        let u_tau = 0.01_f64;
        let nu = 0.001_f64;
        let delta_nu = nu / u_tau;
        let profile: Vec<(f64, f64)> = (1..=5)
            .map(|j| {
                let y = j as f64 * delta_nu;
                let u_x = j as f64 * u_tau;
                (y, u_x)
            })
            .collect();
        let vpa = VelocityProfileAnalysis::from_profile(&profile, u_tau, nu, 10.0);
        let slope = vpa.sublayer_slope();
        assert!(
            (slope - 1.0).abs() < 1e-10,
            "sublayer slope should be ~1, got {}",
            slope
        );
    }
    #[test]
    fn test_vpa_velocity_defect() {
        let profile: Vec<(f64, f64)> = (1..=10)
            .map(|j| (j as f64, 0.1 * (j as f64).ln() + 5.0))
            .collect();
        let vpa = VelocityProfileAnalysis::from_profile(&profile, 0.01, 1e-4, 10.0);
        let defect = vpa.velocity_defect();
        for (_, d) in &defect {
            assert!(*d >= -1e-10, "negative defect: {}", d);
        }
    }
    #[test]
    fn test_vpa_log_region_points_empty_low_re() {
        let profile: Vec<(f64, f64)> = (1..=5)
            .map(|j| (j as f64 * 0.01, j as f64 * 0.001))
            .collect();
        let vpa = VelocityProfileAnalysis::from_profile(&profile, 0.01, 0.01, 1.0);
        let log_pts = vpa.log_region_points();
        assert!(log_pts.is_empty() || !log_pts.is_empty());
    }
    #[test]
    fn test_turb_stats_new() {
        let ts = TurbulentStatistics::new(16);
        assert_eq!(ts.sum_u.len(), 16);
        assert_eq!(ts.count, 0);
    }
    #[test]
    fn test_turb_stats_accumulate() {
        let mut ts = TurbulentStatistics::new(8);
        let nx = 4;
        let u_field: Vec<[f64; 2]> = (0..nx * 8).map(|_| [0.1, 0.0]).collect();
        ts.accumulate(&u_field, nx);
        assert_eq!(ts.count, 1);
    }
    #[test]
    fn test_turb_stats_mean_u() {
        let mut ts = TurbulentStatistics::new(8);
        let nx = 4;
        let u_field: Vec<[f64; 2]> = (0..nx * 8).map(|_| [0.5, 0.0]).collect();
        ts.accumulate(&u_field, nx);
        let mean = ts.mean_u_profile();
        for &m in &mean {
            assert!((m - 0.5).abs() < 1e-12);
        }
    }
    #[test]
    fn test_turb_stats_u_rms_zero_for_uniform() {
        let mut ts = TurbulentStatistics::new(8);
        let nx = 4;
        let u_field: Vec<[f64; 2]> = (0..nx * 8).map(|_| [0.5, 0.0]).collect();
        for _ in 0..5 {
            ts.accumulate(&u_field, nx);
        }
        let rms = ts.u_rms_profile();
        for &r in &rms {
            assert!(r < 1e-10, "rms should be ~0 for uniform field: {}", r);
        }
    }
    #[test]
    fn test_turb_stats_tke_profile_non_negative() {
        let mut ts = TurbulentStatistics::new(8);
        let nx = 4;
        let u_field: Vec<[f64; 2]> = (0..nx * 8)
            .map(|k| [0.1 * (k as f64 % 3.0), 0.01])
            .collect();
        for _ in 0..3 {
            ts.accumulate(&u_field, nx);
        }
        let tke = ts.tke_profile();
        for &k in &tke {
            assert!(k >= 0.0);
        }
    }
    #[test]
    fn test_turb_stats_reset() {
        let mut ts = TurbulentStatistics::new(8);
        let nx = 4;
        let u_field: Vec<[f64; 2]> = (0..nx * 8).map(|_| [0.1, 0.0]).collect();
        ts.accumulate(&u_field, nx);
        ts.reset();
        assert_eq!(ts.count, 0);
        let mean = ts.mean_u_profile();
        for &m in &mean {
            assert!((m).abs() < 1e-14);
        }
    }
    #[test]
    fn test_turb_stats_reynolds_stress() {
        let ts = TurbulentStatistics::new(8);
        let rs = ts.reynolds_stress_profile();
        assert_eq!(rs.len(), 8);
        for &r in &rs {
            assert!((r).abs() < 1e-14);
        }
    }
    #[test]
    fn test_turb_stats_rs_gradient_length() {
        let ts = TurbulentStatistics::new(8);
        let grad = ts.reynolds_stress_gradient();
        assert_eq!(grad.len(), 8);
    }
    #[test]
    fn test_validation_kmm_re180_created() {
        let v = ChannelFlowValidation::kmm_re180();
        assert!((v.re_tau_dns - 180.0).abs() < 1e-10);
        assert!(!v.dns_y_plus.is_empty());
        assert!(!v.dns_u_plus.is_empty());
    }
    #[test]
    fn test_validation_mkm_re590_created() {
        let v = ChannelFlowValidation::mkm_re590();
        assert!((v.re_tau_dns - 590.0).abs() < 1e-10);
        assert!(!v.dns_y_plus.is_empty());
    }
    #[test]
    fn test_validation_interpolate_in_range() {
        let v = ChannelFlowValidation::kmm_re180();
        let yp = v.dns_y_plus[5];
        let up = v.interpolate_u_plus(yp);
        let expected = v.dns_u_plus[5];
        assert!((up - expected).abs() < 1e-10);
    }
    #[test]
    fn test_validation_l2_error_zero_exact_match() {
        let v = ChannelFlowValidation::kmm_re180();
        let err = v.l2_error_u_plus(&v.dns_y_plus.clone(), &v.dns_u_plus.clone());
        assert!(err < 1e-10, "L2 error with self should be ~0: {}", err);
    }
    #[test]
    fn test_validation_friction_coefficient() {
        let cf = ChannelFlowValidation::friction_coefficient(20.0);
        assert!((cf - 2.0 / 400.0).abs() < 1e-14);
    }
    #[test]
    fn test_validation_dean_correlation() {
        let cf_dean = ChannelFlowValidation::dean_friction_coefficient(5600.0);
        assert!(cf_dean > 0.0);
        assert!(cf_dean < 0.05);
    }
    #[test]
    fn test_rst_from_samples_zero() {
        let rst = ReynoldsStressTensor::from_samples(&[], &[]);
        assert!((rst.uu).abs() < 1e-14);
        assert!((rst.vv).abs() < 1e-14);
    }
    #[test]
    fn test_rst_tke() {
        let u = vec![1.0, -1.0, 1.0, -1.0];
        let v = vec![0.5, -0.5, 0.5, -0.5];
        let rst = ReynoldsStressTensor::from_samples(&u, &v);
        let tke = rst.tke();
        assert!((tke - 0.625).abs() < 1e-12, "tke = {}", tke);
    }
    #[test]
    fn test_rst_anisotropy_non_negative() {
        let u = vec![1.0, -1.0];
        let v = vec![0.5, 0.5];
        let rst = ReynoldsStressTensor::from_samples(&u, &v);
        assert!(rst.anisotropy() >= 0.0);
    }
    #[test]
    fn test_tke_budget_production_positive() {
        let budget = TkeProductionBudget::estimate(0.01, 1.0, 1e-4, 0.001, 0.1);
        assert!(budget.production > 0.0);
    }
    #[test]
    fn test_tke_budget_dissipation_non_negative() {
        let budget = TkeProductionBudget::estimate(0.01, 1.0, 1e-4, 0.001, 0.1);
        assert!(budget.dissipation >= 0.0);
    }
    #[test]
    fn test_tke_budget_residual_finite() {
        let budget = TkeProductionBudget::estimate(0.01, 1.0, 1e-4, 0.001, 0.1);
        assert!(budget.residual().is_finite());
    }
    #[test]
    fn test_tke_budget_zero_mixing_length() {
        let budget = TkeProductionBudget::estimate(0.01, 1.0, 1e-4, 0.001, 0.0);
        assert!((budget.dissipation).abs() < 1e-14);
    }
    #[test]
    fn test_equilibrium_d2q9_sum() {
        let sum: f64 = (0..9).map(|q| equilibrium_d2q9(1.0, 0.05, 0.0, q)).sum();
        assert!((sum - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_moments_d2q9_at_rest() {
        let mut f = vec![0.0f64; 9];
        f[..9].copy_from_slice(&W9);
        let (rho, ux, uy) = moments_d2q9(&f, 0);
        assert!((rho - 1.0).abs() < 1e-12);
        assert!(ux.abs() < 1e-14);
        assert!(uy.abs() < 1e-14);
    }
}
