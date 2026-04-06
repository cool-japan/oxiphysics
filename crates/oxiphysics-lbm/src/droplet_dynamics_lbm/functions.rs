//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::manual_range_contains)]
/// Ranz-Marshall Sherwood number correlation.
///
/// Sh = 2 + 0.6 Re^0.5 Sc^(1/3)
pub fn sherwood_number(reynolds: f64, schmidt: f64) -> f64 {
    2.0 + 0.6 * reynolds.sqrt() * schmidt.powf(1.0 / 3.0)
}
/// Probit function (inverse normal CDF) — Rational approximation (Abramowitz & Stegun).
pub(super) fn probit(p: f64) -> f64 {
    let p = p.clamp(1e-10, 1.0 - 1e-10);
    let c = [2.515517, 0.802853, 0.010328];
    let d = [1.432788, 0.189269, 0.001308];
    let sign = if p >= 0.5 { 1.0 } else { -1.0 };
    let q = if p >= 0.5 { 1.0 - p } else { p };
    let t = (-2.0 * q.ln()).sqrt();
    let num = c[0] + c[1] * t + c[2] * t * t;
    let den = 1.0 + d[0] * t + d[1] * t * t + d[2] * t * t * t;
    sign * (t - num / den)
}
/// Laplace pressure jump ΔP = 2σ/R for a spherical droplet.
pub fn laplace_pressure(radius: f64, surface_tension: f64) -> f64 {
    if radius <= 0.0 {
        0.0
    } else {
        2.0 * surface_tension / radius
    }
}
/// Ohnesorge number Oh = μ / sqrt(ρ σ R).
pub fn ohnesorge_number(viscosity: f64, density: f64, surface_tension: f64, radius: f64) -> f64 {
    let denom = (density * surface_tension * radius).sqrt();
    if denom == 0.0 {
        f64::INFINITY
    } else {
        viscosity / denom
    }
}
/// Morton number Mo = g μ⁴ / (ρ σ³).
pub fn morton_number(gravity: f64, viscosity: f64, density: f64, surface_tension: f64) -> f64 {
    if surface_tension == 0.0 || density == 0.0 {
        return f64::INFINITY;
    }
    gravity * viscosity.powi(4) / (density * surface_tension.powi(3))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::droplet_dynamics_lbm::types::*;
    use std::f64::consts::PI;
    #[test]
    fn test_droplet_volume() {
        let d = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let expected = 4.0 / 3.0 * PI;
        assert!((d.volume() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_droplet_speed() {
        let d = Droplet::new([0.0; 3], 1.0, [3.0, 4.0, 0.0], 0.1);
        assert!((d.speed() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_droplet_weber_number() {
        let d = Droplet::new([0.0; 3], 2.0, [1.0, 0.0, 0.0], 1.0);
        assert!((d.weber_number(1.0) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_droplet_weber_zero_sigma() {
        let d = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 0.0);
        assert!(d.weber_number(1.0).is_infinite());
    }
    #[test]
    fn test_droplet_capillary_number() {
        let d = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 1.0);
        assert!((d.capillary_number(0.5) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_droplet_distance() {
        let a = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([3.0, 4.0, 0.0], 1.0, [0.0; 3], 0.1);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_droplet_surface_gap_positive() {
        let a = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([5.0, 0.0, 0.0], 1.0, [0.0; 3], 0.1);
        assert!((a.surface_gap(&b) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_droplet_surface_gap_negative() {
        let a = Droplet::new([0.0; 3], 2.0, [0.0; 3], 0.1);
        let b = Droplet::new([1.0, 0.0, 0.0], 2.0, [0.0; 3], 0.1);
        assert!(a.surface_gap(&b) < 0.0);
    }
    #[test]
    fn test_vdw_pressure_positive_gap() {
        let crit = CoalescenceCriterion::new(0.5, 10.0, 1e-19);
        let p = crit.van_der_waals_pressure(1e-9);
        assert!(p < 0.0);
    }
    #[test]
    fn test_vdw_pressure_zero_gap() {
        let crit = CoalescenceCriterion::new(0.5, 10.0, 1e-19);
        assert_eq!(crit.van_der_waals_pressure(0.0), 0.0);
    }
    #[test]
    fn test_coalescence_should_coalesce_true() {
        let crit = CoalescenceCriterion::new(1.0, 10.0, 1e-19);
        let a = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([2.5, 0.0, 0.0], 1.0, [0.0; 3], 0.1);
        assert!(crit.should_coalesce(&a, &b));
    }
    #[test]
    fn test_coalescence_should_not_coalesce() {
        let crit = CoalescenceCriterion::new(0.1, 10.0, 1e-19);
        let a = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.1);
        let b = Droplet::new([10.0, 0.0, 0.0], 1.0, [0.0; 3], 0.1);
        assert!(!crit.should_coalesce(&a, &b));
    }
    #[test]
    fn test_coalescence_merge_volume_conservation() {
        let crit = CoalescenceCriterion::new(1.0, 10.0, 1e-19);
        let a = Droplet::new([0.0; 3], 1.0, [1.0, 0.0, 0.0], 0.1);
        let b = Droplet::new([3.0, 0.0, 0.0], 2.0, [0.0; 3], 0.1);
        let merged = crit.merge(&a, &b);
        let v_total = a.volume() + b.volume();
        assert!((merged.volume() - v_total).abs() < 1e-8);
    }
    #[test]
    fn test_coalescence_effective_drain_time() {
        let crit = CoalescenceCriterion::new(0.5, 10.0, 0.0);
        assert!((crit.effective_drain_time(0.5) - 15.0).abs() < 1e-10);
    }
    #[test]
    fn test_breakup_weber_criterion_true() {
        let bu = DropletBreakup::new(1.0, 0.1, 9.81);
        let d = Droplet::new([0.0; 3], 1.0, [10.0, 0.0, 0.0], 0.01);
        assert!(bu.weber_breakup(&d, 1.0));
    }
    #[test]
    fn test_breakup_weber_criterion_false() {
        let bu = DropletBreakup::new(1000.0, 0.1, 9.81);
        let d = Droplet::new([0.0; 3], 1.0, [0.1, 0.0, 0.0], 1.0);
        assert!(!bu.weber_breakup(&d, 1.0));
    }
    #[test]
    fn test_breakup_child_radius() {
        let r_child = DropletBreakup::child_radius(2.0_f64.cbrt());
        assert!((r_child - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_breakup_rt_growth_zero_surface_tension() {
        let bu = DropletBreakup::new(12.0, 0.5, 9.81);
        let d = Droplet::new([0.0; 3], 1.0, [0.0; 3], 0.0);
        assert_eq!(bu.rt_growth_rate(&d, 1000.0), 0.0);
    }
    #[test]
    fn test_breakup_rt_positive_with_atwood() {
        let bu = DropletBreakup::new(12.0, 0.5, 100.0);
        let d = Droplet::new([0.0; 3], 10.0, [0.0; 3], 0.001);
        let rate = bu.rt_growth_rate(&d, 1000.0);
        assert!(rate >= 0.0);
    }
    #[test]
    fn test_schiller_naumann_low_re() {
        let lf = LiftDragDroplet::new(1.0, 0.001, [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        let cd = lf.schiller_naumann_cd(0.0005);
        assert!(cd > 20.0);
    }
    #[test]
    fn test_schiller_naumann_high_re() {
        let lf = LiftDragDroplet::new(1000.0, 0.001, [2.0, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        let cd = lf.schiller_naumann_cd(0.01);
        assert!((cd - 0.44).abs() < 1e-10);
    }
    #[test]
    fn test_drag_force_direction() {
        let lf = LiftDragDroplet::new(1.0, 0.01, [1.0, 0.0, 0.0], [0.0; 3], [0.0; 3]);
        let f = lf.drag_force(0.1);
        assert!(f[0] > 0.0);
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn test_magnus_lift_nonzero() {
        let lf = LiftDragDroplet::new(1.0, 0.01, [1.0, 0.0, 0.0], [0.0; 3], [0.0, 0.0, 1.0]);
        let fm = lf.magnus_lift(0.1);
        assert!(fm[1].abs() > 0.0);
    }
    #[test]
    fn test_total_force_nonzero() {
        let lf = LiftDragDroplet::new(1.0, 0.01, [1.0, 0.5, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 0.1]);
        let ft = lf.total_force(0.1);
        let mag: f64 = ft.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!(mag > 0.0);
    }
    #[test]
    fn test_d2_at_time_zero() {
        let ev = EvaporationModel::new(4.0, 0.1, 1e-5, 0.01, 0.1, 800.0);
        assert!((ev.d2_at_time(0.0) - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_d2_at_time_partial() {
        let ev = EvaporationModel::new(4.0, 0.1, 1e-5, 0.01, 0.1, 800.0);
        assert!((ev.d2_at_time(10.0) - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_d2_at_time_clamped() {
        let ev = EvaporationModel::new(1.0, 1.0, 1e-5, 0.0, 0.1, 800.0);
        assert_eq!(ev.d2_at_time(1000.0), 0.0);
    }
    #[test]
    fn test_lifetime() {
        let ev = EvaporationModel::new(4.0, 0.5, 1e-5, 0.0, 0.1, 800.0);
        assert!((ev.lifetime() - 8.0).abs() < 1e-10);
    }
    #[test]
    fn test_spalding_number() {
        let ev = EvaporationModel::new(4.0, 0.1, 1e-5, 0.0, 0.5, 800.0);
        assert!((ev.spalding_number() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mass_transfer_rate_nonzero() {
        let ev = EvaporationModel::new(4.0, 0.1, 1e-5, 0.0, 0.5, 800.0);
        let rate = ev.mass_transfer_rate(0.1, 10.0, 0.7);
        assert!(rate <= 0.0);
    }
    #[test]
    fn test_vapor_pressure_positive() {
        let p = EvaporationModel::vapor_pressure(300.0, 2260.0, 373.0, 101325.0);
        assert!(p > 0.0);
    }
    #[test]
    fn test_sherwood_number_quiescent() {
        assert!((sherwood_number(0.0, 0.0) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_lbm_mass_conservation() {
        let mut lbm = DropletLBM::new(20, 1.0, 0.01, 5, 15);
        let mass_before = lbm.total_red_mass();
        for _ in 0..10 {
            lbm.step();
        }
        let mass_after = lbm.total_red_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "mass_before={mass_before} mass_after={mass_after}"
        );
    }
    #[test]
    fn test_lbm_phase_field_range() {
        let lbm = DropletLBM::new(10, 1.0, 0.01, 2, 7);
        for i in 0..10 {
            let phi = lbm.phase_field(i);
            assert!(phi >= -1.0 && phi <= 1.0, "phi[{i}]={phi}");
        }
    }
    #[test]
    fn test_lbm_density_positive() {
        let lbm = DropletLBM::new(10, 1.0, 0.01, 2, 7);
        for i in 0..10 {
            assert!(lbm.density(i) >= 0.0);
        }
    }
    #[test]
    fn test_log_normal_pdf_positive() {
        let spray = SprayModel::new(
            1.0,
            0.5,
            1.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::LogNormal,
        );
        assert!(spray.log_normal_pdf(1.0) > 0.0);
    }
    #[test]
    fn test_log_normal_pdf_zero_diameter() {
        let spray = SprayModel::new(
            1.0,
            0.5,
            1.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::LogNormal,
        );
        assert_eq!(spray.log_normal_pdf(0.0), 0.0);
    }
    #[test]
    fn test_rosin_rammler_cdf_zero() {
        let spray = SprayModel::new(
            1.0,
            2.0,
            1.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::RosinRammler,
        );
        assert_eq!(spray.rosin_rammler_cdf(0.0), 0.0);
    }
    #[test]
    fn test_rosin_rammler_cdf_approaches_one() {
        let spray = SprayModel::new(
            1.0,
            2.0,
            1.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::RosinRammler,
        );
        let cdf = spray.rosin_rammler_cdf(100.0);
        assert!(cdf > 0.999);
    }
    #[test]
    fn test_sample_diameter_positive() {
        let spray = SprayModel::new(
            50e-6,
            0.5,
            10.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::LogNormal,
        );
        for i in 0..10 {
            let d = spray.sample_diameter(i);
            assert!(d > 0.0, "d[{i}]={d}");
        }
    }
    #[test]
    fn test_sample_diameter_rr_positive() {
        let spray = SprayModel::new(
            50e-6,
            2.0,
            10.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::RosinRammler,
        );
        for i in 0..5 {
            let d = spray.sample_diameter(i);
            assert!(d > 0.0, "d[{i}]={d}");
        }
    }
    #[test]
    fn test_injection_velocity_nonzero() {
        let spray = SprayModel::new(
            50e-6,
            0.5,
            10.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::LogNormal,
        );
        let v = spray.injection_velocity(0, 0.1);
        let mag = (v[0].powi(2) + v[1].powi(2) + v[2].powi(2)).sqrt();
        assert!(mag > 0.0);
    }
    #[test]
    fn test_sauter_mean_diameter() {
        let spray = SprayModel::new(
            1.0,
            0.5,
            1.0,
            [0.0, 0.0, 1.0],
            100.0,
            SizeDistribution::LogNormal,
        );
        let smd = spray.sauter_mean_diameter();
        let mean = spray.log_normal_mean();
        assert!(smd >= mean, "SMD should be >= mean");
    }
    #[test]
    fn test_laplace_pressure() {
        assert!((laplace_pressure(1.0, 0.5) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_laplace_pressure_zero_radius() {
        assert_eq!(laplace_pressure(0.0, 0.5), 0.0);
    }
    #[test]
    fn test_ohnesorge_number_positive() {
        let oh = ohnesorge_number(0.001, 1000.0, 0.072, 1e-3);
        assert!(oh > 0.0);
    }
    #[test]
    fn test_morton_number_zero_sigma() {
        assert!(morton_number(9.81, 0.001, 1000.0, 0.0).is_infinite());
    }
    #[test]
    fn test_morton_number_positive() {
        let mo = morton_number(9.81, 0.001, 1000.0, 0.072);
        assert!(mo > 0.0);
    }
    #[test]
    fn test_probit_midpoint() {
        assert!(probit(0.5).abs() < 0.01);
    }
}
/// Union-find `find` with path compression.
pub(super) fn find(parent: &mut Vec<usize>, x: usize) -> usize {
    if parent[x] != x {
        parent[x] = find(parent, parent[x]);
    }
    parent[x]
}
/// Union-find `union` by label.
pub(super) fn union_uf(parent: &mut Vec<usize>, a: usize, b: usize) {
    let ra = find(parent, a);
    let rb = find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}
#[cfg(test)]
mod tests_new {

    use crate::droplet_dynamics_lbm::types::*;
    use std::f64::consts::PI;
    #[test]
    fn test_params_density_ratio() {
        let p = DropletLbmParams::new(2.0, 1.0, 0.01, 0.005, 0.1, 2.0);
        assert!((p.density_ratio() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_params_viscosity_ratio() {
        let p = DropletLbmParams::new(2.0, 1.0, 0.02, 0.01, 0.1, 2.0);
        assert!((p.viscosity_ratio() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_params_tau_droplet_positive() {
        let p = DropletLbmParams::new(1.0, 1.0, 0.1, 0.05, 0.01, 1.0);
        assert!(p.tau_droplet() > 0.5);
    }
    #[test]
    fn test_params_tau_carrier_positive() {
        let p = DropletLbmParams::new(1.0, 1.0, 0.1, 0.05, 0.01, 1.0);
        assert!(p.tau_carrier() > 0.5);
    }
    #[test]
    fn test_params_cahn_number() {
        let p = DropletLbmParams::new(1.0, 1.0, 0.01, 0.01, 0.01, 4.0);
        assert!((p.cahn_number(100.0) - 0.04).abs() < 1e-12);
    }
    #[test]
    fn test_params_capillary_number() {
        let p = DropletLbmParams::new(1.0, 1.0, 0.01, 0.01, 0.1, 2.0);
        assert!((p.capillary_number(1.0) - 0.1).abs() < 1e-12);
    }
    #[test]
    fn test_params_weber_number() {
        let p = DropletLbmParams::new(1.0, 2.0, 0.01, 0.01, 1.0, 2.0);
        assert!((p.weber_number(1.0, 1.0) - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_ch_lbm_initial_droplet_volume_fraction() {
        let ch = CahnHilliardLbm::new_circular_droplet(32, 32, 8.0, 1.0, 1e-4);
        let vf = ch.droplet_volume_fraction();
        assert!(vf > 0.1 && vf < 0.4, "volume fraction {vf}");
    }
    #[test]
    fn test_ch_lbm_total_order_parameter_finite() {
        let ch = CahnHilliardLbm::new_circular_droplet(16, 16, 4.0, 1.0, 1e-4);
        let phi_sum = ch.total_order_parameter();
        assert!(phi_sum.is_finite());
    }
    #[test]
    fn test_ch_lbm_step_no_panic() {
        let mut ch = CahnHilliardLbm::new_circular_droplet(16, 16, 4.0, 1.0, 1e-6);
        ch.step(0.1);
    }
    #[test]
    fn test_ch_lbm_bulk_free_energy_positive() {
        let ch = CahnHilliardLbm::new_circular_droplet(16, 16, 4.0, 1.0, 1e-6);
        let f = ch.bulk_free_energy();
        assert!(f >= 0.0, "free energy should be non-negative: {f}");
    }
    #[test]
    fn test_ch_lbm_idx_mapping() {
        let ch = CahnHilliardLbm::new_circular_droplet(8, 8, 2.0, 0.5, 1e-5);
        assert_eq!(ch.idx(0, 0), 0);
        assert_eq!(ch.idx(7, 0), 7);
        assert_eq!(ch.idx(0, 7), 56);
        assert_eq!(ch.idx(7, 7), 63);
    }
    #[test]
    fn test_ch_lbm_phi_range() {
        let ch = CahnHilliardLbm::new_circular_droplet(16, 16, 4.0, 1.0, 1e-5);
        for &phi in &ch.phi {
            assert!(phi >= -2.0 && phi <= 2.0, "phi out of range: {phi}");
        }
    }
    #[test]
    fn test_ch_lbm_compute_mu_no_nan() {
        let mut ch = CahnHilliardLbm::new_circular_droplet(16, 16, 4.0, 1.0, 1e-5);
        ch.compute_mu();
        for &mu in &ch.mu {
            assert!(mu.is_finite(), "mu NaN detected");
        }
    }
    #[test]
    fn test_tracking_single_droplet() {
        let mut phi = vec![-1.0_f64; 16];
        phi[5] = 1.0;
        phi[6] = 1.0;
        phi[9] = 1.0;
        phi[10] = 1.0;
        let dt = DropletTracking::label_droplets(&phi, 4, 4, 0.0);
        assert_eq!(dt.n_droplets, 1);
    }
    #[test]
    fn test_tracking_two_droplets() {
        let mut phi = vec![-1.0_f64; 50];
        phi[11] = 1.0;
        phi[12] = 1.0;
        phi[2 * 10 + 1] = 1.0;
        phi[2 * 10 + 2] = 1.0;
        phi[17] = 1.0;
        phi[18] = 1.0;
        phi[2 * 10 + 7] = 1.0;
        phi[2 * 10 + 8] = 1.0;
        let dt = DropletTracking::label_droplets(&phi, 10, 5, 0.0);
        assert_eq!(
            dt.n_droplets, 2,
            "should find 2 droplets, got {}",
            dt.n_droplets
        );
    }
    #[test]
    fn test_tracking_droplet_volume() {
        let mut phi = vec![-1.0_f64; 25];
        for j in 1..=3 {
            for i in 1..=3 {
                phi[j * 5 + i] = 1.0;
            }
        }
        let dt = DropletTracking::label_droplets(&phi, 5, 5, 0.0);
        assert_eq!(dt.n_droplets, 1);
        assert_eq!(dt.droplet_volume(1), 9);
    }
    #[test]
    fn test_tracking_centroid() {
        let mut phi = vec![-1.0_f64; 25];
        phi[2 * 5 + 2] = 1.0;
        let dt = DropletTracking::label_droplets(&phi, 5, 5, 0.0);
        let c = dt.centroid(1).unwrap();
        assert!((c[0] - 2.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_tracking_equivalent_radius_positive() {
        let mut phi = vec![-1.0_f64; 100];
        for j in 3..7 {
            for i in 3..7 {
                phi[j * 10 + i] = 1.0;
            }
        }
        let dt = DropletTracking::label_droplets(&phi, 10, 10, 0.0);
        let r = dt.equivalent_radius(1);
        assert!(r > 0.0, "radius should be positive: {r}");
    }
    #[test]
    fn test_tracking_circularity_range() {
        let mut phi = vec![-1.0_f64; 100];
        for j in 3..7 {
            for i in 3..7 {
                phi[j * 10 + i] = 1.0;
            }
        }
        let dt = DropletTracking::label_droplets(&phi, 10, 10, 0.0);
        let c = dt.circularity(1);
        assert!(
            c > 0.0 && c.is_finite(),
            "circularity should be positive and finite: {c}"
        );
    }
    #[test]
    fn test_coalescence_merged_radius_volume_conservation() {
        let cm = CoalescenceModel::new(0.01, 0.5, 0.1, 1e-20);
        let r1 = 2.0_f64;
        let r2 = 3.0_f64;
        let r_merged = cm.merged_radius(r1, r2);
        let v1 = 4.0 / 3.0 * PI * r1.powi(3);
        let v2 = 4.0 / 3.0 * PI * r2.powi(3);
        let vm = 4.0 / 3.0 * PI * r_merged.powi(3);
        assert!((vm - v1 - v2).abs() / (v1 + v2) < 1e-10);
    }
    #[test]
    fn test_coalescence_neck_growth_rate_inertial_positive() {
        let cm = CoalescenceModel::new(0.01, 0.5, 0.1, 1e-20);
        let rate = cm.neck_growth_rate_inertial(0.1, 1.0);
        assert!(
            rate > 0.0,
            "inertial neck growth should be positive: {rate}"
        );
    }
    #[test]
    fn test_coalescence_neck_growth_rate_viscous_positive() {
        let cm = CoalescenceModel::new(0.01, 0.5, 0.1, 1e-20);
        let rate = cm.neck_growth_rate_viscous(0.01, 1.0, 0.001);
        assert!(rate > 0.0, "viscous neck growth should be positive: {rate}");
    }
    #[test]
    fn test_coalescence_disjoining_pressure_negative() {
        let cm = CoalescenceModel::new(0.01, 0.5, 0.1, 1e-15);
        let pi = cm.disjoining_pressure(0.01);
        assert!(pi < 0.0, "vdW disjoining pressure should be negative: {pi}");
    }
    #[test]
    fn test_coalescence_neck_is_critical() {
        let cm = CoalescenceModel::new(0.01, 0.5, 0.1, 0.0);
        assert!(cm.neck_is_critical(0.4));
        assert!(!cm.neck_is_critical(0.6));
    }
    #[test]
    fn test_coalescence_drainage_timescale_positive() {
        let cm = CoalescenceModel::new(0.01, 0.5, 0.1, 0.0);
        let tau = cm.drainage_timescale(0.01);
        assert!(tau > 0.0 && tau.is_finite());
    }
    #[test]
    fn test_breakup_rp_growth_stable_high_k() {
        let bm = BreakupModel::new(12.0, 0.01, 1.0);
        let r0 = 1.0_f64;
        let k = 2.0 / r0;
        let sigma = bm.rayleigh_plateau_growth_rate(r0, k);
        assert_eq!(sigma, 0.0, "short-wave modes should be stable");
    }
    #[test]
    fn test_breakup_rp_growth_unstable() {
        let bm = BreakupModel::new(12.0, 0.01, 1.0);
        let r0 = 5.0_f64;
        let k_star = bm.most_unstable_wavenumber(r0);
        let sigma = bm.rayleigh_plateau_growth_rate(r0, k_star);
        assert!(
            sigma > 0.0,
            "most-unstable mode should have positive growth rate: {sigma}"
        );
    }
    #[test]
    fn test_breakup_breakup_time_finite() {
        let bm = BreakupModel::new(12.0, 0.01, 1.0);
        let t = bm.breakup_time(5.0);
        assert!(
            t > 0.0 && t.is_finite(),
            "breakup time should be finite: {t}"
        );
    }
    #[test]
    fn test_breakup_child_radius_larger_than_jet_radius() {
        let bm = BreakupModel::new(12.0, 0.01, 1.0);
        let r_child = bm.child_radius_plateau(5.0);
        assert!(r_child > 0.0);
    }
    #[test]
    fn test_breakup_weber_criterion_true() {
        let bm = BreakupModel::new(1.0, 0.01, 1.0);
        assert!(bm.weber_breakup(1.0, 10.0));
    }
    #[test]
    fn test_breakup_weber_criterion_false() {
        let bm = BreakupModel::new(1000.0, 0.01, 1.0);
        assert!(!bm.weber_breakup(1.0, 0.01));
    }
    #[test]
    fn test_contact_angle_young_equation() {
        let ca = ContactAngle::from_energies(0.072, 0.040, 0.030);
        let cos_expected = (0.040_f64 - 0.030_f64) / 0.072_f64;
        assert!((ca.theta_young.cos() - cos_expected).abs() < 1e-10);
    }
    #[test]
    fn test_contact_angle_spreading_coefficient() {
        let ca = ContactAngle::from_energies(0.072, 0.040, 0.030);
        let s = ca.spreading_coefficient();
        assert!((s - (0.040 - 0.030 - 0.072)).abs() < 1e-12);
    }
    #[test]
    fn test_contact_angle_work_of_adhesion_positive() {
        let ca = ContactAngle::from_energies(0.072, 0.070, 0.030);
        assert!(ca.work_of_adhesion() > 0.0);
    }
    #[test]
    fn test_contact_angle_hydrophilic() {
        let ca = ContactAngle::from_energies(0.072, 0.070, 0.030);
        assert!(ca.is_hydrophilic(), "should be hydrophilic");
    }
    #[test]
    fn test_contact_angle_cox_voinov_advancing() {
        let ca = ContactAngle::from_energies(0.072, 0.040, 0.030);
        let theta_dyn = ca.cox_voinov_angle(0.1, 10.0);
        assert!(theta_dyn >= 0.0);
    }
    #[test]
    fn test_contact_angle_wenzel_more_hydrophilic() {
        let mut ca = ContactAngle::from_energies(0.072, 0.070, 0.030);
        ca.roughness = 1.5;
        let theta_w = ca.wenzel_angle();
        assert!(
            theta_w <= ca.theta_young + 1e-10,
            "Wenzel should enhance wettability"
        );
    }
    #[test]
    fn test_contact_angle_cassie_baxter_computed() {
        let mut ca = ContactAngle::from_energies(0.072, 0.040, 0.030);
        ca.solid_fraction = 0.5;
        let theta_cb = ca.cassie_baxter_angle();
        assert!(theta_cb >= 0.0 && theta_cb <= PI);
    }
    #[test]
    fn test_contact_angle_capillary_length_positive() {
        let ca = ContactAngle::from_energies(0.072, 0.040, 0.030);
        let lc = ca.capillary_length(1000.0, 9.81);
        assert!(lc > 0.0 && lc < 0.01, "capillary length ~2.7mm: {lc}");
    }
    #[test]
    fn test_contact_angle_superhydrophobic() {
        let cos_160 = 160.0_f64.to_radians().cos();
        let sigma_lg = 1.0_f64;
        let sigma_sl = 0.0_f64;
        let sigma_sg = sigma_sl + sigma_lg * cos_160;
        let ca = ContactAngle::from_energies(sigma_lg, sigma_sg, sigma_sl);
        assert!(
            ca.is_superhydrophobic(),
            "θ≈160° should be superhydrophobic"
        );
    }
}
