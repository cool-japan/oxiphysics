//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{LayerModel, ProcessWindow, ScanStrategy};

/// Return a [`LayerModel`] preset for typical Ti-6Al-4V LPBF parameters.
pub fn ti64_layer_preset() -> LayerModel {
    LayerModel::new(
        30.0e-6_f64,
        120.0e-6_f64,
        55.0e-6_f64,
        1.2_f64,
        200.0_f64,
        ScanStrategy::Rotating67,
    )
}
/// Return a [`ProcessWindow`] preset for Ti-6Al-4V.
pub fn ti64_process_window() -> ProcessWindow {
    ProcessWindow::new(
        4430.0_f64,
        40.0e9_f64,
        120.0e9_f64,
        360.0_f64,
        500.0_f64,
        2.0e-6_f64,
    )
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::additive_manufacturing::AnisotropicAM;
    use crate::additive_manufacturing::BinderJetting;
    use crate::additive_manufacturing::BuildDirection;
    use crate::additive_manufacturing::DefectModels;
    use crate::additive_manufacturing::DirectedEnergyDeposition;
    use crate::additive_manufacturing::DistortionCompensation;
    use crate::additive_manufacturing::GoldakHeatSource;
    use crate::additive_manufacturing::MeltPoolGeometry;
    use crate::additive_manufacturing::MultiMaterialAM;
    use crate::additive_manufacturing::ProcessWindowOptimizer;
    use crate::additive_manufacturing::ResidualStress;
    use crate::additive_manufacturing::RosenthalSolution;
    use crate::additive_manufacturing::SupportStructure;
    use crate::additive_manufacturing::ThermalHistory;
    use std::f64::consts::PI;
    fn default_layer() -> LayerModel {
        LayerModel::new(
            30.0e-6_f64,
            120.0e-6_f64,
            55.0e-6_f64,
            1.0_f64,
            200.0_f64,
            ScanStrategy::Rotating67,
        )
    }
    #[test]
    fn test_energy_density_positive() {
        let layer = default_layer();
        let ev = layer.volumetric_energy_density();
        assert!(
            ev.is_finite() && ev > 0.0_f64,
            "Energy density must be positive: {ev}"
        );
    }
    #[test]
    fn test_linear_energy_density() {
        let layer = default_layer();
        let led = layer.linear_energy_density();
        assert!(
            (led - 200.0_f64).abs() < 1.0e-10_f64,
            "LED should be 200 J/m, got {led}"
        );
    }
    #[test]
    fn test_hatch_overlap_range() {
        let layer = default_layer();
        let ov = layer.hatch_overlap();
        assert!(ov > -1.0_f64 && ov < 1.0_f64, "Overlap out of range: {ov}");
    }
    #[test]
    fn test_num_tracks_min_one() {
        let layer = default_layer();
        assert_eq!(layer.num_tracks(1.0e-8_f64), 1);
        assert!(layer.num_tracks(1.0e-2_f64) > 1);
    }
    #[test]
    fn test_melt_pool_aspect_ratio() {
        let mp = MeltPoolGeometry {
            half_length: 200.0e-6_f64,
            half_width: 80.0e-6_f64,
            depth: 60.0e-6_f64,
        };
        let ar = mp.aspect_ratio();
        assert!((ar - 60.0e-6_f64 / (2.0_f64 * 80.0e-6_f64)).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_melt_pool_volume() {
        let mp = MeltPoolGeometry {
            half_length: 200.0e-6_f64,
            half_width: 80.0e-6_f64,
            depth: 60.0e-6_f64,
        };
        assert!(mp.volume() > 0.0_f64);
    }
    #[test]
    fn test_thermal_history_from_process() {
        let layer = default_layer();
        let th = ThermalHistory::from_process(
            &layer, 0.35_f64, 6.7_f64, 4430.0_f64, 560.0_f64, 1933.0_f64, 300.0_f64,
        );
        assert!(
            th.peak_temperature > 300.0_f64,
            "Peak temperature must exceed ambient"
        );
        assert!(th.cooling_rate > 0.0_f64, "Cooling rate must be positive");
    }
    #[test]
    fn test_solidification_rate_nonneg() {
        let layer = default_layer();
        let th = ThermalHistory::from_process(
            &layer, 0.35_f64, 6.7_f64, 4430.0_f64, 560.0_f64, 1933.0_f64, 300.0_f64,
        );
        assert!(th.solidification_rate() >= 0.0_f64);
    }
    #[test]
    fn test_pdas_positive() {
        let layer = default_layer();
        let th = ThermalHistory::from_process(
            &layer, 0.35_f64, 6.7_f64, 4430.0_f64, 560.0_f64, 1933.0_f64, 300.0_f64,
        );
        assert!(th.primary_dendrite_arm_spacing() > 0.0_f64);
    }
    #[test]
    fn test_goldak_zero_far() {
        let g = GoldakHeatSource::new(
            200.0_f64,
            0.35_f64,
            200.0e-6_f64,
            600.0e-6_f64,
            100.0e-6_f64,
            80.0e-6_f64,
        );
        let q = g.heat_flux(0.1_f64, 0.1_f64, 0.1_f64);
        assert!(
            q < 1.0e-100_f64,
            "Heat flux should be negligible far away: {q}"
        );
    }
    #[test]
    fn test_goldak_peak_positive() {
        let g = GoldakHeatSource::new(
            200.0_f64,
            0.35_f64,
            200.0e-6_f64,
            600.0e-6_f64,
            100.0e-6_f64,
            80.0e-6_f64,
        );
        let q = g.peak_heat_flux();
        assert!(q > 0.0_f64, "Peak heat flux must be positive: {q}");
    }
    #[test]
    fn test_goldak_front_rear_boundary() {
        let g = GoldakHeatSource::new(
            200.0_f64,
            0.35_f64,
            200.0e-6_f64,
            600.0e-6_f64,
            100.0e-6_f64,
            80.0e-6_f64,
        );
        let q_front = g.heat_flux(1.0e-9_f64, 0.0_f64, 0.0_f64);
        let q_rear = g.heat_flux(-1.0e-9_f64, 0.0_f64, 0.0_f64);
        assert!(q_front > 0.0_f64 && q_rear > 0.0_f64);
    }
    #[test]
    fn test_rosenthal_decreases_with_distance() {
        let rs =
            RosenthalSolution::new(200.0_f64, 0.35_f64, 1.2_f64, 6.7_f64, 2.9e-6_f64, 300.0_f64);
        let t1 = rs.temperature(-100.0e-6_f64, 0.0_f64, 0.0_f64);
        let t2 = rs.temperature(-500.0e-6_f64, 0.0_f64, 0.0_f64);
        assert!(
            t1 > t2,
            "Temperature should decrease further from source: t1={t1}, t2={t2}"
        );
    }
    #[test]
    fn test_rosenthal_ambient_at_infinity() {
        let rs =
            RosenthalSolution::new(200.0_f64, 0.35_f64, 1.2_f64, 6.7_f64, 2.9e-6_f64, 300.0_f64);
        let t_near = rs.temperature(-0.01_f64, 0.0_f64, 0.0_f64);
        let t_far = rs.temperature(-0.1_f64, 0.0_f64, 0.0_f64);
        assert!(
            t_far >= 300.0_f64,
            "Far-field temperature should be at or above ambient: {t_far}"
        );
        assert!(
            t_far < t_near,
            "Temperature should decrease further from source: t_near={t_near}, t_far={t_far}"
        );
    }
    #[test]
    fn test_rosenthal_melt_pool_width_positive() {
        let rs =
            RosenthalSolution::new(200.0_f64, 0.35_f64, 1.2_f64, 6.7_f64, 2.9e-6_f64, 300.0_f64);
        let w = rs.melt_pool_width(1933.0_f64);
        assert!(w > 0.0_f64, "Melt pool width must be positive: {w}");
    }
    #[test]
    fn test_residual_stress_capped_at_yield() {
        let rs = ResidualStress::new(110.0e9_f64, 8.6e-6_f64, 900.0e6_f64, 3000.0_f64, 300.0_f64);
        let sigma = rs.longitudinal_stress();
        assert!(
            sigma <= rs.yield_strength,
            "Stress must not exceed yield: {sigma}"
        );
    }
    #[test]
    fn test_residual_stress_ratio_range() {
        let rs = ResidualStress::new(110.0e9_f64, 8.6e-6_f64, 900.0e6_f64, 3000.0_f64, 300.0_f64);
        let ratio = rs.stress_ratio();
        assert!(
            (0.0_f64..=1.0_f64).contains(&ratio),
            "Stress ratio out of range: {ratio}"
        );
    }
    #[test]
    fn test_accumulated_stress_capped() {
        let rs = ResidualStress::new(110.0e9_f64, 8.6e-6_f64, 900.0e6_f64, 1500.0_f64, 300.0_f64);
        let acc = rs.accumulated_stress(1000, 0.1_f64);
        assert!(
            acc <= rs.yield_strength,
            "Accumulated stress must not exceed yield: {acc}"
        );
    }
    #[test]
    fn test_delamination_risk_high_stress() {
        let rs = ResidualStress::new(110.0e9_f64, 8.6e-6_f64, 200.0e6_f64, 5000.0_f64, 300.0_f64);
        assert!(
            rs.delamination_risk(),
            "Delamination should be flagged at high stress"
        );
    }
    #[test]
    fn test_anisotropy_index() {
        let am = AnisotropicAM::new(
            BuildDirection::PlusZ,
            110.0e9_f64,
            120.0e9_f64,
            900.0e6_f64,
            950.0e6_f64,
            0.7_f64,
            8.0_f64,
        );
        let ai = am.anisotropy_index();
        assert!((ai - (120.0e9_f64 - 110.0e9_f64) / 110.0e9_f64).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_effective_modulus_theta_zero() {
        let am = AnisotropicAM::new(
            BuildDirection::PlusZ,
            110.0e9_f64,
            120.0e9_f64,
            900.0e6_f64,
            950.0e6_f64,
            0.7_f64,
            8.0_f64,
        );
        let e = am.effective_modulus(0.0_f64);
        assert!(
            (e - 110.0e9_f64).abs() < 1.0_f64,
            "At theta=0 effective modulus should be e_parallel: {e}"
        );
    }
    #[test]
    fn test_effective_modulus_theta_90() {
        let am = AnisotropicAM::new(
            BuildDirection::PlusZ,
            110.0e9_f64,
            120.0e9_f64,
            900.0e6_f64,
            950.0e6_f64,
            0.7_f64,
            8.0_f64,
        );
        let e = am.effective_modulus(PI / 2.0_f64);
        assert!(
            (e - 120.0e9_f64).abs() < 1.0_f64,
            "At theta=90 effective modulus should be e_transverse: {e}"
        );
    }
    #[test]
    fn test_lof_below_threshold() {
        let dm = DefectModels::new(
            0.96_f64,
            30.0e9_f64,
            40.0e9_f64,
            120.0e9_f64,
            1.0e9_f64,
            50.0e3_f64,
            50.0e-6_f64,
        );
        assert!(dm.has_lack_of_fusion());
    }
    #[test]
    fn test_keyhole_above_threshold() {
        let dm = DefectModels::new(
            0.98_f64,
            150.0e9_f64,
            40.0e9_f64,
            120.0e9_f64,
            1.0e9_f64,
            50.0e3_f64,
            50.0e-6_f64,
        );
        assert!(dm.has_keyhole_porosity());
    }
    #[test]
    fn test_in_process_window() {
        let dm = DefectModels::new(
            0.998_f64,
            80.0e9_f64,
            40.0e9_f64,
            120.0e9_f64,
            1.0e9_f64,
            50.0e3_f64,
            10.0e-6_f64,
        );
        assert!(dm.in_process_window());
    }
    #[test]
    fn test_process_window_low_density() {
        let pw = ti64_process_window();
        let rho = pw.relative_density(1.0e6_f64);
        assert!(
            rho < 0.1_f64,
            "Relative density should be near 0 at very low energy: {rho}"
        );
    }
    #[test]
    fn test_process_window_high_density() {
        let pw = ti64_process_window();
        let rho = pw.relative_density(1.0e13_f64);
        assert!(
            rho > 0.9_f64,
            "Relative density should be near 1 at very high energy: {rho}"
        );
    }
    #[test]
    fn test_process_window_in_window() {
        let pw = ti64_process_window();
        assert!(
            pw.is_in_window(70.0e9_f64),
            "70 GJ/m³ should be in the Ti64 window"
        );
        assert!(
            !pw.is_in_window(1.0e9_f64),
            "1 GJ/m³ should be below window"
        );
    }
    #[test]
    fn test_ti64_preset_energy_density() {
        let layer = ti64_layer_preset();
        let ev = layer.volumetric_energy_density();
        assert!(
            ev > 1.0e9_f64 && ev < 1.0e12_f64,
            "Ti64 energy density should be realistic: {ev}"
        );
    }
    #[test]
    fn test_hardness_prediction_positive() {
        let pw = ti64_process_window();
        let hv = pw.hardness_prediction(70.0e9_f64, 2.0e-6_f64);
        assert!(hv > 0.0_f64, "Hardness must be positive: {hv}");
    }
    #[test]
    fn test_support_required_downward() {
        let ss = SupportStructure::new(45.0_f64, 4430.0_f64, 0.1_f64);
        assert!(ss.requires_support([0.0_f64, 0.0_f64, -1.0_f64]));
    }
    #[test]
    fn test_no_support_upward() {
        let ss = SupportStructure::new(45.0_f64, 4430.0_f64, 0.1_f64);
        assert!(!ss.requires_support([0.0_f64, 0.0_f64, 1.0_f64]));
    }
    #[test]
    fn test_support_volume_nonneg() {
        let mut ss = SupportStructure::new(45.0_f64, 4430.0_f64, 0.1_f64);
        ss.estimate_support(1.0e-5_f64, 0.3_f64, 0.05_f64);
        assert!(ss.support_volume_m3 >= 0.0_f64);
    }
    #[test]
    fn test_support_mass_scales_with_fill() {
        let mut ss_sparse = SupportStructure::new(45.0_f64, 4430.0_f64, 0.05_f64);
        let mut ss_dense = SupportStructure::new(45.0_f64, 4430.0_f64, 0.20_f64);
        ss_sparse.estimate_support(1.0e-5_f64, 0.3_f64, 0.05_f64);
        ss_dense.estimate_support(1.0e-5_f64, 0.3_f64, 0.05_f64);
        assert!(ss_dense.support_mass_kg() > ss_sparse.support_mass_kg());
    }
    #[test]
    fn test_waste_fraction_range() {
        let mut ss = SupportStructure::new(45.0_f64, 4430.0_f64, 0.1_f64);
        ss.estimate_support(1.0e-5_f64, 0.3_f64, 0.05_f64);
        let wf = ss.waste_fraction(0.05_f64);
        assert!(
            (0.0_f64..=1.0_f64).contains(&wf),
            "Waste fraction out of range: {wf}"
        );
    }
    #[test]
    fn test_distortion_compensation_invert() {
        let dc = DistortionCompensation::new(-1.0_f64, 1.0_f64, 1.0e-9_f64);
        let disp = vec![0.001_f64, -0.002_f64, 0.003_f64];
        let comp = dc.apply(&disp);
        assert!((comp[0] - (-0.001_f64)).abs() < 1.0e-12_f64);
        assert!((comp[1] - 0.002_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_distortion_compensation_clamped() {
        let dc = DistortionCompensation::new(-2.0_f64, 0.5e-3_f64, 1.0e-9_f64);
        let large = vec![1.0_f64];
        let comp = dc.apply(&large);
        assert!(comp[0].abs() <= 0.5e-3_f64 + 1.0e-12_f64);
    }
    #[test]
    fn test_rms_displacement_empty() {
        assert_eq!(DistortionCompensation::rms_displacement(&[]), 0.0_f64);
    }
    #[test]
    fn test_distortion_iterate_length() {
        let mut dc = DistortionCompensation::new(-1.0_f64, 0.01_f64, 1.0e-9_f64);
        let input = vec![0.001_f64; 10];
        let output = dc.iterate(&input, 5);
        assert_eq!(output.len(), input.len());
    }
    #[test]
    fn test_fgm_vf_at_zero() {
        let fgm = MultiMaterialAM::new(
            110.0e9_f64,
            200.0e9_f64,
            4430.0_f64,
            7900.0_f64,
            6.7_f64,
            15.0_f64,
            1.0_f64,
            0.01_f64,
        );
        assert_eq!(fgm.volume_fraction_b(0.0_f64), 0.0_f64);
    }
    #[test]
    fn test_fgm_vf_at_top() {
        let fgm = MultiMaterialAM::new(
            110.0e9_f64,
            200.0e9_f64,
            4430.0_f64,
            7900.0_f64,
            6.7_f64,
            15.0_f64,
            1.0_f64,
            0.01_f64,
        );
        assert!((fgm.volume_fraction_b(0.01_f64) - 1.0_f64).abs() < 1.0e-10_f64);
    }
    #[test]
    fn test_fgm_effective_modulus_bounded() {
        let fgm = MultiMaterialAM::new(
            110.0e9_f64,
            200.0e9_f64,
            4430.0_f64,
            7900.0_f64,
            6.7_f64,
            15.0_f64,
            2.0_f64,
            0.05_f64,
        );
        for i in 0..=10 {
            let z = 0.05_f64 * i as f64 / 10.0_f64;
            let e = fgm.effective_modulus(z);
            assert!((110.0e9_f64 - 1.0_f64..=200.0e9_f64 + 1.0_f64).contains(&e));
        }
    }
    #[test]
    fn test_fgm_average_modulus_analytical() {
        let fgm = MultiMaterialAM::new(
            100.0e9_f64,
            200.0e9_f64,
            4000.0_f64,
            8000.0_f64,
            10.0_f64,
            20.0_f64,
            1.0_f64,
            1.0_f64,
        );
        let expected = (100.0e9_f64 + 200.0e9_f64) / 2.0_f64;
        let avg = fgm.average_modulus();
        assert!(
            (avg - expected).abs() / expected < 1.0e-10_f64,
            "avg={avg}, expected={expected}"
        );
    }
    #[test]
    fn test_binder_jetting_void_fraction() {
        let bj = BinderJetting::new(0.1e-3_f64, 0.6_f64, 0.8_f64, 50.0e-6_f64, 0.18_f64);
        assert!((bj.void_fraction() + bj.packing_density - 1.0_f64).abs() < 1.0e-12_f64);
    }
    #[test]
    fn test_binder_jetting_sintered_smaller() {
        let bj = BinderJetting::new(0.1e-3_f64, 0.6_f64, 0.8_f64, 50.0e-6_f64, 0.18_f64);
        let green = 0.1_f64;
        let sintered = bj.sintered_dimension(green);
        assert!(sintered < green, "sintered < green: {sintered} < {green}");
    }
    #[test]
    fn test_binder_jetting_volumetric_shrinkage() {
        let bj = BinderJetting::new(0.1e-3_f64, 0.6_f64, 0.8_f64, 50.0e-6_f64, 0.18_f64);
        let vs = bj.volumetric_shrinkage();
        assert!(vs > 0.0_f64 && vs < 1.0_f64, "vs={vs}");
    }
    #[test]
    fn test_binder_jetting_spreading_speed() {
        let bj = BinderJetting::new(0.1e-3_f64, 0.6_f64, 0.8_f64, 50.0e-6_f64, 0.18_f64);
        assert!(bj.max_spreading_speed_m_s() > 0.0_f64);
    }
    #[test]
    fn test_ded_effective_power() {
        let ded = DirectedEnergyDeposition::new(
            3000.0_f64, 0.002_f64, 0.01_f64, 2.0e-3_f64, 0.7_f64, 15.0_f64, 0.35_f64,
        );
        assert!(ded.effective_power() <= ded.laser_power);
    }
    #[test]
    fn test_ded_clad_height_positive() {
        let ded = DirectedEnergyDeposition::new(
            3000.0_f64, 0.002_f64, 0.01_f64, 2.0e-3_f64, 0.7_f64, 15.0_f64, 0.35_f64,
        );
        let h = ded.clad_height_m(7900.0_f64);
        assert!(h > 0.0_f64, "Clad height must be positive: {h}");
    }
    #[test]
    fn test_ded_dilution_ratio_range() {
        let ded = DirectedEnergyDeposition::new(
            3000.0_f64, 0.002_f64, 0.01_f64, 2.0e-3_f64, 0.7_f64, 15.0_f64, 0.35_f64,
        );
        let d = ded.dilution_ratio();
        assert!(
            (0.0_f64..=1.0_f64).contains(&d),
            "Dilution ratio out of range: {d}"
        );
    }
    #[test]
    fn test_ded_higher_power_more_dilution() {
        let ded_low = DirectedEnergyDeposition::new(
            500.0_f64, 0.002_f64, 0.01_f64, 2.0e-3_f64, 0.7_f64, 15.0_f64, 0.35_f64,
        );
        let ded_high = DirectedEnergyDeposition::new(
            5000.0_f64, 0.002_f64, 0.01_f64, 2.0e-3_f64, 0.7_f64, 15.0_f64, 0.35_f64,
        );
        assert!(ded_high.dilution_ratio() > ded_low.dilution_ratio());
    }
    #[test]
    fn test_pwo_optimal_within_range() {
        let opt = ProcessWindowOptimizer::new(
            100.0_f64,
            400.0_f64,
            0.5_f64,
            3.0_f64,
            30.0e-6_f64,
            120.0e-6_f64,
            40.0e9_f64,
            120.0e9_f64,
        );
        let (p, v) = opt.optimise(10);
        assert!((100.0_f64..=400.0_f64).contains(&p), "power={p}");
        assert!((0.5_f64..=3.0_f64).contains(&v), "speed={v}");
    }
    #[test]
    fn test_pwo_porosity_score_zero_inside() {
        let opt = ProcessWindowOptimizer::new(
            100.0_f64,
            400.0_f64,
            0.5_f64,
            3.0_f64,
            30.0e-6_f64,
            120.0e-6_f64,
            40.0e9_f64,
            120.0e9_f64,
        );
        let ev_mid = 0.5_f64 * (40.0e9_f64 + 120.0e9_f64);
        assert_eq!(opt.porosity_score(ev_mid), 0.0_f64);
    }
    #[test]
    fn test_pwo_porosity_score_below_low() {
        let opt = ProcessWindowOptimizer::new(
            100.0_f64,
            400.0_f64,
            0.5_f64,
            3.0_f64,
            30.0e-6_f64,
            120.0e-6_f64,
            40.0e9_f64,
            120.0e9_f64,
        );
        assert!(opt.porosity_score(1.0e9_f64) > 0.0_f64);
    }
    #[test]
    fn test_pwo_cost_lower_in_window() {
        let opt = ProcessWindowOptimizer::new(
            50.0_f64,
            500.0_f64,
            0.2_f64,
            4.0_f64,
            30.0e-6_f64,
            120.0e-6_f64,
            40.0e9_f64,
            120.0e9_f64,
        );
        let cost_out = opt.cost(50.0_f64, 4.0_f64);
        let cost_in = opt.cost(200.0_f64, 1.0_f64);
        assert!(
            cost_in < cost_out,
            "cost inside window ({cost_in}) should be lower than outside ({cost_out})"
        );
    }
}
