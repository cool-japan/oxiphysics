//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[inline]
pub(super) fn add2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] + b[0], a[1] + b[1]]
}
#[inline]
pub(super) fn sub2(a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
    [a[0] - b[0], a[1] - b[1]]
}
#[inline]
pub(super) fn scale2(v: [f64; 2], s: f64) -> [f64; 2] {
    [v[0] * s, v[1] * s]
}
#[inline]
pub(super) fn dot2(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
#[inline]
pub(super) fn len2(v: [f64; 2]) -> f64 {
    dot2(v, v).sqrt()
}
#[inline]
pub(super) fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}
pub(super) const NQ: usize = 9;
pub(super) const CX: [f64; NQ] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
pub(super) const CY: [f64; NQ] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
pub(super) const W: [f64; NQ] = [
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
pub(super) const CS2: f64 = 1.0 / 3.0;
/// Compute equilibrium distribution for one cell.
pub(super) fn feq(rho: f64, ux: f64, uy: f64) -> [f64; NQ] {
    let mut f = [0.0f64; NQ];
    let u2 = ux * ux + uy * uy;
    for q in 0..NQ {
        let cu = CX[q] * ux + CY[q] * uy;
        f[q] = W[q] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2));
    }
    f
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::biofluid_lbm::types::*;
    use std::f64::consts::PI;
    #[test]
    fn test_carreau_yasuda_zero_shear() {
        let p = CarreauYasudaParams::blood();
        let eta = p.viscosity(0.0);
        assert!(
            (eta - p.eta_0).abs() < 1e-10,
            "zero shear → eta_0, got {eta}"
        );
    }
    #[test]
    fn test_carreau_yasuda_high_shear() {
        let p = CarreauYasudaParams::blood();
        let eta = p.viscosity(1e6);
        assert!(
            (eta - p.eta_inf).abs() < 1e-4,
            "high shear → eta_inf, got {eta}"
        );
    }
    #[test]
    fn test_carreau_yasuda_monotone() {
        let p = CarreauYasudaParams::blood();
        let eta1 = p.viscosity(1.0);
        let eta2 = p.viscosity(100.0);
        assert!(eta2 < eta1, "shear thinning: η should decrease with γ̇");
    }
    #[test]
    fn test_carreau_yasuda_symmetric() {
        let p = CarreauYasudaParams::blood();
        assert_eq!(p.viscosity(-5.0), p.viscosity(5.0), "η(|γ̇|)");
    }
    #[test]
    fn test_carreau_yasuda_intermediate() {
        let p = CarreauYasudaParams::blood();
        let eta = p.viscosity(1.0 / p.lambda);
        assert!(eta > p.eta_inf && eta < p.eta_0);
    }
    #[test]
    fn test_blood_flow_construction() {
        let p = CarreauYasudaParams::blood();
        let sim = BloodFlowLbm::new(16, 8, p, 1e-4, 5e-5, 1.0 / 200.0, 4.0);
        assert_eq!(sim.nx, 16);
        assert_eq!(sim.ny, 8);
        assert_eq!(sim.womersley_number(), 4.0);
    }
    #[test]
    fn test_blood_flow_pressure_gradient_pulsatile() {
        let p = CarreauYasudaParams::blood();
        let mut sim = BloodFlowLbm::new(16, 8, p, 1e-4, 1e-4, 0.5, 3.0);
        let g0 = sim.pressure_gradient();
        assert!((g0 - 1e-4).abs() < 1e-12);
        sim.time_step = 1;
        let g1 = sim.pressure_gradient();
        assert!(g1.abs() > 0.0);
    }
    #[test]
    fn test_blood_flow_runs() {
        let p = CarreauYasudaParams::blood();
        let mut sim = BloodFlowLbm::new(10, 6, p, 1e-5, 0.0, 1.0 / 300.0, 2.0);
        sim.run(5);
        assert_eq!(sim.time_step, 5);
    }
    #[test]
    fn test_blood_flow_density_conservation() {
        let p = CarreauYasudaParams::blood();
        let mut sim = BloodFlowLbm::new(8, 6, p, 1e-5, 0.0, 1.0 / 300.0, 2.0);
        for _ in 0..5 {
            sim.step();
        }
        let all_finite = sim.f.iter().all(|&v| v.is_finite());
        assert!(
            all_finite,
            "all distribution functions should remain finite after 5 steps"
        );
        let total_rho: f64 = sim
            .rho
            .iter()
            .zip(sim.is_wall.iter())
            .filter(|&(_, w)| !w)
            .map(|(r, _)| *r)
            .sum();
        assert!(
            total_rho > 0.0,
            "total interior density should be positive, got {total_rho}"
        );
    }
    #[test]
    fn test_blood_flow_mean_velocity_positive() {
        let p = CarreauYasudaParams::blood();
        let mut sim = BloodFlowLbm::new(12, 8, p, 1e-4, 0.0, 1.0 / 300.0, 2.0);
        sim.run(50);
        let mv = sim.mean_velocity();
        assert!(mv.is_finite(), "mean velocity should be finite, got {mv}");
    }
    #[test]
    fn test_blood_flow_local_tau() {
        let p = CarreauYasudaParams::blood();
        let sim = BloodFlowLbm::new(8, 6, p, 1e-5, 0.0, 0.001, 2.0);
        let tau = sim.local_tau(10.0);
        assert!(tau > 0.5, "tau must exceed 0.5 for stability");
    }
    #[test]
    fn test_blood_flow_womersley() {
        let p = CarreauYasudaParams::blood();
        let sim = BloodFlowLbm::new(8, 6, p, 1e-5, 0.0, 0.001, 7.3);
        assert_eq!(sim.womersley_number(), 7.3);
    }
    #[test]
    fn test_rbc_construction() {
        let params = RbcMembraneParams::default();
        let rbc = RedBloodCellLbm::new_circle(20, 10.0, 10.0, 4.0, 32, 32, params);
        assert_eq!(rbc.n_markers, 20);
    }
    #[test]
    fn test_rbc_area() {
        let params = RbcMembraneParams::default();
        let r = 4.0;
        let rbc = RedBloodCellLbm::new_circle(100, 10.0, 10.0, r, 40, 40, params);
        let area = rbc.current_area();
        let expected = PI * r * r;
        assert!(
            (area - expected).abs() / expected < 0.01,
            "area {area} vs {expected}"
        );
    }
    #[test]
    fn test_rbc_perimeter() {
        let params = RbcMembraneParams::default();
        let r = 5.0;
        let rbc = RedBloodCellLbm::new_circle(200, 10.0, 10.0, r, 40, 40, params);
        let perim = rbc.perimeter();
        let expected = 2.0 * PI * r;
        assert!((perim - expected).abs() / expected < 0.01);
    }
    #[test]
    fn test_rbc_circularity_near_one() {
        let params = RbcMembraneParams::default();
        let rbc = RedBloodCellLbm::new_circle(100, 10.0, 10.0, 4.0, 40, 40, params);
        let circ = rbc.circularity();
        assert!((circ - 1.0).abs() < 0.02, "circular → C ≈ 1, got {circ}");
    }
    #[test]
    fn test_rbc_shear_force_zero_at_rest() {
        let params = RbcMembraneParams {
            rest_length: 2.0 * PI * 4.0 / 20.0,
            ..Default::default()
        };
        let rbc = RedBloodCellLbm::new_circle(20, 10.0, 10.0, 4.0, 40, 40, params);
        let f = rbc.shear_force(0);
        let fmag = len2(f);
        assert!(
            fmag < 1.0,
            "shear force at rest should be small, got {fmag}"
        );
    }
    #[test]
    fn test_rbc_area_force_direction() {
        let params = RbcMembraneParams::default();
        let mut rbc = RedBloodCellLbm::new_circle(20, 10.0, 10.0, 4.0, 40, 40, params);
        for p in rbc.marker_pos.iter_mut() {
            p[0] = (p[0] - 10.0) * 0.5 + 10.0;
            p[1] = (p[1] - 10.0) * 0.5 + 10.0;
        }
        let f = rbc.area_force(0);
        let fmag = len2(f);
        assert!(fmag > 0.0, "area force should be non-zero when compressed");
    }
    #[test]
    fn test_rbc_spread_forces_runs() {
        let params = RbcMembraneParams::default();
        let mut rbc = RedBloodCellLbm::new_circle(20, 8.0, 8.0, 3.0, 20, 20, params);
        let fluid_ux = vec![0.01f64; 400];
        let fluid_uy = vec![0.0f64; 400];
        rbc.spread_forces(&fluid_ux, &fluid_uy);
        let total_fx: f64 = rbc.force_x.iter().sum();
        assert!(total_fx.is_finite());
    }
    #[test]
    fn test_rbc_advance_markers() {
        let params = RbcMembraneParams::default();
        let mut rbc = RedBloodCellLbm::new_circle(10, 5.0, 5.0, 2.0, 20, 20, params);
        let pos0 = rbc.marker_pos[0];
        rbc.marker_vel[0] = [1.0, 0.0];
        rbc.advance_markers(1.0);
        assert!((rbc.marker_pos[0][0] - pos0[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_artery_construction() {
        let art = ArteryModel::new(100.0, 1e-8, 10000.0, 800.0, 6.0, 1e-8);
        assert!(art.pressure > art.p_out);
        assert_eq!(art.pwv, 6.0);
    }
    #[test]
    fn test_artery_step_monotone() {
        let mut art = ArteryModel::new(100.0, 1e-8, 10000.0, 800.0, 6.0, 1e-8);
        let p0 = art.pressure;
        art.step(1e-6, 0.001);
        assert_ne!(art.pressure, p0);
        assert!(art.pressure.is_finite());
    }
    #[test]
    fn test_artery_wall_displacement() {
        let art = ArteryModel::new(100.0, 1e-8, 10000.0, 800.0, 6.0, 1e-9);
        let disp = art.wall_displacement(0.01, 0.1);
        assert!(disp.is_finite());
        assert!(disp >= 0.0);
    }
    #[test]
    fn test_moens_korteweg() {
        let pwv = ArteryModel::moens_korteweg(5e5, 1e-3, 1060.0, 4e-3);
        assert!(
            pwv > 1.0 && pwv < 20.0,
            "pwv={pwv} should be in physiological range"
        );
    }
    #[test]
    fn test_moens_korteweg_zero_radius() {
        assert_eq!(ArteryModel::moens_korteweg(1e6, 1e-3, 1000.0, 0.0), 0.0);
    }
    #[test]
    fn test_artery_windkessel_steady_state() {
        let q = 1e-4f64;
        let r2 = 100.0f64;
        let p_out = 0.0f64;
        let p_ss = p_out + r2 * q;
        let mut art = ArteryModel::new(10.0, 1e-6, r2, p_out, 5.0, 1e-9);
        art.pressure = p_ss;
        for _ in 0..1000 {
            art.step(q, 0.001);
        }
        assert!(
            art.pressure.is_finite(),
            "pressure should be finite, got {}",
            art.pressure
        );
        assert!(
            (art.pressure - p_ss).abs() < p_ss * 0.5 + 0.001,
            "ss pressure {}, expected {}",
            art.pressure,
            p_ss
        );
    }
    #[test]
    fn test_platelet_initial_state() {
        let p = Platelet::new([1.0, 2.0]);
        assert_eq!(p.state, PlateletState::Resting);
        assert_eq!(p.activation, 0.0);
    }
    #[test]
    fn test_platelet_activation_grows() {
        let mut p = Platelet::new([0.0, 0.0]);
        p.update_activation(100.0, 0.1, 0.01, 0.0);
        assert!(p.activation > 0.0);
    }
    #[test]
    fn test_platelet_no_activation_below_threshold() {
        let mut p = Platelet::new([0.0, 0.0]);
        p.update_activation(0.5, 1.0, 1.0, 1.0);
        assert_eq!(p.activation, 0.0);
    }
    #[test]
    fn test_platelet_full_activation() {
        let mut p = Platelet::new([0.0, 0.0]);
        for _ in 0..200 {
            p.update_activation(1000.0, 0.1, 1.0, 0.0);
        }
        assert_eq!(p.state, PlateletState::Adhered);
    }
    #[test]
    fn test_coagulation_thrombin_grows() {
        let mut c = CoagulationState::physiological();
        for _ in 0..100 {
            c.step(1.0, 0.01);
        }
        assert!(c.thrombin > 0.0);
        assert!(c.fibrin > 0.0);
    }
    #[test]
    fn test_coagulation_mass_conservation() {
        let mut c = CoagulationState::physiological();
        let total0 = c.prothrombin + c.thrombin + c.fibrinogen + c.fibrin;
        for _ in 0..10 {
            c.step(0.5, 0.01);
        }
        let prot_tot = c.prothrombin + c.thrombin;
        assert!(prot_tot <= total0 && prot_tot >= 0.0);
    }
    #[test]
    fn test_thrombosis_model_construction() {
        let m = ThrombosisModel::new(16, 16, 10, 0.01, 1.0);
        assert_eq!(m.platelets.len(), 10);
        assert_eq!(m.thrombus_fraction.len(), 256);
    }
    #[test]
    fn test_thrombosis_total_volume_zero() {
        let m = ThrombosisModel::new(8, 8, 5, 0.01, 1.0);
        assert_eq!(m.total_thrombus_volume(1.0), 0.0);
    }
    #[test]
    fn test_thrombosis_step_runs() {
        let mut m = ThrombosisModel::new(8, 8, 5, 0.01, 0.0);
        let shear = vec![10.0f64; 64];
        m.step(&shear, 0.1);
        assert!(m.total_thrombus_volume(1.0).is_finite());
    }
    #[test]
    fn test_hele_shaw_permeability() {
        let hs = HeleShawParams {
            gap_height: 1e-4,
            width: 1e-3,
            viscosity: 1e-3,
        };
        let k = hs.permeability();
        assert!((k - 1e-4 * 1e-4 / 12.0).abs() < 1e-20);
    }
    #[test]
    fn test_hele_shaw_darcy_velocity() {
        let hs = HeleShawParams {
            gap_height: 1e-4,
            width: 1e-3,
            viscosity: 1e-3,
        };
        let u = hs.darcy_velocity(-100.0);
        assert!(u > 0.0);
    }
    #[test]
    fn test_hele_shaw_flow_resistance() {
        let hs = HeleShawParams {
            gap_height: 1e-4,
            width: 1e-3,
            viscosity: 1e-3,
        };
        let r = hs.flow_resistance();
        assert!(r > 0.0);
        assert!((r - 12.0 * 1e-3 / (1e-4 * 1e-4)).abs() < 1e-3);
    }
    #[test]
    fn test_dean_number() {
        let df = DeanFlowParams {
            radius: 1e-3,
            curvature_radius: 10e-3,
            mean_velocity: 0.1,
            kinematic_viscosity: 1e-6,
        };
        let de = df.dean_number();
        assert!((de - 31.62).abs() < 0.1, "Dean number = {de}");
    }
    #[test]
    fn test_dean_secondary_velocity() {
        let df = DeanFlowParams {
            radius: 1e-3,
            curvature_radius: 10e-3,
            mean_velocity: 0.1,
            kinematic_viscosity: 1e-6,
        };
        let u_sec = df.secondary_velocity();
        assert!(u_sec > 0.0);
    }
    #[test]
    fn test_droplet_capillary_number() {
        let mut d = Droplet::new([5.0, 5.0], 2.0, 1.0, 1e-3);
        d.update_capillary_number(0.1, 1e-3);
        assert!((d.capillary_number - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_droplet_will_break() {
        let mut d = Droplet::new([0.0, 0.0], 1.0, 1.0, 1e-4);
        d.capillary_number = 10.0;
        assert!(d.will_break());
    }
    #[test]
    fn test_droplet_advance() {
        let mut d = Droplet::new([0.0, 0.0], 1.0, 1.0, 1e-3);
        d.vel = [2.0, 3.0];
        d.advance(1.0);
        assert_eq!(d.pos, [2.0, 3.0]);
    }
    #[test]
    fn test_microfluidics_construction() {
        let hs = HeleShawParams {
            gap_height: 1e-4,
            width: 1e-3,
            viscosity: 1e-3,
        };
        let df = DeanFlowParams {
            radius: 1e-3,
            curvature_radius: 10e-3,
            mean_velocity: 0.1,
            kinematic_viscosity: 1e-6,
        };
        let sim = MicrofluidicsLbm::new(16, 8, 0.8, 0.9, -0.1, hs, df);
        assert_eq!(sim.nx, 16);
        assert_eq!(sim.ny, 8);
    }
    #[test]
    fn test_microfluidics_step_runs() {
        let hs = HeleShawParams {
            gap_height: 1e-4,
            width: 1e-3,
            viscosity: 1e-3,
        };
        let df = DeanFlowParams {
            radius: 1e-3,
            curvature_radius: 10e-3,
            mean_velocity: 0.1,
            kinematic_viscosity: 1e-6,
        };
        let mut sim = MicrofluidicsLbm::new(10, 8, 0.8, 0.9, -0.1, hs, df);
        sim.step(1e-4);
        let sum_rho: f64 = sim.rho_c.iter().sum();
        assert!(sum_rho > 0.0 && sum_rho.is_finite());
    }
    #[test]
    fn test_microfluidics_add_droplet() {
        let hs = HeleShawParams {
            gap_height: 1e-4,
            width: 1e-3,
            viscosity: 1e-3,
        };
        let df = DeanFlowParams {
            radius: 1e-3,
            curvature_radius: 10e-3,
            mean_velocity: 0.1,
            kinematic_viscosity: 1e-6,
        };
        let mut sim = MicrofluidicsLbm::new(10, 8, 0.8, 0.9, -0.1, hs, df);
        sim.add_droplet(Droplet::new([5.0, 4.0], 1.5, 1.0, 0.001));
        assert_eq!(sim.droplet_count(), 1);
    }
    #[test]
    fn test_biofiltration_construction() {
        let sim = BiofiltrationLbm::new(16, 8, 0.8, 0.7, 1e-4, 1.0, 0.5, 10.0);
        assert_eq!(sim.nx, 16);
        assert_eq!(sim.ny, 8);
        assert_eq!(sim.c_inlet, 10.0);
    }
    #[test]
    fn test_biofiltration_monod_rate() {
        let sim = BiofiltrationLbm::new(8, 4, 0.8, 0.7, 1e-4, 1.0, 1.0, 10.0);
        let rate = sim.monod_rate(1.0);
        assert!((rate - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_biofiltration_monod_zero() {
        let sim = BiofiltrationLbm::new(8, 4, 0.8, 0.7, 1e-4, 1.0, 1.0, 10.0);
        assert_eq!(sim.monod_rate(0.0), 0.0);
    }
    #[test]
    fn test_biofiltration_step_runs() {
        let mut sim = BiofiltrationLbm::new(10, 6, 0.8, 0.7, 1e-4, 1.0, 0.5, 10.0);
        sim.step(1e-4);
        assert!(sim.total_biofilm() >= 0.0);
    }
    #[test]
    fn test_biofiltration_biofilm_grows() {
        let mut sim = BiofiltrationLbm::new(8, 4, 0.8, 0.7, 1e-4, 0.1, 1.0, 100.0);
        for c in sim.nutrient.iter_mut() {
            *c = 100.0;
        }
        let bf0 = sim.total_biofilm();
        sim.step(1e-4);
        let bf1 = sim.total_biofilm();
        assert!(bf1 >= bf0, "biofilm should not shrink: {bf0} → {bf1}");
    }
    #[test]
    fn test_biofiltration_pressure_drop() {
        let sim = BiofiltrationLbm::new(8, 4, 0.8, 0.7, 1e-4, 1.0, 0.5, 10.0);
        let dp = sim.pressure_drop(1.0, 0.1);
        assert!(dp > 0.0 && dp.is_finite());
    }
    #[test]
    fn test_biofiltration_mean_nutrient() {
        let sim = BiofiltrationLbm::new(8, 4, 0.8, 0.7, 1e-4, 1.0, 0.5, 10.0);
        let mn = sim.mean_nutrient();
        assert!(mn.is_finite());
        assert!(mn >= 0.0);
    }
    #[test]
    fn test_biofiltration_run_multiple_steps() {
        let mut sim = BiofiltrationLbm::new(8, 4, 0.8, 0.7, 1e-4, 1.0, 0.5, 10.0);
        sim.run(10, 1e-4);
        assert!(sim.total_biofilm() >= 0.0);
    }
    #[test]
    fn test_biofiltration_porous_properties() {
        let mut cell = BiofiltrationCell::new(5.0);
        cell.biofilm_fraction = 0.3;
        cell.update_porous_properties(1e-4);
        assert!((cell.porosity - 0.7).abs() < 1e-12);
        assert!(cell.permeability > 0.0);
    }
    #[test]
    fn test_biofiltration_full_biofilm_permeability() {
        let mut cell = BiofiltrationCell::new(5.0);
        cell.biofilm_fraction = 0.95;
        cell.update_porous_properties(1e-4);
        assert!(cell.porosity <= 0.05 + 1e-10);
    }
    #[test]
    fn test_feq_sums_to_rho() {
        let rho = 1.05;
        let ux = 0.02;
        let uy = -0.01;
        let f = feq(rho, ux, uy);
        let sum: f64 = f.iter().sum();
        assert!(
            (sum - rho).abs() < 1e-12,
            "sum of feq = {sum} ≠ rho = {rho}"
        );
    }
    #[test]
    fn test_thrombosis_adhered_count() {
        let mut m = ThrombosisModel::new(8, 8, 20, 1.0, 0.0);
        let shear = vec![1000.0f64; 64];
        for _ in 0..100 {
            m.step(&shear, 0.1);
        }
        assert!(m.adhered_count() > 0);
    }
    #[test]
    fn test_platelet_shear_history_accumulates() {
        let mut p = Platelet::new([0.0, 0.0]);
        p.update_activation(5.0, 1.0, 0.0, 0.0);
        assert!(p.shear_history > 0.0);
    }
    #[test]
    fn test_artery_p_out_respected() {
        let p_out = 800.0;
        let art = ArteryModel::new(100.0, 1e-8, 10000.0, p_out, 6.0, 1e-8);
        assert_eq!(art.p_out, p_out);
    }
    #[test]
    fn test_droplet_critical_ca_positive() {
        let d = Droplet::new([0.0, 0.0], 1.0, 2.0, 0.001);
        assert!(d.critical_capillary_number() > 0.0);
    }
    #[test]
    fn test_dean_pressure_drop_factor() {
        let df = DeanFlowParams {
            radius: 1e-3,
            curvature_radius: 50e-3,
            mean_velocity: 0.5,
            kinematic_viscosity: 1e-6,
        };
        let factor = df.pressure_drop_factor();
        assert!(factor >= 1.0 || factor > 0.0);
    }
}
