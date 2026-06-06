//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use crate::porous::types::*;

/// D2Q9 lattice velocities: (cx, cy) for directions 0..9.
pub(super) const CX: [i32; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
pub(super) const CY: [i32; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
/// D2Q9 weights.
pub(super) const W: [f64; 9] = [
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
/// Opposite direction indices for bounce-back.
pub(super) const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
/// Compute the Kozeny-Carman permeability for a packed bed.
///
/// K = phi^3 * d_p^2 / (c_kc * (1-phi)^2)
///
/// where:
///   - phi is the porosity
///   - d_p is the particle diameter
///   - c_kc is the Kozeny-Carman constant (typically 150-180, default 180)
pub fn kozeny_carman_permeability(porosity: f64, particle_diameter: f64, c_kc: f64) -> f64 {
    let phi = porosity.clamp(1e-10, 1.0 - 1e-10);
    let one_minus_phi = 1.0 - phi;
    phi * phi * phi * particle_diameter * particle_diameter / (c_kc * one_minus_phi * one_minus_phi)
}
/// Compute the Kozeny-Carman permeability with the standard constant (180).
pub fn kozeny_carman_standard(porosity: f64, particle_diameter: f64) -> f64 {
    kozeny_carman_permeability(porosity, particle_diameter, 180.0)
}
/// Compute the Ergun equation permeability (c_kc = 150).
pub fn ergun_permeability(porosity: f64, particle_diameter: f64) -> f64 {
    kozeny_carman_permeability(porosity, particle_diameter, 150.0)
}
#[cfg(test)]
mod tests {
    use super::*;
    /// All-solid domain: after many steps, all velocities remain zero.
    #[test]
    fn test_all_solid_no_flow() {
        let nx = 4;
        let ny = 4;
        let mut pm = PorousMediumD2Q9::new(nx, ny, 1.0 / 6.0);
        for j in 0..ny {
            for i in 0..nx {
                pm.set_solid(i, j);
            }
        }
        pm.set_uniform(1.0, [0.1, 0.0]);
        for _ in 0..20 {
            pm.step();
        }
        for idx in 0..(nx * ny) {
            let (_, u) = pm.macros(idx);
            assert!(
                u[0].abs() < 1e-10 && u[1].abs() < 1e-10,
                "Non-zero velocity in all-solid domain at idx={idx}: u={u:?}"
            );
        }
    }
    /// Darcy velocity is proportional to permeability.
    #[test]
    fn test_darcy_velocity_proportional_to_k() {
        let mu = 1e-3;
        let dp = -1.0;
        let df1 = DarcyFlow::new(1e-12, mu);
        let df2 = DarcyFlow::new(2e-12, mu);
        let v1 = df1.velocity(dp);
        let v2 = df2.velocity(dp);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-10,
            "Velocity ratio = {}, expected 2.0",
            v2 / v1
        );
    }
    /// Pressure drop is proportional to length.
    #[test]
    fn test_darcy_pressure_drop_proportional_to_length() {
        let df = DarcyFlow::new(1e-12, 1e-3);
        let v = 1e-4;
        let dp1 = df.pressure_drop(1.0, v);
        let dp2 = df.pressure_drop(2.0, v);
        assert!(
            (dp2 / dp1 - 2.0).abs() < 1e-10,
            "Pressure drop ratio = {}, expected 2.0",
            dp2 / dp1
        );
    }
    /// Kozeny-Carman estimate: all-fluid domain yields high permeability.
    #[test]
    fn test_permeability_all_fluid_high() {
        let pm_fluid = PorousMediumD2Q9::new(8, 8, 1.0 / 6.0);
        let k_fluid = pm_fluid.permeability_estimate();
        assert!(
            k_fluid > 1.0,
            "Expected high K for all-fluid, got {k_fluid}"
        );
        let mut pm_half = PorousMediumD2Q9::new(8, 8, 1.0 / 6.0);
        for j in 0..8 {
            for i in 0..4 {
                pm_half.set_solid(i, j);
            }
        }
        let k_half = pm_half.permeability_estimate();
        assert!(
            k_half < k_fluid,
            "Half-solid K={k_half} should be less than all-fluid K={k_fluid}"
        );
    }
    /// Mass is conserved after one step in a fully-fluid periodic domain.
    #[test]
    fn test_mass_conserved_after_step() {
        let nx = 6;
        let ny = 6;
        let mut pm = PorousMediumD2Q9::new(nx, ny, 1.0 / 6.0);
        let idx22 = pm.idx(2, 2);
        pm.f[idx22] = PorousMediumD2Q9::equilibrium(1.2, [0.05, -0.02]);
        let mass_before: f64 = pm.f.iter().map(|fi| fi.iter().sum::<f64>()).sum();
        pm.step();
        let mass_after: f64 = pm.f.iter().map(|fi| fi.iter().sum::<f64>()).sum();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    /// DBF: Darcy drag is proportional to velocity.
    #[test]
    fn test_dbf_darcy_drag_proportional() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-3, 0.1, 1e-3, 1.0, 0.4);
        let d1 = dbf.darcy_drag(1.0);
        let d2 = dbf.darcy_drag(2.0);
        assert!(
            (d2 / d1 - 2.0).abs() < 1e-10,
            "Darcy drag should be proportional to velocity"
        );
    }
    /// DBF: Forchheimer drag is proportional to velocity squared.
    #[test]
    fn test_dbf_forchheimer_drag_quadratic() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-3, 0.1, 1e-3, 1.0, 0.4);
        let d1 = dbf.forchheimer_drag(1.0);
        let d2 = dbf.forchheimer_drag(2.0);
        assert!(
            (d2 / d1 - 4.0).abs() < 1e-10,
            "Forchheimer drag ratio = {}, expected 4.0",
            d2 / d1
        );
    }
    /// DBF: total resistance is sum of Darcy and Forchheimer.
    #[test]
    fn test_dbf_total_resistance() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-3, 0.1, 1e-3, 1.0, 0.4);
        let v = 0.5;
        let total = dbf.total_resistance(v);
        let darcy = dbf.darcy_drag(v);
        let forch = dbf.forchheimer_drag(v);
        assert!(
            (total - darcy - forch).abs() < 1e-14,
            "Total resistance should equal Darcy + Forchheimer"
        );
    }
    /// DBF: Forchheimer number is zero at zero velocity.
    #[test]
    fn test_dbf_forchheimer_number_zero() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-3, 0.1, 1e-3, 1.0, 0.4);
        let fo = dbf.forchheimer_number(0.0);
        assert!(fo.abs() < 1e-14, "Fo should be zero at u=0: {fo}");
    }
    /// DBF: effective viscosity ratio is 1/porosity.
    #[test]
    fn test_dbf_effective_viscosity_ratio() {
        let phi = 0.4;
        let dbf = DarcyBrinkmanForchheimer::new(1e-3, 0.1, 1e-3, 1.0, phi);
        let ratio = dbf.effective_viscosity_ratio();
        assert!(
            (ratio - 1.0 / phi).abs() < 1e-14,
            "Effective viscosity ratio = {ratio}, expected {}",
            1.0 / phi
        );
    }
    /// DBF: LBM body force opposes flow direction.
    #[test]
    fn test_dbf_body_force_opposes_flow() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-3, 0.1, 1e-3, 1.0, 0.4);
        let [fx, fy] = dbf.lbm_body_force(0.1, 0.05);
        assert!(fx < 0.0, "fx should be negative for positive ux: fx={fx}");
        assert!(fy < 0.0, "fy should be negative for positive uy: fy={fy}");
    }
    /// KC: higher porosity gives higher permeability.
    #[test]
    fn test_kozeny_carman_porosity_effect() {
        let k_low = kozeny_carman_standard(0.3, 1e-3);
        let k_high = kozeny_carman_standard(0.5, 1e-3);
        assert!(
            k_high > k_low,
            "Higher porosity should give higher K: {k_high} vs {k_low}"
        );
    }
    /// KC: permeability scales with d_p^2.
    #[test]
    fn test_kozeny_carman_diameter_scaling() {
        let k1 = kozeny_carman_standard(0.4, 1e-3);
        let k2 = kozeny_carman_standard(0.4, 2e-3);
        assert!(
            (k2 / k1 - 4.0).abs() < 1e-10,
            "K should scale as dp^2: ratio = {}",
            k2 / k1
        );
    }
    /// Ergun permeability is higher than KC standard (smaller constant).
    #[test]
    fn test_ergun_vs_kc() {
        let k_kc = kozeny_carman_standard(0.4, 1e-3);
        let k_ergun = ergun_permeability(0.4, 1e-3);
        assert!(
            k_ergun > k_kc,
            "Ergun K should exceed KC standard: {k_ergun} vs {k_kc}"
        );
    }
    /// PBB: sigma=0 means fully transparent (fluid).
    #[test]
    fn test_pbb_fully_transparent() {
        let pbb = PartialBounceBack::new(0.0);
        let f_pre = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let f_streamed = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0];
        let mut f_out = [0.0_f64; 9];
        pbb.apply(&f_pre, &f_streamed, &mut f_out);
        for k in 0..9 {
            assert!(
                (f_out[k] - f_streamed[k]).abs() < 1e-14,
                "sigma=0 should give streamed distribution"
            );
        }
    }
    /// PBB: sigma=1 means fully solid (bounce-back).
    #[test]
    fn test_pbb_fully_solid() {
        let pbb = PartialBounceBack::new(1.0);
        let f_pre = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let f_streamed = [10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0];
        let mut f_out = [0.0_f64; 9];
        pbb.apply(&f_pre, &f_streamed, &mut f_out);
        for k in 0..9 {
            let expected = f_pre[OPP[k]];
            assert!(
                (f_out[k] - expected).abs() < 1e-14,
                "sigma=1 should give bounce-back: f_out[{k}]={}, expected {expected}",
                f_out[k]
            );
        }
    }
    /// PBB: intermediate sigma gives linear interpolation.
    #[test]
    fn test_pbb_intermediate_sigma() {
        let sigma = 0.3;
        let pbb = PartialBounceBack::new(sigma);
        let f_pre = [1.0; 9];
        let f_streamed = [2.0; 9];
        let mut f_out = [0.0_f64; 9];
        pbb.apply(&f_pre, &f_streamed, &mut f_out);
        for (k, &fo) in f_out.iter().enumerate() {
            let expected = (1.0 - sigma) * 2.0 + sigma * 1.0;
            assert!(
                (fo - expected).abs() < 1e-14,
                "PBB intermediate: f_out[{k}]={fo}, expected {expected}",
            );
        }
    }
    /// PBB: from_permeability produces valid sigma.
    #[test]
    fn test_pbb_from_permeability() {
        let tau = 1.0;
        let k = 0.1;
        let pbb = PartialBounceBack::from_permeability(k, tau);
        assert!(
            pbb.sigma >= 0.0 && pbb.sigma <= 1.0,
            "sigma out of range: {}",
            pbb.sigma
        );
        let k_eff = pbb.effective_permeability(tau);
        assert!(
            (k_eff - k).abs() < 1e-10,
            "Effective K = {k_eff}, expected {k}"
        );
    }
    /// Gray LB: ns=0 everywhere behaves like standard BGK.
    #[test]
    fn test_gray_lb_all_fluid() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let mut glb = GrayLatticeBoltzmann::new(nx, ny, nu);
        let idx = glb.idx(3, 3);
        glb.f[idx] = PorousMediumD2Q9::equilibrium(1.2, [0.05, -0.02]);
        let mass_before = glb.total_mass();
        glb.step();
        let mass_after = glb.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Gray LB mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    /// Gray LB: ns=1 everywhere gives bounce-back-like behaviour.
    #[test]
    fn test_gray_lb_all_solid() {
        let nx = 4;
        let ny = 4;
        let mut glb = GrayLatticeBoltzmann::new(nx, ny, 1.0 / 6.0);
        for j in 0..ny {
            for i in 0..nx {
                glb.set_grayness(i, j, 1.0);
            }
        }
        glb.set_uniform(1.0, [0.1, 0.0]);
        for _ in 0..20 {
            glb.step();
        }
        let avg_u = glb.average_velocity();
        assert!(
            avg_u < 1e-6,
            "All-solid gray LB should have near-zero velocity: {avg_u}"
        );
    }
    /// Gray LB: intermediate ns reduces flow velocity.
    #[test]
    fn test_gray_lb_intermediate_ns() {
        let nx = 10;
        let ny = 10;
        let nu = 1.0 / 6.0;
        let mut glb_fluid = GrayLatticeBoltzmann::new(nx, ny, nu);
        glb_fluid.set_uniform(1.0, [0.05, 0.0]);
        for _ in 0..50 {
            glb_fluid.step();
        }
        let u_fluid = glb_fluid.average_velocity();
        let mut glb_gray = GrayLatticeBoltzmann::new(nx, ny, nu);
        for j in 0..ny {
            for i in 0..nx {
                glb_gray.set_grayness(i, j, 0.3);
            }
        }
        glb_gray.set_uniform(1.0, [0.05, 0.0]);
        for _ in 0..50 {
            glb_gray.step();
        }
        let u_gray = glb_gray.average_velocity();
        assert!(
            u_gray < u_fluid + 1e-6,
            "Gray medium should have lower velocity: u_gray={u_gray}, u_fluid={u_fluid}"
        );
    }
    /// Gray LB: set_grayness clamps values.
    #[test]
    fn test_gray_lb_grayness_clamping() {
        let mut glb = GrayLatticeBoltzmann::new(4, 4, 1.0 / 6.0);
        glb.set_grayness(0, 0, -0.5);
        glb.set_grayness(1, 1, 1.5);
        assert!(
            (glb.ns[glb.idx(0, 0)] - 0.0).abs() < 1e-14,
            "Negative ns should clamp to 0"
        );
        assert!(
            (glb.ns[glb.idx(1, 1)] - 1.0).abs() < 1e-14,
            "ns > 1 should clamp to 1"
        );
    }
    /// Average velocity is zero for all-solid domain.
    #[test]
    fn test_average_velocity_all_solid() {
        let mut pm = PorousMediumD2Q9::new(4, 4, 1.0 / 6.0);
        for j in 0..4 {
            for i in 0..4 {
                pm.set_solid(i, j);
            }
        }
        pm.set_uniform(1.0, [0.1, 0.0]);
        for _ in 0..20 {
            pm.step();
        }
        let avg = pm.average_velocity();
        assert!(
            avg < 1e-10,
            "All-solid average velocity should be ~0: {avg}"
        );
    }
    /// Porous nodes reduce velocity relative to fluid nodes.
    #[test]
    fn test_porous_reduces_velocity() {
        let nx = 10;
        let ny = 10;
        let nu = 1.0 / 6.0;
        let mut pm_fluid = PorousMediumD2Q9::new(nx, ny, nu);
        pm_fluid.set_uniform(1.0, [0.05, 0.0]);
        for _ in 0..50 {
            pm_fluid.step();
        }
        let u_fluid = pm_fluid.average_velocity();
        let mut pm_porous = PorousMediumD2Q9::new(nx, ny, nu);
        for j in 0..ny {
            for i in 0..nx {
                pm_porous.set_porous(i, j, 0.5);
            }
        }
        pm_porous.set_uniform(1.0, [0.05, 0.0]);
        for _ in 0..50 {
            pm_porous.step();
        }
        let u_porous = pm_porous.average_velocity();
        assert!(
            u_porous < u_fluid + 1e-6,
            "Porous should have lower velocity: u_porous={u_porous}, u_fluid={u_fluid}"
        );
    }
    /// Reynolds number is proportional to velocity.
    #[test]
    fn test_darcy_reynolds_proportional_to_velocity() {
        let df = DarcyFlow::new(1e-12, 1e-3);
        let re1 = df.reynolds_number(1.0, 1e-4, 1000.0);
        let re2 = df.reynolds_number(2.0, 1e-4, 1000.0);
        assert!(
            (re2 / re1 - 2.0).abs() < 1e-10,
            "Re should be proportional to velocity"
        );
    }
}
#[cfg(test)]
mod tests_extended {
    use super::*;
    #[test]
    fn test_rev_volume() {
        let rev = RepresentativeElementaryVolume::new(1e-3, 0.4, 100, 50e-6);
        let v = rev.volume();
        assert!((v - 1e-9).abs() < 1e-20, "REV volume = {v}");
    }
    #[test]
    fn test_rev_pore_and_solid_volumes_sum_to_total() {
        let rev = RepresentativeElementaryVolume::new(1e-3, 0.35, 80, 40e-6);
        let total = rev.pore_volume() + rev.solid_volume();
        let v = rev.volume();
        assert!(
            (total - v).abs() < 1e-25,
            "Pore + solid should equal REV volume"
        );
    }
    #[test]
    fn test_rev_hydraulic_diameter_positive() {
        let rev = RepresentativeElementaryVolume::new(1e-3, 0.4, 1000, 20e-6);
        let dh = rev.hydraulic_diameter();
        assert!(dh > 0.0, "Hydraulic diameter should be positive: {dh}");
    }
    #[test]
    fn test_rev_is_sufficient() {
        let rev_small = RepresentativeElementaryVolume::new(100e-6, 0.4, 10, 50e-6);
        let rev_large = RepresentativeElementaryVolume::new(1e-3, 0.4, 10, 50e-6);
        assert!(
            !rev_small.is_sufficient(10.0),
            "Small REV should not be sufficient"
        );
        assert!(
            rev_large.is_sufficient(10.0),
            "Large REV should be sufficient"
        );
    }
    #[test]
    fn test_brinkman_inverse_porosity() {
        let mu = 1e-3;
        let phi = 0.5;
        let mu_eff = BrinkmanEffectiveViscosity::inverse_porosity(mu, phi);
        assert!(
            (mu_eff - 2e-3).abs() < 1e-16,
            "mu_eff should be 2x mu: {mu_eff}"
        );
    }
    #[test]
    fn test_brinkman_viscosity_ratio_gt_one() {
        let ratio = BrinkmanEffectiveViscosity::viscosity_ratio(0.4);
        assert!(ratio > 1.0, "Ratio should exceed 1 for phi < 1: {ratio}");
    }
    #[test]
    fn test_lundgren_increases_for_lower_porosity() {
        let mu_high_phi = BrinkmanEffectiveViscosity::lundgren(1e-3, 0.8);
        let mu_low_phi = BrinkmanEffectiveViscosity::lundgren(1e-3, 0.3);
        assert!(
            mu_low_phi > mu_high_phi,
            "Lower porosity should give higher mu_eff"
        );
    }
    #[test]
    fn test_forchheimer_number_zero_at_rest() {
        let ff = ForchheimerFlow::new(1e-8, 0.1, 1e-3, 1000.0);
        let fo = ff.forchheimer_number(0.0);
        assert!(fo.abs() < 1e-20, "Fo should be 0 at rest: {fo}");
    }
    #[test]
    fn test_forchheimer_pressure_gradient_sign() {
        let ff = ForchheimerFlow::new(1e-8, 0.1, 1e-3, 1000.0);
        let dp = ff.pressure_gradient(0.01);
        assert!(
            dp > 0.0,
            "Pressure gradient should be positive for positive u: {dp}"
        );
    }
    #[test]
    fn test_forchheimer_pressure_drop_proportional_to_length() {
        let ff = ForchheimerFlow::new(1e-8, 0.1, 1e-3, 1000.0);
        let dp1 = ff.pressure_drop(0.01, 1.0);
        let dp2 = ff.pressure_drop(0.01, 2.0);
        assert!(
            (dp2 / dp1 - 2.0).abs() < 1e-10,
            "Pressure drop proportional to length"
        );
    }
    #[test]
    fn test_forchheimer_ergun_permeability_scaling() {
        let k1 = ForchheimerFlow::ergun_permeability(0.4, 1e-3);
        let k2 = ForchheimerFlow::ergun_permeability(0.4, 2e-3);
        assert!(
            (k2 / k1 - 4.0).abs() < 1e-10,
            "Ergun K scales as d_p^2: {}",
            k2 / k1
        );
    }
    #[test]
    fn test_forchheimer_body_force_opposes_flow() {
        let bf = ForchheimerBodyForce::new(1e-3, 0.1, 1e-6, 1.0);
        let [fx, fy] = bf.force(0.1, 0.05);
        assert!(fx < 0.0, "Force x should oppose positive ux: {fx}");
        assert!(fy < 0.0, "Force y should oppose positive uy: {fy}");
    }
    #[test]
    fn test_forchheimer_body_force_zero_at_rest() {
        let bf = ForchheimerBodyForce::new(1e-3, 0.1, 1e-6, 1.0);
        let [fx, fy] = bf.force(0.0, 0.0);
        assert!(fx.abs() < 1e-20, "Force x should be zero at rest: {fx}");
        assert!(fy.abs() < 1e-20, "Force y should be zero at rest: {fy}");
    }
    #[test]
    fn test_guo_forcing_sums_to_net_force() {
        let bf = ForchheimerBodyForce::new(1e-3, 0.1, 1e-6, 1.0);
        let ux = 0.05;
        let uy = 0.02;
        let omega = 1.0;
        let fi_force = bf.guo_forcing_d2q9(ux, uy, omega);
        let [fx, fy] = bf.force(ux, uy);
        let sum_x: f64 = fi_force
            .iter()
            .enumerate()
            .map(|(k, &f)| f * CX[k] as f64)
            .sum();
        let sum_y: f64 = fi_force
            .iter()
            .enumerate()
            .map(|(k, &f)| f * CY[k] as f64)
            .sum();
        assert!(
            (sum_x - fx * (1.0 - 0.5 * omega)).abs() < 1e-12,
            "Guo x-force mismatch: {sum_x} vs {}",
            fx * (1.0 - 0.5 * omega)
        );
        assert!(
            (sum_y - fy * (1.0 - 0.5 * omega)).abs() < 1e-12,
            "Guo y-force mismatch: {sum_y} vs {}",
            fy * (1.0 - 0.5 * omega)
        );
    }
    #[test]
    fn test_dbf_lbm_mass_conservation() {
        let nx = 6;
        let ny = 6;
        let bf = ForchheimerBodyForce::new(0.01, 0.05, 1.0 / 6.0, 1.0);
        let mut lbm = DbfLbmD2Q9::new(nx, ny, 1.0 / 6.0).with_body_force(bf);
        lbm.set_uniform(1.0, [0.02, 0.0]);
        for j in 2..4 {
            for i in 2..4 {
                lbm.set_porous(i, j, 0.6);
            }
        }
        let mass_before = lbm.total_mass();
        for _ in 0..10 {
            lbm.step();
        }
        let mass_after = lbm.total_mass();
        assert!(
            (mass_before - mass_after).abs() / mass_before < 1e-6,
            "DBF mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_dbf_lbm_average_velocity_positive() {
        let nx = 8;
        let ny = 8;
        let mut lbm = DbfLbmD2Q9::new(nx, ny, 1.0 / 6.0);
        lbm.set_uniform(1.0, [0.05, 0.0]);
        for _ in 0..5 {
            lbm.step();
        }
        let u_avg = lbm.average_velocity();
        assert!(u_avg > 0.0, "Average velocity should be positive: {u_avg}");
    }
    #[test]
    fn test_permeability_tensor_isotropic_darcy() {
        let k = PermeabilityTensor2D::isotropic(1e-12);
        let mu = 1e-3;
        let [ux, uy] = k.darcy_velocity(-1000.0, 0.0, mu);
        let expected = 1e-12 / 1e-3 * 1000.0;
        assert!(
            (ux - expected).abs() < 1e-20,
            "ux = {ux}, expected {expected}"
        );
        assert!(uy.abs() < 1e-20, "uy should be 0: {uy}");
    }
    #[test]
    fn test_permeability_tensor_determinant() {
        let k = PermeabilityTensor2D::isotropic(2.0);
        let det = k.determinant();
        assert!(
            (det - 4.0).abs() < 1e-14,
            "Isotropic det should be k^2: {det}"
        );
    }
    #[test]
    fn test_permeability_tensor_harmonic_mean_isotropic() {
        let k = PermeabilityTensor2D::isotropic(3.0);
        let hm = k.harmonic_mean();
        assert!(
            (hm - 3.0).abs() < 1e-10,
            "Harmonic mean of isotropic should equal k: {hm}"
        );
    }
    #[test]
    fn test_porous_channel_velocity_zero_at_wall() {
        let pcf = PorousChannelFlow::new(1e-3, 1e-3, 1e-8, 0.01);
        let dp_dx = -1000.0;
        let u_wall = pcf.velocity_profile(pcf.half_width, dp_dx);
        assert!(
            u_wall.abs() < 1e-10,
            "Velocity at wall should be ~0: {u_wall}"
        );
    }
    #[test]
    fn test_porous_channel_average_velocity_positive() {
        let pcf = PorousChannelFlow::new(1e-3, 1e-3, 1e-8, 0.01);
        let dp_dx = -1000.0;
        let u_avg = pcf.average_velocity(dp_dx);
        assert!(
            u_avg > 0.0,
            "Average velocity should be positive for dp<0: {u_avg}"
        );
    }
    #[test]
    fn test_porous_channel_effective_permeability_bounded() {
        let k = 1e-8;
        let pcf = PorousChannelFlow::new(1e-3, 1e-3, k, 0.01);
        let k_eff = pcf.effective_permeability(-1000.0);
        assert!(
            k_eff <= k * 1.01,
            "Effective K should be <= K: {k_eff} vs {k}"
        );
    }
    #[test]
    fn test_bruggeman_tortuosity_greater_than_one() {
        let tau = TortuosityModel::bruggeman(0.4);
        assert!(tau > 1.0, "Bruggeman tortuosity should be > 1: {tau}");
    }
    #[test]
    fn test_millington_quirk_tortuosity() {
        let tau = TortuosityModel::millington_quirk(0.4);
        assert!(
            tau > 1.0,
            "Millington-Quirk tortuosity should be > 1: {tau}"
        );
    }
    #[test]
    fn test_weissberg_tortuosity() {
        let tau = TortuosityModel::weissberg(0.4);
        assert!(tau > 1.0, "Weissberg tortuosity should be > 1: {tau}");
    }
    #[test]
    fn test_effective_diffusivity_less_than_free() {
        let d_free = 1e-9;
        let phi = 0.4;
        let tau = TortuosityModel::bruggeman(phi);
        let d_eff = TortuosityModel::effective_diffusivity(d_free, phi, tau);
        assert!(
            d_eff < d_free,
            "Effective diffusivity should be less than free: {d_eff} vs {d_free}"
        );
    }
    #[test]
    fn test_bruggeman_diffusivity_scales_with_phi() {
        let d = 1e-9;
        let d_eff_high = TortuosityModel::bruggeman_diffusivity(d, 0.8);
        let d_eff_low = TortuosityModel::bruggeman_diffusivity(d, 0.3);
        assert!(d_eff_high > d_eff_low, "Higher porosity gives higher D_eff");
    }
    #[test]
    fn test_bruggeman_diffusivity_at_phi_one_equals_d() {
        let d = 1.5e-9;
        let d_eff = TortuosityModel::bruggeman_diffusivity(d, 1.0);
        assert!(
            (d_eff - d).abs() < 1e-20,
            "At phi=1, D_eff should equal D: {d_eff}"
        );
    }
}
#[cfg(test)]
mod tests_porous_extended {
    use super::*;
    #[test]
    fn test_rev_volume_cubic() {
        let rev = RepresentativeElementaryVolume::new(0.01, 0.4, 100, 1e-4);
        let v = rev.volume();
        assert!((v - 1e-6).abs() < 1e-18, "REV volume = {v}");
    }
    #[test]
    fn test_rev_pore_solid_sum_to_total() {
        let rev = RepresentativeElementaryVolume::new(0.01, 0.4, 100, 1e-4);
        let total = rev.pore_volume() + rev.solid_volume();
        assert!(
            (total - rev.volume()).abs() < 1e-20,
            "pore + solid = total: {total}"
        );
    }
    #[test]
    fn test_rev_pore_volume_proportional_to_porosity() {
        let rev1 = RepresentativeElementaryVolume::new(0.01, 0.3, 100, 1e-4);
        let rev2 = RepresentativeElementaryVolume::new(0.01, 0.6, 100, 1e-4);
        assert!(
            rev2.pore_volume() > rev1.pore_volume(),
            "Higher porosity → larger pore volume"
        );
    }
    #[test]
    fn test_rev_is_sufficient_true() {
        let rev = RepresentativeElementaryVolume::new(0.01, 0.4, 100, 1e-4);
        assert!(
            rev.is_sufficient(10.0),
            "REV with L/r=100 should be sufficient"
        );
    }
    #[test]
    fn test_rev_is_sufficient_false_for_small_rev() {
        let rev = RepresentativeElementaryVolume::new(1e-4, 0.4, 100, 1e-4);
        assert!(
            !rev.is_sufficient(10.0),
            "REV with L/r=1 should not be sufficient"
        );
    }
    #[test]
    fn test_rev_hydraulic_diameter_positive() {
        let rev = RepresentativeElementaryVolume::new(0.01, 0.4, 1000, 1e-4);
        let dh = rev.hydraulic_diameter();
        assert!(dh > 0.0, "Hydraulic diameter = {dh}");
    }
    #[test]
    fn test_brinkman_inverse_porosity_at_phi_one() {
        let mu = 1e-3;
        let mu_eff = BrinkmanEffectiveViscosity::inverse_porosity(mu, 1.0);
        assert!((mu_eff - mu).abs() < 1e-20, "mu_eff at phi=1 = {mu_eff}");
    }
    #[test]
    fn test_brinkman_inverse_porosity_increases_for_lower_phi() {
        let mu_eff1 = BrinkmanEffectiveViscosity::inverse_porosity(1e-3, 0.8);
        let mu_eff2 = BrinkmanEffectiveViscosity::inverse_porosity(1e-3, 0.4);
        assert!(
            mu_eff2 > mu_eff1,
            "mu_eff increases as phi decreases: {mu_eff1}, {mu_eff2}"
        );
    }
    #[test]
    fn test_brinkman_lundgren_at_phi_one_is_mu() {
        let mu = 1e-3;
        let mu_eff = BrinkmanEffectiveViscosity::lundgren(mu, 1.0);
        assert!((mu_eff - mu).abs() < 1e-20, "Lundgren at phi=1: {mu_eff}");
    }
    #[test]
    fn test_brinkman_viscosity_ratio_at_phi_one() {
        let ratio = BrinkmanEffectiveViscosity::viscosity_ratio(1.0);
        assert!(
            (ratio - 1.0).abs() < 1e-14,
            "viscosity ratio at phi=1 = {ratio}"
        );
    }
    #[test]
    fn test_forchheimer_number_zero_at_rest() {
        let ff = ForchheimerFlow::new(1e-10, 0.1, 1e-3, 1000.0);
        let fo = ff.forchheimer_number(0.0);
        assert!(fo.abs() < 1e-30, "Forchheimer number at rest = {fo}");
    }
    #[test]
    fn test_forchheimer_pressure_gradient_positive_for_positive_u() {
        let ff = ForchheimerFlow::new(1e-10, 0.1, 1e-3, 1000.0);
        let dp_dx = ff.pressure_gradient(0.01);
        assert!(dp_dx > 0.0, "Pressure gradient should be positive: {dp_dx}");
    }
    #[test]
    fn test_forchheimer_pressure_drop_proportional_to_length() {
        let ff = ForchheimerFlow::new(1e-10, 0.1, 1e-3, 1000.0);
        let dp1 = ff.pressure_drop(0.01, 1.0);
        let dp2 = ff.pressure_drop(0.01, 2.0);
        assert!(
            (dp2 / dp1 - 2.0).abs() < 1e-10,
            "Pressure drop ∝ L: {dp1}, {dp2}"
        );
    }
    #[test]
    fn test_forchheimer_ergun_permeability_increases_with_porosity() {
        let k1 = ForchheimerFlow::ergun_permeability(0.3, 1e-3);
        let k2 = ForchheimerFlow::ergun_permeability(0.6, 1e-3);
        assert!(k2 > k1, "Ergun permeability increases with phi: {k1}, {k2}");
    }
    #[test]
    fn test_forchheimer_inertial_coefficient_ward_decreases_with_phi() {
        let c1 = ForchheimerFlow::inertial_coefficient_ward(0.2);
        let c2 = ForchheimerFlow::inertial_coefficient_ward(0.8);
        assert!(c1 > c2, "Ward coefficient decreases with phi: {c1}, {c2}");
    }
    #[test]
    fn test_permeability_tensor_isotropic_determinant() {
        let k = 1e-10;
        let pt = PermeabilityTensor2D::isotropic(k);
        let det = pt.determinant();
        assert!((det - k * k).abs() < 1e-30, "Isotropic tensor det = {det}");
    }
    #[test]
    fn test_permeability_tensor_isotropic_kxx_equals_kyy() {
        let pt = PermeabilityTensor2D::isotropic(1e-10);
        assert!(
            (pt.kxx - pt.kyy).abs() < 1e-30,
            "Isotropic: kxx should equal kyy"
        );
    }
    #[test]
    fn test_permeability_tensor_darcy_velocity_direction() {
        let pt = PermeabilityTensor2D::isotropic(1e-10);
        let u = pt.darcy_velocity(-1000.0, 0.0, 1e-3);
        assert!(u[0] > 0.0, "Darcy velocity in +x for -dP/dx: {}", u[0]);
    }
    #[test]
    fn test_permeability_tensor_geometric_mean_isotropic() {
        let k = 1e-10;
        let pt = PermeabilityTensor2D::isotropic(k);
        let gm = pt.geometric_mean();
        assert!(
            (gm - k).abs() < 1e-30,
            "Geometric mean of isotropic tensor = {gm}"
        );
    }
    #[test]
    fn test_porous_channel_velocity_zero_at_wall() {
        let h = 0.005;
        let ch = PorousChannelFlow::new(1e-3, 1e-3, 1e-10, h);
        let u_wall = ch.velocity_profile(h, -1000.0);
        assert!(u_wall.abs() < 1e-10, "Velocity at wall = {u_wall}");
    }
    #[test]
    fn test_porous_channel_velocity_maximum_at_centerline() {
        let h = 0.005;
        let ch = PorousChannelFlow::new(1e-3, 1e-3, 1e-10, h);
        let dp_dx = -1000.0;
        let u_center = ch.velocity_profile(0.0, dp_dx);
        let u_quarter = ch.velocity_profile(h / 2.0, dp_dx);
        assert!(
            u_center >= u_quarter,
            "Max velocity at centerline: {u_center}, {u_quarter}"
        );
    }
    #[test]
    fn test_porous_channel_average_velocity_positive() {
        let h = 0.005;
        let ch = PorousChannelFlow::new(1e-3, 1e-3, 1e-10, h);
        let u_avg = ch.average_velocity(-1000.0);
        assert!(u_avg > 0.0, "Average velocity should be positive: {u_avg}");
    }
    #[test]
    fn test_kozeny_carman_zero_permeability_at_zero_porosity() {
        let k = kozeny_carman_permeability(0.0, 1e-3, 180.0);
        assert!(k < 1e-30, "Kozeny-Carman at phi=0 should be ~0: {k}");
    }
    #[test]
    fn test_kozeny_carman_standard_equals_c180() {
        let k1 = kozeny_carman_standard(0.4, 1e-3);
        let k2 = kozeny_carman_permeability(0.4, 1e-3, 180.0);
        assert!((k1 - k2).abs() < 1e-30, "standard == c_kc=180: {k1}, {k2}");
    }
    #[test]
    fn test_ergun_permeability_uses_c150() {
        let k_ergun = ergun_permeability(0.4, 1e-3);
        let k_c150 = kozeny_carman_permeability(0.4, 1e-3, 150.0);
        assert!(
            (k_ergun - k_c150).abs() < 1e-30,
            "Ergun == c_kc=150: {k_ergun}, {k_c150}"
        );
    }
    #[test]
    fn test_gray_lb_mass_conserved_step() {
        let mut glb = GrayLatticeBoltzmann::new(4, 4, 0.1);
        glb.set_uniform(1.0, [0.0; 2]);
        let mass_before = glb.total_mass();
        glb.step();
        let mass_after = glb.total_mass();
        assert!(
            (mass_after - mass_before).abs() < 1e-10,
            "GLB mass conservation: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_gray_lb_grayness_zero_all_solid() {
        let mut glb = GrayLatticeBoltzmann::new(4, 4, 0.1);
        for j in 0..4 {
            for i in 0..4 {
                glb.set_grayness(i, j, 1.0);
            }
        }
        let u_avg = glb.average_velocity();
        assert!(u_avg.abs() < 1e-10, "All-solid GLB avg velocity = {u_avg}");
    }
    #[test]
    fn test_effective_diffusivity_less_than_free_diffusivity() {
        let d_free = 1e-9;
        let phi = 0.4;
        let tau = TortuosityModel::bruggeman(phi);
        let d_eff = TortuosityModel::effective_diffusivity(d_free, phi, tau);
        assert!(
            d_eff < d_free,
            "Effective diffusivity < free diffusivity: {d_eff}"
        );
    }
    #[test]
    fn test_weissberg_tortuosity_bounded() {
        let tau = TortuosityModel::weissberg(0.4);
        assert!(tau >= 1.0, "Weissberg tortuosity >= 1: {tau}");
        assert!(tau < 10.0, "Weissberg tortuosity < 10 for phi=0.4: {tau}");
    }
    #[test]
    fn test_millington_quirk_tortuosity_monotone() {
        let tau1 = TortuosityModel::millington_quirk(0.7);
        let tau2 = TortuosityModel::millington_quirk(0.3);
        assert!(
            tau2 > tau1,
            "MQ tortuosity increases with lower phi: {tau1}, {tau2}"
        );
    }
    #[test]
    fn test_porous_media_kozeny_carman_positive() {
        let pm = PorousMedia::new(0.4, 1e-3, 1e-9);
        let k = pm.compute_kozeny_carman(0.0);
        assert!(k > 0.0, "Kozeny-Carman permeability must be positive: {k}");
    }
    #[test]
    fn test_porous_media_kozeny_carman_increases_with_porosity() {
        let pm1 = PorousMedia::new(0.3, 1e-3, 1e-9);
        let pm2 = PorousMedia::new(0.5, 1e-3, 1e-9);
        let k1 = pm1.compute_kozeny_carman(0.0);
        let k2 = pm2.compute_kozeny_carman(0.0);
        assert!(k2 > k1, "Higher porosity → higher K: k1={k1}, k2={k2}");
    }
    #[test]
    fn test_porous_media_kozeny_carman_scales_with_dp_squared() {
        let pm1 = PorousMedia::new(0.4, 1e-3, 1e-9);
        let pm2 = PorousMedia::new(0.4, 2e-3, 1e-9);
        let k1 = pm1.compute_kozeny_carman(0.0);
        let k2 = pm2.compute_kozeny_carman(0.0);
        assert!(
            (k2 / k1 - 4.0).abs() < 1e-10,
            "Doubling dp should give 4x permeability: ratio = {}",
            k2 / k1
        );
    }
    #[test]
    fn test_porous_media_kozeny_carman_zero_at_zero_porosity() {
        let pm = PorousMedia::new(0.0, 1e-3, 1e-9);
        let k = pm.compute_kozeny_carman(0.0);
        assert_eq!(k, 0.0, "Zero porosity → zero permeability");
    }
    #[test]
    fn test_porous_media_forchheimer_correction_exceeds_darcy() {
        let pm = PorousMedia::new(0.4, 1e-3, 1e-9);
        let u = 0.1;
        let mu = 1e-3;
        let rho = 1000.0;
        let dp_forchheimer = pm.compute_forchheimer_correction(u, mu, rho, 0.0);
        let k = pm.compute_kozeny_carman(0.0);
        let dp_darcy = mu / k * u;
        assert!(
            dp_forchheimer > dp_darcy,
            "Forchheimer correction should exceed Darcy: {dp_forchheimer} vs {dp_darcy}"
        );
    }
    #[test]
    fn test_porous_media_forchheimer_correction_nonnegative() {
        let pm = PorousMedia::new(0.4, 1e-3, 1e-9);
        for &u in &[0.0, 0.01, 0.1, 1.0] {
            let dp = pm.compute_forchheimer_correction(u, 1e-3, 1000.0, 0.0);
            assert!(
                dp >= 0.0,
                "Forchheimer correction must be >= 0 for u={u}: {dp}"
            );
        }
    }
    #[test]
    fn test_porous_media_effective_diffusivity_bruggeman() {
        let pm = PorousMedia::new(0.4, 1e-3, 1e-9);
        let d_eff = pm.compute_effective_diffusivity("bruggeman", 1.0);
        let expected = 1e-9 * 0.4f64.powf(1.5);
        assert!(
            (d_eff - expected).abs() < 1e-22,
            "Bruggeman D_eff mismatch: {d_eff} vs {expected}"
        );
    }
    #[test]
    fn test_porous_media_effective_diffusivity_less_than_free() {
        let pm = PorousMedia::new(0.4, 1e-3, 1e-9);
        for model in &["bruggeman", "millington_quirk", "weissberg"] {
            let d_eff = pm.compute_effective_diffusivity(model, 1.0);
            assert!(
                d_eff < 1e-9,
                "D_eff ({model}) should be < free diffusivity: {d_eff}"
            );
        }
    }
    #[test]
    fn test_porous_media_effective_diffusivity_constant_tortuosity() {
        let pm = PorousMedia::new(0.4, 1e-3, 1e-9);
        let tau = 2.5;
        let d_eff = pm.compute_effective_diffusivity("constant", tau);
        let expected = TortuosityModel::effective_diffusivity(1e-9, 0.4, tau);
        assert!(
            (d_eff - expected).abs() < 1e-25,
            "Constant tortuosity D_eff mismatch: {d_eff} vs {expected}"
        );
    }
    #[test]
    fn test_porous_media_effective_diffusivity_increases_with_porosity() {
        let pm1 = PorousMedia::new(0.2, 1e-3, 1e-9);
        let pm2 = PorousMedia::new(0.6, 1e-3, 1e-9);
        let d1 = pm1.compute_effective_diffusivity("bruggeman", 1.0);
        let d2 = pm2.compute_effective_diffusivity("bruggeman", 1.0);
        assert!(d2 > d1, "Higher porosity → higher D_eff: d1={d1}, d2={d2}");
    }
}
