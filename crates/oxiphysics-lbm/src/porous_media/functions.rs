//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::PorousLbmGrid;

/// D2Q9 velocity components (cx, cy) for directions 0..9
pub(super) const CX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
pub(super) const CY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];
/// D2Q9 weights
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
/// Speed of sound squared for D2Q9: cs^2 = 1/3
pub(super) const CS2: f64 = 1.0 / 3.0;
/// Gravity body force in a partially-saturated porous medium.
///
/// F_z = -water_sat * rho * g  (downward gravity component only)
///
/// Returns a 3-D vector with the gravity force in the z-direction.
pub fn partial_saturation_force(water_sat: f64, rho: f64, g: f64) -> [f64; 3] {
    [0.0, 0.0, -water_sat * rho * g]
}
/// Kozeny-Carman permeability correlation.
///
/// K = d_p^2 * eps^3 / (180 * (1 - eps)^2)
///
/// # Arguments
/// * `d_p`     - particle diameter (m)
/// * `epsilon` - porosity (dimensionless, in (0, 1))
pub fn kozeny_carman_permeability(d_p: f64, epsilon: f64) -> f64 {
    let one_minus_eps = 1.0 - epsilon;
    if one_minus_eps.abs() < 1e-12 {
        return f64::MAX;
    }
    d_p * d_p * epsilon.powi(3) / (180.0 * one_minus_eps * one_minus_eps)
}
/// Ergun equation pressure drop over a packed bed.
///
/// ΔP/L = 150 * mu * (1-eps)^2 / (d_p^2 * eps^3) * u
///       + 1.75 * rho * (1-eps) / (d_p * eps^3) * u^2
///
/// # Arguments
/// * `l`       - bed length (m)
/// * `u`       - superficial velocity (m/s)
/// * `d_p`     - particle diameter (m)
/// * `epsilon` - porosity (dimensionless)
/// * `mu`      - dynamic viscosity (Pa·s)
/// * `rho`     - fluid density (kg/m^3)
pub fn ergun_pressure_drop(l: f64, u: f64, d_p: f64, epsilon: f64, mu: f64, rho: f64) -> f64 {
    let one_minus_eps = 1.0 - epsilon;
    let eps3 = epsilon.powi(3);
    let viscous = 150.0 * mu * one_minus_eps * one_minus_eps / (d_p * d_p * eps3) * u;
    let inertial = 1.75 * rho * one_minus_eps / (d_p * eps3) * u * u;
    (viscous + inertial) * l
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::porous_media::types::*;
    #[test]
    fn test_darcy_velocity_proportional_to_pressure_gradient() {
        let k = 1e-10_f64;
        let mu = 1e-3_f64;
        let grad1 = [1000.0_f64, 0.0_f64];
        let grad2 = [2000.0_f64, 0.0_f64];
        let u1 = DarcyFlow::darcy_velocity(k, mu, grad1);
        let u2 = DarcyFlow::darcy_velocity(k, mu, grad2);
        assert!(
            (u2[0] - 2.0 * u1[0]).abs() < 1e-20,
            "Darcy velocity not proportional: u1={}, u2={}",
            u1[0],
            u2[0]
        );
        assert!(
            u1[0] < 0.0,
            "Darcy velocity should be negative for positive gradient"
        );
    }
    #[test]
    fn test_kozeny_carman_permeability_increases_with_porosity() {
        let model = KozenyCarmanModel::new(1e-3);
        let k1 = model.permeability(0.3);
        let k2 = model.permeability(0.5);
        let k3 = model.permeability(0.7);
        assert!(k1 < k2 && k2 < k3, "k(0.3)={k1}, k(0.5)={k2}, k(0.7)={k3}");
    }
    #[test]
    fn test_kozeny_carman_permeability_high_at_unit_porosity() {
        let model = KozenyCarmanModel::new(1e-3);
        let k_high = model.permeability(0.9999);
        let k_low = model.permeability(0.3);
        assert!(k_high > k_low * 1000.0, "k_high={k_high}, k_low={k_low}");
    }
    #[test]
    fn test_porous_lbm_grid_new_correct_cell_count() {
        let nx = 8;
        let ny = 6;
        let grid = PorousLbmGrid::new(nx, ny);
        assert_eq!(grid.cells.len(), nx * ny);
    }
    #[test]
    fn test_set_porosity_field_assigns_values() {
        let nx = 4;
        let ny = 4;
        let mut grid = PorousLbmGrid::new(nx, ny);
        grid.set_porosity_field(|x, y| if (x + y) % 2 == 0 { 0.6 } else { 0.4 });
        for y in 0..ny {
            for x in 0..nx {
                let idx = grid.idx(x, y);
                let expected = if (x + y) % 2 == 0 { 0.6 } else { 0.4 };
                assert!(
                    (grid.cells[idx].porosity - expected).abs() < 1e-15,
                    "at ({x},{y}): expected {expected}, got {}",
                    grid.cells[idx].porosity
                );
            }
        }
    }
    #[test]
    fn test_effective_thermal_conductivity_between_fluid_and_solid() {
        let k_f = 0.6_f64;
        let k_s = 2.0_f64;
        let eps = 0.4_f64;
        let k_eff = EffectiveMediumProperties::effective_thermal_conductivity(k_f, k_s, eps);
        assert!(
            k_eff >= k_f.min(k_s) - 1e-12 && k_eff <= k_f.max(k_s) + 1e-12,
            "k_eff={k_eff} between k_f={k_f} and k_s={k_s}"
        );
    }
    #[test]
    fn test_effective_diffusivity_le_free_diffusivity() {
        let d_free = 1e-9_f64;
        let porosity = 0.4_f64;
        let tortuosity = EffectiveMediumProperties::tortuosity_estimate(porosity);
        let d_eff = EffectiveMediumProperties::effective_diffusivity(d_free, porosity, tortuosity);
        assert!(d_eff <= d_free + 1e-20, "d_eff={d_eff} > d_free={d_free}");
    }
    #[test]
    fn test_isotropic_tensor_determinant() {
        let t = PermeabilityTensor::isotropic(1e-10);
        let det = t.determinant();
        let expected = 1e-10 * 1e-10;
        assert!((det - expected).abs() < 1e-30, "det = {det}");
    }
    #[test]
    fn test_isotropic_tensor_darcy_velocity() {
        let t = PermeabilityTensor::isotropic(1e-10);
        let u = t.darcy_velocity(1e-3, [1000.0, 0.0]);
        let expected = DarcyFlow::darcy_velocity(1e-10, 1e-3, [1000.0, 0.0]);
        assert!((u[0] - expected[0]).abs() < 1e-15, "ux mismatch");
        assert!((u[1] - expected[1]).abs() < 1e-15, "uy mismatch");
    }
    #[test]
    fn test_anisotropic_tensor_principal_perms() {
        let t = PermeabilityTensor::anisotropic(2e-10, 0.0, 1e-10);
        let (k_min, k_max) = t.principal_permeabilities();
        assert!((k_min - 1e-10).abs() < 1e-22, "k_min = {k_min}");
        assert!((k_max - 2e-10).abs() < 1e-22, "k_max = {k_max}");
    }
    #[test]
    fn test_anisotropy_ratio_isotropic() {
        let t = PermeabilityTensor::isotropic(1e-10);
        let ratio = t.anisotropy_ratio();
        assert!((ratio - 1.0).abs() < 1e-10, "ratio = {ratio}");
    }
    #[test]
    fn test_anisotropy_ratio_anisotropic() {
        let t = PermeabilityTensor::anisotropic(4e-10, 0.0, 1e-10);
        let ratio = t.anisotropy_ratio();
        assert!((ratio - 4.0).abs() < 1e-6, "ratio = {ratio}");
    }
    #[test]
    fn test_tensor_rotation_preserves_eigenvalues() {
        let t = PermeabilityTensor::anisotropic(3e-10, 0.0, 1e-10);
        let t_rot = t.rotate(std::f64::consts::PI / 4.0);
        let (k1_orig, k2_orig) = t.principal_permeabilities();
        let (k1_rot, k2_rot) = t_rot.principal_permeabilities();
        assert!(
            (k1_orig - k1_rot).abs() < 1e-22,
            "k_min changed after rotation"
        );
        assert!(
            (k2_orig - k2_rot).abs() < 1e-22,
            "k_max changed after rotation"
        );
    }
    #[test]
    fn test_rev_average_porosity_uniform() {
        let mut grid = PorousLbmGrid::new(10, 10);
        grid.set_porosity_field(|_, _| 0.5);
        let avg = RevAveraging::average_porosity(&grid, 2, 8, 2, 8);
        assert!((avg - 0.5).abs() < 1e-12, "avg porosity = {avg}");
    }
    #[test]
    fn test_rev_average_velocity_zero() {
        let grid = PorousLbmGrid::new(10, 10);
        let avg = RevAveraging::average_velocity(&grid, 0, 10, 0, 10);
        assert!(avg[0].abs() < 1e-15 && avg[1].abs() < 1e-15);
    }
    #[test]
    fn test_rev_average_density() {
        let grid = PorousLbmGrid::new(5, 5);
        let avg = RevAveraging::average_density(&grid, 0, 5, 0, 5);
        assert!((avg - 1.0).abs() < 1e-12, "avg density = {avg}");
    }
    #[test]
    fn test_interphase_heat_rate() {
        let ht = PorousHeatTransfer::new(0.6, 2.0, 1e4, 4.18e6, 2.0e6);
        let q = ht.interphase_heat_rate(300.0, 350.0);
        assert!((q - 1e4 * 50.0).abs() < 1e-6, "q = {q}");
    }
    #[test]
    fn test_interphase_heat_rate_zero_when_equal() {
        let ht = PorousHeatTransfer::new(0.6, 2.0, 1e4, 4.18e6, 2.0e6);
        let q = ht.interphase_heat_rate(300.0, 300.0);
        assert!(q.abs() < 1e-10, "q = {q}");
    }
    #[test]
    fn test_fluid_temp_rate_heating() {
        let ht = PorousHeatTransfer::new(0.6, 2.0, 1e4, 4.18e6, 2.0e6);
        let rate = ht.fluid_temp_rate(0.4, 0.0, 300.0, 400.0);
        assert!(
            rate > 0.0,
            "Fluid should heat up when T_s > T_f: rate={rate}"
        );
    }
    #[test]
    fn test_solid_temp_rate_cooling() {
        let ht = PorousHeatTransfer::new(0.6, 2.0, 1e4, 4.18e6, 2.0e6);
        let rate = ht.solid_temp_rate(0.4, 0.0, 300.0, 400.0);
        assert!(rate < 0.0, "Solid should cool when T_s > T_f: rate={rate}");
    }
    #[test]
    fn test_nusselt_wakao_zero_re() {
        let nu = PorousHeatTransfer::nusselt_wakao_kaguei(0.0, 0.7);
        assert!((nu - 2.0).abs() < 1e-10, "Nu at Re=0 should be 2: {nu}");
    }
    #[test]
    fn test_nusselt_increases_with_re() {
        let nu1 = PorousHeatTransfer::nusselt_wakao_kaguei(10.0, 0.7);
        let nu2 = PorousHeatTransfer::nusselt_wakao_kaguei(100.0, 0.7);
        assert!(
            nu2 > nu1,
            "Nu should increase with Re: nu1={nu1}, nu2={nu2}"
        );
    }
    #[test]
    fn test_darcy_velocity_field_uniform_pressure_zero() {
        let nx = 5;
        let ny = 5;
        let pressure = vec![100.0; nx * ny];
        let perm = vec![1e-10; nx * ny];
        let vel = DarcyFlow::darcy_velocity_field(&pressure, &perm, 1e-3, nx, ny, 1.0);
        for v in &vel {
            assert!(
                v[0].abs() < 1e-20 && v[1].abs() < 1e-20,
                "Uniform pressure: v = [{}, {}]",
                v[0],
                v[1]
            );
        }
    }
    #[test]
    fn test_darcy_velocity_field_linear_pressure() {
        let nx = 10;
        let ny = 5;
        let dx = 1.0;
        let pressure: Vec<f64> = (0..nx * ny)
            .map(|k| {
                let x = k % nx;
                100.0 - 10.0 * x as f64
            })
            .collect();
        let perm = vec![1e-10; nx * ny];
        let vel = DarcyFlow::darcy_velocity_field(&pressure, &perm, 1e-3, nx, ny, dx);
        let mid = ny / 2 * nx + nx / 2;
        assert!(
            vel[mid][0] > 0.0,
            "Flow should be in +x direction: {}",
            vel[mid][0]
        );
    }
    #[test]
    fn test_velocity_magnitude() {
        let mag = DarcyFlow::velocity_magnitude([3.0, 4.0]);
        assert!((mag - 5.0).abs() < 1e-14);
    }
    #[test]
    fn test_porous_reynolds_number() {
        let re = DarcyFlow::porous_reynolds_number(0.01, 1e-3, 1000.0, 1e-3);
        assert!((re - 10.0).abs() < 1e-10, "Re = {re}");
    }
    #[test]
    fn test_anisotropic_porous_cell_creation() {
        let cell = AnisotropicPorousCell::new(0.4, 1e-10);
        assert!((cell.porosity - 0.4).abs() < 1e-15);
        assert!((cell.t_fluid - 300.0).abs() < 1e-10);
        assert!((cell.t_solid - 300.0).abs() < 1e-10);
        let ratio = cell.permeability.anisotropy_ratio();
        assert!((ratio - 1.0).abs() < 1e-10, "Should be isotropic: {ratio}");
    }
    #[test]
    fn test_brinkman_force_opposes_velocity() {
        let bf = BrinkmanForce::new(1e-10, 1e-3);
        let u = [0.01, 0.0, 0.0];
        let f = bf.body_force(u);
        assert!(
            f[0] < 0.0,
            "Brinkman force should oppose velocity: f[0] = {}",
            f[0]
        );
        assert!(f[1].abs() < 1e-30);
        assert!(f[2].abs() < 1e-30);
    }
    #[test]
    fn test_brinkman_force_scales_linearly() {
        let bf = BrinkmanForce::new(1e-10, 1e-3);
        let f1 = bf.body_force([0.01, 0.0, 0.0]);
        let f2 = bf.body_force([0.02, 0.0, 0.0]);
        assert!(
            (f2[0] - 2.0 * f1[0]).abs() < 1e-20,
            "Force should scale linearly"
        );
    }
    #[test]
    fn test_darcy_pressure_drop_linear_in_length() {
        let dr = DarcyResistance::new(1e-10, 1e-3);
        let dp1 = dr.pressure_drop(1.0, 0.01);
        let dp2 = dr.pressure_drop(2.0, 0.01);
        assert!(
            (dp2 - 2.0 * dp1).abs() < 1e-15,
            "Pressure drop should be linear in L"
        );
    }
    #[test]
    fn test_darcy_pressure_drop_positive() {
        let dr = DarcyResistance::new(1e-10, 1e-3);
        let dp = dr.pressure_drop(0.1, 0.01);
        assert!(dp > 0.0, "Pressure drop should be positive: {dp}");
    }
    #[test]
    fn test_kozeny_carman_increases_with_porosity() {
        let k1 = kozeny_carman_permeability(1e-3, 0.3);
        let k2 = kozeny_carman_permeability(1e-3, 0.5);
        let k3 = kozeny_carman_permeability(1e-3, 0.7);
        assert!(k1 < k2 && k2 < k3, "k(0.3)={k1}, k(0.5)={k2}, k(0.7)={k3}");
    }
    #[test]
    fn test_ergun_pressure_drop_positive() {
        let dp = ergun_pressure_drop(0.1, 0.01, 1e-3, 0.4, 1e-3, 1000.0);
        assert!(dp > 0.0, "Ergun pressure drop should be positive: {dp}");
    }
    #[test]
    fn test_ergun_pressure_drop_increases_with_velocity() {
        let dp1 = ergun_pressure_drop(0.1, 0.01, 1e-3, 0.4, 1e-3, 1000.0);
        let dp2 = ergun_pressure_drop(0.1, 0.05, 1e-3, 0.4, 1e-3, 1000.0);
        assert!(
            dp2 > dp1,
            "Higher velocity should give higher pressure drop"
        );
    }
    #[test]
    fn test_porous_medium_type_enum() {
        let t = PorousMediumType::Brinkman;
        assert_eq!(t, PorousMediumType::Brinkman);
        assert_ne!(t, PorousMediumType::Darcy);
    }
    #[test]
    fn test_partial_saturation_force_direction() {
        let f = partial_saturation_force(0.8, 1000.0, 9.81);
        assert!(
            f[2] < 0.0,
            "Gravity force should be downward: f[2] = {}",
            f[2]
        );
        assert!(f[0].abs() < 1e-30);
        assert!(f[1].abs() < 1e-30);
    }
}
#[cfg(test)]
mod tests_porous_media_extended {
    use super::*;
    use crate::porous_media::types::*;
    #[test]
    fn test_porous_cell_new_solid_fraction() {
        let cell = PorousCell::new(0.4, 1e-10);
        assert!(
            (cell.solid_fraction - 0.6).abs() < 1e-15,
            "solid_fraction = {}",
            cell.solid_fraction
        );
    }
    #[test]
    fn test_porous_cell_default_fluid_velocity_zero() {
        let cell = PorousCell::new(0.4, 1e-10);
        assert_eq!(cell.fluid_velocity, [0.0; 2]);
    }
    #[test]
    fn test_porous_cell_default_forchheimer_zero() {
        let cell = PorousCell::new(0.4, 1e-10);
        assert!(
            cell.forchheimer_coeff.abs() < 1e-30,
            "Default Forchheimer coeff = {}",
            cell.forchheimer_coeff
        );
    }
    #[test]
    fn test_darcy_flow_velocity_zero_at_zero_gradient() {
        let u = DarcyFlow::darcy_velocity(1e-10, 1e-3, [0.0, 0.0]);
        assert!(
            u[0].abs() < 1e-30 && u[1].abs() < 1e-30,
            "Zero gradient → zero velocity"
        );
    }
    #[test]
    fn test_darcy_flow_velocity_direction_opposes_gradient() {
        let u = DarcyFlow::darcy_velocity(1e-10, 1e-3, [1000.0, 0.0]);
        assert!(
            u[0] < 0.0,
            "Positive gradient → negative velocity: {}",
            u[0]
        );
    }
    #[test]
    fn test_darcy_flow_velocity_proportional_to_permeability() {
        let u1 = DarcyFlow::darcy_velocity(1e-10, 1e-3, [-1000.0, 0.0]);
        let u2 = DarcyFlow::darcy_velocity(2e-10, 1e-3, [-1000.0, 0.0]);
        assert!(
            (u2[0] / u1[0] - 2.0).abs() < 1e-10,
            "u ∝ K: {}, {}",
            u1[0],
            u2[0]
        );
    }
    #[test]
    fn test_darcy_resistance_force_opposes_velocity() {
        let f = DarcyFlow::darcy_resistance_force(0.4, 1e-10, 1e-3, [0.01, 0.0]);
        assert!(
            f[0] < 0.0,
            "Resistance force should oppose velocity: {}",
            f[0]
        );
    }
    #[test]
    fn test_forchheimer_correction_zero_at_rest() {
        let f = DarcyFlow::forchheimer_correction([0.0, 0.0], 1e-10, 0.1, 1000.0);
        assert!(f[0].abs() < 1e-30 && f[1].abs() < 1e-30, "No force at rest");
    }
    #[test]
    fn test_forchheimer_correction_opposes_velocity() {
        let f = DarcyFlow::forchheimer_correction([0.01, 0.0], 1e-10, 0.1, 1000.0);
        assert!(
            f[0] < 0.0,
            "Forchheimer correction should oppose flow: {}",
            f[0]
        );
    }
    #[test]
    fn test_kozeny_carman_permeability_increases_with_porosity() {
        let model = KozenyCarmanModel::new(1e-3);
        let k1 = model.permeability(0.3);
        let k2 = model.permeability(0.6);
        assert!(k2 > k1, "K increases with phi: {k1}, {k2}");
    }
    #[test]
    fn test_kozeny_carman_specific_surface_area_positive() {
        let model = KozenyCarmanModel::new(1e-3);
        let s = model.specific_surface_area(0.4);
        assert!(s > 0.0, "Specific surface area = {s}");
    }
    #[test]
    fn test_kozeny_carman_specific_surface_area_zero_porosity() {
        let model = KozenyCarmanModel::new(1e-3);
        let s = model.specific_surface_area(1.0);
        assert!(s.abs() < 1e-30, "S at phi=1 should be 0: {s}");
    }
    #[test]
    fn test_kozeny_constant_value() {
        let model = KozenyCarmanModel::new(1e-3);
        assert!(
            (model.kozeny_constant() - 5.0).abs() < 1e-14,
            "Kozeny constant = {}",
            model.kozeny_constant()
        );
    }
    #[test]
    fn test_brinkman_force_scales_linearly_with_velocity() {
        let brink = BrinkmanExtension::new(1.5);
        let f1 = brink.brinkman_force([0.01, 0.0], 1e-10, 1e-3);
        let f2 = brink.brinkman_force([0.02, 0.0], 1e-10, 1e-3);
        assert!(
            (f2[0] / f1[0] - 2.0).abs() < 1e-10,
            "Brinkman force ∝ u: {}, {}",
            f1[0],
            f2[0]
        );
    }
    #[test]
    fn test_brinkman_force_zero_at_rest() {
        let brink = BrinkmanExtension::new(1.5);
        let f = brink.brinkman_force([0.0, 0.0], 1e-10, 1e-3);
        assert!(
            f[0].abs() < 1e-30 && f[1].abs() < 1e-30,
            "Zero force at rest"
        );
    }
    #[test]
    fn test_effective_thermal_conductivity_bounds() {
        let k_f: f64 = 0.6;
        let k_s: f64 = 2.5;
        let k_eff = EffectiveMediumProperties::effective_thermal_conductivity(k_f, k_s, 0.4);
        assert!(
            k_eff >= k_f.min(k_s) && k_eff <= k_f.max(k_s),
            "k_eff out of bounds: {k_eff}"
        );
    }
    #[test]
    fn test_effective_diffusivity_less_than_free() {
        let d_free = 1e-9;
        let d_eff = EffectiveMediumProperties::effective_diffusivity(d_free, 0.4, 1.5);
        assert!(
            d_eff < d_free,
            "Effective diffusivity should be less than free: {d_eff}"
        );
    }
    #[test]
    fn test_tortuosity_estimate_greater_than_one() {
        let tau = EffectiveMediumProperties::tortuosity_estimate(0.4);
        assert!(tau >= 1.0, "Tortuosity should be >= 1: {tau}");
    }
    #[test]
    fn test_porous_lbm_grid_uniform_porosity_one() {
        let mut grid = PorousLbmGrid::new(4, 4);
        grid.set_porosity_field(|_, _| 1.0);
        grid.compute_macroscopic();
        for j in 0..4 {
            for i in 0..4 {
                let idx = grid.idx(i, j);
                assert!(
                    (grid.cells[idx].density - 1.0).abs() < 1e-12,
                    "Density at ({i},{j}) = {}",
                    grid.cells[idx].density
                );
            }
        }
    }
    #[test]
    fn test_porous_lbm_equilibrium_sums_to_rho() {
        let rho = 1.1;
        let u = [0.0, 0.0];
        let sum: f64 = (0..9).map(|i| PorousLbmGrid::equilibrium(rho, u, i)).sum();
        assert!(
            (sum - rho).abs() < 1e-12,
            "Equilibrium sum = {sum}, expected {rho}"
        );
    }
    #[test]
    fn test_permeability_tensor_principal_permeabilities_isotropic() {
        let pt = PermeabilityTensor::isotropic(1e-10);
        let (k1, k2) = pt.principal_permeabilities();
        assert!((k1 - 1e-10).abs() < 1e-25, "k1 = {k1}");
        assert!((k2 - 1e-10).abs() < 1e-25, "k2 = {k2}");
    }
    #[test]
    fn test_permeability_tensor_anisotropy_ratio_isotropic_is_one() {
        let pt = PermeabilityTensor::isotropic(1e-10);
        let ratio = pt.anisotropy_ratio();
        assert!((ratio - 1.0).abs() < 1e-10, "Isotropic ratio = {ratio}");
    }
    #[test]
    fn test_rev_averaging_average_porosity() {
        let mut grid = PorousLbmGrid::new(3, 3);
        grid.set_porosity_field(|_, _| 0.3);
        let avg = RevAveraging::average_porosity(&grid, 0, 3, 0, 3);
        assert!((avg - 0.3).abs() < 1e-14, "Average porosity = {avg}");
    }
}
/// Computes the effective BGK relaxation frequency for a Brinkman porous-medium
/// cell.  The Brinkman extension modifies the viscosity:
///
/// `nu_eff = nu_fluid * effective_viscosity_ratio`
///
/// which shifts the relaxation time:
///
/// `tau_eff = nu_eff / cs^2 + 0.5`
///
/// # Arguments
/// * `nu_fluid`                - kinematic fluid viscosity (lattice units)
/// * `effective_viscosity_ratio` - `mu_eff / mu_fluid` (≥ 1 for most models)
pub fn brinkman_effective_omega(nu_fluid: f64, effective_viscosity_ratio: f64) -> f64 {
    let nu_eff = nu_fluid * effective_viscosity_ratio;
    let tau_eff = nu_eff / CS2 + 0.5;
    1.0 / tau_eff
}
/// Blake-Kozeny coefficient for the viscous term of the Ergun equation (= 150).
pub const BLAKE_KOZENY_COEFF: f64 = 150.0;
/// Burke-Plummer coefficient for the inertial term of the Ergun equation (= 1.75).
pub const BURKE_PLUMMER_COEFF: f64 = 1.75;
/// Compute the Ergun drag coefficient alpha (1/m^2) and beta (kg/m^3/m) such
/// that the pressure-drop per unit length is:
///
/// `dP/dz = alpha * mu * u + beta * rho * u^2`
///
/// # Arguments
/// * `d_p`     - particle diameter (m)
/// * `epsilon` - porosity
///
/// Returns `(alpha, beta)`.
pub fn ergun_drag_coefficients(d_p: f64, epsilon: f64) -> (f64, f64) {
    let one_minus_eps = 1.0 - epsilon;
    let eps3 = epsilon.powi(3);
    let alpha = BLAKE_KOZENY_COEFF * one_minus_eps * one_minus_eps / (d_p * d_p * eps3);
    let beta = BURKE_PLUMMER_COEFF * one_minus_eps / (d_p * eps3);
    (alpha, beta)
}
/// Pressure drop across a packed bed using the modified Ergun equation
/// with a user-supplied Kozeny-Carman constant `k_kc` (default is 5).
///
/// `ΔP/L = k_kc * 36 * mu * (1-eps)^2 / (dp^2 * eps^3) * u`
///       ` + C_F * rho * (1-eps) / (dp * eps^3) * u^2`
///
/// where `C_F = 0.55` (Forchheimer constant for sphere packings).
///
/// # Arguments
/// * `l`    - bed length (m)
/// * `u`    - superficial velocity (m/s)
/// * `d_p`  - particle diameter (m)
/// * `eps`  - porosity
/// * `mu`   - dynamic viscosity (Pa·s)
/// * `rho`  - fluid density (kg/m^3)
/// * `k_kc` - Kozeny-Carman constant (default 5)
pub fn packed_bed_pressure_drop(
    l: f64,
    u: f64,
    d_p: f64,
    eps: f64,
    mu: f64,
    rho: f64,
    k_kc: f64,
) -> f64 {
    let one_minus_eps = 1.0 - eps;
    let eps3 = eps.powi(3);
    let viscous = k_kc * 36.0 * mu * one_minus_eps * one_minus_eps / (d_p * d_p * eps3) * u;
    let inertial = 0.55 * rho * one_minus_eps / (d_p * eps3) * u * u;
    (viscous + inertial) * l
}
/// Compute the Darcy number: `Da = K / L^2`.
///
/// The Darcy number measures the relative permeability of a porous medium.
/// `Da << 1` indicates Darcy-dominated flow; `Da ~ 1` indicates Brinkman effects.
pub fn darcy_number(permeability: f64, char_length: f64) -> f64 {
    permeability / (char_length * char_length)
}
/// Compute the Forchheimer number: `Fo = C_F * Re_p / sqrt(Da)`.
///
/// The Forchheimer number quantifies the importance of inertial effects.
pub fn forchheimer_number(cf: f64, re_p: f64, da: f64) -> f64 {
    if da < 1e-300 {
        return f64::MAX;
    }
    cf * re_p / da.sqrt()
}
/// Porous Reynolds number using superficial velocity:
///
/// `Re = rho * u * d_p / (mu * (1 - eps))`
pub fn superficial_reynolds_number(rho: f64, u: f64, d_p: f64, mu: f64, eps: f64) -> f64 {
    rho * u * d_p / (mu * (1.0 - eps).max(1e-12))
}
/// Compute the effective lattice permeability of a `PorousLbmGrid` by running
/// a short Poiseuille-like simulation with a known pressure gradient and
/// fitting Darcy's law.
///
/// `K_eff = nu * u_mean / |grad_p|`
///
/// The grid is driven by a body force `f_x` in the x-direction for `n_steps`
/// and the volume-averaged x-velocity is used.
///
/// Returns the estimated permeability in lattice units.
pub fn estimate_lattice_permeability(
    grid: &mut PorousLbmGrid,
    tau: f64,
    f_x: f64,
    n_steps: usize,
) -> f64 {
    let nu = CS2 * (tau - 0.5);
    for _ in 0..n_steps {
        grid.compute_macroscopic();
        let nx = grid.nx;
        let ny = grid.ny;
        for y in 0..ny {
            for x in 0..nx {
                let k = grid.idx(x, y);
                let rho = grid.cells[k].density;
                let omega = 1.0 / tau;
                let u = grid.cells[k].velocity;
                let sigma = grid.cells[k].resistance;
                let denom = 1.0 + sigma * tau;
                let u_eff = [u[0] / denom, u[1] / denom];
                for i in 0..9 {
                    let feq = PorousLbmGrid::equilibrium(rho, u_eff, i);
                    grid.cells[k].f[i] += omega * (feq - grid.cells[k].f[i]);
                    let cu = CX[i] * u_eff[0] + CY[i] * u_eff[1];
                    let fi_term =
                        W[i] * (1.0 - 0.5 * omega) * (CX[i] / CS2 + cu * CX[i] / (CS2 * CS2)) * f_x;
                    grid.cells[k].f[i] += fi_term;
                }
            }
        }
        grid.stream();
        grid.compute_macroscopic();
    }
    let n = grid.nx * grid.ny;
    let u_mean: f64 = (0..n).map(|k| grid.cells[k].velocity[0]).sum::<f64>() / n as f64;
    if f_x.abs() < 1e-300 {
        return 0.0;
    }
    nu * u_mean.abs() / f_x.abs()
}
#[cfg(test)]
mod tests_porous_extra {
    use super::*;
    use crate::porous_media::types::*;
    #[test]
    fn test_ergun_drag_coefficients_positive() {
        let (alpha, beta) = ergun_drag_coefficients(1e-3, 0.4);
        assert!(alpha > 0.0, "alpha = {alpha}");
        assert!(beta > 0.0, "beta = {beta}");
    }
    #[test]
    fn test_ergun_drag_alpha_increases_with_decreasing_porosity() {
        let (a1, _) = ergun_drag_coefficients(1e-3, 0.6);
        let (a2, _) = ergun_drag_coefficients(1e-3, 0.3);
        assert!(
            a2 > a1,
            "More compact bed should have higher viscous coeff: {a1} vs {a2}"
        );
    }
    #[test]
    fn test_packed_bed_pressure_drop_positive() {
        let dp = packed_bed_pressure_drop(0.1, 0.01, 1e-3, 0.4, 1e-3, 1000.0, 5.0);
        assert!(dp > 0.0, "packed bed pressure drop = {dp}");
    }
    #[test]
    fn test_packed_bed_pressure_drop_vs_ergun_order_of_magnitude() {
        let ergun = ergun_pressure_drop(0.1, 0.01, 1e-3, 0.4, 1e-3, 1000.0);
        let packed = packed_bed_pressure_drop(0.1, 0.01, 1e-3, 0.4, 1e-3, 1000.0, 5.0);
        assert!(ergun > 0.0 && packed > 0.0);
        let ratio = ergun / packed;
        assert!(
            ratio > 0.01 && ratio < 100.0,
            "ratio Ergun/packed = {ratio}"
        );
    }
    #[test]
    fn test_darcy_number_positive() {
        let da = darcy_number(1e-10, 0.01);
        assert!(da > 0.0, "Da = {da}");
    }
    #[test]
    fn test_darcy_number_scales_with_permeability() {
        let da1 = darcy_number(1e-10, 0.01);
        let da2 = darcy_number(2e-10, 0.01);
        assert!((da2 / da1 - 2.0).abs() < 1e-10, "Da ∝ K: {da1}, {da2}");
    }
    #[test]
    fn test_forchheimer_number_positive() {
        let fo = forchheimer_number(0.55, 10.0, 1e-4);
        assert!(fo > 0.0 && fo.is_finite(), "Fo = {fo}");
    }
    #[test]
    fn test_superficial_reynolds_number() {
        let re = superficial_reynolds_number(1000.0, 0.01, 1e-3, 1e-3, 0.4);
        assert!(re > 0.0 && re.is_finite(), "Re = {re}");
    }
    #[test]
    fn test_kozeny_carman_extended_permeability_positive() {
        let model = KozenyCarmanExtended::new(1e-3, 5.0, 1.5);
        let k = model.permeability(0.4);
        assert!(k > 0.0 && k.is_finite(), "K = {k}");
    }
    #[test]
    fn test_kozeny_carman_extended_increases_with_porosity() {
        let model = KozenyCarmanExtended::new(1e-3, 5.0, 1.5);
        let k1 = model.permeability(0.3);
        let k2 = model.permeability(0.5);
        assert!(k2 > k1, "K should increase with porosity: {k1} < {k2}");
    }
    #[test]
    fn test_kozeny_carman_extended_hydraulic_diameter_positive() {
        let model = KozenyCarmanExtended::new(1e-3, 5.0, 1.5);
        let dh = model.hydraulic_diameter(0.4);
        assert!(dh > 0.0, "Hydraulic diameter = {dh}");
    }
    #[test]
    fn test_kozeny_carman_extended_vs_standard() {
        let standard = KozenyCarmanModel::new(1e-3);
        let extended = KozenyCarmanExtended::new(1e-3, 5.0, 1.0);
        let k_std = standard.permeability(0.4);
        let k_ext = extended.permeability(0.4);
        let ratio = k_std / k_ext;
        assert!(ratio > 0.01 && ratio < 100.0, "ratio std/ext = {ratio}");
    }
    #[test]
    fn test_dbf_darcy_force_opposes_flow() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-10, 0.4, 1e-3, 0.55, 1000.0);
        let f = dbf.darcy_force([0.01, 0.0]);
        assert!(f[0] < 0.0, "Darcy force should oppose +x flow: {}", f[0]);
    }
    #[test]
    fn test_dbf_forchheimer_force_opposes_flow() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-10, 0.4, 1e-3, 0.55, 1000.0);
        let f = dbf.forchheimer_force([0.01, 0.0]);
        assert!(
            f[0] < 0.0,
            "Forchheimer force should oppose +x flow: {}",
            f[0]
        );
    }
    #[test]
    fn test_dbf_total_force_at_rest_zero() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-10, 0.4, 1e-3, 0.55, 1000.0);
        let f = dbf.total_force([0.0, 0.0]);
        assert!(f[0].abs() < 1e-30 && f[1].abs() < 1e-30, "No force at rest");
    }
    #[test]
    fn test_dbf_brinkman_omega_positive() {
        let dbf = DarcyBrinkmanForchheimer::new(1e-10, 0.4, 1e-3, 0.55, 1000.0);
        let omega = dbf.brinkman_omega();
        assert!(omega > 0.0 && omega < 2.0, "omega = {omega}");
    }
    #[test]
    fn test_brooks_corey_capillary_at_residual_is_max() {
        let bc = BrooksCoreyCapillary::new(1000.0, 2.0, 0.1);
        let pc = bc.capillary_pressure(0.1);
        assert!(
            pc == f64::MAX || pc > 1e10,
            "Pc at residual should be very large: {pc}"
        );
    }
    #[test]
    fn test_brooks_corey_capillary_decreases_with_saturation() {
        let bc = BrooksCoreyCapillary::new(1000.0, 2.0, 0.1);
        let pc1 = bc.capillary_pressure(0.3);
        let pc2 = bc.capillary_pressure(0.7);
        assert!(
            pc2 < pc1,
            "Pc should decrease with saturation: {pc1} vs {pc2}"
        );
    }
    #[test]
    fn test_brooks_corey_water_relperm_at_full_saturation() {
        let bc = BrooksCoreyCapillary::new(1000.0, 2.0, 0.1);
        let krw = bc.water_relative_permeability(1.0);
        assert!((krw - 1.0).abs() < 1e-10, "k_rw at full saturation = {krw}");
    }
    #[test]
    fn test_brooks_corey_water_relperm_increases_with_saturation() {
        let bc = BrooksCoreyCapillary::new(1000.0, 2.0, 0.1);
        let krw1 = bc.water_relative_permeability(0.3);
        let krw2 = bc.water_relative_permeability(0.7);
        assert!(
            krw2 > krw1,
            "k_rw should increase with S_w: {krw1} vs {krw2}"
        );
    }
    #[test]
    fn test_brooks_corey_nonwetting_relperm_at_residual_saturation() {
        let bc = BrooksCoreyCapillary::new(1000.0, 2.0, 0.1);
        let krnw = bc.nonwetting_relative_permeability(0.1);
        assert!((krnw - 1.0).abs() < 1e-10, "k_rnw at residual = {krnw}");
    }
    #[test]
    fn test_porous_media_driver_creates_grid() {
        let drv = PorousMediaDriver::new(8, 6, 1.0 / 6.0, 0.4, 1e-10);
        assert_eq!(drv.grid.cells.len(), 48, "Grid should have 48 cells");
    }
    #[test]
    fn test_porous_media_driver_step_increments_count() {
        let mut drv = PorousMediaDriver::new(6, 6, 1.0 / 6.0, 0.4, 1e-3);
        drv.run(5);
        assert_eq!(drv.step_count, 5, "step_count should be 5");
    }
    #[test]
    fn test_porous_media_driver_body_force_generates_flow() {
        let mut drv = PorousMediaDriver::new(8, 8, 1.0 / 6.0, 0.6, 1e-3);
        drv.set_body_force(1e-4, 0.0);
        drv.run(100);
        let u_mean = drv.mean_superficial_velocity();
        assert!(
            u_mean[0] > 0.0,
            "Body force should drive flow in +x: {}",
            u_mean[0]
        );
    }
    #[test]
    fn test_porous_media_driver_mean_velocity_finite() {
        let mut drv = PorousMediaDriver::new(4, 4, 1.0 / 6.0, 0.4, 1e-3);
        drv.run(10);
        let u = drv.mean_superficial_velocity();
        assert!(
            u[0].is_finite() && u[1].is_finite(),
            "Velocity must be finite"
        );
    }
    #[test]
    fn test_brinkman_effective_omega_ratio_one() {
        let nu = 1.0 / 6.0;
        let omega_standard = 1.0 / (nu / CS2 + 0.5);
        let omega_brinkman = brinkman_effective_omega(nu, 1.0);
        assert!(
            (omega_brinkman - omega_standard).abs() < 1e-12,
            "omega mismatch: {omega_brinkman} vs {omega_standard}"
        );
    }
    #[test]
    fn test_brinkman_effective_omega_larger_ratio_smaller_omega() {
        let nu = 1.0 / 6.0;
        let omega1 = brinkman_effective_omega(nu, 1.0);
        let omega2 = brinkman_effective_omega(nu, 2.0);
        assert!(
            omega2 < omega1,
            "Larger mu_eff ratio → smaller omega: {omega1}, {omega2}"
        );
    }
    #[test]
    fn test_estimate_lattice_permeability_positive() {
        let mut grid = PorousLbmGrid::new(6, 6);
        let k = estimate_lattice_permeability(&mut grid, 1.0, 1e-5, 50);
        assert!(k >= 0.0, "Permeability must be non-negative: {k}");
    }
    #[test]
    fn test_estimate_lattice_permeability_larger_force_same_order() {
        let mut grid1 = PorousLbmGrid::new(4, 4);
        let mut grid2 = PorousLbmGrid::new(4, 4);
        let k1 = estimate_lattice_permeability(&mut grid1, 1.0, 1e-5, 30);
        let k2 = estimate_lattice_permeability(&mut grid2, 1.0, 2e-5, 30);
        if k1 > 1e-20 && k2 > 1e-20 {
            let ratio = k1 / k2;
            assert!(
                ratio > 0.1 && ratio < 10.0,
                "Permeability estimates diverge: {k1}, {k2}"
            );
        }
    }
}
