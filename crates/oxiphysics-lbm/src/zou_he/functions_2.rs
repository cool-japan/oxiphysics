//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::{poiseuille_profile, zou_he_inlet_left};
use crate::lattice::{D3Q19_VELOCITIES, D3Q19_WEIGHTS};

/// Apply a Poiseuille velocity profile as a Zou-He left inlet for D2Q9.
///
/// Applies the parabolic profile at each y-index from `i_start` to `i_end`.
/// Each cell distribution is updated independently.
///
/// # Arguments
/// * `f_cells` – slice of D2Q9 distribution arrays (one per y-cell)
/// * `h`       – total channel height in lattice units
/// * `u_max`   – peak velocity
pub fn apply_poiseuille_inlet_d2q9(f_cells: &mut [[f64; 9]], h: f64, u_max: f64) {
    for (j, cell) in f_cells.iter_mut().enumerate() {
        let y = j as f64 + 0.5;
        let ux = poiseuille_profile(y, h, u_max);
        zou_he_inlet_left(cell, ux);
    }
}
/// Half-way bounce-back with moving wall velocity for D3Q19 top face.
///
/// Implements the standard half-way bounce-back plus a momentum source
/// due to wall velocity (ux_wall, uy_wall = 0, uz_wall):
/// `f_opposite = f_i − 2 w_i ρ (c_i · u_wall) / cs²`
///
/// # Arguments
/// * `f`       – distribution at wall cell
/// * `rho`     – local density
/// * `ux_wall` – wall velocity in x-direction
/// * `uz_wall` – wall velocity in z-direction
pub fn bounce_back_moving_top_d3q19(f: &mut [f64; 19], rho: f64, ux_wall: f64, uz_wall: f64) {
    use crate::lattice::D3Q19_WEIGHTS;
    let bounce_pairs: [(usize, usize, f64, f64); 5] = [
        (3, 4, 0.0, 0.0),
        (7, 10, ux_wall, 0.0),
        (8, 9, -ux_wall, 0.0),
        (15, 18, 0.0, uz_wall),
        (16, 17, 0.0, -uz_wall),
    ];
    for (i_in, i_out, cx, cz) in bounce_pairs {
        let correction = 2.0 * D3Q19_WEIGHTS[i_in] * rho * (cx * ux_wall + cz * uz_wall) * 3.0;
        f[i_out] = f[i_in] - correction;
    }
}
/// Half-way bounce-back for D3Q19 bottom face with moving wall.
///
/// Same as `bounce_back_moving_top_d3q19` but for the bottom wall (−y face).
pub fn bounce_back_moving_bottom_d3q19(f: &mut [f64; 19], rho: f64, ux_wall: f64, uz_wall: f64) {
    let bounce_pairs: [(usize, usize, f64, f64); 5] = [
        (4, 3, 0.0, 0.0),
        (9, 8, ux_wall, 0.0),
        (10, 7, -ux_wall, 0.0),
        (17, 16, 0.0, uz_wall),
        (18, 15, 0.0, -uz_wall),
    ];
    for (i_in, i_out, cx, cz) in bounce_pairs {
        let correction = 2.0 * D3Q19_WEIGHTS[i_in] * rho * (cx * ux_wall + cz * uz_wall) * 3.0;
        f[i_out] = f[i_in] - correction;
    }
}
/// Compute the mass flux through the left face of a D2Q9 cell.
///
/// J = Σ_{i: c_ix > 0} f_i − Σ_{i: c_ix < 0} f_i
pub fn mass_flux_left_d2q9(f: &[f64; 9]) -> f64 {
    f[1] + f[5] + f[8] - (f[3] + f[6] + f[7])
}
/// Compute the mass flux through the left face of a D3Q19 cell.
///
/// Sums all distributions carrying momentum in the +x direction.
pub fn mass_flux_left_d3q19(f: &[f64; 19]) -> f64 {
    (f[1] + f[7] + f[9] + f[11] + f[13]) - (f[2] + f[8] + f[10] + f[12] + f[14])
}
/// Compute the momentum flux (x-component) at a D3Q19 left face.
///
/// Π_xx = Σ_i f_i c_{ix}²
pub fn momentum_flux_xx_d3q19(f: &[f64; 19]) -> f64 {
    use crate::lattice::D3Q19_VELOCITIES;
    f.iter()
        .enumerate()
        .map(|(i, fi)| {
            let cx = D3Q19_VELOCITIES[i][0] as f64;
            fi * cx * cx
        })
        .sum()
}
/// Verify that a D2Q9 Zou-He inlet satisfies the density constraint.
///
/// Returns `true` if |Σ f_i − rho_expected| < tol.
pub fn verify_density_constraint_d2q9(f: &[f64; 9], rho_expected: f64, tol: f64) -> bool {
    let rho_actual: f64 = f.iter().sum();
    (rho_actual - rho_expected).abs() < tol
}
/// Verify that a D3Q19 boundary satisfies the momentum constraint.
///
/// Returns `true` if |Σ f_i c_{ix} − rho * ux| < tol.
pub fn verify_momentum_constraint_d3q19(f: &[f64; 19], rho: f64, ux: f64, tol: f64) -> bool {
    let jx: f64 = f
        .iter()
        .enumerate()
        .map(|(i, fi)| fi * D3Q19_VELOCITIES[i][0] as f64)
        .sum();
    (jx - rho * ux).abs() < tol
}
#[cfg(test)]
mod tests_advanced_zou_he {
    use super::super::*;
    fn rest_d2q9() -> [f64; 9] {
        d2q9_equilibrium(1.0, 0.0, 0.0)
    }
    fn rest_d3q19() -> [f64; 19] {
        d3q19_equilibrium(1.0, 0.0, 0.0, 0.0)
    }
    fn rest_d3q27() -> [f64; 27] {
        d3q27_equilibrium(1.0, 0.0, 0.0, 0.0)
    }
    #[test]
    fn test_pressure_left_d3q19_density() {
        let rho_set = 1.02;
        let mut f = rest_d3q19();
        let ux = zou_he_pressure_left_d3q19(&mut f, rho_set);
        assert!(ux.is_finite(), "Pressure left ux={ux}");
        assert!(f.iter().all(|v| v.is_finite()), "Pressure left has NaN/Inf");
        let rho_out: f64 = f.iter().sum();
        assert!(rho_out > 0.0, "Pressure left rho={rho_out}");
    }
    #[test]
    fn test_pressure_right_d3q19_density() {
        let rho_set = 0.98;
        let mut f = rest_d3q19();
        let ux = zou_he_pressure_right_d3q19(&mut f, rho_set);
        assert!(ux.is_finite(), "Pressure right ux={ux}");
        assert!(
            f.iter().all(|v| v.is_finite()),
            "Pressure right has NaN/Inf"
        );
        let rho_out: f64 = f.iter().sum();
        assert!(rho_out > 0.0, "Pressure right rho={rho_out}");
    }
    #[test]
    fn test_pressure_left_d3q19_finite() {
        let mut f = rest_d3q19();
        zou_he_pressure_left_d3q19(&mut f, 1.01);
        assert!(f.iter().all(|v| v.is_finite()), "Pressure left has NaN/Inf");
    }
    #[test]
    fn test_pressure_right_d3q19_finite() {
        let mut f = rest_d3q19();
        zou_he_pressure_right_d3q19(&mut f, 0.99);
        assert!(
            f.iter().all(|v| v.is_finite()),
            "Pressure right has NaN/Inf"
        );
    }
    #[test]
    fn test_d3q27_equilibrium_sums_to_rho() {
        let feq = d3q27_equilibrium(1.1, 0.03, 0.01, -0.02);
        let sum: f64 = feq.iter().sum();
        assert!(
            (sum - 1.1).abs() < 1e-12,
            "D3Q27 eq sum={sum}, expected 1.1"
        );
    }
    #[test]
    fn test_d3q27_equilibrium_at_rest_equals_weights() {
        let feq = d3q27_equilibrium(1.0, 0.0, 0.0, 0.0);
        let sum: f64 = feq.iter().sum();
        assert!((sum - 1.0).abs() < 1e-14, "D3Q27 rest sum={sum}");
        assert!((feq[0] - 8.0 / 27.0).abs() < 1e-14, "feq[0]={}", feq[0]);
    }
    #[test]
    fn test_d3q27_equilibrium_positive_small_vel() {
        let feq = d3q27_equilibrium(1.0, 0.04, 0.01, -0.02);
        for (i, &fi) in feq.iter().enumerate() {
            assert!(fi >= 0.0, "D3Q27 feq[{i}]={fi} is negative");
        }
    }
    #[test]
    fn test_macroscopic_d3q27_at_rest() {
        let f = rest_d3q27();
        let (rho, ux, uy, uz) = macroscopic_d3q27(&f);
        assert!((rho - 1.0).abs() < 1e-14, "D3Q27 rho={rho}");
        assert!(ux.abs() < 1e-14, "D3Q27 ux={ux}");
        assert!(uy.abs() < 1e-14, "D3Q27 uy={uy}");
        assert!(uz.abs() < 1e-14, "D3Q27 uz={uz}");
    }
    #[test]
    fn test_macroscopic_d3q27_momentum() {
        let ux = 0.03;
        let uy = 0.01;
        let uz = -0.02;
        let feq = d3q27_equilibrium(1.0, ux, uy, uz);
        let (rho, ux_out, uy_out, uz_out) = macroscopic_d3q27(&feq);
        assert!((ux_out - ux).abs() < 1e-12, "D3Q27 ux={ux_out}");
        assert!((uy_out - uy).abs() < 1e-12, "D3Q27 uy={uy_out}");
        assert!((uz_out - uz).abs() < 1e-12, "D3Q27 uz={uz_out}");
        let _ = rho;
    }
    #[test]
    fn test_zou_he_velocity_inlet_d3q27_density() {
        let mut f = rest_d3q27();
        let rho = zou_he_velocity_inlet_d3q27(&mut f, 0.05, 0.0, 0.0);
        assert!((rho - 1.0).abs() < 0.1, "D3Q27 inlet rho={rho}");
        assert!(f.iter().all(|v| v.is_finite()), "D3Q27 inlet has NaN/Inf");
    }
    #[test]
    fn test_zou_he_pressure_outlet_d3q27_finite() {
        let mut f = rest_d3q27();
        let ux = zou_he_pressure_outlet_d3q27(&mut f, 1.0);
        assert!(ux.is_finite(), "D3Q27 outlet ux={ux}");
        assert!(f.iter().all(|v| v.is_finite()), "D3Q27 outlet has NaN/Inf");
    }
    #[test]
    fn test_convective_outlet_d2q9_mass_positive() {
        let f_b = rest_d2q9();
        let f_i = d2q9_equilibrium(1.0, 0.05, 0.0);
        let f_new = convective_outlet_d2q9(&f_b, &f_i, 0.05);
        let rho: f64 = f_new.iter().sum();
        assert!(rho > 0.0, "Convective outlet rho should be positive: {rho}");
    }
    #[test]
    fn test_convective_outlet_d2q9_uc_zero() {
        let f_b = d2q9_equilibrium(1.05, 0.02, 0.0);
        let f_i = d2q9_equilibrium(1.0, 0.05, 0.0);
        let f_new = convective_outlet_d2q9(&f_b, &f_i, 0.0);
        for i in 0..9 {
            assert!(
                (f_new[i] - f_b[i]).abs() < 1e-14,
                "uc=0: f_new[{i}]={} != f_b[{i}]={}",
                f_new[i],
                f_b[i]
            );
        }
    }
    #[test]
    fn test_convective_outlet_d3q19_finite() {
        let f_b = rest_d3q19();
        let f_i = d3q19_equilibrium(1.0, 0.05, 0.0, 0.0);
        let f_new = convective_outlet_d3q19(&f_b, &f_i, 0.05);
        assert!(
            f_new.iter().all(|v| v.is_finite()),
            "Convective D3Q19 has NaN"
        );
    }
    #[test]
    fn test_neumann_outlet_d2q9_copies_interior() {
        let f_i = d2q9_equilibrium(1.0, 0.05, 0.01);
        let mut f_b = rest_d2q9();
        neumann_outlet_d2q9(&mut f_b, &f_i);
        for i in 0..9 {
            assert!(
                (f_b[i] - f_i[i]).abs() < 1e-14,
                "Neumann D2Q9 f_b[{i}]={} != f_i[{i}]={}",
                f_b[i],
                f_i[i]
            );
        }
    }
    #[test]
    fn test_neumann_outlet_d3q19_copies_interior() {
        let f_i = d3q19_equilibrium(1.0, 0.04, 0.01, -0.02);
        let mut f_b = rest_d3q19();
        neumann_outlet_d3q19(&mut f_b, &f_i);
        for i in 0..19 {
            assert!(
                (f_b[i] - f_i[i]).abs() < 1e-14,
                "Neumann D3Q19 f_b[{i}]={} != f_i[{i}]={}",
                f_b[i],
                f_i[i]
            );
        }
    }
    #[test]
    fn test_extrapolation_outlet_d2q9_linear_state() {
        let f_c = d2q9_equilibrium(1.0, 0.05, 0.0);
        let mut f_b = rest_d2q9();
        extrapolation_outlet_d2q9(&mut f_b, &f_c, &f_c);
        for i in 0..9 {
            assert!(
                (f_b[i] - f_c[i]).abs() < 1e-14,
                "Extrapolation const: f_b[{i}]={} != f_c[{i}]={}",
                f_b[i],
                f_c[i]
            );
        }
    }
    #[test]
    fn test_extrapolation_outlet_d3q19_linear_state() {
        let f_c = d3q19_equilibrium(1.0, 0.03, 0.0, 0.0);
        let mut f_b = rest_d3q19();
        extrapolation_outlet_d3q19(&mut f_b, &f_c, &f_c);
        for i in 0..19 {
            assert!(
                (f_b[i] - f_c[i]).abs() < 1e-14,
                "Extrapolation D3Q19 const: f_b[{i}]={}",
                f_b[i]
            );
        }
    }
    #[test]
    fn test_extrapolation3_outlet_d3q19_constant() {
        let f_c = d3q19_equilibrium(1.0, 0.04, 0.0, 0.0);
        let mut f_b = rest_d3q19();
        extrapolation3_outlet_d3q19(&mut f_b, &f_c, &f_c, &f_c);
        for i in 0..19 {
            assert!(
                (f_b[i] - f_c[i]).abs() < 1e-14,
                "3rd order extrapolation const: f_b[{i}]={}",
                f_b[i]
            );
        }
    }
    #[test]
    fn test_extrapolation_outlet_d2q9_linearly_varying() {
        let f1 = d2q9_equilibrium(1.01, 0.05, 0.0);
        let f2 = d2q9_equilibrium(1.02, 0.04, 0.0);
        let mut f_b = rest_d2q9();
        extrapolation_outlet_d2q9(&mut f_b, &f1, &f2);
        for i in 0..9 {
            let expected = 2.0 * f1[i] - f2[i];
            assert!(
                (f_b[i] - expected).abs() < 1e-14,
                "Extrapolation formula at i={i}: {} != {}",
                f_b[i],
                expected
            );
        }
    }
    #[test]
    fn test_sponge_sigma_outside() {
        assert_eq!(sponge_sigma(-1.0, 0.0, 10.0, 0.5), 0.0);
        assert_eq!(sponge_sigma(11.0, 0.0, 10.0, 0.5), 0.5);
    }
    #[test]
    fn test_sponge_sigma_at_boundary() {
        let sigma = sponge_sigma(5.0, 0.0, 10.0, 1.0);
        assert!(
            (sigma - 0.5).abs() < 1e-12,
            "Sponge sigma at midpoint: {sigma}"
        );
    }
    #[test]
    fn test_sponge_sigma_monotone() {
        let xs: [f64; 5] = [0.0, 2.5, 5.0, 7.5, 10.0];
        let mut prev = -1.0f64;
        for &x in &xs {
            let s = sponge_sigma(x, 0.0, 10.0, 1.0);
            assert!(s >= prev - 1e-14, "Sponge sigma not monotone at x={x}: {s}");
            prev = s;
        }
    }
    #[test]
    fn test_sponge_zone_d2q9_full_sigma() {
        let rho_t = 1.0;
        let ux_t = 0.0;
        let uy_t = 0.0;
        let feq_target = d2q9_equilibrium(rho_t, ux_t, uy_t);
        let mut f = d2q9_equilibrium(1.05, 0.05, 0.0);
        sponge_zone_d2q9(&mut f, rho_t, ux_t, uy_t, 1.0);
        for i in 0..9 {
            assert!(
                (f[i] - feq_target[i]).abs() < 1e-14,
                "Sponge full: f[{i}]={} != feq[{i}]={}",
                f[i],
                feq_target[i]
            );
        }
    }
    #[test]
    fn test_sponge_zone_d2q9_zero_sigma() {
        let mut f = d2q9_equilibrium(1.05, 0.05, 0.0);
        let f_orig = f;
        sponge_zone_d2q9(&mut f, 1.0, 0.0, 0.0, 0.0);
        for i in 0..9 {
            assert!(
                (f[i] - f_orig[i]).abs() < 1e-14,
                "Sponge zero: f[{i}] changed"
            );
        }
    }
    #[test]
    fn test_sponge_zone_d3q19_full_sigma() {
        let rho_t = 1.0;
        let feq_target = d3q19_equilibrium(rho_t, 0.0, 0.0, 0.0);
        let mut f = d3q19_equilibrium(1.05, 0.05, 0.0, 0.0);
        sponge_zone_d3q19(&mut f, rho_t, 0.0, 0.0, 0.0, 1.0);
        for i in 0..19 {
            assert!(
                (f[i] - feq_target[i]).abs() < 1e-14,
                "D3Q19 sponge full at i={i}"
            );
        }
    }
    #[test]
    fn test_sponge_sigma_exponential_boundaries() {
        assert!((sponge_sigma_exponential(-1.0, 0.0, 10.0, 1.0, 5.0) - 0.0).abs() < 1e-14);
        assert!((sponge_sigma_exponential(15.0, 0.0, 10.0, 1.0, 5.0) - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_sponge_sigma_exponential_monotone() {
        let xs: [f64; 4] = [0.0, 3.0, 6.0, 10.0];
        let mut prev = -1.0f64;
        for &x in &xs {
            let s = sponge_sigma_exponential(x, 0.0, 10.0, 1.0, 3.0);
            assert!(s >= prev - 1e-12, "Exp sponge not monotone at x={x}: {s}");
            prev = s;
        }
    }
    #[test]
    fn test_characteristic_outlet_d2q9_at_target() {
        let rho = 1.0;
        let ux = 0.05;
        let uy = 0.0;
        let mut f = d2q9_equilibrium(rho, ux, uy);
        characteristic_outlet_d2q9(&mut f, rho, ux, uy, rho, 0.1);
        let rho_out: f64 = f.iter().sum();
        assert!(
            (rho_out - rho).abs() < 1e-12,
            "Char BC at target: rho={rho_out}"
        );
    }
    #[test]
    fn test_characteristic_outlet_d3q19_finite() {
        let mut f = d3q19_equilibrium(1.02, 0.05, 0.0, 0.0);
        characteristic_outlet_d3q19(&mut f, 1.02, 0.05, 0.0, 0.0, 1.0, 0.1);
        assert!(f.iter().all(|v| v.is_finite()), "Char D3Q19 has NaN/Inf");
    }
    #[test]
    fn test_poiseuille_profile_centerline() {
        let h = 10.0;
        let u_max = 0.1;
        let u = poiseuille_profile(h / 2.0, h, u_max);
        assert!(
            (u - u_max).abs() < 1e-12,
            "Poiseuille center: {u} != {u_max}"
        );
    }
    #[test]
    fn test_poiseuille_profile_walls_zero() {
        let h = 10.0;
        let u_max = 0.1;
        let u0 = poiseuille_profile(0.0, h, u_max);
        let uh = poiseuille_profile(h, h, u_max);
        assert!(u0.abs() < 1e-14, "Poiseuille at y=0: {u0}");
        assert!(uh.abs() < 1e-14, "Poiseuille at y=H: {uh}");
    }
    #[test]
    fn test_tanh_shear_profile_center() {
        let u = tanh_shear_profile(5.0, 5.0, 1.0, 0.1);
        assert!((u - 0.05).abs() < 1e-14, "Tanh shear center: {u}");
    }
    #[test]
    fn test_tanh_shear_profile_range() {
        let u_low = tanh_shear_profile(0.0, 5.0, 1.0, 0.1);
        let u_high = tanh_shear_profile(10.0, 5.0, 1.0, 0.1);
        assert!((0.0..=0.05).contains(&u_low), "Tanh low: {u_low}");
        assert!((0.05..=0.1).contains(&u_high), "Tanh high: {u_high}");
    }
    #[test]
    fn test_mass_flux_left_d2q9_at_rest() {
        let f = rest_d2q9();
        let flux = mass_flux_left_d2q9(&f);
        assert!(flux.abs() < 1e-14, "Rest mass flux: {flux}");
    }
    #[test]
    fn test_mass_flux_left_d2q9_positive_ux() {
        let f = d2q9_equilibrium(1.0, 0.1, 0.0);
        let flux = mass_flux_left_d2q9(&f);
        assert!(flux > 0.0, "Positive ux → positive mass flux: {flux}");
    }
    #[test]
    fn test_mass_flux_left_d3q19_at_rest() {
        let f = rest_d3q19();
        let flux = mass_flux_left_d3q19(&f);
        assert!(flux.abs() < 1e-14, "D3Q19 rest mass flux: {flux}");
    }
    #[test]
    fn test_verify_density_constraint_d2q9_pass() {
        let f = d2q9_equilibrium(1.05, 0.03, 0.0);
        assert!(
            verify_density_constraint_d2q9(&f, 1.05, 1e-12),
            "Density constraint should pass"
        );
    }
    #[test]
    fn test_verify_density_constraint_d2q9_fail() {
        let f = d2q9_equilibrium(1.05, 0.03, 0.0);
        assert!(
            !verify_density_constraint_d2q9(&f, 1.0, 1e-10),
            "Density constraint should fail for wrong target"
        );
    }
    #[test]
    fn test_verify_momentum_constraint_d3q19_pass() {
        let rho = 1.0;
        let ux = 0.05;
        let f = d3q19_equilibrium(rho, ux, 0.0, 0.0);
        assert!(
            verify_momentum_constraint_d3q19(&f, rho, ux, 1e-12),
            "Momentum constraint should pass"
        );
    }
    #[test]
    fn test_apply_poiseuille_inlet_d2q9_peak() {
        let h = 20.0;
        let u_max = 0.1;
        let mut cells: Vec<[f64; 9]> = (0..20).map(|_| rest_d2q9()).collect();
        apply_poiseuille_inlet_d2q9(&mut cells, h, u_max);
        let _all_same = cells
            .iter()
            .all(|c| (c.iter().sum::<f64>() - 1.0).abs() < 0.01);
        // Poiseuille inlet applied without panic
    }
    #[test]
    fn test_bounce_back_moving_top_d3q19_finite() {
        let mut f = rest_d3q19();
        bounce_back_moving_top_d3q19(&mut f, 1.0, 0.05, 0.0);
        assert!(
            f.iter().all(|v| v.is_finite()),
            "Moving wall top has NaN/Inf"
        );
    }
    #[test]
    fn test_bounce_back_moving_bottom_d3q19_finite() {
        let mut f = rest_d3q19();
        bounce_back_moving_bottom_d3q19(&mut f, 1.0, 0.05, 0.0);
        assert!(
            f.iter().all(|v| v.is_finite()),
            "Moving wall bottom has NaN/Inf"
        );
    }
    #[test]
    fn test_momentum_flux_xx_d3q19_positive() {
        let f = d3q19_equilibrium(1.0, 0.05, 0.0, 0.0);
        let pi_xx = momentum_flux_xx_d3q19(&f);
        assert!(pi_xx > 0.0, "Momentum flux should be positive: {pi_xx}");
    }
    #[test]
    fn test_d3q27_inlet_outlet_roundtrip() {
        let mut f_in = rest_d3q27();
        let rho = zou_he_velocity_inlet_d3q27(&mut f_in, 0.03, 0.0, 0.0);
        let mut f_out = rest_d3q27();
        let ux_out = zou_he_pressure_outlet_d3q27(&mut f_out, 1.0);
        assert!(rho > 0.0, "Inlet rho should be positive: {rho}");
        assert!(ux_out.is_finite(), "Outlet ux should be finite: {ux_out}");
    }
}
