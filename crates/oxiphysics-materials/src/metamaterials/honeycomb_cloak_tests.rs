//! Regression tests for the re-entrant honeycomb flexure model
//! ([`ReentrantHoneycomb`]) and the ideal acoustic cloak shell
//! ([`AcousticCloakShell`]).

use super::{
    AcousticCloakShell, ElasticMetamaterial, ReentrantHoneycomb, auxetic_poissons_from_geometry,
};
use std::f64::consts::PI;

/// Honeycomb cell with `l = 1`, thin ribs (`t/l = 0.1`) and a solid modulus of 1 GPa.
fn cell(theta_deg: f64, h_over_l: f64) -> ReentrantHoneycomb {
    ReentrantHoneycomb::new(theta_deg.to_radians(), 1.0, h_over_l, 0.1, 1.0e9)
}

/// Asserts `|actual - expected| <= rel * |expected|`.
fn assert_close(actual: f64, expected: f64, rel: f64) {
    let scale = expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        (actual - expected).abs() <= rel * scale,
        "actual {actual} differs from expected {expected} by more than a relative {rel}"
    );
}

/// Valid (closing) cells covering both signs of theta and several rib ratios.
fn valid_cells() -> Vec<ReentrantHoneycomb> {
    let mut cells = Vec::new();
    for &theta_deg in &[-40.0, -30.0, -20.0, -10.0, 10.0, 20.0, 30.0, 45.0, 60.0] {
        for &h_over_l in &[1.0, 1.5, 2.0, 3.0] {
            cells.push(cell(theta_deg, h_over_l));
        }
    }
    cells
}

// ---------------------------------------------------------------------------
// ReentrantHoneycomb
// ---------------------------------------------------------------------------

#[test]
fn regular_hexagon_has_unit_poisson_ratios() {
    // theta = +30 deg, h = l: Gibson-Ashby nu = +1 (not -1).
    let hexagon = cell(30.0, 1.0);
    assert_close(hexagon.poisson_ratio_xy(), 1.0, 1e-12);
    assert_close(hexagon.poisson_ratio_yx(), 1.0, 1e-12);
}

#[test]
fn reentrant_cell_has_minus_one_poisson_ratios() {
    // theta = -30 deg, h/l = 2 closes (h/l + sin theta = 1.5) and is auxetic.
    let reentrant = cell(-30.0, 2.0);
    assert_close(reentrant.poisson_ratio_xy(), -1.0, 1e-12);
    assert_close(reentrant.poisson_ratio_yx(), -1.0, 1e-12);
}

#[test]
fn asymmetric_cells_distinguish_the_two_loading_directions() {
    // theta = -20 deg, h/l = 2: nu_xy = cos^2 / ((h/l + sin) sin) = -1.5571867...,
    // nu_yx = (h/l + sin) sin / cos^2 = -0.6421837...
    let reentrant = cell(-20.0, 2.0);
    assert_close(
        reentrant.poisson_ratio_xy(),
        -1.557_186_745_330_010_1,
        1e-12,
    );
    assert_close(reentrant.poisson_ratio_yx(), -0.642_183_734_866_092, 1e-12);
    // theta = +20 deg, h/l = 1.5: both positive but different.
    let conventional = cell(20.0, 1.5);
    assert_close(
        conventional.poisson_ratio_xy(),
        1.401_604_790_366_811_8,
        1e-12,
    );
    assert_close(
        conventional.poisson_ratio_yx(),
        0.713_467_881_155_208_9,
        1e-12,
    );
}

#[test]
fn reentrant_sign_follows_the_cell_angle() {
    for &h_over_l in &[1.0, 1.5, 2.0, 3.0] {
        for &theta_deg in &[10.0, 20.0, 30.0, 45.0] {
            let conventional = cell(theta_deg, h_over_l);
            let reentrant = cell(-theta_deg.min(30.0), h_over_l);
            assert!(conventional.poisson_ratio_xy() > 0.0);
            assert!(conventional.poisson_ratio_yx() > 0.0);
            assert!(reentrant.poisson_ratio_xy() < 0.0);
            assert!(reentrant.poisson_ratio_yx() < 0.0);
        }
    }
}

#[test]
fn poisson_ratios_are_mutually_reciprocal() {
    for c in valid_cells() {
        let product = c.poisson_ratio_xy() * c.poisson_ratio_yx();
        assert_close(product, 1.0, 1e-12);
    }
}

#[test]
fn poisson_ratios_obey_reciprocity_with_the_moduli() {
    // nu_xy / E_x == nu_yx / E_y for every closing cell.
    for c in valid_cells() {
        let lhs = c.poisson_ratio_xy() / c.young_modulus_x();
        let rhs = c.poisson_ratio_yx() / c.young_modulus_y();
        assert_close(lhs, rhs, 1e-12);
    }
}

#[test]
fn all_three_poisson_ratio_sites_agree() {
    // The honeycomb struct, the geometry helper and the elastic metamaterial all
    // evaluate the Gibson-Ashby nu_21 = (h/l + sin) sin / cos^2 (loading along the
    // vertical ribs); none of them may carry the former sign error.
    for c in valid_cells() {
        let nu_yx = c.poisson_ratio_yx();
        let from_geometry = auxetic_poissons_from_geometry(c.theta.to_degrees(), c.l, c.h);
        assert_close(from_geometry, nu_yx, 1e-10);
        let mut metamaterial =
            ElasticMetamaterial::new(1.0e9, -0.3, 800.0, 3.0, c.theta, 200.0e9, 7850.0);
        let from_metamaterial = metamaterial.compute_auxetic_poisson(c.h / c.l);
        assert_close(from_metamaterial, nu_yx, 1e-10);
        assert_close(metamaterial.effective_poisson, nu_yx, 1e-10);
    }
}

#[test]
fn regular_hexagon_moduli_match_gibson_ashby() {
    // E* = (4/sqrt 3) Es (t/l)^3 in both directions, G* = E* / 4 = Es (t/l)^3 / sqrt 3.
    let hexagon = cell(30.0, 1.0);
    let t_over_l_cubed = 0.1_f64.powi(3);
    let expected_e = 4.0 / 3.0_f64.sqrt() * 1.0e9 * t_over_l_cubed;
    let expected_g = 1.0 / 3.0_f64.sqrt() * 1.0e9 * t_over_l_cubed;
    assert_close(hexagon.young_modulus_x(), expected_e, 1e-12);
    assert_close(hexagon.young_modulus_y(), expected_e, 1e-12);
    assert_close(hexagon.shear_modulus_xy(), expected_g, 1e-12);
    // In-plane isotropy: G = E / (2 (1 + nu)) with nu = 1.
    assert_close(
        hexagon.shear_modulus_xy(),
        hexagon.young_modulus_x() / (2.0 * (1.0 + hexagon.poisson_ratio_xy())),
        1e-12,
    );
}

#[test]
fn honeycomb_moduli_scale_with_the_cube_of_rib_thickness() {
    let thin = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 2.0, 0.05, 1.0e9);
    let thick = ReentrantHoneycomb::new(-PI / 6.0, 1.0, 2.0, 0.10, 1.0e9);
    assert_close(thick.young_modulus_x() / thin.young_modulus_x(), 8.0, 1e-12);
    assert_close(thick.young_modulus_y() / thin.young_modulus_y(), 8.0, 1e-12);
    assert_close(
        thick.shear_modulus_xy() / thin.shear_modulus_xy(),
        8.0,
        1e-12,
    );
    // The Poisson ratios are independent of the rib thickness.
    assert_close(thick.poisson_ratio_xy(), thin.poisson_ratio_xy(), 1e-14);
    assert_close(thick.poisson_ratio_yx(), thin.poisson_ratio_yx(), 1e-14);
}

#[test]
fn cells_that_cannot_close_are_reported_as_nan() {
    // theta = -30 deg with h/l = 0.3: h/l + sin(theta) = -0.2 < 0.
    let open = cell(-30.0, 0.3);
    assert!(open.poisson_ratio_xy().is_nan());
    assert!(open.poisson_ratio_yx().is_nan());
    assert!(open.young_modulus_x().is_nan());
    assert!(open.young_modulus_y().is_nan());
    assert!(open.shear_modulus_xy().is_nan());
}

#[test]
fn non_physical_honeycomb_inputs_are_reported_as_nan() {
    let bad_cells = [
        ReentrantHoneycomb::new(PI / 2.0, 1.0, 1.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(-PI / 2.0, 1.0, 2.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(PI / 6.0, 0.0, 1.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(PI / 6.0, -1.0, 1.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(PI / 6.0, 1.0, 0.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(PI / 6.0, 1.0, -2.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(f64::NAN, 1.0, 1.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(PI / 6.0, f64::INFINITY, 1.0, 0.1, 1.0e9),
        ReentrantHoneycomb::new(PI / 6.0, 1.0, f64::NAN, 0.1, 1.0e9),
    ];
    for c in &bad_cells {
        assert!(c.poisson_ratio_xy().is_nan(), "{c:?}");
        assert!(c.poisson_ratio_yx().is_nan(), "{c:?}");
        assert!(c.young_modulus_x().is_nan(), "{c:?}");
        assert!(c.young_modulus_y().is_nan(), "{c:?}");
        assert!(c.shear_modulus_xy().is_nan(), "{c:?}");
    }
}

#[test]
fn rectangular_cell_is_singular_only_in_the_x_direction() {
    // theta = 0: the horizontal ribs carry the x-load axially, so the pure-flexure
    // E_x and nu_xy are undefined, while nu_yx, E_y and G_xy stay finite.
    let rectangular = cell(0.0, 1.0);
    assert!(rectangular.poisson_ratio_xy().is_nan());
    assert!(rectangular.young_modulus_x().is_nan());
    assert_eq!(rectangular.poisson_ratio_yx(), 0.0);
    let t_over_l_cubed = 0.1_f64.powi(3);
    assert_close(rectangular.young_modulus_y(), 1.0e9 * t_over_l_cubed, 1e-12);
    assert_close(
        rectangular.shear_modulus_xy(),
        1.0e9 * t_over_l_cubed / 3.0,
        1e-12,
    );
}

// ---------------------------------------------------------------------------
// AcousticCloakShell
// ---------------------------------------------------------------------------

const R_INNER: f64 = 0.05;
const R_OUTER: f64 = 0.10;
const RHO_0: f64 = 1000.0;
const KAPPA_0: f64 = 2.2e9;

fn shell() -> AcousticCloakShell {
    AcousticCloakShell::new(R_INNER, R_OUTER, RHO_0, KAPPA_0, 10)
}

/// Radii strictly inside the shell.
fn shell_radii() -> Vec<f64> {
    (1..40)
        .map(|i| R_INNER + (R_OUTER - R_INNER) * f64::from(i) / 40.0)
        .collect()
}

#[test]
fn cloak_parameters_at_mid_shell_match_cummer_schurig() {
    let c = shell();
    // r = 0.075: r/(r - r1) = 3, (r - r1)/r = 1/3, ((r2 - r1)/r2)^2 = 1/4.
    assert_close(c.radial_density(0.075), 3000.0, 1e-12);
    assert_close(c.angular_density(0.075), 1000.0 / 3.0, 1e-12);
    assert_close(c.effective_modulus(0.075), 1.65e9, 1e-12);
}

#[test]
fn cloak_parameters_follow_the_transformation_acoustics_closed_forms() {
    let c = shell();
    let compression = (R_OUTER - R_INNER) / R_OUTER;
    for r in shell_radii() {
        assert_close(c.radial_density(r), RHO_0 * r / (r - R_INNER), 1e-12);
        assert_close(c.angular_density(r), RHO_0 * (r - R_INNER) / r, 1e-12);
        assert_close(
            c.effective_modulus(r),
            KAPPA_0 * compression * compression * r / (r - R_INNER),
            1e-12,
        );
    }
}

#[test]
fn cloak_density_product_is_background_density_squared() {
    let c = shell();
    for r in shell_radii() {
        assert_close(
            c.radial_density(r) * c.angular_density(r),
            RHO_0 * RHO_0,
            1e-12,
        );
    }
}

#[test]
fn cloak_radial_wave_speed_is_constant_and_compressed() {
    // sqrt(kappa / rho_r) = k c0 with k = (r2 - r1) / r2 at every radius.
    let c = shell();
    let compression = (R_OUTER - R_INNER) / R_OUTER;
    let c_0 = (KAPPA_0 / RHO_0).sqrt();
    for r in shell_radii() {
        let ratio = c.effective_modulus(r) / c.radial_density(r);
        assert_close(ratio, (KAPPA_0 / RHO_0) * compression * compression, 1e-12);
        assert_close(ratio.sqrt(), compression * c_0, 1e-12);
    }
}

#[test]
fn cloak_azimuthal_transit_time_equals_that_of_the_virtual_circle() {
    // Going once around radius r takes 2 pi r / c_theta with c_theta = sqrt(kappa / rho_theta);
    // in virtual space the same ray circles the radius s = (r - r1) / k at speed c0.
    let c = shell();
    let compression = (R_OUTER - R_INNER) / R_OUTER;
    let c_0 = (KAPPA_0 / RHO_0).sqrt();
    for r in shell_radii() {
        let c_theta = (c.effective_modulus(r) / c.angular_density(r)).sqrt();
        let physical = 2.0 * PI * r / c_theta;
        let virtual_radius = (r - R_INNER) / compression;
        let virtual_time = 2.0 * PI * virtual_radius / c_0;
        assert_close(physical, virtual_time, 1e-12);
    }
}

#[test]
fn cloak_radial_impedance_is_matched_at_the_outer_boundary() {
    // Normal-incidence impedance sqrt(kappa rho_r) must equal sqrt(kappa0 rho0) at r2.
    let c = shell();
    let just_inside = R_OUTER * (1.0 - 1e-12);
    let impedance = (c.effective_modulus(just_inside) * c.radial_density(just_inside)).sqrt();
    assert_close(impedance, (KAPPA_0 * RHO_0).sqrt(), 1e-9);
}

#[test]
fn cloak_parameters_are_singular_towards_the_inner_boundary() {
    let c = shell();
    let radii = shell_radii();
    // rho_r and kappa fall monotonically outwards, rho_theta rises from zero.
    for pair in radii.windows(2) {
        assert!(c.radial_density(pair[0]) > c.radial_density(pair[1]));
        assert!(c.effective_modulus(pair[0]) > c.effective_modulus(pair[1]));
        assert!(c.angular_density(pair[0]) < c.angular_density(pair[1]));
    }
    // r = r1 (1 + 1e-6): r / (r - r1) is about 1e6 / 1.
    let close_to_core = R_INNER * (1.0 + 1.0e-6);
    assert!(c.radial_density(close_to_core) > 9.0e8);
    assert!(c.effective_modulus(close_to_core) > 5.0e14);
    assert!(c.angular_density(close_to_core) < 2.0e-3);
}

#[test]
fn cloak_core_and_background_regions() {
    let c = shell();
    // Outside the shell the medium is the background fluid.
    for r in [R_OUTER, 0.11, 1.0, f64::INFINITY] {
        assert_eq!(c.radial_density(r), RHO_0);
        assert_eq!(c.angular_density(r), RHO_0);
        assert_eq!(c.effective_modulus(r), KAPPA_0);
    }
    // At and inside r1 the ideal parameters are replaced by their r -> r1+ limits.
    for r in [R_INNER, 0.02, 0.0, -1.0] {
        assert_eq!(c.radial_density(r), f64::INFINITY);
        assert_eq!(c.angular_density(r), 0.0);
        assert_eq!(c.effective_modulus(r), f64::INFINITY);
    }
}

#[test]
fn cloak_with_zero_inner_radius_is_the_identity_transformation() {
    // r1 = 0: nothing is hidden, the map is the identity and the shell is background fluid.
    let c = AcousticCloakShell::new(0.0, R_OUTER, RHO_0, KAPPA_0, 4);
    for r in [0.01, 0.05, 0.099] {
        assert_close(c.radial_density(r), RHO_0, 1e-12);
        assert_close(c.angular_density(r), RHO_0, 1e-12);
        assert_close(c.effective_modulus(r), KAPPA_0, 1e-12);
    }
}

#[test]
fn cloak_with_invalid_geometry_returns_nan() {
    let invalid = [
        AcousticCloakShell::new(0.10, 0.05, RHO_0, KAPPA_0, 4),
        AcousticCloakShell::new(0.05, 0.05, RHO_0, KAPPA_0, 4),
        AcousticCloakShell::new(-0.01, 0.10, RHO_0, KAPPA_0, 4),
        AcousticCloakShell::new(f64::NAN, 0.10, RHO_0, KAPPA_0, 4),
        AcousticCloakShell::new(0.05, f64::INFINITY, RHO_0, KAPPA_0, 4),
    ];
    for c in &invalid {
        assert!(c.radial_density(0.07).is_nan(), "{c:?}");
        assert!(c.angular_density(0.07).is_nan(), "{c:?}");
        assert!(c.effective_modulus(0.07).is_nan(), "{c:?}");
    }
    let c = shell();
    assert!(c.radial_density(f64::NAN).is_nan());
    assert!(c.angular_density(f64::NAN).is_nan());
    assert!(c.effective_modulus(f64::NAN).is_nan());
}

#[test]
fn cloak_layers_sample_the_ideal_parameters_at_mid_radius() {
    let c = shell();
    let layers = c.layer_properties();
    assert_eq!(layers.len(), 10);
    let dr = (R_OUTER - R_INNER) / 10.0;
    for (i, &(r, rho_r, rho_theta, kappa)) in layers.iter().enumerate() {
        assert_close(r, R_INNER + (i as f64 + 0.5) * dr, 1e-12);
        assert_close(rho_r, c.radial_density(r), 1e-14);
        assert_close(rho_theta, c.angular_density(r), 1e-14);
        assert_close(kappa, c.effective_modulus(r), 1e-14);
    }
    // First layer: r = 0.0525, r/(r - r1) = 21.
    assert_close(layers[0].1, 21_000.0, 1e-12);
    // Radial density and modulus decrease outwards, the angular density increases.
    for pair in layers.windows(2) {
        assert!(pair[0].1 > pair[1].1);
        assert!(pair[0].2 < pair[1].2);
        assert!(pair[0].3 > pair[1].3);
    }
}

#[test]
fn cloak_without_layers_has_no_layer_properties() {
    let c = AcousticCloakShell::new(R_INNER, R_OUTER, RHO_0, KAPPA_0, 0);
    assert!(c.layer_properties().is_empty());
}
