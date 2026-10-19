// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the J2 radial return (plastic strain bookkeeping in
//! engineering shear), the J2 consistent tangent (checked against a finite
//! difference of the stress update) and the Drucker-Prager apex return.

use oxiphysics_materials::constitutive::{ConstitutiveModel, J2State};
use oxiphysics_materials::plasticity::{
    DruckerPragerApexReturn, DruckerPragerModel, J2ConsistentTangent, J2ReturnMapping,
};

const YOUNG: f64 = 200.0e9;
const POISSON: f64 = 0.3;
const YIELD: f64 = 250.0e6;
const HARDENING: f64 = 2.0e9;

fn steel() -> J2ReturnMapping {
    J2ReturnMapping::from_young_poisson(YOUNG, POISSON, YIELD, HARDENING)
}

fn max_abs(values: &[f64]) -> f64 {
    values.iter().fold(0.0_f64, |acc, v| acc.max(v.abs()))
}

// ---------------------------------------------------------------------------
// F3027: plastic shear strain is engineering shear
// ---------------------------------------------------------------------------

#[test]
fn return_map_pure_shear_stores_engineering_shear_plastic_strain() {
    let rm = steel();
    let g = rm.shear_modulus;
    let tau = 300.0e6;
    let trial = [0.0, 0.0, 0.0, 0.0, 0.0, tau];
    let (stress, dep, dg) = rm.return_map(&trial, 0.0);

    let q = 3.0_f64.sqrt() * tau;
    let dg_expected = (q - YIELD) / (3.0 * g + HARDENING);
    assert!(
        (dg - dg_expected).abs() <= 1.0e-12 * dg_expected,
        "plastic multiplier {dg} vs {dg_expected}"
    );
    // gamma_p = 2 * eps_p with eps_p = dg * (3/2) * tau / q  ->  sqrt(3) * dg.
    let gamma_p = 3.0_f64.sqrt() * dg_expected;
    assert!(
        (dep[5] - gamma_p).abs() <= 1.0e-12 * gamma_p,
        "engineering plastic shear {} vs {gamma_p}",
        dep[5]
    );
    // The shear stress drops by exactly G * gamma_p.
    assert!(
        (stress[5] - (tau - g * dep[5])).abs() <= 1.0e-9 * tau,
        "tau_new = {} vs tau - G gamma_p = {}",
        stress[5],
        tau - g * dep[5]
    );
    // And the stress sits on the hardened yield surface.
    let vm = J2ReturnMapping::von_mises(&stress);
    let yield_new = YIELD + HARDENING * dg;
    assert!((vm - yield_new).abs() <= 1.0e-9 * yield_new);
}

#[test]
fn return_map_uniaxial_plastic_flow_is_isochoric_with_axial_strain_equal_to_the_multiplier() {
    let rm = steel();
    let trial = [500.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
    let (_, dep, dg) = rm.return_map(&trial, 0.0);
    assert!((dep[0] - dg).abs() <= 1.0e-12 * dg);
    assert!((dep[1] + 0.5 * dg).abs() <= 1.0e-12 * dg);
    assert!((dep[2] + 0.5 * dg).abs() <= 1.0e-12 * dg);
    assert!((dep[0] + dep[1] + dep[2]).abs() <= 1.0e-12 * dg);
    assert!(dep[3..].iter().all(|&s| s == 0.0));
}

#[test]
fn stress_update_on_a_converged_plastic_state_is_elastic_and_reproduces_the_stress() {
    let rm = steel();
    let g = rm.shear_modulus;
    // Shear dominated strain, far beyond yield.
    let strain = [1.0e-3, -4.0e-4, -2.0e-4, 2.0e-3, -1.0e-3, 6.0e-3];
    let first = rm.stress_update(&strain, &J2State::default(), 1.0e-3);
    assert!(
        first.state.equiv_plastic_strain > 1.0e-4,
        "the strain path must yield"
    );

    let second = rm.stress_update(&strain, &first.state, 1.0e-3);
    let increment = second.state.equiv_plastic_strain - first.state.equiv_plastic_strain;
    assert!(
        increment.abs() < 1.0e-12,
        "re-evaluating a converged state must not flow again: {increment}"
    );
    let scale = max_abs(&first.stress);
    for k in 0..6 {
        assert!(
            (second.stress[k] - first.stress[k]).abs() <= 1.0e-10 * scale,
            "stress component {k} changed: {} -> {}",
            first.stress[k],
            second.stress[k]
        );
        assert!(
            (second.state.plastic_strain[k] - first.state.plastic_strain[k]).abs() < 1.0e-12,
            "plastic strain component {k} changed"
        );
    }
    // sigma_shear = G * (gamma - gamma_p) for each engineering shear component.
    for (k, &gamma) in strain.iter().enumerate().skip(3) {
        let expected = g * (gamma - first.state.plastic_strain[k]);
        assert!(
            (first.stress[k] - expected).abs() <= 1.0e-9 * scale,
            "shear {k}: {} vs G (gamma - gamma_p) = {expected}",
            first.stress[k]
        );
    }
}

#[test]
fn stress_update_is_independent_of_the_step_count_for_proportional_loading() {
    // Radial return is exact for proportional strain paths with linear isotropic
    // hardening, so one big step and six small steps must agree.
    let rm = steel();
    let direction = [0.4e-3, -0.1e-3, -0.3e-3, 0.5e-3, -0.2e-3, 1.0e-3];
    let scaled = |factor: f64| direction.map(|d| d * factor);

    let single = rm.stress_update(&scaled(6.0), &J2State::default(), 1.0e-3);
    let mut state = J2State::default();
    let mut last = None;
    for step in 1..=6 {
        let response = rm.stress_update(&scaled(f64::from(step)), &state, 1.0e-3);
        state = response.state.clone();
        last = Some(response);
    }
    let stepped = last.expect("six steps were taken");

    assert!(single.state.equiv_plastic_strain > 1.0e-4);
    let scale = max_abs(&single.stress);
    for k in 0..6 {
        assert!(
            (single.stress[k] - stepped.stress[k]).abs() <= 1.0e-9 * scale,
            "stress {k}: {} vs {}",
            single.stress[k],
            stepped.stress[k]
        );
        assert!(
            (single.state.plastic_strain[k] - stepped.state.plastic_strain[k]).abs() <= 1.0e-12,
            "plastic strain {k}: {} vs {}",
            single.state.plastic_strain[k],
            stepped.state.plastic_strain[k]
        );
    }
    assert!(
        (single.state.equiv_plastic_strain - stepped.state.equiv_plastic_strain).abs() <= 1.0e-12
    );
}

// ---------------------------------------------------------------------------
// F3028: consistent tangent vs finite differences
// ---------------------------------------------------------------------------

/// Central finite difference of the stress with respect to the strain vector.
fn finite_difference_tangent(
    rm: &J2ReturnMapping,
    strain: &[f64; 6],
    state: &J2State,
    h: f64,
) -> [[f64; 6]; 6] {
    let mut d = [[0.0_f64; 6]; 6];
    for j in 0..6 {
        let mut plus = *strain;
        let mut minus = *strain;
        plus[j] += h;
        minus[j] -= h;
        let resp_plus = rm.stress_update(&plus, state, 1.0e-3);
        let resp_minus = rm.stress_update(&minus, state, 1.0e-3);
        assert!(
            resp_plus.state.equiv_plastic_strain > state.equiv_plastic_strain
                && resp_minus.state.equiv_plastic_strain > state.equiv_plastic_strain,
            "both perturbed states must stay plastic (column {j})"
        );
        for (i, row) in d.iter_mut().enumerate() {
            row[j] = (resp_plus.stress[i] - resp_minus.stress[i]) / (2.0 * h);
        }
    }
    d
}

fn assert_tangent_matches_finite_difference(strain: &[f64; 6], state: &J2State) {
    let rm = steel();
    let response = rm.stress_update(strain, state, 1.0e-3);
    assert!(
        response.state.equiv_plastic_strain > state.equiv_plastic_strain,
        "the step must be plastic"
    );
    let fd = finite_difference_tangent(&rm, strain, state, 1.0e-7);
    let flat: Vec<f64> = response.tangent.iter().flatten().copied().collect();
    let scale = max_abs(&flat);
    let worst = fd
        .iter()
        .flatten()
        .zip(response.tangent.iter().flatten())
        .fold(0.0_f64, |acc, (a, b)| acc.max((a - b).abs()));
    assert!(
        worst <= 1.0e-6 * scale,
        "consistent tangent deviates from the finite difference: {worst:.3e} vs scale {scale:.3e}"
    );
    // The algorithmic tangent of associative J2 flow is symmetric.
    for (i, row) in response.tangent.iter().enumerate() {
        for (j, &entry) in row.iter().enumerate().take(i) {
            assert!(
                (entry - response.tangent[j][i]).abs() <= 1.0e-12 * scale,
                "tangent not symmetric at ({i},{j})"
            );
        }
    }
}

#[test]
fn consistent_tangent_matches_finite_difference_from_the_virgin_state() {
    // All six strain components non-zero: exercises the normal/shear cross terms.
    let strain = [2.0e-3, -1.1e-3, 0.4e-3, 1.5e-3, -2.5e-3, 3.5e-3];
    assert_tangent_matches_finite_difference(&strain, &J2State::default());
}

#[test]
fn consistent_tangent_matches_finite_difference_after_hardening() {
    let rm = steel();
    let first_strain = [1.0e-3, -0.4e-3, -0.2e-3, 1.2e-3, 0.6e-3, -2.0e-3];
    let first = rm.stress_update(&first_strain, &J2State::default(), 1.0e-3);
    assert!(first.state.equiv_plastic_strain > 0.0);
    let second_strain = [2.4e-3, -0.9e-3, -0.3e-3, 2.2e-3, 1.1e-3, -3.6e-3];
    assert_tangent_matches_finite_difference(&second_strain, &first.state);
}

#[test]
fn consistent_tangent_for_pure_shear_has_the_hardened_shear_modulus() {
    // Pure shear: n_55^2 = 1/2 and the plastic-multiplier terms cancel, leaving
    // D_55 = G - 3 G^2 / (3 G + H) = G H / (3 G + H) = d tau / d gamma, which also
    // follows from tau = G (gamma - gamma_p) with the hardened yield condition.
    let rm = steel();
    let g = rm.shear_modulus;
    let tan = J2ConsistentTangent::from_young_poisson(YOUNG, POISSON, HARDENING);
    let trial = [0.0, 0.0, 0.0, 0.0, 0.0, 300.0e6];
    let (_, _, dg) = rm.return_map(&trial, 0.0);
    let c = tan.compute_consistent_tangent(&trial, dg);
    let expected = g * HARDENING / (3.0 * g + HARDENING);
    assert!(
        (c[5 * 6 + 5] - expected).abs() <= 1.0e-9 * expected,
        "pure-shear tangent {} vs G H / (3G + H) = {expected}",
        c[5 * 6 + 5]
    );
}

// ---------------------------------------------------------------------------
// F3029: Drucker-Prager apex return
// ---------------------------------------------------------------------------

fn apex_return() -> DruckerPragerApexReturn {
    DruckerPragerApexReturn::new(0.2, 100.0e3, 30.0e6, 15.0e6)
}

#[test]
fn apex_return_is_needed_for_hydrostatic_tension_beyond_the_apex() {
    let apr = apex_return();
    let hydrostatic_tension = [1.0e6, 1.0e6, 1.0e6, 0.0, 0.0, 0.0];
    assert!(apr.needs_apex_return(&hydrostatic_tension));
}

#[test]
fn apex_return_is_not_needed_for_compressive_states() {
    let apr = apex_return();
    let uniaxial_compression = [-1.0e6, 0.0, 0.0, 0.0, 0.0, 0.0];
    assert!(
        !apr.needs_apex_return(&uniaxial_compression),
        "uniaxial compression returns to the smooth cone"
    );
    let deep_hydrostatic_compression = [-5.0e6, -5.0e6, -5.0e6, 0.0, 0.0, 0.0];
    assert!(!apr.needs_apex_return(&deep_hydrostatic_compression));
    // Pure shear at zero mean stress is far from the apex as well.
    let shear = [0.0, 0.0, 0.0, 2.0e6, 0.0, 0.0];
    assert!(!apr.needs_apex_return(&shear));
}

#[test]
fn apex_return_is_not_needed_for_elastic_states_even_in_tension() {
    let apr = apex_return();
    // p = 100 kPa < p_apex = 166.7 kPa: f_tr = 0.2 * 300 kPa - 100 kPa < 0.
    let weak_tension = [1.0e5, 1.0e5, 1.0e5, 0.0, 0.0, 0.0];
    assert!(!apr.needs_apex_return(&weak_tension));
    assert_eq!(apr.compute_apex_multiplier(&weak_tension), 0.0);
}

#[test]
fn apex_region_boundary_matches_the_closed_form() {
    // sqrt(J2)* = G (alpha I1 - k) / (9 alpha^2 K): the cone return would leave
    // exactly sqrt(J2) = 0 on this line.
    let apr = apex_return();
    let p = 1.0e6;
    let i1 = 3.0 * p;
    let boundary = apr.shear_modulus * (apr.alpha * i1 - apr.k_dp)
        / (9.0 * apr.alpha * apr.alpha * apr.bulk_modulus);
    assert!((boundary - 694_444.444).abs() < 1.0, "boundary {boundary}");
    // For a pure shear stress tau on top of the pressure, sqrt(J2) = |tau|.
    let state = |tau: f64| [p, p, p, tau, 0.0, 0.0];
    assert!(apr.needs_apex_return(&state(boundary * (1.0 - 1.0e-6))));
    assert!(!apr.needs_apex_return(&state(boundary * (1.0 + 1.0e-6))));
    assert!(apr.needs_apex_return(&state(600.0e3)));
    assert!(!apr.needs_apex_return(&state(800.0e3)));
}

#[test]
fn apex_multiplier_uses_the_cone_consistency_denominator() {
    let apr = apex_return();
    let hydrostatic_tension = [1.0e6, 1.0e6, 1.0e6, 0.0, 0.0, 0.0];
    let f_tr = 0.2 * 3.0e6 - 100.0e3;
    let denominator = 15.0e6 + 9.0 * 0.2 * 0.2 * 30.0e6;
    let expected = f_tr / denominator;
    let multiplier = apr.compute_apex_multiplier(&hydrostatic_tension);
    assert!(
        (multiplier - expected).abs() <= 1.0e-12 * expected,
        "multiplier {multiplier} vs {expected}"
    );
    assert!((expected - 0.019_379_8).abs() < 1.0e-6);
    // The cone return would drive sqrt(J2) negative here, hence the apex return.
    assert!(15.0e6 * multiplier > 0.0);
}

#[test]
fn apex_volumetric_plastic_strain_resets_the_mean_stress_to_the_apex() {
    let apr = apex_return();
    let trial = [1.0e6, 1.0e6, 1.0e6, 0.0, 0.0, 0.0];
    let plastic_volumetric = apr.apex_volumetric_plastic_strain(&trial);
    let expected = (1.0e6 - 100.0e3 / 0.6) / 30.0e6;
    assert!((plastic_volumetric - expected).abs() <= 1.0e-12 * expected);
    let p_after = 1.0e6 - apr.bulk_modulus * plastic_volumetric;
    assert!((p_after - apr.apex_pressure()).abs() <= 1.0e-6);
    let returned = apr.compute_apex_return(&trial);
    assert!((returned[0] - apr.apex_pressure()).abs() <= 1.0e-9);
    assert_eq!(&returned[3..], &[0.0, 0.0, 0.0]);
}

#[test]
fn apex_return_from_a_mohr_coulomb_matched_model() {
    let phi = 30.0_f64.to_radians();
    let (cohesion, young, nu) = (50.0e3, 50.0e6, 0.3);
    let dp = DruckerPragerModel::new(
        phi,
        cohesion,
        young / (2.0 * (1.0 + nu)),
        young / (3.0 * (1.0 - 2.0 * nu)),
    );
    let apr = DruckerPragerApexReturn::from_dp_model(&dp);
    let p_apex = apr.apex_pressure();
    assert!(
        p_apex > 0.0,
        "the tensile apex lies at positive mean stress"
    );
    let beyond_apex = [2.0 * p_apex, 2.0 * p_apex, 2.0 * p_apex, 0.0, 0.0, 0.0];
    assert!(apr.needs_apex_return(&beyond_apex));
    let compressed = [-p_apex, -p_apex, -p_apex, 0.0, 0.0, 0.0];
    assert!(!apr.needs_apex_return(&compressed));
}
