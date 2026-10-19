//! Regression tests for the bainite, CALPHAD and phase-field models.
//!
//! Each test pins a behaviour that an earlier implementation got wrong:
//! `BainiteTransformation::time_to_fraction` returned NaN, `CalphadModel`
//! added the mixing-entropy term instead of subtracting it, and
//! `PhaseFieldModel2D` had no melting point (a liquid always solidified).

use super::functions::R_GAS;
use super::types::{BainiteTransformation, CalphadModel, PhaseFieldModel2D};

// ---------------------------------------------------------------------------
// BainiteTransformation::time_to_fraction
// ---------------------------------------------------------------------------

fn bainite() -> BainiteTransformation {
    BainiteTransformation::new(2.0, 1.0e-3, 80_000.0, 750.0, 500.0)
}

#[test]
fn bainite_time_to_fraction_is_finite_and_inverts_fraction_at() {
    let bt = bainite();
    let t = bt.time_to_fraction(0.5, 650.0);
    assert!(
        t.is_finite() && t > 0.0,
        "time to 50 % bainite must be finite and positive, got {t}"
    );
    let x = bt.fraction_at(t, 650.0);
    assert!(
        (x - 0.5).abs() < 1.0e-9,
        "fraction_at(time_to_fraction(0.5)) must be 0.5, got {x}"
    );
}

#[test]
fn bainite_time_to_fraction_round_trips_over_targets_and_temperatures() {
    let bt = bainite();
    for &temperature in &[500.0, 600.0, 700.0, 750.0] {
        for &target in &[0.01, 0.1, 0.25, 0.5, 0.9, 0.99, 0.999] {
            let t = bt.time_to_fraction(target, temperature);
            assert!(t.is_finite() && t > 0.0, "t({target}, {temperature}) = {t}");
            let x = bt.fraction_at(t, temperature);
            assert!(
                (x - target).abs() < 1.0e-9,
                "round trip failed at T={temperature}, X={target}: got {x}"
            );
        }
    }
}

#[test]
fn bainite_time_to_fraction_matches_closed_form_for_avrami_exponent_two() {
    let bt = bainite();
    let temperature = 650.0;
    let k = 1.0e-3 * (-80_000.0 / (R_GAS * temperature)).exp();
    for &target in &[0.1, 0.5, 0.9] {
        let expected = (-(1.0_f64 - target).ln() / k).sqrt();
        let got = bt.time_to_fraction(target, temperature);
        assert!(
            (got - expected).abs() <= 1.0e-12 * expected,
            "closed form {expected} vs {got} at X={target}"
        );
    }
}

#[test]
fn bainite_time_to_fraction_increases_with_target_and_decreases_with_temperature() {
    let bt = bainite();
    let t_10 = bt.time_to_fraction(0.1, 650.0);
    let t_50 = bt.time_to_fraction(0.5, 650.0);
    let t_90 = bt.time_to_fraction(0.9, 650.0);
    assert!(
        t_10 < t_50 && t_50 < t_90,
        "{t_10} < {t_50} < {t_90} expected"
    );
    // The Arrhenius rate constant grows with temperature, so the time shrinks.
    assert!(bt.time_to_fraction(0.5, 700.0) < bt.time_to_fraction(0.5, 600.0));
}

#[test]
fn bainite_time_to_fraction_edge_cases() {
    let bt = bainite();
    assert_eq!(
        bt.time_to_fraction(0.0, 650.0),
        0.0,
        "an empty target is reached immediately"
    );
    assert_eq!(bt.time_to_fraction(-0.2, 650.0), 0.0);
    assert!(
        bt.time_to_fraction(0.5, 800.0).is_infinite(),
        "no bainite forms above B_s"
    );
    assert!(
        bt.time_to_fraction(0.5, 1.0).is_infinite(),
        "vanishing rate constant never reaches the target"
    );
    let complete = bt.time_to_fraction(1.0, 650.0);
    assert!(
        complete.is_finite() && complete > bt.time_to_fraction(0.999, 650.0),
        "a target of 1 is clamped just below unity, got {complete}"
    );
}

// ---------------------------------------------------------------------------
// CalphadModel::gibbs_energy
// ---------------------------------------------------------------------------

fn ideal_solution() -> CalphadModel {
    CalphadModel::new(0.0, 0.0, 0.0, 0.0)
}

/// Central second difference of the Gibbs energy in `x_b`.
fn gibbs_second_derivative(model: &CalphadModel, x_b: f64, temperature: f64) -> f64 {
    let h = 1.0e-4;
    (model.gibbs_energy(x_b + h, temperature) - 2.0 * model.gibbs_energy(x_b, temperature)
        + model.gibbs_energy(x_b - h, temperature))
        / (h * h)
}

#[test]
fn calphad_ideal_solution_gibbs_energy_is_minus_t_s_mix() {
    let model = ideal_solution();
    let g = model.gibbs_energy(0.5, 1000.0);
    let expected = R_GAS * 1000.0 * 0.5_f64.ln();
    assert!(
        (g - expected).abs() < 1.0e-9,
        "G(x_B = 0.5, 1000 K) = R T ln 0.5 = {expected} J/mol, got {g}"
    );
    assert!(
        (g + 5763.2).abs() < 0.5,
        "G should be about -5763 J/mol, got {g}"
    );
}

#[test]
fn calphad_ideal_entropy_term_is_non_negative_and_lowers_g() {
    let model = ideal_solution();
    for &x_b in &[0.05, 0.3, 0.5, 0.8] {
        let ts = model.ideal_entropy(x_b, 800.0);
        assert!(
            ts > 0.0,
            "T S_mix must be positive at x_B = {x_b}, got {ts}"
        );
        assert!(
            (model.gibbs_energy(x_b, 800.0) + ts).abs() < 1.0e-9,
            "G = -T S_mix for the ideal solution at x_B = {x_b}"
        );
    }
    // At fixed composition a hotter ideal solution is more stable (lower G).
    assert!(model.gibbs_energy(0.5, 1000.0) < model.gibbs_energy(0.5, 500.0));
}

#[test]
fn calphad_pure_components_return_reference_energies() {
    let model = CalphadModel::new(-5000.0, -8000.0, 2000.0, 500.0);
    assert!((model.gibbs_energy(0.0, 900.0) - (-5000.0)).abs() < 1.0e-12);
    assert!((model.gibbs_energy(1.0, 900.0) - (-8000.0)).abs() < 1.0e-12);
}

#[test]
fn calphad_ideal_solution_is_convex_with_curvature_rt_over_x_a_x_b() {
    let model = ideal_solution();
    for &(x_b, temperature) in &[(0.5, 1000.0), (0.3, 700.0), (0.8, 1200.0)] {
        let curvature = gibbs_second_derivative(&model, x_b, temperature);
        let expected = R_GAS * temperature / (x_b * (1.0 - x_b));
        assert!(
            curvature > 0.0,
            "an ideal solution is stable (G'' > 0), got {curvature}"
        );
        assert!(
            (curvature - expected).abs() < 1.0e-5 * expected,
            "G'' = R T / (x_A x_B): expected {expected}, got {curvature}"
        );
    }
}

#[test]
fn calphad_regular_solution_has_a_miscibility_gap_below_the_critical_temperature() {
    // L0 = 20 kJ/mol  ->  T_c = L0 / (2 R) ~ 1203 K.
    let model = CalphadModel::new(0.0, 0.0, 20_000.0, 0.0);
    let spinodal = gibbs_second_derivative(&model, 0.5, 600.0);
    let stable = gibbs_second_derivative(&model, 0.5, 1500.0);
    let expected_spinodal = 4.0 * R_GAS * 600.0 - 2.0 * 20_000.0;
    let expected_stable = 4.0 * R_GAS * 1500.0 - 2.0 * 20_000.0;
    assert!(
        spinodal < 0.0,
        "G'' < 0 inside the spinodal, got {spinodal}"
    );
    assert!(stable > 0.0, "G'' > 0 above T_c, got {stable}");
    assert!((spinodal - expected_spinodal).abs() < 1.0e-5 * expected_spinodal.abs());
    assert!((stable - expected_stable).abs() < 1.0e-5 * expected_stable.abs());
}

#[test]
fn calphad_chemical_potential_of_ideal_solution_is_rt_ln_x_a() {
    let model = ideal_solution();
    let mu = model.chemical_potential_a(0.3, 1000.0);
    let expected = R_GAS * 1000.0 * 0.7_f64.ln();
    assert!(
        (mu - expected).abs() < 1.0e-3,
        "mu_A = R T ln x_A = {expected} J/mol, got {mu}"
    );
}

// ---------------------------------------------------------------------------
// PhaseFieldModel2D melting point
// ---------------------------------------------------------------------------

const T_MELT: f64 = 1000.0;

fn uniform_model(nx: usize, ny: usize, temperature: f64) -> PhaseFieldModel2D {
    let mut pf = PhaseFieldModel2D::new(nx, ny, 2.0, 1.0, 100.0, T_MELT, 1.0);
    pf.temperature.iter_mut().for_each(|t| *t = temperature);
    pf
}

fn seeded_model(temperature: f64) -> PhaseFieldModel2D {
    let mut pf = uniform_model(48, 48, temperature);
    pf.initialize_seed(24.0, 24.0, 6.0);
    pf
}

/// Solid fraction sampled every `chunk` steps of size `dt`.
fn solid_fraction_history(
    pf: &mut PhaseFieldModel2D,
    dt: f64,
    chunk: usize,
    samples: usize,
) -> Vec<f64> {
    let mut history = vec![pf.solid_fraction()];
    for _ in 0..samples {
        for _ in 0..chunk {
            pf.step(dt, 1.0);
        }
        history.push(pf.solid_fraction());
    }
    history
}

#[test]
fn phase_field_pure_liquid_never_solidifies_spontaneously() {
    for &temperature in &[300.0, T_MELT - 1.0, T_MELT, T_MELT + 200.0] {
        let mut pf = uniform_model(16, 16, temperature);
        for _ in 0..200 {
            pf.step(1.0e-3, 1.0);
        }
        assert!(
            pf.phi.iter().all(|&p| p == 0.0),
            "pure liquid must stay liquid at {temperature} K"
        );
        assert_eq!(pf.solid_fraction(), 0.0);
    }
}

#[test]
fn phase_field_pure_solid_never_melts_spontaneously() {
    for &temperature in &[300.0, T_MELT, T_MELT + 500.0] {
        let mut pf = uniform_model(16, 16, temperature);
        pf.phi.iter_mut().for_each(|p| *p = 1.0);
        for _ in 0..200 {
            pf.step(1.0e-3, 1.0);
        }
        assert!(
            pf.phi.iter().all(|&p| p == 1.0),
            "pure solid must stay solid at {temperature} K"
        );
    }
}

#[test]
fn phase_field_seed_grows_below_and_shrinks_above_the_melting_point() {
    let growing = solid_fraction_history(&mut seeded_model(900.0), 1.0e-3, 100, 3);
    assert!(
        growing.windows(2).all(|w| w[1] > w[0]),
        "an undercooled seed must grow monotonically: {growing:?}"
    );
    assert!(
        growing[3] > 2.0 * growing[0],
        "growth should be substantial"
    );

    let shrinking = solid_fraction_history(&mut seeded_model(1100.0), 1.0e-3, 100, 3);
    assert!(
        shrinking.windows(2).all(|w| w[1] < w[0]),
        "a superheated seed must shrink monotonically: {shrinking:?}"
    );
    assert!(
        shrinking[3] < 0.5 * shrinking[0],
        "melting should be substantial"
    );
}

#[test]
fn phase_field_seed_only_shrinks_slowly_by_curvature_at_the_melting_point() {
    let at_melting = solid_fraction_history(&mut seeded_model(T_MELT), 1.0e-3, 100, 3);
    let superheated = solid_fraction_history(&mut seeded_model(1100.0), 1.0e-3, 100, 3);
    assert!(
        at_melting[3] < at_melting[0],
        "curvature makes a seed recede even without a driving force"
    );
    assert!(
        at_melting[3] > 0.9 * at_melting[0],
        "no thermal driving force: the recession must stay small, got {at_melting:?}"
    );
    assert!(at_melting[3] > superheated[3]);
}

#[test]
fn phase_field_melting_temperature_decides_the_direction_at_fixed_temperature() {
    // Same 1000 K temperature field; only the melting point changes.
    let history = |t_melt: f64| {
        let mut pf = PhaseFieldModel2D::new(48, 48, 2.0, 1.0, 100.0, t_melt, 1.0);
        pf.temperature.iter_mut().for_each(|t| *t = 1000.0);
        pf.initialize_seed(24.0, 24.0, 6.0);
        solid_fraction_history(&mut pf, 1.0e-3, 100, 2)
    };
    let below_melting = history(1100.0);
    let above_melting = history(900.0);
    assert!(below_melting[2] > below_melting[0], "{below_melting:?}");
    assert!(above_melting[2] < above_melting[0], "{above_melting:?}");
}

#[test]
fn phase_field_non_positive_melting_temperature_disables_the_driving_force() {
    let mut pf = PhaseFieldModel2D::new(48, 48, 2.0, 1.0, 100.0, 0.0, 1.0);
    pf.initialize_seed(24.0, 24.0, 6.0);
    let history = solid_fraction_history(&mut pf, 1.0e-3, 100, 2);
    assert!(
        history[2] < history[0],
        "without a thermal driving force only curvature acts: {history:?}"
    );
}

/// A solid slab with two flat interfaces (periodic in both directions).  The
/// slab is wide enough (64 cells, interface width about 3 cells) that the
/// non-periodic tails of the tanh profile are of order 1e-10.
fn slab_model(t_melt: f64, temperature: f64, latent_heat: f64, coupling: f64) -> PhaseFieldModel2D {
    let (nx, ny) = (128, 4);
    let width = 2.0;
    let mut pf = PhaseFieldModel2D::new(nx, ny, width, 1.0, latent_heat, t_melt, coupling);
    pf.temperature.iter_mut().for_each(|t| *t = temperature);
    // Equilibrium profile of the double well f = phi^2 (1 - phi)^2 with the
    // gradient coefficient W^2: phi = (1 + tanh(x / (sqrt(2) W))) / 2.
    let delta = std::f64::consts::SQRT_2 * width;
    for j in 0..ny {
        for i in 0..nx {
            let x = i as f64;
            pf.phi[j * nx + i] = 0.5 * (((x - 32.0) / delta).tanh() - ((x - 96.0) / delta).tanh());
        }
    }
    pf
}

/// Mean interface speed (cells per unit time) of the slab's two interfaces.
fn slab_interface_speed(t_melt: f64, temperature: f64, latent_heat: f64, coupling: f64) -> f64 {
    let mut pf = slab_model(t_melt, temperature, latent_heat, coupling);
    let f0 = pf.solid_fraction();
    let (dt, steps) = (0.02, 500);
    for _ in 0..steps {
        pf.step(dt, 1.0);
    }
    (pf.solid_fraction() - f0) * pf.nx as f64 / 2.0 / (dt * steps as f64)
}

#[test]
fn phase_field_flat_interface_is_stationary_at_the_melting_temperature() {
    let speed = slab_interface_speed(T_MELT, T_MELT, 10.0, 1.0);
    assert!(
        speed.abs() < 1.0e-9,
        "a flat interface at T_m must not move, speed = {speed}"
    );
}

#[test]
fn phase_field_flat_interface_speed_matches_the_sharp_interface_limit() {
    // Small driving force dG = lambda Q (T_m - T) / T_m: the Allen-Cahn front
    // speed is v = L dG / int(phi'^2 dx) = 3 sqrt(2) L W dG.
    let cases = [(10.0, 999.0, 1.0), (10.0, 998.0, 2.0)];
    for &(q, temperature, lambda) in &cases {
        let d_g = lambda * q * (T_MELT - temperature) / T_MELT;
        let predicted = 3.0 * std::f64::consts::SQRT_2 * 1.0 * 2.0 * d_g;
        let measured = slab_interface_speed(T_MELT, temperature, q, lambda);
        assert!(
            (measured - predicted).abs() < 0.06 * predicted,
            "front speed {measured} vs sharp-interface {predicted} (dG = {d_g})"
        );
        // Above the melting point the same magnitude, opposite direction.
        let melting = slab_interface_speed(T_MELT, 2.0 * T_MELT - temperature, q, lambda);
        assert!(
            (melting + measured).abs() < 1.0e-6 * measured,
            "melting speed {melting} must mirror the freezing speed {measured}"
        );
    }
}
