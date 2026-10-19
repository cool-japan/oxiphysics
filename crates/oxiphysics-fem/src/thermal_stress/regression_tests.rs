//! Regression tests for the 1-D thermoelastic analysis driver.
//!
//! The mechanical problem must be driven by the temperature rise over the
//! stress-free reference temperature and the stresses must be recovered from the
//! displacement solution:
//!
//! * a rod that can expand freely carries no stress and grows by `α ΔT L`,
//! * a rod clamped at both ends carries `σ = −E α ΔT`.

use super::*;
use crate::analysis::FemError;

const E: f64 = 200.0e9;
const ALPHA: f64 = 12.0e-6;
const AREA: f64 = 1.0e-4;
const CONDUCTIVITY: f64 = 50.0;

/// Rod of `n_elems` elements and length `length`, both ends held at the given
/// temperatures, reference temperature `t_ref`.
fn rod(
    n_elems: usize,
    length: f64,
    t_ends: (f64, f64),
    mechanical_bc: Vec<BoundaryCondition>,
    t_ref: f64,
) -> ThermalStressAnalysis {
    let mesh = ThermoelasticMesh::uniform_1d(n_elems, length, AREA);
    let thermal_bc = vec![
        BoundaryCondition::new(0, t_ends.0),
        BoundaryCondition::new(n_elems, t_ends.1),
    ];
    ThermalStressAnalysis::new(
        mesh,
        thermal_bc,
        mechanical_bc,
        E,
        0.3,
        ALPHA,
        CONDUCTIVITY,
        0.0,
        t_ref,
    )
}

fn fixed_left() -> Vec<BoundaryCondition> {
    vec![BoundaryCondition::new(0, 0.0)]
}

fn fixed_both(n_elems: usize) -> Vec<BoundaryCondition> {
    vec![
        BoundaryCondition::new(0, 0.0),
        BoundaryCondition::new(n_elems, 0.0),
    ]
}

fn assert_close(got: f64, want: f64, rel: f64, abs: f64, what: &str) {
    assert!(
        (got - want).abs() <= abs + rel * want.abs(),
        "{what}: got {got}, expected {want}"
    );
}

fn expect_invalid(result: Result<ThermoelasticResult, FemError>, what: &str) -> String {
    match result {
        Err(FemError::InvalidInput(message)) => message,
        other => panic!("{what}: expected InvalidInput, got {other:?}"),
    }
}

#[test]
fn free_bar_with_uniform_temperature_rise_expands_without_stress() {
    let (n, length, t_ref, t) = (4, 2.0, 300.0, 350.0);
    let result = rod(n, length, (t, t), fixed_left(), t_ref)
        .solve()
        .expect("well-posed analysis");
    let d_t = t - t_ref;
    let sigma_constrained = E * ALPHA * d_t;
    for (i, &u) in result.displacements.iter().enumerate() {
        let x = length * i as f64 / n as f64;
        assert_close(u, ALPHA * d_t * x, 1e-10, 1e-15, "displacement");
    }
    let tip = result.displacements[n];
    assert_close(tip, ALPHA * d_t * length, 1e-10, 0.0, "tip displacement");
    for (i, stress) in result.stresses.iter().enumerate() {
        assert!(
            stress.direct[0].abs() < 1e-8 * sigma_constrained,
            "element {i}: free expansion must be stress free, got {} Pa (constrained value {sigma_constrained})",
            stress.direct[0]
        );
        assert_close(stress.temperature, t, 1e-12, 0.0, "stored temperature");
    }
}

#[test]
fn bar_clamped_at_both_ends_carries_minus_e_alpha_delta_t() {
    let (n, t_ref, t) = (5, 300.0, 350.0);
    let result = rod(n, 1.0, (t, t), fixed_both(n), t_ref)
        .solve()
        .expect("well-posed analysis");
    let expected = -E * ALPHA * (t - t_ref);
    for (i, &u) in result.displacements.iter().enumerate() {
        assert!(u.abs() < 1e-12, "node {i} must not move, got {u}");
    }
    for (i, stress) in result.stresses.iter().enumerate() {
        assert_close(stress.direct[0], expected, 1e-9, 0.0, "clamped stress");
        assert!(
            stress.direct[0] < 0.0,
            "heating a clamped bar compresses it"
        );
        assert_close(
            stress.von_mises,
            expected.abs(),
            1e-9,
            0.0,
            &format!("von Mises of element {i}"),
        );
    }
    // Cooling reverses the sign.
    let cooled = rod(n, 1.0, (250.0, 250.0), fixed_both(n), t_ref)
        .solve()
        .expect("well-posed analysis");
    for stress in &cooled.stresses {
        assert_close(stress.direct[0], E * ALPHA * 50.0, 1e-9, 0.0, "cooling");
    }
}

#[test]
fn stresses_are_measured_from_the_reference_temperature_not_from_zero() {
    let n = 4;
    // Temperature equals the reference everywhere: nothing happens, however
    // large the absolute temperature is.
    let at_ref = rod(n, 1.0, (300.0, 300.0), fixed_both(n), 300.0)
        .solve()
        .expect("well-posed analysis");
    for (i, stress) in at_ref.stresses.iter().enumerate() {
        assert!(
            stress.direct[0].abs() < 1e-3,
            "element {i}: {}",
            stress.direct[0]
        );
    }
    assert!(at_ref.displacements.iter().all(|u| u.abs() < 1e-15));
    // With a zero reference the same field is a 300 K rise.
    let from_zero = rod(n, 1.0, (300.0, 300.0), fixed_both(n), 0.0)
        .solve()
        .expect("well-posed analysis");
    for stress in &from_zero.stresses {
        assert_close(stress.direct[0], -E * ALPHA * 300.0, 1e-9, 0.0, "σ");
    }
    // A free rod with a non-zero reference expands by α (T − T_ref) L.
    let free = rod(n, 1.0, (320.0, 320.0), fixed_left(), 300.0)
        .solve()
        .expect("well-posed analysis");
    assert_close(free.displacements[n], ALPHA * 20.0, 1e-10, 0.0, "tip");
}

#[test]
fn linear_temperature_gradient_free_expansion_is_stress_free() {
    // T(x) = 300 + 100 x / L, free at the right end: u(x) = α · 100 · x² / (2 L).
    let (n, length) = (8, 2.0);
    let result = rod(n, length, (300.0, 400.0), fixed_left(), 300.0)
        .solve()
        .expect("well-posed analysis");
    for (i, &t) in result.temperatures.iter().enumerate() {
        let x = length * i as f64 / n as f64;
        assert_close(t, 300.0 + 100.0 * x / length, 1e-10, 1e-9, "temperature");
        assert_close(
            result.displacements[i],
            ALPHA * 100.0 * x * x / (2.0 * length),
            1e-9,
            1e-15,
            "displacement",
        );
    }
    let scale = E * ALPHA * 100.0;
    for stress in &result.stresses {
        assert!(
            stress.direct[0].abs() < 1e-8 * scale,
            "σ = {}",
            stress.direct[0]
        );
    }
}

#[test]
fn clamped_bar_with_temperature_gradient_has_uniform_mean_stress() {
    // Both ends fixed, ΔT rising linearly from 0 to 100: compatibility gives a
    // uniform stress −E α mean(ΔT) = −E α · 50.
    let n = 10;
    let result = rod(n, 1.0, (300.0, 400.0), fixed_both(n), 300.0)
        .solve()
        .expect("well-posed analysis");
    for (i, stress) in result.stresses.iter().enumerate() {
        assert_close(
            stress.direct[0],
            -E * ALPHA * 50.0,
            1e-9,
            0.0,
            &format!("element {i}"),
        );
    }
    assert!(result.displacements.iter().all(|u| u.abs() < 1e-12));
}

#[test]
fn prescribed_end_displacement_is_honoured() {
    // Left end fixed, right end pushed by 1 mm against a 50 K rise (free
    // expansion 0.6 mm for L = 1 m): the mechanical strain is the difference.
    let n = 4;
    let mechanical_bc = vec![
        BoundaryCondition::new(0, 0.0),
        BoundaryCondition::new(n, 1.0e-3),
    ];
    let result = rod(n, 1.0, (350.0, 350.0), mechanical_bc, 300.0)
        .solve()
        .expect("well-posed analysis");
    assert_close(result.displacements[n], 1.0e-3, 1e-12, 0.0, "prescribed u");
    let expected = E * (1.0e-3 - ALPHA * 50.0);
    for stress in &result.stresses {
        assert_close(stress.direct[0], expected, 1e-9, 0.0, "σ");
    }
}

#[test]
fn internal_heat_generation_gives_the_parabolic_profile() {
    // k T'' + q = 0, T(0) = T(L) = 300:  T(x) = 300 + q x (L − x) / (2 k).
    let (n, length, q) = (8, 1.0, 1.0e6);
    let mut analysis = rod(n, length, (300.0, 300.0), fixed_left(), 300.0);
    analysis.heat_source = q;
    let result = analysis.solve().expect("well-posed analysis");
    for (i, &t) in result.temperatures.iter().enumerate() {
        let x = length * i as f64 / n as f64;
        let exact = 300.0 + q * x * (length - x) / (2.0 * CONDUCTIVITY);
        assert_close(t, exact, 1e-9, 1e-6, "temperature");
    }
    // Free expansion of the heated rod: tip = α Σ ΔT_mean h.
    let h = length / n as f64;
    let expected_tip: f64 = (0..n)
        .map(|e| 0.5 * (result.temperatures[e] + result.temperatures[e + 1]) - 300.0)
        .sum::<f64>()
        * ALPHA
        * h;
    assert_close(result.displacements[n], expected_tip, 1e-9, 0.0, "tip");
    let scale = E * ALPHA * 2500.0;
    for stress in &result.stresses {
        assert!(stress.direct[0].abs() < 1e-8 * scale);
    }
}

#[test]
fn reversed_element_orientation_gives_the_same_solution() {
    let n = 6;
    let normal = rod(n, 1.5, (310.0, 390.0), fixed_left(), 300.0);
    let mut reversed = normal.clone();
    for connectivity in reversed.mesh.connectivity.iter_mut() {
        connectivity.swap(0, 1);
    }
    let a = normal.solve().expect("normal");
    let b = reversed.solve().expect("reversed");
    for i in 0..=n {
        assert_close(b.displacements[i], a.displacements[i], 1e-10, 1e-15, "u");
        assert_close(b.temperatures[i], a.temperatures[i], 1e-10, 1e-12, "T");
    }
    for (sa, sb) in a.stresses.iter().zip(b.stresses.iter()) {
        assert!((sa.direct[0] - sb.direct[0]).abs() < 1e-3);
    }
}

#[test]
fn two_way_coupling_is_rejected_because_no_feedback_path_exists() {
    let n = 3;
    let mut analysis = rod(n, 1.0, (300.0, 350.0), fixed_left(), 300.0);
    assert_eq!(analysis.coupling_mode, CouplingMode::OneWay);
    assert!(analysis.solve().is_ok());
    analysis.coupling_mode = CouplingMode::TwoWay;
    let message = expect_invalid(analysis.solve(), "two-way coupling");
    assert!(message.contains("TwoWay"), "message: {message}");
}

#[test]
fn ill_posed_inputs_are_reported_instead_of_returning_garbage() {
    let n = 4;
    // no thermal / mechanical boundary condition: singular steady problems
    let mut no_thermal = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    no_thermal.thermal_bc.clear();
    let message = expect_invalid(no_thermal.solve(), "no thermal BC");
    assert!(message.contains("thermal"), "message: {message}");
    let no_mech = rod(n, 1.0, (300.0, 300.0), Vec::new(), 300.0);
    let message = expect_invalid(no_mech.solve(), "no mechanical BC");
    assert!(message.contains("mechanical"), "message: {message}");

    // node out of range, non-finite value, conflicting duplicates
    let out_of_range = rod(
        n,
        1.0,
        (300.0, 300.0),
        vec![BoundaryCondition::new(99, 0.0)],
        300.0,
    );
    expect_invalid(out_of_range.solve(), "BC outside the mesh");
    let non_finite = rod(
        n,
        1.0,
        (300.0, 300.0),
        vec![BoundaryCondition::new(0, f64::NAN)],
        300.0,
    );
    expect_invalid(non_finite.solve(), "NaN BC");
    let conflicting = rod(
        n,
        1.0,
        (300.0, 300.0),
        vec![
            BoundaryCondition::new(0, 0.0),
            BoundaryCondition::new(0, 1.0),
        ],
        300.0,
    );
    expect_invalid(conflicting.solve(), "conflicting BCs");
    // identical duplicates are harmless
    let duplicate = rod(
        n,
        1.0,
        (300.0, 300.0),
        vec![
            BoundaryCondition::new(0, 0.0),
            BoundaryCondition::new(0, 0.0),
        ],
        300.0,
    );
    assert!(duplicate.solve().is_ok());

    // unusable mesh / material data
    let mut zero_length = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    zero_length.mesh.node_positions[2] = zero_length.mesh.node_positions[1];
    expect_invalid(zero_length.solve(), "zero-length element");
    let mut bad_connectivity = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    bad_connectivity.mesh.connectivity[1] = [1, 42];
    expect_invalid(bad_connectivity.solve(), "connectivity out of range");
    let mut bad_area = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    bad_area.mesh.area = -1.0;
    expect_invalid(bad_area.solve(), "negative area");
    let mut bad_modulus = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    bad_modulus.e_modulus = 0.0;
    expect_invalid(bad_modulus.solve(), "zero modulus");
    let mut bad_conductivity = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    bad_conductivity.conductivity = f64::NAN;
    expect_invalid(bad_conductivity.solve(), "NaN conductivity");
    let mut bad_reference = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    bad_reference.reference_temp = f64::INFINITY;
    expect_invalid(bad_reference.solve(), "infinite reference temperature");
    let mut empty = rod(n, 1.0, (300.0, 300.0), fixed_left(), 300.0);
    empty.mesh.connectivity.clear();
    expect_invalid(empty.solve(), "no elements");
}
