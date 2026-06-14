//! Integration tests for [`GearConstraint`].

use oxiphysics_constraints::mechanism::{BodyState, GearConstraint};

// ---------------------------------------------------------------------------
// Test 1 — rigid gear ratio correctness
// ---------------------------------------------------------------------------

/// Rigid gear (ratio=2, no backlash): after 1000 TGS iterations with angle
/// integration, the position drift |θ₁ − 2·θ₂| must be < 1e-4 and the
/// velocity ratio |ω_Bz − 2·ω_Dz| < 1e-6.
#[test]
fn rigid_gear_ratio_correctness() {
    let gc = GearConstraint::new(2.0, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);

    // A and C are static (ground / housing)
    let inv_ia = [0.0f64; 3];
    let inv_ic = [0.0f64; 3];
    // B and D are unit flywheels
    let inv_ib = [1.0f64; 3];
    let inv_id = [1.0f64; 3];

    let mut omega_a = [0.0f64; 3];
    let mut omega_b = [0.0f64, 0.0, 1.0]; // ω_Bz = 1 rad/s
    let mut omega_c = [0.0f64; 3];
    let mut omega_d = [0.0f64, 0.0, 0.3]; // ω_Dz = 0.3 rad/s (not at ratio yet)

    let dt = 1.0 / 60.0;
    let mut theta1 = 0.0f64;
    let mut theta2 = 0.0f64;

    for _ in 0..1000 {
        gc.solve_iteration(
            &mut BodyState::new(&mut omega_a, &inv_ia),
            &mut BodyState::new(&mut omega_b, &inv_ib),
            &mut BodyState::new(&mut omega_c, &inv_ic),
            &mut BodyState::new(&mut omega_d, &inv_id),
            theta1,
            theta2,
            dt,
        );
        // Integrate angles from the corrected velocities
        theta1 += omega_b[2] * dt;
        theta2 += omega_d[2] * dt;
    }

    let drift = (theta1 - 2.0 * theta2).abs();
    let vel_err = (omega_b[2] - 2.0 * omega_d[2]).abs();

    assert!(
        drift < 1e-4,
        "Position drift too large: |θ₁ - 2·θ₂| = {drift} >= 1e-4"
    );
    assert!(
        vel_err < 1e-6,
        "Velocity ratio error too large: |ω_Bz - 2·ω_Dz| = {vel_err} >= 1e-6"
    );
}

// ---------------------------------------------------------------------------
// Test 2 — backlash dead zone
// ---------------------------------------------------------------------------

/// With backlash=0.1, body D must remain stationary while the drift is inside
/// the dead zone, then begin moving once the threshold is crossed.
#[test]
fn backlash_dead_zone() {
    let gc = GearConstraint::new(2.0, 0.1, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);

    let inv_ia = [0.0f64; 3]; // A static
    let inv_ic = [0.0f64; 3]; // C static
    let inv_ib = [1.0f64; 3]; // B flywheel
    let inv_id = [1.0f64; 3]; // D flywheel

    let mut omega_a = [0.0f64; 3];
    let mut omega_b = [0.0f64, 0.0, 0.5]; // drives theta1 upward
    let mut omega_c = [0.0f64; 3];
    let mut omega_d = [0.0f64; 3]; // D starts at rest

    let dt = 1.0 / 60.0;
    let mut theta1 = 0.0f64;
    let mut theta2 = 0.0f64;

    // With dt=1/60 and omega_b=0.5, each iteration: Δθ₁ ≈ 0.5/60 ≈ 0.00833
    // Dead zone: |θ₁| < 0.1, so approximately 12 iterations to cross.
    //
    // We check after 10 iterations (still inside dead zone):
    let iters_in_zone = 10usize;
    for _ in 0..iters_in_zone {
        gc.solve_iteration(
            &mut BodyState::new(&mut omega_a, &inv_ia),
            &mut BodyState::new(&mut omega_b, &inv_ib),
            &mut BodyState::new(&mut omega_c, &inv_ic),
            &mut BodyState::new(&mut omega_d, &inv_id),
            theta1,
            theta2,
            dt,
        );
        theta1 += omega_b[2] * dt;
        theta2 += omega_d[2] * dt;
    }

    // Inside dead zone: D must not have moved
    assert!(
        theta1.abs() < 0.1,
        "Expected theta1 still in dead zone, got {theta1}"
    );
    assert!(
        theta2.abs() < 1e-10,
        "D should not rotate inside dead zone, theta2 = {theta2}"
    );

    // Continue until theta1 > 0.1 (crosses dead-zone boundary)
    let max_extra = 10000usize;
    let mut crossed = false;
    for _ in 0..max_extra {
        gc.solve_iteration(
            &mut BodyState::new(&mut omega_a, &inv_ia),
            &mut BodyState::new(&mut omega_b, &inv_ib),
            &mut BodyState::new(&mut omega_c, &inv_ic),
            &mut BodyState::new(&mut omega_d, &inv_id),
            theta1,
            theta2,
            dt,
        );
        theta1 += omega_b[2] * dt;
        theta2 += omega_d[2] * dt;
        if theta1 > 0.1 {
            crossed = true;
            break;
        }
    }

    assert!(crossed, "theta1 never exceeded 0.1 in {max_extra} iterations");

    // After crossing, run a few more iterations and confirm D is moving
    for _ in 0..100 {
        gc.solve_iteration(
            &mut BodyState::new(&mut omega_a, &inv_ia),
            &mut BodyState::new(&mut omega_b, &inv_ib),
            &mut BodyState::new(&mut omega_c, &inv_ic),
            &mut BodyState::new(&mut omega_d, &inv_id),
            theta1,
            theta2,
            dt,
        );
        theta1 += omega_b[2] * dt;
        theta2 += omega_d[2] * dt;
    }

    assert!(
        theta2 > 0.0,
        "D should start moving after backlash is overcome, theta2 = {theta2}"
    );
}

// ---------------------------------------------------------------------------
// Test 3 — ratio=1 identity
// ---------------------------------------------------------------------------

/// With ratio=1 and no backlash, B and D must converge to equal angular speed.
#[test]
fn ratio_one_identity() {
    let gc = GearConstraint::new(1.0, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);

    let inv_ia = [0.0f64; 3];
    let inv_ic = [0.0f64; 3];
    let inv_ib = [1.0f64; 3];
    let inv_id = [1.0f64; 3];

    let mut omega_a = [0.0f64; 3];
    let mut omega_b = [0.0f64, 0.0, 2.0];
    let mut omega_c = [0.0f64; 3];
    let mut omega_d = [0.0f64, 0.0, 0.5];

    let dt = 1.0 / 60.0;
    let mut theta1 = 0.0f64;
    let mut theta2 = 0.0f64;

    for _ in 0..500 {
        gc.solve_iteration(
            &mut BodyState::new(&mut omega_a, &inv_ia),
            &mut BodyState::new(&mut omega_b, &inv_ib),
            &mut BodyState::new(&mut omega_c, &inv_ic),
            &mut BodyState::new(&mut omega_d, &inv_id),
            theta1,
            theta2,
            dt,
        );
        theta1 += omega_b[2] * dt;
        theta2 += omega_d[2] * dt;
    }

    let vel_err = (omega_b[2] - omega_d[2]).abs();
    assert!(
        vel_err < 1e-6,
        "velocity should match at ratio=1: |ω_Bz - ω_Dz| = {vel_err}"
    );
}

// ---------------------------------------------------------------------------
// Test 4 — no energy injection
// ---------------------------------------------------------------------------

/// The constraint must never add energy to the system (it can only remove it).
#[test]
fn no_energy_injection() {
    let gc = GearConstraint::new(2.0, 0.0, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]);

    let inv_ia = [0.0f64; 3];
    let inv_ic = [0.0f64; 3];
    let inv_ib = [1.0f64; 3]; // I=1 → KE = 0.5·ω²
    let inv_id = [1.0f64; 3];

    let mut omega_a = [0.0f64; 3];
    let mut omega_b = [0.0f64, 0.0, 1.0];
    let mut omega_c = [0.0f64; 3];
    let mut omega_d = [0.0f64, 0.0, 0.3];

    let ke_initial = 0.5 * omega_b[2] * omega_b[2] + 0.5 * omega_d[2] * omega_d[2];

    let dt = 1.0 / 60.0;

    for _ in 0..500 {
        gc.solve_iteration(
            &mut BodyState::new(&mut omega_a, &inv_ia),
            &mut BodyState::new(&mut omega_b, &inv_ib),
            &mut BodyState::new(&mut omega_c, &inv_ic),
            &mut BodyState::new(&mut omega_d, &inv_id),
            0.0,
            0.0, // no angle integration needed for this test
            dt,
        );
    }

    let ke_final = 0.5 * omega_b[2] * omega_b[2] + 0.5 * omega_d[2] * omega_d[2];

    assert!(
        ke_final <= ke_initial + 1e-10,
        "Constraint injected energy: KE_initial={ke_initial}, KE_final={ke_final}"
    );
}
