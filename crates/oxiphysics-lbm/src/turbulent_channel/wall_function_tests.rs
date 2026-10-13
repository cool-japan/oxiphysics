// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the log-law wall function of [`WallFunctionLbm`]:
//! the modelled wall shear stress must act on the populations (it used to be
//! written into a diagnostic velocity that the next collision overwrote), the
//! lattice walls must be free-slip so the model is the only wall drag, and the
//! steady state must satisfy the law of the wall.

use super::functions::*;
use super::types::*;

/// Deterministic pseudo-random numbers in [0, 1) (64-bit LCG).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Set every fluid node to the equilibrium of density `rho` and velocity `(ux, 0)`.
fn set_uniform_flow(wf: &mut WallFunctionLbm, rho: f64, ux: f64) {
    let ny = wf.ny;
    for i in 0..wf.nx {
        for j in 1..ny - 1 {
            for q in 0..9 {
                wf.f_dist[(i * ny + j) * 9 + q] = equilibrium_d2q9(rho, ux, 0.0, q);
            }
            wf.u[i * ny + j] = [ux, 0.0];
        }
    }
}

/// Total streamwise momentum and mass of the fluid rows.
fn momentum_and_mass(wf: &WallFunctionLbm) -> (f64, f64) {
    let mut momentum = 0.0;
    let mut mass = 0.0;
    for i in 0..wf.nx {
        for j in 1..wf.ny - 1 {
            let base = (i * wf.ny + j) * 9;
            for (q, &cx) in CX.iter().enumerate() {
                momentum += cx as f64 * wf.f_dist[base + q];
                mass += wf.f_dist[base + q];
            }
        }
    }
    (momentum, mass)
}

fn max_speed(wf: &WallFunctionLbm) -> f64 {
    let ny = wf.ny;
    let mut worst = 0.0_f64;
    for i in 0..wf.nx {
        for j in 1..ny - 1 {
            let [ux, uy] = wf.u[i * ny + j];
            worst = worst.max(ux.hypot(uy));
        }
    }
    worst
}

/// x-averaged streamwise velocity of row `j`.
fn row_velocity(wf: &WallFunctionLbm, j: usize) -> f64 {
    let (nx, ny) = (wf.nx, wf.ny);
    (0..nx).map(|i| wf.u[i * ny + j][0]).sum::<f64>() / nx as f64
}

// ─── Law of the wall and its inversion ──────────────────────────────────────

#[test]
fn law_of_the_wall_follows_sublayer_and_log_branches() {
    assert_eq!(law_of_the_wall_u_plus(0.0), 0.0);
    assert_eq!(law_of_the_wall_u_plus(-3.0), 0.0);
    assert_eq!(law_of_the_wall_u_plus(5.0), 5.0);
    let expected = 100.0_f64.ln() / KAPPA_VK + B_LL;
    assert!((law_of_the_wall_u_plus(100.0) - expected).abs() < 1e-14);
    let mut previous = 0.0;
    for k in 13..=2000 {
        let u_plus = law_of_the_wall_u_plus(k as f64);
        assert!(
            u_plus > previous,
            "u+ must grow with y+ above the transition"
        );
        previous = u_plus;
    }
}

#[test]
fn wall_friction_velocity_inverts_the_law_in_every_regime() {
    let nu = 1.0e-4;
    let u_tau = 0.01;
    for y_plus in [0.5, 3.0, 8.0, 11.0, 12.5, 20.0, 30.0, 100.0, 300.0, 1000.0] {
        let y = y_plus * nu / u_tau;
        let speed = u_tau * law_of_the_wall_u_plus(y_plus);
        for sign in [1.0, -1.0] {
            let recovered = wall_friction_velocity(sign * speed, y, nu);
            assert!(
                (recovered - u_tau).abs() < 1e-10 * u_tau,
                "y+ = {y_plus}: recovered u_tau {recovered}, expected {u_tau}"
            );
        }
    }
}

#[test]
fn wall_friction_velocity_is_monotone_in_the_velocity() {
    let (y, nu) = (0.5, 1.6e-4);
    let mut previous = 0.0;
    let mut speed = 1.0e-6;
    while speed < 2.0 {
        let u_tau = wall_friction_velocity(speed, y, nu);
        assert!(u_tau >= previous, "u_tau fell at |u| = {speed}");
        // Whichever root is returned reproduces the input velocity.
        let implied = u_tau * law_of_the_wall_u_plus(y * u_tau / nu);
        assert!((implied - speed).abs() < 1e-9 * speed, "|u| = {speed}");
        previous = u_tau;
        speed *= 1.07;
    }
}

#[test]
fn wall_friction_velocity_rejects_degenerate_input() {
    assert_eq!(wall_friction_velocity(0.0, 0.5, 1e-4), 0.0);
    assert_eq!(wall_friction_velocity(0.1, 0.0, 1e-4), 0.0);
    assert_eq!(wall_friction_velocity(0.1, -1.0, 1e-4), 0.0);
    assert_eq!(wall_friction_velocity(0.1, 0.5, 0.0), 0.0);
    assert_eq!(wall_friction_velocity(f64::NAN, 0.5, 1e-4), 0.0);
    assert_eq!(wall_friction_velocity(f64::INFINITY, 0.5, 1e-4), 0.0);
}

// ─── Geometry ──────────────────────────────────────────────────────────────

#[test]
fn geometry_forcing_and_viscosity_are_consistent() {
    let wf = WallFunctionLbm::new(4, 10, 240.0);
    // Eight fluid rows: half-width 4, first row centre half a cell off the wall.
    assert_eq!(wf.params.half_width, 4.0);
    assert_eq!(wf.y_first, 0.5);
    assert!((wf.params.nu - 0.01 * 4.0 / 240.0).abs() < 1e-18);
    assert!((wf.y_plus_first - 30.0).abs() < 1e-9, "{}", wf.y_plus_first);
    assert!((wf.y_plus_at(1) - wf.y_plus_first).abs() < 1e-12);
    assert!((wf.y_plus_at(8) - wf.y_plus_first).abs() < 1e-12);
    assert!((wf.y_plus_at(2) - 3.0 * wf.y_plus_first).abs() < 1e-9);
    // Momentum balance of the discrete channel: body force x half-width = u_tau^2.
    assert!((wf.body_force * wf.params.half_width - wf.params.u_tau.powi(2)).abs() < 1e-18);
    assert!((wf.friction_velocity() - wf.params.u_tau).abs() < 1e-15);
    // Distances to the wall.
    assert_eq!(wf.wall_distance(0), 0.0);
    assert_eq!(wf.wall_distance(9), 0.0);
    assert_eq!(wf.wall_distance(1), 0.5);
    assert_eq!(wf.wall_distance(8), 0.5);
    assert_eq!(wf.wall_distance(2), 1.5);
    assert_eq!(wf.wall_distance(4), 3.5);
    assert_eq!(wf.wall_distance(5), 3.5);
    assert_eq!(wf.y_plus_at(0), 0.0);
}

#[test]
fn constructor_with_friction_velocity_scales_viscosity_and_force() {
    let slow = WallFunctionLbm::with_friction_velocity(4, 8, 180.0, 0.004);
    assert_eq!(slow.params.u_tau, 0.004);
    assert!((slow.params.nu - 0.004 * 3.0 / 180.0).abs() < 1e-18);
    assert!((slow.body_force - 0.004_f64.powi(2) / 3.0).abs() < 1e-18);
    // `new` is the u_tau = 0.01 special case.
    let standard = WallFunctionLbm::new(4, 8, 180.0);
    assert_eq!(standard.params.u_tau, 0.01);
}

// ─── Free-slip walls ───────────────────────────────────────────────────────

#[test]
fn streaming_reflects_populations_specularly_at_both_walls() {
    let (nx, ny) = (5, 7);
    // (row, source direction, expected direction after streaming, column shift)
    let cases = [
        (1, 4, 2, 0),  // (0,-1) -> (0,+1) at the lower wall
        (1, 7, 6, -1), // (-1,-1) -> (-1,+1)
        (1, 8, 5, 1),  // (1,-1) -> (1,+1)
        (5, 2, 4, 0),  // (0,+1) -> (0,-1) at the upper wall
        (5, 5, 8, 1),  // (1,1) -> (1,-1)
        (5, 6, 7, -1), // (-1,1) -> (-1,-1)
        (1, 2, 2, 0),  // free flight away from the wall: lands one row up
        (3, 1, 1, 1),  // purely streamwise population in the bulk
    ];
    for (row, q, expected_q, shift) in cases {
        let mut wf = WallFunctionLbm::new(nx, ny, 180.0);
        wf.f_dist.iter_mut().for_each(|v| *v = 0.0);
        let src_col = 2usize;
        wf.f_dist[(src_col * ny + row) * 9 + q] = 1.0;
        wf.stream();
        let dst_col = (src_col as i64 + shift).rem_euclid(nx as i64) as usize;
        // Only the population flying away from the wall changes its row.
        let dst_row = if (row, q) == (1, 2) { row + 1 } else { row };
        let hit = wf.f_dist[(dst_col * ny + dst_row) * 9 + expected_q];
        assert_eq!(hit, 1.0, "row {row}, direction {q}");
        let total: f64 = wf.f_dist.iter().sum();
        assert_eq!(
            total, 1.0,
            "row {row}, direction {q}: population lost or duplicated"
        );
    }
}

#[test]
fn streaming_conserves_mass_and_streamwise_momentum_and_spares_solid_rows() {
    let (nx, ny) = (6, 9);
    let mut wf = WallFunctionLbm::new(nx, ny, 180.0);
    let mut rng = Lcg(11);
    for i in 0..nx {
        for j in 1..ny - 1 {
            for q in 0..9 {
                wf.f_dist[(i * ny + j) * 9 + q] = rng.next();
            }
        }
    }
    let solid_before: Vec<f64> = (0..nx)
        .flat_map(|i| {
            [0, ny - 1]
                .into_iter()
                .flat_map(move |j| (0..9).map(move |q| (i * ny + j) * 9 + q))
        })
        .map(|k| wf.f_dist[k])
        .collect();
    let (p0, m0) = momentum_and_mass(&wf);
    for _ in 0..25 {
        wf.stream();
    }
    let (p1, m1) = momentum_and_mass(&wf);
    assert!((m1 - m0).abs() < 1e-12 * m0, "mass {m0} -> {m1}");
    assert!((p1 - p0).abs() < 1e-12 * m0, "x momentum {p0} -> {p1}");
    let solid_after: Vec<f64> = (0..nx)
        .flat_map(|i| {
            [0, ny - 1]
                .into_iter()
                .flat_map(move |j| (0..9).map(move |q| (i * ny + j) * 9 + q))
        })
        .map(|k| wf.f_dist[k])
        .collect();
    assert_eq!(solid_before, solid_after, "solid rows must not be touched");
}

#[test]
fn free_slip_walls_without_the_wall_model_exert_no_drag() {
    // collide + stream only: with the stored wall force at zero, neither the
    // lattice walls nor the interior may remove streamwise momentum.
    let mut wf = WallFunctionLbm::new(4, 10, 240.0);
    wf.body_force = 0.0;
    set_uniform_flow(&mut wf, 1.0, 0.05);
    let (p0, m0) = momentum_and_mass(&wf);
    for _ in 0..200 {
        wf.collide();
        wf.stream();
    }
    let (p1, m1) = momentum_and_mass(&wf);
    assert!((m1 - m0).abs() < 1e-12 * m0);
    assert!(
        (p1 - p0).abs() < 1e-12 * p0,
        "momentum {p0} -> {p1}: the lattice walls must be free-slip"
    );
    assert!(wf.wall_force.iter().all(|&f| f == 0.0));
}

// ─── The wall model acts on the dynamics ───────────────────────────────────

#[test]
fn wall_force_follows_the_law_of_the_wall_and_opposes_the_flow() {
    let mut wf = WallFunctionLbm::new(3, 10, 240.0);
    let nu = wf.params.nu;
    let rho = 1.02;
    for ux in [0.1, -0.1, 0.002] {
        set_uniform_flow(&mut wf, rho, ux);
        wf.apply_wall_function();
        assert_eq!(wf.wall_force.len(), 2 * wf.nx);
        assert_eq!(wf.wall_u_tau.len(), 2 * wf.nx);
        for k in 0..2 * wf.nx {
            let u_tau = wf.wall_u_tau[k];
            // The friction velocity reproduces the matching-point velocity ...
            let implied = u_tau * law_of_the_wall_u_plus(wf.y_first * u_tau / nu);
            assert!((implied - ux.abs()).abs() < 1e-12, "ux = {ux}: {implied}");
            // ... and the force is the wall shear stress, against the flow.
            let expected = -ux.signum() * rho * u_tau * u_tau;
            assert!(
                (wf.wall_force[k] - expected).abs() < 1e-15,
                "ux = {ux}, entry {k}"
            );
            assert!(wf.wall_force[k] * ux < 0.0);
        }
        assert!((wf.model_friction_velocity() - wf.wall_u_tau[0]).abs() < 1e-15);
    }
    // No flow, no stress.
    set_uniform_flow(&mut wf, 1.0, 0.0);
    wf.apply_wall_function();
    assert!(wf.wall_force.iter().all(|&f| f == 0.0));
    assert!(wf.wall_u_tau.iter().all(|&u| u == 0.0));
}

#[test]
fn apply_wall_function_records_the_force_and_leaves_the_fields_alone() {
    // The wall model used to overwrite the diagnostic velocity of the first
    // fluid rows; it must now only record the force for the next collision.
    let mut wf = WallFunctionLbm::new(3, 10, 240.0);
    set_uniform_flow(&mut wf, 1.0, 0.1);
    let velocity_before = wf.u.clone();
    let populations_before = wf.f_dist.clone();
    wf.apply_wall_function();
    assert_eq!(wf.u, velocity_before);
    assert_eq!(wf.f_dist, populations_before);
    assert!(wf.wall_force.iter().all(|&f| f < 0.0));
}

#[test]
fn a_step_removes_exactly_the_modelled_wall_shear_momentum() {
    for ny in [3, 4, 10] {
        let mut wf = WallFunctionLbm::new(4, ny, 240.0);
        wf.body_force = 0.0;
        set_uniform_flow(&mut wf, 1.0, 0.1);
        wf.apply_wall_function();
        let u_tau = wf.model_friction_velocity();
        assert!(u_tau > 0.0);
        let (p0, m0) = momentum_and_mass(&wf);
        wf.collide();
        wf.stream();
        let (p1, m1) = momentum_and_mass(&wf);
        // Two walls, one lattice time unit, unit density, columns nx.
        let expected = -2.0 * wf.nx as f64 * u_tau * u_tau;
        assert!(
            ((p1 - p0) - expected).abs() < 1e-13,
            "ny = {ny}: momentum change {} vs {expected}",
            p1 - p0
        );
        assert!((m1 - m0).abs() < 1e-12);
    }
}

#[test]
fn the_wall_shear_force_is_the_only_momentum_sink() {
    // No driving force and no eddy-viscosity momentum loss: over many steps
    // the streamwise momentum falls by exactly the wall forces that the
    // successive collisions applied.
    let mut wf = WallFunctionLbm::new(4, 12, 240.0);
    wf.body_force = 0.0;
    set_uniform_flow(&mut wf, 1.0, 0.12);
    wf.apply_wall_function();
    let (p_start, m_start) = momentum_and_mass(&wf);
    let mut applied = 0.0;
    for _ in 0..50 {
        applied += wf.wall_force.iter().sum::<f64>();
        wf.step();
    }
    let (p_end, m_end) = momentum_and_mass(&wf);
    assert!(
        applied < -1e-3,
        "the walls must decelerate the flow: {applied}"
    );
    assert!(
        (p_end - p_start - applied).abs() < 1e-11,
        "momentum change {} vs applied wall force {applied}",
        p_end - p_start
    );
    assert!((m_end - m_start).abs() < 1e-11 * m_start);
    assert!(p_end / m_end < 0.12, "the mean velocity must have dropped");
}

// ─── Eddy viscosity ────────────────────────────────────────────────────────

#[test]
fn eddy_viscosity_vanishes_for_uniform_flow_and_follows_the_mixing_length() {
    let mut wf = WallFunctionLbm::new(3, 12, 240.0);
    set_uniform_flow(&mut wf, 1.0, 0.08);
    wf.update_eddy_viscosity();
    for j in 2..wf.ny - 2 {
        for i in 0..wf.nx {
            assert!(wf.nu_t[i * wf.ny + j].abs() < 1e-15, "row {j}");
        }
    }
    // Linear profile u = 0.01 j: du/dy = 0.01 between fluid rows.
    for i in 0..wf.nx {
        for j in 1..wf.ny - 1 {
            for q in 0..9 {
                wf.f_dist[(i * wf.ny + j) * 9 + q] = equilibrium_d2q9(1.0, 0.01 * j as f64, 0.0, q);
            }
        }
    }
    wf.update_eddy_viscosity();
    for j in 2..wf.ny - 2 {
        let l_m = wf.mixing_length(j);
        let expected = l_m * l_m * 0.01;
        assert!((wf.nu_t[j] - expected).abs() < 1e-14, "row {j}");
        assert!(wf.nu_t[j] > 0.0);
    }
    // Mixing length: kappa d near the wall, capped at 0.09 h away from it.
    assert!((wf.mixing_length(1) - KAPPA_VK * 0.5).abs() < 1e-15);
    assert!((wf.mixing_length(5) - 0.09 * wf.params.half_width).abs() < 1e-15);
    assert_eq!(wf.mixing_length(0), 0.0);
}

#[test]
fn eddy_viscosity_next_to_the_wall_uses_the_law_of_the_wall_gradient() {
    let mut wf = WallFunctionLbm::new(2, 12, 240.0);
    set_uniform_flow(&mut wf, 1.0, 0.1);
    wf.apply_wall_function();
    wf.update_eddy_viscosity();
    let u_tau = wf.wall_u_tau[0];
    let d = wf.wall_distance(1);
    let y_plus = d * u_tau / wf.params.nu;
    assert!(
        y_plus > Y_PLUS_TRANSITION,
        "set-up must be in the log layer: {y_plus}"
    );
    let l_m = wf.mixing_length(1);
    let expected = l_m * l_m * u_tau / (KAPPA_VK * d);
    assert!((wf.nu_t[1] - expected).abs() < 1e-15);
    assert!((wf.nu_t[wf.ny - 2] - expected).abs() < 1e-15);
    // In the viscous sublayer the stress is molecular.
    set_uniform_flow(&mut wf, 1.0, 1.0e-4);
    wf.apply_wall_function();
    wf.update_eddy_viscosity();
    assert!(wf.wall_u_tau[0] * wf.wall_distance(1) / wf.params.nu < Y_PLUS_TRANSITION);
    assert_eq!(wf.nu_t[1], 0.0);
}

// ─── Steady state ──────────────────────────────────────────────────────────

/// Run to (near) steady state and return the largest speed seen on the way.
fn run_and_track_speed(wf: &mut WallFunctionLbm, steps: usize) -> f64 {
    let mut worst = 0.0_f64;
    for step in 0..steps {
        wf.step();
        if step.is_multiple_of(50) {
            worst = worst.max(max_speed(wf));
        }
    }
    worst
}

#[test]
fn channel_at_y_plus_30_converges_to_the_log_law() {
    // Re_tau = 180 with six fluid rows: half-width 3, first row at y+ = 30.
    let mut wf = WallFunctionLbm::new(4, 8, 180.0);
    assert!((wf.y_plus_first - 30.0).abs() < 1e-9);
    let (_, mass0) = momentum_and_mass(&wf);
    let worst = run_and_track_speed(&mut wf, 40_000);
    let (_, mass1) = momentum_and_mass(&wf);

    let design_u_tau = wf.params.u_tau;
    let measured_u_tau = wf.model_friction_velocity();
    assert!(
        (measured_u_tau - design_u_tau).abs() < 0.02 * design_u_tau,
        "the wall model must carry the driving force: u_tau {measured_u_tau} vs {design_u_tau}"
    );

    // u(j = 1) / u_tau reproduces the log-law value at y+ = 30.
    let u_plus_law = law_of_the_wall_u_plus(wf.y_plus_first);
    let u_plus_measured = row_velocity(&wf, 1) / design_u_tau;
    assert!(
        (u_plus_measured - u_plus_law).abs() < 0.05 * u_plus_law,
        "u+ {u_plus_measured} vs log-law {u_plus_law}"
    );

    // The velocity profile is symmetric, monotone and stays low-Mach.
    for j in 1..wf.ny - 1 {
        let mirror = wf.ny - 1 - j;
        assert!((row_velocity(&wf, j) - row_velocity(&wf, mirror)).abs() < 1e-9);
    }
    // Rows ny/2 - 1 and ny/2 are mirror images across the centreline.
    for j in 1..wf.ny / 2 - 1 {
        assert!(row_velocity(&wf, j + 1) > row_velocity(&wf, j), "row {j}");
    }
    assert!(worst < 0.3, "peak speed {worst}");
    let centre = row_velocity(&wf, wf.ny / 2) / design_u_tau;
    assert!(
        centre > u_plus_law && centre < 1.5 * u_plus_law,
        "centreline u+ {centre} (log-law value at the first row: {u_plus_law})"
    );
    assert!(
        (mass1 - mass0).abs() < 1e-10 * mass0,
        "mass {mass0} -> {mass1}"
    );
}

#[test]
fn channel_converges_to_the_same_state_from_above_and_below() {
    let steps = 40_000;
    let mut from_rest = WallFunctionLbm::new(4, 8, 180.0);
    run_and_track_speed(&mut from_rest, steps);
    let mut from_above = WallFunctionLbm::new(4, 8, 180.0);
    set_uniform_flow(&mut from_above, 1.0, 0.25);
    let worst = run_and_track_speed(&mut from_above, steps);
    assert!(worst < 0.3, "peak speed {worst}");
    for j in 1..from_rest.ny - 1 {
        let a = row_velocity(&from_rest, j);
        let b = row_velocity(&from_above, j);
        assert!((a - b).abs() < 0.01 * a, "row {j}: {a} vs {b}");
    }
    let u_tau = from_above.params.u_tau;
    assert!((from_above.model_friction_velocity() - u_tau).abs() < 0.02 * u_tau);
}

#[test]
fn smaller_friction_velocity_gives_a_lower_mach_number_at_the_same_y_plus() {
    // Same y+ = 30 at the first row (Re_tau / h fixed), a third of the speed.
    let mut slow = WallFunctionLbm::with_friction_velocity(4, 8, 180.0, 0.0033);
    assert!((slow.y_plus_first - 30.0).abs() < 1e-9);
    let worst = run_and_track_speed(&mut slow, 60_000);
    assert!(worst < 0.1, "peak speed {worst}");
    let u_tau = slow.params.u_tau;
    assert!((slow.model_friction_velocity() - u_tau).abs() < 0.02 * u_tau);
    let u_plus = row_velocity(&slow, 1) / u_tau;
    let law = law_of_the_wall_u_plus(slow.y_plus_first);
    assert!((u_plus - law).abs() < 0.05 * law, "u+ {u_plus} vs {law}");
}

#[test]
fn solid_rows_stay_untouched_and_degenerate_sizes_do_not_panic() {
    let mut wf = WallFunctionLbm::new(4, 8, 180.0);
    for _ in 0..200 {
        wf.step();
    }
    for i in 0..wf.nx {
        for j in [0, wf.ny - 1] {
            for (q, &weight) in W9.iter().enumerate() {
                assert_eq!(wf.f_dist[(i * wf.ny + j) * 9 + q], weight);
            }
            assert_eq!(wf.u[i * wf.ny + j], [0.0, 0.0]);
            assert_eq!(wf.nu_t[i * wf.ny + j], 0.0);
        }
    }
    // No fluid rows at all.
    for ny in [0, 1, 2] {
        let mut tiny = WallFunctionLbm::new(3, ny, 180.0);
        for _ in 0..5 {
            tiny.step();
        }
        assert_eq!(tiny.bulk_velocity(), 0.0);
        assert!(tiny.f_dist.iter().all(|v| v.is_finite()));
    }
    // A single fluid row between the two walls.
    let mut single = WallFunctionLbm::new(3, 3, 180.0);
    for _ in 0..500 {
        single.step();
    }
    assert!(single.bulk_velocity() > 0.0);
    assert!(single.f_dist.iter().all(|v| v.is_finite()));
}
