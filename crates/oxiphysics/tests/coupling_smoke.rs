// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoke and regression tests for the OxiCAR domain coupling layer
//! (Blueprint KF-1).
//!
//! The first six tests are the original smoke tests:
//! 1. Empty runtime — step produces an empty report list.
//! 2. Copy semantics — `InterfaceSite`, `InterfaceState`, `InterfaceForce`,
//!    and `CouplingReport` are all `Copy`.
//! 3. Zero-velocity exchange — identical velocities on both sides produce
//!    zero momentum exchange when `stiffness > 0, damping > 0`.
//! 4. Nonzero velocity exchange — confirms the penalty formula is correct.
//! 5. MD–continuum smoke — `MdContinuumAdapter` links two mock domains;
//!    action–reaction is satisfied.
//! 6. Two-pair schedule — runtime with two coupler pairs processes both
//!    pairs in one `step` call and returns two reports.
//!
//! The remaining tests are regressions for the penalty law and the runtime:
//! * the force is `k (x_a - x_b) + c (v_a - v_b)` (stiffness on position,
//!   damping on velocity, in newtons) and pulls side A towards side B,
//! * the runtime samples both domains at the *coupler's* interface sites, so
//!   forces are actually exchanged,
//! * the mock domains integrate the forces with the step length.

use oxiphysics::coupling::fem_sph::{ConstantVelocityDomain, FemSphCoupler, ImpulseDomain};
use oxiphysics::coupling::md_continuum_adapter::{
    MdContinuumAdapter, MockContinuumDomain, MockMdDomain,
};
use oxiphysics::coupling::{
    CouplingDomain, CouplingReport, CouplingRuntime, DomainCoupler, DomainKind, InterfaceForce,
    InterfaceForceVec, InterfaceSite, InterfaceState, InterfaceStateVec,
};

// ─────────────────────────────────────────────────────────────────────────────
// Test 1: empty runtime produces empty reports
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn empty_runtime_produces_no_reports() {
    let mut rt = CouplingRuntime::new();
    let reports = rt.step(0.01);
    assert!(
        reports.is_empty(),
        "expected no reports for empty schedule, got {}",
        reports.len()
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 2: payload types are Copy
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn payload_types_are_copy() {
    // If these compile, the types implement Copy.
    let site = InterfaceSite {
        id: 0,
        position: [1.0, 2.0, 3.0],
    };
    let site2 = site; // copy, not move
    let _ = site; // original still usable
    let _ = site2;

    let state = InterfaceState {
        id: 0,
        position: [0.0; 3],
        velocity: [1.0, 0.0, 0.0],
    };
    let state2 = state;
    let _ = state;
    let _ = state2;

    let force = InterfaceForce {
        id: 0,
        force: [0.0, -9.81, 0.0],
    };
    let force2 = force;
    let _ = force;
    let _ = force2;

    let report = CouplingReport {
        site_count: 3,
        momentum_a_to_b: [0.0; 3],
        momentum_b_to_a: [0.0; 3],
        rms_position_error: 0.0,
    };
    let report2 = report;
    let _ = report;
    let _ = report2;
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 3: zero velocity / position mismatch → zero momentum exchange
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn zero_velocity_mismatch_produces_zero_momentum() {
    // Both domains move at the same velocity; positions are identical.
    // stiffness > 0, damping > 0 — but all differences are zero.
    let vel = [1.0_f64, 2.0, 3.0];
    let sites = vec![
        InterfaceSite {
            id: 0,
            position: [0.0, 0.0, 0.0],
        },
        InterfaceSite {
            id: 1,
            position: [1.0, 0.0, 0.0],
        },
    ];

    let mut coupler = FemSphCoupler::new(sites.clone(), 1000.0, 50.0);

    // Build matching state vecs manually — same pos, same vel on both sides.
    let make_states = |s: &[InterfaceSite]| InterfaceStateVec {
        states: s
            .iter()
            .map(|site| InterfaceState {
                id: site.id,
                position: site.position,
                velocity: vel,
            })
            .collect(),
    };

    let state_a = make_states(&sites);
    let state_b = make_states(&sites);

    let (_fa, _fb, report) = coupler.step(0.01, &state_a, &state_b);

    assert_eq!(report.site_count, 2);
    for j in 0..3 {
        assert!(
            report.momentum_a_to_b[j].abs() < 1e-15,
            "expected zero momentum_a_to_b[{j}], got {}",
            report.momentum_a_to_b[j]
        );
    }
    assert!(
        report.rms_position_error < 1e-15,
        "expected zero rms_position_error, got {}",
        report.rms_position_error
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 4: nonzero velocity mismatch → expected momentum exchange
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn nonzero_velocity_mismatch_produces_correct_momentum() {
    // One site: FEM moves at [1,0,0], SPH at [0,0,0]; positions coincide.
    // f = k*(x_fem - x_sph) + c*(v_fem - v_sph)
    //   = 0 + c*[1,0,0]
    let c = 500.0_f64;
    let dt = 0.02_f64;
    let expected_f_x = c * 1.0; // newtons, along X
    let expected_momentum_x = expected_f_x * dt; // f*dt

    let sites = vec![InterfaceSite {
        id: 0,
        position: [0.0, 0.0, 0.0],
    }];
    let mut coupler = FemSphCoupler::new(sites, 0.0, 0.0);
    coupler.stiffness = 0.0;
    coupler.damping = c;

    let state_a = InterfaceStateVec {
        states: vec![InterfaceState {
            id: 0,
            position: [0.0, 0.0, 0.0],
            velocity: [1.0, 0.0, 0.0],
        }],
    };
    let state_b = InterfaceStateVec {
        states: vec![InterfaceState {
            id: 0,
            position: [0.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
        }],
    };

    let (fa, fb, report) = coupler.step(dt, &state_a, &state_b);

    assert_eq!(report.site_count, 1);
    let tol = 1e-12;
    // The FEM side moves faster, so damping drags it back (-f) and the SPH side
    // along (+f).
    assert!(
        (fa.forces[0].force[0] + expected_f_x).abs() < tol,
        "FEM force[0]: expected {}, got {}",
        -expected_f_x,
        fa.forces[0].force[0]
    );
    assert!(
        (fb.forces[0].force[0] - expected_f_x).abs() < tol,
        "SPH force[0]: expected {expected_f_x}, got {}",
        fb.forces[0].force[0]
    );
    assert!(
        (report.momentum_a_to_b[0] - expected_momentum_x).abs() < tol,
        "momentum_a_to_b[0]: expected {expected_momentum_x}, got {}",
        report.momentum_a_to_b[0]
    );
    // Action–reaction: momentum_b_to_a should be −momentum_a_to_b
    assert!(
        (report.momentum_b_to_a[0] + expected_momentum_x).abs() < tol,
        "momentum_b_to_a[0]: expected {}, got {}",
        -expected_momentum_x,
        report.momentum_b_to_a[0]
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 5: MD–continuum smoke — action–reaction satisfied
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn md_continuum_adapter_action_reaction() {
    // MD particles at [1,0,0] with velocity [2,0,0].
    // Continuum nodes at [0,0,0] with velocity [0,0,0].
    // stiffness=100, damping=0 → f = k*(x_md - x_cont) = [100,0,0] N, which
    // pulls the MD side (-f) towards the continuum and the continuum (+f)
    // towards the MD side.
    let stiffness = 100.0_f64;
    let damping = 0.0_f64;
    let dt = 0.01_f64;

    let state_a = InterfaceStateVec {
        states: vec![InterfaceState {
            id: 0,
            position: [1.0, 0.0, 0.0],
            velocity: [2.0, 0.0, 0.0],
        }],
    };
    let state_b = InterfaceStateVec {
        states: vec![InterfaceState {
            id: 0,
            position: [0.0, 0.0, 0.0],
            velocity: [0.0, 0.0, 0.0],
        }],
    };

    let sites = vec![InterfaceSite {
        id: 0,
        position: [0.5, 0.0, 0.0],
    }];
    let mut adapter = MdContinuumAdapter::new(sites, stiffness, damping);
    assert_eq!(adapter.interface_sites().len(), 1);
    let (fa, fb, report) = adapter.step(dt, &state_a, &state_b);

    // f = 100 * [1,0,0] = [100,0,0]
    let expected_fx = 100.0_f64;
    let tol = 1e-10;

    assert_eq!(report.site_count, 1);
    assert!(
        (fa.forces[0].force[0] + expected_fx).abs() < tol,
        "MD force[0]: expected {}, got {}",
        -expected_fx,
        fa.forces[0].force[0]
    );
    assert!(
        (fb.forces[0].force[0] - expected_fx).abs() < tol,
        "Continuum force[0]: expected {expected_fx}, got {}",
        fb.forces[0].force[0]
    );

    // Total impulse (a_to_b + b_to_a) should sum to zero.
    let total = report.momentum_a_to_b[0] + report.momentum_b_to_a[0];
    assert!(
        total.abs() < tol,
        "total impulse should be zero, got {total}"
    );

    // Mock domain round-trip: create domains and run through runtime.
    let mut md =
        MockMdDomain::from_particles(vec![[1.0, 0.0, 0.0]], vec![[2.0, 0.0, 0.0]], vec![1.0]);
    let mut cont =
        MockContinuumDomain::from_nodes(vec![[0.0, 0.0, 0.0]], vec![[0.0, 0.0, 0.0]], vec![1.0]);

    // Verify sampling works.
    let sites_md = [InterfaceSite {
        id: 0,
        position: [1.0, 0.0, 0.0],
    }];
    let sampled_md = md.sample_interface(&sites_md);
    assert_eq!(sampled_md.states.len(), 1);
    assert!((sampled_md.states[0].velocity[0] - 2.0).abs() < tol);

    let sites_cont = [InterfaceSite {
        id: 0,
        position: [0.0, 0.0, 0.0],
    }];
    let sampled_cont = cont.sample_interface(&sites_cont);
    assert_eq!(sampled_cont.states.len(), 1);
    assert!(sampled_cont.states[0].velocity[0].abs() < tol);

    // apply_interface_force: a force of [10,0,0] N on a 1 kg particle for
    // dt = 0.5 s is an impulse of 5 N·s: v = 2 + 10*0.5/1 = 7, and the position
    // advances with the new velocity: x = 1 + 7*0.5 = 4.5.
    let kick = InterfaceForceVec {
        forces: vec![InterfaceForce {
            id: 0,
            force: [10.0, 0.0, 0.0],
        }],
    };
    md.apply_interface_force(&kick, 0.5);
    assert!(
        (md.velocities[0][0] - 7.0).abs() < tol,
        "MD velocity after kick: expected 7, got {}",
        md.velocities[0][0]
    );
    assert!(
        (md.positions[0][0] - 4.5).abs() < tol,
        "MD position after kick: expected 4.5, got {}",
        md.positions[0][0]
    );

    // Continuum gets the negative kick: v = 0 - 10*0.5/1 = -5, x = 0 - 5*0.5.
    let neg_kick = InterfaceForceVec {
        forces: vec![InterfaceForce {
            id: 0,
            force: [-10.0, 0.0, 0.0],
        }],
    };
    cont.apply_interface_force(&neg_kick, 0.5);
    assert!(
        (cont.node_velocities[0][0] - (-5.0)).abs() < tol,
        "Continuum velocity after kick: expected -5, got {}",
        cont.node_velocities[0][0]
    );
    assert!(
        (cont.node_positions[0][0] - (-2.5)).abs() < tol,
        "Continuum position after kick: expected -2.5, got {}",
        cont.node_positions[0][0]
    );

    // Unknown ids and massless particles are ignored.
    let stray = InterfaceForceVec {
        forces: vec![InterfaceForce {
            id: 99,
            force: [1.0, 1.0, 1.0],
        }],
    };
    md.apply_interface_force(&stray, 0.5);
    assert!((md.velocities[0][0] - 7.0).abs() < tol);
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 6: two-pair schedule — runtime processes both pairs
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn two_pair_schedule_returns_two_reports() {
    // Set up: 4 domains (2 FEM, 2 SPH) linked by 2 FemSphCouplers.  The
    // runtime samples both domains of a pair at the coupler's own sites, so the
    // reports cover 3 and 2 sites respectively.
    let mut rt = CouplingRuntime::new();

    let dom_a = ConstantVelocityDomain::new([1.0, 0.0, 0.0], DomainKind::Fem);
    let dom_b = ConstantVelocityDomain::new([0.0, 0.0, 0.0], DomainKind::Sph);
    let dom_c = ConstantVelocityDomain::new([0.5, 0.0, 0.0], DomainKind::Fem);
    let dom_d = ConstantVelocityDomain::new([0.25, 0.0, 0.0], DomainKind::Sph);

    let idx_a = rt.add_domain(Box::new(dom_a));
    let idx_b = rt.add_domain(Box::new(dom_b));
    let idx_c = rt.add_domain(Box::new(dom_c));
    let idx_d = rt.add_domain(Box::new(dom_d));

    // Coupler 1: links dom_a ↔ dom_b with 3 sites.
    let sites_1 = vec![
        InterfaceSite {
            id: 0,
            position: [0.0, 0.0, 0.0],
        },
        InterfaceSite {
            id: 1,
            position: [1.0, 0.0, 0.0],
        },
        InterfaceSite {
            id: 2,
            position: [2.0, 0.0, 0.0],
        },
    ];
    let coupler_1 = FemSphCoupler::new(sites_1, 100.0, 10.0);

    // Coupler 2: links dom_c ↔ dom_d with 2 sites.
    let sites_2 = vec![
        InterfaceSite {
            id: 0,
            position: [0.5, 0.0, 0.0],
        },
        InterfaceSite {
            id: 1,
            position: [1.5, 0.0, 0.0],
        },
    ];
    let coupler_2 = FemSphCoupler::new(sites_2, 200.0, 20.0);

    let coup_idx_1 = rt.add_coupler(Box::new(coupler_1));
    let coup_idx_2 = rt.add_coupler(Box::new(coupler_2));

    rt.schedule_pair(coup_idx_1, idx_a, idx_b);
    rt.schedule_pair(coup_idx_2, idx_c, idx_d);

    assert_eq!(rt.domain_count(), 4);
    assert_eq!(rt.coupler_count(), 2);
    assert_eq!(rt.schedule_len(), 2);

    let reports = rt.step(0.01);

    assert_eq!(
        reports.len(),
        2,
        "expected 2 reports for 2 scheduled pairs, got {}",
        reports.len()
    );

    // Each report covers every site of its coupler (3 and 2).
    assert_eq!(
        reports[0].site_count, 3,
        "report[0].site_count should be 3, got {}",
        reports[0].site_count
    );
    assert_eq!(
        reports[1].site_count, 2,
        "report[1].site_count should be 2, got {}",
        reports[1].site_count
    );

    // Pair 1: the FEM side moves at 1 m/s and the SPH side is at rest, so the
    // damping force per site is c * (1 - 0) = 10 N and the exchanged impulse
    // c * dv * dt summed over three sites is 3 * 10 * 0.01 = 0.3 N·s along +x.
    let impulse_1 = reports[0].momentum_a_to_b[0];
    assert!(
        (impulse_1 - 0.3).abs() < 1e-12,
        "pair 1 impulse: expected 0.3, got {impulse_1}"
    );
    assert!(
        (reports[0].momentum_b_to_a[0] + impulse_1).abs() < 1e-12,
        "action-reaction violated for pair 1"
    );
    // Pair 2: c = 20, dv = 0.5 - 0.25 = 0.25, two sites: 2 * 20 * 0.25 * 0.01.
    let impulse_2 = reports[1].momentum_a_to_b[0];
    assert!(
        (impulse_2 - 0.1).abs() < 1e-12,
        "pair 2 impulse: expected 0.1, got {impulse_2}"
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Penalty law regressions
// ─────────────────────────────────────────────────────────────────────────────

/// One-site state pair helper: `(position, velocity)` for side A and side B.
fn single_site_states(
    a: ([f64; 3], [f64; 3]),
    b: ([f64; 3], [f64; 3]),
) -> (InterfaceStateVec, InterfaceStateVec) {
    let make = |(position, velocity): ([f64; 3], [f64; 3])| InterfaceStateVec {
        states: vec![InterfaceState {
            id: 0,
            position,
            velocity,
        }],
    };
    (make(a), make(b))
}

fn one_site() -> Vec<InterfaceSite> {
    vec![InterfaceSite {
        id: 0,
        position: [0.5, 0.0, 0.0],
    }]
}

#[test]
fn penalty_force_is_restoring_for_both_couplers() {
    // Side A is 1 m to the right of side B, both at rest.  The spring must pull
    // A to the left (-x) and B to the right (+x); the old law pushed A further
    // away.
    let (state_a, state_b) =
        single_site_states(([1.0, 0.0, 0.0], [0.0; 3]), ([0.0, 0.0, 0.0], [0.0; 3]));

    let mut md = MdContinuumAdapter::new(one_site(), 40.0, 0.0);
    let (fa, fb, _) = md.step(0.01, &state_a, &state_b);
    assert!(
        (fa.forces[0].force[0] + 40.0).abs() < 1e-12,
        "MD force on A"
    );
    assert!(
        (fb.forces[0].force[0] - 40.0).abs() < 1e-12,
        "MD force on B"
    );

    let mut fem = FemSphCoupler::new(one_site(), 40.0, 0.0);
    let (fa, fb, _) = fem.step(0.01, &state_a, &state_b);
    assert!(
        (fa.forces[0].force[0] + 40.0).abs() < 1e-12,
        "FEM force on A"
    );
    assert!(
        (fb.forces[0].force[0] - 40.0).abs() < 1e-12,
        "SPH force on B"
    );
}

#[test]
fn penalty_force_damps_relative_velocity_for_both_couplers() {
    // Same positions, A moves at +2 m/s relative to B: damping must brake A
    // (-x) and accelerate B (+x).
    let (state_a, state_b) = single_site_states(([0.0; 3], [2.0, 0.0, 0.0]), ([0.0; 3], [0.0; 3]));

    let mut md = MdContinuumAdapter::new(one_site(), 0.0, 5.0);
    let (fa, fb, _) = md.step(0.01, &state_a, &state_b);
    assert!((fa.forces[0].force[0] + 10.0).abs() < 1e-12);
    assert!((fb.forces[0].force[0] - 10.0).abs() < 1e-12);

    let mut fem = FemSphCoupler::new(one_site(), 0.0, 5.0);
    let (fa, fb, _) = fem.step(0.01, &state_a, &state_b);
    assert!((fa.forces[0].force[0] + 10.0).abs() < 1e-12);
    assert!((fb.forces[0].force[0] - 10.0).abs() < 1e-12);
}

#[test]
fn stiffness_acts_on_position_and_damping_on_velocity_only() {
    // The old law multiplied the stiffness by the velocity mismatch and the
    // damping by the position mismatch, i.e. it swapped the two roles (and
    // mixed N with N·s).
    let mut md = MdContinuumAdapter::new(one_site(), 1000.0, 0.0);

    // Velocity mismatch only, no damping: a pure spring exerts no force.
    let (state_a, state_b) = single_site_states(([0.0; 3], [3.0, 0.0, 0.0]), ([0.0; 3], [0.0; 3]));
    let (fa, _, report) = md.step(0.01, &state_a, &state_b);
    assert_eq!(fa.forces[0].force, [0.0; 3]);
    assert_eq!(report.momentum_a_to_b, [0.0; 3]);

    // Position mismatch only, no spring: a pure damper exerts no force.
    let mut damper = MdContinuumAdapter::new(one_site(), 0.0, 1000.0);
    let (state_a, state_b) = single_site_states(([2.0, 0.0, 0.0], [0.0; 3]), ([0.0; 3], [0.0; 3]));
    let (fa, _, report) = damper.step(0.01, &state_a, &state_b);
    assert_eq!(fa.forces[0].force, [0.0; 3]);
    assert_eq!(report.momentum_a_to_b, [0.0; 3]);
}

#[test]
fn penalty_force_does_not_depend_on_the_step_length() {
    // f is a force in newtons: only the impulse f*dt scales with dt.
    let (state_a, state_b) =
        single_site_states(([1.0, 0.0, 0.0], [0.5, 0.0, 0.0]), ([0.0; 3], [0.0; 3]));
    let mut md = MdContinuumAdapter::new(one_site(), 30.0, 4.0);
    let (fa_small, _, report_small) = md.step(0.001, &state_a, &state_b);
    let (fa_large, _, report_large) = md.step(0.1, &state_a, &state_b);
    assert_eq!(fa_small.forces[0].force, fa_large.forces[0].force);
    assert!((fa_small.forces[0].force[0] + 32.0).abs() < 1e-12);
    assert!((report_small.momentum_a_to_b[0] - 32.0 * 0.001).abs() < 1e-12);
    assert!((report_large.momentum_a_to_b[0] - 32.0 * 0.1).abs() < 1e-12);
}

// ─────────────────────────────────────────────────────────────────────────────
// Runtime regressions
// ─────────────────────────────────────────────────────────────────────────────

/// Read `(position, velocity)` of node `id` from a registered domain.
fn probe(rt: &CouplingRuntime, domain_idx: usize, id: u64) -> InterfaceState {
    let domain = rt.domain(domain_idx).expect("domain is registered");
    let sites = [InterfaceSite {
        id,
        position: [0.0; 3],
    }];
    domain.sample_interface(&sites).states[0]
}

#[test]
fn runtime_exchanges_forces_between_impulse_domains() {
    // Two 1 kg nodes at the same place; A moves at +2 m/s, B at -2 m/s.
    let mut a = ImpulseDomain::new(2, 1.0, DomainKind::Fem);
    let mut b = ImpulseDomain::new(2, 1.0, DomainKind::Sph);
    for node in 0..2 {
        a.velocities[node] = [2.0, 0.0, 0.0];
        b.velocities[node] = [-2.0, 0.0, 0.0];
    }

    let sites = vec![
        InterfaceSite {
            id: 0,
            position: [0.0; 3],
        },
        InterfaceSite {
            id: 1,
            position: [0.0, 1.0, 0.0],
        },
    ];
    let mut rt = CouplingRuntime::new();
    let ia = rt.add_domain(Box::new(a));
    let ib = rt.add_domain(Box::new(b));
    let coupler = rt.add_coupler(Box::new(FemSphCoupler::new(sites, 0.0, 3.0)));
    rt.schedule_pair(coupler, ia, ib);

    let dt = 0.01;
    let reports = rt.step(dt);
    assert_eq!(reports.len(), 1);
    assert_eq!(
        reports[0].site_count, 2,
        "both coupler sites must be sampled"
    );

    // f = c * (v_a - v_b) = 3 * 4 = 12 N per node: A is braked, B accelerated,
    // dv = f*dt/m = 0.12 m/s.
    for node in 0..2_u64 {
        let state_a = probe(&rt, ia, node);
        let state_b = probe(&rt, ib, node);
        assert!(
            (state_a.velocity[0] - (2.0 - 0.12)).abs() < 1e-12,
            "A node {node}: velocity {}",
            state_a.velocity[0]
        );
        assert!(
            (state_b.velocity[0] - (-2.0 + 0.12)).abs() < 1e-12,
            "B node {node}: velocity {}",
            state_b.velocity[0]
        );
        // Momentum is conserved node by node (equal masses).
        assert!((state_a.velocity[0] + state_b.velocity[0]).abs() < 1e-12);
    }

    // The report carries the impulse f*dt delivered to B and its reaction.
    let impulse = 2.0 * 12.0 * dt;
    assert!((reports[0].momentum_a_to_b[0] - impulse).abs() < 1e-12);
    assert!((reports[0].momentum_b_to_a[0] + impulse).abs() < 1e-12);

    // The relative velocity keeps decaying over further steps.
    let mut previous = 4.0_f64;
    for _ in 0..20 {
        rt.step(dt);
        let relative = probe(&rt, ia, 0).velocity[0] - probe(&rt, ib, 0).velocity[0];
        assert!(relative < previous, "relative velocity must decrease");
        assert!(relative > 0.0, "a damper must not reverse the motion");
        previous = relative;
    }
}

#[test]
fn runtime_md_continuum_gap_shrinks_monotonically() {
    // MD particle at x = 1, continuum node at x = 0, both at rest, 1 kg each.
    // An over-damped spring-damper (k = 100 N/m, c = 40 N·s/m, reduced mass
    // 0.5 kg => damping ratio ~2.8) must pull the two together without ever
    // overshooting.
    let md = MockMdDomain::from_particles(vec![[1.0, 0.0, 0.0]], vec![[0.0; 3]], vec![1.0]);
    let continuum = MockContinuumDomain::from_nodes(vec![[0.0; 3]], vec![[0.0; 3]], vec![1.0]);

    let mut rt = CouplingRuntime::new();
    let i_md = rt.add_domain(Box::new(md));
    let i_cont = rt.add_domain(Box::new(continuum));
    let coupler = rt.add_coupler(Box::new(MdContinuumAdapter::new(one_site(), 100.0, 40.0)));
    rt.schedule_pair(coupler, i_md, i_cont);

    let dt = 0.002;
    let mut gap = probe(&rt, i_md, 0).position[0] - probe(&rt, i_cont, 0).position[0];
    assert!((gap - 1.0).abs() < 1e-15);
    let initial_gap = gap;

    for step in 0..50 {
        let reports = rt.step(dt);
        assert_eq!(reports[0].site_count, 1);

        let md_state = probe(&rt, i_md, 0);
        let cont_state = probe(&rt, i_cont, 0);
        let new_gap = md_state.position[0] - cont_state.position[0];
        assert!(
            new_gap < gap,
            "step {step}: gap must shrink monotonically ({gap} -> {new_gap})"
        );
        assert!(new_gap > 0.0, "step {step}: over-damped, no overshoot");
        // MD moves towards the continuum (-x), the continuum towards MD (+x).
        assert!(md_state.velocity[0] < 0.0, "step {step}: MD velocity");
        assert!(
            cont_state.velocity[0] > 0.0,
            "step {step}: continuum velocity"
        );
        // Equal masses and equal-and-opposite forces: total momentum stays zero.
        assert!(
            (md_state.velocity[0] + cont_state.velocity[0]).abs() < 1e-12,
            "step {step}: momentum not conserved"
        );
        if step == 0 {
            // Both start at rest with 1 kg: the report's impulses are exactly the
            // momenta the continuum (A -> B) and the MD side (B -> A) received.
            assert!((cont_state.velocity[0] - reports[0].momentum_a_to_b[0]).abs() < 1e-12);
            assert!((md_state.velocity[0] - reports[0].momentum_b_to_a[0]).abs() < 1e-12);
            assert!((cont_state.velocity[0] - 100.0 * dt).abs() < 1e-12);
        }
        gap = new_gap;
    }
    assert!(
        gap < 0.9 * initial_gap,
        "the gap must have closed noticeably after 50 steps, got {gap}"
    );
}

#[test]
fn runtime_conserves_momentum_between_unequal_masses() {
    let md = MockMdDomain::from_particles(vec![[1.0, 0.0, 0.0]], vec![[0.0, 0.0, 0.0]], vec![4.0]);
    let continuum =
        MockContinuumDomain::from_nodes(vec![[0.0; 3]], vec![[0.0, 0.0, 0.0]], vec![1.0]);
    let mut rt = CouplingRuntime::new();
    let i_md = rt.add_domain(Box::new(md));
    let i_cont = rt.add_domain(Box::new(continuum));
    let coupler = rt.add_coupler(Box::new(MdContinuumAdapter::new(one_site(), 50.0, 8.0)));
    rt.schedule_pair(coupler, i_md, i_cont);

    for _ in 0..100 {
        rt.step(0.001);
        let p_total = 4.0 * probe(&rt, i_md, 0).velocity[0] + probe(&rt, i_cont, 0).velocity[0];
        assert!(p_total.abs() < 1e-12, "total momentum drifted to {p_total}");
    }
    // The light continuum node moved four times as fast as the heavy MD particle.
    let v_md = probe(&rt, i_md, 0).velocity[0];
    let v_cont = probe(&rt, i_cont, 0).velocity[0];
    assert!(v_md < 0.0 && v_cont > 0.0);
    assert!((v_cont / v_md + 4.0).abs() < 1e-9);
}

#[test]
fn runtime_reversed_schedule_order_gives_the_same_physics() {
    // Scheduling (coupler, B, A) instead of (coupler, A, B) swaps which domain is
    // "side A"; the penalty is symmetric under that swap, so both domains still
    // approach each other.
    let md = MockMdDomain::from_particles(vec![[1.0, 0.0, 0.0]], vec![[0.0; 3]], vec![1.0]);
    let continuum = MockContinuumDomain::from_nodes(vec![[0.0; 3]], vec![[0.0; 3]], vec![1.0]);
    let mut rt = CouplingRuntime::new();
    let i_md = rt.add_domain(Box::new(md));
    let i_cont = rt.add_domain(Box::new(continuum));
    let coupler = rt.add_coupler(Box::new(MdContinuumAdapter::new(one_site(), 100.0, 0.0)));
    // Domain order reversed: continuum is side A, MD is side B.
    rt.schedule_pair(coupler, i_cont, i_md);

    rt.step(0.01);
    assert!(
        probe(&rt, i_md, 0).velocity[0] < 0.0,
        "MD moves towards the continuum"
    );
    assert!(
        probe(&rt, i_cont, 0).velocity[0] > 0.0,
        "continuum moves towards MD"
    );
}

#[test]
fn runtime_constant_velocity_domain_advances_its_sites() {
    // A kinematic FEM wall moving at +1 m/s drags an SPH node through the
    // spring-damper; the wall itself is unaffected by the reaction force.
    let wall = ConstantVelocityDomain::new([1.0, 0.0, 0.0], DomainKind::Fem);
    let fluid = ImpulseDomain::new(1, 2.0, DomainKind::Sph);
    let mut rt = CouplingRuntime::new();
    let i_wall = rt.add_domain(Box::new(wall));
    let i_fluid = rt.add_domain(Box::new(fluid));
    let coupler = rt.add_coupler(Box::new(FemSphCoupler::new(one_site(), 0.0, 6.0)));
    rt.schedule_pair(coupler, i_wall, i_fluid);

    let dt = 0.01;
    for _ in 0..10 {
        rt.step(dt);
    }
    let wall_state = probe(&rt, i_wall, 0);
    assert!(
        (wall_state.velocity[0] - 1.0).abs() < 1e-15,
        "wall is kinematic"
    );
    // 10 steps of 10 ms: the wall site moved 0.1 m from where it was sampled.
    let site_x = probe_at(&rt, i_wall, [0.5, 0.0, 0.0]).position[0];
    assert!((site_x - 0.6).abs() < 1e-12, "wall site at x = {site_x}");

    let fluid_state = probe(&rt, i_fluid, 0);
    assert!(
        fluid_state.velocity[0] > 0.0 && fluid_state.velocity[0] < 1.0,
        "the fluid node is dragged along but not beyond the wall speed: {}",
        fluid_state.velocity[0]
    );
}

/// Sample node 0 of a registered domain at an explicit site position.
fn probe_at(rt: &CouplingRuntime, domain_idx: usize, position: [f64; 3]) -> InterfaceState {
    let domain = rt.domain(domain_idx).expect("domain is registered");
    let sites = [InterfaceSite { id: 0, position }];
    domain.sample_interface(&sites).states[0]
}

#[test]
fn runtime_self_coupling_is_a_no_op_for_antisymmetric_couplers() {
    // Both sides sample the same domain, so the penalty force is exactly zero.
    let node = ImpulseDomain::new(1, 1.0, DomainKind::Custom("loop".to_string()));
    let mut rt = CouplingRuntime::new();
    let idx = rt.add_domain(Box::new(node));
    let coupler = rt.add_coupler(Box::new(FemSphCoupler::new(one_site(), 10.0, 10.0)));
    rt.schedule_pair(coupler, idx, idx);

    let reports = rt.step(0.01);
    assert_eq!(reports.len(), 1);
    assert_eq!(reports[0].site_count, 1);
    assert_eq!(reports[0].momentum_a_to_b, [0.0; 3]);
    assert_eq!(probe(&rt, idx, 0).velocity, [0.0; 3]);
}

#[test]
fn runtime_domain_accessor_is_bounds_checked() {
    let mut rt = CouplingRuntime::new();
    assert!(rt.domain(0).is_none());
    let idx = rt.add_domain(Box::new(ConstantVelocityDomain::new(
        [0.0; 3],
        DomainKind::Lbm,
    )));
    assert_eq!(rt.domain(idx).map(|d| d.kind()), Some(DomainKind::Lbm));
    assert!(rt.domain(idx + 1).is_none());
}
