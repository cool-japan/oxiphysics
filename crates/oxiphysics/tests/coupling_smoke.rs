// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Smoke tests for the OxiCAR domain coupling layer (Blueprint KF-1).
//!
//! Six tests cover:
//! 1. Empty runtime — step produces an empty report list.
//! 2. Copy semantics — `InterfaceSite`, `InterfaceState`, `InterfaceForce`,
//!    and `CouplingReport` are all `Copy`.
//! 3. Zero-velocity exchange — identical velocities on both sides produce
//!    zero momentum exchange when `stiffness > 0, damping = 0`.
//! 4. Nonzero velocity exchange — confirms the penalty formula is correct.
//! 5. MD–continuum smoke — `MdContinuumAdapter` links two mock domains;
//!    action–reaction is satisfied.
//! 6. Two-pair schedule — runtime with two coupler pairs processes both
//!    pairs in one `step` call and returns two reports.

use oxiphysics::coupling::fem_sph::{ConstantVelocityDomain, FemSphCoupler};
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
    // f = k*(v_fem - v_sph)*dt + c*(x_fem - x_sph)
    //   = k*[1,0,0]*dt + 0
    let k = 500.0_f64;
    let dt = 0.02_f64;
    let expected_f_x = k * 1.0 * dt; // along X
    let expected_momentum_x = expected_f_x * dt; // f*dt

    let sites = vec![InterfaceSite {
        id: 0,
        position: [0.0, 0.0, 0.0],
    }];
    let mut coupler = FemSphCoupler::new(sites, 0.0, 0.0);
    coupler.stiffness = k;
    coupler.damping = 0.0;

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

    let (_fa, _fb, report) = coupler.step(dt, &state_a, &state_b);

    assert_eq!(report.site_count, 1);
    let tol = 1e-12;
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
    // stiffness=0, damping=100 → f = c*(x_md - x_cont) = [100,0,0]
    let stiffness = 0.0_f64;
    let damping = 100.0_f64;
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

    let mut adapter = MdContinuumAdapter::new(stiffness, damping);
    let (fa, fb, report) = adapter.step(dt, &state_a, &state_b);

    // f = 100 * [1,0,0] = [100,0,0]
    let expected_fx = 100.0_f64;
    let tol = 1e-10;

    assert_eq!(report.site_count, 1);
    assert!(
        (fa.forces[0].force[0] - expected_fx).abs() < tol,
        "MD force[0]: expected {expected_fx}, got {}",
        fa.forces[0].force[0]
    );
    assert!(
        (fb.forces[0].force[0] + expected_fx).abs() < tol,
        "Continuum force[0]: expected {}, got {}",
        -expected_fx,
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

    // apply_interface_force: give MD a kick of [10,0,0], mass=1 → v += 10.
    let kick = InterfaceForceVec {
        forces: vec![InterfaceForce {
            id: 0,
            force: [10.0, 0.0, 0.0],
        }],
    };
    md.apply_interface_force(&kick);
    assert!(
        (md.velocities[0][0] - 12.0).abs() < tol,
        "MD velocity after kick: expected 12, got {}",
        md.velocities[0][0]
    );

    // Continuum gets negative kick → v -= 10/1 = -10.
    let neg_kick = InterfaceForceVec {
        forces: vec![InterfaceForce {
            id: 0,
            force: [-10.0, 0.0, 0.0],
        }],
    };
    cont.apply_interface_force(&neg_kick);
    assert!(
        (cont.node_velocities[0][0] - (-10.0)).abs() < tol,
        "Continuum velocity after kick: expected -10, got {}",
        cont.node_velocities[0][0]
    );
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 6: two-pair schedule — runtime processes both pairs
// ─────────────────────────────────────────────────────────────────────────────

#[test]
fn two_pair_schedule_returns_two_reports() {
    // Set up: 4 domains (2 FEM, 2 SPH) linked by 2 FemSphCouplers.
    // All sites are trivial (empty); both couplers should produce a report
    // with site_count == 0 (no sites in coupler's list).
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

    // First report: 0 sites (runtime passes empty site slices to domain.sample_interface,
    // then coupler computes min(sites.len()=3, state_a.len()=0, state_b.len()=0) = 0).
    // The ConstantVelocityDomain returns 0 states when given 0 sites.
    assert_eq!(
        reports[0].site_count, 0,
        "report[0].site_count should be 0, got {}",
        reports[0].site_count
    );
    assert_eq!(
        reports[1].site_count, 0,
        "report[1].site_count should be 0, got {}",
        reports[1].site_count
    );
}
