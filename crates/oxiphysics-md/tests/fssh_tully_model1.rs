// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fewest-switches surface hopping on Tully's model I (simple avoided
//! crossing, J. Chem. Phys. 93, 1061 (1990)).
//!
//! A nucleus of mass 2000 m_e crosses the avoided crossing on the lower adiabatic
//! surface.  Without propagation of the electronic amplitudes the hopping
//! probability is identically zero and no trajectory ever reaches the upper
//! surface; with it, a substantial fraction is transmitted on the upper surface
//! as soon as the total energy exceeds the upper asymptote, and none below it.
//!
//! The nuclei are propagated in atomic units (velocity Verlet); the electronic
//! quantities handed to `SurfaceHoppingMd` are converted to its eV / fs units.

use oxiphysics_md::photochemistry_md::SurfaceHoppingMd;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// Hartree in eV.
const HARTREE_EV: f64 = 27.211_386_245_988;
/// Atomic unit of time in fs.
const AU_TIME_FS: f64 = 0.024_188_843_265_857;

/// Nuclear mass (electron masses).
const MASS: f64 = 2000.0;
const A: f64 = 0.01;
const B: f64 = 1.6;
const C: f64 = 0.005;
const D: f64 = 1.0;

/// Diabatic matrix elements `(V11, V12, dV11/dx, dV12/dx)`; `V22 = -V11`.
fn diabatic(x: f64) -> (f64, f64, f64, f64) {
    let (v11, dv11) = if x > 0.0 {
        (A * (1.0 - (-B * x).exp()), A * B * (-B * x).exp())
    } else {
        (-A * (1.0 - (B * x).exp()), A * B * (B * x).exp())
    };
    let v12 = C * (-D * x * x).exp();
    (v11, v12, dv11, -2.0 * D * x * v12)
}

/// Adiabatic energy magnitude `R = sqrt(V11^2 + V12^2)` (the surfaces are `-R`, `+R`),
/// its derivative and the non-adiabatic coupling `<lower| d/dx |upper>`.
fn adiabatic(x: f64) -> (f64, f64, f64) {
    let (v11, v12, dv11, dv12) = diabatic(x);
    let r = (v11 * v11 + v12 * v12).sqrt();
    let dr = (v11 * dv11 + v12 * dv12) / r;
    let coupling = 0.5 * (v11 * dv12 - v12 * dv11) / (r * r);
    (r, dr, coupling)
}

/// Force on the nucleus on adiabatic surface `state` (0 = lower, 1 = upper).
fn force(x: f64, state: usize) -> f64 {
    let (_, dr, _) = adiabatic(x);
    if state == 0 { dr } else { -dr }
}

/// Outcome of one trajectory.
struct Outcome {
    transmitted: bool,
    final_state: usize,
    /// `|E_total(end) - E_total(start)|` in Hartree.
    energy_drift: f64,
}

/// One FSSH trajectory from `x = -5` with momentum `k0` on the lower surface.
fn run_trajectory(k0: f64, rng: &mut StdRng) -> Outcome {
    let dt = 1.0; // atomic units
    let mut x = -5.0_f64;
    let mut v = k0 / MASS;
    let mut shmd = SurfaceHoppingMd::new(2, 0, dt * AU_TIME_FS).expect("valid FSSH state");
    let (r0, _, _) = adiabatic(x);
    let total_start = 0.5 * MASS * v * v - r0;
    let mut accel = force(x, shmd.current_state) / MASS;
    for _ in 0..20_000 {
        x += v * dt + 0.5 * accel * dt * dt;
        let accel_new = force(x, shmd.current_state) / MASS;
        v += 0.5 * (accel + accel_new) * dt;

        let (r, _, coupling) = adiabatic(x);
        shmd.set_hamiltonian(0, 0, -r * HARTREE_EV);
        shmd.set_hamiltonian(1, 1, r * HARTREE_EV);
        shmd.set_nac_velocity(0, 1, v * coupling / AU_TIME_FS);

        let kinetic_ev = 0.5 * MASS * v * v * HARTREE_EV;
        let energies = [-r * HARTREE_EV, r * HARTREE_EV];
        let draw: f64 = rng.random();
        if let Some((_, scale)) = shmd.attempt_hop(draw, kinetic_ev, &energies) {
            v *= scale;
        }
        shmd.propagate_electronic(4);
        accel = force(x, shmd.current_state) / MASS;
        if x.abs() > 5.0 {
            break;
        }
    }
    let (r_end, _, _) = adiabatic(x);
    let potential = if shmd.current_state == 0 {
        -r_end
    } else {
        r_end
    };
    Outcome {
        transmitted: x > 5.0,
        final_state: shmd.current_state,
        energy_drift: (0.5 * MASS * v * v + potential - total_start).abs(),
    }
}

/// Final-state statistics of an ensemble of trajectories.
struct Ensemble {
    upper_transmitted: f64,
    lower_transmitted: f64,
    upper_reflected: f64,
    lower_reflected: f64,
    worst_energy_drift: f64,
}

fn ensemble(k0: f64, n_traj: usize, seed: u64) -> Ensemble {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut counts = [[0_usize; 2]; 2]; // [transmitted][final state]
    let mut worst_energy_drift = 0.0_f64;
    for _ in 0..n_traj {
        let outcome = run_trajectory(k0, &mut rng);
        worst_energy_drift = worst_energy_drift.max(outcome.energy_drift);
        counts[usize::from(outcome.transmitted)][outcome.final_state] += 1;
    }
    let n = n_traj as f64;
    Ensemble {
        upper_transmitted: counts[1][1] as f64 / n,
        lower_transmitted: counts[1][0] as f64 / n,
        upper_reflected: counts[0][1] as f64 / n,
        lower_reflected: counts[0][0] as f64 / n,
        worst_energy_drift,
    }
}

#[test]
fn upper_surface_transmission_is_nonzero_above_the_energy_threshold() {
    // k0 = 20 a.u. (E = 0.1 Eh) is far above the 2A = 0.02 Eh needed to end on
    // the upper surface, and the avoided crossing is passed non-adiabatically.
    let e = ensemble(20.0, 300, 1990);
    assert!(
        e.upper_transmitted > 0.2 && e.upper_transmitted < 0.8,
        "upper-state transmission {} (lower {}, reflected {})",
        e.upper_transmitted,
        e.lower_transmitted,
        e.upper_reflected + e.lower_reflected
    );
    assert!(e.lower_transmitted > 0.2, "lower: {}", e.lower_transmitted);
    let total = e.upper_transmitted + e.lower_transmitted + e.upper_reflected + e.lower_reflected;
    assert!((total - 1.0).abs() < 1e-12);
    // The velocity rescaling at every accepted hop conserves the total energy.
    assert!(
        e.worst_energy_drift < 2e-4,
        "total-energy drift {} Eh",
        e.worst_energy_drift
    );
}

#[test]
fn nothing_is_transmitted_on_the_upper_surface_below_the_energy_threshold() {
    // sqrt(4 M A) = 8.94 a.u.: at k0 = 8.9 the total energy (0.0098 Eh) lies above
    // the upper surface at the crossing (R(0) = 0.005 Eh), so hops are energetically
    // possible, but below the upper asymptote (0.01 Eh): a trajectory that hopped up
    // is turned back, and none is transmitted on the upper surface.
    let e = ensemble(8.9, 300, 85);
    assert_eq!(e.upper_transmitted, 0.0);
    assert!(e.lower_transmitted > 0.8, "lower: {}", e.lower_transmitted);
    assert!(
        e.worst_energy_drift < 2e-4,
        "total-energy drift {} Eh",
        e.worst_energy_drift
    );
}
