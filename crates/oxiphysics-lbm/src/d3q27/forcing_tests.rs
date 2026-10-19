//! Regression tests for the D3Q27 two-relaxation-time collision (`collide_mrt`) and for the
//! Guo body force.
//!
//! The previous TRT kernel relaxed the populations in place, so the second population of every
//! opposite pair was relaxed from an already-updated partner and the result depended on the
//! direction order. The previous Guo routine modified the populations after collision with the
//! `(1 - omega/2)` prefactor but without the `F/2` velocity shift, which injected only
//! `(1 - omega/2) F` of momentum per step. These tests pin the snapshot semantics of TRT, the
//! exact momentum injection `F` per step for every collision operator, the corrected velocity
//! and the Guo momentum flux.

use super::*;

/// Deterministic pseudo-random populations around a moving equilibrium.
fn perturbed_cell(seed: u64, amplitude: f64) -> [f64; Q] {
    let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
    let feq = D3Q27Lattice::equilibrium(1.1, 0.03, -0.02, 0.025);
    std::array::from_fn(|i| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (state >> 11) as f64 / (1u64 << 53) as f64;
        feq[i] * (1.0 + amplitude * (2.0 * unit - 1.0))
    })
}

/// Mass and momentum `[rho, jx, jy, jz]` of one cell.
fn mass_and_momentum(cell: &[f64; Q]) -> [f64; 4] {
    let mut m = [0.0_f64; 4];
    for (&fi, c) in cell.iter().zip(VELOCITIES.iter()) {
        m[0] += fi;
        m[1] += fi * f64::from(c[0]);
        m[2] += fi * f64::from(c[1]);
        m[3] += fi * f64::from(c[2]);
    }
    m
}

/// Second moment `Pi_ab = sum_i f_i c_ia c_ib` of one cell.
fn momentum_flux(cell: &[f64; Q]) -> [[f64; 3]; 3] {
    let mut pi = [[0.0_f64; 3]; 3];
    for (&fi, c) in cell.iter().zip(VELOCITIES.iter()) {
        let c = [f64::from(c[0]), f64::from(c[1]), f64::from(c[2])];
        for (a, row) in pi.iter_mut().enumerate() {
            for (b, entry) in row.iter_mut().enumerate() {
                *entry += fi * c[a] * c[b];
            }
        }
    }
    pi
}

/// A single-cell lattice holding `cell`.
fn single_cell_lattice(cell: [f64; Q], omega: f64) -> D3Q27Lattice {
    let mut lat = D3Q27Lattice::new(1, 1, 1, omega);
    lat.f[0] = cell;
    lat
}

/// Reference TRT operator evaluated from a snapshot, visiting the directions in the given order.
fn reference_trt(cell: &[f64; Q], s_plus: f64, s_minus: f64, reverse: bool) -> [f64; Q] {
    let (rho, ux, uy, uz) = single_cell_lattice(*cell, 1.0).compute_macroscopic(cell);
    let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
    let order: Vec<usize> = if reverse {
        (0..Q).rev().collect()
    } else {
        (0..Q).collect()
    };
    let mut out = [0.0_f64; Q];
    for i in order {
        let opp = OPPOSITES[i];
        let neq_i = cell[i] - feq[i];
        let neq_opp = cell[opp] - feq[opp];
        out[i] = cell[i] - s_plus * 0.5 * (neq_i + neq_opp) - s_minus * 0.5 * (neq_i - neq_opp);
    }
    out
}

fn rates(s_plus: f64, s_minus: f64) -> [f64; Q] {
    let mut s = [1.0; Q];
    s[0] = s_plus;
    s[1] = s_minus;
    s
}

#[test]
fn trt_conserves_mass_and_momentum_for_unequal_rates() {
    for seed in 1..=8u64 {
        let cell = perturbed_cell(seed, 0.08);
        let mut lat = single_cell_lattice(cell, 1.0);
        lat.collide_mrt(&rates(1.2, 1.7));
        let before = mass_and_momentum(&cell);
        let after = mass_and_momentum(&lat.f[0]);
        for (c, name) in ["rho", "jx", "jy", "jz"].iter().enumerate() {
            assert!(
                (before[c] - after[c]).abs() < 1e-13,
                "seed {seed}: {name} changed from {} to {}",
                before[c],
                after[c]
            );
        }
    }
}

#[test]
fn trt_result_does_not_depend_on_the_direction_order() {
    for seed in 1..=6u64 {
        let cell = perturbed_cell(seed, 0.08);
        let mut lat = single_cell_lattice(cell, 1.0);
        lat.collide_mrt(&rates(1.2, 1.7));
        let forward = reference_trt(&cell, 1.2, 1.7, false);
        let backward = reference_trt(&cell, 1.2, 1.7, true);
        for i in 0..Q {
            assert!(
                (forward[i] - backward[i]).abs() < 1e-15,
                "seed {seed}: reference depends on the order at direction {i}"
            );
            assert!(
                (lat.f[0][i] - forward[i]).abs() < 1e-14,
                "seed {seed}: direction {i}: {} vs snapshot reference {}",
                lat.f[0][i],
                forward[i]
            );
        }
    }
}

#[test]
fn trt_with_equal_rates_is_bgk() {
    let omega = 1.35;
    let cell = perturbed_cell(11, 0.05);
    let mut trt = single_cell_lattice(cell, omega);
    let mut bgk = single_cell_lattice(cell, omega);
    trt.collide_mrt(&rates(omega, omega));
    bgk.collide_bgk();
    for i in 0..Q {
        assert!(
            (trt.f[0][i] - bgk.f[0][i]).abs() < 1e-14,
            "direction {i}: TRT {} vs BGK {}",
            trt.f[0][i],
            bgk.f[0][i]
        );
    }
}

/// Total momentum of a lattice.
fn total_momentum(lat: &D3Q27Lattice) -> [f64; 4] {
    let mut total = [0.0_f64; 4];
    for cell in &lat.f {
        let m = mass_and_momentum(cell);
        for (t, v) in total.iter_mut().zip(m.iter()) {
            *t += v;
        }
    }
    total
}

/// The three collision operators of the lattice as `(name, closure)` pairs.
fn collide_with(lat: &mut D3Q27Lattice, operator: usize, omega: f64) {
    match operator {
        0 => lat.collide_bgk(),
        1 => lat.collide_mrt(&rates(omega, 1.9 - 0.5 * omega)),
        _ => {
            let r = MrtRelaxationRates::for_viscosity(CS2 * (1.0 / omega - 0.5), 0.02);
            collide_all_mrt_full(lat, &r);
        }
    }
}

#[test]
fn every_collision_operator_injects_exactly_the_body_force_per_step() {
    let force = [2.0e-4, -1.0e-4, 3.0e-4];
    for &omega in &[0.6, 1.0, 1.4, 1.9] {
        for operator in 0..3 {
            let mut lat = D3Q27Lattice::new(3, 2, 2, omega);
            let n_cells = lat.f.len() as f64;
            for (k, cell) in lat.f.iter_mut().enumerate() {
                *cell = perturbed_cell(100 + k as u64, 0.04);
            }
            lat.apply_guo_forcing(force[0], force[1], force[2]);
            let before = total_momentum(&lat);
            collide_with(&mut lat, operator, omega);
            let after = total_momentum(&lat);
            assert!(
                (after[0] - before[0]).abs() < 1e-12,
                "operator {operator}, omega {omega}: mass changed"
            );
            for a in 0..3 {
                let injected = after[a + 1] - before[a + 1];
                assert!(
                    (injected - n_cells * force[a]).abs() < 1e-12,
                    "operator {operator}, omega {omega}: injected momentum {injected} \
                     along axis {a}, expected {}",
                    n_cells * force[a]
                );
            }
        }
    }
}

#[test]
fn reported_velocity_is_the_guo_corrected_one() {
    let (omega, fx) = (1.25, 4.0e-4);
    let mut lat = D3Q27Lattice::new(4, 3, 3, omega);
    lat.apply_guo_forcing(fx, 0.0, 0.0);
    let (rho, ux, uy, uz) = lat.compute_macroscopic(&lat.f[0]);
    assert!((rho - 1.0).abs() < 1e-14);
    assert!(
        (ux - 0.5 * fx).abs() < 1e-14,
        "the rest state with a body force reports F/(2 rho): {ux}"
    );
    assert!(uy.abs() < 1e-14 && uz.abs() < 1e-14);
    let steps = 40;
    for _ in 0..steps {
        lat.collide_and_stream();
    }
    let expected = (steps as f64 + 0.5) * fx;
    for (k, cell) in lat.f.iter().enumerate() {
        let (rho, ux, uy, uz) = lat.compute_macroscopic(cell);
        assert!((rho - 1.0).abs() < 1e-12, "cell {k}: rho = {rho}");
        assert!(
            (ux - expected).abs() < 1e-12,
            "cell {k}: ux = {ux}, expected (n + 1/2) F = {expected}"
        );
        assert!(uy.abs() < 1e-12 && uz.abs() < 1e-12, "cell {k}");
    }
    let (rho_avg, mean_ux, mean_uy, mean_uz) = lat.global_macroscopic();
    assert!((rho_avg - 1.0).abs() < 1e-12);
    assert!(
        (mean_ux - expected).abs() < 1e-12 && mean_uy.abs() < 1e-12 && mean_uz.abs() < 1e-12,
        "global velocity ({mean_ux}, {mean_uy}, {mean_uz})"
    );
}

#[test]
fn global_velocity_of_a_resting_forced_lattice_is_half_the_force_per_density() {
    let mut lat = D3Q27Lattice::new(2, 2, 2, 1.0);
    lat.apply_guo_forcing(0.0, -6.0e-4, 2.0e-4);
    let (rho, ux, uy, uz) = lat.global_macroscopic();
    assert!((rho - 1.0).abs() < 1e-14);
    assert!(ux.abs() < 1e-14);
    assert!((uy + 3.0e-4).abs() < 1e-14, "uy = {uy}");
    assert!((uz - 1.0e-4).abs() < 1e-14, "uz = {uz}");
}

/// Expected post-collision momentum flux of the Guo scheme started from a uniform equilibrium
/// state with density `rho`, velocity `u0` and force `F`, for the symmetric-part rate `omega`:
/// `Pi = rho cs^2 I + rho u0 u0 + (u0 F + F u0) + (1 - omega/4) F F / rho`.
fn expected_guo_flux(rho: f64, u0: [f64; 3], force: [f64; 3], omega: f64) -> [[f64; 3]; 3] {
    let mut pi = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let delta = if a == b { 1.0 } else { 0.0 };
            pi[a][b] = rho * CS2 * delta
                + rho * u0[a] * u0[b]
                + (u0[a] * force[b] + force[a] * u0[b])
                + (1.0 - 0.25 * omega) * force[a] * force[b] / rho;
        }
    }
    pi
}

#[test]
fn guo_forcing_reproduces_the_exact_momentum_flux() {
    let (rho, u0) = (1.15, [0.04, -0.02, 0.03]);
    let force = [0.02, 0.01, -0.015];
    for &(omega, s_minus) in &[(0.8, 0.8), (1.3, 1.3), (1.6, 0.9), (0.7, 1.8)] {
        // BGK (and TRT with equal rates) and TRT with a different anti-symmetric rate: the
        // momentum flux is an even moment, so it depends on the symmetric rate only.
        let mut lat =
            single_cell_lattice(D3Q27Lattice::equilibrium(rho, u0[0], u0[1], u0[2]), omega);
        lat.apply_guo_forcing(force[0], force[1], force[2]);
        if (omega - s_minus).abs() < 1e-15 {
            lat.collide_bgk();
        } else {
            lat.collide_mrt(&rates(omega, s_minus));
        }
        let momentum = mass_and_momentum(&lat.f[0]);
        for a in 0..3 {
            assert!(
                (momentum[a + 1] - (rho * u0[a] + force[a])).abs() < 1e-13,
                "omega {omega}, s_minus {s_minus}: momentum[{a}] = {}",
                momentum[a + 1]
            );
        }
        let pi = momentum_flux(&lat.f[0]);
        let expected = expected_guo_flux(rho, u0, force, omega);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (pi[a][b] - expected[a][b]).abs() < 1e-13,
                    "omega {omega}, s_minus {s_minus}: Pi[{a}][{b}] = {}, expected {}",
                    pi[a][b],
                    expected[a][b]
                );
            }
        }
    }
}

#[test]
fn mrt_with_uniform_rates_equals_bgk_with_a_body_force() {
    let omega = 1.2;
    let force = [3.0e-4, 1.0e-4, -2.0e-4];
    let mut bgk = D3Q27Lattice::new(2, 2, 2, omega);
    let mut mrt = D3Q27Lattice::new(2, 2, 2, omega);
    for k in 0..bgk.f.len() {
        let cell = perturbed_cell(40 + k as u64, 0.06);
        bgk.f[k] = cell;
        mrt.f[k] = cell;
    }
    bgk.apply_guo_forcing(force[0], force[1], force[2]);
    mrt.apply_guo_forcing(force[0], force[1], force[2]);
    bgk.collide_bgk();
    collide_all_mrt_full(&mut mrt, &MrtRelaxationRates::uniform(omega));
    for (k, (cell_bgk, cell_mrt)) in bgk.f.iter().zip(mrt.f.iter()).enumerate() {
        for i in 0..Q {
            assert!(
                (cell_bgk[i] - cell_mrt[i]).abs() < 1e-12,
                "cell {k}, direction {i}: BGK {} vs MRT {}",
                cell_bgk[i],
                cell_mrt[i]
            );
        }
    }
}

#[test]
fn zero_body_force_leaves_every_operator_unchanged() {
    let omega = 1.1;
    let cell = perturbed_cell(5, 0.05);
    // BGK reference computed directly from the definition.
    let (rho, ux, uy, uz) = single_cell_lattice(cell, omega).compute_macroscopic(&cell);
    let feq = D3Q27Lattice::equilibrium(rho, ux, uy, uz);
    let mut lat = single_cell_lattice(cell, omega);
    lat.collide_bgk();
    for i in 0..Q {
        let expected = cell[i] - omega * (cell[i] - feq[i]);
        assert!(
            (lat.f[0][i] - expected).abs() < 1e-15,
            "direction {i}: {} vs {expected}",
            lat.f[0][i]
        );
    }
    assert_eq!(lat.body_force, [0.0; 3]);
}

#[test]
fn guo_force_term_with_the_corrected_velocity_carries_the_full_force() {
    let omega = 1.3;
    let (rho, force) = (1.2, [1.0e-3, -2.0e-3, 5.0e-4]);
    let raw = [0.03, 0.01, -0.02];
    let (ux, uy, uz) =
        guo_corrected_velocity(raw[0], raw[1], raw[2], force[0], force[1], force[2], rho);
    let source = guo_force_term(ux, uy, uz, force[0], force[1], force[2], omega);
    let mut momentum = [0.0_f64; 3];
    let mut mass = 0.0_f64;
    for (&s, c) in source.iter().zip(VELOCITIES.iter()) {
        mass += s;
        for a in 0..3 {
            momentum[a] += s * f64::from(c[a]);
        }
    }
    assert!(mass.abs() < 1e-15, "the source creates mass: {mass}");
    for a in 0..3 {
        assert!(
            (momentum[a] - (1.0 - 0.5 * omega) * force[a]).abs() < 1e-15,
            "source momentum[{a}] = {}",
            momentum[a]
        );
    }
}
