// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the tribochemistry MD module: Morse force sign, bond
//! virial, virial pressure, shear stress, the periodic / Lees-Edwards
//! minimum-image convention and the stress-driven mechanochemical drivers.

use super::*;

const INV_SQRT2: f64 = std::f64::consts::FRAC_1_SQRT_2;

/// `true` when `a` and `b` agree to the relative tolerance `rel` or to the
/// absolute tolerance `abs_tol`.
fn close(a: f64, b: f64, rel: f64, abs_tol: f64) -> bool {
    let diff = (a - b).abs();
    diff <= abs_tol || diff <= rel * a.abs().max(b.abs())
}

/// Periodic box with lengths `[lx, ly, lz]` and Lees-Edwards offset `tilt`.
fn sheared_box(lx: f64, ly: f64, lz: f64, tilt: f64) -> TriboBox {
    TriboBox {
        lengths: [lx, ly, lz],
        shear_rate: 0.0,
        strain: tilt,
        tilt,
    }
}

// ---------------------------------------------------------------------------
// Pair potentials: sign convention of the force
// ---------------------------------------------------------------------------

#[test]
fn morse_force_is_attractive_when_stretched_and_repulsive_when_compressed() {
    for morse in [MorseParams::cc(), MorseParams::co()] {
        assert!(
            morse.force_magnitude(1.1 * morse.r_eq) < 0.0,
            "a stretched Morse bond must attract (force = -dU/dr < 0)"
        );
        assert!(
            morse.force_magnitude(0.9 * morse.r_eq) > 0.0,
            "a compressed Morse bond must repel (force = -dU/dr > 0)"
        );
        assert_eq!(morse.force_magnitude(morse.r_eq), 0.0);
    }
}

#[test]
fn morse_force_vector_pulls_a_stretched_bond_together() {
    let morse = MorseParams::cc();
    let pos_i = [0.0, 0.0, 0.0];
    let stretched_j = [1.1 * morse.r_eq, 0.0, 0.0];
    let f_i = morse.force_vector(pos_i, stretched_j);
    let f_j = morse.force_vector(stretched_j, pos_i);
    assert!(f_i[0] > 0.0, "atom i must be pulled towards j, got {f_i:?}");
    assert!(f_j[0] < 0.0, "atom j must be pulled towards i, got {f_j:?}");
    assert!(close(f_i[0], -f_j[0], 1.0e-12, 0.0), "Newton's third law");
    assert_eq!(f_i[1], 0.0);
    assert_eq!(f_i[2], 0.0);
    let compressed_j = [0.9 * morse.r_eq, 0.0, 0.0];
    assert!(
        morse.force_vector(pos_i, compressed_j)[0] < 0.0,
        "a compressed bond must push atom i away from j"
    );
}

#[test]
fn morse_force_matches_negative_energy_gradient() {
    for morse in [MorseParams::cc(), MorseParams::co()] {
        let h = 1.0e-5 * morse.r_eq;
        for k in 0..=24 {
            let r = morse.r_eq * (0.8 + 1.2 * f64::from(k) / 24.0);
            if (r - morse.r_eq).abs() < 1.0e-3 * morse.r_eq {
                continue;
            }
            let numeric = -(morse.energy(r + h) - morse.energy(r - h)) / (2.0 * h);
            let analytic = morse.force_magnitude(r);
            assert!(
                close(analytic, numeric, 1.0e-6, 1.0e-14),
                "r/r_eq = {}: F = {analytic:e}, -dU/dr = {numeric:e}",
                r / morse.r_eq
            );
        }
    }
}

#[test]
fn lj_force_matches_negative_energy_gradient() {
    for lj in [LjParams::cc(), LjParams::c_fe(), LjParams::fe_fe()] {
        let h = 1.0e-5 * lj.sigma;
        for k in 0..=21 {
            let r = lj.sigma * (0.9 + 0.1 * f64::from(k));
            let r_min = 2.0_f64.powf(1.0 / 6.0) * lj.sigma;
            if (r - r_min).abs() < 0.02 * lj.sigma {
                continue;
            }
            let numeric = -(lj.energy(r + h) - lj.energy(r - h)) / (2.0 * h);
            let analytic = lj.force_magnitude(r);
            assert!(
                close(analytic, numeric, 1.0e-6, 1.0e-16),
                "r/sigma = {}: F = {analytic:e}, -dU/dr = {numeric:e}",
                r / lj.sigma
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Minimum image in a periodic / Lees-Edwards box
// ---------------------------------------------------------------------------

#[test]
fn min_image_without_tilt_removes_whole_periods() {
    let b = TriboBox::cubic(10.0);
    let d = b.min_image([9.0, -9.5, 6.0]);
    assert!(close(d[0], -1.0, 1.0e-14, 1.0e-14), "dx = {}", d[0]);
    assert!(close(d[1], 0.5, 1.0e-14, 1.0e-14), "dy = {}", d[1]);
    assert!(close(d[2], -4.0, 1.0e-14, 1.0e-14), "dz = {}", d[2]);
}

#[test]
fn min_image_applies_the_lees_edwards_offset_across_z() {
    // The image one box above is displaced by +0.3 along x.
    let b = sheared_box(1.0, 1.0, 1.0, 0.3);
    let above = b.min_image([0.1, 0.0, 0.9]);
    assert!(close(above[0], -0.2, 1.0e-14, 1.0e-14), "dx = {}", above[0]);
    assert!(close(above[2], -0.1, 1.0e-14, 1.0e-14), "dz = {}", above[2]);
    let below = b.min_image([0.1, 0.0, -0.9]);
    assert!(close(below[0], 0.4, 1.0e-14, 1.0e-14), "dx = {}", below[0]);
    assert!(close(below[2], 0.1, 1.0e-14, 1.0e-14), "dz = {}", below[2]);
}

#[test]
fn min_image_inverts_every_lattice_translation() {
    // Lx = 2, Ly = 1.5, Lz = 1 with an image offset of 0.8 <= Lx / 2.
    let b = sheared_box(2.0, 1.5, 1.0, 0.8);
    let offset = b.tilt * b.lengths[2];
    let shorts = [
        [0.2, -0.1, 0.15],
        [-0.3, 0.2, -0.25],
        [0.1, 0.0, 0.0],
        [0.0, 0.3, -0.2],
    ];
    for s in shorts {
        for nx in -3..=3 {
            for ny in -3..=3 {
                for nz in -3..=3 {
                    let raw = [
                        s[0] + f64::from(nx) * b.lengths[0] + f64::from(nz) * offset,
                        s[1] + f64::from(ny) * b.lengths[1],
                        s[2] + f64::from(nz) * b.lengths[2],
                    ];
                    let d = b.min_image(raw);
                    for k in 0..3 {
                        assert!(
                            (d[k] - s[k]).abs() < 1.0e-12,
                            "images ({nx},{ny},{nz}) axis {k}: {} vs {}",
                            d[k],
                            s[k]
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn min_image_agrees_with_lees_edwards_wrapping() {
    let mut b = TriboBox::cubic(1.0e-9);
    b.tilt = 0.27;
    b.shear_rate = 3.0e11;
    // Leaving through the top face ...
    let above = [0.4e-9, 0.3e-9, 1.05e-9];
    let (wrapped_above, v_above) = b.wrap_lees_edwards(above, [0.0; 3]);
    assert!((0.0..b.lengths[2]).contains(&wrapped_above[2]));
    let sep = b.min_image(sub3(above, wrapped_above));
    assert!(
        norm3(sep) < 1.0e-20,
        "an atom and its Lees-Edwards image must coincide, got {sep:?}"
    );
    assert!(
        (v_above[0] + b.shear_rate * b.lengths[2]).abs() < 1.0e-6,
        "x velocity must drop by shear_rate * Lz, got {}",
        v_above[0]
    );
    // ... and through the bottom face.
    let below = [0.4e-9, 0.3e-9, -0.05e-9];
    let (wrapped_below, _) = b.wrap_lees_edwards(below, [0.0; 3]);
    let sep_below = b.min_image(sub3(below, wrapped_below));
    assert!(norm3(sep_below) < 1.0e-20, "got {sep_below:?}");
}

#[test]
fn shear_step_folds_the_offset_modulo_lx() {
    // Lx = 2 Lz: an offset tilt * Lz is periodic with period Lx, i.e. tilt with period 2.
    let mut b = TriboBox {
        lengths: [2.0, 2.0, 1.0],
        shear_rate: 1.3,
        strain: 0.0,
        tilt: 0.0,
    };
    let unfolded = TriboBox {
        tilt: 1.3,
        ..b.clone()
    };
    b.shear_step(1.0);
    assert!(close(b.tilt, -0.7, 1.0e-14, 1.0e-14), "tilt = {}", b.tilt);
    assert!(
        close(b.strain, 1.3, 1.0e-14, 1.0e-14),
        "the accumulated strain must not be folded, got {}",
        b.strain
    );
    // The folded and unfolded boxes describe the same lattice.
    for dr in [[0.3, 0.1, 0.4], [-0.5, 0.2, 0.9], [1.1, -0.4, -0.8]] {
        let a = b.min_image(dr);
        let c = unfolded.min_image(dr);
        for k in 0..3 {
            assert!(
                (a[k] - c[k]).abs() < 1.0e-12,
                "axis {k}: {} vs {}",
                a[k],
                c[k]
            );
        }
    }
    // A cubic box keeps the classic (-0.5, 0.5] window.
    let mut cubic = TriboBox::cubic(1.0);
    cubic.shear_rate = 0.7;
    cubic.shear_step(1.0);
    assert!(
        close(cubic.tilt, -0.3, 1.0e-14, 1.0e-14),
        "tilt = {}",
        cubic.tilt
    );
}

#[test]
fn stretching_force_helpers_are_consistent() {
    let bond = TriboBond::new_cc(0, 1);
    let pos_i = [1.0e-10, 2.0e-10, 0.5e-10];
    let pos_j = [2.3e-10, 2.4e-10, 0.1e-10];
    let from_positions = bond.stretching_force(pos_i, pos_j);
    let from_separation = bond.stretching_force_from_separation(sub3(pos_i, pos_j));
    assert_eq!(from_positions, from_separation);
    assert!(close(
        bond.stretching_force_from_separation([2.0 * bond.length_eq, 0.0, 0.0]),
        bond.k_spring * bond.length_eq,
        1.0e-12,
        0.0
    ));
}

// ---------------------------------------------------------------------------
// Bonds across a periodic boundary
// ---------------------------------------------------------------------------

#[test]
fn bond_straddling_a_boundary_stays_at_its_equilibrium_length_while_crossing_it() {
    let lz = 5.0e-9;
    let config = TriboConfig {
        dt: 1.0e-15,
        temperature: 300.0,
        ..Default::default()
    };
    let mut sim = TriboSimulation::new(config, lz);
    let r0 = BOND_CC;
    let mut a = TriboAtom::new([2.5e-9, 2.5e-9, 0.5 * r0], TriboChem::C);
    let mut b = TriboAtom::new([2.5e-9, 2.5e-9, lz - 0.5 * r0], TriboChem::C);
    a.vel = [0.0, 0.0, -1.5e3];
    b.vel = [0.0, 0.0, -1.5e3];
    sim.add_atom(a);
    sim.add_atom(b);
    sim.add_bond(TriboBond::new_cc(0, 1));

    // The raw separation is ~Lz but the minimum image is exactly the bond
    // length, so the bond is force free.  (The old code produced ~2e-6 N.)
    sim.compute_bond_forces();
    let f_max = sim
        .atoms
        .iter()
        .flat_map(|atom| atom.force)
        .fold(0.0_f64, |m, f| m.max(f.abs()));
    assert!(
        f_max < 1.0e-12,
        "spurious boundary-spanning force: {f_max:e} N"
    );

    // Tiny random numbers: a spurious box-length tension (breaking probability
    // ~1 % per step) would break the bond, a relaxed bond (~1e-40) never does.
    for _ in 0..100 {
        sim.step_forward(&[1.0e-6]);
    }
    assert!(
        !sim.bonds[0].broken,
        "bond broke spuriously across the boundary"
    );
    let dr = sim.box_.min_image(sub3(sim.atoms[0].pos, sim.atoms[1].pos));
    assert!(
        (norm3(dr) - r0).abs() < 1.0e-4 * r0,
        "minimum-image bond length {} m, equilibrium {r0} m",
        norm3(dr)
    );
    assert!(
        sim.atoms.iter().any(|atom| atom.pos[2] > 0.5 * lz),
        "the leading atom should have been wrapped through z = 0"
    );
}

// ---------------------------------------------------------------------------
// Virial, pressure and shear stress
// ---------------------------------------------------------------------------

#[test]
fn pressure_of_a_bond_free_gas_is_n_kb_t_over_v() {
    let mut sim = TriboSimulation::new(TriboConfig::default(), 5.0e-9);
    let velocities = [
        [120.0, -40.0, 15.0],
        [-300.0, 220.0, 90.0],
        [10.0, 5.0, -450.0],
        [260.0, -310.0, 75.0],
        [-95.0, 180.0, -220.0],
    ];
    for (k, v) in velocities.iter().enumerate() {
        let species = if k % 2 == 0 {
            TriboChem::C
        } else {
            TriboChem::Fe
        };
        let mut atom = TriboAtom::new([1.0e-9 * (k as f64 + 1.0), 2.0e-9, 3.0e-9], species);
        atom.vel = *v;
        sim.add_atom(atom);
    }
    let n = sim.atoms.len() as f64;
    let expected = n * K_B * sim.temperature() / sim.box_.volume();
    let p = sim.pressure();
    assert!(expected > 0.0);
    assert!(
        (p - expected).abs() <= 1.0e-12 * expected,
        "P = {p} Pa, N k_B T / V = {expected} Pa (the old formula divided by N)"
    );
}

/// Dimer stretched along (1, 0, 1)/sqrt(2) inside a 4 nm box.
fn stretched_xz_dimer(stretch: f64) -> (TriboSimulation, f64) {
    let lbox = 4.0e-9;
    let mut sim = TriboSimulation::new(TriboConfig::default(), lbox);
    let d = stretch * BOND_CC;
    sim.add_atom(TriboAtom::new([2.0e-9, 2.0e-9, 2.0e-9], TriboChem::C));
    sim.add_atom(TriboAtom::new(
        [2.0e-9 + d * INV_SQRT2, 2.0e-9, 2.0e-9 + d * INV_SQRT2],
        TriboChem::C,
    ));
    sim.add_bond(TriboBond::new_cc(0, 1));
    (sim, d)
}

fn total_virial(sim: &TriboSimulation) -> [f64; 6] {
    let mut total = [0.0; 6];
    for atom in &sim.atoms {
        for (t, w) in total.iter_mut().zip(atom.virial.iter()) {
            *t += w;
        }
    }
    total
}

#[test]
fn stretched_bond_accumulates_virial_pressure_and_shear_stress() {
    let (mut sim, d) = stretched_xz_dimer(1.25);
    sim.compute_bond_forces();
    let tension = 400.0 * (d - BOND_CC);
    assert!(tension > 0.0);
    // The force on atom 0 points towards atom 1, i.e. along +(1, 0, 1)/sqrt(2).
    let f0 = sim.atoms[0].force;
    assert!(
        close(f0[0], tension * INV_SQRT2, 1.0e-12, 0.0),
        "f0 = {f0:?}"
    );
    assert!(
        close(f0[2], tension * INV_SQRT2, 1.0e-12, 0.0),
        "f0 = {f0:?}"
    );
    assert_eq!(f0[1], 0.0);
    let f1 = sim.atoms[1].force;
    for k in 0..3 {
        assert!(close(f0[k], -f1[k], 1.0e-12, 0.0), "Newton's third law");
    }
    // W = r_01 (x) f_0 with r_01 = -d n and f_0 = +T n: every in-plane
    // component equals -d T / 2, [xx, yy, zz, xy, xz, yz].
    let w = total_virial(&sim);
    let half = -0.5 * d * tension;
    let expected = [half, 0.0, half, 0.0, half, 0.0];
    for k in 0..6 {
        assert!(
            close(w[k], expected[k], 1.0e-9, 0.0),
            "virial component {k}: {} vs {}",
            w[k],
            expected[k]
        );
    }
    // The pair virial is shared half and half between the two atoms.
    for k in 0..6 {
        assert!(close(
            sim.atoms[0].virial[k],
            sim.atoms[1].virial[k],
            1.0e-12,
            0.0
        ));
    }
    let vol = sim.box_.volume();
    // Atoms at rest: P = tr W / (3 V) < 0 for a bond in tension.
    let pressure = sim.pressure();
    assert!(
        close(pressure, -d * tension / (3.0 * vol), 1.0e-9, 0.0),
        "P = {pressure} Pa"
    );
    // sigma_xz = -W_xz / V > 0 for a bond in tension.
    let tau = sim.shear_stress();
    assert!(tau > 0.0, "shear stress must be non-zero, got {tau}");
    assert!(
        close(tau, 0.5 * d * tension / vol, 1.0e-9, 0.0),
        "sigma_xz = {tau} Pa, expected {}",
        0.5 * d * tension / vol
    );
}

#[test]
fn virial_is_rebuilt_not_accumulated_and_ignores_broken_bonds() {
    let (mut sim, _) = stretched_xz_dimer(1.25);
    sim.compute_bond_forces();
    let first = total_virial(&sim);
    sim.compute_bond_forces();
    sim.compute_bond_forces();
    let again = total_virial(&sim);
    for k in 0..6 {
        assert_eq!(
            first[k], again[k],
            "component {k} drifted between evaluations"
        );
    }
    sim.bonds[0].broken = true;
    sim.compute_bond_forces();
    assert!(total_virial(&sim).iter().all(|w| *w == 0.0));
    assert!(sim.atoms.iter().all(|a| a.force == [0.0; 3]));
    assert_eq!(sim.shear_stress(), 0.0);
}

#[test]
fn shear_stress_kinetic_term_uses_the_peculiar_velocity() {
    let lbox = 2.0e-9;
    let mut sim = TriboSimulation::new(TriboConfig::default(), lbox);
    sim.box_.shear_rate = 4.0e10;
    let mut atom = TriboAtom::new([1.0e-9, 1.0e-9, 1.5e-9], TriboChem::C);
    atom.vel = [5.0, 0.0, 200.0];
    sim.add_atom(atom);
    // The streaming speed at z = 1.5 nm is 4e10 * (1.5 - 1.0) nm = 20 m/s,
    // so the peculiar x velocity is 5 - 20 = -15 m/s.
    let peculiar_vx = 5.0 - 4.0e10 * (1.5e-9 - 1.0e-9);
    let expected = -(TriboChem::C.mass() * peculiar_vx * 200.0) / sim.box_.volume();
    let tau = sim.shear_stress();
    assert!(tau > 0.0);
    assert!(
        close(tau, expected, 1.0e-12, 0.0),
        "tau = {tau}, expected {expected}"
    );
}

#[test]
fn atoms_moving_with_the_streaming_flow_carry_no_kinetic_shear_stress() {
    let mut sim = TriboSimulation::new(TriboConfig::default(), 2.0e-9);
    sim.box_.shear_rate = 5.0e10;
    for (k, z) in [0.3e-9, 0.9e-9, 1.4e-9, 1.8e-9].into_iter().enumerate() {
        let mut atom = TriboAtom::new([0.5e-9 * (k as f64 + 1.0), 1.0e-9, z], TriboChem::C);
        atom.vel = [
            sim.box_.streaming_velocity_x(z),
            0.0,
            300.0 - 80.0 * k as f64,
        ];
        sim.add_atom(atom);
    }
    assert_eq!(sim.shear_stress(), 0.0);
}

// ---------------------------------------------------------------------------
// Stress-driven mechanochemical drivers
// ---------------------------------------------------------------------------

/// A C-C dimer whose bond crosses the z face and lies along -(1, 0, 1)/sqrt(2)
/// in the minimum-image convention, in a 1 nm box sheared at `sliding_velocity`.
/// The Lees-Edwards offset grows with time and stretches that bond.
fn sheared_dimer_sim(sliding_velocity: f64) -> TriboSimulation {
    let config = TriboConfig {
        dt: 1.0e-15,
        temperature: 300.0,
        normal_stress: 1.0e9,
        sliding_velocity,
        reactive: false,
        ..Default::default()
    };
    let lbox = 1.0e-9;
    let a = BOND_CC * INV_SQRT2;
    let mut sim = TriboSimulation::new(config, lbox);
    sim.add_atom(TriboAtom::new(
        [0.6e-9 - a, 0.5e-9, lbox - 0.5 * a],
        TriboChem::C,
    ));
    sim.add_atom(TriboAtom::new([0.6e-9, 0.5e-9, 0.5 * a], TriboChem::C));
    sim.add_bond(TriboBond::new_cc(0, 1));
    sim.add_additive(LubricantAdditive::zddp(1000));
    sim.phase.tau_crit = 1.0e6;
    sim.phase.k_transform = 1.0e-3;
    sim
}

fn peak_abs(history: &[f64]) -> f64 {
    history.iter().fold(0.0_f64, |m, s| m.max(s.abs()))
}

#[test]
fn shearing_a_bonded_system_produces_a_shear_stress_history() {
    let mut sims = [0.0, 100.0, 400.0].map(sheared_dimer_sim);
    for sim in &mut sims {
        for _ in 0..50 {
            sim.step_forward(&[]);
        }
    }
    // Before the fix the history was identically zero (nothing wrote the virial).
    let still = peak_abs(&sims[0].shear_stress_history);
    let slow = peak_abs(&sims[1].shear_stress_history);
    let fast = peak_abs(&sims[2].shear_stress_history);
    assert_eq!(sims[1].shear_stress_history.len(), 50);
    assert!(still < 1.0, "unsheared dimer carries {still} Pa");
    assert!(slow > 1.0e5, "sheared dimer carries only {slow} Pa");
    assert!(fast > 1.0e5, "sheared dimer carries only {fast} Pa");
}

#[test]
fn mechanochemical_drivers_receive_the_measured_stress_magnitude() {
    let mut sim = sheared_dimer_sim(200.0);
    for _ in 0..50 {
        sim.step_forward(&[]);
    }
    let last = *sim
        .shear_stress_history
        .last()
        .expect("fifty steps were recorded");
    assert!(last != 0.0);
    let add = &sim.additives[0];
    assert_eq!(
        add.shear_stress,
        last.abs(),
        "the additive must see |sigma_xz|, not a fraction of the normal stress"
    );
    // |tau| >= 0, so the mechanochemical factor exp(gamma * tau / kT) can only
    // raise the decomposition rate; a signed negative stress would lower it.
    assert!(add.total_rate() >= add.k_thermal);
    assert!(
        sim.phase.amorphous_fraction > 0.0,
        "the phase transformer must be driven by the measured stress"
    );
}

#[test]
fn additive_rate_depends_on_the_sliding_velocity_through_the_stress() {
    let rate_after = |sliding_velocity: f64| -> f64 {
        let mut sim = sheared_dimer_sim(sliding_velocity);
        for _ in 0..50 {
            sim.step_forward(&[]);
        }
        sim.additives[0].total_rate()
    };
    let k_thermal = LubricantAdditive::zddp(1).k_thermal;
    let still = rate_after(0.0);
    let slow = rate_after(150.0);
    let fast = rate_after(400.0);
    assert!(close(still, k_thermal, 1.0e-9, 0.0), "static rate {still}");
    assert!(slow > k_thermal * (1.0 + 1.0e-9), "slow shear rate {slow}");
    assert!(fast > k_thermal * (1.0 + 1.0e-9), "fast shear rate {fast}");
    assert!(
        (fast - slow).abs() > 1.0e-6 * slow,
        "rate must change with the sliding velocity: {slow} vs {fast}"
    );
}

#[test]
fn drivers_ignore_the_applied_normal_stress_when_nothing_is_sheared() {
    let mut sim = sheared_dimer_sim(0.0);
    sim.config.normal_stress = 5.0e9;
    for _ in 0..20 {
        sim.step_forward(&[]);
    }
    // The former code fed 0.1 * normal_stress (5e8 Pa) to the additives and
    // 0.15 * normal_stress to the phase transformer regardless of the state.
    assert!(
        sim.additives[0].shear_stress < 1.0,
        "{}",
        sim.additives[0].shear_stress
    );
    assert_eq!(sim.phase.amorphous_fraction, 0.0);
    assert!(peak_abs(&sim.shear_stress_history) < 1.0);
}
