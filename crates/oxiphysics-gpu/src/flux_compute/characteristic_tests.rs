// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`characteristic_limited_reconstruct`].
//!
//! The reconstruction used to bind the result of
//! [`characteristic_decompose_2wave`] (which returns `(w_plus, w_minus)`) as
//! `(w_minus, w_plus)`, so [`characteristic_recompose_2wave`] received the two
//! waves swapped and returned `-u1`: a constant state with `u1 = 0.5` came back
//! with `u1 = -0.5`.  It also indexed `arr[j.min(n - 1)]`, which underflows for
//! an empty grid.  The tests below pin the wave bookkeeping with exact
//! solutions (constant states, linear profiles, pure travelling waves, a
//! per-wave limiting pattern, a step) and the guards for degenerate input.

use super::*;

/// Signature of the slope limiters accepted by the reconstruction.
type Limiter = fn(f64, f64) -> f64;

/// Every slope limiter that can be passed to the reconstruction.
const LIMITERS: [(&str, Limiter); 5] = [
    ("minmod", minmod),
    ("superbee", superbee),
    ("van_leer", van_leer),
    ("mc_limiter", mc_limiter),
    ("van_albada", van_albada),
];

/// Conserved components `(u0, u1)` of the system with wave speed `a` for the
/// characteristic amplitudes `wp` (speed `+a`) and `wm` (speed `-a`), written
/// out explicitly: `u0 = w+ + w-`, `u1 = a (w+ - w-)`.
fn from_waves(wp: &[f64], wm: &[f64], a: f64) -> (Vec<f64>, Vec<f64>) {
    let u0 = wp.iter().zip(wm).map(|(p, m)| p + m).collect();
    let u1 = wp.iter().zip(wm).map(|(p, m)| a * (p - m)).collect();
    (u0, u1)
}

fn assert_state(actual: [f64; 2], expected: [f64; 2], what: &str) {
    assert!(
        (actual[0] - expected[0]).abs() < 1e-12 && (actual[1] - expected[1]).abs() < 1e-12,
        "{what}: {actual:?} vs expected {expected:?}"
    );
}

#[test]
fn constant_state_is_reconstructed_exactly() {
    let u0 = [1.0_f64; 6];
    let u1 = [0.5_f64; 6];
    for (name, limiter) in LIMITERS {
        for a in [1.0, 2.0, 0.3, -1.5] {
            let (left, right) = characteristic_limited_reconstruct(&u0, &u1, 2, a, limiter);
            assert_state(left, [1.0, 0.5], &format!("{name} a={a} left"));
            assert_state(right, [1.0, 0.5], &format!("{name} a={a} right"));
        }
    }
}

#[test]
fn linear_profile_is_reproduced_exactly_in_both_components() {
    // u0 = 1 + 0.2 j, u1 = 0.5 - 0.1 j: every characteristic amplitude is
    // linear too, so every TVD limiter returns the exact slope and both faces
    // of the interface 3.5 carry the profile value there.
    let u0: Vec<f64> = (0..8).map(|j| 1.0 + 0.2 * j as f64).collect();
    let u1: Vec<f64> = (0..8).map(|j| 0.5 - 0.1 * j as f64).collect();
    for (name, limiter) in LIMITERS {
        for a in [2.0, 0.7, -1.3] {
            let (left, right) = characteristic_limited_reconstruct(&u0, &u1, 3, a, limiter);
            let expected = [1.0 + 0.2 * 3.5, 0.5 - 0.1 * 3.5];
            assert_state(left, expected, &format!("{name} a={a} left"));
            assert_state(right, expected, &format!("{name} a={a} right"));
        }
    }
}

#[test]
fn travelling_waves_keep_the_sign_of_their_momentum() {
    let a = 3.0;
    let ramp: Vec<f64> = (0..6).map(|j| 0.1 + 0.05 * j as f64).collect();
    // u1 = a u0 is a pure right-moving wave (w- = 0); u1 = -a u0 a pure
    // left-moving one (w+ = 0).  The interface 2.5 carries u0 = 0.225.
    for (direction, sign) in [("right-moving", 1.0), ("left-moving", -1.0)] {
        let u1: Vec<f64> = ramp.iter().map(|u| sign * a * u).collect();
        let (left, right) = characteristic_limited_reconstruct(&ramp, &u1, 2, a, minmod);
        let expected = [0.225, sign * a * 0.225];
        assert_state(left, expected, &format!("{direction} left"));
        assert_state(right, expected, &format!("{direction} right"));
    }
}

#[test]
fn each_wave_is_limited_on_its_own() {
    // w+ is the ramp 0.1 j (exact slope kept); w- alternates 0.4, 0, 0.4, ...
    // so every limiter flattens it to its cell values.  Interface i = 2:
    //   w+ = 0.25 on both faces, w- = w-(2) = 0.4 left and w-(3) = 0 right.
    let a = 2.0;
    let wp: Vec<f64> = (0..8).map(|j| 0.1 * j as f64).collect();
    let wm: Vec<f64> = (0..8).map(|j| if j % 2 == 0 { 0.4 } else { 0.0 }).collect();
    let (u0, u1) = from_waves(&wp, &wm, a);
    for (name, limiter) in LIMITERS {
        let (left, right) = characteristic_limited_reconstruct(&u0, &u1, 2, a, limiter);
        assert_state(
            left,
            [0.25 + 0.4, a * (0.25 - 0.4)],
            &format!("{name} left"),
        );
        assert_state(right, [0.25, a * 0.25], &format!("{name} right"));
    }
}

#[test]
fn step_discontinuity_keeps_the_adjacent_cell_values() {
    let u0 = [1.0, 1.0, 1.0, 1.0, 3.0, 3.0, 3.0, 3.0];
    let u1 = [0.2, 0.2, 0.2, 0.2, -0.6, -0.6, -0.6, -0.6];
    for (name, limiter) in LIMITERS {
        let (left, right) = characteristic_limited_reconstruct(&u0, &u1, 3, 1.5, limiter);
        assert_state(left, [1.0, 0.2], &format!("{name} left of the jump"));
        assert_state(right, [3.0, -0.6], &format!("{name} right of the jump"));
    }
}

/// With `minmod` each limited slope is at most half the neighbouring jump, so
/// both faces of every interface stay between the two adjacent cell values in
/// each characteristic variable (the TVD property).
#[test]
fn minmod_faces_stay_between_the_adjacent_cells_in_every_wave() {
    let n = 40;
    let a = 1.3;
    let mut state = 0x9E37_79B9_7F4A_7C15_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    };
    let u0: Vec<f64> = (0..n).map(|_| next()).collect();
    let u1: Vec<f64> = (0..n).map(|_| next()).collect();
    for i in 0..n - 1 {
        let (left, right) = characteristic_limited_reconstruct(&u0, &u1, i, a, minmod);
        let (wp_i, wm_i) = characteristic_decompose_2wave(u0[i], u1[i], a);
        let (wp_n, wm_n) = characteristic_decompose_2wave(u0[i + 1], u1[i + 1], a);
        for (side, face) in [("left", left), ("right", right)] {
            let (wp, wm) = characteristic_decompose_2wave(face[0], face[1], a);
            assert!(
                wp >= wp_i.min(wp_n) - 1e-12 && wp <= wp_i.max(wp_n) + 1e-12,
                "interface {i} {side}: w+ = {wp} outside [{wp_i}, {wp_n}]"
            );
            assert!(
                wm >= wm_i.min(wm_n) - 1e-12 && wm <= wm_i.max(wm_n) + 1e-12,
                "interface {i} {side}: w- = {wm} outside [{wm_i}, {wm_n}]"
            );
        }
    }
}

#[test]
fn empty_and_mismatched_input_returns_zero_states() {
    let zero = ([0.0; 2], [0.0; 2]);
    assert_eq!(
        characteristic_limited_reconstruct(&[], &[], 0, 1.0, minmod),
        zero
    );
    assert_eq!(
        characteristic_limited_reconstruct(&[], &[], 7, 1.0, minmod),
        zero
    );
    // `u1` shorter or longer than `u0`.
    let u0 = [1.0, 2.0, 3.0];
    assert_eq!(
        characteristic_limited_reconstruct(&u0, &[0.5, 0.5], 1, 1.0, minmod),
        zero
    );
    assert_eq!(
        characteristic_limited_reconstruct(&u0, &[0.5; 4], 1, 1.0, minmod),
        zero
    );
}

#[test]
fn singular_wave_speed_returns_zero_states() {
    let u0 = [1.0, 2.0, 3.0, 4.0];
    let u1 = [0.5; 4];
    for a in [0.0, -0.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            characteristic_limited_reconstruct(&u0, &u1, 1, a, minmod),
            ([0.0; 2], [0.0; 2]),
            "a = {a}"
        );
    }
}

#[test]
fn indices_outside_the_grid_are_clamped_instead_of_panicking() {
    // u0 = 0.2 j, u1 = 0.1 - 0.05 j on five cells, a = 1.
    let u0: Vec<f64> = (0..5).map(|j| 0.2 * j as f64).collect();
    let u1: Vec<f64> = (0..5).map(|j| 0.1 - 0.05 * j as f64).collect();
    let last = [0.8, 0.1 - 0.05 * 4.0];

    // Last cell, far beyond the grid, and the largest index: zero-gradient
    // boundary, so both faces are the last cell.
    for i in [4_usize, 5, 100, usize::MAX] {
        let (left, right) = characteristic_limited_reconstruct(&u0, &u1, i, 1.0, minmod);
        assert_state(left, last, &format!("i = {i} left"));
        assert_state(right, last, &format!("i = {i} right"));
    }

    // First cell: its left slope is flattened by the boundary, the right face
    // is the exact value of the linear profile at j = 0.5.
    let (left, right) = characteristic_limited_reconstruct(&u0, &u1, 0, 1.0, minmod);
    assert_state(left, [0.0, 0.1], "i = 0 left");
    assert_state(right, [0.1, 0.1 - 0.025], "i = 0 right");
}
