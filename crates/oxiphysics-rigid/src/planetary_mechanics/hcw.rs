// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Closed-form Hill-Clohessy-Wiltshire (HCW) two-impulse targeting.
//!
//! The relative motion of a chaser about a circular reference orbit with mean
//! motion `n` obeys, in the rotating frame with `x` radial, `y` along-track
//! and `z` cross-track,
//!
//! ```text
//! x'' =  3 n² x + 2 n y'
//! y'' = -2 n x'
//! z'' =  -n² z
//! ```
//!
//! whose exact state-transition matrix over a time `t` (with `φ = n t`,
//! `s = sin φ`, `c = cos φ`) is
//!
//! ```text
//! Φrr = [[4-3c, 0, 0], [6(s-φ), 1, 0], [0, 0, c]]
//! Φrv = [[s/n, 2(1-c)/n, 0], [-2(1-c)/n, (4s-3φ)/n, 0], [0, 0, s/n]]
//! Φvr = [[3ns, 0, 0], [6n(c-1), 0, 0], [0, 0, -ns]]
//! Φvv = [[c, 2s, 0], [-2s, 4c-3, 0], [0, 0, c]]
//! ```
//!
//! This is the same solution used by
//! [`HillClohessyWiltshire::propagate`](super::HillClohessyWiltshire::propagate)
//! and by `SpacecraftFormation::propagate_exact`; [`hcw_two_impulse`] inverts
//! it to find the two burns that carry a chaser from one relative state to
//! another in a prescribed time.

/// Relative threshold under which a targeting block is considered singular.
const SINGULARITY_TOL: f64 = 1e-10;

fn is_finite_vec(v: [f64; 3]) -> bool {
    v.iter().all(|c| c.is_finite())
}

/// Two-impulse HCW targeting between two relative states.
///
/// A chaser at relative position `r0` with relative velocity `v0` (frame:
/// `x` radial, `y` along-track, `z` cross-track) is to reach position `rf`
/// with velocity `vf` after `tf` seconds.  `n` is the mean motion of the
/// circular reference orbit in rad/s.
///
/// The first burn `dv1` at `t = 0` selects the departure velocity
/// `v0⁺ = Φrv⁻¹ (rf - Φrr r0)` so that the coasting arc passes through `rf` at
/// `tf`; the second burn `dv2` at `tf` removes the difference between the
/// arrival velocity `Φvr r0 + Φvv v0⁺` and `vf`.  For a given `tf` this
/// transfer is unique, so the pair is the only two-impulse solution.
///
/// Returns `Some((dv1, dv2))`, or `None` when the problem is singular or
/// ill-posed: a non-positive or non-finite `tf` or `n`, non-finite states, a
/// singular in-plane block or a vanishing cross-track block.  With
/// `φ = n·tf` the in-plane determinant is `8 (1 - cos φ) - 3 φ sin φ`; it
/// vanishes at `φ = 2πk` (the chaser returns to its starting point, so no
/// in-plane burn can change the arrival position) and at the other roots of
/// `tan(φ/2) = 3φ/8`, the first of which is `φ ≈ 8.8387`.  The cross-track
/// block `sin φ / n` vanishes whenever `φ` is a multiple of `π`.
pub fn hcw_two_impulse(
    n: f64,
    r0: [f64; 3],
    v0: [f64; 3],
    rf: [f64; 3],
    vf: [f64; 3],
    tf: f64,
) -> Option<([f64; 3], [f64; 3])> {
    if !(n.is_finite() && n > 0.0 && tf.is_finite() && tf > 0.0) {
        return None;
    }
    if !(is_finite_vec(r0) && is_finite_vec(v0) && is_finite_vec(rf) && is_finite_vec(vf)) {
        return None;
    }

    let phi = n * tf;
    let (s, c) = phi.sin_cos();

    // Position that the departure velocity has to make up: rf - Φrr·r0.
    let wx = rf[0] - (4.0 - 3.0 * c) * r0[0];
    let wy = rf[1] - (6.0 * (s - phi) * r0[0] + r0[1]);
    let wz = rf[2] - c * r0[2];

    // In-plane block Φrv = [[a, b], [-b, d]] / n, with n·Φrv = [[s, 2(1-c)], [-2(1-c), 4s-3φ]].
    let a = s;
    let b = 2.0 * (1.0 - c);
    let d = 4.0 * s - 3.0 * phi;
    let det = a * d + b * b;
    let scale = a.abs().max(b.abs()).max(d.abs());
    if det.abs() <= SINGULARITY_TOL * scale * scale || s.abs() <= SINGULARITY_TOL {
        return None;
    }

    // v0⁺ = Φrv⁻¹ w;  Φrv⁻¹ = n/det · [[d, -b], [b, a]].
    let vx_plus = n * (d * wx - b * wy) / det;
    let vy_plus = n * (b * wx + a * wy) / det;
    let vz_plus = n * wz / s;
    let v_plus = [vx_plus, vy_plus, vz_plus];

    // Arrival velocity on the coasting arc: Φvr·r0 + Φvv·v0⁺.
    let arrival = [
        3.0 * n * s * r0[0] + c * v_plus[0] + 2.0 * s * v_plus[1],
        6.0 * n * (c - 1.0) * r0[0] - 2.0 * s * v_plus[0] + (4.0 * c - 3.0) * v_plus[1],
        -n * s * r0[2] + c * v_plus[2],
    ];

    let dv1 = [v_plus[0] - v0[0], v_plus[1] - v0[1], v_plus[2] - v0[2]];
    let dv2 = [vf[0] - arrival[0], vf[1] - arrival[1], vf[2] - arrival[2]];
    if is_finite_vec(dv1) && is_finite_vec(dv2) {
        Some((dv1, dv2))
    } else {
        None
    }
}
