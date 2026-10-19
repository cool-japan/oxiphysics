// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! D1Q3 lattice kernel shared by the multiscale LBM solvers.
//!
//! A node carries the three populations `[f_rest, f_plus, f_minus]` belonging
//! to the lattice velocities `c = [0, +1, -1]` (lattice units, `dx = dt = 1`)
//! with the weights `w = [2/3, 1/6, 1/6]` and the lattice speed of sound
//! `cs² = 1/3`.  The moments used throughout the module are
//!
//! - `Σ f_i` (zeroth moment) and
//! - `Σ c_i f_i = f_plus - f_minus` (first moment).
//!
//! Two equilibria are provided:
//!
//! - [`hydrodynamic_equilibrium`]: isothermal, momentum conserving,
//!   `f_i^eq = w_i ρ (1 + 3 c_i u + 4.5 c_i² u² - 1.5 u²)`, whose second
//!   moment is `ρ (cs² + u²)`;
//! - [`advection_equilibrium`]: linear in the transported scalar,
//!   `g_i^eq = w_i φ (1 + 3 c_i u)`, for passive advection-diffusion with the
//!   prescribed velocity `u`.
//!
//! Collision is BGK, `f_i ← f_i - ω (f_i - f_i^eq)`, and the Chapman-Enskog
//! diffusivity of either model is `D = cs² (τ - 1/2)` with `τ = 1/ω`.

use crate::error::{Error, Result};

/// Number of lattice velocities.
pub(super) const Q: usize = 3;

/// Lattice velocities in lattice units, ordered rest, `+x`, `-x`.
pub(super) const VELOCITIES: [f64; Q] = [0.0, 1.0, -1.0];

/// Lattice weights, ordered like [`VELOCITIES`].
pub(super) const WEIGHTS: [f64; Q] = [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0];

/// Squared lattice speed of sound.
pub(super) const CS2: f64 = 1.0 / 3.0;

/// Zeroth moment `Σ f_i` of a node.
///
/// The two moving populations are added first so that the rest state
/// `(2/3, 1/6, 1/6)` sums to exactly `1` in floating point.
#[inline]
pub(super) fn zeroth_moment(f: &[f64; Q]) -> f64 {
    f[0] + (f[1] + f[2])
}

/// First moment `Σ c_i f_i` of a node.
#[inline]
pub(super) fn first_moment(f: &[f64; Q]) -> f64 {
    f[1] - f[2]
}

/// Velocity `u = (Σ c_i f_i) / ρ` of a node, `0` where `ρ` is not positive.
#[inline]
pub(super) fn velocity(f: &[f64; Q]) -> f64 {
    let rho = zeroth_moment(f);
    if rho > 0.0 {
        first_moment(f) / rho
    } else {
        0.0
    }
}

/// Isothermal momentum-conserving equilibrium of the density `rho` moving with
/// the velocity `u`.
///
/// `f_i^eq = w_i ρ (1 + c_i u/cs² + (c_i u)²/(2 cs⁴) - u²/(2 cs²))`.
///
/// The moving populations follow the formula; the rest population is its exact
/// complement `ρ - (f_+ + f_-) = (2/3) ρ (1 - 3u²/2)`.  The double-precision
/// weights do not sum to one (`2/3 + 1/6 + 1/6 = 1 - 5.6e-17`), so evaluating
/// the rest population from its own weight would remove a fixed fraction of the
/// mass at every collision; the complement keeps `Σ f_i^eq = ρ` to within one
/// rounding error.
pub(super) fn hydrodynamic_equilibrium(rho: f64, u: f64) -> [f64; Q] {
    let u2 = u * u;
    let moving = |c: f64| {
        let cu = c * u;
        WEIGHTS[1] * rho * (1.0 + cu / CS2 + cu * cu / (2.0 * CS2 * CS2) - u2 / (2.0 * CS2))
    };
    let plus = moving(VELOCITIES[1]);
    let minus = moving(VELOCITIES[2]);
    [rho - (plus + minus), plus, minus]
}

/// Equilibrium of a passive scalar `phi` advected with the velocity `u`.
///
/// `g_i^eq = w_i φ (1 + c_i u/cs²)`; it is linear in `phi`, so
/// `advection_equilibrium(a, u) + advection_equilibrium(b, u)` is
/// `advection_equilibrium(a + b, u)`.  As for
/// [`hydrodynamic_equilibrium`] the rest population is the exact complement
/// `φ - (g_+ + g_-)` so that `Σ g_i^eq = φ` to within one rounding error.
pub(super) fn advection_equilibrium(phi: f64, u: f64) -> [f64; Q] {
    let moving = |c: f64| WEIGHTS[1] * phi * (1.0 + c * u / CS2);
    let plus = moving(VELOCITIES[1]);
    let minus = moving(VELOCITIES[2]);
    [phi - (plus + minus), plus, minus]
}

/// BGK collision `f ← f - ω (f - f_eq)` of one node.
#[inline]
pub(super) fn bgk_relax(f: &mut [f64; Q], feq: &[f64; Q], omega: f64) {
    for (fi, &eq) in f.iter_mut().zip(feq.iter()) {
        *fi -= omega * (*fi - eq);
    }
}

/// Non-equilibrium part `f - f_eq` of a node.
#[inline]
pub(super) fn nonequilibrium(f: &[f64; Q], feq: &[f64; Q]) -> [f64; Q] {
    [f[0] - feq[0], f[1] - feq[1], f[2] - feq[2]]
}

/// Periodic push streaming of a row of nodes, in place.
///
/// `f_plus` moves one node towards `+x` and `f_minus` one node towards `-x`;
/// the populations leaving through one end re-enter through the other.
pub(super) fn stream_periodic(nodes: &mut [[f64; Q]]) {
    let n = nodes.len();
    if n < 2 {
        return;
    }
    let wrapped_plus = nodes[n - 1][1];
    for i in (1..n).rev() {
        nodes[i][1] = nodes[i - 1][1];
    }
    nodes[0][1] = wrapped_plus;
    let wrapped_minus = nodes[0][2];
    for i in 0..n - 1 {
        nodes[i][2] = nodes[i + 1][2];
    }
    nodes[n - 1][2] = wrapped_minus;
}

/// Open (non-periodic) push streaming of a row of nodes, in place.
///
/// The two populations that would enter through the ends are not defined by
/// the streaming step and are set to zero; the caller must reconstruct them
/// (boundary condition).
pub(super) fn stream_open(nodes: &mut [[f64; Q]]) {
    let n = nodes.len();
    if n == 0 {
        return;
    }
    for i in (1..n).rev() {
        nodes[i][1] = nodes[i - 1][1];
    }
    nodes[0][1] = 0.0;
    for i in 0..n - 1 {
        nodes[i][2] = nodes[i + 1][2];
    }
    nodes[n - 1][2] = 0.0;
}

/// Check that the BGK relaxation rate `omega` is finite and gives a stable
/// relaxation time, i.e. `0 < ω < 2` (`τ = 1/ω > 1/2`).  `name` identifies the
/// parameter in the error.
pub(super) fn validate_relaxation_rate(omega: f64, name: &str) -> Result<()> {
    if !omega.is_finite() || omega <= 0.0 {
        return Err(Error::ParameterOutOfBounds {
            name: name.to_string(),
            value: omega,
            min: 0.0,
            max: 2.0,
        });
    }
    Error::check_tau(1.0 / omega)
}
