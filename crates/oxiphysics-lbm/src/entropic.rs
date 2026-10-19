// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Entropic Lattice Boltzmann Method (H-theorem stabilized) for D2Q9.
//!
//! This module implements two entropic LBM variants:
//!
//! - **EntropicD2Q9**: Standard entropic LBM with adaptive relaxation parameter α
//!   determined by enforcing the discrete H-theorem (H-stability).
//! - **KbcD2Q9**: Karlin–Bösch–Chikatamarla (KBC) multi-relaxation variant that
//!   splits the non-equilibrium part into shear (s) and higher-order (h)
//!   natural moments, relaxes the shear part with the requested viscosity and
//!   the higher-order part with an entropic stabiliser `γ` (see [`kbc_split`]
//!   and [`kbc_stabiliser`]).
//!
//! Also provides:
//! - ELBM diagnostics (entropy monitoring, alpha statistics)
//! - Newton-Raphson alpha optimization (faster convergence)
//! - Entropy function variants (Boltzmann H, relative entropy)
//!
//! ## Relaxation and viscosity
//!
//! The entropic collisions ([`EntropicD2Q9`], [`EntropyOverRelaxation`]) relax
//! the populations towards the local
//! equilibrium `f_eq` as `f' = f + α β (f_eq − f)` with the physical
//! relaxation parameter `β = 1/(2τ) = 1/(6ν + 1)` (`τ = 3ν + 1/2` for the
//! kinematic viscosity `ν` in lattice units).  The entropic parameter `α` is
//! the non-trivial root of `H(f + α (f_eq − f)) = H(f)`; for a state close to
//! equilibrium `α = 2 − O(|f − f_eq|)`, so the scheme reduces to BGK with
//! `ω = 2β = 1/τ` and reproduces exactly the requested viscosity `ν`.  Far
//! from equilibrium `α` adapts so that the H-theorem `H(f') ≤ H(f)` holds.
//! The KBC collision ([`KbcD2Q9`]) uses the same `β` for its shear moments,
//! `f' = f − β (2 Δs + γ Δh)`, so it too reproduces the requested shear
//! viscosity.
//!
//! The collision target is the *entropic* equilibrium, i.e. the minimiser of
//! `H` for the cell's density and momentum (see [`d2q9_equilibrium`]); with a
//! merely polynomial equilibrium the entropy estimate would be dominated by
//! the polynomial's O(u³) misfit instead of by the non-equilibrium part.
//!
//! ## References
//!
//! - Ansumali & Karlin (2002), Phys. Rev. E 65, 056312.
//! - Karlin, Bösch & Chikatamarla (2014), Phys. Rev. E 90, 031302.

use crate::error::{Error, Result};

mod kbc;

pub use kbc::{
    KBC_H_NOISE_FLOOR, entropic_scalar_product, kbc_collide_populations, kbc_shear_part, kbc_split,
    kbc_stabiliser,
};

// ---------------------------------------------------------------------------
// D2Q9 constants
// ---------------------------------------------------------------------------

/// D2Q9 lattice weights: \[centre, axis×4, diagonal×4\]
///
/// Ordering convention (standard):
/// ```text
/// 0: (0,0)  w = 4/9
/// 1: (1,0)  w = 1/9
/// 2: (0,1)  w = 1/9
/// 3: (-1,0) w = 1/9
/// 4: (0,-1) w = 1/9
/// 5: (1,1)  w = 1/36
/// 6: (-1,1) w = 1/36
/// 7: (-1,-1)w = 1/36
/// 8: (1,-1) w = 1/36
/// ```
const W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];

/// D2Q9 lattice velocities: (cx, cy) for each direction.
const C: [(i32, i32); 9] = [
    (0, 0),
    (1, 0),
    (0, 1),
    (-1, 0),
    (0, -1),
    (1, 1),
    (-1, 1),
    (-1, -1),
    (1, -1),
];

/// D2Q9 opposite-direction mapping.
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Relative non-equilibrium deviation `max_i |f_i − f_eq,i| / f_eq,i` below
/// which the entropic root search is skipped and `α = 2` (BGK
/// over-relaxation) is returned.
///
/// Near equilibrium the exact root is `α = 2 − O(δ)` for a deviation `δ`, so
/// `α = 2` is accurate to `O(δ)`; at the same time the entropy change
/// `φ(α) = O(δ²)` drops below what double-precision evaluation of `H` can
/// resolve for `δ ≲ 1e-6`.
const NEAR_EQUILIBRIUM_REL_TOL: f64 = 1e-5;

/// Maximum number of geometric (×1.5) growth steps when bracketing the
/// entropic root from above.
const MAX_BRACKET_GROWTH: usize = 64;

/// Tolerance on the entropy change `H(f') − H(f)` accepted by the Newton
/// solvers when they land on the root from the wrong side.
const ENTROPY_TOL: f64 = 1e-12;

/// Relative back-off from the positivity boundary of the mirror line (see
/// [`back_off_positivity`]).
const POSITIVITY_MARGIN: f64 = 1e-6;

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the second-order polynomial (Taylor-expanded Maxwell–Boltzmann)
/// D2Q9 equilibrium:
///
/// f_eq_i = w_i * ρ * (1 + (c·u)/cs² + (c·u)²/(2cs⁴) − u²/(2cs²)),  cs² = 1/3.
///
/// This is the standard LBM equilibrium.  It is *not* the minimiser of the
/// H-function (its misfit is O(u³)), which is why the entropic collisions of
/// this module use [`d2q9_equilibrium`] instead.
pub fn d2q9_polynomial_equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
    let ux = u[0];
    let uy = u[1];
    let u2 = ux * ux + uy * uy;
    let mut feq = [0.0; 9];
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        let cu = cx * ux + cy * uy;
        feq[i] = W[i] * rho * (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
    }
    feq
}

/// Compute the D2Q9 *entropic* equilibrium: the minimiser of
/// `H = Σ f_i ln(f_i / w_i)` at fixed density `ρ` and momentum `ρ u`.
///
/// For D2Q9 the minimiser has the closed product form of Ansumali, Karlin &
/// Öttinger (2003),
///
/// ```text
/// f_eq_i = ρ w_i Π_{a ∈ {x,y}} (2 − √(1+3u_a²)) · ((2u_a + √(1+3u_a²)) / (1 − u_a))^{c_ia}
/// ```
///
/// It conserves mass and momentum exactly and agrees with
/// [`d2q9_polynomial_equilibrium`] up to `O(u³)`.  Because `H(f_eq)` is the
/// minimum of `H` over all populations with the same `ρ` and `ρ u`, the
/// entropic root `α` of `H(f + α (f_eq − f)) = H(f)` equals `2 − O(|f − f_eq|)`
/// for every non-equilibrium magnitude.
///
/// The product form is only defined for `|u_a| < 1`; outside that range (far
/// beyond the validity of any LBM) the polynomial equilibrium is returned.
pub fn d2q9_equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
    if !(u[0].abs() < 1.0 && u[1].abs() < 1.0) {
        return d2q9_polynomial_equilibrium(rho, u);
    }
    // One-dimensional factors for c = -1, 0, +1.
    let factors = |ua: f64| -> [f64; 3] {
        let s = (1.0 + 3.0 * ua * ua).sqrt();
        let xi = 2.0 - s;
        let phi = (2.0 * ua + s) / (1.0 - ua);
        [xi / phi, xi, xi * phi]
    };
    let fx = factors(u[0]);
    let fy = factors(u[1]);
    let mut feq = [0.0; 9];
    for i in 0..9 {
        let ix = (C[i].0 + 1) as usize;
        let iy = (C[i].1 + 1) as usize;
        feq[i] = W[i] * rho * fx[ix] * fy[iy];
    }
    feq
}

/// Relaxation parameter `β = 1/(2τ) = 1/(6ν + 1)` of an entropic collision for
/// the kinematic viscosity `ν` (lattice units, `τ = 3ν + 1/2`).
///
/// With `α = 2` the collision `f' = f + α β (f_eq − f)` is exactly BGK with
/// `ω = 1/τ`.  A non-positive or non-finite `ν` is outside the physical range
/// (see the `try_new` constructors, which reject it) and is treated as the
/// inviscid limit `β = 1`.
pub fn relaxation_beta(nu: f64) -> f64 {
    if nu.is_finite() && nu > 0.0 {
        1.0 / (6.0 * nu + 1.0)
    } else {
        1.0
    }
}

/// Validate a kinematic viscosity for the entropic collisions: it must be
/// finite and positive, i.e. `τ = 3ν + 1/2 > 1/2`.
fn validate_viscosity(nu: f64) -> Result<()> {
    Error::check_range(nu, 0.0, f64::MAX, "nu")?;
    Error::check_tau(3.0 * nu + 0.5)
}

/// Compute the discrete H-function: H = Σ_i f_i * ln(f_i / w_i).
///
/// This is the negative of the Boltzmann entropy. The H-theorem guarantees
/// H decreases (or stays constant) during collision.
///
/// Returns `f64::INFINITY` if any `f_i ≤ 0` (unphysical state).
pub fn h_function(f: &[f64; 9]) -> f64 {
    let mut h = 0.0;
    for i in 0..9 {
        if f[i] <= 0.0 {
            return f64::INFINITY;
        }
        h += f[i] * (f[i] / W[i]).ln();
    }
    h
}

/// Compute the relative entropy (Kullback-Leibler divergence) of f from f_eq:
///
/// ```text
/// D_KL(f || f_eq) = Σ_i f_i * ln(f_i / f_eq_i)
/// ```
///
/// Returns `f64::INFINITY` if any `f_i ≤ 0` or `f_eq_i ≤ 0`.
pub fn relative_entropy(f: &[f64; 9], f_eq: &[f64; 9]) -> f64 {
    let mut dkl = 0.0;
    for i in 0..9 {
        if f[i] <= 0.0 || f_eq[i] <= 0.0 {
            return f64::INFINITY;
        }
        dkl += f[i] * (f[i] / f_eq[i]).ln();
    }
    dkl
}

/// Compute the Boltzmann entropy: S = -Σ_i f_i * ln(f_i / w_i) = -H.
pub fn boltzmann_entropy(f: &[f64; 9]) -> f64 {
    -h_function(f)
}

/// Compute the mirror state: f_mirror = 2 * f_eq − f.
///
/// The mirror state is used in the entropic overrelaxation to find the
/// α value that keeps the post-collision state on the H = H(f) level set.
pub fn mirror_state(f: &[f64; 9], f_eq: &[f64; 9]) -> [f64; 9] {
    let mut mirror = [0.0; 9];
    for i in 0..9 {
        mirror[i] = 2.0 * f_eq[i] - f[i];
    }
    mirror
}

/// Deviation `d = f_eq − f` along which the entropic mirror line runs.
fn mirror_direction(f: &[f64; 9], f_eq: &[f64; 9]) -> [f64; 9] {
    let mut d = [0.0; 9];
    for i in 0..9 {
        d[i] = f_eq[i] - f[i];
    }
    d
}

/// `true` if every population is within [`NEAR_EQUILIBRIUM_REL_TOL`] (relative)
/// of the equilibrium.
fn is_near_equilibrium(f: &[f64; 9], f_eq: &[f64; 9]) -> bool {
    f.iter()
        .zip(f_eq.iter())
        .all(|(&fi, &ei)| ei > 0.0 && (fi - ei).abs() <= NEAR_EQUILIBRIUM_REL_TOL * ei)
}

/// Entropy change `φ(α) = H(f + α d) − H(f)` along the mirror line.
///
/// Evaluated term by term as `Σ_i [α d_i ln((f_i + α d_i)/w_i) + f_i ln(1 + α d_i/f_i)]`,
/// which is algebraically identical to the difference of the two H values but
/// avoids the catastrophic cancellation of subtracting two O(1) numbers whose
/// difference is O(|d|²).  Returns `+∞` as soon as a population `f_i + α d_i`
/// (or `f_i`) is not strictly positive, where `H` is undefined.
fn entropy_change(f: &[f64; 9], d: &[f64; 9], alpha: f64) -> f64 {
    let mut total = 0.0;
    for i in 0..9 {
        let x = f[i] + alpha * d[i];
        let admissible = f[i] > 0.0 && x > 0.0;
        if !admissible {
            return f64::INFINITY;
        }
        total += alpha * d[i] * (x / W[i]).ln() + f[i] * (alpha * d[i] / f[i]).ln_1p();
    }
    total
}

/// Derivative `dφ/dα = Σ_i d_i (1 + ln((f_i + α d_i)/w_i))` of
/// [`entropy_change`]; `+∞` where `H` is undefined.
fn entropy_change_derivative(f: &[f64; 9], d: &[f64; 9], alpha: f64) -> f64 {
    let mut total = 0.0;
    for i in 0..9 {
        let x = f[i] + alpha * d[i];
        if x > 0.0 {
            total += d[i] * (1.0 + (x / W[i]).ln());
        } else {
            return f64::INFINITY;
        }
    }
    total
}

/// Bracket the non-trivial root of `φ(α) = 0`.
///
/// `φ` is convex with `φ(0) = 0`, so a root `α* > 0` exists iff `φ(1) < 0`
/// (the equilibrium has lower entropy than the pre-collision state) and then
/// `α* > 1`.  Starting from `[1, 2]` the upper end is grown geometrically
/// until `φ(hi) ≥ 0`; `φ = +∞` once a population would turn non-positive, so
/// the growth always terminates.  Returns `None` when there is no such root.
fn bracket_entropic_root(f: &[f64; 9], d: &[f64; 9]) -> Option<(f64, f64)> {
    if entropy_change(f, d, 1.0) < 0.0 {
        let mut lo = 1.0_f64;
        let mut hi = 2.0_f64;
        for _ in 0..MAX_BRACKET_GROWTH {
            if entropy_change(f, d, hi) < 0.0 {
                lo = hi;
                hi *= 1.5;
            } else {
                return Some((lo, hi));
            }
        }
    }
    None
}

/// Back off from the positivity boundary of the mirror line.
///
/// For extremely strong non-equilibrium the line `f + α d` leaves the positive
/// orthant before the entropy returns to `H(f)`; there is then no iso-entropy
/// root and the largest admissible `α` sits on the boundary, where a
/// population would be exactly zero (and `H` undefined at the next step).
/// If `alpha` is within a relative `1e-8` of that boundary it is reduced by
/// [`POSITIVITY_MARGIN`], which keeps every over-relaxed population at least
/// that fraction of its pre-collision value.  The H-theorem is unaffected
/// because `φ < 0` on the whole admissible range.
fn back_off_positivity(f: &[f64; 9], d: &[f64; 9], alpha: f64) -> f64 {
    if entropy_change(f, d, alpha * (1.0 + 1e-8)).is_finite() {
        alpha
    } else {
        alpha * (1.0 - POSITIVITY_MARGIN)
    }
}

/// Safeguarded Newton iteration for `φ(α) = 0` inside the bracket `[lo, hi]`
/// (`φ(lo) < 0 ≤ φ(hi)`), started at `alpha0`.
///
/// A Newton step that leaves the current bracket (or hits a vanishing /
/// infinite derivative) is replaced by bisection.  The result is the root
/// estimate if it satisfies the H-theorem to [`ENTROPY_TOL`], otherwise the
/// safe lower end of the bracket.
fn safeguarded_newton(
    f: &[f64; 9],
    d: &[f64; 9],
    bracket: (f64, f64),
    alpha0: f64,
    tol: f64,
    max_iter: usize,
) -> f64 {
    let (mut lo, mut hi) = bracket;
    let mut alpha = if alpha0 > lo && alpha0 < hi {
        alpha0
    } else {
        0.5 * (lo + hi)
    };
    for _ in 0..max_iter {
        let phi = entropy_change(f, d, alpha);
        if phi < 0.0 {
            lo = alpha;
        } else {
            hi = alpha;
        }
        let slope = entropy_change_derivative(f, d, alpha);
        let newton = if phi.is_finite() && slope.is_finite() && slope != 0.0 {
            alpha - phi / slope
        } else {
            f64::NAN
        };
        let next = if newton > lo && newton < hi {
            newton
        } else {
            0.5 * (lo + hi)
        };
        let converged = (next - alpha).abs() <= tol * alpha.max(1.0) || hi - lo <= tol;
        alpha = next;
        if converged {
            break;
        }
    }
    let root = if entropy_change(f, d, alpha) <= ENTROPY_TOL {
        alpha
    } else {
        lo
    };
    back_off_positivity(f, d, root)
}

/// Find the entropic relaxation parameter α via bracketed bisection.
///
/// With `d = f_eq − f` the mirror line is `f + α d`, and α is the non-trivial
/// root of
///
/// ```text
/// φ(α) = H(f + α d) − H(f) = 0,
/// ```
///
/// i.e. the over-relaxed state has exactly the entropy of the pre-collision
/// state and the discrete H-theorem holds for every state between `f` and
/// `f + α d`.  `φ` is convex with `φ(0) = 0` and `φ(1) = H(f_eq) − H(f)`, so
/// for `H(f_eq) < H(f)` the root lies above 1; it equals `2 − O(|d|)` near
/// equilibrium and is below or above 2 for strongly non-equilibrium states.
/// The returned value is the lower end of the final bracket, so
/// `φ(α) ≤ 0` always holds.  In the extreme case where the mirror line leaves
/// the positive orthant before the entropy returns to `H(f)` (populations
/// would turn negative first), no iso-entropy root exists and the result is
/// the largest positivity-preserving value, backed off by a relative `1e-6`.
///
/// `α = 2` (standard BGK over-relaxation) is returned when
/// - the state is within a relative `1e-5` of
///   equilibrium, where the root is `2 − O(δ)` and cannot be resolved in
///   floating point, or
/// - `H(f_eq) ≥ H(f)` or `H` is undefined (non-positive populations), where
///   no entropy-based correction exists.
///
/// # Arguments
/// * `f`        – pre-collision distribution
/// * `f_eq`     – local equilibrium distribution (use [`d2q9_equilibrium`])
/// * `max_iter` – maximum bisection iterations (20–50 is typically sufficient)
pub fn find_entropic_alpha(f: &[f64; 9], f_eq: &[f64; 9], max_iter: usize) -> f64 {
    if is_near_equilibrium(f, f_eq) {
        return 2.0;
    }
    let d = mirror_direction(f, f_eq);
    let Some((mut lo, mut hi)) = bracket_entropic_root(f, &d) else {
        return 2.0;
    };
    for _ in 0..max_iter {
        let mid = 0.5 * (lo + hi);
        if entropy_change(f, &d, mid) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-14 {
            break;
        }
    }
    back_off_positivity(f, &d, lo)
}

/// Find the entropic α using safeguarded Newton–Raphson iteration (faster
/// convergence than [`find_entropic_alpha`]).
///
/// Solves the same equation `H(f + α (f_eq − f)) = H(f)` using the analytic
/// derivative `dH/dα = Σ_i (f_eq_i − f_i)(1 + ln((f_i + α (f_eq_i − f_i))/w_i))`.
/// The iteration starts at `α = 2`, is kept inside the bracket of the root
/// (a Newton step leaving the bracket falls back to bisection) and is not
/// clamped to `[0.5, 2]`, so roots above 2 are found as well.  The same
/// `α = 2` fall-backs as in [`find_entropic_alpha`] apply.
pub fn find_entropic_alpha_newton(f: &[f64; 9], f_eq: &[f64; 9], max_iter: usize) -> f64 {
    if is_near_equilibrium(f, f_eq) {
        return 2.0;
    }
    let d = mirror_direction(f, f_eq);
    let Some(bracket) = bracket_entropic_root(f, &d) else {
        return 2.0;
    };
    safeguarded_newton(f, &d, bracket, 2.0, 1e-14, max_iter)
}

/// Compute non-equilibrium stress tensor components from f and f_eq.
///
/// Returns the three independent components of the symmetric
/// non-equilibrium stress tensor in 2D: `(σ_xx, σ_xy, σ_yy)`.
pub fn non_equilibrium_stress(f: &[f64; 9], f_eq: &[f64; 9]) -> (f64, f64, f64) {
    let mut s_xx = 0.0;
    let mut s_xy = 0.0;
    let mut s_yy = 0.0;
    for i in 0..9 {
        let cx = C[i].0 as f64;
        let cy = C[i].1 as f64;
        let f_neq = f[i] - f_eq[i];
        s_xx += cx * cx * f_neq;
        s_xy += cx * cy * f_neq;
        s_yy += cy * cy * f_neq;
    }
    (s_xx, s_xy, s_yy)
}

// ---------------------------------------------------------------------------
// EntropicD2Q9
// ---------------------------------------------------------------------------

/// Entropic LBM solver for D2Q9 with adaptive relaxation (H-theorem stable).
///
/// At every cell the collision is `f' = f + α β (f_eq − f)` with the physical
/// relaxation parameter `β = 1/(6ν + 1)` of the kinematic viscosity `nu` and
/// the entropic parameter `α` of [`find_entropic_alpha`], so that the
/// post-collision state never violates the discrete H-theorem.  Near
/// equilibrium `α = 2` and the scheme is exactly BGK with `τ = 3ν + 1/2`.
pub struct EntropicD2Q9 {
    /// Domain width (number of cells in x direction).
    pub nx: usize,
    /// Domain height (number of cells in y direction).
    pub ny: usize,
    /// Distribution functions stored as `f[cell_index][direction]`.
    pub f: Vec<[f64; 9]>,
    /// Kinematic viscosity ν (lattice units); must be positive.
    pub nu: f64,
}

impl EntropicD2Q9 {
    /// Create a new solver with uniform rest equilibrium (ρ=1, u=0).
    ///
    /// `nu` is the kinematic viscosity in lattice units and must be positive
    /// (`τ = 3ν + 1/2 > 1/2`); use [`try_new`](Self::try_new) to have an
    /// invalid value reported instead of being treated as the inviscid limit.
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let f_rest = Self::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![f_rest; n],
            nu,
        }
    }

    /// Validated constructor: rejects an empty grid and a viscosity that is
    /// not finite and positive (`τ = 3ν + 1/2 ≤ 1/2`, checked with
    /// [`Error::check_tau`]).
    pub fn try_new(nx: usize, ny: usize, nu: f64) -> Result<Self> {
        Error::check_grid(nx, ny)?;
        validate_viscosity(nu)?;
        Ok(Self::new(nx, ny, nu))
    }

    /// Relaxation parameter `β = 1/(2τ) = 1/(6ν + 1)` of this solver.
    pub fn beta(&self) -> f64 {
        relaxation_beta(self.nu)
    }

    /// Compute the D2Q9 entropic (H-minimising) equilibrium distribution;
    /// see [`d2q9_equilibrium`].
    pub fn equilibrium(rho: f64, u: [f64; 2]) -> [f64; 9] {
        d2q9_equilibrium(rho, u)
    }

    /// Compute macroscopic density and velocity at cell `idx`.
    ///
    /// Returns `(rho, [ux, uy])`.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let f = &self.f[idx];
        let rho: f64 = f.iter().sum();
        let mut mx = 0.0;
        let mut my = 0.0;
        for i in 0..9 {
            mx += f[i] * C[i].0 as f64;
            my += f[i] * C[i].1 as f64;
        }
        let inv_rho = if rho > 0.0 { 1.0 / rho } else { 0.0 };
        (rho, [mx * inv_rho, my * inv_rho])
    }

    /// Perform entropic collision at every cell.
    ///
    /// For each cell the adaptive α is found by bisection such that
    /// H(f_new) ≤ H(f).  The update rule is `f_new = f + α β (f_eq − f)` with
    /// `β = 1/(6ν + 1)`; α is capped at `2/β` so that the effective
    /// relaxation frequency `α β` never exceeds the linear stability limit 2.
    pub fn collide(&mut self) {
        self.collide_with(|f, f_eq| find_entropic_alpha(f, f_eq, 50));
    }

    /// Perform entropic collision using Newton-Raphson alpha optimization.
    pub fn collide_newton(&mut self) {
        self.collide_with(|f, f_eq| find_entropic_alpha_newton(f, f_eq, 20));
    }

    /// Shared collision loop: `alpha_of(f, f_eq)` supplies the entropic α.
    fn collide_with(&mut self, alpha_of: impl Fn(&[f64; 9], &[f64; 9]) -> f64) {
        let beta = self.beta();
        let alpha_cap = 2.0 / beta;
        let n = self.nx * self.ny;
        for k in 0..n {
            let (rho, u) = self.macros(k);
            let f_eq = Self::equilibrium(rho, u);
            let alpha = alpha_of(&self.f[k], &f_eq).min(alpha_cap);
            let f = &mut self.f[k];
            for i in 0..9 {
                f[i] += alpha * beta * (f_eq[i] - f[i]);
            }
        }
    }

    /// Perform pull-scheme streaming with periodic boundary conditions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k_dst = y * nx + x;
                for i in 0..9 {
                    let cx = C[i].0;
                    let cy = C[i].1;
                    // Pull from the source cell (periodic wrap).
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let k_src = sy * nx + sx;
                    self.f[k_dst][i] = f_old[k_src][i];
                }
            }
        }
    }

    /// Perform one full LBM step: collide then stream.
    pub fn step(&mut self) {
        self.collide();
        self.stream();
    }

    /// Initialize all cells to a uniform density and velocity.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = Self::equilibrium(rho, u);
        for cell in self.f.iter_mut() {
            *cell = feq;
        }
    }

    /// Compute the total H-function (entropy) over the entire domain.
    pub fn total_entropy(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut total_h = 0.0;
        for k in 0..n {
            total_h += h_function(&self.f[k]);
        }
        total_h
    }

    /// Compute the maximum velocity magnitude in the domain.
    pub fn max_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        let mut max_u = 0.0_f64;
        for k in 0..n {
            let (_, u) = self.macros(k);
            let u_mag = (u[0] * u[0] + u[1] * u[1]).sqrt();
            max_u = max_u.max(u_mag);
        }
        max_u
    }

    /// Compute the minimum and maximum density in the domain.
    pub fn density_range(&self) -> (f64, f64) {
        let n = self.nx * self.ny;
        let mut rho_min = f64::MAX;
        let mut rho_max = f64::MIN;
        for k in 0..n {
            let (rho, _) = self.macros(k);
            rho_min = rho_min.min(rho);
            rho_max = rho_max.max(rho);
        }
        (rho_min, rho_max)
    }
}

// ---------------------------------------------------------------------------
// ElbmDiagnostics
// ---------------------------------------------------------------------------

/// Diagnostics for the entropic LBM.
///
/// Tracks entropy evolution, alpha statistics, and stability indicators.
pub struct ElbmDiagnostics {
    /// History of total domain entropy (H function).
    pub entropy_history: Vec<f64>,
    /// History of minimum α values per step.
    pub alpha_min_history: Vec<f64>,
    /// History of average α values per step.
    pub alpha_avg_history: Vec<f64>,
    /// History of maximum velocity magnitudes per step.
    pub max_velocity_history: Vec<f64>,
}

impl ElbmDiagnostics {
    /// Create empty diagnostics.
    pub fn new() -> Self {
        Self {
            entropy_history: Vec::new(),
            alpha_min_history: Vec::new(),
            alpha_avg_history: Vec::new(),
            max_velocity_history: Vec::new(),
        }
    }
}

impl Default for ElbmDiagnostics {
    fn default() -> Self {
        Self::new()
    }
}

impl ElbmDiagnostics {
    /// Record diagnostics for one step from an EntropicD2Q9 solver.
    pub fn record(&mut self, solver: &EntropicD2Q9) {
        let n = solver.nx * solver.ny;

        // Total entropy
        let mut total_h = 0.0;
        let mut alpha_min = 2.0_f64;
        let mut alpha_sum = 0.0;
        let mut max_u = 0.0_f64;

        for k in 0..n {
            let (rho, u) = solver.macros(k);
            let f_eq = EntropicD2Q9::equilibrium(rho, u);
            let alpha = find_entropic_alpha(&solver.f[k], &f_eq, 50);
            total_h += h_function(&solver.f[k]);
            alpha_min = alpha_min.min(alpha);
            alpha_sum += alpha;
            let u_mag = (u[0] * u[0] + u[1] * u[1]).sqrt();
            max_u = max_u.max(u_mag);
        }

        self.entropy_history.push(total_h);
        self.alpha_min_history.push(alpha_min);
        self.alpha_avg_history.push(alpha_sum / n as f64);
        self.max_velocity_history.push(max_u);
    }

    /// Check if entropy is monotonically non-increasing (H-theorem).
    pub fn is_h_theorem_satisfied(&self) -> bool {
        if self.entropy_history.len() < 2 {
            return true;
        }
        for i in 1..self.entropy_history.len() {
            if self.entropy_history[i] > self.entropy_history[i - 1] + 1e-10 {
                return false;
            }
        }
        true
    }

    /// Number of recorded steps.
    pub fn step_count(&self) -> usize {
        self.entropy_history.len()
    }
}

// ---------------------------------------------------------------------------
// KbcD2Q9
// ---------------------------------------------------------------------------

/// Karlin–Bösch–Chikatamarla (KBC) entropic multi-relaxation LBM for D2Q9.
///
/// The non-equilibrium deviation `Δ = f − f_eq` is split in the natural-moment
/// basis into its shear part `Δs` (deviatoric stress `N = M_20 − M_02` and
/// `Π_xy`) and its higher-order part `Δh` (trace `T`, third-order moments
/// `Q_xyy`, `Q_yxx` and fourth-order moment `A`); see [`kbc_split`].  The
/// collision is
///
/// ```text
/// f' = f − β (2 Δs + γ Δh),   β = 1/(2τ) = 1/(6ν + 1),
/// γ  = 1/β − (2 − 1/β) ⟨Δs|Δh⟩ / ⟨Δh|Δh⟩,   ⟨x|y⟩ = Σ_i x_i y_i / f_eq,i.
/// ```
///
/// The shear moments are relaxed with the BGK frequency `1/τ`, so the shear
/// viscosity is exactly the requested `nu`; the stabiliser `γ` (see
/// [`kbc_stabiliser`]) minimises the quadratic entropy estimate of the
/// post-collision state over the relaxation of the higher-order moments,
/// which damps the ghost modes responsible for the instability of BGK at
/// `τ → 1/2`.  Since the trace `T` is part of `Δh`, the bulk viscosity is
/// governed by `γ` (bulk relaxation frequency `βγ`, equal to 1 for orthogonal
/// shear and higher-order parts).
pub struct KbcD2Q9 {
    /// Domain width.
    pub nx: usize,
    /// Domain height.
    pub ny: usize,
    /// Distribution functions stored as `f[cell_index][direction]`.
    pub f: Vec<[f64; 9]>,
    /// Kinematic viscosity ν (lattice units); must be positive.
    pub nu: f64,
}

impl KbcD2Q9 {
    /// Parity decomposition of `f_eq` into its even and odd parts under the
    /// velocity reversal `i ↔ opp(i)`:
    /// - s_i = (f_eq_i + f_eq_{opp(i)}) / 2   (even: density, stress, `A`)
    /// - h_i = f_eq_i − s_i                    (odd: momentum, third-order moments)
    ///
    /// Returns `(s, h)` where `s + h = f_eq`.  This is *not* the KBC
    /// natural-moment split used by [`collide_kbc`](Self::collide_kbc); that
    /// one is [`kbc_split`].
    pub fn decompose_feq(f_eq: &[f64; 9]) -> ([f64; 9], [f64; 9]) {
        let mut s = [0.0; 9];
        let mut h = [0.0; 9];
        for i in 0..9 {
            s[i] = 0.5 * (f_eq[i] + f_eq[OPP[i]]);
            h[i] = f_eq[i] - s[i];
        }
        (s, h)
    }

    /// Parity decomposition of an arbitrary distribution into its even and odd
    /// parts under `i ↔ opp(i)` (see [`decompose_feq`](Self::decompose_feq)).
    pub fn decompose_general(f: &[f64; 9]) -> ([f64; 9], [f64; 9]) {
        let mut s = [0.0; 9];
        let mut h = [0.0; 9];
        for i in 0..9 {
            s[i] = 0.5 * (f[i] + f[OPP[i]]);
            h[i] = f[i] - s[i];
        }
        (s, h)
    }

    /// Create a new KBC solver with uniform rest equilibrium (ρ=1, u=0).
    ///
    /// `nu` must be positive (`τ = 3ν + 1/2 > 1/2`); use
    /// [`try_new`](Self::try_new) to have an invalid value reported.
    pub fn new(nx: usize, ny: usize, nu: f64) -> Self {
        let n = nx * ny;
        let f_rest = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
        Self {
            nx,
            ny,
            f: vec![f_rest; n],
            nu,
        }
    }

    /// Validated constructor: rejects an empty grid and a viscosity that is
    /// not finite and positive (checked with [`Error::check_tau`]).
    pub fn try_new(nx: usize, ny: usize, nu: f64) -> Result<Self> {
        Error::check_grid(nx, ny)?;
        validate_viscosity(nu)?;
        Ok(Self::new(nx, ny, nu))
    }

    /// Relaxation parameter `β = 1/(2τ) = 1/(6ν + 1)` of this solver.
    pub fn beta(&self) -> f64 {
        relaxation_beta(self.nu)
    }

    /// Compute macroscopic density and velocity at cell `idx`.
    pub fn macros(&self, idx: usize) -> (f64, [f64; 2]) {
        let f = &self.f[idx];
        let rho: f64 = f.iter().sum();
        let mut mx = 0.0;
        let mut my = 0.0;
        for i in 0..9 {
            mx += f[i] * C[i].0 as f64;
            my += f[i] * C[i].1 as f64;
        }
        let inv_rho = if rho > 0.0 { 1.0 / rho } else { 0.0 };
        (rho, [mx * inv_rho, my * inv_rho])
    }

    /// Perform the KBC collision at every cell.
    ///
    /// With `Δ = f − f_eq` split into shear and higher-order parts
    /// (`Δ = Δs + Δh`, see [`kbc_split`]) the update is
    ///
    /// ```text
    /// f_new = f − β (2 Δs + γ Δh),   β = 1/(6ν + 1),
    /// ```
    ///
    /// where `γ` is the entropic stabiliser of [`kbc_stabiliser`].  The shear
    /// stress relaxes as `Π' − Π_eq = (1 − 1/τ)(Π − Π_eq)` with
    /// `τ = 3ν + 1/2`; density and momentum are conserved.
    pub fn collide_kbc(&mut self) {
        let beta = self.beta();
        let n = self.nx * self.ny;
        for k in 0..n {
            let (rho, u) = self.macros(k);
            let f_eq = EntropicD2Q9::equilibrium(rho, u);
            let (post, _gamma) = kbc_collide_populations(&self.f[k], &f_eq, beta);
            self.f[k] = post;
        }
    }

    /// Perform pull-scheme streaming with periodic boundary conditions.
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let f_old = self.f.clone();
        for y in 0..ny {
            for x in 0..nx {
                let k_dst = y * nx + x;
                for i in 0..9 {
                    let cx = C[i].0;
                    let cy = C[i].1;
                    let sx = (x as isize - cx as isize).rem_euclid(nx as isize) as usize;
                    let sy = (y as isize - cy as isize).rem_euclid(ny as isize) as usize;
                    let k_src = sy * nx + sx;
                    self.f[k_dst][i] = f_old[k_src][i];
                }
            }
        }
    }

    /// One full KBC step: collide then stream.
    pub fn step(&mut self) {
        self.collide_kbc();
        self.stream();
    }

    /// Set all cells to uniform density and velocity.
    pub fn set_uniform(&mut self, rho: f64, u: [f64; 2]) {
        let feq = EntropicD2Q9::equilibrium(rho, u);
        for cell in self.f.iter_mut() {
            *cell = feq;
        }
    }

    /// Compute total density in the domain.
    pub fn total_density(&self) -> f64 {
        let n = self.nx * self.ny;
        (0..n).map(|k| self.macros(k).0).sum()
    }

    /// Compute the maximum velocity magnitude in the domain.
    pub fn max_velocity(&self) -> f64 {
        let n = self.nx * self.ny;
        (0..n)
            .map(|k| {
                let (_, u) = self.macros(k);
                (u[0] * u[0] + u[1] * u[1]).sqrt()
            })
            .fold(0.0_f64, f64::max)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// H-theorem: H(f_eq) ≤ H(f) for perturbed f.
    #[test]
    fn test_h_theorem_at_equilibrium() {
        let rho = 1.0;
        let u = [0.05, 0.02];
        let f_eq = EntropicD2Q9::equilibrium(rho, u);

        // Perturb f away from equilibrium.
        let mut f_perturbed = f_eq;
        f_perturbed[1] += 0.02;
        f_perturbed[3] -= 0.01;
        f_perturbed[5] += 0.005;

        let h_eq = h_function(&f_eq);
        let h_f = h_function(&f_perturbed);

        assert!(
            h_eq <= h_f + 1e-12,
            "H-theorem violated: H(f_eq)={h_eq} > H(f)={h_f}"
        );
    }

    /// find_entropic_alpha returns a positive over-relaxation parameter whose
    /// mirror state does not raise the entropy (the root may exceed 2).
    #[test]
    fn test_find_entropic_alpha_range() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.03;
        f[3] -= 0.02;

        let alpha = find_entropic_alpha(&f, &f_eq, 50);
        assert!(alpha.is_finite() && alpha > 0.0, "alpha invalid: {alpha}");
        let mut mirrored = f;
        for i in 0..9 {
            mirrored[i] += alpha * (f_eq[i] - f[i]);
        }
        assert!(
            h_function(&mirrored) <= h_function(&f) + 1e-12,
            "alpha = {alpha} violates the H-theorem"
        );
    }

    /// find_entropic_alpha returns 2.0 when f is already at equilibrium.
    #[test]
    fn test_find_entropic_alpha_at_equilibrium() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
        let alpha = find_entropic_alpha(&f_eq, &f_eq, 50);
        assert!(
            (alpha - 2.0).abs() < 1e-10,
            "alpha at equilibrium should be 2.0, got {alpha}"
        );
    }

    /// After collide+stream on a periodic domain, total density must be conserved.
    #[test]
    fn test_density_conservation_after_step() {
        let nx = 8;
        let ny = 8;
        let nu = 1.0 / 6.0;
        let mut solver = EntropicD2Q9::new(nx, ny, nu);

        // Perturb a few cells.
        solver.f[5][1] += 0.02;
        solver.f[5][3] -= 0.01;
        solver.f[10][2] += 0.015;

        let rho_before: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();

        for _ in 0..10 {
            solver.step();
        }

        let rho_after: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();
        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "Density not conserved: before={rho_before}, after={rho_after}"
        );
    }

    /// Uniform flow initialized with set_uniform stays uniform after steps.
    #[test]
    fn test_uniform_flow_stays_uniform() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let rho0 = 1.0;
        let u0 = [0.05, 0.0];

        let mut solver = EntropicD2Q9::new(nx, ny, nu);
        solver.set_uniform(rho0, u0);

        for _ in 0..20 {
            solver.step();
        }

        let n = nx * ny;
        for k in 0..n {
            let (rho, u) = solver.macros(k);
            assert!(
                (rho - rho0).abs() < 1e-12,
                "Density changed at cell {k}: {rho} != {rho0}"
            );
            assert!(
                (u[0] - u0[0]).abs() < 1e-12,
                "ux changed at cell {k}: {} != {}",
                u[0],
                u0[0]
            );
            assert!(
                u[1].abs() < 1e-12,
                "uy became non-zero at cell {k}: {}",
                u[1]
            );
        }
    }

    // -----------------------------------------------------------------------
    // Newton-Raphson alpha: must return α ∈ [0, 2]
    // -----------------------------------------------------------------------
    #[test]
    fn test_newton_alpha_range() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.03;
        f[3] -= 0.02;

        let alpha = find_entropic_alpha_newton(&f, &f_eq, 20);
        assert!(
            alpha.is_finite() && alpha > 0.0,
            "Newton alpha invalid: {alpha}"
        );
        let mut mirrored = f;
        for i in 0..9 {
            mirrored[i] += alpha * (f_eq[i] - f[i]);
        }
        assert!(
            h_function(&mirrored) <= h_function(&f) + 1e-12,
            "Newton alpha = {alpha} violates the H-theorem"
        );
    }

    // -----------------------------------------------------------------------
    // Newton-Raphson: at equilibrium returns 2.0
    // -----------------------------------------------------------------------
    #[test]
    fn test_newton_alpha_at_equilibrium() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.0, 0.0]);
        let alpha = find_entropic_alpha_newton(&f_eq, &f_eq, 20);
        assert!(
            (alpha - 2.0).abs() < 1e-8,
            "Newton alpha at equilibrium should be ~2.0, got {alpha}"
        );
    }

    // -----------------------------------------------------------------------
    // Newton collision preserves density
    // -----------------------------------------------------------------------
    #[test]
    fn test_newton_collision_density_conservation() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let mut solver = EntropicD2Q9::new(nx, ny, nu);

        solver.f[5][1] += 0.01;
        solver.f[5][3] -= 0.005;

        let rho_before: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();
        solver.collide_newton();
        solver.stream();
        let rho_after: f64 = (0..nx * ny).map(|k| solver.macros(k).0).sum();

        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "Density not conserved with Newton: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // Relative entropy: zero at equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_relative_entropy_zero_at_eq() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let dkl = relative_entropy(&f_eq, &f_eq);
        assert!(
            dkl.abs() < 1e-14,
            "D_KL(f_eq || f_eq) should be 0, got {dkl}"
        );
    }

    // -----------------------------------------------------------------------
    // Relative entropy: positive for non-equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_relative_entropy_positive() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.02;
        f[3] -= 0.01;

        let dkl = relative_entropy(&f, &f_eq);
        assert!(
            dkl > 0.0,
            "D_KL(f || f_eq) should be positive for f != f_eq, got {dkl}"
        );
    }

    // -----------------------------------------------------------------------
    // Boltzmann entropy = -H
    // -----------------------------------------------------------------------
    #[test]
    fn test_boltzmann_entropy_is_neg_h() {
        let f = EntropicD2Q9::equilibrium(1.0, [0.03, -0.01]);
        let h = h_function(&f);
        let s = boltzmann_entropy(&f);
        assert!((s + h).abs() < 1e-14, "S should be -H: S={s}, H={h}");
    }

    // -----------------------------------------------------------------------
    // Non-equilibrium stress: zero at equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_neq_stress_zero_at_eq() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let (sxx, sxy, syy) = non_equilibrium_stress(&f_eq, &f_eq);
        assert!(sxx.abs() < 1e-14, "σ_xx neq should be 0: {sxx}");
        assert!(sxy.abs() < 1e-14, "σ_xy neq should be 0: {sxy}");
        assert!(syy.abs() < 1e-14, "σ_yy neq should be 0: {syy}");
    }

    // -----------------------------------------------------------------------
    // Non-equilibrium stress: non-zero when perturbed
    // -----------------------------------------------------------------------
    #[test]
    fn test_neq_stress_nonzero_when_perturbed() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        // Asymmetric perturbation: only perturb direction 1 (cx=1)
        f[1] += 0.05;

        let (sxx, _sxy, _syy) = non_equilibrium_stress(&f, &f_eq);
        assert!(sxx.abs() > 1e-10, "s_xx should be non-zero when perturbed");
    }

    // -----------------------------------------------------------------------
    // ELBM Diagnostics: record and check
    // -----------------------------------------------------------------------
    #[test]
    fn test_elbm_diagnostics() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let solver = EntropicD2Q9::new(nx, ny, nu);

        let mut diag = ElbmDiagnostics::new();
        assert_eq!(diag.step_count(), 0);

        diag.record(&solver);
        assert_eq!(diag.step_count(), 1);
        assert!(diag.entropy_history[0].is_finite());
        assert!(diag.alpha_min_history[0] >= 0.0);
        assert!(diag.alpha_avg_history[0] >= 0.0);
    }

    // -----------------------------------------------------------------------
    // ELBM Diagnostics: H-theorem check
    // -----------------------------------------------------------------------
    #[test]
    fn test_elbm_h_theorem_satisfied() {
        let diag = ElbmDiagnostics {
            entropy_history: vec![10.0, 9.5, 9.0, 8.5],
            alpha_min_history: vec![],
            alpha_avg_history: vec![],
            max_velocity_history: vec![],
        };
        assert!(diag.is_h_theorem_satisfied());

        let diag_bad = ElbmDiagnostics {
            entropy_history: vec![10.0, 9.5, 11.0, 8.5],
            alpha_min_history: vec![],
            alpha_avg_history: vec![],
            max_velocity_history: vec![],
        };
        assert!(!diag_bad.is_h_theorem_satisfied());
    }

    // -----------------------------------------------------------------------
    // Total entropy: finite for equilibrium
    // -----------------------------------------------------------------------
    #[test]
    fn test_total_entropy_finite() {
        let solver = EntropicD2Q9::new(6, 6, 1.0 / 6.0);
        let h = solver.total_entropy();
        assert!(h.is_finite(), "Total entropy should be finite: {h}");
    }

    // -----------------------------------------------------------------------
    // Max velocity at rest = 0
    // -----------------------------------------------------------------------
    #[test]
    fn test_max_velocity_at_rest() {
        let solver = EntropicD2Q9::new(6, 6, 1.0 / 6.0);
        let max_u = solver.max_velocity();
        assert!(max_u < 1e-14, "Max velocity at rest should be ~0: {max_u}");
    }

    // -----------------------------------------------------------------------
    // Density range at uniform state
    // -----------------------------------------------------------------------
    #[test]
    fn test_density_range_uniform() {
        let solver = EntropicD2Q9::new(6, 6, 1.0 / 6.0);
        let (rho_min, rho_max) = solver.density_range();
        assert!((rho_min - 1.0).abs() < 1e-14);
        assert!((rho_max - 1.0).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // KBC: decompose_general round-trip
    // -----------------------------------------------------------------------
    #[test]
    fn test_kbc_decompose_general_round_trip() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.03, -0.02]);
        let (s, h) = KbcD2Q9::decompose_general(&f_eq);
        for i in 0..9 {
            let reconstructed = s[i] + h[i];
            assert!(
                (reconstructed - f_eq[i]).abs() < 1e-14,
                "s+h != f at i={i}: {} vs {}",
                reconstructed,
                f_eq[i]
            );
        }
    }

    // -----------------------------------------------------------------------
    // KBC: total density conservation
    // -----------------------------------------------------------------------
    #[test]
    fn test_kbc_density_conservation() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let mut kbc = KbcD2Q9::new(nx, ny, nu);

        kbc.f[5][1] += 0.02;
        kbc.f[5][3] -= 0.01;

        let rho_before = kbc.total_density();
        kbc.step();
        let rho_after = kbc.total_density();

        assert!(
            (rho_before - rho_after).abs() < 1e-10,
            "KBC density not conserved: before={rho_before}, after={rho_after}"
        );
    }

    // -----------------------------------------------------------------------
    // KBC: uniform flow stays uniform
    // -----------------------------------------------------------------------
    #[test]
    fn test_kbc_uniform_flow() {
        let nx = 6;
        let ny = 6;
        let nu = 1.0 / 6.0;
        let rho0 = 1.0;
        let u0 = [0.03, 0.0];

        let mut kbc = KbcD2Q9::new(nx, ny, nu);
        kbc.set_uniform(rho0, u0);

        for _ in 0..10 {
            kbc.step();
        }

        let n = nx * ny;
        for k in 0..n {
            let (rho, u) = kbc.macros(k);
            assert!(
                (rho - rho0).abs() < 1e-12,
                "KBC density changed at cell {k}"
            );
            assert!((u[0] - u0[0]).abs() < 1e-12, "KBC ux changed at cell {k}");
        }
    }

    // -----------------------------------------------------------------------
    // d2q9_equilibrium free function agrees with method
    // -----------------------------------------------------------------------
    #[test]
    fn test_d2q9_equilibrium_free_fn() {
        let rho = 1.5;
        let u = [0.04, -0.03];
        let feq1 = d2q9_equilibrium(rho, u);
        let feq2 = EntropicD2Q9::equilibrium(rho, u);
        for i in 0..9 {
            assert!(
                (feq1[i] - feq2[i]).abs() < 1e-14,
                "Free fn vs method mismatch at i={i}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // Mirror state: f_mirror + f = 2 * f_eq
    // -----------------------------------------------------------------------
    #[test]
    fn test_mirror_state_identity() {
        let f_eq = EntropicD2Q9::equilibrium(1.0, [0.05, 0.02]);
        let mut f = f_eq;
        f[1] += 0.01;
        f[3] -= 0.01;

        let fm = mirror_state(&f, &f_eq);
        for i in 0..9 {
            let sum = f[i] + fm[i];
            let expected = 2.0 * f_eq[i];
            assert!(
                (sum - expected).abs() < 1e-14,
                "f + f_mirror != 2*f_eq at i={i}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Entropy production and H-theorem diagnostics
// ---------------------------------------------------------------------------

/// Compute the symmetrised Kullback–Leibler divergence of f from f_eq:
///
/// D_sym(f || f_eq) = D_KL(f||f_eq) + D_KL(f_eq||f)
///
/// This is symmetric: D_sym(f,feq) = D_sym(feq,f).
/// Returns `f64::INFINITY` if any distribution has non-positive entries.
pub fn symmetric_kl_divergence(f: &[f64; 9], feq: &[f64; 9]) -> f64 {
    let mut d = 0.0;
    for i in 0..9 {
        if f[i] <= 0.0 || feq[i] <= 0.0 {
            return f64::INFINITY;
        }
        d += (f[i] - feq[i]) * (f[i] / feq[i]).ln();
    }
    d
}

/// Compute the entropy production per step: ΔH = H(f) - H(f*).
///
/// According to the H-theorem this must be ≥ 0.
pub fn entropy_production(f_pre: &[f64; 9], f_post: &[f64; 9]) -> f64 {
    h_function(f_pre) - h_function(f_post)
}

// ---------------------------------------------------------------------------
// Entropic stabiliser — Newton–Raphson alpha solver (improved)
// ---------------------------------------------------------------------------

/// Solve the entropy equation `H(f + α (feq − f)) = H(f)` for the entropic
/// parameter α with a bracketed Newton–Raphson iteration.
///
/// This is the same equation as in [`find_entropic_alpha`], with an explicit
/// stopping tolerance on the Newton step and a caller-supplied starting
/// point: the iteration starts at the α for which `f + α (feq − f)` is the
/// least-squares fit of `f_mirror` (the usual choice is the mirror state
/// `2 feq − f`, i.e. α = 2).  The iteration is safeguarded inside the bracket
/// of the root, is not clamped to `(0, 2]`, and the result obeys the discrete
/// H-theorem (`H(f + α (feq − f)) ≤ H(f)` to 1e-12).
///
/// `α = 2` is returned when the state is within the near-equilibrium
/// tolerance of `feq` or when no entropic root exists (see
/// [`find_entropic_alpha`]).
pub fn solve_alpha_newton_raphson(
    f: &[f64; 9],
    feq: &[f64; 9],
    f_mirror: &[f64; 9],
    tol: f64,
    max_iter: usize,
) -> f64 {
    if is_near_equilibrium(f, feq) {
        return 2.0;
    }
    let d = mirror_direction(f, feq);
    let Some(bracket) = bracket_entropic_root(f, &d) else {
        return 2.0;
    };
    let norm2: f64 = d.iter().map(|x| x * x).sum();
    let alpha0 = if norm2 > 0.0 {
        let proj: f64 = (0..9).map(|i| (f_mirror[i] - f[i]) * d[i]).sum();
        proj / norm2
    } else {
        2.0
    };
    safeguarded_newton(f, &d, bracket, alpha0, tol, max_iter)
}

// ---------------------------------------------------------------------------
// Entropic over-relaxation (EOR) scheme
// ---------------------------------------------------------------------------

/// Entropic over-relaxation (EOR) for D2Q9.
///
/// Uses α = 2 (full over-relaxation toward the mirror state) when the
/// H-theorem is satisfied, and otherwise the largest α in \[1, 2\] that still
/// satisfies it.  The collision is `f' = f + α β (f_eq − f)` with the physical
/// relaxation parameter `β = 1/(6ν + 1)`, so α = 2 is exactly BGK with
/// `τ = 3ν + 1/2`.
///
/// This provides second-order accuracy in time while maintaining the
/// discrete H-theorem unconditionally.
#[derive(Debug, Clone)]
pub struct EntropyOverRelaxation {
    /// Base kinematic viscosity (lattice units); must be positive.
    pub nu: f64,
    /// Whether to enforce H-theorem (true) or always use α=2.
    pub enforce_h: bool,
}

impl EntropyOverRelaxation {
    /// Construct EOR from kinematic viscosity.
    ///
    /// `nu` must be positive (`τ = 3ν + 1/2 > 1/2`); use
    /// [`try_new`](Self::try_new) to have an invalid value reported.
    pub fn new(nu: f64) -> Self {
        Self {
            nu,
            enforce_h: true,
        }
    }

    /// Validated constructor: rejects a viscosity that is not finite and
    /// positive (checked with [`Error::check_tau`]).
    pub fn try_new(nu: f64) -> Result<Self> {
        validate_viscosity(nu)?;
        Ok(Self::new(nu))
    }

    /// Relaxation parameter `β = 1/(2τ) = 1/(6ν + 1)`.
    pub fn beta(&self) -> f64 {
        relaxation_beta(self.nu)
    }

    /// Compute the EOR relaxation parameter for a single cell.
    ///
    /// Returns α ∈ (0, 2] — typically close to 2 for smooth flows.
    pub fn compute_alpha(&self, f: &[f64; 9], rho: f64, u: [f64; 2]) -> f64 {
        let feq = d2q9_equilibrium(rho, u);
        if !self.enforce_h {
            return 2.0;
        }
        let fm = mirror_state(f, &feq);
        // Quick check: does α=2 satisfy H-theorem?
        let h_pre = h_function(f);
        let h_post = h_function(&fm);
        if h_post.is_finite() && h_post <= h_pre + 1e-12 {
            return 2.0;
        }
        // Fallback: binary search for valid α
        let mut lo = 1.0_f64;
        let mut hi = 2.0_f64;
        for _ in 0..32 {
            let mid = 0.5 * (lo + hi);
            let mut ftry = [0.0f64; 9];
            for i in 0..9 {
                ftry[i] = f[i] + mid * (feq[i] - f[i]);
            }
            if h_function(&ftry) <= h_pre + 1e-12 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        lo
    }

    /// Apply EOR collision to a single cell:
    /// `f' = f + α β (f_eq − f)` with `β = 1/(6ν + 1)`.
    pub fn collide(&self, f: &[f64; 9], rho: f64, u: [f64; 2]) -> [f64; 9] {
        let feq = d2q9_equilibrium(rho, u);
        let alpha = self.compute_alpha(f, rho, u);
        let beta = self.beta();
        let mut f_out = [0.0f64; 9];
        for i in 0..9 {
            f_out[i] = f[i] + alpha * beta * (feq[i] - f[i]);
        }
        f_out
    }
}

// ---------------------------------------------------------------------------
// Lyapunov stability indicator
// ---------------------------------------------------------------------------

/// Compute the normalised non-equilibrium entropy indicator.
///
/// δ = H(f) - H(feq) >= 0 (H-theorem).
///
/// This provides a local measure of how far the system is from equilibrium
/// and can be used as a refinement indicator in adaptive LBM.
pub fn neq_entropy_indicator(f: &[f64; 9], rho: f64, u: [f64; 2]) -> f64 {
    let feq = d2q9_equilibrium(rho, u);
    let hf = h_function(f);
    let hfeq = h_function(&feq);
    if hf.is_finite() && hfeq.is_finite() {
        (hf - hfeq).max(0.0)
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Additional entropic tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod extra_entropic_tests {
    use super::*;

    /// Symmetric KL divergence is zero when f == feq.
    #[test]
    fn test_symmetric_kl_at_eq() {
        let feq = d2q9_equilibrium(1.0, [0.05, -0.02]);
        let d = symmetric_kl_divergence(&feq, &feq);
        assert!(d.abs() < 1e-12, "Symmetric KL at eq = {d}");
    }

    /// Symmetric KL divergence ≥ 0 by convexity.
    #[test]
    fn test_symmetric_kl_nonneg() {
        let feq = d2q9_equilibrium(1.0, [0.03, 0.01]);
        let mut f = feq;
        f[1] += 0.005;
        f[3] -= 0.005;
        let d = symmetric_kl_divergence(&f, &feq);
        assert!(d >= 0.0, "Symmetric KL should be non-negative: {d}");
    }

    /// H-theorem: entropy production ≥ 0 after BGK step (α=1).
    #[test]
    fn test_entropy_production_nonneg() {
        let rho = 1.0;
        let u = [0.04, -0.02];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[1] += 0.01;
        f[3] -= 0.01;
        let feq2 = d2q9_equilibrium(rho, u);
        // BGK post-collision: f* = f + (feq - f) [alpha=1, BGK]
        let mut f_post = [0.0f64; 9];
        for i in 0..9 {
            f_post[i] = f[i] + (feq2[i] - f[i]);
        }
        let dp = entropy_production(&f, &f_post);
        assert!(dp >= -1e-13, "Entropy production negative: {dp}");
    }

    /// EOR: mass is conserved after collision.
    #[test]
    fn test_eor_mass_conservation() {
        let eor = EntropyOverRelaxation::new(1.0 / 6.0);
        let rho = 1.0_f64;
        let u = [0.03, 0.01];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[0] += 0.005;
        f[1] -= 0.003;
        f[3] -= 0.002;
        let rho_in: f64 = f.iter().sum();
        let f_out = eor.collide(&f, rho_in, u);
        let rho_out: f64 = f_out.iter().sum();
        assert!(
            (rho_out - rho_in).abs() < 1e-12,
            "EOR mass: {rho_in} -> {rho_out}"
        );
    }

    /// EOR: alpha is in (0, 2].
    #[test]
    fn test_eor_alpha_range() {
        let eor = EntropyOverRelaxation::new(0.05);
        let rho = 1.0_f64;
        let u = [0.05, -0.02];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[1] += 0.005;
        f[3] -= 0.005;
        let alpha = eor.compute_alpha(&f, rho, u);
        assert!(
            alpha > 0.0 && alpha <= 2.0 + 1e-10,
            "alpha={alpha} out of range"
        );
    }

    /// solve_alpha_newton_raphson: at equilibrium alpha → 2.
    #[test]
    fn test_alpha_nr_at_equilibrium() {
        let feq = d2q9_equilibrium(1.0, [0.04, -0.01]);
        let fm = mirror_state(&feq, &feq);
        let alpha = solve_alpha_newton_raphson(&feq, &feq, &fm, 1e-12, 50);
        assert!(alpha > 0.0 && alpha <= 2.0 + 1e-10, "alpha={alpha}");
    }

    /// neq_entropy_indicator: at equilibrium indicator ≈ 0.
    #[test]
    fn test_neq_entropy_indicator_at_eq() {
        let rho = 1.0_f64;
        let u = [0.03, 0.02];
        let feq = d2q9_equilibrium(rho, u);
        let ind = neq_entropy_indicator(&feq, rho, u);
        assert!(ind < 1e-12, "indicator at eq = {ind}");
    }

    /// neq_entropy_indicator: perturbed state has larger indicator.
    #[test]
    fn test_neq_entropy_indicator_perturbed() {
        let rho = 1.0_f64;
        let u = [0.03, 0.02];
        let feq = d2q9_equilibrium(rho, u);
        let mut f = feq;
        f[1] += 0.01;
        f[3] -= 0.01;
        let ind = neq_entropy_indicator(&f, rho, u);
        assert!(ind >= 0.0, "indicator = {ind}");
    }

    /// entropy_production with α=2 (mirror) should be exactly zero.
    #[test]
    fn test_entropy_production_mirror() {
        let feq = d2q9_equilibrium(1.0, [0.04, 0.01]);
        let mut f = feq;
        f[1] += 0.005;
        f[3] -= 0.005;
        let fm = mirror_state(&f, &feq);
        // After mirror step H(fm) should equal H(f_mirror)
        let dp = entropy_production(&f, &fm);
        // dp = H(f) - H(fm); mirror should satisfy H-theorem
        assert!(dp >= -1e-12, "entropy production negative: {dp}");
    }

    /// h_function is lower for equilibrium than for perturbed state.
    #[test]
    fn test_h_function_minimum_at_eq() {
        let feq = d2q9_equilibrium(1.0, [0.02, -0.01]);
        let mut f = feq;
        f[2] += 0.01;
        f[4] -= 0.01;
        let h_eq = h_function(&feq);
        let h_f = h_function(&f);
        assert!(h_f >= h_eq - 1e-12, "H(f)={h_f} < H(feq)={h_eq}");
    }
}

#[cfg(test)]
mod solver_tests;

#[cfg(test)]
mod kbc_tests;
