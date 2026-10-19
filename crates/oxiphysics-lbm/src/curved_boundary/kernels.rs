// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Per-link curved-wall kernels.
//!
//! Notation (shared by every function in this file): `i` is the lattice
//! direction pointing from the fluid node `x_f` into the wall, `ī` its
//! opposite, `x_ff = x_f − c_i` and `x_fff = x_f − 2 c_i` the upstream fluid
//! nodes on the same link line, `f̃` post-collision (pre-streaming)
//! populations, and `q ∈ (0, 1]` the fraction of the link `x_f → x_f + c_i`
//! that lies in the fluid.  Every kernel returns the unknown post-streaming
//! population `f_ī(x_f, t + 1)`.

use super::lattice::{CurvedLattice, dot3};
use crate::error::{Error, Result};
use crate::lattice::CS2;

/// Post-collision populations sampled along one boundary link line.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct LinkPopulations {
    /// `f̃_i(x_f)`: population at the fluid node heading into the wall.
    pub f_i_xf: f64,
    /// `f̃_i(x_f − c_i)`.
    pub f_i_xff: f64,
    /// `f̃_i(x_f − 2 c_i)`.
    pub f_i_xfff: f64,
    /// `f̃_ī(x_f)`: population at the fluid node heading away from the wall.
    pub f_opp_xf: f64,
    /// `f̃_ī(x_f − c_i)`.
    pub f_opp_xff: f64,
}

/// Geometry and kinematics of one wall-cut lattice link.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WallLink {
    /// Fraction of the link lying in the fluid, `q ∈ (0, 1]`.
    pub q: f64,
    /// Lattice weight `w_i` of the link direction.
    pub w_i: f64,
    /// Link velocity `c_i` (pointing from the fluid node into the wall).
    pub c_i: [f64; 3],
    /// Wall velocity at the intersection point `x_f + q c_i`.
    pub u_wall: [f64; 3],
    /// Fluid density used in the moving-wall momentum term.
    pub rho: f64,
}

impl WallLink {
    /// Build a link description, validating `q ∈ (0, 1]`, `w_i > 0` and
    /// finite inputs.
    pub fn new(q: f64, w_i: f64, c_i: [f64; 3], u_wall: [f64; 3], rho: f64) -> Result<Self> {
        if !q.is_finite() || q <= 0.0 || q > 1.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "q".to_string(),
                value: q,
                min: 0.0,
                max: 1.0,
            });
        }
        if !w_i.is_finite() || w_i <= 0.0 {
            return Err(Error::ParameterOutOfBounds {
                name: "w_i".to_string(),
                value: w_i,
                min: 0.0,
                max: 1.0,
            });
        }
        let all_finite = c_i
            .iter()
            .chain(u_wall.iter())
            .chain(std::iter::once(&rho))
            .all(|v| v.is_finite());
        if !all_finite {
            return Err(Error::General(
                "wall link velocities and density must be finite".to_string(),
            ));
        }
        Ok(Self {
            q,
            w_i,
            c_i,
            u_wall,
            rho,
        })
    }

    /// Ladd moving-wall momentum term `2 w_i ρ (c_i · u_w) / c_s²`.
    pub fn moving_wall_term(&self) -> f64 {
        bouzidi_wall_correction(self.rho, self.w_i, self.c_i, self.u_wall)
    }
}

// ---------------------------------------------------------------------------
// Moving-wall term
// ---------------------------------------------------------------------------

/// Moving-wall momentum term `2 w_i ρ (c_i · u_w) / c_s²` (Ladd 1994).
///
/// It is *subtracted* from the bounced population when `c_i` points from the
/// fluid into the wall: `f_ī(x_f, t+1) = f̃_i(x_f, t) − 2 w_i ρ (c_i·u_w)/c_s²`.
pub fn bouzidi_wall_correction(rho: f64, w_i: f64, c_i: [f64; 3], u_wall: [f64; 3]) -> f64 {
    2.0 * w_i * rho * dot3(c_i, u_wall) / CS2
}

// ---------------------------------------------------------------------------
// Bouzidi schemes
// ---------------------------------------------------------------------------

/// Linear Bouzidi–Firdaouss–Lallemand interpolated bounce-back.
///
/// ```text
/// q <  1/2:  f_ī(x_f) = 2q f̃_i(x_f) + (1 − 2q) f̃_i(x_ff) − M
/// q ≥ 1/2:  f_ī(x_f) = f̃_i(x_f)/(2q) + (2q − 1)/(2q) f̃_ī(x_f) − M/(2q)
/// ```
///
/// with `M = 2 w_i ρ (c_i·u_w)/c_s²`.  For `q ≥ 1/2` the bounced particle
/// lands at `x_f + (2q − 1) c_i`, whose interpolation weight is `1/(2q)`, so
/// the moving-wall term carries the same factor.  At `q = 1/2` both branches
/// reduce to the half-way bounce-back with the Ladd correction.
pub fn bouzidi_bc(pops: &LinkPopulations, link: &WallLink) -> f64 {
    let q = link.q;
    let m = link.moving_wall_term();
    if q < 0.5 {
        2.0 * q * pops.f_i_xf + (1.0 - 2.0 * q) * pops.f_i_xff - m
    } else {
        let inv_2q = 1.0 / (2.0 * q);
        inv_2q * pops.f_i_xf + (2.0 * q - 1.0) * inv_2q * pops.f_opp_xf - inv_2q * m
    }
}

/// Quadratic Bouzidi–Firdaouss–Lallemand interpolated bounce-back.
///
/// ```text
/// q <  1/2:  f_ī = q(1+2q) f̃_i(x_f) + (1−4q²) f̃_i(x_ff) − q(1−2q) f̃_i(x_fff) − M
/// q ≥ 1/2:  f_ī = [f̃_i(x_f) − M] / (q(2q+1)) + (2q−1)/q f̃_ī(x_f)
///                 − (2q−1)/(2q+1) f̃_ī(x_ff)
/// ```
///
/// The coefficients are the Lagrange weights of the three interpolation
/// nodes and sum to one in both branches.
pub fn bouzidi_quadratic_bc(pops: &LinkPopulations, link: &WallLink) -> f64 {
    let q = link.q;
    let m = link.moving_wall_term();
    if q < 0.5 {
        q * (1.0 + 2.0 * q) * pops.f_i_xf + (1.0 - 4.0 * q * q) * pops.f_i_xff
            - q * (1.0 - 2.0 * q) * pops.f_i_xfff
            - m
    } else {
        let l_wall = 1.0 / (q * (2.0 * q + 1.0));
        l_wall * (pops.f_i_xf - m) + (2.0 * q - 1.0) / q * pops.f_opp_xf
            - (2.0 * q - 1.0) / (2.0 * q + 1.0) * pops.f_opp_xff
    }
}

/// Lagrange quadratic through `(0, f_0)`, `(1, f_1)`, `(2, f_2)` evaluated at `t`.
///
/// Generic helper for interpolating a quantity sampled at three equally spaced
/// nodes along a link line.
pub fn quadratic_interpolation_bc(f_0: f64, f_1: f64, f_2: f64, t: f64) -> f64 {
    let l0 = (t - 1.0) * (t - 2.0) / 2.0;
    let l1 = -t * (t - 2.0);
    let l2 = t * (t - 1.0) / 2.0;
    f_0 * l0 + f_1 * l1 + f_2 * l2
}

// ---------------------------------------------------------------------------
// Unified interpolated bounce-back (Yu, Mei, Luo & Shyy 2003)
// ---------------------------------------------------------------------------

/// Unified linear interpolated bounce-back of Yu, Mei, Luo & Shyy (2003).
///
/// The streamed population `f_i` is interpolated to the wall point
/// (`f_i(x_w) = (1−q) f̃_i(x_ff) + q f̃_i(x_f)`), bounced there with the Ladd
/// term, and the result is interpolated back to `x_f` together with
/// `f_ī(x_ff, t+1) = f̃_ī(x_f)`:
///
/// ```text
/// f_ī(x_f) = [(1 − q) f̃_i(x_ff) + q f̃_i(x_f) + q f̃_ī(x_f) − M] / (1 + q)
/// ```
///
/// A single formula covers the whole range `q ∈ (0, 1]`.
pub fn interpolated_bounce_back(pops: &LinkPopulations, link: &WallLink) -> f64 {
    let q = link.q;
    ((1.0 - q) * pops.f_i_xff + q * pops.f_i_xf + q * pops.f_opp_xf - link.moving_wall_term())
        / (1.0 + q)
}

/// Standard half-way bounce-back.
///
/// Returns `f_post[outgoing]`, i.e. `f_ī(x_f, t+1) = f̃_i(x_f, t)` for a
/// stationary wall located half-way along the link; `0` if `outgoing` is out
/// of range.
pub fn half_way_bounce_back(f_post: &[f64], outgoing: usize) -> f64 {
    f_post.get(outgoing).copied().unwrap_or(0.0)
}

/// Half-way bounce-back on a moving wall (Ladd correction).
///
/// `f_post` holds the node's post-collision populations and `f` its
/// post-streaming populations.  The unknown `f[ī]` (`ī = opposite(dir)`) is set
/// to `f_post[dir] − 2 w_dir ρ (c_dir · u_w)/c_s²` and returned.
///
/// # Errors
/// Returns an error if `dir` is not a valid direction of `lattice` or the
/// population slices do not have `Q` entries.
pub fn moving_wall_bc(
    f_post: &[f64],
    f: &mut [f64],
    lattice: &CurvedLattice,
    dir: usize,
    u_wall: [f64; 3],
    rho: f64,
) -> Result<f64> {
    let nq = lattice.num_directions();
    if f_post.len() != nq || f.len() != nq {
        return Err(Error::ComponentCountMismatch {
            expected: nq,
            found: f_post.len().min(f.len()),
            field: "moving_wall_bc populations".to_string(),
        });
    }
    let (c_i, w_i, opp) = match (
        lattice.velocity_f64(dir),
        lattice.weight(dir),
        lattice.opposite(dir),
    ) {
        (Some(c), Some(w), Some(o)) => (c, w, o),
        _ => {
            return Err(Error::ParameterOutOfBounds {
                name: "dir".to_string(),
                value: dir as f64,
                min: 0.0,
                max: (nq - 1) as f64,
            });
        }
    };
    let outgoing = f_post.get(dir).copied().unwrap_or(0.0);
    let value = outgoing - bouzidi_wall_correction(rho, w_i, c_i, u_wall);
    if let Some(slot) = f.get_mut(opp) {
        *slot = value;
    }
    Ok(value)
}

// ---------------------------------------------------------------------------
// Filippova–Hänel / Mei–Luo–Shyy
// ---------------------------------------------------------------------------

/// Weighting factor `χ` of the Filippova–Hänel scheme.
///
/// With `ω = 1/τ` (Mei, Shyy, Yu & Luo, *J. Comput. Phys.* 161, 2000):
///
/// ```text
/// q ≥ 1/2:                χ = (2q − 1)/(τ + 1/2) = ω (2q − 1) / (1 + ω/2)
/// q <  1/2, MLS (u_ff):    χ = (2q − 1)/(τ − 2)   = ω (2q − 1) / (1 − 2ω)
/// q <  1/2, original FH:   χ = (2q − 1)/(τ − 1)   = ω (2q − 1) / (1 − ω)
/// ```
///
/// The `q ≥ 1/2` weight is the one consistent with the extrapolated
/// `u_bf = (1 − 3/(2q)) u_f + 3/(2q) u_w` used by [`filippova_hanel`].
///
/// # Errors
/// Returns an error if `ω ∉ (0, 2)` or the `q < 1/2` denominator vanishes
/// (`τ = 2` for MLS, `τ = 1` for the original scheme).
pub fn mei_luo_shyy_chi(q: f64, omega: f64, has_upstream: bool) -> Result<f64> {
    if !omega.is_finite() || omega <= 0.0 || omega >= 2.0 {
        return Err(Error::ParameterOutOfBounds {
            name: "omega".to_string(),
            value: omega,
            min: 0.0,
            max: 2.0,
        });
    }
    if q >= 0.5 {
        return Ok(omega * (2.0 * q - 1.0) / (1.0 + 0.5 * omega));
    }
    let denom = if has_upstream {
        1.0 - 2.0 * omega
    } else {
        1.0 - omega
    };
    if denom.abs() < 1e-12 {
        return Err(Error::IncompatibleParameters {
            message: format!(
                "Filippova-Hanel weighting is singular for q={q} at omega={omega} \
                 (tau={})",
                1.0 / omega
            ),
        });
    }
    Ok(omega * (2.0 * q - 1.0) / denom)
}

/// Fictitious equilibrium `f*_i(x_w)` of the Filippova–Hänel scheme:
///
/// ```text
/// f*_i = w_i ρ [1 + c_i·u_bf/c_s² + (c_i·u_f)²/(2c_s⁴) − u_f·u_f/(2c_s²)]
/// ```
pub fn fictitious_equilibrium(link: &WallLink, u_f: [f64; 3], u_bf: [f64; 3]) -> f64 {
    let cu_bf = dot3(link.c_i, u_bf);
    let cu_f = dot3(link.c_i, u_f);
    let uu = dot3(u_f, u_f);
    link.w_i * link.rho * (1.0 + cu_bf / CS2 + cu_f * cu_f / (2.0 * CS2 * CS2) - uu / (2.0 * CS2))
}

/// Filippova–Hänel curved bounce-back with the Mei–Luo–Shyy improvement.
///
/// ```text
/// f_ī(x_f, t+1) = (1 − χ) f̃_i(x_f) + χ f*_i(x_w) − 2 w_i ρ (c_i·u_w)/c_s²
/// ```
///
/// `u_bf = (1 − 3/(2q)) u_f + 3/(2q) u_w` for `q ≥ 1/2`; for `q < 1/2` it is
/// `u_ff = u(x_f − c_i)` (Mei–Luo–Shyy) when `u_ff` is available, otherwise
/// `u_f` (original Filippova–Hänel).  `χ` is given by [`mei_luo_shyy_chi`];
/// at `q = 1/2` the scheme reduces to half-way bounce-back.  Linear (Couette)
/// profiles are reproduced exactly for every `q` and `τ`.
///
/// # Arguments
/// * `f_i_xf` — post-collision `f̃_i(x_f)`
/// * `link`   — link geometry; `link.rho` is the fluid-node density `ρ_f`
/// * `u_f`    — fluid velocity at `x_f`
/// * `u_ff`   — fluid velocity at `x_f − c_i`, if that node is fluid
/// * `omega`  — BGK relaxation frequency `1/τ`
///
/// # Errors
/// Propagates the errors of [`mei_luo_shyy_chi`].
pub fn filippova_hanel(
    f_i_xf: f64,
    link: &WallLink,
    u_f: [f64; 3],
    u_ff: Option<[f64; 3]>,
    omega: f64,
) -> Result<f64> {
    let q = link.q;
    let chi = mei_luo_shyy_chi(q, omega, u_ff.is_some())?;
    let u_bf = if q >= 0.5 {
        let a = 3.0 / (2.0 * q);
        [
            (1.0 - a) * u_f[0] + a * link.u_wall[0],
            (1.0 - a) * u_f[1] + a * link.u_wall[1],
            (1.0 - a) * u_f[2] + a * link.u_wall[2],
        ]
    } else {
        u_ff.unwrap_or(u_f)
    };
    let f_star = fictitious_equilibrium(link, u_f, u_bf);
    Ok((1.0 - chi) * f_i_xf + chi * f_star - link.moving_wall_term())
}

// ---------------------------------------------------------------------------
// Immersed boundary direct forcing
// ---------------------------------------------------------------------------

/// Immersed-boundary direct-forcing velocity relaxation.
///
/// ```text
/// u_new = u + ε (u_target − u)
/// ```
///
/// `ε = 1` imposes the target velocity exactly; `ε ∈ (0, 1)` relaxes towards it.
pub fn direct_forcing(u: &mut [f64; 3], u_target: [f64; 3], epsilon: f64) {
    for d in 0..3 {
        u[d] += epsilon * (u_target[d] - u[d]);
    }
}

// ---------------------------------------------------------------------------
// D2Q9 equilibrium and signed-distance helpers
// ---------------------------------------------------------------------------

/// D2Q9 equilibrium for density `rho` and velocity `(ux, uy)`.
///
/// Direction order: rest, E, N, W, S, NE, NW, SW, SE (as [`CurvedLattice::d2q9`]).
pub fn d2q9_equilibrium(rho: f64, ux: f64, uy: f64) -> [f64; 9] {
    let lattice = CurvedLattice::d2q9();
    let mut feq = [0.0f64; 9];
    for (i, v) in feq.iter_mut().enumerate() {
        *v = lattice.equilibrium(i, rho, [ux, uy, 0.0]);
    }
    feq
}

/// Signed distance to a sphere (circle in 2-D): positive outside, negative inside.
pub fn sphere_sdf(pos: [f64; 3], center: [f64; 3], radius: f64) -> f64 {
    let d = (0..3)
        .map(|i| (pos[i] - center[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    d - radius
}

/// Signed distance to an axis-aligned box with half-extents `half`:
/// positive outside, negative inside.
pub fn box_sdf(pos: [f64; 3], center: [f64; 3], half: [f64; 3]) -> f64 {
    let dx = (pos[0] - center[0]).abs() - half[0];
    let dy = (pos[1] - center[1]).abs() - half[1];
    let dz = (pos[2] - center[2]).abs() - half[2];
    let outside = [dx.max(0.0), dy.max(0.0), dz.max(0.0)];
    let outside_dist = (outside[0].powi(2) + outside[1].powi(2) + outside[2].powi(2)).sqrt();
    let inside_dist = dx.max(dy).max(dz).min(0.0);
    outside_dist + inside_dist
}
