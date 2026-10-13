// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ferrari-Canny grasp quality: the radius of the largest origin-centred ball
//! contained in the convex hull of the primitive wrenches.
//!
//! For unit directions `u` let `h(u) = max_i u . w_i` be the support function
//! of `conv(W)`.  The ball of radius `rho` lies in the hull iff `h(u) >= rho`
//! for every unit `u`, so the metric is
//!
//! ```text
//! Q1 = min_{|u| = 1} h(u)          (Q1 > 0  <=>  the origin is interior,
//!                                    i.e. the grasp has force closure)
//! ```
//!
//! This minimisation of a convex function over a sphere is not convex, but it
//! has a closed combinatorial solution.  `h(y) <= 1` describes the *polar*
//! polytope `P = { y : w_i . y <= 1 }`, and the boundary point of `P` in
//! direction `u` is `u / h(u)`, hence
//!
//! ```text
//! Q1 = 1 / max_{y in P} |y| = 1 / max over the vertices y of P of |y| ,
//! ```
//!
//! because a convex function attains its maximum over a polytope at a vertex.
//! The vertices of `P` (equivalently the facets of `conv(W)`) are the extreme
//! rays with `s > 0` (scaled to `s = 1`) of the homogenised cone
//! `{ (y, s) : w_i . y - s <= 0, s >= 0 }`, which the double-description
//! method of [`ConeEnumerator`] maintains while constraints are added.  An
//! extreme ray with `s = 0` that no constraint cuts off means `P` is
//! unbounded, i.e. the origin is not an interior point and the quality is zero.
//!
//! Only the farthest vertex is needed, so the constraints are added lazily by a
//! cutting-plane loop: the vertex of the current outer approximation that is
//! farthest from the origin is tested against all constraints, and the most
//! violated one is intersected next.  When the farthest vertex satisfies every
//! constraint it is a vertex of `P` and, because the approximation contains
//! `P`, the maximiser of `|y|`: the result is exact.  Typically only a small
//! fraction of the constraints is ever enumerated, which keeps the vertex sets
//! small even for finely discretised friction cones.  If the work budget of the
//! enumeration runs out first, the farthest vertex of the current approximation
//! still proves a lower bound of Q1, which is returned and flagged as inexact.
//!
//! The wrench components mix forces and torques, so the value depends on the
//! units (and on the torque reference point); it is the usual convention-based
//! metric of the grasping literature.

use super::polar_cone::{Budget, ConeEnumerator, ZERO_TOL};
use super::wrench_set::{distinct_nonzero, distinct_unit};
use super::{Wrench, force_closure_test};

/// Dimension of the homogenised polar cone: six wrench coordinates plus `s`.
const CONE_DIM: usize = 7;

/// A polar vertex whose scale `s` is below this lies (numerically) at
/// infinity: the origin is not an interior point of the hull.
const BOUNDARY_TOL: f64 = 1e-9;

/// A wrench set needs at least seven members for an interior origin in `R^6`.
const MIN_WRENCHES: usize = 7;

/// Result of an inscribed-ball computation.
#[derive(Debug, Clone, Copy)]
pub(super) struct InscribedBall {
    /// Radius of the ball, `0.0` if the origin is not an interior point.
    pub(super) radius: f64,
    /// Unit direction `u*` in which the support function attains `radius`
    /// (the direction of the farthest polar vertex); zero if `radius = 0`.  For
    /// an inexact result it is the direction of the farthest vertex of the
    /// approximation that proved the bound.
    pub(super) direction: [f64; 6],
    /// `false` if the enumeration budget ran out; `radius` is then a lower
    /// bound of the true value.
    pub(super) exact: bool,
}

impl InscribedBall {
    const EMPTY: Self = Self {
        radius: 0.0,
        direction: [0.0; 6],
        exact: true,
    };
}

/// Inscribed ball of the convex hull of the wrenches themselves (Q1).
pub(super) fn q1_ball(wrenches: &[Wrench], budget: Budget) -> InscribedBall {
    match distinct_nonzero(wrenches) {
        Some(set) => ball_of_clean_set(&set, budget),
        None => InscribedBall::EMPTY,
    }
}

/// Inscribed ball of the convex hull of the unit-normalised wrenches (the
/// epsilon metric).
pub(super) fn epsilon_ball(wrenches: &[Wrench], budget: Budget) -> InscribedBall {
    match distinct_unit(wrenches) {
        Some(set) => ball_of_clean_set(&set, budget),
        None => InscribedBall::EMPTY,
    }
}

/// The farthest known vertex of the polar polytope, as a Q1 bound.
#[derive(Debug, Clone, Copy)]
struct Candidate {
    radius: f64,
    direction: [f64; 6],
}

impl Candidate {
    /// An unbounded polar polytope: no positive ball fits.
    const UNBOUNDED: Self = Self {
        radius: 0.0,
        direction: [0.0; 6],
    };

    /// The candidate described by a (finite, `s > BOUNDARY_TOL`) polar ray.
    fn from_ray(ray: &[f64]) -> Option<Self> {
        let scale = *ray.get(6)?;
        let y = ray.get(..6)?;
        let y_norm = y.iter().map(|v| v * v).sum::<f64>().sqrt();
        if scale <= BOUNDARY_TOL || y_norm <= f64::MIN_POSITIVE {
            return None;
        }
        let mut direction = [0.0; 6];
        for (d, v) in direction.iter_mut().zip(y.iter()) {
            *d = v / y_norm;
        }
        Some(Self {
            radius: scale / y_norm,
            direction,
        })
    }
}

/// What the cutting-plane loop does next.
enum Step {
    /// The polar polytope is unbounded, i.e. the origin is not an interior
    /// point: Q1 is zero.
    Unbounded,
    /// The farthest vertex of the outer approximation satisfies every
    /// remaining constraint: it is the farthest vertex of the polar polytope.
    Finished(Candidate),
    /// Intersect this constraint next; `bound` is the Q1 bound the current
    /// outer approximation already proves.
    Cut { row: usize, bound: Candidate },
}

/// The unprocessed constraint violated most by `ray`, with its violation.
fn most_violated(cone: &ConeEnumerator, ray: &[f64]) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64)> = None;
    for index in 0..cone.row_count() {
        if cone.is_processed(index) {
            continue;
        }
        let violation: f64 = cone
            .row(index)
            .iter()
            .zip(ray.iter())
            .map(|(a, x)| a * x)
            .sum();
        if violation > ZERO_TOL && best.is_none_or(|(_, v)| violation > v) {
            best = Some((index, violation));
        }
    }
    best
}

/// Decide the next step from the current outer approximation of the polar
/// polytope.
fn next_step(cone: &ConeEnumerator) -> Step {
    let mut any_at_infinity = false;
    let mut infinite_cut: Option<(usize, f64)> = None;
    let mut farthest: Option<(Candidate, &[f64])> = None;
    for ray in cone.rays() {
        let at_infinity = ray.get(6).is_none_or(|&scale| scale <= BOUNDARY_TOL);
        if at_infinity {
            any_at_infinity = true;
            if let Some((row, violation)) = most_violated(cone, ray)
                && infinite_cut.is_none_or(|(_, v)| violation > v)
            {
                infinite_cut = Some((row, violation));
            }
        } else if let Some(candidate) = Candidate::from_ray(ray)
            && farthest.is_none_or(|(best, _)| candidate.radius < best.radius)
        {
            farthest = Some((candidate, ray));
        }
    }
    if any_at_infinity {
        // A direction of the outer approximation that no constraint cuts off
        // is a direction of the polar polytope itself.
        return match infinite_cut {
            Some((row, _)) => Step::Cut {
                row,
                bound: Candidate::UNBOUNDED,
            },
            None => Step::Unbounded,
        };
    }
    match farthest {
        None => Step::Unbounded,
        Some((candidate, ray)) => match most_violated(cone, ray) {
            None => Step::Finished(candidate),
            Some((row, _)) => Step::Cut {
                row,
                bound: candidate,
            },
        },
    }
}

/// Q1 of a set of finite, distinct, non-zero wrenches (see the module
/// documentation for the cutting-plane algorithm).
fn ball_of_clean_set(set: &[Wrench], budget: Budget) -> InscribedBall {
    // The per-axis sign test is a necessary condition for an interior origin
    // and rules out most non-closing sets without any enumeration.
    if set.len() < MIN_WRENCHES || !force_closure_test(set) {
        return InscribedBall::EMPTY;
    }
    // Homogenised polar cone rows a . (y, s) <= 0: `s >= 0` first, then
    // `w_i . y - s <= 0` for every wrench.
    let mut rows: Vec<Vec<f64>> = Vec::with_capacity(set.len() + 1);
    rows.push(vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0]);
    for w in set {
        rows.push(vec![w[0], w[1], w[2], w[3], w[4], w[5], -1.0]);
    }
    // Rank below seven: the wrenches do not span `R^6`.
    let Some(mut cone) = ConeEnumerator::new(&rows, CONE_DIM, &[0], budget) else {
        return InscribedBall::EMPTY;
    };
    loop {
        match next_step(&cone) {
            Step::Unbounded => return InscribedBall::EMPTY,
            Step::Finished(found) => {
                return InscribedBall {
                    radius: found.radius,
                    direction: found.direction,
                    exact: true,
                };
            }
            Step::Cut { row, bound } => {
                if cone.insert(row).is_err() {
                    // Out of budget: the outer approximation proves a lower
                    // bound of Q1.
                    return InscribedBall {
                        radius: bound.radius,
                        direction: bound.direction,
                        exact: false,
                    };
                }
            }
        }
    }
}

/// Q1 of a wrench set with the standard enumeration budget.
pub(super) fn q1_radius(wrenches: &[Wrench]) -> f64 {
    q1_ball(wrenches, Budget::STANDARD).radius
}

/// Ferrari-Canny Q1 together with the wrench direction in which the grasp is
/// weakest.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraspQuality {
    /// Q1: the radius of the largest origin-centred ball inside the convex hull
    /// of the primitive wrenches; zero when the contacts do not achieve force
    /// closure.
    pub q1: f64,
    /// Unit direction in wrench space in which the convex hull of the primitive
    /// wrenches is closest to the origin: along it the largest wrench in the hull
    /// has magnitude exactly `q1`, along every other direction at least `q1`.
    /// The zero vector when `q1` is zero.  If `exact` is `false` it is the
    /// direction of the approximation that proved the lower bound.
    pub weakest_direction: [f64; 6],
    /// `true` if `q1` is exact; `false` if the facet enumeration was cut off by
    /// its work budget, in which case `q1` is a lower bound of the true value.
    pub exact: bool,
}

/// Q1 grasp quality with its weakest wrench direction and an exactness flag.
///
/// [`grasp_quality_q1`](super::grasp_quality_q1) returns only `.q1` of this
/// result; see there for the definition and the algorithm.
pub fn grasp_quality_q1_detailed(wrenches: &[Wrench]) -> GraspQuality {
    let ball = q1_ball(wrenches, Budget::STANDARD);
    GraspQuality {
        q1: ball.radius,
        weakest_direction: ball.direction,
        exact: ball.exact,
    }
}

/// Epsilon metric of a wrench set with the standard enumeration budget.
pub(super) fn epsilon_radius(wrenches: &[Wrench]) -> f64 {
    epsilon_ball(wrenches, Budget::STANDARD).radius
}
