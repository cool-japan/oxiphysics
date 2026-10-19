// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Signed dihedral angle of a cloth hinge, its gradient and the XPBD
//! projection of a dihedral bending constraint.
//!
//! # Hinge
//!
//! A hinge is four particles `[p0, p1, pk, pl]`.  `(p0, p1)` is the shared
//! edge, `pk` is the wing particle of the face `A = (p0, p1, pk)` and `pl` the
//! wing particle of the face `B = (p1, p0, pl)`.  With the edge vector
//! `e = p1 - p0` the (unnormalised) face normals are
//!
//! ```text
//! N_A = e x (pk - p0)        N_B = (pl - p0) x e
//! ```
//!
//! For a flat sheet whose wings lie on opposite sides of the edge the two
//! normals coincide, and the *signed dihedral angle*
//!
//! ```text
//! theta = atan2( (N_A x N_B) . e_hat , N_A . N_B )
//! ```
//!
//! is the rotation about the edge that carries `N_A` onto `N_B`.  It is `0` for
//! a flat hinge and lies in `(-pi, pi]`.  Folding the wings towards the side
//! the face normals point to makes `theta` negative, folding them away from
//! that side makes it positive.  Unlike `acos(n_A . n_B)` the `atan2` form is
//! smooth through the flat state, so its gradient never blows up there.
//!
//! # Gradient
//!
//! With `l = |e|` the gradient of `theta` (Bridson, Marino and Fedkiw, 2003)
//! is
//!
//! ```text
//! d theta / d pk = -l N_A / |N_A|^2
//! d theta / d pl = -l N_B / |N_B|^2
//! d theta / d p0 = -((pk - p1) . e_hat) N_A / |N_A|^2
//!                  -((pl - p1) . e_hat) N_B / |N_B|^2
//! d theta / d p1 = -(d theta / d p0 + d theta / d pk + d theta / d pl)
//! ```
//!
//! The wing terms rotate a face about the edge; the hinge terms follow from
//! translation and rotation invariance of `theta` and from the fact that
//! sliding a hinge end along the edge does not change `theta`.  Consequently
//! the four gradients sum to zero and carry no net torque, which is what makes
//! the XPBD projection conserve linear and angular momentum.

use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

use super::functions::{add, cross, dot, length, scale, sub};

/// Smallest hinge edge length (m) that still defines a rotation axis.
const MIN_EDGE_LENGTH: f64 = 1e-12;

/// Smallest squared length of an unnormalised face normal (m^4).
///
/// A face whose normal is shorter has (numerically) zero area: the wing
/// particle lies on the hinge line and the face has no orientation.
const MIN_NORMAL_SQ: f64 = 1e-24;

/// Smallest generalised inverse mass `sum_i w_i |grad_i|^2` that can still move
/// a hinge.
const MIN_GENERALISED_INVERSE_MASS: f64 = 1e-30;

/// Signed dihedral angle of a hinge together with its gradient.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct HingeState {
    /// Signed dihedral angle in `(-pi, pi]`, `0` for a flat hinge.
    pub(super) angle: f64,
    /// Gradient of [`HingeState::angle`] with respect to the positions of
    /// `[p0, p1, pk, pl]` (rad/m).
    pub(super) gradients: [[f64; 3]; 4],
}

/// Result of one XPBD projection of a dihedral constraint.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct DihedralProjection {
    /// Change of the constraint multiplier, to be added to the accumulated
    /// multiplier of the constraint.
    pub(super) delta_lambda: f64,
    /// Position corrections for `[p0, p1, pk, pl]`.
    pub(super) corrections: [[f64; 3]; 4],
}

/// Signed dihedral angle of the hinge `[p0, p1, pk, pl]` and its gradient.
///
/// Returns `None` when the hinge is degenerate: the edge `(p0, p1)` has no
/// length, a wing particle lies on the hinge line (so a face has no
/// orientation), or an input coordinate is not finite.
pub(super) fn hinge_state(
    p0: [f64; 3],
    p1: [f64; 3],
    pk: [f64; 3],
    pl: [f64; 3],
) -> Option<HingeState> {
    let edge = sub(p1, p0);
    let edge_len = length(edge);
    if !edge_len.is_finite() || edge_len <= MIN_EDGE_LENGTH {
        return None;
    }
    let edge_hat = scale(edge, 1.0 / edge_len);
    let normal_a = cross(edge, sub(pk, p0));
    let normal_b = cross(sub(pl, p0), edge);
    let norm_a_sq = dot(normal_a, normal_a);
    let norm_b_sq = dot(normal_b, normal_b);
    if !norm_a_sq.is_finite()
        || !norm_b_sq.is_finite()
        || norm_a_sq <= MIN_NORMAL_SQ
        || norm_b_sq <= MIN_NORMAL_SQ
    {
        return None;
    }
    let angle = dot(cross(normal_a, normal_b), edge_hat).atan2(dot(normal_a, normal_b));
    let weighted_a = scale(normal_a, 1.0 / norm_a_sq);
    let weighted_b = scale(normal_b, 1.0 / norm_b_sq);
    let grad_k = scale(weighted_a, -edge_len);
    let grad_l = scale(weighted_b, -edge_len);
    let along_k = dot(sub(pk, p1), edge_hat);
    let along_l = dot(sub(pl, p1), edge_hat);
    let grad_0 = add(scale(weighted_a, -along_k), scale(weighted_b, -along_l));
    let grad_1 = scale(add(add(grad_0, grad_k), grad_l), -1.0);
    Some(HingeState {
        angle,
        gradients: [grad_0, grad_1, grad_k, grad_l],
    })
}

/// One hinge `[a, b, c, d]` for every interior edge of a triangle list.
///
/// An edge shared by exactly two triangles becomes the hinge `[a, b, c, d]`:
/// `(a, b)` is the shared edge, directed as it appears in the first triangle
/// that holds it, `c` is the vertex opposite to the edge in that triangle and
/// `d` is the vertex opposite to it in the other triangle.  For consistently
/// wound triangles `c` therefore lies to the left of `a -> b` and `d` to the
/// right, the faces `(a, b, c)` and `(b, a, d)` are the two triangles with
/// their own winding, and a flat mesh has signed dihedral angle `0` at every
/// hinge.  Boundary edges (one triangle) and non-manifold edges (more than two
/// triangles) get no hinge.  The hinges are ordered by the sorted vertex pair
/// of their edge, which makes the result deterministic.
pub(super) fn interior_edge_hinges(triangles: &[[usize; 3]]) -> Vec<[usize; 4]> {
    let mut edges: BTreeMap<(usize, usize), ([usize; 2], Vec<usize>)> = BTreeMap::new();
    for &[v0, v1, v2] in triangles {
        for (a, b, wing) in [(v0, v1, v2), (v1, v2, v0), (v2, v0, v1)] {
            edges
                .entry((a.min(b), a.max(b)))
                .or_insert_with(|| ([a, b], Vec::new()))
                .1
                .push(wing);
        }
    }
    edges
        .into_values()
        .filter_map(|(edge, wings)| match wings.as_slice() {
            [c, d] => Some([edge[0], edge[1], *c, *d]),
            _ => None,
        })
        .collect()
}

/// Map an angle to `(-pi, pi]`.
///
/// Angles already in that range are returned unchanged (bit for bit), so small
/// angles keep their full precision.
pub(super) fn wrap_angle(angle: f64) -> f64 {
    if angle > -PI && angle <= PI {
        return angle;
    }
    let wrapped = angle - TAU * ((angle + PI) / TAU).floor();
    if wrapped <= -PI { PI } else { wrapped }
}

/// XPBD compliance `alpha_tilde = 1 / (k dt^2)` of a bending constraint with
/// stiffness `k` (N m / rad) for the time step `dt` (s).
///
/// Returns `None` when the constraint exerts no force: a stiffness that is not
/// strictly positive and finite means "no bending resistance", and a time step
/// that is not strictly positive and finite has no compliance.
pub(super) fn bending_compliance(stiffness: f64, dt: f64) -> Option<f64> {
    let usable = stiffness.is_finite() && stiffness > 0.0 && dt.is_finite() && dt > 0.0;
    usable.then(|| 1.0 / (stiffness * dt * dt))
}

/// One XPBD projection step of the dihedral constraint
/// `C = wrap(theta - rest_angle)` on the hinge `[p0, p1, pk, pl]`.
///
/// `lambda` is the multiplier accumulated for this constraint during the
/// current time step (`0` before the first iteration) and `alpha_tilde` the
/// time-scaled compliance (see [`bending_compliance`]; `0` makes the constraint
/// rigid).  The multiplier increment is
///
/// ```text
/// d_lambda = (-C - alpha_tilde * lambda) / (sum_i w_i |grad_i|^2 + alpha_tilde)
/// ```
///
/// and particle `i` moves by `w_i * grad_i * d_lambda`, where `w_i` is its
/// inverse mass.  Because the gradients sum to zero and carry no torque,
/// `sum_i m_i dx_i = 0` and `sum_i x_i x (m_i dx_i) = 0`.
///
/// Returns `None` when the hinge is degenerate, an input is not finite, or no
/// particle of the hinge can move (all inverse masses are zero).
pub(super) fn project_dihedral_xpbd(
    positions: [[f64; 3]; 4],
    inv_masses: [f64; 4],
    rest_angle: f64,
    alpha_tilde: f64,
    lambda: f64,
) -> Option<DihedralProjection> {
    let state = hinge_state(positions[0], positions[1], positions[2], positions[3])?;
    let violation = wrap_angle(state.angle - rest_angle);
    if !violation.is_finite() || !alpha_tilde.is_finite() || !lambda.is_finite() {
        return None;
    }
    let generalised_inverse_mass: f64 = state
        .gradients
        .iter()
        .zip(inv_masses)
        .map(|(gradient, inv_mass)| inv_mass * dot(*gradient, *gradient))
        .sum();
    if generalised_inverse_mass <= MIN_GENERALISED_INVERSE_MASS {
        return None;
    }
    let delta_lambda =
        (-violation - alpha_tilde * lambda) / (generalised_inverse_mass + alpha_tilde);
    let mut corrections = [[0.0_f64; 3]; 4];
    for ((correction, gradient), inv_mass) in
        corrections.iter_mut().zip(state.gradients).zip(inv_masses)
    {
        *correction = scale(gradient, inv_mass * delta_lambda);
    }
    Some(DihedralProjection {
        delta_lambda,
        corrections,
    })
}
