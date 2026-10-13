// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Fold angle of a bending hinge and its gradient.
//!
//! A hinge consists of the two triangles `(x0, x1, x2)` and `(x1, x0, x3)` that
//! share the edge `x0 x1`; `x2` and `x3` are the *wing* vertices.  With
//! `e = x1 - x0` the face normals
//!
//! ```text
//! n_a = e x (x2 - x0)          n_b = (x3 - x0) x e = (x0 - x1) x (x3 - x1)
//! ```
//!
//! point to the same side of a flat hinge, whose wings lie on opposite sides of
//! the edge.  The *signed fold angle* is the rotation about `e` that takes the
//! unit normal of the first face onto the unit normal of the second:
//!
//! ```text
//! phi = atan2( (n_a x n_b) . e_hat , n_a . n_b ),     phi in (-pi, pi]
//! ```
//!
//! It is `0` for a flat hinge, negative when both wings are raised towards the
//! side of `n_a` and positive when they are lowered, and `+-pi` when the wings
//! lie on top of each other.  Because `atan2` is used instead of `acos`, the
//! angle keeps full precision close to the flat state, which is the rest state
//! of a cloth.
//!
//! The gradient follows from rotating a wing about the hinge line: moving `x2`
//! by `t` along `n_a_hat` rotates the first face by `t / h_a`, where `h_a` is the
//! distance of `x2` from the hinge line, and the second wing behaves alike.  The
//! reactions that keep the net force and the net torque at zero act on the hinge
//! endpoints, split by the barycentric coordinate of the wing's projection onto
//! the hinge line:
//!
//! ```text
//! grad_x2 = -n_a_hat / h_a                   grad_x3 = -n_b_hat / h_b
//! grad_x0 = (1 - s_a) n_a_hat / h_a + (1 - s_b) n_b_hat / h_b
//! grad_x1 =      s_a  n_a_hat / h_a +      s_b  n_b_hat / h_b
//! ```
//!
//! with `s = (x_wing - x0) . e / |e|^2`.  The four gradients sum to zero
//! (translation invariance) and have no net torque (rotation invariance).

use super::types::BendingConstraint;
use oxiphysics_core::math::{Real, Vec3};

/// Shortest hinge edge (m) that still defines a fold axis.
const MIN_HINGE_LENGTH: Real = 1.0e-12;

/// Smallest accepted distance of a wing vertex from the hinge line, relative to
/// the length of the hinge edge.  A wing closer than this lies on the hinge line
/// and its face has no usable normal.
const MIN_RELATIVE_WING_HEIGHT: Real = 1.0e-9;

/// Signed fold angle of a hinge together with its gradient.
#[derive(Debug, Clone, Copy)]
pub(super) struct HingeFold {
    /// Signed fold angle `phi` in `(-pi, pi]`, `0` for a flat hinge.
    pub(super) angle: Real,
    /// `d phi / d x_k` for the hinge particles `[x0, x1, x2, x3]`.
    pub(super) gradient: [Vec3; 4],
}

impl BendingConstraint {
    /// Signed fold angle and gradient of the hinge `(x0, x1; x2, x3)`.
    ///
    /// Returns `None` when the hinge is degenerate: a zero-length or non-finite
    /// edge, or a wing vertex that lies on the hinge line.
    pub(super) fn hinge_fold(x0: &Vec3, x1: &Vec3, x2: &Vec3, x3: &Vec3) -> Option<HingeFold> {
        let edge = x1 - x0;
        let edge_len = edge.norm();
        let edge_usable = edge_len.is_finite() && edge_len > MIN_HINGE_LENGTH;
        if !edge_usable {
            return None;
        }
        let to_x2 = x2 - x0;
        let to_x3 = x3 - x0;
        let normal_a = edge.cross(&to_x2);
        let normal_b = to_x3.cross(&edge);
        let len_a = normal_a.norm();
        let len_b = normal_b.norm();
        let min_len = MIN_RELATIVE_WING_HEIGHT * edge_len * edge_len;
        let wings_usable =
            len_a.is_finite() && len_b.is_finite() && len_a > min_len && len_b > min_len;
        if !wings_usable {
            return None;
        }
        let unit_a = normal_a / len_a;
        let unit_b = normal_b / len_b;
        let axis = edge / edge_len;
        let angle = unit_a.cross(&unit_b).dot(&axis).atan2(unit_a.dot(&unit_b));

        let edge_len_sq = edge_len * edge_len;
        let along_a = to_x2.dot(&edge) / edge_len_sq;
        let along_b = to_x3.dot(&edge) / edge_len_sq;
        let push_a = unit_a * (edge_len / len_a);
        let push_b = unit_b * (edge_len / len_b);
        let gradient = [
            push_a * (1.0 - along_a) + push_b * (1.0 - along_b),
            push_a * along_a + push_b * along_b,
            -push_a,
            -push_b,
        ];
        Some(HingeFold { angle, gradient })
    }
}
