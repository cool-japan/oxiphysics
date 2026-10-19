// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact edge-edge continuous collision detection for deformable meshes.
//!
//! Every vertex of an edge moves linearly over the step, so the four end points
//! of an edge pair become coplanar at the roots of a cubic in time (Provot 1997,
//! Bridson et al. 2002). The root finding, the geometric validation of each root
//! and the fallback for pairs that only graze the contact shell are provided by
//! [`ccd_edge_edge_analytic`]; this module adds the mesh-level bookkeeping
//! (index validation, broad-phase pruning), the rounding floor of the contact
//! radius, the rescaling of every pair to a unit-length local frame (the
//! shared routines use absolute tolerances) and the contact normal.

use super::{add3, cross3, dot3, edge_edge_distance, len3, normalize3, scale3, sub3};
use crate::deformable_collision::ccd_edge_edge_analytic;
use crate::soft_body_collision::EdgeCcdResult;

/// Relative floor of the contact radius used to validate a contact.
///
/// Two edges that cross exactly are at distance zero only in exact arithmetic:
/// evaluated at the rounded root time of the coplanarity cubic the distance
/// comes out around `1e-16` times the magnitude of the coordinates, so a
/// contact radius of exactly zero would reject the root of every generic
/// crossing. The radius of a pair is therefore never taken smaller than this
/// fraction of `1 + max |coordinate|` of its end points. The floor is orders of
/// magnitude below any physical contact radius, which it leaves unchanged.
const RADIUS_FLOOR: f64 = 1e-10;

/// An edge whose two end points move linearly over one step.
#[derive(Clone, Copy)]
struct MovingEdge {
    /// End point positions at the start of the step.
    start: [[f64; 3]; 2],
    /// End point positions at the end of the step.
    end: [[f64; 3]; 2],
}

impl MovingEdge {
    /// Build the edge `(i0, i1)`; `None` if an index is outside `verts0` or `verts1`.
    fn new(i0: usize, i1: usize, verts0: &[[f64; 3]], verts1: &[[f64; 3]]) -> Option<Self> {
        Some(MovingEdge {
            start: [*verts0.get(i0)?, *verts0.get(i1)?],
            end: [*verts1.get(i0)?, *verts1.get(i1)?],
        })
    }

    /// End points at normalised step time `t`.
    fn at(&self, t: f64) -> [[f64; 3]; 2] {
        let lerp = |a: [f64; 3], b: [f64; 3]| add3(scale3(a, 1.0 - t), scale3(b, t));
        [
            lerp(self.start[0], self.end[0]),
            lerp(self.start[1], self.end[1]),
        ]
    }

    /// Velocity (per unit step time) of the point at parameter `s` along the edge.
    fn velocity_at(&self, s: f64) -> [f64; 3] {
        let v0 = sub3(self.end[0], self.start[0]);
        let v1 = sub3(self.end[1], self.start[1]);
        add3(scale3(v0, 1.0 - s), scale3(v1, s))
    }

    /// Axis-aligned box enclosing the whole trajectory of the edge.
    ///
    /// Each end point moves on a straight segment, so the box of the four
    /// positions at `t = 0` and `t = 1` contains the edge at every time.
    fn swept_bounds(&self) -> ([f64; 3], [f64; 3]) {
        let mut lo = self.start[0];
        let mut hi = self.start[0];
        for p in self.start.iter().chain(self.end.iter()) {
            for ((l, h), &c) in lo.iter_mut().zip(hi.iter_mut()).zip(p.iter()) {
                *l = l.min(c);
                *h = h.max(c);
            }
        }
        (lo, hi)
    }

    /// Largest absolute coordinate over the four end point positions.
    fn magnitude(&self) -> f64 {
        self.start
            .iter()
            .chain(self.end.iter())
            .flatten()
            .fold(0.0_f64, |largest, c| largest.max(c.abs()))
    }

    /// Length of the longer of the two edge vectors (at the start and at the
    /// end of the step).
    fn longest_vector(&self) -> f64 {
        len3(sub3(self.start[1], self.start[0])).max(len3(sub3(self.end[1], self.end[0])))
    }

    /// The edge with every position mapped to `(position - origin) * inv_scale`.
    fn rescaled(&self, origin: [f64; 3], inv_scale: f64) -> MovingEdge {
        let map = |p: [f64; 3]| scale3(sub3(p, origin), inv_scale);
        MovingEdge {
            start: [map(self.start[0]), map(self.start[1])],
            end: [map(self.end[0]), map(self.end[1])],
        }
    }
}

/// A moving edge together with the data the broad phase needs.
struct PreparedEdge {
    edge: MovingEdge,
    /// Lower corner of the box enclosing the whole trajectory.
    lo: [f64; 3],
    /// Upper corner of the box enclosing the whole trajectory.
    hi: [f64; 3],
    /// Largest absolute coordinate of the end points.
    magnitude: f64,
}

impl PreparedEdge {
    /// Build the edge `(i0, i1)`; `None` if an index is outside `verts0` or `verts1`.
    fn new(i0: usize, i1: usize, verts0: &[[f64; 3]], verts1: &[[f64; 3]]) -> Option<Self> {
        let edge = MovingEdge::new(i0, i1, verts0, verts1)?;
        let (lo, hi) = edge.swept_bounds();
        Some(PreparedEdge {
            magnitude: edge.magnitude(),
            edge,
            lo,
            hi,
        })
    }
}

/// Contact radius used to validate the contacts of a pair of edges whose
/// coordinates are at most `magnitude` in absolute value: `edge_radius`, raised
/// to the rounding floor (see [`RADIUS_FLOOR`]) when it is smaller. Coordinates
/// that are not finite leave the radius unchanged.
fn contact_radius(edge_radius: f64, magnitude: f64) -> f64 {
    if magnitude.is_finite() {
        edge_radius.max(RADIUS_FLOOR * (1.0 + magnitude))
    } else {
        edge_radius
    }
}

/// Closest-point data of an edge pair frozen at one time.
struct PairGeometry {
    /// Closest point on edge A minus closest point on edge B.
    separation: [f64; 3],
    /// Direction vector of edge A (end minus start).
    dir_a: [f64; 3],
    /// Direction vector of edge B (end minus start).
    dir_b: [f64; 3],
    /// Parameter of the closest point on edge A.
    s: f64,
    /// Parameter of the closest point on edge B.
    t: f64,
}

fn pair_geometry(a: &MovingEdge, b: &MovingEdge, time: f64) -> PairGeometry {
    let [pa0, pa1] = a.at(time);
    let [pb0, pb1] = b.at(time);
    let dir_a = sub3(pa1, pa0);
    let dir_b = sub3(pb1, pb0);
    let (_sq_dist, s, t) = edge_edge_distance(pa0, dir_a, pb0, dir_b);
    let closest_a = add3(pa0, scale3(dir_a, s));
    let closest_b = add3(pb0, scale3(dir_b, t));
    PairGeometry {
        separation: sub3(closest_a, closest_b),
        dir_a,
        dir_b,
        s,
        t,
    }
}

/// Contact normal (pointing from edge B towards edge A) of an edge pair at the
/// time of impact `toi`.
///
/// When the edges are still apart at `toi` (a graze of the contact shell) the
/// normal is the unit separation of the closest points. At an exact
/// coplanarity root the closest points coincide and the separation carries no
/// direction, so the common perpendicular of the two edges is used instead,
/// oriented towards the side on which edge A started the step. For parallel
/// touching edges the pre-contact separation, then the approach direction,
/// provide the normal.
fn contact_normal(a: &MovingEdge, b: &MovingEdge, toi: f64) -> [f64; 3] {
    let at_toi = pair_geometry(a, b, toi);
    let length_scale = len3(at_toi.dir_a) + len3(at_toi.dir_b);
    let degenerate = 1e-9 * length_scale;
    let gap = len3(at_toi.separation);
    if gap > degenerate {
        return scale3(at_toi.separation, 1.0 / gap);
    }

    let reference = pair_geometry(a, b, 0.0).separation;
    let reference_len = len3(reference);
    let relative_velocity = sub3(a.velocity_at(at_toi.s), b.velocity_at(at_toi.t));
    let cross = cross3(at_toi.dir_a, at_toi.dir_b);
    let cross_len = len3(cross);
    if cross_len > 1e-12 * len3(at_toi.dir_a) * len3(at_toi.dir_b) {
        let normal = scale3(cross, 1.0 / cross_len);
        let orientation = if reference_len > degenerate {
            dot3(normal, reference)
        } else {
            -dot3(normal, relative_velocity)
        };
        return if orientation < 0.0 {
            scale3(normal, -1.0)
        } else {
            normal
        };
    }
    if reference_len > degenerate {
        return scale3(reference, 1.0 / reference_len);
    }
    normalize3(scale3(relative_velocity, -1.0))
}

/// Exact edge-edge CCD of two edge sets over one step.
///
/// Edges with a vertex index outside the position arrays of their mesh are
/// skipped. The reported `edge_a` / `edge_b` are indices into `edges_a` /
/// `edges_b`. A negative or NaN `edge_radius` can never be reached, so nothing
/// is reported for it; a radius of zero detects exact crossings (see
/// [`RADIUS_FLOOR`]). Each pair is evaluated in a local frame whose unit of
/// length is the longest edge vector of the pair, so the result does not depend
/// on the units of the mesh.
pub(super) fn edge_edge_ccd(
    edge_radius: f64,
    edges_a: &[(usize, usize)],
    verts0_a: &[[f64; 3]],
    verts1_a: &[[f64; 3]],
    edges_b: &[(usize, usize)],
    verts0_b: &[[f64; 3]],
    verts1_b: &[[f64; 3]],
) -> Vec<EdgeCcdResult> {
    if edge_radius.is_nan() || edge_radius < 0.0 {
        return Vec::new();
    }
    let prepared_b: Vec<Option<PreparedEdge>> = edges_b
        .iter()
        .map(|&(b0, b1)| PreparedEdge::new(b0, b1, verts0_b, verts1_b))
        .collect();
    let mut results = Vec::new();
    for (ia, &(a0, a1)) in edges_a.iter().enumerate() {
        let Some(a) = PreparedEdge::new(a0, a1, verts0_a, verts1_a) else {
            continue;
        };
        for (ib, prepared) in prepared_b.iter().enumerate() {
            let Some(b) = prepared else {
                continue;
            };
            let radius = contact_radius(edge_radius, a.magnitude.max(b.magnitude));
            // Broad phase: two edges can only come within `radius` of each other
            // if their swept boxes, grown by `radius`, overlap.
            let separated =
                (0..3).any(|k| a.lo[k] > b.hi[k] + radius || b.lo[k] > a.hi[k] + radius);
            if separated {
                continue;
            }
            // Narrow phase in a local frame: the origin is moved to a vertex of
            // the pair and the longest edge vector becomes the unit of length.
            // Time of impact and normal do not depend on either, and the
            // distance and parallelism tolerances of the shared edge-edge
            // routines are absolute, so they only behave on well-scaled input.
            let length = a.edge.longest_vector().max(b.edge.longest_vector());
            let length = if length.is_finite() && length > 0.0 {
                length
            } else {
                1.0
            };
            let inv_length = 1.0 / length;
            let origin = a.edge.start[0];
            let local_a = a.edge.rescaled(origin, inv_length);
            let local_b = b.edge.rescaled(origin, inv_length);
            let Some(toi) = ccd_edge_edge_analytic(
                local_a.start[0],
                local_a.end[0],
                local_a.start[1],
                local_a.end[1],
                local_b.start[0],
                local_b.end[0],
                local_b.start[1],
                local_b.end[1],
                radius * inv_length,
            ) else {
                continue;
            };
            results.push(EdgeCcdResult {
                edge_a: ia,
                edge_b: ib,
                toi,
                normal: contact_normal(&local_a, &local_b, toi),
            });
        }
    }
    results
}
