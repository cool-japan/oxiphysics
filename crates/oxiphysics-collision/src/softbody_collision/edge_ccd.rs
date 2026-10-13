// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact edge-edge continuous collision detection for deformable meshes.
//!
//! Every vertex of an edge moves linearly over the step, so the four end points
//! of an edge pair become coplanar at the roots of a cubic in time (Provot 1997,
//! Bridson et al. 2002). The root finding, the geometric validation of each root
//! and the fallback for pairs that only graze the contact shell are provided by
//! [`ccd_edge_edge_analytic`]; this module adds the mesh-level bookkeeping
//! (index validation, broad-phase pruning) and the contact normal.

use super::{add3, cross3, dot3, edge_edge_distance, len3, normalize3, scale3, sub3};
use crate::deformable_collision::ccd_edge_edge_analytic;
use crate::soft_body_collision::EdgeCcdResult;

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

    /// Axis-aligned box enclosing the whole trajectory of the edge, grown by `pad`.
    ///
    /// Each end point moves on a straight segment, so the box of the four
    /// positions at `t = 0` and `t = 1` contains the edge at every time.
    fn swept_bounds(&self, pad: f64) -> ([f64; 3], [f64; 3]) {
        let mut lo = self.start[0];
        let mut hi = self.start[0];
        for p in self.start.iter().chain(self.end.iter()) {
            for ((l, h), &c) in lo.iter_mut().zip(hi.iter_mut()).zip(p.iter()) {
                *l = l.min(c);
                *h = h.max(c);
            }
        }
        (
            [lo[0] - pad, lo[1] - pad, lo[2] - pad],
            [hi[0] + pad, hi[1] + pad, hi[2] + pad],
        )
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
/// `edges_b`.
pub(super) fn edge_edge_ccd(
    edge_radius: f64,
    edges_a: &[(usize, usize)],
    verts0_a: &[[f64; 3]],
    verts1_a: &[[f64; 3]],
    edges_b: &[(usize, usize)],
    verts0_b: &[[f64; 3]],
    verts1_b: &[[f64; 3]],
) -> Vec<EdgeCcdResult> {
    let moving_b: Vec<Option<MovingEdge>> = edges_b
        .iter()
        .map(|&(b0, b1)| MovingEdge::new(b0, b1, verts0_b, verts1_b))
        .collect();
    let pad = edge_radius.max(0.0);
    let mut results = Vec::new();
    for (ia, &(a0, a1)) in edges_a.iter().enumerate() {
        let Some(edge_a) = MovingEdge::new(a0, a1, verts0_a, verts1_a) else {
            continue;
        };
        let (lo_a, hi_a) = edge_a.swept_bounds(pad);
        for (ib, moving) in moving_b.iter().enumerate() {
            let Some(edge_b) = moving else {
                continue;
            };
            // Broad phase: the swept boxes (grown by the contact radius) must overlap.
            let (lo_b, hi_b) = edge_b.swept_bounds(0.0);
            let separated = (0..3).any(|k| lo_a[k] > hi_b[k] || lo_b[k] > hi_a[k]);
            if separated {
                continue;
            }
            let Some(toi) = ccd_edge_edge_analytic(
                edge_a.start[0],
                edge_a.end[0],
                edge_a.start[1],
                edge_a.end[1],
                edge_b.start[0],
                edge_b.end[0],
                edge_b.start[1],
                edge_b.end[1],
                edge_radius,
            ) else {
                continue;
            };
            results.push(EdgeCcdResult {
                edge_a: ia,
                edge_b: ib,
                toi,
                normal: contact_normal(&edge_a, edge_b, toi),
            });
        }
    }
    results
}
