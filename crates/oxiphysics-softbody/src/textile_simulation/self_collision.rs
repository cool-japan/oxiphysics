// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Vertex-face self-collision resolution of [`FabricCollision`].
//!
//! Every vertex is tested against the triangles it does not belong to. The
//! exhaustive scan costs `iterations * vertices * triangles` narrow-phase tests
//! per call; a uniform grid over the triangles ([`TriangleGrid`]) restricts
//! each vertex to the triangles in its neighbourhood, which makes the cost
//! proportional to the number of vertices for a mesh of roughly uniform
//! resolution. The corrections are the same as those of the exhaustive scan.

use oxiphysics_core::math::{Real, Vec3};

use super::FabricCollision;
use super::triangle_grid::TriangleGrid;

/// Work done by one call to [`FabricCollision::resolve_self_collisions`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SelfCollisionStats {
    /// Vertex-triangle pairs that reached the narrow-phase test
    /// ([`FabricCollision::vertex_face_test`]).
    ///
    /// The exhaustive scan would test `iterations * vertices * triangles`
    /// pairs.
    pub pairs_tested: usize,
    /// Penetrating pairs that were pushed apart.
    pub contacts_resolved: usize,
    /// Number of times the broad-phase grid had to be rebuilt in the middle of
    /// the call because particles had moved farther than its safety margin.
    /// After too many rebuilds the grid is dropped and the rest of the call
    /// scans all triangles for every vertex.
    pub grid_rebuilds: usize,
}

impl FabricCollision {
    /// Resolve vertex-face penetrations using position-based correction.
    ///
    /// For each of [`FabricCollision::collision_iterations`] iterations every
    /// vertex is tested against the triangles it is not part of, in ascending
    /// triangle order. A vertex that is closer than
    /// [`FabricCollision::collision_thickness`] to the plane of a triangle it
    /// projects into is pushed out along the triangle normal, and the triangle
    /// corners are pushed the opposite way, distributed by inverse mass. Pinned
    /// particles do not move. Iterations stop early once one of them finds
    /// nothing to correct, since the following ones would see the same state.
    ///
    /// A uniform grid over the triangles limits the triangles each vertex is
    /// tested against to its neighbourhood, so the cost grows linearly with the
    /// mesh size rather than with the product of vertex and triangle counts.
    /// The result is identical to testing every vertex against every triangle.
    pub fn resolve_self_collisions(&mut self) {
        self.resolve_self_collisions_with_stats();
    }

    /// Like [`FabricCollision::resolve_self_collisions`], and also reports how
    /// much work the call did.
    pub fn resolve_self_collisions_with_stats(&mut self) -> SelfCollisionStats {
        let mut stats = SelfCollisionStats::default();
        let thickness = self.collision_thickness;
        // `|distance| < thickness` can never hold for a non-positive or NaN
        // thickness, so nothing can be in contact.
        if self.collision_iterations == 0
            || self.particles.is_empty()
            || self.triangles.is_empty()
            || thickness.is_nan()
            || thickness <= 0.0
        {
            return stats;
        }
        let mut grid = TriangleGrid::build(&self.particles, &self.triangles, thickness);
        for _ in 0..self.collision_iterations {
            let resolved_before = stats.contacts_resolved;
            for vertex in 0..self.particles.len() {
                self.resolve_vertex_contacts(vertex, &mut grid, &mut stats);
            }
            if stats.contacts_resolved == resolved_before {
                break;
            }
        }
        stats
    }

    /// Test `vertex` against its candidate triangles and push it out of the
    /// ones it penetrates.
    ///
    /// The candidates are visited in ascending triangle order and every test
    /// sees the positions left by the previous corrections, exactly like the
    /// exhaustive scan. When a correction moves a particle beyond the grid's
    /// safety margin the grid is rebuilt and the scan of this vertex resumes
    /// with the next triangle.
    fn resolve_vertex_contacts(
        &mut self,
        vertex: usize,
        grid: &mut TriangleGrid,
        stats: &mut SelfCollisionStats,
    ) {
        let mut next_triangle = 0;
        loop {
            let mut grid_is_stale = false;
            let candidates = grid.candidates(vertex);
            let first = candidates.partition_point(|&triangle| triangle < next_triangle);
            for &triangle in &candidates[first..] {
                next_triangle = triangle + 1;
                let corners = self.triangles[triangle];
                stats.pairs_tested += 1;
                let Some((depth, normal)) = self.vertex_face_test(vertex, corners) else {
                    continue;
                };
                if !self.push_vertex_off_triangle(vertex, corners, depth, normal) {
                    continue;
                }
                stats.contacts_resolved += 1;
                let moved = [vertex, corners[0], corners[1], corners[2]];
                if grid.has_drifted(&self.particles, &moved) {
                    grid_is_stale = true;
                    break;
                }
            }
            if !grid_is_stale {
                return;
            }
            if grid.rebuild(&self.particles, &self.triangles, self.collision_thickness) {
                stats.grid_rebuilds += 1;
            }
        }
    }

    /// Separate `vertex` and the triangle `corners` along `normal` by `depth`,
    /// distributed by inverse mass. Returns `false` (and moves nothing) when
    /// all the particles involved are effectively immovable.
    fn push_vertex_off_triangle(
        &mut self,
        vertex: usize,
        corners: [usize; 3],
        depth: Real,
        normal: Vec3,
    ) -> bool {
        let w_v = self.particles[vertex].inv_mass;
        let w_tri: Real = corners.iter().map(|&i| self.particles[i].inv_mass).sum();
        let total_w = w_v + w_tri / 3.0;
        if total_w < 1e-12 {
            return false;
        }
        let correction = depth / total_w;
        if !self.particles[vertex].is_pinned() {
            self.particles[vertex].position += normal * (correction * w_v);
        }
        for &corner in corners.iter() {
            if !self.particles[corner].is_pinned() {
                let w_i = self.particles[corner].inv_mass;
                self.particles[corner].position -= normal * (correction * w_i / 3.0);
            }
        }
        true
    }
}
