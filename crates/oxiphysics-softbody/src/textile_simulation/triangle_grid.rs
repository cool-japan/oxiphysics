// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Uniform grid over the triangles of a fabric mesh: the broad phase of the
//! vertex-face self-collision test of [`super::FabricCollision`].
//!
//! # Contract
//!
//! [`TriangleGrid::candidates`] returns, for a vertex, triangle indices in
//! **ascending order** that contain every triangle the vertex can be in contact
//! with. The exhaustive scan visits the triangles in the same ascending order,
//! so running the narrow phase on the candidates only produces exactly the same
//! corrections, in the same order, as running it on all triangles.
//!
//! The guarantee rests on a *safety margin* (`skin`). The grid is built from the
//! particle positions at one instant (the *reference* positions):
//!
//! * a triangle is stored in every cell its reference bounding box overlaps
//!   after being inflated by `reach = thickness + 2 * skin`, and
//! * a vertex looks up the cell of its **reference** position.
//!
//! A vertex can only be in contact with a triangle when it lies within
//! `thickness` of a point of that triangle. If no particle has moved farther
//! than `skin` from its reference position, the vertex and every corner of the
//! triangle have moved by at most `skin` each, so such a pair still has its
//! reference vertex position inside the inflated reference box of the triangle,
//! i.e. in one of the cells the triangle is stored in. As soon as a particle
//! exceeds the margin ([`TriangleGrid::has_drifted`]) the caller rebuilds the
//! grid ([`TriangleGrid::rebuild`]).
//!
//! The margin is the larger of the collision thickness (one correction never
//! moves a particle by more than that) and a fraction of the mean triangle
//! size, so that contact-rich meshes, where particles are pushed many times,
//! rarely have to rebuild. The number of rebuilds per grid is capped; once the
//! cap is hit the grid is dropped in favour of the exhaustive scan, which keeps
//! the worst case within a small factor of the exhaustive cost.
//!
//! The narrow phase itself is evaluated in floating point; the box inflation
//! carries a relative slack of a few dozen ulps on top of `reach`, which is
//! far below the margin that `thickness` and `skin` provide.
//!
//! When the geometry cannot be indexed (non-finite positions, a triangle that
//! refers to a missing particle, cell indices that do not fit, a non-positive or
//! infinite thickness) the grid degrades to the exhaustive scan: every vertex
//! gets every triangle as candidate. A uniform grid also degrades towards the
//! exhaustive scan for meshes whose largest triangle is much bigger than the
//! typical one, because the cell size follows the largest triangle.
//!
//! The grid is private to this module rather than built on the point-based
//! spatial hash of `oxiphysics-geometry`: triangles are boxes that must be
//! stored in every cell they overlap, the candidate lists must stay in
//! ascending triangle order, and the conservativeness argument above needs
//! control over the exact cell arithmetic.

use std::collections::HashMap;

use oxiphysics_core::math::{Real, Vec3};

use super::FabricParticle;

/// Integer coordinates of a grid cell.
type CellKey = [i64; 3];

/// Cell indices must stay below this magnitude (comfortably below `2^52`) so
/// that the `Real -> i64` conversion is exact and cannot saturate.
const MAX_CELL_INDEX: Real = 4.0e15;

/// A triangle that would have to be stored in more cells than this is not worth
/// indexing. With a cell edge of at least the largest triangle extent plus the
/// inflation a triangle overlaps at most two cells per axis, so this is only
/// reached through degenerate floating-point input.
const MAX_CELLS_PER_TRIANGLE: usize = 64;

/// Relative slack, in units of `f64::EPSILON`, added to the box inflation to
/// absorb the rounding of the box and cell arithmetic.
const ROUNDING_SLACK: Real = 64.0 * Real::EPSILON;

/// Fraction of the mean triangle extent used as safety margin when it exceeds
/// the collision thickness. A margin of a quarter of a triangle keeps the number
/// of candidates per vertex in the tens while making rebuilds rare.
const SKIN_FRACTION: Real = 0.25;

/// Rebuilds a grid may perform regardless of the mesh size ...
pub(super) const REBUILD_BUDGET_BASE: usize = 4;

/// ... plus one per this many particles. A rebuild costs about as much as ten
/// narrow-phase tests per triangle, so the budget caps the total rebuild work
/// at roughly one exhaustive pass (`vertices * triangles` tests).
pub(super) const REBUILD_BUDGET_PARTICLES: usize = 16;

/// Broad-phase index of the triangles of a fabric mesh, see the
/// [module documentation](self).
#[derive(Debug, Clone)]
pub(super) struct TriangleGrid {
    /// Every triangle index in ascending order: the candidates of the
    /// exhaustive scan.
    all_triangles: Vec<usize>,
    /// `None` when the mesh could not be indexed or the grid was given up.
    index: Option<GridIndex>,
    /// Rebuilds still allowed before the exhaustive scan takes over.
    rebuilds_left: usize,
}

/// The spatial index proper.
#[derive(Debug, Clone)]
struct GridIndex {
    /// Largest distance a particle may move from its reference position before
    /// the candidate lists stop being conservative.
    skin: Real,
    /// Edge length of the cubic cells.
    cell_size: Real,
    /// Triangle indices (ascending) per cell.
    cells: HashMap<CellKey, Vec<usize>>,
    /// Particle positions the grid was built from; every one of them lies in an
    /// indexable cell.
    reference: Vec<Vec3>,
}

impl TriangleGrid {
    /// Index `triangles` for contact tests with a collision `thickness`.
    ///
    /// The current particle positions become the reference positions.
    pub(super) fn build(
        particles: &[FabricParticle],
        triangles: &[[usize; 3]],
        thickness: Real,
    ) -> Self {
        Self {
            all_triangles: (0..triangles.len()).collect(),
            index: GridIndex::build(particles, triangles, thickness),
            rebuilds_left: REBUILD_BUDGET_BASE + particles.len() / REBUILD_BUDGET_PARTICLES,
        }
    }

    /// Triangle indices (ascending) that `vertex` may be in contact with as
    /// long as no particle has drifted farther than the safety margin
    /// ([`TriangleGrid::has_drifted`]).
    pub(super) fn candidates(&self, vertex: usize) -> &[usize] {
        match &self.index {
            Some(index) => index.candidates(vertex),
            None => &self.all_triangles,
        }
    }

    /// `true` if any of the particles `indices` has moved farther from its
    /// reference position than the safety margin, so that the candidate lists
    /// may no longer be complete and the grid has to be rebuilt.
    ///
    /// Always `false` for the exhaustive fallback, whose lists are complete by
    /// construction. A particle whose position became NaN never counts as
    /// drifted: it cannot be in contact with anything.
    pub(super) fn has_drifted(&self, particles: &[FabricParticle], indices: &[usize]) -> bool {
        self.index
            .as_ref()
            .is_some_and(|index| index.has_drifted(particles, indices))
    }

    /// Re-index the mesh from the current particle positions.
    ///
    /// Returns `true` if the grid was rebuilt. Returns `false` when the rebuild
    /// budget is used up or the mesh can no longer be indexed; the grid then
    /// falls back to the exhaustive scan for good.
    pub(super) fn rebuild(
        &mut self,
        particles: &[FabricParticle],
        triangles: &[[usize; 3]],
        thickness: Real,
    ) -> bool {
        if self.rebuilds_left == 0 {
            self.index = None;
            return false;
        }
        self.rebuilds_left -= 1;
        self.index = GridIndex::build(particles, triangles, thickness);
        self.index.is_some()
    }

    /// `true` if a spatial index exists (otherwise every vertex is tested
    /// against every triangle).
    #[cfg(test)]
    pub(super) fn is_indexed(&self) -> bool {
        self.index.is_some()
    }

    /// Edge length of the grid cells, `None` for the exhaustive fallback.
    #[cfg(test)]
    pub(super) fn cell_size(&self) -> Option<Real> {
        self.index.as_ref().map(|index| index.cell_size)
    }
}

impl GridIndex {
    fn build(
        particles: &[FabricParticle],
        triangles: &[[usize; 3]],
        thickness: Real,
    ) -> Option<Self> {
        if !thickness.is_finite() || thickness <= 0.0 {
            return None;
        }
        if particles
            .iter()
            .any(|p| !p.position.iter().all(|c| c.is_finite()))
        {
            return None;
        }
        if triangles.is_empty() {
            return None;
        }

        let mut boxes: Vec<([Real; 3], [Real; 3])> = Vec::with_capacity(triangles.len());
        let mut max_extent: Real = 0.0;
        let mut total_extent: Real = 0.0;
        for tri in triangles {
            let corners = [
                particles.get(tri[0])?.position,
                particles.get(tri[1])?.position,
                particles.get(tri[2])?.position,
            ];
            let lo: [Real; 3] = std::array::from_fn(|axis| {
                corners
                    .iter()
                    .map(|corner| corner[axis])
                    .fold(Real::INFINITY, Real::min)
            });
            let hi: [Real; 3] = std::array::from_fn(|axis| {
                corners
                    .iter()
                    .map(|corner| corner[axis])
                    .fold(Real::NEG_INFINITY, Real::max)
            });
            let extent = lo
                .iter()
                .zip(hi.iter())
                .map(|(low, high)| high - low)
                .fold(0.0, Real::max);
            max_extent = max_extent.max(extent);
            total_extent += extent;
            boxes.push((lo, hi));
        }
        let mean_extent = total_extent / triangles.len() as Real;
        let skin = thickness.max(SKIN_FRACTION * mean_extent);
        let reach = thickness + 2.0 * skin;

        // A cell at least as wide as the largest inflated triangle box means a
        // triangle overlaps at most two cells per axis.
        let cell_size = max_extent + 2.0 * reach;
        if !cell_size.is_finite() || cell_size <= 0.0 {
            return None;
        }

        let mut cells: HashMap<CellKey, Vec<usize>> = HashMap::new();
        for (triangle, (lo, hi)) in boxes.iter().enumerate() {
            let mut first = [0_i64; 3];
            let mut last = [0_i64; 3];
            let mut cell_count: usize = 1;
            for axis in 0..3 {
                first[axis] = inflated_cell(lo[axis], -1.0, reach, cell_size)?;
                last[axis] = inflated_cell(hi[axis], 1.0, reach, cell_size)?;
                let span = usize::try_from(last[axis] - first[axis] + 1).ok()?;
                cell_count = cell_count.checked_mul(span)?;
            }
            if cell_count > MAX_CELLS_PER_TRIANGLE {
                return None;
            }
            for x in first[0]..=last[0] {
                for y in first[1]..=last[1] {
                    for z in first[2]..=last[2] {
                        cells.entry([x, y, z]).or_default().push(triangle);
                    }
                }
            }
        }

        if particles
            .iter()
            .any(|p| cell_key(&p.position, cell_size).is_none())
        {
            return None;
        }
        Some(Self {
            skin,
            cell_size,
            cells,
            reference: particles.iter().map(|p| p.position).collect(),
        })
    }

    fn candidates(&self, vertex: usize) -> &[usize] {
        self.reference
            .get(vertex)
            .and_then(|position| cell_key(position, self.cell_size))
            .and_then(|key| self.cells.get(&key))
            .map_or(&[], Vec::as_slice)
    }

    fn has_drifted(&self, particles: &[FabricParticle], indices: &[usize]) -> bool {
        indices
            .iter()
            .any(|&i| match (particles.get(i), self.reference.get(i)) {
                (Some(particle), Some(reference)) => {
                    (particle.position - reference).norm() > self.skin
                }
                _ => false,
            })
    }
}

/// Cell index of `value` moved by `direction * reach` (`direction` is `-1` for
/// the lower and `+1` for the upper end of a box), with the rounding slack
/// applied.
fn inflated_cell(value: Real, direction: Real, reach: Real, cell_size: Real) -> Option<i64> {
    let slack = reach + ROUNDING_SLACK * (value.abs() + reach);
    cell_coord(value + direction * slack, cell_size)
}

/// Cell index of one coordinate, `None` when it is not finite or too large to
/// index.
fn cell_coord(value: Real, cell_size: Real) -> Option<i64> {
    let cell = (value / cell_size).floor();
    if cell.abs() < MAX_CELL_INDEX {
        Some(cell as i64)
    } else {
        None
    }
}

/// Cell of a position, `None` when it cannot be indexed.
fn cell_key(position: &Vec3, cell_size: Real) -> Option<CellKey> {
    Some([
        cell_coord(position.x, cell_size)?,
        cell_coord(position.y, cell_size)?,
        cell_coord(position.z, cell_size)?,
    ])
}
