// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Approximate convex decomposition of triangle meshes.
//!
//! The algorithm is described on [`ConvexDecomposition`]; this module also holds
//! the QuickHull helpers shared with the rest of `compound_shapes`.

use super::{Aabb3, MeshShapeQuery, add3, cross3, dot3, norm3, scale3, sub3};
use crate::broadphase::BvhBroadphase;
use crate::contact_generation::clip_polygon_by_plane;
use oxiphysics_core::Aabb;
use oxiphysics_core::math::Vec3;
use oxiphysics_geometry::quickhull::ConvexHull3DVec;

/// A mesh triangle given by its three corner positions.
type Triangle = [[f64; 3]; 3];

/// Fewest subdivisions per hull-face edge used when probing the hull surface.
/// Three subdivisions already place a sample at every face centre and two on
/// every edge.
const MIN_SAMPLING: usize = 3;

/// Most subdivisions per hull-face edge used when probing the hull surface.
const MAX_SAMPLING: usize = 8;

/// Concavity below this fraction of the mesh diameter is numerical noise and
/// never triggers a split, even if the caller asks for a threshold of zero.
const RELATIVE_NOISE: f64 = 1e-9;

/// Cut positions closer than this fraction of the part's extent to its
/// boundary would produce a negligible sliver and are not considered.
const CUT_MARGIN: f64 = 1e-6;

/// A hull vertex that rises less than this fraction of the hull diagonal above
/// its neighbours along its averaged face normal is a point inside a face or on
/// an edge of the hull, not a corner.
const FLAT_VERTEX_TOLERANCE: f64 = 1e-9;

/// Convex hull of `points` computed with QuickHull.
///
/// Non-finite points are ignored. Returns `None` when the remaining points do
/// not span three dimensions (fewer than four points, or all coplanar).
pub(super) fn quickhull(points: &[[f64; 3]]) -> Option<ConvexHull3DVec> {
    let cloud: Vec<Vec3> = points
        .iter()
        .filter(|p| p.iter().all(|c| c.is_finite()))
        .map(|p| Vec3::new(p[0], p[1], p[2]))
        .collect();
    ConvexHull3DVec::build(&cloud)
}

/// The finite points of `points`, sorted lexicographically and deduplicated.
pub(super) fn sorted_unique_points(points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    let mut unique: Vec<[f64; 3]> = points
        .iter()
        .filter(|p| p.iter().all(|c| c.is_finite()))
        // Adding 0.0 turns -0.0 into +0.0 so equal points compare equal.
        .map(|p| [p[0] + 0.0, p[1] + 0.0, p[2] + 0.0])
        .collect();
    unique.sort_by(|a, b| {
        a[0].total_cmp(&b[0])
            .then(a[1].total_cmp(&b[1]))
            .then(a[2].total_cmp(&b[2]))
    });
    unique.dedup();
    unique
}

/// Corner points of the convex hull of `points`.
///
/// When the points do not span three dimensions (no volumetric hull exists) all
/// distinct finite points are returned instead; a support function evaluated
/// over them is exact.
pub(super) fn hull_vertices(points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    match quickhull(points) {
        Some(hull) => extreme_vertices(&hull),
        None => sorted_unique_points(points),
    }
}

/// The corner points of `hull`.
///
/// QuickHull keeps every point it was seeded with, even a seed that turns out
/// to lie inside a face or on an edge of the finished hull. Such a point adds
/// nothing to the shape, so it is dropped here: a vertex is a corner when the
/// direction that averages the normals of its faces is a direction in which
/// that vertex alone is the farthest, i.e. when every neighbouring vertex lies
/// lower along it by more than `FLAT_VERTEX_TOLERANCE` times the hull diagonal.
/// For a vertex inside a face all neighbours lie in the face plane and for a
/// vertex on an edge the edge neighbours do, so neither passes.
fn extreme_vertices(hull: &ConvexHull3DVec) -> Vec<[f64; 3]> {
    let count = hull.vertices.len();
    let mut normal_sums = vec![Vec3::zeros(); count];
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); count];
    for face in &hull.faces {
        let [a, b, c] = *face;
        let (pa, pb, pc) = (hull.vertices[a], hull.vertices[b], hull.vertices[c]);
        let normal = (pb - pa).cross(&(pc - pa));
        let length = normal.norm();
        if length <= 0.0 {
            continue;
        }
        let unit = normal / length;
        for (vertex, others) in [(a, [b, c]), (b, [c, a]), (c, [a, b])] {
            normal_sums[vertex] += unit;
            neighbours[vertex].extend(others);
        }
    }

    let mut lowest = hull.vertices.first().copied().unwrap_or_else(Vec3::zeros);
    let mut highest = lowest;
    for v in &hull.vertices {
        lowest = lowest.inf(v);
        highest = highest.sup(v);
    }
    let tolerance = FLAT_VERTEX_TOLERANCE * (highest - lowest).norm();

    hull.vertices
        .iter()
        .enumerate()
        .filter(|&(index, vertex)| {
            let Some(direction) = normal_sums[index].try_normalize(1e-12) else {
                return true;
            };
            neighbours[index]
                .iter()
                .all(|&other| direction.dot(&(vertex - hull.vertices[other])) > tolerance)
        })
        .map(|(_, v)| [v.x, v.y, v.z])
        .collect()
}

/// Volume of the convex hull of `points` (zero when the hull is degenerate).
pub(super) fn hull_volume(points: &[[f64; 3]]) -> f64 {
    quickhull(points).map_or(0.0, |hull| hull.volume())
}

/// Approximate convex decomposition of a triangle mesh.
///
/// Repeatedly cuts the most concave part of the mesh with an axis-aligned plane
/// until every part is within `max_concavity` of its own convex hull or
/// `max_hulls` parts exist. The output is the convex hull of each part.
///
/// **Concavity.** The concavity of a part is the largest distance between the
/// surface of its convex hull and the surface of the part itself (hull to
/// mesh). Sample points on the faces of the hull are matched against the
/// nearest point of the part's triangles through a bounding-volume hierarchy.
/// A distance is used rather than "how deep do the vertices sit inside the
/// hull" because the latter is blind to many concave shapes: every vertex of an
/// extruded L-shape lies on the top or bottom face of its hull, so its depth
/// below the hull surface is zero although the shape is clearly concave.
///
/// **Splitting.** The part with the largest concavity above the threshold is
/// cut by an axis-aligned plane. Candidate planes are the bounding-box middle,
/// the quartiles, and the planes through the deepest point of the pocket that
/// produced the concavity; the candidate whose two halves have the smallest
/// total hull volume wins. Triangles that cross the plane are clipped, so each
/// half keeps exactly its own geometry. Hull faces that lie on a cut plane are
/// artificial caps and are ignored when measuring the concavity of the halves.
#[derive(Debug, Clone)]
pub struct ConvexDecomposition {
    /// Maximum number of convex hulls (at least one hull is always produced).
    pub max_hulls: usize,
    /// Maximum concavity: a part whose hull surface is further than this from
    /// the part's own surface (in mesh units) is split.
    pub max_concavity: f64,
    /// Sampling resolution: the number of subdivisions per hull-face edge used
    /// to probe the hull surface, clamped to `3..=8`. Higher values find
    /// narrower pockets at a higher cost.
    pub resolution: usize,
    /// Output convex hulls (vertex lists), one per part.
    pub hulls: Vec<Vec<[f64; 3]>>,
}

impl ConvexDecomposition {
    /// Create a new convex decomposition.
    pub fn new(max_hulls: usize, max_concavity: f64, resolution: usize) -> Self {
        Self {
            max_hulls,
            max_concavity,
            resolution,
            hulls: Vec::new(),
        }
    }

    /// Decompose a mesh into convex parts.
    ///
    /// The mesh surface is given by the triangles `indices` into `vertices`.
    /// Triangles with an out-of-range index or a non-finite corner are
    /// skipped, and vertices that no triangle references are ignored. The
    /// result replaces [`ConvexDecomposition::hulls`]; it is empty when no valid
    /// triangle exists.
    pub fn decompose(&mut self, vertices: &[[f64; 3]], indices: &[[usize; 3]]) {
        self.hulls.clear();
        let triangles = gather_triangles(vertices, indices);
        if triangles.is_empty() {
            return;
        }

        let bounds = triangle_bounds(&triangles);
        let scale = diagonal(&bounds);
        let context = Context {
            sampling: self.resolution.clamp(MIN_SAMPLING, MAX_SAMPLING),
            orientation: orientation_sign(&triangles),
            noise: RELATIVE_NOISE * scale,
        };
        let tolerance = self.max_concavity.max(context.noise);
        let max_parts = self.max_hulls.max(1);

        let mut parts = vec![Part::evaluate(triangles, Vec::new(), &context)];
        while parts.len() < max_parts {
            let target = parts
                .iter()
                .enumerate()
                .filter(|(_, part)| part.splittable && part.concavity > tolerance)
                .max_by(|(_, a), (_, b)| a.concavity.total_cmp(&b.concavity))
                .map(|(index, _)| index);
            let Some(target) = target else {
                break;
            };
            match parts[target].split(&context) {
                Some((lower, upper)) => {
                    parts[target] = lower;
                    parts.insert(target + 1, upper);
                }
                None => parts[target].splittable = false,
            }
        }

        self.hulls = parts.into_iter().map(|part| part.hull_points).collect();
    }

    /// Number of output hulls.
    pub fn num_hulls(&self) -> usize {
        self.hulls.len()
    }
}

/// Parameters shared by every step of one decomposition.
struct Context {
    /// Subdivisions per hull-face edge when probing the hull surface.
    sampling: usize,
    /// `+1` when the triangle winding of the mesh is counter-clockwise seen
    /// from outside, `-1` when it is clockwise.
    orientation: f64,
    /// Length below which two coordinates are considered equal.
    noise: f64,
}

/// A candidate cutting plane together with the geometry on both sides.
struct SplitCandidate {
    /// Total hull volume of the two halves (smaller is better).
    score: f64,
    axis: usize,
    position: f64,
    below: Vec<Triangle>,
    above: Vec<Triangle>,
}

/// The widest gap between a part and its convex hull.
#[derive(Debug, Clone, Copy)]
struct Pocket {
    /// Sample point on the hull surface.
    hull_point: [f64; 3],
    /// Nearest point of the part surface: the bottom of the pocket.
    mesh_point: [f64; 3],
}

/// A part of the mesh together with its convex hull and concavity.
struct Part {
    /// Surface triangles (clipped fragments of the input triangles).
    triangles: Vec<Triangle>,
    /// Planes `(axis, position)` that cut this part out of the mesh. Hull faces
    /// lying on them are artificial caps, not mesh surface.
    cuts: Vec<(usize, f64)>,
    /// Vertices of the convex hull (all distinct points when it is degenerate).
    hull_points: Vec<[f64; 3]>,
    /// Largest distance from the hull surface to the part surface.
    concavity: f64,
    /// Where the largest gap between hull and part is, if any.
    pocket: Option<Pocket>,
    /// `false` once no cutting plane could split the part.
    splittable: bool,
}

impl Part {
    /// Build a part and measure its concavity.
    fn evaluate(triangles: Vec<Triangle>, cuts: Vec<(usize, f64)>, context: &Context) -> Part {
        let corners: Vec<[f64; 3]> = triangles.iter().flatten().copied().collect();
        let points = sorted_unique_points(&corners);
        let (hull_points, concavity, pocket) = match quickhull(&points) {
            Some(hull) => {
                let gap = hull_gap(&hull, &triangles, &cuts, context);
                (extreme_vertices(&hull), gap.0, gap.1)
            }
            None => (points, 0.0, None),
        };
        Part {
            triangles,
            cuts,
            hull_points,
            concavity,
            pocket,
            splittable: true,
        }
    }

    /// Cut the part with the best axis-aligned plane. `None` when no plane
    /// leaves geometry on both sides.
    fn split(&self, context: &Context) -> Option<(Part, Part)> {
        let bounds = triangle_bounds(&self.triangles);
        let extents = [
            bounds.max[0] - bounds.min[0],
            bounds.max[1] - bounds.min[1],
            bounds.max[2] - bounds.min[2],
        ];
        let mut axes = [0usize, 1, 2];
        axes.sort_by(|&a, &b| extents[b].total_cmp(&extents[a]));

        let mut best: Option<SplitCandidate> = None;
        for &axis in &axes {
            let low = bounds.min[axis];
            let extent = extents[axis];
            if extent <= context.noise {
                continue;
            }
            let mut positions = vec![low + 0.5 * extent];
            if let Some(pocket) = self.pocket {
                positions.push(pocket.mesh_point[axis]);
                positions.push(pocket.hull_point[axis]);
            }
            positions.push(low + 0.25 * extent);
            positions.push(low + 0.75 * extent);

            let mut tried: Vec<f64> = Vec::new();
            for position in positions {
                let inside = position > low + CUT_MARGIN * extent
                    && position < bounds.max[axis] - CUT_MARGIN * extent;
                if !inside || tried.iter().any(|t| (t - position).abs() <= context.noise) {
                    continue;
                }
                tried.push(position);
                let (below, above) = partition_triangles(
                    &self.triangles,
                    axis,
                    position,
                    context.orientation,
                    context.noise,
                );
                if below.is_empty() || above.is_empty() {
                    continue;
                }
                let score = triangles_hull_volume(&below) + triangles_hull_volume(&above);
                let improves = best
                    .as_ref()
                    .is_none_or(|b| score < b.score - 1e-12 * (1.0 + b.score));
                if improves {
                    best = Some(SplitCandidate {
                        score,
                        axis,
                        position,
                        below,
                        above,
                    });
                }
            }
        }

        let chosen = best?;
        let mut cuts = self.cuts.clone();
        cuts.push((chosen.axis, chosen.position));
        Some((
            Part::evaluate(chosen.below, cuts.clone(), context),
            Part::evaluate(chosen.above, cuts, context),
        ))
    }
}

/// Collect the usable triangles of an indexed mesh.
fn gather_triangles(vertices: &[[f64; 3]], indices: &[[usize; 3]]) -> Vec<Triangle> {
    indices
        .iter()
        .filter_map(|tri| {
            let a = *vertices.get(tri[0])?;
            let b = *vertices.get(tri[1])?;
            let c = *vertices.get(tri[2])?;
            let finite = [a, b, c].iter().flatten().all(|v| v.is_finite());
            finite.then_some([a, b, c])
        })
        .collect()
}

/// Axis-aligned bounds of all triangle corners.
fn triangle_bounds(triangles: &[Triangle]) -> Aabb3 {
    let mut bounds = Aabb3::empty();
    for corner in triangles.iter().flatten() {
        bounds.include_point(*corner);
    }
    bounds
}

/// Length of the diagonal of `bounds`.
fn diagonal(bounds: &Aabb3) -> f64 {
    norm3(sub3(bounds.max, bounds.min))
}

/// Volume of the convex hull of the corners of `triangles`.
fn triangles_hull_volume(triangles: &[Triangle]) -> f64 {
    let corners: Vec<[f64; 3]> = triangles.iter().flatten().copied().collect();
    hull_volume(&sorted_unique_points(&corners))
}

/// `+1` when the mesh winding is counter-clockwise seen from outside (positive
/// signed volume), `-1` otherwise. Only meaningful for (nearly) closed meshes;
/// the outward normal of a triangle is `orientation * (b - a) x (c - a)`.
fn orientation_sign(triangles: &[Triangle]) -> f64 {
    let count = (triangles.len() * 3) as f64;
    let mut centroid = [0.0_f64; 3];
    for corner in triangles.iter().flatten() {
        centroid = add3(centroid, *corner);
    }
    let centroid = scale3(centroid, 1.0 / count);
    let six_volume: f64 = triangles
        .iter()
        .map(|[a, b, c]| {
            let (a, b, c) = (sub3(*a, centroid), sub3(*b, centroid), sub3(*c, centroid));
            dot3(a, cross3(b, c))
        })
        .sum();
    if six_volume < 0.0 { -1.0 } else { 1.0 }
}

/// Split triangles by the plane `x[axis] = position`.
///
/// Triangles that cross the plane are clipped. A triangle lying in the plane
/// goes to the side its outward normal points away from (the side of the
/// solid).
fn partition_triangles(
    triangles: &[Triangle],
    axis: usize,
    position: f64,
    orientation: f64,
    noise: f64,
) -> (Vec<Triangle>, Vec<Triangle>) {
    let mut plane_point = [0.0_f64; 3];
    plane_point[axis] = position;
    let mut up = [0.0_f64; 3];
    up[axis] = 1.0;
    let down = scale3(up, -1.0);

    let mut below = Vec::new();
    let mut above = Vec::new();
    for tri in triangles {
        let heights = [
            tri[0][axis] - position,
            tri[1][axis] - position,
            tri[2][axis] - position,
        ];
        let highest = heights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let lowest = heights.iter().copied().fold(f64::INFINITY, f64::min);
        if highest <= noise && lowest >= -noise {
            let normal = scale3(
                cross3(sub3(tri[1], tri[0]), sub3(tri[2], tri[0])),
                orientation,
            );
            if normal[axis] > 0.0 {
                below.push(*tri);
            } else {
                above.push(*tri);
            }
        } else if highest <= noise {
            below.push(*tri);
        } else if lowest >= -noise {
            above.push(*tri);
        } else {
            let polygon = [tri[0], tri[1], tri[2]];
            fan_triangulate(
                &clip_polygon_by_plane(&polygon, plane_point, down),
                noise,
                &mut below,
            );
            fan_triangulate(
                &clip_polygon_by_plane(&polygon, plane_point, up),
                noise,
                &mut above,
            );
        }
    }
    (below, above)
}

/// Append the fan triangulation of a convex polygon to `out`.
///
/// Neighbouring corners closer than `noise` are merged first: clipping a
/// triangle through one of its corners yields such duplicates, and a fan over
/// them would emit zero-area triangles.
fn fan_triangulate(polygon: &[[f64; 3]], noise: f64, out: &mut Vec<Triangle>) {
    let mut corners: Vec<[f64; 3]> = Vec::with_capacity(polygon.len());
    for &corner in polygon {
        if corners
            .last()
            .is_none_or(|&last| norm3(sub3(corner, last)) > noise)
        {
            corners.push(corner);
        }
    }
    while corners.len() > 1
        && corners
            .first()
            .zip(corners.last())
            .is_some_and(|(&first, &last)| norm3(sub3(first, last)) <= noise)
    {
        corners.pop();
    }
    for pair in corners.windows(2).skip(1) {
        out.push([corners[0], pair[0], pair[1]]);
    }
}

/// Largest distance between the hull surface and the part surface.
///
/// Returns the distance and the [`Pocket`] that realises it. Samples on a cut
/// plane belong to an artificial cap and are skipped.
fn hull_gap(
    hull: &ConvexHull3DVec,
    triangles: &[Triangle],
    cuts: &[(usize, f64)],
    context: &Context,
) -> (f64, Option<Pocket>) {
    let index = TriangleIndex::new(triangles);
    let mut worst = 0.0_f64;
    let mut pocket = None;
    for face in &hull.faces {
        let corner = |k: usize| {
            let v = hull.vertices[face[k]];
            [v.x, v.y, v.z]
        };
        let corners = [corner(0), corner(1), corner(2)];
        for sample in lattice_samples(&corners, context.sampling) {
            if cuts
                .iter()
                .any(|&(axis, position)| (sample[axis] - position).abs() <= context.noise)
            {
                continue;
            }
            if let Some((distance, nearest)) = index.nearest(sample)
                && distance > worst
            {
                worst = distance;
                pocket = Some(Pocket {
                    hull_point: sample,
                    mesh_point: nearest,
                });
            }
        }
    }
    (worst, pocket)
}

/// Barycentric lattice with `n` subdivisions per edge over a triangle, without
/// the three corner points.
fn lattice_samples(corners: &[[f64; 3]; 3], n: usize) -> Vec<[f64; 3]> {
    let mut samples = Vec::with_capacity((n + 1) * (n + 2) / 2);
    let step = 1.0 / n.max(1) as f64;
    for i in 0..=n {
        for j in 0..=(n - i) {
            let k = n - i - j;
            if i == n || j == n || k == n {
                continue;
            }
            samples.push(add3(
                add3(
                    scale3(corners[0], k as f64 * step),
                    scale3(corners[1], i as f64 * step),
                ),
                scale3(corners[2], j as f64 * step),
            ));
        }
    }
    samples
}

/// Nearest-triangle queries over a set of triangles, accelerated by a BVH over
/// the triangle bounding boxes.
struct TriangleIndex<'a> {
    triangles: &'a [Triangle],
    bvh: BvhBroadphase,
    /// Radius of the first search box; doubled until the nearest triangle is
    /// certain.
    first_radius: f64,
    /// Search radius at which every triangle is certainly inside the box.
    covering_radius: f64,
}

impl<'a> TriangleIndex<'a> {
    fn new(triangles: &'a [Triangle]) -> Self {
        let boxes: Vec<Aabb> = triangles
            .iter()
            .map(|tri| {
                let mut bounds = Aabb3::empty();
                for corner in tri {
                    bounds.include_point(*corner);
                }
                Aabb::new(
                    Vec3::new(bounds.min[0], bounds.min[1], bounds.min[2]),
                    Vec3::new(bounds.max[0], bounds.max[1], bounds.max[2]),
                )
            })
            .collect();
        let mut bvh = BvhBroadphase::new();
        bvh.build(&boxes);
        let diameter = diagonal(&triangle_bounds(triangles));
        Self {
            triangles,
            bvh,
            first_radius: (diameter / 64.0).max(1e-12),
            covering_radius: 2.0 * diameter + 1e-9,
        }
    }

    /// Distance from `point` to the nearest triangle and the nearest point.
    fn nearest(&self, point: [f64; 3]) -> Option<(f64, [f64; 3])> {
        if self.triangles.is_empty() {
            return None;
        }
        let mut radius = self.first_radius;
        loop {
            let probe = Aabb::new(
                Vec3::new(point[0] - radius, point[1] - radius, point[2] - radius),
                Vec3::new(point[0] + radius, point[1] + radius, point[2] + radius),
            );
            let mut best: Option<(f64, [f64; 3])> = None;
            for candidate in self.bvh.query(&probe) {
                let Some([a, b, c]) = self.triangles.get(candidate) else {
                    continue;
                };
                let closest = MeshShapeQuery::closest_point_triangle(point, *a, *b, *c);
                let distance = norm3(sub3(point, closest));
                // A degenerate (zero-area) triangle has no defined closest
                // point; it must not poison the comparison.
                if !distance.is_finite() {
                    continue;
                }
                if best.is_none_or(|(d, _)| distance < d) {
                    best = Some((distance, closest));
                }
            }
            // Every triangle within `radius` of the point touches the probe
            // box, so a hit that close is the true nearest one.
            if best.is_some_and(|(d, _)| d <= radius) || radius >= self.covering_radius {
                return best;
            }
            radius *= 2.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compound_shapes::test_support::{
        Mesh, cube_mesh, extrude_polygon, flipped, icosphere, signed_volume, torus,
    };
    use oxiphysics_geometry::quickhull::is_point_set_convex;

    fn l_prism() -> Mesh {
        extrude_polygon(
            &[
                [0.0, 0.0],
                [2.0, 0.0],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ],
            1.0,
        )
    }

    fn u_prism() -> Mesh {
        extrude_polygon(
            &[
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 2.0],
                [2.0, 2.0],
                [2.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [0.0, 2.0],
            ],
            1.0,
        )
    }

    fn staircase_prism() -> Mesh {
        extrude_polygon(
            &[
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 1.0],
                [2.0, 1.0],
                [2.0, 2.0],
                [1.0, 2.0],
                [1.0, 3.0],
                [0.0, 3.0],
            ],
            1.0,
        )
    }

    fn c_prism() -> Mesh {
        extrude_polygon(
            &[
                [0.0, 0.0],
                [3.0, 0.0],
                [3.0, 1.0],
                [1.0, 1.0],
                [1.0, 2.0],
                [3.0, 2.0],
                [3.0, 3.0],
                [0.0, 3.0],
            ],
            1.0,
        )
    }

    fn decompose(mesh: &Mesh, max_hulls: usize, max_concavity: f64) -> ConvexDecomposition {
        let mut decomposition = ConvexDecomposition::new(max_hulls, max_concavity, 4);
        decomposition.decompose(&mesh.0, &mesh.1);
        decomposition
    }

    fn total_volume(decomposition: &ConvexDecomposition) -> f64 {
        decomposition.hulls.iter().map(|h| hull_volume(h)).sum()
    }

    fn as_cloud(points: &[[f64; 3]]) -> Vec<Vec3> {
        points.iter().map(|p| Vec3::new(p[0], p[1], p[2])).collect()
    }

    /// Every vertex of `mesh` lies inside (or on) one of the hulls.
    fn hulls_cover_vertices(decomposition: &ConvexDecomposition, mesh: &Mesh) -> bool {
        let hulls: Vec<ConvexHull3DVec> = decomposition
            .hulls
            .iter()
            .filter_map(|h| quickhull(h))
            .collect();
        mesh.0.iter().all(|v| {
            let p = Vec3::new(v[0], v[1], v[2]);
            hulls.iter().any(|hull| hull.contains_point(p))
        })
    }

    fn merge_meshes(a: &Mesh, b: &Mesh) -> Mesh {
        let offset = a.0.len();
        let mut vertices = a.0.clone();
        vertices.extend_from_slice(&b.0);
        let mut indices = a.1.clone();
        indices.extend(
            b.1.iter()
                .map(|t| [t[0] + offset, t[1] + offset, t[2] + offset]),
        );
        (vertices, indices)
    }

    /// The same surface as `mesh`, as a triangle soup: every triangle owns its
    /// three vertices instead of sharing them with its neighbours.
    fn soup(mesh: &Mesh) -> Mesh {
        let mut vertices = Vec::with_capacity(mesh.1.len() * 3);
        let mut indices = Vec::with_capacity(mesh.1.len());
        for triangle in &mesh.1 {
            let base = vertices.len();
            vertices.extend(triangle.iter().map(|&corner| mesh.0[corner]));
            indices.push([base, base + 1, base + 2]);
        }
        (vertices, indices)
    }

    #[test]
    fn triangle_soup_decomposes_like_the_indexed_mesh() {
        let indexed = l_prism();
        let soup_mesh = soup(&indexed);
        assert_eq!(soup_mesh.0.len(), 3 * indexed.1.len());
        let from_soup = decompose(&soup_mesh, 8, 0.05);
        let from_indexed = decompose(&indexed, 8, 0.05);
        assert_eq!(from_soup.num_hulls(), 2);
        let difference = (total_volume(&from_soup) - total_volume(&from_indexed)).abs();
        assert!(difference < 1e-9, "volumes differ by {difference}");
        assert!(hulls_cover_vertices(&from_soup, &indexed));
    }

    #[test]
    fn convex_cube_stays_a_single_hull() {
        let mesh = cube_mesh(1.0, [0.0; 3]);
        let decomposition = decompose(&mesh, 8, 0.01);
        assert_eq!(decomposition.num_hulls(), 1);
        assert_eq!(decomposition.hulls[0].len(), 8);
        assert!((total_volume(&decomposition) - 8.0).abs() < 1e-9);
    }

    #[test]
    fn convex_icosphere_stays_a_single_hull() {
        let mesh = icosphere(2.0);
        let decomposition = decompose(&mesh, 8, 0.01);
        assert_eq!(decomposition.num_hulls(), 1);
        assert_eq!(decomposition.hulls[0].len(), 42);
    }

    #[test]
    fn l_shaped_prism_is_split_into_convex_parts() {
        // Every vertex of an extruded L lies on the top or bottom face of its
        // hull, so a "depth of the vertices below the hull" measure sees no
        // concavity at all; the hull-to-surface distance does.
        let mesh = l_prism();
        let decomposition = decompose(&mesh, 8, 0.05);
        assert_eq!(decomposition.num_hulls(), 2, "{:?}", decomposition.hulls);
        for hull in &decomposition.hulls {
            assert!(is_point_set_convex(&as_cloud(hull)));
            assert!(hull_volume(hull) > 0.5);
        }
        assert!(hulls_cover_vertices(&decomposition, &mesh));
        // Cutting at the notch yields the two boxes of the L: volume 2 + 1.
        let volume = total_volume(&decomposition);
        assert!((volume - 3.0).abs() < 1e-6, "total volume {volume}");
        // The hull of the whole L has volume 3.5.
        assert!(volume < hull_volume(&mesh.0) - 0.4);
    }

    #[test]
    fn mesh_winding_does_not_matter() {
        let mesh = flipped(&l_prism());
        let decomposition = decompose(&mesh, 8, 0.05);
        assert_eq!(decomposition.num_hulls(), 2);
        assert!((total_volume(&decomposition) - 3.0).abs() < 1e-6);
        assert!(hulls_cover_vertices(&decomposition, &mesh));
    }

    #[test]
    fn separate_components_get_their_own_hull() {
        let mesh = merge_meshes(&cube_mesh(1.0, [0.0; 3]), &cube_mesh(1.0, [5.0, 0.0, 0.0]));
        let decomposition = decompose(&mesh, 8, 0.05);
        assert_eq!(decomposition.num_hulls(), 2);
        for hull in &decomposition.hulls {
            assert_eq!(hull.len(), 8);
            assert!((hull_volume(hull) - 8.0).abs() < 1e-9);
        }
    }

    #[test]
    fn u_shaped_prism_needs_three_convex_parts() {
        let mesh = u_prism();
        let decomposition = decompose(&mesh, 16, 0.05);
        assert!(
            decomposition.num_hulls() >= 3,
            "{}",
            decomposition.num_hulls()
        );
        for hull in &decomposition.hulls {
            assert!(is_point_set_convex(&as_cloud(hull)));
        }
        assert!(hulls_cover_vertices(&decomposition, &mesh));
        let volume = total_volume(&decomposition);
        // The U has volume 5, its hull 6.
        assert!(volume >= 5.0 - 1e-6, "volume {volume}");
        assert!(volume < 5.6, "volume {volume}");
    }

    #[test]
    fn staircase_prism_is_cut_into_boxes_through_its_corners() {
        // Every cut plane passes through vertices of the staircase, so the
        // clipped triangles contain duplicate corners; they must neither
        // produce NaN distances nor hide the remaining concavity.
        let mesh = staircase_prism();
        let decomposition = decompose(&mesh, 16, 0.05);
        assert!(
            decomposition.num_hulls() >= 3,
            "{}",
            decomposition.num_hulls()
        );
        for hull in &decomposition.hulls {
            assert!(is_point_set_convex(&as_cloud(hull)));
        }
        assert!(hulls_cover_vertices(&decomposition, &mesh));
        // The staircase has volume 6, its hull 6.5 + 1 = 7.5.
        let volume = total_volume(&decomposition);
        assert!(volume >= 6.0 - 1e-6, "volume {volume}");
        assert!(volume < 6.3, "volume {volume}");
    }

    #[test]
    fn c_shaped_prism_is_split_around_its_pocket() {
        let mesh = c_prism();
        let decomposition = decompose(&mesh, 16, 0.05);
        assert!(
            decomposition.num_hulls() >= 3,
            "{}",
            decomposition.num_hulls()
        );
        for hull in &decomposition.hulls {
            assert!(is_point_set_convex(&as_cloud(hull)));
        }
        assert!(hulls_cover_vertices(&decomposition, &mesh));
        // The C has volume 7, its hull 9.
        let volume = total_volume(&decomposition);
        assert!(volume >= 7.0 - 1e-6, "volume {volume}");
        assert!(volume < 7.5, "volume {volume}");
    }

    #[test]
    fn hull_vertices_are_exactly_the_corners() {
        // Corners of [-1, 1]^3 plus an interior point, a face centre, a point
        // inside a face and an edge midpoint. The face centre comes last in the
        // x direction, which makes QuickHull seed with it.
        let (mut points, _) = cube_mesh(1.0, [0.0; 3]);
        let corners = points.clone();
        points.push([0.0, 0.0, 0.0]);
        points.push([1.0, 0.0, 0.0]);
        points.push([0.3, -0.2, 1.0]);
        points.push([1.0, 1.0, 0.0]);
        let mut reversed = points.clone();
        reversed.reverse();
        let mut by_x = sorted_unique_points(&points);
        by_x.reverse();
        for order in [points, reversed, by_x] {
            let found = hull_vertices(&order);
            assert_eq!(found.len(), 8, "{found:?}");
            for corner in &corners {
                assert!(found.contains(corner), "corner {corner:?} missing");
            }
        }
    }

    #[test]
    fn fan_triangulation_merges_duplicate_corners() {
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [1.0, 1.0, 0.0];
        let d = [0.0, 1.0, 0.0];
        let mut out = Vec::new();
        // Duplicate in the middle, at the end (wrap-around) and nearly equal.
        fan_triangulate(
            &[a, b, [1.0 + 1e-13, 0.0, 0.0], c, d, [1e-14, 0.0, 0.0]],
            1e-9,
            &mut out,
        );
        assert_eq!(out.len(), 2, "{out:?}");
        for [p, q, r] in &out {
            let area = norm3(cross3(sub3(*q, *p), sub3(*r, *p))) * 0.5;
            assert!(area > 0.4, "degenerate fragment {p:?} {q:?} {r:?}");
        }
        // A polygon that collapses to a point or a segment gives no triangles.
        let mut out = Vec::new();
        fan_triangulate(&[a, a, a], 1e-9, &mut out);
        fan_triangulate(&[a, b, a], 1e-9, &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn nearest_triangle_ignores_degenerate_triangles() {
        // Two coincident corners make the closest-point formula divide 0 by 0.
        let degenerate = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        assert!(
            MeshShapeQuery::closest_point_triangle(
                [0.25, 0.25, 3.0],
                degenerate[0],
                degenerate[1],
                degenerate[2]
            )
            .iter()
            .any(|c| c.is_nan()),
            "the scene must really produce a NaN closest point"
        );
        let flat_line = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let real = [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]];
        let triangles = [degenerate, flat_line, real];
        let index = TriangleIndex::new(&triangles);
        let (distance, point) = index
            .nearest([0.25, 0.25, 3.0])
            .expect("the real triangle is found");
        assert!((distance - 2.0).abs() < 1e-12, "distance {distance}");
        assert!((point[2] - 1.0).abs() < 1e-12);
        // Only a degenerate triangle: no answer instead of NaN.
        let only_bad = [degenerate];
        assert!(
            TriangleIndex::new(&only_bad)
                .nearest([0.5, 0.5, 0.5])
                .is_none()
        );
    }

    #[test]
    fn max_hulls_caps_the_output() {
        let mesh = u_prism();
        assert_eq!(decompose(&mesh, 2, 0.01).num_hulls(), 2);
        assert_eq!(decompose(&mesh, 1, 0.01).num_hulls(), 1);
        assert_eq!(decompose(&mesh, 0, 0.01).num_hulls(), 1);
    }

    #[test]
    fn concavity_threshold_controls_splitting() {
        let mesh = l_prism();
        assert_eq!(decompose(&mesh, 8, 5.0).num_hulls(), 1);
        assert_eq!(decompose(&mesh, 8, 0.05).num_hulls(), 2);
        // A threshold of zero still terminates, and a convex part is never cut.
        assert_eq!(decompose(&mesh, 8, 0.0).num_hulls(), 2);
        assert_eq!(decompose(&cube_mesh(1.0, [0.0; 3]), 8, 0.0).num_hulls(), 1);
    }

    #[test]
    fn sampling_resolution_is_clamped() {
        let mesh = l_prism();
        for resolution in [0, 1, 3, 8, 1000] {
            let mut decomposition = ConvexDecomposition::new(8, 0.05, resolution);
            decomposition.decompose(&mesh.0, &mesh.1);
            assert_eq!(decomposition.num_hulls(), 2, "resolution {resolution}");
        }
    }

    #[test]
    fn invalid_input_is_skipped() {
        let (mut vertices, mut indices) = cube_mesh(1.0, [0.0; 3]);
        // A stray unreferenced vertex, a NaN vertex used by one triangle and a
        // triangle with an out-of-range index.
        vertices.push([100.0, 100.0, 100.0]);
        vertices.push([f64::NAN, 0.0, 0.0]);
        indices.push([0, 1, 9]);
        indices.push([0, 1, 1000]);
        let mut decomposition = ConvexDecomposition::new(8, 0.01, 4);
        decomposition.decompose(&vertices, &indices);
        assert_eq!(decomposition.num_hulls(), 1);
        assert_eq!(decomposition.hulls[0].len(), 8);
        assert!((hull_volume(&decomposition.hulls[0]) - 8.0).abs() < 1e-9);
    }

    #[test]
    fn empty_input_gives_no_hulls() {
        let mut decomposition = ConvexDecomposition::new(8, 0.01, 4);
        decomposition.decompose(&[], &[]);
        assert_eq!(decomposition.num_hulls(), 0);
        decomposition.decompose(&[[0.0; 3], [1.0, 0.0, 0.0]], &[]);
        assert_eq!(decomposition.num_hulls(), 0);
    }

    #[test]
    fn flat_mesh_yields_its_own_points() {
        let vertices = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ];
        let indices = vec![[0, 1, 2], [0, 2, 3]];
        let mut decomposition = ConvexDecomposition::new(8, 0.01, 4);
        decomposition.decompose(&vertices, &indices);
        assert_eq!(decomposition.num_hulls(), 1);
        assert_eq!(decomposition.hulls[0].len(), 4);
    }

    #[test]
    fn torus_is_cut_into_arcs_that_hug_the_surface() {
        // A curved, genuinely non-convex mesh of 576 triangles whose cut planes
        // clip many triangles. The hole is the only concavity: the hull of the
        // whole torus (volume ~87.7) is a poor fit for the solid (volume ~55.9).
        let mesh = torus(3.0, 1.0, 24, 12);
        let mesh_volume = signed_volume(&mesh);
        let whole = hull_volume(&mesh.0);

        // Four parts: the 4-fold symmetry of the mesh gives four equal quarters.
        let coarse = decompose(&mesh, 4, 0.1);
        assert_eq!(coarse.num_hulls(), 4);
        let quarter = hull_volume(&coarse.hulls[0]);
        for hull in &coarse.hulls {
            assert!((hull_volume(hull) - quarter).abs() < 1e-4 * quarter);
        }

        // With a tight threshold the splitting stops by itself, below the cap.
        let fine = decompose(&mesh, 32, 0.05);
        assert!(
            (16..32).contains(&fine.num_hulls()),
            "{} hulls",
            fine.num_hulls()
        );

        for decomposition in [&coarse, &fine] {
            for hull in &decomposition.hulls {
                assert!(is_point_set_convex(&as_cloud(hull)));
            }
            assert!(hulls_cover_vertices(decomposition, &mesh));
            // The hulls enclose the whole solid ...
            assert!(total_volume(decomposition) >= mesh_volume - 1e-6);
        }
        // ... and the more parts there are, the tighter they fit.
        assert!(total_volume(&coarse) < 0.8 * whole);
        assert!(total_volume(&fine) < total_volume(&coarse));
        assert!(total_volume(&fine) < 1.05 * mesh_volume);
    }

    #[test]
    fn decomposition_is_deterministic_and_replaces_previous_output() {
        let mesh = u_prism();
        let mut decomposition = ConvexDecomposition::new(16, 0.05, 4);
        decomposition.decompose(&mesh.0, &mesh.1);
        let first = decomposition.hulls.clone();
        decomposition.decompose(&mesh.0, &mesh.1);
        assert_eq!(decomposition.hulls, first);
        decomposition.decompose(&cube_mesh(1.0, [0.0; 3]).0, &cube_mesh(1.0, [0.0; 3]).1);
        assert_eq!(decomposition.num_hulls(), 1);
    }
}
