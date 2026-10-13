// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Mesh-driven approximate convex decomposition (ACD).
//!
//! The input is a closed triangle mesh that bounds a solid `S`.  A *cluster*
//! is the part `S ∩ H` of the solid cut out by the half-spaces `H` chosen by
//! earlier splits; it is represented by the mesh triangles clipped to `H`.
//! The convex hull of the clipped vertices is the convex part the cluster
//! would be replaced by (the hull of `S ∩ H` is spanned by the vertices of its
//! boundary, and the boundary of `S ∩ H` consists of clipped mesh polygons and
//! planar cut caps whose corners are mesh-edge/plane intersections that the
//! clipped polygons already contain).
//!
//! # Concavity
//!
//! The concavity of a cluster is how far its convex hull sticks out beyond the
//! real part: the largest distance from a point of the hull surface to the
//! part (the one-sided Hausdorff distance from the hull to the part).  It is
//! sampled on a barycentric grid over every hull facet.  A probe point that
//! lies on the mesh surface or inside the solid `S` is at distance zero: a
//! point of the hull boundary inside `S` belongs to `S ∩ H` and therefore lies
//! on the real boundary of the part, namely on a planar cap created by an
//! earlier cut.  This is what keeps cut faces from being mistaken for concave
//! regions.  A probe in empty space is at the distance of the nearest point of
//! the part's clipped mesh polygons: the segment from the probe to the nearest
//! point of the part first enters the solid on the surface, inside `H` (which
//! is convex), so the nearest point of the part is always a surface point.
//!
//! # Splitting
//!
//! Clusters are split best-first (largest concavity first) until every
//! cluster is below `max_concavity` or `max_parts` is reached.  The split
//! planes tried pass through the *most concave point* (the surface point
//! nearest to the worst probe) and through the cluster centroid, with the
//! coordinate axes and the principal axes of the cluster as normals, plus the
//! planes of the mesh face nearest to the worst probe and of the faces sharing
//! a corner with it (the extensions of the walls of a notch).  The plane for
//! which the concavities of the two children add up to the least wins (as in
//! V-HACD; ties go to the smaller total hull volume).

use oxiphysics_core::math::Mat3;

use super::{AcdConfig, ConvexPart, ConvexParts, Vec3, add, cross, dot, neg, norm, scale, sub};
use crate::convex_hull::ConvexHull3D;
use crate::{Error, Result};

/// Each hull facet is probed on the barycentric grid `(i, j, k) / N` with
/// `i + j + k = N` and `N` this value.
const FACET_SAMPLE_RESOLUTION: usize = 4;
/// Distances below this fraction of the mesh diagonal count as zero.
const RELATIVE_TOLERANCE: f64 = 1e-9;
/// Children whose hull volume is not at least this fraction smaller than the
/// parent's do not make progress.
const MIN_RELATIVE_VOLUME_DROP: f64 = 1e-6;
/// Offset, as a fraction of the mesh diagonal, of the probes that decide on
/// which side of a coplanar face the solid lies.
const SIDE_PROBE_FRACTION: f64 = 1e-6;
/// Directions of the rays used for the inside/outside parity test.  They are
/// deliberately not aligned with the axes of typical meshes.
const PARITY_DIRECTIONS: [Vec3; 3] = [
    [0.572_357_329_1, 0.405_431_870_3, 0.711_926_105_2],
    [-0.318_190_447_7, 0.883_277_219_4, -0.343_201_558_3],
    [0.241_503_982_6, -0.512_338_196_4, -0.823_947_025_1],
];

// ─── Small geometry helpers ───────────────────────────────────────────────

fn unit(v: Vec3) -> Vec3 {
    let n = norm(v);
    if n > 0.0 { scale(v, 1.0 / n) } else { v }
}

/// Closest point of triangle `abc` to `p` (Ericson, *Real-Time Collision
/// Detection*, 5.1.5).
pub(super) fn closest_point_on_triangle(p: Vec3, a: Vec3, b: Vec3, c: Vec3) -> Vec3 {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return a;
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return b;
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        return add(a, scale(ab, d1 / (d1 - d3)));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return c;
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        return add(a, scale(ac, d2 / (d2 - d6)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return add(b, scale(sub(c, b), w));
    }
    let denom = 1.0 / (va + vb + vc);
    add(add(a, scale(ab, vb * denom)), scale(ac, vc * denom))
}

/// Whether the ray `origin + t · direction` crosses the interior of the
/// triangle at some `t > 0` (Möller–Trumbore, two-sided, strict interior).
fn ray_crosses_triangle(origin: Vec3, direction: Vec3, tri: &[Vec3; 3]) -> bool {
    let e1 = sub(tri[1], tri[0]);
    let e2 = sub(tri[2], tri[0]);
    let p = cross(direction, e2);
    let det = dot(e1, p);
    if det.abs() <= 1e-14 * norm(e1) * norm(e2) {
        return false;
    }
    let inv = 1.0 / det;
    let t_vec = sub(origin, tri[0]);
    let u = dot(t_vec, p) * inv;
    if u <= 0.0 || u >= 1.0 {
        return false;
    }
    let q = cross(t_vec, e1);
    let v = dot(direction, q) * inv;
    v > 0.0 && u + v < 1.0 && dot(e2, q) * inv > 0.0
}

// ─── Reference surface ────────────────────────────────────────────────────

/// The input mesh: triangle corners and the bounding-box diagonal.
struct Surface {
    triangles: Vec<[Vec3; 3]>,
    diagonal: f64,
}

impl Surface {
    /// Validate the mesh and resolve the triangle corners; triangles with
    /// zero area are dropped.
    fn new(vertices: &[Vec3], triangles: &[[usize; 3]]) -> Result<Self> {
        if vertices.is_empty() || triangles.is_empty() {
            return Err(Error::invalid_mesh(
                "convex decomposition needs at least one vertex and one triangle",
            ));
        }
        if let Some(bad) = vertices
            .iter()
            .position(|v| v.iter().any(|c| !c.is_finite()))
        {
            return Err(Error::invalid_mesh(format!(
                "vertex {bad} has a non-finite coordinate"
            )));
        }
        let mut resolved = Vec::with_capacity(triangles.len());
        for tri in triangles {
            for &index in tri {
                if index >= vertices.len() {
                    return Err(Error::index_out_of_bounds(index, vertices.len()));
                }
            }
            let corners = [vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]];
            let area2 = norm(cross(
                sub(corners[1], corners[0]),
                sub(corners[2], corners[0]),
            ));
            if area2 > 0.0 {
                resolved.push(corners);
            }
        }
        if resolved.is_empty() {
            return Err(Error::degenerate_geometry(
                "every triangle of the mesh has zero area",
            ));
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for corner in resolved.iter().flatten() {
            for ((low, high), &c) in lo.iter_mut().zip(hi.iter_mut()).zip(corner.iter()) {
                *low = low.min(c);
                *high = high.max(c);
            }
        }
        let diagonal = norm(sub(hi, lo));
        if !(diagonal > 0.0 && diagonal.is_finite()) {
            return Err(Error::degenerate_geometry("mesh has no extent"));
        }
        Ok(Self {
            triangles: resolved,
            diagonal,
        })
    }

    fn tolerance(&self) -> f64 {
        RELATIVE_TOLERANCE * self.diagonal
    }

    /// Whether `p` lies inside the solid bounded by the mesh (ray-crossing
    /// parity, majority vote over three oblique rays).
    fn contains(&self, p: Vec3) -> bool {
        let votes = PARITY_DIRECTIONS
            .iter()
            .filter(|&&direction| {
                let crossings = self
                    .triangles
                    .iter()
                    .filter(|t| ray_crosses_triangle(p, direction, t))
                    .count();
                crossings % 2 == 1
            })
            .count();
        votes >= 2
    }
}

// ─── Convex hull of a cluster ─────────────────────────────────────────────

/// Hull of a cluster: vertices, outward triangles and mass properties.
struct HullData {
    vertices: Vec<Vec3>,
    faces: Vec<[usize; 3]>,
    volume: f64,
    centroid: Vec3,
}

impl HullData {
    /// Convex hull of `points`, `None` if they span no volume.
    fn build(points: &[Vec3]) -> Option<Self> {
        let hull = ConvexHull3D::build(points)?;
        let count = hull.vertices.len() as f64;
        let interior = scale(
            hull.vertices.iter().fold([0.0; 3], |acc, &v| add(acc, v)),
            1.0 / count,
        );
        // Outward triangles (flip any that face the interior), skipping
        // zero-area ones.
        let mut faces = Vec::with_capacity(hull.faces.len());
        for face in &hull.faces {
            let (a, b, c) = (
                hull.vertices[face[0]],
                hull.vertices[face[1]],
                hull.vertices[face[2]],
            );
            let normal = cross(sub(b, a), sub(c, a));
            let length = norm(normal);
            if !(length > 0.0 && length.is_finite()) {
                continue;
            }
            let mut tri = *face;
            if dot(normal, sub(interior, a)) > 0.0 {
                tri.swap(1, 2);
            }
            faces.push(tri);
        }
        let mut six_volume = 0.0;
        let mut weighted = [0.0; 3];
        for face in &faces {
            let a = sub(hull.vertices[face[0]], interior);
            let b = sub(hull.vertices[face[1]], interior);
            let c = sub(hull.vertices[face[2]], interior);
            let det = dot(a, cross(b, c));
            six_volume += det;
            weighted = add(weighted, scale(add(add(a, b), c), det / 4.0));
        }
        let volume = six_volume / 6.0;
        if !(volume.is_finite() && volume > 0.0) {
            return None;
        }
        let centroid = add(interior, scale(weighted, 1.0 / six_volume));
        Some(Self {
            vertices: hull.vertices,
            faces,
            volume,
            centroid,
        })
    }
}

// ─── Clusters ─────────────────────────────────────────────────────────────

/// Closest point of a cluster's clipped polygons to a query point.
struct Nearest {
    distance: f64,
    point: Vec3,
    polygon: usize,
}

/// The polygon (as a triangle fan) point of `polygons` closest to `p`.
fn nearest_on_polygons(polygons: &[Vec<Vec3>], p: Vec3) -> Nearest {
    let mut best = Nearest {
        distance: f64::INFINITY,
        point: p,
        polygon: 0,
    };
    for (index, polygon) in polygons.iter().enumerate() {
        for k in 1..polygon.len().saturating_sub(1) {
            let q = closest_point_on_triangle(p, polygon[0], polygon[k], polygon[k + 1]);
            let distance = norm(sub(p, q));
            if distance < best.distance {
                best = Nearest {
                    distance,
                    point: q,
                    polygon: index,
                };
            }
        }
    }
    best
}

/// Unit normal (Newell's method) and a point of a planar polygon.
fn polygon_plane(polygon: &[Vec3]) -> Option<(Vec3, Vec3)> {
    let mut normal = [0.0_f64; 3];
    for (i, a) in polygon.iter().enumerate() {
        let b = polygon[(i + 1) % polygon.len()];
        normal[0] += (a[1] - b[1]) * (a[2] + b[2]);
        normal[1] += (a[2] - b[2]) * (a[0] + b[0]);
        normal[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    let length = norm(normal);
    (length > 0.0 && length.is_finite()).then(|| (scale(normal, 1.0 / length), polygon[0]))
}

/// Remove points closer than `tol` to an earlier point (sweep over x).
fn dedupe_points(points: &[Vec3], tol: f64) -> Vec<Vec3> {
    let mut sorted: Vec<Vec3> = points.to_vec();
    sorted.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut kept: Vec<Vec3> = Vec::with_capacity(sorted.len());
    for p in sorted {
        let duplicate = kept
            .iter()
            .rev()
            .take_while(|q| p[0] - q[0] <= tol)
            .any(|q| norm(sub(p, *q)) <= tol);
        if !duplicate {
            kept.push(p);
        }
    }
    kept
}

/// The worst probe of a cluster: where the hull sticks out furthest beyond
/// the real part.
#[derive(Clone, Copy)]
struct Witness {
    /// The surface point nearest to the worst probe: the most concave point.
    nearest: Vec3,
    /// Index of the cluster polygon that point lies on.
    polygon: usize,
}

struct Cluster {
    /// Mesh triangles clipped to the cluster's half-spaces (convex polygons).
    polygons: Vec<Vec<Vec3>>,
    hull: HullData,
    concavity: f64,
    witness: Option<Witness>,
    /// Set when no split could improve the cluster.
    frozen: bool,
}

impl Cluster {
    /// Build a cluster from clipped polygons; `None` if its vertices span no
    /// volume.
    fn new(polygons: Vec<Vec<Vec3>>, surface: &Surface) -> Option<Self> {
        let all: Vec<Vec3> = polygons.iter().flatten().copied().collect();
        let points = dedupe_points(&all, surface.tolerance());
        let hull = HullData::build(&points)?;
        let (concavity, witness) = measure_concavity(&hull, &polygons, surface);
        Some(Self {
            polygons,
            hull,
            concavity,
            witness,
            frozen: false,
        })
    }
}

/// Largest distance from a probe on the hull surface to the part, with the
/// worst probe's nearest surface point.
fn measure_concavity(
    hull: &HullData,
    polygons: &[Vec<Vec3>],
    surface: &Surface,
) -> (f64, Option<Witness>) {
    let n = FACET_SAMPLE_RESOLUTION;
    let tolerance = surface.tolerance();
    let mut best = 0.0_f64;
    let mut witness = None;
    for face in &hull.faces {
        let (a, b, c) = (
            hull.vertices[face[0]],
            hull.vertices[face[1]],
            hull.vertices[face[2]],
        );
        // Interior grid points only: facet vertices are mesh vertices.
        for i in 0..=n {
            for j in 0..=(n - i) {
                let k = n - i - j;
                if i == n || j == n || k == n {
                    continue;
                }
                let weights = (
                    i as f64 / n as f64,
                    j as f64 / n as f64,
                    k as f64 / n as f64,
                );
                let probe = add(
                    add(scale(a, weights.0), scale(b, weights.1)),
                    scale(c, weights.2),
                );
                let nearest = nearest_on_polygons(polygons, probe);
                // Probes on the surface, or not worse than the current worst,
                // need no inside test.
                if nearest.distance <= tolerance.max(best) {
                    continue;
                }
                // A hull-boundary point inside the solid lies on a cut cap.
                if surface.contains(probe) {
                    continue;
                }
                best = nearest.distance;
                witness = Some(Witness {
                    nearest: nearest.point,
                    polygon: nearest.polygon,
                });
            }
        }
    }
    (best, witness)
}

// ─── Splitting ────────────────────────────────────────────────────────────

/// Clip a convex polygon to the half-space `normal · x <= offset`.
///
/// Vertices within `eps` of the plane count as lying on it and are kept.
fn clip_polygon(polygon: &[Vec3], normal: Vec3, offset: f64, eps: f64) -> Vec<Vec3> {
    let n = polygon.len();
    let mut out: Vec<Vec3> = Vec::with_capacity(n + 1);
    for i in 0..n {
        let a = polygon[i];
        let b = polygon[(i + 1) % n];
        let da = dot(normal, a) - offset;
        let db = dot(normal, b) - offset;
        let a_in = da <= eps;
        let b_in = db <= eps;
        if a_in {
            out.push(a);
        }
        if a_in != b_in {
            let t = (da / (da - db)).clamp(0.0, 1.0);
            out.push(add(a, scale(sub(b, a), t)));
        }
    }
    // Drop repeated vertices produced by clipping through existing ones.
    let mut cleaned: Vec<Vec3> = Vec::with_capacity(out.len());
    for p in out {
        if cleaned.last().is_none_or(|&q| norm(sub(p, q)) > eps) {
            cleaned.push(p);
        }
    }
    while cleaned.len() > 1 && norm(sub(cleaned[0], cleaned[cleaned.len() - 1])) <= eps {
        cleaned.pop();
    }
    cleaned
}

/// Split the cluster polygons by the plane `normal · x = offset` into the
/// parts below (`normal · x <= offset`) and above it.
///
/// A polygon that lies in the plane is a face of the solid, and the solid is
/// on one side of it.  Only the child on that side owns the polygon: handing
/// it to the other child as well would drag zero-thickness slivers of
/// unrelated geometry into that child's hull.  The side is found by probing
/// the solid just in front of and behind the polygon, so it does not depend on
/// the winding of the input mesh.  If the probes cannot tell (a degenerate
/// mesh), the polygon goes to both children.
fn split_polygons(
    polygons: &[Vec<Vec3>],
    normal: Vec3,
    offset: f64,
    surface: &Surface,
) -> (Vec<Vec<Vec3>>, Vec<Vec<Vec3>>) {
    let eps = surface.tolerance();
    let probe_distance = SIDE_PROBE_FRACTION * surface.diagonal;
    let mut below = Vec::with_capacity(polygons.len());
    let mut above = Vec::with_capacity(polygons.len());
    for polygon in polygons {
        let coplanar = polygon
            .iter()
            .all(|&p| (dot(normal, p) - offset).abs() <= eps);
        if coplanar {
            let centre = scale(
                polygon.iter().fold([0.0; 3], |acc, &p| add(acc, p)),
                1.0 / polygon.len() as f64,
            );
            let solid_below = surface.contains(sub(centre, scale(normal, probe_distance)));
            let solid_above = surface.contains(add(centre, scale(normal, probe_distance)));
            if solid_below || !solid_above {
                below.push(polygon.clone());
            }
            if solid_above || !solid_below {
                above.push(polygon.clone());
            }
            continue;
        }
        let lower = clip_polygon(polygon, normal, offset, eps);
        if lower.len() >= 3 {
            below.push(lower);
        }
        let upper = clip_polygon(polygon, neg(normal), -offset, eps);
        if upper.len() >= 3 {
            above.push(upper);
        }
    }
    (below, above)
}

/// Principal axes of the point set (eigenvectors of its covariance).
fn principal_axes(points: &[Vec3]) -> Vec<Vec3> {
    if points.is_empty() {
        return Vec::new();
    }
    let count = points.len() as f64;
    let mean = scale(
        points.iter().fold([0.0; 3], |acc, &p| add(acc, p)),
        1.0 / count,
    );
    let mut covariance = Mat3::zeros();
    for p in points {
        let d = sub(*p, mean);
        for (row, &dr) in d.iter().enumerate() {
            for (col, &dc) in d.iter().enumerate() {
                covariance[(row, col)] += dr * dc / count;
            }
        }
    }
    let eigen = covariance.symmetric_eigen();
    (0..3)
        .map(|k| {
            let column = eigen.eigenvectors.column(k);
            unit([column[0], column[1], column[2]])
        })
        .collect()
}

/// Add the plane `normal · x = offset` unless an equivalent one is present.
fn push_plane(planes: &mut Vec<(Vec3, f64)>, normal: Vec3, offset: f64, tolerance: f64) {
    if norm(normal) < 0.5 {
        return;
    }
    let duplicate = planes.iter().any(|&(n, d)| {
        let alignment = dot(n, normal);
        (alignment >= 1.0 - 1e-9 && (d - offset).abs() <= tolerance)
            || (alignment <= -1.0 + 1e-9 && (d + offset).abs() <= tolerance)
    });
    if !duplicate {
        planes.push((normal, offset));
    }
}

/// Split planes `(unit normal, offset)` worth trying for `cluster`:
///
/// * the coordinate axes and the principal axes of the cluster as normals,
///   through the most concave point and through the hull centroid;
/// * the planes of the mesh face nearest to the worst probe and of the faces
///   that share a corner with it.  Cutting along the extension of a wall of
///   a notch separates the notch from the rest of the solid, which is how a
///   concave corner is resolved exactly.
fn candidate_planes(cluster: &Cluster, witness: &Witness, surface: &Surface) -> Vec<(Vec3, f64)> {
    let tolerance = surface.tolerance();
    let points: Vec<Vec3> = cluster.polygons.iter().flatten().copied().collect();
    let mut normals: Vec<Vec3> = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    normals.extend(principal_axes(&points));
    let mut planes: Vec<(Vec3, f64)> = Vec::new();
    for normal in normals {
        for anchor in [witness.nearest, cluster.hull.centroid] {
            push_plane(&mut planes, normal, dot(normal, anchor), tolerance);
        }
    }
    // Faces of the cluster that share a corner with the polygon nearest to the
    // worst probe (this includes that polygon itself).
    if let Some(hit) = cluster.polygons.get(witness.polygon) {
        for polygon in &cluster.polygons {
            let touches = polygon
                .iter()
                .any(|&a| hit.iter().any(|&b| norm(sub(a, b)) <= tolerance));
            if !touches {
                continue;
            }
            if let Some((normal, corner)) = polygon_plane(polygon) {
                push_plane(&mut planes, normal, dot(normal, corner), tolerance);
            }
        }
    }
    planes
}

/// Best split of `cluster`: the two children whose concavities add up to the
/// least (ties: smallest total hull volume), or `None` if no candidate plane
/// yields two non-trivial children that both shrink the hull.
fn best_split(cluster: &Cluster, surface: &Surface) -> Option<(Cluster, Cluster)> {
    let witness = cluster.witness?;
    let eps = surface.tolerance();
    let parent_volume = cluster.hull.volume;
    let mut best: Option<(f64, f64, Cluster, Cluster)> = None;
    for (normal, offset) in candidate_planes(cluster, &witness, surface) {
        // The plane must cut through the hull.
        let (lowest, highest) = cluster
            .hull
            .vertices
            .iter()
            .map(|&v| dot(normal, v))
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), d| {
                (lo.min(d), hi.max(d))
            });
        if !(lowest < offset - eps && highest > offset + eps) {
            continue;
        }
        let (below, above) = split_polygons(&cluster.polygons, normal, offset, surface);
        let (Some(low), Some(high)) = (Cluster::new(below, surface), Cluster::new(above, surface))
        else {
            continue;
        };
        let shrinks =
            |child: &Cluster| child.hull.volume < parent_volume * (1.0 - MIN_RELATIVE_VOLUME_DROP);
        if !(shrinks(&low) && shrinks(&high)) {
            continue;
        }
        let score = low.concavity + high.concavity;
        let total_volume = low.hull.volume + high.hull.volume;
        let better = match &best {
            None => true,
            Some((best_score, best_volume, _, _)) => {
                score < best_score - eps
                    || (score <= best_score + eps && total_volume < *best_volume)
            }
        };
        if better {
            best = Some((score, total_volume, low, high));
        }
    }
    best.map(|(_, _, low, high)| (low, high))
}

fn into_part(cluster: Cluster) -> ConvexPart {
    let index = |i: usize| u32::try_from(i).unwrap_or(u32::MAX);
    let indices = cluster
        .hull
        .faces
        .iter()
        .map(|f| [index(f[0]), index(f[1]), index(f[2])])
        .collect();
    ConvexPart {
        vertices: cluster.hull.vertices,
        concavity: cluster.concavity,
        indices,
        volume: cluster.hull.volume,
        centroid: cluster.hull.centroid,
    }
}

/// Decompose the closed mesh `(vertices, triangles)`; see the module docs.
pub(super) fn decompose_mesh(
    config: &AcdConfig,
    vertices: &[Vec3],
    triangles: &[[usize; 3]],
) -> Result<ConvexParts> {
    if !(config.max_concavity.is_finite() && config.max_concavity >= 0.0) {
        return Err(Error::out_of_range(
            "max_concavity",
            config.max_concavity,
            0.0,
            f64::MAX,
        ));
    }
    if !(0.0..=1.0).contains(&config.min_volume_fraction) {
        return Err(Error::out_of_range(
            "min_volume_fraction",
            config.min_volume_fraction,
            0.0,
            1.0,
        ));
    }
    if config.max_parts == 0 {
        return Err(Error::out_of_range(
            "max_parts",
            0.0,
            1.0,
            usize::MAX as f64,
        ));
    }
    let surface = Surface::new(vertices, triangles)?;
    let polygons: Vec<Vec<Vec3>> = surface.triangles.iter().map(|t| t.to_vec()).collect();
    let root = Cluster::new(polygons, &surface)
        .ok_or_else(|| Error::degenerate_geometry("the mesh spans no volume"))?;
    let reference_volume = root.hull.volume;

    let mut clusters = vec![root];
    while clusters.len() < config.max_parts {
        let next = clusters
            .iter()
            .enumerate()
            .filter(|(_, c)| !c.frozen && c.concavity > config.max_concavity)
            .max_by(|(_, a), (_, b)| a.concavity.total_cmp(&b.concavity))
            .map(|(i, _)| i);
        let Some(index) = next else { break };
        match best_split(&clusters[index], &surface) {
            Some((low, high)) => {
                clusters[index] = low;
                clusters.push(high);
            }
            None => clusters[index].frozen = true,
        }
    }
    let threshold = config.min_volume_fraction * reference_volume;
    let parts = clusters
        .into_iter()
        .filter(|c| c.hull.volume >= threshold)
        .map(into_part)
        .collect();
    Ok(ConvexParts { parts })
}
