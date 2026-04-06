#![allow(clippy::type_complexity)]
// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// A view frustum defined by 6 half-space planes.
///
/// Each plane is `[nx, ny, nz, d]` where `n·p + d >= 0` is the "inside" half-space.
/// Planes should be: near, far, left, right, top, bottom.
#[allow(dead_code)]
pub struct ViewFrustum {
    /// The 6 frustum planes.
    pub planes: [[f64; 4]; 6],
}
impl ViewFrustum {
    /// Create a frustum from 6 plane coefficients.
    #[allow(dead_code)]
    pub fn new(planes: [[f64; 4]; 6]) -> Self {
        Self { planes }
    }
    /// Test whether a sphere is potentially visible (intersects the frustum).
    ///
    /// Returns `false` only if the sphere is fully outside at least one plane.
    #[allow(dead_code)]
    pub fn test_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        for plane in &self.planes {
            let n = [plane[0], plane[1], plane[2]];
            let d = plane[3];
            let sd = dot(center, n) + d;
            if sd < -radius {
                return false;
            }
        }
        true
    }
    /// Test whether an AABB is potentially visible.
    ///
    /// Uses the p-vertex / n-vertex method for conservative rejection.
    #[allow(dead_code)]
    pub fn test_aabb(&self, min: [f64; 3], max: [f64; 3]) -> bool {
        for plane in &self.planes {
            let n = [plane[0], plane[1], plane[2]];
            let d = plane[3];
            let p = [
                if n[0] >= 0.0 { max[0] } else { min[0] },
                if n[1] >= 0.0 { max[1] } else { min[1] },
                if n[2] >= 0.0 { max[2] } else { min[2] },
            ];
            if dot(p, n) + d < 0.0 {
                return false;
            }
        }
        true
    }
}
/// Result of a ray cast in the all-pairs query.
#[allow(dead_code)]
pub struct AllPairsRayResult {
    /// Index into the input shapes array of the hit shape.
    pub shape_index: usize,
    /// Time of impact.
    pub toi: f64,
    /// Surface normal at the hit point.
    pub normal: [f64; 3],
    /// Hit point in world space.
    pub point: [f64; 3],
}
/// Result of a sphere cast against the scene.
#[allow(dead_code)]
pub struct SphereCastResult {
    /// Index of the first shape hit.
    pub shape_index: usize,
    /// Time of impact.
    pub toi: f64,
    /// Contact normal at the hit (pointing from the hit shape toward the sphere centre).
    pub normal: [f64; 3],
}
/// Result of a continuous time-of-impact (TOI) query between two moving spheres.
#[allow(dead_code)]
pub struct ToiResult {
    /// Time of impact in `[0, max_t]`.
    pub toi: f64,
    /// Contact normal at impact (from B toward A).
    pub normal: [f64; 3],
    /// Contact point in world space.
    pub point: [f64; 3],
}
/// A shape and its distance from the query point.
#[allow(dead_code)]
pub struct NearestShape {
    /// Index of the shape in the input array.
    pub index: usize,
    /// Signed distance from the query point to the shape surface.
    /// Negative if the point is inside the shape.
    pub distance: f64,
    /// Closest point on the shape surface to the query point.
    pub closest_point: [f64; 3],
}
/// Compound shape composed of transformed child shapes.
pub struct CompoundQuery {
    /// Child shapes together with their world transforms.
    pub shapes: Vec<(Box<dyn ShapeQuery>, [[f64; 4]; 4])>,
}
impl CompoundQuery {
    /// Create an empty compound shape.
    pub fn new() -> Self {
        CompoundQuery { shapes: Vec::new() }
    }
    /// Add a child shape with the given 4×4 world transform.
    pub fn add_shape(&mut self, shape: Box<dyn ShapeQuery>, transform: [[f64; 4]; 4]) {
        self.shapes.push((shape, transform));
    }
    /// Cast a ray against all child shapes and return the first (closest) hit.
    pub fn ray_cast_any(
        &self,
        ray_origin: [f64; 3],
        ray_dir: [f64; 3],
        max_toi: f64,
    ) -> Option<RayCastResult> {
        let mut best: Option<RayCastResult> = None;
        for (shape, transform) in &self.shapes {
            let (local_origin, local_dir) = transform_ray(ray_origin, ray_dir, *transform);
            let current_max = best.as_ref().map_or(max_toi, |b| b.toi);
            if let Some(local_res) = shape.ray_cast(local_origin, local_dir, current_max) {
                let world_point = transform_point(local_res.point, *transform);
                let n = local_res.normal;
                let t = transform;
                let world_normal = normalize([
                    t[0][0] * n[0] + t[0][1] * n[1] + t[0][2] * n[2],
                    t[1][0] * n[0] + t[1][1] * n[1] + t[1][2] * n[2],
                    t[2][0] * n[0] + t[2][1] * n[1] + t[2][2] * n[2],
                ]);
                let candidate = RayCastResult {
                    toi: local_res.toi,
                    normal: world_normal,
                    point: world_point,
                };
                let replace = match &best {
                    None => true,
                    Some(b) => candidate.toi < b.toi,
                };
                if replace {
                    best = Some(candidate);
                }
            }
        }
        best
    }
}
/// Descriptor for a shape in the all-pairs query.
#[allow(dead_code)]
pub enum QueryShapeRef<'a> {
    /// A sphere.
    Sphere(&'a QuerySphere),
    /// An axis-aligned box.
    Box(&'a QueryBox),
    /// A capsule.
    Capsule(&'a QueryCapsule),
}
/// Result of a successful ray cast against a shape.
pub struct RayCastResult {
    /// Time of impact along the ray (parametric distance).
    pub toi: f64,
    /// Surface normal at the hit point (world-space unit vector).
    pub normal: [f64; 3],
    /// Hit point in world space.
    pub point: [f64; 3],
}
/// Sphere shape for queries.
pub struct QuerySphere {
    /// Centre of the sphere.
    pub center: [f64; 3],
    /// Radius.
    pub radius: f64,
}
/// Result of a segment-to-segment closest-point query.
#[allow(dead_code)]
pub struct SegmentDistResult {
    /// Closest point on segment A.
    pub point_a: [f64; 3],
    /// Closest point on segment B.
    pub point_b: [f64; 3],
    /// Distance between the two closest points.
    pub distance: f64,
}
/// Result of a closest-point query on a shape.
pub struct PointQueryResult {
    /// Closest point on (or inside) the shape surface.
    pub point: [f64; 3],
    /// Signed distance: negative when the query point is inside the shape.
    pub distance: f64,
    /// `true` when the query point lies inside the shape.
    pub inside: bool,
}
/// Capsule shape for queries (line segment + radius).
pub struct QueryCapsule {
    /// First endpoint of the central segment.
    pub p0: [f64; 3],
    /// Second endpoint of the central segment.
    pub p1: [f64; 3],
    /// Radius of the capsule.
    pub radius: f64,
}
impl QueryCapsule {
    /// Closest point on the segment p0→p1 to `point`, returning the point and
    /// the parameter t ∈ \[0,1\].
    pub(super) fn closest_point_on_segment(&self, point: [f64; 3]) -> ([f64; 3], f64) {
        let ab = sub(self.p1, self.p0);
        let ap = sub(point, self.p0);
        let len_sq = dot(ab, ab);
        if len_sq < 1e-30 {
            return (self.p0, 0.0);
        }
        let t = clamp(dot(ap, ab) / len_sq, 0.0, 1.0);
        let closest = add(self.p0, scale(ab, t));
        (closest, t)
    }
}
/// Axis-aligned box shape for queries.
pub struct QueryBox {
    /// Centre of the box.
    pub center: [f64; 3],
    /// Half-extents along each axis.
    pub half_extents: [f64; 3],
}
