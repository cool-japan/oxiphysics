//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
#[allow(unused_imports)]
use super::functions_2::*;
pub use super::specialized::*;

/// Result of a ray cast against a shape.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RayCastResult {
    /// Whether the ray hit the shape.
    pub hit: bool,
    /// Time of intersection (units of ray direction magnitude).
    pub toi: f64,
    /// Hit point in world space.
    pub hit_point: [f64; 3],
    /// Outward surface normal at the hit point.
    pub normal: [f64; 3],
}
/// A unified contact result produced by the routing layer.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NarrowPhaseContact {
    /// Contact normal pointing from shape A toward shape B (unit length).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlapping).
    pub depth: f64,
    /// World-space witness point on shape A's surface.
    pub point_a: [f64; 3],
    /// World-space witness point on shape B's surface.
    pub point_b: [f64; 3],
}
impl NarrowPhaseContact {
    /// Return a copy with the normal flipped and witness points swapped.
    #[allow(dead_code)]
    pub fn flipped(&self) -> Self {
        Self {
            normal: scale3(self.normal, -1.0),
            depth: self.depth,
            point_a: self.point_b,
            point_b: self.point_a,
        }
    }
    /// Midpoint between the two witness points.
    #[allow(dead_code)]
    pub fn midpoint(&self) -> [f64; 3] {
        scale3(add3(self.point_a, self.point_b), 0.5)
    }
}
/// A compound shape: a collection of child shapes with local offsets.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CompoundShape {
    /// Child shapes in world space (pre-transformed by the caller).
    pub children: Vec<ShapeKind>,
}
impl CompoundShape {
    /// Create a compound shape from world-space children.
    #[allow(dead_code)]
    pub fn new(children: Vec<ShapeKind>) -> Self {
        CompoundShape { children }
    }
    /// Compute the AABB enclosing all children.
    #[allow(dead_code)]
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for child in &self.children {
            let (cmin, cmax) = child.aabb();
            for i in 0..3 {
                mn[i] = mn[i].min(cmin[i]);
                mx[i] = mx[i].max(cmax[i]);
            }
        }
        (mn, mx)
    }
}
/// Geometric feature that produced a contact (for incremental warm-starting).
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ContactFeature {
    /// Contact between two faces.
    FaceFace {
        /// Index of face A.
        face_a: u32,
        /// Index of face B.
        face_b: u32,
    },
    /// Contact between a face and an edge.
    FaceEdge {
        /// Index of the face.
        face: u32,
        /// Index of the edge.
        edge: u32,
    },
    /// Contact between two edges.
    EdgeEdge {
        /// Index of edge A.
        edge_a: u32,
        /// Index of edge B.
        edge_b: u32,
    },
    /// Contact between a vertex and a face.
    VertexFace {
        /// Index of the vertex.
        vertex: u32,
        /// Index of the face.
        face: u32,
    },
    /// Unknown / unclassified feature (e.g. from GJK fallback).
    Unknown,
}
/// Post-processing options applied to a [`NarrowPhaseContact`] after it is
/// produced by the routing layer.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ContactFilter {
    /// Minimum depth to report a contact. Contacts shallower than this are
    /// discarded (default: 0.0).
    pub min_depth: f64,
    /// Maximum depth that will be reported; deeper contacts are clamped.
    /// Set to `f64::INFINITY` to disable clamping.
    pub max_depth: f64,
    /// If `true`, negate the contact normal before returning it.
    pub flip_normal: bool,
}
impl ContactFilter {
    /// Apply this filter to a contact result.
    ///
    /// Returns `None` if the contact should be discarded.
    #[allow(dead_code)]
    pub fn apply(&self, mut c: NarrowPhaseContact) -> Option<NarrowPhaseContact> {
        if c.depth < self.min_depth {
            return None;
        }
        if c.depth > self.max_depth {
            c.depth = self.max_depth;
        }
        if self.flip_normal {
            c.normal = scale3(c.normal, -1.0);
        }
        Some(c)
    }
}
/// A narrowphase contact enriched with feature information.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct FeatureContact {
    /// Underlying geometry contact.
    pub contact: NarrowPhaseContact,
    /// The geometric feature responsible for this contact.
    pub feature: ContactFeature,
}
impl FeatureContact {
    /// Wrap a plain contact with an `Unknown` feature tag.
    #[allow(dead_code)]
    pub fn from_plain(c: NarrowPhaseContact) -> Self {
        FeatureContact {
            contact: c,
            feature: ContactFeature::Unknown,
        }
    }
}
/// Processes multiple broadphase pairs in a single call.
///
/// Pairs are described as indices into a shared slice of [`ShapeKind`] values.
/// Results are returned in the same order as the input pairs; pairs that did
/// not produce a contact get `None`.
#[allow(dead_code)]
#[derive(Default)]
pub struct BatchNarrowPhase {
    /// Post-processing filter applied to every contact.
    pub filter: ContactFilter,
}
impl BatchNarrowPhase {
    /// Create a `BatchNarrowPhase` with default settings.
    #[allow(dead_code)]
    pub fn new() -> Self {
        Self::default()
    }
    /// Run the narrow phase for all `pairs`.
    ///
    /// Returns one `Option<NarrowPhaseContact>` per input pair.
    #[allow(dead_code)]
    pub fn run(
        &self,
        shapes: &[ShapeKind],
        pairs: &[(usize, usize)],
    ) -> Vec<Option<NarrowPhaseContact>> {
        pairs
            .iter()
            .map(|(i, j)| {
                if *i >= shapes.len() || *j >= shapes.len() {
                    return None;
                }
                let contact = shape_shape_contact(&shapes[*i], &shapes[*j])?;
                self.filter.apply(contact)
            })
            .collect()
    }
    /// Run and collect only pairs that produced a contact.
    #[allow(dead_code)]
    pub fn run_compact(
        &self,
        shapes: &[ShapeKind],
        pairs: &[(usize, usize)],
    ) -> Vec<((usize, usize), NarrowPhaseContact)> {
        pairs
            .iter()
            .filter_map(|(i, j)| {
                if *i >= shapes.len() || *j >= shapes.len() {
                    return None;
                }
                let contact = shape_shape_contact(&shapes[*i], &shapes[*j])?;
                let filtered = self.filter.apply(contact)?;
                Some(((*i, *j), filtered))
            })
            .collect()
    }
}
/// Result of a point-in-shape query.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PointQueryResult {
    /// Whether the point is inside the shape.
    pub is_inside: bool,
    /// Closest point on the shape's surface to the query point.
    pub closest_surface_point: [f64; 3],
    /// Signed distance from the query point to the surface (negative = inside).
    pub signed_distance: f64,
    /// Outward surface normal at the closest point.
    pub normal: [f64; 3],
}
/// A triangle mesh (concave shape) represented as a flat list of triangles.
///
/// Each triangle is three consecutive vertices: `[v0, v1, v2]`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TriangleMesh {
    /// Flat list of triangle vertices; length must be a multiple of 3.
    pub triangles: Vec<[f64; 3]>,
}
impl TriangleMesh {
    /// Create a new triangle mesh.
    #[allow(dead_code)]
    pub fn new(triangles: Vec<[f64; 3]>) -> Self {
        debug_assert!(
            triangles.len().is_multiple_of(3),
            "triangle count must be a multiple of 3"
        );
        TriangleMesh { triangles }
    }
    /// Number of triangles.
    #[allow(dead_code)]
    pub fn tri_count(&self) -> usize {
        self.triangles.len() / 3
    }
    /// Get the three vertices of triangle `i`.
    #[allow(dead_code)]
    pub fn triangle(&self, i: usize) -> [[f64; 3]; 3] {
        let base = i * 3;
        [
            self.triangles[base],
            self.triangles[base + 1],
            self.triangles[base + 2],
        ]
    }
    /// Compute the face normal for triangle `i` (not normalized).
    #[allow(dead_code)]
    pub fn face_normal(&self, i: usize) -> [f64; 3] {
        let [v0, v1, v2] = self.triangle(i);
        let e0 = sub3(v1, v0);
        let e1 = sub3(v2, v0);
        cross3(e0, e1)
    }
}
/// Plain-data descriptor for a convex shape, used by the routing layer.
///
/// All geometry is expressed in world space so no transform arithmetic is
/// needed by the dispatcher.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum ShapeKind {
    /// A sphere with centre and radius.
    Sphere {
        /// World-space centre.
        center: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// An axis-aligned box with centre and half-extents.
    Box {
        /// World-space centre.
        center: [f64; 3],
        /// Half-extents along X, Y, Z.
        half_extents: [f64; 3],
    },
    /// A capsule defined by two endpoint centres and a radius.
    Capsule {
        /// First endpoint centre.
        p0: [f64; 3],
        /// Second endpoint centre.
        p1: [f64; 3],
        /// Radius.
        radius: f64,
    },
    /// An infinite plane through the origin with a given normal.
    Plane {
        /// Unit normal.
        normal: [f64; 3],
        /// Distance from origin along the normal (plane equation: n·x = d).
        offset: f64,
    },
    /// A convex polyhedron described by its support vertices.
    Convex {
        /// Vertices in world space.
        vertices: Vec<[f64; 3]>,
    },
}
impl ShapeKind {
    /// Returns an axis-aligned bounding box `(min, max)` for this shape.
    #[allow(dead_code)]
    pub fn aabb(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            ShapeKind::Sphere { center, radius } => (
                [center[0] - radius, center[1] - radius, center[2] - radius],
                [center[0] + radius, center[1] + radius, center[2] + radius],
            ),
            ShapeKind::Box {
                center,
                half_extents,
            } => (sub3(*center, *half_extents), add3(*center, *half_extents)),
            ShapeKind::Capsule { p0, p1, radius } => {
                let mn = [
                    p0[0].min(p1[0]) - radius,
                    p0[1].min(p1[1]) - radius,
                    p0[2].min(p1[2]) - radius,
                ];
                let mx = [
                    p0[0].max(p1[0]) + radius,
                    p0[1].max(p1[1]) + radius,
                    p0[2].max(p1[2]) + radius,
                ];
                (mn, mx)
            }
            ShapeKind::Plane { normal, offset } => {
                let n = *normal;
                let d = *offset;
                let _ = (n, d);
                ([-1e15; 3], [1e15; 3])
            }
            ShapeKind::Convex { vertices } => {
                let mut mn = [f64::INFINITY; 3];
                let mut mx = [f64::NEG_INFINITY; 3];
                for v in vertices {
                    for i in 0..3 {
                        mn[i] = mn[i].min(v[i]);
                        mx[i] = mx[i].max(v[i]);
                    }
                }
                (mn, mx)
            }
        }
    }
    /// The largest sphere that bounds this shape (bounding sphere).
    #[allow(dead_code)]
    pub fn bounding_radius(&self) -> f64 {
        match self {
            ShapeKind::Sphere { radius, .. } => *radius,
            ShapeKind::Box { half_extents, .. } => len3(*half_extents),
            ShapeKind::Capsule { p0, p1, radius } => len3(sub3(*p1, *p0)) * 0.5 + radius,
            ShapeKind::Plane { .. } => f64::INFINITY,
            ShapeKind::Convex { vertices } => {
                vertices.iter().map(|v| len3(*v)).fold(0.0_f64, f64::max)
            }
        }
    }
    /// Support function: point on the shape furthest in direction `dir`.
    #[allow(dead_code)]
    pub fn support(&self, dir: [f64; 3]) -> [f64; 3] {
        match self {
            ShapeKind::Sphere { center, radius } => {
                let d = normalize3(dir);
                add3(*center, scale3(d, *radius))
            }
            ShapeKind::Box {
                center,
                half_extents,
            } => [
                center[0] + half_extents[0] * dir[0].signum(),
                center[1] + half_extents[1] * dir[1].signum(),
                center[2] + half_extents[2] * dir[2].signum(),
            ],
            ShapeKind::Capsule { p0, p1, radius } => {
                let d0 = dot3(*p0, dir);
                let d1 = dot3(*p1, dir);
                let base = if d0 > d1 { *p0 } else { *p1 };
                add3(base, scale3(normalize3(dir), *radius))
            }
            ShapeKind::Plane { normal, offset } => scale3(*normal, *offset + 1e9),
            ShapeKind::Convex { vertices } => vertices
                .iter()
                .max_by(|a, b| {
                    dot3(**a, dir)
                        .partial_cmp(&dot3(**b, dir))
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .copied()
                .unwrap_or([0.0; 3]),
        }
    }
}
/// Result of a segment cast against a shape.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SegmentCastResult {
    /// Whether the segment hit the shape.
    pub hit: bool,
    /// Parameter along the segment `[0, 1]` where the hit occurred.
    pub t: f64,
    /// Hit point in world space.
    pub hit_point: [f64; 3],
    /// Outward surface normal at the hit point.
    pub normal: [f64; 3],
}
