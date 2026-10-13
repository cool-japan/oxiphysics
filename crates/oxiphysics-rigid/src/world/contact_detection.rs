// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Collider-driven collision detection for [`super::PhysicsWorld`].
//!
//! Every non-sensor [`Collider`] attached to a body is flattened into one or
//! more *leaves* (a compound shape contributes one leaf per child, recursively),
//! each carrying its world transform and world-space AABB. Candidate pairs come
//! from the `oxiphysics-collision` sweep-and-prune broadphase over those AABBs,
//! and contact manifolds come from its [`NarrowPhaseDispatcher`] (analytic
//! sphere / box / capsule algorithms, convex face-clip manifolds and GJK/EPA for
//! every other convex pair).
//!
//! Concave shapes ([`TriangleMesh`], [`HeightField`] and [`Torus`]) are never
//! sent through convex GJK (their support mapping describes the convex hull,
//! which would fill holes and valleys). They are decomposed into triangles,
//! cached per shape together with a BVH over the triangle AABBs, and each
//! overlapping triangle is collided against the convex partner. Height-field and
//! torus triangles are one-sided (their outward orientation is known), so the
//! contact normal always pushes the partner out of the surface. Pairs in which
//! *both* leaves are concave produce no contacts: at least one side of every
//! contact must be convex (or a compound of convex children).
//!
//! All produced [`ContactPair`] normals point FROM body B TOWARD body A, with a
//! positive `depth` for penetration (negative for speculative contacts).

use super::types::ContactPair;
use crate::body::{BodyState, BodyType, RigidBody};
use crate::collider::Collider;
use crate::sets::{ColliderSet, RigidBodySet};
use oxiphysics_collision::narrowphase::dispatch::ShapeType;
use oxiphysics_collision::narrowphase::gjk::{Gjk, GjkResult};
use oxiphysics_collision::types::CollisionPair;
use oxiphysics_collision::{BroadPhase, BvhBroadphase, NarrowPhaseDispatcher, SweepAndPrune};
use oxiphysics_core::math::{Mat3, Real, Vec3};
use oxiphysics_core::{Aabb, BodyHandle, ColliderHandle, Transform};
use oxiphysics_geometry::{
    BoxShape, Capsule, Compound, Cone, ConvexHull, Cylinder, HeightField, RayHit, Shape, Sphere,
    Torus, TriangleMesh,
};
use std::collections::HashMap;
use std::sync::Arc;

/// Dispatcher tag for user-defined convex shapes (anything implementing
/// [`Shape`] that is not one of the geometry crate's primitives). They are
/// collided through their support mapping by GJK/EPA.
const USER_CONVEX_TAG: u32 = 1;

/// Dispatcher tag for the per-triangle proxies of concave shapes.
const TRIANGLE_TAG: u32 = 2;

/// Number of segments around the main ring of a tessellated [`Torus`].
const TORUS_MAJOR_SEGMENTS: usize = 64;

/// Number of segments around the tube of a tessellated [`Torus`].
const TORUS_MINOR_SEGMENTS: usize = 32;

/// Overlap used when probing a full manifold for a speculative (separated)
/// pair; small enough to stay in the contact-feature regime of the pair.
const SPECULATIVE_PROBE_OVERLAP: Real = 1e-4;

/// Triangles whose doubled area is below this are degenerate and skipped.
const DEGENERATE_AREA2: Real = 1e-18;

// ── Leaves ────────────────────────────────────────────────────────────────────

/// Concave shape families supported by the triangle path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ConcaveKind {
    /// Arbitrary triangle soup; winding is user-defined so contacts are two-sided.
    Mesh,
    /// Regular height grid; triangles face local `+y` (one-sided).
    HeightField,
    /// Tessellated torus surface with outward winding (one-sided).
    Torus,
}

impl ConcaveKind {
    fn one_sided(self) -> bool {
        matches!(self, ConcaveKind::HeightField | ConcaveKind::Torus)
    }
}

/// How a leaf is collided.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum LeafKind {
    /// Convex leaf dispatched with the given shape type.
    Convex(ShapeType),
    /// Concave leaf decomposed into cached triangles.
    Concave(ConcaveKind),
}

/// One collision primitive of a collider, positioned for the current step.
#[derive(Debug, Clone)]
pub(crate) struct ColliderLeaf {
    /// Collider the leaf belongs to.
    pub(crate) collider: ColliderHandle,
    /// Body the collider is attached to.
    pub(crate) body: BodyHandle,
    /// Leaf geometry (a compound child, or the collider shape itself).
    pub(crate) shape: Arc<dyn Shape>,
    /// Dispatch classification.
    pub(crate) kind: LeafKind,
    /// World transform of the leaf.
    pub(crate) world: Transform,
    /// World-space bounding box of the leaf.
    pub(crate) aabb: Aabb,
    /// Collider friction coefficient.
    pub(crate) friction: Real,
    /// Collider coefficient of restitution.
    pub(crate) restitution: Real,
}

/// Classification result for a shape.
enum Classified<'a> {
    Leaf(LeafKind),
    Compound(&'a Compound),
}

/// Classify a shape by its concrete type.
fn classify(shape: &dyn Shape) -> Classified<'_> {
    let any = shape.as_any();
    if any.is::<Sphere>() {
        Classified::Leaf(LeafKind::Convex(ShapeType::Sphere))
    } else if any.is::<BoxShape>() {
        Classified::Leaf(LeafKind::Convex(ShapeType::Box))
    } else if any.is::<Capsule>() {
        Classified::Leaf(LeafKind::Convex(ShapeType::Capsule))
    } else if any.is::<Cylinder>() {
        Classified::Leaf(LeafKind::Convex(ShapeType::Cylinder))
    } else if any.is::<Cone>() {
        Classified::Leaf(LeafKind::Convex(ShapeType::Cone))
    } else if any.is::<ConvexHull>() {
        Classified::Leaf(LeafKind::Convex(ShapeType::ConvexHull))
    } else if any.is::<TriangleMesh>() {
        Classified::Leaf(LeafKind::Concave(ConcaveKind::Mesh))
    } else if any.is::<HeightField>() {
        Classified::Leaf(LeafKind::Concave(ConcaveKind::HeightField))
    } else if any.is::<Torus>() {
        Classified::Leaf(LeafKind::Concave(ConcaveKind::Torus))
    } else if let Some(compound) = any.downcast_ref::<Compound>() {
        Classified::Compound(compound)
    } else {
        Classified::Leaf(LeafKind::Convex(ShapeType::Custom(USER_CONVEX_TAG)))
    }
}

/// World-space AABB of `shape` placed at `tf` (exact for the local box: the
/// centre is transformed and the half-extents are mapped through `|R|`).
pub(crate) fn world_aabb(shape: &dyn Shape, tf: &Transform) -> Aabb {
    let local = shape.bounding_box();
    transformed_box(&local, tf)
}

/// Bound the box `local` (expressed in the frame of `tf`) in the parent frame.
fn transformed_box(local: &Aabb, tf: &Transform) -> Aabb {
    let centre = (local.min + local.max) * 0.5;
    let half = (local.max - local.min) * 0.5;
    let rot = tf.rotation.to_rotation_matrix();
    let abs_rot = rot.matrix().abs();
    let world_centre = tf.transform_point(&centre);
    let world_half = abs_rot * half;
    Aabb::new(world_centre - world_half, world_centre + world_half)
}

/// Flatten every non-sensor, body-attached collider into leaves.
///
/// Sensors exert no forces and the world exposes no trigger events, so they are
/// excluded. A collider without a parent body cannot be represented in a
/// [`ContactPair`] (which addresses bodies) and is skipped as well.
pub(crate) fn gather_leaves(bodies: &RigidBodySet, colliders: &ColliderSet) -> Vec<ColliderLeaf> {
    let mut leaves = Vec::new();
    for (handle, collider) in colliders.iter() {
        if collider.is_sensor {
            continue;
        }
        let Some(body_handle) = collider.body_handle else {
            continue;
        };
        let Some(body) = bodies.get(body_handle) else {
            continue;
        };
        let world = collider.world_transform(&body.transform);
        push_leaves(
            &collider.shape,
            world,
            handle,
            body_handle,
            collider,
            &mut leaves,
        );
    }
    leaves
}

fn push_leaves(
    shape: &Arc<dyn Shape>,
    world: Transform,
    handle: ColliderHandle,
    body: BodyHandle,
    collider: &Collider,
    out: &mut Vec<ColliderLeaf>,
) {
    match classify(shape.as_ref()) {
        Classified::Compound(compound) => {
            for (child_tf, child) in &compound.children {
                let child_world = world.compose(child_tf);
                push_leaves(child, child_world, handle, body, collider, out);
            }
        }
        Classified::Leaf(kind) => {
            let aabb = world_aabb(shape.as_ref(), &world);
            out.push(ColliderLeaf {
                collider: handle,
                body,
                shape: Arc::clone(shape),
                kind,
                world,
                aabb,
                friction: collider.friction,
                restitution: collider.restitution,
            });
        }
    }
}

/// Whether the body is simulated this step (awake dynamic, or kinematic).
fn is_moving(body: &RigidBody) -> bool {
    match body.body_type {
        BodyType::Dynamic => body.state != BodyState::Sleeping,
        BodyType::Kinematic => true,
        BodyType::Static => false,
    }
}

/// Pick the sweep axis with the largest spread of AABB centres.
fn sweep_axis(aabbs: &[Aabb]) -> usize {
    if aabbs.len() < 2 {
        return 0;
    }
    let n = aabbs.len() as Real;
    let mut mean = Vec3::zeros();
    for b in aabbs {
        mean += (b.min + b.max) * 0.5;
    }
    mean /= n;
    let mut var = Vec3::zeros();
    for b in aabbs {
        let d = (b.min + b.max) * 0.5 - mean;
        var += d.component_mul(&d);
    }
    if var.x >= var.y && var.x >= var.z {
        0
    } else if var.y >= var.z {
        1
    } else {
        2
    }
}

// ── Concave triangle cache ────────────────────────────────────────────────────

/// Triangle decomposition of one concave shape plus a BVH over the triangles.
struct ConcaveEntry {
    /// Keeps the shape alive so its address (the cache key) cannot be reused.
    _shape: Arc<dyn Shape>,
    /// Triangles in the shape's local frame.
    triangles: Vec<[Vec3; 3]>,
    /// BVH over the local triangle AABBs.
    bvh: BvhBroadphase,
    /// Whether the shape was referenced during the current step.
    seen: bool,
}

/// Per-world cache of concave-shape triangle decompositions, keyed by the
/// shape allocation. Entries for shapes no longer referenced by any collider
/// are evicted at the end of each detection pass.
#[derive(Default)]
pub(crate) struct ConcaveCache {
    entries: HashMap<usize, ConcaveEntry>,
}

fn shape_key(shape: &Arc<dyn Shape>) -> usize {
    Arc::as_ptr(shape).cast::<()>() as usize
}

impl ConcaveCache {
    /// Ensure every concave leaf has a cache entry and evict stale entries.
    fn prepare(&mut self, leaves: &[ColliderLeaf]) {
        for entry in self.entries.values_mut() {
            entry.seen = false;
        }
        for leaf in leaves {
            let LeafKind::Concave(kind) = leaf.kind else {
                continue;
            };
            let key = shape_key(&leaf.shape);
            if let Some(entry) = self.entries.get_mut(&key) {
                entry.seen = true;
                continue;
            }
            let triangles = triangulate(leaf.shape.as_ref(), kind);
            let tri_aabbs: Vec<Aabb> = triangles.iter().map(triangle_aabb).collect();
            let mut bvh = BvhBroadphase::new();
            bvh.build(&tri_aabbs);
            self.entries.insert(
                key,
                ConcaveEntry {
                    _shape: Arc::clone(&leaf.shape),
                    triangles,
                    bvh,
                    seen: true,
                },
            );
        }
        self.entries.retain(|_, e| e.seen);
    }

    fn get(&self, shape: &Arc<dyn Shape>) -> Option<&ConcaveEntry> {
        self.entries.get(&shape_key(shape))
    }

    /// Number of cached concave shapes.
    #[cfg(test)]
    pub(crate) fn len(&self) -> usize {
        self.entries.len()
    }
}

fn triangle_aabb(t: &[Vec3; 3]) -> Aabb {
    let min = t[0].inf(&t[1]).inf(&t[2]);
    let max = t[0].sup(&t[1]).sup(&t[2]);
    Aabb::new(min, max)
}

fn push_triangle(out: &mut Vec<[Vec3; 3]>, a: Vec3, b: Vec3, c: Vec3) {
    if (b - a).cross(&(c - a)).norm_squared() > DEGENERATE_AREA2 {
        out.push([a, b, c]);
    }
}

/// Decompose a concave shape into local-frame triangles.
fn triangulate(shape: &dyn Shape, kind: ConcaveKind) -> Vec<[Vec3; 3]> {
    let mut out = Vec::new();
    match kind {
        ConcaveKind::Mesh => {
            if let Some(mesh) = shape.as_any().downcast_ref::<TriangleMesh>() {
                for idx in &mesh.indices {
                    let (Some(a), Some(b), Some(c)) = (
                        mesh.vertices.get(idx[0]),
                        mesh.vertices.get(idx[1]),
                        mesh.vertices.get(idx[2]),
                    ) else {
                        continue;
                    };
                    push_triangle(&mut out, *a, *b, *c);
                }
            }
        }
        ConcaveKind::HeightField => {
            if let Some(hf) = shape.as_any().downcast_ref::<HeightField>() {
                triangulate_height_field(hf, &mut out);
            }
        }
        ConcaveKind::Torus => {
            if let Some(torus) = shape.as_any().downcast_ref::<Torus>() {
                triangulate_torus(torus, &mut out);
            }
        }
    }
    out
}

/// Two triangles per grid cell, wound so their normals face local `+y`.
fn triangulate_height_field(hf: &HeightField, out: &mut Vec<[Vec3; 3]>) {
    if hf.rows < 2 || hf.cols < 2 {
        return;
    }
    let point = |row: usize, col: usize| -> Option<Vec3> {
        let h = *hf.heights.get(row * hf.cols + col)?;
        Some(Vec3::new(
            col as Real * hf.scale_x,
            h,
            row as Real * hf.scale_z,
        ))
    };
    for row in 0..hf.rows - 1 {
        for col in 0..hf.cols - 1 {
            let (Some(p00), Some(p10), Some(p01), Some(p11)) = (
                point(row, col),
                point(row + 1, col),
                point(row, col + 1),
                point(row + 1, col + 1),
            ) else {
                continue;
            };
            push_triangle(out, p00, p10, p01);
            push_triangle(out, p10, p11, p01);
        }
    }
}

/// Tessellate the torus surface (ring in the local XZ plane, axis `+y`) with
/// outward-facing triangles. The vertices lie on the exact surface, so the
/// facets deviate from it by at most
/// `R(1 − cos(π/M)) + r(1 − cos(π/N))` for `M`/`N` major/minor segments.
fn triangulate_torus(torus: &Torus, out: &mut Vec<[Vec3; 3]>) {
    let big_r = torus.major_radius;
    let small_r = torus.minor_radius;
    if !(big_r > 0.0 && small_r > 0.0) {
        return;
    }
    let m = TORUS_MAJOR_SEGMENTS;
    let n = TORUS_MINOR_SEGMENTS;
    let vertex = |i: usize, j: usize| -> (Vec3, Vec3) {
        let theta = std::f64::consts::TAU * (i % m) as Real / m as Real;
        let phi = std::f64::consts::TAU * (j % n) as Real / n as Real;
        let radial = Vec3::new(theta.cos(), 0.0, theta.sin());
        let ring = radial * big_r;
        let offset = radial * (small_r * phi.cos()) + Vec3::new(0.0, small_r * phi.sin(), 0.0);
        (ring + offset, ring)
    };
    for i in 0..m {
        for j in 0..n {
            let (a, ring_a) = vertex(i, j);
            let (b, _) = vertex(i + 1, j);
            let (c, _) = vertex(i + 1, j + 1);
            let (d, _) = vertex(i, j + 1);
            for (p, q, r) in [(a, b, c), (a, c, d)] {
                let normal = (q - p).cross(&(r - p));
                // Orient outward: away from the tube centre line.
                let centroid = (p + q + r) / 3.0;
                if normal.dot(&(centroid - ring_a)) >= 0.0 {
                    push_triangle(out, p, q, r);
                } else {
                    push_triangle(out, p, r, q);
                }
            }
        }
    }
}

// ── Triangle proxy shape ──────────────────────────────────────────────────────

/// A single triangle exposed as a [`Shape`] so the GJK/EPA machinery of the
/// collision crate can collide it against a convex partner.
#[derive(Debug, Clone)]
pub(crate) struct TriangleShape {
    vertices: [Vec3; 3],
}

impl TriangleShape {
    fn new(vertices: [Vec3; 3]) -> Self {
        Self { vertices }
    }

    fn centroid(&self) -> Vec3 {
        (self.vertices[0] + self.vertices[1] + self.vertices[2]) / 3.0
    }
}

impl Shape for TriangleShape {
    fn bounding_box(&self) -> Aabb {
        triangle_aabb(&self.vertices)
    }

    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn support_point(&self, direction: &Vec3) -> Vec3 {
        let mut best = self.vertices[0];
        let mut best_dot = best.dot(direction);
        for v in &self.vertices[1..] {
            let d = v.dot(direction);
            if d > best_dot {
                best_dot = d;
                best = *v;
            }
        }
        best
    }

    fn volume(&self) -> Real {
        0.0
    }

    fn center_of_mass(&self) -> Vec3 {
        self.centroid()
    }

    /// Inertia of a uniform triangular lamina about its centroid:
    /// `C = m/12 · Σ vᵢvᵢᵀ` (vertices relative to the centroid), `I = tr(C)·1 − C`.
    fn inertia_tensor(&self, mass: Real) -> Mat3 {
        let c = self.centroid();
        let mut cov = Mat3::zeros();
        for v in &self.vertices {
            let r = v - c;
            cov += r * r.transpose();
        }
        cov *= mass / 12.0;
        Mat3::identity() * cov.trace() - cov
    }

    /// Möller–Trumbore ray/triangle intersection (two-sided).
    fn ray_cast(&self, ray_origin: &Vec3, ray_direction: &Vec3, max_toi: Real) -> Option<RayHit> {
        let [a, b, c] = self.vertices;
        let e1 = b - a;
        let e2 = c - a;
        let p = ray_direction.cross(&e2);
        let det = e1.dot(&p);
        if det.abs() < 1e-14 {
            return None;
        }
        let inv_det = 1.0 / det;
        let s = ray_origin - a;
        let u = s.dot(&p) * inv_det;
        if !(0.0..=1.0).contains(&u) {
            return None;
        }
        let q = s.cross(&e1);
        let v = ray_direction.dot(&q) * inv_det;
        if v < 0.0 || u + v > 1.0 {
            return None;
        }
        let toi = e2.dot(&q) * inv_det;
        if toi < 0.0 || toi > max_toi {
            return None;
        }
        let mut normal = e1.cross(&e2).normalize();
        if normal.dot(ray_direction) > 0.0 {
            normal = -normal;
        }
        Some(RayHit {
            point: ray_origin + ray_direction * toi,
            normal,
            toi,
        })
    }
}

// ── Narrowphase ───────────────────────────────────────────────────────────────

/// A raw contact before it is tagged with body/collider metadata.
struct RawContact {
    /// Unit normal pointing from leaf B toward leaf A.
    normal: Vec3,
    /// Penetration depth (negative = speculative gap).
    depth: Real,
    /// World-space contact point.
    point: Vec3,
}

/// Output of [`detect_contacts`].
pub(crate) struct DetectionOutput {
    /// All contacts of the step (normal from B toward A).
    pub(crate) contacts: Vec<ContactPair>,
    /// Sleeping dynamic bodies that are touched by a moving body.
    pub(crate) wake: Vec<BodyHandle>,
}

/// Shared narrowphase context for one candidate pair.
struct NarrowContext<'a> {
    dispatcher: &'a NarrowPhaseDispatcher,
    pair: CollisionPair,
    speculative_margin: Option<Real>,
}

/// A convex operand of a narrowphase query.
struct ConvexOperand<'a> {
    shape: &'a dyn Shape,
    shape_type: ShapeType,
    world: &'a Transform,
}

/// The concave operand of a concave–convex query.
struct ConcaveOperand<'a> {
    entry: &'a ConcaveEntry,
    leaf: &'a ColliderLeaf,
    kind: ConcaveKind,
    /// Whether the concave leaf is pair member A.
    is_a: bool,
}

/// Convex leaf pair → raw contacts (manifold, or speculative GJK distance).
fn convex_convex(
    ctx: &NarrowContext<'_>,
    a: &ConvexOperand<'_>,
    b: &ConvexOperand<'_>,
    out: &mut Vec<RawContact>,
) {
    let (shape_a, tf_a, shape_b, tf_b) = (a.shape, a.world, b.shape, b.world);
    let result = ctx.dispatcher.dispatch(
        shape_a,
        a.shape_type,
        tf_a,
        shape_b,
        b.shape_type,
        tf_b,
        ctx.pair,
    );
    if let Some(manifold) = result.manifold {
        let before = out.len();
        for c in &manifold.contacts {
            let len = c.normal.norm();
            if !(len.is_finite() && len > 0.5 && c.depth.is_finite()) {
                continue;
            }
            out.push(RawContact {
                normal: c.normal / len,
                depth: c.depth.max(0.0),
                point: (c.point_a + c.point_b) * 0.5,
            });
        }
        if out.len() > before {
            return;
        }
    }
    let Some(margin) = ctx.speculative_margin else {
        return;
    };
    if let GjkResult::Separated {
        point_a,
        point_b,
        distance,
    } = Gjk::query(shape_a, tf_a, shape_b, tf_b)
        && distance > 1e-12
        && distance <= margin
    {
        let gap = point_a - point_b;
        let gap_len = gap.norm();
        if gap_len <= 1e-12 {
            return;
        }
        let normal = gap / gap_len;
        // Probe for a full speculative manifold: translate A along -n until it
        // overlaps B by `SPECULATIVE_PROBE_OVERLAP`, collect the manifold
        // there, and shift every point back. A face-parallel pair then gets its
        // whole contact patch (not one GJK witness that would add spin).
        let shift = distance + SPECULATIVE_PROBE_OVERLAP;
        let probe_tf = Transform::new(tf_a.position - normal * shift, tf_a.rotation);
        let probe = ctx.dispatcher.dispatch(
            shape_a,
            a.shape_type,
            &probe_tf,
            shape_b,
            b.shape_type,
            tf_b,
            ctx.pair,
        );
        let before = out.len();
        if let Some(manifold) = probe.manifold {
            for c in &manifold.contacts {
                let len = c.normal.norm();
                if !(len.is_finite() && len > 0.5 && c.depth.is_finite()) {
                    continue;
                }
                let depth = c.depth - shift;
                if depth >= 0.0 || depth < -margin {
                    continue;
                }
                out.push(RawContact {
                    normal: c.normal / len,
                    depth,
                    point: (c.point_a + normal * shift + c.point_b) * 0.5,
                });
            }
        }
        if out.len() == before {
            out.push(RawContact {
                normal,
                depth: -distance,
                point: (point_a + point_b) * 0.5,
            });
        }
    }
}

/// Deepest point of a convex leaf below the plane `n·x = offset` (world frame).
fn plane_penetration(shape: &dyn Shape, tf: &Transform, n: &Vec3, offset: Real) -> (Real, Vec3) {
    let local_dir = tf.rotation.inverse() * (-n);
    let deepest = tf.transform_point(&shape.support_point(&local_dir));
    (offset - n.dot(&deepest), deepest)
}

/// Concave leaf vs convex leaf → raw contacts, oriented for the original pair
/// order (`concave.is_a` says whether the concave leaf is pair member A).
fn concave_convex(
    ctx: &NarrowContext<'_>,
    concave: &ConcaveOperand<'_>,
    convex: &ColliderLeaf,
    convex_type: ShapeType,
    out: &mut Vec<RawContact>,
) {
    let margin = ctx.speculative_margin.unwrap_or(0.0);
    let concave_world = &concave.leaf.world;
    let query = transformed_box(&convex.aabb, &concave_world.inverse()).expand(margin);
    let rotation = concave_world.rotation.to_rotation_matrix();
    let convex_op = ConvexOperand {
        shape: convex.shape.as_ref(),
        shape_type: convex_type,
        world: &convex.world,
    };
    for tri_index in concave.entry.bvh.query(&query) {
        let Some(tri) = concave.entry.triangles.get(tri_index) else {
            continue;
        };
        let proxy = TriangleShape::new(*tri);
        let tri_op = ConvexOperand {
            shape: &proxy,
            shape_type: ShapeType::Custom(TRIANGLE_TAG),
            world: concave_world,
        };
        let before = out.len();
        if concave.is_a {
            convex_convex(ctx, &tri_op, &convex_op, out);
        } else {
            convex_convex(ctx, &convex_op, &tri_op, out);
        }
        if !concave.kind.one_sided() || out.len() == before {
            continue;
        }
        // One-sided surface: the convex partner must always be pushed out along
        // the outward face normal.
        let face = (tri[1] - tri[0]).cross(&(tri[2] - tri[0]));
        let face_len = face.norm();
        if face_len <= 0.0 {
            continue;
        }
        let face_n = rotation * (face / face_len);
        let expected = if concave.is_a { -face_n } else { face_n };
        let offset = face_n.dot(&concave_world.transform_point(&tri[0]));
        let mut keep = Vec::with_capacity(out.len() - before);
        for raw in out.drain(before..) {
            if raw.normal.dot(&expected) > 0.0 {
                keep.push(raw);
                continue;
            }
            if raw.depth < 0.0 {
                // A speculative gap measured to the back side is not a contact.
                continue;
            }
            let (depth, deepest) =
                plane_penetration(convex.shape.as_ref(), &convex.world, &face_n, offset);
            if depth > 0.0 {
                keep.push(RawContact {
                    normal: expected,
                    depth,
                    point: deepest + face_n * (0.5 * depth),
                });
            }
        }
        out.extend(keep);
    }
}

/// Run broadphase + narrowphase over all colliders.
///
/// `speculative_margin` (when `Some`) inflates the broadphase AABBs and admits
/// separated convex pairs closer than the margin as negative-depth contacts.
pub(crate) fn detect_contacts(
    bodies: &RigidBodySet,
    colliders: &ColliderSet,
    dispatcher: &NarrowPhaseDispatcher,
    cache: &mut ConcaveCache,
    speculative_margin: Option<Real>,
) -> DetectionOutput {
    let leaves = gather_leaves(bodies, colliders);
    cache.prepare(&leaves);
    let margin = speculative_margin.unwrap_or(0.0).max(0.0);
    let aabbs: Vec<Aabb> = leaves
        .iter()
        .map(|l| {
            if margin > 0.0 {
                l.aabb.expand(margin * 0.5)
            } else {
                l.aabb.clone()
            }
        })
        .collect();
    let pairs = SweepAndPrune::new(sweep_axis(&aabbs)).find_pairs(&aabbs);

    let mut contacts = Vec::new();
    let mut wake = Vec::new();
    let mut raw = Vec::new();
    for pair in pairs {
        let (Some(leaf_a), Some(leaf_b)) = (leaves.get(pair.a), leaves.get(pair.b)) else {
            continue;
        };
        if leaf_a.body == leaf_b.body {
            continue;
        }
        let (Some(body_a), Some(body_b)) = (bodies.get(leaf_a.body), bodies.get(leaf_b.body))
        else {
            continue;
        };
        let any_dynamic =
            body_a.body_type == BodyType::Dynamic || body_b.body_type == BodyType::Dynamic;
        if !any_dynamic || !(is_moving(body_a) || is_moving(body_b)) {
            continue;
        }
        raw.clear();
        let ctx = NarrowContext {
            dispatcher,
            pair,
            speculative_margin,
        };
        match (leaf_a.kind, leaf_b.kind) {
            (LeafKind::Convex(type_a), LeafKind::Convex(type_b)) => convex_convex(
                &ctx,
                &ConvexOperand {
                    shape: leaf_a.shape.as_ref(),
                    shape_type: type_a,
                    world: &leaf_a.world,
                },
                &ConvexOperand {
                    shape: leaf_b.shape.as_ref(),
                    shape_type: type_b,
                    world: &leaf_b.world,
                },
                &mut raw,
            ),
            (LeafKind::Concave(kind), LeafKind::Convex(type_b)) => {
                if let Some(entry) = cache.get(&leaf_a.shape) {
                    let concave = ConcaveOperand {
                        entry,
                        leaf: leaf_a,
                        kind,
                        is_a: true,
                    };
                    concave_convex(&ctx, &concave, leaf_b, type_b, &mut raw);
                }
            }
            (LeafKind::Convex(type_a), LeafKind::Concave(kind)) => {
                if let Some(entry) = cache.get(&leaf_b.shape) {
                    let concave = ConcaveOperand {
                        entry,
                        leaf: leaf_b,
                        kind,
                        is_a: false,
                    };
                    concave_convex(&ctx, &concave, leaf_a, type_a, &mut raw);
                }
            }
            (LeafKind::Concave(_), LeafKind::Concave(_)) => {}
        }
        if raw.is_empty() {
            continue;
        }
        if raw.iter().any(|c| c.depth >= 0.0) {
            for (handle, body, other) in
                [(leaf_a.body, body_a, body_b), (leaf_b.body, body_b, body_a)]
            {
                let asleep =
                    body.body_type == BodyType::Dynamic && body.state == BodyState::Sleeping;
                if asleep && is_moving(other) && !wake.contains(&handle) {
                    wake.push(handle);
                }
            }
        }
        let friction = (leaf_a.friction.max(0.0) * leaf_b.friction.max(0.0)).sqrt();
        let restitution = leaf_a.restitution.max(leaf_b.restitution).clamp(0.0, 1.0);
        for c in raw.drain(..) {
            contacts.push(ContactPair {
                body_a: leaf_a.body,
                body_b: leaf_b.body,
                collider_a: leaf_a.collider,
                collider_b: leaf_b.collider,
                normal: [c.normal.x, c.normal.y, c.normal.z],
                depth: c.depth,
                contact_point: [c.point.x, c.point.y, c.point.z],
                friction,
                restitution,
            });
        }
    }
    DetectionOutput { contacts, wake }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::body::RigidBody;
    use oxiphysics_core::math::Quat;

    fn unit_axis(v: Vec3) -> nalgebra::Unit<Vec3> {
        nalgebra::Unit::new_normalize(v)
    }

    #[test]
    fn world_aabb_of_rotated_box_swaps_extents() {
        let shape = BoxShape::new(Vec3::new(1.0, 0.5, 0.25));
        let tf = Transform::new(
            Vec3::new(1.0, 2.0, 3.0),
            Quat::from_axis_angle(&unit_axis(Vec3::z()), std::f64::consts::FRAC_PI_2),
        );
        let aabb = world_aabb(&shape, &tf);
        let half = (aabb.max - aabb.min) * 0.5;
        let centre = (aabb.max + aabb.min) * 0.5;
        assert!((half - Vec3::new(0.5, 1.0, 0.25)).norm() < 1e-12, "{half}");
        assert!((centre - Vec3::new(1.0, 2.0, 3.0)).norm() < 1e-12);
    }

    #[test]
    fn height_field_triangles_face_up() {
        let hf = HeightField::new((0..12).map(|i| i as Real * 0.1).collect(), 3, 4, 1.0, 2.0);
        let tris = triangulate(&hf, ConcaveKind::HeightField);
        assert_eq!(tris.len(), 2 * 2 * 3);
        for t in &tris {
            let n = (t[1] - t[0]).cross(&(t[2] - t[0]));
            assert!(n.y > 0.0, "height-field facet must face +y: {n}");
        }
    }

    #[test]
    fn torus_triangles_lie_on_the_surface_and_face_outward() {
        let torus = Torus::new(2.0, 0.5);
        let tris = triangulate(&torus, ConcaveKind::Torus);
        assert_eq!(tris.len(), 2 * TORUS_MAJOR_SEGMENTS * TORUS_MINOR_SEGMENTS);
        for t in &tris {
            for v in t {
                let ring = Vec3::new(v.x, 0.0, v.z).normalize() * 2.0;
                assert!(((v - ring).norm() - 0.5).abs() < 1e-9, "vertex off surface");
            }
            let centroid = (t[0] + t[1] + t[2]) / 3.0;
            let ring = Vec3::new(centroid.x, 0.0, centroid.z).normalize() * 2.0;
            let n = (t[1] - t[0]).cross(&(t[2] - t[0]));
            assert!(n.dot(&(centroid - ring)) > 0.0, "facet must face outward");
        }
    }

    #[test]
    fn mesh_triangulation_skips_invalid_and_degenerate_triangles() {
        let mesh = TriangleMesh::new(
            vec![
                Vec3::zeros(),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
                Vec3::new(2.0, 0.0, 0.0),
            ],
            vec![[0, 1, 2], [0, 1, 3], [0, 1, 9]],
        );
        // [0, 1, 3] is collinear (zero area); [0, 1, 9] indexes past the end.
        assert_eq!(triangulate(&mesh, ConcaveKind::Mesh).len(), 1);
    }

    #[test]
    fn triangle_proxy_support_inertia_and_ray_cast() {
        let tri = TriangleShape::new([
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]);
        assert_eq!(tri.support_point(&Vec3::x()), Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(tri.support_point(&Vec3::z()), Vec3::new(0.0, 0.0, 1.0));
        assert_eq!(tri.volume(), 0.0);
        let inertia = tri.inertia_tensor(6.0);
        assert!((inertia - inertia.transpose()).norm() < 1e-12);
        // The lamina lies in y = 0, so I_yy = I_xx + I_zz (perpendicular axis).
        assert!((inertia[(1, 1)] - inertia[(0, 0)] - inertia[(2, 2)]).abs() < 1e-12);
        let hit = tri
            .ray_cast(&Vec3::new(0.2, 1.0, 0.2), &Vec3::new(0.0, -1.0, 0.0), 10.0)
            .expect("ray hits the triangle");
        assert!((hit.toi - 1.0).abs() < 1e-12);
        assert!(hit.normal.y > 0.99, "normal faces the ray origin");
        assert!(
            tri.ray_cast(&Vec3::new(0.8, 1.0, 0.8), &Vec3::new(0.0, -1.0, 0.0), 10.0)
                .is_none(),
            "ray outside the hypotenuse misses"
        );
    }

    #[test]
    fn leaves_flatten_compounds_and_skip_sensors_and_orphans() {
        let mut bodies = RigidBodySet::new();
        let body = bodies.insert(RigidBody::new(1.0));
        let mut colliders = ColliderSet::new();
        let compound = Compound::new(vec![
            (
                Transform::from_position(Vec3::new(1.0, 0.0, 0.0)),
                Arc::new(Sphere::new(0.5)) as Arc<dyn Shape>,
            ),
            (
                Transform::from_position(Vec3::new(-1.0, 0.0, 0.0)),
                Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5))) as Arc<dyn Shape>,
            ),
        ]);
        colliders.insert(Collider::new(Arc::new(compound)).with_body(body));
        colliders.insert(
            Collider::new(Arc::new(Sphere::new(1.0)))
                .with_body(body)
                .as_sensor(),
        );
        colliders.insert(Collider::new(Arc::new(Sphere::new(1.0))));
        let leaves = gather_leaves(&bodies, &colliders);
        assert_eq!(leaves.len(), 2, "two compound children, no sensor/orphan");
        assert_eq!(leaves[0].kind, LeafKind::Convex(ShapeType::Sphere));
        assert_eq!(leaves[1].kind, LeafKind::Convex(ShapeType::Box));
        assert!((leaves[0].world.position.x - 1.0).abs() < 1e-12);
        assert!((leaves[1].world.position.x + 1.0).abs() < 1e-12);
    }

    #[test]
    fn concave_cache_builds_once_and_evicts_unused_shapes() {
        let mut bodies = RigidBodySet::new();
        let body = bodies.insert(RigidBody::new_static());
        let mut colliders = ColliderSet::new();
        let mesh = TriangleMesh::new(
            vec![
                Vec3::zeros(),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ],
            vec![[0, 1, 2]],
        );
        let handle = colliders.insert(Collider::new(Arc::new(mesh)).with_body(body));
        let mut cache = ConcaveCache::default();
        cache.prepare(&gather_leaves(&bodies, &colliders));
        assert_eq!(cache.len(), 1);
        cache.prepare(&gather_leaves(&bodies, &colliders));
        assert_eq!(cache.len(), 1, "rebuilt entries must not accumulate");
        colliders.remove(handle);
        cache.prepare(&gather_leaves(&bodies, &colliders));
        assert_eq!(cache.len(), 0, "entry evicted once the collider is gone");
    }
}
