// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Interactive ray-picking and object selection for CPU-side rasterizers.
//!
//! This module provides [`Ray`], a world-space ray emitted from the camera via
//! mouse coordinates, and [`Picker`], which tests that ray against a scene of
//! pickable objects using Möller–Trumbore triangle intersection.
//!
//! ## Usage
//!
//! ```
//! use oxiphysics_viz::picking::{Picker, PickScene, PickableObject, Ray};
//!
//! // Build a scene with one triangle
//! let mut scene = PickScene::new();
//! scene.add_object(PickableObject {
//!     id:       0,
//!     vertices: vec![[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
//!     indices:  vec![0, 1, 2],
//! });
//!
//! let ray = Ray::new([0.0, 0.5, 5.0], [0.0, 0.0, -1.0]);
//! let picker = Picker::new();
//! if let Some(result) = picker.pick(&ray, &scene) {
//!     println!("Hit object {} at distance {:.3}", result.object_id, result.distance);
//! }
//! ```

// ── Ray ──────────────────────────────────────────────────────────────────────

/// A world-space ray with an origin and a unit-length direction.
///
/// Construct from a screen position with [`Ray::from_screen`] or directly
/// from world-space origin and direction with [`Ray::new`].
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin in world space.
    pub origin: [f32; 3],
    /// Ray direction (should be normalised).
    pub direction: [f32; 3],
}

impl Ray {
    /// Create a new ray.  `direction` is normalised internally.
    pub fn new(origin: [f32; 3], direction: [f32; 3]) -> Self {
        Self {
            origin,
            direction: normalize3(direction),
        }
    }

    /// Construct a ray through the pixel at normalised device coordinates
    /// (NDC) of a camera.
    ///
    /// * `ndc_x`, `ndc_y` — in `[-1, 1]`, x to the right and y up
    /// * `inv_proj`        — inverse projection matrix (column-major 4×4,
    ///   OpenGL clip space: NDC depth `-1` is the near plane, `+1` the far plane)
    /// * `inv_view`        — inverse view matrix, i.e. camera-to-world
    ///   (column-major 4×4, translation in elements 12..15)
    ///
    /// The points on the near and the far plane under the pixel are both
    /// unprojected through `inv_proj` (including the perspective divide) and
    /// carried to world space by `inv_view`; the ray starts at the near-plane
    /// point and heads towards the far-plane point.  This is exact for
    /// perspective, orthographic and infinite-far-plane projections alike.
    ///
    /// Because the origin lies on the near clip plane, distances reported by
    /// the [`Picker`] are measured from the near plane, and geometry in front
    /// of it (which is clipped on screen) cannot be picked.  Matrices that
    /// cannot be unprojected (for example a zero `inv_proj`) yield a ray with
    /// a zero direction at the camera position, which never hits anything.
    pub fn from_screen(ndc_x: f32, ndc_y: f32, inv_proj: &[f32; 16], inv_view: &[f32; 16]) -> Self {
        const MIN_W: f32 = 1e-9;
        let usable_w = |w: f32| w.is_finite() && w.abs() > MIN_W;
        let camera_position = [inv_view[12], inv_view[13], inv_view[14]];
        let no_ray = Self::new(camera_position, [0.0; 3]);

        // Homogeneous view-space points under the pixel on the near and far planes.
        let near_h = mat4_mul_vec4(inv_proj, [ndc_x, ndc_y, -1.0_f32, 1.0_f32]);
        let far_h = mat4_mul_vec4(inv_proj, [ndc_x, ndc_y, 1.0_f32, 1.0_f32]);
        if !usable_w(near_h[3]) {
            return no_ray;
        }
        let near_view = [
            near_h[0] / near_h[3],
            near_h[1] / near_h[3],
            near_h[2] / near_h[3],
        ];
        // The direction is formed in view space, where coordinates are small
        // and f32 keeps its precision however far the camera is from the origin.
        let dir_view = if usable_w(far_h[3]) {
            let far_view = [
                far_h[0] / far_h[3],
                far_h[1] / far_h[3],
                far_h[2] / far_h[3],
            ];
            sub3(far_view, near_view)
        } else {
            // Far plane at infinity: the homogeneous point with w = 0 is the
            // direction itself.  Orient it away from the eye, which sits at the
            // view-space origin of a perspective camera.
            let dir = [far_h[0], far_h[1], far_h[2]];
            if dot3(dir, near_view) < 0.0 {
                scale3(dir, -1.0)
            } else {
                dir
            }
        };
        // Carry the origin (w = 1) and the direction (w = 0) into world space.
        let origin = mat4_mul_vec4(inv_view, [near_view[0], near_view[1], near_view[2], 1.0]);
        let direction = mat4_mul_vec4(inv_view, [dir_view[0], dir_view[1], dir_view[2], 0.0]);
        let finite = origin[..3]
            .iter()
            .chain(direction[..3].iter())
            .all(|v| v.is_finite());
        if !finite {
            return no_ray;
        }
        Self::new(
            [origin[0], origin[1], origin[2]],
            [direction[0], direction[1], direction[2]],
        )
    }

    /// Evaluate a point on the ray: `origin + t * direction`.
    pub fn at(&self, t: f32) -> [f32; 3] {
        add3(self.origin, scale3(self.direction, t))
    }
}

// ── PickResult ───────────────────────────────────────────────────────────────

/// Result of a successful ray–object intersection.
#[derive(Debug, Clone)]
pub struct PickResult {
    /// ID of the hit object (from [`PickableObject::id`]).
    pub object_id: u32,
    /// Index of the first vertex of the hit triangle in `object.indices`.
    pub triangle_index: usize,
    /// Distance from the ray origin to the hit point.
    pub distance: f32,
    /// Barycentric coordinates (u, v) within the hit triangle.
    ///
    /// The third barycentric coordinate is `w = 1 - u - v`.
    pub barycentric: [f32; 2],
    /// World-space position of the hit point.
    pub hit_point: [f32; 3],
    /// Geometric normal of the hit triangle (unnormalised).
    pub normal: [f32; 3],
}

// ── PickableObject ───────────────────────────────────────────────────────────

/// A named triangle mesh that can be tested by the [`Picker`].
#[derive(Debug, Clone)]
pub struct PickableObject {
    /// Unique identifier.
    pub id: u32,
    /// World-space vertex positions.
    pub vertices: Vec<[f32; 3]>,
    /// Triangle indices (groups of 3).
    pub indices: Vec<u32>,
}

impl PickableObject {
    /// Create a single-triangle object.
    pub fn single_triangle(id: u32, a: [f32; 3], b: [f32; 3], c: [f32; 3]) -> Self {
        Self {
            id,
            vertices: vec![a, b, c],
            indices: vec![0, 1, 2],
        }
    }

    /// Return the axis-aligned bounding box of this object.
    pub fn aabb(&self) -> ([f32; 3], [f32; 3]) {
        let mut lo = [f32::INFINITY; 3];
        let mut hi = [f32::NEG_INFINITY; 3];
        for v in &self.vertices {
            for k in 0..3 {
                lo[k] = lo[k].min(v[k]);
                hi[k] = hi[k].max(v[k]);
            }
        }
        (lo, hi)
    }
}

// ── PickScene ────────────────────────────────────────────────────────────────

/// A collection of [`PickableObject`]s that form the pickable scene.
#[derive(Debug, Default)]
pub struct PickScene {
    /// All registered objects.
    pub objects: Vec<PickableObject>,
}

impl PickScene {
    /// Create an empty scene.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an object to the scene.
    pub fn add_object(&mut self, obj: PickableObject) {
        self.objects.push(obj);
    }

    /// Remove an object by ID.  Returns `true` if found.
    pub fn remove_object(&mut self, id: u32) -> bool {
        if let Some(pos) = self.objects.iter().position(|o| o.id == id) {
            self.objects.remove(pos);
            true
        } else {
            false
        }
    }
}

// ── Picker ───────────────────────────────────────────────────────────────────

/// Ray-casting picker that tests a [`Ray`] against all objects in a [`PickScene`].
///
/// Uses Möller–Trumbore triangle intersection with front-face culling disabled
/// (both sides are pickable).
#[derive(Debug, Default)]
pub struct Picker {
    /// Minimum distance from the ray origin to consider a hit.
    pub min_distance: f32,
    /// Maximum distance from the ray origin to consider a hit.
    pub max_distance: f32,
    /// Whether to cull back-facing triangles.
    pub backface_cull: bool,
}

impl Picker {
    const EPSILON: f32 = 1e-7;

    /// Create a picker with default parameters.
    pub fn new() -> Self {
        Self {
            min_distance: 0.001,
            max_distance: 10_000.0,
            backface_cull: false,
        }
    }

    /// Find the **nearest** object that the ray hits.
    ///
    /// Returns `None` if no object is hit within `[min_distance, max_distance]`.
    pub fn pick(&self, ray: &Ray, scene: &PickScene) -> Option<PickResult> {
        let mut best: Option<PickResult> = None;
        for obj in &scene.objects {
            if let Some(r) = self.pick_object(ray, obj) {
                let closer = best.as_ref().is_none_or(|b| r.distance < b.distance);
                if closer {
                    best = Some(r);
                }
            }
        }
        best
    }

    /// Find **all** objects that the ray hits, sorted by ascending distance.
    pub fn pick_all(&self, ray: &Ray, scene: &PickScene) -> Vec<PickResult> {
        let mut results: Vec<PickResult> = scene
            .objects
            .iter()
            .filter_map(|obj| self.pick_object(ray, obj))
            .collect();
        results.sort_by(|a, b| {
            a.distance
                .partial_cmp(&b.distance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        results
    }

    /// Test a ray against a single object.
    pub fn pick_object(&self, ray: &Ray, obj: &PickableObject) -> Option<PickResult> {
        // AABB pre-test for culling
        let (lo, hi) = obj.aabb();
        if !ray_aabb(ray, lo, hi) {
            return None;
        }

        let tri_count = obj.indices.len() / 3;
        let mut best: Option<PickResult> = None;

        for tri in 0..tri_count {
            let i0 = obj.indices[tri * 3] as usize;
            let i1 = obj.indices[tri * 3 + 1] as usize;
            let i2 = obj.indices[tri * 3 + 2] as usize;
            let v0 = obj.vertices[i0];
            let v1 = obj.vertices[i1];
            let v2 = obj.vertices[i2];

            if let Some((t, u, v)) = moller_trumbore(ray, v0, v1, v2, self.backface_cull) {
                if t < self.min_distance || t > self.max_distance {
                    continue;
                }
                let hit = ray.at(t);
                let normal = cross3(sub3(v1, v0), sub3(v2, v0));
                let closer = best.as_ref().is_none_or(|b| t < b.distance);
                if closer {
                    best = Some(PickResult {
                        object_id: obj.id,
                        triangle_index: tri,
                        distance: t,
                        barycentric: [u, v],
                        hit_point: hit,
                        normal,
                    });
                }
            }
        }
        best
    }
}

// ── SelectionSet ─────────────────────────────────────────────────────────────

/// Tracks which objects are currently selected.
///
/// Supports single selection, multi-selection (additive), and box selection.
#[derive(Debug, Default)]
pub struct SelectionSet {
    /// IDs of currently selected objects.
    pub selected: std::collections::HashSet<u32>,
}

impl SelectionSet {
    /// Create an empty selection.
    pub fn new() -> Self {
        Self::default()
    }

    /// Select exactly one object, clearing all previous selections.
    pub fn select_single(&mut self, id: u32) {
        self.selected.clear();
        self.selected.insert(id);
    }

    /// Toggle the selection state of an object.
    pub fn toggle(&mut self, id: u32) {
        if self.selected.contains(&id) {
            self.selected.remove(&id);
        } else {
            self.selected.insert(id);
        }
    }

    /// Add objects inside an AABB to the selection.
    pub fn select_in_aabb(
        &mut self,
        scene: &PickScene,
        lo: [f32; 3],
        hi: [f32; 3],
        additive: bool,
    ) {
        if !additive {
            self.selected.clear();
        }
        for obj in &scene.objects {
            let (olo, ohi) = obj.aabb();
            // Object's AABB overlaps the selection box?
            let overlap = (0..3).all(|k| olo[k] <= hi[k] && ohi[k] >= lo[k]);
            if overlap {
                self.selected.insert(obj.id);
            }
        }
    }

    /// Return `true` if `id` is selected.
    pub fn is_selected(&self, id: u32) -> bool {
        self.selected.contains(&id)
    }

    /// Clear all selections.
    pub fn clear(&mut self) {
        self.selected.clear();
    }
}

// ── Möller–Trumbore intersection ─────────────────────────────────────────────

/// Möller–Trumbore ray–triangle intersection.
///
/// Returns `Some((t, u, v))` on hit where `t` is the ray parameter and
/// `(u, v)` are barycentric coordinates, or `None` on miss.
pub fn moller_trumbore(
    ray: &Ray,
    v0: [f32; 3],
    v1: [f32; 3],
    v2: [f32; 3],
    backface_cull: bool,
) -> Option<(f32, f32, f32)> {
    let edge1 = sub3(v1, v0);
    let edge2 = sub3(v2, v0);
    let h = cross3(ray.direction, edge2);
    let a = dot3(edge1, h);

    if backface_cull && a < Picker::EPSILON {
        return None;
    }
    if a.abs() < Picker::EPSILON {
        return None;
    }

    let f = 1.0 / a;
    let s = sub3(ray.origin, v0);
    let u = f * dot3(s, h);
    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = cross3(s, edge1);
    let v = f * dot3(ray.direction, q);
    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * dot3(edge2, q);
    if t < 0.0 { None } else { Some((t, u, v)) }
}

// ── Ray–AABB intersection ─────────────────────────────────────────────────────

/// Ray–AABB slab test. Returns `true` if the ray intersects the box.
fn ray_aabb(ray: &Ray, lo: [f32; 3], hi: [f32; 3]) -> bool {
    let mut t_min = 0.0_f32;
    let mut t_max = f32::INFINITY;
    for k in 0..3 {
        let inv_d = 1.0 / ray.direction[k];
        let t0 = (lo[k] - ray.origin[k]) * inv_d;
        let t1 = (hi[k] - ray.origin[k]) * inv_d;
        let (t_near, t_far) = if inv_d >= 0.0 { (t0, t1) } else { (t1, t0) };
        t_min = t_min.max(t_near);
        t_max = t_max.min(t_far);
        if t_max < t_min {
            return false;
        }
    }
    true
}

// ── Vector helpers ────────────────────────────────────────────────────────────

#[inline(always)]
fn dot3(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline(always)]
fn cross3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline(always)]
fn sub3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline(always)]
fn add3(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline(always)]
fn scale3(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline(always)]
fn normalize3(a: [f32; 3]) -> [f32; 3] {
    let len = dot3(a, a).sqrt();
    if len < 1e-12 {
        [0.0; 3]
    } else {
        [a[0] / len, a[1] / len, a[2] / len]
    }
}

fn mat4_mul_vec4(m: &[f32; 16], v: [f32; 4]) -> [f32; 4] {
    [
        m[0] * v[0] + m[4] * v[1] + m[8] * v[2] + m[12] * v[3],
        m[1] * v[0] + m[5] * v[1] + m[9] * v[2] + m[13] * v[3],
        m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3],
        m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3],
    ]
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn unit_triangle() -> PickableObject {
        PickableObject::single_triangle(1, [-1.0, -1.0, 0.0], [1.0, -1.0, 0.0], [0.0, 1.0, 0.0])
    }

    #[test]
    fn ray_hits_triangle_through_centre() {
        let mut scene = PickScene::new();
        scene.add_object(unit_triangle());
        let ray = Ray::new([0.0, 0.0, 5.0], [0.0, 0.0, -1.0]);
        let result = Picker::new().pick(&ray, &scene);
        assert!(result.is_some(), "ray should hit triangle");
        let r = result.unwrap();
        assert_eq!(r.object_id, 1);
        assert!((r.hit_point[2]).abs() < 0.01, "hit z should be ~0");
    }

    #[test]
    fn ray_misses_triangle_offset() {
        let mut scene = PickScene::new();
        scene.add_object(unit_triangle());
        let ray = Ray::new([5.0, 5.0, 5.0], [0.0, 0.0, -1.0]);
        assert!(Picker::new().pick(&ray, &scene).is_none());
    }

    #[test]
    fn pick_returns_nearest_hit() {
        let mut scene = PickScene::new();
        scene.add_object(PickableObject::single_triangle(
            1,
            [-1.0, -1.0, 2.0],
            [1.0, -1.0, 2.0],
            [0.0, 1.0, 2.0],
        ));
        scene.add_object(PickableObject::single_triangle(
            2,
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [0.0, 1.0, 1.0],
        ));
        let ray = Ray::new([0.0, 0.0, 5.0], [0.0, 0.0, -1.0]);
        let r = Picker::new().pick(&ray, &scene).unwrap();
        assert_eq!(r.object_id, 1, "closer object at z=2 should be hit first");
    }

    #[test]
    fn pick_all_sorted_ascending() {
        let mut scene = PickScene::new();
        scene.add_object(PickableObject::single_triangle(
            1,
            [-1.0, -1.0, 2.0],
            [1.0, -1.0, 2.0],
            [0.0, 1.0, 2.0],
        ));
        scene.add_object(PickableObject::single_triangle(
            2,
            [-1.0, -1.0, 1.0],
            [1.0, -1.0, 1.0],
            [0.0, 1.0, 1.0],
        ));
        let ray = Ray::new([0.0, 0.0, 5.0], [0.0, 0.0, -1.0]);
        let all = Picker::new().pick_all(&ray, &scene);
        assert_eq!(all.len(), 2);
        assert!(
            all[0].distance < all[1].distance,
            "results should be sorted"
        );
    }

    #[test]
    fn selection_set_toggle() {
        let mut sel = SelectionSet::new();
        sel.toggle(5);
        assert!(sel.is_selected(5));
        sel.toggle(5);
        assert!(!sel.is_selected(5));
    }

    #[test]
    fn selection_set_aabb() {
        let mut scene = PickScene::new();
        scene.add_object(PickableObject::single_triangle(
            10,
            [-0.5, -0.5, -0.5],
            [0.5, -0.5, -0.5],
            [0.0, 0.5, -0.5],
        ));
        let mut sel = SelectionSet::new();
        sel.select_in_aabb(&scene, [-2.0; 3], [2.0; 3], false);
        assert!(sel.is_selected(10));
    }

    // ── Ray::from_screen ──────────────────────────────────────────────────────

    const IDENTITY: [f32; 16] = [
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ];

    /// Column-major product `a * b`.
    fn mat4_mul(a: &[f32; 16], b: &[f32; 16]) -> [f32; 16] {
        let mut out = [0.0_f32; 16];
        for col in 0..4 {
            for row in 0..4 {
                out[col * 4 + row] = (0..4).map(|k| a[k * 4 + row] * b[col * 4 + k]).sum();
            }
        }
        out
    }

    fn assert_identity(m: &[f32; 16], tol: f32) {
        for col in 0..4 {
            for row in 0..4 {
                let expected = if col == row { 1.0 } else { 0.0 };
                assert!(
                    (m[col * 4 + row] - expected).abs() < tol,
                    "element ({row},{col}) = {}, expected {expected}",
                    m[col * 4 + row]
                );
            }
        }
    }

    /// OpenGL perspective projection and its analytic inverse (column-major).
    fn gl_perspective(fov_y: f64, aspect: f64, near: f64, far: f64) -> ([f32; 16], [f32; 16]) {
        let f = 1.0 / (0.5 * fov_y).tan();
        let mut p = [0.0_f64; 16];
        p[0] = f / aspect;
        p[5] = f;
        p[10] = (far + near) / (near - far);
        p[11] = -1.0;
        p[14] = 2.0 * far * near / (near - far);
        let mut inv = [0.0_f64; 16];
        inv[0] = aspect / f;
        inv[5] = 1.0 / f;
        inv[11] = (near - far) / (2.0 * far * near);
        inv[14] = -1.0;
        inv[15] = (far + near) / (2.0 * far * near);
        (p.map(|v| v as f32), inv.map(|v| v as f32))
    }

    /// OpenGL perspective with the far plane at infinity, and its inverse.
    fn gl_infinite_perspective(fov_y: f64, aspect: f64, near: f64) -> ([f32; 16], [f32; 16]) {
        let f = 1.0 / (0.5 * fov_y).tan();
        let mut p = [0.0_f64; 16];
        p[0] = f / aspect;
        p[5] = f;
        p[10] = -1.0;
        p[11] = -1.0;
        p[14] = -2.0 * near;
        let mut inv = [0.0_f64; 16];
        inv[0] = aspect / f;
        inv[5] = 1.0 / f;
        inv[11] = -1.0 / (2.0 * near);
        inv[14] = -1.0;
        inv[15] = 1.0 / (2.0 * near);
        (p.map(|v| v as f32), inv.map(|v| v as f32))
    }

    /// OpenGL orthographic projection and its analytic inverse.
    fn gl_orthographic(bounds: [f64; 4], near: f64, far: f64) -> ([f32; 16], [f32; 16]) {
        let [l, r, b, t] = bounds;
        let mut p = [0.0_f64; 16];
        p[0] = 2.0 / (r - l);
        p[5] = 2.0 / (t - b);
        p[10] = -2.0 / (far - near);
        p[12] = -(r + l) / (r - l);
        p[13] = -(t + b) / (t - b);
        p[14] = -(far + near) / (far - near);
        p[15] = 1.0;
        let mut inv = [0.0_f64; 16];
        inv[0] = (r - l) / 2.0;
        inv[5] = (t - b) / 2.0;
        inv[10] = -(far - near) / 2.0;
        inv[12] = (r + l) / 2.0;
        inv[13] = (t + b) / 2.0;
        inv[14] = -(far + near) / 2.0;
        inv[15] = 1.0;
        (p.map(|v| v as f32), inv.map(|v| v as f32))
    }

    /// Inverse of a rigid transform `[R | t]` (column-major).
    fn rigid_inverse(m: &[f32; 16]) -> [f32; 16] {
        let mut out = [0.0_f32; 16];
        for col in 0..3 {
            for row in 0..3 {
                out[col * 4 + row] = m[row * 4 + col];
            }
        }
        for row in 0..3 {
            out[12 + row] = -(out[row] * m[12] + out[4 + row] * m[13] + out[8 + row] * m[14]);
        }
        out[15] = 1.0;
        out
    }

    /// Camera-to-world matrix of a camera at `position` turned by 90 degrees
    /// about +y, so that it looks along world -x.
    fn camera_looking_down_minus_x(position: [f32; 3]) -> [f32; 16] {
        [
            0.0,
            0.0,
            -1.0,
            0.0, // image of +x
            0.0,
            1.0,
            0.0,
            0.0, // image of +y
            1.0,
            0.0,
            0.0,
            0.0, // image of +z
            position[0],
            position[1],
            position[2],
            1.0,
        ]
    }

    fn assert_vec3_close(actual: [f32; 3], expected: [f32; 3], tol: f32, what: &str) {
        for k in 0..3 {
            assert!(
                (actual[k] - expected[k]).abs() < tol,
                "{what}: {actual:?} vs expected {expected:?}"
            );
        }
    }

    #[test]
    fn test_projection_helpers_are_inverses() {
        let (p, inv) = gl_perspective(std::f64::consts::FRAC_PI_2, 1.5, 0.1, 100.0);
        assert_identity(&mat4_mul(&p, &inv), 1e-4);
        assert_identity(&mat4_mul(&inv, &p), 1e-4);
        let (p, inv) = gl_infinite_perspective(1.0, 1.0, 0.1);
        assert_identity(&mat4_mul(&p, &inv), 1e-5);
        let (p, inv) = gl_orthographic([-2.0, 2.0, -1.0, 1.0], 0.1, 10.0);
        assert_identity(&mat4_mul(&p, &inv), 1e-5);
        let cam = camera_looking_down_minus_x([3.0, 2.0, 5.0]);
        assert_identity(&mat4_mul(&cam, &rigid_inverse(&cam)), 1e-6);
    }

    #[test]
    fn test_from_screen_perspective_edge_rays_use_the_true_field_of_view() {
        // fov 90 degrees, aspect 1, near 0.1: the right edge of the screen is
        // a 45 degree ray.  Forcing view z = -1 (instead of z = -near) used to
        // give atan(0.1) = 5.7 degrees.
        let (_, inv_proj) = gl_perspective(std::f64::consts::FRAC_PI_2, 1.0, 0.1, 100.0);
        let right = Ray::from_screen(1.0, 0.0, &inv_proj, &IDENTITY);
        let s = std::f32::consts::FRAC_1_SQRT_2;
        assert_vec3_close(right.direction, [s, 0.0, -s], 1e-4, "right edge");
        let centre = Ray::from_screen(0.0, 0.0, &inv_proj, &IDENTITY);
        assert_vec3_close(centre.direction, [0.0, 0.0, -1.0], 1e-5, "centre");
        let corner = Ray::from_screen(1.0, 1.0, &inv_proj, &IDENTITY);
        let c = 1.0 / 3.0_f32.sqrt();
        assert_vec3_close(corner.direction, [c, c, -c], 1e-4, "corner");
        // The origin is the near-plane point under the pixel.
        assert_vec3_close(right.origin, [0.1, 0.0, -0.1], 1e-5, "near-plane origin");
    }

    #[test]
    fn test_from_screen_follows_the_camera_pose() {
        let (_, inv_proj) = gl_perspective(std::f64::consts::FRAC_PI_2, 1.0, 0.1, 100.0);
        let inv_view = camera_looking_down_minus_x([3.0, 2.0, 5.0]);
        let centre = Ray::from_screen(0.0, 0.0, &inv_proj, &inv_view);
        assert_vec3_close(centre.direction, [-1.0, 0.0, 0.0], 1e-5, "centre direction");
        assert_vec3_close(centre.origin, [2.9, 2.0, 5.0], 1e-5, "centre origin");
        // Screen-right is the camera's +x, which points along world -z.
        let right = Ray::from_screen(1.0, 0.0, &inv_proj, &inv_view);
        let s = std::f32::consts::FRAC_1_SQRT_2;
        assert_vec3_close(right.direction, [-s, 0.0, -s], 1e-4, "right direction");
    }

    #[test]
    fn test_from_screen_ray_passes_through_the_projected_point() {
        // Project a world point to NDC with the real matrices, then unproject
        // the pixel: the ray must pass through the original point.
        let (proj, inv_proj) = gl_perspective(1.1, 1.6, 0.25, 200.0);
        let inv_view = camera_looking_down_minus_x([3.0, 2.0, 5.0]);
        let view = rigid_inverse(&inv_view);
        let view_proj = mat4_mul(&proj, &view);
        for world in [[-4.0_f32, 2.5, 4.0], [-9.0, 0.0, 7.5], [-20.0, 5.0, 1.0]] {
            let clip = mat4_mul_vec4(&view_proj, [world[0], world[1], world[2], 1.0]);
            let ndc = [clip[0] / clip[3], clip[1] / clip[3]];
            assert!(
                ndc[0].abs() < 1.0 && ndc[1].abs() < 1.0,
                "{world:?} must be on screen"
            );
            let ray = Ray::from_screen(ndc[0], ndc[1], &inv_proj, &inv_view);
            let to_point = sub3(world, ray.origin);
            let along = dot3(to_point, ray.direction);
            let closest = add3(ray.origin, scale3(ray.direction, along));
            let miss = sub3(world, closest);
            let miss_len = dot3(miss, miss).sqrt();
            assert!(along > 0.0, "{world:?} must lie ahead of the ray origin");
            assert!(miss_len < 1e-3, "ray misses {world:?} by {miss_len}");
        }
    }

    #[test]
    fn test_from_screen_orthographic_rays_are_parallel() {
        let (_, inv_proj) = gl_orthographic([-2.0, 2.0, -1.0, 1.0], 0.1, 10.0);
        let a = Ray::from_screen(0.5, 0.25, &inv_proj, &IDENTITY);
        let b = Ray::from_screen(-1.0, -1.0, &inv_proj, &IDENTITY);
        assert_vec3_close(a.direction, [0.0, 0.0, -1.0], 1e-5, "ortho direction a");
        assert_vec3_close(b.direction, [0.0, 0.0, -1.0], 1e-5, "ortho direction b");
        assert_vec3_close(a.origin, [1.0, 0.25, -0.1], 1e-5, "ortho origin a");
        assert_vec3_close(b.origin, [-2.0, -1.0, -0.1], 1e-5, "ortho origin b");
    }

    #[test]
    fn test_from_screen_infinite_far_plane() {
        let (_, inv_proj) = gl_infinite_perspective(std::f64::consts::FRAC_PI_2, 1.0, 0.1);
        let right = Ray::from_screen(1.0, 0.0, &inv_proj, &IDENTITY);
        let s = std::f32::consts::FRAC_1_SQRT_2;
        assert_vec3_close(right.direction, [s, 0.0, -s], 1e-4, "infinite right edge");
        assert_vec3_close(right.origin, [0.1, 0.0, -0.1], 1e-5, "infinite origin");
    }

    #[test]
    fn test_from_screen_picks_the_pixel_under_the_cursor() {
        // A large triangle in the plane z = -2 in front of an identity camera.
        let mut scene = PickScene::new();
        scene.add_object(PickableObject::single_triangle(
            7,
            [-10.0, -10.0, -2.0],
            [10.0, -10.0, -2.0],
            [0.0, 10.0, -2.0],
        ));
        let (proj, inv_proj) = gl_perspective(std::f64::consts::FRAC_PI_2, 1.0, 0.1, 100.0);
        let target = [0.8_f32, 0.3, -2.0];
        let clip = mat4_mul_vec4(&proj, [target[0], target[1], target[2], 1.0]);
        let ray = Ray::from_screen(clip[0] / clip[3], clip[1] / clip[3], &inv_proj, &IDENTITY);
        let hit = Picker::new()
            .pick(&ray, &scene)
            .expect("pixel must hit the triangle");
        assert_eq!(hit.object_id, 7);
        assert_vec3_close(hit.hit_point, target, 1e-3, "picked point");
        // Distance is measured from the near-plane origin (0.04, 0.015, -0.1)
        // to the hit point (0.8, 0.3, -2).
        let expected = (0.76_f32 * 0.76 + 0.285 * 0.285 + 1.9 * 1.9).sqrt();
        assert!(
            (hit.distance - expected).abs() < 1e-3,
            "distance {}",
            hit.distance
        );
    }

    #[test]
    fn test_from_screen_degenerate_matrices_give_a_harmless_ray() {
        let zero = [0.0_f32; 16];
        let ray = Ray::from_screen(0.3, -0.2, &zero, &IDENTITY);
        assert_eq!(ray.direction, [0.0; 3]);
        assert_eq!(ray.origin, [0.0; 3]);
        let nan_proj = [f32::NAN; 16];
        let ray = Ray::from_screen(0.0, 0.0, &nan_proj, &IDENTITY);
        assert_eq!(ray.direction, [0.0; 3]);
        let mut scene = PickScene::new();
        scene.add_object(unit_triangle());
        assert!(Picker::new().pick(&ray, &scene).is_none());
    }
}
