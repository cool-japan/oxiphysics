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

    /// Construct a ray from normalised device coordinates (NDC) and a camera.
    ///
    /// * `ndc_x`, `ndc_y` — in `[-1, 1]` (right-handed NDC)
    /// * `inv_proj`        — inverse projection matrix (row-major 4×4)
    /// * `inv_view`        — inverse view matrix (row-major 4×4)
    pub fn from_screen(ndc_x: f32, ndc_y: f32, inv_proj: &[f32; 16], inv_view: &[f32; 16]) -> Self {
        // Unproject from clip to view space
        let clip = [ndc_x, ndc_y, -1.0_f32, 1.0_f32];
        let view = mat4_mul_vec4(inv_proj, clip);
        let view_dir = [view[0] / view[3], view[1] / view[3], -1.0, 0.0];
        // Rotate into world space (w=0 → direction only)
        let world_dir = mat4_mul_vec4(inv_view, view_dir);
        let origin = [inv_view[12], inv_view[13], inv_view[14]];
        Self::new(origin, [world_dir[0], world_dir[1], world_dir[2]])
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
}
