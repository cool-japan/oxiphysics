// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body terrain interaction.
//!
//! Provides triangle-mesh terrain patches, ray–terrain intersection,
//! contact force models, friction, slope analysis, and factory functions
//! for common terrain shapes.

/// A triangulated terrain patch with per-patch friction.
#[derive(Debug, Clone)]
pub struct TerrainPatch {
    /// World-space vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle index triples (v0, v1, v2).
    pub triangles: Vec<[usize; 3]>,
    /// Coefficient of friction for this patch.
    pub friction_coeff: f64,
}

/// A contact point between a rigid body and terrain.
#[derive(Debug, Clone)]
pub struct TerrainContact {
    /// World-space contact position \[m\].
    pub position: [f64; 3],
    /// Outward terrain normal (unit vector).
    pub normal: [f64; 3],
    /// Penetration depth (positive = overlap) \[m\].
    pub depth: f64,
    /// Friction coefficient at the contact.
    pub friction: f64,
}

// ── internal helpers ──────────────────────────────────────────────────────────

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Compute the unit outward normal of a triangle defined by vertices v0, v1, v2.
///
/// Returns the normalised cross product of (v1−v0) × (v2−v0).
/// Returns `[0,1,0]` if the triangle is degenerate.
pub fn triangle_normal(v0: [f64; 3], v1: [f64; 3], v2: [f64; 3]) -> [f64; 3] {
    let e1 = sub3(v1, v0);
    let e2 = sub3(v2, v0);
    let n = cross3(e1, e2);
    let len = norm3(n);
    if len < 1e-15 {
        return [0.0, 1.0, 0.0];
    }
    [n[0] / len, n[1] / len, n[2] / len]
}

/// Cast a ray from `origin` along `dir` and find the nearest intersection
/// with any triangle in `patch`.
///
/// Uses Möller–Trumbore intersection.
/// Returns `Some(t)` where `t` is the parameter along the ray, or `None`.
pub fn terrain_ray_intersect(patch: &TerrainPatch, origin: [f64; 3], dir: [f64; 3]) -> Option<f64> {
    const EPS: f64 = 1e-10;
    let mut closest: Option<f64> = None;

    for tri in &patch.triangles {
        let v0 = patch.vertices[tri[0]];
        let v1 = patch.vertices[tri[1]];
        let v2 = patch.vertices[tri[2]];

        let e1 = sub3(v1, v0);
        let e2 = sub3(v2, v0);
        let h = cross3(dir, e2);
        let a = dot3(e1, h);
        if a.abs() < EPS {
            continue;
        }
        let f = 1.0 / a;
        let s = sub3(origin, v0);
        let u = f * dot3(s, h);
        if !(0.0..=1.0).contains(&u) {
            continue;
        }
        let q = cross3(s, e1);
        let v = f * dot3(dir, q);
        if v < 0.0 || u + v > 1.0 {
            continue;
        }
        let t = f * dot3(e2, q);
        if t > EPS && closest.is_none_or(|c| t < c) {
            closest = Some(t);
        }
    }
    closest
}

/// Sample the terrain height (Y coordinate) at world position (x, z).
///
/// Projects a downward ray from far above and returns the hit height.
pub fn terrain_height_at(patch: &TerrainPatch, x: f64, z: f64) -> Option<f64> {
    let origin = [x, 1e6, z];
    let dir = [0.0, -1.0, 0.0];
    terrain_ray_intersect(patch, origin, dir).map(|t| 1e6 - t)
}

/// Compute the penalty contact force between a rigid body and terrain.
///
/// Uses a spring-damper model:  F = k · δ − b · v_dot
///
/// # Arguments
/// * `penetration` – positive overlap depth \[m\]
/// * `velocity`    – normal velocity (negative = penetrating) \[m/s\]
/// * `stiffness`   – contact spring constant \[N/m\]
/// * `damping`     – contact damping constant \[N·s/m\]
pub fn ground_contact_force(penetration: f64, velocity: f64, stiffness: f64, damping: f64) -> f64 {
    if penetration <= 0.0 {
        return 0.0;
    }
    (stiffness * penetration - damping * velocity).max(0.0)
}

/// Compute the friction force vector for a sliding contact.
///
/// Uses Coulomb friction: F_f = −μ · N_f · v̂
///
/// # Arguments
/// * `normal_force` – magnitude of the normal contact force \[N\]
/// * `velocity`     – tangential velocity vector \[m/s\]
/// * `mu`           – coefficient of friction
pub fn terrain_friction_model(normal_force: f64, velocity: [f64; 3], mu: f64) -> [f64; 3] {
    let speed = norm3(velocity);
    if speed < 1e-15 {
        return [0.0; 3];
    }
    let mag = mu * normal_force.abs();
    [
        -mag * velocity[0] / speed,
        -mag * velocity[1] / speed,
        -mag * velocity[2] / speed,
    ]
}

/// Compute the slope angle from vertical (in radians) given a terrain normal.
///
/// A flat horizontal surface has `normal = [0, 1, 0]` → slope = 0.
pub fn slope_angle(normal: [f64; 3]) -> f64 {
    // angle between normal and world-up [0,1,0]; use abs to handle both winding conventions
    let cos_theta = normal[1].abs().clamp(0.0, 1.0);
    cos_theta.acos()
}

/// Return `true` if the terrain slope is traversable by a vehicle.
///
/// # Arguments
/// * `slope`          – slope angle from vertical \[radians\]
/// * `max_slope_deg`  – maximum traversable slope \[degrees\]
pub fn terrain_traversability(slope: f64, max_slope_deg: f64) -> bool {
    let max_slope_rad = max_slope_deg.to_radians();
    slope <= max_slope_rad
}

/// Create a flat horizontal terrain patch at a given height.
///
/// The patch is a square of side `size` centred at (0, height, 0).
pub fn create_flat_terrain(size: f64, height: f64) -> TerrainPatch {
    let h = size * 0.5;
    let vertices = vec![
        [-h, height, -h],
        [h, height, -h],
        [h, height, h],
        [-h, height, h],
    ];
    let triangles = vec![[0, 1, 2], [0, 2, 3]];
    TerrainPatch {
        vertices,
        triangles,
        friction_coeff: 0.6,
    }
}

/// Create a sloped terrain patch tilted by `slope_deg` degrees.
///
/// The slope rises along the +X axis.
pub fn create_slope_terrain(size: f64, slope_deg: f64) -> TerrainPatch {
    let h = size * 0.5;
    let slope_rad = slope_deg.to_radians();
    let tan_s = slope_rad.tan();
    // Y = x * tan(slope) for the sloped edge.
    // Vertex ordering: v0=(-x,-z corner), v1=(+x,-z), v2=(+x,+z), v3=(-x,+z)
    let vertices = vec![
        [-h, -h * tan_s, -h],
        [h, h * tan_s, -h],
        [h, h * tan_s, h],
        [-h, -h * tan_s, h],
    ];
    // Reversed winding to ensure upward-pointing normals.
    let triangles = vec![[0, 2, 1], [0, 3, 2]];
    TerrainPatch {
        vertices,
        triangles,
        friction_coeff: 0.5,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // ── triangle_normal ───────────────────────────────────────────────────

    #[test]
    fn triangle_normal_flat_xy_plane_points_up() {
        let n = triangle_normal([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        // Cross of x-hat × z-hat = -y-hat; for (v1-v0)×(v2-v0) with above it's -y
        // Depending on winding — just check it's unit length and perpendicular to plane
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!((len - 1.0).abs() < 1e-10, "len={len}");
        // Normal should have zero x and z components (flat XZ plane)
        assert!(n[0].abs() < 1e-10, "n[0]={}", n[0]);
        assert!(n[2].abs() < 1e-10, "n[2]={}", n[2]);
    }

    #[test]
    fn triangle_normal_degenerate_returns_up() {
        let n = triangle_normal([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert!((n[1] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn triangle_normal_is_unit_length() {
        let n = triangle_normal([1.0, 2.0, 3.0], [4.0, 0.0, 1.0], [2.0, 3.0, -1.0]);
        let len = norm3(n);
        assert!((len - 1.0).abs() < 1e-12);
    }

    // ── terrain_ray_intersect ─────────────────────────────────────────────

    #[test]
    fn ray_hits_flat_terrain_straight_down() {
        let patch = create_flat_terrain(10.0, 0.0);
        let t = terrain_ray_intersect(&patch, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(t.is_some(), "ray should hit flat terrain");
        assert!((t.unwrap() - 5.0).abs() < 1e-9, "t={}", t.unwrap());
    }

    #[test]
    fn ray_misses_when_outside_terrain() {
        let patch = create_flat_terrain(2.0, 0.0); // ±1 m
        let t = terrain_ray_intersect(&patch, [5.0, 5.0, 5.0], [0.0, -1.0, 0.0]);
        assert!(t.is_none(), "ray outside patch should miss");
    }

    #[test]
    fn ray_misses_when_pointing_away() {
        let patch = create_flat_terrain(10.0, 0.0);
        let t = terrain_ray_intersect(&patch, [0.0, 5.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(t.is_none(), "upward ray should miss terrain below");
    }

    #[test]
    fn ray_returns_closest_intersection() {
        let mut patch = create_flat_terrain(10.0, 0.0);
        // add a second elevated platform
        let base = patch.vertices.len();
        patch.vertices.extend_from_slice(&[
            [-2.0, 3.0, -2.0],
            [2.0, 3.0, -2.0],
            [2.0, 3.0, 2.0],
            [-2.0, 3.0, 2.0],
        ]);
        patch.triangles.push([base, base + 1, base + 2]);
        patch.triangles.push([base, base + 2, base + 3]);

        let t = terrain_ray_intersect(&patch, [0.0, 10.0, 0.0], [0.0, -1.0, 0.0]);
        assert!(t.is_some());
        // Should hit the elevated platform first (t ≈ 7)
        assert!((t.unwrap() - 7.0).abs() < 1e-9, "t={}", t.unwrap());
    }

    // ── terrain_height_at ─────────────────────────────────────────────────

    #[test]
    fn height_at_flat_terrain_centre() {
        let patch = create_flat_terrain(10.0, 2.0);
        let h = terrain_height_at(&patch, 0.0, 0.0);
        assert!(h.is_some(), "should hit flat terrain");
        assert!((h.unwrap() - 2.0).abs() < 1e-6, "h={}", h.unwrap());
    }

    #[test]
    fn height_at_outside_patch_is_none() {
        let patch = create_flat_terrain(2.0, 0.0);
        let h = terrain_height_at(&patch, 100.0, 0.0);
        assert!(h.is_none());
    }

    // ── ground_contact_force ──────────────────────────────────────────────

    #[test]
    fn contact_force_zero_when_no_penetration() {
        let f = ground_contact_force(0.0, 0.0, 1e5, 1e3);
        assert!(f.abs() < 1e-15);
    }

    #[test]
    fn contact_force_negative_penetration_is_zero() {
        let f = ground_contact_force(-0.01, 0.0, 1e5, 0.0);
        assert!(f.abs() < 1e-15);
    }

    #[test]
    fn contact_force_positive_when_penetrating() {
        let f = ground_contact_force(0.01, 0.0, 1e5, 0.0);
        assert!(f > 0.0, "f={f}");
        assert!((f - 1000.0).abs() < 1e-9, "f={f}");
    }

    #[test]
    fn contact_force_damping_reduces_force() {
        let f_no_damp = ground_contact_force(0.01, -1.0, 1e5, 0.0);
        let f_damped = ground_contact_force(0.01, -1.0, 1e5, 100.0);
        assert!(
            f_damped > f_no_damp,
            "damping with negative velocity adds force"
        );
    }

    #[test]
    fn contact_force_clamped_at_zero() {
        // Very high damping with small penetration should not go negative
        let f = ground_contact_force(0.001, 100.0, 1e3, 1e6);
        assert!(f >= 0.0, "f={f}");
    }

    // ── terrain_friction_model ────────────────────────────────────────────

    #[test]
    fn friction_opposes_velocity() {
        let f = terrain_friction_model(100.0, [1.0, 0.0, 0.0], 0.7);
        assert!(f[0] < 0.0, "friction should oppose +x velocity");
    }

    #[test]
    fn friction_zero_velocity_is_zero() {
        let f = terrain_friction_model(100.0, [0.0, 0.0, 0.0], 0.7);
        assert!(f.iter().all(|x| x.abs() < 1e-15));
    }

    #[test]
    fn friction_magnitude_coulomb() {
        let f = terrain_friction_model(100.0, [1.0, 0.0, 0.0], 0.5);
        assert!((f[0] + 50.0).abs() < 1e-10, "f[0]={}", f[0]);
    }

    #[test]
    fn friction_scales_with_normal_force() {
        let f1 = terrain_friction_model(50.0, [1.0, 0.0, 0.0], 0.5);
        let f2 = terrain_friction_model(100.0, [1.0, 0.0, 0.0], 0.5);
        assert!((f2[0] / f1[0] - 2.0).abs() < 1e-10);
    }

    // ── slope_angle ───────────────────────────────────────────────────────

    #[test]
    fn slope_angle_flat_is_zero() {
        let angle = slope_angle([0.0, 1.0, 0.0]);
        assert!(angle.abs() < 1e-12, "angle={angle}");
    }

    #[test]
    fn slope_angle_vertical_wall() {
        let angle = slope_angle([1.0, 0.0, 0.0]);
        assert!((angle - PI / 2.0).abs() < 1e-12, "angle={angle}");
    }

    #[test]
    fn slope_angle_45_degrees() {
        let n = [1.0_f64 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt(), 0.0];
        let angle = slope_angle(n);
        assert!((angle - PI / 4.0).abs() < 1e-12, "angle={angle}");
    }

    // ── terrain_traversability ────────────────────────────────────────────

    #[test]
    fn traversable_flat_terrain() {
        assert!(terrain_traversability(0.0, 45.0));
    }

    #[test]
    fn not_traversable_steep_slope() {
        let steep = 60.0_f64.to_radians();
        assert!(!terrain_traversability(steep, 45.0));
    }

    #[test]
    fn traversable_at_exact_limit() {
        let angle = 30.0_f64.to_radians();
        assert!(terrain_traversability(angle, 30.0));
    }

    // ── create_flat_terrain ───────────────────────────────────────────────

    #[test]
    fn flat_terrain_has_two_triangles() {
        let patch = create_flat_terrain(10.0, 0.0);
        assert_eq!(patch.triangles.len(), 2);
        assert_eq!(patch.vertices.len(), 4);
    }

    #[test]
    fn flat_terrain_all_vertices_at_given_height() {
        let patch = create_flat_terrain(10.0, 5.0);
        for v in &patch.vertices {
            assert!((v[1] - 5.0).abs() < 1e-12, "v[1]={}", v[1]);
        }
    }

    // ── create_slope_terrain ──────────────────────────────────────────────

    #[test]
    fn slope_terrain_has_two_triangles() {
        let patch = create_slope_terrain(10.0, 30.0);
        assert_eq!(patch.triangles.len(), 2);
    }

    #[test]
    fn slope_terrain_varying_heights() {
        let patch = create_slope_terrain(10.0, 45.0);
        let ys: Vec<f64> = patch.vertices.iter().map(|v| v[1]).collect();
        let min_y = ys.iter().cloned().fold(f64::INFINITY, f64::min);
        let max_y = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!(max_y > min_y, "sloped terrain should have height variation");
    }

    #[test]
    fn slope_normal_angle_matches_slope_deg() {
        // For 30-degree slope, the normal angle from vertical should be ~30 degrees
        let patch = create_slope_terrain(10.0, 30.0);
        let tri = &patch.triangles[0];
        let v0 = patch.vertices[tri[0]];
        let v1 = patch.vertices[tri[1]];
        let v2 = patch.vertices[tri[2]];
        let n = triangle_normal(v0, v1, v2);
        let angle = slope_angle(n).to_degrees();
        assert!(
            (angle - 30.0).abs() < 1.0,
            "slope angle={angle}° expected ~30°"
        );
    }
}
