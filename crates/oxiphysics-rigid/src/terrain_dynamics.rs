// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body dynamics on terrain: height maps, ground contact, Coulomb
//! friction on uneven surfaces, slope resistance, terrain compliance, and
//! vehicle tip-over stability.

// ---------------------------------------------------------------------------
// Helper math (no nalgebra)
// ---------------------------------------------------------------------------

fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

fn vec3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = vec3_norm(a);
    if n < 1e-12 {
        [0.0, 1.0, 0.0]
    } else {
        vec3_scale(a, 1.0 / n)
    }
}

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// TerrainPatch
// ---------------------------------------------------------------------------

/// A rectangular terrain patch described by a height map grid.
///
/// The patch is axis-aligned in the XZ plane.  `width` and `depth` define
/// the size in the X and Z directions respectively.  The height map stores
/// one height value per grid cell `(row, col)` in row-major order.
#[derive(Debug, Clone)]
pub struct TerrainPatch {
    /// World-space width of the patch along X (m).
    pub width: f64,
    /// World-space depth of the patch along Z (m).
    pub depth: f64,
    /// Number of sample columns along X.
    pub cols: usize,
    /// Number of sample rows along Z.
    pub rows: usize,
    /// Height samples in row-major order: `heights[row * cols + col]`.
    pub heights: Vec<f64>,
    /// Coulomb friction coefficient (dimensionless).
    pub friction: f64,
    /// Terrain compliance (m/N) — how much the ground deforms per unit force.
    pub compliance: f64,
}

impl TerrainPatch {
    /// Create a flat terrain patch with uniform height `h`.
    pub fn flat(width: f64, depth: f64, cols: usize, rows: usize, h: f64, friction: f64) -> Self {
        Self {
            width,
            depth,
            cols,
            rows,
            heights: vec![h; cols * rows],
            friction,
            compliance: 0.0,
        }
    }

    /// Sample the height at world coordinates `(wx, wz)` using bilinear
    /// interpolation.  Clamps to the patch boundary.
    pub fn sample_height(&self, wx: f64, wz: f64) -> f64 {
        if self.cols < 2 || self.rows < 2 {
            return self.heights[0];
        }
        let u = (wx / self.width).clamp(0.0, 1.0) * (self.cols as f64 - 1.0);
        let v = (wz / self.depth).clamp(0.0, 1.0) * (self.rows as f64 - 1.0);
        let c0 = u.floor() as usize;
        let r0 = v.floor() as usize;
        let c1 = (c0 + 1).min(self.cols - 1);
        let r1 = (r0 + 1).min(self.rows - 1);
        let fc = u - c0 as f64;
        let fr = v - r0 as f64;
        let h00 = self.heights[r0 * self.cols + c0];
        let h10 = self.heights[r0 * self.cols + c1];
        let h01 = self.heights[r1 * self.cols + c0];
        let h11 = self.heights[r1 * self.cols + c1];
        h00 * (1.0 - fc) * (1.0 - fr)
            + h10 * fc * (1.0 - fr)
            + h01 * (1.0 - fc) * fr
            + h11 * fc * fr
    }
}

// ---------------------------------------------------------------------------
// terrain_normal
// ---------------------------------------------------------------------------

/// Compute the surface normal of a terrain patch at world coordinates
/// `(wx, wz)` using finite-difference gradients of the height map.
///
/// The normal is guaranteed to point upward (positive Y component).
///
/// # Arguments
/// * `patch` - The terrain patch.
/// * `wx` - World X coordinate.
/// * `wz` - World Z coordinate.
///
/// Returns a unit normal vector `[nx, ny, nz]`.
pub fn terrain_normal(patch: &TerrainPatch, wx: f64, wz: f64) -> [f64; 3] {
    let eps = patch.width / (patch.cols as f64).max(2.0);
    let h_px = patch.sample_height(wx + eps, wz);
    let h_mx = patch.sample_height(wx - eps, wz);
    let h_pz = patch.sample_height(wx, wz + eps);
    let h_mz = patch.sample_height(wx, wz - eps);
    // Tangent vectors along X and Z
    let tx: [f64; 3] = [2.0 * eps, h_px - h_mx, 0.0];
    let tz: [f64; 3] = [0.0, h_pz - h_mz, 2.0 * eps];
    let n = vec3_cross(tx, tz);
    // Ensure upward-facing
    let n = if n[1] < 0.0 { vec3_scale(n, -1.0) } else { n };
    vec3_normalise(n)
}

// ---------------------------------------------------------------------------
// ground_contact_point
// ---------------------------------------------------------------------------

/// Find the lowest world-space point of a body that may be touching the
/// terrain.
///
/// The body is modelled as an axis-aligned box centred at `centre` with
/// half-extents `half_extents`.  The function returns the vertex of the box
/// with the smallest world Y coordinate (i.e., the lowest point), and the
/// terrain height at its XZ position.
///
/// Returns `(contact_point, terrain_height)`.
pub fn ground_contact_point(
    centre: [f64; 3],
    half_extents: [f64; 3],
    patch: &TerrainPatch,
) -> ([f64; 3], f64) {
    let [cx, cy, cz] = centre;
    let [hx, _hy, hz] = half_extents;
    // Lowest corners in XZ
    let candidates = [
        [cx - hx, cz - hz],
        [cx + hx, cz - hz],
        [cx - hx, cz + hz],
        [cx + hx, cz + hz],
    ];
    let mut best_y = f64::INFINITY;
    let mut best_pt = [cx, cy, cz];
    let mut best_th = 0.0;
    for [wx, wz] in candidates {
        let th = patch.sample_height(wx, wz);
        let pt_y = cy - half_extents[1];
        if pt_y < best_y {
            best_y = pt_y;
            best_pt = [wx, pt_y, wz];
            best_th = th;
        }
    }
    (best_pt, best_th)
}

// ---------------------------------------------------------------------------
// terrain_friction_force
// ---------------------------------------------------------------------------

/// Compute the Coulomb friction force on a body sliding over uneven terrain.
///
/// # Arguments
/// * `normal` - Terrain surface normal `[nx, ny, nz]` (unit vector).
/// * `velocity` - Body velocity `[vx, vy, vz]` (m/s).
/// * `normal_force` - Magnitude of the normal contact force (N).
/// * `mu` - Coulomb friction coefficient.
///
/// Returns the friction force vector `[fx, fy, fz]` opposing the tangential
/// velocity.
pub fn terrain_friction_force(
    normal: [f64; 3],
    velocity: [f64; 3],
    normal_force: f64,
    mu: f64,
) -> [f64; 3] {
    // Tangential component of velocity
    let v_dot_n = vec3_dot(velocity, normal);
    let v_normal = vec3_scale(normal, v_dot_n);
    let v_tan = vec3_sub(velocity, v_normal);
    let v_tan_mag = vec3_norm(v_tan);
    if v_tan_mag < 1e-12 {
        return [0.0; 3];
    }
    let friction_mag = mu * normal_force.max(0.0);
    let direction = vec3_scale(v_tan, -1.0 / v_tan_mag);
    vec3_scale(direction, friction_mag)
}

// ---------------------------------------------------------------------------
// slope_resistance
// ---------------------------------------------------------------------------

/// Compute the gravitational resistance force opposing motion up a slope.
///
/// Projects the gravity vector onto the terrain surface, then takes its
/// component along `motion_dir`.  A positive result means the slope
/// opposes the motion (body going uphill); negative means gravity aids it
/// (downhill).
///
/// # Arguments
/// * `mass` - Body mass (kg).
/// * `normal` - Terrain normal `[nx, ny, nz]` (unit vector).
/// * `gravity` - Gravitational acceleration magnitude (m/s²), positive.
/// * `motion_dir` - Unit vector of the desired motion direction `[dx, dy, dz]`.
///
/// Returns the resistance force scalar (N).
pub fn slope_resistance(mass: f64, normal: [f64; 3], gravity: f64, motion_dir: [f64; 3]) -> f64 {
    // Gravity vector (world-space, pointing downward)
    let g_vec = [0.0_f64, -mass * gravity, 0.0];
    // Remove the normal component → tangential gravity on the surface
    let g_dot_n = vec3_dot(g_vec, normal);
    let g_normal_part = vec3_scale(normal, g_dot_n);
    let g_tangential = vec3_sub(g_vec, g_normal_part);
    // Component along motion_dir: positive if gravity pushes in motion direction (downhill)
    let g_along_motion = vec3_dot(g_tangential, motion_dir);
    // Resistance is positive when going uphill (gravity opposes motion),
    // negative when going downhill (gravity aids motion).
    -g_along_motion
}

// ---------------------------------------------------------------------------
// terrain_compliance
// ---------------------------------------------------------------------------

/// Compute the penetration depth and restoring force for soft ground contact.
///
/// # Arguments
/// * `body_bottom_y` - World Y of the body's lowest point (m).
/// * `terrain_height` - Terrain height at the contact XZ position (m).
/// * `compliance` - Terrain compliance in m/N.
/// * `contact_area` - Effective contact area (m²) for pressure calculation.
///
/// Returns `(penetration_depth, restoring_force_n)`.
pub fn terrain_compliance(
    body_bottom_y: f64,
    terrain_height: f64,
    compliance: f64,
    contact_area: f64,
) -> (f64, f64) {
    let penetration = (terrain_height - body_bottom_y).max(0.0);
    if penetration < 1e-12 || compliance < 1e-15 {
        return (penetration, penetration / compliance.max(1e-15));
    }
    // Pressure = penetration / (compliance × contact_area)
    let stiffness = contact_area / compliance;
    let force = stiffness * penetration;
    (penetration, force)
}

// ---------------------------------------------------------------------------
// vehicle_stability
// ---------------------------------------------------------------------------

/// Compute whether a vehicle is at risk of tipping over on sloped terrain.
///
/// Uses the Static Stability Factor (SSF): the ratio of half the track width
/// to the height of the centre of mass above the roll axis.  If the lateral
/// acceleration exceeds the tip-over threshold, the vehicle is unstable.
///
/// # Arguments
/// * `track_width` - Distance between left and right wheel contact points (m).
/// * `cg_height` - Height of the centre of gravity above the ground (m).
/// * `lateral_acceleration` - Lateral acceleration in m/s².
/// * `gravity` - Gravitational acceleration (m/s²).
///
/// Returns `(ssf, is_stable)` where `ssf` is the Static Stability Factor and
/// `is_stable` is `true` when the vehicle is not tipping.
pub fn vehicle_stability(
    track_width: f64,
    cg_height: f64,
    lateral_acceleration: f64,
    gravity: f64,
) -> (f64, bool) {
    let ssf = track_width / (2.0 * cg_height.max(1e-6));
    let tip_acceleration = ssf * gravity;
    let is_stable = lateral_acceleration.abs() < tip_acceleration;
    (ssf, is_stable)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── TerrainPatch helpers ───────────────────────────────────────────────

    fn flat_patch(h: f64) -> TerrainPatch {
        TerrainPatch::flat(10.0, 10.0, 5, 5, h, 0.8)
    }

    // ── sample_height ──────────────────────────────────────────────────────

    #[test]
    fn flat_terrain_uniform_height() {
        let p = flat_patch(3.0);
        assert!((p.sample_height(5.0, 5.0) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn sample_height_clamps_to_boundary() {
        let p = flat_patch(1.0);
        let h = p.sample_height(-999.0, -999.0);
        assert!((h - 1.0).abs() < 1e-10);
    }

    #[test]
    fn sample_height_corner_exact() {
        let mut p = flat_patch(0.0);
        // Set top-left corner to 5.0
        p.heights[0] = 5.0;
        let h = p.sample_height(0.0, 0.0);
        assert!((h - 5.0).abs() < 1e-6);
    }

    #[test]
    fn sample_height_bilinear_midpoint() {
        let mut p = TerrainPatch::flat(2.0, 2.0, 2, 2, 0.0, 0.5);
        p.heights[0] = 0.0;
        p.heights[1] = 4.0;
        p.heights[2] = 0.0;
        p.heights[3] = 4.0;
        // Midpoint should be 2.0
        let h = p.sample_height(1.0, 1.0);
        assert!((h - 2.0).abs() < 1e-6, "h={h}");
    }

    // ── terrain_normal ─────────────────────────────────────────────────────

    #[test]
    fn flat_terrain_normal_points_up() {
        let p = flat_patch(0.0);
        let n = terrain_normal(&p, 5.0, 5.0);
        assert!(n[1] > 0.9, "ny={}", n[1]);
    }

    #[test]
    fn terrain_normal_is_unit_length() {
        let p = flat_patch(0.0);
        let n = terrain_normal(&p, 5.0, 5.0);
        let len = vec3_norm(n);
        assert!((len - 1.0).abs() < 1e-6, "len={len}");
    }

    #[test]
    fn sloped_terrain_normal_tilts() {
        let mut p = TerrainPatch::flat(10.0, 10.0, 3, 3, 0.0, 0.5);
        // Create a slope rising in X: heights increase with column index
        for row in 0..3 {
            for col in 0..3 {
                p.heights[row * 3 + col] = col as f64 * 2.0;
            }
        }
        let n = terrain_normal(&p, 5.0, 5.0);
        // Normal should tilt toward negative X
        assert!(
            n[0] < 0.0,
            "nx should be negative for upslope in +X, got {}",
            n[0]
        );
        assert!(n[1] > 0.0, "ny should be positive");
    }

    #[test]
    fn terrain_normal_boundary_does_not_panic() {
        let p = flat_patch(1.0);
        let _n = terrain_normal(&p, 0.0, 0.0);
        let _n2 = terrain_normal(&p, 10.0, 10.0);
    }

    // ── ground_contact_point ───────────────────────────────────────────────

    #[test]
    fn contact_point_above_flat_terrain() {
        let p = flat_patch(0.0);
        let (pt, th) = ground_contact_point([5.0, 1.0, 5.0], [0.5, 0.5, 0.5], &p);
        assert!((th - 0.0).abs() < 1e-10);
        assert!(pt[1] <= 1.0);
    }

    #[test]
    fn contact_point_y_is_body_bottom() {
        let p = flat_patch(0.0);
        let cy = 2.0;
        let hy = 0.5;
        let (pt, _th) = ground_contact_point([5.0, cy, 5.0], [1.0, hy, 1.0], &p);
        assert!((pt[1] - (cy - hy)).abs() < 1e-10, "pt.y={}", pt[1]);
    }

    #[test]
    fn contact_point_no_panic_at_edge() {
        let p = flat_patch(0.0);
        let _ = ground_contact_point([0.0, 1.0, 0.0], [0.5, 0.5, 0.5], &p);
    }

    // ── terrain_friction_force ─────────────────────────────────────────────

    #[test]
    fn friction_zero_velocity_zero_force() {
        let n = [0.0, 1.0, 0.0];
        let f = terrain_friction_force(n, [0.0, 0.0, 0.0], 100.0, 0.8);
        assert!(vec3_norm(f) < 1e-10);
    }

    #[test]
    fn friction_flat_surface_opposes_motion() {
        let n = [0.0, 1.0, 0.0];
        let v = [1.0, 0.0, 0.0];
        let f = terrain_friction_force(n, v, 100.0, 0.8);
        assert!(f[0] < 0.0, "friction should oppose +x motion");
        assert!((f[0] + 80.0).abs() < 1e-8, "f[0]={}", f[0]);
    }

    #[test]
    fn friction_magnitude_proportional_to_normal_force() {
        let n = [0.0, 1.0, 0.0];
        let v = [1.0, 0.0, 0.0];
        let f1 = terrain_friction_force(n, v, 100.0, 0.5);
        let f2 = terrain_friction_force(n, v, 200.0, 0.5);
        assert!((vec3_norm(f2) - 2.0 * vec3_norm(f1)).abs() < 1e-8);
    }

    #[test]
    fn friction_normal_velocity_no_tangential_friction() {
        let n = [0.0, 1.0, 0.0];
        // Velocity purely in normal direction
        let v = [0.0, 5.0, 0.0];
        let f = terrain_friction_force(n, v, 100.0, 0.8);
        assert!(vec3_norm(f) < 1e-10);
    }

    #[test]
    fn friction_zero_normal_force_zero_friction() {
        let n = [0.0, 1.0, 0.0];
        let f = terrain_friction_force(n, [3.0, 0.0, 0.0], 0.0, 0.8);
        assert!(vec3_norm(f) < 1e-10);
    }

    // ── slope_resistance ───────────────────────────────────────────────────

    #[test]
    fn slope_resistance_flat_terrain_zero() {
        let n = [0.0, 1.0, 0.0]; // flat
        let motion = [1.0, 0.0, 0.0];
        let r = slope_resistance(1000.0, n, 9.81, motion);
        assert!(r.abs() < 1e-6, "r={r}");
    }

    #[test]
    fn slope_resistance_downhill_negative() {
        // Slope: surface tilted so +X is downhill.
        // For +X downhill the normal tilts toward +X: n = [sin θ, cos θ, 0].
        // Downhill tangent on the slope: t = [cos θ, -sin θ, 0].
        // motion_dir along downhill tangent → gravity aids → resistance < 0.
        let angle: f64 = 15.0_f64.to_radians();
        let n = [angle.sin(), angle.cos(), 0.0];
        let motion = [angle.cos(), -angle.sin(), 0.0]; // downhill tangent
        let r = slope_resistance(100.0, n, 9.81, motion);
        assert!(r < 0.0, "r={r}");
    }

    #[test]
    fn slope_resistance_uphill_positive() {
        // Same slope, motion in uphill direction (negative of downhill tangent).
        let angle: f64 = 15.0_f64.to_radians();
        let n = [angle.sin(), angle.cos(), 0.0];
        let motion = [-angle.cos(), angle.sin(), 0.0]; // uphill tangent
        let r = slope_resistance(100.0, n, 9.81, motion);
        assert!(r > 0.0, "r={r}");
    }

    #[test]
    fn slope_resistance_scales_with_mass() {
        // On a flat surface (normal = up), both should be ~0 regardless of mass.
        let n = [0.0, 1.0, 0.0];
        let motion = [1.0, 0.0, 0.0];
        let r1 = slope_resistance(100.0, n, 9.81, motion);
        let r2 = slope_resistance(200.0, n, 9.81, motion);
        assert!((r1 - r2).abs() < 1e-6, "flat terrain: both should be ~0");
    }

    // ── terrain_compliance ─────────────────────────────────────────────────

    #[test]
    fn compliance_no_penetration_zero_force() {
        let (pen, force) = terrain_compliance(1.0, 0.5, 1e-3, 1.0);
        assert!(pen < 1e-10, "pen={pen}");
        assert!(force < 1e-6, "force={force}");
    }

    #[test]
    fn compliance_penetration_increases_force() {
        let (_, f1) = terrain_compliance(0.0, 0.01, 1e-5, 1.0);
        let (_, f2) = terrain_compliance(0.0, 0.02, 1e-5, 1.0);
        assert!(f2 > f1, "deeper penetration should give more force");
    }

    #[test]
    fn compliance_rigid_ground_large_force() {
        // Very low compliance → large stiffness → large force
        let (pen, force) = terrain_compliance(0.0, 0.01, 1e-10, 1.0);
        assert!(pen > 0.0);
        assert!(force > 1.0, "force={force}");
    }

    #[test]
    fn compliance_area_scales_force() {
        let (_, f1) = terrain_compliance(0.0, 0.01, 1e-3, 1.0);
        let (_, f2) = terrain_compliance(0.0, 0.01, 1e-3, 2.0);
        assert!((f2 - 2.0 * f1).abs() < 1e-6, "f1={f1} f2={f2}");
    }

    // ── vehicle_stability ──────────────────────────────────────────────────

    #[test]
    fn stability_wide_low_cg_stable() {
        let (ssf, stable) = vehicle_stability(2.0, 0.5, 2.0, 9.81);
        assert!(ssf > 1.0);
        assert!(stable);
    }

    #[test]
    fn stability_narrow_high_cg_unstable() {
        let (_, stable) = vehicle_stability(1.0, 2.0, 5.0, 9.81);
        assert!(
            !stable,
            "narrow, high CG at high lateral accel should be unstable"
        );
    }

    #[test]
    fn stability_ssf_formula() {
        let track = 1.6;
        let cg = 0.8;
        let (ssf, _) = vehicle_stability(track, cg, 0.0, 9.81);
        let expected = track / (2.0 * cg);
        assert!((ssf - expected).abs() < 1e-10);
    }

    #[test]
    fn stability_zero_lateral_always_stable() {
        let (_, stable) = vehicle_stability(1.2, 1.5, 0.0, 9.81);
        assert!(stable);
    }

    #[test]
    fn stability_returns_ssf_positive() {
        let (ssf, _) = vehicle_stability(2.0, 1.0, 1.0, 9.81);
        assert!(ssf > 0.0);
    }

    #[test]
    fn flat_patch_construction() {
        let p = TerrainPatch::flat(10.0, 20.0, 4, 6, 2.5, 0.7);
        assert_eq!(p.cols, 4);
        assert_eq!(p.rows, 6);
        assert!((p.friction - 0.7).abs() < 1e-10);
        assert_eq!(p.heights.len(), 24);
    }

    #[test]
    fn terrain_compliance_penetration_depth_returned() {
        let (pen, _) = terrain_compliance(0.0, 0.05, 1e-3, 1.0);
        assert!((pen - 0.05).abs() < 1e-10);
    }

    #[test]
    fn friction_tilt_normal_gives_smaller_tangential() {
        // Normal tilted 45°
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let n = [s, s, 0.0];
        let v = [1.0, 0.0, 0.0]; // horizontal motion
        let f_tilt = terrain_friction_force(n, v, 100.0, 1.0);
        // Compare with flat normal
        let nf = [0.0, 1.0, 0.0];
        let f_flat = terrain_friction_force(nf, v, 100.0, 1.0);
        // Tilted surface has smaller tangential component for same velocity
        assert!(vec3_norm(f_tilt) < vec3_norm(f_flat) + 1e-6);
    }
}
