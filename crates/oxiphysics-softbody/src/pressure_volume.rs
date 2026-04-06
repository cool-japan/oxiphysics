#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Pressure-volume soft body dynamics.
//!
//! Models a closed surface mesh inflated by internal pressure. Suitable for
//! balloons, air bags, pneumatic actuators, and biological organs.
//!
//! Uses a triangle-mesh surface representation. Volume is computed using the
//! divergence theorem (signed tetrahedral decomposition).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (no nalgebra; use [f64; 3])
// ─────────────────────────────────────────────────────────────────────────────

/// Vector addition.
#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Vector subtraction.
#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar multiplication.
#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product.
#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm.
#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

/// Normalize. Returns zero vector if near-zero.
#[inline]
fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        v3_scale(a, 1.0 / n)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Mesh vertex
// ─────────────────────────────────────────────────────────────────────────────

/// A vertex on the surface mesh of a pressurized body.
#[derive(Debug, Clone)]
pub struct SurfaceVertex {
    /// World-space position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Accumulated force (N).
    pub force: [f64; 3],
    /// Inverse mass (kg⁻¹). 0.0 for static vertices.
    pub inv_mass: f64,
}

impl SurfaceVertex {
    /// Create a dynamic vertex.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        let inv_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass,
        }
    }

    /// Create a static (pinned) vertex.
    pub fn new_static(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass: 0.0,
        }
    }

    /// Returns true if this vertex is pinned.
    pub fn is_static(&self) -> bool {
        self.inv_mass == 0.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pressurized body
// ─────────────────────────────────────────────────────────────────────────────

/// A closed triangle-mesh soft body inflated by internal pressure.
///
/// The mesh must be a closed, consistently-oriented manifold for the volume
/// computation to be correct. Triangle winding should be outward-facing.
#[derive(Debug, Clone)]
pub struct PressurizedBody {
    /// Surface vertices.
    pub vertices: Vec<SurfaceVertex>,
    /// Triangle indices. Each entry is \[i0, i1, i2\] with outward normals.
    pub triangles: Vec<[usize; 3]>,
    /// Internal (gauge) pressure above ambient (Pa).
    pub pressure: f64,
    /// Rest (target) volume (m³). Used for incompressibility constraint.
    pub rest_volume: f64,
    /// Surface tensile strength (Pa). Burst criterion: max stress > this value.
    pub tensile_strength: f64,
    /// Valve flow rate coefficient (m³/(s·Pa)). Controls pressure regulation.
    pub valve_coefficient: f64,
    /// Ambient pressure (Pa). Typically 101325 Pa.
    pub ambient_pressure: f64,
}

impl PressurizedBody {
    /// Create a new `PressurizedBody`.
    ///
    /// # Arguments
    /// * `vertices` — surface vertex list
    /// * `triangles` — triangle index list (outward-facing winding)
    /// * `pressure` — initial gauge pressure (Pa)
    /// * `tensile_strength` — burst threshold (Pa)
    /// * `valve_coefficient` — flow coefficient (m³/(s·Pa))
    pub fn new(
        vertices: Vec<SurfaceVertex>,
        triangles: Vec<[usize; 3]>,
        pressure: f64,
        tensile_strength: f64,
        valve_coefficient: f64,
    ) -> Self {
        let rest_volume = compute_volume_from_verts(&vertices, &triangles);
        Self {
            vertices,
            triangles,
            pressure,
            rest_volume,
            tensile_strength,
            valve_coefficient,
            ambient_pressure: 101_325.0,
        }
    }

    /// Create a simple axis-aligned box pressurized body.
    ///
    /// Returns a box centered at the origin with given half-extents.
    pub fn unit_box(half_size: f64, pressure: f64, mass_per_vertex: f64) -> Self {
        let s = half_size;
        // 8 corners.
        let positions = [
            [-s, -s, -s],
            [s, -s, -s],
            [s, s, -s],
            [-s, s, -s],
            [-s, -s, s],
            [s, -s, s],
            [s, s, s],
            [-s, s, s],
        ];
        let vertices: Vec<SurfaceVertex> = positions
            .iter()
            .map(|&p| SurfaceVertex::new(p, mass_per_vertex))
            .collect();
        // 12 triangles (2 per face, outward normals).
        let triangles = vec![
            // -Z face
            [0, 2, 1],
            [0, 3, 2],
            // +Z face
            [4, 5, 6],
            [4, 6, 7],
            // -X face
            [0, 4, 7],
            [0, 7, 3],
            // +X face
            [1, 2, 6],
            [1, 6, 5],
            // -Y face
            [0, 1, 5],
            [0, 5, 4],
            // +Y face
            [3, 6, 2],
            [3, 7, 6],
        ];
        Self::new(vertices, triangles, pressure, 1e6, 1e-9)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Volume computation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the signed volume of a closed triangle mesh using the divergence theorem.
///
/// V = (1/6) Σ (v0 · (v1 × v2))
///
/// For a consistently-outward-wound mesh this returns the positive volume.
///
/// # Arguments
/// * `vertices` — vertex position list
/// * `triangles` — triangle index list
pub fn compute_volume_from_verts(vertices: &[SurfaceVertex], triangles: &[[usize; 3]]) -> f64 {
    let positions: Vec<[f64; 3]> = vertices.iter().map(|v| v.position).collect();
    compute_volume(&positions, triangles)
}

/// Compute signed volume from raw position arrays.
///
/// V = (1/6) Σ_{triangles} v0 · (v1 × v2)
pub fn compute_volume(positions: &[[f64; 3]], triangles: &[[usize; 3]]) -> f64 {
    let mut vol = 0.0;
    for tri in triangles {
        let v0 = positions[tri[0]];
        let v1 = positions[tri[1]];
        let v2 = positions[tri[2]];
        vol += v3_dot(v0, v3_cross(v1, v2));
    }
    vol / 6.0
}

// ─────────────────────────────────────────────────────────────────────────────
// Pressure forces
// ─────────────────────────────────────────────────────────────────────────────

/// Accumulate outward pressure forces on all surface vertices.
///
/// For each triangle, the pressure force is distributed equally to its three
/// vertices:  ΔF = p × A_triangle × n̂ / 3
///
/// Forces are added to `body.vertices[i].force`.
pub fn pressure_force(body: &mut PressurizedBody) {
    let p = body.pressure;
    if p == 0.0 {
        return;
    }
    let n_tris = body.triangles.len();
    for t in 0..n_tris {
        let tri = body.triangles[t];
        let v0 = body.vertices[tri[0]].position;
        let v1 = body.vertices[tri[1]].position;
        let v2 = body.vertices[tri[2]].position;
        let e1 = v3_sub(v1, v0);
        let e2 = v3_sub(v2, v0);
        let area_normal = v3_cross(e1, e2); // magnitude = 2 × area
        // Force = p × (area × normal) / 3 per vertex.
        let f_each = v3_scale(area_normal, p / 6.0); // p × 2A/3 / 2 = p×A/3
        for &vi in &[tri[0], tri[1], tri[2]] {
            body.vertices[vi].force = v3_add(body.vertices[vi].force, f_each);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Incompressibility constraint
// ─────────────────────────────────────────────────────────────────────────────

/// Apply incompressibility (volume-preservation) constraint.
///
/// Adjusts all dynamic vertex positions to restore the body's rest volume
/// using a global position correction along the volume gradient.
///
/// Returns the volume constraint violation C = V - V_rest before correction.
pub fn incompressibility_constraint(body: &mut PressurizedBody) -> f64 {
    let current_vol = compute_volume_from_verts(&body.vertices, &body.triangles);
    let c = current_vol - body.rest_volume;
    if c.abs() < 1e-15 {
        return c;
    }
    // Compute gradient of volume w.r.t. each vertex position.
    // dV/dx_i = (1/6) Σ (n_k) where k ranges over triangles containing vertex i.
    let nv = body.vertices.len();
    let mut grad = vec![[0.0_f64; 3]; nv];
    for tri in &body.triangles {
        let v0 = body.vertices[tri[0]].position;
        let v1 = body.vertices[tri[1]].position;
        let v2 = body.vertices[tri[2]].position;
        // dV/dv0 = (1/6)(v1 × v2), etc. (cyclic)
        let g0 = v3_scale(v3_cross(v1, v2), 1.0 / 6.0);
        let g1 = v3_scale(v3_cross(v2, v0), 1.0 / 6.0);
        let g2 = v3_scale(v3_cross(v0, v1), 1.0 / 6.0);
        grad[tri[0]] = v3_add(grad[tri[0]], g0);
        grad[tri[1]] = v3_add(grad[tri[1]], g1);
        grad[tri[2]] = v3_add(grad[tri[2]], g2);
    }
    // XPBD-style correction: Δx = -C × w_i × ∇C_i / Σ w_j |∇C_j|²
    let w_sum: f64 = (0..nv)
        .filter(|&i| !body.vertices[i].is_static())
        .map(|i| body.vertices[i].inv_mass * v3_dot(grad[i], grad[i]))
        .sum();
    if w_sum < 1e-20 {
        return c;
    }
    let lambda = -c / w_sum;
    for i in 0..nv {
        if body.vertices[i].is_static() {
            continue;
        }
        let dx = v3_scale(grad[i], lambda * body.vertices[i].inv_mass);
        body.vertices[i].position = v3_add(body.vertices[i].position, dx);
    }
    c
}

// ─────────────────────────────────────────────────────────────────────────────
// Inflation curve
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the pressure-volume relationship for a neo-Hookean membrane.
///
/// Uses the Mooney-Rivlin thin-shell approximation for a spherical balloon:
///   P(λ) = (2 t₀ / R₀) × G × (1/λ - 1/λ⁷)
///
/// where λ = R/R₀ is the stretch ratio.
///
/// # Arguments
/// * `stretch_ratio` — λ = current_radius / rest_radius
/// * `shear_modulus` — G (Pa), material shear modulus
/// * `wall_thickness` — t₀ (m), undeformed wall thickness
/// * `rest_radius` — R₀ (m), undeformed radius
///
/// Returns gauge pressure in Pa.
pub fn inflation_curve(
    stretch_ratio: f64,
    shear_modulus: f64,
    wall_thickness: f64,
    rest_radius: f64,
) -> f64 {
    if stretch_ratio <= 0.0 || rest_radius <= 0.0 {
        return 0.0;
    }
    let lambda = stretch_ratio;
    (2.0 * wall_thickness / rest_radius) * shear_modulus * (1.0 / lambda - 1.0 / lambda.powi(7))
}

// ─────────────────────────────────────────────────────────────────────────────
// Burst criterion
// ─────────────────────────────────────────────────────────────────────────────

/// Evaluate the maximum membrane stress for a pressurized body.
///
/// For a sphere: σ = p R / (2 t)
/// For a general mesh we use the per-triangle average hoop stress:
///   σ_max = max over triangles of (p × R_eff / (2 t))
///
/// # Arguments
/// * `body` — the pressurized body
/// * `wall_thickness` — membrane thickness (m)
///
/// Returns (has_burst, max_stress_pa).
pub fn burst_criterion(body: &PressurizedBody, wall_thickness: f64) -> (bool, f64) {
    // Compute approximate radius as cube root of 3V/4π for a sphere.
    let vol = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
    let r_eff = (3.0 * vol / (4.0 * std::f64::consts::PI)).cbrt();
    let sigma_max = body.pressure * r_eff / (2.0 * wall_thickness.max(1e-15));
    let has_burst = sigma_max > body.tensile_strength;
    (has_burst, sigma_max)
}

// ─────────────────────────────────────────────────────────────────────────────
// Valve flow / pressure regulation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the pressure change due to valve flow in one time step.
///
/// Models laminar (linear) flow through a valve:
///   dP/dt = -k_v × (P_internal - P_ambient)
///
/// # Arguments
/// * `body` — the pressurized body (pressure modified in place)
/// * `dt` — time step (s)
///
/// Returns the pressure delta applied.
pub fn valve_flow(body: &mut PressurizedBody, dt: f64) -> f64 {
    let delta_p = -(body.pressure - body.ambient_pressure) * body.valve_coefficient * dt;
    body.pressure += delta_p;
    delta_p
}

/// Compute the volumetric flow rate through a valve orifice (Hagen-Poiseuille approximation).
///
/// Q = k_v × ΔP  (m³/s)
///
/// # Arguments
/// * `pressure_diff` — P_internal - P_ambient (Pa)
/// * `valve_coefficient` — flow coefficient k_v (m³/(s·Pa))
pub fn volumetric_flow_rate(pressure_diff: f64, valve_coefficient: f64) -> f64 {
    valve_coefficient * pressure_diff
}

// ─────────────────────────────────────────────────────────────────────────────
// Simple integration step
// ─────────────────────────────────────────────────────────────────────────────

/// Advance the pressurized body one explicit Euler time step.
///
/// 1. Accumulate pressure forces.
/// 2. Integrate velocities and positions.
/// 3. Apply incompressibility constraint.
/// 4. Clear forces.
///
/// # Arguments
/// * `body` — the body to step
/// * `gravity` — gravitational acceleration (m/s²)
/// * `dt` — time step (s)
pub fn step_body(body: &mut PressurizedBody, gravity: [f64; 3], dt: f64) {
    pressure_force(body);
    let nv = body.vertices.len();
    for i in 0..nv {
        if body.vertices[i].is_static() {
            continue;
        }
        let m_inv = body.vertices[i].inv_mass;
        // Gravity force = m × g → impulse = g × dt (since accel = g).
        let gravity_impulse = v3_scale(gravity, dt);
        // Force impulse.
        let force_impulse = v3_scale(body.vertices[i].force, m_inv * dt);
        body.vertices[i].velocity = v3_add(
            v3_add(body.vertices[i].velocity, gravity_impulse),
            force_impulse,
        );
        body.vertices[i].position = v3_add(
            body.vertices[i].position,
            v3_scale(body.vertices[i].velocity, dt),
        );
        body.vertices[i].force = [0.0; 3];
    }
    incompressibility_constraint(body);
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-8;

    // ── Helpers ──────────────────────────────────────────────────────────────

    /// A tetrahedron whose signed volume is easy to compute by hand.
    fn tet_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        // Vertices of a regular right-corner tet with unit edges.
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        // Four outward-facing triangles.
        let triangles = vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]];
        (positions, triangles)
    }

    fn unit_box_body() -> PressurizedBody {
        PressurizedBody::unit_box(0.5, 1000.0, 0.01)
    }

    // 1. Volume of a unit tetrahedron is 1/6.
    #[test]
    fn test_tet_volume() {
        let (pos, tris) = tet_mesh();
        let vol = compute_volume(&pos, &tris).abs();
        assert!(
            (vol - 1.0 / 6.0).abs() < 1e-12,
            "Tet volume should be 1/6, got {vol}"
        );
    }

    // 2. Volume of a unit cube (side=1) is 1.0.
    #[test]
    fn test_unit_box_volume() {
        let body = PressurizedBody::unit_box(0.5, 0.0, 0.01);
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
        assert!(
            (vol - 1.0).abs() < 1e-10,
            "Unit box volume should be 1, got {vol}"
        );
    }

    // 3. Volume of box with half_size=2 is 64.
    #[test]
    fn test_box_volume_scale() {
        let body = PressurizedBody::unit_box(2.0, 0.0, 0.01);
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
        assert!(
            (vol - 64.0).abs() < 1e-8,
            "2×box volume should be 64, got {vol}"
        );
    }

    // 4. Pressure force with p=0 produces zero forces.
    #[test]
    fn test_pressure_force_zero_pressure() {
        let mut body = PressurizedBody::unit_box(0.5, 0.0, 0.01);
        pressure_force(&mut body);
        for v in &body.vertices {
            assert_eq!(v.force, [0.0; 3], "Zero pressure → zero force");
        }
    }

    // 5. Pressure force with positive pressure produces non-zero forces.
    #[test]
    fn test_pressure_force_nonzero() {
        let mut body = unit_box_body();
        pressure_force(&mut body);
        let total_force: [f64; 3] = body
            .vertices
            .iter()
            .fold([0.0; 3], |acc, v| v3_add(acc, v.force));
        // For a symmetric mesh the net force sums to zero, but individual nodes have non-zero force.
        let _net = v3_norm(total_force); // net is ~0 for symmetric body
        let any_nonzero = body.vertices.iter().any(|v| v3_norm(v.force) > EPS);
        assert!(
            any_nonzero,
            "Pressure should produce non-zero vertex forces"
        );
    }

    // 6. Net pressure force on a closed symmetric body is zero.
    #[test]
    fn test_pressure_net_force_zero_for_closed_mesh() {
        let mut body = unit_box_body();
        pressure_force(&mut body);
        let net = body
            .vertices
            .iter()
            .fold([0.0; 3], |acc, v| v3_add(acc, v.force));
        let net_mag = v3_norm(net);
        assert!(
            net_mag < 1e-8,
            "Net pressure force on closed mesh should be 0, got {net_mag}"
        );
    }

    // 7. Incompressibility constraint reduces volume error.
    #[test]
    fn test_incompressibility_reduces_error() {
        let mut body = unit_box_body();
        let rest_vol = body.rest_volume;
        // Compress the body by scaling all vertices inward.
        for v in body.vertices.iter_mut() {
            v.position = v3_scale(v.position, 0.9);
        }
        let vol_before = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
        let error_before = (vol_before - rest_vol).abs();
        incompressibility_constraint(&mut body);
        let vol_after = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
        let error_after = (vol_after - rest_vol).abs();
        assert!(
            error_after < error_before,
            "Incompressibility should reduce volume error: before={error_before}, after={error_after}"
        );
    }

    // 8. Incompressibility constraint on undeformed body returns small C.
    #[test]
    fn test_incompressibility_zero_violation() {
        let mut body = unit_box_body();
        let c = incompressibility_constraint(&mut body);
        assert!(
            c.abs() < 1e-10,
            "Undeformed body has zero volume violation, got C={c}"
        );
    }

    // 9. Inflation curve at lambda=1 is zero (no stretch → no pressure).
    #[test]
    fn test_inflation_curve_no_stretch() {
        // Actually at lambda=1: (1/1 - 1/1) = 0
        let p = inflation_curve(1.0, 1e5, 0.001, 0.05);
        assert!(
            p.abs() < EPS,
            "Inflation pressure at lambda=1 should be 0, got {p}"
        );
    }

    // 10. Inflation curve is positive for lambda > 1 (inflated balloon).
    #[test]
    fn test_inflation_curve_positive_for_inflation() {
        let p = inflation_curve(1.5, 1e5, 0.001, 0.05);
        assert!(
            p > 0.0,
            "Inflation pressure should be positive for lambda>1, got {p}"
        );
    }

    // 11. Inflation curve shows pressure peak near lambda=0.77 for neo-Hookean.
    // (known to peak around lambda ≈ 0.77 → test that p(1.2) < p(0.9) is not needed;
    //  instead test p at zero stretch ratio returns 0.)
    #[test]
    fn test_inflation_curve_zero_stretch_ratio() {
        let p = inflation_curve(0.0, 1e5, 0.001, 0.05);
        assert!(
            p.abs() < EPS,
            "Inflation at lambda=0 should return 0 (guard), got {p}"
        );
    }

    // 12. Burst criterion: no burst when pressure is low.
    #[test]
    fn test_burst_criterion_no_burst() {
        let body = PressurizedBody::unit_box(0.5, 100.0, 0.01);
        let (burst, sigma) = burst_criterion(&body, 0.01);
        assert!(
            !burst,
            "Low pressure should not burst: sigma={sigma}, strength=1e6"
        );
    }

    // 13. Burst criterion: burst when pressure exceeds tensile strength.
    #[test]
    fn test_burst_criterion_burst() {
        let mut body = unit_box_body(); // pressure=1000 Pa, strength=1e6
        // Set extreme pressure.
        body.pressure = 1e10;
        body.tensile_strength = 1.0;
        let (burst, sigma) = burst_criterion(&body, 0.001);
        assert!(
            burst,
            "Very high pressure should trigger burst: sigma={sigma}"
        );
    }

    // 14. Max stress scales with pressure.
    #[test]
    fn test_burst_criterion_stress_scales_with_pressure() {
        let mut body1 = unit_box_body();
        let mut body2 = unit_box_body();
        body1.pressure = 100.0;
        body2.pressure = 200.0;
        let (_, s1) = burst_criterion(&body1, 0.01);
        let (_, s2) = burst_criterion(&body2, 0.01);
        assert!(s2 > s1, "Higher pressure → higher stress: s1={s1}, s2={s2}");
    }

    // 15. Valve flow with zero coefficient produces no pressure change.
    #[test]
    fn test_valve_flow_zero_coefficient() {
        let mut body = unit_box_body();
        body.valve_coefficient = 0.0;
        let initial_p = body.pressure;
        valve_flow(&mut body, 0.01);
        assert_eq!(
            body.pressure, initial_p,
            "Zero valve coefficient → no pressure change"
        );
    }

    // 16. Valve flow reduces internal pressure toward ambient.
    #[test]
    fn test_valve_flow_reduces_pressure() {
        let mut body = unit_box_body(); // p=1000 > ambient=101325? No: 1000 < 101325
        // Set pressure above ambient.
        body.pressure = 200_000.0;
        body.valve_coefficient = 1e-6;
        let p_before = body.pressure;
        valve_flow(&mut body, 1.0);
        assert!(
            body.pressure < p_before,
            "Valve flow should reduce high pressure toward ambient"
        );
    }

    // 17. Volumetric flow rate is proportional to pressure difference.
    #[test]
    fn test_volumetric_flow_proportional() {
        let q1 = volumetric_flow_rate(100.0, 1e-6);
        let q2 = volumetric_flow_rate(200.0, 1e-6);
        assert!(
            (q2 - 2.0 * q1).abs() < EPS,
            "Flow should be proportional to ΔP"
        );
    }

    // 18. Volumetric flow rate with zero coefficient is zero.
    #[test]
    fn test_volumetric_flow_zero_coefficient() {
        let q = volumetric_flow_rate(1000.0, 0.0);
        assert!(q.abs() < EPS, "Zero valve coefficient → zero flow");
    }

    // 19. Step body under gravity lowers vertices (zero pressure to avoid explosion).
    #[test]
    fn test_step_body_gravity() {
        let mut body = PressurizedBody::unit_box(0.5, 0.0, 0.01); // p=0
        let initial_y: f64 =
            body.vertices.iter().map(|v| v.position[1]).sum::<f64>() / body.vertices.len() as f64;
        for _ in 0..10 {
            step_body(&mut body, [0.0, -9.81, 0.0], 1.0 / 60.0);
        }
        let final_y: f64 =
            body.vertices.iter().map(|v| v.position[1]).sum::<f64>() / body.vertices.len() as f64;
        assert!(
            final_y < initial_y - 1e-6,
            "Body should fall under gravity: initial_y={initial_y}, final_y={final_y}"
        );
    }

    // 20. SurfaceVertex static has inv_mass = 0.
    #[test]
    fn test_surface_vertex_static() {
        let v = SurfaceVertex::new_static([0.0; 3]);
        assert!(v.is_static());
        assert_eq!(v.inv_mass, 0.0);
    }

    // 21. SurfaceVertex dynamic has correct inv_mass.
    #[test]
    fn test_surface_vertex_dynamic() {
        let v = SurfaceVertex::new([0.0; 3], 0.5);
        assert!(!v.is_static());
        assert!((v.inv_mass - 2.0).abs() < EPS);
    }

    // 22. PressurizedBody::new sets rest_volume from initial geometry.
    #[test]
    fn test_pressurized_body_rest_volume() {
        let body = unit_box_body();
        assert!(
            (body.rest_volume.abs() - 1.0).abs() < 1e-8,
            "rest_volume should be 1.0 for unit box"
        );
    }

    // 23. compute_volume returns positive value for outward-wound mesh.
    #[test]
    fn test_compute_volume_positive() {
        let (pos, tris) = tet_mesh();
        let vol = compute_volume(&pos, &tris);
        assert!(
            vol > 0.0,
            "Outward-wound mesh should give positive volume, got {vol}"
        );
    }

    // 24. Incompressibility applied multiple times converges to near-zero C.
    #[test]
    fn test_incompressibility_idempotent() {
        let mut body = unit_box_body();
        for v in body.vertices.iter_mut() {
            v.position = v3_scale(v.position, 0.85);
        }
        // Run several iterations to converge (XPBD constraint may need multiple passes).
        for _ in 0..20 {
            incompressibility_constraint(&mut body);
        }
        let c_final = incompressibility_constraint(&mut body);
        assert!(
            c_final.abs() < 1e-4,
            "After many iterations, C should be tiny, got {c_final}"
        );
    }

    // 25. Inflation curve: larger shear modulus → higher pressure.
    #[test]
    fn test_inflation_curve_shear_modulus_effect() {
        let p1 = inflation_curve(1.5, 1e4, 0.001, 0.05);
        let p2 = inflation_curve(1.5, 2e4, 0.001, 0.05);
        assert!(p2 > p1, "Higher shear modulus → higher inflation pressure");
    }

    // 26. Valve flow with negative pressure (below ambient) increases pressure.
    #[test]
    fn test_valve_flow_increases_below_ambient() {
        let mut body = unit_box_body();
        body.pressure = 50_000.0; // below ambient 101325
        body.valve_coefficient = 1e-6;
        let p_before = body.pressure;
        valve_flow(&mut body, 1.0);
        assert!(
            body.pressure > p_before,
            "Pressure below ambient should increase toward ambient"
        );
    }

    // 27. Burst criterion: stress increases with reduced wall thickness.
    #[test]
    fn test_burst_thin_wall_higher_stress() {
        let body = unit_box_body();
        let (_, s_thick) = burst_criterion(&body, 0.01);
        let (_, s_thin) = burst_criterion(&body, 0.001);
        assert!(
            s_thin > s_thick,
            "Thinner wall → higher stress: s_thick={s_thick}, s_thin={s_thin}"
        );
    }

    // 28. pressure_force: forces point outward from center.
    #[test]
    fn test_pressure_force_outward() {
        let mut body = unit_box_body();
        pressure_force(&mut body);
        // For each vertex of the unit box, the force should point away from origin
        // (same sign as position component for at least one axis).
        for v in &body.vertices {
            let outward = v3_dot(v.force, v.position);
            // At least loosely outward.
            assert!(
                outward >= -1e-8,
                "Pressure force should be outward, got dot={outward} for pos={:?}",
                v.position
            );
        }
    }

    // 29. Volume doubles when all positions scaled by cbrt(2).
    #[test]
    fn test_volume_scales_with_positions() {
        let body = unit_box_body();
        let scale = 2.0_f64.cbrt();
        let scaled_verts: Vec<SurfaceVertex> = body
            .vertices
            .iter()
            .map(|v| {
                let mut nv = v.clone();
                nv.position = v3_scale(v.position, scale);
                nv
            })
            .collect();
        let vol = compute_volume_from_verts(&scaled_verts, &body.triangles).abs();
        assert!(
            (vol - 2.0).abs() < 1e-8,
            "Scaled volume should be 2.0, got {vol}"
        );
    }

    // 30. Step body preserves approximate volume (incompressibility), zero pressure.
    #[test]
    fn test_step_body_volume_preservation() {
        let mut body = PressurizedBody::unit_box(0.5, 0.0, 0.01); // zero pressure
        let rest_vol = body.rest_volume;
        for _ in 0..5 {
            step_body(&mut body, [0.0, 0.0, 0.0], 1.0 / 60.0);
        }
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
        let rel_err = (vol - rest_vol).abs() / rest_vol;
        assert!(
            rel_err < 0.05,
            "Volume should stay close to rest after steps: rel_err={rel_err}"
        );
    }
}
