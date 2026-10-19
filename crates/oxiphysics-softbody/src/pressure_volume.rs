// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Pressure-volume soft body dynamics.
//!
//! Models a closed surface mesh inflated by internal pressure. Suitable for
//! balloons, air bags, pneumatic actuators, and biological organs.
//!
//! Uses a triangle-mesh surface representation. Volume is computed using the
//! divergence theorem (signed tetrahedral decomposition).
//!
//! # Dynamics
//!
//! [`step_body`] integrates the surface vertices with semi-implicit Euler under
//! gravity, the membrane edge springs ([`membrane_force`]) and, for an
//! inflatable body, the pressure load ([`pressure_force`]). Two regimes exist,
//! selected by [`PressurizedBody::enforce_incompressible`]:
//!
//! * **Incompressible** (default): the enclosed volume is held at
//!   [`PressurizedBody::rest_volume`] by iterated position projection followed
//!   by the PBD velocity update `v = (x − x_start)/Δt`. On a closed mesh the
//!   nodal pressure force is exactly `p ∂V/∂xᵢ`, i.e. it is parallel to the
//!   volume-constraint gradient; the constraint reaction therefore absorbs it
//!   completely and it is not integrated as a load.
//! * **Inflatable**: the volume is free; the gauge pressure inflates the
//!   membrane until the edge-spring tension balances it
//!   (`p ∂V/∂xᵢ + F_spring,i = 0`).

use std::collections::BTreeSet;

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
#[cfg(test)]
#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
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

/// An elastic edge of the membrane: a linear spring between two surface
/// vertices with a fixed rest length.
#[derive(Debug, Clone, PartialEq)]
pub struct MembraneEdge {
    /// Indices of the two end vertices (`vertices[0] < vertices[1]`).
    pub vertices: [usize; 2],
    /// Rest length (m).
    pub rest_length: f64,
}

/// Collect the unique undirected edges of `triangles` with their current
/// lengths as rest lengths. Edges referencing missing vertices are skipped.
fn build_membrane_edges(vertices: &[SurfaceVertex], triangles: &[[usize; 3]]) -> Vec<MembraneEdge> {
    let mut unique: BTreeSet<(usize, usize)> = BTreeSet::new();
    for tri in triangles {
        for (a, b) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            unique.insert((a.min(b), a.max(b)));
        }
    }
    unique
        .into_iter()
        .filter_map(|(a, b)| {
            let pa = vertices.get(a)?.position;
            let pb = vertices.get(b)?.position;
            let d = v3_sub(pb, pa);
            Some(MembraneEdge {
                vertices: [a, b],
                rest_length: v3_dot(d, d).sqrt(),
            })
        })
        .collect()
}

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
    /// Membrane edge-spring stiffness k (N/m). `0.0` (the default) disables
    /// membrane elasticity.
    pub membrane_stiffness: f64,
    /// Membrane edge springs; built by [`PressurizedBody::new`] from the
    /// triangle list with the construction-time edge lengths as rest lengths.
    /// Call [`PressurizedBody::rebuild_membrane_edges`] after editing
    /// `vertices` or `triangles` to re-capture the rest lengths.
    pub membrane_edges: Vec<MembraneEdge>,
    /// If `true` (the default) [`step_body`] holds the enclosed volume at
    /// `rest_volume` (incompressible fill); if `false` the volume is free and
    /// `pressure` inflates the membrane against `membrane_stiffness`.
    pub enforce_incompressible: bool,
    /// Viscous velocity damping rate (1/s): each step multiplies the vertex
    /// velocities by `1 / (1 + damping · dt)`. `0.0` (the default) is
    /// undamped; an inflating membrane needs some damping to settle.
    pub damping: f64,
}

impl PressurizedBody {
    /// Create a new `PressurizedBody`.
    ///
    /// The body starts incompressible (`enforce_incompressible = true`) with
    /// no membrane elasticity and no damping; edge rest lengths are captured
    /// from the initial geometry.
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
        let membrane_edges = build_membrane_edges(&vertices, &triangles);
        Self {
            vertices,
            triangles,
            pressure,
            rest_volume,
            tensile_strength,
            valve_coefficient,
            ambient_pressure: 101_325.0,
            membrane_stiffness: 0.0,
            membrane_edges,
            enforce_incompressible: true,
            damping: 0.0,
        }
    }

    /// Re-capture the membrane edge table and rest lengths from the current
    /// vertex positions and triangle list.
    pub fn rebuild_membrane_edges(&mut self) {
        self.membrane_edges = build_membrane_edges(&self.vertices, &self.triangles);
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

/// Accumulate the membrane edge-spring forces on all surface vertices.
///
/// Each [`MembraneEdge`] is a linear spring `F = k (ℓ − ℓ₀) d̂` acting along
/// the edge (`k` = [`PressurizedBody::membrane_stiffness`]); the two end
/// vertices receive equal and opposite forces, so the net force is zero.
/// Zero-length edges (undefined direction) and edges referencing missing
/// vertices are skipped. Does nothing when the stiffness is not positive.
///
/// Forces are added to `body.vertices[i].force`.
pub fn membrane_force(body: &mut PressurizedBody) {
    let k = body.membrane_stiffness;
    if k.is_nan() || k <= 0.0 {
        return;
    }
    for edge in &body.membrane_edges {
        let [a, b] = edge.vertices;
        let (Some(va), Some(vb)) = (body.vertices.get(a), body.vertices.get(b)) else {
            continue;
        };
        let d = v3_sub(vb.position, va.position);
        let len = v3_dot(d, d).sqrt();
        if len <= f64::MIN_POSITIVE {
            continue;
        }
        // Tension (ℓ > ℓ₀) pulls the end vertices together.
        let f = v3_scale(d, k * (len - edge.rest_length) / len);
        body.vertices[a].force = v3_add(body.vertices[a].force, f);
        body.vertices[b].force = v3_sub(body.vertices[b].force, f);
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
    for (i, &gi) in grad.iter().enumerate().take(nv) {
        if body.vertices[i].is_static() {
            continue;
        }
        let dx = v3_scale(gi, lambda * body.vertices[i].inv_mass);
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
/// `pressure` is the **gauge** pressure above ambient, so the valve relaxes it
/// toward zero. The volumetric flow through the valve is linear in the gauge
/// pressure, `Q = k_v · p` (m³/s, positive = gas leaving), and the enclosed
/// gas is isothermal, `P_abs · V = const`, hence
///
/// ```text
/// dp/dt = −(P_abs / V) · Q = −(P_abs / V) · k_v · p,    P_abs = p + P_ambient
/// ```
///
/// `V` is the current enclosed volume (rest volume as a fallback for a
/// degenerate mesh). The step is explicit Euler, clamped so the gauge pressure
/// never crosses zero: a body above ambient relaxes down to — never below —
/// ambient, a body below ambient (`p < 0`, gas flowing in) relaxes up to it,
/// and `p = 0` stays at zero. A non-positive/NaN `valve_coefficient`, a
/// non-positive/NaN `dt`, a non-finite pressure or a body without volume leave
/// the pressure unchanged.
///
/// # Arguments
/// * `body` — the pressurized body (pressure modified in place)
/// * `dt` — time step (s)
///
/// Returns the pressure delta applied.
pub fn valve_flow(body: &mut PressurizedBody, dt: f64) -> f64 {
    let p = body.pressure;
    let k_v = body.valve_coefficient;
    if !p.is_finite() || p == 0.0 || k_v.is_nan() || k_v <= 0.0 || dt.is_nan() || dt <= 0.0 {
        return 0.0;
    }
    let mut volume = compute_volume_from_verts(&body.vertices, &body.triangles).abs();
    if volume.is_nan() || volume <= f64::MIN_POSITIVE {
        volume = body.rest_volume.abs();
    }
    if !volume.is_finite() || volume <= f64::MIN_POSITIVE {
        return 0.0;
    }
    let p_abs = (p + body.ambient_pressure).max(0.0);
    let flow = volumetric_flow_rate(p, k_v);
    let raw_delta = -(p_abs / volume) * flow * dt;
    // Relax toward zero gauge pressure without overshooting it.
    let delta_p = if p > 0.0 {
        raw_delta.max(-p)
    } else {
        raw_delta.min(-p)
    };
    body.pressure = p + delta_p;
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

/// Maximum number of volume-projection passes per [`step_body`] call.
const MAX_VOLUME_PROJECTIONS: usize = 32;

/// Relative volume tolerance `|V − V_rest| / |V_rest|` at which the projection
/// loop of [`step_body`] stops.
const VOLUME_PROJECTION_TOLERANCE: f64 = 1.0e-12;

/// Advance the pressurized body one semi-implicit Euler time step.
///
/// 1. Accumulate forces on top of whatever is already in `vertex.force`:
///    the membrane edge springs ([`membrane_force`]) and, **unless**
///    `enforce_incompressible` is set, the pressure load
///    ([`pressure_force`]).
/// 2. Integrate velocities (gravity + forces, then viscous damping) and
///    positions of the dynamic vertices.
/// 3. If `enforce_incompressible`: iterate [`incompressibility_constraint`]
///    until the volume matches `rest_volume` (relative tolerance 1e-12, at
///    most 32 passes), then apply the PBD velocity update
///    `v = (x_final − x_start)/Δt` to the dynamic vertices, so the velocity
///    field carries no volume-changing component.
/// 4. Clear the forces of **all** vertices, static ones included.
///
/// In the incompressible regime the pressure load is not integrated: on a
/// closed mesh `pressure_force` equals `p ∂V/∂xᵢ`, which is parallel to the
/// constraint gradient and would be removed again by the projection (while a
/// large `p` would first throw the vertices far outside the linearisation
/// range of the constraint). In the inflatable regime the pressure inflates
/// the membrane against its edge springs.
///
/// A non-positive or NaN `dt` leaves the body untouched.
///
/// # Arguments
/// * `body` — the body to step
/// * `gravity` — gravitational acceleration (m/s²)
/// * `dt` — time step (s)
pub fn step_body(body: &mut PressurizedBody, gravity: [f64; 3], dt: f64) {
    if dt.is_nan() || dt <= 0.0 {
        return;
    }
    if !body.enforce_incompressible {
        pressure_force(body);
    }
    membrane_force(body);

    let enforce = body.enforce_incompressible;
    let x_start: Vec<[f64; 3]> = if enforce {
        body.vertices.iter().map(|v| v.position).collect()
    } else {
        Vec::new()
    };
    let damping_factor = if body.damping > 0.0 {
        1.0 / (1.0 + body.damping * dt)
    } else {
        1.0
    };

    let gravity_impulse = v3_scale(gravity, dt);
    for v in body.vertices.iter_mut().filter(|v| !v.is_static()) {
        let force_impulse = v3_scale(v.force, v.inv_mass * dt);
        v.velocity = v3_scale(
            v3_add(v3_add(v.velocity, gravity_impulse), force_impulse),
            damping_factor,
        );
        v.position = v3_add(v.position, v3_scale(v.velocity, dt));
    }

    if enforce {
        let scale = body.rest_volume.abs().max(f64::MIN_POSITIVE);
        for _ in 0..MAX_VOLUME_PROJECTIONS {
            let violation = incompressibility_constraint(body);
            if violation.abs() <= VOLUME_PROJECTION_TOLERANCE * scale {
                break;
            }
        }
        let inv_dt = 1.0 / dt;
        for (v, x0) in body.vertices.iter_mut().zip(&x_start) {
            if !v.is_static() {
                v.velocity = v3_scale(v3_sub(v.position, *x0), inv_dt);
            }
        }
    }

    for v in &mut body.vertices {
        v.force = [0.0; 3];
    }
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

    // 16. Valve flow relaxes a positive gauge pressure toward zero (ambient).
    #[test]
    fn test_valve_flow_reduces_pressure() {
        let mut body = unit_box_body();
        body.pressure = 1000.0; // gauge pressure above ambient
        body.valve_coefficient = 1e-6;
        let p_before = body.pressure;
        let dp = valve_flow(&mut body, 1.0);
        assert!(dp < 0.0, "outflow must lower the pressure, got dp={dp}");
        assert!(
            body.pressure < p_before && body.pressure > 0.0,
            "gauge pressure must decrease toward (but not below) 0: {}",
            body.pressure
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

    // 19. Step body under gravity lowers vertices.
    #[test]
    fn test_step_body_gravity() {
        let mut body = PressurizedBody::unit_box(0.5, 0.0, 0.01);
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

    // 26. Valve flow with negative gauge pressure (below ambient) increases pressure.
    #[test]
    fn test_valve_flow_increases_below_ambient() {
        let mut body = unit_box_body();
        body.pressure = -500.0; // 500 Pa below ambient
        body.valve_coefficient = 1e-6;
        let p_before = body.pressure;
        let dp = valve_flow(&mut body, 1.0);
        assert!(dp > 0.0, "inflow must raise the pressure, got dp={dp}");
        assert!(
            body.pressure > p_before && body.pressure < 0.0,
            "gauge pressure must rise toward (but not above) 0: {}",
            body.pressure
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

    // 30. Step body preserves the volume (incompressibility).
    #[test]
    fn test_step_body_volume_preservation() {
        let mut body = PressurizedBody::unit_box(0.5, 0.0, 0.01);
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

    // ── valve_flow (gauge-pressure semantics) ───────────────────────────────

    // 31. Gauge pressure 0 is the fixed point of the valve.
    #[test]
    fn test_valve_flow_zero_gauge_stays_zero() {
        let mut body = unit_box_body();
        body.pressure = 0.0;
        body.valve_coefficient = 1e-3;
        let dp = valve_flow(&mut body, 1.0);
        assert_eq!(dp, 0.0);
        assert_eq!(body.pressure, 0.0, "gauge 0 must not drift toward ambient");
    }

    // 32. The pressure delta is exactly −(P_abs/V)·k_v·p·dt (isothermal, linear valve).
    #[test]
    fn test_valve_flow_matches_isothermal_formula() {
        let mut body = unit_box_body();
        body.pressure = 2000.0;
        body.valve_coefficient = 2e-6;
        let volume = compute_volume_from_verts(&body.vertices, &body.triangles);
        let expected = -((2000.0 + 101_325.0) / volume) * 2e-6 * 2000.0 * 0.5;
        let dp = valve_flow(&mut body, 0.5);
        assert!((dp - expected).abs() < 1e-9, "dp={dp}, expected={expected}");
        assert!((body.pressure - (2000.0 + expected)).abs() < 1e-9);
    }

    // 33. A smaller enclosed volume vents faster (dp ∝ 1/V).
    #[test]
    fn test_valve_flow_rate_scales_inversely_with_volume() {
        let mut small = PressurizedBody::unit_box(0.5, 1000.0, 0.01);
        let mut large = PressurizedBody::unit_box(1.0, 1000.0, 0.01);
        small.valve_coefficient = 1e-6;
        large.valve_coefficient = 1e-6;
        let dp_small = valve_flow(&mut small, 1.0);
        let dp_large = valve_flow(&mut large, 1.0);
        // Volumes differ by a factor 8.
        assert!(
            (dp_small / dp_large - 8.0).abs() < 1e-9,
            "dp_small={dp_small}, dp_large={dp_large}"
        );
    }

    // 34. A step that would overshoot is clamped at zero gauge pressure.
    #[test]
    fn test_valve_flow_never_crosses_zero() {
        let mut high = unit_box_body();
        high.pressure = 1000.0;
        high.valve_coefficient = 1e-3;
        valve_flow(&mut high, 10.0);
        assert_eq!(high.pressure, 0.0, "overshoot must be clamped at 0");

        let mut low = unit_box_body();
        low.pressure = -500.0;
        low.valve_coefficient = 1e-3;
        valve_flow(&mut low, 10.0);
        assert_eq!(low.pressure, 0.0, "overshoot must be clamped at 0");
    }

    // 35. Repeated venting converges monotonically to zero gauge pressure
    //     (the previous model drove the gauge pressure toward +101325 Pa).
    #[test]
    fn test_valve_flow_relaxes_to_zero_gauge_pressure() {
        for start in [5000.0_f64, -5000.0] {
            let mut body = unit_box_body();
            body.pressure = start;
            body.valve_coefficient = 1e-5;
            let mut previous = body.pressure.abs();
            for _ in 0..400 {
                valve_flow(&mut body, 0.1);
                let now = body.pressure.abs();
                assert!(now <= previous, "|p| must not grow: {previous} -> {now}");
                assert!(
                    body.pressure * start >= 0.0,
                    "pressure must keep its sign, got {}",
                    body.pressure
                );
                previous = now;
            }
            assert!(
                body.pressure.abs() < 1e-3,
                "gauge pressure should have relaxed to 0, got {}",
                body.pressure
            );
        }
    }

    // 36. Invalid parameters leave the pressure unchanged.
    #[test]
    fn test_valve_flow_invalid_inputs_are_noops() {
        let mut body = unit_box_body();
        body.valve_coefficient = 1e-6;
        for dt in [0.0, -1.0, f64::NAN] {
            assert_eq!(valve_flow(&mut body, dt), 0.0);
            assert_eq!(body.pressure, 1000.0);
        }
        body.valve_coefficient = -1e-6;
        assert_eq!(valve_flow(&mut body, 1.0), 0.0);
        body.valve_coefficient = f64::NAN;
        assert_eq!(valve_flow(&mut body, 1.0), 0.0);
        assert_eq!(body.pressure, 1000.0);
        body.valve_coefficient = 1e-6;
        body.pressure = f64::NAN;
        assert_eq!(valve_flow(&mut body, 1.0), 0.0);
        assert!(body.pressure.is_nan());
    }

    // 37. A mesh without enclosed volume (and rest volume) has no gas to vent.
    #[test]
    fn test_valve_flow_degenerate_volume_is_noop() {
        let flat = vec![
            SurfaceVertex::new([0.0, 0.0, 0.0], 1.0),
            SurfaceVertex::new([1.0, 0.0, 0.0], 1.0),
            SurfaceVertex::new([0.0, 1.0, 0.0], 1.0),
        ];
        let mut body = PressurizedBody::new(flat, vec![[0, 1, 2]], 500.0, 1e6, 1e-3);
        assert_eq!(valve_flow(&mut body, 1.0), 0.0);
        assert_eq!(body.pressure, 500.0);
    }

    // ── pressure force / membrane ──────────────────────────────────────────

    // 38. On a closed mesh the nodal pressure force equals p ∂V/∂xᵢ.
    #[test]
    fn test_pressure_force_equals_pressure_times_volume_gradient() {
        let mut body = unit_box_body();
        // Perturb the box so the check is not a symmetry accident.
        for (i, v) in body.vertices.iter_mut().enumerate() {
            let f = i as f64;
            v.position[0] += 0.03 * (f * 1.3).sin();
            v.position[1] += 0.03 * (f * 0.7).cos();
            v.position[2] += 0.03 * (f * 2.1).sin();
        }
        pressure_force(&mut body);
        let forces: Vec<[f64; 3]> = body.vertices.iter().map(|v| v.force).collect();
        let positions: Vec<[f64; 3]> = body.vertices.iter().map(|v| v.position).collect();
        let delta = 1e-6;
        for (i, force) in forces.iter().enumerate() {
            for axis in 0..3 {
                let mut plus = positions.clone();
                let mut minus = positions.clone();
                plus[i][axis] += delta;
                minus[i][axis] -= delta;
                let dv = (compute_volume(&plus, &body.triangles)
                    - compute_volume(&minus, &body.triangles))
                    / (2.0 * delta);
                let expected = body.pressure * dv;
                assert!(
                    (force[axis] - expected).abs() < 1e-6 * (1.0 + expected.abs()),
                    "vertex {i} axis {axis}: F={}, p dV/dx={expected}",
                    force[axis]
                );
            }
        }
    }

    // 39. A box has 18 membrane edges (12 box edges + 6 face diagonals).
    #[test]
    fn test_membrane_edges_of_box() {
        let body = unit_box_body();
        assert_eq!(body.membrane_edges.len(), 18);
        let unit = body
            .membrane_edges
            .iter()
            .filter(|e| (e.rest_length - 1.0).abs() < 1e-12)
            .count();
        let diagonal = body
            .membrane_edges
            .iter()
            .filter(|e| (e.rest_length - 2.0_f64.sqrt()).abs() < 1e-12)
            .count();
        assert_eq!((unit, diagonal), (12, 6));
        assert!(
            body.membrane_edges
                .iter()
                .all(|e| e.vertices[0] < e.vertices[1])
        );
    }

    // 40. Undeformed membrane exerts no force; a stretched one pulls back
    //     with zero net force.
    #[test]
    fn test_membrane_force_restoring_and_balanced() {
        let mut body = unit_box_body();
        body.membrane_stiffness = 500.0;
        membrane_force(&mut body);
        assert!(body.vertices.iter().all(|v| v3_norm(v.force) < 1e-9));

        // Uniform 10 % stretch: every vertex is pulled back toward the centre.
        for v in &mut body.vertices {
            v.position = v3_scale(v.position, 1.1);
        }
        membrane_force(&mut body);
        let net = body
            .vertices
            .iter()
            .fold([0.0; 3], |acc, v| v3_add(acc, v.force));
        assert!(v3_norm(net) < 1e-9, "net membrane force {net:?}");
        for v in &body.vertices {
            assert!(
                v3_dot(v.force, v.position) < 0.0,
                "force must point inward at {:?}",
                v.position
            );
        }
    }

    // 41. Zero stiffness disables the membrane.
    #[test]
    fn test_membrane_force_disabled_without_stiffness() {
        let mut body = unit_box_body();
        for v in &mut body.vertices {
            v.position = v3_scale(v.position, 1.5);
        }
        membrane_force(&mut body);
        assert!(body.vertices.iter().all(|v| v.force == [0.0; 3]));
    }

    // ── step_body: incompressible regime ───────────────────────────────────

    fn kinetic_energy(body: &PressurizedBody) -> f64 {
        body.vertices
            .iter()
            .filter(|v| !v.is_static())
            .map(|v| 0.5 * v3_dot(v.velocity, v.velocity) / v.inv_mass)
            .sum()
    }

    // 42. A pressurised, volume-constrained body stays at rest: the pressure
    //     impulse is removed by the projection *and* the velocity update.
    #[test]
    fn test_step_body_pressurised_incompressible_stays_at_rest() {
        let mut body = PressurizedBody::unit_box(0.5, 1.0e3, 0.01);
        assert!(body.enforce_incompressible);
        let rest = body.rest_volume;
        let start: Vec<[f64; 3]> = body.vertices.iter().map(|v| v.position).collect();
        for _ in 0..1000 {
            step_body(&mut body, [0.0; 3], 1.0 / 60.0);
            assert!(kinetic_energy(&body) < 1e-6, "KE={}", kinetic_energy(&body));
        }
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles);
        assert!(((vol - rest) / rest).abs() < 1e-9, "volume drifted: {vol}");
        for (v, x0) in body.vertices.iter().zip(&start) {
            assert!(v3_norm(v3_sub(v.position, *x0)) < 1e-6, "vertex moved");
        }
    }

    // 43. A body under gravity *and* pressure still falls rigidly with its
    //     volume intact.
    #[test]
    fn test_step_body_pressurised_falls_under_gravity() {
        let mut body = PressurizedBody::unit_box(0.5, 5.0e3, 0.01);
        let rest = body.rest_volume;
        for _ in 0..30 {
            step_body(&mut body, [0.0, -9.81, 0.0], 1.0 / 60.0);
        }
        let mean_y: f64 =
            body.vertices.iter().map(|v| v.position[1]).sum::<f64>() / body.vertices.len() as f64;
        // 0.5 s of free fall ≈ 1.2 m (semi-implicit Euler: slightly more).
        assert!(mean_y < -1.0 && mean_y > -1.5, "mean_y={mean_y}");
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles);
        assert!(((vol - rest) / rest).abs() < 1e-9, "volume drifted: {vol}");
    }

    // 44. PBD velocity update: the volume-changing velocity component is
    //     removed and the rest (the tangential part) is kept.
    #[test]
    fn test_step_body_velocity_update_projects_onto_constraint_tangent() {
        let mut body = PressurizedBody::unit_box(0.5, 0.0, 0.01);
        // Explosive initial velocity field: radial, 2·x (≈ 1.7 m/s per corner).
        let v0: Vec<[f64; 3]> = body
            .vertices
            .iter()
            .map(|v| v3_scale(v.position, 2.0))
            .collect();
        for (v, v_init) in body.vertices.iter_mut().zip(&v0) {
            v.velocity = *v_init;
        }
        // ∂V/∂xᵢ of the initial geometry (= pressure_force at p = 1).
        let mut probe = body.clone();
        probe.pressure = 1.0;
        pressure_force(&mut probe);
        let grad: Vec<[f64; 3]> = probe.vertices.iter().map(|v| v.force).collect();
        let along: f64 = grad.iter().zip(&v0).map(|(g, v)| v3_dot(*g, *v)).sum();
        let norm_sq: f64 = grad.iter().map(|g| v3_dot(*g, *g)).sum();
        assert!(along > 0.5, "the test velocity must expand the body");
        let mu = along / norm_sq;

        step_body(&mut body, [0.0; 3], 1.0 / 240.0);

        // Volume rate of the new velocity field vanishes ...
        let rate_after: f64 = grad
            .iter()
            .zip(&body.vertices)
            .map(|(g, v)| v3_dot(*g, v.velocity))
            .sum();
        assert!(
            rate_after.abs() < 0.02 * along,
            "dV/dt = {rate_after} (was {along})"
        );
        // ... and every vertex velocity equals the tangent-space projection
        // v − μ ∇V of the initial one (to first order in Δt).
        for (i, ((v, g), v_init)) in body.vertices.iter().zip(&grad).zip(&v0).enumerate() {
            let expected = v3_sub(*v_init, v3_scale(*g, mu));
            let err = v3_norm(v3_sub(v.velocity, expected));
            assert!(
                err < 0.05 * v3_norm(*v_init),
                "vertex {i}: v={:?}, expected={expected:?}",
                v.velocity
            );
        }
    }

    // 45. Forces are cleared on every vertex, static ones included.
    #[test]
    fn test_step_body_clears_forces_on_static_vertices() {
        for enforce in [true, false] {
            let mut body = unit_box_body();
            body.enforce_incompressible = enforce;
            body.membrane_stiffness = 100.0;
            let pinned = body.vertices[0].position;
            body.vertices[0] = SurfaceVertex::new_static(pinned);
            body.vertices[0].force = [1.0, 2.0, 3.0];
            body.vertices[3].force = [0.5, 0.0, 0.0];
            step_body(&mut body, [0.0, -9.81, 0.0], 1e-4);
            for (i, v) in body.vertices.iter().enumerate() {
                assert_eq!(
                    v.force, [0.0; 3],
                    "vertex {i} kept a force (enforce={enforce})"
                );
            }
            assert_eq!(body.vertices[0].position, pinned, "pinned vertex moved");
            assert_eq!(body.vertices[0].velocity, [0.0; 3]);
        }
    }

    // 46. Invalid time steps are ignored.
    #[test]
    fn test_step_body_invalid_dt_is_noop() {
        let mut body = unit_box_body();
        let before: Vec<[f64; 3]> = body.vertices.iter().map(|v| v.position).collect();
        for dt in [0.0, -0.01, f64::NAN] {
            step_body(&mut body, [0.0, -9.81, 0.0], dt);
            let after: Vec<[f64; 3]> = body.vertices.iter().map(|v| v.position).collect();
            assert_eq!(before, after);
            assert!(body.vertices.iter().all(|v| v.velocity == [0.0; 3]));
        }
    }

    // ── step_body: inflatable regime ───────────────────────────────────────

    /// Maximum nodal residual |p ∂V/∂xᵢ + F_membrane,i|.
    fn max_force_residual(body: &PressurizedBody) -> f64 {
        let mut probe = body.clone();
        for v in &mut probe.vertices {
            v.force = [0.0; 3];
        }
        pressure_force(&mut probe);
        membrane_force(&mut probe);
        probe
            .vertices
            .iter()
            .map(|v| v3_norm(v.force))
            .fold(0.0, f64::max)
    }

    fn inflatable_box(pressure: f64) -> PressurizedBody {
        let mut body = PressurizedBody::unit_box(0.5, pressure, 0.01);
        body.enforce_incompressible = false;
        body.membrane_stiffness = 1.0e4;
        body.damping = 200.0;
        body
    }

    // 47. With the volume free, pressure inflates the membrane to a state of
    //     force balance (pressure load = spring tension), with KE → 0.
    #[test]
    fn test_inflatable_box_inflates_to_force_balance() {
        let mut body = inflatable_box(1000.0);
        let rest = body.rest_volume;
        let load_scale = {
            let mut probe = body.clone();
            pressure_force(&mut probe);
            probe
                .vertices
                .iter()
                .map(|v| v3_norm(v.force))
                .fold(0.0, f64::max)
        };
        for _ in 0..4000 {
            step_body(&mut body, [0.0; 3], 1.0e-4);
        }
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles);
        assert!(
            vol > 1.02 * rest,
            "should have inflated: V={vol}, rest={rest}"
        );
        assert!(vol < 1.15 * rest, "inflation should be bounded: V={vol}");
        assert!(kinetic_energy(&body) < 1e-9, "KE={}", kinetic_energy(&body));
        let residual = max_force_residual(&body);
        assert!(
            residual < 1e-3 * load_scale,
            "residual {residual} vs load scale {load_scale}"
        );
    }

    // 48. A higher pressure inflates the same membrane further.
    #[test]
    fn test_inflatable_box_volume_grows_with_pressure() {
        let mut volumes = Vec::new();
        for pressure in [500.0, 1000.0, 2000.0] {
            let mut body = inflatable_box(pressure);
            for _ in 0..4000 {
                step_body(&mut body, [0.0; 3], 1.0e-4);
            }
            volumes.push(compute_volume_from_verts(&body.vertices, &body.triangles));
        }
        assert!(
            volumes[0] > 1.0 && volumes[1] > volumes[0] && volumes[2] > volumes[1],
            "volumes {volumes:?}"
        );
    }

    // 49. A stiffer membrane inflates less.
    #[test]
    fn test_inflatable_box_stiffer_membrane_inflates_less() {
        let mut volumes = Vec::new();
        for k in [5.0e3, 2.0e4] {
            let mut body = inflatable_box(1000.0);
            body.membrane_stiffness = k;
            for _ in 0..4000 {
                step_body(&mut body, [0.0; 3], 1.0e-4);
            }
            volumes.push(compute_volume_from_verts(&body.vertices, &body.triangles));
        }
        assert!(volumes[0] > volumes[1], "volumes {volumes:?}");
    }

    // 50. Without membrane stiffness damping is the only dissipation and an
    //     unpressurised free body just keeps its shape.
    #[test]
    fn test_unpressurised_inflatable_keeps_shape() {
        let mut body = inflatable_box(0.0);
        let rest = body.rest_volume;
        for _ in 0..500 {
            step_body(&mut body, [0.0; 3], 1.0e-4);
        }
        let vol = compute_volume_from_verts(&body.vertices, &body.triangles);
        assert!((vol - rest).abs() < 1e-12, "V={vol}");
    }

    // 51. Rebuilding the edge table re-captures the rest lengths.
    #[test]
    fn test_rebuild_membrane_edges_recaptures_rest_lengths() {
        let mut body = unit_box_body();
        for v in &mut body.vertices {
            v.position = v3_scale(v.position, 2.0);
        }
        body.rebuild_membrane_edges();
        body.membrane_stiffness = 100.0;
        membrane_force(&mut body);
        assert!(body.vertices.iter().all(|v| v3_norm(v.force) < 1e-9));
        assert!(body.membrane_edges.iter().all(|e| e.rest_length > 1.9));
    }
}
