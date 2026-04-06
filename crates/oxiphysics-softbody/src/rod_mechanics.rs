#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Elastic rod mechanics based on Cosserat / Kirchhoff rod theory.
//!
//! Implements stretch-shear energy, bending-twist energy, Bishop
//! (rotation-minimizing) frames, writhe integrals, and self-contact
//! repulsion forces for slender elastic rods such as DNA filaments, cables,
//! and hair strands.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// Types
// ─────────────────────────────────────────────────────────────────────────────

/// A 3-D vector stored as `[f64; 3]` (x, y, z).
pub type Vec3 = [f64; 3];

/// A 3×3 matrix stored in row-major order.
pub type Mat3 = [f64; 9];

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers
// ─────────────────────────────────────────────────────────────────────────────

fn dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(v: Vec3) -> f64 {
    dot(v, v).sqrt()
}

fn sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn normalize(v: Vec3) -> Vec3 {
    let n = norm(v);
    if n < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        scale(v, 1.0 / n)
    }
}

/// Project `v` onto the plane perpendicular to unit vector `n`, then normalize.
fn project_and_normalize(v: Vec3, n: Vec3) -> Vec3 {
    let proj = dot(v, n);
    let out = [v[0] - proj * n[0], v[1] - proj * n[1], v[2] - proj * n[2]];
    normalize(out)
}

// ─────────────────────────────────────────────────────────────────────────────
// CosseratRod
// ─────────────────────────────────────────────────────────────────────────────

/// A discretized Cosserat rod with a centerline and attached material frames.
///
/// The rod is represented as `n` nodes with positions `centerline[i]`, and
/// `n - 1` edges, each carrying an orthonormal material frame
/// `(d1[i], d2[i], d3[i])` where `d3` is the tangent direction.
#[derive(Debug, Clone)]
pub struct CosseratRod {
    /// Centerline node positions \[m\].
    pub centerline: Vec<Vec3>,
    /// First material director (normal to tangent).
    pub d1: Vec<Vec3>,
    /// Second material director (binormal direction).
    pub d2: Vec<Vec3>,
    /// Tangent director (along the rod edge).
    pub d3: Vec<Vec3>,
    /// Stretching stiffness EA \[N\].
    pub stiffness_stretch: f64,
    /// Shear stiffness GA \[N\].
    pub stiffness_shear: f64,
    /// Bending stiffness EI \[N·m²\].
    pub stiffness_bend: f64,
    /// Torsional stiffness GJ \[N·m²\].
    pub stiffness_twist: f64,
    /// Rest length of each segment \[m\].
    pub rest_lengths: Vec<f64>,
}

impl CosseratRod {
    /// Create a straight rod along the +Z axis.
    ///
    /// Nodes are placed at `z = i * segment_length` for i = 0..=n_nodes-1.
    ///
    /// # Arguments
    /// * `n_nodes` – Number of nodes (must be ≥ 2).
    /// * `segment_length` – Rest length of each segment \[m\].
    /// * `stiffness_stretch` – EA \[N\].
    /// * `stiffness_shear` – GA \[N\].
    /// * `stiffness_bend` – EI \[N·m²\].
    /// * `stiffness_twist` – GJ \[N·m²\].
    pub fn new_straight(
        n_nodes: usize,
        segment_length: f64,
        stiffness_stretch: f64,
        stiffness_shear: f64,
        stiffness_bend: f64,
        stiffness_twist: f64,
    ) -> Self {
        assert!(n_nodes >= 2, "Rod needs at least 2 nodes");
        let n_segs = n_nodes - 1;
        let centerline: Vec<Vec3> = (0..n_nodes)
            .map(|i| [0.0, 0.0, i as f64 * segment_length])
            .collect();
        let d1 = vec![[1.0_f64, 0.0, 0.0]; n_segs];
        let d2 = vec![[0.0_f64, 1.0, 0.0]; n_segs];
        let d3 = vec![[0.0_f64, 0.0, 1.0]; n_segs];
        let rest_lengths = vec![segment_length; n_segs];
        Self {
            centerline,
            d1,
            d2,
            d3,
            stiffness_stretch,
            stiffness_shear,
            stiffness_bend,
            stiffness_twist,
            rest_lengths,
        }
    }

    /// Number of centerline nodes.
    pub fn n_nodes(&self) -> usize {
        self.centerline.len()
    }

    /// Number of rod segments (edges).
    pub fn n_segments(&self) -> usize {
        self.n_nodes() - 1
    }

    /// Compute the current tangent of segment `i` (un-normalized).
    pub fn edge_tangent(&self, i: usize) -> Vec3 {
        sub(self.centerline[i + 1], self.centerline[i])
    }

    /// Compute the current length of segment `i`.
    pub fn edge_length(&self, i: usize) -> f64 {
        norm(self.edge_tangent(i))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Stretch-shear energy
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the total stretch-shear elastic energy \[J\] of the rod.
///
/// Each segment contributes:
///   E_stretch = (EA/2) * (ε)² * L₀
///   E_shear   = (GA/2) * (γ₁² + γ₂²) * L₀
///
/// where ε = (L − L₀)/L₀ is the axial strain, and γ₁, γ₂ are the shear
/// strains (components of the edge tangent not aligned with d3).
///
/// For a Kirchhoff rod (no shear) set `stiffness_shear = 0`.
///
/// # Arguments
/// * `rod` – The Cosserat rod.
pub fn stretch_shear_energy(rod: &CosseratRod) -> f64 {
    let mut energy = 0.0;
    for i in 0..rod.n_segments() {
        let l0 = rod.rest_lengths[i];
        let edge = rod.edge_tangent(i);
        let l = norm(edge);
        // Axial stretch strain.
        let eps = (l - l0) / l0;
        energy += 0.5 * rod.stiffness_stretch * eps * eps * l0;
        // Shear strains: components of (edge/L) in d1 and d2 directions.
        if rod.stiffness_shear > 0.0 && l > 1e-300 {
            let t = scale(edge, 1.0 / l);
            let gamma1 = dot(t, rod.d1[i]);
            let gamma2 = dot(t, rod.d2[i]);
            energy += 0.5 * rod.stiffness_shear * (gamma1 * gamma1 + gamma2 * gamma2) * l0;
        }
    }
    energy
}

// ─────────────────────────────────────────────────────────────────────────────
// Bending-twist energy
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the total bending and twist energy \[J\] of the rod (Kirchhoff model).
///
/// At each interior node the discrete curvature vector κ is approximated via
/// the Darboux vector between adjacent frames:
///
///   κ ≈ (d3\[i\] × d3\[i-1\]) / (L_avg)
///
/// E_bend = (EI/2) * (κ₁² + κ₂²) * L_avg  (bending about d1, d2)
/// E_twist = (GJ/2) * κ₃² * L_avg           (twist about d3)
///
/// # Arguments
/// * `rod` – The Cosserat rod.
pub fn bending_twist_energy(rod: &CosseratRod) -> f64 {
    let mut energy = 0.0;
    for i in 1..rod.n_segments() {
        let l_avg = 0.5 * (rod.rest_lengths[i - 1] + rod.rest_lengths[i]);
        let t_prev = normalize(rod.edge_tangent(i - 1));
        let t_curr = normalize(rod.edge_tangent(i));
        // Discrete curvature from tangent difference.
        let kappa_vec = cross(t_prev, t_curr);
        // Resolve into material frame components.
        let kappa1 = dot(kappa_vec, rod.d1[i]);
        let kappa2 = dot(kappa_vec, rod.d2[i]);
        // Twist: rotation of d1 frame around d3.
        let d1_prev_in_plane = project_and_normalize(rod.d1[i - 1], t_curr);
        let twist_cos = dot(d1_prev_in_plane, rod.d1[i]).clamp(-1.0, 1.0);
        let twist_sin_vec = cross(d1_prev_in_plane, rod.d1[i]);
        let twist_sign = if dot(twist_sin_vec, t_curr) >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let kappa3 = twist_sign * twist_cos.acos() / l_avg;
        energy += 0.5 * rod.stiffness_bend * (kappa1 * kappa1 + kappa2 * kappa2) * l_avg;
        energy += 0.5 * rod.stiffness_twist * kappa3 * kappa3 * l_avg;
    }
    energy
}

// ─────────────────────────────────────────────────────────────────────────────
// Material frame transport (parallel transport)
// ─────────────────────────────────────────────────────────────────────────────

/// Parallel-transport a vector `v` from tangent `t0` to tangent `t1`.
///
/// Uses the double-reflection (Bishop) formula: one reflection maps `t0` → `t1`,
/// the second reflection brings the transported vector back into the plane.
///
/// # Arguments
/// * `v` – Vector attached to the frame at `t0`.
/// * `t0` – Source tangent (unit vector).
/// * `t1` – Target tangent (unit vector).
pub fn material_frame_transport(v: Vec3, t0: Vec3, t1: Vec3) -> Vec3 {
    // b = t0 × t1 (rotation axis).
    let b = cross(t0, t1);
    let sin_theta = norm(b);
    let cos_theta = dot(t0, t1).clamp(-1.0, 1.0);
    if sin_theta < 1e-12 {
        // Nearly parallel tangents: no rotation needed.
        return v;
    }
    let b_hat = scale(b, 1.0 / sin_theta);
    // Rodrigues' rotation formula: v' = v*cos + (b̂ × v)*sin + b̂*(b̂·v)*(1-cos)
    let bdotv = dot(b_hat, v);
    let bxv = cross(b_hat, v);
    [
        v[0] * cos_theta + bxv[0] * sin_theta + b_hat[0] * bdotv * (1.0 - cos_theta),
        v[1] * cos_theta + bxv[1] * sin_theta + b_hat[1] * bdotv * (1.0 - cos_theta),
        v[2] * cos_theta + bxv[2] * sin_theta + b_hat[2] * bdotv * (1.0 - cos_theta),
    ]
}

// ─────────────────────────────────────────────────────────────────────────────
// Bishop (rotation-minimizing) frame
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Bishop (rotation-minimizing) frame along the rod.
///
/// Returns `(u, v, t)` triples for each node, where `t` is the unit tangent,
/// and `(u, v)` are the two normal directors propagated by parallel transport
/// with no twist.
///
/// # Arguments
/// * `centerline` – Sequence of node positions (at least 2).
/// * `u0` – Initial normal director at node 0 (must be perpendicular to the
///   first tangent).
///
/// Returns a `Vec` of `(u, v, t)` tuples, one per node.
pub fn bishop_frame(centerline: &[Vec3], u0: Vec3) -> Vec<(Vec3, Vec3, Vec3)> {
    assert!(centerline.len() >= 2);
    let n = centerline.len();
    let mut frames: Vec<(Vec3, Vec3, Vec3)> = Vec::with_capacity(n);

    // First tangent.
    let t0 = normalize(sub(centerline[1], centerline[0]));
    // Ensure u0 is perpendicular to t0.
    let u0 = project_and_normalize(u0, t0);
    let v0 = cross(t0, u0);
    frames.push((u0, v0, t0));

    for i in 1..n {
        let t_prev = frames[i - 1].2;
        let t_curr = if i + 1 < n {
            normalize(sub(centerline[i + 1], centerline[i]))
        } else {
            normalize(sub(centerline[i], centerline[i - 1]))
        };
        let u_prev = frames[i - 1].0;
        let u_curr = material_frame_transport(u_prev, t_prev, t_curr);
        let v_curr = cross(t_curr, u_curr);
        frames.push((u_curr, v_curr, t_curr));
    }
    frames
}

// ─────────────────────────────────────────────────────────────────────────────
// Writhe number
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the writhe number Wr of the rod centerline using the
/// Gauss double-integral approximation (discrete version).
///
/// Wr = (1/4π) ∑_{i≠j} (r_ij × r_ij') · (dr_i × dr_j) / |r_ij|³
///
/// For efficiency the sum is computed only for non-adjacent segments.
/// Returns the signed writhe (positive = right-handed crossings).
///
/// # Arguments
/// * `centerline` – Ordered sequence of node positions (at least 3).
pub fn writhe_number(centerline: &[Vec3]) -> f64 {
    let n = centerline.len();
    if n < 3 {
        return 0.0;
    }
    let mut wr = 0.0;
    for i in 0..n - 1 {
        let dr_i = sub(centerline[i + 1], centerline[i]);
        let mid_i = scale(add(centerline[i], centerline[i + 1]), 0.5);
        for j in (i + 2)..n - 1 {
            let dr_j = sub(centerline[j + 1], centerline[j]);
            let mid_j = scale(add(centerline[j], centerline[j + 1]), 0.5);
            let r = sub(mid_i, mid_j);
            let r2 = dot(r, r);
            if r2 < 1e-30 {
                continue;
            }
            let r3 = r2 * r2.sqrt();
            let cross_dr = cross(dr_i, dr_j);
            wr += dot(cross_dr, r) / r3;
        }
    }
    wr / (4.0 * PI)
}

// ─────────────────────────────────────────────────────────────────────────────
// Self-contact detection and repulsion
// ─────────────────────────────────────────────────────────────────────────────

/// Compute pairwise self-contact repulsion forces on the rod centerline.
///
/// For each pair of non-adjacent nodes closer than `contact_radius`, a
/// soft repulsion force is applied:
///   F = k * (contact_radius - d) * r̂ / d
///
/// Returns a `Vec` of force vectors, one per node.
///
/// # Arguments
/// * `centerline` – Node positions.
/// * `contact_radius` – Radius at which repulsion activates \[m\].
/// * `stiffness` – Repulsion stiffness \[N/m\].
pub fn rod_contact_force(centerline: &[Vec3], contact_radius: f64, stiffness: f64) -> Vec<Vec3> {
    let n = centerline.len();
    let mut forces = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        for j in (i + 2)..n {
            let r = sub(centerline[j], centerline[i]);
            let d = norm(r);
            if d < contact_radius && d > 1e-300 {
                let penetration = contact_radius - d;
                let f_mag = stiffness * penetration;
                let r_hat = scale(r, 1.0 / d);
                // Push apart along the connecting vector.
                forces[i] = [
                    forces[i][0] - f_mag * r_hat[0],
                    forces[i][1] - f_mag * r_hat[1],
                    forces[i][2] - f_mag * r_hat[2],
                ];
                forces[j] = [
                    forces[j][0] + f_mag * r_hat[0],
                    forces[j][1] + f_mag * r_hat[1],
                    forces[j][2] + f_mag * r_hat[2],
                ];
            }
        }
    }
    forces
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. Straight rod at rest has zero stretch energy.
    #[test]
    fn test_stretch_energy_rest() {
        let rod = CosseratRod::new_straight(5, 1.0, 1000.0, 0.0, 100.0, 50.0);
        let e = stretch_shear_energy(&rod);
        assert!(
            e.abs() < EPS,
            "Straight rod at rest should have zero stretch energy, got {e}"
        );
    }

    // 2. Stretch energy is positive after elongation.
    #[test]
    fn test_stretch_energy_positive_after_elongation() {
        let mut rod = CosseratRod::new_straight(3, 1.0, 1000.0, 0.0, 100.0, 50.0);
        // Move the last node to stretch the rod.
        rod.centerline[2] = [0.0, 0.0, 2.5]; // 2.5 instead of 2.0
        let e = stretch_shear_energy(&rod);
        assert!(e > 0.0, "Stretched rod must have positive energy, got {e}");
    }

    // 3. Stretch energy scales as k * ε² / 2.
    #[test]
    fn test_stretch_energy_formula() {
        let ea = 1000.0;
        let l0 = 1.0;
        let mut rod = CosseratRod::new_straight(2, l0, ea, 0.0, 100.0, 50.0);
        let stretch = 0.1; // 10% strain
        rod.centerline[1] = [0.0, 0.0, l0 * (1.0 + stretch)];
        let e = stretch_shear_energy(&rod);
        let expected = 0.5 * ea * stretch * stretch * l0;
        assert!(
            (e - expected).abs() < EPS,
            "Stretch energy formula mismatch: {e} vs {expected}"
        );
    }

    // 4. Straight rod has zero bending energy.
    #[test]
    fn test_bending_energy_straight_rod() {
        let rod = CosseratRod::new_straight(6, 0.5, 1000.0, 0.0, 100.0, 50.0);
        let e = bending_twist_energy(&rod);
        assert!(
            e.abs() < EPS,
            "Straight rod must have zero bending energy, got {e}"
        );
    }

    // 5. Bent rod has positive bending energy.
    #[test]
    fn test_bending_energy_bent_rod() {
        let mut rod = CosseratRod::new_straight(4, 1.0, 1000.0, 0.0, 100.0, 50.0);
        // Introduce a kink at node 2.
        rod.centerline[2] = [0.5, 0.0, 2.0];
        rod.centerline[3] = [0.5, 0.0, 3.0];
        // Update d3 for segment 1 (now angled).
        rod.d3[1] = normalize(sub(rod.centerline[2], rod.centerline[1]));
        let e = bending_twist_energy(&rod);
        assert!(
            e > 0.0,
            "Bent rod must have positive bending energy, got {e}"
        );
    }

    // 6. Rod has correct number of nodes after construction.
    #[test]
    fn test_rod_n_nodes() {
        let rod = CosseratRod::new_straight(7, 0.3, 500.0, 0.0, 50.0, 25.0);
        assert_eq!(rod.n_nodes(), 7);
        assert_eq!(rod.n_segments(), 6);
    }

    // 7. Edge length of straight rod matches segment_length.
    #[test]
    fn test_edge_length_straight() {
        let rod = CosseratRod::new_straight(4, 2.0, 100.0, 0.0, 10.0, 5.0);
        for i in 0..rod.n_segments() {
            let l = rod.edge_length(i);
            assert!(
                (l - 2.0).abs() < EPS,
                "Edge {i} length should be 2.0, got {l}"
            );
        }
    }

    // 8. Parallel transport preserves vector length.
    #[test]
    fn test_parallel_transport_length() {
        let t0 = normalize([1.0, 0.0, 0.0]);
        let t1 = normalize([0.0, 1.0, 0.0]);
        let v = [0.0, 0.0, 1.0_f64];
        let vt = material_frame_transport(v, t0, t1);
        let len_in = norm(v);
        let len_out = norm(vt);
        assert!(
            (len_out - len_in).abs() < EPS,
            "Parallel transport must preserve vector length"
        );
    }

    // 9. Parallel transport with identical tangents returns unchanged vector.
    #[test]
    fn test_parallel_transport_identity() {
        let t = normalize([1.0, 1.0, 0.0]);
        let v = [0.0, 0.0, 1.0_f64];
        let vt = material_frame_transport(v, t, t);
        for k in 0..3 {
            assert!(
                (vt[k] - v[k]).abs() < EPS,
                "Transport with same tangent should be identity"
            );
        }
    }

    // 10. Parallel transport maps t0 to t1 correctly.
    #[test]
    fn test_parallel_transport_tangent_maps() {
        let t0 = [1.0_f64, 0.0, 0.0];
        let t1 = [0.0_f64, 0.0, 1.0];
        // Transporting t0 itself should give t1 (since it's in the rotation plane).
        let vt = material_frame_transport(t0, t0, t1);
        for k in 0..3 {
            assert!(
                (vt[k] - t1[k]).abs() < EPS,
                "Transported tangent mismatch at component {k}"
            );
        }
    }

    // 11. Bishop frame has unit tangents.
    #[test]
    fn test_bishop_frame_unit_tangents() {
        let cl: Vec<Vec3> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let frames = bishop_frame(&cl, [0.0, 1.0, 0.0]);
        for (i, (_u, _v, t)) in frames.iter().enumerate() {
            let n = norm(*t);
            assert!(
                (n - 1.0).abs() < EPS,
                "Bishop frame tangent {i} must be unit, got {n}"
            );
        }
    }

    // 12. Bishop frame normal is perpendicular to tangent.
    #[test]
    fn test_bishop_frame_orthogonality() {
        let cl: Vec<Vec3> = (0..5).map(|i| [i as f64, 0.0, 0.0]).collect();
        let frames = bishop_frame(&cl, [0.0, 1.0, 0.0]);
        for (i, (u, _v, t)) in frames.iter().enumerate() {
            let d = dot(*u, *t).abs();
            assert!(
                d < EPS,
                "Bishop u must be perpendicular to t at node {i}, dot={d}"
            );
        }
    }

    // 13. Bishop frame returns correct number of frames.
    #[test]
    fn test_bishop_frame_count() {
        let cl: Vec<Vec3> = (0..8).map(|i| [0.0, 0.0, i as f64]).collect();
        let frames = bishop_frame(&cl, [1.0, 0.0, 0.0]);
        assert_eq!(frames.len(), 8);
    }

    // 14. Writhe of straight rod is zero.
    #[test]
    fn test_writhe_straight_zero() {
        let cl: Vec<Vec3> = (0..10).map(|i| [0.0, 0.0, i as f64]).collect();
        let wr = writhe_number(&cl);
        assert!(
            wr.abs() < EPS,
            "Writhe of straight rod must be zero, got {wr}"
        );
    }

    // 15. Writhe of planar closed loop is zero.
    #[test]
    fn test_writhe_planar_loop_zero() {
        // Create a planar ring; writhe = 0 for planar curves.
        let n = 20;
        let cl: Vec<Vec3> = (0..n)
            .map(|i| {
                let theta = 2.0 * PI * i as f64 / n as f64;
                [theta.cos(), theta.sin(), 0.0]
            })
            .collect();
        let wr = writhe_number(&cl);
        assert!(wr.abs() < 0.1, "Planar loop writhe should be ~0, got {wr}");
    }

    // 16. Writhe with fewer than 3 nodes returns zero.
    #[test]
    fn test_writhe_too_short() {
        let cl: Vec<Vec3> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let wr = writhe_number(&cl);
        assert_eq!(wr, 0.0);
    }

    // 17. Rod contact force: no force when nodes are far apart.
    #[test]
    fn test_contact_force_no_force_far() {
        let cl: Vec<Vec3> = (0..5).map(|i| [0.0, 0.0, i as f64 * 2.0]).collect();
        let forces = rod_contact_force(&cl, 0.5, 1000.0);
        for (i, f) in forces.iter().enumerate() {
            let fn_ = norm(*f);
            assert!(
                fn_ < EPS,
                "No contact force expected for far nodes at {i}, got {fn_}"
            );
        }
    }

    // 18. Rod contact force: force is generated when nodes overlap.
    #[test]
    fn test_contact_force_repulsion_generated() {
        // Two non-adjacent nodes at the same position.
        let cl: Vec<Vec3> = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1e-6, 0.0, 0.0]];
        let forces = rod_contact_force(&cl, 0.1, 1000.0);
        // Nodes 0 and 2 are within contact radius.
        let f_total = norm(forces[0]);
        assert!(
            f_total > 0.0,
            "Contact force should be non-zero for overlapping nodes"
        );
    }

    // 19. Contact forces are antisymmetric (Newton's third law).
    #[test]
    fn test_contact_force_newton_third_law() {
        let cl: Vec<Vec3> = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.05, 0.0, 0.0]];
        let forces = rod_contact_force(&cl, 0.1, 500.0);
        for k in 0..3 {
            assert!(
                (forces[0][k] + forces[2][k]).abs() < EPS,
                "Forces must be antisymmetric (Newton III) at component {k}"
            );
        }
    }

    // 20. Stretch energy is symmetric: equal stretch and compression give equal energy.
    #[test]
    fn test_stretch_energy_symmetry() {
        let ea = 500.0;
        let l0 = 1.0;
        let eps_val = 0.2;

        let mut rod_stretch = CosseratRod::new_straight(2, l0, ea, 0.0, 10.0, 5.0);
        rod_stretch.centerline[1] = [0.0, 0.0, l0 * (1.0 + eps_val)];
        let e_stretch = stretch_shear_energy(&rod_stretch);

        let mut rod_compress = CosseratRod::new_straight(2, l0, ea, 0.0, 10.0, 5.0);
        rod_compress.centerline[1] = [0.0, 0.0, l0 * (1.0 - eps_val)];
        let e_compress = stretch_shear_energy(&rod_compress);

        assert!(
            (e_stretch - e_compress).abs() < EPS,
            "Stretch energy must be symmetric: {e_stretch} vs {e_compress}"
        );
    }

    // 21. Bending energy increases with larger bend angle.
    #[test]
    fn test_bending_energy_increases_with_bend() {
        let make_bent_rod = |angle: f64| {
            let mut rod = CosseratRod::new_straight(4, 1.0, 1000.0, 0.0, 100.0, 50.0);
            rod.centerline[2] = [angle.sin(), 0.0, 2.0 * angle.cos()];
            rod.centerline[3] = [2.0 * angle.sin(), 0.0, 3.0 * angle.cos()];
            rod.d3[1] = normalize(sub(rod.centerline[2], rod.centerline[1]));
            rod
        };
        let e_small = bending_twist_energy(&make_bent_rod(0.1));
        let e_large = bending_twist_energy(&make_bent_rod(0.5));
        assert!(
            e_large > e_small,
            "More bent rod must have higher bending energy"
        );
    }

    // 22. Bishop frame v = t × u (right-handed frame check).
    #[test]
    fn test_bishop_frame_right_handed() {
        let cl: Vec<Vec3> = (0..4).map(|i| [i as f64, 0.0, 0.0]).collect();
        let frames = bishop_frame(&cl, [0.0, 1.0, 0.0]);
        for (i, (u, v, t)) in frames.iter().enumerate() {
            let v_expected = cross(*t, *u);
            for k in 0..3 {
                assert!(
                    (v[k] - v_expected[k]).abs() < EPS,
                    "Frame {i} must be right-handed at component {k}"
                );
            }
        }
    }

    // 23. Rod rest_lengths vector has correct size.
    #[test]
    fn test_rest_lengths_size() {
        let rod = CosseratRod::new_straight(10, 0.5, 200.0, 0.0, 20.0, 10.0);
        assert_eq!(
            rod.rest_lengths.len(),
            9,
            "rest_lengths should have n_nodes - 1 entries"
        );
    }

    // 24. Contact force magnitude scales with stiffness.
    #[test]
    fn test_contact_force_scales_with_stiffness() {
        let cl: Vec<Vec3> = vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.05, 0.0, 0.0]];
        let f1 = rod_contact_force(&cl, 0.1, 100.0);
        let f2 = rod_contact_force(&cl, 0.1, 200.0);
        let n1 = norm(f1[0]);
        let n2 = norm(f2[0]);
        assert!(
            (n2 / n1 - 2.0).abs() < EPS,
            "Contact force must scale with stiffness: {n1}, {n2}"
        );
    }

    // 25. Bishop frame for a straight rod always produces the same u and v.
    #[test]
    fn test_bishop_frame_straight_constant() {
        let cl: Vec<Vec3> = (0..6).map(|i| [0.0, 0.0, i as f64]).collect();
        let u0 = [1.0_f64, 0.0, 0.0];
        let frames = bishop_frame(&cl, u0);
        for (i, (u, _v, _t)) in frames.iter().enumerate() {
            // Along +Z, the u director should remain [1, 0, 0].
            assert!(
                (u[0] - 1.0).abs() < EPS && u[1].abs() < EPS && u[2].abs() < EPS,
                "Bishop u at node {i} should be [1,0,0], got {:?}",
                u
            );
        }
    }
}
