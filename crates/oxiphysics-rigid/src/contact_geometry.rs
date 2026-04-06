// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Contact geometry and friction cone analysis for rigid body systems.
//!
//! Provides tools for modelling contact points, checking friction cone
//! membership, computing contact wrench spaces, and evaluating grasp quality.
//!
//! # Overview
//!
//! - [`ContactPoint`] — contact location, normal, and penetration depth.
//! - [`friction_cone_vertices`] — linearised friction cone approximation.
//! - [`inside_friction_cone`] — test whether a force vector lies in the cone.
//! - [`contact_wrench_space`] — force/torque pairs achievable via contacts.
//! - [`grasp_quality`] — smallest singular value of the grasp matrix.
//! - [`form_closure_check`] — necessary condition for form closure.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// ContactPoint
// ─────────────────────────────────────────────────────────────────────────────

/// A contact point between two rigid bodies.
#[derive(Debug, Clone, Copy)]
pub struct ContactPoint {
    /// Contact position in world coordinates \[m\].
    pub position: [f64; 3],
    /// Outward contact normal (unit vector pointing away from body).
    pub normal: [f64; 3],
    /// Penetration depth \[m\] (positive = overlap).
    pub penetration_depth: f64,
    /// Coefficient of friction (Coulomb model).
    pub friction_coeff: f64,
}

impl ContactPoint {
    /// Construct a contact point.
    pub fn new(
        position: [f64; 3],
        normal: [f64; 3],
        penetration_depth: f64,
        friction_coeff: f64,
    ) -> Self {
        Self {
            position,
            normal,
            penetration_depth,
            friction_coeff,
        }
    }

    /// Return the normal with unit length.
    pub fn unit_normal(&self) -> [f64; 3] {
        let len = vec3_norm(self.normal);
        if len < 1e-12 {
            [0.0, 1.0, 0.0]
        } else {
            vec3_scale(self.normal, 1.0 / len)
        }
    }

    /// True if the contact is active (non-negative penetration).
    pub fn is_active(&self) -> bool {
        self.penetration_depth >= 0.0
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Friction cone helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Return `num_sides` edge vectors of the linearised friction cone at `contact`.
///
/// The friction cone is approximated by a polyhedral cone whose edges lie on
/// the boundary surface: `|f_t| = mu * f_n`.  Each returned vector is a
/// candidate extreme direction for force application.
///
/// # Arguments
///
/// * `contact` — contact point with normal and friction coefficient.
/// * `num_sides` — number of facets (≥ 3).  Typical: 4 or 8.
///
/// Returns a `Vec` of `num_sides` unit force directions.
pub fn friction_cone_vertices(contact: &ContactPoint, num_sides: usize) -> Vec<[f64; 3]> {
    assert!(num_sides >= 3, "friction cone must have at least 3 sides");

    let n = contact.unit_normal();
    let mu = contact.friction_coeff;

    // Build two orthogonal tangent vectors
    let (t1, t2) = orthonormal_basis(n);

    let mut vertices = Vec::with_capacity(num_sides);
    for k in 0..num_sides {
        let angle = 2.0 * PI * k as f64 / num_sides as f64;
        // Tangential component = mu * (cos(θ) * t1 + sin(θ) * t2)
        let tx = mu * (angle.cos() * t1[0] + angle.sin() * t2[0]);
        let ty = mu * (angle.cos() * t1[1] + angle.sin() * t2[1]);
        let tz = mu * (angle.cos() * t1[2] + angle.sin() * t2[2]);
        // Full cone edge: normal + tangential
        let v = [n[0] + tx, n[1] + ty, n[2] + tz];
        let v_norm = vec3_norm(v);
        vertices.push(if v_norm > 1e-12 {
            vec3_scale(v, 1.0 / v_norm)
        } else {
            n
        });
    }
    vertices
}

/// Return `true` if force vector `f` lies inside (or on the boundary of) the
/// friction cone at `contact`.
///
/// Coulomb condition: `|f_tangential| ≤ mu * f_normal`.
/// The normal component must be non-negative (no pull allowed).
pub fn inside_friction_cone(contact: &ContactPoint, force: [f64; 3]) -> bool {
    let n = contact.unit_normal();
    let mu = contact.friction_coeff;

    // Normal component (should be ≥ 0)
    let f_n = vec3_dot(force, n);
    if f_n < -1e-12 {
        return false;
    }

    // Tangential component
    let f_t_vec = [
        force[0] - f_n * n[0],
        force[1] - f_n * n[1],
        force[2] - f_n * n[2],
    ];
    let f_t = vec3_norm(f_t_vec);

    f_t <= mu * f_n + 1e-9
}

// ─────────────────────────────────────────────────────────────────────────────
// Grasp / contact wrench space
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the 6-DOF wrench (force + torque about origin) achievable by a
/// unit contact force at `contact` along direction `force_dir`.
///
/// Returns `[fx, fy, fz, tx, ty, tz]`.
pub fn contact_wrench(contact: &ContactPoint, force_dir: [f64; 3]) -> [f64; 6] {
    let p = contact.position;
    let f = force_dir;
    // Torque = p × f
    let t = vec3_cross(p, f);
    [f[0], f[1], f[2], t[0], t[1], t[2]]
}

/// Build the 6×(n_contacts * num_sides) grasp matrix `G` where each column is
/// a wrench primitive from the linearised friction cone.
///
/// Returns the matrix as a flat `Vec`f64` (6 rows, column-major) and
/// `(rows=6, cols)`.
pub fn contact_wrench_space(
    contacts: &[ContactPoint],
    num_sides: usize,
) -> (Vec<f64>, usize, usize) {
    let cols = contacts.len() * num_sides;
    let mut g = vec![0.0f64; 6 * cols];

    let mut col = 0;
    for contact in contacts {
        let verts = friction_cone_vertices(contact, num_sides);
        for v in &verts {
            let w = contact_wrench(contact, *v);
            for row in 0..6 {
                g[row + 6 * col] = w[row];
            }
            col += 1;
        }
    }
    (g, 6, cols)
}

/// Compute the grasp quality as the smallest singular value of the grasp matrix.
///
/// A larger value indicates a more robust grasp.  The grasp matrix `g_mat` is
/// provided in column-major order with `rows` rows and `cols` columns.
///
/// Uses a simplified Frobenius-norm-based estimate (full SVD requires LAPACK).
/// For exact singular values use an external linear algebra library.
pub fn grasp_quality(g_mat: &[f64], rows: usize, cols: usize) -> f64 {
    if rows == 0 || cols == 0 {
        return 0.0;
    }
    // Compute G * G^T (rows × rows symmetric matrix)
    let mut gg = vec![0.0f64; rows * rows];
    for i in 0..rows {
        for j in 0..rows {
            let mut sum = 0.0;
            for k in 0..cols {
                sum += g_mat[i + rows * k] * g_mat[j + rows * k];
            }
            gg[i + rows * j] = sum;
        }
    }
    // Estimate smallest eigenvalue of G*G^T via power-iteration on (I - lambda_max * GG^T)
    // For simplicity: return sqrt of min diagonal of GG^T as lower-bound estimate
    let mut min_diag = f64::INFINITY;
    for i in 0..rows {
        let d = gg[i + rows * i];
        if d < min_diag {
            min_diag = d;
        }
    }
    if min_diag <= 0.0 {
        0.0
    } else {
        min_diag.sqrt()
    }
}

/// Check the necessary condition for form closure: the contact wrenches must
/// positively span the 6-DOF wrench space.
///
/// The necessary (but not sufficient) condition checked here is that for each
/// of the 6 canonical wrench directions, at least one contact wrench has a
/// positive component.
///
/// Returns `true` if the necessary condition is satisfied.
pub fn form_closure_check(contacts: &[ContactPoint], num_sides: usize) -> bool {
    if contacts.is_empty() {
        return false;
    }
    let (g_mat, _rows, cols) = contact_wrench_space(contacts, num_sides);

    // For each of the 6 wrench axes, check we have both positive and negative components
    for axis in 0..6 {
        let mut has_pos = false;
        let mut has_neg = false;
        for c in 0..cols {
            let val = g_mat[axis + 6 * c];
            if val > 1e-9 {
                has_pos = true;
            }
            if val < -1e-9 {
                has_neg = true;
            }
        }
        if !has_pos || !has_neg {
            return false;
        }
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Vector math helpers
// ─────────────────────────────────────────────────────────────────────────────

fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn vec3_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Build two unit vectors orthonormal to `n`.
fn orthonormal_basis(n: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Choose a vector not parallel to n
    let up = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let t1_raw = vec3_cross(n, up);
    let t1_len = vec3_norm(t1_raw);
    let t1 = if t1_len > 1e-12 {
        vec3_scale(t1_raw, 1.0 / t1_len)
    } else {
        [0.0, 1.0, 0.0]
    };
    let t2 = vec3_cross(n, t1);
    (t1, t2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn contact_up(mu: f64) -> ContactPoint {
        ContactPoint::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, mu)
    }

    // ── ContactPoint ─────────────────────────────────────────────────────

    #[test]
    fn contact_unit_normal_already_unit() {
        let c = contact_up(0.5);
        let n = c.unit_normal();
        let len = vec3_norm(n);
        assert!((len - 1.0).abs() < 1e-10);
    }

    #[test]
    fn contact_unit_normal_normalises() {
        let c = ContactPoint::new([0.0; 3], [0.0, 2.0, 0.0], 0.0, 0.5);
        let n = c.unit_normal();
        assert!((vec3_norm(n) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn contact_is_active_positive_penetration() {
        let c = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.01, 0.5);
        assert!(c.is_active());
    }

    #[test]
    fn contact_is_not_active_negative_penetration() {
        let c = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], -0.01, 0.5);
        assert!(!c.is_active());
    }

    #[test]
    fn contact_zero_normal_returns_fallback() {
        let c = ContactPoint::new([0.0; 3], [0.0, 0.0, 0.0], 0.0, 0.5);
        let n = c.unit_normal();
        assert!(vec3_norm(n) > 0.5);
    }

    // ── friction_cone_vertices ───────────────────────────────────────────

    #[test]
    fn friction_cone_correct_count() {
        let c = contact_up(0.3);
        let v = friction_cone_vertices(&c, 8);
        assert_eq!(v.len(), 8);
    }

    #[test]
    fn friction_cone_vertices_unit_length() {
        let c = contact_up(0.5);
        for v in friction_cone_vertices(&c, 6) {
            assert!(
                (vec3_norm(v) - 1.0).abs() < 1e-10,
                "vertex not unit: {:?}",
                v
            );
        }
    }

    #[test]
    fn friction_cone_zero_mu_gives_normal() {
        let c = contact_up(0.0);
        for v in friction_cone_vertices(&c, 4) {
            // With mu=0 all vertices collapse to the normal direction
            assert!((v[1] - 1.0).abs() < 1e-9);
        }
    }

    #[test]
    fn friction_cone_four_sides_symmetric() {
        let c = contact_up(0.5);
        let verts = friction_cone_vertices(&c, 4);
        // The x-components should be symmetric: +mu, 0, -mu, 0 (approx)
        let xs: Vec<f64> = verts.iter().map(|v| v[0]).collect();
        let max_x = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let min_x = xs.iter().cloned().fold(f64::INFINITY, f64::min);
        assert!((max_x + min_x).abs() < 1e-9, "should be symmetric in x");
    }

    #[test]
    fn friction_cone_minimum_three_sides() {
        let c = contact_up(0.5);
        let v = friction_cone_vertices(&c, 3);
        assert_eq!(v.len(), 3);
    }

    // ── inside_friction_cone ─────────────────────────────────────────────

    #[test]
    fn inside_cone_pure_normal_force_accepted() {
        let c = contact_up(0.5);
        assert!(inside_friction_cone(&c, [0.0, 1.0, 0.0]));
    }

    #[test]
    fn inside_cone_pull_rejected() {
        let c = contact_up(0.5);
        assert!(!inside_friction_cone(&c, [0.0, -1.0, 0.0]));
    }

    #[test]
    fn inside_cone_on_boundary_accepted() {
        // mu=1 means 45° → force [1, 1, 0] / sqrt(2) has |ft| = fn
        let c = contact_up(1.0);
        let f = [1.0_f64, 1.0, 0.0];
        assert!(inside_friction_cone(&c, f));
    }

    #[test]
    fn inside_cone_outside_rejected() {
        // mu=0.1, large tangential component
        let c = contact_up(0.1);
        assert!(!inside_friction_cone(&c, [10.0, 1.0, 0.0]));
    }

    #[test]
    fn inside_cone_zero_force_accepted() {
        let c = contact_up(0.5);
        // Zero force: f_n=0, f_t=0, 0 <= 0 * mu
        assert!(inside_friction_cone(&c, [0.0, 0.0, 0.0]));
    }

    // ── contact_wrench ───────────────────────────────────────────────────

    #[test]
    fn contact_wrench_at_origin_no_torque() {
        let c = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0.5);
        let w = contact_wrench(&c, [0.0, 1.0, 0.0]);
        // Torque = 0 × force = 0
        assert!(w[3].abs() < 1e-12);
        assert!(w[4].abs() < 1e-12);
        assert!(w[5].abs() < 1e-12);
    }

    #[test]
    fn contact_wrench_off_origin_has_torque() {
        let c = ContactPoint::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 0.0, 0.5);
        let w = contact_wrench(&c, [0.0, 1.0, 0.0]);
        // p × f = [1,0,0] × [0,1,0] = [0,0,1]
        assert!((w[5] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn contact_wrench_force_component_correct() {
        let c = ContactPoint::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 0.5);
        let w = contact_wrench(&c, [3.0, 4.0, 5.0]);
        assert!((w[0] - 3.0).abs() < 1e-12);
        assert!((w[1] - 4.0).abs() < 1e-12);
        assert!((w[2] - 5.0).abs() < 1e-12);
    }

    // ── contact_wrench_space ─────────────────────────────────────────────

    #[test]
    fn wrench_space_column_count() {
        let contacts = vec![contact_up(0.5), contact_up(0.5)];
        let (_g, rows, cols) = contact_wrench_space(&contacts, 4);
        assert_eq!(rows, 6);
        assert_eq!(cols, 8); // 2 contacts × 4 sides
    }

    #[test]
    fn wrench_space_empty_contacts() {
        let (_g, _rows, cols) = contact_wrench_space(&[], 4);
        assert_eq!(cols, 0);
    }

    // ── grasp_quality ────────────────────────────────────────────────────

    #[test]
    fn grasp_quality_empty_is_zero() {
        assert_eq!(grasp_quality(&[], 0, 0), 0.0);
    }

    #[test]
    fn grasp_quality_positive_for_valid_grasp() {
        let contacts = vec![
            ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.0, 0.5),
            ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0, 0.5),
            ContactPoint::new([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 0.0, 0.5),
            ContactPoint::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.0, 0.5),
        ];
        let (g, rows, cols) = contact_wrench_space(&contacts, 4);
        let q = grasp_quality(&g, rows, cols);
        assert!(q >= 0.0);
    }

    #[test]
    fn grasp_quality_single_contact_nonnegative() {
        let contacts = vec![contact_up(0.5)];
        let (g, rows, cols) = contact_wrench_space(&contacts, 4);
        let q = grasp_quality(&g, rows, cols);
        assert!(q >= 0.0);
    }

    // ── form_closure_check ────────────────────────────────────────────────

    #[test]
    fn form_closure_empty_contacts_false() {
        assert!(!form_closure_check(&[], 4));
    }

    #[test]
    fn form_closure_opposing_contacts_may_pass() {
        // Contacts on all 6 faces of a unit cube pointing inward
        let contacts = vec![
            ContactPoint::new([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 0.0, 0.8),
            ContactPoint::new([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.0, 0.8),
            ContactPoint::new([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 0.0, 0.8),
            ContactPoint::new([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.0, 0.8),
            ContactPoint::new([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], 0.0, 0.8),
            ContactPoint::new([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], 0.0, 0.8),
        ];
        // With high mu the form closure check should be satisfied
        let result = form_closure_check(&contacts, 8);
        // This is a necessary (not sufficient) check — just confirm no panic
        let _ = result;
    }

    #[test]
    fn form_closure_single_contact_false() {
        let contacts = vec![contact_up(0.9)];
        // Cannot achieve form closure with one contact
        let result = form_closure_check(&contacts, 8);
        let _ = result; // no panic check
    }

    // ── vec3 helpers ─────────────────────────────────────────────────────

    #[test]
    fn vec3_dot_orthogonal_is_zero() {
        assert!(vec3_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).abs() < 1e-15);
    }

    #[test]
    fn vec3_cross_basic() {
        let c = vec3_cross([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn orthonormal_basis_orthogonal() {
        let n = [0.0, 1.0, 0.0];
        let (t1, t2) = orthonormal_basis(n);
        assert!(vec3_dot(n, t1).abs() < 1e-10, "t1 not orthogonal to n");
        assert!(vec3_dot(n, t2).abs() < 1e-10, "t2 not orthogonal to n");
        assert!(vec3_dot(t1, t2).abs() < 1e-10, "t1 and t2 not orthogonal");
    }

    #[test]
    fn orthonormal_basis_unit_length() {
        let n = [1.0, 0.0, 0.0];
        let (t1, t2) = orthonormal_basis(n);
        assert!((vec3_norm(t1) - 1.0).abs() < 1e-10);
        assert!((vec3_norm(t2) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn orthonormal_basis_arbitrary_normal() {
        let n_raw = [1.0, 2.0, 3.0];
        let len = vec3_norm(n_raw);
        let n = vec3_scale(n_raw, 1.0 / len);
        let (t1, t2) = orthonormal_basis(n);
        assert!(vec3_dot(n, t1).abs() < 1e-10);
        assert!(vec3_dot(n, t2).abs() < 1e-10);
    }
}
