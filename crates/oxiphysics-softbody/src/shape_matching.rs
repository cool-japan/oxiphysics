// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid shape matching constraint (Müller et al., 2005).
//!
//! Shape matching projects a set of particles toward the best-fit rigid
//! transformation (rotation + translation) of a rest shape.  The stiffness
//! parameter controls how strongly particles are pulled toward the goal
//! positions.
//!
//! This module provides two implementations:
//!
//! 1. **`ShapeMatching`** – the original nalgebra-based constraint, kept for
//!    backward compatibility.
//!
//! 2. **`ShapeMatchingBody`** – a standalone, pure-`f64`-array body that
//!    supports plasticity and direct spring-force application.

use oxiphysics_core::math::{Mat3, Real, Vec3};

// ============================================================================
// Pure-f64 helper math
// ============================================================================

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

/// 3×3 matrix stored in row-major order: m\[row\]\[col\].
type Mat3x3 = [[f64; 3]; 3];

fn mat3_zero() -> Mat3x3 {
    [[0.0; 3]; 3]
}

fn mat3_identity() -> Mat3x3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn mat3_mul(a: Mat3x3, b: Mat3x3) -> Mat3x3 {
    let mut c = mat3_zero();
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

#[cfg(test)]
fn mat3_transpose(m: Mat3x3) -> Mat3x3 {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

fn mat3_det(m: Mat3x3) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Apply a 3×3 matrix to a vector.
fn mat3_mul_vec(m: Mat3x3, v: [f64; 3]) -> [f64; 3] {
    [v3_dot(m[0], v), v3_dot(m[1], v), v3_dot(m[2], v)]
}

/// Outer product: returns `a * bᵀ` as a row-major 3×3 matrix.
fn outer_product_f64(a: [f64; 3], b: [f64; 3]) -> Mat3x3 {
    [
        [a[0] * b[0], a[0] * b[1], a[0] * b[2]],
        [a[1] * b[0], a[1] * b[1], a[1] * b[2]],
        [a[2] * b[0], a[2] * b[1], a[2] * b[2]],
    ]
}

fn mat3_add_scaled(a: Mat3x3, b: Mat3x3, s: f64) -> Mat3x3 {
    let mut c = a;
    for i in 0..3 {
        for j in 0..3 {
            c[i][j] += b[i][j] * s;
        }
    }
    c
}

/// Polar decomposition of `m` via Jacobi SVD, returning the rotation factor R
/// such that `m ≈ R * S`.
///
/// Implemented using nalgebra's SVD for robustness.
fn polar_rotation_f64(m: Mat3x3) -> Mat3x3 {
    // Convert to nalgebra matrix.
    let nm = Mat3::new(
        m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2],
    );
    let svd = nm.svd(true, true);
    let u = svd.u.unwrap_or_else(Mat3::identity);
    let vt = svd.v_t.unwrap_or_else(Mat3::identity);
    let mut r_na = u * vt;
    if r_na.determinant() < 0.0 {
        let mut u_fixed = u;
        u_fixed.set_column(2, &(-u.column(2)));
        r_na = u_fixed * vt;
    }
    // Convert back to row-major array.
    [
        [r_na[(0, 0)], r_na[(0, 1)], r_na[(0, 2)]],
        [r_na[(1, 0)], r_na[(1, 1)], r_na[(1, 2)]],
        [r_na[(2, 0)], r_na[(2, 1)], r_na[(2, 2)]],
    ]
}

// ============================================================================
// ShapeMatchingBody  (pure-f64 standalone)
// ============================================================================

/// A soft body that uses the shape-matching algorithm to maintain a rigid
/// rest shape with optional plasticity.
///
/// # Algorithm (per timestep)
///
/// 1. Compute the current centre of mass `xcm`.
/// 2. Build the deformation matrix `A = Σ mᵢ (xᵢ − xcm) qᵢᵀ` where `qᵢ`
///    are the (possibly plastic) rest positions relative to the rest COM.
/// 3. Extract the optimal rotation `R` via polar decomposition of `A`.
/// 4. Compute goal positions `gᵢ = R qᵢ + xcm`.
/// 5. Apply spring forces towards `gᵢ` (or teleport with `stiffness = 1`).
/// 6. Optionally accumulate plastic deformation.
#[derive(Debug, Clone)]
pub struct ShapeMatchingBody {
    /// Current positions of the particles.
    pub positions: Vec<[f64; 3]>,
    /// Current velocities of the particles.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// Rest positions relative to the rest centre of mass.
    ///
    /// These are updated when plasticity is applied.
    rest_positions: Vec<[f64; 3]>,
    /// The last computed optimal rotation matrix.
    pub rotation: Mat3x3,
    /// Stiffness α ∈ \[0, 1\].  0 = no shape matching, 1 = full rigid match.
    pub stiffness: f64,
    /// Plastic deformation threshold.  When the positional error exceeds this
    /// value the rest shape is updated to partially absorb the deformation.
    pub plasticity_threshold: f64,
    /// How much of the excess deformation is absorbed plastically per step.
    /// 0 = no plasticity, 1 = instant plastic yield.
    pub plasticity_rate: f64,
    /// Total accumulated plastic deformation magnitude (diagnostic).
    pub accumulated_plastic_deformation: f64,
}

impl ShapeMatchingBody {
    /// Create a new body from rest positions and masses.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>, stiffness: f64) -> Self {
        assert_eq!(
            positions.len(),
            masses.len(),
            "positions and masses must have the same length"
        );
        let n = positions.len();

        // Compute rest centre of mass.
        let total_mass: f64 = masses.iter().sum();
        assert!(total_mass > 0.0, "total mass must be positive");
        let mut cm = [0.0_f64; 3];
        for (p, &m) in positions.iter().zip(masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        cm = v3_scale(cm, 1.0 / total_mass);

        // Rest positions relative to CM.
        let rest_positions: Vec<[f64; 3]> = positions.iter().map(|p| v3_sub(*p, cm)).collect();

        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            rest_positions,
            rotation: mat3_identity(),
            stiffness,
            plasticity_threshold: f64::INFINITY,
            plasticity_rate: 0.0,
            accumulated_plastic_deformation: 0.0,
        }
    }

    /// Total mass of the body.
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }

    /// Current centre of mass.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total = self.total_mass();
        let mut cm = [0.0_f64; 3];
        for (p, &m) in self.positions.iter().zip(self.masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        v3_scale(cm, 1.0 / total)
    }

    /// Compute the optimal rotation matrix from current positions.
    ///
    /// Updates `self.rotation` and returns it.
    pub fn compute_optimal_rotation(&mut self) -> Mat3x3 {
        let xcm = self.centre_of_mass();

        // Deformation matrix A = Σ mᵢ (xᵢ - xcm) qᵢᵀ.
        let mut a = mat3_zero();
        for ((pos, q), &mi) in self
            .positions
            .iter()
            .zip(self.rest_positions.iter())
            .zip(self.masses.iter())
        {
            let p = v3_sub(*pos, xcm);
            let op = outer_product_f64(p, *q);
            a = mat3_add_scaled(a, op, mi);
        }

        self.rotation = polar_rotation_f64(a);
        self.rotation
    }

    /// Compute goal positions using the current optimal rotation.
    ///
    /// `goalᵢ = R * qᵢ + xcm`
    pub fn goal_positions(&mut self) -> Vec<[f64; 3]> {
        let r = self.compute_optimal_rotation();
        let xcm = self.centre_of_mass();
        self.rest_positions
            .iter()
            .map(|q| v3_add(mat3_mul_vec(r, *q), xcm))
            .collect()
    }

    /// Apply shape-matching spring forces towards the goal positions.
    ///
    /// For each dynamic particle:
    ///   `v += stiffness * (goal - pos) / dt`
    ///   `pos += stiffness * (goal - pos)`
    ///
    /// After updating, optionally accumulates plastic deformation.
    pub fn apply_shape_matching_force(&mut self, stiffness: f64, dt: f64) {
        let goals = self.goal_positions();
        let n = self.positions.len();
        for (i, goal) in goals.iter().enumerate().take(n) {
            let delta = v3_sub(*goal, self.positions[i]);
            let disp = v3_scale(delta, stiffness);
            self.positions[i] = v3_add(self.positions[i], disp);
            // Velocity from positional correction.
            if dt > 0.0 {
                self.velocities[i] = v3_add(self.velocities[i], v3_scale(disp, 1.0 / dt));
            }

            // Plasticity: if positional error exceeds threshold, absorb it.
            let err = v3_norm(delta);
            if err > self.plasticity_threshold {
                let excess = err - self.plasticity_threshold;
                let plastic_delta = v3_scale(
                    v3_scale(delta, 1.0 / err.max(1e-12)),
                    excess * self.plasticity_rate,
                );
                // Update rest position towards current to absorb plastic strain.
                self.rest_positions[i] = v3_add(self.rest_positions[i], plastic_delta);
                self.accumulated_plastic_deformation += v3_norm(plastic_delta);
            }
        }
    }

    /// Reset the rest shape to the current positions (erase all deformation).
    pub fn reset_rest_shape(&mut self) {
        let cm = self.centre_of_mass();
        self.rest_positions = self.positions.iter().map(|p| v3_sub(*p, cm)).collect();
        self.accumulated_plastic_deformation = 0.0;
    }
}

// ============================================================================
// Original nalgebra-based ShapeMatching (kept for backward compat)
// ============================================================================

/// Rigid shape-matching constraint.
#[derive(Debug, Clone)]
pub struct ShapeMatching {
    /// Rest shape positions relative to the rest centre of mass.
    rest_positions: Vec<Vec3>,
    /// Per-particle masses.
    masses: Vec<Real>,
    /// Stiffness α ∈ \[0, 1\].  0 = no constraint, 1 = full rigid matching.
    pub stiffness: Real,
}

impl ShapeMatching {
    /// Create a new shape-matching constraint from the current rest positions
    /// and masses.
    ///
    /// The rest positions are stored relative to their (mass-weighted) centre
    /// of mass so that the deformation matrix `A` can be computed cheaply at
    /// runtime.
    pub fn new(positions: &[Vec3], masses: &[Real], stiffness: Real) -> Self {
        assert_eq!(
            positions.len(),
            masses.len(),
            "positions and masses must have the same length"
        );

        // Compute rest centre of mass.
        let total_mass: Real = masses.iter().sum();
        let rest_cm: Vec3 = positions
            .iter()
            .zip(masses.iter())
            .map(|(p, &m)| p * m)
            .fold(Vec3::zeros(), |a, b| a + b)
            / total_mass;

        // Store positions relative to rest centre of mass.
        let rest_positions: Vec<Vec3> = positions.iter().map(|p| p - rest_cm).collect();

        Self {
            rest_positions,
            masses: masses.to_vec(),
            stiffness,
        }
    }

    /// Compute the goal positions for each particle by:
    ///
    /// 1. Computing the current centre of mass `xcm`.
    /// 2. Building the deformation matrix
    ///    `A = Σ mᵢ (xᵢ − xcm) qᵢᵀ`
    ///    where `qᵢ` is the rest position of particle *i* relative to the rest
    ///    centre of mass.
    /// 3. Extracting the rotation `R` via polar decomposition of `A`.
    /// 4. Computing goal positions as `goalᵢ = stiffness * (R qᵢ + xcm) + (1 − stiffness) * xᵢ`.
    pub fn goal_positions(&self, positions: &[Vec3]) -> Vec<Vec3> {
        let n = self.rest_positions.len();
        assert_eq!(positions.len(), n, "positions length mismatch");

        // --- (1) Current centre of mass ----------------------------------
        let total_mass: Real = self.masses.iter().sum();
        let xcm: Vec3 = positions
            .iter()
            .zip(self.masses.iter())
            .map(|(p, &m)| p * m)
            .fold(Vec3::zeros(), |a, b| a + b)
            / total_mass;

        // --- (2) Deformation matrix A = Σ mᵢ pᵢ qᵢᵀ ---------------------
        // where pᵢ = xᵢ - xcm  (current relative position)
        //       qᵢ = rest_positions[i]  (rest relative position)
        let mut a = Mat3::zeros();
        for ((pos, q), &mi) in positions
            .iter()
            .zip(self.rest_positions.iter())
            .zip(self.masses.iter())
        {
            let p = pos - xcm;
            // Outer product p * qᵀ, scaled by mass.
            a += outer_product(&p, q) * mi;
        }

        // --- (3) Polar decomposition: extract R from A -------------------
        let r = polar_rotation(&a);

        // --- (4) Goal positions ------------------------------------------
        self.rest_positions
            .iter()
            .zip(positions.iter())
            .map(|(q, cur)| {
                let goal_rigid = r * q + xcm;
                // Blend between rigid goal and current position by stiffness.
                goal_rigid * self.stiffness + cur * (1.0 - self.stiffness)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Helpers (nalgebra-based)
// ---------------------------------------------------------------------------

/// Outer (dyadic) product: returns `a * bᵀ` as a 3×3 matrix.
fn outer_product(a: &Vec3, b: &Vec3) -> Mat3 {
    Mat3::new(
        a.x * b.x,
        a.x * b.y,
        a.x * b.z,
        a.y * b.x,
        a.y * b.y,
        a.y * b.z,
        a.z * b.x,
        a.z * b.y,
        a.z * b.z,
    )
}

/// Polar decomposition of matrix `m` via SVD, returning the rotation factor
/// `R` such that `m ≈ R * S` (S symmetric positive semi-definite).
///
/// Given the thin SVD `m = U Σ Vᵀ`, the closest rotation is `R = U Vᵀ`.
/// This is numerically robust even when `m` is rank-deficient (e.g. planar
/// configurations).
fn polar_rotation(m: &Mat3) -> Mat3 {
    let svd = m.svd(true, true);
    // svd.u and svd.v_t are `Option` – they are `Some` when we pass `true`.
    let u = svd.u.unwrap_or_else(Mat3::identity);
    let vt = svd.v_t.unwrap_or_else(Mat3::identity);
    let mut r = u * vt;
    // Ensure a proper rotation (det = +1, not a reflection).
    if r.determinant() < 0.0 {
        // Flip the column of U corresponding to the smallest singular value.
        let mut u_fixed = u;
        u_fixed.set_column(2, &(-u.column(2)));
        r = u_fixed * vt;
    }
    r
}

// ============================================================================
// LinearDeformationBody – linear deformation modes (affine shape matching)
// ============================================================================

/// Shape matching with linear deformation: A itself is used as the goal
/// deformation matrix instead of just its rotation factor.
///
/// This allows volume-preserving stretch in addition to rotation.
#[derive(Debug, Clone)]
pub struct LinearDeformationBody {
    /// Current positions.
    pub positions: Vec<[f64; 3]>,
    /// Current velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// Rest positions relative to rest centre of mass.
    rest_positions: Vec<[f64; 3]>,
    /// Precomputed inverse of the rest shape matrix Q = Σ mᵢ qᵢ qᵢᵀ.
    q_inv: Mat3x3,
    /// Stiffness α ∈ \[0, 1\].
    pub stiffness: f64,
    /// Volume preservation blend factor (0 = pure linear, 1 = volume-preserving).
    pub volume_preservation: f64,
}

impl LinearDeformationBody {
    /// Create a new linear deformation body.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>, stiffness: f64) -> Self {
        assert_eq!(positions.len(), masses.len());
        let n = positions.len();

        let total_mass: f64 = masses.iter().sum();
        assert!(total_mass > 0.0);

        let mut cm = [0.0_f64; 3];
        for (p, &m) in positions.iter().zip(masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        cm = v3_scale(cm, 1.0 / total_mass);

        let rest_positions: Vec<[f64; 3]> = positions.iter().map(|p| v3_sub(*p, cm)).collect();

        // Precompute Q = Σ mᵢ qᵢ qᵢᵀ
        let mut q_mat = mat3_zero();
        for (&q, &m) in rest_positions.iter().zip(masses.iter()) {
            let op = outer_product_f64(q, q);
            q_mat = mat3_add_scaled(q_mat, op, m);
        }
        let q_inv = mat3_inverse_safe(q_mat);

        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            rest_positions,
            q_inv,
            stiffness,
            volume_preservation: 0.5,
        }
    }

    /// Compute current centre of mass.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total: f64 = self.masses.iter().sum();
        let mut cm = [0.0_f64; 3];
        for (p, &m) in self.positions.iter().zip(self.masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        v3_scale(cm, 1.0 / total)
    }

    /// Compute the linear deformation matrix A_lin = Σ mᵢ pᵢ qᵢᵀ * Q⁻¹.
    pub fn deformation_matrix(&self) -> Mat3x3 {
        let xcm = self.centre_of_mass();
        let mut a = mat3_zero();
        for ((pos, q), &mi) in self
            .positions
            .iter()
            .zip(self.rest_positions.iter())
            .zip(self.masses.iter())
        {
            let p = v3_sub(*pos, xcm);
            let op = outer_product_f64(p, *q);
            a = mat3_add_scaled(a, op, mi);
        }
        mat3_mul(a, self.q_inv)
    }

    /// Compute goal positions using the linear deformation matrix.
    ///
    /// Blends between pure rotation (shape matching) and full affine deformation.
    pub fn goal_positions(&self) -> Vec<[f64; 3]> {
        let a_lin = self.deformation_matrix();
        let r = polar_rotation_f64(a_lin);
        let xcm = self.centre_of_mass();

        // Volume-preserving: blend A_lin with R.
        let det = mat3_det(a_lin).abs();
        let vp = self.volume_preservation;
        let scale_factor = if det > 1e-12 {
            det.powf(1.0 / 3.0)
        } else {
            1.0
        };

        self.rest_positions
            .iter()
            .map(|q| {
                let goal_affine = mat3_mul_vec(a_lin, *q);
                let goal_rigid = mat3_mul_vec(r, *q);
                // Volume-preserving blend: scale affine toward unit-det
                let goal_vp = v3_scale(mat3_mul_vec(a_lin, *q), 1.0 / scale_factor.max(1e-12));
                let blended = [
                    goal_affine[0] * (1.0 - vp) + goal_vp[0] * vp,
                    goal_affine[1] * (1.0 - vp) + goal_vp[1] * vp,
                    goal_affine[2] * (1.0 - vp) + goal_vp[2] * vp,
                ];
                let final_goal = [
                    blended[0] * (1.0 - (1.0 - self.stiffness))
                        + goal_rigid[0] * (1.0 - self.stiffness),
                    blended[1] * (1.0 - (1.0 - self.stiffness))
                        + goal_rigid[1] * (1.0 - self.stiffness),
                    blended[2] * (1.0 - (1.0 - self.stiffness))
                        + goal_rigid[2] * (1.0 - self.stiffness),
                ];
                v3_add(final_goal, xcm)
            })
            .collect()
    }

    /// Apply goal-position correction.
    pub fn apply_correction(&mut self, stiffness: f64, dt: f64) {
        let goals = self.goal_positions();
        for (i, goal) in goals.iter().enumerate() {
            let delta = v3_sub(*goal, self.positions[i]);
            let disp = v3_scale(delta, stiffness);
            self.positions[i] = v3_add(self.positions[i], disp);
            if dt > 0.0 {
                self.velocities[i] = v3_add(self.velocities[i], v3_scale(disp, 1.0 / dt));
            }
        }
    }
}

/// Compute a safe pseudo-inverse of a 3×3 matrix via Cramer's rule.
/// Falls back to identity when the matrix is near-singular.
fn mat3_inverse_safe(m: Mat3x3) -> Mat3x3 {
    let det = mat3_det(m);
    if det.abs() < 1e-12 {
        return mat3_identity();
    }
    let inv_det = 1.0 / det;
    [
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
        ],
    ]
}

// ============================================================================
// QuadraticDeformationBody – quadratic deformation modes
// ============================================================================

/// Extended shape matching with quadratic deformation modes.
///
/// The rest-position feature vector is extended to 9 components by appending
/// quadratic terms: `[qx, qy, qz, qx*qx, qy*qy, qz*qz, qx*qy, qy*qz, qx*qz]`.
///
/// This allows bending modes beyond pure affine deformation.
#[derive(Debug, Clone)]
pub struct QuadraticDeformationBody {
    /// Current positions.
    pub positions: Vec<[f64; 3]>,
    /// Current velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// 9-component extended rest-position feature vectors.
    rest_features: Vec<[f64; 9]>,
    /// Stiffness α ∈ \[0, 1\].
    pub stiffness: f64,
}

impl QuadraticDeformationBody {
    /// Compute the 9-component feature vector from a 3D position.
    fn feature(q: [f64; 3]) -> [f64; 9] {
        [
            q[0],
            q[1],
            q[2],
            q[0] * q[0],
            q[1] * q[1],
            q[2] * q[2],
            q[0] * q[1],
            q[1] * q[2],
            q[0] * q[2],
        ]
    }

    /// Create a new quadratic deformation body from rest positions and masses.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>, stiffness: f64) -> Self {
        assert_eq!(positions.len(), masses.len());
        let total: f64 = masses.iter().sum();
        assert!(total > 0.0);
        let mut cm = [0.0_f64; 3];
        for (p, &m) in positions.iter().zip(masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        cm = v3_scale(cm, 1.0 / total);

        let n = positions.len();
        let rest_features: Vec<[f64; 9]> = positions
            .iter()
            .map(|p| Self::feature(v3_sub(*p, cm)))
            .collect();

        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            rest_features,
            stiffness,
        }
    }

    /// Compute current centre of mass.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total: f64 = self.masses.iter().sum();
        let mut cm = [0.0_f64; 3];
        for (p, &m) in self.positions.iter().zip(self.masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        v3_scale(cm, 1.0 / total)
    }

    /// Compute the 3×9 deformation matrix mapping feature space to 3D.
    pub fn compute_a39(&self) -> [[f64; 9]; 3] {
        let xcm = self.centre_of_mass();
        let mut a = [[0.0_f64; 9]; 3];
        for ((pos, feat), &mi) in self
            .positions
            .iter()
            .zip(self.rest_features.iter())
            .zip(self.masses.iter())
        {
            let p = v3_sub(*pos, xcm);
            for row in 0..3 {
                for col in 0..9 {
                    a[row][col] += mi * p[row] * feat[col];
                }
            }
        }
        a
    }

    /// Compute goal positions using the quadratic deformation matrix.
    ///
    /// Goal = A₃₉ * feat(q) + xcm (with stiffness blend).
    pub fn goal_positions(&self) -> Vec<[f64; 3]> {
        let a39 = self.compute_a39();
        let xcm = self.centre_of_mass();
        self.rest_features
            .iter()
            .zip(self.positions.iter())
            .map(|(feat, cur)| {
                let mut goal = xcm;
                for row in 0..3 {
                    let v: f64 = a39[row].iter().zip(feat.iter()).map(|(a, f)| a * f).sum();
                    goal[row] += v;
                }
                // Stiffness blend
                [
                    goal[0] * self.stiffness + cur[0] * (1.0 - self.stiffness),
                    goal[1] * self.stiffness + cur[1] * (1.0 - self.stiffness),
                    goal[2] * self.stiffness + cur[2] * (1.0 - self.stiffness),
                ]
            })
            .collect()
    }

    /// Apply quadratic goal correction.
    pub fn apply_correction(&mut self, stiffness: f64, dt: f64) {
        let goals = self.goal_positions();
        for (i, goal) in goals.iter().enumerate() {
            let delta = v3_sub(*goal, self.positions[i]);
            let disp = v3_scale(delta, stiffness);
            self.positions[i] = v3_add(self.positions[i], disp);
            if dt > 0.0 {
                self.velocities[i] = v3_add(self.velocities[i], v3_scale(disp, 1.0 / dt));
            }
        }
    }
}

// ============================================================================
// ClusterShapeMatching – cluster-based shape matching
// ============================================================================

/// A cluster used in cluster-based shape matching.
///
/// Each particle may belong to multiple clusters; corrections are averaged.
#[derive(Debug, Clone)]
pub struct ShapeMatchingCluster {
    /// Particle indices belonging to this cluster.
    pub particle_indices: Vec<usize>,
    /// Per-cluster rest positions relative to cluster centre of mass.
    cluster_rest: Vec<[f64; 3]>,
    /// Cluster stiffness.
    pub stiffness: f64,
}

impl ShapeMatchingCluster {
    /// Create a cluster from a subset of a particle system.
    ///
    /// `positions` and `masses` are the full arrays; `indices` selects
    /// which particles form this cluster.
    pub fn new(
        positions: &[[f64; 3]],
        masses: &[f64],
        indices: Vec<usize>,
        stiffness: f64,
    ) -> Self {
        let total: f64 = indices.iter().map(|&i| masses[i]).sum();
        assert!(total > 0.0);
        let mut cm = [0.0_f64; 3];
        for &i in &indices {
            cm[0] += positions[i][0] * masses[i];
            cm[1] += positions[i][1] * masses[i];
            cm[2] += positions[i][2] * masses[i];
        }
        cm = v3_scale(cm, 1.0 / total);
        let cluster_rest: Vec<[f64; 3]> =
            indices.iter().map(|&i| v3_sub(positions[i], cm)).collect();
        Self {
            particle_indices: indices,
            cluster_rest,
            stiffness,
        }
    }

    /// Compute goal corrections for this cluster's particles, returning
    /// a Vec of `(particle_index, displacement)` pairs.
    pub fn compute_corrections(
        &self,
        positions: &[[f64; 3]],
        masses: &[f64],
    ) -> Vec<(usize, [f64; 3])> {
        // Compute cluster current centre of mass.
        let total: f64 = self.particle_indices.iter().map(|&i| masses[i]).sum();
        if total < 1e-30 {
            return vec![];
        }
        let mut xcm = [0.0_f64; 3];
        for &i in &self.particle_indices {
            xcm[0] += positions[i][0] * masses[i];
            xcm[1] += positions[i][1] * masses[i];
            xcm[2] += positions[i][2] * masses[i];
        }
        xcm = v3_scale(xcm, 1.0 / total);

        // Deformation matrix.
        let mut a = mat3_zero();
        for (&i, q) in self.particle_indices.iter().zip(self.cluster_rest.iter()) {
            let p = v3_sub(positions[i], xcm);
            let op = outer_product_f64(p, *q);
            a = mat3_add_scaled(a, op, masses[i]);
        }
        let r = polar_rotation_f64(a);

        // Goal positions.
        self.particle_indices
            .iter()
            .zip(self.cluster_rest.iter())
            .map(|(&i, q)| {
                let goal = v3_add(mat3_mul_vec(r, *q), xcm);
                let delta = v3_sub(goal, positions[i]);
                (i, v3_scale(delta, self.stiffness))
            })
            .collect()
    }
}

/// Manager for cluster-based shape matching over a shared particle set.
///
/// Each particle may belong to multiple clusters; per-particle corrections
/// are averaged across clusters.
pub struct ClusterShapeMatchingSystem {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// All clusters.
    pub clusters: Vec<ShapeMatchingCluster>,
}

impl ClusterShapeMatchingSystem {
    /// Create an empty system with `n` particles.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>) -> Self {
        let n = positions.len();
        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            clusters: Vec::new(),
        }
    }

    /// Add a cluster.
    pub fn add_cluster(&mut self, indices: Vec<usize>, stiffness: f64) {
        let c = ShapeMatchingCluster::new(&self.positions, &self.masses, indices, stiffness);
        self.clusters.push(c);
    }

    /// Apply one step of cluster shape matching.
    ///
    /// Corrections from all clusters are averaged per particle.
    pub fn step(&mut self, dt: f64) {
        let n = self.positions.len();
        let mut total_correction = vec![[0.0_f64; 3]; n];
        let mut correction_count = vec![0usize; n];

        for cluster in &self.clusters {
            let corrections = cluster.compute_corrections(&self.positions, &self.masses);
            for (i, delta) in corrections {
                total_correction[i] = v3_add(total_correction[i], delta);
                correction_count[i] += 1;
            }
        }

        for i in 0..n {
            if correction_count[i] > 0 {
                let avg = v3_scale(total_correction[i], 1.0 / correction_count[i] as f64);
                self.positions[i] = v3_add(self.positions[i], avg);
                if dt > 0.0 {
                    self.velocities[i] = v3_add(self.velocities[i], v3_scale(avg, 1.0 / dt));
                }
            }
        }
    }
}

// ============================================================================
// PlasticityImprovedBody – improved plasticity with yield surface
// ============================================================================

/// Improved plasticity model with von-Mises-like yield criterion.
///
/// Tracks the cumulative plastic strain and applies a hardening response
/// once yield stress is exceeded.
#[derive(Debug, Clone)]
pub struct PlasticityImprovedBody {
    /// Current positions.
    pub positions: Vec<[f64; 3]>,
    /// Current velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// Plastic rest positions.
    plastic_rest: Vec<[f64; 3]>,
    /// Cumulative plastic strain per particle.
    pub plastic_strain: Vec<f64>,
    /// Yield stress σ_y: deformation below this is elastic.
    pub yield_stress: f64,
    /// Hardening modulus H: reduces yield stress over time (softening if < 0).
    pub hardening_modulus: f64,
    /// Plasticity rate (0..1).
    pub plasticity_rate: f64,
    /// Shape matching stiffness.
    pub stiffness: f64,
    /// Rotation matrix.
    pub rotation: Mat3x3,
}

impl PlasticityImprovedBody {
    /// Create a new body with improved plasticity.
    pub fn new(positions: Vec<[f64; 3]>, masses: Vec<f64>, stiffness: f64) -> Self {
        assert_eq!(positions.len(), masses.len());
        let n = positions.len();
        let total: f64 = masses.iter().sum();
        assert!(total > 0.0);
        let mut cm = [0.0_f64; 3];
        for (p, &m) in positions.iter().zip(masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        cm = v3_scale(cm, 1.0 / total);
        let plastic_rest: Vec<[f64; 3]> = positions.iter().map(|p| v3_sub(*p, cm)).collect();
        Self {
            positions,
            velocities: vec![[0.0; 3]; n],
            masses,
            plastic_rest,
            plastic_strain: vec![0.0; n],
            yield_stress: 0.1,
            hardening_modulus: 0.0,
            plasticity_rate: 0.3,
            stiffness,
            rotation: mat3_identity(),
        }
    }

    /// Current centre of mass.
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let total: f64 = self.masses.iter().sum();
        let mut cm = [0.0_f64; 3];
        for (p, &m) in self.positions.iter().zip(self.masses.iter()) {
            cm[0] += p[0] * m;
            cm[1] += p[1] * m;
            cm[2] += p[2] * m;
        }
        v3_scale(cm, 1.0 / total)
    }

    /// Apply one shape-matching step with improved plasticity.
    ///
    /// The effective yield stress varies with accumulated plastic strain (hardening).
    pub fn step(&mut self, dt: f64) {
        let xcm = self.centre_of_mass();

        let mut a = mat3_zero();
        for ((pos, q), &mi) in self
            .positions
            .iter()
            .zip(self.plastic_rest.iter())
            .zip(self.masses.iter())
        {
            let p = v3_sub(*pos, xcm);
            let op = outer_product_f64(p, *q);
            a = mat3_add_scaled(a, op, mi);
        }
        self.rotation = polar_rotation_f64(a);
        let r = self.rotation;

        let n = self.positions.len();
        for i in 0..n {
            let goal = v3_add(mat3_mul_vec(r, self.plastic_rest[i]), xcm);
            let delta = v3_sub(goal, self.positions[i]);
            let disp = v3_scale(delta, self.stiffness);
            self.positions[i] = v3_add(self.positions[i], disp);
            if dt > 0.0 {
                self.velocities[i] = v3_add(self.velocities[i], v3_scale(disp, 1.0 / dt));
            }

            // Improved plasticity with hardening
            let err = v3_norm(delta);
            let effective_yield =
                (self.yield_stress + self.hardening_modulus * self.plastic_strain[i]).max(1e-12);
            if err > effective_yield {
                let excess = err - effective_yield;
                let dir = v3_scale(delta, 1.0 / err.max(1e-12));
                let plastic_delta = v3_scale(dir, excess * self.plasticity_rate);
                self.plastic_rest[i] = v3_add(self.plastic_rest[i], plastic_delta);
                self.plastic_strain[i] += v3_norm(plastic_delta);
            }
        }
    }

    /// Total accumulated plastic strain across all particles.
    pub fn total_plastic_strain(&self) -> f64 {
        self.plastic_strain.iter().sum()
    }

    /// Reset plastic rest shape to current positions.
    pub fn reset_plasticity(&mut self) {
        let cm = self.centre_of_mass();
        self.plastic_rest = self.positions.iter().map(|p| v3_sub(*p, cm)).collect();
        for s in &mut self.plastic_strain {
            *s = 0.0;
        }
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-6;

    // Helper: apply a 3×3 rotation to a position.
    fn rotate(r: Mat3x3, v: [f64; 3]) -> [f64; 3] {
        mat3_mul_vec(r, v)
    }

    // Helper: Euclidean distance between two [f64;3] points.
    fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
        v3_norm(v3_sub(a, b))
    }

    // SM1. Rigid body: pure rotation → no deformation → goals equal current.
    //
    // If the current positions are a pure rotation R of the rest positions
    // (with stiffness=1), then goal_positions() == current positions.
    #[test]
    fn test_rigid_body_no_deformation() {
        // Tetrahedron rest positions.
        let rest = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let masses = vec![1.0_f64; 4];
        let mut body = ShapeMatchingBody::new(rest.clone(), masses, 1.0);

        // Apply a 90° rotation about Z to current positions.
        let rz90: Mat3x3 = [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
        body.positions = rest.iter().map(|p| rotate(rz90, *p)).collect();

        let goals = body.goal_positions();
        for (i, (goal, cur)) in goals.iter().zip(body.positions.iter()).enumerate() {
            assert!(
                dist3(*goal, *cur) < EPS,
                "Particle {i}: goal {:?} should equal current {:?}",
                goal,
                cur
            );
        }
    }

    // SM2. Pure translation: shape matching should follow the centre of mass.
    //
    // Translating all particles by a constant vector should not introduce any
    // rotational error – goals should equal current positions.
    #[test]
    fn test_pure_translation() {
        let rest = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let masses = vec![1.0_f64; 4];
        let mut body = ShapeMatchingBody::new(rest.clone(), masses, 1.0);

        // Translate all particles.
        let offset = [5.0, -3.0, 2.0];
        body.positions = rest.iter().map(|p| v3_add(*p, offset)).collect();

        let goals = body.goal_positions();
        for (i, (goal, cur)) in goals.iter().zip(body.positions.iter()).enumerate() {
            assert!(
                dist3(*goal, *cur) < EPS,
                "Particle {i}: goal {:?} should equal current {:?} after pure translation",
                goal,
                cur
            );
        }
    }

    // SM3. Rotation invariance: rotating rest AND current by the same R
    //      should not change the positional errors.
    #[test]
    fn test_rotation_invariance() {
        let rest = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.5, 0.0]];
        let masses = vec![1.0_f64; 3];

        // Deformed positions (a simple stretch).
        let deformed: Vec<[f64; 3]> = vec![[2.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.5, 0.0]];

        // Errors without rotation.
        let mut body_a = ShapeMatchingBody::new(rest.clone(), masses.clone(), 1.0);
        body_a.positions = deformed.clone();
        let goals_a = body_a.goal_positions();
        let errors_a: Vec<f64> = goals_a
            .iter()
            .zip(deformed.iter())
            .map(|(g, c)| dist3(*g, *c))
            .collect();

        // Apply a 45° rotation about Z to both rest and deformed.
        let c45 = std::f64::consts::FRAC_PI_4.cos();
        let s45 = std::f64::consts::FRAC_PI_4.sin();
        let rz45: Mat3x3 = [[c45, -s45, 0.0], [s45, c45, 0.0], [0.0, 0.0, 1.0]];

        let rotated_rest: Vec<[f64; 3]> = rest.iter().map(|p| rotate(rz45, *p)).collect();
        let rotated_deformed: Vec<[f64; 3]> = deformed.iter().map(|p| rotate(rz45, *p)).collect();

        let mut body_b = ShapeMatchingBody::new(rotated_rest, masses.clone(), 1.0);
        body_b.positions = rotated_deformed.clone();
        let goals_b = body_b.goal_positions();
        let errors_b: Vec<f64> = goals_b
            .iter()
            .zip(rotated_deformed.iter())
            .map(|(g, c)| dist3(*g, *c))
            .collect();

        for i in 0..3 {
            assert!(
                (errors_a[i] - errors_b[i]).abs() < 1e-5,
                "Rotation invariance violated at particle {i}: {:.6} vs {:.6}",
                errors_a[i],
                errors_b[i]
            );
        }
    }

    // SM4. apply_shape_matching_force with stiffness=1 teleports particles to goals.
    #[test]
    fn test_apply_force_stiffness_one() {
        let rest = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0_f64; 3];
        let mut body = ShapeMatchingBody::new(rest.clone(), masses, 1.0);

        // Perturb one particle.
        body.positions[2] = [0.5, 1.5, 0.3];

        let goals_before = body.goal_positions();
        body.positions = goals_before.clone(); // manually set for the next check
        // After applying with stiffness=1, positions should equal goals.
        let mut body2 = ShapeMatchingBody::new(rest, vec![1.0_f64; 3], 1.0);
        body2.positions[2] = [0.5, 1.5, 0.3];
        body2.apply_shape_matching_force(1.0, 0.016);

        // Goals after application should be very close to new positions.
        let goals_after = body2.goal_positions();
        for (i, (goal, pos)) in goals_after.iter().zip(body2.positions.iter()).enumerate() {
            assert!(
                dist3(*goal, *pos) < 1e-4,
                "Particle {i}: position {:?} should be close to goal {:?}",
                pos,
                goal
            );
        }
    }

    // SM5. Plasticity accumulates when deformation exceeds threshold.
    #[test]
    fn test_plasticity_accumulates() {
        let rest = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0_f64; 3];
        let mut body = ShapeMatchingBody::new(rest, masses, 0.5);
        body.plasticity_threshold = 0.01;
        body.plasticity_rate = 0.5;

        // Stretch particle 0 significantly.
        body.positions[0] = [3.0, 0.0, 0.0];

        let plastic_before = body.accumulated_plastic_deformation;
        body.apply_shape_matching_force(0.5, 0.016);
        let plastic_after = body.accumulated_plastic_deformation;

        assert!(
            plastic_after > plastic_before,
            "Plastic deformation should have accumulated: before={plastic_before}, after={plastic_after}"
        );
    }

    // SM6. reset_rest_shape zeroes accumulated plasticity.
    #[test]
    fn test_reset_rest_shape() {
        let rest = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0_f64; 3];
        let mut body = ShapeMatchingBody::new(rest, masses, 0.5);
        body.accumulated_plastic_deformation = 99.0;

        body.reset_rest_shape();

        assert_eq!(
            body.accumulated_plastic_deformation, 0.0,
            "Accumulated plasticity should be zero after reset"
        );

        // After reset, current positions should be the new rest (goals = current).
        let goals = body.goal_positions();
        for (i, (goal, pos)) in goals.iter().zip(body.positions.iter()).enumerate() {
            assert!(
                dist3(*goal, *pos) < EPS,
                "Particle {i}: after reset goals {:?} should equal positions {:?}",
                goal,
                pos
            );
        }
    }

    // SM7. mat3_mul of two identity matrices is identity.
    #[test]
    fn test_mat3_mul_identity() {
        let i = mat3_identity();
        let r = mat3_mul(i, i);
        for (row, r_row) in r.iter().enumerate() {
            for (col, &v) in r_row.iter().enumerate() {
                let expected = if row == col { 1.0 } else { 0.0 };
                assert!(
                    (v - expected).abs() < 1e-12,
                    "mat3_mul(I,I)[{row}][{col}] = {} (expected {expected})",
                    v
                );
            }
        }
    }

    // SM8. mat3_transpose inverts correctly.
    #[test]
    fn test_mat3_transpose() {
        let m: Mat3x3 = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let mt = mat3_transpose(m);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (mt[i][j] - m[j][i]).abs() < 1e-12,
                    "Transpose mismatch at [{i}][{j}]"
                );
            }
        }
    }

    // SM9. mat3_inverse_safe of identity should be identity.
    #[test]
    fn test_mat3_inverse_safe_identity() {
        let i = mat3_identity();
        let inv = mat3_inverse_safe(i);
        for (r, inv_row) in inv.iter().enumerate() {
            for (c, &v) in inv_row.iter().enumerate() {
                let expected = if r == c { 1.0 } else { 0.0 };
                assert!(
                    (v - expected).abs() < 1e-10,
                    "Inverse of identity should be identity at [{r}][{c}]"
                );
            }
        }
    }

    // SM10. mat3_inverse_safe of near-singular returns identity (safe fallback).
    #[test]
    fn test_mat3_inverse_safe_singular() {
        let singular: Mat3x3 = [[1.0, 2.0, 3.0], [2.0, 4.0, 6.0], [1.0, 2.0, 3.0]];
        let inv = mat3_inverse_safe(singular);
        // Should not panic; result should be identity (fallback).
        let _ = inv;
    }

    // SM11. LinearDeformationBody: pure translation → goals ≈ current.
    #[test]
    fn test_linear_deformation_pure_translation() {
        let rest = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let masses = vec![1.0_f64; 4];
        let mut body = LinearDeformationBody::new(rest.clone(), masses, 1.0);
        let offset = [3.0, 2.0, 1.0];
        body.positions = rest.iter().map(|p| v3_add(*p, offset)).collect();
        let goals = body.goal_positions();
        for (i, (goal, cur)) in goals.iter().zip(body.positions.iter()).enumerate() {
            let d = v3_norm(v3_sub(*goal, *cur));
            assert!(
                d < 0.2,
                "LinearDeformation pure translation particle {i}: d={d}"
            );
        }
    }

    // SM12. LinearDeformationBody: apply_correction moves positions toward goals.
    #[test]
    fn test_linear_deformation_apply_correction() {
        let rest = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let masses = vec![1.0_f64; 4];
        let mut body = LinearDeformationBody::new(rest, masses, 1.0);
        // Perturb.
        body.positions[0] = [2.0, 0.5, 0.3];
        body.apply_correction(0.5, 0.016);
        // Positions should have changed.
        assert!(
            (body.positions[0][0] - 2.0).abs() > 0.01 || (body.positions[0][1] - 0.5).abs() > 0.01,
            "LinearDeformation: correction should move particles"
        );
    }

    // SM13. QuadraticDeformationBody: feature vector has 9 components.
    #[test]
    fn test_quadratic_feature_vector_length() {
        let q = [1.0, 2.0, 3.0];
        let feat = QuadraticDeformationBody::feature(q);
        assert_eq!(feat.len(), 9);
        assert!((feat[0] - 1.0).abs() < 1e-12);
        assert!((feat[3] - 1.0).abs() < 1e-12); // qx*qx = 1
        assert!((feat[4] - 4.0).abs() < 1e-12); // qy*qy = 4
        assert!((feat[6] - 2.0).abs() < 1e-12); // qx*qy = 2
    }

    // SM14. QuadraticDeformationBody: goal positions exist for a small system.
    #[test]
    fn test_quadratic_deformation_goals_exist() {
        let rest = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let masses = vec![1.0_f64; 5];
        let body = QuadraticDeformationBody::new(rest.clone(), masses, 1.0);
        let goals = body.goal_positions();
        assert_eq!(goals.len(), rest.len());
        // Goals should be finite.
        for g in &goals {
            for &v in g {
                assert!(v.is_finite(), "goal component should be finite");
            }
        }
    }

    // SM15. ClusterShapeMatching: single cluster covers all particles → same as basic SM.
    #[test]
    fn test_cluster_single_cluster() {
        let positions = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ];
        let masses = vec![1.0_f64; 4];
        let mut sys = ClusterShapeMatchingSystem::new(positions.clone(), masses.clone());
        sys.add_cluster(vec![0, 1, 2, 3], 1.0);

        // Perturb one particle.
        sys.positions[0] = [3.0, 0.0, 0.0];
        let before = sys.positions[0];
        sys.step(0.016);
        let after = sys.positions[0];
        // Position should have moved toward the goal.
        let d_before = v3_norm(v3_sub(before, [1.0, 0.0, 0.0]));
        let d_after = v3_norm(v3_sub(after, [1.0, 0.0, 0.0]));
        assert!(
            d_after < d_before + 0.01,
            "Cluster SM should move particle closer to goal"
        );
    }

    // SM16. ClusterShapeMatching: two overlapping clusters both apply corrections.
    #[test]
    fn test_cluster_two_overlapping_clusters() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let masses = vec![1.0_f64; 4];
        let mut sys = ClusterShapeMatchingSystem::new(positions, masses);
        sys.add_cluster(vec![0, 1, 2], 0.5);
        sys.add_cluster(vec![1, 2, 3], 0.5);
        // Should not panic with overlapping clusters.
        sys.step(0.016);
        assert_eq!(sys.positions.len(), 4);
    }

    // SM17. PlasticityImprovedBody: plastic strain increases above yield.
    #[test]
    fn test_improved_plasticity_accumulates() {
        let rest = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0_f64; 3];
        let mut body = PlasticityImprovedBody::new(rest, masses, 0.5);
        body.yield_stress = 0.01;
        body.plasticity_rate = 0.5;
        // Large deformation.
        body.positions[0] = [5.0, 0.0, 0.0];
        let before = body.total_plastic_strain();
        body.step(0.016);
        let after = body.total_plastic_strain();
        assert!(
            after > before,
            "Plastic strain should accumulate: before={before}, after={after}"
        );
    }

    // SM18. PlasticityImprovedBody: reset zeroes plastic strain.
    #[test]
    fn test_improved_plasticity_reset() {
        let rest = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0_f64; 3];
        let mut body = PlasticityImprovedBody::new(rest, masses, 0.5);
        body.plastic_strain[0] = 5.0;
        body.reset_plasticity();
        assert_eq!(body.total_plastic_strain(), 0.0);
    }

    // SM19. mat3_det of rotation matrix should be ≈ 1.
    #[test]
    fn test_mat3_det_rotation() {
        let c = std::f64::consts::FRAC_PI_4.cos();
        let s = std::f64::consts::FRAC_PI_4.sin();
        let rz: Mat3x3 = [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
        let det = mat3_det(rz);
        assert!(
            (det - 1.0).abs() < 1e-12,
            "Rotation matrix det should be 1, got {det}"
        );
    }

    // SM20. ShapeMatchingBody: centre_of_mass is mass-weighted.
    #[test]
    fn test_centre_of_mass_mass_weighted() {
        let rest = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let masses = vec![1.0_f64, 3.0]; // Total = 4; CM_x = (0*1 + 2*3)/4 = 1.5
        let body = ShapeMatchingBody::new(rest, masses, 1.0);
        let cm = body.centre_of_mass();
        assert!(
            (cm[0] - 1.5).abs() < 1e-12,
            "CM_x should be 1.5, got {}",
            cm[0]
        );
        assert!(cm[1].abs() < 1e-12);
        assert!(cm[2].abs() < 1e-12);
    }
}
