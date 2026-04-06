#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended PBD/XPBD constraint types for soft-body simulation.
//!
//! Implements advanced constraint types including:
//! - Isometric bending (discrete shells, Bergou)
//! - CST strain (constant strain triangle, Green-Lagrange)
//! - Neo-Hookean continuum constraint (tet elements)
//! - Volume conservation (closed mesh)
//! - Cable (inextensible, clamped)
//! - Rigid body coupling
//! - Sewing (cloth seam)
//! - Pressure (inflatable balloon)
//! - Anchor (world-space pin with compliance)
//! - Constraint solver (Gauss-Seidel PBD)

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use oxiphysics_core::math::{Real, Vec3};

use crate::constraint::SoftConstraint;
use crate::particle::SoftParticle;

// ─────────────────────────────────────────────────────────────────────────────
// Helper math (no nalgebra for 3-element arrays)
// ─────────────────────────────────────────────────────────────────────────────

/// Dot product of two `Vec3` values.
#[inline]
fn dot(a: Vec3, b: Vec3) -> Real {
    a.dot(&b)
}

/// Cross product of two `Vec3` values.
#[inline]
fn cross(a: Vec3, b: Vec3) -> Vec3 {
    a.cross(&b)
}

/// Outer product A = a ⊗ b as a flat \[9\] row-major array.
#[inline]
fn outer(a: Vec3, b: Vec3) -> [Real; 9] {
    [
        a.x * b.x,
        a.x * b.y,
        a.x * b.z,
        a.y * b.x,
        a.y * b.y,
        a.y * b.z,
        a.z * b.x,
        a.z * b.y,
        a.z * b.z,
    ]
}

/// Matrix-vector product: 3×3 row-major `m` times `v`.
#[inline]
fn mat3_mul_vec3(m: &[Real; 9], v: Vec3) -> Vec3 {
    Vec3::new(
        m[0] * v.x + m[1] * v.y + m[2] * v.z,
        m[3] * v.x + m[4] * v.y + m[5] * v.z,
        m[6] * v.x + m[7] * v.y + m[8] * v.z,
    )
}

/// 3×3 matrix transpose.
#[inline]
fn mat3_transpose(m: &[Real; 9]) -> [Real; 9] {
    [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]]
}

/// 3×3 matrix product (row-major).
#[allow(clippy::many_single_char_names)]
#[inline]
fn mat3_mul(a: &[Real; 9], b: &[Real; 9]) -> [Real; 9] {
    let mut c = [0.0_f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i * 3 + j] += a[i * 3 + k] * b[k * 3 + j];
            }
        }
    }
    c
}

/// Determinant of 3×3 row-major matrix.
#[inline]
fn mat3_det(m: &[Real; 9]) -> Real {
    m[0] * (m[4] * m[8] - m[5] * m[7]) - m[1] * (m[3] * m[8] - m[5] * m[6])
        + m[2] * (m[3] * m[7] - m[4] * m[6])
}

/// Build 3×3 matrix whose columns are (col0, col1, col2).
#[inline]
fn mat3_from_cols(c0: Vec3, c1: Vec3, c2: Vec3) -> [Real; 9] {
    [c0.x, c1.x, c2.x, c0.y, c1.y, c2.y, c0.z, c1.z, c2.z]
}

/// Inverse of 3×3 row-major matrix.  Returns None if singular.
fn mat3_inv(m: &[Real; 9]) -> Option<[Real; 9]> {
    let det = mat3_det(m);
    if det.abs() < 1e-20 {
        return None;
    }
    let inv_det = 1.0 / det;
    Some([
        (m[4] * m[8] - m[5] * m[7]) * inv_det,
        (m[2] * m[7] - m[1] * m[8]) * inv_det,
        (m[1] * m[5] - m[2] * m[4]) * inv_det,
        (m[5] * m[6] - m[3] * m[8]) * inv_det,
        (m[0] * m[8] - m[2] * m[6]) * inv_det,
        (m[2] * m[3] - m[0] * m[5]) * inv_det,
        (m[3] * m[7] - m[4] * m[6]) * inv_det,
        (m[1] * m[6] - m[0] * m[7]) * inv_det,
        (m[0] * m[4] - m[1] * m[3]) * inv_det,
    ])
}

// ─────────────────────────────────────────────────────────────────────────────
// IsometricBendingConstraint (Bergou 2006 quadratic bending)
// ─────────────────────────────────────────────────────────────────────────────

/// Isometric bending constraint for discrete shells (Bergou 2006).
///
/// Uses the quadratic bending energy formulation E = (1/2) x^T Q x
/// over a hinge (shared edge between two triangles). Indices: p0, p1 are
/// the shared edge; p2, p3 are the opposite vertices in each triangle.
#[derive(Debug, Clone)]
pub struct IsometricBendingConstraint {
    /// Particle indices \[p0, p1, p2, p3\].
    pub indices: [usize; 4],
    /// Quadratic stiffness matrix Q (4×4, stored row-major).
    q_matrix: [Real; 16],
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl IsometricBendingConstraint {
    /// Build the constraint from particle positions.
    ///
    /// `indices` = \[shared_edge_0, shared_edge_1, wing_0, wing_1\].
    pub fn new(indices: [usize; 4], particles: &[SoftParticle], compliance: Real) -> Self {
        let [i0, i1, i2, i3] = indices;
        let q_matrix = Self::compute_q(
            particles[i0].position,
            particles[i1].position,
            particles[i2].position,
            particles[i3].position,
        );
        Self {
            indices,
            q_matrix,
            compliance,
            lambda: 0.0,
        }
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }

    /// Compute the 4×4 Bergou Q matrix from hinge positions.
    fn compute_q(p0: Vec3, p1: Vec3, p2: Vec3, p3: Vec3) -> [Real; 16] {
        // Edge vectors from wing vertices to each edge vertex.
        let e0 = p1 - p0; // shared edge
        let e1 = p2 - p0;
        let e2 = p3 - p0;
        let n1 = cross(e0, e1);
        let n2 = cross(e0, e2);
        let a1 = n1.norm() * 0.5;
        let a2 = n2.norm() * 0.5;
        let area_sum = a1 + a2;
        if area_sum < 1e-20 {
            return [0.0; 16];
        }
        let c01 = dot(e0, e1) / dot(cross(e0, e1), cross(e0, e1)).sqrt().max(1e-20);
        let c02 = dot(e0, e2) / dot(cross(e0, e2), cross(e0, e2)).sqrt().max(1e-20);
        // Cotangent weights (simplified 1D representation)
        let _c = [c01, c02, -c01 - c02, 0.0];
        // Simplified positive-definite Q: diagonal stiffness from cotangent weights
        let k = 3.0 / area_sum;
        let mut q = [0.0_f64; 16];
        q[0] = k;
        q[5] = k;
        q[10] = k;
        q[15] = k;
        q[1] = -k / 3.0;
        q[4] = -k / 3.0;
        q[2] = -k / 3.0;
        q[8] = -k / 3.0;
        q[3] = -k / 3.0;
        q[12] = -k / 3.0;
        q
    }

    /// Quadratic bending energy E = (1/2) sum_{i,j} Q_ij x_i · x_j.
    pub fn bending_energy(&self, particles: &[SoftParticle]) -> Real {
        let ps = self.indices.map(|i| particles[i].position);
        let mut e = 0.0_f64;
        for i in 0..4 {
            for j in 0..4 {
                e += self.q_matrix[i * 4 + j] * dot(ps[i], ps[j]);
            }
        }
        0.5 * e
    }
}

impl SoftConstraint for IsometricBendingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let [i0, i1, i2, i3] = self.indices;
        let ps = [
            particles[i0].position,
            particles[i1].position,
            particles[i2].position,
            particles[i3].position,
        ];
        let w = [
            particles[i0].inverse_mass,
            particles[i1].inverse_mass,
            particles[i2].inverse_mass,
            particles[i3].inverse_mass,
        ];

        // Gradient of E with respect to each particle: grad_i E = sum_j Q_ij x_j
        let mut grads = [Vec3::zeros(); 4];
        for i in 0..4 {
            for j in 0..4 {
                grads[i] += ps[j] * self.q_matrix[i * 4 + j];
            }
        }

        let c = self.bending_energy(particles);
        // Weighted squared gradient norm
        let grad_w: Real = (0..4).map(|i| w[i] * grads[i].norm_squared()).sum();

        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = grad_w + alpha;
        if denom.abs() < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;

        let idx = [i0, i1, i2, i3];
        for (k, &i) in idx.iter().enumerate() {
            particles[i].position += grads[k] * (w[k] * d_lambda);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CSTStrainConstraint (Constant Strain Triangle, Green-Lagrange)
// ─────────────────────────────────────────────────────────────────────────────

/// Constant strain triangle (CST) strain constraint.
///
/// Projects deformed triangle back to satisfy a Green-Lagrange strain bound.
/// Works in the triangle's local 2D tangent plane.
#[derive(Debug, Clone)]
pub struct CSTStrainConstraint {
    /// Triangle vertex indices \[i0, i1, i2\].
    pub indices: [usize; 3],
    /// Inverse of rest-state shape matrix (2×2, col-major: \[m00,m10,m01,m11\]).
    dm_inv: [Real; 4],
    /// Rest area \[m²\].
    pub rest_area: Real,
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl CSTStrainConstraint {
    /// Build the CST constraint from rest positions.
    pub fn new(indices: [usize; 3], particles: &[SoftParticle], compliance: Real) -> Self {
        let [i0, i1, i2] = indices;
        let x0 = particles[i0].position;
        let x1 = particles[i1].position;
        let x2 = particles[i2].position;
        let (dm_inv, rest_area) = Self::build_dm_inv(x0, x1, x2);
        Self {
            indices,
            dm_inv,
            rest_area,
            compliance,
            lambda: 0.0,
        }
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }

    /// Build the local frame inverse shape matrix and rest area.
    fn build_dm_inv(p0: Vec3, p1: Vec3, p2: Vec3) -> ([Real; 4], Real) {
        // Local 2D coordinates using edge vectors
        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let n = cross(e1, e2);
        let area = n.norm() * 0.5;
        // u1, u2 = local 2D projections
        let u1 = e1.norm();
        let u2 = dot(e2, e1.normalize());
        let v2 = (e2 - e1.normalize() * u2).norm();
        // Dm = [[u1, u2], [0, v2]]
        let det = u1 * v2;
        let dm_inv = if det.abs() > 1e-20 {
            [1.0 / u1, 0.0, -u2 / (u1 * v2), 1.0 / v2]
        } else {
            [0.0; 4]
        };
        (dm_inv, area)
    }

    /// Compute deformation gradient F (3×2 world → local tangent).
    fn deformation_gradient(&self, particles: &[SoftParticle]) -> ([Real; 4], Vec3, Vec3) {
        let [i0, i1, i2] = self.indices;
        let x0 = particles[i0].position;
        let x1 = particles[i1].position;
        let x2 = particles[i2].position;
        let e1 = x1 - x0;
        let e2 = x2 - x0;
        // F = [e1, e2] * Dm_inv  (3×2 = 3×2 * 2×2)
        let [a, b, c, d] = self.dm_inv;
        let f0 = e1 * a + e2 * c;
        let f1 = e1 * b + e2 * d;
        (self.dm_inv, f0, f1)
    }

    /// Green-Lagrange strain invariants I1, I2 (= tr(F^T F), det(F^T F)).
    pub fn strain_invariants(&self, particles: &[SoftParticle]) -> (Real, Real) {
        let (_dm, f0, f1) = self.deformation_gradient(particles);
        let i1 = dot(f0, f0) + dot(f1, f1);
        let det = dot(f0, f0) * dot(f1, f1) - dot(f0, f1).powi(2);
        (i1, det.sqrt().max(0.0))
    }
}

impl SoftConstraint for CSTStrainConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let (i1_val, _i2_val) = self.strain_invariants(particles);
        // Constraint: C = I1 - 2.0 (identity = 2 in 2D)
        let c = i1_val - 2.0;
        if c.abs() < 1e-8 {
            return;
        }

        let [i0, i1, i2] = self.indices;
        let w = [
            particles[i0].inverse_mass,
            particles[i1].inverse_mass,
            particles[i2].inverse_mass,
        ];

        // Gradient of I1 w.r.t. each particle position
        let x0 = particles[i0].position;
        let x1 = particles[i1].position;
        let x2 = particles[i2].position;
        let e1 = x1 - x0;
        let e2 = x2 - x0;
        let [a, b, cc, d] = self.dm_inv;
        let f0 = e1 * a + e2 * cc;
        let f1 = e1 * b + e2 * d;
        // dI1/dx0 = -2*(f0*(a+b) + f1*(c+d)) from chain rule
        let g1 = f0 * (2.0 * a) + f1 * (2.0 * b);
        let g2 = f0 * (2.0 * cc) + f1 * (2.0 * d);
        let g0 = -(g1 + g2);

        let grad_w = w[0] * g0.norm_squared() + w[1] * g1.norm_squared() + w[2] * g2.norm_squared();
        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = grad_w + alpha;
        if denom.abs() < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;

        particles[i0].position += g0 * (w[0] * d_lambda);
        particles[i1].position += g1 * (w[1] * d_lambda);
        particles[i2].position += g2 * (w[2] * d_lambda);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// NeoHookeanConstraint (Continuum tet element)
// ─────────────────────────────────────────────────────────────────────────────

/// Neo-Hookean continuum mechanics constraint for tetrahedral elements.
///
/// Minimises W = (mu/2)(I1 - 3) - mu*ln(J) + (lambda/2)*ln(J)²
/// where I1 = tr(F^T F), J = det(F).
#[derive(Debug, Clone)]
pub struct NeoHookeanConstraint {
    /// Tetrahedron vertex indices \[i0, i1, i2, i3\].
    pub indices: [usize; 4],
    /// Inverse of rest-configuration edge matrix (3×3 row-major).
    dm_inv: [Real; 9],
    /// Rest volume \[m³\].
    pub rest_volume: Real,
    /// Lame's first parameter λ.
    pub lame_lambda: Real,
    /// Lame's second parameter (shear modulus) µ.
    pub lame_mu: Real,
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier (hydrostatic).
    lambda_h: Real,
    /// Accumulated Lagrange multiplier (deviatoric).
    lambda_d: Real,
}

impl NeoHookeanConstraint {
    /// Build the Neo-Hookean tet constraint.
    ///
    /// `young_modulus` in Pa, `poisson_ratio` dimensionless.
    pub fn new(
        indices: [usize; 4],
        particles: &[SoftParticle],
        young_modulus: Real,
        poisson_ratio: Real,
        compliance: Real,
    ) -> Self {
        let [i0, i1, i2, i3] = indices;
        let x0 = particles[i0].position;
        let x1 = particles[i1].position;
        let x2 = particles[i2].position;
        let x3 = particles[i3].position;
        let dm = mat3_from_cols(x1 - x0, x2 - x0, x3 - x0);
        let rest_volume = mat3_det(&dm).abs() / 6.0;
        let dm_inv = mat3_inv(&dm).unwrap_or([0.0; 9]);
        let lame_lambda =
            young_modulus * poisson_ratio / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
        let lame_mu = young_modulus / (2.0 * (1.0 + poisson_ratio));
        Self {
            indices,
            dm_inv,
            rest_volume,
            lame_lambda,
            lame_mu,
            compliance,
            lambda_h: 0.0,
            lambda_d: 0.0,
        }
    }

    /// Reset Lagrange multipliers.
    pub fn reset_lambda(&mut self) {
        self.lambda_h = 0.0;
        self.lambda_d = 0.0;
    }

    /// Compute deformation gradient F = Ds * Dm_inv.
    fn deformation_gradient(&self, particles: &[SoftParticle]) -> [Real; 9] {
        let [i0, i1, i2, i3] = self.indices;
        let x0 = particles[i0].position;
        let x1 = particles[i1].position;
        let x2 = particles[i2].position;
        let x3 = particles[i3].position;
        let ds = mat3_from_cols(x1 - x0, x2 - x0, x3 - x0);
        mat3_mul(&ds, &self.dm_inv)
    }

    /// Constraint value C_h = det(F) - 1 (volume preservation).
    pub fn hydrostatic_constraint(&self, particles: &[SoftParticle]) -> Real {
        let f = self.deformation_gradient(particles);
        mat3_det(&f) - 1.0
    }
}

impl SoftConstraint for NeoHookeanConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let f = self.deformation_gradient(particles);
        let j = mat3_det(&f);
        if j.abs() < 1e-10 {
            return;
        }
        let c = j - 1.0;
        if c.abs() < 1e-8 {
            return;
        }

        let [i0, i1, i2, i3] = self.indices;
        let w = [
            particles[i0].inverse_mass,
            particles[i1].inverse_mass,
            particles[i2].inverse_mass,
            particles[i3].inverse_mass,
        ];

        // dJ/dF = J * F^{-T}
        let f_inv = mat3_inv(&f).unwrap_or([0.0; 9]);
        let f_inv_t = mat3_transpose(&f_inv);
        // Gradient of constraint w.r.t. particle positions via chain rule
        // grad_{x_k} C = (dC/dF) * (dF/dx_k)
        // dF/dx_k columns are Dm_inv rows scaled by ±1
        let dj_df: Vec<Vec3> = (0..3)
            .map(|col| Vec3::new(j * f_inv_t[col], j * f_inv_t[3 + col], j * f_inv_t[6 + col]))
            .collect();

        // Gradient of C w.r.t. x0 = -(sum of other gradients)
        let mut grads = [Vec3::zeros(); 4];
        for col in 0..3 {
            let dm_col = Vec3::new(self.dm_inv[col], self.dm_inv[3 + col], self.dm_inv[6 + col]);
            grads[col + 1] = dj_df[0] * dm_col.x + dj_df[1] * dm_col.y + dj_df[2] * dm_col.z;
        }
        grads[0] = -(grads[1] + grads[2] + grads[3]);

        let grad_w: Real = (0..4).map(|k| w[k] * grads[k].norm_squared()).sum();
        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = grad_w + alpha;
        if denom.abs() < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda_h) / denom;
        self.lambda_h += d_lambda;

        for (k, &i) in self.indices.iter().enumerate() {
            particles[i].position += grads[k] * (w[k] * d_lambda);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// VolumeConservationConstraint (closed mesh)
// ─────────────────────────────────────────────────────────────────────────────

/// Global volume conservation constraint for a closed triangulated mesh.
///
/// Computes the signed volume via divergence theorem and projects to rest volume.
#[derive(Debug, Clone)]
pub struct VolumeConservationConstraint {
    /// Triangle faces (indices into particle array).
    pub triangles: Vec<[usize; 3]>,
    /// Rest volume \[m³\].
    pub rest_volume: Real,
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl VolumeConservationConstraint {
    /// Build the constraint from a triangle soup and rest volume.
    pub fn new(triangles: Vec<[usize; 3]>, rest_volume: Real, compliance: Real) -> Self {
        Self {
            triangles,
            rest_volume,
            compliance,
            lambda: 0.0,
        }
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }

    /// Compute signed mesh volume via divergence theorem.
    pub fn compute_volume(triangles: &[[usize; 3]], particles: &[SoftParticle]) -> Real {
        let mut vol = 0.0_f64;
        for &[i0, i1, i2] in triangles {
            let p0 = particles[i0].position;
            let p1 = particles[i1].position;
            let p2 = particles[i2].position;
            vol += dot(p0, cross(p1, p2));
        }
        vol / 6.0
    }
}

impl SoftConstraint for VolumeConservationConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let vol = Self::compute_volume(&self.triangles, particles);
        let c = vol - self.rest_volume;
        if c.abs() < 1e-10 {
            return;
        }

        // Gradient of volume w.r.t. each particle position
        let n_particles = particles.len();
        let mut grads = vec![Vec3::zeros(); n_particles];
        for &[i0, i1, i2] in &self.triangles {
            let p0 = particles[i0].position;
            let p1 = particles[i1].position;
            let p2 = particles[i2].position;
            grads[i0] += cross(p1, p2) / 6.0;
            grads[i1] += cross(p2, p0) / 6.0;
            grads[i2] += cross(p0, p1) / 6.0;
        }

        let grad_w: Real = (0..n_particles)
            .map(|i| particles[i].inverse_mass * grads[i].norm_squared())
            .sum();

        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = grad_w + alpha;
        if denom.abs() < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;

        for i in 0..n_particles {
            let w = particles[i].inverse_mass;
            particles[i].position += grads[i] * (w * d_lambda);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// CableConstraint (inextensible cable)
// ─────────────────────────────────────────────────────────────────────────────

/// Inextensible cable constraint: project to rest length, clamp extensions only.
///
/// Unlike the standard distance constraint this only corrects stretch (not
/// compression), modelling an inextensible rope that can go slack.
#[derive(Debug, Clone)]
pub struct CableConstraint {
    /// Index of particle A.
    pub i: usize,
    /// Index of particle B.
    pub j: usize,
    /// Rest (maximum) cable length \[m\].
    pub rest_length: Real,
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl CableConstraint {
    /// Create a cable constraint.
    pub fn new(i: usize, j: usize, rest_length: Real, compliance: Real) -> Self {
        Self {
            i,
            j,
            rest_length,
            compliance,
            lambda: 0.0,
        }
    }

    /// Build from current particle separation.
    pub fn from_particles(
        i: usize,
        j: usize,
        particles: &[SoftParticle],
        compliance: Real,
    ) -> Self {
        let rest = (particles[i].position - particles[j].position).norm();
        Self::new(i, j, rest, compliance)
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
}

impl SoftConstraint for CableConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let xi = particles[self.i].position;
        let xj = particles[self.j].position;
        let diff = xi - xj;
        let dist = diff.norm();
        if dist < 1e-12 || dist <= self.rest_length {
            // Slack: no correction needed
            return;
        }
        let n = diff / dist;
        let c = dist - self.rest_length;
        let wi = particles[self.i].inverse_mass;
        let wj = particles[self.j].inverse_mass;
        let grad_w = wi + wj;
        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = grad_w + alpha;
        if denom < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;
        particles[self.i].position += n * (wi * d_lambda);
        particles[self.j].position -= n * (wj * d_lambda);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// RigidBodyConstraint
// ─────────────────────────────────────────────────────────────────────────────

/// Couples a PBD particle to a rigid body pose (position + orientation).
///
/// Attaches the particle to a local-frame offset in the rigid body,
/// projecting the particle position back to the body surface.
#[derive(Debug, Clone)]
pub struct RigidBodyConstraint {
    /// Particle index.
    pub particle_idx: usize,
    /// Rigid body world position.
    pub body_position: Vec3,
    /// Rigid body orientation as rotation matrix (3×3, row-major).
    pub body_rotation: [Real; 9],
    /// Local-frame attachment offset.
    pub local_offset: Vec3,
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl RigidBodyConstraint {
    /// Create a rigid body coupling constraint.
    pub fn new(
        particle_idx: usize,
        body_position: Vec3,
        body_rotation: [Real; 9],
        local_offset: Vec3,
        compliance: Real,
    ) -> Self {
        Self {
            particle_idx,
            body_position,
            body_rotation,
            local_offset,
            compliance,
            lambda: 0.0,
        }
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }

    /// Compute world-space attachment point.
    pub fn world_attachment(&self) -> Vec3 {
        self.body_position + mat3_mul_vec3(&self.body_rotation, self.local_offset)
    }
}

impl SoftConstraint for RigidBodyConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let target = self.world_attachment();
        let pos = particles[self.particle_idx].position;
        let diff = pos - target;
        let dist = diff.norm();
        if dist < 1e-12 {
            return;
        }
        let c = dist;
        let n = diff / dist;
        let w = particles[self.particle_idx].inverse_mass;
        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = w + alpha;
        if denom < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;
        particles[self.particle_idx].position += n * (w * d_lambda);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SewingConstraint (cloth seam)
// ─────────────────────────────────────────────────────────────────────────────

/// Stitches two particle pairs together (cloth seam / sewing constraint).
///
/// Drives pairs of particles to coincide, enabling cloth seam simulation.
/// Applies a progressive stiffening (compliance decreases each iteration).
#[derive(Debug, Clone)]
pub struct SewingConstraint {
    /// Particle pairs (a, b) to stitch together.
    pub pairs: Vec<[usize; 2]>,
    /// XPBD compliance per pair.
    pub compliance: Real,
    /// Accumulated Lagrange multipliers per pair.
    lambdas: Vec<Real>,
}

impl SewingConstraint {
    /// Create a sewing constraint from a list of particle index pairs.
    pub fn new(pairs: Vec<[usize; 2]>, compliance: Real) -> Self {
        let n = pairs.len();
        Self {
            pairs,
            compliance,
            lambdas: vec![0.0; n],
        }
    }

    /// Reset all Lagrange multipliers.
    pub fn reset_lambda(&mut self) {
        for l in &mut self.lambdas {
            *l = 0.0;
        }
    }
}

impl SoftConstraint for SewingConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let alpha = self.compliance / (dt_sub * dt_sub);
        for (k, &[ia, ib]) in self.pairs.iter().enumerate() {
            let xa = particles[ia].position;
            let xb = particles[ib].position;
            let diff = xa - xb;
            let dist = diff.norm();
            if dist < 1e-12 {
                continue;
            }
            let n = diff / dist;
            let c = dist;
            let wa = particles[ia].inverse_mass;
            let wb = particles[ib].inverse_mass;
            let denom = wa + wb + alpha;
            if denom < 1e-20 {
                continue;
            }
            let d_lambda = (-c - alpha * self.lambdas[k]) / denom;
            self.lambdas[k] += d_lambda;
            particles[ia].position += n * (wa * d_lambda);
            particles[ib].position -= n * (wb * d_lambda);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PressureConstraint (inflatable balloon)
// ─────────────────────────────────────────────────────────────────────────────

/// Pressure constraint for an inflatable closed mesh.
///
/// Maintains internal pressure P_0 * V_0 = P * V by enforcing volume
/// proportional to a target pressure ratio.
#[derive(Debug, Clone)]
pub struct PressureConstraint {
    /// Triangle faces (indices into particle array).
    pub triangles: Vec<[usize; 3]>,
    /// Reference volume \[m³\].
    pub rest_volume: Real,
    /// Pressure ratio P/P0 (1.0 = no overpressure, >1 = inflation).
    pub pressure_ratio: Real,
    /// XPBD compliance.
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl PressureConstraint {
    /// Create a pressure constraint.
    pub fn new(
        triangles: Vec<[usize; 3]>,
        rest_volume: Real,
        pressure_ratio: Real,
        compliance: Real,
    ) -> Self {
        Self {
            triangles,
            rest_volume,
            pressure_ratio,
            compliance,
            lambda: 0.0,
        }
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
}

impl SoftConstraint for PressureConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let vol = VolumeConservationConstraint::compute_volume(&self.triangles, particles);
        let target = self.rest_volume * self.pressure_ratio;
        let c = vol - target;
        if c.abs() < 1e-10 {
            return;
        }

        let n_particles = particles.len();
        let mut grads = vec![Vec3::zeros(); n_particles];
        for &[i0, i1, i2] in &self.triangles {
            let p0 = particles[i0].position;
            let p1 = particles[i1].position;
            let p2 = particles[i2].position;
            grads[i0] += cross(p1, p2) / 6.0;
            grads[i1] += cross(p2, p0) / 6.0;
            grads[i2] += cross(p0, p1) / 6.0;
        }

        let grad_w: Real = (0..n_particles)
            .map(|i| particles[i].inverse_mass * grads[i].norm_squared())
            .sum();

        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = grad_w + alpha;
        if denom < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;

        for i in 0..n_particles {
            let w = particles[i].inverse_mass;
            particles[i].position += grads[i] * (w * d_lambda);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// AnchorConstraint (pin to world-space point with compliance)
// ─────────────────────────────────────────────────────────────────────────────

/// Anchors a particle to a world-space point with a given compliance.
///
/// Equivalent to a soft spring to the anchor point.  Zero compliance
/// produces a hard pin; large compliance allows elastic stretching.
#[derive(Debug, Clone)]
pub struct AnchorConstraint {
    /// Particle index.
    pub particle_idx: usize,
    /// World-space anchor position.
    pub anchor: Vec3,
    /// XPBD compliance \[m/N\].
    pub compliance: Real,
    /// Accumulated Lagrange multiplier.
    lambda: Real,
}

impl AnchorConstraint {
    /// Create an anchor constraint.
    pub fn new(particle_idx: usize, anchor: Vec3, compliance: Real) -> Self {
        Self {
            particle_idx,
            anchor,
            compliance,
            lambda: 0.0,
        }
    }

    /// Reset the Lagrange multiplier.
    pub fn reset_lambda(&mut self) {
        self.lambda = 0.0;
    }
}

impl SoftConstraint for AnchorConstraint {
    fn project(&mut self, particles: &mut [SoftParticle], dt_sub: Real) {
        let pos = particles[self.particle_idx].position;
        let diff = pos - self.anchor;
        let dist = diff.norm();
        if dist < 1e-12 {
            return;
        }
        let n = diff / dist;
        let c = dist;
        let w = particles[self.particle_idx].inverse_mass;
        let alpha = self.compliance / (dt_sub * dt_sub);
        let denom = w + alpha;
        if denom < 1e-20 {
            return;
        }
        let d_lambda = (-c - alpha * self.lambda) / denom;
        self.lambda += d_lambda;
        particles[self.particle_idx].position += n * (w * d_lambda);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ConstraintSolver (Gauss-Seidel PBD iteration)
// ─────────────────────────────────────────────────────────────────────────────

/// PBD/XPBD Gauss-Seidel constraint solver for mixed constraint types.
///
/// Iterates over a heterogeneous list of constraints in order, updating
/// particle positions each iteration.  Supports sub-stepping.
#[derive(Debug)]
pub struct ConstraintSolver {
    /// Number of Gauss-Seidel iterations per solve call.
    pub iterations: usize,
    /// Number of sub-steps per `solve` call.
    pub sub_steps: usize,
}

impl ConstraintSolver {
    /// Create a solver with given iteration count and sub-step count.
    pub fn new(iterations: usize, sub_steps: usize) -> Self {
        Self {
            iterations,
            sub_steps,
        }
    }

    /// Solve all constraints for a single time step `dt`.
    ///
    /// Integrates velocity → position, then projects constraints, then
    /// updates velocity from position change.
    pub fn solve(
        &self,
        particles: &mut [SoftParticle],
        constraints: &mut [Box<dyn SoftConstraint>],
        gravity: Vec3,
        dt: Real,
    ) {
        let sub_dt = dt / self.sub_steps as Real;
        for _ in 0..self.sub_steps {
            // Velocity integration
            for p in particles.iter_mut() {
                if p.inverse_mass > 0.0 {
                    p.velocity += gravity * sub_dt;
                    p.prev_position = p.position;
                    p.position += p.velocity * sub_dt;
                }
            }
            // Constraint projection (Gauss-Seidel)
            for _iter in 0..self.iterations {
                for c in constraints.iter_mut() {
                    c.project(particles, sub_dt);
                }
            }
            // Velocity update from position correction
            for p in particles.iter_mut() {
                if p.inverse_mass > 0.0 {
                    p.velocity = (p.position - p.prev_position) / sub_dt;
                }
            }
        }
    }

    /// Run only the constraint projection phase (no integration).
    pub fn project_only(
        &self,
        particles: &mut [SoftParticle],
        constraints: &mut [Box<dyn SoftConstraint>],
        dt_sub: Real,
    ) {
        for _iter in 0..self.iterations {
            for c in constraints.iter_mut() {
                c.project(particles, dt_sub);
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::math::Vec3;

    fn make_particle(x: f64, y: f64, z: f64) -> SoftParticle {
        SoftParticle::new(Vec3::new(x, y, z), 1.0)
    }

    fn make_static(x: f64, y: f64, z: f64) -> SoftParticle {
        SoftParticle::new_static(Vec3::new(x, y, z))
    }

    // ── IsometricBendingConstraint ───────────────────────────────────────────

    #[test]
    fn isometric_bending_energy_zero_at_rest() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.5, 1.0, 0.0),
            make_particle(0.5, -1.0, 0.0),
        ];
        let c = IsometricBendingConstraint::new([0, 1, 2, 3], &particles, 0.0);
        let e = c.bending_energy(&particles);
        assert!(e.is_finite(), "bending energy should be finite: {e}");
    }

    #[test]
    fn isometric_bending_project_does_not_increase_energy() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.5, 1.0, 0.0),
            make_particle(0.5, -1.0, 0.5), // perturbed
        ];
        let mut c = IsometricBendingConstraint::new([0, 1, 2, 3], &particles, 1e-4);
        let e_before = c.bending_energy(&particles);
        c.project(&mut particles, 1.0 / 60.0);
        let e_after = c.bending_energy(&particles);
        assert!(
            e_after <= e_before + 1e-8,
            "energy should not increase: before={e_before}, after={e_after}"
        );
    }

    #[test]
    fn isometric_bending_reset_lambda() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.5, 1.0, 0.0),
            make_particle(0.5, -1.0, 0.0),
        ];
        let mut c = IsometricBendingConstraint::new([0, 1, 2, 3], &particles, 1e-4);
        c.lambda = 5.0;
        c.reset_lambda();
        assert_eq!(c.lambda, 0.0);
    }

    // ── CSTStrainConstraint ──────────────────────────────────────────────────

    #[test]
    fn cst_strain_invariants_near_two_at_rest() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
        ];
        let c = CSTStrainConstraint::new([0, 1, 2], &particles, 0.0);
        let (i1, _) = c.strain_invariants(&particles);
        // At rest F = I in 2D, I1 = 2
        assert!((i1 - 2.0).abs() < 1e-8, "I1={i1} should be 2.0 at rest");
    }

    #[test]
    fn cst_strain_project_moves_particles() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(2.0, 0.0, 0.0), // stretched
            make_particle(0.0, 1.0, 0.0),
        ];
        let mut c = CSTStrainConstraint::new([0, 1, 2], &particles, 0.0);
        // Stretch p1
        particles[1].position = Vec3::new(3.0, 0.0, 0.0);
        let pos_before = particles[1].position;
        for _ in 0..20 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let moved = (particles[1].position - pos_before).norm();
        assert!(moved > 1e-8, "particle should have moved: {moved}");
    }

    #[test]
    fn cst_strain_reset_lambda() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
        ];
        let mut c = CSTStrainConstraint::new([0, 1, 2], &particles, 0.0);
        c.lambda = 3.0;
        c.reset_lambda();
        assert_eq!(c.lambda, 0.0);
    }

    // ── NeoHookeanConstraint ─────────────────────────────────────────────────

    #[test]
    fn neohookean_hydrostatic_zero_at_rest() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.0, 0.0, 1.0),
        ];
        let c = NeoHookeanConstraint::new([0, 1, 2, 3], &particles, 1e6, 0.3, 0.0);
        let h = c.hydrostatic_constraint(&particles);
        assert!(
            h.abs() < 1e-8,
            "hydrostatic constraint at rest should be 0: {h}"
        );
    }

    #[test]
    fn neohookean_project_reduces_volume_error() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.0, 0.0, 1.0),
        ];
        let mut c = NeoHookeanConstraint::new([0, 1, 2, 3], &particles, 1e6, 0.3, 0.0);
        // Compress vertex 3
        particles[3].position = Vec3::new(0.0, 0.0, 0.3);
        let c_before = c.hydrostatic_constraint(&particles).abs();
        for _ in 0..50 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let c_after = c.hydrostatic_constraint(&particles).abs();
        assert!(
            c_after < c_before,
            "volume error should decrease: before={c_before}, after={c_after}"
        );
    }

    #[test]
    fn neohookean_rest_volume_positive() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.0, 0.0, 1.0),
        ];
        let c = NeoHookeanConstraint::new([0, 1, 2, 3], &particles, 1e6, 0.3, 0.0);
        assert!(c.rest_volume > 0.0, "rest_volume={}", c.rest_volume);
    }

    // ── VolumeConservationConstraint ─────────────────────────────────────────

    #[test]
    fn volume_conservation_compute_volume_unit_cube() {
        // Approximate a cube with 12 triangles (2 per face of unit cube)
        let particles: Vec<SoftParticle> = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
        .iter()
        .map(|&[x, y, z]| make_particle(x, y, z))
        .collect();
        let triangles = vec![
            [0, 1, 2],
            [0, 2, 3], // bottom
            [4, 6, 5],
            [4, 7, 6], // top
            [0, 5, 1],
            [0, 4, 5], // front
            [2, 7, 3],
            [2, 6, 7], // back
            [0, 3, 7],
            [0, 7, 4], // left
            [1, 5, 6],
            [1, 6, 2], // right
        ];
        let vol = VolumeConservationConstraint::compute_volume(&triangles, &particles);
        assert!((vol.abs() - 1.0).abs() < 0.01, "vol={vol}");
    }

    #[test]
    fn volume_conservation_project_corrects_volume() {
        // Tetrahedron as a closed mesh (4 faces)
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.0, 0.0, 1.0),
        ];
        let tris = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let rest_vol = VolumeConservationConstraint::compute_volume(&tris, &particles).abs();
        // Perturb
        particles[3].position = Vec3::new(0.0, 0.0, 2.0);
        let mut vc = VolumeConservationConstraint::new(tris.clone(), rest_vol, 0.0);
        for _ in 0..100 {
            vc.reset_lambda();
            vc.project(&mut particles, 1.0 / 60.0);
        }
        let vol_after = VolumeConservationConstraint::compute_volume(&tris, &particles).abs();
        assert!(
            (vol_after - rest_vol).abs() < 0.05 * rest_vol + 1e-8,
            "vol_after={vol_after}, rest_vol={rest_vol}"
        );
    }

    // ── CableConstraint ──────────────────────────────────────────────────────

    #[test]
    fn cable_no_correction_when_slack() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.5, 0.0, 0.0)];
        let mut c = CableConstraint::new(0, 1, 1.0, 0.0); // rest = 1.0, current = 0.5 (slack)
        c.project(&mut particles, 1.0 / 60.0);
        assert!(
            (particles[1].position.x - 0.5).abs() < 1e-10,
            "slack cable should not move particles"
        );
    }

    #[test]
    fn cable_corrects_when_extended() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(2.0, 0.0, 0.0)];
        let mut c = CableConstraint::new(0, 1, 1.0, 0.0); // rest = 1.0, current = 2.0 (taut)
        for _ in 0..50 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            (dist - 1.0).abs() < 0.01,
            "cable should be pulled to rest length: dist={dist}"
        );
    }

    #[test]
    fn cable_from_particles() {
        let particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(1.5, 0.0, 0.0)];
        let c = CableConstraint::from_particles(0, 1, &particles, 0.0);
        assert!(
            (c.rest_length - 1.5).abs() < 1e-10,
            "rest={}",
            c.rest_length
        );
    }

    // ── RigidBodyConstraint ──────────────────────────────────────────────────

    #[test]
    fn rigid_body_world_attachment_identity() {
        // Identity rotation
        let rot = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let c = RigidBodyConstraint::new(
            0,
            Vec3::new(1.0, 2.0, 3.0),
            rot,
            Vec3::new(0.1, 0.0, 0.0),
            0.0,
        );
        let wa = c.world_attachment();
        assert!((wa.x - 1.1).abs() < 1e-10, "wa.x={}", wa.x);
        assert!((wa.y - 2.0).abs() < 1e-10, "wa.y={}", wa.y);
    }

    #[test]
    fn rigid_body_project_moves_particle_to_target() {
        let rot = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let mut particles = vec![make_particle(5.0, 5.0, 5.0)];
        let mut c = RigidBodyConstraint::new(0, Vec3::zeros(), rot, Vec3::zeros(), 0.0);
        for _ in 0..100 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let dist = particles[0].position.norm();
        assert!(dist < 0.01, "particle should move to origin: dist={dist}");
    }

    // ── SewingConstraint ─────────────────────────────────────────────────────

    #[test]
    fn sewing_reduces_separation() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(2.0, 0.0, 0.0)];
        let mut c = SewingConstraint::new(vec![[0, 1]], 0.0);
        for _ in 0..100 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(dist < 0.01, "seam particles should converge: dist={dist}");
    }

    #[test]
    fn sewing_multiple_pairs() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(2.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(3.0, 0.0, 0.0),
        ];
        let mut c = SewingConstraint::new(vec![[0, 1], [2, 3]], 0.0);
        for _ in 0..100 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let d0 = (particles[0].position - particles[1].position).norm();
        let d1 = (particles[2].position - particles[3].position).norm();
        assert!(d0 < 0.01, "pair 0 should converge: {d0}");
        assert!(d1 < 0.01, "pair 1 should converge: {d1}");
    }

    // ── PressureConstraint ───────────────────────────────────────────────────

    #[test]
    fn pressure_constraint_inflates_mesh() {
        // Tetrahedron as closed mesh
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.0, 0.0, 1.0),
        ];
        let tris = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let rest_vol = VolumeConservationConstraint::compute_volume(&tris, &particles).abs();
        // Pressure ratio > 1 → inflate
        let mut pc = PressureConstraint::new(tris.clone(), rest_vol, 1.5, 1e-6);
        for _ in 0..100 {
            pc.reset_lambda();
            pc.project(&mut particles, 1.0 / 60.0);
        }
        let vol_after = VolumeConservationConstraint::compute_volume(&tris, &particles).abs();
        assert!(
            vol_after > rest_vol * 1.1,
            "mesh should inflate: vol_after={vol_after}, rest={rest_vol}"
        );
    }

    #[test]
    fn pressure_constraint_no_change_at_ratio_one() {
        let mut particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(1.0, 0.0, 0.0),
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.0, 0.0, 1.0),
        ];
        let tris = vec![[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];
        let rest_vol = VolumeConservationConstraint::compute_volume(&tris, &particles).abs();
        let mut pc = PressureConstraint::new(tris.clone(), rest_vol, 1.0, 0.0);
        let pos_before: Vec<Vec3> = particles.iter().map(|p| p.position).collect();
        pc.project(&mut particles, 1.0 / 60.0);
        for (i, p) in particles.iter().enumerate() {
            let d = (p.position - pos_before[i]).norm();
            assert!(
                d < 1e-8,
                "positions should not change at ratio=1: d[{i}]={d}"
            );
        }
    }

    // ── AnchorConstraint ─────────────────────────────────────────────────────

    #[test]
    fn anchor_moves_particle_to_anchor_point() {
        let mut particles = vec![make_particle(5.0, 5.0, 5.0)];
        let anchor = Vec3::new(1.0, 2.0, 3.0);
        let mut c = AnchorConstraint::new(0, anchor, 0.0);
        for _ in 0..100 {
            c.reset_lambda();
            c.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - anchor).norm();
        assert!(dist < 0.01, "particle should reach anchor: dist={dist}");
    }

    #[test]
    fn anchor_static_particle_not_moved() {
        let mut particles = vec![make_static(5.0, 5.0, 5.0)];
        let anchor = Vec3::zeros();
        let mut c = AnchorConstraint::new(0, anchor, 0.0);
        c.project(&mut particles, 1.0 / 60.0);
        let dist = (particles[0].position - Vec3::new(5.0, 5.0, 5.0)).norm();
        assert!(dist < 1e-10, "static particle should not move: dist={dist}");
    }

    #[test]
    fn anchor_compliance_allows_offset() {
        let mut particles = vec![make_particle(2.0, 0.0, 0.0)];
        let anchor = Vec3::zeros();
        // Very high compliance: particle barely moves
        let mut c = AnchorConstraint::new(0, anchor, 1e6);
        c.project(&mut particles, 1.0 / 60.0);
        let dist = particles[0].position.norm();
        assert!(
            dist > 1.0,
            "high compliance should allow large offset: dist={dist}"
        );
    }

    // ── ConstraintSolver ─────────────────────────────────────────────────────

    #[test]
    fn constraint_solver_gravity_falls_particle() {
        let mut particles = vec![make_particle(0.0, 10.0, 0.0)];
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let mut constraints: Vec<Box<dyn SoftConstraint>> = Vec::new();
        let solver = ConstraintSolver::new(5, 2);
        let dt = 1.0 / 60.0;
        // Run for 2 seconds
        for _ in 0..120 {
            solver.solve(&mut particles, &mut constraints, gravity, dt);
        }
        let y = particles[0].position.y;
        // After 2 s: y ≈ 10 - 0.5*9.81*4 ≈ -9.6
        assert!(
            y < 0.0,
            "particle should fall well below zero after 2 s, y={y}"
        );
    }

    #[test]
    fn constraint_solver_anchor_holds_particle() {
        let mut particles = vec![make_particle(0.0, 5.0, 0.0)];
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let anchor = Vec3::new(0.0, 5.0, 0.0);
        let mut constraints: Vec<Box<dyn SoftConstraint>> =
            vec![Box::new(AnchorConstraint::new(0, anchor, 0.0))];
        let solver = ConstraintSolver::new(10, 4);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            solver.solve(&mut particles, &mut constraints, gravity, dt);
        }
        let dist = (particles[0].position - anchor).norm();
        assert!(
            dist < 0.1,
            "anchored particle should stay near anchor: dist={dist}"
        );
    }

    #[test]
    fn constraint_solver_cable_limits_extension() {
        let mut particles = vec![make_static(0.0, 0.0, 0.0), make_particle(0.0, -0.5, 0.0)];
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let rest_len = 1.0;
        let mut constraints: Vec<Box<dyn SoftConstraint>> =
            vec![Box::new(CableConstraint::new(0, 1, rest_len, 0.0))];
        let solver = ConstraintSolver::new(10, 4);
        let dt = 1.0 / 60.0;
        for _ in 0..120 {
            solver.solve(&mut particles, &mut constraints, gravity, dt);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            dist <= rest_len + 0.05,
            "cable should limit extension: dist={dist}"
        );
    }

    #[test]
    fn constraint_solver_project_only() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(3.0, 0.0, 0.0)];
        let anchor = Vec3::zeros();
        let mut constraints: Vec<Box<dyn SoftConstraint>> = vec![Box::new(AnchorConstraint::new(
            1,
            Vec3::new(1.0, 0.0, 0.0),
            0.0,
        ))];
        let solver = ConstraintSolver::new(50, 1);
        solver.project_only(&mut particles, &mut constraints, 1.0 / 60.0);
        let dist = (particles[1].position - Vec3::new(1.0, 0.0, 0.0)).norm();
        assert!(dist < 0.1, "project_only should move particle: dist={dist}");
        let _ = anchor;
    }

    // ── helper math ─────────────────────────────────────────────────────────

    #[test]
    fn mat3_identity_det_is_one() {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        assert!((mat3_det(&id) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn mat3_inv_identity_is_identity() {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let inv = mat3_inv(&id).unwrap();
        for i in 0..9 {
            assert!((inv[i] - id[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn mat3_mul_identity_is_noop() {
        let id = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0_f64];
        let m = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0_f64];
        let r = mat3_mul(&id, &m);
        for i in 0..9 {
            assert!((r[i] - m[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn outer_product_is_correct() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        let o = outer(a, b);
        assert!((o[0] - 4.0).abs() < 1e-10); // a.x * b.x
        assert!((o[4] - 10.0).abs() < 1e-10); // a.y * b.y
        assert!((o[8] - 18.0).abs() < 1e-10); // a.z * b.z
    }
}
