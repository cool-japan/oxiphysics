// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spatial vectors (6D) for Featherstone articulated-body dynamics.
//!
//! Following Featherstone 2008 "Rigid Body Dynamics Algorithms" notation.
//!
//! Spatial velocity:  v = [ω; v_linear]   (angular first, then linear)
//! Spatial force:     f = [n; f_linear]   (moment first, then force)
//! Spatial inertia:   I (6×6 symmetric positive-definite)
//!
//! The [`SpatialInertia`] stores the rotational inertia **about the body-frame
//! origin** (not about the centre of mass), so the spatial-inertia multiply
//! formula `I·v` takes the simple block-matrix form without extra cross-product
//! terms. Use [`SpatialInertia::from_com`] to construct from per-COM data; it
//! applies the parallel-axis theorem internally.

use std::ops::{Add, Mul, Neg, Sub};

// ─── helpers ─────────────────────────────────────────────────────────────────

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(s: f64, a: [f64; 3]) -> [f64; 3] {
    [s * a[0], s * a[1], s * a[2]]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// 3×3 matrix multiply: result[i][j] = sum_k A[i][k] * B[k][j]
fn mat3_mul(a: [[f64; 3]; 3], b: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut out = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                out[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    out
}

// 3×3 matrix times 3-vector
fn mat3_vec(m: [[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

// Transpose of a 3×3 matrix
fn mat3_transpose(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Skew-symmetric matrix of a vector: cross(a, v) = skew(a) * v
fn skew(a: [f64; 3]) -> [[f64; 3]; 3] {
    [[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]]
}

// ─── SpatialVec ──────────────────────────────────────────────────────────────

/// 6D spatial vector (angular part first, then linear part).
///
/// Used for both motion vectors (velocity, acceleration) and force vectors
/// (moment + force). Context determines the interpretation.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpatialVec {
    /// Angular component [ωx, ωy, ωz] (or moment [nx, ny, nz] for forces).
    pub angular: [f64; 3],
    /// Linear component [vx, vy, vz] (or force [fx, fy, fz] for forces).
    pub linear: [f64; 3],
}

impl SpatialVec {
    /// The zero spatial vector.
    pub const ZERO: Self = Self {
        angular: [0.0, 0.0, 0.0],
        linear: [0.0, 0.0, 0.0],
    };

    /// Construct from angular and linear parts.
    #[inline]
    pub fn new(angular: [f64; 3], linear: [f64; 3]) -> Self {
        Self { angular, linear }
    }

    /// Dot product of two spatial vectors: ω·ω' + v·v'.
    #[inline]
    pub fn dot(&self, other: &Self) -> f64 {
        dot3(self.angular, other.angular) + dot3(self.linear, other.linear)
    }

    /// Spatial cross product of *motion* vectors: `self × other`
    ///
    /// Used for velocity propagation. Given `v = [ω; v_lin]` and `u = [ω'; v']`:
    ///
    /// ```text
    /// v × u = [ ω × ω' ;  ω × v' + v_lin × ω' ]
    /// ```
    #[inline]
    pub fn cross_motion(&self, other: &Self) -> Self {
        let ang = cross3(self.angular, other.angular);
        let lin = add3(
            cross3(self.angular, other.linear),
            cross3(self.linear, other.angular),
        );
        Self::new(ang, lin)
    }

    /// Spatial cross product for *force* vectors: `self ×* force`
    ///
    /// Given `v = [ω; v_lin]` and `f = [n; f_lin]`:
    ///
    /// ```text
    /// v ×* f = [ ω × n + v_lin × f_lin ;  ω × f_lin ]
    /// ```
    #[inline]
    pub fn cross_force(&self, force: &Self) -> Self {
        let ang = add3(
            cross3(self.angular, force.angular),
            cross3(self.linear, force.linear),
        );
        let lin = cross3(self.angular, force.linear);
        Self::new(ang, lin)
    }

    /// Scale by a scalar.
    #[inline]
    pub fn scale(&self, s: f64) -> Self {
        Self::new(scale3(s, self.angular), scale3(s, self.linear))
    }
}

impl Add for SpatialVec {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::new(
            add3(self.angular, rhs.angular),
            add3(self.linear, rhs.linear),
        )
    }
}

impl Sub for SpatialVec {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::new(
            sub3(self.angular, rhs.angular),
            sub3(self.linear, rhs.linear),
        )
    }
}

impl Neg for SpatialVec {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::new(
            [-self.angular[0], -self.angular[1], -self.angular[2]],
            [-self.linear[0], -self.linear[1], -self.linear[2]],
        )
    }
}

impl Mul<f64> for SpatialVec {
    type Output = Self;
    #[inline]
    fn mul(self, s: f64) -> Self {
        self.scale(s)
    }
}

// ─── SpatialInertia ──────────────────────────────────────────────────────────

/// Spatial rigid-body inertia (6×6 symmetric positive-definite).
///
/// Internally stored as the block-matrix form:
///
/// ```text
/// I = [ I_origin    m·c× ]
///     [ m·c×ᵀ       m·E  ]
/// ```
///
/// where `I_origin` is the 3×3 rotational inertia **about the body-frame
/// origin** (not the COM), `c` is the COM position in the body frame, and `m`
/// is the mass. This representation makes `I·v` a simple block-multiply.
///
/// Use `from_com` to construct from (mass, com, inertia-about-com).
#[derive(Debug, Clone, Copy)]
pub struct SpatialInertia {
    /// Mass \[kg\].
    pub mass: f64,
    /// Centre-of-mass position in body frame \[m\].
    pub com: [f64; 3],
    /// Rotational inertia about the **body-frame origin** (3×3, row-major).
    pub rot_inertia_origin: [[f64; 3]; 3],
}

impl SpatialInertia {
    /// Construct from physical parameters (inertia **about COM**).
    ///
    /// Applies the parallel-axis theorem to convert to origin-frame storage:
    /// `I_origin = I_com + m (c·c E − c cᵀ)`.
    pub fn from_com(mass: f64, com: [f64; 3], rot_inertia_com: [[f64; 3]; 3]) -> Self {
        // Parallel-axis: I_origin = I_com + m * (|c|² I - c cᵀ)
        let cc = dot3(com, com); // c · c
        let mut io = rot_inertia_com;
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                io[i][j] += mass * (cc * delta - com[i] * com[j]);
            }
        }
        Self {
            mass,
            com,
            rot_inertia_origin: io,
        }
    }

    /// Construct directly from (mass, com, inertia-about-origin) — lower-level.
    pub fn from_origin(mass: f64, com: [f64; 3], rot_inertia_origin: [[f64; 3]; 3]) -> Self {
        Self {
            mass,
            com,
            rot_inertia_origin,
        }
    }

    /// Multiply: I · v → spatial force.
    ///
    /// Block-matrix formula (using `h = m·c`):
    ///
    /// ```text
    /// [n]   [ I_origin   h× ] [ω]
    /// [f] = [ h×ᵀ        mE ] [v]
    ///
    ///  n = I_origin · ω + h × v
    ///  f = m · v − h × ω      (h×ᵀ · ω = −h × ω, but sign flips for forces)
    /// ```
    ///
    /// More explicitly: `f = m·v + m·(−c × ω) = m·(v − c × ω)`.
    pub fn mul_vec(&self, v: &SpatialVec) -> SpatialVec {
        let h = scale3(self.mass, self.com); // h = m·c
        // angular (moment): I_origin · ω + h × v_lin
        let i_omega = mat3_vec(self.rot_inertia_origin, v.angular);
        let h_cross_v = cross3(h, v.linear);
        let ang = add3(i_omega, h_cross_v);
        // linear (force): m·v_lin - h × ω  (= m·v_lin + m·c × ω, but note
        //   h×ᵀ · ω ≠ h × ω; block (2,1) = −h× so term = −h × ω = h × ω with minus)
        // Correct: f = m*v_lin + (-c) × (m*ω) = m*v_lin - h×ω
        let h_cross_omega = cross3(h, v.angular);
        let lin = sub3(scale3(self.mass, v.linear), h_cross_omega);
        SpatialVec::new(ang, lin)
    }

    /// Add two spatial inertias (used in ABA articulated inertia accumulation).
    pub fn add(&self, other: &Self) -> Self {
        let mass = self.mass + other.mass;
        // Combined COM (weighted)
        let com = if mass < 1e-30 {
            [0.0; 3]
        } else {
            let a = scale3(self.mass, self.com);
            let b = scale3(other.mass, other.com);
            scale3(1.0 / mass, add3(a, b))
        };
        // Combined inertia about origin (simple addition in this representation)
        let mut io = self.rot_inertia_origin;
        for (io_row, other_row) in io.iter_mut().zip(other.rot_inertia_origin.iter()) {
            for (io_cell, other_cell) in io_row.iter_mut().zip(other_row.iter()) {
                *io_cell += other_cell;
            }
        }
        Self {
            mass,
            com,
            rot_inertia_origin: io,
        }
    }
}

// ─── SpatialTransform ────────────────────────────────────────────────────────

/// Spatial coordinate transform (Plücker transform) between two frames.
///
/// A `SpatialTransform` stores the **pose** of a child frame `C` relative to a
/// parent frame `P`:
///
/// - `rot` (`R`) maps child coordinates to parent coordinates
///   (`x_P = R · x_C + r`),
/// - `trans` (`r`) is the child origin expressed in parent coordinates.
///
/// The corresponding Plücker *motion* transform from parent to child
/// coordinates (`C X_P` in Featherstone's notation) uses `E = Rᵀ`:
///
/// ```text
/// C X_P · [ω; v] = [ Rᵀ·ω ;  Rᵀ·(v − r×ω) ]
/// ```
///
/// and the matching *force* transform from child to parent coordinates
/// (`P X_C^* = (C X_P)ᵀ`) is
///
/// ```text
/// (C X_P)ᵀ · [n; f] = [ R·n + r×(R·f) ;  R·f ]
/// ```
///
/// Composition ([`SpatialTransform::compose`]) chains poses, so
/// `world_pose_of_child = world_pose_of_parent.compose(&parent_pose_of_child)`.
#[derive(Debug, Clone, Copy)]
pub struct SpatialTransform {
    /// Rotation matrix from child to parent coordinates (row-major, R\[i\]\[j\]).
    pub rot: [[f64; 3]; 3],
    /// Translation from parent origin to child origin, in parent coordinates \[m\].
    pub trans: [f64; 3],
}

impl SpatialTransform {
    /// Identity transform (same frame).
    pub const IDENTITY: Self = Self {
        rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        trans: [0.0, 0.0, 0.0],
    };

    /// Construct from rotation matrix and translation.
    pub fn from_rotation_translation(rot: [[f64; 3]; 3], trans: [f64; 3]) -> Self {
        Self { rot, trans }
    }

    /// Apply the Plücker motion transform `C X_P`: express a parent-frame
    /// motion vector in child coordinates.
    ///
    /// `ω' = Rᵀ · ω`
    /// `v' = Rᵀ · (v − r × ω)`
    pub fn apply_velocity(&self, v: &SpatialVec) -> SpatialVec {
        let rt = mat3_transpose(self.rot);
        let omega_child = mat3_vec(rt, v.angular);
        let r_cross_omega = cross3(self.trans, v.angular);
        let v_minus = sub3(v.linear, r_cross_omega);
        let v_child = mat3_vec(rt, v_minus);
        SpatialVec::new(omega_child, v_child)
    }

    /// Apply the force transform `(C X_P)ᵀ = P X_C^*`: express a child-frame
    /// force vector in parent coordinates (used by the inward RNEA/ABA passes).
    ///
    /// `f' = R · f`
    /// `n' = R · n + r × (R · f)`
    pub fn apply_force(&self, f: &SpatialVec) -> SpatialVec {
        let f_parent = mat3_vec(self.rot, f.linear);
        let n_rotated = mat3_vec(self.rot, f.angular);
        let r_cross_f = cross3(self.trans, f_parent);
        let n_parent = add3(n_rotated, r_cross_f);
        SpatialVec::new(n_parent, f_parent)
    }

    /// Compose transforms: `self ∘ other` (apply `other` first, then `self`).
    ///
    /// For chained joints: parent←child1 composed with child1←child2.
    pub fn compose(&self, other: &Self) -> Self {
        // Rotation: R = R_self · R_other
        let rot = mat3_mul(self.rot, other.rot);
        // Translation: r = r_self + R_self · r_other
        let r_other_in_self = mat3_vec(self.rot, other.trans);
        let trans = add3(self.trans, r_other_in_self);
        Self { rot, trans }
    }

    /// Compute the inverse transform (from child back to parent reversed).
    ///
    /// If X transforms from B to A, then X.inverse() transforms from A to B:
    ///
    /// `R_inv = R^T`
    /// `r_inv = −R^T · r`
    pub fn inverse(&self) -> Self {
        let rt = mat3_transpose(self.rot);
        let neg_r = [-self.trans[0], -self.trans[1], -self.trans[2]];
        let trans_inv = mat3_vec(rt, neg_r);
        Self {
            rot: rt,
            trans: trans_inv,
        }
    }

    /// Build a rotation-only transform from axis-angle (Rodrigues' formula).
    pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if n < 1e-14 {
            return Self::IDENTITY;
        }
        let u = [axis[0] / n, axis[1] / n, axis[2] / n];
        let c = angle.cos();
        let s = angle.sin();
        let t = 1.0 - c;
        let rot = [
            [
                t * u[0] * u[0] + c,
                t * u[0] * u[1] - s * u[2],
                t * u[0] * u[2] + s * u[1],
            ],
            [
                t * u[0] * u[1] + s * u[2],
                t * u[1] * u[1] + c,
                t * u[1] * u[2] - s * u[0],
            ],
            [
                t * u[0] * u[2] - s * u[1],
                t * u[1] * u[2] + s * u[0],
                t * u[2] * u[2] + c,
            ],
        ];
        Self {
            rot,
            trans: [0.0; 3],
        }
    }

    /// Build a translation-only transform.
    pub fn from_translation(trans: [f64; 3]) -> Self {
        Self {
            rot: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            trans,
        }
    }

    /// Express a child-frame spatial inertia in parent coordinates.
    ///
    /// The mass is unchanged, the COM becomes `c' = R·c + r`, and the
    /// rotational inertia about the parent origin follows from the
    /// parallel-axis theorem applied about the COM:
    ///
    /// ```text
    /// I'_O = R·I_O·Rᵀ + m·[ (|r|² + 2 r·c')·E − r rᵀ − r c'ᵀ − c' rᵀ ]
    /// ```
    ///
    /// This equals `(C X_P)ᵀ · I · C X_P` (see
    /// [`SpatialInertia6x6::transform_to_parent`]).
    pub fn transform_inertia(&self, inertia: &SpatialInertia) -> SpatialInertia {
        let m = inertia.mass;
        let c_rot = mat3_vec(self.rot, inertia.com);
        let new_com = add3(c_rot, self.trans);
        let rt = mat3_transpose(self.rot);
        let mut io = mat3_mul(mat3_mul(self.rot, inertia.rot_inertia_origin), rt);
        let t = self.trans;
        let diag = dot3(t, t) + 2.0 * dot3(t, c_rot);
        for (i, io_row) in io.iter_mut().enumerate() {
            for (j, io_ij) in io_row.iter_mut().enumerate() {
                let delta = if i == j { 1.0 } else { 0.0 };
                *io_ij += m * (diag * delta - t[i] * t[j] - t[i] * c_rot[j] - c_rot[i] * t[j]);
            }
        }
        SpatialInertia {
            mass: m,
            com: new_com,
            rot_inertia_origin: io,
        }
    }
}

// ─── Full 6×6 Spatial Inertia for ABA ────────────────────────────────────────

/// Full 6×6 spatial inertia matrix for use in the ABA algorithm.
///
/// ABA requires rank-1 updates and direct linear-system solves on the full
/// 6×6 block. This type stores all 36 entries (symmetric, so only 21 unique)
/// and supports the required operations.
#[derive(Debug, Clone, Copy)]
pub struct SpatialInertia6x6 {
    /// Row-major 6×6 symmetric matrix.
    pub data: [[f64; 6]; 6],
}

impl SpatialInertia6x6 {
    /// Construct from a `SpatialInertia` (packed form → full 6×6).
    pub fn from_packed(si: &SpatialInertia) -> Self {
        let m = si.mass;
        let c = si.com;
        let io = si.rot_inertia_origin;
        // Build h = m·c and h× (skew)
        let h = scale3(m, c);
        let hx = skew(h);
        // Full 6×6:
        // [ I_origin   h×  ]
        // [ -h×        m·I ]
        // Note: (h×)^T = -h× so bottom-left = -h× = (h×)^T
        let mut d = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = io[i][j];
                d[i][j + 3] = hx[i][j];
                d[i + 3][j] = -hx[i][j]; // -h× = (h×)^T
                d[i + 3][j + 3] = if i == j { m } else { 0.0 };
            }
        }
        Self { data: d }
    }

    /// Multiply 6×6 inertia by a spatial vector (returns spatial force).
    pub fn mul_vec6(&self, v: &SpatialVec) -> SpatialVec {
        let vv = [
            v.angular[0],
            v.angular[1],
            v.angular[2],
            v.linear[0],
            v.linear[1],
            v.linear[2],
        ];
        let mut result = [0.0f64; 6];
        for (result_i, data_row) in result.iter_mut().zip(self.data.iter()) {
            for (data_ij, vv_j) in data_row.iter().zip(vv.iter()) {
                *result_i += data_ij * vv_j;
            }
        }
        SpatialVec::new(
            [result[0], result[1], result[2]],
            [result[3], result[4], result[5]],
        )
    }

    /// Add two 6×6 spatial inertias.
    pub fn add(&self, other: &Self) -> Self {
        let mut d = [[0.0f64; 6]; 6];
        for (d_row, (self_row, other_row)) in
            d.iter_mut().zip(self.data.iter().zip(other.data.iter()))
        {
            for (d_cell, (self_cell, other_cell)) in
                d_row.iter_mut().zip(self_row.iter().zip(other_row.iter()))
            {
                *d_cell = self_cell + other_cell;
            }
        }
        Self { data: d }
    }

    /// Subtract outer product: `I -= u * u^T * d_inv` (rank-1 update for ABA).
    pub fn sub_rank1(&self, u: &SpatialVec, d_inv: f64) -> Self {
        let uu = [
            u.angular[0],
            u.angular[1],
            u.angular[2],
            u.linear[0],
            u.linear[1],
            u.linear[2],
        ];
        let mut d = self.data;
        for i in 0..6 {
            for j in 0..6 {
                d[i][j] -= uu[i] * uu[j] * d_inv;
            }
        }
        Self { data: d }
    }

    /// Express a child-frame articulated/spatial inertia in parent coordinates.
    ///
    /// `x` is the pose of the child frame relative to the parent (see
    /// [`SpatialTransform`]). With the parent→child motion transform
    ///
    /// ```text
    /// X = C X_P = [  Rᵀ       0  ]
    ///             [ −Rᵀ·r×    Rᵀ ]
    /// ```
    ///
    /// the inertia seen from the parent is `I_parent = Xᵀ · I_child · X`.
    pub fn transform_to_parent(&self, x: &SpatialTransform) -> Self {
        let e = mat3_transpose(x.rot); // E = Rᵀ (parent → child rotation)
        let rx = skew(x.trans); // r× = skew(r)

        // neg_e_rx = −E · r× (lower-left block of X)
        let mut neg_e_rx = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                neg_e_rx[i][j] = -(e[i][0] * rx[0][j] + e[i][1] * rx[1][j] + e[i][2] * rx[2][j]);
            }
        }

        let mut x6 = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                x6[i][j] = e[i][j];
                x6[i + 3][j] = neg_e_rx[i][j];
                x6[i + 3][j + 3] = e[i][j];
            }
        }

        // Compute I_parent = X^T * I_child * X
        // Step 1: tmp = I_child * X
        let mut tmp = [[0.0f64; 6]; 6];
        for (tmp_row, self_row) in tmp.iter_mut().zip(self.data.iter()) {
            for j in 0..6 {
                for (self_ik, x6_kj) in self_row.iter().zip(x6.iter().map(|x6k| &x6k[j])) {
                    tmp_row[j] += self_ik * x6_kj;
                }
            }
        }
        // Step 2: out = X^T * tmp
        let mut out = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    out[i][j] += x6[k][i] * tmp[k][j]; // x6[k][i] = X^T[i][k]
                }
            }
        }
        Self { data: out }
    }
}

// ─── Unit tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    #[test]
    fn test_cross_motion() {
        // v = [1,0,0; 0,0,0], u = [0,1,0; 0,0,0]
        // v × u = [1×[0,1,0]; 1×[0,0,0] + [0,0,0]×[0,1,0]]
        //       = [ω×ω'; ω×v' + v×ω']
        //       = [[0,0,-1] × [0,1,0]; [0,0,-1]..] -- recalc with correct values
        let v = SpatialVec::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let u = SpatialVec::new([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        let result = v.cross_motion(&u);
        // ω × ω' = [1,0,0] × [0,1,0] = [0,0,1]
        assert!(
            approx_eq(result.angular[2], 1.0, 1e-12),
            "angular z should be 1, got {:.6}",
            result.angular[2]
        );
        assert!(
            approx_eq(result.linear[0], 0.0, 1e-12),
            "linear should be zero"
        );
    }

    #[test]
    fn test_cross_force() {
        // v = [0,0,1; 0,0,0], f = [0,0,0; 1,0,0]
        // v ×* f = [ω×n + v×f; ω×f]
        //        = [[0,0,1]×[0,0,0] + [0,0,0]×[1,0,0]; [0,0,1]×[1,0,0]]
        //        = [[0;0;0]; [0,1,0]]
        // Wait: [0,0,1]×[1,0,0] = [0*0-1*0, 1*1-0*0, 0*0-0*1] = [0,1,0]
        let v = SpatialVec::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let f = SpatialVec::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let result = v.cross_force(&f);
        assert!(
            approx_eq(result.angular[0], 0.0, 1e-12),
            "ang[0] should be 0"
        );
        assert!(
            approx_eq(result.linear[1], 1.0, 1e-12),
            "lin[1] should be 1, got {:.6}",
            result.linear[1]
        );
    }

    #[test]
    fn test_inertia_point_mass() {
        // Point mass m=1 at c=[L,0,0] — verify angular inertia is m*L^2 about z-axis
        // I_com = 0 (point mass), so I_origin = m*(|c|^2 I - c·c^T)
        // For rotation about z (ω=[0,0,1]), v_lin=0:
        // I_origin * ω → should give moment = m*L^2 * [0,0,1]
        let l = 2.0;
        let m = 3.0;
        let si = SpatialInertia::from_com(m, [l, 0.0, 0.0], [[0.0; 3]; 3]);
        let v = SpatialVec::new([0.0, 0.0, 1.0], [0.0, 0.0, 0.0]);
        let fv = si.mul_vec(&v);
        let expected_ang_z = m * l * l;
        assert!(
            approx_eq(fv.angular[2], expected_ang_z, 1e-10),
            "point mass I*ω_z should be m*L^2={:.4}, got {:.4}",
            expected_ang_z,
            fv.angular[2]
        );
        assert!(
            approx_eq(fv.linear[0], 0.0, 1e-10),
            "linear force for pure angular velocity and point mass at x should be zero"
        );
    }

    #[test]
    fn test_transform_identity() {
        let v = SpatialVec::new([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        let xt = SpatialTransform::IDENTITY;
        let vt = xt.apply_velocity(&v);
        assert!(
            approx_eq(vt.angular[0], v.angular[0], 1e-12),
            "identity should preserve angular"
        );
        assert!(
            approx_eq(vt.linear[2], v.linear[2], 1e-12),
            "identity should preserve linear"
        );
    }
}
#[cfg(test)]
mod transform_tests {
    use super::*;

    fn approx_eq_6x6(a: &SpatialInertia6x6, b: &SpatialInertia6x6, eps: f64, label: &str) {
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (a.data[i][j] - b.data[i][j]).abs() < eps,
                    "{label}: data[{i}][{j}] expected {:.8}, got {:.8}",
                    b.data[i][j],
                    a.data[i][j]
                );
            }
        }
    }

    /// transform_to_parent on identity should return itself.
    #[test]
    fn test_transform_to_parent_identity() {
        let si = SpatialInertia::from_com(2.0, [0.3, 0.1, 0.0], [[0.0; 3]; 3]);
        let i6 = SpatialInertia6x6::from_packed(&si);
        let result = i6.transform_to_parent(&SpatialTransform::IDENTITY);
        approx_eq_6x6(&result, &i6, 1e-10, "identity transform should be no-op");
    }

    /// For a point mass m at body-frame origin, translating by r should give
    /// from_com(m, r, 0).
    #[test]
    fn test_transform_to_parent_pure_translation() {
        let m = 3.0;
        let r = [0.5, 0.0, 0.0];
        // Body inertia: mass at origin (com=[0,0,0])
        let si_body = SpatialInertia::from_com(m, [0.0; 3], [[0.0; 3]; 3]);
        let i6 = SpatialInertia6x6::from_packed(&si_body);

        // Transform with translation r (child origin at r in parent frame)
        let x = SpatialTransform::from_translation(r);
        let result = i6.transform_to_parent(&x);

        // Expected: same mass at com=r in parent frame
        let si_expected = SpatialInertia::from_com(m, r, [[0.0; 3]; 3]);
        let expected_6x6 = SpatialInertia6x6::from_packed(&si_expected);

        approx_eq_6x6(&result, &expected_6x6, 1e-10, "pure translation transform");
    }

    /// Check that from_packed round-trips: from_packed(si).mul_vec6(v) == si.mul_vec(v)
    #[test]
    fn test_from_packed_roundtrip() {
        let si = SpatialInertia::from_com(
            2.0,
            [0.3, 0.1, 0.05],
            [[0.1, 0.0, 0.0], [0.0, 0.2, 0.0], [0.0, 0.0, 0.15]],
        );
        let i6 = SpatialInertia6x6::from_packed(&si);

        let v = SpatialVec::new([1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        let f1 = si.mul_vec(&v);
        let f2 = i6.mul_vec6(&v);

        let eps = 1e-10;
        assert!(
            (f1.angular[0] - f2.angular[0]).abs() < eps,
            "ang[0] mismatch: {} vs {}",
            f1.angular[0],
            f2.angular[0]
        );
        assert!(
            (f1.angular[1] - f2.angular[1]).abs() < eps,
            "ang[1] mismatch: {} vs {}",
            f1.angular[1],
            f2.angular[1]
        );
        assert!(
            (f1.angular[2] - f2.angular[2]).abs() < eps,
            "ang[2] mismatch: {} vs {}",
            f1.angular[2],
            f2.angular[2]
        );
        assert!(
            (f1.linear[0] - f2.linear[0]).abs() < eps,
            "lin[0] mismatch: {} vs {}",
            f1.linear[0],
            f2.linear[0]
        );
        assert!(
            (f1.linear[1] - f2.linear[1]).abs() < eps,
            "lin[1] mismatch: {} vs {}",
            f1.linear[1],
            f2.linear[1]
        );
        assert!(
            (f1.linear[2] - f2.linear[2]).abs() < eps,
            "lin[2] mismatch: {} vs {}",
            f1.linear[2],
            f2.linear[2]
        );
    }
}

#[cfg(test)]
mod transform_tests2 {
    use super::*;

    /// Verify that transform_to_parent with rotation+translation gives correct result.
    ///
    /// For a point mass m at body-frame origin (com=[0,0,0]) expressed in a child frame
    /// that is rotated and translated from parent, the transformed inertia should
    /// equal from_com(m, translated_com, 0) in parent frame.
    #[test]
    fn test_transform_to_parent_rot_trans() {
        let m = 1.0_f64;
        let r = [1.0_f64, 0.3_f64, 0.0_f64]; // Translation: child origin off-axis in parent
        let angle = -0.7_f64; // Rotation by -0.7 rad about z

        // Body inertia: point mass at com=[0.4,0.1,0] in child frame (non-axis-aligned)
        let com_child = [0.4_f64, 0.1_f64, 0.0_f64];
        let si_child = SpatialInertia::from_com(m, com_child, [[0.0_f64; 3]; 3]);
        let i6_child = SpatialInertia6x6::from_packed(&si_child);

        // Pose of the child in the parent: X_T (translate r) then X_J (rotate angle)
        let x_t = SpatialTransform::from_translation(r);
        let x_j = SpatialTransform::from_axis_angle([0.0, 0.0, 1.0], angle);
        let x_full = x_t.compose(&x_j);

        // Transform the inertia to parent frame
        let i6_parent = i6_child.transform_to_parent(&x_full);

        // `x_full` is the pose of the child frame in the parent frame, so a
        // child-frame point maps to the parent frame as `r + R·com_child`.
        let com_parent = add3(x_full.trans, mat3_vec(x_full.rot, com_child));

        // Expected inertia: mass at com_parent in parent frame (with zero I_com)
        let si_expected = SpatialInertia::from_com(m, com_parent, [[0.0_f64; 3]; 3]);
        let i6_expected = SpatialInertia6x6::from_packed(&si_expected);

        // Check
        let eps = 1e-8;
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (i6_parent.data[i][j] - i6_expected.data[i][j]).abs() < eps,
                    "data[{i}][{j}]: got {:.10}, expected {:.10}",
                    i6_parent.data[i][j],
                    i6_expected.data[i][j]
                );
            }
        }
    }
}

#[cfg(test)]
mod convention_tests {
    use super::*;

    fn sample_pose() -> SpatialTransform {
        SpatialTransform::from_translation([0.4, -0.7, 1.1])
            .compose(&SpatialTransform::from_axis_angle([0.3, -0.5, 0.8], 1.234))
    }

    fn i6_close(a: &SpatialInertia6x6, b: &SpatialInertia6x6, eps: f64) {
        for i in 0..6 {
            for j in 0..6 {
                assert!(
                    (a.data[i][j] - b.data[i][j]).abs() < eps,
                    "[{i}][{j}]: {:.12} vs {:.12}",
                    a.data[i][j],
                    b.data[i][j]
                );
            }
        }
    }

    /// Power is frame-invariant: `(X·v)·f_child == v·(Xᵀ·f_child)`.
    #[test]
    fn test_motion_and_force_transforms_are_dual() {
        let x = sample_pose();
        let v_parent = SpatialVec::new([0.3, -1.2, 0.5], [2.0, 0.1, -0.7]);
        let f_child = SpatialVec::new([-0.4, 0.9, 1.3], [0.6, -2.2, 0.8]);
        let lhs = x.apply_velocity(&v_parent).dot(&f_child);
        let rhs = v_parent.dot(&x.apply_force(&f_child));
        assert!((lhs - rhs).abs() < 1e-12, "power mismatch {lhs} vs {rhs}");
    }

    /// `apply_velocity` re-expresses a rigid motion: the child-frame linear
    /// part is the velocity of the child origin, in child coordinates.
    #[test]
    fn test_apply_velocity_matches_point_velocity() {
        let x = sample_pose();
        let omega = [0.2, -0.4, 0.9];
        let v_o = [1.0, 0.5, -0.3]; // velocity of the parent origin (parent coords)
        let v_child = x.apply_velocity(&SpatialVec::new(omega, v_o));
        // Child origin velocity in parent coords: v_o + ω × r.
        let v_point = add3(v_o, cross3(omega, x.trans));
        let rt = mat3_transpose(x.rot);
        let expected_lin = mat3_vec(rt, v_point);
        let expected_ang = mat3_vec(rt, omega);
        for k in 0..3 {
            assert!((v_child.linear[k] - expected_lin[k]).abs() < 1e-12);
            assert!((v_child.angular[k] - expected_ang[k]).abs() < 1e-12);
        }
    }

    /// `apply_velocity` composes like the poses it is built from.
    #[test]
    fn test_apply_velocity_respects_composition() {
        let a = sample_pose();
        let b = SpatialTransform::from_translation([-0.2, 0.3, 0.5])
            .compose(&SpatialTransform::from_axis_angle([1.0, 0.2, 0.0], -0.6));
        let v = SpatialVec::new([0.7, 0.1, -0.3], [0.2, -0.9, 1.4]);
        let chained = b.apply_velocity(&a.apply_velocity(&v));
        let direct = a.compose(&b).apply_velocity(&v);
        for k in 0..3 {
            assert!((chained.angular[k] - direct.angular[k]).abs() < 1e-12);
            assert!((chained.linear[k] - direct.linear[k]).abs() < 1e-12);
        }
    }

    /// Full (non-diagonal) inertia: `transform_to_parent` equals the explicit
    /// parent-frame inertia built from the rotated COM and `R·I_com·Rᵀ`, and
    /// `transform_inertia` (packed form) agrees with both.
    #[test]
    fn test_transform_to_parent_rotates_full_inertia() {
        let x = sample_pose();
        let m = 1.8;
        let com = [0.15, -0.25, 0.4];
        let i_com = [[0.21, 0.03, -0.05], [0.03, 0.17, 0.02], [-0.05, 0.02, 0.11]];
        let si = SpatialInertia::from_com(m, com, i_com);

        let got = SpatialInertia6x6::from_packed(&si).transform_to_parent(&x);

        let com_parent = add3(x.trans, mat3_vec(x.rot, com));
        let i_com_parent = mat3_mul(mat3_mul(x.rot, i_com), mat3_transpose(x.rot));
        let expected =
            SpatialInertia6x6::from_packed(&SpatialInertia::from_com(m, com_parent, i_com_parent));
        i6_close(&got, &expected, 1e-12);

        let packed = SpatialInertia6x6::from_packed(&x.transform_inertia(&si));
        i6_close(&packed, &expected, 1e-12);
    }

    /// `I·v` transforms consistently: `Xᵀ·(I_child·(X·v)) == I_parent·v`.
    #[test]
    fn test_inertia_transform_consistent_with_motion_and_force() {
        let x = sample_pose();
        let si = SpatialInertia::from_com(
            2.3,
            [0.3, 0.1, -0.2],
            [[0.3, 0.01, 0.0], [0.01, 0.2, -0.03], [0.0, -0.03, 0.25]],
        );
        let i_child = SpatialInertia6x6::from_packed(&si);
        let i_parent = i_child.transform_to_parent(&x);
        let v = SpatialVec::new([0.4, -0.8, 0.3], [-1.1, 0.6, 0.9]);
        let via_child = x.apply_force(&i_child.mul_vec6(&x.apply_velocity(&v)));
        let direct = i_parent.mul_vec6(&v);
        for k in 0..3 {
            assert!((via_child.angular[k] - direct.angular[k]).abs() < 1e-12);
            assert!((via_child.linear[k] - direct.linear[k]).abs() < 1e-12);
        }
    }
}
