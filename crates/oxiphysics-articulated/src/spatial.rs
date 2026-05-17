// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
#![allow(clippy::needless_range_loop)]

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
        for i in 0..3 {
            for j in 0..3 {
                io[i][j] += other.rot_inertia_origin[i][j];
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
/// Transforms spatial vectors from one coordinate frame to another using the
/// Plücker coordinate representation:
///
/// ```text
/// X * v_motion = [ R·ω ;  R·(v - r×ω) ]
/// X^T * f_force = [ R^T·(n - r×f) ;  R^T·f ]    (force transform)
/// ```
///
/// where `R` is the rotation matrix and `r` is the translation from the
/// parent frame origin to the child frame origin (expressed in parent frame).
#[derive(Debug, Clone, Copy)]
pub struct SpatialTransform {
    /// Rotation matrix from child to parent (row-major, R\[i\]\[j\]).
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

    /// Apply Plücker motion transform: `X * v_motion`.
    ///
    /// `v'_ang = R · ω`
    /// `v'_lin = R · (v_lin − r × ω)`
    pub fn apply_velocity(&self, v: &SpatialVec) -> SpatialVec {
        let r_omega = mat3_vec(self.rot, v.angular);
        let r_cross_omega = cross3(self.trans, v.angular);
        let v_minus = sub3(v.linear, r_cross_omega);
        let r_v_minus = mat3_vec(self.rot, v_minus);
        SpatialVec::new(r_omega, r_v_minus)
    }

    /// Apply force transform: `X^{−T} * f = X^* · f` (co-vector transform).
    ///
    /// For a force vector `f = [n; f_lin]`:
    ///
    /// `f'_ang = R · n + (R · r) × (R · f_lin)`
    /// `f'_lin = R · f_lin`
    ///
    /// Equivalently (same as transpose of inverse of motion transform):
    ///
    /// `n' = R · n + r × R · f_lin`  (in parent frame)
    /// `f' = R · f_lin`
    ///
    /// Wait — the correct Featherstone force transform for X^T is:
    ///
    /// `[n']   [  R    0 ] [n ]   [   R·n          ]`
    /// `[f'] = [ -r×R  R ] [f ] = [ -r × (R·n) + R·f ]`
    ///
    /// But `X^T` applied to a *covector* (force) from child to parent is:
    ///
    /// `n' = R^T · n − (R^T · r) × (R^T · f)`
    ///
    /// The direction depends on convention. Here we implement the transform that
    /// carries a force **from child frame to parent frame**, which matches the
    /// inward pass of RNEA:
    ///
    /// `n'_parent = R^T · n_child + r × R^T · f_child`
    /// `f'_parent = R^T · f_child`
    pub fn apply_force(&self, f: &SpatialVec) -> SpatialVec {
        let rt = mat3_transpose(self.rot);
        let f_parent = mat3_vec(rt, f.linear);
        let n_rotated = mat3_vec(rt, f.angular);
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

    /// Apply the inertia transform: `X * I * X^T` (transform inertia to new frame).
    ///
    /// Used to express body inertia in the parent frame.
    pub fn transform_inertia(&self, inertia: &SpatialInertia) -> SpatialInertia {
        // New mass is unchanged.
        // New com in parent frame: R · com_child + trans
        let new_com = add3(mat3_vec(self.rot, inertia.com), self.trans);
        // New inertia about parent origin (apply parallel-axis via frame change).
        // I'_origin_parent = R · I_origin_child · R^T + parallel_axis_correction_for_translation
        let rt = mat3_transpose(self.rot);
        let ri = mat3_mul(self.rot, inertia.rot_inertia_origin);
        let rirt = mat3_mul(ri, rt);
        // Now add parallel-axis for the translation from child_origin → parent_origin
        // I'_parent = R·I_child_origin·R^T + m*(|trans|^2 I - trans ⊗ trans)
        // where trans is the vector from parent origin to child origin (self.trans)
        let m = inertia.mass;
        let t = self.trans;
        let tt = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
        let mut io = rirt;
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                io[i][j] += m * (tt * delta - t[i] * t[j]);
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
        for i in 0..6 {
            for j in 0..6 {
                result[i] += self.data[i][j] * vv[j];
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
        for i in 0..6 {
            for j in 0..6 {
                d[i][j] = self.data[i][j] + other.data[i][j];
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

    /// Transform inertia from child frame to parent frame using the motion transform `X`.
    ///
    /// Given `X` as the **parent-to-child** motion transform (i.e., `v_child = X * v_parent`),
    /// this computes the child's inertia expressed in the parent frame:
    ///
    /// ```text
    /// I_parent = X^T * I_child * X
    /// ```
    ///
    /// where `X` is the 6×6 Plücker motion matrix:
    ///
    /// ```text
    /// X = [  R      0  ]
    ///     [ -R·r×   R  ]
    /// ```
    ///
    /// and `r` is the translation from parent origin to child origin in parent coordinates.
    pub fn transform_to_parent(&self, x: &SpatialTransform) -> Self {
        // Build 6×6 Plücker motion transform X (parent→child)
        // v_child = X * v_parent: [R·ω_p; R·(v_p - r×ω_p)]
        let r = x.rot;
        let rx = skew(x.trans); // r× = skew(r)

        // neg_r_rx = -R · r× (lower-left block of X)
        let mut neg_r_rx = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                neg_r_rx[i][j] = -(r[i][0] * rx[0][j] + r[i][1] * rx[1][j] + r[i][2] * rx[2][j]);
            }
        }

        // x6[i][j] = X[i][j]
        // [R  0 ]  → top-left (0..3, 0..3) = R, top-right (0..3, 3..6) = 0
        // [-Rrx R] → bottom-left (3..6, 0..3) = -R·r×, bottom-right (3..6, 3..6) = R
        let mut x6 = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                x6[i][j] = r[i][j];
                x6[i + 3][j] = neg_r_rx[i][j];
                x6[i + 3][j + 3] = r[i][j];
            }
        }

        // Compute I_parent = X^T * I_child * X
        // Step 1: tmp = I_child * X
        let mut tmp = [[0.0f64; 6]; 6];
        for i in 0..6 {
            for j in 0..6 {
                for k in 0..6 {
                    tmp[i][j] += self.data[i][k] * x6[k][j];
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

        // Build the parent-to-child transform: X_T(translate r) then X_J(rotate angle)
        let x_t = SpatialTransform::from_translation(r);
        let x_j = SpatialTransform::from_axis_angle([0.0, 0.0, 1.0], angle);
        let x_full = x_t.compose(&x_j);

        // Transform the inertia to parent frame
        let i6_parent = i6_child.transform_to_parent(&x_full);

        // Expected: com in parent frame
        // Child origin is at r=[1,0,0] in parent frame.
        // com_child=[0.4,0,0] in child frame, rotated by R_j2, translated by r.
        // com_parent = R_j2^T * com_child + r  -- wait, which direction?
        //
        // Actually: com_parent = r + R_j2^T * com_child
        //  No! The transform is parent-to-child, so:
        //  v_child = R * v_parent + ...
        //  To go from child frame to parent frame:
        //  v_parent = R^T * v_child (for rotation)
        //
        // A point at com_child in child frame is at com_parent = r + R^T * com_child
        // where R = x_full.rot (parent→child rotation).
        let r_mat = x_full.rot;
        let r_t = mat3_transpose(r_mat);
        let com_parent = add3(r, mat3_vec(r_t, com_child));

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
