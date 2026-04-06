// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Dual quaternions for rigid body transformations and screw motions.
//!
//! Dual quaternions combine a rotation quaternion (real part) with a
//! translation-encoding quaternion (dual part) to represent general rigid
//! body motions in a single algebraic object.  They are particularly well
//! suited for screw-motion interpolation (ScLERP).

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Quaternion
// ─────────────────────────────────────────────────────────────────────────────

/// A unit or general quaternion `w + xi + yj + zk`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Quaternion {
    /// Scalar (real) component.
    pub w: f64,
    /// `i` component.
    pub x: f64,
    /// `j` component.
    pub y: f64,
    /// `k` component.
    pub z: f64,
}

impl Quaternion {
    /// Creates a new quaternion from its four components.
    pub fn new(w: f64, x: f64, y: f64, z: f64) -> Self {
        Self { w, x, y, z }
    }

    /// Returns the multiplicative identity quaternion `(1, 0, 0, 0)`.
    pub fn identity() -> Self {
        Self::new(1.0, 0.0, 0.0, 0.0)
    }

    /// Constructs a unit quaternion representing a rotation of `angle` radians
    /// around `axis`.  `axis` need not be normalised; it is normalised
    /// internally.
    ///
    /// # Panics
    /// Panics (debug) if `axis` has zero length.
    pub fn from_axis_angle(axis: [f64; 3], angle: f64) -> Self {
        let len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let inv = 1.0 / len;
        let nx = axis[0] * inv;
        let ny = axis[1] * inv;
        let nz = axis[2] * inv;
        let half = angle * 0.5;
        let s = half.sin();
        Self::new(half.cos(), nx * s, ny * s, nz * s)
    }

    /// Hamilton product `self * other`.
    pub fn mul(&self, other: &Quaternion) -> Quaternion {
        Quaternion::new(
            self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
            self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
        )
    }

    /// Returns the conjugate `(w, -x, -y, -z)`.
    pub fn conjugate(&self) -> Quaternion {
        Quaternion::new(self.w, -self.x, -self.y, -self.z)
    }

    /// Returns the Euclidean norm.
    pub fn norm(&self) -> f64 {
        (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Returns a normalised copy.  If the norm is zero the result is
    /// `identity`.
    pub fn normalize(&self) -> Quaternion {
        let n = self.norm();
        if n < f64::EPSILON {
            return Quaternion::identity();
        }
        let inv = 1.0 / n;
        Quaternion::new(self.w * inv, self.x * inv, self.y * inv, self.z * inv)
    }

    /// Rotates a 3-D vector by the unit quaternion using `q v q*`.
    pub fn rotate_vector(&self, v: [f64; 3]) -> [f64; 3] {
        let qv = Quaternion::new(0.0, v[0], v[1], v[2]);
        let res = self.mul(&qv).mul(&self.conjugate());
        [res.x, res.y, res.z]
    }

    /// Euclidean dot product in 4-D.
    pub fn dot(&self, other: &Quaternion) -> f64 {
        self.w * other.w + self.x * other.x + self.y * other.y + self.z * other.z
    }

    /// Spherical linear interpolation between `self` and `other` at parameter
    /// `t ∈ [0, 1]`.  Both quaternions should be unit quaternions.
    pub fn slerp(&self, other: &Quaternion, t: f64) -> Quaternion {
        let mut d = self.dot(other);
        // Take the shorter arc.
        let other_adj = if d < 0.0 {
            d = -d;
            Quaternion::new(-other.w, -other.x, -other.y, -other.z)
        } else {
            *other
        };

        if d > 1.0 - 1e-10 {
            // Nearly identical — linear blend and renormalize.
            let w = self.w + t * (other_adj.w - self.w);
            let x = self.x + t * (other_adj.x - self.x);
            let y = self.y + t * (other_adj.y - self.y);
            let z = self.z + t * (other_adj.z - self.z);
            return Quaternion::new(w, x, y, z).normalize();
        }

        let theta = d.acos();
        let sin_theta = theta.sin();
        let s0 = ((1.0 - t) * theta).sin() / sin_theta;
        let s1 = (t * theta).sin() / sin_theta;
        Quaternion::new(
            s0 * self.w + s1 * other_adj.w,
            s0 * self.x + s1 * other_adj.x,
            s0 * self.y + s1 * other_adj.y,
            s0 * self.z + s1 * other_adj.z,
        )
    }

    /// Converts to a 3×3 rotation matrix.
    pub fn to_rotation_matrix(&self) -> [[f64; 3]; 3] {
        let q = self.normalize();
        let (w, x, y, z) = (q.w, q.x, q.y, q.z);
        [
            [
                1.0 - 2.0 * (y * y + z * z),
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                1.0 - 2.0 * (x * x + z * z),
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                1.0 - 2.0 * (x * x + y * y),
            ],
        ]
    }

    /// Constructs a quaternion from a 3×3 rotation matrix using Shepperd's
    /// method.
    pub fn from_rotation_matrix(m: &[[f64; 3]; 3]) -> Quaternion {
        let trace = m[0][0] + m[1][1] + m[2][2];
        if trace > 0.0 {
            let s = 0.5 / (trace + 1.0).sqrt();
            Quaternion::new(
                0.25 / s,
                (m[2][1] - m[1][2]) * s,
                (m[0][2] - m[2][0]) * s,
                (m[1][0] - m[0][1]) * s,
            )
        } else if m[0][0] > m[1][1] && m[0][0] > m[2][2] {
            let s = 2.0 * (1.0 + m[0][0] - m[1][1] - m[2][2]).sqrt();
            Quaternion::new(
                (m[2][1] - m[1][2]) / s,
                0.25 * s,
                (m[0][1] + m[1][0]) / s,
                (m[0][2] + m[2][0]) / s,
            )
        } else if m[1][1] > m[2][2] {
            let s = 2.0 * (1.0 + m[1][1] - m[0][0] - m[2][2]).sqrt();
            Quaternion::new(
                (m[0][2] - m[2][0]) / s,
                (m[0][1] + m[1][0]) / s,
                0.25 * s,
                (m[1][2] + m[2][1]) / s,
            )
        } else {
            let s = 2.0 * (1.0 + m[2][2] - m[0][0] - m[1][1]).sqrt();
            Quaternion::new(
                (m[1][0] - m[0][1]) / s,
                (m[0][2] + m[2][0]) / s,
                (m[1][2] + m[2][1]) / s,
                0.25 * s,
            )
        }
    }

    /// Quaternion exponential.  For a pure quaternion `(0, v)` this produces
    /// `cos(|v|) + sin(|v|)/|v| * v`.
    pub fn exp(&self) -> Quaternion {
        let v_norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        let ew = self.w.exp();
        if v_norm < f64::EPSILON {
            return Quaternion::new(ew, 0.0, 0.0, 0.0);
        }
        let s = ew * v_norm.sin() / v_norm;
        Quaternion::new(ew * v_norm.cos(), self.x * s, self.y * s, self.z * s)
    }

    /// Quaternion natural logarithm.  For a unit quaternion this produces a
    /// pure quaternion encoding the rotation axis and half-angle.
    pub fn ln(&self) -> Quaternion {
        let n = self.norm();
        let v_norm = (self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if v_norm < f64::EPSILON {
            return Quaternion::new(n.ln(), 0.0, 0.0, 0.0);
        }
        let theta = (self.w / n).acos();
        let s = theta / v_norm;
        Quaternion::new(n.ln(), self.x * s, self.y * s, self.z * s)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DualNumber
// ─────────────────────────────────────────────────────────────────────────────

/// A dual number `a + ε b` where `ε² = 0`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualNumber {
    /// Real part.
    pub real: f64,
    /// Dual (infinitesimal) part.
    pub dual: f64,
}

impl DualNumber {
    /// Creates a new dual number.
    pub fn new(real: f64, dual: f64) -> Self {
        Self { real, dual }
    }

    /// Component-wise addition.
    pub fn add(&self, other: &DualNumber) -> DualNumber {
        DualNumber::new(self.real + other.real, self.dual + other.dual)
    }

    /// Dual-number multiplication: `(a + εb)(c + εd) = ac + ε(ad + bc)`.
    pub fn mul(&self, other: &DualNumber) -> DualNumber {
        DualNumber::new(
            self.real * other.real,
            self.real * other.dual + self.dual * other.real,
        )
    }

    /// Dual conjugate `(a, -b)`.
    pub fn conjugate(&self) -> DualNumber {
        DualNumber::new(self.real, -self.dual)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DualQuaternion
// ─────────────────────────────────────────────────────────────────────────────

/// A dual quaternion `q_r + ε q_d` encoding a rigid body transformation.
///
/// The *real* part `q_r` is a unit quaternion encoding the rotation.  The
/// *dual* part `q_d = ½ t̂ q_r` encodes the translation, where `t̂` is the
/// pure quaternion formed from the translation vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualQuaternion {
    /// Rotation part (unit quaternion).
    pub real: Quaternion,
    /// Translation-encoding part.
    pub dual: Quaternion,
}

impl DualQuaternion {
    /// Returns the identity dual quaternion (no rotation, no translation).
    pub fn identity() -> Self {
        Self {
            real: Quaternion::identity(),
            dual: Quaternion::new(0.0, 0.0, 0.0, 0.0),
        }
    }

    /// Constructs a dual quaternion from a unit rotation quaternion `q` and a
    /// translation vector `t`.
    ///
    /// The dual part is `½ * (0, t) * q`.
    pub fn from_rotation_translation(q: Quaternion, t: [f64; 3]) -> Self {
        let t_quat = Quaternion::new(0.0, t[0], t[1], t[2]);
        let dual = t_quat.mul(&q).scale(0.5);
        Self { real: q, dual }
    }

    /// Dual-quaternion product `self * other`.
    pub fn mul(&self, other: &DualQuaternion) -> DualQuaternion {
        DualQuaternion {
            real: self.real.mul(&other.real),
            dual: self.real.mul(&other.dual).add(&self.dual.mul(&other.real)),
        }
    }

    /// Dual-quaternion conjugate: conjugate both the real and dual parts.
    pub fn conjugate(&self) -> DualQuaternion {
        DualQuaternion {
            real: self.real.conjugate(),
            dual: self.dual.conjugate(),
        }
    }

    /// Normalises so that the real part has unit norm.
    pub fn normalize(&self) -> DualQuaternion {
        let n = self.real.norm();
        if n < f64::EPSILON {
            return DualQuaternion::identity();
        }
        let inv = 1.0 / n;
        DualQuaternion {
            real: self.real.scale(inv),
            dual: self.dual.scale(inv),
        }
    }

    /// Returns the norm (magnitude) of the real part.
    pub fn norm(&self) -> f64 {
        self.real.norm()
    }

    /// Transforms a 3-D point using the rigid body motion encoded in `self`.
    ///
    /// Equivalent to `dq * (1 + ε p̂) * dq*` in the dual-quaternion algebra.
    pub fn transform_point(&self, p: [f64; 3]) -> [f64; 3] {
        let dq = self.normalize();
        // Rotation
        let rotated = dq.real.rotate_vector(p);
        // Translation: t = 2 * dual * real_conjugate  (imaginary parts)
        let t = dq.to_translation();
        [rotated[0] + t[0], rotated[1] + t[1], rotated[2] + t[2]]
    }

    /// Returns the rotation as a normalised quaternion.
    pub fn to_rotation(&self) -> Quaternion {
        self.real.normalize()
    }

    /// Recovers the translation vector from the dual quaternion.
    ///
    /// `t = 2 * q_d * q_r*` (imaginary parts only).
    pub fn to_translation(&self) -> [f64; 3] {
        let dq = self.normalize();
        let t = dq.dual.mul(&dq.real.conjugate()).scale(2.0);
        [t.x, t.y, t.z]
    }

    /// Screw linear interpolation (ScLERP) between `self` and `other` at `t`.
    ///
    /// Uses the formula `(other * self⁻¹)^t * self`.
    pub fn sclerp(&self, other: &DualQuaternion, t: f64) -> DualQuaternion {
        let self_inv = self.conjugate();
        let diff = self_inv.mul(other).normalize();
        self.mul(&diff.pow(t))
    }

    /// Constructs a dual quaternion from a general screw motion.
    ///
    /// # Parameters
    /// * `axis`  – unit direction of the screw axis.
    /// * `angle` – rotation angle in radians.
    /// * `pitch` – translation per radian.
    /// * `point` – a point on the screw axis.
    pub fn from_twist(axis: [f64; 3], angle: f64, pitch: f64, point: [f64; 3]) -> DualQuaternion {
        let translation = pitch * angle;
        // Pure rotation about `axis` through origin.
        let q_rot = Quaternion::from_axis_angle(axis, angle);
        // Translation along axis plus moment contribution.
        let d = translation;
        let t = [
            axis[0] * d
                + (point[1] * axis[2] - point[2] * axis[1])
                    * 2.0
                    * (angle * 0.5).sin()
                    * (angle * 0.5).sin()
                    * 0.0, // moment term placeholder — see full formula below
            axis[1] * d,
            axis[2] * d,
        ];
        // Full screw: translation along axis = d * axis_dir,
        // plus a rotation of `point` around the axis.
        let q_t = Quaternion::from_axis_angle(axis, angle);
        // Translate to origin, rotate, translate back.
        let t_to_origin = [
            -point[0] + point[0],
            -point[1] + point[1],
            -point[2] + point[2],
        ];
        // The proper way: DQ(rotation) composed with DQ(translation along axis).
        let t_vec = [axis[0] * d, axis[1] * d, axis[2] * d];
        // Apply point-on-axis shift.
        let moment_rotated = q_rot.rotate_vector(point);
        let full_t = [
            t_vec[0] + (point[0] - moment_rotated[0]),
            t_vec[1] + (point[1] - moment_rotated[1]),
            t_vec[2] + (point[2] - moment_rotated[2]),
        ];
        let _ = t;
        let _ = t_to_origin;
        let _ = q_t;
        DualQuaternion::from_rotation_translation(q_rot, full_t)
    }

    // ── Internal helpers ──────────────────────────────────────────────────────

    /// Raises a dual quaternion to the real power `t` via the ScLERP formula.
    fn pow(&self, t: f64) -> DualQuaternion {
        // Decompose into screw parameters and re-compose.
        let sm = ScrewMotion::from_dual_quaternion(self);
        DualQuaternion::from_twist(sm.axis_dir, sm.angle * t, sm.pitch, sm.axis_point)
    }
}

// ── Private arithmetic helpers on Quaternion ─────────────────────────────────

impl Quaternion {
    fn scale(&self, s: f64) -> Quaternion {
        Quaternion::new(self.w * s, self.x * s, self.y * s, self.z * s)
    }

    fn add(&self, other: &Quaternion) -> Quaternion {
        Quaternion::new(
            self.w + other.w,
            self.x + other.x,
            self.y + other.y,
            self.z + other.z,
        )
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ScrewMotion
// ─────────────────────────────────────────────────────────────────────────────

/// A rigid body screw motion described by a Chasles-style parameterisation.
///
/// Every rigid body motion (up to pure translation) can be represented as a
/// rotation by `angle` radians around an axis, combined with a translation
/// of `pitch * angle` metres along that same axis.
#[derive(Debug, Clone, Copy)]
pub struct ScrewMotion {
    /// Unit direction vector of the screw axis.
    pub axis_dir: [f64; 3],
    /// A point on the screw axis.
    pub axis_point: [f64; 3],
    /// Rotation angle in radians.
    pub angle: f64,
    /// Translation per radian along the screw axis.
    pub pitch: f64,
}

impl ScrewMotion {
    /// Converts this screw motion to its dual-quaternion representation.
    pub fn to_dual_quaternion(&self) -> DualQuaternion {
        DualQuaternion::from_twist(self.axis_dir, self.angle, self.pitch, self.axis_point)
    }

    /// Extracts the screw-motion parameters from a (normalised) dual
    /// quaternion.
    pub fn from_dual_quaternion(dq: &DualQuaternion) -> ScrewMotion {
        let dq = dq.normalize();
        // Real part → rotation axis & angle.
        let q = dq.real;
        // Clamp for numerical safety.
        let w_clamped = q.w.clamp(-1.0, 1.0);
        let angle = 2.0 * w_clamped.acos();
        let sin_half = (1.0 - w_clamped * w_clamped).sqrt();

        let (axis_dir, axis_point, pitch) = if sin_half < 1e-10 {
            // Pure translation (no rotation).
            let t = dq.to_translation();
            let t_norm = (t[0] * t[0] + t[1] * t[1] + t[2] * t[2]).sqrt();
            let dir = if t_norm < f64::EPSILON {
                [0.0, 0.0, 1.0]
            } else {
                [t[0] / t_norm, t[1] / t_norm, t[2] / t_norm]
            };
            ([dir[0], dir[1], dir[2]], [0.0_f64, 0.0, 0.0], t_norm)
        } else {
            let inv = 1.0 / sin_half;
            let axis = [q.x * inv, q.y * inv, q.z * inv];

            // Dual part → moment and pitch.
            // d = 2 * q_d * q_r* (imaginary parts = translation).
            let t = dq.to_translation();
            // Pitch = translation along axis.
            let pitch_val = t[0] * axis[0] + t[1] * axis[1] + t[2] * axis[2];
            // Moment = axis × translation / |axis|² (moment of the screw).
            // axis_point = axis × (t - pitch*axis) / angle (simplified).
            let t_perp = [
                t[0] - pitch_val * axis[0],
                t[1] - pitch_val * axis[1],
                t[2] - pitch_val * axis[2],
            ];
            // axis_point from (I-R)*p = t_perp:
            // p = [(1-cos)*t_perp + sin*cross] / (2*(1-cos))
            // where cross = axis × t_perp, 1-cos = 2*sin²(θ/2)
            let cross = [
                axis[1] * t_perp[2] - axis[2] * t_perp[1],
                axis[2] * t_perp[0] - axis[0] * t_perp[2],
                axis[0] * t_perp[1] - axis[1] * t_perp[0],
            ];
            let one_minus_cos = 2.0 * sin_half * sin_half; // 1 - cos(θ)
            let sin_theta = 2.0 * sin_half * w_clamped; // sin(θ)
            let denom = 2.0 * one_minus_cos;
            let p = if denom.abs() < 1e-10 {
                [0.0, 0.0, 0.0]
            } else {
                [
                    (one_minus_cos * t_perp[0] + sin_theta * cross[0]) / denom,
                    (one_minus_cos * t_perp[1] + sin_theta * cross[1]) / denom,
                    (one_minus_cos * t_perp[2] + sin_theta * cross[2]) / denom,
                ]
            };
            let pitch_per_rad = if angle.abs() < 1e-10 {
                0.0
            } else {
                pitch_val / angle
            };
            (axis, p, pitch_per_rad)
        };

        ScrewMotion {
            axis_dir,
            axis_point,
            angle,
            pitch,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DLB — Dual Quaternion Linear Blending
// ─────────────────────────────────────────────────────────────────────────────

/// Dual Quaternion Linear Blending (DLB) — blend N dual quaternions with given weights.
///
/// The result is normalised so the real part has unit norm.
/// This is the standard blending used in dual quaternion skinning (GPU-friendly).
///
/// If all weights sum to zero the identity dual quaternion is returned.
#[allow(dead_code)]
pub fn dual_quaternion_blend(dqs: &[DualQuaternion], weights: &[f64]) -> DualQuaternion {
    let n = dqs.len().min(weights.len());
    if n == 0 {
        return DualQuaternion::identity();
    }
    let mut sum_real = Quaternion::new(0.0, 0.0, 0.0, 0.0);
    let mut sum_dual = Quaternion::new(0.0, 0.0, 0.0, 0.0);
    let first_real = dqs[0].real;
    for i in 0..n {
        // Ensure the real part is on the same hemisphere as the first quaternion
        let dq = if dqs[i].real.dot(&first_real) < 0.0 {
            DualQuaternion {
                real: dqs[i].real.scale(-1.0),
                dual: dqs[i].dual.scale(-1.0),
            }
        } else {
            dqs[i]
        };
        let w = weights[i];
        sum_real = sum_real.add(&dq.real.scale(w));
        sum_dual = sum_dual.add(&dq.dual.scale(w));
    }
    DualQuaternion {
        real: sum_real,
        dual: sum_dual,
    }
    .normalize()
}

// ─────────────────────────────────────────────────────────────────────────────
// Dual quaternion interpolation helpers
// ─────────────────────────────────────────────────────────────────────────────

impl DualQuaternion {
    /// Linear interpolation (nLERP) between two dual quaternions.
    ///
    /// Cheaper than ScLERP; less accurate for large rotations.
    /// Both dual quaternions are normalised before blending.
    #[allow(dead_code)]
    pub fn nlerp(&self, other: &DualQuaternion, t: f64) -> DualQuaternion {
        let a = self.normalize();
        let b = other.normalize();
        // Ensure shortest path on real part
        let b = if a.real.dot(&b.real) < 0.0 {
            DualQuaternion {
                real: b.real.scale(-1.0),
                dual: b.dual.scale(-1.0),
            }
        } else {
            b
        };
        let real = a.real.scale(1.0 - t).add(&b.real.scale(t));
        let dual = a.dual.scale(1.0 - t).add(&b.dual.scale(t));
        DualQuaternion { real, dual }.normalize()
    }

    /// Rigid transform composition: compose `self` with `other`.
    ///
    /// The resulting transformation first applies `self`, then `other`.
    /// Equivalent to `other.mul(self)`.
    #[allow(dead_code)]
    pub fn compose_transforms(&self, other: &DualQuaternion) -> DualQuaternion {
        other.mul(self)
    }

    /// Invert a unit dual quaternion (conjugate of the real and dual parts).
    #[allow(dead_code)]
    pub fn invert(&self) -> DualQuaternion {
        self.conjugate()
    }

    /// Apply this rigid transform to a list of points.
    #[allow(dead_code)]
    pub fn transform_points(&self, points: &[[f64; 3]]) -> Vec<[f64; 3]> {
        points.iter().map(|&p| self.transform_point(p)).collect()
    }

    /// Extract the rotation angle (in radians) from the real quaternion.
    #[allow(dead_code)]
    pub fn rotation_angle(&self) -> f64 {
        let q = self.real.normalize();
        2.0 * q.w.clamp(-1.0, 1.0).acos()
    }

    /// Extract the rotation axis (unit vector) from the real quaternion.
    ///
    /// Returns \[0, 0, 1\] for the identity rotation.
    #[allow(dead_code)]
    pub fn rotation_axis(&self) -> [f64; 3] {
        let q = self.real.normalize();
        let s = (1.0 - q.w * q.w).sqrt();
        if s < 1e-12 {
            [0.0, 0.0, 1.0]
        } else {
            [q.x / s, q.y / s, q.z / s]
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Screw interpolation (ScLERP) utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Interpolate a sequence of `n` dual quaternion poses evenly spaced between
/// `start` and `end` (including both endpoints).
///
/// Uses ScLERP for each intermediate pose.
#[allow(dead_code)]
pub fn sclerp_sequence(
    start: &DualQuaternion,
    end: &DualQuaternion,
    n: usize,
) -> Vec<DualQuaternion> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![*start];
    }
    (0..n)
        .map(|i| {
            let t = i as f64 / (n - 1) as f64;
            start.sclerp(end, t)
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Plücker Coordinates
// ─────────────────────────────────────────────────────────────────────────────

/// A line in 3D represented by Plücker coordinates (direction, moment).
///
/// The direction vector `l` points along the line, and the moment `m = p × l`
/// where `p` is any point on the line.  Together `(l, m)` uniquely describe
/// a 3D line up to scale.
#[derive(Debug, Clone, Copy)]
pub struct PluckerLine {
    /// Unit direction vector of the line.
    pub direction: [f64; 3],
    /// Moment vector `m = p × direction` for any point `p` on the line.
    pub moment: [f64; 3],
}

impl PluckerLine {
    /// Construct a Plücker line passing through two distinct points.
    pub fn from_two_points(p: [f64; 3], q: [f64; 3]) -> Self {
        let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let direction = if len < f64::EPSILON {
            [0.0, 0.0, 1.0]
        } else {
            [d[0] / len, d[1] / len, d[2] / len]
        };
        // Moment = p × direction
        let moment = [
            p[1] * direction[2] - p[2] * direction[1],
            p[2] * direction[0] - p[0] * direction[2],
            p[0] * direction[1] - p[1] * direction[0],
        ];
        Self { direction, moment }
    }

    /// Construct from a point on the line and a direction vector.
    pub fn from_point_direction(p: [f64; 3], d: [f64; 3]) -> Self {
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let direction = if len < f64::EPSILON {
            [0.0, 0.0, 1.0]
        } else {
            [d[0] / len, d[1] / len, d[2] / len]
        };
        let moment = [
            p[1] * direction[2] - p[2] * direction[1],
            p[2] * direction[0] - p[0] * direction[2],
            p[0] * direction[1] - p[1] * direction[0],
        ];
        Self { direction, moment }
    }

    /// Returns the shortest distance between two skew lines (or 0 if they intersect).
    pub fn distance_to_line(&self, other: &PluckerLine) -> f64 {
        // For two lines l1=(d1,m1), l2=(d2,m2):
        // distance = |d1·m2 + d2·m1| / |d1 × d2|
        let cross_d = [
            self.direction[1] * other.direction[2] - self.direction[2] * other.direction[1],
            self.direction[2] * other.direction[0] - self.direction[0] * other.direction[2],
            self.direction[0] * other.direction[1] - self.direction[1] * other.direction[0],
        ];
        let cross_len =
            (cross_d[0] * cross_d[0] + cross_d[1] * cross_d[1] + cross_d[2] * cross_d[2]).sqrt();

        let dot1 = self.direction[0] * other.moment[0]
            + self.direction[1] * other.moment[1]
            + self.direction[2] * other.moment[2];
        let dot2 = other.direction[0] * self.moment[0]
            + other.direction[1] * self.moment[1]
            + other.direction[2] * self.moment[2];

        if cross_len < 1e-12 {
            // Parallel lines: pick points on each line (d × m for unit direction)
            // and compute |(p1 - p2) × d| / |d|
            let p1 = [
                self.direction[1] * self.moment[2] - self.direction[2] * self.moment[1],
                self.direction[2] * self.moment[0] - self.direction[0] * self.moment[2],
                self.direction[0] * self.moment[1] - self.direction[1] * self.moment[0],
            ];
            let p2 = [
                other.direction[1] * other.moment[2] - other.direction[2] * other.moment[1],
                other.direction[2] * other.moment[0] - other.direction[0] * other.moment[2],
                other.direction[0] * other.moment[1] - other.direction[1] * other.moment[0],
            ];
            let dp = [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]];
            // cross dp × self.direction
            let cx = dp[1] * self.direction[2] - dp[2] * self.direction[1];
            let cy = dp[2] * self.direction[0] - dp[0] * self.direction[2];
            let cz = dp[0] * self.direction[1] - dp[1] * self.direction[0];
            (cx * cx + cy * cy + cz * cz).sqrt()
        } else {
            (dot1 + dot2).abs() / cross_len
        }
    }

    /// Returns the foot of the perpendicular from a point to the line.
    pub fn closest_point_to(&self, p: [f64; 3]) -> [f64; 3] {
        // Project p onto the line: p - ((p - q) · d_hat) * d_hat + q_on_line
        // Reconstruct a point on the line: q = d × m / |d|² (for unit d, q = d × m)
        let q_on_line = [
            self.direction[1] * self.moment[2] - self.direction[2] * self.moment[1],
            self.direction[2] * self.moment[0] - self.direction[0] * self.moment[2],
            self.direction[0] * self.moment[1] - self.direction[1] * self.moment[0],
        ];
        let ap = [
            p[0] - q_on_line[0],
            p[1] - q_on_line[1],
            p[2] - q_on_line[2],
        ];
        let t = ap[0] * self.direction[0] + ap[1] * self.direction[1] + ap[2] * self.direction[2];
        [
            q_on_line[0] + t * self.direction[0],
            q_on_line[1] + t * self.direction[1],
            q_on_line[2] + t * self.direction[2],
        ]
    }

    /// Converts this Plücker line to a dual quaternion representation.
    ///
    /// The DQ encodes the line as a pure rotation at the moment direction.
    pub fn to_dual_quaternion_rotation(&self, angle: f64) -> DualQuaternion {
        // Rotate around the direction axis by `angle`, with the axis passing through
        // the point closest to the origin (given by d × m for unit d)
        let axis_point = [
            self.direction[1] * self.moment[2] - self.direction[2] * self.moment[1],
            self.direction[2] * self.moment[0] - self.direction[0] * self.moment[2],
            self.direction[0] * self.moment[1] - self.direction[1] * self.moment[0],
        ];
        DualQuaternion::from_twist(self.direction, angle, 0.0, axis_point)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Velocity Screw
// ─────────────────────────────────────────────────────────────────────────────

/// A velocity screw in body or spatial form, combining angular and linear velocity.
///
/// In screw theory, any instantaneous rigid body motion can be described by a
/// screw velocity: an angular velocity `ω` and a linear velocity `v` of the origin.
/// The velocity of a point `p` is `v + ω × p`.
#[derive(Debug, Clone, Copy)]
pub struct VelocityScrew {
    /// Angular velocity vector (ω).
    pub angular: [f64; 3],
    /// Linear velocity of the reference point (v).
    pub linear: [f64; 3],
}

impl VelocityScrew {
    /// Returns the zero velocity screw.
    pub fn zero() -> Self {
        Self {
            angular: [0.0; 3],
            linear: [0.0; 3],
        }
    }

    /// Computes the velocity of a point `p` attached to the body.
    ///
    /// `vel(p) = v + ω × p`
    pub fn velocity_at_point(&self, p: [f64; 3]) -> [f64; 3] {
        let w = self.angular;
        let v = self.linear;
        // ω × p
        let cross = [
            w[1] * p[2] - w[2] * p[1],
            w[2] * p[0] - w[0] * p[2],
            w[0] * p[1] - w[1] * p[0],
        ];
        [v[0] + cross[0], v[1] + cross[1], v[2] + cross[2]]
    }

    /// Returns the screw pitch: v · ω / |ω|² (translation per radian along axis).
    ///
    /// Returns 0 if angular velocity is zero.
    pub fn pitch(&self) -> f64 {
        let w_sq = self.angular[0] * self.angular[0]
            + self.angular[1] * self.angular[1]
            + self.angular[2] * self.angular[2];
        if w_sq < f64::EPSILON {
            return 0.0;
        }
        let v_dot_w = self.linear[0] * self.angular[0]
            + self.linear[1] * self.angular[1]
            + self.linear[2] * self.angular[2];
        v_dot_w / w_sq
    }

    /// Returns the screw axis direction (unit angular velocity), or z-axis if zero.
    pub fn axis_direction(&self) -> [f64; 3] {
        let len = (self.angular[0] * self.angular[0]
            + self.angular[1] * self.angular[1]
            + self.angular[2] * self.angular[2])
            .sqrt();
        if len < f64::EPSILON {
            [0.0, 0.0, 1.0]
        } else {
            [
                self.angular[0] / len,
                self.angular[1] / len,
                self.angular[2] / len,
            ]
        }
    }

    /// Adds two velocity screws component-wise.
    pub fn add(&self, other: &VelocityScrew) -> VelocityScrew {
        VelocityScrew {
            angular: [
                self.angular[0] + other.angular[0],
                self.angular[1] + other.angular[1],
                self.angular[2] + other.angular[2],
            ],
            linear: [
                self.linear[0] + other.linear[0],
                self.linear[1] + other.linear[1],
                self.linear[2] + other.linear[2],
            ],
        }
    }

    /// Scales the velocity screw by a scalar.
    pub fn scale(&self, s: f64) -> VelocityScrew {
        VelocityScrew {
            angular: [
                self.angular[0] * s,
                self.angular[1] * s,
                self.angular[2] * s,
            ],
            linear: [self.linear[0] * s, self.linear[1] * s, self.linear[2] * s],
        }
    }

    /// Converts this instantaneous velocity screw to a dual quaternion
    /// for the given time step `dt` (first-order integration).
    pub fn to_dual_quaternion_step(&self, dt: f64) -> DualQuaternion {
        let angle = {
            let w = self.angular;
            (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt() * dt
        };
        let axis = self.axis_direction();
        let pitch = self.pitch();
        DualQuaternion::from_twist(axis, angle, pitch, [0.0; 3])
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Denavit-Hartenberg Transform
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the Denavit-Hartenberg (DH) rigid body transform as a dual quaternion.
///
/// The standard DH convention applies four elementary transformations:
/// 1. Rotate `theta` about z-axis
/// 2. Translate `d` along z-axis
/// 3. Translate `a` along x-axis
/// 4. Rotate `alpha` about x-axis
///
/// # Arguments
/// * `a`     – link length (translation along x)
/// * `alpha` – link twist (rotation about x)
/// * `d`     – link offset (translation along z)
/// * `theta` – joint angle (rotation about z)
///
/// # Returns
/// A dual quaternion representing the combined rigid transform.
#[allow(clippy::too_many_arguments)]
pub fn dh_transform(a: f64, alpha: f64, d: f64, theta: f64) -> DualQuaternion {
    // Step 1+2: rotate theta about z then translate d along z
    let q_theta = Quaternion::from_axis_angle([0.0, 0.0, 1.0], theta);
    let dq_theta = DualQuaternion::from_rotation_translation(q_theta, [0.0, 0.0, d]);

    // Step 3+4: translate a along x then rotate alpha about x
    let q_alpha = Quaternion::from_axis_angle([1.0, 0.0, 0.0], alpha);
    let dq_alpha = DualQuaternion::from_rotation_translation(q_alpha, [a, 0.0, 0.0]);

    dq_theta.mul(&dq_alpha)
}

/// Builds a full DH kinematic chain from a list of DH parameter sets.
///
/// Each entry is `(a, alpha, d, theta)`.  The transforms are composed in order:
/// `T_total = T_0 * T_1 * ... * T_{n-1}`.
pub fn dh_chain(params: &[(f64, f64, f64, f64)]) -> DualQuaternion {
    params
        .iter()
        .fold(DualQuaternion::identity(), |acc, &(a, al, d, th)| {
            acc.mul(&dh_transform(a, al, d, th))
        })
}

// ─────────────────────────────────────────────────────────────────────────────
// Rigid body pose blending
// ─────────────────────────────────────────────────────────────────────────────

/// Blend N rigid body poses (dual quaternions) with given weights using DLB.
///
/// This is a wrapper around `dual_quaternion_blend` with a more descriptive name
/// for the rigid body pose blending use case (e.g. character skinning, IK blending).
pub fn rigid_body_pose_blend(poses: &[DualQuaternion], weights: &[f64]) -> DualQuaternion {
    dual_quaternion_blend(poses, weights)
}

/// Geodesic interpolation between two poses using ScLERP.
///
/// Equivalent to `start.sclerp(end, t)` but exposed as a free function.
pub fn geodesic_interpolate(
    start: &DualQuaternion,
    end: &DualQuaternion,
    t: f64,
) -> DualQuaternion {
    start.sclerp(end, t)
}

/// Compute the relative transform from `from` to `to`: `from⁻¹ * to`.
pub fn relative_transform(from: &DualQuaternion, to: &DualQuaternion) -> DualQuaternion {
    from.invert().mul(to)
}

/// Extract angular velocity from two adjacent dual quaternions separated by time `dt`.
///
/// Uses finite differences on the rotation part: `ω ≈ 2 * ln(q1* q2).xyz / dt`.
pub fn finite_difference_angular_velocity(
    dq1: &DualQuaternion,
    dq2: &DualQuaternion,
    dt: f64,
) -> [f64; 3] {
    if dt.abs() < f64::EPSILON {
        return [0.0; 3];
    }
    let q1 = dq1.to_rotation();
    let q2 = dq2.to_rotation();
    let dq = q1.conjugate().mul(&q2);
    let log = dq.ln();
    // log of unit quaternion: (0, θ/2 * axis_xyz)
    let scale = 2.0 / dt;
    [log.x * scale, log.y * scale, log.z * scale]
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    fn approx_eq(a: f64, b: f64, eps: f64) -> bool {
        (a - b).abs() < eps
    }

    fn vec3_approx_eq(a: [f64; 3], b: [f64; 3], eps: f64) -> bool {
        approx_eq(a[0], b[0], eps) && approx_eq(a[1], b[1], eps) && approx_eq(a[2], b[2], eps)
    }

    fn quat_approx_eq(a: Quaternion, b: Quaternion, eps: f64) -> bool {
        // Account for double-cover: q and -q represent the same rotation.
        let direct = approx_eq(a.w, b.w, eps)
            && approx_eq(a.x, b.x, eps)
            && approx_eq(a.y, b.y, eps)
            && approx_eq(a.z, b.z, eps);
        let negated = approx_eq(a.w, -b.w, eps)
            && approx_eq(a.x, -b.x, eps)
            && approx_eq(a.y, -b.y, eps)
            && approx_eq(a.z, -b.z, eps);
        direct || negated
    }

    fn dq_approx_eq(a: DualQuaternion, b: DualQuaternion, eps: f64) -> bool {
        quat_approx_eq(a.real, b.real, eps) && quat_approx_eq(a.dual, b.dual, eps)
    }

    // ── Quaternion tests ──────────────────────────────────────────────────────

    #[test]
    fn test_rotate_z90() {
        let q = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let v = q.rotate_vector([1.0, 0.0, 0.0]);
        assert!(vec3_approx_eq(v, [0.0, 1.0, 0.0], 1e-12));
    }

    #[test]
    fn test_slerp_endpoints() {
        let q1 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], 0.0);
        let q2 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let s0 = q1.slerp(&q2, 0.0);
        let s1 = q1.slerp(&q2, 1.0);
        assert!(quat_approx_eq(s0, q1, 1e-10));
        assert!(quat_approx_eq(s1, q2, 1e-10));
    }

    #[test]
    fn test_to_rotation_matrix_orthogonal() {
        let q = Quaternion::from_axis_angle([1.0, 1.0, 0.0], PI / 3.0);
        let m = q.to_rotation_matrix();
        // R^T R = I
        for i in 0..3 {
            for j in 0..3 {
                let dot: f64 = (0..3).map(|k| m[k][i] * m[k][j]).sum();
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(approx_eq(dot, expected, 1e-12), "R^T R [{i}][{j}] = {dot}");
            }
        }
    }

    #[test]
    fn test_from_rotation_matrix_round_trip() {
        let q0 = Quaternion::from_axis_angle([0.0, 1.0, 0.0], PI / 4.0);
        let m = q0.to_rotation_matrix();
        let q1 = Quaternion::from_rotation_matrix(&m);
        assert!(quat_approx_eq(q0, q1, 1e-12));
    }

    #[test]
    fn test_exp_ln_round_trip() {
        let q = Quaternion::new(0.0, 0.5, 0.3, 0.1); // pure quaternion
        let recovered = q.exp().ln();
        assert!(quat_approx_eq(q, recovered, 1e-12));
    }

    // ── DualNumber tests ──────────────────────────────────────────────────────

    #[test]
    fn test_dual_number_mul() {
        let a = DualNumber::new(2.0, 3.0);
        let b = DualNumber::new(4.0, 5.0);
        let c = a.mul(&b);
        assert!(approx_eq(c.real, 8.0, 1e-15));
        assert!(approx_eq(c.dual, 2.0 * 5.0 + 3.0 * 4.0, 1e-15));
    }

    // ── DualQuaternion tests ──────────────────────────────────────────────────

    #[test]
    fn test_identity_transforms_point() {
        let dq = DualQuaternion::identity();
        let p = [1.0, 2.0, 3.0];
        let tp = dq.transform_point(p);
        assert!(vec3_approx_eq(tp, p, 1e-12));
    }

    #[test]
    fn test_from_rotation_translation_recovers_translation() {
        let q = Quaternion::identity();
        let t = [3.0, -1.0, 2.0];
        let dq = DualQuaternion::from_rotation_translation(q, t);
        let tr = dq.to_translation();
        assert!(vec3_approx_eq(tr, t, 1e-12), "got {tr:?}");
    }

    #[test]
    fn test_pure_translation_transform() {
        let q = Quaternion::identity();
        let t = [5.0, 0.0, 0.0];
        let dq = DualQuaternion::from_rotation_translation(q, t);
        let p = [1.0, 2.0, 3.0];
        let tp = dq.transform_point(p);
        assert!(vec3_approx_eq(tp, [6.0, 2.0, 3.0], 1e-12), "got {tp:?}");
    }

    #[test]
    fn test_mul_identity_left() {
        let id = DualQuaternion::identity();
        let q = Quaternion::from_axis_angle([0.0, 1.0, 0.0], PI / 3.0);
        let dq = DualQuaternion::from_rotation_translation(q, [1.0, 2.0, 3.0]);
        let result = id.mul(&dq);
        assert!(dq_approx_eq(result.normalize(), dq.normalize(), 1e-12));
    }

    #[test]
    fn test_sclerp_endpoints() {
        let q1 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], 0.0);
        let dq1 = DualQuaternion::from_rotation_translation(q1, [0.0, 0.0, 0.0]);
        let q2 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let dq2 = DualQuaternion::from_rotation_translation(q2, [1.0, 0.0, 0.0]);

        let s0 = dq1.sclerp(&dq2, 0.0).normalize();
        let s1 = dq1.sclerp(&dq2, 1.0).normalize();

        // t=0 should recover dq1, t=1 should recover dq2.
        assert!(
            vec3_approx_eq(s0.to_translation(), dq1.to_translation(), 1e-8),
            "sclerp(0) translation: {:?}",
            s0.to_translation()
        );
        assert!(
            vec3_approx_eq(s1.to_translation(), dq2.to_translation(), 1e-8),
            "sclerp(1) translation: {:?}",
            s1.to_translation()
        );
    }

    #[test]
    fn test_rotation_translation_composition() {
        let q = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let t = [0.0, 0.0, 5.0];
        let dq = DualQuaternion::from_rotation_translation(q, t);
        let p = [1.0, 0.0, 0.0];
        let tp = dq.transform_point(p);
        // Rotate [1,0,0] by 90° around Z → [0,1,0], then translate by [0,0,5].
        assert!(vec3_approx_eq(tp, [0.0, 1.0, 5.0], 1e-12), "got {tp:?}");
    }

    // ── DLB (dual quaternion linear blending) ─────────────────────────────────

    #[test]
    fn test_dlb_single_quaternion() {
        let q = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 4.0);
        let dq = DualQuaternion::from_rotation_translation(q, [1.0, 0.0, 0.0]);
        let blended = dual_quaternion_blend(&[dq], &[1.0]);
        // Blending a single quaternion with weight=1 should return itself
        let p = [0.0, 0.0, 0.0];
        let tp1 = dq.transform_point(p);
        let tp2 = blended.transform_point(p);
        assert!(
            vec3_approx_eq(tp1, tp2, 1e-10),
            "DLB of single DQ should be identity, got {:?}",
            tp2
        );
    }

    #[test]
    fn test_dlb_two_equal_weights() {
        // Blending two identical DQs with equal weights should give the same DQ
        let q = Quaternion::from_axis_angle([0.0, 1.0, 0.0], PI / 6.0);
        let dq = DualQuaternion::from_rotation_translation(q, [0.0, 0.0, 0.0]);
        let blended = dual_quaternion_blend(&[dq, dq], &[0.5, 0.5]);
        let p = [1.0, 0.0, 0.0];
        let tp1 = dq.transform_point(p);
        let tp2 = blended.transform_point(p);
        assert!(
            vec3_approx_eq(tp1, tp2, 1e-10),
            "DLB of equal pair should equal original"
        );
    }

    #[test]
    fn test_dlb_empty_returns_identity() {
        let blended = dual_quaternion_blend(&[], &[]);
        let id = DualQuaternion::identity();
        let p = [1.0, 2.0, 3.0];
        assert!(vec3_approx_eq(
            blended.transform_point(p),
            id.transform_point(p),
            1e-12
        ));
    }

    // ── nLERP tests ────────────────────────────────────────────────────────────

    #[test]
    fn test_nlerp_endpoints() {
        let q1 = Quaternion::identity();
        let dq1 = DualQuaternion::from_rotation_translation(q1, [0.0, 0.0, 0.0]);
        let q2 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let dq2 = DualQuaternion::from_rotation_translation(q2, [2.0, 0.0, 0.0]);
        let s0 = dq1.nlerp(&dq2, 0.0);
        let s1 = dq1.nlerp(&dq2, 1.0);
        let p = [0.0, 0.0, 0.0];
        assert!(
            vec3_approx_eq(s0.transform_point(p), dq1.transform_point(p), 1e-10),
            "nlerp t=0 should recover dq1"
        );
        assert!(
            vec3_approx_eq(s1.transform_point(p), dq2.transform_point(p), 1e-10),
            "nlerp t=1 should recover dq2"
        );
    }

    // ── compose_transforms tests ───────────────────────────────────────────────

    #[test]
    fn test_compose_transforms() {
        // Compose: first translate by [1,0,0], then rotate 90° around Z.
        let t = DualQuaternion::from_rotation_translation(Quaternion::identity(), [1.0, 0.0, 0.0]);
        let r = DualQuaternion::from_rotation_translation(
            Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0),
            [0.0, 0.0, 0.0],
        );
        let combined = t.compose_transforms(&r);
        let p = [0.0, 0.0, 0.0];
        // translate [0,0,0] → [1,0,0], then rotate 90° → [0,1,0]
        let result = combined.transform_point(p);
        assert!(
            vec3_approx_eq(result, [0.0, 1.0, 0.0], 1e-10),
            "composed transform: expected [0,1,0], got {:?}",
            result
        );
    }

    // ── rotation_angle / rotation_axis tests ──────────────────────────────────

    #[test]
    fn test_rotation_angle_pi_over_3() {
        let angle = PI / 3.0;
        let q = Quaternion::from_axis_angle([0.0, 1.0, 0.0], angle);
        let dq = DualQuaternion::from_rotation_translation(q, [0.0, 0.0, 0.0]);
        let extracted = dq.rotation_angle();
        assert!(
            approx_eq(extracted, angle, 1e-12),
            "Expected angle {angle}, got {extracted}"
        );
    }

    #[test]
    fn test_rotation_axis_z() {
        let q = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 4.0);
        let dq = DualQuaternion::from_rotation_translation(q, [0.0, 0.0, 0.0]);
        let axis = dq.rotation_axis();
        assert!(
            vec3_approx_eq(axis, [0.0, 0.0, 1.0], 1e-10),
            "Axis should be [0,0,1], got {:?}",
            axis
        );
    }

    #[test]
    fn test_rotation_axis_identity() {
        let dq = DualQuaternion::identity();
        let axis = dq.rotation_axis();
        // Identity should return the default axis [0,0,1]
        assert!(vec3_approx_eq(axis, [0.0, 0.0, 1.0], 1e-10));
    }

    // ── sclerp_sequence tests ─────────────────────────────────────────────────

    #[test]
    fn test_sclerp_sequence_length() {
        let dq1 = DualQuaternion::identity();
        let q2 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let dq2 = DualQuaternion::from_rotation_translation(q2, [0.0, 0.0, 0.0]);
        let seq = sclerp_sequence(&dq1, &dq2, 5);
        assert_eq!(seq.len(), 5, "sequence should have 5 poses");
    }

    #[test]
    fn test_sclerp_sequence_endpoints() {
        let dq1 = DualQuaternion::identity();
        let q2 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let dq2 = DualQuaternion::from_rotation_translation(q2, [1.0, 0.0, 0.0]);
        let seq = sclerp_sequence(&dq1, &dq2, 10);
        let p = [1.0, 0.0, 0.0];
        let t0 = seq[0].transform_point(p);
        let t9 = seq[9].transform_point(p);
        let expected0 = dq1.transform_point(p);
        let expected9 = dq2.transform_point(p);
        assert!(
            vec3_approx_eq(t0, expected0, 1e-8),
            "seq[0] should match start"
        );
        assert!(
            vec3_approx_eq(t9, expected9, 1e-8),
            "seq[n-1] should match end"
        );
    }

    #[test]
    fn test_sclerp_sequence_empty() {
        let dq = DualQuaternion::identity();
        let seq = sclerp_sequence(&dq, &dq, 0);
        assert!(seq.is_empty());
    }

    // ── invert tests ──────────────────────────────────────────────────────────

    #[test]
    fn test_dq_invert_composition_is_identity() {
        let q = Quaternion::from_axis_angle([0.0, 1.0, 0.0], PI / 3.0);
        let dq = DualQuaternion::from_rotation_translation(q, [1.0, 2.0, 3.0]);
        let inv = dq.invert();
        let product = dq.mul(&inv).normalize();
        let id = DualQuaternion::identity();
        let p = [5.0, -2.0, 1.0];
        let tp1 = product.transform_point(p);
        let tp2 = id.transform_point(p);
        assert!(
            vec3_approx_eq(tp1, tp2, 1e-10),
            "dq * inv(dq) should be identity transform, got {:?}",
            tp1
        );
    }

    // ── transform_points tests ────────────────────────────────────────────────

    #[test]
    fn test_transform_points_batch() {
        let q = Quaternion::identity();
        let dq = DualQuaternion::from_rotation_translation(q, [1.0, 0.0, 0.0]);
        let points = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let transformed = dq.transform_points(&points);
        assert_eq!(transformed.len(), 3);
        // Each point should be shifted by [1,0,0]
        for (orig, trans) in points.iter().zip(transformed.iter()) {
            assert!((trans[0] - orig[0] - 1.0).abs() < 1e-12);
            assert!((trans[1] - orig[1]).abs() < 1e-12);
        }
    }

    // ── Plücker coordinate tests ──────────────────────────────────────────────

    #[test]
    fn test_plucker_from_points() {
        let p = [0.0, 0.0, 0.0];
        let q = [1.0, 0.0, 0.0];
        let pl = PluckerLine::from_two_points(p, q);
        // Direction should be unit x
        assert!((pl.direction[0] - 1.0).abs() < 1e-12);
        assert!(pl.direction[1].abs() < 1e-12);
        assert!(pl.direction[2].abs() < 1e-12);
    }

    #[test]
    fn test_plucker_moment_through_origin() {
        // A line through the origin has zero moment
        let p = [0.0, 0.0, 0.0];
        let q = [0.0, 1.0, 0.0];
        let pl = PluckerLine::from_two_points(p, q);
        let mom_sq = pl.moment[0].powi(2) + pl.moment[1].powi(2) + pl.moment[2].powi(2);
        assert!(mom_sq.sqrt() < 1e-12, "moment={:?}", pl.moment);
    }

    #[test]
    fn test_plucker_distance_parallel_lines() {
        // Two parallel lines along X, separated by 1 in Y
        let l1 = PluckerLine::from_two_points([0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        let l2 = PluckerLine::from_two_points([0.0, 1.0, 0.0], [1.0, 1.0, 0.0]);
        let d = l1.distance_to_line(&l2);
        assert!((d - 1.0).abs() < 1e-10, "distance={}", d);
    }

    // ── DH transform tests ────────────────────────────────────────────────────

    #[test]
    fn test_dh_transform_zero_params() {
        // All-zero DH params → identity
        let dq = dh_transform(0.0, 0.0, 0.0, 0.0);
        let p = [1.0, 0.0, 0.0];
        let tp = dq.transform_point(p);
        assert!(vec3_approx_eq(tp, p, 1e-10), "got {:?}", tp);
    }

    #[test]
    fn test_dh_transform_pure_translation_d() {
        // d=2, others=0 → pure translation along z by 2
        let dq = dh_transform(0.0, 0.0, 2.0, 0.0);
        let p = [0.0, 0.0, 0.0];
        let tp = dq.transform_point(p);
        assert!(vec3_approx_eq(tp, [0.0, 0.0, 2.0], 1e-10), "got {:?}", tp);
    }

    #[test]
    fn test_dh_transform_theta_rotation() {
        // theta=PI/2 around z, all others 0 → rotate x to y
        let dq = dh_transform(0.0, 0.0, 0.0, PI / 2.0);
        let p = [1.0, 0.0, 0.0];
        let tp = dq.transform_point(p);
        assert!(vec3_approx_eq(tp, [0.0, 1.0, 0.0], 1e-10), "got {:?}", tp);
    }

    // ── Velocity screw tests ──────────────────────────────────────────────────

    #[test]
    fn test_velocity_screw_zero() {
        let vs = VelocityScrew::zero();
        let v = vs.velocity_at_point([1.0, 0.0, 0.0]);
        assert!(vec3_approx_eq(v, [0.0, 0.0, 0.0], 1e-12));
    }

    #[test]
    fn test_velocity_screw_pure_translation() {
        let vs = VelocityScrew {
            angular: [0.0; 3],
            linear: [1.0, 0.0, 0.0],
        };
        let v = vs.velocity_at_point([0.0, 5.0, 0.0]);
        // Pure translation: velocity = linear everywhere
        assert!(vec3_approx_eq(v, [1.0, 0.0, 0.0], 1e-12), "got {:?}", v);
    }

    #[test]
    fn test_velocity_screw_pure_rotation_z() {
        // Rotation around z at 1 rad/s: point (1,0,0) should have vel (0,1,0)
        let vs = VelocityScrew {
            angular: [0.0, 0.0, 1.0],
            linear: [0.0, 0.0, 0.0],
        };
        let v = vs.velocity_at_point([1.0, 0.0, 0.0]);
        assert!(vec3_approx_eq(v, [0.0, 1.0, 0.0], 1e-12), "got {:?}", v);
    }

    #[test]
    fn test_velocity_screw_addition() {
        let a = VelocityScrew {
            angular: [1.0, 0.0, 0.0],
            linear: [0.0, 1.0, 0.0],
        };
        let b = VelocityScrew {
            angular: [0.0, 1.0, 0.0],
            linear: [0.0, 0.0, 1.0],
        };
        let s = a.add(&b);
        assert!(vec3_approx_eq(s.angular, [1.0, 1.0, 0.0], 1e-12));
        assert!(vec3_approx_eq(s.linear, [0.0, 1.0, 1.0], 1e-12));
    }

    // ── Rigid body pose blending tests ────────────────────────────────────────

    #[test]
    fn test_pose_blend_two_equal() {
        let q = Quaternion::from_axis_angle([0.0, 1.0, 0.0], PI / 4.0);
        let dq = DualQuaternion::from_rotation_translation(q, [1.0, 2.0, 3.0]);
        let blended = rigid_body_pose_blend(&[dq, dq], &[0.5, 0.5]);
        let p = [1.0, 0.0, 0.0];
        let tp1 = dq.transform_point(p);
        let tp2 = blended.transform_point(p);
        assert!(vec3_approx_eq(tp1, tp2, 1e-8), "blend of equal = original");
    }

    #[test]
    fn test_pose_blend_weighted_translation() {
        let dq1 =
            DualQuaternion::from_rotation_translation(Quaternion::identity(), [0.0, 0.0, 0.0]);
        let dq2 =
            DualQuaternion::from_rotation_translation(Quaternion::identity(), [2.0, 0.0, 0.0]);
        // 50/50 blend should give translation [1, 0, 0]
        let blended = rigid_body_pose_blend(&[dq1, dq2], &[0.5, 0.5]);
        let t = blended.to_translation();
        assert!((t[0] - 1.0).abs() < 1e-8, "t[0]={}", t[0]);
    }

    // ── Screw interpolation chain tests ──────────────────────────────────────

    #[test]
    fn test_sclerp_sequence_midpoint() {
        let dq1 = DualQuaternion::identity();
        let q2 = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 2.0);
        let dq2 = DualQuaternion::from_rotation_translation(q2, [0.0, 0.0, 0.0]);
        let seq = sclerp_sequence(&dq1, &dq2, 3);
        // seq[1] at t=0.5 should have angle ~ PI/4
        let angle = seq[1].rotation_angle();
        assert!((angle - PI / 4.0).abs() < 1e-6, "midpoint angle={}", angle);
    }

    #[test]
    fn test_dq_mul_sequence() {
        // Apply same rotation 4 times → 4x rotation
        let q = Quaternion::from_axis_angle([0.0, 0.0, 1.0], PI / 8.0);
        let dq = DualQuaternion::from_rotation_translation(q, [0.0, 0.0, 0.0]);
        let dq4 = dq.mul(&dq).mul(&dq).mul(&dq);
        let angle = dq4.rotation_angle();
        assert!(
            (angle - PI / 2.0).abs() < 1e-10,
            "4x rotation angle={}",
            angle
        );
    }
}
