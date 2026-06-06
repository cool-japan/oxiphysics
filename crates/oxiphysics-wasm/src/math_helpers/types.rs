//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::functions::*;
use crate::types::{QuatWasm, TransformWasm, Vec3Wasm};
use serde::{Deserialize, Serialize};

/// A flat `[x, y, z]` representation of a 3D vector for JavaScript interop.
///
/// Use `JsVec3::from_vec3` / `JsVec3::to_vec3` to convert between this and
/// the richer [`Vec3Wasm`] type.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::math_helpers::JsVec3;
///
/// let js = JsVec3::new(1.0, 2.0, 3.0);
/// assert_eq!(js.to_array(), [1.0, 2.0, 3.0]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsVec3 {
    /// X component.
    pub x: f64,
    /// Y component.
    pub y: f64,
    /// Z component.
    pub z: f64,
}
impl JsVec3 {
    /// Create a new `JsVec3` from components.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
    /// Create the zero vector.
    pub fn zero() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
    /// Create from a `[f64; 3]` array.
    pub fn from_array(arr: [f64; 3]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
            z: arr[2],
        }
    }
    /// Convert to a `[f64; 3]` array (for JavaScript `Float64Array` slices).
    pub fn to_array(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
    /// Convert from a [`Vec3Wasm`].
    pub fn from_vec3(v: &Vec3Wasm) -> Self {
        Self {
            x: v.x,
            y: v.y,
            z: v.z,
        }
    }
    /// Convert to a [`Vec3Wasm`].
    pub fn to_vec3(&self) -> Vec3Wasm {
        Vec3Wasm::new(self.x, self.y, self.z)
    }
    /// Return the squared length.
    pub fn length_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }
    /// Return the length.
    pub fn length(&self) -> f64 {
        self.length_squared().sqrt()
    }
    /// Return a normalized copy, or zero if length is zero.
    pub fn normalized(&self) -> Self {
        let len = self.length();
        if len < 1e-15 {
            Self::zero()
        } else {
            Self {
                x: self.x / len,
                y: self.y / len,
                z: self.z / len,
            }
        }
    }
    /// Dot product.
    pub fn dot(&self, other: &Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
    /// Cross product.
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            x: self.y * other.z - self.z * other.y,
            y: self.z * other.x - self.x * other.z,
            z: self.x * other.y - self.y * other.x,
        }
    }
    /// Component-wise addition.
    pub fn add(&self, other: &Self) -> Self {
        Self {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
    /// Component-wise subtraction.
    pub fn sub(&self, other: &Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
    /// Scalar multiply.
    pub fn scale(&self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
    /// Negate.
    pub fn neg(&self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
    /// Linear interpolation (t clamped to \[0, 1\]).
    pub fn lerp(&self, other: &Self, t: f64) -> Self {
        let t = t.clamp(0.0, 1.0);
        Self {
            x: self.x + t * (other.x - self.x),
            y: self.y + t * (other.y - self.y),
            z: self.z + t * (other.z - self.z),
        }
    }
    /// Distance to another vector.
    pub fn distance_to(&self, other: &Self) -> f64 {
        self.sub(other).length()
    }
    /// Pack multiple `JsVec3` values into a flat `Vec`f64` interleaved `\[x0,y0,z0, x1,y1,z1,...\]`.
    pub fn pack_flat(vecs: &[Self]) -> Vec<f64> {
        let mut out = Vec::with_capacity(vecs.len() * 3);
        for v in vecs {
            out.push(v.x);
            out.push(v.y);
            out.push(v.z);
        }
        out
    }
    /// Unpack a flat `\[x0,y0,z0, x1,y1,z1,...\]` slice into a `Vec`JsVec3`.
    ///
    /// Panics if `data.len()` is not divisible by 3.
    pub fn unpack_flat(data: &[f64]) -> Vec<Self> {
        assert!(
            data.len().is_multiple_of(3),
            "flat vec3 data length must be divisible by 3"
        );
        data.chunks_exact(3)
            .map(|c| Self::new(c[0], c[1], c[2]))
            .collect()
    }
}
/// A 3-D ray defined by an origin and a (non-normalized) direction.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Ray origin.
    pub origin: [f64; 3],
    /// Ray direction (need not be unit-length).
    pub direction: [f64; 3],
}
impl Ray {
    /// Create a new ray.
    pub fn new(origin: [f64; 3], direction: [f64; 3]) -> Self {
        Self { origin, direction }
    }
    /// Test intersection with an AABB `(min, max)`.
    ///
    /// Returns `Some(t_min)` where `t_min >= 0` is the entry distance, or
    /// `None` if the ray misses.
    pub fn intersect_aabb(&self, min: [f64; 3], max: [f64; 3]) -> Option<f64> {
        let mut t_min = 0.0f64;
        let mut t_max = f64::INFINITY;
        for i in 0..3 {
            let inv_d = 1.0 / self.direction[i];
            let mut t0 = (min[i] - self.origin[i]) * inv_d;
            let mut t1 = (max[i] - self.origin[i]) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }
            t_min = t_min.max(t0);
            t_max = t_max.min(t1);
            if t_max < t_min {
                return None;
            }
        }
        if t_min >= 0.0 {
            Some(t_min)
        } else if t_max >= 0.0 {
            Some(0.0)
        } else {
            None
        }
    }
    /// Test intersection with a sphere at `center` with `radius`.
    /// Returns `Some(t)` (entry distance), or `None` if miss.
    pub fn intersect_sphere(&self, center: [f64; 3], radius: f64) -> Option<f64> {
        let oc = [
            self.origin[0] - center[0],
            self.origin[1] - center[1],
            self.origin[2] - center[2],
        ];
        let d = self.direction;
        let a = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        let b = 2.0 * (oc[0] * d[0] + oc[1] * d[1] + oc[2] * d[2]);
        let c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - radius * radius;
        let disc = b * b - 4.0 * a * c;
        if disc < 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        let t0 = (-b - sq) / (2.0 * a);
        let t1 = (-b + sq) / (2.0 * a);
        if t0 >= 0.0 {
            Some(t0)
        } else if t1 >= 0.0 {
            Some(t1)
        } else {
            None
        }
    }
    /// Evaluate a point along the ray at parameter `t`.
    pub fn at(&self, t: f64) -> [f64; 3] {
        [
            self.origin[0] + t * self.direction[0],
            self.origin[1] + t * self.direction[1],
            self.origin[2] + t * self.direction[2],
        ]
    }
}
/// A JS-friendly rigid body transform: flat `[px, py, pz, qx, qy, qz, qw]`.
///
/// Can be converted to/from a column-major 4×4 matrix (for WebGL / Three.js).
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::math_helpers::JsTransform;
///
/// let t = JsTransform::from_position(1.0, 2.0, 3.0);
/// let arr = t.to_array7();
/// assert!((arr[0] - 1.0).abs() < 1e-10);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsTransform {
    /// Position.
    pub position: JsVec3,
    /// Orientation (unit quaternion).
    pub rotation: JsQuat,
}
impl JsTransform {
    /// Create from position and rotation.
    pub fn new(position: JsVec3, rotation: JsQuat) -> Self {
        Self { position, rotation }
    }
    /// Create at the given position with identity rotation.
    pub fn from_position(x: f64, y: f64, z: f64) -> Self {
        Self {
            position: JsVec3::new(x, y, z),
            rotation: JsQuat::identity(),
        }
    }
    /// The identity transform (origin, no rotation).
    pub fn identity() -> Self {
        Self {
            position: JsVec3::zero(),
            rotation: JsQuat::identity(),
        }
    }
    /// Convert from a [`TransformWasm`].
    pub fn from_transform(t: &TransformWasm) -> Self {
        Self {
            position: JsVec3::from_vec3(&t.position),
            rotation: JsQuat::from_quat(&t.rotation),
        }
    }
    /// Convert to a [`TransformWasm`].
    pub fn to_transform(&self) -> TransformWasm {
        TransformWasm::new(self.position.to_vec3(), self.rotation.to_quat())
    }
    /// Encode as a flat `[px, py, pz, qx, qy, qz, qw]` array.
    pub fn to_array7(&self) -> [f64; 7] {
        [
            self.position.x,
            self.position.y,
            self.position.z,
            self.rotation.x,
            self.rotation.y,
            self.rotation.z,
            self.rotation.w,
        ]
    }
    /// Decode from a flat `[px, py, pz, qx, qy, qz, qw]` array.
    pub fn from_array7(arr: [f64; 7]) -> Self {
        Self {
            position: JsVec3::new(arr[0], arr[1], arr[2]),
            rotation: JsQuat::new(arr[3], arr[4], arr[5], arr[6]),
        }
    }
    /// Encode as a column-major 4×4 matrix `[f64; 16]` (WebGL / Three.js compatible).
    pub fn to_matrix4(&self) -> [f64; 16] {
        let q = &self.rotation;
        let (x, y, z, w) = (q.x, q.y, q.z, q.w);
        let (tx, ty, tz) = (self.position.x, self.position.y, self.position.z);
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y + w * z),
            2.0 * (x * z - w * y),
            0.0,
            2.0 * (x * y - w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z + w * x),
            0.0,
            2.0 * (x * z + w * y),
            2.0 * (y * z - w * x),
            1.0 - 2.0 * (x * x + y * y),
            0.0,
            tx,
            ty,
            tz,
            1.0,
        ]
    }
    /// Transform a point from local space to world space.
    pub fn transform_point(&self, p: &JsVec3) -> JsVec3 {
        self.rotation.rotate_vec(p).add(&self.position)
    }
    /// Transform a direction vector (ignores translation).
    pub fn transform_vector(&self, v: &JsVec3) -> JsVec3 {
        self.rotation.rotate_vec(v)
    }
    /// Compute the inverse transform.
    pub fn inverse(&self) -> Self {
        let inv_rot = self.rotation.conjugate();
        let neg_pos = self.position.neg();
        let inv_pos = inv_rot.rotate_vec(&neg_pos);
        Self {
            position: inv_pos,
            rotation: inv_rot,
        }
    }
    /// Pack multiple `JsTransform` values into a flat `Vec`f64` of 7 floats each.
    pub fn pack_flat(transforms: &[Self]) -> Vec<f64> {
        let mut out = Vec::with_capacity(transforms.len() * 7);
        for t in transforms {
            let arr = t.to_array7();
            out.extend_from_slice(&arr);
        }
        out
    }
    /// Unpack a flat slice (7 floats per transform) into a `Vec<JsTransform>`.
    ///
    /// Panics if `data.len()` is not divisible by 7.
    pub fn unpack_flat(data: &[f64]) -> Vec<Self> {
        assert!(
            data.len().is_multiple_of(7),
            "flat transform data length must be divisible by 7"
        );
        data.chunks_exact(7)
            .map(|c| Self::from_array7([c[0], c[1], c[2], c[3], c[4], c[5], c[6]]))
            .collect()
    }
}
/// A flat `\[x, y, z, w\]` representation of a unit quaternion for JavaScript interop.
///
/// Follows the convention where `w` is the scalar part and `(x, y, z)` is the
/// imaginary vector part. The identity quaternion is `\[0, 0, 0, 1\]`.
///
/// # Example
///
/// ```no_run
/// use oxiphysics_wasm::math_helpers::JsQuat;
///
/// let q = JsQuat::identity();
/// assert!((q.w - 1.0).abs() < 1e-10);
/// assert_eq!(q.to_array(), [0.0, 0.0, 0.0, 1.0]);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct JsQuat {
    /// X imaginary component.
    pub x: f64,
    /// Y imaginary component.
    pub y: f64,
    /// Z imaginary component.
    pub z: f64,
    /// W real (scalar) component.
    pub w: f64,
}
impl JsQuat {
    /// Create a quaternion from components.
    pub fn new(x: f64, y: f64, z: f64, w: f64) -> Self {
        Self { x, y, z, w }
    }
    /// The identity quaternion `(0, 0, 0, 1)`.
    pub fn identity() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            w: 1.0,
        }
    }
    /// Create from a `\[f64; 4\]` array `\[x, y, z, w\]`.
    pub fn from_array(arr: [f64; 4]) -> Self {
        Self {
            x: arr[0],
            y: arr[1],
            z: arr[2],
            w: arr[3],
        }
    }
    /// Convert to a `\[f64; 4\]` array `\[x, y, z, w\]`.
    pub fn to_array(&self) -> [f64; 4] {
        [self.x, self.y, self.z, self.w]
    }
    /// Convert from a [`QuatWasm`].
    pub fn from_quat(q: &QuatWasm) -> Self {
        Self {
            x: q.x,
            y: q.y,
            z: q.z,
            w: q.w,
        }
    }
    /// Convert to a [`QuatWasm`].
    pub fn to_quat(&self) -> QuatWasm {
        QuatWasm::new(self.x, self.y, self.z, self.w)
    }
    /// Normalize the quaternion (so it has unit length).
    pub fn normalized(&self) -> Self {
        let len = (self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w).sqrt();
        if len < 1e-15 {
            Self::identity()
        } else {
            Self {
                x: self.x / len,
                y: self.y / len,
                z: self.z / len,
                w: self.w / len,
            }
        }
    }
    /// Conjugate (inverse for unit quaternions).
    pub fn conjugate(&self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            w: self.w,
        }
    }
    /// Hamilton (quaternion) product.
    pub fn mul(&self, other: &Self) -> Self {
        Self {
            x: self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            y: self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            z: self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
            w: self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
        }
    }
    /// Rotate a 3D vector by this quaternion.
    pub fn rotate_vec(&self, v: &JsVec3) -> JsVec3 {
        let qv = JsVec3::new(self.x, self.y, self.z);
        let uv = qv.cross(v);
        let uuv = qv.cross(&uv);
        v.add(&uv.scale(2.0 * self.w)).add(&uuv.scale(2.0))
    }
    /// Build a rotation quaternion from an axis and angle (radians).
    ///
    /// `axis` need not be normalized (this function normalizes internally).
    pub fn from_axis_angle(axis: &JsVec3, angle_rad: f64) -> Self {
        let n = axis.normalized();
        let half = angle_rad * 0.5;
        let s = half.sin();
        Self {
            x: n.x * s,
            y: n.y * s,
            z: n.z * s,
            w: half.cos(),
        }
    }
    /// Build a quaternion from Euler angles (ZYX convention, all in radians).
    pub fn from_euler_zyx(roll: f64, pitch: f64, yaw: f64) -> Self {
        let (sr, cr) = (roll * 0.5).sin_cos();
        let (sp, cp) = (pitch * 0.5).sin_cos();
        let (sy, cy) = (yaw * 0.5).sin_cos();
        Self {
            x: sr * cp * cy - cr * sp * sy,
            y: cr * sp * cy + sr * cp * sy,
            z: cr * cp * sy - sr * sp * cy,
            w: cr * cp * cy + sr * sp * sy,
        }
        .normalized()
    }
    /// Spherical linear interpolation toward `other` by factor `t` (clamped [0, 1]).
    pub fn slerp(&self, other: &Self, t: f64) -> Self {
        let t = t.clamp(0.0, 1.0);
        let mut dot = self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w;
        let other_adj = if dot < 0.0 {
            dot = -dot;
            Self::new(-other.x, -other.y, -other.z, -other.w)
        } else {
            *other
        };
        if dot > 0.9995 {
            return Self {
                x: self.x + t * (other_adj.x - self.x),
                y: self.y + t * (other_adj.y - self.y),
                z: self.z + t * (other_adj.z - self.z),
                w: self.w + t * (other_adj.w - self.w),
            }
            .normalized();
        }
        let theta_0 = dot.acos();
        let theta = theta_0 * t;
        let sin_theta = theta.sin();
        let sin_theta_0 = theta_0.sin();
        let s1 = theta.cos() - dot * sin_theta / sin_theta_0;
        let s2 = sin_theta / sin_theta_0;
        Self {
            x: s1 * self.x + s2 * other_adj.x,
            y: s1 * self.y + s2 * other_adj.y,
            z: s1 * self.z + s2 * other_adj.z,
            w: s1 * self.w + s2 * other_adj.w,
        }
    }
    /// Pack multiple `JsQuat` values into a flat `Vec`f64` `[x0,y0,z0,w0, x1,...]`.
    pub fn pack_flat(quats: &[Self]) -> Vec<f64> {
        let mut out = Vec::with_capacity(quats.len() * 4);
        for q in quats {
            out.push(q.x);
            out.push(q.y);
            out.push(q.z);
            out.push(q.w);
        }
        out
    }
    /// Unpack a flat `[x0,y0,z0,w0, x1,...]` slice into a `Vec`JsQuat`.
    ///
    /// Panics if `data.len()` is not divisible by 4.
    pub fn unpack_flat(data: &[f64]) -> Vec<Self> {
        assert!(
            data.len().is_multiple_of(4),
            "flat quat data length must be divisible by 4"
        );
        data.chunks_exact(4)
            .map(|c| Self::new(c[0], c[1], c[2], c[3]))
            .collect()
    }
}
impl JsQuat {
    /// Convert this quaternion to Euler angles `\[roll, pitch, yaw\]` (ZYX convention, radians).
    pub fn to_euler_zyx(&self) -> [f64; 3] {
        let q = self.normalized();
        let sinr_cosp = 2.0 * (q.w * q.x + q.y * q.z);
        let cosr_cosp = 1.0 - 2.0 * (q.x * q.x + q.y * q.y);
        let roll = sinr_cosp.atan2(cosr_cosp);
        let sinp = 2.0 * (q.w * q.y - q.z * q.x);
        let pitch = if sinp.abs() >= 1.0 {
            std::f64::consts::FRAC_PI_2.copysign(sinp)
        } else {
            sinp.asin()
        };
        let siny_cosp = 2.0 * (q.w * q.z + q.x * q.y);
        let cosy_cosp = 1.0 - 2.0 * (q.y * q.y + q.z * q.z);
        let yaw = siny_cosp.atan2(cosy_cosp);
        [roll, pitch, yaw]
    }
    /// Compute the angular difference (angle in radians) between this and `other`.
    pub fn angle_to(&self, other: &Self) -> f64 {
        let dot = (self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w)
            .abs()
            .min(1.0);
        2.0 * dot.acos()
    }
    /// Compute the squared norm (for sanity checks).
    pub fn norm_squared(&self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z + self.w * self.w
    }
    /// Build a quaternion that rotates `from` direction to `to` direction.
    pub fn from_rotation_between(from: &JsVec3, to: &JsVec3) -> Self {
        let f = from.normalized();
        let t = to.normalized();
        let dot = f.dot(&t).clamp(-1.0, 1.0);
        if dot > 0.9999 {
            return Self::identity();
        }
        if dot < -0.9999 {
            let axis = if f.x.abs() < 0.9 {
                f.cross(&JsVec3::new(1.0, 0.0, 0.0)).normalized()
            } else {
                f.cross(&JsVec3::new(0.0, 1.0, 0.0)).normalized()
            };
            return Self::from_axis_angle(&axis, std::f64::consts::PI);
        }
        let axis = f.cross(&t);
        let s = ((1.0 + dot) * 2.0).sqrt();
        let inv_s = 1.0 / s;
        Self {
            x: axis.x * inv_s,
            y: axis.y * inv_s,
            z: axis.z * inv_s,
            w: s * 0.5,
        }
        .normalized()
    }
    /// Convert to a 3×3 rotation matrix (row-major `\[f64; 9\]`).
    pub fn to_rotation_matrix3(&self) -> [f64; 9] {
        let q = self.normalized();
        let (x, y, z, w) = (q.x, q.y, q.z, q.w);
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
            2.0 * (x * y + w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - w * x),
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            1.0 - 2.0 * (x * x + y * y),
        ]
    }
    /// Nlerp (normalized linear interpolation) — cheaper than slerp, good for small angles.
    pub fn nlerp(&self, other: &Self, t: f64) -> Self {
        let t = t.clamp(0.0, 1.0);
        let dot = self.x * other.x + self.y * other.y + self.z * other.z + self.w * other.w;
        let flip = if dot < 0.0 { -1.0 } else { 1.0 };
        Self {
            x: self.x + t * (flip * other.x - self.x),
            y: self.y + t * (flip * other.y - self.y),
            z: self.z + t * (flip * other.z - self.z),
            w: self.w + t * (flip * other.w - self.w),
        }
        .normalized()
    }
}
/// Collection of common easing functions mapping `t ∈ \[0, 1\]` to `\[0, 1\]`.
pub struct Easing;
impl Easing {
    /// Linear (no easing).
    pub fn linear(t: f64) -> f64 {
        t
    }
    /// Quadratic ease-in.
    pub fn ease_in_quad(t: f64) -> f64 {
        t * t
    }
    /// Quadratic ease-out.
    pub fn ease_out_quad(t: f64) -> f64 {
        t * (2.0 - t)
    }
    /// Quadratic ease-in-out.
    pub fn ease_in_out_quad(t: f64) -> f64 {
        if t < 0.5 {
            2.0 * t * t
        } else {
            -1.0 + (4.0 - 2.0 * t) * t
        }
    }
    /// Cubic ease-in.
    pub fn ease_in_cubic(t: f64) -> f64 {
        t * t * t
    }
    /// Cubic ease-out.
    pub fn ease_out_cubic(t: f64) -> f64 {
        let t1 = t - 1.0;
        t1 * t1 * t1 + 1.0
    }
    /// Cubic ease-in-out.
    pub fn ease_in_out_cubic(t: f64) -> f64 {
        if t < 0.5 {
            4.0 * t * t * t
        } else {
            let t1 = 2.0 * t - 2.0;
            0.5 * t1 * t1 * t1 + 1.0
        }
    }
    /// Quartic ease-in.
    pub fn ease_in_quart(t: f64) -> f64 {
        t * t * t * t
    }
    /// Quartic ease-out.
    pub fn ease_out_quart(t: f64) -> f64 {
        let t1 = t - 1.0;
        1.0 - t1 * t1 * t1 * t1
    }
    /// Sinusoidal ease-in.
    pub fn ease_in_sine(t: f64) -> f64 {
        1.0 - ((t * std::f64::consts::PI * 0.5).cos())
    }
    /// Sinusoidal ease-out.
    pub fn ease_out_sine(t: f64) -> f64 {
        (t * std::f64::consts::PI * 0.5).sin()
    }
    /// Sinusoidal ease-in-out.
    pub fn ease_in_out_sine(t: f64) -> f64 {
        -0.5 * ((std::f64::consts::PI * t).cos() - 1.0)
    }
    /// Exponential ease-in.
    pub fn ease_in_expo(t: f64) -> f64 {
        if t == 0.0 {
            0.0
        } else {
            (2.0f64).powf(10.0 * t - 10.0)
        }
    }
    /// Exponential ease-out.
    pub fn ease_out_expo(t: f64) -> f64 {
        if t == 1.0 {
            1.0
        } else {
            1.0 - (2.0f64).powf(-10.0 * t)
        }
    }
    /// Elastic ease-out.
    pub fn ease_out_elastic(t: f64) -> f64 {
        if t == 0.0 {
            return 0.0;
        }
        if t == 1.0 {
            return 1.0;
        }
        let c4 = (2.0 * std::f64::consts::PI) / 3.0;
        (2.0f64).powf(-10.0 * t) * ((t * 10.0 - 0.75) * c4).sin() + 1.0
    }
    /// Bounce ease-out.
    pub fn ease_out_bounce(t: f64) -> f64 {
        let n1 = 7.5625;
        let d1 = 2.75;
        if t < 1.0 / d1 {
            n1 * t * t
        } else if t < 2.0 / d1 {
            let t2 = t - 1.5 / d1;
            n1 * t2 * t2 + 0.75
        } else if t < 2.5 / d1 {
            let t2 = t - 2.25 / d1;
            n1 * t2 * t2 + 0.9375
        } else {
            let t2 = t - 2.625 / d1;
            n1 * t2 * t2 + 0.984375
        }
    }
    /// Back ease-in (slight overshoot backward).
    pub fn ease_in_back(t: f64) -> f64 {
        let c1 = 1.70158;
        let c3 = c1 + 1.0;
        c3 * t * t * t - c1 * t * t
    }
}
/// A column-major 4×4 matrix (compatible with WebGL / Three.js / WGSL).
///
/// Stored as `\[f64; 16\]` in column-major order:
/// `m\[col * 4 + row\]`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Mat4(pub [f64; 16]);
impl Mat4 {
    /// Identity matrix.
    pub fn identity() -> Self {
        Self([
            1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
        ])
    }
    /// Zero matrix.
    pub fn zero() -> Self {
        Self([0.0; 16])
    }
    /// Build a translation matrix from `\[tx, ty, tz\]`.
    pub fn from_translation(t: [f64; 3]) -> Self {
        let mut m = Self::identity();
        m.0[12] = t[0];
        m.0[13] = t[1];
        m.0[14] = t[2];
        m
    }
    /// Build a uniform scale matrix.
    pub fn from_scale(s: f64) -> Self {
        let mut m = Self::zero();
        m.0[0] = s;
        m.0[5] = s;
        m.0[10] = s;
        m.0[15] = 1.0;
        m
    }
    /// Build a non-uniform scale matrix.
    pub fn from_scale_xyz(sx: f64, sy: f64, sz: f64) -> Self {
        let mut m = Self::zero();
        m.0[0] = sx;
        m.0[5] = sy;
        m.0[10] = sz;
        m.0[15] = 1.0;
        m
    }
    /// Build a rotation matrix from a quaternion `\[x, y, z, w\]`.
    pub fn from_quat(qx: f64, qy: f64, qz: f64, qw: f64) -> Self {
        let x2 = qx + qx;
        let y2 = qy + qy;
        let z2 = qz + qz;
        let xx = qx * x2;
        let xy = qx * y2;
        let xz = qx * z2;
        let yy = qy * y2;
        let yz = qy * z2;
        let zz = qz * z2;
        let wx = qw * x2;
        let wy = qw * y2;
        let wz = qw * z2;
        Self([
            1.0 - (yy + zz),
            xy + wz,
            xz - wy,
            0.0,
            xy - wz,
            1.0 - (xx + zz),
            yz + wx,
            0.0,
            xz + wy,
            yz - wx,
            1.0 - (xx + yy),
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
        ])
    }
    /// Multiply two 4×4 matrices (this * other).
    pub fn mul(&self, other: &Self) -> Self {
        let a = &self.0;
        let b = &other.0;
        let mut c = [0.0f64; 16];
        for col in 0..4 {
            for row in 0..4 {
                let mut sum = 0.0;
                for k in 0..4 {
                    sum += a[k * 4 + row] * b[col * 4 + k];
                }
                c[col * 4 + row] = sum;
            }
        }
        Self(c)
    }
    /// Transform a 4-component vector `\[x, y, z, w\]` by this matrix.
    pub fn transform_vec4(&self, v: [f64; 4]) -> [f64; 4] {
        let m = &self.0;
        [
            m[0] * v[0] + m[4] * v[1] + m[8] * v[2] + m[12] * v[3],
            m[1] * v[0] + m[5] * v[1] + m[9] * v[2] + m[13] * v[3],
            m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3],
            m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3],
        ]
    }
    /// Transform a point `\[x, y, z\]` (w=1).
    pub fn transform_point3(&self, p: [f64; 3]) -> [f64; 3] {
        let r = self.transform_vec4([p[0], p[1], p[2], 1.0]);
        let w_inv = if r[3].abs() > 1e-15 { 1.0 / r[3] } else { 1.0 };
        [r[0] * w_inv, r[1] * w_inv, r[2] * w_inv]
    }
    /// Transform a direction `\[x, y, z\]` (w=0, ignores translation).
    pub fn transform_dir3(&self, d: [f64; 3]) -> [f64; 3] {
        let r = self.transform_vec4([d[0], d[1], d[2], 0.0]);
        [r[0], r[1], r[2]]
    }
    /// Transpose this matrix.
    pub fn transpose(&self) -> Self {
        let m = &self.0;
        Self([
            m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13], m[2], m[6], m[10], m[14], m[3], m[7],
            m[11], m[15],
        ])
    }
    /// Build a perspective projection matrix (right-handed, depth 0..1 NDC).
    ///
    /// - `fov_y_rad`: vertical field of view in radians.
    /// - `aspect`: width / height.
    /// - `near`, `far`: near/far clip distances.
    pub fn perspective(fov_y_rad: f64, aspect: f64, near: f64, far: f64) -> Self {
        let f = 1.0 / (fov_y_rad * 0.5).tan();
        let range_inv = 1.0 / (near - far);
        Self([
            f / aspect,
            0.0,
            0.0,
            0.0,
            0.0,
            f,
            0.0,
            0.0,
            0.0,
            0.0,
            (far) * range_inv,
            -1.0,
            0.0,
            0.0,
            near * far * range_inv,
            0.0,
        ])
    }
    /// Build an orthographic projection matrix (right-handed, depth 0..1 NDC).
    pub fn orthographic(left: f64, right: f64, bottom: f64, top: f64, near: f64, far: f64) -> Self {
        let rl = 1.0 / (right - left);
        let tb = 1.0 / (top - bottom);
        let nf = 1.0 / (near - far);
        Self([
            2.0 * rl,
            0.0,
            0.0,
            0.0,
            0.0,
            2.0 * tb,
            0.0,
            0.0,
            0.0,
            0.0,
            nf,
            0.0,
            -(right + left) * rl,
            -(top + bottom) * tb,
            near * nf,
            1.0,
        ])
    }
    /// Build a look-at view matrix (right-handed).
    pub fn look_at(eye: [f64; 3], target: [f64; 3], up: [f64; 3]) -> Self {
        let e = JsVec3::from_array(eye);
        let t = JsVec3::from_array(target);
        let u = JsVec3::from_array(up);
        let f = t.sub(&e).normalized();
        let r = f.cross(&u.normalized()).normalized();
        let u2 = r.cross(&f);
        Self([
            r.x,
            u2.x,
            -f.x,
            0.0,
            r.y,
            u2.y,
            -f.y,
            0.0,
            r.z,
            u2.z,
            -f.z,
            0.0,
            -r.dot(&e),
            -u2.dot(&e),
            f.dot(&e),
            1.0,
        ])
    }
    /// Return the raw `\[f64; 16\]` backing array.
    pub fn to_array(&self) -> [f64; 16] {
        self.0
    }
}
/// A view frustum defined by six half-space planes (left, right, bottom, top, near, far).
///
/// Each plane is represented as `\[nx, ny, nz, d\]` where the plane equation is
/// `nx*x + ny*y + nz*z + d >= 0` for points inside.
#[derive(Debug, Clone, Copy)]
pub struct Frustum {
    /// Six planes: [left, right, bottom, top, near, far].
    pub planes: [[f64; 4]; 6],
}
impl Frustum {
    /// Extract frustum planes from a combined view-projection matrix.
    pub fn from_view_proj(vp: &Mat4) -> Self {
        let m = &vp.0;
        let row = |r: usize| -> [f64; 4] { [m[r], m[4 + r], m[8 + r], m[12 + r]] };
        let r0 = row(0);
        let r1 = row(1);
        let r2 = row(2);
        let r3 = row(3);
        let planes = [
            normalize_plane([r3[0] + r0[0], r3[1] + r0[1], r3[2] + r0[2], r3[3] + r0[3]]),
            normalize_plane([r3[0] - r0[0], r3[1] - r0[1], r3[2] - r0[2], r3[3] - r0[3]]),
            normalize_plane([r3[0] + r1[0], r3[1] + r1[1], r3[2] + r1[2], r3[3] + r1[3]]),
            normalize_plane([r3[0] - r1[0], r3[1] - r1[1], r3[2] - r1[2], r3[3] - r1[3]]),
            normalize_plane([r3[0] + r2[0], r3[1] + r2[1], r3[2] + r2[2], r3[3] + r2[3]]),
            normalize_plane([r3[0] - r2[0], r3[1] - r2[1], r3[2] - r2[2], r3[3] - r2[3]]),
        ];
        Self { planes }
    }
    /// Test if an AABB `(min, max)` intersects the frustum.
    /// Returns `true` if the AABB is at least partially inside.
    pub fn intersects_aabb(&self, min: [f64; 3], max: [f64; 3]) -> bool {
        for plane in &self.planes {
            let [nx, ny, nz, d] = *plane;
            let px = if nx >= 0.0 { max[0] } else { min[0] };
            let py = if ny >= 0.0 { max[1] } else { min[1] };
            let pz = if nz >= 0.0 { max[2] } else { min[2] };
            if nx * px + ny * py + nz * pz + d < 0.0 {
                return false;
            }
        }
        true
    }
    /// Test if a sphere is fully or partially inside the frustum.
    pub fn intersects_sphere(&self, center: [f64; 3], radius: f64) -> bool {
        for plane in &self.planes {
            let [nx, ny, nz, d] = *plane;
            let dist = nx * center[0] + ny * center[1] + nz * center[2] + d;
            if dist < -radius {
                return false;
            }
        }
        true
    }
}
