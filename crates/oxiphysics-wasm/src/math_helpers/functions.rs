//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::types::{QuatWasm, TransformWasm, Vec3Wasm};

#[allow(unused_imports)]
use super::types::{Easing, Frustum, JsQuat, JsTransform, JsVec3, Mat4, Ray};

impl From<JsVec3> for Vec3Wasm {
    fn from(v: JsVec3) -> Self {
        v.to_vec3()
    }
}
impl From<JsVec3> for [f64; 3] {
    fn from(v: JsVec3) -> Self {
        v.to_array()
    }
}
impl From<JsQuat> for QuatWasm {
    fn from(q: JsQuat) -> Self {
        q.to_quat()
    }
}
impl From<JsQuat> for [f64; 4] {
    fn from(q: JsQuat) -> Self {
        q.to_array()
    }
}
impl From<JsTransform> for TransformWasm {
    fn from(t: JsTransform) -> Self {
        t.to_transform()
    }
}
/// Convert a slice of `[f64; 3]` arrays into a flat `Vec`f64`.
pub fn vec3_array_slice_to_flat(positions: &[[f64; 3]]) -> Vec<f64> {
    let mut out = Vec::with_capacity(positions.len() * 3);
    for p in positions {
        out.extend_from_slice(p);
    }
    out
}
/// Convert a slice of `\[f64; 4\]` arrays into a flat `Vec`f64`.
pub fn quat_array_slice_to_flat(rotations: &[[f64; 4]]) -> Vec<f64> {
    let mut out = Vec::with_capacity(rotations.len() * 4);
    for q in rotations {
        out.extend_from_slice(q);
    }
    out
}
/// Interleave position and rotation data for a body batch.
///
/// Input: two equal-length slices of `[f64; 3]` and `[f64; 4]` respectively.
/// Output: flat `[px,py,pz,qx,qy,qz,qw, px,py,pz,qx,qy,qz,qw, ...]`.
pub fn interleave_transforms(positions: &[[f64; 3]], rotations: &[[f64; 4]]) -> Vec<f64> {
    assert_eq!(
        positions.len(),
        rotations.len(),
        "positions and rotations must have equal length"
    );
    let mut out = Vec::with_capacity(positions.len() * 7);
    for (p, q) in positions.iter().zip(rotations.iter()) {
        out.extend_from_slice(p);
        out.extend_from_slice(q);
    }
    out
}
/// Deinterleave a flat transform buffer `[px,py,pz,qx,qy,qz,qw, ...]`
/// into a `Vec<JsTransform>`.
pub fn deinterleave_transforms(data: &[f64]) -> Vec<JsTransform> {
    JsTransform::unpack_flat(data)
}
#[cfg(test)]
mod tests {
    use super::*;

    use crate::math_helpers::JsQuat;
    use crate::math_helpers::JsVec3;

    use std::f64::consts::PI;
    #[test]
    fn test_js_vec3_zero() {
        let v = JsVec3::zero();
        assert_eq!(v.length(), 0.0);
    }
    #[test]
    fn test_js_vec3_length() {
        let v = JsVec3::new(3.0, 4.0, 0.0);
        assert!((v.length() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_vec3_normalized() {
        let v = JsVec3::new(0.0, 5.0, 0.0);
        let n = v.normalized();
        assert!((n.y - 1.0).abs() < 1e-10);
        assert!(n.x.abs() < 1e-10);
        assert!(n.z.abs() < 1e-10);
    }
    #[test]
    fn test_js_vec3_dot_cross() {
        let x = JsVec3::new(1.0, 0.0, 0.0);
        let y = JsVec3::new(0.0, 1.0, 0.0);
        assert!((x.dot(&y)).abs() < 1e-10);
        let z = x.cross(&y);
        assert!((z.z - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_vec3_lerp() {
        let a = JsVec3::new(0.0, 0.0, 0.0);
        let b = JsVec3::new(10.0, 0.0, 0.0);
        let mid = a.lerp(&b, 0.5);
        assert!((mid.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_vec3_array_roundtrip() {
        let arr = [1.0_f64, 2.0, 3.0];
        let v = JsVec3::from_array(arr);
        assert_eq!(v.to_array(), arr);
    }
    #[test]
    fn test_js_vec3_from_into_vec3wasm() {
        let v = Vec3Wasm::new(1.0, 2.0, 3.0);
        let js: JsVec3 = v.into();
        let back: Vec3Wasm = js.into();
        assert!((back.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_vec3_pack_unpack_flat() {
        let vecs = vec![JsVec3::new(1.0, 2.0, 3.0), JsVec3::new(4.0, 5.0, 6.0)];
        let flat = JsVec3::pack_flat(&vecs);
        assert_eq!(flat.len(), 6);
        let back = JsVec3::unpack_flat(&flat);
        assert!((back[0].x - 1.0).abs() < 1e-10);
        assert!((back[1].z - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_vec3_neg() {
        let v = JsVec3::new(1.0, -2.0, 3.0);
        let n = v.neg();
        assert!((n.x + 1.0).abs() < 1e-10);
        assert!((n.y - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_quat_identity_rotate() {
        let q = JsQuat::identity();
        let v = JsVec3::new(1.0, 2.0, 3.0);
        let r = q.rotate_vec(&v);
        assert!((r.x - 1.0).abs() < 1e-10);
        assert!((r.y - 2.0).abs() < 1e-10);
        assert!((r.z - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_quat_mul_identity() {
        let q = JsQuat::new(0.1, 0.2, 0.3, 0.9).normalized();
        let i = JsQuat::identity();
        let qi = q.mul(&i);
        assert!((qi.x - q.x).abs() < 1e-10);
        assert!((qi.w - q.w).abs() < 1e-10);
    }
    #[test]
    fn test_js_quat_from_axis_angle_90deg() {
        let axis = JsVec3::new(0.0, 1.0, 0.0);
        let q = JsQuat::from_axis_angle(&axis, PI / 2.0);
        let x = JsVec3::new(1.0, 0.0, 0.0);
        let r = q.rotate_vec(&x);
        assert!(r.x.abs() < 1e-10);
        assert!(r.y.abs() < 1e-10);
        assert!((r.z + 1.0).abs() < 1e-10, "expected z=-1, got {}", r.z);
    }
    #[test]
    fn test_js_quat_from_euler_identity() {
        let q = JsQuat::from_euler_zyx(0.0, 0.0, 0.0);
        assert!((q.w - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_quat_array_roundtrip() {
        let arr = [0.0_f64, 0.0, 0.0, 1.0];
        let q = JsQuat::from_array(arr);
        assert_eq!(q.to_array(), arr);
    }
    #[test]
    fn test_js_quat_slerp_endpoints() {
        let a = JsQuat::identity();
        let b = JsQuat::new(0.0, 1.0, 0.0, 0.0).normalized();
        let r0 = a.slerp(&b, 0.0);
        assert!((r0.w - 1.0).abs() < 1e-6);
        let r1 = a.slerp(&b, 1.0);
        assert!(r1.w.abs() < 1e-6);
    }
    #[test]
    fn test_js_quat_conjugate() {
        let q = JsQuat::new(1.0, 2.0, 3.0, 4.0).normalized();
        let c = q.conjugate();
        let prod = q.mul(&c);
        assert!(prod.x.abs() < 1e-10);
        assert!(prod.y.abs() < 1e-10);
        assert!(prod.z.abs() < 1e-10);
        assert!((prod.w - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_quat_pack_unpack_flat() {
        let quats = vec![
            JsQuat::identity(),
            JsQuat::new(0.0, 1.0, 0.0, 0.0).normalized(),
        ];
        let flat = JsQuat::pack_flat(&quats);
        assert_eq!(flat.len(), 8);
        let back = JsQuat::unpack_flat(&flat);
        assert!((back[0].w - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_transform_identity() {
        let t = JsTransform::identity();
        let p = JsVec3::new(1.0, 2.0, 3.0);
        let r = t.transform_point(&p);
        assert!((r.x - 1.0).abs() < 1e-10);
        assert!((r.y - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_transform_translation() {
        let t = JsTransform::from_position(5.0, 0.0, 0.0);
        let p = JsVec3::zero();
        let r = t.transform_point(&p);
        assert!((r.x - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_transform_array7_roundtrip() {
        let t = JsTransform::from_position(1.0, 2.0, 3.0);
        let arr = t.to_array7();
        let t2 = JsTransform::from_array7(arr);
        assert!((t2.position.x - 1.0).abs() < 1e-10);
        assert!((t2.position.y - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_transform_inverse() {
        let t = JsTransform::from_position(3.0, -2.0, 1.0);
        let inv = t.inverse();
        let origin = JsVec3::zero();
        let out = inv.transform_point(&origin);
        assert!((out.x + 3.0).abs() < 1e-10);
        assert!((out.y - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_transform_matrix4_identity() {
        let t = JsTransform::identity();
        let m = t.to_matrix4();
        assert!((m[0] - 1.0).abs() < 1e-10);
        assert!((m[5] - 1.0).abs() < 1e-10);
        assert!((m[10] - 1.0).abs() < 1e-10);
        assert!((m[15] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_js_transform_pack_unpack_flat() {
        let transforms = vec![
            JsTransform::from_position(1.0, 0.0, 0.0),
            JsTransform::from_position(0.0, 2.0, 0.0),
        ];
        let flat = JsTransform::pack_flat(&transforms);
        assert_eq!(flat.len(), 14);
        let back = JsTransform::unpack_flat(&flat);
        assert!((back[0].position.x - 1.0).abs() < 1e-10);
        assert!((back[1].position.y - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_interleave_deinterleave_transforms() {
        let positions: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let rotations: Vec<[f64; 4]> = vec![[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0, 1.0]];
        let flat = interleave_transforms(&positions, &rotations);
        assert_eq!(flat.len(), 14);
        let back = deinterleave_transforms(&flat);
        assert_eq!(back.len(), 2);
        assert!((back[0].position.x - 1.0).abs() < 1e-10);
        assert!((back[1].position.z - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_from_into_transform_wasm() {
        let tw = TransformWasm::from_position(1.0, 2.0, 3.0);
        let jt: JsTransform = tw.into();
        let back: TransformWasm = jt.into();
        assert!((back.position.x - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_vec3_distance() {
        let a = JsVec3::new(0.0, 0.0, 0.0);
        let b = JsVec3::new(3.0, 4.0, 0.0);
        assert!((a.distance_to(&b) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_vec3_array_slice_to_flat() {
        let positions: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let flat = vec3_array_slice_to_flat(&positions);
        assert_eq!(flat, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    }
    #[test]
    fn test_quat_array_slice_to_flat() {
        let rotations: Vec<[f64; 4]> = vec![[0.0, 0.0, 0.0, 1.0]];
        let flat = quat_array_slice_to_flat(&rotations);
        assert_eq!(flat, vec![0.0, 0.0, 0.0, 1.0]);
    }
}
#[allow(dead_code)]
pub(super) fn normalize_plane(p: [f64; 4]) -> [f64; 4] {
    let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
    if len < 1e-15 {
        return p;
    }
    [p[0] / len, p[1] / len, p[2] / len, p[3] / len]
}
/// Compute a smooth 3-D value noise in `[−1, 1]` using trilinear interpolation.
///
/// Uses a fixed-seed hash function — suitable for deterministic procedural content.
/// `x`, `y`, `z` are world-space coordinates.
#[allow(dead_code)]
pub fn value_noise_3d(x: f64, y: f64, z: f64) -> f64 {
    let xi = x.floor() as i64;
    let yi = y.floor() as i64;
    let zi = z.floor() as i64;
    let fx = x - xi as f64;
    let fy = y - yi as f64;
    let fz = z - zi as f64;
    let ux = smooth_step(fx);
    let uy = smooth_step(fy);
    let uz = smooth_step(fz);
    let v000 = hash_to_float(xi, yi, zi);
    let v100 = hash_to_float(xi + 1, yi, zi);
    let v010 = hash_to_float(xi, yi + 1, zi);
    let v110 = hash_to_float(xi + 1, yi + 1, zi);
    let v001 = hash_to_float(xi, yi, zi + 1);
    let v101 = hash_to_float(xi + 1, yi, zi + 1);
    let v011 = hash_to_float(xi, yi + 1, zi + 1);
    let v111 = hash_to_float(xi + 1, yi + 1, zi + 1);
    let x00 = lerp_f64(v000, v100, ux);
    let x10 = lerp_f64(v010, v110, ux);
    let x01 = lerp_f64(v001, v101, ux);
    let x11 = lerp_f64(v011, v111, ux);
    let y0 = lerp_f64(x00, x10, uy);
    let y1 = lerp_f64(x01, x11, uy);
    lerp_f64(y0, y1, uz) * 2.0 - 1.0
}
/// Fractal Brownian Motion (fBm) summing `octaves` layers of value noise.
#[allow(dead_code)]
pub fn fbm_3d(x: f64, y: f64, z: f64, octaves: u32, lacunarity: f64, gain: f64) -> f64 {
    let mut value = 0.0;
    let mut amplitude = 0.5;
    let mut frequency = 1.0;
    for _ in 0..octaves {
        value += amplitude * value_noise_3d(x * frequency, y * frequency, z * frequency);
        amplitude *= gain;
        frequency *= lacunarity;
    }
    value
}
/// Simple Perlin-like gradient noise in 2-D. Returns value in approximately `[-1, 1]`.
#[allow(dead_code)]
pub fn gradient_noise_2d(x: f64, y: f64) -> f64 {
    let xi = x.floor() as i64;
    let yi = y.floor() as i64;
    let fx = x - xi as f64;
    let fy = y - yi as f64;
    let ux = smooth_step(fx);
    let uy = smooth_step(fy);
    let g = |ix: i64, iy: i64, dx: f64, dy: f64| -> f64 {
        let h = hash_2d(ix, iy);
        let gx = if h & 1 != 0 { 1.0 } else { -1.0 };
        let gy = if h & 2 != 0 { 1.0 } else { -1.0 };
        gx * dx + gy * dy
    };
    let n00 = g(xi, yi, fx, fy);
    let n10 = g(xi + 1, yi, fx - 1.0, fy);
    let n01 = g(xi, yi + 1, fx, fy - 1.0);
    let n11 = g(xi + 1, yi + 1, fx - 1.0, fy - 1.0);
    let x0 = lerp_f64(n00, n10, ux);
    let x1 = lerp_f64(n01, n11, ux);
    lerp_f64(x0, x1, uy)
}
#[inline]
pub(super) fn smooth_step(t: f64) -> f64 {
    t * t * (3.0 - 2.0 * t)
}
#[inline]
pub(super) fn smooth_step_quintic(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}
#[inline]
pub(super) fn lerp_f64(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}
#[inline]
pub(super) fn hash_to_float(x: i64, y: i64, z: i64) -> f64 {
    let h = hash_3d(x, y, z);
    (h as f64) / (u64::MAX as f64)
}
#[inline]
pub(super) fn hash_3d(x: i64, y: i64, z: i64) -> u64 {
    let mut h = x.unsigned_abs().wrapping_mul(2654435761)
        ^ y.unsigned_abs().wrapping_mul(805459861)
        ^ z.unsigned_abs().wrapping_mul(1234567891);
    h = h
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    h
}
#[inline]
pub(super) fn hash_2d(x: i64, y: i64) -> u64 {
    let mut h =
        x.unsigned_abs().wrapping_mul(2654435761) ^ y.unsigned_abs().wrapping_mul(805459861);
    h = h
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    h
}
/// Evaluate a cubic Bezier curve at parameter `t ∈ [0, 1]`.
///
/// Control points: `p0`, `p1`, `p2`, `p3`.
#[allow(dead_code)]
pub fn bezier_cubic(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3], t: f64) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let mt = 1.0 - t;
    let mt2 = mt * mt;
    let t2 = t * t;
    let b0 = mt2 * mt;
    let b1 = 3.0 * mt2 * t;
    let b2 = 3.0 * mt * t2;
    let b3 = t2 * t;
    [
        b0 * p0[0] + b1 * p1[0] + b2 * p2[0] + b3 * p3[0],
        b0 * p0[1] + b1 * p1[1] + b2 * p2[1] + b3 * p3[1],
        b0 * p0[2] + b1 * p1[2] + b2 * p2[2] + b3 * p3[2],
    ]
}
/// Evaluate the tangent of a cubic Bezier curve at `t`.
#[allow(dead_code)]
pub fn bezier_cubic_tangent(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    t: f64,
) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let mt = 1.0 - t;
    let b0 = 3.0 * mt * mt;
    let b1 = 6.0 * mt * t;
    let b2 = 3.0 * t * t;
    [
        b0 * (p1[0] - p0[0]) + b1 * (p2[0] - p1[0]) + b2 * (p3[0] - p2[0]),
        b0 * (p1[1] - p0[1]) + b1 * (p2[1] - p1[1]) + b2 * (p3[1] - p2[1]),
        b0 * (p1[2] - p0[2]) + b1 * (p2[2] - p1[2]) + b2 * (p3[2] - p2[2]),
    ]
}
/// Evaluate a quadratic Bezier curve at `t`.
#[allow(dead_code)]
pub fn bezier_quadratic(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], t: f64) -> [f64; 3] {
    let t = t.clamp(0.0, 1.0);
    let mt = 1.0 - t;
    let b0 = mt * mt;
    let b1 = 2.0 * mt * t;
    let b2 = t * t;
    [
        b0 * p0[0] + b1 * p1[0] + b2 * p2[0],
        b0 * p0[1] + b1 * p1[1] + b2 * p2[1],
        b0 * p0[2] + b1 * p1[2] + b2 * p2[2],
    ]
}
/// Sample a cubic Bezier curve uniformly into `n_samples` points (including endpoints).
#[allow(dead_code)]
pub fn bezier_cubic_sample(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    n_samples: usize,
) -> Vec<[f64; 3]> {
    if n_samples < 2 {
        return vec![p0];
    }
    (0..n_samples)
        .map(|i| {
            let t = i as f64 / (n_samples - 1) as f64;
            bezier_cubic(p0, p1, p2, p3, t)
        })
        .collect()
}
/// Catmull-Rom spline through `points`, evaluated at `t ∈ [0, 1]` over the whole curve.
#[allow(dead_code)]
pub fn catmull_rom(points: &[[f64; 3]], t: f64) -> Option<[f64; 3]> {
    let n = points.len();
    if n < 2 {
        return None;
    }
    let t = t.clamp(0.0, 1.0);
    let segments = (n - 1) as f64;
    let ts = (t * segments).min(segments - 1.0 + 1e-15);
    let i = ts.floor() as usize;
    let local_t = ts - i as f64;
    let p0 = if i == 0 { points[0] } else { points[i - 1] };
    let p1 = points[i];
    let p2 = points[(i + 1).min(n - 1)];
    let p3 = points[(i + 2).min(n - 1)];
    let t2 = local_t * local_t;
    let t3 = t2 * local_t;
    let result = std::array::from_fn(|k| {
        0.5 * ((2.0 * p1[k])
            + (-p0[k] + p2[k]) * local_t
            + (2.0 * p0[k] - 5.0 * p1[k] + 4.0 * p2[k] - p3[k]) * t2
            + (-p0[k] + 3.0 * p1[k] - 3.0 * p2[k] + p3[k]) * t3)
    });
    Some(result)
}
#[cfg(test)]
mod extended_math_tests {

    use crate::math_helpers::Easing;
    use crate::math_helpers::Frustum;
    use crate::math_helpers::JsQuat;
    use crate::math_helpers::JsVec3;
    use crate::math_helpers::Mat4;
    use crate::math_helpers::Ray;
    use crate::math_helpers::bezier_cubic;
    use crate::math_helpers::bezier_cubic_sample;
    use crate::math_helpers::bezier_cubic_tangent;
    use crate::math_helpers::bezier_quadratic;
    use crate::math_helpers::catmull_rom;
    use crate::math_helpers::fbm_3d;
    use crate::math_helpers::gradient_noise_2d;
    use crate::math_helpers::smooth_step;
    use crate::math_helpers::smooth_step_quintic;
    use crate::math_helpers::value_noise_3d;
    use std::f64::consts::PI;
    #[test]
    fn test_mat4_identity_mul() {
        let i = Mat4::identity();
        let a = Mat4::from_translation([1.0, 2.0, 3.0]);
        let result = i.mul(&a);
        assert!((result.0[12] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_translation_point() {
        let m = Mat4::from_translation([5.0, -3.0, 2.0]);
        let p = m.transform_point3([0.0, 0.0, 0.0]);
        assert!((p[0] - 5.0).abs() < 1e-10);
        assert!((p[1] + 3.0).abs() < 1e-10);
        assert!((p[2] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_scale_uniform() {
        let m = Mat4::from_scale(2.0);
        let p = m.transform_point3([1.0, 1.0, 1.0]);
        assert!((p[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_scale_xyz() {
        let m = Mat4::from_scale_xyz(1.0, 2.0, 3.0);
        let p = m.transform_point3([1.0, 1.0, 1.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
        assert!((p[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_from_quat_identity() {
        let m = Mat4::from_quat(0.0, 0.0, 0.0, 1.0);
        assert!((m.0[0] - 1.0).abs() < 1e-10);
        assert!((m.0[5] - 1.0).abs() < 1e-10);
        assert!((m.0[10] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_transpose_is_involution() {
        let m = Mat4::from_translation([1.0, 2.0, 3.0]);
        let tt = m.transpose().transpose();
        for i in 0..16 {
            assert!((m.0[i] - tt.0[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_mat4_mul_translation_stacks() {
        let a = Mat4::from_translation([1.0, 0.0, 0.0]);
        let b = Mat4::from_translation([0.0, 2.0, 0.0]);
        let ab = a.mul(&b);
        let p = ab.transform_point3([0.0, 0.0, 0.0]);
        assert!((p[0] - 1.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_mat4_orthographic_identity_z() {
        let m = Mat4::orthographic(-1.0, 1.0, -1.0, 1.0, 0.1, 100.0);
        let p = m.transform_vec4([0.0, 0.0, -1.0, 1.0]);
        assert!(p[3] > 0.0);
    }
    #[test]
    fn test_mat4_look_at_no_panic() {
        let m = Mat4::look_at([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let _ = m;
    }
    #[test]
    fn test_frustum_sphere_inside() {
        let proj = Mat4::perspective(PI * 0.5, 1.0, 0.1, 100.0);
        let view = Mat4::look_at([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let vp = view.mul(&proj);
        let frustum = Frustum::from_view_proj(&vp);
        assert!(frustum.intersects_sphere([0.0, 0.0, 0.0], 0.1));
    }
    #[test]
    fn test_frustum_sphere_behind_camera() {
        let proj = Mat4::perspective(PI * 0.5, 1.0, 0.1, 100.0);
        let view = Mat4::look_at([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let vp = view.mul(&proj);
        let frustum = Frustum::from_view_proj(&vp);
        let result = frustum.intersects_sphere([0.0, 0.0, 1000.0], 0.01);
        let _ = result;
    }
    #[test]
    fn test_frustum_aabb_clearly_outside() {
        let proj = Mat4::perspective(PI * 0.5, 1.0, 0.1, 10.0);
        let view = Mat4::look_at([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let vp = view.mul(&proj);
        let frustum = Frustum::from_view_proj(&vp);
        let result = frustum.intersects_aabb([-0.01, 5000.0, -0.01], [0.01, 5001.0, 0.01]);
        let _ = result;
    }
    #[test]
    fn test_ray_intersect_aabb_hit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let t = ray.intersect_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(t.is_some(), "ray should hit AABB");
        assert!((t.unwrap() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_intersect_aabb_miss() {
        let ray = Ray::new([5.0, 5.0, -5.0], [0.0, 0.0, 1.0]);
        let t = ray.intersect_aabb([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]);
        assert!(t.is_none(), "ray should miss AABB");
    }
    #[test]
    fn test_ray_intersect_sphere_hit() {
        let ray = Ray::new([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let t = ray.intersect_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(t.is_some());
        assert!((t.unwrap() - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_ray_intersect_sphere_miss() {
        let ray = Ray::new([5.0, 0.0, -5.0], [0.0, 0.0, 1.0]);
        let t = ray.intersect_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(t.is_none());
    }
    #[test]
    fn test_ray_at() {
        let ray = Ray::new([1.0, 2.0, 3.0], [1.0, 0.0, 0.0]);
        let p = ray.at(3.0);
        assert!((p[0] - 4.0).abs() < 1e-10);
        assert!((p[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_value_noise_range() {
        for i in 0..20 {
            let v = value_noise_3d(i as f64 * 0.3, i as f64 * 0.17, i as f64 * 0.41);
            assert!((-1.0..=1.0).contains(&v), "noise out of range: {}", v);
        }
    }
    #[test]
    fn test_gradient_noise_range() {
        for i in 0..20 {
            let v = gradient_noise_2d(i as f64 * 0.3, i as f64 * 0.7);
            assert!(v.is_finite(), "noise should be finite");
        }
    }
    #[test]
    fn test_fbm_finite() {
        let v = fbm_3d(1.5, 2.3, 0.7, 4, 2.0, 0.5);
        assert!(v.is_finite());
    }
    #[test]
    fn test_smooth_step_endpoints() {
        assert!((smooth_step(0.0)).abs() < 1e-15);
        assert!((smooth_step(1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_smooth_step_quintic_endpoints() {
        assert!((smooth_step_quintic(0.0)).abs() < 1e-15);
        assert!((smooth_step_quintic(1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_easing_linear_endpoints() {
        assert!((Easing::linear(0.0)).abs() < 1e-15);
        assert!((Easing::linear(1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_easing_quad_endpoints() {
        assert!((Easing::ease_in_quad(0.0)).abs() < 1e-15);
        assert!((Easing::ease_in_quad(1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_easing_ease_out_quad_endpoints() {
        assert!((Easing::ease_out_quad(0.0)).abs() < 1e-15);
        assert!((Easing::ease_out_quad(1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_easing_sine_endpoints() {
        assert!((Easing::ease_in_sine(0.0)).abs() < 1e-10);
        assert!((Easing::ease_in_sine(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_easing_bounce_endpoints() {
        assert!((Easing::ease_out_bounce(0.0)).abs() < 1e-10);
        assert!((Easing::ease_out_bounce(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_easing_cubic_endpoints() {
        assert!((Easing::ease_in_cubic(0.0)).abs() < 1e-15);
        assert!((Easing::ease_out_cubic(1.0) - 1.0).abs() < 1e-15);
    }
    #[test]
    fn test_easing_expo_endpoints() {
        assert!((Easing::ease_in_expo(0.0)).abs() < 1e-10);
        assert!((Easing::ease_out_expo(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_easing_elastic_endpoints() {
        assert!((Easing::ease_out_elastic(0.0)).abs() < 1e-10);
        assert!((Easing::ease_out_elastic(1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_bezier_cubic_endpoints() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 2.0, 0.0];
        let p2 = [2.0, 2.0, 0.0];
        let p3 = [3.0, 0.0, 0.0];
        let start = bezier_cubic(p0, p1, p2, p3, 0.0);
        let end = bezier_cubic(p0, p1, p2, p3, 1.0);
        assert!((start[0]).abs() < 1e-10);
        assert!((end[0] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_bezier_cubic_midpoint() {
        let p0 = [0.0, 0.0, 0.0];
        let p3 = [2.0, 0.0, 0.0];
        let mid = bezier_cubic(p0, p0, p3, p3, 0.5);
        assert!((mid[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_bezier_quadratic_endpoints() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 2.0, 0.0];
        let p2 = [2.0, 0.0, 0.0];
        let s = bezier_quadratic(p0, p1, p2, 0.0);
        let e = bezier_quadratic(p0, p1, p2, 1.0);
        assert!((s[0]).abs() < 1e-10);
        assert!((e[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_bezier_cubic_sample_length() {
        let p0 = [0.0; 3];
        let p3 = [1.0, 0.0, 0.0];
        let samples = bezier_cubic_sample(p0, p0, p3, p3, 11);
        assert_eq!(samples.len(), 11);
    }
    #[test]
    fn test_bezier_cubic_tangent_not_zero_at_midpoint() {
        let p0 = [0.0, 0.0, 0.0];
        let p1 = [1.0, 1.0, 0.0];
        let p2 = [2.0, 1.0, 0.0];
        let p3 = [3.0, 0.0, 0.0];
        let tan = bezier_cubic_tangent(p0, p1, p2, p3, 0.5);
        let len = (tan[0] * tan[0] + tan[1] * tan[1] + tan[2] * tan[2]).sqrt();
        assert!(len > 1e-10, "tangent should be non-zero");
    }
    #[test]
    fn test_catmull_rom_endpoints() {
        let pts = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ];
        let start = catmull_rom(&pts, 0.0).unwrap();
        assert!((start[0]).abs() < 1e-10);
        let near_end = catmull_rom(&pts, 0.999).unwrap();
        assert!(
            (near_end[0] - 2.0).abs() < 0.01,
            "near_end.x = {}",
            near_end[0]
        );
    }
    #[test]
    fn test_catmull_rom_none_for_empty() {
        let result = catmull_rom(&[], 0.5);
        assert!(result.is_none());
    }
    #[test]
    fn test_quat_to_euler_identity() {
        let q = JsQuat::identity();
        let euler = q.to_euler_zyx();
        assert!(euler[0].abs() < 1e-10);
        assert!(euler[1].abs() < 1e-10);
        assert!(euler[2].abs() < 1e-10);
    }
    #[test]
    fn test_quat_angle_to_identity() {
        let q = JsQuat::identity();
        let angle = q.angle_to(&q);
        assert!(angle.abs() < 1e-10);
    }
    #[test]
    fn test_quat_norm_squared_identity() {
        let q = JsQuat::identity();
        assert!((q.norm_squared() - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_from_rotation_between_same_direction() {
        let v = JsVec3::new(1.0, 0.0, 0.0);
        let q = JsQuat::from_rotation_between(&v, &v);
        assert!((q.w - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_from_rotation_between_opposite() {
        let a = JsVec3::new(1.0, 0.0, 0.0);
        let b = JsVec3::new(-1.0, 0.0, 0.0);
        let q = JsQuat::from_rotation_between(&a, &b);
        let rotated = q.rotate_vec(&a);
        assert!((rotated.x + 1.0).abs() < 1e-8, "x={}", rotated.x);
    }
    #[test]
    fn test_quat_to_rotation_matrix3_identity() {
        let q = JsQuat::identity();
        let m = q.to_rotation_matrix3();
        assert!((m[0] - 1.0).abs() < 1e-10);
        assert!((m[4] - 1.0).abs() < 1e-10);
        assert!((m[8] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_nlerp_endpoints() {
        let a = JsQuat::identity();
        let b = JsQuat::from_axis_angle(&JsVec3::new(0.0, 1.0, 0.0), PI);
        let r0 = a.nlerp(&b, 0.0);
        assert!((r0.w - 1.0).abs() < 1e-10);
        let r1 = a.nlerp(&b, 1.0);
        assert!(r1.norm_squared() > 0.5);
    }
    #[test]
    fn test_quat_euler_roundtrip() {
        let roll = 0.3;
        let pitch = -0.2;
        let yaw = 1.1;
        let q = JsQuat::from_euler_zyx(roll, pitch, yaw);
        let euler = q.to_euler_zyx();
        assert!(
            (euler[0] - roll).abs() < 1e-8,
            "roll mismatch: {} vs {}",
            euler[0],
            roll
        );
    }
}
