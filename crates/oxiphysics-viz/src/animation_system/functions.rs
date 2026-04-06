//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{AnimLod, Easing, IkJoint, SkinWeight};

/// 3-component vector.
#[allow(dead_code)]
pub type Vec3 = [f64; 3];
/// 4-component quaternion stored as `[x, y, z, w]`.
#[allow(dead_code)]
pub type Quat = [f64; 4];
/// 4×4 column-major matrix (indices: col * 4 + row).
#[allow(dead_code)]
pub type Mat4 = [f64; 16];
/// Add two Vec3s.
#[allow(dead_code)]
pub fn vec3_add(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Subtract two Vec3s.
#[allow(dead_code)]
pub fn vec3_sub(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scale a Vec3.
#[allow(dead_code)]
pub fn vec3_scale(a: Vec3, s: f64) -> Vec3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
/// Dot product.
#[allow(dead_code)]
pub fn vec3_dot(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product.
#[allow(dead_code)]
pub fn vec3_cross(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Euclidean length.
#[allow(dead_code)]
pub fn vec3_len(a: Vec3) -> f64 {
    vec3_dot(a, a).sqrt()
}
/// Normalise (returns zero-vector if degenerate).
#[allow(dead_code)]
pub fn vec3_normalize(a: Vec3) -> Vec3 {
    let l = vec3_len(a);
    if l < 1e-300 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / l)
    }
}
/// Linear interpolation of Vec3.
#[allow(dead_code)]
pub fn vec3_lerp(a: Vec3, b: Vec3, t: f64) -> Vec3 {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}
/// Quaternion normalise.
#[allow(dead_code)]
pub fn quat_normalize(q: Quat) -> Quat {
    let len2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if len2 < 1e-300 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    let inv = 1.0 / len2.sqrt();
    [q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv]
}
/// Quaternion multiplication `p * q`.
#[allow(dead_code)]
pub fn quat_mul(p: Quat, q: Quat) -> Quat {
    [
        p[3] * q[0] + p[0] * q[3] + p[1] * q[2] - p[2] * q[1],
        p[3] * q[1] - p[0] * q[2] + p[1] * q[3] + p[2] * q[0],
        p[3] * q[2] + p[0] * q[1] - p[1] * q[0] + p[2] * q[3],
        p[3] * q[3] - p[0] * q[0] - p[1] * q[1] - p[2] * q[2],
    ]
}
/// Quaternion SLERP.
#[allow(dead_code)]
pub fn quat_slerp(a: Quat, b: Quat, t: f64) -> Quat {
    let mut cos_half = a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
    let b_adj = if cos_half < 0.0 {
        cos_half = -cos_half;
        [-b[0], -b[1], -b[2], -b[3]]
    } else {
        b
    };
    let (scale_a, scale_b) = if 1.0 - cos_half > 1e-6 {
        let half = cos_half.acos();
        let sin_half = half.sin();
        (
            ((1.0 - t) * half).sin() / sin_half,
            (t * half).sin() / sin_half,
        )
    } else {
        (1.0 - t, t)
    };
    quat_normalize([
        scale_a * a[0] + scale_b * b_adj[0],
        scale_a * a[1] + scale_b * b_adj[1],
        scale_a * a[2] + scale_b * b_adj[2],
        scale_a * a[3] + scale_b * b_adj[3],
    ])
}
/// Convert a unit quaternion to a 4×4 rotation matrix (column-major).
#[allow(dead_code)]
pub fn quat_to_mat4(q: Quat) -> Mat4 {
    let [x, y, z, w] = q;
    let x2 = x * x;
    let y2 = y * y;
    let z2 = z * z;
    let xy = x * y;
    let xz = x * z;
    let yz = y * z;
    let wx = w * x;
    let wy = w * y;
    let wz = w * z;
    [
        1.0 - 2.0 * (y2 + z2),
        2.0 * (xy + wz),
        2.0 * (xz - wy),
        0.0,
        2.0 * (xy - wz),
        1.0 - 2.0 * (x2 + z2),
        2.0 * (yz + wx),
        0.0,
        2.0 * (xz + wy),
        2.0 * (yz - wx),
        1.0 - 2.0 * (x2 + y2),
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
    ]
}
/// Compose a TRS matrix (translation * rotation * scale).
#[allow(dead_code)]
pub fn trs_matrix(translation: Vec3, rotation: Quat, scale: Vec3) -> Mat4 {
    let r = quat_to_mat4(rotation);
    [
        r[0] * scale[0],
        r[1] * scale[0],
        r[2] * scale[0],
        0.0,
        r[4] * scale[1],
        r[5] * scale[1],
        r[6] * scale[1],
        0.0,
        r[8] * scale[2],
        r[9] * scale[2],
        r[10] * scale[2],
        0.0,
        translation[0],
        translation[1],
        translation[2],
        1.0,
    ]
}
/// Multiply two 4×4 column-major matrices.
#[allow(dead_code)]
pub fn mat4_mul(a: Mat4, b: Mat4) -> Mat4 {
    let mut out = [0.0_f64; 16];
    for col in 0..4 {
        for row in 0..4 {
            let mut s = 0.0_f64;
            for k in 0..4 {
                s += a[k * 4 + row] * b[col * 4 + k];
            }
            out[col * 4 + row] = s;
        }
    }
    out
}
/// Transform a point by a Mat4.
#[allow(dead_code)]
pub fn mat4_transform_point(m: Mat4, p: Vec3) -> Vec3 {
    let w = m[3] * p[0] + m[7] * p[1] + m[11] * p[2] + m[15];
    let inv_w = if w.abs() < 1e-300 { 1.0 } else { 1.0 / w };
    [
        (m[0] * p[0] + m[4] * p[1] + m[8] * p[2] + m[12]) * inv_w,
        (m[1] * p[0] + m[5] * p[1] + m[9] * p[2] + m[13]) * inv_w,
        (m[2] * p[0] + m[6] * p[1] + m[10] * p[2] + m[14]) * inv_w,
    ]
}
/// Apply the chosen easing function to normalised time `t ∈ [0, 1]`.
#[allow(dead_code)]
pub fn apply_easing(easing: Easing, t: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    match easing {
        Easing::Linear => t,
        Easing::EaseInQuad => t * t,
        Easing::EaseOutQuad => t * (2.0 - t),
        Easing::EaseInOutQuad => {
            if t < 0.5 {
                2.0 * t * t
            } else {
                -1.0 + (4.0 - 2.0 * t) * t
            }
        }
        Easing::EaseInCubic => t * t * t,
        Easing::EaseOutCubic => {
            let t1 = t - 1.0;
            t1 * t1 * t1 + 1.0
        }
        Easing::EaseInOutCubic => {
            if t < 0.5 {
                4.0 * t * t * t
            } else {
                let t1 = 2.0 * t - 2.0;
                0.5 * t1 * t1 * t1 + 1.0
            }
        }
        Easing::EaseInSine => 1.0 - (t * PI * 0.5).cos(),
        Easing::EaseOutSine => (t * PI * 0.5).sin(),
        Easing::EaseInOutSine => 0.5 * (1.0 - (PI * t).cos()),
        Easing::EaseInExpo => {
            if t <= 0.0 {
                0.0
            } else {
                (2.0_f64).powf(10.0 * t - 10.0)
            }
        }
        Easing::EaseOutExpo => {
            if t >= 1.0 {
                1.0
            } else {
                1.0 - (2.0_f64).powf(-10.0 * t)
            }
        }
        Easing::EaseOutElastic => {
            if t == 0.0 {
                0.0
            } else if t == 1.0 {
                1.0
            } else {
                let c4 = TAU / 3.0;
                (2.0_f64).powf(-10.0 * t) * ((t * 10.0 - 0.75) * c4).sin() + 1.0
            }
        }
        Easing::EaseOutBounce => ease_out_bounce(t),
    }
}
/// Bounce ease-out helper.
#[allow(dead_code)]
pub(super) fn ease_out_bounce(t: f64) -> f64 {
    pub(super) const N: f64 = 7.5625;
    pub(super) const D: f64 = 2.75;
    if t < 1.0 / D {
        N * t * t
    } else if t < 2.0 / D {
        let t2 = t - 1.5 / D;
        N * t2 * t2 + 0.75
    } else if t < 2.5 / D {
        let t2 = t - 2.25 / D;
        N * t2 * t2 + 0.9375
    } else {
        let t2 = t - 2.625 / D;
        N * t2 * t2 + 0.984375
    }
}
/// Evaluate a cubic Bézier curve at parameter `t ∈ [0, 1]`.
/// Control points: `p0` (start), `p1` (ctrl1), `p2` (ctrl2), `p3` (end).
#[allow(dead_code)]
pub fn cubic_bezier(p0: f64, p1: f64, p2: f64, p3: f64, t: f64) -> f64 {
    let mt = 1.0 - t;
    mt * mt * mt * p0 + 3.0 * mt * mt * t * p1 + 3.0 * mt * t * t * p2 + t * t * t * p3
}
/// Evaluate a cubic Hermite spline at `t ∈ [0, 1]`.
/// `p0`, `p1` are endpoint values; `m0`, `m1` are tangents.
#[allow(dead_code)]
pub fn cubic_hermite(p0: f64, m0: f64, p1: f64, m1: f64, t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    (2.0 * t3 - 3.0 * t2 + 1.0) * p0
        + (t3 - 2.0 * t2 + t) * m0
        + (-2.0 * t3 + 3.0 * t2) * p1
        + (t3 - t2) * m1
}
/// Apply skinning to a vertex position given per-bone skinning matrices.
#[allow(dead_code)]
pub fn skin_vertex(rest_pos: Vec3, weights: &[SkinWeight], skin_matrices: &[Mat4]) -> Vec3 {
    let mut result = [0.0_f64; 3];
    for sw in weights {
        if sw.bone_index < skin_matrices.len() {
            let p = mat4_transform_point(skin_matrices[sw.bone_index], rest_pos);
            result[0] += p[0] * sw.weight;
            result[1] += p[1] * sw.weight;
            result[2] += p[2] * sw.weight;
        }
    }
    result
}
/// Linear interpolation of two Mat4 (element-wise).
#[allow(dead_code)]
pub(super) fn lerp_mat4(a: Mat4, b: Mat4, t: f64) -> Mat4 {
    let mut out = [0.0_f64; 16];
    for i in 0..16 {
        out[i] = a[i] + (b[i] - a[i]) * t;
    }
    out
}
/// Identity 4×4 matrix.
#[allow(dead_code)]
pub(super) fn identity_mat4() -> Mat4 {
    [
        1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0,
    ]
}
/// CCD (Cyclic Coordinate Descent) IK solver.
/// Returns the updated joint positions after `max_iter` iterations.
#[allow(dead_code)]
pub fn ccd_ik(
    joints: &[IkJoint],
    bone_lengths: &[f64],
    target: Vec3,
    max_iter: usize,
    tolerance: f64,
) -> Vec<Vec3> {
    let n = joints.len();
    if n == 0 {
        return Vec::new();
    }
    let mut positions: Vec<Vec3> = joints.iter().map(|j| j.position).collect();
    for _iter in 0..max_iter {
        let end = positions[n - 1];
        if vec3_len(vec3_sub(end, target)) < tolerance {
            break;
        }
        for i in (0..n - 1).rev() {
            let dir_to_end = vec3_normalize(vec3_sub(positions[n - 1], positions[i]));
            let dir_to_target = vec3_normalize(vec3_sub(target, positions[i]));
            let cos_a = vec3_dot(dir_to_end, dir_to_target).clamp(-1.0, 1.0);
            let angle = cos_a.acos();
            if angle.abs() < 1e-9 {
                continue;
            }
            let axis = vec3_normalize(vec3_cross(dir_to_end, dir_to_target));
            let clamped = angle.min(joints[i].angle_limit);
            let s = (clamped * 0.5).sin();
            let c = (clamped * 0.5).cos();
            let q: Quat = [axis[0] * s, axis[1] * s, axis[2] * s, c];
            for k in i + 1..n {
                let rel = vec3_sub(positions[k], positions[i]);
                let qv: Quat = [rel[0], rel[1], rel[2], 0.0];
                let q_conj: Quat = [-q[0], -q[1], -q[2], q[3]];
                let rotated = quat_mul(quat_mul(q, qv), q_conj);
                positions[k] = vec3_add(positions[i], [rotated[0], rotated[1], rotated[2]]);
            }
        }
        for i in 0..n - 1 {
            let dir = vec3_normalize(vec3_sub(positions[i + 1], positions[i]));
            let bl = if i < bone_lengths.len() {
                bone_lengths[i]
            } else {
                1.0
            };
            positions[i + 1] = vec3_add(positions[i], vec3_scale(dir, bl));
        }
    }
    positions
}
/// Choose animation LOD based on distance from camera.
#[allow(dead_code)]
pub fn lod_for_distance(distance: f64) -> AnimLod {
    if distance < 10.0 {
        AnimLod::Full
    } else if distance < 30.0 {
        AnimLod::Medium
    } else if distance < 80.0 {
        AnimLod::Low
    } else {
        AnimLod::Disabled
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::animation_system::AnimLodManager;
    use crate::animation_system::AnimState;
    use crate::animation_system::AnimStateMachine;
    use crate::animation_system::AnimationClip;
    use crate::animation_system::Animator;
    use crate::animation_system::BlendNode;
    use crate::animation_system::Bone;
    use crate::animation_system::MocapFrame;
    use crate::animation_system::MocapSequence;
    use crate::animation_system::NodeAnimation;
    use crate::animation_system::QuatCurve;
    use crate::animation_system::QuatKey;
    use crate::animation_system::ScalarCurve;
    use crate::animation_system::ScalarKey;
    use crate::animation_system::Skeleton;
    use crate::animation_system::Transition;
    use crate::animation_system::Vec3Curve;
    use crate::animation_system::Vec3Key;
    use std::collections::HashMap;
    #[test]
    fn test_vec3_add() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let c = vec3_add(a, b);
        assert!((c[0] - 5.0).abs() < 1e-12);
        assert!((c[1] - 7.0).abs() < 1e-12);
        assert!((c[2] - 9.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_normalize_unit() {
        let a = [3.0_f64, 0.0, 0.0];
        let n = vec3_normalize(a);
        assert!((n[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_normalize_degenerate() {
        let n = vec3_normalize([0.0, 0.0, 0.0]);
        assert!((vec3_len(n)).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_cross() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = vec3_cross(x, y);
        assert!((z[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_vec3_lerp_endpoints() {
        let a = [0.0; 3];
        let b = [1.0; 3];
        let mid = vec3_lerp(a, b, 0.5);
        assert!((mid[0] - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_quat_normalize_unit() {
        let q = quat_normalize([0.0, 0.0, 0.0, 2.0]);
        assert!((q[3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_quat_slerp_at_zero() {
        let a = [0.0, 0.0, 0.0, 1.0];
        let b = [0.0, 1.0, 0.0, 0.0];
        let r = quat_slerp(a, b, 0.0);
        assert!((r[3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_slerp_at_one() {
        let a = [0.0, 0.0, 0.0, 1.0];
        let b = [0.0, 1.0, 0.0, 0.0];
        let r = quat_slerp(a, b, 1.0);
        assert!((r[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_mul_identity() {
        let id = [0.0, 0.0, 0.0, 1.0];
        let q = [0.5, 0.5, 0.5, 0.5];
        let r = quat_mul(id, q);
        for (a, b) in r.iter().zip(q.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
    #[test]
    fn test_easing_endpoints() {
        let easings = [
            Easing::Linear,
            Easing::EaseInQuad,
            Easing::EaseOutQuad,
            Easing::EaseInOutQuad,
            Easing::EaseInCubic,
            Easing::EaseOutCubic,
            Easing::EaseInOutCubic,
            Easing::EaseInSine,
            Easing::EaseOutSine,
            Easing::EaseInOutSine,
        ];
        for e in &easings {
            assert!(apply_easing(*e, 0.0).abs() < 1e-10, "{e:?} at 0");
            assert!((apply_easing(*e, 1.0) - 1.0).abs() < 1e-10, "{e:?} at 1");
        }
    }
    #[test]
    fn test_ease_out_bounce_endpoints() {
        assert!(apply_easing(Easing::EaseOutBounce, 0.0).abs() < 1e-10);
        assert!((apply_easing(Easing::EaseOutBounce, 1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_ease_out_elastic_endpoints() {
        assert!(apply_easing(Easing::EaseOutElastic, 0.0).abs() < 1e-10);
        assert!((apply_easing(Easing::EaseOutElastic, 1.0) - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_cubic_bezier_linear() {
        for i in 0..=10 {
            let t = i as f64 / 10.0;
            let v = cubic_bezier(0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0, t);
            assert!((v - t).abs() < 1e-10, "bezier {v} vs {t}");
        }
    }
    #[test]
    fn test_cubic_hermite_endpoints() {
        assert!((cubic_hermite(0.0, 1.0, 5.0, 1.0, 0.0)).abs() < 1e-12);
        assert!((cubic_hermite(0.0, 1.0, 5.0, 1.0, 1.0) - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_curve_linear_interp() {
        let mut c = ScalarCurve::new();
        c.push_key(ScalarKey::linear(0.0, 0.0));
        c.push_key(ScalarKey::linear(1.0, 10.0));
        assert!((c.sample(0.5) - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_scalar_curve_clamp() {
        let mut c = ScalarCurve::new();
        c.push_key(ScalarKey::linear(0.0, 3.0));
        c.push_key(ScalarKey::linear(1.0, 7.0));
        assert!((c.sample(-1.0) - 3.0).abs() < 1e-12);
        assert!((c.sample(2.0) - 7.0).abs() < 1e-12);
    }
    #[test]
    fn test_scalar_curve_hermite() {
        let mut c = ScalarCurve::new();
        c.push_key(ScalarKey::hermite(0.0, 0.0, 0.0, 0.0));
        c.push_key(ScalarKey::hermite(1.0, 1.0, 0.0, 0.0));
        let v = c.sample(0.5);
        assert!(v > 0.3 && v < 0.7, "hermite midpoint {v}");
    }
    #[test]
    fn test_vec3_curve_lerp() {
        let mut c = Vec3Curve::new();
        c.push_key(Vec3Key::linear(0.0, [0.0, 0.0, 0.0]));
        c.push_key(Vec3Key::linear(2.0, [2.0, 4.0, 6.0]));
        let v = c.sample(1.0);
        assert!((v[0] - 1.0).abs() < 1e-10);
        assert!((v[1] - 2.0).abs() < 1e-10);
        assert!((v[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_quat_curve_slerp_midpoint() {
        let mut c = QuatCurve::new();
        c.push_key(QuatKey::linear(0.0, [0.0, 0.0, 0.0, 1.0]));
        c.push_key(QuatKey::linear(1.0, [0.0, 0.0, 0.0, 1.0]));
        let q = c.sample(0.5);
        assert!((q[3] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_clip_normalise_time_loop() {
        let clip = AnimationClip::new("walk", 2.0, 30.0, true);
        assert!((clip.normalise_time(2.5) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_clip_normalise_time_clamp() {
        let clip = AnimationClip::new("idle", 1.0, 24.0, false);
        assert!((clip.normalise_time(5.0) - 1.0).abs() < 1e-12);
        assert!((clip.normalise_time(-1.0) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_clip_sample_empty() {
        let clip = AnimationClip::new("empty", 1.0, 30.0, false);
        let pose = clip.sample(0.5);
        assert!(pose.is_empty());
    }
    #[test]
    fn test_skeleton_find_bone() {
        let mut skel = Skeleton::new();
        skel.add_bone(Bone::root(0, "hip"));
        assert!(skel.find_bone("hip").is_some());
        assert!(skel.find_bone("head").is_none());
    }
    #[test]
    fn test_skeleton_world_matrices_identity() {
        let mut skel = Skeleton::new();
        skel.add_bone(Bone::root(0, "root"));
        let local = [identity_mat4()];
        let world = skel.compute_world_matrices(&local);
        for (a, b) in world[0].iter().zip(identity_mat4().iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
    #[test]
    fn test_skin_vertex_single_bone() {
        let pos = [1.0, 0.0, 0.0];
        let weights = [SkinWeight::new(0, 1.0)];
        let m = identity_mat4();
        let skinned = skin_vertex(pos, &weights, &[m]);
        assert!((skinned[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_skin_vertex_two_bones_equal() {
        let pos = [2.0, 0.0, 0.0];
        let weights = [SkinWeight::new(0, 0.5), SkinWeight::new(1, 0.5)];
        let m = identity_mat4();
        let skinned = skin_vertex(pos, &weights, &[m, m]);
        assert!((skinned[0] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_ccd_ik_already_at_target() {
        let joints = vec![IkJoint::new([0.0, 0.0, 0.0]), IkJoint::new([1.0, 0.0, 0.0])];
        let target = [1.0, 0.0, 0.0];
        let result = ccd_ik(&joints, &[1.0], target, 10, 1e-6);
        assert_eq!(result.len(), 2);
        assert!((result[1][0] - 1.0).abs() < 0.01);
    }
    #[test]
    fn test_ccd_ik_empty_joints() {
        let result = ccd_ik(&[], &[], [1.0, 0.0, 0.0], 5, 1e-4);
        assert!(result.is_empty());
    }
    #[test]
    fn test_mocap_empty_returns_none() {
        let seq = MocapSequence::new(30.0, false);
        assert!(seq.sample(0.5).is_none());
    }
    #[test]
    fn test_mocap_single_frame_clamp() {
        let mut seq = MocapSequence::new(30.0, false);
        let mut f = MocapFrame::new(0.0, [0.0; 3]);
        f.set_bone(0, [0.0, 0.0, 0.0, 1.0]);
        seq.push_frame(f);
        let result = seq.sample(5.0).unwrap();
        assert!((result.time).abs() < 1e-12);
    }
    #[test]
    fn test_mocap_interpolation_root() {
        let mut seq = MocapSequence::new(30.0, false);
        seq.push_frame(MocapFrame::new(0.0, [0.0, 0.0, 0.0]));
        seq.push_frame(MocapFrame::new(1.0, [1.0, 0.0, 0.0]));
        let result = seq.sample(0.5).unwrap();
        assert!((result.root_position[0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_lod_for_distance() {
        assert_eq!(lod_for_distance(5.0), AnimLod::Full);
        assert_eq!(lod_for_distance(20.0), AnimLod::Medium);
        assert_eq!(lod_for_distance(50.0), AnimLod::Low);
        assert_eq!(lod_for_distance(100.0), AnimLod::Disabled);
    }
    #[test]
    fn test_lod_update_rate() {
        assert!((AnimLodManager::update_rate(AnimLod::Full) - 1.0).abs() < 1e-12);
        assert!((AnimLodManager::update_rate(AnimLod::Disabled)).abs() < 1e-12);
    }
    #[test]
    fn test_lod_manager_update_entity() {
        let mut mgr = AnimLodManager::new();
        mgr.update_entity(42, 5.0);
        assert_eq!(mgr.get_lod(42), AnimLod::Full);
        mgr.update_entity(42, 100.0);
        assert_eq!(mgr.get_lod(42), AnimLod::Disabled);
    }
    #[test]
    fn test_state_machine_entry() {
        let mut sm = AnimStateMachine::new();
        sm.add_state(AnimState::new("idle", "idle_clip", 1.0));
        sm.set_entry("idle");
        assert_eq!(sm.current.as_deref(), Some("idle"));
    }
    #[test]
    fn test_state_machine_trigger_transition() {
        let mut sm = AnimStateMachine::new();
        sm.add_state(AnimState::new("idle", "idle_clip", 1.0));
        sm.add_state(AnimState::new("run", "run_clip", 1.0));
        sm.add_transition(Transition::triggered("idle", "run", 0.2, "start_run"));
        sm.set_entry("idle");
        let mut idle_clip = AnimationClip::new("idle_clip", 2.0, 30.0, true);
        idle_clip.add_node(NodeAnimation::new(0));
        let mut clips = HashMap::new();
        clips.insert("idle_clip".to_string(), idle_clip);
        sm.fire_trigger("start_run");
        sm.update(0.016, &clips);
        assert_eq!(sm.current.as_deref(), Some("run"));
    }
    #[test]
    fn test_blend_node_missing_clip() {
        let node = BlendNode::Clip {
            clip_name: "nonexistent".to_string(),
            time: 0.0,
        };
        let clips = HashMap::new();
        let result = node.evaluate(&clips);
        assert!(result.is_empty());
    }
    #[test]
    fn test_blend_node_lerp_empty() {
        let node = BlendNode::Lerp {
            child_a: Box::new(BlendNode::Clip {
                clip_name: "a".into(),
                time: 0.0,
            }),
            child_b: Box::new(BlendNode::Clip {
                clip_name: "b".into(),
                time: 0.0,
            }),
            weight: 0.5,
        };
        let clips = HashMap::new();
        let result = node.evaluate(&clips);
        assert!(result.is_empty());
    }
    #[test]
    fn test_trs_matrix_identity_rotation() {
        let m = trs_matrix([0.0; 3], [0.0, 0.0, 0.0, 1.0], [1.0; 3]);
        let id = identity_mat4();
        for (a, b) in m.iter().zip(id.iter()) {
            assert!((a - b).abs() < 1e-12, "trs vs identity: {a} vs {b}");
        }
    }
    #[test]
    fn test_trs_matrix_translation() {
        let m = trs_matrix([5.0, 3.0, 1.0], [0.0, 0.0, 0.0, 1.0], [1.0; 3]);
        assert!((m[12] - 5.0).abs() < 1e-12);
        assert!((m[13] - 3.0).abs() < 1e-12);
        assert!((m[14] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_mat4_mul_identity() {
        let id = identity_mat4();
        let m = trs_matrix([1.0, 2.0, 3.0], [0.0, 0.0, 0.0, 1.0], [1.0; 3]);
        let result = mat4_mul(id, m);
        for (a, b) in result.iter().zip(m.iter()) {
            assert!((a - b).abs() < 1e-12);
        }
    }
    #[test]
    fn test_animator_update_no_crash() {
        let mut skel = Skeleton::new();
        skel.add_bone(Bone::root(0, "root"));
        let mut anim = Animator::new(skel);
        let clip = AnimationClip::new("idle", 1.0, 30.0, true);
        anim.add_clip(clip);
        let idle_state = AnimState::new("idle", "idle", 1.0);
        anim.state_machine.add_state(idle_state);
        anim.state_machine.set_entry("idle");
        anim.update(0.033);
        assert_eq!(anim.skeleton.bones.len(), 1);
    }
}
