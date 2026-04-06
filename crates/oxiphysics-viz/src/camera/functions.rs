//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::CameraKeyframe;

/// Linear interpolation for camera paths.
///
/// Given a list of keyframes sorted by time, interpolates position and
/// target at a given time `t`.
#[allow(dead_code)]
pub fn interpolate_camera_path(keyframes: &[CameraKeyframe], t: f64) -> ([f64; 3], [f64; 3]) {
    if keyframes.is_empty() {
        return ([0.0; 3], [0.0, 0.0, -1.0]);
    }
    if keyframes.len() == 1 || t <= keyframes[0].time {
        return (keyframes[0].position, keyframes[0].target);
    }
    let last = keyframes.len() - 1;
    if t >= keyframes[last].time {
        return (keyframes[last].position, keyframes[last].target);
    }
    let mut idx = 0;
    for i in 0..last {
        if keyframes[i + 1].time >= t {
            idx = i;
            break;
        }
    }
    let a = &keyframes[idx];
    let b = &keyframes[idx + 1];
    let alpha = if (b.time - a.time).abs() > 1e-12 {
        (t - a.time) / (b.time - a.time)
    } else {
        0.0
    };
    let pos = lerp3(a.position, b.position, alpha);
    let tgt = lerp3(a.target, b.target, alpha);
    (pos, tgt)
}
/// Catmull-Rom interpolation for camera paths (smoother than linear).
#[allow(dead_code)]
pub fn interpolate_camera_path_catmull_rom(
    keyframes: &[CameraKeyframe],
    t: f64,
) -> ([f64; 3], [f64; 3]) {
    if keyframes.len() < 4 {
        return interpolate_camera_path(keyframes, t);
    }
    let last = keyframes.len() - 1;
    if t <= keyframes[0].time {
        return (keyframes[0].position, keyframes[0].target);
    }
    if t >= keyframes[last].time {
        return (keyframes[last].position, keyframes[last].target);
    }
    let mut seg = 1;
    for i in 1..last {
        if keyframes[i + 1].time >= t {
            seg = i;
            break;
        }
    }
    let p0_idx = if seg > 0 { seg - 1 } else { 0 };
    let p3_idx = (seg + 2).min(last);
    let a = &keyframes[seg];
    let b = &keyframes[seg + 1];
    let alpha = if (b.time - a.time).abs() > 1e-12 {
        (t - a.time) / (b.time - a.time)
    } else {
        0.0
    };
    let pos = catmull_rom_3(
        keyframes[p0_idx].position,
        a.position,
        b.position,
        keyframes[p3_idx].position,
        alpha,
    );
    let tgt = catmull_rom_3(
        keyframes[p0_idx].target,
        a.target,
        b.target,
        keyframes[p3_idx].target,
        alpha,
    );
    (pos, tgt)
}
/// Smooth noise function on a 1-D lattice using cosine interpolation.
pub(super) fn smooth_noise_1d(t: f64) -> f64 {
    let i = t.floor() as i64;
    let f = t - t.floor();
    let u = f * f * (3.0 - 2.0 * f);
    let a = lattice_noise(i);
    let b = lattice_noise(i + 1);
    a + (b - a) * u
}
/// Deterministic pseudo-random noise value in \[-1, 1\] for lattice index `i`.
pub(super) fn lattice_noise(i: i64) -> f64 {
    let x = i.wrapping_mul(1664525).wrapping_add(1013904223) as u32;
    let x = x ^ (x >> 16);
    let x = x.wrapping_mul(0x45d9f3b);
    let x = x ^ (x >> 16);
    (x as f64 / u32::MAX as f64) * 2.0 - 1.0
}
/// Compute the Halton sequence value for index `index` in base `base`.
#[allow(dead_code)]
pub fn halton(mut index: usize, base: usize) -> f64 {
    let mut result = 0.0_f64;
    let mut denom = 1.0_f64;
    while index > 0 {
        denom *= base as f64;
        result += (index % base) as f64 / denom;
        index /= base;
    }
    result
}
pub(super) fn add4(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]]
}
pub(super) fn sub4(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]]
}
/// Linear interpolation between two 3-vectors.
#[allow(dead_code)]
pub fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}
/// Add two 3-vectors.
#[allow(dead_code)]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector by a scalar.
#[allow(dead_code)]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Catmull-Rom interpolation for a single 3-vector component.
#[allow(dead_code)]
pub(super) fn catmull_rom_3(
    p0: [f64; 3],
    p1: [f64; 3],
    p2: [f64; 3],
    p3: [f64; 3],
    t: f64,
) -> [f64; 3] {
    let t2 = t * t;
    let t3 = t2 * t;
    let mut result = [0.0; 3];
    for i in 0..3 {
        result[i] = 0.5
            * ((2.0 * p1[i])
                + (-p0[i] + p2[i]) * t
                + (2.0 * p0[i] - 5.0 * p1[i] + 4.0 * p2[i] - p3[i]) * t2
                + (-p0[i] + 3.0 * p1[i] - 3.0 * p2[i] + p3[i]) * t3);
    }
    result
}
/// Length of a 3-vector.
#[allow(dead_code)]
pub fn length3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
/// 4x4 matrix multiply (column-major).
#[allow(dead_code)]
pub fn mat4_mul(a: [[f64; 4]; 4], b: [[f64; 4]; 4]) -> [[f64; 4]; 4] {
    let mut out = [[0.0f64; 4]; 4];
    for col in 0..4 {
        for row in 0..4 {
            out[col][row] = (0..4).map(|k| a[k][row] * b[col][k]).sum();
        }
    }
    out
}
/// Normalize a 3-vector.
#[allow(dead_code)]
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len > 1e-12 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        v
    }
}
/// Dot product of two 3-vectors.
#[allow(dead_code)]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Cross product of two 3-vectors.
#[allow(dead_code)]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Subtract two 3-vectors.
#[allow(dead_code)]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::camera::Camera;
    use crate::camera::CameraRig;
    use crate::camera::CameraShake;
    use crate::camera::CinematicMove;
    use crate::camera::CinematicMoveType;
    use crate::camera::DepthOfField;
    use crate::camera::FlyCamera;
    use crate::camera::Frustum;
    use crate::camera::MultiTargetCamera;
    use crate::camera::PerlinShake;
    use crate::camera::PhysicalDoF;
    use crate::camera::TaaHistory;
    use std::f64::consts::PI;
    fn default_camera() -> Camera {
        Camera::new(
            [0.0, 0.0, 5.0],
            [0.0, 0.0, 0.0],
            60.0,
            16.0 / 9.0,
            0.1,
            100.0,
        )
    }
    #[test]
    fn test_normalize3_unit_length() {
        let v = normalize3([3.0, 4.0, 0.0]);
        let len = dot3(v, v).sqrt();
        assert!((len - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalize3_zero_passthrough() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert_eq!(v, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_dot3_orthogonal() {
        assert!((dot3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])).abs() < 1e-12);
    }
    #[test]
    fn test_cross3_standard_basis() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_sub3() {
        let r = sub3([3.0, 2.0, 1.0], [1.0, 1.0, 1.0]);
        assert_eq!(r, [2.0, 1.0, 0.0]);
    }
    #[test]
    fn test_mat4_mul_identity() {
        let id = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let m = default_camera().view_matrix();
        let result = mat4_mul(id, m);
        for col in 0..4 {
            for row in 0..4 {
                assert!((result[col][row] - m[col][row]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn test_camera_fov_converted_to_radians() {
        let cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 90.0, 1.0, 0.1, 100.0);
        assert!((cam.fov_y - PI / 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_camera_default_up_is_y() {
        let cam = default_camera();
        assert_eq!(cam.up, [0.0, 1.0, 0.0]);
    }
    #[test]
    fn test_view_matrix_column_major_homogeneous() {
        let m = default_camera().view_matrix();
        assert!((m[3][3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_view_matrix_orthonormal_rows() {
        let m = default_camera().view_matrix();
        let r: [f64; 3] = [m[0][0], m[1][0], m[2][0]];
        let u: [f64; 3] = [m[0][1], m[1][1], m[2][1]];
        let f: [f64; 3] = [m[0][2], m[1][2], m[2][2]];
        assert!((dot3(r, r).sqrt() - 1.0).abs() < 1e-12);
        assert!((dot3(u, u).sqrt() - 1.0).abs() < 1e-12);
        assert!((dot3(f, f).sqrt() - 1.0).abs() < 1e-12);
        assert!(dot3(r, u).abs() < 1e-12);
    }
    #[test]
    fn test_projection_matrix_near_maps_to_minus_one() {
        let cam = Camera::new([0.0, 0.0, 0.0], [0.0, 0.0, -1.0], 90.0, 1.0, 1.0, 10.0);
        let p = cam.projection_matrix();
        let n = cam.near;
        let clip_z = p[2][2] * (-n) + p[3][2];
        let clip_w = p[2][3] * (-n) + p[3][3];
        let ndc_z = clip_z / clip_w;
        assert!((ndc_z - (-1.0)).abs() < 1e-10, "ndc_z={}", ndc_z);
    }
    #[test]
    fn test_projection_matrix_far_maps_to_plus_one() {
        let cam = Camera::new([0.0, 0.0, 0.0], [0.0, 0.0, -1.0], 90.0, 1.0, 1.0, 10.0);
        let p = cam.projection_matrix();
        let fa = cam.far;
        let clip_z = p[2][2] * (-fa) + p[3][2];
        let clip_w = p[2][3] * (-fa) + p[3][3];
        let ndc_z = clip_z / clip_w;
        assert!((ndc_z - 1.0).abs() < 1e-10, "ndc_z={}", ndc_z);
    }
    #[test]
    fn test_forward_is_toward_target() {
        let cam = default_camera();
        let fwd = cam.forward();
        assert!((fwd[0]).abs() < 1e-12);
        assert!((fwd[1]).abs() < 1e-12);
        assert!((fwd[2] + 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_right_perpendicular_to_forward() {
        let cam = default_camera();
        assert!(dot3(cam.right(), cam.forward()).abs() < 1e-12);
    }
    #[test]
    fn test_up_vector_perpendicular_to_forward() {
        let cam = default_camera();
        assert!(dot3(cam.up_vector(), cam.forward()).abs() < 1e-12);
    }
    #[test]
    fn test_up_vector_perpendicular_to_right() {
        let cam = default_camera();
        assert!(dot3(cam.up_vector(), cam.right()).abs() < 1e-12);
    }
    #[test]
    fn test_orbit_preserves_distance_to_target() {
        let mut cam = default_camera();
        let dist_before = sub3(cam.position, cam.target);
        let r_before = dot3(dist_before, dist_before).sqrt();
        cam.orbit(0.3, 0.1);
        let dist_after = sub3(cam.position, cam.target);
        let r_after = dot3(dist_after, dist_after).sqrt();
        assert!((r_before - r_after).abs() < 1e-10);
    }
    #[test]
    fn test_orbit_changes_position() {
        let mut cam = default_camera();
        let pos_before = cam.position;
        cam.orbit(0.5, 0.2);
        let changed = pos_before
            .iter()
            .zip(cam.position.iter())
            .any(|(a, b)| (a - b).abs() > 1e-10);
        assert!(changed);
    }
    #[test]
    fn test_zoom_moves_position_forward() {
        let mut cam = default_camera();
        let pos_before = cam.position;
        cam.zoom(1.0);
        assert!(cam.position[2] < pos_before[2]);
    }
    #[test]
    fn test_zoom_preserves_forward_direction() {
        let mut cam = default_camera();
        let fwd_before = cam.forward();
        cam.zoom(2.0);
        let fwd_after = cam.forward();
        for i in 0..3 {
            assert!((fwd_before[i] - fwd_after[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_frustum_contains_target() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        assert!(frustum.contains_point(cam.target));
    }
    #[test]
    fn test_frustum_excludes_behind_camera() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        assert!(!frustum.contains_point([0.0, 0.0, 1000.0]));
    }
    #[test]
    fn test_frustum_sphere_large_enough_always_visible() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        assert!(frustum.contains_sphere([0.0, 0.0, 0.0], 1000.0));
    }
    #[test]
    fn test_frustum_sphere_behind_camera_excluded() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        assert!(!frustum.contains_sphere([0.0, 0.0, 1000.0], 0.1));
    }
    #[test]
    fn test_frustum_aabb_at_origin_visible() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        assert!(frustum.contains_aabb([-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]));
    }
    #[test]
    fn test_frustum_aabb_behind_camera_excluded() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        assert!(!frustum.contains_aabb([0.0, 0.0, 200.0], [1.0, 1.0, 201.0]));
    }
    #[test]
    fn test_distance_to_target() {
        let cam = default_camera();
        assert!((cam.distance_to_target() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_pan_moves_both_position_and_target() {
        let mut cam = default_camera();
        let pos0 = cam.position;
        let tgt0 = cam.target;
        cam.pan(1.0, 0.0);
        let dp = sub3(cam.position, pos0);
        let dt = sub3(cam.target, tgt0);
        for i in 0..3 {
            assert!(
                (dp[i] - dt[i]).abs() < 1e-12,
                "pan should move pos and target equally"
            );
        }
    }
    #[test]
    fn test_pan_preserves_distance() {
        let mut cam = default_camera();
        let d_before = cam.distance_to_target();
        cam.pan(2.0, 1.5);
        let d_after = cam.distance_to_target();
        assert!((d_before - d_after).abs() < 1e-10);
    }
    #[test]
    fn test_look_at_distance() {
        let mut cam = default_camera();
        cam.look_at_distance([1.0, 2.0, 3.0], 10.0);
        assert_eq!(cam.target, [1.0, 2.0, 3.0]);
        assert!((cam.position[2] - 13.0).abs() < 1e-12);
        assert!((cam.distance_to_target() - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_fly_camera_initial_forward() {
        let cam = FlyCamera::new([0.0, 0.0, 0.0], 60.0, 16.0 / 9.0);
        let fwd = cam.forward();
        assert!((fwd[0]).abs() < 1e-12);
        assert!((fwd[1]).abs() < 1e-12);
        assert!((fwd[2] + 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fly_camera_move_forward() {
        let mut cam = FlyCamera::new([0.0, 0.0, 0.0], 60.0, 1.0);
        cam.speed = 1.0;
        cam.move_forward(1.0);
        assert!(cam.position[2] < -0.5);
    }
    #[test]
    fn test_fly_camera_strafe() {
        let mut cam = FlyCamera::new([0.0, 0.0, 0.0], 60.0, 1.0);
        cam.speed = 1.0;
        let pos_before = cam.position;
        cam.move_right(1.0);
        let changed = pos_before
            .iter()
            .zip(cam.position.iter())
            .any(|(a, b)| (a - b).abs() > 1e-6);
        assert!(changed, "strafe should change position");
    }
    #[test]
    fn test_fly_camera_move_up() {
        let mut cam = FlyCamera::new([0.0, 0.0, 0.0], 60.0, 1.0);
        cam.speed = 1.0;
        cam.move_up(1.0);
        assert!((cam.position[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fly_camera_rotate_pitch_clamp() {
        let mut cam = FlyCamera::new([0.0, 0.0, 0.0], 60.0, 1.0);
        cam.sensitivity = 1.0;
        cam.rotate(0.0, -100.0);
        assert!(cam.pitch < PI / 2.0, "pitch should be clamped below PI/2");
        assert!(cam.pitch > -PI / 2.0, "pitch should be clamped above -PI/2");
    }
    #[test]
    fn test_fly_camera_view_matrix_valid() {
        let cam = FlyCamera::new([0.0, 0.0, 5.0], 60.0, 1.0);
        let m = cam.view_matrix();
        assert!((m[3][3] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_fly_camera_projection_matrix_valid() {
        let cam = FlyCamera::new([0.0, 0.0, 5.0], 60.0, 1.0);
        let p = cam.projection_matrix();
        assert!(p[0][0].abs() > 0.1);
    }
    #[test]
    fn test_interpolate_empty_keyframes() {
        let (pos, tgt) = interpolate_camera_path(&[], 0.5);
        assert_eq!(pos, [0.0; 3]);
        assert_eq!(tgt, [0.0, 0.0, -1.0]);
    }
    #[test]
    fn test_interpolate_single_keyframe() {
        let kf = vec![CameraKeyframe {
            time: 0.0,
            position: [1.0, 2.0, 3.0],
            target: [0.0, 0.0, 0.0],
        }];
        let (pos, tgt) = interpolate_camera_path(&kf, 5.0);
        assert_eq!(pos, [1.0, 2.0, 3.0]);
        assert_eq!(tgt, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn test_interpolate_midpoint() {
        let kf = vec![
            CameraKeyframe {
                time: 0.0,
                position: [0.0, 0.0, 0.0],
                target: [0.0, 0.0, -1.0],
            },
            CameraKeyframe {
                time: 1.0,
                position: [10.0, 0.0, 0.0],
                target: [10.0, 0.0, -1.0],
            },
        ];
        let (pos, tgt) = interpolate_camera_path(&kf, 0.5);
        assert!((pos[0] - 5.0).abs() < 1e-10);
        assert!((tgt[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_interpolate_before_first() {
        let kf = vec![
            CameraKeyframe {
                time: 1.0,
                position: [1.0, 0.0, 0.0],
                target: [0.0; 3],
            },
            CameraKeyframe {
                time: 2.0,
                position: [2.0, 0.0, 0.0],
                target: [0.0; 3],
            },
        ];
        let (pos, _) = interpolate_camera_path(&kf, 0.0);
        assert_eq!(pos, [1.0, 0.0, 0.0]);
    }
    #[test]
    fn test_interpolate_after_last() {
        let kf = vec![
            CameraKeyframe {
                time: 0.0,
                position: [0.0; 3],
                target: [0.0; 3],
            },
            CameraKeyframe {
                time: 1.0,
                position: [10.0, 0.0, 0.0],
                target: [0.0; 3],
            },
        ];
        let (pos, _) = interpolate_camera_path(&kf, 5.0);
        assert_eq!(pos, [10.0, 0.0, 0.0]);
    }
    #[test]
    fn test_catmull_rom_falls_back_for_few_keyframes() {
        let kf = vec![
            CameraKeyframe {
                time: 0.0,
                position: [0.0; 3],
                target: [0.0; 3],
            },
            CameraKeyframe {
                time: 1.0,
                position: [10.0, 0.0, 0.0],
                target: [0.0; 3],
            },
        ];
        let (pos, _) = interpolate_camera_path_catmull_rom(&kf, 0.5);
        assert!((pos[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_catmull_rom_with_four_keyframes() {
        let kf = vec![
            CameraKeyframe {
                time: 0.0,
                position: [0.0, 0.0, 0.0],
                target: [0.0; 3],
            },
            CameraKeyframe {
                time: 1.0,
                position: [1.0, 0.0, 0.0],
                target: [0.0; 3],
            },
            CameraKeyframe {
                time: 2.0,
                position: [2.0, 0.0, 0.0],
                target: [0.0; 3],
            },
            CameraKeyframe {
                time: 3.0,
                position: [3.0, 0.0, 0.0],
                target: [0.0; 3],
            },
        ];
        let (pos, _) = interpolate_camera_path_catmull_rom(&kf, 1.5);
        assert!(
            (pos[0] - 1.5).abs() < 0.2,
            "CR interpolation at 1.5 got {}",
            pos[0]
        );
    }
    #[test]
    fn test_camera_shake_initial_offset() {
        let shake = CameraShake::new(1.0, 2.0, 10.0);
        let off = shake.offset();
        assert!((off[0]).abs() < 1e-12);
    }
    #[test]
    fn test_camera_shake_decays() {
        let mut shake = CameraShake::new(1.0, 5.0, 10.0);
        shake.update(0.1);
        assert!(shake.intensity() < 1.0, "intensity should decrease");
    }
    #[test]
    fn test_camera_shake_impulse() {
        let mut shake = CameraShake::new(0.0, 2.0, 10.0);
        assert!(shake.is_done());
        shake.impulse(1.0);
        assert!(!shake.is_done());
        assert!((shake.intensity() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_camera_shake_eventually_done() {
        let mut shake = CameraShake::new(1.0, 10.0, 5.0);
        for _ in 0..100 {
            shake.update(0.05);
        }
        assert!(shake.is_done(), "shake should fade out");
    }
    #[test]
    fn test_dof_in_focus() {
        let dof = DepthOfField::new(10.0, 4.0, 5.0);
        assert!(
            (dof.blur_factor(10.0)).abs() < 1e-12,
            "focal plane should be in focus"
        );
        assert!((dof.blur_factor(9.0)).abs() < 1e-12, "within focal range");
        assert!((dof.blur_factor(11.0)).abs() < 1e-12, "within focal range");
    }
    #[test]
    fn test_dof_blur_increases_with_distance() {
        let dof = DepthOfField::new(10.0, 2.0, 5.0);
        let b_near = dof.blur_factor(5.0);
        let b_far = dof.blur_factor(20.0);
        assert!(b_near > 0.0, "outside focal range should blur");
        assert!(b_far > 0.0, "far from focal should blur");
    }
    #[test]
    fn test_dof_coc_radius() {
        let dof = DepthOfField::new(10.0, 2.0, 8.0);
        let coc = dof.coc_radius(10.0);
        assert!((coc).abs() < 1e-12, "at focus, CoC should be 0");
        let coc_far = dof.coc_radius(100.0);
        assert!(coc_far > 0.0);
    }
    #[test]
    fn test_lerp3() {
        let r = lerp3([0.0, 0.0, 0.0], [10.0, 20.0, 30.0], 0.5);
        assert!((r[0] - 5.0).abs() < 1e-12);
        assert!((r[1] - 10.0).abs() < 1e-12);
        assert!((r[2] - 15.0).abs() < 1e-12);
    }
    #[test]
    fn test_add3() {
        let r = add3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }
    #[test]
    fn test_scale3() {
        let r = scale3([1.0, 2.0, 3.0], 2.0);
        assert_eq!(r, [2.0, 4.0, 6.0]);
    }
    #[test]
    fn test_length3() {
        let l = length3([3.0, 4.0, 0.0]);
        assert!((l - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_frustum_signed_distance() {
        let cam = default_camera();
        let frustum = Frustum::from_camera(&cam);
        let d = frustum.signed_distance(4, cam.target);
        assert!(d > 0.0, "target should be in front of near plane, got {d}");
    }
    #[test]
    fn test_perlin_shake_initial_amplitude() {
        let shake = PerlinShake::new(2.0, 1.0, 5.0);
        assert!((shake.amplitude() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_perlin_shake_decays_over_time() {
        let mut shake = PerlinShake::new(1.0, 2.0, 5.0);
        shake.update(0.5);
        assert!(shake.amplitude() < 1.0, "amplitude should decrease");
    }
    #[test]
    fn test_perlin_shake_impulse_increases_amplitude() {
        let mut shake = PerlinShake::new(0.0, 1.0, 5.0);
        assert!(shake.is_done());
        shake.impulse(0.5);
        assert!(!shake.is_done());
    }
    #[test]
    fn test_perlin_shake_offset_bounds() {
        let shake = PerlinShake::new(1.0, 0.1, 10.0);
        let off = shake.offset();
        for v in &off {
            assert!(v.abs() <= 1.01, "offset out of range: {v}");
        }
    }
    #[test]
    fn test_perlin_shake_eventually_done() {
        let mut shake = PerlinShake::new(1.0, 10.0, 5.0);
        for _ in 0..200 {
            shake.update(0.05);
        }
        assert!(shake.is_done(), "perlin shake should fade out");
    }
    #[test]
    fn test_taa_history_starts_empty() {
        let h = TaaHistory::new(8);
        assert!(h.is_empty());
        assert_eq!(h.len(), 0);
    }
    #[test]
    fn test_taa_history_push_and_len() {
        let mut h = TaaHistory::new(4);
        h.push([1.0, 0.0, 0.0], [0.0; 3]);
        h.push([2.0, 0.0, 0.0], [0.0; 3]);
        assert_eq!(h.len(), 2);
    }
    #[test]
    fn test_taa_history_cap_overflow() {
        let mut h = TaaHistory::new(3);
        for i in 0..10 {
            h.push([i as f64, 0.0, 0.0], [0.0; 3]);
        }
        assert_eq!(h.len(), 3, "history should not exceed capacity");
    }
    #[test]
    fn test_taa_history_blended_position_single() {
        let mut h = TaaHistory::new(4);
        h.push([5.0, 0.0, 0.0], [0.0; 3]);
        let pos = h.blended_position(0.5).unwrap();
        assert!((pos[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_taa_history_blended_position_empty_is_none() {
        let h = TaaHistory::new(4);
        assert!(h.blended_position(0.1).is_none());
    }
    #[test]
    fn test_taa_halton_jitter_bounded() {
        for i in 0..32 {
            let (jx, jy) = TaaHistory::halton_jitter(i, 1920, 1080);
            assert!(jx.abs() < 1.0, "jx out of range: {jx}");
            assert!(jy.abs() < 1.0, "jy out of range: {jy}");
        }
    }
    #[test]
    fn test_halton_sequence_base2() {
        let h1 = halton(1, 2);
        assert!((h1 - 0.5).abs() < 1e-12);
        let h2 = halton(2, 2);
        assert!((h2 - 0.25).abs() < 1e-12);
    }
    #[test]
    fn test_smooth_noise_bounded() {
        for i in 0..100 {
            let v = smooth_noise_1d(i as f64 * 0.1);
            assert!(v.abs() <= 1.01, "smooth_noise_1d out of range: {v}");
        }
    }
    #[test]
    fn test_cinematic_move_progress_starts_zero() {
        let m = CinematicMove::new(CinematicMoveType::Dolly, 5.0, 2.0);
        assert!((m.progress()).abs() < 1e-12);
        assert!(!m.is_done());
    }
    #[test]
    fn test_cinematic_dolly_moves_camera() {
        let mut cam = Camera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut m = CinematicMove::new(CinematicMoveType::Dolly, 2.0, 1.0);
        m.advance(&mut cam, 1.0);
        assert!(cam.position[2] < 10.0, "dolly should move camera forward");
    }
    #[test]
    fn test_cinematic_crane_moves_up() {
        let mut cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut m = CinematicMove::new(CinematicMoveType::Crane, 3.0, 1.0);
        m.advance(&mut cam, 1.0);
        assert!(cam.position[1] > 0.0, "crane should lift camera");
    }
    #[test]
    fn test_cinematic_jib_pans() {
        let mut cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let pos0 = cam.position;
        let mut m = CinematicMove::new(CinematicMoveType::Jib, 2.0, 1.0);
        m.advance(&mut cam, 1.0);
        let changed = pos0
            .iter()
            .zip(cam.position.iter())
            .any(|(a, b)| (a - b).abs() > 1e-6);
        assert!(changed, "jib should move camera laterally");
    }
    #[test]
    fn test_cinematic_arc_orbits() {
        let mut cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let pos0 = cam.position;
        let mut m = CinematicMove::new(CinematicMoveType::Arc, 0.3, 1.0);
        m.advance(&mut cam, 1.0);
        let changed = pos0
            .iter()
            .zip(cam.position.iter())
            .any(|(a, b)| (a - b).abs() > 1e-6);
        assert!(changed, "arc should orbit camera");
    }
    #[test]
    fn test_cinematic_move_is_done_after_full_duration() {
        let mut cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut m = CinematicMove::new(CinematicMoveType::Dolly, 1.0, 0.5);
        m.advance(&mut cam, 0.5);
        assert!(m.is_done());
    }
    #[test]
    fn test_camera_rig_starts_idle() {
        let cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let rig = CameraRig::new(cam);
        assert!(rig.is_idle());
        assert_eq!(rig.pending_count(), 0);
    }
    #[test]
    fn test_camera_rig_enqueue_and_count() {
        let cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut rig = CameraRig::new(cam);
        rig.enqueue(CinematicMove::new(CinematicMoveType::Dolly, 1.0, 1.0));
        rig.enqueue(CinematicMove::new(CinematicMoveType::Crane, 2.0, 1.0));
        assert_eq!(rig.pending_count(), 2);
    }
    #[test]
    fn test_camera_rig_update_advances_move() {
        let cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut rig = CameraRig::new(cam);
        rig.enqueue(CinematicMove::new(CinematicMoveType::Crane, 3.0, 0.5));
        rig.update(0.5);
        assert!(rig.is_idle());
    }
    #[test]
    fn test_camera_rig_sequential_moves() {
        let cam = Camera::new([0.0, 0.0, 5.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut rig = CameraRig::new(cam);
        rig.enqueue(CinematicMove::new(CinematicMoveType::Dolly, 1.0, 0.1));
        rig.enqueue(CinematicMove::new(CinematicMoveType::Crane, 1.0, 0.1));
        rig.update(0.1);
        assert_eq!(rig.pending_count(), 1, "first move should have completed");
        rig.update(0.1);
        assert!(rig.is_idle(), "both moves should be done");
    }
    #[test]
    fn test_multi_target_camera_no_targets() {
        let cam = Camera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mt = MultiTargetCamera::new(cam);
        assert!(mt.bounding_sphere().is_none());
    }
    #[test]
    fn test_multi_target_camera_single_target() {
        let cam = Camera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut mt = MultiTargetCamera::new(cam);
        mt.add_target([3.0, 0.0, 0.0]);
        let (centroid, radius) = mt.bounding_sphere().unwrap();
        assert!((centroid[0] - 3.0).abs() < 1e-10);
        assert!(radius >= mt.padding);
    }
    #[test]
    fn test_multi_target_camera_centroid() {
        let cam = Camera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut mt = MultiTargetCamera::new(cam);
        mt.add_target([-1.0, 0.0, 0.0]);
        mt.add_target([1.0, 0.0, 0.0]);
        let (centroid, _) = mt.bounding_sphere().unwrap();
        assert!((centroid[0]).abs() < 1e-10, "centroid should be at origin");
    }
    #[test]
    fn test_multi_target_camera_update_changes_position() {
        let cam = Camera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut mt = MultiTargetCamera::new(cam);
        mt.add_target([50.0, 0.0, 0.0]);
        mt.smoothing = 0.0;
        mt.update(1.0 / 60.0);
        assert!(mt.camera.position[0].abs() > 0.1 || mt.camera.position[2].abs() > 5.0);
    }
    #[test]
    fn test_multi_target_camera_clear() {
        let cam = Camera::new([0.0, 0.0, 10.0], [0.0, 0.0, 0.0], 60.0, 1.0, 0.1, 100.0);
        let mut mt = MultiTargetCamera::new(cam);
        mt.add_target([1.0, 2.0, 3.0]);
        mt.clear_targets();
        assert_eq!(mt.target_count(), 0);
        assert!(mt.bounding_sphere().is_none());
    }
    #[test]
    fn test_physical_dof_in_focus_is_zero() {
        let dof = PhysicalDoF::portrait();
        let coc = dof.coc_mm(dof.focus_distance);
        assert!(coc.abs() < 0.01, "CoC at focus should be ~0, got {coc}");
    }
    #[test]
    fn test_physical_dof_blur_increases_with_distance() {
        let dof = PhysicalDoF::portrait();
        let b_near = dof.blur_factor(1.0);
        let b_far = dof.blur_factor(10.0);
        assert!(
            b_near > 0.0 || b_far > 0.0,
            "at least one should be blurred"
        );
    }
    #[test]
    fn test_physical_dof_range_valid() {
        let dof = PhysicalDoF::portrait();
        let (near, far) = dof.dof_range();
        assert!(near >= 0.0, "near DoF must be non-negative");
        assert!(far >= near, "far DoF must be >= near");
    }
    #[test]
    fn test_physical_dof_zero_depth_is_zero() {
        let dof = PhysicalDoF::portrait();
        assert!((dof.coc_mm(0.0)).abs() < 1e-10);
    }
}
