//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// A 3-D point or vector as `[x, y, z]` in metres.
pub type Vec3 = [f64; 3];
/// A unit quaternion `[w, x, y, z]`.
pub type Quat = [f64; 4];
/// A homogeneous 4×4 transform stored row-major.
pub type Mat4 = [f64; 16];
/// Compute the dot product of two 3-vectors.
#[inline]
pub fn dot3(a: Vec3, b: Vec3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Compute the Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(v: Vec3) -> f64 {
    dot3(v, v).sqrt()
}
/// Subtract two 3-vectors: `a - b`.
#[inline]
pub fn sub3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Add two 3-vectors: `a + b`.
#[inline]
pub fn add3(a: Vec3, b: Vec3) -> Vec3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Scale a 3-vector by scalar `s`.
#[inline]
pub fn scale3(v: Vec3, s: f64) -> Vec3 {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: Vec3, b: Vec3) -> Vec3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Normalise a 3-vector. Returns `[0,0,0]` for zero vectors.
#[inline]
pub fn normalize3(v: Vec3) -> Vec3 {
    let n = norm3(v);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / n)
    }
}
/// Build a 4×4 identity matrix.
#[inline]
pub fn mat4_identity() -> Mat4 {
    let mut m = [0.0f64; 16];
    m[0] = 1.0;
    m[5] = 1.0;
    m[10] = 1.0;
    m[15] = 1.0;
    m
}
/// Build a 4×4 translation matrix.
#[inline]
pub fn mat4_translation(t: Vec3) -> Mat4 {
    let mut m = mat4_identity();
    m[3] = t[0];
    m[7] = t[1];
    m[11] = t[2];
    m
}
/// Multiply two 4×4 matrices (row-major).
pub fn mat4_mul(a: Mat4, b: Mat4) -> Mat4 {
    let mut c = [0.0f64; 16];
    for row in 0..4 {
        for col in 0..4 {
            let mut s = 0.0;
            for k in 0..4 {
                s += a[row * 4 + k] * b[k * 4 + col];
            }
            c[row * 4 + col] = s;
        }
    }
    c
}
/// Convert a quaternion `[w,x,y,z]` to a rotation matrix embedded in a 4×4.
pub fn quat_to_mat4(q: Quat) -> Mat4 {
    let [w, x, y, z] = q;
    let mut m = mat4_identity();
    m[0] = 1.0 - 2.0 * (y * y + z * z);
    m[1] = 2.0 * (x * y - z * w);
    m[2] = 2.0 * (x * z + y * w);
    m[4] = 2.0 * (x * y + z * w);
    m[5] = 1.0 - 2.0 * (x * x + z * z);
    m[6] = 2.0 * (y * z - x * w);
    m[8] = 2.0 * (x * z - y * w);
    m[9] = 2.0 * (y * z + x * w);
    m[10] = 1.0 - 2.0 * (x * x + y * y);
    m
}
/// Convert roll-pitch-yaw (intrinsic ZYX) to a 4×4 rotation matrix.
pub fn rpy_to_mat4(rpy: Vec3) -> Mat4 {
    let [roll, pitch, yaw] = rpy;
    let cr = roll.cos();
    let sr = roll.sin();
    let cp = pitch.cos();
    let sp = pitch.sin();
    let cy = yaw.cos();
    let sy = yaw.sin();
    let mut m = mat4_identity();
    m[0] = cy * cp;
    m[1] = cy * sp * sr - sy * cr;
    m[2] = cy * sp * cr + sy * sr;
    m[4] = sy * cp;
    m[5] = sy * sp * sr + cy * cr;
    m[6] = sy * sp * cr - cy * sr;
    m[8] = -sp;
    m[9] = cp * sr;
    m[10] = cp * cr;
    m
}
pub(super) const ROS_MAGIC: [u8; 4] = [0x52, 0x4F, 0x53, 0x31];
/// Apply a 4×4 homogeneous transform (row-major) to a 3-D point.
pub fn apply_transform(m: &Mat4, p: Vec3) -> Vec3 {
    let x = m[0] * p[0] + m[1] * p[1] + m[2] * p[2] + m[3];
    let y = m[4] * p[0] + m[5] * p[1] + m[6] * p[2] + m[7];
    let z = m[8] * p[0] + m[9] * p[1] + m[10] * p[2] + m[11];
    let w = m[12] * p[0] + m[13] * p[1] + m[14] * p[2] + m[15];
    if (w - 1.0).abs() > 1e-9 && w.abs() > 1e-9 {
        [x / w, y / w, z / w]
    } else {
        [x, y, z]
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::robotics_io::CameraIntrinsics;
    use crate::robotics_io::GripperState;
    use crate::robotics_io::IkSolution;
    use crate::robotics_io::IkSolutionLog;
    use crate::robotics_io::JointLimits;
    use crate::robotics_io::JointStateMsg;
    use crate::robotics_io::JointType;
    use crate::robotics_io::JointWaypoint;
    use crate::robotics_io::MsgType;
    use crate::robotics_io::RobotPointCloud;
    use crate::robotics_io::RobotTrajectory;
    use crate::robotics_io::RosMsg;
    use crate::robotics_io::ScanPoint;
    use crate::robotics_io::UrdfGeometry;
    use crate::robotics_io::UrdfJoint;
    use crate::robotics_io::UrdfLink;
    use crate::robotics_io::UrdfRobot;
    use crate::robotics_io::UrdfVisualElement;
    use crate::robotics_io::WorkspaceGrid;
    use crate::robotics_io::Wrench;
    #[test]
    fn test_dot3_orthogonal() {
        assert!((dot3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0])).abs() < 1e-12);
    }
    #[test]
    fn test_norm3_unit_vector() {
        let v = normalize3([3.0, 4.0, 0.0]);
        assert!((norm3(v) - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_cross3_right_hand() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[0]).abs() < 1e-12);
        assert!((c[1]).abs() < 1e-12);
        assert!((c[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_add3_sub3() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let s = add3(a, b);
        assert_eq!(s, [5.0, 7.0, 9.0]);
        let d = sub3(b, a);
        assert_eq!(d, [3.0, 3.0, 3.0]);
    }
    #[test]
    fn test_mat4_identity_mul() {
        let id = mat4_identity();
        let t = mat4_translation([1.0, 2.0, 3.0]);
        let result = mat4_mul(id, t);
        for i in 0..16 {
            assert!((result[i] - t[i]).abs() < 1e-12);
        }
    }
    #[test]
    fn test_apply_transform_identity() {
        let id = mat4_identity();
        let p = [1.0, 2.0, 3.0];
        let q = apply_transform(&id, p);
        assert!((q[0] - 1.0).abs() < 1e-12);
        assert!((q[1] - 2.0).abs() < 1e-12);
        assert!((q[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_apply_transform_translation() {
        let t = mat4_translation([5.0, -3.0, 1.0]);
        let p = [0.0, 0.0, 0.0];
        let q = apply_transform(&t, p);
        assert!((q[0] - 5.0).abs() < 1e-12);
        assert!((q[1] + 3.0).abs() < 1e-12);
        assert!((q[2] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_urdf_geometry_bounding_radius_sphere() {
        let elem = UrdfVisualElement::new(UrdfGeometry::Sphere { radius: 0.5 });
        assert!((elem.bounding_radius() - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_urdf_geometry_bounding_radius_box() {
        let elem = UrdfVisualElement::new(UrdfGeometry::Box {
            half_extents: [1.0, 0.0, 0.0],
        });
        assert!((elem.bounding_radius() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_urdf_link_approximate_volume_empty() {
        let link = UrdfLink::new("base");
        assert!((link.approximate_visual_volume()).abs() < 1e-12);
    }
    #[test]
    fn test_urdf_link_with_sphere() {
        let mut link = UrdfLink::new("link1");
        link.visuals
            .push(UrdfVisualElement::new(UrdfGeometry::Sphere { radius: 1.0 }));
        let vol = link.approximate_visual_volume();
        let expected = (4.0 / 3.0) * std::f64::consts::PI;
        assert!((vol - expected).abs() < 1e-6);
    }
    #[test]
    fn test_joint_limits_in_range() {
        let lim = JointLimits {
            lower: -1.0,
            upper: 1.0,
            ..JointLimits::default()
        };
        assert!(lim.in_range(0.0));
        assert!(!lim.in_range(2.0));
    }
    #[test]
    fn test_joint_limits_clamp() {
        let lim = JointLimits {
            lower: 0.0,
            upper: 3.125,
            ..JointLimits::default()
        };
        assert!((lim.clamp(-1.0) - 0.0).abs() < 1e-12);
        assert!((lim.clamp(10.0) - 3.125).abs() < 1e-12);
    }
    #[test]
    fn test_joint_limits_range() {
        let lim = JointLimits {
            lower: -1.5,
            upper: 1.5,
            ..JointLimits::default()
        };
        assert!((lim.range() - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_urdf_robot_dof() {
        let mut robot = UrdfRobot::new("arm");
        robot.add_link(UrdfLink::new("base"));
        robot.add_link(UrdfLink::new("link1"));
        let mut j = UrdfJoint::revolute("joint1", "base", "link1");
        j.joint_type = JointType::Revolute;
        robot.add_joint(j);
        assert_eq!(robot.dof(), 1);
    }
    #[test]
    fn test_urdf_robot_total_mass() {
        let mut robot = UrdfRobot::new("test");
        let mut l1 = UrdfLink::new("l1");
        l1.inertial.mass = 2.5;
        let mut l2 = UrdfLink::new("l2");
        l2.inertial.mass = 1.5;
        robot.add_link(l1);
        robot.add_link(l2);
        assert!((robot.total_mass() - 4.0).abs() < 1e-12);
    }
    #[test]
    fn test_urdf_robot_to_text_round_trip() {
        let mut robot = UrdfRobot::new("my_robot");
        let mut link = UrdfLink::new("base_link");
        link.inertial.mass = 3.0;
        robot.add_link(link);
        let text = robot.to_text();
        assert!(text.contains("my_robot"));
        assert!(text.contains("base_link"));
    }
    #[test]
    fn test_joint_state_position_of() {
        let msg = JointStateMsg::new(
            1.0,
            vec!["j0".to_string(), "j1".to_string()],
            vec![0.1, 0.2],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        assert!((msg.position_of("j1").unwrap() - 0.2).abs() < 1e-12);
        assert!(msg.position_of("j99").is_none());
    }
    #[test]
    fn test_joint_state_unit_kinetic_energy() {
        let msg = JointStateMsg::new(0.0, vec!["a".to_string()], vec![0.0], vec![2.0], vec![0.0]);
        assert!((msg.unit_kinetic_energy() - 2.0).abs() < 1e-12);
    }
    #[test]
    fn test_joint_state_to_bytes_length() {
        let msg = JointStateMsg::new(
            0.0,
            vec!["a".to_string(), "b".to_string()],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
            vec![0.0, 0.0],
        );
        let bytes = msg.to_bytes();
        assert_eq!(bytes.len(), 8 + 2 * 24);
    }
    #[test]
    fn test_joint_state_csv_line_contains_timestamp() {
        let msg = JointStateMsg::new(
            42.0,
            vec!["q0".to_string()],
            vec![1.0],
            vec![0.5],
            vec![0.1],
        );
        let line = msg.to_csv_line();
        assert!(line.starts_with("42.000000"));
    }
    #[test]
    fn test_trajectory_duration() {
        let mut traj = RobotTrajectory::new(vec!["j0".to_string()]);
        traj.push(JointWaypoint::new(0.0, vec![0.0]));
        traj.push(JointWaypoint::new(5.0, vec![1.0]));
        assert!((traj.duration() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_trajectory_interpolate_endpoints() {
        let mut traj = RobotTrajectory::new(vec!["j0".to_string()]);
        traj.push(JointWaypoint::new(0.0, vec![0.0]));
        traj.push(JointWaypoint::new(1.0, vec![1.0]));
        let q0 = traj.interpolate(0.0).unwrap();
        let q1 = traj.interpolate(1.0).unwrap();
        assert!((q0[0] - 0.0).abs() < 1e-12);
        assert!((q1[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_trajectory_interpolate_midpoint_monotone() {
        let mut traj = RobotTrajectory::new(vec!["j0".to_string()]);
        traj.push(JointWaypoint::new(0.0, vec![0.0]));
        traj.push(JointWaypoint::new(1.0, vec![1.0]));
        let qm = traj.interpolate(0.5).unwrap()[0];
        assert!(qm > 0.0 && qm < 1.0, "midpoint {} out of range", qm);
    }
    #[test]
    fn test_trajectory_sample_uniform_count() {
        let mut traj = RobotTrajectory::new(vec!["j0".to_string()]);
        traj.push(JointWaypoint::new(0.0, vec![0.0]));
        traj.push(JointWaypoint::new(2.0, vec![2.0]));
        let samples = traj.sample_uniform(11);
        assert_eq!(samples.len(), 11);
    }
    #[test]
    fn test_trajectory_csv_header() {
        let traj = RobotTrajectory::new(vec!["shoulder".to_string(), "elbow".to_string()]);
        let csv = traj.to_csv();
        assert!(csv.starts_with("time,shoulder,elbow"));
    }
    #[test]
    fn test_gripper_openness() {
        let g = GripperState {
            position: 0.0425,
            ..Default::default()
        };
        let o = g.openness();
        assert!((o - 0.5).abs() < 0.01);
    }
    #[test]
    fn test_gripper_is_closed() {
        let mut g = GripperState {
            position: 0.0,
            ..Default::default()
        };
        assert!(g.is_closed());
        g.position = 0.05;
        assert!(!g.is_closed());
    }
    #[test]
    fn test_gripper_round_trip_bytes() {
        let g = GripperState {
            position: 0.04,
            force: 30.0,
            max_width: 0.085,
            object_detected: true,
            temperature: 28.5,
        };
        let bytes = g.to_bytes();
        let g2 = GripperState::from_bytes(&bytes);
        assert!((g2.position - g.position).abs() < 1e-12);
        assert!((g2.force - g.force).abs() < 1e-12);
        assert!(g2.object_detected);
    }
    #[test]
    fn test_wrench_force_magnitude() {
        let w = Wrench::new([3.0, 4.0, 0.0], [0.0; 3]);
        assert!((w.force_magnitude() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_wrench_add() {
        let a = Wrench::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let b = Wrench::new([2.0, 0.0, 0.0], [0.0, 2.0, 0.0]);
        let c = a.add(&b);
        assert!((c.force[0] - 3.0).abs() < 1e-12);
        assert!((c.torque[1] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_wrench_to_bytes_length() {
        let w = Wrench::default();
        assert_eq!(w.to_bytes().len(), 48);
    }
    #[test]
    fn test_ros_msg_round_trip() {
        let payload = vec![1u8, 2, 3, 4, 5];
        let msg = RosMsg::new(MsgType::Generic, 42, payload.clone());
        let bytes = msg.to_bytes();
        let msg2 = RosMsg::from_bytes(&bytes).expect("should parse");
        assert_eq!(msg2.msg_type, MsgType::Generic);
        assert_eq!(msg2.seq, 42);
        assert_eq!(msg2.payload, payload);
    }
    #[test]
    fn test_ros_msg_bad_magic_rejected() {
        let mut bytes = RosMsg::new(MsgType::Generic, 0, vec![]).to_bytes();
        bytes[0] = 0xFF;
        assert!(RosMsg::from_bytes(&bytes).is_none());
    }
    #[test]
    fn test_ros_msg_bad_checksum_rejected() {
        let mut bytes = RosMsg::new(MsgType::Generic, 0, vec![0xAA, 0xBB]).to_bytes();
        let last = bytes.len() - 1;
        bytes[last] ^= 0xFF;
        assert!(RosMsg::from_bytes(&bytes).is_none());
    }
    #[test]
    fn test_workspace_grid_set_and_get() {
        let mut grid = WorkspaceGrid::new(10, 10, 10, 0.1, [0.0; 3]);
        let pos = [0.05, 0.05, 0.05];
        grid.set_at(pos, 0.9);
        let val = grid.get_at(pos).unwrap();
        assert!((val - 0.9).abs() < 1e-6);
    }
    #[test]
    fn test_workspace_grid_out_of_bounds() {
        let grid = WorkspaceGrid::new(5, 5, 5, 0.1, [0.0; 3]);
        assert!(grid.get_at([10.0, 10.0, 10.0]).is_none());
    }
    #[test]
    fn test_workspace_grid_reachability_fraction() {
        let mut grid = WorkspaceGrid::new(2, 2, 2, 1.0, [0.0; 3]);
        for v in grid.data.iter_mut().take(4) {
            *v = 1.0;
        }
        let frac = grid.reachability_fraction(0.5);
        assert!((frac - 0.5).abs() < 1e-6);
    }
    #[test]
    fn test_ik_solution_total_error() {
        let sol = IkSolution::success(vec![0.0; 6], 0.003, 0.004, 100);
        assert!((sol.total_error() - 0.005).abs() < 1e-9);
    }
    #[test]
    fn test_ik_solution_failure_is_not_success() {
        let sol = IkSolution::failure(50);
        assert!(!sol.success);
        assert!(sol.position_error.is_infinite());
    }
    #[test]
    fn test_ik_solution_log_success_rate() {
        let mut log = IkSolutionLog::new();
        log.push(IkSolution::success(vec![], 0.0, 0.0, 1));
        log.push(IkSolution::success(vec![], 0.0, 0.0, 1));
        log.push(IkSolution::failure(10));
        let rate = log.success_rate();
        assert!((rate - 2.0 / 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_camera_project_origin() {
        let cam = CameraIntrinsics::new(500.0, 500.0, 320.0, 240.0, 640, 480);
        let px = cam.project([0.0, 0.0, 1.0]).unwrap();
        assert!((px[0] - 320.0).abs() < 1e-9);
        assert!((px[1] - 240.0).abs() < 1e-9);
    }
    #[test]
    fn test_camera_project_negative_z_returns_none() {
        let cam = CameraIntrinsics::new(500.0, 500.0, 320.0, 240.0, 640, 480);
        assert!(cam.project([0.0, 0.0, -1.0]).is_none());
    }
    #[test]
    fn test_camera_fov_reasonable() {
        let cam = CameraIntrinsics::new(500.0, 500.0, 320.0, 240.0, 640, 480);
        let fov = cam.fov_x();
        assert!(fov > 0.5 && fov < 2.0, "fov_x = {fov}");
    }
    #[test]
    fn test_point_cloud_centroid() {
        let mut cloud = RobotPointCloud::new(0.0, "base_link");
        cloud.points.push(ScanPoint::from_position([1.0, 0.0, 0.0]));
        cloud
            .points
            .push(ScanPoint::from_position([-1.0, 0.0, 0.0]));
        let c = cloud.centroid().unwrap();
        assert!(c[0].abs() < 1e-12);
    }
    #[test]
    fn test_point_cloud_aabb() {
        let mut cloud = RobotPointCloud::new(0.0, "f");
        cloud.points.push(ScanPoint::from_position([1.0, 2.0, 3.0]));
        cloud
            .points
            .push(ScanPoint::from_position([-1.0, -2.0, -3.0]));
        let (mn, mx) = cloud.aabb().unwrap();
        assert!((mn[0] + 1.0).abs() < 1e-12);
        assert!((mx[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_point_cloud_round_trip_bytes() {
        let mut cloud = RobotPointCloud::new(1.5, "base");
        for i in 0..5 {
            cloud
                .points
                .push(ScanPoint::from_position([i as f64, 0.0, 0.0]));
        }
        let bytes = cloud.to_bytes();
        let cloud2 = RobotPointCloud::from_bytes(&bytes).unwrap();
        assert_eq!(cloud2.points.len(), 5);
        assert!((cloud2.points[2].position[0] - 2.0).abs() < 1e-4);
    }
    #[test]
    fn test_point_cloud_downsample() {
        let mut cloud = RobotPointCloud::new(0.0, "f");
        for i in 0..10 {
            cloud
                .points
                .push(ScanPoint::from_position([i as f64, 0.0, 0.0]));
        }
        let down = cloud.downsample(2);
        assert_eq!(down.len(), 5);
    }
    #[test]
    fn test_point_cloud_transform() {
        let mut cloud = RobotPointCloud::new(0.0, "f");
        cloud.points.push(ScanPoint::from_position([0.0, 0.0, 0.0]));
        let t = mat4_translation([5.0, 0.0, 0.0]);
        cloud.transform_in_place(&t);
        assert!((cloud.points[0].position[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_quat_identity() {
        let m = quat_to_mat4([1.0, 0.0, 0.0, 0.0]);
        let id = mat4_identity();
        for i in 0..16 {
            assert!((m[i] - id[i]).abs() < 1e-12, "diff at {i}");
        }
    }
    #[test]
    fn test_quat_180_deg_around_z() {
        use std::f64::consts::FRAC_PI_2;
        let m = quat_to_mat4([0.0, 0.0, 0.0, 1.0]);
        let p = apply_transform(&m, [1.0, 0.0, 0.0]);
        assert!((p[0] + 1.0).abs() < 1e-10, "x should flip: {}", p[0]);
        let _ = FRAC_PI_2;
    }
    #[test]
    fn test_rpy_zero_is_identity() {
        let m = rpy_to_mat4([0.0, 0.0, 0.0]);
        let id = mat4_identity();
        for i in 0..16 {
            assert!((m[i] - id[i]).abs() < 1e-12);
        }
    }
}
