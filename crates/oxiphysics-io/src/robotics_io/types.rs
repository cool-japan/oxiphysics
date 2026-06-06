//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;

use super::functions::*;
use super::functions::{Mat4, ROS_MAGIC, Vec3};

/// A ROS-like joint state message.
#[derive(Debug, Clone, Default)]
pub struct JointStateMsg {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// Joint names.
    pub names: Vec<String>,
    /// Joint positions (rad or m).
    pub positions: Vec<f64>,
    /// Joint velocities (rad/s or m/s).
    pub velocities: Vec<f64>,
    /// Joint efforts (N·m or N).
    pub efforts: Vec<f64>,
}
impl JointStateMsg {
    /// Construct a new joint state with the given timestamp and equal-length vectors.
    ///
    /// Panics in debug mode if the vectors have different lengths.
    pub fn new(
        timestamp: f64,
        names: Vec<String>,
        positions: Vec<f64>,
        velocities: Vec<f64>,
        efforts: Vec<f64>,
    ) -> Self {
        debug_assert_eq!(names.len(), positions.len());
        debug_assert_eq!(names.len(), velocities.len());
        debug_assert_eq!(names.len(), efforts.len());
        Self {
            timestamp,
            names,
            positions,
            velocities,
            efforts,
        }
    }
    /// Return the number of joints.
    pub fn len(&self) -> usize {
        self.names.len()
    }
    /// Return `true` if there are no joints.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }
    /// Look up the index of a joint by name.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.names.iter().position(|n| n == name)
    }
    /// Get the position of a joint by name.
    pub fn position_of(&self, name: &str) -> Option<f64> {
        self.index_of(name).map(|i| self.positions[i])
    }
    /// Kinetic energy assuming unit inertia for each joint.
    pub fn unit_kinetic_energy(&self) -> f64 {
        self.velocities.iter().map(|v| 0.5 * v * v).sum()
    }
    /// Serialise to CSV line: `timestamp,name0:pos0:vel0:eff0,...`.
    pub fn to_csv_line(&self) -> String {
        let fields: Vec<String> = self
            .names
            .iter()
            .enumerate()
            .map(|(i, n)| {
                format!(
                    "{}:{:.6}:{:.6}:{:.6}",
                    n, self.positions[i], self.velocities[i], self.efforts[i]
                )
            })
            .collect();
        format!("{:.6},{}", self.timestamp, fields.join(","))
    }
    /// Serialise to a binary buffer (little-endian f64 per value).
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = self.names.len();
        let mut buf = Vec::with_capacity(8 + n * 24);
        buf.extend_from_slice(&(self.timestamp).to_le_bytes());
        for i in 0..n {
            buf.extend_from_slice(&self.positions[i].to_le_bytes());
            buf.extend_from_slice(&self.velocities[i].to_le_bytes());
            buf.extend_from_slice(&self.efforts[i].to_le_bytes());
        }
        buf
    }
}
/// A 3-D reachability grid for a robot workspace.
///
/// Each voxel stores a reachability score in `[0, 1]`: `1` = fully reachable,
/// `0` = unreachable or not evaluated.
#[derive(Debug, Clone)]
pub struct WorkspaceGrid {
    /// Number of cells along X.
    pub nx: usize,
    /// Number of cells along Y.
    pub ny: usize,
    /// Number of cells along Z.
    pub nz: usize,
    /// Grid cell size in metres.
    pub cell_size: f64,
    /// World-space origin `[x0, y0, z0]` of the (0,0,0) voxel.
    pub origin: Vec3,
    /// Reachability scores, linearised as `data[ix + nx*(iy + ny*iz)]`.
    pub data: Vec<f32>,
}
impl WorkspaceGrid {
    /// Construct an all-zero workspace grid.
    pub fn new(nx: usize, ny: usize, nz: usize, cell_size: f64, origin: Vec3) -> Self {
        Self {
            nx,
            ny,
            nz,
            cell_size,
            origin,
            data: vec![0.0; nx * ny * nz],
        }
    }
    /// Convert a world position to a voxel index `(ix, iy, iz)`.
    pub fn world_to_index(&self, pos: Vec3) -> Option<(usize, usize, usize)> {
        let ix = ((pos[0] - self.origin[0]) / self.cell_size) as isize;
        let iy = ((pos[1] - self.origin[1]) / self.cell_size) as isize;
        let iz = ((pos[2] - self.origin[2]) / self.cell_size) as isize;
        if ix < 0 || iy < 0 || iz < 0 {
            return None;
        }
        let (ix, iy, iz) = (ix as usize, iy as usize, iz as usize);
        if ix >= self.nx || iy >= self.ny || iz >= self.nz {
            return None;
        }
        Some((ix, iy, iz))
    }
    /// Get the reachability score at world position `pos`.
    pub fn get_at(&self, pos: Vec3) -> Option<f32> {
        let (ix, iy, iz) = self.world_to_index(pos)?;
        Some(self.data[ix + self.nx * (iy + self.ny * iz)])
    }
    /// Set the reachability score at world position `pos`.
    pub fn set_at(&mut self, pos: Vec3, value: f32) -> bool {
        if let Some((ix, iy, iz)) = self.world_to_index(pos) {
            self.data[ix + self.nx * (iy + self.ny * iz)] = value;
            true
        } else {
            false
        }
    }
    /// Fraction of voxels with reachability above `threshold`.
    pub fn reachability_fraction(&self, threshold: f32) -> f64 {
        let count = self.data.iter().filter(|&&v| v >= threshold).count();
        count as f64 / self.data.len() as f64
    }
    /// Maximum reachability score in the grid.
    pub fn max_reachability(&self) -> f32 {
        self.data.iter().cloned().fold(0.0f32, f32::max)
    }
    /// Serialise to bytes (header + data).
    pub fn to_bytes(&self) -> Vec<u8> {
        let header_size = 4 * 3 + 8 + 8 * 3;
        let mut buf = Vec::with_capacity(header_size + self.data.len() * 4);
        buf.extend_from_slice(&(self.nx as u32).to_le_bytes());
        buf.extend_from_slice(&(self.ny as u32).to_le_bytes());
        buf.extend_from_slice(&(self.nz as u32).to_le_bytes());
        buf.extend_from_slice(&self.cell_size.to_le_bytes());
        for v in &self.origin {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &self.data {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }
}
/// A single joint in a URDF robot description.
#[derive(Debug, Clone)]
pub struct UrdfJoint {
    /// Unique joint name.
    pub name: String,
    /// Joint type.
    pub joint_type: JointType,
    /// Name of the parent link.
    pub parent: String,
    /// Name of the child link.
    pub child: String,
    /// Origin translation `[x,y,z]` relative to parent.
    pub origin_xyz: Vec3,
    /// Origin rotation in roll-pitch-yaw (radians).
    pub origin_rpy: Vec3,
    /// Joint axis (unit vector, default z-axis).
    pub axis: Vec3,
    /// Position/velocity/effort limits.
    pub limits: JointLimits,
    /// Damping coefficient.
    pub damping: f64,
    /// Friction coefficient.
    pub friction: f64,
}
impl UrdfJoint {
    /// Construct a revolute joint with default settings.
    pub fn revolute(
        name: impl Into<String>,
        parent: impl Into<String>,
        child: impl Into<String>,
    ) -> Self {
        Self {
            name: name.into(),
            joint_type: JointType::Revolute,
            parent: parent.into(),
            child: child.into(),
            origin_xyz: [0.0; 3],
            origin_rpy: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            limits: JointLimits::default(),
            damping: 0.1,
            friction: 0.01,
        }
    }
    /// Return `true` if this is a moveable (non-fixed) joint.
    pub fn is_moveable(&self) -> bool {
        !matches!(self.joint_type, JointType::Fixed)
    }
    /// Build the local transform of this joint at zero position.
    pub fn local_transform_zero(&self) -> Mat4 {
        let t = mat4_translation(self.origin_xyz);
        let r = rpy_to_mat4(self.origin_rpy);
        mat4_mul(t, r)
    }
}
/// Result of an inverse kinematics solve.
#[derive(Debug, Clone)]
pub struct IkSolution {
    /// Whether a solution was found.
    pub success: bool,
    /// Joint configuration at the solution (radians or metres).
    pub joint_positions: Vec<f64>,
    /// Residual Cartesian position error (metres).
    pub position_error: f64,
    /// Residual orientation error (radians).
    pub orientation_error: f64,
    /// Number of iterations taken.
    pub iterations: u32,
    /// Name of the IK solver used.
    pub solver_name: String,
}
impl IkSolution {
    /// Construct a successful IK solution.
    pub fn success(
        joint_positions: Vec<f64>,
        position_error: f64,
        orientation_error: f64,
        iterations: u32,
    ) -> Self {
        Self {
            success: true,
            joint_positions,
            position_error,
            orientation_error,
            iterations,
            solver_name: "default".to_string(),
        }
    }
    /// Construct a failed IK solution.
    pub fn failure(iterations: u32) -> Self {
        Self {
            success: false,
            joint_positions: Vec::new(),
            position_error: f64::INFINITY,
            orientation_error: f64::INFINITY,
            iterations,
            solver_name: "default".to_string(),
        }
    }
    /// Total error (Euclidean sum of position and orientation errors).
    pub fn total_error(&self) -> f64 {
        (self.position_error * self.position_error
            + self.orientation_error * self.orientation_error)
            .sqrt()
    }
    /// Serialise to a JSON-like string.
    pub fn to_json_string(&self) -> String {
        let joints: Vec<String> = self
            .joint_positions
            .iter()
            .map(|q| format!("{q:.6}"))
            .collect();
        format!(
            "{{\"success\":{},\"joints\":[{}],\"pos_err\":{:.8},\"ori_err\":{:.8},\"iter\":{}}}",
            self.success,
            joints.join(","),
            self.position_error,
            self.orientation_error,
            self.iterations
        )
    }
}
/// A log of multiple IK solutions, e.g. for benchmarking.
#[derive(Debug, Default, Clone)]
pub struct IkSolutionLog {
    /// Stored solutions.
    pub solutions: Vec<IkSolution>,
}
impl IkSolutionLog {
    /// Create an empty log.
    pub fn new() -> Self {
        Self::default()
    }
    /// Append a solution.
    pub fn push(&mut self, sol: IkSolution) {
        self.solutions.push(sol);
    }
    /// Success rate of all logged solutions.
    pub fn success_rate(&self) -> f64 {
        if self.solutions.is_empty() {
            return 0.0;
        }
        let n = self.solutions.iter().filter(|s| s.success).count();
        n as f64 / self.solutions.len() as f64
    }
    /// Mean position error of successful solutions.
    pub fn mean_position_error(&self) -> f64 {
        let succ: Vec<_> = self.solutions.iter().filter(|s| s.success).collect();
        if succ.is_empty() {
            return f64::NAN;
        }
        succ.iter().map(|s| s.position_error).sum::<f64>() / succ.len() as f64
    }
    /// Mean number of iterations.
    pub fn mean_iterations(&self) -> f64 {
        if self.solutions.is_empty() {
            return 0.0;
        }
        self.solutions
            .iter()
            .map(|s| s.iterations as f64)
            .sum::<f64>()
            / self.solutions.len() as f64
    }
}
/// Camera intrinsic parameters.
#[derive(Debug, Clone, Copy)]
pub struct CameraIntrinsics {
    /// Focal length in x (pixels).
    pub fx: f64,
    /// Focal length in y (pixels).
    pub fy: f64,
    /// Principal point x (pixels).
    pub cx: f64,
    /// Principal point y (pixels).
    pub cy: f64,
    /// Image width in pixels.
    pub width: u32,
    /// Image height in pixels.
    pub height: u32,
}
impl CameraIntrinsics {
    /// Construct standard intrinsics.
    pub fn new(fx: f64, fy: f64, cx: f64, cy: f64, width: u32, height: u32) -> Self {
        Self {
            fx,
            fy,
            cx,
            cy,
            width,
            height,
        }
    }
    /// Field of view in x (radians).
    pub fn fov_x(&self) -> f64 {
        2.0 * (self.width as f64 / (2.0 * self.fx)).atan()
    }
    /// Field of view in y (radians).
    pub fn fov_y(&self) -> f64 {
        2.0 * (self.height as f64 / (2.0 * self.fy)).atan()
    }
    /// Project a 3-D point `[x,y,z]` (camera frame) to pixel `[u,v]`.
    ///
    /// Returns `None` if `z <= 0`.
    pub fn project(&self, point: Vec3) -> Option<[f64; 2]> {
        let z = point[2];
        if z <= 0.0 {
            return None;
        }
        let u = self.fx * point[0] / z + self.cx;
        let v = self.fy * point[1] / z + self.cy;
        Some([u, v])
    }
}
/// A single point in a robot-mounted scanner point cloud.
#[derive(Debug, Clone, Copy)]
pub struct ScanPoint {
    /// 3-D position in robot base frame `[x, y, z]` (metres).
    pub position: Vec3,
    /// Intensity (if available, else `0`).
    pub intensity: f32,
    /// Scan ring index (for multi-layer lidars).
    pub ring: u16,
    /// Return number (1-based, for multi-return lidars).
    pub return_num: u8,
}
impl ScanPoint {
    /// Construct a point from position only.
    pub fn from_position(pos: Vec3) -> Self {
        Self {
            position: pos,
            intensity: 0.0,
            ring: 0,
            return_num: 1,
        }
    }
}
/// Coefficients `[a0, a1, a2, a3]` for a cubic polynomial `a0 + a1*t + a2*t² + a3*t³`.
#[derive(Debug, Clone, Copy)]
pub(super) struct CubicCoeffs {
    pub(super) a0: f64,
    pub(super) a1: f64,
    pub(super) a2: f64,
    pub(super) a3: f64,
}
impl CubicCoeffs {
    fn eval(&self, t: f64) -> f64 {
        self.a0 + self.a1 * t + self.a2 * t * t + self.a3 * t * t * t
    }
}
/// A visual or collision element attached to a URDF link.
#[derive(Debug, Clone)]
pub struct UrdfVisualElement {
    /// Optional name of the visual element.
    pub name: Option<String>,
    /// Origin pose: translation `[x,y,z]` in metres.
    pub origin_xyz: Vec3,
    /// Origin pose: rotation in roll-pitch-yaw (radians).
    pub origin_rpy: Vec3,
    /// Geometry description.
    pub geometry: UrdfGeometry,
    /// Optional material name.
    pub material: Option<String>,
}
impl UrdfVisualElement {
    /// Create a new visual element at the identity origin.
    pub fn new(geometry: UrdfGeometry) -> Self {
        Self {
            name: None,
            origin_xyz: [0.0; 3],
            origin_rpy: [0.0; 3],
            geometry,
            material: None,
        }
    }
    /// Approximate bounding-sphere radius of the geometry.
    pub fn bounding_radius(&self) -> f64 {
        match &self.geometry {
            UrdfGeometry::Box { half_extents } => norm3(*half_extents),
            UrdfGeometry::Sphere { radius } => *radius,
            UrdfGeometry::Cylinder {
                radius,
                half_length,
            } => (radius * radius + half_length * half_length).sqrt(),
            UrdfGeometry::Mesh { scale, .. } => norm3(*scale),
        }
    }
}
/// A stamped wrench measurement (time + value).
#[derive(Debug, Clone, Copy)]
pub struct WrenchStamped {
    /// Timestamp in seconds.
    pub timestamp: f64,
    /// The wrench measurement.
    pub wrench: Wrench,
}
impl WrenchStamped {
    /// Construct a stamped wrench.
    pub fn new(timestamp: f64, wrench: Wrench) -> Self {
        Self { timestamp, wrench }
    }
}
/// Camera-to-robot extrinsic calibration (a rigid-body transform).
#[derive(Debug, Clone)]
pub struct CameraRobotCalibration {
    /// Transform from camera frame to robot base frame.
    pub camera_to_robot: Mat4,
    /// Transform from robot base frame to camera frame.
    pub robot_to_camera: Mat4,
    /// Camera intrinsics.
    pub intrinsics: CameraIntrinsics,
    /// Calibration residual (reprojection error in pixels).
    pub residual_pixels: f64,
    /// Calibration date (ISO 8601 string, informational).
    pub calibration_date: String,
}
impl CameraRobotCalibration {
    /// Construct a calibration with identity transform.
    pub fn identity(intrinsics: CameraIntrinsics) -> Self {
        Self {
            camera_to_robot: mat4_identity(),
            robot_to_camera: mat4_identity(),
            intrinsics,
            residual_pixels: 0.0,
            calibration_date: "2026-01-01".to_string(),
        }
    }
    /// Transform a 3-D point from camera frame to robot base frame.
    pub fn camera_to_robot_point(&self, p: Vec3) -> Vec3 {
        apply_transform(&self.camera_to_robot, p)
    }
    /// Transform a 3-D point from robot base frame to camera frame.
    pub fn robot_to_camera_point(&self, p: Vec3) -> Vec3 {
        apply_transform(&self.robot_to_camera, p)
    }
}
/// A 6-DOF wrench: force `[fx, fy, fz]` (N) and torque `[tx, ty, tz]` (N·m).
#[derive(Debug, Clone, Copy, Default)]
pub struct Wrench {
    /// Force vector `[fx, fy, fz]` in Newtons.
    pub force: Vec3,
    /// Torque vector `[tx, ty, tz]` in Newton-metres.
    pub torque: Vec3,
}
impl Wrench {
    /// Create a wrench from force and torque arrays.
    pub fn new(force: Vec3, torque: Vec3) -> Self {
        Self { force, torque }
    }
    /// Magnitude of the force vector.
    pub fn force_magnitude(&self) -> f64 {
        norm3(self.force)
    }
    /// Magnitude of the torque vector.
    pub fn torque_magnitude(&self) -> f64 {
        norm3(self.torque)
    }
    /// Combine two wrenches by addition.
    pub fn add(&self, other: &Wrench) -> Wrench {
        Wrench {
            force: add3(self.force, other.force),
            torque: add3(self.torque, other.torque),
        }
    }
    /// Scale the wrench by a scalar.
    pub fn scale(&self, s: f64) -> Wrench {
        Wrench {
            force: scale3(self.force, s),
            torque: scale3(self.torque, s),
        }
    }
    /// Serialise to bytes (little-endian f64).
    pub fn to_bytes(&self) -> [u8; 48] {
        let mut buf = [0u8; 48];
        for (i, v) in self.force.iter().chain(self.torque.iter()).enumerate() {
            buf[i * 8..(i + 1) * 8].copy_from_slice(&v.to_le_bytes());
        }
        buf
    }
}
/// Message type identifier used in the simplified ROS-like protocol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MsgType {
    /// Joint state message.
    JointState = 1,
    /// Wrench (force/torque) message.
    Wrench = 2,
    /// Gripper state message.
    GripperState = 3,
    /// Point cloud message.
    PointCloud = 4,
    /// Generic data blob.
    Generic = 255,
}
/// Inertial properties of a URDF link.
#[derive(Debug, Clone)]
pub struct UrdfInertial {
    /// Mass in kilograms.
    pub mass: f64,
    /// Centre of mass position `[x,y,z]`.
    pub com: Vec3,
    /// Inertia tensor diagonal `[ixx, iyy, izz]` (kg·m²).
    pub inertia_diag: Vec3,
    /// Off-diagonal inertia `[ixy, ixz, iyz]` (kg·m²).
    pub inertia_off: Vec3,
}
/// Geometry kind attached to a URDF link's visual or collision element.
#[derive(Debug, Clone, PartialEq)]
pub enum UrdfGeometry {
    /// Box primitive with half-extents `[hx, hy, hz]` in metres.
    Box {
        /// Half-extents (x, y, z) in metres.
        half_extents: Vec3,
    },
    /// Sphere primitive.
    Sphere {
        /// Radius in metres.
        radius: f64,
    },
    /// Cylinder primitive.
    Cylinder {
        /// Radius in metres.
        radius: f64,
        /// Half-length in metres.
        half_length: f64,
    },
    /// Mesh reference (path to file, optional scale).
    Mesh {
        /// Path to the mesh file.
        filename: String,
        /// Scale factor applied to the mesh.
        scale: Vec3,
    },
}
/// Joint type as found in URDF.
#[derive(Debug, Clone, PartialEq)]
pub enum JointType {
    /// Fixed joint (no movement).
    Fixed,
    /// Revolute (rotational, with limits).
    Revolute,
    /// Prismatic (translational, with limits).
    Prismatic,
    /// Continuous (revolute, no limits).
    Continuous,
    /// Planar (movement in a plane).
    Planar,
    /// Floating (6-DOF).
    Floating,
}
/// A single link in a URDF robot description.
#[derive(Debug, Clone)]
pub struct UrdfLink {
    /// Unique link name.
    pub name: String,
    /// Inertial properties.
    pub inertial: UrdfInertial,
    /// Visual geometry elements.
    pub visuals: Vec<UrdfVisualElement>,
    /// Collision geometry elements.
    pub collisions: Vec<UrdfVisualElement>,
}
impl UrdfLink {
    /// Construct a new link with the given name and default inertial.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            inertial: UrdfInertial::default(),
            visuals: Vec::new(),
            collisions: Vec::new(),
        }
    }
    /// Total approximate volume of all visual bounding spheres.
    pub fn approximate_visual_volume(&self) -> f64 {
        self.visuals
            .iter()
            .map(|v| {
                let r = v.bounding_radius();
                (4.0 / 3.0) * std::f64::consts::PI * r * r * r
            })
            .sum()
    }
}
/// A single waypoint in a robot trajectory.
#[derive(Debug, Clone)]
pub struct JointWaypoint {
    /// Time from the start of the trajectory (seconds).
    pub time: f64,
    /// Joint positions at this waypoint.
    pub positions: Vec<f64>,
    /// Optional joint velocities.
    pub velocities: Option<Vec<f64>>,
    /// Optional joint accelerations.
    pub accelerations: Option<Vec<f64>>,
}
impl JointWaypoint {
    /// Construct a waypoint with only positions specified.
    pub fn new(time: f64, positions: Vec<f64>) -> Self {
        Self {
            time,
            positions,
            velocities: None,
            accelerations: None,
        }
    }
}
/// A framed message in the simplified binary protocol.
///
/// Wire format: `[magic:4][msg_type:1][seq:4][payload_len:4][payload:N][crc:4]`
#[derive(Debug, Clone)]
pub struct RosMsg {
    /// Message type.
    pub msg_type: MsgType,
    /// Sequence number.
    pub seq: u32,
    /// Raw payload bytes.
    pub payload: Vec<u8>,
}
impl RosMsg {
    /// Construct a new message.
    pub fn new(msg_type: MsgType, seq: u32, payload: Vec<u8>) -> Self {
        Self {
            msg_type,
            seq,
            payload,
        }
    }
    /// Simple 32-bit checksum (sum of all payload bytes, modulo 2³²).
    fn checksum(data: &[u8]) -> u32 {
        data.iter().fold(0u32, |acc, &b| acc.wrapping_add(b as u32))
    }
    /// Serialise the message to wire bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let len = self.payload.len() as u32;
        let crc = Self::checksum(&self.payload);
        let mut buf = Vec::with_capacity(17 + self.payload.len());
        buf.extend_from_slice(&ROS_MAGIC);
        buf.push(self.msg_type as u8);
        buf.extend_from_slice(&self.seq.to_le_bytes());
        buf.extend_from_slice(&len.to_le_bytes());
        buf.extend_from_slice(&self.payload);
        buf.extend_from_slice(&crc.to_le_bytes());
        buf
    }
    /// Deserialise from wire bytes.  Returns `None` if the frame is invalid.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 17 {
            return None;
        }
        if data[0..4] != ROS_MAGIC {
            return None;
        }
        let msg_type = MsgType::try_from(data[4]).ok()?;
        let seq = u32::from_le_bytes(data[5..9].try_into().ok()?);
        let len = u32::from_le_bytes(data[9..13].try_into().ok()?) as usize;
        if data.len() < 13 + len + 4 {
            return None;
        }
        let payload = data[13..13 + len].to_vec();
        let crc_stored = u32::from_le_bytes(data[13 + len..17 + len].try_into().ok()?);
        let crc_computed = Self::checksum(&payload);
        if crc_stored != crc_computed {
            return None;
        }
        Some(Self {
            msg_type,
            seq,
            payload,
        })
    }
}
/// A point cloud from a robot-mounted scanner.
#[derive(Debug, Default, Clone)]
pub struct RobotPointCloud {
    /// Timestamp in seconds (scan start).
    pub timestamp: f64,
    /// Points in robot base frame.
    pub points: Vec<ScanPoint>,
    /// Frame ID (e.g. "base_link").
    pub frame_id: String,
    /// Scanner model name.
    pub scanner_model: String,
}
impl RobotPointCloud {
    /// Create an empty point cloud.
    pub fn new(timestamp: f64, frame_id: impl Into<String>) -> Self {
        Self {
            timestamp,
            points: Vec::new(),
            frame_id: frame_id.into(),
            scanner_model: "unknown".to_string(),
        }
    }
    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }
    /// Return `true` if the cloud is empty.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// Axis-aligned bounding box of the cloud: `(min, max)`.
    pub fn aabb(&self) -> Option<(Vec3, Vec3)> {
        if self.points.is_empty() {
            return None;
        }
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for p in &self.points {
            for k in 0..3 {
                if p.position[k] < mn[k] {
                    mn[k] = p.position[k];
                }
                if p.position[k] > mx[k] {
                    mx[k] = p.position[k];
                }
            }
        }
        Some((mn, mx))
    }
    /// Mean point position.
    pub fn centroid(&self) -> Option<Vec3> {
        if self.points.is_empty() {
            return None;
        }
        let n = self.points.len() as f64;
        let mut sum = [0.0; 3];
        for p in &self.points {
            sum[0] += p.position[0];
            sum[1] += p.position[1];
            sum[2] += p.position[2];
        }
        Some([sum[0] / n, sum[1] / n, sum[2] / n])
    }
    /// Apply an in-place transform to all points.
    pub fn transform_in_place(&mut self, m: &Mat4) {
        for p in &mut self.points {
            p.position = apply_transform(m, p.position);
        }
    }
    /// Downsample by keeping every `stride`-th point.
    pub fn downsample(&self, stride: usize) -> Self {
        let stride = stride.max(1);
        let points: Vec<ScanPoint> = self
            .points
            .iter()
            .enumerate()
            .filter(|(i, _)| i % stride == 0)
            .map(|(_, p)| *p)
            .collect();
        Self {
            timestamp: self.timestamp,
            points,
            frame_id: self.frame_id.clone(),
            scanner_model: self.scanner_model.clone(),
        }
    }
    /// Serialise to a minimal binary PCD-inspired format.
    ///
    /// Format: `[n_points:u32][point:4×f32]×n`  (x,y,z,intensity as f32).
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = self.points.len();
        let mut buf = Vec::with_capacity(4 + n * 16);
        buf.extend_from_slice(&(n as u32).to_le_bytes());
        for p in &self.points {
            for k in 0..3 {
                buf.extend_from_slice(&(p.position[k] as f32).to_le_bytes());
            }
            buf.extend_from_slice(&p.intensity.to_le_bytes());
        }
        buf
    }
    /// Deserialise from the binary format produced by `to_bytes`.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 4 {
            return None;
        }
        let n = u32::from_le_bytes(data[0..4].try_into().ok()?) as usize;
        if data.len() < 4 + n * 16 {
            return None;
        }
        let mut points = Vec::with_capacity(n);
        for i in 0..n {
            let base = 4 + i * 16;
            let read_f32 = |b: &[u8]| f32::from_le_bytes(b.try_into().unwrap_or([0u8; 4]));
            let x = read_f32(&data[base..base + 4]) as f64;
            let y = read_f32(&data[base + 4..base + 8]) as f64;
            let z = read_f32(&data[base + 8..base + 12]) as f64;
            let intensity = read_f32(&data[base + 12..base + 16]);
            points.push(ScanPoint {
                position: [x, y, z],
                intensity,
                ring: 0,
                return_num: 1,
            });
        }
        Some(Self {
            timestamp: 0.0,
            points,
            frame_id: String::new(),
            scanner_model: String::new(),
        })
    }
}
/// Limits for a URDF joint.
#[derive(Debug, Clone)]
pub struct JointLimits {
    /// Lower position limit (rad or m).
    pub lower: f64,
    /// Upper position limit (rad or m).
    pub upper: f64,
    /// Maximum effort (N or N·m).
    pub effort: f64,
    /// Maximum velocity (rad/s or m/s).
    pub velocity: f64,
}
impl JointLimits {
    /// Return `true` if `q` lies within `[lower, upper]`.
    pub fn in_range(&self, q: f64) -> bool {
        q >= self.lower && q <= self.upper
    }
    /// Clamp `q` to `[lower, upper]`.
    pub fn clamp(&self, q: f64) -> f64 {
        q.clamp(self.lower, self.upper)
    }
    /// Range of motion.
    pub fn range(&self) -> f64 {
        self.upper - self.lower
    }
}
/// State of a parallel-jaw gripper.
#[derive(Debug, Clone, Copy)]
pub struct GripperState {
    /// Current jaw opening width in metres (`0` = fully closed, `max_width` = fully open).
    pub position: f64,
    /// Current gripping force in Newtons.
    pub force: f64,
    /// Maximum jaw opening width in metres.
    pub max_width: f64,
    /// Whether an object is detected in the gripper.
    pub object_detected: bool,
    /// Gripper temperature (°C), if available.
    pub temperature: f64,
}
impl GripperState {
    /// Normalised opening (`0` = closed, `1` = fully open).
    pub fn openness(&self) -> f64 {
        if self.max_width < 1e-9 {
            0.0
        } else {
            (self.position / self.max_width).clamp(0.0, 1.0)
        }
    }
    /// Return `true` if the gripper is effectively closed (within 1 mm).
    pub fn is_closed(&self) -> bool {
        self.position < 1e-3
    }
    /// Serialise to bytes.
    pub fn to_bytes(&self) -> [u8; 40] {
        let mut buf = [0u8; 40];
        buf[0..8].copy_from_slice(&self.position.to_le_bytes());
        buf[8..16].copy_from_slice(&self.force.to_le_bytes());
        buf[16..24].copy_from_slice(&self.max_width.to_le_bytes());
        buf[24..32].copy_from_slice(&self.temperature.to_le_bytes());
        buf[32] = self.object_detected as u8;
        buf
    }
    /// Deserialise from bytes.
    pub fn from_bytes(b: &[u8; 40]) -> Self {
        let read = |s: &[u8]| f64::from_le_bytes(s.try_into().unwrap_or([0u8; 8]));
        Self {
            position: read(&b[0..8]),
            force: read(&b[8..16]),
            max_width: read(&b[16..24]),
            temperature: read(&b[24..32]),
            object_detected: b[32] != 0,
        }
    }
}
/// A complete robot trajectory consisting of [`JointWaypoint`]s.
#[derive(Debug, Default, Clone)]
pub struct RobotTrajectory {
    /// Ordered waypoints (should be sorted by time).
    pub waypoints: Vec<JointWaypoint>,
    /// Joint names (informational).
    pub joint_names: Vec<String>,
}
impl RobotTrajectory {
    /// Create an empty trajectory for `joint_names`.
    pub fn new(joint_names: Vec<String>) -> Self {
        Self {
            waypoints: Vec::new(),
            joint_names,
        }
    }
    /// Append a waypoint.
    pub fn push(&mut self, wp: JointWaypoint) {
        self.waypoints.push(wp);
    }
    /// Total duration of the trajectory (last waypoint time).
    pub fn duration(&self) -> f64 {
        self.waypoints.last().map(|w| w.time).unwrap_or(0.0)
    }
    /// Number of degrees of freedom (from first waypoint).
    pub fn dof(&self) -> usize {
        self.waypoints
            .first()
            .map(|w| w.positions.len())
            .unwrap_or(0)
    }
    /// Cubic spline interpolation at time `t`.
    ///
    /// Uses natural boundary conditions (zero velocity at endpoints if not specified).
    /// Returns `None` if the trajectory has fewer than 2 waypoints.
    pub fn interpolate(&self, t: f64) -> Option<Vec<f64>> {
        let wps = &self.waypoints;
        if wps.len() < 2 {
            return None;
        }
        let dof = self.dof();
        let t_clamped = t.clamp(wps[0].time, wps[wps.len() - 1].time);
        let seg = wps
            .windows(2)
            .enumerate()
            .find(|(_, pair)| t_clamped >= pair[0].time && t_clamped <= pair[1].time)
            .map(|(i, _)| i)
            .unwrap_or(wps.len() - 2);
        let t0 = wps[seg].time;
        let t1 = wps[seg + 1].time;
        let dt = (t1 - t0).max(1e-12);
        let u = (t_clamped - t0) / dt;
        let mut result = vec![0.0; dof];
        for j in 0..dof {
            let p0 = wps[seg].positions[j];
            let p1 = wps[seg + 1].positions[j];
            let v0 = wps[seg]
                .velocities
                .as_ref()
                .map(|v| v[j] * dt)
                .unwrap_or(0.0);
            let v1 = wps[seg + 1]
                .velocities
                .as_ref()
                .map(|v| v[j] * dt)
                .unwrap_or(0.0);
            let c = CubicCoeffs {
                a0: p0,
                a1: v0,
                a2: 3.0 * (p1 - p0) - 2.0 * v0 - v1,
                a3: -2.0 * (p1 - p0) + v0 + v1,
            };
            result[j] = c.eval(u);
        }
        Some(result)
    }
    /// Sample the trajectory at `n_samples` uniformly spaced times.
    pub fn sample_uniform(&self, n_samples: usize) -> Vec<(f64, Vec<f64>)> {
        if n_samples == 0 || self.waypoints.is_empty() {
            return Vec::new();
        }
        let dur = self.duration();
        (0..n_samples)
            .filter_map(|i| {
                let t = if n_samples == 1 {
                    0.0
                } else {
                    i as f64 * dur / (n_samples - 1) as f64
                };
                self.interpolate(t).map(|q| (t, q))
            })
            .collect()
    }
    /// Serialise to CSV with header.
    pub fn to_csv(&self) -> String {
        let header: Vec<String> = std::iter::once("time".to_string())
            .chain(self.joint_names.iter().cloned())
            .collect();
        let mut lines = vec![header.join(",")];
        for wp in &self.waypoints {
            let row: Vec<String> = std::iter::once(format!("{:.6}", wp.time))
                .chain(wp.positions.iter().map(|p| format!("{:.6}", p)))
                .collect();
            lines.push(row.join(","));
        }
        lines.join("\n")
    }
}
/// A simplified URDF robot model: a collection of links and joints.
#[derive(Debug, Default, Clone)]
pub struct UrdfRobot {
    /// Robot name.
    pub name: String,
    /// All links, keyed by name.
    pub links: HashMap<String, UrdfLink>,
    /// All joints, keyed by name.
    pub joints: HashMap<String, UrdfJoint>,
    /// Ordered list of joint names (insertion order).
    pub joint_order: Vec<String>,
}
impl UrdfRobot {
    /// Create an empty robot with the given name.
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            ..Default::default()
        }
    }
    /// Add a link.  Overwrites if a link with the same name exists.
    pub fn add_link(&mut self, link: UrdfLink) {
        self.links.insert(link.name.clone(), link);
    }
    /// Add a joint.  Overwrites if a joint with the same name exists.
    pub fn add_joint(&mut self, joint: UrdfJoint) {
        if !self.joint_order.contains(&joint.name) {
            self.joint_order.push(joint.name.clone());
        }
        self.joints.insert(joint.name.clone(), joint);
    }
    /// Return a list of moveable joint names in insertion order.
    pub fn moveable_joint_names(&self) -> Vec<&str> {
        self.joint_order
            .iter()
            .filter(|n| {
                self.joints
                    .get(n.as_str())
                    .map(|j| j.is_moveable())
                    .unwrap_or(false)
            })
            .map(|s| s.as_str())
            .collect()
    }
    /// Total number of degrees of freedom (moveable joints).
    pub fn dof(&self) -> usize {
        self.moveable_joint_names().len()
    }
    /// Total mass of the robot (sum of all link masses).
    pub fn total_mass(&self) -> f64 {
        self.links.values().map(|l| l.inertial.mass).sum()
    }
    /// Serialise the robot description to a human-readable text form.
    pub fn to_text(&self) -> String {
        let mut out = format!("robot name=\"{}\"\n", self.name);
        for (name, link) in &self.links {
            out += &format!("  link \"{name}\" mass={:.4} kg\n", link.inertial.mass);
        }
        for name in &self.joint_order {
            if let Some(j) = self.joints.get(name) {
                out += &format!(
                    "  joint \"{}\" type={} parent=\"{}\" child=\"{}\"\n",
                    j.name, j.joint_type, j.parent, j.child
                );
            }
        }
        out
    }
    /// Parse a minimal URDF-like text representation produced by [`UrdfRobot::to_text`].
    ///
    /// This is intentionally a simplified round-trip parser, not a full XML URDF parser.
    pub fn from_text(text: &str) -> Option<Self> {
        let mut robot = UrdfRobot::default();
        for line in text.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix("robot name=\"") {
                robot.name = rest.trim_end_matches('"').to_string();
            } else if let Some(rest) = line.strip_prefix("link \"") {
                let parts: Vec<&str> = rest.splitn(2, '"').collect();
                if let Some(name) = parts.first() {
                    let mass_str = rest
                        .find("mass=")
                        .map(|i| &rest[i + 5..])
                        .and_then(|s| s.split_whitespace().next())
                        .and_then(|s| s.trim_end_matches(" kg").parse::<f64>().ok())
                        .unwrap_or(1.0);
                    let mut link = UrdfLink::new(*name);
                    link.inertial.mass = mass_str;
                    robot.add_link(link);
                }
            }
        }
        if robot.name.is_empty() {
            None
        } else {
            Some(robot)
        }
    }
}
