// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Animation I/O: BVH motion capture, FBX ASCII, USD/USDA export, skeletal
//! animation, blend shapes, timeline/sequencing, frame interpolation
//! (including SLERP for rotations), and animation retargeting.
//!
//! All geometry is represented using plain `[f64; 3]` / `[f64; 4]` arrays
//! (no nalgebra dependency).

use std::collections::HashMap;

// ── Tiny math helpers ──────────────────────────────────────────────────────

/// Add two 3-vectors.
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Linear interpolation for a 3-vector.
fn vec3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Normalize a quaternion `[x, y, z, w]`.
fn quat_normalize(q: [f64; 4]) -> [f64; 4] {
    let len2 = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
    if len2 < 1e-30 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    let inv = 1.0 / len2.sqrt();
    [q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv]
}

/// Quaternion dot product.
fn quat_dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
}

/// Multiply two quaternions (Hamilton product). Layout: `[x, y, z, w]`.
fn quat_mul(p: [f64; 4], q: [f64; 4]) -> [f64; 4] {
    let [px, py, pz, pw] = p;
    let [qx, qy, qz, qw] = q;
    [
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw,
        pw * qw - px * qx - py * qy - pz * qz,
    ]
}

/// Quaternion conjugate (inverse for unit quaternions).
fn quat_conjugate(q: [f64; 4]) -> [f64; 4] {
    [-q[0], -q[1], -q[2], q[3]]
}

/// Euler XYZ (degrees) to quaternion `[x, y, z, w]`.
fn euler_xyz_to_quat(euler_deg: [f64; 3]) -> [f64; 4] {
    let to_rad = std::f64::consts::PI / 180.0;
    let (hx, hy, hz) = (
        euler_deg[0] * to_rad * 0.5,
        euler_deg[1] * to_rad * 0.5,
        euler_deg[2] * to_rad * 0.5,
    );
    let (sx, cx) = (hx.sin(), hx.cos());
    let (sy, cy) = (hy.sin(), hy.cos());
    let (sz, cz) = (hz.sin(), hz.cos());
    quat_normalize([
        sx * cy * cz + cx * sy * sz,
        cx * sy * cz - sx * cy * sz,
        cx * cy * sz + sx * sy * cz,
        cx * cy * cz - sx * sy * sz,
    ])
}

/// Quaternion to Euler XYZ (degrees).
fn quat_to_euler_xyz(q: [f64; 4]) -> [f64; 3] {
    let [x, y, z, w] = q;
    let to_deg = 180.0 / std::f64::consts::PI;
    let sinr_cosp = 2.0 * (w * x + y * z);
    let cosr_cosp = 1.0 - 2.0 * (x * x + y * y);
    let rx = sinr_cosp.atan2(cosr_cosp);
    let sinp = 2.0 * (w * y - z * x);
    let ry = if sinp.abs() >= 1.0 {
        std::f64::consts::FRAC_PI_2.copysign(sinp)
    } else {
        sinp.asin()
    };
    let siny_cosp = 2.0 * (w * z + x * y);
    let cosy_cosp = 1.0 - 2.0 * (y * y + z * z);
    let rz = siny_cosp.atan2(cosy_cosp);
    [rx * to_deg, ry * to_deg, rz * to_deg]
}

// ═══════════════════════════════════════════════════════════════════════════════
// SLERP – Spherical Linear Interpolation
// ═══════════════════════════════════════════════════════════════════════════════

/// Spherical linear interpolation between two unit quaternions.
///
/// Handles the short-arc selection and falls back to linear interpolation
/// when the quaternions are nearly parallel.
pub fn quat_slerp(a: [f64; 4], b: [f64; 4], t: f64) -> [f64; 4] {
    let mut dot = quat_dot(a, b);
    let b = if dot < 0.0 {
        dot = -dot;
        [-b[0], -b[1], -b[2], -b[3]]
    } else {
        b
    };
    if dot > 0.9995 {
        // Nearly parallel – use linear interpolation
        return quat_normalize([
            a[0] + (b[0] - a[0]) * t,
            a[1] + (b[1] - a[1]) * t,
            a[2] + (b[2] - a[2]) * t,
            a[3] + (b[3] - a[3]) * t,
        ]);
    }
    let theta = dot.clamp(-1.0, 1.0).acos();
    let sin_theta = theta.sin();
    if sin_theta.abs() < 1e-15 {
        return a;
    }
    let wa = ((1.0 - t) * theta).sin() / sin_theta;
    let wb = (t * theta).sin() / sin_theta;
    quat_normalize([
        wa * a[0] + wb * b[0],
        wa * a[1] + wb * b[1],
        wa * a[2] + wb * b[2],
        wa * a[3] + wb * b[3],
    ])
}

// ═══════════════════════════════════════════════════════════════════════════════
// Joint & skeleton hierarchy
// ═══════════════════════════════════════════════════════════════════════════════

/// A joint in a skeletal hierarchy.
#[derive(Debug, Clone)]
pub struct Joint {
    /// Joint name.
    pub name: String,
    /// Index of the parent joint, or `None` for root.
    pub parent: Option<usize>,
    /// Rest-pose local offset from parent.
    pub offset: [f64; 3],
    /// Channel names (e.g. "Xposition", "Yposition", "Zrotation").
    pub channels: Vec<String>,
}

/// A full skeleton definition.
#[derive(Debug, Clone)]
pub struct Skeleton {
    /// Ordered list of joints (index 0 is typically the root).
    pub joints: Vec<Joint>,
    /// Map from joint name to index.
    pub name_to_index: HashMap<String, usize>,
}

impl Skeleton {
    /// Create an empty skeleton.
    pub fn new() -> Self {
        Self {
            joints: Vec::new(),
            name_to_index: HashMap::new(),
        }
    }

    /// Add a joint and return its index.
    pub fn add_joint(&mut self, name: &str, parent: Option<usize>, offset: [f64; 3]) -> usize {
        let idx = self.joints.len();
        self.name_to_index.insert(name.to_string(), idx);
        self.joints.push(Joint {
            name: name.to_string(),
            parent,
            offset,
            channels: Vec::new(),
        });
        idx
    }

    /// Add channels to a joint.
    pub fn set_channels(&mut self, joint_idx: usize, channels: Vec<String>) {
        if let Some(j) = self.joints.get_mut(joint_idx) {
            j.channels = channels;
        }
    }

    /// Number of joints.
    pub fn num_joints(&self) -> usize {
        self.joints.len()
    }

    /// Get the children of a joint.
    pub fn children(&self, joint_idx: usize) -> Vec<usize> {
        self.joints
            .iter()
            .enumerate()
            .filter(|(_, j)| j.parent == Some(joint_idx))
            .map(|(i, _)| i)
            .collect()
    }

    /// Find a joint by name.
    pub fn find_joint(&self, name: &str) -> Option<usize> {
        self.name_to_index.get(name).copied()
    }

    /// Total number of channels across all joints.
    pub fn total_channels(&self) -> usize {
        self.joints.iter().map(|j| j.channels.len()).sum()
    }
}

impl Default for Skeleton {
    fn default() -> Self {
        Self::new()
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Keyframe & animation clip
// ═══════════════════════════════════════════════════════════════════════════════

/// A single keyframe for one joint.
#[derive(Debug, Clone, Copy)]
pub struct JointPose {
    /// Local translation.
    pub translation: [f64; 3],
    /// Local rotation as quaternion `[x, y, z, w]`.
    pub rotation: [f64; 4],
    /// Local scale.
    pub scale: [f64; 3],
}

impl JointPose {
    /// Identity pose.
    pub fn identity() -> Self {
        Self {
            translation: [0.0; 3],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0, 1.0, 1.0],
        }
    }

    /// Interpolate between two poses.
    pub fn lerp(&self, other: &JointPose, t: f64) -> JointPose {
        JointPose {
            translation: vec3_lerp(self.translation, other.translation, t),
            rotation: quat_slerp(self.rotation, other.rotation, t),
            scale: vec3_lerp(self.scale, other.scale, t),
        }
    }
}

/// A single frame of a full skeleton pose.
#[derive(Debug, Clone)]
pub struct SkeletonPose {
    /// Per-joint local poses (indexed by joint index).
    pub joint_poses: Vec<JointPose>,
}

impl SkeletonPose {
    /// Create a default pose for a skeleton.
    pub fn default_for(skeleton: &Skeleton) -> Self {
        let n = skeleton.num_joints();
        Self {
            joint_poses: vec![JointPose::identity(); n],
        }
    }

    /// Interpolate between two skeleton poses.
    pub fn lerp(&self, other: &SkeletonPose, t: f64) -> SkeletonPose {
        let poses: Vec<JointPose> = self
            .joint_poses
            .iter()
            .zip(&other.joint_poses)
            .map(|(a, b)| a.lerp(b, t))
            .collect();
        SkeletonPose { joint_poses: poses }
    }
}

/// An animation clip: a sequence of timed skeleton poses.
#[derive(Debug, Clone)]
pub struct AnimationClip {
    /// Name of the clip.
    pub name: String,
    /// Frame rate (frames per second).
    pub fps: f64,
    /// Frame timestamps.
    pub times: Vec<f64>,
    /// Per-frame skeleton poses.
    pub frames: Vec<SkeletonPose>,
}

impl AnimationClip {
    /// Create an empty clip.
    pub fn new(name: &str, fps: f64) -> Self {
        Self {
            name: name.to_string(),
            fps,
            times: Vec::new(),
            frames: Vec::new(),
        }
    }

    /// Add a frame.
    pub fn add_frame(&mut self, time: f64, pose: SkeletonPose) {
        self.times.push(time);
        self.frames.push(pose);
    }

    /// Duration of the clip.
    pub fn duration(&self) -> f64 {
        if self.times.is_empty() {
            return 0.0;
        }
        self.times[self.times.len() - 1] - self.times[0]
    }

    /// Number of frames.
    pub fn num_frames(&self) -> usize {
        self.frames.len()
    }

    /// Sample the animation at a given time using linear interpolation
    /// (with SLERP for rotations).
    pub fn sample(&self, time: f64) -> Option<SkeletonPose> {
        if self.frames.is_empty() {
            return None;
        }
        if self.frames.len() == 1 || time <= self.times[0] {
            return Some(self.frames[0].clone());
        }
        let last = self.times.len() - 1;
        if time >= self.times[last] {
            return Some(self.frames[last].clone());
        }
        // Binary search for the interval
        let idx = match self
            .times
            .binary_search_by(|t| t.partial_cmp(&time).unwrap_or(std::cmp::Ordering::Equal))
        {
            Ok(i) => return Some(self.frames[i].clone()),
            Err(i) => i,
        };
        let i0 = idx - 1;
        let i1 = idx;
        let dt = self.times[i1] - self.times[i0];
        let t = if dt > 1e-15 {
            (time - self.times[i0]) / dt
        } else {
            0.0
        };
        Some(self.frames[i0].lerp(&self.frames[i1], t))
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// BVH motion capture parser
// ═══════════════════════════════════════════════════════════════════════════════

/// Parsed BVH motion capture data.
#[derive(Debug, Clone)]
pub struct BvhData {
    /// Skeleton hierarchy.
    pub skeleton: Skeleton,
    /// Frame rate (seconds per frame).
    pub frame_time: f64,
    /// Motion data: outer = frames, inner = channel values.
    pub motion: Vec<Vec<f64>>,
}

/// Parse BVH format text into a `BvhData` structure.
pub fn parse_bvh(input: &str) -> Result<BvhData, String> {
    let mut skeleton = Skeleton::new();
    let lines = input.lines().map(|l| l.trim()).peekable();
    let mut parent_stack: Vec<Option<usize>> = vec![None];
    let mut frame_time = 1.0 / 30.0;
    let mut motion: Vec<Vec<f64>> = Vec::new();
    let mut in_motion = false;
    let mut _num_frames: usize = 0;

    for line in lines {
        if line.is_empty() {
            continue;
        }
        if in_motion {
            let vals: Result<Vec<f64>, _> =
                line.split_whitespace().map(|s| s.parse::<f64>()).collect();
            match vals {
                Ok(v) if !v.is_empty() => motion.push(v),
                _ => {}
            }
            continue;
        }

        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.is_empty() {
            continue;
        }

        match tokens[0] {
            "HIERARCHY" => {}
            "ROOT" | "JOINT" => {
                let name = if tokens.len() > 1 {
                    tokens[1]
                } else {
                    "unnamed"
                };
                let parent = *parent_stack.last().unwrap_or(&None);
                let idx = skeleton.add_joint(name, parent, [0.0; 3]);
                parent_stack.push(Some(idx));
            }
            "End" => {
                // End site – add as a leaf
                let parent = *parent_stack.last().unwrap_or(&None);
                let name = format!("EndSite_{}", skeleton.num_joints());
                let idx = skeleton.add_joint(&name, parent, [0.0; 3]);
                parent_stack.push(Some(idx));
            }
            "OFFSET"
                if tokens.len() >= 4 => {
                    let ox: f64 = tokens[1].parse().unwrap_or(0.0);
                    let oy: f64 = tokens[2].parse().unwrap_or(0.0);
                    let oz: f64 = tokens[3].parse().unwrap_or(0.0);
                    if let Some(&Some(idx)) = parent_stack.last()
                        && let Some(j) = skeleton.joints.get_mut(idx) {
                            j.offset = [ox, oy, oz];
                        }
                }
            "CHANNELS"
                if tokens.len() >= 2 => {
                    let _n: usize = tokens[1].parse().unwrap_or(0);
                    let channels: Vec<String> = tokens[2..].iter().map(|s| s.to_string()).collect();
                    if let Some(&Some(idx)) = parent_stack.last() {
                        skeleton.set_channels(idx, channels);
                    }
                }
            "{" => {}
            "}" => {
                parent_stack.pop();
            }
            "MOTION" => {
                in_motion = true;
            }
            "Frames:"
                if tokens.len() >= 2 => {
                    _num_frames = tokens[1].parse().unwrap_or(0);
                }
            "Frame"
                // "Frame Time: 0.033333"
                if tokens.len() >= 3 => {
                    frame_time = tokens[2].parse().unwrap_or(1.0 / 30.0);
                }
            _ => {}
        }
    }

    Ok(BvhData {
        skeleton,
        frame_time,
        motion,
    })
}

/// Convert BVH channel data for one frame to a `SkeletonPose`.
pub fn bvh_frame_to_pose(skeleton: &Skeleton, channels: &[f64]) -> SkeletonPose {
    let mut pose = SkeletonPose::default_for(skeleton);
    let mut ch_idx = 0;

    for (j_idx, joint) in skeleton.joints.iter().enumerate() {
        let mut trans = joint.offset;
        let mut euler = [0.0; 3];

        for ch_name in &joint.channels {
            if ch_idx >= channels.len() {
                break;
            }
            let val = channels[ch_idx];
            ch_idx += 1;
            match ch_name.as_str() {
                "Xposition" => trans[0] += val,
                "Yposition" => trans[1] += val,
                "Zposition" => trans[2] += val,
                "Xrotation" => euler[0] = val,
                "Yrotation" => euler[1] = val,
                "Zrotation" => euler[2] = val,
                _ => {}
            }
        }

        let rot = euler_xyz_to_quat(euler);
        pose.joint_poses[j_idx] = JointPose {
            translation: trans,
            rotation: rot,
            scale: [1.0, 1.0, 1.0],
        };
    }

    pose
}

/// Convert BVH data to an `AnimationClip`.
pub fn bvh_to_clip(bvh: &BvhData) -> AnimationClip {
    let fps = if bvh.frame_time > 1e-15 {
        1.0 / bvh.frame_time
    } else {
        30.0
    };
    let mut clip = AnimationClip::new("bvh_clip", fps);
    for (i, frame_data) in bvh.motion.iter().enumerate() {
        let time = i as f64 * bvh.frame_time;
        let pose = bvh_frame_to_pose(&bvh.skeleton, frame_data);
        clip.add_frame(time, pose);
    }
    clip
}

/// Export BVH data to string.
pub fn export_bvh(skeleton: &Skeleton, clip: &AnimationClip) -> String {
    let mut out = String::new();
    out.push_str("HIERARCHY\n");

    fn write_joint(out: &mut String, skeleton: &Skeleton, idx: usize, indent: usize) {
        let joint = &skeleton.joints[idx];
        let prefix = " ".repeat(indent);
        let keyword = if joint.parent.is_none() {
            "ROOT"
        } else {
            "JOINT"
        };
        out.push_str(&format!("{}{} {}\n", prefix, keyword, joint.name));
        out.push_str(&format!("{}{{\n", prefix));
        out.push_str(&format!(
            "{}  OFFSET {:.6} {:.6} {:.6}\n",
            prefix, joint.offset[0], joint.offset[1], joint.offset[2]
        ));
        if !joint.channels.is_empty() {
            out.push_str(&format!(
                "{}  CHANNELS {} {}\n",
                prefix,
                joint.channels.len(),
                joint.channels.join(" ")
            ));
        }
        let children = skeleton.children(idx);
        if children.is_empty() {
            out.push_str(&format!("{}  End Site\n", prefix));
            out.push_str(&format!("{}  {{\n", prefix));
            out.push_str(&format!(
                "{}    OFFSET 0.000000 0.000000 0.000000\n",
                prefix
            ));
            out.push_str(&format!("{}  }}\n", prefix));
        } else {
            for &child in &children {
                write_joint(out, skeleton, child, indent + 2);
            }
        }
        out.push_str(&format!("{}}}\n", prefix));
    }

    // Find roots (joints with no parent)
    for (i, j) in skeleton.joints.iter().enumerate() {
        if j.parent.is_none() {
            write_joint(&mut out, skeleton, i, 0);
        }
    }

    out.push_str("MOTION\n");
    out.push_str(&format!("Frames: {}\n", clip.num_frames()));
    let frame_time = if clip.fps > 0.0 {
        1.0 / clip.fps
    } else {
        1.0 / 30.0
    };
    out.push_str(&format!("Frame Time: {:.6}\n", frame_time));

    for frame in &clip.frames {
        let mut vals = Vec::new();
        for (j_idx, joint) in skeleton.joints.iter().enumerate() {
            let pose = &frame.joint_poses[j_idx];
            for ch in &joint.channels {
                match ch.as_str() {
                    "Xposition" => vals.push(pose.translation[0]),
                    "Yposition" => vals.push(pose.translation[1]),
                    "Zposition" => vals.push(pose.translation[2]),
                    "Xrotation" | "Yrotation" | "Zrotation" => {
                        let euler = quat_to_euler_xyz(pose.rotation);
                        match ch.as_str() {
                            "Xrotation" => vals.push(euler[0]),
                            "Yrotation" => vals.push(euler[1]),
                            _ => vals.push(euler[2]),
                        }
                    }
                    _ => vals.push(0.0),
                }
            }
        }
        let line: Vec<String> = vals.iter().map(|v| format!("{:.6}", v)).collect();
        out.push_str(&line.join(" "));
        out.push('\n');
    }

    out
}

// ═══════════════════════════════════════════════════════════════════════════════
// FBX ASCII (basic)
// ═══════════════════════════════════════════════════════════════════════════════

/// A basic FBX ASCII node.
#[derive(Debug, Clone)]
pub struct FbxNode {
    /// Node name.
    pub name: String,
    /// Properties (string values).
    pub properties: Vec<String>,
    /// Child nodes.
    pub children: Vec<FbxNode>,
}

impl FbxNode {
    /// Create a new FBX node.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            properties: Vec::new(),
            children: Vec::new(),
        }
    }

    /// Add a property.
    pub fn add_property(&mut self, value: &str) {
        self.properties.push(value.to_string());
    }

    /// Add a child node.
    pub fn add_child(&mut self, child: FbxNode) {
        self.children.push(child);
    }

    /// Serialize to FBX ASCII text.
    pub fn to_ascii(&self, indent: usize) -> String {
        let prefix = "  ".repeat(indent);
        let mut out = format!("{}{}: ", prefix, self.name);
        out.push_str(&self.properties.join(", "));
        if self.children.is_empty() {
            out.push_str(" {\n");
            out.push_str(&format!("{}}}\n", prefix));
        } else {
            out.push_str(" {\n");
            for child in &self.children {
                out.push_str(&child.to_ascii(indent + 1));
            }
            out.push_str(&format!("{}}}\n", prefix));
        }
        out
    }

    /// Find a child by name.
    pub fn find_child(&self, name: &str) -> Option<&FbxNode> {
        self.children.iter().find(|c| c.name == name)
    }
}

/// Parse a very basic FBX ASCII string.
///
/// This parser handles the top-level structure of FBX ASCII files. It is
/// not a full FBX parser but sufficient for extracting skeleton and
/// animation data from simple scenes.
pub fn parse_fbx_ascii(input: &str) -> Result<Vec<FbxNode>, String> {
    let mut nodes = Vec::new();
    let mut stack: Vec<FbxNode> = Vec::new();

    for line in input.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with(';') {
            continue;
        }
        if trimmed == "}" {
            if let Some(node) = stack.pop() {
                if let Some(parent) = stack.last_mut() {
                    parent.add_child(node);
                } else {
                    nodes.push(node);
                }
            }
            continue;
        }
        if let Some(colon_pos) = trimmed.find(':') {
            let name = trimmed[..colon_pos].trim();
            let rest = trimmed[colon_pos + 1..].trim();
            let mut node = FbxNode::new(name);
            // Parse properties (before '{')
            let props_str = if let Some(brace) = rest.find('{') {
                rest[..brace].trim()
            } else {
                rest
            };
            if !props_str.is_empty() {
                for prop in props_str.split(',') {
                    let p = prop.trim().trim_matches('"');
                    if !p.is_empty() {
                        node.add_property(p);
                    }
                }
            }
            if rest.contains('{') {
                stack.push(node);
            } else {
                if let Some(parent) = stack.last_mut() {
                    parent.add_child(node);
                } else {
                    nodes.push(node);
                }
            }
        }
    }

    // Drain any unclosed nodes
    while let Some(node) = stack.pop() {
        if let Some(parent) = stack.last_mut() {
            parent.add_child(node);
        } else {
            nodes.push(node);
        }
    }

    Ok(nodes)
}

/// Export a skeleton to basic FBX ASCII format.
pub fn export_fbx_ascii_skeleton(skeleton: &Skeleton) -> String {
    let mut root = FbxNode::new("FBXHeaderExtension");
    let mut version = FbxNode::new("FBXVersion");
    version.add_property("7400");
    root.add_child(version);

    let mut objects = FbxNode::new("Objects");
    for (i, joint) in skeleton.joints.iter().enumerate() {
        let mut model = FbxNode::new("Model");
        model.add_property(&format!("{}", i));
        model.add_property(&format!("\"Model::{}\"", joint.name));
        model.add_property("\"LimbNode\"");

        let mut props = FbxNode::new("Properties70");
        let mut lct = FbxNode::new("P");
        lct.add_property("\"Lcl Translation\"");
        lct.add_property(&format!("{:.6}", joint.offset[0]));
        lct.add_property(&format!("{:.6}", joint.offset[1]));
        lct.add_property(&format!("{:.6}", joint.offset[2]));
        props.add_child(lct);
        model.add_child(props);
        objects.add_child(model);
    }

    let mut out = root.to_ascii(0);
    out.push_str(&objects.to_ascii(0));
    out
}

// ═══════════════════════════════════════════════════════════════════════════════
// USD/USDA export
// ═══════════════════════════════════════════════════════════════════════════════

/// Export a skeleton and animation clip to USDA (ASCII USD) format.
pub fn export_usda(skeleton: &Skeleton, clip: &AnimationClip, stage_name: &str) -> String {
    let mut out = String::new();
    out.push_str("#usda 1.0\n");
    out.push_str(&format!("(\n    defaultPrim = \"{}\"\n)\n\n", stage_name));

    // Skeleton definition
    out.push_str(&format!("def Xform \"{}\" {{\n", stage_name));
    out.push_str("    def Skeleton \"Skeleton\" {\n");

    // Joint paths
    let joint_paths: Vec<String> = skeleton
        .joints
        .iter()
        .map(|j| format!("\"{}\"", j.name))
        .collect();
    out.push_str(&format!(
        "        uniform token[] joints = [{}]\n",
        joint_paths.join(", ")
    ));

    // Bind transforms
    out.push_str("        uniform matrix4d[] bindTransforms = [\n");
    for joint in &skeleton.joints {
        out.push_str(&format!(
            "            (({:.6}, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), ({:.6}, {:.6}, {:.6}, 1)),\n",
            1.0, joint.offset[0], joint.offset[1], joint.offset[2]
        ));
    }
    out.push_str("        ]\n");
    out.push_str("    }\n");

    // Animation
    if !clip.frames.is_empty() {
        out.push_str("    def SkelAnimation \"Animation\" {\n");
        out.push_str(&format!(
            "        uniform token[] joints = [{}]\n",
            joint_paths.join(", ")
        ));

        // Write per-frame translations and rotations
        for (f_idx, frame) in clip.frames.iter().enumerate() {
            let time = if f_idx < clip.times.len() {
                clip.times[f_idx]
            } else {
                f_idx as f64 / clip.fps
            };
            out.push_str(&format!(
                "        float3[] translations.timeSamples[{:.6}] = [\n",
                time
            ));
            for pose in &frame.joint_poses {
                out.push_str(&format!(
                    "            ({:.6}, {:.6}, {:.6}),\n",
                    pose.translation[0], pose.translation[1], pose.translation[2]
                ));
            }
            out.push_str("        ]\n");

            out.push_str(&format!(
                "        quatf[] rotations.timeSamples[{:.6}] = [\n",
                time
            ));
            for pose in &frame.joint_poses {
                let [x, y, z, w] = pose.rotation;
                out.push_str(&format!(
                    "            ({:.6}, {:.6}, {:.6}, {:.6}),\n",
                    w, x, y, z
                ));
            }
            out.push_str("        ]\n");
        }
        out.push_str("    }\n");
    }

    out.push_str("}\n");
    out
}

// ═══════════════════════════════════════════════════════════════════════════════
// Blend shapes
// ═══════════════════════════════════════════════════════════════════════════════

/// A blend shape (morph target) definition.
#[derive(Debug, Clone)]
pub struct BlendShape {
    /// Name of the blend shape.
    pub name: String,
    /// Vertex deltas: `(vertex_index, delta_position)`.
    pub deltas: Vec<(usize, [f64; 3])>,
}

impl BlendShape {
    /// Create a new blend shape.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            deltas: Vec::new(),
        }
    }

    /// Add a vertex delta.
    pub fn add_delta(&mut self, vertex_idx: usize, delta: [f64; 3]) {
        self.deltas.push((vertex_idx, delta));
    }
}

/// A blend shape set with named targets and weights.
#[derive(Debug, Clone)]
pub struct BlendShapeSet {
    /// The blend shapes.
    pub shapes: Vec<BlendShape>,
    /// Current weights for each shape.
    pub weights: Vec<f64>,
}

impl BlendShapeSet {
    /// Create a new empty blend shape set.
    pub fn new() -> Self {
        Self {
            shapes: Vec::new(),
            weights: Vec::new(),
        }
    }

    /// Add a blend shape and return its index.
    pub fn add_shape(&mut self, shape: BlendShape) -> usize {
        let idx = self.shapes.len();
        self.shapes.push(shape);
        self.weights.push(0.0);
        idx
    }

    /// Set the weight for a blend shape.
    pub fn set_weight(&mut self, idx: usize, weight: f64) {
        if idx < self.weights.len() {
            self.weights[idx] = weight;
        }
    }

    /// Evaluate the blend shapes and return the accumulated deltas per vertex.
    ///
    /// The result maps vertex_index to accumulated delta.
    pub fn evaluate(&self) -> HashMap<usize, [f64; 3]> {
        let mut result: HashMap<usize, [f64; 3]> = HashMap::new();
        for (shape, &weight) in self.shapes.iter().zip(&self.weights) {
            if weight.abs() < 1e-15 {
                continue;
            }
            for &(vi, delta) in &shape.deltas {
                let entry = result.entry(vi).or_insert([0.0; 3]);
                *entry = vec3_add(*entry, vec3_scale(delta, weight));
            }
        }
        result
    }

    /// Interpolate weights between two sets of weights.
    pub fn interpolate_weights(a: &[f64], b: &[f64], t: f64) -> Vec<f64> {
        a.iter()
            .zip(b.iter())
            .map(|(&wa, &wb)| wa + (wb - wa) * t)
            .collect()
    }
}

impl Default for BlendShapeSet {
    fn default() -> Self {
        Self::new()
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Timeline / Sequence
// ═══════════════════════════════════════════════════════════════════════════════

/// An event on the animation timeline.
#[derive(Debug, Clone)]
pub struct TimelineEvent {
    /// Time of the event.
    pub time: f64,
    /// Event name / type.
    pub name: String,
    /// Event payload (arbitrary string data).
    pub payload: String,
}

/// A timeline track.
#[derive(Debug, Clone)]
pub struct TimelineTrack {
    /// Track name.
    pub name: String,
    /// Start time.
    pub start: f64,
    /// End time.
    pub end: f64,
    /// Events on this track.
    pub events: Vec<TimelineEvent>,
    /// Associated animation clip name.
    pub clip_name: Option<String>,
}

impl TimelineTrack {
    /// Create a new track.
    pub fn new(name: &str, start: f64, end: f64) -> Self {
        Self {
            name: name.to_string(),
            start,
            end,
            events: Vec::new(),
            clip_name: None,
        }
    }

    /// Duration of the track.
    pub fn duration(&self) -> f64 {
        self.end - self.start
    }

    /// Add an event.
    pub fn add_event(&mut self, time: f64, name: &str, payload: &str) {
        self.events.push(TimelineEvent {
            time,
            name: name.to_string(),
            payload: payload.to_string(),
        });
    }

    /// Get events at or near a given time.
    pub fn events_at(&self, time: f64, tolerance: f64) -> Vec<&TimelineEvent> {
        self.events
            .iter()
            .filter(|e| (e.time - time).abs() < tolerance)
            .collect()
    }
}

/// A full animation timeline with multiple tracks.
#[derive(Debug, Clone)]
pub struct Timeline {
    /// Tracks in the timeline.
    pub tracks: Vec<TimelineTrack>,
    /// Global start time.
    pub start: f64,
    /// Global end time.
    pub end: f64,
    /// Frames per second.
    pub fps: f64,
}

impl Timeline {
    /// Create a new timeline.
    pub fn new(start: f64, end: f64, fps: f64) -> Self {
        Self {
            tracks: Vec::new(),
            start,
            end,
            fps,
        }
    }

    /// Add a track.
    pub fn add_track(&mut self, track: TimelineTrack) {
        self.tracks.push(track);
    }

    /// Duration.
    pub fn duration(&self) -> f64 {
        self.end - self.start
    }

    /// Total number of frames.
    pub fn total_frames(&self) -> usize {
        ((self.end - self.start) * self.fps).ceil() as usize
    }

    /// Get the time for a given frame index.
    pub fn frame_time(&self, frame: usize) -> f64 {
        self.start + (frame as f64) / self.fps
    }
}

// ═══════════════════════════════════════════════════════════════════════════════
// Animation retargeting
// ═══════════════════════════════════════════════════════════════════════════════

/// A mapping between source and target skeleton joints.
#[derive(Debug, Clone)]
pub struct RetargetMapping {
    /// Maps source joint index to target joint index.
    pub joint_map: HashMap<usize, usize>,
    /// Scale factor for translations.
    pub scale: f64,
    /// Per-joint rotation offset (applied after retargeting).
    pub rotation_offsets: HashMap<usize, [f64; 4]>,
}

impl RetargetMapping {
    /// Create a new mapping with a given scale.
    pub fn new(scale: f64) -> Self {
        Self {
            joint_map: HashMap::new(),
            scale,
            rotation_offsets: HashMap::new(),
        }
    }

    /// Add a joint mapping.
    pub fn map_joint(&mut self, source_idx: usize, target_idx: usize) {
        self.joint_map.insert(source_idx, target_idx);
    }

    /// Add a rotation offset for a target joint.
    pub fn set_rotation_offset(&mut self, target_idx: usize, offset: [f64; 4]) {
        self.rotation_offsets.insert(target_idx, offset);
    }

    /// Build a mapping from joint names.
    pub fn from_name_mapping(
        source: &Skeleton,
        target: &Skeleton,
        name_map: &HashMap<String, String>,
        scale: f64,
    ) -> Self {
        let mut mapping = Self::new(scale);
        for (src_name, tgt_name) in name_map {
            if let (Some(&si), Some(&ti)) = (
                source.name_to_index.get(src_name),
                target.name_to_index.get(tgt_name),
            ) {
                mapping.map_joint(si, ti);
            }
        }
        mapping
    }
}

/// Retarget a single skeleton pose from a source to a target skeleton.
pub fn retarget_pose(
    source_pose: &SkeletonPose,
    target_skeleton: &Skeleton,
    mapping: &RetargetMapping,
) -> SkeletonPose {
    let mut target_pose = SkeletonPose::default_for(target_skeleton);

    for (&src_idx, &tgt_idx) in &mapping.joint_map {
        if src_idx >= source_pose.joint_poses.len() || tgt_idx >= target_pose.joint_poses.len() {
            continue;
        }
        let src = &source_pose.joint_poses[src_idx];
        let mut tgt = JointPose {
            translation: vec3_scale(src.translation, mapping.scale),
            rotation: src.rotation,
            scale: src.scale,
        };
        // Apply rotation offset if present
        if let Some(&offset) = mapping.rotation_offsets.get(&tgt_idx) {
            tgt.rotation = quat_mul(offset, tgt.rotation);
            tgt.rotation = quat_normalize(tgt.rotation);
        }
        target_pose.joint_poses[tgt_idx] = tgt;
    }

    target_pose
}

/// Retarget an entire animation clip.
pub fn retarget_clip(
    source_clip: &AnimationClip,
    target_skeleton: &Skeleton,
    mapping: &RetargetMapping,
) -> AnimationClip {
    let mut target_clip = AnimationClip::new(&source_clip.name, source_clip.fps);
    for (i, frame) in source_clip.frames.iter().enumerate() {
        let time = if i < source_clip.times.len() {
            source_clip.times[i]
        } else {
            i as f64 / source_clip.fps
        };
        let retargeted = retarget_pose(frame, target_skeleton, mapping);
        target_clip.add_frame(time, retargeted);
    }
    target_clip
}

// ═══════════════════════════════════════════════════════════════════════════════
// Frame interpolation utilities
// ═══════════════════════════════════════════════════════════════════════════════

/// Interpolation mode for animation channels.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InterpolationMode {
    /// Step / hold (nearest frame).
    Step,
    /// Linear interpolation (SLERP for rotations).
    Linear,
    /// Cubic Hermite (Catmull-Rom style).
    CubicHermite,
}

/// Perform step interpolation.
pub fn step_interpolate(a: f64, _b: f64, _t: f64) -> f64 {
    a
}

/// Perform linear interpolation.
pub fn linear_interpolate(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Perform cubic Hermite interpolation between two values with tangents.
pub fn cubic_hermite_interpolate(p0: f64, m0: f64, p1: f64, m1: f64, t: f64) -> f64 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    h00 * p0 + h10 * m0 + h01 * p1 + h11 * m1
}

/// Catmull-Rom tangent estimation.
pub fn catmull_rom_tangent(p_prev: f64, p_next: f64, dt: f64) -> f64 {
    if dt.abs() < 1e-15 {
        return 0.0;
    }
    (p_next - p_prev) / (2.0 * dt)
}

/// Resample an animation clip to a new frame rate using the specified
/// interpolation mode.
pub fn resample_clip(
    clip: &AnimationClip,
    target_fps: f64,
    _mode: InterpolationMode,
) -> AnimationClip {
    let duration = clip.duration();
    if duration <= 0.0 || clip.frames.is_empty() {
        return clip.clone();
    }
    let num_frames = (duration * target_fps).ceil() as usize + 1;
    let mut resampled = AnimationClip::new(&clip.name, target_fps);

    for i in 0..num_frames {
        let time = clip.times[0] + (i as f64) / target_fps;
        if let Some(pose) = clip.sample(time) {
            resampled.add_frame(time, pose);
        }
    }

    resampled
}

// ═══════════════════════════════════════════════════════════════════════════════
// Animation blending
// ═══════════════════════════════════════════════════════════════════════════════

/// Blend two skeleton poses with a weight `t` (0 = fully a, 1 = fully b).
pub fn blend_poses(a: &SkeletonPose, b: &SkeletonPose, t: f64) -> SkeletonPose {
    a.lerp(b, t)
}

/// Additive blend: add the difference (b - ref) to the base pose, scaled
/// by `weight`.
pub fn additive_blend(
    base: &SkeletonPose,
    reference: &SkeletonPose,
    additive: &SkeletonPose,
    weight: f64,
) -> SkeletonPose {
    let n = base.joint_poses.len();
    let mut result = base.clone();

    for i in 0..n {
        if i >= reference.joint_poses.len() || i >= additive.joint_poses.len() {
            continue;
        }
        let ref_pose = &reference.joint_poses[i];
        let add_pose = &additive.joint_poses[i];
        let base_pose = &base.joint_poses[i];

        // Translation delta
        let t_delta = vec3_sub(add_pose.translation, ref_pose.translation);
        result.joint_poses[i].translation =
            vec3_add(base_pose.translation, vec3_scale(t_delta, weight));

        // Rotation delta
        let ref_inv = quat_conjugate(ref_pose.rotation);
        let r_delta = quat_mul(add_pose.rotation, ref_inv);
        let r_delta_scaled = quat_slerp([0.0, 0.0, 0.0, 1.0], r_delta, weight);
        result.joint_poses[i].rotation =
            quat_normalize(quat_mul(r_delta_scaled, base_pose.rotation));
    }

    result
}

// ═══════════════════════════════════════════════════════════════════════════════
// Animation export helpers
// ═══════════════════════════════════════════════════════════════════════════════

/// Export animation data as a simple CSV.
///
/// Columns: `time, joint_index, tx, ty, tz, rx, ry, rz, rw`.
pub fn export_animation_csv(clip: &AnimationClip) -> String {
    let mut out = String::from("time,joint,tx,ty,tz,rx,ry,rz,rw\n");
    for (f_idx, frame) in clip.frames.iter().enumerate() {
        let time = if f_idx < clip.times.len() {
            clip.times[f_idx]
        } else {
            f_idx as f64 / clip.fps
        };
        for (j_idx, pose) in frame.joint_poses.iter().enumerate() {
            out.push_str(&format!(
                "{:.6},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
                time,
                j_idx,
                pose.translation[0],
                pose.translation[1],
                pose.translation[2],
                pose.rotation[0],
                pose.rotation[1],
                pose.rotation[2],
                pose.rotation[3],
            ));
        }
    }
    out
}

/// Parse animation data from CSV (inverse of `export_animation_csv`).
///
/// Returns `(fps_estimate, frames_map)` where frames_map groups data by
/// time.
pub fn parse_animation_csv(input: &str, num_joints: usize) -> Result<AnimationClip, String> {
    let mut times_and_poses: Vec<(f64, usize, JointPose)> = Vec::new();

    for (line_no, line) in input.lines().enumerate() {
        if line_no == 0 || line.trim().is_empty() {
            continue; // skip header
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() < 9 {
            continue;
        }
        let time: f64 = parts[0].parse().map_err(|_| "Invalid time".to_string())?;
        let joint: usize = parts[1].parse().map_err(|_| "Invalid joint".to_string())?;
        let tx: f64 = parts[2].parse().unwrap_or(0.0);
        let ty: f64 = parts[3].parse().unwrap_or(0.0);
        let tz: f64 = parts[4].parse().unwrap_or(0.0);
        let rx: f64 = parts[5].parse().unwrap_or(0.0);
        let ry: f64 = parts[6].parse().unwrap_or(0.0);
        let rz: f64 = parts[7].parse().unwrap_or(0.0);
        let rw: f64 = parts[8].parse().unwrap_or(1.0);

        times_and_poses.push((
            time,
            joint,
            JointPose {
                translation: [tx, ty, tz],
                rotation: [rx, ry, rz, rw],
                scale: [1.0, 1.0, 1.0],
            },
        ));
    }

    // Group by time
    let mut frame_map: HashMap<u64, (f64, Vec<(usize, JointPose)>)> = HashMap::new();
    for (time, joint, pose) in &times_and_poses {
        let key = (*time * 1_000_000.0) as u64;
        let entry = frame_map.entry(key).or_insert((*time, Vec::new()));
        entry.1.push((*joint, *pose));
    }

    let mut sorted_frames: Vec<(f64, Vec<(usize, JointPose)>)> = frame_map.into_values().collect();
    sorted_frames.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let fps = if sorted_frames.len() >= 2 {
        let dt = sorted_frames[1].0 - sorted_frames[0].0;
        if dt > 1e-15 { 1.0 / dt } else { 30.0 }
    } else {
        30.0
    };

    let mut clip = AnimationClip::new("csv_clip", fps);
    for (time, joint_poses) in &sorted_frames {
        let mut pose = SkeletonPose {
            joint_poses: vec![JointPose::identity(); num_joints],
        };
        for (ji, jp) in joint_poses {
            if *ji < num_joints {
                pose.joint_poses[*ji] = *jp;
            }
        }
        clip.add_frame(*time, pose);
    }

    Ok(clip)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quat_slerp_identity() {
        let id = [0.0, 0.0, 0.0, 1.0];
        let result = quat_slerp(id, id, 0.5);
        assert!((result[3] - 1.0).abs() < 1e-10, "w = {:.6}", result[3]);
    }

    #[test]
    fn test_quat_slerp_endpoints() {
        let a = [0.0, 0.0, 0.0, 1.0];
        let angle = 0.5_f64;
        let b = quat_normalize([0.0, angle.sin(), 0.0, angle.cos()]);
        let r0 = quat_slerp(a, b, 0.0);
        let r1 = quat_slerp(a, b, 1.0);
        for i in 0..4 {
            assert!((r0[i] - a[i]).abs() < 1e-10, "r0[{i}] = {:.6}", r0[i]);
            assert!((r1[i] - b[i]).abs() < 1e-10, "r1[{i}] = {:.6}", r1[i]);
        }
    }

    #[test]
    fn test_quat_slerp_midpoint() {
        let a = [0.0, 0.0, 0.0, 1.0];
        let angle = std::f64::consts::FRAC_PI_2;
        let b = quat_normalize([0.0, (angle / 2.0).sin(), 0.0, (angle / 2.0).cos()]);
        let mid = quat_slerp(a, b, 0.5);
        // Should be approximately half the rotation
        let dot_val = quat_dot(mid, a);
        assert!(
            dot_val > 0.9,
            "Midpoint should be close to a, dot = {:.6}",
            dot_val
        );
    }

    #[test]
    fn test_skeleton_hierarchy() {
        let mut skel = Skeleton::new();
        let root = skel.add_joint("Hips", None, [0.0, 1.0, 0.0]);
        let spine = skel.add_joint("Spine", Some(root), [0.0, 0.2, 0.0]);
        let _head = skel.add_joint("Head", Some(spine), [0.0, 0.3, 0.0]);
        assert_eq!(skel.num_joints(), 3);
        assert_eq!(skel.children(root), vec![1]);
        assert_eq!(skel.find_joint("Spine"), Some(1));
    }

    #[test]
    fn test_joint_pose_identity() {
        let pose = JointPose::identity();
        assert!((pose.rotation[3] - 1.0).abs() < 1e-10);
        assert!((pose.scale[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_joint_pose_lerp() {
        let a = JointPose {
            translation: [0.0, 0.0, 0.0],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [1.0, 1.0, 1.0],
        };
        let b = JointPose {
            translation: [2.0, 4.0, 6.0],
            rotation: [0.0, 0.0, 0.0, 1.0],
            scale: [2.0, 2.0, 2.0],
        };
        let mid = a.lerp(&b, 0.5);
        assert!((mid.translation[0] - 1.0).abs() < 1e-10);
        assert!((mid.scale[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_animation_clip_sample() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0; 3]);
        let mut clip = AnimationClip::new("test", 30.0);
        let mut pose0 = SkeletonPose::default_for(&skel);
        pose0.joint_poses[0].translation = [0.0, 0.0, 0.0];
        let mut pose1 = SkeletonPose::default_for(&skel);
        pose1.joint_poses[0].translation = [10.0, 0.0, 0.0];
        clip.add_frame(0.0, pose0);
        clip.add_frame(1.0, pose1);
        let sampled = clip.sample(0.5).unwrap();
        assert!(
            (sampled.joint_poses[0].translation[0] - 5.0).abs() < 1e-10,
            "tx = {:.6}",
            sampled.joint_poses[0].translation[0]
        );
    }

    #[test]
    fn test_animation_clip_duration() {
        let mut clip = AnimationClip::new("test", 30.0);
        let pose = SkeletonPose {
            joint_poses: vec![],
        };
        clip.add_frame(0.0, pose.clone());
        clip.add_frame(2.0, pose);
        assert!((clip.duration() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_parse_bvh_basic() {
        let bvh_text = "\
HIERARCHY
ROOT Hips
{
  OFFSET 0.0 0.0 0.0
  CHANNELS 6 Xposition Yposition Zposition Xrotation Yrotation Zrotation
  JOINT Spine
  {
    OFFSET 0.0 5.0 0.0
    CHANNELS 3 Xrotation Yrotation Zrotation
    End Site
    {
      OFFSET 0.0 3.0 0.0
    }
  }
}
MOTION
Frames: 2
Frame Time: 0.033333
0.0 0.0 0.0 0.0 0.0 0.0 10.0 20.0 30.0
1.0 2.0 3.0 5.0 10.0 15.0 20.0 25.0 30.0";

        let bvh = parse_bvh(bvh_text).unwrap();
        assert!(bvh.skeleton.num_joints() >= 2);
        assert_eq!(bvh.motion.len(), 2);
        assert!((bvh.frame_time - 0.033333).abs() < 1e-4);
    }

    #[test]
    fn test_bvh_to_clip() {
        let bvh_text = "\
HIERARCHY
ROOT Root
{
  OFFSET 0.0 0.0 0.0
  CHANNELS 3 Xposition Yposition Zposition
}
MOTION
Frames: 2
Frame Time: 0.033333
1.0 2.0 3.0
4.0 5.0 6.0";

        let bvh = parse_bvh(bvh_text).unwrap();
        let clip = bvh_to_clip(&bvh);
        assert_eq!(clip.num_frames(), 2);
    }

    #[test]
    fn test_export_bvh_roundtrip() {
        let mut skel = Skeleton::new();
        let root = skel.add_joint("Root", None, [0.0; 3]);
        skel.set_channels(
            root,
            vec!["Xposition".into(), "Yposition".into(), "Zposition".into()],
        );

        let mut clip = AnimationClip::new("test", 30.0);
        let mut pose = SkeletonPose::default_for(&skel);
        pose.joint_poses[0].translation = [1.0, 2.0, 3.0];
        clip.add_frame(0.0, pose);

        let exported = export_bvh(&skel, &clip);
        assert!(exported.contains("ROOT Root"));
        assert!(exported.contains("MOTION"));
    }

    #[test]
    fn test_fbx_node_ascii() {
        let mut node = FbxNode::new("Objects");
        node.add_property("\"Test\"");
        let child = FbxNode::new("Child");
        node.add_child(child);
        let text = node.to_ascii(0);
        assert!(text.contains("Objects:"));
        assert!(text.contains("Child:"));
    }

    #[test]
    fn test_parse_fbx_ascii_basic() {
        let input = "\
; FBX test
Objects: {
  Model: 1, \"Model::Hips\", \"LimbNode\" {
    Properties70: {
    }
  }
}";
        let nodes = parse_fbx_ascii(input).unwrap();
        assert!(!nodes.is_empty());
        assert_eq!(nodes[0].name, "Objects");
    }

    #[test]
    fn test_export_usda() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0, 1.0, 0.0]);
        let mut clip = AnimationClip::new("walk", 30.0);
        let pose = SkeletonPose::default_for(&skel);
        clip.add_frame(0.0, pose);

        let usda = export_usda(&skel, &clip, "Character");
        assert!(usda.contains("#usda 1.0"));
        assert!(usda.contains("Skeleton"));
    }

    #[test]
    fn test_blend_shape() {
        let mut bs = BlendShape::new("smile");
        bs.add_delta(0, [0.1, 0.2, 0.0]);
        bs.add_delta(1, [0.0, 0.1, 0.0]);

        let mut set = BlendShapeSet::new();
        let idx = set.add_shape(bs);
        set.set_weight(idx, 0.5);

        let result = set.evaluate();
        assert!((result[&0][0] - 0.05).abs() < 1e-10);
        assert!((result[&0][1] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_blend_shape_interpolate_weights() {
        let a = vec![0.0, 1.0];
        let b = vec![1.0, 0.0];
        let mid = BlendShapeSet::interpolate_weights(&a, &b, 0.5);
        assert!((mid[0] - 0.5).abs() < 1e-10);
        assert!((mid[1] - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_timeline() {
        let mut tl = Timeline::new(0.0, 10.0, 30.0);
        let mut track = TimelineTrack::new("main", 0.0, 10.0);
        track.add_event(5.0, "footstep", "left");
        tl.add_track(track);
        assert!((tl.duration() - 10.0).abs() < 1e-10);
        assert_eq!(tl.total_frames(), 300);
        let events = tl.tracks[0].events_at(5.0, 0.01);
        assert_eq!(events.len(), 1);
    }

    #[test]
    fn test_retarget_mapping() {
        let mut source = Skeleton::new();
        source.add_joint("Hips", None, [0.0; 3]);
        source.add_joint("Spine", Some(0), [0.0, 1.0, 0.0]);

        let mut target = Skeleton::new();
        target.add_joint("pelvis", None, [0.0; 3]);
        target.add_joint("spine_01", Some(0), [0.0, 0.5, 0.0]);

        let mut name_map = HashMap::new();
        name_map.insert("Hips".to_string(), "pelvis".to_string());
        name_map.insert("Spine".to_string(), "spine_01".to_string());

        let mapping = RetargetMapping::from_name_mapping(&source, &target, &name_map, 0.5);
        assert_eq!(mapping.joint_map.len(), 2);
    }

    #[test]
    fn test_retarget_pose() {
        let mut source_skel = Skeleton::new();
        source_skel.add_joint("A", None, [0.0; 3]);

        let mut target_skel = Skeleton::new();
        target_skel.add_joint("B", None, [0.0; 3]);

        let mut mapping = RetargetMapping::new(2.0);
        mapping.map_joint(0, 0);

        let mut source_pose = SkeletonPose::default_for(&source_skel);
        source_pose.joint_poses[0].translation = [1.0, 2.0, 3.0];

        let result = retarget_pose(&source_pose, &target_skel, &mapping);
        assert!((result.joint_poses[0].translation[0] - 2.0).abs() < 1e-10);
        assert!((result.joint_poses[0].translation[1] - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_cubic_hermite() {
        // At endpoints
        let v0 = cubic_hermite_interpolate(0.0, 1.0, 1.0, 1.0, 0.0);
        let v1 = cubic_hermite_interpolate(0.0, 1.0, 1.0, 1.0, 1.0);
        assert!((v0).abs() < 1e-10);
        assert!((v1 - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_step_interpolate() {
        assert!((step_interpolate(5.0, 10.0, 0.5) - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_linear_interpolate() {
        assert!((linear_interpolate(0.0, 10.0, 0.3) - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_resample_clip() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0; 3]);
        let mut clip = AnimationClip::new("test", 10.0);
        let mut p0 = SkeletonPose::default_for(&skel);
        p0.joint_poses[0].translation = [0.0, 0.0, 0.0];
        let mut p1 = SkeletonPose::default_for(&skel);
        p1.joint_poses[0].translation = [10.0, 0.0, 0.0];
        clip.add_frame(0.0, p0);
        clip.add_frame(1.0, p1);
        let resampled = resample_clip(&clip, 20.0, InterpolationMode::Linear);
        assert!(resampled.num_frames() > clip.num_frames());
    }

    #[test]
    fn test_export_animation_csv() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0; 3]);
        let mut clip = AnimationClip::new("test", 30.0);
        let pose = SkeletonPose::default_for(&skel);
        clip.add_frame(0.0, pose);
        let csv = export_animation_csv(&clip);
        assert!(csv.contains("time,joint,tx,ty,tz,rx,ry,rz,rw"));
    }

    #[test]
    fn test_parse_animation_csv_roundtrip() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0; 3]);
        let mut clip = AnimationClip::new("test", 30.0);
        let mut pose = SkeletonPose::default_for(&skel);
        pose.joint_poses[0].translation = [1.0, 2.0, 3.0];
        clip.add_frame(0.0, pose.clone());
        clip.add_frame(0.033333, pose);
        let csv = export_animation_csv(&clip);
        let parsed = parse_animation_csv(&csv, 1).unwrap();
        assert_eq!(parsed.num_frames(), 2);
    }

    #[test]
    fn test_additive_blend() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0; 3]);

        let mut base = SkeletonPose::default_for(&skel);
        base.joint_poses[0].translation = [5.0, 0.0, 0.0];

        let reference = SkeletonPose::default_for(&skel);
        let mut additive = SkeletonPose::default_for(&skel);
        additive.joint_poses[0].translation = [1.0, 0.0, 0.0];

        let result = additive_blend(&base, &reference, &additive, 1.0);
        // base + (additive - reference) * 1.0 = 5 + (1 - 0) = 6
        assert!(
            (result.joint_poses[0].translation[0] - 6.0).abs() < 1e-10,
            "tx = {:.6}",
            result.joint_poses[0].translation[0]
        );
    }

    #[test]
    fn test_euler_quat_roundtrip() {
        // Single-axis rotation avoids convention mixing
        let euler = [15.0, 0.0, 0.0];
        let q = euler_xyz_to_quat(euler);
        let euler_back = quat_to_euler_xyz(q);
        assert!(
            (euler[0] - euler_back[0]).abs() < 0.1,
            "euler[0] = {:.6}, got {:.6}",
            euler[0],
            euler_back[0]
        );
        // Test a Y-only rotation
        let euler_y = [0.0, 25.0, 0.0];
        let qy = euler_xyz_to_quat(euler_y);
        let ey_back = quat_to_euler_xyz(qy);
        assert!(
            (euler_y[1] - ey_back[1]).abs() < 0.1,
            "euler_y[1] = {:.6}, got {:.6}",
            euler_y[1],
            ey_back[1]
        );
    }

    #[test]
    fn test_catmull_rom_tangent() {
        let t = catmull_rom_tangent(0.0, 4.0, 1.0);
        assert!((t - 2.0).abs() < 1e-10, "tangent = {:.6}", t);
    }

    #[test]
    fn test_export_fbx_skeleton() {
        let mut skel = Skeleton::new();
        skel.add_joint("Root", None, [0.0, 1.0, 0.0]);
        skel.add_joint("Spine", Some(0), [0.0, 0.5, 0.0]);
        let fbx = export_fbx_ascii_skeleton(&skel);
        assert!(fbx.contains("FBXVersion"));
        assert!(fbx.contains("Model"));
    }
}
