// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! BVH motion-capture import and export.
//!
//! # Conventions
//!
//! - **Rotation order.** The rotation channels of a joint are applied in the
//!   order the `CHANNELS` line lists them.  With `Zrotation Xrotation
//!   Yrotation` and the values `(a, b, c)` the local rotation is
//!   `R = Rz(a) · Rx(b) · Ry(c)`.  [`bvh_frame_to_pose`] composes the channels
//!   that way and [`export_bvh`] decomposes a rotation in the same order, so
//!   the two are inverse to each other (up to the usual ambiguity of Euler
//!   angles at gimbal lock).
//! - **Translation.** A position channel adds its value to the `OFFSET` of the
//!   joint, so the pose translation is the offset plus the channel values;
//!   the exporter writes the pose translation minus the offset.
//! - **End sites.** An `End Site` block becomes a leaf joint flagged
//!   [`Joint::is_end_site`], and the exporter writes such a joint back as an
//!   `End Site` block, so `parse → export → parse` keeps the hierarchy.
//! - **Channel values** are consumed in joint-index order, which for a
//!   skeleton produced by [`parse_bvh`] is the order of the file.

use super::{AnimationClip, Joint, JointPose, Skeleton, SkeletonPose, quat_mul, quat_normalize};
use std::fmt::Write as _;

/// Frame time assumed when a file does not contain a `Frame Time:` line.
const DEFAULT_FRAME_TIME: f64 = 1.0 / 30.0;

/// Smallest `cos(beta)` for which the outer Euler angles are still separable;
/// below it the decomposition is at gimbal lock.
const GIMBAL_COS_EPSILON: f64 = 1e-9;

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

/// Coordinate axis of a BVH channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Axis {
    X,
    Y,
    Z,
}

impl Axis {
    const ALL: [Axis; 3] = [Axis::X, Axis::Y, Axis::Z];

    fn index(self) -> usize {
        match self {
            Axis::X => 0,
            Axis::Y => 1,
            Axis::Z => 2,
        }
    }
}

/// A BVH channel: a translation along, or a rotation about, one axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Channel {
    Position(Axis),
    Rotation(Axis),
}

impl Channel {
    /// Parse a channel name such as `Zrotation` (the case is ignored).
    fn parse(name: &str) -> Option<Self> {
        let lower = name.to_ascii_lowercase();
        let axis = match lower.get(..1)? {
            "x" => Axis::X,
            "y" => Axis::Y,
            "z" => Axis::Z,
            _ => return None,
        };
        match lower.get(1..)? {
            "position" => Some(Channel::Position(axis)),
            "rotation" => Some(Channel::Rotation(axis)),
            _ => None,
        }
    }
}

// ── Rotations ──────────────────────────────────────────────────────────────

/// Rotation about a coordinate axis by `degrees`, as a quaternion
/// `[x, y, z, w]`.
fn axis_angle_quat(axis: Axis, degrees: f64) -> [f64; 4] {
    let (sin_half, cos_half) = (degrees.to_radians() * 0.5).sin_cos();
    let mut q = [0.0, 0.0, 0.0, cos_half];
    q[axis.index()] = sin_half;
    q
}

/// Row-major rotation matrix of the quaternion `q = [x, y, z, w]` (which is
/// normalised first).
fn quat_to_matrix(q: [f64; 4]) -> [[f64; 3]; 3] {
    let [x, y, z, w] = quat_normalize(q);
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// Decompose the rotation `q` into three rotations about the distinct axes of
/// `order`, in degrees: the result `[a, b, c]` satisfies
/// `R(q) = R_order[0](a) · R_order[1](b) · R_order[2](c)`.
///
/// The middle angle lies in `[-90, 90]` degrees and the outer ones in
/// `(-180, 180]`.  At gimbal lock (middle angle of ±90 degrees) only the sum
/// or difference of the outer angles is defined; the last one is then zero.
///
/// All six Tait-Bryan orders share one formula.  Relabelling the axes as
/// `order` turns the product into the `X·Y·Z` product, whose matrix elements
/// give the angles directly; an odd permutation of the axes mirrors the
/// coordinate system, which flips the sign of every angle.
fn quat_to_euler(q: [f64; 4], order: [Axis; 3]) -> [f64; 3] {
    let m = quat_to_matrix(q);
    let axes = order.map(Axis::index);
    let sign = if axes[1] == (axes[0] + 1) % 3 {
        1.0
    } else {
        -1.0
    };
    // Element (a, b) of the matrix in the relabelled coordinate system.
    let r = |a: usize, b: usize| m[axes[a]][axes[b]];
    let cos_beta = r(1, 2).hypot(r(2, 2));
    let beta = r(0, 2).atan2(cos_beta);
    let (alpha, gamma) = if cos_beta > GIMBAL_COS_EPSILON {
        ((-r(1, 2)).atan2(r(2, 2)), (-r(0, 1)).atan2(r(0, 0)))
    } else {
        (r(2, 1).atan2(r(1, 1)), 0.0)
    };
    [
        sign * alpha.to_degrees(),
        sign * beta.to_degrees(),
        sign * gamma.to_degrees(),
    ]
}

// ── Parsing ────────────────────────────────────────────────────────────────

/// Meaning of one line of the `MOTION` section.
enum MotionLine {
    /// The `Frames:` line (its number is validated, not used).
    FrameCount,
    /// The `Frame Time:` line.
    FrameTime(f64),
    /// A row of channel values.
    Row(Vec<f64>),
}

/// Interpret one non-blank line of the `MOTION` section.
fn parse_motion_line(
    tokens: &[&str],
    line_no: usize,
    channel_count: usize,
) -> Result<MotionLine, String> {
    let starts_with_letter = tokens
        .first()
        .and_then(|t| t.as_bytes().first())
        .is_some_and(u8::is_ascii_alphabetic);
    if starts_with_letter {
        let header = tokens.join(" ").to_ascii_lowercase();
        if let Some(value) = header.strip_prefix("frames:") {
            let value = value.trim();
            value
                .parse::<usize>()
                .map_err(|_| format!("line {line_no}: invalid frame count `{value}`"))?;
            return Ok(MotionLine::FrameCount);
        }
        if let Some(value) = header.strip_prefix("frame time:") {
            let value = value.trim();
            let seconds = value
                .parse::<f64>()
                .ok()
                .filter(|s| s.is_finite())
                .ok_or_else(|| format!("line {line_no}: invalid frame time `{value}`"))?;
            return Ok(MotionLine::FrameTime(seconds));
        }
    }
    let row = tokens
        .iter()
        .map(|token| {
            token
                .parse::<f64>()
                .ok()
                .filter(|value| value.is_finite())
                .ok_or_else(|| format!("line {line_no}: invalid motion value `{token}`"))
        })
        .collect::<Result<Vec<f64>, String>>()?;
    if row.len() != channel_count {
        return Err(format!(
            "line {line_no}: motion row holds {} values but the hierarchy declares {channel_count} channels",
            row.len()
        ));
    }
    Ok(MotionLine::Row(row))
}

/// Parse an `OFFSET x y z` line.
fn parse_offset(tokens: &[&str], line_no: usize) -> Result<[f64; 3], String> {
    let components = tokens
        .get(1..4)
        .ok_or_else(|| format!("line {line_no}: OFFSET needs three numbers"))?;
    let mut offset = [0.0; 3];
    for (slot, token) in offset.iter_mut().zip(components) {
        *slot = token
            .parse::<f64>()
            .ok()
            .filter(|value| value.is_finite())
            .ok_or_else(|| format!("line {line_no}: invalid OFFSET component `{token}`"))?;
    }
    Ok(offset)
}

/// Parse a `CHANNELS n name...` line into the channel names.
fn parse_channels(tokens: &[&str], line_no: usize) -> Result<Vec<String>, String> {
    let declared = tokens
        .get(1)
        .and_then(|t| t.parse::<usize>().ok())
        .ok_or_else(|| format!("line {line_no}: CHANNELS needs a channel count"))?;
    let names = tokens.get(2..).unwrap_or(&[]);
    if names.len() != declared {
        return Err(format!(
            "line {line_no}: CHANNELS announces {declared} channels but lists {}",
            names.len()
        ));
    }
    if let Some(unknown) = names.iter().find(|name| Channel::parse(name).is_none()) {
        return Err(format!("line {line_no}: unknown channel `{unknown}`"));
    }
    Ok(names.iter().map(|name| (*name).to_string()).collect())
}

/// Parse BVH format text into a `BvhData` structure.
///
/// The `Frames:` and `Frame Time:` lines of the `MOTION` section are read
/// (`Frame Time` defaults to 1/30 s when absent).  The number of frames
/// announced by `Frames:` is not enforced, so a truncated recording still
/// loads with the rows it has.
///
/// # Errors
///
/// Returns a message naming the offending (1-based) line when
///
/// - a `CHANNELS` line lists an unknown channel or announces a different
///   number of channels than it lists,
/// - an `OFFSET` line does not hold three numbers,
/// - `Frames:` or `Frame Time:` has no valid number,
/// - a motion row contains something that is not a number, or does not hold
///   exactly one value per channel of the hierarchy.
pub fn parse_bvh(input: &str) -> Result<BvhData, String> {
    let mut skeleton = Skeleton::new();
    // Joint whose block is open; the bottom entry stands for "no joint yet".
    let mut parent_stack: Vec<Option<usize>> = vec![None];
    let mut frame_time = DEFAULT_FRAME_TIME;
    let mut motion: Vec<Vec<f64>> = Vec::new();
    let mut in_motion = false;
    let mut channel_count = 0usize;

    for (line_index, raw_line) in input.lines().enumerate() {
        let line_no = line_index + 1;
        let tokens: Vec<&str> = raw_line.split_whitespace().collect();
        let Some(&keyword) = tokens.first() else {
            continue;
        };

        if in_motion {
            match parse_motion_line(&tokens, line_no, channel_count)? {
                MotionLine::FrameCount => {}
                MotionLine::FrameTime(seconds) => frame_time = seconds,
                MotionLine::Row(row) => motion.push(row),
            }
            continue;
        }

        match keyword {
            "ROOT" | "JOINT" => {
                let name = tokens.get(1).copied().unwrap_or("unnamed");
                let parent = parent_stack.last().copied().flatten();
                let idx = skeleton.add_joint(name, parent, [0.0; 3]);
                parent_stack.push(Some(idx));
            }
            "End" => {
                let parent = parent_stack.last().copied().flatten();
                let idx = skeleton.add_end_site(parent, [0.0; 3]);
                parent_stack.push(Some(idx));
            }
            "OFFSET" => {
                let offset = parse_offset(&tokens, line_no)?;
                if let Some(&Some(idx)) = parent_stack.last()
                    && let Some(joint) = skeleton.joints.get_mut(idx)
                {
                    joint.offset = offset;
                }
            }
            "CHANNELS" => {
                let channels = parse_channels(&tokens, line_no)?;
                if let Some(&Some(idx)) = parent_stack.last() {
                    skeleton.set_channels(idx, channels);
                }
            }
            "}" => {
                parent_stack.pop();
            }
            "MOTION" => {
                in_motion = true;
                channel_count = skeleton.total_channels();
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
///
/// The values are consumed joint by joint in joint-index order, channel by
/// channel.  A position channel adds its value to the joint offset.  The
/// rotation channels are composed in the order the joint lists them: starting
/// from the identity, each rotation channel post-multiplies the running
/// rotation, so `Zrotation Xrotation Yrotation` gives `Rz · Rx · Ry`.  If
/// `channels` is shorter than the skeleton needs, the remaining joints keep
/// their rest pose.
pub fn bvh_frame_to_pose(skeleton: &Skeleton, channels: &[f64]) -> SkeletonPose {
    let mut pose = SkeletonPose::default_for(skeleton);
    let mut values = channels.iter().copied();

    for (joint, joint_pose) in skeleton.joints.iter().zip(pose.joint_poses.iter_mut()) {
        let mut translation = joint.offset;
        let mut rotation = [0.0, 0.0, 0.0, 1.0];

        for name in &joint.channels {
            let Some(value) = values.next() else {
                break;
            };
            match Channel::parse(name) {
                Some(Channel::Position(axis)) => translation[axis.index()] += value,
                Some(Channel::Rotation(axis)) => {
                    rotation = quat_mul(rotation, axis_angle_quat(axis, value));
                }
                None => {}
            }
        }

        *joint_pose = JointPose {
            translation,
            rotation: quat_normalize(rotation),
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

// ── Export ─────────────────────────────────────────────────────────────────

/// What the exporter needs to know about one joint to write its channels.
struct JointExport<'a> {
    /// Index of the joint in the skeleton (and in every `SkeletonPose`).
    index: usize,
    joint: &'a Joint,
    /// The joint's channels, in file order.
    channels: Vec<Channel>,
    /// The three rotation axes in which the rotation is decomposed: the
    /// joint's rotation channels in file order, completed by the missing axes.
    euler_order: Option<[Axis; 3]>,
}

/// One step of the depth-first walk that writes the hierarchy.
enum Visit {
    Open { index: usize, depth: usize },
    Close { depth: usize },
}

/// Roots (in index order) and the children of every joint (in index order).
fn joint_tree(skeleton: &Skeleton) -> (Vec<usize>, Vec<Vec<usize>>) {
    let mut roots = Vec::new();
    let mut children: Vec<Vec<usize>> = vec![Vec::new(); skeleton.joints.len()];
    for (index, joint) in skeleton.joints.iter().enumerate() {
        match joint.parent {
            None => roots.push(index),
            Some(parent) => {
                if let Some(list) = children.get_mut(parent) {
                    list.push(index);
                }
            }
        }
    }
    (roots, children)
}

/// Parse the channel names of `joint` and derive its Euler order.
fn plan_channels(joint: &Joint) -> Result<(Vec<Channel>, Option<[Axis; 3]>), String> {
    let mut channels = Vec::with_capacity(joint.channels.len());
    for name in &joint.channels {
        let channel = Channel::parse(name)
            .ok_or_else(|| format!("joint `{}` has the unknown channel `{name}`", joint.name))?;
        if channels.contains(&channel) {
            return Err(format!(
                "joint `{}` lists the channel `{name}` twice",
                joint.name
            ));
        }
        channels.push(channel);
    }
    let mut axes: Vec<Axis> = channels
        .iter()
        .filter_map(|channel| match channel {
            Channel::Rotation(axis) => Some(*axis),
            Channel::Position(_) => None,
        })
        .collect();
    if axes.is_empty() {
        return Ok((channels, None));
    }
    for axis in Axis::ALL {
        if !axes.contains(&axis) {
            axes.push(axis);
        }
    }
    let order = <[Axis; 3]>::try_from(axes)
        .map_err(|_| format!("joint `{}` has an invalid rotation order", joint.name))?;
    Ok((channels, Some(order)))
}

/// Write the opening of `joint` (keyword, braces, offset, channels) and return
/// what is needed to write its motion values.  `None` is returned for an
/// `End Site`, which has no channels.
fn open_joint<'a>(
    out: &mut String,
    joint: &'a Joint,
    index: usize,
    depth: usize,
    is_leaf: bool,
) -> Result<Option<JointExport<'a>>, String> {
    let pad = "  ".repeat(depth);
    let [ox, oy, oz] = joint.offset;
    let fmt = |e: std::fmt::Error| format!("formatting the BVH hierarchy failed: {e}");

    if joint.is_end_site {
        if joint.parent.is_none() {
            return Err(format!("End Site `{}` has no parent joint", joint.name));
        }
        if !joint.channels.is_empty() || !is_leaf {
            return Err(format!(
                "End Site `{}` must have neither channels nor children",
                joint.name
            ));
        }
        writeln!(out, "{pad}End Site").map_err(fmt)?;
        writeln!(out, "{pad}{{").map_err(fmt)?;
        writeln!(out, "{pad}  OFFSET {ox:.6} {oy:.6} {oz:.6}").map_err(fmt)?;
        return Ok(None);
    }

    if joint.name.is_empty() || joint.name.contains(char::is_whitespace) {
        return Err(format!(
            "joint {index} has the name `{}`, which is empty or contains whitespace",
            joint.name
        ));
    }
    let (channels, euler_order) = plan_channels(joint)?;
    let keyword = if joint.parent.is_none() {
        "ROOT"
    } else {
        "JOINT"
    };
    writeln!(out, "{pad}{keyword} {}", joint.name).map_err(fmt)?;
    writeln!(out, "{pad}{{").map_err(fmt)?;
    writeln!(out, "{pad}  OFFSET {ox:.6} {oy:.6} {oz:.6}").map_err(fmt)?;
    if !channels.is_empty() {
        writeln!(
            out,
            "{pad}  CHANNELS {} {}",
            channels.len(),
            joint.channels.join(" ")
        )
        .map_err(fmt)?;
    }
    if is_leaf {
        // BVH requires every chain of joints to end in an End Site.
        writeln!(out, "{pad}  End Site").map_err(fmt)?;
        writeln!(out, "{pad}  {{").map_err(fmt)?;
        writeln!(out, "{pad}    OFFSET 0.000000 0.000000 0.000000").map_err(fmt)?;
        writeln!(out, "{pad}  }}").map_err(fmt)?;
    }
    Ok(Some(JointExport {
        index,
        joint,
        channels,
        euler_order,
    }))
}

/// Append the values of the channels of one joint for the pose `pose`.
fn push_pose_values(values: &mut Vec<String>, export: &JointExport<'_>, pose: &JointPose) {
    // Angles of the three axes, by axis; zero for an axis the joint does not
    // rotate about.
    let mut angle_by_axis = [0.0; 3];
    if let Some(order) = export.euler_order {
        let angles = quat_to_euler(pose.rotation, order);
        for (axis, angle) in order.iter().zip(angles) {
            angle_by_axis[axis.index()] = angle;
        }
    }
    for channel in &export.channels {
        let value = match *channel {
            Channel::Position(axis) => {
                pose.translation[axis.index()] - export.joint.offset[axis.index()]
            }
            Channel::Rotation(axis) => angle_by_axis[axis.index()],
        };
        // `+ 0.0` turns a negative zero into a plain zero.
        values.push(format!("{:.6}", value + 0.0));
    }
}

/// Export a clip as BVH text.
///
/// The hierarchy is written depth first (roots and children in joint-index
/// order) and the motion values follow the same order, so the file is
/// self-consistent for any parent/child layout.  Every rotation is decomposed
/// into Euler angles in the order of the joint's rotation channels (all six
/// Tait-Bryan orders are supported), the inverse of what [`parse_bvh`] and
/// [`bvh_frame_to_pose`] apply.  A joint with fewer than three rotation
/// channels is decomposed with the missing axes appended to its order, and
/// only its own channels are written.  Scale is not representable in BVH and
/// is dropped.  A joint that is a leaf but not an end site gets an `End Site`
/// with a zero offset, because BVH requires every chain to end in one.
///
/// # Errors
///
/// Returns a message when the skeleton or the clip cannot be expressed as BVH:
///
/// - a joint is unreachable from a root (a cycle or a missing parent), has an
///   empty or blank-containing name, an unknown channel, a repeated channel,
///   or is an end site with channels, children or no parent,
/// - a frame does not hold exactly one pose per joint.
pub fn export_bvh(skeleton: &Skeleton, clip: &AnimationClip) -> Result<String, String> {
    let fmt = |e: std::fmt::Error| format!("formatting the BVH file failed: {e}");
    let (roots, children) = joint_tree(skeleton);
    let mut out = String::from("HIERARCHY\n");
    let mut exports: Vec<JointExport<'_>> = Vec::with_capacity(skeleton.joints.len());
    let mut visited = 0usize;

    let mut stack: Vec<Visit> = roots
        .iter()
        .rev()
        .map(|&index| Visit::Open { index, depth: 0 })
        .collect();
    while let Some(visit) = stack.pop() {
        match visit {
            Visit::Open { index, depth } => {
                let joint = skeleton
                    .joints
                    .get(index)
                    .ok_or_else(|| format!("joint index {index} is out of range"))?;
                let kids = children.get(index).map(Vec::as_slice).unwrap_or_default();
                visited += 1;
                if let Some(export) = open_joint(&mut out, joint, index, depth, kids.is_empty())? {
                    exports.push(export);
                }
                stack.push(Visit::Close { depth });
                stack.extend(kids.iter().rev().map(|&child| Visit::Open {
                    index: child,
                    depth: depth + 1,
                }));
            }
            Visit::Close { depth } => {
                writeln!(out, "{}}}", "  ".repeat(depth)).map_err(fmt)?;
            }
        }
    }
    if visited != skeleton.joints.len() {
        return Err(format!(
            "{} of {} joints are not reachable from a root joint (cycle or missing parent)",
            skeleton.joints.len() - visited,
            skeleton.joints.len()
        ));
    }

    out.push_str("MOTION\n");
    writeln!(out, "Frames: {}", clip.num_frames()).map_err(fmt)?;
    let frame_time = if clip.fps.is_finite() && clip.fps > 0.0 {
        1.0 / clip.fps
    } else {
        DEFAULT_FRAME_TIME
    };
    writeln!(out, "Frame Time: {frame_time:.6}").map_err(fmt)?;

    for (frame_index, frame) in clip.frames.iter().enumerate() {
        if frame.joint_poses.len() != skeleton.joints.len() {
            return Err(format!(
                "frame {frame_index} holds {} joint poses but the skeleton has {} joints",
                frame.joint_poses.len(),
                skeleton.joints.len()
            ));
        }
        let mut values: Vec<String> = Vec::new();
        for export in &exports {
            let pose = frame.joint_poses.get(export.index).ok_or_else(|| {
                format!("frame {frame_index} has no pose for joint {}", export.index)
            })?;
            push_pose_values(&mut values, export, pose);
        }
        out.push_str(&values.join(" "));
        out.push('\n');
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    type Matrix = [[f64; 3]; 3];

    fn rot(axis: Axis, degrees: f64) -> Matrix {
        let (s, c) = degrees.to_radians().sin_cos();
        match axis {
            Axis::X => [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]],
            Axis::Y => [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]],
            Axis::Z => [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
        }
    }

    fn mat_mul(a: &Matrix, b: &Matrix) -> Matrix {
        let mut out = [[0.0; 3]; 3];
        for (i, row) in out.iter_mut().enumerate() {
            for (j, cell) in row.iter_mut().enumerate() {
                *cell = (0..3).map(|k| a[i][k] * b[k][j]).sum();
            }
        }
        out
    }

    /// `R_order[0](angles[0]) · R_order[1](angles[1]) · R_order[2](angles[2])`
    /// built from explicit 3x3 matrices, independent of the quaternion code.
    fn product(order: [Axis; 3], angles: [f64; 3]) -> Matrix {
        let first = mat_mul(&rot(order[0], angles[0]), &rot(order[1], angles[1]));
        mat_mul(&first, &rot(order[2], angles[2]))
    }

    fn assert_matrix_close(actual: &Matrix, expected: &Matrix, tolerance: f64) {
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (actual[i][j] - expected[i][j]).abs() < tolerance,
                    "element ({i}, {j}): {} vs {}\n{actual:?}\n{expected:?}",
                    actual[i][j],
                    expected[i][j]
                );
            }
        }
    }

    const ALL_ORDERS: [[Axis; 3]; 6] = [
        [Axis::X, Axis::Y, Axis::Z],
        [Axis::X, Axis::Z, Axis::Y],
        [Axis::Y, Axis::X, Axis::Z],
        [Axis::Y, Axis::Z, Axis::X],
        [Axis::Z, Axis::X, Axis::Y],
        [Axis::Z, Axis::Y, Axis::X],
    ];

    fn channel_names(order: [Axis; 3]) -> String {
        order
            .iter()
            .map(|axis| match axis {
                Axis::X => "Xrotation",
                Axis::Y => "Yrotation",
                Axis::Z => "Zrotation",
            })
            .collect::<Vec<_>>()
            .join(" ")
    }

    /// A one-joint BVH file with a single frame.
    fn single_joint_bvh(channels: &str, count: usize, values: &str) -> String {
        format!(
            "HIERARCHY\nROOT Hips\n{{\n  OFFSET 0.0 0.0 0.0\n  CHANNELS {count} {channels}\n  End Site\n  {{\n    OFFSET 0.0 1.0 0.0\n  }}\n}}\nMOTION\nFrames: 1\nFrame Time: 0.008333\n{values}\n"
        )
    }

    fn root_rotation(bvh_text: &str) -> Matrix {
        let bvh = parse_bvh(bvh_text).unwrap();
        let clip = bvh_to_clip(&bvh);
        quat_to_matrix(clip.frames[0].joint_poses[0].rotation)
    }

    // ── parsing ────────────────────────────────────────────────────────────

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
        assert_eq!(bvh.skeleton.num_joints(), 3);
        assert_eq!(bvh.skeleton.joints[1].offset, [0.0, 5.0, 0.0]);
        assert!(bvh.skeleton.joints[2].is_end_site);
        assert_eq!(bvh.skeleton.joints[2].offset, [0.0, 3.0, 0.0]);
        assert_eq!(bvh.skeleton.joints[2].parent, Some(1));
        assert_eq!(bvh.motion.len(), 2);
        assert!((bvh.frame_time - 0.033333).abs() < 1e-9);
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
        assert_eq!(clip.frames[1].joint_poses[0].translation, [4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_frame_time_line_of_the_motion_section_is_read() {
        let text = single_joint_bvh("Zrotation Xrotation Yrotation", 3, "0 0 0");
        let bvh = parse_bvh(&text).unwrap();
        assert_eq!(bvh.frame_time, 0.008333);
        let clip = bvh_to_clip(&bvh);
        assert!((clip.fps - 1.0 / 0.008333).abs() < 1e-9, "fps {}", clip.fps);
    }

    #[test]
    fn test_frame_time_defaults_to_one_thirtieth_and_accepts_loose_layout() {
        let without =
            "HIERARCHY\nROOT R\n{\n OFFSET 0 0 0\n CHANNELS 1 Xposition\n}\nMOTION\n1.5\n";
        let bvh = parse_bvh(without).unwrap();
        assert_eq!(bvh.frame_time, 1.0 / 30.0);
        assert_eq!(bvh.motion, vec![vec![1.5]]);
        let loose = "HIERARCHY\nROOT R\n{\n OFFSET 0 0 0\n CHANNELS 1 Xposition\n}\nMOTION\nFrames:\t2\nFrame   Time:\t0.0625\n1\n2\n";
        let bvh = parse_bvh(loose).unwrap();
        assert_eq!(bvh.frame_time, 0.0625);
        assert_eq!(bvh.motion, vec![vec![1.0], vec![2.0]]);
    }

    #[test]
    fn test_motion_row_with_wrong_length_is_rejected() {
        let text = single_joint_bvh("Zrotation Xrotation Yrotation", 3, "10 20");
        let message = parse_bvh(&text).unwrap_err();
        assert!(message.contains("2 values"), "{message}");
        assert!(message.contains("3 channels"), "{message}");
        assert!(message.contains("line 14"), "{message}");
        let text = single_joint_bvh("Zrotation Xrotation Yrotation", 3, "10 20 30 40");
        assert!(parse_bvh(&text).unwrap_err().contains("4 values"));
    }

    #[test]
    fn test_non_numeric_motion_value_is_rejected_instead_of_skipped() {
        let text = single_joint_bvh("Zrotation Xrotation Yrotation", 3, "10 twenty 30");
        let message = parse_bvh(&text).unwrap_err();
        assert!(message.contains("twenty"), "{message}");
    }

    #[test]
    fn test_non_finite_numbers_are_rejected() {
        for bad in ["NaN", "inf", "-inf"] {
            let text =
                single_joint_bvh("Zrotation Xrotation Yrotation", 3, &format!("10 {bad} 30"));
            let message = parse_bvh(&text).unwrap_err();
            assert!(message.contains(bad), "{message}");
            let offset = format!("HIERARCHY\nROOT R\n{{\n OFFSET 0 {bad} 0\n}}\n");
            let message = parse_bvh(&offset).unwrap_err();
            assert!(
                message.contains(bad) && message.contains("OFFSET"),
                "{message}"
            );
        }
    }

    #[test]
    fn test_malformed_header_lines_are_rejected() {
        let hierarchy = "HIERARCHY\nROOT R\n{\n";
        let unknown = format!("{hierarchy} OFFSET 0 0 0\n CHANNELS 1 Wrotation\n}}\n");
        assert!(parse_bvh(&unknown).unwrap_err().contains("Wrotation"));
        let count = format!("{hierarchy} OFFSET 0 0 0\n CHANNELS 3 Xposition Yposition\n}}\n");
        assert!(
            parse_bvh(&count)
                .unwrap_err()
                .contains("announces 3 channels but lists 2")
        );
        let offset = format!("{hierarchy} OFFSET 0 zero 0\n}}\n");
        assert!(parse_bvh(&offset).unwrap_err().contains("zero"));
        let short = format!("{hierarchy} OFFSET 0 0\n}}\n");
        assert!(parse_bvh(&short).unwrap_err().contains("three numbers"));
        let frames = format!("{hierarchy}}}\nMOTION\nFrames: many\n");
        assert!(parse_bvh(&frames).unwrap_err().contains("many"));
        let time = format!("{hierarchy}}}\nMOTION\nFrame Time: soon\n");
        assert!(parse_bvh(&time).unwrap_err().contains("soon"));
    }

    // ── rotation order ─────────────────────────────────────────────────────

    #[test]
    fn test_zxy_channels_compose_as_rz_rx_ry() {
        let text = single_joint_bvh("Zrotation Xrotation Yrotation", 3, "30 45 60");
        let expected = product([Axis::Z, Axis::X, Axis::Y], [30.0, 45.0, 60.0]);
        assert_matrix_close(&root_rotation(&text), &expected, 1e-12);
    }

    #[test]
    fn test_every_channel_order_composes_in_file_order() {
        for order in ALL_ORDERS {
            let text = single_joint_bvh(&channel_names(order), 3, "30 45 60");
            let expected = product(order, [30.0, 45.0, 60.0]);
            assert_matrix_close(&root_rotation(&text), &expected, 1e-12);
        }
    }

    #[test]
    fn test_channel_order_matters() {
        let zxy = root_rotation(&single_joint_bvh(
            "Zrotation Xrotation Yrotation",
            3,
            "30 45 60",
        ));
        let xyz = root_rotation(&single_joint_bvh(
            "Xrotation Yrotation Zrotation",
            3,
            "30 45 60",
        ));
        let difference: f64 = (0..3)
            .flat_map(|i| (0..3).map(move |j| (i, j)))
            .map(|(i, j)| (zxy[i][j] - xyz[i][j]).abs())
            .sum();
        assert!(difference > 0.1, "difference {difference}");
    }

    #[test]
    fn test_single_and_double_rotation_channels() {
        let text = single_joint_bvh("Yrotation", 1, "25");
        assert_matrix_close(&root_rotation(&text), &rot(Axis::Y, 25.0), 1e-12);
        let text = single_joint_bvh("Zrotation Xrotation", 2, "40 -10");
        let expected = mat_mul(&rot(Axis::Z, 40.0), &rot(Axis::X, -10.0));
        assert_matrix_close(&root_rotation(&text), &expected, 1e-12);
    }

    #[test]
    fn test_position_channels_add_to_the_offset() {
        let text = "HIERARCHY\nROOT R\n{\n OFFSET 1 2 3\n CHANNELS 2 Zposition Xposition\n}\nMOTION\nFrames: 1\nFrame Time: 0.1\n10 20\n";
        let bvh = parse_bvh(text).unwrap();
        let pose = bvh_frame_to_pose(&bvh.skeleton, &bvh.motion[0]);
        assert_eq!(pose.joint_poses[0].translation, [21.0, 2.0, 13.0]);
    }

    // ── Euler extraction ───────────────────────────────────────────────────

    fn quat_of(order: [Axis; 3], angles: [f64; 3]) -> [f64; 4] {
        order
            .iter()
            .zip(angles)
            .fold([0.0, 0.0, 0.0, 1.0], |q, (axis, degrees)| {
                quat_mul(q, axis_angle_quat(*axis, degrees))
            })
    }

    #[test]
    fn test_euler_extraction_inverts_the_composition_for_all_six_orders() {
        let outer = [-170.0, -100.0, -45.0, 0.0, 30.0, 90.0, 150.0, 180.0];
        let middle = [-85.0, -30.0, 0.0, 45.0, 80.0];
        for order in ALL_ORDERS {
            for &a in &outer {
                for &b in &middle {
                    for &c in &outer {
                        let q = quat_of(order, [a, b, c]);
                        let got = quat_to_euler(q, order);
                        // 180 and -180 designate the same angle.
                        let wrap = |x: f64, y: f64| {
                            let d = (x - y).rem_euclid(360.0);
                            d.min(360.0 - d)
                        };
                        assert!(
                            wrap(got[0], a) < 1e-8
                                && (got[1] - b).abs() < 1e-8
                                && wrap(got[2], c) < 1e-8,
                            "order {order:?}: input ({a}, {b}, {c}) gave {got:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_euler_extraction_reproduces_the_rotation_at_gimbal_lock() {
        for order in ALL_ORDERS {
            for middle in [90.0, -90.0] {
                for (a, c) in [(20.0, 35.0), (-70.0, 10.0), (0.0, 0.0), (120.0, -150.0)] {
                    let q = quat_of(order, [a, middle, c]);
                    let got = quat_to_euler(q, order);
                    assert!(
                        (got[1].abs() - 90.0).abs() < 1e-6,
                        "order {order:?}: middle angle {}",
                        got[1]
                    );
                    assert_eq!(got[2], 0.0, "order {order:?}");
                    let rebuilt = product(order, got);
                    let expected = product(order, [a, middle, c]);
                    assert_matrix_close(&rebuilt, &expected, 1e-9);
                }
            }
        }
    }

    // ── export ─────────────────────────────────────────────────────────────

    #[test]
    fn test_export_bvh_basic() {
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

        let exported = export_bvh(&skel, &clip).unwrap();
        assert!(exported.contains("ROOT Root"));
        assert!(exported.contains("MOTION"));
        assert!(exported.contains("Frames: 1"));
        assert!(exported.contains("Frame Time: 0.033333"));
        assert!(
            exported.ends_with("1.000000 2.000000 3.000000\n"),
            "{exported}"
        );
    }

    const CMU_LIKE: &str = "\
HIERARCHY
ROOT Hips
{
  OFFSET 0.0 0.0 0.0
  CHANNELS 6 Xposition Yposition Zposition Zrotation Xrotation Yrotation
  JOINT Chest
  {
    OFFSET 0.0 5.21 0.0
    CHANNELS 3 Zrotation Xrotation Yrotation
    JOINT Head
    {
      OFFSET 0.0 5.45 0.0
      CHANNELS 3 Zrotation Xrotation Yrotation
      End Site
      {
        OFFSET 0.0 3.87 0.0
      }
    }
  }
  JOINT LeftHip
  {
    OFFSET 3.91 0.0 0.0
    CHANNELS 3 Zrotation Yrotation Xrotation
    End Site
    {
      OFFSET 0.0 -10.5 0.0
    }
  }
}
MOTION
Frames: 2
Frame Time: 0.008333
1.5 -2.5 3.25 30.0 45.0 60.0 -20.0 10.0 5.0 15.0 -30.0 40.0 100.0 -45.0 -170.0
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0
";

    fn export_text(bvh: &BvhData) -> String {
        export_bvh(&bvh.skeleton, &bvh_to_clip(bvh)).unwrap()
    }

    #[test]
    fn test_parse_export_parse_round_trips() {
        let first = parse_bvh(CMU_LIKE).unwrap();
        let exported = export_text(&first);
        let second = parse_bvh(&exported).unwrap();

        assert_eq!(second.skeleton.num_joints(), first.skeleton.num_joints());
        for (a, b) in first.skeleton.joints.iter().zip(&second.skeleton.joints) {
            assert_eq!(a.name, b.name);
            assert_eq!(a.parent, b.parent);
            assert_eq!(a.channels, b.channels);
            assert_eq!(a.is_end_site, b.is_end_site);
            for (x, y) in a.offset.iter().zip(&b.offset) {
                assert!((x - y).abs() < 1e-6, "offset {x} vs {y}");
            }
        }
        assert!((second.frame_time - first.frame_time).abs() < 1e-6);
        assert_eq!(second.motion.len(), first.motion.len());
        for (row_a, row_b) in first.motion.iter().zip(&second.motion) {
            assert_eq!(row_a.len(), row_b.len());
            for (a, b) in row_a.iter().zip(row_b) {
                assert!((a - b).abs() < 1e-5, "{a} vs {b}\n{exported}");
            }
        }
    }

    #[test]
    fn test_second_export_is_identical_to_the_first() {
        let first = parse_bvh(CMU_LIKE).unwrap();
        let text_one = export_text(&first);
        let text_two = export_text(&parse_bvh(&text_one).unwrap());
        assert_eq!(text_one, text_two);
        // The end sites survive as End Site blocks, not as extra joints.
        assert_eq!(text_one.matches("End Site").count(), 2);
        assert_eq!(text_one.matches("JOINT").count(), 3);
    }

    #[test]
    fn test_round_trip_for_every_channel_order() {
        for order in ALL_ORDERS {
            let text = single_joint_bvh(&channel_names(order), 3, "30 45 60");
            let first = parse_bvh(&text).unwrap();
            let second = parse_bvh(&export_text(&first)).unwrap();
            for (a, b) in first.motion[0].iter().zip(&second.motion[0]) {
                assert!((a - b).abs() < 1e-5, "order {order:?}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn test_export_writes_translation_minus_offset() {
        let text = "HIERARCHY\nROOT R\n{\n OFFSET 1 2 3\n CHANNELS 3 Xposition Yposition Zposition\n}\nMOTION\nFrames: 1\nFrame Time: 0.1\n10 20 30\n";
        let first = parse_bvh(text).unwrap();
        let exported = export_text(&first);
        let second = parse_bvh(&exported).unwrap();
        assert_eq!(second.motion, vec![vec![10.0, 20.0, 30.0]]);
    }

    #[test]
    fn test_export_decomposes_in_the_channel_order_of_each_joint() {
        let mut skel = Skeleton::new();
        let root = skel.add_joint("Root", None, [0.0; 3]);
        skel.set_channels(
            root,
            vec!["Zrotation".into(), "Xrotation".into(), "Yrotation".into()],
        );
        let mut clip = AnimationClip::new("rot", 60.0);
        let mut pose = SkeletonPose::default_for(&skel);
        pose.joint_poses[0].rotation = quat_of([Axis::Z, Axis::X, Axis::Y], [30.0, 45.0, 60.0]);
        clip.add_frame(0.0, pose);
        let exported = export_bvh(&skel, &clip).unwrap();
        assert!(
            exported.ends_with("30.000000 45.000000 60.000000\n"),
            "{exported}"
        );
    }

    #[test]
    fn test_motion_values_follow_the_written_hierarchy_order() {
        // Joints are added breadth first: Root, A, B, then C below A.  The
        // hierarchy is written depth first (Root, A, C, B), so the motion
        // values must be in that order too.
        let mut skel = Skeleton::new();
        let root = skel.add_joint("Root", None, [0.0; 3]);
        let a = skel.add_joint("A", Some(root), [0.0; 3]);
        let b = skel.add_joint("B", Some(root), [0.0; 3]);
        let c = skel.add_joint("C", Some(a), [0.0; 3]);
        for joint in [root, a, b, c] {
            skel.set_channels(joint, vec!["Xposition".into()]);
        }
        let mut pose = SkeletonPose::default_for(&skel);
        for (joint, x) in [(root, 1.0), (a, 2.0), (b, 3.0), (c, 4.0)] {
            pose.joint_poses[joint].translation = [x, 0.0, 0.0];
        }
        let mut clip = AnimationClip::new("order", 30.0);
        clip.add_frame(0.0, pose);

        let exported = export_bvh(&skel, &clip).unwrap();
        let reparsed = parse_bvh(&exported).unwrap();
        let names: Vec<&str> = reparsed
            .skeleton
            .joints
            .iter()
            .filter(|j| !j.is_end_site)
            .map(|j| j.name.as_str())
            .collect();
        assert_eq!(names, vec!["Root", "A", "C", "B"]);
        assert_eq!(reparsed.motion, vec![vec![1.0, 2.0, 4.0, 3.0]]);
    }

    #[test]
    fn test_leaf_joint_without_end_site_gets_one() {
        let mut skel = Skeleton::new();
        let root = skel.add_joint("Root", None, [0.0; 3]);
        skel.set_channels(root, vec!["Xposition".into()]);
        let mut clip = AnimationClip::new("leaf", 30.0);
        clip.add_frame(0.0, SkeletonPose::default_for(&skel));
        let exported = export_bvh(&skel, &clip).unwrap();
        assert_eq!(exported.matches("End Site").count(), 1, "{exported}");
        let reparsed = parse_bvh(&exported).unwrap();
        assert_eq!(reparsed.skeleton.num_joints(), 2);
        assert!(reparsed.skeleton.joints[1].is_end_site);
    }

    #[test]
    fn test_missing_poses_are_an_error_instead_of_a_panic() {
        let mut skel = Skeleton::new();
        let root = skel.add_joint("Root", None, [0.0; 3]);
        skel.add_joint("Child", Some(root), [0.0, 1.0, 0.0]);
        let mut clip = AnimationClip::new("short", 30.0);
        clip.add_frame(
            0.0,
            SkeletonPose {
                joint_poses: vec![JointPose::identity()],
            },
        );
        let message = export_bvh(&skel, &clip).unwrap_err();
        assert!(message.contains("frame 0"), "{message}");
        assert!(message.contains("1 joint poses"), "{message}");
    }

    #[test]
    fn test_unexportable_skeletons_are_rejected() {
        let clip = AnimationClip::new("empty", 30.0);

        let mut cyclic = Skeleton::new();
        cyclic.add_joint("A", Some(1), [0.0; 3]);
        cyclic.add_joint("B", Some(0), [0.0; 3]);
        assert!(
            export_bvh(&cyclic, &clip)
                .unwrap_err()
                .contains("not reachable")
        );

        let mut repeated = Skeleton::new();
        let r = repeated.add_joint("R", None, [0.0; 3]);
        repeated.set_channels(
            r,
            vec!["Zrotation".into(), "Xrotation".into(), "Zrotation".into()],
        );
        assert!(export_bvh(&repeated, &clip).unwrap_err().contains("twice"));

        let mut unknown = Skeleton::new();
        let r = unknown.add_joint("R", None, [0.0; 3]);
        unknown.set_channels(r, vec!["Wrotation".into()]);
        assert!(
            export_bvh(&unknown, &clip)
                .unwrap_err()
                .contains("Wrotation")
        );

        let mut spaced = Skeleton::new();
        spaced.add_joint("left arm", None, [0.0; 3]);
        assert!(
            export_bvh(&spaced, &clip)
                .unwrap_err()
                .contains("whitespace")
        );

        let mut bad_end = Skeleton::new();
        let r = bad_end.add_joint("R", None, [0.0; 3]);
        let e = bad_end.add_end_site(Some(r), [0.0; 3]);
        bad_end.add_joint("Below", Some(e), [0.0; 3]);
        assert!(
            export_bvh(&bad_end, &clip)
                .unwrap_err()
                .contains("End Site")
        );
    }
}
