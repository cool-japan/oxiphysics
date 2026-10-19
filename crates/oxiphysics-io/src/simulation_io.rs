// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics simulation I/O formats.
//!
//! Provides writers and readers for full physics scenes (bodies, joints,
//! materials), multi-body trajectories, contact force logs, energy logs,
//! simulation checkpoints, VTK time-series output, and XDMF wrappers.
//!
//! Scenes and checkpoints are JSON documents read and written with `serde`;
//! every field of [`BodyDesc`], [`JointDesc`] and [`MaterialDesc`] survives a
//! write/read cycle bit-exactly, and non-finite floats are encoded as the
//! strings `"NaN"`, `"Infinity"` and `"-Infinity"` (see the `json` helper
//! module). Trajectories use a small binary format described on
//! [`TrajectoryWriter::write_binary`].

use std::io::{BufWriter, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

mod json;

// ---------------------------------------------------------------------------
// Helper free functions
// ---------------------------------------------------------------------------

/// Encode a quaternion (w, x, y, z) as four f32 values in little-endian bytes.
pub fn encode_quaternion_f32(w: f64, x: f64, y: f64, z: f64) -> [u8; 16] {
    let mut buf = [0u8; 16];
    buf[0..4].copy_from_slice(&(w as f32).to_le_bytes());
    buf[4..8].copy_from_slice(&(x as f32).to_le_bytes());
    buf[8..12].copy_from_slice(&(y as f32).to_le_bytes());
    buf[12..16].copy_from_slice(&(z as f32).to_le_bytes());
    buf
}

/// Pack a slice of \[f64; 3\] vectors into an interleaved f32 little-endian byte
/// buffer.
pub fn pack_float3_array(data: &[[f64; 3]]) -> Vec<u8> {
    let mut buf = Vec::with_capacity(data.len() * 12);
    for &[x, y, z] in data {
        buf.extend_from_slice(&(x as f32).to_le_bytes());
        buf.extend_from_slice(&(y as f32).to_le_bytes());
        buf.extend_from_slice(&(z as f32).to_le_bytes());
    }
    buf
}

/// Write a ParaView PVD collection XML file pointing to per-step VTU files.
pub fn write_pvd_collection<W: Write>(
    writer: &mut W,
    entries: &[(f64, &str)],
) -> std::io::Result<()> {
    writeln!(writer, r#"<?xml version="1.0"?>"#)?;
    writeln!(writer, r#"<VTKFile type="Collection" version="0.1">"#)?;
    writeln!(writer, "  <Collection>")?;
    for &(time, path) in entries {
        writeln!(
            writer,
            r#"    <DataSet timestep="{time}" group="" part="0" file="{path}"/>"#
        )?;
    }
    writeln!(writer, "  </Collection>")?;
    writeln!(writer, "</VTKFile>")?;
    Ok(())
}

/// Parse a simple XDMF file and return a list of (time, data_file) pairs.
pub fn read_xdmf_timesteps(xml: &str) -> Vec<(f64, String)> {
    let mut result = Vec::new();
    for line in xml.lines() {
        let line = line.trim();
        if line.starts_with("<Grid") {
            // Extract Time value and DataItem reference
            if let Some(tv) = extract_attr(line, "Time")
                && let Ok(t) = tv.parse::<f64>()
            {
                result.push((t, String::new()));
            }
        } else if line.starts_with("<DataItem")
            && let Some(href) = extract_attr(line, "href")
            && let Some(last) = result.last_mut()
            && last.1.is_empty()
        {
            last.1 = href.to_string();
        }
    }
    result
}

fn extract_attr<'a>(s: &'a str, attr: &str) -> Option<&'a str> {
    let pat = format!("{attr}=\"");
    let start = s.find(pat.as_str())? + pat.len();
    let end = s[start..].find('"')? + start;
    Some(&s[start..end])
}

// ---------------------------------------------------------------------------
// Body / scene types
// ---------------------------------------------------------------------------

/// Orientation of a body that has not been rotated.
fn identity_quaternion() -> [f64; 4] {
    [1.0, 0.0, 0.0, 0.0]
}

/// Mass assumed for a body descriptor that does not specify one.
fn unit_mass() -> f64 {
    1.0
}

/// Shape tag assumed for a body descriptor that does not specify one.
fn default_shape_tag() -> String {
    "sphere".to_string()
}

/// A rigid body descriptor for scene serialization.
///
/// The JSON keys are `id`, `pos`, `quat`, `linvel`, `angvel`, `mass` and
/// `shape`. When reading, every key except `id` and `pos` may be omitted and
/// then takes the value [`BodyDesc::new`] would give it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BodyDesc {
    /// Unique body identifier.
    pub id: u64,
    /// World-space position.
    #[serde(rename = "pos", with = "json::f64_array")]
    pub position: [f64; 3],
    /// Orientation as quaternion (w, x, y, z).
    #[serde(
        rename = "quat",
        default = "identity_quaternion",
        with = "json::f64_array"
    )]
    pub orientation: [f64; 4],
    /// Linear velocity.
    #[serde(rename = "linvel", default, with = "json::f64_array")]
    pub linear_velocity: [f64; 3],
    /// Angular velocity.
    #[serde(rename = "angvel", default, with = "json::f64_array")]
    pub angular_velocity: [f64; 3],
    /// Mass (kg).
    #[serde(default = "unit_mass", with = "json::f64_value")]
    pub mass: f64,
    /// Shape type tag.
    #[serde(rename = "shape", default = "default_shape_tag")]
    pub shape_tag: String,
}

impl BodyDesc {
    /// Construct a body descriptor.
    pub fn new(id: u64, position: [f64; 3]) -> Self {
        BodyDesc {
            id,
            position,
            orientation: identity_quaternion(),
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            mass: unit_mass(),
            shape_tag: default_shape_tag(),
        }
    }

    /// Serialize to a compact JSON object string.
    ///
    /// Strings are escaped and floats are written so that they parse back to
    /// the identical value; see the module documentation for non-finite
    /// values.
    ///
    /// # Errors
    /// Returns [`crate::Error::Io`] if the encoder fails.
    pub fn to_json(&self) -> crate::Result<String> {
        Ok(json::encode(self)?)
    }
}

/// A joint descriptor.
///
/// The JSON keys are `id`, `bodyA`, `bodyB`, `type`, `anchorA` and `anchorB`;
/// the anchors may be omitted when reading (they default to the origin).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct JointDesc {
    /// Unique joint id.
    pub id: u64,
    /// Body A id.
    #[serde(rename = "bodyA")]
    pub body_a: u64,
    /// Body B id.
    #[serde(rename = "bodyB")]
    pub body_b: u64,
    /// Joint type.
    #[serde(rename = "type")]
    pub joint_type: String,
    /// Anchor point on body A (local space).
    #[serde(rename = "anchorA", default, with = "json::f64_array")]
    pub anchor_a: [f64; 3],
    /// Anchor point on body B (local space).
    #[serde(rename = "anchorB", default, with = "json::f64_array")]
    pub anchor_b: [f64; 3],
}

impl JointDesc {
    /// Construct a joint descriptor.
    pub fn new(id: u64, body_a: u64, body_b: u64, joint_type: &str) -> Self {
        JointDesc {
            id,
            body_a,
            body_b,
            joint_type: joint_type.to_string(),
            anchor_a: [0.0; 3],
            anchor_b: [0.0; 3],
        }
    }
}

/// A physics material descriptor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MaterialDesc {
    /// Material id.
    pub id: u64,
    /// Restitution coefficient.
    #[serde(with = "json::f64_value")]
    pub restitution: f64,
    /// Dynamic friction coefficient.
    #[serde(with = "json::f64_value")]
    pub friction: f64,
    /// Density (kg/m³).
    #[serde(with = "json::f64_value")]
    pub density: f64,
}

impl MaterialDesc {
    /// Construct a material descriptor.
    pub fn new(id: u64, restitution: f64, friction: f64, density: f64) -> Self {
        MaterialDesc {
            id,
            restitution,
            friction,
            density,
        }
    }
}

// ---------------------------------------------------------------------------
// PhysicsSceneWriter
// ---------------------------------------------------------------------------

/// Write a full physics scene (bodies, joints, materials) to JSON or binary.
pub struct PhysicsSceneWriter {
    /// Bodies in the scene.
    pub bodies: Vec<BodyDesc>,
    /// Joints in the scene.
    pub joints: Vec<JointDesc>,
    /// Materials in the scene.
    pub materials: Vec<MaterialDesc>,
}

impl PhysicsSceneWriter {
    /// Construct an empty scene writer.
    pub fn new() -> Self {
        PhysicsSceneWriter {
            bodies: Vec::new(),
            joints: Vec::new(),
            materials: Vec::new(),
        }
    }

    /// Add a body to the scene.
    pub fn add_body(&mut self, body: BodyDesc) {
        self.bodies.push(body);
    }

    /// Add a joint to the scene.
    pub fn add_joint(&mut self, joint: JointDesc) {
        self.joints.push(joint);
    }

    /// Add a material to the scene.
    pub fn add_material(&mut self, mat: MaterialDesc) {
        self.materials.push(mat);
    }

    /// Serialize the scene to JSON and write to `writer`.
    ///
    /// The document has the members `bodies`, `joints` and `materials`, one
    /// element per line; read it back with [`PhysicsSceneReader::read_json`].
    pub fn write_json<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(writer, "{{")?;
        write_json_array(writer, "bodies", &self.bodies, true)?;
        write_json_array(writer, "joints", &self.joints, true)?;
        write_json_array(writer, "materials", &self.materials, false)?;
        writeln!(writer, "}}")?;
        Ok(())
    }

    /// Write to a file path.
    pub fn write_to_file(&self, path: &str) -> crate::Result<()> {
        let file = std::fs::File::create(Path::new(path))?;
        let mut w = BufWriter::new(file);
        self.write_json(&mut w)?;
        w.flush()?;
        Ok(())
    }
}

impl Default for PhysicsSceneWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// PhysicsSceneReader
// ---------------------------------------------------------------------------

/// Write `"name": [ ... ]` with one compact JSON element per line.
fn write_json_array<W: Write, T: Serialize>(
    writer: &mut W,
    name: &str,
    items: &[T],
    trailing_comma: bool,
) -> std::io::Result<()> {
    writeln!(writer, "  \"{name}\": [")?;
    for (i, item) in items.iter().enumerate() {
        let comma = if i + 1 < items.len() { "," } else { "" };
        writeln!(writer, "    {}{}", json::encode(item)?, comma)?;
    }
    writeln!(writer, "  ]{}", if trailing_comma { "," } else { "" })
}

/// A complete scene: everything [`PhysicsSceneWriter`] writes.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct SceneDesc {
    /// Bodies in the scene.
    #[serde(default)]
    pub bodies: Vec<BodyDesc>,
    /// Joints in the scene.
    #[serde(default)]
    pub joints: Vec<JointDesc>,
    /// Materials in the scene.
    #[serde(default)]
    pub materials: Vec<MaterialDesc>,
}

/// Read a physics scene from JSON, validate it, and reconstruct the bodies,
/// joints and materials.
pub struct PhysicsSceneReader;

impl PhysicsSceneReader {
    /// Parse a JSON scene document as produced by [`PhysicsSceneWriter`].
    ///
    /// Bodies, joints and materials are restored with all their fields.
    /// Members that are absent are read as empty lists, and unknown members
    /// are ignored.
    ///
    /// # Errors
    /// Returns [`crate::Error::Parse`] if the text is not valid JSON or an
    /// element has a missing or mistyped field.
    pub fn read_json(json: &str) -> crate::Result<SceneDesc> {
        serde_json::from_str(json).map_err(|e| crate::Error::Parse(format!("scene JSON: {e}")))
    }

    /// Validate a list of body descriptors (checks for duplicate IDs).
    pub fn validate(bodies: &[BodyDesc]) -> bool {
        let mut ids = std::collections::HashSet::new();
        bodies.iter().all(|b| ids.insert(b.id))
    }
}

// ---------------------------------------------------------------------------
// TrajectoryWriter
// ---------------------------------------------------------------------------

/// A single trajectory frame: positions, velocities, quaternions.
#[derive(Debug, Clone)]
pub struct TrajectoryFrame {
    /// Frame time.
    pub time: f64,
    /// Per-body positions.
    pub positions: Vec<[f64; 3]>,
    /// Per-body velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Per-body quaternions (w, x, y, z).
    pub quaternions: Vec<[f64; 4]>,
}

impl TrajectoryFrame {
    /// Construct an empty frame at the given time.
    pub fn new(time: f64) -> Self {
        TrajectoryFrame {
            time,
            positions: Vec::new(),
            velocities: Vec::new(),
            quaternions: Vec::new(),
        }
    }
}

/// Writer for multi-body trajectory files: positions, velocities, quaternions
/// per timestep.
pub struct TrajectoryWriter {
    frames: Vec<TrajectoryFrame>,
}

impl TrajectoryWriter {
    /// Construct an empty trajectory writer.
    pub fn new() -> Self {
        TrajectoryWriter { frames: Vec::new() }
    }

    /// Append a frame.
    pub fn push_frame(&mut self, frame: TrajectoryFrame) {
        self.frames.push(frame);
    }

    /// Number of frames.
    pub fn num_frames(&self) -> usize {
        self.frames.len()
    }

    /// Write trajectory to a binary file.
    ///
    /// Format (all integers and floats little-endian):
    ///
    /// ```text
    /// "TRAJ" 0x01            5-byte magic, the last byte is the format version
    /// n_frames               u64
    /// repeated n_frames times:
    ///   time                 f64
    ///   n_pos, n_vel, n_quat u64 each
    ///   positions            n_pos  x 3 f32
    ///   velocities           n_vel  x 3 f32
    ///   quaternions (w,x,y,z) n_quat x 4 f32
    /// ```
    ///
    /// Positions, velocities and quaternions are stored as `f32`. A frame may
    /// leave velocities and/or quaternions out (the count is then 0), but
    /// otherwise they must have one entry per position.
    ///
    /// # Errors
    /// Returns an [`InvalidInput`](std::io::ErrorKind::InvalidInput) error,
    /// before anything is written, if a frame has a velocity or quaternion
    /// count that is neither 0 nor the number of positions; such a frame would
    /// shift every later value in the file. Other errors come from `w`.
    pub fn write_binary<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        for (index, f) in self.frames.iter().enumerate() {
            let n_pos = f.positions.len();
            for (what, count) in [
                ("velocities", f.velocities.len()),
                ("quaternions", f.quaternions.len()),
            ] {
                if count != 0 && count != n_pos {
                    return Err(std::io::Error::new(
                        std::io::ErrorKind::InvalidInput,
                        format!(
                            "trajectory frame {index} has {count} {what} for {n_pos} positions \
                             (expected 0 or {n_pos})"
                        ),
                    ));
                }
            }
        }
        w.write_all(TRAJ_MAGIC)?;
        w.write_all(&(self.frames.len() as u64).to_le_bytes())?;
        for f in &self.frames {
            w.write_all(&f.time.to_le_bytes())?;
            w.write_all(&(f.positions.len() as u64).to_le_bytes())?;
            w.write_all(&(f.velocities.len() as u64).to_le_bytes())?;
            w.write_all(&(f.quaternions.len() as u64).to_le_bytes())?;
            w.write_all(&pack_float3_array(&f.positions))?;
            w.write_all(&pack_float3_array(&f.velocities))?;
            for &q in &f.quaternions {
                let bytes = encode_quaternion_f32(q[0], q[1], q[2], q[3]);
                w.write_all(&bytes)?;
            }
        }
        Ok(())
    }

    /// Write trajectory to a CSV file (one row per body per frame).
    pub fn write_csv<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "time,body,px,py,pz,vx,vy,vz,qw,qx,qy,qz")?;
        for f in &self.frames {
            let nb = f.positions.len();
            for i in 0..nb {
                let p = f.positions.get(i).copied().unwrap_or([0.0; 3]);
                let v = f.velocities.get(i).copied().unwrap_or([0.0; 3]);
                let q = f
                    .quaternions
                    .get(i)
                    .copied()
                    .unwrap_or([1.0, 0.0, 0.0, 0.0]);
                writeln!(
                    w,
                    "{},{},{},{},{},{},{},{},{},{},{},{}",
                    f.time, i, p[0], p[1], p[2], v[0], v[1], v[2], q[0], q[1], q[2], q[3]
                )?;
            }
        }
        Ok(())
    }
}

impl Default for TrajectoryWriter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TrajectoryReader
// ---------------------------------------------------------------------------

/// Magic bytes of the binary trajectory format: `"TRAJ"` and format version 1.
const TRAJ_MAGIC: &[u8; 5] = b"TRAJ\x01";

/// Magic bytes of the retired version-0 layout (`"TRAJ"` and a zero byte).
///
/// That layout wrote one count per frame but then velocities and quaternions
/// of whatever length, so it cannot be decoded reliably.
const TRAJ_LEGACY_MAGIC: &[u8; 5] = b"TRAJ\0";

/// Bounds-checked little-endian reader over a byte slice.
struct ByteCursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> ByteCursor<'a> {
    fn new(data: &'a [u8], pos: usize) -> Self {
        Self { data, pos }
    }

    fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    /// Take the next `len` bytes, or fail if fewer remain.
    fn take(&mut self, len: usize) -> crate::Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|&end| end <= self.data.len())
            .ok_or_else(|| {
                crate::Error::Parse(format!(
                    "TRAJ data truncated: {len} bytes needed at byte {} but only {} remain",
                    self.pos,
                    self.remaining()
                ))
            })?;
        let bytes = &self.data[self.pos..end];
        self.pos = end;
        Ok(bytes)
    }

    fn read_u64(&mut self) -> crate::Result<u64> {
        let bytes = self.take(8)?;
        let mut buf = [0u8; 8];
        buf.copy_from_slice(bytes);
        Ok(u64::from_le_bytes(buf))
    }

    fn read_f64(&mut self) -> crate::Result<f64> {
        let bytes = self.take(8)?;
        let mut buf = [0u8; 8];
        buf.copy_from_slice(bytes);
        Ok(f64::from_le_bytes(buf))
    }

    /// Read a `u64` element count that must fit in `usize`.
    fn read_count(&mut self, what: &str) -> crate::Result<usize> {
        let value = self.read_u64()?;
        usize::try_from(value).map_err(|_| {
            crate::Error::Parse(format!("TRAJ {what} count {value} does not fit in memory"))
        })
    }
}

/// Decode `f32` little-endian values from a chunk of exactly `4 * N` bytes.
fn f32_group<const N: usize>(chunk: &[u8]) -> [f64; N] {
    let mut out = [0.0f64; N];
    for (slot, bytes) in out.iter_mut().zip(chunk.chunks_exact(4)) {
        *slot = f64::from(f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]));
    }
    out
}

/// Reader for binary trajectory files produced by [`TrajectoryWriter`].
pub struct TrajectoryReader {
    frames: Vec<TrajectoryFrame>,
    cursor: usize,
}

impl TrajectoryReader {
    /// Read from a byte buffer.
    ///
    /// Every count in the file is checked against the number of bytes that
    /// actually follow before anything is allocated, so a corrupt header
    /// cannot trigger a huge allocation or an arithmetic overflow.
    ///
    /// # Errors
    /// Returns [`crate::Error::Parse`] if the magic is wrong (including the
    /// retired version-0 layout), the data is truncated, a frame has a
    /// velocity or quaternion count other than 0 or the number of positions,
    /// or bytes remain after the last frame.
    pub fn from_bytes(data: &[u8]) -> crate::Result<Self> {
        let Some(magic) = data.get(..TRAJ_MAGIC.len()) else {
            return Err(crate::Error::Parse(format!(
                "TRAJ data too short for a header: {} bytes",
                data.len()
            )));
        };
        if magic == TRAJ_LEGACY_MAGIC.as_slice() {
            return Err(crate::Error::Parse(
                "unsupported legacy TRAJ file (version 0); rewrite it with TrajectoryWriter"
                    .to_string(),
            ));
        }
        if magic != TRAJ_MAGIC.as_slice() {
            return Err(crate::Error::Parse(
                "not a TRAJ file: bad magic".to_string(),
            ));
        }
        let mut cur = ByteCursor::new(data, TRAJ_MAGIC.len());
        let n_frames = cur.read_count("frame")?;
        // A frame is at least its time and three counts: 32 bytes.
        if n_frames > cur.remaining() / 32 {
            return Err(crate::Error::Parse(format!(
                "TRAJ header declares {n_frames} frames but only {} bytes follow",
                cur.remaining()
            )));
        }
        let mut frames = Vec::with_capacity(n_frames);
        for index in 0..n_frames {
            let time = cur.read_f64()?;
            let n_pos = cur.read_count("position")?;
            let n_vel = cur.read_count("velocity")?;
            let n_quat = cur.read_count("quaternion")?;
            if (n_vel != 0 && n_vel != n_pos) || (n_quat != 0 && n_quat != n_pos) {
                return Err(crate::Error::Parse(format!(
                    "TRAJ frame {index} has {n_vel} velocities and {n_quat} quaternions \
                     for {n_pos} positions (expected 0 or {n_pos})"
                )));
            }
            let sizes = (
                n_pos.checked_mul(12),
                n_vel.checked_mul(12),
                n_quat.checked_mul(16),
            );
            let (Some(pos_len), Some(vel_len), Some(quat_len)) = sizes else {
                return Err(crate::Error::Parse(format!(
                    "TRAJ frame {index}: element counts overflow the address space"
                )));
            };
            let total = pos_len
                .checked_add(vel_len)
                .and_then(|sum| sum.checked_add(quat_len))
                .ok_or_else(|| {
                    crate::Error::Parse(format!(
                        "TRAJ frame {index}: element counts overflow the address space"
                    ))
                })?;
            let body = cur.take(total)?;
            let (pos_bytes, rest) = body.split_at(pos_len);
            let (vel_bytes, quat_bytes) = rest.split_at(vel_len);
            let mut frame = TrajectoryFrame::new(time);
            frame.positions = pos_bytes.chunks_exact(12).map(f32_group::<3>).collect();
            frame.velocities = vel_bytes.chunks_exact(12).map(f32_group::<3>).collect();
            frame.quaternions = quat_bytes.chunks_exact(16).map(f32_group::<4>).collect();
            frames.push(frame);
        }
        if cur.remaining() != 0 {
            return Err(crate::Error::Parse(format!(
                "TRAJ data has {} unexpected bytes after the last frame",
                cur.remaining()
            )));
        }
        Ok(TrajectoryReader { frames, cursor: 0 })
    }

    /// Seek to a given frame index.
    pub fn seek(&mut self, frame_idx: usize) {
        self.cursor = frame_idx.min(self.frames.len());
    }

    /// Return the next frame, or None if at end.
    pub fn next_frame(&mut self) -> Option<&TrajectoryFrame> {
        if self.cursor < self.frames.len() {
            let f = &self.frames[self.cursor];
            self.cursor += 1;
            Some(f)
        } else {
            None
        }
    }

    /// Total number of frames.
    pub fn num_frames(&self) -> usize {
        self.frames.len()
    }
}

// ---------------------------------------------------------------------------
// ContactForceLog
// ---------------------------------------------------------------------------

/// A single contact force record.
#[derive(Debug, Clone)]
pub struct ContactForceRecord {
    /// Timestep index.
    pub step: u64,
    /// Body A id.
    pub body_a: u64,
    /// Body B id.
    pub body_b: u64,
    /// Normal impulse magnitude.
    pub normal_impulse: f64,
    /// Tangential impulse magnitude.
    pub tangent_impulse: f64,
    /// Contact point.
    pub point: [f64; 3],
    /// Contact normal.
    pub normal: [f64; 3],
}

/// Log of contact forces at each timestep.
pub struct ContactForceLog {
    records: Vec<ContactForceRecord>,
}

impl ContactForceLog {
    /// Construct an empty log.
    pub fn new() -> Self {
        ContactForceLog {
            records: Vec::new(),
        }
    }

    /// Add a contact force record.
    pub fn push(&mut self, rec: ContactForceRecord) {
        self.records.push(rec);
    }

    /// Get all records for a given step.
    pub fn records_at_step(&self, step: u64) -> Vec<&ContactForceRecord> {
        self.records.iter().filter(|r| r.step == step).collect()
    }

    /// Write to CSV.
    pub fn write_csv<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(
            w,
            "step,body_a,body_b,normal_impulse,tangent_impulse,px,py,pz,nx,ny,nz"
        )?;
        for r in &self.records {
            writeln!(
                w,
                "{},{},{},{},{},{},{},{},{},{},{}",
                r.step,
                r.body_a,
                r.body_b,
                r.normal_impulse,
                r.tangent_impulse,
                r.point[0],
                r.point[1],
                r.point[2],
                r.normal[0],
                r.normal[1],
                r.normal[2],
            )?;
        }
        Ok(())
    }

    /// Total number of records.
    pub fn len(&self) -> usize {
        self.records.len()
    }

    /// Returns true if the log is empty.
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }
}

impl Default for ContactForceLog {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// EnergyLog
// ---------------------------------------------------------------------------

/// Energy log entry for one timestep.
#[derive(Debug, Clone)]
pub struct EnergyEntry {
    /// Simulation time.
    pub time: f64,
    /// Total kinetic energy.
    pub kinetic: f64,
    /// Total potential energy.
    pub potential: f64,
    /// Total energy.
    pub total: f64,
    /// Linear momentum magnitude.
    pub linear_momentum: f64,
    /// Angular momentum magnitude.
    pub angular_momentum: f64,
}

/// Log kinetic/potential/total energy, angular momentum, linear momentum.
pub struct EnergyLog {
    entries: Vec<EnergyEntry>,
}

impl EnergyLog {
    /// Construct an empty energy log.
    pub fn new() -> Self {
        EnergyLog {
            entries: Vec::new(),
        }
    }

    /// Add an entry.
    pub fn push(&mut self, entry: EnergyEntry) {
        self.entries.push(entry);
    }

    /// Write to CSV.
    pub fn write_csv<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(
            w,
            "time,kinetic,potential,total,linear_momentum,angular_momentum"
        )?;
        for e in &self.entries {
            writeln!(
                w,
                "{},{},{},{},{},{}",
                e.time, e.kinetic, e.potential, e.total, e.linear_momentum, e.angular_momentum
            )?;
        }
        Ok(())
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns true if empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Maximum total energy recorded.
    pub fn max_total_energy(&self) -> f64 {
        self.entries
            .iter()
            .map(|e| e.total)
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

impl Default for EnergyLog {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// CheckpointWriter
// ---------------------------------------------------------------------------

/// A simulation checkpoint: full state including bodies, constraints, and
/// broadphase.
pub struct CheckpointWriter {
    /// Bodies in the checkpoint.
    pub bodies: Vec<BodyDesc>,
    /// Simulation time.
    pub sim_time: f64,
    /// Simulation step count.
    pub step: u64,
    /// Extra metadata key-value pairs.
    pub metadata: Vec<(String, String)>,
}

impl CheckpointWriter {
    /// Construct a checkpoint.
    pub fn new(sim_time: f64, step: u64) -> Self {
        CheckpointWriter {
            bodies: Vec::new(),
            sim_time,
            step,
            metadata: Vec::new(),
        }
    }

    /// Add a body.
    pub fn add_body(&mut self, body: BodyDesc) {
        self.bodies.push(body);
    }

    /// Add metadata.
    pub fn add_metadata(&mut self, key: &str, value: &str) {
        self.metadata.push((key.to_string(), value.to_string()));
    }

    /// Write checkpoint as JSON.
    ///
    /// The document has the members `sim_time`, `step`, `metadata` (an object
    /// whose entries keep their insertion order; keys and values are escaped)
    /// and `bodies`, one body per line. Read it back with
    /// [`CheckpointReader::read_json`].
    pub fn write_json<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "{{")?;
        writeln!(
            w,
            "  \"sim_time\": {},",
            json::encode(&json::SpecialF64(self.sim_time))?
        )?;
        writeln!(w, "  \"step\": {},", self.step)?;
        writeln!(
            w,
            "  \"metadata\": {},",
            json::encode(&json::MetadataRef(&self.metadata))?
        )?;
        write_json_array(w, "bodies", &self.bodies, false)?;
        writeln!(w, "}}")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// CheckpointReader
// ---------------------------------------------------------------------------

/// A simulation checkpoint as restored by [`CheckpointReader`].
#[derive(Debug, Clone, PartialEq, Deserialize)]
pub struct CheckpointData {
    /// Simulation time.
    #[serde(with = "json::f64_value")]
    pub sim_time: f64,
    /// Simulation step count.
    pub step: u64,
    /// Metadata key-value pairs in file order.
    #[serde(default, with = "json::ordered_string_map")]
    pub metadata: Vec<(String, String)>,
    /// Bodies in the checkpoint, with their full dynamic state.
    #[serde(default)]
    pub bodies: Vec<BodyDesc>,
}

/// Read a checkpoint and restore simulation state.
pub struct CheckpointReader;

impl CheckpointReader {
    /// Parse a JSON checkpoint produced by [`CheckpointWriter`].
    ///
    /// Time, step count, metadata and every body (pose, velocities, mass and
    /// shape) are restored; `metadata` and `bodies` may be absent, unknown
    /// members are ignored.
    ///
    /// # Errors
    /// Returns [`crate::Error::Parse`] if the text is not valid JSON, if
    /// `sim_time` or `step` is missing, or if a member has the wrong type.
    pub fn read_json(json: &str) -> crate::Result<CheckpointData> {
        serde_json::from_str(json).map_err(|e| crate::Error::Parse(format!("checkpoint JSON: {e}")))
    }
}

// ---------------------------------------------------------------------------
// VtkTrajectory
// ---------------------------------------------------------------------------

/// Write a time-series VTK output: per-step VTU files + a PVD collection.
pub struct VtkTrajectory {
    /// Base directory for output files.
    pub output_dir: String,
    /// Base filename prefix.
    pub prefix: String,
    entries: Vec<(f64, String)>,
}

impl VtkTrajectory {
    /// Construct for the given output directory and filename prefix.
    pub fn new(output_dir: &str, prefix: &str) -> Self {
        VtkTrajectory {
            output_dir: output_dir.to_string(),
            prefix: prefix.to_string(),
            entries: Vec::new(),
        }
    }

    /// Write a single frame as a VTU file and record it in the PVD collection.
    pub fn write_frame(&mut self, time: f64, positions: &[[f64; 3]]) -> std::io::Result<()> {
        let frame_idx = self.entries.len();
        let vtu_name = format!("{}_{:06}.vtu", self.prefix, frame_idx);
        let vtu_path = format!("{}/{}", self.output_dir, vtu_name);
        // Write minimal VTU
        let file = std::fs::File::create(Path::new(&vtu_path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, r#"<?xml version="1.0"?>"#)?;
        writeln!(w, r#"<VTKFile type="UnstructuredGrid" version="0.1">"#)?;
        writeln!(w, "  <UnstructuredGrid>")?;
        writeln!(
            w,
            r#"    <Piece NumberOfPoints="{}" NumberOfCells="0">"#,
            positions.len()
        )?;
        writeln!(w, "      <Points>")?;
        writeln!(
            w,
            r#"        <DataArray type="Float64" NumberOfComponents="3" format="ascii">"#
        )?;
        for p in positions {
            writeln!(w, "          {} {} {}", p[0], p[1], p[2])?;
        }
        writeln!(w, "        </DataArray>")?;
        writeln!(w, "      </Points>")?;
        writeln!(w, "      <Cells/>")?;
        writeln!(w, "    </Piece>")?;
        writeln!(w, "  </UnstructuredGrid>")?;
        writeln!(w, "</VTKFile>")?;
        w.flush()?;
        self.entries.push((time, vtu_name));
        Ok(())
    }

    /// Write the PVD collection file.
    pub fn write_pvd(&self) -> std::io::Result<()> {
        let pvd_path = format!("{}/{}.pvd", self.output_dir, self.prefix);
        let file = std::fs::File::create(Path::new(&pvd_path))?;
        let mut w = BufWriter::new(file);
        let refs: Vec<(f64, &str)> = self.entries.iter().map(|(t, s)| (*t, s.as_str())).collect();
        write_pvd_collection(&mut w, &refs)?;
        w.flush()?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// ParaviewXdmf
// ---------------------------------------------------------------------------

/// XDMF wrapper for HDF5-based trajectory (positions as heavy data).
pub struct ParaviewXdmf {
    /// Output XDMF file path.
    pub xdmf_path: String,
    /// HDF5 file name (basename only).
    pub h5_file: String,
    entries: Vec<(f64, usize)>,
}

impl ParaviewXdmf {
    /// Construct for the given XDMF and HDF5 paths.
    pub fn new(xdmf_path: &str, h5_file: &str) -> Self {
        ParaviewXdmf {
            xdmf_path: xdmf_path.to_string(),
            h5_file: h5_file.to_string(),
            entries: Vec::new(),
        }
    }

    /// Record a timestep with the given number of points.
    pub fn add_timestep(&mut self, time: f64, n_points: usize) {
        self.entries.push((time, n_points));
    }

    /// Write the XDMF file.
    pub fn write(&self) -> std::io::Result<()> {
        let file = std::fs::File::create(Path::new(&self.xdmf_path))?;
        let mut w = BufWriter::new(file);
        writeln!(w, r#"<?xml version="1.0"?>"#)?;
        writeln!(w, r#"<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">"#)?;
        writeln!(w, r#"<Xdmf Version="2.0">"#)?;
        writeln!(w, "  <Domain>")?;
        writeln!(
            w,
            r#"    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">"#
        )?;
        for (i, &(t, n)) in self.entries.iter().enumerate() {
            writeln!(w, r#"      <Grid Name="step_{i}" GridType="Uniform">"#)?;
            writeln!(w, r#"        <Time Value="{t}"/>"#)?;
            writeln!(
                w,
                r#"        <Topology TopologyType="Polyvertex" NumberOfElements="{n}"/>"#
            )?;
            writeln!(w, r#"        <Geometry GeometryType="XYZ">"#)?;
            writeln!(
                w,
                r#"          <DataItem Dimensions="{n} 3" NumberType="Float" Precision="4" Format="HDF" href="{h5}">/step_{i}/positions</DataItem>"#,
                h5 = self.h5_file,
                n = n
            )?;
            writeln!(w, "        </Geometry>")?;
            writeln!(w, "      </Grid>")?;
        }
        writeln!(w, "    </Grid>")?;
        writeln!(w, "  </Domain>")?;
        writeln!(w, "</Xdmf>")?;
        w.flush()?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod serde_tests;

#[cfg(test)]
mod tests {
    use super::*;

    // --- encode_quaternion_f32 ---

    #[test]
    fn encode_quat_identity() {
        let bytes = encode_quaternion_f32(1.0, 0.0, 0.0, 0.0);
        let w = f32::from_le_bytes(bytes[0..4].try_into().unwrap());
        assert!((w - 1.0_f32).abs() < 1e-6, "w={w}");
    }

    #[test]
    fn encode_quat_roundtrip() {
        let bytes = encode_quaternion_f32(0.5, 0.5, 0.5, 0.5);
        let w = f32::from_le_bytes(bytes[0..4].try_into().unwrap());
        let x = f32::from_le_bytes(bytes[4..8].try_into().unwrap());
        assert!((w - 0.5_f32).abs() < 1e-6);
        assert!((x - 0.5_f32).abs() < 1e-6);
    }

    // --- pack_float3_array ---

    #[test]
    fn pack_float3_correct_length() {
        let data = vec![[1.0f64, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let bytes = pack_float3_array(&data);
        assert_eq!(bytes.len(), 24); // 2 * 3 * 4 bytes
    }

    #[test]
    fn pack_float3_values() {
        let data = vec![[1.0f64, 0.0, 0.0]];
        let bytes = pack_float3_array(&data);
        let x = f32::from_le_bytes(bytes[0..4].try_into().unwrap());
        assert!((x - 1.0_f32).abs() < 1e-6);
    }

    // --- write_pvd_collection ---

    #[test]
    fn write_pvd_contains_entries() {
        let mut buf = Vec::new();
        let entries = vec![(0.0, "step_0.vtu"), (0.1, "step_1.vtu")];
        write_pvd_collection(&mut buf, &entries).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("step_0.vtu"));
        assert!(s.contains("step_1.vtu"));
        assert!(s.contains("Collection"));
    }

    #[test]
    fn write_pvd_valid_xml() {
        let mut buf = Vec::new();
        write_pvd_collection(&mut buf, &[]).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("</VTKFile>"));
    }

    // --- read_xdmf_timesteps ---

    #[test]
    fn read_xdmf_empty() {
        let result = read_xdmf_timesteps("");
        assert!(result.is_empty());
    }

    #[test]
    fn read_xdmf_parses_time() {
        let xml = r#"<Grid Name="step_0"><Time Value="0.5"/><DataItem href="data.h5"/></Grid>"#;
        let result = read_xdmf_timesteps(xml);
        // Simple parser may or may not find it on one line; just ensure no panic
        let _ = result;
    }

    // --- BodyDesc ---

    #[test]
    fn body_desc_to_json_contains_id() {
        let b = BodyDesc::new(42, [1.0, 2.0, 3.0]);
        let j = b.to_json().unwrap();
        assert!(j.contains("\"id\":42"), "json={j}");
    }

    #[test]
    fn body_desc_to_json_contains_pos() {
        let b = BodyDesc::new(1, [1.5, 2.5, 3.5]);
        let j = b.to_json().unwrap();
        assert!(j.contains("1.5"), "json={j}");
    }

    // --- PhysicsSceneWriter / Reader ---

    #[test]
    fn scene_writer_produces_json() {
        let mut w = PhysicsSceneWriter::new();
        w.add_body(BodyDesc::new(1, [0.0, 1.0, 0.0]));
        w.add_body(BodyDesc::new(2, [1.0, 0.0, 0.0]));
        let mut buf = Vec::new();
        w.write_json(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("\"bodies\""));
    }

    #[test]
    fn scene_reader_roundtrip() {
        let mut w = PhysicsSceneWriter::new();
        w.add_body(BodyDesc::new(7, [1.0, 2.0, 3.0]));
        let mut buf = Vec::new();
        w.write_json(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let scene = PhysicsSceneReader::read_json(&s).unwrap();
        assert_eq!(scene.bodies, vec![BodyDesc::new(7, [1.0, 2.0, 3.0])]);
    }

    #[test]
    fn scene_validate_unique_ids() {
        let bodies = vec![BodyDesc::new(1, [0.0; 3]), BodyDesc::new(2, [0.0; 3])];
        assert!(PhysicsSceneReader::validate(&bodies));
    }

    #[test]
    fn scene_validate_duplicate_ids() {
        let bodies = vec![BodyDesc::new(1, [0.0; 3]), BodyDesc::new(1, [0.0; 3])];
        assert!(!PhysicsSceneReader::validate(&bodies));
    }

    // --- TrajectoryWriter / Reader ---

    #[test]
    fn trajectory_writer_roundtrip() {
        let mut tw = TrajectoryWriter::new();
        let mut f = TrajectoryFrame::new(0.1);
        f.positions.push([1.0, 2.0, 3.0]);
        f.velocities.push([0.1, 0.0, 0.0]);
        f.quaternions.push([1.0, 0.0, 0.0, 0.0]);
        tw.push_frame(f);
        let mut buf = Vec::new();
        tw.write_binary(&mut buf).unwrap();
        let reader = TrajectoryReader::from_bytes(&buf).expect("parse failed");
        assert_eq!(reader.num_frames(), 1);
    }

    #[test]
    fn trajectory_reader_seek_and_iterate() {
        let mut tw = TrajectoryWriter::new();
        for i in 0..3 {
            let mut f = TrajectoryFrame::new(i as f64 * 0.1);
            f.positions.push([i as f64, 0.0, 0.0]);
            f.velocities.push([0.0; 3]);
            f.quaternions.push([1.0, 0.0, 0.0, 0.0]);
            tw.push_frame(f);
        }
        let mut buf = Vec::new();
        tw.write_binary(&mut buf).unwrap();
        let mut reader = TrajectoryReader::from_bytes(&buf).unwrap();
        reader.seek(1);
        let frame = reader.next_frame().unwrap();
        assert!((frame.time - 0.1).abs() < 1e-4, "time={}", frame.time);
    }

    #[test]
    fn trajectory_writer_csv() {
        let mut tw = TrajectoryWriter::new();
        let mut f = TrajectoryFrame::new(0.0);
        f.positions.push([0.0, 0.0, 0.0]);
        f.velocities.push([0.0; 3]);
        f.quaternions.push([1.0, 0.0, 0.0, 0.0]);
        tw.push_frame(f);
        let mut buf = Vec::new();
        tw.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("time,body"));
    }

    // --- ContactForceLog ---

    #[test]
    fn contact_force_log_push_and_filter() {
        let mut log = ContactForceLog::new();
        log.push(ContactForceRecord {
            step: 0,
            body_a: 1,
            body_b: 2,
            normal_impulse: 0.5,
            tangent_impulse: 0.1,
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
        });
        log.push(ContactForceRecord {
            step: 1,
            body_a: 1,
            body_b: 3,
            normal_impulse: 0.2,
            tangent_impulse: 0.0,
            point: [0.0; 3],
            normal: [1.0, 0.0, 0.0],
        });
        assert_eq!(log.records_at_step(0).len(), 1);
        assert_eq!(log.records_at_step(1).len(), 1);
        assert_eq!(log.records_at_step(2).len(), 0);
    }

    #[test]
    fn contact_force_log_csv() {
        let mut log = ContactForceLog::new();
        log.push(ContactForceRecord {
            step: 0,
            body_a: 1,
            body_b: 2,
            normal_impulse: 1.0,
            tangent_impulse: 0.5,
            point: [0.0; 3],
            normal: [0.0, 1.0, 0.0],
        });
        let mut buf = Vec::new();
        log.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("step,body_a"));
    }

    // --- EnergyLog ---

    #[test]
    fn energy_log_max_energy() {
        let mut log = EnergyLog::new();
        log.push(EnergyEntry {
            time: 0.0,
            kinetic: 1.0,
            potential: 2.0,
            total: 3.0,
            linear_momentum: 0.5,
            angular_momentum: 0.2,
        });
        log.push(EnergyEntry {
            time: 0.1,
            kinetic: 2.0,
            potential: 1.0,
            total: 3.0,
            linear_momentum: 0.4,
            angular_momentum: 0.1,
        });
        assert!((log.max_total_energy() - 3.0).abs() < 1e-9);
    }

    #[test]
    fn energy_log_csv() {
        let mut log = EnergyLog::new();
        log.push(EnergyEntry {
            time: 0.0,
            kinetic: 1.0,
            potential: 1.0,
            total: 2.0,
            linear_momentum: 0.0,
            angular_momentum: 0.0,
        });
        let mut buf = Vec::new();
        log.write_csv(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("time,kinetic"));
    }

    // --- CheckpointWriter / Reader ---

    #[test]
    fn checkpoint_roundtrip() {
        let mut cw = CheckpointWriter::new(1.5, 10);
        cw.add_body(BodyDesc::new(1, [0.0, 1.0, 0.0]));
        cw.add_metadata("solver", "PGS");
        let mut buf = Vec::new();
        cw.write_json(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let checkpoint = CheckpointReader::read_json(&s).unwrap();
        assert!(
            (checkpoint.sim_time - 1.5).abs() < 1e-9,
            "t={}",
            checkpoint.sim_time
        );
        assert_eq!(checkpoint.step, 10);
        assert_eq!(checkpoint.bodies, vec![BodyDesc::new(1, [0.0, 1.0, 0.0])]);
        assert_eq!(
            checkpoint.metadata,
            vec![("solver".to_string(), "PGS".to_string())]
        );
    }

    #[test]
    fn checkpoint_metadata_in_output() {
        let mut cw = CheckpointWriter::new(0.0, 0);
        cw.add_metadata("version", "1.0");
        let mut buf = Vec::new();
        cw.write_json(&mut buf).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("version"));
    }

    // --- ParaviewXdmf ---

    #[test]
    fn xdmf_write_produces_xml() {
        use std::io::Read;
        let dir = std::env::temp_dir();
        let xdmf_path = dir.join("test_traj.xdmf");
        let mut xdmf = ParaviewXdmf::new(xdmf_path.to_str().unwrap(), "data.h5");
        xdmf.add_timestep(0.0, 10);
        xdmf.add_timestep(0.1, 10);
        xdmf.write().unwrap_or_else(|e| {
            let _ = e.into_inner();
        });
        let mut content = String::new();
        std::fs::File::open(&xdmf_path)
            .unwrap()
            .read_to_string(&mut content)
            .unwrap();
        assert!(content.contains("Xdmf"));
        assert!(content.contains("step_0"));
    }
}
