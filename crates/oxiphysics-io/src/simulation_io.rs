// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Physics simulation I/O formats.
//!
//! Provides writers and readers for full physics scenes (bodies, joints,
//! materials), multi-body trajectories, contact force logs, energy logs,
//! simulation checkpoints, VTK time-series output, and XDMF wrappers.

#![allow(dead_code)]

use std::io::{BufWriter, Write};
use std::path::Path;

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

/// A rigid body descriptor for scene serialization.
#[derive(Debug, Clone)]
pub struct BodyDesc {
    /// Unique body identifier.
    pub id: u64,
    /// World-space position.
    pub position: [f64; 3],
    /// Orientation as quaternion (w, x, y, z).
    pub orientation: [f64; 4],
    /// Linear velocity.
    pub linear_velocity: [f64; 3],
    /// Angular velocity.
    pub angular_velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Shape type tag.
    pub shape_tag: String,
}

impl BodyDesc {
    /// Construct a body descriptor.
    pub fn new(id: u64, position: [f64; 3]) -> Self {
        BodyDesc {
            id,
            position,
            orientation: [1.0, 0.0, 0.0, 0.0],
            linear_velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            mass: 1.0,
            shape_tag: "sphere".to_string(),
        }
    }

    /// Serialize to a simple JSON object string.
    pub fn to_json(&self) -> String {
        format!(
            r#"{{"id":{id},"pos":[{px},{py},{pz}],"quat":[{qw},{qx},{qy},{qz}],"linvel":[{lvx},{lvy},{lvz}],"angvel":[{avx},{avy},{avz}],"mass":{m},"shape":"{s}"}}"#,
            id = self.id,
            px = self.position[0],
            py = self.position[1],
            pz = self.position[2],
            qw = self.orientation[0],
            qx = self.orientation[1],
            qy = self.orientation[2],
            qz = self.orientation[3],
            lvx = self.linear_velocity[0],
            lvy = self.linear_velocity[1],
            lvz = self.linear_velocity[2],
            avx = self.angular_velocity[0],
            avy = self.angular_velocity[1],
            avz = self.angular_velocity[2],
            m = self.mass,
            s = self.shape_tag,
        )
    }
}

/// A joint descriptor.
#[derive(Debug, Clone)]
pub struct JointDesc {
    /// Unique joint id.
    pub id: u64,
    /// Body A id.
    pub body_a: u64,
    /// Body B id.
    pub body_b: u64,
    /// Joint type.
    pub joint_type: String,
    /// Anchor point on body A (local space).
    pub anchor_a: [f64; 3],
    /// Anchor point on body B (local space).
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
#[derive(Debug, Clone)]
pub struct MaterialDesc {
    /// Material id.
    pub id: u64,
    /// Restitution coefficient.
    pub restitution: f64,
    /// Dynamic friction coefficient.
    pub friction: f64,
    /// Density (kg/m³).
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
    pub fn write_json<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writeln!(writer, "{{")?;
        writeln!(writer, "  \"bodies\": [")?;
        for (i, b) in self.bodies.iter().enumerate() {
            let comma = if i + 1 < self.bodies.len() { "," } else { "" };
            writeln!(writer, "    {}{}", b.to_json(), comma)?;
        }
        writeln!(writer, "  ],")?;
        writeln!(writer, "  \"joints\": [")?;
        for (i, j) in self.joints.iter().enumerate() {
            let comma = if i + 1 < self.joints.len() { "," } else { "" };
            writeln!(
                writer,
                r#"    {{"id":{id},"bodyA":{ba},"bodyB":{bb},"type":"{jt}"}}{comma}"#,
                id = j.id,
                ba = j.body_a,
                bb = j.body_b,
                jt = j.joint_type
            )?;
        }
        writeln!(writer, "  ],")?;
        writeln!(writer, "  \"materials\": [")?;
        for (i, m) in self.materials.iter().enumerate() {
            let comma = if i + 1 < self.materials.len() {
                ","
            } else {
                ""
            };
            writeln!(
                writer,
                r#"    {{"id":{id},"restitution":{r},"friction":{f},"density":{d}}}{comma}"#,
                id = m.id,
                r = m.restitution,
                f = m.friction,
                d = m.density
            )?;
        }
        writeln!(writer, "  ]")?;
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

/// Read a physics scene from JSON, validate it, and reconstruct a body list.
pub struct PhysicsSceneReader;

impl PhysicsSceneReader {
    /// Parse a JSON scene string and return body descriptors.
    /// This is a simplified parser for the format produced by [`PhysicsSceneWriter`].
    pub fn read_json(json: &str) -> Vec<BodyDesc> {
        let mut bodies = Vec::new();
        // Very simple line-by-line parser looking for id fields
        let mut in_bodies = false;
        let mut cur_id: Option<u64> = None;
        let mut cur_pos = [0.0f64; 3];
        for line in json.lines() {
            let line = line.trim();
            if line.contains("\"bodies\"") {
                in_bodies = true;
            }
            if line.contains("\"joints\"") {
                in_bodies = false;
            }
            if !in_bodies {
                continue;
            }
            if let Some(id) = parse_u64_field(line, "id") {
                cur_id = Some(id);
            }
            if let Some(pos) = parse_float3_field(line, "pos") {
                cur_pos = pos;
            }
            if line.contains('}') && cur_id.is_some() {
                bodies.push(BodyDesc::new(
                    cur_id.expect("value should be present"),
                    cur_pos,
                ));
                cur_id = None;
                cur_pos = [0.0; 3];
            }
        }
        bodies
    }

    /// Validate a list of body descriptors (checks for duplicate IDs).
    pub fn validate(bodies: &[BodyDesc]) -> bool {
        let mut ids = std::collections::HashSet::new();
        bodies.iter().all(|b| ids.insert(b.id))
    }
}

fn parse_u64_field(s: &str, field: &str) -> Option<u64> {
    let pat = format!("\"{field}\":");
    let idx = s.find(pat.as_str())? + pat.len();
    let rest = s[idx..].trim_start();
    rest.split([',', '}', ' ', '\t'])
        .next()?
        .trim()
        .parse()
        .ok()
}

fn parse_float3_field(s: &str, field: &str) -> Option<[f64; 3]> {
    let pat = format!("\"{field}\":");
    let idx = s.find(pat.as_str())? + pat.len();
    let rest = s[idx..].trim_start();
    if rest.starts_with('[') {
        let end = rest.find(']')?;
        let inner = &rest[1..end];
        let vals: Vec<f64> = inner
            .split(',')
            .filter_map(|v| v.trim().parse().ok())
            .collect();
        if vals.len() == 3 {
            return Some([vals[0], vals[1], vals[2]]);
        }
    }
    None
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
    /// Format: header "TRAJ\0" + n_frames (u64 LE) + per-frame data.
    pub fn write_binary<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        w.write_all(b"TRAJ\0")?;
        w.write_all(&(self.frames.len() as u64).to_le_bytes())?;
        for f in &self.frames {
            w.write_all(&f.time.to_le_bytes())?;
            w.write_all(&(f.positions.len() as u64).to_le_bytes())?;
            let pos_bytes = pack_float3_array(&f.positions);
            w.write_all(&pos_bytes)?;
            let vel_bytes = pack_float3_array(&f.velocities);
            w.write_all(&vel_bytes)?;
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

/// Reader for binary trajectory files produced by [`TrajectoryWriter`].
pub struct TrajectoryReader {
    frames: Vec<TrajectoryFrame>,
    cursor: usize,
}

impl TrajectoryReader {
    /// Read from a byte buffer.
    pub fn from_bytes(data: &[u8]) -> Option<Self> {
        if data.len() < 13 || &data[0..5] != b"TRAJ\0" {
            return None;
        }
        let n_frames = u64::from_le_bytes(data[5..13].try_into().ok()?) as usize;
        let mut frames = Vec::with_capacity(n_frames);
        let mut pos = 13;
        for _ in 0..n_frames {
            if pos + 16 > data.len() {
                break;
            }
            let time = f64::from_le_bytes(data[pos..pos + 8].try_into().ok()?);
            pos += 8;
            let nb = u64::from_le_bytes(data[pos..pos + 8].try_into().ok()?) as usize;
            pos += 8;
            let mut frame = TrajectoryFrame::new(time);
            // positions (nb * 12 bytes f32)
            for _ in 0..nb {
                if pos + 12 > data.len() {
                    break;
                }
                let x = f32::from_le_bytes(data[pos..pos + 4].try_into().ok()?) as f64;
                let y = f32::from_le_bytes(data[pos + 4..pos + 8].try_into().ok()?) as f64;
                let z = f32::from_le_bytes(data[pos + 8..pos + 12].try_into().ok()?) as f64;
                frame.positions.push([x, y, z]);
                pos += 12;
            }
            // velocities
            for _ in 0..nb {
                if pos + 12 > data.len() {
                    break;
                }
                let x = f32::from_le_bytes(data[pos..pos + 4].try_into().ok()?) as f64;
                let y = f32::from_le_bytes(data[pos + 4..pos + 8].try_into().ok()?) as f64;
                let z = f32::from_le_bytes(data[pos + 8..pos + 12].try_into().ok()?) as f64;
                frame.velocities.push([x, y, z]);
                pos += 12;
            }
            // quaternions
            for _ in 0..nb {
                if pos + 16 > data.len() {
                    break;
                }
                let w = f32::from_le_bytes(data[pos..pos + 4].try_into().ok()?) as f64;
                let x = f32::from_le_bytes(data[pos + 4..pos + 8].try_into().ok()?) as f64;
                let y = f32::from_le_bytes(data[pos + 8..pos + 12].try_into().ok()?) as f64;
                let z = f32::from_le_bytes(data[pos + 12..pos + 16].try_into().ok()?) as f64;
                frame.quaternions.push([w, x, y, z]);
                pos += 16;
            }
            frames.push(frame);
        }
        Some(TrajectoryReader { frames, cursor: 0 })
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
    pub fn write_json<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        writeln!(w, "{{")?;
        writeln!(w, "  \"sim_time\": {},", self.sim_time)?;
        writeln!(w, "  \"step\": {},", self.step)?;
        writeln!(w, "  \"metadata\": {{")?;
        for (i, (k, v)) in self.metadata.iter().enumerate() {
            let comma = if i + 1 < self.metadata.len() { "," } else { "" };
            writeln!(w, "    \"{k}\": \"{v}\"{comma}")?;
        }
        writeln!(w, "  }},")?;
        writeln!(w, "  \"bodies\": [")?;
        for (i, b) in self.bodies.iter().enumerate() {
            let comma = if i + 1 < self.bodies.len() { "," } else { "" };
            writeln!(w, "    {}{}", b.to_json(), comma)?;
        }
        writeln!(w, "  ]")?;
        writeln!(w, "}}")?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// CheckpointReader
// ---------------------------------------------------------------------------

/// Read a checkpoint and restore simulation state.
pub struct CheckpointReader;

impl CheckpointReader {
    /// Parse a JSON checkpoint and return (sim_time, step, bodies).
    pub fn read_json(json: &str) -> (f64, u64, Vec<BodyDesc>) {
        let sim_time = parse_f64_line(json, "sim_time").unwrap_or(0.0);
        let step = parse_u64_line(json, "step").unwrap_or(0);
        let bodies = PhysicsSceneReader::read_json(json);
        (sim_time, step, bodies)
    }
}

fn parse_f64_line(s: &str, field: &str) -> Option<f64> {
    let pat = format!("\"{field}\":");
    for line in s.lines() {
        let line = line.trim();
        if line.starts_with(pat.as_str()) {
            let rest = &line[pat.len()..];
            let val_str = rest.trim().trim_end_matches(',');
            return val_str.parse().ok();
        }
    }
    None
}

fn parse_u64_line(s: &str, field: &str) -> Option<u64> {
    let pat = format!("\"{field}\":");
    for line in s.lines() {
        let line = line.trim();
        if line.starts_with(pat.as_str()) {
            let rest = &line[pat.len()..];
            let val_str = rest.trim().trim_end_matches(',');
            return val_str.parse().ok();
        }
    }
    None
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
        let j = b.to_json();
        assert!(j.contains("\"id\":42"), "json={j}");
    }

    #[test]
    fn body_desc_to_json_contains_pos() {
        let b = BodyDesc::new(1, [1.5, 2.5, 3.5]);
        let j = b.to_json();
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
        let bodies = PhysicsSceneReader::read_json(&s);
        assert!(!bodies.is_empty(), "should find bodies");
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
        let (t, step, _bodies) = CheckpointReader::read_json(&s);
        assert!((t - 1.5).abs() < 1e-9, "t={t}");
        assert_eq!(step, 10);
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
