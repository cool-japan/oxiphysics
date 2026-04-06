//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions_2::*;
use crate::Error;
use crate::types::{PyContactResult, PyVec3};
use serde::{Deserialize, Serialize};

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{PICKLE_MAGIC, PICKLE_VERSION};

use std::collections::HashMap;

/// Serializable snapshot of the entire world state (legacy format).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WorldState {
    /// Gravity vector.
    pub gravity: PyVec3,
    /// Simulation time.
    pub time: f64,
    /// Positions of all bodies.
    pub positions: Vec<PyVec3>,
    /// Number of bodies.
    pub num_bodies: usize,
}
/// Configuration for incremental state export.
#[derive(Debug, Clone)]
pub struct IncrementalExportConfig {
    /// Only export bodies whose speed exceeds this threshold (m/s).
    pub min_speed_threshold: f64,
    /// Maximum number of bodies per export batch.
    pub max_batch_size: usize,
    /// Whether to include sleeping bodies in the export.
    pub include_sleeping: bool,
}
/// Serializable state of a single rigid body within a snapshot.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct SimBodyState {
    /// Slot handle (u32 index in the world).
    pub handle: u32,
    /// Position `[x, y, z]` in world space.
    pub position: [f64; 3],
    /// Linear velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub orientation: [f64; 4],
    /// Angular velocity `[wx, wy, wz]`.
    pub angular_velocity: [f64; 3],
    /// Whether this body is currently sleeping.
    pub is_sleeping: bool,
    /// Whether this body is static.
    pub is_static: bool,
    /// Optional user tag.
    pub tag: Option<String>,
}
impl SimBodyState {
    /// Create a `SimBodyState` with default motion (at rest, identity orientation).
    pub fn at_rest(handle: u32, position: [f64; 3]) -> Self {
        Self {
            handle,
            position,
            velocity: [0.0; 3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            is_sleeping: false,
            is_static: false,
            tag: None,
        }
    }
    /// Speed (magnitude of linear velocity).
    pub fn speed(&self) -> f64 {
        let v = &self.velocity;
        (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
    }
    /// Angular speed (magnitude of angular velocity).
    pub fn angular_speed(&self) -> f64 {
        let w = &self.angular_velocity;
        (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt()
    }
    /// Kinetic energy estimate using a unit-mass approximation.
    pub fn kinetic_energy_proxy(&self) -> f64 {
        let v = self.speed();
        0.5 * v * v
    }
    /// Distance from origin.
    pub fn distance_from_origin(&self) -> f64 {
        let p = &self.position;
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt()
    }
    /// Whether the body is effectively at rest (speed < threshold).
    pub fn is_at_rest(&self, linear_threshold: f64, angular_threshold: f64) -> bool {
        self.speed() < linear_threshold && self.angular_speed() < angular_threshold
    }
}
/// A complete snapshot of simulation state at a particular moment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    /// Snapshot format version (currently 1).
    pub version: u32,
    /// Simulation time when the snapshot was taken.
    pub time: f64,
    /// Gravity vector `[gx, gy, gz]` at snapshot time.
    pub gravity: [f64; 3],
    /// State of every live body at snapshot time.
    pub bodies: Vec<SimBodyState>,
    /// Contacts active at snapshot time (informational).
    pub contacts: Vec<PyContactResult>,
    /// Number of sleeping bodies at snapshot time.
    pub sleeping_count: usize,
    /// Optional human-readable description.
    pub description: Option<String>,
    /// Arbitrary key-value metadata.
    pub metadata: std::collections::HashMap<String, String>,
}
impl SimulationSnapshot {
    /// Current snapshot format version.
    pub const FORMAT_VERSION: u32 = 1;
    /// Create an empty snapshot at time zero.
    pub fn empty() -> Self {
        Self {
            version: Self::FORMAT_VERSION,
            time: 0.0,
            gravity: [0.0, -9.81, 0.0],
            bodies: Vec::new(),
            contacts: Vec::new(),
            sleeping_count: 0,
            description: None,
            metadata: std::collections::HashMap::new(),
        }
    }
    /// Number of active bodies in this snapshot.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
    /// Number of sleeping bodies in this snapshot.
    pub fn sleeping_count(&self) -> usize {
        self.sleeping_count
    }
    /// Find a body by its handle.
    pub fn find_body(&self, handle: u32) -> Option<&SimBodyState> {
        self.bodies.iter().find(|b| b.handle == handle)
    }
    /// Find a body by its tag.
    pub fn find_by_tag(&self, tag: &str) -> Option<&SimBodyState> {
        self.bodies.iter().find(|b| b.tag.as_deref() == Some(tag))
    }
    /// Total kinetic-energy proxy across all non-sleeping bodies.
    pub fn total_kinetic_energy_proxy(&self) -> f64 {
        self.bodies
            .iter()
            .filter(|b| !b.is_sleeping)
            .map(|b| b.kinetic_energy_proxy())
            .sum()
    }
    /// Return the snapshot as a pretty-printed JSON string.
    pub fn to_pretty_json(&self) -> String {
        serde_json::to_string_pretty(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Serialize to compact JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from JSON.
    pub fn from_json(json: &str) -> Result<Self, crate::Error> {
        serde_json::from_str(json)
            .map_err(|e| crate::Error::General(format!("snapshot deserialization failed: {e}")))
    }
    /// Add a metadata key-value pair.
    pub fn with_metadata(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.metadata.insert(key.into(), value.into());
        self
    }
    /// Set a human-readable description.
    pub fn with_description(mut self, desc: impl Into<String>) -> Self {
        self.description = Some(desc.into());
        self
    }
    /// Serialize to MessagePack bytes.
    pub fn to_msgpack(&self) -> Vec<u8> {
        let json = self.to_json();
        let json_bytes = json.as_bytes();
        let mut result = Vec::with_capacity(8 + json_bytes.len());
        result.extend_from_slice(b"OXIP");
        result.extend_from_slice(&(json_bytes.len() as u32).to_le_bytes());
        result.extend_from_slice(json_bytes);
        result
    }
    /// Deserialize from MessagePack bytes.
    pub fn from_msgpack(data: &[u8]) -> Result<Self, Error> {
        if data.len() < 8 {
            return Err(Error::General("msgpack data too short".to_string()));
        }
        if &data[0..4] != b"OXIP" {
            return Err(Error::General("invalid msgpack magic bytes".to_string()));
        }
        let len = u32::from_le_bytes([data[4], data[5], data[6], data[7]]) as usize;
        if data.len() < 8 + len {
            return Err(Error::General("msgpack data truncated".to_string()));
        }
        let json = std::str::from_utf8(&data[8..8 + len])
            .map_err(|e| Error::General(format!("invalid UTF-8: {e}")))?;
        Self::from_json(json)
    }
    /// Count of static bodies.
    pub fn static_body_count(&self) -> usize {
        self.bodies.iter().filter(|b| b.is_static).count()
    }
    /// Count of dynamic (non-static) bodies.
    pub fn dynamic_body_count(&self) -> usize {
        self.bodies.iter().filter(|b| !b.is_static).count()
    }
    /// All body handles.
    pub fn handles(&self) -> Vec<u32> {
        self.bodies.iter().map(|b| b.handle).collect()
    }
    /// Filter bodies by tag prefix.
    pub fn find_by_tag_prefix(&self, prefix: &str) -> Vec<&SimBodyState> {
        self.bodies
            .iter()
            .filter(|b| b.tag.as_deref().is_some_and(|t| t.starts_with(prefix)))
            .collect()
    }
}
impl SimulationSnapshot {
    /// Compute the diff between `self` (snapshot A) and `other` (snapshot B).
    ///
    /// A body is considered "moved" if its position changed by more than
    /// `position_threshold`.
    pub fn diff(&self, other: &SimulationSnapshot, position_threshold: f64) -> SnapshotDiff {
        let a_map: HashMap<u32, &SimBodyState> =
            self.bodies.iter().map(|b| (b.handle, b)).collect();
        let b_map: HashMap<u32, &SimBodyState> =
            other.bodies.iter().map(|b| (b.handle, b)).collect();
        let removed: Vec<u32> = a_map
            .keys()
            .filter(|h| !b_map.contains_key(h))
            .copied()
            .collect();
        let added: Vec<u32> = b_map
            .keys()
            .filter(|h| !a_map.contains_key(h))
            .copied()
            .collect();
        let mut moved = Vec::new();
        let mut max_disp = 0.0_f64;
        for (handle, a_body) in &a_map {
            if let Some(b_body) = b_map.get(handle) {
                let dx = b_body.position[0] - a_body.position[0];
                let dy = b_body.position[1] - a_body.position[1];
                let dz = b_body.position[2] - a_body.position[2];
                let disp = (dx * dx + dy * dy + dz * dz).sqrt();
                if disp > max_disp {
                    max_disp = disp;
                }
                if disp > position_threshold {
                    moved.push(*handle);
                }
            }
        }
        SnapshotDiff {
            removed,
            added,
            moved,
            max_displacement: max_disp,
            time_delta: other.time - self.time,
        }
    }
}
/// Python-dict-compatible representation of a single body state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BodyDict {
    /// Body handle.
    pub handle: u32,
    /// Position as list.
    pub pos: [f64; 3],
    /// Velocity as list.
    pub vel: [f64; 3],
    /// Orientation quaternion \[x, y, z, w\].
    pub quat: [f64; 4],
    /// Angular velocity.
    pub omega: [f64; 3],
    /// Whether sleeping.
    pub sleeping: bool,
    /// Whether static.
    pub static_body: bool,
    /// Optional tag.
    pub tag: Option<String>,
}
impl BodyDict {
    /// Convert from `SimBodyState`.
    #[allow(dead_code)]
    pub fn from_sim_body(b: &SimBodyState) -> Self {
        Self {
            handle: b.handle,
            pos: b.position,
            vel: b.velocity,
            quat: b.orientation,
            omega: b.angular_velocity,
            sleeping: b.is_sleeping,
            static_body: b.is_static,
            tag: b.tag.clone(),
        }
    }
    /// Convert back to `SimBodyState`.
    #[allow(dead_code)]
    pub fn to_sim_body(&self) -> SimBodyState {
        SimBodyState {
            handle: self.handle,
            position: self.pos,
            velocity: self.vel,
            orientation: self.quat,
            angular_velocity: self.omega,
            is_sleeping: self.sleeping,
            is_static: self.static_body,
            tag: self.tag.clone(),
        }
    }
}
/// A Python-pickle-compatible binary envelope for a `SimulationSnapshot`.
///
/// The real `pickle` protocol uses a complex op-code stream. Here we implement
/// a self-describing binary envelope that can be round-tripped without Python:
///
/// ```text
/// [4]  magic  "OXPK"
/// [1]  version  (currently 2)
/// [4]  payload_len  (u32 LE)
/// [payload_len]  JSON-encoded snapshot
/// [1]  terminal  0x2e  ('.')
/// ```
#[derive(Debug, Clone)]
pub struct PickleEnvelope {
    /// The snapshot stored in this envelope.
    pub snapshot: SimulationSnapshot,
}
impl PickleEnvelope {
    /// Wrap a snapshot in a pickle envelope.
    #[allow(dead_code)]
    pub fn new(snapshot: SimulationSnapshot) -> Self {
        Self { snapshot }
    }
    /// Serialize to bytes.
    #[allow(dead_code)]
    pub fn to_bytes(&self) -> Vec<u8> {
        let json = self.snapshot.to_json();
        let payload = json.as_bytes();
        let mut buf = Vec::with_capacity(10 + payload.len());
        buf.extend_from_slice(PICKLE_MAGIC);
        buf.push(PICKLE_VERSION);
        buf.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        buf.extend_from_slice(payload);
        buf.push(b'.');
        buf
    }
    /// Deserialize from bytes.
    #[allow(dead_code)]
    pub fn from_bytes(data: &[u8]) -> Result<Self, Error> {
        if data.len() < 10 {
            return Err(Error::General("pickle envelope too short".to_string()));
        }
        if &data[0..4] != PICKLE_MAGIC {
            return Err(Error::General("invalid pickle magic bytes".to_string()));
        }
        let _version = data[4];
        let payload_len = u32::from_le_bytes([data[5], data[6], data[7], data[8]]) as usize;
        if data.len() < 9 + payload_len + 1 {
            return Err(Error::General("pickle envelope truncated".to_string()));
        }
        let json = std::str::from_utf8(&data[9..9 + payload_len])
            .map_err(|e| Error::General(format!("invalid UTF-8 in pickle: {e}")))?;
        let snapshot = SimulationSnapshot::from_json(json)?;
        Ok(Self { snapshot })
    }
    /// Serialize to a hex string (for embedding in Python source).
    #[allow(dead_code)]
    pub fn to_hex(&self) -> String {
        self.to_bytes()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect::<Vec<_>>()
            .join("")
    }
}
/// Result of JSON schema validation.
#[derive(Debug, Clone)]
pub struct SchemaValidationResult {
    /// Whether the JSON conforms to the expected schema.
    pub is_valid: bool,
    /// List of schema violation descriptions.
    pub errors: Vec<String>,
}
impl SchemaValidationResult {
    /// Create a passing result.
    #[allow(dead_code)]
    pub fn ok() -> Self {
        Self {
            is_valid: true,
            errors: Vec::new(),
        }
    }
    /// Create a failing result with a single error.
    #[allow(dead_code)]
    pub fn err(msg: impl Into<String>) -> Self {
        Self {
            is_valid: false,
            errors: vec![msg.into()],
        }
    }
}
/// A NumPy-compatible flat array descriptor for body positions.
///
/// Provides the data, shape, and dtype needed to reconstruct a NumPy array on
/// the Python side via `np.frombuffer` or `np.array`.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumpyPositionArray {
    /// Flat f64 data, row-major: `[x0, y0, z0, x1, y1, z1, ...]`.
    pub data: Vec<f64>,
    /// Shape: `[n_bodies, 3]`.
    pub shape: [usize; 2],
    /// NumPy dtype string.
    pub dtype: String,
    /// Whether data is in C (row-major) order.
    pub c_order: bool,
}
impl NumpyPositionArray {
    /// Build a position array from a snapshot.
    #[allow(dead_code)]
    pub fn from_snapshot(snap: &SimulationSnapshot) -> Self {
        let n = snap.bodies.len();
        let mut data = Vec::with_capacity(n * 3);
        for b in &snap.bodies {
            data.extend_from_slice(&b.position);
        }
        Self {
            data,
            shape: [n, 3],
            dtype: "float64".to_string(),
            c_order: true,
        }
    }
    /// Build a velocity array from a snapshot.
    #[allow(dead_code)]
    pub fn velocity_array(snap: &SimulationSnapshot) -> Self {
        let n = snap.bodies.len();
        let mut data = Vec::with_capacity(n * 3);
        for b in &snap.bodies {
            data.extend_from_slice(&b.velocity);
        }
        Self {
            data,
            shape: [n, 3],
            dtype: "float64".to_string(),
            c_order: true,
        }
    }
    /// Get the position of body at row `i` as `[x, y, z]`.
    #[allow(dead_code)]
    pub fn get_row(&self, i: usize) -> Option<[f64; 3]> {
        if i >= self.shape[0] {
            return None;
        }
        let base = i * 3;
        Some([self.data[base], self.data[base + 1], self.data[base + 2]])
    }
    /// Number of rows (bodies).
    #[allow(dead_code)]
    pub fn n_rows(&self) -> usize {
        self.shape[0]
    }
    /// Total number of f64 elements.
    #[allow(dead_code)]
    pub fn size(&self) -> usize {
        self.data.len()
    }
    /// Serialize to JSON.
    #[allow(dead_code)]
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Serialize data to raw bytes (f64 LE).
    #[allow(dead_code)]
    pub fn to_raw_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(self.data.len() * 8);
        for &v in &self.data {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }
}
/// A serialized representation of a single body's state as a JSON-ready struct.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct BodyStateJson {
    /// Body handle.
    pub handle: u32,
    /// World-space position `[x, y, z]`.
    pub position: [f64; 3],
    /// Linear velocity `[vx, vy, vz]`.
    pub velocity: [f64; 3],
    /// Orientation quaternion `[x, y, z, w]`.
    pub orientation: [f64; 4],
    /// Angular velocity `[wx, wy, wz]`.
    pub angular_velocity: [f64; 3],
    /// Whether body is sleeping.
    pub is_sleeping: bool,
    /// Whether body is static.
    pub is_static: bool,
    /// Optional user tag.
    pub tag: Option<String>,
    /// Schema version string.
    pub schema_version: String,
}
impl BodyStateJson {
    /// Create a `BodyStateJson` from a `SimBodyState`.
    pub fn from_sim_body(body: &SimBodyState) -> Self {
        Self {
            handle: body.handle,
            position: body.position,
            velocity: body.velocity,
            orientation: body.orientation,
            angular_velocity: body.angular_velocity,
            is_sleeping: body.is_sleeping,
            is_static: body.is_static,
            tag: body.tag.clone(),
            schema_version: "1.0.0".to_string(),
        }
    }
    /// Convert back to a `SimBodyState`.
    pub fn to_sim_body(&self) -> SimBodyState {
        SimBodyState {
            handle: self.handle,
            position: self.position,
            velocity: self.velocity,
            orientation: self.orientation,
            angular_velocity: self.angular_velocity,
            is_sleeping: self.is_sleeping,
            is_static: self.is_static,
            tag: self.tag.clone(),
        }
    }
}
/// A simulation checkpoint with metadata for save/restore workflows.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationCheckpoint {
    /// Checkpoint format version.
    pub version: u32,
    /// Human-readable checkpoint label.
    pub label: String,
    /// Wall-clock timestamp (Unix seconds, approximate).
    pub timestamp: f64,
    /// Simulation time at checkpoint.
    pub sim_time: f64,
    /// Number of steps taken.
    pub step_count: u64,
    /// Gravity vector at checkpoint.
    pub gravity: [f64; 3],
    /// All body states.
    pub bodies: Vec<BodyStateJson>,
    /// Arbitrary metadata key-value pairs.
    pub metadata: std::collections::HashMap<String, String>,
}
impl SimulationCheckpoint {
    /// Current checkpoint format version.
    pub const FORMAT_VERSION: u32 = 1;
    /// Create an empty checkpoint.
    pub fn empty(label: impl Into<String>) -> Self {
        Self {
            version: Self::FORMAT_VERSION,
            label: label.into(),
            timestamp: 0.0,
            sim_time: 0.0,
            step_count: 0,
            gravity: [0.0, -9.81, 0.0],
            bodies: Vec::new(),
            metadata: std::collections::HashMap::new(),
        }
    }
    /// Build a checkpoint from a `SimulationSnapshot`.
    pub fn from_snapshot(
        snap: &SimulationSnapshot,
        label: impl Into<String>,
        step_count: u64,
        timestamp: f64,
    ) -> Self {
        let bodies = snap
            .bodies
            .iter()
            .map(BodyStateJson::from_sim_body)
            .collect();
        Self {
            version: Self::FORMAT_VERSION,
            label: label.into(),
            timestamp,
            sim_time: snap.time,
            step_count,
            gravity: snap.gravity,
            bodies,
            metadata: snap.metadata.clone(),
        }
    }
    /// Serialize to JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from JSON.
    pub fn from_json(json: &str) -> Result<Self, crate::Error> {
        serde_json::from_str(json)
            .map_err(|e| crate::Error::General(format!("checkpoint deserialization failed: {e}")))
    }
    /// Convert back to a `SimulationSnapshot`.
    pub fn to_snapshot(&self) -> SimulationSnapshot {
        let bodies = self.bodies.iter().map(|b| b.to_sim_body()).collect();
        SimulationSnapshot {
            version: SimulationSnapshot::FORMAT_VERSION,
            time: self.sim_time,
            gravity: self.gravity,
            bodies,
            contacts: Vec::new(),
            sleeping_count: 0,
            description: Some(format!("Restored from checkpoint '{}'", self.label)),
            metadata: self.metadata.clone(),
        }
    }
    /// Number of bodies in the checkpoint.
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }
}
/// An incremental state update (delta) for efficient streaming.
///
/// Instead of sending the full snapshot every frame, only changed bodies
/// are included in the delta.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IncrementalUpdate {
    /// Sequence number for ordering.
    pub sequence: u64,
    /// Timestamp of this update.
    pub time: f64,
    /// Only bodies that changed since the last update.
    pub changed_bodies: Vec<SimBodyState>,
    /// Handles of bodies that were removed.
    pub removed_handles: Vec<u32>,
    /// Handles of bodies that were added.
    pub added_handles: Vec<u32>,
}
impl IncrementalUpdate {
    /// Create an empty incremental update.
    pub fn empty(sequence: u64, time: f64) -> Self {
        Self {
            sequence,
            time,
            changed_bodies: Vec::new(),
            removed_handles: Vec::new(),
            added_handles: Vec::new(),
        }
    }
    /// Serialize to JSON.
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from JSON.
    pub fn from_json(json: &str) -> Result<Self, Error> {
        serde_json::from_str(json)
            .map_err(|e| Error::General(format!("incremental update deserialization: {e}")))
    }
    /// Whether this update is empty (no changes).
    pub fn is_empty(&self) -> bool {
        self.changed_bodies.is_empty()
            && self.removed_handles.is_empty()
            && self.added_handles.is_empty()
    }
    /// Number of changes in this update.
    pub fn change_count(&self) -> usize {
        self.changed_bodies.len() + self.removed_handles.len() + self.added_handles.len()
    }
}
/// A single incremental export batch.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExportBatch {
    /// Batch index (0-based).
    pub batch_index: usize,
    /// Whether this is the final batch.
    pub is_last: bool,
    /// Total number of batches.
    pub total_batches: usize,
    /// Simulation time.
    pub time: f64,
    /// Body states in this batch.
    pub bodies: Vec<SimBodyState>,
}
impl ExportBatch {
    /// Serialize to JSON.
    #[allow(dead_code)]
    pub fn to_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from JSON.
    #[allow(dead_code)]
    pub fn from_json(json: &str) -> Result<Self, Error> {
        serde_json::from_str(json)
            .map_err(|e| Error::General(format!("ExportBatch deserialization: {e}")))
    }
}
/// Schema version information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SchemaVersion {
    /// Major version (breaking changes).
    pub major: u32,
    /// Minor version (backward-compatible additions).
    pub minor: u32,
    /// Patch version.
    pub patch: u32,
}
impl SchemaVersion {
    /// Current schema version.
    pub fn current() -> Self {
        Self {
            major: 1,
            minor: 0,
            patch: 0,
        }
    }
    /// Check if this version is compatible with another.
    pub fn is_compatible_with(&self, other: &SchemaVersion) -> bool {
        self.major == other.major
    }
    /// Version string (e.g. "1.0.0").
    pub fn to_string_version(&self) -> String {
        format!("{}.{}.{}", self.major, self.minor, self.patch)
    }
}
/// A Python-dict-compatible representation of a `SimulationSnapshot`.
///
/// The `to_dict_json` method returns a JSON string structured like a Python dict
/// that the Python side can `eval()` or `json.loads()` directly.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnapshotDict {
    /// Snapshot version.
    pub version: u32,
    /// Simulation time.
    pub time: f64,
    /// Gravity as list.
    pub gravity: [f64; 3],
    /// Body states as list of dicts.
    pub bodies: Vec<BodyDict>,
    /// Number of contacts.
    pub n_contacts: usize,
    /// Optional description.
    pub description: Option<String>,
}
impl SnapshotDict {
    /// Build from a `SimulationSnapshot`.
    #[allow(dead_code)]
    pub fn from_snapshot(snap: &SimulationSnapshot) -> Self {
        Self {
            version: snap.version,
            time: snap.time,
            gravity: snap.gravity,
            bodies: snap.bodies.iter().map(BodyDict::from_sim_body).collect(),
            n_contacts: snap.contacts.len(),
            description: snap.description.clone(),
        }
    }
    /// Convert back to `SimulationSnapshot`.
    #[allow(dead_code)]
    pub fn to_snapshot(&self) -> SimulationSnapshot {
        let bodies: Vec<SimBodyState> = self.bodies.iter().map(|b| b.to_sim_body()).collect();
        let sleeping_count = bodies.iter().filter(|b| b.is_sleeping).count();
        SimulationSnapshot {
            version: self.version,
            time: self.time,
            gravity: self.gravity,
            bodies,
            contacts: Vec::new(),
            sleeping_count,
            description: self.description.clone(),
            metadata: std::collections::HashMap::new(),
        }
    }
    /// Serialize to JSON string (like Python's `json.dumps`).
    #[allow(dead_code)]
    pub fn to_dict_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }
    /// Deserialize from JSON string.
    #[allow(dead_code)]
    pub fn from_dict_json(json: &str) -> Result<Self, Error> {
        serde_json::from_str(json)
            .map_err(|e| Error::General(format!("SnapshotDict deserialization: {e}")))
    }
}
/// Difference report between two snapshots.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SnapshotDiff {
    /// Handles present in A but not in B (removed bodies).
    pub removed: Vec<u32>,
    /// Handles present in B but not in A (added bodies).
    pub added: Vec<u32>,
    /// Handles present in both with changed position (beyond threshold).
    pub moved: Vec<u32>,
    /// Maximum position displacement among all common bodies.
    pub max_displacement: f64,
    /// Time difference between snapshots.
    pub time_delta: f64,
}
impl SnapshotDiff {
    /// Returns `true` if there are no differences.
    pub fn is_identical(&self) -> bool {
        self.removed.is_empty() && self.added.is_empty() && self.moved.is_empty()
    }
    /// Total number of changed bodies.
    pub fn change_count(&self) -> usize {
        self.removed.len() + self.added.len() + self.moved.len()
    }
}
/// Validation result for a deserialized snapshot.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ValidationResult {
    /// Whether the snapshot is valid.
    pub is_valid: bool,
    /// List of validation warnings/errors.
    pub issues: Vec<String>,
}
