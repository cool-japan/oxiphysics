// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! WebAssembly I/O bridge.
//!
//! Provides types for exporting particle data, mesh geometry, state snapshots,
//! trajectory buffers, and simulation configurations across the WASM boundary.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::collections::VecDeque;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// ExportFormat
// ---------------------------------------------------------------------------

/// Enumeration of supported data-export formats.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExportFormat {
    /// Comma-separated values.
    Csv,
    /// JavaScript Object Notation.
    Json,
    /// Raw binary data.
    Binary,
    /// Extended XYZ format (chemistry/molecular dynamics convention).
    Xyz,
}

// ---------------------------------------------------------------------------
// ParticleRecord (internal helper)
// ---------------------------------------------------------------------------

/// Internal record storing per-particle data.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct ParticleRecord {
    /// Position in world space (m).
    pos: [f64; 3],
    /// Velocity (m/s).
    vel: [f64; 3],
    /// Mass (kg).
    mass: f64,
}

// ---------------------------------------------------------------------------
// ParticleExporter
// ---------------------------------------------------------------------------

/// Collects particle data and exports it in various formats.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParticleExporter {
    /// The format that was chosen at construction time.
    pub format: ExportFormat,
    /// Accumulated particle records.
    particles: Vec<ParticleRecord>,
}

impl ParticleExporter {
    /// Create a new `ParticleExporter` targeting the given `format`.
    pub fn new(format: ExportFormat) -> Self {
        ParticleExporter {
            format,
            particles: Vec::new(),
        }
    }

    /// Append a particle with the given `pos`ition, `vel`ocity, and `mass`.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64) {
        self.particles.push(ParticleRecord { pos, vel, mass });
    }

    /// Export all particles as a CSV string.
    ///
    /// Header: `px,py,pz,vx,vy,vz,mass`
    pub fn export_csv(&self) -> String {
        let mut out = String::from("px,py,pz,vx,vy,vz,mass\n");
        for p in &self.particles {
            out.push_str(&format!(
                "{},{},{},{},{},{},{}\n",
                p.pos[0], p.pos[1], p.pos[2], p.vel[0], p.vel[1], p.vel[2], p.mass
            ));
        }
        out
    }

    /// Export all particles as a JSON string.
    pub fn export_json(&self) -> String {
        serde_json::to_string(&self.particles).unwrap_or_else(|_| "[]".to_string())
    }

    /// Export all particles in XYZ format (count line, comment line, then one
    /// atom per line as `X px py pz`).
    pub fn export_xyz(&self) -> String {
        let mut out = format!("{}\n", self.particles.len());
        out.push_str("OxiPhysics XYZ export\n");
        for p in &self.particles {
            out.push_str(&format!("X {} {} {}\n", p.pos[0], p.pos[1], p.pos[2]));
        }
        out
    }

    /// Number of particles currently stored.
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// Returns `true` if no particles have been added.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }
}

// ---------------------------------------------------------------------------
// MeshExporter
// ---------------------------------------------------------------------------

/// Collects mesh vertices and triangular faces and exports them as OBJ or STL.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshExporter {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangular face indices (zero-based).
    pub faces: Vec<[usize; 3]>,
}

impl MeshExporter {
    /// Create a new empty `MeshExporter`.
    pub fn new() -> Self {
        MeshExporter {
            vertices: Vec::new(),
            faces: Vec::new(),
        }
    }

    /// Append a vertex at position `v`.
    pub fn add_vertex(&mut self, v: [f64; 3]) {
        self.vertices.push(v);
    }

    /// Append a triangular face defined by vertex indices `i` (zero-based).
    pub fn add_face(&mut self, i: [usize; 3]) {
        self.faces.push(i);
    }

    /// Export the mesh as a Wavefront OBJ string.
    ///
    /// Vertex indices in the OBJ file are one-based as per the format
    /// specification.
    pub fn export_obj(&self) -> String {
        let mut out = String::from("# OxiPhysics OBJ export\n");
        for v in &self.vertices {
            out.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
        }
        for f in &self.faces {
            // OBJ indices are 1-based
            out.push_str(&format!("f {} {} {}\n", f[0] + 1, f[1] + 1, f[2] + 1));
        }
        out
    }

    /// Export the mesh as an ASCII STL string.
    pub fn export_stl_ascii(&self) -> String {
        let mut out = String::from("solid oxiphysics\n");
        for f in &self.faces {
            let v0 = self.vertices[f[0]];
            let v1 = self.vertices[f[1]];
            let v2 = self.vertices[f[2]];
            // Compute face normal via cross product
            let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
            let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
            let nx = e1[1] * e2[2] - e1[2] * e2[1];
            let ny = e1[2] * e2[0] - e1[0] * e2[2];
            let nz = e1[0] * e2[1] - e1[1] * e2[0];
            let len = (nx * nx + ny * ny + nz * nz).sqrt().max(1e-30);
            out.push_str(&format!(
                "  facet normal {} {} {}\n    outer loop\n",
                nx / len,
                ny / len,
                nz / len
            ));
            out.push_str(&format!(
                "      vertex {} {} {}\n      vertex {} {} {}\n      vertex {} {} {}\n",
                v0[0], v0[1], v0[2], v1[0], v1[1], v1[2], v2[0], v2[1], v2[2]
            ));
            out.push_str("    endloop\n  endfacet\n");
        }
        out.push_str("endsolid oxiphysics\n");
        out
    }

    /// Number of vertices stored.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of faces stored.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }
}

impl Default for MeshExporter {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// StateSnapshot
// ---------------------------------------------------------------------------

/// A complete serialisable snapshot of simulation state at a single step.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StateSnapshot {
    /// Simulation step index.
    pub step: u64,
    /// Simulation time (seconds).
    pub time: f64,
    /// Number of rigid bodies in the snapshot.
    pub n_bodies: usize,
    /// Body positions `[x, y, z]` for each body.
    pub positions: Vec<[f64; 3]>,
    /// Body velocities `[vx, vy, vz]` for each body.
    pub velocities: Vec<[f64; 3]>,
}

impl StateSnapshot {
    /// Create a new `StateSnapshot`.
    pub fn new(step: u64, time: f64, positions: Vec<[f64; 3]>, velocities: Vec<[f64; 3]>) -> Self {
        let n_bodies = positions.len();
        StateSnapshot {
            step,
            time,
            n_bodies,
            positions,
            velocities,
        }
    }

    /// Serialise this snapshot to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(self).unwrap_or_else(|_| "{}".to_string())
    }

    /// Deserialise a `StateSnapshot` from a JSON string.
    ///
    /// Returns `None` if parsing fails.
    pub fn deserialize_json(s: &str) -> Option<Self> {
        serde_json::from_str(s).ok()
    }
}

// ---------------------------------------------------------------------------
// TrajectoryBuffer
// ---------------------------------------------------------------------------

/// A bounded ring-buffer of `StateSnapshot` frames.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryBuffer {
    /// Maximum number of frames to retain.
    pub max_frames: usize,
    /// Stored frames in insertion order.
    pub frames: VecDeque<StateSnapshot>,
}

impl TrajectoryBuffer {
    /// Create a new `TrajectoryBuffer` with the given capacity.
    pub fn new(max_frames: usize) -> Self {
        TrajectoryBuffer {
            max_frames,
            frames: VecDeque::new(),
        }
    }

    /// Push a new frame, evicting the oldest if the buffer is full.
    pub fn push_frame(&mut self, s: StateSnapshot) {
        if self.frames.len() >= self.max_frames {
            self.frames.pop_front();
        }
        self.frames.push_back(s);
    }

    /// Return a reference to the frame at index `i`, or `None` if out of
    /// bounds.
    pub fn get_frame(&self, i: usize) -> Option<&StateSnapshot> {
        self.frames.get(i)
    }

    /// Export the full trajectory as a JSON string.
    pub fn export_trajectory_json(&self) -> String {
        let frames: Vec<&StateSnapshot> = self.frames.iter().collect();
        serde_json::to_string(&frames).unwrap_or_else(|_| "[]".to_string())
    }

    /// Number of frames currently stored.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Returns `true` if no frames have been stored.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }
}

// ---------------------------------------------------------------------------
// ConfigSerializer
// ---------------------------------------------------------------------------

/// Serialises and deserialises arbitrary simulation configuration maps
/// (key-value pairs of strings) to and from JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConfigSerializer {
    /// Internal key-value store for configuration entries.
    pub entries: std::collections::HashMap<String, String>,
}

impl ConfigSerializer {
    /// Create a new empty `ConfigSerializer`.
    pub fn new() -> Self {
        ConfigSerializer {
            entries: std::collections::HashMap::new(),
        }
    }

    /// Insert or update a configuration key-value pair.
    pub fn set(&mut self, key: &str, value: &str) {
        self.entries.insert(key.to_string(), value.to_string());
    }

    /// Retrieve a configuration value by key, or `None` if absent.
    pub fn get(&self, key: &str) -> Option<&str> {
        self.entries.get(key).map(String::as_str)
    }

    /// Serialise all configuration entries to a JSON string.
    pub fn serialize_json(&self) -> String {
        serde_json::to_string(&self.entries).unwrap_or_else(|_| "{}".to_string())
    }

    /// Deserialise configuration from a JSON string, replacing the current
    /// entries.
    ///
    /// Returns `false` if parsing fails.
    pub fn deserialize_json(&mut self, s: &str) -> bool {
        match serde_json::from_str::<std::collections::HashMap<String, String>>(s) {
            Ok(map) => {
                self.entries = map;
                true
            }
            Err(_) => false,
        }
    }

    /// Number of configuration entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if no entries are stored.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

impl Default for ConfigSerializer {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---- ExportFormat ----

    #[test]
    fn test_export_format_variants() {
        assert_eq!(ExportFormat::Csv, ExportFormat::Csv);
        assert_ne!(ExportFormat::Csv, ExportFormat::Json);
        assert_ne!(ExportFormat::Binary, ExportFormat::Xyz);
    }

    // ---- ParticleExporter ----

    #[test]
    fn test_particle_exporter_new_empty() {
        let pe = ParticleExporter::new(ExportFormat::Csv);
        assert!(pe.is_empty());
        assert_eq!(pe.len(), 0);
    }

    #[test]
    fn test_particle_exporter_add_particle() {
        let mut pe = ParticleExporter::new(ExportFormat::Json);
        pe.add_particle([1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 1.0);
        assert_eq!(pe.len(), 1);
        assert!(!pe.is_empty());
    }

    #[test]
    fn test_particle_exporter_csv_header() {
        let pe = ParticleExporter::new(ExportFormat::Csv);
        let csv = pe.export_csv();
        assert!(csv.starts_with("px,py,pz,vx,vy,vz,mass\n"));
    }

    #[test]
    fn test_particle_exporter_csv_single_particle() {
        let mut pe = ParticleExporter::new(ExportFormat::Csv);
        pe.add_particle([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], 5.0);
        let csv = pe.export_csv();
        assert!(csv.contains("1,2,3"));
        assert!(csv.contains("5"));
    }

    #[test]
    fn test_particle_exporter_csv_multiple_particles() {
        let mut pe = ParticleExporter::new(ExportFormat::Csv);
        pe.add_particle([0.0; 3], [0.0; 3], 1.0);
        pe.add_particle([1.0; 3], [1.0; 3], 2.0);
        let csv = pe.export_csv();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 data rows
    }

    #[test]
    fn test_particle_exporter_json_valid() {
        let mut pe = ParticleExporter::new(ExportFormat::Json);
        pe.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0);
        let json = pe.export_json();
        assert!(json.starts_with('['));
        assert!(json.contains("pos"));
    }

    #[test]
    fn test_particle_exporter_json_empty() {
        let pe = ParticleExporter::new(ExportFormat::Json);
        assert_eq!(pe.export_json(), "[]");
    }

    #[test]
    fn test_particle_exporter_xyz_count_line() {
        let mut pe = ParticleExporter::new(ExportFormat::Xyz);
        pe.add_particle([0.0; 3], [0.0; 3], 1.0);
        pe.add_particle([1.0; 3], [0.0; 3], 1.0);
        let xyz = pe.export_xyz();
        let first_line = xyz.lines().next().unwrap();
        assert_eq!(first_line.trim(), "2");
    }

    #[test]
    fn test_particle_exporter_xyz_atom_lines() {
        let mut pe = ParticleExporter::new(ExportFormat::Xyz);
        pe.add_particle([1.0, 2.0, 3.0], [0.0; 3], 1.0);
        let xyz = pe.export_xyz();
        assert!(xyz.contains("X 1 2 3"));
    }

    // ---- MeshExporter ----

    #[test]
    fn test_mesh_exporter_empty() {
        let me = MeshExporter::new();
        assert_eq!(me.vertex_count(), 0);
        assert_eq!(me.face_count(), 0);
    }

    #[test]
    fn test_mesh_exporter_add_vertex() {
        let mut me = MeshExporter::new();
        me.add_vertex([1.0, 2.0, 3.0]);
        assert_eq!(me.vertex_count(), 1);
    }

    #[test]
    fn test_mesh_exporter_add_face() {
        let mut me = MeshExporter::new();
        me.add_vertex([0.0, 0.0, 0.0]);
        me.add_vertex([1.0, 0.0, 0.0]);
        me.add_vertex([0.0, 1.0, 0.0]);
        me.add_face([0, 1, 2]);
        assert_eq!(me.face_count(), 1);
    }

    #[test]
    fn test_mesh_exporter_obj_vertices() {
        let mut me = MeshExporter::new();
        me.add_vertex([1.0, 2.0, 3.0]);
        let obj = me.export_obj();
        assert!(obj.contains("v 1 2 3"));
    }

    #[test]
    fn test_mesh_exporter_obj_faces_one_based() {
        let mut me = MeshExporter::new();
        me.add_vertex([0.0, 0.0, 0.0]);
        me.add_vertex([1.0, 0.0, 0.0]);
        me.add_vertex([0.0, 1.0, 0.0]);
        me.add_face([0, 1, 2]);
        let obj = me.export_obj();
        assert!(obj.contains("f 1 2 3"));
    }

    #[test]
    fn test_mesh_exporter_stl_header_footer() {
        let me = MeshExporter::new();
        let stl = me.export_stl_ascii();
        assert!(stl.starts_with("solid oxiphysics"));
        assert!(stl.ends_with("endsolid oxiphysics\n"));
    }

    #[test]
    fn test_mesh_exporter_stl_has_facets() {
        let mut me = MeshExporter::new();
        me.add_vertex([0.0, 0.0, 0.0]);
        me.add_vertex([1.0, 0.0, 0.0]);
        me.add_vertex([0.0, 1.0, 0.0]);
        me.add_face([0, 1, 2]);
        let stl = me.export_stl_ascii();
        assert!(stl.contains("facet normal"));
        assert!(stl.contains("outer loop"));
        assert!(stl.contains("endloop"));
        assert!(stl.contains("endfacet"));
    }

    #[test]
    fn test_mesh_exporter_default() {
        let me = MeshExporter::default();
        assert_eq!(me.vertex_count(), 0);
    }

    // ---- StateSnapshot ----

    #[test]
    fn test_state_snapshot_new() {
        let snap = StateSnapshot::new(10, 1.0, vec![[0.0, 1.0, 2.0]], vec![[0.1, 0.2, 0.3]]);
        assert_eq!(snap.step, 10);
        assert!((snap.time - 1.0).abs() < 1e-10);
        assert_eq!(snap.n_bodies, 1);
    }

    #[test]
    fn test_state_snapshot_serialize_deserialize() {
        let snap = StateSnapshot::new(5, 0.5, vec![[1.0, 2.0, 3.0]], vec![[0.0; 3]]);
        let json = snap.serialize_json();
        let restored = StateSnapshot::deserialize_json(&json).expect("deserialise failed");
        assert_eq!(restored.step, 5);
        assert!((restored.time - 0.5).abs() < 1e-10);
        assert_eq!(restored.n_bodies, 1);
    }

    #[test]
    fn test_state_snapshot_deserialize_invalid() {
        assert!(StateSnapshot::deserialize_json("not json at all").is_none());
    }

    #[test]
    fn test_state_snapshot_n_bodies_matches_positions() {
        let snap = StateSnapshot::new(
            0,
            0.0,
            vec![[0.0; 3], [1.0; 3], [2.0; 3]],
            vec![[0.0; 3], [0.0; 3], [0.0; 3]],
        );
        assert_eq!(snap.n_bodies, 3);
    }

    // ---- TrajectoryBuffer ----

    #[test]
    fn test_trajectory_buffer_empty() {
        let tb = TrajectoryBuffer::new(100);
        assert!(tb.is_empty());
        assert_eq!(tb.len(), 0);
    }

    #[test]
    fn test_trajectory_buffer_push_frame() {
        let mut tb = TrajectoryBuffer::new(10);
        let snap = StateSnapshot::new(0, 0.0, vec![], vec![]);
        tb.push_frame(snap);
        assert_eq!(tb.len(), 1);
    }

    #[test]
    fn test_trajectory_buffer_evicts_oldest() {
        let mut tb = TrajectoryBuffer::new(3);
        for i in 0u64..5 {
            tb.push_frame(StateSnapshot::new(i, i as f64, vec![], vec![]));
        }
        assert_eq!(tb.len(), 3);
        // Oldest remaining should be step 2
        assert_eq!(tb.get_frame(0).unwrap().step, 2);
    }

    #[test]
    fn test_trajectory_buffer_get_frame_valid() {
        let mut tb = TrajectoryBuffer::new(10);
        tb.push_frame(StateSnapshot::new(7, 0.7, vec![], vec![]));
        let f = tb.get_frame(0).unwrap();
        assert_eq!(f.step, 7);
    }

    #[test]
    fn test_trajectory_buffer_get_frame_out_of_bounds() {
        let tb = TrajectoryBuffer::new(10);
        assert!(tb.get_frame(0).is_none());
    }

    #[test]
    fn test_trajectory_buffer_export_json_valid() {
        let mut tb = TrajectoryBuffer::new(10);
        tb.push_frame(StateSnapshot::new(1, 0.1, vec![], vec![]));
        let json = tb.export_trajectory_json();
        assert!(json.starts_with('['));
    }

    #[test]
    fn test_trajectory_buffer_export_json_empty() {
        let tb = TrajectoryBuffer::new(10);
        assert_eq!(tb.export_trajectory_json(), "[]");
    }

    // ---- ConfigSerializer ----

    #[test]
    fn test_config_serializer_new() {
        let cs = ConfigSerializer::new();
        assert!(cs.is_empty());
        assert_eq!(cs.len(), 0);
    }

    #[test]
    fn test_config_serializer_set_get() {
        let mut cs = ConfigSerializer::new();
        cs.set("gravity", "-9.81");
        assert_eq!(cs.get("gravity"), Some("-9.81"));
    }

    #[test]
    fn test_config_serializer_get_missing() {
        let cs = ConfigSerializer::new();
        assert!(cs.get("nonexistent").is_none());
    }

    #[test]
    fn test_config_serializer_overwrite() {
        let mut cs = ConfigSerializer::new();
        cs.set("key", "old");
        cs.set("key", "new");
        assert_eq!(cs.get("key"), Some("new"));
    }

    #[test]
    fn test_config_serializer_serialize_json() {
        let mut cs = ConfigSerializer::new();
        cs.set("dt", "0.016");
        let json = cs.serialize_json();
        assert!(json.contains("dt"));
        assert!(json.contains("0.016"));
    }

    #[test]
    fn test_config_serializer_roundtrip() {
        let mut cs = ConfigSerializer::new();
        cs.set("solver", "pgs");
        cs.set("iterations", "10");
        let json = cs.serialize_json();
        let mut cs2 = ConfigSerializer::new();
        assert!(cs2.deserialize_json(&json));
        assert_eq!(cs2.get("solver"), Some("pgs"));
        assert_eq!(cs2.get("iterations"), Some("10"));
    }

    #[test]
    fn test_config_serializer_deserialize_invalid() {
        let mut cs = ConfigSerializer::new();
        assert!(!cs.deserialize_json("bad json"));
    }

    #[test]
    fn test_config_serializer_default() {
        let cs = ConfigSerializer::default();
        assert!(cs.is_empty());
    }

    #[test]
    fn test_config_serializer_len() {
        let mut cs = ConfigSerializer::new();
        cs.set("a", "1");
        cs.set("b", "2");
        assert_eq!(cs.len(), 2);
    }

    #[test]
    fn test_particle_exporter_format_stored() {
        let pe = ParticleExporter::new(ExportFormat::Binary);
        assert_eq!(pe.format, ExportFormat::Binary);
    }

    #[test]
    fn test_mesh_exporter_obj_comment() {
        let me = MeshExporter::new();
        let obj = me.export_obj();
        assert!(obj.contains("# OxiPhysics OBJ export"));
    }
}
