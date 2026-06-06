//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::Error;
use crate::types::{PySimConfig, PyVec3};
use crate::world_api::PyPhysicsWorld;

use super::types::{
    BodyStateJson, ExportBatch, IncrementalExportConfig, IncrementalUpdate, SchemaValidationResult,
    SimBodyState, SimulationSnapshot, ValidationResult, WorldState,
};

/// Compute an incremental update between two snapshots.
pub fn compute_incremental_update(
    old: &SimulationSnapshot,
    new: &SimulationSnapshot,
    sequence: u64,
) -> IncrementalUpdate {
    let old_handles: std::collections::HashSet<u32> = old.bodies.iter().map(|b| b.handle).collect();
    let new_handles: std::collections::HashSet<u32> = new.bodies.iter().map(|b| b.handle).collect();
    let removed: Vec<u32> = old_handles.difference(&new_handles).copied().collect();
    let added: Vec<u32> = new_handles.difference(&old_handles).copied().collect();
    let mut changed = Vec::new();
    for new_body in &new.bodies {
        if let Some(old_body) = old.find_body(new_body.handle) {
            let pos_changed =
                (0..3).any(|k| (new_body.position[k] - old_body.position[k]).abs() > 1e-12);
            let vel_changed =
                (0..3).any(|k| (new_body.velocity[k] - old_body.velocity[k]).abs() > 1e-12);
            if pos_changed || vel_changed {
                changed.push(new_body.clone());
            }
        } else {
            changed.push(new_body.clone());
        }
    }
    IncrementalUpdate {
        sequence,
        time: new.time,
        changed_bodies: changed,
        removed_handles: removed,
        added_handles: added,
    }
}
/// Apply an incremental update to a snapshot.
pub fn apply_incremental_update(snap: &mut SimulationSnapshot, update: &IncrementalUpdate) {
    snap.bodies
        .retain(|b| !update.removed_handles.contains(&b.handle));
    for changed in &update.changed_bodies {
        if let Some(existing) = snap.bodies.iter_mut().find(|b| b.handle == changed.handle) {
            *existing = changed.clone();
        } else {
            snap.bodies.push(changed.clone());
        }
    }
    snap.time = update.time;
}
/// Validate a deserialized snapshot for consistency.
pub fn validate_snapshot(snap: &SimulationSnapshot) -> ValidationResult {
    let mut issues = Vec::new();
    if snap.version != SimulationSnapshot::FORMAT_VERSION {
        issues.push(format!(
            "version mismatch: expected {}, got {}",
            SimulationSnapshot::FORMAT_VERSION,
            snap.version
        ));
    }
    let mut handles = std::collections::HashSet::new();
    for body in &snap.bodies {
        if !handles.insert(body.handle) {
            issues.push(format!("duplicate handle: {}", body.handle));
        }
    }
    for body in &snap.bodies {
        for k in 0..3 {
            if !body.position[k].is_finite() {
                issues.push(format!("body {}: non-finite position[{k}]", body.handle));
            }
            if !body.velocity[k].is_finite() {
                issues.push(format!("body {}: non-finite velocity[{k}]", body.handle));
            }
            if !body.angular_velocity[k].is_finite() {
                issues.push(format!(
                    "body {}: non-finite angular_velocity[{k}]",
                    body.handle
                ));
            }
        }
        for k in 0..4 {
            if !body.orientation[k].is_finite() {
                issues.push(format!("body {}: non-finite orientation[{k}]", body.handle));
            }
        }
        let qlen = body.orientation.iter().map(|x| x * x).sum::<f64>().sqrt();
        if (qlen - 1.0).abs() > 0.01 {
            issues.push(format!(
                "body {}: quaternion not normalized (len={})",
                body.handle, qlen
            ));
        }
    }
    let actual_sleeping = snap.bodies.iter().filter(|b| b.is_sleeping).count();
    if actual_sleeping != snap.sleeping_count {
        issues.push(format!(
            "sleeping_count mismatch: field says {}, actual is {}",
            snap.sleeping_count, actual_sleeping
        ));
    }
    if !snap.time.is_finite() {
        issues.push("non-finite simulation time".to_string());
    }
    for k in 0..3 {
        if !snap.gravity[k].is_finite() {
            issues.push(format!("non-finite gravity[{k}]"));
        }
    }
    ValidationResult {
        is_valid: issues.is_empty(),
        issues,
    }
}
/// Serialize the current state of `world` into a compact JSON string.
pub fn save_snapshot(world: &PyPhysicsWorld) -> String {
    build_snapshot(world).to_json()
}
/// Serialize the current state of `world` into a pretty-printed JSON string.
pub fn save_snapshot_pretty(world: &PyPhysicsWorld) -> String {
    build_snapshot(world).to_pretty_json()
}
/// Deserialize a `SimulationSnapshot` from a JSON string.
pub fn load_snapshot(json: &str) -> Result<SimulationSnapshot, Error> {
    SimulationSnapshot::from_json(json)
}
/// Build a `SimulationSnapshot` from `world` without serializing it.
pub fn build_snapshot(world: &PyPhysicsWorld) -> SimulationSnapshot {
    let handles = world.active_handles();
    let bodies: Vec<SimBodyState> = handles
        .iter()
        .map(|&h| SimBodyState {
            handle: h,
            position: world.get_position(h).unwrap_or([0.0; 3]),
            velocity: world.get_velocity(h).unwrap_or([0.0; 3]),
            orientation: world.get_orientation(h).unwrap_or([0.0, 0.0, 0.0, 1.0]),
            angular_velocity: world.get_angular_velocity(h).unwrap_or([0.0; 3]),
            is_sleeping: world.is_sleeping(h),
            is_static: false,
            tag: world.get_tag(h),
        })
        .collect();
    SimulationSnapshot {
        version: SimulationSnapshot::FORMAT_VERSION,
        time: world.time(),
        gravity: world.gravity(),
        bodies,
        contacts: world.get_contacts(),
        sleeping_count: world.sleeping_count(),
        description: None,
        metadata: std::collections::HashMap::new(),
    }
}
/// Restore world state from a snapshot (positions and velocities only).
pub fn apply_snapshot(world: &mut PyPhysicsWorld, snap: &SimulationSnapshot) -> usize {
    let mut updated = 0;
    for body_state in &snap.bodies {
        let h = body_state.handle;
        if world.get_position(h).is_some() {
            world.set_position(h, body_state.position);
            world.set_velocity(h, body_state.velocity);
            world.set_orientation(h, body_state.orientation);
            world.set_angular_velocity(h, body_state.angular_velocity);
            updated += 1;
        }
    }
    updated
}
/// Serialize a `PySimConfig` to JSON.
pub fn config_to_json(config: &PySimConfig) -> String {
    serde_json::to_string(config).unwrap_or_else(|_| "{}".to_string())
}
/// Deserialize a `PySimConfig` from JSON.
pub fn config_from_json(json: &str) -> Result<PySimConfig, Error> {
    serde_json::from_str(json)
        .map_err(|e| Error::General(format!("config deserialization failed: {e}")))
}
/// Serialize a `PyPhysicsWorld` to a JSON string (legacy).
pub fn to_json(world: &PyPhysicsWorld) -> String {
    let state = WorldState {
        gravity: PyVec3::from_array(world.gravity()),
        time: world.time(),
        positions: world
            .all_positions()
            .into_iter()
            .map(PyVec3::from_array)
            .collect(),
        num_bodies: world.body_count(),
    };
    serde_json::to_string(&state).unwrap_or_else(|_| "{}".to_string())
}
/// Deserialize a `WorldState` from a JSON string (legacy).
pub fn from_json(json: &str) -> Option<WorldState> {
    serde_json::from_str(json).ok()
}
/// Magic bytes that identify the OxiPhysics pickle envelope.
pub(super) const PICKLE_MAGIC: &[u8; 4] = b"OXPK";
/// Protocol version stored in the pickle envelope.
pub(super) const PICKLE_VERSION: u8 = 2;
/// Compute pairwise distances between all body pairs from a snapshot.
/// Returns a flat upper-triangular distance array.
pub fn compute_pairwise_distances(snap: &SimulationSnapshot) -> Vec<f64> {
    let n = snap.bodies.len();
    let mut dists = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let pi = &snap.bodies[i].position;
            let pj = &snap.bodies[j].position;
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            dists.push((dx * dx + dy * dy + dz * dz).sqrt());
        }
    }
    dists
}
/// Validate that a JSON string is a well-formed `SimulationSnapshot`.
///
/// Checks:
/// - Parses as valid JSON object
/// - Contains required top-level keys
/// - `bodies` is an array
/// - All body entries have `handle`, `position`, `velocity` fields
/// - `time` and `gravity` are finite numbers / arrays
pub fn validate_snapshot_json(json: &str) -> SchemaValidationResult {
    let value: serde_json::Value = match serde_json::from_str(json) {
        Ok(v) => v,
        Err(e) => return SchemaValidationResult::err(format!("JSON parse error: {e}")),
    };
    let obj = match value.as_object() {
        Some(o) => o,
        None => return SchemaValidationResult::err("root is not an object"),
    };
    let mut errors = Vec::new();
    for key in &["version", "time", "gravity", "bodies"] {
        if !obj.contains_key(*key) {
            errors.push(format!("missing required key: {key}"));
        }
    }
    if let Some(t) = obj.get("time").and_then(|v| v.as_f64())
        && !t.is_finite()
    {
        errors.push("time is not finite".to_string());
    }
    if let Some(g) = obj.get("gravity").and_then(|v| v.as_array()) {
        if g.len() != 3 {
            errors.push(format!("gravity must have 3 components, got {}", g.len()));
        }
        for (i, gi) in g.iter().enumerate() {
            if let Some(f) = gi.as_f64() {
                if !f.is_finite() {
                    errors.push(format!("gravity[{i}] is not finite"));
                }
            } else {
                errors.push(format!("gravity[{i}] is not a number"));
            }
        }
    } else if obj.contains_key("gravity") {
        errors.push("gravity is not an array".to_string());
    }
    if let Some(bodies) = obj.get("bodies").and_then(|v| v.as_array()) {
        for (idx, body) in bodies.iter().enumerate() {
            let body_obj = match body.as_object() {
                Some(o) => o,
                None => {
                    errors.push(format!("bodies[{idx}] is not an object"));
                    continue;
                }
            };
            for field in &["handle", "position", "velocity"] {
                if !body_obj.contains_key(*field) {
                    errors.push(format!("bodies[{idx}] missing field: {field}"));
                }
            }
            if let Some(pos) = body_obj.get("position").and_then(|v| v.as_array())
                && pos.len() != 3
            {
                errors.push(format!("bodies[{idx}].position must have 3 components"));
            }
        }
    } else if obj.contains_key("bodies") {
        errors.push("bodies is not an array".to_string());
    }
    SchemaValidationResult {
        is_valid: errors.is_empty(),
        errors,
    }
}
/// Stream-export a snapshot as a series of `ExportBatch` objects.
pub fn export_snapshot_incremental(
    snap: &SimulationSnapshot,
    config: &IncrementalExportConfig,
) -> Vec<ExportBatch> {
    let filtered: Vec<&SimBodyState> = snap
        .bodies
        .iter()
        .filter(|b| {
            if !config.include_sleeping && b.is_sleeping {
                return false;
            }
            b.speed() >= config.min_speed_threshold
        })
        .collect();
    let chunk_size = config.max_batch_size.max(1);
    let total_batches = (filtered.len() + chunk_size - 1).max(1) / chunk_size;
    filtered
        .chunks(chunk_size)
        .enumerate()
        .map(|(i, chunk)| ExportBatch {
            batch_index: i,
            is_last: i + 1 >= total_batches,
            total_batches,
            time: snap.time,
            bodies: chunk.iter().map(|b| (*b).clone()).collect(),
        })
        .collect()
}
/// Reconstruct a snapshot by merging export batches (in order).
pub fn merge_export_batches(batches: &[ExportBatch]) -> SimulationSnapshot {
    let time = batches.first().map(|b| b.time).unwrap_or(0.0);
    let bodies: Vec<SimBodyState> = batches.iter().flat_map(|b| b.bodies.clone()).collect();
    let sleeping_count = bodies.iter().filter(|b| b.is_sleeping).count();
    SimulationSnapshot {
        version: SimulationSnapshot::FORMAT_VERSION,
        time,
        gravity: [0.0, -9.81, 0.0],
        bodies,
        contacts: Vec::new(),
        sleeping_count,
        description: None,
        metadata: std::collections::HashMap::new(),
    }
}
/// Serialize a `SimBodyState` to a JSON string via `BodyStateJson`.
pub fn serialize_body_state_json(body: &SimBodyState) -> String {
    let bsj = BodyStateJson::from_sim_body(body);
    serde_json::to_string(&bsj).unwrap_or_else(|_| "{}".to_string())
}
/// Deserialize a `SimBodyState` from a JSON string produced by
/// `serialize_body_state_json`.
pub fn deserialize_body_state_json(json: &str) -> Result<SimBodyState, crate::Error> {
    let bsj: BodyStateJson = serde_json::from_str(json)
        .map_err(|e| crate::Error::General(format!("body state deserialization failed: {e}")))?;
    Ok(bsj.to_sim_body())
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::serialization::BodyDict;
    use crate::serialization::NumpyPositionArray;
    use crate::serialization::PickleEnvelope;
    use crate::serialization::SchemaVersion;
    use crate::serialization::SnapshotDict;
    use crate::types::{PyRigidBodyConfig, PyRigidBodyDesc, PySimConfig, PyVec3};
    fn make_world_with_bodies() -> PyPhysicsWorld {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        let cfg = PyRigidBodyConfig::dynamic(1.0, [1.0, 2.0, 3.0]).with_tag("body_0");
        world.add_rigid_body(cfg);
        let cfg2 = PyRigidBodyConfig::static_body([0.0, 0.0, 0.0]);
        world.add_rigid_body(cfg2);
        world
    }
    #[test]
    fn test_legacy_serialization_roundtrip() {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        let desc = PyRigidBodyDesc {
            mass: 1.0,
            position: PyVec3::new(1.0, 2.0, 3.0),
            is_static: false,
        };
        world.add_body_legacy(&desc);
        let json = to_json(&world);
        let state = from_json(&json).expect("should deserialize");
        assert_eq!(state.num_bodies, 1);
        assert!((state.positions[0].x - 1.0).abs() < 1e-10);
        assert!((state.gravity.y + 9.81).abs() < 1e-10);
    }
    #[test]
    fn test_save_load_snapshot_roundtrip() {
        let world = make_world_with_bodies();
        let json = save_snapshot(&world);
        let snap = load_snapshot(&json).expect("deserialize snapshot");
        assert_eq!(snap.bodies.len(), 2);
        assert_eq!(snap.version, SimulationSnapshot::FORMAT_VERSION);
    }
    #[test]
    fn test_snapshot_gravity_preserved() {
        let world = PyPhysicsWorld::new(PySimConfig::moon_gravity());
        let json = save_snapshot(&world);
        let snap = load_snapshot(&json).expect("deserialize");
        assert!((snap.gravity[1] + 1.62).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_time_preserved() {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        world.step(0.1);
        world.step(0.1);
        let json = save_snapshot(&world);
        let snap = load_snapshot(&json).expect("deserialize");
        assert!((snap.time - 0.2).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_body_position_preserved() {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        let cfg = PyRigidBodyConfig::dynamic(1.0, [5.0, 10.0, -3.0]);
        world.add_rigid_body(cfg);
        let snap = build_snapshot(&world);
        assert_eq!(snap.bodies.len(), 1);
        assert!((snap.bodies[0].position[0] - 5.0).abs() < 1e-10);
        assert!((snap.bodies[0].position[1] - 10.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_find_by_tag() {
        let world = make_world_with_bodies();
        let snap = build_snapshot(&world);
        let found = snap.find_by_tag("body_0");
        assert!(found.is_some());
        assert!((found.unwrap().position[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_find_by_handle() {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [7.0, 0.0, 0.0]));
        let snap = build_snapshot(&world);
        let found = snap.find_body(h);
        assert!(found.is_some());
        assert!((found.unwrap().position[0] - 7.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_total_kinetic_energy_proxy() {
        let mut world = PyPhysicsWorld::new(PySimConfig::zero_gravity());
        let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
        world.set_velocity(h, [3.0, 4.0, 0.0]);
        let snap = build_snapshot(&world);
        let ke = snap.total_kinetic_energy_proxy();
        assert!((ke - 12.5).abs() < 1e-10);
    }
    #[test]
    fn test_sim_body_state_speed() {
        let mut s = SimBodyState::at_rest(0, [0.0; 3]);
        s.velocity = [3.0, 4.0, 0.0];
        assert!((s.speed() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_apply_snapshot_updates_body() {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        let h = world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0; 3]));
        let mut snap = build_snapshot(&world);
        snap.bodies[0].position = [99.0, 0.0, 0.0];
        let updated = apply_snapshot(&mut world, &snap);
        assert_eq!(updated, 1);
        let pos = world.get_position(h).unwrap();
        assert!((pos[0] - 99.0).abs() < 1e-10);
    }
    #[test]
    fn test_config_json_roundtrip() {
        let cfg = PySimConfig::earth_gravity();
        let json = config_to_json(&cfg);
        let restored = config_from_json(&json).expect("restore config");
        assert!((restored.gravity[1] + 9.81).abs() < 1e-10);
        assert_eq!(restored.solver_iterations, cfg.solver_iterations);
    }
    #[test]
    fn test_snapshot_metadata() {
        let snap = SimulationSnapshot::empty()
            .with_metadata("author", "test")
            .with_description("unit test snapshot");
        let json = snap.to_json();
        let back = SimulationSnapshot::from_json(&json).expect("deserialize");
        assert_eq!(
            back.metadata.get("author").map(String::as_str),
            Some("test")
        );
        assert_eq!(back.description.as_deref(), Some("unit test snapshot"));
    }
    #[test]
    fn test_snapshot_pretty_json_is_valid() {
        let world = make_world_with_bodies();
        let pretty = save_snapshot_pretty(&world);
        let parsed: serde_json::Value =
            serde_json::from_str(&pretty).expect("pretty JSON must be valid");
        assert!(parsed.is_object());
    }
    #[test]
    fn test_load_snapshot_invalid_json() {
        let result = load_snapshot("this is not json");
        assert!(result.is_err());
    }
    #[test]
    fn test_snapshot_after_step() {
        let mut world = PyPhysicsWorld::new(PySimConfig::earth_gravity());
        world.add_rigid_body(PyRigidBodyConfig::dynamic(1.0, [0.0, 100.0, 0.0]));
        world.step(1.0);
        let snap = build_snapshot(&world);
        assert!(snap.bodies[0].position[1] < 100.0);
        assert!(snap.bodies[0].velocity[1] < 0.0);
    }
    #[test]
    fn test_sim_body_state_serde_roundtrip() {
        let mut s = SimBodyState::at_rest(42, [1.0, 2.0, 3.0]);
        s.velocity = [0.1, -0.2, 0.3];
        s.tag = Some("player".to_string());
        let json = serde_json::to_string(&s).expect("serialize");
        let back: SimBodyState = serde_json::from_str(&json).expect("deserialize");
        assert_eq!(back, s);
    }
    #[test]
    fn test_msgpack_roundtrip() {
        let snap = SimulationSnapshot::empty()
            .with_metadata("key", "value")
            .with_description("msgpack test");
        let bytes = snap.to_msgpack();
        let restored = SimulationSnapshot::from_msgpack(&bytes).expect("msgpack roundtrip");
        assert_eq!(restored.description.as_deref(), Some("msgpack test"));
        assert_eq!(
            restored.metadata.get("key").map(String::as_str),
            Some("value")
        );
    }
    #[test]
    fn test_msgpack_invalid_magic() {
        let bad_data = b"BAAD\x00\x00\x00\x00";
        let result = SimulationSnapshot::from_msgpack(bad_data);
        assert!(result.is_err());
    }
    #[test]
    fn test_msgpack_truncated() {
        let result = SimulationSnapshot::from_msgpack(b"OXI");
        assert!(result.is_err());
    }
    #[test]
    fn test_schema_version_current() {
        let v = SchemaVersion::current();
        assert_eq!(v.major, 1);
        assert_eq!(v.to_string_version(), "1.0.0");
    }
    #[test]
    fn test_schema_version_compatibility() {
        let v1 = SchemaVersion {
            major: 1,
            minor: 0,
            patch: 0,
        };
        let v2 = SchemaVersion {
            major: 1,
            minor: 1,
            patch: 0,
        };
        let v3 = SchemaVersion {
            major: 2,
            minor: 0,
            patch: 0,
        };
        assert!(v1.is_compatible_with(&v2));
        assert!(!v1.is_compatible_with(&v3));
    }
    #[test]
    fn test_incremental_update_empty() {
        let update = IncrementalUpdate::empty(0, 0.0);
        assert!(update.is_empty());
        assert_eq!(update.change_count(), 0);
    }
    #[test]
    fn test_incremental_update_json_roundtrip() {
        let mut update = IncrementalUpdate::empty(1, 0.5);
        update
            .changed_bodies
            .push(SimBodyState::at_rest(0, [1.0, 0.0, 0.0]));
        update.removed_handles.push(5);
        let json = update.to_json();
        let restored = IncrementalUpdate::from_json(&json).expect("roundtrip");
        assert_eq!(restored.sequence, 1);
        assert_eq!(restored.changed_bodies.len(), 1);
        assert_eq!(restored.removed_handles.len(), 1);
    }
    #[test]
    fn test_compute_incremental_update() {
        let mut old = SimulationSnapshot::empty();
        old.bodies.push(SimBodyState::at_rest(0, [0.0, 0.0, 0.0]));
        old.bodies.push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        let mut new = SimulationSnapshot::empty();
        new.time = 1.0;
        new.bodies.push(SimBodyState::at_rest(0, [0.0, 0.0, 0.0]));
        let mut moved_body = SimBodyState::at_rest(1, [2.0, 0.0, 0.0]);
        moved_body.velocity = [1.0, 0.0, 0.0];
        new.bodies.push(moved_body);
        new.bodies.push(SimBodyState::at_rest(2, [3.0, 0.0, 0.0]));
        let update = compute_incremental_update(&old, &new, 1);
        assert!(!update.is_empty());
        assert!(!update.changed_bodies.is_empty());
        assert!(update.added_handles.contains(&2));
    }
    #[test]
    fn test_apply_incremental_update() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [0.0, 0.0, 0.0]));
        snap.bodies.push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        let mut update = IncrementalUpdate::empty(1, 1.0);
        update
            .changed_bodies
            .push(SimBodyState::at_rest(0, [5.0, 0.0, 0.0]));
        update.removed_handles.push(1);
        apply_incremental_update(&mut snap, &update);
        assert_eq!(snap.bodies.len(), 1);
        assert!((snap.bodies[0].position[0] - 5.0).abs() < 1e-10);
        assert!((snap.time - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_validate_valid_snapshot() {
        let snap = SimulationSnapshot::empty();
        let result = validate_snapshot(&snap);
        assert!(
            result.is_valid,
            "empty snapshot should be valid: {:?}",
            result.issues
        );
    }
    #[test]
    fn test_validate_duplicate_handles() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [0.0; 3]));
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 0.0, 0.0]));
        let result = validate_snapshot(&snap);
        assert!(!result.is_valid);
        assert!(result.issues.iter().any(|s| s.contains("duplicate")));
    }
    #[test]
    fn test_validate_nan_position() {
        let mut snap = SimulationSnapshot::empty();
        let mut body = SimBodyState::at_rest(0, [0.0; 3]);
        body.position[0] = f64::NAN;
        snap.bodies.push(body);
        let result = validate_snapshot(&snap);
        assert!(!result.is_valid);
        assert!(result.issues.iter().any(|s| s.contains("non-finite")));
    }
    #[test]
    fn test_validate_quaternion_normalization() {
        let mut snap = SimulationSnapshot::empty();
        let mut body = SimBodyState::at_rest(0, [0.0; 3]);
        body.orientation = [0.0, 0.0, 0.0, 0.0];
        snap.bodies.push(body);
        let result = validate_snapshot(&snap);
        assert!(!result.is_valid);
    }
    #[test]
    fn test_distance_from_origin() {
        let s = SimBodyState::at_rest(0, [3.0, 4.0, 0.0]);
        assert!((s.distance_from_origin() - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_is_at_rest() {
        let mut s = SimBodyState::at_rest(0, [0.0; 3]);
        assert!(s.is_at_rest(0.1, 0.1));
        s.velocity = [1.0, 0.0, 0.0];
        assert!(!s.is_at_rest(0.1, 0.1));
    }
    #[test]
    fn test_snapshot_static_dynamic_counts() {
        let mut snap = SimulationSnapshot::empty();
        let mut b1 = SimBodyState::at_rest(0, [0.0; 3]);
        b1.is_static = true;
        snap.bodies.push(b1);
        snap.bodies.push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        assert_eq!(snap.static_body_count(), 1);
        assert_eq!(snap.dynamic_body_count(), 1);
    }
    #[test]
    fn test_snapshot_handles() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(5, [0.0; 3]));
        snap.bodies.push(SimBodyState::at_rest(10, [0.0; 3]));
        let handles = snap.handles();
        assert_eq!(handles, vec![5, 10]);
    }
    #[test]
    fn test_find_by_tag_prefix() {
        let mut snap = SimulationSnapshot::empty();
        let mut b1 = SimBodyState::at_rest(0, [0.0; 3]);
        b1.tag = Some("player_1".to_string());
        let mut b2 = SimBodyState::at_rest(1, [0.0; 3]);
        b2.tag = Some("player_2".to_string());
        let mut b3 = SimBodyState::at_rest(2, [0.0; 3]);
        b3.tag = Some("enemy_1".to_string());
        snap.bodies.push(b1);
        snap.bodies.push(b2);
        snap.bodies.push(b3);
        let players = snap.find_by_tag_prefix("player");
        assert_eq!(players.len(), 2);
    }
    #[test]
    fn test_pickle_envelope_roundtrip() {
        let snap = SimulationSnapshot::empty()
            .with_metadata("source", "test")
            .with_description("pickle test");
        let env = PickleEnvelope::new(snap.clone());
        let bytes = env.to_bytes();
        let restored = PickleEnvelope::from_bytes(&bytes).expect("pickle roundtrip");
        assert_eq!(
            restored.snapshot.description.as_deref(),
            Some("pickle test")
        );
        assert_eq!(
            restored.snapshot.metadata.get("source").map(String::as_str),
            Some("test")
        );
    }
    #[test]
    fn test_pickle_envelope_invalid_magic() {
        let bad = b"BAAD\x02\x00\x00\x00\x00{}";
        assert!(PickleEnvelope::from_bytes(bad).is_err());
    }
    #[test]
    fn test_pickle_envelope_truncated() {
        let result = PickleEnvelope::from_bytes(b"OX");
        assert!(result.is_err());
    }
    #[test]
    fn test_pickle_envelope_to_hex_non_empty() {
        let snap = SimulationSnapshot::empty();
        let env = PickleEnvelope::new(snap);
        let hex = env.to_hex();
        assert!(!hex.is_empty());
        assert!(
            hex.starts_with("4f58504b"),
            "expected OXPK prefix, got: {hex}"
        );
    }
    #[test]
    fn test_pickle_with_bodies() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 2.0, 3.0]));
        snap.bodies.push(SimBodyState::at_rest(1, [4.0, 5.0, 6.0]));
        let env = PickleEnvelope::new(snap);
        let bytes = env.to_bytes();
        let restored = PickleEnvelope::from_bytes(&bytes).unwrap();
        assert_eq!(restored.snapshot.bodies.len(), 2);
        assert!((restored.snapshot.bodies[1].position[2] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_numpy_position_array_from_snapshot() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 2.0, 3.0]));
        snap.bodies.push(SimBodyState::at_rest(1, [4.0, 5.0, 6.0]));
        let arr = NumpyPositionArray::from_snapshot(&snap);
        assert_eq!(arr.shape, [2, 3]);
        assert_eq!(arr.size(), 6);
        assert_eq!(arr.n_rows(), 2);
        assert!(arr.c_order);
        assert_eq!(arr.dtype, "float64");
        let row0 = arr.get_row(0).unwrap();
        assert!((row0[0] - 1.0).abs() < 1e-10);
        assert!((row0[2] - 3.0).abs() < 1e-10);
        let row1 = arr.get_row(1).unwrap();
        assert!((row1[0] - 4.0).abs() < 1e-10);
    }
    #[test]
    fn test_numpy_velocity_array() {
        let mut snap = SimulationSnapshot::empty();
        let mut b = SimBodyState::at_rest(0, [0.0; 3]);
        b.velocity = [1.0, 2.0, 3.0];
        snap.bodies.push(b);
        let arr = NumpyPositionArray::velocity_array(&snap);
        assert_eq!(arr.n_rows(), 1);
        let row = arr.get_row(0).unwrap();
        assert!((row[0] - 1.0).abs() < 1e-10);
        assert!((row[1] - 2.0).abs() < 1e-10);
    }
    #[test]
    fn test_numpy_position_array_get_row_out_of_range() {
        let snap = SimulationSnapshot::empty();
        let arr = NumpyPositionArray::from_snapshot(&snap);
        assert!(arr.get_row(0).is_none());
    }
    #[test]
    fn test_numpy_position_array_to_raw_bytes() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 2.0, 3.0]));
        let arr = NumpyPositionArray::from_snapshot(&snap);
        let bytes = arr.to_raw_bytes();
        assert_eq!(bytes.len(), 24);
        let first = f64::from_le_bytes(bytes[0..8].try_into().unwrap());
        assert!((first - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_numpy_position_array_json_roundtrip() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [5.0, 6.0, 7.0]));
        let arr = NumpyPositionArray::from_snapshot(&snap);
        let json = arr.to_json();
        let back: NumpyPositionArray = serde_json::from_str(&json).unwrap();
        assert_eq!(back.shape, [1, 3]);
        assert!((back.data[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_compute_pairwise_distances() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [0.0, 0.0, 0.0]));
        snap.bodies.push(SimBodyState::at_rest(1, [3.0, 4.0, 0.0]));
        snap.bodies.push(SimBodyState::at_rest(2, [0.0, 0.0, 1.0]));
        let dists = compute_pairwise_distances(&snap);
        assert_eq!(dists.len(), 3);
        assert!(
            (dists[0] - 5.0).abs() < 1e-10,
            "expected 5.0, got {}",
            dists[0]
        );
        assert!(
            (dists[1] - 1.0).abs() < 1e-10,
            "expected 1.0, got {}",
            dists[1]
        );
    }
    #[test]
    fn test_body_dict_from_and_to_sim_body() {
        let mut b = SimBodyState::at_rest(7, [1.0, 2.0, 3.0]);
        b.velocity = [0.1, 0.2, 0.3];
        b.tag = Some("test".to_string());
        b.is_sleeping = true;
        let bd = BodyDict::from_sim_body(&b);
        assert_eq!(bd.handle, 7);
        assert!((bd.pos[0] - 1.0).abs() < 1e-10);
        assert!(bd.sleeping);
        assert_eq!(bd.tag.as_deref(), Some("test"));
        let back = bd.to_sim_body();
        assert_eq!(back, b);
    }
    #[test]
    fn test_snapshot_dict_roundtrip() {
        let mut snap = SimulationSnapshot::empty();
        snap.time = 5.0;
        snap.gravity = [0.0, -9.81, 0.0];
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 2.0, 3.0]));
        let dict = SnapshotDict::from_snapshot(&snap);
        assert_eq!(dict.time, 5.0);
        assert_eq!(dict.bodies.len(), 1);
        let json = dict.to_dict_json();
        let back = SnapshotDict::from_dict_json(&json).unwrap();
        assert_eq!(back.bodies.len(), 1);
        assert!((back.time - 5.0).abs() < 1e-10);
        let back_snap = back.to_snapshot();
        assert_eq!(back_snap.bodies.len(), 1);
        assert!((back_snap.bodies[0].position[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_dict_n_contacts() {
        let snap = SimulationSnapshot::empty();
        let dict = SnapshotDict::from_snapshot(&snap);
        assert_eq!(dict.n_contacts, 0);
    }
    #[test]
    fn test_snapshot_dict_invalid_json() {
        let result = SnapshotDict::from_dict_json("not json");
        assert!(result.is_err());
    }
    #[test]
    fn test_validate_snapshot_json_valid() {
        let snap = SimulationSnapshot::empty();
        let json = snap.to_json();
        let result = validate_snapshot_json(&json);
        assert!(
            result.is_valid,
            "empty snapshot JSON should be valid: {:?}",
            result.errors
        );
    }
    #[test]
    fn test_validate_snapshot_json_invalid_json() {
        let result = validate_snapshot_json("this is not json");
        assert!(!result.is_valid);
        assert!(result.errors.iter().any(|e| e.contains("JSON parse")));
    }
    #[test]
    fn test_validate_snapshot_json_missing_key() {
        let json = r#"{"version":1,"time":0.0}"#;
        let result = validate_snapshot_json(json);
        assert!(!result.is_valid);
        assert!(result.errors.iter().any(|e| e.contains("gravity")));
    }
    #[test]
    fn test_validate_snapshot_json_invalid_bodies() {
        let json = r#"{"version":1,"time":0.0,"gravity":[0,0,0],"bodies":"not_an_array"}"#;
        let result = validate_snapshot_json(json);
        assert!(!result.is_valid);
        assert!(
            result
                .errors
                .iter()
                .any(|e| e.contains("bodies is not an array"))
        );
    }
    #[test]
    fn test_validate_snapshot_json_body_missing_field() {
        let json = r#"{"version":1,"time":0.0,"gravity":[0,0,0],"bodies":[{"handle":0}]}"#;
        let result = validate_snapshot_json(json);
        assert!(!result.is_valid);
        assert!(result.errors.iter().any(|e| e.contains("position")));
    }
    #[test]
    fn test_validate_schema_validation_result_ok() {
        let r = SchemaValidationResult::ok();
        assert!(r.is_valid);
        assert!(r.errors.is_empty());
    }
    #[test]
    fn test_validate_schema_validation_result_err() {
        let r = SchemaValidationResult::err("something wrong");
        assert!(!r.is_valid);
        assert_eq!(r.errors.len(), 1);
    }
    #[test]
    fn test_export_snapshot_incremental_basic() {
        let mut snap = SimulationSnapshot::empty();
        for i in 0..10 {
            snap.bodies
                .push(SimBodyState::at_rest(i, [i as f64, 0.0, 0.0]));
        }
        let config = IncrementalExportConfig {
            max_batch_size: 3,
            ..Default::default()
        };
        let batches = export_snapshot_incremental(&snap, &config);
        assert!(!batches.is_empty());
        assert_eq!(batches.len(), 4);
        assert!(batches.last().unwrap().is_last);
        assert!(!batches[0].is_last);
        let total_bodies: usize = batches.iter().map(|b| b.bodies.len()).sum();
        assert_eq!(total_bodies, 10);
    }
    #[test]
    fn test_export_snapshot_incremental_single_batch() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [0.0; 3]));
        let config = IncrementalExportConfig::default();
        let batches = export_snapshot_incremental(&snap, &config);
        assert_eq!(batches.len(), 1);
        assert!(batches[0].is_last);
        assert_eq!(batches[0].batch_index, 0);
    }
    #[test]
    fn test_export_snapshot_incremental_speed_filter() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [0.0; 3]));
        let mut moving = SimBodyState::at_rest(1, [0.0; 3]);
        moving.velocity = [10.0, 0.0, 0.0];
        snap.bodies.push(moving);
        let config = IncrementalExportConfig {
            min_speed_threshold: 5.0,
            ..Default::default()
        };
        let batches = export_snapshot_incremental(&snap, &config);
        let total_bodies: usize = batches.iter().map(|b| b.bodies.len()).sum();
        assert_eq!(total_bodies, 1, "only the moving body should be exported");
    }
    #[test]
    fn test_export_snapshot_incremental_exclude_sleeping() {
        let mut snap = SimulationSnapshot::empty();
        let mut sleeping = SimBodyState::at_rest(0, [0.0; 3]);
        sleeping.is_sleeping = true;
        snap.bodies.push(sleeping);
        snap.bodies.push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        let config = IncrementalExportConfig {
            include_sleeping: false,
            ..Default::default()
        };
        let batches = export_snapshot_incremental(&snap, &config);
        let total_bodies: usize = batches.iter().map(|b| b.bodies.len()).sum();
        assert_eq!(total_bodies, 1);
    }
    #[test]
    fn test_merge_export_batches() {
        let mut snap = SimulationSnapshot::empty();
        snap.time = 3.0;
        for i in 0..6 {
            snap.bodies
                .push(SimBodyState::at_rest(i, [i as f64, 0.0, 0.0]));
        }
        let config = IncrementalExportConfig {
            max_batch_size: 2,
            ..Default::default()
        };
        let batches = export_snapshot_incremental(&snap, &config);
        let merged = merge_export_batches(&batches);
        assert_eq!(merged.bodies.len(), 6);
        assert!((merged.time - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_export_batch_json_roundtrip() {
        let batch = ExportBatch {
            batch_index: 0,
            is_last: true,
            total_batches: 1,
            time: 1.0,
            bodies: vec![SimBodyState::at_rest(0, [1.0, 2.0, 3.0])],
        };
        let json = batch.to_json();
        let back = ExportBatch::from_json(&json).unwrap();
        assert_eq!(back.batch_index, 0);
        assert!(back.is_last);
        assert_eq!(back.bodies.len(), 1);
    }
    #[test]
    fn test_merge_empty_batches() {
        let merged = merge_export_batches(&[]);
        assert!(merged.bodies.is_empty());
        assert!((merged.time - 0.0).abs() < 1e-10);
    }
}
