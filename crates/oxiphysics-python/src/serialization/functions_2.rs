//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)
use super::types::{SimulationCheckpoint, SimulationSnapshot};

/// Serialize a complete simulation checkpoint from a snapshot.
///
/// Convenience wrapper around `SimulationCheckpoint::from_snapshot`.
pub fn serialize_simulation_checkpoint(
    snap: &SimulationSnapshot,
    label: impl Into<String>,
    step_count: u64,
    timestamp: f64,
) -> String {
    let cp = SimulationCheckpoint::from_snapshot(snap, label, step_count, timestamp);
    cp.to_json()
}
#[cfg(test)]
mod json_checkpoint_tests {
    use super::super::functions::{deserialize_body_state_json, serialize_body_state_json};
    use super::*;
    use crate::serialization::BodyStateJson;
    use crate::serialization::SimBodyState;
    #[test]
    fn test_serialize_body_state_json_roundtrip() {
        let body = SimBodyState {
            handle: 42,
            position: [1.0, 2.0, 3.0],
            velocity: [0.1, 0.2, 0.3],
            orientation: [0.0, 0.0, 0.0, 1.0],
            angular_velocity: [0.0; 3],
            is_sleeping: false,
            is_static: false,
            tag: Some("ball".to_string()),
        };
        let json = serialize_body_state_json(&body);
        let back = deserialize_body_state_json(&json).unwrap();
        assert_eq!(back.handle, 42);
        assert!((back.position[0] - 1.0).abs() < 1e-12);
        assert_eq!(back.tag.as_deref(), Some("ball"));
    }
    #[test]
    fn test_deserialize_body_state_json_invalid() {
        let result = deserialize_body_state_json("not valid json");
        assert!(result.is_err());
    }
    #[test]
    fn test_body_state_json_schema_version() {
        let body = SimBodyState::at_rest(0, [0.0; 3]);
        let bsj = BodyStateJson::from_sim_body(&body);
        assert_eq!(bsj.schema_version, "1.0.0");
    }
    #[test]
    fn test_simulation_checkpoint_from_snapshot() {
        let mut snap = SimulationSnapshot::empty();
        snap.time = 5.0;
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 2.0, 3.0]));
        let cp = SimulationCheckpoint::from_snapshot(&snap, "test_checkpoint", 300, 12345.0);
        assert_eq!(cp.label, "test_checkpoint");
        assert!((cp.sim_time - 5.0).abs() < 1e-10);
        assert_eq!(cp.step_count, 300);
        assert_eq!(cp.body_count(), 1);
    }
    #[test]
    fn test_simulation_checkpoint_json_roundtrip() {
        let mut snap = SimulationSnapshot::empty();
        snap.time = 2.5;
        snap.gravity = [0.0, -9.81, 0.0];
        snap.bodies.push(SimBodyState::at_rest(1, [0.0, 5.0, 0.0]));
        let json = serialize_simulation_checkpoint(&snap, "save_1", 100, 0.0);
        let cp = SimulationCheckpoint::from_json(&json).unwrap();
        assert_eq!(cp.body_count(), 1);
        assert!((cp.sim_time - 2.5).abs() < 1e-10);
        let restored = cp.to_snapshot();
        assert_eq!(restored.bodies.len(), 1);
        assert!((restored.bodies[0].position[1] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_diff_identical() {
        let mut snap = SimulationSnapshot::empty();
        snap.bodies.push(SimBodyState::at_rest(0, [1.0, 0.0, 0.0]));
        snap.bodies.push(SimBodyState::at_rest(1, [2.0, 0.0, 0.0]));
        let diff = snap.diff(&snap, 1e-6);
        assert!(
            diff.is_identical(),
            "identical snapshots should yield empty diff"
        );
        assert!((diff.max_displacement).abs() < 1e-12);
    }
    #[test]
    fn test_snapshot_diff_added_removed() {
        let mut snap_a = SimulationSnapshot::empty();
        snap_a.bodies.push(SimBodyState::at_rest(0, [0.0; 3]));
        snap_a
            .bodies
            .push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        let mut snap_b = SimulationSnapshot::empty();
        snap_b
            .bodies
            .push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        snap_b
            .bodies
            .push(SimBodyState::at_rest(2, [5.0, 0.0, 0.0]));
        let diff = snap_a.diff(&snap_b, 1e-6);
        assert_eq!(diff.removed.len(), 1);
        assert!(diff.removed.contains(&0));
        assert_eq!(diff.added.len(), 1);
        assert!(diff.added.contains(&2));
    }
    #[test]
    fn test_snapshot_diff_moved() {
        let mut snap_a = SimulationSnapshot::empty();
        snap_a.time = 0.0;
        snap_a
            .bodies
            .push(SimBodyState::at_rest(0, [0.0, 0.0, 0.0]));
        let mut snap_b = SimulationSnapshot::empty();
        snap_b.time = 1.0;
        let mut moved = SimBodyState::at_rest(0, [3.0, 4.0, 0.0]);
        moved.velocity = [1.0, 0.0, 0.0];
        snap_b.bodies.push(moved);
        let diff = snap_a.diff(&snap_b, 1e-6);
        assert_eq!(diff.moved.len(), 1);
        assert_eq!(diff.moved[0], 0);
        assert!(
            (diff.max_displacement - 5.0).abs() < 1e-10,
            "disp = {}",
            diff.max_displacement
        );
        assert!((diff.time_delta - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_snapshot_diff_change_count() {
        let snap_a = SimulationSnapshot::empty();
        let mut snap_b = SimulationSnapshot::empty();
        snap_b.bodies.push(SimBodyState::at_rest(0, [0.0; 3]));
        snap_b
            .bodies
            .push(SimBodyState::at_rest(1, [1.0, 0.0, 0.0]));
        let diff = snap_a.diff(&snap_b, 1e-6);
        assert_eq!(diff.change_count(), 2);
        assert!(!diff.is_identical());
    }
    #[test]
    fn test_checkpoint_empty() {
        let cp = SimulationCheckpoint::empty("init");
        assert_eq!(cp.body_count(), 0);
        assert_eq!(cp.label, "init");
        assert_eq!(cp.version, SimulationCheckpoint::FORMAT_VERSION);
    }
}
