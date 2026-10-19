//! Round-trip and robustness tests for the scene, checkpoint and trajectory
//! formats of [`simulation_io`](super).

use super::*;
use std::f64::consts::FRAC_1_SQRT_2;

/// Round a value through `f32`, the precision the trajectory format stores.
fn f32_round(value: f64) -> f64 {
    f64::from(value as f32)
}

fn scene_to_string(scene: &PhysicsSceneWriter) -> String {
    let mut buf = Vec::new();
    scene.write_json(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

fn checkpoint_to_string(checkpoint: &CheckpointWriter) -> String {
    let mut buf = Vec::new();
    checkpoint.write_json(&mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

fn dynamic_body(id: u64) -> BodyDesc {
    BodyDesc {
        id,
        position: [1.5, -2.25, 1.0e-3],
        orientation: [FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2, 0.0],
        linear_velocity: [1.0, 2.0, 3.0],
        angular_velocity: [-0.5, 0.25, 4.0],
        mass: 5.0,
        shape_tag: "box".to_string(),
    }
}

// ── Scene round trips ───────────────────────────────────────────────────────

#[test]
fn scene_reader_restores_orientation_velocities_mass_and_shape() {
    // Regression: the reader parsed only `id` and `pos` and rebuilt every body
    // with BodyDesc::new, so a restart put all bodies at rest with unit mass.
    let mut writer = PhysicsSceneWriter::new();
    writer.add_body(dynamic_body(7));
    let scene = PhysicsSceneReader::read_json(&scene_to_string(&writer)).unwrap();
    assert_eq!(scene.bodies, vec![dynamic_body(7)]);
    let body = &scene.bodies[0];
    assert_eq!(body.orientation, [FRAC_1_SQRT_2, 0.0, FRAC_1_SQRT_2, 0.0]);
    assert_eq!(body.linear_velocity, [1.0, 2.0, 3.0]);
    assert_eq!(body.mass, 5.0);
    assert_eq!(body.shape_tag, "box");
}

#[test]
fn scene_reader_restores_joints_and_materials() {
    let mut writer = PhysicsSceneWriter::new();
    writer.add_body(dynamic_body(1));
    writer.add_body(BodyDesc::new(2, [0.0, 3.0, 0.0]));
    let mut hinge = JointDesc::new(10, 1, 2, "hinge");
    hinge.anchor_a = [0.5, 0.0, -0.5];
    hinge.anchor_b = [-1.0, 2.0, 0.125];
    writer.add_joint(hinge.clone());
    let steel = MaterialDesc::new(3, 0.3, 0.6, 7850.0);
    writer.add_material(steel.clone());
    let scene = PhysicsSceneReader::read_json(&scene_to_string(&writer)).unwrap();
    assert_eq!(scene.bodies, writer.bodies);
    assert_eq!(scene.joints, vec![hinge]);
    assert_eq!(scene.materials, vec![steel]);
}

#[test]
fn scene_roundtrip_is_bit_exact_for_awkward_floats() {
    let awkward = [
        0.1,
        1.0 / 3.0,
        f64::MIN_POSITIVE,
        5.0e-324,
        f64::MAX,
        -1.0e-200,
        123_456_789.123_456_79,
    ];
    for &value in &awkward {
        let mut body = BodyDesc::new(1, [value, -value, value * 0.5]);
        body.mass = value.abs();
        let mut writer = PhysicsSceneWriter::new();
        writer.add_body(body.clone());
        let scene = PhysicsSceneReader::read_json(&scene_to_string(&writer)).unwrap();
        for k in 0..3 {
            assert_eq!(
                scene.bodies[0].position[k].to_bits(),
                body.position[k].to_bits()
            );
        }
        assert_eq!(scene.bodies[0].mass.to_bits(), body.mass.to_bits());
    }
}

#[test]
fn scene_roundtrip_preserves_signed_zero_and_non_finite_values() {
    // JSON cannot spell these; they travel as strings so that a checkpoint of
    // a diverged run can still be read back.
    let mut body = BodyDesc::new(4, [f64::NAN, f64::INFINITY, f64::NEG_INFINITY]);
    body.mass = -0.0;
    body.angular_velocity = [0.0, -0.0, f64::INFINITY];
    let mut writer = PhysicsSceneWriter::new();
    writer.add_body(body);
    let text = scene_to_string(&writer);
    let value: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(value["bodies"][0]["pos"][0], "NaN");
    assert_eq!(value["bodies"][0]["pos"][1], "Infinity");
    assert_eq!(value["bodies"][0]["pos"][2], "-Infinity");
    let back = &PhysicsSceneReader::read_json(&text).unwrap().bodies[0];
    assert!(back.position[0].is_nan());
    assert_eq!(back.position[1], f64::INFINITY);
    assert_eq!(back.position[2], f64::NEG_INFINITY);
    assert!(back.mass == 0.0 && back.mass.is_sign_negative());
    assert!(back.angular_velocity[0].is_sign_positive());
    assert!(back.angular_velocity[1] == 0.0 && back.angular_velocity[1].is_sign_negative());
    assert_eq!(back.angular_velocity[2], f64::INFINITY);
}

#[test]
fn strings_are_escaped_in_scene_documents() {
    let mut body = BodyDesc::new(1, [0.0; 3]);
    body.shape_tag = "tri\"mesh\\v2\n\tünï".to_string();
    let mut writer = PhysicsSceneWriter::new();
    writer.add_body(body.clone());
    writer.add_joint(JointDesc::new(2, 1, 1, "ball \"socket\""));
    let text = scene_to_string(&writer);
    // Still a valid JSON document ...
    serde_json::from_str::<serde_json::Value>(&text).unwrap();
    // ... that reads back unchanged.
    let scene = PhysicsSceneReader::read_json(&text).unwrap();
    assert_eq!(scene.bodies[0].shape_tag, body.shape_tag);
    assert_eq!(scene.joints[0].joint_type, "ball \"socket\"");
}

#[test]
fn hand_written_minimal_bodies_use_the_documented_defaults() {
    let json = r#"{"bodies":[{"id":3,"pos":[1,2,3]}]}"#;
    let scene = PhysicsSceneReader::read_json(json).unwrap();
    assert_eq!(scene.bodies, vec![BodyDesc::new(3, [1.0, 2.0, 3.0])]);
    assert!(scene.joints.is_empty());
    assert!(scene.materials.is_empty());
}

#[test]
fn joints_written_without_anchors_are_still_readable() {
    let json = r#"{"joints":[{"id":1,"bodyA":4,"bodyB":5,"type":"hinge"}],"extra":true}"#;
    let scene = PhysicsSceneReader::read_json(json).unwrap();
    assert_eq!(scene.joints, vec![JointDesc::new(1, 4, 5, "hinge")]);
}

#[test]
fn scene_reader_reports_malformed_documents() {
    let bad = [
        "not json at all",
        "",
        r#"{"bodies":[{"pos":[0,0,0]}]}"#,
        r#"{"bodies":[{"id":1,"pos":[1,2]}]}"#,
        r#"{"bodies":[{"id":1,"pos":[1,2,3,4]}]}"#,
        r#"{"bodies":[{"id":1,"pos":[1,2,"x"]}]}"#,
        r#"{"bodies":[{"id":-1,"pos":[1,2,3]}]}"#,
        r#"{"joints":[{"id":1,"bodyA":1,"type":"hinge"}]}"#,
        r#"{"materials":[{"id":1,"restitution":0.5}]}"#,
        r#"{"bodies":[]} trailing"#,
    ];
    for json in bad {
        let err = PhysicsSceneReader::read_json(json).unwrap_err();
        assert!(matches!(err, crate::Error::Parse(_)), "{json}: {err}");
    }
}

#[test]
fn scene_file_roundtrip_through_a_temp_file() {
    let mut writer = PhysicsSceneWriter::new();
    writer.add_body(dynamic_body(9));
    let path = std::env::temp_dir().join(format!(
        "oxiphysics_scene_roundtrip_{}.json",
        std::process::id()
    ));
    writer.write_to_file(path.to_str().unwrap_or("")).unwrap();
    let text = std::fs::read_to_string(&path).unwrap();
    let _ = std::fs::remove_file(&path);
    let scene = PhysicsSceneReader::read_json(&text).unwrap();
    assert_eq!(scene.bodies, vec![dynamic_body(9)]);
}

// ── Checkpoints ─────────────────────────────────────────────────────────────

#[test]
fn checkpoint_restores_dynamic_state_and_ordered_metadata() {
    let mut writer = CheckpointWriter::new(12.625, 4096);
    writer.add_body(dynamic_body(1));
    writer.add_body(dynamic_body(2));
    // Deliberately not in alphabetical order, with a repeated key.
    writer.add_metadata("solver", "PGS");
    writer.add_metadata("build", "release \"fast\"");
    writer.add_metadata("author", "ünï\tcode");
    writer.add_metadata("build", "second");
    let text = checkpoint_to_string(&writer);
    serde_json::from_str::<serde_json::Value>(&text).unwrap();
    let checkpoint = CheckpointReader::read_json(&text).unwrap();
    assert_eq!(checkpoint.sim_time, 12.625);
    assert_eq!(checkpoint.step, 4096);
    assert_eq!(checkpoint.bodies, writer.bodies);
    assert_eq!(checkpoint.metadata, writer.metadata);
}

#[test]
fn checkpoint_reader_keeps_file_order_and_duplicates() {
    let json = r#"{"sim_time":0,"step":0,"metadata":{"z":"1","a":"2","z":"3"}}"#;
    let checkpoint = CheckpointReader::read_json(json).unwrap();
    let pairs: Vec<(&str, &str)> = checkpoint
        .metadata
        .iter()
        .map(|(k, v)| (k.as_str(), v.as_str()))
        .collect();
    assert_eq!(pairs, vec![("z", "1"), ("a", "2"), ("z", "3")]);
    assert!(checkpoint.bodies.is_empty());
}

#[test]
fn checkpoint_with_diverged_time_can_still_be_read() {
    let writer = CheckpointWriter::new(f64::NAN, 5);
    let checkpoint = CheckpointReader::read_json(&checkpoint_to_string(&writer)).unwrap();
    assert!(checkpoint.sim_time.is_nan());
    assert_eq!(checkpoint.step, 5);
}

#[test]
fn checkpoint_reader_requires_time_and_step() {
    for json in [
        r#"{"step":3}"#,
        r#"{"sim_time":1.5}"#,
        r#"{"sim_time":"fast","step":3}"#,
        r#"{"sim_time":1.5,"step":-3}"#,
        r#"{"sim_time":1.5,"step":3,"metadata":{"k":5}}"#,
        "[]",
        "",
    ] {
        let err = CheckpointReader::read_json(json).unwrap_err();
        assert!(matches!(err, crate::Error::Parse(_)), "{json}: {err}");
    }
}

// ── Trajectories ────────────────────────────────────────────────────────────

fn write_trajectory(frames: Vec<TrajectoryFrame>) -> Vec<u8> {
    let mut writer = TrajectoryWriter::new();
    for frame in frames {
        writer.push_frame(frame);
    }
    let mut buf = Vec::new();
    writer.write_binary(&mut buf).unwrap();
    buf
}

fn positions_only_frame(time: f64, n: usize) -> TrajectoryFrame {
    let mut frame = TrajectoryFrame::new(time);
    for i in 0..n {
        frame.positions.push([i as f64 + 0.5, -(i as f64), time]);
    }
    frame
}

fn assert_frame_matches(got: &TrajectoryFrame, want: &TrajectoryFrame) {
    assert_eq!(got.time.to_bits(), want.time.to_bits());
    assert_eq!(got.positions.len(), want.positions.len());
    assert_eq!(got.velocities.len(), want.velocities.len());
    assert_eq!(got.quaternions.len(), want.quaternions.len());
    for (g, w) in got.positions.iter().zip(&want.positions) {
        assert_eq!(*g, w.map(f32_round));
    }
    for (g, w) in got.velocities.iter().zip(&want.velocities) {
        assert_eq!(*g, w.map(f32_round));
    }
    for (g, w) in got.quaternions.iter().zip(&want.quaternions) {
        assert_eq!(*g, w.map(f32_round));
    }
}

#[test]
fn trajectory_positions_only_frames_roundtrip() {
    // Regression: the writer emitted velocities and quaternions of whatever
    // length but a single count, so a positions-only frame (the
    // TrajectoryFrame::new default) shifted every later frame.
    let frames = vec![
        positions_only_frame(0.0, 3),
        positions_only_frame(0.5, 3),
        positions_only_frame(1.0, 2),
    ];
    let buf = write_trajectory(frames.clone());
    let mut reader = TrajectoryReader::from_bytes(&buf).unwrap();
    assert_eq!(reader.num_frames(), 3);
    for want in &frames {
        let got = reader.next_frame().unwrap();
        assert_frame_matches(got, want);
        assert!(got.velocities.is_empty() && got.quaternions.is_empty());
    }
    assert!(reader.next_frame().is_none());
}

#[test]
fn trajectory_frames_with_different_optional_blocks_roundtrip() {
    let mut full = positions_only_frame(0.0, 2);
    full.velocities = vec![[0.1, 0.2, 0.3], [-1.0, 0.0, 1.0e-3]];
    full.quaternions = vec![[1.0, 0.0, 0.0, 0.0], [0.5, 0.5, 0.5, 0.5]];
    let mut with_velocities = positions_only_frame(1.0, 2);
    with_velocities.velocities = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
    let mut with_quaternions = positions_only_frame(2.0, 1);
    with_quaternions.quaternions = vec![[0.0, 1.0, 0.0, 0.0]];
    let empty = TrajectoryFrame::new(3.0);
    let frames = vec![
        full,
        positions_only_frame(0.5, 4),
        with_velocities,
        with_quaternions,
        empty,
    ];
    let buf = write_trajectory(frames.clone());
    let mut reader = TrajectoryReader::from_bytes(&buf).unwrap();
    for want in &frames {
        assert_frame_matches(reader.next_frame().unwrap(), want);
    }
}

#[test]
fn trajectory_writer_rejects_inconsistent_frames_without_writing() {
    let mut frame = positions_only_frame(0.0, 2);
    frame.velocities = vec![[0.0; 3]];
    let mut writer = TrajectoryWriter::new();
    writer.push_frame(positions_only_frame(0.0, 1));
    writer.push_frame(frame);
    let mut buf = Vec::new();
    let err = writer.write_binary(&mut buf).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
    assert!(buf.is_empty(), "nothing may be written for a bad frame");

    let mut frame = positions_only_frame(0.0, 2);
    frame.quaternions = vec![[1.0, 0.0, 0.0, 0.0]; 3];
    let mut writer = TrajectoryWriter::new();
    writer.push_frame(frame);
    let err = writer.write_binary(&mut Vec::new()).unwrap_err();
    assert_eq!(err.kind(), std::io::ErrorKind::InvalidInput);
}

fn header_with_one_frame(n_pos: u64, n_vel: u64, n_quat: u64) -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend_from_slice(b"TRAJ\x01");
    buf.extend_from_slice(&1u64.to_le_bytes());
    buf.extend_from_slice(&0.25f64.to_le_bytes());
    buf.extend_from_slice(&n_pos.to_le_bytes());
    buf.extend_from_slice(&n_vel.to_le_bytes());
    buf.extend_from_slice(&n_quat.to_le_bytes());
    buf
}

#[test]
fn trajectory_reader_rejects_huge_header_counts() {
    // Regression: `Vec::with_capacity(n_frames)` with an untrusted u64 panicked
    // on capacity overflow, and per-frame counts drove unchecked arithmetic.
    let mut frames_max = Vec::new();
    frames_max.extend_from_slice(b"TRAJ\x01");
    frames_max.extend_from_slice(&u64::MAX.to_le_bytes());
    assert!(TrajectoryReader::from_bytes(&frames_max).is_err());

    for n_pos in [u64::MAX, 1u64 << 61, 1u64 << 40, 1u64 << 33] {
        let buf = header_with_one_frame(n_pos, 0, 0);
        let err = TrajectoryReader::from_bytes(&buf).err();
        assert!(
            matches!(err, Some(crate::Error::Parse(_))),
            "n_pos = {n_pos}: {err:?}"
        );
    }
    for n_quat in [u64::MAX, u64::MAX / 16 + 1] {
        let buf = header_with_one_frame(n_quat, 0, n_quat);
        assert!(TrajectoryReader::from_bytes(&buf).is_err());
    }
}

#[test]
fn trajectory_reader_rejects_every_truncation() {
    let mut full = positions_only_frame(0.0, 2);
    full.velocities = vec![[0.5; 3]; 2];
    full.quaternions = vec![[1.0, 0.0, 0.0, 0.0]; 2];
    let buf = write_trajectory(vec![full, positions_only_frame(1.0, 3)]);
    assert!(TrajectoryReader::from_bytes(&buf).is_ok());
    for len in 0..buf.len() {
        assert!(
            TrajectoryReader::from_bytes(&buf[..len]).is_err(),
            "a {len}-byte prefix of a {}-byte file was accepted",
            buf.len()
        );
    }
}

#[test]
fn trajectory_reader_rejects_trailing_bytes_and_bad_magic() {
    let mut buf = write_trajectory(vec![positions_only_frame(0.0, 1)]);
    buf.push(0);
    assert!(TrajectoryReader::from_bytes(&buf).is_err());

    let mut legacy = write_trajectory(vec![positions_only_frame(0.0, 1)]);
    legacy[4] = 0;
    let message = TrajectoryReader::from_bytes(&legacy)
        .err()
        .map(|e| e.to_string());
    assert!(
        message.as_deref().is_some_and(|m| m.contains("legacy")),
        "{message:?}"
    );

    let mut garbage = write_trajectory(vec![positions_only_frame(0.0, 1)]);
    garbage[0] = b'X';
    assert!(TrajectoryReader::from_bytes(&garbage).is_err());
}

#[test]
fn trajectory_reader_rejects_inconsistent_counts() {
    let mut buf = header_with_one_frame(2, 1, 0);
    buf.extend_from_slice(&[0u8; 2 * 12 + 12]);
    assert!(TrajectoryReader::from_bytes(&buf).is_err());
    let mut buf = header_with_one_frame(2, 0, 3);
    buf.extend_from_slice(&[0u8; 2 * 12 + 3 * 16]);
    assert!(TrajectoryReader::from_bytes(&buf).is_err());
}

#[test]
fn empty_trajectory_roundtrips() {
    let buf = write_trajectory(Vec::new());
    assert_eq!(buf.len(), 13);
    let reader = TrajectoryReader::from_bytes(&buf).unwrap();
    assert_eq!(reader.num_frames(), 0);
}
