//! Integration tests for the URDF XML parser, writer, and articulated mapping.

use std::path::{Path, PathBuf};

fn fixture_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

#[test]
fn test_xml_parser_basic() {
    let xml = r#"<?xml version="1.0"?>
    <!-- a comment -->
    <root attr1="hello" attr2='world'>
      <child key="val&amp;ue"/>
      <child2>text</child2>
    </root>"#;
    let elem = oxiphysics_io::robotics_io::xml::parse_xml(xml).expect("parse_xml failed");
    assert_eq!(elem.name, "root");
    assert_eq!(elem.attr("attr1"), Some("hello"));
    assert_eq!(elem.attr("attr2"), Some("world"));
    assert!(elem.child("child").is_some());
    let child = elem.child("child").expect("child element present");
    assert_eq!(child.attr("key"), Some("val&ue"));
    let child2 = elem.child("child2").expect("child2 element present");
    assert_eq!(child2.text.trim(), "text");
}

#[test]
fn test_urdf_parse_two_link() {
    let path = fixture_path("two_link_pendulum.urdf");
    let xml_str = std::fs::read_to_string(&path).expect("read fixture");
    let robot = oxiphysics_io::robotics_io::parse_urdf(&xml_str).expect("parse_urdf");
    assert_eq!(robot.name, "two_link_pendulum");
    assert!(robot.links.contains_key("link1"));
    assert!(robot.links.contains_key("link2"));
    assert!(robot.joints.contains_key("joint1"));
    let j = robot.joints.get("joint1").expect("joint1 present");
    use oxiphysics_io::robotics_io::types::JointType;
    assert!(matches!(j.joint_type, JointType::Revolute));
    assert_eq!(j.axis, [0.0, 1.0, 0.0]);
    assert!((j.limits.effort - 100.0).abs() < 1e-12);
    let l1 = robot.links.get("link1").expect("link1 present");
    assert!((l1.inertial.mass - 1.0).abs() < 1e-12);
    assert_eq!(l1.inertial.com, [0.0, 0.0, 0.5]);
}

#[test]
fn test_urdf_roundtrip() {
    let path = fixture_path("two_link_pendulum.urdf");
    let xml_str = std::fs::read_to_string(&path).expect("read fixture");
    let robot = oxiphysics_io::robotics_io::parse_urdf(&xml_str).expect("parse first");
    let written = oxiphysics_io::robotics_io::write_urdf(&robot);
    let robot2 = oxiphysics_io::robotics_io::parse_urdf(&written).expect("parse second");
    assert_eq!(robot.name, robot2.name);
    assert_eq!(robot.links.len(), robot2.links.len());
    assert_eq!(robot.joints.len(), robot2.joints.len());
    for name in &robot.joint_order {
        let a = robot.joints.get(name).expect("joint in first");
        let b = robot2.joints.get(name).expect("joint in second");
        assert_eq!(a.joint_type, b.joint_type, "joint type for {name}");
        assert_eq!(a.parent, b.parent);
        assert_eq!(a.child, b.child);
    }
}

#[test]
fn test_urdf_to_articulated_two_link() {
    let path = fixture_path("two_link_pendulum.urdf");
    let xml_str = std::fs::read_to_string(&path).expect("read fixture");
    let robot = oxiphysics_io::robotics_io::parse_urdf(&xml_str).expect("parse_urdf");
    let root = oxiphysics_io::robotics_io::validate_tree(&robot).expect("validate_tree");
    assert_eq!(root, "world");
    let gravity = [0.0, 0.0, -9.81];
    let model = oxiphysics_io::robotics_io::urdf_to_articulated(&robot, &root, gravity)
        .expect("urdf_to_articulated");

    assert_eq!(model.num_bodies(), 4, "expected 4 bodies");

    let n_dof = model.total_dof();
    assert_eq!(n_dof, 2, "expected 2 DOFs from 2 revolute joints, got {n_dof}");

    let q = vec![0.0; n_dof];
    let qd = vec![0.0; n_dof];
    let qdd = vec![0.0; n_dof];
    let tau = oxiphysics_articulated::rnea::rnea(&model, &q, &qd, &qdd);
    assert_eq!(tau.len(), 2);

    assert!(
        tau[0].abs() < 1e-6,
        "joint1 gravity torque at vertical should be ~0, got {}",
        tau[0]
    );
    assert!(
        tau[1].abs() < 1e-6,
        "joint2 gravity torque at vertical should be ~0, got {}",
        tau[1]
    );
}

#[test]
fn test_urdf_to_articulated_horizontal_torque() {
    let path = fixture_path("two_link_pendulum.urdf");
    let xml_str = std::fs::read_to_string(&path).expect("read fixture");
    let robot = oxiphysics_io::robotics_io::parse_urdf(&xml_str).expect("parse_urdf");
    let root = oxiphysics_io::robotics_io::validate_tree(&robot).expect("validate_tree");
    let gravity = [0.0, 0.0, -9.81];
    let model = oxiphysics_io::robotics_io::urdf_to_articulated(&robot, &root, gravity)
        .expect("urdf_to_articulated");

    let q = vec![std::f64::consts::FRAC_PI_2, 0.0];
    let qd = vec![0.0; 2];
    let qdd = vec![0.0; 2];
    let tau = oxiphysics_articulated::rnea::rnea(&model, &q, &qd, &qdd);

    assert!(tau[0].is_finite() && tau[1].is_finite());
    assert!(
        tau[0].abs() > 1.0,
        "joint1 horizontal holding torque should be substantial, got {}",
        tau[0]
    );
}
