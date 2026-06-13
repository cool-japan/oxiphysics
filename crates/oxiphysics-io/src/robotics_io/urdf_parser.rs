//! Parse URDF XML into the in-memory [`UrdfRobot`] data model.

use super::types::{
    JointLimits, JointType, UrdfGeometry, UrdfInertial, UrdfJoint, UrdfLink, UrdfRobot,
    UrdfVisualElement,
};
use super::urdf_error::{UrdfError, UrdfResult};
use super::xml::{parse_float_str, parse_vec3_str, parse_xml, XmlElement};
use std::collections::{HashMap, HashSet};

/// Parse a URDF XML document into a [`UrdfRobot`].
///
/// The root element must be `<robot>` and carry a `name` attribute. All
/// `<link>` and `<joint>` children are parsed in document order.
pub fn parse_urdf(xml_str: &str) -> UrdfResult<UrdfRobot> {
    let root = parse_xml(xml_str).map_err(UrdfError::XmlParse)?;
    if root.name != "robot" {
        return Err(UrdfError::XmlParse(format!(
            "expected root <robot>, found <{}>",
            root.name
        )));
    }
    let name = root.attr("name").ok_or_else(|| UrdfError::MissingAttribute {
        element: "robot".into(),
        attr: "name".into(),
    })?;
    let mut robot = UrdfRobot::new(name);
    for link_el in root.children_named("link") {
        robot.add_link(parse_link(link_el)?);
    }
    for joint_el in root.children_named("joint") {
        robot.add_joint(parse_joint(joint_el)?);
    }
    Ok(robot)
}

/// Parse a `<link>` element into a [`UrdfLink`].
fn parse_link(el: &XmlElement) -> UrdfResult<UrdfLink> {
    let name = el.attr("name").ok_or_else(|| UrdfError::MissingAttribute {
        element: "link".into(),
        attr: "name".into(),
    })?;
    let mut link = UrdfLink::new(name);
    if let Some(inertial_el) = el.child("inertial") {
        link.inertial = parse_inertial(inertial_el)?;
    }
    for v in el.children_named("visual") {
        link.visuals.push(parse_visual(v)?);
    }
    for c in el.children_named("collision") {
        link.collisions.push(parse_visual(c)?);
    }
    Ok(link)
}

/// Parse an `<inertial>` element into a [`UrdfInertial`].
fn parse_inertial(el: &XmlElement) -> UrdfResult<UrdfInertial> {
    let mass = match el.child("mass") {
        Some(mass_el) => req_float(mass_el, "value")?,
        None => 0.0,
    };
    // The `rpy` of an inertial origin is ignored; `UrdfInertial` has no rpy field.
    let com = opt_vec3(el.child("origin"), "xyz", [0.0; 3])?;
    let (inertia_diag, inertia_off) = match el.child("inertia") {
        Some(in_el) => {
            let ixx = req_float(in_el, "ixx")?;
            let iyy = req_float(in_el, "iyy")?;
            let izz = req_float(in_el, "izz")?;
            let ixy = opt_float(Some(in_el), "ixy", 0.0)?;
            let ixz = opt_float(Some(in_el), "ixz", 0.0)?;
            let iyz = opt_float(Some(in_el), "iyz", 0.0)?;
            ([ixx, iyy, izz], [ixy, ixz, iyz])
        }
        None => ([0.0; 3], [0.0; 3]),
    };
    Ok(UrdfInertial {
        mass,
        com,
        inertia_diag,
        inertia_off,
    })
}

/// Parse a `<visual>` or `<collision>` element into a [`UrdfVisualElement`].
fn parse_visual(el: &XmlElement) -> UrdfResult<UrdfVisualElement> {
    let geom_el = el.child("geometry").ok_or_else(|| UrdfError::MissingAttribute {
        element: "visual/collision".into(),
        attr: "geometry".into(),
    })?;
    let geometry = parse_geometry(geom_el)?;
    let mut ve = UrdfVisualElement::new(geometry);
    ve.name = el.attr("name").map(String::from);
    ve.origin_xyz = opt_vec3(el.child("origin"), "xyz", [0.0; 3])?;
    ve.origin_rpy = opt_vec3(el.child("origin"), "rpy", [0.0; 3])?;
    ve.material = el.child("material").and_then(|m| m.attr("name")).map(String::from);
    Ok(ve)
}

/// Parse a `<geometry>` element into a [`UrdfGeometry`].
fn parse_geometry(el: &XmlElement) -> UrdfResult<UrdfGeometry> {
    if let Some(b) = el.child("box") {
        let size = req_vec3(b, "size")?;
        Ok(UrdfGeometry::Box {
            half_extents: [size[0] / 2.0, size[1] / 2.0, size[2] / 2.0],
        })
    } else if let Some(s) = el.child("sphere") {
        Ok(UrdfGeometry::Sphere {
            radius: req_float(s, "radius")?,
        })
    } else if let Some(c) = el.child("cylinder") {
        let radius = req_float(c, "radius")?;
        let length = req_float(c, "length")?;
        Ok(UrdfGeometry::Cylinder {
            radius,
            half_length: length / 2.0,
        })
    } else if let Some(m) = el.child("mesh") {
        let filename = m
            .attr("filename")
            .ok_or_else(|| UrdfError::MissingAttribute {
                element: "mesh".into(),
                attr: "filename".into(),
            })?
            .to_string();
        let scale = opt_vec3(Some(m), "scale", [1.0, 1.0, 1.0])?;
        Ok(UrdfGeometry::Mesh { filename, scale })
    } else {
        Err(UrdfError::XmlParse(
            "geometry element has no known shape child".into(),
        ))
    }
}

/// Parse a `<joint>` element into a [`UrdfJoint`].
fn parse_joint(el: &XmlElement) -> UrdfResult<UrdfJoint> {
    let name = el
        .attr("name")
        .ok_or_else(|| UrdfError::MissingAttribute {
            element: "joint".into(),
            attr: "name".into(),
        })?
        .to_string();
    let type_str = el.attr("type").ok_or_else(|| UrdfError::MissingAttribute {
        element: "joint".into(),
        attr: "type".into(),
    })?;
    let joint_type = match type_str {
        "fixed" => JointType::Fixed,
        "revolute" => JointType::Revolute,
        "prismatic" => JointType::Prismatic,
        "continuous" => JointType::Continuous,
        "planar" => JointType::Planar,
        "floating" => JointType::Floating,
        other => return Err(UrdfError::UnsupportedJointType(other.to_string())),
    };
    let parent = el
        .child("parent")
        .and_then(|p| p.attr("link"))
        .ok_or_else(|| UrdfError::MissingAttribute {
            element: "joint/parent".into(),
            attr: "link".into(),
        })?
        .to_string();
    let child = el
        .child("child")
        .and_then(|c| c.attr("link"))
        .ok_or_else(|| UrdfError::MissingAttribute {
            element: "joint/child".into(),
            attr: "link".into(),
        })?
        .to_string();
    let origin_xyz = opt_vec3(el.child("origin"), "xyz", [0.0; 3])?;
    let origin_rpy = opt_vec3(el.child("origin"), "rpy", [0.0; 3])?;
    let axis = opt_vec3(el.child("axis"), "xyz", [0.0, 0.0, 1.0])?;
    let limit_el = el.child("limit");
    let lower = opt_float(limit_el, "lower", 0.0)?;
    let upper = opt_float(limit_el, "upper", 0.0)?;
    let effort = opt_float(limit_el, "effort", 0.0)?;
    let velocity = opt_float(limit_el, "velocity", 0.0)?;
    let limits = JointLimits {
        lower,
        upper,
        effort,
        velocity,
    };
    let damping = opt_float(el.child("dynamics"), "damping", 0.0)?;
    let friction = opt_float(el.child("dynamics"), "friction", 0.0)?;
    Ok(UrdfJoint {
        name,
        joint_type,
        parent,
        child,
        origin_xyz,
        origin_rpy,
        axis,
        limits,
        damping,
        friction,
    })
}

/// Read a required float attribute, erroring if absent or malformed.
fn req_float(el: &XmlElement, attr: &str) -> UrdfResult<f64> {
    let v = el.attr(attr).ok_or_else(|| UrdfError::MissingAttribute {
        element: el.name.clone(),
        attr: attr.to_string(),
    })?;
    parse_float_str(v).map_err(|reason| UrdfError::InvalidValue {
        attr: attr.to_string(),
        value: v.to_string(),
        reason,
    })
}

/// Read a required 3-vector attribute, erroring if absent or malformed.
fn req_vec3(el: &XmlElement, attr: &str) -> UrdfResult<[f64; 3]> {
    let v = el.attr(attr).ok_or_else(|| UrdfError::MissingAttribute {
        element: el.name.clone(),
        attr: attr.to_string(),
    })?;
    parse_vec3_str(v).map_err(|reason| UrdfError::InvalidValue {
        attr: attr.to_string(),
        value: v.to_string(),
        reason,
    })
}

/// Read an optional float attribute, falling back to `default` when absent.
fn opt_float(el: Option<&XmlElement>, attr: &str, default: f64) -> UrdfResult<f64> {
    match el {
        None => Ok(default),
        Some(e) => match e.attr(attr) {
            None => Ok(default),
            Some(v) => parse_float_str(v).map_err(|reason| UrdfError::InvalidValue {
                attr: attr.to_string(),
                value: v.to_string(),
                reason,
            }),
        },
    }
}

/// Read an optional 3-vector attribute, falling back to `default` when absent.
fn opt_vec3(el: Option<&XmlElement>, attr: &str, default: [f64; 3]) -> UrdfResult<[f64; 3]> {
    match el {
        None => Ok(default),
        Some(e) => match e.attr(attr) {
            None => Ok(default),
            Some(v) => parse_vec3_str(v).map_err(|reason| UrdfError::InvalidValue {
                attr: attr.to_string(),
                value: v.to_string(),
                reason,
            }),
        },
    }
}

/// Validate the kinematic tree of `robot` and return the name of its root link.
///
/// Checks that every joint references existing parent/child links, that exactly
/// one root link exists (a link that is never a child), and that the tree is
/// acyclic. Root candidates and adjacency are sorted so results are
/// deterministic regardless of `HashMap` iteration order.
pub fn validate_tree(robot: &UrdfRobot) -> UrdfResult<String> {
    let mut child_set: HashSet<&str> = HashSet::new();
    for joint in robot.joints.values() {
        if !robot.links.contains_key(&joint.parent) {
            return Err(UrdfError::UnknownLink(
                joint.parent.clone(),
                joint.name.clone(),
            ));
        }
        if !robot.links.contains_key(&joint.child) {
            return Err(UrdfError::UnknownLink(joint.child.clone(), joint.name.clone()));
        }
        child_set.insert(joint.child.as_str());
    }

    let mut candidates: Vec<String> = robot
        .links
        .keys()
        .filter(|name| !child_set.contains(name.as_str()))
        .cloned()
        .collect();
    candidates.sort();
    if candidates.is_empty() {
        return Err(UrdfError::NoRoot);
    }
    if candidates.len() >= 2 {
        return Err(UrdfError::MultipleRoots(
            candidates[0].clone(),
            candidates[1].clone(),
        ));
    }
    let root = candidates[0].clone();

    let mut adj: HashMap<&str, Vec<&str>> = HashMap::new();
    for joint in robot.joints.values() {
        adj.entry(joint.parent.as_str())
            .or_default()
            .push(joint.child.as_str());
    }
    for children in adj.values_mut() {
        children.sort_unstable();
    }

    #[derive(Clone, Copy, PartialEq)]
    enum Color {
        White,
        Gray,
        Black,
    }
    let mut color: HashMap<&str, Color> = HashMap::new();
    for name in robot.links.keys() {
        color.insert(name.as_str(), Color::White);
    }
    let mut stack: Vec<(&str, usize)> = Vec::new();
    color.insert(root.as_str(), Color::Gray);
    stack.push((root.as_str(), 0));
    while let Some(&(node, idx)) = stack.last() {
        let children = adj.get(node).map(|v| v.as_slice()).unwrap_or(&[]);
        if idx < children.len() {
            if let Some(frame) = stack.last_mut() {
                frame.1 += 1;
            }
            let next = children[idx];
            match color.get(next).copied().unwrap_or(Color::White) {
                Color::Gray => return Err(UrdfError::KinematicCycle(next.to_string())),
                Color::Black => { /* already fully explored, skip */ }
                Color::White => {
                    color.insert(next, Color::Gray);
                    stack.push((next, 0));
                }
            }
        } else {
            color.insert(node, Color::Black);
            stack.pop();
        }
    }

    Ok(root)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_small_robot() {
        let urdf = r#"<?xml version="1.0"?>
<robot name="test_bot">
  <link name="base"/>
  <link name="arm"/>
  <joint name="j1" type="revolute">
    <parent link="base"/>
    <child link="arm"/>
    <axis xyz="0 0 1"/>
    <limit lower="-1.0" upper="1.0" effort="50" velocity="5"/>
  </joint>
</robot>"#;
        let robot = parse_urdf(urdf).unwrap();
        assert_eq!(robot.name, "test_bot");
        assert!(robot.links.contains_key("base"));
        assert!(robot.links.contains_key("arm"));
        assert!(robot.joints.contains_key("j1"));
        let root = validate_tree(&robot).unwrap();
        assert_eq!(root, "base");
    }

    #[test]
    fn parses_joint_fields() {
        let urdf = r#"<robot name="r">
  <link name="base"/>
  <link name="arm"/>
  <joint name="j1" type="revolute">
    <parent link="base"/>
    <child link="arm"/>
    <axis xyz="0 0 1"/>
    <limit lower="-1.0" upper="1.0" effort="50" velocity="5"/>
    <dynamics damping="0.1" friction="0.2"/>
  </joint>
</robot>"#;
        let robot = parse_urdf(urdf).unwrap();
        let j1 = robot.joints.get("j1").expect("j1 present");
        assert_eq!(j1.joint_type, JointType::Revolute);
        assert_eq!(j1.parent, "base");
        assert_eq!(j1.child, "arm");
        assert_eq!(j1.axis, [0.0, 0.0, 1.0]);
        assert_eq!(j1.limits.lower, -1.0);
        assert_eq!(j1.limits.upper, 1.0);
        assert_eq!(j1.limits.effort, 50.0);
        assert_eq!(j1.limits.velocity, 5.0);
        assert_eq!(j1.damping, 0.1);
        assert_eq!(j1.friction, 0.2);
    }

    #[test]
    fn parses_geometry_and_inertial() {
        let urdf = r#"<robot name="g">
  <link name="base">
    <inertial>
      <mass value="2.0"/>
      <origin xyz="0.1 0.2 0.3"/>
      <inertia ixx="1.0" iyy="2.0" izz="3.0" ixy="0.5" ixz="0.6" iyz="0.7"/>
    </inertial>
    <visual name="v0">
      <geometry>
        <box size="2.0 4.0 6.0"/>
      </geometry>
      <material name="red"/>
    </visual>
    <collision>
      <geometry>
        <cylinder radius="0.5" length="2.0"/>
      </geometry>
    </collision>
  </link>
</robot>"#;
        let robot = parse_urdf(urdf).unwrap();
        let base = robot.links.get("base").expect("base present");
        assert_eq!(base.inertial.mass, 2.0);
        assert_eq!(base.inertial.com, [0.1, 0.2, 0.3]);
        assert_eq!(base.inertial.inertia_diag, [1.0, 2.0, 3.0]);
        assert_eq!(base.inertial.inertia_off, [0.5, 0.6, 0.7]);
        assert_eq!(base.visuals.len(), 1);
        match base.visuals[0].geometry {
            UrdfGeometry::Box { half_extents } => {
                assert_eq!(half_extents, [1.0, 2.0, 3.0]);
            }
            _ => panic!("expected box geometry"),
        }
        assert_eq!(base.visuals[0].name.as_deref(), Some("v0"));
        assert_eq!(base.visuals[0].material.as_deref(), Some("red"));
        assert_eq!(base.collisions.len(), 1);
        match base.collisions[0].geometry {
            UrdfGeometry::Cylinder {
                radius,
                half_length,
            } => {
                assert_eq!(radius, 0.5);
                assert_eq!(half_length, 1.0);
            }
            _ => panic!("expected cylinder geometry"),
        }
    }

    #[test]
    fn missing_robot_name_errors() {
        let err = parse_urdf(r#"<robot><link name="a"/></robot>"#);
        assert!(matches!(err, Err(UrdfError::MissingAttribute { .. })));
    }

    #[test]
    fn unsupported_joint_type_errors() {
        let urdf = r#"<robot name="r">
  <link name="base"/>
  <link name="arm"/>
  <joint name="j1" type="screw">
    <parent link="base"/>
    <child link="arm"/>
  </joint>
</robot>"#;
        assert!(matches!(
            parse_urdf(urdf),
            Err(UrdfError::UnsupportedJointType(_))
        ));
    }

    #[test]
    fn detects_kinematic_cycle() {
        // A pure cycle a->b->c->a has no root (every node is a child), so we
        // attach a real root link and a joint root->a to make the cycle
        // reachable and thus detectable by the DFS.
        let mut robot = UrdfRobot::new("cyclic");
        for name in ["root", "a", "b", "c"] {
            robot.add_link(UrdfLink::new(name));
        }
        let make_joint = |name: &str, parent: &str, child: &str| UrdfJoint {
            name: name.to_string(),
            joint_type: JointType::Revolute,
            parent: parent.to_string(),
            child: child.to_string(),
            origin_xyz: [0.0; 3],
            origin_rpy: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            limits: JointLimits {
                lower: 0.0,
                upper: 0.0,
                effort: 0.0,
                velocity: 0.0,
            },
            damping: 0.0,
            friction: 0.0,
        };
        robot.add_joint(make_joint("jr", "root", "a"));
        robot.add_joint(make_joint("j1", "a", "b"));
        robot.add_joint(make_joint("j2", "b", "c"));
        robot.add_joint(make_joint("j3", "c", "a"));
        assert!(matches!(
            validate_tree(&robot),
            Err(UrdfError::KinematicCycle(_))
        ));
    }

    #[test]
    fn unknown_link_errors() {
        let mut robot = UrdfRobot::new("u");
        robot.add_link(UrdfLink::new("base"));
        robot.add_joint(UrdfJoint {
            name: "j1".to_string(),
            joint_type: JointType::Fixed,
            parent: "base".to_string(),
            child: "missing".to_string(),
            origin_xyz: [0.0; 3],
            origin_rpy: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            limits: JointLimits {
                lower: 0.0,
                upper: 0.0,
                effort: 0.0,
                velocity: 0.0,
            },
            damping: 0.0,
            friction: 0.0,
        });
        assert!(matches!(
            validate_tree(&robot),
            Err(UrdfError::UnknownLink(_, _))
        ));
    }
}
