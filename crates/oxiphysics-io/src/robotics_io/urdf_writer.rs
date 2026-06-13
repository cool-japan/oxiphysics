//! Serialise a [`UrdfRobot`] back to URDF XML text.

use super::types::{UrdfGeometry, UrdfRobot, UrdfVisualElement};

/// Serialise a [`UrdfRobot`] to URDF XML text.
///
/// Links are emitted in sorted order by name for deterministic output, and
/// joints are emitted in the robot's `joint_order`. The result round-trips
/// through the URDF parser.
pub fn write_urdf(robot: &UrdfRobot) -> String {
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
    out.push_str(&format!("<robot name=\"{}\">\n", xml_escape(&robot.name)));

    // Links (sorted by name for deterministic output).
    let mut link_names: Vec<&String> = robot.links.keys().collect();
    link_names.sort();
    for name in link_names {
        if let Some(link) = robot.links.get(name) {
            out.push_str(&format!("  <link name=\"{}\">\n", xml_escape(&link.name)));
            out.push_str("    <inertial>\n");
            out.push_str(&format!(
                "      <origin xyz=\"{}\" rpy=\"0 0 0\"/>\n",
                fmt_vec3(link.inertial.com)
            ));
            out.push_str(&format!(
                "      <mass value=\"{}\"/>\n",
                fmt_f64(link.inertial.mass)
            ));
            out.push_str(&format!(
                "      <inertia ixx=\"{}\" iyy=\"{}\" izz=\"{}\" ixy=\"{}\" ixz=\"{}\" iyz=\"{}\"/>\n",
                fmt_f64(link.inertial.inertia_diag[0]),
                fmt_f64(link.inertial.inertia_diag[1]),
                fmt_f64(link.inertial.inertia_diag[2]),
                fmt_f64(link.inertial.inertia_off[0]),
                fmt_f64(link.inertial.inertia_off[1]),
                fmt_f64(link.inertial.inertia_off[2])
            ));
            out.push_str("    </inertial>\n");
            for v in &link.visuals {
                write_visual(&mut out, v, "visual");
            }
            for c in &link.collisions {
                write_visual(&mut out, c, "collision");
            }
            out.push_str("  </link>\n");
        }
    }

    for name in &robot.joint_order {
        if let Some(joint) = robot.joints.get(name) {
            out.push_str(&format!(
                "  <joint name=\"{}\" type=\"{}\">\n",
                xml_escape(&joint.name),
                joint.joint_type
            ));
            out.push_str(&format!(
                "    <parent link=\"{}\"/>\n",
                xml_escape(&joint.parent)
            ));
            out.push_str(&format!(
                "    <child link=\"{}\"/>\n",
                xml_escape(&joint.child)
            ));
            out.push_str(&format!(
                "    <origin xyz=\"{}\" rpy=\"{}\"/>\n",
                fmt_vec3(joint.origin_xyz),
                fmt_vec3(joint.origin_rpy)
            ));
            out.push_str(&format!("    <axis xyz=\"{}\"/>\n", fmt_vec3(joint.axis)));
            out.push_str(&format!(
                "    <limit lower=\"{}\" upper=\"{}\" effort=\"{}\" velocity=\"{}\"/>\n",
                fmt_f64(joint.limits.lower),
                fmt_f64(joint.limits.upper),
                fmt_f64(joint.limits.effort),
                fmt_f64(joint.limits.velocity)
            ));
            out.push_str(&format!(
                "    <dynamics damping=\"{}\" friction=\"{}\"/>\n",
                fmt_f64(joint.damping),
                fmt_f64(joint.friction)
            ));
            out.push_str("  </joint>\n");
        }
    }
    out.push_str("</robot>\n");
    out
}

/// Write a `<visual>` or `<collision>` block (selected by `tag`).
fn write_visual(out: &mut String, ve: &UrdfVisualElement, tag: &str) {
    match &ve.name {
        Some(n) => out.push_str(&format!("    <{} name=\"{}\">\n", tag, xml_escape(n))),
        None => out.push_str(&format!("    <{}>\n", tag)),
    }
    out.push_str(&format!(
        "      <origin xyz=\"{}\" rpy=\"{}\"/>\n",
        fmt_vec3(ve.origin_xyz),
        fmt_vec3(ve.origin_rpy)
    ));
    out.push_str("      <geometry>\n");
    match &ve.geometry {
        UrdfGeometry::Box { half_extents } => out.push_str(&format!(
            "        <box size=\"{}\"/>\n",
            fmt_vec3([
                half_extents[0] * 2.0,
                half_extents[1] * 2.0,
                half_extents[2] * 2.0
            ])
        )),
        UrdfGeometry::Sphere { radius } => out.push_str(&format!(
            "        <sphere radius=\"{}\"/>\n",
            fmt_f64(*radius)
        )),
        UrdfGeometry::Cylinder {
            radius,
            half_length,
        } => out.push_str(&format!(
            "        <cylinder radius=\"{}\" length=\"{}\"/>\n",
            fmt_f64(*radius),
            fmt_f64(half_length * 2.0)
        )),
        UrdfGeometry::Mesh { filename, scale } => out.push_str(&format!(
            "        <mesh filename=\"{}\" scale=\"{}\"/>\n",
            xml_escape(filename),
            fmt_vec3(*scale)
        )),
    }
    out.push_str("      </geometry>\n");
    if let Some(m) = &ve.material {
        out.push_str(&format!("      <material name=\"{}\"/>\n", xml_escape(m)));
    }
    out.push_str(&format!("    </{}>\n", tag));
}

/// Format an `f64` for URDF output.
fn fmt_f64(v: f64) -> String {
    v.to_string()
}

/// Format a 3-vector as space-separated values.
fn fmt_vec3(v: [f64; 3]) -> String {
    format!("{} {} {}", v[0], v[1], v[2])
}

/// Escape the five predefined XML entities (`&` first).
fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

#[cfg(test)]
mod tests {
    use super::super::types::{JointLimits, JointType, UrdfJoint, UrdfLink};
    use super::super::urdf_parser::parse_urdf;
    use super::*;

    #[test]
    fn roundtrip_tiny_robot() {
        let mut robot = UrdfRobot::new("tiny");
        robot.add_link(UrdfLink::new("base"));
        robot.add_link(UrdfLink::new("arm"));
        robot.add_joint(UrdfJoint {
            name: "j1".to_string(),
            joint_type: JointType::Revolute,
            parent: "base".to_string(),
            child: "arm".to_string(),
            origin_xyz: [0.0, 0.0, 1.0],
            origin_rpy: [0.0, 0.0, 0.0],
            axis: [0.0, 0.0, 1.0],
            limits: JointLimits {
                lower: -1.0,
                upper: 1.0,
                effort: 50.0,
                velocity: 5.0,
            },
            damping: 0.1,
            friction: 0.01,
        });

        let xml = write_urdf(&robot);
        let parsed = parse_urdf(&xml).unwrap();
        assert_eq!(parsed.name, "tiny");
        assert_eq!(parsed.links.len(), 2);
        assert_eq!(parsed.joints.len(), 1);
        assert!(parsed.links.contains_key("base"));
        assert!(parsed.links.contains_key("arm"));
        let j = parsed.joints.get("j1").unwrap();
        assert_eq!(j.parent, "base");
        assert_eq!(j.child, "arm");
        assert!(matches!(j.joint_type, JointType::Revolute));
    }
}
