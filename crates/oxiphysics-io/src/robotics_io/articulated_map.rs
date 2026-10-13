//! Map a validated [`UrdfRobot`] onto an `ArticulatedModel` (Featherstone tree).

use std::collections::HashMap;

use oxiphysics_articulated::body::RigidBody;
use oxiphysics_articulated::joint::{
    FixedJoint, FreeFloatingJoint, Joint, PrismaticJoint, RevoluteJoint,
};
use oxiphysics_articulated::model::ArticulatedModel;
use oxiphysics_articulated::spatial::{SpatialInertia, SpatialTransform};

use super::functions::rpy_to_mat3;
use super::types::{JointType, UrdfInertial, UrdfJoint, UrdfRobot};
use super::urdf_error::{UrdfError, UrdfResult};

/// Build an [`ArticulatedModel`] (Featherstone kinematic tree) from a validated
/// [`UrdfRobot`].
///
/// The tree is rooted at `root_link` and expanded breadth-first. Children of
/// each link are visited in lexicographic order so that body indices are fully
/// deterministic regardless of the `HashMap` iteration order. The root link is
/// attached to the world with a [`FixedJoint`]; every other link is attached to
/// its parent via the joint whose `child` field names it.
///
/// `gravity` is forwarded to [`ArticulatedModel::new`].
///
/// # Errors
///
/// Returns [`UrdfError::UnknownLink`] if a link or its parent joint cannot be
/// resolved, and [`UrdfError::UnsupportedJointType`] for joint kinds that the
/// articulated mapping does not yet support (currently `planar`).
pub fn urdf_to_articulated(
    robot: &UrdfRobot,
    root_link: &str,
    gravity: [f64; 3],
) -> UrdfResult<ArticulatedModel> {
    // 1. Map each child link to the joint that connects it to its parent.
    //    Iterating `joint_order` keeps the mapping deterministic.
    let mut child_to_joint: HashMap<&str, &UrdfJoint> = HashMap::new();
    for name in &robot.joint_order {
        if let Some(j) = robot.joints.get(name) {
            child_to_joint.insert(j.child.as_str(), j);
        }
    }

    // 2. Build a parent -> children adjacency list, sorted for determinism.
    let mut adj: HashMap<&str, Vec<&str>> = HashMap::new();
    for name in &robot.joint_order {
        if let Some(j) = robot.joints.get(name) {
            adj.entry(j.parent.as_str())
                .or_default()
                .push(j.child.as_str());
        }
    }
    for children in adj.values_mut() {
        children.sort_unstable();
    }

    // 3. Breadth-first traversal from the root, recording the visit order.
    let mut order: Vec<String> = Vec::new();
    let mut link_to_index: HashMap<String, usize> = HashMap::new();
    let mut queue: std::collections::VecDeque<&str> = std::collections::VecDeque::new();
    queue.push_back(root_link);
    while let Some(name) = queue.pop_front() {
        if link_to_index.contains_key(name) {
            continue;
        }
        let idx = order.len();
        link_to_index.insert(name.to_string(), idx);
        order.push(name.to_string());
        if let Some(children) = adj.get(name) {
            for c in children {
                queue.push_back(*c);
            }
        }
    }

    // 4. Create the model with the requested gravity.
    let mut model = ArticulatedModel::new(gravity);

    // 5. Add bodies in BFS order so parent indices always precede children.
    for (idx, name) in order.iter().enumerate() {
        let link = robot
            .links
            .get(name)
            .ok_or_else(|| UrdfError::UnknownLink(name.clone(), "model".into()))?;
        let inertia = build_spatial_inertia(&link.inertial);
        if idx == 0 {
            let joint: Box<dyn Joint> = Box::new(FixedJoint);
            model.add_body(
                RigidBody::new(name.clone(), inertia, None, SpatialTransform::IDENTITY),
                joint,
            );
        } else {
            let pj = child_to_joint
                .get(name.as_str())
                .ok_or_else(|| UrdfError::UnknownLink(name.clone(), "no parent joint".into()))?;
            let parent_id = *link_to_index
                .get(&pj.parent)
                .ok_or_else(|| UrdfError::UnknownLink(pj.parent.clone(), pj.name.clone()))?;
            let parent_transform = spatial_transform_from_origin(pj.origin_xyz, pj.origin_rpy);
            let joint: Box<dyn Joint> = match pj.joint_type {
                JointType::Fixed => Box::new(FixedJoint),
                JointType::Revolute | JointType::Continuous => {
                    Box::new(RevoluteJoint::new(pj.axis))
                }
                JointType::Prismatic => Box::new(PrismaticJoint::new(pj.axis)),
                JointType::Floating => Box::new(FreeFloatingJoint),
                JointType::Planar => return Err(UrdfError::UnsupportedJointType("planar".into())),
            };
            model.add_body(
                RigidBody::new(name.clone(), inertia, Some(parent_id), parent_transform),
                joint,
            );
        }
    }

    // 6. The fully-populated Featherstone tree.
    Ok(model)
}

/// Build a child->parent [`SpatialTransform`] from a URDF `<origin xyz rpy>`.
///
/// `rpy` is the extrinsic-XYZ / intrinsic-ZYX Euler convention used by URDF and
/// is converted with [`rpy_to_mat3`], the same routine the URDF parser uses to
/// bring `<inertial>` tensors into the link frame. `xyz` becomes the
/// translation directly.
fn spatial_transform_from_origin(xyz: [f64; 3], rpy: [f64; 3]) -> SpatialTransform {
    SpatialTransform::from_rotation_translation(rpy_to_mat3(rpy), xyz)
}

/// Build a [`SpatialInertia`] (about the COM) from URDF inertial data.
///
/// `UrdfInertial::inertia_diag` / `inertia_off` are already expressed in the
/// link frame (the parser folds the `<inertial><origin rpy>` rotation into
/// them), so they are passed to [`SpatialInertia::from_com`] unchanged.
fn build_spatial_inertia(inert: &UrdfInertial) -> SpatialInertia {
    let [ixx, iyy, izz] = inert.inertia_diag;
    let [ixy, ixz, iyz] = inert.inertia_off;
    let i_com = [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]];
    SpatialInertia::from_com(inert.mass, inert.com, i_com)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_transform_from_zero_rpy() {
        let t = spatial_transform_from_origin([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        // rpy(0,0,0) -> identity rotation
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for (r, (rot_row, id_row)) in t.rot.iter().zip(id.iter()).enumerate() {
            for (c, (rot_val, id_val)) in rot_row.iter().zip(id_row.iter()).enumerate() {
                assert!((rot_val - id_val).abs() < 1e-12, "rot[{r}][{c}]");
            }
        }
        assert!((t.trans[0] - 1.0).abs() < 1e-12);
        assert!((t.trans[1] - 2.0).abs() < 1e-12);
        assert!((t.trans[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn yaw_rotation_is_about_z() {
        // A 90-degree yaw rotates +x onto +y (column 0 of the rotation matrix).
        let t =
            spatial_transform_from_origin([0.0, 0.0, 0.0], [0.0, 0.0, std::f64::consts::FRAC_PI_2]);
        assert!((t.rot[0][0] - 0.0).abs() < 1e-12);
        assert!((t.rot[1][0] - 1.0).abs() < 1e-12);
        assert!((t.rot[2][2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn inertial_origin_rpy_reaches_the_articulated_body_inertia() {
        // `<inertial><origin rpy>` used to be ignored, so the body kept the
        // tensor in the (differently oriented) inertial frame.
        let urdf = r#"<robot name="r">
  <link name="base"/>
  <joint name="j" type="revolute">
    <parent link="base"/>
    <child link="arm"/>
    <axis xyz="0 0 1"/>
    <limit lower="-1.0" upper="1.0" effort="1" velocity="1"/>
  </joint>
  <link name="arm">
    <inertial>
      <origin xyz="0.1 0.2 0.3" rpy="0 0 1.5707963267948966"/>
      <mass value="2.0"/>
      <inertia ixx="1.0" iyy="2.0" izz="3.0" ixy="0.5" ixz="0.0" iyz="0.0"/>
    </inertial>
  </link>
</robot>"#;
        let robot = super::super::urdf_parser::parse_urdf(urdf).unwrap();
        let model = urdf_to_articulated(&robot, "base", [0.0, 0.0, -9.81]).unwrap();
        let arm = model
            .bodies
            .iter()
            .find(|b| b.name == "arm")
            .expect("arm body present");
        assert_eq!(arm.inertia.mass, 2.0);
        assert_eq!(arm.inertia.com, [0.1, 0.2, 0.3]);

        // Rz(90 deg) * [[1, .5, 0], [.5, 2, 0], [0, 0, 3]] * Rz(90 deg)^T
        //   = [[2, -.5, 0], [-.5, 1, 0], [0, 0, 3]]   (inertia about the COM)
        // plus the parallel-axis term m * (|c|^2 E - c c^T) with m = 2,
        // c = (0.1, 0.2, 0.3), |c|^2 = 0.14.
        let expected = [
            [2.26, -0.54, -0.06],
            [-0.54, 1.2, -0.12],
            [-0.06, -0.12, 3.1],
        ];
        for (r, (got_row, want_row)) in arm
            .inertia
            .rot_inertia_origin
            .iter()
            .zip(expected.iter())
            .enumerate()
        {
            for (c, (got, want)) in got_row.iter().zip(want_row.iter()).enumerate() {
                assert!(
                    (got - want).abs() < 1e-12,
                    "rot_inertia_origin[{r}][{c}] = {got}, expected {want}"
                );
            }
        }
    }

    #[test]
    fn spatial_inertia_from_urdf_does_not_panic() {
        // Degenerate inertia (izz = 0.0) is intentionally accepted here.
        let inert = UrdfInertial {
            mass: 2.5,
            com: [0.1, 0.0, -0.2],
            inertia_diag: [1.0, 1.0, 0.0],
            inertia_off: [0.0, 0.0, 0.0],
        };
        let _ = build_spatial_inertia(&inert);
    }
}
