//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::types::*;

#[inline]
pub(super) fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub(super) fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[inline]
pub(super) fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
#[inline]
pub(super) fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let len = length3(v);
    if len < 1e-30 {
        return [0.0; 3];
    }
    scale3(v, 1.0 / len)
}
#[inline]
pub(super) fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub(super) fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub(super) fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Compute the squared length of a 3-D vector.
#[inline]
pub(super) fn len2_3(v: [f64; 3]) -> f64 {
    dot3(v, v)
}
/// Rotate a point `p` by unit quaternion `q = [w, x, y, z]`.
pub(super) fn quat_rotate(q: [f64; 4], p: [f64; 3]) -> [f64; 3] {
    let [w, x, y, z] = q;
    let qv = [x, y, z];
    let t = scale3(cross3(qv, p), 2.0);
    add3(add3(p, scale3(t, w)), cross3(qv, t))
}
#[cfg(test)]
mod tests {
    use super::*;
    fn make_simple_mesh(n: usize) -> DeformableGpuMesh {
        let mut mesh = DeformableGpuMesh::new();
        for i in 0..n {
            let x = i as f64;
            mesh.add_vertex(DeformableVertex::new([x, 0.0, 0.0]));
        }
        mesh
    }
    #[test]
    fn test_mesh_add_vertex() {
        let mut mesh = DeformableGpuMesh::new();
        let idx = mesh.add_vertex(DeformableVertex::new([1.0, 2.0, 3.0]));
        assert_eq!(idx, 0);
        assert_eq!(mesh.vertices.len(), 1);
    }
    #[test]
    fn test_mesh_add_face() {
        let mut mesh = make_simple_mesh(4);
        mesh.add_face(0, 1, 2);
        assert_eq!(mesh.faces.len(), 1);
        assert_eq!(mesh.faces[0].a, 0);
    }
    #[test]
    fn test_mesh_aabb() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([0.0, 0.0, 0.0]));
        mesh.add_vertex(DeformableVertex::new([1.0, 2.0, 3.0]));
        let (lo, hi) = mesh.aabb();
        assert!((lo[0]).abs() < 1e-10);
        assert!((hi[2] - 3.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_aabb_single_vertex() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([5.0, 6.0, 7.0]));
        let (lo, hi) = mesh.aabb();
        assert!((lo[0] - 5.0).abs() < 1e-10);
        assert!((hi[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_mesh_face_normals() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([0.0, 0.0, 0.0]));
        mesh.add_vertex(DeformableVertex::new([1.0, 0.0, 0.0]));
        mesh.add_vertex(DeformableVertex::new([0.0, 1.0, 0.0]));
        mesh.add_face(0, 1, 2);
        let normals = mesh.compute_face_normals();
        assert_eq!(normals.len(), 1);
        assert!((normals[0][2] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_blend_shape_apply() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([0.0, 0.0, 0.0]));
        let bs = BlendShape::new("test", vec![[0.0, 1.0, 0.0]]);
        let idx = mesh.add_blend_shape(bs);
        mesh.blend_shapes[idx].weight = 1.0;
        mesh.apply_blend_shapes();
        assert!((mesh.vertices[0].position[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_blend_shape_zero_weight() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([0.0, 0.0, 0.0]));
        let bs = BlendShape::new("test", vec![[0.0, 1.0, 0.0]]);
        mesh.add_blend_shape(bs);
        mesh.apply_blend_shapes();
        assert!((mesh.vertices[0].position[1]).abs() < 1e-10);
    }
    #[test]
    fn test_joint_transform_identity() {
        let jt = JointTransform::identity();
        let p = [1.0, 2.0, 3.0];
        let tp = jt.transform_point(p);
        for k in 0..3 {
            assert!((tp[k] - p[k]).abs() < 1e-10, "k={k}");
        }
    }
    #[test]
    fn test_joint_transform_translation() {
        let mut jt = JointTransform::identity();
        jt.row0[3] = 5.0;
        let tp = jt.transform_point([0.0, 0.0, 0.0]);
        assert!((tp[0] - 5.0).abs() < 1e-10);
    }
    #[test]
    fn test_dualquat_identity_transform() {
        let dq = DualQuat::identity();
        let p = [1.0, 2.0, 3.0];
        let tp = dq.transform_point(p);
        for k in 0..3 {
            assert!((tp[k] - p[k]).abs() < 1e-9, "k={k}, got {}", tp[k]);
        }
    }
    #[test]
    fn test_dualquat_pure_translation() {
        let dq = DualQuat::from_translation([3.0, 0.0, 0.0]);
        let p = [0.0, 0.0, 0.0];
        let tp = dq.transform_point(p);
        assert!((tp[0] - 3.0).abs() < 1e-9, "got {}", tp[0]);
    }
    #[test]
    fn test_dualquat_normalize_identity() {
        let dq = DualQuat::identity();
        let dn = dq.normalize();
        assert!((dn.real[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_skinning_lbs_identity() {
        let kernel = SkinningKernel::new(SkinningMode::LinearBlend, 1);
        let v = DeformableVertex::new([1.0, 2.0, 3.0]);
        let result = kernel.apply_lbs(&v);
        for (result_val, rest_val) in result.iter().zip(v.rest_position.iter()) {
            assert!((result_val - rest_val).abs() < 1e-10);
        }
    }
    #[test]
    fn test_skinning_dqs_identity() {
        let kernel = SkinningKernel::new(SkinningMode::DualQuaternion, 1);
        let v = DeformableVertex::new([1.0, 2.0, 3.0]);
        let result = kernel.apply_dqs(&v);
        for (k, (result_val, rest_val)) in result.iter().zip(v.rest_position.iter()).enumerate() {
            assert!((result_val - rest_val).abs() < 1e-9, "k={k}");
        }
    }
    #[test]
    fn test_skinning_dispatch_lbs() {
        let mut kernel = SkinningKernel::new(SkinningMode::LinearBlend, 1);
        kernel.joints[0].row1[3] = 5.0;
        let mut mesh = make_simple_mesh(3);
        kernel.dispatch(&mut mesh);
        for v in &mesh.vertices {
            assert!((v.position[1] - 5.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_skinning_dispatch_dqs() {
        let mut kernel = SkinningKernel::new(SkinningMode::DualQuaternion, 1);
        kernel.dq_joints[0] = DualQuat::from_translation([0.0, 5.0, 0.0]);
        let mut mesh = make_simple_mesh(3);
        kernel.dispatch(&mut mesh);
        for v in &mesh.vertices {
            assert!((v.position[1] - 5.0).abs() < 1e-9, "got {}", v.position[1]);
        }
    }
    #[test]
    fn test_xpbd_no_constraint() {
        let kernel = XPBDGpuKernel::new(4);
        let mut pos = vec![[0.0f64; 3], [1.0, 0.0, 0.0]];
        let inv_m = vec![1.0, 1.0];
        kernel.solve_step(&mut pos, &inv_m, 0.01);
        assert!((pos[0][0]).abs() < 1e-10);
        assert!((pos[1][0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_xpbd_distance_constraint_satisfied() {
        let mut kernel = XPBDGpuKernel::new(20);
        kernel.add_constraint(XpbdDistConstraint::new(0, 1, 1.0, 0.0));
        let mut pos = vec![[0.0f64; 3], [2.0, 0.0, 0.0]];
        let inv_m = vec![1.0, 1.0];
        kernel.solve_step(&mut pos, &inv_m, 0.01);
        let dist = length3(sub3(pos[1], pos[0]));
        assert!((dist - 1.0).abs() < 1e-6, "dist = {dist}");
    }
    #[test]
    fn test_xpbd_fixed_particle() {
        let mut kernel = XPBDGpuKernel::new(10);
        kernel.add_constraint(XpbdDistConstraint::new(0, 1, 1.0, 0.0));
        let mut pos = vec![[0.0f64; 3], [5.0, 0.0, 0.0]];
        let inv_m = vec![0.0, 1.0];
        kernel.solve_step(&mut pos, &inv_m, 0.01);
        assert!((pos[0][0]).abs() < 1e-10);
        let dist = length3(sub3(pos[1], pos[0]));
        assert!((dist - 1.0).abs() < 1e-6, "dist = {dist}");
    }
    #[test]
    fn test_xpbd_integrate_gravity() {
        let kernel = XPBDGpuKernel::new(1);
        let mut pos = vec![[0.0f64; 3]];
        let mut vel = vec![[0.0f64; 3]];
        let inv_m = vec![1.0];
        let gravity = [0.0, -9.81, 0.0];
        let dt = 0.1;
        kernel.integrate(&mut pos, &mut vel, &inv_m, gravity, dt);
        assert!((vel[0][1] - (-9.81 * dt)).abs() < 1e-10);
        assert!((pos[0][1] - (-9.81 * dt * dt)).abs() < 1e-10);
    }
    #[test]
    fn test_xpbd_integrate_fixed_particle() {
        let kernel = XPBDGpuKernel::new(1);
        let mut pos = vec![[0.0f64; 3]];
        let mut vel = vec![[0.0f64; 3]];
        let inv_m = vec![0.0];
        kernel.integrate(&mut pos, &mut vel, &inv_m, [0.0, -9.81, 0.0], 0.1);
        assert!((pos[0][1]).abs() < 1e-10);
    }
    #[test]
    fn test_fem_forces_at_rest() {
        let mut fem = FEMGpuKernel::new(4);
        let elem = TetElement::new([0, 1, 2, 3], 1.0 / 6.0, 1e6, 0.3);
        fem.add_element(elem);
        let rest = vec![
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let forces = fem.compute_forces(&rest, &rest);
        for f in &forces {
            let mag = length3(*f);
            assert!(mag < 1e-10, "mag = {mag}");
        }
    }
    #[test]
    fn test_fem_forces_nonzero_when_stretched() {
        let mut fem = FEMGpuKernel::new(4);
        let elem = TetElement::new([0, 1, 2, 3], 1.0 / 6.0, 1e6, 0.3);
        fem.add_element(elem);
        let rest = vec![
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut deformed = rest.clone();
        deformed[1] = [2.0, 0.0, 0.0];
        let forces = fem.compute_forces(&deformed, &rest);
        let total_mag: f64 = forces.iter().map(|f| length3(*f)).sum();
        assert!(total_mag > 1e-6, "expected nonzero forces");
    }
    #[test]
    fn test_fem_strain_energy_zero_at_rest() {
        let mut fem = FEMGpuKernel::new(4);
        let elem = TetElement::new([0, 1, 2, 3], 1.0, 1e5, 0.3);
        fem.add_element(elem);
        let rest = vec![
            [0.0f64, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let energy = fem.compute_strain_energy(&rest, &rest);
        assert!(energy < 1e-10, "energy = {energy}");
    }
    #[test]
    fn test_fem_lame_params() {
        let elem = TetElement::new([0, 1, 2, 3], 1.0, 1e6, 0.25);
        let (lambda, mu) = elem.lame_params();
        assert!(lambda > 0.0);
        assert!(mu > 0.0);
    }
    #[test]
    fn test_collision_signed_distance_above() {
        let cr = CollisionResponseGpu::new(CollisionResponseParams::default());
        let sd = cr.signed_distance_to_triangle(
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        assert!(sd > 0.0, "sd = {sd}");
    }
    #[test]
    fn test_collision_signed_distance_below() {
        let cr = CollisionResponseGpu::new(CollisionResponseParams::default());
        let sd = cr.signed_distance_to_triangle(
            [0.0, -0.5, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
        );
        assert!(sd < 0.0, "sd = {sd}");
    }
    #[test]
    fn test_collision_penalty_forces_below_threshold() {
        let cr = CollisionResponseGpu::new(CollisionResponseParams::new(1e4, 10.0, 1.0));
        let positions = vec![[0.5, -0.001, 0.5f64]];
        let velocities = vec![[0.0f64; 3]];
        let faces = vec![([0.0, 0.0, 0.0f64], [1.0, 0.0, 0.0f64], [0.0, 0.0, 1.0f64])];
        let mut forces = vec![[0.0f64; 3]];
        cr.apply_penalty_forces(&positions, &velocities, &faces, &mut forces);
        assert!(forces[0][1] > 0.0, "expected upward penalty force");
    }
    #[test]
    fn test_collision_penalty_forces_above_plane() {
        let cr = CollisionResponseGpu::new(CollisionResponseParams::default());
        let positions = vec![[0.5, 1.0, 0.5f64]];
        let velocities = vec![[0.0f64; 3]];
        let faces = vec![([0.0, 0.0, 0.0f64], [1.0, 0.0, 0.0f64], [0.0, 0.0, 1.0f64])];
        let mut forces = vec![[0.0f64; 3]];
        cr.apply_penalty_forces(&positions, &velocities, &faces, &mut forces);
        let mag = length3(forces[0]);
        assert!(mag < 1e-10, "expected no force, got mag = {mag}");
    }
    #[test]
    fn test_collision_sphere_penalty_inside() {
        let cr = CollisionResponseGpu::new(CollisionResponseParams::new(1e4, 0.0, 0.01));
        let positions = vec![[0.0, 0.5, 0.0f64]];
        let velocities = vec![[0.0f64; 3]];
        let mut forces = vec![[0.0f64; 3]];
        cr.apply_sphere_penalty(&positions, &velocities, [0.0, 0.0, 0.0], 1.0, &mut forces);
        let mag = length3(forces[0]);
        assert!(mag > 0.0, "expected sphere penalty force");
    }
    #[test]
    fn test_collision_sphere_penalty_outside() {
        let cr = CollisionResponseGpu::new(CollisionResponseParams::default());
        let positions = vec![[0.0, 2.0, 0.0f64]];
        let velocities = vec![[0.0f64; 3]];
        let mut forces = vec![[0.0f64; 3]];
        cr.apply_sphere_penalty(&positions, &velocities, [0.0, 0.0, 0.0], 1.0, &mut forces);
        let mag = length3(forces[0]);
        assert!(mag < 1e-10, "expected no sphere force");
    }
    #[test]
    fn test_pipeline_creation() {
        let mesh = make_simple_mesh(4);
        let config = PhysicsPipelineConfig::default();
        let pipeline = PhysicsGpuPipeline::new(config, mesh);
        assert_eq!(pipeline.mesh.vertices.len(), 4);
        assert!((pipeline.time).abs() < 1e-10);
    }
    #[test]
    fn test_pipeline_step_advances_time() {
        let mesh = make_simple_mesh(4);
        let config = PhysicsPipelineConfig::default();
        let mut pipeline = PhysicsGpuPipeline::new(config, mesh);
        pipeline.step();
        assert!(pipeline.time > 0.0);
    }
    #[test]
    fn test_pipeline_gravity_falls() {
        let mut mesh = DeformableGpuMesh::new();
        let v = DeformableVertex::new([0.0, 10.0, 0.0]);
        mesh.add_vertex(v);
        let config = PhysicsPipelineConfig {
            dt: 0.1,
            substeps: 1,
            xpbd_enabled: true,
            fem_enabled: false,
            collision_enabled: false,
            ..Default::default()
        };
        let mut pipeline = PhysicsGpuPipeline::new(config, mesh);
        pipeline.step();
        assert!(
            pipeline.mesh.vertices[0].position[1] < 10.0,
            "y = {}",
            pipeline.mesh.vertices[0].position[1]
        );
    }
    #[test]
    fn test_pipeline_centre_of_mass() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([0.0, 0.0, 0.0]));
        mesh.add_vertex(DeformableVertex::new([2.0, 0.0, 0.0]));
        let config = PhysicsPipelineConfig::default();
        let pipeline = PhysicsGpuPipeline::new(config, mesh);
        let com = pipeline.centre_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_pipeline_kinetic_energy_zero_at_rest() {
        let mesh = make_simple_mesh(4);
        let config = PhysicsPipelineConfig::default();
        let pipeline = PhysicsGpuPipeline::new(config, mesh);
        assert!((pipeline.kinetic_energy()).abs() < 1e-10);
    }
    #[test]
    fn test_pipeline_xpbd_constraint_spring() {
        let mut mesh = DeformableGpuMesh::new();
        let mut v0 = DeformableVertex::new([0.0, 0.0, 0.0]);
        v0.mass = 1e30;
        mesh.add_vertex(v0);
        let v1 = DeformableVertex::new([5.0, 0.0, 0.0]);
        mesh.add_vertex(v1);
        let config = PhysicsPipelineConfig {
            gravity: [0.0; 3],
            dt: 0.01,
            substeps: 10,
            collision_enabled: false,
            fem_enabled: false,
            xpbd_enabled: true,
        };
        let mut pipeline = PhysicsGpuPipeline::new(config, mesh);
        pipeline
            .xpbd
            .add_constraint(XpbdDistConstraint::new(0, 1, 1.0, 0.0));
        pipeline.step();
        let dx = pipeline.mesh.vertices[1].position[0];
        assert!(dx < 5.0, "particle should move, dx = {dx}");
    }
    #[test]
    fn test_xpbd_multiple_constraints() {
        let mut kernel = XPBDGpuKernel::new(50);
        kernel.add_constraint(XpbdDistConstraint::new(0, 1, 1.0, 0.0));
        kernel.add_constraint(XpbdDistConstraint::new(1, 2, 1.0, 0.0));
        let mut pos = vec![[0.0f64; 3], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
        let inv_m = vec![0.0, 1.0, 1.0];
        kernel.solve_step(&mut pos, &inv_m, 0.01);
        let d01 = length3(sub3(pos[1], pos[0]));
        let d12 = length3(sub3(pos[2], pos[1]));
        assert!((d01 - 1.0).abs() < 1e-5, "d01 = {d01}");
        assert!((d12 - 1.0).abs() < 1e-5, "d12 = {d12}");
    }
    #[test]
    fn test_blend_shape_partial_weight() {
        let mut mesh = DeformableGpuMesh::new();
        mesh.add_vertex(DeformableVertex::new([0.0, 0.0, 0.0]));
        let bs = BlendShape::new("half", vec![[0.0, 2.0, 0.0]]);
        let idx = mesh.add_blend_shape(bs);
        mesh.blend_shapes[idx].weight = 0.5;
        mesh.apply_blend_shapes();
        assert!((mesh.vertices[0].position[1] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_fem_add_element_updates_num_nodes() {
        let mut fem = FEMGpuKernel::new(0);
        fem.add_element(TetElement::new([0, 1, 2, 10], 1.0, 1e5, 0.3));
        assert!(fem.num_nodes >= 11);
    }
    #[test]
    fn test_collision_params_default() {
        let p = CollisionResponseParams::default();
        assert!(p.penalty_stiffness > 0.0);
        assert!(p.damping >= 0.0);
        assert!(p.contact_threshold > 0.0);
    }
    #[test]
    fn test_deformable_vertex_default_weights() {
        let v = DeformableVertex::new([1.0, 2.0, 3.0]);
        assert!((v.joint_weights[0] - 1.0).abs() < 1e-10);
        assert!((v.joint_weights[1]).abs() < 1e-10);
    }
    #[test]
    fn test_quat_rotate_identity() {
        let q = [1.0, 0.0, 0.0, 0.0];
        let p = [3.0, 4.0, 5.0];
        let r = quat_rotate(q, p);
        for k in 0..3 {
            assert!((r[k] - p[k]).abs() < 1e-9);
        }
    }
    #[test]
    fn test_skinning_mode_eq() {
        assert_eq!(SkinningMode::LinearBlend, SkinningMode::LinearBlend);
        assert_ne!(SkinningMode::LinearBlend, SkinningMode::DualQuaternion);
    }
    #[test]
    fn test_pipeline_multiple_steps() {
        let mesh = make_simple_mesh(2);
        let config = PhysicsPipelineConfig {
            gravity: [0.0, -9.81, 0.0],
            dt: 0.016,
            substeps: 2,
            collision_enabled: false,
            fem_enabled: false,
            xpbd_enabled: true,
        };
        let mut pipeline = PhysicsGpuPipeline::new(config, mesh);
        for _ in 0..10 {
            pipeline.step();
        }
        assert!((pipeline.time - 0.16).abs() < 1e-6);
    }
    #[test]
    fn test_gpu_node_new() {
        let node = DeformableGpuNode::new([1.0, 2.0, 3.0], 2.5_f64);
        assert!((node.mass - 2.5_f64).abs() < 1e-10);
        assert_eq!(node.vel, [0.0; 3]);
        assert_eq!(node.force, [0.0; 3]);
    }
    #[test]
    fn test_gpu_node_apply_force() {
        let mut node = DeformableGpuNode::new([0.0; 3], 1.0_f64);
        node.apply_force([1.0, 0.0, 0.0]);
        node.apply_force([0.0, 2.0, 0.0]);
        assert!((node.force[0] - 1.0_f64).abs() < 1e-10);
        assert!((node.force[1] - 2.0_f64).abs() < 1e-10);
    }
    #[test]
    fn test_gpu_node_reset_force() {
        let mut node = DeformableGpuNode::new([0.0; 3], 1.0_f64);
        node.apply_force([5.0, 5.0, 5.0]);
        node.reset_force();
        assert_eq!(node.force, [0.0; 3]);
    }
    #[test]
    fn test_gpu_node_kinetic_energy() {
        let mut node = DeformableGpuNode::new([0.0; 3], 2.0_f64);
        node.vel = [3.0, 4.0, 0.0];
        assert!((node.kinetic_energy() - 25.0_f64).abs() < 1e-10);
    }
    #[test]
    fn test_gpu_buffer_push_len() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 1.0_f64);
        buf.push([1.0, 0.0, 0.0], [0.0; 3], 1.0_f64);
        assert_eq!(buf.len(), 2);
        assert!(!buf.is_empty());
    }
    #[test]
    fn test_gpu_buffer_from_nodes() {
        let nodes = vec![
            DeformableGpuNode::new([0.0; 3], 1.0_f64),
            DeformableGpuNode::new([1.0, 0.0, 0.0], 2.0_f64),
        ];
        let buf = DeformableGpuBuffer::from_nodes(&nodes);
        assert_eq!(buf.len(), 2);
        assert!((buf.masses[1] - 2.0_f64).abs() < 1e-10);
        assert!((buf.inv_masses[1] - 0.5_f64).abs() < 1e-10);
    }
    #[test]
    fn test_gpu_buffer_reset_forces() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 1.0_f64);
        buf.forces[0] = [5.0, 5.0, 5.0];
        buf.reset_forces();
        assert_eq!(buf.forces[0], [0.0; 3]);
    }
    #[test]
    fn test_gpu_buffer_apply_gravity() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 2.0_f64);
        buf.apply_gravity([0.0, -9.81_f64, 0.0]);
        assert!((buf.forces[0][1] - (-19.62_f64)).abs() < 1e-5);
    }
    #[test]
    fn test_gpu_buffer_integrate() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 1.0_f64);
        buf.forces[0] = [0.0, -9.81_f64, 0.0];
        buf.integrate(0.1_f64);
        assert!((buf.velocities[0][1] - (-0.981_f64)).abs() < 1e-6);
    }
    #[test]
    fn test_gpu_buffer_total_ke() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [1.0, 0.0, 0.0], 2.0_f64);
        assert!((buf.total_kinetic_energy() - 1.0_f64).abs() < 1e-10);
    }
    #[test]
    fn test_gpu_buffer_flatten_positions() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([1.0, 2.0, 3.0], [0.0; 3], 1.0_f64);
        let flat = buf.flatten_positions();
        assert_eq!(flat, vec![1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_gpu_buffer_fixed_particle_inv_mass() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 0.0_f64);
        assert!((buf.inv_masses[0]).abs() < 1e-15);
    }
    #[test]
    fn test_xpbd_solver_distance_constraint() {
        let mut solver = XPBDGpuSolver::new(20);
        solver.add_distance(0, 1, 1.0_f64, 0.0_f64);
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 1.0_f64);
        buf.push([3.0, 0.0, 0.0], [0.0; 3], 1.0_f64);
        solver.solve(&mut buf, 0.01_f64);
        let dist = length3(sub3(buf.positions[1], buf.positions[0]));
        assert!((dist - 1.0_f64).abs() < 0.01_f64, "dist = {dist}");
    }
    #[test]
    fn test_xpbd_solver_no_constraint_unchanged() {
        let solver = XPBDGpuSolver::new(5);
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 1.0_f64);
        buf.push([2.0, 0.0, 0.0], [0.0; 3], 1.0_f64);
        solver.solve(&mut buf, 0.01_f64);
        assert!((buf.positions[0][0]).abs() < 1e-10);
        assert!((buf.positions[1][0] - 2.0_f64).abs() < 1e-10);
    }
    #[test]
    fn test_xpbd_solver_shape_matching() {
        let mut solver = XPBDGpuSolver::new(20);
        let rest = vec![[0.0; 3], [1.0, 0.0, 0.0]];
        solver.add_shape(ShapeMatchConstraint::new(vec![0, 1], rest, 1.0_f64));
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 1.0_f64);
        buf.push([3.0, 0.0, 0.0], [0.0; 3], 1.0_f64);
        solver.solve(&mut buf, 0.01_f64);
        let gap = (buf.positions[1][0] - buf.positions[0][0]).abs();
        assert!(gap < 3.0_f64, "gap = {gap}");
    }
    #[test]
    fn test_fem_solver_rest_zero_energy() {
        let rest = [
            [0.0, 0.0, 0.0f64],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let elem = CorotationalElement::from_rest([0, 1, 2, 3], &rest, 1e6_f64, 0.3_f64);
        let mut solver = FEMGpuSolver::new(4);
        solver.add_element(elem);
        let mut buf = DeformableGpuBuffer::new();
        for &p in &rest {
            buf.push(p, [0.0; 3], 1.0_f64);
        }
        let energy = solver.strain_energy(&buf);
        assert!(energy < 1e-5_f64, "energy = {energy}");
    }
    #[test]
    fn test_fem_solver_assemble_forces_at_rest() {
        let rest = [
            [0.0, 0.0, 0.0f64],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let elem = CorotationalElement::from_rest([0, 1, 2, 3], &rest, 1e5_f64, 0.3_f64);
        let mut solver = FEMGpuSolver::new(4);
        solver.add_element(elem);
        let mut buf = DeformableGpuBuffer::new();
        for &p in &rest {
            buf.push(p, [0.0; 3], 1.0_f64);
        }
        buf.reset_forces();
        solver.assemble_forces(&mut buf);
        for i in 0..4 {
            let mag = length3(buf.forces[i]);
            assert!(mag < 1e-5_f64, "node {i} force mag = {mag}");
        }
    }
    #[test]
    fn test_fem_corotational_element_lame() {
        let rest = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let elem = CorotationalElement::from_rest([0, 1, 2, 3], &rest, 1e6_f64, 0.25_f64);
        let (lambda, mu) = elem.lame();
        assert!(lambda > 0.0_f64);
        assert!(mu > 0.0_f64);
    }
    #[test]
    fn test_fem_solver_num_nodes_update() {
        let rest = [[0.0; 3]; 4];
        let elem = CorotationalElement::from_rest([0, 1, 9, 3], &rest, 1e5_f64, 0.3_f64);
        let mut solver = FEMGpuSolver::new(0);
        solver.add_element(elem);
        assert!(solver.num_nodes >= 10);
    }
    #[test]
    fn test_gpu_collision_plane_detect() {
        let cr = GPUCollisionResponse::default_params();
        let positions = vec![[0.0, -0.01_f64, 0.0]];
        let contacts = cr.detect_plane_contacts(&positions, [0.0, 1.0, 0.0], 0.0_f64);
        assert!(!contacts.is_empty());
        assert_eq!(contacts[0].particle_idx, 0);
    }
    #[test]
    fn test_gpu_collision_plane_no_contact_above() {
        let cr = GPUCollisionResponse::default_params();
        let positions = vec![[0.0, 1.0_f64, 0.0]];
        let contacts = cr.detect_plane_contacts(&positions, [0.0, 1.0, 0.0], 0.0_f64);
        assert!(contacts.is_empty());
    }
    #[test]
    fn test_gpu_collision_sphere_detect() {
        let cr = GPUCollisionResponse::default_params();
        let positions = vec![[0.5_f64, 0.0, 0.0]];
        let contacts = cr.detect_sphere_contacts(&positions, [0.0; 3], 1.0_f64);
        assert!(!contacts.is_empty());
    }
    #[test]
    fn test_gpu_collision_resolve_positions() {
        let cr = GPUCollisionResponse::default_params();
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0, -0.05_f64, 0.0], [0.0; 3], 1.0_f64);
        let contacts = cr.detect_plane_contacts(&buf.positions, [0.0, 1.0, 0.0], 0.0_f64);
        cr.resolve_positions(&mut buf, &contacts);
        assert!(
            buf.positions[0][1] >= 0.0_f64 - 1e-10,
            "y = {}",
            buf.positions[0][1]
        );
    }
    #[test]
    fn test_gpu_collision_resolve_velocities() {
        let cr = GPUCollisionResponse::new(1e4_f64, 10.0_f64, 0.01_f64, 0.5_f64);
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0, -0.001_f64, 0.0], [0.0, -2.0_f64, 0.0], 1.0_f64);
        let contacts = cr.detect_plane_contacts(&buf.positions, [0.0, 1.0, 0.0], 0.0_f64);
        cr.resolve_velocities(&mut buf, &contacts);
        assert!(
            buf.velocities[0][1] > 0.0_f64,
            "vel = {}",
            buf.velocities[0][1]
        );
    }
    #[test]
    fn test_deformable_pipeline_creation() {
        let pipeline = DeformableGpuPipeline::new([0.0, -9.81_f64, 0.0], 1.0 / 60.0_f64, 4);
        assert_eq!(pipeline.buf.len(), 0);
        assert!((pipeline.time).abs() < 1e-15);
    }
    #[test]
    fn test_deformable_pipeline_add_particle() {
        let mut pipeline = DeformableGpuPipeline::new([0.0; 3], 0.016_f64, 2);
        let idx = pipeline.add_particle([0.0, 5.0_f64, 0.0], 1.0_f64);
        assert_eq!(idx, 0);
        assert_eq!(pipeline.buf.len(), 1);
    }
    #[test]
    fn test_deformable_pipeline_step_advances_time() {
        let mut pipeline = DeformableGpuPipeline::new([0.0, -9.81_f64, 0.0], 0.016_f64, 2);
        pipeline.add_particle([0.0, 10.0_f64, 0.0], 1.0_f64);
        pipeline.step();
        assert!(pipeline.time > 0.0_f64);
    }
    #[test]
    fn test_deformable_pipeline_gravity_falls() {
        let mut pipeline = DeformableGpuPipeline::new([0.0, -9.81_f64, 0.0], 0.1_f64, 1);
        pipeline.collision_enabled = false;
        pipeline.xpbd_enabled = false;
        pipeline.fem_enabled = false;
        pipeline.add_particle([0.0, 10.0_f64, 0.0], 1.0_f64);
        pipeline.step();
        assert!(
            pipeline.buf.positions[0][1] < 10.0_f64,
            "y = {}",
            pipeline.buf.positions[0][1]
        );
    }
    #[test]
    fn test_deformable_pipeline_kinetic_energy() {
        let mut pipeline = DeformableGpuPipeline::new([0.0; 3], 0.016_f64, 1);
        pipeline.add_particle([0.0; 3], 1.0_f64);
        assert!((pipeline.kinetic_energy()).abs() < 1e-15);
    }
    #[test]
    fn test_deformable_pipeline_centre_of_mass() {
        let mut pipeline = DeformableGpuPipeline::new([0.0; 3], 0.016_f64, 1);
        pipeline.add_particle([0.0; 3], 1.0_f64);
        pipeline.add_particle([2.0, 0.0, 0.0], 1.0_f64);
        let com = pipeline.centre_of_mass();
        assert!((com[0] - 1.0_f64).abs() < 1e-10, "com.x = {}", com[0]);
    }
    #[test]
    fn test_deformable_pipeline_empty_com() {
        let pipeline = DeformableGpuPipeline::new([0.0; 3], 0.016_f64, 1);
        let com = pipeline.centre_of_mass();
        assert_eq!(com, [0.0; 3]);
    }
    #[test]
    fn test_xpbd_solver_with_buffer_distance() {
        let mut solver = XPBDGpuSolver::new(30);
        solver.add_distance(0, 1, 1.0_f64, 0.0_f64);
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 0.0_f64);
        buf.push([5.0, 0.0, 0.0], [0.0; 3], 1.0_f64);
        solver.solve(&mut buf, 0.01_f64);
        let dist = length3(sub3(buf.positions[1], buf.positions[0]));
        assert!((dist - 1.0_f64).abs() < 0.01_f64, "dist = {dist}");
    }
    #[test]
    fn test_gpu_buffer_fixed_particle_doesnt_move() {
        let mut buf = DeformableGpuBuffer::new();
        buf.push([0.0; 3], [0.0; 3], 0.0_f64);
        buf.forces[0] = [100.0, 100.0, 100.0];
        buf.integrate(0.1_f64);
        assert_eq!(buf.positions[0], [0.0; 3]);
    }
}
