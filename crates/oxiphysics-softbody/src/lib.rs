// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soft body simulation for the OxiPhysics engine.
//!
//! Provides PBD / XPBD based soft-body simulation including:
//!
//! - Distance, bending, volume, and collision constraints
//! - Cloth simulation with wind forces
//! - Rope / hair chains
//! - Corotational FEM soft bodies
#![warn(missing_docs)]
#![allow(ambiguous_glob_reexports)]
#![allow(dead_code)]

mod error;
pub use error::*;

pub mod aerodynamics;
pub mod cloth;
pub mod constraint;
pub mod crack_propagation;
pub mod fem_soft;
pub mod fracture;
pub mod fracture_dynamics;
pub mod inflatable;
pub mod particle;
pub mod pbd_system;
pub mod rope;
pub mod shape_matching;
pub mod solver;
pub mod volume;
pub mod xpbd;

pub use aerodynamics::AerodynamicsModel;
pub use cloth::*;
pub use constraint::{
    BendingConstraint, CollisionConstraint, DistanceConstraint, SoftConstraint, VolumeConstraint,
};
pub use fem_soft::{CorotationalElement, FemSoftBody};
pub use fracture_dynamics::*;
pub use inflatable::*;
pub use particle::{SoftBody, SoftParticle};
pub use rope::{HairStrand, HairSystem, Rope, RopeSolver, XpbdRope};
pub use shape_matching::ShapeMatching;
pub use solver::XpbdSolver;

/// Trait for soft body solvers.
pub trait SoftBodySolver {
    /// Initialize this component.
    fn init(&mut self);
}

#[cfg(test)]
mod tests {
    use super::*;
    use oxiphysics_core::math::Mat3;
    use oxiphysics_core::math::{Real, Vec3};

    const EPS: Real = 1e-6;

    // 1. Distance constraint maintains rest length.
    #[test]
    fn test_distance_constraint_maintains_rest_length() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(2.0, 0.0, 0.0), 1.0),
        ];
        let rest = 1.0;
        let mut dc = DistanceConstraint::new(0, 1, rest, 0.0);
        // Project many times to converge.
        for _ in 0..100 {
            dc.reset_lambda();
            dc.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            (dist - rest).abs() < 0.01,
            "Distance {dist} should be close to rest length {rest}"
        );
    }

    // 2. Two particles connected by distance constraint converge.
    #[test]
    fn test_two_particles_converge() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(5.0, 0.0, 0.0), 1.0),
        ];
        let rest = 1.0;
        let mut dc = DistanceConstraint::new(0, 1, rest, 0.0);
        for _ in 0..200 {
            dc.reset_lambda();
            dc.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            (dist - rest).abs() < 0.05,
            "Particles should converge to rest length"
        );
    }

    // 3. Volume constraint preserves volume.
    #[test]
    fn test_volume_constraint_preserves_volume() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 1.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let rest_vol = VolumeConstraint::compute_tet_volume(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        // Perturb vertex 3.
        particles[3].position = Vec3::new(0.0, 0.0, 1.5);
        let mut vc = VolumeConstraint::new([0, 1, 2, 3], rest_vol, 0.0);
        for _ in 0..100 {
            vc.reset_lambda();
            vc.project(&mut particles, 1.0 / 60.0);
        }
        let new_vol = VolumeConstraint::compute_tet_volume(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        assert!(
            (new_vol - rest_vol).abs() < 0.01,
            "Volume {new_vol} should be close to rest volume {rest_vol}"
        );
    }

    // 4. XPBD solver step doesn't crash.
    #[test]
    fn test_xpbd_solver_runs() {
        let mut body = SoftBody::from_particles(vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
        ]);
        body.apply_force(&Vec3::new(0.0, -9.81, 0.0));
        let mut constraints: Vec<Box<dyn SoftConstraint>> = vec![Box::new(
            DistanceConstraint::from_particles(0, 1, &body.particles, 0.0),
        )];
        let mut solver = XpbdSolver::new(5);
        solver.solve(&mut body, &mut constraints, 1.0 / 60.0);
        // Just verify no panic and positions have changed.
        let y0 = body.particles[0].position.y;
        assert!(y0 < 0.0, "Particle should have fallen under gravity");
    }

    // 5. Cloth mesh creation has correct topology.
    #[test]
    fn test_cloth_mesh_topology() {
        let cloth = XpbdClothMesh::new(4, 3, 3.0, 2.0, 1.0, 0.001);
        assert_eq!(cloth.num_particles(), 12); // 4 * 3
        // Triangles: (nx-1)*(ny-1)*2 = 3*2*2 = 12
        assert_eq!(cloth.num_triangles(), 12);
        assert!(!cloth.distance_constraints.is_empty());
        assert!(!cloth.bending_constraints.is_empty());
    }

    // 6. Wind force applies in correct direction.
    #[test]
    fn test_wind_force_direction() {
        // A triangle lying in the XZ plane, normal pointing up (+Y).
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let tris = vec![[0, 1, 2]];
        let wind = Vec3::new(0.0, -1.0, 0.0); // Wind blowing downward.
        apply_wind(&mut particles, &wind, &tris);
        // The normal is (0, -1, 0) (depends on winding).  Either way, force
        // should have a Y component.
        let total_fy: Real = particles.iter().map(|p| p.external_force.y).sum();
        assert!(
            total_fy.abs() > EPS,
            "Wind should produce a Y-direction force, got {total_fy}"
        );
    }

    // 7. Rope creation has correct number of particles/constraints.
    #[test]
    fn test_rope_creation() {
        let rope = XpbdRope::new(
            Vec3::new(0.0, 10.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
            10,
            0.5,
            0.0,
        );
        assert_eq!(rope.num_particles(), 11);
        assert_eq!(rope.num_segments(), 10);
        // First particle should be static.
        assert!(rope.body.particles[0].is_static());
    }

    // 8. Corotational element identity gives zero force.
    #[test]
    fn test_corotational_identity_zero_force() {
        let particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 1.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let elem = CorotationalElement::new([0, 1, 2, 3], &particles, 1000.0, 0.3);
        let forces = elem.compute_forces(&particles);
        for (k, f) in forces.iter().enumerate() {
            assert!(
                f.norm() < 1e-8,
                "Force on vertex {k} should be zero at rest, got {f:?}"
            );
        }
    }

    // 9. Bending constraint projection.
    #[test]
    fn test_bending_constraint_projection() {
        // Four particles: shared edge (0,1), wing vertices (2,3).
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.5, 1.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.5, -1.0, 0.0), 1.0),
        ];
        let mut bc = BendingConstraint::from_particles([0, 1, 2, 3], &particles, 0.0);
        // Perturb wing vertex.
        particles[2].position = Vec3::new(0.5, 1.0, 0.5);
        let old_pos = particles[2].position;
        for _ in 0..50 {
            bc.reset_lambda();
            bc.project(&mut particles, 1.0 / 60.0);
        }
        // The wing vertex should have moved.
        let moved = (particles[2].position - old_pos).norm();
        assert!(
            moved > 1e-8,
            "Bending constraint should move perturbed wing vertex"
        );
    }

    // 10. Particle free-fall under gravity.
    #[test]
    fn test_particle_free_fall() {
        let mut body =
            SoftBody::from_particles(vec![SoftParticle::new(Vec3::new(0.0, 10.0, 0.0), 1.0)]);
        body.apply_force(&Vec3::new(0.0, -9.81, 0.0));
        let mut constraints: Vec<Box<dyn SoftConstraint>> = Vec::new();
        let mut solver = XpbdSolver::new(1);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            solver.solve(&mut body, &mut constraints, dt);
        }
        // After ~1 second of gravity, y should be significantly lower than 10.
        let y = body.particles[0].position.y;
        assert!(
            y < 9.0,
            "Particle should have fallen under gravity, y = {y}"
        );
    }

    // 11. Collision constraint prevents penetration.
    #[test]
    fn test_collision_constraint_plane() {
        let mut particles = vec![SoftParticle::new(Vec3::new(0.0, -1.0, 0.0), 1.0)];
        let mut cc =
            CollisionConstraint::new(0, Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 0.0));
        cc.project(&mut particles, 1.0 / 60.0);
        assert!(
            particles[0].position.y >= -EPS,
            "Particle should be pushed above the plane"
        );
    }

    // 12. FEM soft body step runs.
    #[test]
    fn test_fem_soft_body_step() {
        let particles = vec![
            SoftParticle::new_static(Vec3::new(0.0, 0.0, 0.0)),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 1.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let elem = CorotationalElement::new([0, 1, 2, 3], &particles, 1000.0, 0.3);
        let mut fem = FemSoftBody::new(particles, vec![elem], 0.01);
        let gravity = Vec3::new(0.0, -9.81, 0.0);
        for _ in 0..10 {
            fem.step(1.0 / 60.0, &gravity);
        }
        // Dynamic particles should have moved.
        assert!(
            fem.particles[1].position != Vec3::new(1.0, 0.0, 0.0),
            "Dynamic particle should have moved"
        );
    }

    // 13. Corotational polar decomposition extracts pure rotation correctly.
    //
    // When the deformed tetrahedron is a pure Rz(45°) rotation of the rest
    // configuration, the deformation gradient F = R (no stretch).  The
    // corotational model should therefore produce zero elastic forces because
    // R^T * F - I = R^T * R - I = 0.
    #[test]
    fn test_corotational_polar_decomp() {
        use oxiphysics_core::math::Mat3;

        // Rest tetrahedron.
        let rest_particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(1.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 1.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let elem = CorotationalElement::new([0, 1, 2, 3], &rest_particles, 1000.0, 0.3);

        // Apply a pure Rz(45°) rotation to all vertices.
        let angle: Real = std::f64::consts::FRAC_PI_4; // 45 degrees
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let rz = Mat3::new(cos_a, -sin_a, 0.0, sin_a, cos_a, 0.0, 0.0, 0.0, 1.0);
        let rotated_particles: Vec<SoftParticle> = rest_particles
            .iter()
            .map(|p| SoftParticle::new(rz * p.position, 1.0))
            .collect();

        // Forces on a purely-rotated (unstretched) element must be zero.
        let forces = elem.compute_forces(&rotated_particles);
        for (k, f) in forces.iter().enumerate() {
            assert!(
                f.norm() < 1e-8,
                "Force on vertex {k} should be zero for pure rotation, got {f:?}"
            );
        }
    }

    // 14. Cloth sags under gravity (center vertex moves downward).
    //
    // A 3x3 grid cloth is pinned at its four corners.  After 100 XPBD steps
    // with gravity applied, the center particle (index 4) must have a negative
    // Y position (sag below the initial flat XZ plane).
    #[test]
    fn test_cloth_gravity_sag() {
        let mut cloth = XpbdClothMesh::new(3, 3, 2.0, 2.0, 1.0, 1e-4);
        // Pin the four corners.
        cloth.pin(0, 0);
        cloth.pin(2, 0);
        cloth.pin(0, 2);
        cloth.pin(2, 2);

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 1.0 / 60.0;
        let mut solver = XpbdSolver::new(10);

        // Collect all constraints into a single boxed list.
        let mut constraints: Vec<Box<dyn SoftConstraint>> = cloth
            .distance_constraints
            .iter()
            .cloned()
            .map(|c| Box::new(c) as Box<dyn SoftConstraint>)
            .chain(
                cloth
                    .bending_constraints
                    .iter()
                    .cloned()
                    .map(|c| Box::new(c) as Box<dyn SoftConstraint>),
            )
            .collect();

        for _ in 0..100 {
            cloth.body.apply_force(&gravity);
            solver.solve(&mut cloth.body, &mut constraints, dt);
            cloth.body.clear_forces();
        }

        // Center particle is index nx/2 * nx + nx/2 = 1 * 3 + 1 = 4 for a 3x3 grid.
        let center_y = cloth.body.particles[4].position.y;
        assert!(
            center_y < -0.01,
            "Center particle should sag below y=0 under gravity, got y={center_y}"
        );
    }

    // -----------------------------------------------------------------------
    // Aerodynamics tests
    // -----------------------------------------------------------------------

    // A1. Zero wind velocity → zero forces on every vertex.
    #[test]
    fn test_aero_zero_wind() {
        let aero = AerodynamicsModel::new(1.225, 1.0, 0.5);
        let zero = Vec3::zeros();
        let (f0, f1, f2) = aero.triangle_force(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            zero,
            zero,
            zero,
            zero, // zero wind
        );
        assert!(
            f0.norm() < 1e-10,
            "f0 should be zero with zero wind, got {:?}",
            f0
        );
        assert!(
            f1.norm() < 1e-10,
            "f1 should be zero with zero wind, got {:?}",
            f1
        );
        assert!(
            f2.norm() < 1e-10,
            "f2 should be zero with zero wind, got {:?}",
            f2
        );
    }

    // A2. Wind perpendicular to triangle (head-on) → purely drag force.
    //
    // Triangle in XZ plane (normal along +Y).  Wind blows in +Y direction.
    // The lift term vanishes because v_rel × n̂ = 0 when v_rel ∥ n̂.
    #[test]
    fn test_aero_head_on_wind() {
        let cd = 1.5;
        let cl = 0.8;
        let aero = AerodynamicsModel::new(1.225, cd, cl);
        let zero = Vec3::zeros();
        let wind = Vec3::new(0.0, 10.0, 0.0); // perpendicular to XZ plane
        let (f0, f1, f2) = aero.triangle_force(
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            zero,
            zero,
            zero,
            wind,
        );
        let f_total = f0 + f1 + f2;
        // Force must be entirely in the Y direction (within tolerance).
        assert!(
            f_total.x.abs() < 1e-10 && f_total.z.abs() < 1e-10,
            "Head-on wind should produce only Y force, got {:?}",
            f_total
        );
        assert!(
            f_total.y > 1e-10,
            "Head-on wind should produce positive Y drag, got {:?}",
            f_total
        );
    }

    // A3. Doubling rho_air doubles the total aerodynamic force.
    #[test]
    fn test_aero_force_scales_with_rho() {
        let aero1 = AerodynamicsModel::new(1.0, 1.0, 0.5);
        let aero2 = AerodynamicsModel::new(2.0, 1.0, 0.5);
        let zero = Vec3::zeros();
        let wind = Vec3::new(5.0, 3.0, 1.0);
        let v0 = Vec3::new(0.0, 0.0, 0.0);
        let v1 = Vec3::new(1.0, 0.0, 0.0);
        let v2 = Vec3::new(0.0, 0.0, 1.0);

        let (a0, a1, a2) = aero1.triangle_force(v0, v1, v2, zero, zero, zero, wind);
        let (b0, b1, b2) = aero2.triangle_force(v0, v1, v2, zero, zero, zero, wind);

        let fa = (a0 + a1 + a2).norm();
        let fb = (b0 + b1 + b2).norm();

        assert!(
            (fb - 2.0 * fa).abs() < 1e-10,
            "Doubling rho should double force: fa={fa}, fb={fb}"
        );
    }

    // -----------------------------------------------------------------------
    // Shape matching tests
    // -----------------------------------------------------------------------

    // S1. Pure rotation: shape-matching goal positions = rotated rest positions.
    //
    // If the current positions are a pure rotation R of the rest positions
    // (with stiffness=1), the goals must equal the current positions.
    #[test]
    fn test_shape_matching_rigid_body() {
        let rest = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
        ];
        let masses = vec![1.0f64; 4];
        let sm = ShapeMatching::new(&rest, &masses, 1.0);

        // Apply a 90° rotation about Z.
        let rz90 = Mat3::new(0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0);
        let current: Vec<Vec3> = rest.iter().map(|p| rz90 * p).collect();
        let goals = sm.goal_positions(&current);

        for (i, (goal, cur)) in goals.iter().zip(current.iter()).enumerate() {
            assert!(
                (goal - cur).norm() < 1e-6,
                "Particle {i}: goal {:?} should equal current {:?}",
                goal,
                cur
            );
        }
    }

    // S2. stiffness=0 → goals equal current positions (no constraint).
    #[test]
    fn test_shape_matching_stiffness_zero() {
        let rest = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
        ];
        let masses = vec![1.0f64; 3];
        let sm = ShapeMatching::new(&rest, &masses, 0.0); // stiffness = 0

        // Arbitrary deformed positions.
        let current = vec![
            Vec3::new(3.0, 1.0, -2.0),
            Vec3::new(-2.0, 0.5, 1.0),
            Vec3::new(0.5, 4.0, 0.3),
        ];
        let goals = sm.goal_positions(&current);

        for (i, (goal, cur)) in goals.iter().zip(current.iter()).enumerate() {
            assert!(
                (goal - cur).norm() < 1e-10,
                "With stiffness=0 goal {i} {:?} should equal current {:?}",
                goal,
                cur
            );
        }
    }

    // S3. Identity: undeformed positions → goals equal current positions.
    #[test]
    fn test_shape_matching_identity() {
        let rest = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let masses = vec![1.0f64; 4];
        let sm = ShapeMatching::new(&rest, &masses, 1.0);

        // Current positions = rest positions (no deformation).
        let goals = sm.goal_positions(&rest);

        for (i, (goal, cur)) in goals.iter().zip(rest.iter()).enumerate() {
            assert!(
                (goal - cur).norm() < 1e-6,
                "Particle {i}: goal {:?} should equal rest {:?}",
                goal,
                cur
            );
        }
    }

    // 15. Distance constraint pushes two particles that are too close to rest length.
    //
    // Place two particles at the same position (distance = 0).  After many
    // XPBD iterations the distance must converge to the rest length.
    #[test]
    fn test_distance_constraint_correction() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0), // same position
        ];
        let rest = 1.5;
        // To avoid the degenerate (dist < 1e-12) guard, give them a tiny offset.
        particles[1].position = Vec3::new(1e-6, 0.0, 0.0);

        let mut dc = DistanceConstraint::new(0, 1, rest, 0.0);
        for _ in 0..200 {
            dc.reset_lambda();
            dc.project(&mut particles, 1.0 / 60.0);
        }
        let dist = (particles[0].position - particles[1].position).norm();
        assert!(
            (dist - rest).abs() < 0.1,
            "Particles too close should be pushed to rest length {rest}, got dist={dist}"
        );
    }

    // 16. Volume constraint keeps tetrahedron volume close to rest after perturbation.
    //
    // Form a unit tetrahedron, record its volume, then compress one vertex
    // inward.  After running the volume constraint the volume should return
    // to within 5 % of the original.
    #[test]
    fn test_volume_constraint_preserve() {
        let mut particles = vec![
            SoftParticle::new(Vec3::new(0.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(2.0, 0.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 2.0, 0.0), 1.0),
            SoftParticle::new(Vec3::new(0.0, 0.0, 2.0), 1.0),
        ];
        let rest_vol = VolumeConstraint::compute_tet_volume(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );

        // Compress: move vertex 3 so volume becomes roughly half.
        particles[3].position = Vec3::new(0.0, 0.0, 0.5);

        let mut vc = VolumeConstraint::new([0, 1, 2, 3], rest_vol, 0.0);
        for _ in 0..200 {
            vc.reset_lambda();
            vc.project(&mut particles, 1.0 / 60.0);
        }

        let final_vol = VolumeConstraint::compute_tet_volume(
            &particles[0].position,
            &particles[1].position,
            &particles[2].position,
            &particles[3].position,
        );
        let relative_error = ((final_vol - rest_vol) / rest_vol).abs();
        assert!(
            relative_error < 0.05,
            "Volume after constraint should be close to rest ({rest_vol}), got {final_vol} (error={relative_error:.3})"
        );
    }
}

pub mod active_matter;
pub mod active_origami;
pub mod bio_softbody;
pub mod bioinspired;
pub mod biomech_simulation;
pub mod cable_nets;
pub mod cloth_advanced;
pub mod cloth_simulation;
pub mod collision_response;
pub mod constraint_visualization;
pub mod cosserat_rods;
pub mod cosserat_softbody;
pub mod elastic_wave;
pub mod electroactive_softbody;
pub mod fluid_film_softbody;
pub mod food_physics;
pub mod granular_softbody;
pub mod growth_mechanics;
pub mod hair;
pub mod hair_fur;
pub mod hair_sim;
pub mod hair_softbody;
pub mod haptic_softbody;
pub mod liquid_sim;
pub mod material_point;
pub mod membrane_biophysics;
pub mod membrane_softbody;
pub mod metamaterial_softbody;
pub mod metamorphic_mesh;
pub mod morphogenesis_softbody;
pub mod mud_snow;
pub mod muscle_sim;
pub mod muscle_simulation;
pub mod muscle_softbody;
pub mod muscle_tendon;
pub mod neural_deform;
pub mod neural_soft;
pub mod neural_softbody;
pub mod numerical_softbody;
pub mod origami_folding;
pub mod origami_mech;
pub mod origami_mechanics;
pub mod origami_softbody;
pub mod pbd_constraints;
pub mod peridynamics;
pub mod position_based_fluids;
pub mod pressure_volume;
pub mod rigid_soft_coupling;
pub mod rod_mechanics;
pub mod rope_simulation;
pub mod sand_sim;
pub mod smart_materials;
pub mod soft_gripper;
pub mod surgery_simulation;
pub mod surgical_simulation;
pub mod tearing;
pub mod tendon_softbody;
pub mod textile_simulation;
pub mod tissue_sim;
pub mod topology_opt;
pub mod topology_optimization;
pub mod underwater_softbody;
pub mod wrinkling;
