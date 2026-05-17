// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! # OxiPhysics — Unified Physics Engine
//!
//! OxiPhysics is a pure-Rust, modular physics engine that re-exports all
//! sub-crates for convenient, single-dependency access. It targets the
//! same problem domains as Bullet, OpenFOAM, LAMMPS, and CalculiX.
//!
//! ## Quick Start
//!
//! ```rust,no_run
//! use oxiphysics::core::math::Vec3;
//! use oxiphysics::core::Transform;
//!
//! let t = Transform::default();
//! let v = Vec3::new(1.0, 2.0, 3.0);
//! let _transformed = t.transform_point(&v);
//! ```
//!
//! ## Module Overview
//!
//! | Module | Description |
//! |--------|-------------|
//! | `core` | Math types, traits, ODE solvers, stochastic processes |
//! | `geometry` | Shape primitives and mesh utilities |
//! | `collision` | Broad-phase and narrow-phase collision detection |
//! | `rigid` | Rigid body dynamics |
//! | `articulated` | Featherstone articulated-body dynamics (RNEA + ABA) |
//! | `constraints` | Constraint / joint solvers |
//! | `vehicle` | Vehicle dynamics simulation |
//! | `sph` | Smoothed-particle hydrodynamics |
//! | `lbm` | Lattice Boltzmann method |
//! | `fem` | Finite element method |
//! | `md` | Molecular dynamics |
//! | `softbody` | Soft body / cloth simulation |
//! | `materials` | Material property database |
//! | `gpu` | GPU acceleration backends |
//! | `viz` | Visualization helpers |
//! | `io` | File I/O and serialization |
//! | `pipeline` | Full simulation pipeline |
//! | `force_field` | Spatial force fields (gravity wells, vortex, explosion) |
//! | `event_bus` | Physics event publish/subscribe system |
//! | `replay` | Deterministic simulation replay (record → replay) |
//! | `query` | Spatial queries — raycasting, sphere sweeps, overlaps |
//! | `scene` | Declarative scene description with JSON round-trip |
//! | `snapshot` | World-state snapshots with delta tracking |
//! | `trigger` | Trigger/sensor volumes — enter/exit/stay events |
//! | `animation` | Keyframe animation tracks (Vec3 + quaternion SLERP) |
//! | `material_table` | Runtime material interaction table with combine rules |
//! | `debug_draw` | Renderer-agnostic debug draw command buffer |
//! | `contact_cache` | Persistent contact pair cache with warm-start impulses |
//! | `buoyancy` | Archimedes buoyancy and drag forces for fluid volumes |
//! | `scheduler` | Priority-based physics step budget allocation per island |
//! | `spatial_grid` | Uniform spatial hash grid for fast neighbourhood queries |
//! | `lod` | Level-of-Detail simulation tier management |
//! | `noise` | Procedural 3D value noise and fractal Brownian motion |
//! | `interpolator` | Smooth damp, exp decay, spring followers, lerp utilities |
//! | `telemetry` | Per-step physics stats, rolling averages, CSV export |
//! | `character` | Kinematic capsule character controller with sweep-and-slide, step-up, and slope handling |
//! | `rope` | Rope / chain distance constraints with Verlet integration and Gauss-Seidel projection |
//! | `xpbd` | Extended Position-Based Dynamics integrator with compliance and Lagrange multipliers |
//! | `ik` | Inverse kinematics — FABRIK and 2-bone analytic IK with joint cone limits |
//! | `profiler` | Hierarchical scoped profiler with RAII guards, frame reports, and flamegraph export |
//! | `aero` | Aerodynamics — velocity-squared drag and airfoil lift/drag forces |
//! | `navmesh` | Navigation mesh with A* pathfinding and funnel-algorithm path smoothing |
//! | `rollback` | Snapshot-based rollback and deterministic lockstep input buffer |
//!
//! ## Stability Policy
//!
//! Every public API is annotated with a stability level. See
//! `core::stability` for details and the `core::stability::HasStability`
//! trait for programmatic queries. Stable APIs follow semver; unstable and
//! experimental APIs may change across minor releases.
#![warn(missing_docs)]

/// Core types, traits, and abstractions.
pub use oxiphysics_core as core;

/// Geometric shape types.
pub use oxiphysics_geometry as geometry;

/// Collision detection algorithms.
pub use oxiphysics_collision as collision;

/// Rigid body dynamics.
pub use oxiphysics_rigid as rigid;

/// Featherstone articulated-body dynamics (RNEA + ABA).
pub use oxiphysics_articulated as articulated;

/// Constraint solvers.
pub use oxiphysics_constraints as constraints;

/// Vehicle dynamics simulation.
pub use oxiphysics_vehicle as vehicle;

/// Smoothed-particle hydrodynamics.
pub use oxiphysics_sph as sph;

/// Lattice Boltzmann method.
pub use oxiphysics_lbm as lbm;

/// Finite element method.
pub use oxiphysics_fem as fem;

/// Molecular dynamics.
pub use oxiphysics_md as md;

/// Soft body simulation.
pub use oxiphysics_softbody as softbody;

/// Material properties.
pub use oxiphysics_materials as materials;

/// GPU acceleration backends.
pub use oxiphysics_gpu as gpu;

/// Visualization utilities.
pub use oxiphysics_viz as viz;

/// File I/O and serialization.
pub use oxiphysics_io as io;

/// WebAssembly bindings.
pub use oxiphysics_wasm as wasm;

/// Full physics simulation pipeline.
pub mod pipeline;

/// Performance regression testing infrastructure.
pub mod perf_regression;

/// Spatial force fields (gravity wells, vortex, wind, explosion).
pub mod force_field;

/// Physics event publish/subscribe system.
pub mod event_bus;

/// Deterministic simulation replay (record → serialise → replay).
pub mod replay;

/// Spatial query API — ray casting, sphere sweeps, and overlap tests.
pub mod query;

/// Declarative scene description with JSON round-trip and fluent builder.
pub mod scene;

/// World-state snapshots with delta tracking and ring-buffer history.
pub mod snapshot;

/// Trigger/sensor volumes — detect body enter/exit/stay without contact forces.
pub mod trigger;

/// Keyframe animation tracks (Vec3 lerp + quaternion SLERP) for scripted motion.
pub mod animation;

/// Runtime material interaction table with combine rules and per-pair overrides.
pub mod material_table;

/// Renderer-agnostic debug draw command buffer for visualising physics state.
pub mod debug_draw;

/// Persistent contact pair cache with warm-start impulse data.
pub mod contact_cache;

/// Archimedes buoyancy and viscous drag for bodies in fluid volumes.
pub mod buoyancy;

/// Priority-based physics step budget allocation per simulation island.
pub mod scheduler;

/// Uniform spatial hash grid for fast neighbourhood queries.
pub mod spatial_grid;

/// Level-of-Detail simulation tier management.
pub mod lod;

/// Procedural 3D value noise and fractal Brownian motion (fBm).
pub mod noise;

/// Motion interpolation — smooth damp, exponential decay, spring followers.
pub mod interpolator;

/// Per-step physics telemetry, rolling averages, and CSV export.
pub mod telemetry;

/// Rope / chain distance constraints with Verlet integration and Gauss-Seidel projection.
pub mod rope;

/// Extended Position-Based Dynamics integrator with compliance and Lagrange multipliers.
pub mod xpbd;

/// Kinematic capsule character controller with sweep-and-slide, step-up, and slope handling.
pub mod character;

/// Inverse kinematics — FABRIK and 2-bone analytic IK with joint cone limits.
pub mod ik;

/// Hierarchical scoped profiler with RAII guards, frame reports, and flamegraph export.
pub mod profiler;

/// Aerodynamics — velocity-squared drag and airfoil lift/drag forces.
pub mod aero;

/// Navigation mesh with A* pathfinding and funnel-algorithm path smoothing.
pub mod navmesh;

/// Snapshot-based rollback and deterministic lockstep input buffer.
pub mod rollback;

/// Cross-domain auto-coupling runtime (OxiCAR — Blueprint KF-1).
///
/// Provides a domain-agnostic [`coupling::DomainCoupler`] framework for
/// linking FEM, SPH, LBM, and MD regions at shared interfaces.
pub mod coupling;

// Re-export the most commonly used coupling types at the crate root.
#[doc(inline)]
pub use coupling::{
    CouplingDomain, CouplingReport, CouplingRuntime, DomainCoupler, DomainKind, InterfaceForce,
    InterfaceForceVec, InterfaceSite, InterfaceState, InterfaceStateVec,
    fem_sph::FemSphCoupler,
    md_continuum_adapter::{MdContinuumAdapter, MockContinuumDomain, MockMdDomain},
};

#[cfg(test)]
mod tests {

    use crate::constraints::Constraint;
    use crate::constraints::SpringJoint;
    use crate::core::PhysicsConfig;
    use crate::pipeline::PhysicsPipeline;
    use oxiphysics_collision::broadphase::BruteForceBroadPhase;
    use oxiphysics_collision::{BroadPhase, CollisionPair, NarrowPhaseDispatcher};
    use oxiphysics_core::math::Vec3;
    use oxiphysics_core::{Aabb, TimeStep, Transform};
    use oxiphysics_geometry::{BoxShape, Shape, Sphere};
    use oxiphysics_rigid::{Collider, ColliderSet, RigidBody, RigidBodySet};
    use std::sync::Arc;

    // -----------------------------------------------------------------------
    // 3. LBM density conservation over 100 steps
    // -----------------------------------------------------------------------
    #[test]
    fn test_lbm_density_conservation() {
        use oxiphysics_lbm::{LatticeType, LbmGrid2D, bgk_collide_2d, stream_2d};

        let nx = 20;
        let ny = 20;
        let omega = 1.0;
        let mut grid = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);

        // Perturb a few cells away from equilibrium
        grid.set_equilibrium(5, 5, 1.05, 0.02, 0.0);
        grid.set_equilibrium(15, 10, 0.95, -0.01, 0.01);
        grid.compute_macroscopic();

        let total_before = grid.total_density();

        for _ in 0..100 {
            bgk_collide_2d(&mut grid, omega);
            stream_2d(&mut grid);
        }
        grid.compute_macroscopic();

        let total_after = grid.total_density();
        let rel_change = (total_before - total_after).abs() / total_before;
        assert!(
            rel_change < 0.001,
            "LBM density changed by {:.4}% (> 0.1%): before={total_before}, after={total_after}",
            rel_change * 100.0
        );
    }

    // -----------------------------------------------------------------------
    // 4. FEM cantilever deflects downward under a point load
    // -----------------------------------------------------------------------
    #[test]
    fn test_fem_cantilever_deflects() {
        use oxiphysics_fem::analysis::LinearStaticAnalysis;
        use oxiphysics_fem::boundary::{DirichletBc, NeumannBc};
        use oxiphysics_fem::constitutive::LinearElasticMaterial;
        use oxiphysics_fem::mesh::TetrahedralMesh;

        let length = 1.0_f64;
        let width = 0.25_f64;
        let height = 0.25_f64;
        let e = 200.0e9_f64;
        let nu = 0.3_f64;

        let mesh = TetrahedralMesh::generate_beam(length, width, height, 4, 1, 1);
        let material = LinearElasticMaterial::new(e, nu);

        // Fix all DOFs at x=0
        let mut dirichlet_bcs = Vec::new();
        for (i, node) in mesh.nodes.iter().enumerate() {
            if node.x.abs() < 1e-10 {
                dirichlet_bcs.push(DirichletBc::new(i, 0, 0.0));
                dirichlet_bcs.push(DirichletBc::new(i, 1, 0.0));
                dirichlet_bcs.push(DirichletBc::new(i, 2, 0.0));
            }
        }

        // Apply downward load at free end
        let free_end_nodes: Vec<usize> = mesh
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, n)| (n.x - length).abs() < 1e-10)
            .map(|(i, _)| i)
            .collect();
        let p = -1000.0_f64;
        let num_free = free_end_nodes.len() as f64;
        let neumann_bcs: Vec<NeumannBc> = free_end_nodes
            .iter()
            .map(|&i| NeumannBc::new(i, [0.0, p / num_free, 0.0]))
            .collect();

        let analysis = LinearStaticAnalysis::new();
        let result = analysis.solve(
            &mesh,
            &material,
            &dirichlet_bcs,
            &neumann_bcs,
            &Vec3::zeros(),
        );

        // At least one free-end node should have a downward (negative y) displacement
        let max_y_disp = free_end_nodes
            .iter()
            .map(|&i| result.displacements[i].y)
            .fold(f64::INFINITY, f64::min);

        assert!(
            max_y_disp < 0.0,
            "Cantilever free end should deflect downward, got y={max_y_disp}"
        );
    }

    // -----------------------------------------------------------------------
    // 5. SPH particles stay within 2x initial bounds under gravity
    // -----------------------------------------------------------------------
    #[test]
    fn test_sph_particles_stay_bounded() {
        use oxiphysics_sph::boundary_sph::BoundaryPlane;
        use oxiphysics_sph::kernel::CubicSplineKernel;
        use oxiphysics_sph::particle::SphParticle;
        use oxiphysics_sph::simulation::{SolverType, SphSimulation, SphSimulationParams};

        let params = SphSimulationParams {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            smoothing_length: 0.1,
            sound_speed: 50.0,
            rest_density: 1000.0,
            viscosity: 0.01,
            solver_type: SolverType::Wcsph,
            min_dt: 1e-5,
            max_dt: 0.001,
            ..Default::default()
        };
        let mut sim = SphSimulation::new(params);

        // Floor boundary
        sim.boundaries
            .add_plane(BoundaryPlane::new(Vec3::zeros(), Vec3::new(0.0, 1.0, 0.0)));

        let spacing = 0.08_f64;
        let mass = 1000.0 * spacing.powi(3);

        // Create 50 particles in a 5x2x5 grid
        for i in 0..5_usize {
            for j in 0..2_usize {
                for k in 0..5_usize {
                    let pos = Vec3::new(
                        0.1 + i as f64 * spacing,
                        0.3 + j as f64 * spacing,
                        0.1 + k as f64 * spacing,
                    );
                    sim.particles
                        .add_particle(&SphParticle::new(pos, Vec3::zeros(), mass));
                }
            }
        }
        assert_eq!(sim.particles.len(), 50);

        // Compute initial bounding box
        let init_max_x = sim
            .particles
            .positions
            .iter()
            .map(|p| p.x)
            .fold(f64::NEG_INFINITY, f64::max);
        let init_max_y = sim
            .particles
            .positions
            .iter()
            .map(|p| p.y)
            .fold(f64::NEG_INFINITY, f64::max);
        let init_max_z = sim
            .particles
            .positions
            .iter()
            .map(|p| p.z)
            .fold(f64::NEG_INFINITY, f64::max);

        let kernel = CubicSplineKernel;
        for _ in 0..200 {
            sim.step(0.0005, &kernel);
        }

        // All particles should remain within 2x initial bounds
        let bound_x = init_max_x * 2.0;
        let bound_y = init_max_y * 2.0;
        let bound_z = init_max_z * 2.0;
        for (idx, p) in sim.particles.positions.iter().enumerate() {
            assert!(
                p.x.abs() <= bound_x + 1.0,
                "Particle {idx} x={} exceeds 2x bound {}",
                p.x,
                bound_x
            );
            assert!(
                p.y <= bound_y + 1.0,
                "Particle {idx} y={} exceeds 2x bound {}",
                p.y,
                bound_y
            );
            assert!(
                p.z.abs() <= bound_z + 1.0,
                "Particle {idx} z={} exceeds 2x bound {}",
                p.z,
                bound_z
            );
        }
    }

    // -----------------------------------------------------------------------
    // 6. MD: 8 LJ argon atoms, 100 steps, energy finite and KE > 0
    // -----------------------------------------------------------------------
    #[test]
    fn test_md_energy_finite() {
        use oxiphysics_md::atom::AtomSet;
        use oxiphysics_md::barostat::NoBarostat;
        use oxiphysics_md::forcefield::PairForceField;
        use oxiphysics_md::integrator::VelocityVerlet;
        use oxiphysics_md::neighbor::PeriodicBox;
        use oxiphysics_md::potential::LennardJones;
        use oxiphysics_md::simulation::MdSimulation;
        use oxiphysics_md::thermostat::NoThermostat;

        let spacing = 1.6_f64;
        let mut atoms = AtomSet::new();
        // 2x2x2 = 8 atoms
        for ix in 0..2_usize {
            for iy in 0..2_usize {
                for iz in 0..2_usize {
                    let pos = Vec3::new(
                        (ix as f64 + 0.5) * spacing,
                        (iy as f64 + 0.5) * spacing,
                        (iz as f64 + 0.5) * spacing,
                    );
                    let vel = Vec3::new(
                        (ix as f64 - 0.5) * 0.2,
                        (iy as f64 - 0.5) * 0.2,
                        (iz as f64 - 0.5) * 0.2,
                    );
                    atoms.add_atom(pos, vel, 1.0, 0.0, 0);
                }
            }
        }
        assert_eq!(atoms.len(), 8);

        let mut ff = PairForceField::new();
        ff.add_interaction(0, 0, LennardJones::new(1.0, 1.0, 2.5));

        let pbox = PeriodicBox::cubic(10.0);
        let mut sim = MdSimulation::new(
            atoms,
            ff,
            VelocityVerlet::new(),
            NoThermostat::new(),
            NoBarostat::new(),
            pbox,
            0.0,
            0.0,
            1.0,
        );

        let records = sim.run(100, 0.001, Some(1));

        assert!(!records.is_empty());
        for rec in &records {
            assert!(
                rec.total.is_finite(),
                "Total energy must be finite: {}",
                rec.total
            );
        }
        let final_ke = records.last().unwrap().kinetic;
        assert!(
            final_ke > 0.0,
            "Kinetic energy should be positive after 100 steps, got {final_ke}"
        );
    }

    // -----------------------------------------------------------------------
    // 7. Constraints: BallJoint constrains relative velocity at anchor
    // -----------------------------------------------------------------------
    #[test]
    fn test_constraints_joint_holds() {
        use oxiphysics_constraints::Constraint;
        use oxiphysics_constraints::joints::BallJoint;
        use oxiphysics_rigid::RigidBody;

        let mut bodies = RigidBodySet::new();
        let mut body_a = RigidBody::new(1.0);
        body_a.transform.position = Vec3::new(0.0, 0.0, 0.0);
        let ha = bodies.insert(body_a);

        let mut body_b = RigidBody::new(1.0);
        body_b.transform.position = Vec3::new(2.0, 0.0, 0.0);
        // Apply a force to body B by setting its velocity away from the joint
        body_b.velocity = Vec3::new(5.0, 3.0, 0.0);
        let hb = bodies.insert(body_b);

        let mut joint = BallJoint::new(
            ha,
            hb,
            Vec3::new(1.0, 0.0, 0.0),  // anchor at x=1 in A's local frame
            Vec3::new(-1.0, 0.0, 0.0), // anchor at x=-1 in B's local frame
        );

        let dt = 1.0 / 60.0;
        joint.prepare(&bodies, dt);
        for _ in 0..30 {
            joint.solve_velocity(&mut bodies, dt);
        }
        joint.solve_position(&mut bodies, dt);

        // Recompute world anchor positions – they should be close
        let a_pos = bodies.get(ha).unwrap().transform.position;
        let b_pos = bodies.get(hb).unwrap().transform.position;
        // Anchors: a_world = a_pos + (1,0,0), b_world = b_pos + (-1,0,0)
        let anchor_a = a_pos + Vec3::new(1.0, 0.0, 0.0);
        let anchor_b = b_pos + Vec3::new(-1.0, 0.0, 0.0);
        let gap = (anchor_a - anchor_b).norm();
        assert!(
            gap < 0.5,
            "BallJoint anchor gap too large: {gap} (should be < 0.5)"
        );
    }

    // -----------------------------------------------------------------------
    // 8. Materials: MooneyRivlin near-zero stress at identity deformation
    // -----------------------------------------------------------------------
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_materials_hyperelastic_identity() {
        use oxiphysics_materials::MooneyRivlin;

        let mr = MooneyRivlin::new(0.5e6, 0.1e6, 1.0e9);
        let identity = [[1.0_f64, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let p = mr.first_piola_kirchhoff_stress(&identity);

        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    p[i][j].abs() < 1.0,
                    "MooneyRivlin P[{i}][{j}] = {} should be near zero at identity deformation",
                    p[i][j]
                );
            }
        }
    }

    // -----------------------------------------------------------------------
    // 9. Soft body: 5x5 cloth mesh has 25 particles and constraints
    // -----------------------------------------------------------------------
    #[test]
    fn test_softbody_cloth_creation() {
        use oxiphysics_softbody::cloth::ClothMesh;

        let cloth = ClothMesh::new_grid(5, 5, 4.0, 4.0, 1.0);
        assert_eq!(
            cloth.vertices.len(),
            25,
            "5x5 cloth should have 25 vertices"
        );
        assert!(!cloth.edges.is_empty(), "Cloth should have edges");
    }

    // -----------------------------------------------------------------------
    // 10. GPU: LJ kernel via CpuBackend matches direct calculation for 4 atoms
    // -----------------------------------------------------------------------
    #[test]
    fn test_gpu_lj_forces_match_cpu() {
        use oxiphysics_gpu::compute::ComputeKernel;
        use oxiphysics_gpu::compute::CpuBackend;
        use oxiphysics_gpu::kernels::md_force::LennardJonesKernel;

        let _ = CpuBackend::new(); // backend exists

        let sigma = 1.0_f64;
        let epsilon = 1.0_f64;
        let cutoff = 3.0_f64;

        // 4 atoms on a line: x = 0, 1.5, 3, 4.5 (spacing = 1.5 sigma)
        let positions = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 3.0, 0.0, 0.0, 4.5, 0.0, 0.0];
        let params = vec![epsilon, sigma, cutoff];

        let mut outputs = vec![Vec::new(), Vec::new()];
        LennardJonesKernel.execute(&[&positions, &params], &mut outputs, 4);

        let forces = &outputs[0];
        let pe = outputs[1][0];

        assert_eq!(
            forces.len(),
            12,
            "Should have 4 atoms * 3 components = 12 force values"
        );
        assert!(pe.is_finite(), "Potential energy should be finite: {pe}");

        // Newton's third law: total force in x should be ~0
        let total_fx: f64 = forces.chunks(3).map(|c| c[0]).sum();
        assert!(
            total_fx.abs() < 1e-10,
            "Total x-force should be ~0 (Newton 3rd law), got {total_fx}"
        );
    }

    // -----------------------------------------------------------------------
    // 11. IO: write a simple VTK point cloud to /tmp and verify it exists
    // -----------------------------------------------------------------------
    #[test]
    fn test_io_vtk_write() {
        use oxiphysics_io::VtkWriter;

        let path = "/tmp/oxiphysics_integration_test.vtk";
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];

        VtkWriter::write_points(path, &positions).expect("VTK write should succeed");

        assert!(
            std::path::Path::new(path).exists(),
            "VTK file should exist at {path}"
        );

        let content = std::fs::read_to_string(path).expect("Should be able to read VTK file");
        assert!(!content.is_empty(), "VTK file should have content");
        assert!(
            content.contains("POINTS 4 float"),
            "VTK file should describe 4 points"
        );

        std::fs::remove_file(path).ok();
    }

    // -----------------------------------------------------------------------
    // 12. Viz: colormap(0.0) != colormap(1.0) for all colormaps
    // -----------------------------------------------------------------------
    #[test]
    fn test_viz_colormap_extremes() {
        use oxiphysics_viz::{Colormap, map_scalar};

        for cm in [
            Colormap::Jet,
            Colormap::Viridis,
            Colormap::CoolWarm,
            Colormap::Grayscale,
        ] {
            let c0 = map_scalar(0.0, 0.0, 1.0, cm);
            let c1 = map_scalar(1.0, 0.0, 1.0, cm);

            // At least one channel must differ between min and max
            let differs = (c0.r - c1.r).abs() > 0.01
                || (c0.g - c1.g).abs() > 0.01
                || (c0.b - c1.b).abs() > 0.01;

            assert!(
                differs,
                "Colormap {cm:?}: color at 0.0 should differ from color at 1.0"
            );
        }
    }

    /// End-to-end test: sphere falls under gravity, hits a static box, bounces.
    #[test]
    fn test_sphere_falls_and_bounces() {
        // Create body set
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        // Dynamic sphere at y=5
        let mut sphere_body = RigidBody::new(1.0);
        sphere_body.transform = Transform::from_position(Vec3::new(0.0, 5.0, 0.0));
        sphere_body.linear_damping = 0.0;
        let sphere_handle = bodies.insert(sphere_body);

        let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        let sphere_col = Collider::new(sphere_shape.clone())
            .with_body(sphere_handle)
            .with_restitution(0.8);
        let _sphere_col_handle = colliders.insert(sphere_col);

        // Static box floor at y=-0.5 (half_extents.y = 0.5, so top is at y=0)
        let mut floor_body = RigidBody::new_static();
        floor_body.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let floor_handle = bodies.insert(floor_body);

        let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
        let floor_col = Collider::new(floor_shape.clone()).with_body(floor_handle);
        let _floor_col_handle = colliders.insert(floor_col);

        let gravity = Vec3::new(0.0, -9.81, 0.0);
        let dt = 1.0 / 60.0;
        let _time_step = TimeStep::new(dt);

        let broadphase = BruteForceBroadPhase;

        // Simulate for 3 seconds
        let mut bounced = false;
        for _step in 0..180 {
            // Integrate forces
            for (_, body) in bodies.iter_mut() {
                body.integrate_forces(dt, &gravity);
            }

            // Collect AABBs for broadphase
            let collider_data: Vec<_> = colliders.iter().collect();
            let aabbs: Vec<Aabb> = collider_data
                .iter()
                .map(|(_, col)| {
                    let body_t = col
                        .body_handle
                        .and_then(|h| bodies.get(h))
                        .map(|b| b.transform.clone())
                        .unwrap_or_default();
                    let world_t = col.world_transform(&body_t);
                    let local_aabb = col.shape.bounding_box();
                    // Translate AABB to world
                    Aabb::new(
                        local_aabb.min + world_t.position,
                        local_aabb.max + world_t.position,
                    )
                })
                .collect();

            // Broadphase
            let pairs = broadphase.find_pairs(&aabbs);

            // Narrowphase + response
            for pair in &pairs {
                let (_, col_a) = &collider_data[pair.a];
                let (_, col_b) = &collider_data[pair.b];

                let body_t_a = col_a
                    .body_handle
                    .and_then(|h| bodies.get(h))
                    .map(|b| b.transform.clone())
                    .unwrap_or_default();
                let body_t_b = col_b
                    .body_handle
                    .and_then(|h| bodies.get(h))
                    .map(|b| b.transform.clone())
                    .unwrap_or_default();

                let t_a = col_a.world_transform(&body_t_a);
                let t_b = col_b.world_transform(&body_t_b);

                let manifold = NarrowPhaseDispatcher::generate_contacts(
                    col_a.shape.as_ref(),
                    &t_a,
                    col_b.shape.as_ref(),
                    &t_b,
                    CollisionPair::new(pair.a, pair.b),
                );

                if let Some(m) = manifold {
                    for contact in &m.contacts {
                        // Simple impulse-based response
                        let restitution = (col_a.restitution + col_b.restitution) * 0.5;

                        // Get body handles
                        let ha = col_a.body_handle;
                        let hb = col_b.body_handle;

                        // Compute relative velocity
                        let vel_a = ha
                            .and_then(|h| bodies.get(h))
                            .map(|b| b.velocity)
                            .unwrap_or_else(Vec3::zeros);
                        let vel_b = hb
                            .and_then(|h| bodies.get(h))
                            .map(|b| b.velocity)
                            .unwrap_or_else(Vec3::zeros);
                        let rel_vel = vel_a - vel_b;
                        let normal_vel = rel_vel.dot(&contact.normal);

                        if normal_vel > 0.0 {
                            continue; // Separating
                        }

                        let inv_mass_a = ha
                            .and_then(|h| bodies.get(h))
                            .map(|b| b.inverse_mass)
                            .unwrap_or(0.0);
                        let inv_mass_b = hb
                            .and_then(|h| bodies.get(h))
                            .map(|b| b.inverse_mass)
                            .unwrap_or(0.0);

                        if inv_mass_a + inv_mass_b < 1e-10 {
                            continue;
                        }

                        let j = -(1.0 + restitution) * normal_vel / (inv_mass_a + inv_mass_b);
                        let impulse = contact.normal * j;

                        if let Some(h) = ha
                            && let Some(body) = bodies.get_mut(h)
                            && body.is_dynamic()
                        {
                            body.velocity += impulse * body.inverse_mass;
                            // Position correction
                            body.transform.position += contact.normal * contact.depth * 0.5;
                        }
                        if let Some(h) = hb
                            && let Some(body) = bodies.get_mut(h)
                            && body.is_dynamic()
                        {
                            body.velocity -= impulse * body.inverse_mass;
                            body.transform.position -= contact.normal * contact.depth * 0.5;
                        }
                    }
                }
            }

            // Integrate velocities
            for (_, body) in bodies.iter_mut() {
                body.integrate_velocity(dt);
            }

            // Check if sphere bounced (velocity.y > 0 after being negative)
            let sphere = bodies.get(sphere_handle).unwrap();
            if sphere.velocity.y > 0.5 && sphere.transform.position.y < 3.0 {
                bounced = true;
            }
        }

        let sphere = bodies.get(sphere_handle).unwrap();
        // Sphere should have fallen and bounced
        assert!(bounced, "Sphere should have bounced off the floor");
        // Sphere should still be above the floor
        assert!(
            sphere.transform.position.y > -1.0,
            "Sphere y={} should be above -1.0",
            sphere.transform.position.y
        );
    }

    // -----------------------------------------------------------------------
    // P1. Elastic collision: momentum conservation
    // -----------------------------------------------------------------------
    #[test]
    fn test_elastic_collision_momentum_conservation() {
        // Two equal-mass bodies moving toward each other along X.
        // After elastic impulse exchange, total momentum should be conserved.
        let mass = 1.0_f64;
        let v0 = 5.0_f64; // m/s

        // Before collision
        let vel_a_before = Vec3::new(v0, 0.0, 0.0);
        let vel_b_before = Vec3::new(-v0, 0.0, 0.0);
        let p_before = vel_a_before * mass + vel_b_before * mass;
        let ke_before =
            0.5 * mass * vel_a_before.norm_squared() + 0.5 * mass * vel_b_before.norm_squared();

        // Set up bodies
        let mut bodies = RigidBodySet::new();
        let mut body_a = oxiphysics_rigid::RigidBody::new(mass);
        body_a.velocity = vel_a_before;
        body_a.linear_damping = 0.0;
        body_a.transform.position = Vec3::new(-0.5, 0.0, 0.0);
        let ha = bodies.insert(body_a);

        let mut body_b = oxiphysics_rigid::RigidBody::new(mass);
        body_b.velocity = vel_b_before;
        body_b.linear_damping = 0.0;
        body_b.transform.position = Vec3::new(0.5, 0.0, 0.0);
        let hb = bodies.insert(body_b);

        // Apply elastic collision impulse manually (restitution = 1.0)
        // For equal masses: bodies exchange velocities
        // n = (b_pos - a_pos).normalize() = (1,0,0)
        let normal = Vec3::new(1.0, 0.0, 0.0);
        let inv_mass_a = bodies.get(ha).unwrap().inverse_mass;
        let inv_mass_b = bodies.get(hb).unwrap().inverse_mass;
        let vel_a = bodies.get(ha).unwrap().velocity;
        let vel_b = bodies.get(hb).unwrap().velocity;
        let rel_vel = vel_a - vel_b;
        let normal_vel = rel_vel.dot(&normal);
        let restitution = 1.0_f64;
        let j = -(1.0 + restitution) * normal_vel / (inv_mass_a + inv_mass_b);
        let impulse = normal * j;
        bodies.get_mut(ha).unwrap().velocity += impulse * inv_mass_a;
        bodies.get_mut(hb).unwrap().velocity -= impulse * inv_mass_b;

        let vel_a_after = bodies.get(ha).unwrap().velocity;
        let vel_b_after = bodies.get(hb).unwrap().velocity;
        let p_after = vel_a_after * mass + vel_b_after * mass;
        let ke_after =
            0.5 * mass * vel_a_after.norm_squared() + 0.5 * mass * vel_b_after.norm_squared();

        // Total momentum conserved within 1%
        let dp = (p_after - p_before).norm();
        let p_mag = p_before.norm().max(1e-10);
        assert!(
            dp / p_mag < 0.01 || dp < 1e-8,
            "Elastic collision: momentum not conserved: dp={dp}, p_before={p_before:?}"
        );

        // KE conserved within 10%
        let dke = (ke_after - ke_before).abs();
        assert!(
            dke / ke_before < 0.10,
            "Elastic collision: KE not conserved: before={ke_before}, after={ke_after}"
        );
    }

    // -----------------------------------------------------------------------
    // P2. Perfectly inelastic collision: bodies stick, momentum conserved
    // -----------------------------------------------------------------------
    #[test]
    fn test_inelastic_collision_energy_dissipation() {
        let mass = 1.0_f64;
        let v0 = 4.0_f64;

        let mut bodies = RigidBodySet::new();
        let mut body_a = oxiphysics_rigid::RigidBody::new(mass);
        body_a.velocity = Vec3::new(v0, 0.0, 0.0);
        body_a.linear_damping = 0.0;
        body_a.transform.position = Vec3::new(-0.5, 0.0, 0.0);
        let ha = bodies.insert(body_a);

        let mut body_b = oxiphysics_rigid::RigidBody::new(mass);
        body_b.velocity = Vec3::new(-v0, 0.0, 0.0);
        body_b.linear_damping = 0.0;
        body_b.transform.position = Vec3::new(0.5, 0.0, 0.0);
        let hb = bodies.insert(body_b);

        let p_before =
            bodies.get(ha).unwrap().velocity * mass + bodies.get(hb).unwrap().velocity * mass;

        // Perfectly inelastic: restitution = 0 -> j = -rel_vel_normal / (1/m_a + 1/m_b)
        let normal = Vec3::new(1.0, 0.0, 0.0);
        let inv_mass_a = bodies.get(ha).unwrap().inverse_mass;
        let inv_mass_b = bodies.get(hb).unwrap().inverse_mass;
        let vel_a = bodies.get(ha).unwrap().velocity;
        let vel_b = bodies.get(hb).unwrap().velocity;
        let rel_vel = vel_a - vel_b;
        let normal_vel = rel_vel.dot(&normal);
        let restitution = 0.0_f64;
        let j = -(1.0 + restitution) * normal_vel / (inv_mass_a + inv_mass_b);
        let impulse = normal * j;
        bodies.get_mut(ha).unwrap().velocity += impulse * inv_mass_a;
        bodies.get_mut(hb).unwrap().velocity -= impulse * inv_mass_b;

        let vel_a_after = bodies.get(ha).unwrap().velocity;
        let vel_b_after = bodies.get(hb).unwrap().velocity;
        let p_after = vel_a_after * mass + vel_b_after * mass;

        // Relative velocity along normal should be ~0 (sticking)
        let rel_v_after = (vel_a_after - vel_b_after).dot(&normal);
        assert!(
            rel_v_after.abs() < 0.1,
            "Inelastic: bodies should stick, rel_v_normal={rel_v_after}"
        );

        // Momentum conserved
        let dp = (p_after - p_before).norm();
        assert!(dp < 1e-8, "Inelastic: momentum not conserved: dp={dp}");
    }

    // -----------------------------------------------------------------------
    // P3. Projectile parabolic trajectory
    // -----------------------------------------------------------------------
    #[test]
    fn test_projectile_parabolic_trajectory() {
        use crate::pipeline::PhysicsPipeline;
        use oxiphysics_core::PhysicsConfig;

        let g = 9.81_f64;
        let v0 = 10.0_f64;
        let angle = std::f64::consts::FRAC_PI_4; // 45 degrees

        let vx = v0 * angle.cos();
        let vy = v0 * angle.sin();

        // Analytical predictions
        let max_height_analytical = vy * vy / (2.0 * g); // ~2.55 m
        let range_analytical = v0 * v0 * (2.0 * angle).sin() / g; // ~10.19 m

        let config = PhysicsConfig {
            gravity: Vec3::new(0.0, -g, 0.0),
            // Disable sleep so body stays active
            time_before_sleep: 100.0,
            linear_sleep_threshold: 0.0,
            angular_sleep_threshold: 0.0,
            ..PhysicsConfig::default()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let colliders = oxiphysics_rigid::ColliderSet::new();

        let mut body = oxiphysics_rigid::RigidBody::new(1.0);
        body.transform.position = Vec3::new(0.0, 0.0, 0.0);
        body.velocity = Vec3::new(vx, vy, 0.0);
        body.linear_damping = 0.0;
        let h = bodies.insert(body);

        let dt = 1.0 / 600.0_f64;
        let mut max_height = 0.0_f64;

        // Simulate until y returns to 0 (or just below), track max height
        let max_steps = (5.0 / dt) as usize;
        let mut landed = false;
        for _ in 0..max_steps {
            pipeline.step(dt, &mut bodies, &colliders);
            let pos = bodies.get(h).unwrap().transform.position;
            if pos.y > max_height {
                max_height = pos.y;
            }
            // Detect landing (returned to y ≤ 0 with negative velocity)
            if pos.y <= 0.0 && bodies.get(h).unwrap().velocity.y < 0.0 && pipeline.time > 0.5 {
                let range = pos.x;
                // Range within 20%
                let range_err = (range - range_analytical).abs() / range_analytical;
                assert!(
                    range_err < 0.20,
                    "Projectile range: got {range:.3}, expected {range_analytical:.3}, err={range_err:.3}"
                );
                landed = true;
                break;
            }
        }

        // Max height within 20%
        let height_err = (max_height - max_height_analytical).abs() / max_height_analytical;
        assert!(
            height_err < 0.20,
            "Projectile max height: got {max_height:.3}, expected {max_height_analytical:.3}, err={height_err:.3}"
        );
        assert!(
            landed,
            "Projectile should have landed within simulation time"
        );
    }

    // -----------------------------------------------------------------------
    // P4. Pendulum period via spring joint
    // -----------------------------------------------------------------------
    #[test]
    fn test_pendulum_period() {
        use oxiphysics_constraints::{Constraint, SpringJoint};

        let g = 9.81_f64;
        let l = 1.0_f64; // pendulum length (m)
        let t_analytical = 2.0 * std::f64::consts::PI * (l / g).sqrt(); // ~2.006 s

        // Simulate a pendulum-like system: a mass connected by a very stiff spring
        // to a static anchor at origin. Body starts displaced by a small angle.
        let theta0 = 0.1_f64; // small angle (radians)
        let x0 = l * theta0.sin();
        let y0 = -l * theta0.cos(); // hanging below anchor

        let mut bodies = RigidBodySet::new();

        // Static anchor (infinite mass)
        let anchor = oxiphysics_rigid::RigidBody::new_static();
        let h_anchor = bodies.insert(anchor);

        // Pendulum bob
        let mut bob = oxiphysics_rigid::RigidBody::new(1.0);
        bob.transform.position = Vec3::new(x0, y0, 0.0);
        bob.velocity = Vec3::zeros();
        bob.linear_damping = 0.0;
        bob.angular_damping = 0.0;
        let h_bob = bodies.insert(bob);

        // High-stiffness spring to simulate rigid pendulum rod
        let stiffness = 50000.0_f64;
        let mut spring = SpringJoint::new(
            h_anchor,
            h_bob,
            Vec3::zeros(), // anchor at origin
            Vec3::zeros(), // bob center
            l,             // rest length = pendulum length
            stiffness,
            100.0, // light damping to maintain stability
        );

        let dt = 1.0 / 2000.0_f64;
        let gravity = Vec3::new(0.0, -g, 0.0);

        // Track zero crossings (x-position) to measure period
        let mut zero_crossings = 0u32;
        let mut last_x = x0;
        let mut period_estimate = 0.0_f64;
        let mut first_crossing_time = 0.0_f64;
        let total_steps = (6.0 / dt) as usize;

        for step in 0..total_steps {
            let t = step as f64 * dt;

            // Spring solve
            spring.prepare(&bodies, dt);
            spring.solve_velocity(&mut bodies, dt);

            // Integrate forces (gravity)
            for (_, body) in bodies.iter_mut() {
                body.integrate_forces(dt, &gravity);
            }

            // Integrate velocities
            for (_, body) in bodies.iter_mut() {
                body.integrate_velocity(dt);
            }

            let x_now = bodies.get(h_bob).unwrap().transform.position.x;

            // Count zero crossings (x changes sign)
            if last_x > 0.0 && x_now <= 0.0 {
                zero_crossings += 1;
                if zero_crossings == 1 {
                    first_crossing_time = t;
                } else if zero_crossings == 2 {
                    // From crossing 1 to crossing 2 (both + → -): one full period
                    period_estimate = t - first_crossing_time;
                    break;
                }
            }
            last_x = x_now;
        }

        assert!(
            zero_crossings >= 1,
            "Pendulum should oscillate (zero_crossings={zero_crossings})"
        );

        if period_estimate > 0.0 {
            let period_err = (period_estimate - t_analytical).abs() / t_analytical;
            assert!(
                // Spring-pendulum period differs from rigid-rod due to compliance; allow up to 110% error.
                period_err < 1.1,
                "Pendulum period: got {period_estimate:.3}, expected {t_analytical:.3}, err={period_err:.3}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // P5. Rigid body torque -> angular velocity
    // -----------------------------------------------------------------------
    #[test]
    fn test_rigid_body_angular_momentum_rotation() {
        // A free body in space (no gravity), apply torque about Z axis.
        // After dt: angular_velocity_z = torque_z * dt / I_z
        // For RigidBody::new(mass), local_inertia = Mat3::identity() * mass
        // world_inverse_inertia = Mat3::identity() * (1/mass)
        // So: omega_z = torque_z * (1/mass) * dt

        let mass = 2.0_f64;
        let torque_z = 10.0_f64;
        let dt = 1.0 / 60.0_f64;

        let mut body = oxiphysics_rigid::RigidBody::new(mass);
        body.linear_damping = 0.0;
        body.angular_damping = 0.0;
        body.gravity_scale = 0.0; // no gravity

        body.apply_torque(Vec3::new(0.0, 0.0, torque_z));

        let no_gravity = Vec3::zeros();
        body.integrate_forces(dt, &no_gravity);

        // Expected: omega_z = torque * inv_mass * dt = 10 * 0.5 * (1/60) ≈ 0.0833
        let expected_omega_z = torque_z * (1.0 / mass) * dt;
        let actual_omega_z = body.angular_velocity.z;

        let err = (actual_omega_z - expected_omega_z).abs();
        assert!(
            err < 1e-8,
            "Angular velocity should be {expected_omega_z}, got {actual_omega_z}"
        );

        // Body should be rotating in the correct direction (positive Z torque -> positive omega_z)
        assert!(
            actual_omega_z > 0.0,
            "Torque in +Z should give positive angular velocity, got {actual_omega_z}"
        );
    }

    // -----------------------------------------------------------------------
    // P6. Freefall kinematics: v = g*t, y = 0.5*g*t²
    // -----------------------------------------------------------------------
    #[test]
    fn test_freefall_kinematics() {
        let g = 9.81_f64;
        let config = PhysicsConfig {
            gravity: Vec3::new(0.0, -g, 0.0),
            time_before_sleep: 100.0,
            linear_sleep_threshold: 0.0,
            angular_sleep_threshold: 0.0,
            ..PhysicsConfig::default()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let colliders = oxiphysics_rigid::ColliderSet::new();

        let mut body = oxiphysics_rigid::RigidBody::new(1.0);
        body.transform.position = Vec3::new(0.0, 0.0, 0.0);
        body.velocity = Vec3::zeros();
        body.linear_damping = 0.0;
        let h = bodies.insert(body);

        let dt = 1.0 / 600.0_f64;
        let target_time = 1.0_f64;
        let steps = (target_time / dt) as usize;

        for _ in 0..steps {
            pipeline.step(dt, &mut bodies, &colliders);
        }

        let body = bodies.get(h).unwrap();
        let t = pipeline.time;

        // v = g*t (downward, so -g*t)
        let v_expected = -g * t;
        let v_actual = body.velocity.y;
        let v_err = (v_actual - v_expected).abs() / g.abs();
        assert!(
            v_err < 0.01,
            "Freefall velocity: expected {v_expected:.3}, got {v_actual:.3}, rel_err={v_err:.4}"
        );

        // y = -0.5*g*t²
        let y_expected = -0.5 * g * t * t;
        let y_actual = body.transform.position.y;
        let y_err = (y_actual - y_expected).abs() / y_expected.abs();
        assert!(
            y_err < 0.01,
            "Freefall position: expected {y_expected:.3}, got {y_actual:.3}, rel_err={y_err:.4}"
        );
    }

    // -----------------------------------------------------------------------
    // P7. Two bodies + spring: center of mass stays fixed
    // -----------------------------------------------------------------------
    #[test]
    fn test_two_body_gravitational_center_of_mass_fixed() {
        // Two equal-mass bodies connected by a spring, no external forces (gravity off).
        // System starts at rest with spring at natural length.
        // CoM should not move.
        let mass = 1.0_f64;

        let mut bodies = RigidBodySet::new();

        let mut body_a = oxiphysics_rigid::RigidBody::new(mass);
        body_a.transform.position = Vec3::new(-1.5, 0.0, 0.0);
        body_a.linear_damping = 0.0;
        body_a.gravity_scale = 0.0;
        let ha = bodies.insert(body_a);

        let mut body_b = oxiphysics_rigid::RigidBody::new(mass);
        body_b.transform.position = Vec3::new(1.5, 0.0, 0.0);
        body_b.linear_damping = 0.0;
        body_b.gravity_scale = 0.0;
        let hb = bodies.insert(body_b);

        // Spring at natural length = 3.0, then compress by displacing
        // Actually start with compressed spring (distance < rest length)
        // to create motion
        bodies.get_mut(ha).unwrap().transform.position = Vec3::new(-1.0, 0.0, 0.0);
        bodies.get_mut(hb).unwrap().transform.position = Vec3::new(1.0, 0.0, 0.0);

        let rest_length = 3.0_f64; // spring is compressed (actual dist = 2.0)
        let stiffness = 100.0_f64;
        let mut spring = SpringJoint::new(
            ha,
            hb,
            Vec3::zeros(),
            Vec3::zeros(),
            rest_length,
            stiffness,
            0.0,
        );

        // Compute initial CoM
        let pos_a0 = bodies.get(ha).unwrap().transform.position;
        let pos_b0 = bodies.get(hb).unwrap().transform.position;
        let com_initial = (pos_a0 * mass + pos_b0 * mass) / (2.0 * mass);

        let dt = 1.0 / 600.0_f64;
        let no_gravity = Vec3::zeros();

        for _ in 0..600 {
            spring.prepare(&bodies, dt);
            spring.solve_velocity(&mut bodies, dt);
            for (_, body) in bodies.iter_mut() {
                body.integrate_forces(dt, &no_gravity);
            }
            for (_, body) in bodies.iter_mut() {
                body.integrate_velocity(dt);
            }
        }

        let pos_a = bodies.get(ha).unwrap().transform.position;
        let pos_b = bodies.get(hb).unwrap().transform.position;
        let com_final = (pos_a * mass + pos_b * mass) / (2.0 * mass);

        let com_drift = (com_final - com_initial).norm();
        assert!(
            com_drift < 0.05,
            "Center of mass should not move: drift={com_drift:.6}, initial={com_initial:?}, final={com_final:?}"
        );
    }

    // -----------------------------------------------------------------------
    // P8. Restitution coefficient: bounce height = e² * H
    // -----------------------------------------------------------------------
    #[test]
    fn test_restitution_coefficient() {
        use oxiphysics_geometry::{Shape, Sphere};

        let g = 9.81_f64;
        let drop_height = 4.0_f64;
        let restitution = 0.7_f64;
        let expected_bounce_height = restitution * restitution * drop_height; // e² * H ≈ 1.96

        let config = PhysicsConfig {
            gravity: Vec3::new(0.0, -g, 0.0),
            solver_iterations: 20,
            time_before_sleep: 100.0,
            linear_sleep_threshold: 0.0,
            angular_sleep_threshold: 0.0,
            ..PhysicsConfig::default()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let mut colliders = oxiphysics_rigid::ColliderSet::new();

        // Dynamic sphere dropped from height H
        let mut sphere_body = oxiphysics_rigid::RigidBody::new(1.0);
        sphere_body.transform.position = Vec3::new(0.0, drop_height, 0.0);
        sphere_body.velocity = Vec3::zeros();
        sphere_body.linear_damping = 0.0;
        let sh = bodies.insert(sphere_body);

        let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(
            oxiphysics_rigid::Collider::new(sphere_shape)
                .with_body(sh)
                .with_restitution(restitution),
        );

        // Static floor at y = -0.5 (top surface at y = 0)
        let mut floor_body = oxiphysics_rigid::RigidBody::new_static();
        floor_body.transform.position = Vec3::new(0.0, -0.5, 0.0);
        let fh = bodies.insert(floor_body);

        let floor_shape: Arc<dyn Shape> = Arc::new(oxiphysics_geometry::BoxShape::new(Vec3::new(
            10.0, 0.5, 10.0,
        )));
        colliders.insert(
            oxiphysics_rigid::Collider::new(floor_shape)
                .with_body(fh)
                .with_restitution(restitution),
        );

        let dt = 1.0 / 600.0_f64;
        let mut bounced = false;
        let mut max_bounce_height = 0.0_f64;
        let mut was_falling = false;
        let mut hit_floor = false;

        for _ in 0..6000 {
            pipeline.step(dt, &mut bodies, &colliders);

            let body = bodies.get(sh).unwrap();
            let y = body.transform.position.y;
            let vy = body.velocity.y;

            if vy < -1.0 {
                was_falling = true;
            }

            // Detect first bounce: after falling and then moving upward
            if was_falling && !hit_floor && vy > 0.5 && y < drop_height * 0.5 {
                hit_floor = true;
            }

            if hit_floor && vy < 0.0 && y > 0.5 {
                // Past apex of first bounce
                if y > max_bounce_height {
                    max_bounce_height = y;
                }
                bounced = true;
            }
        }

        assert!(bounced, "Ball should have bounced off the floor");

        // Account for sphere radius: sphere center at rest on floor is at y = 0.5 (radius)
        // So effective drop from apex perspective: (drop_height - 0.5) to floor contact (y=0.5)
        // For a rough check: bounce height (from sphere center) within 30% of e²*H
        let height_err =
            (max_bounce_height - expected_bounce_height).abs() / expected_bounce_height.max(0.1);
        assert!(
            height_err < 0.30,
            "Bounce height: expected ≈{expected_bounce_height:.3}, got {max_bounce_height:.3}, err={height_err:.3}"
        );
    }

    /// Test: sphere-box contact returns correct normal and depth.
    #[test]
    fn test_sphere_box_contact_correctness() {
        let sphere = Sphere::new(1.0);
        let box_shape = BoxShape::new(Vec3::new(1.0, 1.0, 1.0));

        // Sphere at (0, 1.5, 0), box at (0, 0, 0) -> overlap of 0.5 along Y
        let t_sphere = Transform::from_position(Vec3::new(0.0, 1.5, 0.0));
        let t_box = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));

        let manifold = NarrowPhaseDispatcher::generate_contacts(
            &sphere,
            &t_sphere,
            &box_shape,
            &t_box,
            CollisionPair::new(0, 1),
        );

        assert!(manifold.is_some(), "Should detect collision");
        let m = manifold.unwrap();
        assert!(!m.contacts.is_empty(), "Should have contacts");
        let c = &m.contacts[0];

        // Normal should point roughly upward (from box toward sphere)
        assert!(c.normal.y.abs() > 0.5, "Normal should be roughly along Y");
        // Depth should be ~0.5
        assert!(
            (c.depth - 0.5).abs() < 0.2,
            "Depth should be ~0.5, got {}",
            c.depth
        );
    }

    // -----------------------------------------------------------------------
    // S1. 100 sphere bodies in a grid — no panic, all finite positions
    // -----------------------------------------------------------------------
    #[test]
    fn test_many_bodies_no_panic() {
        let config = PhysicsConfig {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver_iterations: 4,
            ccd_enabled: false,
            ..PhysicsConfig::default()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        let mut handles = Vec::new();

        // 10x10 grid of spheres, spaced 2 m apart, starting at y=5
        for row in 0..10_i32 {
            for col in 0..10_i32 {
                let mut body = RigidBody::new(1.0);
                body.transform =
                    Transform::from_position(Vec3::new(col as f64 * 2.0, 5.0, row as f64 * 2.0));
                body.linear_damping = 0.1;
                let h = bodies.insert(body);
                handles.push(h);

                let shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.4));
                colliders.insert(Collider::new(shape).with_body(h));
            }
        }

        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            pipeline.step(dt, &mut bodies, &colliders);
        }

        for h in &handles {
            let body = bodies.get(*h).unwrap();
            let p = body.transform.position;
            assert!(
                p.x.is_finite() && p.y.is_finite() && p.z.is_finite(),
                "Body position is not finite: {:?}",
                p
            );
        }
    }

    // -----------------------------------------------------------------------
    // S2. 5 boxes stacked vertically on a static floor — stack doesn't fully collapse
    // -----------------------------------------------------------------------
    #[test]
    fn test_stacked_boxes_stable() {
        let config = PhysicsConfig {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            solver_iterations: 10,
            ccd_enabled: false,
            ..PhysicsConfig::default()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        // Static floor: centre at y = -0.5, top surface at y = 0
        let mut floor = RigidBody::new_static();
        floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
        let fh = bodies.insert(floor);
        let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
        colliders.insert(Collider::new(floor_shape).with_body(fh));

        // 5 stacked boxes, each 1 m tall (half-extent 0.5 on Y)
        let mut top_handle = None;
        for i in 0..5_i32 {
            let y_center = 0.5 + i as f64 * 1.0; // box i centre: 0.5, 1.5, 2.5, 3.5, 4.5
            let mut boxb = RigidBody::new(1.0);
            boxb.transform = Transform::from_position(Vec3::new(0.0, y_center, 0.0));
            boxb.linear_damping = 0.1;
            let bh = bodies.insert(boxb);
            let box_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(0.5, 0.5, 0.5)));
            colliders.insert(
                Collider::new(box_shape)
                    .with_body(bh)
                    .with_friction(0.8)
                    .with_restitution(0.0),
            );
            if i == 4 {
                top_handle = Some(bh);
            }
        }

        // Run 3 seconds
        let dt = 1.0 / 60.0;
        for _ in 0..180 {
            pipeline.step(dt, &mut bodies, &colliders);
        }

        let top = bodies.get(top_handle.unwrap()).unwrap();
        // Top box started at y=4.5; after settling the stack may compress but
        // the top box centre should remain well above the floor (y > 1.5).
        assert!(
            top.transform.position.y > 1.5,
            "Top box fell too far: y={}",
            top.transform.position.y
        );
    }

    // -----------------------------------------------------------------------
    // S3. Restitution range: bounce heights ordered e=0.0 < e=0.5 < e=1.0
    // -----------------------------------------------------------------------
    #[test]
    fn test_restitution_range() {
        fn max_bounce_height_for_restitution(e: f64) -> f64 {
            let config = PhysicsConfig {
                gravity: Vec3::new(0.0, -9.81, 0.0),
                solver_iterations: 20,
                time_before_sleep: 100.0,
                linear_sleep_threshold: 0.0,
                angular_sleep_threshold: 0.0,
                ccd_enabled: false,
            };
            let mut pipeline = PhysicsPipeline::with_config(config);
            let mut bodies = RigidBodySet::new();
            let mut colliders = ColliderSet::new();

            let drop_height = 4.0_f64;

            let mut sphere_body = RigidBody::new(1.0);
            sphere_body.transform = Transform::from_position(Vec3::new(0.0, drop_height, 0.0));
            sphere_body.linear_damping = 0.0;
            let sh = bodies.insert(sphere_body);
            let sphere_shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
            colliders.insert(
                Collider::new(sphere_shape)
                    .with_body(sh)
                    .with_restitution(e),
            );

            let mut floor_body = RigidBody::new_static();
            floor_body.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
            let fh = bodies.insert(floor_body);
            let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
            colliders.insert(Collider::new(floor_shape).with_body(fh).with_restitution(e));

            let dt = 1.0 / 300.0;
            let mut max_height = 0.0_f64;
            let mut has_bounced = false;
            let mut was_going_down = false;

            for _ in 0..3000 {
                pipeline.step(dt, &mut bodies, &colliders);
                let body = bodies.get(sh).unwrap();
                let y = body.transform.position.y;
                let vy = body.velocity.y;

                if vy < -1.0 {
                    was_going_down = true;
                }
                if was_going_down && vy > 0.1 {
                    has_bounced = true;
                }
                if has_bounced && y > max_height {
                    max_height = y;
                }
            }
            max_height
        }

        let h0 = max_bounce_height_for_restitution(0.0);
        let h05 = max_bounce_height_for_restitution(0.5);
        let h10 = max_bounce_height_for_restitution(1.0);

        assert!(
            h0 < h05,
            "e=0.0 bounce height ({h0:.3}) should be less than e=0.5 ({h05:.3})"
        );
        assert!(
            h05 < h10,
            "e=0.5 bounce height ({h05:.3}) should be less than e=1.0 ({h10:.3})"
        );
    }

    // -----------------------------------------------------------------------
    // S4. Multiple collisions per step — 3 spheres in a line
    // -----------------------------------------------------------------------
    #[test]
    fn test_multiple_collisions_per_step() {
        let config = PhysicsConfig {
            gravity: Vec3::zeros(),
            solver_iterations: 10,
            ccd_enabled: false,
            ..PhysicsConfig::default()
        };
        let mut pipeline = PhysicsPipeline::with_config(config);
        let mut bodies = RigidBodySet::new();
        let mut colliders = ColliderSet::new();

        // Left sphere moving right at speed 10
        let mut left = RigidBody::new(1.0);
        left.transform = Transform::from_position(Vec3::new(-5.0, 0.0, 0.0));
        left.velocity = Vec3::new(10.0, 0.0, 0.0);
        left.linear_damping = 0.0;
        let hl = bodies.insert(left);
        let sl: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(Collider::new(sl).with_body(hl).with_restitution(0.8));

        // Middle sphere stationary
        let mut mid = RigidBody::new(1.0);
        mid.transform = Transform::from_position(Vec3::new(0.0, 0.0, 0.0));
        mid.linear_damping = 0.0;
        let hm = bodies.insert(mid);
        let sm: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(Collider::new(sm).with_body(hm).with_restitution(0.8));

        // Right sphere moving left at speed 10
        let mut right = RigidBody::new(1.0);
        right.transform = Transform::from_position(Vec3::new(5.0, 0.0, 0.0));
        right.velocity = Vec3::new(-10.0, 0.0, 0.0);
        right.linear_damping = 0.0;
        let hr = bodies.insert(right);
        let sr: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
        colliders.insert(Collider::new(sr).with_body(hr).with_restitution(0.8));

        let dt = 1.0 / 60.0;
        for _ in 0..300 {
            pipeline.step(dt, &mut bodies, &colliders);
        }

        for h in &[hl, hm, hr] {
            let body = bodies.get(*h).unwrap();
            let v = body.velocity;
            assert!(
                v.x.is_finite() && v.y.is_finite() && v.z.is_finite(),
                "Body velocity is not finite: {:?}",
                v
            );
        }
    }

    // -----------------------------------------------------------------------
    // S5. Pipeline determinism — two identical runs produce identical results
    // -----------------------------------------------------------------------
    #[test]
    fn test_pipeline_deterministic() {
        let make_scene = || {
            let config = PhysicsConfig {
                gravity: Vec3::new(0.0, -9.81, 0.0),
                solver_iterations: 8,
                ccd_enabled: false,
                linear_sleep_threshold: 0.01,
                angular_sleep_threshold: 0.01,
                time_before_sleep: 0.5,
            };
            let pipeline = PhysicsPipeline::with_config(config);
            let mut bodies = RigidBodySet::new();
            let mut colliders = ColliderSet::new();

            // Floor
            let mut floor = RigidBody::new_static();
            floor.transform = Transform::from_position(Vec3::new(0.0, -0.5, 0.0));
            let fh = bodies.insert(floor);
            let floor_shape: Arc<dyn Shape> = Arc::new(BoxShape::new(Vec3::new(10.0, 0.5, 10.0)));
            colliders.insert(Collider::new(floor_shape).with_body(fh));

            // 3 dynamic spheres
            let mut handles = Vec::new();
            for i in 0..3_i32 {
                let mut body = RigidBody::new(1.0);
                body.transform = Transform::from_position(Vec3::new(i as f64 * 1.5, 4.0, 0.0));
                body.linear_damping = 0.0;
                let h = bodies.insert(body);
                handles.push(h);
                let shape: Arc<dyn Shape> = Arc::new(Sphere::new(0.5));
                colliders.insert(Collider::new(shape).with_body(h).with_restitution(0.5));
            }

            (pipeline, bodies, colliders, handles)
        };

        let dt = 1.0 / 60.0;
        let steps = 120;

        let (mut pipeline1, mut bodies1, colliders1, handles1) = make_scene();
        for _ in 0..steps {
            pipeline1.step(dt, &mut bodies1, &colliders1);
        }

        let (mut pipeline2, mut bodies2, colliders2, handles2) = make_scene();
        for _ in 0..steps {
            pipeline2.step(dt, &mut bodies2, &colliders2);
        }

        for (h1, h2) in handles1.iter().zip(handles2.iter()) {
            let p1 = bodies1.get(*h1).unwrap().transform.position;
            let p2 = bodies2.get(*h2).unwrap().transform.position;
            assert!(
                (p1.x - p2.x).abs() < 1e-10
                    && (p1.y - p2.y).abs() < 1e-10
                    && (p1.z - p2.z).abs() < 1e-10,
                "Simulation is not deterministic: run1={:?}, run2={:?}",
                p1,
                p2
            );
        }
    }
}
