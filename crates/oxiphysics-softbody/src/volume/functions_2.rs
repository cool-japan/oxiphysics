//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod extended_volume_tests {
    use super::*;
    use crate::volume::CollisionPlane;
    use crate::volume::CollisionSphere;
    use crate::volume::TetrahedralSoftBody;
    use crate::volume::VolumeParticle;
    fn make_p(pos: [f64; 3]) -> VolumeParticle {
        VolumeParticle {
            position: pos,
            prev_position: pos,
            mass: 1.0,
            fixed: false,
        }
    }
    fn make_p_fixed(pos: [f64; 3]) -> VolumeParticle {
        VolumeParticle {
            position: pos,
            prev_position: pos,
            mass: 1.0,
            fixed: true,
        }
    }
    fn unit_tet_triangles() -> Vec<[usize; 3]> {
        vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
    }
    fn unit_tet_positions() -> Vec<[f64; 3]> {
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }
    #[test]
    fn test_signed_surface_volume_matches_tet() {
        let positions = unit_tet_positions();
        let triangles = unit_tet_triangles();
        let sv = signed_surface_volume(&triangles, &positions).abs();
        assert!((sv - 1.0 / 6.0).abs() < 1e-10, "sv={sv}");
    }
    #[test]
    fn test_signed_surface_volume_scaling() {
        let positions = unit_tet_positions();
        let triangles = unit_tet_triangles();
        let v1 = signed_surface_volume(&triangles, &positions).abs();
        let scale = 2.0_f64;
        let scaled: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| [p[0] * scale, p[1] * scale, p[2] * scale])
            .collect();
        let v2 = signed_surface_volume(&triangles, &scaled).abs();
        assert!(
            (v2 - v1 * scale.powi(3)).abs() < 1e-10,
            "v1={v1}, v2={v2}, expected={}",
            v1 * scale.powi(3)
        );
    }
    #[test]
    fn test_signed_surface_volume_gradient_nonzero() {
        let positions = unit_tet_positions();
        let triangles = unit_tet_triangles();
        for i in 0..4 {
            let g = signed_surface_volume_gradient(i, &triangles, &positions);
            let mag = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            assert!(mag > 1e-10, "gradient for vertex {i} is zero");
        }
    }
    #[test]
    fn test_pressurized_volume_constraint_converges() {
        let mut positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.5, 0.0, 0.0],
            [0.0, 1.5, 0.0],
            [0.0, 0.0, 1.5],
        ];
        let triangles = unit_tet_triangles();
        let rest_volume = 1.0 / 6.0;
        let inv_masses = vec![1.0_f64; 4];
        let v_before = signed_surface_volume(&triangles, &positions).abs();
        let err_before = (v_before - rest_volume).abs();
        for _ in 0..50 {
            solve_pressurized_volume_constraint(
                &triangles,
                rest_volume,
                1.0,
                0.5,
                &inv_masses,
                &mut positions,
            );
        }
        let v_after = signed_surface_volume(&triangles, &positions).abs();
        let err_after = (v_after - rest_volume).abs();
        assert!(
            err_after < err_before,
            "err_before={err_before}, err_after={err_after}"
        );
    }
    #[test]
    fn test_volume_preservation_energy_at_rest() {
        let positions = unit_tet_positions();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let rest_volumes = vec![1.0 / 6.0_f64];
        let e = volume_preservation_energy(&tetrahedra, &rest_volumes, &positions, 100.0);
        assert!(e < 1e-20, "energy at rest should be ~0, got {e}");
    }
    #[test]
    fn test_volume_preservation_energy_positive_when_deformed() {
        let mut positions = unit_tet_positions();
        positions[3] = [0.0, 0.0, 2.0];
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let rest_volumes = vec![1.0 / 6.0_f64];
        let e = volume_preservation_energy(&tetrahedra, &rest_volumes, &positions, 100.0);
        assert!(e > 0.0, "energy should be positive when deformed, got {e}");
    }
    #[test]
    fn test_volume_preservation_gradient_size() {
        let positions = unit_tet_positions();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let rest_volumes = vec![1.0 / 6.0_f64];
        let g = volume_preservation_gradient(&tetrahedra, &rest_volumes, &positions, 100.0);
        assert_eq!(g.len(), 4);
    }
    #[test]
    fn test_boyle_pressure_law() {
        let p0 = 101325.0_f64;
        let v0 = 1.0_f64;
        let v1 = 0.5_f64;
        let p1 = boyle_pressure(p0, v0, v1);
        assert!((p1 * v1 - p0 * v0).abs() < 1e-6, "Boyle's law violated");
    }
    #[test]
    fn test_boyle_pressure_zero_volume() {
        let p = boyle_pressure(101325.0, 1.0, 0.0);
        assert_eq!(p, 0.0);
    }
    #[test]
    fn test_adiabatic_pressure_law() {
        let p0 = 101325.0_f64;
        let v0 = 1.0_f64;
        let v1 = 0.8_f64;
        let gamma = 1.4_f64;
        let p1 = adiabatic_pressure(p0, v0, v1, gamma);
        let lhs = p1 * v1.powf(gamma);
        let rhs = p0 * v0.powf(gamma);
        assert!((lhs - rhs).abs() / rhs < 1e-10, "adiabatic law violated");
    }
    #[test]
    fn test_adiabatic_vs_boyle_compression() {
        let p0 = 1.0_f64;
        let v0 = 1.0_f64;
        let v1 = 0.5_f64;
        let p_boyle = boyle_pressure(p0, v0, v1);
        let p_adiabatic = adiabatic_pressure(p0, v0, v1, 1.4);
        assert!(
            p_adiabatic > p_boyle,
            "adiabatic pressure should exceed isothermal under compression"
        );
    }
    #[test]
    fn test_gas_pressure_modifies_velocities() {
        let mut positions = unit_tet_positions();
        let triangles = unit_tet_triangles();
        let inv_masses = vec![1.0_f64; 4];
        let mut velocities = vec![[0.0_f64; 3]; 4];
        apply_gas_pressure_force(
            &triangles,
            1.0,
            1.0 / 6.0,
            &inv_masses,
            &mut positions,
            &mut velocities,
            0.01,
        );
        let total_delta: f64 = velocities
            .iter()
            .flat_map(|v| v.iter())
            .map(|x| x.abs())
            .sum();
        assert!(total_delta > 0.0, "gas pressure should change velocities");
    }
    #[test]
    fn test_plane_collision_pushes_up() {
        let plane = CollisionPlane {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
        };
        let mut positions = vec![[0.0_f64, 0.0, -0.1]];
        let mut velocities = vec![[0.0_f64, 0.0, -1.0]];
        let inv_masses = vec![1.0_f64];
        resolve_plane_collision(
            &plane,
            100.0,
            0.5,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        assert!(
            positions[0][2] >= 0.0,
            "particle should be above plane after collision"
        );
    }
    #[test]
    fn test_plane_collision_no_effect_above() {
        let plane = CollisionPlane {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
        };
        let mut positions = vec![[0.0_f64, 0.0, 1.0]];
        let mut velocities = vec![[0.0_f64, 0.0, -1.0]];
        let inv_masses = vec![1.0_f64];
        let pos_before = positions[0];
        resolve_plane_collision(
            &plane,
            100.0,
            0.5,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        assert_eq!(positions[0], pos_before);
    }
    #[test]
    fn test_sphere_collision_pushes_out() {
        let sphere = CollisionSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let mut positions = vec![[0.5_f64, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]];
        let inv_masses = vec![1.0_f64];
        resolve_sphere_collision(
            &sphere,
            100.0,
            0.5,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        let dist =
            (positions[0][0].powi(2) + positions[0][1].powi(2) + positions[0][2].powi(2)).sqrt();
        assert!(
            dist >= 1.0 - 1e-10,
            "particle should be outside sphere, dist={dist}"
        );
    }
    #[test]
    fn test_sphere_collision_no_effect_outside() {
        let sphere = CollisionSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let mut positions = vec![[2.0_f64, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]];
        let inv_masses = vec![1.0_f64];
        let pos_before = positions[0];
        resolve_sphere_collision(
            &sphere,
            100.0,
            0.5,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        assert_eq!(positions[0], pos_before);
    }
    #[test]
    fn test_particle_collision_separates() {
        let mut positions = vec![[0.0_f64, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]; 2];
        let inv_masses = vec![1.0_f64; 2];
        resolve_particle_collisions(
            1.0,
            100.0,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        let dx = positions[0][0] - positions[1][0];
        let dist = dx.abs();
        assert!(
            dist >= 2.0 - 1e-10,
            "particles should be separated, dist={dist}"
        );
    }
    #[test]
    fn test_compute_current_volumes() {
        let pts = unit_tet_positions();
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_p(p)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let body = TetrahedralSoftBody::new(particles, tetrahedra);
        let cvs = body.compute_current_volumes();
        assert_eq!(cvs.len(), 1);
        assert!(
            (cvs[0] - 1.0 / 6.0).abs() < 1e-10,
            "current volume={}",
            cvs[0]
        );
    }
    #[test]
    fn test_free_particle_count() {
        let mut particles: Vec<VolumeParticle> =
            unit_tet_positions().iter().map(|&p| make_p(p)).collect();
        particles[0] = make_p_fixed([0.0, 0.0, 0.0]);
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let body = TetrahedralSoftBody::new(particles, tetrahedra);
        assert_eq!(body.free_particle_count(), 3);
    }
    #[test]
    fn test_kinetic_energy_zero_at_rest() {
        let pts = unit_tet_positions();
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_p(p)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let body = TetrahedralSoftBody::new(particles, tetrahedra);
        assert!((body.approx_kinetic_energy()).abs() < 1e-30);
    }
    #[test]
    fn test_kinetic_energy_positive_after_step() {
        let pts = unit_tet_positions();
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_p(p)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        body.step(0.01, [0.0, -9.8, 0.0]);
        assert!(body.approx_kinetic_energy() > 0.0);
    }
    #[test]
    fn test_apply_volume_preservation_moves_positions() {
        let mut pts = unit_tet_positions();
        pts[3] = [0.0, 0.0, 3.0];
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_p(p)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let rest_volumes = vec![1.0 / 6.0_f64];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        let pos_before = body.particles[3].position;
        body.apply_volume_preservation(&rest_volumes, 1000.0);
        let pos_after = body.particles[3].position;
        let changed = (pos_before[0] - pos_after[0]).abs()
            + (pos_before[1] - pos_after[1]).abs()
            + (pos_before[2] - pos_after[2]).abs();
        assert!(
            changed > 0.0,
            "positions should change after volume preservation"
        );
    }
    #[test]
    fn test_body_plane_collision() {
        let mut pts = unit_tet_positions();
        for p in pts.iter_mut() {
            p[2] -= 0.5;
        }
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_p(p)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        let plane = CollisionPlane {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
        };
        body.apply_plane_collision(&plane);
        for p in &body.particles {
            assert!(
                p.position[2] >= -1e-10,
                "particle below floor: z={}",
                p.position[2]
            );
        }
    }
    #[test]
    fn test_step_with_pressure_positions_change() {
        let pts = unit_tet_positions();
        let particles: Vec<VolumeParticle> = pts.iter().map(|&p| make_p(p)).collect();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let triangles = unit_tet_triangles();
        let mut body = TetrahedralSoftBody::new(particles, tetrahedra);
        let pos_before: Vec<[f64; 3]> = body.particles.iter().map(|p| p.position).collect();
        body.step_with_pressure(&triangles, 10.0, 1.0 / 6.0, [0.0, 0.0, -9.8], 0.01);
        let pos_after: Vec<[f64; 3]> = body.particles.iter().map(|p| p.position).collect();
        let total_change: f64 = pos_before
            .iter()
            .zip(&pos_after)
            .flat_map(|(a, b)| a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()))
            .sum();
        assert!(
            total_change > 0.0,
            "positions should change after step_with_pressure"
        );
    }
    #[test]
    fn test_fixed_particle_unaffected_by_collision() {
        let plane = CollisionPlane {
            point: [0.0, 0.0, 0.0],
            normal: [0.0, 0.0, 1.0],
        };
        let mut positions = vec![[0.0_f64, 0.0, -1.0]];
        let mut velocities = vec![[0.0_f64; 3]];
        let inv_masses = vec![0.0_f64];
        let pos_before = positions[0];
        resolve_plane_collision(
            &plane,
            100.0,
            0.5,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        assert_eq!(positions[0], pos_before, "fixed particle should not move");
    }
    #[test]
    fn test_fixed_particles_no_collision() {
        let mut positions = vec![[0.0_f64, 0.0, 0.0], [0.1, 0.0, 0.0]];
        let mut velocities = vec![[0.0_f64; 3]; 2];
        let inv_masses = vec![0.0_f64, 0.0];
        let p_before = positions.clone();
        resolve_particle_collisions(
            1.0,
            100.0,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        assert_eq!(positions, p_before, "fixed particles should not move");
    }
    #[test]
    fn test_pressure_volume_work_compression() {
        let w = pressure_volume_work(1.0, 2.0, 1.0);
        assert!((w - 1.0).abs() < 1e-10, "work={w}");
    }
    #[test]
    fn test_adiabatic_pressure_at_rest() {
        let p0 = 101325.0_f64;
        let v0 = 1.0_f64;
        let p = adiabatic_pressure(p0, v0, v0, 1.4);
        assert!(
            (p - p0).abs() < 1e-6,
            "pressure at rest volume should equal p0, got {p}"
        );
    }
    #[test]
    fn test_volume_preservation_gradient_zero_at_rest() {
        let positions = unit_tet_positions();
        let tetrahedra = vec![[0usize, 1, 2, 3]];
        let rest_volumes = vec![1.0 / 6.0_f64];
        let g = volume_preservation_gradient(&tetrahedra, &rest_volumes, &positions, 100.0);
        for (i, gi) in g.iter().enumerate() {
            let mag = (gi[0] * gi[0] + gi[1] * gi[1] + gi[2] * gi[2]).sqrt();
            assert!(
                mag < 1e-10,
                "gradient at rest for vertex {i} should be zero, mag={mag}"
            );
        }
    }
    #[test]
    fn test_adiabatic_gamma_1_equals_boyle() {
        let p0 = 101325.0_f64;
        let v0 = 1.0_f64;
        let v1 = 0.7_f64;
        let p_boyle = boyle_pressure(p0, v0, v1);
        let p_adiabatic = adiabatic_pressure(p0, v0, v1, 1.0);
        assert!(
            (p_boyle - p_adiabatic).abs() < 1e-6,
            "gamma=1 adiabatic should equal Boyle's law"
        );
    }
    #[test]
    fn test_sphere_collision_velocity_response() {
        let sphere = CollisionSphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        };
        let mut positions = vec![[0.5_f64, 0.0, 0.0]];
        let mut velocities = vec![[-1.0_f64, 0.0, 0.0]];
        let inv_masses = vec![1.0_f64];
        resolve_sphere_collision(
            &sphere,
            0.0,
            1.0,
            &mut positions,
            &mut velocities,
            &inv_masses,
            0.01,
        );
        assert!(
            velocities[0][0] >= 0.0,
            "outward velocity should be non-negative after elastic collision"
        );
    }
    #[test]
    fn test_pressurized_volume_monotone_convergence() {
        let mut positions: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        ];
        let triangles = unit_tet_triangles();
        let rest_volume = 1.0 / 6.0_f64;
        let inv_masses = vec![1.0_f64; 4];
        let mut prev_err = f64::MAX;
        for _ in 0..30 {
            solve_pressurized_volume_constraint(
                &triangles,
                rest_volume,
                1.0,
                0.5,
                &inv_masses,
                &mut positions,
            );
            let v = signed_surface_volume(&triangles, &positions).abs();
            let err = (v - rest_volume).abs();
            assert!(err <= prev_err + 1e-10, "error should not increase");
            prev_err = err;
        }
    }
    #[test]
    fn test_signed_surface_volume_gradient_finite_diff() {
        let positions = unit_tet_positions();
        let triangles = unit_tet_triangles();
        let eps = 1e-6_f64;
        let v0 = signed_surface_volume(&triangles, &positions);
        let mut pos_plus = positions.clone();
        pos_plus[0][0] += eps;
        let v_plus = signed_surface_volume(&triangles, &pos_plus);
        let fd = (v_plus - v0) / eps;
        let g = signed_surface_volume_gradient(0, &triangles, &positions);
        assert!((fd - g[0]).abs() < 1e-5, "FD={fd}, analytic={}", g[0]);
    }
}
