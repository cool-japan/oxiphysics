//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use super::types::{MpmGrid, MpmParticle};

/// G2P with explicit dt parameter (same as g2p_transfer above).
pub fn g2p_transfer_with_dt(grid: &MpmGrid, particles: &mut [MpmParticle], _dt: f64) {
    let h = grid.h;
    let inv_h = 1.0 / h;
    for p in particles.iter_mut() {
        if !p.active {
            continue;
        }
        let (base, frac) = grid.base_and_frac(p.position);
        let (w, _) = compute_weights(base, frac, h);
        let mut new_vel = [0.0f64; 3];
        let mut new_c = mat3_zero();
        for di in 0..3_usize {
            for dj in 0..3_usize {
                for dk in 0..3_usize {
                    let ix = base[0] + di as i64 - 1;
                    let iy = base[1] + dj as i64 - 1;
                    let iz = base[2] + dk as i64 - 1;
                    if !grid.in_bounds(ix, iy, iz) {
                        continue;
                    }
                    let w_ijk = w[0][di] * w[1][dj] * w[2][dk];
                    let node_idx = grid.idx(ix as usize, iy as usize, iz as usize);
                    let vi = grid.nodes[node_idx].velocity;
                    new_vel = add3(new_vel, scale3(vi, w_ijk));
                    let xip = sub3(
                        grid.node_pos(ix as usize, iy as usize, iz as usize),
                        p.position,
                    );
                    let outer = outer3(vi, xip);
                    new_c = mat3_add(new_c, mat3_scale(outer, w_ijk * 4.0 * inv_h * inv_h));
                }
            }
        }
        p.velocity = new_vel;
        p.affine_c = new_c;
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::material_point::MaterialType;
    use crate::material_point::MpmSimulation;
    use crate::material_point::MultiMaterialGrid;
    use crate::material_point::SandParams;
    use crate::material_point::SnowParams;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_add3() {
        let r = add3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert!((r[0] - 5.0).abs() < EPS);
        assert!((r[1] - 7.0).abs() < EPS);
        assert!((r[2] - 9.0).abs() < EPS);
    }
    #[test]
    fn test_dot3() {
        let d = dot3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(d.abs() < EPS);
        let d2 = dot3([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
        assert!((d2 - 14.0).abs() < EPS);
    }
    #[test]
    fn test_cross3() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(c[0].abs() < EPS);
        assert!(c[1].abs() < EPS);
        assert!((c[2] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_mat3_identity_mul() {
        let i = mat3_identity();
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let r = mat3_mul(i, m);
        for row in 0..3 {
            for col in 0..3 {
                assert!((r[row][col] - m[row][col]).abs() < EPS);
            }
        }
    }
    #[test]
    fn test_mat3_det_identity() {
        assert!((mat3_det(mat3_identity()) - 1.0).abs() < EPS);
    }
    #[test]
    fn test_mat3_det_zero() {
        let m = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        assert!(mat3_det(m).abs() < 1e-8);
    }
    #[test]
    fn test_mat3_inv_identity() {
        let i = mat3_identity();
        let inv = mat3_inv(i);
        for row in 0..3 {
            for col in 0..3 {
                let expected = if row == col { 1.0 } else { 0.0 };
                assert!((inv[row][col] - expected).abs() < 1e-8);
            }
        }
    }
    #[test]
    fn test_mat3_inv_round_trip() {
        let m = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]];
        let inv = mat3_inv(m);
        let prod = mat3_mul(m, inv);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (prod[i][j] - expected).abs() < 1e-8,
                    "prod[{i}][{j}]={}",
                    prod[i][j]
                );
            }
        }
    }
    #[test]
    fn test_outer3() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let o = outer3(a, b);
        assert!(o[0][1].abs() - 1.0 < EPS.abs());
        assert!(o[0][0].abs() < EPS);
    }
    #[test]
    fn test_polar_decomp_identity() {
        let (r, s) = polar_decomp(mat3_identity());
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((r[i][j] - expected).abs() < 1e-6, "r[{i}][{j}]={}", r[i][j]);
                assert!((s[i][j] - expected).abs() < 1e-6, "s[{i}][{j}]={}", s[i][j]);
            }
        }
    }
    #[test]
    fn test_bspline_quadratic_at_zero() {
        let w = bspline_quadratic(0.0);
        assert!((w - 0.75).abs() < EPS);
    }
    #[test]
    fn test_bspline_quadratic_partition_of_unity() {
        let frac = 0.3;
        let w0 = bspline_quadratic(frac + 1.0);
        let w1 = bspline_quadratic(frac);
        let w2 = bspline_quadratic(frac - 1.0);
        assert!((w0 + w1 + w2 - 1.0).abs() < 1e-6, "sum={}", w0 + w1 + w2);
    }
    #[test]
    fn test_bspline_quadratic_zero_outside_support() {
        assert!(bspline_quadratic(2.0).abs() < EPS);
        assert!(bspline_quadratic(-2.0).abs() < EPS);
    }
    #[test]
    fn test_bspline_cubic_at_zero() {
        let w = bspline_cubic(0.0);
        assert!((w - 2.0 / 3.0).abs() < EPS);
    }
    #[test]
    fn test_grid_creation() {
        let g = MpmGrid::new([8, 8, 8], 0.1, [0.0; 3]);
        assert_eq!(g.num_nodes(), 512);
    }
    #[test]
    fn test_grid_idx() {
        let g = MpmGrid::new([4, 5, 6], 0.1, [0.0; 3]);
        assert_eq!(g.idx(0, 0, 0), 0);
        assert_eq!(g.idx(1, 0, 0), 5 * 6);
        assert_eq!(g.idx(0, 1, 0), 6);
    }
    #[test]
    fn test_grid_node_pos() {
        let g = MpmGrid::new([8, 8, 8], 0.5, [1.0, 2.0, 3.0]);
        let pos = g.node_pos(2, 1, 0);
        assert!((pos[0] - 2.0).abs() < EPS);
        assert!((pos[1] - 2.5).abs() < EPS);
        assert!((pos[2] - 3.0).abs() < EPS);
    }
    #[test]
    fn test_grid_in_bounds() {
        let g = MpmGrid::new([8, 8, 8], 0.1, [0.0; 3]);
        assert!(g.in_bounds(0, 0, 0));
        assert!(g.in_bounds(7, 7, 7));
        assert!(!g.in_bounds(8, 0, 0));
        assert!(!g.in_bounds(-1, 0, 0));
    }
    #[test]
    fn test_grid_reset() {
        let mut g = MpmGrid::new([4, 4, 4], 0.1, [0.0; 3]);
        g.nodes[0].mass = 5.0;
        g.nodes[0].velocity = [1.0, 2.0, 3.0];
        g.reset();
        assert!(g.nodes[0].mass.abs() < EPS);
        for c in g.nodes[0].velocity {
            assert!(c.abs() < EPS);
        }
    }
    #[test]
    fn test_grid_base_and_frac() {
        let g = MpmGrid::new([16, 16, 16], 1.0, [0.0; 3]);
        let (base, frac) = g.base_and_frac([3.7, 2.1, 0.9]);
        assert_eq!(base[0], 3);
        assert_eq!(base[1], 2);
        assert_eq!(base[2], 0);
        assert!((frac[0] - 0.7).abs() < 1e-10);
    }
    #[test]
    fn test_particle_creation() {
        let p = MpmParticle::new(
            [0.5, 0.5, 0.5],
            [0.0; 3],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        );
        assert!((p.mass - 1.0).abs() < EPS);
        assert!(p.active);
    }
    #[test]
    fn test_particle_jacobian_identity() {
        let p = MpmParticle::new(
            [0.0; 3],
            [0.0; 3],
            1.0,
            1.0,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        );
        assert!((p.jacobian() - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_particle_volume() {
        let vol0 = 0.01;
        let p = MpmParticle::new(
            [0.0; 3],
            [0.0; 3],
            1.0,
            vol0,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        );
        assert!((p.volume() - vol0).abs() < EPS);
    }
    #[test]
    fn test_neo_hookean_pk1_identity_zero() {
        let f = mat3_identity();
        let pk1 = neo_hookean_pk1(f, 1.0, 1.0);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    pk1[i][j].abs() < 1e-8,
                    "pk1[{i}][{j}]={} should be ~0",
                    pk1[i][j]
                );
            }
        }
    }
    #[test]
    fn test_neo_hookean_energy_identity() {
        let f = mat3_identity();
        let e = neo_hookean_energy(f, 1.0, 1.0);
        assert!(e.abs() < 1e-10, "energy at identity should be 0, got {e}");
    }
    #[test]
    fn test_fixed_corotated_energy_identity() {
        let f = mat3_identity();
        let e = fixed_corotated_energy(f, 1.0, 1.0);
        assert!(e.abs() < 1e-10, "corotated energy at identity = {e}");
    }
    #[test]
    fn test_neo_hookean_energy_positive_under_stretch() {
        let mut f = mat3_identity();
        f[0][0] = 1.5;
        let e = neo_hookean_energy(f, 1.0e3, 1.0e3);
        assert!(e > 0.0, "energy should be positive under stretch");
    }
    #[test]
    fn test_init_particle_block_count() {
        let particles = init_particle_block(
            [0.0; 3],
            [1.0, 1.0, 1.0],
            0.5,
            1.0,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        );
        assert_eq!(particles.len(), 27);
    }
    #[test]
    fn test_init_particle_sphere_all_inside() {
        let center = [0.5, 0.5, 0.5];
        let radius = 0.4;
        let particles =
            init_particle_sphere(center, radius, 0.1, 1.0, MaterialType::NeoHookean, 1e5, 2e5);
        for p in &particles {
            let d = len3(sub3(p.position, center));
            assert!(d <= radius + 1e-10, "particle outside sphere: d={d}");
        }
    }
    #[test]
    fn test_init_particle_sphere_nonempty() {
        let center = [0.5, 0.5, 0.5];
        let radius = 0.3;
        let particles =
            init_particle_sphere(center, radius, 0.1, 1.0, MaterialType::NeoHookean, 1e5, 2e5);
        assert!(!particles.is_empty(), "sphere should have particles");
    }
    #[test]
    fn test_mpm_step_runs() {
        let mut grid = MpmGrid::new([16, 16, 16], 0.1, [0.0; 3]);
        let mut particles = init_particle_block(
            [0.4, 0.4, 0.4],
            [0.6, 0.6, 0.6],
            0.05,
            0.001,
            MaterialType::NeoHookean,
            1e4,
            2e4,
        );
        let dt = 1e-4;
        let gravity = [0.0, -9.81, 0.0];
        p2g_transfer(&mut grid, &particles, dt, gravity);
        grid.apply_boundary_conditions(2);
        g2p_transfer(&grid, &mut particles);
        update_deformation_gradient(&mut particles, dt, 0.025, 0.0075);
        advect_particles(&mut particles, dt);
        let avg_vy: f64 =
            particles.iter().map(|p| p.velocity[1]).sum::<f64>() / particles.len() as f64;
        assert!(avg_vy < 0.0, "avg vy={avg_vy} should be negative (gravity)");
    }
    #[test]
    fn test_mpm_simulation_context() {
        let mut sim = MpmSimulation::new([16, 16, 16], 0.1, [0.0; 3], [0.0, -9.81, 0.0]);
        let particles = init_particle_block(
            [0.4, 0.4, 0.4],
            [0.6, 0.6, 0.6],
            0.05,
            0.001,
            MaterialType::NeoHookean,
            1e4,
            2e4,
        );
        let n = particles.len();
        sim.add_particles(particles);
        assert_eq!(sim.num_active_particles(), n);
        let mass0 = sim.total_mass();
        sim.step(1e-4);
        let mass1 = sim.total_mass();
        assert!((mass0 - mass1).abs() < 1e-10, "mass should be conserved");
    }
    #[test]
    fn test_cfl_dt() {
        let mut particles = vec![MpmParticle::new(
            [0.5; 3],
            [10.0, 0.0, 0.0],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        )];
        particles.push(MpmParticle::new(
            [0.6; 3],
            [5.0, 0.0, 0.0],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        ));
        let dt = compute_cfl_dt(&particles, 0.1, 0.5);
        assert!((dt - 0.005).abs() < 1e-10, "dt={dt}");
    }
    #[test]
    fn test_snow_particle_creation() {
        let params = SnowParams::default();
        let p = make_snow_particle([0.5, 0.5, 0.5], [0.0; 3], 0.001, &params);
        assert_eq!(p.material, MaterialType::Snow);
        assert!((p.hardening - 10.0).abs() < EPS);
    }
    #[test]
    fn test_sand_particle_creation() {
        let params = SandParams::default();
        let p = make_sand_particle([0.5, 0.5, 0.5], [0.0; 3], 0.001, &params);
        assert_eq!(p.material, MaterialType::Sand);
        assert!((p.yield_param - 30.0).abs() < EPS);
    }
    #[test]
    fn test_drucker_prager_elastic_return() {
        let fe = mat3_identity();
        let fe_out = drucker_prager_return_mapping(fe, 1e5, 2e5, 30.0);
        for i in 0..3 {
            let expected = if i < 3 { 1.0 } else { 0.0 };
            assert!((fe_out[i][i] - expected).abs() < 1e-3);
        }
    }
    #[test]
    fn test_von_mises_elastic_return() {
        let fe = mat3_identity();
        let fe_out = von_mises_return_mapping(fe, 1e5, 2e5, 1e10);
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((fe_out[i][j] - expected).abs() < 1e-6);
            }
        }
    }
    #[test]
    fn test_snow_hardening_params_scaling() {
        let (mu_e, lambda_e) = snow_hardening_params(
            mat3_identity(),
            mat3_identity(),
            1.0,
            1.0,
            0.0,
            0.025,
            0.0075,
        );
        assert!((mu_e - 1.0).abs() < 1e-8);
        assert!((lambda_e - 1.0).abs() < 1e-8);
    }
    #[test]
    fn test_p2g_mass_conservation() {
        let mut grid = MpmGrid::new([16, 16, 16], 0.1, [0.0; 3]);
        let particles = init_particle_block(
            [0.4, 0.4, 0.4],
            [0.6, 0.6, 0.6],
            0.05,
            1.0,
            MaterialType::NeoHookean,
            1e4,
            2e4,
        );
        let total_particle_mass: f64 = particles.iter().map(|p| p.mass).sum();
        p2g_transfer(&mut grid, &particles, 1e-4, [0.0; 3]);
        let total_grid_mass: f64 = grid.nodes.iter().map(|n| n.mass).sum();
        let rel_err = (total_grid_mass - total_particle_mass).abs() / total_particle_mass;
        assert!(rel_err < 1e-10, "P2G mass conservation: rel_err={rel_err}");
    }
    #[test]
    fn test_advect_particles_zero_velocity() {
        let mut particles = vec![MpmParticle::new(
            [1.0, 2.0, 3.0],
            [0.0; 3],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        )];
        advect_particles(&mut particles, 0.01);
        let pos = particles[0].position;
        assert!((pos[0] - 1.0).abs() < EPS);
        assert!((pos[1] - 2.0).abs() < EPS);
        assert!((pos[2] - 3.0).abs() < EPS);
    }
    #[test]
    fn test_advect_particles_constant_velocity() {
        let mut particles = vec![MpmParticle::new(
            [0.0; 3],
            [1.0, 2.0, 3.0],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        )];
        let dt = 0.5;
        advect_particles(&mut particles, dt);
        let pos = particles[0].position;
        assert!((pos[0] - 0.5).abs() < EPS);
        assert!((pos[1] - 1.0).abs() < EPS);
        assert!((pos[2] - 1.5).abs() < EPS);
    }
    #[test]
    fn test_cull_out_of_bounds() {
        let mut sim = MpmSimulation::new([10, 10, 10], 0.1, [0.0; 3], [0.0, -9.81, 0.0]);
        let mut p = MpmParticle::new(
            [5.0, 5.0, 5.0],
            [0.0; 3],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        );
        p.active = true;
        sim.particles.push(p);
        sim.cull_out_of_bounds_particles();
        assert_eq!(sim.num_active_particles(), 0);
    }
    #[test]
    fn test_deformation_gradient_update_identity() {
        let mut particles = vec![MpmParticle::new(
            [0.5; 3],
            [0.0; 3],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        )];
        update_deformation_gradient(&mut particles, 0.01, 0.025, 0.0075);
        let f = particles[0].deformation_gradient;
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (f[i][j] - expected).abs() < 1e-10,
                    "F[{i}][{j}]={}",
                    f[i][j]
                );
            }
        }
    }
    #[test]
    fn test_multi_material_grid_creation() {
        let mmg = MultiMaterialGrid::new([8, 8, 8], 0.1, [0.0; 3], 3);
        assert_eq!(mmg.num_materials, 3);
        let n = 8 * 8 * 8;
        assert_eq!(mmg.velocities.len(), 3 * n);
    }
    #[test]
    fn test_level_set_all_max_for_empty_particles() {
        let grid = MpmGrid::new([8, 8, 8], 0.1, [0.0; 3]);
        let phi = compute_level_set(&[], &grid, 0.1);
        for &v in &phi {
            assert_eq!(v, f64::MAX);
        }
    }
    #[test]
    fn test_kinetic_energy_zero() {
        let mut sim = MpmSimulation::new([8, 8, 8], 0.1, [0.0; 3], [0.0; 3]);
        sim.particles.push(MpmParticle::new(
            [0.5; 3],
            [0.0; 3],
            1.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        ));
        assert!(sim.kinetic_energy().abs() < EPS);
    }
    #[test]
    fn test_kinetic_energy_nonzero() {
        let mut sim = MpmSimulation::new([8, 8, 8], 0.1, [0.0; 3], [0.0; 3]);
        sim.particles.push(MpmParticle::new(
            [0.5; 3],
            [1.0, 0.0, 0.0],
            2.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        ));
        assert!((sim.kinetic_energy() - 1.0).abs() < EPS);
    }
    #[test]
    fn test_total_momentum() {
        let mut sim = MpmSimulation::new([8, 8, 8], 0.1, [0.0; 3], [0.0; 3]);
        sim.particles.push(MpmParticle::new(
            [0.5; 3],
            [3.0, 0.0, 0.0],
            2.0,
            0.001,
            MaterialType::NeoHookean,
            1e5,
            2e5,
        ));
        let mom = sim.total_momentum();
        assert!((mom[0] - 6.0).abs() < EPS);
    }
}
