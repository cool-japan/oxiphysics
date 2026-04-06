//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_sim_gpu::FlipParticle;
    use crate::fluid_sim_gpu::FluidSimConfig;
    use crate::fluid_sim_gpu::FluidSimulation;
    use crate::fluid_sim_gpu::GpuBoundaryBox;
    use crate::fluid_sim_gpu::GpuNeighborList;
    use crate::fluid_sim_gpu::LbmCellType;
    use crate::fluid_sim_gpu::LbmD2Q9;
    use crate::fluid_sim_gpu::MacGrid;
    use crate::fluid_sim_gpu::MultiGpuDomain;
    use crate::fluid_sim_gpu::SphConfig;
    use crate::fluid_sim_gpu::SphKernels;
    use crate::fluid_sim_gpu::SphParticle;
    #[test]
    fn test_poly6_inside() {
        let w = SphKernels::poly6(0.0, 0.1);
        assert!(w > 0.0);
    }
    #[test]
    fn test_poly6_outside() {
        let w = SphKernels::poly6(0.2, 0.1);
        assert!((w).abs() < 1e-15);
    }
    #[test]
    fn test_spiky_grad_inside() {
        let r_vec = [0.05, 0.0, 0.0];
        let g = SphKernels::spiky_grad(r_vec, 0.05, 0.1);
        assert!(g[0] != 0.0);
    }
    #[test]
    fn test_viscosity_laplacian() {
        let lap = SphKernels::viscosity_laplacian(0.05, 0.1);
        assert!(lap > 0.0);
    }
    #[test]
    fn test_sph_compute_density_single_particle() {
        let mut particles = vec![SphParticle::new([0.0; 3], 1.0)];
        let config = SphConfig::default();
        sph_compute_density(&mut particles, &config);
        assert!(particles[0].density > 0.0);
    }
    #[test]
    fn test_sph_step_does_not_panic() {
        let mut particles = vec![
            SphParticle::new([0.0, 0.0, 0.0], 0.001),
            SphParticle::new([0.05, 0.0, 0.0], 0.001),
        ];
        let config = SphConfig::default();
        sph_step(&mut particles, &config);
    }
    #[test]
    fn test_lbm_d2q9_new() {
        let lbm = LbmD2Q9::new(8, 8, 1.0);
        assert_eq!(lbm.nx, 8);
        assert_eq!(lbm.ny, 8);
    }
    #[test]
    fn test_lbm_density_initial() {
        let lbm = LbmD2Q9::new(4, 4, 1.0);
        let rho = lbm.density(2, 2);
        assert!((rho - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_lbm_velocity_initial_zero() {
        let lbm = LbmD2Q9::new(4, 4, 1.0);
        let [ux, uy] = lbm.velocity(2, 2);
        assert!(ux.abs() < 1e-10);
        assert!(uy.abs() < 1e-10);
    }
    #[test]
    fn test_lbm_f_equilibrium_rest() {
        let f_eq = LbmD2Q9::f_equilibrium(1.0, 0.0, 0.0, 0);
        assert!((f_eq - D2Q9_W[0]).abs() < 1e-12);
    }
    #[test]
    fn test_lbm_step_does_not_panic() {
        let mut lbm = LbmD2Q9::new(8, 8, 0.6);
        lbm.set_inlet_velocity(0.05, 0.0);
        lbm.step();
    }
    #[test]
    fn test_lbm_bounce_back_solid() {
        let mut lbm = LbmD2Q9::new(4, 4, 1.0);
        lbm.cell_type[0] = LbmCellType::Solid;
        lbm.apply_bounce_back();
    }
    #[test]
    fn test_mac_grid_new() {
        let g = MacGrid::new(4, 4, 4, 0.1);
        assert_eq!(g.nx, 4);
        assert_eq!(g.u.len(), 5 * 4 * 4);
    }
    #[test]
    fn test_mac_grid_set_get_u() {
        let mut g = MacGrid::new(4, 4, 4, 0.1);
        g.set_u(2, 1, 1, 3.125);
        assert!((g.get_u(2, 1, 1) - 3.125).abs() < 1e-12);
    }
    #[test]
    fn test_mac_grid_compute_divergence_zero() {
        let mut g = MacGrid::new(4, 4, 4, 0.1);
        g.flags = vec![1u8; 4 * 4 * 4];
        g.compute_divergence();
        for &d in &g.div {
            assert!(d.abs() < 1e-12);
        }
    }
    #[test]
    fn test_jacobi_pressure_solve() {
        let mut g = MacGrid::new(4, 4, 4, 0.1);
        g.flags = vec![1u8; 4 * 4 * 4];
        g.jacobi_pressure_solve(1000.0, 0.01, 5);
    }
    #[test]
    fn test_compute_vorticity_zero_vel() {
        let g = MacGrid::new(8, 8, 8, 0.1);
        let vort = compute_vorticity(&g);
        assert_eq!(vort.len(), 8 * 8 * 8);
        for v in &vort {
            assert!(v[0].abs() < 1e-12);
        }
    }
    #[test]
    fn test_surface_tension_csf() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let phi = vec![-1.0f64; nx * ny * nz];
        let forces = surface_tension_csf(&phi, nx, ny, nz, 0.1, 0.0728);
        assert_eq!(forces.len(), nx * ny * nz);
    }
    #[test]
    fn test_flip_particle_new() {
        let p = FlipParticle::new([1.0, 2.0, 3.0]);
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.velocity, [0.0; 3]);
    }
    #[test]
    fn test_p2g_g2p_roundtrip() {
        let mut particles = vec![FlipParticle::new([0.25, 0.25, 0.25])];
        particles[0].velocity = [1.0, 0.0, 0.0];
        let mut grid = MacGrid::new(4, 4, 4, 0.1);
        p2g_transfer(&particles, &mut grid);
    }
    #[test]
    fn test_fluid_simulation_new() {
        let config = FluidSimConfig::default();
        let sim = FluidSimulation::new(config);
        assert_eq!(sim.particles.len(), 0);
        assert!((sim.time).abs() < 1e-15);
    }
    #[test]
    fn test_fluid_simulation_add_block() {
        let config = FluidSimConfig::default();
        let mut sim = FluidSimulation::new(config);
        sim.add_fluid_block([0.1; 3], [0.9; 3], 8);
        assert!(sim.particles.len() >= 8);
    }
    #[test]
    fn test_fluid_simulation_step() {
        let config = FluidSimConfig {
            grid_size: [8, 8, 8],
            pressure_iters: 5,
            ..Default::default()
        };
        let mut sim = FluidSimulation::new(config);
        sim.add_fluid_block([0.1; 3], [0.5; 3], 4);
        sim.step();
        assert_eq!(sim.step_count, 1);
        assert!(sim.time > 0.0);
    }
    #[test]
    fn test_fluid_kinetic_energy_increases_under_gravity() {
        let config = FluidSimConfig {
            grid_size: [8, 8, 8],
            pressure_iters: 5,
            ..Default::default()
        };
        let mut sim = FluidSimulation::new(config);
        sim.add_fluid_block([0.4; 3], [0.6; 3], 4);
        let ke0 = sim.kinetic_energy();
        sim.step();
        let ke1 = sim.kinetic_energy();
        assert!(ke1 >= ke0);
    }
    #[test]
    fn test_advect_scalar_uniform_field() {
        let nx = 4;
        let ny = 4;
        let nz = 4;
        let grid = MacGrid::new(nx, ny, nz, 0.1);
        let phi = vec![1.0f64; nx * ny * nz];
        let result = advect_scalar(&phi, &grid, 0.01, [0.0, -9.81, 0.0]);
        assert_eq!(result.len(), nx * ny * nz);
        for v in &result {
            assert!((v - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_sph_particle_gravity_descent() {
        let mut p = vec![SphParticle::new([0.0, 1.0, 0.0], 0.001)];
        let config = SphConfig {
            gravity: [0.0, -9.81, 0.0],
            ..Default::default()
        };
        let y0 = p[0].position[1];
        for _ in 0..100 {
            sph_step(&mut p, &config);
        }
        assert!(p[0].position[1] < y0);
    }
    #[test]
    fn test_lbm_outlet_pressure() {
        let mut lbm = LbmD2Q9::new(8, 8, 0.8);
        lbm.set_outlet_pressure(1.0);
        assert_eq!(lbm.cell_type[7], LbmCellType::Outlet);
    }
    #[test]
    fn test_mac_grid_pressure_project_no_panic() {
        let mut g = MacGrid::new(4, 4, 4, 0.1);
        g.flags = vec![1u8; 64];
        g.pressure_project(1000.0, 0.01);
    }
    #[test]
    fn test_gpu_sph_density_summation_parallel() {
        let mut particles = vec![
            SphParticle::new([0.0, 0.0, 0.0], 0.01),
            SphParticle::new([0.05, 0.0, 0.0], 0.01),
        ];
        let config = SphConfig::default();
        gpu_sph_density_parallel(&mut particles, &config);
        assert!(particles[0].density > 0.0);
        assert!(particles[1].density > 0.0);
    }
    #[test]
    fn test_gpu_pressure_jacobi_single_iter() {
        let mut g = MacGrid::new(4, 4, 4, 0.1);
        g.flags = vec![1u8; 64];
        gpu_jacobi_pressure_solve(&mut g, 1000.0, 0.01, 3);
        assert_eq!(g.p.len(), 64);
    }
    #[test]
    fn test_gpu_neighbor_list_update() {
        let positions = vec![[0.0f64, 0.0, 0.0], [0.05, 0.0, 0.0], [0.3, 0.0, 0.0]];
        let nl = GpuNeighborList::build(&positions, 0.1, [0.5, 0.5, 0.5]);
        let neighbors = nl.neighbors_of(0);
        assert!(neighbors.contains(&1));
        assert!(!neighbors.contains(&2));
    }
    #[test]
    fn test_gpu_lbm_bgk_collision() {
        let mut lbm = LbmD2Q9::new(8, 8, 0.8);
        gpu_lbm_bgk_collide(&mut lbm);
        let rho = lbm.density(4, 4);
        assert!((rho - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_morton_sort_particles() {
        let mut particles = vec![
            SphParticle::new([0.9, 0.9, 0.9], 1.0),
            SphParticle::new([0.0, 0.0, 0.0], 1.0),
            SphParticle::new([0.5, 0.5, 0.5], 1.0),
        ];
        morton_sort_particles(&mut particles, [1.0, 1.0, 1.0]);
        assert!(particles[0].position[0] < 0.2);
    }
    #[test]
    fn test_gpu_particle_integrate_euler() {
        let mut particles = vec![SphParticle::new([0.0, 0.0, 0.0], 1.0)];
        particles[0].force = [1.0, 0.0, 0.0];
        particles[0].density = 1.0;
        gpu_particle_integrate_euler(&mut particles, 0.01);
        assert!(particles[0].velocity[0] > 0.0);
        assert!(particles[0].position[0] > 0.0);
    }
    #[test]
    fn test_gpu_particle_integrate_verlet() {
        let mut particles = vec![SphParticle::new([0.0, 1.0, 0.0], 1.0)];
        particles[0].force = [0.0, -9.81, 0.0];
        particles[0].density = 1.0;
        let dt = 0.01;
        let y0 = particles[0].position[1];
        gpu_particle_integrate_verlet(&mut particles, dt, dt);
        assert!(particles[0].position[1] < y0);
    }
    #[test]
    fn test_gpu_boundary_condition_box() {
        let mut particles = vec![SphParticle::new([1.5, 0.5, 0.5], 1.0)];
        particles[0].velocity = [1.0, 0.0, 0.0];
        let bounds = GpuBoundaryBox {
            min: [0.0; 3],
            max: [1.0; 3],
            restitution: 0.5,
        };
        gpu_apply_boundary_box(&mut particles, &bounds);
        assert!(particles[0].position[0] <= 1.0);
        assert!(particles[0].velocity[0] <= 0.0);
    }
    #[test]
    fn test_multi_gpu_domain_decomp_2_devices() {
        let positions: Vec<[f64; 3]> = (0..8).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let domains = MultiGpuDomain::decompose_x(&positions, 2);
        assert_eq!(domains.len(), 2);
        let total: usize = domains.iter().map(|d| d.particle_indices.len()).sum();
        assert_eq!(total, 8);
    }
    #[test]
    fn test_gpu_reduce_energy() {
        let particles = [
            SphParticle::new([0.0; 3], 1.0),
            SphParticle::new([0.0; 3], 1.0),
        ];
        let mut p0 = particles[0].clone();
        let mut p1 = particles[1].clone();
        p0.velocity = [2.0, 0.0, 0.0];
        p1.velocity = [0.0, 1.0, 0.0];
        let particles2 = vec![p0, p1];
        let ke = gpu_reduce_kinetic_energy(&particles2, 0.5);
        assert!((ke - 1.25).abs() < 1e-10);
    }
    #[test]
    fn test_gpu_reduce_momentum() {
        let mut p = SphParticle::new([0.0; 3], 1.0);
        p.velocity = [3.0, 0.0, 0.0];
        let particles = vec![p];
        let mom = gpu_reduce_momentum(&particles, 2.0);
        assert!((mom[0] - 6.0).abs() < 1e-10);
    }
    #[test]
    fn test_gpu_advect_field() {
        let nx = 4usize;
        let ny = 4;
        let field: Vec<f64> = (0..nx * ny).map(|i| i as f64).collect();
        let vel = vec![[0.0f64; 2]; nx * ny];
        let result = gpu_advect_2d(&field, &vel, nx, ny, 0.1, 0.01);
        assert_eq!(result.len(), nx * ny);
    }
    #[test]
    fn test_gpu_pressure_poisson_staggered() {
        let mut p_grid = vec![0.0f64; 4 * 4];
        let div = vec![0.0f64; 4 * 4];
        gpu_pressure_poisson_jacobi_2d(&mut p_grid, &div, 4, 4, 0.1, 5);
        for &v in &p_grid {
            assert!(v.abs() < 1e-10);
        }
    }
    #[test]
    fn test_neighbor_list_cell_counts() {
        let positions: Vec<[f64; 3]> = (0..4).map(|i| [i as f64 * 0.02, 0.0, 0.0]).collect();
        let nl = GpuNeighborList::build(&positions, 0.1, [0.5, 0.5, 0.5]);
        let n0 = nl.neighbors_of(0);
        assert_eq!(n0.len(), 3);
    }
    #[test]
    fn test_morton_code_ordering() {
        let a = morton_encode_3d(0, 0, 0);
        let b = morton_encode_3d(1, 0, 0);
        let c = morton_encode_3d(0, 1, 0);
        assert!(a < b);
        assert!(a < c);
    }
    #[test]
    fn test_gpu_sph_verlet_energy_conservation_approx() {
        let mut particles = vec![SphParticle::new([0.0, 0.5, 0.0], 0.01)];
        let config = SphConfig {
            gravity: [0.0, 0.0, 0.0],
            ..Default::default()
        };
        sph_compute_density(&mut particles, &config);
        sph_compute_pressure(&mut particles, &config);
        let v_before = particles[0].velocity;
        gpu_particle_integrate_verlet(&mut particles, 0.001, 0.001);
        let v_after = particles[0].velocity;
        assert!((v_before[1] - v_after[1]).abs() < 1e-8);
    }
    #[test]
    fn test_domain_decomp_particles_go_to_correct_half() {
        let positions: Vec<[f64; 3]> = vec![
            [0.1, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.6, 0.0, 0.0],
            [0.8, 0.0, 0.0],
        ];
        let domains = MultiGpuDomain::decompose_x(&positions, 2);
        assert_eq!(domains[0].particle_indices.len(), 2);
        assert_eq!(domains[1].particle_indices.len(), 2);
    }
}
