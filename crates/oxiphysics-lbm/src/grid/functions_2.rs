//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
mod extended_grid_tests {

    use crate::grid::functions::{
        count_fluid_cells, enstrophy_2d, mean_density, mean_kinetic_energy, outflow_neumann_east,
        pressure_inlet_west, pressure_outlet_east, velocity_l2_norm_diff,
    };
    use crate::grid::grid_extended::{CellularGrid3D, TrtGrid2D};

    use crate::grid::types::{CellularGrid2D, ChannelFlow2D};
    use crate::lattice::{equilibrium_d2q9, macros_from_d2q9};
    #[test]
    fn test_cellular_grid2d_initial_mass() {
        let g = CellularGrid2D::new(6, 4, 1.5);
        let mass = g.total_mass();
        assert!((mass - 24.0).abs() < 1e-12, "mass = {mass}");
    }
    #[test]
    fn test_cellular_grid2d_initial_ke_zero() {
        let g = CellularGrid2D::new(4, 4, 1.5);
        assert!(g.kinetic_energy() < 1e-20, "ke = {}", g.kinetic_energy());
    }
    #[test]
    fn test_cellular_grid2d_idx_correct() {
        let g = CellularGrid2D::new(5, 3, 1.5);
        assert_eq!(g.idx(2, 1), 7);
    }
    #[test]
    fn test_cellular_grid2d_set_wall_marks_correctly() {
        let mut g = CellularGrid2D::new(4, 4, 1.5);
        g.set_wall(2, 3);
        assert!(g.wall[g.idx(2, 3)]);
        assert!(!g.wall[g.idx(1, 1)]);
    }
    #[test]
    fn test_cellular_grid2d_clear_wall() {
        let mut g = CellularGrid2D::new(4, 4, 1.5);
        g.set_wall(1, 1);
        g.clear_wall(1, 1);
        assert!(!g.wall[g.idx(1, 1)]);
    }
    #[test]
    fn test_cellular_grid2d_collide_conserves_mass() {
        let mut g = CellularGrid2D::new(6, 4, 1.5);
        g.pop[3][1] += 0.01;
        g.pop[3][3] -= 0.01;
        let mass_before = g.total_mass();
        g.collide();
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-13);
    }
    #[test]
    fn test_cellular_grid2d_stream_conserves_mass() {
        let mut g = CellularGrid2D::new(6, 4, 1.5);
        g.pop[5][1] += 0.05;
        g.pop[5][3] -= 0.05;
        let mass_before = g.total_mass();
        g.stream();
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-12);
    }
    #[test]
    fn test_cellular_grid2d_step_conserves_mass_no_walls() {
        let mut g = CellularGrid2D::new(8, 8, 1.2);
        let mass_before = g.total_mass();
        for _ in 0..20 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-10);
    }
    #[test]
    fn test_cellular_grid2d_step_with_walls_conserves_mass() {
        let mut g = CellularGrid2D::new(6, 6, 1.5);
        for x in 0..6 {
            g.set_wall(x, 0);
            g.set_wall(x, 5);
        }
        let mass_before = g.total_mass();
        for _ in 0..10 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-10);
    }
    #[test]
    fn test_cellular_grid2d_set_uniform_equilibrium() {
        let mut g = CellularGrid2D::new(4, 4, 1.5);
        g.set_uniform_equilibrium(1.2, 0.05, 0.0);
        for node in &g.pop {
            let (rho, ux, _) = macros_from_d2q9(node);
            assert!((rho - 1.2).abs() < 1e-12);
            assert!((ux - 0.05).abs() < 1e-12);
        }
    }
    #[test]
    fn test_cellular_grid2d_body_force_increases_momentum() {
        let mut g = CellularGrid2D::new(8, 4, 1.5);
        let jx_before = g.total_momentum_x();
        for _ in 0..20 {
            g.step_with_force(1e-4, 0.0);
        }
        let jx_after = g.total_momentum_x();
        assert!(jx_after > jx_before, "force should add x-momentum");
    }
    #[test]
    fn test_channel_flow2d_develops_positive_mean_ux() {
        let mut ch = ChannelFlow2D::new(4, 8, 1.0, 1e-4);
        for _ in 0..500 {
            ch.step();
        }
        let mean = ch.mean_ux();
        assert!(mean > 0.0, "Channel flow should be positive: {mean}");
    }
    #[test]
    fn test_channel_flow2d_poiseuille_umax_positive() {
        let ch = ChannelFlow2D::new(4, 8, 1.0, 1e-5);
        assert!(ch.poiseuille_umax() > 0.0);
    }
    #[test]
    fn test_cellular_grid3d_initial_mass() {
        let g = CellularGrid3D::new(4, 4, 4, 1.5);
        let mass = g.total_mass();
        assert!((mass - 64.0).abs() < 1e-11, "3D mass = {mass}");
    }
    #[test]
    fn test_cellular_grid3d_idx_correct() {
        let g = CellularGrid3D::new(4, 4, 4, 1.5);
        assert_eq!(g.idx(1, 2, 3), 3 * 16 + 2 * 4 + 1);
    }
    #[test]
    fn test_cellular_grid3d_step_conserves_mass() {
        let mut g = CellularGrid3D::new(4, 4, 4, 1.5);
        let mass_before = g.total_mass();
        for _ in 0..10 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-9);
    }
    #[test]
    fn test_cellular_grid3d_set_wall() {
        let mut g = CellularGrid3D::new(4, 4, 4, 1.5);
        g.set_wall(2, 2, 2);
        assert!(g.wall[g.idx(2, 2, 2)]);
    }
    #[test]
    fn test_trt_grid2d_initial_mass() {
        let g = TrtGrid2D::new_magic(6, 4, 1.5);
        let mass = g.total_mass();
        assert!((mass - 24.0).abs() < 1e-12, "TRT mass = {mass}");
    }
    #[test]
    fn test_trt_grid2d_step_conserves_mass() {
        let mut g = TrtGrid2D::new_magic(6, 6, 1.5);
        let mass_before = g.total_mass();
        for _ in 0..20 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-10);
    }
    #[test]
    fn test_trt_grid2d_magic_lambda() {
        let g = TrtGrid2D::new_magic(4, 4, 1.2);
        let tau_s: f64 = 1.0 / g.omega_sym;
        let tau_a: f64 = 1.0 / g.omega_anti;
        let lambda = (tau_s - 0.5) * (tau_a - 0.5);
        assert!((lambda - 3.0 / 16.0).abs() < 1e-12, "λ = {lambda}");
    }
    #[test]
    fn test_trt_grid2d_collide_conserves_mass() {
        let mut g = TrtGrid2D::new_magic(4, 4, 1.5);
        g.pop[3][1] += 0.02;
        g.pop[3][3] -= 0.02;
        let mass_before = g.total_mass();
        g.collide();
        let mass_after = g.total_mass();
        assert!((mass_before - mass_after).abs() < 1e-13);
    }
    #[test]
    fn test_mean_kinetic_energy_zero_at_rest() {
        let pop = vec![equilibrium_d2q9(1.0, 0.0, 0.0); 16];
        let ke = mean_kinetic_energy(&pop);
        assert!(ke < 1e-20, "ke at rest = {ke}");
    }
    #[test]
    fn test_mean_kinetic_energy_positive_with_flow() {
        let pop = vec![equilibrium_d2q9(1.0, 0.1, 0.0); 16];
        let ke = mean_kinetic_energy(&pop);
        assert!(ke > 0.0, "ke with flow = {ke}");
    }
    #[test]
    fn test_count_fluid_cells() {
        let mut g = CellularGrid2D::new(4, 4, 1.5);
        g.set_wall(0, 0);
        g.set_wall(1, 1);
        assert_eq!(count_fluid_cells(&g), 14);
    }
    #[test]
    fn test_mean_density_uniform() {
        let g = CellularGrid2D::new(4, 4, 1.5);
        let rho = mean_density(&g);
        assert!((rho - 1.0).abs() < 1e-13, "rho = {rho}");
    }
    #[test]
    fn test_velocity_l2_norm_diff_identical_is_zero() {
        let pop = vec![equilibrium_d2q9(1.0, 0.05, 0.0); 16];
        let diff = velocity_l2_norm_diff(&pop, &pop);
        assert!(diff < 1e-14, "diff = {diff}");
    }
    #[test]
    fn test_velocity_l2_norm_diff_nonzero_for_different() {
        let a = vec![equilibrium_d2q9(1.0, 0.1, 0.0); 4];
        let b = vec![equilibrium_d2q9(1.0, 0.0, 0.0); 4];
        let diff = velocity_l2_norm_diff(&a, &b);
        assert!(diff > 0.0, "diff = {diff}");
    }
    #[test]
    fn test_pressure_inlet_west_sets_density() {
        let nx = 4usize;
        let ny = 2usize;
        let mut pop = vec![equilibrium_d2q9(1.0, 0.0, 0.0); nx * ny];
        pressure_inlet_west(&mut pop, nx, ny, 1.05);
        for y in 0..ny {
            let (rho, _, _) = macros_from_d2q9(&pop[y * nx]);
            assert!((rho - 1.05).abs() < 1e-12, "inlet rho = {rho}");
        }
    }
    #[test]
    fn test_pressure_outlet_east_sets_density() {
        let nx = 4usize;
        let ny = 2usize;
        let mut pop = vec![equilibrium_d2q9(1.0, 0.05, 0.0); nx * ny];
        pressure_outlet_east(&mut pop, nx, ny, 0.98);
        for y in 0..ny {
            let (rho, _, _) = macros_from_d2q9(&pop[y * nx + (nx - 1)]);
            assert!((rho - 0.98).abs() < 1e-12, "outlet rho = {rho}");
        }
    }
    #[test]
    fn test_outflow_neumann_east_copies_neighbor() {
        let nx = 4usize;
        let ny = 2usize;
        let mut pop = vec![equilibrium_d2q9(1.0, 0.0, 0.0); nx * ny];
        for y in 0..ny {
            pop[y * nx + (nx - 2)][0] = 99.0;
        }
        outflow_neumann_east(&mut pop, nx, ny);
        for y in 0..ny {
            assert!((pop[y * nx + (nx - 1)][0] - 99.0).abs() < 1e-14);
        }
    }
    #[test]
    fn test_enstrophy_2d_zero_at_rest() {
        let pop = vec![equilibrium_d2q9(1.0, 0.0, 0.0); 16];
        let enst = enstrophy_2d(&pop, 4, 4);
        assert!(enst < 1e-20, "enstrophy at rest = {enst}");
    }
    #[test]
    fn test_enstrophy_2d_positive_with_shear() {
        let nx = 4usize;
        let ny = 4usize;
        let mut pop = vec![equilibrium_d2q9(1.0, 0.0, 0.0); nx * ny];
        for x in 0..nx {
            for y in 0..ny {
                let u = if y < ny / 2 { 0.05 } else { -0.05 };
                pop[y * nx + x] = equilibrium_d2q9(1.0, u, 0.0);
            }
        }
        let enst = enstrophy_2d(&pop, nx, ny);
        assert!(enst > 0.0, "enstrophy with shear = {enst}");
    }
    #[test]
    fn test_ux_profile_all_zero_at_rest() {
        let g = CellularGrid2D::new(4, 6, 1.5);
        let profile = g.ux_profile(2);
        for (y, &ux) in profile.iter().enumerate() {
            assert!(ux.abs() < 1e-14, "ux[{y}] = {ux}");
        }
    }
    #[test]
    fn test_max_speed_zero_at_rest() {
        let g = CellularGrid2D::new(4, 4, 1.5);
        assert!(g.max_speed() < 1e-14);
    }
    #[test]
    fn test_max_speed_positive_with_flow() {
        let mut g = CellularGrid2D::new(4, 4, 1.5);
        g.set_uniform_equilibrium(1.0, 0.1, 0.0);
        assert!(g.max_speed() > 0.0);
    }
    #[test]
    fn test_cellular_grid3d_macros_at_rest() {
        let g = CellularGrid3D::new(4, 4, 4, 1.5);
        let (rho, ux, uy, uz) = g.macros_at(2, 2, 2);
        assert!((rho - 1.0).abs() < 1e-13);
        assert!(ux.abs() < 1e-13);
        assert!(uy.abs() < 1e-13);
        assert!(uz.abs() < 1e-13);
    }
    #[test]
    fn test_cellular_grid3d_max_speed_zero_at_rest() {
        let g = CellularGrid3D::new(4, 4, 4, 1.5);
        assert!(g.max_speed() < 1e-14);
    }
    #[test]
    fn test_flagged_step_mass_conservation() {
        use crate::grid::grid_extended::FlaggedGrid3D;
        let nx = 6usize;
        let ny = 6usize;
        let nz = 6usize;
        let omega = 1.0;
        let rho0 = 1.0;
        let mut g = FlaggedGrid3D::new(nx, ny, nz, omega, rho0);
        // Perturb a few cells with a non-zero uniform velocity to create dynamics.
        g.initialize_uniform(rho0, 0.05, 0.01, 0.0);
        let mass_before: f64 = g.total_mass();
        for _ in 0..10 {
            g.step();
        }
        let mass_after: f64 = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "mass not conserved: before={mass_before}, after={mass_after}, diff={}",
            (mass_before - mass_after).abs()
        );
    }
    #[test]
    fn test_flagged_step_solid_bounce_back() {
        use crate::grid::grid_extended::{FlaggedGrid3D, NodeFlag};
        use crate::lattice::D3Q19_VELOCITIES;
        let nx = 6usize;
        let ny = 4usize;
        let nz = 4usize;
        let omega = 1.0;
        let rho0 = 1.0;
        let mut g = FlaggedGrid3D::new(nx, ny, nz, omega, rho0);
        // Add a solid wall layer at x=0.
        for z in 0..nz {
            for y in 0..ny {
                g.set_wall(0, y, z);
            }
        }
        // Drive flow in +x direction on fluid nodes only via initialize_uniform
        // (it skips non-Fluid nodes).
        g.initialize_uniform(rho0, 0.05, 0.0, 0.0);
        // Run many steps so the wall bounce-back has fully acted.
        for _ in 0..200 {
            g.step();
        }
        // Verify that wall flags are preserved.
        for z in 0..nz {
            for y in 0..ny {
                let widx = g.cell_idx(0, y, z);
                assert_eq!(
                    g.flags[widx],
                    NodeFlag::Wall,
                    "flag at wall (0,{y},{z}) changed"
                );
            }
        }
        // The net x-momentum stored at a wall node must be near zero because
        // full-way bounce-back exactly reverses all x-velocity contributions.
        // After streaming, opposite populations pair up so the sum of
        // f[a] * cx[a] cancels to within floating-point noise.
        for z in 0..nz {
            for y in 0..ny {
                let base = g.cell_idx(0, y, z) * 19;
                let mx: f64 = g.f[base..base + 19]
                    .iter()
                    .zip(D3Q19_VELOCITIES.iter())
                    .map(|(&fa, vel)| fa * vel[0] as f64)
                    .sum();
                assert!(
                    mx.abs() < 1e-6,
                    "wall node (0,{y},{z}) net x-momentum = {mx}, expected ~0"
                );
            }
        }
    }
    #[test]
    fn test_flagged_step_reduces_to_full_grid() {
        use crate::grid::grid_extended::{FlaggedGrid3D, FullGrid3D};
        let nx = 4usize;
        let ny = 4usize;
        let nz = 4usize;
        let omega = 1.2;
        let rho0 = 1.0;
        let ux0 = 0.05;
        let uy0 = 0.02;
        let uz0 = 0.0;
        // Build both grids with same initial conditions.
        let mut fg = FlaggedGrid3D::new(nx, ny, nz, omega, rho0);
        fg.initialize_uniform(rho0, ux0, uy0, uz0);
        let mut full = FullGrid3D::new(nx, ny, nz, omega, rho0);
        full.initialize_uniform(rho0, ux0, uy0, uz0);
        // Run 5 steps on each.
        for _ in 0..5 {
            fg.step();
            full.step();
        }
        // Per-cell f values must agree within floating-point tolerance.
        let n = nx * ny * nz * 19;
        let max_diff =
            fg.f.iter()
                .zip(full.f.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
        assert!(
            max_diff < 1e-10,
            "FlaggedGrid3D diverged from FullGrid3D: max_diff={max_diff} over {n} entries"
        );
    }
}
