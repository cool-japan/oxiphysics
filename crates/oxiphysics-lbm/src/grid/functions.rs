//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::lattice::{CS2, equilibrium_d2q9, macros_from_d2q9};

use super::grid_extended::{FullGrid3D, NodeFlag};
use super::types::{CellularGrid2D, FlaggedGrid2D};

/// Compute equilibrium distribution for a single direction (2D).
///
/// feq_i = w_i * rho * (1 + (e_i · u)/cs² + (e_i · u)²/(2 cs⁴) − u·u/(2 cs²))
pub fn equilibrium_2d(w: f64, rho: f64, ux: f64, uy: f64, cx: f64, cy: f64) -> f64 {
    let eu = cx * ux + cy * uy;
    let u_sq = ux * ux + uy * uy;
    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2))
}
/// Compute equilibrium distribution for a single direction (3D).
///
/// feq_i = w_i * rho * (1 + (e_i · u)/cs² + (e_i · u)²/(2 cs⁴) − u·u/(2 cs²))
pub fn equilibrium_3d(
    w: f64,
    rho: f64,
    ux: f64,
    uy: f64,
    uz: f64,
    cx: f64,
    cy: f64,
    cz: f64,
) -> f64 {
    let eu = cx * ux + cy * uy + cz * uz;
    let u_sq = ux * ux + uy * uy + uz * uz;
    w * rho * (1.0 + eu / CS2 + eu * eu / (2.0 * CS2 * CS2) - u_sq / (2.0 * CS2))
}
/// D2Q9 weights (local copy for self-contained module).
pub(crate) const W9: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
/// D2Q9 velocities (local copy).
pub(crate) const C9: [[i32; 2]; 9] = [
    [0, 0],
    [1, 0],
    [0, 1],
    [-1, 0],
    [0, -1],
    [1, 1],
    [-1, 1],
    [-1, -1],
    [1, -1],
];
/// D2Q9 opposite directions.
pub(crate) const OPP9: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
#[cfg(test)]
mod tests {

    use crate::grid::grid_extended::*;
    use crate::grid::types::*;
    #[test]
    fn test_flagged_grid_new_uniform_density() {
        let g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        for y in 0..4 {
            for x in 0..4 {
                let rho = g.density_at(x, y);
                assert!(
                    (rho - 1.0).abs() < 1e-14,
                    "density at ({x},{y}) = {rho}, expected 1.0"
                );
            }
        }
    }
    #[test]
    fn test_flagged_grid_total_mass() {
        let g = FlaggedGrid2D::new(5, 5, 1.0, 2.0);
        let mass = g.total_mass();
        assert!(
            (mass - 50.0).abs() < 1e-12,
            "total mass = {mass}, expected 50.0"
        );
    }
    #[test]
    fn test_flagged_grid_set_wall() {
        let mut g = FlaggedGrid2D::new(5, 5, 1.0, 1.0);
        g.set_wall(0, 0);
        assert_eq!(g.flags[0], NodeFlag::Wall);
    }
    #[test]
    fn test_flagged_grid_set_inlet() {
        let mut g = FlaggedGrid2D::new(5, 5, 1.0, 1.0);
        g.set_inlet(0, 2, 1.0, 0.05, 0.0);
        match g.flags[2 * 5] {
            NodeFlag::Inlet { .. } => {}
            other => panic!("Expected Inlet, got {other:?}"),
        }
    }
    #[test]
    fn test_flagged_grid_set_outlet() {
        let mut g = FlaggedGrid2D::new(5, 5, 1.0, 1.0);
        g.set_outlet(4, 2);
        assert_eq!(g.flags[2 * 5 + 4], NodeFlag::Outlet);
    }
    #[test]
    fn test_flagged_grid_step_mass_conservation() {
        let nx = 8;
        let ny = 8;
        let mut g = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..20 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Step did not conserve mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_flagged_grid_wall_bounce_back() {
        let mut g = FlaggedGrid2D::new(5, 5, 1.0, 1.0);
        g.set_wall(2, 2);
        let idx = 2 * 5 + 2;
        g.f[idx * 9 + 1] = 0.5;
        g.f[idx * 9 + 3] = 0.1;
        g.step();
        let rho = g.density_at(2, 2);
        assert!(rho > 0.0, "Wall cell density should be positive");
    }
    #[test]
    fn test_flagged_grid_velocity_at_rest() {
        let g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        for y in 0..4 {
            for x in 0..4 {
                let vel = g.velocity_at(x, y);
                assert!(
                    vel[0].abs() < 1e-14 && vel[1].abs() < 1e-14,
                    "velocity at ({x},{y}) should be zero"
                );
            }
        }
    }
    #[test]
    fn test_flagged_grid3d_stub() {
        let mut g = FlaggedGrid3D::new(3, 3, 3, 1.0, 1.0);
        g.step();
        assert_eq!(g.flags.len(), 27);
    }
    #[test]
    fn test_node_flag_equality() {
        assert_eq!(NodeFlag::Fluid, NodeFlag::Fluid);
        assert_ne!(NodeFlag::Wall, NodeFlag::Fluid);
        assert_ne!(NodeFlag::Outlet, NodeFlag::Symmetry);
    }
}
/// D3Q19 weights (local copy for the full 3D step).
pub(crate) const W19: [f64; 19] = crate::lattice::D3Q19_WEIGHTS;
/// D3Q19 velocity vectors.
pub(crate) const C19: [[i32; 3]; 19] = crate::lattice::D3Q19_VELOCITIES;
/// D3Q19 opposite direction indices.
pub(crate) const OPP19: [usize; 19] = crate::lattice::D3Q19_OPPOSITES;
#[cfg(test)]
mod full_grid3d_tests {

    use crate::grid::grid_extended::*;

    #[test]
    fn test_full_grid3d_initial_density() {
        let g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..4 {
                    let rho = g.density_at(x, y, z);
                    assert!((rho - 1.0).abs() < 1e-13, "rho at ({x},{y},{z}) = {rho}");
                }
            }
        }
    }
    #[test]
    fn test_full_grid3d_total_mass() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 2.0);
        let mass = g.total_mass();
        assert!((mass - 54.0).abs() < 1e-10, "total mass = {mass}");
    }
    #[test]
    fn test_full_grid3d_step_mass_conservation() {
        let mut g = FullGrid3D::new(5, 5, 5, 1.0, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..10 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "mass before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_full_grid3d_wall_bounce_back() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        g.set_wall(1, 1, 1);
        assert_eq!(g.flags[g.cell_idx(1, 1, 1)], NodeFlag::Wall);
        g.step();
    }
    #[test]
    fn test_full_grid3d_velocity_at_rest() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    let vel = g.velocity_at(x, y, z);
                    for &v in &vel {
                        assert!(v.abs() < 1e-13, "velocity component should be zero");
                    }
                }
            }
        }
    }
    #[test]
    fn test_full_grid3d_checkpoint_save_load() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        let checkpoint = g.checkpoint_save();
        g.f[0] += 0.5;
        g.checkpoint_load(&checkpoint);
        assert!(
            (g.f[0] - checkpoint[0]).abs() < 1e-15,
            "checkpoint restore failed"
        );
    }
    #[test]
    fn test_full_grid3d_mass_conservation_error_zero() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let ck = g.checkpoint_save();
        let err = g.mass_conservation_error(&ck);
        assert!(err.abs() < 1e-15, "error = {err}");
    }
    #[test]
    fn test_full_grid3d_initialize_uniform() {
        let mut g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        g.initialize_uniform(1.2, 0.05, 0.0, 0.0);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    let rho = g.density_at(x, y, z);
                    assert!((rho - 1.2).abs() < 1e-10, "rho = {rho}");
                }
            }
        }
    }
    #[test]
    fn test_full_grid3d_set_inlet() {
        let mut g = FullGrid3D::new(5, 5, 5, 1.0, 1.0);
        g.set_inlet(0, 2, 2, 1.0, 0.05, 0.0);
        let idx = g.cell_idx(0, 2, 2);
        match g.flags[idx] {
            NodeFlag::Inlet { .. } => {}
            other => panic!("Expected Inlet, got {other:?}"),
        }
    }
    #[test]
    fn test_full_grid3d_max_velocity_at_rest() {
        let g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        let max_v = g.max_velocity_magnitude();
        assert!(max_v < 1e-13, "max_v at rest = {max_v}");
    }
}
#[cfg(test)]
mod cell_grid_tests {

    use crate::grid::types::*;

    #[test]
    fn test_cell_grid2d_init_rho_uniform() {
        let rho0 = 1.5;
        let g = CellGrid2D::new(4, 4, rho0);
        for (k, cell) in g.cells.iter().enumerate() {
            assert!(
                (cell.rho - rho0).abs() < 1e-13,
                "cell[{k}].rho = {}, expected {rho0}",
                cell.rho
            );
        }
    }
    #[test]
    fn test_cell_grid2d_step_mass_conservation() {
        let mut g = CellGrid2D::new(6, 6, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..20 {
            g.step(1.0);
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-9,
            "Mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_cell_grid2d_obstacle_bounce_back() {
        let mut g = CellGrid2D::new(5, 5, 1.0);
        g.cells[2 * 5 + 2].obstacle = true;
        g.apply_bounce_back();
        let rho: f64 = g.cells[2 * 5 + 2].f.iter().sum();
        assert!(rho >= 0.0, "Obstacle cell rho should be >= 0");
    }
    #[test]
    fn test_cell_grid2d_mean_velocity_at_rest() {
        let g = CellGrid2D::new(4, 4, 1.0);
        let (mux, muy) = g.mean_velocity();
        assert!(
            mux.abs() < 1e-14 && muy.abs() < 1e-14,
            "Mean velocity at rest should be zero"
        );
    }
    #[test]
    fn test_cell_grid3d_init_rho_uniform() {
        let rho0 = 2.0;
        let g = CellGrid3D::new(3, 3, 3, rho0);
        for (k, cell) in g.cells.iter().enumerate() {
            assert!(
                (cell.rho - rho0).abs() < 1e-12,
                "cell3d[{k}].rho = {}, expected {rho0}",
                cell.rho
            );
        }
    }
    #[test]
    fn test_cell_grid3d_step_mass_conservation() {
        let mut g = CellGrid3D::new(4, 4, 4, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..10 {
            g.step(1.0);
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "3D mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_cell_grid3d_mean_velocity_at_rest() {
        let g = CellGrid3D::new(3, 3, 3, 1.0);
        let mean = g.mean_velocity();
        for &v in &mean {
            assert!(
                v.abs() < 1e-13,
                "3D mean velocity component should be zero, got {v}"
            );
        }
    }
}
#[cfg(test)]
mod lbm_grid2d_full_tests {

    use crate::grid::grid_extended::*;
    use crate::grid::types::*;

    #[test]
    fn test_lbm_grid2d_full_init_density() {
        let g = LbmGrid2DFull::new(8, 8, 1.0, 1.0);
        for y in 0..8 {
            for x in 0..8 {
                let rho = g.density_at(x, y);
                assert!((rho - 1.0).abs() < 1e-13, "rho at ({x},{y}) = {rho}");
            }
        }
    }
    #[test]
    fn test_lbm_grid2d_full_init_velocity_zero() {
        let g = LbmGrid2DFull::new(6, 6, 1.0, 1.0);
        for y in 0..6 {
            for x in 0..6 {
                let vel = g.velocity_at(x, y);
                assert!(
                    vel[0].abs() < 1e-14 && vel[1].abs() < 1e-14,
                    "velocity at ({x},{y}) = {:?}",
                    vel
                );
            }
        }
    }
    #[test]
    fn test_lbm_grid2d_full_mass_conservation() {
        let nx = 16;
        let ny = 8;
        let mut g = LbmGrid2DFull::new(nx, ny, 1.0, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..50 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-9,
            "Mass not conserved: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_lbm_grid2d_full_solid_flag() {
        let mut g = LbmGrid2DFull::new(8, 8, 1.0, 1.0);
        g.set_solid(0, 3);
        assert_eq!(g.node_type[g.cell_idx(0, 3)], BoundaryNodeType::Solid);
    }
    #[test]
    fn test_lbm_grid2d_full_inlet_flag() {
        let mut g = LbmGrid2DFull::new(8, 8, 1.0, 1.0);
        g.set_inlet(0, 4, 1.0, 0.05, 0.0);
        match g.node_type[g.cell_idx(0, 4)] {
            BoundaryNodeType::Inlet { ux, uy, rho } => {
                assert!((ux - 0.05).abs() < 1e-14);
                assert!(uy.abs() < 1e-14);
                assert!((rho - 1.0).abs() < 1e-14);
            }
            other => panic!("Expected Inlet, got {other:?}"),
        }
    }
    #[test]
    fn test_lbm_grid2d_full_outlet_flag() {
        let mut g = LbmGrid2DFull::new(8, 8, 1.0, 1.0);
        g.set_outlet(7, 4);
        assert_eq!(g.node_type[g.cell_idx(7, 4)], BoundaryNodeType::Outlet);
    }
    #[test]
    fn test_lbm_grid2d_full_body_force() {
        let nx = 8;
        let ny = 8;
        let mut g = LbmGrid2DFull::new(nx, ny, 1.0, 1.0);
        for _ in 0..20 {
            g.apply_body_force(1e-4, 0.0);
            g.step();
        }
        let n = nx * ny;
        let mean_ux: f64 = (0..n)
            .map(|idx| {
                let x = idx % nx;
                let y = idx / nx;
                g.velocity_at(x, y)[0]
            })
            .sum::<f64>()
            / n as f64;
        assert!(
            mean_ux > 0.0,
            "Body force should drive positive mean ux, got {mean_ux}"
        );
    }
    #[test]
    fn test_poiseuille_flow_convergence() {
        let nx = 4_usize;
        let ny = 20_usize;
        let omega = 1.0_f64;
        let fx = 1e-5_f64;
        let mut g = LbmGrid2DFull::new(nx, ny, omega, 1.0);
        for x in 0..nx {
            g.set_solid(x, 0);
            g.set_solid(x, ny - 1);
        }
        for _ in 0..5_000 {
            g.apply_body_force(fx, 0.0);
            g.step();
        }
        let tau = 1.0 / omega;
        let nu = (tau - 0.5) * (1.0_f64 / 3.0);
        let h = (ny - 2) as f64;
        let u_peak_analytic = fx * h * h / (8.0 * nu);
        let y_mid = ny / 2;
        let u_peak_lbm: f64 = (0..nx).map(|x| g.velocity_at(x, y_mid)[0]).sum::<f64>() / nx as f64;
        let rel_err = (u_peak_lbm - u_peak_analytic).abs() / u_peak_analytic.abs().max(1e-30);
        assert!(
            rel_err < 0.05,
            "Poiseuille: LBM={u_peak_lbm:.6e}, analytic={u_peak_analytic:.6e}, rel_err={rel_err:.4}"
        );
    }
}
#[cfg(test)]
mod flat_grid_tests {
    use super::*;
    use crate::grid::grid_extended::*;

    #[test]
    fn test_flat_lbm_grid2d_initial_density() {
        let g = FlatLbmGrid2D::new(8, 8, 1.0);
        for y in 0..8 {
            for x in 0..8 {
                let rho = g.density(x, y);
                assert!((rho - 1.0).abs() < 1e-13, "density at ({x},{y})={rho}");
            }
        }
    }
    #[test]
    fn test_flat_lbm_grid2d_initial_velocity_zero() {
        let g = FlatLbmGrid2D::new(6, 6, 1.0);
        for y in 0..6 {
            for x in 0..6 {
                let vel = g.velocity(x, y);
                assert!(
                    vel[0].abs() < 1e-14 && vel[1].abs() < 1e-14,
                    "velocity at ({x},{y}) = {vel:?}"
                );
            }
        }
    }
    #[test]
    fn test_flat_lbm_grid2d_mass_conservation() {
        let mut g = FlatLbmGrid2D::new(10, 10, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..30 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-9,
            "mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_flat_lbm_grid2d_set_inlet_velocity() {
        let mut g = FlatLbmGrid2D::new(8, 6, 1.0);
        g.set_inlet_velocity(0.05, 0.0);
        for y in 0..6 {
            let vel = g.velocity(0, y);
            assert!(
                vel[0].abs() > 0.0,
                "inlet velocity at y={y} should be nonzero"
            );
        }
    }
    #[test]
    fn test_flat_lbm_grid2d_apply_bounce_back() {
        let mut g = FlatLbmGrid2D::new(8, 8, 1.0);
        g.apply_bounce_back();
        assert!(
            g.total_mass() > 0.0,
            "mass should remain positive after bounce-back"
        );
    }
    #[test]
    fn test_flat_lbm_grid2d_inlet_and_step() {
        let mut g = FlatLbmGrid2D::new(8, 8, 1.0);
        g.set_inlet_velocity(0.05, 0.0);
        g.step();
        for y in 0..8 {
            for x in 0..8 {
                let rho = g.density(x, y);
                assert!(
                    rho > 0.0 && rho.is_finite(),
                    "density blow-up at ({x},{y})={rho}"
                );
            }
        }
    }
    #[test]
    fn test_flat_lbm_grid2d_poiseuille_density_conservation() {
        let nx = 4;
        let ny = 12;
        let tau = 1.0;
        let mut g = FlatLbmGrid2D::new(nx, ny, tau);
        let mass_before = g.total_mass();
        for _ in 0..50 {
            let omega = 1.0 / tau;
            let fx = 1e-5_f64;
            let n = nx * ny;
            for idx in 0..n {
                let base = idx * 9;
                for a in 0..9 {
                    let cx = C9[a][0] as f64;
                    g.f[base + a] += W9[a] * cx * fx / CS2;
                }
            }
            let _ = omega;
            g.step();
            g.apply_bounce_back();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() / mass_before < 0.01,
            "Poiseuille mass drift: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_flat_lbm_grid3d_initial_density() {
        let g = FlatLbmGrid3D::new(4, 4, 4, 1.0);
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..4 {
                    let rho = g.density(x, y, z);
                    assert!(
                        (rho - 1.0).abs() < 1e-13,
                        "3D density at ({x},{y},{z})={rho}"
                    );
                }
            }
        }
    }
    #[test]
    fn test_flat_lbm_grid3d_mass_conservation() {
        let mut g = FlatLbmGrid3D::new(4, 4, 4, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..10 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "3D mass: before={mass_before}, after={mass_after}"
        );
    }
}
/// Compute grid-wide statistics over distribution functions stored flat as
/// `f[cell_idx * Q + alpha]` with Q velocities and `n = nx * ny` cells.
///
/// Returns `(min_rho, max_rho, min_umag, max_umag)`.
pub fn grid_stats_2d(f: &[f64], nx: usize, ny: usize) -> (f64, f64, f64, f64) {
    let n = nx * ny;
    let mut min_rho = f64::INFINITY;
    let mut max_rho = f64::NEG_INFINITY;
    let mut min_umag = f64::INFINITY;
    let mut max_umag = f64::NEG_INFINITY;
    for idx in 0..n {
        let base = idx * 9;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for a in 0..9 {
            let fi = f[base + a];
            rho += fi;
            mx += fi * C9[a][0] as f64;
            my += fi * C9[a][1] as f64;
        }
        if rho < min_rho {
            min_rho = rho;
        }
        if rho > max_rho {
            max_rho = rho;
        }
        let umag = if rho.abs() > 1e-15 {
            let ux = mx / rho;
            let uy = my / rho;
            (ux * ux + uy * uy).sqrt()
        } else {
            0.0
        };
        if umag < min_umag {
            min_umag = umag;
        }
        if umag > max_umag {
            max_umag = umag;
        }
    }
    (min_rho, max_rho, min_umag, max_umag)
}
/// Compute grid-wide statistics for a 3D D3Q19 grid stored flat as
/// `f[cell_idx * 19 + alpha]` with `n = nx * ny * nz` cells.
///
/// Returns `(min_rho, max_rho, min_umag, max_umag)`.
pub fn grid_stats_3d(f: &[f64], nx: usize, ny: usize, nz: usize) -> (f64, f64, f64, f64) {
    let n = nx * ny * nz;
    let mut min_rho = f64::INFINITY;
    let mut max_rho = f64::NEG_INFINITY;
    let mut min_umag = f64::INFINITY;
    let mut max_umag = f64::NEG_INFINITY;
    for idx in 0..n {
        let base = idx * 19;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for a in 0..19 {
            let fi = f[base + a];
            rho += fi;
            mx += fi * C19[a][0] as f64;
            my += fi * C19[a][1] as f64;
            mz += fi * C19[a][2] as f64;
        }
        if rho < min_rho {
            min_rho = rho;
        }
        if rho > max_rho {
            max_rho = rho;
        }
        let umag = if rho.abs() > 1e-15 {
            let ux = mx / rho;
            let uy = my / rho;
            let uz = mz / rho;
            (ux * ux + uy * uy + uz * uz).sqrt()
        } else {
            0.0
        };
        if umag < min_umag {
            min_umag = umag;
        }
        if umag > max_umag {
            max_umag = umag;
        }
    }
    (min_rho, max_rho, min_umag, max_umag)
}
/// Save the state of a 2D flagged grid to a `Vec`f64` checkpoint buffer.
///
/// Format: `\[nx as f64, ny as f64, omega, f\[0\\], f\[1\], ...]`
pub fn checkpoint_save_2d(grid: &FlaggedGrid2D) -> Vec<f64> {
    let mut buf = Vec::with_capacity(3 + grid.f.len());
    buf.push(grid.nx as f64);
    buf.push(grid.ny as f64);
    buf.push(grid.omega);
    buf.extend_from_slice(&grid.f);
    buf
}
/// Restore a 2D flagged grid's distribution functions from a checkpoint buffer.
///
/// Only restores `f`; boundary flags and geometry are unchanged.
///
/// # Panics
///
/// Panics if the buffer size does not match the grid's expected distribution size.
pub fn checkpoint_restore_2d(grid: &mut FlaggedGrid2D, buf: &[f64]) {
    let expected = 3 + grid.nx * grid.ny * 9;
    assert_eq!(
        buf.len(),
        expected,
        "Checkpoint buffer size mismatch: got {}, expected {expected}",
        buf.len()
    );
    grid.f.copy_from_slice(&buf[3..]);
}
/// Save the state of a 3D full grid to a `Vec`f64` checkpoint buffer.
///
/// Format: `[nx as f64, ny as f64, nz as f64, omega, f[0\], ...]`
pub fn checkpoint_save_3d(grid: &FullGrid3D) -> Vec<f64> {
    let mut buf = Vec::with_capacity(4 + grid.f.len());
    buf.push(grid.nx as f64);
    buf.push(grid.ny as f64);
    buf.push(grid.nz as f64);
    buf.push(grid.omega);
    buf.extend_from_slice(&grid.f);
    buf
}
/// Restore a 3D full grid's distribution functions from a checkpoint buffer.
///
/// # Panics
///
/// Panics if the buffer size does not match the grid's expected distribution size.
pub fn checkpoint_restore_3d(grid: &mut FullGrid3D, buf: &[f64]) {
    let expected = 4 + grid.nx * grid.ny * grid.nz * 19;
    assert_eq!(
        buf.len(),
        expected,
        "3D checkpoint buffer size mismatch: got {}, expected {expected}",
        buf.len()
    );
    grid.f.copy_from_slice(&buf[4..]);
}
/// Extract the x-velocity profile along a vertical slice at column `x`.
///
/// Returns `Vec`f64` of length `ny` where entry `j` is `ux` at `(x, j)`.
pub fn extract_ux_profile_2d(grid: &FlaggedGrid2D, x: usize) -> Vec<f64> {
    let nx = grid.nx;
    let ny = grid.ny;
    let mut profile = Vec::with_capacity(ny);
    for y in 0..ny {
        let idx = y * nx + x;
        let base = idx * 9;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        for (a, c9a) in C9.iter().enumerate() {
            let fi = grid.f[base + a];
            rho += fi;
            mx += fi * c9a[0] as f64;
        }
        profile.push(if rho.abs() > 1e-15 { mx / rho } else { 0.0 });
    }
    profile
}
/// Extract the y-velocity profile along a horizontal slice at row `y`.
///
/// Returns `Vec`f64` of length `nx` where entry `j` is `uy` at `(j, y)`.
pub fn extract_uy_profile_2d(grid: &FlaggedGrid2D, y: usize) -> Vec<f64> {
    let nx = grid.nx;
    let mut profile = Vec::with_capacity(nx);
    for x in 0..nx {
        let idx = y * nx + x;
        let base = idx * 9;
        let mut rho = 0.0_f64;
        let mut my = 0.0_f64;
        for (a, c9a) in C9.iter().enumerate() {
            let fi = grid.f[base + a];
            rho += fi;
            my += fi * c9a[1] as f64;
        }
        profile.push(if rho.abs() > 1e-15 { my / rho } else { 0.0 });
    }
    profile
}
/// Compute total kinetic energy of the flow: `KE = 0.5 * sum_k rho_k * |u_k|^2`.
pub fn kinetic_energy_2d(grid: &FlaggedGrid2D) -> f64 {
    let n = grid.nx * grid.ny;
    let mut ke = 0.0_f64;
    for idx in 0..n {
        let base = idx * 9;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for (a, c9a) in C9.iter().enumerate() {
            let fi = grid.f[base + a];
            rho += fi;
            mx += fi * c9a[0] as f64;
            my += fi * c9a[1] as f64;
        }
        if rho.abs() > 1e-15 {
            let ux = mx / rho;
            let uy = my / rho;
            ke += 0.5 * rho * (ux * ux + uy * uy);
        }
    }
    ke
}
/// Compute total kinetic energy of a 3D FullGrid3D.
pub fn kinetic_energy_3d(grid: &FullGrid3D) -> f64 {
    let n = grid.nx * grid.ny * grid.nz;
    let mut ke = 0.0_f64;
    for idx in 0..n {
        let base = idx * 19;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        let mut mz = 0.0_f64;
        for (a, c19a) in C19.iter().enumerate() {
            let fi = grid.f[base + a];
            rho += fi;
            mx += fi * c19a[0] as f64;
            my += fi * c19a[1] as f64;
            mz += fi * c19a[2] as f64;
        }
        if rho.abs() > 1e-15 {
            let ux = mx / rho;
            let uy = my / rho;
            let uz = mz / rho;
            ke += 0.5 * rho * (ux * ux + uy * uy + uz * uz);
        }
    }
    ke
}
/// Compute maximum velocity magnitude across all fluid cells (D2Q9 grid).
pub fn max_velocity_magnitude_2d(grid: &FlaggedGrid2D) -> f64 {
    let n = grid.nx * grid.ny;
    let mut max_mag = 0.0_f64;
    for idx in 0..n {
        if grid.flags[idx] == NodeFlag::Wall {
            continue;
        }
        let base = idx * 9;
        let mut rho = 0.0_f64;
        let mut mx = 0.0_f64;
        let mut my = 0.0_f64;
        for (a, c9a) in C9.iter().enumerate() {
            let fi = grid.f[base + a];
            rho += fi;
            mx += fi * c9a[0] as f64;
            my += fi * c9a[1] as f64;
        }
        if rho.abs() > 1e-15 {
            let ux = mx / rho;
            let uy = my / rho;
            let mag = (ux * ux + uy * uy).sqrt();
            if mag > max_mag {
                max_mag = mag;
            }
        }
    }
    max_mag
}
/// Apply a uniform body force `[fx, fy]` to all fluid cells using the
/// simple momentum correction approach:
/// `f_i += w_i * (c_i · F) / cs^2 * dt`
///
/// This is the simplest (Ladd-style) body force scheme.
pub fn apply_body_force_2d(grid: &mut FlaggedGrid2D, fx: f64, fy: f64) {
    use crate::lattice::CS2;
    let n = grid.nx * grid.ny;
    for idx in 0..n {
        if grid.flags[idx] != NodeFlag::Fluid {
            continue;
        }
        let base = idx * 9;
        for a in 0..9 {
            let cx = C9[a][0] as f64;
            let cy = C9[a][1] as f64;
            let c_dot_f = cx * fx + cy * fy;
            grid.f[base + a] += W9[a] * c_dot_f / CS2;
        }
    }
}
/// Apply a uniform body force `[fx, fy, fz]` to all fluid cells (3D D3Q19).
pub fn apply_body_force_3d(grid: &mut FullGrid3D, fx: f64, fy: f64, fz: f64) {
    let n = grid.nx * grid.ny * grid.nz;
    for idx in 0..n {
        if grid.flags[idx] != NodeFlag::Fluid {
            continue;
        }
        let base = idx * 19;
        for a in 0..19 {
            let cx = C19[a][0] as f64;
            let cy = C19[a][1] as f64;
            let cz = C19[a][2] as f64;
            let c_dot_f = cx * fx + cy * fy + cz * fz;
            grid.f[base + a] += W19[a] * c_dot_f / CS2;
        }
    }
}
/// Compute the L2 norm of the density difference between two snapshots.
///
/// Useful for convergence monitoring: `||rho_new - rho_old||_2 / n`.
pub fn density_l2_diff_2d(f_new: &[f64], f_old: &[f64], nx: usize, ny: usize) -> f64 {
    let n = nx * ny;
    let mut sum_sq = 0.0_f64;
    for idx in 0..n {
        let base = idx * 9;
        let mut rho_new = 0.0_f64;
        let mut rho_old = 0.0_f64;
        for a in 0..9 {
            rho_new += f_new[base + a];
            rho_old += f_old[base + a];
        }
        let diff = rho_new - rho_old;
        sum_sq += diff * diff;
    }
    (sum_sq / n as f64).sqrt()
}
/// Compute L2 norm of velocity magnitude difference between two 2D snapshots.
pub fn velocity_l2_diff_2d(f_new: &[f64], f_old: &[f64], nx: usize, ny: usize) -> f64 {
    let n = nx * ny;
    let mut sum_sq = 0.0_f64;
    for idx in 0..n {
        let base = idx * 9;
        let (mut rn, mut mxn, mut myn) = (0.0_f64, 0.0_f64, 0.0_f64);
        let (mut ro, mut mxo, mut myo) = (0.0_f64, 0.0_f64, 0.0_f64);
        for a in 0..9 {
            let fn_ = f_new[base + a];
            let fo = f_old[base + a];
            rn += fn_;
            mxn += fn_ * C9[a][0] as f64;
            myn += fn_ * C9[a][1] as f64;
            ro += fo;
            mxo += fo * C9[a][0] as f64;
            myo += fo * C9[a][1] as f64;
        }
        let (uxn, uyn) = if rn.abs() > 1e-15 {
            (mxn / rn, myn / rn)
        } else {
            (0.0, 0.0)
        };
        let (uxo, uyo) = if ro.abs() > 1e-15 {
            (mxo / ro, myo / ro)
        } else {
            (0.0, 0.0)
        };
        let du = (uxn - uxo).hypot(uyn - uyo);
        sum_sq += du * du;
    }
    (sum_sq / n as f64).sqrt()
}
#[cfg(test)]
mod grid_physics_tests {
    use super::*;
    use crate::grid::types::*;
    #[test]
    fn test_grid_stats_2d_uniform_density() {
        let g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        let (min_rho, max_rho, _min_umag, _max_umag) = grid_stats_2d(&g.f, g.nx, g.ny);
        assert!((min_rho - 1.0).abs() < 1e-13, "min_rho={min_rho}");
        assert!((max_rho - 1.0).abs() < 1e-13, "max_rho={max_rho}");
    }
    #[test]
    fn test_grid_stats_2d_at_rest_zero_velocity() {
        let g = FlaggedGrid2D::new(6, 6, 1.0, 1.0);
        let (_min_rho, _max_rho, min_umag, max_umag) = grid_stats_2d(&g.f, g.nx, g.ny);
        assert!(min_umag.abs() < 1e-14, "min_umag at rest={min_umag}");
        assert!(max_umag.abs() < 1e-14, "max_umag at rest={max_umag}");
    }
    #[test]
    fn test_grid_stats_2d_nonzero_velocity() {
        let mut g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        g.set_inlet(1, 1, 1.0, 0.05, 0.0);
        let (min_rho, max_rho, _min_u, max_umag) = grid_stats_2d(&g.f, g.nx, g.ny);
        assert!(min_rho > 0.0, "min_rho should be positive");
        assert!(max_rho.is_finite(), "max_rho should be finite");
        assert!(
            max_umag >= 0.0,
            "max_umag should be non-negative: {max_umag}"
        );
    }
    #[test]
    fn test_grid_stats_3d_uniform_density() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let (min_rho, max_rho, _min_u, _max_u) = grid_stats_3d(&g.f, g.nx, g.ny, g.nz);
        assert!((min_rho - 1.0).abs() < 1e-12, "3D min_rho={min_rho}");
        assert!((max_rho - 1.0).abs() < 1e-12, "3D max_rho={max_rho}");
    }
    #[test]
    fn test_grid_stats_3d_at_rest_zero_velocity() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let (_min_rho, _max_rho, min_u, max_u) = grid_stats_3d(&g.f, g.nx, g.ny, g.nz);
        assert!(min_u.abs() < 1e-14, "3D min_u at rest={min_u}");
        assert!(max_u.abs() < 1e-14, "3D max_u at rest={max_u}");
    }
    #[test]
    fn test_checkpoint_save_restore_2d_roundtrip() {
        let nx = 5;
        let ny = 5;
        let mut g1 = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        g1.set_inlet(0, 2, 1.0, 0.05, 0.0);
        g1.set_outlet(nx - 1, 2);
        for _ in 0..5 {
            g1.step();
        }
        let buf = checkpoint_save_2d(&g1);
        let mut g2 = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        checkpoint_restore_2d(&mut g2, &buf);
        for i in 0..g1.f.len() {
            assert!(
                (g1.f[i] - g2.f[i]).abs() < 1e-14,
                "Distribution mismatch at f[{i}]: {} vs {}",
                g1.f[i],
                g2.f[i]
            );
        }
    }
    #[test]
    fn test_checkpoint_save_2d_header() {
        let nx = 4;
        let ny = 6;
        let omega = 1.2;
        let g = FlaggedGrid2D::new(nx, ny, omega, 1.0);
        let buf = checkpoint_save_2d(&g);
        assert!(
            (buf[0] - nx as f64).abs() < 1e-15,
            "buf[0] should be nx={nx}"
        );
        assert!(
            (buf[1] - ny as f64).abs() < 1e-15,
            "buf[1] should be ny={ny}"
        );
        assert!(
            (buf[2] - omega).abs() < 1e-15,
            "buf[2] should be omega={omega}"
        );
    }
    #[test]
    fn test_checkpoint_save_restore_3d_roundtrip() {
        let nx = 3;
        let ny = 3;
        let nz = 3;
        let mut g1 = FullGrid3D::new(nx, ny, nz, 1.0, 1.0);
        for _ in 0..3 {
            g1.step();
        }
        let buf = checkpoint_save_3d(&g1);
        let mut g2 = FullGrid3D::new(nx, ny, nz, 1.0, 1.0);
        checkpoint_restore_3d(&mut g2, &buf);
        for i in 0..g1.f.len() {
            assert!(
                (g1.f[i] - g2.f[i]).abs() < 1e-14,
                "3D distribution mismatch at f[{i}]: {} vs {}",
                g1.f[i],
                g2.f[i]
            );
        }
    }
    #[test]
    fn test_checkpoint_save_3d_header() {
        let nx = 4;
        let ny = 5;
        let nz = 3;
        let omega = 1.5;
        let g = FullGrid3D::new(nx, ny, nz, omega, 1.0);
        let buf = checkpoint_save_3d(&g);
        assert!((buf[0] - nx as f64).abs() < 1e-15, "buf[0]={}", buf[0]);
        assert!((buf[1] - ny as f64).abs() < 1e-15, "buf[1]={}", buf[1]);
        assert!((buf[2] - nz as f64).abs() < 1e-15, "buf[2]={}", buf[2]);
        assert!((buf[3] - omega).abs() < 1e-15, "buf[3]={}", buf[3]);
    }
    #[test]
    fn test_extract_ux_profile_at_rest() {
        let g = FlaggedGrid2D::new(8, 8, 1.0, 1.0);
        let profile = extract_ux_profile_2d(&g, 4);
        for (j, &ux) in profile.iter().enumerate() {
            assert!(
                ux.abs() < 1e-14,
                "ux profile at rest should be zero, got {} at y={j}",
                ux
            );
        }
    }
    #[test]
    fn test_extract_ux_profile_length() {
        let g = FlaggedGrid2D::new(8, 12, 1.0, 1.0);
        let profile = extract_ux_profile_2d(&g, 3);
        assert_eq!(profile.len(), 12, "ux profile length should equal ny");
    }
    #[test]
    fn test_extract_uy_profile_length() {
        let g = FlaggedGrid2D::new(10, 6, 1.0, 1.0);
        let profile = extract_uy_profile_2d(&g, 3);
        assert_eq!(profile.len(), 10, "uy profile length should equal nx");
    }
    #[test]
    fn test_extract_uy_profile_at_rest() {
        let g = FlaggedGrid2D::new(8, 8, 1.0, 1.0);
        let profile = extract_uy_profile_2d(&g, 4);
        for (j, &uy) in profile.iter().enumerate() {
            assert!(
                uy.abs() < 1e-14,
                "uy profile at rest should be zero, got {} at x={j}",
                uy
            );
        }
    }
    #[test]
    fn test_kinetic_energy_2d_at_rest_zero() {
        let g = FlaggedGrid2D::new(6, 6, 1.0, 1.0);
        let ke = kinetic_energy_2d(&g);
        assert!(ke.abs() < 1e-14, "KE at rest should be ~0: {ke}");
    }
    #[test]
    fn test_kinetic_energy_3d_at_rest_zero() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let ke = kinetic_energy_3d(&g);
        assert!(ke.abs() < 1e-14, "3D KE at rest should be ~0: {ke}");
    }
    #[test]
    fn test_kinetic_energy_2d_positive_after_inlet() {
        let mut g = FlaggedGrid2D::new(8, 8, 1.0, 1.0);
        g.set_inlet(0, 4, 1.0, 0.05, 0.0);
        for _ in 0..5 {
            g.step();
        }
        let ke = kinetic_energy_2d(&g);
        assert!(ke >= 0.0, "Kinetic energy should be non-negative: {ke}");
    }
    #[test]
    fn test_max_velocity_magnitude_2d_at_rest() {
        let g = FlaggedGrid2D::new(6, 6, 1.0, 1.0);
        let max_u = max_velocity_magnitude_2d(&g);
        assert!(
            max_u.abs() < 1e-14,
            "max velocity at rest should be ~0: {max_u}"
        );
    }
    #[test]
    fn test_max_velocity_magnitude_2d_after_step() {
        let mut g = FlaggedGrid2D::new(8, 8, 1.0, 1.0);
        g.set_inlet(0, 4, 1.0, 0.05, 0.0);
        g.step();
        let max_u = max_velocity_magnitude_2d(&g);
        assert!(max_u >= 0.0, "max velocity should be non-negative");
        assert!(max_u.is_finite(), "max velocity should be finite");
    }
    #[test]
    fn test_apply_body_force_2d_increases_momentum() {
        let mut g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        let fx = 1e-4_f64;
        let fy = 0.0;
        let n = g.nx * g.ny;
        let mx_before: f64 = (0..n)
            .map(|idx| {
                let base = idx * 9;
                (0..9).map(|a| g.f[base + a] * C9[a][0] as f64).sum::<f64>()
            })
            .sum();
        apply_body_force_2d(&mut g, fx, fy);
        let mx_after: f64 = (0..n)
            .map(|idx| {
                let base = idx * 9;
                (0..9).map(|a| g.f[base + a] * C9[a][0] as f64).sum::<f64>()
            })
            .sum();
        assert!(
            mx_after > mx_before,
            "Body force in x should increase x-momentum: before={mx_before}, after={mx_after}"
        );
    }
    #[test]
    fn test_apply_body_force_2d_conserves_mass() {
        let mut g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        let mass_before = g.total_mass();
        apply_body_force_2d(&mut g, 1e-4, 0.0);
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "Body force should conserve mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_apply_body_force_3d_conserves_mass() {
        let mut g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let mass_before = g.total_mass();
        apply_body_force_3d(&mut g, 1e-4, 0.0, 0.0);
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "3D body force mass: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_apply_body_force_3d_increases_x_momentum() {
        let mut g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let n = g.nx * g.ny * g.nz;
        let mx_before: f64 = (0..n)
            .map(|idx| {
                let base = idx * 19;
                (0..19)
                    .map(|a| g.f[base + a] * C19[a][0] as f64)
                    .sum::<f64>()
            })
            .sum();
        apply_body_force_3d(&mut g, 1e-4, 0.0, 0.0);
        let mx_after: f64 = (0..n)
            .map(|idx| {
                let base = idx * 19;
                (0..19)
                    .map(|a| g.f[base + a] * C19[a][0] as f64)
                    .sum::<f64>()
            })
            .sum();
        assert!(
            mx_after > mx_before,
            "3D body force should increase x-momentum: before={mx_before}, after={mx_after}"
        );
    }
    #[test]
    fn test_density_l2_diff_2d_identical_grids_zero() {
        let g = FlaggedGrid2D::new(6, 6, 1.0, 1.0);
        let diff = density_l2_diff_2d(&g.f, &g.f, g.nx, g.ny);
        assert!(
            diff.abs() < 1e-14,
            "L2 diff of identical grids should be 0: {diff}"
        );
    }
    #[test]
    fn test_density_l2_diff_2d_different_grids_positive() {
        let g1 = FlaggedGrid2D::new(6, 6, 1.0, 1.0);
        let g2 = FlaggedGrid2D::new(6, 6, 1.0, 1.5);
        let diff = density_l2_diff_2d(&g1.f, &g2.f, g1.nx, g1.ny);
        assert!(
            diff > 0.0,
            "L2 diff of different grids should be positive: {diff}"
        );
    }
    #[test]
    fn test_velocity_l2_diff_2d_identical_zero() {
        let g = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        let diff = velocity_l2_diff_2d(&g.f, &g.f, g.nx, g.ny);
        assert!(
            diff.abs() < 1e-14,
            "velocity L2 diff of identical grids should be 0: {diff}"
        );
    }
    #[test]
    fn test_poiseuille_body_force_develops_flow() {
        let nx = 4;
        let ny = 10;
        let omega = 1.2;
        let mut g = FlaggedGrid2D::new(nx, ny, omega, 1.0);
        for x in 0..nx {
            g.set_wall(x, 0);
            g.set_wall(x, ny - 1);
        }
        let fx = 1e-5_f64;
        for _ in 0..200 {
            apply_body_force_2d(&mut g, fx, 0.0);
            g.step();
        }
        let max_u = max_velocity_magnitude_2d(&g);
        assert!(
            max_u > 0.0,
            "Poiseuille flow should develop with body force"
        );
    }
    #[test]
    fn test_full_grid_3d_step_mass_conservation_with_body_force() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..10 {
            apply_body_force_3d(&mut g, 1e-5, 0.0, 0.0);
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "3D mass with body force: before={mass_before}, after={mass_after}"
        );
    }
    #[test]
    fn test_checkpoint_restore_identical_step_outcome() {
        let nx = 4;
        let ny = 4;
        let mut g1 = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        for _ in 0..3 {
            g1.step();
        }
        let buf = checkpoint_save_2d(&g1);
        let mut g2 = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        checkpoint_restore_2d(&mut g2, &buf);
        g1.step();
        g2.step();
        for i in 0..g1.f.len() {
            assert!(
                (g1.f[i] - g2.f[i]).abs() < 1e-13,
                "After restore+step, distributions should match at f[{i}]"
            );
        }
    }
    #[test]
    fn test_l2_convergence_decreases_with_body_force() {
        let nx = 4;
        let ny = 8;
        let omega = 1.2;
        let mut g = FlaggedGrid2D::new(nx, ny, omega, 1.0);
        for x in 0..nx {
            g.set_wall(x, 0);
            g.set_wall(x, ny - 1);
        }
        let fx = 1e-5_f64;
        let mut f_prev = g.f.clone();
        let mut _l2_early = 0.0_f64;
        let mut l2_late = 0.0_f64;
        for step in 0..500 {
            apply_body_force_2d(&mut g, fx, 0.0);
            g.step();
            if step == 50 {
                _l2_early = density_l2_diff_2d(&g.f, &f_prev, nx, ny);
                f_prev = g.f.clone();
            }
            if step == 490 {
                l2_late = density_l2_diff_2d(&g.f, &f_prev, nx, ny);
            }
        }
        assert!(
            l2_late < 0.01,
            "L2 density change near steady state should be small: {l2_late}"
        );
    }
}
#[cfg(test)]
mod grid_extended_tests {
    use super::*;
    use crate::grid::grid_extended::*;
    use crate::grid::types::*;
    use crate::lattice::LatticeType;
    #[test]
    fn test_lbm_grid2d_initial_density_uniform() {
        let g = LbmGrid2D::new(5, 4, LatticeType::D2Q9);
        for y in 0..4 {
            for x in 0..5 {
                let rho = g.density_at(x, y);
                assert!((rho - 1.0).abs() < 1e-14, "rho({x},{y})={rho}");
            }
        }
    }
    #[test]
    fn test_lbm_grid2d_initial_velocity_zero() {
        let g = LbmGrid2D::new(5, 4, LatticeType::D2Q9);
        for y in 0..4 {
            for x in 0..5 {
                let (ux, uy) = g.velocity_at(x, y);
                assert!(ux.abs() < 1e-15 && uy.abs() < 1e-15);
            }
        }
    }
    #[test]
    fn test_lbm_grid2d_set_equilibrium_recovers_values() {
        let mut g = LbmGrid2D::new(6, 6, LatticeType::D2Q9);
        g.set_equilibrium(2, 3, 1.2, 0.05, -0.03);
        g.compute_macroscopic();
        let (ux, uy) = g.velocity_at(2, 3);
        let rho = g.density_at(2, 3);
        assert!((rho - 1.2).abs() < 1e-13, "rho={rho}");
        assert!((ux - 0.05).abs() < 1e-13, "ux={ux}");
        assert!((uy + 0.03).abs() < 1e-13, "uy={uy}");
    }
    #[test]
    fn test_lbm_grid2d_total_density_after_set_equilibrium() {
        let nx = 4;
        let ny = 4;
        let mut g = LbmGrid2D::new(nx, ny, LatticeType::D2Q9);
        g.set_equilibrium(1, 1, 2.0, 0.0, 0.0);
        g.compute_macroscopic();
        let total = g.total_density();
        let expected = (nx * ny - 1) as f64 + 2.0;
        assert!((total - expected).abs() < 1e-13, "total rho={total}");
    }
    #[test]
    fn test_lbm_grid2d_idx_linearizes_correctly() {
        let g = LbmGrid2D::new(8, 6, LatticeType::D2Q9);
        assert_eq!(g.idx(0, 0), 0);
        assert_eq!(g.idx(1, 0), 1);
        assert_eq!(g.idx(0, 1), 8);
        assert_eq!(g.idx(7, 5), 5 * 8 + 7);
    }
    #[test]
    fn test_lbm_grid3d_initial_density_uniform() {
        let g = LbmGrid3D::new(3, 3, 3, LatticeType::D3Q19);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    let rho = g.density_at(x, y, z);
                    assert!((rho - 1.0).abs() < 1e-14, "rho({x},{y},{z})={rho}");
                }
            }
        }
    }
    #[test]
    fn test_lbm_grid3d_initial_velocity_zero() {
        let g = LbmGrid3D::new(3, 3, 3, LatticeType::D3Q19);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    let (ux, uy, uz) = g.velocity_at(x, y, z);
                    assert!(ux.abs() < 1e-15 && uy.abs() < 1e-15 && uz.abs() < 1e-15);
                }
            }
        }
    }
    #[test]
    fn test_lbm_grid3d_set_equilibrium_recovers_values() {
        let mut g = LbmGrid3D::new(4, 4, 4, LatticeType::D3Q19);
        g.set_equilibrium(1, 2, 3, 1.1, 0.04, -0.02, 0.01);
        g.compute_macroscopic();
        let (ux, uy, uz) = g.velocity_at(1, 2, 3);
        let rho = g.density_at(1, 2, 3);
        assert!((rho - 1.1).abs() < 1e-12, "rho={rho}");
        assert!((ux - 0.04).abs() < 1e-12, "ux={ux}");
        assert!((uy + 0.02).abs() < 1e-12, "uy={uy}");
        assert!((uz - 0.01).abs() < 1e-12, "uz={uz}");
    }
    #[test]
    fn test_lbm_grid3d_idx_linearizes_correctly() {
        let g = LbmGrid3D::new(4, 5, 6, LatticeType::D3Q19);
        assert_eq!(g.idx(0, 0, 0), 0);
        assert_eq!(g.idx(1, 0, 0), 1);
        assert_eq!(g.idx(0, 1, 0), 4);
        assert_eq!(g.idx(0, 0, 1), 4 * 5);
        assert_eq!(g.idx(3, 4, 5), 5 * 4 * 5 + 4 * 4 + 3);
    }
    #[test]
    fn test_lbm_grid3d_total_density() {
        let g = LbmGrid3D::new(3, 3, 3, LatticeType::D3Q19);
        let total = g.total_density();
        assert!((total - 27.0).abs() < 1e-12, "total 3D rho={total}");
    }
    #[test]
    fn test_flagged_grid2d_wall_stays_at_equilibrium_mass() {
        let nx = 6;
        let ny = 6;
        let mut g = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        for x in 0..nx {
            g.set_wall(x, 0);
            g.set_wall(x, ny - 1);
        }
        for y in 0..ny {
            g.set_wall(0, y);
            g.set_wall(nx - 1, y);
        }
        let mass_before = g.total_mass();
        for _ in 0..10 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Wall BC mass: {mass_before} vs {mass_after}"
        );
    }
    #[test]
    fn test_flagged_grid2d_set_outlet_flag() {
        let mut g = FlaggedGrid2D::new(8, 8, 1.0, 1.0);
        g.set_outlet(7, 3);
        let idx = g.cell_idx(7, 3);
        assert_eq!(g.flags[idx], NodeFlag::Outlet);
    }
    #[test]
    fn test_flagged_grid2d_inlet_density_after_step() {
        let mut g = FlaggedGrid2D::new(8, 4, 1.0, 1.0);
        g.set_inlet(0, 2, 1.0, 0.05, 0.0);
        g.step();
        let rho = g.density_at(0, 2);
        assert!(rho > 0.0, "Inlet cell should have positive density: {rho}");
    }
    #[test]
    fn test_full_grid3d_initial_uniform_density() {
        let g = FullGrid3D::new(4, 4, 4, 1.0, 1.5);
        let mass = (0..4 * 4 * 4)
            .map(|idx| {
                let base = idx * 19;
                g.f[base..base + 19].iter().sum::<f64>()
            })
            .sum::<f64>();
        assert!(
            (mass - 1.5 * 64.0).abs() < 1e-10,
            "FullGrid3D mass = {mass}"
        );
    }
    #[test]
    fn test_full_grid3d_step_runs_without_panic() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        for _ in 0..5 {
            g.step();
        }
    }
    #[test]
    fn test_full_grid3d_wall_set() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        g.set_wall(1, 2, 3);
        let idx = g.cell_idx(1, 2, 3);
        assert_eq!(g.flags[idx], NodeFlag::Wall);
    }
    #[test]
    fn test_full_grid3d_density_at() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 2.0);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    let rho = g.density_at(x, y, z);
                    assert!(
                        (rho - 2.0).abs() < 1e-12,
                        "FullGrid3D rho({x},{y},{z})={rho}"
                    );
                }
            }
        }
    }
    #[test]
    fn test_full_grid3d_velocity_at_zero_initially() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        for z in 0..3 {
            for y in 0..3 {
                for x in 0..3 {
                    let v = g.velocity_at(x, y, z);
                    assert!(v[0].abs() < 1e-14 && v[1].abs() < 1e-14 && v[2].abs() < 1e-14);
                }
            }
        }
    }
    #[test]
    fn test_full_grid3d_total_mass() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let mass = g.total_mass();
        assert!(
            (mass - 27.0).abs() < 1e-10,
            "FullGrid3D total mass = {mass}"
        );
    }
    #[test]
    fn test_full_grid3d_max_velocity_zero_at_rest() {
        let g = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        let v_max = g.max_velocity_magnitude();
        assert!(v_max < 1e-14, "Max velocity at rest = {v_max}");
    }
    #[test]
    fn test_checkpoint_save_restore_2d_identical() {
        let mut g1 = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        g1.f[0] += 0.1;
        let buf = checkpoint_save_2d(&g1);
        let mut g2 = FlaggedGrid2D::new(4, 4, 1.0, 1.0);
        checkpoint_restore_2d(&mut g2, &buf);
        for i in 0..g1.f.len() {
            assert!(
                (g1.f[i] - g2.f[i]).abs() < 1e-15,
                "Checkpoint mismatch at [{i}]"
            );
        }
    }
    #[test]
    fn test_checkpoint_save_restore_3d_identical() {
        let mut g1 = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        g1.f[5] += 0.2;
        let buf = checkpoint_save_3d(&g1);
        let mut g2 = FullGrid3D::new(3, 3, 3, 1.0, 1.0);
        checkpoint_restore_3d(&mut g2, &buf);
        for i in 0..g1.f.len() {
            assert!(
                (g1.f[i] - g2.f[i]).abs() < 1e-15,
                "3D Checkpoint mismatch at [{i}]"
            );
        }
    }
    #[test]
    fn test_apply_body_force_2d_increases_x_momentum() {
        let nx = 4;
        let ny = 4;
        let mut g = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        let mom_before: f64 = (0..nx * ny)
            .map(|idx| {
                let base = idx * 9;
                (0..9).map(|a| g.f[base + a] * C9[a][0] as f64).sum::<f64>()
            })
            .sum();
        apply_body_force_2d(&mut g, 1e-4, 0.0);
        let mom_after: f64 = (0..nx * ny)
            .map(|idx| {
                let base = idx * 9;
                (0..9).map(|a| g.f[base + a] * C9[a][0] as f64).sum::<f64>()
            })
            .sum();
        assert!(
            mom_after > mom_before,
            "Body force should increase x-momentum: {mom_before} -> {mom_after}"
        );
    }
    #[test]
    fn test_apply_body_force_3d_does_not_change_mass() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.0, 1.0);
        let mass_before = g.total_mass();
        apply_body_force_3d(&mut g, 1e-4, 0.0, 0.0);
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "Body force 3D changed mass: {mass_before} vs {mass_after}"
        );
    }
    #[test]
    fn test_max_velocity_magnitude_2d_zero_at_rest() {
        let g = FlaggedGrid2D::new(5, 5, 1.0, 1.0);
        let v_max = max_velocity_magnitude_2d(&g);
        assert!(v_max < 1e-14, "Max velocity at rest = {v_max}");
    }
    #[test]
    fn test_density_l2_diff_2d_zero_for_identical() {
        let nx = 4;
        let ny = 4;
        let g = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        let diff = density_l2_diff_2d(&g.f, &g.f, nx, ny);
        assert!(diff < 1e-15, "L2 diff of identical = {diff}");
    }
    #[test]
    fn test_density_l2_diff_2d_positive_for_different() {
        let nx = 4;
        let ny = 4;
        let g1 = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        let g2 = FlaggedGrid2D::new(nx, ny, 1.0, 2.0);
        let diff = density_l2_diff_2d(&g1.f, &g2.f, nx, ny);
        assert!(
            diff > 0.0,
            "L2 diff should be positive for different grids: {diff}"
        );
    }
    #[test]
    fn test_node_flag_fluid_is_default() {
        let flags = vec![NodeFlag::Fluid; 10];
        for f in &flags {
            assert_eq!(*f, NodeFlag::Fluid);
        }
    }
    #[test]
    fn test_node_flag_wall_inequality_with_fluid() {
        assert_ne!(NodeFlag::Wall, NodeFlag::Fluid);
        assert_ne!(NodeFlag::Wall, NodeFlag::Outlet);
    }
    #[test]
    fn test_equilibrium_2d_positive_weight_result() {
        let val = equilibrium_2d(4.0 / 9.0, 1.0, 0.0, 0.0, 0.0, 0.0);
        assert!((val - 4.0 / 9.0).abs() < 1e-14, "eq 2d at rest = {val}");
    }
    #[test]
    fn test_equilibrium_3d_positive_weight_result() {
        let val = equilibrium_3d(1.0 / 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!((val - 1.0 / 3.0).abs() < 1e-14, "eq 3d at rest = {val}");
    }
    #[test]
    fn test_flagged_grid3d_step_does_not_panic() {
        let mut g = FlaggedGrid3D::new(3, 3, 3, 1.0, 1.0);
        for _ in 0..5 {
            g.step();
        }
    }
    #[test]
    fn test_flagged_grid3d_initial_mass() {
        let nx = 3;
        let ny = 3;
        let nz = 3;
        let rho0 = 1.5;
        let g = FlaggedGrid3D::new(nx, ny, nz, 1.0, rho0);
        let mass: f64 = (0..nx * ny * nz)
            .map(|idx| g.f[idx * 19..idx * 19 + 19].iter().sum::<f64>())
            .sum();
        let expected = rho0 * (nx * ny * nz) as f64;
        assert!(
            (mass - expected).abs() < 1e-10,
            "FlaggedGrid3D mass={mass}, expected={expected}"
        );
    }
    #[test]
    fn test_lbm_grid3d_d3q27_initial_density() {
        let g = LbmGrid3D::new(3, 3, 3, LatticeType::D3Q27);
        let total = g.total_density();
        assert!((total - 27.0).abs() < 1e-12, "D3Q27 total density={total}");
    }
    #[test]
    fn test_full_grid3d_periodic_steps_mass_conservation() {
        let mut g = FullGrid3D::new(4, 4, 4, 1.2, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..20 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-8,
            "3D periodic mass: {mass_before} vs {mass_after}"
        );
    }
    #[test]
    fn test_flagged_grid2d_omega_1_8_conserves_mass() {
        let nx = 6;
        let ny = 6;
        let mut g = FlaggedGrid2D::new(nx, ny, 1.8, 1.0);
        let mass_before = g.total_mass();
        for _ in 0..30 {
            g.step();
        }
        let mass_after = g.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-10,
            "omega=1.8 mass: {mass_before} vs {mass_after}"
        );
    }
    #[test]
    fn test_body_force_produces_flow_2d() {
        let nx = 8;
        let ny = 6;
        let mut g = FlaggedGrid2D::new(nx, ny, 1.0, 1.0);
        for x in 0..nx {
            g.set_wall(x, 0);
            g.set_wall(x, ny - 1);
        }
        for _ in 0..100 {
            apply_body_force_2d(&mut g, 1e-5, 0.0);
            g.step();
        }
        let v_max = max_velocity_magnitude_2d(&g);
        assert!(v_max > 0.0, "Body force should produce flow: v_max={v_max}");
    }
}
/// Compute L2 norm of velocity magnitude difference between two population arrays.
pub fn velocity_l2_norm_diff(a: &[[f64; 9]], b: &[[f64; 9]]) -> f64 {
    let mut sum = 0.0f64;
    for (na, nb) in a.iter().zip(b.iter()) {
        let (_, ua_x, ua_y) = macros_from_d2q9(na);
        let (_, ub_x, ub_y) = macros_from_d2q9(nb);
        let dx = ua_x - ub_x;
        let dy = ua_y - ub_y;
        sum += dx * dx + dy * dy;
    }
    sum.sqrt()
}
/// Compute the mean kinetic energy over a population array.
pub fn mean_kinetic_energy(pop: &[[f64; 9]]) -> f64 {
    let ke: f64 = pop
        .iter()
        .map(|node| {
            let (rho, ux, uy) = macros_from_d2q9(node);
            0.5 * rho * (ux * ux + uy * uy)
        })
        .sum();
    ke / pop.len() as f64
}
/// Compute enstrophy (sum of squared vorticity) on a 2D grid from velocity fields.
///
/// Uses centered finite differences with periodic wrapping.
pub fn enstrophy_2d(pop: &[[f64; 9]], nx: usize, ny: usize) -> f64 {
    let mut enst = 0.0f64;
    for y in 0..ny {
        for x in 0..nx {
            let xp = (x + 1) % nx;
            let xm = (x + nx - 1) % nx;
            let yp = (y + 1) % ny;
            let ym = (y + ny - 1) % ny;
            let (_, _, uy_xp) = macros_from_d2q9(&pop[y * nx + xp]);
            let (_, _, uy_xm) = macros_from_d2q9(&pop[y * nx + xm]);
            let (_, ux_yp, _) = macros_from_d2q9(&pop[yp * nx + x]);
            let (_, ux_ym, _) = macros_from_d2q9(&pop[ym * nx + x]);
            let omega_z = 0.5 * (uy_xp - uy_xm) - 0.5 * (ux_yp - ux_ym);
            enst += omega_z * omega_z;
        }
    }
    enst
}
/// Count fluid cells (non-wall) in a `CellularGrid2D`.
pub fn count_fluid_cells(grid: &CellularGrid2D) -> usize {
    grid.wall.iter().filter(|&&w| !w).count()
}
/// Compute mean density over all fluid cells.
pub fn mean_density(grid: &CellularGrid2D) -> f64 {
    let fluid: Vec<f64> = grid
        .pop
        .iter()
        .enumerate()
        .filter(|(i, _)| !grid.wall[*i])
        .map(|(_, node)| macros_from_d2q9(node).0)
        .collect();
    if fluid.is_empty() {
        return 0.0;
    }
    fluid.iter().sum::<f64>() / fluid.len() as f64
}
/// Apply zero-gradient (Neumann) outflow at the east boundary `x = nx-1`.
///
/// Copies populations from `x = nx-2` to `x = nx-1`.
pub fn outflow_neumann_east(pop: &mut [[f64; 9]], nx: usize, ny: usize) {
    for y in 0..ny {
        pop[y * nx + (nx - 1)] = pop[y * nx + (nx - 2)];
    }
}
/// Apply constant-density (Dirichlet) pressure inlet at the west boundary.
pub fn pressure_inlet_west(pop: &mut [[f64; 9]], nx: usize, ny: usize, rho_in: f64) {
    for y in 0..ny {
        let idx = y * nx;
        let (_, ux, uy) = macros_from_d2q9(&pop[idx]);
        pop[idx] = equilibrium_d2q9(rho_in, ux, uy);
    }
}
/// Apply constant-density outlet at the east boundary.
pub fn pressure_outlet_east(pop: &mut [[f64; 9]], nx: usize, ny: usize, rho_out: f64) {
    for y in 0..ny {
        let idx = y * nx + (nx - 1);
        let (_, ux, uy) = macros_from_d2q9(&pop[idx]);
        pop[idx] = equilibrium_d2q9(rho_out, ux, uy);
    }
}
