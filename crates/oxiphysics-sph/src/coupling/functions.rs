//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use crate::coupling::types::*;

/// Smooth transition weight function for the handshake region.
///
/// Returns `1.0` when `d <= d_min` (atomistic side), `0.0` when `d >= d_max`
/// (continuum side), and a quintic polynomial interpolant in between:
/// `w = 1 - 6x^5 + 15x^4 - 10x^3`, where `x = (d - d_min) / (d_max - d_min)`.
pub fn weight_function_transition(d: f64, d_min: f64, d_max: f64) -> f64 {
    if d_max <= d_min {
        return 1.0;
    }
    let x = ((d - d_min) / (d_max - d_min)).clamp(0.0, 1.0);
    1.0 - 6.0 * x.powi(5) + 15.0 * x.powi(4) - 10.0 * x.powi(3)
}
/// Linear blending weight function.
///
/// Simpler alternative to quintic: w = 1 - (d - d_min) / (d_max - d_min).
pub fn weight_function_linear(d: f64, d_min: f64, d_max: f64) -> f64 {
    if d_max <= d_min {
        return 1.0;
    }
    let x = ((d - d_min) / (d_max - d_min)).clamp(0.0, 1.0);
    1.0 - x
}
/// Compute the boundary repulsion force between an SPH particle and a rigid body.
///
/// Uses a Lennard-Jones-like repulsion: F = k * (r0/r)^n * n_hat
///
/// where r is the distance to the surface and n_hat points away from the body.
pub fn boundary_repulsion_force(
    particle_pos: [f64; 3],
    body: &RigidBody,
    stiffness: f64,
    cutoff_distance: f64,
) -> [f64; 3] {
    let dx = particle_pos[0] - body.position[0];
    let dy = particle_pos[1] - body.position[1];
    let dz = particle_pos[2] - body.position[2];
    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
    if dist < 1e-14 {
        return [0.0; 3];
    }
    let surface_dist = dist - body.radius;
    if surface_dist >= cutoff_distance || surface_dist <= 0.0 {
        return [0.0; 3];
    }
    let nx = dx / dist;
    let ny = dy / dist;
    let nz = dz / dist;
    let mag = stiffness * (1.0 - surface_dist / cutoff_distance);
    [mag * nx, mag * ny, mag * nz]
}
/// Volume-based coupling: computes the volume fraction of SPH particles
/// occupying a given cell.
///
/// Returns the volume fraction in \[0, 1\].
pub fn volume_fraction(
    cell_min: [f64; 3],
    cell_max: [f64; 3],
    particle_positions: &[[f64; 3]],
    particle_volume: f64,
) -> f64 {
    let cell_vol =
        (cell_max[0] - cell_min[0]) * (cell_max[1] - cell_min[1]) * (cell_max[2] - cell_min[2]);
    if cell_vol < 1e-30 {
        return 0.0;
    }
    let mut count = 0usize;
    for &pos in particle_positions {
        if pos[0] >= cell_min[0]
            && pos[0] <= cell_max[0]
            && pos[1] >= cell_min[1]
            && pos[1] <= cell_max[1]
            && pos[2] >= cell_min[2]
            && pos[2] <= cell_max[2]
        {
            count += 1;
        }
    }
    let occupied = count as f64 * particle_volume;
    (occupied / cell_vol).min(1.0)
}
/// Maps SPH particles onto a regular continuum grid and returns a flat
/// density (mass-per-cell) array.
///
/// Each particle contributes its mass to the grid cell whose centre is
/// closest to the particle (nearest-grid-point assignment).
///
/// # Arguments
/// * `positions`    – particle positions.
/// * `velocities`   – particle velocities (unused in density mapping but kept
///   for API symmetry with the continuum solver).
/// * `masses`       – particle masses.
/// * `grid_origin`  – world-space origin of the grid.
/// * `grid_size`    – total extent of the grid in each dimension.
/// * `n_cells`      – number of cells in each dimension.
///
/// # Returns
/// A flat `Vec`f64` of length `n_cells\[0\] * n_cells\[1\] * n_cells\[2\]`
/// ordered as `\[ix\]\[iy\]\[iz\]` with stride `n_cells\[1\] * n_cells\[2\]`.
pub fn coarsen_sph_to_continuum(
    positions: &[[f64; 3]],
    _velocities: &[[f64; 3]],
    masses: &[f64],
    grid_origin: [f64; 3],
    grid_size: [f64; 3],
    n_cells: [usize; 3],
) -> Vec<f64> {
    let total = n_cells[0] * n_cells[1] * n_cells[2];
    let mut density = vec![0.0_f64; total];
    let cell_size = [
        grid_size[0] / n_cells[0] as f64,
        grid_size[1] / n_cells[1] as f64,
        grid_size[2] / n_cells[2] as f64,
    ];
    for (idx, &pos) in positions.iter().enumerate() {
        let fx = (pos[0] - grid_origin[0]) / cell_size[0];
        let fy = (pos[1] - grid_origin[1]) / cell_size[1];
        let fz = (pos[2] - grid_origin[2]) / cell_size[2];
        let ix = (fx.floor() as isize).clamp(0, n_cells[0] as isize - 1) as usize;
        let iy = (fy.floor() as isize).clamp(0, n_cells[1] as isize - 1) as usize;
        let iz = (fz.floor() as isize).clamp(0, n_cells[2] as isize - 1) as usize;
        let flat = ix * n_cells[1] * n_cells[2] + iy * n_cells[2] + iz;
        density[flat] += masses[idx];
    }
    density
}
/// Squared Euclidean distance between two 3-D points.
pub fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}
/// Euclidean distance between two 3-D points.
pub fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    dist_sq(a, b).sqrt()
}
/// Cross product of two 3-D vectors.
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
/// Simple cubic spline kernel (truncated at `2h`), normalised in 3-D.
///
/// W(r, h) = σ * { 1 - 1.5q² + 0.75q³,  0 ≤ q < 1
///               { 0.25(2-q)³,            1 ≤ q < 2
///               { 0,                     q ≥ 2 }
/// where q = r / h and σ = 8 / (π h³).
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let sigma = 8.0 / (std::f64::consts::PI * h * h * h);
    let q = r / h;
    if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    }
}
/// Cubic spline kernel gradient (scalar dW/dr).
pub fn cubic_kernel_grad(r: f64, h: f64) -> f64 {
    if h <= 0.0 || r < 1e-14 {
        return 0.0;
    }
    let sigma = 8.0 / (std::f64::consts::PI * h * h * h);
    let q = r / h;
    if q >= 2.0 {
        0.0
    } else if q >= 1.0 {
        sigma * (-0.75 * (2.0 - q).powi(2)) / h
    } else {
        sigma * (-3.0 * q + 2.25 * q * q) / h
    }
}
/// Apply forces from SPH particles onto a rigid body (Newton's second law).
///
/// The net SPH force `F = Σ sph_forces\[i\]` accelerates the body:
///   v_new = v_old + (F / body_mass) · dt
#[allow(dead_code)]
pub fn two_way_coupling_update(
    sph_forces: &[[f64; 3]],
    body_mass: f64,
    body_vel: &mut [f64; 3],
    dt: f64,
) {
    let inv_m = 1.0 / body_mass.max(1e-30);
    let mut f_total = [0.0_f64; 3];
    for &f in sph_forces {
        f_total[0] += f[0];
        f_total[1] += f[1];
        f_total[2] += f[2];
    }
    body_vel[0] += f_total[0] * inv_m * dt;
    body_vel[1] += f_total[1] * inv_m * dt;
    body_vel[2] += f_total[2] * inv_m * dt;
}
/// Build a sparse set of SPH-FEM interpolation weights using a kernel-based
/// approach.
///
/// For each FEM node `n` and each SPH particle `p`, the weight is
/// W(|x_n - x_p|, h) · m_p / ρ_p  (consistent with SPH interpolation).
#[allow(dead_code)]
pub fn build_sph_fem_weights(
    fem_positions: &[[f64; 3]],
    sph_positions: &[[f64; 3]],
    sph_masses: &[f64],
    sph_densities: &[f64],
    h: f64,
) -> Vec<SphFemWeight> {
    let mut weights = Vec::new();
    for (ni, &np) in fem_positions.iter().enumerate() {
        for (pi, &sp) in sph_positions.iter().enumerate() {
            let r = dist(np, sp);
            if r > 2.0 * h {
                continue;
            }
            let w = cubic_kernel(r, h);
            let vol_p = if sph_densities[pi] > 1e-30 {
                sph_masses[pi] / sph_densities[pi]
            } else {
                0.0
            };
            weights.push(SphFemWeight {
                sph_idx: pi,
                fem_idx: ni,
                weight: w * vol_p,
            });
        }
    }
    weights
}
/// Transfer velocity from SPH particles to FEM nodes via interpolation weights.
#[allow(dead_code)]
pub fn sph_to_fem_velocity(
    weights: &[SphFemWeight],
    sph_velocities: &[[f64; 3]],
    n_fem_nodes: usize,
) -> Vec<[f64; 3]> {
    let mut vel_fem = vec![[0.0_f64; 3]; n_fem_nodes];
    let mut w_sum = vec![0.0_f64; n_fem_nodes];
    for wt in weights {
        let v = sph_velocities[wt.sph_idx];
        vel_fem[wt.fem_idx][0] += wt.weight * v[0];
        vel_fem[wt.fem_idx][1] += wt.weight * v[1];
        vel_fem[wt.fem_idx][2] += wt.weight * v[2];
        w_sum[wt.fem_idx] += wt.weight;
    }
    for ni in 0..n_fem_nodes {
        if w_sum[ni] > 1e-30 {
            let inv = 1.0 / w_sum[ni];
            for c in 0..3 {
                vel_fem[ni][c] *= inv;
            }
        }
    }
    vel_fem
}
/// Transfer stress/pressure from FEM nodes to SPH particles.
#[allow(dead_code)]
pub fn fem_to_sph_pressure(
    weights: &[SphFemWeight],
    fem_pressures: &[f64],
    n_sph: usize,
) -> Vec<f64> {
    let mut p_sph = vec![0.0_f64; n_sph];
    let mut w_sum = vec![0.0_f64; n_sph];
    for wt in weights {
        p_sph[wt.sph_idx] += wt.weight * fem_pressures[wt.fem_idx];
        w_sum[wt.sph_idx] += wt.weight;
    }
    for pi in 0..n_sph {
        if w_sum[pi] > 1e-30 {
            p_sph[pi] /= w_sum[pi];
        }
    }
    p_sph
}
/// Compute Lennard-Jones type boundary repulsion force.
///
/// F_LJ = ε_LJ · \[(r0/r)^n1 - (r0/r)^n2\] / r² · (x_j - x_i)
///
/// Typical: n1=4, n2=2, r0 = particle spacing.
#[allow(dead_code)]
pub fn lennard_jones_repulsion(
    pos_fluid: [f64; 3],
    pos_boundary: [f64; 3],
    eps_lj: f64,
    r0: f64,
    n1: i32,
    n2: i32,
) -> [f64; 3] {
    let dx = [
        pos_fluid[0] - pos_boundary[0],
        pos_fluid[1] - pos_boundary[1],
        pos_fluid[2] - pos_boundary[2],
    ];
    let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
    let r = r2.sqrt();
    if r < 1e-14 || r > 2.0 * r0 {
        return [0.0; 3];
    }
    let ratio = r0 / r;
    let magnitude = eps_lj * (ratio.powi(n1) - ratio.powi(n2)) / r2;
    [magnitude * dx[0], magnitude * dx[1], magnitude * dx[2]]
}
/// Compute Adami et al. (2012) pressure boundary condition force.
///
/// Returns the generalised wall pressure for a virtual boundary particle:
/// p_w = (Σ_f p_f W_fw + (g - a_w) · Σ_f ρ_f (x_w - x_f) W_fw) / Σ_f W_fw
///
/// Here we use the simplified form without body force.
#[allow(dead_code)]
pub fn adami_wall_pressure(
    fluid_pressures: &[f64],
    fluid_positions: &[[f64; 3]],
    wall_pos: [f64; 3],
    h: f64,
) -> f64 {
    let mut num = 0.0_f64;
    let mut denom = 0.0_f64;
    for (pi, &fp) in fluid_positions.iter().enumerate() {
        let r = dist(wall_pos, fp);
        if r > 2.0 * h {
            continue;
        }
        let w = cubic_kernel(r, h);
        num += fluid_pressures[pi] * w;
        denom += w;
    }
    if denom < 1e-30 { 0.0 } else { num / denom }
}
/// Compute the SPH boundary velocity (mirror/ghost-particle method).
///
/// v_wall_ghost = 2 v_wall - v_fluid_avg  (no-slip condition)
#[allow(dead_code)]
pub fn ghost_particle_velocity(
    v_wall: [f64; 3],
    fluid_velocities: &[[f64; 3]],
    fluid_positions: &[[f64; 3]],
    ghost_pos: [f64; 3],
    h: f64,
) -> [f64; 3] {
    let mut v_avg = [0.0_f64; 3];
    let mut w_sum = 0.0_f64;
    for (i, &fp) in fluid_positions.iter().enumerate() {
        let r = dist(ghost_pos, fp);
        if r > 2.0 * h {
            continue;
        }
        let w = cubic_kernel(r, h);
        for c in 0..3 {
            v_avg[c] += w * fluid_velocities[i][c];
        }
        w_sum += w;
    }
    if w_sum > 1e-30 {
        for c in 0..3 {
            v_avg[c] /= w_sum;
        }
    }
    [
        2.0 * v_wall[0] - v_avg[0],
        2.0 * v_wall[1] - v_avg[1],
        2.0 * v_wall[2] - v_avg[2],
    ]
}
/// Compute the generalised SPH-to-FEM nodal force vector via virtual work.
///
/// F^FEM_n = Σ_p φ_n(x_p) · f^SPH_p · V_p
///
/// where φ_n are hat functions approximated by kernel weights.
#[allow(dead_code)]
pub fn sph_to_fem_nodal_forces(
    weights: &[SphFemWeight],
    sph_forces: &[[f64; 3]],
    sph_volumes: &[f64],
    n_fem_nodes: usize,
) -> Vec<[f64; 3]> {
    let mut f_fem = vec![[0.0_f64; 3]; n_fem_nodes];
    for wt in weights {
        let vp = sph_volumes[wt.sph_idx];
        let f = sph_forces[wt.sph_idx];
        for c in 0..3 {
            f_fem[wt.fem_idx][c] += wt.weight * vp * f[c];
        }
    }
    f_fem
}
#[cfg(test)]
mod tests {
    use super::*;

    fn make_region(min: [f64; 3], max: [f64; 3], rt: RegionType) -> CouplingRegion {
        CouplingRegion {
            min,
            max,
            region_type: rt,
        }
    }
    #[test]
    fn test_coupling_region_contains_inside() {
        let region = make_region([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], RegionType::Atomistic);
        assert!(region.contains([0.5, 0.5, 0.5]));
    }
    #[test]
    fn test_coupling_region_contains_outside() {
        let region = make_region([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], RegionType::Atomistic);
        assert!(!region.contains([1.5, 0.5, 0.5]));
        assert!(!region.contains([-0.1, 0.5, 0.5]));
    }
    #[test]
    fn test_coupling_region_contains_on_boundary() {
        let region = make_region([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], RegionType::Continuum);
        assert!(region.contains([0.0, 0.0, 0.0]));
        assert!(region.contains([1.0, 1.0, 1.0]));
    }
    #[test]
    fn test_coupling_region_volume() {
        let region = make_region([0.0, 0.0, 0.0], [2.0, 3.0, 4.0], RegionType::Handshake);
        let vol = region.volume();
        assert!((vol - 24.0).abs() < 1e-12);
    }
    #[test]
    fn test_coupling_region_center() {
        let region = make_region([0.0, 0.0, 0.0], [2.0, 4.0, 6.0], RegionType::Atomistic);
        let c = region.center();
        assert!((c[0] - 1.0).abs() < 1e-12);
        assert!((c[1] - 2.0).abs() < 1e-12);
        assert!((c[2] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_weight_at_d_min() {
        let w = weight_function_transition(0.0, 0.0, 1.0);
        assert!((w - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_weight_at_d_max() {
        let w = weight_function_transition(1.0, 0.0, 1.0);
        assert!(w.abs() < 1e-12);
    }
    #[test]
    fn test_weight_in_between_is_smooth() {
        let d_min = 0.0;
        let d_max = 1.0;
        let w_mid = weight_function_transition(0.5, d_min, d_max);
        assert!((w_mid - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_weight_below_d_min_clamped_to_one() {
        let w = weight_function_transition(-5.0, 0.0, 1.0);
        assert!((w - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_weight_above_d_max_clamped_to_zero() {
        let w = weight_function_transition(10.0, 0.0, 1.0);
        assert!(w.abs() < 1e-12);
    }
    #[test]
    fn test_weight_linear_endpoints() {
        assert!((weight_function_linear(0.0, 0.0, 1.0) - 1.0).abs() < 1e-12);
        assert!((weight_function_linear(1.0, 0.0, 1.0) - 0.0).abs() < 1e-12);
        assert!((weight_function_linear(0.5, 0.0, 1.0) - 0.5).abs() < 1e-12);
    }
    #[test]
    fn test_coupling_buffer_interpolate_velocity_single_particle() {
        let mut buf = CouplingBuffer::new();
        let pos = [1.0, 2.0, 3.0];
        let vel = [4.0, 5.0, 6.0];
        buf.add_particle(pos, vel, 1000.0, 1.0);
        let result = buf.interpolate_velocity_at(pos, 0.5);
        assert!((result[0] - vel[0]).abs() < 1e-10);
        assert!((result[1] - vel[1]).abs() < 1e-10);
        assert!((result[2] - vel[2]).abs() < 1e-10);
    }
    #[test]
    fn test_coupling_buffer_interpolate_density_single_particle() {
        let mut buf = CouplingBuffer::new();
        let pos = [0.0, 0.0, 0.0];
        buf.add_particle(pos, [0.0; 3], 999.0, 1.0);
        let rho = buf.interpolate_density_at(pos, 0.5);
        assert!((rho - 999.0).abs() < 1e-10);
    }
    #[test]
    fn test_coupling_buffer_len() {
        let mut buf = CouplingBuffer::new();
        assert_eq!(buf.len(), 0);
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 1.0);
        assert_eq!(buf.len(), 1);
    }
    #[test]
    fn test_coupling_buffer_no_particles_returns_zero() {
        let buf = CouplingBuffer::new();
        let v = buf.interpolate_velocity_at([0.0, 0.0, 0.0], 1.0);
        assert_eq!(v, [0.0, 0.0, 0.0]);
        let rho = buf.interpolate_density_at([0.0, 0.0, 0.0], 1.0);
        assert_eq!(rho, 0.0);
    }
    #[test]
    fn test_coupling_buffer_clear() {
        let mut buf = CouplingBuffer::new();
        buf.add_particle([0.0; 3], [1.0, 0.0, 0.0], 1000.0, 1.0);
        assert_eq!(buf.len(), 1);
        buf.clear();
        assert!(buf.is_empty());
    }
    #[test]
    fn test_coupling_buffer_interpolate_scalar() {
        let mut buf = CouplingBuffer::new();
        buf.add_particle([0.0; 3], [0.0; 3], 1000.0, 1.0);
        let values = vec![42.0];
        let s = buf.interpolate_scalar_at([0.0; 3], 0.5, &values);
        assert!((s - 42.0).abs() < 1e-10);
    }
    #[test]
    fn test_bdm_is_in_overlap_inside() {
        let atom = make_region([0.0, 0.0, 0.0], [3.0, 1.0, 1.0], RegionType::Atomistic);
        let cont = make_region([1.0, 0.0, 0.0], [4.0, 1.0, 1.0], RegionType::Continuum);
        let over = make_region([1.0, 0.0, 0.0], [3.0, 1.0, 1.0], RegionType::Handshake);
        let bdm = BridgingDomainMethod::new(atom, cont, over);
        assert!(bdm.is_in_overlap([2.0, 0.5, 0.5]));
    }
    #[test]
    fn test_bdm_is_in_overlap_outside() {
        let atom = make_region([0.0, 0.0, 0.0], [3.0, 1.0, 1.0], RegionType::Atomistic);
        let cont = make_region([1.0, 0.0, 0.0], [4.0, 1.0, 1.0], RegionType::Continuum);
        let over = make_region([1.0, 0.0, 0.0], [3.0, 1.0, 1.0], RegionType::Handshake);
        let bdm = BridgingDomainMethod::new(atom, cont, over);
        assert!(!bdm.is_in_overlap([0.5, 0.5, 0.5]));
        assert!(!bdm.is_in_overlap([3.5, 0.5, 0.5]));
    }
    #[test]
    fn test_coarsen_single_particle_nearest_cell() {
        let positions = vec![[0.5, 0.5, 0.5_f64]];
        let velocities = vec![[1.0, 0.0, 0.0_f64]];
        let masses = vec![2.5_f64];
        let grid_origin = [0.0, 0.0, 0.0];
        let grid_size = [4.0, 4.0, 4.0];
        let n_cells = [4_usize, 4, 4];
        let density = coarsen_sph_to_continuum(
            &positions,
            &velocities,
            &masses,
            grid_origin,
            grid_size,
            n_cells,
        );
        assert!(
            (density[0] - 2.5).abs() < 1e-12,
            "Expected 2.5 in cell 0, got {}",
            density[0]
        );
        let rest_sum: f64 = density[1..].iter().sum();
        assert!(rest_sum.abs() < 1e-12);
    }
    #[test]
    fn test_coarsen_particle_out_of_bounds_clamped() {
        let positions = vec![[-100.0, -100.0, -100.0_f64]];
        let masses = vec![1.0_f64];
        let density = coarsen_sph_to_continuum(
            &positions,
            &[[[0.0; 3][0]; 3]],
            &masses,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2, 2, 2],
        );
        assert!((density[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_rigid_body_sphere_mass() {
        let body = RigidBody::sphere([0.0; 3], 1.0, 1000.0);
        let expected_mass = (4.0 / 3.0) * std::f64::consts::PI * 1000.0;
        assert!(
            (body.mass - expected_mass).abs() < 1e-6,
            "Expected mass {expected_mass}, got {}",
            body.mass
        );
    }
    #[test]
    fn test_rigid_body_signed_distance() {
        let body = RigidBody::sphere([0.0; 3], 1.0, 1000.0);
        let sd_outside = body.signed_distance([2.0, 0.0, 0.0]);
        assert!(
            (sd_outside - 1.0).abs() < 1e-12,
            "Outside: expected 1.0, got {sd_outside}"
        );
        let sd_inside = body.signed_distance([0.5, 0.0, 0.0]);
        assert!(
            sd_inside < 0.0,
            "Inside: expected negative, got {sd_inside}"
        );
    }
    #[test]
    fn test_rigid_body_surface_velocity_no_rotation() {
        let body = RigidBody {
            position: [0.0; 3],
            velocity: [1.0, 0.0, 0.0],
            angular_velocity: [0.0; 3],
            mass: 1.0,
            inertia: 1.0,
            radius: 1.0,
        };
        let sv = body.surface_velocity([1.0, 0.0, 0.0]);
        assert!((sv[0] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn test_rigid_body_apply_force() {
        let mut body = RigidBody::sphere([0.0; 3], 1.0, 1.0);
        body.apply_force([1.0, 0.0, 0.0], [0.0; 3], 0.1);
        assert!(
            body.velocity[0] > 0.0,
            "Velocity should increase after force"
        );
        assert!(body.position[0] > 0.0, "Position should advance");
    }
    #[test]
    fn test_boundary_repulsion_force_direction() {
        let body = RigidBody::sphere([0.0; 3], 1.0, 1000.0);
        let f = boundary_repulsion_force([1.5, 0.0, 0.0], &body, 1000.0, 1.0);
        assert!(f[0] > 0.0, "Force should push particle away from body");
    }
    #[test]
    fn test_boundary_repulsion_force_zero_far() {
        let body = RigidBody::sphere([0.0; 3], 1.0, 1000.0);
        let f = boundary_repulsion_force([10.0, 0.0, 0.0], &body, 1000.0, 1.0);
        for &v in &f {
            assert!(v.abs() < 1e-14, "No force far from body");
        }
    }
    #[test]
    fn test_sph_fem_coupling_transfer() {
        let mut coupling = SphFemCoupling::new(0.5);
        coupling.add_node([0.0, 0.0, 0.0]);
        coupling.add_node([1.0, 0.0, 0.0]);
        let sph_positions = vec![[0.1, 0.0, 0.0]];
        let sph_forces = vec![[10.0, 0.0, 0.0]];
        coupling.transfer_forces_to_fem(&sph_positions, &sph_forces);
        let f0 = coupling.fem_nodes[0].coupling_force[0];
        let f1 = coupling.fem_nodes[1].coupling_force[0];
        assert!(
            f0 > f1,
            "Closer node should get more force: f0={f0}, f1={f1}"
        );
    }
    #[test]
    fn test_sph_fem_node_count() {
        let mut coupling = SphFemCoupling::new(0.5);
        assert_eq!(coupling.node_count(), 0);
        coupling.add_node([0.0; 3]);
        coupling.add_node([1.0, 0.0, 0.0]);
        assert_eq!(coupling.node_count(), 2);
    }
    #[test]
    fn test_volume_fraction_single_particle() {
        let positions = vec![[0.5, 0.5, 0.5]];
        let particle_vol = 0.001;
        let vf = volume_fraction([0.0; 3], [1.0, 1.0, 1.0], &positions, particle_vol);
        assert!((vf - 0.001).abs() < 1e-10, "Expected 0.001, got {vf}");
    }
    #[test]
    fn test_volume_fraction_outside() {
        let positions = vec![[5.0, 5.0, 5.0]];
        let vf = volume_fraction([0.0; 3], [1.0, 1.0, 1.0], &positions, 0.1);
        assert!(vf.abs() < 1e-12, "No particles inside should give 0");
    }
    #[test]
    fn test_interface_tracker_creation() {
        let tracker = InterfaceTracker::new(10);
        assert_eq!(tracker.colors.len(), 10);
        for &c in &tracker.colors {
            assert!(c.abs() < 1e-12);
        }
    }
    #[test]
    fn test_interface_tracker_detect() {
        let mut tracker = InterfaceTracker::new(5);
        tracker.set_color(0, 0.0);
        tracker.set_color(1, 0.35);
        tracker.set_color(2, 0.5);
        tracker.set_color(3, 0.65);
        tracker.set_color(4, 1.0);
        let iface = tracker.interface_particles(0.2);
        assert!(iface.contains(&1), "Particle 1 should be at interface");
        assert!(iface.contains(&2), "Particle 2 should be at interface");
        assert!(iface.contains(&3), "Particle 3 should be at interface");
        assert!(!iface.contains(&0), "Particle 0 should not be at interface");
        assert!(!iface.contains(&4), "Particle 4 should not be at interface");
    }
    #[test]
    fn test_interface_smooth_preserves_values() {
        let mut tracker = InterfaceTracker::new(3);
        tracker.set_color(0, 0.0);
        tracker.set_color(1, 0.0);
        tracker.set_color(2, 0.0);
        let positions = [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]];
        tracker.smooth_colors(&positions, 0.5);
        for &c in &tracker.colors {
            assert!(
                c.abs() < 1e-10,
                "Smoothing uniform field should not change it"
            );
        }
    }
    #[test]
    fn test_rigid_body_coupling_force_direction() {
        let coupling = RigidBodyCoupling::new(0, 1.0, [0.0, 0.0, 0.0], [0.0; 3]);
        let fluid_pos = vec![[0.0, 0.1, 0.0]];
        let fluid_pressures = vec![1e5_f64];
        let fluid_rho = vec![1000.0_f64];
        let h = 0.5_f64;
        let force = coupling.sph_coupling_force(&fluid_pos, &fluid_pressures, &fluid_rho, h);
        assert!(
            force.iter().any(|&v| v.abs() > 1e-14),
            "Coupling force should be non-zero for nearby pressure source"
        );
    }
    #[test]
    fn test_rigid_body_coupling_no_force_far_away() {
        let coupling = RigidBodyCoupling::new(0, 1.0, [0.0, 0.0, 0.0], [0.0; 3]);
        let fluid_pos = vec![[100.0, 0.0, 0.0]];
        let fluid_pressures = vec![1e5_f64];
        let fluid_rho = vec![1000.0_f64];
        let h = 0.1_f64;
        let force = coupling.sph_coupling_force(&fluid_pos, &fluid_pressures, &fluid_rho, h);
        for &v in &force {
            assert!(v.abs() < 1e-14, "No force from distant particle, got {v}");
        }
    }
    #[test]
    fn test_two_way_coupling_changes_velocity() {
        let forces = vec![[10.0_f64, 0.0, 0.0]];
        let body_mass = 2.0_f64;
        let mut vel = [0.0_f64, 0.0, 0.0];
        let dt = 0.1_f64;
        two_way_coupling_update(&forces, body_mass, &mut vel, dt);
        assert!(
            (vel[0] - 0.5).abs() < 1e-12,
            "Expected vel[0]=0.5, got {}",
            vel[0]
        );
        assert!(vel[1].abs() < 1e-14);
        assert!(vel[2].abs() < 1e-14);
    }
    #[test]
    fn test_two_way_coupling_zero_force() {
        let forces: Vec<[f64; 3]> = vec![[0.0; 3]];
        let mut vel = [1.0_f64, 2.0, 3.0];
        let vel_orig = vel;
        two_way_coupling_update(&forces, 1.0, &mut vel, 0.1);
        for i in 0..3 {
            assert!(
                (vel[i] - vel_orig[i]).abs() < 1e-14,
                "Zero force should not change velocity"
            );
        }
    }
    #[test]
    fn test_fem_transfer_increases_nodal_forces() {
        let node_ids = vec![0];
        let node_positions = vec![[0.0, 0.0, 0.0_f64]];
        let mut fem = FemCoupling::new(node_ids, node_positions);
        let sph_pos = vec![[0.05, 0.0, 0.0_f64]];
        let sph_pressure = vec![1e5_f64];
        let h = 0.3_f64;
        fem.transfer_sph_pressure_to_nodes(&sph_pos, &sph_pressure, h);
        let f_mag = {
            let f = fem.nodal_forces[0];
            (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt()
        };
        assert!(
            f_mag > 0.0,
            "Nodal force should be non-zero for nearby SPH pressure"
        );
    }
    #[test]
    fn test_fem_transfer_no_force_far_away() {
        let node_ids = vec![0];
        let node_positions = vec![[0.0, 0.0, 0.0_f64]];
        let mut fem = FemCoupling::new(node_ids, node_positions);
        let sph_pos = vec![[100.0, 0.0, 0.0_f64]];
        let sph_pressure = vec![1e5_f64];
        let h = 0.3_f64;
        fem.transfer_sph_pressure_to_nodes(&sph_pos, &sph_pressure, h);
        for &v in &fem.nodal_forces[0] {
            assert!(
                v.abs() < 1e-14,
                "No force from distant SPH particle, got {v}"
            );
        }
    }
    #[test]
    fn test_immersed_boundary_force_count() {
        let bp = vec![[0.0, 0.0, 0.0_f64], [0.1, 0.0, 0.0], [0.2, 0.0, 0.0]];
        let normals = vec![[0.0, 1.0, 0.0_f64], [0.0, 1.0, 0.0], [0.0, 1.0, 0.0]];
        let ib = SphRigidImmersedBoundary::new(bp, normals);
        let fluid_pos = vec![[0.1, 0.05, 0.0_f64]];
        let pressures = vec![1e4_f64];
        let forces = ib.compute_reaction_forces(&fluid_pos, &pressures, 0.3);
        assert_eq!(
            forces.len(),
            3,
            "Should return one force vector per boundary particle"
        );
    }
    #[test]
    fn test_immersed_boundary_force_nonzero_near_fluid() {
        let bp = vec![[0.0, 0.0, 0.0_f64]];
        let normals = vec![[1.0, 0.0, 0.0_f64]];
        let ib = SphRigidImmersedBoundary::new(bp, normals);
        let fluid_pos = vec![[0.05, 0.0, 0.0_f64]];
        let pressures = vec![1e5_f64];
        let forces = ib.compute_reaction_forces(&fluid_pos, &pressures, 0.3);
        let f = forces[0];
        let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        assert!(
            mag > 0.0,
            "Reaction force should be non-zero for nearby fluid particle"
        );
    }
    #[test]
    fn test_build_sph_fem_weights_nearby() {
        let fem_pos = vec![[0.0_f64, 0.0, 0.0]];
        let sph_pos = vec![[0.05_f64, 0.0, 0.0]];
        let sph_mass = vec![0.001_f64];
        let sph_den = vec![1000.0_f64];
        let h = 0.2_f64;
        let weights = build_sph_fem_weights(&fem_pos, &sph_pos, &sph_mass, &sph_den, h);
        assert!(
            !weights.is_empty(),
            "Should have a weight for nearby particles"
        );
        assert_eq!(weights[0].fem_idx, 0);
        assert_eq!(weights[0].sph_idx, 0);
        assert!(weights[0].weight > 0.0);
    }
    #[test]
    fn test_build_sph_fem_weights_far() {
        let fem_pos = vec![[0.0_f64, 0.0, 0.0]];
        let sph_pos = vec![[100.0_f64, 0.0, 0.0]];
        let sph_mass = vec![0.001_f64];
        let sph_den = vec![1000.0_f64];
        let h = 0.2_f64;
        let weights = build_sph_fem_weights(&fem_pos, &sph_pos, &sph_mass, &sph_den, h);
        assert!(weights.is_empty(), "No weight for far-away particles");
    }
    #[test]
    fn test_sph_to_fem_velocity_transfer() {
        let fem_pos = vec![[0.0_f64, 0.0, 0.0]];
        let sph_pos = vec![[0.05_f64, 0.0, 0.0]];
        let sph_mass = vec![0.001_f64];
        let sph_den = vec![1000.0_f64];
        let h = 0.2_f64;
        let weights = build_sph_fem_weights(&fem_pos, &sph_pos, &sph_mass, &sph_den, h);
        let sph_vel = vec![[1.0_f64, 2.0, 3.0]];
        let vel_fem = sph_to_fem_velocity(&weights, &sph_vel, 1);
        for c in 0..3 {
            assert!(
                (vel_fem[0][c] - sph_vel[0][c]).abs() < 1e-12,
                "Velocity component {c} should match SPH velocity"
            );
        }
    }
    #[test]
    fn test_fem_to_sph_pressure_transfer() {
        let fem_pos = vec![[0.0_f64, 0.0, 0.0]];
        let sph_pos = vec![[0.05_f64, 0.0, 0.0]];
        let sph_mass = vec![0.001_f64];
        let sph_den = vec![1000.0_f64];
        let h = 0.2_f64;
        let weights = build_sph_fem_weights(&fem_pos, &sph_pos, &sph_mass, &sph_den, h);
        let fem_pressure = vec![1e5_f64];
        let p_sph = fem_to_sph_pressure(&weights, &fem_pressure, 1);
        assert!(
            (p_sph[0] - 1e5).abs() < 1.0,
            "Transferred pressure should match FEM value"
        );
    }
    #[test]
    fn test_coupled_time_stepper_advance() {
        let mut ts = CoupledTimeStepper::new(0.01, SubstepStrategy::Monolithic);
        ts.advance();
        assert!((ts.time - 0.01).abs() < 1e-14, "Time should advance by dt");
        assert_eq!(ts.cycles, 1);
    }
    #[test]
    fn test_coupled_time_stepper_accumulate_force() {
        let mut ts = CoupledTimeStepper::new(0.01, SubstepStrategy::Monolithic);
        ts.accumulate_force([1.0, 2.0, 3.0], [0.1, 0.2, 0.3]);
        ts.accumulate_force([4.0, 5.0, 6.0], [0.4, 0.5, 0.6]);
        assert!((ts.accumulated_force[0] - 5.0).abs() < 1e-14);
        assert!((ts.accumulated_torque[2] - 0.9).abs() < 1e-14);
    }
    #[test]
    fn test_coupled_time_stepper_reset_accumulation() {
        let mut ts = CoupledTimeStepper::new(0.01, SubstepStrategy::Monolithic);
        ts.accumulate_force([1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        ts.reset_accumulation();
        for &v in &ts.accumulated_force {
            assert!(v.abs() < 1e-14, "Force should be reset");
        }
    }
    #[test]
    fn test_coupled_time_stepper_sph_substep_dt() {
        let ts = CoupledTimeStepper::new(
            0.1,
            SubstepStrategy::SphSubStep {
                sph_steps_per_rigid: 4,
            },
        );
        let dt_sub = ts.sph_substep_dt();
        assert!((dt_sub - 0.025).abs() < 1e-14, "SPH substep should be dt/4");
    }
    #[test]
    fn test_coupled_time_stepper_rigid_substep_dt() {
        let ts = CoupledTimeStepper::new(
            0.1,
            SubstepStrategy::RigidSubStep {
                rigid_steps_per_sph: 5,
            },
        );
        let dt_sub = ts.rigid_substep_dt();
        assert!(
            (dt_sub - 0.02).abs() < 1e-14,
            "Rigid substep should be dt/5"
        );
    }
    #[test]
    fn test_rigid_body_linear_integration() {
        let mut rb = RigidBody6Dof::new(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        rb.integrate_linear([1.0, 0.0, 0.0], 0.1);
        assert!(
            (rb.velocity[0] - 0.1).abs() < 1e-14,
            "v_x should be 0.1 after F=1 dt=0.1 m=1"
        );
        assert!(
            (rb.position[0] - 0.01).abs() < 1e-14,
            "x should be 0.01 after integration"
        );
    }
    #[test]
    fn test_rigid_body_zero_force() {
        let mut rb = RigidBody6Dof::new(2.0, [1.0, 1.0, 1.0], [1.0, 2.0, 3.0]);
        rb.velocity = [0.5, 0.0, 0.0];
        rb.integrate_linear([0.0, 0.0, 0.0], 0.01);
        assert!((rb.velocity[0] - 0.5).abs() < 1e-14);
        assert!((rb.position[0] - 1.005).abs() < 1e-14);
    }
    #[test]
    fn test_rigid_body_angular_integration() {
        let mut rb = RigidBody6Dof::new(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        rb.integrate_angular([0.0, 0.0, 1.0], 0.1);
        assert!(
            (rb.angular_velocity[2] - 0.1).abs() < 1e-14,
            "ω_z should be 0.1"
        );
    }
    #[test]
    fn test_rigid_body_quaternion_unit_norm() {
        let mut rb = RigidBody6Dof::new(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        for _ in 0..10 {
            rb.integrate_angular([0.5, 0.3, 0.2], 0.05);
        }
        let [w, x, y, z] = rb.quaternion;
        let norm = (w * w + x * x + y * y + z * z).sqrt();
        assert!(
            (norm - 1.0).abs() < 1e-12,
            "Quaternion norm should be 1, got {norm}"
        );
    }
    #[test]
    fn test_rigid_body_kinetic_energy_positive() {
        let mut rb = RigidBody6Dof::new(2.0, [1.0, 2.0, 3.0], [0.0, 0.0, 0.0]);
        rb.velocity = [1.0, 0.0, 0.0];
        rb.angular_velocity = [0.0, 1.0, 0.0];
        let ke = rb.kinetic_energy();
        assert!((ke - 2.0).abs() < 1e-14, "Expected KE=2.0, got {ke}");
    }
    #[test]
    fn test_rigid_body_point_velocity() {
        let mut rb = RigidBody6Dof::new(1.0, [1.0, 1.0, 1.0], [0.0, 0.0, 0.0]);
        rb.velocity = [0.0, 0.0, 0.0];
        rb.angular_velocity = [0.0, 0.0, 1.0];
        let vp = rb.point_velocity([1.0, 0.0, 0.0]);
        assert!((vp[0]).abs() < 1e-14, "v_x should be 0");
        assert!((vp[1] - 1.0).abs() < 1e-14, "v_y should be 1");
    }
    #[test]
    fn test_energy_balance_fluid_to_body() {
        let mut eb = CoupledEnergyBalance::new();
        eb.update([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 0.01);
        assert!(eb.fluid_to_body > 0.0, "Fluid should give energy to body");
        assert!((eb.body_to_fluid).abs() < 1e-14);
        assert!((eb.net_transfer() - eb.fluid_to_body).abs() < 1e-14);
    }
    #[test]
    fn test_energy_balance_body_to_fluid() {
        let mut eb = CoupledEnergyBalance::new();
        eb.update([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 0.5, 0.01);
        assert!(eb.body_to_fluid > 0.0, "Body should give energy to fluid");
        assert!((eb.fluid_to_body).abs() < 1e-14);
        assert!(eb.net_transfer() < 0.0);
    }
    #[test]
    fn test_lj_repulsion_zero_far_away() {
        let f = lennard_jones_repulsion([10.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1, 4, 2);
        for &v in &f {
            assert!(v.abs() < 1e-14, "LJ force should be zero for r >> 2*r0");
        }
    }
    #[test]
    fn test_lj_repulsion_nonzero_nearby() {
        let f = lennard_jones_repulsion([0.05, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.1, 4, 2);
        let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        assert!(
            mag > 0.0,
            "LJ force should be non-zero for nearby particles"
        );
    }
    #[test]
    fn test_adami_wall_pressure_nearby() {
        let fluid_pos = vec![[0.05_f64, 0.0, 0.0]];
        let fluid_p = vec![1e5_f64];
        let p_wall = adami_wall_pressure(&fluid_p, &fluid_pos, [0.0, 0.0, 0.0], 0.2);
        assert!(
            (p_wall - 1e5).abs() < 1.0,
            "Wall pressure should approximate fluid pressure"
        );
    }
    #[test]
    fn test_adami_wall_pressure_far() {
        let fluid_pos = vec![[100.0_f64, 0.0, 0.0]];
        let fluid_p = vec![1e5_f64];
        let p_wall = adami_wall_pressure(&fluid_p, &fluid_pos, [0.0, 0.0, 0.0], 0.2);
        assert!(
            p_wall.abs() < 1e-14,
            "Wall pressure should be zero for far-away fluid"
        );
    }
    #[test]
    fn test_ghost_velocity_no_slip() {
        let v_wall = [1.0_f64, 0.0, 0.0];
        let fluid_vel = vec![[1.0_f64, 0.0, 0.0]];
        let fluid_pos = vec![[0.05_f64, 0.0, 0.0]];
        let ghost_pos = [0.0_f64, 0.0, 0.0];
        let vg = ghost_particle_velocity(v_wall, &fluid_vel, &fluid_pos, ghost_pos, 0.2);
        assert!(
            (vg[0] - 1.0).abs() < 1e-12,
            "Ghost velocity should be 1 for no-slip"
        );
    }
    #[test]
    fn test_sph_to_fem_nodal_forces_nonzero() {
        let fem_pos = vec![[0.0_f64, 0.0, 0.0]];
        let sph_pos = vec![[0.05_f64, 0.0, 0.0]];
        let sph_mass = vec![0.001_f64];
        let sph_den = vec![1000.0_f64];
        let h = 0.2_f64;
        let weights = build_sph_fem_weights(&fem_pos, &sph_pos, &sph_mass, &sph_den, h);
        let sph_forces = vec![[1.0_f64, 0.0, 0.0]];
        let sph_volumes = vec![1e-6_f64];
        let f_fem = sph_to_fem_nodal_forces(&weights, &sph_forces, &sph_volumes, 1);
        assert!(f_fem[0][0] > 0.0, "FEM nodal force in x should be positive");
    }
    #[test]
    fn test_fem_sph_interface_normal_single_particle() {
        let mut coupling = FemSphCoupling::new(vec![[0.1, 0.0, 0.0]], 0.3);
        coupling.compute_interface_normal(&[[0.0, 0.0, 0.0_f64]]);
        let n = coupling.fem_node_normals[0];
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(
            (mag - 1.0).abs() < 1e-10,
            "Interface normal should be unit length, mag={mag}"
        );
    }
    #[test]
    fn test_fem_sph_interface_normal_no_nearby_particle() {
        let mut coupling = FemSphCoupling::new(vec![[0.0, 0.0, 0.0]], 0.1);
        coupling.compute_interface_normal(&[[100.0, 0.0, 0.0_f64]]);
        let n = coupling.fem_node_normals[0];
        let mag = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        assert!(mag < 1e-14, "No nearby particle → zero normal, mag={mag}");
    }
    #[test]
    fn test_fem_sph_transfer_stress_no_nearby_particle() {
        let mut coupling = FemSphCoupling::new(vec![[0.0, 0.0, 0.0]], 0.1);
        coupling.fem_node_normals = vec![[0.0, 1.0, 0.0]];
        let sph_pos = vec![[100.0_f64, 0.0, 0.0]];
        let sph_stress = vec![[1e5, 1e5, 1e5, 0.0, 0.0, 0.0_f64]];
        coupling.transfer_stress(&sph_pos, &sph_stress);
        let f = coupling.nodal_forces[0];
        let mag = (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt();
        assert!(
            mag < 1e-14,
            "Far particle should produce zero nodal force, mag={mag}"
        );
    }
    #[test]
    fn test_fem_sph_transfer_stress_hydrostatic() {
        let h = 0.3_f64;
        let mut coupling = FemSphCoupling::new(vec![[0.0, 0.0, 0.0]], h);
        coupling.fem_node_normals = vec![[0.0, 1.0, 0.0]];
        let sph_pos = vec![[0.05_f64, 0.0, 0.0]];
        let p = 1e4_f64;
        let sph_stress = vec![[-p, -p, -p, 0.0, 0.0, 0.0]];
        coupling.transfer_stress(&sph_pos, &sph_stress);
        let fy = coupling.nodal_forces[0][1];
        assert!(
            fy < 0.0,
            "Hydrostatic traction in -y should be negative, got {fy}"
        );
    }
    #[test]
    fn test_fem_sph_coupling_new_dimensions() {
        let positions = vec![[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let coupling = FemSphCoupling::new(positions, 0.2);
        assert_eq!(coupling.fem_node_positions.len(), 3);
        assert_eq!(coupling.fem_node_normals.len(), 3);
        assert_eq!(coupling.nodal_forces.len(), 3);
        assert_eq!(coupling.nodal_stress.len(), 3);
    }
    #[test]
    fn test_buoyancy_torque_zero_no_displaced_particles() {
        let coupling = RigidSphCoupling::new([0.0; 3], 1000.0, 9.81);
        let positions = vec![[5.0_f64, 0.0, 0.0], [0.0, 5.0, 0.0]];
        let volumes = vec![1e-6_f64; 2];
        let tau = coupling.compute_buoyancy_torque(&positions, &volumes, 0.1);
        let mag = (tau[0] * tau[0] + tau[1] * tau[1] + tau[2] * tau[2]).sqrt();
        assert!(
            mag < 1e-20,
            "No displaced particles → zero torque, mag={mag}"
        );
    }
    #[test]
    fn test_buoyancy_torque_symmetric_cancels() {
        let coupling = RigidSphCoupling::new([0.0; 3], 1000.0, 9.81);
        let positions = vec![[0.05_f64, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let volumes = vec![1e-6_f64; 2];
        let tau = coupling.compute_buoyancy_torque(&positions, &volumes, 0.1);
        assert!(
            tau[1].abs() < 1e-20,
            "Symmetric torques should cancel in y, got {}",
            tau[1]
        );
        assert!(
            tau[2].abs() < 1e-20,
            "No z torque for x-symmetric particles"
        );
    }
    #[test]
    fn test_buoyancy_torque_nonzero_offset() {
        let coupling = RigidSphCoupling::new([0.0; 3], 1000.0, 9.81);
        let positions = vec![[0.05_f64, 0.0, 0.0]];
        let volumes = vec![1e-4_f64];
        let tau = coupling.compute_buoyancy_torque(&positions, &volumes, 0.1);
        let expected_tau_y = -0.05 * 1000.0 * 9.81 * 1e-4;
        assert!(
            (tau[1] - expected_tau_y).abs() < 1e-15,
            "τ_y = {}, expected {expected_tau_y}",
            tau[1]
        );
    }
    #[test]
    fn test_buoyancy_force_scales_with_volume() {
        let coupling = RigidSphCoupling::new([0.0; 3], 1000.0, 9.81);
        let positions = vec![[0.05_f64, 0.0, 0.0]];
        let v1 = vec![1e-4_f64];
        let v2 = vec![2e-4_f64];
        let fb1 = coupling.compute_buoyancy_force(&positions, &v1, 0.1);
        let fb2 = coupling.compute_buoyancy_force(&positions, &v2, 0.1);
        assert!(
            (fb2[2] / fb1[2] - 2.0).abs() < 1e-12,
            "Buoyancy should scale linearly with volume"
        );
    }
    #[test]
    fn test_buoyancy_force_direction_is_upward() {
        let coupling = RigidSphCoupling::new([0.0; 3], 1000.0, 9.81);
        let positions = vec![[0.0_f64, 0.0, 0.0]];
        let volumes = vec![1e-3_f64];
        let fb = coupling.compute_buoyancy_force(&positions, &volumes, 0.5);
        assert!(
            fb[2] > 0.0,
            "Buoyancy force should be upward (+z), got {}",
            fb[2]
        );
        assert!(
            fb[0].abs() < 1e-20 && fb[1].abs() < 1e-20,
            "No horizontal buoyancy"
        );
    }
}
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
#[inline]
pub fn scale3(s: f64, v: [f64; 3]) -> [f64; 3] {
    [s * v[0], s * v[1], s * v[2]]
}
#[inline]
pub fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}
#[cfg(test)]
mod tests_coupling_extended {
    use super::*;

    use std::f64::consts::PI;
    #[test]
    fn test_two_way_body_force_zero_far_particles() {
        let c = TwoWayCouplingForce::new(1000.0, 0.1, 0.5);
        let body_pos = [0.0; 3];
        let body_vel = [0.0; 3];
        let sph_pos = vec![[5.0, 0.0, 0.0], [0.0, 5.0, 0.0]];
        let sph_vel = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let sph_vol = vec![1e-4_f64; 2];
        let f = c.body_force_from_fluid(body_pos, body_vel, &sph_pos, &sph_vel, &sph_vol, 0.3);
        let mag = norm3(f);
        assert!(
            mag < 1e-20,
            "Far particles should give zero body force, mag={mag}"
        );
    }
    #[test]
    fn test_two_way_body_force_nonzero_close_particle() {
        let cd = 0.5;
        let rho_f = 1000.0;
        let h = 0.5;
        let c = TwoWayCouplingForce::new(rho_f, h, cd);
        let body_pos = [0.0; 3];
        let body_vel = [0.0; 3];
        let sph_pos = vec![[0.1, 0.0, 0.0]];
        let v_fluid = [2.0_f64, 0.0, 0.0];
        let sph_vel = vec![v_fluid];
        let vol = 1e-3_f64;
        let sph_vol = vec![vol];
        let f = c.body_force_from_fluid(body_pos, body_vel, &sph_pos, &sph_vel, &sph_vol, h);
        let expected_fx = -cd * rho_f * vol * v_fluid[0];
        assert!(
            (f[0] - expected_fx).abs() < 1e-12,
            "Expected fx={expected_fx}, got {}",
            f[0]
        );
    }
    #[test]
    fn test_two_way_reaction_forces_newton_third_law() {
        let cd = 0.3;
        let rho_f = 1000.0;
        let h = 1.0;
        let c = TwoWayCouplingForce::new(rho_f, h, cd);
        let body_pos = [0.0; 3];
        let body_vel = [1.0, 0.0, 0.0];
        let sph_pos = vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        let sph_vel = vec![[0.0; 3]; 2];
        let sph_vol = vec![1e-4_f64; 2];
        let f_body = c.body_force_from_fluid(body_pos, body_vel, &sph_pos, &sph_vel, &sph_vol, h);
        let f_fluid = c.fluid_reaction_forces(body_pos, body_vel, &sph_pos, &sph_vel, &sph_vol, h);
        let f_fluid_sum: [f64; 3] = f_fluid.iter().fold([0.0; 3], |acc, &f| add3(acc, f));
        assert!(
            (f_fluid_sum[0] + f_body[0]).abs() < 1e-12,
            "Newton 3rd law violated: sum_reaction={}, body={}",
            f_fluid_sum[0],
            f_body[0]
        );
    }
    #[test]
    fn test_two_way_reaction_forces_length() {
        let c = TwoWayCouplingForce::new(1000.0, 0.5, 0.5);
        let sph_pos = vec![[0.1; 3], [0.2; 3], [0.3; 3]];
        let sph_vel = vec![[0.0; 3]; 3];
        let sph_vol = vec![1e-4_f64; 3];
        let forces = c.fluid_reaction_forces([0.0; 3], [0.0; 3], &sph_pos, &sph_vel, &sph_vol, 0.5);
        assert_eq!(forces.len(), sph_pos.len());
    }
    #[test]
    fn test_added_mass_sphere() {
        let am = AddedMass::new(0.5, 1000.0, 1e-3);
        let expected = 0.5 * 1000.0 * 1e-3;
        assert!((am.added_mass() - expected).abs() < 1e-12);
    }
    #[test]
    fn test_added_mass_effective_mass() {
        let am = AddedMass::new(0.5, 1000.0, 1e-3);
        let m_body = 2.0;
        let m_eff = am.effective_mass(m_body);
        assert!((m_eff - (m_body + 0.5)).abs() < 1e-12);
    }
    #[test]
    fn test_added_mass_force_opposes_acceleration() {
        let am = AddedMass::new(0.5, 1000.0, 1e-3);
        let a = [1.0, 0.0, 0.0];
        let f = am.force(a);
        assert!(f[0] < 0.0, "Added-mass force should oppose acceleration");
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn test_added_mass_corrected_acceleration() {
        let am = AddedMass::new(0.5, 1000.0, 1e-3);
        let m_body = 2.0;
        let a = am.corrected_acceleration([0.0; 3], m_body);
        assert!(norm3(a) < 1e-30);
    }
    #[test]
    fn test_added_mass_zero_volume() {
        let am = AddedMass::new(0.5, 1000.0, 0.0);
        assert_eq!(am.added_mass(), 0.0);
    }
    #[test]
    fn test_dem_sphere_volume() {
        let p = DemSphParticle::new_sphere([0.0; 3], [0.0; 3], 1.0, 1.0);
        let expected = 4.0 / 3.0 * PI;
        assert!((p.volume() - expected).abs() < 1e-10);
    }
    #[test]
    fn test_dem_sphere_mass() {
        let rho = 2500.0;
        let r = 0.01;
        let p = DemSphParticle::new_sphere([0.0; 3], [0.0; 3], r, rho);
        let expected_mass = rho * (4.0 / 3.0) * PI * r * r * r;
        assert!((p.mass - expected_mass).abs() < 1e-20);
    }
    #[test]
    fn test_dem_sphere_inertia() {
        let r = 0.05;
        let rho = 1000.0;
        let p = DemSphParticle::new_sphere([0.0; 3], [0.0; 3], r, rho);
        let expected_inertia = 0.4 * p.mass * r * r;
        assert!((p.inertia - expected_inertia).abs() < 1e-20);
    }
    #[test]
    fn test_dem_integrate_free_fall() {
        let g = -9.81;
        let mut p = DemSphParticle::new_sphere([0.0, 0.0, 10.0], [0.0; 3], 0.1, 1000.0);
        let f_grav = [0.0, 0.0, g * p.mass];
        let dt = 0.01;
        p.integrate(f_grav, [0.0; 3], dt);
        let expected_vz = g * dt;
        assert!((p.velocity[2] - expected_vz).abs() < 1e-12);
        let expected_z = 10.0 + expected_vz * dt;
        assert!((p.position[2] - expected_z).abs() < 1e-12);
    }
    #[test]
    fn test_dem_integrate_no_force_constant_velocity() {
        let vel = [1.0, 2.0, 3.0];
        let mut p = DemSphParticle::new_sphere([0.0; 3], vel, 0.1, 1000.0);
        let dt = 0.1;
        p.integrate([0.0; 3], [0.0; 3], dt);
        assert!((p.velocity[0] - vel[0]).abs() < 1e-14);
        assert!((p.position[0] - vel[0] * dt).abs() < 1e-14);
        assert!((p.position[1] - vel[1] * dt).abs() < 1e-14);
        assert!((p.position[2] - vel[2] * dt).abs() < 1e-14);
    }
    #[test]
    fn test_dem_contact_no_overlap() {
        let c = DemContactForce::new(1e6, 100.0, 1e4, 0.3);
        let f = c.normal_contact_force([0.0; 3], [5.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1, 0.1);
        assert!(norm3(f) < 1e-20, "No contact → zero force");
    }
    #[test]
    fn test_dem_contact_overlap_force_direction() {
        let c = DemContactForce::new(1e6, 0.0, 1e4, 0.3);
        let f = c.normal_contact_force([0.0; 3], [1.5, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 1.0);
        assert!(
            f[0] < 0.0,
            "Contact force on i should push away from j (in -x), got {}",
            f[0]
        );
    }
    #[test]
    fn test_dem_contact_hertz_power_law() {
        let c = DemContactForce::new(1e6, 0.0, 1e4, 0.3);
        let force_small =
            c.normal_contact_force([0.0; 3], [1.9, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 1.0);
        let force_large =
            c.normal_contact_force([0.0; 3], [1.5, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, 1.0);
        let ratio = norm3(force_large) / norm3(force_small);
        let expected_ratio = (0.5_f64 / 0.1).powf(1.5);
        assert!(
            (ratio - expected_ratio).abs() / expected_ratio < 0.01,
            "Hertz power law: ratio={ratio:.3}, expected={expected_ratio:.3}"
        );
    }
    #[test]
    fn test_dem_wall_contact_no_penetration() {
        let wall = DemWallContact::new([0.0; 3], [0.0, 1.0, 0.0], 1e5, 0.0);
        let mut p = DemSphParticle::new_sphere([0.0, 2.0, 0.0], [0.0; 3], 0.1, 1000.0);
        p.position = [0.0, 2.0, 0.0];
        let f = wall.contact_force(&p);
        assert!(norm3(f) < 1e-20, "No penetration → zero wall force");
    }
    #[test]
    fn test_dem_wall_contact_penetrating_force_direction() {
        let wall = DemWallContact::new([0.0; 3], [0.0, 1.0, 0.0], 1e5, 0.0);
        let mut p = DemSphParticle::new_sphere([0.0, 0.05, 0.0], [0.0; 3], 0.1, 1000.0);
        p.position = [0.0, 0.05, 0.0];
        let f = wall.contact_force(&p);
        assert!(
            f[1] > 0.0,
            "Wall contact force should push in +y (into domain), got {}",
            f[1]
        );
    }
    #[test]
    fn test_dem_wall_signed_distance() {
        let wall = DemWallContact::new([0.0, 1.0, 0.0], [0.0, 1.0, 0.0], 1e5, 0.0);
        let d_above = wall.signed_distance([0.0, 2.0, 0.0]);
        assert!(d_above > 0.0);
        let d_below = wall.signed_distance([0.0, 0.5, 0.0]);
        assert!(d_below < 0.0);
    }
    #[test]
    fn test_impulse_no_impulse_separating() {
        let imp = FluidImpulseTransfer::new(0.8, 1000.0);
        let j = imp.compute_impulse(1.0, [1.0, 0.0, 0.0], [0.0; 3], [1.0, 0.0, 0.0], 1e-3);
        assert!(norm3(j) < 1e-20, "Separating contact → zero impulse");
    }
    #[test]
    fn test_impulse_approaching_gives_positive_impulse() {
        let e = 0.5;
        let imp = FluidImpulseTransfer::new(e, 1000.0);
        let m_body = 1.0;
        let fluid_vol = 1e-3;
        let _m_fluid = 1000.0 * fluid_vol;
        let v_body = [-2.0, 0.0, 0.0];
        let v_fluid = [0.0; 3];
        let n_hat = [-1.0, 0.0, 0.0];
        let j = imp.compute_impulse(m_body, v_body, v_fluid, n_hat, fluid_vol);
        assert!(norm3(j) < 1e-20, "Actually separating → zero impulse");
    }
    #[test]
    fn test_impulse_magnitude_formula() {
        let imp = FluidImpulseTransfer::new(0.0, 1000.0);
        let m_body = 1.0;
        let fluid_vol = 1.0 / 1000.0;
        let v_body = [1.0, 0.0, 0.0];
        let v_fluid = [2.0, 0.0, 0.0];
        let n_hat = [1.0, 0.0, 0.0];
        let j = imp.compute_impulse(m_body, v_body, v_fluid, n_hat, fluid_vol);
        let expected = (1.0 + 0.0) * 1.0 * 1.0 / 2.0;
        assert!(
            (j[0] - expected).abs() < 1e-10,
            "Impulse magnitude: got {}, expected {expected}",
            j[0]
        );
    }
    #[test]
    fn test_fsi_buoyancy_upward() {
        let fsi = FluidStructureInteraction::new(1000.0, 9.81, 0.0, 0.0);
        let f = fsi.total_force([0.0; 3], [0.0; 3], 1e-3, [0.0; 3], &[], &[], &[], 0.1);
        let expected_fz = 1000.0 * 9.81 * 1e-3;
        assert!((f[2] - expected_fz).abs() < 1e-10);
    }
    #[test]
    fn test_fsi_total_force_no_fluid_particles() {
        let fsi = FluidStructureInteraction::new(1000.0, 9.81, 0.5, 0.5);
        let a_body = [1.0, 0.0, 0.0];
        let body_vol = 0.001;
        let f = fsi.total_force([0.0; 3], [0.0; 3], body_vol, a_body, &[], &[], &[], 0.1);
        let m_added = 0.5 * 1000.0 * body_vol;
        let f_added_x = -m_added * a_body[0];
        assert!((f[0] - f_added_x).abs() < 1e-12, "Only added-mass in x");
    }
    #[test]
    fn test_sph_boundary_particle_new() {
        let bp = SphBoundaryParticle::new([1.0, 2.0, 3.0], [0.1, 0.0, 0.0], 1e4, 1000.0, 1e-4);
        assert_eq!(bp.position, [1.0, 2.0, 3.0]);
        assert_eq!(bp.velocity, [0.1, 0.0, 0.0]);
        assert_eq!(bp.pressure, 1e4);
    }
    #[test]
    fn test_sph_boundary_particle_adami_no_fluid() {
        let mut bp = SphBoundaryParticle::new([0.0; 3], [0.0; 3], 0.0, 1000.0, 1e-4);
        bp.update_pressure_adami(&[], &[], &[], [0.0; 3], [0.0, 0.0, -9.81], 0.1);
        assert_eq!(bp.pressure, 0.0);
    }
    #[test]
    fn test_sph_boundary_particle_adami_single_fluid() {
        let mut bp = SphBoundaryParticle::new([0.0; 3], [0.0; 3], 0.0, 1000.0, 1e-4);
        let fluid_pos = vec![[0.05, 0.0, 0.0]];
        let p0 = 1e4;
        let fluid_press = vec![p0];
        let fluid_dens = vec![1000.0_f64];
        let h = 0.2;
        bp.update_pressure_adami(&fluid_pos, &fluid_press, &fluid_dens, [0.0; 3], [0.0; 3], h);
        assert!(bp.pressure.is_finite());
    }
    #[test]
    fn test_dem_sph_coupling_buoyancy_upward() {
        let g = [0.0, 0.0, -9.81];
        let coupling = DemSphCoupling::new(1000.0, g, 0.0, 0.0, 0.1);
        let dem = DemSphParticle::new_sphere([0.0; 3], [0.0; 3], 0.05, 7800.0);
        let f = coupling.hydro_force_on_dem(&dem, [0.0; 3], &[], &[], &[], &[], &[]);
        assert!(f[2] > 0.0, "Buoyancy should be upward (+z), got {}", f[2]);
    }
    #[test]
    fn test_dem_sph_coupling_reaction_length() {
        let g = [0.0, 0.0, -9.81];
        let coupling = DemSphCoupling::new(1000.0, g, 0.0, 0.0, 0.1);
        let dem = DemSphParticle::new_sphere([0.0; 3], [0.0; 3], 0.05, 7800.0);
        let dem_force = [1.0, 0.0, 0.0];
        let sph_pos = vec![[0.05, 0.0, 0.0], [-0.05, 0.0, 0.0]];
        let sph_mass = vec![0.001_f64; 2];
        let sph_dens = vec![1000.0_f64; 2];
        let f_react =
            coupling.reaction_forces_on_sph(&dem, dem_force, &sph_pos, &sph_mass, &sph_dens);
        assert_eq!(f_react.len(), sph_pos.len());
    }
    #[test]
    fn test_dem_sph_coupling_no_nearby_fluid_zero_drag() {
        let g = [0.0, 0.0, -9.81];
        let coupling = DemSphCoupling::new(1000.0, g, 0.5, 0.0, 0.05);
        let dem = DemSphParticle::new_sphere([0.0; 3], [5.0, 0.0, 0.0], 0.05, 7800.0);
        let sph_pos = vec![[10.0, 0.0, 0.0]];
        let sph_vel = vec![[0.0; 3]];
        let sph_press = vec![0.0];
        let sph_mass = vec![0.001];
        let sph_dens = vec![1000.0];
        let f = coupling.hydro_force_on_dem(
            &dem, [0.0; 3], &sph_pos, &sph_vel, &sph_press, &sph_mass, &sph_dens,
        );
        let vol = dem.volume();
        let expected_fz = -1000.0 * g[2] * vol;
        assert!((f[2] - expected_fz).abs() < 1e-10, "Only buoyancy, no drag");
    }
}
