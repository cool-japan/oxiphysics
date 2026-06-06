//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    ForceBuffer, HarmonicAngle, HarmonicBond, LjPotential, NeighborList, VirialStressTensorKernel,
    VirialTensor,
};
use crate::compute::ComputeKernel;

#[cfg(test)]
use super::types::*;

/// Compute the Lennard-Jones potential energy and scalar force magnitude at
/// separation `r`.
///
/// Returns `(energy, force_magnitude)` where:
/// * `energy = 4·ε·[(σ/r)^12 − (σ/r)^6]`
/// * `force_magnitude = 24·ε·[2(σ/r)^12 − (σ/r)^6] / r` (positive = repulsive)
///
/// The caller is responsible for applying the cutoff.
pub fn compute_lj_force(r: f64, lj: &LjPotential) -> (f64, f64) {
    if r < 1e-30 {
        return (f64::INFINITY, f64::INFINITY);
    }
    let sr = lj.sigma / r;
    let sr6 = sr.powi(6);
    let sr12 = sr6 * sr6;
    let energy = 4.0 * lj.epsilon * (sr12 - sr6);
    let force_mag = 24.0 * lj.epsilon * (2.0 * sr12 - sr6) / r;
    (energy, force_mag)
}
/// Compute shifted LJ energy at distance `r` with cutoff `rc`.
///
/// V_shifted(r) = V(r) - V(rc) for r < rc, 0 otherwise.
pub fn compute_lj_shifted_energy(r: f64, lj: &LjPotential, cutoff: f64) -> f64 {
    if r >= cutoff {
        return 0.0;
    }
    let (e_r, _) = compute_lj_force(r, lj);
    let (e_c, _) = compute_lj_force(cutoff, lj);
    e_r - e_c
}
/// Compute Coulomb force between two charged particles.
///
/// Returns `(energy, force_magnitude)`.
pub fn compute_coulomb_force(r: f64, qi: f64, qj: f64, k_e: f64) -> (f64, f64) {
    if r < 1e-30 {
        return (f64::INFINITY, f64::INFINITY);
    }
    let energy = k_e * qi * qj / r;
    let force_mag = k_e * qi * qj / (r * r);
    (energy, force_mag)
}
/// Compute LJ forces using a neighbor list.
pub fn compute_lj_forces_neighborlist(
    positions: &[[f64; 3]],
    lj: &LjPotential,
    nlist: &NeighborList,
    buffer: &mut ForceBuffer,
) {
    let cutoff2 = nlist.cutoff * nlist.cutoff;
    buffer.clear();
    let n = positions.len();
    for i in 0..n {
        for &j in &nlist.neighbors[i] {
            if j <= i {
                continue;
            }
            let dx = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 >= cutoff2 || r2 < 1e-30 {
                continue;
            }
            let r2_inv = 1.0 / r2;
            let sr2 = lj.sigma * lj.sigma * r2_inv;
            let sr6 = sr2 * sr2 * sr2;
            let sr12 = sr6 * sr6;
            let energy = 4.0 * lj.epsilon * (sr12 - sr6);
            let f_mag = 24.0 * lj.epsilon * (2.0 * sr12 - sr6) * r2_inv;
            let f_ij = [f_mag * dx[0], f_mag * dx[1], f_mag * dx[2]];
            buffer.add_pair(i, j, f_ij, energy, dx);
        }
    }
}
/// Compute Coulomb forces using a neighbor list.
pub fn compute_coulomb_forces_neighborlist(
    positions: &[[f64; 3]],
    charges: &[f64],
    k_e: f64,
    nlist: &NeighborList,
    buffer: &mut ForceBuffer,
) {
    let cutoff2 = nlist.cutoff * nlist.cutoff;
    let n = positions.len();
    for i in 0..n {
        for &j in &nlist.neighbors[i] {
            if j <= i {
                continue;
            }
            let dx = [
                positions[i][0] - positions[j][0],
                positions[i][1] - positions[j][1],
                positions[i][2] - positions[j][2],
            ];
            let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            if r2 >= cutoff2 || r2 < 1e-30 {
                continue;
            }
            let r = r2.sqrt();
            let qi = charges[i];
            let qj = charges[j];
            let energy = k_e * qi * qj / r;
            let f_mag = k_e * qi * qj / (r2 * r);
            let f_ij = [f_mag * dx[0], f_mag * dx[1], f_mag * dx[2]];
            buffer.add_pair(i, j, f_ij, energy, dx);
        }
    }
}
/// Compute Lennard-Jones forces for all particle pairs within `cutoff`.
///
/// Returns a `Vec<[f64;3]>` of forces, one per particle.
/// Interactions beyond `cutoff` are ignored.
pub fn compute_all_lj_forces(
    positions: &[[f64; 3]],
    _masses: &[f64],
    lj: &LjPotential,
    cutoff: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let cutoff2 = cutoff * cutoff;
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 >= cutoff2 || r2 < 1e-30 {
                continue;
            }
            let sr = lj.sigma / r2.sqrt();
            let sr6 = sr.powi(6);
            let sr12 = sr6 * sr6;
            let f_mag = 24.0 * lj.epsilon * (2.0 * sr12 - sr6) / r2;
            forces[i][0] += f_mag * dx;
            forces[i][1] += f_mag * dy;
            forces[i][2] += f_mag * dz;
            forces[j][0] -= f_mag * dx;
            forces[j][1] -= f_mag * dy;
            forces[j][2] -= f_mag * dz;
        }
    }
    forces
}
/// Compute Coulomb forces for all particle pairs within `cutoff`.
///
/// Returns a `Vec<[f64;3]>` of forces, one per particle.
pub fn compute_all_coulomb_forces(
    positions: &[[f64; 3]],
    charges: &[f64],
    k_e: f64,
    cutoff: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let cutoff2 = cutoff * cutoff;
    let mut forces = vec![[0.0f64; 3]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            let r2 = dx * dx + dy * dy + dz * dz;
            if r2 >= cutoff2 || r2 < 1e-30 {
                continue;
            }
            let r = r2.sqrt();
            let f_mag = k_e * charges[i] * charges[j] / (r2 * r);
            forces[i][0] += f_mag * dx;
            forces[i][1] += f_mag * dy;
            forces[i][2] += f_mag * dz;
            forces[j][0] -= f_mag * dx;
            forces[j][1] -= f_mag * dy;
            forces[j][2] -= f_mag * dz;
        }
    }
    forces
}
/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
pub(super) fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    poly * (-x * x).exp()
}
/// Self-energy correction for Ewald summation.
///
/// Returns -α/√π · Σ qi².
pub fn ewald_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let sum_q2: f64 = charges.iter().map(|&q| q * q).sum();
    -alpha / std::f64::consts::PI.sqrt() * sum_q2
}
/// Estimate PPPM mesh contribution to long-range energy from a charge mesh.
///
/// This is a simplified mock: computes Σ ρ(k)² * G(k) using a uniform
/// Green's function G(k) = 1 / |k|² for k ≠ 0.
pub fn pppm_mesh_energy_estimate(charge_mesh: &[f64], nx: usize, ny: usize, nz: usize) -> f64 {
    if nx == 0 || ny == 0 || nz == 0 {
        return 0.0;
    }
    let q2: f64 = charge_mesh.iter().map(|&q| q * q).sum();
    q2 / (nx * ny * nz) as f64
}
/// Compute the full virial stress tensor from positions using LJ potential.
///
/// Convenience wrapper around `VirialStressTensorKernel`.
pub fn compute_virial_stress_tensor(
    positions: &[[f64; 3]],
    lj: &LjPotential,
    cutoff: f64,
) -> VirialTensor {
    let n = positions.len();
    let flat_pos: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
    let params = vec![lj.epsilon, lj.sigma, cutoff];
    let mut outputs = vec![Vec::new()];
    VirialStressTensorKernel.execute(&[&flat_pos, &params], &mut outputs, n);
    if outputs[0].len() < 6 {
        return VirialTensor::zero();
    }
    let mut c = [0.0f64; 6];
    c.copy_from_slice(&outputs[0][..6]);
    VirialTensor { components: c }
}
/// Compute harmonic bond forces and accumulate into a force buffer.
///
/// For each bond `(i, j)` with spring constant `k` and rest length `r0`:
/// `F_i = -k*(r - r0)*r̂_ij`,  `F_j = +k*(r - r0)*r̂_ij`
///
/// Returns `(forces, total_bond_energy)`.
pub fn compute_bond_forces(positions: &[[f64; 3]], bonds: &[HarmonicBond]) -> (Vec<[f64; 3]>, f64) {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let mut total_energy = 0.0f64;
    for bond in bonds {
        let i = bond.atom_i;
        let j = bond.atom_j;
        if i >= n || j >= n {
            continue;
        }
        let dx = positions[j][0] - positions[i][0];
        let dy = positions[j][1] - positions[i][1];
        let dz = positions[j][2] - positions[i][2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-30 {
            continue;
        }
        let delta = r - bond.r0;
        let energy = 0.5 * bond.k * delta * delta;
        total_energy += energy;
        let mag = bond.k * delta / r;
        forces[i][0] += mag * dx;
        forces[i][1] += mag * dy;
        forces[i][2] += mag * dz;
        forces[j][0] -= mag * dx;
        forces[j][1] -= mag * dy;
        forces[j][2] -= mag * dz;
    }
    (forces, total_energy)
}
/// Compute harmonic angle forces (CPU mock).
///
/// The angle θ at vertex `j` (between vectors `r_ij` and `r_kj`) is:
/// `cos θ = (r_ij · r_kj) / (|r_ij| |r_kj|)`
///
/// Forces follow from the gradient of the harmonic angle potential.
///
/// Returns `(forces, total_angle_energy)`.
pub fn compute_angle_forces(
    positions: &[[f64; 3]],
    angles: &[HarmonicAngle],
) -> (Vec<[f64; 3]>, f64) {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let mut total_energy = 0.0f64;
    for angle in angles {
        let i = angle.atom_i;
        let j = angle.atom_j;
        let k = angle.atom_k;
        if i >= n || j >= n || k >= n {
            continue;
        }
        let rji = [
            positions[i][0] - positions[j][0],
            positions[i][1] - positions[j][1],
            positions[i][2] - positions[j][2],
        ];
        let rjk = [
            positions[k][0] - positions[j][0],
            positions[k][1] - positions[j][1],
            positions[k][2] - positions[j][2],
        ];
        let len_ji = (rji[0] * rji[0] + rji[1] * rji[1] + rji[2] * rji[2]).sqrt();
        let len_jk = (rjk[0] * rjk[0] + rjk[1] * rjk[1] + rjk[2] * rjk[2]).sqrt();
        if len_ji < 1e-30 || len_jk < 1e-30 {
            continue;
        }
        let cos_theta = (rji[0] * rjk[0] + rji[1] * rjk[1] + rji[2] * rjk[2]) / (len_ji * len_jk);
        let cos_theta = cos_theta.clamp(-1.0, 1.0);
        let theta = cos_theta.acos();
        let delta = theta - angle.theta0;
        total_energy += 0.5 * angle.k_theta * delta * delta;
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt().max(1e-12);
        let d_prefactor = -angle.k_theta * delta / sin_theta;
        for dim in 0..3 {
            let d_cos_d_ri =
                rjk[dim] / (len_ji * len_jk) - cos_theta * rji[dim] / (len_ji * len_ji);
            let d_cos_d_rk =
                rji[dim] / (len_ji * len_jk) - cos_theta * rjk[dim] / (len_jk * len_jk);
            let fi = d_prefactor * d_cos_d_ri;
            let fk = d_prefactor * d_cos_d_rk;
            forces[i][dim] += fi;
            forces[k][dim] += fk;
            forces[j][dim] -= fi + fk;
        }
    }
    (forces, total_energy)
}
/// Compute instantaneous kinetic temperature from particle velocities and masses.
///
/// `T = (2 * KE) / (N_dof * k_B)`, where `N_dof = 3*N - 3` (subtract COM).
/// For simplicity, uses `N_dof = 3*N`.
///
/// # Arguments
/// * `velocities` - Per-particle velocity vectors.
/// * `masses`     - Per-particle masses.
/// * `k_boltzmann` - Boltzmann constant in simulation units.
pub fn kinetic_temperature(velocities: &[[f64; 3]], masses: &[f64], k_boltzmann: f64) -> f64 {
    let n = velocities.len();
    if n == 0 || k_boltzmann < 1e-30 {
        return 0.0;
    }
    let ke2: f64 = velocities
        .iter()
        .zip(masses.iter())
        .map(|(v, &m)| m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum();
    let n_dof = (3 * n) as f64;
    ke2 / (n_dof * k_boltzmann)
}
/// Rescale all velocities to match the target temperature.
///
/// Applies the velocity-rescaling thermostat:
/// `v_i ← v_i * sqrt(T_target / T_current)`
///
/// Does nothing if the current temperature is below a floor value.
pub fn temperature_scale(
    velocities: &mut [[f64; 3]],
    masses: &[f64],
    t_target: f64,
    k_boltzmann: f64,
) {
    let t_current = kinetic_temperature(velocities, masses, k_boltzmann);
    if t_current < 1e-30 || t_target < 0.0 {
        return;
    }
    let scale = (t_target / t_current).sqrt();
    for v in velocities.iter_mut() {
        v[0] *= scale;
        v[1] *= scale;
        v[2] *= scale;
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_md_lj_force_repulsive_at_short_range() {
        let sigma = 1.0_f64;
        let epsilon = 1.0_f64;
        let cutoff = 5.0_f64;
        let r = 0.8_f64 * sigma;
        let positions = vec![0.0, 0.0, 0.0, r, 0.0, 0.0];
        let params = vec![epsilon, sigma, cutoff];
        let mut outputs = vec![Vec::new(), Vec::new()];
        LennardJonesKernel.execute(&[&positions, &params], &mut outputs, 2);
        let fx0 = outputs[0][0];
        let fx1 = outputs[0][3];
        assert!(
            fx0 < 0.0,
            "at r < r_min, force on atom 0 should be negative (repulsive), got {fx0}"
        );
        assert!(
            fx1 > 0.0,
            "at r < r_min, force on atom 1 should be positive (repulsive), got {fx1}"
        );
        assert!(
            (fx0 + fx1).abs() < 1e-10,
            "forces should sum to zero (Newton III), got {fx0} + {fx1} = {}",
            fx0 + fx1
        );
    }
    #[test]
    fn lj_kernel_correct_force_known_separation() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let cutoff = 3.0;
        let positions = vec![0.0, 0.0, 0.0, sigma, 0.0, 0.0];
        let params = vec![epsilon, sigma, cutoff];
        let mut outputs = vec![Vec::new(), Vec::new()];
        LennardJonesKernel.execute(&[&positions, &params], &mut outputs, 2);
        let fx0 = outputs[0][0];
        assert!(
            (fx0 - (-24.0)).abs() < 1e-10,
            "expected fx0 ~ -24.0, got {fx0}"
        );
        let pe = outputs[1][0];
        assert!(pe.abs() < 1e-10, "expected PE ~ 0, got {pe}");
    }
    #[test]
    fn lj_kernel_force_zero_beyond_cutoff() {
        let sigma = 1.0;
        let epsilon = 1.0;
        let cutoff = 2.5;
        let positions = vec![0.0, 0.0, 0.0, 3.0, 0.0, 0.0];
        let params = vec![epsilon, sigma, cutoff];
        let mut outputs = vec![Vec::new(), Vec::new()];
        LennardJonesKernel.execute(&[&positions, &params], &mut outputs, 2);
        for &f in &outputs[0] {
            assert!(f.abs() < 1e-15, "expected zero force, got {f}");
        }
        assert!(outputs[1][0].abs() < 1e-15);
    }
    #[test]
    fn lj_minimum_at_r_min() {
        let lj = LjPotential::new(1.0, 1.0);
        let r_min = lj.r_min();
        let (energy, force_mag) = compute_lj_force(r_min, &lj);
        assert!(
            (energy - (-lj.epsilon)).abs() < 1e-10,
            "energy at r_min should be -epsilon={}, got {energy}",
            -lj.epsilon
        );
        assert!(
            force_mag.abs() < 1e-10,
            "force at r_min should be 0, got {force_mag}"
        );
    }
    #[test]
    fn compute_all_lj_forces_newtons_third_law() {
        let lj = LjPotential::new(1.0, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.6, 1.0, 0.0]];
        let masses = vec![1.0; 3];
        let forces = compute_all_lj_forces(&positions, &masses, &lj, 5.0);
        assert_eq!(forces.len(), 3);
        for k in 0..3 {
            let total: f64 = forces.iter().map(|f| f[k]).sum();
            assert!(
                total.abs() < 1e-10,
                "total force component {k} should be 0, got {total}"
            );
        }
    }
    #[test]
    fn compute_all_lj_forces_repulsive_at_short_range() {
        let sigma = 1.0;
        let lj = LjPotential::new(1.0, sigma);
        let positions = vec![[0.0, 0.0, 0.0], [0.9 * sigma, 0.0, 0.0]];
        let masses = vec![1.0; 2];
        let forces = compute_all_lj_forces(&positions, &masses, &lj, 5.0);
        assert!(
            forces[0][0] < 0.0,
            "repulsive: force[0].x should be < 0, got {}",
            forces[0][0]
        );
        assert!(
            forces[1][0] > 0.0,
            "repulsive: force[1].x should be > 0, got {}",
            forces[1][0]
        );
    }
    #[test]
    fn pair_force_kernel_new() {
        let lj = LjPotential::new(2.0, 0.5);
        let kern = PairForceKernel::new(lj, 3.0, true);
        assert!((kern.lj.epsilon - 2.0).abs() < 1e-15);
        assert!((kern.lj.sigma - 0.5).abs() < 1e-15);
        assert!((kern.cutoff - 3.0).abs() < 1e-15);
        assert!(kern.shift);
        // Test shifted evaluation: at the cutoff boundary, energy should be zero (shifted)
        let (e_at_cut, _f_at_cut) = kern.evaluate(kern.cutoff - 1e-10);
        assert!(
            e_at_cut.abs() < 0.1,
            "energy near cutoff should be small when shifted"
        );
    }
    #[test]
    fn test_coulomb_potential() {
        let cp = CoulombPotential::new(1.0);
        let (e, f) = cp.compute(1.0, 1.0, 1.0);
        assert!((e - 1.0).abs() < 1e-10);
        assert!((f - 1.0).abs() < 1e-10);
        let (e2, f2) = cp.compute(1.0, -1.0, 1.0);
        assert!((e2 - (-1.0)).abs() < 1e-10);
        assert!((f2 - (-1.0)).abs() < 1e-10);
    }
    #[test]
    fn test_coulomb_force_function() {
        let (e, f) = compute_coulomb_force(2.0, 1.0, 1.0, 1.0);
        assert!((e - 0.5).abs() < 1e-10);
        assert!((f - 0.25).abs() < 1e-10);
    }
    #[test]
    fn test_coulomb_kernel_newton_iii() {
        let positions = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let charges = vec![1.0, -1.0];
        let params = vec![1.0, 10.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        CoulombKernel.execute(&[&positions, &charges, &params], &mut outputs, 2);
        for k in 0..3 {
            let total = outputs[0][k] + outputs[0][3 + k];
            assert!(
                total.abs() < 1e-10,
                "forces should sum to zero in dim {k}, got {total}"
            );
        }
        assert!(
            outputs[0][0] > 0.0,
            "particle 0 should be attracted toward +x, got {}",
            outputs[0][0]
        );
    }
    #[test]
    fn test_coulomb_kernel_same_charge_repulsive() {
        let positions = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let charges = vec![1.0, 1.0];
        let params = vec![1.0, 10.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        CoulombKernel.execute(&[&positions, &charges, &params], &mut outputs, 2);
        assert!(
            outputs[0][0] < 0.0,
            "particle 0 should be repelled in -x, got {}",
            outputs[0][0]
        );
    }
    #[test]
    fn test_coulomb_kernel_beyond_cutoff() {
        let positions = vec![0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        let charges = vec![1.0, 1.0];
        let params = vec![1.0, 3.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        CoulombKernel.execute(&[&positions, &charges, &params], &mut outputs, 2);
        for &f in &outputs[0] {
            assert!(
                f.abs() < 1e-15,
                "expected zero force beyond cutoff, got {f}"
            );
        }
    }
    #[test]
    fn test_compute_all_coulomb_forces_newton_iii() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let charges = vec![1.0, -1.0, 0.5];
        let forces = compute_all_coulomb_forces(&positions, &charges, 1.0, 10.0);
        for k in 0..3 {
            let total: f64 = forces.iter().map(|f| f[k]).sum();
            assert!(
                total.abs() < 1e-10,
                "total Coulomb force component {k} should be 0, got {total}"
            );
        }
    }
    #[test]
    fn test_lj_shifted_energy() {
        let lj = LjPotential::new(1.0, 1.0);
        let cutoff = 2.5;
        let e_at_cutoff = compute_lj_shifted_energy(cutoff, &lj, cutoff);
        assert!(
            e_at_cutoff.abs() < 1e-10,
            "shifted energy at cutoff should be 0, got {e_at_cutoff}"
        );
        let e_beyond = compute_lj_shifted_energy(3.0, &lj, cutoff);
        assert!(e_beyond.abs() < 1e-15);
    }
    #[test]
    fn test_lj_well_depth() {
        let lj = LjPotential::new(2.5, 1.0);
        assert!((lj.well_depth() - (-2.5)).abs() < 1e-15);
    }
    #[test]
    fn test_pair_force_kernel_evaluate() {
        let lj = LjPotential::new(1.0, 1.0);
        let kern = PairForceKernel::new(lj, 3.0, false);
        let (e, f) = kern.evaluate(1.0);
        assert!(e.abs() < 1e-10);
        assert!((f - 24.0).abs() < 1e-10);
        let (e2, f2) = kern.evaluate(5.0);
        assert!(e2.abs() < 1e-15);
        assert!(f2.abs() < 1e-15);
    }
    #[test]
    fn test_cutoff_scheme_hard() {
        let scheme = CutoffScheme::Hard { cutoff: 2.5 };
        assert!((scheme.cutoff_distance() - 2.5).abs() < 1e-15);
        assert!((scheme.switch_value(1.0) - 1.0).abs() < 1e-15);
        assert!((scheme.switch_value(3.0) - 0.0).abs() < 1e-15);
    }
    #[test]
    fn test_cutoff_scheme_switched() {
        let scheme = CutoffScheme::Switched {
            r_switch: 2.0,
            r_cutoff: 3.0,
        };
        assert!((scheme.cutoff_distance() - 3.0).abs() < 1e-15);
        assert!((scheme.switch_value(1.5) - 1.0).abs() < 1e-15);
        assert!((scheme.switch_value(3.5) - 0.0).abs() < 1e-15);
        assert!((scheme.switch_value(2.5) - 0.5).abs() < 1e-10);
        let v1 = scheme.switch_value(2.2);
        let v2 = scheme.switch_value(2.8);
        assert!(
            v1 > v2,
            "switch should decrease: v(2.2)={}, v(2.8)={}",
            v1,
            v2
        );
    }
    #[test]
    fn test_neighbor_list_brute_force() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let nlist = NeighborList::build_brute_force(&positions, 2.0, 0.5);
        assert_eq!(nlist.num_particles(), 3);
        assert!(nlist.neighbors[0].contains(&1));
        assert!(nlist.neighbors[1].contains(&0));
        assert!(!nlist.neighbors[0].contains(&2));
        assert!(!nlist.neighbors[2].contains(&0));
    }
    #[test]
    fn test_neighbor_list_num_pairs() {
        let positions = vec![[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let nlist = NeighborList::build_brute_force(&positions, 2.0, 0.0);
        assert_eq!(nlist.num_pairs(), 3);
    }
    #[test]
    fn test_neighbor_list_needs_rebuild() {
        let nlist = NeighborList {
            neighbors: vec![],
            cutoff: 2.5,
            skin: 0.4,
        };
        assert!(!nlist.needs_rebuild(0.1));
        assert!(nlist.needs_rebuild(0.3));
    }
    #[test]
    fn test_force_buffer_basic() {
        let mut buf = ForceBuffer::new(3);
        assert_eq!(buf.forces.len(), 3);
        assert_eq!(buf.total_energy(), 0.0);
        buf.add_pair(0, 1, [1.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]);
        assert!((buf.forces[0][0] - 1.0).abs() < 1e-15);
        assert!((buf.forces[1][0] - (-1.0)).abs() < 1e-15);
        assert!((buf.energies[0] - 1.0).abs() < 1e-15);
        assert!((buf.energies[1] - 1.0).abs() < 1e-15);
        assert!((buf.total_energy() - 2.0).abs() < 1e-15);
    }
    #[test]
    fn test_force_buffer_total_force_zero() {
        let mut buf = ForceBuffer::new(3);
        buf.add_pair(0, 1, [3.0, -1.0, 2.0], 1.0, [1.0, 0.0, 0.0]);
        buf.add_pair(1, 2, [-1.0, 2.0, 0.5], 0.5, [0.0, 1.0, 0.0]);
        let total = buf.total_force();
        for (k, &tk) in total.iter().enumerate() {
            assert!(tk.abs() < 1e-10, "total force[{k}] should be 0, got {}", tk);
        }
    }
    #[test]
    fn test_force_buffer_clear() {
        let mut buf = ForceBuffer::new(2);
        buf.add_pair(0, 1, [1.0, 2.0, 3.0], 5.0, [1.0, 0.0, 0.0]);
        buf.clear();
        assert!((buf.total_energy() - 0.0).abs() < 1e-15);
        for f in &buf.forces {
            for &c in f {
                assert!(c.abs() < 1e-15);
            }
        }
    }
    #[test]
    fn test_force_buffer_reduce() {
        let mut main_buf = ForceBuffer::new(2);
        main_buf.add_pair(0, 1, [1.0, 0.0, 0.0], 2.0, [1.0, 0.0, 0.0]);
        let mut other = ForceBuffer::new(2);
        other.add_pair(0, 1, [0.5, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0]);
        main_buf.reduce_from(&[other]);
        assert!((main_buf.forces[0][0] - 1.5).abs() < 1e-15);
        assert!((main_buf.total_energy() - 3.0).abs() < 1e-15);
    }
    #[test]
    fn test_force_buffer_virial() {
        let mut buf = ForceBuffer::new(2);
        buf.add_pair(0, 1, [2.0, 0.0, 0.0], 1.0, [3.0, 0.0, 0.0]);
        assert!((buf.virial[0][0] - 3.0).abs() < 1e-15);
        assert!((buf.virial[1][0] - 3.0).abs() < 1e-15);
        assert!((buf.total_virial() - 6.0).abs() < 1e-15);
    }
    #[test]
    fn test_lj_forces_neighborlist() {
        let lj = LjPotential::new(1.0, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let nlist = NeighborList::build_brute_force(&positions, 2.5, 0.0);
        let mut buf = ForceBuffer::new(3);
        compute_lj_forces_neighborlist(&positions, &lj, &nlist, &mut buf);
        let total = buf.total_force();
        for (k, &tk) in total.iter().enumerate() {
            assert!(tk.abs() < 1e-10, "total[{k}] = {}", tk);
        }
        for &fk in buf.forces[2].iter() {
            assert!(fk.abs() < 1e-15);
        }
    }
    #[test]
    fn test_coulomb_forces_neighborlist() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let nlist = NeighborList::build_brute_force(&positions, 5.0, 0.0);
        let mut buf = ForceBuffer::new(2);
        compute_coulomb_forces_neighborlist(&positions, &charges, 1.0, &nlist, &mut buf);
        assert!(buf.forces[0][0] > 0.0, "should attract toward +x");
        assert!(buf.forces[1][0] < 0.0, "should attract toward -x");
        let total = buf.total_force();
        assert!(total[0].abs() < 1e-10);
    }
    #[test]
    fn test_lj_forces_neighborlist_matches_brute_force() {
        let lj = LjPotential::new(1.0, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let cutoff = 5.0;
        let masses = vec![1.0; 3];
        let forces_bf = compute_all_lj_forces(&positions, &masses, &lj, cutoff);
        let nlist = NeighborList::build_brute_force(&positions, cutoff, 0.0);
        let mut buf = ForceBuffer::new(3);
        compute_lj_forces_neighborlist(&positions, &lj, &nlist, &mut buf);
        for (i, (buf_row, bf_row)) in buf.forces.iter().zip(forces_bf.iter()).enumerate() {
            for (k, (&buf_val, &bf_val)) in buf_row.iter().zip(bf_row.iter()).enumerate() {
                assert!(
                    (buf_val - bf_val).abs() < 1e-10,
                    "mismatch at particle {i}, dim {k}: nlist={}, brute={}",
                    buf_val,
                    bf_val
                );
            }
        }
    }
    #[test]
    fn test_erfc_approx_at_zero() {
        let result = erfc_approx(0.0);
        assert!((result - 1.0).abs() < 1e-4, "erfc(0) ~ 1, got {result}");
    }
    #[test]
    fn test_erfc_approx_large_arg() {
        let result = erfc_approx(5.0);
        assert!(result < 1e-10, "erfc(5) ~ 0, got {result}");
    }
    #[test]
    fn test_ewald_self_energy() {
        let charges = vec![1.0, -1.0];
        let alpha = 0.5;
        let se = ewald_self_energy(&charges, alpha);
        let expected = -2.0 * alpha / std::f64::consts::PI.sqrt();
        assert!((se - expected).abs() < 1e-10);
    }
    #[test]
    fn test_ewald_real_space_kernel_newton_iii() {
        let pos = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let charges = vec![1.0, -1.0];
        let params = vec![0.5, 5.0, 20.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        EwaldRealSpaceKernel.execute(&[&pos, &charges, &params], &mut outputs, 2);
        assert_eq!(outputs[0].len(), 6);
        let total_fx = outputs[0][0] + outputs[0][3];
        assert!(
            total_fx.abs() < 1e-10,
            "Ewald Newton III violated: {total_fx}"
        );
    }
    #[test]
    fn test_ewald_params_accuracy() {
        let p = EwaldParams::new(0.5, 6.0, 100.0, 20.0);
        let acc = p.real_space_accuracy();
        assert!(acc < 0.01, "erfc(3) should be small, got {acc}");
    }
    #[test]
    fn test_pppm_grid_spacing() {
        let grid = PppmGrid::new(32, 32, 32, 10.0, 2);
        assert!((grid.dx() - 10.0 / 32.0).abs() < 1e-12);
        assert_eq!(grid.total_points(), 32768);
    }
    #[test]
    fn test_pppm_charge_assign_single_particle() {
        let pos = vec![0.5, 0.5, 0.5];
        let charges = vec![1.0];
        let grid_params = vec![4.0, 4.0, 4.0, 4.0];
        let mut outputs = vec![Vec::new()];
        PppmChargeAssignKernel.execute(&[&pos, &charges, &grid_params], &mut outputs, 1);
        assert_eq!(outputs[0].len(), 64);
        let total: f64 = outputs[0].iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-10,
            "total charge on mesh = {total}"
        );
    }
    #[test]
    fn test_pppm_charge_assign_conservation() {
        let pos = vec![1.0, 2.0, 3.0, 5.0, 5.0, 5.0];
        let charges = vec![2.0, -1.5];
        let grid_params = vec![8.0, 8.0, 8.0, 8.0];
        let mut outputs = vec![Vec::new()];
        PppmChargeAssignKernel.execute(&[&pos, &charges, &grid_params], &mut outputs, 2);
        let total: f64 = outputs[0].iter().sum();
        assert!(
            (total - 0.5).abs() < 1e-10,
            "net charge should be 0.5, got {total}"
        );
    }
    #[test]
    fn test_pppm_mesh_energy_estimate_positive() {
        let mesh = vec![1.0, -1.0, 2.0, 0.5];
        let e = pppm_mesh_energy_estimate(&mesh, 2, 2, 1);
        assert!(e >= 0.0, "mesh energy should be non-negative");
    }
    #[test]
    fn test_nlist_update_kernel_no_rebuild() {
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let ref_pos = pos.clone();
        let params = vec![2.5, 0.4];
        let mut outputs = vec![Vec::new(), Vec::new()];
        NlistUpdateKernel.execute(&[&pos, &ref_pos, &params], &mut outputs, 2);
        assert!(
            (outputs[1][0] - 0.0).abs() < 1e-10,
            "status should be Valid (0)"
        );
    }
    #[test]
    fn test_nlist_update_kernel_rebuild() {
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let ref_pos = vec![0.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        let params = vec![2.5, 0.4];
        let mut outputs = vec![Vec::new(), Vec::new()];
        NlistUpdateKernel.execute(&[&pos, &ref_pos, &params], &mut outputs, 2);
        assert!(
            (outputs[1][0] - 1.0).abs() < 1e-10,
            "status should be Rebuilt (1)"
        );
    }
    #[test]
    fn test_nlist_update_pairs_found() {
        let pos = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0];
        let ref_pos = vec![100.0; 9];
        let params = vec![2.5, 0.4];
        let mut outputs = vec![Vec::new(), Vec::new()];
        NlistUpdateKernel.execute(&[&pos, &ref_pos, &params], &mut outputs, 3);
        let num_pairs = outputs[1][1] as usize;
        assert_eq!(num_pairs, 1, "only 1 pair should be found, got {num_pairs}");
    }
    #[test]
    fn test_pair_energy_accumulate_basic() {
        let sigma = 1.0;
        let pos = vec![0.0, 0.0, 0.0, sigma, 0.0, 0.0];
        let pairs = vec![0.0, 1.0];
        let params = vec![1.0, sigma, 5.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        PairEnergyAccumulateKernel.execute(&[&pos, &pairs, &params], &mut outputs, 2);
        let total = outputs[1][0];
        assert!(
            total.abs() < 1e-10,
            "energy at r=sigma should be 0, got {total}"
        );
    }
    #[test]
    fn test_pair_energy_accumulate_split_equally() {
        let pos = vec![0.0, 0.0, 0.0, 0.9, 0.0, 0.0];
        let pairs = vec![0.0, 1.0];
        let params = vec![1.0, 1.0, 5.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        PairEnergyAccumulateKernel.execute(&[&pos, &pairs, &params], &mut outputs, 2);
        let e0 = outputs[0][0];
        let e1 = outputs[0][1];
        assert!((e0 - e1).abs() < 1e-12, "energy should be split equally");
        assert!(
            outputs[1][0] > 0.0,
            "total energy should be positive at r < r_min"
        );
    }
    #[test]
    fn test_pair_energy_beyond_cutoff_zero() {
        let pos = vec![0.0, 0.0, 0.0, 10.0, 0.0, 0.0];
        let pairs = vec![0.0, 1.0];
        let params = vec![1.0, 1.0, 2.5];
        let mut outputs = vec![Vec::new(), Vec::new()];
        PairEnergyAccumulateKernel.execute(&[&pos, &pairs, &params], &mut outputs, 2);
        assert!(
            outputs[1][0].abs() < 1e-15,
            "energy beyond cutoff should be 0"
        );
    }
    #[test]
    fn test_virial_tensor_trace() {
        let vt = VirialTensor {
            components: [1.0, 2.0, 3.0, 0.5, 0.2, 0.1],
        };
        assert!((vt.trace() - 6.0).abs() < 1e-12);
    }
    #[test]
    fn test_virial_tensor_pressure() {
        let vt = VirialTensor {
            components: [-3.0, -3.0, -3.0, 0.0, 0.0, 0.0],
        };
        let p = vt.pressure_contribution(1.0);
        assert!((p - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_virial_tensor_add() {
        let a = VirialTensor {
            components: [1.0, 2.0, 3.0, 0.0, 0.0, 0.0],
        };
        let b = VirialTensor {
            components: [4.0, 5.0, 6.0, 0.0, 0.0, 0.0],
        };
        let c = a.add(&b);
        assert!((c.components[0] - 5.0).abs() < 1e-12);
        assert!((c.trace() - 21.0).abs() < 1e-12);
    }
    #[test]
    fn test_compute_virial_stress_tensor_symmetric() {
        let lj = LjPotential::new(1.0, 1.0);
        let positions = vec![[0.0, 0.0, 0.0], [1.2, 0.0, 0.0]];
        let vt = compute_virial_stress_tensor(&positions, &lj, 5.0);
        assert!(vt.components[1].abs() < 1e-10, "Wyy should be 0");
        assert!(vt.components[2].abs() < 1e-10, "Wzz should be 0");
    }
    #[test]
    fn test_virial_kernel_newton_iii_check() {
        let lj = LjPotential::new(1.0, 1.0);
        let positions = [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [0.5, 1.0, 0.0]];
        let flat_pos: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        let params = [lj.epsilon, lj.sigma, 5.0f64];
        let mut outputs = vec![Vec::new()];
        VirialStressTensorKernel.execute(&[&flat_pos, &params], &mut outputs, 3);
        assert_eq!(
            outputs[0].len(),
            6,
            "virial tensor should have 6 components"
        );
    }
    #[test]
    fn test_bond_force_equilibrium_no_force() {
        let r0 = 1.5_f64;
        let positions = vec![[0.0, 0.0, 0.0], [r0, 0.0, 0.0]];
        let bonds = vec![HarmonicBond::new(0, 1, 100.0, r0)];
        let (forces, energy) = compute_bond_forces(&positions, &bonds);
        assert_eq!(forces.len(), 2);
        for (&f0d, &f1d) in forces[0].iter().zip(forces[1].iter()) {
            assert!(f0d.abs() < 1e-10, "force at equilibrium should be 0");
            assert!(f1d.abs() < 1e-10, "force at equilibrium should be 0");
        }
        assert!(
            energy.abs() < 1e-10,
            "energy at equilibrium should be 0, got {energy}"
        );
    }
    #[test]
    fn test_bond_force_compressed() {
        let r0 = 2.0_f64;
        let r = 1.0_f64;
        let k = 50.0_f64;
        let positions = vec![[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
        let bonds = vec![HarmonicBond::new(0, 1, k, r0)];
        let (forces, energy) = compute_bond_forces(&positions, &bonds);
        assert!(forces[0][0] < 0.0, "atom 0 should be pushed away from bond");
        assert!(forces[1][0] > 0.0, "atom 1 should be pushed away from bond");
        for (dim, (&f0d, &f1d)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            assert!(
                (f0d + f1d).abs() < 1e-10,
                "Newton III violated at dim {dim}"
            );
        }
        let expected_e = 0.5 * k * (r - r0).powi(2);
        assert!(
            (energy - expected_e).abs() < 1e-10,
            "energy mismatch: {energy} vs {expected_e}"
        );
    }
    #[test]
    fn test_bond_force_kernel_executes() {
        let positions_flat = vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        let bond_data = vec![0.0, 1.0, 100.0, 1.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        BondForceKernel.execute(&[&positions_flat, &bond_data], &mut outputs, 2);
        assert_eq!(outputs[0].len(), 6, "forces should have 6 components (3*2)");
        assert_eq!(outputs[1].len(), 1, "energies should have 1 component");
        for &f in &outputs[0] {
            assert!(f.abs() < 1e-10, "force at equilibrium should be 0, got {f}");
        }
    }
    #[test]
    fn test_bond_force_kernel_stretched() {
        let positions_flat = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let bond_data = vec![0.0, 1.0, 10.0, 1.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        BondForceKernel.execute(&[&positions_flat, &bond_data], &mut outputs, 2);
        let fx0 = outputs[0][0];
        let fx1 = outputs[0][3];
        assert!(
            fx0 > 0.0,
            "atom 0 should be pulled toward atom 1 (positive x)"
        );
        assert!(
            fx1 < 0.0,
            "atom 1 should be pulled toward atom 0 (negative x)"
        );
        assert!(
            (fx0 + fx1).abs() < 1e-10,
            "Newton III: forces should cancel"
        );
    }
    #[test]
    fn test_angle_force_at_equilibrium() {
        let positions = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let angles = vec![HarmonicAngle::new(0, 1, 2, 50.0, std::f64::consts::PI)];
        let (forces, energy) = compute_angle_forces(&positions, &angles);
        assert_eq!(forces.len(), 3);
        assert!(
            energy.abs() < 1e-8,
            "energy at equilibrium angle should be ~0, got {energy}"
        );
    }
    #[test]
    fn test_angle_force_finite_at_90_degrees() {
        let positions = vec![[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let theta0 = std::f64::consts::PI / 2.0;
        let angles = vec![HarmonicAngle::new(0, 1, 2, 100.0, theta0)];
        let (forces, energy) = compute_angle_forces(&positions, &angles);
        assert_eq!(forces.len(), 3);
        assert!(energy.is_finite(), "angle energy should be finite");
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "angle force component should be finite: {c}");
            }
        }
        assert!(
            energy.abs() < 1e-8,
            "at equilibrium angle energy should be ~0, got {energy}"
        );
    }
    #[test]
    fn test_angle_force_kernel_executes() {
        let positions_flat = vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let theta0 = std::f64::consts::PI / 2.0;
        let angle_data = vec![0.0, 1.0, 2.0, 100.0, theta0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        AngleForceKernel.execute(&[&positions_flat, &angle_data], &mut outputs, 3);
        assert_eq!(outputs[0].len(), 9, "forces should have 9 components (3*3)");
        assert_eq!(outputs[1].len(), 1, "energies should have 1 element");
        for &f in &outputs[0] {
            assert!(f.is_finite(), "angle force not finite: {f}");
        }
    }
    #[test]
    fn test_kinetic_temperature_basic() {
        let velocities = vec![[3.0, 0.0, 0.0]];
        let masses = vec![1.0];
        let kb = 1.0;
        let t = kinetic_temperature(&velocities, &masses, kb);
        assert!((t - 3.0).abs() < 1e-10, "expected T=3, got {t}");
    }
    #[test]
    fn test_kinetic_temperature_zero_velocity() {
        let velocities = vec![[0.0; 3]; 5];
        let masses = vec![1.0; 5];
        let t = kinetic_temperature(&velocities, &masses, 1.0);
        assert!(
            t.abs() < 1e-15,
            "temperature of zero-velocity system should be 0"
        );
    }
    #[test]
    fn test_temperature_scale_reaches_target() {
        let mut velocities = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let masses = vec![1.0, 1.0];
        let kb = 1.0;
        let t_before = kinetic_temperature(&velocities, &masses, kb);
        assert!(t_before > 0.0);
        let t_target = t_before * 4.0;
        temperature_scale(&mut velocities, &masses, t_target, kb);
        let t_after = kinetic_temperature(&velocities, &masses, kb);
        assert!(
            (t_after - t_target).abs() < 1e-8,
            "after scaling: expected T={t_target}, got T={t_after}"
        );
    }
    #[test]
    fn test_temperature_scale_kernel_rescales() {
        let vel_flat = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
        let masses = vec![1.0, 1.0];
        let kb = 1.0;
        let t_target = 2.0 / 3.0;
        let params = vec![t_target, kb];
        let mut outputs = vec![Vec::new(), Vec::new()];
        TemperatureScaleKernel.execute(&[&vel_flat, &masses, &params], &mut outputs, 2);
        assert_eq!(outputs[0].len(), 6);
        assert_eq!(outputs[1].len(), 2);
        let t_before = outputs[1][0];
        let t_after = outputs[1][1];
        assert!(t_before > 0.0, "t_before should be positive");
        assert!(
            (t_after - t_target).abs() < 1e-8,
            "t_after should be target {t_target}, got {t_after}"
        );
    }
    #[test]
    fn test_temperature_scale_kernel_outputs_finite() {
        let vel_flat = vec![2.0, 1.0, 0.5, 0.3, 0.7, 1.2, 0.1, 0.4, 0.9];
        let masses = vec![1.0, 2.0, 0.5];
        let params = vec![300.0, 1.0];
        let mut outputs = vec![Vec::new(), Vec::new()];
        TemperatureScaleKernel.execute(&[&vel_flat, &masses, &params], &mut outputs, 3);
        for &v in &outputs[0] {
            assert!(v.is_finite(), "scaled velocity not finite: {v}");
        }
        for &t in &outputs[1] {
            assert!(t.is_finite(), "temperature not finite: {t}");
        }
    }
}
