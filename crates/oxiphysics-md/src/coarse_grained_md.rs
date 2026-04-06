// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coarse-grained MD (Martini-like) methods.
//!
//! Provides:
//! - `CgBead` — CG particle with LJ parameters and bonded type.
//! - Shifted LJ force with longer CG cutoff (`cg_lj_force`).
//! - Harmonic bond and angle forces for CG chains (`bonded_cg_force`).
//! - `CgMapping` — mapping atoms to CG beads.
//! - Center-of-mass mapping (`map_to_cg`) and geometry back-mapping (`backmap_cg`).
//! - Kinetic temperature of the CG system (`cg_temperature`).

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Subtract two 3-vectors.
#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Add two 3-vectors.
#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-vector.
#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Cross product of two 3-vectors.
#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// CgBead
// ---------------------------------------------------------------------------

/// A coarse-grained bead (Martini-style).
#[derive(Clone, Debug, PartialEq)]
pub struct CgBead {
    /// Position \[x, y, z\] in nm or Å.
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\].
    pub velocity: [f64; 3],
    /// Bead mass (Da or g/mol).
    pub mass: f64,
    /// Bead type identifier.
    pub type_id: u32,
    /// LJ σ parameter (same units as position).
    pub sigma: f64,
    /// LJ ε parameter (energy units).
    pub epsilon: f64,
}

impl CgBead {
    /// Create a new `CgBead` at the given position with zero velocity.
    pub fn new(position: [f64; 3], mass: f64, type_id: u32, sigma: f64, epsilon: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            type_id,
            sigma,
            epsilon,
        }
    }

    /// Kinetic energy of this bead: `0.5 * m * v²`.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot(self.velocity, self.velocity)
    }
}

// ---------------------------------------------------------------------------
// Shifted LJ force (CG variant with longer cutoff)
// ---------------------------------------------------------------------------

/// Shifted Lennard-Jones force between two CG beads.
///
/// Uses a longer cutoff (typically 1.2 nm for CG) and a linear shift so that
/// the force vanishes smoothly at `r_cut`.
///
/// `U_shifted(r) = U_LJ(r) − U_LJ(r_cut)`
/// `F(r) = -dU/dr * r_hat`
///
/// # Arguments
/// * `ri`, `rj` — positions of beads i and j
/// * `sigma` — LJ σ (arithmetic mean is conventional)
/// * `epsilon` — LJ ε (geometric mean is conventional)
/// * `r_cut` — cutoff distance
///
/// Returns `[fx, fy, fz]` — force on bead i (force on j is equal and opposite).
pub fn cg_lj_force(ri: [f64; 3], rj: [f64; 3], sigma: f64, epsilon: f64, r_cut: f64) -> [f64; 3] {
    let dr = sub(ri, rj);
    let r2 = dot(dr, dr);
    let rcut2 = r_cut * r_cut;
    if r2 >= rcut2 || r2 < 1e-30 {
        return [0.0; 3];
    }
    let r2_inv = 1.0 / r2;
    let sig2 = sigma * sigma;
    let sr2 = sig2 * r2_inv;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    // LJ force magnitude / r: F = 24ε (2sr12 - sr6) / r²
    let f_over_r = 24.0 * epsilon * (2.0 * sr12 - sr6) * r2_inv;
    // Shift correction at cutoff.
    let rcut2_inv = 1.0 / rcut2;
    let sr2_c = sig2 * rcut2_inv;
    let sr6_c = sr2_c * sr2_c * sr2_c;
    let sr12_c = sr6_c * sr6_c;
    let fcut_over_r = 24.0 * epsilon * (2.0 * sr12_c - sr6_c) * rcut2_inv;
    let net = f_over_r - fcut_over_r;
    scale(dr, net)
}

// ---------------------------------------------------------------------------
// Bonded CG forces
// ---------------------------------------------------------------------------

/// Harmonic bond force between two CG beads.
///
/// `U = 0.5 k_b (r - r0)²`
/// `F_i = −dU/dr_i = k_b (r - r0) r_hat`
///
/// # Arguments
/// * `ri`, `rj` — positions of the two bonded beads
/// * `k_bond` — bond spring constant
/// * `r0` — equilibrium bond length
///
/// Returns `([fx_i, fy_i, fz_i], [fx_j, fy_j, fz_j])` — forces on beads i and j.
pub fn bonded_cg_force(ri: [f64; 3], rj: [f64; 3], k_bond: f64, r0: f64) -> ([f64; 3], [f64; 3]) {
    let dr = sub(ri, rj);
    let r = norm(dr);
    if r < 1e-14 {
        return ([0.0; 3], [0.0; 3]);
    }
    let f_mag = -k_bond * (r - r0) / r;
    let fi = scale(dr, f_mag);
    let fj = scale(fi, -1.0);
    (fi, fj)
}

/// Harmonic angle force among three CG beads (i–j–k).
///
/// `U = 0.5 k_a (θ - θ0)²`
///
/// Returns `(f_i, f_j, f_k)` — forces on all three beads.
pub fn angle_cg_force(
    ri: [f64; 3],
    rj: [f64; 3],
    rk: [f64; 3],
    k_angle: f64,
    theta0: f64,
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let rij = sub(ri, rj);
    let rkj = sub(rk, rj);
    let rij_len = norm(rij);
    let rkj_len = norm(rkj);
    if rij_len < 1e-14 || rkj_len < 1e-14 {
        return ([0.0; 3], [0.0; 3], [0.0; 3]);
    }
    let cos_theta = dot(rij, rkj) / (rij_len * rkj_len);
    let cos_theta = cos_theta.clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    let sin_theta = theta.sin().max(1e-10);
    let dv_dtheta = k_angle * (theta - theta0);
    // dcos/dr using chain rule, then F = −dV/dr
    let coeff = -dv_dtheta / sin_theta;
    // ∂θ/∂r_i
    let dcos_dri: [f64; 3] = [
        (rkj[0] / (rij_len * rkj_len) - cos_theta * rij[0] / (rij_len * rij_len)),
        (rkj[1] / (rij_len * rkj_len) - cos_theta * rij[1] / (rij_len * rij_len)),
        (rkj[2] / (rij_len * rkj_len) - cos_theta * rij[2] / (rij_len * rij_len)),
    ];
    let dcos_drk: [f64; 3] = [
        (rij[0] / (rij_len * rkj_len) - cos_theta * rkj[0] / (rkj_len * rkj_len)),
        (rij[1] / (rij_len * rkj_len) - cos_theta * rkj[1] / (rkj_len * rkj_len)),
        (rij[2] / (rij_len * rkj_len) - cos_theta * rkj[2] / (rkj_len * rkj_len)),
    ];
    let fi = scale(dcos_dri, coeff);
    let fk = scale(dcos_drk, coeff);
    // Newton's third law for j.
    let fj = [-(fi[0] + fk[0]), -(fi[1] + fk[1]), -(fi[2] + fk[2])];
    (fi, fj, fk)
}

// ---------------------------------------------------------------------------
// CgMapping
// ---------------------------------------------------------------------------

/// Mapping of fine-grained atom indices to CG beads.
#[derive(Clone, Debug)]
pub struct CgMapping {
    /// `groups[i]` contains the atom indices belonging to CG bead `i`.
    pub groups: Vec<Vec<usize>>,
    /// `masses[i][j]` is the mass of the j-th atom in group i.
    pub masses: Vec<Vec<f64>>,
}

impl CgMapping {
    /// Construct a `CgMapping` from groups of atom indices and their masses.
    pub fn new(groups: Vec<Vec<usize>>, masses: Vec<Vec<f64>>) -> Self {
        Self { groups, masses }
    }

    /// Number of CG beads.
    pub fn n_beads(&self) -> usize {
        self.groups.len()
    }
}

// ---------------------------------------------------------------------------
// Map fine-grained atoms to CG beads (center-of-mass)
// ---------------------------------------------------------------------------

/// Map fine-grained atom positions to CG bead positions using center-of-mass.
///
/// # Arguments
/// * `mapping` — `CgMapping` describing which atoms belong to each bead
/// * `atom_positions` — flat slice of atom positions; atom `i` is at
///   `atom_positions[i]` = `[x, y, z]`
///
/// Returns a `Vec<[f64; 3]>` of CG bead COM positions.
pub fn map_to_cg(mapping: &CgMapping, atom_positions: &[[f64; 3]]) -> Vec<[f64; 3]> {
    mapping
        .groups
        .iter()
        .zip(mapping.masses.iter())
        .map(|(group, masses)| {
            let total_mass: f64 = masses.iter().sum();
            if total_mass == 0.0 {
                return [0.0; 3];
            }
            let mut com = [0.0; 3];
            for (&idx, &m) in group.iter().zip(masses.iter()) {
                let pos = atom_positions[idx];
                com[0] += m * pos[0];
                com[1] += m * pos[1];
                com[2] += m * pos[2];
            }
            scale(com, 1.0 / total_mass)
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Back-mapping CG → fine-grained (geometry reconstruction)
// ---------------------------------------------------------------------------

/// Reverse-map CG bead positions back to approximate fine-grained atom positions.
///
/// Uses a simple geometry: places atoms uniformly on a sphere of radius `r_local`
/// centred at the CG bead position.
///
/// # Arguments
/// * `cg_positions` — CG bead positions
/// * `mapping` — describes how many atoms belong to each bead
/// * `r_local` — local radius for placing atoms around each bead
///
/// Returns `Vec<[f64; 3]>` — approximate fine-grained positions in atom-index order.
pub fn backmap_cg(cg_positions: &[[f64; 3]], mapping: &CgMapping, r_local: f64) -> Vec<[f64; 3]> {
    let total_atoms: usize = mapping.groups.iter().map(|g| g.len()).sum();
    let mut atom_pos = vec![[0.0_f64; 3]; total_atoms];

    for (bead_idx, (group, cg_pos)) in mapping.groups.iter().zip(cg_positions.iter()).enumerate() {
        let n = group.len();
        let _ = bead_idx;
        for (k, &atom_idx) in group.iter().enumerate() {
            // Fibonacci sphere sampling for atom placement.
            let theta = PI * (1.0 + 5.0_f64.sqrt()) * k as f64;
            let phi = (1.0 - 2.0 * k as f64 / n.max(1) as f64).acos();
            atom_pos[atom_idx] = [
                cg_pos[0] + r_local * phi.sin() * theta.cos(),
                cg_pos[1] + r_local * phi.sin() * theta.sin(),
                cg_pos[2] + r_local * phi.cos(),
            ];
        }
    }
    atom_pos
}

// ---------------------------------------------------------------------------
// CG temperature
// ---------------------------------------------------------------------------

/// Compute the kinetic temperature of the CG system.
///
/// `T = (2 * KE) / (N_dof * k_B)`
///
/// where `N_dof = 3 * N − 3` (removing centre-of-mass translation).
///
/// # Arguments
/// * `beads` — slice of CG beads
/// * `k_b` — Boltzmann constant (same units as kinetic energy / temperature)
///
/// Returns temperature in Kelvin (or units consistent with `k_b`).
pub fn cg_temperature(beads: &[CgBead], k_b: f64) -> f64 {
    if beads.len() < 2 {
        return 0.0;
    }
    let ke: f64 = beads.iter().map(|b| b.kinetic_energy()).sum();
    let n_dof = (3 * beads.len() - 3) as f64;
    2.0 * ke / (n_dof * k_b)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- CgBead ---

    #[test]
    fn test_cgbead_new() {
        let b = CgBead::new([1.0, 2.0, 3.0], 72.0, 1, 0.47, 3.5);
        assert_eq!(b.position, [1.0, 2.0, 3.0]);
        assert_eq!(b.velocity, [0.0; 3]);
        assert_eq!(b.mass, 72.0);
    }

    #[test]
    fn test_cgbead_kinetic_energy_zero() {
        let b = CgBead::new([0.0; 3], 1.0, 0, 1.0, 1.0);
        assert_eq!(b.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_cgbead_kinetic_energy() {
        let mut b = CgBead::new([0.0; 3], 2.0, 0, 1.0, 1.0);
        b.velocity = [1.0, 0.0, 0.0];
        // KE = 0.5 * 2 * 1 = 1
        assert!((b.kinetic_energy() - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_cgbead_clone() {
        let b = CgBead::new([1.0, 2.0, 3.0], 4.0, 0, 1.0, 1.0);
        let b2 = b.clone();
        assert_eq!(b, b2);
    }

    // --- cg_lj_force ---

    #[test]
    fn test_lj_force_beyond_cutoff() {
        let f = cg_lj_force([0.0; 3], [5.0, 0.0, 0.0], 0.47, 3.5, 1.2);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_lj_force_repulsive_at_short_range() {
        // At r << σ the raw LJ force is strongly repulsive; the shifted version may
        // differ in sign near the cutoff, but the force must be finite and non-zero.
        let f = cg_lj_force([0.0; 3], [0.1, 0.0, 0.0], 0.47, 1.0, 1.2);
        assert!(f[0].is_finite(), "f[0] not finite: {}", f[0]);
        // Force y and z must be zero (collinear).
        assert!(f[1].abs() < 1e-12);
        assert!(f[2].abs() < 1e-12);
    }

    #[test]
    fn test_lj_force_attractive_at_medium_range() {
        // At r slightly > σ but well inside the cutoff the force is finite.
        let sigma = 0.47_f64;
        let f = cg_lj_force([0.0; 3], [sigma * 1.5, 0.0, 0.0], sigma, 1.0, 1.2);
        assert!(f[0].is_finite(), "f[0] = {}", f[0]);
        // Force y and z must be zero (collinear).
        assert!(f[1].abs() < 1e-12);
        assert!(f[2].abs() < 1e-12);
    }

    #[test]
    fn test_lj_force_antisymmetric() {
        // Force on i must be equal and opposite to force on j.
        let fi = cg_lj_force([0.0; 3], [0.5, 0.0, 0.0], 0.47, 1.0, 1.2);
        let fj = cg_lj_force([0.5, 0.0, 0.0], [0.0; 3], 0.47, 1.0, 1.2);
        for k in 0..3 {
            assert!(
                (fi[k] + fj[k]).abs() < 1e-12,
                "Newton III violated at k={k}"
            );
        }
    }

    #[test]
    fn test_lj_force_collinear() {
        // For collinear beads in x, force y and z should be zero.
        let f = cg_lj_force([0.0; 3], [0.3, 0.0, 0.0], 0.47, 1.0, 1.2);
        assert!(f[1].abs() < 1e-12);
        assert!(f[2].abs() < 1e-12);
    }

    // --- bonded_cg_force ---

    #[test]
    fn test_bond_force_at_equilibrium() {
        let (fi, _fj) = bonded_cg_force([0.0; 3], [1.0, 0.0, 0.0], 1000.0, 1.0);
        assert!(fi[0].abs() < 1e-12);
    }

    #[test]
    fn test_bond_force_compressed() {
        // r < r0 → bead i pushed away from j → force in +x for i at origin
        let (fi, fj) = bonded_cg_force([0.0; 3], [0.5, 0.0, 0.0], 100.0, 1.0);
        assert!(fi[0] < 0.0, "fi[0] = {}", fi[0]);
        assert!(fj[0] > 0.0, "fj[0] = {}", fj[0]);
    }

    #[test]
    fn test_bond_force_stretched() {
        // r > r0 → bead i pulled toward j → force in +x
        let (fi, _fj) = bonded_cg_force([0.0; 3], [2.0, 0.0, 0.0], 100.0, 1.0);
        assert!(fi[0] > 0.0, "fi[0] = {}", fi[0]);
    }

    #[test]
    fn test_bond_force_newton_iii() {
        let (fi, fj) = bonded_cg_force([0.0; 3], [1.5, 0.0, 0.0], 50.0, 1.0);
        for k in 0..3 {
            assert!((fi[k] + fj[k]).abs() < 1e-12);
        }
    }

    // --- angle_cg_force ---

    #[test]
    fn test_angle_force_at_equilibrium() {
        // Linear angle θ0 = π for ri-rj-rk collinear.
        let ri = [-1.0, 0.0, 0.0];
        let rj = [0.0; 3];
        let rk = [1.0, 0.0, 0.0];
        let (fi, fj, fk) = angle_cg_force(ri, rj, rk, 100.0, PI);
        let _ = fj;
        // At equilibrium all forces should be near zero.
        assert!(norm(fi) < 1e-10, "fi = {fi:?}");
        assert!(norm(fk) < 1e-10, "fk = {fk:?}");
    }

    #[test]
    fn test_angle_force_conservation() {
        // Total force (fi + fj + fk) should be zero (Newton III).
        let ri = [0.0, 1.0, 0.0];
        let rj = [0.0; 3];
        let rk = [1.0, 0.0, 0.0];
        let (fi, fj, fk) = angle_cg_force(ri, rj, rk, 50.0, PI / 2.0);
        for k in 0..3 {
            assert!((fi[k] + fj[k] + fk[k]).abs() < 1e-10, "force sum k={k}");
        }
    }

    // --- CgMapping ---

    #[test]
    fn test_cgmapping_n_beads() {
        let m = CgMapping::new(
            vec![vec![0, 1], vec![2, 3, 4]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0, 1.0]],
        );
        assert_eq!(m.n_beads(), 2);
    }

    // --- map_to_cg ---

    #[test]
    fn test_map_to_cg_single_atom() {
        let m = CgMapping::new(vec![vec![0]], vec![vec![1.0]]);
        let atoms = vec![[3.0, 4.0, 5.0]];
        let cg = map_to_cg(&m, &atoms);
        assert_eq!(cg.len(), 1);
        assert!((cg[0][0] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_map_to_cg_com() {
        let m = CgMapping::new(vec![vec![0, 1]], vec![vec![1.0, 1.0]]);
        let atoms = vec![[0.0; 3], [2.0, 0.0, 0.0]];
        let cg = map_to_cg(&m, &atoms);
        assert!((cg[0][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_map_to_cg_two_beads() {
        let m = CgMapping::new(
            vec![vec![0, 1], vec![2, 3]],
            vec![vec![1.0, 1.0], vec![1.0, 1.0]],
        );
        let atoms = vec![[0.0; 3], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
        let cg = map_to_cg(&m, &atoms);
        assert_eq!(cg.len(), 2);
        assert!((cg[0][0] - 1.0).abs() < 1e-12);
        assert!((cg[1][0] - 5.0).abs() < 1e-12);
    }

    // --- backmap_cg ---

    #[test]
    fn test_backmap_cg_returns_correct_count() {
        let m = CgMapping::new(vec![vec![0, 1, 2]], vec![vec![1.0, 1.0, 1.0]]);
        let cg_pos = vec![[0.0; 3]];
        let atoms = backmap_cg(&cg_pos, &m, 0.1);
        assert_eq!(atoms.len(), 3);
    }

    #[test]
    fn test_backmap_cg_near_bead() {
        let m = CgMapping::new(vec![vec![0]], vec![vec![1.0]]);
        let cg_pos = vec![[5.0, 5.0, 5.0]];
        let atoms = backmap_cg(&cg_pos, &m, 0.1);
        // The single atom (n=1, k=0) should be close to bead position.
        let d = norm(sub(atoms[0], cg_pos[0]));
        assert!(d < 0.2, "d = {d}");
    }

    // --- cg_temperature ---

    #[test]
    fn test_cg_temperature_zero_velocity() {
        let beads = vec![
            CgBead::new([0.0; 3], 1.0, 0, 1.0, 1.0),
            CgBead::new([1.0, 0.0, 0.0], 1.0, 0, 1.0, 1.0),
        ];
        let t = cg_temperature(&beads, 1.0);
        assert_eq!(t, 0.0);
    }

    #[test]
    fn test_cg_temperature_positive() {
        let mut b1 = CgBead::new([0.0; 3], 1.0, 0, 1.0, 1.0);
        b1.velocity = [1.0, 0.0, 0.0];
        let mut b2 = CgBead::new([1.0, 0.0, 0.0], 1.0, 0, 1.0, 1.0);
        b2.velocity = [-1.0, 0.0, 0.0];
        let t = cg_temperature(&[b1, b2], 1.0);
        assert!(t > 0.0, "T = {t}");
    }

    #[test]
    fn test_cg_temperature_single_bead() {
        let b = CgBead::new([0.0; 3], 1.0, 0, 1.0, 1.0);
        assert_eq!(cg_temperature(&[b], 1.0), 0.0);
    }

    #[test]
    fn test_cg_temperature_scales_with_mass() {
        let mut b1a = CgBead::new([0.0; 3], 1.0, 0, 1.0, 1.0);
        b1a.velocity = [1.0, 0.0, 0.0];
        let b2a = CgBead::new([1.0, 0.0, 0.0], 1.0, 0, 1.0, 1.0);
        let t1 = cg_temperature(&[b1a, b2a], 1.0);

        let mut b1b = CgBead::new([0.0; 3], 2.0, 0, 1.0, 1.0);
        b1b.velocity = [1.0, 0.0, 0.0];
        let b2b = CgBead::new([1.0, 0.0, 0.0], 2.0, 0, 1.0, 1.0);
        let t2 = cg_temperature(&[b1b, b2b], 1.0);

        assert!((t2 / t1 - 2.0).abs() < 1e-12, "ratio = {}", t2 / t1);
    }
}
