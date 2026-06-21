// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core Coulomb interactions, electric field/potential, PME lattice, and dipole utilities.

use oxifft;

// ---------------------------------------------------------------------------
// Physical constant
// ---------------------------------------------------------------------------

/// Coulomb constant in GROMACS-compatible units: kJ*angstrom*mol^-1*e^-2.
///
/// Value: 138.935458 kJ*angstrom/(mol*e^2).
pub const COULOMB_K: f64 = 138.935_458;

// ---------------------------------------------------------------------------
// Basic Coulomb functions
// ---------------------------------------------------------------------------

/// Direct Coulomb energy (kJ mol^-1) between charges `q_i` and `q_j`
/// separated by distance `r` (angstrom).
///
/// ```text
/// E = K * q_i * q_j / r
/// ```
pub fn coulomb_energy(q_i: f64, q_j: f64, r: f64) -> f64 {
    COULOMB_K * q_i * q_j / r
}

/// Coulomb force magnitude (kJ mol^-1 angstrom^-1) between charges `q_i` and `q_j`
/// separated by distance `r` (angstrom).
///
/// Positive value -> repulsive (same-sign charges).
///
/// ```text
/// F = K * q_i * q_j / r^2
/// ```
pub fn coulomb_force_scalar(q_i: f64, q_j: f64, r: f64) -> f64 {
    COULOMB_K * q_i * q_j / (r * r)
}

// ---------------------------------------------------------------------------
// erfc approximation (Abramowitz & Stegun 7.1.26)
// ---------------------------------------------------------------------------

/// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
///
/// Maximum absolute error ~ 1.5e-7.
pub fn erfc_approx(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc_approx(-x);
    }
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let poly = t
        * (0.254_829_592
            + t * (-0.284_496_736
                + t * (1.421_413_741 + t * (-1.453_152_027 + t * 1.061_405_429))));
    poly * (-x * x).exp()
}

// ---------------------------------------------------------------------------
// Ewald functions
// ---------------------------------------------------------------------------

/// Ewald real-space (direct) pair energy (kJ mol^-1) using erfc damping.
///
/// ```text
/// E = K * q_i * q_j * erfc(alpha * r) / r
/// ```
pub fn ewald_direct_energy(q_i: f64, q_j: f64, r: f64, alpha: f64) -> f64 {
    COULOMB_K * q_i * q_j * erfc_approx(alpha * r) / r
}

/// Ewald self-interaction energy (kJ mol^-1).
///
/// Removes the spurious self-interaction introduced by the Ewald splitting:
/// ```text
/// E_self = -K * alpha / sqrt(pi) * sum qi^2
/// ```
pub fn ewald_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
    -COULOMB_K * alpha / std::f64::consts::PI.sqrt() * sum_q2
}

// ---------------------------------------------------------------------------
// PmeLattice
// ---------------------------------------------------------------------------

/// A simple PME charge grid for spreading charges onto a regular lattice.
///
/// Provides nearest-grid-point (NGP) charge spreading and an FFT-based
/// reciprocal-space Ewald energy estimate
/// ([`PmeLattice::reciprocal_energy_approx`]). The zeroth-order NGP assignment
/// makes the estimate approximate; the smooth-PME variant with B-spline
/// interpolation lives in [`crate::electrostatics::pme`].
pub struct PmeLattice {
    /// Number of grid points along x.
    pub nx: usize,
    /// Number of grid points along y.
    pub ny: usize,
    /// Number of grid points along z.
    pub nz: usize,
    /// Flattened grid values (row-major: index = ix*ny*nz + iy*nz + iz).
    pub grid: Vec<f64>,
    /// Box lengths \[Lx, Ly, Lz\] in angstrom.
    pub box_lengths: [f64; 3],
}

impl PmeLattice {
    /// Create a new zero-initialised [`PmeLattice`].
    pub fn new(nx: usize, ny: usize, nz: usize, box_lengths: [f64; 3]) -> Self {
        Self {
            nx,
            ny,
            nz,
            grid: vec![0.0; nx * ny * nz],
            box_lengths,
        }
    }

    /// Map a real-space position to a grid-cell index using fractional coordinates.
    ///
    /// Coordinates are wrapped periodically into \[0, L) before mapping.
    pub fn cell_index(&self, pos: [f64; 3]) -> (usize, usize, usize) {
        let ix = {
            let frac = pos[0].rem_euclid(self.box_lengths[0]) / self.box_lengths[0];
            ((frac * self.nx as f64) as usize).min(self.nx - 1)
        };
        let iy = {
            let frac = pos[1].rem_euclid(self.box_lengths[1]) / self.box_lengths[1];
            ((frac * self.ny as f64) as usize).min(self.ny - 1)
        };
        let iz = {
            let frac = pos[2].rem_euclid(self.box_lengths[2]) / self.box_lengths[2];
            ((frac * self.nz as f64) as usize).min(self.nz - 1)
        };
        (ix, iy, iz)
    }

    /// Spread charges onto the grid using nearest-grid-point (NGP) assignment.
    ///
    /// Each charge is assigned entirely to the single nearest grid cell.
    pub fn spread_charges(&mut self, positions: &[[f64; 3]], charges: &[f64]) {
        for (pos, &q) in positions.iter().zip(charges.iter()) {
            let (ix, iy, iz) = self.cell_index(*pos);
            self.grid[ix * self.ny * self.nz + iy * self.nz + iz] += q;
        }
    }

    /// Zero all grid values.
    pub fn clear(&mut self) {
        for v in self.grid.iter_mut() {
            *v = 0.0;
        }
    }

    /// Reciprocal-space Ewald energy (kJ mol⁻¹) from a 3-D FFT of the charge grid.
    ///
    /// Performs an OxiFFT 3-D FFT of the charge density already spread onto this
    /// lattice and evaluates the standard Ewald reciprocal-space sum
    ///
    /// ```text
    /// E_recip = K/(2V) Σ_{k≠0} (4π/k²) exp(−k²/4α²) |ρ̂(k)|²
    /// ```
    ///
    /// where ρ̂(k) is the discrete Fourier transform of the grid and the
    /// splitting parameter α defaults to 5/min(L) as a safe heuristic.
    ///
    /// This is an *approximate* PME: charges are assigned with nearest-grid-point
    /// (zeroth-order) interpolation by [`Self::spread_charges`] and no B-spline
    /// influence-function deconvolution is applied, so it converges to the true
    /// reciprocal energy only on a sufficiently fine grid. For the smooth-PME
    /// variant with B-spline interpolation and deconvolution use
    /// [`crate::electrostatics::pme::pme_reciprocal_energy`]. The returned value
    /// is a genuine FFT evaluation, not the former `sum(grid²)·V/N²` proxy.
    pub fn reciprocal_energy_approx(&self) -> f64 {
        let n_total = self.nx * self.ny * self.nz;
        let volume = self.box_lengths[0] * self.box_lengths[1] * self.box_lengths[2];
        // Default alpha: 5.0 / min(box_length) as a safe heuristic
        let alpha = 5.0
            / self.box_lengths[0]
                .min(self.box_lengths[1])
                .min(self.box_lengths[2]);
        let four_alpha_sq = 4.0 * alpha * alpha;
        let two_pi = 2.0 * std::f64::consts::PI;

        // 3D FFT of the pre-spread charge grid
        let rho_imag = vec![0.0f64; n_total];
        let (rho_re, rho_im) =
            oxifft::fft3d_split::<f64>(&self.grid, &rho_imag, self.nx, self.ny, self.nz);

        let prefactor = COULOMB_K / (2.0 * volume);
        let mut energy = 0.0f64;

        for ix in 0..self.nx {
            let nx = if ix <= self.nx / 2 {
                ix as i64
            } else {
                ix as i64 - self.nx as i64
            };
            let kx = two_pi * nx as f64 / self.box_lengths[0];
            for iy in 0..self.ny {
                let ny = if iy <= self.ny / 2 {
                    iy as i64
                } else {
                    iy as i64 - self.ny as i64
                };
                let ky = two_pi * ny as f64 / self.box_lengths[1];
                for iz in 0..self.nz {
                    if ix == 0 && iy == 0 && iz == 0 {
                        continue;
                    }
                    let nz = if iz <= self.nz / 2 {
                        iz as i64
                    } else {
                        iz as i64 - self.nz as i64
                    };
                    let kz = two_pi * nz as f64 / self.box_lengths[2];
                    let k2 = kx * kx + ky * ky + kz * kz;
                    if k2 < 1e-20 {
                        continue;
                    }
                    let gaussian = (-k2 / four_alpha_sq).exp();
                    let flat = ix * self.ny * self.nz + iy * self.nz + iz;
                    let rho_sq = rho_re[flat] * rho_re[flat] + rho_im[flat] * rho_im[flat];
                    energy += (4.0 * std::f64::consts::PI / k2) * gaussian * rho_sq;
                }
            }
        }
        prefactor * energy
    }

    /// Total charge on the grid.
    pub fn total_charge(&self) -> f64 {
        self.grid.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// reaction_field_correction — standalone function
// ---------------------------------------------------------------------------

/// Reaction-field correction energy (kJ mol^-1) for a pair of charges.
///
/// Adds a continuum dielectric correction beyond the cutoff radius using
/// the Onsager reaction-field model:
///
/// ```text
/// E_rf = K * qi * qj * (k_rf * r^2 - c_rf)
/// ```
///
/// where:
/// - k_rf = (epsilon_r - 1) / (2*epsilon_r + 1) / cutoff^3
/// - c_rf = 3*epsilon_r / (2*epsilon_r + 1) / cutoff
///
/// Returns 0 if `r >= cutoff`.
pub fn reaction_field_correction(qi: f64, qj: f64, r: f64, cutoff: f64, epsilon_r: f64) -> f64 {
    if r >= cutoff || r <= 0.0 {
        return 0.0;
    }
    let k_rf = (epsilon_r - 1.0) / (2.0 * epsilon_r + 1.0) / (cutoff * cutoff * cutoff);
    let c_rf = 3.0 * epsilon_r / (2.0 * epsilon_r + 1.0) / cutoff;
    COULOMB_K * qi * qj * (k_rf * r * r - c_rf)
}

// ---------------------------------------------------------------------------
// coulomb_energy_direct — O(N²) direct Coulomb sum
// ---------------------------------------------------------------------------

/// O(N²) direct Coulomb energy (kJ mol^-1) without any cutoff or periodic
/// boundary conditions.
///
/// ```text
/// E = COULOMB_K * Σ_{i<j} qᵢ qⱼ / rᵢⱼ
/// ```
///
/// # Arguments
/// * `positions` - atom positions (angstrom).
/// * `charges`   - atom charges (e).
pub fn coulomb_energy_direct(positions: &[[f64; 3]], charges: &[f64]) -> f64 {
    let n = positions.len();
    assert_eq!(charges.len(), n, "positions and charges length mismatch");
    let mut energy = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r > 1e-20 {
                energy += COULOMB_K * charges[i] * charges[j] / r;
            }
        }
    }
    energy
}

/// O(N²) direct Coulomb forces (kJ mol^-1 angstrom^-1) without cutoff.
///
/// Returns a `Vec` of force vectors, one per atom.
/// Sign convention: for same-sign charges, the force on atom i points away
/// from atom j (repulsive); for opposite-sign charges, force on atom i
/// points toward atom j (attractive).
pub fn coulomb_forces_direct(positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
    let n = positions.len();
    assert_eq!(charges.len(), n, "positions and charges length mismatch");
    let mut forces = vec![[0.0_f64; 3]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            // dr points from i to j
            let dr = [
                positions[j][0] - positions[i][0],
                positions[j][1] - positions[i][1],
                positions[j][2] - positions[i][2],
            ];
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if r2 < 1e-40 {
                continue;
            }
            let r = r2.sqrt();
            // F_mag = K * qi * qj / r^2 (positive = repulsive for same-sign)
            // Force on i: F_mag acts along -dr/r for repulsion (away from j)
            // F_i = -K*qi*qj/r^2 * dr/r  => F_i[a] = -f_mag * dr[a]/r
            let f_mag = COULOMB_K * charges[i] * charges[j] / r2;
            for a in 0..3 {
                // Positive f_mag & same sign: force on i in -dr direction (away from j)
                let f_a = -f_mag * dr[a] / r;
                forces[i][a] += f_a;
                forces[j][a] -= f_a; // Newton III: force on j in +dr direction
            }
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// PointDipoleInteraction — pair interaction between two point dipoles
// ---------------------------------------------------------------------------

/// Compute the interaction energy between two point dipoles `d1` and `d2`
/// separated by displacement `r_vec` (Å).
///
/// E = (1/r³) * \[d1·d2 - 3*(d1·r̂)*(d2·r̂)\] * COULOMB_K
pub fn dipole_dipole_energy(d1: [f64; 3], d2: [f64; 3], r_vec: [f64; 3]) -> f64 {
    let r2: f64 = r_vec.iter().map(|&x| x * x).sum();
    if r2 < 1e-20 {
        return 0.0;
    }
    let r = r2.sqrt();
    let r3 = r2 * r;
    let d1_dot_d2: f64 = (0..3).map(|a| d1[a] * d2[a]).sum();
    let d1_dot_r: f64 = (0..3).map(|a| d1[a] * r_vec[a]).sum::<f64>() / r;
    let d2_dot_r: f64 = (0..3).map(|a| d2[a] * r_vec[a]).sum::<f64>() / r;
    COULOMB_K * (d1_dot_d2 - 3.0_f64 * d1_dot_r * d2_dot_r) / r3
}

/// Compute the force on dipole 1 due to dipole 2 (kJ mol⁻¹ Å⁻¹).
///
/// Returns `[f64; 3]` force vector on dipole 1.
pub fn dipole_dipole_force(d1: [f64; 3], d2: [f64; 3], r_vec: [f64; 3]) -> [f64; 3] {
    let r2: f64 = r_vec.iter().map(|&x| x * x).sum();
    if r2 < 1e-20 {
        return [0.0; 3];
    }
    let r = r2.sqrt();
    let r5 = r2 * r2 * r;
    let d1_dot_d2: f64 = (0..3).map(|a| d1[a] * d2[a]).sum();
    let d1_dot_r: f64 = (0..3).map(|a| d1[a] * r_vec[a]).sum();
    let d2_dot_r: f64 = (0..3).map(|a| d2[a] * r_vec[a]).sum();
    let mut f = [0.0_f64; 3];
    for a in 0..3 {
        f[a] = COULOMB_K / r5
            * (3.0_f64 * (d1[a] * d2_dot_r + d2[a] * d1_dot_r) + 3.0_f64 * d1_dot_d2 * r_vec[a]
                - 15.0_f64 * d1_dot_r * d2_dot_r * r_vec[a] / r2);
    }
    f
}

// ---------------------------------------------------------------------------
// ElectricField — field from a collection of point charges
// ---------------------------------------------------------------------------

/// Compute the electric field (kJ mol⁻¹ Å⁻¹ e⁻¹) at a set of evaluation points
/// due to a collection of point charges.
///
/// Each field vector is `E_a = COULOMB_K * q / r³ * r_a`.
pub fn electric_field(
    source_positions: &[[f64; 3]],
    charges: &[f64],
    eval_points: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let ns = source_positions.len().min(charges.len());
    eval_points
        .iter()
        .map(|&ep| {
            let mut e_field = [0.0_f64; 3];
            for j in 0..ns {
                let dr = [
                    ep[0] - source_positions[j][0],
                    ep[1] - source_positions[j][1],
                    ep[2] - source_positions[j][2],
                ];
                let r2: f64 = dr.iter().map(|&x| x * x).sum();
                if r2 < 1e-20 {
                    continue;
                }
                let r3 = r2 * r2.sqrt();
                for a in 0..3 {
                    e_field[a] += COULOMB_K * charges[j] * dr[a] / r3;
                }
            }
            e_field
        })
        .collect()
}

/// Compute the electric potential (kJ mol⁻¹ e⁻¹) at a set of evaluation points.
pub fn electric_potential(
    source_positions: &[[f64; 3]],
    charges: &[f64],
    eval_points: &[[f64; 3]],
) -> Vec<f64> {
    let ns = source_positions.len().min(charges.len());
    eval_points
        .iter()
        .map(|&ep| {
            (0..ns)
                .map(|j| {
                    let r2: f64 = (0..3)
                        .map(|a| (ep[a] - source_positions[j][a]).powi(2))
                        .sum();
                    if r2 < 1e-20 {
                        0.0
                    } else {
                        COULOMB_K * charges[j] / r2.sqrt()
                    }
                })
                .sum()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// ScreenedInteractionMatrix — pairwise screened interaction table
// ---------------------------------------------------------------------------

/// Precomputes a pairwise screened interaction matrix for fast energy lookups.
///
/// `W[i][j] = COULOMB_K * q_i * q_j * exp(-kappa * r_ij) / r_ij`
#[derive(Debug, Clone)]
pub struct ScreenedInteractionMatrix {
    /// Interaction energies (upper triangle stored).
    pub w: Vec<Vec<f64>>,
    /// Number of particles.
    pub n: usize,
}

impl ScreenedInteractionMatrix {
    /// Build the interaction matrix for `n` particles.
    ///
    /// `kappa` is the inverse Debye length (Å⁻¹).
    pub fn build(positions: &[[f64; 3]], charges: &[f64], kappa: f64) -> Self {
        let n = positions.len().min(charges.len());
        let mut w = vec![vec![0.0_f64; n]; n];
        for i in 0..n {
            for j in (i + 1)..n {
                let r2: f64 = (0..3)
                    .map(|a| (positions[j][a] - positions[i][a]).powi(2))
                    .sum();
                if r2 < 1e-20 {
                    continue;
                }
                let r = r2.sqrt();
                let e = COULOMB_K * charges[i] * charges[j] * (-kappa * r).exp() / r;
                w[i][j] = e;
                w[j][i] = e;
            }
        }
        Self { w, n }
    }

    /// Total interaction energy (kJ mol⁻¹).
    pub fn total_energy(&self) -> f64 {
        let mut e = 0.0_f64;
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                e += self.w[i][j];
            }
        }
        e
    }

    /// Interaction energy of atom `i` with all other atoms.
    pub fn atom_energy(&self, i: usize) -> f64 {
        if i >= self.n {
            return 0.0;
        }
        self.w[i].iter().sum()
    }
}

// ---------------------------------------------------------------------------
// Additional Ewald utilities
// ---------------------------------------------------------------------------

/// Compute the Ewald reciprocal-space energy for a single k-vector.
///
/// E_k = (2*pi*K / V) * |S(k)|^2 * exp(-k^2/(4*alpha^2)) / k^2
///
/// where S(k) = sum_i q_i * exp(i * k . r_i) is the structure factor.
///
/// Since we work with real arithmetic, we compute |S(k)|^2 = Re(S)^2 + Im(S)^2.
pub fn ewald_reciprocal_k_energy(
    positions: &[[f64; 3]],
    charges: &[f64],
    k_vec: [f64; 3],
    alpha: f64,
    volume: f64,
) -> f64 {
    let k2 = k_vec[0] * k_vec[0] + k_vec[1] * k_vec[1] + k_vec[2] * k_vec[2];
    if k2 < 1e-20 {
        return 0.0;
    }

    let mut s_re = 0.0_f64;
    let mut s_im = 0.0_f64;
    for (pos, &q) in positions.iter().zip(charges.iter()) {
        let phase = k_vec[0] * pos[0] + k_vec[1] * pos[1] + k_vec[2] * pos[2];
        s_re += q * phase.cos();
        s_im += q * phase.sin();
    }
    let s_sq = s_re * s_re + s_im * s_im;
    let exp_factor = (-k2 / (4.0 * alpha * alpha)).exp();
    2.0 * std::f64::consts::PI * COULOMB_K * s_sq * exp_factor / (k2 * volume)
}

/// Ewald reciprocal-space energy summed over a grid of k-vectors.
///
/// Sums over vectors k = (nx, ny, nz) * 2*pi/L for -n_max <= n <= n_max.
pub fn ewald_reciprocal_energy_grid(
    positions: &[[f64; 3]],
    charges: &[f64],
    box_lengths: [f64; 3],
    alpha: f64,
    n_max: i32,
) -> f64 {
    let volume = box_lengths[0] * box_lengths[1] * box_lengths[2];
    let mut energy = 0.0;
    for nx in -n_max..=n_max {
        for ny in -n_max..=n_max {
            for nz in -n_max..=n_max {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                let kx = 2.0 * std::f64::consts::PI * nx as f64 / box_lengths[0];
                let ky = 2.0 * std::f64::consts::PI * ny as f64 / box_lengths[1];
                let kz = 2.0 * std::f64::consts::PI * nz as f64 / box_lengths[2];
                energy +=
                    ewald_reciprocal_k_energy(positions, charges, [kx, ky, kz], alpha, volume);
            }
        }
    }
    energy
}

// ---------------------------------------------------------------------------
// Additional electrostatic utilities
// ---------------------------------------------------------------------------

/// Compute electric potential energy of a dipole mu in external field E.
///
/// U = -mu . E
pub fn dipole_in_field_energy(mu: [f64; 3], e_field: [f64; 3]) -> f64 {
    -(mu[0] * e_field[0] + mu[1] * e_field[1] + mu[2] * e_field[2])
}

/// Torque on dipole mu in electric field E: tau = mu × E.
pub fn dipole_torque(mu: [f64; 3], e_field: [f64; 3]) -> [f64; 3] {
    [
        mu[1] * e_field[2] - mu[2] * e_field[1],
        mu[2] * e_field[0] - mu[0] * e_field[2],
        mu[0] * e_field[1] - mu[1] * e_field[0],
    ]
}

/// Compute the total dipole moment of a charge distribution.
pub fn total_dipole_moment(positions: &[[f64; 3]], charges: &[f64]) -> [f64; 3] {
    let mut mu = [0.0_f64; 3];
    for (pos, &q) in positions.iter().zip(charges.iter()) {
        mu[0] += q * pos[0];
        mu[1] += q * pos[1];
        mu[2] += q * pos[2];
    }
    mu
}

/// Lennard-Jones + Coulomb combined pair energy.
///
/// E = eps * \[(sigma/r)^12 - 2*(sigma/r)^6\] + K*qi*qj/r
pub fn lj_coulomb_pair_energy(qi: f64, qj: f64, r: f64, eps: f64, sigma: f64) -> f64 {
    if r < 1e-10 {
        return 0.0;
    }
    let sr = sigma / r;
    let sr6 = sr.powi(6);
    let lj = eps * (sr6 * sr6 - 2.0 * sr6);
    let coul = COULOMB_K * qi * qj / r;
    lj + coul
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electrostatics::*;

    #[test]
    fn test_coulomb_energy_unit() {
        let e = coulomb_energy(1.0, 1.0, 1.0);
        assert!(
            (e - COULOMB_K).abs() < 1e-6,
            "coulomb_energy(1,1,1) = {e}, expected {COULOMB_K}"
        );
    }

    #[test]
    fn test_erfc_approx_zero() {
        let v = erfc_approx(0.0);
        assert!(
            (v - 1.0).abs() < 1e-6,
            "erfc_approx(0) = {v}, expected ~1.0"
        );
    }

    #[test]
    fn test_erfc_approx_large() {
        let v = erfc_approx(3.0);
        assert!(v < 1e-4, "erfc_approx(3) = {v}, expected ~0");
    }

    #[test]
    fn test_ewald_self_energy_single_charge() {
        let alpha = 0.3;
        let e = ewald_self_energy(&[1.0], alpha);
        assert!(e < 0.0, "ewald_self_energy should be negative, got {e}");
        let expected = -COULOMB_K * alpha / std::f64::consts::PI.sqrt();
        assert!(
            (e - expected).abs() < 1e-8,
            "ewald_self_energy = {e}, expected {expected}"
        );
    }

    #[test]
    fn test_pme_lattice_spread_and_clear() {
        let mut lattice = PmeLattice::new(4, 4, 4, [10.0, 10.0, 10.0]);
        lattice.spread_charges(&[[2.5, 2.5, 2.5]], &[1.0]);
        let any_nonzero = lattice.grid.iter().any(|&v| v != 0.0);
        assert!(
            any_nonzero,
            "grid should have a non-zero cell after spreading"
        );
        lattice.clear();
        let all_zero = lattice.grid.iter().all(|&v| v == 0.0);
        assert!(all_zero, "grid should be all-zero after clear");
    }

    #[test]
    fn test_pme_reciprocal_energy_real_fft() {
        // The reciprocal energy must be a genuine FFT evaluation: finite,
        // non-zero for a charged distribution, and distinct from the former
        // sum(grid^2)*V/N^2 proxy that it replaced.
        let mut lattice = PmeLattice::new(8, 8, 8, [12.0, 12.0, 12.0]);
        lattice.spread_charges(&[[2.0, 2.0, 2.0], [8.0, 8.0, 8.0]], &[1.0, -1.0]);
        let e = lattice.reciprocal_energy_approx();
        assert!(e.is_finite());
        assert!(e.abs() > 1e-12, "reciprocal energy should be non-zero");

        let n = (lattice.nx * lattice.ny * lattice.nz) as f64;
        let volume = 12.0 * 12.0 * 12.0;
        let proxy: f64 = lattice.grid.iter().map(|g| g * g).sum::<f64>() * volume / (n * n);
        assert!(
            (e - proxy).abs() > 1e-9,
            "must not reproduce the old proxy formula"
        );

        // A cleared (neutral, empty) grid yields exactly zero reciprocal energy.
        lattice.clear();
        assert_eq!(lattice.reciprocal_energy_approx(), 0.0);
    }

    #[test]
    fn test_debye_potential_decays() {
        let dh = DebyeHuckelModel::new(0.5, 78.4);
        let v1 = dh.potential(2.0, 1.0);
        let v2 = dh.potential(4.0, 1.0);
        assert!(
            v1.abs() > v2.abs(),
            "potential should decay with distance: |V(2)| = {} > |V(4)| = {}",
            v1.abs(),
            v2.abs()
        );
    }

    #[test]
    fn test_debye_force_repulsive_same_sign() {
        let dh = DebyeHuckelModel::new(0.3, 78.4);
        let f = dh.force([3.0, 0.0, 0.0], 1.0, 1.0);
        assert!(
            f[0] > 0.0,
            "same-sign force should be repulsive (+x), got {}",
            f[0]
        );
    }

    #[test]
    fn test_pme_lattice_total_charge() {
        let mut lattice = PmeLattice::new(4, 4, 4, [10.0, 10.0, 10.0]);
        lattice.spread_charges(&[[2.5, 2.5, 2.5], [7.5, 7.5, 7.5]], &[1.0, -1.0]);
        let total = lattice.total_charge();
        assert!(
            total.abs() < 1e-12,
            "total charge should be 0 for neutral system, got {total}"
        );
    }

    #[test]
    fn test_coulomb_energy_direct_unit_charges() {
        // Two unit charges at 1 angstrom: E = COULOMB_K
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let e = coulomb_energy_direct(&positions, &charges);
        assert!(
            (e - COULOMB_K).abs() < 1e-6,
            "Direct Coulomb at r=1: got {e}, expected {COULOMB_K}"
        );
    }

    #[test]
    fn test_coulomb_energy_direct_attractive() {
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = coulomb_energy_direct(&positions, &charges);
        assert!(
            e < 0.0,
            "opposite charges should have negative direct energy, got {e}"
        );
    }

    #[test]
    fn test_coulomb_energy_direct_scales_with_1_over_r() {
        // At r=2 the energy should be half that at r=1 (for same charges)
        let pos1 = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let pos2 = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let e1 = coulomb_energy_direct(&pos1, &charges);
        let e2 = coulomb_energy_direct(&pos2, &charges);
        assert!(
            (e2 * 2.0 - e1).abs() < 1e-8,
            "1/r scaling: e1={e1}, 2*e2={}",
            2.0 * e2
        );
    }

    #[test]
    fn test_coulomb_energy_direct_three_charges() {
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let charges = [1.0, 1.0, 1.0];
        let e = coulomb_energy_direct(&positions, &charges);
        // Pairs (0,1), (0,2), (1,2) at r=1,2,1
        // E = K*(1/1 + 1/2 + 1/1) = K * 2.5
        let expected = COULOMB_K * 2.5;
        assert!(
            (e - expected).abs() < 1e-6,
            "three-charge direct energy: got {e}, expected {expected}"
        );
    }

    #[test]
    fn test_coulomb_forces_direct_newton_third() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(
                sum.abs() < 1e-10,
                "Newton III violated on axis {a}: sum = {sum}"
            );
        }
    }

    #[test]
    fn test_coulomb_forces_direct_attractive_direction() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        // Atom 0 attracted to atom 1 (+x)
        assert!(
            forces[0][0] > 0.0,
            "atom 0 should be attracted in +x, got {}",
            forces[0][0]
        );
        // Atom 1 attracted to atom 0 (-x)
        assert!(
            forces[1][0] < 0.0,
            "atom 1 should be attracted in -x, got {}",
            forces[1][0]
        );
    }

    #[test]
    fn test_coulomb_forces_direct_repulsive_direction() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        // Atom 0 repelled away from atom 1 (-x)
        assert!(
            forces[0][0] < 0.0,
            "atom 0 should be repelled in -x, got {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_coulomb_direct_vs_function() {
        // coulomb_energy_direct for a pair at r should match coulomb_energy(1,1,r)
        let r = 4.0;
        let positions = [[0.0, 0.0, 0.0], [r, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let e_direct = coulomb_energy_direct(&positions, &charges);
        let e_fn = coulomb_energy(1.0, 1.0, r);
        assert!(
            (e_direct - e_fn).abs() < 1e-10,
            "coulomb_energy_direct: {e_direct}, coulomb_energy: {e_fn}"
        );
    }

    #[test]
    fn test_coulomb_energy_scales_with_r() {
        let e1 = coulomb_energy(1.0, 1.0, 2.0);
        let e2 = coulomb_energy(1.0, 1.0, 4.0);
        assert!(
            (e1 - 2.0 * e2).abs() < 1e-10,
            "Coulomb energy scales as 1/r"
        );
    }

    #[test]
    fn test_coulomb_force_scales_with_r() {
        let f1 = coulomb_force_scalar(1.0, 1.0, 2.0);
        let f2 = coulomb_force_scalar(1.0, 1.0, 4.0);
        assert!(
            (f1 - 4.0 * f2).abs() < 1e-10,
            "Coulomb force scales as 1/r^2"
        );
    }

    #[test]
    fn test_coulomb_force_positive_same_sign() {
        let f = coulomb_force_scalar(1.0, 1.0, 3.0);
        assert!(
            f > 0.0,
            "same-sign repulsive force should be positive, got {f}"
        );
    }

    #[test]
    fn test_coulomb_force_negative_opposite_sign() {
        let f = coulomb_force_scalar(1.0, -1.0, 3.0);
        assert!(
            f < 0.0,
            "opposite-sign attractive force should be negative, got {f}"
        );
    }

    #[test]
    fn test_coulomb_energy_direct_single_particle_zero() {
        let positions = [[0.0, 0.0, 0.0]];
        let charges = [1.0];
        let e = coulomb_energy_direct(&positions, &charges);
        assert_eq!(e, 0.0, "single particle has no pair interaction");
    }

    #[test]
    fn test_coulomb_forces_direct_newton_third_v2() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        for (a, (&f0a, &f1a)) in forces[0].iter().zip(forces[1].iter()).enumerate() {
            let sum = f0a + f1a;
            assert!(sum.abs() < 1e-10, "Newton III violated on axis {a}: {sum}");
        }
    }

    #[test]
    fn test_coulomb_forces_direct_repulsion_direction() {
        // Charge +1 at x=0, charge +1 at x=3 -> force on particle 0 in -x direction
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        assert!(
            forces[0][0] < 0.0,
            "repulsive force on particle 0 should be -x, got {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_coulomb_forces_direct_attraction_direction() {
        // Charge +1 at x=0, charge -1 at x=3 -> force on particle 0 in +x direction
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let forces = coulomb_forces_direct(&positions, &charges);
        assert!(
            forces[0][0] > 0.0,
            "attractive force on particle 0 should be +x, got {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_coulomb_forces_direct_three_charges_sum_finite() {
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]];
        let charges = [1.0, -1.0, 0.5];
        let forces = coulomb_forces_direct(&positions, &charges);
        assert_eq!(forces.len(), 3);
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "force component must be finite");
            }
        }
    }

    #[test]
    fn test_dipole_dipole_energy_collinear_attractive() {
        // Two parallel dipoles along x, separation along x:
        // head-to-tail → attractive (negative energy)
        let d1 = [1.0, 0.0, 0.0];
        let d2 = [1.0, 0.0, 0.0];
        let r_vec = [3.0, 0.0, 0.0]; // along x
        let e = dipole_dipole_energy(d1, d2, r_vec);
        assert!(e < 0.0, "collinear head-to-tail dipoles → attractive: {e}");
    }

    #[test]
    fn test_dipole_dipole_energy_antiparallel_repulsive() {
        // Anti-parallel dipoles side by side → repulsive
        let d1 = [1.0, 0.0, 0.0];
        let d2 = [-1.0, 0.0, 0.0];
        let r_vec = [3.0, 0.0, 0.0];
        let e = dipole_dipole_energy(d1, d2, r_vec);
        assert!(e > 0.0, "antiparallel collinear dipoles → repulsive: {e}");
    }

    #[test]
    fn test_dipole_dipole_energy_zero_separation_zero() {
        let d1 = [1.0, 0.0, 0.0];
        let d2 = [1.0, 0.0, 0.0];
        let r_vec = [0.0; 3];
        let e = dipole_dipole_energy(d1, d2, r_vec);
        assert_eq!(e, 0.0, "zero separation → 0 (guard)");
    }

    #[test]
    fn test_dipole_dipole_force_finite() {
        let d1 = [1.0, 0.0, 0.0];
        let d2 = [0.0, 1.0, 0.0];
        let r_vec = [3.0, 2.0, 1.0];
        let f = dipole_dipole_force(d1, d2, r_vec);
        for &fc in &f {
            assert!(fc.is_finite(), "force component must be finite: {fc}");
        }
    }

    #[test]
    fn test_electric_field_single_charge_direction() {
        let src = [[0.0; 3]];
        let q = [1.0];
        let eval = vec![[1.0, 0.0, 0.0]];
        let e_field = electric_field(&src, &q, &eval);
        // E_x should be positive (pointing away from positive charge)
        assert!(
            e_field[0][0] > 0.0,
            "field from +q should point away: {}",
            e_field[0][0]
        );
        assert!(e_field[0][1].abs() < 1e-12);
        assert!(e_field[0][2].abs() < 1e-12);
    }

    #[test]
    fn test_electric_field_falls_off_as_r_squared() {
        let src = [[0.0; 3]];
        let q = [1.0];
        let eval1 = vec![[1.0, 0.0, 0.0]];
        let eval2 = vec![[2.0, 0.0, 0.0]];
        let e1 = electric_field(&src, &q, &eval1)[0][0];
        let e2 = electric_field(&src, &q, &eval2)[0][0];
        let ratio = e1 / e2;
        assert!((ratio - 4.0).abs() < 1e-10, "field ∝ 1/r², ratio = {ratio}");
    }

    #[test]
    fn test_electric_potential_single_charge() {
        let src = [[0.0; 3]];
        let q = [1.0];
        let eval = vec![[5.0, 0.0, 0.0]];
        let v = electric_potential(&src, &q, &eval)[0];
        assert!((v - COULOMB_K / 5.0).abs() < 1e-10, "V = {v}");
    }

    #[test]
    fn test_electric_potential_superposition() {
        let src = [[0.0; 3], [10.0, 0.0, 0.0]];
        let q = [1.0, 1.0];
        let eval = vec![[5.0, 0.0, 0.0]]; // midpoint
        let v = electric_potential(&src, &q, &eval)[0];
        let expected = 2.0 * COULOMB_K / 5.0;
        assert!(
            (v - expected).abs() < 1e-8,
            "V at midpoint = {v}, expected {expected}"
        );
    }

    #[test]
    fn test_screened_matrix_total_energy_negative_for_opposite_charges() {
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let mat = ScreenedInteractionMatrix::build(&positions, &charges, 0.1);
        let e = mat.total_energy();
        assert!(e < 0.0, "opposite charges → negative total energy: {e}");
    }

    #[test]
    fn test_screened_matrix_atom_energy_finite() {
        let positions = [[0.0; 3], [3.0, 0.0, 0.0], [6.0, 0.0, 0.0]];
        let charges = [1.0, -1.0, 0.5];
        let mat = ScreenedInteractionMatrix::build(&positions, &charges, 0.1);
        let e = mat.atom_energy(0);
        assert!(e.is_finite(), "atom energy must be finite: {e}");
    }

    #[test]
    fn test_screened_matrix_same_charge_repulsive() {
        let positions = [[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = [1.0, 1.0];
        let mat = ScreenedInteractionMatrix::build(&positions, &charges, 0.0);
        let e = mat.total_energy();
        assert!(e > 0.0, "same charges → positive energy: {e}");
    }

    #[test]
    fn test_screened_matrix_n_correct() {
        let positions = vec![[0.0; 3]; 5];
        let charges = vec![1.0; 5];
        let mat = ScreenedInteractionMatrix::build(&positions, &charges, 0.1);
        assert_eq!(mat.n, 5);
        assert_eq!(mat.w.len(), 5);
    }

    #[test]
    fn test_screened_matrix_larger_kappa_smaller_energy() {
        // Stronger screening → smaller interaction energy magnitude
        let positions = [[0.0; 3], [2.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e_low = ScreenedInteractionMatrix::build(&positions, &charges, 0.01)
            .total_energy()
            .abs();
        let e_high = ScreenedInteractionMatrix::build(&positions, &charges, 2.0)
            .total_energy()
            .abs();
        assert!(
            e_high < e_low,
            "stronger screening → weaker interaction: {e_high} < {e_low}"
        );
    }

    #[test]
    fn test_ewald_reciprocal_k_energy_finite() {
        let positions = [[0.0; 3], [2.5, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = ewald_reciprocal_k_energy(&positions, &charges, [1.0, 0.0, 0.0], 0.3, 125.0);
        assert!(e.is_finite(), "reciprocal energy = {e}");
    }

    #[test]
    fn test_ewald_reciprocal_k_zero_for_zero_k() {
        let positions = [[0.0; 3]];
        let charges = [1.0];
        let e = ewald_reciprocal_k_energy(&positions, &charges, [0.0; 3], 0.3, 125.0);
        assert_eq!(e, 0.0, "k=0 returns 0");
    }

    #[test]
    fn test_ewald_grid_energy_finite() {
        let positions = [[0.0; 3], [2.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = ewald_reciprocal_energy_grid(&positions, &charges, [10.0; 3], 0.3, 2);
        assert!(e.is_finite(), "grid energy = {e}");
    }

    #[test]
    fn test_dipole_in_field_energy_formula() {
        let mu = [1.0, 0.0, 0.0];
        let e = [2.0, 0.0, 0.0];
        let u = dipole_in_field_energy(mu, e);
        assert!((u + 2.0).abs() < 1e-12, "U = -mu.E = -2, got {u}");
    }

    #[test]
    fn test_dipole_torque_perpendicular() {
        let mu = [1.0, 0.0, 0.0];
        let e = [0.0, 1.0, 0.0];
        let tau = dipole_torque(mu, e);
        // mu × E = (0*0 - 0*1, 0*0 - 1*0, 1*1 - 0*0) = (0, 0, 1)
        assert!((tau[2] - 1.0).abs() < 1e-12, "tau_z = 1, got {}", tau[2]);
    }

    #[test]
    fn test_total_dipole_moment() {
        let pos = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [-1.0, 1.0];
        let mu = total_dipole_moment(&pos, &charges);
        // mu = -1*[0,0,0] + 1*[1,0,0] = [1,0,0]
        assert!((mu[0] - 1.0).abs() < 1e-12, "mu_x = 1, got {}", mu[0]);
    }

    #[test]
    fn test_lj_coulomb_pair_energy_finite() {
        let e = lj_coulomb_pair_energy(0.5, -0.5, 3.0, 1.0, 3.0);
        assert!(e.is_finite(), "LJ+Coulomb energy = {e}");
    }

    #[test]
    fn test_lj_coulomb_pair_energy_zero_at_zero_r() {
        let e = lj_coulomb_pair_energy(1.0, -1.0, 0.0, 1.0, 3.0);
        assert_eq!(e, 0.0, "zero r guard");
    }
}
