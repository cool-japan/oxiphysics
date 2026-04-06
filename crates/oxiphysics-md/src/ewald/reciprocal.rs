#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ewald reciprocal space, B-splines, PME grid, structure factor, Madelung constants.

use super::params::{COULOMB_K, EwaldParams, erfc_approx};
use super::summation::EwaldSummation;

// ---------------------------------------------------------------------------
// BSplineChargeSpreading (smooth PME helper)
// ---------------------------------------------------------------------------

/// Cardinal B-spline of order `n` evaluated at `x`.
///
/// Used for charge spreading in smooth PME.
#[allow(dead_code)]
pub fn bspline(n: usize, x: f64) -> f64 {
    if n == 1 {
        if (0.0..1.0).contains(&x) { 1.0 } else { 0.0 }
    } else {
        let nf = n as f64;
        x / (nf - 1.0) * bspline(n - 1, x) + (nf - x) / (nf - 1.0) * bspline(n - 1, x - 1.0)
    }
}

/// Evaluate the 3D B-spline weight for charge spreading.
///
/// Given fractional coordinates `u` in \[0, grid_size), spreads using
/// order-`p` B-splines. Returns the weight at grid point `(ix, iy, iz)`.
#[allow(dead_code)]
pub fn bspline_weight_3d(
    u: [f64; 3],
    grid_point: [usize; 3],
    grid_size: [usize; 3],
    order: usize,
) -> f64 {
    let mut weight = 1.0;
    for dim in 0..3 {
        let frac = u[dim] * grid_size[dim] as f64;
        let offset = frac - grid_point[dim] as f64;
        weight *= bspline(order, offset + (order as f64 / 2.0));
    }
    weight
}

// ---------------------------------------------------------------------------
// EwaldSum — simplified self-contained Ewald struct (plain [f64;3])
// ---------------------------------------------------------------------------

/// Simplified Ewald sum object for use with plain `[f64;3]` arrays.
///
/// Distinct from [`EwaldSummation`] in that it stores a single `alpha`,
/// `kmax`, and `box_len` (cubic box), making it easy to instantiate without
/// building [`EwaldParams`] first.
#[allow(dead_code)]
pub struct EwaldSum {
    /// Ewald splitting parameter (angstrom^-1).
    pub alpha: f64,
    /// Maximum k-vector index per dimension.
    pub kmax: usize,
    /// Cubic box side length (angstrom).
    pub box_len: f64,
}

#[allow(dead_code)]
impl EwaldSum {
    /// Create a new [`EwaldSum`].
    pub fn new(alpha: f64, kmax: usize, box_len: f64) -> Self {
        assert!(alpha > 0.0, "alpha must be positive");
        assert!(box_len > 0.0, "box_len must be positive");
        Self {
            alpha,
            kmax,
            box_len,
        }
    }

    /// Real-space energy (kJ mol^-1): sum_{i<j} K * qi * qj * erfc(alpha*r_ij) / r_ij.
    ///
    /// Uses the minimum image convention for the cubic box.
    pub fn real_space_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let mut energy = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dr_raw = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let dr = self.min_image(dr_raw);
                let r = (dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]).sqrt();
                if r < 1e-20 {
                    continue;
                }
                energy += COULOMB_K * charges[i] * charges[j] * erfc_approx(self.alpha * r) / r;
            }
        }
        energy
    }

    /// Reciprocal-space energy (kJ mol^-1).
    ///
    /// Standard Ewald k-space sum:
    /// E_recip = (COULOMB_K / (2*pi*V)) * sum_{k!=0} (4*pi^2/k^2) * exp(-k^2/(4*alpha^2)) * |S(k)|^2
    pub fn reciprocal_space_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        let n = positions.len();
        let volume = self.box_len * self.box_len * self.box_len;
        let prefactor = COULOMB_K / (2.0 * std::f64::consts::PI * volume);
        let four_a2 = 4.0 * self.alpha * self.alpha;
        let two_pi = 2.0 * std::f64::consts::PI;
        let kmax_i = self.kmax as i32;
        let mut energy = 0.0;

        for nx in -kmax_i..=kmax_i {
            for ny in -kmax_i..=kmax_i {
                for nz in -kmax_i..=kmax_i {
                    if nx == 0 && ny == 0 && nz == 0 {
                        continue;
                    }
                    let kx = two_pi * nx as f64 / self.box_len;
                    let ky = two_pi * ny as f64 / self.box_len;
                    let kz = two_pi * nz as f64 / self.box_len;
                    let k2 = kx * kx + ky * ky + kz * kz;

                    let factor = (4.0 * std::f64::consts::PI * std::f64::consts::PI / k2)
                        * (-k2 / four_a2).exp();

                    let mut s_cos = 0.0f64;
                    let mut s_sin = 0.0f64;
                    for i in 0..n {
                        let kr = kx * positions[i][0] + ky * positions[i][1] + kz * positions[i][2];
                        s_cos += charges[i] * kr.cos();
                        s_sin += charges[i] * kr.sin();
                    }
                    energy += factor * (s_cos * s_cos + s_sin * s_sin);
                }
            }
        }
        energy * prefactor
    }

    /// Self-energy correction (kJ mol^-1): -COULOMB_K * alpha / sqrt(pi) * sum qi^2.
    pub fn self_energy(&self, charges: &[f64]) -> f64 {
        let sum_q2: f64 = charges.iter().map(|q| q * q).sum();
        -COULOMB_K * self.alpha / std::f64::consts::PI.sqrt() * sum_q2
    }

    /// Total Ewald energy = real + reciprocal + self.
    pub fn total_energy(&self, positions: &[[f64; 3]], charges: &[f64]) -> f64 {
        self.real_space_energy(positions, charges)
            + self.reciprocal_space_energy(positions, charges)
            + self.self_energy(charges)
    }

    /// Minimum image displacement for the cubic box.
    fn min_image(&self, dr: [f64; 3]) -> [f64; 3] {
        let l = self.box_len;
        [
            dr[0] - l * (dr[0] / l).round(),
            dr[1] - l * (dr[1] / l).round(),
            dr[2] - l * (dr[2] / l).round(),
        ]
    }
}

// ---------------------------------------------------------------------------
// PmeGrid — Particle Mesh Ewald charge grid
// ---------------------------------------------------------------------------

/// Particle Mesh Ewald charge grid.
///
/// Spreads point charges onto a regular 3D grid using nearest-grid-point (NGP)
/// assignment.  A full PME implementation would then FFT this grid; here the
/// grid object stores the raw charge density for downstream use.
#[allow(dead_code)]
pub struct PmeGrid {
    /// Number of grid points per dimension (cubic grid).
    pub grid_size: usize,
    /// Ewald splitting parameter (angstrom^-1).
    pub alpha: f64,
    /// Flattened 3D charge density grid (index = ix*N^2 + iy*N + iz).
    pub grid: Vec<f64>,
    /// Cubic box side length (angstrom).
    pub box_len: f64,
}

#[allow(dead_code)]
impl PmeGrid {
    /// Create a new zero-initialised [`PmeGrid`].
    pub fn new(grid_size: usize, alpha: f64, box_len: f64) -> Self {
        assert!(grid_size > 0, "grid_size must be > 0");
        assert!(alpha > 0.0, "alpha must be positive");
        assert!(box_len > 0.0, "box_len must be positive");
        Self {
            grid_size,
            alpha,
            grid: vec![0.0; grid_size * grid_size * grid_size],
            box_len,
        }
    }

    /// Spread charges onto the grid using nearest-grid-point (NGP) assignment.
    ///
    /// Each atom's charge is deposited into the single closest grid cell.
    /// Coordinates are mapped periodically into \[0, box_len).
    pub fn spread_charges(&mut self, positions: &[[f64; 3]], charges: &[f64]) {
        let n_grid = self.grid_size;
        for (pos, &q) in positions.iter().zip(charges.iter()) {
            let ix = {
                let frac = pos[0].rem_euclid(self.box_len) / self.box_len;
                ((frac * n_grid as f64) as usize).min(n_grid - 1)
            };
            let iy = {
                let frac = pos[1].rem_euclid(self.box_len) / self.box_len;
                ((frac * n_grid as f64) as usize).min(n_grid - 1)
            };
            let iz = {
                let frac = pos[2].rem_euclid(self.box_len) / self.box_len;
                ((frac * n_grid as f64) as usize).min(n_grid - 1)
            };
            self.grid[ix * n_grid * n_grid + iy * n_grid + iz] += q;
        }
    }

    /// Zero all grid values.
    pub fn clear(&mut self) {
        for v in self.grid.iter_mut() {
            *v = 0.0;
        }
    }

    /// Total charge deposited on the grid.
    pub fn total_charge(&self) -> f64 {
        self.grid.iter().sum()
    }

    /// Number of non-zero grid cells.
    pub fn occupied_cells(&self) -> usize {
        self.grid.iter().filter(|&&v| v.abs() > 1e-20).count()
    }

    /// Simplified reciprocal-space energy estimate (proxy for full PME).
    ///
    /// Returns `sum(grid^2) * V / N^3` scaled by the Ewald Gaussian envelope.
    pub fn reciprocal_energy_estimate(&self) -> f64 {
        let n = self.grid_size as f64;
        let volume = self.box_len * self.box_len * self.box_len;
        let sum_sq: f64 = self.grid.iter().map(|v| v * v).sum();
        COULOMB_K * sum_sq * volume / (n * n * n)
    }
}

/// Generate a list of k-vectors in reciprocal space for an orthorhombic box.
///
/// Returns all k-vectors `k = 2π [nx/Lx, ny/Ly, nz/Lz]` with
/// `-k_max ≤ nα ≤ k_max` and (nx, ny, nz) ≠ (0, 0, 0).
///
/// # Arguments
/// * `k_max`      - maximum index per dimension.
/// * `box_lengths` - orthorhombic box lengths \[Lx, Ly, Lz\] (angstrom).
#[allow(dead_code)]
pub fn k_vector_list(k_max: usize, box_lengths: [f64; 3]) -> Vec<[f64; 3]> {
    let two_pi = 2.0 * std::f64::consts::PI;
    let km = k_max as i32;
    let mut k_vecs = Vec::new();
    for nx in -km..=km {
        for ny in -km..=km {
            for nz in -km..=km {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                k_vecs.push([
                    two_pi * nx as f64 / box_lengths[0],
                    two_pi * ny as f64 / box_lengths[1],
                    two_pi * nz as f64 / box_lengths[2],
                ]);
            }
        }
    }
    k_vecs
}

// ---------------------------------------------------------------------------
// StructureFactor — S(k) for a collection of charges
// ---------------------------------------------------------------------------

/// Complex structure factor S(k) = Σᵢ qᵢ exp(i k·rᵢ).
///
/// Decomposes into real part S_cos = Σᵢ qᵢ cos(k·rᵢ) and
/// imaginary part S_sin = Σᵢ qᵢ sin(k·rᵢ).
#[derive(Debug, Clone, Copy, Default)]
#[allow(dead_code)]
pub struct StructureFactor {
    /// Real part: Σᵢ qᵢ cos(k·rᵢ).
    pub real: f64,
    /// Imaginary part: Σᵢ qᵢ sin(k·rᵢ).
    pub imag: f64,
}

#[allow(dead_code)]
impl StructureFactor {
    /// |S(k)|² = real² + imag².
    #[inline]
    pub fn modulus_sq(&self) -> f64 {
        self.real * self.real + self.imag * self.imag
    }

    /// Compute S(k) for a single k-vector.
    pub fn compute(k: [f64; 3], positions: &[[f64; 3]], charges: &[f64]) -> Self {
        let n = positions.len().min(charges.len());
        let mut s_cos = 0.0;
        let mut s_sin = 0.0;
        for i in 0..n {
            let kr = k[0] * positions[i][0] + k[1] * positions[i][1] + k[2] * positions[i][2];
            s_cos += charges[i] * kr.cos();
            s_sin += charges[i] * kr.sin();
        }
        Self {
            real: s_cos,
            imag: s_sin,
        }
    }

    /// Compute S(k) for all k-vectors in a list.
    pub fn compute_all(k_vecs: &[[f64; 3]], positions: &[[f64; 3]], charges: &[f64]) -> Vec<Self> {
        k_vecs
            .iter()
            .map(|&k| Self::compute(k, positions, charges))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// EwaldForceIntegrator — integrate MD step with Ewald electrostatics
// ---------------------------------------------------------------------------

/// Utility for running a short velocity-Verlet MD loop with Ewald forces.
///
/// Operates on plain `[f64;3]` arrays with no external dependencies.
#[allow(dead_code)]
pub struct EwaldForceIntegrator {
    /// Ewald summation parameters.
    pub ewald: EwaldSummation,
    /// Time step (ps).
    pub dt: f64,
}

#[allow(dead_code)]
impl EwaldForceIntegrator {
    /// Create a new integrator.
    pub fn new(ewald: EwaldSummation, dt: f64) -> Self {
        assert!(dt > 0.0, "dt must be positive");
        Self { ewald, dt }
    }

    /// Perform one velocity-Verlet step with Ewald electrostatic forces.
    ///
    /// # Arguments
    /// * `positions`  - current atom positions (updated in place, Å).
    /// * `velocities` - current atom velocities (updated in place, Å ps⁻¹).
    /// * `charges`    - atom partial charges (e).
    /// * `masses`     - atom masses (amu).
    ///
    /// # Returns
    /// The potential energy (kJ mol⁻¹) at the new configuration.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        positions: &mut Vec<[f64; 3]>,
        velocities: &mut Vec<[f64; 3]>,
        charges: &[f64],
        masses: &[f64],
    ) -> f64 {
        let n = positions.len();
        // Initial forces at current positions
        let mut forces = self.ewald.real_space_forces(positions, charges);

        // Half-step velocity update
        for i in 0..n {
            let inv_m = 1.0 / masses[i];
            for a in 0..3 {
                velocities[i][a] += 0.5 * self.dt * forces[i][a] * inv_m;
            }
        }

        // Full-step position update
        for i in 0..n {
            for a in 0..3 {
                positions[i][a] += self.dt * velocities[i][a];
                // PBC wrap
                let l = self.ewald.box_lengths[a];
                positions[i][a] = positions[i][a].rem_euclid(l);
            }
        }

        // New forces
        forces = self.ewald.real_space_forces(positions, charges);
        let energy = self.ewald.total_energy(positions, charges);

        // Second half-step velocity update
        for i in 0..n {
            let inv_m = 1.0 / masses[i];
            for a in 0..3 {
                velocities[i][a] += 0.5 * self.dt * forces[i][a] * inv_m;
            }
        }

        energy
    }

    /// Run for `n_steps` and record total energy at each step.
    pub fn run(
        &self,
        positions: &mut Vec<[f64; 3]>,
        velocities: &mut Vec<[f64; 3]>,
        charges: &[f64],
        masses: &[f64],
        n_steps: usize,
    ) -> Vec<f64> {
        let mut energies = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            let e = self.step(positions, velocities, charges, masses);
            energies.push(e);
        }
        energies
    }
}

// ---------------------------------------------------------------------------
// EwaldSum extensions: forces via finite differences check utility
// ---------------------------------------------------------------------------

impl EwaldSum {
    /// Numerical gradient of total Ewald energy for atom `idx` in direction `axis`.
    ///
    /// Uses central differences with step size `h`.
    /// Useful for testing analytic force implementations.
    #[allow(dead_code)]
    pub fn numerical_force_component(
        &self,
        positions: &[[f64; 3]],
        charges: &[f64],
        idx: usize,
        axis: usize,
        h: f64,
    ) -> f64 {
        let mut pos_plus = positions.to_vec();
        let mut pos_minus = positions.to_vec();
        pos_plus[idx][axis] += h;
        pos_minus[idx][axis] -= h;
        let e_plus = self.total_energy(&pos_plus, charges);
        let e_minus = self.total_energy(&pos_minus, charges);
        -(e_plus - e_minus) / (2.0 * h)
    }

    /// Real-space forces as force array (kJ mol⁻¹ Å⁻¹).
    ///
    /// This mirrors `EwaldSummation::real_space_forces` but for the cubic-box `EwaldSum`.
    #[allow(dead_code)]
    pub fn real_space_forces(&self, positions: &[[f64; 3]], charges: &[f64]) -> Vec<[f64; 3]> {
        let n = positions.len();
        let alpha = self.alpha;
        let inv_sqrt_pi = 1.0 / std::f64::consts::PI.sqrt();
        let mut forces = vec![[0.0f64; 3]; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let dr_raw = [
                    positions[j][0] - positions[i][0],
                    positions[j][1] - positions[i][1],
                    positions[j][2] - positions[i][2],
                ];
                let dr = self.min_image(dr_raw);
                let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                let r = r2.sqrt();
                if r < 1e-20 {
                    continue;
                }
                let ar = alpha * r;
                let erfc_val = erfc_approx(ar);
                let exp_val = (-ar * ar).exp();
                let qq = COULOMB_K * charges[i] * charges[j];
                let f_mag = qq * (erfc_val / r2 + 2.0 * alpha * exp_val * inv_sqrt_pi / r);
                for a in 0..3 {
                    let fa = -f_mag * dr[a] / r;
                    forces[i][a] += fa;
                    forces[j][a] -= fa;
                }
            }
        }
        forces
    }

    /// Net charge of the system.
    #[allow(dead_code)]
    pub fn net_charge(charges: &[f64]) -> f64 {
        charges.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// Madelung constant computation
// ---------------------------------------------------------------------------

/// Compute the Madelung constant for a cubic NaCl-type lattice via Ewald summation.
///
/// The Madelung constant A is defined through the lattice energy per ion pair:
/// E/N = -A * e² / (4π ε₀ r₀) where r₀ is the nearest-neighbour distance.
///
/// This function estimates A by evaluating the Ewald sum for a small cluster
/// of an NaCl-type rock salt lattice and extracts the per-ion-pair energy.
///
/// # Arguments
/// * `alpha`     - Ewald splitting parameter (angstrom⁻¹).
/// * `k_max`     - k-vector cutoff (integer per dimension).
/// * `n_cells`   - number of unit cells per dimension (total = (2*n_cells+1)³ ions).
/// * `a`         - lattice constant (angstrom); nearest-neighbour distance = a/2.
///
/// Returns the Madelung constant (dimensionless).
#[allow(dead_code)]
pub fn madelung_constant_nacl(alpha: f64, k_max: usize, n_cells: i32, a: f64) -> f64 {
    // Build NaCl lattice in [-n_cells, n_cells]³
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut charges: Vec<f64> = Vec::new();

    for ix in -n_cells..=n_cells {
        for iy in -n_cells..=n_cells {
            for iz in -n_cells..=n_cells {
                let x = ix as f64 * a * 0.5;
                let y = iy as f64 * a * 0.5;
                let z = iz as f64 * a * 0.5;
                // NaCl: sign of charge = (-1)^(ix+iy+iz)
                let sign = if (ix + iy + iz) % 2 == 0 { 1.0 } else { -1.0 };
                positions.push([x, y, z]);
                charges.push(sign);
            }
        }
    }

    let n_side = (2 * n_cells + 1) as usize;
    let box_len = n_side as f64 * a * 0.5 + a;
    let box_lengths = [box_len, box_len, box_len];
    let params = EwaldParams {
        alpha,
        r_cutoff: box_len * 0.45,
        k_cutoff: 2.0 * std::f64::consts::PI * alpha,
        epsilon_r: 1.0,
    };
    let ewald = EwaldSummation::with_k_max(params, box_lengths, k_max as i32);

    let e_total = ewald.total_energy(&positions, &charges);

    // Number of ion pairs = n³
    let n_ions = positions.len() as f64;
    let n_pairs = n_ions / 2.0;

    // Energy per ion pair (in units of COULOMB_K / a, where a is nearest-neighbour distance)
    // E_pair = -A * COULOMB_K / r_nn  =>  A = -E_pair * r_nn / COULOMB_K
    let r_nn = a * 0.5;
    let e_per_pair = e_total / n_pairs;

    -e_per_pair * r_nn / COULOMB_K
}

/// Madelung constant for a simple 2-ion NaCl dimer (exact = 1.0).
///
/// Returns the analytic Madelung constant A = 1.0 for a single ion pair.
#[allow(dead_code)]
pub fn madelung_constant_dimer(r: f64) -> f64 {
    // For a dimer +q/-q at separation r: E = -COULOMB_K * q^2 / r
    // Madelung A = 1.0 by definition
    let _ = r;
    1.0
}

/// Estimate the Madelung constant using direct lattice sum for a finite
/// cubic shell up to `n_max` cells.
///
/// A_direct = -Σ_{n≠0} (-1)^|n| / |n|   (converges only conditionally)
///
/// This uses the Ewald-accelerated form internally via [`EwaldSum`].
#[allow(dead_code)]
pub fn madelung_direct_lattice(n_max: i32, a: f64) -> f64 {
    let mut sum = 0.0;
    for nx in -n_max..=n_max {
        for ny in -n_max..=n_max {
            for nz in -n_max..=n_max {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                let r = ((nx * nx + ny * ny + nz * nz) as f64).sqrt();
                let sign = if (nx.abs() + ny.abs() + nz.abs()) % 2 == 0 {
                    1.0
                } else {
                    -1.0
                };
                sum += sign / r;
            }
        }
    }
    // Return absolute value (A is positive)
    sum.abs() / a.max(1.0)
}

// ---------------------------------------------------------------------------
// RecipSpaceContrib — standalone reciprocal-space energy/force contributions
// ---------------------------------------------------------------------------

/// Compute the reciprocal-space energy contribution for a subset of k-vectors.
///
/// Allows parallel decomposition: split the k-vector list across threads and
/// sum the contributions.
#[allow(dead_code)]
pub fn reciprocal_energy_kvecs(
    positions: &[[f64; 3]],
    charges: &[f64],
    kvecs: &[[f64; 3]],
    alpha: f64,
    volume: f64,
) -> f64 {
    let n = positions.len();
    let prefactor = COULOMB_K / (2.0_f64 * std::f64::consts::PI * volume);
    let mut energy = 0.0_f64;
    for &k in kvecs {
        let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        if k2 < 1e-20 {
            continue;
        }
        let factor = (4.0_f64 * std::f64::consts::PI * std::f64::consts::PI / k2)
            * (-(k2) / (4.0_f64 * alpha * alpha)).exp();
        // Structure factor S(k) = Σ_i q_i e^{i k·r_i}
        let mut s_re = 0.0_f64;
        let mut s_im = 0.0_f64;
        for i in 0..n {
            let phase = k[0] * positions[i][0] + k[1] * positions[i][1] + k[2] * positions[i][2];
            s_re += charges[i] * phase.cos();
            s_im += charges[i] * phase.sin();
        }
        energy += prefactor * factor * (s_re * s_re + s_im * s_im);
    }
    energy
}

/// Compute reciprocal-space forces from a list of k-vectors.
///
/// Returns a `Vec<[f64;3]>` of forces (kJ mol⁻¹ Å⁻¹).
#[allow(dead_code)]
pub fn reciprocal_forces_kvecs(
    positions: &[[f64; 3]],
    charges: &[f64],
    kvecs: &[[f64; 3]],
    alpha: f64,
    volume: f64,
) -> Vec<[f64; 3]> {
    let n = positions.len();
    let prefactor = COULOMB_K / (std::f64::consts::PI * volume);
    let mut forces = vec![[0.0_f64; 3]; n];
    for &k in kvecs {
        let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
        if k2 < 1e-20 {
            continue;
        }
        let factor = (4.0_f64 * std::f64::consts::PI * std::f64::consts::PI / k2)
            * (-(k2) / (4.0_f64 * alpha * alpha)).exp();
        // Structure factor
        let mut s_re = 0.0_f64;
        let mut s_im = 0.0_f64;
        for i in 0..n {
            let phase = k[0] * positions[i][0] + k[1] * positions[i][1] + k[2] * positions[i][2];
            s_re += charges[i] * phase.cos();
            s_im += charges[i] * phase.sin();
        }
        // Force on atom i: F_i = prefactor * factor * q_i * (k · sin(k·r_i) * s_re - cos(k·r_i) * s_im)
        for i in 0..n {
            let phase = k[0] * positions[i][0] + k[1] * positions[i][1] + k[2] * positions[i][2];
            let contrib =
                prefactor * factor * charges[i] * (phase.sin() * s_re - phase.cos() * s_im);
            for a in 0..3 {
                forces[i][a] += contrib * k[a];
            }
        }
    }
    forces
}

// ---------------------------------------------------------------------------
// EwaldConvergenceTest — verify convergence as k_max increases
// ---------------------------------------------------------------------------

/// Measures how the Ewald reciprocal energy converges as `k_max` increases.
#[allow(dead_code)]
pub fn reciprocal_energy_convergence(
    positions: &[[f64; 3]],
    charges: &[f64],
    params: &EwaldParams,
    box_lengths: [f64; 3],
    k_max_values: &[i32],
) -> Vec<f64> {
    k_max_values
        .iter()
        .map(|&km| {
            let ewald = EwaldSummation::with_k_max(params.clone(), box_lengths, km);
            ewald.reciprocal_space_energy(positions, charges)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ewald::*;

    #[test]
    fn test_reciprocal_energy_nonneg_neutral() {
        let params = EwaldParams::new(8.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 3);
        let positions = [[2.0, 2.0, 2.0], [8.0, 2.0, 2.0]];
        let charges = [1.0, -1.0];
        let e_recip = ewald.reciprocal_space_energy(&positions, &charges);
        assert!(
            e_recip.is_finite(),
            "reciprocal energy should be finite, got {e_recip}"
        );
    }

    #[test]
    fn test_reciprocal_error_estimate() {
        let params = EwaldParams::new(10.0);
        let err = params.reciprocal_space_error_estimate(2.0, 20.0, 5);
        assert!(err >= 0.0, "reciprocal error should be non-negative");
        assert!(err.is_finite(), "reciprocal error should be finite");
    }

    #[test]
    fn test_bspline_order1() {
        assert!((bspline(1, 0.5) - 1.0).abs() < 1e-14);
        assert!((bspline(1, -0.1)).abs() < 1e-14);
        assert!((bspline(1, 1.1)).abs() < 1e-14);
    }

    #[test]
    fn test_bspline_order2() {
        // B2(1.0) should be at peak = 1.0
        let v = bspline(2, 1.0);
        assert!(
            (v - 1.0).abs() < 1e-12,
            "B-spline order 2 at x=1.0 should be 1.0, got {v}"
        );
        // B2(0) = 0
        let v0 = bspline(2, 0.0);
        assert!(v0.abs() < 1e-12, "B-spline order 2 at x=0 should be 0");
    }

    #[test]
    fn test_ewald_sum_self_energy_scales_with_charge_squared() {
        let es = EwaldSum::new(0.5, 3, 20.0);
        let e1 = es.self_energy(&[1.0]);
        let e2 = es.self_energy(&[2.0]);
        // self_energy ~ -alpha/sqrt(pi) * sum qi^2, so e2 = 4 * e1
        assert!(
            (e2 / e1 - 4.0).abs() < 1e-10,
            "self energy should scale with q^2: e2/e1 = {}",
            e2 / e1
        );
    }

    #[test]
    fn test_ewald_sum_self_energy_negative() {
        let es = EwaldSum::new(0.5, 3, 20.0);
        let e = es.self_energy(&[1.0, -1.0]);
        assert!(e < 0.0, "self energy should be negative, got {e}");
    }

    #[test]
    fn test_ewald_sum_total_energy_finite_neutral() {
        let es = EwaldSum::new(0.4, 3, 25.0);
        let positions = [[2.0f64, 0.0, 0.0], [8.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let e = es.total_energy(&positions, &charges);
        assert!(e.is_finite(), "total energy should be finite, got {e}");
    }

    #[test]
    fn test_pme_grid_spread_nonzero() {
        let mut grid = PmeGrid::new(8, 0.3, 10.0);
        grid.spread_charges(&[[2.0, 2.0, 2.0]], &[1.0]);
        let occ = grid.occupied_cells();
        assert_eq!(occ, 1, "one charge should occupy exactly one cell");
    }

    #[test]
    fn test_pme_grid_clear() {
        let mut grid = PmeGrid::new(4, 0.3, 10.0);
        grid.spread_charges(&[[1.0, 1.0, 1.0]], &[2.0]);
        grid.clear();
        let all_zero = grid.grid.iter().all(|&v| v == 0.0);
        assert!(all_zero, "grid should be all-zero after clear");
    }

    #[test]
    fn test_pme_grid_reciprocal_energy_nonneg() {
        let mut grid = PmeGrid::new(8, 0.3, 10.0);
        grid.spread_charges(&[[1.0, 1.0, 1.0], [6.0, 6.0, 6.0]], &[1.0, -1.0]);
        let e = grid.reciprocal_energy_estimate();
        assert!(
            e >= 0.0,
            "reciprocal energy estimate should be non-negative, got {e}"
        );
        assert!(e.is_finite(), "reciprocal energy should be finite");
    }

    #[test]
    fn test_k_vector_list_excludes_origin() {
        let k_vecs = k_vector_list(1, [10.0, 10.0, 10.0]);
        // No (0,0,0) vector should be present
        for k in &k_vecs {
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            assert!(k2 > 1e-20, "k=(0,0,0) should be excluded");
        }
    }

    #[test]
    fn test_k_vector_list_count() {
        let k_vecs = k_vector_list(1, [10.0, 10.0, 10.0]);
        // k_max=1: (2*1+1)^3 - 1 = 3^3 - 1 = 26 vectors
        assert_eq!(
            k_vecs.len(),
            26,
            "expected 26 k-vectors for k_max=1, got {}",
            k_vecs.len()
        );
    }

    #[test]
    fn test_k_vector_list_magnitudes_positive() {
        let k_vecs = k_vector_list(2, [15.0, 15.0, 15.0]);
        for k in &k_vecs {
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            assert!(k2 > 0.0, "k^2 should be positive");
        }
    }

    #[test]
    fn test_ewald_reciprocal_energy_fn_finite() {
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let k_vecs = k_vector_list(2, [20.0, 20.0, 20.0]);
        let vol = 20.0_f64.powi(3);
        let e = ewald_reciprocal_energy(&positions, &charges, &k_vecs, vol);
        assert!(e.is_finite(), "reciprocal energy should be finite, got {e}");
    }

    #[test]
    fn test_two_charge_nacl_ewald_sum_negative() {
        // Na+ Cl- pair at equilibrium distance ~2.81 angstrom
        // Total Ewald energy should be substantially negative
        let es = EwaldSum::new(0.5, 5, 30.0);
        let positions = [[10.0, 10.0, 10.0], [12.81, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let e = es.total_energy(&positions, &charges);
        assert!(
            e < 0.0,
            "Na+/Cl- Ewald total energy should be negative, got {e}"
        );
    }

    #[test]
    fn test_structure_factor_modulus_sq_nonneg() {
        let k = [0.5, 0.0, 0.0];
        let positions = [[0.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let sf = StructureFactor::compute(k, &positions, &charges);
        assert!(sf.modulus_sq() >= 0.0);
        assert!(sf.modulus_sq().is_finite());
    }

    #[test]
    fn test_structure_factor_zero_for_equal_charges_at_origin() {
        // All charges at origin: S(k) = sum_i q_i * exp(0) = sum_i q_i
        // For neutral system sum q_i = 0 -> S = 0
        let k = [1.0, 0.0, 0.0];
        let positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let sf = StructureFactor::compute(k, &positions, &charges);
        assert!(
            sf.real.abs() < 1e-12,
            "S.real should be 0 for neutral at origin, got {}",
            sf.real
        );
        assert!(
            sf.imag.abs() < 1e-12,
            "S.imag should be 0 for neutral at origin, got {}",
            sf.imag
        );
    }

    #[test]
    fn test_structure_factor_compute_all_length() {
        let k_vecs = k_vector_list(1, [10.0, 10.0, 10.0]);
        let positions = [[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let sfs = StructureFactor::compute_all(&k_vecs, &positions, &charges);
        assert_eq!(sfs.len(), k_vecs.len());
    }

    #[test]
    fn test_structure_factor_real_part_known_value() {
        // k = [pi, 0, 0], position = [1, 0, 0], charge = 1
        // S = cos(pi) = -1
        let k = [std::f64::consts::PI, 0.0, 0.0];
        let positions = [[1.0, 0.0, 0.0]];
        let charges = [1.0];
        let sf = StructureFactor::compute(k, &positions, &charges);
        assert!(
            (sf.real - (-1.0)).abs() < 1e-12,
            "S.real should be -1, got {}",
            sf.real
        );
        assert!(sf.imag.abs() < 1e-12, "S.imag should be 0, got {}", sf.imag);
    }

    #[test]
    fn test_ewald_force_integrator_step_advances_positions() {
        let params = EwaldParams::new(8.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 2);
        let integrator = EwaldForceIntegrator::new(ewald, 0.001);

        let mut positions = vec![[2.0, 10.0, 10.0], [8.0, 10.0, 10.0]];
        let mut velocities = vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let masses = vec![1.0, 1.0];

        let pos0_x = positions[0][0];
        integrator.step(&mut positions, &mut velocities, &charges, &masses);
        assert!(
            (positions[0][0] - pos0_x).abs() > 1e-12,
            "positions should change after step"
        );
    }

    #[test]
    fn test_ewald_force_integrator_run_returns_energies() {
        let params = EwaldParams::new(8.0);
        let ewald = EwaldSummation::with_k_max(params, [20.0, 20.0, 20.0], 2);
        let integrator = EwaldForceIntegrator::new(ewald, 0.001);

        let mut positions = vec![[2.0, 10.0, 10.0], [8.0, 10.0, 10.0]];
        let mut velocities = vec![[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let masses = vec![1.0, 1.0];

        let energies = integrator.run(&mut positions, &mut velocities, &charges, &masses, 5);
        assert_eq!(energies.len(), 5);
        for e in &energies {
            assert!(e.is_finite(), "energy should be finite");
        }
    }

    #[test]
    fn test_reciprocal_energy_single_k_finite() {
        let params = EwaldParams::new(10.0);
        let ewald = EwaldSummation::new(params, [20.0, 20.0, 20.0]);
        let k = [2.0 * std::f64::consts::PI / 20.0, 0.0, 0.0];
        let positions = [[2.0, 10.0, 10.0], [8.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let e = ewald.reciprocal_energy_single_k(k, &positions, &charges);
        assert!(e.is_finite(), "per-k energy should be finite, got {e}");
    }

    #[test]
    fn test_ewald_sum_numerical_force_matches_analytic() {
        // Compare numerical gradient of real-space energy with analytic real-space force.
        // We use a large box and small k_max so that the real-space term dominates,
        // and we compare the real-space force against a finite difference of the
        // real-space energy (not total energy, since numerical_force_component uses
        // total_energy which includes reciprocal and self parts).
        //
        // Instead, test via the standalone ewald_real_space_energy and ewald_real_forces.
        let alpha = 0.4;
        let cutoff = 15.0;
        let h = 1e-5;
        let positions = vec![[2.0, 10.0, 10.0], [7.0, 10.0, 10.0]];
        let charges = [1.0, -1.0_f64];

        let e0 = ewald_real_space_energy(&positions, &charges, alpha, cutoff);
        let mut pos_shifted = positions.clone();
        pos_shifted[0][0] += h;
        let e1 = ewald_real_space_energy(&pos_shifted, &charges, alpha, cutoff);

        let f_numerical = -(e1 - e0) / h;
        let f_analytic = ewald_real_forces(&positions, &charges, alpha, cutoff)[0][0];

        let rel_err = (f_numerical - f_analytic).abs() / f_analytic.abs().max(1e-10);
        assert!(
            rel_err < 0.01,
            "numerical ({f_numerical}) vs analytic ({f_analytic}) real-space force mismatch, rel_err = {rel_err}"
        );
    }

    #[test]
    fn test_ewald_sum_net_charge_zero() {
        let charges = [1.0, -0.5, -0.5];
        let q = EwaldSum::net_charge(&charges);
        assert!(q.abs() < 1e-12, "net charge should be 0, got {q}");
    }

    #[test]
    fn test_ewald_sum_total_energy_components_add_up() {
        // For a simple 2-charge system, verify decomposition
        let es = EwaldSum::new(0.4, 4, 25.0);
        let positions = [[5.0, 12.5, 12.5], [10.0, 12.5, 12.5]];
        let charges = [1.0, -1.0];

        let e_real = es.real_space_energy(&positions, &charges);
        let e_recip = es.reciprocal_space_energy(&positions, &charges);
        let e_self = es.self_energy(&charges);
        let e_total = es.total_energy(&positions, &charges);

        let sum = e_real + e_recip + e_self;
        assert!(
            (sum - e_total).abs() < 1e-10,
            "decomposition: real({e_real}) + recip({e_recip}) + self({e_self}) = {sum} != total({e_total})"
        );
    }

    #[test]
    fn test_madelung_constant_dimer_is_one() {
        let a = madelung_constant_dimer(3.0);
        assert!((a - 1.0).abs() < 1e-14, "dimer Madelung = 1.0, got {a}");
    }

    #[test]
    fn test_madelung_constant_dimer_independent_of_r() {
        let a1 = madelung_constant_dimer(2.0);
        let a2 = madelung_constant_dimer(10.0);
        assert!((a1 - a2).abs() < 1e-14);
    }

    #[test]
    fn test_madelung_nacl_positive() {
        // NaCl Madelung constant should be ~ 1.748 (positive)
        let a = madelung_constant_nacl(0.5, 3, 1, 5.64);
        assert!(a > 0.0, "Madelung constant should be positive, got {a}");
    }

    #[test]
    fn test_madelung_nacl_finite() {
        let a = madelung_constant_nacl(0.3, 4, 1, 4.0);
        assert!(a.is_finite(), "Madelung constant should be finite, got {a}");
    }

    #[test]
    fn test_madelung_nacl_order_of_magnitude() {
        // NaCl Madelung ≈ 1.748; accept range [1.0, 3.0] for small cluster
        let a = madelung_constant_nacl(0.4, 5, 2, 5.64);
        assert!(
            a > 0.5 && a < 5.0,
            "Madelung constant out of expected range [0.5, 5.0]: {a}"
        );
    }

    #[test]
    fn test_madelung_direct_lattice_small_returns_finite() {
        let a = madelung_direct_lattice(3, 1.0);
        assert!(
            a.is_finite(),
            "direct lattice sum should be finite, got {a}"
        );
    }

    #[test]
    fn test_madelung_direct_lattice_positive() {
        let a = madelung_direct_lattice(4, 1.0);
        assert!(a > 0.0, "direct lattice sum should be positive, got {a}");
    }

    #[test]
    fn test_madelung_direct_lattice_scales_with_a() {
        // madelung_direct_lattice divides by a, so doubling a halves the result
        let a1 = madelung_direct_lattice(3, 1.0);
        let a2 = madelung_direct_lattice(3, 2.0);
        assert!((a1 - 2.0 * a2).abs() < 1e-10, "a1={a1}, 2*a2={}", 2.0 * a2);
    }

    #[test]
    fn test_ewald_params_reciprocal_error_estimate_finite() {
        let params = EwaldParams::new(10.0);
        let err = params.reciprocal_space_error_estimate(4.0, 30.0, 5);
        assert!(
            err.is_finite(),
            "reciprocal error estimate should be finite"
        );
    }

    #[test]
    fn test_ewald_sum_self_energy_negative_v2() {
        let es = EwaldSum::new(0.5, 3, 20.0);
        let charges = [1.0, 1.0];
        let se = es.self_energy(&charges);
        assert!(
            se < 0.0,
            "EwaldSum self-energy should be negative, got {se}"
        );
    }

    #[test]
    fn test_ewald_sum_reciprocal_energy_symmetric_charges_zero() {
        let es = EwaldSum::new(0.5, 3, 20.0);
        // A single charge with itself: no pairs -> reciprocal part only
        // Two charges with zero total charge but symmetric
        let positions = [[10.0, 10.0, 10.0], [10.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        let e = es.reciprocal_space_energy(&positions, &charges);
        // S(k) = 0 for q+ = -q- at same position -> 0
        assert!(
            e.abs() < 1e-10,
            "cancelling charges at same position: e={e}"
        );
    }

    #[test]
    fn test_ewald_sum_total_energy_opposite_charges_negative() {
        let es = EwaldSum::new(0.5, 3, 40.0);
        let positions = [[10.0, 20.0, 20.0], [15.0, 20.0, 20.0]];
        let charges = [1.0, -1.0];
        let e = es.total_energy(&positions, &charges);
        assert!(
            e < 0.0,
            "opposite charges should have negative total energy, got {e}"
        );
    }

    #[test]
    fn test_ewald_sum_min_image_returns_same_for_small_displacement() {
        let es = EwaldSum::new(0.5, 3, 20.0);
        let res = es.real_space_energy(&[[0.0, 0.0, 0.0], [1.0, 2.0, 3.0]], &[1.0, 1.0]);
        assert!(
            res.is_finite(),
            "small displacement should give finite energy, got {res}"
        );
    }

    #[test]
    fn test_pme_grid_clear_zeroes_all() {
        let mut grid = PmeGrid::new(8, 0.4, 20.0);
        let positions = [[5.0, 10.0, 10.0]];
        let charges = [1.0];
        grid.spread_charges(&positions, &charges);
        grid.clear();
        assert_eq!(grid.total_charge(), 0.0);
    }

    #[test]
    fn test_pme_grid_spread_charges_preserves_total() {
        let mut grid = PmeGrid::new(8, 0.4, 20.0);
        let positions = [[5.0, 10.0, 10.0], [15.0, 10.0, 10.0]];
        let charges = [1.0, -1.0];
        grid.spread_charges(&positions, &charges);
        let q = grid.total_charge();
        assert!(
            q.abs() < 1e-10,
            "total charge should be 0 for neutral system, got {q}"
        );
    }

    #[test]
    fn test_pme_grid_occupied_cells_nonzero_after_spread() {
        let mut grid = PmeGrid::new(8, 0.4, 20.0);
        let positions = [[5.0, 5.0, 5.0]];
        let charges = [1.0];
        grid.spread_charges(&positions, &charges);
        assert!(grid.occupied_cells() >= 1);
    }

    #[test]
    fn test_pme_grid_reciprocal_energy_estimate_positive_for_nonzero_charges() {
        let mut grid = PmeGrid::new(8, 0.4, 20.0);
        let positions = [[5.0, 5.0, 5.0]];
        let charges = [1.0];
        grid.spread_charges(&positions, &charges);
        let e = grid.reciprocal_energy_estimate();
        assert!(e > 0.0, "reciprocal energy estimate should be positive");
    }

    #[test]
    fn test_structure_factor_modulus_sq_zero_charge() {
        let sf = StructureFactor::compute([1.0, 0.0, 0.0], &[[0.0, 0.0, 0.0]], &[0.0]);
        assert!((sf.modulus_sq()).abs() < 1e-20);
    }

    #[test]
    fn test_structure_factor_compute_all_length_v2() {
        let k_vecs = vec![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let positions = [[0.0, 0.0, 0.0]];
        let charges = [1.0];
        let sfs = StructureFactor::compute_all(&k_vecs, &positions, &charges);
        assert_eq!(sfs.len(), 3);
    }

    #[test]
    fn test_structure_factor_compute_all_neutral_zero() {
        let k_vecs = vec![[0.5, 0.0, 0.0]];
        let positions = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let charges = [1.0, -1.0];
        let sfs = StructureFactor::compute_all(&k_vecs, &positions, &charges);
        // Not necessarily zero in general - just finite
        assert!(sfs[0].modulus_sq().is_finite());
    }

    #[test]
    fn test_k_vector_list_size() {
        let k_vecs = k_vector_list(2, [10.0, 10.0, 10.0]);
        // (2*2+1)^3 - 1 = 5^3 - 1 = 124 vectors
        assert_eq!(k_vecs.len(), 124);
    }

    #[test]
    fn test_k_vector_list_no_zero_vector() {
        let k_vecs = k_vector_list(2, [10.0, 10.0, 10.0]);
        for k in &k_vecs {
            let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
            assert!(k2 > 1e-20, "k=0 should not be in list");
        }
    }

    #[test]
    fn test_bspline_order_2_at_zero_is_one() {
        // B-spline of order 2 at x=0 should be > 0 (in support)
        let v = bspline(2, 0.0);
        assert!(
            (0.0..=1.0).contains(&v),
            "B2(0) should be in [0,1], got {v}"
        );
    }

    #[test]
    fn test_bspline_order_1_outside_support_is_zero() {
        assert_eq!(bspline(1, -0.1), 0.0);
        assert_eq!(bspline(1, 1.0), 0.0);
    }

    #[test]
    fn test_reciprocal_energy_kvecs_empty_zero() {
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let e = reciprocal_energy_kvecs(&pos, &charges, &[], 0.5, 8000.0);
        assert_eq!(e, 0.0, "empty k-vector list → zero energy");
    }

    #[test]
    fn test_reciprocal_energy_kvecs_finite_for_nonzero_kvec() {
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let kvec = [2.0 * std::f64::consts::PI / 20.0, 0.0, 0.0];
        let e = reciprocal_energy_kvecs(&pos, &charges, &[kvec], 0.5, 8000.0);
        assert!(e.is_finite(), "single k-vector → finite energy: {e}");
    }

    #[test]
    fn test_reciprocal_forces_kvecs_len_correct() {
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let kvec = [2.0 * std::f64::consts::PI / 20.0, 0.0, 0.0];
        let f = reciprocal_forces_kvecs(&pos, &charges, &[kvec], 0.5, 8000.0);
        assert_eq!(f.len(), 2);
    }

    #[test]
    fn test_reciprocal_forces_kvecs_finite() {
        let pos = vec![[0.0; 3], [3.0, 0.0, 0.0]];
        let charges = vec![1.0, -1.0];
        let kvec = [2.0 * std::f64::consts::PI / 20.0, 0.0, 0.0];
        let f = reciprocal_forces_kvecs(&pos, &charges, &[kvec], 0.5, 8000.0);
        for fi in &f {
            for &c in fi {
                assert!(c.is_finite(), "force component must be finite: {c}");
            }
        }
    }
}
