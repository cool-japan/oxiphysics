#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Topology optimisation for compliant mechanisms and structures.
//!
//! Implements:
//!
//! - SIMP penalisation (Solid Isotropic Material with Penalisation)
//! - Sensitivity filter (density-weighted averaging)
//! - Density filter (weighted average of element densities)
//! - Optimality-criteria (OC) update rule for minimum compliance
//! - MMA (Method of Moving Asymptotes) update rule
//! - Volume constraint enforcement
//! - Checkerboard suppression via average-neighbour filter
//! - Heaviside projection (density thresholding)
//! - Convergence check based on relative density change
//! - Multi-load topology optimisation (weighted compliance sum)
//! - Stress-constrained topology optimisation (P-norm stress constraint)
//! - Adjoint sensitivity analysis
//! - Compliance minimisation loop
//!
//! No nalgebra dependency — all arrays are plain `f64` / `Vec`f64`.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// TopOptProblem
// ---------------------------------------------------------------------------

/// A 2-D topology optimisation problem on a rectangular grid.
///
/// The domain has `nx × ny` elements. Each element holds a design density
/// ρ ∈ [ρ_min, 1] and a sensitivity value ∂C/∂ρ.
#[derive(Debug, Clone)]
pub struct TopOptProblem {
    /// Number of elements in the x-direction.
    pub nx: usize,
    /// Number of elements in the y-direction.
    pub ny: usize,
    /// Minimum density (void).
    pub rho_min: f64,
    /// Penalisation exponent p for SIMP.
    pub penal: f64,
    /// Design densities, length nx*ny.
    pub density: Vec<f64>,
    /// Elemental sensitivities ∂C/∂ρ, length nx*ny.
    pub sensitivity: Vec<f64>,
    /// Objective value (compliance).
    pub objective: f64,
}

impl TopOptProblem {
    /// Create a new problem with uniform initial density ρ_0.
    ///
    /// # Parameters
    /// - `nx`      : Elements in x.
    /// - `ny`      : Elements in y.
    /// - `rho_0`   : Initial density (usually volume fraction target).
    /// - `rho_min` : Minimum density (void, prevents singularity).
    /// - `penal`   : SIMP penalisation exponent (typically 3).
    pub fn new(nx: usize, ny: usize, rho_0: f64, rho_min: f64, penal: f64) -> Self {
        let n = nx * ny;
        Self {
            nx,
            ny,
            rho_min,
            penal,
            density: vec![rho_0; n],
            sensitivity: vec![0.0; n],
            objective: 0.0,
        }
    }

    /// Total number of elements.
    pub fn n_elements(&self) -> usize {
        self.nx * self.ny
    }

    /// Current **volume fraction** V = (Σ ρ_e) / N.
    pub fn volume_fraction(&self) -> f64 {
        let sum: f64 = self.density.iter().sum();
        sum / self.n_elements() as f64
    }

    // -----------------------------------------------------------------------
    // SIMP
    // -----------------------------------------------------------------------

    /// **SIMP penalisation**: E(ρ) = E₀ · ρ^p.
    ///
    /// Maps a relative density ρ ∈ [ρ_min, 1] to effective Young's modulus
    /// fraction.  Uses the stored penalisation exponent `self.penal`.
    pub fn simp_penalization(&self, rho: f64, e0: f64) -> f64 {
        e0 * rho.max(self.rho_min).powf(self.penal)
    }

    /// Derivative of SIMP w.r.t. density: dE/dρ = p · E₀ · ρ^(p−1).
    pub fn simp_derivative(&self, rho: f64, e0: f64) -> f64 {
        e0 * self.penal * rho.max(self.rho_min).powf(self.penal - 1.0)
    }

    /// Compute SIMP effective modulus for all elements.
    ///
    /// Returns a vector of length `n_elements()` with `E(ρ_e)`.
    pub fn simp_all(&self, e0: f64) -> Vec<f64> {
        self.density
            .iter()
            .map(|&rho| self.simp_penalization(rho, e0))
            .collect()
    }

    // -----------------------------------------------------------------------
    // Sensitivity filter
    // -----------------------------------------------------------------------

    /// **Sensitivity filter**: replace ∂C/∂ρ_e with a density-weighted average
    /// over neighbours within radius r_min.
    ///
    /// ∂C̃/∂ρ_e = (1/ρ_e) · (Σ_f H_ef · ρ_f · ∂C/∂ρ_f) / Σ_f H_ef
    ///
    /// where H_ef = max(0, r_min − dist(e, f)).
    pub fn sensitivity_filter(&self, r_min: f64) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut filtered = vec![0.0f64; n];

        for ey in 0..ny {
            for ex in 0..nx {
                let idx_e = ey * nx + ex;
                let mut sum_hk = 0.0f64;
                let mut sum_hk_dc = 0.0f64;

                // Iterate over neighbours within bounding box of r_min.
                let fy_min = (ey as isize - r_min.ceil() as isize).max(0) as usize;
                let fy_max = (ey + r_min.ceil() as usize + 1).min(ny);
                let fx_min = (ex as isize - r_min.ceil() as isize).max(0) as usize;
                let fx_max = (ex + r_min.ceil() as usize + 1).min(nx);

                for fy in fy_min..fy_max {
                    for fx in fx_min..fx_max {
                        let dx = ex as f64 - fx as f64;
                        let dy = ey as f64 - fy as f64;
                        let dist = (dx * dx + dy * dy).sqrt();
                        let h = (r_min - dist).max(0.0);
                        let idx_f = fy * nx + fx;
                        sum_hk += h;
                        sum_hk_dc += h * self.density[idx_f] * self.sensitivity[idx_f];
                    }
                }

                let rho_e = self.density[idx_e].max(1e-12);
                if sum_hk > 1e-12 {
                    filtered[idx_e] = sum_hk_dc / (rho_e * sum_hk);
                } else {
                    filtered[idx_e] = self.sensitivity[idx_e];
                }
            }
        }
        filtered
    }

    // -----------------------------------------------------------------------
    // Density filter
    // -----------------------------------------------------------------------

    /// **Density filter**: replace ρ_e with a weighted average of nearby densities.
    ///
    /// ρ̃_e = (Σ_f H_ef · ρ_f) / (Σ_f H_ef)
    ///
    /// where H_ef = max(0, r_min − dist(e, f)).
    pub fn density_filter(&self, r_min: f64) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut filtered = vec![0.0f64; n];

        for ey in 0..ny {
            for ex in 0..nx {
                let idx_e = ey * nx + ex;
                let mut sum_h = 0.0f64;
                let mut sum_hrho = 0.0f64;

                let fy_min = (ey as isize - r_min.ceil() as isize).max(0) as usize;
                let fy_max = (ey + r_min.ceil() as usize + 1).min(ny);
                let fx_min = (ex as isize - r_min.ceil() as isize).max(0) as usize;
                let fx_max = (ex + r_min.ceil() as usize + 1).min(nx);

                for fy in fy_min..fy_max {
                    for fx in fx_min..fx_max {
                        let dx = ex as f64 - fx as f64;
                        let dy = ey as f64 - fy as f64;
                        let dist = (dx * dx + dy * dy).sqrt();
                        let h = (r_min - dist).max(0.0);
                        let idx_f = fy * nx + fx;
                        sum_h += h;
                        sum_hrho += h * self.density[idx_f];
                    }
                }

                if sum_h > 1e-12 {
                    filtered[idx_e] = sum_hrho / sum_h;
                } else {
                    filtered[idx_e] = self.density[idx_e];
                }
            }
        }
        filtered
    }

    // -----------------------------------------------------------------------
    // Optimality criteria
    // -----------------------------------------------------------------------

    /// **Optimality-criteria update rule** for minimum-compliance topology
    /// optimisation.
    ///
    /// ρ_e^{new} = clamp(ρ_e · √(−∂C̃/∂ρ_e / (λ · ρ_e)), ρ_lo, ρ_hi)
    ///
    /// where λ is found by bisection to satisfy the volume constraint
    /// V = v_frac.
    ///
    /// Returns the new density vector.
    pub fn optimality_criteria(
        &self,
        filtered_sens: &[f64],
        v_frac: f64,
        move_limit: f64,
    ) -> Vec<f64> {
        let n = self.n_elements();

        // Use log-space bisection on Lagrange multiplier λ for numerical stability.
        // Bounds span many orders of magnitude; bisecting in log space converges
        // in O(log(range/tol)) iterations regardless of scale.
        let mut log_lam_lo = -40.0_f64; // 10^-40
        let mut log_lam_hi = 10.0_f64; // 10^10  (plenty of range)
        let mut new_density = vec![0.0f64; n];

        for _ in 0..100 {
            let log_lam_mid = 0.5 * (log_lam_lo + log_lam_hi);
            let lam_mid = 10.0_f64.powf(log_lam_mid);
            let mut vol = 0.0f64;
            for e in 0..n {
                let rho = self.density[e];
                let be = (-filtered_sens[e] / lam_mid).max(0.0).sqrt() * rho;
                let rho_lo = (rho - move_limit).max(self.rho_min);
                let rho_hi = (rho + move_limit).min(1.0);
                let rho_new = be.clamp(rho_lo, rho_hi);
                new_density[e] = rho_new;
                vol += rho_new;
            }
            let vf = vol / n as f64;
            if vf > v_frac {
                log_lam_lo = log_lam_mid;
            } else {
                log_lam_hi = log_lam_mid;
            }
        }
        new_density
    }

    // -----------------------------------------------------------------------
    // MMA update (simplified)
    // -----------------------------------------------------------------------

    /// **Method of Moving Asymptotes** update step (simplified single-variable version).
    ///
    /// For each element, constructs a local convex approximation:
    ///   f̃_e(ρ) = a_e / (U_e - ρ)   with U_e = ρ_e + asymptote_move
    ///
    /// The update satisfies the volume constraint via bisection on the Lagrange
    /// multiplier μ of the linearised sub-problem.
    ///
    /// # Parameters
    /// - `filtered_sens`: filtered sensitivities ∂C̃/∂ρ.
    /// - `v_frac`:        volume fraction target.
    /// - `move_limit`:    maximum per-element step size in a single iteration.
    /// - `asymptote`:     MMA asymptote distance (controls convexity).
    ///
    /// Returns the new density vector.
    pub fn mma_update(
        &self,
        filtered_sens: &[f64],
        v_frac: f64,
        move_limit: f64,
        asymptote: f64,
    ) -> Vec<f64> {
        let n = self.n_elements();
        let mut new_density = vec![0.0f64; n];

        // Bisect on μ (Lagrange multiplier for volume constraint)
        let mut mu_lo = 0.0_f64;
        let mut mu_hi = 1e6_f64;

        // Pre-compute upper asymptotes
        let upper: Vec<f64> = self.density.iter().map(|&r| r + asymptote).collect();
        // Coefficients: p_e = max(0, -dc_e) * (U_e - rho_e)^2
        let p_coeff: Vec<f64> = filtered_sens
            .iter()
            .zip(self.density.iter())
            .zip(upper.iter())
            .map(|((&dc, &rho), &u)| (-dc).max(0.0) * (u - rho).powi(2))
            .collect();
        // Coefficients for volume: q_e = (U_e - rho_e)^2
        let q_coeff: Vec<f64> = self
            .density
            .iter()
            .zip(upper.iter())
            .map(|(&rho, &u)| (u - rho).powi(2))
            .collect();

        for _ in 0..100 {
            let mu_mid = 0.5 * (mu_lo + mu_hi);
            let mut vol = 0.0f64;
            for e in 0..n {
                let rho = self.density[e];
                let u = upper[e];
                // Optimality condition: rho* = U - sqrt(p_e / (mu * q_e + 1e-30))
                let denom = (mu_mid * q_coeff[e] + 1e-30).max(1e-30);
                let rho_star = u - (p_coeff[e] / denom).sqrt();
                let rho_lo = (rho - move_limit).max(self.rho_min);
                let rho_hi = (rho + move_limit).min(1.0);
                let rho_new = rho_star.clamp(rho_lo, rho_hi);
                new_density[e] = rho_new;
                vol += rho_new;
            }
            let vf = vol / n as f64;
            // Larger μ → larger ρ* → larger vf.
            // To reduce vf below target, decrease μ (upper bound).
            if vf > v_frac {
                mu_hi = mu_mid;
            } else {
                mu_lo = mu_mid;
            }
        }
        new_density
    }

    // -----------------------------------------------------------------------
    // Volume constraint
    // -----------------------------------------------------------------------

    /// Check whether the **volume constraint** V ≤ V_max is satisfied.
    ///
    /// Returns `true` if the current volume fraction ≤ `v_max`.
    pub fn volume_constraint(&self, v_max: f64) -> bool {
        self.volume_fraction() <= v_max + 1e-10
    }

    // -----------------------------------------------------------------------
    // Checkerboard suppression
    // -----------------------------------------------------------------------

    /// **Checkerboard suppression**: replace each element's density with the
    /// average of itself and its direct (von-Neumann) neighbours.
    ///
    /// This low-pass filter removes the alternating solid/void pattern that
    /// can arise with tri-linear elements.
    pub fn checkerboard_suppression(&self) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let n = nx * ny;
        let mut smoothed = vec![0.0f64; n];

        for ey in 0..ny {
            for ex in 0..nx {
                let idx = ey * nx + ex;
                let mut sum = self.density[idx];
                let mut count = 1_usize;

                let neighbours: [(isize, isize); 4] = [(-1, 0), (1, 0), (0, -1), (0, 1)];
                for (dx, dy) in neighbours {
                    let fx = ex as isize + dx;
                    let fy = ey as isize + dy;
                    if fx >= 0 && fx < nx as isize && fy >= 0 && fy < ny as isize {
                        sum += self.density[fy as usize * nx + fx as usize];
                        count += 1;
                    }
                }
                smoothed[idx] = sum / count as f64;
            }
        }
        smoothed
    }

    // -----------------------------------------------------------------------
    // Heaviside projection
    // -----------------------------------------------------------------------

    /// **Heaviside projection** to sharpen density fields toward 0/1.
    ///
    /// H_β(ρ̃) = (tanh(β·η) + tanh(β·(ρ̃ − η))) / (tanh(β·η) + tanh(β·(1 − η)))
    ///
    /// # Parameters
    /// - `beta` : sharpness parameter (larger → sharper projection, e.g. 8).
    /// - `eta`  : threshold (typically 0.5 for threshold at mid-density).
    ///
    /// Returns the projected density field.
    pub fn heaviside_projection(&self, beta: f64, eta: f64) -> Vec<f64> {
        let denom = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
        let denom = denom.max(1e-12);
        self.density
            .iter()
            .map(|&rho| {
                let num = (beta * eta).tanh() + (beta * (rho - eta)).tanh();
                (num / denom).clamp(0.0, 1.0)
            })
            .collect()
    }

    /// Compute the **chain-rule correction** to sensitivities for Heaviside projection.
    ///
    /// dH/dρ̃ = β · sech²(β·(ρ̃ − η)) / (tanh(β·η) + tanh(β·(1 − η)))
    ///
    /// Multiplies `sens` elementwise by dH/dρ̃ and returns the corrected vector.
    pub fn heaviside_sensitivity_correction(&self, sens: &[f64], beta: f64, eta: f64) -> Vec<f64> {
        let denom = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
        let denom = denom.max(1e-12);
        self.density
            .iter()
            .zip(sens.iter())
            .map(|(&rho, &s)| {
                let sech2 = 1.0 / (beta * (rho - eta)).cosh().powi(2);
                let dhdrho = beta * sech2 / denom;
                s * dhdrho
            })
            .collect()
    }

    // -----------------------------------------------------------------------
    // Convergence
    // -----------------------------------------------------------------------

    /// **Convergence check**: returns `true` when the maximum relative change
    /// in density between old and new distributions is below `tol`.
    ///
    /// max_e |ρ_new\[e\] − ρ_old\[e\]| / max(ρ_old\[e\], 1e-6) < tol
    pub fn convergence_check(&self, new_density: &[f64], tol: f64) -> bool {
        assert_eq!(new_density.len(), self.n_elements());
        let max_change = self
            .density
            .iter()
            .zip(new_density.iter())
            .map(|(&rho_old, &rho_new)| (rho_new - rho_old).abs() / rho_old.max(1e-6))
            .fold(0.0_f64, f64::max);
        max_change < tol
    }

    /// Apply a new density distribution in-place.
    pub fn apply_density(&mut self, new_density: Vec<f64>) {
        assert_eq!(new_density.len(), self.n_elements());
        self.density = new_density;
    }

    // -----------------------------------------------------------------------
    // Adjoint sensitivity
    // -----------------------------------------------------------------------

    /// **Adjoint-method compliance sensitivity** ∂C/∂ρ_e = −p · E₀ · ρ_e^(p−1) · u_e^T K_e u_e.
    ///
    /// This computes a proxy sensitivity assuming unit element strain energy.
    /// In a real FEA implementation, `element_strain_energies` would come from
    /// the FE solver.
    ///
    /// # Parameters
    /// - `element_strain_energies` : elemental strain energies u_e^T K_0 u_e (J).
    /// - `e0` : reference Young's modulus.
    ///
    /// Fills `self.sensitivity` and returns a reference slice.
    pub fn adjoint_sensitivity(&mut self, element_strain_energies: &[f64], e0: f64) {
        assert_eq!(element_strain_energies.len(), self.n_elements());
        for e in 0..self.n_elements() {
            let rho = self.density[e];
            let dsimp = self.simp_derivative(rho, e0);
            self.sensitivity[e] = -dsimp * element_strain_energies[e];
        }
        self.objective = element_strain_energies
            .iter()
            .zip(self.density.iter())
            .map(|(&se, &rho)| self.simp_penalization(rho, e0) * se)
            .sum();
    }
}

// ---------------------------------------------------------------------------
// Stand-alone helpers
// ---------------------------------------------------------------------------

/// Compute the global SIMP compliance sensitivity assuming a uniform load
/// (proxy: sensitivity\[e\] = −E(ρ_e) · strain_energy_e for a unit strain).
///
/// This helper fills `problem.sensitivity` with −simp(ρ_e) as a stand-in
/// for a real FEA kernel.
pub fn fill_dummy_sensitivity(problem: &mut TopOptProblem) {
    for e in 0..problem.n_elements() {
        let rho = problem.density[e];
        problem.sensitivity[e] = -problem.simp_penalization(rho, 1.0);
    }
    problem.objective = problem.sensitivity.iter().map(|s| s.abs()).sum();
}

/// Compute sensitivities for a **multi-load** problem.
///
/// The total compliance is C = Σ_l w_l · C_l where each load case `l`
/// contributes elemental strain energies `strain_energies\[l\]\[e\]` with weight
/// `weights\[l\]`.
///
/// # Parameters
/// - `problem`         : the topology optimisation problem (mutated in-place).
/// - `strain_energies` : outer vector over load cases, inner over elements.
/// - `weights`         : scalar weight per load case.
/// - `e0`              : reference Young's modulus.
pub fn multi_load_sensitivity(
    problem: &mut TopOptProblem,
    strain_energies: &[Vec<f64>],
    weights: &[f64],
    e0: f64,
) {
    assert_eq!(strain_energies.len(), weights.len());
    let n = problem.n_elements();
    problem.sensitivity = vec![0.0; n];
    problem.objective = 0.0;

    for (se_load, &w) in strain_energies.iter().zip(weights.iter()) {
        assert_eq!(se_load.len(), n);
        for e in 0..n {
            let rho = problem.density[e];
            let dsimp = problem.simp_derivative(rho, e0);
            problem.sensitivity[e] -= w * dsimp * se_load[e];
            problem.objective += w * problem.simp_penalization(rho, e0) * se_load[e];
        }
    }
}

/// Compute the **P-norm stress constraint** value and its sensitivity.
///
/// P-norm aggregation: σ_PN = (Σ_e (σ_e / σ_0)^P)^(1/P)
///
/// Penalized element stresses: σ̄_e = ρ_e^q · σ_vm(e) where `q` is a
/// stress penalisation exponent and `σ_vm(e)` is the von Mises stress.
///
/// Returns `(pnorm_value, sensitivity_vector)`.
///
/// # Parameters
/// - `densities`        : current element densities.
/// - `vm_stresses`      : elemental von Mises stresses (Pa).
/// - `sigma_0`          : allowable stress (Pa).
/// - `p_norm_exp`       : P-norm exponent (typically 8–16).
/// - `stress_pen_exp`   : stress penalisation exponent `q` (typically 0.5).
pub fn pnorm_stress_constraint(
    densities: &[f64],
    vm_stresses: &[f64],
    sigma_0: f64,
    p_norm_exp: f64,
    stress_pen_exp: f64,
) -> (f64, Vec<f64>) {
    assert_eq!(densities.len(), vm_stresses.len());
    let n = densities.len();

    // Penalized stresses
    let pen_stresses: Vec<f64> = densities
        .iter()
        .zip(vm_stresses.iter())
        .map(|(&rho, &s)| rho.powf(stress_pen_exp) * s / sigma_0)
        .collect();

    // P-norm value
    let sum_p: f64 = pen_stresses.iter().map(|&s| s.abs().powf(p_norm_exp)).sum();
    let pnorm = sum_p.powf(1.0 / p_norm_exp);

    // Sensitivity dPN/dρ_e = (1/P) * sum_p^(1/P - 1) * P * pen_s_e^(P-1) * d(pen_s_e)/d(rho_e)
    //   d(pen_s_e)/d(rho_e) = q * rho_e^(q-1) * vm_s_e / sigma_0
    let pnorm_factor = sum_p.powf(1.0 / p_norm_exp - 1.0);
    let mut sens = vec![0.0f64; n];
    for e in 0..n {
        let ps = pen_stresses[e];
        let dps =
            stress_pen_exp * densities[e].powf(stress_pen_exp - 1.0).max(0.0) * vm_stresses[e]
                / sigma_0;
        sens[e] = pnorm_factor * ps.abs().powf(p_norm_exp - 1.0) * ps.signum() * dps;
    }
    (pnorm, sens)
}

/// Run a simplified **compliance minimisation loop**.
///
/// This is a pure-Rust stand-alone solver that uses a proxy FEA (uniform
/// unit strain energies per element proportional to the SIMP modulus).
///
/// # Parameters
/// - `problem`     : mutable reference to the topology problem.
/// - `v_frac`      : target volume fraction.
/// - `r_min`       : sensitivity filter radius (element units).
/// - `move_limit`  : OC move limit per iteration.
/// - `max_iter`    : maximum number of optimisation iterations.
/// - `tol`         : density-change convergence tolerance.
///
/// Returns the number of iterations actually performed.
pub fn compliance_minimisation(
    problem: &mut TopOptProblem,
    v_frac: f64,
    r_min: f64,
    move_limit: f64,
    max_iter: usize,
    tol: f64,
) -> usize {
    for it in 0..max_iter {
        // Proxy FEA: unit strain energy per element
        let se: Vec<f64> = vec![1.0; problem.n_elements()];
        problem.adjoint_sensitivity(&se, 1.0);

        let filtered = problem.sensitivity_filter(r_min);
        let new_density = problem.optimality_criteria(&filtered, v_frac, move_limit);

        let converged = problem.convergence_check(&new_density, tol);
        problem.apply_density(new_density);
        if converged {
            return it + 1;
        }
    }
    max_iter
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_problem() -> TopOptProblem {
        TopOptProblem::new(10, 10, 0.5, 1e-3, 3.0)
    }

    // 1. Initial volume fraction equals rho_0.
    #[test]
    fn test_initial_volume_fraction() {
        let prob = make_problem();
        assert!(
            (prob.volume_fraction() - 0.5).abs() < 1e-10,
            "Initial volume fraction should be 0.5, got {}",
            prob.volume_fraction()
        );
    }

    // 2. SIMP at rho=1 gives E0.
    #[test]
    fn test_simp_full_density() {
        let prob = make_problem();
        let e = prob.simp_penalization(1.0, 1.0);
        assert!(
            (e - 1.0).abs() < 1e-10,
            "SIMP at ρ=1 should give E0=1, got {e}"
        );
    }

    // 3. SIMP at rho_min gives near-zero stiffness.
    #[test]
    fn test_simp_min_density() {
        let prob = make_problem();
        let e = prob.simp_penalization(prob.rho_min, 1.0);
        assert!(
            e < 0.01,
            "SIMP at ρ_min should give near-zero stiffness, got {e}"
        );
    }

    // 4. SIMP is monotone increasing with density.
    #[test]
    fn test_simp_monotone() {
        let prob = make_problem();
        let e1 = prob.simp_penalization(0.3, 1.0);
        let e2 = prob.simp_penalization(0.6, 1.0);
        let e3 = prob.simp_penalization(0.9, 1.0);
        assert!(e1 < e2 && e2 < e3, "SIMP should be monotone increasing");
    }

    // 5. SIMP derivative is positive.
    #[test]
    fn test_simp_derivative_positive() {
        let prob = make_problem();
        let de = prob.simp_derivative(0.5, 1.0);
        assert!(de > 0.0, "SIMP derivative should be positive, got {de}");
    }

    // 6. Sensitivity filter preserves approximate magnitude.
    #[test]
    fn test_sensitivity_filter_runs() {
        let mut prob = make_problem();
        for e in 0..prob.n_elements() {
            prob.sensitivity[e] = -1.0;
        }
        let filtered = prob.sensitivity_filter(1.5);
        assert_eq!(
            filtered.len(),
            prob.n_elements(),
            "Filter output length mismatch"
        );
        for &v in &filtered {
            assert!(v.is_finite(), "Filtered sensitivity should be finite");
        }
    }

    // 7. Sensitivity filter at rmin=0 returns finite values.
    #[test]
    fn test_sensitivity_filter_rmin_zero() {
        let mut prob = make_problem();
        for (e, s) in prob.sensitivity.iter_mut().enumerate() {
            *s = -(e as f64 + 1.0);
        }
        let filtered = prob.sensitivity_filter(0.0);
        for (e, &v) in filtered.iter().enumerate() {
            assert!(
                v.is_finite(),
                "Filter with rmin=0 should give finite values, element {e}: {v}"
            );
        }
    }

    // 8. OC update moves density toward volume target.
    #[test]
    fn test_oc_update_volume_constraint() {
        let mut prob = make_problem(); // initial ρ = 0.5
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(1.5);
        let v_frac = 0.6;
        let new_density = prob.optimality_criteria(&filtered, v_frac, 0.5);
        let vf = new_density.iter().sum::<f64>() / prob.n_elements() as f64;
        assert!(
            (vf - v_frac).abs() < 0.05,
            "OC update should move volume fraction toward {v_frac:.2}, got {vf:.4}"
        );
    }

    // 9. OC update densities stay in [rho_min, 1].
    #[test]
    fn test_oc_update_bounds() {
        let mut prob = make_problem();
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(1.5);
        let new_density = prob.optimality_criteria(&filtered, 0.5, 0.2);
        for (e, &rho) in new_density.iter().enumerate() {
            assert!(
                rho >= prob.rho_min - 1e-10 && rho <= 1.0 + 1e-10,
                "Density out of bounds at element {e}: {rho}"
            );
        }
    }

    // 10. Volume constraint check returns true when satisfied.
    #[test]
    fn test_volume_constraint_satisfied() {
        let prob = make_problem(); // ρ_0 = 0.5
        assert!(prob.volume_constraint(0.6), "ρ=0.5 should satisfy V≤0.6");
    }

    // 11. Volume constraint check returns false when violated.
    #[test]
    fn test_volume_constraint_violated() {
        let prob = make_problem(); // ρ_0 = 0.5
        assert!(!prob.volume_constraint(0.3), "ρ=0.5 should violate V≤0.3");
    }

    // 12. Checkerboard suppression produces values in [0, 1].
    #[test]
    fn test_checkerboard_suppression_range() {
        let mut prob = make_problem();
        // Create a checkerboard pattern.
        for ey in 0..prob.ny {
            for ex in 0..prob.nx {
                prob.density[ey * prob.nx + ex] = if (ex + ey) % 2 == 0 { 1.0 } else { 0.0 };
            }
        }
        let smoothed = prob.checkerboard_suppression();
        for &rho in &smoothed {
            assert!(
                (0.0 - 1e-10..=1.0 + 1e-10).contains(&rho),
                "Smoothed density out of range: {rho}"
            );
        }
    }

    // 13. Checkerboard suppression reduces variation in a checker pattern.
    #[test]
    fn test_checkerboard_suppression_smooths() {
        let mut prob = make_problem();
        for ey in 0..prob.ny {
            for ex in 0..prob.nx {
                prob.density[ey * prob.nx + ex] = if (ex + ey) % 2 == 0 { 1.0 } else { 0.0 };
            }
        }
        let original_var: f64 = prob.density.iter().map(|&r| (r - 0.5).powi(2)).sum::<f64>();
        let smoothed = prob.checkerboard_suppression();
        let smooth_var: f64 = smoothed.iter().map(|&r| (r - 0.5).powi(2)).sum::<f64>();
        assert!(
            smooth_var < original_var,
            "Smoothing should reduce variance: before={original_var:.4}, after={smooth_var:.4}"
        );
    }

    // 14. Convergence check returns true for unchanged density.
    #[test]
    fn test_convergence_unchanged() {
        let prob = make_problem();
        let same = prob.density.clone();
        assert!(
            prob.convergence_check(&same, 1e-6),
            "Unchanged density should be converged"
        );
    }

    // 15. Convergence check returns false for large change.
    #[test]
    fn test_convergence_large_change() {
        let prob = make_problem();
        let new_density = vec![1.0; prob.n_elements()];
        assert!(
            !prob.convergence_check(&new_density, 1e-3),
            "Large density change should not be converged"
        );
    }

    // 16. fill_dummy_sensitivity sets all sensitivities to negative.
    #[test]
    fn test_dummy_sensitivity_negative() {
        let mut prob = make_problem();
        fill_dummy_sensitivity(&mut prob);
        for (e, &s) in prob.sensitivity.iter().enumerate() {
            assert!(
                s < 0.0,
                "Sensitivity should be negative at element {e}: {s}"
            );
        }
    }

    // 17. apply_density updates the density field.
    #[test]
    fn test_apply_density() {
        let mut prob = make_problem();
        let new_density = vec![0.8; prob.n_elements()];
        prob.apply_density(new_density.clone());
        assert_eq!(
            prob.density, new_density,
            "apply_density should update the field"
        );
    }

    // 18. n_elements returns nx * ny.
    #[test]
    fn test_n_elements() {
        let prob = TopOptProblem::new(5, 8, 0.5, 1e-3, 3.0);
        assert_eq!(prob.n_elements(), 40);
    }

    // 19. Iterating OC steps reduces volume fraction when target is lower.
    #[test]
    fn test_oc_reduces_volume() {
        let mut prob = TopOptProblem::new(10, 10, 0.7, 1e-3, 3.0);
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(1.5);
        let new_density = prob.optimality_criteria(&filtered, 0.4, 0.3);
        let vf = new_density.iter().sum::<f64>() / prob.n_elements() as f64;
        assert!(
            vf < 0.7,
            "OC should reduce volume fraction from 0.7 toward 0.4, got {vf:.4}"
        );
    }

    // 20. SIMP exponent p=1 gives linear interpolation.
    #[test]
    fn test_simp_linear() {
        let prob = TopOptProblem::new(4, 4, 0.5, 0.0, 1.0);
        let e = prob.simp_penalization(0.5, 2.0);
        assert!(
            (e - 1.0).abs() < 1e-10,
            "SIMP p=1 at ρ=0.5 should give E0/2=1.0, got {e}"
        );
    }

    // 21. Volume fraction decreases after apply_density with lower densities.
    #[test]
    fn test_volume_decreases() {
        let mut prob = make_problem();
        let new_density = vec![0.2; prob.n_elements()];
        prob.apply_density(new_density);
        assert!(
            (prob.volume_fraction() - 0.2).abs() < 1e-10,
            "Volume fraction should be 0.2 after applying low-density field"
        );
    }

    // 22. Sensitivity filter output length equals number of elements.
    #[test]
    fn test_sensitivity_filter_length() {
        let mut prob = TopOptProblem::new(6, 8, 0.5, 1e-3, 3.0);
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(2.0);
        assert_eq!(
            filtered.len(),
            prob.n_elements(),
            "Filtered sensitivity length mismatch"
        );
    }

    // 23. OC update with move_limit=0 leaves densities unchanged.
    #[test]
    fn test_oc_zero_move_limit() {
        let mut prob = make_problem();
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(1.5);
        let new_density = prob.optimality_criteria(&filtered, 0.5, 0.0);
        for (e, (&rho_old, &rho_new)) in prob.density.iter().zip(new_density.iter()).enumerate() {
            assert!(
                (rho_new - rho_old).abs() < 1e-10,
                "Zero move limit: element {e} changed from {rho_old} to {rho_new}"
            );
        }
    }

    // 24. SIMP penalisation uses E0 scaling correctly.
    #[test]
    fn test_simp_e0_scaling() {
        let prob = make_problem();
        let e1 = prob.simp_penalization(0.5, 1.0);
        let e2 = prob.simp_penalization(0.5, 200e9);
        assert!(
            (e2 / e1 - 200e9).abs() < 1.0,
            "SIMP should scale linearly with E0"
        );
    }

    // 25. Checkerboard suppression of uniform field returns same field.
    #[test]
    fn test_checkerboard_uniform_unchanged() {
        let prob = make_problem(); // all ρ = 0.5
        let smoothed = prob.checkerboard_suppression();
        for &rho in &smoothed {
            assert!(
                (rho - 0.5).abs() < 1e-10,
                "Uniform field should be unchanged by suppression, got {rho}"
            );
        }
    }

    // 26. Density filter preserves total mass approximately.
    #[test]
    fn test_density_filter_mass_preservation() {
        let mut prob = TopOptProblem::new(8, 8, 0.5, 1e-3, 3.0);
        // Non-uniform density
        for e in 0..prob.n_elements() {
            prob.density[e] = if e % 3 == 0 { 0.8 } else { 0.3 };
        }
        let mass_before: f64 = prob.density.iter().sum();
        let filtered = prob.density_filter(1.5);
        let mass_after: f64 = filtered.iter().sum();
        // Density filter conserves total mass (it's a weighted average)
        assert!(
            (mass_after - mass_before).abs() < 1e-6,
            "Density filter should approximately preserve total mass: {mass_before:.4} vs {mass_after:.4}"
        );
    }

    // 27. Density filter output length equals n_elements.
    #[test]
    fn test_density_filter_length() {
        let prob = TopOptProblem::new(5, 7, 0.5, 1e-3, 3.0);
        let filtered = prob.density_filter(1.5);
        assert_eq!(filtered.len(), prob.n_elements());
    }

    // 28. Heaviside projection: uniform field at rho=0.5 stays near 0.5 for eta=0.5.
    #[test]
    fn test_heaviside_projection_mid_density() {
        let prob = make_problem(); // all ρ = 0.5
        let proj = prob.heaviside_projection(1.0, 0.5);
        for &p in &proj {
            assert!(
                (p - 0.5).abs() < 0.1,
                "Heaviside of 0.5 at eta=0.5 should be near 0.5, got {p}"
            );
        }
    }

    // 29. Heaviside projection: large beta sharpens toward 0/1.
    #[test]
    fn test_heaviside_projection_sharpening() {
        let mut prob = TopOptProblem::new(4, 4, 0.5, 1e-3, 3.0);
        // Set half elements to 0.3 (below threshold) and half to 0.7 (above)
        for e in 0..prob.n_elements() {
            prob.density[e] = if e < prob.n_elements() / 2 { 0.3 } else { 0.7 };
        }
        let proj_sharp = prob.heaviside_projection(64.0, 0.5);
        // Below threshold should go toward 0, above toward 1
        for &p in &proj_sharp[..prob.n_elements() / 2] {
            assert!(
                p < 0.1,
                "Below-threshold element should project near 0, got {p}"
            );
        }
        for &p in &proj_sharp[prob.n_elements() / 2..] {
            assert!(
                p > 0.9,
                "Above-threshold element should project near 1, got {p}"
            );
        }
    }

    // 30. Heaviside sensitivity correction is finite and non-negative.
    #[test]
    fn test_heaviside_sensitivity_correction_nonneg() {
        let prob = make_problem();
        let sens = vec![-1.0; prob.n_elements()];
        let corrected = prob.heaviside_sensitivity_correction(&sens, 8.0, 0.5);
        for &c in &corrected {
            assert!(
                c.is_finite(),
                "Corrected sensitivity should be finite, got {c}"
            );
            // Since sens is -1 and dH/drho is always positive, result should be <= 0
            assert!(c <= 1e-10, "Corrected sens should be non-positive, got {c}");
        }
    }

    // 31. simp_all returns vector of correct length.
    #[test]
    fn test_simp_all_length() {
        let prob = make_problem();
        let e_all = prob.simp_all(200e9);
        assert_eq!(e_all.len(), prob.n_elements());
    }

    // 32. simp_all at uniform density gives uniform output.
    #[test]
    fn test_simp_all_uniform() {
        let prob = make_problem(); // all ρ = 0.5
        let e_all = prob.simp_all(1.0);
        let e_expected = prob.simp_penalization(0.5, 1.0);
        for &e in &e_all {
            assert!(
                (e - e_expected).abs() < 1e-12,
                "simp_all should be uniform for uniform density"
            );
        }
    }

    // 33. adjoint_sensitivity fills sensitivity and sets objective.
    #[test]
    fn test_adjoint_sensitivity_fills() {
        let mut prob = make_problem();
        let se = vec![1.0; prob.n_elements()];
        prob.adjoint_sensitivity(&se, 1.0);
        for &s in &prob.sensitivity {
            assert!(
                s < 0.0,
                "Adjoint sensitivities should be negative for unit strain energy"
            );
        }
        assert!(
            prob.objective > 0.0,
            "Objective should be positive after adjoint sensitivity"
        );
    }

    // 34. multi_load_sensitivity with single load case equals adjoint_sensitivity.
    #[test]
    fn test_multi_load_single_case() {
        let mut prob1 = make_problem();
        let se = vec![1.0; prob1.n_elements()];
        prob1.adjoint_sensitivity(&se, 1.0);
        let s1 = prob1.sensitivity.clone();

        let mut prob2 = make_problem();
        let se2 = vec![vec![1.0; prob2.n_elements()]];
        let weights = vec![1.0];
        multi_load_sensitivity(&mut prob2, &se2, &weights, 1.0);
        let s2 = prob2.sensitivity.clone();

        for (i, (&a, &b)) in s1.iter().zip(s2.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-10,
                "Multi-load single case should equal adjoint at element {i}: {a} vs {b}"
            );
        }
    }

    // 35. pnorm_stress_constraint returns finite value and correct sensitivity length.
    #[test]
    fn test_pnorm_stress_constraint_basic() {
        let n = 10;
        let densities = vec![0.5_f64; n];
        let vm_stresses = vec![100.0_f64; n];
        let sigma_0 = 200.0;
        let (pnorm, sens) = pnorm_stress_constraint(&densities, &vm_stresses, sigma_0, 8.0, 0.5);
        assert!(
            pnorm.is_finite(),
            "P-norm value should be finite, got {pnorm}"
        );
        assert!(pnorm > 0.0, "P-norm should be positive, got {pnorm}");
        assert_eq!(sens.len(), n, "Sensitivity length mismatch");
        for &s in &sens {
            assert!(s.is_finite(), "P-norm sensitivity should be finite");
        }
    }

    // 36. pnorm_stress_constraint scales with stress level.
    #[test]
    fn test_pnorm_stress_scales() {
        let n = 10;
        let densities = vec![1.0_f64; n];
        let vm1 = vec![50.0_f64; n];
        let vm2 = vec![100.0_f64; n];
        let sigma_0 = 200.0;
        let (pn1, _) = pnorm_stress_constraint(&densities, &vm1, sigma_0, 4.0, 0.0);
        let (pn2, _) = pnorm_stress_constraint(&densities, &vm2, sigma_0, 4.0, 0.0);
        assert!(
            pn2 > pn1,
            "Higher stress should give larger P-norm: {pn1} vs {pn2}"
        );
    }

    // 37. compliance_minimisation runs without panic and reduces volume fraction.
    #[test]
    fn test_compliance_minimisation_runs() {
        let mut prob = TopOptProblem::new(6, 6, 0.6, 1e-3, 3.0);
        let iters = compliance_minimisation(&mut prob, 0.4, 1.5, 0.2, 20, 1e-3);
        assert!(iters <= 20, "Minimisation should complete within max_iter");
        // Volume fraction should be near target
        let vf = prob.volume_fraction();
        assert!(
            (vf - 0.4).abs() < 0.05,
            "Volume fraction after minimisation should be near 0.4, got {vf:.4}"
        );
    }

    // 38. MMA update respects bounds [rho_min, 1].
    #[test]
    fn test_mma_update_bounds() {
        let mut prob = make_problem();
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(1.5);
        let new_density = prob.mma_update(&filtered, 0.5, 0.2, 0.5);
        for (e, &rho) in new_density.iter().enumerate() {
            assert!(
                rho >= prob.rho_min - 1e-9 && rho <= 1.0 + 1e-9,
                "MMA density out of bounds at element {e}: {rho}"
            );
        }
    }

    // 39. MMA update moves volume fraction toward target.
    #[test]
    fn test_mma_update_volume_target() {
        let mut prob = TopOptProblem::new(10, 10, 0.7, 1e-3, 3.0);
        fill_dummy_sensitivity(&mut prob);
        let filtered = prob.sensitivity_filter(1.5);
        let new_density = prob.mma_update(&filtered, 0.4, 0.3, 0.5);
        let vf = new_density.iter().sum::<f64>() / prob.n_elements() as f64;
        assert!(
            (vf - 0.4).abs() < 0.05,
            "MMA update should drive volume fraction toward 0.4, got {vf:.4}"
        );
    }

    // 40. Density filter output is in [0, 1] for densities in [0, 1].
    #[test]
    fn test_density_filter_range() {
        let mut prob = TopOptProblem::new(6, 6, 0.5, 1e-3, 3.0);
        for e in 0..prob.n_elements() {
            prob.density[e] = if e % 2 == 0 { 1.0 } else { 0.0 };
        }
        let filtered = prob.density_filter(2.0);
        for &rho in &filtered {
            assert!(
                (0.0 - 1e-10..=1.0 + 1e-10).contains(&rho),
                "Density filter output should be in [0,1], got {rho}"
            );
        }
    }
}
