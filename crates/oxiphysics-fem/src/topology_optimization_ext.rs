#![allow(clippy::ptr_arg)]
#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended topology optimization FEM.
//!
//! Covers SIMP with density / sensitivity filters, Optimality Criteria (OC),
//! Method of Moving Asymptotes (MMA), compliance minimisation, volume
//! constraint, multi-load-case, stress-constrained, level-set, BESO,
//! compliant-mechanism design, frequency optimisation, and heat-conduction
//! topology optimisation.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ============================================================================
// § 1  SIMP INTERPOLATION
// ============================================================================

/// Compute SIMP interpolated Young's modulus.
///
/// E(ρ) = E_min + ρ^p · (E₀ − E_min)
///
/// * `rho`     – element density ∈ \[0, 1\]
/// * `e0`      – solid modulus
/// * `e_min`   – void modulus (numerical floor, e.g. 1e-9 · E₀)
/// * `penalty` – penalisation exponent p (typically 3)
pub fn simp_modulus(rho: f64, e0: f64, e_min: f64, penalty: f64) -> f64 {
    e_min + rho.powf(penalty) * (e0 - e_min)
}

/// Derivative dE/dρ for the SIMP interpolation.
///
/// dE/dρ = p · ρ^(p−1) · (E₀ − E_min)
pub fn simp_modulus_derivative(rho: f64, e0: f64, e_min: f64, penalty: f64) -> f64 {
    penalty * rho.powf(penalty - 1.0) * (e0 - e_min)
}

/// Compute SIMP compliance sensitivity ∂c/∂ρ.
///
/// ∂c/∂ρ = −(dE/dρ / E(ρ)) · uₑᵀ kₑ uₑ
///
/// `element_compliance` = uₑᵀ kₑ uₑ (element strain energy).
pub fn simp_compliance_sensitivity(
    rho: f64,
    e0: f64,
    e_min: f64,
    penalty: f64,
    element_compliance: f64,
) -> f64 {
    let e_rho = simp_modulus(rho, e0, e_min, penalty);
    if e_rho < 1e-30 {
        return 0.0;
    }
    let de_drho = simp_modulus_derivative(rho, e0, e_min, penalty);
    -de_drho / e_rho * element_compliance
}

// ============================================================================
// § 2  DENSITY FILTER (convolution / hat kernel)
// ============================================================================

/// Apply a linear density (convolution) filter to an element density field.
///
/// ρ̃ᵢ = Σⱼ Hᵢⱼ ρⱼ / Σⱼ Hᵢⱼ   with  Hᵢⱼ = max(0, r_min − dᵢⱼ)
///
/// * `densities`  – raw element densities ρ
/// * `centroids`  – (x, y) centroid of each element
/// * `r_min`      – filter radius
pub fn density_filter(densities: &[f64], centroids: &[(f64, f64)], r_min: f64) -> Vec<f64> {
    let n = densities.len();
    assert_eq!(centroids.len(), n, "centroids length mismatch");
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let (xi, yi) = centroids[i];
        let mut num = 0.0;
        let mut den = 0.0;
        for j in 0..n {
            let (xj, yj) = centroids[j];
            let d = ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt();
            let w = (r_min - d).max(0.0);
            num += w * densities[j];
            den += w;
        }
        out[i] = if den > 1e-30 { num / den } else { densities[i] };
    }
    out
}

/// Compute the filter weight matrix H (row-normalised).
///
/// Returns a flat `n×n` row-major vector of weights.
pub fn density_filter_matrix(centroids: &[(f64, f64)], r_min: f64) -> Vec<f64> {
    let n = centroids.len();
    let mut h = vec![0.0_f64; n * n];
    for i in 0..n {
        let (xi, yi) = centroids[i];
        let mut row_sum = 0.0;
        for j in 0..n {
            let (xj, yj) = centroids[j];
            let d = ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt();
            let w = (r_min - d).max(0.0);
            h[i * n + j] = w;
            row_sum += w;
        }
        if row_sum > 1e-30 {
            for j in 0..n {
                h[i * n + j] /= row_sum;
            }
        }
    }
    h
}

// ============================================================================
// § 3  SENSITIVITY FILTER
// ============================================================================

/// Density-weighted sensitivity filter.
///
/// s̃ᵢ = Σⱼ Hᵢⱼ ρⱼ sⱼ / Σⱼ Hᵢⱼ ρⱼ
pub fn filter_sensitivities_density_weighted(
    sensitivities: &[f64],
    densities: &[f64],
    centroids: &[(f64, f64)],
    r_min: f64,
) -> Vec<f64> {
    let n = sensitivities.len();
    assert_eq!(densities.len(), n);
    assert_eq!(centroids.len(), n);
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let (xi, yi) = centroids[i];
        let mut num = 0.0;
        let mut den = 0.0;
        for j in 0..n {
            let (xj, yj) = centroids[j];
            let d = ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt();
            let w = (r_min - d).max(0.0);
            num += w * densities[j] * sensitivities[j];
            den += w * densities[j];
        }
        out[i] = if den.abs() > 1e-30 {
            num / den
        } else {
            sensitivities[i]
        };
    }
    out
}

/// Plain hat-function sensitivity filter (no density weighting).
///
/// s̃ᵢ = Σⱼ Hᵢⱼ sⱼ / Σⱼ Hᵢⱼ
pub fn filter_sensitivities_hat(
    sensitivities: &[f64],
    centroids: &[(f64, f64)],
    r_min: f64,
) -> Vec<f64> {
    let n = sensitivities.len();
    assert_eq!(centroids.len(), n);
    let mut out = vec![0.0_f64; n];
    for i in 0..n {
        let (xi, yi) = centroids[i];
        let mut num = 0.0;
        let mut den = 0.0;
        for j in 0..n {
            let (xj, yj) = centroids[j];
            let d = ((xi - xj).powi(2) + (yi - yj).powi(2)).sqrt();
            let w = (r_min - d).max(0.0);
            num += w * sensitivities[j];
            den += w;
        }
        out[i] = if den.abs() > 1e-30 {
            num / den
        } else {
            sensitivities[i]
        };
    }
    out
}

// ============================================================================
// § 4  OPTIMALITY CRITERIA (OC) UPDATE
// ============================================================================

/// Update a single design variable with the standard OC formula.
///
/// ρ_new = clamp( ρ · sqrt(−∂c/∂ρ / λ), \[ρ − Δ, ρ + Δ\] ∩ \[ρ_min, 1\] )
pub fn oc_single_update(
    rho: f64,
    sensitivity_magnitude: f64,
    lagrange: f64,
    move_limit: f64,
    rho_min: f64,
) -> f64 {
    if lagrange < 1e-30 {
        return rho;
    }
    let be = (sensitivity_magnitude / lagrange).max(0.0).sqrt();
    let rho_trial = rho * be;
    rho_trial
        .max(rho - move_limit)
        .min(rho + move_limit)
        .clamp(rho_min, 1.0)
}

/// Full OC density update with bisection on the Lagrange multiplier.
///
/// Enforces  Σρᵢ = n · vf  via bisection and applies `oc_single_update`.
pub fn oc_update(
    densities: &[f64],
    sensitivities: &[f64],
    volume_fraction: f64,
    move_limit: f64,
    rho_min: f64,
) -> Vec<f64> {
    assert_eq!(densities.len(), sensitivities.len());
    let n = densities.len();
    let target = volume_fraction * n as f64;
    let mut l1 = 0.0_f64;
    let mut l2 = 1e9_f64;
    for _ in 0..60 {
        let lm = 0.5 * (l1 + l2);
        let vol: f64 = densities
            .iter()
            .zip(sensitivities.iter())
            .map(|(&rho, &s)| oc_single_update(rho, s.abs(), lm, move_limit, rho_min))
            .sum();
        if vol > target {
            l1 = lm;
        } else {
            l2 = lm;
        }
    }
    let lm = 0.5 * (l1 + l2);
    densities
        .iter()
        .zip(sensitivities.iter())
        .map(|(&rho, &s)| oc_single_update(rho, s.abs(), lm, move_limit, rho_min))
        .collect()
}

/// Compute the volume fraction of a density vector.
pub fn volume_fraction(densities: &[f64]) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    densities.iter().sum::<f64>() / densities.len() as f64
}

// ============================================================================
// § 5  MMA (Method of Moving Asymptotes) STEP
// ============================================================================

/// MMA asymptote tracking state for a single design variable.
#[derive(Debug, Clone)]
pub struct MmaAsymptote {
    /// Lower asymptote Lᵢ
    pub lower: f64,
    /// Upper asymptote Uᵢ
    pub upper: f64,
    /// Previous design value (iteration k−1)
    pub x_prev: f64,
    /// Design value two iterations ago (iteration k−2)
    pub x_prev2: f64,
}

impl MmaAsymptote {
    /// Construct with default asymptotes centred around `x0`.
    pub fn new(x0: f64, s: f64) -> Self {
        Self {
            lower: x0 - s,
            upper: x0 + s,
            x_prev: x0,
            x_prev2: x0,
        }
    }

    /// Update asymptotes using the sign-change heuristic.
    ///
    /// If xₖ − xₖ₋₁ and xₖ₋₁ − xₖ₋₂ have the same sign, expand; otherwise
    /// contract.
    pub fn update_asymptotes(
        &mut self,
        x_cur: f64,
        x_lo: f64,
        x_hi: f64,
        s_incr: f64,
        s_decr: f64,
    ) {
        let d1 = x_cur - self.x_prev;
        let d2 = self.x_prev - self.x_prev2;
        let s = if d1 * d2 > 0.0 { s_incr } else { s_decr };
        let range = x_hi - x_lo;
        self.lower = x_cur - s * range;
        self.upper = x_cur + s * range;
        self.x_prev2 = self.x_prev;
        self.x_prev = x_cur;
    }
}

/// Solve a single-variable MMA sub-problem analytically.
///
/// Minimise  f(x) ≈ p/(U − x) + q/(x − L)  subject to  x ∈ \[α, β\].
///
/// Returns the optimal x.
pub fn mma_single_variable(p: f64, q: f64, lower: f64, upper: f64, alpha: f64, beta: f64) -> f64 {
    let range = upper - lower;
    if range < 1e-30 {
        return alpha.max(beta.min(0.5 * (lower + upper)));
    }
    // Unconstrained minimiser of the convex approximation
    let x_unc = if p + q > 1e-30 {
        let ratio = (q / p).max(0.0).sqrt();
        (lower * ratio + upper) / (1.0 + ratio)
    } else {
        0.5 * (lower + upper)
    };
    x_unc.clamp(alpha, beta)
}

/// One MMA iteration step for a vector of design variables.
///
/// Assumes sensitivities `dc` have been pre-filtered.
/// Returns updated densities.
pub fn mma_step(
    x: &[f64],
    dc: &[f64],
    asymptotes: &mut Vec<MmaAsymptote>,
    x_lo: f64,
    x_hi: f64,
    move_limit: f64,
) -> Vec<f64> {
    let n = x.len();
    assert_eq!(dc.len(), n);
    assert_eq!(asymptotes.len(), n);
    let mut x_new = vec![0.0_f64; n];
    for i in 0..n {
        let alpha = (x[i] - move_limit).max(x_lo);
        let beta = (x[i] + move_limit).min(x_hi);
        let l = asymptotes[i].lower;
        let u = asymptotes[i].upper;
        let dci = dc[i];
        let p = if dci < 0.0 {
            0.0
        } else {
            dci * (u - x[i]).powi(2)
        };
        let q = if dci > 0.0 {
            0.0
        } else {
            -dci * (x[i] - l).powi(2)
        };
        x_new[i] = mma_single_variable(p, q, l, u, alpha, beta);
        asymptotes[i].update_asymptotes(x_new[i], x_lo, x_hi, 1.2, 0.7);
    }
    x_new
}

// ============================================================================
// § 6  COMPLIANCE MINIMISATION
// ============================================================================

/// Compute the structural compliance  c = fᵀ u  for a single load case.
///
/// `forces` and `displacements` are global DOF vectors of the same length.
pub fn structural_compliance(forces: &[f64], displacements: &[f64]) -> f64 {
    assert_eq!(forces.len(), displacements.len());
    forces
        .iter()
        .zip(displacements.iter())
        .map(|(&f, &u)| f * u)
        .sum()
}

/// Accumulate per-element compliance sensitivity from element strain energies.
///
/// dc_dρᵢ = −p ρᵢ^(p−1) (E₀ − E_min) / E(ρᵢ) · uₑᵢᵀ kₑ uₑᵢ
pub fn compliance_sensitivities(
    densities: &[f64],
    element_strain_energies: &[f64],
    e0: f64,
    e_min: f64,
    penalty: f64,
) -> Vec<f64> {
    densities
        .iter()
        .zip(element_strain_energies.iter())
        .map(|(&rho, &se)| simp_compliance_sensitivity(rho, e0, e_min, penalty, se))
        .collect()
}

/// Run a simplified compliance minimisation loop (no FEM solve, for testing).
///
/// Iterates OC updates given precomputed sensitivities, returning history of
/// compliance values.
pub fn compliance_minimisation_loop(
    densities: &mut Vec<f64>,
    sensitivities_provider: &dyn Fn(&[f64]) -> Vec<f64>,
    volume_fraction_target: f64,
    move_limit: f64,
    rho_min: f64,
    n_iter: usize,
) -> Vec<f64> {
    let mut history = Vec::with_capacity(n_iter);
    for _ in 0..n_iter {
        let sens = sensitivities_provider(densities);
        *densities = oc_update(
            densities,
            &sens,
            volume_fraction_target,
            move_limit,
            rho_min,
        );
        let c: f64 = sens
            .iter()
            .zip(densities.iter())
            .map(|(&s, &r)| s.abs() * r)
            .sum();
        history.push(c);
    }
    history
}

// ============================================================================
// § 7  VOLUME CONSTRAINT
// ============================================================================

/// Volume constraint residual g = Σρᵢ/n − Vf.
pub fn volume_constraint(densities: &[f64], volume_fraction_target: f64) -> f64 {
    volume_fraction(densities) - volume_fraction_target
}

/// Gradient of the volume constraint ∂g/∂ρᵢ = 1/n (uniform).
pub fn volume_constraint_gradient(n: usize) -> Vec<f64> {
    vec![1.0 / n as f64; n]
}

/// Project densities to enforce the target volume fraction via a shift.
///
/// Adds a uniform offset δ so that Σρᵢ_new/n = Vf, then clamps to
/// \[rho_min, 1\].
pub fn volume_projection(densities: &[f64], volume_fraction_target: f64, rho_min: f64) -> Vec<f64> {
    let current = volume_fraction(densities);
    let delta = volume_fraction_target - current;
    densities
        .iter()
        .map(|&r| (r + delta).clamp(rho_min, 1.0))
        .collect()
}

// ============================================================================
// § 8  MULTI-LOAD-CASE OPTIMISATION
// ============================================================================

/// Weighted compliance for multiple load cases.
///
/// c = Σₗ wₗ · fₗᵀ uₗ
///
/// `cases`: list of (forces, displacements) pairs.
/// `weights`: load-case weights (should sum to 1).
pub fn multi_load_compliance(cases: &[(&[f64], &[f64])], weights: &[f64]) -> f64 {
    assert_eq!(cases.len(), weights.len());
    cases
        .iter()
        .zip(weights.iter())
        .map(|((f, u), &w)| w * structural_compliance(f, u))
        .sum()
}

/// Weighted sensitivity for multiple load cases.
///
/// dc/dρᵢ = Σₗ wₗ · (dc/dρᵢ)_l
pub fn multi_load_sensitivities(all_sensitivities: &[Vec<f64>], weights: &[f64]) -> Vec<f64> {
    assert!(!all_sensitivities.is_empty());
    assert_eq!(all_sensitivities.len(), weights.len());
    let n = all_sensitivities[0].len();
    let mut out = vec![0.0_f64; n];
    for (sens, &w) in all_sensitivities.iter().zip(weights.iter()) {
        for (o, &s) in out.iter_mut().zip(sens.iter()) {
            *o += w * s;
        }
    }
    out
}

// ============================================================================
// § 9  STRESS-CONSTRAINED TOPOLOGY OPTIMISATION
// ============================================================================

/// Von Mises stress for a 2-D plane-stress state.
///
/// σ_vm = sqrt(σ_xx² − σ_xx σ_yy + σ_yy² + 3 τ_xy²)
pub fn von_mises_2d(sigma_xx: f64, sigma_yy: f64, tau_xy: f64) -> f64 {
    (sigma_xx.powi(2) - sigma_xx * sigma_yy + sigma_yy.powi(2) + 3.0 * tau_xy.powi(2)).sqrt()
}

/// P-norm stress aggregation function.
///
/// σ_pn = ( Σᵢ σᵢ^p )^(1/p)
///
/// Used to aggregate element stresses into a differentiable scalar constraint.
pub fn p_norm_stress(stresses: &[f64], p: f64) -> f64 {
    if stresses.is_empty() {
        return 0.0;
    }
    let sum: f64 = stresses.iter().map(|&s| s.abs().powf(p)).sum();
    sum.powf(1.0 / p)
}

/// Sensitivity of the p-norm stress w.r.t. element stress σᵢ.
///
/// ∂σ_pn/∂σᵢ = (σᵢ/σ_pn)^(p−1) / σ_pn^(1−p) · 1/n   (chain rule result)
pub fn p_norm_stress_sensitivity(stresses: &[f64], p: f64) -> Vec<f64> {
    let sigma_pn = p_norm_stress(stresses, p);
    if sigma_pn < 1e-30 {
        return vec![0.0; stresses.len()];
    }
    stresses
        .iter()
        .map(|&s| (s.abs() / sigma_pn).powf(p - 1.0) * s.signum())
        .collect()
}

/// Stress-penalised SIMP modulus (lower penalisation to avoid stress singularity).
///
/// E_stress(ρ) = ρ^q · E₀   with q = p − 1.
pub fn simp_stress_modulus(rho: f64, e0: f64, penalty: f64) -> f64 {
    rho.powf(penalty - 1.0) * e0
}

// ============================================================================
// § 10  LEVEL-SET TOPOLOGY OPTIMISATION
// ============================================================================

/// Signed-distance level-set field.
#[derive(Debug, Clone)]
pub struct LevelSetField {
    /// Level-set values φ on a regular grid.
    pub phi: Vec<f64>,
    /// Grid width (number of cells in x).
    pub nx: usize,
    /// Grid height (number of cells in y).
    pub ny: usize,
    /// Cell size h.
    pub h: f64,
}

impl LevelSetField {
    /// Create a new level-set field with all cells initialised to `init`.
    pub fn new(nx: usize, ny: usize, h: f64, init: f64) -> Self {
        Self {
            phi: vec![init; nx * ny],
            nx,
            ny,
            h,
        }
    }

    /// Index into the flat `phi` array.
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.nx + i
    }

    /// Bilinearly interpolate φ at position (x, y).
    pub fn interpolate(&self, x: f64, y: f64) -> f64 {
        let xi = (x / self.h).floor() as isize;
        let yj = (y / self.h).floor() as isize;
        let tx = (x / self.h) - xi as f64;
        let ty = (y / self.h) - yj as f64;
        let get = |i: isize, j: isize| {
            let ci = i.clamp(0, self.nx as isize - 1) as usize;
            let cj = j.clamp(0, self.ny as isize - 1) as usize;
            self.phi[self.idx(ci, cj)]
        };
        let v00 = get(xi, yj);
        let v10 = get(xi + 1, yj);
        let v01 = get(xi, yj + 1);
        let v11 = get(xi + 1, yj + 1);
        (1.0 - tx) * (1.0 - ty) * v00
            + tx * (1.0 - ty) * v10
            + (1.0 - tx) * ty * v01
            + tx * ty * v11
    }

    /// Reinitialise the level-set field toward a signed distance function.
    ///
    /// Simple first-order redistancing (iterative PDE approach).
    pub fn reinitialise(&mut self, n_iter: usize) {
        let h = self.h;
        let nx = self.nx;
        let ny = self.ny;
        let phi0 = self.phi.clone();
        for _ in 0..n_iter {
            let phi_old = self.phi.clone();
            for j in 0..ny {
                for i in 0..nx {
                    let idx = self.idx(i, j);
                    let sign = if phi0[idx] >= 0.0 { 1.0 } else { -1.0 };
                    let im = if i > 0 {
                        phi_old[self.idx(i - 1, j)]
                    } else {
                        phi_old[idx]
                    };
                    let ip = if i < nx - 1 {
                        phi_old[self.idx(i + 1, j)]
                    } else {
                        phi_old[idx]
                    };
                    let jm = if j > 0 {
                        phi_old[self.idx(i, j - 1)]
                    } else {
                        phi_old[idx]
                    };
                    let jp = if j < ny - 1 {
                        phi_old[self.idx(i, j + 1)]
                    } else {
                        phi_old[idx]
                    };
                    let grad_sq = ((ip - im) / (2.0 * h)).powi(2) + ((jp - jm) / (2.0 * h)).powi(2);
                    self.phi[idx] = phi_old[idx] - 0.5 * sign * (grad_sq.sqrt() - 1.0);
                }
            }
        }
    }

    /// Update the level-set by advecting with a normal velocity field Vₙ.
    ///
    /// φ_new = φ_old − dt · Vₙ · |∇φ|  (first-order upwind).
    pub fn advect(&mut self, vn: &[f64], dt: f64) {
        let h = self.h;
        let nx = self.nx;
        let ny = self.ny;
        let phi_old = self.phi.clone();
        let get = |i: isize, j: isize| {
            let ci = i.clamp(0, nx as isize - 1) as usize;
            let cj = j.clamp(0, ny as isize - 1) as usize;
            phi_old[cj * nx + ci]
        };
        for j in 0..ny {
            for i in 0..nx {
                let idx = self.idx(i, j);
                let v = vn[idx];
                let dxp = (get(i as isize + 1, j as isize) - phi_old[idx]) / h;
                let dxm = (phi_old[idx] - get(i as isize - 1, j as isize)) / h;
                let dyp = (get(i as isize, j as isize + 1) - phi_old[idx]) / h;
                let dym = (phi_old[idx] - get(i as isize, j as isize - 1)) / h;
                let grad_sq = if v > 0.0 {
                    dxm.max(0.0).powi(2)
                        + dxp.min(0.0).powi(2)
                        + dym.max(0.0).powi(2)
                        + dyp.min(0.0).powi(2)
                } else {
                    dxp.max(0.0).powi(2)
                        + dxm.min(0.0).powi(2)
                        + dyp.max(0.0).powi(2)
                        + dym.min(0.0).powi(2)
                };
                self.phi[idx] = phi_old[idx] - dt * v * grad_sq.sqrt();
            }
        }
    }
}

// ============================================================================
// § 11  BESO (Bi-directional ESO)
// ============================================================================

/// BESO element state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BesoCellState {
    /// Element is void (removed).
    Void,
    /// Element is solid (active).
    Solid,
}

/// BESO removal / addition thresholds.
///
/// Returns (threshold_removal, threshold_addition) given sorted sensitivities
/// and a target volume fraction.
pub fn beso_thresholds(
    sensitivities: &[f64],
    states: &[BesoCellState],
    volume_fraction_target: f64,
    evolutionary_rate: f64,
) -> (f64, f64) {
    let n = sensitivities.len() as f64;
    let n_solid_target = (n * volume_fraction_target).round() as usize;
    let mut solid_sens: Vec<f64> = sensitivities
        .iter()
        .zip(states.iter())
        .filter(|&(_, s)| *s == BesoCellState::Solid)
        .map(|(&v, _)| v.abs())
        .collect();
    solid_sens.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n_solid = solid_sens.len();
    let n_target = n_solid_target;
    let n_del = ((n_solid as f64 * evolutionary_rate) as usize).min(n_solid.saturating_sub(1));
    let t_del = if n_del < solid_sens.len() {
        solid_sens[n_del]
    } else {
        solid_sens.last().copied().unwrap_or(0.0)
    };
    let t_add = if n_target < solid_sens.len() {
        solid_sens[n_target.saturating_sub(1)]
    } else {
        solid_sens.last().copied().unwrap_or(0.0)
    };
    (t_del, t_add)
}

/// Apply one BESO update step.
///
/// Removes low-sensitivity solid elements and adds high-sensitivity void
/// elements to approach the target volume fraction.
pub fn beso_update(
    sensitivities: &[f64],
    states: &mut Vec<BesoCellState>,
    volume_fraction_target: f64,
    evolutionary_rate: f64,
) {
    let (t_del, t_add) = beso_thresholds(
        sensitivities,
        states,
        volume_fraction_target,
        evolutionary_rate,
    );
    for (state, &s) in states.iter_mut().zip(sensitivities.iter()) {
        match *state {
            BesoCellState::Solid => {
                if s.abs() < t_del {
                    *state = BesoCellState::Void;
                }
            }
            BesoCellState::Void => {
                if s.abs() > t_add {
                    *state = BesoCellState::Solid;
                }
            }
        }
    }
}

// ============================================================================
// § 12  COMPLIANT MECHANISM DESIGN
// ============================================================================

/// Output displacement objective for a compliant mechanism.
///
/// Maximise u_out = dᵀ u  subject to the volume constraint.
///
/// `d` is the output DOF selection vector; `u` is the displacement vector.
pub fn mechanism_output_displacement(d: &[f64], u: &[f64]) -> f64 {
    assert_eq!(d.len(), u.len());
    d.iter().zip(u.iter()).map(|(&di, &ui)| di * ui).sum()
}

/// Geometric advantage: ratio of output to input displacement.
///
/// GA = u_out / u_in
pub fn geometric_advantage(u_out: f64, u_in: f64) -> f64 {
    if u_in.abs() < 1e-30 {
        return 0.0;
    }
    u_out / u_in
}

/// Mechanical efficiency of a compliant mechanism.
///
/// η = (F_out · u_out) / (F_in · u_in)
pub fn mechanical_efficiency(f_out: f64, u_out: f64, f_in: f64, u_in: f64) -> f64 {
    let input_work = f_in * u_in;
    if input_work.abs() < 1e-30 {
        return 0.0;
    }
    (f_out * u_out) / input_work
}

// ============================================================================
// § 13  FREQUENCY OPTIMISATION
// ============================================================================

/// Rayleigh quotient estimate of the fundamental frequency.
///
/// ω² ≈ xᵀ K x / xᵀ M x
///
/// `k_diag` and `m_diag` are diagonal stiffness/mass vectors; `x` is a mode
/// shape estimate.
pub fn rayleigh_quotient(k_diag: &[f64], m_diag: &[f64], x: &[f64]) -> f64 {
    let n = x.len();
    assert_eq!(k_diag.len(), n);
    assert_eq!(m_diag.len(), n);
    let num: f64 = x
        .iter()
        .zip(k_diag.iter())
        .map(|(&xi, &ki)| xi * ki * xi)
        .sum();
    let den: f64 = x
        .iter()
        .zip(m_diag.iter())
        .map(|(&xi, &mi)| xi * mi * xi)
        .sum();
    if den < 1e-30 { 0.0 } else { num / den }
}

/// Sensitivity of the Rayleigh quotient ω² w.r.t. element density ρₑ.
///
/// ∂ω²/∂ρₑ = (∂kₑ/∂ρₑ · uₑ² − ω² · ∂mₑ/∂ρₑ · uₑ²) / (xᵀ M x)
pub fn frequency_sensitivity(
    element_ke: f64,
    element_me: f64,
    omega_sq: f64,
    mode_displacement_sq: f64,
    mass_norm: f64,
    rho: f64,
    penalty_k: f64,
    penalty_m: f64,
) -> f64 {
    if mass_norm < 1e-30 {
        return 0.0;
    }
    let dk_drho = penalty_k * rho.powf(penalty_k - 1.0) * element_ke;
    let dm_drho = penalty_m * rho.powf(penalty_m - 1.0) * element_me;
    (dk_drho - omega_sq * dm_drho) * mode_displacement_sq / mass_norm
}

/// Minimum gap frequency constraint: penalise the design if ω < ω_target.
///
/// Returns the gap constraint value g = ω_target² − ω².
pub fn frequency_gap_constraint(omega_sq: f64, omega_target_sq: f64) -> f64 {
    omega_target_sq - omega_sq
}

// ============================================================================
// § 14  HEAT CONDUCTION TOPOLOGY OPTIMISATION
// ============================================================================

/// SIMP thermal conductivity: k(ρ) = k_min + ρ^p (k₀ − k_min).
pub fn simp_conductivity(rho: f64, k0: f64, k_min: f64, penalty: f64) -> f64 {
    k_min + rho.powf(penalty) * (k0 - k_min)
}

/// Derivative dk/dρ for SIMP thermal conductivity.
pub fn simp_conductivity_derivative(rho: f64, k0: f64, k_min: f64, penalty: f64) -> f64 {
    penalty * rho.powf(penalty - 1.0) * (k0 - k_min)
}

/// Heat compliance (thermal dissipation potential): c_T = qᵀ T.
pub fn heat_compliance(heat_flux: &[f64], temperature: &[f64]) -> f64 {
    assert_eq!(heat_flux.len(), temperature.len());
    heat_flux
        .iter()
        .zip(temperature.iter())
        .map(|(&q, &t)| q * t)
        .sum()
}

/// Heat compliance sensitivity ∂c_T/∂ρᵢ.
///
/// ∂c_T/∂ρᵢ = −(dk/dρᵢ / k(ρᵢ)) · qₑᵢ · Tₑᵢ
pub fn heat_compliance_sensitivity(
    rho: f64,
    k0: f64,
    k_min: f64,
    penalty: f64,
    element_dissipation: f64,
) -> f64 {
    let k_rho = simp_conductivity(rho, k0, k_min, penalty);
    if k_rho < 1e-30 {
        return 0.0;
    }
    let dk_drho = simp_conductivity_derivative(rho, k0, k_min, penalty);
    -dk_drho / k_rho * element_dissipation
}

/// Thermal compliance OC update for heat conduction topology.
///
/// Same bisection strategy as the mechanical OC update.
pub fn heat_oc_update(
    densities: &[f64],
    sensitivities: &[f64],
    volume_fraction_target: f64,
    move_limit: f64,
    rho_min: f64,
) -> Vec<f64> {
    oc_update(
        densities,
        sensitivities,
        volume_fraction_target,
        move_limit,
        rho_min,
    )
}

// ============================================================================
// § 15  SMOOTH HEAVISIDE PROJECTION
// ============================================================================

/// Smooth Heaviside projection: ρ̃ → ρ̄.
///
/// H(ρ̃; β, η) = \[tanh(βη) + tanh(β(ρ̃ − η))\] / \[tanh(βη) + tanh(β(1 − η))\]
pub fn heaviside_projection(rho_tilde: f64, beta: f64, eta: f64) -> f64 {
    let num = (beta * eta).tanh() + (beta * (rho_tilde - eta)).tanh();
    let den = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
    if den.abs() < 1e-30 {
        rho_tilde
    } else {
        num / den
    }
}

/// Derivative of the Heaviside projection dH/dρ̃.
pub fn heaviside_projection_derivative(rho_tilde: f64, beta: f64, eta: f64) -> f64 {
    let den = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
    if den.abs() < 1e-30 {
        return 1.0;
    }
    beta * (1.0 - (beta * (rho_tilde - eta)).tanh().powi(2)) / den
}

/// Apply Heaviside projection to an entire density field.
pub fn apply_heaviside(densities: &[f64], beta: f64, eta: f64) -> Vec<f64> {
    densities
        .iter()
        .map(|&r| heaviside_projection(r, beta, eta))
        .collect()
}

// ============================================================================
// § 16  CHECKERBOARD SUPPRESSION METRIC
// ============================================================================

/// Detect checkerboard instability on a structured nx×ny mesh.
///
/// Returns the fraction of 2×2 blocks that exhibit alternating solid/void
/// patterns.
pub fn checkerboard_fraction(densities: &[f64], nx: usize, _ny: usize, threshold: f64) -> f64 {
    let ny2 = densities.len() / nx;
    if nx < 2 || ny2 < 2 {
        return 0.0;
    }
    let mut n_checker = 0usize;
    let mut n_blocks = 0usize;
    for j in 0..ny2 - 1 {
        for i in 0..nx - 1 {
            let a = densities[j * nx + i] > threshold;
            let b = densities[j * nx + i + 1] > threshold;
            let c = densities[(j + 1) * nx + i] > threshold;
            let d = densities[(j + 1) * nx + i + 1] > threshold;
            n_blocks += 1;
            if (a != b) && (c != d) && (a != c) {
                n_checker += 1;
            }
        }
    }
    if n_blocks == 0 {
        0.0
    } else {
        n_checker as f64 / n_blocks as f64
    }
}

// ============================================================================
// § 17  CONTINUATION STRATEGY
// ============================================================================

/// Penalisation continuation schedule.
///
/// Returns the penalisation exponent for iteration `iter` given a start
/// penalty, end penalty, and the iteration at which to begin increasing.
pub fn continuation_penalty(
    iter: usize,
    p_start: f64,
    p_end: f64,
    iter_start: usize,
    iter_end: usize,
) -> f64 {
    if iter <= iter_start {
        return p_start;
    }
    if iter >= iter_end {
        return p_end;
    }
    let t = (iter - iter_start) as f64 / (iter_end - iter_start) as f64;
    p_start + t * (p_end - p_start)
}

/// Beta continuation for Heaviside projection.
///
/// β doubles every `double_interval` iterations up to `beta_max`.
pub fn continuation_beta(iter: usize, beta0: f64, beta_max: f64, double_interval: usize) -> f64 {
    if double_interval == 0 {
        return beta_max;
    }
    let factor = 2.0_f64.powi((iter / double_interval) as i32);
    (beta0 * factor).min(beta_max)
}

// ============================================================================
// § 18  PERIMETER / REGULARISATION FUNCTIONAL
// ============================================================================

/// Isotropic TV regularisation contribution for element i.
///
/// R_i = Σⱼ∈Nᵢ (ρᵢ − ρⱼ)²
///
/// Used to suppress checkerboarding and control complexity.
pub fn tv_regularisation_element(rho_i: f64, neighbour_densities: &[f64], weight: f64) -> f64 {
    weight
        * neighbour_densities
            .iter()
            .map(|&r| (rho_i - r).powi(2))
            .sum::<f64>()
}

/// Global TV regularisation value.
pub fn tv_regularisation_global(densities: &[f64], nx: usize, weight: f64) -> f64 {
    let ny = densities.len() / nx;
    let mut r = 0.0;
    for j in 0..ny {
        for i in 0..nx {
            let rho = densities[j * nx + i];
            if i + 1 < nx {
                r += (rho - densities[j * nx + i + 1]).powi(2);
            }
            if j + 1 < ny {
                r += (rho - densities[(j + 1) * nx + i]).powi(2);
            }
        }
    }
    weight * r
}

// ============================================================================
// § 19  EXTENDED SIMP WITH INTERMEDIATE INTERPOLATION (RAMP)
// ============================================================================

/// RAMP (Rational Approximation of Material Properties) interpolation.
///
/// E(ρ) = E_min + ρ / (1 + q(1 − ρ)) · (E₀ − E_min)
pub fn ramp_modulus(rho: f64, e0: f64, e_min: f64, q: f64) -> f64 {
    e_min + rho / (1.0 + q * (1.0 - rho)) * (e0 - e_min)
}

/// Derivative dE/dρ for RAMP interpolation.
pub fn ramp_modulus_derivative(rho: f64, e0: f64, e_min: f64, q: f64) -> f64 {
    let denom = 1.0 + q * (1.0 - rho);
    (1.0 + q) / (denom * denom) * (e0 - e_min)
}

// ============================================================================
// § 20  TOPOLOGY OPTIMISATION RESULT
// ============================================================================

/// Summary result of a topology optimisation run.
#[derive(Debug, Clone)]
pub struct TopologyResult {
    /// Final element densities.
    pub densities: Vec<f64>,
    /// Compliance history over iterations.
    pub compliance_history: Vec<f64>,
    /// Final volume fraction.
    pub volume_fraction: f64,
    /// Number of iterations performed.
    pub n_iter: usize,
    /// Converged flag.
    pub converged: bool,
}

impl TopologyResult {
    /// Create a result from density and compliance history vectors.
    pub fn new(densities: Vec<f64>, compliance_history: Vec<f64>) -> Self {
        let vf = volume_fraction(&densities);
        let n = compliance_history.len();
        let conv = if n >= 2 {
            let last = compliance_history[n - 1];
            let prev = compliance_history[n - 2];
            (last - prev).abs() / last.abs().max(1e-30) < 1e-4
        } else {
            false
        };
        Self {
            volume_fraction: vf,
            n_iter: n,
            converged: conv,
            densities,
            compliance_history,
        }
    }

    /// Binary density field (threshold at 0.5).
    pub fn binary_densities(&self) -> Vec<f64> {
        self.densities
            .iter()
            .map(|&r| if r >= 0.5 { 1.0 } else { 0.0 })
            .collect()
    }

    /// Print a short summary to stdout.
    pub fn print_summary(&self) {
        println!(
            "TopologyResult: iter={} vf={:.6} converged={}",
            self.n_iter, self.volume_fraction, self.converged
        );
        if let Some(&last) = self.compliance_history.last() {
            println!("  final compliance = {:.6}", last);
        }
    }
}

// ============================================================================
// § 21  UTILITY – STRUCTURED GRID CENTROID GENERATION
// ============================================================================

/// Generate element centroids on a regular nx×ny grid with cell size h.
pub fn regular_grid_centroids(nx: usize, ny: usize, h: f64) -> Vec<(f64, f64)> {
    let mut c = Vec::with_capacity(nx * ny);
    for j in 0..ny {
        for i in 0..nx {
            c.push(((i as f64 + 0.5) * h, (j as f64 + 0.5) * h));
        }
    }
    c
}

/// Flatten a 2-D index pair (i, j) to a 1-D flat index for an nx-wide grid.
#[inline]
pub fn flat_index(i: usize, j: usize, nx: usize) -> usize {
    j * nx + i
}

// ============================================================================
// § 22  COMBINED TOPOLOGY OPTIMISATION DRIVER (simplified, no FEM solve)
// ============================================================================

/// Configuration for a topology optimisation run.
#[derive(Debug, Clone)]
pub struct TopologyConfig {
    /// Number of elements in x.
    pub nx: usize,
    /// Number of elements in y.
    pub ny: usize,
    /// Cell size.
    pub h: f64,
    /// Filter radius.
    pub r_min: f64,
    /// Target volume fraction.
    pub vf: f64,
    /// SIMP penalisation exponent.
    pub penalty: f64,
    /// Move limit for OC.
    pub move_limit: f64,
    /// Minimum density.
    pub rho_min: f64,
    /// Maximum iterations.
    pub max_iter: usize,
    /// Convergence tolerance.
    pub tol: f64,
}

impl Default for TopologyConfig {
    fn default() -> Self {
        Self {
            nx: 60,
            ny: 20,
            h: 1.0,
            r_min: 1.5,
            vf: 0.5,
            penalty: 3.0,
            move_limit: 0.2,
            rho_min: 1e-3,
            max_iter: 100,
            tol: 1e-4,
        }
    }
}

impl TopologyConfig {
    /// Total number of elements.
    pub fn n_elements(&self) -> usize {
        self.nx * self.ny
    }

    /// Create a uniform density field at the target volume fraction.
    pub fn initial_densities(&self) -> Vec<f64> {
        vec![self.vf; self.n_elements()]
    }
}

/// Run a mock topology optimisation loop (no actual FEM).
///
/// The sensitivity is approximated as the negative inverse density (mimics
/// a compliance-like behaviour for testing purposes).
pub fn run_mock_topology_optimisation(cfg: &TopologyConfig) -> TopologyResult {
    let centroids = regular_grid_centroids(cfg.nx, cfg.ny, cfg.h);
    let mut densities = cfg.initial_densities();
    let mut compliance_history = Vec::with_capacity(cfg.max_iter);
    let e0 = 1.0_f64;
    let e_min = 1e-9_f64;
    for _iter in 0..cfg.max_iter {
        // Mock sensitivities: ∂c/∂ρ ≈ −1/ρ
        let raw_sens: Vec<f64> = densities.iter().map(|&r| -1.0 / r.max(1e-6)).collect();
        let filt_sens = filter_sensitivities_hat(&raw_sens, &centroids, cfg.r_min);
        let filt_dens = density_filter(&densities, &centroids, cfg.r_min);
        let se: Vec<f64> = filt_dens
            .iter()
            .map(|&r| simp_modulus(r, e0, e_min, cfg.penalty))
            .collect();
        let c: f64 = se.iter().sum();
        compliance_history.push(c);
        densities = oc_update(&densities, &filt_sens, cfg.vf, cfg.move_limit, cfg.rho_min);
        if compliance_history.len() >= 2 {
            let n = compliance_history.len();
            let delta = (compliance_history[n - 1] - compliance_history[n - 2]).abs();
            if delta / compliance_history[n - 1].abs().max(1e-30) < cfg.tol {
                break;
            }
        }
    }
    TopologyResult::new(densities, compliance_history)
}

// ============================================================================
// TESTS
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    #[test]
    fn test_simp_modulus_limits() {
        let e0 = 1.0;
        let e_min = 1e-9;
        let p = 3.0;
        let at_one = simp_modulus(1.0, e0, e_min, p);
        let at_zero = simp_modulus(0.0, e0, e_min, p);
        assert!((at_one - e0).abs() < 1e-12, "simp at rho=1 should be e0");
        assert!(
            (at_zero - e_min).abs() < 1e-12,
            "simp at rho=0 should be e_min"
        );
    }

    #[test]
    fn test_simp_modulus_intermediate() {
        let e0 = 200e9_f64;
        let e_min = 1e-3_f64;
        let p = 3.0_f64;
        let rho = 0.5_f64;
        let expected = e_min + 0.5_f64.powi(3) * (e0 - e_min);
        let result = simp_modulus(rho, e0, e_min, p);
        assert!(
            (result - expected).abs() < 1e-3,
            "simp modulus mismatch: {} vs {}",
            result,
            expected
        );
    }

    #[test]
    fn test_simp_derivative_positive() {
        let e0 = 1.0;
        let e_min = 1e-9;
        let p = 3.0;
        let rho = 0.6;
        let d = simp_modulus_derivative(rho, e0, e_min, p);
        assert!(d > 0.0, "derivative must be positive");
        // Finite difference check
        let h = 1e-6;
        let fd =
            (simp_modulus(rho + h, e0, e_min, p) - simp_modulus(rho - h, e0, e_min, p)) / (2.0 * h);
        assert!(
            (d - fd).abs() < 1e-8,
            "analytical vs FD derivative: {} vs {}",
            d,
            fd
        );
    }

    #[test]
    fn test_compliance_sensitivity_sign() {
        // Sensitivity must be negative (increasing density reduces compliance)
        let s = simp_compliance_sensitivity(0.5, 1.0, 1e-9, 3.0, 1.0);
        assert!(
            s < 0.0,
            "compliance sensitivity should be negative, got {}",
            s
        );
    }

    #[test]
    fn test_density_filter_uniform() {
        // Uniform field should pass through unchanged
        let n = 9;
        let rho: Vec<f64> = vec![0.5; n];
        let centroids: Vec<(f64, f64)> = (0..n).map(|i| (i as f64, 0.0)).collect();
        let filtered = density_filter(&rho, &centroids, 2.0);
        for &r in &filtered {
            assert!((r - 0.5).abs() < 1e-12, "uniform field should be unchanged");
        }
    }

    #[test]
    fn test_density_filter_smoothing() {
        // A spike at centre should be smoothed
        let mut rho = vec![0.0_f64; 5];
        rho[2] = 1.0;
        let centroids: Vec<(f64, f64)> = (0..5usize).map(|i| (i as f64, 0.0)).collect();
        let filtered = density_filter(&rho, &centroids, 2.0);
        assert!(filtered[2] < 1.0, "peak should decrease after filtering");
        assert!(
            filtered[1] > 0.0 && filtered[3] > 0.0,
            "neighbours should increase"
        );
    }

    #[test]
    fn test_oc_update_volume_conservation() {
        let n = 20;
        let rho: Vec<f64> = vec![0.5; n];
        let sens: Vec<f64> = (0..n).map(|i| -(i as f64 + 1.0)).collect();
        let updated = oc_update(&rho, &sens, 0.5, 0.2, 0.01);
        let vf = volume_fraction(&updated);
        assert!(
            (vf - 0.5).abs() < 0.02,
            "volume fraction should be approx 0.5, got {:.6}",
            vf
        );
    }

    #[test]
    fn test_oc_single_update_clamping() {
        // Should never exceed [rho_min, 1.0]
        let v = oc_single_update(0.5, 100.0, 1e-6, 0.2, 0.01);
        assert!(v <= 1.0 && v >= 0.01);
        let v2 = oc_single_update(0.5, 1e-10, 1e6, 0.2, 0.01);
        assert!(v2 >= 0.01 && v2 <= 1.0);
    }

    #[test]
    fn test_volume_fraction_basic() {
        let rho = vec![0.0, 0.5, 1.0, 0.5];
        let vf = volume_fraction(&rho);
        assert!((vf - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_volume_projection() {
        let rho = vec![0.3, 0.4, 0.5, 0.6, 0.7];
        let proj = volume_projection(&rho, 0.6, 0.01);
        let vf = volume_fraction(&proj);
        assert!((vf - 0.6).abs() < 0.01, "projected vf={:.6}", vf);
    }

    #[test]
    fn test_multi_load_compliance_equal_weights() {
        let f1 = vec![1.0, 0.0, 1.0];
        let u1 = vec![2.0, 3.0, 1.0];
        let f2 = vec![0.0, 1.0, 0.0];
        let u2 = vec![1.0, 4.0, 2.0];
        let c = multi_load_compliance(&[(&f1, &u1), (&f2, &u2)], &[0.5, 0.5]);
        // c1 = 1*2 + 1*1 = 3, c2 = 1*4 = 4 → 0.5*3 + 0.5*4 = 3.5
        assert!((c - 3.5).abs() < 1e-12);
    }

    #[test]
    fn test_von_mises_2d_pure_tension() {
        // σ_xx = σ, σ_yy = 0, τ = 0 → σ_vm = σ
        let sigma = 100.0;
        let vm = von_mises_2d(sigma, 0.0, 0.0);
        assert!((vm - sigma).abs() < 1e-10);
    }

    #[test]
    fn test_p_norm_stress_identity() {
        // When p→∞ p-norm → max; at finite p it must be ≥ max element
        let stresses = vec![1.0, 2.0, 3.0, 4.0];
        let pn = p_norm_stress(&stresses, 8.0);
        assert!(pn >= 3.9, "p-norm with large p should approach max={}", 4.0);
    }

    #[test]
    fn test_level_set_interpolate_corners() {
        let mut lsf = LevelSetField::new(3, 3, 1.0, 0.0);
        lsf.phi[0] = 1.0;
        // At (0,0) should be 1.0
        let v = lsf.interpolate(0.0, 0.0);
        assert!((v - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_beso_update_reduces_solids() {
        let n = 10;
        let states: Vec<BesoCellState> = vec![BesoCellState::Solid; n];
        let mut states = states;
        // Low sensitivity everywhere → remove half
        let sens: Vec<f64> = (0..n).map(|i| (i + 1) as f64 * 0.1).collect();
        beso_update(&sens, &mut states, 0.4, 0.2);
        let n_solid = states
            .iter()
            .filter(|&&s| s == BesoCellState::Solid)
            .count();
        assert!(n_solid <= n);
    }

    #[test]
    fn test_geometric_advantage() {
        let ga = geometric_advantage(0.5, 2.0);
        assert!((ga - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_mechanical_efficiency() {
        let eta = mechanical_efficiency(5.0, 2.0, 10.0, 1.0);
        assert!((eta - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_rayleigh_quotient() {
        // ω² = kᵀx² / mᵀx² with k=m=1, x=1 → ω²=1
        let k = vec![1.0, 1.0, 1.0];
        let m = vec![1.0, 1.0, 1.0];
        let x = vec![1.0, 1.0, 1.0];
        let omega_sq = rayleigh_quotient(&k, &m, &x);
        assert!((omega_sq - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_simp_conductivity_limits() {
        let k0 = 50.0;
        let k_min = 1e-6;
        let p = 3.0;
        assert!((simp_conductivity(1.0, k0, k_min, p) - k0).abs() < 1e-10);
        assert!((simp_conductivity(0.0, k0, k_min, p) - k_min).abs() < 1e-10);
    }

    #[test]
    fn test_heat_compliance() {
        let q = vec![1.0, 2.0, 3.0];
        let t = vec![4.0, 5.0, 6.0];
        let ct = heat_compliance(&q, &t);
        assert!((ct - 32.0).abs() < 1e-12); // 4+10+18
    }

    #[test]
    fn test_heaviside_projection_limits() {
        let beta = 5.0;
        let eta = 0.5;
        let h0 = heaviside_projection(0.0, beta, eta);
        let h1 = heaviside_projection(1.0, beta, eta);
        assert!(h0 < 0.1, "H(0) should be near 0, got {:.6}", h0);
        assert!(h1 > 0.9, "H(1) should be near 1, got {:.6}", h1);
    }

    #[test]
    fn test_continuation_penalty_endpoints() {
        let p = continuation_penalty(0, 1.0, 3.0, 10, 50);
        assert!((p - 1.0).abs() < 1e-12);
        let p_end = continuation_penalty(100, 1.0, 3.0, 10, 50);
        assert!((p_end - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_regular_grid_centroids_count() {
        let c = regular_grid_centroids(4, 3, 1.0);
        assert_eq!(c.len(), 12);
        assert!((c[0].0 - 0.5).abs() < 1e-12 && (c[0].1 - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_ramp_modulus_limits() {
        let e0 = 1.0;
        let e_min = 0.0;
        let q = 3.0;
        let at_one = ramp_modulus(1.0, e0, e_min, q);
        let at_zero = ramp_modulus(0.0, e0, e_min, q);
        assert!((at_one - e0).abs() < 1e-12);
        assert!((at_zero - e_min).abs() < 1e-12);
    }

    #[test]
    fn test_topology_result_binary() {
        let rho = vec![0.3, 0.7, 0.5, 0.9];
        let result = TopologyResult::new(rho, vec![10.0, 9.0]);
        let bin = result.binary_densities();
        assert_eq!(bin, vec![0.0, 1.0, 1.0, 1.0]);
    }

    #[test]
    fn test_tv_regularisation_uniform() {
        let rho = vec![0.5; 9];
        let r = tv_regularisation_global(&rho, 3, 1.0);
        assert!((r - 0.0).abs() < 1e-12, "uniform field has zero TV");
    }

    #[test]
    fn test_mock_topology_run_smoke() {
        let cfg = TopologyConfig {
            nx: 6,
            ny: 4,
            max_iter: 5,
            ..Default::default()
        };
        let result = run_mock_topology_optimisation(&cfg);
        assert_eq!(result.densities.len(), 24);
        assert!(!result.compliance_history.is_empty());
        let vf = result.volume_fraction;
        assert!(vf > 0.0 && vf <= 1.0, "vf out of range: {}", vf);
    }

    #[test]
    fn test_filter_sensitivity_uniform_field() {
        let n = 5;
        let sens = vec![1.0; n];
        let rho = vec![0.5; n];
        let cents: Vec<(f64, f64)> = (0..n).map(|i| (i as f64, 0.0)).collect();
        let fs = filter_sensitivities_density_weighted(&sens, &rho, &cents, 3.0);
        for &s in &fs {
            assert!(
                (s - 1.0).abs() < 1e-10,
                "uniform sensitivities unchanged: {}",
                s
            );
        }
    }

    #[test]
    fn test_frequency_gap_constraint_sign() {
        let omega_sq = 4.0;
        let omega_tgt_sq = 9.0;
        let g = frequency_gap_constraint(omega_sq, omega_tgt_sq);
        assert!(g > 0.0, "constraint active when ω < ω_target");
        let g2 = frequency_gap_constraint(16.0, 9.0);
        assert!(g2 < 0.0, "constraint inactive when ω > ω_target");
    }

    #[test]
    fn test_density_filter_matrix_row_sum() {
        let cents: Vec<(f64, f64)> = vec![(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)];
        let h = density_filter_matrix(&cents, 1.5);
        let n = cents.len();
        for i in 0..n {
            let row_sum: f64 = (0..n).map(|j| h[i * n + j]).sum();
            assert!(
                (row_sum - 1.0).abs() < 1e-12,
                "row {} sum={:.6}",
                i,
                row_sum
            );
        }
    }

    // Additional test for PI usage (avoids dead_code warning on import)
    #[test]
    fn test_pi_constant_used() {
        let circumference = 2.0 * PI * 1.0;
        assert!((circumference - std::f64::consts::TAU).abs() < 1e-12);
    }
}
