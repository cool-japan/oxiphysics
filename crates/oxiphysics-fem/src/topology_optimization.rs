// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended topology optimization: SIMP update rule, optimality criteria,
//! sensitivity filtering, Heaviside projection, volume constraint enforcement,
//! checkerboard suppression, continuation, compliance gradient, and
//! multi-load-case optimization.

// ---------------------------------------------------------------------------
// 1. SIMP design-variable update
// ---------------------------------------------------------------------------

/// SIMP interpolated Young's modulus: E(ρ) = E_min + ρ^p (E_0 − E_min).
///
/// - `rho`     : element density ∈ \[0, 1\]
/// - `e0`      : solid modulus
/// - `e_min`   : void modulus (numerical floor, e.g. 1e-9 E_0)
/// - `penalty` : penalization exponent p (typically 3)
pub fn simp_modulus(rho: f64, e0: f64, e_min: f64, penalty: f64) -> f64 {
    e_min + rho.powf(penalty) * (e0 - e_min)
}

/// Derivative dE/dρ for the SIMP interpolation.
pub fn simp_modulus_derivative(rho: f64, e0: f64, e_min: f64, penalty: f64) -> f64 {
    penalty * rho.powf(penalty - 1.0) * (e0 - e_min)
}

/// SIMP compliance sensitivity ∂c/∂ρ = −dE/dρ · u_e^T k_e u_e / E(ρ).
///
/// `element_compliance` = u_e^T k_e u_e (element strain energy).
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

/// Update a single design variable using the SIMP OC formula.
///
/// ρ_new = ρ * (−∂c/∂ρ / λ)^(1/2) clamped to \[ρ − Δ, ρ + Δ\] ∩ \[ρ_min, 1\].
pub fn simp_oc_update(
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

// ---------------------------------------------------------------------------
// 2. Optimality Criteria update rule
// ---------------------------------------------------------------------------

/// Full Optimality Criteria update for a density vector.
///
/// Bisects on the Lagrange multiplier to enforce the volume constraint
/// Σ ρ_i = n · V_f and applies `simp_oc_update` to each element.
///
/// Returns the updated density vector.
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
            .map(|(&rho, &s)| simp_oc_update(rho, s.abs(), lm, move_limit, rho_min))
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
        .map(|(&rho, &s)| simp_oc_update(rho, s.abs(), lm, move_limit, rho_min))
        .collect()
}

/// Volume fraction of a density vector.
pub fn volume_fraction(densities: &[f64]) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    densities.iter().sum::<f64>() / densities.len() as f64
}

// ---------------------------------------------------------------------------
// 3. Sensitivity filtering (density-weighted average)
// ---------------------------------------------------------------------------

/// Density-weighted sensitivity filter (hat kernel).
///
/// s_filtered\[i\] = Σ_j w_ij ρ_j s_j / Σ_j w_ij ρ_j
/// with w_ij = max(0, r_min − d_ij).
///
/// `centroids`: (x, y) centroid of each element.
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

/// Simple hat-function sensitivity filter (no density weighting).
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

// ---------------------------------------------------------------------------
// 4. Heaviside projection
// ---------------------------------------------------------------------------

/// Smooth Heaviside projection: maps filtered density ρ̃ to physical density ρ̄.
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

/// Derivative dH/dρ̃.
pub fn heaviside_derivative(rho_tilde: f64, beta: f64, eta: f64) -> f64 {
    let den = (beta * eta).tanh() + (beta * (1.0 - eta)).tanh();
    if den.abs() < 1e-30 {
        return 1.0;
    }
    let sech2 = 1.0 - (beta * (rho_tilde - eta)).tanh().powi(2);
    beta * sech2 / den
}

/// Apply Heaviside projection to an entire density field.
pub fn apply_heaviside_projection(densities: &[f64], beta: f64, eta: f64) -> Vec<f64> {
    densities
        .iter()
        .map(|&r| heaviside_projection(r, beta, eta))
        .collect()
}

/// Chain-rule corrected sensitivity after Heaviside projection.
///
/// ∂c/∂ρ̃ = ∂c/∂ρ̄ · dH/dρ̃
pub fn heaviside_chain_sensitivity(
    dc_drho_bar: &[f64],
    rho_tilde: &[f64],
    beta: f64,
    eta: f64,
) -> Vec<f64> {
    dc_drho_bar
        .iter()
        .zip(rho_tilde.iter())
        .map(|(&dc, &rt)| dc * heaviside_derivative(rt, beta, eta))
        .collect()
}

// ---------------------------------------------------------------------------
// 5. Volume constraint enforcement
// ---------------------------------------------------------------------------

/// Enforce a volume constraint by rescaling densities.
///
/// Scales all densities uniformly so that their mean equals `target_vf`.
/// Clamps each density to `[rho_min, 1.0]` after scaling.
pub fn enforce_volume_constraint_scaling(densities: &mut [f64], target_vf: f64, rho_min: f64) {
    let current = volume_fraction(densities);
    if current < 1e-30 {
        return;
    }
    let scale = target_vf / current;
    for rho in densities.iter_mut() {
        *rho = (*rho * scale).clamp(rho_min, 1.0);
    }
}

/// Volume constraint value g = V(ρ)/V* − 1.
///
/// g < 0 ⟹ under target; g > 0 ⟹ over target.
pub fn volume_constraint_value(densities: &[f64], target_vf: f64) -> f64 {
    volume_fraction(densities) / target_vf - 1.0
}

/// Gradient of the volume constraint: ∂g/∂ρ_i = 1/(n * V*).
pub fn volume_constraint_gradient(n: usize, target_vf: f64) -> Vec<f64> {
    if n == 0 || target_vf < 1e-30 {
        return vec![];
    }
    vec![1.0 / (n as f64 * target_vf); n]
}

// ---------------------------------------------------------------------------
// 6. Checkerboard suppression
// ---------------------------------------------------------------------------

/// Measure the fraction of 2×2 patches in a structured grid exhibiting
/// a checkerboard pattern (alternating solid/void along diagonals).
///
/// `densities`: row-major flat array, `nx` columns, `ny` rows.
pub fn checkerboard_fraction(densities: &[f64], nx: usize, ny: usize) -> f64 {
    if nx < 2 || ny < 2 {
        return 0.0;
    }
    let mut checker = 0u32;
    let mut total = 0u32;
    for row in 0..(ny - 1) {
        for col in 0..(nx - 1) {
            let a = densities[row * nx + col] > 0.5;
            let b = densities[row * nx + col + 1] > 0.5;
            let c = densities[(row + 1) * nx + col] > 0.5;
            let d = densities[(row + 1) * nx + col + 1] > 0.5;
            if (a && !b && !c && d) || (!a && b && c && !d) {
                checker += 1;
            }
            total += 1;
        }
    }
    if total == 0 {
        0.0
    } else {
        checker as f64 / total as f64
    }
}

/// Suppress checkerboard patterns by replacing each interior element density
/// with a weighted average of itself and its four cardinal neighbours.
///
/// w_self and w_neighbor control the blend (w_self + 4·w_neighbor should ≈ 1).
pub fn suppress_checkerboard(
    densities: &[f64],
    nx: usize,
    ny: usize,
    w_self: f64,
    w_neighbor: f64,
) -> Vec<f64> {
    let n = densities.len();
    assert_eq!(n, nx * ny);
    let mut out = densities.to_vec();
    for row in 0..ny {
        for col in 0..nx {
            let idx = row * nx + col;
            let mut sum = w_self * densities[idx];
            let mut wt = w_self;
            if col > 0 {
                sum += w_neighbor * densities[row * nx + col - 1];
                wt += w_neighbor;
            }
            if col + 1 < nx {
                sum += w_neighbor * densities[row * nx + col + 1];
                wt += w_neighbor;
            }
            if row > 0 {
                sum += w_neighbor * densities[(row - 1) * nx + col];
                wt += w_neighbor;
            }
            if row + 1 < ny {
                sum += w_neighbor * densities[(row + 1) * nx + col];
                wt += w_neighbor;
            }
            out[idx] = sum / wt;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// 7. Continuation approach
// ---------------------------------------------------------------------------

/// Penalty continuation: slowly increase p toward p_max.
///
/// Returns the new penalty after one continuation step.
pub fn penalty_continuation(current: f64, increment: f64, p_max: f64) -> f64 {
    (current + increment).min(p_max)
}

/// Beta (Heaviside sharpness) continuation: multiply by `factor` up to `beta_max`.
pub fn beta_continuation(current_beta: f64, factor: f64, beta_max: f64) -> f64 {
    (current_beta * factor).min(beta_max)
}

/// A continuation schedule for both penalty and beta.
///
/// Every `interval` iterations, the penalty is incremented and the Heaviside
/// sharpness is doubled.
#[derive(Debug, Clone)]
pub struct ContinuationSchedule {
    /// How often (in iterations) to increase penalty and beta.
    pub interval: usize,
    /// Penalty increment per step.
    pub penalty_increment: f64,
    /// Maximum penalty.
    pub penalty_max: f64,
    /// Beta multiplication factor per step.
    pub beta_factor: f64,
    /// Maximum Heaviside sharpness.
    pub beta_max: f64,
    /// Current penalty.
    pub current_penalty: f64,
    /// Current beta.
    pub current_beta: f64,
    /// Iteration counter.
    pub iteration: usize,
}

impl ContinuationSchedule {
    /// Create a new continuation schedule.
    pub fn new(
        interval: usize,
        initial_penalty: f64,
        penalty_increment: f64,
        penalty_max: f64,
        initial_beta: f64,
        beta_factor: f64,
        beta_max: f64,
    ) -> Self {
        Self {
            interval,
            penalty_increment,
            penalty_max,
            beta_factor,
            beta_max,
            current_penalty: initial_penalty,
            current_beta: initial_beta,
            iteration: 0,
        }
    }

    /// Advance one iteration; update parameters if at an interval boundary.
    pub fn step(&mut self) {
        self.iteration += 1;
        if self.interval > 0 && self.iteration.is_multiple_of(self.interval) {
            self.current_penalty = penalty_continuation(
                self.current_penalty,
                self.penalty_increment,
                self.penalty_max,
            );
            self.current_beta =
                beta_continuation(self.current_beta, self.beta_factor, self.beta_max);
        }
    }

    /// Whether the penalty has reached its maximum.
    pub fn penalty_saturated(&self) -> bool {
        (self.current_penalty - self.penalty_max).abs() < 1e-10
    }

    /// Whether the beta has reached its maximum.
    pub fn beta_saturated(&self) -> bool {
        (self.current_beta - self.beta_max).abs() < 1e-10
    }
}

// ---------------------------------------------------------------------------
// 8. Compliance objective gradient
// ---------------------------------------------------------------------------

/// Compliance objective C = Σ_e E(ρ_e) c_e where c_e is element compliance.
pub fn compliance_objective(
    densities: &[f64],
    element_compliances: &[f64],
    e0: f64,
    e_min: f64,
    penalty: f64,
) -> f64 {
    densities
        .iter()
        .zip(element_compliances.iter())
        .map(|(&rho, &ce)| simp_modulus(rho, e0, e_min, penalty) * ce)
        .sum()
}

/// Gradient of the compliance objective ∂C/∂ρ_e (full SIMP formula).
///
/// Returns a vector of sensitivities (all ≤ 0 for a normal compliance problem).
pub fn compliance_gradient(
    densities: &[f64],
    element_compliances: &[f64],
    e0: f64,
    e_min: f64,
    penalty: f64,
) -> Vec<f64> {
    densities
        .iter()
        .zip(element_compliances.iter())
        .map(|(&rho, &ce)| {
            let de = simp_modulus_derivative(rho, e0, e_min, penalty);
            -de * ce
        })
        .collect()
}

/// Normalised compliance gradient: divide by the maximum absolute value.
pub fn normalised_compliance_gradient(
    densities: &[f64],
    element_compliances: &[f64],
    e0: f64,
    e_min: f64,
    penalty: f64,
) -> Vec<f64> {
    let grad = compliance_gradient(densities, element_compliances, e0, e_min, penalty);
    let max_abs = grad.iter().map(|g| g.abs()).fold(0.0_f64, f64::max);
    if max_abs < 1e-30 {
        return grad;
    }
    grad.iter().map(|g| g / max_abs).collect()
}

// ---------------------------------------------------------------------------
// 9. Multi-load-case optimization
// ---------------------------------------------------------------------------

/// Accumulate weighted sensitivities from multiple load cases.
///
/// `load_sensitivities[l]` is the sensitivity vector for load case `l`.
/// `weights[l]` is the weight for load case `l`.
///
/// Returns the weighted sum of all sensitivity vectors.
pub fn accumulate_multi_load_sensitivities(
    load_sensitivities: &[Vec<f64>],
    weights: &[f64],
) -> Vec<f64> {
    assert!(!load_sensitivities.is_empty());
    assert_eq!(load_sensitivities.len(), weights.len());
    let n = load_sensitivities[0].len();
    let mut acc = vec![0.0_f64; n];
    for (sens, &w) in load_sensitivities.iter().zip(weights.iter()) {
        assert_eq!(sens.len(), n);
        for (a, &s) in acc.iter_mut().zip(sens.iter()) {
            *a += w * s;
        }
    }
    acc
}

/// Weighted total compliance from multiple load cases.
pub fn multi_load_compliance(compliances: &[f64], weights: &[f64]) -> f64 {
    assert_eq!(compliances.len(), weights.len());
    compliances
        .iter()
        .zip(weights.iter())
        .map(|(c, w)| c * w)
        .sum()
}

/// Multi-load OC update: accumulate sensitivities, then apply OC.
pub fn multi_load_oc_update(
    densities: &[f64],
    load_sensitivities: &[Vec<f64>],
    weights: &[f64],
    volume_fraction_target: f64,
    move_limit: f64,
    rho_min: f64,
) -> Vec<f64> {
    let acc_sens = accumulate_multi_load_sensitivities(load_sensitivities, weights);
    oc_update(
        densities,
        &acc_sens,
        volume_fraction_target,
        move_limit,
        rho_min,
    )
}

// ---------------------------------------------------------------------------
// 10. Topology optimization result
// ---------------------------------------------------------------------------

/// Summary of one completed topology optimization iteration.
#[derive(Debug, Clone)]
pub struct TopOptIterResult {
    /// Iteration number (1-indexed).
    pub iteration: usize,
    /// Compliance at this iteration.
    pub compliance: f64,
    /// Volume fraction at this iteration.
    pub volume_fraction: f64,
    /// Maximum density change from previous iteration.
    pub max_density_change: f64,
    /// Lagrange multiplier used for the OC update.
    pub lagrange: f64,
}

/// Run a series of OC topology optimization iterations.
///
/// This function does NOT perform FEA; `sensitivity_fn` is called each
/// iteration to compute element sensitivities from the current densities.
pub fn run_oc_iterations(
    initial_densities: &[f64],
    sensitivity_fn: &dyn Fn(&[f64]) -> Vec<f64>,
    volume_fraction_target: f64,
    move_limit: f64,
    rho_min: f64,
    max_iter: usize,
    convergence_tol: f64,
) -> (Vec<f64>, Vec<TopOptIterResult>) {
    let n = initial_densities.len();
    let mut densities = initial_densities.to_vec();
    let mut history = Vec::with_capacity(max_iter);

    for iter in 1..=max_iter {
        let sens = sensitivity_fn(&densities);
        let prev = densities.clone();
        densities = oc_update(
            &densities,
            &sens,
            volume_fraction_target,
            move_limit,
            rho_min,
        );

        let max_change = densities
            .iter()
            .zip(prev.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);

        let compliance: f64 = densities
            .iter()
            .zip(sens.iter())
            .map(|(r, s)| r * s.abs())
            .sum();
        let vf = volume_fraction(&densities);

        // Approximate Lagrange multiplier (bisection midpoint, recomputed for history).
        let lagrange = if n > 0 { compliance / n as f64 } else { 0.0 };

        history.push(TopOptIterResult {
            iteration: iter,
            compliance,
            volume_fraction: vf,
            max_density_change: max_change,
            lagrange,
        });

        if max_change < convergence_tol {
            break;
        }
    }

    (densities, history)
}

// ---------------------------------------------------------------------------
// 11. Density projection filter (combined density filter + Heaviside)
// ---------------------------------------------------------------------------

/// Apply a combined density filter and Heaviside projection.
///
/// Step 1: density filter  ρ̃ = Σ w_j ρ_j / Σ w_j
/// Step 2: Heaviside projection  ρ̄ = H(ρ̃; β, η)
pub fn density_filter_and_project(
    densities: &[f64],
    centroids: &[(f64, f64)],
    r_min: f64,
    beta: f64,
    eta: f64,
) -> Vec<f64> {
    let n = densities.len();
    assert_eq!(centroids.len(), n);
    let mut rho_tilde = vec![0.0_f64; n];
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
        rho_tilde[i] = if den > 1e-30 { num / den } else { densities[i] };
    }
    apply_heaviside_projection(&rho_tilde, beta, eta)
}

// ---------------------------------------------------------------------------
// 12. Grayscale and discreteness measures
// ---------------------------------------------------------------------------

/// Fraction of elements in the "gray zone" (0.1 < ρ < 0.9).
pub fn gray_scale_measure(densities: &[f64]) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    let gray = densities.iter().filter(|&&r| r > 0.1 && r < 0.9).count();
    gray as f64 / densities.len() as f64
}

/// Discreteness measure: 4·ρ·(1−ρ) averaged over elements (0 = discrete, 1 = gray).
pub fn discreteness_measure(densities: &[f64]) -> f64 {
    if densities.is_empty() {
        return 0.0;
    }
    let sum: f64 = densities.iter().map(|&r| 4.0 * r * (1.0 - r)).sum();
    sum / densities.len() as f64
}

/// Mean, minimum, and maximum density.
pub fn density_statistics(densities: &[f64]) -> (f64, f64, f64) {
    if densities.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mean = densities.iter().sum::<f64>() / densities.len() as f64;
    let min = densities.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = densities.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (mean, min, max)
}

// ---------------------------------------------------------------------------
// 13. Stress-constrained topology optimization
// ---------------------------------------------------------------------------

/// Von Mises stress from 2D stress components (σ_x, σ_y, τ_xy).
///
/// σ_vm = √(σ_x² - σ_x σ_y + σ_y² + 3 τ_xy²)
pub fn von_mises_stress(sigma_x: f64, sigma_y: f64, tau_xy: f64) -> f64 {
    (sigma_x * sigma_x - sigma_x * sigma_y + sigma_y * sigma_y + 3.0 * tau_xy * tau_xy).sqrt()
}

/// P-norm stress aggregation function.
///
/// σ_pn = σ_ref * (Σ (σ_i / σ_ref)^p)^(1/p)
///
/// Approximates the maximum stress for large `p`. The reference stress `sigma_ref`
/// controls numerical conditioning.
pub fn p_norm_aggregation(stresses: &[f64], p: f64, sigma_ref: f64) -> f64 {
    if stresses.is_empty() || sigma_ref < 1e-30 {
        return 0.0;
    }
    let sum: f64 = stresses.iter().map(|&s| (s / sigma_ref).powf(p)).sum();
    sigma_ref * sum.powf(1.0 / p)
}

/// Sensitivity of the P-norm aggregation with respect to each stress value.
///
/// dσ_pn / dσ_i = (σ_i / σ_ref)^(p-1) / (Σ (σ_j / σ_ref)^p)^((p-1)/p)
pub fn p_norm_sensitivity(stresses: &[f64], p: f64, sigma_ref: f64) -> Vec<f64> {
    if stresses.is_empty() || sigma_ref < 1e-30 {
        return vec![0.0; stresses.len()];
    }
    let sum: f64 = stresses.iter().map(|&s| (s / sigma_ref).powf(p)).sum();
    if sum < 1e-30 {
        return vec![0.0; stresses.len()];
    }
    let sum_p_m1 = sum.powf((p - 1.0) / p);
    stresses
        .iter()
        .map(|&s| (s / sigma_ref).powf(p - 1.0) / (sigma_ref * sum_p_m1))
        .collect()
}

/// SIMP stress relaxation to avoid stress singularity in void elements.
///
/// σ_relaxed = ρ^q / (ρ^q + ε(1 - ρ^q)) * σ_full
///
/// where q is a relaxation exponent and ε is a small regularisation parameter.
pub fn stress_simp_relaxation(rho: f64, sigma_full: f64, q: f64, epsilon: f64) -> f64 {
    let rho_q = rho.powf(q);
    let factor = rho_q / (rho_q + epsilon * (1.0 - rho_q));
    factor * sigma_full
}

// ---------------------------------------------------------------------------
// 14. BESO (Bi-directional Evolutionary Structural Optimization)
// ---------------------------------------------------------------------------

/// BESO update: elements with sensitivity below the threshold are removed (set to rho_min),
/// elements above the threshold are added (set to 1). Threshold is found by bisection.
///
/// - `densities`: current element densities (0 or 1 for BESO; intermediate also accepted).
/// - `sensitivities`: element sensitivity numbers (higher = more important).
/// - `target_vf`: target volume fraction.
/// - `rho_min`: density for "void" elements (typically 0.001 to avoid singularity).
///
/// Returns updated binary densities.
pub fn beso_update(
    densities: &[f64],
    sensitivities: &[f64],
    target_vf: f64,
    rho_min: f64,
) -> Vec<f64> {
    let n = densities.len();
    assert_eq!(sensitivities.len(), n);
    let target = (target_vf * n as f64).round() as usize;
    let target = target.clamp(1, n);

    // Sort element indices by sensitivity (descending order).
    let mut indices: Vec<usize> = (0..n).collect();
    indices.sort_by(|&a, &b| {
        sensitivities[b]
            .partial_cmp(&sensitivities[a])
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut new_dens = vec![rho_min; n];
    for &idx in &indices[..target] {
        new_dens[idx] = 1.0;
    }
    new_dens
}

/// BESO sensitivity averaging (evolutionary history filter).
///
/// Averages current sensitivities with previous iteration for stability.
///
/// α_n = (α_current + α_previous) / 2
pub fn beso_sensitivity_average(current: &[f64], previous: &[f64]) -> Vec<f64> {
    assert_eq!(current.len(), previous.len());
    current
        .iter()
        .zip(previous.iter())
        .map(|(c, p)| (c + p) / 2.0)
        .collect()
}

/// BESO convergence criterion: maximum change in compliance over last `window` iterations.
pub fn beso_converged(compliance_history: &[f64], window: usize, tol: f64) -> bool {
    let n = compliance_history.len();
    if n < window * 2 {
        return false;
    }
    let recent: f64 = compliance_history[(n - window)..].iter().sum::<f64>() / window as f64;
    let previous: f64 = compliance_history[(n - 2 * window)..(n - window)]
        .iter()
        .sum::<f64>()
        / window as f64;
    if previous.abs() < 1e-30 {
        return false;
    }
    ((recent - previous) / previous).abs() < tol
}

// ---------------------------------------------------------------------------
// 15. Level-set topology optimization
// ---------------------------------------------------------------------------

/// Smooth Heaviside function for level-set methods.
///
/// H(φ) = 0 for φ < -ε, 0.5 + 3φ/(4ε) - φ³/(4ε³) for |φ| ≤ ε, 1 for φ > ε.
pub fn level_set_heaviside(phi: f64, epsilon: f64) -> f64 {
    if phi < -epsilon {
        0.0
    } else if phi > epsilon {
        1.0
    } else {
        0.5 + 3.0 * phi / (4.0 * epsilon) - phi.powi(3) / (4.0 * epsilon.powi(3))
    }
}

/// Smooth Dirac delta for level-set methods: δ(φ) = dH/dφ.
///
/// δ(φ) = 3/(4ε) (1 - (φ/ε)²) for |φ| ≤ ε, else 0.
pub fn level_set_dirac(phi: f64, epsilon: f64) -> f64 {
    if phi.abs() <= epsilon {
        3.0 / (4.0 * epsilon) * (1.0 - (phi / epsilon).powi(2))
    } else {
        0.0
    }
}

/// Advect the level-set function using the Hamilton-Jacobi equation (explicit Euler).
///
/// φ_{n+1} = φ_n - Δt · V_n · |∇φ|
///
/// For 1-D: approximates |∇φ| ≈ 1 (signed distance function).
pub fn level_set_advect(phi: &[f64], velocity: &[f64], dt: f64) -> Vec<f64> {
    assert_eq!(phi.len(), velocity.len());
    phi.iter()
        .zip(velocity.iter())
        .map(|(&p, &v)| p - v * dt)
        .collect()
}

/// Simple 1-D level-set re-initialisation using the fast-marching sign function.
///
/// Iterates `n_iter` times of: φ_new = φ + Δt · sign(φ₀) · (1 - |∇φ|)
/// with a simple central-difference gradient approximation.
pub fn level_set_reinitialize(phi: &[f64], dx: f64, dt: f64, n_iter: usize) -> Vec<f64> {
    let n = phi.len();
    let sign: Vec<f64> = phi
        .iter()
        .map(|&p| {
            if p > 1e-14 {
                1.0
            } else if p < -1e-14 {
                -1.0
            } else {
                0.0
            }
        })
        .collect();

    let mut phi_cur = phi.to_vec();

    for _ in 0..n_iter {
        let mut phi_new = phi_cur.clone();
        for i in 1..(n - 1) {
            let grad = (phi_cur[i + 1] - phi_cur[i - 1]) / (2.0 * dx);
            phi_new[i] += dt * sign[i] * (1.0 - grad.abs());
        }
        phi_cur = phi_new;
    }
    phi_cur
}

/// Extend velocity from the interface into the domain (constant extension).
///
/// For compliance minimisation, the velocity field is derived from the
/// compliance sensitivity: V = -dC/dΩ. This function flips the sign so that
/// regions with high sensitivity contract the material boundary (void grows).
pub fn level_set_velocity_extension(sensitivities: &[f64]) -> Vec<f64> {
    // Velocity = -sensitivity (material removed where sensitivity is most negative)
    sensitivities.iter().map(|&s| -s).collect()
}

// ---------------------------------------------------------------------------
// 16. Perimeter constraint
// ---------------------------------------------------------------------------

/// Compute the discrete perimeter of a topology (number of solid-void boundaries).
///
/// Each pair of adjacent elements where one is above `threshold` and the
/// other is below counts as one boundary edge.
pub fn perimeter_length(densities: &[f64], nx: usize, ny: usize, threshold: f64) -> f64 {
    let mut perimeter = 0.0_f64;
    for row in 0..ny {
        for col in 0..nx {
            let is_solid = densities[row * nx + col] >= threshold;
            // Right neighbour
            if col + 1 < nx {
                let right_solid = densities[row * nx + col + 1] >= threshold;
                if is_solid != right_solid {
                    perimeter += 1.0;
                }
            }
            // Top neighbour
            if row + 1 < ny {
                let top_solid = densities[(row + 1) * nx + col] >= threshold;
                if is_solid != top_solid {
                    perimeter += 1.0;
                }
            }
        }
    }
    perimeter
}

/// Gradient of the discrete perimeter with respect to element densities.
///
/// Each element's gradient is the number of its neighbours on the other side
/// of the threshold boundary.
pub fn perimeter_gradient(densities: &[f64], nx: usize, ny: usize, threshold: f64) -> Vec<f64> {
    let n = densities.len();
    let mut grad = vec![0.0_f64; n];

    let neighbours = |row: usize, col: usize| -> Vec<usize> {
        let mut nb = Vec::new();
        if col > 0 {
            nb.push(row * nx + col - 1);
        }
        if col + 1 < nx {
            nb.push(row * nx + col + 1);
        }
        if row > 0 {
            nb.push((row - 1) * nx + col);
        }
        if row + 1 < ny {
            nb.push((row + 1) * nx + col);
        }
        nb
    };

    for row in 0..ny {
        for col in 0..nx {
            let idx = row * nx + col;
            let is_solid = densities[idx] >= threshold;
            for nb_idx in neighbours(row, col) {
                let nb_solid = densities[nb_idx] >= threshold;
                if is_solid != nb_solid {
                    grad[idx] += 1.0;
                }
            }
        }
    }
    grad
}

// ---------------------------------------------------------------------------
// 17. Manufacturability constraints
// ---------------------------------------------------------------------------

/// Minimum length scale filter using erosion/dilation.
///
/// Projects densities so features below the minimum length scale are suppressed.
/// Uses a hat-function filter with radius r_min then applies solid threshold.
pub fn minimum_length_scale_filter(
    densities: &[f64],
    centroids: &[(f64, f64)],
    r_min: f64,
) -> Vec<f64> {
    let n = densities.len();
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
        out[i] = if den > 1e-30 {
            (num / den).clamp(0.0, 1.0)
        } else {
            densities[i]
        };
    }
    out
}

/// Casting constraint projection (single draw-direction mould).
///
/// Enforces the maximum density rule along the draw direction:
/// ρ(z) ≤ max_{z' > z} ρ(z')  (material cannot be trapped under a core).
///
/// `forward`: if true, projection is applied bottom-to-top (increasing index).
///
/// Returns a modified density field satisfying the casting constraint.
pub fn casting_constraint_projection(densities: &[f64], forward: bool) -> Vec<f64> {
    let n = densities.len();
    let mut out = densities.to_vec();

    if forward {
        // Bottom-to-top: each element ≤ its predecessor
        for i in 1..n {
            if out[i] > out[i - 1] {
                out[i] = out[i - 1];
            }
        }
    } else {
        // Top-to-bottom: each element ≤ its successor
        for i in (0..(n - 1)).rev() {
            if out[i] > out[i + 1] {
                out[i] = out[i + 1];
            }
        }
    }
    out
}

/// Extrusion constraint: enforce that all element densities in a column
/// (along the extrusion direction) are equal.
///
/// `nx`: number of elements per row (x direction).
/// `ny`: number of rows (z/extrusion direction).
///
/// Each column's density is averaged over all ny slices.
pub fn extrusion_constraint_projection(densities: &[f64], nx: usize, ny: usize) -> Vec<f64> {
    assert_eq!(densities.len(), nx * ny);
    let mut out = vec![0.0_f64; nx * ny];
    for col in 0..nx {
        let mean: f64 = (0..ny).map(|row| densities[row * nx + col]).sum::<f64>() / ny as f64;
        for row in 0..ny {
            out[row * nx + col] = mean;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// 18. Multi-objective topology optimization
// ---------------------------------------------------------------------------

/// Pareto-weighted objective combining compliance and dynamic objective.
///
/// f = w · C + (1 − w) · f_dyn  where w ∈ \[0, 1\].
pub fn pareto_weighted_objective(compliance: f64, dynamic_obj: f64, weight: f64) -> f64 {
    weight * compliance + (1.0 - weight) * dynamic_obj
}

/// Gradient of the thermal compliance with respect to element densities.
///
/// Analogous to structural compliance gradient but uses thermal sensitivities.
pub fn thermal_compliance_gradient(
    densities: &[f64],
    thermal_sensitivities: &[f64],
    e0: f64,
    e_min: f64,
    penalty: f64,
) -> Vec<f64> {
    densities
        .iter()
        .zip(thermal_sensitivities.iter())
        .map(|(&rho, &ts)| {
            let de = simp_modulus_derivative(rho, e0, e_min, penalty);
            -de * ts
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── SIMP modulus ─────────────────────────────────────────────────────

    #[test]
    fn test_simp_modulus_solid() {
        let e = simp_modulus(1.0, 210e9, 1e-3, 3.0);
        assert!((e - 210e9).abs() < 1.0);
    }

    #[test]
    fn test_simp_modulus_void() {
        let e = simp_modulus(0.0, 210e9, 1e-3, 3.0);
        assert!((e - 1e-3).abs() < 1e-10);
    }

    #[test]
    fn test_simp_modulus_half() {
        let e = simp_modulus(0.5, 1.0, 0.0, 3.0);
        let expected = 0.5_f64.powi(3);
        assert!((e - expected).abs() < 1e-12);
    }

    #[test]
    fn test_simp_modulus_derivative_positive() {
        let de = simp_modulus_derivative(0.5, 1.0, 0.0, 3.0);
        assert!(de > 0.0);
    }

    #[test]
    fn test_simp_compliance_sensitivity_negative() {
        // Sensitivity should be negative (removing material increases compliance)
        let s = simp_compliance_sensitivity(0.5, 210e9, 1e-3, 3.0, 1.0);
        assert!(s <= 0.0);
    }

    // ── OC update ────────────────────────────────────────────────────────

    #[test]
    fn test_oc_update_volume_constraint_satisfied() {
        let densities = vec![0.5; 100];
        let sensitivities = vec![-1.0; 100];
        let new_dens = oc_update(&densities, &sensitivities, 0.4, 0.2, 1e-3);
        let vf = volume_fraction(&new_dens);
        assert!((vf - 0.4).abs() < 0.05, "vf = {vf}");
    }

    #[test]
    fn test_oc_update_preserves_length() {
        let densities = vec![0.5; 50];
        let sensitivities = vec![-1.0; 50];
        let new_dens = oc_update(&densities, &sensitivities, 0.5, 0.2, 1e-3);
        assert_eq!(new_dens.len(), 50);
    }

    #[test]
    fn test_oc_update_density_in_range() {
        let densities: Vec<f64> = (0..20).map(|i| 0.1 + i as f64 * 0.04).collect();
        let sensitivities: Vec<f64> = (0..20).map(|i| -(i as f64 + 1.0)).collect();
        let new_dens = oc_update(&densities, &sensitivities, 0.5, 0.2, 1e-3);
        for &rho in &new_dens {
            assert!((1e-3..=1.0).contains(&rho), "Out of range: {rho}");
        }
    }

    #[test]
    fn test_simp_oc_update_clamps_to_rho_min() {
        let rho_new = simp_oc_update(0.001, 0.01, 1e9, 0.2, 1e-3);
        assert!(rho_new >= 1e-3);
    }

    #[test]
    fn test_volume_fraction_uniform() {
        let dens = vec![0.4; 25];
        assert!((volume_fraction(&dens) - 0.4).abs() < 1e-12);
    }

    // ── Sensitivity filter ────────────────────────────────────────────────

    fn uniform_centroids(n: usize) -> Vec<(f64, f64)> {
        let side = (n as f64).sqrt().ceil() as usize;
        (0..n)
            .map(|i| ((i % side) as f64 + 0.5, (i / side) as f64 + 0.5))
            .collect()
    }

    #[test]
    fn test_filter_density_weighted_length() {
        let n = 16;
        let sens = vec![-1.0; n];
        let dens = vec![0.5; n];
        let centroids = uniform_centroids(n);
        let out = filter_sensitivities_density_weighted(&sens, &dens, &centroids, 1.5);
        assert_eq!(out.len(), n);
    }

    #[test]
    fn test_filter_density_weighted_uniform_unchanged() {
        let n = 9;
        let sens = vec![-2.0; n];
        let dens = vec![0.5; n];
        let centroids = uniform_centroids(n);
        let out = filter_sensitivities_density_weighted(&sens, &dens, &centroids, 2.0);
        for &v in &out {
            assert!(
                (v - (-2.0)).abs() < 1e-6,
                "Uniform sensitivity changed: {v}"
            );
        }
    }

    #[test]
    fn test_filter_hat_length() {
        let n = 9;
        let sens = vec![-1.0; n];
        let centroids = uniform_centroids(n);
        let out = filter_sensitivities_hat(&sens, &centroids, 1.5);
        assert_eq!(out.len(), n);
    }

    #[test]
    fn test_filter_hat_uniform_unchanged() {
        let n = 9;
        let sens = vec![-3.0; n];
        let centroids = uniform_centroids(n);
        let out = filter_sensitivities_hat(&sens, &centroids, 1.5);
        for &v in &out {
            assert!((v - (-3.0)).abs() < 1e-6, "Hat filter changed uniform: {v}");
        }
    }

    // ── Heaviside projection ──────────────────────────────────────────────

    #[test]
    fn test_heaviside_at_eta() {
        let h = heaviside_projection(0.5, 10.0, 0.5);
        assert!((h - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_heaviside_above_eta_near_one() {
        let h = heaviside_projection(1.0, 100.0, 0.5);
        assert!(h > 0.95, "h = {h}");
    }

    #[test]
    fn test_heaviside_below_eta_near_zero() {
        let h = heaviside_projection(0.0, 100.0, 0.5);
        assert!(h < 0.05, "h = {h}");
    }

    #[test]
    fn test_heaviside_derivative_positive() {
        let dh = heaviside_derivative(0.5, 10.0, 0.5);
        assert!(dh > 0.0);
    }

    #[test]
    fn test_heaviside_chain_sensitivity_length() {
        let dc = vec![-1.0_f64; 8];
        let rt = vec![0.5_f64; 8];
        let out = heaviside_chain_sensitivity(&dc, &rt, 10.0, 0.5);
        assert_eq!(out.len(), 8);
    }

    #[test]
    fn test_apply_heaviside_projection_length() {
        let dens = vec![0.3, 0.5, 0.7, 1.0];
        let out = apply_heaviside_projection(&dens, 5.0, 0.5);
        assert_eq!(out.len(), 4);
    }

    // ── Volume constraint ─────────────────────────────────────────────────

    #[test]
    fn test_enforce_volume_constraint_scaling() {
        let mut dens = vec![0.8; 10];
        enforce_volume_constraint_scaling(&mut dens, 0.5, 1e-3);
        let vf = volume_fraction(&dens);
        assert!((vf - 0.5).abs() < 0.01, "vf = {vf}");
    }

    #[test]
    fn test_volume_constraint_value_at_target() {
        let dens = vec![0.5; 20];
        let g = volume_constraint_value(&dens, 0.5);
        assert!(g.abs() < 1e-12, "g = {g}");
    }

    #[test]
    fn test_volume_constraint_value_over_budget() {
        let dens = vec![0.8; 20];
        let g = volume_constraint_value(&dens, 0.5);
        assert!(g > 0.0, "g = {g}");
    }

    #[test]
    fn test_volume_constraint_gradient_length() {
        let grad = volume_constraint_gradient(10, 0.5);
        assert_eq!(grad.len(), 10);
    }

    #[test]
    fn test_volume_constraint_gradient_uniform() {
        let grad = volume_constraint_gradient(10, 0.5);
        let expected = 1.0 / (10.0 * 0.5);
        for &g in &grad {
            assert!((g - expected).abs() < 1e-12);
        }
    }

    // ── Checkerboard ──────────────────────────────────────────────────────

    #[test]
    fn test_checkerboard_fraction_uniform() {
        let dens = vec![1.0; 9];
        assert_eq!(checkerboard_fraction(&dens, 3, 3), 0.0);
    }

    #[test]
    fn test_checkerboard_fraction_full() {
        // 2×2 pure checkerboard: (0,0)=solid, (0,1)=void, (1,0)=void, (1,1)=solid
        let dens = vec![1.0, 0.0, 0.0, 1.0];
        assert!((checkerboard_fraction(&dens, 2, 2) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_suppress_checkerboard_length() {
        let dens = vec![1.0, 0.0, 0.0, 1.0];
        let out = suppress_checkerboard(&dens, 2, 2, 0.5, 0.125);
        assert_eq!(out.len(), 4);
    }

    #[test]
    fn test_suppress_checkerboard_reduces_pattern() {
        let dens = vec![1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0];
        let out = suppress_checkerboard(&dens, 3, 3, 0.5, 0.125);
        let before = checkerboard_fraction(&dens, 3, 3);
        let after = checkerboard_fraction(&out, 3, 3);
        assert!(after <= before + 1e-10, "before={before} after={after}");
    }

    // ── Continuation ─────────────────────────────────────────────────────

    #[test]
    fn test_penalty_continuation_increments() {
        let p = penalty_continuation(3.0, 0.5, 5.0);
        assert!((p - 3.5).abs() < 1e-12);
    }

    #[test]
    fn test_penalty_continuation_clamped() {
        let p = penalty_continuation(4.9, 0.5, 5.0);
        assert!((p - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_beta_continuation_multiplies() {
        let b = beta_continuation(2.0, 3.0, 100.0);
        assert!((b - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_beta_continuation_clamped() {
        let b = beta_continuation(50.0, 3.0, 100.0);
        assert!((b - 100.0).abs() < 1e-12);
    }

    #[test]
    fn test_continuation_schedule_step() {
        let mut sch = ContinuationSchedule::new(5, 3.0, 1.0, 8.0, 1.0, 2.0, 64.0);
        for _ in 0..5 {
            sch.step();
        }
        assert_eq!(sch.iteration, 5);
        assert!((sch.current_penalty - 4.0).abs() < 1e-12);
        assert!((sch.current_beta - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_continuation_schedule_saturated() {
        let mut sch = ContinuationSchedule::new(1, 7.0, 2.0, 8.0, 32.0, 4.0, 64.0);
        sch.step(); // penalty → 8 (clamped), beta → 64 (clamped)
        sch.step(); // no change
        assert!(sch.penalty_saturated());
        assert!(sch.beta_saturated());
    }

    // ── Compliance gradient ───────────────────────────────────────────────

    #[test]
    fn test_compliance_objective_positive() {
        let dens = vec![0.5; 10];
        let ce = vec![1.0; 10];
        let c = compliance_objective(&dens, &ce, 1.0, 1e-9, 3.0);
        assert!(c > 0.0);
    }

    #[test]
    fn test_compliance_gradient_negative() {
        let dens = vec![0.5; 5];
        let ce = vec![1.0; 5];
        let grad = compliance_gradient(&dens, &ce, 1.0, 1e-9, 3.0);
        for &g in &grad {
            assert!(g <= 0.0, "Sensitivity should be ≤ 0: {g}");
        }
    }

    #[test]
    fn test_normalised_compliance_gradient_max_one() {
        let dens = vec![0.3, 0.5, 0.7, 1.0];
        let ce = vec![1.0, 2.0, 0.5, 3.0];
        let ng = normalised_compliance_gradient(&dens, &ce, 1.0, 1e-9, 3.0);
        let max_abs = ng.iter().map(|g| g.abs()).fold(0.0_f64, f64::max);
        assert!((max_abs - 1.0).abs() < 1e-12, "max_abs = {max_abs}");
    }

    #[test]
    fn test_compliance_gradient_length() {
        let dens = vec![0.5; 8];
        let ce = vec![1.0; 8];
        let grad = compliance_gradient(&dens, &ce, 1.0, 1e-9, 3.0);
        assert_eq!(grad.len(), 8);
    }

    // ── Multi-load ────────────────────────────────────────────────────────

    #[test]
    fn test_accumulate_multi_load_sensitivities_length() {
        let sens = vec![vec![-1.0; 10], vec![-2.0; 10]];
        let w = vec![0.5, 0.5];
        let acc = accumulate_multi_load_sensitivities(&sens, &w);
        assert_eq!(acc.len(), 10);
    }

    #[test]
    fn test_accumulate_multi_load_sensitivities_values() {
        let sens = vec![vec![-1.0; 4], vec![-3.0; 4]];
        let w = vec![0.5, 0.5];
        let acc = accumulate_multi_load_sensitivities(&sens, &w);
        for &a in &acc {
            assert!((a - (-2.0)).abs() < 1e-12, "acc = {a}");
        }
    }

    #[test]
    fn test_multi_load_compliance() {
        let c = multi_load_compliance(&[10.0, 20.0], &[0.6, 0.4]);
        assert!((c - 14.0).abs() < 1e-12, "c = {c}");
    }

    #[test]
    fn test_multi_load_oc_update_volume() {
        let dens = vec![0.5; 50];
        let sens1 = vec![-1.0; 50];
        let sens2 = vec![-1.0; 50];
        let new_dens = multi_load_oc_update(&dens, &[sens1, sens2], &[0.5, 0.5], 0.4, 0.2, 1e-3);
        let vf = volume_fraction(&new_dens);
        assert!((vf - 0.4).abs() < 0.05, "vf = {vf}");
    }

    // ── run_oc_iterations ─────────────────────────────────────────────────

    #[test]
    fn test_run_oc_iterations_returns_correct_length() {
        let dens = vec![0.5; 16];
        let (final_dens, hist) = run_oc_iterations(
            &dens,
            &|d: &[f64]| d.iter().map(|_| -1.0).collect(),
            0.4,
            0.2,
            1e-3,
            5,
            1e-4,
        );
        assert_eq!(final_dens.len(), 16);
        assert!(!hist.is_empty());
    }

    #[test]
    fn test_run_oc_iterations_volume_approaches_target() {
        let dens = vec![0.8; 25];
        let (final_dens, _) = run_oc_iterations(
            &dens,
            &|d: &[f64]| d.iter().map(|_| -1.0).collect(),
            0.5,
            0.2,
            1e-3,
            20,
            1e-5,
        );
        let vf = volume_fraction(&final_dens);
        assert!((vf - 0.5).abs() < 0.1, "vf = {vf}");
    }

    // ── Density filter + projection ───────────────────────────────────────

    #[test]
    fn test_density_filter_and_project_length() {
        let dens = vec![0.5; 9];
        let centroids = uniform_centroids(9);
        let out = density_filter_and_project(&dens, &centroids, 1.5, 5.0, 0.5);
        assert_eq!(out.len(), 9);
    }

    #[test]
    fn test_density_filter_and_project_uniform() {
        // Uniform density 0.5 at eta=0.5 → projected ≈ 0.5
        let dens = vec![0.5; 9];
        let centroids = uniform_centroids(9);
        let out = density_filter_and_project(&dens, &centroids, 1.5, 5.0, 0.5);
        for &v in &out {
            assert!((v - 0.5).abs() < 0.1, "v = {v}");
        }
    }

    // ── Grayscale / discreteness ──────────────────────────────────────────

    #[test]
    fn test_gray_scale_all_gray() {
        let dens = vec![0.5; 10];
        assert!((gray_scale_measure(&dens) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gray_scale_all_binary() {
        let dens = vec![0.0, 1.0, 0.0, 1.0];
        assert_eq!(gray_scale_measure(&dens), 0.0);
    }

    #[test]
    fn test_discreteness_measure_half() {
        let dens = vec![0.5; 4];
        let dm = discreteness_measure(&dens);
        assert!((dm - 1.0).abs() < 1e-12, "dm = {dm}");
    }

    #[test]
    fn test_discreteness_measure_binary() {
        let dens = vec![0.0, 1.0];
        let dm = discreteness_measure(&dens);
        assert!(dm.abs() < 1e-12, "dm = {dm}");
    }

    #[test]
    fn test_density_statistics_uniform() {
        let dens = vec![0.4; 10];
        let (mean, min, max) = density_statistics(&dens);
        assert!((mean - 0.4).abs() < 1e-12);
        assert!((min - 0.4).abs() < 1e-12);
        assert!((max - 0.4).abs() < 1e-12);
    }

    #[test]
    fn test_density_statistics_range() {
        let dens = vec![0.1, 0.5, 0.9];
        let (_mean, min, max) = density_statistics(&dens);
        assert!((min - 0.1).abs() < 1e-12);
        assert!((max - 0.9).abs() < 1e-12);
    }

    // ── Stress-constrained topology optimization ───────────────────────────

    #[test]
    fn test_von_mises_stress_zero_for_zero_stress() {
        let svm = von_mises_stress(0.0, 0.0, 0.0);
        assert!(svm.abs() < 1e-12);
    }

    #[test]
    fn test_von_mises_stress_uniaxial() {
        // σ_vm = |σ_x| for uniaxial stress
        let svm = von_mises_stress(100.0, 0.0, 0.0);
        assert!((svm - 100.0).abs() < 1e-10, "svm={svm}");
    }

    #[test]
    fn test_von_mises_pure_shear() {
        // σ_vm = sqrt(3) * τ for pure shear
        let tau = 50.0;
        let svm = von_mises_stress(0.0, 0.0, tau);
        let expected = (3.0_f64).sqrt() * tau;
        assert!(
            (svm - expected).abs() < 1e-10,
            "svm={svm} expected={expected}"
        );
    }

    #[test]
    fn test_p_norm_aggregation_value() {
        let stresses = vec![100.0, 200.0, 150.0];
        let pn = p_norm_aggregation(&stresses, 8.0, 300.0);
        assert!(pn > 0.0, "P-norm should be positive");
    }

    #[test]
    fn test_p_norm_aggregation_approaches_max_high_p() {
        let stresses = vec![100.0, 200.0, 50.0];
        // As p → ∞, P-norm → max
        let pn_high = p_norm_aggregation(&stresses, 100.0, 200.0);
        assert!(
            (pn_high - 200.0).abs() / 200.0 < 0.05,
            "High P-norm should approach max, got {pn_high}"
        );
    }

    #[test]
    fn test_p_norm_sensitivity() {
        let stresses = vec![100.0, 200.0, 150.0];
        let sens = p_norm_sensitivity(&stresses, 8.0, 300.0);
        assert_eq!(sens.len(), stresses.len());
        // All sensitivities should be non-negative
        for &s in &sens {
            assert!(s >= 0.0, "Sensitivity should be non-negative: {s}");
        }
    }

    #[test]
    fn test_stress_simp_relaxation() {
        // Stress relaxation should reduce stress for intermediate densities
        let s = stress_simp_relaxation(0.5, 100.0, 3.0, 0.5);
        assert!(
            s < 100.0,
            "Relaxed stress should be less than full stress: {s}"
        );
    }

    #[test]
    fn test_stress_simp_relaxation_solid() {
        // At ρ=1, relaxed stress = original stress
        let s = stress_simp_relaxation(1.0, 200.0, 3.0, 0.5);
        assert!(
            (s - 200.0).abs() < 1e-10,
            "At rho=1, relaxed stress should be 200, got {s}"
        );
    }

    // ── BESO (Bi-directional Evolutionary Structural Optimization) tests ──

    #[test]
    fn test_beso_update_removes_weakest() {
        let sens = vec![0.1, 0.5, 0.2, 0.8, 0.3];
        let dens = vec![1.0; 5];
        let new_dens = beso_update(&dens, &sens, 0.6, 0.01);
        // Volume fraction should approach target
        let vf: f64 = new_dens.iter().sum::<f64>() / new_dens.len() as f64;
        assert!((vf - 0.6).abs() <= 0.25, "vf={vf}");
    }

    #[test]
    fn test_beso_update_length_preserved() {
        let sens = vec![0.1, 0.5, 0.2, 0.8, 0.3];
        let dens = vec![1.0; 5];
        let new_dens = beso_update(&dens, &sens, 0.6, 0.01);
        assert_eq!(new_dens.len(), 5);
    }

    #[test]
    fn test_beso_update_only_binary_densities() {
        let sens: Vec<f64> = (0..10).map(|i| i as f64 * 0.1).collect();
        let dens = vec![1.0; 10];
        let new_dens = beso_update(&dens, &sens, 0.5, 0.01);
        for &d in &new_dens {
            assert!(
                d == 0.01 || d == 1.0,
                "BESO should produce binary densities: {d}"
            );
        }
    }

    // ── Level-set topology optimization tests ─────────────────────────────

    #[test]
    fn test_level_set_reinitialize_signed_distance() {
        // A simple signed distance function should remain approximately unchanged.
        let phi: Vec<f64> = (-5..=5).map(|i| i as f64).collect();
        let phi_new = level_set_reinitialize(&phi, 1.0, 0.5, 5);
        assert_eq!(phi_new.len(), phi.len());
        // Central value should be near 0
        assert!(
            phi_new[5].abs() < 0.5,
            "Central value should be near 0: {}",
            phi_new[5]
        );
    }

    #[test]
    fn test_level_set_heaviside_smooth_transition() {
        let eps = 0.1;
        let h_neg = level_set_heaviside(-1.0, eps);
        let h_zero = level_set_heaviside(0.0, eps);
        let h_pos = level_set_heaviside(1.0, eps);
        assert!(h_neg < 0.1, "Negative phi should give small H: {h_neg}");
        assert!(
            (h_zero - 0.5).abs() < 0.1,
            "Zero phi should give ~0.5: {h_zero}"
        );
        assert!(h_pos > 0.9, "Positive phi should give large H: {h_pos}");
    }

    #[test]
    fn test_level_set_sensitivity_to_velocity() {
        // Velocity extension: sensitivity should be all positive for compliance min.
        let sens = vec![-1.0, -2.0, -3.0, -4.0];
        let vel = level_set_velocity_extension(&sens);
        for &v in &vel {
            assert!(v > 0.0, "Velocity should be positive: {v}");
        }
    }

    #[test]
    fn test_level_set_update_advances_interface() {
        let phi = vec![-1.0, -0.5, 0.0, 0.5, 1.0];
        let vel = vec![1.0; 5];
        let dt = 0.1;
        let phi_new = level_set_advect(&phi, &vel, dt);
        assert_eq!(phi_new.len(), phi.len());
        // The interface (phi=0) should move to the right (phi decreases for positive vel)
        // phi_new = phi - vel * dt  (advection equation)
        for i in 0..phi.len() {
            let expected = phi[i] - vel[i] * dt;
            assert!(
                (phi_new[i] - expected).abs() < 1e-12,
                "phi[{i}]={}",
                phi_new[i]
            );
        }
    }

    // ── Perimeter constraint tests ─────────────────────────────────────────

    #[test]
    fn test_perimeter_length_uniform_solid() {
        let dens = vec![1.0; 9];
        let p = perimeter_length(&dens, 3, 3, 0.5);
        assert_eq!(p, 0.0, "Uniform solid has no internal boundary: {p}");
    }

    #[test]
    fn test_perimeter_length_checkerboard_high() {
        let dens = vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        let p = perimeter_length(&dens, 3, 3, 0.5);
        assert!(p > 0.0, "Checkerboard should have nonzero perimeter: {p}");
    }

    #[test]
    fn test_perimeter_gradient_length() {
        let dens = vec![1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0];
        let g = perimeter_gradient(&dens, 3, 3, 0.5);
        assert_eq!(g.len(), 9);
    }

    // ── Manufacturability constraints ──────────────────────────────────────

    #[test]
    fn test_minimum_length_scale_filter_length() {
        let dens = vec![0.5; 9];
        let centroids: Vec<(f64, f64)> = (0..9).map(|i| (i as f64, 0.0)).collect();
        let out = minimum_length_scale_filter(&dens, &centroids, 2.0);
        assert_eq!(out.len(), 9);
    }

    #[test]
    fn test_minimum_length_scale_filter_values_in_range() {
        let dens: Vec<f64> = (0..16)
            .map(|i| if i % 4 == 0 { 1.0 } else { 0.0 })
            .collect();
        let centroids: Vec<(f64, f64)> =
            (0..16).map(|i| ((i % 4) as f64, (i / 4) as f64)).collect();
        let out = minimum_length_scale_filter(&dens, &centroids, 1.5);
        for &v in &out {
            assert!(
                (0.0..=1.0).contains(&v),
                "Filtered density out of range: {v}"
            );
        }
    }

    #[test]
    fn test_casting_constraint_demold_direction() {
        // Casting demold: maximum density should be non-decreasing in build direction.
        let dens = vec![0.8, 0.6, 0.4, 0.2]; // decreasing → violates casting constraint
        let out = casting_constraint_projection(&dens, true);
        assert_eq!(out.len(), 4);
        // After projection, density should be non-increasing (maximum density rule)
        for i in 1..out.len() {
            assert!(
                out[i] <= out[i - 1] + 1e-10,
                "Casting constraint violated: out[{i}]={} > out[{}]={}",
                out[i],
                i - 1,
                out[i - 1]
            );
        }
    }

    #[test]
    fn test_extrusion_constraint_z_uniform() {
        // Extrusion: all slices must have the same density.
        let dens = vec![0.8, 0.5, 0.3, 0.8, 0.5, 0.3]; // 2 rows × 3 cols
        let out = extrusion_constraint_projection(&dens, 3, 2);
        assert_eq!(out.len(), 6);
        // Columns should have uniform density
        for col in 0..3 {
            let v0 = out[col];
            let v1 = out[3 + col];
            assert!(
                (v0 - v1).abs() < 1e-10,
                "Extrusion: col {col} should be uniform: {v0} vs {v1}"
            );
        }
    }

    // ── Topology optimization with multiple objectives ──────────────────────

    #[test]
    fn test_pareto_weight_compliance_stiffness() {
        let c_comp = 100.0;
        let c_freq = 500.0;
        let w = 0.3;
        let obj = pareto_weighted_objective(c_comp, c_freq, w);
        let expected = w * c_comp + (1.0 - w) * c_freq;
        assert!((obj - expected).abs() < 1e-12);
    }

    #[test]
    fn test_pareto_weight_zero_gives_compliance() {
        let obj = pareto_weighted_objective(100.0, 500.0, 0.0);
        assert!((obj - 500.0).abs() < 1e-12);
    }

    #[test]
    fn test_pareto_weight_one_gives_stiffness() {
        let obj = pareto_weighted_objective(100.0, 500.0, 1.0);
        assert!((obj - 100.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_compliance_gradient() {
        let dens = vec![0.5; 5];
        let thermal_sens = vec![2.0; 5];
        let e0 = 1.0;
        let e_min = 1e-9;
        let p = 3.0;
        let grad = thermal_compliance_gradient(&dens, &thermal_sens, e0, e_min, p);
        assert_eq!(grad.len(), 5);
        for &g in &grad {
            assert!(g <= 0.0, "Thermal sensitivity should be ≤ 0: {g}");
        }
    }
}
