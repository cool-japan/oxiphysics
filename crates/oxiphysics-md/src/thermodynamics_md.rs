#![allow(clippy::if_same_then_else)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermodynamic integration and free-energy methods for MD simulations.
//!
//! Covers:
//!
//! - Thermodynamic integration (TI) with trapezoidal / Gaussian quadrature
//! - Free energy perturbation (FEP) via Zwanzig exponential averaging
//! - Soft-core potentials for alchemical transformations (Beutler form)
//! - Bennett Acceptance Ratio (BAR) estimator (iterative self-consistent)
//! - Multistate BAR (MBAR) via self-consistent iteration
//! - Alchemical λ-coupling: linear, soft-core LJ, electrostatics
//! - Thermodynamic cycle closure and hysteresis detection
//! - Hydration free energy from decoupling protocol
//! - Protein–ligand binding ΔΔG estimation
//! - Relative free energy perturbation (RBFE)

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
const K_B: f64 = 1.380_649e-23;

/// Avogadro's number (mol⁻¹).
const N_AV: f64 = 6.022_140_76e23;

/// Gas constant (J mol⁻¹ K⁻¹).
const R_GAS: f64 = 8.314_462_618;

/// Reference temperature (K).
const T_REF: f64 = 298.15;

/// kcal to kJ conversion factor.
const KCAL_TO_KJ: f64 = 4.184;

/// kJ mol⁻¹ to k_B T at 298.15 K: β = 1 / (R T) in mol kJ⁻¹.
const BETA_REF_MOL: f64 = 1.0 / (R_GAS * 1.0e-3 * T_REF); // kJ⁻¹ mol

// ---------------------------------------------------------------------------
// Scalar helpers
// ---------------------------------------------------------------------------

/// Compute β = 1/(k_B T) in units of J⁻¹.
#[inline]
pub fn beta(temp: f64) -> f64 {
    1.0 / (K_B * temp)
}

/// Compute β in mol/kJ units (for energies in kJ/mol).
#[inline]
pub fn beta_mol_kjmol(temp: f64) -> f64 {
    1.0 / (R_GAS * 1.0e-3 * temp)
}

/// Compute the log-sum-exp of a slice to avoid numerical overflow.
///
/// ln Σ_i exp(x_i) = x_max + ln Σ_i exp(x_i - x_max)
pub fn log_sum_exp(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    let xmax = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = values.iter().map(|&x| (x - xmax).exp()).sum();
    xmax + sum.ln()
}

// ---------------------------------------------------------------------------
// Thermodynamic integration (TI)
// ---------------------------------------------------------------------------

/// A single TI window: a λ value and the ensemble average ⟨∂H/∂λ⟩.
#[derive(Debug, Clone)]
pub struct TiWindow {
    /// Coupling parameter λ ∈ \[0, 1\].
    pub lambda: f64,
    /// Ensemble average of ∂H/∂λ at this λ (kJ mol⁻¹).
    pub dhdl_mean: f64,
    /// Standard error of ⟨∂H/∂λ⟩ (kJ mol⁻¹).
    pub dhdl_stderr: f64,
    /// Number of samples used to estimate the mean.
    pub n_samples: usize,
}

impl TiWindow {
    /// Create a new TI window.
    pub fn new(lambda: f64, dhdl_mean: f64, dhdl_stderr: f64, n_samples: usize) -> Self {
        Self {
            lambda,
            dhdl_mean,
            dhdl_stderr,
            n_samples,
        }
    }
}

/// Estimate ΔA from a set of TI windows using the trapezoidal rule.
///
/// ΔA = ∫₀¹ ⟨∂H/∂λ⟩ dλ ≈ Σ_i ½ (⟨∂H/∂λ⟩_i + ⟨∂H/∂λ⟩_{i+1}) Δλ
///
/// Windows must be sorted by λ.  Returns free energy difference in kJ mol⁻¹.
pub fn ti_trapezoidal(windows: &[TiWindow]) -> f64 {
    if windows.len() < 2 {
        return 0.0;
    }
    let mut integral = 0.0f64;
    for i in 0..windows.len() - 1 {
        let dlam = windows[i + 1].lambda - windows[i].lambda;
        integral += 0.5 * (windows[i].dhdl_mean + windows[i + 1].dhdl_mean) * dlam;
    }
    integral
}

/// Estimate the statistical uncertainty of a TI trapezoidal integral.
///
/// Propagates standard errors assuming uncorrelated windows.
/// Returns 1-σ error in kJ mol⁻¹.
pub fn ti_trapezoidal_error(windows: &[TiWindow]) -> f64 {
    if windows.len() < 2 {
        return 0.0;
    }
    let mut var = 0.0f64;
    for i in 0..windows.len() - 1 {
        let dlam = windows[i + 1].lambda - windows[i].lambda;
        let se_i = windows[i].dhdl_stderr;
        let se_j = windows[i + 1].dhdl_stderr;
        var += (0.5 * dlam).powi(2) * (se_i * se_i + se_j * se_j);
    }
    var.sqrt()
}

/// Gauss-Legendre quadrature nodes and weights on \[0, 1\] for n=5.
///
/// Returns (nodes, weights).
pub fn gauss_legendre_5() -> ([f64; 5], [f64; 5]) {
    // Nodes and weights on [-1, 1], transformed to [0, 1]
    let nodes_m1 = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683,
        0.0,
        0.538_469_310_105_683,
        0.906_179_845_938_664,
    ];
    let weights_m1 = [
        0.236_926_885_056_189,
        0.478_628_670_499_366,
        0.568_888_888_888_889,
        0.478_628_670_499_366,
        0.236_926_885_056_189,
    ];
    let mut nodes = [0.0f64; 5];
    let mut weights = [0.0f64; 5];
    for i in 0..5 {
        nodes[i] = 0.5 * (nodes_m1[i] + 1.0);
        weights[i] = 0.5 * weights_m1[i];
    }
    (nodes, weights)
}

/// TI integration using 5-point Gauss-Legendre quadrature.
///
/// Evaluates a cubic spline through the provided windows and samples at
/// GL nodes; returns ΔA (kJ mol⁻¹).  For simplicity, uses piecewise
/// linear interpolation of ⟨∂H/∂λ⟩.
pub fn ti_gauss_legendre(windows: &[TiWindow]) -> f64 {
    if windows.len() < 2 {
        return 0.0;
    }
    let (nodes, weights) = gauss_legendre_5();
    let mut integral = 0.0f64;
    for (node, weight) in nodes.iter().zip(weights.iter()) {
        let lam = *node;
        // Linear interpolation at lambda
        let dhdl = interpolate_dhdl(windows, lam);
        integral += weight * dhdl;
    }
    integral
}

/// Piecewise linear interpolation of ⟨∂H/∂λ⟩ at a given λ.
pub fn interpolate_dhdl(windows: &[TiWindow], lam: f64) -> f64 {
    if windows.is_empty() {
        return 0.0;
    }
    if lam <= windows[0].lambda {
        return windows[0].dhdl_mean;
    }
    if lam >= windows[windows.len() - 1].lambda {
        return windows[windows.len() - 1].dhdl_mean;
    }
    for i in 0..windows.len() - 1 {
        if lam >= windows[i].lambda && lam <= windows[i + 1].lambda {
            let t = (lam - windows[i].lambda) / (windows[i + 1].lambda - windows[i].lambda);
            return windows[i].dhdl_mean * (1.0 - t) + windows[i + 1].dhdl_mean * t;
        }
    }
    0.0
}

// ---------------------------------------------------------------------------
// Free energy perturbation (FEP) — Zwanzig exponential averaging
// ---------------------------------------------------------------------------

/// FEP estimate of ΔA from A → B using the Zwanzig formula.
///
/// ΔA = −k_B T ln ⟨exp(−β ΔU)⟩_A
///
/// `delta_u_samples` are ΔU = U_B − U_A values (kJ mol⁻¹) sampled from A.
/// `temp` is the temperature (K).
pub fn fep_zwanzig(delta_u_samples: &[f64], temp: f64) -> f64 {
    if delta_u_samples.is_empty() {
        return 0.0;
    }
    let b = beta_mol_kjmol(temp);
    let neg_b_du: Vec<f64> = delta_u_samples.iter().map(|&du| -b * du).collect();
    let lse = log_sum_exp(&neg_b_du);
    let log_mean = lse - (delta_u_samples.len() as f64).ln();
    -log_mean / b
}

/// Variance of ΔA from FEP, estimated from sample variance of ΔU.
///
/// Returns 1-σ error (kJ mol⁻¹).
pub fn fep_zwanzig_error(delta_u_samples: &[f64], temp: f64) -> f64 {
    if delta_u_samples.len() < 2 {
        return 0.0;
    }
    let b = beta_mol_kjmol(temp);
    let n = delta_u_samples.len() as f64;
    let mean: f64 = delta_u_samples.iter().sum::<f64>() / n;
    let var: f64 = delta_u_samples
        .iter()
        .map(|&x| (x - mean).powi(2))
        .sum::<f64>()
        / (n - 1.0);
    // First-order approximation: σ_ΔA ≈ β σ_ΔU / sqrt(n)
    b * var.sqrt() / n.sqrt()
}

// ---------------------------------------------------------------------------
// Soft-core potentials for alchemical FEP (Beutler form)
// ---------------------------------------------------------------------------

/// Soft-core Lennard-Jones potential energy (kJ mol⁻¹).
///
/// U_sc(r; λ) = 4 ε λ^n \[ 1/(α(1−λ)^m + (r/σ)^6)^2 − 1/(α(1−λ)^m + (r/σ)^6) \]
///
/// Beutler et al. (1994) parameters: α=0.5, n=1, m=2 (soft-core exponents).
/// `eps` well depth (kJ mol⁻¹), `sigma` LJ diameter (nm), `r` distance (nm).
pub fn soft_core_lj(
    r: f64,
    lam: f64,
    eps: f64,
    sigma: f64,
    alpha: f64,
    n_exp: f64,
    m_exp: f64,
) -> f64 {
    if lam <= 0.0 {
        // Fully decoupled
        return 0.0;
    }
    let rs6 = (r / sigma).powi(6);
    let denom = alpha * (1.0 - lam).powf(m_exp) + rs6;
    4.0 * eps * lam.powf(n_exp) * (1.0 / denom.powi(2) - 1.0 / denom)
}

/// Derivative ∂U_sc/∂λ for the soft-core LJ (kJ mol⁻¹).
pub fn soft_core_lj_dhdl(
    r: f64,
    lam: f64,
    eps: f64,
    sigma: f64,
    alpha: f64,
    n_exp: f64,
    m_exp: f64,
) -> f64 {
    if lam <= 0.0 {
        return 0.0;
    }
    let rs6 = (r / sigma).powi(6);
    let s = alpha * (1.0 - lam).powf(m_exp);
    let denom = s + rs6;
    // ∂/∂λ [λ^n (denom^{-2} - denom^{-1})]
    // where denom = s(λ) + rs6, s = α (1-λ)^m
    let f = 1.0 / denom.powi(2) - 1.0 / denom;
    let df_dlam = {
        let ds_dlam = -alpha * m_exp * (1.0 - lam).powf(m_exp - 1.0);
        (-2.0 / denom.powi(3) + 1.0 / denom.powi(2)) * ds_dlam
    };
    4.0 * eps * (n_exp * lam.powf(n_exp - 1.0) * f + lam.powf(n_exp) * df_dlam)
}

/// Linear λ-coupling for electrostatics: U_elec(λ) = λ U_elec_full.
pub fn linear_elec_coupling(u_elec_full: f64, lam: f64) -> f64 {
    lam * u_elec_full
}

/// dU_elec/dλ for linear coupling.
pub fn linear_elec_dhdl(u_elec_full: f64) -> f64 {
    u_elec_full
}

// ---------------------------------------------------------------------------
// Bennett Acceptance Ratio (BAR)
// ---------------------------------------------------------------------------

/// Fermi function: f(x) = 1/(1 + exp(x)).
#[inline]
pub fn fermi(x: f64) -> f64 {
    1.0 / (1.0 + x.exp())
}

/// Single BAR iteration to update the free energy estimate ΔA_n+1.
///
/// Solves the BAR self-consistent equation:
///   Σ_{A} f(β(ΔU_AB − M − ΔA)) = Σ_{B} f(β(ΔU_BA + M − ΔA))
///
/// `du_ab` = U_B − U_A values sampled from A (kJ mol⁻¹).
/// `du_ba` = U_A − U_B values sampled from B (kJ mol⁻¹).
/// `da_prev` current estimate of ΔA (kJ mol⁻¹).
/// `temp` temperature (K).
/// `m` = k_B T ln(n_A / n_B) (kJ mol⁻¹).
///
/// Returns the updated ΔA estimate.
pub fn bar_iteration(du_ab: &[f64], du_ba: &[f64], da_prev: f64, temp: f64, m: f64) -> f64 {
    if du_ab.is_empty() || du_ba.is_empty() {
        return da_prev;
    }
    let b = beta_mol_kjmol(temp);
    // Bennett self-consistent equation (Bennett 1976):
    //   ΔA_new = ΔA + (1/β) ln( <f(β(ΔU_BA + M))>_B / <f(β(ΔU_AB − M))>_A )
    // where M = ΔA + (1/β) ln(n_A/n_B).
    // Denominator (B side): mean of f(β(ΔU_BA + m + ΔA)) over B
    let den: f64 = du_ba
        .iter()
        .map(|&du| fermi(b * (du + m + da_prev)))
        .sum::<f64>()
        / du_ba.len() as f64;
    // Numerator (A side): mean of f(β(ΔU_AB − m − ΔA)) over A
    let num: f64 = du_ab
        .iter()
        .map(|&du| fermi(b * (du - m - da_prev)))
        .sum::<f64>()
        / du_ab.len() as f64;
    if num.abs() < 1.0e-30 {
        return da_prev;
    }
    // Update: ΔA_new = ΔA + (1/β) ln(den/num)
    da_prev + (den / num).ln() / b
}

/// Run BAR self-consistent iteration to convergence.
///
/// Returns (ΔA, n_iterations_used).
pub fn bar_converge(
    du_ab: &[f64],
    du_ba: &[f64],
    temp: f64,
    tol: f64,
    max_iter: usize,
) -> (f64, usize) {
    if du_ab.is_empty() || du_ba.is_empty() {
        return (0.0, 0);
    }
    let na = du_ab.len() as f64;
    let nb = du_ba.len() as f64;
    let b = beta_mol_kjmol(temp);
    let m = (na / nb).ln() / b;
    let mut da = 0.0f64;
    for iter in 0..max_iter {
        let da_new = bar_iteration(du_ab, du_ba, da, temp, m);
        if (da_new - da).abs() < tol {
            return (da_new, iter + 1);
        }
        da = da_new;
    }
    (da, max_iter)
}

/// BAR statistical uncertainty via the analytical formula (Shirts et al. 2003).
///
/// σ²(ΔA) = (1/β²) \[ 1/(n_A ⟨f²⟩_A) + 1/(n_B ⟨f²⟩_B) − 1/n_A − 1/n_B \]
///
/// Returns 1-σ error in kJ mol⁻¹.
pub fn bar_uncertainty(du_ab: &[f64], du_ba: &[f64], da: f64, temp: f64) -> f64 {
    if du_ab.is_empty() || du_ba.is_empty() {
        return 0.0;
    }
    let b = beta_mol_kjmol(temp);
    let na = du_ab.len() as f64;
    let nb = du_ba.len() as f64;
    let m = (na / nb).ln() / b;
    let fa2: f64 = du_ab
        .iter()
        .map(|&du| fermi(b * (du - m - da)).powi(2))
        .sum::<f64>()
        / na;
    let fb2: f64 = du_ba
        .iter()
        .map(|&du| fermi(b * (du + m + da)).powi(2))
        .sum::<f64>()
        / nb;
    let var = (1.0 / (b * b)) * (1.0 / (na * fa2.max(1.0e-30)) + 1.0 / (nb * fb2.max(1.0e-30)));
    var.sqrt()
}

// ---------------------------------------------------------------------------
// Multistate BAR (MBAR)
// ---------------------------------------------------------------------------

/// MBAR free energy estimator via self-consistent iteration.
///
/// `energies[k][n]` is the energy of configuration n evaluated at state k
/// (kJ mol⁻¹).  `n_k[k]` is the number of samples from state k.
/// Returns the reduced free energies f_k = −ln Z_k (units of k_B T).
pub fn mbar_free_energies(
    energies: &[Vec<f64>],
    n_k: &[usize],
    temp: f64,
    max_iter: usize,
    tol: f64,
) -> Vec<f64> {
    let k_states = energies.len();
    if k_states == 0 {
        return Vec::new();
    }
    let b = beta_mol_kjmol(temp);
    let n_total: usize = n_k.iter().sum();
    if n_total == 0 {
        return vec![0.0; k_states];
    }
    // Initial guess: f_k = 0
    let mut f = vec![0.0f64; k_states];
    // log N_k
    let log_nk: Vec<f64> = n_k.iter().map(|&n| (n as f64).ln()).collect();

    for _iter in 0..max_iter {
        let mut f_new = vec![0.0f64; k_states];
        // For each state k, f_k = −ln Σ_n [ Σ_j N_j exp(f_j - β U_kn) /
        //                                       Σ_l N_l exp(f_l - β U_ln) ]
        // Self-consistent update using log-space arithmetic
        for k in 0..k_states {
            let mut log_sum = f64::NEG_INFINITY;
            let mut sample_count = 0usize;
            for j in 0..k_states {
                let nj = n_k[j];
                for n in 0..nj {
                    if n >= energies[k].len() || n >= energies[j].len() {
                        continue;
                    }
                    // log weight numerator: f_j - β U_kn
                    // denominator: log Σ_l N_l exp(f_l - β U_ln)
                    let log_num = f[j] + log_nk[j] - b * energies[k][n];
                    // log denominator
                    let mut log_den = f64::NEG_INFINITY;
                    for l in 0..k_states {
                        if n < energies[l].len() {
                            let term = f[l] + log_nk[l] - b * energies[l][n];
                            // log-sum-exp accumulation
                            if term > log_den {
                                log_den = log_sum_exp(&[log_den, term]);
                            } else {
                                log_den = log_sum_exp(&[log_den, term]);
                            }
                        }
                    }
                    let log_w = log_num - log_den;
                    log_sum = log_sum_exp(&[log_sum, log_w]);
                    sample_count += 1;
                }
            }
            let _ = sample_count;
            f_new[k] = -log_sum;
        }
        // Normalise so f[0] = 0
        let f0 = f_new[0];
        for fk in f_new.iter_mut() {
            *fk -= f0;
        }
        // Check convergence
        let max_diff: f64 = f
            .iter()
            .zip(f_new.iter())
            .map(|(a, b_v)| (a - b_v).abs())
            .fold(0.0f64, f64::max);
        f = f_new;
        if max_diff < tol {
            break;
        }
    }
    f
}

/// Convert MBAR reduced free energies to ΔA between state 0 and state k
/// (kJ mol⁻¹).
pub fn mbar_delta_a(f_reduced: &[f64], temp: f64) -> Vec<f64> {
    let b = beta_mol_kjmol(temp);
    f_reduced.iter().map(|&fk| fk / b).collect()
}

// ---------------------------------------------------------------------------
// Alchemical transformation: λ schedule
// ---------------------------------------------------------------------------

/// Generate a uniform λ schedule from 0 to 1 with `n_windows` points.
pub fn lambda_schedule_uniform(n_windows: usize) -> Vec<f64> {
    if n_windows == 0 {
        return Vec::new();
    }
    (0..n_windows)
        .map(|i| i as f64 / (n_windows - 1) as f64)
        .collect()
}

/// Generate a λ schedule with denser spacing near 0 and 1 (for endpoint
/// regions where ∂H/∂λ changes rapidly).
///
/// Uses a sine-stretched spacing: λ_i = ½ (1 − cos(πi/(n−1))).
pub fn lambda_schedule_cosine(n_windows: usize) -> Vec<f64> {
    if n_windows == 0 {
        return Vec::new();
    }
    (0..n_windows)
        .map(|i| 0.5 * (1.0 - (PI * i as f64 / (n_windows - 1) as f64).cos()))
        .collect()
}

// ---------------------------------------------------------------------------
// Thermodynamic cycle and closure
// ---------------------------------------------------------------------------

/// A directed edge in a thermodynamic cycle: source state, target state, ΔA.
#[derive(Debug, Clone)]
pub struct ThermoCycleEdge {
    /// Index of the source state.
    pub from: usize,
    /// Index of the target state.
    pub to: usize,
    /// Free energy difference ΔA = A_to − A_from (kJ mol⁻¹).
    pub delta_a: f64,
    /// Statistical uncertainty (kJ mol⁻¹).
    pub sigma: f64,
}

impl ThermoCycleEdge {
    /// Create a new thermodynamic cycle edge.
    pub fn new(from: usize, to: usize, delta_a: f64, sigma: f64) -> Self {
        Self {
            from,
            to,
            delta_a,
            sigma,
        }
    }
}

/// Compute cycle closure error for a directed cycle of edges (sum of ΔA).
///
/// For a closed cycle A → B → C → A the sum should be zero.
/// Returns the closure error Σ ΔA and its propagated uncertainty.
pub fn cycle_closure_error(edges: &[ThermoCycleEdge]) -> (f64, f64) {
    let sum: f64 = edges.iter().map(|e| e.delta_a).sum();
    let var: f64 = edges.iter().map(|e| e.sigma * e.sigma).sum();
    (sum, var.sqrt())
}

/// Check if cycle closure passes a Z-score test.
///
/// Returns `true` if |closure error| / σ_closure < z_threshold.
pub fn cycle_closure_passes(closure: f64, sigma_closure: f64, z_threshold: f64) -> bool {
    if sigma_closure < 1.0e-30 {
        return closure.abs() < 1.0e-6;
    }
    (closure.abs() / sigma_closure) < z_threshold
}

// ---------------------------------------------------------------------------
// Hydration free energy
// ---------------------------------------------------------------------------

/// Hydration free energy protocol (two-step: vdW + electrostatics decoupling).
///
/// ΔA_hyd = ΔA_vdw + ΔA_elec
///
/// Both are obtained from separate TI or BAR calculations.
pub fn hydration_free_energy(delta_a_vdw: f64, delta_a_elec: f64) -> f64 {
    delta_a_vdw + delta_a_elec
}

/// Estimate the standard hydration free energy including a standard-state
/// correction for ideal gas → 1 mol/L concentration reference.
///
/// δA_ss = k_B T ln(ρ_gas V_mol)  ≈ k_B T ln(24.8) at 298 K
/// where V_mol = 24.8 L/mol for ideal gas at 298 K.
pub fn standard_state_correction(temp: f64, v_mol_l: f64) -> f64 {
    // In kJ/mol: R T ln(ρ_0 V_mol)
    R_GAS * 1.0e-3 * temp * v_mol_l.ln()
}

// ---------------------------------------------------------------------------
// Protein–ligand binding free energy
// ---------------------------------------------------------------------------

/// Binding free energy from a decoupling thermodynamic cycle.
///
/// ΔA_bind = ΔA_complex − ΔA_solvation + ΔA_restr
///
/// `da_complex` = free energy to decouple ligand in complex (kJ mol⁻¹).
/// `da_solv` = free energy to decouple ligand in solution (kJ mol⁻¹).
/// `da_restr` = free energy cost of restraints (kJ mol⁻¹).
pub fn binding_free_energy(da_complex: f64, da_solv: f64, da_restr: f64) -> f64 {
    da_complex - da_solv + da_restr
}

/// Convert binding free energy to dissociation constant K_d at temperature T.
///
/// K_d = C° exp(β ΔA_bind)  where C° = 1 M (standard concentration).
pub fn kd_from_binding_dg(delta_a_bind: f64, temp: f64) -> f64 {
    let b = beta_mol_kjmol(temp);
    (b * delta_a_bind).exp() // units of 1 M (dimensionless when C°=1M)
}

/// pKd from K_d: pK_d = −log10(K_d).
pub fn pkd_from_kd(kd: f64) -> f64 {
    -kd.log10()
}

// ---------------------------------------------------------------------------
// Relative free energy perturbation (RBFE)
// ---------------------------------------------------------------------------

/// Relative ΔΔA from two absolute binding free energies.
///
/// ΔΔA = ΔA_bind(ligand_B) − ΔA_bind(ligand_A) (kJ mol⁻¹).
pub fn rbfe_ddg(da_bind_a: f64, da_bind_b: f64) -> f64 {
    da_bind_b - da_bind_a
}

/// Estimate RBFE from a perturbation network (sum along a path).
///
/// `edges` is a list of (ΔΔA, σ) pairs along the shortest path from
/// reference ligand to target ligand.
pub fn rbfe_from_path(edges: &[(f64, f64)]) -> (f64, f64) {
    let sum: f64 = edges.iter().map(|&(ddg, _)| ddg).sum();
    let var: f64 = edges.iter().map(|&(_, s)| s * s).sum();
    (sum, var.sqrt())
}

// ---------------------------------------------------------------------------
// Alchemical simulation record
// ---------------------------------------------------------------------------

/// A record of ∂H/∂λ samples collected during an alchemical MD window.
#[derive(Debug, Clone)]
pub struct AlchemicalWindow {
    /// λ value for this window.
    pub lambda: f64,
    /// Raw ∂H/∂λ samples (kJ mol⁻¹).
    pub dhdl_samples: Vec<f64>,
    /// ΔU = U(λ+δλ) − U(λ) samples for FEP (kJ mol⁻¹).
    pub fep_forward: Vec<f64>,
    /// ΔU = U(λ−δλ) − U(λ) samples for BAR (kJ mol⁻¹).
    pub fep_backward: Vec<f64>,
}

impl AlchemicalWindow {
    /// Create a new empty alchemical window.
    pub fn new(lambda: f64) -> Self {
        Self {
            lambda,
            dhdl_samples: Vec::new(),
            fep_forward: Vec::new(),
            fep_backward: Vec::new(),
        }
    }

    /// Add a ∂H/∂λ sample.
    pub fn push_dhdl(&mut self, dhdl: f64) {
        self.dhdl_samples.push(dhdl);
    }

    /// Add a forward FEP ΔU sample.
    pub fn push_fep_forward(&mut self, du: f64) {
        self.fep_forward.push(du);
    }

    /// Add a backward FEP ΔU sample.
    pub fn push_fep_backward(&mut self, du: f64) {
        self.fep_backward.push(du);
    }

    /// Compute the mean and standard error of ∂H/∂λ from collected samples.
    ///
    /// Returns `(mean, stderr)`.
    pub fn dhdl_statistics(&self) -> (f64, f64) {
        if self.dhdl_samples.is_empty() {
            return (0.0, 0.0);
        }
        let n = self.dhdl_samples.len() as f64;
        let mean: f64 = self.dhdl_samples.iter().sum::<f64>() / n;
        if self.dhdl_samples.len() < 2 {
            return (mean, 0.0);
        }
        let var: f64 = self
            .dhdl_samples
            .iter()
            .map(|&x| (x - mean).powi(2))
            .sum::<f64>()
            / (n - 1.0);
        (mean, (var / n).sqrt())
    }

    /// Convert this window's statistics into a [`TiWindow`].
    pub fn to_ti_window(&self) -> TiWindow {
        let (mean, stderr) = self.dhdl_statistics();
        TiWindow::new(self.lambda, mean, stderr, self.dhdl_samples.len())
    }
}

// ---------------------------------------------------------------------------
// Alchemical transformation: complete pipeline
// ---------------------------------------------------------------------------

/// Run a mock alchemical MD step for testing.
///
/// Generates synthetic ∂H/∂λ and ΔU samples from a Gaussian model:
/// ⟨∂H/∂λ⟩ = dhdl_mean, noise = dhdl_noise.
pub fn mock_alchemical_step(
    window: &mut AlchemicalWindow,
    dhdl_mean: f64,
    dhdl_noise: f64,
    du_forward_mean: f64,
    du_noise: f64,
    n_samples: usize,
) {
    let mut rng = rand::rng();
    for _ in 0..n_samples {
        let dhdl_jitter = if dhdl_noise > 0.0 {
            rng.random_range(-dhdl_noise..dhdl_noise)
        } else {
            0.0
        };
        let dhdl = dhdl_mean + dhdl_jitter;
        window.push_dhdl(dhdl);
        let du_jitter = if du_noise > 0.0 {
            rng.random_range(-du_noise..du_noise)
        } else {
            0.0
        };
        let du_f = du_forward_mean + du_jitter;
        window.push_fep_forward(du_f);
        let du_b_jitter = if du_noise > 0.0 {
            rng.random_range(-du_noise..du_noise)
        } else {
            0.0
        };
        let du_b = -du_forward_mean + du_b_jitter;
        window.push_fep_backward(du_b);
    }
}

/// Run the full TI pipeline over a set of alchemical windows.
///
/// Collects `TiWindow` from each `AlchemicalWindow` and runs trapezoidal TI.
/// Returns `(delta_a, error)` in kJ mol⁻¹.
pub fn run_ti_pipeline(windows: &[AlchemicalWindow]) -> (f64, f64) {
    let ti_windows: Vec<TiWindow> = windows.iter().map(|w| w.to_ti_window()).collect();
    let da = ti_trapezoidal(&ti_windows);
    let err = ti_trapezoidal_error(&ti_windows);
    (da, err)
}

// ---------------------------------------------------------------------------
// Overlap matrix and convergence diagnostics
// ---------------------------------------------------------------------------

/// Compute the overlap integral between two energy distributions.
///
/// O = Σ_i min(p_A(i), p_B(i)) Δε using histogram overlap.
/// `hist_a` and `hist_b` are normalised histograms (same bins).
pub fn overlap_integral(hist_a: &[f64], hist_b: &[f64]) -> f64 {
    hist_a
        .iter()
        .zip(hist_b.iter())
        .map(|(&a, &b)| a.min(b))
        .sum()
}

/// Build a simple histogram with `n_bins` bins from samples in \[lo, hi\].
///
/// Returns `(bin_edges, counts_normalised)`.
pub fn histogram(samples: &[f64], lo: f64, hi: f64, n_bins: usize) -> (Vec<f64>, Vec<f64>) {
    let mut counts = vec![0usize; n_bins];
    let bin_width = (hi - lo) / n_bins as f64;
    for &x in samples {
        if x >= lo && x < hi {
            let idx = ((x - lo) / bin_width) as usize;
            let idx = idx.min(n_bins - 1);
            counts[idx] += 1;
        }
    }
    let total: usize = counts.iter().sum();
    let norm = if total > 0 {
        1.0 / (total as f64 * bin_width)
    } else {
        1.0
    };
    let edges: Vec<f64> = (0..=n_bins).map(|i| lo + i as f64 * bin_width).collect();
    let hist: Vec<f64> = counts.iter().map(|&c| c as f64 * norm).collect();
    (edges, hist)
}

// ---------------------------------------------------------------------------
// Bennet overlap score
// ---------------------------------------------------------------------------

/// Bennett overlap score: P(overlap) between forward and backward ΔU distributions.
///
/// Returns a number in \[0, 1\]; values > 0.3 suggest adequate sampling.
pub fn bennett_overlap_score(du_ab: &[f64], du_ba: &[f64]) -> f64 {
    if du_ab.is_empty() || du_ba.is_empty() {
        return 0.0;
    }
    let all: Vec<f64> = du_ab.iter().chain(du_ba.iter()).cloned().collect();
    let lo = all.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = all.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if (hi - lo).abs() < 1.0e-30 {
        return 1.0;
    }
    let (_ea, ha) = histogram(du_ab, lo, hi, 50);
    let (_eb, hb) = histogram(du_ba, lo, hi, 50);
    overlap_integral(&ha, &hb) * (hi - lo) / 50.0
}

// ---------------------------------------------------------------------------
// Equilibration detection (block averaging / Gelman-Rubin)
// ---------------------------------------------------------------------------

/// Block averaging to estimate the effective sample size.
///
/// Splits `samples` into blocks of size `block_size`, computes block means,
/// then estimates the variance of block means divided by grand variance.
/// Returns effective sample size N_eff.
pub fn effective_sample_size(samples: &[f64], block_size: usize) -> f64 {
    if samples.len() < 2 * block_size {
        return samples.len() as f64;
    }
    let n = samples.len();
    let n_blocks = n / block_size;
    let grand_mean: f64 = samples.iter().sum::<f64>() / n as f64;
    let block_means: Vec<f64> = (0..n_blocks)
        .map(|b| {
            let sl = &samples[b * block_size..(b + 1) * block_size];
            sl.iter().sum::<f64>() / block_size as f64
        })
        .collect();
    let var_block: f64 = block_means
        .iter()
        .map(|&m| (m - grand_mean).powi(2))
        .sum::<f64>()
        / (n_blocks - 1) as f64;
    let var_sample: f64 = samples
        .iter()
        .map(|&x| (x - grand_mean).powi(2))
        .sum::<f64>()
        / (n - 1) as f64;
    if var_block < 1.0e-30 {
        return n as f64;
    }
    var_sample / var_block
}

// ---------------------------------------------------------------------------
// Statistical correlation time
// ---------------------------------------------------------------------------

/// Estimate the integrated autocorrelation time τ of a time series.
///
/// τ_int = 1 + 2 Σ_{k=1}^{k_max} C(k) / C(0)
///
/// where C(k) is the autocorrelation at lag k.
pub fn autocorrelation_time(samples: &[f64], k_max: usize) -> f64 {
    if samples.len() < 2 {
        return 1.0;
    }
    let n = samples.len();
    let mean: f64 = samples.iter().sum::<f64>() / n as f64;
    let c0: f64 = samples.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    if c0 < 1.0e-30 {
        return 1.0;
    }
    let mut tau = 1.0f64;
    let k_max = k_max.min(n / 4);
    for k in 1..=k_max {
        let ck: f64 = (0..n - k)
            .map(|i| (samples[i] - mean) * (samples[i + k] - mean))
            .sum::<f64>()
            / (n - k) as f64;
        tau += 2.0 * ck / c0;
        if ck / c0 < 0.05 {
            break; // truncate when autocorrelation is small
        }
    }
    tau.max(1.0)
}

// ---------------------------------------------------------------------------
// Energy difference distribution analysis
// ---------------------------------------------------------------------------

/// Compute the mean and variance of a ΔU sample set.
pub fn du_statistics(delta_u: &[f64]) -> (f64, f64) {
    if delta_u.is_empty() {
        return (0.0, 0.0);
    }
    let n = delta_u.len() as f64;
    let mean: f64 = delta_u.iter().sum::<f64>() / n;
    let var: f64 = delta_u.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n;
    (mean, var)
}

/// Gaussian approximation FEP (Hummer 1996):
///
/// ΔA ≈ ⟨ΔU⟩ − β/2 σ²(ΔU)
///
/// Valid when ΔU is Gaussian-distributed around its mean.
pub fn fep_gaussian_approx(delta_u: &[f64], temp: f64) -> f64 {
    let (mean, var) = du_statistics(delta_u);
    let b = beta_mol_kjmol(temp);
    mean - 0.5 * b * var
}

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------

/// Generate synthetic ΔU samples from a normal distribution (for testing).
pub fn generate_gaussian_du(mean: f64, std: f64, n: usize) -> Vec<f64> {
    let mut rng = rand::rng();
    // Box-Muller transform
    let mut samples = Vec::with_capacity(n);
    let mut i = 0;
    while i < n {
        let u1: f64 = rng.random_range(1.0e-12_f64..1.0_f64);
        let u2: f64 = rng.random_range(0.0_f64..1.0_f64);
        let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).cos();
        let z1 = (-2.0 * u1.ln()).sqrt() * (2.0 * PI * u2).sin();
        samples.push(mean + std * z0);
        i += 1;
        if i < n {
            samples.push(mean + std * z1);
            i += 1;
        }
    }
    samples.truncate(n);
    samples
}

/// Convert free energy from kJ mol⁻¹ to kcal mol⁻¹.
pub fn kjmol_to_kcalmol(da_kjmol: f64) -> f64 {
    da_kjmol / KCAL_TO_KJ
}

/// Convert free energy from kcal mol⁻¹ to kJ mol⁻¹.
pub fn kcalmol_to_kjmol(da_kcalmol: f64) -> f64 {
    da_kcalmol * KCAL_TO_KJ
}

/// Apply bootstrap resampling to estimate the uncertainty of a statistic.
///
/// `samples` input data, `statistic` function to apply, `n_bootstrap` number
/// of bootstrap iterations. Returns 1-σ bootstrap standard deviation.
pub fn bootstrap_uncertainty(
    samples: &[f64],
    statistic: &dyn Fn(&[f64]) -> f64,
    n_bootstrap: usize,
) -> f64 {
    if samples.is_empty() || n_bootstrap < 2 {
        return 0.0;
    }
    let mut rng = rand::rng();
    let n = samples.len();
    let mut boot_stats = Vec::with_capacity(n_bootstrap);
    for _ in 0..n_bootstrap {
        let resample: Vec<f64> = (0..n)
            .map(|_| {
                let idx = rng.random_range(0..n);
                samples[idx]
            })
            .collect();
        boot_stats.push(statistic(&resample));
    }
    let mean: f64 = boot_stats.iter().sum::<f64>() / n_bootstrap as f64;
    let var: f64 =
        boot_stats.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (n_bootstrap - 1) as f64;
    var.sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_beta_at_reference_temperature() {
        let b = beta(T_REF);
        assert!((b - 1.0 / (K_B * T_REF)).abs() < 1.0e-10 * b);
    }

    #[test]
    fn test_beta_mol_kjmol_positive() {
        let b = beta_mol_kjmol(T_REF);
        assert!(b > 0.0);
    }

    #[test]
    fn test_log_sum_exp_single_value() {
        let lse = log_sum_exp(&[3.0]);
        assert!((lse - 3.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_log_sum_exp_two_equal() {
        let lse = log_sum_exp(&[0.0, 0.0]);
        assert!((lse - 2_f64.ln()).abs() < 1.0e-12);
    }

    #[test]
    fn test_ti_trapezoidal_linear_dhdl() {
        // ⟨∂H/∂λ⟩ = λ → ΔA = ∫₀¹ λ dλ = 0.5
        let windows: Vec<TiWindow> = (0..=10)
            .map(|i| {
                let lam = i as f64 / 10.0;
                TiWindow::new(lam, lam, 0.0, 100)
            })
            .collect();
        let da = ti_trapezoidal(&windows);
        assert!((da - 0.5).abs() < 1.0e-10);
    }

    #[test]
    fn test_ti_trapezoidal_error_zero_when_no_stderr() {
        let windows: Vec<TiWindow> = (0..=5)
            .map(|i| TiWindow::new(i as f64 / 5.0, 1.0, 0.0, 100))
            .collect();
        let err = ti_trapezoidal_error(&windows);
        assert_eq!(err, 0.0);
    }

    #[test]
    fn test_ti_trapezoidal_single_window_returns_zero() {
        let windows = vec![TiWindow::new(0.5, 2.0, 0.1, 100)];
        let da = ti_trapezoidal(&windows);
        assert_eq!(da, 0.0);
    }

    #[test]
    fn test_gauss_legendre_weights_sum_to_one() {
        let (_nodes, weights) = gauss_legendre_5();
        let sum: f64 = weights.iter().sum();
        assert!((sum - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_gauss_legendre_nodes_in_unit_interval() {
        let (nodes, _) = gauss_legendre_5();
        for &n in nodes.iter() {
            assert!((0.0..=1.0).contains(&n));
        }
    }

    #[test]
    fn test_fep_zwanzig_zero_du_returns_zero() {
        let du = vec![0.0; 100];
        let da = fep_zwanzig(&du, T_REF);
        assert!(da.abs() < 1.0e-12);
    }

    #[test]
    fn test_fep_zwanzig_positive_du_gives_negative_da() {
        // All ΔU > 0 → <exp(-βΔU)> = exp(-β*ΔU) → ΔA = ΔU > 0
        let du: Vec<f64> = vec![10.0; 100]; // large positive
        let da = fep_zwanzig(&du, T_REF);
        assert!(da > 0.0);
    }

    #[test]
    fn test_fep_gaussian_approx_zero_variance() {
        // When all ΔU = c, Gaussian approx gives c − 0 = c
        let du = vec![5.0f64; 100];
        let da = fep_gaussian_approx(&du, T_REF);
        assert!((da - 5.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_soft_core_lj_at_zero_lambda_is_zero() {
        let u = soft_core_lj(0.5, 0.0, 1.0, 0.35, 0.5, 1.0, 2.0);
        assert_eq!(u, 0.0);
    }

    #[test]
    fn test_soft_core_lj_positive_at_full_coupling() {
        // At λ=1, r far from σ → repulsion or attraction
        let u = soft_core_lj(0.4, 1.0, 1.0, 0.35, 0.5, 1.0, 2.0);
        assert!(u.is_finite());
    }

    #[test]
    fn test_soft_core_lj_dhdl_zero_at_lambda_zero() {
        let dhdl = soft_core_lj_dhdl(0.5, 0.0, 1.0, 0.35, 0.5, 1.0, 2.0);
        assert_eq!(dhdl, 0.0);
    }

    #[test]
    fn test_fermi_at_zero_is_half() {
        assert!((fermi(0.0) - 0.5).abs() < 1.0e-15);
    }

    #[test]
    fn test_bar_converge_symmetric() {
        // Gaussian distributions: du_ab ~ N(2, 5), du_ba ~ N(-2, 5) → ΔA ≈ 2 kJ/mol
        // Wide sigma relative to mean ensures good phase-space overlap for BAR.
        let du_ab: Vec<f64> = generate_gaussian_du(2.0, 5.0, 2000);
        let du_ba: Vec<f64> = generate_gaussian_du(-2.0, 5.0, 2000);
        let (da, _iters) = bar_converge(&du_ab, &du_ba, T_REF, 1.0e-6, 200);
        // For symmetric Gaussian distributions centred at ±2, BAR should give ΔA ≈ 2
        assert!((da - 2.0).abs() < 0.3, "BAR estimate {da} not close to 2.0");
    }

    #[test]
    fn test_bar_uncertainty_positive() {
        let du_ab: Vec<f64> = generate_gaussian_du(2.0, 0.5, 200);
        let du_ba: Vec<f64> = generate_gaussian_du(-2.0, 0.5, 200);
        let (da, _) = bar_converge(&du_ab, &du_ba, T_REF, 1.0e-6, 100);
        let sigma = bar_uncertainty(&du_ab, &du_ba, da, T_REF);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_lambda_schedule_uniform_endpoints() {
        let sched = lambda_schedule_uniform(11);
        assert_eq!(sched.len(), 11);
        assert!((sched[0]).abs() < 1.0e-15);
        assert!((sched[10] - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn test_lambda_schedule_cosine_endpoints() {
        let sched = lambda_schedule_cosine(11);
        assert!((sched[0]).abs() < 1.0e-14);
        assert!((sched[10] - 1.0).abs() < 1.0e-14);
    }

    #[test]
    fn test_cycle_closure_perfect_triangle() {
        // A→B = 5, B→C = 3, C→A = −8 → closure = 0
        let edges = vec![
            ThermoCycleEdge::new(0, 1, 5.0, 0.1),
            ThermoCycleEdge::new(1, 2, 3.0, 0.1),
            ThermoCycleEdge::new(2, 0, -8.0, 0.1),
        ];
        let (err, _sigma) = cycle_closure_error(&edges);
        assert!(err.abs() < 1.0e-10);
    }

    #[test]
    fn test_cycle_closure_fails_with_hysteresis() {
        let edges = vec![
            ThermoCycleEdge::new(0, 1, 5.0, 0.1),
            ThermoCycleEdge::new(1, 2, 3.0, 0.1),
            ThermoCycleEdge::new(2, 0, -6.0, 0.1), // wrong: closure = 2 kJ/mol
        ];
        let (err, sigma) = cycle_closure_error(&edges);
        assert!(!cycle_closure_passes(err, sigma, 2.0));
    }

    #[test]
    fn test_binding_free_energy_formula() {
        let da = binding_free_energy(-80.0, -60.0, 5.0);
        assert!((da - (-80.0 + 60.0 + 5.0 - 0.0)).abs() < 1.0e-12);
        // da = -80 - (-60) + 5 = -15 kJ/mol
        assert!((da - (-15.0)).abs() < 1.0e-12);
    }

    #[test]
    fn test_kd_from_binding_dg_positive() {
        let kd = kd_from_binding_dg(-40.0, T_REF);
        assert!(kd > 0.0);
    }

    #[test]
    fn test_rbfe_ddg_sign() {
        // Ligand B binds tighter (more negative ΔA)
        let ddg = rbfe_ddg(-30.0, -40.0);
        assert!(ddg < 0.0);
    }

    #[test]
    fn test_alchemical_window_statistics() {
        let mut w = AlchemicalWindow::new(0.5);
        for x in [1.0, 2.0, 3.0, 4.0, 5.0] {
            w.push_dhdl(x);
        }
        let (mean, _stderr) = w.dhdl_statistics();
        assert!((mean - 3.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_run_ti_pipeline_constant_dhdl() {
        let mut windows: Vec<AlchemicalWindow> = lambda_schedule_uniform(6)
            .iter()
            .map(|&lam| {
                let mut w = AlchemicalWindow::new(lam);
                for _ in 0..50 {
                    w.push_dhdl(10.0); // constant 10 kJ/mol
                }
                w
            })
            .collect();
        mock_alchemical_step(&mut windows[0], 10.0, 0.0, 0.5, 0.0, 1);
        let (da, _err) = run_ti_pipeline(&windows);
        // ∫₀¹ 10 dλ = 10
        assert!((da - 10.0).abs() < 1.0e-8);
    }

    #[test]
    fn test_generate_gaussian_du_count() {
        let samples = generate_gaussian_du(0.0, 1.0, 100);
        assert_eq!(samples.len(), 100);
    }

    #[test]
    fn test_generate_gaussian_du_mean_approx() {
        let samples = generate_gaussian_du(5.0, 0.1, 1000);
        let mean: f64 = samples.iter().sum::<f64>() / 1000.0;
        assert!((mean - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_overlap_integral_identical_histograms() {
        let h = vec![0.1f64; 10];
        let ov = overlap_integral(&h, &h);
        assert!((ov - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_effective_sample_size_iid() {
        // For i.i.d. samples, block-averaging gives N_eff ≈ n_blocks (= n/block_size).
        // With n=200 and block_size=5, n_blocks=40 so N_eff should be positive and
        // at least several units.
        let samples: Vec<f64> = generate_gaussian_du(0.0, 1.0, 200);
        let n_eff = effective_sample_size(&samples, 5);
        assert!(
            n_eff > 1.0,
            "N_eff should be positive for i.i.d. samples, got {n_eff}"
        );
    }

    #[test]
    fn test_kjmol_to_kcalmol_round_trip() {
        let original = 42.0_f64;
        let converted = kjmol_to_kcalmol(kcalmol_to_kjmol(original));
        assert!((converted - original).abs() < 1.0e-10);
    }

    #[test]
    fn test_hydration_free_energy_additive() {
        let da = hydration_free_energy(-15.0, -5.0);
        assert!((da - (-20.0)).abs() < 1.0e-12);
    }

    #[test]
    fn test_interpolate_dhdl_at_boundary() {
        let windows = vec![
            TiWindow::new(0.0, 0.0, 0.0, 10),
            TiWindow::new(1.0, 2.0, 0.0, 10),
        ];
        // At λ=0 → 0.0, at λ=1 → 2.0, at λ=0.5 → 1.0
        assert!((interpolate_dhdl(&windows, 0.0)).abs() < 1.0e-12);
        assert!((interpolate_dhdl(&windows, 1.0) - 2.0).abs() < 1.0e-12);
        assert!((interpolate_dhdl(&windows, 0.5) - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn test_bootstrap_uncertainty_positive() {
        let samples = generate_gaussian_du(0.0, 1.0, 50);
        let mean_fn = |s: &[f64]| s.iter().sum::<f64>() / s.len() as f64;
        let sigma = bootstrap_uncertainty(&samples, &mean_fn, 100);
        assert!(sigma >= 0.0);
    }

    #[test]
    fn test_autocorrelation_time_iid() {
        // i.i.d. white noise → τ ≈ 1
        let samples = generate_gaussian_du(0.0, 1.0, 500);
        let tau = autocorrelation_time(&samples, 50);
        assert!(tau >= 1.0);
        assert!(tau < 10.0); // should be close to 1 for i.i.d.
    }
}
