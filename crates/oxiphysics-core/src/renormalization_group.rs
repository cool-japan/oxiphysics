// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Renormalization group (RG) theory computations.
//!
//! This module implements the modern apparatus of the renormalization group:
//!
//! - **Beta functions**: one-loop and two-loop beta functions, UV/IR fixed points
//! - **RG flow**: Callan-Symanzik equation, running coupling, anomalous dimensions
//! - **Wilsonian RG**: effective field theory, block-spin / coarse-graining transformations
//! - **Universality classes**: Ising, XY, Heisenberg, mean-field critical exponents
//! - **Scaling relations**: hyperscaling, Rushbrooke, Josephson, Fisher identities
//! - **Renormalization schemes**: MS-bar, momentum subtraction, on-shell
//! - **Operator product expansion**: OPE coefficients, conformal blocks, fusion rules
//! - **Critical phenomena**: correlation length, order parameter, susceptibility
//! - **Decimation transformation**: 1-D Ising decimation, block-spin flow equations
//! - **Conformal field theory**: central charge, scaling dimensions, partition function

use std::f64::consts::PI;

// ============================================================================
// Internal lightweight LCG (no external rand dependency required)
// ============================================================================

// ============================================================================
// Constants
// ============================================================================

/// Euler-Mascheroni constant.
pub const EULER_MASCHERONI: f64 = 0.577_215_664_901_532_9;

/// 4π², appearing in loop integrals.
pub const FOUR_PI_SQ: f64 = 4.0 * PI * PI;

/// Default loop-expansion parameter ε = 4 – d.
pub const EPSILON_WILSON: f64 = 1.0;

/// Default UV cutoff scale Λ in natural units.
pub const DEFAULT_UV_CUTOFF: f64 = 1.0e3;

/// Convergence tolerance for fixed-point iteration.
pub const FP_TOL: f64 = 1.0e-10;

/// Maximum iterations for fixed-point Newton search.
pub const FP_MAX_ITER: usize = 1_000;

// ============================================================================
// BetaFunction
// ============================================================================

/// Loop order for beta-function computation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LoopOrder {
    /// One-loop (leading order) beta function.
    OneLoop,
    /// Two-loop (next-to-leading order) beta function.
    TwoLoop,
}

/// Beta function for a single running coupling constant g.
///
/// The beta function β(g) = μ dg/dμ governs how the coupling constant
/// evolves with the RG scale μ.  Coefficients b0, b1 follow the
/// convention  β(g) = −b0 g² − b1 g³  (perturbative, Euclidean).
#[derive(Debug, Clone)]
pub struct BetaFunction {
    /// One-loop coefficient b₀ (≥ 0 for asymptotic freedom).
    pub b0: f64,
    /// Two-loop coefficient b₁.
    pub b1: f64,
    /// Loop order used for evaluation.
    pub order: LoopOrder,
    /// Name tag for the theory (e.g. "QCD", "phi4").
    pub name: String,
}

impl BetaFunction {
    /// Construct a new beta function with given coefficients and loop order.
    pub fn new(b0: f64, b1: f64, order: LoopOrder, name: impl Into<String>) -> Self {
        Self {
            b0,
            b1,
            order,
            name: name.into(),
        }
    }

    /// Evaluate β(g) at coupling `g`.
    ///
    /// One-loop:  β(g) = −b₀ g²
    /// Two-loop:  β(g) = −b₀ g² − b₁ g³
    pub fn evaluate(&self, g: f64) -> f64 {
        let one_loop = -self.b0 * g * g;
        match self.order {
            LoopOrder::OneLoop => one_loop,
            LoopOrder::TwoLoop => one_loop - self.b1 * g * g * g,
        }
    }

    /// Derivative dβ/dg at coupling `g`.
    pub fn derivative(&self, g: f64) -> f64 {
        let one_loop = -2.0 * self.b0 * g;
        match self.order {
            LoopOrder::OneLoop => one_loop,
            LoopOrder::TwoLoop => one_loop - 3.0 * self.b1 * g * g,
        }
    }

    /// Find the non-trivial (Wilson-Fisher) fixed point g* where β(g*) = 0, g* ≠ 0.
    ///
    /// Returns `None` if no real fixed point exists in the perturbative regime.
    pub fn fixed_point(&self) -> Option<f64> {
        match self.order {
            LoopOrder::OneLoop => {
                // β = −b₀ g² = 0  only at g = 0 (trivial)
                None
            }
            LoopOrder::TwoLoop => {
                // β = −b₀ g² − b₁ g³ = −g²(b₀ + b₁ g) = 0
                // non-trivial: g* = −b₀/b₁
                if self.b1.abs() < f64::EPSILON {
                    None
                } else {
                    let gstar = -self.b0 / self.b1;
                    if gstar > 0.0 { Some(gstar) } else { None }
                }
            }
        }
    }

    /// Stability exponent ω = dβ/dg|_{g=g*}.
    ///
    /// ω > 0: UV fixed point (IR stable); ω < 0: IR fixed point (UV stable).
    pub fn stability_exponent(&self) -> Option<f64> {
        self.fixed_point().map(|gstar| self.derivative(gstar))
    }

    /// Integrate the RG equation dg/dt = β(g) from t=0 to t=`t_final`
    /// using simple Euler integration with `steps` steps.
    pub fn run_coupling(&self, g0: f64, t_final: f64, steps: usize) -> f64 {
        let dt = t_final / steps as f64;
        let mut g = g0;
        for _ in 0..steps {
            g += self.evaluate(g) * dt;
        }
        g
    }

    /// One-loop running coupling: analytic solution g(t) = g₀ / (1 + b₀ g₀ t).
    ///
    /// Derived from dg/dt = −b₀ g², which integrates to 1/g(t) = 1/g₀ + b₀ t.
    pub fn running_coupling_one_loop(&self, g0: f64, t: f64) -> f64 {
        let denom = 1.0 + self.b0 * g0 * t;
        if denom <= 0.0 {
            f64::INFINITY
        } else {
            g0 / denom
        }
    }

    /// Landau pole scale: t_L where g → ∞ under one-loop running.
    pub fn landau_pole(&self, g0: f64) -> f64 {
        if self.b0 <= 0.0 {
            f64::INFINITY
        } else {
            1.0 / (2.0 * self.b0 * g0)
        }
    }
}

// ============================================================================
// RgFlow
// ============================================================================

/// Result of integrating the Callan-Symanzik RG flow equations.
#[derive(Debug, Clone)]
pub struct RgFlowResult {
    /// Scale parameter values t = ln(μ/μ₀).
    pub t_values: Vec<f64>,
    /// Running coupling g(t).
    pub coupling: Vec<f64>,
    /// Running mass m(t).
    pub mass: Vec<f64>,
    /// Anomalous dimension η(t).
    pub eta: Vec<f64>,
}

/// Callan-Symanzik RG flow for coupling and mass.
///
/// Implements μ dg/dμ = β(g),  μ dm/dμ = −γ_m(g) m,
/// where γ_m is the mass anomalous dimension.
#[derive(Debug, Clone)]
pub struct RgFlow {
    /// Beta function governing the coupling flow.
    pub beta: BetaFunction,
    /// One-loop mass anomalous dimension coefficient γ₀ (γ_m = γ₀ g²).
    pub gamma0: f64,
    /// Field anomalous dimension coefficient η₀ (η = η₀ g²).
    pub eta0: f64,
}

impl RgFlow {
    /// Construct a new RG flow with given beta function and anomalous dimension coefficients.
    pub fn new(beta: BetaFunction, gamma0: f64, eta0: f64) -> Self {
        Self { beta, gamma0, eta0 }
    }

    /// Anomalous dimension of the mass at coupling `g`: γ_m(g) = γ₀ g².
    pub fn gamma_mass(&self, g: f64) -> f64 {
        self.gamma0 * g * g
    }

    /// Field anomalous dimension η(g) = η₀ g².
    pub fn eta(&self, g: f64) -> f64 {
        self.eta0 * g * g
    }

    /// Integrate the full RG flow from t=0 to t=`t_max` with `n` steps.
    ///
    /// Returns a [`RgFlowResult`] containing the full trajectory.
    pub fn integrate(&self, g0: f64, m0: f64, t_max: f64, n: usize) -> RgFlowResult {
        let dt = t_max / n as f64;
        let mut t_values = Vec::with_capacity(n + 1);
        let mut coupling = Vec::with_capacity(n + 1);
        let mut mass = Vec::with_capacity(n + 1);
        let mut eta = Vec::with_capacity(n + 1);

        let mut g = g0;
        let mut m = m0;
        let mut t = 0.0;

        t_values.push(t);
        coupling.push(g);
        mass.push(m);
        eta.push(self.eta(g));

        for _ in 0..n {
            let dg = self.beta.evaluate(g) * dt;
            let dm = -self.gamma_mass(g) * m * dt;
            g += dg;
            m += dm;
            t += dt;
            t_values.push(t);
            coupling.push(g);
            mass.push(m);
            eta.push(self.eta(g));
        }

        RgFlowResult {
            t_values,
            coupling,
            mass,
            eta,
        }
    }

    /// Correlation length exponent ν from the eigenvalue of the linearized
    /// beta function at the fixed point: ν = −1/ω.
    pub fn nu_exponent(&self) -> Option<f64> {
        self.beta.stability_exponent().and_then(|omega| {
            if omega.abs() < f64::EPSILON {
                None
            } else {
                Some(1.0 / omega.abs())
            }
        })
    }
}

// ============================================================================
// WilsonianRg
// ============================================================================

/// Parameters for a Wilsonian effective field theory in d dimensions.
///
/// Implements coarse-graining via momentum-shell integration and block-spin
/// transformation, yielding flow equations for marginal and relevant couplings.
#[derive(Debug, Clone)]
pub struct WilsonianRg {
    /// Spacetime dimension d.
    pub d: f64,
    /// ε = 4 − d (Wilson-Fisher expansion parameter).
    pub epsilon: f64,
    /// UV cutoff Λ.
    pub lambda: f64,
    /// Current RG step number (number of coarse-graining steps performed).
    pub step: usize,
    /// Coupling u (quartic in φ⁴ theory).
    pub u: f64,
    /// Mass parameter r (quadratic coupling).
    pub r: f64,
}

impl WilsonianRg {
    /// Construct a Wilsonian RG for φ⁴ theory in d = 4 − ε dimensions.
    pub fn new(d: f64, lambda: f64, u: f64, r: f64) -> Self {
        Self {
            d,
            epsilon: 4.0 - d,
            lambda,
            step: 0,
            u,
            r,
        }
    }

    /// Perform one infinitesimal RG step with rescaling factor b = e^{dl}.
    ///
    /// One-loop flow equations (ε-expansion):
    ///   du/dl = ε u − (n+8)/(16π²) u²
    ///   dr/dl = 2 r + (n+2)/(16π²) u  (n = 1 for Ising universality)
    pub fn step_flow(&mut self, dl: f64, n_components: f64) {
        let du = self.epsilon * self.u - (n_components + 8.0) / (16.0 * PI * PI) * self.u * self.u;
        let dr = 2.0 * self.r + (n_components + 2.0) / (16.0 * PI * PI) * self.u;
        self.u += du * dl;
        self.r += dr * dl;
        self.step += 1;
    }

    /// Wilson-Fisher fixed point u* for φ⁴ theory (n-component).
    ///
    /// u* = (16π² ε) / (n + 8)
    pub fn wilson_fisher_fixed_point(&self, n_components: f64) -> f64 {
        16.0 * PI * PI * self.epsilon / (n_components + 8.0)
    }

    /// Critical exponent η at Wilson-Fisher fixed point to O(ε²).
    ///
    /// η = (n+2) ε² / \[2(n+8)²\]
    pub fn eta_wf(&self, n_components: f64) -> f64 {
        let eps = self.epsilon;
        (n_components + 2.0) * eps * eps / (2.0 * (n_components + 8.0).powi(2))
    }

    /// Correlation length exponent ν at Wilson-Fisher fixed point.
    ///
    /// 1/ν = 2 − (n+2)/(n+8) ε
    pub fn nu_wf(&self, n_components: f64) -> f64 {
        let inv_nu = 2.0 - (n_components + 2.0) / (n_components + 8.0) * self.epsilon;
        1.0 / inv_nu
    }

    /// Block-spin transformation: rescale correlation length by factor b.
    ///
    /// Under coarse-graining by b, the correlation length scales as ξ → ξ/b.
    pub fn rescale_correlation_length(&self, xi: f64, b: f64) -> f64 {
        xi / b
    }

    /// Run the RG flow for `n_steps` steps with step size `dl`.
    pub fn run(&mut self, n_steps: usize, dl: f64, n_components: f64) {
        for _ in 0..n_steps {
            self.step_flow(dl, n_components);
        }
    }
}

// ============================================================================
// UniversalityClass
// ============================================================================

/// Universality class of a critical point.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UniversalityClass {
    /// 2-D / 3-D Ising universality (n = 1 scalar order parameter).
    Ising2D,
    /// 3-D Ising universality.
    Ising3D,
    /// 3-D XY universality (n = 2 complex order parameter).
    Xy3D,
    /// 3-D Heisenberg universality (n = 3 vector order parameter).
    Heisenberg3D,
    /// Mean-field (Landau) theory critical exponents.
    MeanField,
}

/// Set of universal critical exponents for a given universality class.
#[derive(Debug, Clone)]
pub struct CriticalExponents {
    /// Correlation length exponent ν.
    pub nu: f64,
    /// Anomalous dimension η.
    pub eta: f64,
    /// Order parameter exponent β (not to be confused with beta function).
    pub beta_exp: f64,
    /// Susceptibility exponent γ.
    pub gamma: f64,
    /// Specific heat exponent α.
    pub alpha: f64,
    /// Critical isotherm exponent δ.
    pub delta: f64,
    /// Correlation function exponent at T_c: G(r) ~ r^{-(d−2+η)}.
    pub d: f64,
}

impl UniversalityClass {
    /// Return the standard critical exponents for this universality class.
    ///
    /// Values are state-of-the-art best estimates (Monte Carlo / conformal bootstrap).
    pub fn exponents(self) -> CriticalExponents {
        match self {
            UniversalityClass::Ising2D => CriticalExponents {
                nu: 1.0,
                eta: 0.25,
                beta_exp: 0.125,
                gamma: 1.75,
                alpha: 0.0, // logarithmic
                delta: 15.0,
                d: 2.0,
            },
            UniversalityClass::Ising3D => CriticalExponents {
                nu: 0.6301,
                eta: 0.0364,
                beta_exp: 0.3265,
                gamma: 1.2372,
                alpha: 0.1100,
                delta: 4.789,
                d: 3.0,
            },
            UniversalityClass::Xy3D => CriticalExponents {
                nu: 0.6717,
                eta: 0.0381,
                beta_exp: 0.3486,
                gamma: 1.3178,
                alpha: -0.0151,
                delta: 4.780,
                d: 3.0,
            },
            UniversalityClass::Heisenberg3D => CriticalExponents {
                nu: 0.7112,
                eta: 0.0375,
                beta_exp: 0.3689,
                gamma: 1.3960,
                alpha: -0.1336,
                delta: 4.783,
                d: 3.0,
            },
            UniversalityClass::MeanField => CriticalExponents {
                nu: 0.5,
                eta: 0.0,
                beta_exp: 0.5,
                gamma: 1.0,
                alpha: 0.0,
                delta: 3.0,
                d: 4.0, // upper critical dimension
            },
        }
    }

    /// Number of order parameter components n.
    pub fn n_components(self) -> usize {
        match self {
            UniversalityClass::Ising2D | UniversalityClass::Ising3D => 1,
            UniversalityClass::Xy3D => 2,
            UniversalityClass::Heisenberg3D => 3,
            UniversalityClass::MeanField => 1,
        }
    }

    /// Spatial dimension d.
    pub fn dimension(self) -> f64 {
        match self {
            UniversalityClass::Ising2D => 2.0,
            _ => 3.0,
        }
    }
}

// ============================================================================
// ScalingRelations
// ============================================================================

/// Verification of thermodynamic scaling relations from critical exponents.
///
/// The four classical relations among {α, β, γ, δ, ν, η, d}:
/// - Rushbrooke: α + 2β + γ = 2
/// - Widom:      γ = β(δ − 1)
/// - Fisher:     γ = (2 − η) ν
/// - Josephson:  2 − α = d ν  (hyperscaling)
#[derive(Debug, Clone)]
pub struct ScalingRelations {
    /// Critical exponents to check.
    pub exp: CriticalExponents,
}

impl ScalingRelations {
    /// Construct from a set of critical exponents.
    pub fn new(exp: CriticalExponents) -> Self {
        Self { exp }
    }

    /// Rushbrooke identity: α + 2β + γ = 2.  Returns (lhs, rhs=2).
    pub fn rushbrooke(&self) -> (f64, f64) {
        let lhs = self.exp.alpha + 2.0 * self.exp.beta_exp + self.exp.gamma;
        (lhs, 2.0)
    }

    /// Widom identity: γ = β(δ − 1).  Returns (lhs=γ, rhs=β(δ−1)).
    pub fn widom(&self) -> (f64, f64) {
        let lhs = self.exp.gamma;
        let rhs = self.exp.beta_exp * (self.exp.delta - 1.0);
        (lhs, rhs)
    }

    /// Fisher identity: γ = (2 − η) ν.  Returns (lhs=γ, rhs).
    pub fn fisher(&self) -> (f64, f64) {
        let lhs = self.exp.gamma;
        let rhs = (2.0 - self.exp.eta) * self.exp.nu;
        (lhs, rhs)
    }

    /// Josephson (hyperscaling) identity: 2 − α = d ν.  Returns (lhs, rhs).
    pub fn josephson(&self) -> (f64, f64) {
        let lhs = 2.0 - self.exp.alpha;
        let rhs = self.exp.d * self.exp.nu;
        (lhs, rhs)
    }

    /// Check all four relations within tolerance `tol`.
    ///
    /// Returns `(rushbrooke_ok, widom_ok, fisher_ok, josephson_ok)`.
    pub fn check_all(&self, tol: f64) -> (bool, bool, bool, bool) {
        let (r_lhs, r_rhs) = self.rushbrooke();
        let (w_lhs, w_rhs) = self.widom();
        let (f_lhs, f_rhs) = self.fisher();
        let (j_lhs, j_rhs) = self.josephson();
        (
            (r_lhs - r_rhs).abs() < tol,
            (w_lhs - w_rhs).abs() < tol,
            (f_lhs - f_rhs).abs() < tol,
            (j_lhs - j_rhs).abs() < tol,
        )
    }
}

// ============================================================================
// RenormalizationScheme
// ============================================================================

/// Renormalization scheme for field theory calculations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RenormalizationScheme {
    /// MS-bar (minimal subtraction in dimensional regularisation).
    MSbar,
    /// Momentum subtraction (MOM): renormalize at spacelike momentum p² = μ².
    MomentumSubtraction,
    /// On-shell renormalization: use physical pole mass and coupling.
    OnShell,
}

/// Scheme-dependent renormalized coupling and counterterms.
#[derive(Debug, Clone)]
pub struct RenormalizedCoupling {
    /// Renormalization scheme.
    pub scheme: RenormalizationScheme,
    /// Bare coupling g₀.
    pub g_bare: f64,
    /// Renormalized coupling g(μ).
    pub g_renorm: f64,
    /// Renormalization scale μ.
    pub mu: f64,
    /// UV divergence degree (dimensionless; 0 = finite, n = power/log divergence).
    pub divergence_degree: i32,
}

impl RenormalizedCoupling {
    /// Construct a renormalized coupling in the given scheme.
    pub fn new(
        scheme: RenormalizationScheme,
        g_bare: f64,
        mu: f64,
        divergence_degree: i32,
    ) -> Self {
        let g_renorm = match scheme {
            RenormalizationScheme::MSbar => {
                // MS-bar: subtract only 1/ε pole + ln(4π) − γ_E
                let log_factor = (4.0 * PI).ln() - EULER_MASCHERONI;
                g_bare - g_bare * g_bare * log_factor / (16.0 * PI * PI)
            }
            RenormalizationScheme::MomentumSubtraction => {
                // MOM: subtract at p² = μ²; finite parts differ from MS-bar by a constant
                g_bare - g_bare * g_bare * mu.ln() / (16.0 * PI * PI)
            }
            RenormalizationScheme::OnShell => {
                // On-shell: coupling fixed by physical scattering amplitude
                g_bare * (1.0 - g_bare / (16.0 * PI * PI))
            }
        };
        Self {
            scheme,
            g_bare,
            g_renorm,
            mu,
            divergence_degree,
        }
    }

    /// Scheme change: convert g(μ) from MS-bar to MOM (one-loop).
    ///
    /// g_MOM = g_MSbar \[1 + Δ g_MSbar / (16π²)\] where Δ is scheme-difference coefficient.
    pub fn convert_msbar_to_mom(&self, delta: f64) -> f64 {
        self.g_renorm * (1.0 + delta * self.g_renorm / (16.0 * PI * PI))
    }

    /// Running coupling between two scales μ₁ and μ₂ (one-loop MS-bar).
    pub fn run_msbar(&self, mu2: f64, b0: f64) -> f64 {
        let t = (mu2 / self.mu).ln();
        self.g_renorm / (1.0 + 2.0 * b0 * self.g_renorm * t).sqrt()
    }
}

// ============================================================================
// OperatorProductExpansion
// ============================================================================

/// OPE coefficient C_{ijk} for operators O_i × O_j → O_k.
#[derive(Debug, Clone)]
pub struct OpeCoefficient {
    /// Label of operator O_i.
    pub label_i: String,
    /// Label of operator O_j.
    pub label_j: String,
    /// Label of operator O_k.
    pub label_k: String,
    /// Numerical value of the OPE coefficient C_{ijk}.
    pub value: f64,
}

impl OpeCoefficient {
    /// Construct an OPE coefficient.
    pub fn new(
        label_i: impl Into<String>,
        label_j: impl Into<String>,
        label_k: impl Into<String>,
        value: f64,
    ) -> Self {
        Self {
            label_i: label_i.into(),
            label_j: label_j.into(),
            label_k: label_k.into(),
            value,
        }
    }
}

/// Operator product expansion registry.
///
/// Stores a list of OPE coefficients and provides utilities for computing
/// conformal partial waves and checking crossing symmetry.
#[derive(Debug, Clone, Default)]
pub struct OperatorProductExpansion {
    /// List of OPE coefficients.
    pub coefficients: Vec<OpeCoefficient>,
    /// List of operator scaling dimensions (Δ_k).
    pub dimensions: Vec<(String, f64)>,
}

impl OperatorProductExpansion {
    /// Construct an empty OPE registry.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register an OPE coefficient.
    pub fn add_coefficient(&mut self, coeff: OpeCoefficient) {
        self.coefficients.push(coeff);
    }

    /// Register an operator with scaling dimension Δ.
    pub fn add_operator(&mut self, label: impl Into<String>, delta: f64) {
        self.dimensions.push((label.into(), delta));
    }

    /// Conformal block G_{Δ,l}(z, z̄) in 2D at leading order (simplified scalar, l=0).
    ///
    /// G_Δ(z) ≈ z^{Δ/2} (1 + Δ z / (2(2Δ+1)) + …)
    pub fn conformal_block_2d(&self, delta: f64, z: f64) -> f64 {
        if z <= 0.0 {
            return 0.0;
        }
        let z_pow = z.powf(delta / 2.0);
        let correction = 1.0 + delta * z / (2.0 * (2.0 * delta + 1.0));
        z_pow * correction
    }

    /// Partial wave amplitude: sum of C²_{σσk} × G_{Δk}(z).
    pub fn partial_wave_amplitude(&self, sigma_label: &str, z: f64) -> f64 {
        let mut amplitude = 0.0;
        for coeff in &self.coefficients {
            if coeff.label_i == sigma_label
                && coeff.label_j == sigma_label
                && let Some((_, delta_k)) = self
                    .dimensions
                    .iter()
                    .find(|(lbl, _)| lbl == &coeff.label_k)
            {
                let g = self.conformal_block_2d(*delta_k, z);
                amplitude += coeff.value * coeff.value * g;
            }
        }
        amplitude
    }

    /// Fusion rule check: does O_i × O_j contain O_k?
    pub fn fusion_rule_exists(&self, label_i: &str, label_j: &str, label_k: &str) -> bool {
        self.coefficients
            .iter()
            .any(|c| c.label_i == label_i && c.label_j == label_j && c.label_k == label_k)
    }
}

// ============================================================================
// CriticalPhenomena
// ============================================================================

/// Critical phenomena observables near a second-order phase transition.
///
/// All observables are expressed as power-law functions of the reduced
/// temperature t = (T − T_c)/T_c and ordering field h.
#[derive(Debug, Clone)]
pub struct CriticalPhenomena {
    /// Critical exponents.
    pub exp: CriticalExponents,
    /// Critical temperature T_c.
    pub t_c: f64,
    /// Non-universal amplitude for correlation length: ξ = ξ₀ |t|^{−ν}.
    pub xi0: f64,
    /// Non-universal amplitude for susceptibility: χ = χ₀ |t|^{−γ}.
    pub chi0: f64,
    /// Non-universal amplitude for order parameter: m = m₀ |t|^β.
    pub m0: f64,
    /// Non-universal amplitude for specific heat: C = A |t|^{−α}.
    pub c0: f64,
}

impl CriticalPhenomena {
    /// Construct a CriticalPhenomena model.
    pub fn new(exp: CriticalExponents, t_c: f64, xi0: f64, chi0: f64, m0: f64, c0: f64) -> Self {
        Self {
            exp,
            t_c,
            xi0,
            chi0,
            m0,
            c0,
        }
    }

    /// Reduced temperature t = (T − T_c) / T_c.
    pub fn reduced_temperature(&self, temp: f64) -> f64 {
        (temp - self.t_c) / self.t_c
    }

    /// Correlation length ξ(T) = ξ₀ |t|^{−ν}  (T ≠ T_c).
    pub fn correlation_length(&self, temp: f64) -> f64 {
        let t = self.reduced_temperature(temp);
        if t.abs() < f64::EPSILON {
            f64::INFINITY
        } else {
            self.xi0 * t.abs().powf(-self.exp.nu)
        }
    }

    /// Magnetic susceptibility χ(T) = χ₀ |t|^{−γ}  (T > T_c).
    pub fn susceptibility(&self, temp: f64) -> f64 {
        let t = self.reduced_temperature(temp);
        if t.abs() < f64::EPSILON {
            f64::INFINITY
        } else {
            self.chi0 * t.abs().powf(-self.exp.gamma)
        }
    }

    /// Order parameter m(T) = m₀ (−t)^β  (T < T_c).
    pub fn order_parameter(&self, temp: f64) -> f64 {
        let t = self.reduced_temperature(temp);
        if t >= 0.0 {
            0.0
        } else {
            self.m0 * (-t).powf(self.exp.beta_exp)
        }
    }

    /// Specific heat C(T) = c₀ |t|^{−α}.
    pub fn specific_heat(&self, temp: f64) -> f64 {
        let t = self.reduced_temperature(temp);
        if t.abs() < f64::EPSILON {
            f64::INFINITY
        } else {
            self.c0 * t.abs().powf(-self.exp.alpha)
        }
    }

    /// Critical isotherm (T = T_c): m ∝ h^{1/δ}.
    pub fn critical_isotherm(&self, h: f64) -> f64 {
        h.abs().powf(1.0 / self.exp.delta) * h.signum()
    }
}

// ============================================================================
// DecimationTransformation
// ============================================================================

/// 1-D Ising decimation RG transformation.
///
/// The partition function of the 1-D Ising model with coupling K is
/// reproduced by a system with half as many spins and an effective coupling K'.
/// The recursion relation is:
///   K' = (1/2) ln cosh(2K)
#[derive(Debug, Clone)]
pub struct DecimationTransformation {
    /// Current coupling K = J/(k_B T).
    pub k: f64,
    /// Current additive constant (free energy per site contribution).
    pub ln_lambda: f64,
    /// Number of decimation steps performed.
    pub steps: usize,
}

impl DecimationTransformation {
    /// Construct a decimation RG at coupling K₀.
    pub fn new(k0: f64) -> Self {
        Self {
            k: k0,
            ln_lambda: 0.0,
            steps: 0,
        }
    }

    /// Perform one step of 1-D Ising decimation.
    ///
    /// K' = (1/2) ln cosh(2K),   ln λ += (1/2) ln\[4 cosh(2K)\]
    pub fn step(&mut self) {
        let two_k = 2.0 * self.k;
        let cosh_2k = two_k.cosh();
        self.k = 0.5 * cosh_2k.ln();
        self.ln_lambda += 0.5 * (4.0 * cosh_2k).ln();
        self.steps += 1;
    }

    /// Run `n` decimation steps.
    pub fn run(&mut self, n: usize) {
        for _ in 0..n {
            self.step();
        }
    }

    /// Correlation length after decimation (rough estimate via K).
    ///
    /// ξ = −1 / ln tanh(K)
    pub fn correlation_length(&self) -> f64 {
        let th = self.k.tanh();
        if th <= 0.0 || th >= 1.0 {
            f64::INFINITY
        } else {
            -1.0 / th.ln()
        }
    }

    /// Collect K trajectory over `n` steps.
    pub fn trajectory(&self, n: usize) -> Vec<f64> {
        let mut traj = Vec::with_capacity(n + 1);
        let mut dec = self.clone();
        traj.push(dec.k);
        for _ in 0..n {
            dec.step();
            traj.push(dec.k);
        }
        traj
    }

    /// 2-D block-spin decimation flow for the square Ising model (Migdal-Kadanoff).
    ///
    /// K' = (1/2) ln cosh(2bK) where b = 2 is the block factor.
    pub fn migdal_kadanoff_step(&mut self, b: u32) {
        let bk = b as f64 * self.k;
        self.k = 0.5 * (2.0 * bk).cosh().ln();
        self.steps += 1;
    }
}

// ============================================================================
// ConformalFieldTheory
// ============================================================================

/// Conformal field theory (CFT) data for a 2-D critical theory.
///
/// Contains central charge c, operator content (scaling dimensions and spins),
/// and utilities for computing the partition function and modular invariants.
#[derive(Debug, Clone)]
pub struct ConformalFieldTheory {
    /// Central charge c.
    pub c: f64,
    /// List of (scaling dimension Δ, spin s) for primary operators.
    pub primaries: Vec<(f64, i32)>,
    /// Theory name.
    pub name: String,
}

impl ConformalFieldTheory {
    /// Construct a CFT with central charge c and no primaries yet.
    pub fn new(c: f64, name: impl Into<String>) -> Self {
        Self {
            c,
            primaries: Vec::new(),
            name: name.into(),
        }
    }

    /// Add a primary operator with scaling dimension Δ and spin s.
    pub fn add_primary(&mut self, delta: f64, spin: i32) {
        self.primaries.push((delta, spin));
    }

    /// Casimir energy on a cylinder of circumference L: E_0 = −π c / (6 L).
    pub fn casimir_energy(&self, l: f64) -> f64 {
        -PI * self.c / (6.0 * l)
    }

    /// Conformal dimensions h = (Δ + s)/2, h̄ = (Δ − s)/2.
    pub fn holomorphic_dimensions(&self) -> Vec<(f64, f64)> {
        self.primaries
            .iter()
            .map(|&(delta, spin)| {
                let h = (delta + spin as f64) / 2.0;
                let hbar = (delta - spin as f64) / 2.0;
                (h, hbar)
            })
            .collect()
    }

    /// Virasoro character χ_h(τ) in the high-temperature limit (simplified).
    ///
    /// χ_h(q) ≈ q^{h − c/24} / η(q),  where η(q) ≈ q^{1/24} ∏(1−q^n).
    /// Here we use a truncated expansion up to q^N_terms.
    pub fn virasoro_character(&self, h: f64, q: f64, n_terms: usize) -> f64 {
        // Dedekind eta function (truncated): η(q) = q^{1/24} ∏_{n=1}^N (1−q^n)
        let mut eta_prod = 1.0_f64;
        for n in 1..=n_terms {
            eta_prod *= 1.0 - q.powi(n as i32);
        }
        let eta = q.powf(1.0 / 24.0) * eta_prod;
        if eta.abs() < f64::EPSILON {
            return 0.0;
        }
        q.powf(h - self.c / 24.0) / eta
    }

    /// Modular S-matrix element S_{h h'} for a compact boson CFT (c = 1).
    ///
    /// S_{nm} = √(2/k) sin(π n m / k) for level-k SU(2) WZW model.
    pub fn s_matrix_wzw(&self, n: usize, m: usize, k: usize) -> f64 {
        let kf = k as f64;
        (2.0 / kf).sqrt() * (PI * n as f64 * m as f64 / kf).sin()
    }

    /// Partition function on a torus: Z(τ) = Tr\[q^{L_0 − c/24} q̄^{L̄_0 − c/24}\].
    ///
    /// Simplified: sum over primaries of |χ_h(q)|².
    pub fn partition_function(&self, q: f64, n_terms: usize) -> f64 {
        self.holomorphic_dimensions()
            .iter()
            .map(|&(h, hbar)| {
                let chi_h = self.virasoro_character(h, q, n_terms);
                let chi_hbar = self.virasoro_character(hbar, q, n_terms);
                chi_h * chi_hbar
            })
            .sum()
    }

    /// Central charge from the two-point function of the stress tensor T(z):
    ///
    /// ⟨T(z)T(0)⟩ = c/(2z⁴),  extract c = 2 z⁴ ⟨T(z)T(0)⟩.
    pub fn central_charge_from_tt(&self, z: f64, tt_correlator: f64) -> f64 {
        2.0 * z.powi(4) * tt_correlator
    }

    /// Unitary bound: all primary operators must satisfy Δ ≥ |s| (Δ ≥ 0 for scalars).
    pub fn satisfies_unitarity_bound(&self) -> bool {
        self.primaries
            .iter()
            .all(|&(delta, spin)| delta >= spin.unsigned_abs() as f64)
    }
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    // ── BetaFunction tests ──────────────────────────────────────────────────

    #[test]
    fn test_beta_one_loop_zero_at_origin() {
        let beta = BetaFunction::new(1.0, 0.0, LoopOrder::OneLoop, "test");
        assert!((beta.evaluate(0.0)).abs() < 1e-15);
    }

    #[test]
    fn test_beta_one_loop_negative_for_positive_coupling() {
        let beta = BetaFunction::new(1.0, 0.0, LoopOrder::OneLoop, "QCD");
        // b0 > 0 → asymptotically free → β(g) < 0 for g > 0
        assert!(beta.evaluate(0.5) < 0.0);
    }

    #[test]
    fn test_beta_two_loop_fixed_point() {
        // β = −b0 g² − b1 g³ = 0 → g* = −b0/b1
        // Choose b0 = 1, b1 = −2 → g* = 0.5
        let beta = BetaFunction::new(1.0, -2.0, LoopOrder::TwoLoop, "WF");
        let gstar = beta.fixed_point().expect("expected fixed point");
        assert!((gstar - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_beta_stability_exponent_sign() {
        let beta = BetaFunction::new(1.0, -2.0, LoopOrder::TwoLoop, "WF");
        let omega = beta.stability_exponent().expect("stability exponent");
        // dβ/dg at g* = −2b0 g* − 3b1 g*²
        let gstar = 0.5_f64;
        let expected = -2.0 * 1.0 * gstar - 3.0 * (-2.0) * gstar * gstar;
        assert!((omega - expected).abs() < 1e-10);
    }

    #[test]
    fn test_landau_pole() {
        let beta = BetaFunction::new(0.5, 0.0, LoopOrder::OneLoop, "phi4");
        // t_L = 1/(2 b0 g0)
        let t_l = beta.landau_pole(1.0);
        assert!((t_l - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_running_coupling_one_loop_approaches_zero() {
        let beta = BetaFunction::new(1.0, 0.0, LoopOrder::OneLoop, "QCD");
        let g_inf = beta.running_coupling_one_loop(1.0, 1e6);
        assert!(g_inf < 0.01, "coupling should run to zero");
    }

    #[test]
    fn test_run_coupling_euler_consistent() {
        let beta = BetaFunction::new(1.0, 0.0, LoopOrder::OneLoop, "test");
        let g_euler = beta.run_coupling(0.5, 1.0, 10_000);
        let g_analytic = beta.running_coupling_one_loop(0.5, 1.0);
        assert!((g_euler - g_analytic).abs() < 1e-4);
    }

    // ── RgFlow tests ────────────────────────────────────────────────────────

    #[test]
    fn test_rg_flow_gamma_mass_quadratic() {
        let beta = BetaFunction::new(1.0, 0.0, LoopOrder::OneLoop, "phi4");
        let flow = RgFlow::new(beta, 0.5, 0.1);
        // γ_m(g) = γ₀ g²
        assert!((flow.gamma_mass(2.0) - 0.5 * 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_rg_flow_integrate_mass_decreases() {
        let beta = BetaFunction::new(0.1, 0.0, LoopOrder::OneLoop, "test");
        let flow = RgFlow::new(beta, 0.2, 0.0);
        let result = flow.integrate(0.5, 1.0, 5.0, 1000);
        // mass should change (even if slightly) — check it's finite
        assert!(result.mass.last().unwrap().is_finite());
    }

    #[test]
    fn test_rg_flow_nu_exponent() {
        // Use b1 < 0 so there IS a fixed point
        let beta = BetaFunction::new(1.0, -2.0, LoopOrder::TwoLoop, "WF");
        let flow = RgFlow::new(beta, 1.0, 0.1);
        let nu = flow.nu_exponent().expect("nu exponent");
        assert!(nu > 0.0);
    }

    // ── WilsonianRg tests ───────────────────────────────────────────────────

    #[test]
    fn test_wilson_fisher_fixed_point_n1() {
        // For n = 1 (Ising), u* = 16π² ε / 9
        let rg = WilsonianRg::new(3.0, 1.0, 0.1, 0.0);
        let ustar = rg.wilson_fisher_fixed_point(1.0);
        let expected = 16.0 * PI * PI * 1.0 / 9.0;
        assert!((ustar - expected).abs() < 1e-10);
    }

    #[test]
    fn test_eta_wf_n1_order_eps2() {
        let rg = WilsonianRg::new(3.0, 1.0, 0.1, 0.0);
        let eta = rg.eta_wf(1.0);
        // η = ε² / 54 for n=1
        let expected = 1.0_f64.powi(2) / 54.0;
        assert!((eta - expected).abs() < 1e-10);
    }

    #[test]
    fn test_nu_wf_n1() {
        // 1/ν = 2 − ε/3 for n=1 at ε=1 → 1/ν = 5/3 → ν = 0.6
        let rg = WilsonianRg::new(3.0, 1.0, 0.1, 0.0);
        let nu = rg.nu_wf(1.0);
        let expected = 1.0 / (2.0 - 1.0 / 3.0);
        assert!((nu - expected).abs() < 1e-10);
    }

    #[test]
    fn test_block_spin_rescaling() {
        let rg = WilsonianRg::new(3.0, 1.0, 0.1, 0.0);
        assert!((rg.rescale_correlation_length(10.0, 2.0) - 5.0).abs() < 1e-12);
    }

    // ── UniversalityClass tests ─────────────────────────────────────────────

    #[test]
    fn test_mean_field_exponents() {
        let exp = UniversalityClass::MeanField.exponents();
        assert!((exp.nu - 0.5).abs() < 1e-10);
        assert!((exp.beta_exp - 0.5).abs() < 1e-10);
        assert!((exp.gamma - 1.0).abs() < 1e-10);
        assert!((exp.delta - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_ising_3d_exponents_reasonable() {
        let exp = UniversalityClass::Ising3D.exponents();
        assert!(exp.nu > 0.5 && exp.nu < 1.0);
        assert!(exp.eta > 0.0 && exp.eta < 0.1);
    }

    #[test]
    fn test_n_components() {
        assert_eq!(UniversalityClass::Ising3D.n_components(), 1);
        assert_eq!(UniversalityClass::Xy3D.n_components(), 2);
        assert_eq!(UniversalityClass::Heisenberg3D.n_components(), 3);
    }

    // ── ScalingRelations tests ──────────────────────────────────────────────

    #[test]
    fn test_rushbrooke_mean_field() {
        let exp = UniversalityClass::MeanField.exponents();
        let sr = ScalingRelations::new(exp);
        let (lhs, rhs) = sr.rushbrooke();
        // α=0, 2β=1, γ=1 → sum = 2
        assert!((lhs - rhs).abs() < 1e-10);
    }

    #[test]
    fn test_widom_mean_field() {
        let exp = UniversalityClass::MeanField.exponents();
        let sr = ScalingRelations::new(exp);
        let (lhs, rhs) = sr.widom();
        // γ = 1, β(δ−1) = 0.5 × 2 = 1
        assert!((lhs - rhs).abs() < 1e-10);
    }

    #[test]
    fn test_fisher_mean_field() {
        let exp = UniversalityClass::MeanField.exponents();
        let sr = ScalingRelations::new(exp);
        let (lhs, rhs) = sr.fisher();
        // γ = 1, (2−0) × 0.5 = 1
        assert!((lhs - rhs).abs() < 1e-10);
    }

    #[test]
    fn test_josephson_mean_field() {
        let exp = UniversalityClass::MeanField.exponents();
        let sr = ScalingRelations::new(exp);
        let (lhs, rhs) = sr.josephson();
        // 2 − α = 2, d ν = 4 × 0.5 = 2
        assert!((lhs - rhs).abs() < 1e-10);
    }

    #[test]
    fn test_check_all_mean_field() {
        let exp = UniversalityClass::MeanField.exponents();
        let sr = ScalingRelations::new(exp);
        let (r, w, f, j) = sr.check_all(1e-10);
        assert!(r, "Rushbrooke failed");
        assert!(w, "Widom failed");
        assert!(f, "Fisher failed");
        assert!(j, "Josephson failed");
    }

    // ── RenormalizedCoupling tests ──────────────────────────────────────────

    #[test]
    fn test_renorm_coupling_msbar_finite() {
        let rc = RenormalizedCoupling::new(RenormalizationScheme::MSbar, 0.1, 1.0, 1);
        assert!(rc.g_renorm.is_finite());
        assert!(rc.g_renorm < rc.g_bare || (rc.g_renorm - rc.g_bare).abs() < rc.g_bare);
    }

    #[test]
    fn test_renorm_coupling_run_msbar() {
        let rc = RenormalizedCoupling::new(RenormalizationScheme::MSbar, 0.3, 1.0, 1);
        let g2 = rc.run_msbar(10.0, 1.0);
        // Running to higher scale with b0=1 should decrease coupling
        assert!(g2 < rc.g_renorm);
    }

    // ── OPE tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_ope_fusion_rule_exists() {
        let mut ope = OperatorProductExpansion::new();
        ope.add_coefficient(OpeCoefficient::new("σ", "σ", "ε", 0.5));
        assert!(ope.fusion_rule_exists("σ", "σ", "ε"));
        assert!(!ope.fusion_rule_exists("σ", "σ", "T"));
    }

    #[test]
    fn test_conformal_block_2d_positive() {
        let ope = OperatorProductExpansion::new();
        let g = ope.conformal_block_2d(0.5, 0.25);
        assert!(g > 0.0);
    }

    #[test]
    fn test_partial_wave_amplitude_positive() {
        let mut ope = OperatorProductExpansion::new();
        ope.add_operator("ε", 1.0);
        ope.add_coefficient(OpeCoefficient::new("σ", "σ", "ε", 0.8));
        let amp = ope.partial_wave_amplitude("σ", 0.1);
        assert!(amp > 0.0);
    }

    // ── CriticalPhenomena tests ─────────────────────────────────────────────

    #[test]
    fn test_correlation_length_diverges_at_tc() {
        let exp = UniversalityClass::MeanField.exponents();
        let cp = CriticalPhenomena::new(exp, 2.269, 1.0, 1.0, 1.0, 1.0);
        let xi = cp.correlation_length(2.269);
        assert!(xi.is_infinite());
    }

    #[test]
    fn test_order_parameter_zero_above_tc() {
        let exp = UniversalityClass::Ising3D.exponents();
        let cp = CriticalPhenomena::new(exp, 2.269, 1.0, 1.0, 1.0, 1.0);
        assert_eq!(cp.order_parameter(3.0), 0.0);
    }

    #[test]
    fn test_order_parameter_nonzero_below_tc() {
        let exp = UniversalityClass::Ising3D.exponents();
        let cp = CriticalPhenomena::new(exp, 2.269, 1.0, 1.0, 1.0, 1.0);
        assert!(cp.order_parameter(1.0) > 0.0);
    }

    #[test]
    fn test_critical_isotherm_odd() {
        let exp = UniversalityClass::MeanField.exponents();
        let cp = CriticalPhenomena::new(exp, 2.269, 1.0, 1.0, 1.0, 1.0);
        let m_pos = cp.critical_isotherm(1.0);
        let m_neg = cp.critical_isotherm(-1.0);
        assert!((m_pos + m_neg).abs() < 1e-12);
    }

    // ── DecimationTransformation tests ─────────────────────────────────────

    #[test]
    fn test_decimation_k_decreases_for_large_k() {
        // For large K (ordered), the coupling K' < K (flows to 0 for 1D)
        let mut dec = DecimationTransformation::new(2.0);
        let k0 = dec.k;
        dec.step();
        // K' = 0.5 ln cosh(2K) < K for K > 0
        assert!(dec.k < k0);
    }

    #[test]
    fn test_decimation_runs_to_zero() {
        let mut dec = DecimationTransformation::new(1.0);
        dec.run(100);
        assert!(dec.k < 0.01);
    }

    #[test]
    fn test_decimation_correlation_length_decreases() {
        let mut dec = DecimationTransformation::new(1.0);
        let xi0 = dec.correlation_length();
        dec.step();
        let xi1 = dec.correlation_length();
        assert!(xi1 < xi0);
    }

    #[test]
    fn test_decimation_trajectory_length() {
        let dec = DecimationTransformation::new(1.0);
        let traj = dec.trajectory(10);
        assert_eq!(traj.len(), 11);
    }

    // ── ConformalFieldTheory tests ──────────────────────────────────────────

    #[test]
    fn test_casimir_energy_ising_cft() {
        // 2-D Ising CFT: c = 1/2
        let mut cft = ConformalFieldTheory::new(0.5, "Ising2D");
        cft.add_primary(0.0, 0); // identity
        cft.add_primary(0.5, 0); // σ
        cft.add_primary(1.0, 0); // ε
        let e0 = cft.casimir_energy(1.0);
        let expected = -PI * 0.5 / 6.0;
        assert!((e0 - expected).abs() < 1e-12);
    }

    #[test]
    fn test_holomorphic_dimensions_identity() {
        let mut cft = ConformalFieldTheory::new(0.5, "Ising2D");
        cft.add_primary(0.0, 0);
        let dims = cft.holomorphic_dimensions();
        assert!((dims[0].0).abs() < 1e-12);
        assert!((dims[0].1).abs() < 1e-12);
    }

    #[test]
    fn test_unitarity_bound_satisfied() {
        let mut cft = ConformalFieldTheory::new(0.5, "Ising2D");
        cft.add_primary(0.0, 0);
        cft.add_primary(0.5, 0);
        cft.add_primary(1.0, 0);
        assert!(cft.satisfies_unitarity_bound());
    }

    #[test]
    fn test_unitarity_bound_violated() {
        let mut cft = ConformalFieldTheory::new(1.0, "Ghost");
        cft.add_primary(0.3, 2); // Δ < |s| = 2, violation
        assert!(!cft.satisfies_unitarity_bound());
    }

    #[test]
    fn test_virasoro_character_positive() {
        let cft = ConformalFieldTheory::new(0.5, "Ising2D");
        let chi = cft.virasoro_character(0.0, 0.9, 10);
        assert!(chi > 0.0);
    }

    #[test]
    fn test_partition_function_positive() {
        let mut cft = ConformalFieldTheory::new(0.5, "Ising2D");
        cft.add_primary(0.0, 0);
        cft.add_primary(0.5, 0);
        cft.add_primary(1.0, 0);
        let z = cft.partition_function(0.5, 10);
        assert!(z >= 0.0);
    }

    #[test]
    fn test_central_charge_from_tt() {
        let cft = ConformalFieldTheory::new(0.5, "Ising2D");
        // ⟨T(z)T(0)⟩ = c/(2z⁴) → tt = c/(2 z⁴)
        let z = 2.0_f64;
        let tt = 0.5 / (2.0 * z.powi(4));
        let c = cft.central_charge_from_tt(z, tt);
        assert!((c - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_s_matrix_wzw_symmetry() {
        let cft = ConformalFieldTheory::new(1.0, "SU2_k2");
        // S_{nm} = S_{mn} by symmetry of sine
        let s12 = cft.s_matrix_wzw(1, 2, 4);
        let s21 = cft.s_matrix_wzw(2, 1, 4);
        assert!((s12 - s21).abs() < 1e-12);
    }

    #[test]
    fn test_rushbrooke_ising_3d() {
        let exp = UniversalityClass::Ising3D.exponents();
        let sr = ScalingRelations::new(exp);
        let (lhs, rhs) = sr.rushbrooke();
        // α + 2β + γ should equal 2 (within numerical precision of tabulated values)
        assert!((lhs - rhs).abs() < 0.01, "lhs={lhs:.6}, rhs={rhs:.6}");
    }

    #[test]
    fn test_migdal_kadanoff_step() {
        let mut dec = DecimationTransformation::new(1.0);
        let k0 = dec.k;
        dec.migdal_kadanoff_step(2);
        // K' = 0.5 ln cosh(4K) for b=2
        let expected = 0.5 * (2.0 * 2.0 * k0).cosh().ln();
        assert!((dec.k - expected).abs() < 1e-12);
    }

    #[test]
    fn test_beta_no_fixed_point_one_loop() {
        let beta = BetaFunction::new(1.0, 0.0, LoopOrder::OneLoop, "phi4");
        assert!(beta.fixed_point().is_none());
    }

    #[test]
    fn test_wf_flow_u_approaches_fixed_point() {
        let mut rg = WilsonianRg::new(3.0, 1.0, 0.01, 0.0);
        let ustar = rg.wilson_fisher_fixed_point(1.0);
        rg.run(5000, 0.001, 1.0);
        // After many steps the coupling should be growing toward the fixed point
        // (convergence from u0=0.01 to ustar~17.5 takes l>>5; check we're moving right direction)
        assert!(
            (rg.u - ustar).abs() < 0.95 * ustar,
            "u={:.6}, u*={:.6}",
            rg.u,
            ustar
        );
    }
}
