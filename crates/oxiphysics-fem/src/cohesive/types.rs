//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use crate::error::Error;

/// Parameters for the bilinear (linear-softening) traction-separation law.
///
/// The law has two phases:
/// 1. Linear loading:  T = K * delta            for delta < delta_0
/// 2. Linear softening: T = T_max * (delta_f - delta) / (delta_f - delta_0),
///    for delta_0 <= delta < delta_f
/// 3. Complete failure: T = 0                   for delta >= delta_f
#[derive(Debug, Clone, Copy)]
pub struct BilinearCzmParams {
    /// Penalty stiffness (N/m³) – slope of the initial linear branch.
    pub penalty_stiffness: f64,
    /// Peak traction (Pa).
    pub traction_max: f64,
    /// Critical (failure) separation (m) at which traction → 0.
    pub delta_failure: f64,
}
impl BilinearCzmParams {
    /// Create a new `BilinearCzmParams`.
    pub fn new(penalty_stiffness: f64, traction_max: f64, delta_failure: f64) -> Self {
        Self {
            penalty_stiffness,
            traction_max,
            delta_failure,
        }
    }
    /// Onset (damage-initiation) separation: δ₀ = T_max / K.
    pub fn delta_onset(&self) -> f64 {
        self.traction_max / self.penalty_stiffness.max(f64::EPSILON)
    }
    /// Mode-I fracture toughness G_c = ½ T_max * δ_f (area under T-δ curve).
    pub fn fracture_toughness(&self) -> f64 {
        0.5 * self.traction_max * self.delta_failure
    }
    /// Evaluate the traction for a given scalar separation `delta` (≥ 0).
    ///
    /// Returns the traction magnitude (N/m²).
    pub fn traction(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let d0 = self.delta_onset();
        if delta < d0 {
            self.penalty_stiffness * delta
        } else if delta < self.delta_failure {
            let df = self.delta_failure;
            self.traction_max * (df - delta) / (df - d0).max(f64::EPSILON)
        } else {
            0.0
        }
    }
    /// Secant stiffness at separation `delta`: T(delta) / delta.
    ///
    /// Returns 0 when delta ≤ 0.
    pub fn secant_stiffness(&self, delta: f64) -> f64 {
        if delta <= f64::EPSILON {
            return self.penalty_stiffness;
        }
        self.traction(delta) / delta
    }
    /// Scalar damage variable D ∈ \[0, 1\] at separation `delta`.
    ///
    /// D = 0: undamaged, D = 1: fully failed.
    pub fn damage(&self, delta_max: f64) -> f64 {
        let d0 = self.delta_onset();
        let df = self.delta_failure;
        if delta_max <= d0 {
            0.0
        } else if delta_max >= df {
            1.0
        } else {
            (df * (delta_max - d0)) / (delta_max * (df - d0).max(f64::EPSILON))
        }
    }
}
/// Number of Gauss–Legendre panels used by [`PprParams::mode_i_energy_numerical`].
const PPR_ENERGY_PANELS: usize = 64;

/// Nodes and weights of the 5-point Gauss–Legendre rule on `[-1, 1]`.
const GAUSS_LEGENDRE_5: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
    (-0.538_469_310_105_683, 0.478_628_670_499_366_5),
    (0.0, 0.568_888_888_888_889),
    (0.538_469_310_105_683, 0.478_628_670_499_366_5),
    (0.906_179_845_938_664, 0.236_926_885_056_189_1),
];

/// PPR (Park-Paulino-Roesler, 2009) potential-based cohesive zone model.
///
/// The PPR potential is built from the polynomial
///
/// `f(s) = (1 − s)^α · (1 + α·s/m)^m`,  `s = Δ/δ`,
///
/// with the shape exponent `m = α(α − 1)·λ² / (1 − α·λ²)`, which makes `f`
/// fall from 1 at `s = 0` to 0 at `s = 1` and places the traction maximum
/// exactly at the initial-slope indicator `s = λ`.  For pure mode I (no
/// tangential opening) the potential is `Ψ = Γ_n·(1 − f(Δ_n/δ_n))`, so the
/// normal traction is its exact derivative
///
/// `T_n(Δ_n) = Γ_n/δ_n · α·(1 + α/m) · s·(1 − s)^(α−1)·(1 + α·s/m)^(m−1)`
///
/// and the work of separation `∫₀^δn T_n dΔ_n` equals the mode-I fracture
/// energy `Γ_n` for every admissible `(α, λ_n)`.  The peak traction is
/// [`sigma_max`](Self::sigma_max) at `Δ_n = λ_n·δ_n`
/// ([`delta_n_peak`](Self::delta_n_peak)).
///
/// The tangential parameters (`Γ_t`, `β`, `δ_t`, `λ_t`) define the mode-II law
/// of the same form; this type evaluates the mode-I traction and the
/// tangential peak position.
///
/// The model is only defined for `α > 1` and `0 < λ < 1/√α` (otherwise `m` is
/// not a positive finite number); [`validate`](Self::validate) checks the whole
/// parameter set, and the traction functions return `NaN` for parameters
/// outside the domain instead of inventing a value.
///
/// References: Park, Paulino, Roesler (2009) J. Mech. Phys. Solids 57(6).
#[derive(Debug, Clone, Copy)]
pub struct PprParams {
    /// Mode-I fracture energy Γ_n (J/m²).
    pub gamma_n: f64,
    /// Mode-II fracture energy Γ_t (J/m²).
    pub gamma_t: f64,
    /// Shape parameter for normal direction (controls softening slope) α > 1.
    pub alpha: f64,
    /// Shape parameter for tangential direction β > 1.
    pub beta: f64,
    /// Critical normal separation δ_n (m) where traction → 0.
    pub delta_n: f64,
    /// Critical tangential separation δ_t (m) where traction → 0.
    pub delta_t: f64,
    /// Normal initial slope indicator λ_n = δ_nc/δ_n ∈ (0, 1/√α): the peak
    /// traction occurs at Δ_n = λ_n·δ_n.
    pub lambda_n: f64,
    /// Tangential initial slope indicator λ_t ∈ (0, 1/√β).
    pub lambda_t: f64,
}
impl PprParams {
    /// Create a new PPR parameter set.
    ///
    /// The values are stored as given; use [`try_new`](Self::try_new) or
    /// [`validate`](Self::validate) to check that they lie in the PPR domain.
    pub fn new(
        gamma_n: f64,
        gamma_t: f64,
        alpha: f64,
        beta: f64,
        delta_n: f64,
        delta_t: f64,
        lambda_n: f64,
        lambda_t: f64,
    ) -> Self {
        Self {
            gamma_n,
            gamma_t,
            alpha,
            beta,
            delta_n,
            delta_t,
            lambda_n,
            lambda_t,
        }
    }
    /// Create a PPR parameter set and check it with [`validate`](Self::validate).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] for parameters outside the PPR domain.
    pub fn try_new(
        gamma_n: f64,
        gamma_t: f64,
        alpha: f64,
        beta: f64,
        delta_n: f64,
        delta_t: f64,
        lambda_n: f64,
        lambda_t: f64,
    ) -> Result<Self, Error> {
        let params = Self::new(
            gamma_n, gamma_t, alpha, beta, delta_n, delta_t, lambda_n, lambda_t,
        );
        params.validate()?;
        Ok(params)
    }
    /// Check that the parameters lie in the PPR domain: finite positive
    /// fracture energies and critical separations, shape parameters `α, β > 1`
    /// and slope indicators `0 < λ_n < 1/√α`, `0 < λ_t < 1/√β` (so that both
    /// shape exponents `m`, `n` are positive and finite).
    ///
    /// # Errors
    ///
    /// [`Error::InvalidInput`] naming the first offending parameter.
    pub fn validate(&self) -> Result<(), Error> {
        let positive = |name: &str, value: f64| -> Result<(), Error> {
            if value.is_finite() && value > 0.0 {
                Ok(())
            } else {
                Err(Error::InvalidInput(format!(
                    "PPR parameter {name} must be finite and positive, got {value}"
                )))
            }
        };
        positive("gamma_n", self.gamma_n)?;
        positive("gamma_t", self.gamma_t)?;
        positive("delta_n", self.delta_n)?;
        positive("delta_t", self.delta_t)?;
        for (name, shape, lambda) in [
            ("alpha", self.alpha, self.lambda_n),
            ("beta", self.beta, self.lambda_t),
        ] {
            if !(shape.is_finite() && shape > 1.0) {
                return Err(Error::InvalidInput(format!(
                    "PPR shape parameter {name} must be finite and greater than 1, got {shape}"
                )));
            }
            if !(lambda.is_finite() && lambda > 0.0 && shape * lambda * lambda < 1.0) {
                return Err(Error::InvalidInput(format!(
                    "PPR slope indicator for {name} = {shape} must satisfy 0 < λ < 1/√{name} \
                     (= {}), got {lambda}",
                    1.0 / shape.sqrt()
                )));
            }
        }
        Ok(())
    }
    /// Shape exponent `m = α(α − 1)·λ_n² / (1 − α·λ_n²)` of the mode-I law, or
    /// `None` when `α` and `λ_n` are outside the PPR domain.
    pub fn normal_shape_exponent(&self) -> Option<f64> {
        let alpha = self.alpha;
        let lambda = self.lambda_n;
        if !(alpha.is_finite() && alpha > 1.0 && lambda.is_finite() && lambda > 0.0) {
            return None;
        }
        let denominator = 1.0 - alpha * lambda * lambda;
        if denominator > 0.0 {
            Some(alpha * (alpha - 1.0) * lambda * lambda / denominator)
        } else {
            None
        }
    }
    /// Mode-I traction at the normalised opening `s ∈ (0, 1)` for the shape
    /// exponent `m` (see the type-level documentation).
    fn mode_i_shape_traction(&self, s: f64, m: f64) -> f64 {
        let alpha = self.alpha;
        // (1 + α s/m)^(m−1) through ln_1p stays accurate for very large m.
        let growth = ((m - 1.0) * (alpha * s / m).ln_1p()).exp();
        self.gamma_n / self.delta_n
            * alpha
            * (1.0 + alpha / m)
            * s
            * (1.0 - s).powf(alpha - 1.0)
            * growth
    }
    /// Separation at peak normal traction: δ_n^peak = λ_n * δ_n.
    pub fn delta_n_peak(&self) -> f64 {
        self.lambda_n * self.delta_n
    }
    /// Separation at peak tangential traction: δ_t^peak = λ_t * δ_t.
    pub fn delta_t_peak(&self) -> f64 {
        self.lambda_t * self.delta_t
    }
    /// Peak normal traction σ_max (Pa), attained at Δ_n = λ_n·δ_n:
    ///
    /// `σ_max = Γ_n/δ_n · α·(1 + α/m) · λ_n·(1 − λ_n)^(α−1)·(1 + α·λ_n/m)^(m−1)`.
    ///
    /// `NaN` for parameters outside the PPR domain (see [`validate`](Self::validate)).
    pub fn sigma_max(&self) -> f64 {
        match (self.normal_shape_exponent(), self.validate()) {
            (Some(m), Ok(())) => self.mode_i_shape_traction(self.lambda_n, m),
            _ => f64::NAN,
        }
    }
    /// Evaluate the normal traction T_n at separation Δ_n (mode I, Δ_t = 0):
    ///
    /// `T_n = Γ_n/δ_n · α·(1 + α/m) · s·(1 − s)^(α−1)·(1 + α·s/m)^(m−1)`,
    /// `s = Δ_n/δ_n`, for `0 < Δ_n < δ_n`, and 0 outside that range.
    ///
    /// This is the exact derivative of the PPR mode-I potential, so
    /// `∫ T_n dΔ_n = Γ_n`, and its maximum [`sigma_max`](Self::sigma_max) lies at
    /// `Δ_n = λ_n·δ_n`.  `NaN` for parameters outside the PPR domain (see
    /// [`validate`](Self::validate)).
    pub fn normal_traction(&self, delta_n: f64) -> f64 {
        let m = match (self.normal_shape_exponent(), self.validate()) {
            (Some(m), Ok(())) => m,
            _ => return f64::NAN,
        };
        if delta_n.is_nan() {
            return f64::NAN;
        }
        if delta_n <= 0.0 || delta_n >= self.delta_n {
            return 0.0;
        }
        self.mode_i_shape_traction(delta_n / self.delta_n, m)
    }
    /// Evaluate the mode-I fracture energy by integrating T_n dΔ_n numerically
    /// (64 panels of 5-point Gauss–Legendre quadrature); it reproduces the
    /// fracture energy `gamma_n` to better than `1e-4` relative.
    ///
    /// `NaN` for parameters outside the PPR domain.
    pub fn mode_i_energy_numerical(&self) -> f64 {
        let panel_width = self.delta_n / PPR_ENERGY_PANELS as f64;
        let mut sum = 0.0;
        for panel in 0..PPR_ENERGY_PANELS {
            let start = panel as f64 * panel_width;
            for &(node, weight) in &GAUSS_LEGENDRE_5 {
                let d = start + 0.5 * panel_width * (node + 1.0);
                sum += 0.5 * panel_width * weight * self.normal_traction(d);
            }
        }
        sum
    }
}
/// Scalar damage evolution tracking for a cohesive element.
///
/// Tracks the maximum separation `delta_max` experienced and prevents
/// healing (damage is irreversible).
#[derive(Debug, Clone)]
pub struct CohesiveDamageState {
    /// Current scalar damage D ∈ \[0, 1\].
    pub damage: f64,
    /// Maximum effective separation ever experienced (m).
    pub delta_max: f64,
    /// Whether the element is fully failed (D = 1).
    pub failed: bool,
}
impl CohesiveDamageState {
    /// Create a new undamaged state.
    pub fn new() -> Self {
        Self {
            damage: 0.0,
            delta_max: 0.0,
            failed: false,
        }
    }
    /// Update the damage state given the current effective separation `delta`
    /// and the bilinear cohesive parameters.
    ///
    /// Returns the updated damage variable.
    pub fn update(&mut self, delta: f64, params: &BilinearCzmParams) -> f64 {
        if self.failed {
            return 1.0;
        }
        if delta > self.delta_max {
            self.delta_max = delta;
        }
        self.damage = params.damage(self.delta_max);
        if self.damage >= 1.0 {
            self.damage = 1.0;
            self.failed = true;
        }
        self.damage
    }
    /// Compute the current traction given effective separation `delta`.
    ///
    /// Uses the bilinear law with irreversibility (D is monotonically increasing).
    pub fn traction(&self, delta: f64, params: &BilinearCzmParams) -> f64 {
        if self.failed || delta <= 0.0 {
            return 0.0;
        }
        (1.0 - self.damage) * params.penalty_stiffness * delta
    }
}
/// Descriptor for a cohesive interface to be inserted between two element faces.
///
/// Cohesive element insertion duplicates nodes along the interface and
/// inserts zero-thickness cohesive elements between adjacent bulk elements.
#[derive(Debug, Clone)]
pub struct CohesiveInterface {
    /// Global indices of the nodes on face A (first bulk element face).
    pub nodes_a: Vec<usize>,
    /// Global indices of the nodes on face B (second bulk element face, duplicate).
    pub nodes_b: Vec<usize>,
    /// Material parameters for this interface.
    pub params: BilinearCzmParams,
    /// Current damage state for each interface (one per integration point).
    pub damage_states: Vec<CohesiveDamageState>,
}
impl CohesiveInterface {
    /// Create a new cohesive interface from two matching node sets.
    ///
    /// `nodes_a` and `nodes_b` must have the same length (the shared face topology).
    pub fn new(nodes_a: Vec<usize>, nodes_b: Vec<usize>, params: BilinearCzmParams) -> Self {
        let n_ip = nodes_a.len();
        Self {
            nodes_a,
            nodes_b,
            params,
            damage_states: vec![CohesiveDamageState::new(); n_ip],
        }
    }
    /// Number of nodes on each face.
    pub fn n_nodes(&self) -> usize {
        self.nodes_a.len()
    }
    /// Check whether all integration points have failed.
    pub fn is_fully_failed(&self) -> bool {
        self.damage_states.iter().all(|s| s.failed)
    }
    /// Maximum damage across all integration points.
    pub fn max_damage(&self) -> f64 {
        self.damage_states
            .iter()
            .map(|s| s.damage)
            .fold(0.0_f64, f64::max)
    }
}
/// Thermo-mechanical cohesive zone model.
///
/// Accounts for thermal effects on the cohesive strength via:
/// T_eff(δ, T) = T_coh(δ) * f_thermal(T)
///
/// where f_thermal reduces strength at elevated temperature.
#[derive(Debug, Clone, Copy)]
pub struct ThermalCzmParams {
    /// Bilinear cohesive law at reference temperature.
    pub base: BilinearCzmParams,
    /// Reference temperature T_ref (K or °C).
    pub t_ref: f64,
    /// Thermal softening exponent m: f(T) = (1 - (T-T_ref)/T_melt)^m.
    pub softening_exponent: f64,
    /// Melt/degradation temperature T_melt (relative to T_ref).
    pub t_melt: f64,
}
impl ThermalCzmParams {
    /// Create a new thermo-mechanical CZM parameter set.
    pub fn new(base: BilinearCzmParams, t_ref: f64, softening_exponent: f64, t_melt: f64) -> Self {
        Self {
            base,
            t_ref,
            softening_exponent,
            t_melt,
        }
    }
    /// Thermal reduction factor f(T) ∈ \[0, 1\].
    ///
    /// f(T) = max(0, 1 - ((T - T_ref) / T_melt)^m)
    pub fn thermal_factor(&self, temperature: f64) -> f64 {
        let dt = temperature - self.t_ref;
        if dt <= 0.0 {
            return 1.0;
        }
        let ratio = dt / self.t_melt.max(f64::EPSILON);
        (1.0 - ratio.powf(self.softening_exponent)).max(0.0)
    }
    /// Traction at separation `delta` and temperature `temp`.
    pub fn traction(&self, delta: f64, temp: f64) -> f64 {
        self.base.traction(delta) * self.thermal_factor(temp)
    }
    /// Thermal effective fracture toughness G_c(T) = G_c_ref * f(T).
    pub fn fracture_toughness(&self, temp: f64) -> f64 {
        self.base.fracture_toughness() * self.thermal_factor(temp)
    }
    /// Thermal contribution to traction-separation: penalty heat flux through interface.
    ///
    /// Simplified interface conductance model: q = h_int * ΔT
    /// where h_int is the interface thermal conductance (W/m²K) and ΔT is the
    /// temperature jump across the interface.
    pub fn heat_flux(&self, h_int: f64, delta_t: f64) -> f64 {
        h_int * delta_t
    }
}
/// Xu-Needleman exponential traction-separation law (1993).
///
/// The potential (coupling parameters `q = φ_t/φ_n`, `r = 0`) is
///
/// `Φ = φ_n − (1 + Δ_n/δ_n)·exp(−Δ_n/δ_n)·[(φ_n − φ_t) + φ_t·exp(−Δ_t²/δ_t²)]`
///
/// which for equal works of separation `φ_t = φ_n` is the familiar
/// `Φ = φ_n·[1 − (1 + Δ_n/δ_n)·exp(−Δ_n/δ_n)·exp(−Δ_t²/δ_t²)]`.  φ_n and φ_t
/// are the normal and tangential works of separation, δ_n and δ_t the
/// characteristic lengths, and Δ_n, Δ_t the normal and tangential separations.
/// The tractions are the exact derivatives of Φ:
///
/// `T_n = ∂Φ/∂Δ_n = (1/δ_n)·(Δ_n/δ_n)·exp(−Δ_n/δ_n)·[(φ_n − φ_t) + φ_t·exp(−Δ_t²/δ_t²)]`
///
/// `T_t = ∂Φ/∂Δ_t = (2φ_t/δ_t²)·Δ_t·(1 + Δ_n/δ_n)·exp(−Δ_n/δ_n)·exp(−Δ_t²/δ_t²)`
///
/// so the normal traction vanishes at zero opening and peaks at `Δ_n = δ_n`
/// with [`sigma_max`](Self::sigma_max).
#[derive(Debug, Clone, Copy)]
pub struct XuNeedlemanParams {
    /// Normal work of separation φ_n = σ_max * e * δ_n  (J/m²).
    pub phi_n: f64,
    /// Tangential work of separation φ_t (J/m²).
    pub phi_t: f64,
    /// Characteristic normal length δ_n (m).
    pub delta_n: f64,
    /// Characteristic tangential length δ_t (m).
    pub delta_t: f64,
}
impl XuNeedlemanParams {
    /// Create a new `XuNeedlemanParams`.
    pub fn new(phi_n: f64, phi_t: f64, delta_n: f64, delta_t: f64) -> Self {
        Self {
            phi_n,
            phi_t,
            delta_n,
            delta_t,
        }
    }
    /// Peak normal traction σ_max = φ_n / (e * δ_n).
    pub fn sigma_max(&self) -> f64 {
        self.phi_n / (std::f64::consts::E * self.delta_n.max(f64::EPSILON))
    }
    /// Peak tangential traction τ_max = φ_t / (sqrt(e/2) * δ_t).
    pub fn tau_max(&self) -> f64 {
        self.phi_t / ((std::f64::consts::E / 2.0).sqrt() * self.delta_t.max(f64::EPSILON))
    }
    /// Normal traction `T_n = ∂Φ/∂Δ_n` at (Δ_n, Δ_t); zero at zero opening.
    pub fn normal_traction(&self, delta_n: f64, delta_t: f64) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let dt = self.delta_t.max(f64::EPSILON);
        let s = delta_n / dn;
        let exp_t = (-(delta_t * delta_t) / (dt * dt)).exp();
        let coupling = (self.phi_n - self.phi_t) + self.phi_t * exp_t;
        s * (-s).exp() * coupling / dn
    }
    /// Tangential traction `T_t = ∂Φ/∂Δ_t` at (Δ_n, Δ_t) — magnitude in the
    /// tangential plane.
    pub fn tangential_traction(&self, delta_n: f64, delta_t: f64) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let dt = self.delta_t.max(f64::EPSILON);
        let s = delta_n / dn;
        let exp_t = (-(delta_t * delta_t) / (dt * dt)).exp();
        2.0 * (self.phi_t / (dt * dt)) * delta_t * (1.0 + s) * (-s).exp() * exp_t
    }
    /// Evaluate the full cohesive potential Φ at (Δ_n, Δ_t).
    pub fn potential(&self, delta_n: f64, delta_t: f64) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let dt = self.delta_t.max(f64::EPSILON);
        let s = delta_n / dn;
        let exp_t = (-(delta_t * delta_t) / (dt * dt)).exp();
        self.phi_n - (1.0 + s) * (-s).exp() * ((self.phi_n - self.phi_t) + self.phi_t * exp_t)
    }
}
/// Fatigue cohesive zone model parameters (Paris-type interface fatigue).
#[derive(Debug, Clone, Copy)]
pub struct FatigueCzmParams {
    /// Monotonic cohesive parameters.
    pub mono: BilinearCzmParams,
    /// Paris fatigue constant C (relates Δδ^m_fatigue to dD/dN).
    pub c_fatigue: f64,
    /// Paris fatigue exponent m_fatigue.
    pub m_fatigue: f64,
    /// Fatigue threshold: no fatigue damage below this separation range.
    pub delta_threshold: f64,
}
impl FatigueCzmParams {
    /// Create new fatigue CZM parameters.
    pub fn new(
        mono: BilinearCzmParams,
        c_fatigue: f64,
        m_fatigue: f64,
        delta_threshold: f64,
    ) -> Self {
        Self {
            mono,
            c_fatigue,
            m_fatigue,
            delta_threshold,
        }
    }
    /// Number of cycles to failure at a constant separation range Δδ.
    ///
    /// N_f = (1 - D_0) / (C * Δδ^m_f / δ_f)
    pub fn cycles_to_failure(&self, delta_range: f64, d_initial: f64) -> Option<f64> {
        if delta_range <= self.delta_threshold {
            return None;
        }
        let dd_dn = self.c_fatigue * delta_range.powf(self.m_fatigue)
            / self.mono.delta_failure.max(f64::EPSILON);
        if dd_dn < f64::EPSILON {
            return None;
        }
        Some((1.0 - d_initial.clamp(0.0, 1.0)) / dd_dn)
    }
}
/// Frictional cohesive contact model for interfaces that can carry both
/// cohesive traction and frictional shear after crack initiation.
///
/// Uses a Coulomb friction law in the sliding regime:
/// |T_t| = μ * |T_n| + T_cohesive_shear
#[derive(Debug, Clone, Copy)]
pub struct FrictionalCzmParams {
    /// Cohesive mode-I parameters.
    pub normal: BilinearCzmParams,
    /// Cohesive mode-II (shear) parameters.
    pub shear: BilinearCzmParams,
    /// Coulomb friction coefficient μ.
    pub friction_coefficient: f64,
}
impl FrictionalCzmParams {
    /// Create new frictional cohesive zone parameters.
    pub fn new(
        normal: BilinearCzmParams,
        shear: BilinearCzmParams,
        friction_coefficient: f64,
    ) -> Self {
        Self {
            normal,
            shear,
            friction_coefficient,
        }
    }
    /// Normal traction at gap `gap` (compressive contact at gap < 0).
    pub fn normal_traction(&self, gap: f64, k_contact: f64) -> f64 {
        if gap < 0.0 {
            k_contact * gap.abs()
        } else {
            self.normal.traction(gap)
        }
    }
    /// Shear traction under mixed cohesive + frictional contribution.
    ///
    /// # Arguments
    /// * `delta_t` – tangential separation magnitude (m)
    /// * `gap` – normal gap (negative = contact)
    /// * `k_contact` – contact penalty stiffness (N/m³)
    pub fn shear_traction(&self, delta_t: f64, gap: f64, k_contact: f64) -> f64 {
        let t_coh = self.shear.traction(delta_t);
        if gap < 0.0 {
            let t_n = k_contact * gap.abs();
            t_coh + self.friction_coefficient * t_n
        } else {
            t_coh
        }
    }
    /// Maximum shear traction capacity (cohesive + frictional).
    pub fn max_shear_capacity(&self, t_normal_contact: f64) -> f64 {
        self.shear.traction_max + self.friction_coefficient * t_normal_contact
    }
}
/// Mixed-mode interaction criterion for cohesive zone initiation/failure.
///
/// Uses the quadratic interaction:
/// (T_n / T_n^c)² + (T_t / T_t^c)² + (T_3 / T_3^c)² = 1
///
/// where T_n^c, T_t^c, T_3^c are critical tractions in each mode.
#[derive(Debug, Clone, Copy)]
pub struct MixedModeCriterion {
    /// Critical normal traction T_n^c (Pa).
    pub t_n_crit: f64,
    /// Critical tangential (mode II) traction T_t_crit (Pa).
    pub t_t_crit: f64,
    /// Critical mode-III traction T_3_crit (Pa).
    pub t_3_crit: f64,
}
impl MixedModeCriterion {
    /// Create a new `MixedModeCriterion`.
    pub fn new(t_n_crit: f64, t_t_crit: f64, t_3_crit: f64) -> Self {
        Self {
            t_n_crit,
            t_t_crit,
            t_3_crit,
        }
    }
    /// Evaluate the quadratic mixed-mode index f ∈ \[0, 1\].
    ///
    /// f = 1 means the criterion is exactly met; f < 1 means no failure.
    /// Only compressive (negative) normal tractions do not contribute.
    pub fn interaction_index(&self, t_n: f64, t_t: f64, t_3: f64) -> f64 {
        let tn_contrib = if t_n > 0.0 {
            (t_n / self.t_n_crit.max(f64::EPSILON)).powi(2)
        } else {
            0.0
        };
        let tt_contrib = (t_t / self.t_t_crit.max(f64::EPSILON)).powi(2);
        let t3_contrib = (t_3 / self.t_3_crit.max(f64::EPSILON)).powi(2);
        (tn_contrib + tt_contrib + t3_contrib).sqrt()
    }
    /// Check whether the initiation criterion is met (f ≥ 1).
    pub fn is_initiated(&self, t_n: f64, t_t: f64, t_3: f64) -> bool {
        self.interaction_index(t_n, t_t, t_3) >= 1.0
    }
    /// Mixed-mode ratio β = G_II / (G_I + G_II) (2D, mode I+II only).
    ///
    /// β = 0: pure mode I; β = 1: pure mode II.
    pub fn mixed_mode_ratio(g1: f64, g2: f64) -> f64 {
        let g_total = g1 + g2;
        if g_total < f64::EPSILON {
            return 0.0;
        }
        g2 / g_total
    }
}
/// Rate-dependent cohesive law with viscous regularization.
///
/// The effective traction is augmented by a viscous damper in parallel:
/// T_visc = T_coh(δ) + η * dδ/dt
///
/// where η is the viscosity parameter (Pa·s/m) and dδ/dt is the separation rate.
#[derive(Debug, Clone, Copy)]
pub struct ViscousCzmParams {
    /// Underlying bilinear cohesive parameters.
    pub base: BilinearCzmParams,
    /// Viscosity η (Pa·s/m).
    pub viscosity: f64,
}
impl ViscousCzmParams {
    /// Create a new viscous CZM parameter set.
    pub fn new(base: BilinearCzmParams, viscosity: f64) -> Self {
        Self { base, viscosity }
    }
    /// Evaluate the total traction at separation `delta` and separation rate `delta_dot`.
    pub fn traction(&self, delta: f64, delta_dot: f64) -> f64 {
        self.base.traction(delta) + self.viscosity * delta_dot
    }
    /// Evaluate the viscous contribution only.
    pub fn viscous_traction(&self, delta_dot: f64) -> f64 {
        self.viscosity * delta_dot
    }
    /// Compute the viscous energy dissipated over a time step Δt.
    ///
    /// Approximated as: W_v = η * (dδ/dt)² * Δt  (power dissipation * time).
    pub fn viscous_dissipation(&self, delta_dot: f64, dt: f64) -> f64 {
        self.viscosity * delta_dot * delta_dot * dt
    }
    /// Consistent tangent including viscous term: ∂T/∂δ + ∂T/∂δ̇ * (1/Δt).
    ///
    /// For implicit time integration with time step Δt:
    /// K_eff = K_coh + η/Δt
    pub fn consistent_tangent(&self, delta: f64, dt: f64) -> f64 {
        let k_coh = self.base.secant_stiffness(delta);
        k_coh + self.viscosity / dt.max(f64::EPSILON)
    }
}
/// R-curve model: fracture toughness increases with crack extension Δa.
///
/// The exponential R-curve model:
/// K_R(Δa) = K_0 + (K_ss - K_0) * (1 - exp(-Δa / a_0))
///
/// where:
/// - K_0  = initial (initiation) toughness
/// - K_ss = steady-state (plateau) toughness
/// - a_0  = characteristic crack extension scale
#[derive(Debug, Clone, Copy)]
pub struct RCurveModel {
    /// Initial fracture toughness K_0 (Pa√m).
    pub k_init: f64,
    /// Steady-state toughness K_ss (Pa√m).
    pub k_steady: f64,
    /// Characteristic crack extension scale a_0 (m).
    pub a_scale: f64,
}
impl RCurveModel {
    /// Create a new R-curve model.
    pub fn new(k_init: f64, k_steady: f64, a_scale: f64) -> Self {
        Self {
            k_init,
            k_steady,
            a_scale,
        }
    }
    /// Evaluate K_R at crack extension Δa.
    pub fn k_resistance(&self, delta_a: f64) -> f64 {
        let da = delta_a.max(0.0);
        let a0 = self.a_scale.max(f64::EPSILON);
        self.k_init + (self.k_steady - self.k_init) * (1.0 - (-da / a0).exp())
    }
    /// Check whether the applied SIF `k_applied` exceeds the resistance.
    ///
    /// Returns `true` if crack propagation occurs (K ≥ K_R).
    pub fn propagates(&self, k_applied: f64, delta_a: f64) -> bool {
        k_applied >= self.k_resistance(delta_a)
    }
    /// Compute the stable crack extension Δa for which K_applied = K_R.
    ///
    /// Uses a bisection search within `[0, da_max]`.
    /// Returns `None` if no equilibrium is found.
    pub fn stable_crack_extension(&self, k_applied: f64, da_max: f64, tol: f64) -> Option<f64> {
        if k_applied < self.k_init {
            return Some(0.0);
        }
        if k_applied >= self.k_steady {
            return None;
        }
        let mut lo = 0.0_f64;
        let mut hi = da_max;
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            let kr = self.k_resistance(mid);
            if (kr - k_applied).abs() < tol {
                return Some(mid);
            }
            if kr < k_applied {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }
}
/// Fatigue degradation state for cohesive zone model.
///
/// Tracks cyclic damage accumulation via:
/// D_f += ΔD_f per cycle, based on separation range Δδ.
///
/// The total damage variable is:
/// D_total = D_monotonic + D_fatigue (capped at 1)
#[derive(Debug, Clone)]
pub struct FatigueCzmState {
    /// Monotonic damage D_mono ∈ \[0, 1\].
    pub d_mono: f64,
    /// Fatigue damage D_fatigue ∈ \[0, 1\].
    pub d_fatigue: f64,
    /// Number of cycles accumulated.
    pub n_cycles: u64,
    /// Maximum separation in current cycle.
    pub delta_max_cycle: f64,
    /// Minimum separation in current cycle.
    pub delta_min_cycle: f64,
    /// Maximum historic separation.
    pub delta_historic_max: f64,
}
impl FatigueCzmState {
    /// Create a new undamaged fatigue CZM state.
    pub fn new() -> Self {
        Self {
            d_mono: 0.0,
            d_fatigue: 0.0,
            n_cycles: 0,
            delta_max_cycle: 0.0,
            delta_min_cycle: f64::MAX,
            delta_historic_max: 0.0,
        }
    }
    /// Total damage D = D_mono + D_fatigue, capped at 1.
    pub fn total_damage(&self) -> f64 {
        (self.d_mono + self.d_fatigue).min(1.0)
    }
    /// Whether the element is fully failed.
    pub fn is_failed(&self) -> bool {
        self.total_damage() >= 1.0
    }
    /// Record a new separation value within the current half-cycle.
    pub fn record_separation(&mut self, delta: f64) {
        if delta > self.delta_max_cycle {
            self.delta_max_cycle = delta;
        }
        if delta < self.delta_min_cycle {
            self.delta_min_cycle = delta;
        }
        if delta > self.delta_historic_max {
            self.delta_historic_max = delta;
        }
    }
    /// Increment the cycle count and compute fatigue damage for the cycle.
    ///
    /// Uses a Paris-type fatigue law for interfaces:
    /// dD/dN = C * (Δδ)^m / δ_f
    ///
    /// where Δδ = δ_max_cycle - δ_min_cycle is the separation range.
    pub fn advance_cycle(&mut self, c_fatigue: f64, m_fatigue: f64, delta_failure: f64) {
        self.n_cycles += 1;
        let delta_min = if self.delta_min_cycle == f64::MAX {
            0.0
        } else {
            self.delta_min_cycle
        };
        let delta_range = (self.delta_max_cycle - delta_min).max(0.0);
        let dd_dn = c_fatigue * delta_range.powf(m_fatigue) / delta_failure.max(f64::EPSILON);
        self.d_fatigue = (self.d_fatigue + dd_dn).min(1.0);
        self.delta_max_cycle = 0.0;
        self.delta_min_cycle = f64::MAX;
    }
    /// Update monotonic damage from current separation and bilinear params.
    pub fn update_monotonic(&mut self, delta: f64, params: &BilinearCzmParams) {
        let d_new = params.damage(delta);
        if d_new > self.d_mono {
            self.d_mono = d_new;
        }
    }
}
