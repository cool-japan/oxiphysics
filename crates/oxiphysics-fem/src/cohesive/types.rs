//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

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
/// PPR (Park-Paulino-Roesler, 2009) potential-based cohesive zone model.
///
/// The PPR potential unifies mode-I and mode-II with a polynomial form:
///
/// Ψ = Γ_n/m * (m/n * (1 - Δ_n/δ_n)^m - 1)^n * (Γ_t/n * (n/m * (1 - |Δ_t|/δ_t)^n - 1)^m + Γ_n)
///
/// Simplified implementation: independent polynomial traction-separation law.
///
/// References: Park, Paulino, Roesler (2009) J. Mech. Phys. Solids 57(6).
#[derive(Debug, Clone, Copy)]
pub struct PprParams {
    /// Mode-I fracture energy Γ_n (J/m²).
    pub gamma_n: f64,
    /// Mode-II fracture energy Γ_t (J/m²).
    pub gamma_t: f64,
    /// Shape parameter for normal direction (controls softening slope) α.
    pub alpha: f64,
    /// Shape parameter for tangential direction β.
    pub beta: f64,
    /// Critical normal separation δ_n (m) where traction → 0.
    pub delta_n: f64,
    /// Critical tangential separation δ_t (m) where traction → 0.
    pub delta_t: f64,
    /// Normal initial slope indicator λ_n ∈ (0,1).
    pub lambda_n: f64,
    /// Tangential initial slope indicator λ_t ∈ (0,1).
    pub lambda_t: f64,
}
impl PprParams {
    /// Create a new PPR parameter set.
    #[allow(clippy::too_many_arguments)]
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
    /// Separation at peak normal traction: δ_n^peak = λ_n * δ_n.
    pub fn delta_n_peak(&self) -> f64 {
        self.lambda_n * self.delta_n
    }
    /// Separation at peak tangential traction: δ_t^peak = λ_t * δ_t.
    pub fn delta_t_peak(&self) -> f64 {
        self.lambda_t * self.delta_t
    }
    /// Peak normal traction σ_max (Pa).
    ///
    /// Derived from the PPR potential as: σ_max = Γ_n / δ_n * (α/λ_n) * (α-1)^(α-1) ...
    /// Simplified polynomial approximation:
    pub fn sigma_max(&self) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let lam = self.lambda_n.clamp(f64::EPSILON, 1.0 - f64::EPSILON);
        self.gamma_n / dn * self.alpha * (1.0 - lam).powf(self.alpha - 1.0) / lam.powf(self.alpha)
    }
    /// Evaluate the normal traction T_n at separation Δ_n (mode-I, Δ_t=0).
    ///
    /// Uses polynomial traction-separation shape:
    /// T_n = σ_max * (Δ_n/δ_n^p)^(α-1) * (1 - Δ_n/δ_n)^α   for 0 ≤ Δ_n ≤ δ_n
    pub fn normal_traction(&self, delta_n: f64) -> f64 {
        if delta_n <= 0.0 || delta_n >= self.delta_n {
            return 0.0;
        }
        let s = delta_n / self.delta_n;
        let sp = self.lambda_n.clamp(f64::EPSILON, 1.0 - f64::EPSILON);
        let t_max = self.sigma_max();
        let numerator = (s / sp).powf(self.alpha - 1.0) * (1.0 - s).powf(self.alpha);
        let denom = (1.0 - sp).powf(self.alpha);
        t_max * numerator / denom.max(f64::EPSILON)
    }
    /// Evaluate the mode-I fracture energy by integrating T_n dΔ_n numerically (n=200 pts).
    pub fn mode_i_energy_numerical(&self) -> f64 {
        let n = 200usize;
        let dn_max = self.delta_n;
        let h = dn_max / n as f64;
        let mut sum = 0.0;
        for i in 0..n {
            let d = (i as f64 + 0.5) * h;
            sum += self.normal_traction(d) * h;
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
/// The potential is:
/// Φ = φ_n * \[ 1 - (1 + Δ_n/δ_n) * exp(-Δ_n/δ_n) * exp(-Δ_t²/δ_t²) \]
///
/// where φ_n is the normal work of separation, δ_n and δ_t are
/// characteristic lengths, and Δ_n, Δ_t are normal and tangential separations.
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
    /// Normal traction T_n at (Δ_n, Δ_t).
    pub fn normal_traction(&self, delta_n: f64, delta_t: f64) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let dt = self.delta_t.max(f64::EPSILON);
        let exp_n = (-delta_n / dn).exp();
        let exp_t = (-(delta_t * delta_t) / (dt * dt)).exp();
        (self.phi_n / dn) * exp_n * exp_t * (1.0 + delta_n / dn)
    }
    /// Tangential traction T_t at (Δ_n, Δ_t) — magnitude in the tangential plane.
    pub fn tangential_traction(&self, delta_n: f64, delta_t: f64) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let dt = self.delta_t.max(f64::EPSILON);
        let exp_n = (-delta_n / dn).exp();
        let exp_t = (-(delta_t * delta_t) / (dt * dt)).exp();
        2.0 * (self.phi_t / (dt * dt)) * delta_t * exp_n * exp_t
    }
    /// Evaluate the full cohesive potential Φ at (Δ_n, Δ_t).
    pub fn potential(&self, delta_n: f64, delta_t: f64) -> f64 {
        let dn = self.delta_n.max(f64::EPSILON);
        let dt = self.delta_t.max(f64::EPSILON);
        self.phi_n
            * (1.0
                - (1.0 + delta_n / dn)
                    * (-delta_n / dn).exp()
                    * (-(delta_t * delta_t) / (dt * dt)).exp())
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
