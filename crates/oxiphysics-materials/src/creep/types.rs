//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Thermally activated rate model for creep.
///
/// Strain rate: eps_dot = nu_0 * exp(-ΔG / (k_B * T)) * sinh(V * sigma / (k_B * T))
///
/// where nu_0 is the attempt frequency, ΔG is the activation free energy,
/// V is the activation volume, and k_B is Boltzmann's constant.
pub struct ThermallyActivatedCreep {
    /// Attempt frequency ν₀ (1/s).
    pub nu0: f64,
    /// Activation free energy ΔG (J).
    pub delta_g: f64,
    /// Activation volume V (m³).
    pub activation_volume: f64,
    /// Boltzmann constant k_B (J/K).
    pub k_b: f64,
}
impl ThermallyActivatedCreep {
    /// Create a thermally activated creep model.
    pub fn new(nu0: f64, delta_g: f64, activation_volume: f64) -> Self {
        Self {
            nu0,
            delta_g,
            activation_volume,
            k_b: 1.380_649e-23,
        }
    }
    /// Strain rate at given stress and temperature.
    pub fn strain_rate(&self, stress: f64, temperature: f64) -> f64 {
        let kt = self.k_b * temperature;
        let boltzmann = (-self.delta_g / kt).exp();
        let sinh_term = (self.activation_volume * stress / kt).sinh();
        self.nu0 * boltzmann * sinh_term
    }
    /// Activation stress: σ* = k_B * T / V * arcsinh(eps_dot / nu_0 / exp(-ΔG/k_BT))
    pub fn activation_stress(&self, strain_rate: f64, temperature: f64) -> f64 {
        let kt = self.k_b * temperature;
        let boltzmann = (-self.delta_g / kt).exp();
        let arg = strain_rate / (self.nu0 * boltzmann);
        (kt / self.activation_volume) * arg.asinh()
    }
    /// Apparent activation energy at constant stress (from Arrhenius plot slope).
    pub fn apparent_activation_energy(&self, stress: f64) -> f64 {
        (self.delta_g - stress * self.activation_volume).max(0.0)
    }
}
/// Monkman-Grant creep life relation.
///
/// Empirical rule: ε_ss^m * t_f = C_MG
///
/// where ε_ss is the steady-state creep rate, t_f is the rupture life,
/// m is typically close to 1, and C_MG is a material constant.
pub struct MonkmanGrant {
    /// Monkman-Grant constant C_MG.
    pub c_mg: f64,
    /// Exponent m (typically ~1.0).
    pub m_exp: f64,
}
impl MonkmanGrant {
    /// Create a new Monkman-Grant model.
    pub fn new(c_mg: f64, m_exp: f64) -> Self {
        Self { c_mg, m_exp }
    }
    /// Predict rupture life from steady-state creep rate.
    ///
    /// t_f = C_MG / ε_ss^m
    pub fn rupture_life(&self, steady_state_rate: f64) -> f64 {
        if steady_state_rate <= 0.0 {
            return f64::INFINITY;
        }
        self.c_mg / steady_state_rate.powf(self.m_exp)
    }
    /// Back-calculate steady-state rate from known rupture life.
    ///
    /// ε_ss = (C_MG / t_f)^(1/m)
    pub fn steady_state_rate_from_life(&self, rupture_life: f64) -> f64 {
        if rupture_life <= 0.0 {
            return f64::INFINITY;
        }
        (self.c_mg / rupture_life).powf(1.0 / self.m_exp)
    }
    /// Check consistency: given ε_ss, compute t_f, then back-calculate ε_ss.
    pub fn is_consistent(&self, steady_state_rate: f64, tolerance: f64) -> bool {
        let t_f = self.rupture_life(steady_state_rate);
        let rate_back = self.steady_state_rate_from_life(t_f);
        (rate_back - steady_state_rate).abs() / steady_state_rate < tolerance
    }
}
/// Chaboche nonlinear kinematic hardening (back-stress evolution).
pub struct ChabocheKinematicHardening {
    /// Kinematic hardening modulus C
    pub c: f64,
    /// Dynamic recovery coefficient gamma
    pub gamma: f64,
}
impl ChabocheKinematicHardening {
    /// Create a new ChabocheKinematicHardening model.
    pub fn new(c: f64, gamma: f64) -> Self {
        Self { c, gamma }
    }
    /// Rate of back stress evolution.
    ///
    /// dot(X) = C * eps_dot_p - gamma * X * |eps_dot_p|
    pub fn back_stress_rate(&self, plastic_strain_rate: f64, back_stress: f64) -> f64 {
        self.c * plastic_strain_rate - self.gamma * back_stress * plastic_strain_rate.abs()
    }
    /// Explicit Euler update of back stress over a plastic strain increment dep.
    pub fn update_back_stress(&self, x: f64, dep: f64, dt: f64) -> f64 {
        let rate = self.back_stress_rate(dep / dt.max(f64::EPSILON), x);
        x + rate * dt
    }
    /// Saturation back stress: X_sat = C / gamma.
    pub fn saturation_stress(&self) -> f64 {
        if self.gamma.abs() > f64::EPSILON {
            self.c / self.gamma
        } else {
            f64::INFINITY
        }
    }
}
/// Primary creep using Bailey-Norton time-hardening: ε_p = A * σ^n * t^m.
#[allow(dead_code)]
pub struct PrimaryCreepModel {
    /// Constant A
    pub a: f64,
    /// Stress exponent n
    pub n: f64,
    /// Time exponent m (0 < m < 1 for primary hardening)
    pub m: f64,
}
impl PrimaryCreepModel {
    /// Create a primary creep model.
    pub fn new(a: f64, n: f64, m: f64) -> Self {
        Self { a, n, m }
    }
    /// Accumulated primary creep strain at stress `sigma` and time `t`.
    pub fn strain(&self, sigma: f64, t: f64) -> f64 {
        self.a * sigma.powf(self.n) * t.powf(self.m)
    }
    /// Primary creep rate: dε/dt = A * σ^n * m * t^(m-1).
    pub fn strain_rate(&self, sigma: f64, t: f64) -> f64 {
        if t <= 0.0 {
            return f64::INFINITY;
        }
        self.a * sigma.powf(self.n) * self.m * t.powf(self.m - 1.0)
    }
    /// Time at which primary rate equals the steady-state (Norton) rate.
    ///
    /// Returns the transition time from primary to secondary creep.
    pub fn transition_time(&self, norton: &NortonCreep, sigma: f64, temperature: f64) -> f64 {
        let eps_dot_ss = norton.creep_strain_rate(sigma, temperature);
        if eps_dot_ss <= 0.0 || self.m >= 1.0 {
            return f64::INFINITY;
        }
        let numerator = self.a * sigma.powf(self.n) * self.m;
        (numerator / eps_dot_ss).powf(1.0 / (1.0 - self.m))
    }
}
/// Generalised power-law creep: eps_dot = B * (sigma/sigma_0)^n.
///
/// No temperature dependence (use NortonCreep for Arrhenius activation).
pub struct PowerLawCreep {
    /// Reference strain rate B (1/s).
    pub b: f64,
    /// Reference stress sigma_0 (Pa).
    pub sigma_0: f64,
    /// Stress exponent n.
    pub n: f64,
}
impl PowerLawCreep {
    /// Create a new power-law creep model.
    pub fn new(b: f64, sigma_0: f64, n: f64) -> Self {
        Self { b, sigma_0, n }
    }
    /// Strain rate at given stress.
    pub fn strain_rate(&self, stress: f64) -> f64 {
        self.b * (stress / self.sigma_0).powf(self.n)
    }
    /// Creep strain after time t at constant stress.
    pub fn strain_at_time(&self, stress: f64, t: f64) -> f64 {
        self.strain_rate(stress) * t
    }

    /// Creep compliance J(t) = b + sigma_0 · t^n.
    ///
    /// Interprets the power-law parameters as a linear viscoelastic
    /// compliance function where `b` is the instantaneous compliance J₀,
    /// `sigma_0` is the time-dependent coefficient, and `n` is the time
    /// exponent.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        self.b + self.sigma_0 * t.powf(self.n)
    }

    /// Time derivative of creep compliance: dJ/dt = sigma_0 · n · t^(n−1).
    pub fn creep_rate(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return f64::INFINITY;
        }
        self.sigma_0 * self.n * t.powf(self.n - 1.0)
    }

    /// Creep strain at stress σ and time t: ε(t) = σ · J(t).
    pub fn creep_strain(&self, stress: f64, t: f64) -> f64 {
        stress * self.creep_compliance(t)
    }
}
/// Empirical three-stage creep curve.
pub struct CreepCurve {
    /// Primary creep coefficient A
    pub primary_rate: f64,
    /// Secondary (steady-state) creep rate (1/s)
    pub secondary_rate: f64,
    /// Strain at which tertiary stage begins
    pub tertiary_start_strain: f64,
}
impl CreepCurve {
    /// Create a new CreepCurve.
    pub fn new(primary_rate: f64, secondary_rate: f64, tertiary_start_strain: f64) -> Self {
        Self {
            primary_rate,
            secondary_rate,
            tertiary_start_strain,
        }
    }
    /// Total creep strain at time t.
    ///
    /// eps(t) = primary_rate * sqrt(t) + secondary_rate * t
    pub fn strain_at_time(&self, t: f64) -> f64 {
        self.primary_rate * t.sqrt() + self.secondary_rate * t
    }
    /// Instantaneous creep rate at time t.
    ///
    /// d(eps)/dt = primary_rate/(2*sqrt(t)) + secondary_rate
    pub fn strain_rate_at_time(&self, t: f64) -> f64 {
        if t < 1e-30 {
            return f64::INFINITY;
        }
        self.primary_rate / (2.0 * t.sqrt()) + self.secondary_rate
    }
    /// Find the time at which a target strain is first reached using bisection.
    pub fn time_to_strain(&self, target_strain: f64) -> Option<f64> {
        if target_strain <= 0.0 {
            return Some(0.0);
        }
        let t_max = if self.secondary_rate > 0.0 {
            target_strain / self.secondary_rate * 100.0
        } else if self.primary_rate > 0.0 {
            (target_strain / self.primary_rate).powi(2) * 10.0
        } else {
            return None;
        };
        if self.strain_at_time(t_max) < target_strain {
            return None;
        }
        let mut lo = 0.0_f64;
        let mut hi = t_max;
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if self.strain_at_time(mid) < target_strain {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(0.5 * (lo + hi))
    }
    /// Determine which creep stage we are in at a given strain level.
    pub fn stage_at_strain(&self, strain: f64) -> CreepStage {
        if strain < self.primary_rate * 100.0_f64.sqrt() {
            CreepStage::Primary { exponent: 0.5 }
        } else if strain < self.tertiary_start_strain {
            CreepStage::Secondary
        } else {
            CreepStage::Tertiary { acceleration: 2.0 }
        }
    }
}
/// Andrade primary creep model.
///
/// The Andrade equation describes the transient (primary) creep stage:
///
///   ε(t) = β * t^(1/3) + k * t
///
/// where β is the primary creep coefficient and k is the steady-state rate.
/// The first term dominates early (primary stage) and the second at longer times.
pub struct AndradeCreep {
    /// Primary creep coefficient β (dimensionless per s^(1/3)).
    pub beta: f64,
    /// Steady-state creep coefficient k (1/s).
    pub k_steady: f64,
}
impl AndradeCreep {
    /// Create a new Andrade creep model.
    pub fn new(beta: f64, k_steady: f64) -> Self {
        Self { beta, k_steady }
    }
    /// Total creep strain at time t.
    ///
    /// ε(t) = β * t^(1/3) + k * t
    pub fn strain(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        self.beta * t.powf(1.0 / 3.0) + self.k_steady * t
    }
    /// Instantaneous creep rate.
    ///
    /// dε/dt = β/(3*t^(2/3)) + k
    pub fn strain_rate(&self, t: f64) -> f64 {
        if t < 1e-30 {
            return f64::INFINITY;
        }
        self.beta / (3.0 * t.powf(2.0 / 3.0)) + self.k_steady
    }
    /// Time at which the primary and secondary contributions are equal.
    ///
    /// β * t^(1/3) = k * t  →  t_transition = (β/k)^(3/2)
    pub fn transition_time(&self) -> f64 {
        if self.k_steady <= 0.0 {
            return f64::INFINITY;
        }
        (self.beta / self.k_steady).powf(1.5)
    }
}
/// Manson-Haferd time-temperature parameter.
///
/// P_MH = (log10(t_r) - log10(t_a)) / (T - T_a)
///
/// where t_a and T_a are material constants (coordinates of the convergence point).
pub struct MansonHaferdParameter {
    /// Log(ta) — log of the convergence time.
    pub log10_ta: f64,
    /// T_a — convergence temperature (K).
    pub t_a: f64,
}
impl MansonHaferdParameter {
    /// Create a new Manson-Haferd parameter.
    pub fn new(log10_ta: f64, t_a: f64) -> Self {
        Self { log10_ta, t_a }
    }
    /// Compute the MH parameter for given rupture time (hours) and temperature (K).
    pub fn mh_parameter(&self, time_hours: f64, temperature: f64) -> f64 {
        (time_hours.log10() - self.log10_ta) / (temperature - self.t_a)
    }
    /// Compute rupture time (hours) given temperature and parameter P.
    pub fn rupture_time(&self, temperature: f64, p_mh: f64) -> f64 {
        let log10_t = self.log10_ta + p_mh * (temperature - self.t_a);
        10.0_f64.powf(log10_t)
    }
    /// Maximum allowable temperature for a given life and parameter.
    pub fn max_temperature(&self, time_hours: f64, p_mh: f64) -> f64 {
        let log10_t = time_hours.log10();
        self.t_a + (log10_t - self.log10_ta) / p_mh
    }
}
/// Sherby-Dorn parameter for creep life prediction under variable temperature.
///
/// P_SD = ∫ exp(-Q / (R * T)) dt ≈ Σ exp(-Q/(RT_i)) * Δt_i
///
/// For constant conditions: P_SD = t * exp(-Q/(R*T))
pub struct SherbyDornParameter {
    /// Activation energy Q (J/mol).
    pub activation_energy: f64,
}
impl SherbyDornParameter {
    /// Create a new Sherby-Dorn parameter.
    pub fn new(activation_energy: f64) -> Self {
        Self { activation_energy }
    }
    /// Compute P_SD for constant temperature and time.
    ///
    /// P_SD = t_hours * exp(-Q/(R*T))
    pub fn parameter_constant_temp(&self, time_hours: f64, temperature: f64) -> f64 {
        let arrhenius = (-self.activation_energy / (GAS_CONSTANT * temperature)).exp();
        time_hours * arrhenius
    }
    /// Compute temperature-compensated time from variable temperature history.
    ///
    /// P_SD = Σ exp(-Q/(R*T_i)) * Δt_i
    ///
    /// Accepts a slice of (delta_t, temperature) pairs.
    pub fn parameter_variable_temp(&self, segments: &[(f64, f64)]) -> f64 {
        segments
            .iter()
            .map(|&(dt, temp)| {
                let arr = (-self.activation_energy / (GAS_CONSTANT * temp)).exp();
                dt * arr
            })
            .sum()
    }
    /// Equivalent rupture time at a reference temperature that gives the same P_SD.
    ///
    /// t_ref = P_SD / exp(-Q/(R * T_ref))
    pub fn equivalent_time(&self, p_sd: f64, temperature_ref: f64) -> f64 {
        let arr = (-self.activation_energy / (GAS_CONSTANT * temperature_ref)).exp();
        if arr <= 0.0 {
            return f64::INFINITY;
        }
        p_sd / arr
    }
}
/// Coble (grain-boundary diffusion) creep model.
///
/// eps_dot = A_C * D_gb * delta * sigma * Omega / (k_B * T * d^3)
///
/// where D_gb is grain-boundary diffusion coefficient, delta is GB width.
pub struct CobleCreep {
    /// GB diffusion pre-exponential D_gb0 (m^2/s).
    pub d_gb0: f64,
    /// Activation energy for GB diffusion (J/mol).
    pub activation_energy: f64,
    /// Grain boundary width delta (m).
    pub gb_width: f64,
    /// Atomic volume Omega (m^3).
    pub atomic_volume: f64,
    /// Grain size d (m).
    pub grain_size: f64,
    /// Numerical constant A_C (typically ~50).
    pub a_c: f64,
}
impl CobleCreep {
    /// Create a new Coble creep model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        d_gb0: f64,
        activation_energy: f64,
        gb_width: f64,
        atomic_volume: f64,
        grain_size: f64,
    ) -> Self {
        Self {
            d_gb0,
            activation_energy,
            gb_width,
            atomic_volume,
            grain_size,
            a_c: 50.0,
        }
    }
    /// Grain-boundary diffusion coefficient at temperature T.
    pub fn gb_diffusion_coefficient(&self, temperature: f64) -> f64 {
        self.d_gb0 * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp()
    }
    /// Creep strain rate.
    pub fn strain_rate(&self, stress: f64, temperature: f64) -> f64 {
        let k_b = 1.380_649e-23;
        let d_gb = self.gb_diffusion_coefficient(temperature);
        self.a_c * d_gb * self.gb_width * stress * self.atomic_volume
            / (k_b * temperature * self.grain_size.powi(3))
    }
}
impl CobleCreep {
    /// Effective grain-boundary diffusivity at temperature.
    pub fn gb_diffusivity(&self, temperature: f64) -> f64 {
        self.d_gb0 * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp()
    }
    /// Relative contribution of Coble vs Nabarro-Herring at given conditions.
    ///
    /// Returns Coble fraction in \[0, 1\].
    pub fn coble_fraction(&self, nh: &NabarroHerringCreep, stress: f64, temperature: f64) -> f64 {
        let rate_coble = self.strain_rate(stress, temperature);
        let rate_nh = nh.strain_rate(stress, temperature);
        let total = rate_coble + rate_nh;
        if total > 0.0 { rate_coble / total } else { 0.5 }
    }
}
/// Nabarro-Herring (lattice diffusion) creep model.
///
/// eps_dot = A_NH * D_v * sigma * Omega / (k_B * T * d^2)
///
/// where D_v is lattice diffusion coefficient, Omega is atomic volume,
/// d is grain size.
pub struct NabarroHerringCreep {
    /// Diffusion pre-exponential D_0 (m^2/s).
    pub d0: f64,
    /// Activation energy for lattice diffusion (J/mol).
    pub activation_energy: f64,
    /// Atomic volume Omega (m^3).
    pub atomic_volume: f64,
    /// Grain size d (m).
    pub grain_size: f64,
    /// Numerical constant A_NH (typically ~14).
    pub a_nh: f64,
}
impl NabarroHerringCreep {
    /// Create a new Nabarro-Herring creep model.
    pub fn new(d0: f64, activation_energy: f64, atomic_volume: f64, grain_size: f64) -> Self {
        Self {
            d0,
            activation_energy,
            atomic_volume,
            grain_size,
            a_nh: 14.0,
        }
    }
    /// Lattice diffusion coefficient at temperature T.
    pub fn diffusion_coefficient(&self, temperature: f64) -> f64 {
        self.d0 * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp()
    }
    /// Creep strain rate at given stress and temperature.
    pub fn strain_rate(&self, stress: f64, temperature: f64) -> f64 {
        let k_b = 1.380_649e-23;
        let d_v = self.diffusion_coefficient(temperature);
        self.a_nh * d_v * stress * self.atomic_volume
            / (k_b * temperature * self.grain_size * self.grain_size)
    }
}
impl NabarroHerringCreep {
    /// Activation energy Q (J/mol) for diffusion.
    pub fn activation_energy(&self) -> f64 {
        self.activation_energy
    }
    /// Effective diffusivity at temperature (D_eff = D_0 * exp(-Q/RT)).
    pub fn effective_diffusivity(&self, temperature: f64) -> f64 {
        self.d0 * (-self.activation_energy / (GAS_CONSTANT * temperature)).exp()
    }
    /// Strain rate normalised by stress (creep compliance, s^-1 Pa^-1).
    pub fn compliance(&self, temperature: f64) -> f64 {
        let d_eff = self.effective_diffusivity(temperature);
        self.a_nh * d_eff * self.atomic_volume / (self.grain_size * self.grain_size)
    }
}
/// Tertiary creep model using Kachanov-Rabotnov damage mechanics.
///
/// Strain rate: ε̇ = A * σ^n / (1 − ω)^n
/// Damage rate: ω̇ = B * σ^χ / (1 − ω)^φ
#[allow(dead_code)]
pub struct TertiaryCreepModel {
    /// Creep constant A
    pub a: f64,
    /// Creep stress exponent n
    pub n: f64,
    /// Damage constant B
    pub b_damage: f64,
    /// Damage stress exponent χ
    pub chi: f64,
    /// Damage softening exponent φ
    pub phi: f64,
    /// Current damage variable ω ∈ \[0, 1)
    pub omega: f64,
}
impl TertiaryCreepModel {
    /// Create a new tertiary creep model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(a: f64, n: f64, b_damage: f64, chi: f64, phi: f64) -> Self {
        Self {
            a,
            n,
            b_damage,
            chi,
            phi,
            omega: 0.0,
        }
    }
    /// Creep strain rate at current damage state.
    pub fn strain_rate(&self, stress: f64) -> f64 {
        let denom = (1.0 - self.omega).max(1e-9);
        self.a * stress.powf(self.n) / denom.powf(self.n)
    }
    /// Damage rate dω/dt at current damage state.
    pub fn damage_rate(&self, stress: f64) -> f64 {
        let denom = (1.0 - self.omega).max(1e-9);
        self.b_damage * stress.powf(self.chi) / denom.powf(self.phi)
    }
    /// Advance one time step using forward Euler.
    pub fn step(&mut self, stress: f64, dt: f64) {
        let eps_dot = self.strain_rate(stress);
        let omega_dot = self.damage_rate(stress);
        let _ = eps_dot;
        self.omega = (self.omega + omega_dot * dt).min(0.9999);
    }
    /// True when tertiary stage is complete (ω → 1 numerically).
    pub fn is_ruptured(&self) -> bool {
        self.omega >= 0.999
    }
    /// Estimated time to rupture (approximate, constant stress).
    pub fn rupture_time(&self, stress: f64) -> f64 {
        let omega_dot = self.damage_rate(stress);
        if omega_dot > 0.0 {
            (1.0 - self.omega) / omega_dot
        } else {
            f64::INFINITY
        }
    }
    /// Reset damage to zero.
    pub fn reset(&mut self) {
        self.omega = 0.0;
    }
}
/// Kachanov-Rabotnov continuum damage mechanics model.
pub struct CreepDamage {
    /// Damage rate constant A
    pub a: f64,
    /// Stress exponent m
    pub m: f64,
    /// Damage exponent r
    pub r: f64,
}
impl CreepDamage {
    /// Create a new CreepDamage model.
    pub fn new(a: f64, m: f64, r: f64) -> Self {
        Self { a, m, r }
    }
    /// Damage evolution rate.
    ///
    /// dot_omega = A * sigma^m / (1 - omega)^r
    pub fn damage_rate(&self, stress: f64, damage: f64) -> f64 {
        let denominator = (1.0 - damage).powf(self.r);
        if denominator < f64::EPSILON {
            f64::INFINITY
        } else {
            self.a * stress.powf(self.m) / denominator
        }
    }
    /// Integrate damage evolution from 0.
    pub fn integrate(&self, stress: f64, dt: f64, n_steps: usize) -> Vec<f64> {
        let mut omega = 0.0_f64;
        let mut history = Vec::with_capacity(n_steps);
        for _ in 0..n_steps {
            let rate = self.damage_rate(stress, omega);
            omega = (omega + rate * dt).min(1.0);
            history.push(omega);
            if omega >= 0.99 {
                break;
            }
        }
        history
    }
    /// Analytical time to failure.
    ///
    /// t_f = 1 / ((r+1) * A * sigma^m)
    pub fn time_to_failure(&self, stress: f64) -> f64 {
        let base_rate = self.a * stress.powf(self.m);
        if base_rate <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / ((self.r + 1.0) * base_rate)
    }
    /// Effective stress accounting for damage: sigma_eff = sigma / (1 - omega).
    pub fn effective_stress(&self, stress: f64, damage: f64) -> f64 {
        if damage >= 1.0 {
            f64::INFINITY
        } else {
            stress / (1.0 - damage)
        }
    }
}
/// Stress-dependent Larson-Miller parameter using a power-law fit.
///
/// The LM parameter P is a function of stress:
///   P(σ) = A - B * log10(σ/σ_ref)
///
/// where A and B are material constants fit to rupture data.
pub struct LarsonMillerStress {
    /// LM parameter at reference stress: P(σ_ref) = A.
    pub a_lm: f64,
    /// Slope dP/d(log σ): B (negative for decreasing P with stress).
    pub b_lm: f64,
    /// Reference stress (Pa).
    pub sigma_ref: f64,
    /// Larson-Miller constant C.
    pub c_lm: f64,
}
impl LarsonMillerStress {
    /// Create a new stress-dependent LM model.
    pub fn new(a_lm: f64, b_lm: f64, sigma_ref: f64, c_lm: f64) -> Self {
        Self {
            a_lm,
            b_lm,
            sigma_ref,
            c_lm,
        }
    }
    /// Compute the LM parameter for a given stress.
    pub fn lm_parameter(&self, sigma: f64) -> f64 {
        if sigma <= 0.0 {
            return f64::INFINITY;
        }
        self.a_lm - self.b_lm * (sigma / self.sigma_ref).log10()
    }
    /// Predict rupture life (hours) for a given stress and temperature.
    ///
    /// t_r = 10^(P(σ)/T - C)
    pub fn rupture_time(&self, sigma: f64, temperature: f64) -> f64 {
        let p = self.lm_parameter(sigma);
        let log10_t = p / temperature - self.c_lm;
        10.0_f64.powf(log10_t)
    }
    /// Maximum allowable stress for a given life and temperature.
    ///
    /// Inverts P(σ) = T*(C + log10(t)) to give σ.
    pub fn max_stress(&self, temperature: f64, time_hours: f64) -> f64 {
        let p_required = temperature * (self.c_lm + time_hours.log10());
        if self.b_lm.abs() < 1e-30 {
            return self.sigma_ref;
        }
        let log_sigma_ratio = (self.a_lm - p_required) / self.b_lm;
        self.sigma_ref * 10.0_f64.powf(log_sigma_ratio)
    }
}
/// Creep rupture envelope using Larson-Miller master curve.
///
/// Maps (temperature, stress) → predicted rupture life using a polynomial
/// master curve fit to the Larson-Miller parameter.
#[derive(Debug, Clone)]
pub struct RuptureEnvelope {
    /// Larson-Miller constant C.
    pub c: f64,
    /// Polynomial coefficients \[a0, a1, a2\] for LM master curve:
    /// P_LM = a0 + a1*log10(sigma) + a2*log10(sigma)^2
    pub curve_coeffs: [f64; 3],
}
impl RuptureEnvelope {
    /// Create a rupture envelope.
    pub fn new(c: f64, curve_coeffs: [f64; 3]) -> Self {
        Self { c, curve_coeffs }
    }
    /// Evaluate Larson-Miller parameter from stress (Pa).
    pub fn lm_parameter_from_stress(&self, stress: f64) -> f64 {
        if stress <= 0.0 {
            return 0.0;
        }
        let log_s = stress.log10();
        self.curve_coeffs[0] + self.curve_coeffs[1] * log_s + self.curve_coeffs[2] * log_s * log_s
    }
    /// Rupture life (s) at given stress and temperature (K).
    pub fn rupture_life(&self, stress: f64, temperature: f64) -> f64 {
        let p = self.lm_parameter_from_stress(stress);
        if temperature <= 0.0 {
            return f64::INFINITY;
        }
        let exponent = p / temperature - self.c;
        10.0_f64.powf(exponent)
    }
    /// Maximum allowable stress for a given temperature and design life (s).
    ///
    /// Finds σ such that lm_parameter_from_stress(σ) = T*(C + log10(t_r)).
    /// Assumes the LM master curve is monotonically decreasing with stress
    /// (negative first coefficient).
    pub fn allowable_stress(&self, temperature: f64, life_s: f64) -> f64 {
        if life_s <= 0.0 || temperature <= 0.0 {
            return 0.0;
        }
        let p_target = temperature * (self.c + life_s.log10());
        let mut lo = 1.0e3_f64;
        let mut hi = 1.0e9_f64;
        let p_lo = self.lm_parameter_from_stress(lo);
        let p_hi = self.lm_parameter_from_stress(hi);
        if p_lo <= p_hi {
            for _ in 0..80 {
                let mid = (lo + hi) / 2.0;
                if self.lm_parameter_from_stress(mid) < p_target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        } else {
            for _ in 0..80 {
                let mid = (lo + hi) / 2.0;
                if self.lm_parameter_from_stress(mid) > p_target {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
        }
        (lo + hi) / 2.0
    }
}
/// Composite three-stage creep model combining primary, secondary, and
/// tertiary contributions.
#[allow(dead_code)]
pub struct ThreeStageCreep {
    /// Primary creep parameters.
    pub primary: PrimaryCreepModel,
    /// Secondary (steady-state) creep.
    pub secondary: SecondaryCreepModel,
    /// Tertiary damage model.
    pub tertiary: TertiaryCreepModel,
    /// Strain threshold beyond which tertiary stage begins.
    pub tertiary_onset_strain: f64,
    /// Accumulated total strain.
    pub strain: f64,
    /// Elapsed time.
    pub time: f64,
}
impl ThreeStageCreep {
    /// Create a three-stage creep model.
    pub fn new(
        primary: PrimaryCreepModel,
        secondary: SecondaryCreepModel,
        tertiary: TertiaryCreepModel,
        tertiary_onset_strain: f64,
    ) -> Self {
        Self {
            primary,
            secondary,
            tertiary,
            tertiary_onset_strain,
            strain: 0.0,
            time: 0.0,
        }
    }
    /// Advance one time step at given stress and temperature.
    ///
    /// Returns the incremental strain for this step.
    pub fn step(&mut self, stress: f64, temperature: f64, dt: f64) -> f64 {
        let eps_dot_primary = self.primary.strain_rate(stress, self.time + dt);
        let eps_dot_secondary = self.secondary.strain_rate(stress, temperature);
        let eps_dot = if self.strain < self.tertiary_onset_strain {
            eps_dot_primary.max(eps_dot_secondary)
        } else {
            self.tertiary.step(stress, dt);
            self.tertiary.strain_rate(stress)
        };
        let d_eps = eps_dot * dt;
        self.strain += d_eps;
        self.time += dt;
        d_eps
    }
    /// Current creep stage as a string tag.
    pub fn current_stage(&self) -> &'static str {
        if self.tertiary.omega > 0.01 {
            "tertiary"
        } else if self.strain > self.tertiary_onset_strain * 0.5 {
            "secondary"
        } else {
            "primary"
        }
    }
    /// True when material has ruptured.
    pub fn is_ruptured(&self) -> bool {
        self.tertiary.is_ruptured()
    }
}
/// Predict creep rupture life using multiple extrapolation methods.
///
/// Returns a struct with results from different time-temperature parameters.
#[derive(Debug, Clone)]
pub struct RuptureLifePrediction {
    /// Larson-Miller prediction (hours).
    pub larson_miller: f64,
    /// Orr-Sherby-Dorn prediction (hours).
    pub orr_sherby_dorn: f64,
    /// Monkman-Grant prediction (hours).
    pub monkman_grant: f64,
}
/// Integrates Norton creep over a stress-time history at constant temperature.
pub struct ThermalCreepIntegrator {
    /// Norton creep law
    pub norton: NortonCreep,
}
impl ThermalCreepIntegrator {
    /// Create a new ThermalCreepIntegrator.
    pub fn new(norton: NortonCreep) -> Self {
        Self { norton }
    }
    /// Creep strain increment over a time step dt.
    pub fn creep_strain_increment(&self, sigma: f64, temperature: f64, dt: f64) -> f64 {
        self.norton.creep_strain_rate(sigma, temperature) * dt
    }
    /// Integrate stress history (time, stress) pairs at constant temperature.
    pub fn integrate(&self, stress_history: &[(f64, f64)], temperature: f64, dt: f64) -> Vec<f64> {
        let mut accumulated = 0.0_f64;
        let mut strains = Vec::with_capacity(stress_history.len());
        for &(_t, sigma) in stress_history {
            accumulated += self.creep_strain_increment(sigma, temperature, dt);
            strains.push(accumulated);
        }
        strains
    }
}
/// Stress relaxation model (Maxwell model).
///
/// sigma(t) = sigma_0 * exp(-t / tau)
///
/// where tau = eta / E is the relaxation time.
pub struct StressRelaxation {
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Viscosity eta (Pa*s).
    pub viscosity: f64,
}
impl StressRelaxation {
    /// Create a new stress relaxation model.
    pub fn new(young_modulus: f64, viscosity: f64) -> Self {
        Self {
            young_modulus,
            viscosity,
        }
    }
    /// Relaxation time tau = eta / E (s).
    pub fn relaxation_time(&self) -> f64 {
        self.viscosity / self.young_modulus
    }
    /// Stress at time t given initial stress sigma_0.
    pub fn stress_at_time(&self, sigma_0: f64, t: f64) -> f64 {
        let tau = self.relaxation_time();
        sigma_0 * (-t / tau).exp()
    }
    /// Time to relax to a fraction f of the initial stress.
    ///
    /// t = -tau * ln(f)
    pub fn time_to_fraction(&self, fraction: f64) -> f64 {
        let tau = self.relaxation_time();
        -tau * fraction.ln()
    }
    /// Stress history from t=0 to t_end with n_steps.
    pub fn stress_history(&self, sigma_0: f64, t_end: f64, n_steps: usize) -> Vec<(f64, f64)> {
        (0..=n_steps)
            .map(|i| {
                let t = t_end * (i as f64) / (n_steps as f64);
                (t, self.stress_at_time(sigma_0, t))
            })
            .collect()
    }
}
/// Creep-fatigue interaction model using the damage summation approach.
///
/// The interaction rule combines fatigue damage fraction (N/Nf) and creep
/// damage fraction (t/tr) with an interaction factor:
///
///   D_total = Σ(N_i / N_fi) + Σ(t_j / t_rj) + interaction_term
///
/// Failure is predicted when D_total ≥ 1.
///
/// Interaction diagram (Langer-Halford): tri-linear failure boundary.
pub struct CreepFatigueInteraction {
    /// Fatigue intercept on the creep-fatigue diagram (creep damage at zero fatigue).
    pub creep_intercept: f64,
    /// Fatigue intercept (fatigue damage at zero creep).
    pub fatigue_intercept: f64,
    /// Interaction exponent.
    pub interaction_exponent: f64,
}
impl CreepFatigueInteraction {
    /// Create a new creep-fatigue interaction model.
    pub fn new(creep_intercept: f64, fatigue_intercept: f64, interaction_exponent: f64) -> Self {
        Self {
            creep_intercept,
            fatigue_intercept,
            interaction_exponent,
        }
    }
    /// Check whether a given (fatigue_damage, creep_damage) combination is safe.
    ///
    /// Uses the Langer interaction locus:
    ///   (Df/Di_f)^n + (Dc/Di_c)^n <= 1
    pub fn is_safe(&self, fatigue_damage: f64, creep_damage: f64) -> bool {
        let n = self.interaction_exponent;
        let term_f = (fatigue_damage / self.fatigue_intercept).powf(n);
        let term_c = (creep_damage / self.creep_intercept).powf(n);
        term_f + term_c <= 1.0
    }
    /// Total interaction damage index.
    pub fn interaction_damage(&self, fatigue_damage: f64, creep_damage: f64) -> f64 {
        let n = self.interaction_exponent;
        let term_f = (fatigue_damage / self.fatigue_intercept).powf(n);
        let term_c = (creep_damage / self.creep_intercept).powf(n);
        (term_f + term_c).powf(1.0 / n)
    }
    /// Remaining allowable creep damage given a fatigue damage fraction.
    pub fn remaining_creep(&self, fatigue_damage: f64) -> f64 {
        let n = self.interaction_exponent;
        let term_f = (fatigue_damage / self.fatigue_intercept).powf(n);
        let remaining = (1.0 - term_f).max(0.0);
        self.creep_intercept * remaining.powf(1.0 / n)
    }
}
/// Multi-axial creep model using the von Mises equivalent stress.
///
/// For a general stress state (σ1, σ2, σ3), the equivalent (von Mises) stress is:
///   σ_eq = sqrt\[ ((σ1-σ2)² + (σ2-σ3)² + (σ3-σ1)²) / 2 \]
///
/// The creep strain rate in each principal direction is:
///   eps_dot_i = eps_dot_eq * (σ_i - (σ1+σ2+σ3)/3) / (σ_eq / (3/2))
///
/// where eps_dot_eq = A * σ_eq^n * exp(-Q/RT).
pub struct MultiaxialCreep {
    /// Norton model for equivalent strain rate.
    pub norton: NortonCreep,
}
impl MultiaxialCreep {
    /// Create a new multiaxial creep model.
    pub fn new(norton: NortonCreep) -> Self {
        Self { norton }
    }
    /// Von Mises equivalent stress.
    pub fn von_mises_stress(&self, s1: f64, s2: f64, s3: f64) -> f64 {
        (0.5 * ((s1 - s2).powi(2) + (s2 - s3).powi(2) + (s3 - s1).powi(2))).sqrt()
    }
    /// Hydrostatic (mean) stress.
    pub fn mean_stress(&self, s1: f64, s2: f64, s3: f64) -> f64 {
        (s1 + s2 + s3) / 3.0
    }
    /// Principal creep strain rates \[eps1_dot, eps2_dot, eps3_dot\].
    ///
    /// Returns `[0, 0, 0]` if von Mises stress is zero.
    pub fn principal_creep_rates(&self, s1: f64, s2: f64, s3: f64, temperature: f64) -> [f64; 3] {
        let s_eq = self.von_mises_stress(s1, s2, s3);
        if s_eq < 1e-30 {
            return [0.0; 3];
        }
        let s_mean = self.mean_stress(s1, s2, s3);
        let eps_dot_eq = self.norton.creep_strain_rate(s_eq, temperature);
        let scale = 1.5 * eps_dot_eq / s_eq;
        [
            scale * (s1 - s_mean),
            scale * (s2 - s_mean),
            scale * (s3 - s_mean),
        ]
    }
    /// Effective volumetric creep strain rate (should be zero for volume-preserving).
    pub fn volumetric_creep_rate(&self, s1: f64, s2: f64, s3: f64, temperature: f64) -> f64 {
        let rates = self.principal_creep_rates(s1, s2, s3, temperature);
        rates[0] + rates[1] + rates[2]
    }
}
/// Simplified deformation mechanism map for a material.
///
/// Based on homologous temperature (T/T_m) and normalised stress (σ/μ),
/// determines the dominant deformation mechanism.
#[derive(Debug, Clone)]
pub struct DeformationMechanismMap {
    /// Melting temperature T_m (K).
    pub t_melting: f64,
    /// Shear modulus μ (Pa).
    pub shear_modulus: f64,
    /// Peierls stress (normalised, σ/μ threshold for dislocation glide).
    pub peierls_stress_norm: f64,
    /// Boundary for power-law creep vs diffusion creep (normalised stress).
    pub power_law_diff_boundary: f64,
}
impl DeformationMechanismMap {
    /// Create a deformation mechanism map.
    pub fn new(
        t_melting: f64,
        shear_modulus: f64,
        peierls_stress_norm: f64,
        power_law_diff_boundary: f64,
    ) -> Self {
        Self {
            t_melting,
            shear_modulus,
            peierls_stress_norm,
            power_law_diff_boundary,
        }
    }
    /// Nickel approximation.
    pub fn nickel() -> Self {
        Self::new(1728.0, 76.0e9, 3e-3, 1e-4)
    }
    /// Identify dominant creep mechanism.
    ///
    /// Returns: "elastic", "dislocation_glide", "power_law_creep", "diffusion_creep", or "superplastic".
    pub fn dominant_mechanism(&self, stress: f64, temperature: f64) -> &'static str {
        let t_hom = temperature / self.t_melting;
        let sigma_norm = stress / self.shear_modulus;
        if sigma_norm > self.peierls_stress_norm {
            "dislocation_glide"
        } else if t_hom < 0.3 {
            "elastic"
        } else if sigma_norm > self.power_law_diff_boundary {
            "power_law_creep"
        } else if t_hom > 0.85 {
            "superplastic"
        } else {
            "diffusion_creep"
        }
    }
    /// Approximate steady-state strain rate combining power law + diffusion creep.
    ///
    /// eps_dot ≈ A_PL * (σ/μ)^5 * exp(-Q_PL/(R*T)) + A_D * (σ/μ) * exp(-Q_D/(R*T))
    pub fn combined_strain_rate(
        &self,
        stress: f64,
        temperature: f64,
        a_pl: f64,
        q_pl: f64,
        a_d: f64,
        q_d: f64,
    ) -> f64 {
        let r = 8.314;
        let sigma_norm = stress / self.shear_modulus;
        let pl = a_pl * sigma_norm.powi(5) * (-q_pl / (r * temperature)).exp();
        let diff = a_d * sigma_norm * (-q_d / (r * temperature)).exp();
        pl + diff
    }
}
/// Norton power-law creep model.
///
/// Strain rate: eps_dot = A * sigma^n * exp(-Q / (R * T))
pub struct NortonCreep {
    /// Norton constant A (s^-1 Pa^-n)
    pub a: f64,
    /// Stress exponent n
    pub n: f64,
    /// Activation energy Q (J/mol)
    pub activation_energy: f64,
    /// Universal gas constant R (J/mol/K)
    pub gas_constant: f64,
}
impl NortonCreep {
    /// Create a new NortonCreep model.
    pub fn new(a: f64, n: f64, activation_energy: f64) -> Self {
        Self {
            a,
            n,
            activation_energy,
            gas_constant: GAS_CONSTANT,
        }
    }
    /// Compute steady-state creep strain rate at given stress and temperature.
    pub fn creep_strain_rate(&self, stress: f64, temperature: f64) -> f64 {
        let arrhenius = (-self.activation_energy / (self.gas_constant * temperature)).exp();
        self.a * stress.powf(self.n) * arrhenius
    }
    /// Normalized Norton-Bailey creep rate using yield stress.
    pub fn creep_strain_rate_normalized(
        &self,
        stress: f64,
        yield_stress: f64,
        temperature: f64,
    ) -> f64 {
        let sigma_norm = stress / yield_stress;
        let arrhenius = (-self.activation_energy / (self.gas_constant * temperature)).exp();
        self.a * sigma_norm.powf(self.n) * arrhenius
    }
    /// Estimate time to rupture (simplified).
    pub fn time_to_rupture(&self, stress: f64, temperature: f64) -> f64 {
        let rate = self.creep_strain_rate(stress, temperature);
        if rate > 0.0 {
            1.0 / rate
        } else {
            f64::INFINITY
        }
    }
    /// Alias for [`creep_strain_rate`](Self::creep_strain_rate).
    pub fn creep_rate(&self, stress: f64, temperature: f64) -> f64 {
        self.creep_strain_rate(stress, temperature)
    }

    /// Effective viscosity at given stress and temperature.
    ///
    /// eta = sigma / (3 * eps_dot)
    pub fn effective_viscosity(&self, stress: f64, temperature: f64) -> f64 {
        let rate = self.creep_strain_rate(stress, temperature);
        if rate > 1e-30 {
            stress / (3.0 * rate)
        } else {
            f64::INFINITY
        }
    }
}
/// Creep stage classification.
pub enum CreepStage {
    /// Primary (transient) creep: eps ~ t^m with m < 1.
    Primary {
        /// Time exponent (typically 0.3-0.5).
        exponent: f64,
    },
    /// Secondary (steady-state) creep: constant strain rate.
    Secondary,
    /// Tertiary creep: accelerating rate leading to fracture.
    Tertiary {
        /// Acceleration coefficient.
        acceleration: f64,
    },
}
/// Orr-Sherby-Dorn (OSD) time-temperature parameter.
///
/// P_OSD = log10(t_r) - Q_OSD / (2.303 * R * T)
///
/// where Q_OSD is the creep activation energy.
pub struct OrrSherbyDornParameter {
    /// Activation energy Q (J/mol).
    pub q: f64,
}
impl OrrSherbyDornParameter {
    /// Create a new OSD parameter.
    pub fn new(q: f64) -> Self {
        Self { q }
    }
    /// Compute OSD parameter for given rupture time (hours) and temperature (K).
    pub fn osd_parameter(&self, time_hours: f64, temperature: f64) -> f64 {
        time_hours.log10() - self.q / (2.303 * GAS_CONSTANT * temperature)
    }
    /// Compute rupture time (hours) given temperature and parameter.
    pub fn rupture_time(&self, temperature: f64, p_osd: f64) -> f64 {
        let log10_t = p_osd + self.q / (2.303 * GAS_CONSTANT * temperature);
        10.0_f64.powf(log10_t)
    }
}
/// Coupled creep-damage model: damage accelerates the creep rate.
///
/// Effective creep strain rate accounting for damage:
///   eps_dot = A * (sigma / (1 - omega))^n * exp(-Q / RT)
///
/// Damage evolution:
///   omega_dot = B * sigma^m / (1 - omega)^phi
pub struct CoupledCreepDamage {
    /// Norton creep constant A (s^-1 Pa^-n).
    pub a: f64,
    /// Creep stress exponent n.
    pub n: f64,
    /// Activation energy Q (J/mol).
    pub q: f64,
    /// Damage rate constant B.
    pub b_damage: f64,
    /// Damage stress exponent m.
    pub m_damage: f64,
    /// Damage evolution exponent phi.
    pub phi_damage: f64,
}
impl CoupledCreepDamage {
    /// Create a new coupled creep-damage model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(a: f64, n: f64, q: f64, b_damage: f64, m_damage: f64, phi_damage: f64) -> Self {
        Self {
            a,
            n,
            q,
            b_damage,
            m_damage,
            phi_damage,
        }
    }
    /// Effective creep strain rate at given stress, temperature and damage.
    pub fn creep_rate(&self, stress: f64, temperature: f64, damage: f64) -> f64 {
        let eff_stress = stress / (1.0 - damage).max(f64::EPSILON);
        let arr = (-self.q / (GAS_CONSTANT * temperature)).exp();
        self.a * eff_stress.powf(self.n) * arr
    }
    /// Damage evolution rate.
    pub fn damage_rate(&self, stress: f64, damage: f64) -> f64 {
        let denom = (1.0 - damage).powf(self.phi_damage);
        if denom < f64::EPSILON {
            f64::INFINITY
        } else {
            self.b_damage * stress.powf(self.m_damage) / denom
        }
    }
    /// Integrate coupled creep-damage equations using explicit Euler.
    ///
    /// Returns (strain history, damage history) each of length n_steps+1.
    pub fn integrate(
        &self,
        stress: f64,
        temperature: f64,
        dt: f64,
        n_steps: usize,
    ) -> (Vec<f64>, Vec<f64>) {
        let mut eps = 0.0_f64;
        let mut omega = 0.0_f64;
        let mut eps_hist = Vec::with_capacity(n_steps + 1);
        let mut dam_hist = Vec::with_capacity(n_steps + 1);
        eps_hist.push(eps);
        dam_hist.push(omega);
        for _ in 0..n_steps {
            let cr = self.creep_rate(stress, temperature, omega);
            let dr = self.damage_rate(stress, omega);
            eps += cr * dt;
            omega = (omega + dr * dt).min(1.0);
            eps_hist.push(eps);
            dam_hist.push(omega);
            if omega >= 0.99 {
                break;
            }
        }
        (eps_hist, dam_hist)
    }
}
/// Secondary (steady-state) creep stage — delegates to NortonCreep.
///
/// This is a thin wrapper that represents the constant-rate stage.
#[allow(dead_code)]
pub struct SecondaryCreepModel {
    /// Underlying Norton power-law parameters.
    pub norton: NortonCreep,
}
impl SecondaryCreepModel {
    /// Create from Norton parameters.
    pub fn new(a: f64, n: f64, activation_energy: f64) -> Self {
        Self {
            norton: NortonCreep::new(a, n, activation_energy),
        }
    }
    /// Steady-state creep rate (same as NortonCreep).
    pub fn strain_rate(&self, stress: f64, temperature: f64) -> f64 {
        self.norton.creep_strain_rate(stress, temperature)
    }
    /// Accumulated strain after duration `dt` at steady state.
    pub fn accumulated_strain(&self, stress: f64, temperature: f64, dt: f64) -> f64 {
        self.strain_rate(stress, temperature) * dt
    }
    /// Effective viscosity η = σ / (3 ε̇).
    pub fn viscosity(&self, stress: f64, temperature: f64) -> f64 {
        self.norton.effective_viscosity(stress, temperature)
    }
}
/// Perzyna viscoplastic model.
pub struct ViscoplasticModel {
    /// Initial yield stress (Pa)
    pub yield_stress: f64,
    /// Viscosity parameter eta (Pa*s)
    pub viscosity: f64,
    /// Isotropic hardening modulus H (Pa)
    pub hardening: f64,
}
impl ViscoplasticModel {
    /// Create a new ViscoplasticModel.
    pub fn new(yield_stress: f64, viscosity: f64, hardening: f64) -> Self {
        Self {
            yield_stress,
            viscosity,
            hardening,
        }
    }
    /// Overstress function (Macaulay bracket).
    pub fn overstress(&self, total_stress: f64, accumulated_plastic: f64) -> f64 {
        let flow_stress = self.yield_stress + self.hardening * accumulated_plastic;
        let excess = total_stress.abs() - flow_stress;
        if excess > 0.0 {
            excess / self.yield_stress
        } else {
            0.0
        }
    }
    /// Plastic strain rate: eps_dot_p = (1/eta) * phi.
    pub fn plastic_strain_rate(&self, stress: f64, plastic_strain: f64) -> f64 {
        let phi = self.overstress(stress, plastic_strain);
        phi / self.viscosity
    }
    /// Explicit Euler step: returns (new_plastic_strain, new_accumulated_plastic).
    pub fn step(&self, stress: f64, plastic_strain: f64, dt: f64) -> (f64, f64) {
        let rate = self.plastic_strain_rate(stress, plastic_strain);
        let dep = rate * dt * stress.signum();
        let new_plastic = plastic_strain + dep;
        let new_accumulated = plastic_strain.abs() + dep.abs();
        (new_plastic, new_accumulated)
    }
}
/// Larson-Miller parameter for creep rupture life prediction.
pub struct LarsonMillerParameter {
    /// Larson-Miller constant C (typically ~20 for metals)
    pub c_lm: f64,
}
impl LarsonMillerParameter {
    /// Create a new LarsonMillerParameter.
    pub fn new(c_lm: f64) -> Self {
        Self { c_lm }
    }
    /// Compute Larson-Miller parameter P = T * (C + log10(t_r)).
    pub fn lm_parameter(&self, temperature: f64, time_hours: f64) -> f64 {
        temperature * (self.c_lm + time_hours.log10())
    }
    /// Compute rupture time in hours given temperature (K) and LM parameter.
    pub fn rupture_time(&self, temperature: f64, lm_param: f64) -> f64 {
        let log10_t = lm_param / temperature - self.c_lm;
        10_f64.powf(log10_t)
    }
    /// Maximum allowable temperature (K) for a given life and LM parameter.
    pub fn max_temperature(&self, time_hours: f64, lm_param: f64) -> f64 {
        lm_param / (self.c_lm + time_hours.log10())
    }
    /// Equivalent time at a different temperature for the same LM parameter.
    ///
    /// Given a known (T1, t1), compute t2 at temperature T2.
    pub fn equivalent_time(&self, t1: f64, temp1: f64, temp2: f64) -> f64 {
        let p = self.lm_parameter(temp1, t1);
        self.rupture_time(temp2, p)
    }
}
/// Zener-Hollomon parameter Z = ε_dot * exp(Q/(R*T)).
///
/// Combines the effect of strain rate and temperature into a single
/// temperature-compensated strain rate. Useful for hot deformation.
pub struct ZenerHollomon {
    /// Activation energy Q (J/mol).
    pub activation_energy: f64,
}
impl ZenerHollomon {
    /// Create a new Zener-Hollomon parameter.
    pub fn new(activation_energy: f64) -> Self {
        Self { activation_energy }
    }
    /// Compute the Zener-Hollomon parameter Z.
    ///
    /// Z = ε_dot * exp(Q / (R * T))
    pub fn z_parameter(&self, strain_rate: f64, temperature: f64) -> f64 {
        let arr = (self.activation_energy / (GAS_CONSTANT * temperature)).exp();
        strain_rate * arr
    }
    /// Compute the equivalent strain rate at a different temperature
    /// that gives the same Z value.
    ///
    /// ε_dot_2 = Z * exp(-Q/(R*T_2))
    pub fn equivalent_strain_rate(&self, z: f64, temperature: f64) -> f64 {
        let arr = (-self.activation_energy / (GAS_CONSTANT * temperature)).exp();
        z * arr
    }
    /// Compute the temperature required to achieve a given Z at a given strain rate.
    ///
    /// T = Q / (R * ln(Z / ε_dot))
    pub fn temperature_for_z(&self, z: f64, strain_rate: f64) -> f64 {
        if strain_rate <= 0.0 || z <= 0.0 {
            return f64::NAN;
        }
        self.activation_energy / (GAS_CONSTANT * (z / strain_rate).ln())
    }
}
/// Norton-Bailey time-hardening creep model.
///
/// ε(t) = A * σ^n * t^m
///
/// where m is the time hardening exponent (0 < m < 1) and A, n are material
/// constants. The instantaneous strain rate is:
///
/// dε/dt = A * m * σ^n * t^(m-1)
pub struct NortonBaileyTimeHardening {
    /// Creep coefficient A.
    pub a: f64,
    /// Stress exponent n.
    pub n: f64,
    /// Time exponent m (0 < m ≤ 1).
    pub m: f64,
}
impl NortonBaileyTimeHardening {
    /// Create a new Norton-Bailey time-hardening model.
    pub fn new(a: f64, n: f64, m: f64) -> Self {
        Self { a, n, m }
    }
    /// Creep strain at time t.
    pub fn strain(&self, sigma: f64, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        self.a * sigma.powf(self.n) * t.powf(self.m)
    }
    /// Creep strain rate at time t.
    pub fn strain_rate(&self, sigma: f64, t: f64) -> f64 {
        if t < 1e-30 {
            return f64::INFINITY;
        }
        self.a * self.m * sigma.powf(self.n) * t.powf(self.m - 1.0)
    }
    /// Equivalent time method: given accumulated strain ε_0 at stress σ_0,
    /// find the time t_eq on the new-stress creep curve.
    ///
    /// t_eq = (ε_0 / (A * σ_new^n))^(1/m)
    pub fn equivalent_time(&self, accumulated_strain: f64, sigma_new: f64) -> f64 {
        if accumulated_strain <= 0.0 || sigma_new <= 0.0 {
            return 0.0;
        }
        let base = accumulated_strain / (self.a * sigma_new.powf(self.n));
        base.powf(1.0 / self.m)
    }
}
