//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::E;

use super::functions::R_GAS;
#[allow(unused_imports)]
use super::functions::*;

/// Stress-assisted and strain-induced martensitic transformation model.
///
/// The martensite start temperature shifts with applied stress:
/// `Ms(σ) = Ms_0 - k_σ * σ` for stress-assisted transformation.
///
/// For strain-induced transformation, uses the Olson-Cohen model:
/// `f_m = 1 - exp(-β * f_sh^n)`
/// where `f_sh` is the volume fraction of shear bands.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct StressAidedMartensite {
    /// Stress-free martensite start temperature Ms₀ \[K\]
    pub ms_0: f64,
    /// Stress coefficient dMs/dσ \[K/Pa\] (typically negative)
    pub k_sigma: f64,
    /// Olson-Cohen β parameter (probability of martensite formation per shear band intersection)
    pub beta_oc: f64,
    /// Olson-Cohen n exponent (geometrical factor, typically ~4.5)
    pub n_oc: f64,
    /// Shear band volume fraction parameter α_oc (controls shear band density)
    pub alpha_oc: f64,
}
impl StressAidedMartensite {
    /// Create a new stress-aided martensite model.
    pub fn new(ms_0: f64, k_sigma: f64, beta_oc: f64, n_oc: f64, alpha_oc: f64) -> Self {
        Self {
            ms_0,
            k_sigma,
            beta_oc,
            n_oc,
            alpha_oc,
        }
    }
    /// Stress-shifted Ms temperature \[K\].
    ///
    /// `Ms(σ) = Ms₀ - k_σ * σ`
    pub fn ms_temperature_under_stress(&self, stress_pa: f64) -> f64 {
        self.ms_0 - self.k_sigma * stress_pa
    }
    /// Shear band volume fraction at equivalent plastic strain ε_p.
    ///
    /// `f_sh = 1 - exp(-α_oc * ε_p)`
    pub fn shear_band_fraction(&self, strain_p: f64) -> f64 {
        if strain_p <= 0.0 {
            return 0.0;
        }
        (1.0 - (-self.alpha_oc * strain_p).exp()).clamp(0.0, 1.0)
    }
    /// Olson-Cohen martensite fraction from plastic strain.
    ///
    /// `f_m = 1 - exp(-β * f_sh^n)`
    pub fn martensite_fraction_from_strain(&self, strain_p: f64) -> f64 {
        let f_sh = self.shear_band_fraction(strain_p);
        (1.0 - (-self.beta_oc * f_sh.powf(self.n_oc)).exp()).clamp(0.0, 1.0)
    }
    /// Whether transformation is stress-assisted (T > Ms(σ)) or thermal (T < Ms₀).
    pub fn transformation_mode(&self, temperature: f64, stress_pa: f64) -> &'static str {
        let ms_stress = self.ms_temperature_under_stress(stress_pa);
        if temperature < self.ms_0 {
            "thermal"
        } else if temperature < ms_stress {
            "stress-assisted"
        } else {
            "none"
        }
    }
}
/// Multi-phase Allen-Cahn order parameter field.
///
/// Allen-Cahn equation:
///   dφ_i/dt = -M * δF/δφ_i = M * (-dF/dφ_i)
#[derive(Debug, Clone)]
pub struct PhaseFieldOrder {
    /// Order parameter values for each grid point (length = n_phases).
    pub phi: Vec<f64>,
    /// Number of phases tracked.
    pub n_phases: usize,
}
impl PhaseFieldOrder {
    /// Create a new phase field with all order parameters initialised to zero.
    pub fn new(n_phases: usize) -> Self {
        Self {
            phi: vec![0.0; n_phases],
            n_phases,
        }
    }
    /// Advance each φ_i by one Allen-Cahn step.
    ///
    /// dφ_i/dt = M * (-dF/dφ_i)
    ///
    /// # Arguments
    /// * `dt`      - Time step (s).
    /// * `m`       - Mobility coefficient M.
    /// * `df_dphi` - Variational derivative δF/δφ for each phase (length = n_phases).
    #[allow(clippy::needless_range_loop)]
    pub fn update_allen_cahn(&mut self, dt: f64, m: f64, df_dphi: &[f64]) {
        assert_eq!(df_dphi.len(), self.n_phases, "df_dphi length mismatch");
        for i in 0..self.n_phases {
            self.phi[i] += dt * m * (-df_dphi[i]);
        }
    }
    /// Total "free energy" proxy: Σ φ_i^2 (for monitoring purposes).
    #[allow(dead_code)]
    pub fn order_norm(&self) -> f64 {
        self.phi.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}
/// Bainite transformation model using JMAK (Avrami) kinetics.
///
/// The fraction transformed follows:
/// X_b(t, T) = 1 - exp(-k(T) · t^n)
///
/// where k(T) = k0 · exp(-Q_b / (R·T)) is Arrhenius-activated.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BainiteTransformation {
    /// Avrami exponent n (dimensionless, typically 2..4 for bainite).
    pub n: f64,
    /// Pre-exponential rate constant k0 \[s^-n\].
    pub k0: f64,
    /// Activation energy for bainite transformation Q_b \[J/mol\].
    pub q_b: f64,
    /// Bainite start temperature B_s \[K\].
    pub b_s: f64,
    /// Bainite finish temperature B_f \[K\].
    pub b_f: f64,
    /// Current bainite volume fraction.
    pub fraction: f64,
}
impl BainiteTransformation {
    /// Create a new bainite transformation model.
    pub fn new(n: f64, k0: f64, q_b: f64, b_s: f64, b_f: f64) -> Self {
        Self {
            n,
            k0,
            q_b,
            b_s,
            b_f,
            fraction: 0.0,
        }
    }
    /// Rate constant k at temperature T \[K\].
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        self.k0 * (-self.q_b / (R_GAS * temperature)).exp()
    }
    /// Bainite volume fraction at time t \[s\] and temperature T \[K\].
    ///
    /// Returns 0 if T > B_s (no bainite above start temperature).
    pub fn fraction_at(&self, time: f64, temperature: f64) -> f64 {
        if temperature > self.b_s || time <= 0.0 {
            return 0.0;
        }
        let k = self.rate_constant(temperature);
        1.0 - (-k * time.powf(self.n)).exp()
    }
    /// Time to reach a target fraction at temperature T.
    ///
    /// t = (-ln(1 - X) / k)^(1/n)
    pub fn time_to_fraction(&self, target_fraction: f64, temperature: f64) -> f64 {
        if temperature > self.b_s || target_fraction <= 0.0 {
            return f64::INFINITY;
        }
        let target = target_fraction.clamp(0.0, 1.0 - f64::EPSILON);
        let k = self.rate_constant(temperature);
        if k <= 0.0 {
            return f64::INFINITY;
        }
        ((-1.0 - target).ln().abs() / k).powf(1.0 / self.n)
    }
    /// Nose temperature: temperature giving maximum transformation rate.
    ///
    /// Approximated as the midpoint of the B_s–B_f range.
    pub fn nose_temperature(&self) -> f64 {
        0.5 * (self.b_s + self.b_f)
    }
    /// Isothermal transformation rate dX/dt at time t.
    pub fn transformation_rate(&self, time: f64, temperature: f64) -> f64 {
        if temperature > self.b_s || time <= 0.0 {
            return 0.0;
        }
        let k = self.rate_constant(temperature);
        let x = self.fraction_at(time, temperature);
        k * self.n * time.powf(self.n - 1.0) * (1.0 - x)
    }
}
/// Clausius-Clapeyron equation for solid/liquid phase boundary slope.
///
/// dP/dT = L / (T * ΔV)
///
/// where L is the latent heat, T is the transition temperature, and ΔV is the
/// molar volume change on transformation.
#[derive(Debug, Clone)]
pub struct ClausiusClapeyron {
    /// Latent heat of transition (J/mol).
    pub latent_heat: f64,
    /// Reference transition temperature at P₀ (K).
    pub t0: f64,
    /// Reference pressure P₀ (Pa).
    pub p0: f64,
    /// Molar volume change ΔV = V_liquid - V_solid (m³/mol).
    pub delta_v: f64,
}
impl ClausiusClapeyron {
    /// Create a new Clausius-Clapeyron model.
    pub fn new(latent_heat: f64, t0: f64, p0: f64, delta_v: f64) -> Self {
        Self {
            latent_heat,
            t0,
            p0,
            delta_v,
        }
    }
    /// Water/ice approximation: L=6010 J/mol, T₀=273.15 K, P₀=101325 Pa,
    /// ΔV ≈ −1.63×10⁻⁶ m³/mol (ice is less dense than water → negative slope).
    pub fn water_ice() -> Self {
        Self::new(6010.0, 273.15, 101325.0, -1.63e-6)
    }
    /// dP/dT slope of the phase boundary at temperature T.
    pub fn slope(&self, t: f64) -> f64 {
        if t.abs() < 1e-10 || self.delta_v.abs() < 1e-30 {
            return 0.0;
        }
        self.latent_heat / (t * self.delta_v)
    }
    /// Melting temperature at pressure P using linear approximation.
    ///
    /// T_m(P) ≈ T₀ + (P - P₀) / slope(T₀)
    pub fn melting_temperature(&self, pressure: f64) -> f64 {
        let slope0 = self.slope(self.t0);
        if slope0.abs() < 1e-30 {
            return self.t0;
        }
        self.t0 + (pressure - self.p0) / slope0
    }
    /// Equilibrium pressure at a given temperature.
    ///
    /// P_eq(T) ≈ P₀ + slope(T₀) * (T - T₀)
    pub fn equilibrium_pressure(&self, temperature: f64) -> f64 {
        let slope0 = self.slope(self.t0);
        self.p0 + slope0 * (temperature - self.t0)
    }
    /// Latent heat at temperature T using the Kirchhoff approximation.
    ///
    /// L(T) ≈ L(T₀) + ΔCp * (T - T₀)
    pub fn latent_heat_at(&self, temperature: f64, delta_cp: f64) -> f64 {
        self.latent_heat + delta_cp * (temperature - self.t0)
    }
}
/// Scheil-Gulliver non-equilibrium solidification model.
///
/// Assumes complete mixing in the liquid and no diffusion in the solid.
/// The composition of the liquid evolves as solidification proceeds:
///
/// `C_L(f_s) = C_0 * (1 - f_s)^(k - 1)`
///
/// where:
/// - `C_0` is the initial alloy composition (wt% or mol fraction)
/// - `f_s` is the fraction solidified \[0, 1)
/// - `k`   is the equilibrium partition coefficient C_solid/C_liquid
///
/// Ref: Scheil (1942), Gulliver (1913).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ScheilSolidification {
    /// Initial alloy composition C₀ (wt% or mol fraction).
    pub c0: f64,
    /// Equilibrium partition coefficient k = C_solid / C_liquid.
    pub k_partition: f64,
}
impl ScheilSolidification {
    /// Create a new Scheil model.
    pub fn new(c0: f64, k_partition: f64) -> Self {
        Self { c0, k_partition }
    }
    /// Binary Al-Cu alloy (k ≈ 0.17 for Cu in Al at eutectic temperature).
    pub fn al_cu() -> Self {
        Self::new(4.0, 0.17)
    }
    /// Fe-C system (k ≈ 0.17 for C in Fe).
    pub fn fe_c() -> Self {
        Self::new(0.4, 0.17)
    }
    /// Liquid composition \[wt%\] at solid fraction f_s.
    ///
    /// `C_L = C_0 * (1 - f_s)^(k-1)`
    pub fn liquid_composition(&self, f_s: f64) -> f64 {
        let f_s = f_s.clamp(0.0, 1.0 - 1.0e-12);
        self.c0 * (1.0 - f_s).powf(self.k_partition - 1.0)
    }
    /// Solid composition at the solidification front \[wt%\].
    ///
    /// `C_S* = k * C_L`
    pub fn solid_composition_at_front(&self, f_s: f64) -> f64 {
        self.k_partition * self.liquid_composition(f_s)
    }
    /// Mean solid composition up to fraction f_s (integrated Scheil).
    ///
    /// `C_S_bar = C_0 * (1 - (1 - f_s)^k) / f_s`
    pub fn mean_solid_composition(&self, f_s: f64) -> f64 {
        let f_s = f_s.clamp(1.0e-12, 1.0 - 1.0e-12);
        self.c0 * (1.0 - (1.0 - f_s).powf(self.k_partition)) / f_s
    }
    /// Solid fraction at which the eutectic composition C_eut is reached.
    ///
    /// `f_s_eut = 1 - (C_eut / C_0)^(1/(k-1))`
    ///
    /// Returns `None` if k = 1 (no segregation) or C_eut < C_0.
    pub fn eutectic_solid_fraction(&self, c_eutectic: f64) -> Option<f64> {
        if (self.k_partition - 1.0).abs() < 1e-12 {
            return None;
        }
        if c_eutectic <= self.c0 {
            return None;
        }
        let exponent = 1.0 / (self.k_partition - 1.0);
        let f_s = 1.0 - (c_eutectic / self.c0).powf(exponent);
        Some(f_s.clamp(0.0, 1.0))
    }
    /// Solidification range: temperature span between liquidus and eutectic.
    ///
    /// Using the linearised phase diagram with liquidus slope m_L \[°C/wt%\]:
    /// `ΔT = m_L * C_L(f_s=0)` to `m_L * C_eut`
    ///
    /// Returns (T_liquidus, T_eutectic) given melting point T_melt,
    /// liquidus slope m_L \[°C/wt%\], and eutectic composition C_eut.
    pub fn solidification_range(
        &self,
        t_melt: f64,
        m_liquidus: f64,
        c_eutectic: f64,
    ) -> (f64, f64) {
        let t_liquidus = t_melt + m_liquidus * self.c0;
        let t_eutectic = t_melt + m_liquidus * c_eutectic;
        (t_liquidus, t_eutectic)
    }
}
/// Transformation plasticity model (Greenwood-Johnson mechanism).
///
/// The instantaneous transformation plasticity strain rate:
///   dε_tp/dt = K_tp / σ_y * σ_dev * df/dt
///
/// where K_tp is a material constant, σ_y is the yield stress, σ_dev is the
/// deviatoric stress, and df/dt is the phase fraction rate.
#[derive(Debug, Clone)]
pub struct TransformationPlasticity {
    /// Transformation plasticity constant K_tp (dimensionless).
    pub k_tp: f64,
    /// Yield stress of the weaker phase (Pa).
    pub yield_stress: f64,
    /// Accumulated transformation plastic strain.
    pub eps_tp: f64,
}
impl TransformationPlasticity {
    /// Create a new transformation plasticity model.
    pub fn new(k_tp: f64, yield_stress: f64) -> Self {
        Self {
            k_tp,
            yield_stress,
            eps_tp: 0.0,
        }
    }
    /// Increment the transformation plastic strain.
    ///
    /// `sigma_dev`: deviatoric stress magnitude (Pa).
    /// `df`: phase fraction increment in this step.
    pub fn increment(&mut self, sigma_dev: f64, df: f64) -> f64 {
        if df <= 0.0 || self.yield_stress < 1e-15 {
            return 0.0;
        }
        let de = self.k_tp / self.yield_stress * sigma_dev * df;
        self.eps_tp += de;
        de
    }
    /// Total transformation plastic strain accumulated.
    pub fn total_strain(&self) -> f64 {
        self.eps_tp
    }
}
/// Standalone Koistinen-Marburger helper for austenite-to-martensite kinetics.
///
/// f_m(T) = 1 − exp(−k · (Ms − T))  for T < Ms
#[derive(Debug, Clone)]
pub struct KoistinenMarburger {
    /// Martensite start temperature (K).
    pub ms: f64,
    /// Rate constant k (K⁻¹).
    pub k: f64,
}
impl KoistinenMarburger {
    /// Create a new KM model.
    pub fn new(ms: f64, k: f64) -> Self {
        Self { ms, k }
    }
    /// Martensite volume fraction at temperature `t`.
    pub fn fraction(&self, t: f64) -> f64 {
        if t >= self.ms {
            0.0
        } else {
            (1.0 - E.powf(-self.k * (self.ms - t))).min(1.0)
        }
    }
    /// Temperature at which a given fraction `f` is reached.
    pub fn temperature_for_fraction(&self, f: f64) -> f64 {
        let f_clamped = f.clamp(0.0, 1.0 - f64::EPSILON);
        self.ms + f64::ln(1.0 - f_clamped) / self.k
    }
}
/// Simplified TRIP (Transformation-Induced Plasticity) effect model.
///
/// The additional TRIP strain increment is:
///   dε_TRIP = K_TRIP · σ_eq · (1 − f_m) · df_m
///
/// where σ_eq is the equivalent stress and df_m is the martensite increment.
#[derive(Debug, Clone)]
pub struct TripEffect {
    /// TRIP constant K (Pa⁻¹).
    pub k_trip: f64,
    /// Accumulated TRIP strain (dimensionless).
    pub accumulated_strain: f64,
}
impl TripEffect {
    /// Create a new TRIP model.
    pub fn new(k_trip: f64) -> Self {
        Self {
            k_trip,
            accumulated_strain: 0.0,
        }
    }
    /// Compute and accumulate TRIP strain increment.
    ///
    /// `sigma_eq`: equivalent (von Mises) stress (Pa).
    /// `f_m`: current martensite fraction.
    /// `df_m`: martensite fraction increment in this step.
    pub fn trip_strain_increment(&mut self, sigma_eq: f64, f_m: f64, df_m: f64) -> f64 {
        if df_m <= 0.0 {
            return 0.0;
        }
        let de = self.k_trip * sigma_eq * (1.0 - f_m) * df_m;
        self.accumulated_strain += de;
        de
    }
    /// Reset accumulated strain (e.g., after a full cycle).
    pub fn reset(&mut self) {
        self.accumulated_strain = 0.0;
    }
}
/// Tracks the volume fractions of multiple phases (austenite, martensite,
/// bainite, pearlite, ferrite) as functions of time and temperature.
///
/// Fractions are constrained to sum to ≤ 1.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct VolumeFractionTracker {
    /// Phase labels.
    pub phases: Vec<String>,
    /// Current volume fractions (index-aligned with `phases`).
    pub fractions: Vec<f64>,
    /// History: (time, temperature, fractions snapshot).
    pub history: Vec<(f64, f64, Vec<f64>)>,
}
impl VolumeFractionTracker {
    /// Create a tracker with given phase labels.
    /// Initial state: first phase = 1.0, rest = 0.0 (fully austenitic).
    pub fn new(phases: Vec<String>) -> Self {
        let n = phases.len();
        let mut fractions = vec![0.0; n];
        if n > 0 {
            fractions[0] = 1.0;
        }
        Self {
            phases,
            fractions,
            history: Vec::new(),
        }
    }
    /// Update phase fractions (normalised to sum ≤ 1).
    #[allow(clippy::needless_range_loop)]
    pub fn update(&mut self, new_fractions: Vec<f64>, time: f64, temperature: f64) {
        let n = self.fractions.len().min(new_fractions.len());
        for i in 0..n {
            self.fractions[i] = new_fractions[i].clamp(0.0, 1.0);
        }
        self.history
            .push((time, temperature, self.fractions.clone()));
    }
    /// Total volume fraction (should be ≤ 1).
    pub fn total(&self) -> f64 {
        self.fractions.iter().sum()
    }
    /// Volume fraction of phase at index `i`.
    pub fn fraction(&self, i: usize) -> f64 {
        self.fractions.get(i).copied().unwrap_or(0.0)
    }
    /// Find the dominant phase index.
    pub fn dominant_phase(&self) -> usize {
        self.fractions
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
}
/// TTT C-curve for a single transformation using parabolic Avrami.
///
/// The incubation time τ(T) follows a C-curve shape:
/// `τ(T) = τ_nose * exp(A * (T - T_nose)²)`
///
/// where A is the curvature parameter.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TttCcurve {
    /// Nose temperature T_nose \[K\]
    pub t_nose: f64,
    /// Nose time (incubation time at T_nose) \[s\]
    pub tau_nose: f64,
    /// C-curve curvature parameter A \[K⁻²\]
    pub a_curvature: f64,
    /// Avrami exponent n
    pub n: f64,
}
impl TttCcurve {
    /// Create a new TTT C-curve.
    pub fn new(t_nose: f64, tau_nose: f64, a_curvature: f64, n: f64) -> Self {
        Self {
            t_nose,
            tau_nose,
            a_curvature,
            n,
        }
    }
    /// Incubation time \[s\] at temperature T.
    pub fn incubation_time(&self, temperature: f64) -> f64 {
        let dt = temperature - self.t_nose;
        self.tau_nose * E.powf(self.a_curvature * dt * dt)
    }
    /// Avrami rate constant k at temperature T.
    ///
    /// `k = -ln(0.99) / τ(T)^n`   (1% = start of transformation)
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        let tau = self.incubation_time(temperature);
        if tau < 1e-30 {
            return f64::INFINITY;
        }
        (-f64::ln(0.99)) / tau.powf(self.n)
    }
    /// Transformed fraction at time t and temperature T.
    pub fn fraction(&self, time: f64, temperature: f64) -> f64 {
        if time <= 0.0 {
            return 0.0;
        }
        let k = self.rate_constant(temperature);
        (1.0 - E.powf(-k * time.powf(self.n))).clamp(0.0, 1.0)
    }
    /// Time to reach fraction f at temperature T.
    ///
    /// `t = (-ln(1-f) / k)^(1/n)`
    pub fn time_to_fraction(&self, fraction: f64, temperature: f64) -> f64 {
        let f = fraction.clamp(1e-12, 1.0 - 1e-12);
        let k = self.rate_constant(temperature);
        if k <= 0.0 || !k.is_finite() {
            return f64::INFINITY;
        }
        ((-f64::ln(1.0 - f)) / k).powf(1.0 / self.n)
    }
    /// Fastest transformation temperature: minimum of τ(T) curve.
    ///
    /// This is the nose temperature T_nose by construction.
    pub fn nose_temperature(&self) -> f64 {
        self.t_nose
    }
}
/// Phase-dependent mechanical properties interpolation.
///
/// Given the volume fractions of martensite, bainite, ferrite, and austenite,
/// compute the effective Young's modulus, yield strength, and hardness by
/// the rule of mixtures.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PhaseDependentProperties {
    /// Young's moduli \[Pa\] for each phase: \[martensite, bainite, ferrite, austenite\].
    pub youngs_moduli: [f64; 4],
    /// Yield strengths \[Pa\] for each phase.
    pub yield_strengths: [f64; 4],
    /// Vickers hardness \[HV\] for each phase.
    pub hardness: [f64; 4],
}
impl PhaseDependentProperties {
    /// Create phase-dependent properties.
    pub fn new(youngs_moduli: [f64; 4], yield_strengths: [f64; 4], hardness: [f64; 4]) -> Self {
        Self {
            youngs_moduli,
            yield_strengths,
            hardness,
        }
    }
    /// Typical medium-carbon steel phase properties.
    pub fn medium_carbon_steel() -> Self {
        Self {
            youngs_moduli: [205e9, 200e9, 210e9, 195e9],
            yield_strengths: [1500e6, 800e6, 350e6, 250e6],
            hardness: [700.0, 400.0, 150.0, 200.0],
        }
    }
    /// Effective Young's modulus by rule of mixtures.
    pub fn effective_youngs_modulus(&self, fractions: &[f64; 4]) -> f64 {
        self.youngs_moduli
            .iter()
            .zip(fractions.iter())
            .map(|(&e, &f)| e * f)
            .sum()
    }
    /// Effective yield strength by rule of mixtures.
    pub fn effective_yield_strength(&self, fractions: &[f64; 4]) -> f64 {
        self.yield_strengths
            .iter()
            .zip(fractions.iter())
            .map(|(&sy, &f)| sy * f)
            .sum()
    }
    /// Effective Vickers hardness by rule of mixtures.
    pub fn effective_hardness(&self, fractions: &[f64; 4]) -> f64 {
        self.hardness
            .iter()
            .zip(fractions.iter())
            .map(|(&hv, &f)| hv * f)
            .sum()
    }
    /// Estimate ultimate tensile strength from hardness (Tabor relation):
    /// UTS ≈ 3.45 · HV \[MPa\].
    pub fn uts_from_hardness(&self, fractions: &[f64; 4]) -> f64 {
        3.45e6 * self.effective_hardness(fractions)
    }
}
/// Tracks the evolution of multiple phase fractions over a simulation.
///
/// Stores a time-series of each phase fraction and allows querying the
/// current state, computing phase growth/dissolution rates, and exporting data.
#[derive(Debug, Clone)]
pub struct PhaseFractionEvolution {
    /// Number of phases tracked.
    pub n_phases: usize,
    /// Time history (s).
    pub times: Vec<f64>,
    /// Phase fractions history: `fractions[i]` is a Vec over time of phase i fraction.
    pub fractions: Vec<Vec<f64>>,
}
impl PhaseFractionEvolution {
    /// Create a new tracker for `n_phases` phases.
    pub fn new(n_phases: usize) -> Self {
        Self {
            n_phases,
            times: Vec::new(),
            fractions: vec![Vec::new(); n_phases],
        }
    }
    /// Record the current fractions at time `t`.
    ///
    /// `fracs` must have length `n_phases` and sum ≤ 1.
    pub fn record(&mut self, t: f64, fracs: &[f64]) {
        assert_eq!(fracs.len(), self.n_phases);
        self.times.push(t);
        for (i, &f) in fracs.iter().enumerate() {
            self.fractions[i].push(f.clamp(0.0, 1.0));
        }
    }
    /// Current (last recorded) fraction of phase `i`.
    pub fn current(&self, phase: usize) -> f64 {
        self.fractions
            .get(phase)
            .and_then(|v| v.last())
            .cloned()
            .unwrap_or(0.0)
    }
    /// Approximate rate of change of phase `phase` (fraction/s).
    pub fn rate(&self, phase: usize) -> f64 {
        let v = &self.fractions[phase];
        let n = v.len();
        if n < 2 {
            return 0.0;
        }
        let dt = self.times[n - 1] - self.times[n - 2];
        if dt.abs() < 1e-15 {
            return 0.0;
        }
        (v[n - 1] - v[n - 2]) / dt
    }
    /// Sum of all current phase fractions (should be ≤ 1 in a well-defined model).
    pub fn total_fraction(&self) -> f64 {
        (0..self.n_phases).map(|i| self.current(i)).sum()
    }
}
/// Multi-component lever rule for two-phase equilibrium in an alloy.
///
/// Determines the phase fractions and compositions in two-phase coexistence
/// using mass balance (lever rule).
#[derive(Debug, Clone)]
pub struct LeverRule {
    /// Composition of solute in the α phase (mole fraction).
    pub x_alpha: f64,
    /// Composition of solute in the β phase (mole fraction).
    pub x_beta: f64,
}
impl LeverRule {
    /// Create a lever rule model with given phase compositions.
    pub fn new(x_alpha: f64, x_beta: f64) -> Self {
        Self { x_alpha, x_beta }
    }
    /// Phase fraction of β phase from overall composition x.
    ///
    /// f_β = (x - x_α) / (x_β - x_α)
    pub fn beta_fraction(&self, x_overall: f64) -> f64 {
        let denom = self.x_beta - self.x_alpha;
        if denom.abs() < 1e-30 {
            return 0.0;
        }
        ((x_overall - self.x_alpha) / denom).clamp(0.0, 1.0)
    }
    /// Phase fraction of α phase.
    pub fn alpha_fraction(&self, x_overall: f64) -> f64 {
        1.0 - self.beta_fraction(x_overall)
    }
    /// Solidus composition (minimum β-fraction lever endpoint).
    pub fn solidus(&self) -> f64 {
        self.x_alpha
    }
    /// Liquidus composition (maximum β-fraction lever endpoint).
    pub fn liquidus(&self) -> f64 {
        self.x_beta
    }
    /// Check if overall composition is in two-phase region.
    pub fn is_two_phase(&self, x_overall: f64) -> bool {
        let lo = self.x_alpha.min(self.x_beta);
        let hi = self.x_alpha.max(self.x_beta);
        x_overall > lo && x_overall < hi
    }
    /// Tie-line length (compositional separation between phases).
    pub fn tie_line_length(&self) -> f64 {
        (self.x_beta - self.x_alpha).abs()
    }
}
/// Simplified JMatPro-style model for martensite and bainite start/finish.
#[derive(Debug, Clone)]
pub struct JmatProModel {
    /// Martensite start temperature (K).
    pub ms_temp: f64,
    /// Martensite finish temperature (K).
    pub mf_temp: f64,
    /// Bainite start temperature (K).
    pub bs_temp: f64,
}
impl JmatProModel {
    /// Create a new JMatPro model.
    pub fn new(ms_temp: f64, mf_temp: f64, bs_temp: f64) -> Self {
        Self {
            ms_temp,
            mf_temp,
            bs_temp,
        }
    }
    /// Martensite volume fraction via Koistinen-Marburger equation.
    ///
    /// X_m = 1 - exp(-0.011 * (Ms - T))  for T < Ms, else 0.
    ///
    /// # Arguments
    /// * `t`  - Current temperature (K).
    /// * `ms` - Martensite start temperature (K).
    pub fn martensite_fraction_koistinen_marburger(t: f64, ms: f64) -> f64 {
        if t >= ms {
            0.0
        } else {
            (1.0 - E.powf(-0.011 * (ms - t))).clamp(0.0, 1.0)
        }
    }
    /// Whether `t` is in the martensite transformation range \[Mf, Ms\].
    pub fn in_martensite_range(&self, t: f64) -> bool {
        t >= self.mf_temp && t <= self.ms_temp
    }
    /// Whether `t` is above the bainite start temperature.
    pub fn above_bainite_start(&self, t: f64) -> bool {
        t >= self.bs_temp
    }
}
/// Johnson-Mehl-Avrami-Kolmogorov (JMAK) isothermal transformation kinetics.
///
/// f(t) = 1 - exp(-k · t^n)
#[derive(Debug, Clone)]
pub struct JohnsonMehlAvramiKolmogorov {
    /// Rate constant k (s⁻ⁿ).
    pub k: f64,
    /// Avrami exponent n (dimensionless; typically 1–4).
    pub n: f64,
}
impl JohnsonMehlAvramiKolmogorov {
    /// Create a new JMAK model.
    pub fn new(k: f64, n: f64) -> Self {
        Self { k, n }
    }
    /// Transformed fraction at time `t`.
    ///
    /// f = 1 - exp(-k · t^n)
    pub fn fraction_transformed(&self, time: f64) -> f64 {
        if time <= 0.0 {
            return 0.0;
        }
        1.0 - E.powf(-self.k * time.powf(self.n))
    }
    /// Instantaneous transformation rate df/dt.
    ///
    /// df/dt = k · n · t^(n-1) · exp(-k · t^n)
    pub fn transformation_rate(&self, time: f64) -> f64 {
        if time <= 0.0 {
            return 0.0;
        }
        let kt_n = self.k * time.powf(self.n);
        self.k * self.n * time.powf(self.n - 1.0) * E.powf(-kt_n)
    }
    /// Time required to reach a given `fraction` transformed.
    ///
    /// t = (-ln(1 - f) / k)^(1/n)
    pub fn incubation_time(&self, fraction: f64) -> f64 {
        let f = fraction.clamp(0.0, 1.0 - f64::EPSILON);
        ((-f64::ln(1.0 - f)) / self.k).powf(1.0 / self.n)
    }
}
/// Pearlite transformation model using Avrami kinetics.
///
/// Very similar to bainite but with different kinetic parameters and
/// applicable temperature range (above bainite nose).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PearliteTransformation {
    /// Avrami exponent n.
    pub n: f64,
    /// Pre-exponential constant k0 \[s^-n\].
    pub k0: f64,
    /// Activation energy \[J/mol\].
    pub q_p: f64,
    /// Pearlite start temperature P_s \[K\].
    pub p_s: f64,
    /// Eutectoid temperature A_e1 \[K\] — no pearlite above this.
    pub a_e1: f64,
}
impl PearliteTransformation {
    /// Create a new pearlite transformation model.
    pub fn new(n: f64, k0: f64, q_p: f64, p_s: f64, a_e1: f64) -> Self {
        Self {
            n,
            k0,
            q_p,
            p_s,
            a_e1,
        }
    }
    /// Avrami rate constant at temperature T.
    pub fn rate_constant(&self, temperature: f64) -> f64 {
        self.k0 * (-self.q_p / (R_GAS * temperature)).exp()
    }
    /// Pearlite volume fraction at time t and temperature T.
    pub fn fraction_at(&self, time: f64, temperature: f64) -> f64 {
        if temperature > self.a_e1 || temperature < self.p_s || time <= 0.0 {
            return 0.0;
        }
        let k = self.rate_constant(temperature);
        1.0 - (-k * time.powf(self.n)).exp()
    }
    /// Maximum pearlite fraction (constrained to ≤ 1).
    pub fn max_fraction_at(&self, time: f64, temperature: f64) -> f64 {
        self.fraction_at(time, temperature).min(1.0)
    }
    /// Incubation time (1% transformed) at temperature T.
    pub fn incubation_time(&self, temperature: f64) -> f64 {
        let k = self.rate_constant(temperature);
        if k <= 0.0 {
            return f64::INFINITY;
        }
        (0.010536 / k).powf(1.0 / self.n)
    }
}
/// Extended TTT diagram with Avrami kinetics for multiple transformation products.
///
/// Stores separate C-curve parameters for:
/// - Pearlite (diffusion-controlled)
/// - Bainite (mixed mechanism)
/// - Martensite (displacive, instantaneous at Ms)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TttExtended {
    /// Nose temperature for pearlite C-curve \[K\]
    pub t_nose_pearlite: f64,
    /// Nose time for pearlite (1% start) \[s\]
    pub t_nose_pearlite_time: f64,
    /// Nose temperature for bainite C-curve \[K\]
    pub t_nose_bainite: f64,
    /// Nose time for bainite (1% start) \[s\]
    pub t_nose_bainite_time: f64,
    /// Martensite start temperature Ms \[K\]
    pub ms: f64,
    /// Eutectoid temperature A_e1 \[K\]
    pub a_e1: f64,
    /// Avrami exponent for pearlite
    pub n_pearlite: f64,
    /// Avrami exponent for bainite
    pub n_bainite: f64,
}
impl TttExtended {
    /// Create a new extended TTT diagram.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        t_nose_pearlite: f64,
        t_nose_pearlite_time: f64,
        t_nose_bainite: f64,
        t_nose_bainite_time: f64,
        ms: f64,
        a_e1: f64,
        n_pearlite: f64,
        n_bainite: f64,
    ) -> Self {
        Self {
            t_nose_pearlite,
            t_nose_pearlite_time,
            t_nose_bainite,
            t_nose_bainite_time,
            ms,
            a_e1,
            n_pearlite,
            n_bainite,
        }
    }
    /// Avrami rate constant k for pearlite at temperature T.
    ///
    /// Parabolic C-curve: `k = ln(0.99) / -(t_start(T)^n)`
    /// where `t_start` is estimated from the C-curve shape.
    fn pearlite_rate_constant(&self, temperature: f64) -> f64 {
        let dt = (temperature - self.t_nose_pearlite).powi(2);
        let t_start = self.t_nose_pearlite_time * E.powf(0.001 * dt);
        if t_start < 1e-30 {
            return f64::INFINITY;
        }
        (-f64::ln(0.99)) / t_start.powf(self.n_pearlite)
    }
    /// Fraction of pearlite at time t and temperature T.
    pub fn pearlite_fraction(&self, time: f64, temperature: f64) -> f64 {
        if temperature >= self.a_e1 || temperature <= self.ms || time <= 0.0 {
            return 0.0;
        }
        let k = self.pearlite_rate_constant(temperature);
        (1.0 - E.powf(-k * time.powf(self.n_pearlite))).clamp(0.0, 1.0)
    }
    /// Fraction of bainite at time t and temperature T.
    pub fn bainite_fraction(&self, time: f64, temperature: f64) -> f64 {
        if temperature >= self.t_nose_pearlite || temperature <= self.ms || time <= 0.0 {
            return 0.0;
        }
        let dt = (temperature - self.t_nose_bainite).powi(2);
        let t_start = self.t_nose_bainite_time * E.powf(0.001 * dt);
        if t_start < 1e-30 {
            return 1.0;
        }
        let k = (-f64::ln(0.99)) / t_start.powf(self.n_bainite);
        (1.0 - E.powf(-k * time.powf(self.n_bainite))).clamp(0.0, 1.0)
    }
    /// Martensite fraction via Koistinen-Marburger at temperature T.
    pub fn martensite_fraction(&self, temperature: f64) -> f64 {
        if temperature >= self.ms {
            return 0.0;
        }
        (1.0 - E.powf(-0.011 * (self.ms - temperature))).min(1.0)
    }
    /// Dominant transformation product at (time, temperature).
    ///
    /// Returns: `"austenite"`, `"pearlite"`, `"bainite"`, `"martensite"`, or `"mixed"`.
    pub fn dominant_product(&self, time: f64, temperature: f64) -> &'static str {
        if temperature > self.a_e1 {
            return "austenite";
        }
        if temperature < self.ms {
            return "martensite";
        }
        let fp = self.pearlite_fraction(time, temperature);
        let fb = self.bainite_fraction(time, temperature);
        if fp > 0.5 {
            return "pearlite";
        }
        if fb > 0.5 {
            return "bainite";
        }
        if fp + fb < 0.05 {
            return "austenite";
        }
        "mixed"
    }
}
/// Continuous Cooling Transformation (CCT) diagram data.
///
/// Stores the start and finish temperatures and cooling rate boundaries
/// for each transformation product (ferrite, pearlite, bainite, martensite).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CctDiagram {
    /// Cooling rate at which bainite transformation starts \[K/s\].
    pub bainite_start_rate: f64,
    /// Cooling rate at which bainite transformation ends \[K/s\].
    pub bainite_end_rate: f64,
    /// Cooling rate above which martensite forms \[K/s\].
    pub martensite_critical_rate: f64,
    /// Martensite start temperature \[K\].
    pub ms_temperature: f64,
    /// Pearlite nose temperature \[K\].
    pub pearlite_nose_temperature: f64,
    /// Pearlite nose time \[s\].
    pub pearlite_nose_time: f64,
}
impl CctDiagram {
    /// Create a new CCT diagram.
    pub fn new(
        bainite_start_rate: f64,
        bainite_end_rate: f64,
        martensite_critical_rate: f64,
        ms_temperature: f64,
        pearlite_nose_temperature: f64,
        pearlite_nose_time: f64,
    ) -> Self {
        Self {
            bainite_start_rate,
            bainite_end_rate,
            martensite_critical_rate,
            ms_temperature,
            pearlite_nose_temperature,
            pearlite_nose_time,
        }
    }
    /// Predict dominant microstructure for a given cooling rate \[K/s\].
    ///
    /// Returns: `"martensite"`, `"bainite"`, `"pearlite"`, or `"mixed"`.
    pub fn predict_microstructure(&self, cooling_rate: f64) -> &'static str {
        if cooling_rate >= self.martensite_critical_rate {
            "martensite"
        } else if cooling_rate >= self.bainite_end_rate {
            "bainite"
        } else if cooling_rate >= self.bainite_start_rate {
            "mixed"
        } else {
            "pearlite"
        }
    }
    /// Estimate martensite fraction for a given cooling rate (simplified).
    pub fn martensite_fraction(&self, cooling_rate: f64) -> f64 {
        if cooling_rate >= self.martensite_critical_rate {
            return 1.0;
        }
        if cooling_rate <= self.bainite_end_rate {
            return 0.0;
        }
        (cooling_rate - self.bainite_end_rate)
            / (self.martensite_critical_rate - self.bainite_end_rate)
    }
}
/// Simplified Brinson shape memory alloy (SMA) model.
///
/// Tracks stress-induced (`ξs`) and temperature-induced (`ξT`) martensite
/// fractions and computes the SMA stress-strain response.
///
/// Phase fractions: ξ = ξs + ξT (total martensite fraction).
#[derive(Debug, Clone)]
pub struct BrinsonModel {
    /// Young's modulus of austenite (Pa).
    pub e_austenite: f64,
    /// Young's modulus of martensite (Pa).
    pub e_martensite: f64,
    /// Martensite start temperature on cooling (K).
    pub ms: f64,
    /// Martensite finish temperature on cooling (K).
    pub mf: f64,
    /// Austenite start temperature on heating (K).
    pub as_temp: f64,
    /// Austenite finish temperature on heating (K).
    pub af_temp: f64,
    /// Transformation strain (maximum).
    pub eps_l: f64,
    /// Stress-induced martensite fraction.
    pub xi_s: f64,
    /// Temperature-induced martensite fraction.
    pub xi_t: f64,
}
impl BrinsonModel {
    /// Create a new Brinson SMA model.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        e_austenite: f64,
        e_martensite: f64,
        ms: f64,
        mf: f64,
        as_temp: f64,
        af_temp: f64,
        eps_l: f64,
    ) -> Self {
        Self {
            e_austenite,
            e_martensite,
            ms,
            mf,
            as_temp,
            af_temp,
            eps_l,
            xi_s: 0.0,
            xi_t: 0.0,
        }
    }
    /// Total martensite fraction.
    pub fn xi(&self) -> f64 {
        (self.xi_s + self.xi_t).clamp(0.0, 1.0)
    }
    /// Effective Young's modulus (linear mixture).
    pub fn effective_modulus(&self) -> f64 {
        let xi = self.xi();
        self.e_martensite * xi + self.e_austenite * (1.0 - xi)
    }
    /// Update temperature-induced martensite at temperature `t`.
    pub fn update_temperature(&mut self, t: f64) {
        if t <= self.mf {
            self.xi_t = 1.0;
        } else if t < self.ms {
            self.xi_t =
                0.5 * (1.0 - ((std::f64::consts::PI * (t - self.ms) / (self.mf - self.ms)).cos()));
        } else if t >= self.af_temp {
            self.xi_t = 0.0;
        } else if t > self.as_temp {
            self.xi_t = (self.xi_t
                * 0.5
                * (1.0
                    + ((std::f64::consts::PI * (t - self.as_temp)
                        / (self.af_temp - self.as_temp))
                        .cos())))
            .max(0.0);
        }
    }
    /// Compute stress given mechanical strain `eps` and temperature `t`.
    pub fn stress(&self, eps: f64) -> f64 {
        let e = self.effective_modulus();
        let xi = self.xi();
        e * (eps - self.eps_l * self.xi_s - self.eps_l * xi * 0.1)
    }

    /// Create a Brinson model with typical Nitinol (NiTi) parameters.
    ///
    /// * E_A = 75 GPa, E_M = 28 GPa
    /// * Ms = 291 K, Mf = 273 K, As = 307 K, Af = 325 K
    /// * ε_L = 0.067 (maximum transformation strain)
    pub fn nitinol() -> Self {
        Self::new(
            75.0e9, // E_austenite
            28.0e9, // E_martensite
            291.0,  // Ms
            273.0,  // Mf
            307.0,  // As
            325.0,  // Af
            0.067,  // eps_l
        )
    }
}
/// Time-Temperature-Transformation (TTT) diagram for isothermal
/// transformation of supercooled austenite.
///
/// The C-curve is parameterised by three points:
/// - nose (T_nose, t_nose): fastest transformation
/// - high-temperature arm
/// - martensite start temperature M_s
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct TttDiagram {
    /// Nose temperature \[K\].
    pub t_nose: f64,
    /// Nose time (1% transformation start) \[s\].
    pub t_nose_time: f64,
    /// Eutectoid temperature A_e1 \[K\].
    pub a_e1: f64,
    /// Martensite start temperature M_s \[K\].
    pub ms: f64,
    /// Bainite start temperature B_s \[K\].
    pub b_s: f64,
}
impl TttDiagram {
    /// Create a new TTT diagram.
    pub fn new(t_nose: f64, t_nose_time: f64, a_e1: f64, ms: f64, b_s: f64) -> Self {
        Self {
            t_nose,
            t_nose_time,
            a_e1,
            ms,
            b_s,
        }
    }
    /// Estimate the start time of transformation at temperature T (1% fraction).
    ///
    /// Uses a parabolic C-curve: ln(t) = ln(t_nose) + c·(T - T_nose)^2
    pub fn start_time(&self, temperature: f64) -> f64 {
        if temperature >= self.a_e1 || temperature <= self.ms {
            return f64::INFINITY;
        }
        let c =
            4.0 * (temperature - self.t_nose).powi(2) / ((self.a_e1 - self.t_nose).powi(2) + 1.0);
        self.t_nose_time * E.powf(c * 0.01)
    }
    /// Check if transformation occurs before quench at given cooling rate.
    pub fn transformation_avoided(&self, cooling_rate: f64) -> bool {
        let t_at_nose = self.t_nose_time;
        let time_to_cool_through_nose = (self.a_e1 - self.t_nose) / cooling_rate;
        time_to_cool_through_nose < t_at_nose
    }
    /// Temperature range of the bainite transformation region.
    pub fn bainite_range(&self) -> (f64, f64) {
        (self.ms, self.b_s)
    }
    /// Temperature range of the pearlite transformation region.
    pub fn pearlite_range(&self) -> (f64, f64) {
        (self.t_nose, self.a_e1)
    }
}
/// Martensitic phase transformation model with hysteresis.
#[derive(Debug, Clone)]
pub struct MartensiteTransformation {
    /// Martensite-start temperature (K).
    pub ms_temperature: f64,
    /// Martensite-finish temperature (K).
    pub mf_temperature: f64,
    /// Austenite-start temperature on reverse transformation (K).
    pub as_temperature: f64,
    /// Austenite-finish temperature (K).
    pub af_temperature: f64,
    /// Koistinen-Marburger rate constant (typically ~0.011 K⁻¹).
    pub k_km: f64,
    /// Maximum transformation strain (dimensionless).
    pub eps_max: f64,
}
impl MartensiteTransformation {
    /// Create a new `MartensiteTransformation` model.
    pub fn new(ms: f64, mf: f64, as_temp: f64, af: f64, k_km: f64, eps_max: f64) -> Self {
        Self {
            ms_temperature: ms,
            mf_temperature: mf,
            as_temperature: as_temp,
            af_temperature: af,
            k_km,
            eps_max,
        }
    }
    /// Martensite volume fraction via the Koistinen-Marburger equation.
    ///
    /// f_m = 1 - exp(-k·(Ms - T))  for T < Ms, else 0.
    pub fn martensite_fraction_koistinen_marburger(&self, t: f64) -> f64 {
        if t >= self.ms_temperature {
            0.0
        } else {
            let driving = self.ms_temperature - t;
            (1.0 - E.powf(-self.k_km * driving)).min(1.0)
        }
    }
    /// Austenite volume fraction on heating (linear between As and Af).
    pub fn austenite_fraction(&self, t: f64) -> f64 {
        if t <= self.as_temperature {
            0.0
        } else if t >= self.af_temperature {
            1.0
        } else {
            (t - self.as_temperature) / (self.af_temperature - self.as_temperature)
        }
    }
    /// Whether temperature `t` falls inside the hysteresis loop.
    ///
    /// On cooling (`cooling = true`) the active window is \[Mf, Ms\].
    /// On heating (`cooling = false`) the active window is \[As, Af\].
    pub fn is_in_hysteresis(&self, t: f64, cooling: bool) -> bool {
        if cooling {
            t >= self.mf_temperature && t <= self.ms_temperature
        } else {
            t >= self.as_temperature && t <= self.af_temperature
        }
    }
    /// Transformation strain for a given martensite fraction.
    ///
    /// eps_trans = eps_max * f_m
    pub fn transformation_strain(&self, f_m: f64) -> f64 {
        self.eps_max * f_m.clamp(0.0, 1.0)
    }
}
/// 1-D Cahn-Hilliard concentration field with finite-difference update.
///
/// ∂c/∂t = M ∇²(μ)  where  μ = dF/dc - κ ∇²c
/// and F(c) is the Landau-style bulk free energy.
#[derive(Debug, Clone)]
pub struct CahnHilliardField {
    /// Concentration field (dimensionless, typically 0..1).
    pub c: Vec<f64>,
    /// Number of grid points.
    pub nx: usize,
}
impl CahnHilliardField {
    /// Create a new Cahn-Hilliard field with uniform initial concentration.
    pub fn new(nx: usize, c0: f64) -> Self {
        Self {
            c: vec![c0; nx],
            nx,
        }
    }
    /// Advance the field by one explicit time step using finite differences.
    ///
    /// Uses the standard 1-D explicit scheme with periodic boundary conditions.
    ///
    /// # Arguments
    /// * `dt`    - Time step (s).
    /// * `m`     - Mobility M.
    /// * `kappa` - Gradient energy coefficient κ.
    pub fn update(&mut self, dt: f64, m: f64, kappa: f64) {
        let nx = self.nx;
        let dx = 1.0;
        let c = self.c.clone();
        let mut mu = vec![0.0; nx];
        for i in 0..nx {
            let ip = if i + 1 < nx { i + 1 } else { 0 };
            let im = if i > 0 { i - 1 } else { nx - 1 };
            let laplacian_c = (c[ip] + c[im] - 2.0 * c[i]) / (dx * dx);
            let df_dc = 2.0 * c[i] * (1.0 - c[i]) * (1.0 - 2.0 * c[i]);
            mu[i] = df_dc - kappa * laplacian_c;
        }
        let mut c_new = self.c.clone();
        for i in 0..nx {
            let ip = if i + 1 < nx { i + 1 } else { 0 };
            let im = if i > 0 { i - 1 } else { nx - 1 };
            let laplacian_mu = (mu[ip] + mu[im] - 2.0 * mu[i]) / (dx * dx);
            c_new[i] = c[i] + dt * m * laplacian_mu;
        }
        self.c = c_new;
    }
    /// Total concentration (should be conserved).
    #[allow(dead_code)]
    pub fn total_concentration(&self) -> f64 {
        self.c.iter().sum()
    }
}
/// Simplified CALPHAD (Calculation of Phase Diagrams) model for a binary
/// A–B alloy using a Redlich-Kister excess interaction.
#[derive(Debug, Clone)]
pub struct CalphadModel {
    /// Molar Gibbs energy of pure component A (J/mol).
    pub g_a: f64,
    /// Molar Gibbs energy of pure component B (J/mol).
    pub g_b: f64,
    /// Zeroth Redlich-Kister interaction parameter L0 (J/mol).
    pub l0: f64,
    /// First Redlich-Kister interaction parameter L1 (J/mol).
    pub l1: f64,
}
impl CalphadModel {
    /// Create a new `CalphadModel`.
    pub fn new(g_a: f64, g_b: f64, l0: f64, l1: f64) -> Self {
        Self { g_a, g_b, l0, l1 }
    }
    /// Molar mixing enthalpy.
    ///
    /// H_mix = x_A · x_B · (L0 + L1·(x_A - x_B))
    pub fn mixing_enthalpy(&self, x_b: f64) -> f64 {
        let x_a = 1.0 - x_b;
        x_a * x_b * (self.l0 + self.l1 * (x_a - x_b))
    }
    /// Ideal configurational entropy contribution: -R·T·(x_A·ln x_A + x_B·ln x_B).
    pub fn ideal_entropy(&self, x_b: f64, temperature: f64) -> f64 {
        let x_a = 1.0 - x_b;
        let s_mix = if x_a > 0.0 && x_b > 0.0 {
            x_a * x_a.ln() + x_b * x_b.ln()
        } else {
            0.0
        };
        -R_GAS * temperature * s_mix
    }
    /// Molar Gibbs energy of mixing.
    ///
    /// G = x_A·G_A + x_B·G_B + H_mix - T·S_mix
    pub fn gibbs_energy(&self, x_b: f64, temperature: f64) -> f64 {
        let x_a = 1.0 - x_b;
        let g_ref = x_a * self.g_a + x_b * self.g_b;
        let h_mix = self.mixing_enthalpy(x_b);
        let ts_mix = self.ideal_entropy(x_b, temperature);
        g_ref + h_mix + ts_mix
    }
    /// Chemical potential of component A: μ_A = G - x_B · (dG/dx_B).
    ///
    /// Computed via a central finite difference on `gibbs_energy`.
    pub fn chemical_potential_a(&self, x_b: f64, temperature: f64) -> f64 {
        let h = 1.0e-7_f64;
        let x_b_hi = (x_b + h).min(1.0 - 1.0e-12);
        let x_b_lo = (x_b - h).max(1.0e-12);
        let dg_dx_b = (self.gibbs_energy(x_b_hi, temperature)
            - self.gibbs_energy(x_b_lo, temperature))
            / (x_b_hi - x_b_lo);
        self.gibbs_energy(x_b, temperature) - x_b * dg_dx_b
    }
}
/// Full austenite ↔ martensite phase fraction evolution with temperature
/// history tracking.  On each `step` call the model determines whether the
/// system is cooling (forming martensite) or heating (reverting to austenite)
/// and updates the phase fractions accordingly.
#[derive(Debug, Clone)]
pub struct AusteniteMartensiteKinetics {
    pub(super) km: KoistinenMarburger,
    /// Austenite-start temperature for reverse transformation (K).
    pub as_temp: f64,
    /// Austenite-finish temperature (K).
    pub af_temp: f64,
    /// Current martensite fraction \[0, 1\].
    pub f_martensite: f64,
    /// Last temperature (K) – used to detect cooling/heating.
    pub(super) last_temp: f64,
}
impl AusteniteMartensiteKinetics {
    /// Create a new kinetics model.
    pub fn new(ms: f64, k: f64, as_temp: f64, af_temp: f64, initial_temp: f64) -> Self {
        Self {
            km: KoistinenMarburger::new(ms, k),
            as_temp,
            af_temp,
            f_martensite: 0.0,
            last_temp: initial_temp,
        }
    }
    /// Advance the model to a new temperature `t`.
    ///
    /// On cooling below Ms, the KM equation is applied.
    /// On heating above As, martensite dissolves linearly.
    pub fn step(&mut self, t: f64) {
        let cooling = t < self.last_temp;
        self.last_temp = t;
        if cooling {
            let km_frac = self.km.fraction(t);
            if km_frac > self.f_martensite {
                self.f_martensite = km_frac;
            }
        } else if t >= self.af_temp {
            self.f_martensite = 0.0;
        } else if t > self.as_temp {
            let reverted = (t - self.as_temp) / (self.af_temp - self.as_temp);
            self.f_martensite = (self.f_martensite * (1.0 - reverted)).max(0.0);
        }
    }
    /// Current austenite fraction.
    pub fn f_austenite(&self) -> f64 {
        1.0 - self.f_martensite
    }
}
/// Type of solid-state phase transformation.
#[derive(Debug, Clone, PartialEq)]
pub enum PhaseTransformType {
    /// Diffusionless (displacive) martensitic transformation.
    Martensitic,
    /// Diffusion-controlled pearlitic transformation.
    Pearlitic,
    /// Intermediate bainitic transformation.
    Bainitic,
    /// Reversion to austenite (high-temperature parent phase).
    Austenitic,
}
/// Cahn-Hilliard spinodal decomposition model.
///
/// Determines the spinodal boundary and amplification factor for
/// concentration fluctuations in a binary alloy.
///
/// The free energy of mixing is approximated as:
/// f(c) = N*k_B*T * \[c*ln(c) + (1-c)*ln(1-c)\] + Ω*c*(1-c)
///
/// where Ω is the interaction parameter.
#[derive(Debug, Clone)]
pub struct SpinocalDecomposition {
    /// Interaction (mixing) parameter Ω (J/mol).
    pub omega: f64,
    /// Gradient energy coefficient κ (J·m²/mol).
    pub kappa: f64,
    /// Molar volume V_m (m³/mol).
    pub molar_volume: f64,
}
impl SpinocalDecomposition {
    /// Create a spinodal decomposition model.
    pub fn new(omega: f64, kappa: f64, molar_volume: f64) -> Self {
        Self {
            omega,
            kappa,
            molar_volume,
        }
    }
    /// Second derivative of free energy f''(c, T) = d²f/dc².
    ///
    /// f''(c,T) = R*T / (c*(1-c)) − 2*Ω
    pub fn f_double_prime(&self, c: f64, temperature: f64) -> f64 {
        let r = 8.314_462_618;
        if c <= 0.0 || c >= 1.0 {
            return f64::INFINITY;
        }
        r * temperature / (c * (1.0 - c)) - 2.0 * self.omega
    }
    /// Within the spinodal region: f''(c, T) < 0.
    ///
    /// Returns true if the composition c is inside the spinodal at temperature T.
    pub fn is_spinodal(&self, c: f64, temperature: f64) -> bool {
        self.f_double_prime(c, temperature) < 0.0
    }
    /// Spinodal boundary compositions at temperature T.
    ///
    /// At the spinodal: R*T / (c*(1-c)) = 2*Ω.
    /// Solving: c² - c + R*T/(2*Ω) = 0.
    ///
    /// Returns (c_low, c_high) or None if temperature is above spinodal limit.
    pub fn spinodal_boundary(&self, temperature: f64) -> Option<(f64, f64)> {
        let r = 8.314_462_618;
        if self.omega <= 0.0 {
            return None;
        }
        let discriminant = 0.25 - r * temperature / (2.0 * self.omega);
        if discriminant <= 0.0 {
            return None;
        }
        let sqrt_d = discriminant.sqrt();
        Some((0.5 - sqrt_d, 0.5 + sqrt_d))
    }
    /// Critical temperature T_c above which no spinodal decomposition occurs.
    ///
    /// T_c = Ω / (2 * R)
    pub fn critical_temperature(&self) -> f64 {
        let r = 8.314_462_618;
        self.omega / (2.0 * r)
    }
    /// Amplification factor R(q) for a perturbation with wavenumber q.
    ///
    /// R(q) = -M * (f''(c,T) * q² + 2*κ*q^4)
    ///
    /// where M is the mobility. Decomposition occurs where R(q) > 0.
    pub fn amplification_factor(&self, q: f64, c: f64, temperature: f64, mobility: f64) -> f64 {
        let f2 = self.f_double_prime(c, temperature);
        -mobility * (f2 * q * q + 2.0 * self.kappa * q.powi(4))
    }
    /// Fastest growing wavelength λ* at which R(q) is maximized.
    ///
    /// dR/dq² = 0 → q*² = -f''/(4κ)  → λ* = 2π/q*
    pub fn fastest_growing_wavelength(&self, c: f64, temperature: f64) -> Option<f64> {
        let f2 = self.f_double_prime(c, temperature);
        if f2 >= 0.0 || self.kappa <= 0.0 {
            return None;
        }
        let q_star_sq = -f2 / (4.0 * self.kappa);
        if q_star_sq <= 0.0 {
            return None;
        }
        let q_star = q_star_sq.sqrt();
        Some(2.0 * std::f64::consts::PI / q_star)
    }
}
/// Landau free-energy model extended with a φ^6 term for first-order
/// (discontinuous) transitions.
///
/// F(eta, T) = a*(T - Tc)/2 * eta^2 + b/4 * eta^4 + c/6 * eta^6
#[derive(Debug, Clone)]
pub struct FirstOrderTransition {
    /// Coefficient of the phi^2 term (analogous to a0 above).
    pub a: f64,
    /// Coefficient of the phi^4 term (may be negative for first-order).
    pub b: f64,
    /// Coefficient of the phi^6 term (must be > 0 for stability).
    pub c: f64,
    /// Spinodal temperature (limit of metastability of disordered phase).
    pub tc: f64,
    /// First-order transition temperature (free-energy degeneracy).
    pub t0: f64,
}
impl FirstOrderTransition {
    /// Create a new `FirstOrderTransition` model.
    pub fn new(a: f64, b: f64, c: f64, tc: f64, t0: f64) -> Self {
        Self { a, b, c, tc, t0 }
    }
    /// Free energy at order parameter `eta` and temperature `t`.
    pub fn free_energy(&self, eta: f64, t: f64) -> f64 {
        let a_t = self.a * (t - self.tc);
        a_t / 2.0 * eta * eta + self.b / 4.0 * eta.powi(4) + self.c / 6.0 * eta.powi(6)
    }
    /// Find local minima (eta values where dF/deta = 0) using Newton-Raphson.
    ///
    /// dF/deta = a*(T-Tc)*eta + b*eta^3 + c*eta^5
    ///         = eta * (a*(T-Tc) + b*eta^2 + c*eta^4)
    ///
    /// Non-trivial roots satisfy: c*u^2 + b*u + a*(T-Tc) = 0, where u = eta^2.
    pub fn local_minima(&self, t: f64) -> Vec<f64> {
        let mut minima = vec![0.0f64];
        let alpha = self.a * (t - self.tc);
        let discriminant = self.b * self.b - 4.0 * self.c * alpha;
        if discriminant >= 0.0 {
            let sqrt_d = discriminant.sqrt();
            let u1 = (-self.b + sqrt_d) / (2.0 * self.c);
            let u2 = (-self.b - sqrt_d) / (2.0 * self.c);
            for u in [u1, u2] {
                if u > 0.0 {
                    let eta = u.sqrt();
                    minima.push(eta);
                    minima.push(-eta);
                }
            }
        }
        minima.retain(|&eta| {
            let d2 = self.a * (t - self.tc) + 3.0 * self.b * eta * eta + 5.0 * self.c * eta.powi(4);
            d2 > 0.0
        });
        minima
    }
    /// Approximate latent heat at the first-order transition temperature t0.
    ///
    /// Uses the discontinuity in entropy: L ≈ T0 * a * eta_eq² where eta_eq
    /// is the non-trivial equilibrium order parameter at t0.
    pub fn latent_heat_approx(&self) -> f64 {
        let alpha = self.a * (self.t0 - self.tc);
        let discriminant = self.b * self.b - 4.0 * self.c * alpha;
        if discriminant < 0.0 {
            return 0.0;
        }
        let sqrt_d = discriminant.sqrt();
        let u = (-self.b + sqrt_d) / (2.0 * self.c);
        if u <= 0.0 {
            return 0.0;
        }
        self.t0 * self.a * u
    }
}
/// Landau free-energy model for a second-order (continuous) phase transition.
///
/// F(eta, T) = a0*(T - Tc)/2 * eta^2 + b/4 * eta^4
#[derive(Debug, Clone)]
pub struct LandauFreeEnergy {
    /// Temperature-dependent coefficient at critical point: a = a0*(T-Tc).
    pub a: f64,
    /// Quartic (stabilising) coefficient (must be > 0).
    pub b: f64,
    /// Critical temperature (K).
    pub tc: f64,
    /// Coefficient for the linear temperature dependence of `a`.
    pub a0: f64,
}
impl LandauFreeEnergy {
    /// Create a new `LandauFreeEnergy` model.
    pub fn new(a0: f64, b: f64, tc: f64) -> Self {
        Self { a: 0.0, b, tc, a0 }
    }
    /// Landau free energy at order parameter `eta` and `temperature`.
    pub fn free_energy(&self, eta: f64, temperature: f64) -> f64 {
        let a_t = self.a0 * (temperature - self.tc);
        a_t / 2.0 * eta * eta + self.b / 4.0 * eta.powi(4)
    }
    /// Equilibrium order parameter at `temperature`.
    ///
    /// Returns 0 above Tc, sqrt(-a/b) below Tc.
    pub fn equilibrium_order_parameter(&self, temperature: f64) -> f64 {
        if temperature >= self.tc {
            0.0
        } else {
            let a_t = self.a0 * (temperature - self.tc);
            (-a_t / self.b).sqrt()
        }
    }
    /// Heat-capacity jump at Tc: ΔCp = a0² / (2·b).
    pub fn heat_capacity_jump(&self) -> f64 {
        self.a0 * self.a0 / (2.0 * self.b)
    }
    /// Latent heat – 0 for a second-order (continuous) transition.
    pub fn latent_heat(&self) -> f64 {
        0.0
    }
}
/// Extended JMAK (Avrami) model with impingement correction.
///
/// f(t) = 1 - exp(-b * t^n * impingement_factor)
#[derive(Debug, Clone)]
pub struct JmakExtended {
    /// Pre-exponential rate constant b (s⁻ⁿ).
    pub b: f64,
    /// Avrami exponent n.
    pub n: f64,
    /// Activation energy Q (J/mol).
    pub activation_energy: f64,
    /// Impingement correction (1.0 = random nucleation, <1 for site saturation).
    pub impingement: f64,
}
impl JmakExtended {
    /// Create an extended JMAK model.
    pub fn new(b: f64, n: f64, activation_energy: f64, impingement: f64) -> Self {
        Self {
            b,
            n,
            activation_energy,
            impingement,
        }
    }
    /// Temperature-dependent rate constant.
    pub fn rate_at_temperature(&self, temperature: f64) -> f64 {
        let r = 8.314_462_618;
        self.b * (-self.activation_energy / (r * temperature)).exp()
    }
    /// Transformed fraction at time t and temperature T.
    pub fn fraction(&self, t: f64, temperature: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let k = self.rate_at_temperature(temperature);
        (1.0 - (-k * t.powf(self.n) * self.impingement).exp()).clamp(0.0, 1.0)
    }
    /// Time to reach a target fraction at given temperature.
    ///
    /// From f = 1 - exp(-k * t^n): t = (-ln(1-f) / (k * impingement))^(1/n)
    pub fn time_to_fraction(&self, fraction: f64, temperature: f64) -> f64 {
        if fraction <= 0.0 {
            return 0.0;
        }
        if fraction >= 1.0 {
            return f64::INFINITY;
        }
        let k = self.rate_at_temperature(temperature);
        if k <= 0.0 {
            return f64::INFINITY;
        }
        let arg = (-(1.0 - fraction).ln() / (k * self.impingement)).max(0.0);
        arg.powf(1.0 / self.n)
    }
    /// Transformation rate df/dt at time t and temperature T.
    pub fn transformation_rate(&self, t: f64, temperature: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        let k = self.rate_at_temperature(temperature);
        let f = self.fraction(t, temperature);
        self.n * k * t.powf(self.n - 1.0) * self.impingement * (1.0 - f)
    }
    /// Isothermal TTT start time (1% transformation).
    pub fn ttt_start_time(&self, temperature: f64) -> f64 {
        self.time_to_fraction(0.01, temperature)
    }
    /// Isothermal TTT finish time (99% transformation).
    pub fn ttt_finish_time(&self, temperature: f64) -> f64 {
        self.time_to_fraction(0.99, temperature)
    }
}
/// Eutectic transformation: simultaneous solidification of two phases.
///
/// For a binary eutectic alloy, the eutectic composition solidifies
/// at a constant temperature (eutectic temperature) as a lamellar or
/// rod microstructure.
///
/// The interlamellar spacing is given by the Jackson-Hunt theory:
/// `λ² * v = K` (constant for a given alloy)
/// where `λ` is lamellar spacing \[m\] and `v` is growth velocity \[m/s\].
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct EutecticTransformation {
    /// Eutectic temperature \[K\]
    pub t_eutectic: f64,
    /// Eutectic composition (wt% or mol fraction of component B)
    pub c_eutectic: f64,
    /// Composition of alpha phase at eutectic \[wt%\]
    pub c_alpha: f64,
    /// Composition of beta phase at eutectic \[wt%\]
    pub c_beta: f64,
    /// Jackson-Hunt constant K = λ² * v \[m³/s\]
    pub k_jh: f64,
}
impl EutecticTransformation {
    /// Create a new eutectic transformation model.
    pub fn new(t_eutectic: f64, c_eutectic: f64, c_alpha: f64, c_beta: f64, k_jh: f64) -> Self {
        Self {
            t_eutectic,
            c_eutectic,
            c_alpha,
            c_beta,
            k_jh,
        }
    }
    /// Al-Si eutectic (12.6 wt% Si at 577°C = 850 K).
    pub fn al_si() -> Self {
        Self::new(850.0, 12.6, 1.5, 99.0, 1.0e-18)
    }
    /// Volume fraction of alpha phase from lever rule at eutectic.
    ///
    /// `f_α = (C_β - C_eut) / (C_β - C_α)`
    pub fn alpha_fraction(&self) -> f64 {
        let denom = self.c_beta - self.c_alpha;
        if denom.abs() < 1e-30 {
            return 0.5;
        }
        ((self.c_beta - self.c_eutectic) / denom).clamp(0.0, 1.0)
    }
    /// Volume fraction of beta phase at eutectic.
    pub fn beta_fraction(&self) -> f64 {
        1.0 - self.alpha_fraction()
    }
    /// Equilibrium interlamellar spacing \[m\] at growth velocity v \[m/s\].
    ///
    /// `λ = sqrt(K / v)`
    pub fn lamellar_spacing(&self, velocity_m_s: f64) -> f64 {
        if velocity_m_s < f64::EPSILON {
            return f64::INFINITY;
        }
        (self.k_jh / velocity_m_s).sqrt()
    }
    /// Growth velocity \[m/s\] for a given interlamellar spacing λ \[m\].
    ///
    /// `v = K / λ²`
    pub fn growth_velocity(&self, spacing_m: f64) -> f64 {
        if spacing_m < f64::EPSILON {
            return f64::INFINITY;
        }
        self.k_jh / (spacing_m * spacing_m)
    }
    /// Undercooling required to nucleate eutectic \[K\] (Gibbs-Thomson estimate).
    ///
    /// `ΔT = 2 * Γ / λ`  where Γ = γ * T_eut / ΔH_f is the capillarity parameter.
    ///
    /// `Γ` ≈ 1.0e-7 K·m for metals (default).
    pub fn nucleation_undercooling(&self, spacing_m: f64) -> f64 {
        let capillarity = 1.0e-7_f64;
        if spacing_m < f64::EPSILON {
            return f64::INFINITY;
        }
        2.0 * capillarity / spacing_m
    }
    /// Check if a given composition is at (or near) the eutectic.
    pub fn is_eutectic_composition(&self, c: f64, tolerance: f64) -> bool {
        (c - self.c_eutectic).abs() < tolerance
    }
}
/// 2-D phase-field (Allen-Cahn) model for solidification.
///
/// dphi/dt = L · \[W²·∇²phi - f'(phi) + λ·Q/T·(1-phi²)²\]
///
/// where f(phi) = phi²·(1-phi)² is the double-well potential and
/// f'(phi) = 2·phi·(1-phi)·(1-2·phi).
#[derive(Debug, Clone)]
pub struct PhaseFieldModel2D {
    /// Grid size in x direction.
    pub nx: usize,
    /// Grid size in y direction.
    pub ny: usize,
    /// Order parameter field (phi=1 solid, phi=0 liquid).
    pub phi: Vec<f64>,
    /// Temperature field (K).
    pub temperature: Vec<f64>,
    /// Interface width W (m).
    pub interface_width: f64,
    /// Interface mobility L (m² s⁻¹ J⁻¹).
    pub mobility: f64,
    /// Latent heat Q (J m⁻³).
    pub latent_heat: f64,
}
impl PhaseFieldModel2D {
    /// Create a new flat liquid domain at uniform temperature 0 K.
    pub fn new(nx: usize, ny: usize, w: f64, l: f64, q: f64) -> Self {
        Self {
            nx,
            ny,
            phi: vec![0.0; nx * ny],
            temperature: vec![0.0; nx * ny],
            interface_width: w,
            mobility: l,
            latent_heat: q,
        }
    }
    /// Initialise a spherical solid nucleus centred at `(cx, cy)` with the
    /// given `radius` (in grid cells).
    pub fn initialize_seed(&mut self, cx: f64, cy: f64, radius: f64) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let dx = i as f64 - cx;
                let dy = j as f64 - cy;
                let r = (dx * dx + dy * dy).sqrt();
                let phi = 0.5 * (1.0 - ((r - radius) / (self.interface_width + 1.0e-12)).tanh());
                self.phi[j * self.nx + i] = phi.clamp(0.0, 1.0);
            }
        }
    }
    /// Advance the phase field by one time step `dt` with grid spacing `dx`.
    pub fn step(&mut self, dt: f64, dx: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let w2 = self.interface_width * self.interface_width;
        let l = self.mobility;
        let q = self.latent_heat;
        let lambda = 1.0;
        let mut phi_new = self.phi.clone();
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                let phi_c = self.phi[idx];
                let t_c = self.temperature[idx].max(1.0);
                let ip = if i + 1 < nx { i + 1 } else { 0 };
                let im = if i > 0 { i - 1 } else { nx - 1 };
                let jp = if j + 1 < ny { j + 1 } else { 0 };
                let jm = if j > 0 { j - 1 } else { ny - 1 };
                let lap = (self.phi[j * nx + ip]
                    + self.phi[j * nx + im]
                    + self.phi[jp * nx + i]
                    + self.phi[jm * nx + i]
                    - 4.0 * phi_c)
                    / (dx * dx);
                let df = 2.0 * phi_c * (1.0 - phi_c) * (1.0 - 2.0 * phi_c);
                let coupling = lambda * q / t_c * (1.0 - phi_c * phi_c).powi(2);
                let dphi = l * (w2 * lap - df + coupling);
                phi_new[idx] = (phi_c + dt * dphi).clamp(0.0, 1.0);
            }
        }
        self.phi = phi_new;
    }
    /// Volume fraction of solid (mean of phi field).
    pub fn solid_fraction(&self) -> f64 {
        let n = (self.nx * self.ny) as f64;
        self.phi.iter().sum::<f64>() / n
    }
    /// Approximate interface length: Σ |∇phi| · dx.
    pub fn interface_length(&self, dx: f64) -> f64 {
        let nx = self.nx;
        let ny = self.ny;
        let mut total = 0.0_f64;
        for j in 0..ny {
            for i in 0..nx {
                let ip = if i + 1 < nx { i + 1 } else { 0 };
                let jp = if j + 1 < ny { j + 1 } else { 0 };
                let grad_x = (self.phi[j * nx + ip] - self.phi[j * nx + i]) / dx;
                let grad_y = (self.phi[jp * nx + i] - self.phi[j * nx + i]) / dx;
                total += (grad_x * grad_x + grad_y * grad_y).sqrt();
            }
        }
        total * dx
    }
}
/// Latent heat release model tracking enthalpy change during phase transformation.
///
/// H_latent = L * Δf  where Δf is the increment in transformed fraction.
#[derive(Debug, Clone)]
pub struct LatentHeatModel {
    /// Latent heat per unit mass (J/kg).
    pub latent_heat_per_mass: f64,
    /// Density of material (kg/m³).
    pub density: f64,
}
impl LatentHeatModel {
    /// Create a latent heat model.
    pub fn new(latent_heat_per_mass: f64, density: f64) -> Self {
        Self {
            latent_heat_per_mass,
            density,
        }
    }
    /// Steel solidification: L ≈ 271 kJ/kg, ρ ≈ 7800 kg/m³.
    pub fn steel_solidification() -> Self {
        Self::new(271_000.0, 7800.0)
    }
    /// Aluminium solidification: L ≈ 397 kJ/kg, ρ ≈ 2700 kg/m³.
    pub fn aluminum_solidification() -> Self {
        Self::new(397_000.0, 2700.0)
    }
    /// Volumetric latent heat (J/m³).
    pub fn volumetric_latent_heat(&self) -> f64 {
        self.latent_heat_per_mass * self.density
    }
    /// Enthalpy release for a fractional transformation df (J/m³).
    pub fn enthalpy_release(&self, df: f64) -> f64 {
        self.volumetric_latent_heat() * df.abs()
    }
    /// Temperature rise due to latent heat release in an adiabatic system.
    ///
    /// ΔT = L * df / Cp_volumetric
    ///
    /// where Cp_volumetric = Cp_specific * density.
    pub fn adiabatic_temperature_rise(&self, df: f64, cp_specific: f64) -> f64 {
        if cp_specific <= 0.0 {
            return 0.0;
        }
        self.latent_heat_per_mass * df.abs() / cp_specific
    }
    /// Effective latent heat fraction released between two fractions f1 and f2.
    pub fn fraction_released(&self, f1: f64, f2: f64) -> f64 {
        (f2 - f1).clamp(-1.0, 1.0).abs()
    }
}
