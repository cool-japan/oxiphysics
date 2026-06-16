//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use crate::constitutive::{ConstitutiveModel, ConstitutiveResponse, ViscoelasticState};

/// Kelvin-Voigt chain (series of KV elements).
///
/// Creep compliance J(t) = sum_i (1/E_i) * (1 - exp(-t*E_i/eta_i))
#[derive(Debug, Clone)]
pub struct KelvinVoigtChain {
    /// Individual KV elements.
    pub elements: Vec<KelvinVoigtModel>,
}
impl KelvinVoigtChain {
    /// Create a new KV chain.
    pub fn new() -> Self {
        Self {
            elements: Vec::new(),
        }
    }
}
impl Default for KelvinVoigtChain {
    fn default() -> Self {
        Self::new()
    }
}
impl KelvinVoigtChain {
    /// Add a Kelvin-Voigt element with elastic modulus `e` and viscosity `eta`.
    pub fn add_element(&mut self, e: f64, eta: f64) {
        self.elements.push(KelvinVoigtModel::new(e, eta));
    }
    /// Creep compliance J(t) = sum_i J_i(t).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        self.elements.iter().map(|kv| kv.creep_compliance(t)).sum()
    }
    /// Number of elements in the chain.
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }
    /// Retardation times (tau_i = eta_i / E_i) for all elements.
    pub fn retardation_times(&self) -> Vec<f64> {
        self.elements
            .iter()
            .map(|kv| kv.viscosity / kv.elastic_modulus)
            .collect()
    }
}
/// Fractional Standard Linear Solid (FSLS).
///
/// Replaces the dashpot in the classical SLS with a springpot:
///
///   E*(ω) = E∞ + E₁(iωτ)^α / (1 + (iωτ)^α)
///
/// where τ = (C/E₁)^(1/α) is the characteristic relaxation time.
#[derive(Debug, Clone)]
pub struct FractionalSls {
    /// Equilibrium modulus E∞ (Pa).
    pub e_inf: f64,
    /// Non-equilibrium spring modulus E₁ (Pa).
    pub e1: f64,
    /// Springpot element.
    pub springpot: SpringPot,
}
impl FractionalSls {
    /// Create a new Fractional SLS.
    pub fn new(e_inf: f64, e1: f64, springpot: SpringPot) -> Self {
        Self {
            e_inf,
            e1,
            springpot,
        }
    }
    /// Characteristic relaxation time τ = (C / E₁)^(1/α).
    pub fn relaxation_time(&self) -> f64 {
        let alpha = self.springpot.alpha;
        (self.springpot.quasi_property / self.e1).powf(1.0 / alpha)
    }
    /// Storage modulus E'(ω).
    ///
    /// E'(ω) = E∞ + E₁ (ωτ)^α \[(ωτ)^α + cos(απ/2)\] / \[1 + 2(ωτ)^α cos(απ/2) + (ωτ)^(2α)\]
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let pi = std::f64::consts::PI;
        let tau = self.relaxation_time();
        let wt = (omega * tau).powf(self.springpot.alpha);
        let cos_ap2 = (self.springpot.alpha * pi / 2.0).cos();
        let denom = 1.0 + 2.0 * wt * cos_ap2 + wt * wt;
        self.e_inf + self.e1 * wt * (wt + cos_ap2) / denom
    }
    /// Loss modulus E''(ω).
    ///
    /// E''(ω) = E₁ (ωτ)^α sin(απ/2) / \[1 + 2(ωτ)^α cos(απ/2) + (ωτ)^(2α)\]
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let pi = std::f64::consts::PI;
        let tau = self.relaxation_time();
        let wt = (omega * tau).powf(self.springpot.alpha);
        let cos_ap2 = (self.springpot.alpha * pi / 2.0).cos();
        let sin_ap2 = (self.springpot.alpha * pi / 2.0).sin();
        let denom = 1.0 + 2.0 * wt * cos_ap2 + wt * wt;
        self.e1 * wt * sin_ap2 / denom
    }
    /// Loss tangent tan δ = E'' / E'.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let ep = self.storage_modulus(omega);
        let epp = self.loss_modulus(omega);
        if ep.abs() < 1e-30 { 0.0 } else { epp / ep }
    }
    /// Glassy (short-time) modulus E₀ = E∞ + E₁.
    pub fn glassy_modulus(&self) -> f64 {
        self.e_inf + self.e1
    }
}
/// Prony series representation for relaxation modulus.
///
/// E(t) = E_inf + sum_i g_i * exp(-t / tau_i)
#[derive(Debug, Clone)]
pub struct PronySeries {
    /// Equilibrium modulus E_inf.
    pub e_inf: f64,
    /// Prony coefficients (g_i, tau_i).
    pub terms: Vec<(f64, f64)>,
}
impl PronySeries {
    /// Create a new Prony series.
    pub fn new(e_inf: f64) -> Self {
        Self {
            e_inf,
            terms: Vec::new(),
        }
    }
    /// Add a term (g_i, tau_i).
    pub fn add_term(&mut self, g: f64, tau: f64) {
        self.terms.push((g, tau));
    }
    /// Relaxation modulus E(t).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let sum: f64 = self
            .terms
            .iter()
            .map(|&(g, tau)| g * (-t / tau).exp())
            .sum();
        self.e_inf + sum
    }
    /// Storage modulus E'(omega).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let sum: f64 = self
            .terms
            .iter()
            .map(|&(g, tau)| {
                let wt = omega * tau;
                g * wt * wt / (1.0 + wt * wt)
            })
            .sum();
        self.e_inf + sum
    }
    /// Loss modulus E''(omega).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.terms
            .iter()
            .map(|&(g, tau)| {
                let wt = omega * tau;
                g * wt / (1.0 + wt * wt)
            })
            .sum()
    }
    /// Convert to a GeneralizedMaxwell model (each term becomes a Maxwell element with weight 1).
    pub fn to_generalized_maxwell(&self) -> GeneralizedMaxwell {
        let mut gm = GeneralizedMaxwell::new(self.e_inf);
        for &(g, tau) in &self.terms {
            gm.add_element(MaxwellModel::new(g, g * tau), 1.0);
        }
        gm
    }
    /// Compute the relaxation modulus G(t) from the Prony series.
    ///
    /// The shear relaxation modulus (often denoted G(t) for rubbery / soft materials)
    /// is expressed as a Prony series of decaying exponentials:
    ///
    /// G(t) = G_inf + Σᵢ gᵢ · exp(-t / τᵢ)
    ///
    /// This is equivalent to `relaxation_modulus` but the naming emphasises the
    /// *shear* relaxation context used in polymer viscoelasticity.
    ///
    /// # Arguments
    /// * `t` - Time \[s\] (t ≥ 0)
    ///
    /// # Returns
    /// Shear relaxation modulus G(t) \[Pa\]
    pub fn compute_relaxation_modulus(&self, t: f64) -> f64 {
        let sum: f64 = self
            .terms
            .iter()
            .map(|&(g, tau)| g * (-t / tau).exp())
            .sum();
        self.e_inf + sum
    }
}
/// Scott-Blair (springpot) element: intermediate between spring and dashpot.
///
/// The springpot constitutive law uses the fractional Caputo derivative:
///
///   σ(t) = C · D^α ε(t)
///
/// where α ∈ \[0, 1\] interpolates between a purely elastic spring (α=0)
/// and a purely viscous dashpot (α=1).
///
/// For sinusoidal excitation ε = ε₀ exp(iωt) the complex modulus is:
///
///   E*(ω) = C · (iω)^α
///
/// giving storage and loss moduli:
///
///   E'(ω) = C ω^α cos(α π/2)
///   E''(ω) = C ω^α sin(α π/2)
#[derive(Debug, Clone)]
pub struct SpringPot {
    /// Quasi-property C (Pa·s^α).
    pub quasi_property: f64,
    /// Fractional order α ∈ \[0, 1\].
    pub alpha: f64,
}
impl SpringPot {
    /// Create a new Scott-Blair springpot.
    ///
    /// # Panics
    /// Panics (debug only) if `alpha` is outside `[0, 1]`.
    pub fn new(quasi_property: f64, alpha: f64) -> Self {
        debug_assert!((0.0..=1.0).contains(&alpha), "alpha must be in [0, 1]");
        Self {
            quasi_property,
            alpha,
        }
    }
    /// Storage modulus E'(ω) = C ω^α cos(απ/2).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let pi = std::f64::consts::PI;
        self.quasi_property * omega.powf(self.alpha) * (self.alpha * pi / 2.0).cos()
    }
    /// Loss modulus E''(ω) = C ω^α sin(απ/2).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let pi = std::f64::consts::PI;
        self.quasi_property * omega.powf(self.alpha) * (self.alpha * pi / 2.0).sin()
    }
    /// Loss tangent tan δ = E'' / E' = tan(απ/2).
    ///
    /// Note: independent of frequency — a hallmark of fractional models.
    pub fn loss_tangent(&self) -> f64 {
        let pi = std::f64::consts::PI;
        (self.alpha * pi / 2.0).tan()
    }
    /// Complex modulus magnitude |E*(ω)| = C ω^α.
    pub fn complex_modulus_magnitude(&self, omega: f64) -> f64 {
        self.quasi_property * omega.powf(self.alpha)
    }
    /// Phase angle δ = α π/2 (radians).
    pub fn phase_angle(&self) -> f64 {
        self.alpha * std::f64::consts::FRAC_PI_2
    }
}
/// Generalized Maxwell model (Prony series).
///
/// E(t) = E_inf + sum w_i * E_i * exp(-t / tau_i)
#[derive(Debug, Clone)]
pub struct GeneralizedMaxwell {
    /// Long-time (equilibrium) spring modulus E_inf (Pa)
    pub spring_inf: f64,
    /// Individual Maxwell elements
    pub maxwell_elements: Vec<MaxwellModel>,
    /// Weights for each Maxwell element
    pub weights: Vec<f64>,
}
impl GeneralizedMaxwell {
    /// Create a new Generalized Maxwell model with equilibrium modulus `spring_inf`.
    pub fn new(spring_inf: f64) -> Self {
        Self {
            spring_inf,
            maxwell_elements: Vec::new(),
            weights: Vec::new(),
        }
    }
    /// Add a Maxwell element with the given weight.
    pub fn add_element(&mut self, model: MaxwellModel, weight: f64) {
        self.maxwell_elements.push(model);
        self.weights.push(weight);
    }
    /// Relaxation modulus E(t) = E_inf + sum w_i * E_i * exp(-t / tau_i).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let sum: f64 = self
            .maxwell_elements
            .iter()
            .zip(self.weights.iter())
            .map(|(m, &w)| w * m.relaxation_modulus(t))
            .sum();
        self.spring_inf + sum
    }
    /// Return Prony coefficients as (E_i * w_i, tau_i) pairs.
    pub fn prony_coefficients(&self) -> Vec<(f64, f64)> {
        self.maxwell_elements
            .iter()
            .zip(self.weights.iter())
            .map(|(m, &w)| (m.elastic_modulus * w, m.relaxation_time()))
            .collect()
    }
    /// Storage modulus E'(omega) = E_inf + sum g_i * (omega*tau_i)^2 / (1 + (omega*tau_i)^2).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let sum: f64 = self
            .maxwell_elements
            .iter()
            .zip(self.weights.iter())
            .map(|(m, &w)| {
                let tau = m.relaxation_time();
                let wt = omega * tau;
                w * m.elastic_modulus * wt * wt / (1.0 + wt * wt)
            })
            .sum();
        self.spring_inf + sum
    }
    /// Loss modulus E''(omega) = sum g_i * omega*tau_i / (1 + (omega*tau_i)^2).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.maxwell_elements
            .iter()
            .zip(self.weights.iter())
            .map(|(m, &w)| {
                let tau = m.relaxation_time();
                let wt = omega * tau;
                w * m.elastic_modulus * wt / (1.0 + wt * wt)
            })
            .sum()
    }
    /// Loss tangent tan(delta) = E''(omega) / E'(omega).
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let ep = self.storage_modulus(omega);
        let epp = self.loss_modulus(omega);
        if ep.abs() < 1e-30 { 0.0 } else { epp / ep }
    }
    /// Compute tan(δ) — the material loss tangent at angular frequency ω.
    ///
    /// tan(δ) = G''(ω) / G'(ω)
    ///
    /// This is the ratio of energy dissipated per cycle to the peak stored
    /// energy, and is a direct measure of the material damping capacity.
    /// Values of tan(δ) >> 1 indicate a predominantly viscous response;
    /// values << 1 indicate predominantly elastic behaviour.
    ///
    /// # Arguments
    /// * `omega` - Angular frequency ω \[rad/s\]
    ///
    /// # Returns
    /// Loss tangent tan(δ) \[dimensionless, ≥ 0\]
    pub fn compute_tan_delta(&self, omega: f64) -> f64 {
        let g_prime = self.storage_modulus(omega);
        let g_double_prime = self.loss_modulus(omega);
        if g_prime.abs() < 1e-30 {
            f64::INFINITY
        } else {
            g_double_prime / g_prime
        }
    }
    /// Short-time (glassy) modulus E0 = E_inf + sum w_i * E_i.
    pub fn glassy_modulus(&self) -> f64 {
        let sum: f64 = self
            .maxwell_elements
            .iter()
            .zip(self.weights.iter())
            .map(|(m, &w)| w * m.elastic_modulus)
            .sum();
        self.spring_inf + sum
    }
    /// Number of Maxwell elements.
    pub fn n_elements(&self) -> usize {
        self.maxwell_elements.len()
    }
}
/// WLF (Williams-Landel-Ferry) shift factor.
///
/// log10(a_T) = -C1 * (T - T_ref) / (C2 + (T - T_ref))
#[derive(Debug, Clone)]
pub struct WlfShift {
    /// WLF constant C1.
    pub c1: f64,
    /// WLF constant C2 (K or degC).
    pub c2: f64,
    /// Reference temperature T_ref (K or degC).
    pub t_ref: f64,
}
impl WlfShift {
    /// Create a new WLF shift factor model.
    pub fn new(c1: f64, c2: f64, t_ref: f64) -> Self {
        Self { c1, c2, t_ref }
    }
    /// Compute log10(a_T) at temperature T.
    pub fn log_shift_factor(&self, t: f64) -> f64 {
        let dt = t - self.t_ref;
        -self.c1 * dt / (self.c2 + dt)
    }
    /// Compute the shift factor a_T at temperature T.
    pub fn shift_factor(&self, t: f64) -> f64 {
        10.0_f64.powf(self.log_shift_factor(t))
    }
    /// Shift a frequency by the temperature factor: omega_reduced = omega * a_T.
    pub fn shifted_frequency(&self, omega: f64, t: f64) -> f64 {
        omega * self.shift_factor(t)
    }
}
/// Maxwell element (spring + dashpot in series) -- new API.
///
/// Uses field names `elastic_modulus` and `viscosity` to distinguish from the
/// legacy [`Maxwell`] struct.
#[derive(Debug, Clone)]
pub struct MaxwellModel {
    /// Young's modulus E (Pa)
    pub elastic_modulus: f64,
    /// Dynamic viscosity eta (Pa*s)
    pub viscosity: f64,
}
impl MaxwellModel {
    /// Create a new Maxwell model with elastic modulus `e` and viscosity `eta`.
    pub fn new(e: f64, eta: f64) -> Self {
        Self {
            elastic_modulus: e,
            viscosity: eta,
        }
    }
    /// Relaxation time tau = eta / E.
    pub fn relaxation_time(&self) -> f64 {
        self.viscosity / self.elastic_modulus
    }
    /// Relaxation modulus E(t) = E * exp(-t / tau).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let tau = self.relaxation_time();
        self.elastic_modulus * (-t / tau).exp()
    }
    /// Creep compliance J(t) = 1/E + t/eta.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        1.0 / self.elastic_modulus + t / self.viscosity
    }
    /// Stress update using exact integration of the Maxwell ODE.
    ///
    /// Given old stress `sigma_old`, strain rate `epsilon_dot`, and time step `dt`:
    ///
    /// sigma_new = sigma_old * exp(-dt/tau) + E * epsilon_dot * tau * (1 - exp(-dt/tau))
    pub fn stress_update(&self, _epsilon: f64, epsilon_dot: f64, sigma_old: f64, dt: f64) -> f64 {
        let tau = self.relaxation_time();
        let exp_term = (-dt / tau).exp();
        sigma_old * exp_term + self.elastic_modulus * epsilon_dot * tau * (1.0 - exp_term)
    }
    /// Storage modulus E'(omega) = E * (omega*tau)^2 / (1 + (omega*tau)^2).
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let tau = self.relaxation_time();
        let wt = omega * tau;
        self.elastic_modulus * wt * wt / (1.0 + wt * wt)
    }
    /// Loss modulus E''(omega) = E * omega*tau / (1 + (omega*tau)^2).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let tau = self.relaxation_time();
        let wt = omega * tau;
        self.elastic_modulus * wt / (1.0 + wt * wt)
    }
    /// Loss tangent tan(delta) = E'' / E' = 1 / (omega * tau).
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        let tau = self.relaxation_time();
        1.0 / (omega * tau)
    }
    /// Complex modulus magnitude |E*| = sqrt(E'^2 + E''^2).
    pub fn complex_modulus_magnitude(&self, omega: f64) -> f64 {
        let ep = self.storage_modulus(omega);
        let epp = self.loss_modulus(omega);
        (ep * ep + epp * epp).sqrt()
    }
    /// Compute the dynamic (complex) modulus components (G', G'') for a Maxwell element.
    ///
    /// For a single Maxwell element in shear:
    ///
    /// G'(ω)  = G · (ωτ)² / (1 + (ωτ)²)   \[storage / elastic part\]
    /// G''(ω) = G · ωτ    / (1 + (ωτ)²)   \[loss / viscous part\]
    ///
    /// # Arguments
    /// * `omega` - Angular frequency ω \[rad/s\]
    ///
    /// # Returns
    /// `(G_prime, G_double_prime)` — storage and loss moduli \[Pa\]
    pub fn compute_dynamic_modulus(&self, omega: f64) -> (f64, f64) {
        let tau = self.relaxation_time();
        let wt = omega * tau;
        let denom = 1.0 + wt * wt;
        let g_prime = self.elastic_modulus * wt * wt / denom;
        let g_double_prime = self.elastic_modulus * wt / denom;
        (g_prime, g_double_prime)
    }
}
/// Power-law creep compliance model.
///
/// J(t) = J₀ + J₁ t^n
///
/// where J₀ = 1/E_g is the instantaneous (glassy) compliance, J₁ is the
/// creep coefficient, and n ∈ (0, 1] is the power-law exponent.
///
/// This form is often used for polymers and metals at elevated temperature.
#[derive(Debug, Clone)]
pub struct PowerLawCreep {
    /// Instantaneous compliance J₀ = 1/E_g (1/Pa).
    pub j0: f64,
    /// Power-law coefficient J₁ (1/(Pa·s^n)).
    pub j1: f64,
    /// Power-law exponent n ∈ (0, 1].
    pub n: f64,
}
impl PowerLawCreep {
    /// Create a new power-law creep model.
    pub fn new(j0: f64, j1: f64, n: f64) -> Self {
        Self { j0, j1, n }
    }
    /// Creep compliance J(t) = J₀ + J₁ t^n.
    pub fn creep_compliance(&self, t: f64) -> f64 {
        self.j0 + self.j1 * t.powf(self.n)
    }
    /// Creep strain under constant stress σ₀.
    pub fn creep_strain(&self, sigma0: f64, t: f64) -> f64 {
        sigma0 * self.creep_compliance(t)
    }
    /// Creep rate d J/dt = J₁ n t^(n-1).
    pub fn creep_rate(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return if self.n < 1.0 {
                f64::INFINITY
            } else {
                self.j1 * self.n
            };
        }
        self.j1 * self.n * t.powf(self.n - 1.0)
    }
    /// Retardation spectrum approximation (Alfrey approximation).
    ///
    /// L(τ) ≈ -J₁ n (τ^n) / (d ln J / d ln t) ≈ J₁ n τ^(n-1) / τ
    ///      = J₁ n τ^(n-1)  (valid for n << 1)
    ///
    /// Returns L at retardation time τ.
    pub fn retardation_spectrum(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return 0.0;
        }
        self.j1 * self.n * tau.powf(self.n - 1.0)
    }
}
/// Kelvin-Voigt element (spring + dashpot in parallel) -- new API.
///
/// Uses field names `elastic_modulus` and `viscosity`.
#[derive(Debug, Clone)]
pub struct KelvinVoigtModel {
    /// Young's modulus E (Pa)
    pub elastic_modulus: f64,
    /// Dynamic viscosity eta (Pa*s)
    pub viscosity: f64,
}
impl KelvinVoigtModel {
    /// Create a new Kelvin-Voigt model with elastic modulus `e` and viscosity `eta`.
    pub fn new(e: f64, eta: f64) -> Self {
        Self {
            elastic_modulus: e,
            viscosity: eta,
        }
    }
    /// Creep compliance J(t) = (1/E) * (1 - exp(-t*E/eta)).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        (1.0 / self.elastic_modulus) * (1.0 - (-t * self.elastic_modulus / self.viscosity).exp())
    }
    /// Relaxation modulus for KV model (t > 0).
    ///
    /// Strictly E(t) = E + eta*delta(t); for t > 0 this returns E.
    pub fn relaxation_modulus(&self, _t: f64) -> f64 {
        self.elastic_modulus
    }
    /// Strain update using exact integration of the KV ODE under constant stress.
    ///
    /// epsilon_new = epsilon_old * exp(-E*dt/eta) + (sigma/E) * (1 - exp(-E*dt/eta))
    pub fn strain_update(&self, sigma: f64, epsilon_old: f64, dt: f64) -> f64 {
        let exp_term = (-self.elastic_modulus * dt / self.viscosity).exp();
        epsilon_old * exp_term + (sigma / self.elastic_modulus) * (1.0 - exp_term)
    }
    /// Storage modulus E'(omega) = E.
    pub fn storage_modulus(&self, _omega: f64) -> f64 {
        self.elastic_modulus
    }
    /// Loss modulus E''(omega) = eta * omega.
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        self.viscosity * omega
    }
    /// Loss tangent tan(delta) = eta*omega / E.
    pub fn loss_tangent(&self, omega: f64) -> f64 {
        self.viscosity * omega / self.elastic_modulus
    }
}
/// Kelvin-Voigt viscoelastic model.
///
/// sigma = E*epsilon + eta*epsilon_dot  (elastic + dashpot in parallel)
#[derive(Debug, Clone)]
pub struct KelvinVoigt {
    /// Young's modulus E (Pa)
    pub young_modulus: f64,
    /// Dynamic viscosity eta (Pa*s)
    pub viscosity: f64,
}
impl KelvinVoigt {
    /// Create a new Kelvin-Voigt model.
    pub fn new(young_modulus: f64, viscosity: f64) -> Self {
        Self {
            young_modulus,
            viscosity,
        }
    }
    /// Compute stress given strain and strain rate.
    ///
    /// sigma = E*epsilon + eta*epsilon_dot
    pub fn stress(&self, strain: f64, strain_rate: f64) -> f64 {
        self.young_modulus * strain + self.viscosity * strain_rate
    }
    /// Compute creep strain as a function of time under constant stress sigma0.
    ///
    /// epsilon(t) = (sigma0/E) * (1 - exp(-E*t/eta))
    pub fn creep_strain(&self, sigma0: f64, t: f64) -> f64 {
        let tau = self.relaxation_time();
        (sigma0 / self.young_modulus) * (1.0 - (-t / tau).exp())
    }
    /// Time constant tau = eta/E.
    pub fn relaxation_time(&self) -> f64 {
        self.viscosity / self.young_modulus
    }
    /// Compute the creep compliance J(t) for the Kelvin-Voigt model.
    ///
    /// Under a step stress σ₀ applied at t = 0, the creep compliance is:
    ///
    /// J(t) = (1/E) · (1 − exp(−t/τ))
    ///
    /// where τ = η/E is the retardation time. The compliance represents the
    /// strain per unit applied stress at time t.
    ///
    /// Note: unlike the Maxwell model, J(t) saturates to 1/E rather than
    /// growing without bound, reflecting the restoring elastic spring.
    ///
    /// # Arguments
    /// * `t` - Time elapsed since load application \[s\] (t ≥ 0)
    ///
    /// # Returns
    /// Creep compliance J(t) \[1/Pa\]
    pub fn compute_creep_compliance(&self, t: f64) -> f64 {
        let tau = self.relaxation_time();
        (1.0 / self.young_modulus) * (1.0 - (-t / tau).exp())
    }
}
/// Arrhenius shift factor.
///
/// ln(a_T) = (E_a / R) * (1/T - 1/T_ref)
#[derive(Debug, Clone)]
pub struct ArrheniusShift {
    /// Activation energy E_a (J/mol).
    pub activation_energy: f64,
    /// Gas constant R = 8.314 J/(mol*K).
    pub r: f64,
    /// Reference temperature T_ref (K).
    pub t_ref: f64,
}
impl ArrheniusShift {
    /// Create a new Arrhenius shift model.
    pub fn new(activation_energy: f64, t_ref: f64) -> Self {
        Self {
            activation_energy,
            r: 8.314,
            t_ref,
        }
    }
    /// Compute ln(a_T) at temperature T.
    pub fn ln_shift_factor(&self, t: f64) -> f64 {
        (self.activation_energy / self.r) * (1.0 / t - 1.0 / self.t_ref)
    }
    /// Compute the shift factor a_T at temperature T (must be in Kelvin).
    pub fn shift_factor(&self, t: f64) -> f64 {
        self.ln_shift_factor(t).exp()
    }
    /// Shift a frequency by the temperature factor.
    pub fn shifted_frequency(&self, omega: f64, t: f64) -> f64 {
        omega * self.shift_factor(t)
    }
}
/// Nonlinear Kelvin-Voigt model with power-law dashpot.
///
/// Constitutive law: σ = E·ε + η |dε/dt|^(m-1) · dε/dt
///
/// For m=1 this reduces to the classical (linear) Kelvin-Voigt model.
/// For m≠1 the dashpot is nonlinear (strain-rate hardening/softening).
#[derive(Debug, Clone)]
pub struct NonlinearKelvinVoigt {
    /// Elastic stiffness E (Pa).
    pub young_modulus: f64,
    /// Viscosity coefficient η (Pa·s^m).
    pub viscosity: f64,
    /// Rate sensitivity exponent m (m=1 → linear).
    pub m: f64,
}
impl NonlinearKelvinVoigt {
    /// Create a new nonlinear KV model.
    pub fn new(young_modulus: f64, viscosity: f64, m: f64) -> Self {
        Self {
            young_modulus,
            viscosity,
            m,
        }
    }
    /// Compute stress.
    ///
    /// σ = E·ε + η · sign(dε/dt) · |dε/dt|^m
    pub fn stress(&self, strain: f64, strain_rate: f64) -> f64 {
        let sign = if strain_rate >= 0.0 { 1.0 } else { -1.0 };
        self.young_modulus * strain + self.viscosity * sign * strain_rate.abs().powf(self.m)
    }
    /// Linear Kelvin-Voigt relaxation time τ = η / E (valid for m=1 only).
    pub fn linear_relaxation_time(&self) -> f64 {
        self.viscosity / self.young_modulus
    }
    /// Effective dynamic modulus at frequency ω and amplitude ε₀.
    ///
    /// For a sinusoidal input ε = ε₀ sin(ωt) the apparent viscous stress
    /// amplitude is proportional to (ω ε₀)^m so:
    ///
    ///   |E*(ω, ε₀)| ≈ sqrt(E² + (η (ω ε₀)^(m-1) ω)²)
    ///
    /// This reduces to the linear result for m=1.
    pub fn effective_dynamic_modulus(&self, omega: f64, epsilon0: f64) -> f64 {
        let viscous_term = self.viscosity * (omega * epsilon0).powf(self.m - 1.0) * omega;
        (self.young_modulus.powi(2) + viscous_term.powi(2)).sqrt()
    }
}
/// Burgers model: KV element in series with Maxwell element.
///
/// Total strain: ε = ε_spring + ε_KV + ε_dashpot
///
/// Under constant stress σ₀:
///
///   ε(t) = σ₀/E_M + (σ₀/E_KV)(1 - exp(-E_KV t / η_KV)) + σ₀ t / η_M
///
/// where the three terms are instantaneous elastic, delayed elastic
/// (retarded), and steady-state viscous flow.
#[derive(Debug, Clone)]
pub struct Burgers {
    /// Maxwell spring modulus E_M (Pa).
    pub e_maxwell: f64,
    /// Maxwell dashpot viscosity η_M (Pa·s).
    pub eta_maxwell: f64,
    /// Kelvin-Voigt spring modulus E_KV (Pa).
    pub e_kv: f64,
    /// Kelvin-Voigt dashpot viscosity η_KV (Pa·s).
    pub eta_kv: f64,
}
impl Burgers {
    /// Create a new Burgers model.
    pub fn new(e_maxwell: f64, eta_maxwell: f64, e_kv: f64, eta_kv: f64) -> Self {
        Self {
            e_maxwell,
            eta_maxwell,
            e_kv,
            eta_kv,
        }
    }
    /// Retardation time of the KV element τ_KV = η_KV / E_KV.
    pub fn retardation_time(&self) -> f64 {
        self.eta_kv / self.e_kv
    }
    /// Creep compliance J(t) under constant stress (three-term expression).
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let tau_kv = self.retardation_time();
        1.0 / self.e_maxwell
            + (1.0 / self.e_kv) * (1.0 - (-t / tau_kv).exp())
            + t / self.eta_maxwell
    }
    /// Creep strain under constant stress σ₀.
    pub fn creep_strain(&self, sigma0: f64, t: f64) -> f64 {
        sigma0 * self.creep_compliance(t)
    }
    /// Instantaneous compliance J(0) = 1 / E_M.
    pub fn instantaneous_compliance(&self) -> f64 {
        1.0 / self.e_maxwell
    }
    /// Long-time creep rate dε/dt → σ₀ / η_M (Newtonian viscous flow).
    pub fn steady_state_creep_rate(&self, sigma0: f64) -> f64 {
        sigma0 / self.eta_maxwell
    }
    /// Storage modulus E'(ω) (linearised, valid for small deformations).
    ///
    /// For the Burgers model in the frequency domain:
    ///
    ///   1/E*(ω) = 1/(iω η_M) + 1/\[E_M + iω η_M_spring\] + 1/\[E_KV + iω η_KV\]
    ///
    /// Here we use the simpler series combination of a Maxwell and a KV element.
    pub fn storage_modulus(&self, omega: f64) -> f64 {
        let tau_m = self.eta_maxwell / self.e_maxwell;
        let wt_m = omega * tau_m;
        let ep_maxwell = self.e_maxwell * wt_m * wt_m / (1.0 + wt_m * wt_m);
        let ep_kv = self.e_kv;
        if ep_maxwell < 1.0 {
            return ep_kv;
        }
        1.0 / (1.0 / ep_maxwell + 1.0 / ep_kv)
    }
    /// Loss modulus E''(ω).
    pub fn loss_modulus(&self, omega: f64) -> f64 {
        let tau_m = self.eta_maxwell / self.e_maxwell;
        let wt_m = omega * tau_m;
        let epp_maxwell = self.e_maxwell * wt_m / (1.0 + wt_m * wt_m);
        let epp_kv = self.eta_kv * omega;
        epp_maxwell + epp_kv
    }
}
/// Kohlrausch-Williams-Watts (KWW) stretched-exponential relaxation.
///
/// Relaxation modulus:
///
///   E(t) = (E₀ - E∞) · exp(-(t/τ)^β) + E∞
///
/// where β ∈ (0, 1] is the stretching exponent and τ is the characteristic
/// relaxation time.  β=1 recovers the classical Maxwell form.
///
/// This model is widely used for glassy and amorphous materials.
#[derive(Debug, Clone)]
pub struct Kww {
    /// Glassy modulus E₀ (Pa).
    pub e0: f64,
    /// Equilibrium modulus E∞ (Pa).
    pub e_inf: f64,
    /// Characteristic relaxation time τ (s).
    pub tau: f64,
    /// Stretching exponent β ∈ (0, 1].
    pub beta: f64,
}
impl Kww {
    /// Create a new KWW model.
    pub fn new(e0: f64, e_inf: f64, tau: f64, beta: f64) -> Self {
        Self {
            e0,
            e_inf,
            tau,
            beta,
        }
    }
    /// Relaxation modulus E(t) = (E₀-E∞)·exp(-(t/τ)^β) + E∞.
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let phi = (-(t / self.tau).powf(self.beta)).exp();
        (self.e0 - self.e_inf) * phi + self.e_inf
    }
    /// Mean relaxation time ⟨τ⟩ = (τ/β) Γ(1/β).
    ///
    /// Uses Stirling/Lanczos approximation of Γ.
    pub fn mean_relaxation_time(&self) -> f64 {
        self.tau / self.beta * gamma_approx(1.0 / self.beta)
    }
    /// Normalised relaxation function φ(t) ∈ \[0, 1\].
    pub fn phi(&self, t: f64) -> f64 {
        (-(t / self.tau).powf(self.beta)).exp()
    }
    /// Half-relaxation time t_{1/2} where φ(t) = 0.5.
    ///
    /// t_{1/2} = τ (ln 2)^(1/β)
    pub fn half_relaxation_time(&self) -> f64 {
        self.tau * (2.0_f64.ln()).powf(1.0 / self.beta)
    }
}
/// Maxwell viscoelastic model.
///
/// depsilon/dt = (1/E)*dsigma/dt + sigma/eta  (elastic + dashpot in series)
#[derive(Debug, Clone)]
pub struct Maxwell {
    /// Young's modulus E (Pa)
    pub young_modulus: f64,
    /// Dynamic viscosity eta (Pa*s)
    pub viscosity: f64,
}
impl Maxwell {
    /// Create a new Maxwell model.
    pub fn new(young_modulus: f64, viscosity: f64) -> Self {
        Self {
            young_modulus,
            viscosity,
        }
    }
    /// Time constant tau = eta/E.
    pub fn relaxation_time(&self) -> f64 {
        self.viscosity / self.young_modulus
    }
    /// Stress relaxation under constant strain epsilon0.
    ///
    /// sigma(t) = E * epsilon0 * exp(-t/tau)
    pub fn stress_relaxation(&self, epsilon0: f64, t: f64) -> f64 {
        let tau = self.relaxation_time();
        self.young_modulus * epsilon0 * (-t / tau).exp()
    }
}
/// Standard Linear Solid (Zener model): springs E1, E2 and dashpot eta.
#[derive(Debug, Clone)]
pub struct StandardLinearSolid {
    /// Equilibrium spring modulus (Pa)
    pub e1: f64,
    /// Non-equilibrium spring modulus (Pa)
    pub e2: f64,
    /// Dashpot viscosity eta (Pa*s)
    pub eta: f64,
}
impl StandardLinearSolid {
    /// Create a new Standard Linear Solid model.
    pub fn new(e1: f64, e2: f64, eta: f64) -> Self {
        Self { e1, e2, eta }
    }
    /// Time constant tau = eta/E2.
    pub fn relaxation_time(&self) -> f64 {
        self.eta / self.e2
    }
    /// Long-time (equilibrium) modulus E_inf = E1.
    pub fn long_time_modulus(&self) -> f64 {
        self.e1
    }
    /// Short-time (glassy) modulus E0 = E1 + E2.
    pub fn short_time_modulus(&self) -> f64 {
        self.e1 + self.e2
    }
    /// Relaxation modulus E(t) = E1 + E2 * exp(-t/tau).
    pub fn relaxation_modulus(&self, t: f64) -> f64 {
        let tau = self.relaxation_time();
        self.e1 + self.e2 * (-t / tau).exp()
    }
    /// Creep compliance J(t) for the SLS model.
    ///
    /// J(t) = 1/E1 - (E2 / (E1*(E1+E2))) * exp(-E1*t / (eta*(1 + E1/E2)))
    /// Simplified form: J(t) = 1/(E1+E2) + (E2/(E1*(E1+E2))) * (1 - exp(-t/tau_c))
    /// where tau_c = eta * (E1 + E2) / (E1 * E2)
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let e_sum = self.e1 + self.e2;
        let tau_c = self.eta * e_sum / (self.e1 * self.e2);
        1.0 / e_sum + (self.e2 / (self.e1 * e_sum)) * (1.0 - (-t / tau_c).exp())
    }
}
/// Time-Temperature Superposition master curve builder.
///
/// Given a set of isothermal relaxation modulus curves and a shift-factor
/// model, this struct accumulates (reduced_time, modulus) pairs to form
/// a master curve at the reference temperature.
#[derive(Debug, Clone)]
pub struct MasterCurveBuilder {
    /// Reference temperature T_ref (K or °C — must be consistent).
    pub t_ref: f64,
    /// Accumulated (log10 t_reduced, modulus) data points.
    pub data: Vec<(f64, f64)>,
}
impl MasterCurveBuilder {
    /// Create a new empty master-curve builder.
    pub fn new(t_ref: f64) -> Self {
        Self {
            t_ref,
            data: Vec::new(),
        }
    }
    /// Add an isothermal data point measured at temperature `t_meas`,
    /// time `t_meas_time` (s), and relaxation modulus `e_t` (Pa).
    ///
    /// Uses a WLF shift factor to map to reduced time.
    pub fn add_point_wlf(&mut self, wlf: &WlfShift, t_meas: f64, t_meas_time: f64, e_t: f64) {
        let log_at = wlf.log_shift_factor(t_meas);
        let log_t_reduced = t_meas_time.log10() + log_at;
        self.data.push((log_t_reduced, e_t));
    }
    /// Add an isothermal data point using an Arrhenius shift factor.
    pub fn add_point_arrhenius(
        &mut self,
        arr: &ArrheniusShift,
        t_meas: f64,
        t_meas_time: f64,
        e_t: f64,
    ) {
        let ln_at = arr.ln_shift_factor(t_meas);
        let log_at = ln_at / std::f64::consts::LN_10;
        let log_t_reduced = t_meas_time.log10() + log_at;
        self.data.push((log_t_reduced, e_t));
    }
    /// Sort the master-curve data by reduced time.
    pub fn sort(&mut self) {
        self.data
            .sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }
    /// Linearly interpolate the master-curve modulus at log10(t_reduced) = log_t.
    pub fn interpolate(&self, log_t: f64) -> Option<f64> {
        if self.data.len() < 2 {
            return None;
        }
        for w in self.data.windows(2) {
            let (x0, y0) = w[0];
            let (x1, y1) = w[1];
            if log_t >= x0 && log_t <= x1 {
                let frac = (log_t - x0) / (x1 - x0);
                return Some(y0 + frac * (y1 - y0));
            }
        }
        None
    }
    /// Number of data points accumulated.
    pub fn n_points(&self) -> usize {
        self.data.len()
    }
}
/// Continuous retardation spectrum L(τ) modelled as a log-normal distribution.
///
/// Creep compliance: J(t) = J₀ + ∫ L(τ)/τ · (1 - exp(-t/τ)) d ln τ
///
/// For the log-normal spectrum:
///
///   L(τ) = J_lm / (σ √(2π)) · exp(-(ln(τ/τ_c))² / (2σ²))
///
/// The integral is evaluated numerically over a wide range of τ.
#[derive(Debug, Clone)]
pub struct LogNormalRetardation {
    /// Instantaneous compliance J₀ (1/Pa).
    pub j0: f64,
    /// Total creep compliance amplitude J_lm (1/Pa).
    pub j_lm: f64,
    /// Central retardation time τ_c (s).
    pub tau_c: f64,
    /// Log-width σ (dimensionless, > 0).
    pub sigma: f64,
}
impl LogNormalRetardation {
    /// Create a new log-normal retardation model.
    pub fn new(j0: f64, j_lm: f64, tau_c: f64, sigma: f64) -> Self {
        Self {
            j0,
            j_lm,
            tau_c,
            sigma,
        }
    }
    /// L(τ): log-normal retardation spectrum evaluated at τ.
    pub fn spectrum(&self, tau: f64) -> f64 {
        if tau <= 0.0 {
            return 0.0;
        }
        let x = (tau / self.tau_c).ln() / self.sigma;
        let pi = std::f64::consts::PI;
        self.j_lm / (self.sigma * (2.0 * pi).sqrt()) * (-0.5 * x * x).exp()
    }
    /// Creep compliance J(t) evaluated by Gauss-Legendre quadrature (20 pts)
    /// over ln τ ∈ \[ln τ_c - 5σ, ln τ_c + 5σ\].
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let ln_lo = self.tau_c.ln() - 5.0 * self.sigma;
        let ln_hi = self.tau_c.ln() + 5.0 * self.sigma;
        let n = 200usize;
        let d_ln = (ln_hi - ln_lo) / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let ln_tau = ln_lo + (i as f64 + 0.5) * d_ln;
                let tau = ln_tau.exp();
                let l_tau = self.spectrum(tau);
                l_tau * (1.0 - (-t / tau).exp()) * d_ln
            })
            .sum();
        self.j0 + integral
    }
    /// Storage compliance J'(ω) = J₀ + ∫ L(τ)/(1 + (ωτ)²) d ln τ.
    pub fn storage_compliance(&self, omega: f64) -> f64 {
        let ln_lo = self.tau_c.ln() - 5.0 * self.sigma;
        let ln_hi = self.tau_c.ln() + 5.0 * self.sigma;
        let n = 200usize;
        let d_ln = (ln_hi - ln_lo) / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let ln_tau = ln_lo + (i as f64 + 0.5) * d_ln;
                let tau = ln_tau.exp();
                let l_tau = self.spectrum(tau);
                let wt2 = (omega * tau).powi(2);
                l_tau / (1.0 + wt2) * d_ln
            })
            .sum();
        self.j0 + integral
    }
    /// Loss compliance J''(ω) = ∫ L(τ) ωτ/(1 + (ωτ)²) d ln τ.
    pub fn loss_compliance(&self, omega: f64) -> f64 {
        let ln_lo = self.tau_c.ln() - 5.0 * self.sigma;
        let ln_hi = self.tau_c.ln() + 5.0 * self.sigma;
        let n = 200usize;
        let d_ln = (ln_hi - ln_lo) / n as f64;
        (0..n)
            .map(|i| {
                let ln_tau = ln_lo + (i as f64 + 0.5) * d_ln;
                let tau = ln_tau.exp();
                let l_tau = self.spectrum(tau);
                let wt = omega * tau;
                l_tau * wt / (1.0 + wt * wt) * d_ln
            })
            .sum()
    }
}

impl GeneralizedMaxwell {
    /// Tensorial Voigt-6 viscoelastic stress update (component-wise scalar
    /// relaxation modulus, exact exponential integration per Maxwell branch).
    ///
    /// Each branch partial stress relaxes from its previous value toward the
    /// instantaneous branch stress `E_a · ε` over the step via the discrete
    /// recurrence `h_new = exp(−dt/τ_a)·h_old + E_a·(1−exp(−dt/τ_a))·ε`. The
    /// total stress is `σ = E_inf·ε + Σ_a h_new` and the consistent tangent is
    /// the diagonal `D·I₆` with `D = E_inf + Σ_a E_a·(1−exp(−dt/τ_a))`.
    ///
    /// * `strain` — total strain at t_{n+1} (Voigt-6).
    /// * `branch_stress` — per-branch history (partial) stresses at t_n; a
    ///   shorter/empty slice is treated as all-zero virgin history.
    /// * `dt` — time-step size (s).
    ///
    /// Returns `(stress, tangent, new_branch_stress)`.
    pub fn voigt_stress_update(
        &self,
        strain: &[f64; 6],
        branch_stress: &[[f64; 6]],
        dt: f64,
    ) -> ([f64; 6], [[f64; 6]; 6], Vec<[f64; 6]>) {
        let n = self.n_elements();
        let mut new_hist: Vec<[f64; 6]> = Vec::with_capacity(n);
        let mut d_scalar = self.spring_inf;
        let mut branch_sum = [0.0_f64; 6];
        for (a, (m, &w)) in self
            .maxwell_elements
            .iter()
            .zip(self.weights.iter())
            .enumerate()
        {
            let e_a = w * m.elastic_modulus;
            let tau = m.relaxation_time();
            let exp_a = if tau > 0.0 { (-dt / tau).exp() } else { 0.0 };
            let old = branch_stress.get(a).copied().unwrap_or([0.0; 6]);
            let mut h_new = [0.0_f64; 6];
            for (k, (hk, (&ok, &sk))) in h_new
                .iter_mut()
                .zip(old.iter().zip(strain.iter()))
                .enumerate()
            {
                *hk = exp_a * ok + e_a * (1.0 - exp_a) * sk;
                branch_sum[k] += *hk;
            }
            d_scalar += e_a * (1.0 - exp_a);
            new_hist.push(h_new);
        }
        let mut stress = [0.0_f64; 6];
        for (sk, (&bk, &ek)) in stress
            .iter_mut()
            .zip(branch_sum.iter().zip(strain.iter()))
        {
            *sk = self.spring_inf * ek + bk;
        }
        let mut tangent = [[0.0_f64; 6]; 6];
        for (k, row) in tangent.iter_mut().enumerate() {
            row[k] = d_scalar;
        }
        (stress, tangent, new_hist)
    }
}

impl ConstitutiveModel for GeneralizedMaxwell {
    type State = ViscoelasticState;
    fn stress_update(
        &self,
        strain: &[f64; 6],
        state: &ViscoelasticState,
        dt: f64,
    ) -> ConstitutiveResponse<ViscoelasticState> {
        let (stress, tangent, new_hist) =
            self.voigt_stress_update(strain, &state.branch_stress, dt);
        ConstitutiveResponse {
            stress,
            tangent,
            state: ViscoelasticState {
                branch_stress: new_hist,
            },
        }
    }
    fn n_state_vars(&self) -> usize {
        6 * self.n_elements()
    }
}
