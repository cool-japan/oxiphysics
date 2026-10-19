//! Brinson one-dimensional constitutive model for shape memory alloys (SMA).
//!
//! The martensite volume fraction is split into a stress-induced (detwinned)
//! part `ξ_S` and a temperature-induced (twinned) part `ξ_T`, `ξ = ξ_S + ξ_T`
//! (Brinson, *J. Intell. Mater. Syst. Struct.* 4 (1993) 229-242).  With the
//! rule-of-mixtures modulus `E(ξ) = E_A + ξ (E_M − E_A)` the constitutive law is
//!
//! ```text
//! σ = E(ξ)·(ε − ε_L·ξ_S) + Θ·(T − T₀)
//! ```
//!
//! so only the stress-induced martensite produces transformation strain.  The
//! fractions evolve by cosine kinetics in three transformation windows:
//!
//! * **forward** (loading, A → detwinned M and twinned M → detwinned M): the
//!   stress window `σ_s^cr + C_M·(T − M_s)⁺ … σ_f^cr + C_M·(T − M_s)⁺`,
//!   `ξ_S = 1 − (1 − ξ_S0)·(1 + cos πs)/2`, `ξ_T = ξ_T0·(1 + cos πs)/2`;
//! * **reverse** (unloading or heating, M → A): the window of
//!   `T − σ/C_A` between `A_s` and `A_f`,
//!   `ξ = ξ₀·(1 + cos πs)/2` (both `ξ_S` and `ξ_T` scale with it);
//! * **cooling** (A → twinned M): the window `M_s … M_f`, the remaining
//!   austenite becomes twinned martensite.
//!
//! `s ∈ [0, 1]` is the progress through the window, and `ξ_S0, ξ_T0, ξ₀` are the
//! fractions at the start of the transformation.  The tensile part of the
//! stress, `max(σ, 0)`, drives the stress windows.
//!
//! ## Transformation memory
//!
//! Every window remembers the largest progress reached since the fractions last
//! changed, and a branch transforms only when the progress exceeds that
//! frontier: the untransformed material then decays by the survival ratio
//! `(1 + cos π s_new) / (1 + cos π s_front)`.  Because the cosine ramp
//! factorises, applying this ratio reproduces Brinson's closed form from the
//! start of the transformation exactly and independently of the step size, while
//! elastic unloading and reloading inside a window retrace the same path
//! without jumps.  Whenever any branch transforms, the frontiers of all branches
//! restart from the current state.
//!
//! The stress is not known before the kinetics are evaluated (the windows are
//! stress windows) and the strain is the driving variable, so [`BrinsonModel::update`]
//! solves the scalar equation `σ = E(ξ(σ))·(ε − ε_L ξ_S(σ)) + Θ (T − T₀)` by
//! bisection.  A plain fixed-point iteration is not an option: in the
//! transformation window the map has the slope `E ε_L dξ_S/dσ`, which is of the
//! order of 100 for NiTi.  The bracket `±(E_max (|ε| + ε_L) + |Θ (T − T₀)|)` always
//! contains a root, because the right-hand side of the constitutive law is
//! bounded by it for any fractions in `[0, 1]`, so the residual `σ − (…)` changes
//! sign across it; bisection needs nothing more.  For physical parameters the
//! residual also increases with σ (the transformation only softens the
//! response), which makes the root unique.

use std::f64::consts::PI;

use crate::{
    Error, Result, require_finite, require_in_range, require_non_negative, require_positive,
};

/// Default stress-influence coefficient `C_M` (Pa/K): typical NiTi.
const DEFAULT_C_M: f64 = 8.0e6;
/// Default stress-influence coefficient `C_A` (Pa/K): typical NiTi.
const DEFAULT_C_A: f64 = 13.0e6;
/// Default critical stress for the start of detwinning (Pa): typical NiTi.
const DEFAULT_SIGMA_S_CR: f64 = 100.0e6;
/// Default critical stress for the finish of detwinning (Pa): typical NiTi.
const DEFAULT_SIGMA_F_CR: f64 = 170.0e6;

/// Relative width of the stress bracket at which the bisection stops: about
/// four units in the last place, which is always reachable (the bracket cannot
/// get narrower than one ulp).  A tight bracket keeps the returned stress
/// `E(ξ)(ε − ε_L ξ_S) + Θ(T − T₀)` consistent with the root even in a very
/// steep transformation window.
const SOLVE_TOLERANCE: f64 = 1.0e-15;
/// Upper limit of bisection steps (the bracket halves every step, ~85 are needed).
const MAX_BISECTIONS: usize = 200;

/// Fraction of the transformable material that is still untransformed after a
/// cosine ramp has advanced to progress `s ∈ [0, 1]` of its window:
/// `(1 + cos πs)/2 = sin²(π(1 − s)/2)`.  The sine form is exact at both ends
/// (1 at `s = 0`, 0 at `s = 1`) and free of cancellation near `s = 1`.
fn cosine_remaining(s: f64) -> f64 {
    let x = 0.5 * PI * (1.0 - s);
    x.sin() * x.sin()
}

/// Fraction of the material that survives untransformed when a ramp moves from
/// progress `s_old` to `s_new ≥ s_old`.
fn survival_ratio(s_old: f64, s_new: f64) -> f64 {
    let remaining_old = cosine_remaining(s_old);
    if remaining_old <= f64::MIN_POSITIVE {
        return 0.0;
    }
    (cosine_remaining(s_new) / remaining_old).clamp(0.0, 1.0)
}

/// Position of `value` in the window `[start, finish]`, clamped to `[0, 1]`.  A
/// vanishing or inverted window acts as a step at `finish`.
fn window_progress(value: f64, start: f64, finish: f64) -> f64 {
    let width = finish - start;
    if width <= 0.0 {
        return if value >= finish { 1.0 } else { 0.0 };
    }
    ((value - start) / width).clamp(0.0, 1.0)
}

/// Fractions after a kinetic update, with the window progress it was evaluated at.
#[derive(Debug, Clone, Copy)]
struct Outcome {
    xi_s: f64,
    xi_t: f64,
    /// Progress `[forward, reverse, cooling]` at the evaluated state.
    progress: [f64; 3],
    /// Whether any branch transformed.
    advanced: bool,
}

/// Transformation memory: the fractions it belongs to and the frontier
/// `[forward, reverse, cooling]` of the three windows.
#[derive(Debug, Clone, Copy)]
struct Memory {
    xi_s: f64,
    xi_t: f64,
    frontier: [f64; 3],
}

/// Brinson shape memory alloy model with stress-induced and temperature-induced
/// martensite.
///
/// See the [module documentation](self) for the constitutive law and the
/// kinetics.  Strain-controlled loading is advanced with
/// [`BrinsonModel::update`]; a stress-free thermal history with
/// [`BrinsonModel::update_temperature`].
///
/// `ξ_S` and `ξ_T` start at zero (austenite).  The kinetic constants of
/// [`BrinsonModel::new`] default to typical NiTi values; set them for another
/// alloy with [`BrinsonModel::with_stress_kinetics`].
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
    /// Stress-influence coefficient of the martensitic transformation
    /// `C_M = dσ_cr/dT` (Pa/K).
    pub c_m: f64,
    /// Stress-influence coefficient of the austenitic transformation
    /// `C_A = dσ/dT` of the reverse windows (Pa/K); must be positive.
    pub c_a: f64,
    /// Critical stress at the start of the stress-induced transformation below
    /// `M_s` (Pa).
    pub sigma_s_cr: f64,
    /// Critical stress at the finish of the stress-induced transformation below
    /// `M_s` (Pa); larger than `sigma_s_cr`.
    pub sigma_f_cr: f64,
    /// Thermoelastic coefficient Θ of the constitutive law (Pa/K).
    pub theta: f64,
    /// Reference temperature `T₀` at which `Θ·(T − T₀)` vanishes (K).
    pub t0: f64,
    /// Temperature of the last update (K).
    temperature: f64,
    /// Stress of the last update (Pa).
    stress: f64,
    /// Transformation memory; `None` until the first update.
    memory: Option<Memory>,
}

impl BrinsonModel {
    /// Create a new Brinson SMA model in the austenitic state.
    ///
    /// The kinetic constants default to typical NiTi values (`C_M = 8 MPa/K`,
    /// `C_A = 13 MPa/K`, `σ_s^cr = 100 MPa`, `σ_f^cr = 170 MPa`); the
    /// thermoelastic coefficient `Θ` is zero and the reference temperature
    /// `T₀` is `A_f`.  The first update treats the material as having been
    /// cooled from `A_f` to the new temperature.
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
            c_m: DEFAULT_C_M,
            c_a: DEFAULT_C_A,
            sigma_s_cr: DEFAULT_SIGMA_S_CR,
            sigma_f_cr: DEFAULT_SIGMA_F_CR,
            theta: 0.0,
            t0: af_temp,
            temperature: af_temp,
            stress: 0.0,
            memory: None,
        }
    }

    /// Return the model with the stress-influence coefficients `c_m`, `c_a`
    /// (Pa/K) and the critical detwinning stresses `sigma_s_cr`, `sigma_f_cr`
    /// (Pa) set.  The values are checked by [`BrinsonModel::validate`].
    pub fn with_stress_kinetics(
        mut self,
        c_m: f64,
        c_a: f64,
        sigma_s_cr: f64,
        sigma_f_cr: f64,
    ) -> Self {
        self.c_m = c_m;
        self.c_a = c_a;
        self.sigma_s_cr = sigma_s_cr;
        self.sigma_f_cr = sigma_f_cr;
        self
    }

    /// Return the model with the thermoelastic coefficient `theta` (Pa/K) and
    /// the reference temperature `t0` (K) of `Θ·(T − T₀)` set.
    pub fn with_thermal_stress(mut self, theta: f64, t0: f64) -> Self {
        self.theta = theta;
        self.t0 = t0;
        self
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

    /// Temperature of the last update (K).
    pub fn temperature(&self) -> f64 {
        self.temperature
    }

    /// Stress of the last update (Pa); zero before the first update.
    pub fn current_stress(&self) -> f64 {
        self.stress
    }

    /// Stress window `(σ_start, σ_finish)` (Pa) of the forward (loading)
    /// transformation at temperature `t`: `σ_s^cr … σ_f^cr`, shifted up by
    /// `C_M·(T − M_s)` above `M_s`.
    pub fn forward_window(&self, t: f64) -> (f64, f64) {
        let shift = self.c_m * (t - self.ms).max(0.0);
        (self.sigma_s_cr + shift, self.sigma_f_cr + shift)
    }

    /// Stress window `(σ_start, σ_finish)` (Pa) of the reverse (unloading)
    /// transformation at temperature `t`: it starts when the stress falls to
    /// `C_A·(T − A_s)` and finishes at `C_A·(T − A_f)`.
    pub fn reverse_window(&self, t: f64) -> (f64, f64) {
        (self.c_a * (t - self.as_temp), self.c_a * (t - self.af_temp))
    }

    /// Check the parameters and the state; the error names the offending field.
    ///
    /// Required: positive moduli and `C_A`, non-negative `ε_L`, `C_M` and
    /// `σ_s^cr`, `M_s > M_f`, `A_f > A_s`, `σ_f^cr > σ_s^cr`, fractions in
    /// `[0, 1]` with `ξ_S + ξ_T ≤ 1`, and finite values throughout.
    pub fn validate(&self) -> Result<()> {
        let finite = [
            ("e_austenite", self.e_austenite),
            ("e_martensite", self.e_martensite),
            ("ms", self.ms),
            ("mf", self.mf),
            ("as_temp", self.as_temp),
            ("af_temp", self.af_temp),
            ("eps_l", self.eps_l),
            ("c_m", self.c_m),
            ("c_a", self.c_a),
            ("sigma_s_cr", self.sigma_s_cr),
            ("sigma_f_cr", self.sigma_f_cr),
            ("theta", self.theta),
            ("t0", self.t0),
            ("xi_s", self.xi_s),
            ("xi_t", self.xi_t),
        ];
        for (name, value) in finite {
            require_finite(name, value)?;
        }
        require_positive("e_austenite", self.e_austenite)?;
        require_positive("e_martensite", self.e_martensite)?;
        require_positive("c_a", self.c_a)?;
        require_non_negative("eps_l", self.eps_l)?;
        require_non_negative("c_m", self.c_m)?;
        require_non_negative("sigma_s_cr", self.sigma_s_cr)?;
        if self.sigma_f_cr <= self.sigma_s_cr {
            return Err(Error::below_minimum(
                "sigma_f_cr",
                self.sigma_f_cr,
                self.sigma_s_cr,
            ));
        }
        if self.ms <= self.mf {
            return Err(Error::below_minimum("ms", self.ms, self.mf));
        }
        if self.af_temp <= self.as_temp {
            return Err(Error::below_minimum("af_temp", self.af_temp, self.as_temp));
        }
        require_in_range("xi_s", self.xi_s, 0.0, 1.0)?;
        require_in_range("xi_t", self.xi_t, 0.0, 1.0)?;
        if self.xi_s + self.xi_t > 1.0 + 1.0e-12 {
            return Err(Error::above_maximum(
                "xi_s + xi_t",
                self.xi_s + self.xi_t,
                1.0,
            ));
        }
        Ok(())
    }

    /// Set the martensite fractions and the temperature, restarting the
    /// transformation memory at this stress-free state.
    ///
    /// Use it to start a simulation from a state other than fully austenitic
    /// at the reference temperature, e.g. twinned martensite
    /// (`xi_s = 0, xi_t = 1`) or austenite already at the test temperature.
    pub fn set_state(&mut self, xi_s: f64, xi_t: f64, temperature: f64) -> Result<()> {
        require_in_range("xi_s", xi_s, 0.0, 1.0)?;
        require_in_range("xi_t", xi_t, 0.0, 1.0)?;
        if xi_s + xi_t > 1.0 + 1.0e-12 {
            return Err(Error::above_maximum("xi_s + xi_t", xi_s + xi_t, 1.0));
        }
        require_finite("temperature", temperature)?;
        self.xi_s = xi_s;
        self.xi_t = xi_t;
        self.temperature = temperature;
        self.stress = 0.0;
        self.memory = Some(Memory {
            xi_s,
            xi_t,
            frontier: self.progress(0.0, temperature),
        });
        Ok(())
    }

    /// Progress of the forward transformation through its stress window.
    fn forward_progress(&self, sigma: f64, t: f64) -> f64 {
        let (start, finish) = self.forward_window(t);
        window_progress(sigma.max(0.0), start, finish)
    }

    /// Progress of the reverse transformation: the position of `T − σ/C_A`
    /// between `A_s` and `A_f` (without the stress shift if `C_A` is not
    /// positive, which [`BrinsonModel::validate`] rejects).
    fn reverse_progress(&self, sigma: f64, t: f64) -> f64 {
        let stress_shift = if self.c_a > 0.0 {
            sigma.max(0.0) / self.c_a
        } else {
            0.0
        };
        window_progress(t - stress_shift, self.as_temp, self.af_temp)
    }

    /// Progress of the stress-free cooling transformation through `M_s … M_f`.
    fn cooling_progress(&self, t: f64) -> f64 {
        window_progress(self.ms - t, 0.0, self.ms - self.mf)
    }

    /// Progress `[forward, reverse, cooling]` at stress `sigma` and temperature `t`.
    fn progress(&self, sigma: f64, t: f64) -> [f64; 3] {
        [
            self.forward_progress(sigma, t),
            self.reverse_progress(sigma, t),
            self.cooling_progress(t),
        ]
    }

    /// Frontier of the transformation memory valid for the current fractions;
    /// when the memory is missing or the fractions were edited directly, the
    /// progress at the last stress and temperature.
    fn frontier(&self) -> [f64; 3] {
        match self.memory {
            Some(memory) if memory.xi_s == self.xi_s && memory.xi_t == self.xi_t => memory.frontier,
            _ => self.progress(self.stress, self.temperature),
        }
    }

    /// Fractions after the branches whose progress at `(sigma, t)` exceeds the
    /// `frontier` have transformed.  Pure: the model is not modified.
    fn advance(&self, sigma: f64, t: f64, frontier: [f64; 3]) -> Outcome {
        let progress = self.progress(sigma, t);
        let mut xi_s = self.xi_s.clamp(0.0, 1.0);
        let mut xi_t = self.xi_t.clamp(0.0, 1.0);
        let mut advanced = false;
        if progress[0] > frontier[0] {
            // Austenite and twinned martensite are converted to detwinned
            // martensite: what is not yet detwinned decays by the survival ratio.
            let ratio = survival_ratio(frontier[0], progress[0]);
            xi_s = 1.0 - (1.0 - xi_s) * ratio;
            xi_t *= ratio;
            advanced = true;
        }
        if progress[1] > frontier[1] {
            let ratio = survival_ratio(frontier[1], progress[1]);
            xi_s *= ratio;
            xi_t *= ratio;
            advanced = true;
        }
        if progress[2] > frontier[2] {
            // The remaining austenite becomes twinned martensite.
            let ratio = survival_ratio(frontier[2], progress[2]);
            let austenite = (1.0 - xi_s - xi_t).max(0.0);
            xi_t += austenite * (1.0 - ratio);
            advanced = true;
        }
        Outcome {
            xi_s,
            xi_t,
            progress,
            advanced,
        }
    }

    /// Constitutive law `σ = E(ξ)(ε − ε_L ξ_S) + Θ (T − T₀)`.
    fn law(&self, eps: f64, xi_s: f64, xi_t: f64, t: f64) -> f64 {
        let xi = (xi_s + xi_t).clamp(0.0, 1.0);
        let modulus = self.e_austenite + xi * (self.e_martensite - self.e_austenite);
        modulus * (eps - self.eps_l * xi_s) + self.theta * (t - self.t0)
    }

    /// Store the outcome of an update as the new state.
    fn commit(&mut self, outcome: &Outcome, frontier: [f64; 3], stress: f64, t: f64) {
        self.xi_s = outcome.xi_s;
        self.xi_t = outcome.xi_t;
        self.temperature = t;
        self.stress = stress;
        self.memory = Some(Memory {
            xi_s: outcome.xi_s,
            xi_t: outcome.xi_t,
            frontier: if outcome.advanced {
                outcome.progress
            } else {
                frontier
            },
        });
    }

    /// Advance the model to the total strain `eps` at temperature `t` (K) and
    /// return the stress (Pa).
    ///
    /// The stress and the fractions are solved together: the fractions follow
    /// the transformation kinetics at the stress, and the stress follows the
    /// constitutive law at the fractions.  The strain history therefore matters
    /// only through the sequence of calls; the result for a monotone path does
    /// not depend on the step size.
    ///
    /// # Errors
    /// Returns the [`BrinsonModel::validate`] error for invalid parameters or
    /// fractions, an error for a non-finite `eps` or `t`, and
    /// [`Error::NumericalDivergence`] if the stress cannot be bracketed or the
    /// bisection does not converge.  The model is left unchanged on error.
    pub fn update(&mut self, eps: f64, t: f64) -> Result<f64> {
        self.validate()?;
        require_finite("eps", eps)?;
        require_finite("temperature", t)?;
        let frontier = self.frontier();

        // The right-hand side of the constitutive law is bounded for any
        // fractions in [0, 1]; the stress lies strictly inside this bracket.
        let modulus_max = self.e_austenite.max(self.e_martensite);
        let bound =
            modulus_max * (eps.abs() + self.eps_l) + (self.theta * (t - self.t0)).abs() + 1.0;
        if !bound.is_finite() {
            return Err(Error::diverged("brinson_stress_bracket", 0, bound));
        }
        let residual = |sigma: f64| {
            let outcome = self.advance(sigma, t, frontier);
            sigma - self.law(eps, outcome.xi_s, outcome.xi_t, t)
        };
        let (mut lo, mut hi) = (-bound, bound);
        let (r_lo, r_hi) = (residual(lo), residual(hi));
        if r_lo.is_nan() || r_hi.is_nan() || r_lo > 0.0 || r_hi < 0.0 {
            return Err(Error::diverged(
                "brinson_stress_bracket",
                0,
                r_lo.abs().max(r_hi.abs()),
            ));
        }
        let mut iterations = 0;
        loop {
            let mid = 0.5 * (lo + hi);
            let r_mid = residual(mid);
            if !r_mid.is_finite() {
                return Err(Error::diverged("brinson_stress_solve", iterations, r_mid));
            }
            if r_mid > 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
            iterations += 1;
            if hi - lo <= SOLVE_TOLERANCE * (1.0 + mid.abs()) {
                break;
            }
            if iterations >= MAX_BISECTIONS {
                return Err(Error::diverged("brinson_stress_solve", iterations, hi - lo));
            }
        }

        let sigma = 0.5 * (lo + hi);
        let outcome = self.advance(sigma, t, frontier);
        let stress = self.law(eps, outcome.xi_s, outcome.xi_t, t);
        self.commit(&outcome, frontier, stress, t);
        Ok(stress)
    }

    /// Update the temperature-induced martensite for a stress-free change of the
    /// temperature to `t` (K).
    ///
    /// Cooling below `M_s` turns the remaining austenite into twinned
    /// martensite (`ξ_T → 1 − ξ_S` at `M_f`); heating above `A_s` reverts the
    /// martensite (`ξ_T` and `ξ_S` → 0 at `A_f`).  Repeating the call at the
    /// same temperature changes nothing.  A non-finite `t` is ignored.
    pub fn update_temperature(&mut self, t: f64) {
        if !t.is_finite() {
            return;
        }
        let frontier = self.frontier();
        let outcome = self.advance(0.0, t, frontier);
        self.commit(&outcome, frontier, 0.0, t);
    }

    /// Stress (Pa) for the total strain `eps` at the current fractions and the
    /// temperature of the last update: `E(ξ)(ε − ε_L ξ_S) + Θ (T − T₀)`.
    pub fn stress(&self, eps: f64) -> f64 {
        self.law(eps, self.xi_s, self.xi_t, self.temperature)
    }

    /// Create a Brinson model with typical Nitinol (NiTi) parameters.
    ///
    /// * E_A = 75 GPa, E_M = 28 GPa
    /// * Ms = 291 K, Mf = 273 K, As = 307 K, Af = 325 K
    /// * ε_L = 0.067 (maximum transformation strain)
    /// * C_M = 8 MPa/K, C_A = 13 MPa/K, σ_s^cr = 100 MPa, σ_f^cr = 170 MPa
    ///   (the defaults of [`BrinsonModel::new`]), Θ = 0
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Test temperature 60 °C, above `A_f = 325 K`: the alloy is superelastic.
    const T_TEST: f64 = 333.0;

    /// Superelastic NiTi at `T_TEST` whose stress vanishes there for zero strain.
    fn superelastic() -> BrinsonModel {
        let mut model = BrinsonModel::nitinol().with_thermal_stress(0.0, T_TEST);
        model
            .set_state(0.0, 0.0, T_TEST)
            .expect("fractions and temperature are valid");
        model
    }

    #[derive(Debug, Clone, Copy)]
    struct Sample {
        eps: f64,
        sigma: f64,
        xi_s: f64,
        xi_t: f64,
    }

    /// Strain path `0 → peak → 0` in `steps` equal increments per branch (the
    /// initial zero is not included).
    fn triangle_path(peak: f64, steps: usize) -> Vec<f64> {
        let up = (1..=steps).map(|k| peak * k as f64 / steps as f64);
        let down = (0..steps).rev().map(|k| peak * k as f64 / steps as f64);
        up.chain(down).collect()
    }

    /// Drive `model` through `path` at constant temperature.
    fn drive(model: &mut BrinsonModel, path: &[f64], t: f64) -> Vec<Sample> {
        path.iter()
            .map(|&eps| {
                let sigma = model.update(eps, t).expect("update must converge");
                Sample {
                    eps,
                    sigma,
                    xi_s: model.xi_s,
                    xi_t: model.xi_t,
                }
            })
            .collect()
    }

    fn sample_at(samples: &[Sample], eps: f64) -> Sample {
        *samples
            .iter()
            .find(|s| (s.eps - eps).abs() < 1.0e-12)
            .expect("the path contains this strain")
    }

    /// Brinson (1993) forward closed form of `ξ_S` for `ξ_S0 = 0` above `M_s`:
    /// `ξ_S = ½·cos{π/(σ_s − σ_f)·[σ − σ_f − C_M (T − M_s)]} + ½`.
    fn brinson_forward_xi_s(model: &BrinsonModel, sigma: f64, t: f64) -> f64 {
        let arg = PI / (model.sigma_s_cr - model.sigma_f_cr)
            * (sigma - model.sigma_f_cr - model.c_m * (t - model.ms));
        0.5 * arg.cos() + 0.5
    }

    /// Brinson (1993) reverse closed form
    /// `ξ = ξ₀/2·{cos[π/(A_f − A_s)·(T − A_s − σ/C_A)] + 1}`.
    fn brinson_reverse_xi(model: &BrinsonModel, sigma: f64, t: f64, xi0: f64) -> f64 {
        let arg = PI / (model.af_temp - model.as_temp) * (t - model.as_temp - sigma / model.c_a);
        0.5 * xi0 * (arg.cos() + 1.0)
    }

    // ── Superelasticity ─────────────────────────────────────────────────────

    #[test]
    fn transformation_windows_follow_the_phase_diagram() {
        let model = BrinsonModel::nitinol();
        // Above M_s the forward window shifts up by C_M (T − M_s) = 8 MPa/K · 42 K.
        let (start, finish) = model.forward_window(T_TEST);
        assert!((start - 436.0e6).abs() < 1.0 && (finish - 506.0e6).abs() < 1.0);
        // Below M_s it is the plain detwinning window.
        let (cold_start, cold_finish) = model.forward_window(250.0);
        assert!((cold_start - 100.0e6).abs() < 1.0 && (cold_finish - 170.0e6).abs() < 1.0);
        // Reverse window: C_A (T − A_s) down to C_A (T − A_f) = 13 MPa/K · (26 K, 8 K).
        let (rev_start, rev_finish) = model.reverse_window(T_TEST);
        assert!((rev_start - 338.0e6).abs() < 1.0 && (rev_finish - 104.0e6).abs() < 1.0);
        // Hysteresis: the reverse window lies below the forward window.
        assert!(rev_start < start);
    }

    #[test]
    fn superelastic_loop_closes_with_hysteresis() {
        let mut model = superelastic();
        let steps = 1000;
        let samples = drive(&mut model, &triangle_path(0.10, steps), T_TEST);
        let (loading, unloading) = samples.split_at(steps);

        // Initial elastic branch: σ = E_A ε, no transformation.
        let elastic = sample_at(loading, 0.002);
        assert!((elastic.sigma - 75.0e9 * 0.002).abs() < 1.0e-6 * elastic.sigma);
        assert_eq!((elastic.xi_s, elastic.xi_t), (0.0, 0.0));

        // Loading plateau: inside the forward window, ξ_S follows Brinson's
        // closed form from ξ_S0 = 0, and no twinned martensite appears.
        let (start, finish) = model.forward_window(T_TEST);
        let plateau: Vec<&Sample> = loading
            .iter()
            .filter(|s| s.xi_s > 1.0e-6 && s.xi_s < 1.0 - 1.0e-6)
            .collect();
        assert!(plateau.len() > 500, "plateau samples: {}", plateau.len());
        for s in &plateau {
            assert!(
                s.sigma > start - 1.0 && s.sigma < finish + 1.0,
                "plateau stress {} outside [{start}, {finish}]",
                s.sigma
            );
            let expected = brinson_forward_xi_s(&model, s.sigma, T_TEST);
            assert!(
                (s.xi_s - expected).abs() < 1.0e-9,
                "ξ_S = {} vs Brinson {expected} at σ = {}",
                s.xi_s,
                s.sigma
            );
            assert_eq!(s.xi_t, 0.0);
        }
        assert!(
            loading.windows(2).all(|w| w[1].xi_s >= w[0].xi_s),
            "ξ_S must not decrease while loading"
        );

        // End of loading: fully transformed, σ = E_M (ε − ε_L).
        let peak = loading[steps - 1];
        assert!((peak.xi_s - 1.0).abs() < 1.0e-12 && peak.xi_t == 0.0);
        assert!((peak.sigma - 28.0e9 * (0.10 - 0.067)).abs() < 1.0);

        // Unloading: elastic along E_M first, then the reverse plateau, whose
        // fractions follow Brinson's reverse closed form with ξ₀ = 1.
        let (rev_start, rev_finish) = model.reverse_window(T_TEST);
        let reverse: Vec<&Sample> = unloading
            .iter()
            .filter(|s| s.xi_s > 1.0e-6 && s.xi_s < 1.0 - 1.0e-6)
            .collect();
        assert!(
            reverse.len() > 500,
            "reverse plateau samples: {}",
            reverse.len()
        );
        for s in &reverse {
            assert!(
                s.sigma < rev_start + 1.0 && s.sigma > rev_finish - 1.0,
                "reverse stress {} outside [{rev_finish}, {rev_start}]",
                s.sigma
            );
            let expected = brinson_reverse_xi(&model, s.sigma, T_TEST, 1.0);
            assert!(
                (s.xi_s + s.xi_t - expected).abs() < 1.0e-9,
                "ξ = {} vs Brinson {expected} at σ = {}",
                s.xi_s + s.xi_t,
                s.sigma
            );
            assert_eq!(s.xi_t, 0.0);
        }
        assert!(
            unloading.windows(2).all(|w| w[1].xi_s <= w[0].xi_s),
            "ξ_S must not increase while unloading"
        );

        // The loop closes: zero strain, zero stress, zero martensite.
        let end = unloading[steps - 1];
        assert_eq!(end.eps, 0.0);
        assert!(end.sigma.abs() < 1.0e-3, "residual stress {}", end.sigma);
        assert_eq!((end.xi_s, end.xi_t), (0.0, 0.0));

        // Hysteresis: at the same strain the loading stress is far above the
        // unloading stress.
        let up = sample_at(loading, 0.04);
        let down = sample_at(unloading, 0.04);
        assert!(up.sigma > 436.0e6 - 1.0 && down.sigma < 338.0e6 + 1.0);
        assert!(up.sigma - down.sigma > 100.0e6);
        assert!(up.xi_s < 0.5 + 0.1 && down.xi_s > up.xi_s);

        // Dissipated energy per cycle (loop area) is of the order of
        // (Δσ_plateau · ε_L) ≈ 17 MJ/m³.
        let mut area = 0.0;
        let mut previous = (0.0, 0.0);
        for s in &samples {
            area += 0.5 * (s.sigma + previous.1) * (s.eps - previous.0);
            previous = (s.eps, s.sigma);
        }
        assert!(
            (15.0e6..22.0e6).contains(&area),
            "dissipated energy {area:e} J/m³"
        );
    }

    #[test]
    fn superelastic_response_does_not_depend_on_the_step_size() {
        let mut fine = superelastic();
        let fine_samples = drive(&mut fine, &triangle_path(0.10, 1000), T_TEST);
        let mut coarse = superelastic();
        let coarse_samples = drive(&mut coarse, &triangle_path(0.10, 10), T_TEST);
        for (branch_fine, branch_coarse) in [
            (&fine_samples[..1000], &coarse_samples[..10]),
            (&fine_samples[1000..], &coarse_samples[10..]),
        ] {
            for eps in [0.04, 0.06, 0.08] {
                let a = sample_at(branch_fine, eps);
                let b = sample_at(branch_coarse, eps);
                assert!(
                    (a.sigma - b.sigma).abs() < 1.0e-9 * a.sigma.abs().max(1.0),
                    "σ({eps}): {} vs {}",
                    a.sigma,
                    b.sigma
                );
                assert!((a.xi_s - b.xi_s).abs() < 1.0e-9, "ξ_S({eps})");
            }
        }
        let (a, b) = (fine_samples[1999], coarse_samples[19]);
        assert_eq!((a.xi_s, a.xi_t), (b.xi_s, b.xi_t));
    }

    #[test]
    fn stress_induced_martensite_evolves_under_load() {
        // The earlier model never wrote ξ_S: loading did not transform anything.
        let mut model = superelastic();
        assert_eq!(model.xi_s, 0.0);
        model.update(0.06, T_TEST).expect("update");
        assert!(model.xi_s > 0.5, "ξ_S = {}", model.xi_s);
        assert_eq!(model.xi_t, 0.0);
    }

    #[test]
    fn no_transformation_below_the_critical_stress() {
        let mut model = superelastic();
        for k in 1..=20 {
            let eps = 0.0002 * f64::from(k);
            let sigma = model.update(eps, T_TEST).expect("update");
            assert!((sigma - 75.0e9 * eps).abs() < 1.0e-6 * sigma);
            assert_eq!((model.xi_s, model.xi_t), (0.0, 0.0));
        }
    }

    // ── Inner loops: partial transformations ────────────────────────────────

    fn load_to(model: &mut BrinsonModel, from: f64, to: f64, steps: usize) {
        for k in 1..=steps {
            let eps = from + (to - from) * k as f64 / steps as f64;
            model.update(eps, T_TEST).expect("update");
        }
    }

    #[test]
    fn partial_unloading_is_elastic_and_reloading_retraces_without_a_jump() {
        // Reference: monotone loading to ε = 0.08.
        let mut reference = superelastic();
        load_to(&mut reference, 0.0, 0.08, 800);

        let mut model = superelastic();
        load_to(&mut model, 0.0, 0.05, 500);
        let (sigma_top, xi_top) = (model.current_stress(), model.xi_s);
        assert!(
            xi_top > 0.3 && xi_top < 0.9,
            "partial transformation: {xi_top}"
        );

        // Unload by 0.001 (stays above the reverse window): purely elastic with
        // the current modulus, ξ unchanged.
        let modulus = model.effective_modulus();
        load_to(&mut model, 0.05, 0.049, 10);
        assert_eq!(
            model.xi_s, xi_top,
            "no transformation while unloading elastically"
        );
        let drop = sigma_top - model.current_stress();
        assert!(
            (drop - modulus * 0.001).abs() < 1.0e-9 * drop,
            "elastic drop {drop} vs E ε = {}",
            modulus * 0.001
        );

        // Reloading returns to the same point and continues on the original path.
        load_to(&mut model, 0.049, 0.05, 10);
        assert!((model.current_stress() - sigma_top).abs() < 1.0e-9 * sigma_top);
        assert!((model.xi_s - xi_top).abs() < 1.0e-12);
        load_to(&mut model, 0.05, 0.08, 300);
        assert!(
            (model.current_stress() - reference.current_stress()).abs()
                < 1.0e-9 * reference.current_stress()
        );
        assert!((model.xi_s - reference.xi_s).abs() < 1.0e-12);
    }

    #[test]
    fn partial_reverse_transformation_followed_by_reloading_reaches_full_martensite() {
        let mut model = superelastic();
        load_to(&mut model, 0.0, 0.10, 1000);
        load_to(&mut model, 0.10, 0.074, 260);
        let xi_partial = model.xi_s;
        assert!(
            xi_partial > 0.5 && xi_partial < 1.0,
            "partial reverse transformation: {xi_partial}"
        );
        load_to(&mut model, 0.074, 0.10, 260);
        assert!((model.xi_s - 1.0).abs() < 1.0e-12 && model.xi_t == 0.0);
        assert!((model.current_stress() - 28.0e9 * (0.10 - 0.067)).abs() < 1.0);
    }

    // ── Shape memory effect (low temperature, thermal cycling) ──────────────

    #[test]
    fn cooling_forms_twinned_martensite_which_detwins_under_load_and_recovers_on_heating() {
        let mut model = BrinsonModel::nitinol();
        model.update_temperature(250.0);
        assert_eq!(
            (model.xi_s, model.xi_t),
            (0.0, 1.0),
            "fully twinned below M_f"
        );

        // Detwinning at 250 K < M_s: plain window σ_s^cr … σ_f^cr, ξ_S + ξ_T = 1.
        let samples = drive(&mut model, &triangle_path(0.10, 1000)[..1000], 250.0);
        let (start, finish) = model.forward_window(250.0);
        let plateau: Vec<&Sample> = samples
            .iter()
            .filter(|s| s.xi_s > 1.0e-6 && s.xi_s < 1.0 - 1.0e-6)
            .collect();
        assert!(plateau.len() > 100);
        for s in &plateau {
            assert!(s.sigma > start - 1.0 && s.sigma < finish + 1.0);
            assert!(
                (s.xi_s + s.xi_t - 1.0).abs() < 1.0e-12,
                "no austenite present"
            );
        }
        assert!((model.xi_s - 1.0).abs() < 1.0e-12 && model.xi_t == 0.0);

        // Unloading at 250 K < A_s does not transform back: residual strain ε_L.
        let unloaded = drive(&mut model, &[0.05, 0.0], 250.0);
        assert_eq!(unloaded[1].xi_s, 1.0);
        let residual = model.update(0.067, 250.0).expect("update");
        assert!(residual.abs() < 1.0e-3, "stress at ε = ε_L: {residual}");

        // Stress-free heating recovers the shape between A_s and A_f:
        // free strain ε_L ξ_S = ε_L at A_s, ε_L/2 at the middle, 0 at A_f.
        model.update_temperature(307.0);
        assert!((model.xi_s - 1.0).abs() < 1.0e-12);
        model.update_temperature(316.0);
        assert!((model.xi_s - 0.5).abs() < 1.0e-12, "ξ_S = {}", model.xi_s);
        model.update_temperature(325.0);
        assert_eq!((model.xi_s, model.xi_t), (0.0, 0.0));
    }

    #[test]
    fn heating_does_not_compound_and_repeated_updates_are_idempotent() {
        // Mid-way through the austenite window the twinned martensite is halved
        // once, however often the temperature is applied (the old update
        // multiplied by 0.5 on every call).
        let mut model = BrinsonModel::new(70.0e9, 30.0e9, 300.0, 200.0, 340.0, 380.0, 0.05);
        model.update_temperature(150.0);
        assert_eq!(model.xi_t, 1.0);
        for _ in 0..5 {
            model.update_temperature(360.0);
            assert!((model.xi_t - 0.5).abs() < 1.0e-12, "ξ_T = {}", model.xi_t);
        }
        // A quarter of the way: (1 + cos π/4)/2.
        let mut quarter = BrinsonModel::new(70.0e9, 30.0e9, 300.0, 200.0, 340.0, 380.0, 0.05);
        quarter.update_temperature(150.0);
        quarter.update_temperature(350.0);
        let expected = 0.5 * (1.0 + (PI * 0.25).cos());
        assert!((quarter.xi_t - expected).abs() < 1.0e-12);
        // Cooling back above M_s does not re-form martensite.
        quarter.update_temperature(320.0);
        assert!((quarter.xi_t - expected).abs() < 1.0e-12);
    }

    #[test]
    fn cooling_is_independent_of_the_temperature_step() {
        let reference = |temperatures: &[f64]| {
            let mut model = BrinsonModel::new(70.0e9, 30.0e9, 300.0, 200.0, 340.0, 380.0, 0.05);
            for &t in temperatures {
                model.update_temperature(t);
            }
            model.xi_t
        };
        let direct = reference(&[250.0]);
        let stepped = reference(&[390.0, 310.0, 290.0, 280.0, 270.0, 260.0, 250.0]);
        // 0.5 (1 − cos(π (T − M_s)/(M_f − M_s))) at T = 250: 0.5.
        assert!((direct - 0.5).abs() < 1.0e-12);
        assert!((stepped - direct).abs() < 1.0e-12);
    }

    // ── Constitutive law ────────────────────────────────────────────────────

    #[test]
    fn stress_has_no_invented_twinned_martensite_term() {
        // σ = E(ξ)(ε − ε_L ξ_S): twinned martensite (ξ_T) produces no strain.
        let mut model = BrinsonModel::nitinol();
        model.set_state(0.0, 0.5, 300.0).expect("valid state");
        let modulus = 75.0e9 + 0.5 * (28.0e9 - 75.0e9);
        assert!((model.stress(0.01) - modulus * 0.01).abs() < 1.0e-3);
        // Stress-induced martensite does: ε_L ξ_S is subtracted from the strain.
        model.set_state(0.4, 0.0, 300.0).expect("valid state");
        let modulus = 75.0e9 + 0.4 * (28.0e9 - 75.0e9);
        assert!((model.stress(0.05) - modulus * (0.05 - 0.067 * 0.4)).abs() < 1.0e-3);
    }

    #[test]
    fn thermoelastic_term_enters_the_stress() {
        let mut model = BrinsonModel::nitinol().with_thermal_stress(0.55e6, 300.0);
        let sigma = model.update(0.0, 340.0).expect("update");
        assert!((sigma - 0.55e6 * 40.0).abs() < 1.0e-3, "σ = {sigma}");
        assert_eq!(model.temperature(), 340.0);
        assert!((model.stress(0.001) - (75.0e9 * 0.001 + 22.0e6)).abs() < 1.0e-3);
    }

    #[test]
    fn steep_transformation_window_is_solved_without_divergence() {
        // A 1 kPa wide window: the plateau is almost flat, the map σ ↦ E(ξ(σ))(ε − ε_L ξ_S)
        // has a huge slope and would defeat a fixed-point iteration.
        let mut model =
            superelastic().with_stress_kinetics(8.0e6, 13.0e6, 100.0e6, 100.0e6 + 1.0e3);
        let (start, finish) = model.forward_window(T_TEST);
        let mut previous = 0.0;
        for k in 1..=200 {
            let eps = 0.0005 * f64::from(k);
            let sigma = model.update(eps, T_TEST).expect("update");
            assert!(sigma.is_finite());
            if model.xi_s > 1.0e-9 && model.xi_s < 1.0 - 1.0e-9 {
                // The window is 1 kPa wide, so the map's slope is about 5e6 and a
                // stress error of a few ulp (~1e-7 Pa) moves the law by ~1 Pa.
                assert!(
                    (start - 50.0..=finish + 50.0).contains(&sigma),
                    "σ = {sigma} outside the 1 kPa window [{start}, {finish}]"
                );
            }
            assert!(model.xi_s >= previous);
            previous = model.xi_s;
        }
        assert_eq!(model.xi_s, 1.0);
    }

    // ── Validation and robustness ───────────────────────────────────────────

    /// A function that breaks one parameter or state field of a model.
    type Corruption = fn(&mut BrinsonModel);

    #[test]
    fn invalid_parameters_are_rejected_and_leave_the_model_unchanged() {
        let reference = superelastic();
        let cases: [(&'static str, Corruption); 6] = [
            ("c_a", |m| m.c_a = 0.0),
            ("sigma_f_cr", |m| m.sigma_f_cr = m.sigma_s_cr),
            ("ms", |m| m.ms = m.mf),
            ("af_temp", |m| m.af_temp = m.as_temp),
            ("e_martensite", |m| m.e_martensite = -1.0),
            ("xi_s", |m| m.xi_s = 1.5),
        ];
        for (field, corrupt) in cases {
            let mut model = reference.clone();
            corrupt(&mut model);
            let before = (
                model.xi_s,
                model.xi_t,
                model.current_stress(),
                model.temperature(),
            );
            let error = model.update(0.05, T_TEST).expect_err("must be rejected");
            match error {
                Error::ParameterOutOfRange { param, .. } => assert_eq!(param, field),
                other => panic!("unexpected error for {field}: {other}"),
            }
            assert_eq!(
                before,
                (
                    model.xi_s,
                    model.xi_t,
                    model.current_stress(),
                    model.temperature()
                ),
                "state changed despite the error for {field}"
            );
        }
    }

    #[test]
    fn non_finite_inputs_and_unbracketable_stresses_are_errors() {
        let mut model = superelastic();
        assert!(model.update(f64::NAN, T_TEST).is_err());
        assert!(model.update(0.01, f64::INFINITY).is_err());
        model.update_temperature(f64::NAN);
        assert_eq!(model.temperature(), T_TEST, "a NaN temperature is ignored");

        let mut huge = superelastic();
        huge.e_austenite = f64::MAX;
        match huge.update(10.0, T_TEST) {
            Err(Error::NumericalDivergence { solver, .. }) => {
                assert_eq!(solver, "brinson_stress_bracket");
            }
            other => panic!("expected NumericalDivergence, got {other:?}"),
        }
    }

    #[test]
    fn directly_edited_fractions_restart_the_transformation_memory() {
        let mut model = BrinsonModel::nitinol();
        model.update_temperature(250.0);
        assert_eq!(model.xi_t, 1.0);
        // Overwrite the fraction behind the model's back: the cooling memory of
        // the old state must not apply to the new one.
        model.xi_t = 0.6;
        model.update_temperature(250.0);
        assert_eq!(model.xi_t, 0.6);
        model.update_temperature(316.0);
        assert!((model.xi_t - 0.3).abs() < 1.0e-12, "ξ_T = {}", model.xi_t);
    }

    #[test]
    fn set_state_validates_and_anchors_the_memory() {
        let mut model = BrinsonModel::nitinol();
        assert!(model.set_state(0.7, 0.5, 300.0).is_err(), "ξ_S + ξ_T > 1");
        assert!(model.set_state(-0.1, 0.0, 300.0).is_err());
        assert!(model.set_state(0.0, 0.0, f64::NAN).is_err());
        // Austenite already at the (cold) test temperature does not "cool".
        model.set_state(0.0, 0.0, 280.0).expect("valid state");
        model.update_temperature(280.0);
        assert_eq!((model.xi_s, model.xi_t), (0.0, 0.0));
        // Further cooling does form martensite.
        model.update_temperature(273.0);
        assert!((model.xi_t - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn nitinol_parameters_are_valid() {
        BrinsonModel::nitinol()
            .validate()
            .expect("preset must validate");
        let custom = BrinsonModel::new(70.0e9, 30.0e9, 300.0, 200.0, 340.0, 380.0, 0.05)
            .with_stress_kinetics(6.5e6, 6.5e6, 90.0e6, 160.0e6);
        custom.validate().expect("custom parameters validate");
        assert_eq!(custom.c_m, 6.5e6);
        assert_eq!(custom.sigma_f_cr, 160.0e6);
    }
}
