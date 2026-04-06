//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{normal_quantile, standard_normal_cdf};

/// Coffin-Manson plastic strain-life equation.
///
/// Δε_p/2 = ε_f' * (2N_f)^c
///
/// Relates the plastic strain amplitude to the number of reversals to failure.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CoffinManson {
    /// Fatigue ductility coefficient ε_f'
    pub epsilon_f: f64,
    /// Fatigue ductility exponent c (negative, typically -0.5 to -0.7)
    pub c: f64,
}
#[allow(dead_code)]
impl CoffinManson {
    /// Create a new Coffin-Manson model.
    pub fn new(epsilon_f: f64, c: f64) -> Self {
        Self { epsilon_f, c }
    }
    /// Plastic strain amplitude at 2N reversals.
    pub fn plastic_strain_amplitude(&self, two_n: f64) -> f64 {
        self.epsilon_f * two_n.powf(self.c)
    }
    /// Reversals to failure at a given plastic strain amplitude.
    ///
    /// 2N_f = (Δε_p/2 / ε_f')^(1/c)
    pub fn reversals_to_failure(&self, plastic_strain_amp: f64) -> f64 {
        (plastic_strain_amp / self.epsilon_f).powf(1.0 / self.c)
    }
    /// Cycles to failure (N_f = 2N_f / 2).
    pub fn cycles_to_failure(&self, plastic_strain_amp: f64) -> f64 {
        self.reversals_to_failure(plastic_strain_amp) / 2.0
    }
    /// Transition life: number of reversals where elastic and plastic strains are equal.
    ///
    /// At the transition life 2N_t:
    /// (σ_f'/E) * (2N_t)^b = ε_f' * (2N_t)^c
    /// → 2N_t = (ε_f' * E / σ_f')^(1/(b-c))
    pub fn transition_life(&self, sigma_f: f64, b: f64, e_modulus: f64) -> f64 {
        let ratio = self.epsilon_f * e_modulus / sigma_f;
        ratio.powf(1.0 / (b - self.c))
    }
}
/// Palmgren-Miner cumulative damage rule with pre-computed damage ratios.
///
/// Each element of `damages` is the ratio n_i/N_i for one loading block.
/// Failure is predicted when the sum reaches 1.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MinerRule {
    /// Individual damage contributions D_i = n_i / N_i.
    pub damages: Vec<f64>,
}
#[allow(dead_code)]
impl MinerRule {
    /// Create a new Miner rule accumulator.
    pub fn new(damages: Vec<f64>) -> Self {
        Self { damages }
    }
    /// Total accumulated damage D = Σ D_i.
    pub fn total_damage(&self) -> f64 {
        self.damages.iter().sum()
    }
    /// Returns `true` when total damage D ≥ 1 (failure predicted).
    pub fn is_failed(&self) -> bool {
        self.total_damage() >= 1.0
    }
}
/// Biaxial fatigue failure surface (Sines-type ellipse).
///
/// Failure when: (σ_a / σ_e)^2 + a*σ_a*σ_b/σ_e^2 + (σ_b / σ_e)^2 = 1
/// A simplified form using parameters `a` (interaction) and `b` (biaxiality ratio).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FatigueSurface {
    /// Interaction coefficient a.
    pub a: f64,
    /// Biaxiality coefficient b.
    pub b: f64,
    /// Ultimate tensile strength σ_ult (Pa).
    pub sigma_ult: f64,
    /// Endurance limit σ_e (Pa).
    pub sigma_endurance: f64,
}
#[allow(dead_code)]
impl FatigueSurface {
    /// Create a new biaxial fatigue surface.
    pub fn new(a: f64, b: f64, sigma_ult: f64, sigma_endurance: f64) -> Self {
        Self {
            a,
            b,
            sigma_ult,
            sigma_endurance,
        }
    }
    /// Evaluate the biaxial damage index for stress amplitudes (σ_a1, σ_a2).
    ///
    /// Returns a value < 1 for safe, ≥ 1 for failure.
    pub fn damage_index(&self, sigma_a1: f64, sigma_a2: f64) -> f64 {
        let se = self.sigma_endurance;
        if se <= 0.0 {
            return f64::INFINITY;
        }
        let r1 = sigma_a1 / se;
        let r2 = sigma_a2 / se;
        r1 * r1 + self.a * r1 * r2 + self.b * r2 * r2
    }
    /// Returns `true` when the stress state is outside the failure surface.
    pub fn is_failed(&self, sigma_a1: f64, sigma_a2: f64) -> bool {
        self.damage_index(sigma_a1, sigma_a2) >= 1.0
    }
}
/// Fatigue life predictor combining a Basquin S-N curve with a Goodman
/// mean-stress correction.
///
/// N = (σ_ar / A)^(1/B)  where σ_ar = σ_a / (1 - σ_m / σ_ult)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct FatigueLifePredictor {
    /// Basquin S-N curve.
    pub sn: BasquinCurve,
    /// Ultimate tensile strength for Goodman correction (Pa).
    pub sigma_ult: f64,
}
#[allow(dead_code)]
impl FatigueLifePredictor {
    /// Create a new predictor.
    pub fn new(sn: BasquinCurve, sigma_ult: f64) -> Self {
        Self { sn, sigma_ult }
    }
    /// Predict cycles to failure for a given stress amplitude and mean stress.
    ///
    /// Uses Goodman correction: σ_ar = σ_a / (1 - σ_m / σ_ult)
    /// then applies the Basquin S-N curve.
    pub fn predict_cycles(&self, sigma_a: f64, sigma_m: f64) -> f64 {
        let denom = 1.0 - sigma_m / self.sigma_ult;
        if denom <= 0.0 {
            return 0.0;
        }
        let sigma_ar = sigma_a / denom;
        self.sn.cycles_to_failure(sigma_ar)
    }
}
/// Morrow mean stress correction for strain-life fatigue.
///
/// Modified Basquin equation:
/// Δε_e/2 = ((σ_f' - σ_m) / E) * (2N_f)^b
///
/// The mean stress σ_m reduces the effective fatigue strength coefficient.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MorrowCorrection {
    /// Fatigue strength coefficient σ_f' (Pa)
    pub sigma_f: f64,
    /// Fatigue strength exponent b
    pub b: f64,
    /// Young's modulus E (Pa)
    pub e_modulus: f64,
}
#[allow(dead_code)]
impl MorrowCorrection {
    /// Create a new Morrow correction model.
    pub fn new(sigma_f: f64, b: f64, e_modulus: f64) -> Self {
        Self {
            sigma_f,
            b,
            e_modulus,
        }
    }
    /// Corrected elastic strain amplitude with mean stress.
    ///
    /// Δε_e/2 = ((σ_f' - σ_m) / E) * (2N)^b
    pub fn corrected_strain_amplitude(&self, mean_stress: f64, two_n: f64) -> f64 {
        ((self.sigma_f - mean_stress) / self.e_modulus) * two_n.powf(self.b)
    }
    /// Effective fatigue strength coefficient with mean stress correction.
    ///
    /// σ_f_eff = σ_f' - σ_m
    pub fn effective_sigma_f(&self, mean_stress: f64) -> f64 {
        self.sigma_f - mean_stress
    }
    /// Cycles to failure with Morrow correction (bisection).
    ///
    /// Solves: strain_amp = ((σ_f' - σ_m) / E) * (2N)^b for N.
    pub fn cycles_to_failure(&self, strain_amplitude: f64, mean_stress: f64) -> f64 {
        let sigma_eff = self.effective_sigma_f(mean_stress);
        if sigma_eff <= 0.0 {
            return 0.0;
        }
        let two_n = (strain_amplitude * self.e_modulus / sigma_eff).powf(1.0 / self.b);
        two_n / 2.0
    }
}
/// Soderberg diagram defined by yield strength and endurance limit.
///
/// Allowable amplitude: σ_a = σ_e * (1 - σ_m / σ_yield)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SoderbergDiagram {
    /// Yield strength σ_yield (Pa).
    pub sigma_yield: f64,
    /// Fully-reversed endurance limit σ_e (Pa).
    pub sigma_endurance: f64,
}
#[allow(dead_code)]
impl SoderbergDiagram {
    /// Create a new Soderberg diagram.
    pub fn new(sigma_yield: f64, sigma_endurance: f64) -> Self {
        Self {
            sigma_yield,
            sigma_endurance,
        }
    }
    /// Allowable stress amplitude for a given mean stress.
    ///
    /// σ_a_allow = σ_e * (1 - σ_m / σ_yield)  \[clamped to 0\]
    pub fn soderberg_allowed_amplitude(&self, mean_stress: f64) -> f64 {
        (self.sigma_endurance * (1.0 - mean_stress / self.sigma_yield)).max(0.0)
    }
}
/// Miner's rule with a user-specified critical damage threshold D_crit.
///
/// The standard Palmgren-Miner rule uses D_crit = 1.0, but experimental
/// evidence shows that actual failure may occur between D = 0.3 and 3.0.
/// This struct allows a custom D_crit.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MinerWithCrit {
    /// Accumulated damage D.
    pub damage: f64,
    /// Critical damage threshold D_crit (failure when D >= D_crit).
    pub d_crit: f64,
}
#[allow(dead_code)]
impl MinerWithCrit {
    /// Create a new Miner rule accumulator with a custom critical damage.
    pub fn new(d_crit: f64) -> Self {
        Self {
            damage: 0.0,
            d_crit,
        }
    }
    /// Add damage from one loading block.
    pub fn add_block(&mut self, n_applied: f64, n_failure: f64) {
        if n_failure > 0.0 {
            self.damage += n_applied / n_failure;
        }
    }
    /// Returns true if accumulated damage >= D_crit.
    pub fn is_failed(&self) -> bool {
        self.damage >= self.d_crit
    }
    /// Remaining life fraction: (D_crit - D) / D_crit.
    pub fn remaining_life_fraction(&self) -> f64 {
        ((self.d_crit - self.damage) / self.d_crit).max(0.0)
    }
    /// Number of cycles to failure at a constant stress level N_f.
    pub fn cycles_to_failure(&self, n_f: f64) -> f64 {
        if self.damage >= self.d_crit {
            return 0.0;
        }
        (self.d_crit - self.damage) * n_f
    }
}
/// Cyclic Ramberg-Osgood stress-strain model.
///
/// ε = σ/E + (σ/K')^(1/n')
///
/// Describes the cyclic (stabilised hysteresis loop) stress-strain behaviour.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RambergOsgood {
    /// Young's modulus E (Pa).
    pub e_modulus: f64,
    /// Cyclic strength coefficient K' (Pa).
    pub k_prime: f64,
    /// Cyclic strain hardening exponent n' (dimensionless, 0 < n' < 1).
    pub n_prime: f64,
}
#[allow(dead_code)]
impl RambergOsgood {
    /// Create a new Ramberg-Osgood model.
    pub fn new(e_modulus: f64, k_prime: f64, n_prime: f64) -> Self {
        Self {
            e_modulus,
            k_prime,
            n_prime,
        }
    }
    /// Total strain ε for a given stress σ.
    pub fn strain(&self, sigma: f64) -> f64 {
        sigma / self.e_modulus
            + (sigma.abs() / self.k_prime).powf(1.0 / self.n_prime) * sigma.signum()
    }
    /// Elastic strain σ/E.
    pub fn elastic_strain(&self, sigma: f64) -> f64 {
        sigma / self.e_modulus
    }
    /// Plastic strain (σ/K')^(1/n').
    pub fn plastic_strain(&self, sigma: f64) -> f64 {
        (sigma.abs() / self.k_prime).powf(1.0 / self.n_prime) * sigma.signum()
    }
    /// Solve for stress given total strain using bisection.
    pub fn stress_from_strain(&self, epsilon: f64, max_iter: usize) -> f64 {
        if epsilon.abs() < 1e-30 {
            return 0.0;
        }
        let sign = epsilon.signum();
        let eps_abs = epsilon.abs();
        let sigma_max = eps_abs * self.e_modulus * 10.0;
        let mut lo = 0.0_f64;
        let mut hi = sigma_max;
        for _ in 0..max_iter {
            let mid = 0.5 * (lo + hi);
            let eps_mid = self.strain(mid * sign).abs();
            if eps_mid < eps_abs {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) / sigma_max < 1e-12 {
                break;
            }
        }
        sign * 0.5 * (lo + hi)
    }
    /// Hysteresis loop width (plastic strain range) for a stress range Δσ.
    ///
    /// Δε_p = 2 * (Δσ / (2*K'))^(1/n')  (Masing hypothesis)
    pub fn plastic_strain_range(&self, delta_sigma: f64) -> f64 {
        2.0 * (delta_sigma / (2.0 * self.k_prime)).powf(1.0 / self.n_prime)
    }
    /// Hysteresis loop area (energy per cycle per unit volume) ≈ Δσ * Δε_p * 4/π.
    ///
    /// For a Ramberg-Osgood loop the exact energy per cycle is:
    /// W = Δσ * Δε_p * 2*n'/(n'+1)
    pub fn energy_per_cycle(&self, delta_sigma: f64) -> f64 {
        let dep = self.plastic_strain_range(delta_sigma);
        delta_sigma * dep * 2.0 * self.n_prime / (self.n_prime + 1.0)
    }
}
/// Palmgren-Miner cumulative damage model.
///
/// Linear damage accumulation: D = Σ (n_i / N_i)
/// Failure is predicted when D ≥ 1.
#[derive(Debug, Clone)]
pub struct PalmgrenMinor {
    /// Cumulative damage D
    pub damage: f64,
}
impl PalmgrenMinor {
    /// Create a new Palmgren-Miner damage accumulator with D = 0.
    pub fn new() -> Self {
        Self { damage: 0.0 }
    }
    /// Apply a block of n_applied cycles at a stress level with n_failure cycles to failure.
    ///
    /// D += n_applied / n_failure
    ///
    /// # Arguments
    /// * `n_applied` - Number of applied cycles at this stress level
    /// * `n_failure` - Cycles to failure at this stress level
    pub fn apply_cycle_block(&mut self, n_applied: f64, n_failure: f64) {
        self.damage += n_applied / n_failure;
    }
    /// Return the total accumulated damage D.
    pub fn total_damage(&self) -> f64 {
        self.damage
    }
    /// Returns true if D ≥ 1.0 (failure predicted by Miner's rule).
    pub fn is_failed(&self) -> bool {
        self.damage >= 1.0
    }
    /// Reset the damage accumulator to zero.
    pub fn reset(&mut self) {
        self.damage = 0.0;
    }
    /// Remaining life fraction: 1.0 - D (clamped to \[0, 1\]).
    #[allow(dead_code)]
    pub fn remaining_life(&self) -> f64 {
        (1.0 - self.damage).clamp(0.0, 1.0)
    }
    /// Apply damage from a list of (n_applied, n_failure) pairs.
    #[allow(dead_code)]
    pub fn apply_spectrum(&mut self, blocks: &[(f64, f64)]) {
        for &(n_applied, n_failure) in blocks {
            self.apply_cycle_block(n_applied, n_failure);
        }
    }
    /// Compute the number of additional cycles at a given stress level
    /// before failure (D reaches 1.0).
    #[allow(dead_code)]
    pub fn remaining_cycles(&self, n_failure_at_stress: f64) -> f64 {
        if self.damage >= 1.0 {
            return 0.0;
        }
        (1.0 - self.damage) * n_failure_at_stress
    }
}
/// Explicit Basquin stress-life solver.
///
/// Stress amplitude: σ_a = σ_f_prime * (2N)^b
/// Inverted to give N_f = 0.5 * (σ_a / σ_f_prime)^(1/b)
#[allow(dead_code)]
pub struct BasquinStressLife {
    /// Fatigue strength coefficient σ_f' (Pa)
    pub sigma_f_prime: f64,
    /// Fatigue strength exponent b (negative)
    pub b_exp: f64,
    /// Endurance limit (Pa) — cycles beyond which no damage is counted
    pub endurance_limit: f64,
}
impl BasquinStressLife {
    /// Create a new Basquin stress-life model.
    pub fn new(sigma_f_prime: f64, b_exp: f64, endurance_limit: f64) -> Self {
        Self {
            sigma_f_prime,
            b_exp,
            endurance_limit,
        }
    }
    /// Stress amplitude at a given number of reversals 2N.
    pub fn stress_amplitude(&self, n_reversals: f64) -> f64 {
        self.sigma_f_prime * n_reversals.powf(self.b_exp)
    }
    /// Cycles to failure for a given stress amplitude.
    ///
    /// Returns `f64::INFINITY` when `sigma_a <= endurance_limit`.
    pub fn cycles_to_failure(&self, sigma_a: f64) -> f64 {
        if sigma_a <= self.endurance_limit {
            return f64::INFINITY;
        }
        let two_n_f = (sigma_a / self.sigma_f_prime).powf(1.0 / self.b_exp);
        two_n_f / 2.0
    }
    /// Fatigue damage per cycle at given stress amplitude.
    ///
    /// D = 1 / N_f.  Returns 0 when below endurance limit.
    pub fn damage_per_cycle(&self, sigma_a: f64) -> f64 {
        let n_f = self.cycles_to_failure(sigma_a);
        if n_f.is_infinite() { 0.0 } else { 1.0 / n_f }
    }
    /// Endurance ratio: σ_endurance / σ_f_prime.
    pub fn endurance_ratio(&self) -> f64 {
        self.endurance_limit / self.sigma_f_prime
    }
    /// Cycles at which the S-N curve crosses the endurance limit.
    pub fn transition_cycles(&self) -> f64 {
        self.cycles_to_failure(self.endurance_limit * 1.000_000_001)
    }
    /// Apply a safety factor to the stress amplitude before computing life.
    pub fn cycles_with_safety_factor(&self, sigma_a: f64, safety_factor: f64) -> f64 {
        self.cycles_to_failure(sigma_a * safety_factor)
    }
}
/// A counted fatigue cycle from rainflow analysis.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct RainflowCycle {
    /// Cycle range (peak − valley).
    pub range: f64,
    /// Cycle mean ((peak + valley) / 2).
    pub mean: f64,
    /// Counting weight (1.0 for full cycle, 0.5 for half-cycle).
    pub count: f64,
}
impl RainflowCycle {
    /// Amplitude = range / 2.
    pub fn amplitude(&self) -> f64 {
        self.range / 2.0
    }
    /// Stress ratio R = (mean − amplitude) / (mean + amplitude).
    pub fn stress_ratio(&self) -> f64 {
        let s_min = self.mean - self.amplitude();
        let s_max = self.mean + self.amplitude();
        if s_max.abs() > 1e-30 {
            s_min / s_max
        } else {
            0.0
        }
    }
}
/// Basquin (stress-life) S-N curve.
///
/// σ_a = A * N^B
///
/// where A is the fatigue strength coefficient and B is the exponent.
/// Alternatively: N = (σ_a / A)^(1/B)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct BasquinCurve {
    /// Fatigue strength coefficient A (Pa)
    pub a: f64,
    /// Fatigue strength exponent B (negative, typically -0.05 to -0.15)
    pub b_exp: f64,
    /// Endurance limit σ_e (Pa). Below this stress, infinite life is assumed.
    pub endurance_limit: f64,
}
#[allow(dead_code)]
impl BasquinCurve {
    /// Create a new Basquin S-N curve.
    pub fn new(a: f64, b_exp: f64, endurance_limit: f64) -> Self {
        Self {
            a,
            b_exp,
            endurance_limit,
        }
    }
    /// Create from two known (N, σ) data points.
    ///
    /// σ_1 = A * N_1^B and σ_2 = A * N_2^B
    /// → B = ln(σ_1/σ_2) / ln(N_1/N_2)
    /// → A = σ_1 / N_1^B
    pub fn from_two_points(
        n1: f64,
        sigma1: f64,
        n2: f64,
        sigma2: f64,
        endurance_limit: f64,
    ) -> Self {
        let b_exp = (sigma1 / sigma2).ln() / (n1 / n2).ln();
        let a = sigma1 / n1.powf(b_exp);
        Self {
            a,
            b_exp,
            endurance_limit,
        }
    }
    /// Stress amplitude at N cycles: σ_a = A * N^B
    pub fn stress_at_n(&self, n: f64) -> f64 {
        self.a * n.powf(self.b_exp)
    }
    /// Cycles to failure at stress amplitude σ_a: N = (σ_a / A)^(1/B)
    ///
    /// Returns `f64::INFINITY` if σ_a <= endurance_limit.
    pub fn cycles_to_failure(&self, stress_amplitude: f64) -> f64 {
        if stress_amplitude <= self.endurance_limit {
            return f64::INFINITY;
        }
        (stress_amplitude / self.a).powf(1.0 / self.b_exp)
    }
}
/// Coffin-Manson-Basquin strain-life fatigue parameters.
///
/// The total strain amplitude is decomposed into elastic and plastic parts:
/// - Δε_e/2 = (σ_f/E) * (2N)^b   (Basquin elastic term)
/// - Δε_p/2 = ε_f * (2N)^c        (Coffin-Manson plastic term)
/// - Δε/2   = Δε_e/2 + Δε_p/2    (total)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SNcurve {
    /// Fatigue strength coefficient σ_f (Pa)
    pub sigma_f: f64,
    /// Fatigue strength exponent b (typically -0.05 to -0.12)
    pub b: f64,
    /// Fatigue ductility coefficient ε_f (dimensionless)
    pub epsilon_f: f64,
    /// Fatigue ductility exponent c (typically -0.5 to -0.7)
    pub c: f64,
    /// Young's modulus E (Pa)
    pub e_modulus: f64,
}
impl SNcurve {
    /// Create a new Coffin-Manson-Basquin S-N curve.
    ///
    /// # Arguments
    /// * `sigma_f`   - Fatigue strength coefficient (Pa)
    /// * `b`         - Fatigue strength exponent (dimensionless, typically negative)
    /// * `epsilon_f` - Fatigue ductility coefficient (dimensionless)
    /// * `c`         - Fatigue ductility exponent (dimensionless, typically negative)
    /// * `e_modulus` - Young's modulus (Pa)
    pub fn new(sigma_f: f64, b: f64, epsilon_f: f64, c: f64, e_modulus: f64) -> Self {
        Self {
            sigma_f,
            b,
            epsilon_f,
            c,
            e_modulus,
        }
    }
    /// Elastic strain amplitude: Δε_e/2 = (σ_f/E) * (2N)^b
    ///
    /// # Arguments
    /// * `n_reversals` - Number of reversals 2N (not cycles)
    pub fn elastic_strain_amplitude(&self, n_reversals: f64) -> f64 {
        (self.sigma_f / self.e_modulus) * n_reversals.powf(self.b)
    }
    /// Plastic strain amplitude: Δε_p/2 = ε_f * (2N)^c
    ///
    /// # Arguments
    /// * `n_reversals` - Number of reversals 2N (not cycles)
    pub fn plastic_strain_amplitude(&self, n_reversals: f64) -> f64 {
        self.epsilon_f * n_reversals.powf(self.c)
    }
    /// Total strain amplitude: Δε/2 = elastic + plastic
    ///
    /// # Arguments
    /// * `n_reversals` - Number of reversals 2N (not cycles)
    pub fn total_strain_amplitude(&self, n_reversals: f64) -> f64 {
        self.elastic_strain_amplitude(n_reversals) + self.plastic_strain_amplitude(n_reversals)
    }
    /// Cycles to failure for a given strain amplitude (bisection method).
    ///
    /// Solves Δε/2 = (σ_f/E)*(2N)^b + ε_f*(2N)^c for N using bisection.
    ///
    /// # Arguments
    /// * `strain_amplitude` - Total strain amplitude Δε/2
    pub fn cycles_to_failure_strain(&self, strain_amplitude: f64) -> f64 {
        let mut lo = 2.0_f64;
        let mut hi = 2.0e12_f64;
        let f_lo = self.total_strain_amplitude(lo) - strain_amplitude;
        let f_hi = self.total_strain_amplitude(hi) - strain_amplitude;
        if f_lo <= 0.0 {
            return lo / 2.0;
        }
        if f_hi >= 0.0 {
            return hi / 2.0;
        }
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            let f_mid = self.total_strain_amplitude(mid) - strain_amplitude;
            if f_mid > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) / hi < 1.0e-12 {
                break;
            }
        }
        (lo + hi) / 2.0 / 2.0
    }
    /// Stress amplitude from Basquin equation: σ_a = σ_f * (2N)^b
    ///
    /// # Arguments
    /// * `n_cycles` - Number of cycles N (not reversals)
    pub fn stress_amplitude_from_n(&self, n_cycles: f64) -> f64 {
        let two_n = 2.0 * n_cycles;
        self.sigma_f * two_n.powf(self.b)
    }
}
/// Neuber's rule for notch analysis.
///
/// Neuber's rule relates the theoretical (elastic) stress/strain at a notch
/// to the actual local elastic-plastic stress and strain:
///
/// (Kt * σ_nom)² / E = σ_local * ε_local
///
/// i.e., the geometric mean of local stress and strain energy equals the
/// nominal elastic value scaled by Kt².
///
/// Given nominal stress σ_nom, stress concentration factor Kt, and
/// Young's modulus E, the Neuber hyperbola is:
/// σ_local * ε_local = (Kt * σ_nom)² / E
///
/// When combined with a cyclic stress-strain curve (e.g., Ramberg-Osgood)
/// the intersection gives the actual local stress and strain.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct NeuberRule {
    /// Stress concentration factor Kt.
    pub kt: f64,
    /// Young's modulus E (Pa).
    pub e_modulus: f64,
    /// Ramberg-Osgood cyclic strength coefficient K' (Pa).
    pub k_prime: f64,
    /// Ramberg-Osgood cyclic strain hardening exponent n'.
    pub n_prime: f64,
}
#[allow(dead_code)]
impl NeuberRule {
    /// Create a new Neuber rule model.
    pub fn new(kt: f64, e_modulus: f64, k_prime: f64, n_prime: f64) -> Self {
        Self {
            kt,
            e_modulus,
            k_prime,
            n_prime,
        }
    }
    /// Neuber hyperbola product: (Kt * σ_nom)² / E
    pub fn neuber_product(&self, sigma_nom: f64) -> f64 {
        (self.kt * sigma_nom).powi(2) / self.e_modulus
    }
    /// Ramberg-Osgood total strain for a given local stress.
    ///
    /// ε = σ/E + (σ/K')^(1/n')
    pub fn ramberg_osgood_strain(&self, sigma_local: f64) -> f64 {
        sigma_local / self.e_modulus + (sigma_local / self.k_prime).powf(1.0 / self.n_prime)
    }
    /// Solve for local stress using Neuber's rule + Ramberg-Osgood (bisection).
    ///
    /// Finds σ_local such that σ_local * ε(σ_local) = (Kt * σ_nom)² / E.
    pub fn local_stress(&self, sigma_nom: f64, max_iter: usize) -> f64 {
        let target = self.neuber_product(sigma_nom);
        if target <= 0.0 {
            return 0.0;
        }
        let sigma_el = self.kt * sigma_nom.abs();
        let mut lo = 0.0_f64;
        let mut hi = sigma_el * 3.0;
        for _ in 0..max_iter {
            let mid = 0.5 * (lo + hi);
            let eps_mid = self.ramberg_osgood_strain(mid);
            let prod = mid * eps_mid;
            if prod < target {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) / sigma_el.max(1.0) < 1e-12 {
                break;
            }
        }
        0.5 * (lo + hi)
    }
    /// Solve for local strain from local stress (direct Ramberg-Osgood).
    pub fn local_strain(&self, sigma_local: f64) -> f64 {
        self.ramberg_osgood_strain(sigma_local)
    }
    /// Cyclic Neuber analysis: compute local stress and strain amplitudes
    /// for a given nominal stress amplitude using Masing's hypothesis.
    ///
    /// For cyclic loading, the cyclic Ramberg-Osgood uses factor 2:
    /// Δε = Δσ/E + 2*(Δσ/2K')^(1/n')
    ///
    /// Returns (sigma_a_local, epsilon_a_local).
    pub fn cyclic_local_stress_strain(&self, sigma_nom_amp: f64, max_iter: usize) -> (f64, f64) {
        let target = self.neuber_product(sigma_nom_amp);
        if target <= 0.0 {
            return (0.0, 0.0);
        }
        let sigma_el = self.kt * sigma_nom_amp;
        let mut lo = 0.0_f64;
        let mut hi = sigma_el * 3.0;
        for _ in 0..max_iter {
            let mid = 0.5 * (lo + hi);
            let eps_mid = self.ramberg_osgood_strain(mid);
            let prod = mid * eps_mid;
            if prod < target {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) / sigma_el.max(1.0) < 1e-12 {
                break;
            }
        }
        let sigma_a = 0.5 * (lo + hi);
        (sigma_a, self.ramberg_osgood_strain(sigma_a))
    }
}
/// Coffin-Manson low-cycle fatigue model (extended variant with primed coefficients).
///
/// Total strain amplitude:
/// Δε/2 = ε_f' * (2N_f)^c
///
/// where c is the ductility exponent (typically −0.5 to −0.7).
#[allow(dead_code)]
pub struct CoffinMansonLcf {
    /// Fatigue ductility coefficient ε_f' (dimensionless)
    pub epsilon_f_prime: f64,
    /// Fatigue ductility exponent c (negative)
    pub c_exp: f64,
}
impl CoffinMansonLcf {
    /// Create a new Coffin-Manson low-cycle fatigue model.
    pub fn new(epsilon_f_prime: f64, c_exp: f64) -> Self {
        Self {
            epsilon_f_prime,
            c_exp,
        }
    }
    /// Plastic strain amplitude at a given number of reversals 2N.
    pub fn plastic_strain_amplitude(&self, n_reversals: f64) -> f64 {
        self.epsilon_f_prime * n_reversals.powf(self.c_exp)
    }
    /// Reversals to failure for a given plastic strain amplitude.
    pub fn reversals_to_failure(&self, delta_epsilon_p_half: f64) -> f64 {
        (delta_epsilon_p_half / self.epsilon_f_prime).powf(1.0 / self.c_exp)
    }
    /// Cycles to failure (= reversals / 2).
    pub fn cycles_to_failure(&self, delta_epsilon_p_half: f64) -> f64 {
        self.reversals_to_failure(delta_epsilon_p_half) / 2.0
    }
    /// Transition reversals between low- and high-cycle fatigue.
    ///
    /// At the transition, plastic strain amplitude = elastic strain amplitude.
    /// Requires the Basquin parameters (σ_f/E, b) to solve simultaneously.
    pub fn transition_reversals(&self, sigma_f: f64, e_modulus: f64, b_exp: f64) -> f64 {
        let ratio = self.epsilon_f_prime * e_modulus / sigma_f;
        ratio.powf(1.0 / (b_exp - self.c_exp))
    }
    /// Damage per cycle for Palmgren-Miner accumulation.
    pub fn damage_per_cycle(&self, delta_epsilon_p_half: f64) -> f64 {
        let n_f = self.cycles_to_failure(delta_epsilon_p_half);
        if n_f > 0.0 && n_f.is_finite() {
            1.0 / n_f
        } else {
            0.0
        }
    }
}
/// Goodman diagram defined by ultimate strength and endurance limit.
///
/// Allowable amplitude: σ_a = σ_e * (1 - σ_m / σ_ult)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GoodmanDiagramNew {
    /// Ultimate tensile strength σ_ult (Pa).
    pub sigma_ult: f64,
    /// Fully-reversed endurance limit σ_e (Pa).
    pub sigma_endurance: f64,
}
#[allow(dead_code)]
impl GoodmanDiagramNew {
    /// Create a new Goodman diagram.
    pub fn new(sigma_ult: f64, sigma_endurance: f64) -> Self {
        Self {
            sigma_ult,
            sigma_endurance,
        }
    }
    /// Allowable stress amplitude for a given mean stress.
    ///
    /// σ_a_allow = σ_e * (1 - σ_m / σ_ult)  \[clamped to 0\]
    pub fn goodman_allowed_amplitude(&self, mean_stress: f64) -> f64 {
        (self.sigma_endurance * (1.0 - mean_stress / self.sigma_ult)).max(0.0)
    }
}
/// Smith-Watson-Topper (SWT) damage parameter.
///
/// SWT = σ_max * Δε/2
///
/// Used as a fatigue damage indicator that accounts for both mean stress
/// and strain amplitude. Failure occurs when SWT equals the material
/// SWT capacity.
///
/// The SWT parameter for a Coffin-Manson-Basquin material:
/// SWT = (σ_f')² / E * (2N)^(2b) + σ_f' * ε_f' * (2N)^(b+c)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SwtParameter {
    /// Fatigue strength coefficient σ_f' (Pa)
    pub sigma_f: f64,
    /// Fatigue strength exponent b
    pub b: f64,
    /// Fatigue ductility coefficient ε_f'
    pub epsilon_f: f64,
    /// Fatigue ductility exponent c
    pub c: f64,
    /// Young's modulus E (Pa)
    pub e_modulus: f64,
}
#[allow(dead_code)]
impl SwtParameter {
    /// Create a new SWT parameter model.
    pub fn new(sigma_f: f64, b: f64, epsilon_f: f64, c: f64, e_modulus: f64) -> Self {
        Self {
            sigma_f,
            b,
            epsilon_f,
            c,
            e_modulus,
        }
    }
    /// Compute the SWT damage parameter from max stress and strain amplitude.
    ///
    /// SWT = σ_max * Δε/2
    pub fn compute(max_stress: f64, strain_amplitude: f64) -> f64 {
        if max_stress <= 0.0 {
            return 0.0;
        }
        max_stress * strain_amplitude
    }
    /// Material SWT capacity at a given number of reversals.
    ///
    /// SWT = (σ_f')²/E * (2N)^(2b) + σ_f' * ε_f' * (2N)^(b+c)
    pub fn material_swt(&self, two_n: f64) -> f64 {
        let term1 = self.sigma_f.powi(2) / self.e_modulus * two_n.powf(2.0 * self.b);
        let term2 = self.sigma_f * self.epsilon_f * two_n.powf(self.b + self.c);
        term1 + term2
    }
    /// Cycles to failure from SWT parameter (bisection).
    pub fn cycles_to_failure(&self, swt_value: f64) -> f64 {
        if swt_value <= 0.0 {
            return f64::INFINITY;
        }
        let mut lo = 2.0_f64;
        let mut hi = 2.0e12_f64;
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            let swt_mid = self.material_swt(mid);
            if swt_mid > swt_value {
                lo = mid;
            } else {
                hi = mid;
            }
            if (hi - lo) / hi < 1e-12 {
                break;
            }
        }
        (lo + hi) / 4.0
    }
}
/// Log-normal S-N scatter band model.
///
/// At a given stress amplitude, the fatigue life follows a log-normal
/// distribution with parameters (log_mean_n, log_std_n) derived from
/// a reference life N_ref and a scatter factor T_n.
///
/// T_n = N_p10 / N_p90 (ratio of 10%ile to 90%ile life).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SNScatterBand {
    /// Basquin S-N curve (median/mean life).
    pub sn: BasquinCurve,
    /// Scatter factor T_n = N_90% / N_10% (ratio of 90th to 10th percentile lives).
    pub t_n: f64,
}
#[allow(dead_code)]
impl SNScatterBand {
    /// Create a new S-N scatter band.
    ///
    /// # Arguments
    /// * `sn`  - Median Basquin S-N curve.
    /// * `t_n` - Scatter factor (ratio N_10 / N_90, typically 3..30 for metals).
    pub fn new(sn: BasquinCurve, t_n: f64) -> Self {
        Self { sn, t_n }
    }
    /// Log-normal standard deviation of log10(N).
    ///
    /// log10(T_n) = 2 * 1.281 * s_log10_n
    /// → s_log10_n = log10(T_n) / (2 * 1.281)
    pub fn log_std(&self) -> f64 {
        self.t_n.log10() / (2.0 * 1.2816)
    }
    /// Life at a given probability of survival p_s (0..1).
    ///
    /// log10(N_ps) = log10(N_median) + z_p * s_log10
    /// where z_p is the normal quantile (< 0 for p_s < 0.5).
    ///
    /// Uses the rational approximation for the normal quantile.
    pub fn life_at_survival(&self, sigma_amp: f64, p_survival: f64) -> f64 {
        if sigma_amp <= self.sn.endurance_limit {
            return f64::INFINITY;
        }
        let n_median = self.sn.cycles_to_failure(sigma_amp);
        if !n_median.is_finite() {
            return f64::INFINITY;
        }
        let log_n_med = n_median.log10();
        let s = self.log_std();
        let z = normal_quantile(1.0 - p_survival);
        10.0_f64.powf(log_n_med + z * s)
    }
    /// Probability of survival at a given life N for a given stress amplitude.
    ///
    /// P_s = Φ((log10(N) - log10(N_median)) / s_log10)
    pub fn survival_probability(&self, sigma_amp: f64, n_cycles: f64) -> f64 {
        if sigma_amp <= self.sn.endurance_limit {
            return 1.0;
        }
        let n_median = self.sn.cycles_to_failure(sigma_amp);
        if !n_median.is_finite() {
            return 1.0;
        }
        let s = self.log_std();
        if s <= 0.0 {
            return if n_cycles < n_median { 1.0 } else { 0.0 };
        }
        let z = (n_cycles.log10() - n_median.log10()) / s;
        standard_normal_cdf(z)
    }
}
/// Cycle-by-cycle Palmgren-Miner damage accumulator.
///
/// Accumulates `n / N_f` for each loading block.  A block represents one or
/// more applied cycles at a given stress amplitude.
#[allow(dead_code)]
pub struct CycleDamageAccumulator {
    /// Accumulated damage D = Σ (n_i / N_fi).
    pub damage: f64,
    /// Critical damage at failure (Miner's rule: 1.0).
    pub d_crit: f64,
    /// Log of (stress_amplitude, N_f, n_applied) for each block.
    pub history: Vec<(f64, f64, f64)>,
}
impl CycleDamageAccumulator {
    /// Create a new accumulator with given failure criterion.
    pub fn new(d_crit: f64) -> Self {
        Self {
            damage: 0.0,
            d_crit,
            history: Vec::new(),
        }
    }
    /// Apply `n_applied` cycles at `sigma_a` using the given S-N curve.
    ///
    /// Cycles below the endurance limit contribute zero damage.
    pub fn apply_cycles(&mut self, n_applied: f64, sigma_a: f64, sn: &BasquinStressLife) {
        let n_f = sn.cycles_to_failure(sigma_a);
        let d = if n_f.is_infinite() {
            0.0
        } else {
            n_applied / n_f
        };
        self.damage += d;
        self.history.push((sigma_a, n_f, n_applied));
    }
    /// Apply a rainflow-counted cycle set using a Basquin S-N curve.
    pub fn apply_rainflow(&mut self, cycles: &[RainflowCycle], sn: &BasquinStressLife) {
        for cyc in cycles {
            self.apply_cycles(cyc.count, cyc.amplitude(), sn);
        }
    }
    /// True when accumulated damage meets or exceeds `d_crit`.
    pub fn is_failed(&self) -> bool {
        self.damage >= self.d_crit
    }
    /// Remaining life fraction: (D_crit − D) / D_crit.
    pub fn remaining_life_fraction(&self) -> f64 {
        ((self.d_crit - self.damage) / self.d_crit).max(0.0)
    }
    /// Estimate remaining cycles at constant amplitude `sigma_a`.
    pub fn remaining_cycles(&self, sigma_a: f64, sn: &BasquinStressLife) -> f64 {
        let n_f = sn.cycles_to_failure(sigma_a);
        if n_f.is_infinite() {
            return f64::INFINITY;
        }
        let remaining_d = (self.d_crit - self.damage).max(0.0);
        remaining_d * n_f
    }
    /// Reset accumulator state.
    pub fn reset(&mut self) {
        self.damage = 0.0;
        self.history.clear();
    }
    /// Number of blocks applied.
    pub fn block_count(&self) -> usize {
        self.history.len()
    }
}
/// Paris law crack growth parameters.
///
/// da/dN = C * (ΔK)^m  (valid between ΔK_th and K_c)
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ParisLaw {
    /// Paris coefficient C.
    pub c: f64,
    /// Paris exponent m.
    pub m: f64,
    /// Threshold stress intensity factor range ΔK_th (below this, no growth).
    pub delta_k_threshold: f64,
    /// Fracture toughness K_c (above this, fast fracture).
    pub k_fracture: f64,
}
#[allow(dead_code)]
impl ParisLaw {
    /// Create a new Paris law model.
    pub fn new(c: f64, m: f64, delta_k_threshold: f64, k_fracture: f64) -> Self {
        Self {
            c,
            m,
            delta_k_threshold,
            k_fracture,
        }
    }
    /// Crack growth rate da/dN for a given ΔK.
    ///
    /// Returns 0.0 if ΔK < ΔK_th, returns f64::INFINITY if ΔK ≥ K_c.
    pub fn crack_growth_rate(&self, delta_k: f64) -> f64 {
        if delta_k < self.delta_k_threshold {
            return 0.0;
        }
        if delta_k >= self.k_fracture {
            return f64::INFINITY;
        }
        self.c * delta_k.powf(self.m)
    }
    /// Integrate crack growth from initial crack size a_0 to final a_f.
    ///
    /// Uses simple forward Euler integration:
    /// a_new = a + (da/dN) * Δcycles
    ///
    /// Returns the number of cycles to grow from a_0 to a_f,
    /// or `None` if fast fracture occurs.
    pub fn cycles_to_grow(
        &self,
        a_initial: f64,
        a_final: f64,
        delta_k_fn: &impl Fn(f64) -> f64,
        da_step: f64,
    ) -> Option<f64> {
        let mut a = a_initial;
        let mut n_cycles = 0.0_f64;
        while a < a_final {
            let delta_k = delta_k_fn(a);
            let da_dn = self.crack_growth_rate(delta_k);
            if da_dn.is_infinite() || da_dn <= 0.0 {
                return if da_dn.is_infinite() {
                    None
                } else {
                    Some(f64::INFINITY)
                };
            }
            let dn = da_step / da_dn;
            n_cycles += dn;
            a += da_step;
        }
        Some(n_cycles)
    }
    /// Stress intensity factor range for a center crack in an infinite plate.
    ///
    /// ΔK = Δσ * √(π * a) * Y
    /// where Y is the geometry factor (usually ≈ 1.0 for simple cases).
    pub fn delta_k_center_crack(delta_sigma: f64, a: f64, y: f64) -> f64 {
        delta_sigma * (std::f64::consts::PI * a).sqrt() * y
    }
    /// NASGRO modified Paris law including crack closure.
    ///
    /// da/dN = C * ((1-f)/(1-R))^m * (ΔK)^m
    /// where f is the crack opening function and R is the stress ratio.
    pub fn nasgro_rate(&self, delta_k: f64, r_ratio: f64, crack_opening_f: f64) -> f64 {
        let effective_dk = delta_k * (1.0 - crack_opening_f) / (1.0 - r_ratio).max(0.01);
        self.crack_growth_rate(effective_dk)
    }
}
/// Modified Goodman diagram for mean stress correction.
///
/// Accounts for the effect of non-zero mean stress on fatigue life.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GoodmanDiagram {
    /// Ultimate tensile strength σ_u (Pa)
    pub ultimate_strength: f64,
    /// Yield strength σ_ys (Pa)
    pub yield_strength: f64,
}
impl GoodmanDiagram {
    /// Create a new Goodman diagram.
    ///
    /// # Arguments
    /// * `uts` - Ultimate tensile strength (Pa)
    /// * `ys`  - Yield strength (Pa)
    pub fn new(uts: f64, ys: f64) -> Self {
        Self {
            ultimate_strength: uts,
            yield_strength: ys,
        }
    }
    /// Allowable stress amplitude via modified Goodman criterion.
    ///
    /// σ_a = σ_e * (1 - σ_m/σ_u) / FS
    ///
    /// where σ_e is the fully-reversed endurance limit.
    ///
    /// # Arguments
    /// * `mean_stress`      - Mean stress σ_m (Pa)
    /// * `factor_of_safety` - Safety factor FS (> 0)
    pub fn allowable_amplitude(
        &self,
        mean_stress: f64,
        endurance: f64,
        factor_of_safety: f64,
    ) -> f64 {
        endurance * (1.0 - mean_stress / self.ultimate_strength) / factor_of_safety
    }
    /// Allowable stress amplitude via Soderberg criterion (uses yield strength).
    ///
    /// σ_a = σ_e * (1 - σ_m/σ_ys)
    ///
    /// # Arguments
    /// * `mean_stress` - Mean stress σ_m (Pa)
    /// * `endurance`   - Fully-reversed endurance limit σ_e (Pa)
    pub fn soderberg_amplitude(&self, mean_stress: f64, endurance: f64) -> f64 {
        endurance * (1.0 - mean_stress / self.yield_strength)
    }
    /// Check whether a (mean_stress, amplitude) point is safe by Goodman criterion.
    ///
    /// Safe if: σ_a/σ_e + σ_m/σ_u ≤ 1
    ///
    /// # Arguments
    /// * `mean_stress` - Mean stress σ_m (Pa)
    /// * `amplitude`   - Stress amplitude σ_a (Pa)
    /// * `endurance`   - Fully-reversed endurance limit σ_e (Pa)
    pub fn is_safe(&self, mean_stress: f64, amplitude: f64, endurance: f64) -> bool {
        amplitude / endurance + mean_stress / self.ultimate_strength <= 1.0
    }
    /// Compute the Goodman equivalent fully-reversed stress amplitude.
    ///
    /// σ_ar = σ_a / (1 - σ_m / σ_u)
    ///
    /// This is the equivalent fully-reversed stress amplitude that produces
    /// the same fatigue damage as the given (σ_a, σ_m) combination.
    #[allow(dead_code)]
    pub fn equivalent_amplitude(&self, mean_stress: f64, amplitude: f64) -> f64 {
        let denom = 1.0 - mean_stress / self.ultimate_strength;
        if denom <= 0.0 {
            return f64::INFINITY;
        }
        amplitude / denom
    }
    /// Factor of safety for a given (mean, amplitude) point.
    ///
    /// FS = 1 / (σ_a/σ_e + σ_m/σ_u)
    #[allow(dead_code)]
    pub fn factor_of_safety(&self, mean_stress: f64, amplitude: f64, endurance: f64) -> f64 {
        let ratio = amplitude / endurance + mean_stress / self.ultimate_strength;
        if ratio <= 0.0 {
            return f64::INFINITY;
        }
        1.0 / ratio
    }
}
/// Damage tolerance analysis result.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DamageToleranceResult {
    /// Initial crack size.
    pub a_initial: f64,
    /// Critical crack size (fracture).
    pub a_critical: f64,
    /// Inspection interval (cycles).
    pub inspection_interval: f64,
    /// Cycles to critical crack size.
    pub cycles_to_critical: Option<f64>,
    /// Safety factor on life.
    pub safety_factor: f64,
}
#[allow(dead_code)]
impl DamageToleranceResult {
    /// Whether the component is safe (cycles_to_critical > inspection_interval * safety_factor).
    pub fn is_safe(&self) -> bool {
        if let Some(cycles) = self.cycles_to_critical {
            cycles > self.inspection_interval * self.safety_factor
        } else {
            false
        }
    }
    /// Remaining useful life (cycles) given current cycle count.
    pub fn remaining_life(&self, current_cycles: f64) -> f64 {
        match self.cycles_to_critical {
            Some(total) => (total - current_cycles).max(0.0),
            None => 0.0,
        }
    }
}
/// Gerber diagram defined by ultimate strength and endurance limit.
///
/// Allowable amplitude: σ_a = σ_e * (1 - (σ_m / σ_ult)^2)
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct GerberDiagram {
    /// Ultimate tensile strength σ_ult (Pa).
    pub sigma_ult: f64,
    /// Fully-reversed endurance limit σ_e (Pa).
    pub sigma_endurance: f64,
}
#[allow(dead_code)]
impl GerberDiagram {
    /// Create a new Gerber diagram.
    pub fn new(sigma_ult: f64, sigma_endurance: f64) -> Self {
        Self {
            sigma_ult,
            sigma_endurance,
        }
    }
    /// Allowable stress amplitude for a given mean stress (Gerber parabola).
    ///
    /// σ_a_allow = σ_e * (1 - (σ_m / σ_ult)^2)  \[clamped to 0\]
    pub fn gerber_allowed_amplitude(&self, mean_stress: f64) -> f64 {
        let ratio = mean_stress / self.sigma_ult;
        (self.sigma_endurance * (1.0 - ratio * ratio)).max(0.0)
    }
}
