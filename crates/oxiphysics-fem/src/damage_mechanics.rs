// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Continuum damage mechanics for FEM: CDM, Lemaitre, Gurson-Tvergaard-Needleman.
//!
//! This module implements:
//! - [`DamageVariable`] - scalar damage D ∈ \[0,1\]
//! - [`LemaitreDamage`] - Lemaitre CDM model
//! - [`GursonModel`] - Gurson-Tvergaard-Needleman porous plasticity
//! - [`MazarsModel`] - Mazars damage model for concrete
//! - [`CoupledDamagePlasticity`] - Chaboche coupled model
//! - [`WeibullModel`] - probabilistic failure (Weibull statistics)
//! - [`FatigueLifeModel`] - Basquin + Coffin-Manson combined
//! - [`DamageLocalization`] - damage zone width and non-local averaging

/// A scalar damage variable D ∈ \[0, 1\].
///
/// D = 0 means intact material, D = 1 means complete failure.
#[derive(Debug, Clone)]
pub struct DamageVariable {
    /// Current damage value, clamped to \[0, 1\].
    pub d: f64,
}

impl DamageVariable {
    /// Creates a new undamaged variable (D = 0).
    pub fn new() -> Self {
        Self { d: 0.0 }
    }

    /// Creates a damage variable with the given initial value.
    ///
    /// # Panics
    /// Panics if `d` is not in \[0, 1\].
    pub fn with_damage(d: f64) -> Self {
        assert!((0.0..=1.0).contains(&d), "damage must be in [0, 1]");
        Self { d }
    }

    /// Computes effective stress: σ_eff = σ / (1 - D).
    ///
    /// Returns `f64::INFINITY` when D = 1 (total failure).
    pub fn effective_stress(&self, nominal_stress: f64) -> f64 {
        let denom = 1.0 - self.d;
        if denom.abs() < 1e-15 {
            f64::INFINITY
        } else {
            nominal_stress / denom
        }
    }

    /// Increases damage by `delta`, clamped to 1.0.
    pub fn accumulate(&mut self, delta: f64) {
        self.d = (self.d + delta).min(1.0);
    }

    /// Returns true if the material has completely failed (D >= 1).
    pub fn is_failed(&self) -> bool {
        self.d >= 1.0
    }
}

impl Default for DamageVariable {
    fn default() -> Self {
        Self::new()
    }
}

/// Lemaitre continuum damage mechanics model.
///
/// Damage evolution law:
/// dD/dN = (Y / S)^s
/// where Y is the strain energy release rate, S is a material constant, and s is the exponent.
#[derive(Debug, Clone)]
pub struct LemaitreDamage {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Damage threshold Y_0 (strain energy release rate threshold) \[Pa\].
    pub y_threshold: f64,
    /// Material parameter S in the damage evolution law \[Pa\].
    pub s_param: f64,
    /// Exponent s in the damage evolution law (dimensionless).
    pub s_exponent: f64,
    /// Current damage value D ∈ \[0, 1\].
    pub damage: f64,
}

impl LemaitreDamage {
    /// Creates a new Lemaitre damage model.
    pub fn new(young_modulus: f64, y_threshold: f64, s_param: f64, s_exponent: f64) -> Self {
        Self {
            young_modulus,
            y_threshold,
            s_param,
            s_exponent,
            damage: 0.0,
        }
    }

    /// Computes the strain energy release rate Y = σ²/(2E(1-D)²) for 1D case.
    pub fn strain_energy_release_rate(&self, stress: f64) -> f64 {
        let denom = 2.0 * self.young_modulus * (1.0 - self.damage).powi(2);
        if denom < 1e-15 {
            f64::INFINITY
        } else {
            stress * stress / denom
        }
    }

    /// Computes the damage increment for the given stress state.
    ///
    /// Returns 0 if Y < Y_0 (below threshold).
    pub fn damage_increment(&self, stress: f64) -> f64 {
        let y = self.strain_energy_release_rate(stress);
        if y <= self.y_threshold {
            return 0.0;
        }
        (y / self.s_param).powf(self.s_exponent)
    }

    /// Updates the damage variable by one load cycle.
    pub fn update(&mut self, stress: f64) {
        let delta = self.damage_increment(stress);
        self.damage = (self.damage + delta).min(1.0);
    }

    /// Computes the effective strain: ε̃ = ε / (1 - D).
    pub fn effective_strain(&self, nominal_strain: f64) -> f64 {
        let denom = 1.0 - self.damage;
        if denom < 1e-15 {
            f64::INFINITY
        } else {
            nominal_strain / denom
        }
    }

    /// Computes the effective stress: σ̃ = E * ε̃.
    pub fn effective_stress_from_strain(&self, nominal_strain: f64) -> f64 {
        self.young_modulus * self.effective_strain(nominal_strain)
    }

    /// Returns the secant modulus E*(1-D).
    pub fn secant_modulus(&self) -> f64 {
        self.young_modulus * (1.0 - self.damage)
    }
}

/// Gurson-Tvergaard-Needleman (GTN) porous plasticity model.
///
/// Yield function:
/// Φ = (σ_eq/σ_y)² + 2*q1*f*cosh(3*q2*σ_h/(2*σ_y)) - (1 + q3*f²) = 0
#[derive(Debug, Clone)]
pub struct GursonModel {
    /// Initial void volume fraction.
    pub void_fraction: f64,
    /// Critical void fraction f_c (onset of coalescence).
    pub f_critical: f64,
    /// Failure void fraction f_F.
    pub f_failure: f64,
    /// GTN parameter q1 (≈ 1.5 for metals).
    pub q1: f64,
    /// GTN parameter q2 (≈ 1.0 for metals).
    pub q2: f64,
    /// GTN parameter q3 (= q1² for metals).
    pub q3: f64,
    /// Yield stress of the matrix material \[Pa\].
    pub yield_stress: f64,
    /// Nucleation strain ε_N.
    pub nucleation_strain: f64,
    /// Nucleation standard deviation s_N.
    pub nucleation_std: f64,
    /// Nucleation volume fraction f_N.
    pub nucleation_fraction: f64,
}

impl GursonModel {
    /// Creates a new GTN model with standard Tvergaard parameters.
    pub fn new(
        void_fraction: f64,
        f_critical: f64,
        f_failure: f64,
        yield_stress: f64,
        nucleation_strain: f64,
        nucleation_std: f64,
        nucleation_fraction: f64,
    ) -> Self {
        Self {
            void_fraction,
            f_critical,
            f_failure,
            q1: 1.5,
            q2: 1.0,
            q3: 2.25, // q1²
            yield_stress,
            nucleation_strain,
            nucleation_std,
            nucleation_fraction,
        }
    }

    /// Evaluates the GTN yield function.
    ///
    /// Returns 0 when on the yield surface, < 0 in elastic domain, > 0 if yielded.
    pub fn yield_function(&self, sigma_eq: f64, sigma_h: f64) -> f64 {
        let f = self.effective_void_fraction();
        let ratio_eq = sigma_eq / self.yield_stress;
        let arg = 3.0 * self.q2 * sigma_h / (2.0 * self.yield_stress);
        ratio_eq * ratio_eq + 2.0 * self.q1 * f * arg.cosh() - (1.0 + self.q3 * f * f)
    }

    /// Computes the effective void fraction accounting for coalescence.
    pub fn effective_void_fraction(&self) -> f64 {
        let f = self.void_fraction;
        if f <= self.f_critical {
            f
        } else if f >= self.f_failure {
            1.0 / self.q1
        } else {
            let f_star_c = self.f_critical;
            let f_u = 1.0 / self.q1;
            f_star_c + (f_u - f_star_c) * (f - self.f_critical) / (self.f_failure - self.f_critical)
        }
    }

    /// Void nucleation rate (Chu-Needleman model).
    ///
    /// df_nucleation/dε_p = (f_N / (s_N * √(2π))) * exp(-0.5*((ε_p - ε_N)/s_N)²)
    pub fn nucleation_rate(&self, plastic_strain: f64) -> f64 {
        let two_pi: f64 = std::f64::consts::TAU;
        let denom = self.nucleation_std * two_pi.sqrt();
        let exponent =
            -0.5 * ((plastic_strain - self.nucleation_strain) / self.nucleation_std).powi(2);
        self.nucleation_fraction / denom * exponent.exp()
    }

    /// Void growth rate: df_growth/dε_p = (1-f)*tr(dε_p).
    ///
    /// `volumetric_plastic_strain_rate` is tr(ε̇_p).
    pub fn growth_rate(&self, volumetric_plastic_strain_rate: f64) -> f64 {
        (1.0 - self.void_fraction) * volumetric_plastic_strain_rate
    }

    /// Total void evolution rate.
    pub fn total_void_rate(&self, plastic_strain: f64, volumetric_plastic_strain_rate: f64) -> f64 {
        self.nucleation_rate(plastic_strain) + self.growth_rate(volumetric_plastic_strain_rate)
    }

    /// Updates the void fraction by one step.
    pub fn update_void_fraction(
        &mut self,
        plastic_strain: f64,
        volumetric_plastic_strain_rate: f64,
        dt: f64,
    ) {
        let rate = self.total_void_rate(plastic_strain, volumetric_plastic_strain_rate);
        self.void_fraction = (self.void_fraction + rate * dt).clamp(0.0, 1.0);
    }

    /// Returns true if the material has failed (void fraction >= f_failure).
    pub fn is_failed(&self) -> bool {
        self.void_fraction >= self.f_failure
    }

    /// Triaxiality ratio: T = σ_h / σ_eq.
    pub fn triaxiality(sigma_h: f64, sigma_eq: f64) -> f64 {
        if sigma_eq.abs() < 1e-15 {
            0.0
        } else {
            sigma_h / sigma_eq
        }
    }
}

/// Mazars damage model for concrete.
///
/// Uses tension/compression split with equivalent strain:
/// ε̃ = √(Σ<ε_i>²) where <·> is the Macaulay bracket.
#[derive(Debug, Clone)]
pub struct MazarsModel {
    /// Damage threshold ε_0 (strain below which damage does not occur).
    pub epsilon_0: f64,
    /// Tension damage parameter A_t.
    pub a_tension: f64,
    /// Tension damage parameter B_t.
    pub b_tension: f64,
    /// Compression damage parameter A_c.
    pub a_compression: f64,
    /// Compression damage parameter B_c.
    pub b_compression: f64,
    /// Current damage variable D ∈ \[0, 1\].
    pub damage: f64,
    /// Maximum historical equivalent strain κ.
    pub kappa: f64,
}

impl MazarsModel {
    /// Creates a new Mazars model with typical concrete parameters.
    pub fn new(
        epsilon_0: f64,
        a_tension: f64,
        b_tension: f64,
        a_compression: f64,
        b_compression: f64,
    ) -> Self {
        Self {
            epsilon_0,
            a_tension,
            b_tension,
            a_compression,
            b_compression,
            damage: 0.0,
            kappa: epsilon_0,
        }
    }

    /// Computes the equivalent strain from principal strains.
    ///
    /// ε̃ = √(Σ max(ε_i, 0)²)
    pub fn equivalent_strain(principal_strains: &[f64; 3]) -> f64 {
        let sum_sq: f64 = principal_strains
            .iter()
            .map(|&e| {
                let pos = e.max(0.0);
                pos * pos
            })
            .sum();
        sum_sq.sqrt()
    }

    /// Computes tension damage D_t.
    pub fn tension_damage(&self, kappa: f64) -> f64 {
        if kappa <= self.epsilon_0 {
            return 0.0;
        }
        1.0 - self.epsilon_0 * (1.0 - self.a_tension) / kappa
            - self.a_tension * (-self.b_tension * (kappa - self.epsilon_0)).exp()
    }

    /// Computes compression damage D_c.
    pub fn compression_damage(&self, kappa: f64) -> f64 {
        if kappa <= self.epsilon_0 {
            return 0.0;
        }
        1.0 - self.epsilon_0 * (1.0 - self.a_compression) / kappa
            - self.a_compression * (-self.b_compression * (kappa - self.epsilon_0)).exp()
    }

    /// Combines tension and compression damage.
    ///
    /// `alpha_t` is the tension weight (0 to 1).
    pub fn combined_damage(&self, kappa: f64, alpha_t: f64) -> f64 {
        let d_t = self.tension_damage(kappa);
        let d_c = self.compression_damage(kappa);
        let alpha_c = 1.0 - alpha_t;
        (alpha_t * d_t + alpha_c * d_c).clamp(0.0, 1.0)
    }

    /// Updates the model given new principal strains.
    ///
    /// Returns the new damage value.
    pub fn update(&mut self, principal_strains: &[f64; 3], alpha_t: f64) -> f64 {
        let epsilon_tilde = Self::equivalent_strain(principal_strains);
        if epsilon_tilde > self.kappa {
            self.kappa = epsilon_tilde;
        }
        self.damage = self.combined_damage(self.kappa, alpha_t);
        self.damage
    }
}

/// Chaboche coupled damage-plasticity model.
///
/// Incorporates kinematic (back-stress) and isotropic hardening,
/// with damage coupled to the plastic strain accumulation.
#[derive(Debug, Clone)]
pub struct CoupledDamagePlasticity {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Initial yield stress \[Pa\].
    pub sigma_y0: f64,
    /// Isotropic hardening modulus Q \[Pa\].
    pub q_hardening: f64,
    /// Isotropic hardening rate b.
    pub b_rate: f64,
    /// Kinematic hardening modulus C \[Pa\].
    pub c_kinematic: f64,
    /// Kinematic hardening recall coefficient γ.
    pub gamma_recall: f64,
    /// Damage parameter S (energy release rate scale) \[Pa\].
    pub s_damage: f64,
    /// Damage exponent s.
    pub s_exponent: f64,
    /// Damage threshold Y_0 \[Pa\].
    pub y_threshold: f64,
    /// Current accumulated plastic strain p.
    pub plastic_strain: f64,
    /// Current isotropic hardening variable R.
    pub r_isotropic: f64,
    /// Current back-stress X \[Pa\].
    pub back_stress: f64,
    /// Current damage D ∈ \[0, 1\].
    pub damage: f64,
}

impl CoupledDamagePlasticity {
    /// Creates a new coupled damage-plasticity model.
    pub fn new(
        young_modulus: f64,
        sigma_y0: f64,
        q_hardening: f64,
        b_rate: f64,
        c_kinematic: f64,
        gamma_recall: f64,
        s_damage: f64,
        s_exponent: f64,
        y_threshold: f64,
    ) -> Self {
        Self {
            young_modulus,
            sigma_y0,
            q_hardening,
            b_rate,
            c_kinematic,
            gamma_recall,
            s_damage,
            s_exponent,
            y_threshold,
            plastic_strain: 0.0,
            r_isotropic: 0.0,
            back_stress: 0.0,
            damage: 0.0,
        }
    }

    /// Current yield stress including isotropic hardening.
    pub fn current_yield_stress(&self) -> f64 {
        self.sigma_y0 + self.r_isotropic
    }

    /// Strain energy release rate Y = σ²/(2E(1-D)²).
    pub fn strain_energy_release_rate(&self, stress: f64) -> f64 {
        let denom = 2.0 * self.young_modulus * (1.0 - self.damage).powi(2);
        if denom < 1e-15 {
            f64::INFINITY
        } else {
            stress * stress / denom
        }
    }

    /// Performs one return-mapping step (1D simplified).
    ///
    /// Returns (updated_stress, delta_p).
    pub fn return_mapping(&mut self, trial_stress: f64) -> (f64, f64) {
        let eff_modulus = self.young_modulus * (1.0 - self.damage);
        let overstress = (trial_stress - self.back_stress).abs() - self.current_yield_stress();
        if overstress <= 0.0 {
            return (trial_stress, 0.0);
        }
        let denominator = eff_modulus + self.q_hardening * self.b_rate + self.c_kinematic;
        let delta_p = overstress / denominator;
        let sign = if trial_stress - self.back_stress >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let stress = trial_stress - sign * eff_modulus * delta_p;

        // Update hardening variables
        self.plastic_strain += delta_p;
        self.r_isotropic += self.q_hardening * (1.0 - (-self.b_rate * delta_p).exp());
        self.back_stress +=
            self.c_kinematic * delta_p * sign - self.gamma_recall * self.back_stress * delta_p;

        // Update damage
        let y = self.strain_energy_release_rate(stress);
        if y > self.y_threshold {
            let dd = ((y / self.s_damage).powf(self.s_exponent)) * delta_p;
            self.damage = (self.damage + dd).min(1.0);
        }

        (stress, delta_p)
    }
}

/// Weibull probabilistic failure model.
///
/// Failure probability: P_f = 1 - exp(-(σ/σ_0)^m)
#[derive(Debug, Clone)]
pub struct WeibullModel {
    /// Scale parameter σ_0 (characteristic strength) \[Pa\].
    pub sigma_0: f64,
    /// Shape parameter m (Weibull modulus).
    pub m_modulus: f64,
    /// Reference volume V_0 \[m³\] for volume-scaled probability.
    pub v_ref: f64,
}

impl WeibullModel {
    /// Creates a new Weibull model.
    pub fn new(sigma_0: f64, m_modulus: f64, v_ref: f64) -> Self {
        Self {
            sigma_0,
            m_modulus,
            v_ref,
        }
    }

    /// Computes the failure probability P_f = 1 - exp(-(σ/σ_0)^m).
    pub fn failure_probability(&self, stress: f64) -> f64 {
        if stress <= 0.0 {
            return 0.0;
        }
        let arg = (stress / self.sigma_0).powf(self.m_modulus);
        1.0 - (-arg).exp()
    }

    /// Computes the survival probability S = 1 - P_f.
    pub fn survival_probability(&self, stress: f64) -> f64 {
        1.0 - self.failure_probability(stress)
    }

    /// Volume-scaled failure probability for a volume V.
    ///
    /// P_f(V) = 1 - exp(-V/V_0 * (σ/σ_0)^m)
    pub fn volume_scaled_failure_probability(&self, stress: f64, volume: f64) -> f64 {
        if stress <= 0.0 {
            return 0.0;
        }
        let arg = (volume / self.v_ref) * (stress / self.sigma_0).powf(self.m_modulus);
        1.0 - (-arg).exp()
    }

    /// Returns the stress at which failure probability equals `p`.
    ///
    /// Inverse CDF: σ = σ_0 * (-ln(1-p))^(1/m)
    pub fn quantile(&self, p: f64) -> f64 {
        assert!(p > 0.0 && p < 1.0, "probability must be in (0, 1)");
        self.sigma_0 * (-((1.0 - p).ln())).powf(1.0 / self.m_modulus)
    }

    /// Mean strength: σ_mean = σ_0 * Γ(1 + 1/m).
    ///
    /// Uses Lanczos approximation for the gamma function.
    pub fn mean_strength(&self) -> f64 {
        self.sigma_0 * gamma_function(1.0 + 1.0 / self.m_modulus)
    }
}

/// Computes Γ(x) using Lanczos approximation (valid for x > 0.5).
fn gamma_function(x: f64) -> f64 {
    // Lanczos coefficients g=7, n=9
    let coeffs = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * gamma_function(1.0 - x))
    } else {
        let z = x - 1.0;
        let mut a = coeffs[0];
        let t = z + 8.5; // g + 0.5
        for (i, &c) in coeffs[1..].iter().enumerate() {
            a += c / (z + (i + 1) as f64);
        }
        let two_pi: f64 = std::f64::consts::TAU;
        (two_pi).sqrt() * t.powf(z + 0.5) * (-t).exp() * a
    }
}

/// Fatigue life model combining Basquin and Coffin-Manson laws.
///
/// Total strain amplitude: Δε/2 = σ_f'/E * (2N)^b + ε_f' * (2N)^c
#[derive(Debug, Clone)]
pub struct FatigueLifeModel {
    /// Young's modulus E \[Pa\].
    pub young_modulus: f64,
    /// Fatigue strength coefficient σ_f' \[Pa\].
    pub sigma_f_prime: f64,
    /// Fatigue strength exponent b (Basquin, typically -0.05 to -0.12).
    pub b_exponent: f64,
    /// Fatigue ductility coefficient ε_f'.
    pub epsilon_f_prime: f64,
    /// Fatigue ductility exponent c (Coffin-Manson, typically -0.5 to -0.7).
    pub c_exponent: f64,
    /// Current accumulated damage (cycle fraction).
    pub miner_damage: f64,
}

impl FatigueLifeModel {
    /// Creates a new fatigue life model.
    pub fn new(
        young_modulus: f64,
        sigma_f_prime: f64,
        b_exponent: f64,
        epsilon_f_prime: f64,
        c_exponent: f64,
    ) -> Self {
        Self {
            young_modulus,
            sigma_f_prime,
            b_exponent,
            epsilon_f_prime,
            c_exponent,
            miner_damage: 0.0,
        }
    }

    /// Computes total strain amplitude for given life 2N (reversals).
    ///
    /// Δε/2 = σ_f'/E * (2N)^b + ε_f' * (2N)^c
    pub fn strain_amplitude(&self, reversals_2n: f64) -> f64 {
        let elastic = self.sigma_f_prime / self.young_modulus * reversals_2n.powf(self.b_exponent);
        let plastic = self.epsilon_f_prime * reversals_2n.powf(self.c_exponent);
        elastic + plastic
    }

    /// Computes the stress amplitude from Basquin's law: σ_a = σ_f' * (2N)^b.
    pub fn stress_amplitude_basquin(&self, reversals_2n: f64) -> f64 {
        self.sigma_f_prime * reversals_2n.powf(self.b_exponent)
    }

    /// Estimates life in reversals (2N) given a strain amplitude Δε/2.
    ///
    /// Uses Newton-Raphson iteration.
    pub fn life_from_strain(&self, strain_amplitude: f64) -> f64 {
        let mut two_n = 1000.0f64;
        for _ in 0..100 {
            let f = self.strain_amplitude(two_n) - strain_amplitude;
            let df = self.sigma_f_prime / self.young_modulus
                * self.b_exponent
                * two_n.powf(self.b_exponent - 1.0)
                + self.epsilon_f_prime * self.c_exponent * two_n.powf(self.c_exponent - 1.0);
            if df.abs() < 1e-30 {
                break;
            }
            let two_n_new = (two_n - f / df).max(1.0);
            if (two_n_new - two_n).abs() / two_n < 1e-10 {
                two_n = two_n_new;
                break;
            }
            two_n = two_n_new;
        }
        two_n
    }

    /// Applies Miner's rule for one block of `n_cycles` at a given strain amplitude.
    ///
    /// Returns true if cumulative damage >= 1.0 (failure).
    pub fn apply_miner_cycle(&mut self, strain_amplitude: f64, n_cycles: f64) -> bool {
        let life_reversals = self.life_from_strain(strain_amplitude);
        let life_cycles = life_reversals / 2.0;
        if life_cycles > 0.0 {
            self.miner_damage += n_cycles / life_cycles;
        }
        self.miner_damage >= 1.0
    }
}

/// Non-local damage localization model.
///
/// Regularizes damage to prevent mesh-dependent solutions.
#[derive(Debug, Clone)]
pub struct DamageLocalization {
    /// Non-local averaging radius R (internal length scale) \[m\].
    pub radius: f64,
    /// Regularization type (0 = Gaussian, 1 = linear cone).
    pub regularization_type: u8,
    /// Fracture energy G_f \[J/m²\].
    pub fracture_energy: f64,
    /// Tensile strength f_t \[Pa\].
    pub tensile_strength: f64,
}

impl DamageLocalization {
    /// Creates a new damage localization regularization.
    pub fn new(radius: f64, fracture_energy: f64, tensile_strength: f64) -> Self {
        Self {
            radius,
            regularization_type: 0,
            fracture_energy,
            tensile_strength,
        }
    }

    /// Characteristic length l_c = G_f * E / f_t².
    pub fn characteristic_length(&self, young_modulus: f64) -> f64 {
        self.fracture_energy * young_modulus / (self.tensile_strength * self.tensile_strength)
    }

    /// Gaussian weight function W(r) = exp(-r²/(2*R²)) / normalization.
    pub fn gaussian_weight(&self, r: f64) -> f64 {
        let two_r_sq = 2.0 * self.radius * self.radius;
        (-r * r / two_r_sq).exp()
    }

    /// Linear cone weight W(r) = max(1 - r/R, 0).
    pub fn cone_weight(&self, r: f64) -> f64 {
        (1.0 - r / self.radius).max(0.0)
    }

    /// Computes the non-local equivalent strain at position x from neighboring points.
    ///
    /// `positions` and `local_strains` must have the same length.
    pub fn nonlocal_strain(&self, x: f64, positions: &[f64], local_strains: &[f64]) -> f64 {
        let mut weight_sum = 0.0;
        let mut weighted_strain = 0.0;
        for (&pos, &strain) in positions.iter().zip(local_strains.iter()) {
            let r = (x - pos).abs();
            let w = match self.regularization_type {
                0 => self.gaussian_weight(r),
                _ => self.cone_weight(r),
            };
            weight_sum += w;
            weighted_strain += w * strain;
        }
        if weight_sum < 1e-15 {
            0.0
        } else {
            weighted_strain / weight_sum
        }
    }

    /// Width of the localization band: w_b = π * R (Gaussian).
    pub fn band_width(&self) -> f64 {
        std::f64::consts::PI * self.radius
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // --- DamageVariable tests ---

    #[test]
    fn test_damage_zero_effective_stress_equals_nominal() {
        let dv = DamageVariable::new();
        assert!((dv.effective_stress(100.0) - 100.0).abs() < 1e-12);
    }

    #[test]
    fn test_damage_one_effective_stress_is_infinity() {
        let dv = DamageVariable::with_damage(1.0);
        assert!(dv.effective_stress(100.0).is_infinite());
    }

    #[test]
    fn test_damage_half_doubles_stress() {
        let dv = DamageVariable::with_damage(0.5);
        let eff = dv.effective_stress(50.0);
        assert!((eff - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_damage_accumulate_clamped() {
        let mut dv = DamageVariable::with_damage(0.9);
        dv.accumulate(0.2);
        assert!((dv.d - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_damage_is_failed() {
        let mut dv = DamageVariable::new();
        assert!(!dv.is_failed());
        dv.accumulate(1.0);
        assert!(dv.is_failed());
    }

    #[test]
    fn test_damage_default() {
        let dv = DamageVariable::default();
        assert!((dv.d - 0.0).abs() < 1e-12);
    }

    // --- LemaitreDamage tests ---

    #[test]
    fn test_lemaitre_below_threshold_no_damage() {
        let model = LemaitreDamage::new(200e9, 1e4, 1e5, 2.0);
        // Y = σ²/(2E) < Y_0 for small stress
        let delta = model.damage_increment(1e6);
        // Y = (1e6)^2 / (2 * 200e9) = 2500 < 1e4
        assert!((delta - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_lemaitre_above_threshold_positive_damage() {
        let model = LemaitreDamage::new(200e9, 100.0, 1e5, 1.0);
        // Y = (1e8)^2 / (2 * 200e9 * 1) = 25000 >> 100
        let delta = model.damage_increment(1e8);
        assert!(delta > 0.0);
    }

    #[test]
    fn test_lemaitre_damage_increases_monotonically() {
        let mut model = LemaitreDamage::new(70e9, 100.0, 1e4, 1.0);
        let mut prev = 0.0;
        for i in 1..20 {
            model.update(i as f64 * 1e7);
            assert!(model.damage >= prev);
            prev = model.damage;
        }
    }

    #[test]
    fn test_lemaitre_effective_strain_undamaged() {
        let model = LemaitreDamage::new(200e9, 1e4, 1e5, 2.0);
        assert!((model.effective_strain(0.001) - 0.001).abs() < 1e-12);
    }

    #[test]
    fn test_lemaitre_secant_modulus() {
        let mut model = LemaitreDamage::new(200e9, 100.0, 1e4, 1.0);
        model.damage = 0.5;
        assert!((model.secant_modulus() - 100e9).abs() < 1.0);
    }

    #[test]
    fn test_lemaitre_damage_clamped_at_one() {
        let mut model = LemaitreDamage::new(200e9, 1.0, 1.0, 0.5);
        for _ in 0..1000 {
            model.update(1e9);
        }
        assert!(model.damage <= 1.0);
    }

    // --- GursonModel tests ---

    #[test]
    fn test_gtn_yield_function_zero_void_von_mises() {
        // When f=0, GTN reduces to von Mises: Φ = (σ_eq/σ_y)² - 1
        let mut model = GursonModel::new(0.0, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        model.void_fraction = 0.0;
        let phi = model.yield_function(500e6, 0.0);
        // (1)² + 0 - 1 = 0
        assert!(phi.abs() < 1e-10);
    }

    #[test]
    fn test_gtn_yield_function_positive_void() {
        let model = GursonModel::new(0.01, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        // With nonzero void fraction and hydrostatic stress, yield surface shrinks
        let phi_no_void = (500e6_f64 / 500e6_f64).powi(2) - 1.0;
        let phi_void = model.yield_function(500e6, 0.0);
        // phi_void > phi_no_void due to additional void terms
        assert!(phi_void > phi_no_void);
    }

    #[test]
    fn test_gtn_nucleation_rate_peak_at_epsilon_n() {
        let model = GursonModel::new(0.001, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        let rate_at_peak = model.nucleation_rate(0.3);
        let rate_far = model.nucleation_rate(0.8);
        assert!(rate_at_peak > rate_far);
    }

    #[test]
    fn test_gtn_void_growth_positive_triaxiality() {
        let model = GursonModel::new(0.01, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        let rate = model.growth_rate(0.01);
        assert!(rate > 0.0);
    }

    #[test]
    fn test_gtn_void_growth_zero_for_zero_volumetric() {
        let model = GursonModel::new(0.01, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        let rate = model.growth_rate(0.0);
        assert!((rate - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_gtn_failure_detection() {
        let mut model = GursonModel::new(0.24, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        assert!(!model.is_failed());
        model.void_fraction = 0.25;
        assert!(model.is_failed());
    }

    #[test]
    fn test_gtn_void_update_increases() {
        let mut model = GursonModel::new(0.001, 0.1, 0.25, 500e6, 0.3, 0.1, 0.04);
        let f0 = model.void_fraction;
        model.update_void_fraction(0.3, 0.01, 0.001);
        assert!(model.void_fraction >= f0);
    }

    #[test]
    fn test_gtn_triaxiality() {
        let t = GursonModel::triaxiality(100e6, 300e6);
        assert!((t - 1.0 / 3.0).abs() < 1e-12);
    }

    // --- MazarsModel tests ---

    #[test]
    fn test_mazars_below_threshold_no_damage() {
        let model = MazarsModel::new(1e-4, 0.8, 20000.0, 1.4, 5000.0);
        let eps = [1e-5_f64, 0.0, 0.0];
        assert!((model.tension_damage(MazarsModel::equivalent_strain(&eps)) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_mazars_above_threshold_positive_damage() {
        let model = MazarsModel::new(1e-4, 0.8, 20000.0, 1.4, 5000.0);
        let eps_tilde = 5e-4;
        let d = model.tension_damage(eps_tilde);
        assert!(d > 0.0);
    }

    #[test]
    fn test_mazars_equivalent_strain_positive_only() {
        // Compressive strains should not contribute
        let eps = [-0.01_f64, -0.005, 0.001];
        let e_tilde = MazarsModel::equivalent_strain(&eps);
        assert!((e_tilde - 0.001).abs() < 1e-12);
    }

    #[test]
    fn test_mazars_equivalent_strain_all_positive() {
        let eps = [0.003_f64, 0.004, 0.0];
        let e_tilde = MazarsModel::equivalent_strain(&eps);
        assert!((e_tilde - 0.005).abs() < 1e-12);
    }

    #[test]
    fn test_mazars_damage_monotonic_with_loading() {
        let mut model = MazarsModel::new(1e-4, 0.8, 20000.0, 1.4, 5000.0);
        let mut prev_d = 0.0;
        for i in 1..15 {
            let eps = [i as f64 * 5e-5, 0.0, 0.0];
            let d = model.update(&eps, 1.0);
            assert!(d >= prev_d);
            prev_d = d;
        }
    }

    #[test]
    fn test_mazars_damage_bounded() {
        let mut model = MazarsModel::new(1e-4, 0.8, 20000.0, 1.4, 5000.0);
        for i in 1..50 {
            let eps = [i as f64 * 1e-3, 0.0, 0.0];
            let d = model.update(&eps, 1.0);
            assert!((0.0..=1.0).contains(&d));
        }
    }

    // --- CoupledDamagePlasticity tests ---

    #[test]
    fn test_coupled_elastic_return() {
        let mut model =
            CoupledDamagePlasticity::new(200e9, 250e6, 10e6, 5.0, 50e6, 0.5, 1e4, 1.5, 1e3);
        let (stress, dp) = model.return_mapping(100e6);
        assert!((stress - 100e6).abs() < 1.0);
        assert!((dp - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_coupled_plastic_yielding() {
        let mut model =
            CoupledDamagePlasticity::new(200e9, 250e6, 10e6, 5.0, 50e6, 0.5, 1e4, 1.5, 1e3);
        let (stress, dp) = model.return_mapping(300e6);
        assert!(dp > 0.0);
        assert!(stress < 300e6);
    }

    #[test]
    fn test_coupled_damage_accumulates() {
        let mut model =
            CoupledDamagePlasticity::new(200e9, 50e6, 10e6, 5.0, 50e6, 0.5, 1.0, 1.0, 1.0);
        for _ in 0..10 {
            model.return_mapping(200e6);
        }
        assert!(model.damage > 0.0);
    }

    // --- WeibullModel tests ---

    #[test]
    fn test_weibull_pf_zero_at_zero_stress() {
        let model = WeibullModel::new(100e6, 10.0, 1.0);
        assert!((model.failure_probability(0.0) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn test_weibull_pf_at_characteristic_stress() {
        let model = WeibullModel::new(100e6, 10.0, 1.0);
        // P_f(σ_0) = 1 - exp(-1) ≈ 0.6321
        let pf = model.failure_probability(100e6);
        assert!((pf - (1.0 - (-1.0_f64).exp())).abs() < 1e-12);
    }

    #[test]
    fn test_weibull_pf_approaches_one() {
        let model = WeibullModel::new(100e6, 5.0, 1.0);
        let pf = model.failure_probability(1e12);
        assert!(pf > 0.9999);
    }

    #[test]
    fn test_weibull_survival_plus_failure_equals_one() {
        let model = WeibullModel::new(100e6, 8.0, 1.0);
        let stress = 80e6;
        assert!(
            (model.failure_probability(stress) + model.survival_probability(stress) - 1.0).abs()
                < 1e-12
        );
    }

    #[test]
    fn test_weibull_quantile_roundtrip() {
        let model = WeibullModel::new(100e6, 5.0, 1.0);
        let p = 0.5;
        let s = model.quantile(p);
        let pf = model.failure_probability(s);
        assert!((pf - p).abs() < 1e-10);
    }

    #[test]
    fn test_weibull_volume_scaled_larger_volume_higher_pf() {
        let model = WeibullModel::new(100e6, 5.0, 1.0);
        let pf_small = model.volume_scaled_failure_probability(80e6, 0.5);
        let pf_large = model.volume_scaled_failure_probability(80e6, 2.0);
        assert!(pf_large > pf_small);
    }

    // --- FatigueLifeModel tests ---

    #[test]
    fn test_fatigue_strain_amplitude_at_reference() {
        // Δε/2 = σ_f'/E * (2N)^b + ε_f' * (2N)^c
        let model = FatigueLifeModel::new(200e9, 1000e6, -0.1, 0.5, -0.6);
        let amp = model.strain_amplitude(1000.0);
        assert!(amp > 0.0);
    }

    #[test]
    fn test_fatigue_life_roundtrip() {
        let model = FatigueLifeModel::new(200e9, 1000e6, -0.1, 0.5, -0.6);
        let amp = model.strain_amplitude(10000.0);
        let life = model.life_from_strain(amp);
        assert!((life - 10000.0).abs() / 10000.0 < 0.01);
    }

    #[test]
    fn test_fatigue_basquin_stress() {
        let model = FatigueLifeModel::new(200e9, 1000e6, -0.1, 0.5, -0.6);
        let s = model.stress_amplitude_basquin(1000.0);
        assert!(s > 0.0 && s < 1000e6);
    }

    #[test]
    fn test_fatigue_miner_rule_accumulates() {
        let mut model = FatigueLifeModel::new(200e9, 1000e6, -0.1, 0.5, -0.6);
        let amp = model.strain_amplitude(10000.0);
        model.apply_miner_cycle(amp, 1000.0);
        assert!(model.miner_damage > 0.0);
    }

    // --- DamageLocalization tests ---

    #[test]
    fn test_localization_gaussian_weight_max_at_zero() {
        let loc = DamageLocalization::new(0.01, 100.0, 3e6);
        let w0 = loc.gaussian_weight(0.0);
        let w1 = loc.gaussian_weight(0.005);
        assert!(w0 > w1);
    }

    #[test]
    fn test_localization_cone_weight_zero_outside_radius() {
        let loc = DamageLocalization::new(0.01, 100.0, 3e6);
        assert!((loc.cone_weight(0.02) - 0.0).abs() < 1e-15);
    }

    #[test]
    fn test_localization_nonlocal_strain_uniform_field() {
        let loc = DamageLocalization::new(0.05, 100.0, 3e6);
        let positions: Vec<f64> = (0..10).map(|i| i as f64 * 0.01).collect();
        let strains: Vec<f64> = vec![1e-3; 10];
        let eps_nl = loc.nonlocal_strain(0.05, &positions, &strains);
        assert!((eps_nl - 1e-3).abs() < 1e-10);
    }

    #[test]
    fn test_localization_characteristic_length() {
        let loc = DamageLocalization::new(0.01, 100.0, 3e6);
        let lc = loc.characteristic_length(30e9);
        // G_f * E / f_t² = 100 * 30e9 / (3e6)² = 100 * 30e9 / 9e12 = 0.333...
        assert!((lc - 100.0 * 30e9 / (3e6_f64 * 3e6)).abs() < 1e-3);
    }

    #[test]
    fn test_localization_band_width() {
        let loc = DamageLocalization::new(1.0, 100.0, 3e6);
        assert!((loc.band_width() - std::f64::consts::PI).abs() < 1e-12);
    }
}
