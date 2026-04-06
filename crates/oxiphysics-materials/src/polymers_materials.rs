// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! Polymer material models: rubber elasticity, viscoelasticity, degradation.
//!
//! Provides models for:
//! - Freely Jointed Chain (FJC) polymer model
//! - Worm-Like Chain (WLC) for semi-flexible polymers
//! - Neo-Hookean rubber elasticity
//! - Viscoelastic Maxwell and Kelvin-Voigt models
//! - Glass transition WLF equation
//! - Polymer degradation kinetics
//! - Diblock copolymer microphase separation
//! - Entanglement and reptation dynamics

/// Boltzmann constant in J/K
pub const K_BOLTZMANN: f64 = 1.380649e-23;

/// Avogadro's number
pub const N_AVOGADRO: f64 = 6.02214076e23;

/// Freely Jointed Chain polymer model.
///
/// Models a polymer as N rigid segments of length b (Kuhn length),
/// with free rotation about each joint.
#[allow(dead_code)]
pub struct FreelyJointedChain {
    /// Number of Kuhn segments
    pub n_segments: f64,
    /// Kuhn segment length \[m\]
    pub kuhn_length: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}

impl FreelyJointedChain {
    /// Create a new Freely Jointed Chain model.
    ///
    /// # Arguments
    /// * `n_segments` - Number of Kuhn segments N
    /// * `kuhn_length` - Kuhn length b \[m\]
    /// * `temperature` - Temperature \[K\]
    pub fn new(n_segments: f64, kuhn_length: f64, temperature: f64) -> Self {
        Self {
            n_segments,
            kuhn_length,
            temperature,
        }
    }

    /// Mean square end-to-end distance: `R²` = N * b².
    pub fn mean_square_end_to_end(&self) -> f64 {
        self.n_segments * self.kuhn_length * self.kuhn_length
    }

    /// RMS end-to-end distance: `R²`^(1/2) = b * sqrt(N).
    pub fn rms_end_to_end(&self) -> f64 {
        self.mean_square_end_to_end().sqrt()
    }

    /// Contour length: L_c = N * b.
    pub fn contour_length(&self) -> f64 {
        self.n_segments * self.kuhn_length
    }

    /// Inverse Langevin function approximation (Padé approximant).
    ///
    /// L_inv(x) ≈ x * (3 - x²) / (1 - x²) for |x| < 1.
    pub fn inverse_langevin(x: f64) -> f64 {
        let x = x.clamp(-0.9999, 0.9999);
        x * (3.0 - x * x) / (1.0 - x * x)
    }

    /// Force-extension relationship: F = (kT/b) * L_inv(R/L_c).
    ///
    /// # Arguments
    /// * `extension` - End-to-end distance R \[m\], must be < contour length
    pub fn force_extension(&self, extension: f64) -> f64 {
        let l_c = self.contour_length();
        let ratio = (extension / l_c).clamp(-0.9999, 0.9999);
        let kbt = K_BOLTZMANN * self.temperature;
        (kbt / self.kuhn_length) * Self::inverse_langevin(ratio)
    }

    /// Entropic spring constant at small extensions: k = 3kT / (N*b²).
    pub fn spring_constant(&self) -> f64 {
        3.0 * K_BOLTZMANN * self.temperature
            / (self.n_segments * self.kuhn_length * self.kuhn_length)
    }
}

/// Worm-Like Chain polymer model for semi-flexible polymers.
///
/// Uses the Marko-Siggia interpolation formula for force-extension.
#[allow(dead_code)]
pub struct WormLikeChainPolymer {
    /// Persistence length L_p \[m\]
    pub persistence_length: f64,
    /// Contour length L_c \[m\]
    pub contour_length: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}

impl WormLikeChainPolymer {
    /// Create a new Worm-Like Chain model.
    ///
    /// # Arguments
    /// * `persistence_length` - Persistence length L_p \[m\]
    /// * `contour_length` - Contour length L_c \[m\]
    /// * `temperature` - Temperature \[K\]
    pub fn new(persistence_length: f64, contour_length: f64, temperature: f64) -> Self {
        Self {
            persistence_length,
            contour_length,
            temperature,
        }
    }

    /// Force-extension using Marko-Siggia interpolation:
    /// F = (kT/L_p) * \[1/(4*(1-x)²) - 1/4 + x\]
    /// where x = R/L_c.
    ///
    /// # Arguments
    /// * `extension` - End-to-end distance R \[m\], must be < contour length
    pub fn force_extension(&self, extension: f64) -> f64 {
        let x = (extension / self.contour_length).clamp(0.0, 0.9999);
        let kbt = K_BOLTZMANN * self.temperature;
        (kbt / self.persistence_length) * (0.25 / (1.0 - x).powi(2) - 0.25 + x)
    }

    /// Extension at a given force (numerical inversion via bisection).
    ///
    /// # Arguments
    /// * `force` - Applied force \[N\]
    pub fn extension_at_force(&self, force: f64) -> f64 {
        if force <= 0.0 {
            return 0.0;
        }
        let mut lo = 0.0_f64;
        let mut hi = 0.9999 * self.contour_length;
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.force_extension(mid) < force {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// Mean square end-to-end distance in the WLC model:
    /// `R²` = 2 * L_p * L_c * \[1 - (L_p/L_c)*(1 - exp(-L_c/L_p))\]
    pub fn mean_square_end_to_end(&self) -> f64 {
        let ratio = self.contour_length / self.persistence_length;
        2.0 * self.persistence_length
            * self.contour_length
            * (1.0 - (1.0 / ratio) * (1.0 - (-ratio).exp()))
    }
}

/// Neo-Hookean rubber elasticity model.
///
/// Strain energy density: W = (μ/2)*(I_1 - 3) - μ*ln(J) + (λ/2)*ln²(J)
/// where I_1 is the first invariant of the left Cauchy-Green tensor,
/// J is the volumetric deformation, μ is the shear modulus,
/// λ is the Lamé first parameter.
#[allow(dead_code)]
pub struct RubberElasticityNeoHookean {
    /// Shear modulus μ = nkT \[Pa\]
    pub shear_modulus: f64,
    /// Bulk modulus K \[Pa\]
    pub bulk_modulus: f64,
    /// Network chain density n \[chains/m³\]
    pub network_density: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}

impl RubberElasticityNeoHookean {
    /// Create Neo-Hookean rubber model.
    ///
    /// # Arguments
    /// * `network_density` - Number of network chains per unit volume \[1/m³\]
    /// * `bulk_modulus` - Bulk modulus K \[Pa\]
    /// * `temperature` - Temperature \[K\]
    pub fn new(network_density: f64, bulk_modulus: f64, temperature: f64) -> Self {
        let shear_modulus = network_density * K_BOLTZMANN * temperature;
        Self {
            shear_modulus,
            bulk_modulus,
            network_density,
            temperature,
        }
    }

    /// Create with explicit shear modulus.
    ///
    /// # Arguments
    /// * `shear_modulus` - Shear modulus μ \[Pa\]
    /// * `bulk_modulus` - Bulk modulus K \[Pa\]
    pub fn from_moduli(shear_modulus: f64, bulk_modulus: f64) -> Self {
        let temperature = 298.15;
        let network_density = shear_modulus / (K_BOLTZMANN * temperature);
        Self {
            shear_modulus,
            bulk_modulus,
            network_density,
            temperature,
        }
    }

    /// Lamé first parameter: λ = K - (2/3)*μ.
    pub fn lame_lambda(&self) -> f64 {
        self.bulk_modulus - (2.0 / 3.0) * self.shear_modulus
    }

    /// Strain energy density W = (μ/2)*(I_1 - 3) - μ*ln(J) + (λ/2)*ln²(J).
    ///
    /// # Arguments
    /// * `i1` - First invariant of left Cauchy-Green tensor (= λ₁² + λ₂² + λ₃²)
    /// * `j` - Determinant of deformation gradient J = λ₁*λ₂*λ₃
    pub fn strain_energy_density(&self, i1: f64, j: f64) -> f64 {
        let mu = self.shear_modulus;
        let lambda = self.lame_lambda();
        let ln_j = j.ln();
        (mu / 2.0) * (i1 - 3.0) - mu * ln_j + (lambda / 2.0) * ln_j * ln_j
    }

    /// Cauchy stress in uniaxial tension: σ = μ*(λ² - 1/λ) / V
    /// For incompressible rubber (J=1): σ = μ*(λ² - λ^(-1)).
    ///
    /// # Arguments
    /// * `stretch` - Principal stretch ratio λ (= 1 at no deformation)
    pub fn uniaxial_stress(&self, stretch: f64) -> f64 {
        self.shear_modulus * (stretch * stretch - 1.0 / stretch)
    }

    /// Shear stress: τ = μ * γ (for small shear γ).
    pub fn shear_stress(&self, shear_strain: f64) -> f64 {
        self.shear_modulus * shear_strain
    }
}

/// Viscoelastic polymer model (Maxwell + Kelvin-Voigt).
///
/// Maxwell model: spring E in series with dashpot η.
/// Kelvin-Voigt model: spring E in parallel with dashpot η.
#[allow(dead_code)]
pub struct ViscoelasticPolymer {
    /// Elastic modulus E \[Pa\]
    pub elastic_modulus: f64,
    /// Viscosity η \[Pa·s\]
    pub viscosity: f64,
    /// Relaxation time τ = η/E \[s\]
    pub relaxation_time: f64,
}

impl ViscoelasticPolymer {
    /// Create viscoelastic polymer model.
    ///
    /// # Arguments
    /// * `elastic_modulus` - Young's modulus E \[Pa\]
    /// * `viscosity` - Dynamic viscosity η \[Pa·s\]
    pub fn new(elastic_modulus: f64, viscosity: f64) -> Self {
        let relaxation_time = viscosity / elastic_modulus;
        Self {
            elastic_modulus,
            viscosity,
            relaxation_time,
        }
    }

    /// Maxwell stress relaxation: σ(t) = E * ε_0 * exp(-t/τ).
    ///
    /// # Arguments
    /// * `initial_strain` - Applied strain ε_0
    /// * `time` - Time \[s\]
    pub fn maxwell_stress_relaxation(&self, initial_strain: f64, time: f64) -> f64 {
        self.elastic_modulus * initial_strain * (-time / self.relaxation_time).exp()
    }

    /// Maxwell creep compliance: J(t) = (1/E) + t/η.
    pub fn maxwell_creep_compliance(&self, time: f64) -> f64 {
        1.0 / self.elastic_modulus + time / self.viscosity
    }

    /// Kelvin-Voigt creep: ε(t) = (σ/E) * (1 - exp(-t/τ)).
    ///
    /// # Arguments
    /// * `applied_stress` - Applied stress σ \[Pa\]
    /// * `time` - Time \[s\]
    pub fn kelvin_voigt_creep(&self, applied_stress: f64, time: f64) -> f64 {
        (applied_stress / self.elastic_modulus) * (1.0 - (-time / self.relaxation_time).exp())
    }

    /// Kelvin-Voigt stress: σ = E*ε + η*dε/dt.
    ///
    /// # Arguments
    /// * `strain` - Current strain ε
    /// * `strain_rate` - Strain rate dε/dt \[1/s\]
    pub fn kelvin_voigt_stress(&self, strain: f64, strain_rate: f64) -> f64 {
        self.elastic_modulus * strain + self.viscosity * strain_rate
    }

    /// Storage modulus E'(ω) = E * ω²τ² / (1 + ω²τ²).
    pub fn storage_modulus(&self, angular_freq: f64) -> f64 {
        let wt = angular_freq * self.relaxation_time;
        self.elastic_modulus * wt * wt / (1.0 + wt * wt)
    }

    /// Loss modulus E''(ω) = E * ωτ / (1 + ω²τ²).
    pub fn loss_modulus(&self, angular_freq: f64) -> f64 {
        let wt = angular_freq * self.relaxation_time;
        self.elastic_modulus * wt / (1.0 + wt * wt)
    }

    /// Loss tangent tan(δ) = E''/E'.
    pub fn loss_tangent(&self, angular_freq: f64) -> f64 {
        let wt = angular_freq * self.relaxation_time;
        1.0 / wt
    }
}

/// Glass transition temperature model using the WLF equation.
///
/// log₁₀(a_T) = -C₁*(T - T_ref) / (C₂ + (T - T_ref))
#[allow(dead_code)]
pub struct GlassTransition {
    /// Reference temperature T_ref \[K\]
    pub t_ref: f64,
    /// WLF constant C₁ (dimensionless)
    pub c1: f64,
    /// WLF constant C₂ \[K\]
    pub c2: f64,
    /// Glass transition temperature T_g \[K\]
    pub t_g: f64,
}

impl GlassTransition {
    /// Create glass transition model with WLF parameters.
    ///
    /// Universal constants at T_g: C₁ = 17.44, C₂ = 51.6 K.
    ///
    /// # Arguments
    /// * `t_g` - Glass transition temperature \[K\]
    /// * `t_ref` - Reference temperature \[K\] (often = T_g)
    /// * `c1` - WLF C₁ constant (default 17.44)
    /// * `c2` - WLF C₂ constant \[K\] (default 51.6)
    pub fn new(t_g: f64, t_ref: f64, c1: f64, c2: f64) -> Self {
        Self { t_ref, c1, c2, t_g }
    }

    /// WLF shift factor log₁₀(a_T) = -C₁*(T-T_ref)/(C₂+(T-T_ref)).
    pub fn log_shift_factor(&self, temperature: f64) -> f64 {
        let dt = temperature - self.t_ref;
        -self.c1 * dt / (self.c2 + dt)
    }

    /// Shift factor a_T = 10^(log₁₀(a_T)).
    pub fn shift_factor(&self, temperature: f64) -> f64 {
        10.0_f64.powf(self.log_shift_factor(temperature))
    }

    /// Check if temperature is above glass transition.
    pub fn is_rubbery(&self, temperature: f64) -> bool {
        temperature > self.t_g
    }

    /// Relaxation time at temperature T relative to reference.
    pub fn relaxation_time(&self, t_ref_tau: f64, temperature: f64) -> f64 {
        t_ref_tau * self.shift_factor(temperature)
    }
}

/// Polymer degradation model (hydrolysis kinetics).
///
/// First-order: M_n(t) = M_n0 * exp(-k*t)
/// Autocatalytic: dM/dt = -k * M * \[H⁺\]
#[allow(dead_code)]
pub struct PolymerDegradation {
    /// Initial number-average molecular weight M_n0 \[g/mol\]
    pub initial_mw: f64,
    /// First-order degradation rate constant k \[1/s\]
    pub rate_constant: f64,
    /// Autocatalytic coefficient α (dimensionless)
    pub autocatalytic_coeff: f64,
}

impl PolymerDegradation {
    /// Create polymer degradation model.
    ///
    /// # Arguments
    /// * `initial_mw` - Initial molecular weight M_n0 \[g/mol\]
    /// * `rate_constant` - Degradation rate k \[1/s\]
    /// * `autocatalytic_coeff` - Autocatalytic coefficient α
    pub fn new(initial_mw: f64, rate_constant: f64, autocatalytic_coeff: f64) -> Self {
        Self {
            initial_mw,
            rate_constant,
            autocatalytic_coeff,
        }
    }

    /// First-order hydrolysis: M_n(t) = M_n0 * exp(-k*t).
    pub fn molecular_weight(&self, time: f64) -> f64 {
        self.initial_mw * (-self.rate_constant * time).exp()
    }

    /// Degree of degradation: α = 1 - M_n(t)/M_n0.
    pub fn degree_of_degradation(&self, time: f64) -> f64 {
        1.0 - self.molecular_weight(time) / self.initial_mw
    }

    /// Half-life: t_{1/2} = ln(2) / k.
    pub fn half_life(&self) -> f64 {
        std::f64::consts::LN_2 / self.rate_constant
    }

    /// Autocatalytic degradation rate: dM/dt = -k * M * (1 + α * (1 - M/M0)).
    ///
    /// # Arguments
    /// * `current_mw` - Current molecular weight \[g/mol\]
    pub fn autocatalytic_rate(&self, current_mw: f64) -> f64 {
        let degradation = 1.0 - current_mw / self.initial_mw;
        -self.rate_constant * current_mw * (1.0 + self.autocatalytic_coeff * degradation)
    }

    /// Remaining mechanical strength ratio ~ (M_n/M_n0)^(1/2).
    pub fn strength_retention(&self, time: f64) -> f64 {
        (self.molecular_weight(time) / self.initial_mw).sqrt()
    }
}

/// Morphology type in diblock copolymer phase diagram.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub enum CopolymerMorphology {
    /// Disordered (mixed) phase
    Disordered,
    /// Body-centered cubic spheres
    Spheres,
    /// Hexagonally packed cylinders
    Cylinders,
    /// Gyroid bicontinuous network
    Gyroid,
    /// Lamellar (alternating layers)
    Lamellar,
}

/// Diblock copolymer microphase separation model.
///
/// Based on Flory-Huggins theory and self-consistent field theory (SCFT).
#[allow(dead_code)]
pub struct DiblockCopolymer {
    /// Degree of polymerization N
    pub degree_of_polymerization: f64,
    /// Volume fraction of block A, f_A ∈ (0, 1)
    pub volume_fraction_a: f64,
    /// Flory-Huggins interaction parameter χ
    pub chi_parameter: f64,
    /// Monomer reference volume \[m³\]
    pub monomer_volume: f64,
}

impl DiblockCopolymer {
    /// Create diblock copolymer model.
    ///
    /// # Arguments
    /// * `degree_of_polymerization` - Total chain length N
    /// * `volume_fraction_a` - Volume fraction of A block, f_A
    /// * `chi_parameter` - Flory-Huggins χ parameter
    /// * `monomer_volume` - Reference volume per monomer \[m³\]
    pub fn new(
        degree_of_polymerization: f64,
        volume_fraction_a: f64,
        chi_parameter: f64,
        monomer_volume: f64,
    ) -> Self {
        Self {
            degree_of_polymerization,
            volume_fraction_a,
            chi_parameter,
            monomer_volume,
        }
    }

    /// χN segregation parameter.
    pub fn chi_n(&self) -> f64 {
        self.chi_parameter * self.degree_of_polymerization
    }

    /// Mean-field critical χN = 10.495 for symmetric diblock (f_A = 0.5).
    pub fn critical_chi_n() -> f64 {
        10.495
    }

    /// Predict equilibrium morphology based on χN and f_A.
    ///
    /// Phase boundaries from SCFT (approximate):
    /// - Disordered: χN < 10.5
    /// - Spheres: f_A < 0.15 or f_A > 0.85
    /// - Cylinders: 0.15 ≤ f_A < 0.30 or 0.70 < f_A ≤ 0.85
    /// - Gyroid: 0.30 ≤ f_A < 0.35 or 0.65 < f_A ≤ 0.70
    /// - Lamellar: 0.35 ≤ f_A ≤ 0.65
    pub fn morphology(&self) -> CopolymerMorphology {
        if self.chi_n() < 10.5 {
            return CopolymerMorphology::Disordered;
        }
        let f = self.volume_fraction_a;
        if !(0.15..=0.85).contains(&f) {
            CopolymerMorphology::Spheres
        } else if (0.15..0.30).contains(&f) || (f > 0.70 && f <= 0.85) {
            CopolymerMorphology::Cylinders
        } else if (0.30..0.35).contains(&f) || (f > 0.65 && f <= 0.70) {
            CopolymerMorphology::Gyroid
        } else {
            CopolymerMorphology::Lamellar
        }
    }

    /// Flory-Huggins free energy of mixing per monomer:
    /// ΔG_mix / (nkT) = (f_A/N_A)*ln(f_A) + (f_B/N_B)*ln(f_B) + χ*f_A*f_B
    pub fn free_energy_of_mixing(&self) -> f64 {
        let f_a = self.volume_fraction_a;
        let f_b = 1.0 - f_a;
        let n = self.degree_of_polymerization;
        (f_a / n) * f_a.ln() + (f_b / n) * f_b.ln() + self.chi_parameter * f_a * f_b
    }

    /// Lamellar period d ~ b*N^(2/3) * (χN)^(1/6) (strong segregation limit).
    ///
    /// # Arguments
    /// * `segment_length` - Statistical segment length b \[m\]
    pub fn lamellar_period(&self, segment_length: f64) -> f64 {
        let n = self.degree_of_polymerization;
        segment_length * n.powf(2.0 / 3.0) * self.chi_n().powf(1.0 / 6.0)
    }
}

/// Entanglement and reptation dynamics model.
///
/// Describes chain diffusion, plateau modulus, and Rouse-reptation crossover.
#[allow(dead_code)]
pub struct EntanglementModel {
    /// Polymer density ρ \[kg/m³\]
    pub density: f64,
    /// Number-average molecular weight M_n \[kg/mol\]
    pub molecular_weight: f64,
    /// Entanglement molecular weight M_e \[kg/mol\]
    pub entanglement_mw: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Monomeric friction coefficient ζ_0 \[N·s/m\]
    pub friction_coeff: f64,
    /// Segment length b \[m\]
    pub segment_length: f64,
}

impl EntanglementModel {
    /// Create entanglement model.
    ///
    /// # Arguments
    /// * `density` - Polymer density \[kg/m³\]
    /// * `molecular_weight` - Polymer M_n \[kg/mol\]
    /// * `entanglement_mw` - Entanglement M_e \[kg/mol\]
    /// * `temperature` - Temperature \[K\]
    /// * `friction_coeff` - Monomeric friction coefficient \[N·s/m\]
    /// * `segment_length` - Statistical segment length \[m\]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        density: f64,
        molecular_weight: f64,
        entanglement_mw: f64,
        temperature: f64,
        friction_coeff: f64,
        segment_length: f64,
    ) -> Self {
        Self {
            density,
            molecular_weight,
            entanglement_mw,
            temperature,
            friction_coeff,
            segment_length,
        }
    }

    /// Plateau modulus G_N = ρ*R*T / M_e \[Pa\].
    pub fn plateau_modulus(&self) -> f64 {
        let r_gas = 8.314; // J/(mol·K)
        self.density * r_gas * self.temperature / self.entanglement_mw
    }

    /// Number of entanglements per chain Z = M_n / M_e.
    pub fn entanglement_number(&self) -> f64 {
        self.molecular_weight / self.entanglement_mw
    }

    /// Rouse relaxation time: τ_R = ζ_0 * N² * b² / (6π² * kT).
    pub fn rouse_time(&self) -> f64 {
        let n = self.molecular_weight / 0.1e-3; // approximate monomer count
        self.friction_coeff * n * n * self.segment_length * self.segment_length
            / (6.0 * std::f64::consts::PI.powi(2) * K_BOLTZMANN * self.temperature)
    }

    /// Reptation (terminal) time: τ_d = 3 * Z³ * τ_e
    /// (simplified: τ_d ~ Z³ * τ_R).
    pub fn reptation_time(&self) -> f64 {
        let z = self.entanglement_number();
        z * z * z * self.rouse_time() / (z * z)
    }

    /// Self-diffusion coefficient in reptation regime: D ~ kT*M_e / (ζ*N²).
    pub fn diffusion_coefficient(&self) -> f64 {
        let n = self.molecular_weight / 0.1e-3;
        K_BOLTZMANN * self.temperature * self.entanglement_mw
            / (self.friction_coeff * n * n * 0.1e-3)
    }

    /// Is the chain entangled? (M_n > 2 * M_e is the common criterion)
    pub fn is_entangled(&self) -> bool {
        self.molecular_weight > 2.0 * self.entanglement_mw
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. FJC end-to-end scales as sqrt(N)
    #[test]
    fn test_fjc_rms_scales_sqrt_n() {
        let fjc1 = FreelyJointedChain::new(100.0, 1e-9, 300.0);
        let fjc2 = FreelyJointedChain::new(400.0, 1e-9, 300.0);
        let r1 = fjc1.rms_end_to_end();
        let r2 = fjc2.rms_end_to_end();
        // r2/r1 should be sqrt(400/100) = 2
        assert!(
            (r2 / r1 - 2.0).abs() < 1e-10,
            "RMS end-to-end should scale as sqrt(N): ratio={:.6}",
            r2 / r1
        );
    }

    // 2. FJC contour length = N*b
    #[test]
    fn test_fjc_contour_length() {
        let b = 3e-10;
        let n = 200.0;
        let fjc = FreelyJointedChain::new(n, b, 300.0);
        let l_c = fjc.contour_length();
        assert!(
            (l_c - n * b).abs() < EPS,
            "Contour length should be N*b={:.6e}, got {:.6e}",
            n * b,
            l_c
        );
    }

    // 3. FJC mean square end-to-end = N*b²
    #[test]
    fn test_fjc_mean_square_end_to_end() {
        let b = 2e-10;
        let n = 50.0;
        let fjc = FreelyJointedChain::new(n, b, 300.0);
        let r2 = fjc.mean_square_end_to_end();
        let expected = n * b * b;
        assert!(
            (r2 - expected).abs() < EPS,
            "Mean square end-to-end should be {:.6e}, got {:.6e}",
            expected,
            r2
        );
    }

    // 4. FJC spring constant scales as 1/N
    #[test]
    fn test_fjc_spring_constant_scales_inv_n() {
        let fjc1 = FreelyJointedChain::new(100.0, 1e-9, 300.0);
        let fjc2 = FreelyJointedChain::new(200.0, 1e-9, 300.0);
        let k1 = fjc1.spring_constant();
        let k2 = fjc2.spring_constant();
        assert!(
            (k1 / k2 - 2.0).abs() < 1e-10,
            "Spring constant should scale as 1/N"
        );
    }

    // 5. FJC force extension is positive and finite
    #[test]
    fn test_fjc_force_extension_positive() {
        let fjc = FreelyJointedChain::new(100.0, 1e-9, 300.0);
        let f = fjc.force_extension(50e-9); // half contour length
        assert!(
            f > 0.0,
            "Force should be positive at extension, got {:.6e}",
            f
        );
    }

    // 6. WLC extension approaches L as F → ∞
    #[test]
    fn test_wlc_extension_limit() {
        let lp = 50e-9;
        let lc = 1e-6;
        let wlc = WormLikeChainPolymer::new(lp, lc, 300.0);
        // Very large force (100 pN >> kT/L_p): extension should be > 0.98 * L_c
        let ext = wlc.extension_at_force(1e-10);
        assert!(
            ext / lc > 0.98,
            "Extension should approach L_c at large force: {:.6}",
            ext / lc
        );
    }

    // 7. WLC force increases with extension
    #[test]
    fn test_wlc_force_monotone() {
        let wlc = WormLikeChainPolymer::new(50e-9, 1e-6, 300.0);
        let f1 = wlc.force_extension(0.3e-6);
        let f2 = wlc.force_extension(0.7e-6);
        assert!(f2 > f1, "WLC force should increase with extension");
    }

    // 8. WLC mean square end-to-end is positive
    #[test]
    fn test_wlc_mean_square_positive() {
        let wlc = WormLikeChainPolymer::new(50e-9, 1e-6, 300.0);
        let r2 = wlc.mean_square_end_to_end();
        assert!(r2 > 0.0, "Mean square end-to-end should be positive");
    }

    // 9. Neo-Hookean stress proportional to μ
    #[test]
    fn test_neohookean_stress_proportional_mu() {
        let rubber1 = RubberElasticityNeoHookean::from_moduli(1e6, 2e9);
        let rubber2 = RubberElasticityNeoHookean::from_moduli(2e6, 2e9);
        let stretch = 1.5;
        let s1 = rubber1.uniaxial_stress(stretch);
        let s2 = rubber2.uniaxial_stress(stretch);
        assert!(
            (s2 / s1 - 2.0).abs() < 1e-10,
            "Stress should be proportional to μ: ratio={:.6}",
            s2 / s1
        );
    }

    // 10. Neo-Hookean stress is zero at stretch = 1
    #[test]
    fn test_neohookean_zero_stress_at_rest() {
        let rubber = RubberElasticityNeoHookean::from_moduli(1e6, 2e9);
        let stress = rubber.uniaxial_stress(1.0);
        assert!(
            stress.abs() < 1e-9,
            "Stress at rest (λ=1) should be zero, got {:.6}",
            stress
        );
    }

    // 11. Neo-Hookean strain energy is zero at identity
    #[test]
    fn test_neohookean_zero_energy_at_identity() {
        let rubber = RubberElasticityNeoHookean::from_moduli(1e6, 2e9);
        let w = rubber.strain_energy_density(3.0, 1.0); // I1=3 for no deformation
        assert!(
            w.abs() < 1e-9,
            "Strain energy at I1=3,J=1 should be zero, got {:.6}",
            w
        );
    }

    // 12. Maxwell stress relaxes to zero
    #[test]
    fn test_maxwell_relaxes_to_zero() {
        let poly = ViscoelasticPolymer::new(1e6, 1e3);
        let stress_long = poly.maxwell_stress_relaxation(0.01, 100.0 * poly.relaxation_time);
        assert!(
            stress_long < 1e-30,
            "Maxwell stress should relax to zero, got {:.6e}",
            stress_long
        );
    }

    // 13. Maxwell relaxation at t=0 equals E*ε₀
    #[test]
    fn test_maxwell_initial_stress() {
        let poly = ViscoelasticPolymer::new(2e6, 5e3);
        let eps0 = 0.02;
        let s0 = poly.maxwell_stress_relaxation(eps0, 0.0);
        assert!(
            (s0 - poly.elastic_modulus * eps0).abs() < 1e-9,
            "Initial Maxwell stress should be E*ε₀"
        );
    }

    // 14. Kelvin-Voigt creep saturates at σ/E
    #[test]
    fn test_kelvin_voigt_creep_saturation() {
        let poly = ViscoelasticPolymer::new(1e6, 1e4);
        let sigma = 1e4;
        let eps_sat = poly.kelvin_voigt_creep(sigma, 100.0 * poly.relaxation_time);
        let expected = sigma / poly.elastic_modulus;
        assert!(
            (eps_sat - expected).abs() / expected < 1e-4,
            "KV creep should saturate at σ/E={:.6e}, got {:.6e}",
            expected,
            eps_sat
        );
    }

    // 15. Kelvin-Voigt creep at t=0 is zero
    #[test]
    fn test_kelvin_voigt_creep_zero_initial() {
        let poly = ViscoelasticPolymer::new(1e6, 1e4);
        let eps0 = poly.kelvin_voigt_creep(1e4, 0.0);
        assert!(
            eps0.abs() < EPS,
            "KV creep at t=0 should be zero, got {:.6e}",
            eps0
        );
    }

    // 16. WLF shift factor at T_ref = 1
    #[test]
    fn test_wlf_shift_factor_at_tref() {
        let gt = GlassTransition::new(200.0, 200.0, 17.44, 51.6);
        let at = gt.shift_factor(200.0); // T = T_ref
        assert!(
            (at - 1.0).abs() < 1e-10,
            "Shift factor at T_ref should be 1.0, got {:.6}",
            at
        );
    }

    // 17. WLF log shift factor = 0 at T_ref
    #[test]
    fn test_wlf_log_shift_zero_at_tref() {
        let gt = GlassTransition::new(200.0, 200.0, 17.44, 51.6);
        let log_at = gt.log_shift_factor(200.0);
        assert!(
            log_at.abs() < EPS,
            "Log shift factor at T_ref should be 0, got {:.6}",
            log_at
        );
    }

    // 18. WLF: above T_ref, shift factor < 1 (polymer flows faster)
    #[test]
    fn test_wlf_above_tref_shift_less_than_one() {
        let gt = GlassTransition::new(200.0, 200.0, 17.44, 51.6);
        let at = gt.shift_factor(250.0);
        assert!(
            at < 1.0,
            "Above T_ref, shift factor should be < 1 (faster relaxation)"
        );
    }

    // 19. Polymer degradation first-order kinetics
    #[test]
    fn test_degradation_first_order() {
        let deg = PolymerDegradation::new(100e3, 0.001, 0.0);
        let mw_half = deg.molecular_weight(deg.half_life());
        assert!(
            (mw_half - 50e3).abs() / 50e3 < 1e-10,
            "Molecular weight at half-life should be M_n0/2"
        );
    }

    // 20. Degradation at t=0 gives M_n0
    #[test]
    fn test_degradation_initial_mw() {
        let deg = PolymerDegradation::new(50e3, 0.01, 1.0);
        let mw0 = deg.molecular_weight(0.0);
        assert!(
            (mw0 - 50e3).abs() < EPS,
            "Molecular weight at t=0 should be M_n0"
        );
    }

    // 21. Degradation degree of degradation increases with time
    #[test]
    fn test_degradation_degree_increases() {
        let deg = PolymerDegradation::new(100e3, 0.001, 0.5);
        let d1 = deg.degree_of_degradation(100.0);
        let d2 = deg.degree_of_degradation(1000.0);
        assert!(d2 > d1, "Degree of degradation should increase with time");
    }

    // 22. Strength retention decreases with degradation
    #[test]
    fn test_strength_retention_decreases() {
        let deg = PolymerDegradation::new(100e3, 0.001, 0.0);
        let s0 = deg.strength_retention(0.0);
        let s1 = deg.strength_retention(1000.0);
        assert!(
            (s0 - 1.0).abs() < EPS,
            "Strength retention at t=0 should be 1"
        );
        assert!(s1 < s0, "Strength retention should decrease over time");
    }

    // 23. Diblock: χN < 10.5 → disordered
    #[test]
    fn test_diblock_disordered_below_critical() {
        let copol = DiblockCopolymer::new(50.0, 0.5, 0.1, 1e-28);
        // χN = 0.1 * 50 = 5 < 10.5
        assert_eq!(copol.morphology(), CopolymerMorphology::Disordered);
    }

    // 24. Diblock: χN = 10.5, f=0.5 → lamellar
    #[test]
    fn test_diblock_lamellar_symmetric() {
        let copol = DiblockCopolymer::new(100.0, 0.5, 0.105, 1e-28);
        // χN = 10.5
        assert_eq!(copol.morphology(), CopolymerMorphology::Lamellar);
    }

    // 25. Diblock: χN = 20, f=0.1 → spheres
    #[test]
    fn test_diblock_spheres_asymmetric() {
        let copol = DiblockCopolymer::new(100.0, 0.1, 0.2, 1e-28);
        assert_eq!(copol.morphology(), CopolymerMorphology::Spheres);
    }

    // 26. Diblock: χN = 20, f=0.22 → cylinders
    #[test]
    fn test_diblock_cylinders() {
        let copol = DiblockCopolymer::new(100.0, 0.22, 0.2, 1e-28);
        assert_eq!(copol.morphology(), CopolymerMorphology::Cylinders);
    }

    // 27. Diblock: χN = 20, f=0.32 → gyroid
    #[test]
    fn test_diblock_gyroid() {
        let copol = DiblockCopolymer::new(100.0, 0.32, 0.2, 1e-28);
        assert_eq!(copol.morphology(), CopolymerMorphology::Gyroid);
    }

    // 28. PDMS plateau modulus ~ 0.2 MPa
    #[test]
    fn test_pdms_plateau_modulus() {
        // PDMS: density=970 kg/m³, M_e=12000 g/mol = 12 kg/mol, T=298K
        // G_N = ρ*R*T/M_e = 970 * 8.314 * 298 / 12 ≈ 200 kPa
        let model = EntanglementModel::new(970.0, 1.0, 12.0, 298.0, 1e-11, 4.3e-10);
        let g_n = model.plateau_modulus();
        let expected = 970.0 * 8.314 * 298.0 / 12.0;
        assert!(
            (g_n - expected).abs() / expected < 0.01,
            "PDMS G_N should be ~{:.0} Pa, got {:.6e}",
            expected,
            g_n
        );
    }

    // 29. Entanglement: Z = M_n / M_e
    #[test]
    fn test_entanglement_number() {
        let model = EntanglementModel::new(900.0, 0.1, 0.01, 300.0, 1e-11, 5e-10);
        let z = model.entanglement_number();
        assert!(
            (z - 10.0).abs() < EPS,
            "Z should be M/M_e = 10, got {:.6}",
            z
        );
    }

    // 30. Entangled check: M > 2*M_e
    #[test]
    fn test_entanglement_check() {
        let entangled = EntanglementModel::new(900.0, 0.1, 0.01, 300.0, 1e-11, 5e-10);
        let unentangled = EntanglementModel::new(900.0, 0.015, 0.01, 300.0, 1e-11, 5e-10);
        assert!(
            entangled.is_entangled(),
            "M=0.1 > 2*M_e=0.02 should be entangled"
        );
        assert!(
            !unentangled.is_entangled(),
            "M=0.015 < 2*M_e=0.02 should be unentangled"
        );
    }

    // 31. Diblock free energy of mixing is finite
    #[test]
    fn test_diblock_free_energy_finite() {
        let copol = DiblockCopolymer::new(100.0, 0.5, 0.2, 1e-28);
        let dg = copol.free_energy_of_mixing();
        assert!(dg.is_finite(), "Free energy of mixing should be finite");
    }

    // 32. Viscoelastic loss tangent at ω*τ = 1 equals 1
    #[test]
    fn test_loss_tangent_at_wt1() {
        let poly = ViscoelasticPolymer::new(1e6, 1e4);
        let tau = poly.relaxation_time;
        let tan_d = poly.loss_tangent(1.0 / tau);
        assert!(
            (tan_d - 1.0).abs() < 1e-10,
            "Loss tangent at ω*τ=1 should be 1, got {:.6}",
            tan_d
        );
    }

    // 33. Maxwell creep compliance increases with time
    #[test]
    fn test_maxwell_creep_compliance_increases() {
        let poly = ViscoelasticPolymer::new(1e6, 1e4);
        let j1 = poly.maxwell_creep_compliance(1.0);
        let j2 = poly.maxwell_creep_compliance(2.0);
        assert!(
            j2 > j1,
            "Maxwell creep compliance should increase with time"
        );
    }

    // 34. WLC Marko-Siggia at low force gives approximately linear response
    #[test]
    fn test_wlc_linear_at_low_force() {
        let lp = 50e-9;
        let lc = 1e-6;
        let wlc = WormLikeChainPolymer::new(lp, lc, 300.0);
        let ext1 = wlc.extension_at_force(1e-15);
        let ext2 = wlc.extension_at_force(2e-15);
        // At very low force, extension is small and roughly proportional to force
        assert!(ext2 > ext1, "Extension should increase with force");
    }

    // 35. Plateau modulus scales with temperature
    #[test]
    fn test_plateau_modulus_scales_temperature() {
        let model1 = EntanglementModel::new(970.0, 1.0, 0.012, 300.0, 1e-11, 4.3e-10);
        let model2 = EntanglementModel::new(970.0, 1.0, 0.012, 600.0, 1e-11, 4.3e-10);
        let g1 = model1.plateau_modulus();
        let g2 = model2.plateau_modulus();
        assert!(
            (g2 / g1 - 2.0).abs() < 1e-10,
            "Plateau modulus should double when T doubles: ratio={:.6}",
            g2 / g1
        );
    }
}
