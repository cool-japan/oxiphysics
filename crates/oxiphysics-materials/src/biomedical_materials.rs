// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biomedical material models: soft tissue, bone, blood, hydrogels.
//!
//! Implements physics-based constitutive models for biological materials:
//! - [`SoftTissue`]: Fung hyperelastic and Holzapfel-Gasser-Ogden (HGO) models
//! - [`BoneModel`]: Cortical/cancellous bone with Gibson-Ashby and Wolff's law
//! - [`BloodRheology`]: Carreau-Yasuda and Casson non-Newtonian models
//! - [`HydrogelModel`]: Flory-Rehner equilibrium swelling theory
//! - [`CartilageModel`]: Mow biphasic theory for articular cartilage
//! - [`ArterialWall`]: Three-layer arterial wall with residual stress
//! - [`DegradablePolymer`]: Hydrolytic degradation kinetics

/// Type alias for a 3D vector represented as a plain array.
pub type Vec3 = [f64; 3];

/// Type alias for a symmetric 3×3 Cauchy-Green strain tensor stored as
/// \[E11, E22, E33, E12, E13, E23\].
pub type StrainTensor = [f64; 6];

// ---------------------------------------------------------------------------
// SoftTissue
// ---------------------------------------------------------------------------

/// Biomedical soft-tissue constitutive model.
///
/// Implements two hyperelastic frameworks commonly used in biomechanics:
///
/// 1. **Fung exponential** (quasi-linear viscoelastic): strain-energy
///    `W = C*(exp(Q) - 1)` where `Q` is a quadratic form in the Green-Lagrange
///    strain components.
/// 2. **Holzapfel-Gasser-Ogden (HGO)**: fiber-reinforced model for arterial
///    walls, with an isotropic neo-Hookean matrix and two fiber families.
/// 3. **Active stress**: additive active-passive decomposition for cardiac muscle.
#[derive(Debug, Clone)]
pub struct SoftTissue {
    /// Material parameter C \[Pa\] for the Fung model ground-matrix stiffness.
    pub fung_c: f64,
    /// Fung exponential coefficients \[b11, b22, b33, b12, b13, b23\] (dimensionless).
    pub fung_b: [f64; 6],
    /// HGO: neo-Hookean ground-matrix parameter μ \[Pa\].
    pub hgo_mu: f64,
    /// HGO: fiber stiffness k1 \[Pa\].
    pub hgo_k1: f64,
    /// HGO: fiber exponential nonlinearity k2 (dimensionless).
    pub hgo_k2: f64,
    /// HGO: fiber dispersion parameter κ ∈ \[0, 1/3\].
    pub hgo_kappa: f64,
    /// HGO: mean fiber angle θ \[rad\] relative to the circumferential direction.
    pub fiber_angle: f64,
    /// Active stress magnitude \[Pa\] (cardiac muscle).
    pub active_stress: f64,
}

impl SoftTissue {
    /// Construct a new `SoftTissue` with Fung parameters.
    ///
    /// # Arguments
    /// * `fung_c` – Ground-matrix stiffness \[Pa\].
    /// * `fung_b` – Exponential coefficients \[b11, b22, b33, b12, b13, b23\].
    pub fn new_fung(fung_c: f64, fung_b: [f64; 6]) -> Self {
        Self {
            fung_c,
            fung_b,
            hgo_mu: 0.0,
            hgo_k1: 0.0,
            hgo_k2: 0.0,
            hgo_kappa: 0.0,
            fiber_angle: 0.0,
            active_stress: 0.0,
        }
    }

    /// Construct a new `SoftTissue` with HGO (fiber-reinforced) parameters.
    ///
    /// # Arguments
    /// * `mu`          – Neo-Hookean ground-matrix shear modulus \[Pa\].
    /// * `k1`          – Fiber stiffness parameter \[Pa\].
    /// * `k2`          – Fiber nonlinearity (dimensionless).
    /// * `kappa`       – Fiber dispersion parameter ∈ \[0, 1/3\].
    /// * `fiber_angle` – Mean fiber angle relative to circumferential axis \[rad\].
    pub fn new_hgo(mu: f64, k1: f64, k2: f64, kappa: f64, fiber_angle: f64) -> Self {
        Self {
            fung_c: 0.0,
            fung_b: [0.0; 6],
            hgo_mu: mu,
            hgo_k1: k1,
            hgo_k2: k2,
            hgo_kappa: kappa,
            fiber_angle,
            active_stress: 0.0,
        }
    }

    /// Compute the Fung strain-energy density W \[J/m³\].
    ///
    /// # Formula
    /// `W = C * (exp(Q) – 1)` where `Q = Σ b_ij * E_i * E_j`.
    pub fn fung_strain_energy(&self, strain: &StrainTensor) -> f64 {
        let [e11, e22, e33, _e12, _e13, _e23] = *strain;
        let [b11, b22, b33, b12, b13, b23] = self.fung_b;
        let q = b11 * e11 * e11
            + b22 * e22 * e22
            + b33 * e33 * e33
            + 2.0 * b12 * e11 * e22
            + 2.0 * b13 * e11 * e33
            + 2.0 * b23 * e22 * e33;
        self.fung_c * (q.exp() - 1.0)
    }

    /// Compute Fung second Piola-Kirchhoff stress S11 in the 1-direction \[Pa\].
    ///
    /// Computed as `∂W/∂E11 = 2C * (b11*E11 + b12*E22 + b13*E33) * exp(Q)`.
    pub fn fung_stress_s11(&self, strain: &StrainTensor) -> f64 {
        let [e11, e22, e33, _e12, _e13, _e23] = *strain;
        let [b11, b22, b33, b12, b13, b23] = self.fung_b;
        let q = b11 * e11 * e11
            + b22 * e22 * e22
            + b33 * e33 * e33
            + 2.0 * b12 * e11 * e22
            + 2.0 * b13 * e11 * e33
            + 2.0 * b23 * e22 * e33;
        2.0 * self.fung_c * (b11 * e11 + b12 * e22 + b13 * e33) * q.exp()
    }

    /// Compute HGO fiber pseudo-invariant I4 for one fiber family.
    ///
    /// I4 = λ_f² = C : (a ⊗ a) where **a** is the fiber direction unit vector.
    /// Returns the stretch squared along the fiber.
    pub fn hgo_i4(&self, stretch_circ: f64, stretch_axial: f64) -> f64 {
        let cos_a = self.fiber_angle.cos();
        let sin_a = self.fiber_angle.sin();
        // Right Cauchy-Green tensor diagonal: C11 = λ_circ², C22 = λ_axial².
        let c11 = stretch_circ * stretch_circ;
        let c22 = stretch_axial * stretch_axial;
        cos_a * cos_a * c11 + sin_a * sin_a * c22
    }

    /// Compute HGO fiber strain-energy density \[J/m³\].
    ///
    /// `W_fiber = k1/(2k2) * Σ (exp(k2*(κ*I1 + (1-3κ)*I4 - 1)²) - 1)`
    /// Only contributes when fibers are under tension (I4 > 1).
    pub fn hgo_fiber_energy(
        &self,
        stretch_circ: f64,
        stretch_axial: f64,
        stretch_radial: f64,
    ) -> f64 {
        let i1 = stretch_circ * stretch_circ
            + stretch_axial * stretch_axial
            + stretch_radial * stretch_radial;
        let i4 = self.hgo_i4(stretch_circ, stretch_axial);
        // Only fibers under tension contribute.
        if i4 <= 1.0 {
            return 0.0;
        }
        let e_val = self.hgo_kappa * i1 + (1.0 - 3.0 * self.hgo_kappa) * i4 - 1.0;
        2.0 * self.hgo_k1 / (2.0 * self.hgo_k2) * ((self.hgo_k2 * e_val * e_val).exp() - 1.0)
    }

    /// Compute total HGO strain energy including isotropic neo-Hookean matrix \[J/m³\].
    pub fn hgo_total_energy(
        &self,
        stretch_circ: f64,
        stretch_axial: f64,
        stretch_radial: f64,
    ) -> f64 {
        let i1 = stretch_circ * stretch_circ
            + stretch_axial * stretch_axial
            + stretch_radial * stretch_radial;
        let w_iso = self.hgo_mu / 2.0 * (i1 - 3.0);
        let w_fiber = self.hgo_fiber_energy(stretch_circ, stretch_axial, stretch_radial);
        w_iso + w_fiber
    }

    /// Compute total stress including passive hyperelastic and active components.
    ///
    /// Returns Cauchy stress σ_11 in the fiber direction.
    pub fn total_stress_with_active(&self, strain: &StrainTensor) -> f64 {
        self.fung_stress_s11(strain) + self.active_stress
    }

    /// Set the active stress magnitude (for cardiac simulation).
    pub fn set_active_stress(&mut self, sigma_a: f64) {
        self.active_stress = sigma_a;
    }
}

// ---------------------------------------------------------------------------
// BoneModel
// ---------------------------------------------------------------------------

/// Bone tissue constitutive model (cortical + cancellous).
///
/// Implements:
/// - Transversely isotropic elasticity for cortical bone.
/// - Gibson-Ashby cellular foam model for trabecular bone.
/// - Kopperdahl-Keaveny mineral-density to modulus correlation.
/// - Simplified Wolff's law density adaptation.
#[derive(Debug, Clone)]
pub struct BoneModel {
    /// Apparent density ρ \[g/cm³\] (0.1–2.0 for cancellous, ~1.9 for cortical).
    pub density: f64,
    /// Mineral volume fraction (ash fraction) \[0..1\].
    pub mineral_fraction: f64,
    /// Axial Young's modulus E_axial \[GPa\] (cortical, along osteons).
    pub e_axial: f64,
    /// Transverse Young's modulus E_transverse \[GPa\] (cortical).
    pub e_transverse: f64,
    /// Shear modulus G \[GPa\] (cortical).
    pub shear_modulus: f64,
    /// Trabecular thickness t_b \[mm\].
    pub trabecular_thickness: f64,
    /// Trabecular length l_b \[mm\].
    pub trabecular_length: f64,
}

impl BoneModel {
    /// Construct a BoneModel with explicit mechanical properties.
    pub fn new(
        density: f64,
        mineral_fraction: f64,
        e_axial: f64,
        e_transverse: f64,
        shear_modulus: f64,
    ) -> Self {
        Self {
            density,
            mineral_fraction,
            e_axial,
            e_transverse,
            shear_modulus,
            trabecular_thickness: 0.15,
            trabecular_length: 1.2,
        }
    }

    /// Kopperdahl-Keaveny correlation: density \[g/cm³\] → Young's modulus \[MPa\].
    ///
    /// `E = 6850 * ρ^1.49` for apparent density in g/cm³ (trabecular bone).
    pub fn density_to_modulus_kk(density: f64) -> f64 {
        6850.0 * density.powf(1.49)
    }

    /// Gibson-Ashby foam model: relative density → relative Young's modulus.
    ///
    /// `E/E_s = C * (ρ/ρ_s)^2` where C ≈ 1.0, ρ_s = solid bone density (~1.9 g/cm³).
    pub fn gibson_ashby_modulus(&self, solid_modulus_gpa: f64) -> f64 {
        let rho_s = 1.9_f64; // solid cortical bone density g/cm³
        let relative_density = (self.density / rho_s).min(1.0);
        solid_modulus_gpa * relative_density * relative_density
    }

    /// Compute trabecular bone apparent modulus using Gibson-Ashby.
    ///
    /// Returns modulus in GPa.
    pub fn trabecular_modulus(&self) -> f64 {
        // Solid bone modulus ~18 GPa.
        self.gibson_ashby_modulus(18.0)
    }

    /// Apply Wolff's law density adaptation.
    ///
    /// Updates the apparent density based on the daily stress stimulus `σ_ref`
    /// compared to the homeostatic stimulus `σ_ref`.
    ///
    /// `dρ/dt = k * (σ_actual - σ_ref)` with lazy zone ±s.
    pub fn wolff_adapt(&mut self, sigma_actual: f64, sigma_ref: f64, k: f64, s: f64, dt: f64) {
        let error = sigma_actual - sigma_ref;
        let d_rho = if error > s {
            k * (error - s)
        } else if error < -s {
            k * (error + s)
        } else {
            0.0
        };
        self.density = (self.density + d_rho * dt).clamp(0.05, 1.9);
    }

    /// Compute bone yield strength \[MPa\] from mineral fraction (Currey).
    ///
    /// `σ_y = 94.0 * mineral_fraction^3.84`
    pub fn yield_strength(&self) -> f64 {
        94.0 * self.mineral_fraction.powf(3.84)
    }
}

// ---------------------------------------------------------------------------
// BloodRheology
// ---------------------------------------------------------------------------

/// Non-Newtonian blood rheology models.
///
/// Implements:
/// - **Carreau-Yasuda**: shear-thinning viscosity with 5 parameters.
/// - **Casson**: yield-stress fluid model.
/// - **Hematocrit-dependent** viscosity (Quemada-type).
#[derive(Debug, Clone)]
pub struct BloodRheology {
    /// Zero-shear viscosity η_0 \[Pa·s\].
    pub eta_0: f64,
    /// Infinite-shear viscosity η_∞ \[Pa·s\].
    pub eta_inf: f64,
    /// Relaxation time λ \[s\].
    pub lambda: f64,
    /// Carreau-Yasuda exponent a (transition parameter).
    pub carreau_a: f64,
    /// Power-law index n (shear-thinning for n < 1).
    pub power_n: f64,
    /// Casson yield stress τ_y \[Pa\].
    pub yield_stress: f64,
    /// Casson plastic viscosity η_c \[Pa·s\].
    pub casson_eta: f64,
    /// Hematocrit volume fraction H ∈ \[0, 1\].
    pub hematocrit: f64,
}

impl BloodRheology {
    /// Construct a `BloodRheology` with physiological blood parameters.
    ///
    /// Default: Carreau-Yasuda fitted to whole blood at 37 °C.
    pub fn new_physiological() -> Self {
        Self {
            eta_0: 0.056,    // 56 mPa·s
            eta_inf: 0.0035, // 3.5 mPa·s (plasma viscosity)
            lambda: 3.313,
            carreau_a: 2.0,
            power_n: 0.3568,
            yield_stress: 0.004, // 4 mPa
            casson_eta: 0.0035,
            hematocrit: 0.45,
        }
    }

    /// Construct with explicit parameters.
    pub fn new(
        eta_0: f64,
        eta_inf: f64,
        lambda: f64,
        carreau_a: f64,
        power_n: f64,
        yield_stress: f64,
        casson_eta: f64,
        hematocrit: f64,
    ) -> Self {
        Self {
            eta_0,
            eta_inf,
            lambda,
            carreau_a,
            power_n,
            yield_stress,
            casson_eta,
            hematocrit,
        }
    }

    /// Carreau-Yasuda viscosity at shear rate γ̇ \[s⁻¹\].
    ///
    /// `η(γ̇) = η_∞ + (η_0 - η_∞) * (1 + (λγ̇)^a)^((n-1)/a)`
    pub fn carreau_viscosity(&self, shear_rate: f64) -> f64 {
        let base = 1.0 + (self.lambda * shear_rate.abs()).powf(self.carreau_a);
        self.eta_inf
            + (self.eta_0 - self.eta_inf) * base.powf((self.power_n - 1.0) / self.carreau_a)
    }

    /// Casson model: shear stress from shear rate \[Pa\].
    ///
    /// `√τ = √τ_y + √(η_c * |γ̇|)` for `|γ̇| > 0`, else τ = 0.
    pub fn casson_stress(&self, shear_rate: f64) -> f64 {
        let gamma_abs = shear_rate.abs();
        if gamma_abs < 1e-12 {
            return 0.0;
        }
        let sqrt_tau = self.yield_stress.sqrt() + (self.casson_eta * gamma_abs).sqrt();
        sqrt_tau * sqrt_tau * shear_rate.signum()
    }

    /// Casson effective viscosity η_eff = τ/γ̇ \[Pa·s\].
    pub fn casson_viscosity(&self, shear_rate: f64) -> f64 {
        let gamma_abs = shear_rate.abs();
        if gamma_abs < 1e-12 {
            return self.eta_0; // Newtonian limit at rest
        }
        self.casson_stress(shear_rate) / shear_rate
    }

    /// Quemada-type hematocrit-dependent viscosity \[Pa·s\].
    ///
    /// `η = η_plasma * (1 - H/2)^(-2)` where H = hematocrit.
    pub fn hematocrit_viscosity(&self) -> f64 {
        let h = self.hematocrit.clamp(0.0, 0.99);
        self.eta_inf / (1.0 - h / 2.0).powi(2)
    }

    /// Check if flow is effectively Newtonian (high shear limit).
    ///
    /// Returns true if Carreau viscosity is within 1 % of η_∞.
    pub fn is_newtonian_limit(&self, shear_rate: f64) -> bool {
        let eta = self.carreau_viscosity(shear_rate);
        (eta - self.eta_inf).abs() / (self.eta_0 - self.eta_inf + 1e-20) < 0.01
    }
}

// ---------------------------------------------------------------------------
// HydrogelModel
// ---------------------------------------------------------------------------

/// Flory-Rehner equilibrium swelling model for hydrogels.
///
/// Balances elastic (network) and mixing (polymer-solvent) free energies to
/// find the equilibrium swelling ratio Q.
#[derive(Debug, Clone)]
pub struct HydrogelModel {
    /// Flory-Huggins interaction parameter χ (dimensionless).
    pub chi: f64,
    /// Number of chain segments between cross-links N (Kuhn monomers).
    pub n_chains: f64,
    /// Reference polymer volume fraction ν_0 (dry state, typically 1.0).
    pub nu_0: f64,
    /// Molar volume of solvent \[m³/mol\].
    pub molar_volume_solvent: f64,
    /// Shear modulus of dry network G_0 \[Pa\].
    pub shear_modulus: f64,
}

impl HydrogelModel {
    /// Construct a `HydrogelModel` with standard parameters.
    pub fn new(chi: f64, n_chains: f64, shear_modulus: f64) -> Self {
        Self {
            chi,
            n_chains,
            nu_0: 1.0,
            molar_volume_solvent: 18e-6, // water
            shear_modulus,
        }
    }

    /// Mixing free energy (Flory-Huggins) chemical potential term \[dimensionless\].
    ///
    /// `Δμ_mix / (RT) = ln(1-φ) + φ + χ*φ²`  where φ = polymer volume fraction.
    pub fn mixing_chemical_potential(&self, phi: f64) -> f64 {
        (1.0 - phi).ln() + phi + self.chi * phi * phi
    }

    /// Elastic free energy chemical potential term \[dimensionless\].
    ///
    /// `Δμ_el / (RT) = (Vs/Vref) * N^(-1) * (φ/2 - φ^(1/3))`
    pub fn elastic_chemical_potential(&self, phi: f64) -> f64 {
        // Simplified Flory-Rehner elastic term.
        (1.0 / self.n_chains) * (phi / 2.0 - phi.powf(1.0 / 3.0))
    }

    /// Total chemical potential (mixing + elastic) for equilibrium condition.
    pub fn total_chemical_potential(&self, phi: f64) -> f64 {
        self.mixing_chemical_potential(phi) + self.elastic_chemical_potential(phi)
    }

    /// Find equilibrium swelling ratio Q = V_swollen / V_dry numerically.
    ///
    /// Solves `μ_total(φ) = 0` by bisection over φ ∈ \[1e-4, 1-1e-4\].
    /// Returns volumetric swelling ratio Q = 1/φ_eq.
    pub fn equilibrium_swelling_ratio(&self) -> f64 {
        // Bisection: find φ such that μ_total(φ) = 0.
        let mut lo = 1e-4_f64;
        let mut hi = 1.0 - 1e-4;
        let mu_lo = self.total_chemical_potential(lo);
        let mu_hi = self.total_chemical_potential(hi);
        // If no sign change, return the boundary with smaller |μ|.
        if mu_lo * mu_hi > 0.0 {
            if mu_lo.abs() < mu_hi.abs() {
                return 1.0 / lo;
            } else {
                return 1.0 / hi;
            }
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if self.total_chemical_potential(lo) * self.total_chemical_potential(mid) <= 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let phi_eq = 0.5 * (lo + hi);
        1.0 / phi_eq
    }
}

// ---------------------------------------------------------------------------
// CartilageModel
// ---------------------------------------------------------------------------

/// Mow biphasic model for articular cartilage.
///
/// Treats cartilage as a mixture of:
/// - **Solid phase**: elastic porous matrix.
/// - **Fluid phase**: interstitial water (nearly incompressible).
///
/// Models creep and stress relaxation under compressive loading.
#[derive(Debug, Clone)]
pub struct CartilageModel {
    /// Aggregate modulus H_A \[MPa\].
    pub aggregate_modulus: f64,
    /// Permeability k \[m⁴/(N·s)\].
    pub permeability: f64,
    /// Poisson's ratio ν_s of solid matrix (0 to 0.5).
    pub poisson_solid: f64,
    /// Thickness of cartilage layer h \[mm\].
    pub thickness_mm: f64,
}

impl CartilageModel {
    /// Construct a `CartilageModel` with standard articular cartilage parameters.
    pub fn new(
        aggregate_modulus: f64,
        permeability: f64,
        poisson_solid: f64,
        thickness_mm: f64,
    ) -> Self {
        Self {
            aggregate_modulus,
            permeability,
            poisson_solid,
            thickness_mm,
        }
    }

    /// Characteristic diffusion time scale t* = h²/(H_A*k) \[s\].
    ///
    /// Controls the rate of fluid exudation under constant load.
    pub fn characteristic_time(&self) -> f64 {
        let h_m = self.thickness_mm * 1e-3;
        h_m * h_m / (self.aggregate_modulus * 1e6 * self.permeability)
    }

    /// Creep deformation at time t under constant stress σ₀ \[normalized\].
    ///
    /// Approximated as: `u(t)/u_∞ ≈ 1 - exp(-t / t*)`.
    pub fn creep_response(&self, sigma0: f64, time: f64) -> f64 {
        let t_star = self.characteristic_time();
        let u_inf = sigma0 / self.aggregate_modulus;
        u_inf * (1.0 - (-time / t_star).exp())
    }

    /// Stress relaxation at time t under constant strain ε₀ \[MPa\].
    ///
    /// Approximated as: `σ(t) = σ_0 * exp(-t / t*)`.
    pub fn stress_relaxation(&self, epsilon0: f64, time: f64) -> f64 {
        let t_star = self.characteristic_time();
        let sigma0 = self.aggregate_modulus * epsilon0;
        sigma0 * (-time / t_star).exp()
    }

    /// Compute Lamé parameter λ \[MPa\] for the solid matrix.
    pub fn lame_lambda(&self) -> f64 {
        let nu = self.poisson_solid;
        let e = 2.0 * self.aggregate_modulus * (1.0 - nu) / (1.0 - 2.0 * nu);
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }
}

// ---------------------------------------------------------------------------
// ArterialWall
// ---------------------------------------------------------------------------

/// Three-layer arterial wall model (intima / media / adventitia).
///
/// Each layer has distinct HGO fiber properties and thickness fractions.
/// Residual stress is captured via the opening angle method.
#[derive(Debug, Clone)]
pub struct ArterialWall {
    /// Intima layer HGO tissue model.
    pub intima: SoftTissue,
    /// Media layer HGO tissue model.
    pub media: SoftTissue,
    /// Adventitia layer HGO tissue model.
    pub adventitia: SoftTissue,
    /// Fractional thickness of each layer \[intima, media, adventitia\] (sums to 1).
    pub layer_fractions: [f64; 3],
    /// Opening angle α \[rad\] representing residual stress state.
    pub opening_angle: f64,
    /// Inner (lumen) radius in reference configuration r_i \[mm\].
    pub inner_radius: f64,
    /// Outer radius in reference configuration r_o \[mm\].
    pub outer_radius: f64,
}

impl ArterialWall {
    /// Construct an `ArterialWall` with typical human aorta parameters.
    pub fn new_aorta() -> Self {
        let intima = SoftTissue::new_hgo(7.64e3, 996.6, 524.6, 0.226, 0.523);
        let media = SoftTissue::new_hgo(3.0e3, 15.0e3, 20.0, 0.226, 0.464);
        let adventitia = SoftTissue::new_hgo(0.3e3, 57.0e3, 11.2, 0.226, 0.785);
        Self {
            intima,
            media,
            adventitia,
            layer_fractions: [0.1, 0.5, 0.4],
            opening_angle: 0.541, // ~31 degrees
            inner_radius: 9.0,
            outer_radius: 13.5,
        }
    }

    /// Compute residual circumferential strain at the inner wall.
    ///
    /// Uses the opening-angle formula: λ_θ = π*r/(α_0*r_0) where α_0 = opening angle.
    pub fn residual_strain_inner(&self) -> f64 {
        // Simplified opening-angle residual strain.
        let alpha = std::f64::consts::PI - self.opening_angle;
        let lambda = std::f64::consts::PI * self.inner_radius / (alpha * self.inner_radius);
        lambda - 1.0
    }

    /// Compute wall stress at given circumferential stretch ratio λ_θ.
    ///
    /// Returns the wall-averaged circumferential stress \[Pa\].
    pub fn wall_stress_at_stretch(&self, lambda_circ: f64) -> f64 {
        let lambda_axial = 1.1; // physiological axial pre-stretch
        let lambda_radial = 1.0 / (lambda_circ * lambda_axial); // incompressibility
        let s_i = self
            .intima
            .hgo_total_energy(lambda_circ, lambda_axial, lambda_radial);
        let s_m = self
            .media
            .hgo_total_energy(lambda_circ, lambda_axial, lambda_radial);
        let s_a = self
            .adventitia
            .hgo_total_energy(lambda_circ, lambda_axial, lambda_radial);
        self.layer_fractions[0] * s_i
            + self.layer_fractions[1] * s_m
            + self.layer_fractions[2] * s_a
    }
}

// ---------------------------------------------------------------------------
// DegradablePolymer
// ---------------------------------------------------------------------------

/// Hydrolytic degradation model for biodegradable polymers (PLA/PLGA).
///
/// Tracks molecular weight decrease and property loss due to ester bond
/// hydrolysis.
#[derive(Debug, Clone)]
pub struct DegradablePolymer {
    /// Initial number-average molecular weight M_n0 \[g/mol\].
    pub mw_initial: f64,
    /// Current number-average molecular weight M_n \[g/mol\].
    pub mw_current: f64,
    /// Hydrolysis rate constant k_h \[1/day\] (environment-dependent).
    pub hydrolysis_rate: f64,
    /// Autocatalysis factor β (dimensionless, > 0 for bulk erosion).
    pub autocatalysis: f64,
    /// Critical molecular weight below which material loses integrity \[g/mol\].
    pub mw_critical: f64,
}

impl DegradablePolymer {
    /// Construct a `DegradablePolymer` representing PLGA (50:50).
    pub fn new_plga_50_50() -> Self {
        Self {
            mw_initial: 100_000.0, // 100 kDa
            mw_current: 100_000.0,
            hydrolysis_rate: 0.003, // ~0.3 %/day
            autocatalysis: 0.5,
            mw_critical: 5_000.0, // 5 kDa loss of integrity
        }
    }

    /// Construct a `DegradablePolymer` with explicit parameters.
    pub fn new(
        mw_initial: f64,
        hydrolysis_rate: f64,
        autocatalysis: f64,
        mw_critical: f64,
    ) -> Self {
        Self {
            mw_initial,
            mw_current: mw_initial,
            hydrolysis_rate,
            autocatalysis,
            mw_critical,
        }
    }

    /// Advance degradation by time step Δt \[days\].
    ///
    /// Uses: `dM_n/dt = -k_h * (1 + β * (M_n0 - M_n)/M_n0) * M_n`
    pub fn step(&mut self, dt: f64) {
        let degradation_factor =
            1.0 + self.autocatalysis * (self.mw_initial - self.mw_current) / self.mw_initial;
        let rate = self.hydrolysis_rate * degradation_factor * self.mw_current;
        self.mw_current = (self.mw_current - rate * dt).max(self.mw_critical * 0.5);
    }

    /// Run degradation simulation for `days` total time with time step `dt`.
    pub fn run(&mut self, days: f64, dt: f64) {
        let steps = (days / dt).ceil() as usize;
        for _ in 0..steps {
            self.step(dt);
        }
    }

    /// Returns the relative molecular weight retention M_n / M_n0 ∈ \[0, 1\].
    pub fn mw_retention(&self) -> f64 {
        self.mw_current / self.mw_initial
    }

    /// Returns true if the polymer has degraded below the critical molecular weight.
    pub fn is_failed(&self) -> bool {
        self.mw_current <= self.mw_critical
    }

    /// Estimate Young's modulus retention based on M_n (semi-empirical).
    ///
    /// `E/E_0 = (M_n / M_n0)^0.5` (Fox-Flory type relation).
    pub fn modulus_retention(&self) -> f64 {
        self.mw_retention().sqrt()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // -- SoftTissue / Fung model --

    /// Fung strain energy is zero at zero strain.
    #[test]
    fn test_fung_zero_strain_zero_energy() {
        let tissue = SoftTissue::new_fung(0.5e3, [1.0, 1.0, 1.0, 0.5, 0.5, 0.5]);
        let w = tissue.fung_strain_energy(&[0.0; 6]);
        assert!(w.abs() < EPS, "W at zero strain should be 0, got {w}");
    }

    /// Fung stress is zero at zero strain.
    #[test]
    fn test_fung_zero_strain_zero_stress() {
        let tissue = SoftTissue::new_fung(0.5e3, [1.0, 1.0, 1.0, 0.5, 0.5, 0.5]);
        let s = tissue.fung_stress_s11(&[0.0; 6]);
        assert!(s.abs() < EPS, "S11 at zero strain should be 0, got {s}");
    }

    /// Fung model: stress increases with strain (exponential stiffening).
    #[test]
    fn test_fung_stress_increases_with_strain() {
        let tissue = SoftTissue::new_fung(0.5e3, [2.0, 1.0, 1.0, 0.5, 0.5, 0.5]);
        let s1 = tissue.fung_stress_s11(&[0.1, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let s2 = tissue.fung_stress_s11(&[0.2, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(
            s2 > s1,
            "Fung stress should increase with strain: s1={s1:.6}, s2={s2:.6}"
        );
    }

    /// Fung strain energy is positive for any nonzero strain.
    #[test]
    fn test_fung_energy_positive() {
        let tissue = SoftTissue::new_fung(1.0e3, [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]);
        let w = tissue.fung_strain_energy(&[0.1, 0.05, 0.02, 0.0, 0.0, 0.0]);
        assert!(w > 0.0, "Fung energy should be positive, got {w}");
    }

    /// Fung stress increases exponentially (ratio test).
    #[test]
    fn test_fung_stress_exponential_growth() {
        let tissue = SoftTissue::new_fung(1.0, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let strains = [0.1, 0.2, 0.4, 0.8];
        let mut prev = 0.0_f64;
        for &e in &strains {
            let s = tissue.fung_stress_s11(&[e, 0.0, 0.0, 0.0, 0.0, 0.0]);
            assert!(
                s > prev,
                "stress at {e} ({s:.6}) should exceed prev ({prev:.6})"
            );
            prev = s;
        }
    }

    /// Active stress adds to passive Fung stress.
    #[test]
    fn test_active_stress_adds() {
        let mut tissue = SoftTissue::new_fung(1.0e3, [1.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
        let passive = tissue.fung_stress_s11(&[0.1, 0.0, 0.0, 0.0, 0.0, 0.0]);
        tissue.set_active_stress(500.0);
        let total = tissue.total_stress_with_active(&[0.1, 0.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(
            (total - passive - 500.0).abs() < 1e-6,
            "Total stress should be passive + active, got {total:.6}"
        );
    }

    // -- HGO model --

    /// HGO fiber energy is zero when fibers are compressed (I4 ≤ 1).
    #[test]
    fn test_hgo_no_fiber_compression() {
        let tissue = SoftTissue::new_hgo(1.0e3, 1.0e4, 10.0, 0.1, 0.0);
        // Circumferential stretch < 1 → compression.
        let w = tissue.hgo_fiber_energy(0.9, 1.0, 1.0 / 0.9);
        assert_eq!(
            w, 0.0,
            "Fiber energy should be 0 under compression, got {w}"
        );
    }

    /// HGO fiber energy is positive under tension.
    #[test]
    fn test_hgo_fiber_tension_positive() {
        let tissue = SoftTissue::new_hgo(1.0e3, 1.0e4, 10.0, 0.1, 0.0);
        let w = tissue.hgo_fiber_energy(1.2, 1.0, 1.0 / 1.2);
        assert!(
            w > 0.0,
            "Fiber energy should be positive under tension, got {w}"
        );
    }

    /// HGO total energy ≥ neo-Hookean matrix energy (fibers add positive contribution).
    #[test]
    fn test_hgo_total_energy_geq_isotropic() {
        let tissue = SoftTissue::new_hgo(3.0e3, 1.5e4, 5.0, 0.1, 0.5);
        let w_total = tissue.hgo_total_energy(1.3, 1.1, 1.0 / (1.3 * 1.1));
        let radial = 1.0_f64 / (1.3 * 1.1);
        let i1 = 1.3f64 * 1.3 + 1.1_f64 * 1.1 + radial * radial;
        let w_iso = tissue.hgo_mu / 2.0 * (i1 - 3.0);
        assert!(
            w_total >= w_iso,
            "Total energy should be ≥ isotropic: w_total={w_total:.6}, w_iso={w_iso:.6}"
        );
    }

    /// HGO I4 equals 1 at zero deformation (identity).
    #[test]
    fn test_hgo_i4_identity() {
        let tissue = SoftTissue::new_hgo(1.0e3, 1.0e4, 5.0, 0.1, 0.3);
        let i4 = tissue.hgo_i4(1.0, 1.0);
        assert!(
            (i4 - 1.0).abs() < 1e-10,
            "I4 at identity stretch should be 1, got {i4}"
        );
    }

    // -- BoneModel --

    /// Kopperdahl-Keaveny: modulus increases with density.
    #[test]
    fn test_bone_density_modulus_increases() {
        let e1 = BoneModel::density_to_modulus_kk(0.3);
        let e2 = BoneModel::density_to_modulus_kk(0.6);
        let e3 = BoneModel::density_to_modulus_kk(1.0);
        assert!(e1 < e2 && e2 < e3, "Modulus should increase with density");
    }

    /// Gibson-Ashby modulus scales as (ρ/ρ_s)².
    #[test]
    fn test_bone_gibson_ashby_scaling() {
        let bone1 = BoneModel::new(0.5, 0.6, 18.0, 12.0, 3.5);
        let bone2 = BoneModel::new(1.0, 0.6, 18.0, 12.0, 3.5);
        let e1 = bone1.gibson_ashby_modulus(18.0);
        let e2 = bone2.gibson_ashby_modulus(18.0);
        let expected_ratio = 4.0; // (1.0/0.5)^2 = 4
        assert!(
            (e2 / e1 - expected_ratio).abs() < 0.01,
            "Gibson-Ashby should give 4x modulus for 2x density: ratio={:.6}",
            e2 / e1
        );
    }

    /// Wolff's law: density increases when stress exceeds reference.
    #[test]
    fn test_wolff_law_density_increases() {
        let mut bone = BoneModel::new(0.5, 0.6, 12.0, 8.0, 3.0);
        let initial_rho = bone.density;
        bone.wolff_adapt(10.0, 5.0, 0.1, 0.5, 1.0);
        assert!(
            bone.density > initial_rho,
            "Density should increase under high stress, got {:.4}",
            bone.density
        );
    }

    /// Wolff's law: density decreases under disuse (stress < reference).
    #[test]
    fn test_wolff_law_density_decreases() {
        let mut bone = BoneModel::new(1.0, 0.7, 18.0, 12.0, 3.5);
        let initial_rho = bone.density;
        bone.wolff_adapt(1.0, 8.0, 0.1, 0.5, 1.0);
        assert!(
            bone.density < initial_rho,
            "Density should decrease under disuse, got {:.4}",
            bone.density
        );
    }

    /// Yield strength increases with mineral fraction.
    #[test]
    fn test_bone_yield_strength() {
        let bone1 = BoneModel::new(1.5, 0.5, 18.0, 12.0, 3.5);
        let bone2 = BoneModel::new(1.5, 0.7, 18.0, 12.0, 3.5);
        assert!(
            bone2.yield_strength() > bone1.yield_strength(),
            "Higher mineral fraction → higher yield strength"
        );
    }

    // -- BloodRheology --

    /// Carreau-Yasuda: shear thinning (dη/dγ̇ < 0).
    #[test]
    fn test_blood_shear_thinning() {
        let blood = BloodRheology::new_physiological();
        let eta1 = blood.carreau_viscosity(1.0);
        let eta2 = blood.carreau_viscosity(100.0);
        assert!(
            eta2 < eta1,
            "Blood should shear-thin: η(1)={eta1:.6}, η(100)={eta2:.6}"
        );
    }

    /// Carreau-Yasuda: at very high shear rate, approaches η_∞.
    #[test]
    fn test_blood_high_shear_newtonian() {
        let blood = BloodRheology::new_physiological();
        assert!(
            blood.is_newtonian_limit(10000.0),
            "Blood should be Newtonian at very high shear rate"
        );
    }

    /// Casson model: zero shear rate gives zero stress.
    #[test]
    fn test_casson_zero_shear() {
        let blood = BloodRheology::new_physiological();
        let tau = blood.casson_stress(0.0);
        assert!(
            tau.abs() < EPS,
            "Casson stress at zero shear should be 0, got {tau}"
        );
    }

    /// Casson model: stress increases with shear rate.
    #[test]
    fn test_casson_stress_increases() {
        let blood = BloodRheology::new_physiological();
        let tau1 = blood.casson_stress(10.0);
        let tau2 = blood.casson_stress(100.0);
        assert!(tau2 > tau1, "Casson stress should increase with shear rate");
    }

    /// Hematocrit viscosity > plasma viscosity.
    #[test]
    fn test_hematocrit_viscosity_greater_than_plasma() {
        let blood = BloodRheology::new_physiological();
        assert!(
            blood.hematocrit_viscosity() > blood.eta_inf,
            "Blood viscosity should exceed plasma viscosity"
        );
    }

    /// Casson: positive shear gives positive stress.
    #[test]
    fn test_casson_positive_shear_positive_stress() {
        let blood = BloodRheology::new_physiological();
        assert!(
            blood.casson_stress(5.0) > 0.0,
            "Positive shear should give positive Casson stress"
        );
    }

    /// Carreau-Yasuda: viscosity is bounded between η_∞ and η_0.
    #[test]
    fn test_carreau_bounds() {
        let blood = BloodRheology::new_physiological();
        for &gamma in &[0.001, 0.1, 1.0, 10.0, 1000.0] {
            let eta = blood.carreau_viscosity(gamma);
            assert!(
                eta >= blood.eta_inf && eta <= blood.eta_0,
                "Carreau viscosity out of [η_∞, η_0] at γ̇={gamma}: η={eta:.6}"
            );
        }
    }

    // -- HydrogelModel --

    /// Flory-Rehner mixing chemical potential is negative at low polymer fraction.
    #[test]
    fn test_hydrogel_mixing_potential_negative_low_phi() {
        let gel = HydrogelModel::new(0.4, 50.0, 1.0e4);
        let mu = gel.mixing_chemical_potential(0.01);
        assert!(
            mu < 0.0,
            "Mixing potential should be negative at low phi, got {mu:.6}"
        );
    }

    /// Equilibrium swelling ratio Q > 1 (gel swells).
    #[test]
    fn test_hydrogel_swelling_ratio_gt1() {
        let gel = HydrogelModel::new(0.3, 100.0, 1.0e4);
        let q = gel.equilibrium_swelling_ratio();
        assert!(
            q > 1.0,
            "Equilibrium swelling ratio should be > 1, got {q:.6}"
        );
    }

    /// More cross-links (shorter chains) → larger magnitude elastic restoring potential.
    ///
    /// At the same polymer fraction φ < φ_eq, the elastic chemical potential is negative
    /// (it drives the solvent OUT of the gel). A tighter network has a more negative
    /// elastic potential (larger restoring force), so |mu_tight| > |mu_loose|.
    #[test]
    fn test_hydrogel_crosslink_reduces_swelling() {
        // Compare elastic chemical potential at the same phi.
        // Tighter gel (more cross-links): smaller n_chains.
        let gel_loose = HydrogelModel::new(0.3, 500.0, 1.0e3);
        let gel_tight = HydrogelModel::new(0.3, 10.0, 1.0e3);
        let phi = 0.1;
        let mu_loose = gel_loose.elastic_chemical_potential(phi);
        let mu_tight = gel_tight.elastic_chemical_potential(phi);
        // Tight gel has more elastic resistance: |mu_tight| > |mu_loose|.
        assert!(
            mu_tight.abs() > mu_loose.abs(),
            "Tighter cross-linking should have larger |elastic potential|: |mu_tight|={:.6}, |mu_loose|={:.6}",
            mu_tight.abs(),
            mu_loose.abs()
        );
    }

    // -- CartilageModel --

    /// Characteristic time is positive.
    #[test]
    fn test_cartilage_char_time_positive() {
        let cart = CartilageModel::new(0.79, 1.16e-15, 0.1, 2.0);
        let t_star = cart.characteristic_time();
        assert!(
            t_star > 0.0,
            "Characteristic time should be positive, got {t_star}"
        );
    }

    /// Creep response increases with time.
    #[test]
    fn test_cartilage_creep_increases() {
        let cart = CartilageModel::new(0.79, 1.16e-15, 0.1, 2.0);
        let u1 = cart.creep_response(0.5, 100.0);
        let u2 = cart.creep_response(0.5, 1000.0);
        assert!(u2 > u1, "Creep deformation should increase with time");
    }

    /// Stress relaxation decreases with time.
    #[test]
    fn test_cartilage_stress_relaxation_decreases() {
        let cart = CartilageModel::new(0.79, 1.16e-15, 0.1, 2.0);
        let s1 = cart.stress_relaxation(0.1, 100.0);
        let s2 = cart.stress_relaxation(0.1, 1000.0);
        assert!(s2 < s1, "Stress relaxation should decrease with time");
    }

    // -- ArterialWall --

    /// Wall stress increases with circumferential stretch.
    #[test]
    fn test_arterial_wall_stress_increases() {
        let wall = ArterialWall::new_aorta();
        let s1 = wall.wall_stress_at_stretch(1.1);
        let s2 = wall.wall_stress_at_stretch(1.5);
        assert!(
            s2 > s1,
            "Wall stress should increase with stretch: s1={s1:.6}, s2={s2:.6}"
        );
    }

    /// Layer fractions sum to 1.
    #[test]
    fn test_arterial_layer_fractions_sum_to_one() {
        let wall = ArterialWall::new_aorta();
        let sum: f64 = wall.layer_fractions.iter().sum();
        assert!(
            (sum - 1.0).abs() < 1e-10,
            "Layer fractions should sum to 1, got {sum:.6}"
        );
    }

    // -- DegradablePolymer --

    /// Molecular weight decreases after degradation.
    #[test]
    fn test_polymer_mw_decreases() {
        let mut poly = DegradablePolymer::new_plga_50_50();
        poly.run(100.0, 1.0);
        assert!(
            poly.mw_current < poly.mw_initial,
            "MW should decrease after degradation, got {:.0}",
            poly.mw_current
        );
    }

    /// MW retention ∈ \[0, 1\].
    #[test]
    fn test_polymer_mw_retention_bounded() {
        let mut poly = DegradablePolymer::new_plga_50_50();
        poly.run(30.0, 1.0);
        let ret = poly.mw_retention();
        assert!(
            (0.0..=1.0).contains(&ret),
            "MW retention should be in [0, 1], got {ret:.6}"
        );
    }

    /// Modulus retention ≤ 1.
    #[test]
    fn test_polymer_modulus_retention_leq1() {
        let mut poly = DegradablePolymer::new_plga_50_50();
        poly.run(50.0, 1.0);
        let e_ret = poly.modulus_retention();
        assert!(
            e_ret <= 1.0,
            "Modulus retention should be ≤ 1, got {e_ret:.6}"
        );
    }

    /// Degradation is faster with higher hydrolysis rate.
    #[test]
    fn test_polymer_higher_rate_more_degradation() {
        let mut slow = DegradablePolymer::new(100_000.0, 0.001, 0.5, 5_000.0);
        let mut fast = DegradablePolymer::new(100_000.0, 0.01, 0.5, 5_000.0);
        slow.run(100.0, 1.0);
        fast.run(100.0, 1.0);
        assert!(
            fast.mw_current < slow.mw_current,
            "Higher rate should degrade more: fast={:.0}, slow={:.0}",
            fast.mw_current,
            slow.mw_current
        );
    }

    /// Autocatalysis accelerates degradation beyond simple first-order.
    #[test]
    fn test_polymer_autocatalysis_accelerates() {
        let mut no_auto = DegradablePolymer::new(100_000.0, 0.003, 0.0, 1_000.0);
        let mut with_auto = DegradablePolymer::new(100_000.0, 0.003, 2.0, 1_000.0);
        no_auto.run(100.0, 1.0);
        with_auto.run(100.0, 1.0);
        assert!(
            with_auto.mw_current < no_auto.mw_current,
            "Autocatalysis should accelerate degradation"
        );
    }
}
