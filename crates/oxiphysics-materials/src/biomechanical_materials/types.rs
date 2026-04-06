//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
/// Biphasic model for articular cartilage.
///
/// Cartilage is modeled as a mixture of a solid matrix (proteoglycan-collagen
/// network) and an interstitial fluid phase. The model includes:
///
/// - Aggregate modulus H_A for the solid phase
/// - Hydraulic permeability k for fluid flow
/// - Poisson's ratio for the solid matrix
/// - Depth-dependent material properties (surface/middle/deep zones)
#[derive(Debug, Clone)]
pub struct CartilageBiphasic {
    /// Aggregate modulus (Pa) – stiffness in confined compression.
    pub aggregate_modulus: f64,
    /// Hydraulic permeability (m^4 / N s).
    pub permeability: f64,
    /// Solid matrix Poisson's ratio.
    pub poisson_ratio: f64,
    /// Solid volume fraction (porosity = 1 - phi_s).
    pub solid_fraction: f64,
    /// Young's modulus of the solid matrix (Pa).
    pub youngs_modulus: f64,
    /// Strain-dependent permeability parameter M.
    pub permeability_strain_coeff: f64,
    /// Current fluid pressure (Pa).
    pub(super) fluid_pressure: f64,
}
impl CartilageBiphasic {
    /// Create a new biphasic cartilage model.
    pub fn new(
        aggregate_modulus: f64,
        permeability: f64,
        poisson_ratio: f64,
        solid_fraction: f64,
    ) -> Self {
        let nu = poisson_ratio;
        let ha = aggregate_modulus;
        let youngs = ha * (1.0 + nu) * (1.0 - 2.0 * nu) / (1.0 - nu);
        Self {
            aggregate_modulus,
            permeability,
            poisson_ratio,
            solid_fraction,
            youngs_modulus: youngs,
            permeability_strain_coeff: 5.0,
            fluid_pressure: 0.0,
        }
    }
    /// Create a model for the superficial zone of cartilage.
    pub fn superficial_zone() -> Self {
        Self::new(0.42e6, 4.0e-15, 0.1, 0.15)
    }
    /// Create a model for the middle zone of cartilage.
    pub fn middle_zone() -> Self {
        Self::new(0.60e6, 2.0e-15, 0.25, 0.20)
    }
    /// Create a model for the deep zone of cartilage.
    pub fn deep_zone() -> Self {
        Self::new(0.85e6, 1.0e-15, 0.35, 0.25)
    }
    /// Compute the strain-dependent permeability.
    ///
    /// k(e) = k0 * exp(M * e_vol) where e_vol is volumetric strain
    pub fn current_permeability(&self, volumetric_strain: f64) -> f64 {
        self.permeability * (self.permeability_strain_coeff * volumetric_strain).exp()
    }
    /// Compute the effective (drained) stress of the solid matrix.
    pub fn effective_stress(&self, strain: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let e_mod = self.youngs_modulus;
        let nu = self.poisson_ratio;
        let lambda = e_mod * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e_mod / (2.0 * (1.0 + nu));
        let tr_e = trace3(strain);
        let mut sigma = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                sigma[i][j] = 2.0 * mu * strain[i][j];
                if i == j {
                    sigma[i][j] += lambda * tr_e;
                }
            }
        }
        sigma
    }
    /// Compute the total stress (effective + fluid pressure).
    #[allow(clippy::needless_range_loop)]
    pub fn total_stress(&self, strain: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut sigma = self.effective_stress(strain);
        for i in 0..3 {
            sigma[i][i] -= self.fluid_pressure;
        }
        sigma
    }
    /// Update the fluid pressure based on volumetric strain rate.
    ///
    /// dp/dt = -H_A * d(e_vol)/dt + (k * H_A / L^2) * p
    /// Simplified update for a given characteristic length L.
    pub fn update_pressure(
        &mut self,
        volumetric_strain_rate: f64,
        characteristic_length: f64,
        dt: f64,
    ) -> f64 {
        let k = self.permeability;
        let ha = self.aggregate_modulus;
        let l2 = characteristic_length * characteristic_length;
        let diffusivity = k * ha / l2;
        let dp = (-ha * volumetric_strain_rate - diffusivity * self.fluid_pressure) * dt;
        self.fluid_pressure += dp;
        self.fluid_pressure
    }
    /// Compute the consolidation time constant.
    ///
    /// tau = L^2 / (H_A * k)
    pub fn consolidation_time(&self, characteristic_length: f64) -> f64 {
        let l2 = characteristic_length * characteristic_length;
        l2 / (self.aggregate_modulus * self.permeability)
    }
    /// Compute the porosity (fluid volume fraction).
    pub fn porosity(&self) -> f64 {
        1.0 - self.solid_fraction
    }
    /// Get the current fluid pressure.
    pub fn fluid_pressure(&self) -> f64 {
        self.fluid_pressure
    }
    /// Set the fluid pressure directly.
    pub fn set_fluid_pressure(&mut self, p: f64) {
        self.fluid_pressure = p;
    }
}
/// Fung quasi-linear viscoelastic (QLV) model for soft tissues.
///
/// The stress response is a convolution of the elastic response with a
/// reduced relaxation function:
///   sigma(t) = integral_0^t G(t-tau) * dS^e/dtau dtau
///
/// The relaxation function uses a continuous spectrum:
///   G(t) = 1 + c * \[E1(t/tau2) - E1(t/tau1)\] / ln(tau2/tau1)
///
/// where E1 is the exponential integral. For computational efficiency, we
/// use a Prony series approximation.
#[derive(Debug, Clone)]
pub struct QuasiLinearViscoelastic {
    /// Elastic response model (Fung).
    pub elastic: FungHyperelastic,
    /// Prony series coefficients (relaxation moduli).
    pub g_coeffs: Vec<f64>,
    /// Prony series time constants (seconds).
    pub tau_coeffs: Vec<f64>,
    /// Long-term (equilibrium) relaxation ratio.
    pub g_inf: f64,
    /// Internal state variables for each Prony term.
    pub(super) state: Vec<f64>,
    /// Previous elastic stress (scalar representation).
    pub(super) prev_stress: f64,
}
impl QuasiLinearViscoelastic {
    /// Create a new QLV model.
    pub fn new(
        elastic: FungHyperelastic,
        g_coeffs: Vec<f64>,
        tau_coeffs: Vec<f64>,
        g_inf: f64,
    ) -> Self {
        let n = g_coeffs.len();
        Self {
            elastic,
            g_coeffs,
            tau_coeffs,
            g_inf,
            state: vec![0.0; n],
            prev_stress: 0.0,
        }
    }
    /// Create a default QLV model for tendon.
    pub fn tendon_default() -> Self {
        let elastic = FungHyperelastic::isotropic(200.0, 30.0);
        Self::new(elastic, vec![0.3, 0.15, 0.05], vec![1.0, 10.0, 100.0], 0.5)
    }
    /// Reset all internal state variables.
    pub fn reset(&mut self) {
        for s in &mut self.state {
            *s = 0.0;
        }
        self.prev_stress = 0.0;
    }
    /// Update the viscoelastic stress for a given deformation gradient and time step.
    ///
    /// Returns the total (viscoelastic) stress scalar.
    pub fn update(&mut self, f: &[[f64; 3]; 3], dt: f64) -> f64 {
        let e = green_lagrange(f);
        let q = self.elastic.compute_q(&e);
        let elastic_stress = self.elastic.c * q.exp() * q;
        let d_stress = elastic_stress - self.prev_stress;
        for (i, state_i) in self.state.iter_mut().enumerate() {
            let tau = self.tau_coeffs[i];
            let exp_dt = (-dt / tau).exp();
            *state_i = exp_dt * (*state_i) + self.g_coeffs[i] * d_stress;
        }
        let viscoelastic_stress = self.g_inf * elastic_stress + self.state.iter().sum::<f64>();
        self.prev_stress = elastic_stress;
        viscoelastic_stress
    }
    /// Compute the relaxation function G(t) at a given time.
    pub fn relaxation_function(&self, t: f64) -> f64 {
        let mut g = self.g_inf;
        for (i, &gi) in self.g_coeffs.iter().enumerate() {
            g += gi * (-t / self.tau_coeffs[i]).exp();
        }
        g
    }
    /// Compute the instantaneous (glassy) modulus ratio.
    pub fn instantaneous_modulus_ratio(&self) -> f64 {
        self.g_inf + self.g_coeffs.iter().sum::<f64>()
    }
}
/// Multi-layer tissue model (e.g., skin + subcutaneous fat + muscle).
#[derive(Debug, Clone)]
pub struct MultiLayerTissue {
    /// The tissue layers from outermost to innermost.
    pub layers: Vec<TissueLayer>,
}
impl MultiLayerTissue {
    /// Create a new multi-layer tissue model.
    pub fn new(layers: Vec<TissueLayer>) -> Self {
        Self { layers }
    }
    /// Create a typical skin-fat-muscle layered model.
    pub fn skin_fat_muscle() -> Self {
        Self::new(vec![
            TissueLayer::new(TissueLayerType::Epithelial, 0.002, 1100.0, 0.6e6, 0.49),
            TissueLayer::new(TissueLayerType::Adipose, 0.010, 900.0, 0.003e6, 0.49),
            TissueLayer::new(TissueLayerType::Muscular, 0.030, 1060.0, 0.5e6, 0.45),
        ])
    }
    /// Total thickness of all layers.
    pub fn total_thickness(&self) -> f64 {
        self.layers.iter().map(|l| l.thickness).sum()
    }
    /// Total mass per unit area (kg/m^2).
    pub fn mass_per_area(&self) -> f64 {
        self.layers.iter().map(|l| l.density * l.thickness).sum()
    }
    /// Effective Young's modulus (Voigt upper bound for in-plane loading).
    pub fn effective_youngs_modulus_voigt(&self) -> f64 {
        let total_t = self.total_thickness();
        if total_t < 1e-15 {
            return 0.0;
        }
        self.layers
            .iter()
            .map(|l| l.youngs_modulus * l.thickness)
            .sum::<f64>()
            / total_t
    }
    /// Effective Young's modulus (Reuss lower bound for through-thickness loading).
    pub fn effective_youngs_modulus_reuss(&self) -> f64 {
        let total_t = self.total_thickness();
        if total_t < 1e-15 {
            return 0.0;
        }
        let inv_sum: f64 = self
            .layers
            .iter()
            .map(|l| {
                if l.youngs_modulus > 1e-30 {
                    l.thickness / l.youngs_modulus
                } else {
                    0.0
                }
            })
            .sum();
        if inv_sum > 1e-30 {
            total_t / inv_sum
        } else {
            0.0
        }
    }
    /// Find the depth of a specific layer type (distance from surface).
    pub fn layer_depth(&self, layer_type: TissueLayerType) -> Option<f64> {
        let mut depth = 0.0;
        for l in &self.layers {
            if l.layer_type == layer_type {
                return Some(depth);
            }
            depth += l.thickness;
        }
        None
    }
}
/// Bone remodeling model based on Wolff's law.
///
/// Bone adapts its density in response to mechanical loading. The model
/// follows the strain energy density stimulus approach:
///   d(rho)/dt = B * (S/S_ref - 1) when |S/S_ref - 1| > dead_zone
///
/// where S is the strain energy density stimulus, S_ref is the reference
/// (homeostatic) stimulus, and B is the remodeling rate.
#[derive(Debug, Clone)]
pub struct BoneRemodeling {
    /// Current bone density (kg/m^3).
    pub density: f64,
    /// Minimum density (fully resorbed cortical bone).
    pub density_min: f64,
    /// Maximum density (dense cortical bone).
    pub density_max: f64,
    /// Remodeling rate constant (kg/m^3/s).
    pub remodeling_rate: f64,
    /// Reference (homeostatic) strain energy density (Pa).
    pub stimulus_ref: f64,
    /// Dead zone half-width (fractional).
    pub dead_zone: f64,
    /// Power-law exponent relating density to elastic modulus.
    pub density_exponent: f64,
    /// Reference elastic modulus at reference density (Pa).
    pub e_ref: f64,
    /// Reference density for modulus scaling (kg/m^3).
    pub rho_ref: f64,
}
impl BoneRemodeling {
    /// Create a new bone remodeling model.
    pub fn new(
        density: f64,
        density_min: f64,
        density_max: f64,
        remodeling_rate: f64,
        stimulus_ref: f64,
        dead_zone: f64,
    ) -> Self {
        Self {
            density,
            density_min,
            density_max,
            remodeling_rate,
            stimulus_ref,
            dead_zone,
            density_exponent: 2.0,
            e_ref: 17.0e9,
            rho_ref: 1900.0,
        }
    }
    /// Create a default cortical bone model.
    pub fn cortical_bone() -> Self {
        Self::new(1900.0, 500.0, 2100.0, 0.1, 0.004, 0.1)
    }
    /// Create a default trabecular (spongy) bone model.
    pub fn trabecular_bone() -> Self {
        Self::new(800.0, 100.0, 1500.0, 0.15, 0.002, 0.15)
    }
    /// Compute the current elastic modulus based on density.
    ///
    /// E = E_ref * (rho/rho_ref)^n
    pub fn elastic_modulus(&self) -> f64 {
        self.e_ref * (self.density / self.rho_ref).powf(self.density_exponent)
    }
    /// Compute the strain energy density for a given strain state.
    pub fn strain_energy_density(&self, strain: &[[f64; 3]; 3]) -> f64 {
        let e_mod = self.elastic_modulus();
        let nu = 0.3;
        let lambda = e_mod * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e_mod / (2.0 * (1.0 + nu));
        let tr_e = trace3(strain);
        let e2 = double_contract(strain, strain);
        0.5 * lambda * tr_e * tr_e + mu * e2
    }
    /// Update bone density based on the current loading stimulus.
    ///
    /// Returns the new density.
    pub fn update_density(&mut self, stimulus: f64, dt: f64) -> f64 {
        let s_ratio = stimulus / self.stimulus_ref;
        let deviation = s_ratio - 1.0;
        if deviation.abs() > self.dead_zone {
            let effective_deviation = if deviation > 0.0 {
                deviation - self.dead_zone
            } else {
                deviation + self.dead_zone
            };
            let d_rho = self.remodeling_rate * effective_deviation * dt;
            self.density = clamp(self.density + d_rho, self.density_min, self.density_max);
        }
        self.density
    }
    /// Compute the bone quality index (0 = fully resorbed, 1 = maximum density).
    pub fn quality_index(&self) -> f64 {
        (self.density - self.density_min) / (self.density_max - self.density_min)
    }
    /// Check if the bone is osteoporotic (density below threshold).
    pub fn is_osteoporotic(&self) -> bool {
        self.density < 0.5 * self.rho_ref
    }
    /// Compute the yield stress based on density (empirical relation).
    pub fn yield_stress(&self) -> f64 {
        let rho_gcc = self.density / 1000.0;
        137.0e6 * rho_gcc.powf(1.88)
    }
}
/// Holzapfel-Gasser-Ogden anisotropic model for arterial walls.
///
/// The strain energy is:
///   W = (mu/2)*(I1-3) + sum_i (k1/(2*k2)) * (exp(k2*`E_i`^2) - 1)
///
/// where E_i = kappa*(I1-3) + (1-3*kappa)*(I4i-1) and `x` = max(0,x).
///
/// Two fiber families are supported, defined by direction vectors a1, a2.
#[derive(Debug, Clone)]
pub struct HolzapfelGasserOgden {
    /// Ground substance shear modulus (Pa).
    pub mu: f64,
    /// Fiber stiffness parameter k1 (Pa).
    pub k1: f64,
    /// Fiber exponential parameter k2 (dimensionless).
    pub k2: f64,
    /// Dispersion parameter kappa in \[0, 1/3\].
    pub kappa: f64,
    /// First fiber family direction (unit vector).
    pub a1: [f64; 3],
    /// Second fiber family direction (unit vector).
    pub a2: [f64; 3],
    /// Bulk modulus for volumetric penalty (Pa).
    pub bulk_modulus: f64,
}
impl HolzapfelGasserOgden {
    /// Create a new HGO model.
    pub fn new(
        mu: f64,
        k1: f64,
        k2: f64,
        kappa: f64,
        a1: [f64; 3],
        a2: [f64; 3],
        bulk_modulus: f64,
    ) -> Self {
        let n1 = norm3(&a1);
        let n2 = norm3(&a2);
        let a1n = if n1 > 1e-15 {
            scale3(&a1, 1.0 / n1)
        } else {
            a1
        };
        let a2n = if n2 > 1e-15 {
            scale3(&a2, 1.0 / n2)
        } else {
            a2
        };
        Self {
            mu,
            k1,
            k2,
            kappa: kappa.clamp(0.0, 1.0 / 3.0),
            a1: a1n,
            a2: a2n,
            bulk_modulus,
        }
    }
    /// Human aortic adventitia (Holzapfel 2000, approximate).
    ///
    /// c₁ = 7.64 kPa, k₁ = 996.6 Pa, k₂ = 524.6, κ = 0.226
    pub fn aorta_adventitia() -> Self {
        Self::new(
            7640.0,
            996.6,
            524.6,
            0.226,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            100.0e3,
        )
    }
    /// Human aortic media (approximate).
    pub fn aorta_media() -> Self {
        Self::new(
            3000.0,
            2362.0,
            100.7,
            0.1,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            100.0e3,
        )
    }
    /// Create a model for a typical arterial wall.
    pub fn arterial_wall() -> Self {
        let angle = 39.0_f64.to_radians();
        Self::new(
            15.0e3,
            2.3632e3,
            0.8393,
            0.226,
            [angle.cos(), angle.sin(), 0.0],
            [angle.cos(), -angle.sin(), 0.0],
            100.0e3,
        )
    }
    /// Compute pseudo-invariant I4 = a . C . a for a fiber direction.
    fn compute_i4(&self, c: &[[f64; 3]; 3], a: &[f64; 3]) -> f64 {
        let mut i4 = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                i4 += a[i] * c[i][j] * a[j];
            }
        }
        i4
    }
    /// Compute the strain energy density (alias for `strain_energy`).
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        self.strain_energy(f)
    }
    /// Compute the anisotropic strain-like quantity E_i.
    fn fiber_strain(&self, i1: f64, i4: f64) -> f64 {
        let e = self.kappa * (i1 - 3.0) + (1.0 - 3.0 * self.kappa) * (i4 - 1.0);
        e.max(0.0)
    }
    /// Compute the strain energy density.
    pub fn strain_energy(&self, f: &[[f64; 3]; 3]) -> f64 {
        let c = right_cauchy_green(f);
        let j = jacobian(f);
        let i1 = invariant_i1(&c);
        let i4_1 = self.compute_i4(&c, &self.a1);
        let i4_2 = self.compute_i4(&c, &self.a2);
        let j_23 = j.powf(-2.0 / 3.0);
        let i1_bar = j_23 * i1;
        let w_iso = 0.5 * self.mu * (i1_bar - 3.0);
        let e1 = self.fiber_strain(i1_bar, i4_1);
        let e2 = self.fiber_strain(i1_bar, i4_2);
        let w_fib = if self.k2.abs() > 1e-30 {
            (self.k1 / (2.0 * self.k2))
                * ((self.k2 * e1 * e1).exp() - 1.0 + (self.k2 * e2 * e2).exp() - 1.0)
        } else {
            self.k1 * (e1 * e1 + e2 * e2)
        };
        let w_vol = 0.5 * self.bulk_modulus * (j - 1.0) * (j - 1.0);
        w_iso + w_fib + w_vol
    }
    /// Compute the isotropic contribution to the Cauchy stress.
    pub fn isotropic_cauchy_stress(&self, f: &[[f64; 3]; 3]) -> f64 {
        let c = right_cauchy_green(f);
        let j = jacobian(f);
        let j_23 = j.powf(-2.0 / 3.0);
        let _i1_bar = j_23 * invariant_i1(&c);
        self.mu * j_23 / j
    }
    /// Compute the fiber stress contribution for a single fiber family.
    pub fn fiber_stress(&self, f: &[[f64; 3]; 3], fiber_index: usize) -> f64 {
        let c = right_cauchy_green(f);
        let j = jacobian(f);
        let j_23 = j.powf(-2.0 / 3.0);
        let i1_bar = j_23 * invariant_i1(&c);
        let a = if fiber_index == 0 { &self.a1 } else { &self.a2 };
        let i4 = self.compute_i4(&c, a);
        let e = self.fiber_strain(i1_bar, i4);
        if e > 0.0 {
            2.0 * self.k1 * e * (self.k2 * e * e).exp()
        } else {
            0.0
        }
    }
    /// Compute the volumetric penalty stress (hydrostatic pressure).
    pub fn volumetric_stress(&self, f: &[[f64; 3]; 3]) -> f64 {
        let j = jacobian(f);
        self.bulk_modulus * (j - 1.0)
    }
}
/// A single tissue layer with thickness and material properties.
#[derive(Debug, Clone)]
pub struct TissueLayer {
    /// Layer type.
    pub layer_type: TissueLayerType,
    /// Thickness (m).
    pub thickness: f64,
    /// Density (kg/m^3).
    pub density: f64,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}
impl TissueLayer {
    /// Create a new tissue layer.
    pub fn new(
        layer_type: TissueLayerType,
        thickness: f64,
        density: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
    ) -> Self {
        Self {
            layer_type,
            thickness,
            density,
            youngs_modulus,
            poisson_ratio,
        }
    }
    /// Compute the shear modulus.
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
    /// Compute the bulk modulus.
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
    /// Compute the P-wave modulus (constrained modulus).
    pub fn constrained_modulus(&self) -> f64 {
        let nu = self.poisson_ratio;
        self.youngs_modulus * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }
}
/// Tissue damage model using a Kachanov-type continuum damage variable.
#[derive(Debug, Clone)]
pub struct TissueDamage {
    /// Damage variable D in \[0, 1\] (0 = intact, 1 = fully damaged).
    pub damage: f64,
    /// Damage threshold strain.
    pub threshold_strain: f64,
    /// Damage evolution rate.
    pub damage_rate: f64,
    /// Maximum historical strain (Weibull-type).
    pub(super) max_strain: f64,
}
impl TissueDamage {
    /// Create a new tissue damage model.
    pub fn new(threshold_strain: f64, damage_rate: f64) -> Self {
        Self {
            damage: 0.0,
            threshold_strain,
            damage_rate,
            max_strain: 0.0,
        }
    }
    /// Default parameters for soft tissue rupture.
    pub fn soft_tissue_default() -> Self {
        Self::new(0.3, 2.0)
    }
    /// Update damage based on current effective strain.
    pub fn update(&mut self, effective_strain: f64) -> f64 {
        if effective_strain > self.max_strain {
            self.max_strain = effective_strain;
        }
        if self.max_strain > self.threshold_strain {
            let excess = self.max_strain - self.threshold_strain;
            self.damage = 1.0 - (-self.damage_rate * excess).exp();
        }
        self.damage
    }
    /// Compute the effective (damaged) stress from an undamaged stress.
    pub fn effective_stress(&self, undamaged_stress: f64) -> f64 {
        (1.0 - self.damage) * undamaged_stress
    }
    /// Check if the tissue is ruptured (damage > 0.99).
    pub fn is_ruptured(&self) -> bool {
        self.damage > 0.99
    }
    /// Reset the damage state.
    pub fn reset(&mut self) {
        self.damage = 0.0;
        self.max_strain = 0.0;
    }
}
/// Multi-term Ogden model for skin mechanics.
///
/// The strain energy density is:
///   W = sum_p (mu_p / alpha_p) * (lambda1^alpha_p + lambda2^alpha_p + lambda3^alpha_p - 3)
///     + (K/2)*(J-1)^2
///
/// where lambda_i are the principal stretches.
///
/// The Ogden model is particularly effective for skin because it captures:
/// - Toe region at low strains (collagen recruitment)
/// - Stiffening at large strains
/// - Anisotropy through directional parameters
#[derive(Debug, Clone)]
pub struct SkinMechanicsOgden {
    /// Shear moduli mu_p (Pa).
    pub mu: Vec<f64>,
    /// Exponents alpha_p (dimensionless).
    pub alpha: Vec<f64>,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
    /// Langer's line direction (unit vector in reference configuration).
    pub langer_direction: [f64; 3],
    /// Anisotropy ratio (>1 means stiffer along Langer's lines).
    pub anisotropy_ratio: f64,
    /// Pre-stress (residual stress in skin, Pa).
    pub pre_stress: f64,
}
impl SkinMechanicsOgden {
    /// Create a new Ogden skin model.
    pub fn new(mu: Vec<f64>, alpha: Vec<f64>, bulk_modulus: f64) -> Self {
        assert_eq!(mu.len(), alpha.len(), "mu and alpha must have same length");
        Self {
            mu,
            alpha,
            bulk_modulus,
            langer_direction: [1.0, 0.0, 0.0],
            anisotropy_ratio: 1.0,
            pre_stress: 0.0,
        }
    }
    /// Create a model for adult human skin (forearm region).
    pub fn human_skin() -> Self {
        Self {
            mu: vec![0.12e6, -0.05e6],
            alpha: vec![10.0, -2.0],
            bulk_modulus: 20.0e6,
            langer_direction: [1.0, 0.0, 0.0],
            anisotropy_ratio: 1.3,
            pre_stress: 1.0e3,
        }
    }
    /// Create a model for aged skin (reduced stiffness).
    pub fn aged_skin() -> Self {
        Self {
            mu: vec![0.08e6, -0.03e6],
            alpha: vec![8.0, -1.5],
            bulk_modulus: 15.0e6,
            langer_direction: [1.0, 0.0, 0.0],
            anisotropy_ratio: 1.1,
            pre_stress: 0.5e3,
        }
    }
    /// Compute the strain energy density from principal stretches.
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let j_inv_third = j.powf(-1.0 / 3.0);
        let l1 = lambda1 * j_inv_third;
        let l2 = lambda2 * j_inv_third;
        let l3 = lambda3 * j_inv_third;
        let mut w_dev = 0.0;
        for p in 0..self.mu.len() {
            let ap = self.alpha[p];
            let mp = self.mu[p];
            w_dev += (mp / ap) * (l1.powf(ap) + l2.powf(ap) + l3.powf(ap) - 3.0);
        }
        let w_vol = 0.5 * self.bulk_modulus * (j - 1.0) * (j - 1.0);
        w_dev + w_vol
    }
    /// Compute the strain energy density from a deformation gradient.
    pub fn strain_energy(&self, f: &[[f64; 3]; 3]) -> f64 {
        let stretches = self.principal_stretches(f);
        self.strain_energy_from_stretches(stretches[0], stretches[1], stretches[2])
    }
    /// Compute the principal Cauchy stress from principal stretches.
    pub fn principal_cauchy_stress(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> [f64; 3] {
        let j = lambda1 * lambda2 * lambda3;
        let j_inv_third = j.powf(-1.0 / 3.0);
        let l = [
            lambda1 * j_inv_third,
            lambda2 * j_inv_third,
            lambda3 * j_inv_third,
        ];
        let mut sigma_dev = [0.0; 3];
        for p in 0..self.mu.len() {
            let ap = self.alpha[p];
            let mp = self.mu[p];
            let mean_pow: f64 = (l[0].powf(ap) + l[1].powf(ap) + l[2].powf(ap)) / 3.0;
            for i in 0..3 {
                sigma_dev[i] += (mp / j) * (l[i].powf(ap) - mean_pow);
            }
        }
        let p = self.bulk_modulus * (j - 1.0);
        [
            sigma_dev[0] + p + self.pre_stress,
            sigma_dev[1] + p + self.pre_stress,
            sigma_dev[2] + p,
        ]
    }
    /// Compute principal stretches from deformation gradient (simple SVD-like).
    pub fn principal_stretches(&self, f: &[[f64; 3]; 3]) -> [f64; 3] {
        let c = right_cauchy_green(f);
        let i1 = invariant_i1(&c);
        let i2 = invariant_i2(&c);
        let i3 = invariant_i3(&c);
        let eigenvalues = solve_cubic_symmetric(i1, i2, i3);
        [
            eigenvalues[0].max(1e-30).sqrt(),
            eigenvalues[1].max(1e-30).sqrt(),
            eigenvalues[2].max(1e-30).sqrt(),
        ]
    }
    /// Compute the initial shear modulus.
    pub fn initial_shear_modulus(&self) -> f64 {
        let mut mu_total = 0.0;
        for p in 0..self.mu.len() {
            mu_total += 0.5 * self.mu[p] * self.alpha[p];
        }
        mu_total
    }
    /// Compute the skin tension along a given direction.
    pub fn directional_tension(&self, f: &[[f64; 3]; 3], direction: &[f64; 3]) -> f64 {
        let c = right_cauchy_green(f);
        let mut lambda_sq = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                lambda_sq += direction[i] * c[i][j] * direction[j];
            }
        }
        let lambda = lambda_sq.max(1e-30).sqrt();
        let mut stress = 0.0;
        for p in 0..self.mu.len() {
            let ap = self.alpha[p];
            stress += self.mu[p] * (lambda.powf(ap - 1.0) - lambda.powf(-ap / 2.0 - 1.0));
        }
        stress + self.pre_stress
    }
}
/// Growth and remodeling model for biological tissues.
///
/// Implements volumetric growth via a growth tensor Fg, so that
/// F = Fe * Fg where Fe is the elastic deformation gradient.
#[derive(Debug, Clone)]
pub struct TissueGrowth {
    /// Growth tensor (starts as identity).
    pub growth_tensor: [[f64; 3]; 3],
    /// Growth rate constant.
    pub growth_rate: f64,
    /// Maximum growth stretch.
    pub max_growth: f64,
    /// Homeostatic stress target (Pa).
    pub target_stress: f64,
}
impl TissueGrowth {
    /// Create a new tissue growth model.
    pub fn new(growth_rate: f64, max_growth: f64, target_stress: f64) -> Self {
        Self {
            growth_tensor: identity3(),
            growth_rate,
            max_growth,
            target_stress,
        }
    }
    /// Default arterial growth model.
    pub fn arterial_growth() -> Self {
        Self::new(0.01, 1.5, 100.0e3)
    }
    /// Update the growth tensor based on the current stress state.
    ///
    /// Isotropic growth: Fg = theta * I where theta is the growth stretch.
    pub fn update_isotropic(&mut self, current_stress: f64, dt: f64) {
        let stress_ratio = current_stress / self.target_stress;
        let d_theta = self.growth_rate * (stress_ratio - 1.0) * dt;
        let theta = self.growth_tensor[0][0];
        let new_theta = clamp(theta + d_theta, 1.0, self.max_growth);
        self.growth_tensor = identity3();
        for i in 0..3 {
            self.growth_tensor[i][i] = new_theta;
        }
    }
    /// Update the growth tensor for anisotropic (directional) growth.
    #[allow(clippy::needless_range_loop)]
    pub fn update_anisotropic(&mut self, current_stress: f64, direction: &[f64; 3], dt: f64) {
        let stress_ratio = current_stress / self.target_stress;
        let d_theta = self.growth_rate * (stress_ratio - 1.0) * dt;
        let nn = outer3(direction, direction);
        for i in 0..3 {
            for j in 0..3 {
                let new_val = self.growth_tensor[i][j] + d_theta * nn[i][j];
                self.growth_tensor[i][j] = if i == j {
                    clamp(new_val, 1.0, self.max_growth)
                } else {
                    new_val
                };
            }
        }
    }
    /// Compute the elastic part of the deformation gradient: Fe = F * Fg^{-1}.
    pub fn elastic_deformation(&self, f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let fg_inv = inv3(&self.growth_tensor);
        mat_mul3(f, &fg_inv)
    }
    /// Compute the current growth volume ratio det(Fg).
    pub fn growth_volume_ratio(&self) -> f64 {
        det3(&self.growth_tensor)
    }
}
/// Mooney-Rivlin model adapted for biological soft tissues.
///
/// W = c10*(I1-3) + c01*(I2-3) + c11*(I1-3)*(I2-3) + (K/2)*(J-1)^2
///
/// The c11 coupling term captures the nonlinear stiffening seen in
/// biological tissues at large strains.
#[derive(Debug, Clone)]
pub struct MooneyRivlinBiological {
    /// First material constant c10 (Pa).
    pub c10: f64,
    /// Second material constant c01 (Pa).
    pub c01: f64,
    /// Coupling constant c11 (Pa).
    pub c11: f64,
    /// Bulk modulus K (Pa).
    pub bulk_modulus: f64,
}
impl MooneyRivlinBiological {
    /// Create a new Mooney-Rivlin biological model.
    pub fn new(c10: f64, c01: f64, c11: f64, bulk_modulus: f64) -> Self {
        Self {
            c10,
            c01,
            c11,
            bulk_modulus,
        }
    }
    /// Create a model for generic soft tissue.
    pub fn soft_tissue() -> Self {
        Self::new(1.0e3, 0.5e3, 0.1e3, 50.0e3)
    }
    /// Create a model for brain tissue.
    pub fn brain_tissue() -> Self {
        Self::new(0.5e3, 0.25e3, 0.05e3, 30.0e3)
    }
    /// Compute the strain energy density.
    pub fn strain_energy(&self, f: &[[f64; 3]; 3]) -> f64 {
        let c = right_cauchy_green(f);
        let j = jacobian(f);
        let j_23 = j.powf(-2.0 / 3.0);
        let j_43 = j_23 * j_23;
        let i1_bar = j_23 * invariant_i1(&c);
        let i2_bar = j_43 * invariant_i2(&c);
        let w_dev = self.c10 * (i1_bar - 3.0)
            + self.c01 * (i2_bar - 3.0)
            + self.c11 * (i1_bar - 3.0) * (i2_bar - 3.0);
        let w_vol = 0.5 * self.bulk_modulus * (j - 1.0) * (j - 1.0);
        w_dev + w_vol
    }
    /// Compute the initial shear modulus mu = 2*(c10 + c01).
    pub fn shear_modulus(&self) -> f64 {
        2.0 * (self.c10 + self.c01)
    }
    /// Compute the second Piola-Kirchhoff stress for the deviatoric part.
    pub fn deviatoric_stress(&self, f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let c = right_cauchy_green(f);
        let j = jacobian(f);
        let j_23 = j.powf(-2.0 / 3.0);
        let j_43 = j_23 * j_23;
        let i1_bar = j_23 * invariant_i1(&c);
        let i2_bar = j_43 * invariant_i2(&c);
        let dw_di1 = self.c10 + self.c11 * (i2_bar - 3.0);
        let dw_di2 = self.c01 + self.c11 * (i1_bar - 3.0);
        let mut s = [[0.0; 3]; 3];
        let id = identity3();
        for i in 0..3 {
            for j_idx in 0..3 {
                s[i][j_idx] =
                    2.0 * (dw_di1 * id[i][j_idx] + dw_di2 * (i1_bar * id[i][j_idx] - c[i][j_idx]));
            }
        }
        s
    }
    /// Compute the effective tangent stiffness (scalar estimate).
    pub fn tangent_modulus_estimate(&self, f: &[[f64; 3]; 3]) -> f64 {
        let c = right_cauchy_green(f);
        let j = jacobian(f);
        let j_23 = j.powf(-2.0 / 3.0);
        let j_43 = j_23 * j_23;
        let i1_bar = j_23 * invariant_i1(&c);
        let i2_bar = j_43 * invariant_i2(&c);
        4.0 * (self.c10 + self.c01 * i1_bar + self.c11 * (i1_bar + i2_bar - 6.0))
    }
}
/// Tissue layer type for multi-layer models.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TissueLayerType {
    /// Epithelial tissue (skin surface).
    Epithelial,
    /// Connective tissue (dermis, fascia).
    Connective,
    /// Muscular tissue.
    Muscular,
    /// Adipose (fat) tissue.
    Adipose,
    /// Cartilage.
    Cartilage,
    /// Bone.
    Bone,
    /// Vascular tissue (blood vessels).
    Vascular,
    /// Neural tissue.
    Neural,
}
/// Hill-type muscle activation model for active force generation.
///
/// Computes the active stress contribution from muscle contraction,
/// combining force-length and force-velocity relationships with
/// neural activation dynamics.
#[derive(Debug, Clone)]
pub struct HillMuscleActivation {
    /// Maximum isometric force (N).
    pub f_max: f64,
    /// Optimal fiber length (m).
    pub l_opt: f64,
    /// Maximum contraction velocity (lengths/s).
    pub v_max: f64,
    /// Pennation angle at optimal length (rad).
    pub alpha_opt: f64,
    /// Activation time constant (s).
    pub tau_act: f64,
    /// Deactivation time constant (s).
    pub tau_deact: f64,
    /// Current activation level \[0, 1\].
    pub(super) activation: f64,
}
impl HillMuscleActivation {
    /// Create a new Hill muscle activation model.
    pub fn new(
        f_max: f64,
        l_opt: f64,
        v_max: f64,
        alpha_opt: f64,
        tau_act: f64,
        tau_deact: f64,
    ) -> Self {
        Self {
            f_max,
            l_opt,
            v_max,
            alpha_opt,
            tau_act,
            tau_deact,
            activation: 0.0,
        }
    }
    /// Create a model for a generic skeletal muscle.
    pub fn skeletal_muscle(f_max: f64, l_opt: f64) -> Self {
        Self::new(f_max, l_opt, 10.0 * l_opt, 0.0, 0.01, 0.04)
    }
    /// Set the activation level directly.
    pub fn set_activation(&mut self, a: f64) {
        self.activation = clamp(a, 0.0, 1.0);
    }
    /// Get the current activation level.
    pub fn activation(&self) -> f64 {
        self.activation
    }
    /// Update the activation dynamics given a neural excitation u in \[0,1\].
    pub fn update_activation(&mut self, u: f64, dt: f64) {
        let u_clamped = clamp(u, 0.0, 1.0);
        let tau = if u_clamped > self.activation {
            self.tau_act
        } else {
            self.tau_deact
        };
        let da = (u_clamped - self.activation) / tau * dt;
        self.activation = clamp(self.activation + da, 0.0, 1.0);
    }
    /// Normalized force-length relationship (Gaussian).
    ///
    /// f_L(l_norm) = exp(-((l_norm - 1) / gamma)^2)
    pub fn force_length(&self, fiber_length: f64) -> f64 {
        let l_norm = fiber_length / self.l_opt;
        let gamma = 0.45;
        (-(((l_norm - 1.0) / gamma).powi(2))).exp()
    }
    /// Normalized force-velocity relationship (Hill hyperbola).
    ///
    /// Concentric: f_V = (v_max - v) / (v_max + v/a_hill)
    /// Eccentric: f_V = (1.8 - 0.8 * (v_max + v) / (v_max - 7.56/a_hill * v))
    pub fn force_velocity(&self, velocity: f64) -> f64 {
        let v_norm = velocity / self.v_max;
        let a_hill = 0.25;
        if v_norm <= 0.0 {
            let v_abs = (-v_norm).min(0.99);
            (1.0 - v_abs) / (1.0 + v_abs / a_hill)
        } else {
            let f_ecc_max = 1.8;
            let v_clamped = v_norm.min(0.99);
            1.0 + (f_ecc_max - 1.0) * v_clamped / (v_clamped + 0.1)
        }
    }
    /// Passive force-length relationship (exponential).
    pub fn passive_force_length(&self, fiber_length: f64) -> f64 {
        let l_norm = fiber_length / self.l_opt;
        if l_norm > 1.0 {
            let strain = l_norm - 1.0;
            0.05 * ((10.0 * strain).exp() - 1.0)
        } else {
            0.0
        }
    }
    /// Compute the total muscle force.
    pub fn compute_force(&self, fiber_length: f64, velocity: f64) -> f64 {
        let f_l = self.force_length(fiber_length);
        let f_v = self.force_velocity(velocity);
        let f_pe = self.passive_force_length(fiber_length);
        let cos_alpha = self.alpha_opt.cos();
        self.f_max * (self.activation * f_l * f_v + f_pe) * cos_alpha
    }
    /// Compute the active stress (force per unit area) for a given
    /// cross-sectional area.
    pub fn active_stress(&self, fiber_length: f64, velocity: f64, area: f64) -> f64 {
        if area > 1e-15 {
            self.compute_force(fiber_length, velocity) / area
        } else {
            0.0
        }
    }
}
/// Fung hyperelastic model for soft biological tissues.
///
/// The strain energy density is:
///   W = (c/2) * (exp(Q) - 1)
/// where Q = E : D : E (a quadratic form of the Green-Lagrange strain E
/// weighted by the material tensor D).
///
/// This model captures the exponential stiffening observed in arteries,
/// ligaments, and other collagenous tissues.
#[derive(Debug, Clone)]
pub struct FungHyperelastic {
    /// Material constant c (Pa).
    pub c: f64,
    /// Diagonal entries of the material tensor D (Voigt: D11, D22, D33, D12, D13, D23).
    pub d: [f64; 6],
}
impl FungHyperelastic {
    /// Create a new Fung hyperelastic model.
    pub fn new(c: f64, d: [f64; 6]) -> Self {
        Self { c, d }
    }
    /// Create with isotropic D (same stiffness in all directions).
    pub fn isotropic(c: f64, d_val: f64) -> Self {
        Self::new(c, [d_val, d_val, d_val, 0.0, 0.0, 0.0])
    }
    /// Compute the quadratic form Q = E : D : E.
    ///
    /// Using Voigt notation:
    ///   Q = D11*E11^2 + D22*E22^2 + D33*E33^2
    ///     + 2*D12*E11*E22 + 2*D13*E11*E33 + 2*D23*E22*E33
    ///     (off-diagonal shear terms neglected in the simplified form)
    pub fn compute_q(&self, e: &[[f64; 3]; 3]) -> f64 {
        let e11 = e[0][0];
        let e22 = e[1][1];
        let e33 = e[2][2];
        let e12 = e[0][1];
        let e13 = e[0][2];
        let e23 = e[1][2];
        self.d[0] * e11 * e11
            + self.d[1] * e22 * e22
            + self.d[2] * e33 * e33
            + 2.0 * self.d[3] * e11 * e22
            + 2.0 * self.d[4] * e11 * e33
            + 2.0 * self.d[5] * e22 * e33
            + 4.0 * self.d[0].min(self.d[1]) * e12 * e12
            + 4.0 * self.d[0].min(self.d[2]) * e13 * e13
            + 4.0 * self.d[1].min(self.d[2]) * e23 * e23
    }
    /// Compute the strain energy density W(F).
    pub fn strain_energy(&self, f: &[[f64; 3]; 3]) -> f64 {
        let e = green_lagrange(f);
        let q = self.compute_q(&e);
        0.5 * self.c * (q.exp() - 1.0)
    }
    /// Compute the second Piola-Kirchhoff stress S = dW/dE.
    ///
    /// S = c * exp(Q) * (dQ/dE)
    pub fn second_piola_kirchhoff(&self, f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let e = green_lagrange(f);
        let q = self.compute_q(&e);
        let exp_q = q.exp();
        let e11 = e[0][0];
        let e22 = e[1][1];
        let e33 = e[2][2];
        let e12 = e[0][1];
        let e13 = e[0][2];
        let e23 = e[1][2];
        let dq_de11 = 2.0 * self.d[0] * e11 + 2.0 * self.d[3] * e22 + 2.0 * self.d[4] * e33;
        let dq_de22 = 2.0 * self.d[1] * e22 + 2.0 * self.d[3] * e11 + 2.0 * self.d[5] * e33;
        let dq_de33 = 2.0 * self.d[2] * e33 + 2.0 * self.d[4] * e11 + 2.0 * self.d[5] * e22;
        let dq_de12 = 4.0 * self.d[0].min(self.d[1]) * e12;
        let dq_de13 = 4.0 * self.d[0].min(self.d[2]) * e13;
        let dq_de23 = 4.0 * self.d[1].min(self.d[2]) * e23;
        let factor = self.c * exp_q;
        [
            [factor * dq_de11, factor * dq_de12, factor * dq_de13],
            [factor * dq_de12, factor * dq_de22, factor * dq_de23],
            [factor * dq_de13, factor * dq_de23, factor * dq_de33],
        ]
    }
    /// Compute the tangent stiffness at the current deformation.
    pub fn tangent_modulus(&self, f: &[[f64; 3]; 3]) -> f64 {
        let e = green_lagrange(f);
        let q = self.compute_q(&e);
        self.c * q.exp() * (1.0 + q)
    }
}
