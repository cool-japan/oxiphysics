//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Holzapfel-Gasser-Ogden (HGO) model for fiber-reinforced biological tissue.
///
/// Strain energy:
///
///   W = c₁/2 * (Ī₁ - 3) + k₁/(2k₂) * Σᵢ (exp(k₂*(κ*Ī₁ + (1-3κ)*Ī₄ᵢ - 1)²) - 1)
///       + κ_bulk/2 * (J-1)²
///
/// where κ ∈ \[0, 1/3\] is the fiber dispersion parameter:
/// - κ = 0: aligned fibers (no dispersion)
/// - κ = 1/3: fully isotropic distribution
///
/// Two fiber families with mean directions a₁ and a₂.
#[derive(Debug, Clone)]
pub struct HolzapfelGasserOgden {
    /// Neo-Hookean ground matrix stiffness c₁ (Pa).
    pub c1: f64,
    /// Fiber stiffness k₁ (Pa).
    pub k1: f64,
    /// Fiber nonlinearity k₂ (dimensionless).
    pub k2: f64,
    /// Fiber dispersion κ ∈ \[0, 1/3\].
    pub kappa_disp: f64,
    /// Fiber family 1 mean direction (unit vector).
    pub a1: [f64; 3],
    /// Fiber family 2 mean direction (unit vector).
    pub a2: [f64; 3],
    /// Bulk modulus (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl HolzapfelGasserOgden {
    /// Create an HGO model.
    pub fn new(
        c1: f64,
        k1: f64,
        k2: f64,
        kappa_disp: f64,
        a1: [f64; 3],
        a2: [f64; 3],
        bulk_modulus: f64,
    ) -> Self {
        Self {
            c1,
            k1,
            k2,
            kappa_disp,
            a1,
            a2,
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
            1.0e6,
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
            1.0e6,
        )
    }
    /// Generalized fiber invariant: I_4_eff = κ*I₁_bar + (1 - 3κ)*I₄_bar.
    fn generalized_invariant(i1_bar: f64, i4_bar: f64, kappa_disp: f64) -> f64 {
        kappa_disp * i1_bar + (1.0 - 3.0 * kappa_disp) * i4_bar
    }
    /// Isochoric fiber stretch squared: Ī₄ = a · C̄ · a.
    fn isochoric_fiber_stretch_sq(f: &[[f64; 3]; 3], a: &[f64; 3]) -> f64 {
        let j = det3(f);
        let j_13 = j.powf(-1.0 / 3.0);
        let f_bar: [[f64; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|jj| j_13 * f[i][jj]));
        let fa: [f64; 3] =
            std::array::from_fn(|i| f_bar[i].iter().zip(a.iter()).map(|(fb, ak)| fb * ak).sum());
        fa[0] * fa[0] + fa[1] * fa[1] + fa[2] * fa[2]
    }
    /// Strain energy density from deformation gradient F.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let j = det3(f);
        let (i1, _i2, i3) = deformation_invariants(f);
        let j_23 = i3.sqrt().powf(-2.0 / 3.0);
        let i1_bar = j_23 * i1;
        let i4_bar_1 = Self::isochoric_fiber_stretch_sq(f, &self.a1);
        let i4_bar_2 = Self::isochoric_fiber_stretch_sq(f, &self.a2);
        let e1 = Self::generalized_invariant(i1_bar, i4_bar_1, self.kappa_disp) - 1.0;
        let e2 = Self::generalized_invariant(i1_bar, i4_bar_2, self.kappa_disp) - 1.0;
        let w_matrix = self.c1 / 2.0 * (i1_bar - 3.0);
        let w_f1 = if e1 > 0.0 {
            self.k1 / (2.0 * self.k2) * ((self.k2 * e1 * e1).exp() - 1.0)
        } else {
            0.0
        };
        let w_f2 = if e2 > 0.0 {
            self.k1 / (2.0 * self.k2) * ((self.k2 * e2 * e2).exp() - 1.0)
        } else {
            0.0
        };
        let w_vol = self.bulk_modulus / 2.0 * (j - 1.0).powi(2);
        w_matrix + w_f1 + w_f2 + w_vol
    }
    /// Check if fibers are active (strain energy > 0 from fiber families).
    pub fn fibers_active(&self, f: &[[f64; 3]; 3]) -> (bool, bool) {
        let (i1, _i2, i3) = deformation_invariants(f);
        let j_23 = i3.sqrt().powf(-2.0 / 3.0);
        let i1_bar = j_23 * i1;
        let i4_bar_1 = Self::isochoric_fiber_stretch_sq(f, &self.a1);
        let i4_bar_2 = Self::isochoric_fiber_stretch_sq(f, &self.a2);
        let e1 = Self::generalized_invariant(i1_bar, i4_bar_1, self.kappa_disp) - 1.0;
        let e2 = Self::generalized_invariant(i1_bar, i4_bar_2, self.kappa_disp) - 1.0;
        (e1 > 0.0, e2 > 0.0)
    }
}
/// Gent hyperelastic model.
///
/// W = -μ*J_m/2 * ln(1 - (I1 - 3)/J_m) + D1*(J-1)^2
///
/// where J_m is the limiting value of I1-3 (locking stretch parameter).
#[derive(Debug, Clone)]
pub struct Gent {
    /// Shear modulus μ (Pa).
    pub mu: f64,
    /// Limiting parameter J_m (dimensionless, > 0).
    pub j_m: f64,
    /// Bulk modulus (Pa).
    pub bulk_modulus: f64,
}
impl Gent {
    /// Create a new Gent material.
    pub fn new(mu: f64, j_m: f64, bulk_modulus: f64) -> Self {
        Self {
            mu,
            j_m,
            bulk_modulus,
        }
    }
    /// Strain energy density from principal stretches.
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let j13 = j.powf(-1.0 / 3.0);
        let l1 = j13 * lambda1;
        let l2 = j13 * lambda2;
        let l3 = j13 * lambda3;
        let i1_bar = l1 * l1 + l2 * l2 + l3 * l3;
        let arg = 1.0 - (i1_bar - 3.0) / self.j_m;
        if arg <= 0.0 {
            return f64::INFINITY;
        }
        let d1 = self.bulk_modulus / 2.0;
        -self.mu * self.j_m / 2.0 * arg.ln() + d1 * (j - 1.0).powi(2)
    }
    /// Strain energy density from deformation gradient.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
    /// Small-strain shear modulus.
    pub fn shear_modulus(&self) -> f64 {
        self.mu
    }
    /// Initial bulk modulus.
    pub fn bulk_modulus(&self) -> f64 {
        self.bulk_modulus
    }
}
/// Drucker-Prager plasticity model for geomaterials (soil, rock, concrete).
///
/// Yield surface: F = sqrt(J2) + alpha * I1 - k = 0
/// where alpha and k are derived from friction angle and cohesion.
#[derive(Debug, Clone)]
pub struct DruckerPrager {
    /// Friction angle (radians).
    pub friction_angle: f64,
    /// Cohesion (Pa).
    pub cohesion: f64,
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}
impl DruckerPrager {
    /// Create a new Drucker-Prager model.
    pub fn new(friction_angle: f64, cohesion: f64, young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            friction_angle,
            cohesion,
            young_modulus,
            poisson_ratio,
        }
    }
    /// Alpha parameter from friction angle (outer cone match to Mohr-Coulomb).
    pub fn alpha(&self) -> f64 {
        2.0 * self.friction_angle.sin() / (3.0_f64.sqrt() * (3.0 - self.friction_angle.sin()))
    }
    /// k parameter from cohesion and friction angle.
    pub fn k(&self) -> f64 {
        6.0 * self.cohesion * self.friction_angle.cos()
            / (3.0_f64.sqrt() * (3.0 - self.friction_angle.sin()))
    }
    /// Evaluate yield function: F > 0 means yielding.
    pub fn yield_function(&self, stress: &[f64; 6]) -> f64 {
        let i1 = stress[0] + stress[1] + stress[2];
        let j2 = Self::j2_invariant(stress);
        j2.sqrt() + self.alpha() * i1 - self.k()
    }
    /// Second deviatoric stress invariant J2.
    fn j2_invariant(stress: &[f64; 6]) -> f64 {
        let mean = (stress[0] + stress[1] + stress[2]) / 3.0;
        let dev = [stress[0] - mean, stress[1] - mean, stress[2] - mean];
        0.5 * (dev[0].powi(2) + dev[1].powi(2) + dev[2].powi(2))
            + stress[3].powi(2)
            + stress[4].powi(2)
            + stress[5].powi(2)
    }
}
/// JWL (Jones-Wilkins-Lee) equation of state for detonation products.
///
/// P = A*(1 - omega/(R1*V))*exp(-R1*V) + B*(1 - omega/(R2*V))*exp(-R2*V) + omega*E/V
///
/// Where V = rho0/rho (relative volume).
#[derive(Debug, Clone)]
pub struct JwlEos {
    /// First pressure coefficient A (Pa).
    pub a: f64,
    /// Second pressure coefficient B (Pa).
    pub b: f64,
    /// First exponential coefficient R1.
    pub r1: f64,
    /// Second exponential coefficient R2.
    pub r2: f64,
    /// Gruneisen parameter omega.
    pub omega: f64,
    /// Reference density (kg/m^3).
    pub rho0: f64,
}
impl JwlEos {
    /// Create a new JWL EOS.
    pub fn new(a: f64, b: f64, r1: f64, r2: f64, omega: f64, rho0: f64) -> Self {
        Self {
            a,
            b,
            r1,
            r2,
            omega,
            rho0,
        }
    }
    /// Compute pressure from density and specific internal energy.
    pub fn pressure(&self, density: f64, energy: f64) -> f64 {
        let v = self.rho0 / density;
        let p1 = self.a * (1.0 - self.omega / (self.r1 * v)) * (-self.r1 * v).exp();
        let p2 = self.b * (1.0 - self.omega / (self.r2 * v)) * (-self.r2 * v).exp();
        let p3 = self.omega * energy * density;
        p1 + p2 + p3
    }
}
/// Yeoh hyperelastic model (3-parameter).
///
/// W = C10*(I1-3) + C20*(I1-3)^2 + C30*(I1-3)^3 + D1*(J-1)^2
#[derive(Debug, Clone)]
pub struct Yeoh {
    /// First material parameter C10 (Pa).
    pub c10: f64,
    /// Second material parameter C20 (Pa).
    pub c20: f64,
    /// Third material parameter C30 (Pa).
    pub c30: f64,
    /// Bulk modulus (Pa).
    pub bulk_modulus: f64,
}
impl Yeoh {
    /// Create a new Yeoh material.
    pub fn new(c10: f64, c20: f64, c30: f64, bulk_modulus: f64) -> Self {
        Self {
            c10,
            c20,
            c30,
            bulk_modulus,
        }
    }
    /// Initial shear modulus: μ = 2*C10.
    pub fn shear_modulus(&self) -> f64 {
        2.0 * self.c10
    }
    /// Strain energy density from principal stretches.
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let j13 = j.powf(-1.0 / 3.0);
        let l1 = j13 * lambda1;
        let l2 = j13 * lambda2;
        let l3 = j13 * lambda3;
        let i1_bar = l1 * l1 + l2 * l2 + l3 * l3;
        let di = i1_bar - 3.0;
        let d1 = self.bulk_modulus / 2.0;
        self.c10 * di + self.c20 * di.powi(2) + self.c30 * di.powi(3) + d1 * (j - 1.0).powi(2)
    }
    /// Strain energy density from deformation gradient.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
}
/// Mooney-Rivlin hyperelastic material model.
///
/// Strain energy: W = C10 * (I1 - 3) + C01 * (I2 - 3) + D1 * (J - 1)^2
///
/// Where I1, I2 are reduced invariants of the right Cauchy-Green tensor
/// and D1 = K/2 enforces near-incompressibility.
#[derive(Debug, Clone)]
pub struct MooneyRivlin {
    /// First material constant (Pa).
    pub c10: f64,
    /// Second material constant (Pa).
    pub c01: f64,
    /// Bulk modulus (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl MooneyRivlin {
    /// Create a new Mooney-Rivlin material.
    pub fn new(c10: f64, c01: f64, bulk_modulus: f64) -> Self {
        Self {
            c10,
            c01,
            bulk_modulus,
        }
    }
    /// Create from shear modulus and ratio.
    /// mu = 2 * (c10 + c01), ratio = c10 / c01.
    pub fn from_shear_and_ratio(shear_modulus: f64, ratio: f64, bulk_modulus: f64) -> Self {
        let c01 = shear_modulus / (2.0 * (1.0 + ratio));
        let c10 = ratio * c01;
        Self::new(c10, c01, bulk_modulus)
    }
    /// Initial shear modulus mu = 2 * (C10 + C01).
    pub fn shear_modulus(&self) -> f64 {
        2.0 * (self.c10 + self.c01)
    }
    /// Strain energy density.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (i1, i2, i3) = deformation_invariants(f);
        let j = i3.sqrt();
        let j_23 = j.powf(-2.0 / 3.0);
        let i1_bar = j_23 * i1;
        let i2_bar = j_23 * j_23 * i2;
        let d1 = self.bulk_modulus / 2.0;
        self.c10 * (i1_bar - 3.0) + self.c01 * (i2_bar - 3.0) + d1 * (j - 1.0).powi(2)
    }
    /// First Piola-Kirchhoff stress tensor P = dW/dF.
    pub fn first_piola_kirchhoff_stress(&self, f: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let j = det3(f);
        if j.abs() < 1e-30 {
            return [[0.0; 3]; 3];
        }
        let j_23 = j.powf(-2.0 / 3.0);
        let f_inv_t = transpose3(&inv3(f));
        let c: [[f64; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|k| (0..3).map(|l| f[l][i] * f[l][k]).sum())
        });
        let i1 = c[0][0] + c[1][1] + c[2][2];
        let c2_trace: f64 = (0..3)
            .flat_map(|a| (0..3).map(move |b| c[a][b] * c[b][a]))
            .sum();
        let mut p = [[0.0; 3]; 3];
        for (i, pi) in p.iter_mut().enumerate() {
            for (j_idx, pij) in pi.iter_mut().enumerate() {
                let iso1 = 2.0 * self.c10 * j_23 * (f[i][j_idx] - i1 / 3.0 * f_inv_t[i][j_idx]);
                let cfij: f64 = (0..3).map(|k| c[j_idx][k] * f[i][k]).sum();
                let iso2 = 2.0
                    * self.c01
                    * j_23
                    * j_23
                    * (i1 * f[i][j_idx] - cfij - (i1 * i1 - c2_trace) / 3.0 * f_inv_t[i][j_idx]);
                let vol = self.bulk_modulus * (j - 1.0) * j * f_inv_t[i][j_idx];
                *pij = iso1 + iso2 + vol;
            }
        }
        p
    }
}
/// Bilinear traction-separation cohesive zone law.
///
/// Used for modelling delamination in composites and adhesive joints.
///
/// T(δ) = {  K δ                          for δ ≤ δ₀  (linear loading)
///         { T_max (δ_f - δ)/(δ_f - δ₀)  for δ₀ < δ ≤ δ_f  (softening)
///         { 0                            for δ > δ_f  (fully failed)
///
/// where δ₀ is the damage initiation displacement, δ_f is the complete
/// failure displacement, K is the interface stiffness, T_max = K δ₀.
#[derive(Debug, Clone)]
pub struct BilinearCohesiveZone {
    /// Interface stiffness K (Pa/m).
    pub stiffness: f64,
    /// Damage initiation displacement δ₀ (m).
    pub delta0: f64,
    /// Complete failure displacement δ_f (m).
    pub delta_f: f64,
}
impl BilinearCohesiveZone {
    /// Create a new bilinear cohesive zone law.
    pub fn new(stiffness: f64, delta0: f64, delta_f: f64) -> Self {
        assert!(delta_f > delta0, "delta_f must exceed delta0");
        Self {
            stiffness,
            delta0,
            delta_f,
        }
    }
    /// Peak traction T_max = K · δ₀ (Pa).
    pub fn peak_traction(&self) -> f64 {
        self.stiffness * self.delta0
    }
    /// Mode-I fracture toughness G_c = T_max · δ_f / 2 (J/m²).
    pub fn fracture_toughness(&self) -> f64 {
        self.peak_traction() * self.delta_f / 2.0
    }
    /// Traction T(δ) for opening displacement δ ≥ 0.
    pub fn traction(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        if delta <= self.delta0 {
            self.stiffness * delta
        } else if delta < self.delta_f {
            self.peak_traction() * (self.delta_f - delta) / (self.delta_f - self.delta0)
        } else {
            0.0
        }
    }
    /// Damage variable D ∈ \[0, 1\] given current opening δ and maximum
    /// opening δ_max experienced so far.
    pub fn damage_variable(&self, delta_max: f64) -> f64 {
        if delta_max <= self.delta0 {
            0.0
        } else if delta_max >= self.delta_f {
            1.0
        } else {
            self.delta_f * (delta_max - self.delta0) / (delta_max * (self.delta_f - self.delta0))
        }
    }
    /// Energy dissipated up to displacement δ (J/m²).
    pub fn dissipated_energy(&self, delta: f64) -> f64 {
        if delta <= 0.0 {
            0.0
        } else if delta <= self.delta0 {
            0.5 * self.stiffness * delta * delta
        } else if delta < self.delta_f {
            let t_max = self.peak_traction();
            let e0 = 0.5 * t_max * self.delta0;
            let slope = t_max / (self.delta_f - self.delta0);
            let dd = delta - self.delta0;
            let ddf = self.delta_f - self.delta0;
            e0 + t_max * dd - 0.5 * slope * dd * dd - (t_max * ddf - 0.5 * slope * ddf * ddf) * 0.0
        } else {
            self.fracture_toughness()
        }
    }
}
/// J2 (von Mises) plasticity model.
///
/// Combines linear elasticity with yield surface based on deviatoric stress.
/// Uses isotropic hardening.
#[derive(Debug, Clone)]
pub struct J2Plasticity {
    /// Young's modulus (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Initial yield stress (Pa).
    pub yield_stress: f64,
    /// Hardening modulus H (Pa).
    pub hardening_modulus: f64,
}
impl J2Plasticity {
    /// Create a new J2 plasticity model.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        yield_stress: f64,
        hardening_modulus: f64,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            yield_stress,
            hardening_modulus,
        }
    }
    /// Shear modulus G = E / (2(1+nu)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
    /// Bulk modulus K = E / (3(1-2nu)).
    pub fn bulk_modulus(&self) -> f64 {
        self.young_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
    /// Current yield stress with isotropic hardening.
    pub fn current_yield_stress(&self, equivalent_plastic_strain: f64) -> f64 {
        self.yield_stress + self.hardening_modulus * equivalent_plastic_strain
    }
    /// Von Mises stress from Cauchy stress tensor (Voigt: \[xx, yy, zz, xy, yz, xz\]).
    pub fn von_mises_stress(stress: &[f64; 6]) -> f64 {
        let s_xx = stress[0];
        let s_yy = stress[1];
        let s_zz = stress[2];
        let s_xy = stress[3];
        let s_yz = stress[4];
        let s_xz = stress[5];
        let term1 = (s_xx - s_yy).powi(2) + (s_yy - s_zz).powi(2) + (s_zz - s_xx).powi(2);
        let term2 = 6.0 * (s_xy.powi(2) + s_yz.powi(2) + s_xz.powi(2));
        ((term1 + term2) / 2.0).sqrt()
    }
    /// Radial return mapping for stress update.
    /// Returns (updated_stress, delta_plastic_strain).
    pub fn return_mapping(
        &self,
        trial_stress: &[f64; 6],
        equivalent_plastic_strain: f64,
    ) -> ([f64; 6], f64) {
        let vm = Self::von_mises_stress(trial_stress);
        let sigma_y = self.current_yield_stress(equivalent_plastic_strain);
        if vm <= sigma_y {
            return (*trial_stress, 0.0);
        }
        let g = self.shear_modulus();
        let delta_gamma = (vm - sigma_y) / (3.0 * g + self.hardening_modulus);
        let mean = (trial_stress[0] + trial_stress[1] + trial_stress[2]) / 3.0;
        let dev = [
            trial_stress[0] - mean,
            trial_stress[1] - mean,
            trial_stress[2] - mean,
            trial_stress[3],
            trial_stress[4],
            trial_stress[5],
        ];
        let scale = 1.0 - 3.0 * g * delta_gamma / vm;
        let stress = [
            mean + dev[0] * scale,
            mean + dev[1] * scale,
            mean + dev[2] * scale,
            dev[3] * scale,
            dev[4] * scale,
            dev[5] * scale,
        ];
        (stress, delta_gamma)
    }
}
/// Consistent algorithmic tangent modulus for the J2 plasticity return mapping.
///
/// Given the result of a plastic return, compute the consistent tangent
/// operator (scalar) for the uniaxial case.
///
/// C_ep = E · H / (E + H)   (isotropic hardening)
impl J2Plasticity {
    /// Consistent tangent modulus C_ep = E · H / (E + H).
    ///
    /// Useful for predictor-corrector integration.
    pub fn consistent_tangent_uniaxial(&self) -> f64 {
        self.young_modulus * self.hardening_modulus / (self.young_modulus + self.hardening_modulus)
    }
    /// Isotropic hardening: compute flow stress σ_y at equivalent plastic strain.
    pub fn flow_stress(&self, ep_bar: f64) -> f64 {
        self.current_yield_stress(ep_bar)
    }
    /// Kinematic hardening stress update (back-stress evolution).
    ///
    /// Simple Prager rule: Δα = (2/3) C Δε_p  where C is the kinematic hardening modulus.
    pub fn back_stress_increment(&self, c_kinematic: f64, delta_ep_voigt: &[f64; 6]) -> [f64; 6] {
        let scale = 2.0 / 3.0 * c_kinematic;
        [
            scale * delta_ep_voigt[0],
            scale * delta_ep_voigt[1],
            scale * delta_ep_voigt[2],
            scale * delta_ep_voigt[3],
            scale * delta_ep_voigt[4],
            scale * delta_ep_voigt[5],
        ]
    }
    /// Equivalent plastic strain increment Δε̄_p from the return-mapping result.
    ///
    /// Δε̄_p = sqrt(2/3) |Δε_p|  (von Mises equivalent).
    pub fn equivalent_plastic_strain_increment(delta_ep_voigt: &[f64; 6]) -> f64 {
        let sum = delta_ep_voigt[0].powi(2)
            + delta_ep_voigt[1].powi(2)
            + delta_ep_voigt[2].powi(2)
            + 2.0
                * (delta_ep_voigt[3].powi(2)
                    + delta_ep_voigt[4].powi(2)
                    + delta_ep_voigt[5].powi(2));
        (2.0 / 3.0 * sum).sqrt()
    }
    /// Perform a full plane-stress return mapping (σ_zz = 0 enforced iteratively).
    ///
    /// Returns updated (σ_xx, σ_yy, σ_xy) and Δε̄_p.
    pub fn return_mapping_plane_stress(
        &self,
        trial_sxx: f64,
        trial_syy: f64,
        trial_sxy: f64,
        ep_bar: f64,
    ) -> ([f64; 3], f64) {
        let trial_stress_full = [trial_sxx, trial_syy, 0.0, trial_sxy, 0.0, 0.0];
        let vm = Self::von_mises_stress(&trial_stress_full);
        let sigma_y = self.current_yield_stress(ep_bar);
        if vm <= sigma_y {
            return ([trial_sxx, trial_syy, trial_sxy], 0.0);
        }
        let g = self.shear_modulus();
        let delta_gamma = (vm - sigma_y) / (3.0 * g + self.hardening_modulus);
        let scale = 1.0 - 3.0 * g * delta_gamma / vm;
        let mean = (trial_sxx + trial_syy) / 3.0;
        let dev_xx = trial_sxx - mean;
        let dev_yy = trial_syy - mean;
        (
            [
                mean + dev_xx * scale,
                mean + dev_yy * scale,
                trial_sxy * scale,
            ],
            delta_gamma,
        )
    }
}
/// Arruda-Boyce (8-chain) hyperelastic model for rubber.
///
/// Strain energy: W = μ * Σ_{i=1}^{5} C_i / λ_m^{2(i-1)} * (I1^i - 3^i)
/// where C_1=0.5, C_2=1/20, C_3=11/1050, C_4=19/7000, C_5=519/673750.
#[derive(Debug, Clone)]
pub struct ArrudaBoyce {
    /// Initial shear modulus μ (Pa).
    pub mu: f64,
    /// Locking stretch λ_m (≥ 1).
    pub lambda_m: f64,
    /// Bulk modulus (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl ArrudaBoyce {
    /// Arruda-Boyce series coefficients C_i.
    const C: [f64; 5] = [
        0.5,
        1.0 / 20.0,
        11.0 / 1050.0,
        19.0 / 7000.0,
        519.0 / 673750.0,
    ];
    /// Create a new Arruda-Boyce model.
    pub fn new(mu: f64, lambda_m: f64, bulk_modulus: f64) -> Self {
        Self {
            mu,
            lambda_m,
            bulk_modulus,
        }
    }
    /// Strain energy density from principal stretches.
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let j13 = j.powf(-1.0 / 3.0);
        let l1 = j13 * lambda1;
        let l2 = j13 * lambda2;
        let l3 = j13 * lambda3;
        let i1_bar = l1 * l1 + l2 * l2 + l3 * l3;
        let mut w = 0.0;
        let lm2 = self.lambda_m * self.lambda_m;
        let mut lm_pow = 1.0_f64;
        for i in 0..5 {
            let i1_pow_i = i1_bar.powi((i + 1) as i32);
            let three_pow_i = 3.0_f64.powi((i + 1) as i32);
            w += self.mu * Self::C[i] / lm_pow * (i1_pow_i - three_pow_i);
            lm_pow *= lm2;
        }
        let d1 = self.bulk_modulus / 2.0;
        w + d1 * (j - 1.0).powi(2)
    }
    /// Strain energy density from deformation gradient F.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
    /// Initial shear modulus.
    pub fn shear_modulus(&self) -> f64 {
        self.mu
    }
}
/// Mie-Gruneisen equation of state for shock-loaded solids.
///
/// P = P_H + gamma * rho * (E - E_H)
///
/// Where P_H is the Hugoniot pressure and E_H is the Hugoniot energy.
#[derive(Debug, Clone)]
pub struct MieGruneisenEos {
    /// Reference density (kg/m^3).
    pub rho0: f64,
    /// Reference sound speed (m/s).
    pub c0: f64,
    /// Slope of Us-Up Hugoniot (dimensionless).
    pub s: f64,
    /// Gruneisen parameter.
    pub gamma0: f64,
}
impl MieGruneisenEos {
    /// Create a new Mie-Gruneisen EOS.
    pub fn new(rho0: f64, c0: f64, s: f64, gamma0: f64) -> Self {
        Self {
            rho0,
            c0,
            s,
            gamma0,
        }
    }
    /// Compute pressure from density and specific internal energy.
    pub fn pressure(&self, density: f64, energy: f64) -> f64 {
        let eta = 1.0 - self.rho0 / density;
        if eta.abs() < 1e-12 {
            return self.gamma0 * density * energy;
        }
        let denom = 1.0 - self.s * eta;
        if denom.abs() < 1e-12 {
            return self.gamma0 * density * energy;
        }
        let p_h = self.rho0 * self.c0 * self.c0 * eta / (denom * denom);
        let e_h = 0.5 * p_h * eta / self.rho0;
        p_h + self.gamma0 * density * (energy - e_h)
    }
    /// Bulk sound speed at reference state.
    pub fn sound_speed_reference(&self) -> f64 {
        self.c0
    }
}
/// Fung hyperelastic model for soft biological tissue.
///
/// Exponential strain energy density:
///
///   W = c/2 * (exp(Q) - 1)
///
/// where Q = b1*(I1-3) + b2*(I1-3)^2 for isotropic tissue,
/// or a more general anisotropic form including invariants I4.
///
/// The exponential form captures the toe-region (recruitment) and
/// linear stiffening of collagen-reinforced tissues.
#[derive(Debug, Clone)]
pub struct Fung {
    /// Overall stiffness parameter c (Pa).
    pub c: f64,
    /// Linear invariant coefficient b1.
    pub b1: f64,
    /// Quadratic invariant coefficient b2.
    pub b2: f64,
    /// Bulk modulus (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl Fung {
    /// Create a new Fung model.
    pub fn new(c: f64, b1: f64, b2: f64, bulk_modulus: f64) -> Self {
        Self {
            c,
            b1,
            b2,
            bulk_modulus,
        }
    }
    /// Aorta (approximate Fung parameters, Demiray 1972).
    ///
    /// c = 26.95 Pa, b1 = 1.95, b2 = 0.0.
    pub fn aorta() -> Self {
        Self::new(26.95, 1.95, 0.0, 1.0e6)
    }
    /// Myocardium passive (approximate).
    pub fn myocardium() -> Self {
        Self::new(100.0, 4.0, 0.5, 1.0e6)
    }
    /// Strain energy density from principal stretches.
    ///
    /// Uses isochoric I1 (reduced first invariant).
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let j13 = j.powf(-1.0 / 3.0);
        let l1 = j13 * lambda1;
        let l2 = j13 * lambda2;
        let l3 = j13 * lambda3;
        let i1_bar = l1 * l1 + l2 * l2 + l3 * l3;
        let di = i1_bar - 3.0;
        let q = self.b1 * di + self.b2 * di * di;
        let d1 = self.bulk_modulus / 2.0;
        self.c / 2.0 * (q.exp() - 1.0) + d1 * (j - 1.0).powi(2)
    }
    /// Strain energy density from deformation gradient F.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
    /// Cauchy stress estimate for uniaxial incompressible loading via FD.
    ///
    /// σ = λ * dW/dλ (incompressible uniaxial, using finite difference).
    pub fn uniaxial_cauchy_stress(&self, lambda: f64) -> f64 {
        let inv = 1.0 / lambda.sqrt();
        let h = 1e-6;
        let w_plus = self.strain_energy_from_stretches(lambda + h, inv, inv);
        let w_minus = self.strain_energy_from_stretches(lambda - h, inv, inv);
        lambda * (w_plus - w_minus) / (2.0 * h)
    }
    /// Small-strain shear modulus: μ = c * b1.
    pub fn small_strain_shear_modulus(&self) -> f64 {
        self.c * self.b1
    }
}
/// Hencky hyperelastic model (compressible isotropic).
///
/// Strain energy in terms of principal logarithmic (Hencky) strains h_i = ln λ_i:
///
///   W = λ/2 (h₁+h₂+h₃)² + μ (h₁²+h₂²+h₃²)
///
/// where λ = K - 2G/3 (Lamé first parameter) and μ = G (shear modulus).
/// This form recovers St.Venant-Kirchhoff in the small-strain limit.
#[derive(Debug, Clone)]
pub struct Hencky {
    /// Lamé first parameter λ (Pa).
    pub lame_lambda: f64,
    /// Shear modulus μ (Pa).
    pub shear_modulus: f64,
}
impl Hencky {
    /// Create a new Hencky model from Lamé parameters.
    pub fn new(lame_lambda: f64, shear_modulus: f64) -> Self {
        Self {
            lame_lambda,
            shear_modulus,
        }
    }
    /// Create from Young's modulus and Poisson's ratio.
    pub fn from_young_poisson(e: f64, nu: f64) -> Self {
        let lame_lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        Self::new(lame_lambda, mu)
    }
    /// Bulk modulus K = λ + 2G/3.
    pub fn bulk_modulus(&self) -> f64 {
        self.lame_lambda + 2.0 * self.shear_modulus / 3.0
    }
    /// Young's modulus E = μ(3λ + 2μ)/(λ + μ).
    pub fn young_modulus(&self) -> f64 {
        let l = self.lame_lambda;
        let m = self.shear_modulus;
        m * (3.0 * l + 2.0 * m) / (l + m)
    }
    /// Poisson's ratio ν = λ / (2(λ + μ)).
    pub fn poisson_ratio(&self) -> f64 {
        self.lame_lambda / (2.0 * (self.lame_lambda + self.shear_modulus))
    }
    /// Strain energy density from principal stretches.
    ///
    /// W = λ/2 · (Σh_i)² + μ · Σh_i²
    /// where h_i = ln(λ_i).
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let h1 = lambda1.ln();
        let h2 = lambda2.ln();
        let h3 = lambda3.ln();
        let h_sum = h1 + h2 + h3;
        let h_sq = h1 * h1 + h2 * h2 + h3 * h3;
        self.lame_lambda / 2.0 * h_sum * h_sum + self.shear_modulus * h_sq
    }
    /// Strain energy density from deformation gradient F.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
    /// Cauchy stress tensor (principal components) from principal stretches.
    ///
    /// σ_i = (1/J) \[λ ln J + 2μ h_i\]   where J = λ₁λ₂λ₃.
    pub fn principal_cauchy_stress(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> [f64; 3] {
        let h1 = lambda1.ln();
        let h2 = lambda2.ln();
        let h3 = lambda3.ln();
        let j = lambda1 * lambda2 * lambda3;
        let ln_j = j.ln();
        let scale = 1.0 / j;
        [
            scale * (self.lame_lambda * ln_j + 2.0 * self.shear_modulus * h1),
            scale * (self.lame_lambda * ln_j + 2.0 * self.shear_modulus * h2),
            scale * (self.lame_lambda * ln_j + 2.0 * self.shear_modulus * h3),
        ]
    }
}
/// Blatz-Ko hyperelastic model for compressible foam/rubber.
///
/// Strain energy:
///
///   W = (μ/2) \[ (I₁/I₃ - 3) + (I₃^(-1) - 1)/f \]
///
/// where f ∈ (0, 1] is a volume fraction parameter controlling
/// compressibility (f=1 → incompressible rubber, f=0.5 → 50/50 foam).
///
/// The original Blatz-Ko (1962) paper uses two parameter sets; here we
/// implement the general form.
#[derive(Debug, Clone)]
pub struct BlatzKo {
    /// Shear modulus μ (Pa).
    pub shear_modulus: f64,
    /// Volume fraction f ∈ (0, 1].
    pub f: f64,
    /// Poisson's ratio ν (for initial tangent modulus).
    pub poisson_ratio: f64,
}
impl BlatzKo {
    /// Create a new Blatz-Ko model.
    pub fn new(shear_modulus: f64, f: f64, poisson_ratio: f64) -> Self {
        Self {
            shear_modulus,
            f,
            poisson_ratio,
        }
    }
    /// Standard foam model (f=0.5, ν=0.25).
    pub fn foam_default(shear_modulus: f64) -> Self {
        Self::new(shear_modulus, 0.5, 0.25)
    }
    /// Strain energy density from principal stretches.
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let i1 = lambda1 * lambda1 + lambda2 * lambda2 + lambda3 * lambda3;
        let i3 = j * j;
        let w_isochoric = self.f * (i1 / i3.cbrt() - 3.0);
        let w_volumetric =
            (1.0 - self.f) * (1.0 / j.powf(2.0 * self.f / (1.0 - 2.0 * self.poisson_ratio)) - 1.0);
        self.shear_modulus / 2.0 * (w_isochoric + w_volumetric)
    }
    /// Strain energy from deformation gradient.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
    /// Small-strain Young's modulus E = 2μ(1+ν).
    pub fn young_modulus(&self) -> f64 {
        2.0 * self.shear_modulus * (1.0 + self.poisson_ratio)
    }
}
/// Varga hyperelastic model for slightly compressible rubber.
///
/// Strain energy:
///
///   W = 2μ (λ₁ + λ₂ + λ₃ - 3)
///
/// (purely in terms of principal stretches, no volumetric penalty term).
/// Near-incompressibility is enforced by a bulk penalty.
#[derive(Debug, Clone)]
pub struct Varga {
    /// Shear modulus μ (Pa).
    pub shear_modulus: f64,
    /// Bulk modulus (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl Varga {
    /// Create a new Varga model.
    pub fn new(shear_modulus: f64, bulk_modulus: f64) -> Self {
        Self {
            shear_modulus,
            bulk_modulus,
        }
    }
    /// Strain energy density from principal stretches.
    ///
    /// W = 2μ(λ₁ + λ₂ + λ₃ - 3) + K/2 (J-1)²
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let isochoric = 2.0 * self.shear_modulus * (lambda1 + lambda2 + lambda3 - 3.0);
        let volumetric = self.bulk_modulus / 2.0 * (j - 1.0).powi(2);
        isochoric + volumetric
    }
    /// Strain energy from deformation gradient.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
    /// Uniaxial Cauchy stress (incompressible):
    ///
    /// σ = 2μ(1 - λ^(-3/2))
    pub fn uniaxial_stress_incompressible(&self, lambda: f64) -> f64 {
        2.0 * self.shear_modulus * (1.0 - lambda.powf(-1.5))
    }
    /// Small-strain shear modulus.
    pub fn shear_modulus(&self) -> f64 {
        self.shear_modulus
    }
}
/// Holzapfel-Ogden (HO) model for fibre-reinforced biological tissues.
///
/// Strain energy:
///
///   W = a/(2b) exp(b(I₁-3)) - a ln J
///     + Σ_{i=1,2} \[a_i/(2b_i) (exp(b_i(I₄ᵢ-1)²) - 1)\] + κ/2 (J-1)²
///
/// where a, b > 0 are isotropic constants, a_i, b_i > 0 are fibre
/// constants, I₄ᵢ = a_i · C · a_i (fibre stretch squared), and J = det F.
///
/// Two fibre families with directions a₁ and a₂ are supported.
#[derive(Debug, Clone)]
pub struct HolzapfelOgden {
    /// Isotropic stiffness parameter a (Pa).
    pub a: f64,
    /// Isotropic exponential parameter b (dimensionless).
    pub b: f64,
    /// Fibre stiffness a₁ (Pa) for family 1.
    pub a_f1: f64,
    /// Fibre exponential b₁ (dimensionless) for family 1.
    pub b_f1: f64,
    /// Fibre stiffness a₂ (Pa) for family 2.
    pub a_f2: f64,
    /// Fibre exponential b₂ (dimensionless) for family 2.
    pub b_f2: f64,
    /// Fibre direction 1 (unit vector in reference configuration).
    pub a1: [f64; 3],
    /// Fibre direction 2 (unit vector in reference configuration).
    pub a2: [f64; 3],
    /// Bulk modulus κ (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl HolzapfelOgden {
    /// Create a new Holzapfel-Ogden model.
    pub fn new(
        a: f64,
        b: f64,
        a_f1: f64,
        b_f1: f64,
        a_f2: f64,
        b_f2: f64,
        a1: [f64; 3],
        a2: [f64; 3],
        bulk_modulus: f64,
    ) -> Self {
        Self {
            a,
            b,
            a_f1,
            b_f1,
            a_f2,
            b_f2,
            a1,
            a2,
            bulk_modulus,
        }
    }
    /// Create a single-family model (a2 is set to a1 with zero stiffness).
    pub fn single_family(a: f64, b: f64, a_f: f64, b_f: f64, fibre: [f64; 3], bulk: f64) -> Self {
        Self::new(a, b, a_f, b_f, 0.0, 1.0, fibre, fibre, bulk)
    }
    /// Myocardium (heart muscle) default parameters (Holzapfel 2009).
    pub fn myocardium_default() -> Self {
        Self::new(
            59.0,
            8.023,
            18472.0,
            16.026,
            2481.0,
            11.12,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            1.0e6,
        )
    }
    /// Compute fibre stretch squared I₄ = (Fa)·(Fa) for a given F and fibre direction a.
    fn fibre_stretch_squared(f: &[[f64; 3]; 3], a: &[f64; 3]) -> f64 {
        let fa: [f64; 3] =
            std::array::from_fn(|i| f[i].iter().zip(a.iter()).map(|(fij, aj)| fij * aj).sum());
        fa[0] * fa[0] + fa[1] * fa[1] + fa[2] * fa[2]
    }
    /// Strain energy density given deformation gradient F.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let j = det3(f);
        let (i1, _i2, _i3) = deformation_invariants(f);
        let w_iso = self.a / (2.0 * self.b) * ((self.b * (i1 - 3.0)).exp() - 1.0) - self.a * j.ln();
        let i4_1 = Self::fibre_stretch_squared(f, &self.a1);
        let w_f1 = if i4_1 > 1.0 {
            self.a_f1 / (2.0 * self.b_f1) * ((self.b_f1 * (i4_1 - 1.0).powi(2)).exp() - 1.0)
        } else {
            0.0
        };
        let i4_2 = Self::fibre_stretch_squared(f, &self.a2);
        let w_f2 = if i4_2 > 1.0 {
            self.a_f2 / (2.0 * self.b_f2) * ((self.b_f2 * (i4_2 - 1.0).powi(2)).exp() - 1.0)
        } else {
            0.0
        };
        let w_vol = self.bulk_modulus / 2.0 * (j - 1.0).powi(2);
        w_iso + w_f1 + w_f2 + w_vol
    }
    /// Check if a fibre family is active (in tension) given F and fibre index (1 or 2).
    pub fn fibre_active(&self, f: &[[f64; 3]; 3], family: usize) -> bool {
        let a = if family == 1 { &self.a1 } else { &self.a2 };
        Self::fibre_stretch_squared(f, a) > 1.0
    }
}
/// Ogden hyperelastic material model.
///
/// Strain energy: W = sum_p (mu_p / alpha_p) * (lambda1^alpha_p + lambda2^alpha_p + lambda3^alpha_p - 3)
///                  + D1 * (J - 1)^2
///
/// Where lambda_i are the principal stretches.
#[derive(Debug, Clone)]
pub struct Ogden {
    /// Material constants mu_p (Pa).
    pub mu: Vec<f64>,
    /// Exponents alpha_p.
    pub alpha: Vec<f64>,
    /// Bulk modulus (Pa) for volumetric penalty.
    pub bulk_modulus: f64,
}
impl Ogden {
    /// Create a new Ogden material with N terms.
    pub fn new(mu: Vec<f64>, alpha: Vec<f64>, bulk_modulus: f64) -> Self {
        assert_eq!(mu.len(), alpha.len(), "mu and alpha must have same length");
        Self {
            mu,
            alpha,
            bulk_modulus,
        }
    }
    /// Create a 1-term Ogden model (equivalent to Neo-Hookean when alpha=2).
    pub fn one_term(mu: f64, alpha: f64, bulk_modulus: f64) -> Self {
        Self::new(vec![mu], vec![alpha], bulk_modulus)
    }
    /// Create a 3-term Ogden model for rubber.
    pub fn rubber_default() -> Self {
        Self::new(vec![6.3e5, 1.2e3, -1.0e4], vec![1.3, 5.0, -2.0], 1.0e9)
    }
    /// Initial shear modulus: mu = 0.5 * sum(mu_p * alpha_p).
    pub fn shear_modulus(&self) -> f64 {
        0.5 * self
            .mu
            .iter()
            .zip(self.alpha.iter())
            .map(|(m, a)| m * a)
            .sum::<f64>()
    }
    /// Strain energy density given principal stretches.
    pub fn strain_energy_from_stretches(&self, lambda1: f64, lambda2: f64, lambda3: f64) -> f64 {
        let j = lambda1 * lambda2 * lambda3;
        let j_13 = j.powf(-1.0 / 3.0);
        let l1 = j_13 * lambda1;
        let l2 = j_13 * lambda2;
        let l3 = j_13 * lambda3;
        let mut w = 0.0;
        for (mu_p, alpha_p) in self.mu.iter().zip(self.alpha.iter()) {
            w += (mu_p / alpha_p)
                * (l1.powf(*alpha_p) + l2.powf(*alpha_p) + l3.powf(*alpha_p) - 3.0);
        }
        let d1 = self.bulk_modulus / 2.0;
        w += d1 * (j - 1.0).powi(2);
        w
    }
    /// Strain energy density from deformation gradient.
    pub fn strain_energy_density(&self, f: &[[f64; 3]; 3]) -> f64 {
        let (l1, l2, l3) = principal_stretches(f);
        self.strain_energy_from_stretches(l1, l2, l3)
    }
}
