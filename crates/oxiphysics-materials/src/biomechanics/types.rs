//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
use std::f64::consts::PI;

/// Transversely isotropic meniscus model (fibrocartilage).
#[derive(Debug, Clone)]
pub struct MeniscusModel {
    /// Circumferential (hoop) modulus \[Pa\].
    pub e_circumferential: f64,
    /// Radial modulus \[Pa\].
    pub e_radial: f64,
    /// Shear modulus \[Pa\].
    pub shear_modulus: f64,
    /// Permeability \[m^4/Ns\].
    pub permeability: f64,
    /// Fluid fraction.
    pub fluid_fraction: f64,
}
impl MeniscusModel {
    /// Create a new meniscus model.
    pub fn new(
        e_circumferential: f64,
        e_radial: f64,
        shear_modulus: f64,
        permeability: f64,
        fluid_fraction: f64,
    ) -> Self {
        Self {
            e_circumferential,
            e_radial,
            shear_modulus,
            permeability,
            fluid_fraction,
        }
    }
    /// Typical medial meniscus (knee).
    pub fn medial_meniscus() -> Self {
        Self::new(150e6, 1.4e6, 0.12e6, 1.2e-15, 0.72)
    }
    /// Hoop stress from joint reaction force F \[N\] acting radially.
    pub fn hoop_stress(&self, joint_force: f64, area: f64) -> f64 {
        joint_force / area
    }
}
/// Cortical bone model using Carter-Hayes density-elasticity relationships.
#[derive(Debug, Clone)]
pub struct CorticalBone {
    /// Apparent density ρ_app \[g/cm^3\].
    pub density: f64,
    /// Ash fraction (mineral content) \[dimensionless\].
    pub ash_fraction: f64,
}
impl CorticalBone {
    /// Create a new cortical bone model.
    pub fn new(density: f64, ash_fraction: f64) -> Self {
        Self {
            density,
            ash_fraction,
        }
    }
    /// Typical human femur cortical bone.
    pub fn femur_cortical() -> Self {
        Self::new(1.85, 0.62)
    }
    /// Elastic modulus E \[MPa\] via Carter-Hayes relationship.
    ///
    /// E = 2017 * ρ_app^2.5 \[MPa\] (Carter & Hayes 1977, density in g/cm^3)
    pub fn elastic_modulus(&self) -> f64 {
        2017.0 * self.density.powf(2.5)
    }
    /// Compressive strength \[MPa\] via density power law.
    ///
    /// σ_c = 68 * ρ_app^1.6 \[MPa\]
    pub fn strength_compressive(&self) -> f64 {
        68.0 * self.density.powf(1.6)
    }
    /// Tensile strength \[MPa\] via density power law.
    ///
    /// σ_t = 131 * ρ_app \[MPa\] (approximate linear relationship)
    pub fn strength_tensile(&self) -> f64 {
        131.0 * self.density
    }
    /// Mineral content contribution to stiffness (Katz & Ukraincik).
    pub fn mineral_stiffness_contribution(&self) -> f64 {
        self.elastic_modulus() * (0.6 + 0.4 * self.ash_fraction)
    }
}
/// Muscle model using Hill's three-element (active + passive + series) model.
#[derive(Debug, Clone)]
pub struct MuscleModel {
    /// Maximum isometric force \[N\].
    pub f_max: f64,
    /// Optimal fiber length \[m\].
    pub l_opt: f64,
    /// Maximum shortening velocity \[m/s\].
    pub v_max: f64,
    /// Pennation angle at optimal length \[radians\].
    pub pennation_angle: f64,
}
impl MuscleModel {
    /// Create a new muscle model.
    pub fn new(f_max: f64, l_opt: f64, v_max: f64, pennation_angle: f64) -> Self {
        Self {
            f_max,
            l_opt,
            v_max,
            pennation_angle,
        }
    }
    /// Typical human tibialis anterior muscle.
    pub fn tibialis_anterior() -> Self {
        Self::new(600.0, 0.068, 0.51, 0.17)
    }
    /// Typical human quadriceps (rectus femoris).
    pub fn quadriceps() -> Self {
        Self::new(1169.0, 0.084, 0.63, 0.17)
    }
    /// Active force-length relationship f_l(l̃) where l̃ = l/l_opt.
    ///
    /// Gaussian: f_l = exp(-((l̃ - 1)/ω)^2), ω ≈ 0.45
    pub fn active_force_length(&self, l_normalized: f64) -> f64 {
        let omega = 0.45;
        let x = (l_normalized - 1.0) / omega;
        (-x * x).exp()
    }
    /// Passive force-length relationship f_p(l̃).
    ///
    /// Exponential rise for l̃ > 1: f_p = (exp(k_p*(l̃-1)/ε_0) - 1) / (exp(k_p) - 1)
    pub fn passive_force_length(&self, l_normalized: f64) -> f64 {
        if l_normalized <= 1.0 {
            return 0.0;
        }
        let k_p = 5.0;
        let eps_0 = 0.6;
        let numerator = (k_p * (l_normalized - 1.0) / eps_0).exp() - 1.0;
        let denominator = k_p.exp() - 1.0;
        (numerator / denominator).min(1.0)
    }
    /// Force-velocity relationship f_v(ṽ) where ṽ = v/v_max (Hill's equation).
    ///
    /// Concentric: f_v = (1 - ṽ)/(1 + ṽ/a_v), a_v ≈ 0.25
    /// Eccentric:  f_v = (1.8 - 0.8*(1+ṽ)/(1-7.56*ṽ))
    pub fn force_velocity(&self, v_normalized: f64) -> f64 {
        if v_normalized <= 0.0 {
            let a_v = 0.25;
            let v = v_normalized.max(-1.0 + 1e-9);
            ((1.0 + v) / (1.0 - v / a_v)).max(0.0)
        } else {
            let v = v_normalized.min(1.0);
            (1.8 - 0.8 * (1.0 + v) / (1.0 - 7.56 * v / 6.0)).max(0.0)
        }
    }
    /// Compute muscle force \[N\] given fiber length \[m\], velocity \[m/s\], and activation a ∈ \[0,1\].
    pub fn muscle_force(&self, l: f64, v: f64, activation: f64) -> f64 {
        let l_norm = l / self.l_opt;
        let v_norm = v / self.v_max;
        let fl = self.active_force_length(l_norm);
        let fp = self.passive_force_length(l_norm);
        let fv = self.force_velocity(v_norm);
        let a = activation.clamp(0.0, 1.0);
        let cos_psi = self.pennation_angle.cos();
        (a * fl * fv + fp) * self.f_max * cos_psi
    }
    /// Maximum power output \[W\].
    pub fn max_power(&self) -> f64 {
        self.f_max * self.v_max / 4.0
    }
}
/// Blood vessel model using thin-walled vessel theory (Laplace's law).
#[derive(Debug, Clone)]
pub struct BloodVesselModel {
    /// Reference (unloaded) inner radius R0 \[m\].
    pub r0: f64,
    /// Wall thickness h0 \[m\].
    pub h0: f64,
    /// Vessel wall elastic modulus E \[Pa\].
    pub e_wall: f64,
    /// Poisson's ratio of vessel wall.
    pub nu_wall: f64,
}
impl BloodVesselModel {
    /// Create a new blood vessel model.
    pub fn new(r0: f64, h0: f64, e_wall: f64, nu_wall: f64) -> Self {
        Self {
            r0,
            h0,
            e_wall,
            nu_wall,
        }
    }
    /// Typical human aorta.
    pub fn aorta() -> Self {
        Self::new(12.5e-3, 2.0e-3, 500e3, 0.45)
    }
    /// Typical femoral artery.
    pub fn femoral_artery() -> Self {
        Self::new(3.5e-3, 0.6e-3, 800e3, 0.45)
    }
    /// Circumferential (hoop) wall stress via Laplace's law \[Pa\].
    ///
    /// σ_θ = P * R / h
    pub fn laplace_wall_stress(&self, p_transmural: f64) -> f64 {
        p_transmural * self.r0 / self.h0
    }
    /// Pressure-radius compliance C \[m/Pa\] = dR/dP.
    ///
    /// Approximate: C = R0 * (1 - ν^2) / (E * h/R0)
    pub fn compliance(&self, _p: f64) -> f64 {
        self.r0 * (1.0 - self.nu_wall * self.nu_wall) / (self.e_wall * self.h0 / self.r0)
    }
    /// Pulse wave velocity \[m/s\] (Moens-Korteweg equation).
    ///
    /// c = sqrt(E * h / (2 * ρ * R))
    pub fn pulse_wave_velocity(&self, rho: f64) -> f64 {
        (self.e_wall * self.h0 / (2.0 * rho * self.r0)).sqrt()
    }
    /// Pressure-dependent radius (linearized expansion) \[m\].
    pub fn radius_at_pressure(&self, p: f64) -> f64 {
        let c = self.compliance(p);
        self.r0 + c * p
    }
}
/// Corneal biomechanical model (Ogden-type hyperelastic).
#[derive(Debug, Clone)]
pub struct CorneaModel {
    /// Tangent modulus at small strain \[Pa\].
    pub e_small: f64,
    /// Tangent modulus at large strain (linearized) \[Pa\].
    pub e_large: f64,
    /// Transition strain.
    pub transition_strain: f64,
    /// Corneal thickness (central) \[m\].
    pub cct: f64,
    /// Intraocular pressure \[Pa\].
    pub iop: f64,
}
impl CorneaModel {
    /// Create a new cornea model.
    pub fn new(e_small: f64, e_large: f64, transition_strain: f64, cct: f64, iop: f64) -> Self {
        Self {
            e_small,
            e_large,
            transition_strain,
            cct,
            iop,
        }
    }
    /// Typical healthy human cornea.
    pub fn healthy_cornea() -> Self {
        Self::new(0.01e6, 0.5e6, 0.05, 550e-6, 2000.0)
    }
    /// Hoop stress from IOP (thin shell approximation) \[Pa\].
    ///
    /// σ = IOP * R / (2 * t)  (sphere approximation R ≈ 7.8 mm)
    pub fn hoop_stress_iop(&self) -> f64 {
        let r = 7.8e-3;
        self.iop * r / (2.0 * self.cct)
    }
    /// Effective stiffness modulus at a given strain level \[Pa\].
    pub fn tangent_modulus(&self, strain: f64) -> f64 {
        if strain < self.transition_strain {
            self.e_small
        } else {
            self.e_large
        }
    }
}
/// Failure mode of a biological tissue.
#[derive(Debug, Clone, PartialEq)]
pub enum FailureMode {
    /// No failure predicted.
    NoFailure,
    /// Tensile (ultimate tensile strength exceeded).
    TensileFailure,
    /// Compressive (ultimate compressive strength exceeded).
    CompressiveFailure,
    /// Shear (ultimate shear strength exceeded).
    ShearFailure,
    /// Fatigue (cyclic loading below monotonic UTS).
    FatigueFailure,
    /// Fracture (crack propagation, KIc exceeded).
    FractureFailure,
}
/// Anisotropic heart valve leaflet material model.
#[derive(Debug, Clone)]
pub struct HeartValveTissue {
    /// Circumferential (preferred fiber) modulus \[Pa\].
    pub e_circ: f64,
    /// Radial modulus \[Pa\].
    pub e_radial: f64,
    /// Flexural rigidity \[Pa·m^3\].
    pub flexural_rigidity: f64,
    /// Tissue thickness \[m\].
    pub thickness: f64,
}
impl HeartValveTissue {
    /// Create a new heart valve tissue model.
    pub fn new(e_circ: f64, e_radial: f64, flexural_rigidity: f64, thickness: f64) -> Self {
        Self {
            e_circ,
            e_radial,
            flexural_rigidity,
            thickness,
        }
    }
    /// Typical aortic valve leaflet.
    pub fn aortic_valve() -> Self {
        Self::new(15e6, 2e6, 0.001, 0.5e-3)
    }
    /// Bending stiffness D = E*t^3/(12*(1-ν^2)) \[Pa·m^3\].
    pub fn bending_stiffness(&self, nu: f64) -> f64 {
        self.e_circ * self.thickness.powi(3) / (12.0 * (1.0 - nu * nu))
    }
    /// Anisotropy ratio E_circ / E_radial.
    pub fn anisotropy_ratio(&self) -> f64 {
        self.e_circ / self.e_radial
    }
}
/// Soft tissue material using the Fung-type (exponential) hyperelastic model.
///
/// Strain energy: W = c1 * (exp(Q) - 1) where Q = c2*(I1-3) + alpha*(I2-3) + beta*(J-1)^2
#[derive(Debug, Clone)]
pub struct SoftTissueMaterial {
    /// Material parameter c1 \[Pa\] — scaling of exponential term.
    pub c1: f64,
    /// Material parameter c2 \[dimensionless\] — controls I1 contribution.
    pub c2: f64,
    /// Material parameter alpha \[dimensionless\] — controls I2 contribution.
    pub alpha: f64,
    /// Material parameter beta \[Pa\] — volumetric penalty.
    pub beta: f64,
}
impl SoftTissueMaterial {
    /// Create a new Fung-type soft tissue model.
    pub fn new(c1: f64, c2: f64, alpha: f64, beta: f64) -> Self {
        Self {
            c1,
            c2,
            alpha,
            beta,
        }
    }
    /// Typical liver soft tissue parameters.
    pub fn liver() -> Self {
        Self::new(350.0, 0.57, 0.1, 1000.0)
    }
    /// Typical brain tissue parameters (gray matter).
    pub fn brain() -> Self {
        Self::new(264.0, 0.4, 0.05, 500.0)
    }
    /// Compute the Q exponent argument from invariants I1, I2, J.
    pub fn q_exponent(&self, i1: f64, i2: f64, j: f64) -> f64 {
        self.c2 * (i1 - 3.0) + self.alpha * (i2 - 3.0) + self.beta * (j - 1.0).powi(2)
    }
    /// Compute strain energy density W \[Pa\] from invariants I1, I2, J.
    ///
    /// W = c1 * (exp(Q) - 1)
    pub fn strain_energy_density(&self, i1: f64, i2: f64, j: f64) -> f64 {
        let q = self.q_exponent(i1, i2, j);
        self.c1 * (q.exp() - 1.0)
    }
    /// Compute the Cauchy stress tensor \[Pa\] from a Green-Lagrange strain tensor E.
    ///
    /// The deformation gradient is approximated as F = I + 2E (small-to-moderate strain).
    /// Returns the Cauchy stress as \[\[f64;3\];3\].
    pub fn cauchy_stress_green_lagrange(&self, strain: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let id = identity3();
        let mut f = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                f[i][j] = id[i][j] + strain[i][j];
            }
        }
        let c = right_cauchy_green(&f);
        let i1 = trace3(&c);
        let i2 = second_invariant_c(&c);
        let j = det3(&f).abs().max(1e-14);
        let q = self.q_exponent(i1, i2, j);
        let exp_q = q.exp();
        let dw_di1 = self.c1 * self.c2 * exp_q;
        let dw_di2 = self.c1 * self.alpha * exp_q;
        let dw_dj = self.c1 * self.beta * 2.0 * (j - 1.0) * exp_q;
        let ft = transpose3(&f);
        let b = matmul3(&f, &ft);
        let b2 = matmul3(&b, &b);
        let i1_b = trace3(&b);
        let term1 = scale3(2.0 / j * dw_di1, &b);
        let term2_inner = add3(&scale3(i1_b, &b), &scale3(-1.0, &b2));
        let term2 = scale3(2.0 / j * dw_di2, &term2_inner);
        let term3 = scale3(dw_dj, &id);
        let sigma1 = add3(&term1, &term2);
        add3(&sigma1, &term3)
    }
}
/// Tendon and ligament model using quasi-linear viscoelastic (QLV) theory (Fung 1972).
///
/// Models the nonlinear elastic toe-to-linear response and stress relaxation.
#[derive(Debug, Clone)]
pub struct TendonLigamentModel {
    /// Toe-region tangent modulus \[Pa\].
    pub toe_region_modulus: f64,
    /// Linear-region (post-transition) modulus \[Pa\].
    pub linear_modulus: f64,
    /// Transition strain (toe-to-linear) \[dimensionless\].
    pub transition_strain: f64,
    /// QLV reduced relaxation coefficient C \[dimensionless\].
    pub c_relax: f64,
    /// QLV fast time constant τ1 \[s\].
    pub tau1: f64,
    /// QLV slow time constant τ2 \[s\].
    pub tau2: f64,
}
impl TendonLigamentModel {
    /// Create a new tendon/ligament model.
    pub fn new(
        toe_region_modulus: f64,
        linear_modulus: f64,
        transition_strain: f64,
        c_relax: f64,
        tau1: f64,
        tau2: f64,
    ) -> Self {
        Self {
            toe_region_modulus,
            linear_modulus,
            transition_strain,
            c_relax,
            tau1,
            tau2,
        }
    }
    /// Typical human patellar tendon.
    pub fn patellar_tendon() -> Self {
        Self::new(200e6, 1500e6, 0.02, 0.4, 0.01, 10.0)
    }
    /// Typical human ACL (anterior cruciate ligament).
    pub fn acl() -> Self {
        Self::new(50e6, 300e6, 0.03, 0.35, 0.02, 15.0)
    }
    /// Nonlinear stress-strain response σ(ε) \[Pa\].
    ///
    /// Toe region (ε < ε_t): σ = E_toe * ε * (ε / ε_t)
    /// Linear region (ε ≥ ε_t): σ = E_toe*ε_t + E_linear*(ε - ε_t)
    pub fn stress_strain_nonlinear(&self, strain: f64) -> f64 {
        if strain < 0.0 {
            return 0.0;
        }
        if strain < self.transition_strain {
            self.toe_region_modulus * strain * strain / self.transition_strain
        } else {
            let sigma_transition = self.toe_region_modulus * self.transition_strain;
            sigma_transition + self.linear_modulus * (strain - self.transition_strain)
        }
    }
    /// Reduced relaxation function G(t) \[dimensionless\] from QLV theory.
    ///
    /// G(t) = (1 + C * (E1(t/τ1) - E1(t/τ2))) / (1 + C*ln(τ2/τ1))
    /// Approximated using a two-exponential Prony series.
    pub fn relaxation_function(&self, t: f64) -> f64 {
        if t <= 0.0 {
            return 1.0;
        }
        let g1 = (-t / self.tau1).exp();
        let g2 = (-t / self.tau2).exp();
        let denom = 1.0 + self.c_relax * (self.tau2 / self.tau1).ln();
        let amplitude = self.c_relax / denom;
        let g_inf = 1.0 - amplitude;
        g_inf + 0.5 * amplitude * g1 + 0.5 * amplitude * g2
    }
    /// Viscoelastic stress via convolution integral (discrete time step approximation).
    ///
    /// σ(t) = G(t) * σ_e(ε) where σ_e is the instantaneous (elastic) stress.
    pub fn viscoelastic_stress(&self, strain: f64, t: f64) -> f64 {
        let sigma_elastic = self.stress_strain_nonlinear(strain);
        let g_t = self.relaxation_function(t);
        sigma_elastic * g_t
    }
    /// Ultimate tensile strength estimate from linear modulus and failure strain.
    pub fn ultimate_stress(&self, failure_strain: f64) -> f64 {
        self.stress_strain_nonlinear(failure_strain)
    }
}
/// Trabecular (cancellous) bone model.
#[derive(Debug, Clone)]
pub struct TrabecularBone {
    /// Bone volume fraction BV/TV \[dimensionless\] ∈ (0, 1).
    pub bv_tv: f64,
    /// Ash fraction \[dimensionless\].
    pub ash_fraction: f64,
}
impl TrabecularBone {
    /// Create a new trabecular bone model.
    pub fn new(bv_tv: f64, ash_fraction: f64) -> Self {
        Self {
            bv_tv,
            ash_fraction,
        }
    }
    /// Typical vertebral trabecular bone.
    pub fn vertebral() -> Self {
        Self::new(0.15, 0.55)
    }
    /// Elastic modulus \[MPa\] of trabecular bone via BV/TV.
    ///
    /// E = 1904 * BV/TV^1.64 \[MPa\] (Kopperdahl & Keaveny 1998)
    pub fn elastic_modulus_trabecular(&self) -> f64 {
        1904.0 * self.bv_tv.powf(1.64)
    }
    /// Yield stress \[MPa\] in compression.
    ///
    /// σ_y = 24.8 * BV/TV^1.6 \[MPa\]
    pub fn yield_stress_trabecular(&self) -> f64 {
        24.8 * self.bv_tv.powf(1.6)
    }
    /// Fabric tensor isotropy ratio (0 = fully anisotropic, 1 = isotropic).
    pub fn isotropy_ratio(&self) -> f64 {
        (self.bv_tv / 0.35).min(1.0)
    }
}
/// A simpler ligament wrapping model for joint constraint forces.
#[derive(Debug, Clone)]
pub struct LigamentWrapModel {
    /// Rest (zero-force) length \[m\].
    pub rest_length: f64,
    /// Stiffness k \[N/m\].
    pub stiffness: f64,
    /// Slack length (below this: no force) \[m\].
    pub slack_length: f64,
}
impl LigamentWrapModel {
    /// Create a new ligament wrap model.
    pub fn new(rest_length: f64, stiffness: f64, slack_length: f64) -> Self {
        Self {
            rest_length,
            stiffness,
            slack_length,
        }
    }
    /// Restraint force \[N\] for a given current length \[m\].
    pub fn restraint_force(&self, current_length: f64) -> f64 {
        if current_length <= self.slack_length {
            0.0
        } else {
            self.stiffness * (current_length - self.slack_length)
        }
    }
    /// Strain = (l - l0) / l0.
    pub fn strain(&self, current_length: f64) -> f64 {
        (current_length - self.rest_length) / self.rest_length
    }
}
/// Failure criteria for biological soft and hard tissues.
///
/// Implements:
/// - Maximum principal stress criterion (brittle fracture analog)
/// - Tsai-Wu criterion (anisotropic composite analog for fibrous tissues)
/// - Bone fracture criterion based on stress intensity factor KI
/// - Soft tissue tear criterion based on critical energy release rate Gc
#[derive(Debug, Clone)]
pub struct TissueFailureCriteria {
    /// Ultimate tensile strength \[Pa\].
    pub uts: f64,
    /// Ultimate compressive strength \[Pa\] (positive value).
    pub ucs: f64,
    /// Ultimate shear strength \[Pa\].
    pub uss: f64,
    /// Fatigue strength coefficient (Basquin σ'f) \[Pa\].
    pub fatigue_coeff: f64,
    /// Fatigue exponent (Basquin b, negative).
    pub fatigue_exp: f64,
    /// Critical stress intensity factor KIc \[Pa√m\].
    pub k_ic: f64,
    /// Critical energy release rate Gc \[J/m²\].
    pub g_c: f64,
    /// Tsai-Wu interaction coefficient F12* (normalized).
    pub tsai_wu_f12_star: f64,
}
impl TissueFailureCriteria {
    /// Cortical bone failure criteria (Reilly & Burstein 1975).
    pub fn cortical_bone() -> Self {
        Self {
            uts: 133e6,
            ucs: 193e6,
            uss: 68e6,
            fatigue_coeff: 140e6,
            fatigue_exp: -0.08,
            k_ic: 2.2e6,
            g_c: 600.0,
            tsai_wu_f12_star: -0.5,
        }
    }
    /// Articular cartilage failure criteria.
    pub fn articular_cartilage() -> Self {
        Self {
            uts: 15e6,
            ucs: 5e6,
            uss: 4e6,
            fatigue_coeff: 10e6,
            fatigue_exp: -0.10,
            k_ic: 0.2e6,
            g_c: 100.0,
            tsai_wu_f12_star: -0.5,
        }
    }
    /// Ligament/tendon failure criteria.
    pub fn ligament() -> Self {
        Self {
            uts: 38e6,
            ucs: 1e6,
            uss: 10e6,
            fatigue_coeff: 30e6,
            fatigue_exp: -0.09,
            k_ic: 0.5e6,
            g_c: 200.0,
            tsai_wu_f12_star: -0.5,
        }
    }
    /// Maximum principal stress criterion.
    ///
    /// Returns the [`FailureMode`] given the three principal stresses \[Pa\].
    pub fn max_principal_stress(&self, s1: f64, s2: f64, s3: f64) -> FailureMode {
        let max_tensile = s1.max(s2).max(s3);
        let max_compressive = s1.min(s2).min(s3);
        if max_tensile >= self.uts {
            FailureMode::TensileFailure
        } else if max_compressive <= -self.ucs {
            FailureMode::CompressiveFailure
        } else {
            FailureMode::NoFailure
        }
    }
    /// Safety factor against tensile failure: UTS / σ_max.
    pub fn tensile_safety_factor(&self, sigma_max: f64) -> f64 {
        if sigma_max <= 0.0 {
            f64::INFINITY
        } else {
            self.uts / sigma_max
        }
    }
    /// Tsai-Wu criterion for fibrous tissue (fiber direction = 1, transverse = 2).
    ///
    /// F1·σ1 + F2·σ2 + F11·σ1² + F22·σ2² + 2·F12·σ1·σ2 + F66·τ12² = 1 at failure.
    /// Returns the failure index (< 1 = safe, ≥ 1 = failure).
    pub fn tsai_wu_index(&self, sigma1: f64, sigma2: f64, tau12: f64) -> f64 {
        let f1 = 1.0 / self.uts - 1.0 / self.ucs;
        let f2 = f1;
        let f11 = 1.0 / (self.uts * self.ucs);
        let f22 = f11;
        let f66 = 1.0 / (self.uss * self.uss);
        let f12 = self.tsai_wu_f12_star * (f11 * f22).sqrt();
        f1 * sigma1
            + f2 * sigma2
            + f11 * sigma1 * sigma1
            + f22 * sigma2 * sigma2
            + 2.0 * f12 * sigma1 * sigma2
            + f66 * tau12 * tau12
    }
    /// Fatigue life in cycles (Basquin law): σ_a = σ'f · (2N)^b.
    ///
    /// Returns the number of half-cycles to failure (reversals 2N_f).
    pub fn fatigue_life_cycles(&self, stress_amplitude: f64) -> f64 {
        if stress_amplitude <= 0.0 {
            return f64::INFINITY;
        }
        let ratio = stress_amplitude / self.fatigue_coeff;
        if ratio <= 0.0 {
            return f64::INFINITY;
        }
        let two_n = ratio.powf(1.0 / self.fatigue_exp);
        two_n / 2.0
    }
    /// Fracture criterion: compare mode-I stress intensity factor KI \[Pa√m\]
    /// with the critical value KIc.  Returns [`FailureMode::FractureFailure`]
    /// if KI ≥ KIc.
    pub fn fracture_criterion(&self, k_i: f64) -> FailureMode {
        if k_i >= self.k_ic {
            FailureMode::FractureFailure
        } else {
            FailureMode::NoFailure
        }
    }
    /// Energy release rate G \[J/m²\] from stress intensity factor (plane strain).
    ///
    /// G = KI² * (1 - ν²) / E.
    pub fn energy_release_rate(&self, k_i: f64, e_modulus: f64, poisson: f64) -> f64 {
        k_i * k_i * (1.0 - poisson * poisson) / e_modulus
    }
    /// Check whether critical energy release rate Gc is exceeded.
    pub fn toughness_criterion(&self, k_i: f64, e_modulus: f64, poisson: f64) -> FailureMode {
        let g = self.energy_release_rate(k_i, e_modulus, poisson);
        if g >= self.g_c {
            FailureMode::FractureFailure
        } else {
            FailureMode::NoFailure
        }
    }
}
/// Biphasic cartilage model (Mow et al.).
///
/// Treats cartilage as a mixture of solid and fluid phases.
#[derive(Debug, Clone)]
pub struct CartilageModel {
    /// Aggregate modulus H_A \[Pa\] — equilibrium compressive stiffness.
    pub solid_modulus_h_a: f64,
    /// Hydraulic permeability k \[m^4/Ns\].
    pub permeability: f64,
    /// Fluid volume fraction (porosity) φ \[dimensionless\].
    pub fluid_fraction: f64,
    /// Specimen thickness H \[m\].
    pub thickness: f64,
}
impl CartilageModel {
    /// Create a new biphasic cartilage model.
    pub fn new(
        solid_modulus_h_a: f64,
        permeability: f64,
        fluid_fraction: f64,
        thickness: f64,
    ) -> Self {
        Self {
            solid_modulus_h_a,
            permeability,
            fluid_fraction,
            thickness,
        }
    }
    /// Typical articular cartilage (knee).
    pub fn articular_cartilage() -> Self {
        Self::new(0.5e6, 1e-15, 0.75, 2e-3)
    }
    /// Consolidation coefficient c_v = H_A * k \[m^2/s\].
    pub fn consolidation_coefficient(&self) -> f64 {
        self.solid_modulus_h_a * self.permeability
    }
    /// Creep compliance J(t) under constant stress σ_0 \[1/Pa\].
    ///
    /// J(t) = (1/H_A) * \[1 - sum_n (8/(π^2(2n-1)^2)) * exp(-cv*(2n-1)^2*π^2*t/(4H^2))\]
    pub fn creep_compliance(&self, t: f64) -> f64 {
        let cv = self.consolidation_coefficient();
        let h = self.thickness;
        let mut sum = 0.0;
        for n in 1..=20usize {
            let m = (2 * n - 1) as f64;
            let coeff = 8.0 / (PI * PI * m * m);
            let exponent = -cv * m * m * PI * PI * t / (4.0 * h * h);
            sum += coeff * exponent.exp();
        }
        (1.0 / self.solid_modulus_h_a) * (1.0 - sum)
    }
    /// Stress relaxation G(t) normalized by initial stress \[dimensionless\].
    ///
    /// G(t) = sum_n (8/(π^2(2n-1)^2)) * exp(-cv*(2n-1)^2*π^2*t/(4H^2))
    pub fn stress_relaxation(&self, t: f64) -> f64 {
        let cv = self.consolidation_coefficient();
        let h = self.thickness;
        let mut sum = 0.0;
        for n in 1..=20usize {
            let m = (2 * n - 1) as f64;
            let coeff = 8.0 / (PI * PI * m * m);
            let exponent = -cv * m * m * PI * PI * t / (4.0 * h * h);
            sum += coeff * exponent.exp();
        }
        sum
    }
}
/// Two-layer skin model with epidermis and dermis.
#[derive(Debug, Clone)]
pub struct SkinModel {
    /// Young's modulus of epidermis \[Pa\].
    pub e_epidermis: f64,
    /// Young's modulus of dermis \[Pa\].
    pub e_dermis: f64,
    /// Poisson's ratio of epidermis.
    pub nu_epidermis: f64,
    /// Poisson's ratio of dermis.
    pub nu_dermis: f64,
    /// Thickness of epidermis \[m\].
    pub thickness_e: f64,
    /// Thickness of dermis \[m\].
    pub thickness_d: f64,
}
impl SkinModel {
    /// Create a new two-layer skin model.
    pub fn new(
        e_epidermis: f64,
        e_dermis: f64,
        nu_epidermis: f64,
        nu_dermis: f64,
        thickness_e: f64,
        thickness_d: f64,
    ) -> Self {
        Self {
            e_epidermis,
            e_dermis,
            nu_epidermis,
            nu_dermis,
            thickness_e,
            thickness_d,
        }
    }
    /// Typical human forearm skin.
    pub fn human_forearm() -> Self {
        Self::new(1e6, 0.08e6, 0.48, 0.49, 0.1e-3, 1.2e-3)
    }
    /// Effective modulus for a given indentation depth \[Pa\].
    ///
    /// Uses depth-dependent weighting: shallow → epidermis dominated, deep → dermis dominated.
    pub fn effective_modulus(&self, indentation: f64) -> f64 {
        let h_e = self.thickness_e;
        let w = (indentation / h_e).min(1.0);

        (1.0 - w) * self.e_epidermis + w * self.e_dermis
    }
    /// Effective Poisson's ratio (depth-dependent).
    pub fn effective_poisson(&self, indentation: f64) -> f64 {
        let w = (indentation / self.thickness_e).min(1.0);
        (1.0 - w) * self.nu_epidermis + w * self.nu_dermis
    }
    /// Hertz contact force \[N\] for spherical indenter of radius r \[m\],
    /// effective modulus E_eff \[Pa\], indentation δ \[m\].
    ///
    /// F = (4/3) * E_eff / (1-ν^2) * sqrt(r) * δ^(3/2)
    pub fn indentation_force_hertz(&self, r_indenter: f64, e_eff: f64, delta: f64) -> f64 {
        let nu = self.effective_poisson(delta);
        let e_reduced = e_eff / (1.0 - nu * nu);
        (4.0 / 3.0) * e_reduced * r_indenter.sqrt() * delta.abs().powf(1.5)
    }
    /// Total skin thickness \[m\].
    pub fn total_thickness(&self) -> f64 {
        self.thickness_e + self.thickness_d
    }
}
/// Simplified intervertebral disc model (annulus fibrosus + nucleus pulposus).
#[derive(Debug, Clone)]
pub struct IntervertebralDiscModel {
    /// Nucleus pulposus bulk modulus \[Pa\].
    pub nucleus_bulk_modulus: f64,
    /// Annulus fibrosus shear modulus \[Pa\].
    pub annulus_shear_modulus: f64,
    /// Disc height \[m\].
    pub height: f64,
    /// Outer disc radius \[m\].
    pub outer_radius: f64,
    /// Inner (nucleus) radius \[m\].
    pub inner_radius: f64,
}
impl IntervertebralDiscModel {
    /// Create a new intervertebral disc model.
    pub fn new(
        nucleus_bulk_modulus: f64,
        annulus_shear_modulus: f64,
        height: f64,
        outer_radius: f64,
        inner_radius: f64,
    ) -> Self {
        Self {
            nucleus_bulk_modulus,
            annulus_shear_modulus,
            height,
            outer_radius,
            inner_radius,
        }
    }
    /// Typical lumbar L4-L5 disc.
    pub fn lumbar_l4_l5() -> Self {
        Self::new(1.72e6, 0.17e6, 11e-3, 23e-3, 12e-3)
    }
    /// Effective compressive stiffness of the disc \[N/m\].
    pub fn compressive_stiffness(&self) -> f64 {
        let area_n = PI * self.inner_radius * self.inner_radius;
        let k_nucleus = self.nucleus_bulk_modulus * area_n / self.height;
        let area_a =
            PI * (self.outer_radius * self.outer_radius - self.inner_radius * self.inner_radius);
        let k_annulus = self.annulus_shear_modulus * area_a / self.height;
        k_nucleus + k_annulus
    }
    /// Intradiscal pressure under axial load F \[N\] → Pa.
    pub fn intradiscal_pressure(&self, axial_force: f64) -> f64 {
        let area_n = PI * self.inner_radius * self.inner_radius;
        axial_force / area_n
    }
}
/// Fiber-reinforced soft tissue with Holzapfel-Gasser-Ogden style fiber model.
///
/// Models a fiber family with mean orientation `a0` and angular dispersion `kappa`.
#[derive(Debug, Clone)]
pub struct FiberReinforcedTissue {
    /// Ground matrix c \[Pa\].
    pub c_matrix: f64,
    /// Mean fiber direction (unit vector) in reference configuration.
    pub a0: [f64; 3],
    /// Fiber dispersion parameter κ ∈ \[0, 1/3\] (0 = aligned, 1/3 = isotropic).
    pub kappa: f64,
    /// Fiber stiffness k1 \[Pa\].
    pub k1: f64,
    /// Fiber exponential k2 \[dimensionless\].
    pub k2: f64,
}
impl FiberReinforcedTissue {
    /// Create a new fiber-reinforced tissue model.
    pub fn new(c_matrix: f64, a0: [f64; 3], kappa: f64, k1: f64, k2: f64) -> Self {
        Self {
            c_matrix,
            a0,
            kappa,
            k1,
            k2,
        }
    }
    /// Typical aortic tissue model (circumferential fiber family).
    pub fn aortic_tissue() -> Self {
        let a0 = [0.0, 1.0, 0.0];
        Self::new(18.4e3, a0, 0.226, 996.6e3, 524.6)
    }
    /// Compute the fiber strain energy W_f \[Pa\] for a given I4 invariant.
    ///
    /// W_f = (k1/(2k2)) * (exp(k2*(I4-1)^2) - 1)  when I4 >= 1 (tension only).
    pub fn fiber_strain_energy(&self, i4: f64) -> f64 {
        if i4 < 1.0 {
            return 0.0;
        }
        let e = i4 - 1.0;
        (self.k1 / (2.0 * self.k2)) * ((self.k2 * e * e).exp() - 1.0)
    }
    /// Compute total strain energy W \[Pa\] from I1 and I4.
    ///
    /// W = c * (I1 - 3) + W_fiber(I4)
    pub fn total_strain_energy(&self, i1: f64, i4: f64) -> f64 {
        let w_matrix = self.c_matrix * (i1 - 3.0);
        let w_fiber = self.fiber_strain_energy(i4);
        w_matrix + w_fiber
    }
    /// Compute the effective I4 = a0 . C . a0 given deformation gradient F.
    pub fn compute_i4(&self, f: &[[f64; 3]; 3]) -> f64 {
        let c = right_cauchy_green(f);
        let mut ca = [0.0f64; 3];
        for (ca_i, c_row) in ca.iter_mut().zip(c.iter()) {
            *ca_i = c_row
                .iter()
                .zip(self.a0.iter())
                .map(|(&c_ij, &a_j)| c_ij * a_j)
                .sum();
        }
        ca.iter()
            .zip(self.a0.iter())
            .map(|(&ca_i, &a0_i)| a0_i * ca_i)
            .sum()
    }
}
/// Cardiac muscle constitutive model (passive + active).
///
/// Combines:
/// - **Passive**: transversely isotropic hyperelastic material (Guccione et al. 1991)
///   with exponential fiber-sheet description.
/// - **Active**: time-varying elastance model (Sagawa) for systolic force.
///
/// Fiber direction is assumed to be aligned with the local e₁ axis.
#[derive(Debug, Clone)]
pub struct CardiacMuscleModel {
    /// Passive stiffness parameter C \[Pa\].
    pub c_passive: f64,
    /// Fiber exponential coefficient bf.
    pub bf: f64,
    /// Cross-fiber coefficient bt.
    pub bt: f64,
    /// Fiber-sheet coupling coefficient bfs.
    pub bfs: f64,
    /// Maximum active stress Tmax \[Pa\].
    pub t_max: f64,
    /// Calcium sensitivity parameter Ca50 \[μM\].
    pub ca50: f64,
    /// Hill coefficient n for calcium activation.
    pub n_hill: f64,
    /// Sarcomere length at zero active stress l0 \[μm\].
    pub l0: f64,
}
impl CardiacMuscleModel {
    /// Default left ventricular myocardium parameters (Guccione et al. 1995).
    pub fn left_ventricle() -> Self {
        Self {
            c_passive: 880.0,
            bf: 18.48,
            bt: 3.58,
            bfs: 1.627,
            t_max: 135_480.0,
            ca50: 0.805,
            n_hill: 2.0,
            l0: 1.58,
        }
    }
    /// Passive strain energy density \[Pa\] from the Guccione exponential model.
    ///
    /// `e_ff`, `e_ss`, `e_nn` are normal Green strains in fiber, sheet, normal directions.
    /// `e_fs` is the fiber-sheet shear strain.
    pub fn passive_strain_energy(&self, e_ff: f64, e_ss: f64, e_nn: f64, e_fs: f64) -> f64 {
        let q = self.bf * e_ff * e_ff
            + self.bt * (e_ss * e_ss + e_nn * e_nn)
            + self.bfs * (2.0 * e_fs * e_fs);
        self.c_passive / 2.0 * (q.exp() - 1.0)
    }
    /// Active fiber stress \[Pa\] given sarcomere length `l` \[μm\] and
    /// intracellular calcium concentration `ca` \[μM\].
    ///
    /// Uses the steady-state Hill equation for calcium activation.
    pub fn active_fiber_stress(&self, l: f64, ca: f64) -> f64 {
        if l <= self.l0 {
            return 0.0;
        }
        let ca50_n = self.ca50.powf(self.n_hill);
        let ca_n = ca.powf(self.n_hill);
        let activation = ca_n / (ca_n + ca50_n);
        let l_factor = ((l - self.l0) / (2.3 - self.l0)).clamp(0.0, 1.0);
        self.t_max * activation * l_factor
    }
    /// Total (passive + active) fiber stress \[Pa\].
    pub fn total_fiber_stress(
        &self,
        e_ff: f64,
        e_ss: f64,
        e_nn: f64,
        e_fs: f64,
        l: f64,
        ca: f64,
    ) -> f64 {
        let w = self.passive_strain_energy(e_ff, e_ss, e_nn, e_fs);
        let dw = {
            let delta = 1e-6;
            let w2 = self.passive_strain_energy(e_ff + delta, e_ss, e_nn, e_fs);
            (w2 - w) / delta
        };
        dw + self.active_fiber_stress(l, ca)
    }
    /// End-systolic pressure–volume relationship slope Ees \[Pa/m³\] (simplified).
    ///
    /// Uses the time-varying elastance: E(t) = (Emax - Emin)*f(t) + Emin.
    pub fn elastance_at_peak(&self, volume: f64, v0: f64) -> f64 {
        self.t_max / (volume - v0).abs().max(1e-6)
    }
    /// Stroke work \[J\] given end-diastolic volume `edv` \[m³\], end-systolic volume `esv` \[m³\],
    /// and mean arterial pressure `map` \[Pa\].
    pub fn stroke_work(&self, edv: f64, esv: f64, _map: f64) -> f64 {
        let sv = edv - esv;
        let ees = self.elastance_at_peak(edv, esv);
        ees * sv * sv / 2.0
    }
}
/// Composite bone model combining cortical and trabecular bone.
#[derive(Debug, Clone)]
pub struct BoneModel {
    /// Cortical bone layer.
    pub cortical: CorticalBone,
    /// Trabecular bone core.
    pub trabecular: TrabecularBone,
    /// Cortical shell thickness \[m\].
    pub cortical_thickness: f64,
    /// Total bone diameter/width \[m\].
    pub total_width: f64,
}
impl BoneModel {
    /// Create a new composite bone model.
    pub fn new(
        cortical: CorticalBone,
        trabecular: TrabecularBone,
        cortical_thickness: f64,
        total_width: f64,
    ) -> Self {
        Self {
            cortical,
            trabecular,
            cortical_thickness,
            total_width,
        }
    }
    /// Typical femoral neck cross-section.
    pub fn femoral_neck() -> Self {
        Self::new(
            CorticalBone::femur_cortical(),
            TrabecularBone::vertebral(),
            2e-3,
            30e-3,
        )
    }
    /// Estimate effective axial modulus via rule of mixtures \[MPa\].
    pub fn effective_axial_modulus(&self) -> f64 {
        let r_total = self.total_width / 2.0;
        let r_inner = r_total - self.cortical_thickness;
        let a_cortical = PI * (r_total * r_total - r_inner * r_inner);
        let a_trabecular = PI * r_inner * r_inner;
        let a_total = PI * r_total * r_total;
        (self.cortical.elastic_modulus() * a_cortical
            + self.trabecular.elastic_modulus_trabecular() * a_trabecular)
            / a_total
    }
}
/// Huxley (1957) crossbridge model for skeletal muscle contraction.
///
/// Describes active force generation through attachment/detachment kinetics of
/// myosin crossbridges along the actin filament.  The model tracks the fraction
/// `n(x)` of attached crossbridges at extension `x`.
///
/// Rate constants:
/// - Attachment: `f(x)` = `f1` if 0 < x < h, 0 otherwise.
/// - Detachment: `g(x)` = `g1` if 0 < x < h, `g2` if x < 0.
#[derive(Debug, Clone)]
pub struct HuxleyModel {
    /// Stiffness of a single crossbridge \[N/m\].
    pub crossbridge_stiffness: f64,
    /// Maximum stroke distance h \[m\] — range of positive attachment.
    pub h: f64,
    /// Attachment rate constant f1 \[s⁻¹\].
    pub f1: f64,
    /// Detachment rate in positive region g1 \[s⁻¹\].
    pub g1: f64,
    /// Detachment rate in negative region g2 \[s⁻¹\].
    pub g2: f64,
    /// Number density of crossbridges per unit length \[m⁻¹\].
    pub crossbridge_density: f64,
}
impl HuxleyModel {
    /// Construct a new Huxley model with explicit parameters.
    pub fn new(
        crossbridge_stiffness: f64,
        h: f64,
        f1: f64,
        g1: f64,
        g2: f64,
        crossbridge_density: f64,
    ) -> Self {
        Self {
            crossbridge_stiffness,
            h,
            f1,
            g1,
            g2,
            crossbridge_density,
        }
    }
    /// Default fast-twitch skeletal muscle parameters.
    pub fn fast_twitch() -> Self {
        Self::new(2.0e-3, 11e-9, 3.68e9, 1.0e9, 3.68e9, 5.0e14)
    }
    /// Steady-state attached fraction `n_ss(x)` at extension `x` \[m\].
    ///
    /// Derived analytically from `(f+g)*n = f` in the attachment zone.
    pub fn steady_state_fraction(&self, x: f64) -> f64 {
        if x > 0.0 && x < self.h {
            self.f1 / (self.f1 + self.g1)
        } else {
            0.0
        }
    }
    /// Isometric (zero velocity) force per unit cross-sectional area \[Pa\].
    ///
    /// Integrates `n_ss(x) * k * x` over the crossbridge extension range.
    pub fn isometric_force(&self) -> f64 {
        let n_ss = self.f1 / (self.f1 + self.g1);
        self.crossbridge_density * n_ss * self.crossbridge_stiffness * self.h * self.h / 2.0
    }
    /// Force–velocity relationship using the Huxley model at sliding velocity `v` \[m/s\].
    ///
    /// Positive `v` = shortening (crossbridges move from positive to negative extension).
    /// Returns force relative to isometric force (dimensionless).
    pub fn force_velocity_ratio(&self, v: f64) -> f64 {
        if v.abs() < 1e-30 {
            return 1.0;
        }
        let f1 = self.f1;
        let g1 = self.g1;
        let g2 = self.g2;
        let h = self.h;
        if v > 0.0 {
            let phi = v / (h * (f1 + g1));
            let numerator = 1.0 - phi * (g2 / (f1 + g2));
            numerator.max(0.0)
        } else {
            let phi = (-v) / (h * (f1 + g1));
            (1.0 + phi * g2 / (f1 + g2)).min(1.8)
        }
    }
    /// Average power output at shortening velocity `v` \[W per unit area, Pa·m/s\].
    pub fn power_output(&self, v: f64) -> f64 {
        let f0 = self.isometric_force();
        f0 * self.force_velocity_ratio(v) * v
    }
    /// Optimal shortening velocity for maximum power \[m/s\].
    ///
    /// Uses a numerical golden-section search over \[0, v_max\].
    pub fn optimal_velocity(&self, v_max: f64) -> f64 {
        let mut a = 0.0f64;
        let mut b = v_max;
        let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
        let tol = 1e-9;
        let mut c = b - phi * (b - a);
        let mut d = a + phi * (b - a);
        while (b - a).abs() > tol {
            if self.power_output(c) < self.power_output(d) {
                a = c;
            } else {
                b = d;
            }
            c = b - phi * (b - a);
            d = a + phi * (b - a);
        }
        (a + b) / 2.0
    }
}
