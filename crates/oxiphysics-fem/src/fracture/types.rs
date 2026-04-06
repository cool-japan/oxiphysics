//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use std::f64::consts::PI;

/// Common geometry correction factors for standard crack configurations.
#[allow(dead_code)]
pub struct GeometryFactors;
#[allow(dead_code)]
impl GeometryFactors {
    /// Central through-crack in an infinite plate: F = 1.0.
    pub fn infinite_plate() -> f64 {
        1.0
    }
    /// Edge crack in a semi-infinite plate: F = 1.12.
    pub fn edge_crack() -> f64 {
        1.12
    }
    /// Penny-shaped (circular) internal crack: F = 2/pi.
    pub fn penny_crack() -> f64 {
        2.0 / PI
    }
    /// Central crack in a finite-width plate: F(a/W) polynomial approximation.
    ///
    /// F = (1 - 0.025 (a/W)^2 + 0.06 (a/W)^4) * sqrt(sec(pi*a/(2W)))
    pub fn finite_width(a: f64, w: f64) -> f64 {
        let ratio = a / w;
        let poly = 1.0 - 0.025 * ratio * ratio + 0.06 * ratio.powi(4);
        let sec_term = (PI * a / (2.0 * w)).cos();
        if sec_term.abs() < 1e-30 {
            return f64::INFINITY;
        }
        poly * (1.0 / sec_term).sqrt()
    }
}
/// A single point on the J-integral contour path.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct JIntegralPoint {
    /// Stress tensor components \[sigma_xx, sigma_yy, tau_xy\] at this point.
    pub stress: [f64; 3],
    /// Displacement gradient du_x/dx, du_x/dy, du_y/dx, du_y/dy.
    pub displacement_gradient: [f64; 4],
    /// Outward unit normal to the contour \[n_x, n_y\].
    pub normal: [f64; 2],
    /// Arc-length increment ds for this segment.
    pub ds: f64,
}
/// Stress intensity factor mode classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SifMode {
    /// Opening mode -- tensile stress perpendicular to crack plane.
    ModeI,
    /// Sliding mode -- in-plane shear parallel to crack.
    ModeII,
    /// Tearing mode -- out-of-plane shear (anti-plane).
    ModeIII,
}
/// Plane-strain fracture toughness test configuration (compact tension, CT).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct CompactTensionSpec {
    /// Width W (m) — from load line to back face.
    pub width: f64,
    /// Crack length a (m).
    pub crack_length: f64,
    /// Specimen thickness B (m).
    pub thickness: f64,
}
#[allow(dead_code)]
impl CompactTensionSpec {
    /// Create a compact tension specimen.
    pub fn new(width: f64, crack_length: f64, thickness: f64) -> Self {
        Self {
            width,
            crack_length,
            thickness,
        }
    }
    /// Normalised crack length ratio a/W.
    pub fn a_over_w(&self) -> f64 {
        self.crack_length / self.width
    }
    /// Geometry factor f(a/W) for the compact tension specimen.
    ///
    /// ASTM E399 polynomial fit:
    /// f = (2 + a/W) / (1 - a/W)^{3/2} * (0.886 + 4.64*a/W - 13.32*(a/W)^2 + 14.72*(a/W)^3 - 5.6*(a/W)^4)
    pub fn geometry_factor(&self) -> f64 {
        let x = self.a_over_w();
        let poly = 0.886 + 4.64 * x - 13.32 * x * x + 14.72 * x.powi(3) - 5.6 * x.powi(4);
        let factor = (2.0 + x) / (1.0 - x).powf(1.5);
        factor * poly
    }
    /// Compute K_I from applied load P (N).
    ///
    /// K_I = P / (B * sqrt(W)) * f(a/W)
    pub fn k_from_load(&self, load: f64) -> f64 {
        let f = self.geometry_factor();
        load / (self.thickness * self.width.sqrt()) * f
    }
    /// Check plane-strain validity: B >= 2.5 * (K_Ic / sigma_y)^2
    pub fn is_plane_strain_valid(&self, k_ic: f64, yield_stress: f64) -> bool {
        let required = 2.5 * (k_ic / yield_stress).powi(2);
        self.thickness >= required && self.crack_length >= required
    }
}
/// A quadrature point for the domain J-integral (equivalent domain integral).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DomainPoint {
    /// Stress tensor \[sigma_xx, sigma_yy, tau_xy\].
    pub stress: [f64; 3],
    /// Strain tensor \[eps_xx, eps_yy, gamma_xy\].
    pub strain: [f64; 3],
    /// Displacement gradient \[du_x/dx, du_x/dy, du_y/dx, du_y/dy\].
    pub displacement_gradient: [f64; 4],
    /// Weight function q at this point (q = 1 near tip, 0 at outer boundary).
    pub q: f64,
    /// q gradient \[dq/dx, dq/dy\].
    pub dq: [f64; 2],
    /// Area weight dA for this quadrature point.
    pub area: f64,
}
/// Bilinear traction-separation law evaluator.
#[allow(dead_code)]
pub struct CzmBilinear;
impl CzmBilinear {
    /// Normal traction for a given normal opening `delta_n`.
    ///
    /// Returns `t_n_max * (delta_n / delta_n_c)` for `delta_n <= delta_n_c / 2`
    /// (loading branch) and linearly decreasing to zero at `delta_n_c`.
    /// Returns 0 for `delta_n >= delta_n_c`.
    pub fn normal_traction(p: &CzmBilinearParams, delta_n: f64) -> f64 {
        if delta_n <= 0.0 {
            return 0.0;
        }
        let delta_n_0 = p.delta_n_c / 2.0;
        if delta_n <= delta_n_0 {
            p.t_n_max * delta_n / delta_n_0
        } else if delta_n < p.delta_n_c {
            p.t_n_max * (p.delta_n_c - delta_n) / (p.delta_n_c - delta_n_0)
        } else {
            0.0
        }
    }
    /// Shear traction for a given sliding displacement `delta_s`.
    pub fn shear_traction(p: &CzmBilinearParams, delta_s: f64) -> f64 {
        let abs_s = delta_s.abs();
        let delta_s_0 = p.delta_s_c / 2.0;
        let t = if abs_s <= delta_s_0 {
            p.t_s_max * abs_s / delta_s_0
        } else if abs_s < p.delta_s_c {
            p.t_s_max * (p.delta_s_c - abs_s) / (p.delta_s_c - delta_s_0)
        } else {
            0.0
        };
        if delta_s >= 0.0 { t } else { -t }
    }
    /// Mixed-mode effective opening displacement.
    ///
    /// `delta_eff = sqrt(delta_n^2 + delta_s^2)`
    pub fn effective_opening(delta_n: f64, delta_s: f64) -> f64 {
        (delta_n * delta_n + delta_s * delta_s).sqrt()
    }
    /// Mixed-mode criterion (BK criterion).
    ///
    /// Returns the mixed-mode fracture energy
    /// `G_c = G_I + (G_II - G_I) * (G_II / (G_I + G_II))^eta`.
    pub fn bk_fracture_energy(g_c_i: f64, g_c_ii: f64, g_i: f64, g_ii: f64, eta: f64) -> f64 {
        let g_total = g_i + g_ii;
        if g_total == 0.0 {
            return g_c_i;
        }
        let beta = g_ii / g_total;
        g_c_i + (g_c_ii - g_c_i) * beta.powf(eta)
    }
    /// Damage variable `D` for a given effective opening and history variable.
    ///
    /// `D = (delta_eff_max - delta_0) / (delta_f - delta_0)` clamped to \[0, 1\].
    pub fn damage_variable(delta_eff_max: f64, delta_0: f64, delta_f: f64) -> f64 {
        if delta_eff_max <= delta_0 {
            return 0.0;
        }
        ((delta_eff_max - delta_0) / (delta_f - delta_0)).min(1.0)
    }
    /// Secant stiffness `K_n * (1 - D)` for the normal direction.
    pub fn secant_normal_stiffness(k_n: f64, damage: f64) -> f64 {
        k_n * (1.0 - damage)
    }
    /// Xu-Needleman exponential normal traction.
    ///
    /// `t_n = (G_c_I / delta_n_c) * (delta_n / delta_n_c) * exp(-delta_n / delta_n_c)`
    pub fn xu_needleman_normal(g_c_i: f64, delta_n_c: f64, delta_n: f64) -> f64 {
        if delta_n < 0.0 {
            return 0.0;
        }
        let zeta = delta_n / delta_n_c;
        (g_c_i / delta_n_c) * zeta * (-zeta).exp()
    }
    /// Crack opening work per unit area up to opening `delta_n` (trapezoidal integration).
    #[allow(clippy::too_many_arguments)]
    pub fn crack_opening_work(p: &CzmBilinearParams, steps: usize) -> f64 {
        let h = p.delta_n_c / steps as f64;
        let mut work = 0.0;
        for i in 0..steps {
            let d0 = i as f64 * h;
            let d1 = d0 + h;
            let t0 = Self::normal_traction(p, d0);
            let t1 = Self::normal_traction(p, d1);
            work += 0.5 * (t0 + t1) * h;
        }
        work
    }
}
/// Williams crack-tip stress expansion (higher-order terms).
///
/// The full asymptotic expansion of the crack-tip stress field is:
///
/// ```text
/// sigma_ij = sum_{n=1}^{N} A_n * r^{(n/2 - 1)} * f_ij^(n)(theta)
/// ```
///
/// This struct computes contributions up to order N for Mode-I.
#[allow(dead_code)]
pub struct WilliamsExpansion {
    /// Stress intensity factor K_I (Pa sqrt(m)).
    pub k_i: f64,
    /// T-stress (non-singular term, Pa).
    pub t_stress: f64,
    /// Third-order Williams coefficient A_3 (Pa m^{-1/2}).
    pub a3: f64,
    /// Fourth-order Williams coefficient A_4 (Pa m^{-1}).
    pub a4: f64,
}
impl WilliamsExpansion {
    /// Create a new `WilliamsExpansion`.
    pub fn new(k_i: f64, t_stress: f64, a3: f64, a4: f64) -> Self {
        Self {
            k_i,
            t_stress,
            a3,
            a4,
        }
    }
    /// Compute higher-order terms of the Williams expansion at (r, theta).
    ///
    /// Returns `[sigma_xx, sigma_yy, tau_xy]` including terms up to order 4.
    ///
    /// The expansion follows:
    /// - Order 1 (r^{-1/2}): classical singular SIF term
    /// - Order 2 (r^0): T-stress
    /// - Order 3 (r^{1/2}): sub-singular term
    /// - Order 4 (r^1): linear term
    #[allow(clippy::too_many_arguments)]
    pub fn compute_higher_order_terms(&self, r: f64, theta: f64) -> [f64; 3] {
        if r < 1e-20 {
            return [f64::INFINITY, f64::INFINITY, 0.0];
        }
        let half = theta / 2.0;
        let three_half = 3.0 * theta / 2.0;
        let s1 = self.k_i / (2.0 * PI * r).sqrt();
        let cos_h = half.cos();
        let sin_h = half.sin();
        let cos_3h = three_half.cos();
        let sin_3h = three_half.sin();
        let s1_xx = s1 * cos_h * (1.0 - sin_h * sin_3h);
        let s1_yy = s1 * cos_h * (1.0 + sin_h * sin_3h);
        let s1_xy = s1 * sin_h * cos_h * cos_3h;
        let s2_xx = self.t_stress;
        let s2_yy = 0.0;
        let s2_xy = 0.0;
        let r_sqrt = r.sqrt();
        let s3_xx = self.a3 * r_sqrt * (cos_h * (1.0 - 2.0 * sin_h * sin_h));
        let s3_yy = self.a3 * r_sqrt * (cos_h * (1.0 + 2.0 * sin_h * sin_h));
        let s3_xy = self.a3 * r_sqrt * sin_h * (2.0 * cos_h * cos_h - 1.0);
        let x = r * theta.cos();
        let y = r * theta.sin();
        let s4_xx = self.a4 * (x - y);
        let s4_yy = self.a4 * (x + y);
        let s4_xy = self.a4 * y;
        [
            s1_xx + s2_xx + s3_xx + s4_xx,
            s1_yy + s2_yy + s3_yy + s4_yy,
            s1_xy + s2_xy + s3_xy + s4_xy,
        ]
    }
    /// Extract Mode-I SIF from the singular term at given (r, theta=0).
    ///
    /// From the Williams expansion: K_I = sigma_yy * sqrt(2 pi r) at theta=0.
    pub fn extract_ki_from_stress(&self, r: f64, sigma_yy_at_theta0: f64) -> f64 {
        sigma_yy_at_theta0 * (2.0 * PI * r).sqrt()
    }
}
/// Linear elastic fracture mechanics solver.
///
/// All formulae follow Broek, *Elementary Engineering Fracture Mechanics* (4th ed.)
/// and Anderson, *Fracture Mechanics* (3rd ed.).
pub struct LinearFracture;
impl LinearFracture {
    /// Compute the Mode-I stress intensity factor K_I.
    ///
    /// ```text
    /// K_I = sigma_inf * sqrt(pi a) * F
    /// ```
    pub fn stress_intensity_ki(sigma_far: f64, a: f64, geometry_factor: f64) -> f64 {
        sigma_far * (PI * a).sqrt() * geometry_factor
    }
    /// Compute the Mode-II stress intensity factor K_II.
    ///
    /// ```text
    /// K_II = tau * sqrt(pi a) * F
    /// ```
    pub fn stress_intensity_kii(tau: f64, a: f64, geometry_factor: f64) -> f64 {
        tau * (PI * a).sqrt() * geometry_factor
    }
    /// Compute crack-tip stresses using the Williams asymptotic expansion.
    ///
    /// Returns `[sigma_xx, sigma_yy, tau_xy]` in Pa.
    pub fn crack_tip_stress(r: f64, theta: f64, ki: f64, kii: f64) -> [f64; 3] {
        let half_theta = theta / 2.0;
        let prefactor = 1.0 / (2.0 * PI * r).sqrt();
        let cos_h = half_theta.cos();
        let sin_h = half_theta.sin();
        let cos_3h = (3.0 * half_theta).cos();
        let sin_3h = (3.0 * half_theta).sin();
        let sigma_xx = prefactor
            * (ki * cos_h * (1.0 - sin_h * sin_3h) - kii * sin_h * (2.0 + cos_h * cos_3h));
        let sigma_yy =
            prefactor * (ki * cos_h * (1.0 + sin_h * sin_3h) + kii * sin_h * cos_h * cos_3h);
        let tau_xy =
            prefactor * (ki * sin_h * cos_h * cos_3h + kii * cos_h * (1.0 - sin_h * sin_3h));
        [sigma_xx, sigma_yy, tau_xy]
    }
    /// Compute the Griffith critical (fracture) stress for a given crack size.
    pub fn griffith_critical_stress(a: f64, fracture_toughness: f64, geometry_factor: f64) -> f64 {
        fracture_toughness / ((PI * a).sqrt() * geometry_factor)
    }
    /// Estimate the J-integral from stress intensity factors (plane stress).
    ///
    /// J = (K_I^2 + K_II^2) / E*
    pub fn j_integral_estimate(ki: f64, kii: f64, e_star: f64) -> f64 {
        (ki * ki + kii * kii) / e_star
    }
    /// Compute fatigue crack growth rate using the Paris-Erdogan law.
    ///
    /// da/dN = C * (delta_K)^m
    pub fn paris_law_da_dn(delta_k: f64, c: f64, m: f64) -> f64 {
        c * delta_k.powf(m)
    }
    /// Integrate the Paris law to estimate fatigue life in cycles.
    #[allow(clippy::too_many_arguments)]
    pub fn fatigue_life_cycles(
        a_initial: f64,
        a_final: f64,
        delta_sigma: f64,
        geometry_factor: f64,
        c: f64,
        m: f64,
    ) -> f64 {
        let dk_coeff = (delta_sigma * geometry_factor * PI.sqrt()).powf(m);
        let denominator_common = c * dk_coeff;
        let exponent = 1.0 - m / 2.0;
        if exponent.abs() < 1e-12 {
            (a_final / a_initial).ln() / denominator_common
        } else {
            (a_final.powf(exponent) - a_initial.powf(exponent)) / (exponent * denominator_common)
        }
    }
}
/// Mixed-mode SIF extraction via the interaction integral method.
///
/// The interaction integral is:
///
/// ```text
/// I^{(1,2)} = integral_Omega (sigma_ij^{(1)} du_i^{(2)}/dx1
///                            + sigma_ij^{(2)} du_i^{(1)}/dx1
///                            - W^{(1,2)} delta_{1j}) dq/dx_j dV
/// ```
///
/// Using auxiliary fields (superscript 2 = analytical Williams field),
/// the individual SIFs are recovered via:
///
/// K_I = (E* / 2) * I^{(1,a_I)}   (auxiliary = Mode-I unit field)
/// K_II = (E* / 2) * I^{(1,a_II)} (auxiliary = Mode-II unit field)
#[allow(dead_code)]
pub struct SifMixed {
    /// Young's modulus (Pa).
    pub e: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Plane stress (false) or plane strain (true).
    pub plane_strain: bool,
}
impl SifMixed {
    /// Create a new `SifMixed` solver.
    pub fn new(e: f64, nu: f64, plane_strain: bool) -> Self {
        Self {
            e,
            nu,
            plane_strain,
        }
    }
    /// Effective modulus E* (plane stress or plane strain).
    fn e_star(&self) -> f64 {
        if self.plane_strain {
            self.e / (1.0 - self.nu * self.nu)
        } else {
            self.e
        }
    }
    /// Compute the interaction integral I^{(1,2)}.
    ///
    /// Returns the scalar interaction integral value.
    pub fn compute_interaction_integral(&self, points: &[InteractionIntegralPoint]) -> f64 {
        let mut integral = 0.0;
        for pt in points {
            let s1 = pt.stress_fem;
            let s2 = pt.stress_aux;
            let du1 = pt.grad_u_fem;
            let du2 = pt.grad_u_aux;
            let _e1 = pt.strain_fem;
            let e2 = pt.strain_aux;
            let w12 = s1[0] * e2[0] + s1[1] * e2[1] + 2.0 * s1[2] * e2[2];
            let sw1 = s1[0] * du2[0] + s1[2] * du2[2] + s2[0] * du1[0] + s2[2] * du1[2];
            let sw2 = s1[2] * du2[0] + s1[1] * du2[2] + s2[2] * du1[0] + s2[1] * du1[2];
            let dq1 = pt.grad_q[0];
            let dq2 = pt.grad_q[1];
            let integrand = (sw1 - w12) * dq1 + sw2 * dq2;
            integral += integrand * pt.dv;
        }
        integral
    }
    /// Extract K_I from the interaction integral with auxiliary Mode-I field.
    ///
    /// K_I = (E* / 2) * I^{(1, aux_I)}  where K_I_aux = 1 Pa sqrt(m)
    pub fn extract_ki(&self, interaction_integral: f64) -> f64 {
        self.e_star() / 2.0 * interaction_integral
    }
    /// Extract K_II from the interaction integral with auxiliary Mode-II field.
    ///
    /// K_II = (E* / 2) * I^{(1, aux_II)}  where K_II_aux = 1 Pa sqrt(m)
    pub fn extract_kii(&self, interaction_integral: f64) -> f64 {
        self.e_star() / 2.0 * interaction_integral
    }
    /// Compute the effective K_eff = sqrt(K_I^2 + K_II^2).
    pub fn effective_sif(&self, ki: f64, kii: f64) -> f64 {
        (ki * ki + kii * kii).sqrt()
    }
    /// Check fracture using the mixed-mode criterion:
    /// K_I^2 + K_II^2 >= K_Ic^2
    pub fn is_fracture(&self, ki: f64, kii: f64, k_ic: f64) -> bool {
        ki * ki + kii * kii >= k_ic * k_ic
    }
}
/// Parameters for a mixed-mode bilinear cohesive zone model.
#[derive(Debug, Clone)]
pub struct CzmBilinearParams {
    /// Peak normal traction `t_n_max` (Pa).
    pub t_n_max: f64,
    /// Peak shear traction `t_s_max` (Pa).
    pub t_s_max: f64,
    /// Normal critical opening displacement `delta_n_c` (m).
    pub delta_n_c: f64,
    /// Shear critical sliding displacement `delta_s_c` (m).
    pub delta_s_c: f64,
    /// Mode-I fracture energy `G_c_I` (J/m²).
    pub g_c_i: f64,
    /// Mode-II fracture energy `G_c_II` (J/m²).
    pub g_c_ii: f64,
}
impl CzmBilinearParams {
    /// Create parameters from peak tractions and critical displacements.
    /// Fracture energies are set to `0.5 * t_max * delta_c`.
    pub fn from_peak_and_critical(
        t_n_max: f64,
        t_s_max: f64,
        delta_n_c: f64,
        delta_s_c: f64,
    ) -> Self {
        Self {
            t_n_max,
            t_s_max,
            delta_n_c,
            delta_s_c,
            g_c_i: 0.5 * t_n_max * delta_n_c,
            g_c_ii: 0.5 * t_s_max * delta_s_c,
        }
    }
}
/// A single integration point for the domain integral.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DomainIntegralPoint {
    /// Stress tensor \[sxx, syy, sxy\] (Pa).
    pub stress: [f64; 3],
    /// Strain tensor \[exx, eyy, exy\] (consistent with stress).
    pub strain: [f64; 3],
    /// Displacement gradient: \[du1_dx1, du1_dx2, du2_dx1, du2_dx2\].
    pub grad_u: [f64; 4],
    /// Weight function gradient: \[dq_dx1, dq_dx2\].
    pub grad_q: [f64; 2],
    /// Jacobian-weighted volume element (volume weight for this integration point).
    pub dv: f64,
}
/// Three-point bending fracture toughness specimen (SENB).
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ThreePointBendSpec {
    /// Span S (m).
    pub span: f64,
    /// Beam width W (m).
    pub width: f64,
    /// Crack length a (m).
    pub crack_length: f64,
    /// Thickness B (m).
    pub thickness: f64,
}
#[allow(dead_code)]
impl ThreePointBendSpec {
    /// Create a three-point bending specimen.
    pub fn new(span: f64, width: f64, crack_length: f64, thickness: f64) -> Self {
        Self {
            span,
            width,
            crack_length,
            thickness,
        }
    }
    /// Geometry factor for SENB (ASTM E399, simplified).
    ///
    /// f(a/W) = 3*(a/W)^{0.5} * \[1.99 - (a/W)*(1-a/W)*(2.15 - 3.93*a/W + 2.7*(a/W)^2)\]
    ///          / (2*(1+2*a/W)*(1-a/W)^{3/2})
    pub fn geometry_factor(&self) -> f64 {
        let x = self.crack_length / self.width;
        let num = 3.0 * x.sqrt() * (1.99 - x * (1.0 - x) * (2.15 - 3.93 * x + 2.7 * x * x));
        let denom = 2.0 * (1.0 + 2.0 * x) * (1.0 - x).powf(1.5);
        num / denom.max(f64::EPSILON)
    }
    /// K_I from applied load P (N).
    pub fn k_from_load(&self, load: f64) -> f64 {
        let f = self.geometry_factor();
        (load * self.span) / (self.thickness * self.width * self.width) * f
    }
}
/// Mixed-mode fracture criterion type.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MixedModeCriterion {
    /// Maximum Energy Release Rate (MERR).
    MaxEnergyReleaseRate,
    /// Minimum Strain Energy Density (SED, Sih criterion).
    MinStrainEnergyDensity,
    /// Maximum Tangential Stress (MTS, Erdogan-Sih criterion).
    MaxTangentialStress,
}
/// Traction-separation law type for cohesive zone modelling.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum CzmLaw {
    /// Bilinear (linear softening).
    Bilinear,
    /// Exponential (Xu-Needleman).
    Exponential,
    /// Trapezoidal.
    Trapezoidal,
    /// PPR (Park-Paulino-Roesler).
    Ppr,
}
/// A single integration point for the interaction integral.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct InteractionIntegralPoint {
    /// Stress from the actual (FEM) field \[sxx, syy, sxy\] (Pa).
    pub stress_fem: [f64; 3],
    /// Displacement gradient from FEM: \[du1/dx1, du1/dx2, du2/dx1, du2/dx2\].
    pub grad_u_fem: [f64; 4],
    /// Stress from the auxiliary (analytical) field \[sxx, syy, sxy\] (Pa).
    pub stress_aux: [f64; 3],
    /// Displacement gradient from auxiliary field: \[du1/dx1, du1/dx2, du2/dx1, du2/dx2\].
    pub grad_u_aux: [f64; 4],
    /// Strain from FEM field \[exx, eyy, exy\].
    pub strain_fem: [f64; 3],
    /// Strain from auxiliary field \[exx, eyy, exy\].
    pub strain_aux: [f64; 3],
    /// Weight function gradient: \[dq/dx1, dq/dx2\].
    pub grad_q: [f64; 2],
    /// Volume weight for this integration point.
    pub dv: f64,
}
/// Domain integral (EDI) computation of the J-integral.
///
/// Uses the equivalent domain integral which converts the contour integral
/// to a volume integral over a domain Omega surrounding the crack tip.
///
/// J = integral_Omega (sigma_ij du_i/dx1 - W delta_1j) dq/dx_j dV
///
/// where q is a smooth weight function (=1 near crack tip, =0 far away).
#[allow(dead_code)]
pub struct JIntegral {
    /// Young's modulus (Pa).
    pub e: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Plane stress (false) or plane strain (true).
    pub plane_strain: bool,
}
impl JIntegral {
    /// Create a new `JIntegral` solver.
    pub fn new(e: f64, nu: f64, plane_strain: bool) -> Self {
        Self {
            e,
            nu,
            plane_strain,
        }
    }
    /// Compute the domain (equivalent domain) integral form of J.
    ///
    /// For each integration point the integrand is:
    /// `(sigma_i1 * du_i/dx_1 - W * delta_1j) * dq/dx_j`
    ///
    /// where W = 0.5 * sigma:epsilon is the strain energy density.
    pub fn compute_domain_integral(&self, points: &[DomainIntegralPoint]) -> f64 {
        let mut j = 0.0;
        for pt in points {
            let sxx = pt.stress[0];
            let syy = pt.stress[1];
            let sxy = pt.stress[2];
            let exx = pt.strain[0];
            let eyy = pt.strain[1];
            let exy = pt.strain[2];
            let w = 0.5 * (sxx * exx + syy * eyy + 2.0 * sxy * exy);
            let du1_dx1 = pt.grad_u[0];
            let du2_dx1 = pt.grad_u[2];
            let dq_dx1 = pt.grad_q[0];
            let dq_dx2 = pt.grad_q[1];
            let sigma_du_dx1 = sxx * du1_dx1 + sxy * du2_dx1;
            let sigma_du_dx2_part = sxy * du1_dx1 + syy * du2_dx1;
            let integrand = (sigma_du_dx1 - w) * dq_dx1 + sigma_du_dx2_part * dq_dx2;
            j += integrand * pt.dv;
        }
        j
    }
    /// Estimate J from stress intensity factors (plane stress or plane strain).
    pub fn j_from_sif(&self, ki: f64, kii: f64) -> f64 {
        let e_eff = if self.plane_strain {
            self.e / (1.0 - self.nu * self.nu)
        } else {
            self.e
        };
        (ki * ki + kii * kii) / e_eff
    }
    /// Compute Mode-I SIF from J (plane stress or plane strain).
    pub fn ki_from_j(&self, j: f64) -> f64 {
        let e_eff = if self.plane_strain {
            self.e / (1.0 - self.nu * self.nu)
        } else {
            self.e
        };
        (j * e_eff).sqrt()
    }
}
