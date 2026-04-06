// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Geomaterial constitutive models.
//!
//! Implements advanced constitutive models for soils, rocks, and special
//! geomaterials:
//!
//! - [`ModifiedCamClay`]: critical-state elasto-plastic model (Roscoe & Burland).
//! - [`MohrCoulombModel`]: linear Mohr-Coulomb failure with tension cut-off.
//! - [`CapPlasticityModel`]: Drucker-Prager / Cap model for soil compaction.
//! - [`HoekBrownRock`]: empirical strength criterion for jointed rock masses.
//! - [`KozenyCarmanPermeability`]: permeability from void ratio and particle size.
//! - [`TerzaghiConsolidation`]: 1-D consolidation (pore-pressure dissipation).
//! - [`SwellingClayModel`]: expansive clay with suction-dependent swelling strain.
//! - [`LiquefactionCriteria`]: simplified cyclic-stress-ratio screening.
//! - [`RockfillModel`]: nonlinear strength and stiffness for coarse fill.
//! - [`FrozenSoilModel`]: temperature-dependent strength and creep for frozen ground.
//!
//! References
//! ----------
//! Roscoe, K. H. & Burland, J. B. (1968) *On the generalised stress–strain
//! behaviour of "wet" clay*, Cambridge.
//! Hoek, E. & Brown, E. T. (1997) *Int. J. Rock Mech.* 34(8), 1165–1186.
//! Kozeny, J. (1927); Carman, P. C. (1937,1956).
//! Terzaghi, K. (1943) *Theoretical Soil Mechanics*, Wiley.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration \[m s^-2\].
const G_ACC: f64 = 9.81;

/// Unit weight of water \[kN m^-3\].
const GAMMA_W: f64 = 9.81;

/// Absolute zero offset for frozen soil temperature \[C\].
const ABS_ZERO_C: f64 = -273.15;

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

/// Clamp a value to \[lo, hi\].
#[inline]
fn clamp_f64(v: f64, lo: f64, hi: f64) -> f64 {
    if v < lo {
        lo
    } else if v > hi {
        hi
    } else {
        v
    }
}

/// Mean effective stress from principal stresses.
#[inline]
fn mean_stress(s1: f64, s2: f64, s3: f64) -> f64 {
    (s1 + s2 + s3) / 3.0
}

/// Deviatoric stress invariant q from principal stresses.
#[inline]
fn deviatoric_q(s1: f64, s2: f64, s3: f64) -> f64 {
    let t1 = (s1 - s2).powi(2);
    let t2 = (s2 - s3).powi(2);
    let t3 = (s3 - s1).powi(2);
    ((t1 + t2 + t3) / 2.0).sqrt()
}

// ===========================================================================
// Modified Cam-Clay
// ===========================================================================

/// Modified Cam-Clay critical-state model.
///
/// The yield surface in (p', q) space is:
///
///   q^2 / M^2 + p' (p' - p_c) = 0
///
/// where M is the critical-state stress ratio and p_c is the
/// pre-consolidation pressure (hardening variable).
#[derive(Debug, Clone)]
pub struct ModifiedCamClay {
    /// Slope of critical state line M.
    pub m_slope: f64,
    /// Slope of normal compression line lambda.
    pub lambda: f64,
    /// Slope of swelling/recompression line kappa.
    pub kappa: f64,
    /// Initial pre-consolidation pressure p_c0 \[Pa\].
    pub pc0: f64,
    /// Current pre-consolidation pressure p_c \[Pa\].
    pub pc: f64,
    /// Initial specific volume v0.
    pub v0: f64,
    /// Poisson's ratio for elastic unloading.
    pub poisson: f64,
}

impl ModifiedCamClay {
    /// Create a new Modified Cam-Clay model.
    pub fn new(m_slope: f64, lambda: f64, kappa: f64, pc0: f64, v0: f64, poisson: f64) -> Self {
        Self {
            m_slope,
            lambda,
            kappa,
            pc0,
            pc: pc0,
            v0,
            poisson,
        }
    }

    /// Evaluate the yield function f(p', q).
    ///
    /// Returns < 0 inside the yield surface (elastic), 0 on the surface,
    /// > 0 outside (inadmissible before return mapping).
    pub fn yield_function(&self, p: f64, q: f64) -> f64 {
        q * q / (self.m_slope * self.m_slope) + p * (p - self.pc)
    }

    /// Check whether the stress state is on or outside the yield surface.
    pub fn is_yielding(&self, p: f64, q: f64) -> bool {
        self.yield_function(p, q) >= -1e-12
    }

    /// Elastic bulk modulus K from current mean stress.
    pub fn bulk_modulus(&self, p: f64) -> f64 {
        self.v0 * p / self.kappa
    }

    /// Elastic shear modulus G (derived from K and nu).
    pub fn shear_modulus(&self, p: f64) -> f64 {
        let k = self.bulk_modulus(p);
        3.0 * k * (1.0 - 2.0 * self.poisson) / (2.0 * (1.0 + self.poisson))
    }

    /// Harden: update pre-consolidation pressure from plastic volumetric
    /// strain increment d_eps_v_p (positive = compression).
    pub fn harden(&mut self, d_eps_v_p: f64) {
        let factor = self.v0 / (self.lambda - self.kappa);
        self.pc *= (factor * d_eps_v_p).exp();
    }

    /// Stress ratio eta = q / p'.
    pub fn stress_ratio(p: f64, q: f64) -> f64 {
        if p.abs() < 1e-15 { 0.0 } else { q / p }
    }

    /// Over-consolidation ratio OCR = p_c / p'.
    pub fn ocr(&self, p: f64) -> f64 {
        if p.abs() < 1e-15 {
            f64::INFINITY
        } else {
            self.pc / p
        }
    }

    /// Volumetric strain on the normal compression line from p'_0 to p'.
    pub fn ncl_volumetric_strain(&self, p0: f64, p: f64) -> f64 {
        if p0 <= 0.0 || p <= 0.0 {
            return 0.0;
        }
        (self.lambda / self.v0) * (p / p0).ln()
    }

    /// Elastic (swelling-line) volumetric strain from p0 to p.
    pub fn elastic_volumetric_strain(&self, p0: f64, p: f64) -> f64 {
        if p0 <= 0.0 || p <= 0.0 {
            return 0.0;
        }
        (self.kappa / self.v0) * (p / p0).ln()
    }

    /// Size of the yield ellipse along the p-axis (= pc).
    pub fn yield_surface_size(&self) -> f64 {
        self.pc
    }

    /// Critical-state mean stress for given q.
    pub fn critical_state_p(&self, q: f64) -> f64 {
        q / self.m_slope
    }

    /// Reset to initial pre-consolidation.
    pub fn reset(&mut self) {
        self.pc = self.pc0;
    }
}

// ===========================================================================
// Mohr-Coulomb
// ===========================================================================

/// Linear Mohr-Coulomb failure criterion with optional tension cut-off.
///
/// tau_f = c + sigma_n * tan(phi)
#[derive(Debug, Clone)]
pub struct MohrCoulombModel {
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
    /// Friction angle phi \[radians\].
    pub friction_angle: f64,
    /// Dilation angle psi \[radians\].
    pub dilation_angle: f64,
    /// Tensile strength sigma_t \[Pa\] (positive value, used as cut-off).
    pub tensile_strength: f64,
}

impl MohrCoulombModel {
    /// Construct from cohesion \[Pa\], friction angle \[deg\], dilation angle \[deg\],
    /// and tensile strength \[Pa\].
    pub fn new(cohesion: f64, friction_deg: f64, dilation_deg: f64, tensile_strength: f64) -> Self {
        Self {
            cohesion,
            friction_angle: friction_deg.to_radians(),
            dilation_angle: dilation_deg.to_radians(),
            tensile_strength,
        }
    }

    /// Shear strength at given normal stress sigma_n (positive = compressive).
    pub fn shear_strength(&self, sigma_n: f64) -> f64 {
        self.cohesion + sigma_n * self.friction_angle.tan()
    }

    /// Yield function value: f = |tau| - c - sigma_n * tan(phi).
    /// Negative means elastic.
    pub fn yield_function(&self, sigma_n: f64, tau: f64) -> f64 {
        tau.abs() - self.shear_strength(sigma_n)
    }

    /// Check failure given principal stresses sigma_1 >= sigma_3.
    pub fn is_failed_principals(&self, sigma1: f64, sigma3: f64) -> bool {
        let sin_phi = self.friction_angle.sin();
        let lhs = (sigma1 - sigma3) / 2.0;
        let rhs = (sigma1 + sigma3) / 2.0 * sin_phi + self.cohesion * self.friction_angle.cos();
        lhs >= rhs - 1e-12
    }

    /// Maximum principal stress difference at failure (confined by sigma3).
    pub fn max_deviator(&self, sigma3: f64) -> f64 {
        let sin_phi = self.friction_angle.sin();
        let _cos_phi = self.friction_angle.cos();
        let n_phi = (1.0 + sin_phi) / (1.0 - sin_phi);
        2.0 * self.cohesion * n_phi.sqrt() + sigma3 * (n_phi - 1.0)
    }

    /// Check tension cut-off: fails if sigma3 < -tensile_strength.
    pub fn tension_cutoff(&self, sigma3: f64) -> bool {
        sigma3 < -self.tensile_strength
    }

    /// Plastic potential dilation factor d_eps_v / d_gamma.
    pub fn dilation_factor(&self) -> f64 {
        self.dilation_angle.sin()
    }

    /// Equivalent friction angle from cohesion and UCS.
    pub fn from_ucs(ucs: f64, cohesion: f64) -> f64 {
        // sigma_c = 2 c cos(phi) / (1 - sin(phi))
        // Solve iteratively with Newton.
        let mut phi = 30.0_f64.to_radians();
        for _ in 0..20 {
            let sp = phi.sin();
            let cp = phi.cos();
            let f_val = ucs * (1.0 - sp) - 2.0 * cohesion * cp;
            let df = -ucs * cp + 2.0 * cohesion * sp;
            if df.abs() < 1e-30 {
                break;
            }
            phi -= f_val / df;
            phi = clamp_f64(phi, 0.001, PI / 2.0 - 0.001);
        }
        phi
    }
}

// ===========================================================================
// Cap Plasticity
// ===========================================================================

/// Cap plasticity model (Drucker-Prager shear surface + elliptical cap).
///
/// Used for modelling compactive yielding in soils, concrete, and powders.
#[derive(Debug, Clone)]
pub struct CapPlasticityModel {
    /// Drucker-Prager friction parameter alpha.
    pub alpha: f64,
    /// Drucker-Prager cohesion parameter k.
    pub k_cohesion: f64,
    /// Cap eccentricity R (ratio of ellipse axes).
    pub cap_r: f64,
    /// Hardening parameter: initial cap position X_0 (hydrostatic yield).
    pub x0: f64,
    /// Current cap position X.
    pub x_cap: f64,
    /// Hardening modulus W (relates plastic vol strain to cap movement).
    pub hardening_w: f64,
    /// Cap shape parameter D.
    pub hardening_d: f64,
}

impl CapPlasticityModel {
    /// Construct a new cap plasticity model.
    pub fn new(
        alpha: f64,
        k_cohesion: f64,
        cap_r: f64,
        x0: f64,
        hardening_w: f64,
        hardening_d: f64,
    ) -> Self {
        Self {
            alpha,
            k_cohesion,
            cap_r,
            x0,
            x_cap: x0,
            hardening_w,
            hardening_d,
        }
    }

    /// Drucker-Prager shear yield function: f_s = q - alpha * p - k.
    pub fn shear_yield(&self, p: f64, q: f64) -> f64 {
        q - self.alpha * p - self.k_cohesion
    }

    /// Cap yield function (elliptical): f_c = sqrt((p - L)^2 + (Rq)^2) - (X - L).
    /// L = X - R * (alpha * X + k).
    pub fn cap_yield(&self, p: f64, q: f64) -> f64 {
        let l = self.cap_l();
        let rq = self.cap_r * q;
        let dp = p - l;
        (dp * dp + rq * rq).sqrt() - (self.x_cap - l)
    }

    /// Intersection of DP line and hydrostatic axis: L.
    fn cap_l(&self) -> f64 {
        self.x_cap - self.cap_r * (self.alpha * self.x_cap + self.k_cohesion)
    }

    /// Check if stress state is elastic.
    pub fn is_elastic(&self, p: f64, q: f64) -> bool {
        self.shear_yield(p, q) < 1e-12 && self.cap_yield(p, q) < 1e-12
    }

    /// Harden cap from plastic volumetric strain increment.
    pub fn harden_cap(&mut self, d_eps_v_p: f64) {
        // X moves according to exponential hardening
        let eps_total = self.current_plastic_vol_strain() + d_eps_v_p;
        self.x_cap = self.x0
            - (1.0 / self.hardening_d) * (1.0 - eps_total / self.hardening_w).ln().max(-50.0);
    }

    /// Current accumulated plastic volumetric strain from cap position.
    pub fn current_plastic_vol_strain(&self) -> f64 {
        self.hardening_w * (1.0 - (-self.hardening_d * (self.x_cap - self.x0)).exp())
    }

    /// Reset cap to initial position.
    pub fn reset(&mut self) {
        self.x_cap = self.x0;
    }
}

// ===========================================================================
// Hoek-Brown Rock Criterion
// ===========================================================================

/// Generalised Hoek-Brown rock-mass strength criterion.
///
/// sigma_1 = sigma_3 + sigma_ci * (m_b * sigma_3 / sigma_ci + s)^a
#[derive(Debug, Clone)]
pub struct HoekBrownRock {
    /// Intact uniaxial compressive strength sigma_ci \[Pa\].
    pub sigma_ci: f64,
    /// Material constant m_i (intact rock).
    pub m_i: f64,
    /// Geological Strength Index GSI \[0..100\].
    pub gsi: f64,
    /// Disturbance factor D \[0..1\].
    pub disturbance: f64,
    /// Derived: m_b.
    pub m_b: f64,
    /// Derived: s.
    pub s: f64,
    /// Derived: a.
    pub a: f64,
}

impl HoekBrownRock {
    /// Build from intact rock parameters and GSI.
    pub fn new(sigma_ci: f64, m_i: f64, gsi: f64, disturbance: f64) -> Self {
        let m_b = m_i * ((gsi - 100.0) / (28.0 - 14.0 * disturbance)).exp();
        let s = ((gsi - 100.0) / (9.0 - 3.0 * disturbance)).exp();
        let a = 0.5 + ((-gsi / 15.0_f64).exp() - (-20.0 / 3.0_f64).exp()) / 6.0;
        Self {
            sigma_ci,
            m_i,
            gsi,
            disturbance,
            m_b,
            s,
            a,
        }
    }

    /// Major principal stress at failure given confining stress sigma3 >= 0.
    pub fn sigma1_at_failure(&self, sigma3: f64) -> f64 {
        let ratio = sigma3 / self.sigma_ci;
        let inner = self.m_b * ratio + self.s;
        if inner <= 0.0 {
            sigma3
        } else {
            sigma3 + self.sigma_ci * inner.powf(self.a)
        }
    }

    /// Uniaxial compressive strength of the rock mass.
    pub fn ucs_rock_mass(&self) -> f64 {
        self.sigma_ci * self.s.powf(self.a)
    }

    /// Tensile strength estimate (sigma_3 where sigma_1 = sigma_3).
    pub fn tensile_strength(&self) -> f64 {
        -self.s * self.sigma_ci / self.m_b
    }

    /// Equivalent Mohr-Coulomb cohesion and friction angle for a given
    /// confining stress range \[sig3_min, sig3_max\].
    ///
    /// Returns (cohesion \[Pa\], friction_angle \[radians\]).
    pub fn equivalent_mc(&self, sig3_min: f64, sig3_max: f64) -> (f64, f64) {
        let n = 10;
        let ds = (sig3_max - sig3_min) / n as f64;
        let mut sum_tau = 0.0;
        let mut sum_sigma = 0.0;
        let mut sum_st = 0.0;
        let mut sum_s2 = 0.0;
        let nf = (n + 1) as f64;
        for i in 0..=n {
            let s3 = sig3_min + i as f64 * ds;
            let s1 = self.sigma1_at_failure(s3);
            let sigma_n = (s1 + s3) / 2.0;
            let tau = (s1 - s3) / 2.0;
            sum_sigma += sigma_n;
            sum_tau += tau;
            sum_st += sigma_n * tau;
            sum_s2 += sigma_n * sigma_n;
        }
        let denom = nf * sum_s2 - sum_sigma * sum_sigma;
        if denom.abs() < 1e-30 {
            return (0.0, 0.0);
        }
        let tan_phi = (nf * sum_st - sum_sigma * sum_tau) / denom;
        let c = (sum_tau - tan_phi * sum_sigma) / nf;
        (c.max(0.0), tan_phi.atan().max(0.0))
    }

    /// Deformation modulus E_rm \[Pa\] (Hoek-Diederichs 2006).
    pub fn deformation_modulus(&self) -> f64 {
        let e_i = self.sigma_ci * 1000.0; // rough estimate
        e_i * (0.02
            + (1.0 - self.disturbance / 2.0)
                / (1.0 + ((60.0 + 15.0 * self.disturbance - self.gsi) / 11.0).exp()))
    }
}

// ===========================================================================
// Kozeny-Carman Permeability
// ===========================================================================

/// Kozeny-Carman permeability model.
///
/// k = C * (e^3 / (1 + e)) * d10^2
#[derive(Debug, Clone)]
pub struct KozenyCarmanPermeability {
    /// Shape/tortuosity factor C.
    pub shape_factor: f64,
    /// Effective grain diameter d10 \[m\].
    pub d10: f64,
    /// Dynamic viscosity of pore fluid \[Pa s\].
    pub viscosity: f64,
}

impl KozenyCarmanPermeability {
    /// Create a Kozeny-Carman model with given shape factor, d10, viscosity.
    pub fn new(shape_factor: f64, d10: f64, viscosity: f64) -> Self {
        Self {
            shape_factor,
            d10,
            viscosity,
        }
    }

    /// Default model for water at 20 C with standard shape factor.
    pub fn default_water(d10: f64) -> Self {
        Self::new(1.0 / 180.0, d10, 1.002e-3)
    }

    /// Intrinsic permeability k \[m^2\] for void ratio e.
    pub fn intrinsic_permeability(&self, e: f64) -> f64 {
        if e <= 0.0 {
            return 0.0;
        }
        self.shape_factor * (e.powi(3) / (1.0 + e)) * self.d10.powi(2)
    }

    /// Hydraulic conductivity K \[m/s\] for void ratio e.
    pub fn hydraulic_conductivity(&self, e: f64) -> f64 {
        let k_intr = self.intrinsic_permeability(e);
        k_intr * GAMMA_W * 1000.0 / self.viscosity // gamma_w in N/m3
    }

    /// Ratio of permeability at e1 vs e0.
    pub fn permeability_ratio(&self, e0: f64, e1: f64) -> f64 {
        if e0 <= 0.0 {
            return 0.0;
        }
        let k0 = self.intrinsic_permeability(e0);
        let k1 = self.intrinsic_permeability(e1);
        if k0.abs() < 1e-30 { 0.0 } else { k1 / k0 }
    }

    /// Porosity n from void ratio e.
    pub fn porosity(e: f64) -> f64 {
        e / (1.0 + e)
    }

    /// Void ratio from porosity n.
    pub fn void_ratio_from_porosity(n: f64) -> f64 {
        if n >= 1.0 {
            return f64::INFINITY;
        }
        n / (1.0 - n)
    }
}

// ===========================================================================
// Terzaghi 1-D Consolidation
// ===========================================================================

/// Terzaghi one-dimensional consolidation model.
///
/// Solves the 1-D diffusion equation for excess pore-water pressure:
///   du/dt = c_v * d2u/dz2
#[derive(Debug, Clone)]
pub struct TerzaghiConsolidation {
    /// Coefficient of consolidation c_v \[m^2/s\].
    pub cv: f64,
    /// Drainage path length H_dr \[m\].
    pub h_dr: f64,
    /// Initial excess pore pressure u_0 \[Pa\].
    pub u0: f64,
}

impl TerzaghiConsolidation {
    /// Create a Terzaghi consolidation model.
    pub fn new(cv: f64, h_dr: f64, u0: f64) -> Self {
        Self { cv, h_dr, u0 }
    }

    /// Time factor T_v = c_v * t / H_dr^2.
    pub fn time_factor(&self, t: f64) -> f64 {
        self.cv * t / (self.h_dr * self.h_dr)
    }

    /// Degree of consolidation U(T_v) using series approximation.
    pub fn degree_of_consolidation(&self, t: f64) -> f64 {
        let tv = self.time_factor(t);
        if tv <= 0.0 {
            return 0.0;
        }
        // Series solution: U = 1 - sum_{m=0}^{N} (2/M^2) exp(-M^2 Tv)
        // where M = pi/2 * (2m+1)
        let mut u = 0.0;
        for m in 0..50 {
            let big_m = PI / 2.0 * (2.0 * m as f64 + 1.0);
            let term = (2.0 / (big_m * big_m)) * (-big_m * big_m * tv).exp();
            u += term;
            if term.abs() < 1e-15 {
                break;
            }
        }
        clamp_f64(1.0 - u, 0.0, 1.0)
    }

    /// Excess pore pressure at depth z and time t (double drainage).
    pub fn pore_pressure(&self, z: f64, t: f64) -> f64 {
        let tv = self.time_factor(t);
        if tv <= 0.0 {
            return self.u0;
        }
        let mut u = 0.0;
        for m in 0..50 {
            let big_m = PI / 2.0 * (2.0 * m as f64 + 1.0);
            let term = (2.0 * self.u0 / big_m)
                * (big_m * z / self.h_dr).sin()
                * (-big_m * big_m * tv).exp();
            u += term;
            if term.abs() < 1e-15 {
                break;
            }
        }
        clamp_f64(u, 0.0, self.u0)
    }

    /// Settlement at time t: s(t) = U(t) * s_final.
    pub fn settlement(&self, t: f64, s_final: f64) -> f64 {
        self.degree_of_consolidation(t) * s_final
    }

    /// Time for a given degree of consolidation U (0..1).
    /// Uses the approximate formula: t = T_v * H_dr^2 / c_v.
    pub fn time_for_consolidation(&self, u_target: f64) -> f64 {
        let u_c = clamp_f64(u_target, 0.0, 0.999);
        // Approximate Tv from U
        let tv = if u_c < 0.6 {
            PI / 4.0 * u_c * u_c
        } else {
            -0.9332 * (1.0 - u_c).ln() - 0.0851
        };
        tv * self.h_dr * self.h_dr / self.cv
    }

    /// Coefficient of volume compressibility m_v from compression index Cc,
    /// initial void ratio e0, and stress increment.
    pub fn mv_from_cc(cc: f64, e0: f64, sigma0: f64, delta_sigma: f64) -> f64 {
        if delta_sigma <= 0.0 || sigma0 <= 0.0 {
            return 0.0;
        }
        let de = cc * ((sigma0 + delta_sigma) / sigma0).log10();
        de / ((1.0 + e0) * delta_sigma)
    }
}

// ===========================================================================
// Swelling Clay Model
// ===========================================================================

/// Expansive (swelling) clay model with suction-dependent volumetric strain.
///
/// Models swell/shrink behaviour as a function of matric suction change.
#[derive(Debug, Clone)]
pub struct SwellingClayModel {
    /// Swelling index C_s (slope of e vs log(suction) on wetting).
    pub swelling_index: f64,
    /// Shrinkage limit suction \[Pa\].
    pub shrinkage_limit_suction: f64,
    /// Initial void ratio.
    pub e0: f64,
    /// Initial matric suction \[Pa\].
    pub suction0: f64,
    /// Maximum swelling strain (cap).
    pub max_swell_strain: f64,
}

impl SwellingClayModel {
    /// Create a new swelling clay model.
    pub fn new(
        swelling_index: f64,
        shrinkage_limit_suction: f64,
        e0: f64,
        suction0: f64,
        max_swell_strain: f64,
    ) -> Self {
        Self {
            swelling_index,
            shrinkage_limit_suction,
            e0,
            suction0,
            max_swell_strain,
        }
    }

    /// Volumetric swelling strain for suction change from suction0 to s_new.
    /// Negative suction change (wetting) -> positive swell strain.
    pub fn volumetric_strain(&self, s_new: f64) -> f64 {
        if s_new <= 0.0 || self.suction0 <= 0.0 {
            return 0.0;
        }
        let de = -self.swelling_index * (s_new / self.suction0).log10();
        let eps_v = de / (1.0 + self.e0);
        clamp_f64(eps_v, -self.max_swell_strain, self.max_swell_strain)
    }

    /// Void ratio at new suction level.
    pub fn void_ratio_at(&self, s_new: f64) -> f64 {
        if s_new <= 0.0 || self.suction0 <= 0.0 {
            return self.e0;
        }
        let de = -self.swelling_index * (s_new / self.suction0).log10();
        (self.e0 + de).max(0.0)
    }

    /// Swell pressure estimate \[Pa\]: suction at which soil starts to swell
    /// if loaded at given vertical stress.
    pub fn swell_pressure(&self) -> f64 {
        // Simplified: swell pressure ~ suction at which de = 0 starting from
        // shrinkage limit.
        self.shrinkage_limit_suction
    }

    /// Check if shrinkage limit is reached.
    pub fn is_at_shrinkage_limit(&self, suction: f64) -> bool {
        suction >= self.shrinkage_limit_suction
    }

    /// Linear swell potential from plasticity index (empirical).
    /// Returns percent swell.
    pub fn linear_swell_from_pi(plasticity_index: f64) -> f64 {
        // Seed et al. (1962): swell% ~ 2.16e-3 * PI^2.44
        2.16e-3 * plasticity_index.powf(2.44)
    }
}

// ===========================================================================
// Liquefaction Criteria
// ===========================================================================

/// Simplified liquefaction screening based on cyclic stress ratio.
///
/// Uses the Seed-Idriss simplified procedure for evaluating liquefaction
/// potential under earthquake loading.
#[derive(Debug, Clone)]
pub struct LiquefactionCriteria {
    /// Earthquake magnitude M_w.
    pub magnitude: f64,
    /// Peak ground acceleration a_max \[g\].
    pub pga: f64,
    /// Total vertical stress at depth \[Pa\].
    pub sigma_v: f64,
    /// Effective vertical stress at depth \[Pa\].
    pub sigma_v_eff: f64,
    /// SPT blow count N1_60 (corrected).
    pub n1_60: f64,
    /// Fines content \[%\].
    pub fines_content: f64,
}

impl LiquefactionCriteria {
    /// Create liquefaction criteria for a given site.
    pub fn new(
        magnitude: f64,
        pga: f64,
        sigma_v: f64,
        sigma_v_eff: f64,
        n1_60: f64,
        fines_content: f64,
    ) -> Self {
        Self {
            magnitude,
            pga,
            sigma_v,
            sigma_v_eff,
            n1_60,
            fines_content,
        }
    }

    /// Cyclic stress ratio CSR (Seed-Idriss).
    pub fn csr(&self) -> f64 {
        if self.sigma_v_eff <= 0.0 {
            return 0.0;
        }
        let rd = self.depth_reduction_factor();
        0.65 * (self.sigma_v / self.sigma_v_eff) * self.pga * rd
    }

    /// Depth reduction factor r_d (simplified linear approximation).
    fn depth_reduction_factor(&self) -> f64 {
        // Use sigma_v / gamma to estimate depth (assume gamma ~ 18 kN/m3)
        let _depth = self.sigma_v / 18_000.0;
        let z = self.sigma_v / 18_000.0;
        if z <= 9.15 {
            1.0 - 0.00765 * z
        } else if z <= 23.0 {
            1.174 - 0.0267 * z
        } else {
            0.5 // conservative
        }
    }

    /// Cyclic resistance ratio CRR from corrected SPT (clean sand).
    pub fn crr_clean_sand(&self) -> f64 {
        let n = self.n1_60_cs();
        if n >= 30.0 {
            return f64::INFINITY; // non-liquefiable
        }
        // Youd & Idriss (2001) curve
        let a = n / 34.8;
        let b = n / 135.0;
        let c = 50.0 / (10.0 * n + 45.0).powi(2);
        1.0 / (a + b + c) * 0.048 + 0.004 * n
    }

    /// Clean-sand equivalent N1_60_cs (fines correction).
    pub fn n1_60_cs(&self) -> f64 {
        let fc = self.fines_content;
        let alpha = if fc <= 5.0 {
            0.0
        } else if fc >= 35.0 {
            5.0
        } else {
            (fc - 5.0) / 6.0
        };
        let beta = if fc <= 5.0 {
            1.0
        } else if fc >= 35.0 {
            1.2
        } else {
            1.0 + 0.2 * (fc - 5.0) / 30.0
        };
        alpha + beta * self.n1_60
    }

    /// Magnitude scaling factor MSF.
    pub fn magnitude_scaling_factor(&self) -> f64 {
        10.0_f64.powf(2.24) / self.magnitude.powf(2.56)
    }

    /// Factor of safety against liquefaction.
    pub fn factor_of_safety(&self) -> f64 {
        let csr = self.csr();
        if csr <= 0.0 {
            return f64::INFINITY;
        }
        let crr = self.crr_clean_sand();
        let msf = self.magnitude_scaling_factor();
        crr * msf / csr
    }

    /// Is the site liquefiable? (FoS < 1.0).
    pub fn is_liquefiable(&self) -> bool {
        self.factor_of_safety() < 1.0
    }
}

// ===========================================================================
// Rockfill Model
// ===========================================================================

/// Nonlinear rockfill material model.
///
/// Models coarse granular fill with stress-dependent friction angle and
/// stiffness (e.g., for earth/rock dams).
#[derive(Debug, Clone)]
pub struct RockfillModel {
    /// Base friction angle at reference stress \[radians\].
    pub phi_base: f64,
    /// Friction angle reduction parameter delta_phi.
    pub delta_phi: f64,
    /// Reference confining stress for phi reduction \[Pa\].
    pub sigma_ref: f64,
    /// Modulus number K_mod (Janbu).
    pub modulus_number: f64,
    /// Modulus exponent n_mod (Janbu).
    pub modulus_exponent: f64,
    /// Atmospheric pressure for normalization \[Pa\].
    pub p_atm: f64,
    /// Maximum void ratio (loose state).
    pub e_max: f64,
    /// Minimum void ratio (dense state).
    pub e_min: f64,
}

impl RockfillModel {
    /// Create a new rockfill model.
    pub fn new(
        phi_base_deg: f64,
        delta_phi_deg: f64,
        sigma_ref: f64,
        modulus_number: f64,
        modulus_exponent: f64,
    ) -> Self {
        Self {
            phi_base: phi_base_deg.to_radians(),
            delta_phi: delta_phi_deg.to_radians(),
            sigma_ref,
            modulus_number,
            modulus_exponent,
            p_atm: 101_325.0,
            e_max: 0.9,
            e_min: 0.4,
        }
    }

    /// Stress-dependent friction angle \[radians\] at confining pressure sigma3.
    pub fn friction_angle(&self, sigma3: f64) -> f64 {
        if sigma3 <= 0.0 {
            return self.phi_base;
        }
        let reduction = self.delta_phi * (sigma3 / self.sigma_ref).log10().max(0.0);
        (self.phi_base - reduction).max(0.1_f64.to_radians())
    }

    /// Janbu modulus E \[Pa\] at given confining stress.
    pub fn youngs_modulus(&self, sigma3: f64) -> f64 {
        if sigma3 <= 0.0 {
            return self.modulus_number * self.p_atm;
        }
        self.modulus_number * self.p_atm * (sigma3 / self.p_atm).powf(self.modulus_exponent)
    }

    /// Shear strength at confining pressure.
    pub fn shear_strength(&self, sigma3: f64) -> f64 {
        sigma3 * self.friction_angle(sigma3).tan()
    }

    /// Relative density from void ratio.
    pub fn relative_density(&self, e: f64) -> f64 {
        if (self.e_max - self.e_min).abs() < 1e-15 {
            return 0.5;
        }
        clamp_f64((self.e_max - e) / (self.e_max - self.e_min), 0.0, 1.0)
    }

    /// Breakage index: fraction of particles crushed (simplified).
    pub fn breakage_index(&self, sigma3: f64) -> f64 {
        // Empirical: B = a * log(sigma3 / p_atm)
        if sigma3 <= self.p_atm {
            return 0.0;
        }
        let b = 0.15 * (sigma3 / self.p_atm).log10();
        clamp_f64(b, 0.0, 1.0)
    }
}

// ===========================================================================
// Frozen Soil Model
// ===========================================================================

/// Temperature-dependent frozen soil model.
///
/// Accounts for ice-bonding strength and creep at sub-zero temperatures.
#[derive(Debug, Clone)]
pub struct FrozenSoilModel {
    /// Unfrozen soil cohesion \[Pa\].
    pub c_unfrozen: f64,
    /// Unfrozen friction angle \[radians\].
    pub phi_unfrozen: f64,
    /// Ice bonding strength coefficient \[Pa / C\].
    pub ice_strength_coeff: f64,
    /// Creep exponent n.
    pub creep_n: f64,
    /// Creep reference rate A \[1/s\].
    pub creep_a: f64,
    /// Temperature at which freezing starts \[C\].
    pub freeze_temp: f64,
    /// Unfrozen water content model parameter alpha.
    pub unfrozen_alpha: f64,
    /// Unfrozen water content model parameter beta.
    pub unfrozen_beta: f64,
}

impl FrozenSoilModel {
    /// Create a new frozen soil model.
    pub fn new(
        c_unfrozen: f64,
        phi_unfrozen_deg: f64,
        ice_strength_coeff: f64,
        creep_n: f64,
        creep_a: f64,
    ) -> Self {
        Self {
            c_unfrozen,
            phi_unfrozen: phi_unfrozen_deg.to_radians(),
            ice_strength_coeff,
            creep_n,
            creep_a,
            freeze_temp: 0.0,
            unfrozen_alpha: 0.1,
            unfrozen_beta: -0.5,
        }
    }

    /// Total cohesion at temperature T \[C\].
    pub fn cohesion(&self, temp_c: f64) -> f64 {
        if temp_c >= self.freeze_temp {
            return self.c_unfrozen;
        }
        let dt = self.freeze_temp - temp_c;
        self.c_unfrozen + self.ice_strength_coeff * dt
    }

    /// Friction angle \[radians\] at temperature T \[C\].
    /// Ice bonding slightly reduces effective friction at very low temperatures.
    pub fn friction_angle(&self, temp_c: f64) -> f64 {
        if temp_c >= self.freeze_temp {
            return self.phi_unfrozen;
        }
        let dt = self.freeze_temp - temp_c;
        // Slight reduction at extreme cold
        (self.phi_unfrozen - 0.005 * dt).max(5.0_f64.to_radians())
    }

    /// Uniaxial compressive strength at temperature T \[C\].
    pub fn ucs(&self, temp_c: f64) -> f64 {
        let c = self.cohesion(temp_c);
        let phi = self.friction_angle(temp_c);
        2.0 * c * phi.cos() / (1.0 - phi.sin())
    }

    /// Creep strain rate \[1/s\] at deviatoric stress q and temperature T.
    pub fn creep_rate(&self, q: f64, temp_c: f64) -> f64 {
        if temp_c >= self.freeze_temp || q <= 0.0 {
            return 0.0;
        }
        // Arrhenius-type temperature dependence
        let dt = (self.freeze_temp - temp_c).max(0.01);
        let temp_factor = (-0.05 * dt).exp(); // simplified
        self.creep_a * (q / 1e6).powf(self.creep_n) * temp_factor
    }

    /// Unfrozen water content at temperature T \[C\].
    pub fn unfrozen_water_content(&self, temp_c: f64) -> f64 {
        if temp_c >= self.freeze_temp {
            return 1.0; // all water is liquid
        }
        let dt = (self.freeze_temp - temp_c).max(0.01);
        clamp_f64(self.unfrozen_alpha * dt.powf(self.unfrozen_beta), 0.0, 1.0)
    }

    /// Thermal conductivity \[W/(m K)\] using geometric mean of frozen/unfrozen.
    pub fn thermal_conductivity(&self, temp_c: f64, k_soil: f64, k_ice: f64, k_water: f64) -> f64 {
        let wu = self.unfrozen_water_content(temp_c);
        let wi = 1.0 - wu;
        // Geometric mean: k = k_soil^(1-n) * k_water^(n*wu) * k_ice^(n*wi)
        // Simplified: assume porosity n = 0.4
        let n = 0.4;
        k_soil.powf(1.0 - n) * k_water.powf(n * wu) * k_ice.powf(n * wi)
    }

    /// Is the soil frozen at temperature T?
    pub fn is_frozen(&self, temp_c: f64) -> bool {
        temp_c < self.freeze_temp
    }
}

// ===========================================================================
// Additional helpers
// ===========================================================================

/// Compute Drucker-Prager parameters from Mohr-Coulomb (inner cone).
///
/// Returns (alpha, k) for the DP criterion: q = alpha * p + k.
pub fn mohr_coulomb_to_drucker_prager(cohesion: f64, friction_deg: f64) -> (f64, f64) {
    let phi = friction_deg.to_radians();
    let sin_phi = phi.sin();
    let cos_phi = phi.cos();
    let alpha = 6.0 * sin_phi / (3.0 - sin_phi);
    let k = 6.0 * cohesion * cos_phi / (3.0 - sin_phi);
    (alpha, k)
}

/// Compute effective stress from total stress and pore pressure.
pub fn effective_stress(total_stress: f64, pore_pressure: f64) -> f64 {
    total_stress - pore_pressure
}

/// Compute void ratio from dry density and specific gravity.
pub fn void_ratio_from_density(dry_density: f64, specific_gravity: f64) -> f64 {
    let rho_w = 1000.0; // kg/m3
    specific_gravity * rho_w / dry_density - 1.0
}

/// Compute compression index Cc from liquid limit (Terzaghi & Peck).
pub fn compression_index_from_ll(liquid_limit: f64) -> f64 {
    0.009 * (liquid_limit - 10.0)
}

/// Compute coefficient of earth pressure at rest K0 (Jaky).
pub fn k0_jaky(friction_deg: f64) -> f64 {
    1.0 - friction_deg.to_radians().sin()
}

/// Compute bearing capacity (Terzaghi) for strip footing on cohesive soil.
///
/// q_ult = c * Nc + gamma * D * Nq + 0.5 * gamma * B * Ngamma
pub fn terzaghi_bearing_capacity(
    cohesion: f64,
    friction_deg: f64,
    gamma: f64,
    depth: f64,
    width: f64,
) -> f64 {
    let phi = friction_deg.to_radians();
    let tan_phi = phi.tan();
    // Bearing capacity factors (approximate)
    let nq = if phi.abs() < 1e-10 {
        1.0
    } else {
        (PI * tan_phi).exp() * (PI / 4.0 + phi / 2.0).tan().powi(2)
    };
    let nc = if phi.abs() < 1e-10 {
        5.14
    } else {
        (nq - 1.0) / tan_phi
    };
    let ngamma = 2.0 * (nq + 1.0) * tan_phi;
    cohesion * nc + gamma * depth * nq + 0.5 * gamma * width * ngamma
}

/// Rankine active earth pressure coefficient.
pub fn ka_rankine(friction_deg: f64) -> f64 {
    let phi = friction_deg.to_radians();
    let t = (PI / 4.0 - phi / 2.0).tan();
    t * t
}

/// Rankine passive earth pressure coefficient.
pub fn kp_rankine(friction_deg: f64) -> f64 {
    let phi = friction_deg.to_radians();
    let t = (PI / 4.0 + phi / 2.0).tan();
    t * t
}

/// Specific surface area from grain size (spherical particles).
pub fn specific_surface(d: f64) -> f64 {
    if d <= 0.0 {
        return 0.0;
    }
    6.0 / d
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ---- Modified Cam-Clay ----

    #[test]
    fn test_mcc_yield_inside() {
        let mcc = ModifiedCamClay::new(1.0, 0.1, 0.02, 200.0, 2.0, 0.3);
        // p = 100, q = 0 -> well inside
        let f = mcc.yield_function(100.0, 0.0);
        assert!(f < 0.0, "Should be inside yield surface, f = {f}");
    }

    #[test]
    fn test_mcc_yield_on_surface() {
        let mcc = ModifiedCamClay::new(1.0, 0.1, 0.02, 200.0, 2.0, 0.3);
        // On the ellipse apex: p = pc, q = 0 -> f = pc*(pc-pc) = 0
        let f = mcc.yield_function(200.0, 0.0);
        assert!(f.abs() < TOL, "Should be on yield surface, f = {f}");
    }

    #[test]
    fn test_mcc_yield_outside() {
        let mcc = ModifiedCamClay::new(1.0, 0.1, 0.02, 200.0, 2.0, 0.3);
        // Large q at small p
        let f = mcc.yield_function(50.0, 200.0);
        assert!(f > 0.0, "Should be outside yield surface, f = {f}");
    }

    #[test]
    fn test_mcc_hardening() {
        let mut mcc = ModifiedCamClay::new(1.0, 0.1, 0.02, 200.0, 2.0, 0.3);
        let old_pc = mcc.pc;
        mcc.harden(0.01);
        assert!(mcc.pc > old_pc, "pc should increase after compression");
    }

    #[test]
    fn test_mcc_ocr() {
        let mcc = ModifiedCamClay::new(1.0, 0.1, 0.02, 200.0, 2.0, 0.3);
        let ocr = mcc.ocr(100.0);
        assert!((ocr - 2.0).abs() < TOL, "OCR should be 2.0, got {ocr}");
    }

    #[test]
    fn test_mcc_bulk_modulus() {
        let mcc = ModifiedCamClay::new(1.0, 0.1, 0.02, 200.0, 2.0, 0.3);
        let k = mcc.bulk_modulus(100.0);
        // K = v0 * p / kappa = 2.0 * 100.0 / 0.02 = 10000
        assert!((k - 10_000.0).abs() < TOL, "K should be 10000, got {k}");
    }

    // ---- Mohr-Coulomb ----

    #[test]
    fn test_mc_shear_strength() {
        let mc = MohrCoulombModel::new(50_000.0, 30.0, 0.0, 0.0);
        let tau = mc.shear_strength(100_000.0);
        let expected = 50_000.0 + 100_000.0 * 30.0_f64.to_radians().tan();
        assert!(
            (tau - expected).abs() < 1.0,
            "tau = {tau}, expected = {expected}"
        );
    }

    #[test]
    fn test_mc_failure_check() {
        let mc = MohrCoulombModel::new(0.0, 30.0, 0.0, 0.0);
        // At failure: sigma1 = sigma3 * Nφ
        let n_phi = (1.0 + 30.0_f64.to_radians().sin()) / (1.0 - 30.0_f64.to_radians().sin());
        let sigma3 = 100_000.0;
        let sigma1 = sigma3 * n_phi;
        assert!(
            mc.is_failed_principals(sigma1, sigma3),
            "Should be at failure"
        );
    }

    #[test]
    fn test_mc_elastic_state() {
        let mc = MohrCoulombModel::new(100_000.0, 30.0, 0.0, 0.0);
        assert!(
            !mc.is_failed_principals(150_000.0, 100_000.0),
            "Should be elastic with high cohesion"
        );
    }

    // ---- Cap Plasticity ----

    #[test]
    fn test_cap_elastic() {
        // Use positive p convention. alpha=0.5, k=50000, R=0.5, X0=-200000
        // shear: q < alpha*p + k = 0.5*p + 50000
        // cap L = X - R*(alpha*X + k) = -200000 - 0.5*(0.5*(-200000)+50000)
        //       = -200000 - 0.5*(-50000) = -200000 + 25000 = -175000
        // cap: sqrt((p - L)^2 + (Rq)^2) < (X - L) = (-200000 - (-175000)) = -25000
        // That's negative, so let's use more reasonable params.
        let cap = CapPlasticityModel::new(0.3, 20_000.0, 0.5, 200_000.0, 0.1, 0.001);
        // L = 200000 - 0.5*(0.3*200000 + 20000) = 200000 - 0.5*80000 = 160000
        // X - L = 40000
        // At p=180000, q=5000: shear = 5000 - 0.3*180000 - 20000 = 5000 - 54000 - 20000 < 0 OK
        // cap: sqrt((180000-160000)^2 + (0.5*5000)^2) = sqrt(400e6 + 6.25e6) ~ 20156 < 40000 OK
        assert!(cap.is_elastic(180_000.0, 5_000.0), "Should be elastic");
    }

    #[test]
    fn test_cap_shear_yield() {
        let cap = CapPlasticityModel::new(0.5, 50_000.0, 0.5, -200_000.0, 0.1, 0.001);
        let f = cap.shear_yield(-10_000.0, 200_000.0);
        assert!(f > 0.0, "Should exceed shear yield, f = {f}");
    }

    #[test]
    fn test_cap_hardening() {
        let mut cap = CapPlasticityModel::new(0.5, 50_000.0, 0.5, -200_000.0, 0.1, 0.001);
        let old_x = cap.x_cap;
        cap.harden_cap(0.01);
        assert!(
            (cap.x_cap - old_x).abs() > 1e-10,
            "Cap should move after hardening"
        );
    }

    // ---- Hoek-Brown ----

    #[test]
    fn test_hb_ucs_intact() {
        // GSI = 100, D = 0 -> m_b = m_i, s = 1, a ~ 0.5
        let hb = HoekBrownRock::new(100e6, 25.0, 100.0, 0.0);
        let ucs = hb.ucs_rock_mass();
        // s^a ~ 1^0.5 = 1, so UCS_rm ~ sigma_ci
        assert!(
            (ucs - 100e6).abs() / 100e6 < 0.05,
            "UCS_rm should be ~ sigma_ci for GSI=100, got {ucs}"
        );
    }

    #[test]
    fn test_hb_sigma1_at_failure() {
        let hb = HoekBrownRock::new(100e6, 25.0, 75.0, 0.0);
        let s1 = hb.sigma1_at_failure(10e6);
        assert!(s1 > 10e6, "sigma1 at failure must exceed sigma3");
    }

    #[test]
    fn test_hb_tensile_strength() {
        let hb = HoekBrownRock::new(100e6, 25.0, 75.0, 0.0);
        let ts = hb.tensile_strength();
        assert!(ts < 0.0, "Tensile strength should be negative (tension)");
    }

    // ---- Kozeny-Carman ----

    #[test]
    fn test_kc_permeability_increases_with_e() {
        let kc = KozenyCarmanPermeability::default_water(0.001);
        let k1 = kc.intrinsic_permeability(0.5);
        let k2 = kc.intrinsic_permeability(1.0);
        assert!(
            k2 > k1,
            "Permeability should increase with void ratio: k1={k1}, k2={k2}"
        );
    }

    #[test]
    fn test_kc_zero_void_ratio() {
        let kc = KozenyCarmanPermeability::default_water(0.001);
        let k = kc.intrinsic_permeability(0.0);
        assert!(k.abs() < TOL, "Zero void ratio -> zero permeability");
    }

    #[test]
    fn test_kc_porosity_conversion() {
        let n = KozenyCarmanPermeability::porosity(1.0);
        assert!((n - 0.5).abs() < TOL, "e=1 -> n=0.5, got {n}");
        let e = KozenyCarmanPermeability::void_ratio_from_porosity(0.5);
        assert!((e - 1.0).abs() < TOL, "n=0.5 -> e=1.0, got {e}");
    }

    // ---- Terzaghi Consolidation ----

    #[test]
    fn test_terzaghi_initial() {
        let tc = TerzaghiConsolidation::new(1e-7, 5.0, 100_000.0);
        let u = tc.degree_of_consolidation(0.0);
        assert!(u.abs() < TOL, "U at t=0 should be 0, got {u}");
    }

    #[test]
    fn test_terzaghi_long_time() {
        let tc = TerzaghiConsolidation::new(1e-7, 1.0, 100_000.0);
        // Very long time
        let u = tc.degree_of_consolidation(1e12);
        assert!(
            (u - 1.0).abs() < 0.01,
            "U should approach 1.0 at long time, got {u}"
        );
    }

    #[test]
    fn test_terzaghi_settlement() {
        let tc = TerzaghiConsolidation::new(1e-7, 1.0, 100_000.0);
        let s = tc.settlement(1e12, 0.5);
        assert!(
            (s - 0.5).abs() < 0.01,
            "Final settlement should be 0.5 m, got {s}"
        );
    }

    // ---- Swelling Clay ----

    #[test]
    fn test_swell_wetting() {
        let sc = SwellingClayModel::new(0.1, 1e6, 0.8, 500_000.0, 0.2);
        // Wetting: reduce suction from 500 kPa to 100 kPa
        let eps = sc.volumetric_strain(100_000.0);
        assert!(eps > 0.0, "Wetting should cause positive swell, got {eps}");
    }

    #[test]
    fn test_swell_drying() {
        let sc = SwellingClayModel::new(0.1, 1e6, 0.8, 500_000.0, 0.2);
        // Drying: increase suction
        let eps = sc.volumetric_strain(1_000_000.0);
        assert!(eps < 0.0, "Drying should cause shrinkage, got {eps}");
    }

    // ---- Liquefaction ----

    #[test]
    fn test_liquefaction_csr() {
        let liq = LiquefactionCriteria::new(7.5, 0.3, 100_000.0, 60_000.0, 10.0, 5.0);
        let csr = liq.csr();
        assert!(csr > 0.0, "CSR should be positive, got {csr}");
        assert!(csr < 1.0, "CSR should be < 1 for typical conditions");
    }

    #[test]
    fn test_liquefaction_dense_sand() {
        // Dense sand with high N1_60 should not liquefy
        let liq = LiquefactionCriteria::new(7.5, 0.15, 100_000.0, 80_000.0, 35.0, 5.0);
        assert!(
            !liq.is_liquefiable(),
            "Dense sand (N1_60=35) should not liquefy"
        );
    }

    // ---- Rockfill ----

    #[test]
    fn test_rockfill_phi_reduction() {
        let rf = RockfillModel::new(45.0, 8.0, 100_000.0, 500.0, 0.5);
        let phi_low = rf.friction_angle(100_000.0);
        let phi_high = rf.friction_angle(1_000_000.0);
        assert!(
            phi_high < phi_low,
            "Friction should decrease with stress: low={phi_low:.4}, high={phi_high:.4}"
        );
    }

    #[test]
    fn test_rockfill_modulus() {
        let rf = RockfillModel::new(45.0, 8.0, 100_000.0, 500.0, 0.5);
        let e1 = rf.youngs_modulus(100_000.0);
        let e2 = rf.youngs_modulus(400_000.0);
        assert!(e2 > e1, "Modulus should increase with stress");
    }

    // ---- Frozen Soil ----

    #[test]
    fn test_frozen_cohesion_increase() {
        let fs = FrozenSoilModel::new(20_000.0, 25.0, 50_000.0, 3.0, 1e-12);
        let c_warm = fs.cohesion(5.0); // unfrozen
        let c_cold = fs.cohesion(-10.0); // frozen
        assert!(
            c_cold > c_warm,
            "Cohesion should increase when frozen: warm={c_warm}, cold={c_cold}"
        );
    }

    #[test]
    fn test_frozen_ucs() {
        let fs = FrozenSoilModel::new(20_000.0, 25.0, 50_000.0, 3.0, 1e-12);
        let ucs_warm = fs.ucs(5.0);
        let ucs_cold = fs.ucs(-20.0);
        assert!(
            ucs_cold > ucs_warm,
            "UCS should increase when frozen: warm={ucs_warm}, cold={ucs_cold}"
        );
    }

    #[test]
    fn test_frozen_creep_rate() {
        let fs = FrozenSoilModel::new(20_000.0, 25.0, 50_000.0, 3.0, 1e-12);
        let rate = fs.creep_rate(1e6, -5.0);
        assert!(rate > 0.0, "Creep rate should be positive at sub-zero temp");
    }

    #[test]
    fn test_frozen_no_creep_above_zero() {
        let fs = FrozenSoilModel::new(20_000.0, 25.0, 50_000.0, 3.0, 1e-12);
        let rate = fs.creep_rate(1e6, 5.0);
        assert!(rate.abs() < TOL, "No creep above freezing, got {rate}");
    }

    // ---- Helper functions ----

    #[test]
    fn test_dp_conversion() {
        let (alpha, k) = mohr_coulomb_to_drucker_prager(50_000.0, 30.0);
        assert!(alpha > 0.0, "alpha should be positive");
        assert!(k > 0.0, "k should be positive");
    }

    #[test]
    fn test_effective_stress() {
        let es = effective_stress(200_000.0, 50_000.0);
        assert!(
            (es - 150_000.0).abs() < TOL,
            "Effective stress should be 150 kPa, got {es}"
        );
    }

    #[test]
    fn test_k0_jaky_value() {
        let k0 = k0_jaky(30.0);
        let expected = 1.0 - 30.0_f64.to_radians().sin();
        assert!(
            (k0 - expected).abs() < TOL,
            "K0 = {k0}, expected {expected}"
        );
    }

    #[test]
    fn test_bearing_capacity_positive() {
        let q = terzaghi_bearing_capacity(50_000.0, 30.0, 18_000.0, 1.0, 2.0);
        assert!(q > 0.0, "Bearing capacity should be positive, got {q}");
    }

    #[test]
    fn test_rankine_ka_kp_product() {
        let ka = ka_rankine(30.0);
        let kp = kp_rankine(30.0);
        // Ka * Kp = 1.0
        assert!(
            (ka * kp - 1.0).abs() < 0.01,
            "Ka * Kp should be ~1.0, got {:.6}",
            ka * kp
        );
    }

    #[test]
    fn test_mean_deviatoric_stress() {
        let p = mean_stress(300.0, 200.0, 100.0);
        assert!((p - 200.0).abs() < TOL, "Mean stress should be 200");
        let q = deviatoric_q(300.0, 200.0, 100.0);
        // q = sqrt(0.5*((300-200)^2 + (200-100)^2 + (100-300)^2))
        // = sqrt(0.5*(10000 + 10000 + 40000)) = sqrt(30000) ~ 173.2
        let expected = (30000.0_f64).sqrt();
        assert!((q - expected).abs() < 0.1, "q = {q}, expected ~ {expected}");
    }
}
