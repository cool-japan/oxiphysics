// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exchange-correlation kernels for closed-shell (spin-unpolarised) Kohn–Sham DFT.
//!
//! Every functional is written once as an energy density per unit volume
//! `f(ρ, σ, τ)` with `σ = |∇ρ|²` and `τ = ½ Σ_i |∇ψ_i|²`.  The derivatives
//! `∂f/∂ρ`, `∂f/∂σ` and `∂f/∂τ` that enter the Kohn–Sham potential are obtained
//! exactly with forward-mode automatic differentiation ([`XcDual`]), so the
//! potential is always consistent with the energy.
//!
//! Implemented functionals:
//!
//! * **LDA** – Slater exchange + VWN5 correlation (Vosko, Wilk & Nusair 1980).
//! * **PBE** – PBE exchange and correlation on top of PW92 (Perdew, Burke &
//!   Ernzerhof 1996; Perdew & Wang 1992).
//! * **B3LYP** – `0.08 Slater + 0.72 B88 + 0.20 HF` exchange and
//!   `0.19 VWN(RPA) + 0.81 LYP` correlation (Stephens et al. 1994; Becke 1988;
//!   Lee, Yang & Parr 1988 in the Miehlich form).
//! * **TPSS** – meta-GGA exchange and correlation (Tao, Perdew, Staroverov &
//!   Scuseria 2003).
//!
//! The fraction of exact (Hartree–Fock) exchange `a₀` is applied to the
//! semilocal exchange as `(1 − a₀)`; for B3LYP the B88 gradient correction keeps
//! its fixed weight 0.72.

use std::f64::consts::PI;
use std::ops::{Add, Div, Mul, Neg, Sub};

use super::types_3::XcFunctional;

/// Densities below this threshold contribute nothing to the XC energy.
pub const XC_DENSITY_THRESHOLD: f64 = 1e-14;
/// Smallest |∇ρ|² used inside the gradient-corrected kernels (keeps the
/// derivatives finite where the density gradient vanishes exactly).
const SIGMA_FLOOR: f64 = 1e-30;
/// B88 weight in B3LYP.
const B3LYP_B88_WEIGHT: f64 = 0.72;
/// VWN(RPA) correlation weight in B3LYP.
const B3LYP_VWN_WEIGHT: f64 = 0.19;
/// LYP correlation weight in B3LYP.
const B3LYP_LYP_WEIGHT: f64 = 0.81;
/// Default exact-exchange fraction of B3LYP.
pub const B3LYP_EXACT_EXCHANGE: f64 = 0.20;

/// Forward-mode dual number carrying derivatives with respect to `(ρ, σ, τ)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct XcDual {
    /// Value.
    pub(crate) v: f64,
    /// Partial derivatives ∂/∂ρ, ∂/∂σ, ∂/∂τ.
    pub(crate) d: [f64; 3],
}

impl XcDual {
    /// A constant (zero derivatives).
    pub(crate) fn constant(v: f64) -> Self {
        Self { v, d: [0.0; 3] }
    }
    /// The independent variable number `slot` (0 = ρ, 1 = σ, 2 = τ) with value `v`.
    pub(crate) fn variable(v: f64, slot: usize) -> Self {
        let mut d = [0.0; 3];
        if let Some(x) = d.get_mut(slot) {
            *x = 1.0;
        }
        Self { v, d }
    }
    fn chain(self, value: f64, derivative: f64) -> Self {
        Self {
            v: value,
            d: [
                derivative * self.d[0],
                derivative * self.d[1],
                derivative * self.d[2],
            ],
        }
    }
    pub(crate) fn powf(self, n: f64) -> Self {
        let value = self.v.powf(n);
        self.chain(value, n * self.v.powf(n - 1.0))
    }
    pub(crate) fn sqrt(self) -> Self {
        let value = self.v.sqrt();
        self.chain(value, 0.5 / value)
    }
    pub(crate) fn exp(self) -> Self {
        let value = self.v.exp();
        self.chain(value, value)
    }
    pub(crate) fn ln(self) -> Self {
        self.chain(self.v.ln(), 1.0 / self.v)
    }
    pub(crate) fn atan(self) -> Self {
        self.chain(self.v.atan(), 1.0 / (1.0 + self.v * self.v))
    }
    pub(crate) fn asinh(self) -> Self {
        self.chain(self.v.asinh(), 1.0 / (1.0 + self.v * self.v).sqrt())
    }
    /// The larger of two values (derivative of the selected branch).
    pub(crate) fn max(self, other: Self) -> Self {
        if self.v >= other.v { self } else { other }
    }
}

impl Add for XcDual {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d: [self.d[0] + o.d[0], self.d[1] + o.d[1], self.d[2] + o.d[2]],
        }
    }
}
impl Sub for XcDual {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            v: self.v - o.v,
            d: [self.d[0] - o.d[0], self.d[1] - o.d[1], self.d[2] - o.d[2]],
        }
    }
}
impl Mul for XcDual {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self {
            v: self.v * o.v,
            d: [
                self.d[0] * o.v + self.v * o.d[0],
                self.d[1] * o.v + self.v * o.d[1],
                self.d[2] * o.v + self.v * o.d[2],
            ],
        }
    }
}
impl Div for XcDual {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        let inv = 1.0 / o.v;
        let q = self.v * inv;
        Self {
            v: q,
            d: [
                (self.d[0] - q * o.d[0]) * inv,
                (self.d[1] - q * o.d[1]) * inv,
                (self.d[2] - q * o.d[2]) * inv,
            ],
        }
    }
}
impl Neg for XcDual {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: [-self.d[0], -self.d[1], -self.d[2]],
        }
    }
}
impl Add<f64> for XcDual {
    type Output = Self;
    fn add(self, o: f64) -> Self {
        Self {
            v: self.v + o,
            d: self.d,
        }
    }
}
impl Sub<f64> for XcDual {
    type Output = Self;
    fn sub(self, o: f64) -> Self {
        Self {
            v: self.v - o,
            d: self.d,
        }
    }
}
impl Mul<f64> for XcDual {
    type Output = Self;
    fn mul(self, o: f64) -> Self {
        Self {
            v: self.v * o,
            d: [self.d[0] * o, self.d[1] * o, self.d[2] * o],
        }
    }
}
impl Div<f64> for XcDual {
    type Output = Self;
    fn div(self, o: f64) -> Self {
        self * (1.0 / o)
    }
}
impl Add<XcDual> for f64 {
    type Output = XcDual;
    fn add(self, o: XcDual) -> XcDual {
        o + self
    }
}
impl Sub<XcDual> for f64 {
    type Output = XcDual;
    fn sub(self, o: XcDual) -> XcDual {
        -o + self
    }
}
impl Mul<XcDual> for f64 {
    type Output = XcDual;
    fn mul(self, o: XcDual) -> XcDual {
        o * self
    }
}
impl Div<XcDual> for f64 {
    type Output = XcDual;
    fn div(self, o: XcDual) -> XcDual {
        XcDual::constant(self) / o
    }
}

/// XC energy density and its partial derivatives at one point.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct XcPointValue {
    /// Energy per unit volume `f = ρ ε_xc` \[Hartree/Bohr³\].
    pub energy_density: f64,
    /// `∂f/∂ρ`.
    pub v_rho: f64,
    /// `∂f/∂σ` with `σ = |∇ρ|²`.
    pub v_sigma: f64,
    /// `∂f/∂τ` with `τ = ½ Σ_i |∇ψ_i|²`.
    pub v_tau: f64,
}

/// Whether a functional depends on the density gradient.
pub fn xc_needs_gradient(functional: XcFunctional) -> bool {
    !matches!(functional, XcFunctional::Lda)
}

/// Whether a functional depends on the kinetic-energy density τ.
pub fn xc_needs_tau(functional: XcFunctional) -> bool {
    matches!(functional, XcFunctional::Tpss)
}

/// Kinetic-energy density of the uniform electron gas, `(3/10)(3π²)^{2/3} ρ^{5/3}`.
pub fn uniform_gas_tau(rho: f64) -> f64 {
    0.3 * (3.0 * PI * PI).powf(2.0 / 3.0) * rho.powf(5.0 / 3.0)
}

/// Evaluate the semilocal part of `functional` at one point.
///
/// `exact_exchange` is the fraction `a₀` of Hartree–Fock exchange added by the
/// caller; the semilocal exchange is scaled by `1 − a₀`.  `sigma` and `tau` are
/// ignored by functionals that do not use them.
pub fn evaluate_xc(
    functional: XcFunctional,
    exact_exchange: f64,
    rho: f64,
    sigma: f64,
    tau: f64,
) -> XcPointValue {
    if rho.is_nan() || rho <= XC_DENSITY_THRESHOLD {
        return XcPointValue::default();
    }
    let tau = tau.max(0.0);
    // Meta-GGAs require τ ≥ τ_W = σ/(8ρ); like libxc, bound σ by 8ρτ so that
    // numerical noise in τ cannot produce an unphysical z = τ_W/τ > 1.
    let sigma = if xc_needs_tau(functional) {
        sigma.min(8.0 * rho * tau)
    } else {
        sigma
    };
    let r = XcDual::variable(rho, 0);
    let s = XcDual::variable(sigma.max(SIGMA_FLOOR), 1);
    let t = XcDual::variable(tau, 2);
    let x_scale = 1.0 - exact_exchange;
    let f = match functional {
        XcFunctional::Lda => slater_exchange(r) * x_scale + vwn_correlation(r, &VWN5),
        XcFunctional::Pbe => pbe_exchange(r, s) * x_scale + pbe_correlation(r, s),
        XcFunctional::B3lyp => {
            let slater = slater_exchange(r);
            slater * x_scale
                + (b88_exchange(r, s) - slater) * B3LYP_B88_WEIGHT
                + vwn_correlation(r, &VWN_RPA) * B3LYP_VWN_WEIGHT
                + lyp_correlation(r, s) * B3LYP_LYP_WEIGHT
        }
        XcFunctional::Tpss => tpss_exchange(r, s, t) * x_scale + tpss_correlation(r, s, t),
    };
    XcPointValue {
        energy_density: f.v,
        v_rho: f.d[0],
        v_sigma: f.d[1],
        v_tau: f.d[2],
    }
}

// ---------------------------------------------------------------------------
// LDA pieces
// ---------------------------------------------------------------------------

/// Wigner–Seitz radius `r_s = (3/(4πρ))^{1/3}`.
fn wigner_seitz_radius(rho: XcDual) -> XcDual {
    (XcDual::constant(3.0 / (4.0 * PI)) / rho).powf(1.0 / 3.0)
}

/// Slater (Dirac) exchange energy per volume `−(3/4)(3/π)^{1/3} ρ^{4/3}`.
fn slater_exchange(rho: XcDual) -> XcDual {
    rho.powf(4.0 / 3.0) * (-0.75 * (3.0 / PI).powf(1.0 / 3.0))
}

/// Parameters of a paramagnetic VWN interpolation.
pub(crate) struct VwnParams {
    a: f64,
    b: f64,
    c: f64,
    x0: f64,
}

/// VWN functional V (fit to Ceperley–Alder), paramagnetic.
pub(crate) const VWN5: VwnParams = VwnParams {
    a: 0.031_090_7,
    b: 3.727_44,
    c: 12.935_2,
    x0: -0.104_98,
};
/// VWN fit to the RPA correlation energy (used by B3LYP), paramagnetic.
pub(crate) const VWN_RPA: VwnParams = VwnParams {
    a: 0.031_090_7,
    b: 13.072_0,
    c: 42.719_8,
    x0: -0.409_286,
};

/// VWN correlation energy per particle as a function of `r_s`.
pub(crate) fn vwn_epsilon(rs: XcDual, p: &VwnParams) -> XcDual {
    let x = rs.sqrt();
    let big_x = x * x + x * p.b + p.c;
    let big_x0 = p.x0 * p.x0 + p.b * p.x0 + p.c;
    let q = (4.0 * p.c - p.b * p.b).sqrt();
    let at = (XcDual::constant(q) / (x * 2.0 + p.b)).atan();
    let term1 = (x * x / big_x).ln() + at * (2.0 * p.b / q);
    let xm = x - p.x0;
    let term2 = (xm * xm / big_x).ln() + at * (2.0 * (p.b + 2.0 * p.x0) / q);
    (term1 - term2 * (p.b * p.x0 / big_x0)) * p.a
}

/// VWN correlation energy per volume.
fn vwn_correlation(rho: XcDual, p: &VwnParams) -> XcDual {
    rho * vwn_epsilon(wigner_seitz_radius(rho), p)
}

/// Parameters of the PW92 G(r_s) interpolation.
struct Pw92Params {
    a: f64,
    alpha1: f64,
    beta: [f64; 4],
}

/// PW92 paramagnetic correlation (ζ = 0), "PW_MOD" constants as used by PBE.
const PW92_PARAMAGNETIC: Pw92Params = Pw92Params {
    a: 0.031_090_7,
    alpha1: 0.213_70,
    beta: [7.595_7, 3.587_6, 1.638_2, 0.492_94],
};
/// PW92 ferromagnetic correlation (ζ = 1).
const PW92_FERROMAGNETIC: Pw92Params = Pw92Params {
    a: 0.015_545_35,
    alpha1: 0.205_48,
    beta: [14.118_9, 6.197_7, 3.366_2, 0.625_17],
};

/// PW92 correlation energy per particle
/// `G = −2A(1 + α₁r_s) ln[1 + 1/(2A(β₁r_s^{1/2} + β₂r_s + β₃r_s^{3/2} + β₄r_s²))]`.
fn pw92_epsilon(rs: XcDual, p: &Pw92Params) -> XcDual {
    let srs = rs.sqrt();
    let den = (srs * p.beta[0] + rs * p.beta[1] + rs * srs * p.beta[2] + rs * rs * p.beta[3])
        * (2.0 * p.a);
    ((1.0 / den) + 1.0).ln() * ((rs * p.alpha1 + 1.0) * (-2.0 * p.a))
}

// ---------------------------------------------------------------------------
// GGA pieces
// ---------------------------------------------------------------------------

/// PBE parameters.
const PBE_KAPPA: f64 = 0.804;
const PBE_MU: f64 = 0.219_514_972_764_517_1;
const PBE_BETA: f64 = 0.066_724_550_603_149_22;

/// `(3π²)^{1/3}`.
fn three_pi_squared_cbrt() -> f64 {
    (3.0 * PI * PI).powf(1.0 / 3.0)
}

/// Squared reduced gradient `s² = σ / (4 (3π²)^{2/3} ρ^{8/3})`.
fn reduced_gradient_squared(rho: XcDual, sigma: XcDual) -> XcDual {
    let k = three_pi_squared_cbrt();
    sigma / (rho.powf(8.0 / 3.0) * (4.0 * k * k))
}

/// PBE exchange energy per volume.
fn pbe_exchange(rho: XcDual, sigma: XcDual) -> XcDual {
    let s2 = reduced_gradient_squared(rho, sigma);
    let fx = (1.0 + PBE_KAPPA) - PBE_KAPPA / (s2 * (PBE_MU / PBE_KAPPA) + 1.0);
    slater_exchange(rho) * fx
}

/// PBE gradient correction `H(r_s, ζ, t)` per particle for spin-scaling factor
/// `phi`, uniform-gas correlation `eps_unif`, and `t²`.
fn pbe_h(eps_unif: XcDual, phi: f64, t2: XcDual) -> XcDual {
    let gamma = (1.0 - std::f64::consts::LN_2) / (PI * PI);
    let phi3 = phi * phi * phi;
    let a = (PBE_BETA / gamma) / ((-eps_unif / (gamma * phi3)).exp() - 1.0);
    let at2 = a * t2;
    let frac = (at2 + 1.0) / (at2 + 1.0 + at2 * at2);
    (t2 * frac * (PBE_BETA / gamma) + 1.0).ln() * (gamma * phi3)
}

/// PBE correlation energy per particle of a spin-unpolarised density.
fn pbe_correlation_epsilon(rho: XcDual, sigma: XcDual) -> XcDual {
    let rs = wigner_seitz_radius(rho);
    let eps = pw92_epsilon(rs, &PW92_PARAMAGNETIC);
    let kf = (rho * (3.0 * PI * PI)).powf(1.0 / 3.0);
    let ks2 = kf * (4.0 / PI);
    let t2 = sigma / (ks2 * rho * rho * 4.0);
    eps + pbe_h(eps, 1.0, t2)
}

/// PBE correlation energy per particle of a fully spin-polarised density `n`
/// with gradient² `sigma` (ζ = 1, φ = 2^{−1/3}).
fn pbe_correlation_epsilon_polarised(n: XcDual, sigma: XcDual) -> XcDual {
    let phi = 2.0_f64.powf(-1.0 / 3.0);
    let rs = wigner_seitz_radius(n);
    let eps = pw92_epsilon(rs, &PW92_FERROMAGNETIC);
    let kf = (n * (3.0 * PI * PI)).powf(1.0 / 3.0);
    let ks2 = kf * (4.0 / PI);
    let t2 = sigma / (ks2 * n * n * (4.0 * phi * phi));
    eps + pbe_h(eps, phi, t2)
}

/// PBE correlation energy per volume.
fn pbe_correlation(rho: XcDual, sigma: XcDual) -> XcDual {
    rho * pbe_correlation_epsilon(rho, sigma)
}

/// `x asinh x` as a function of `y = x²`, with its Taylor series near zero.
fn x_asinh_x(y: XcDual) -> XcDual {
    if y.v < 1e-10 {
        // x asinh x = y − y²/6 + 3y³/40 − …
        y - y * y / 6.0 + y * y * y * (3.0 / 40.0)
    } else {
        let x = y.sqrt();
        x * x.asinh()
    }
}

/// Becke 1988 exchange energy per volume (closed shell: ρ_σ = ρ/2, σ_σσ = σ/4).
fn b88_exchange(rho: XcDual, sigma: XcDual) -> XcDual {
    const BETA: f64 = 0.0042;
    let rho_s = rho * 0.5;
    let sigma_s = sigma * 0.25;
    let rho_s43 = rho_s.powf(4.0 / 3.0);
    let x2 = sigma_s / (rho_s43 * rho_s43);
    let c_x = 1.5 * (3.0 / (4.0 * PI)).powf(1.0 / 3.0);
    let per_spin = -(rho_s43 * c_x) - rho_s43 * BETA * x2 / (x_asinh_x(x2) * (6.0 * BETA) + 1.0);
    per_spin * 2.0
}

/// Lee–Yang–Parr correlation energy per volume in the Laplacian-free form of
/// Miehlich, Savin, Stoll & Preuss (1989), evaluated for ρ_α = ρ_β = ρ/2.
fn lyp_correlation(rho: XcDual, sigma: XcDual) -> XcDual {
    const A: f64 = 0.049_18;
    const B: f64 = 0.132;
    const C: f64 = 0.253_3;
    const D: f64 = 0.349;
    let c_f = 0.3 * (3.0 * PI * PI).powf(2.0 / 3.0);
    let ra = rho * 0.5;
    let rb = rho * 0.5;
    // |∇ρ_α|² = |∇ρ_β|² = σ/4, |∇ρ|² = σ.
    let gaa = sigma * 0.25;
    let gbb = sigma * 0.25;
    let rm13 = rho.powf(-1.0 / 3.0);
    let denom = rm13 * D + 1.0;
    let omega = (rm13 * (-C)).exp() / denom * rho.powf(-11.0 / 3.0);
    let delta = rm13 * C + rm13 * D / denom;
    let ra83 = ra.powf(8.0 / 3.0);
    let rb83 = rb.powf(8.0 / 3.0);
    let bracket = (ra83 + rb83) * (2.0_f64.powf(11.0 / 3.0) * c_f)
        + sigma * (47.0 / 18.0 - delta * (7.0 / 18.0))
        - (gaa + gbb) * (2.5 - delta / 18.0)
        - (delta - 11.0) / 9.0 * (ra / rho * gaa + rb / rho * gbb);
    let rho2 = rho * rho;
    let rest = ra * rb * bracket - rho2 * sigma * (2.0 / 3.0)
        + (rho2 * (2.0 / 3.0) - ra * ra) * gbb
        + (rho2 * (2.0 / 3.0) - rb * rb) * gaa;
    -(ra * rb / rho * 4.0 * A / denom) - omega * rest * (A * B)
}

// ---------------------------------------------------------------------------
// meta-GGA pieces (TPSS)
// ---------------------------------------------------------------------------

/// Ratio `z = τ_W/τ` with the von Weizsäcker density `τ_W = σ/(8ρ)`
/// (σ has already been bounded by 8ρτ, so z ≤ 1 up to rounding).
fn weizsacker_ratio(rho: XcDual, sigma: XcDual, tau: XcDual) -> XcDual {
    if tau.v <= 0.0 {
        return XcDual::constant(1.0);
    }
    let z = sigma / (rho * 8.0) / tau;
    if z.v > 1.0 + 1e-10 {
        XcDual::constant(1.0)
    } else {
        z
    }
}

/// TPSS exchange energy per volume.
fn tpss_exchange(rho: XcDual, sigma: XcDual, tau: XcDual) -> XcDual {
    const B: f64 = 0.40;
    const C: f64 = 1.590_96;
    const E: f64 = 1.537;
    const KAPPA: f64 = 0.804;
    const MU: f64 = 0.219_51;
    let p = reduced_gradient_squared(rho, sigma);
    let z = weizsacker_ratio(rho, sigma, tau);
    // α = (τ − τ_W)/τ_unif = (5p/3)(1/z − 1)
    let alpha = p * (5.0 / 3.0) * (1.0 / z - 1.0);
    let qb = (alpha - 1.0) * 0.45 / (alpha * (alpha - 1.0) * B + 1.0).sqrt() + p * (2.0 / 3.0);
    let ten_81 = 10.0 / 81.0;
    let z2 = z * z;
    let opz2 = z2 + 1.0;
    let three_fifths_z = z * 0.6;
    let root = (three_fifths_z * three_fifths_z * 0.5 + p * p * 0.5).sqrt();
    let sqrt_e = E.sqrt();
    let numerator = (z2 * C / (opz2 * opz2) + ten_81) * p + qb * qb * (146.0 / 2025.0)
        - qb * root * (73.0 / 405.0)
        + p * p * (ten_81 * ten_81 / KAPPA)
        + three_fifths_z * three_fifths_z * (2.0 * sqrt_e * ten_81)
        + p * p * p * (E * MU);
    let one_plus = p * sqrt_e + 1.0;
    let x = numerator / (one_plus * one_plus);
    let fx = (1.0 + KAPPA) - KAPPA / (x / KAPPA + 1.0);
    slater_exchange(rho) * fx
}

/// TPSS correlation energy per volume (closed shell: ζ = 0, ξ = 0, C = 0.53).
fn tpss_correlation(rho: XcDual, sigma: XcDual, tau: XcDual) -> XcDual {
    const C0: f64 = 0.53;
    const D: f64 = 2.8;
    let z = weizsacker_ratio(rho, sigma, tau);
    let z2 = z * z;
    let eps_pbe = pbe_correlation_epsilon(rho, sigma);
    // ε̃_c^σ = max[ε_c^PBE(ρ_σ, 0, ∇ρ_σ, 0), ε_c^PBE(ρ_α, ρ_β, ∇ρ_α, ∇ρ_β)];
    // both spins are equal, so Σ_σ (ρ_σ/ρ) ε̃_c^σ = ε̃_c.
    let eps_pol = pbe_correlation_epsilon_polarised(rho * 0.5, sigma * 0.25);
    let eps_tilde = eps_pol.max(eps_pbe);
    let rev = eps_pbe * (z2 * C0 + 1.0) - z2 * eps_tilde * (1.0 + C0);
    let eps_c = rev * (rev * z2 * z * D + 1.0);
    rho * eps_c
}

#[cfg(test)]
mod tests {
    use super::*;

    /// (ρ, σ, τ) sample points.
    const POINTS: [(f64, f64, f64); 3] = [(0.1, 0.01, 0.05), (1.0, 0.5, 2.0), (0.01, 0.001, 0.003)];

    fn assert_close(label: &str, got: f64, want: f64) {
        let tol = 1e-9 * want.abs().max(1.0);
        assert!(
            (got - want).abs() < tol,
            "{label}: got {got:.15e}, want {want:.15e}"
        );
    }

    /// Reference values: PySCF 2.14 / libxc `eval_xc` (exc·ρ, vrho, vsigma, vtau).
    fn check(functional: XcFunctional, a0: f64, refs: [[f64; 4]; 3]) {
        for ((rho, sigma, tau), r) in POINTS.iter().zip(refs.iter()) {
            let v = evaluate_xc(functional, a0, *rho, *sigma, *tau);
            assert_close("energy density", v.energy_density, r[0]);
            assert_close("vrho", v.v_rho, r[1]);
            assert_close("vsigma", v.v_sigma, r[2]);
            assert_close("vtau", v.v_tau, r[3]);
        }
    }

    #[test]
    fn lda_matches_libxc() {
        check(
            XcFunctional::Lda,
            0.0,
            [
                [
                    -0.039_620_590_148_651_23,
                    -0.517_890_180_065_344_8,
                    0.0,
                    0.0,
                ],
                [-0.810_151_378_688_813, -1.064_683_405_018_682, 0.0, 0.0],
                [
                    -0.001_967_628_529_542_297,
                    -0.256_029_540_036_804_67,
                    0.0,
                    0.0,
                ],
            ],
        );
    }

    #[test]
    fn pbe_matches_libxc() {
        check(
            XcFunctional::Pbe,
            0.0,
            [
                [
                    -0.039_691_828_163_848_704,
                    -0.514_908_531_647_090_6,
                    -0.015_691_775_511_290_112,
                    0.0,
                ],
                [
                    -0.809_820_406_744_841_7,
                    -1.063_985_166_119_453_5,
                    -0.000_240_665_007_370_997_2,
                    0.0,
                ],
                [
                    -0.002_386_953_539_304_280_3,
                    -0.245_116_537_479_683_4,
                    -0.271_726_197_404_336_6,
                    0.0,
                ],
            ],
        );
    }

    #[test]
    fn b3lyp_semilocal_part_matches_libxc() {
        check(
            XcFunctional::B3lyp,
            B3LYP_EXACT_EXCHANGE,
            [
                [
                    -0.032_191_134_292_220_98,
                    -0.407_004_482_735_042_36,
                    -0.056_429_535_473_983_51,
                    0.0,
                ],
                [
                    -0.648_239_476_693_461_8,
                    -0.846_996_811_452_456_4,
                    -0.003_405_365_221_255_088_5,
                    0.0,
                ],
                [
                    -0.001_916_794_853_971_783_6,
                    -0.222_000_219_633_237_04,
                    -0.106_902_097_335_018_14,
                    0.0,
                ],
            ],
        );
    }

    #[test]
    fn tpss_matches_libxc() {
        check(
            XcFunctional::Tpss,
            0.0,
            [
                [
                    -0.039_693_348_557_112_665,
                    -0.522_824_200_633_744_4,
                    -0.056_474_863_110_709_76,
                    0.022_396_387_447_339_75,
                ],
                [
                    -0.810_311_092_781_997_4,
                    -1.076_105_677_568_709_2,
                    0.000_833_169_916_262_767_8,
                    0.003_009_415_060_717_61,
                ],
                [
                    -0.002_020_857_634_779_253,
                    -0.249_152_891_476_152_6,
                    -0.436_301_980_852_774_5,
                    0.023_236_070_592_688_535,
                ],
            ],
        );
    }

    #[test]
    fn derivatives_are_consistent_with_finite_differences() {
        for f in [
            XcFunctional::Lda,
            XcFunctional::Pbe,
            XcFunctional::B3lyp,
            XcFunctional::Tpss,
        ] {
            let (rho, sigma, tau) = (0.3, 0.07, 0.4);
            let v = evaluate_xc(f, 0.0, rho, sigma, tau);
            let h = 1e-6;
            let num = |dr: f64, ds: f64, dt: f64| {
                (evaluate_xc(f, 0.0, rho + dr, sigma + ds, tau + dt).energy_density
                    - evaluate_xc(f, 0.0, rho - dr, sigma - ds, tau - dt).energy_density)
                    / (2.0 * h)
            };
            assert!((num(h, 0.0, 0.0) - v.v_rho).abs() < 1e-6, "{f:?} vrho");
            assert!((num(0.0, h, 0.0) - v.v_sigma).abs() < 1e-6, "{f:?} vsigma");
            assert!((num(0.0, 0.0, h) - v.v_tau).abs() < 1e-6, "{f:?} vtau");
        }
    }

    #[test]
    fn vanishing_density_gives_zero() {
        let v = evaluate_xc(XcFunctional::Pbe, 0.0, 0.0, 0.0, 0.0);
        assert_eq!(v, XcPointValue::default());
    }

    #[test]
    fn exact_exchange_fraction_scales_semilocal_exchange() {
        let full = evaluate_xc(XcFunctional::Lda, 0.0, 0.5, 0.0, 0.0);
        let hybrid = evaluate_xc(XcFunctional::Lda, 0.25, 0.5, 0.0, 0.0);
        let ex = slater_exchange(XcDual::constant(0.5)).v;
        assert!((full.energy_density - hybrid.energy_density - 0.25 * ex).abs() < 1e-14);
    }
}
