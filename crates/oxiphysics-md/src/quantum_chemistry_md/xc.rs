// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Spin-unpolarised exchange-correlation functionals for the grid KS driver.
//!
//! Every functional is expressed as an energy density per volume
//! f(ρ, σ) with σ = |∇ρ|², so that E_xc = ∫ f dV:
//!
//! * [`XcFunctional::Lda`]: Slater exchange + VWN5 correlation
//!   (Vosko, Wilk & Nusair, Can. J. Phys. 58, 1200 (1980), paramagnetic fit).
//! * [`XcFunctional::Pbe`]: PBE exchange and correlation (Perdew, Burke &
//!   Ernzerhof, PRL 77, 3865 (1996)) on top of PW92 correlation (Perdew &
//!   Wang, PRB 45, 13244 (1992)).
//! * [`XcFunctional::B3lyp`]: semi-local part of B3LYP,
//!   0.8 E_x^Slater + 0.72 ΔE_x^B88 + 0.19 E_c^VWN(RPA) + 0.81 E_c^LYP,
//!   with the VWN RPA parametrisation used by the original Gaussian
//!   implementation (Stephens et al., J. Phys. Chem. 98, 11623 (1994)).  The
//!   remaining 0.2 E_x^HF is orbital dependent and handled by the driver.
//!
//! LDA potentials are analytic.  GGA potentials need ∂f/∂ρ and ∂f/∂σ, which
//! are obtained by fourth-order finite differences of the closed-form f with
//! steps relative to the natural scales ρ and ρ^{8/3}.

use std::f64::consts::PI;

use super::XcFunctional;

/// Densities below this value contribute nothing to E_xc or v_xc.
pub(super) const RHO_MIN: f64 = 1e-12;

/// (3/4)(3/π)^{1/3}: Slater exchange prefactor for the total density.
const SLATER_CX: f64 = 0.738_558_766_382_022_4;

/// Fraction of exact (Hartree-Fock) exchange in B3LYP.
const B3LYP_A0: f64 = 0.20;
/// Weight of the B88 gradient correction in B3LYP.
const B3LYP_AX: f64 = 0.72;
/// Weight of the LYP correlation in B3LYP.
const B3LYP_AC: f64 = 0.81;

/// VWN parameters (paramagnetic) of ε_c(x), x = √r_s, in Hartree.
#[derive(Debug, Clone, Copy)]
struct VwnParams {
    a: f64,
    x0: f64,
    b: f64,
    c: f64,
}

/// VWN5 (Monte-Carlo fit) parameters.
const VWN5: VwnParams = VwnParams {
    a: 0.031_090_7,
    x0: -0.104_98,
    b: 3.727_44,
    c: 12.935_2,
};

/// VWN RPA parameters (the VWN flavour of B3LYP).
const VWN_RPA: VwnParams = VwnParams {
    a: 0.031_090_7,
    x0: -0.409_286,
    b: 13.072_0,
    c: 42.719_8,
};

/// Energy densities (per volume) of the exchange and correlation parts.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct XcEnergyDensity {
    /// Semi-local exchange energy density \[Hartree/Bohr³\].
    pub(super) exchange: f64,
    /// Correlation energy density \[Hartree/Bohr³\].
    pub(super) correlation: f64,
}

/// Partial derivatives of the total semi-local energy density.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(super) struct XcDerivatives {
    /// ∂f/∂ρ \[Hartree\].
    pub(super) d_rho: f64,
    /// ∂f/∂σ \[Hartree·Bohr²\].
    pub(super) d_sigma: f64,
}

/// Whether the functional depends on the density gradient.
pub(super) fn is_gga(functional: &XcFunctional) -> bool {
    !matches!(functional, XcFunctional::Lda)
}

/// Fraction of exact exchange mixed into the functional.
pub(super) fn exact_exchange_fraction(functional: &XcFunctional) -> f64 {
    match functional {
        XcFunctional::B3lyp => B3LYP_A0,
        XcFunctional::Lda | XcFunctional::Pbe => 0.0,
    }
}

/// Exchange and correlation energy densities at density `rho` and squared
/// gradient `sigma`; zero below [`RHO_MIN`].
pub(super) fn energy_density(functional: &XcFunctional, rho: f64, sigma: f64) -> XcEnergyDensity {
    if rho.is_nan() || rho <= RHO_MIN {
        return XcEnergyDensity::default();
    }
    raw_energy_density(functional, rho, sigma.max(0.0))
}

/// ∂f/∂ρ and ∂f/∂σ of the total semi-local energy density; zero below [`RHO_MIN`].
pub(super) fn derivatives(functional: &XcFunctional, rho: f64, sigma: f64) -> XcDerivatives {
    if rho.is_nan() || rho <= RHO_MIN {
        return XcDerivatives::default();
    }
    match functional {
        XcFunctional::Lda => XcDerivatives {
            d_rho: slater_potential(rho) + vwn_potential(&VWN5, rho),
            d_sigma: 0.0,
        },
        XcFunctional::Pbe | XcFunctional::B3lyp => {
            let sigma = sigma.max(0.0);
            let total = |r: f64, s: f64| {
                let e = raw_energy_density(functional, r, s);
                e.exchange + e.correlation
            };
            let d_r = 1e-3 * rho;
            let d_rho = (total(rho - 2.0 * d_r, sigma) - 8.0 * total(rho - d_r, sigma)
                + 8.0 * total(rho + d_r, sigma)
                - total(rho + 2.0 * d_r, sigma))
                / (12.0 * d_r);
            let d_s = 1e-3 * sigma.max(rho.powf(8.0 / 3.0));
            let d_sigma = if sigma >= 2.0 * d_s {
                (total(rho, sigma - 2.0 * d_s) - 8.0 * total(rho, sigma - d_s)
                    + 8.0 * total(rho, sigma + d_s)
                    - total(rho, sigma + 2.0 * d_s))
                    / (12.0 * d_s)
            } else {
                (-25.0 * total(rho, sigma) + 48.0 * total(rho, sigma + d_s)
                    - 36.0 * total(rho, sigma + 2.0 * d_s)
                    + 16.0 * total(rho, sigma + 3.0 * d_s)
                    - 3.0 * total(rho, sigma + 4.0 * d_s))
                    / (12.0 * d_s)
            };
            XcDerivatives { d_rho, d_sigma }
        }
    }
}

fn raw_energy_density(functional: &XcFunctional, rho: f64, sigma: f64) -> XcEnergyDensity {
    match functional {
        XcFunctional::Lda => XcEnergyDensity {
            exchange: slater_exchange(rho),
            correlation: rho * vwn_epsilon(&VWN5, wigner_seitz_radius(rho)).0,
        },
        XcFunctional::Pbe => XcEnergyDensity {
            exchange: pbe_exchange(rho, sigma),
            correlation: pbe_correlation(rho, sigma),
        },
        XcFunctional::B3lyp => XcEnergyDensity {
            exchange: (1.0 - B3LYP_A0) * slater_exchange(rho)
                + B3LYP_AX * b88_gradient_correction(rho, sigma),
            correlation: (1.0 - B3LYP_AC) * rho * vwn_epsilon(&VWN_RPA, wigner_seitz_radius(rho)).0
                + B3LYP_AC * lyp_correlation(rho, sigma),
        },
    }
}

/// Wigner-Seitz radius r_s = (3/(4πρ))^{1/3}.
pub(super) fn wigner_seitz_radius(rho: f64) -> f64 {
    (3.0 / (4.0 * PI * rho)).cbrt()
}

/// Slater exchange energy density −(3/4)(3/π)^{1/3} ρ^{4/3}.
pub(super) fn slater_exchange(rho: f64) -> f64 {
    -SLATER_CX * rho.powf(4.0 / 3.0)
}

/// Slater exchange potential −(3/π)^{1/3} ρ^{1/3}.
fn slater_potential(rho: f64) -> f64 {
    -(4.0 / 3.0) * SLATER_CX * rho.cbrt()
}

/// VWN ε_c(x) and dε_c/dx for x = √r_s.
fn vwn_epsilon(p: &VwnParams, rs: f64) -> (f64, f64) {
    let x = rs.sqrt();
    let big_x = |t: f64| t * t + p.b * t + p.c;
    let q = (4.0 * p.c - p.b * p.b).sqrt();
    let xx = big_x(x);
    let xx0 = big_x(p.x0);
    let atan_term = (q / (2.0 * x + p.b)).atan();
    let eps = p.a
        * ((x * x / xx).ln() + 2.0 * p.b / q * atan_term
            - p.b * p.x0 / xx0
                * (((x - p.x0).powi(2) / xx).ln() + 2.0 * (p.b + 2.0 * p.x0) / q * atan_term));
    let denom = (2.0 * x + p.b).powi(2) + q * q;
    let deps_dx = p.a
        * (2.0 / x
            - (2.0 * x + p.b) / xx
            - 4.0 * p.b / denom
            - p.b * p.x0 / xx0
                * (2.0 / (x - p.x0) - (2.0 * x + p.b) / xx - 4.0 * (p.b + 2.0 * p.x0) / denom));
    (eps, deps_dx)
}

/// VWN correlation potential ε_c − (x/6) dε_c/dx.
fn vwn_potential(p: &VwnParams, rho: f64) -> f64 {
    let rs = wigner_seitz_radius(rho);
    let (eps, deps_dx) = vwn_epsilon(p, rs);
    eps - rs.sqrt() / 6.0 * deps_dx
}

/// PW92 unpolarised correlation energy per electron ε_c(r_s).
pub(super) fn pw92_epsilon(rs: f64) -> f64 {
    const A: f64 = 0.031_091;
    const ALPHA1: f64 = 0.213_70;
    const BETA1: f64 = 7.595_7;
    const BETA2: f64 = 3.587_6;
    const BETA3: f64 = 1.638_2;
    const BETA4: f64 = 0.492_94;
    let srs = rs.sqrt();
    let q1 = 2.0 * A * (BETA1 * srs + BETA2 * rs + BETA3 * rs * srs + BETA4 * rs * rs);
    -2.0 * A * (1.0 + ALPHA1 * rs) * (1.0 + 1.0 / q1).ln()
}

/// PBE exchange energy density ρ ε_x^unif F_x(s).
fn pbe_exchange(rho: f64, sigma: f64) -> f64 {
    const KAPPA: f64 = 0.804;
    const MU: f64 = 0.219_514_972_764_517_1;
    let kf = (3.0 * PI * PI * rho).cbrt();
    let s2 = sigma / (4.0 * kf * kf * rho * rho);
    let fx = 1.0 + KAPPA - KAPPA / (1.0 + MU * s2 / KAPPA);
    slater_exchange(rho) * fx
}

/// PBE correlation energy density ρ (ε_c^PW92 + H(r_s, t)).
fn pbe_correlation(rho: f64, sigma: f64) -> f64 {
    const BETA: f64 = 0.066_724_550_603_149_22;
    let gamma = (1.0 - std::f64::consts::LN_2) / (PI * PI);
    let eps_c = pw92_epsilon(wigner_seitz_radius(rho));
    let kf = (3.0 * PI * PI * rho).cbrt();
    let ks2 = 4.0 * kf / PI;
    let t2 = sigma / (4.0 * ks2 * rho * rho);
    let expm1 = (-eps_c / gamma).exp_m1();
    let h = if expm1 > 0.0 {
        let a = BETA / gamma / expm1;
        let at2 = a * t2;
        gamma * (1.0 + BETA / gamma * t2 * (1.0 + at2) / (1.0 + at2 + at2 * at2)).ln()
    } else {
        // ε_c → 0⁻: A → ∞ and H → γ ln(1 + (β/γ) t²·(1/(A t²))) → 0 for t > 0.
        0.0
    };
    rho * (eps_c + h)
}

/// Becke 88 gradient correction to the exchange energy density (closed shell).
///
/// ΔE_x^B88 = −β Σ_σ ρ_σ^{4/3} x_σ² / (1 + 6β x_σ sinh⁻¹ x_σ),
/// x_σ = |∇ρ_σ| / ρ_σ^{4/3}, β = 0.0042; ρ_σ = ρ/2, |∇ρ_σ|² = σ/4.
fn b88_gradient_correction(rho: f64, sigma: f64) -> f64 {
    const BETA: f64 = 0.004_2;
    let rho_s = 0.5 * rho;
    let rho_s43 = rho_s.powf(4.0 / 3.0);
    let x = (0.25 * sigma).sqrt() / rho_s43;
    -2.0 * BETA * rho_s43 * x * x / (1.0 + 6.0 * BETA * x * x.asinh())
}

/// LYP correlation energy density for a closed shell, in the gradient-only
/// form of Miehlich, Savin, Stoll & Preuss (CPL 157, 200 (1989)).
fn lyp_correlation(rho: f64, sigma: f64) -> f64 {
    const A: f64 = 0.049_18;
    const B: f64 = 0.132;
    const C: f64 = 0.253_3;
    const D: f64 = 0.349;
    let cf = 0.3 * (3.0 * PI * PI).powf(2.0 / 3.0);
    let rm13 = rho.powf(-1.0 / 3.0);
    let denom = 1.0 + D * rm13;
    let omega = (-C * rm13).exp() / denom * rho.powf(-11.0 / 3.0);
    let delta = C * rm13 + D * rm13 / denom;
    let rho_a = 0.5 * rho;
    let sigma_a = 0.25 * sigma;
    let first = -A * 4.0 / denom * rho_a * rho_a / rho;
    let bracket = rho_a
        * rho_a
        * (2.0_f64.powf(11.0 / 3.0) * cf * 2.0 * rho_a.powf(8.0 / 3.0)
            + (47.0 / 18.0 - 7.0 * delta / 18.0) * sigma
            - (2.5 - delta / 18.0) * 2.0 * sigma_a
            - (delta - 11.0) / 9.0 * (2.0 * rho_a / rho * sigma_a))
        - 2.0 / 3.0 * rho * rho * sigma
        + 2.0 * (2.0 / 3.0 * rho * rho - rho_a * rho_a) * sigma_a;
    first - A * B * omega * bracket
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ∫ 4πr² g(r) dr by the composite Simpson rule on [r_min, r_max].
    fn radial_integral(g: impl Fn(f64) -> f64, r_max: f64, n: usize) -> f64 {
        let h = r_max / n as f64;
        let mut sum = 0.0;
        for i in 0..=n {
            let r = (i as f64 * h).max(1e-9);
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            sum += w * 4.0 * PI * r * r * g(r);
        }
        sum * h / 3.0
    }

    #[test]
    fn test_slater_exchange_per_electron_at_rs_one() {
        // ε_x = −0.458165/r_s Hartree.
        let rho = 3.0 / (4.0 * PI);
        assert!((slater_exchange(rho) / rho + 0.458_165_293).abs() < 1e-8);
    }

    #[test]
    fn test_lda_correlation_matches_ceperley_alder_fits() {
        // VWN5 and PW92 both fit the Ceperley-Alder data and agree within 0.5 mHa.
        for rs in [0.5, 1.0, 2.0, 5.0, 10.0, 20.0] {
            let vwn = vwn_epsilon(&VWN5, rs).0;
            let pw = pw92_epsilon(rs);
            assert!((vwn - pw).abs() < 5e-4, "rs {rs}: {vwn} vs {pw}");
        }
        assert!((pw92_epsilon(2.0) - (-0.044_76)).abs() < 1e-4);
        // RPA overestimates correlation: about −0.0793 Ha at r_s = 1.
        assert!((vwn_epsilon(&VWN_RPA, 1.0).0 - (-0.079_31)).abs() < 1e-4);
    }

    #[test]
    fn test_vwn_derivative_is_consistent() {
        for rs in [0.3, 1.0, 4.0, 12.0] {
            let h = 1e-6;
            let x = f64::sqrt(rs);
            let fd = (vwn_epsilon(&VWN5, (x + h).powi(2)).0
                - vwn_epsilon(&VWN5, (x - h).powi(2)).0)
                / (2.0 * h);
            assert!((fd - vwn_epsilon(&VWN5, rs).1).abs() < 1e-8);
        }
    }

    #[test]
    fn test_potentials_are_density_derivatives() {
        for functional in [XcFunctional::Lda, XcFunctional::Pbe, XcFunctional::B3lyp] {
            for &(rho, sigma) in &[(0.3, 0.0), (0.05, 0.002), (1.7, 4.0), (1e-4, 1e-8)] {
                let d = derivatives(&functional, rho, sigma);
                let total = |r: f64, s: f64| {
                    let e = energy_density(&functional, r, s);
                    e.exchange + e.correlation
                };
                let hr = 1e-5 * rho;
                let fd_r = (total(rho + hr, sigma) - total(rho - hr, sigma)) / (2.0 * hr);
                assert!(
                    (fd_r - d.d_rho).abs() < 1e-6 * d.d_rho.abs().max(1e-3),
                    "{functional:?} rho {rho}: {fd_r} vs {}",
                    d.d_rho
                );
                if is_gga(&functional) {
                    let hs = 1e-4 * sigma.max(rho.powf(8.0 / 3.0));
                    let fd_s = if sigma >= hs {
                        (total(rho, sigma + hs) - total(rho, sigma - hs)) / (2.0 * hs)
                    } else {
                        // Second-order one-sided difference at σ = 0 (PBE's
                        // gradient terms cancel there, so f_σ ≈ 0).
                        (-3.0 * total(rho, sigma) + 4.0 * total(rho, sigma + hs)
                            - total(rho, sigma + 2.0 * hs))
                            / (2.0 * hs)
                    };
                    // Round-off of the difference quotient: ε |f| / hs.
                    let noise = 1e-13 * total(rho, sigma).abs() / hs;
                    assert!(
                        (fd_s - d.d_sigma).abs() < 1e-3 * d.d_sigma.abs() + noise,
                        "{functional:?} sigma: {fd_s} vs {}",
                        d.d_sigma
                    );
                } else {
                    assert_eq!(d.d_sigma, 0.0);
                }
            }
        }
    }

    #[test]
    fn test_gga_reduces_to_uniform_gas_without_gradient() {
        let rho = 0.2;
        let pbe = energy_density(&XcFunctional::Pbe, rho, 0.0);
        assert!((pbe.exchange - slater_exchange(rho)).abs() < 1e-14);
        let lda_c = rho * pw92_epsilon(wigner_seitz_radius(rho));
        assert!((pbe.correlation - lda_c).abs() < 1e-14);
        assert_eq!(b88_gradient_correction(rho, 0.0), 0.0);
    }

    #[test]
    fn test_exchange_functionals_for_hydrogenic_helium_density() {
        // ρ = 2ζ³/π e^{−2ζr} with ζ = 27/16: exact exchange is −5ζ/8.
        let zeta: f64 = 27.0 / 16.0;
        let rho = |r: f64| 2.0 * zeta.powi(3) / PI * (-2.0 * zeta * r).exp();
        let sigma = |r: f64| (2.0 * zeta * rho(r)).powi(2);
        let exact = -5.0 * zeta / 8.0;
        let lda = radial_integral(|r| slater_exchange(rho(r)), 25.0, 200_000);
        let b88 = radial_integral(
            |r| slater_exchange(rho(r)) + b88_gradient_correction(rho(r), sigma(r)),
            25.0,
            200_000,
        );
        let pbe = radial_integral(|r| pbe_exchange(rho(r), sigma(r)), 25.0, 200_000);
        // LDA underestimates |E_x| by ~14 %; B88 and PBE are within 1–2 %.
        assert!((lda - (-0.904_627)).abs() < 1e-4, "lda {lda}");
        assert!(((b88 - exact) / exact).abs() < 0.015, "b88 {b88}");
        assert!(((pbe - exact) / exact).abs() < 0.03, "pbe {pbe}");
        let lda_c = radial_integral(
            |r| rho(r) * pw92_epsilon(wigner_seitz_radius(rho(r))),
            25.0,
            200_000,
        );
        let pbe_c = radial_integral(|r| pbe_correlation(rho(r), sigma(r)), 25.0, 200_000);
        // LDA overbinds He correlation (≈ −0.11 Ha); PBE is close to the exact −0.042 Ha.
        assert!((lda_c - (-0.1141)).abs() < 1e-3, "lda_c {lda_c}");
        assert!((pbe_c - (-0.0448)).abs() < 1e-3, "pbe_c {pbe_c}");
    }

    #[test]
    fn test_lyp_gradient_form_matches_original_laplacian_form() {
        // The Miehlich form is obtained from the original LYP expression
        // (with ∇²ρ) by partial integration, so the integrals must agree.
        const A: f64 = 0.049_18;
        const B: f64 = 0.132;
        const C: f64 = 0.253_3;
        const D: f64 = 0.349;
        let cf = 0.3 * (3.0 * PI * PI).powf(2.0 / 3.0);
        for &(n_el, alpha) in &[(2.0, 1.0), (2.0, 3.0), (10.0, 0.7)] {
            let rho = |r: f64| n_el * (alpha / PI).powf(1.5) * (-alpha * r * r).exp();
            let sigma = |r: f64| (2.0 * alpha * r * rho(r)).powi(2);
            let lap = |r: f64| (4.0 * alpha * alpha * r * r - 6.0 * alpha) * rho(r);
            let original = |r: f64| {
                let p = rho(r);
                if p < 1e-30 {
                    return 0.0;
                }
                let rm13 = p.powf(-1.0 / 3.0);
                let t_w = sigma(r) / (8.0 * p) - lap(r) / 8.0;
                -A / (1.0 + D * rm13)
                    * (p + B
                        * p.powf(-2.0 / 3.0)
                        * (cf * p.powf(5.0 / 3.0) - 2.0 * t_w + (t_w / 9.0 + lap(r) / 18.0))
                        * (-C * rm13).exp())
            };
            let miehlich = |r: f64| {
                let p = rho(r);
                if p < 1e-30 {
                    0.0
                } else {
                    lyp_correlation(p, sigma(r))
                }
            };
            let e_m = radial_integral(miehlich, 15.0, 200_000);
            let e_o = radial_integral(original, 15.0, 200_000);
            assert!((e_m - e_o).abs() < 1e-9, "{e_m} vs {e_o}");
            assert!(e_m < 0.0);
        }
    }
}
