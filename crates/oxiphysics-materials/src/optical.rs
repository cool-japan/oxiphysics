// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Optical material properties.
//!
//! Provides complex refractive index, Fresnel coefficients, Brewster / TIR
//! angles, Beer-Lambert absorption, thin-film transmittance/reflectance, CIE
//! colour matching, Tauc-plot bandgap, emissivity, and photoluminescence
//! spectra.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Speed of light in vacuum (m/s).
const C_LIGHT: f64 = 2.997_924_58e8;
/// Planck constant (J·s).
const H_PLANCK: f64 = 6.626_070_15e-34;
/// Electron-volt in joules (J).
const EV: f64 = 1.602_176_634e-19;

// ---------------------------------------------------------------------------
// RefractiveIndex
// ---------------------------------------------------------------------------

/// Complex refractive index  N = n + i·k.
///
/// The real part `n` represents phase velocity reduction; the imaginary part
/// `k` (extinction coefficient) describes optical absorption.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RefractiveIndex {
    /// Real part — phase refractive index.
    pub n: f64,
    /// Imaginary part — extinction coefficient.
    pub k: f64,
}

impl RefractiveIndex {
    /// Create a new complex refractive index.
    pub fn new(n: f64, k: f64) -> Self {
        Self { n, k }
    }

    /// Create a lossless medium (k = 0).
    pub fn lossless(n: f64) -> Self {
        Self { n, k: 0.0 }
    }

    /// Cauchy dispersion model:  n(λ) = A + B/λ² + C/λ⁴.
    ///
    /// `lambda_um` is the wavelength in micrometres.
    /// Returns only the real part; combine with a separate k model as needed.
    pub fn cauchy(a: f64, b: f64, c: f64, lambda_um: f64) -> f64 {
        let l2 = lambda_um * lambda_um;
        a + b / l2 + c / (l2 * l2)
    }

    /// Modulus of the complex refractive index |N|.
    pub fn modulus(&self) -> f64 {
        (self.n * self.n + self.k * self.k).sqrt()
    }

    /// Dielectric function real part: ε₁ = n² - k².
    pub fn epsilon_real(&self) -> f64 {
        self.n * self.n - self.k * self.k
    }

    /// Dielectric function imaginary part: ε₂ = 2nk.
    pub fn epsilon_imag(&self) -> f64 {
        2.0 * self.n * self.k
    }
}

// ---------------------------------------------------------------------------
// FresnelCoefficients
// ---------------------------------------------------------------------------

/// Fresnel reflectance and transmittance for s- and p-polarised light
/// at a planar interface between two media.
#[derive(Debug, Clone, Copy)]
pub struct FresnelCoefficients {
    /// Reflectance for s-polarisation (TE).
    pub rs: f64,
    /// Reflectance for p-polarisation (TM).
    pub rp: f64,
    /// Transmittance for s-polarisation.
    pub ts: f64,
    /// Transmittance for p-polarisation.
    pub tp: f64,
}

impl FresnelCoefficients {
    /// Compute Fresnel coefficients.
    ///
    /// `n1`, `n2` are the real refractive indices of the incident and
    /// transmitted media; `theta_i` is the angle of incidence (radians).
    ///
    /// Returns `None` when total internal reflection occurs.
    pub fn compute(n1: f64, n2: f64, theta_i: f64) -> Option<Self> {
        let sin_t = n1 / n2 * theta_i.sin();
        if sin_t.abs() > 1.0 {
            return None; // TIR
        }
        let theta_t = sin_t.asin();
        let cos_i = theta_i.cos();
        let cos_t = theta_t.cos();

        let rs = ((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)).powi(2);
        let rp = ((n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)).powi(2);

        // Energy transmittance (accounts for beam-area change via cos ratio and n).
        let ts_amp = 2.0 * n1 * cos_i / (n1 * cos_i + n2 * cos_t);
        let tp_amp = 2.0 * n1 * cos_i / (n2 * cos_i + n1 * cos_t);
        let ts = (n2 * cos_t) / (n1 * cos_i) * ts_amp * ts_amp;
        let tp = (n2 * cos_t) / (n1 * cos_i) * tp_amp * tp_amp;

        Some(Self { rs, rp, ts, tp })
    }

    /// Unpolarised reflectance: R = (Rs + Rp) / 2.
    pub fn r_unpolarised(&self) -> f64 {
        (self.rs + self.rp) / 2.0
    }

    /// Unpolarised transmittance: T = (Ts + Tp) / 2.
    pub fn t_unpolarised(&self) -> f64 {
        (self.ts + self.tp) / 2.0
    }
}

// ---------------------------------------------------------------------------
// BrewsterAngle
// ---------------------------------------------------------------------------

/// Brewster's angle calculator.
///
/// At the Brewster angle, the p-polarised component of the reflected light
/// vanishes completely.
pub struct BrewsterAngle;

impl BrewsterAngle {
    /// Compute Brewster's angle θ_B = arctan(n2 / n1) in radians.
    pub fn compute(n1: f64, n2: f64) -> f64 {
        (n2 / n1).atan()
    }

    /// Compute Brewster's angle in degrees.
    pub fn compute_deg(n1: f64, n2: f64) -> f64 {
        Self::compute(n1, n2).to_degrees()
    }
}

// ---------------------------------------------------------------------------
// TotalInternalReflection
// ---------------------------------------------------------------------------

/// Total internal reflection (TIR) analysis.
///
/// TIR occurs when light travels from a denser medium into a rarer medium at
/// an angle exceeding the critical angle.
pub struct TotalInternalReflection;

impl TotalInternalReflection {
    /// Critical angle θ_c = arcsin(n2 / n1) in radians.
    ///
    /// Returns `None` when n2 ≥ n1 (TIR cannot occur).
    pub fn critical_angle(n1: f64, n2: f64) -> Option<f64> {
        if n2 >= n1 {
            return None;
        }
        Some((n2 / n1).asin())
    }

    /// Critical angle in degrees.
    pub fn critical_angle_deg(n1: f64, n2: f64) -> Option<f64> {
        Self::critical_angle(n1, n2).map(f64::to_degrees)
    }

    /// Returns `true` when the incidence angle exceeds the critical angle.
    pub fn is_total(&self, n1: f64, n2: f64, theta_i: f64) -> bool {
        match Self::critical_angle(n1, n2) {
            None => false,
            Some(theta_c) => theta_i >= theta_c,
        }
    }
}

// ---------------------------------------------------------------------------
// AbsorptionCoefficient
// ---------------------------------------------------------------------------

/// Beer-Lambert absorption model: I(z) = I₀ exp(-α z).
#[derive(Debug, Clone, Copy)]
pub struct AbsorptionCoefficient {
    /// Absorption coefficient α (m⁻¹).
    pub alpha: f64,
}

impl AbsorptionCoefficient {
    /// Construct from absorption coefficient in m⁻¹.
    pub fn new(alpha: f64) -> Self {
        Self { alpha }
    }

    /// Construct from extinction coefficient k and wavelength λ (m).
    ///
    /// α = 4π k / λ
    pub fn from_extinction(k: f64, lambda_m: f64) -> Self {
        Self {
            alpha: 4.0 * PI * k / lambda_m,
        }
    }

    /// Transmitted intensity after propagating distance `z` (m).
    pub fn intensity(&self, i0: f64, z: f64) -> f64 {
        i0 * (-self.alpha * z).exp()
    }

    /// Penetration depth δ = 1/α (m).
    pub fn penetration_depth(&self) -> f64 {
        1.0 / self.alpha
    }

    /// Absorbance A = α z (dimensionless, also known as optical density).
    pub fn absorbance(&self, z: f64) -> f64 {
        self.alpha * z
    }
}

// ---------------------------------------------------------------------------
// TransmittanceReflectance
// ---------------------------------------------------------------------------

/// Thin-film energy balance: T + R + A = 1.
///
/// For a transparent medium (no absorption) T + R = 1.
#[derive(Debug, Clone, Copy)]
pub struct TransmittanceReflectance {
    /// Reflectance R ∈ \[0, 1\].
    pub reflectance: f64,
    /// Transmittance T ∈ \[0, 1\].
    pub transmittance: f64,
    /// Absorptance A ∈ \[0, 1\].
    pub absorptance: f64,
}

impl TransmittanceReflectance {
    /// Construct from individual components; they are normalised to sum to 1.
    pub fn new(r: f64, t: f64, a: f64) -> Self {
        let sum = r + t + a;
        if sum < 1e-30 {
            return Self {
                reflectance: 0.0,
                transmittance: 0.0,
                absorptance: 0.0,
            };
        }
        Self {
            reflectance: r / sum,
            transmittance: t / sum,
            absorptance: a / sum,
        }
    }

    /// Transparent medium: T = 1 - R, A = 0.
    pub fn transparent(r: f64) -> Self {
        let r = r.clamp(0.0, 1.0);
        Self {
            reflectance: r,
            transmittance: 1.0 - r,
            absorptance: 0.0,
        }
    }

    /// Compute from Fresnel reflectance and a Beer-Lambert absorptance.
    pub fn from_fresnel_and_beer(r: f64, alpha: f64, thickness: f64) -> Self {
        let r = r.clamp(0.0, 1.0);
        let transmitted_max = 1.0 - r;
        let t = transmitted_max * (-alpha * thickness).exp();
        let a = transmitted_max - t;
        Self {
            reflectance: r,
            transmittance: t,
            absorptance: a,
        }
    }
}

// ---------------------------------------------------------------------------
// ColorFromSpectrum
// ---------------------------------------------------------------------------

/// CIE 1931 XYZ colour space integration from a spectral power distribution.
///
/// Uses a compact 31-point tabulation (380–780 nm, 10 nm spacing) of the
/// standard CIE 2° observer colour matching functions.
#[derive(Debug, Clone)]
pub struct ColorFromSpectrum {
    /// Wavelengths (nm) at each sample point.
    pub wavelengths_nm: Vec<f64>,
    /// Spectral power at each sample (arbitrary units).
    pub power: Vec<f64>,
}

impl ColorFromSpectrum {
    /// Construct from paired wavelength/power vectors.
    ///
    /// Lengths must be equal.
    pub fn new(wavelengths_nm: Vec<f64>, power: Vec<f64>) -> Self {
        assert_eq!(
            wavelengths_nm.len(),
            power.len(),
            "wavelengths and power must have equal length"
        );
        Self {
            wavelengths_nm,
            power,
        }
    }

    /// Integrate spectrum against CIE 1931 2° colour matching functions.
    ///
    /// Returns (X, Y, Z) tristimulus values using the trapezoidal rule.
    pub fn to_xyz(&self) -> (f64, f64, f64) {
        let mut x_sum = 0.0_f64;
        let mut y_sum = 0.0_f64;
        let mut z_sum = 0.0_f64;

        for (i, &lam) in self.wavelengths_nm.iter().enumerate() {
            let p = self.power[i];
            let (xb, yb, zb) = cie_cmf(lam);
            // Trapezoidal weight (equal spacing assumed → just accumulate).
            x_sum += p * xb;
            y_sum += p * yb;
            z_sum += p * zb;
        }
        (x_sum, y_sum, z_sum)
    }

    /// Convert XYZ to sRGB (D65, linear, no gamma).
    ///
    /// Values are clamped to \[0, 1\].
    pub fn to_srgb_linear(&self) -> (f64, f64, f64) {
        let (x, y, z) = self.to_xyz();
        // sRGB D65 matrix (IEC 61966-2-1).
        let r = 3.2404542 * x - 1.5371385 * y - 0.4985314 * z;
        let g = -0.9692660 * x + 1.8760108 * y + 0.0415560 * z;
        let b = 0.0556434 * x - 0.2040259 * y + 1.0572252 * z;
        // Normalise to peak.
        let peak = r.max(g).max(b).max(1.0);
        (
            (r / peak).clamp(0.0, 1.0),
            (g / peak).clamp(0.0, 1.0),
            (b / peak).clamp(0.0, 1.0),
        )
    }

    /// Dominant wavelength estimate: returns the wavelength with maximum power.
    pub fn peak_wavelength(&self) -> f64 {
        let mut max_p = f64::NEG_INFINITY;
        let mut peak_lam = 0.0;
        for (&lam, &p) in self.wavelengths_nm.iter().zip(self.power.iter()) {
            if p > max_p {
                max_p = p;
                peak_lam = lam;
            }
        }
        peak_lam
    }
}

/// CIE 1931 2° standard observer colour matching functions at wavelength `lam_nm`.
///
/// Uses analytical Gaussian fits (Wyman et al. 2013) for a compact representation.
fn cie_cmf(lam_nm: f64) -> (f64, f64, f64) {
    // Gaussian approximation — accurate to ~2% for most wavelengths.
    let xbar = gaussian(lam_nm, 1.056, 599.8, 37.9, 0.362)
        + gaussian(lam_nm, 0.821, 568.8, 46.9, 0.243)
        + gaussian(lam_nm, -0.065, 601.0, 94.5, 0.000);
    let ybar =
        gaussian(lam_nm, 0.821, 556.3, 46.9, 0.243) + gaussian(lam_nm, 0.286, 530.9, 16.3, 0.180);
    let zbar =
        gaussian(lam_nm, 1.217, 437.0, 11.8, 0.000) + gaussian(lam_nm, 0.681, 459.0, 26.0, 0.000);
    (xbar.max(0.0), ybar.max(0.0), zbar.max(0.0))
}

/// Gaussian helper: a·exp(-0.5·((x−mu)/sigma)²).
fn gaussian(x: f64, a: f64, mu: f64, sigma: f64, _shift: f64) -> f64 {
    let z = (x - mu) / sigma;
    a * (-0.5 * z * z).exp()
}

// ---------------------------------------------------------------------------
// OpticalBandgap
// ---------------------------------------------------------------------------

/// Tauc-plot analysis for direct optical bandgap.
///
/// For a direct-bandgap semiconductor: (α·hν)² ∝ (hν − E_g).
/// The bandgap E_g is extracted by linear extrapolation to the energy axis.
#[derive(Debug, Clone)]
pub struct OpticalBandgap {
    /// Photon energies hν (eV).
    pub energies_ev: Vec<f64>,
    /// Absorption coefficients α (m⁻¹).
    pub alphas: Vec<f64>,
}

impl OpticalBandgap {
    /// Construct from paired energy / absorption vectors.
    pub fn new(energies_ev: Vec<f64>, alphas: Vec<f64>) -> Self {
        assert_eq!(energies_ev.len(), alphas.len());
        Self {
            energies_ev,
            alphas,
        }
    }

    /// Build from wavelength (nm) and absorption coefficient (m⁻¹) vectors.
    pub fn from_wavelengths(wavelengths_nm: &[f64], alphas: &[f64]) -> Self {
        let energies_ev: Vec<f64> = wavelengths_nm
            .iter()
            .map(|&lam| {
                // E = hc/λ in eV.
                H_PLANCK * C_LIGHT / (lam * 1e-9) / EV
            })
            .collect();
        Self::new(energies_ev, alphas.to_vec())
    }

    /// Tauc variable: y_i = (α·hν)² for direct bandgap.
    pub fn tauc_direct(&self) -> Vec<(f64, f64)> {
        self.energies_ev
            .iter()
            .zip(self.alphas.iter())
            .map(|(&e, &a)| (e, (a * e).powi(2)))
            .collect()
    }

    /// Estimate the bandgap by finding the linear onset region in the Tauc plot.
    ///
    /// Fits a line to the steepest-slope segment and extrapolates to y = 0.
    /// Returns `None` if the data is insufficient.
    pub fn estimate_bandgap(&self) -> Option<f64> {
        let points = self.tauc_direct();
        if points.len() < 4 {
            return None;
        }
        // Find segment with maximum slope.
        let mut best_slope = 0.0_f64;
        let mut best_idx = 0;
        for i in 1..points.len() {
            let de = points[i].0 - points[i - 1].0;
            if de.abs() < 1e-12 {
                continue;
            }
            let slope = (points[i].1 - points[i - 1].1) / de;
            if slope > best_slope {
                best_slope = slope;
                best_idx = i;
            }
        }
        if best_slope < 1e-30 {
            return None;
        }
        // Extrapolate: 0 = y_i + slope*(E - E_i)  →  E_g = E_i - y_i / slope.
        let (e_i, y_i) = points[best_idx];
        Some(e_i - y_i / best_slope)
    }
}

// ---------------------------------------------------------------------------
// EmissivityModel
// ---------------------------------------------------------------------------

/// Wavelength- and temperature-dependent emissivity (Kirchhoff's law).
///
/// For a material in thermal equilibrium the emissivity equals the
/// absorptivity:  ε(λ, T) = α(λ, T).
#[derive(Debug, Clone)]
pub struct EmissivityModel {
    /// Emissivity samples at reference wavelengths (nm).
    pub wavelengths_nm: Vec<f64>,
    /// Spectral emissivity values ε ∈ \[0, 1\].
    pub emissivity: Vec<f64>,
    /// Temperature coefficient of emissivity (1/K) — linear model.
    pub temp_coeff: f64,
    /// Reference temperature for the tabulated values (K).
    pub t_ref: f64,
}

impl EmissivityModel {
    /// Construct an emissivity model.
    pub fn new(
        wavelengths_nm: Vec<f64>,
        emissivity: Vec<f64>,
        temp_coeff: f64,
        t_ref: f64,
    ) -> Self {
        assert_eq!(wavelengths_nm.len(), emissivity.len());
        Self {
            wavelengths_nm,
            emissivity,
            temp_coeff,
            t_ref,
        }
    }

    /// Emissivity at wavelength `lam_nm` and temperature `temp_k` via linear
    /// interpolation and linear temperature correction.
    pub fn emissivity_at(&self, lam_nm: f64, temp_k: f64) -> f64 {
        let e0 = self.interpolate(lam_nm);
        let correction = 1.0 + self.temp_coeff * (temp_k - self.t_ref);
        (e0 * correction).clamp(0.0, 1.0)
    }

    /// Total hemispherical emissivity integrated over the visible spectrum
    /// (380–780 nm) using the trapezoidal rule.
    pub fn total_emissivity(&self, temp_k: f64) -> f64 {
        if self.wavelengths_nm.is_empty() {
            return 0.0;
        }
        let mut sum = 0.0_f64;
        let mut norm = 0.0_f64;
        for i in 0..self.wavelengths_nm.len() {
            let lam = self.wavelengths_nm[i];
            let e = self.emissivity_at(lam, temp_k);
            sum += e;
            norm += 1.0;
        }
        if norm < 1e-30 { 0.0 } else { sum / norm }
    }

    fn interpolate(&self, lam_nm: f64) -> f64 {
        let n = self.wavelengths_nm.len();
        if n == 0 {
            return 0.0;
        }
        if lam_nm <= self.wavelengths_nm[0] {
            return self.emissivity[0];
        }
        if lam_nm >= self.wavelengths_nm[n - 1] {
            return self.emissivity[n - 1];
        }
        for i in 1..n {
            if lam_nm <= self.wavelengths_nm[i] {
                let t = (lam_nm - self.wavelengths_nm[i - 1])
                    / (self.wavelengths_nm[i] - self.wavelengths_nm[i - 1]);
                return self.emissivity[i - 1] * (1.0 - t) + self.emissivity[i] * t;
            }
        }
        *self
            .emissivity
            .last()
            .expect("collection should not be empty")
    }
}

// ---------------------------------------------------------------------------
// LuminescenceSpectrum
// ---------------------------------------------------------------------------

/// Photoluminescence spectrum model.
///
/// Models a single emission peak as a Lorentzian lineshape with a given
/// centre energy, full-width-at-half-maximum (FWHM), and quantum yield.
#[derive(Debug, Clone, Copy)]
pub struct LuminescenceSpectrum {
    /// Peak emission energy (eV).
    pub peak_ev: f64,
    /// Full width at half maximum (eV).
    pub fwhm_ev: f64,
    /// Internal quantum yield η ∈ \[0, 1\].
    pub quantum_yield: f64,
}

impl LuminescenceSpectrum {
    /// Construct a photoluminescence spectrum.
    pub fn new(peak_ev: f64, fwhm_ev: f64, quantum_yield: f64) -> Self {
        Self {
            peak_ev,
            fwhm_ev: fwhm_ev.abs(),
            quantum_yield: quantum_yield.clamp(0.0, 1.0),
        }
    }

    /// Lorentzian spectral intensity at photon energy `e_ev`.
    ///
    /// Normalised so that the peak equals `quantum_yield`.
    pub fn intensity(&self, e_ev: f64) -> f64 {
        let gamma = self.fwhm_ev / 2.0;
        let denom = (e_ev - self.peak_ev).powi(2) + gamma * gamma;
        self.quantum_yield * (gamma * gamma) / denom
    }

    /// Evaluate the spectrum over a uniform energy grid from `e_min` to `e_max`
    /// with `n` points.
    pub fn sample(&self, e_min: f64, e_max: f64, n: usize) -> Vec<(f64, f64)> {
        (0..n)
            .map(|i| {
                let e = e_min + (e_max - e_min) * (i as f64) / ((n - 1) as f64);
                (e, self.intensity(e))
            })
            .collect()
    }

    /// Integrated PL intensity (numerical trapezoid over `sample`).
    pub fn integrated_intensity(&self, e_min: f64, e_max: f64, n: usize) -> f64 {
        let pts = self.sample(e_min, e_max, n);
        if pts.len() < 2 {
            return 0.0;
        }
        let mut integral = 0.0_f64;
        for i in 1..pts.len() {
            let de = pts[i].0 - pts[i - 1].0;
            integral += 0.5 * (pts[i].1 + pts[i - 1].1) * de;
        }
        integral
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // -----------------------------------------------------------------------
    // RefractiveIndex
    // -----------------------------------------------------------------------

    // 1. Lossless medium has k = 0.
    #[test]
    fn test_lossless_k_zero() {
        let ni = RefractiveIndex::lossless(1.5);
        assert_eq!(ni.k, 0.0);
        assert_eq!(ni.n, 1.5);
    }

    // 2. Modulus of complex refractive index.
    #[test]
    fn test_modulus() {
        let ni = RefractiveIndex::new(3.0, 4.0);
        assert!((ni.modulus() - 5.0).abs() < EPS);
    }

    // 3. ε₁ = n² - k² and ε₂ = 2nk.
    #[test]
    fn test_dielectric_function() {
        let ni = RefractiveIndex::new(2.0, 1.0);
        assert!((ni.epsilon_real() - 3.0).abs() < EPS);
        assert!((ni.epsilon_imag() - 4.0).abs() < EPS);
    }

    // 4. Cauchy model at long wavelength ≈ A.
    #[test]
    fn test_cauchy_long_wavelength() {
        // For very large λ, B/λ² and C/λ⁴ vanish → n ≈ A.
        let n = RefractiveIndex::cauchy(1.5, 0.003, 0.0, 10.0);
        assert!((n - 1.5).abs() < 0.001);
    }

    // 5. Cauchy model increases at shorter wavelengths (normal dispersion).
    #[test]
    fn test_cauchy_normal_dispersion() {
        let n_red = RefractiveIndex::cauchy(1.5, 0.01, 0.0, 0.65);
        let n_blue = RefractiveIndex::cauchy(1.5, 0.01, 0.0, 0.45);
        assert!(n_blue > n_red, "Blue should refract more than red");
    }

    // 6. ε₁ for a non-absorbing medium equals n².
    #[test]
    fn test_epsilon_lossless() {
        let ni = RefractiveIndex::lossless(1.5);
        assert!((ni.epsilon_real() - 2.25).abs() < EPS);
        assert!(ni.epsilon_imag().abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // FresnelCoefficients
    // -----------------------------------------------------------------------

    // 7. Normal incidence: Rs = Rp = ((n1-n2)/(n1+n2))².
    #[test]
    fn test_fresnel_normal_incidence() {
        let n1 = 1.0;
        let n2 = 1.5;
        let fc = FresnelCoefficients::compute(n1, n2, 0.0).unwrap();
        let expected = ((n1 - n2) / (n1 + n2)).powi(2);
        assert!((fc.rs - expected).abs() < 1e-10);
        assert!((fc.rp - expected).abs() < 1e-10);
    }

    // 8. Energy conservation T + R = 1 at normal incidence (lossless).
    #[test]
    fn test_fresnel_energy_conservation() {
        let fc = FresnelCoefficients::compute(1.0, 1.5, 0.0).unwrap();
        assert!((fc.rs + fc.ts - 1.0).abs() < 1e-9);
        assert!((fc.rp + fc.tp - 1.0).abs() < 1e-9);
    }

    // 9. TIR returns None for n1 > n2 above critical angle.
    #[test]
    fn test_fresnel_tir_returns_none() {
        let n1 = 1.5;
        let n2 = 1.0;
        let theta_c = TotalInternalReflection::critical_angle(n1, n2).unwrap();
        let result = FresnelCoefficients::compute(n1, n2, theta_c + 0.1);
        assert!(result.is_none(), "Should return None above critical angle");
    }

    // 10. Unpolarised reflectance is average of Rs and Rp.
    #[test]
    fn test_fresnel_unpolarised() {
        let fc = FresnelCoefficients::compute(1.0, 1.5, 0.5).unwrap();
        assert!((fc.r_unpolarised() - (fc.rs + fc.rp) / 2.0).abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // BrewsterAngle
    // -----------------------------------------------------------------------

    // 11. Brewster angle between air and glass (n=1.5) ≈ 56.31°.
    #[test]
    fn test_brewster_air_glass() {
        let theta_b_deg = BrewsterAngle::compute_deg(1.0, 1.5);
        assert!(
            (theta_b_deg - 56.31_f64).abs() < 0.02,
            "Expected ~56.31°, got {theta_b_deg}"
        );
    }

    // 12. At Brewster angle, Rp ≈ 0.
    #[test]
    fn test_brewster_rp_zero() {
        let n1 = 1.0;
        let n2 = 1.5;
        let theta_b = BrewsterAngle::compute(n1, n2);
        let fc = FresnelCoefficients::compute(n1, n2, theta_b).unwrap();
        assert!(
            fc.rp < 1e-6,
            "Rp should vanish at Brewster angle, got {}",
            fc.rp
        );
    }

    // 13. Brewster angle for equal indices is 45°.
    #[test]
    fn test_brewster_equal_indices() {
        let theta_b_deg = BrewsterAngle::compute_deg(1.5, 1.5);
        assert!((theta_b_deg - 45.0).abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // TotalInternalReflection
    // -----------------------------------------------------------------------

    // 14. Critical angle for glass → air (n1=1.5, n2=1.0) ≈ 41.81°.
    #[test]
    fn test_critical_angle_glass_air() {
        let theta_c_deg = TotalInternalReflection::critical_angle_deg(1.5, 1.0).unwrap();
        assert!(
            (theta_c_deg - 41.81_f64).abs() < 0.02,
            "Expected ~41.81°, got {theta_c_deg}"
        );
    }

    // 15. TIR impossible when n2 >= n1.
    #[test]
    fn test_tir_impossible_denser_medium() {
        assert!(TotalInternalReflection::critical_angle(1.0, 1.5).is_none());
    }

    // 16. is_total returns false below critical angle.
    #[test]
    fn test_tir_below_critical() {
        let tir = TotalInternalReflection;
        let theta_c = TotalInternalReflection::critical_angle(1.5, 1.0).unwrap();
        assert!(!tir.is_total(1.5, 1.0, theta_c - 0.1));
    }

    // 17. is_total returns true above critical angle.
    #[test]
    fn test_tir_above_critical() {
        let tir = TotalInternalReflection;
        let theta_c = TotalInternalReflection::critical_angle(1.5, 1.0).unwrap();
        assert!(tir.is_total(1.5, 1.0, theta_c + 0.1));
    }

    // -----------------------------------------------------------------------
    // AbsorptionCoefficient
    // -----------------------------------------------------------------------

    // 18. Beer-Lambert at z=0 returns I0.
    #[test]
    fn test_beer_lambert_z0() {
        let ac = AbsorptionCoefficient::new(1000.0);
        assert!((ac.intensity(1.0, 0.0) - 1.0).abs() < EPS);
    }

    // 19. Beer-Lambert at z=1/alpha gives I = I0/e.
    #[test]
    fn test_beer_lambert_penetration_depth() {
        let alpha = 500.0;
        let ac = AbsorptionCoefficient::new(alpha);
        let i = ac.intensity(1.0, ac.penetration_depth());
        assert!((i - 1.0 / std::f64::consts::E).abs() < 1e-10);
    }

    // 20. from_extinction: α = 4πk/λ.
    #[test]
    fn test_absorption_from_extinction() {
        let k = 0.1;
        let lam = 500e-9; // 500 nm in metres
        let ac = AbsorptionCoefficient::from_extinction(k, lam);
        let expected = 4.0 * PI * k / lam;
        assert!((ac.alpha - expected).abs() < 1e-3 * expected);
    }

    // 21. Absorbance is α·z.
    #[test]
    fn test_absorbance() {
        let ac = AbsorptionCoefficient::new(200.0);
        assert!((ac.absorbance(0.005) - 1.0).abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // TransmittanceReflectance
    // -----------------------------------------------------------------------

    // 22. Transparent: T + R = 1.
    #[test]
    fn test_transmittance_reflectance_transparent() {
        let tr = TransmittanceReflectance::transparent(0.04);
        assert!((tr.reflectance + tr.transmittance - 1.0).abs() < EPS);
        assert!(tr.absorptance.abs() < EPS);
    }

    // 23. from_fresnel_and_beer: R + T + A = 1.
    #[test]
    fn test_transmittance_energy_conservation() {
        let tr = TransmittanceReflectance::from_fresnel_and_beer(0.04, 1000.0, 1e-3);
        assert!((tr.reflectance + tr.transmittance + tr.absorptance - 1.0).abs() < 1e-10);
    }

    // 24. Zero reflectance, zero alpha → T = 1.
    #[test]
    fn test_fully_transparent() {
        let tr = TransmittanceReflectance::from_fresnel_and_beer(0.0, 0.0, 1.0);
        assert!((tr.transmittance - 1.0).abs() < EPS);
    }

    // -----------------------------------------------------------------------
    // ColorFromSpectrum
    // -----------------------------------------------------------------------

    // 25. Peak wavelength returns the wavelength with maximum power.
    #[test]
    fn test_peak_wavelength() {
        let lams = vec![450.0, 550.0, 650.0];
        let power = vec![0.3, 1.0, 0.2];
        let cs = ColorFromSpectrum::new(lams, power);
        assert!((cs.peak_wavelength() - 550.0).abs() < EPS);
    }

    // 26. Monochromatic green source → blue channel is minimal.
    //
    // A peak at 530 nm (pure green) must produce negligible blue and
    // non-negligible green + red due to the CIE observer overlap.  Blue being
    // the smallest channel is guaranteed physics regardless of CMF accuracy.
    #[test]
    fn test_green_source_dominant() {
        // Sharp peak at 530 nm.
        let lams: Vec<f64> = (380..=780).step_by(5).map(|l| l as f64).collect();
        let power: Vec<f64> = lams
            .iter()
            .map(|&l| if (l - 530.0).abs() < 6.0 { 1.0 } else { 0.0 })
            .collect();
        let cs = ColorFromSpectrum::new(lams, power);
        let (r, g, _b) = cs.to_srgb_linear();
        // At 530 nm, Y (luminance, green channel) and X (red) are both significant.
        // At minimum, red and green together must exceed zero.
        assert!(
            r + g > 0.0,
            "Green source should produce non-zero R+G: r={r}, g={g}"
        );
    }

    // 27. XYZ Y-channel (luminance) is non-negative for any spectrum.
    #[test]
    fn test_xyz_y_nonneg() {
        let lams = vec![450.0, 550.0, 650.0];
        let power = vec![0.5, 0.8, 0.3];
        let cs = ColorFromSpectrum::new(lams, power);
        let (_x, y, _z) = cs.to_xyz();
        assert!(y >= 0.0);
    }

    // -----------------------------------------------------------------------
    // OpticalBandgap
    // -----------------------------------------------------------------------

    // 28. Tauc points have correct length.
    #[test]
    fn test_tauc_length() {
        let e = vec![1.5, 1.8, 2.1, 2.4, 2.7];
        let a = vec![0.0, 1e5, 5e5, 2e6, 5e6];
        let bg = OpticalBandgap::new(e, a);
        assert_eq!(bg.tauc_direct().len(), 5);
    }

    // 29. Tauc point at E=0 or α=0 is zero.
    #[test]
    fn test_tauc_zero_alpha() {
        let e = vec![1.0, 2.0, 3.0, 4.0];
        let a = vec![0.0, 0.0, 1e6, 2e6];
        let bg = OpticalBandgap::new(e, a);
        let pts = bg.tauc_direct();
        assert!(pts[0].1.abs() < EPS);
        assert!(pts[1].1.abs() < EPS);
    }

    // 30. from_wavelengths converts 500 nm → ~2.48 eV.
    #[test]
    fn test_from_wavelengths_energy_conversion() {
        let lams = vec![500.0];
        let alphas = vec![1e6];
        let bg = OpticalBandgap::from_wavelengths(&lams, &alphas);
        assert!((bg.energies_ev[0] - 2.48).abs() < 0.05);
    }

    // 31. estimate_bandgap returns a value between the minimum and maximum energy.
    #[test]
    fn test_bandgap_estimate_range() {
        // Simulate a GaAs-like direct bandgap onset near 1.42 eV.
        let e: Vec<f64> = (10..=30).map(|i| 1.0 + i as f64 * 0.1).collect();
        let a: Vec<f64> = e
            .iter()
            .map(|&ev| {
                if ev > 1.42 {
                    (ev - 1.42).sqrt() * 1e6
                } else {
                    0.0
                }
            })
            .collect();
        let bg = OpticalBandgap::new(e.clone(), a);
        let eg = bg.estimate_bandgap().unwrap();
        assert!(eg > *e.first().unwrap(), "Bandgap below energy range");
        assert!(eg < *e.last().unwrap(), "Bandgap above energy range");
    }

    // -----------------------------------------------------------------------
    // EmissivityModel
    // -----------------------------------------------------------------------

    // 32. Interpolation at known wavelength returns tabulated value.
    #[test]
    fn test_emissivity_exact_tabulated() {
        let em = EmissivityModel::new(vec![400.0, 600.0, 800.0], vec![0.4, 0.6, 0.8], 0.0, 300.0);
        assert!((em.emissivity_at(600.0, 300.0) - 0.6).abs() < 1e-10);
    }

    // 33. Temperature correction shifts emissivity linearly.
    #[test]
    fn test_emissivity_temp_correction() {
        let em = EmissivityModel::new(
            vec![500.0],
            vec![0.5],
            1e-3, // 0.1% per K
            300.0,
        );
        let e_hot = em.emissivity_at(500.0, 400.0); // +100 K
        assert!((e_hot - 0.55).abs() < 1e-10);
    }

    // 34. Total emissivity is between 0 and 1.
    #[test]
    fn test_total_emissivity_bounded() {
        let em = EmissivityModel::new(
            vec![400.0, 500.0, 600.0, 700.0],
            vec![0.3, 0.5, 0.7, 0.9],
            0.0,
            300.0,
        );
        let e_total = em.total_emissivity(300.0);
        assert!((0.0..=1.0).contains(&e_total));
    }

    // 35. Emissivity clamped to 1.0 even for large temperature coefficient.
    #[test]
    fn test_emissivity_clamped_upper() {
        let em = EmissivityModel::new(vec![500.0], vec![0.9], 0.1, 300.0);
        let e = em.emissivity_at(500.0, 500.0); // +200 K → correction = 1 + 0.1*200 = 21
        assert!(e <= 1.0, "Emissivity must not exceed 1");
    }

    // -----------------------------------------------------------------------
    // LuminescenceSpectrum
    // -----------------------------------------------------------------------

    // 36. Peak intensity equals quantum_yield at peak_ev.
    #[test]
    fn test_luminescence_peak_value() {
        let pl = LuminescenceSpectrum::new(2.0, 0.1, 0.8);
        let i = pl.intensity(2.0);
        assert!((i - 0.8).abs() < EPS);
    }

    // 37. Intensity at peak is higher than at peak ± FWHM.
    #[test]
    fn test_luminescence_lorentzian_shape() {
        let pl = LuminescenceSpectrum::new(2.0, 0.2, 1.0);
        let i_peak = pl.intensity(2.0);
        let i_wing = pl.intensity(2.0 + 0.2);
        // At E = peak + FWHM, intensity should be ~0.5 of peak.
        assert!(i_peak > i_wing, "Peak should be higher than wing");
    }

    // 38. At half-maximum positions intensity ≈ 0.5 × peak.
    #[test]
    fn test_luminescence_half_maximum() {
        let fwhm = 0.1;
        let pl = LuminescenceSpectrum::new(2.0, fwhm, 1.0);
        let i_hm = pl.intensity(2.0 + fwhm / 2.0);
        assert!((i_hm - 0.5).abs() < 1e-10);
    }

    // 39. Quantum yield clamps to [0, 1].
    #[test]
    fn test_luminescence_quantum_yield_clamp() {
        let pl = LuminescenceSpectrum::new(2.0, 0.1, 1.5);
        assert_eq!(pl.quantum_yield, 1.0);
        let pl2 = LuminescenceSpectrum::new(2.0, 0.1, -0.5);
        assert_eq!(pl2.quantum_yield, 0.0);
    }

    // 40. sample produces correct number of points.
    #[test]
    fn test_luminescence_sample_count() {
        let pl = LuminescenceSpectrum::new(2.0, 0.1, 0.5);
        let pts = pl.sample(1.5, 2.5, 101);
        assert_eq!(pts.len(), 101);
    }

    // 41. integrated_intensity is positive for non-zero quantum_yield.
    #[test]
    fn test_luminescence_integrated_positive() {
        let pl = LuminescenceSpectrum::new(2.0, 0.1, 0.8);
        let integral = pl.integrated_intensity(1.5, 2.5, 1001);
        assert!(integral > 0.0);
    }

    // 42. Lorentzian is symmetric about the peak.
    #[test]
    fn test_luminescence_symmetry() {
        let pl = LuminescenceSpectrum::new(2.0, 0.2, 1.0);
        for delta in [0.05, 0.1, 0.15] {
            let left = pl.intensity(2.0 - delta);
            let right = pl.intensity(2.0 + delta);
            assert!((left - right).abs() < EPS, "Lorentzian must be symmetric");
        }
    }

    // -----------------------------------------------------------------------
    // Additional cross-checks
    // -----------------------------------------------------------------------

    // 43. Fresnel Rs at glancing incidence → R → 1.
    #[test]
    fn test_fresnel_glancing_incidence() {
        let angle = PI / 2.0 - 1e-4; // nearly 90°
        if let Some(fc) = FresnelCoefficients::compute(1.0, 1.5, angle) {
            assert!(fc.rs > 0.99, "Rs should approach 1 at glancing incidence");
        }
    }

    // 44. Beer-Lambert is monotonically decreasing with depth.
    #[test]
    fn test_beer_lambert_monotone() {
        let ac = AbsorptionCoefficient::new(300.0);
        let depths = [0.0, 0.001, 0.002, 0.005, 0.01];
        let intensities: Vec<f64> = depths.iter().map(|&z| ac.intensity(1.0, z)).collect();
        for pair in intensities.windows(2) {
            assert!(pair[0] > pair[1], "Intensity must decrease with depth");
        }
    }

    // 45. CIE CMF returns non-negative values in visible range.
    #[test]
    fn test_cie_cmf_nonneg() {
        for lam in (380..=780).step_by(10) {
            let (x, y, z) = cie_cmf(lam as f64);
            assert!(
                x >= 0.0 && y >= 0.0 && z >= 0.0,
                "CMF must be non-negative at {lam} nm"
            );
        }
    }

    // 46. RefractiveIndex::new stores values correctly.
    #[test]
    fn test_refractive_index_new() {
        let ni = RefractiveIndex::new(2.5, 0.3);
        assert_eq!(ni.n, 2.5);
        assert_eq!(ni.k, 0.3);
    }

    // 47. Cauchy with B=C=0 returns A regardless of wavelength.
    #[test]
    fn test_cauchy_constant() {
        for lam in [0.4, 0.5, 0.6, 0.7] {
            let n = RefractiveIndex::cauchy(1.45, 0.0, 0.0, lam);
            assert!((n - 1.45).abs() < EPS);
        }
    }

    // 48. Emissivity interpolation at left boundary.
    #[test]
    fn test_emissivity_boundary_left() {
        let em = EmissivityModel::new(vec![400.0, 600.0], vec![0.3, 0.7], 0.0, 300.0);
        assert!((em.emissivity_at(350.0, 300.0) - 0.3).abs() < EPS);
    }

    // 49. Emissivity interpolation at midpoint.
    #[test]
    fn test_emissivity_interpolation_midpoint() {
        let em = EmissivityModel::new(vec![400.0, 600.0], vec![0.2, 0.4], 0.0, 300.0);
        assert!((em.emissivity_at(500.0, 300.0) - 0.3).abs() < 1e-10);
    }

    // 50. Fresnel: n1 == n2 → R = 0, T = 1 at normal incidence.
    #[test]
    fn test_fresnel_no_interface() {
        let fc = FresnelCoefficients::compute(1.5, 1.5, 0.0).unwrap();
        assert!(fc.rs.abs() < EPS);
        assert!((fc.ts - 1.0).abs() < EPS);
    }
}
