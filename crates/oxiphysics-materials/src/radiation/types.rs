//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    BOLTZMANN_K, SIGMA, blackbody_emissive_power, blackbody_spectral_intensity, wien_displacement,
};

/// Simple single-diode solar cell model.
///
/// Computes the current–voltage (I–V) characteristic of an ideal single-diode
/// solar cell using the Shockley equation:
///
/// `I = I_ph - I_0 * (exp(q*V/(n*k*T)) - 1)`
///
/// where `I_ph` is the photocurrent and `I_0` is the dark saturation current.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct SolarCell {
    /// Photocurrent \[A\]
    pub i_ph: f64,
    /// Dark saturation current \[A\]
    pub i_0: f64,
    /// Ideality factor (1 for ideal)
    pub n_ideal: f64,
    /// Temperature \[K\]
    pub temperature: f64,
    /// Series resistance \[Ω\]
    pub r_series: f64,
}
impl SolarCell {
    /// Create a new solar cell model.
    #[allow(dead_code)]
    pub fn new(i_ph: f64, i_0: f64, n_ideal: f64, temperature: f64, r_series: f64) -> Self {
        Self {
            i_ph,
            i_0,
            n_ideal,
            temperature,
            r_series,
        }
    }
    /// Thermal voltage \[V\]: `V_t = k*T/q`.
    #[allow(dead_code)]
    pub fn thermal_voltage(&self) -> f64 {
        BOLTZMANN_K * self.temperature / 1.602e-19
    }
    /// Short-circuit current \[A\] (at V = 0): `I_sc ≈ I_ph`.
    #[allow(dead_code)]
    pub fn short_circuit_current(&self) -> f64 {
        self.i_ph - self.i_0 * (1.0 / (self.n_ideal * self.thermal_voltage())).exp()
    }
    /// Current \[A\] at terminal voltage `v` \[V\].
    ///
    /// Neglects series resistance for simplicity (direct explicit form).
    #[allow(dead_code)]
    pub fn current_at_voltage(&self, v: f64) -> f64 {
        let v_t = self.thermal_voltage();
        (self.i_ph - self.i_0 * ((v / (self.n_ideal * v_t)).exp() - 1.0)).max(0.0)
    }
    /// Power \[W\] at terminal voltage `v` \[V\].
    #[allow(dead_code)]
    pub fn power_at_voltage(&self, v: f64) -> f64 {
        v * self.current_at_voltage(v)
    }
    /// Open-circuit voltage \[V\]: solved numerically (bisection).
    ///
    /// At V_oc: I = 0 → `I_ph = I_0*(exp(V_oc/(n*V_t)) - 1)`.
    /// Exact: `V_oc = n*V_t * ln(I_ph/I_0 + 1)`.
    #[allow(dead_code)]
    pub fn open_circuit_voltage(&self) -> f64 {
        if self.i_0 < f64::EPSILON {
            return 0.0;
        }
        self.n_ideal * self.thermal_voltage() * (self.i_ph / self.i_0 + 1.0).ln()
    }
    /// Fill factor (FF): `FF = P_max / (I_sc * V_oc)`.
    ///
    /// Estimated using the empirical Green formula:
    /// `FF ≈ (v_oc - ln(v_oc + 0.72)) / (v_oc + 1)`
    /// where `v_oc = V_oc / V_t` (normalised open-circuit voltage).
    #[allow(dead_code)]
    pub fn fill_factor(&self) -> f64 {
        let v_t = self.thermal_voltage();
        let v_oc = self.open_circuit_voltage();
        let v_oc_norm = v_oc / (self.n_ideal * v_t);
        if v_oc_norm < 1.0 {
            return 0.0;
        }
        (v_oc_norm - (v_oc_norm + 0.72).ln()) / (v_oc_norm + 1.0)
    }
    /// Efficiency: `η = P_max / P_in` where `P_max = FF * I_sc * V_oc`.
    ///
    /// # Arguments
    /// * `irradiance` — incident irradiance \[W/m²\]
    /// * `area`       — cell area \[m²\]
    #[allow(dead_code)]
    pub fn efficiency(&self, irradiance: f64, area: f64) -> f64 {
        let p_in = irradiance * area;
        if p_in < f64::EPSILON {
            return 0.0;
        }
        let i_sc = self.short_circuit_current();
        let v_oc = self.open_circuit_voltage();
        let ff = self.fill_factor();
        ff * i_sc * v_oc / p_in
    }
}
/// A gray, diffuse radiation surface for network analysis.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct RadiationSurface {
    /// Surface area \[m²\]
    pub area: f64,
    /// Surface emissivity (0–1)
    pub emissivity: f64,
    /// Surface temperature \[K\]
    pub temperature: f64,
    /// Descriptive name
    pub name: String,
}
impl RadiationSurface {
    /// Create a new radiation surface.
    #[allow(dead_code)]
    pub fn new(area: f64, emissivity: f64, temp_k: f64, name: &str) -> Self {
        Self {
            area,
            emissivity,
            temperature: temp_k,
            name: name.to_string(),
        }
    }
    /// Radiosity \[W/m²\]: simplified as J = epsilon * sigma * T^4.
    #[allow(dead_code)]
    pub fn radiosity(&self) -> f64 {
        self.emissivity * SIGMA * self.temperature.powi(4)
    }
    /// Surface (blackbody) resistance \[(W/m²)^{-1}\]: (1 - epsilon) / (epsilon * A).
    ///
    /// Returns infinity for a blackbody (epsilon = 1).
    #[allow(dead_code)]
    pub fn surface_resistance(&self) -> f64 {
        if (self.emissivity - 1.0).abs() < 1e-12 {
            0.0
        } else {
            (1.0 - self.emissivity) / (self.emissivity * self.area)
        }
    }
}
/// Neutron moderation properties for a moderator material.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct NeutronModeration {
    /// Macroscopic scattering cross-section Σ_s \[cm⁻¹\]
    pub sigma_s: f64,
    /// Macroscopic absorption cross-section Σ_a \[cm⁻¹\]
    pub sigma_a: f64,
    /// Average logarithmic energy decrement ξ (dimensionless)
    pub xi: f64,
}
impl NeutronModeration {
    /// Create a new `NeutronModeration` model.
    #[allow(dead_code)]
    pub fn new(sigma_s: f64, sigma_a: f64, xi: f64) -> Self {
        Self {
            sigma_s,
            sigma_a,
            xi,
        }
    }
    /// Slowing-down power \[cm⁻¹\]: `SDP = ξ · Σ_s`.
    #[allow(dead_code)]
    pub fn slowing_down_power(&self) -> f64 {
        self.xi * self.sigma_s
    }
    /// Moderation ratio (figure of merit): `MR = ξ · Σ_s / Σ_a`.
    #[allow(dead_code)]
    pub fn moderation_ratio(&self) -> f64 {
        if self.sigma_a < f64::EPSILON {
            return f64::INFINITY;
        }
        self.xi * self.sigma_s / self.sigma_a
    }
    /// Migration length squared \[cm²\]: `M² = D / Σ_a` where `D = 1/(3·Σ_s)`.
    #[allow(dead_code)]
    pub fn migration_length_sq(&self) -> f64 {
        let d = 1.0 / (3.0 * self.sigma_s);
        d / self.sigma_a
    }
}
/// Blackbody spectrum utilities.
///
/// Groups the Planck function, Stefan-Boltzmann total power, and Wien peak
/// for a surface at temperature `temp_k`.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct BlackbodySpectrum {
    /// Surface temperature \[K\]
    pub temp_k: f64,
}
impl BlackbodySpectrum {
    /// Create a new `BlackbodySpectrum` at temperature `temp_k`.
    #[allow(dead_code)]
    pub fn new(temp_k: f64) -> Self {
        Self { temp_k }
    }
    /// Spectral radiance `B(λ, T)` \[W/(m²·sr·m)\] at wavelength `lambda_m` \[m\].
    #[allow(dead_code)]
    pub fn planck(&self, lambda_m: f64) -> f64 {
        blackbody_spectral_intensity(lambda_m, self.temp_k)
    }
    /// Total emissive power `E_b = σ·T⁴` \[W/m²\].
    #[allow(dead_code)]
    pub fn total_power(&self) -> f64 {
        blackbody_emissive_power(self.temp_k)
    }
    /// Peak wavelength via Wien's displacement law \[m\].
    #[allow(dead_code)]
    pub fn peak_wavelength(&self) -> f64 {
        wien_displacement(self.temp_k)
    }
    /// Fractional emissive power in wavelength range \[lambda_a, lambda_b\] via
    /// a simple trapezoidal quadrature over `n_steps` intervals.
    #[allow(dead_code)]
    pub fn band_fraction(&self, lambda_a: f64, lambda_b: f64, n_steps: usize) -> f64 {
        let eb = self.total_power();
        if eb < f64::EPSILON {
            return 0.0;
        }
        let n = n_steps.max(2);
        let dlam = (lambda_b - lambda_a) / (n as f64);
        let mut sum = 0.0;
        for i in 0..=n {
            let lam = lambda_a + (i as f64) * dlam;
            let w = if i == 0 || i == n { 0.5 } else { 1.0 };
            sum += w * std::f64::consts::PI * self.planck(lam);
        }
        (sum * dlam / eb).clamp(0.0, 1.0)
    }
}
/// Simple Monte Carlo ray tracer for view factor estimation.
///
/// Estimates the view factor F12 between two planar surfaces by firing
/// random rays from surface 1 and counting how many hit surface 2.
///
/// For reproducible tests a seeded "pseudo-random" deterministic sequence
/// is used (simple LCG) so no external rand dependency is needed here.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct MonteCarloViewFactor {
    /// Number of rays per calculation.
    pub n_rays: usize,
    /// LCG seed for reproducibility.
    pub(super) seed: u64,
}
impl MonteCarloViewFactor {
    /// Create a new Monte Carlo view factor calculator.
    pub fn new(n_rays: usize) -> Self {
        Self {
            n_rays,
            seed: 12345678901,
        }
    }
    /// Simple 64-bit LCG random number in \[0, 1).
    fn rand_next(&mut self) -> f64 {
        self.seed = self
            .seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.seed >> 33) as f64 / (u32::MAX as f64)
    }
    /// Estimate the view factor for two coaxial parallel disks using MC.
    ///
    /// Uses the hemisphere cosine sampling method.  Returns an estimate
    /// of F12 (fraction of rays from disk 1 that hit disk 2).
    ///
    /// # Arguments
    /// * `r1` — radius of disk 1
    /// * `r2` — radius of disk 2
    /// * `h`  — separation distance
    pub fn parallel_disks(&mut self, r1: f64, r2: f64, h: f64) -> f64 {
        let mut hits = 0usize;
        for _ in 0..self.n_rays {
            let r_orig = r1 * self.rand_next().sqrt();
            let phi_orig = 2.0 * std::f64::consts::PI * self.rand_next();
            let ox = r_orig * phi_orig.cos();
            let oy = r_orig * phi_orig.sin();
            let theta = (self.rand_next()).sqrt().asin();
            let phi_dir = 2.0 * std::f64::consts::PI * self.rand_next();
            let dz = theta.cos();
            let dx = theta.sin() * phi_dir.cos();
            let dy = theta.sin() * phi_dir.sin();
            if dz < 1e-12 {
                continue;
            }
            let t = h / dz;
            let ix = ox + dx * t;
            let iy = oy + dy * t;
            if ix * ix + iy * iy <= r2 * r2 {
                hits += 1;
            }
        }
        hits as f64 / self.n_rays as f64
    }
}
/// N-surface radiation network (gray, diffuse enclosure).
///
/// Solves the full N×N radiosity system for a closed diffuse gray enclosure.
/// The net heat flow on each surface is computed from the radiosity vector.
///
/// Reference: Incropera et al., "Fundamentals of Heat and Mass Transfer".
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RadiationNetwork {
    /// Number of surfaces.
    pub n: usize,
    /// Surface temperatures \[K\].
    pub temperatures: Vec<f64>,
    /// Surface emissivities (0–1).
    pub emissivities: Vec<f64>,
    /// Surface areas \[m²\].
    pub areas: Vec<f64>,
    /// View factor matrix F\[i\]\[j\] (row-major, n×n).
    pub view_factors: Vec<f64>,
}
impl RadiationNetwork {
    /// Create a new radiation network.
    pub fn new(
        temperatures: Vec<f64>,
        emissivities: Vec<f64>,
        areas: Vec<f64>,
        view_factors: Vec<f64>,
    ) -> Self {
        let n = temperatures.len();
        assert_eq!(emissivities.len(), n);
        assert_eq!(areas.len(), n);
        assert_eq!(view_factors.len(), n * n);
        Self {
            n,
            temperatures,
            emissivities,
            areas,
            view_factors,
        }
    }
    /// View factor F\[i\]\[j\].
    pub fn f(&self, i: usize, j: usize) -> f64 {
        self.view_factors[i * self.n + j]
    }
    /// Net heat flow from surface i \[W\] using the radiosity method.
    ///
    /// Simplified: q_i = ε_i·A_i·σ·(T_i⁴) - A_i·sum_j(F_ij · ε_j·σ·T_j⁴)
    ///
    /// This is the optically-simplified (no inter-reflection) version.
    pub fn net_heat_flow_simple(&self, i: usize) -> f64 {
        let emit_i = self.emissivities[i] * SIGMA * self.temperatures[i].powi(4) * self.areas[i];
        let absorbed: f64 = (0..self.n)
            .map(|j| {
                self.f(i, j)
                    * self.emissivities[j]
                    * SIGMA
                    * self.temperatures[j].powi(4)
                    * self.areas[i]
            })
            .sum();
        emit_i - absorbed
    }
    /// Total radiative power emitted by surface i \[W\].
    pub fn emitted_power(&self, i: usize) -> f64 {
        self.emissivities[i] * SIGMA * self.temperatures[i].powi(4) * self.areas[i]
    }
}
/// Spectrally resolved emissivity model using piecewise linear interpolation.
///
/// Allows emissivity to vary with wavelength, which is important for
/// selective emitters/absorbers in solar thermal and thermophotovoltaic systems.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct SpectralEmissivityModel {
    /// Wavelengths \[m\] at which emissivity is specified (sorted ascending).
    pub wavelengths: Vec<f64>,
    /// Emissivity values (0–1) at each wavelength knot.
    pub emissivities: Vec<f64>,
}
impl SpectralEmissivityModel {
    /// Create from matched wavelength and emissivity vectors.
    ///
    /// Wavelengths must be sorted in ascending order.
    #[allow(dead_code)]
    pub fn new(wavelengths: Vec<f64>, emissivities: Vec<f64>) -> Self {
        assert_eq!(
            wavelengths.len(),
            emissivities.len(),
            "wavelengths and emissivities must have the same length"
        );
        Self {
            wavelengths,
            emissivities,
        }
    }
    /// Linearly interpolate emissivity at wavelength `lambda` \[m\].
    ///
    /// Clamps to the first/last value outside the defined range.
    #[allow(dead_code)]
    pub fn emissivity_at(&self, lambda: f64) -> f64 {
        let n = self.wavelengths.len();
        if n == 0 {
            return 0.0;
        }
        if n == 1 {
            return self.emissivities[0];
        }
        if lambda <= self.wavelengths[0] {
            return self.emissivities[0];
        }
        if lambda >= self.wavelengths[n - 1] {
            return self.emissivities[n - 1];
        }
        for i in 0..n - 1 {
            let w0 = self.wavelengths[i];
            let w1 = self.wavelengths[i + 1];
            if lambda >= w0 && lambda <= w1 {
                let t = (lambda - w0) / (w1 - w0);
                return self.emissivities[i] * (1.0 - t) + self.emissivities[i + 1] * t;
            }
        }
        self.emissivities[n - 1]
    }
    /// Compute the total (hemispherical) emissivity at temperature `temp_k`
    /// by integrating `ε(λ) * E_b(λ, T)` over all wavelengths.
    ///
    /// Uses simple trapezoidal rule with `n_steps` intervals over the range
    /// \[100 nm, 100 µm\].
    #[allow(dead_code)]
    pub fn effective_total_emissivity(&self, temp_k: f64, n_steps: usize) -> f64 {
        let lambda_a = 100e-9_f64;
        let lambda_b = 100e-6_f64;
        let n = n_steps.max(2);
        let dl = (lambda_b - lambda_a) / n as f64;
        let e_b_total = blackbody_emissive_power(temp_k);
        if e_b_total < f64::EPSILON {
            return 0.0;
        }
        let mut sum = 0.0;
        for i in 0..=n {
            let lam = lambda_a + i as f64 * dl;
            let w = if i == 0 || i == n { 0.5 } else { 1.0 };
            let eps = self.emissivity_at(lam);
            let e_b_lam = std::f64::consts::PI * blackbody_spectral_intensity(lam, temp_k);
            sum += w * eps * e_b_lam;
        }
        (sum * dl / e_b_total).clamp(0.0, 1.0)
    }
}
/// Irradiation dosimetry model for neutron/gamma radiation exposure.
///
/// Computes absorbed dose, fluence, kerma, and dose-equivalent quantities
/// for radiation shielding and materials damage analysis.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct IrradiationDosimetry {
    /// Particle flux \[particles/(cm²·s)\].
    pub flux: f64,
    /// Exposure time \[s\].
    pub time_s: f64,
    /// Energy-transfer cross-section σ_kin \[cm²\] (for kerma calculation).
    pub kerma_cross_section: f64,
}
impl IrradiationDosimetry {
    /// Create a new dosimetry model.
    #[allow(dead_code)]
    pub fn new(flux: f64, time_s: f64, kerma_cross_section: f64) -> Self {
        Self {
            flux,
            time_s,
            kerma_cross_section,
        }
    }
    /// Total fluence \[particles/cm²\].
    #[allow(dead_code)]
    pub fn fluence(&self) -> f64 {
        self.flux * self.time_s
    }
    /// Absorbed dose \[Gy = J/kg\] for a target material of density `rho_kg_m3`.
    ///
    /// D \[Gy\] = Φ · σ_kin · E_avg / ρ
    ///
    /// Uses a fixed average energy transfer of 1 MeV = 1.602e-13 J.
    #[allow(dead_code)]
    pub fn absorbed_dose_gray(&self) -> f64 {
        let e_transfer = 1.602e-13;
        self.fluence() * self.kerma_cross_section * 1.0e4 * e_transfer
    }
    /// Kerma rate \[Gy/s\] for a material of density `rho` \[kg/m³\].
    ///
    /// `K = Φ · σ_kin · E_n / ρ_target`
    #[allow(dead_code)]
    pub fn kerma_rate(&self, rho: f64) -> f64 {
        let e_transfer = 1.602e-13;
        self.flux * self.kerma_cross_section * 1.0e4 * e_transfer / rho
    }
    /// Dose equivalent \[Sv\] = absorbed dose \[Gy\] × quality factor `qf`.
    ///
    /// Quality factors: photons/electrons = 1, neutrons = 5–20, alpha = 20.
    #[allow(dead_code)]
    pub fn dose_equivalent_sievert(&self, quality_factor: f64) -> f64 {
        self.absorbed_dose_gray() * quality_factor
    }
    /// Equivalent rem dose (1 rem = 0.01 Sv).
    #[allow(dead_code)]
    pub fn rem_dose(&self, quality_factor: f64) -> f64 {
        self.dose_equivalent_sievert(quality_factor) / 0.01
    }
}
/// Numerical integrator for Stefan-Boltzmann-related spectral calculations.
///
/// Provides accurate band-limited and spectrally weighted integration of
/// the Planck blackbody function using Gaussian quadrature (Simpson's rule).
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct StefanBoltzmannIntegrator {
    /// Surface temperature \[K\].
    pub temp_k: f64,
}
impl StefanBoltzmannIntegrator {
    /// Create a new integrator for a blackbody at `temp_k`.
    #[allow(dead_code)]
    pub fn new(temp_k: f64) -> Self {
        Self { temp_k }
    }
    /// Total integrated spectral irradiance over \[100 nm, 1 mm\] \[W/m²\].
    ///
    /// For a blackbody this should approach `σ·T⁴` when the range is wide enough.
    #[allow(dead_code)]
    pub fn integrate_full_spectrum(&self, n_steps: usize) -> f64 {
        let lambda_a = 100e-9_f64;
        let lambda_b = 1e-3_f64;
        self.integrate_range(lambda_a, lambda_b, n_steps)
    }
    /// Integrate spectral irradiance over \[lambda_a, lambda_b\] \[W/m²\].
    ///
    /// Uses Simpson's rule for accuracy.
    #[allow(dead_code)]
    pub fn integrate_range(&self, lambda_a: f64, lambda_b: f64, n_steps: usize) -> f64 {
        let n = if n_steps.is_multiple_of(2) {
            n_steps
        } else {
            n_steps + 1
        };
        let n = n.max(2);
        let dl = (lambda_b - lambda_a) / n as f64;
        let mut sum = 0.0;
        for i in 0..=n {
            let lam = lambda_a + i as f64 * dl;
            let intensity = std::f64::consts::PI * blackbody_spectral_intensity(lam, self.temp_k);
            let w = if i == 0 || i == n {
                1.0
            } else if i % 2 == 1 {
                4.0
            } else {
                2.0
            };
            sum += w * intensity;
        }
        sum * dl / 3.0
    }
    /// Fraction of total blackbody power emitted in wavelength band \[lambda_a, lambda_b\].
    ///
    /// Returns a value in \[0, 1\].
    #[allow(dead_code)]
    pub fn band_fraction_in_range(&self, lambda_a: f64, lambda_b: f64, n_steps: usize) -> f64 {
        let total = blackbody_emissive_power(self.temp_k);
        if total < f64::EPSILON {
            return 0.0;
        }
        let band = self.integrate_range(lambda_a, lambda_b, n_steps);
        (band / total).clamp(0.0, 1.0)
    }
    /// Effective emissivity of a gray coating integrated over the Planck spectrum.
    ///
    /// `ε_eff = ∫ ε(λ) E_b(λ,T) dλ / E_b(T)`
    ///
    /// For a uniform emissivity this equals the emissivity itself.
    #[allow(dead_code)]
    pub fn weighted_emissivity(&self, emissivity_fn: impl Fn(f64) -> f64, n_steps: usize) -> f64 {
        let total = blackbody_emissive_power(self.temp_k);
        if total < f64::EPSILON {
            return 0.0;
        }
        let lambda_a = 100e-9_f64;
        let lambda_b = 1e-3_f64;
        let n = n_steps.max(2);
        let dl = (lambda_b - lambda_a) / n as f64;
        let mut sum = 0.0;
        for i in 0..=n {
            let lam = lambda_a + i as f64 * dl;
            let w = if i == 0 || i == n { 0.5 } else { 1.0 };
            let eps = emissivity_fn(lam);
            let e_b_lam = std::f64::consts::PI * blackbody_spectral_intensity(lam, self.temp_k);
            sum += w * eps * e_b_lam;
        }
        (sum * dl / total).clamp(0.0, 1.0)
    }
}
/// Radiative surface properties satisfying the energy balance
/// `emissivity + reflectivity + transmissivity ≤ 1`.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct RadiativeProperties {
    /// Emissivity ε (0–1)
    pub emissivity: f64,
    /// Absorptivity α (0–1)
    pub absorptivity: f64,
    /// Transmissivity τ (0–1)
    pub transmissivity: f64,
    /// Reflectivity ρ = 1 − α − τ (derived, ≥ 0)
    pub reflectivity: f64,
}
impl RadiativeProperties {
    /// Create from emissivity `eps`, transmissivity `tau`, with
    /// absorptivity set via Kirchhoff's law (`absorptivity = emissivity`)
    /// and reflectivity derived from the energy balance.
    ///
    /// # Panics
    /// Panics in debug mode if `eps + tau > 1`.
    #[allow(dead_code)]
    pub fn new(emissivity: f64, transmissivity: f64) -> Self {
        let absorptivity = emissivity;
        let reflectivity = (1.0 - absorptivity - transmissivity).max(0.0);
        Self {
            emissivity,
            absorptivity,
            transmissivity,
            reflectivity,
        }
    }
    /// Returns `true` if the energy balance is satisfied within tolerance.
    #[allow(dead_code)]
    pub fn is_consistent(&self) -> bool {
        let sum = self.absorptivity + self.transmissivity + self.reflectivity;
        (sum - 1.0).abs() < 1.0e-9
    }
    /// Effective emissive power \[W/m²\] for a surface at temperature `temp_k`.
    #[allow(dead_code)]
    pub fn emissive_power(&self, temp_k: f64) -> f64 {
        self.emissivity * SIGMA * temp_k.powi(4)
    }
}
/// Participating media (gas radiation) model.
///
/// Models radiative transfer through an absorbing/emitting gas via
/// the mean-beam-length approximation.
///
/// Reference: Modest, "Radiative Heat Transfer", 3rd ed.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ParticipatingMedia {
    /// Absorption coefficient κ \[1/m\].
    pub kappa: f64,
    /// Scattering coefficient σ_s \[1/m\].
    pub sigma_s: f64,
    /// Temperature of the medium \[K\].
    pub temperature: f64,
}
impl ParticipatingMedia {
    /// Create a new participating media model.
    pub fn new(kappa: f64, sigma_s: f64, temperature: f64) -> Self {
        Self {
            kappa,
            sigma_s,
            temperature,
        }
    }
    /// Extinction coefficient β = κ + σ_s \[1/m\].
    pub fn extinction(&self) -> f64 {
        self.kappa + self.sigma_s
    }
    /// Single-scattering albedo ω = σ_s / β.
    pub fn albedo(&self) -> f64 {
        let b = self.extinction();
        if b < f64::EPSILON {
            return 0.0;
        }
        self.sigma_s / b
    }
    /// Optical depth over a path length L: τ = β · L.
    pub fn optical_depth(&self, path_length: f64) -> f64 {
        self.extinction() * path_length
    }
    /// Transmittance through a slab of thickness L: T = exp(-τ).
    pub fn transmittance(&self, path_length: f64) -> f64 {
        (-self.optical_depth(path_length)).exp()
    }
    /// Emitted radiation from a gas column of thickness L \[W/m²\].
    ///
    /// q_emit = κ · σ · T⁴ · L  (optically thin approximation)
    pub fn emission_optically_thin(&self, path_length: f64) -> f64 {
        self.kappa * SIGMA * self.temperature.powi(4) * path_length
    }
    /// Effective emissivity of a gas layer of thickness L.
    ///
    /// ε_eff = 1 - exp(-κ·L)  (gray gas approximation)
    pub fn effective_emissivity(&self, path_length: f64) -> f64 {
        1.0 - (-self.kappa * path_length).exp()
    }
    /// Mean-beam length for a sphere of radius R: L_m = 0.65 * 2R.
    pub fn mean_beam_length_sphere(radius: f64) -> f64 {
        0.65 * 2.0 * radius
    }
    /// Mean-beam length for a cylinder (infinite, radius R): L_m = 0.95 * 2R.
    pub fn mean_beam_length_cylinder(radius: f64) -> f64 {
        0.95 * 2.0 * radius
    }
}
/// Nuclear radiation damage model.
///
/// Computes the displacement-per-atom (dpa) dose from fluence and cross-section.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct RadiationDamage {
    /// Neutron fluence \[n/cm²\]
    pub fluence: f64,
    /// Displacement cross-section σ_d \[cm²\]
    pub displacement_cross_section: f64,
    /// Atomic number density N \[atoms/cm³\]
    pub atomic_density: f64,
}
impl RadiationDamage {
    /// Create a new `RadiationDamage` model.
    #[allow(dead_code)]
    pub fn new(fluence: f64, sigma_d: f64, n_atoms: f64) -> Self {
        Self {
            fluence,
            displacement_cross_section: sigma_d,
            atomic_density: n_atoms,
        }
    }
    /// Displacements per atom (dpa).
    ///
    /// `dpa = fluence · σ_d · N / N = fluence · σ_d`
    /// (simplified NRT model, one dpa per unit fluence × cross-section)
    #[allow(dead_code)]
    pub fn dpa(&self) -> f64 {
        self.fluence * self.displacement_cross_section
    }
    /// Fraction of atoms displaced (saturates at ~1).
    #[allow(dead_code)]
    pub fn displaced_fraction(&self) -> f64 {
        self.dpa().min(1.0)
    }
    /// Effective swelling model \[dimensionless volume increase / volume\].
    ///
    /// Simple empirical relation: `ΔV/V = A · dpa^n`, with `A = 1e-3`, `n = 2`.
    #[allow(dead_code)]
    pub fn swelling_fraction(&self) -> f64 {
        let dpa = self.dpa();
        1.0e-3 * dpa * dpa
    }
}
/// Gray body model with temperature-dependent emissivity.
///
/// Represents a real surface whose emissivity varies with temperature using
/// a linear model: `ε(T) = ε₀ + dε/dT · (T - T_ref)`.
#[derive(Debug, Clone, Copy)]
#[allow(dead_code)]
pub struct GrayBodyModel {
    /// Emissivity at reference temperature
    pub eps_ref: f64,
    /// Rate of change of emissivity with temperature \[1/K\]
    pub deps_dt: f64,
    /// Reference temperature \[K\]
    pub temp_ref: f64,
}
impl GrayBodyModel {
    /// Create a new gray body model.
    #[allow(dead_code)]
    pub fn new(eps_ref: f64, deps_dt: f64, temp_ref: f64) -> Self {
        Self {
            eps_ref,
            deps_dt,
            temp_ref,
        }
    }
    /// Emissivity at temperature `temp_k`.
    #[allow(dead_code)]
    pub fn emissivity(&self, temp_k: f64) -> f64 {
        (self.eps_ref + self.deps_dt * (temp_k - self.temp_ref)).clamp(0.0, 1.0)
    }
    /// Emissive power \[W/m²\] at temperature `temp_k`.
    #[allow(dead_code)]
    pub fn emissive_power(&self, temp_k: f64) -> f64 {
        self.emissivity(temp_k) * SIGMA * temp_k.powi(4)
    }
    /// Effective temperature that a blackbody would need to emit the same
    /// total power as this gray body at `temp_k`.
    ///
    /// `T_eff = T * ε^0.25`
    #[allow(dead_code)]
    pub fn effective_blackbody_temperature(&self, temp_k: f64) -> f64 {
        let eps = self.emissivity(temp_k);
        temp_k * eps.powf(0.25)
    }
}
/// Simple 1-D Monte Carlo radiation transport through an absorbing/scattering slab.
///
/// Uses path-length sampling to simulate photon transport through a homogeneous
/// participating medium of thickness `L` with extinction coefficient `β = κ + σ_s`.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MCRadiationTransport {
    /// Number of photon packets.
    pub n_photons: usize,
    /// LCG seed.
    pub(super) seed: u64,
}
impl MCRadiationTransport {
    /// Create a new MC transport calculator.
    #[allow(dead_code)]
    pub fn new(n_photons: usize, seed: u64) -> Self {
        Self { n_photons, seed }
    }
    /// Advance the LCG.
    fn rand_next(&mut self) -> f64 {
        self.seed = self
            .seed
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.seed >> 33) as f64 / (u32::MAX as f64)
    }
    /// Estimate slab transmittance for a purely absorbing medium.
    ///
    /// # Arguments
    /// * `kappa` — absorption coefficient \[1/m\]
    /// * `thickness` — slab thickness \[m\]
    #[allow(dead_code)]
    pub fn slab_transmittance(&mut self, kappa: f64, thickness: f64) -> f64 {
        self.slab_transmittance_full(kappa, thickness, 0.0)
    }
    /// Estimate slab transmittance for an absorbing + scattering medium.
    ///
    /// # Arguments
    /// * `kappa`     — absorption coefficient \[1/m\]
    /// * `thickness` — slab thickness \[m\]
    /// * `albedo`    — single-scattering albedo ω = σ_s / (κ + σ_s)
    #[allow(dead_code)]
    pub fn slab_transmittance_full(&mut self, kappa: f64, thickness: f64, albedo: f64) -> f64 {
        if self.n_photons == 0 {
            return 0.0;
        }
        let beta = kappa / (1.0 - albedo).max(1e-15);
        let mut transmitted = 0usize;
        for _ in 0..self.n_photons {
            let mut x = 0.0_f64;
            let mut weight = 1.0_f64;
            let mut alive = true;
            while alive && x < thickness {
                let xi = self.rand_next().max(1e-15);
                let s = -xi.ln() / beta;
                x += s;
                if x >= thickness {
                    transmitted += 1;
                    break;
                }
                if albedo < f64::EPSILON {
                    alive = false;
                } else {
                    weight *= albedo;
                    if weight < 0.001 {
                        let rr = self.rand_next();
                        if rr < 0.1 {
                            weight /= 0.1;
                        } else {
                            alive = false;
                        }
                    }
                    let cos_theta = 2.0 * self.rand_next() - 1.0;
                    if cos_theta < 0.0 {
                        alive = false;
                    }
                }
            }
        }
        transmitted as f64 / self.n_photons as f64
    }
    /// Estimate the mean free path \[m\] in a medium with extinction coefficient `beta`.
    ///
    /// Analytical result: `mfp = 1/beta`. This method provides the MC estimate.
    #[allow(dead_code)]
    pub fn estimate_mean_free_path(&mut self, beta: f64) -> f64 {
        let mut total_path = 0.0;
        let n = self.n_photons.max(1);
        for _ in 0..n {
            let xi = self.rand_next().max(1e-15);
            total_path += -xi.ln() / beta;
        }
        total_path / n as f64
    }
}
