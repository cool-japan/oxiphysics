//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    CO2_PREINDUSTRIAL, CP_SEAWATER, E9X, E9Y, EARTH_RADIUS, LATENT_HEAT_ICE, RHO_SEAWATER,
    SOLAR_CONSTANT, STEFAN_BOLTZMANN, W9,
};

/// Charney climate sensitivity and transient climate response.
///
/// Climate sensitivity *S* relates equilibrium global-mean temperature change
/// to CO₂ doubling: S = ΔF₂× / λ_net, where λ_net combines Planck response
/// and all feedback parameters.
#[derive(Debug, Clone)]
pub struct ClimateSensitivity {
    /// Planck feedback parameter (W m⁻² K⁻¹) — intrinsic cooling, negative.
    pub planck_feedback: f64,
    /// Water vapour feedback (W m⁻² K⁻¹) — positive.
    pub water_vapour_feedback: f64,
    /// Lapse-rate feedback (W m⁻² K⁻¹) — negative in tropics.
    pub lapse_rate_feedback: f64,
    /// Cloud feedback (W m⁻² K⁻¹) — uncertain, positive or negative.
    pub cloud_feedback: f64,
    /// Ice-albedo feedback (W m⁻² K⁻¹) — positive.
    pub ice_albedo_feedback: f64,
    /// Ocean heat uptake efficiency κ (W m⁻² K⁻¹).
    pub ocean_heat_efficiency: f64,
}
impl ClimateSensitivity {
    /// Create with best-estimate IPCC AR6 feedback values.
    pub fn new() -> Self {
        Self {
            planck_feedback: -3.22,
            water_vapour_feedback: 1.77,
            lapse_rate_feedback: -0.50,
            cloud_feedback: 0.42,
            ice_albedo_feedback: 0.35,
            ocean_heat_efficiency: 0.60,
        }
    }
    /// Net feedback parameter λ (W m⁻² K⁻¹).
    pub fn net_feedback(&self) -> f64 {
        self.planck_feedback
            + self.water_vapour_feedback
            + self.lapse_rate_feedback
            + self.cloud_feedback
            + self.ice_albedo_feedback
    }
    /// Equilibrium climate sensitivity ECS (K per CO₂ doubling).
    pub fn ecs(&self) -> f64 {
        -Co2FeedbackModel::co2_doubling_forcing() / self.net_feedback()
    }
    /// Transient climate response TCR (K) — at time of doubling.
    ///
    /// Accounts for ocean heat uptake: TCR = ECS · (1 − κ / |λ|).
    pub fn tcr(&self) -> f64 {
        let lam = self.net_feedback().abs();
        self.ecs() * (1.0 - self.ocean_heat_efficiency / lam)
    }
    /// Temperature change for a given forcing (W m⁻²) at equilibrium.
    pub fn delta_t_eq(&self, forcing: f64) -> f64 {
        -forcing / self.net_feedback()
    }
    /// Feedback factor (ratio of response with feedbacks to Planck-only).
    pub fn feedback_factor(&self) -> f64 {
        self.planck_feedback / self.net_feedback()
    }
}
/// Monsoon index based on land-sea thermal contrast.
///
/// A positive index indicates monsoon onset (land warmer than sea).
#[derive(Debug, Clone)]
pub struct MonsoonIndex {
    /// Land temperature (°C).
    pub land_temp: f64,
    /// Sea-surface temperature (°C).
    pub sea_temp: f64,
    /// Critical contrast threshold for onset (K).
    pub onset_threshold: f64,
}
impl MonsoonIndex {
    /// Create with default values.
    pub fn new() -> Self {
        Self {
            land_temp: 30.0,
            sea_temp: 27.0,
            onset_threshold: 2.0,
        }
    }
    /// Land-sea temperature contrast (K).
    pub fn contrast(&self) -> f64 {
        self.land_temp - self.sea_temp
    }
    /// Return true if monsoon conditions are active.
    pub fn is_monsoon_active(&self) -> bool {
        self.contrast() > self.onset_threshold
    }
    /// Monsoon precipitation proxy (mm day⁻¹).
    pub fn precipitation_proxy(&self) -> f64 {
        (self.contrast() - self.onset_threshold).max(0.0) * 3.5
    }
}
/// 2-D LBM grid for simplified global atmospheric circulation.
///
/// Each cell represents a latitude/longitude tile and carries:
/// - `f[i]` — distribution functions for 9 velocity directions
/// - `temp` — surface air temperature (°C)
/// - `precip` — precipitation (mm day⁻¹)
/// - `albedo` — local albedo
/// - `co2_forcing` — local CO₂ radiative forcing (W m⁻²)
#[derive(Debug, Clone)]
pub struct ClimateLbmGrid {
    /// Number of longitude grid points.
    pub nx: usize,
    /// Number of latitude grid points.
    pub ny: usize,
    /// Distribution functions \[ny\]\[nx\]\[9\].
    pub f: Vec<Vec<[f64; 9]>>,
    /// Temperature field (°C).
    pub temp: Vec<Vec<f64>>,
    /// Precipitation field (mm day⁻¹).
    pub precip: Vec<Vec<f64>>,
    /// Albedo field.
    pub albedo: Vec<Vec<f64>>,
    /// Insolation field (W m⁻²).
    pub insolation: Vec<Vec<f64>>,
    /// BGK relaxation parameter ω.
    pub omega: f64,
}
impl ClimateLbmGrid {
    /// Create a grid with `nx` × `ny` cells and uniform initial temperature.
    pub fn new(nx: usize, ny: usize, t0: f64) -> Self {
        let f_eq = Self::f_equilibrium_uniform(t0);
        let f = vec![vec![f_eq; nx]; ny];
        let temp = vec![vec![t0; nx]; ny];
        let precip = vec![vec![0.0; nx]; ny];
        let albedo = vec![vec![0.30; nx]; ny];
        let insolation = Self::compute_insolation(nx, ny);
        Self {
            nx,
            ny,
            f,
            temp,
            precip,
            albedo,
            insolation,
            omega: 1.8,
        }
    }
    fn f_equilibrium_uniform(t: f64) -> [f64; 9] {
        let mut feq = [0.0_f64; 9];
        for (q, feq_q) in feq.iter_mut().enumerate() {
            *feq_q = W9[q] * t;
        }
        feq
    }
    fn compute_insolation(nx: usize, ny: usize) -> Vec<Vec<f64>> {
        let mut ins = vec![vec![0.0_f64; nx]; ny];
        for j in 0..ny {
            let lat = PI * (j as f64 / (ny as f64 - 1.0) - 0.5);
            let q = SOLAR_CONSTANT / 4.0 * lat.cos().max(0.0);
            for i in 0..nx {
                ins[j][i] = q;
            }
        }
        ins
    }
    /// Compute D2Q9 equilibrium distribution for a cell with temperature *t*
    /// and macroscopic velocity (*ux*, *uy*).
    pub fn f_eq(t: f64, ux: f64, uy: f64) -> [f64; 9] {
        let cs2 = 1.0 / 3.0;
        let u2 = ux * ux + uy * uy;
        let mut feq = [0.0_f64; 9];
        for q in 0..9 {
            let eu = E9X[q] as f64 * ux + E9Y[q] as f64 * uy;
            feq[q] = W9[q] * t * (1.0 + eu / cs2 + eu * eu / (2.0 * cs2 * cs2) - u2 / (2.0 * cs2));
        }
        feq
    }
    /// BGK collision at a single cell.
    fn collide_cell(f: &[f64; 9], t: f64, ux: f64, uy: f64, omega: f64) -> [f64; 9] {
        let feq = Self::f_eq(t, ux, uy);
        let mut f_new = [0.0_f64; 9];
        for q in 0..9 {
            f_new[q] = f[q] - omega * (f[q] - feq[q]);
        }
        f_new
    }
    /// Perform one collision step (BGK) across the whole grid.
    pub fn collide(&mut self) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let t = self.temp[j][i];
                let (ux, uy) = self.compute_velocity(j, i);
                self.f[j][i] = Self::collide_cell(&self.f[j][i], t, ux, uy, self.omega);
            }
        }
    }
    /// Perform one streaming step (periodic boundaries).
    pub fn stream(&mut self) {
        let nx = self.nx;
        let ny = self.ny;
        let mut f_new = vec![vec![[0.0_f64; 9]; nx]; ny];
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    let src_i = ((i as isize - E9X[q] as isize).rem_euclid(nx as isize)) as usize;
                    let src_j = ((j as isize - E9Y[q] as isize).rem_euclid(ny as isize)) as usize;
                    f_new[j][i][q] = self.f[src_j][src_i][q];
                }
            }
        }
        self.f = f_new;
    }
    /// Update macroscopic temperature field from zeroth moment.
    pub fn update_temperature(&mut self) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                self.temp[j][i] = self.f[j][i].iter().sum();
            }
        }
    }
    /// Compute macroscopic velocity (thermal wind proxy) at a cell.
    pub fn compute_velocity(&self, j: usize, i: usize) -> (f64, f64) {
        let t = self.temp[j][i];
        if t <= 0.0 {
            return (0.0, 0.0);
        }
        let mut ux = 0.0_f64;
        let mut uy = 0.0_f64;
        for q in 0..9 {
            ux += E9X[q] as f64 * self.f[j][i][q];
            uy += E9Y[q] as f64 * self.f[j][i][q];
        }
        (ux / t, uy / t)
    }
    /// Apply solar forcing and OLR relaxation to temperature.
    ///
    /// `dt` is in units of LBM time steps (each ~1 day).
    pub fn apply_forcing(&mut self, dt: f64, extra_forcing: f64) {
        let ebm = EnergyBalanceModel::new();
        for j in 0..self.ny {
            for i in 0..self.nx {
                let absorbed = self.insolation[j][i] * (1.0 - self.albedo[j][i]);
                let olr = ebm.olr_a + ebm.olr_b * self.temp[j][i];
                let flux = absorbed - olr + extra_forcing;
                self.temp[j][i] += dt * flux / ebm.heat_capacity;
            }
        }
    }
    /// Apply simple precipitation parameterisation:
    /// precipitation ∝ max(T - T_crit, 0).
    pub fn update_precipitation(&mut self, t_crit: f64) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                self.precip[j][i] = (self.temp[j][i] - t_crit).max(0.0) * 0.5;
            }
        }
    }
    /// Global mean temperature (°C).
    pub fn global_mean_temperature(&self) -> f64 {
        let sum: f64 = self.temp.iter().flat_map(|row| row.iter()).sum();
        sum / (self.nx * self.ny) as f64
    }
    /// Global mean albedo.
    pub fn global_mean_albedo(&self) -> f64 {
        let sum: f64 = self.albedo.iter().flat_map(|row| row.iter()).sum();
        sum / (self.nx * self.ny) as f64
    }
    /// Advance by one combined LBM + forcing step.
    pub fn step(&mut self, dt_force: f64, extra_forcing: f64) {
        self.collide();
        self.stream();
        self.update_temperature();
        self.apply_forcing(dt_force, extra_forcing);
        self.update_precipitation(20.0);
    }
}
/// Multi-gas greenhouse effect model.
///
/// Tracks CO₂ (ppm), CH₄ (ppb), and N₂O (ppb) concentrations and computes
/// effective radiative forcing, effective emissivity, and greenhouse
/// temperature enhancement.
#[derive(Debug, Clone)]
pub struct GreenhouseGasEffect {
    /// CO₂ concentration (ppm).
    pub co2_ppm: f64,
    /// CH₄ concentration (ppb).
    pub ch4_ppb: f64,
    /// N₂O concentration (ppb).
    pub n2o_ppb: f64,
    /// Pre-industrial CO₂ (ppm).
    pub co2_ref: f64,
    /// Pre-industrial CH₄ (ppb).
    pub ch4_ref: f64,
    /// Pre-industrial N₂O (ppb).
    pub n2o_ref: f64,
}
impl GreenhouseGasEffect {
    /// Create with pre-industrial concentrations.
    pub fn new() -> Self {
        Self {
            co2_ppm: CO2_PREINDUSTRIAL,
            ch4_ppb: 722.0,
            n2o_ppb: 270.0,
            co2_ref: CO2_PREINDUSTRIAL,
            ch4_ref: 722.0,
            n2o_ref: 270.0,
        }
    }
    /// CO₂ radiative forcing relative to reference (W m⁻²).
    pub fn co2_forcing(&self) -> f64 {
        5.35 * (self.co2_ppm / self.co2_ref).ln()
    }
    /// N₂O simplified forcing (W m⁻²).
    pub fn n2o_forcing(&self) -> f64 {
        0.12 * (self.n2o_ppb.sqrt() - self.n2o_ref.sqrt())
    }
    /// CH₄ simplified forcing (W m⁻²).
    pub fn ch4_forcing(&self) -> f64 {
        RadiativeForcing::from_ch4(self.ch4_ppb, self.ch4_ref, self.n2o_ppb)
    }
    /// Total GHG radiative forcing (W m⁻²).
    pub fn total_forcing(&self) -> f64 {
        self.co2_forcing() + self.ch4_forcing() + self.n2o_forcing()
    }
    /// Effective atmospheric emissivity.
    ///
    /// Approximated as ε = 1 − exp(−τ) where τ is derived from forcing.
    pub fn effective_emissivity(&self) -> f64 {
        let tau_base = 0.78;
        let delta_tau = self.total_forcing() / (4.0 * STEFAN_BOLTZMANN * 255_f64.powi(3));
        (tau_base + delta_tau).clamp(0.0, 1.0)
    }
    /// Greenhouse temperature enhancement above no-atmosphere case (K).
    pub fn greenhouse_enhancement(&self) -> f64 {
        let t_eff = (SOLAR_CONSTANT * (1.0 - 0.3) / (4.0 * STEFAN_BOLTZMANN)).powf(0.25);
        let eps = self.effective_emissivity();
        let t_surface = t_eff / (1.0 - eps / 2.0).powf(0.25);
        t_surface - t_eff
    }
}
/// Ice-albedo positive feedback parameterisation.
///
/// As temperature drops below the freezing threshold, sea-ice / land-ice
/// fraction increases, raising the planetary albedo and amplifying cooling.
#[derive(Debug, Clone)]
pub struct IceAlbedoFeedback {
    /// Ice-free albedo (ocean/land).
    pub albedo_ice_free: f64,
    /// Fully ice-covered albedo.
    pub albedo_ice_covered: f64,
    /// Temperature at which ice begins forming (°C).
    pub freeze_threshold: f64,
    /// Temperature width of the ice-transition zone (°C).
    pub transition_width: f64,
    /// Current surface temperature (°C).
    pub temperature: f64,
}
impl IceAlbedoFeedback {
    /// Create with default Arctic parameters.
    pub fn new() -> Self {
        Self {
            albedo_ice_free: 0.06,
            albedo_ice_covered: 0.80,
            freeze_threshold: -2.0,
            transition_width: 3.0,
            temperature: 10.0,
        }
    }
    /// Ice fraction as a smooth sigmoid of temperature.
    pub fn ice_fraction(&self) -> f64 {
        let x = -(self.temperature - self.freeze_threshold) / self.transition_width;
        1.0 / (1.0 + (-x).exp())
    }
    /// Effective surface albedo including ice fraction.
    pub fn effective_albedo(&self) -> f64 {
        let fi = self.ice_fraction();
        (1.0 - fi) * self.albedo_ice_free + fi * self.albedo_ice_covered
    }
    /// Change in albedo per degree warming (dα/dT).
    pub fn dalpha_dt(&self) -> f64 {
        let fi = self.ice_fraction();
        let dfi_dt = -fi * (1.0 - fi) / self.transition_width;
        dfi_dt * (self.albedo_ice_covered - self.albedo_ice_free)
    }
    /// Ice-albedo feedback parameter \[W m⁻² K⁻¹\].
    pub fn feedback_parameter(&self) -> f64 {
        -(SOLAR_CONSTANT / 4.0) * self.dalpha_dt()
    }
}
/// Simplified Hadley cell parameterisation.
///
/// Relates the Hadley cell width and strength to the tropical SST gradient
/// following Held (2000) scaling.
#[derive(Debug, Clone)]
pub struct HadleyCell {
    /// Equatorial SST (°C).
    pub sst_equator: f64,
    /// Subtropical SST (°C).
    pub sst_subtropics: f64,
    /// Tropopause height (m).
    pub tropopause_height: f64,
    /// Rotation rate (rad s⁻¹).
    pub omega_earth: f64,
}
impl HadleyCell {
    /// Create with typical tropical values.
    pub fn new() -> Self {
        Self {
            sst_equator: 28.0,
            sst_subtropics: 20.0,
            tropopause_height: 15_000.0,
            omega_earth: 7.292e-5,
        }
    }
    /// Meridional temperature gradient Δθ / Δy (K m⁻¹).
    pub fn temperature_gradient(&self) -> f64 {
        let delta_lat = 30.0_f64.to_radians() * EARTH_RADIUS;
        (self.sst_equator - self.sst_subtropics) / delta_lat
    }
    /// Held scaling for Hadley cell extent (radians).
    ///
    /// φ_H ≈ (5 Δh / (18 Ω² a²)) ^ (1/2)
    pub fn cell_width_radians(&self) -> f64 {
        let delta_h =
            9.81 * (self.sst_equator - self.sst_subtropics) / 300.0 * self.tropopause_height;
        let factor = 5.0 * delta_h / (18.0 * self.omega_earth.powi(2) * EARTH_RADIUS.powi(2));
        factor.sqrt().min(PI / 4.0)
    }
    /// Maximum upwelling velocity (m s⁻¹) proxy.
    pub fn upwelling_velocity(&self) -> f64 {
        let grad = self.temperature_gradient();
        2.0e-3 * grad.abs() * EARTH_RADIUS
    }
}
/// Zonal-mean temperature profile on a latitude grid.
///
/// Derived from the LBM grid by averaging over longitude bands.
#[derive(Debug, Clone)]
pub struct ZonalMeanProfile {
    /// Latitude grid (radians, −π/2 to π/2).
    pub latitudes: Vec<f64>,
    /// Zonal mean temperature at each latitude (°C).
    pub temperature: Vec<f64>,
    /// Zonal mean albedo at each latitude.
    pub albedo: Vec<f64>,
    /// Zonal mean precipitation (mm day⁻¹).
    pub precip: Vec<f64>,
}
impl ZonalMeanProfile {
    /// Compute zonal means from a [`ClimateLbmGrid`].
    pub fn from_grid(grid: &ClimateLbmGrid) -> Self {
        let ny = grid.ny;
        let nx = grid.nx as f64;
        let latitudes: Vec<f64> = (0..ny)
            .map(|j| PI * (j as f64 / (ny as f64 - 1.0) - 0.5))
            .collect();
        let temperature: Vec<f64> = (0..ny)
            .map(|j| grid.temp[j].iter().sum::<f64>() / nx)
            .collect();
        let albedo: Vec<f64> = (0..ny)
            .map(|j| grid.albedo[j].iter().sum::<f64>() / nx)
            .collect();
        let precip: Vec<f64> = (0..ny)
            .map(|j| grid.precip[j].iter().sum::<f64>() / nx)
            .collect();
        Self {
            latitudes,
            temperature,
            albedo,
            precip,
        }
    }
    /// Meridional temperature gradient (K rad⁻¹) at index *j*.
    pub fn meridional_gradient(&self, j: usize) -> f64 {
        let n = self.latitudes.len();
        if j == 0 || j >= n - 1 {
            return 0.0;
        }
        (self.temperature[j + 1] - self.temperature[j - 1])
            / (self.latitudes[j + 1] - self.latitudes[j - 1])
    }
}
/// Arctic amplification factor calculation.
///
/// The Arctic warms approximately 2–4× faster than the global mean due to
/// ice-albedo, lapse-rate, and Planck feedbacks.
#[derive(Debug, Clone)]
pub struct ArcticAmplification {
    /// Global mean temperature change (K).
    pub global_delta_t: f64,
    /// Arctic (60–90°N) mean temperature change (K).
    pub arctic_delta_t: f64,
    /// Arctic amplification factor (dimensionless).
    pub amplification_factor: f64,
}
impl ArcticAmplification {
    /// Compute from zonal mean profiles at two times.
    pub fn compute(profile_now: &ZonalMeanProfile, profile_ref: &ZonalMeanProfile) -> Self {
        let n = profile_now.temperature.len();
        let t_now: f64 = profile_now.temperature.iter().sum::<f64>() / n as f64;
        let t_ref: f64 = profile_ref.temperature.iter().sum::<f64>() / n as f64;
        let global_delta_t = t_now - t_ref;
        let arctic_start = 5 * n / 6;
        let na = n - arctic_start;
        let t_arc_now: f64 =
            profile_now.temperature[arctic_start..].iter().sum::<f64>() / na as f64;
        let t_arc_ref: f64 =
            profile_ref.temperature[arctic_start..].iter().sum::<f64>() / na as f64;
        let arctic_delta_t = t_arc_now - t_arc_ref;
        let amplification_factor = if global_delta_t.abs() > 1e-10 {
            arctic_delta_t / global_delta_t
        } else {
            1.0
        };
        Self {
            global_delta_t,
            arctic_delta_t,
            amplification_factor,
        }
    }
    /// Return true if the amplification is "strong" (> 2).
    pub fn is_strong(&self) -> bool {
        self.amplification_factor > 2.0
    }
}
/// Milankovitch orbital forcing.
///
/// Computes eccentricity, obliquity, and precession contributions to
/// insolation at the summer solstice.
#[derive(Debug, Clone)]
pub struct MilankovitchForcing {
    /// Orbital eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Obliquity / axial tilt (radians).
    pub obliquity: f64,
    /// Precession angle (radians).
    pub precession: f64,
}
impl MilankovitchForcing {
    /// Create with present-day orbital parameters.
    pub fn present_day() -> Self {
        Self {
            eccentricity: 0.0167,
            obliquity: 23.44_f64.to_radians(),
            precession: 102.9_f64.to_radians(),
        }
    }
    /// Top-of-atmosphere insolation at latitude φ and day-of-year d (W m⁻²).
    pub fn insolation_toa(&self, phi: f64, day: f64) -> f64 {
        let mean_anomaly = 2.0 * PI * day / 365.25;
        let true_anomaly = mean_anomaly + 2.0 * self.eccentricity * mean_anomaly.sin();
        let r_factor = (1.0 + self.eccentricity * true_anomaly.cos()).powi(2);
        let declination = self.obliquity * (mean_anomaly + self.precession).sin();
        let cos_zen = phi.sin() * declination.sin() + phi.cos() * declination.cos();
        SOLAR_CONSTANT * r_factor * cos_zen.max(0.0)
    }
    /// Summer solstice insolation at 65°N (key Milankovitch diagnostic, W m⁻²).
    pub fn insolation_65n_summer(&self) -> f64 {
        self.insolation_toa(65.0_f64.to_radians(), 172.0)
    }
}
/// Simple cloud feedback parameterisation.
///
/// Low cloud fraction decreases with temperature (positive feedback);
/// high cloud fraction increases (negative feedback).
#[derive(Debug, Clone)]
pub struct CloudFeedback {
    /// Low cloud fraction at reference temperature.
    pub low_cloud_ref: f64,
    /// High cloud fraction at reference temperature.
    pub high_cloud_ref: f64,
    /// Low-cloud sensitivity to temperature (K⁻¹).
    pub low_cloud_sensitivity: f64,
    /// High-cloud sensitivity to temperature (K⁻¹).
    pub high_cloud_sensitivity: f64,
    /// Reference temperature (°C).
    pub t_ref: f64,
    /// Current temperature (°C).
    pub temperature: f64,
}
impl CloudFeedback {
    /// Create with best-estimate parameters.
    pub fn new() -> Self {
        Self {
            low_cloud_ref: 0.35,
            high_cloud_ref: 0.20,
            low_cloud_sensitivity: -0.006,
            high_cloud_sensitivity: 0.003,
            t_ref: 14.0,
            temperature: 14.0,
        }
    }
    /// Current low cloud fraction.
    pub fn low_cloud_fraction(&self) -> f64 {
        let delta_t = self.temperature - self.t_ref;
        (self.low_cloud_ref + self.low_cloud_sensitivity * delta_t).clamp(0.0, 1.0)
    }
    /// Current high cloud fraction.
    pub fn high_cloud_fraction(&self) -> f64 {
        let delta_t = self.temperature - self.t_ref;
        (self.high_cloud_ref + self.high_cloud_sensitivity * delta_t).clamp(0.0, 1.0)
    }
    /// Net cloud radiative effect (W m⁻²).
    pub fn net_cloud_forcing(&self) -> f64 {
        let sw_effect = -50.0 * self.low_cloud_fraction();
        let lw_effect = 25.0 * self.high_cloud_fraction();
        sw_effect + lw_effect
    }
    /// Cloud feedback parameter (W m⁻² K⁻¹).
    pub fn cloud_feedback_parameter(&self) -> f64 {
        let cf0 = -50.0 * self.low_cloud_ref + 25.0 * self.high_cloud_ref;
        let cf1 = -50.0 * (self.low_cloud_ref + self.low_cloud_sensitivity)
            + 25.0 * (self.high_cloud_ref + self.high_cloud_sensitivity);
        cf1 - cf0
    }
}
/// High-level climate simulation driver.
///
/// Orchestrates [`ClimateLbmGrid`], [`Co2FeedbackModel`], [`IceAlbedoFeedback`],
/// [`OceanHeatUptake`], [`SeaLevelRiseModel`], and [`CarbonCycle`] together.
#[derive(Debug)]
pub struct ClimateSimulation {
    /// 2-D LBM atmospheric circulation grid.
    pub grid: ClimateLbmGrid,
    /// CO₂ feedback model.
    pub co2: Co2FeedbackModel,
    /// Ice-albedo feedback.
    pub ice: IceAlbedoFeedback,
    /// Ocean heat uptake.
    pub ocean: OceanHeatUptake,
    /// Sea-level rise model.
    pub sea_level: SeaLevelRiseModel,
    /// Carbon cycle model.
    pub carbon: CarbonCycle,
    /// Climate sensitivity parameters.
    pub sensitivity: ClimateSensitivity,
    /// Greenhouse gas effect.
    pub ghg: GreenhouseGasEffect,
    /// Simulation year.
    pub year: f64,
}
impl ClimateSimulation {
    /// Create a new simulation starting at year `start_year`.
    pub fn new(nx: usize, ny: usize, start_year: f64) -> Self {
        Self {
            grid: ClimateLbmGrid::new(nx, ny, 14.0),
            co2: Co2FeedbackModel::new(),
            ice: IceAlbedoFeedback::new(),
            ocean: OceanHeatUptake::new(),
            sea_level: SeaLevelRiseModel::new(),
            carbon: CarbonCycle::new(),
            sensitivity: ClimateSensitivity::new(),
            ghg: GreenhouseGasEffect::new(),
            year: start_year,
        }
    }
    /// Advance simulation by `dt` years.
    pub fn step(&mut self, dt: f64) {
        self.co2.step(dt);
        self.ghg.co2_ppm = self.co2.co2_ppm;
        let forcing = self.co2.radiative_forcing();
        self.grid.step(dt, forcing);
        let t_global = self.grid.global_mean_temperature();
        self.ice.temperature = t_global;
        let new_albedo = self.ice.effective_albedo();
        for row in &mut self.grid.albedo {
            for a in row.iter_mut() {
                *a = new_albedo;
            }
        }
        self.ocean.atm_temperature = t_global;
        self.ocean.step(dt * 3.156e7);
        let delta_t = t_global - self.sea_level.reference_temperature;
        self.sea_level.step(dt, delta_t);
        self.carbon.emission_rate = self.co2.emission_rate;
        self.carbon.step(dt);
        self.year += dt;
    }
    /// Run from `start_year` to `end_year` with time step `dt`.
    pub fn run(&mut self, end_year: f64, dt: f64) {
        while self.year < end_year {
            self.step(dt);
        }
    }
    /// Summary snapshot of the simulation state.
    pub fn snapshot(&self) -> ClimateSnapshot {
        ClimateSnapshot {
            year: self.year,
            co2_ppm: self.co2.co2_ppm,
            global_mean_temperature: self.grid.global_mean_temperature(),
            sea_level_rise: self.sea_level.sea_level,
            ocean_temperature: self.ocean.temperature,
            radiative_forcing: self.co2.radiative_forcing(),
            ice_fraction: self.ice.ice_fraction(),
            atmospheric_carbon: self.carbon.atmosphere,
        }
    }
}
/// Atmospheric CO₂ feedback model including emissions, ocean uptake, and
/// land-carbon response.
///
/// Tracks the atmospheric CO₂ concentration \[ppm\] and computes the associated
/// logarithmic radiative forcing.
#[derive(Debug, Clone)]
pub struct Co2FeedbackModel {
    /// Current CO₂ concentration (ppm).
    pub co2_ppm: f64,
    /// Annual anthropogenic emission rate (GtC yr⁻¹).
    pub emission_rate: f64,
    /// Ocean uptake fraction (dimensionless, typically 0.25–0.30).
    pub ocean_uptake_fraction: f64,
    /// Land uptake fraction (dimensionless, typically 0.10–0.20).
    pub land_uptake_fraction: f64,
    /// Beta (land-carbon CO₂ fertilisation, ppm⁻¹).
    pub beta_land: f64,
    /// Simulation year.
    pub year: f64,
}
impl Co2FeedbackModel {
    /// Create a model starting at pre-industrial CO₂.
    pub fn new() -> Self {
        Self {
            co2_ppm: CO2_PREINDUSTRIAL,
            emission_rate: 10.0,
            ocean_uptake_fraction: 0.26,
            land_uptake_fraction: 0.30,
            beta_land: 0.0,
            year: 1850.0,
        }
    }
    /// Radiative forcing relative to pre-industrial \[W m⁻²\].
    ///
    /// Uses the IPCC AR6 simplified expression: ΔF = 5.35 ln(C/C₀).
    pub fn radiative_forcing(&self) -> f64 {
        5.35 * (self.co2_ppm / CO2_PREINDUSTRIAL).ln()
    }
    /// Effective CO₂ doubling forcing (W m⁻²) ≈ 3.71.
    pub fn co2_doubling_forcing() -> f64 {
        5.35 * 2_f64.ln()
    }
    /// Advance CO₂ concentration by `dt` years.
    pub fn step(&mut self, dt: f64) {
        let beta = self.beta_land * (self.co2_ppm - CO2_PREINDUSTRIAL);
        let land_frac = (self.land_uptake_fraction + beta).clamp(0.0, 0.50);
        let net_fraction = 1.0 - self.ocean_uptake_fraction - land_frac;
        let dppm = self.emission_rate * net_fraction * 0.4717 * dt;
        self.co2_ppm += dppm;
        self.year += dt;
    }
    /// Return the year at which CO₂ doubles from pre-industrial.
    pub fn doubling_year(&self) -> Option<f64> {
        if self.co2_ppm >= CO2_PREINDUSTRIAL * 2.0 {
            Some(self.year)
        } else {
            None
        }
    }
}
/// Snapshot of climate simulation state at a given time.
#[derive(Debug, Clone)]
pub struct ClimateSnapshot {
    /// Simulation year.
    pub year: f64,
    /// Atmospheric CO₂ (ppm).
    pub co2_ppm: f64,
    /// Global mean surface temperature (°C).
    pub global_mean_temperature: f64,
    /// Cumulative sea-level rise (m).
    pub sea_level_rise: f64,
    /// Ocean mixed-layer temperature (°C).
    pub ocean_temperature: f64,
    /// CO₂ radiative forcing (W m⁻²).
    pub radiative_forcing: f64,
    /// Ice fraction (0–1).
    pub ice_fraction: f64,
    /// Atmospheric carbon stock (GtC).
    pub atmospheric_carbon: f64,
}
/// Permafrost thaw and carbon release model.
///
/// Uses a simple Stefan freezing-depth model to estimate active-layer
/// depth and associated carbon release.
#[derive(Debug, Clone)]
pub struct PermafrostModel {
    /// Mean annual temperature (°C).
    pub mat: f64,
    /// Active layer depth (m).
    pub active_layer_depth: f64,
    /// Soil carbon density (kgC m⁻³).
    pub carbon_density: f64,
    /// Soil thermal diffusivity (m² s⁻¹).
    pub thermal_diffusivity: f64,
    /// Decomposition rate at 0 °C (yr⁻¹).
    pub decomposition_rate: f64,
    /// Q10 temperature coefficient.
    pub q10: f64,
}
impl PermafrostModel {
    /// Create with typical Siberian permafrost parameters.
    pub fn new() -> Self {
        Self {
            mat: -5.0,
            active_layer_depth: 0.5,
            carbon_density: 25.0,
            thermal_diffusivity: 5e-7,
            decomposition_rate: 0.05,
            q10: 2.0,
        }
    }
    /// Stefan formula for active layer depth (m).
    ///
    /// d = √(2 κ DDT / L_f · ρ_ice)
    pub fn stefan_depth(&self) -> f64 {
        let thaw_days = (self.mat + 5.0).max(0.0) * 365.0;
        let ddt = thaw_days * 86400.0;
        let lf_rho = LATENT_HEAT_ICE * 917.0;
        (2.0 * self.thermal_diffusivity * ddt / lf_rho).sqrt()
    }
    /// Carbon release rate (kgC m⁻² yr⁻¹).
    pub fn carbon_release_rate(&self) -> f64 {
        let depth = self.stefan_depth();
        let r = self.decomposition_rate * self.q10.powf(self.mat / 10.0);
        r * self.carbon_density * depth
    }
    /// Update active layer depth for new temperature (°C).
    pub fn update_temperature(&mut self, new_mat: f64) {
        self.mat = new_mat;
        self.active_layer_depth = self.stefan_depth().max(0.01);
    }
}
/// Zero/one-dimensional energy-balance climate model (Budyko-Sellers).
///
/// The model integrates
/// ```text
/// C dT/dt = Q(1 − α) − A − B·T + D·∇²T
/// ```
/// where *T* is surface temperature (°C), *Q* = S₀/4 the mean insolation,
/// *α* the planetary albedo, *A* and *B* the OLR fit coefficients, and
/// *D* the meridional diffusivity.
#[derive(Debug, Clone)]
pub struct EnergyBalanceModel {
    /// Surface temperature (°C).
    pub temperature: f64,
    /// Planetary albedo (dimensionless, 0–1).
    pub albedo: f64,
    /// Thermal inertia / heat capacity (W yr m⁻² K⁻¹).
    pub heat_capacity: f64,
    /// Budyko OLR intercept A (W m⁻²).
    pub olr_a: f64,
    /// Budyko OLR slope B (W m⁻² K⁻¹).
    pub olr_b: f64,
    /// Meridional diffusivity (W m⁻² K⁻¹) — used in 1-D variant.
    pub diffusivity: f64,
    /// Simulation time (years).
    pub time: f64,
}
impl EnergyBalanceModel {
    /// Create a new EBM with typical pre-industrial parameters.
    pub fn new() -> Self {
        Self {
            temperature: 14.0,
            albedo: 0.30,
            heat_capacity: 10.0,
            olr_a: 203.3,
            olr_b: 2.09,
            diffusivity: 0.65,
            time: 0.0,
        }
    }
    /// Compute net absorbed short-wave radiation (W m⁻²).
    pub fn absorbed_shortwave(&self) -> f64 {
        (SOLAR_CONSTANT / 4.0) * (1.0 - self.albedo)
    }
    /// Compute outgoing long-wave radiation using the Budyko linearisation
    /// OLR = A + B·T (W m⁻²).
    pub fn outgoing_longwave(&self) -> f64 {
        self.olr_a + self.olr_b * self.temperature
    }
    /// Net top-of-atmosphere energy imbalance (W m⁻²).
    pub fn net_toa_flux(&self) -> f64 {
        self.absorbed_shortwave() - self.outgoing_longwave()
    }
    /// Integrate one time step of size `dt` years with optional forcing
    /// `extra_forcing` (W m⁻²).
    pub fn step(&mut self, dt: f64, extra_forcing: f64) {
        let flux = self.net_toa_flux() + extra_forcing;
        self.temperature += dt * flux / self.heat_capacity;
        self.time += dt;
    }
    /// Run until equilibrium (|dT/dt| < `tol` K yr⁻¹).
    pub fn run_to_equilibrium(&mut self, dt: f64, tol: f64, max_steps: usize) -> usize {
        for step in 0..max_steps {
            let t_prev = self.temperature;
            self.step(dt, 0.0);
            if (self.temperature - t_prev).abs() / dt < tol {
                return step + 1;
            }
        }
        max_steps
    }
    /// Equilibrium temperature for a given albedo and forcing (°C).
    pub fn equilibrium_temperature(&self, albedo: f64, forcing: f64) -> f64 {
        let asw = (SOLAR_CONSTANT / 4.0) * (1.0 - albedo);
        (asw - self.olr_a + forcing) / self.olr_b
    }
}
/// Slab-ocean heat uptake model.
///
/// The mixed layer with depth *h* exchanges heat with the deep ocean at rate κ,
/// and with the atmosphere via the air-sea flux *F*.
#[derive(Debug, Clone)]
pub struct OceanHeatUptake {
    /// Mixed-layer depth (m).
    pub mixed_layer_depth: f64,
    /// Ocean temperature (°C).
    pub temperature: f64,
    /// Deep-ocean temperature (°C).
    pub deep_temperature: f64,
    /// Ocean heat diffusivity (m² s⁻¹).
    pub diffusivity: f64,
    /// Air-sea heat exchange coefficient (W m⁻² K⁻¹).
    pub air_sea_coefficient: f64,
    /// Atmospheric temperature driving the flux (°C).
    pub atm_temperature: f64,
}
impl OceanHeatUptake {
    /// Create with 50 m mixed layer.
    pub fn new() -> Self {
        Self {
            mixed_layer_depth: 50.0,
            temperature: 15.0,
            deep_temperature: 4.0,
            diffusivity: 1e-4,
            air_sea_coefficient: 20.0,
            atm_temperature: 15.0,
        }
    }
    /// Heat capacity of the mixed layer per unit area (J m⁻² K⁻¹).
    pub fn heat_capacity(&self) -> f64 {
        RHO_SEAWATER * CP_SEAWATER * self.mixed_layer_depth
    }
    /// Air-sea heat flux (W m⁻²) — positive into ocean.
    pub fn air_sea_flux(&self) -> f64 {
        self.air_sea_coefficient * (self.atm_temperature - self.temperature)
    }
    /// Entrainment flux to deep ocean (W m⁻²).
    pub fn deep_ocean_flux(&self) -> f64 {
        RHO_SEAWATER * CP_SEAWATER * self.diffusivity * (self.temperature - self.deep_temperature)
            / self.mixed_layer_depth
    }
    /// Net heat flux into mixed layer (W m⁻²).
    pub fn net_flux(&self) -> f64 {
        self.air_sea_flux() - self.deep_ocean_flux()
    }
    /// Integrate by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let c = self.heat_capacity();
        self.temperature += dt * self.net_flux() / c;
    }
    /// Ocean heat content anomaly relative to initial 15 °C (J m⁻²).
    pub fn heat_content_anomaly(&self, t_ref: f64) -> f64 {
        self.heat_capacity() * (self.temperature - t_ref)
    }
}
/// Sea-level rise model combining thermal expansion and ice-sheet melt.
///
/// Uses simple parameterisations calibrated to IPCC AR6 projections.
#[derive(Debug, Clone)]
pub struct SeaLevelRiseModel {
    /// Thermal expansion coefficient (m K⁻¹ per century).
    pub thermal_expansion_coeff: f64,
    /// Greenland ice-sheet contribution per degree (m K⁻¹ per century).
    pub greenland_contribution: f64,
    /// West Antarctic Ice Sheet contribution per degree (m K⁻¹ per century).
    pub wais_contribution: f64,
    /// Mountain glacier contribution per degree (m K⁻¹ per century).
    pub glacier_contribution: f64,
    /// Cumulative sea-level rise (m).
    pub sea_level: f64,
    /// Reference temperature (°C).
    pub reference_temperature: f64,
}
impl SeaLevelRiseModel {
    /// Create with IPCC AR6 mid-range coefficients.
    pub fn new() -> Self {
        Self {
            thermal_expansion_coeff: 0.38,
            greenland_contribution: 0.10,
            wais_contribution: 0.06,
            glacier_contribution: 0.18,
            sea_level: 0.0,
            reference_temperature: 14.0,
        }
    }
    /// Sea-level rise rate (m yr⁻¹) for a given temperature anomaly (K).
    pub fn rate(&self, delta_t: f64) -> f64 {
        (self.thermal_expansion_coeff
            + self.greenland_contribution
            + self.wais_contribution
            + self.glacier_contribution)
            * delta_t
            / 100.0
    }
    /// Advance sea level by `dt` years with temperature anomaly `delta_t` (K).
    pub fn step(&mut self, dt: f64, delta_t: f64) {
        self.sea_level += self.rate(delta_t) * dt;
    }
    /// Project sea-level rise to end year.
    pub fn project(
        &mut self,
        start_year: f64,
        end_year: f64,
        dt: f64,
        temperature_fn: impl Fn(f64) -> f64,
    ) -> f64 {
        let mut year = start_year;
        while year < end_year {
            let delta_t = temperature_fn(year) - self.reference_temperature;
            self.step(dt, delta_t);
            year += dt;
        }
        self.sea_level
    }
    /// Ice-sheet contribution fraction of total.
    pub fn ice_fraction(&self) -> f64 {
        let ice = self.greenland_contribution + self.wais_contribution + self.glacier_contribution;
        let total = ice + self.thermal_expansion_coeff;
        if total > 0.0 { ice / total } else { 0.0 }
    }
}
/// Simple two-box ENSO proxy (Zebiak-Cane reduced model).
///
/// Tracks west-Pacific (WP) and east-Pacific (EP) ocean temperatures.
/// The SST anomaly in the east Pacific modulates the Walker circulation.
#[derive(Debug, Clone)]
pub struct EnsoProxy {
    /// West-Pacific SST anomaly (°C).
    pub sst_west: f64,
    /// East-Pacific SST anomaly (°C).
    pub sst_east: f64,
    /// Thermocline depth anomaly (m).
    pub thermocline: f64,
    /// Coupling coefficient (K m⁻¹).
    pub coupling: f64,
    /// Damping rate (yr⁻¹).
    pub damping: f64,
    /// Delayed oscillator period (yr).
    pub period: f64,
    /// Internal phase (rad).
    pub phase: f64,
}
impl EnsoProxy {
    /// Create with typical ENSO parameters.
    pub fn new() -> Self {
        Self {
            sst_west: 0.0,
            sst_east: 0.0,
            thermocline: 0.0,
            coupling: 0.1,
            damping: 0.3,
            period: 4.0,
            phase: 0.0,
        }
    }
    /// Advance ENSO state by `dt` years.
    pub fn step(&mut self, dt: f64) {
        let omega = 2.0 * PI / self.period;
        self.phase += omega * dt;
        self.sst_east = 2.5 * self.phase.sin() * (-self.damping * self.phase / omega).exp();
        self.sst_west = -self.sst_east;
        self.thermocline = 30.0 * self.phase.cos();
    }
    /// ENSO index (Niño 3.4): east-Pacific SST anomaly (°C).
    pub fn nino34(&self) -> f64 {
        self.sst_east
    }
    /// Return phase classification.
    pub fn phase_name(&self) -> &'static str {
        let n = self.nino34();
        if n > 0.5 {
            "El Nino"
        } else if n < -0.5 {
            "La Nina"
        } else {
            "Neutral"
        }
    }
}
/// Aggregated radiative forcing from multiple agents (W m⁻²).
///
/// Includes CO₂, CH₄, N₂O, tropospheric ozone, aerosols, solar, and
/// land-use change.
#[derive(Debug, Clone)]
pub struct RadiativeForcing {
    /// CO₂ forcing (W m⁻²).
    pub co2: f64,
    /// CH₄ forcing (W m⁻²).
    pub ch4: f64,
    /// N₂O forcing (W m⁻²).
    pub n2o: f64,
    /// Tropospheric ozone forcing (W m⁻²).
    pub ozone: f64,
    /// Direct aerosol forcing (W m⁻²) — typically negative.
    pub aerosol_direct: f64,
    /// Indirect aerosol (cloud) forcing (W m⁻²) — typically negative.
    pub aerosol_indirect: f64,
    /// Solar irradiance change (W m⁻²).
    pub solar: f64,
    /// Land-use change forcing (W m⁻²).
    pub land_use: f64,
}
impl RadiativeForcing {
    /// Create a zero-forcing baseline.
    pub fn zero() -> Self {
        Self {
            co2: 0.0,
            ch4: 0.0,
            n2o: 0.0,
            ozone: 0.0,
            aerosol_direct: 0.0,
            aerosol_indirect: 0.0,
            solar: 0.0,
            land_use: 0.0,
        }
    }
    /// Compute CO₂ forcing from concentration ratio (dimensionless).
    pub fn from_co2_ratio(ratio: f64) -> f64 {
        5.35 * ratio.ln()
    }
    /// Compute CH₄ forcing from current (ppb) and reference (ppb) concentrations.
    ///
    /// Uses simplified expression: 0.036·(√M − √M₀) with overlap correction.
    pub fn from_ch4(m_ppb: f64, m0_ppb: f64, n2o_ppb: f64) -> f64 {
        let alpha = 0.47
            * (1.0
                + 2.01e-5 * (m_ppb * n2o_ppb).powf(0.75)
                + 5.31e-15 * m_ppb * (m_ppb * n2o_ppb).powf(1.52))
            .ln();
        let alpha0 = 0.47
            * (1.0
                + 2.01e-5 * (m0_ppb * n2o_ppb).powf(0.75)
                + 5.31e-15 * m0_ppb * (m0_ppb * n2o_ppb).powf(1.52))
            .ln();
        0.036 * (m_ppb.sqrt() - m0_ppb.sqrt()) - (alpha - alpha0)
    }
    /// Total effective radiative forcing (W m⁻²).
    pub fn total(&self) -> f64 {
        self.co2
            + self.ch4
            + self.n2o
            + self.ozone
            + self.aerosol_direct
            + self.aerosol_indirect
            + self.solar
            + self.land_use
    }
    /// Anthropogenic forcing only (excludes solar).
    pub fn anthropogenic(&self) -> f64 {
        self.total() - self.solar
    }
}
/// Three-reservoir carbon cycle model (atmosphere, ocean, land biosphere).
///
/// Atmosphere ↔ Ocean exchange follows a first-order piston-velocity model;
/// Atmosphere ↔ Land follows NPP-respiration balance with CO₂ fertilisation.
#[derive(Debug, Clone)]
pub struct CarbonCycle {
    /// Atmospheric carbon stock (GtC).
    pub atmosphere: f64,
    /// Surface ocean dissolved inorganic carbon stock (GtC).
    pub ocean_surface: f64,
    /// Deep ocean DIC stock (GtC).
    pub ocean_deep: f64,
    /// Terrestrial vegetation stock (GtC).
    pub land_vegetation: f64,
    /// Soil / litter stock (GtC).
    pub land_soil: f64,
    /// Air-sea gas exchange piston velocity (GtC yr⁻¹ ppm⁻¹).
    pub piston_velocity: f64,
    /// Ocean vertical mixing coefficient (GtC yr⁻¹ per GtC difference).
    pub ocean_mixing: f64,
    /// Net primary productivity at reference CO₂ (GtC yr⁻¹).
    pub npp_ref: f64,
    /// CO₂ fertilisation factor β (dimensionless).
    pub beta: f64,
    /// Heterotrophic respiration rate constant (yr⁻¹).
    pub respiration_rate: f64,
    /// Anthropogenic emission rate (GtC yr⁻¹).
    pub emission_rate: f64,
}
impl CarbonCycle {
    /// Create with approximate year-2000 stocks.
    pub fn new() -> Self {
        Self {
            atmosphere: 590.0,
            ocean_surface: 900.0,
            ocean_deep: 37_100.0,
            land_vegetation: 550.0,
            land_soil: 1500.0,
            piston_velocity: 0.24,
            ocean_mixing: 3.65,
            npp_ref: 60.0,
            beta: 0.35,
            respiration_rate: 0.04,
            emission_rate: 9.5,
        }
    }
    /// Atmospheric CO₂ in ppm from GtC.
    ///
    /// Conversion: 1 ppm ≈ 2.12 GtC.
    pub fn co2_ppm(&self) -> f64 {
        self.atmosphere / 2.12
    }
    /// Air-sea flux (GtC yr⁻¹) — positive means into ocean.
    pub fn air_sea_flux(&self) -> f64 {
        self.piston_velocity * (self.co2_ppm() - self.ocean_surface / 45.0)
    }
    /// Net primary productivity with CO₂ fertilisation (GtC yr⁻¹).
    pub fn npp(&self) -> f64 {
        let co2_ratio = self.co2_ppm() / CO2_PREINDUSTRIAL;
        self.npp_ref * (1.0 + self.beta * co2_ratio.ln())
    }
    /// Heterotrophic respiration from soil (GtC yr⁻¹).
    pub fn respiration(&self) -> f64 {
        self.respiration_rate * self.land_soil
    }
    /// Integrate carbon cycle by `dt` years.
    pub fn step(&mut self, dt: f64) {
        let f_as = self.air_sea_flux();
        let f_mix = self.ocean_mixing * (self.ocean_surface - self.ocean_deep / 41.2);
        let npp = self.npp();
        let resp = self.respiration();
        let da_atm = (self.emission_rate - f_as - (npp - resp)) * dt;
        let da_os = (f_as - f_mix) * dt;
        let da_od = f_mix * dt;
        let da_veg = (npp - resp) * dt;
        let da_soil = (resp - npp) * dt;
        self.atmosphere += da_atm;
        self.ocean_surface += da_os;
        self.ocean_deep += da_od;
        self.land_vegetation += da_veg;
        self.land_soil += da_soil;
    }
    /// Total carbon in the system (GtC) — should be approximately conserved.
    pub fn total_carbon(&self) -> f64 {
        self.atmosphere
            + self.ocean_surface
            + self.ocean_deep
            + self.land_vegetation
            + self.land_soil
    }
}
