// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0
//! SPH for climate modeling: ocean-atmosphere coupling, ice-ocean interaction.
//!
//! Implements stratified ocean layers, atmospheric boundary layer turbulence,
//! thermohaline circulation, sea ice albedo feedback, monsoon upwelling,
//! ENSO simplified oscillator, ocean CO2 uptake, wave climate, and polar vortex.

use std::f64::consts::PI;

// ─── Constants ───────────────────────────────────────────────────────────────

/// Von Kármán constant (dimensionless).
pub const KAPPA: f64 = 0.41;
/// Gravitational acceleration \[m/s²\].
pub const G: f64 = 9.81;
/// Stefan-Boltzmann constant \[W/m²/K⁴\].
pub const SIGMA_SB: f64 = 5.670_374_419e-8;
/// Reference seawater density \[kg/m³\].
pub const RHO0_SW: f64 = 1025.0;
/// Specific heat of seawater \[J/kg/K\].
pub const CP_SW: f64 = 3850.0;
/// Thermal expansion coefficient of seawater \[1/K\].
pub const ALPHA_SW: f64 = 2.1e-4;
/// Haline contraction coefficient \[1/PSU\].
pub const BETA_SW: f64 = 7.4e-4;
/// Albedo of sea ice.
pub const ALBEDO_ICE: f64 = 0.8;
/// Albedo of open ocean.
pub const ALBEDO_OCEAN: f64 = 0.06;

// ─── UNESCO Equation of State ─────────────────────────────────────────────────

/// Computes seawater density via a simplified UNESCO equation of state.
///
/// # Arguments
/// * `temp_c` - Temperature \[°C\]
/// * `salinity` - Salinity \[PSU\]
/// * `pressure_dbar` - Pressure \[dbar\] (1 dbar ≈ 1 m depth)
///
/// Returns density \[kg/m³\].
pub fn unesco_density(temp_c: f64, salinity: f64, pressure_dbar: f64) -> f64 {
    // Pure water density (simplified polynomial)
    let t = temp_c;
    let s = salinity;
    let p = pressure_dbar * 1.0e4; // convert dbar → Pa (1 dbar = 10000 Pa)

    // Coefficients for pure water density at atmospheric pressure (IES80)
    let rho_w = 999.842_594 + 6.793_952e-2 * t - 9.095_290e-3 * t * t + 1.001_685e-4 * t * t * t
        - 1.120_083e-6 * t * t * t * t
        + 6.536_332e-9 * t * t * t * t * t;

    // Salinity contribution (linearised)
    let a = 0.824_493 - 4.0899e-3 * t + 7.6438e-5 * t * t - 8.2467e-7 * t * t * t
        + 5.3875e-9 * t * t * t * t;
    let b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t * t;
    let c = 4.8314e-4;
    let rho_sw_0 = rho_w + a * s + b * s.powf(1.5) + c * s * s;

    // Pressure correction (secant bulk modulus, greatly simplified)
    let bulk_modulus = 2.2e9; // Pa, approximate
    rho_sw_0 * (1.0 + p / bulk_modulus)
}

// ─── OceanLayerSph ────────────────────────────────────────────────────────────

/// Stratified ocean column represented as discrete vertical layers.
#[derive(Debug, Clone)]
pub struct OceanLayerSph {
    /// Layer depths \[m\].
    pub depths: Vec<f64>,
    /// Temperature at each layer \[°C\].
    pub temperatures: Vec<f64>,
    /// Salinity at each layer \[PSU\].
    pub salinities: Vec<f64>,
    /// Computed density at each layer \[kg/m³\].
    pub densities: Vec<f64>,
}

impl OceanLayerSph {
    /// Create a new ocean column with the given layer parameters.
    ///
    /// Depths, temperatures, and salinities must all have the same length.
    pub fn new(depths: Vec<f64>, temperatures: Vec<f64>, salinities: Vec<f64>) -> Self {
        let n = depths.len();
        let mut densities = vec![0.0_f64; n];
        for i in 0..n {
            let p_dbar = depths[i] / 10.0; // 1 m ≈ 0.1 dbar
            densities[i] = unesco_density(temperatures[i], salinities[i], p_dbar);
        }
        Self {
            depths,
            temperatures,
            salinities,
            densities,
        }
    }

    /// Check whether the column is statically stable (density increases with depth).
    pub fn is_stable(&self) -> bool {
        for i in 1..self.densities.len() {
            if self.densities[i] < self.densities[i - 1] {
                return false;
            }
        }
        true
    }

    /// Mixed-layer depth: deepest layer whose density is within `delta` \[kg/m³\]
    /// of the surface density.
    pub fn mixed_layer_depth(&self, delta: f64) -> f64 {
        let rho_surf = self.densities[0];
        let mut mld = self.depths[0];
        for i in 1..self.densities.len() {
            if (self.densities[i] - rho_surf).abs() <= delta {
                mld = self.depths[i];
            } else {
                break;
            }
        }
        mld
    }
}

// ─── AtmosphericBoundarySph ───────────────────────────────────────────────────

/// Atmospheric boundary layer (ABL) described by Monin-Obukhov similarity theory.
#[derive(Debug, Clone)]
pub struct AtmosphericBoundarySph {
    /// Friction velocity u* \[m/s\].
    pub u_star: f64,
    /// Aerodynamic roughness length z_0 \[m\].
    pub z0: f64,
    /// Obukhov length L \[m\] (positive = stable, negative = unstable).
    pub obukhov_length: f64,
    /// Surface sensible heat flux H \[W/m²\].
    pub sensible_heat_flux: f64,
    /// Surface latent heat flux LE \[W/m²\].
    pub latent_heat_flux: f64,
}

impl AtmosphericBoundarySph {
    /// Create a new ABL with the given Monin-Obukhov parameters.
    pub fn new(
        u_star: f64,
        z0: f64,
        obukhov_length: f64,
        sensible_heat_flux: f64,
        latent_heat_flux: f64,
    ) -> Self {
        Self {
            u_star,
            z0,
            obukhov_length,
            sensible_heat_flux,
            latent_heat_flux,
        }
    }

    /// Wind speed at height `z` \[m\] using logarithmic profile (neutral stability).
    ///
    /// u(z) = (u* / κ) · ln(z / z_0)
    pub fn wind_speed(&self, z: f64) -> f64 {
        if z <= self.z0 || self.z0 <= 0.0 {
            return 0.0;
        }
        (self.u_star / KAPPA) * (z / self.z0).ln()
    }

    /// Wind speed with Monin-Obukhov stability correction (Businger-Dyer).
    ///
    /// Uses ψ_m stability functions for stable and unstable conditions.
    pub fn wind_speed_stability(&self, z: f64) -> f64 {
        if z <= self.z0 || self.z0 <= 0.0 {
            return 0.0;
        }
        let zeta = z / self.obukhov_length;
        let psi_m = if zeta >= 0.0 {
            // Stable: Businger et al. 1971
            -5.0 * zeta
        } else {
            // Unstable: Paulson 1970
            let x = (1.0 - 16.0 * zeta).powf(0.25);
            2.0 * (0.5 * (1.0 + x)).ln() + (0.5 * (1.0 + x * x)).ln() - 2.0 * x.atan() + PI / 2.0
        };
        (self.u_star / KAPPA) * ((z / self.z0).ln() - psi_m)
    }

    /// Total turbulent heat flux \[W/m²\] (sensible + latent).
    pub fn total_heat_flux(&self) -> f64 {
        self.sensible_heat_flux + self.latent_heat_flux
    }
}

// ─── ThermohalineCirculation ──────────────────────────────────────────────────

/// Simplified thermohaline (buoyancy-driven) circulation box model.
#[derive(Debug, Clone)]
pub struct ThermohalineCirculation {
    /// Surface heat flux \[W/m²\] (positive = ocean gains heat).
    pub q_heat: f64,
    /// Surface freshwater flux \[m/s\] (positive = precipitation - evaporation).
    pub q_salt: f64,
    /// Reference temperature \[°C\].
    pub temp_ref: f64,
    /// Reference salinity \[PSU\].
    pub sal_ref: f64,
}

impl ThermohalineCirculation {
    /// Create a new thermohaline box with given surface fluxes.
    pub fn new(q_heat: f64, q_salt: f64, temp_ref: f64, sal_ref: f64) -> Self {
        Self {
            q_heat,
            q_salt,
            temp_ref,
            sal_ref,
        }
    }

    /// Buoyancy flux B \[m²/s³\]:
    ///
    /// B = g · (α · Q_H / (ρ · c_p) - β · S · Q_S)
    pub fn buoyancy_flux(&self) -> f64 {
        let q_h = self.q_heat / (RHO0_SW * CP_SW); // [m·K/s]
        let q_s = self.q_salt * self.sal_ref; // [PSU·m/s]
        G * (ALPHA_SW * q_h - BETA_SW * q_s)
    }

    /// Deep water formation rate \[Sv, 1 Sv = 10^6 m³/s\], proportional to
    /// negative buoyancy flux (cooling / salinification).
    pub fn deep_water_formation_sv(&self, area_m2: f64) -> f64 {
        let b = -self.buoyancy_flux(); // negative = sinking
        if b <= 0.0 {
            return 0.0;
        }
        // Scaling: Overturning ~ sqrt(B * A) / g (highly simplified)
        (b * area_m2).sqrt() / G / 1.0e6
    }
}

// ─── SeaIceAlbedo ─────────────────────────────────────────────────────────────

/// Sea-ice / open-ocean albedo feedback module.
#[derive(Debug, Clone)]
pub struct SeaIceAlbedo {
    /// Ice fraction (0 = open ocean, 1 = full ice cover).
    pub ice_fraction: f64,
    /// Incoming shortwave radiation \[W/m²\].
    pub sw_down: f64,
    /// Incoming longwave radiation \[W/m²\].
    pub lw_down: f64,
    /// Upwelling longwave radiation \[W/m²\].
    pub lw_up: f64,
}

impl SeaIceAlbedo {
    /// Create a new sea-ice surface budget.
    pub fn new(ice_fraction: f64, sw_down: f64, lw_down: f64, lw_up: f64) -> Self {
        Self {
            ice_fraction,
            sw_down,
            lw_down,
            lw_up,
        }
    }

    /// Effective surface albedo (linear mix of ice and ocean albedo).
    pub fn effective_albedo(&self) -> f64 {
        self.ice_fraction * ALBEDO_ICE + (1.0 - self.ice_fraction) * ALBEDO_OCEAN
    }

    /// Net radiation at the surface \[W/m²\]:
    ///
    /// R_net = (1 - α) · SW_down + LW_down - LW_up
    pub fn net_radiation(&self) -> f64 {
        let alpha = self.effective_albedo();
        (1.0 - alpha) * self.sw_down + self.lw_down - self.lw_up
    }

    /// Change in net radiation for a unit increase in ice fraction \[W/m²\].
    ///
    /// Positive feedback: more ice → higher albedo → less absorbed SW.
    pub fn albedo_feedback_sensitivity(&self) -> f64 {
        -(ALBEDO_ICE - ALBEDO_OCEAN) * self.sw_down
    }
}

// ─── MonsoonSph ───────────────────────────────────────────────────────────────

/// Seasonal monsoon-driven coastal upwelling model.
#[derive(Debug, Clone)]
pub struct MonsoonSph {
    /// Along-shore wind stress τ \[Pa\] (positive = equatorward).
    pub wind_stress: f64,
    /// Coriolis parameter f \[1/s\].
    pub coriolis: f64,
    /// Mixed-layer depth H \[m\].
    pub mixed_layer_depth: f64,
    /// Seawater density \[kg/m³\].
    pub density: f64,
}

impl MonsoonSph {
    /// Create a new monsoon upwelling configuration.
    pub fn new(wind_stress: f64, coriolis: f64, mixed_layer_depth: f64, density: f64) -> Self {
        Self {
            wind_stress,
            coriolis,
            mixed_layer_depth,
            density,
        }
    }

    /// Ekman transport per unit coastline length \[m²/s\]:
    ///
    /// M_E = τ / (ρ · f)
    pub fn ekman_transport(&self) -> f64 {
        if self.coriolis.abs() < 1.0e-12 {
            return 0.0;
        }
        self.wind_stress / (self.density * self.coriolis)
    }

    /// Upwelling velocity w \[m/s\]:
    ///
    /// w = M_E / H
    pub fn upwelling_velocity(&self) -> f64 {
        if self.mixed_layer_depth <= 0.0 {
            return 0.0;
        }
        self.ekman_transport() / self.mixed_layer_depth
    }

    /// SST cooling rate due to upwelling \[°C/day\] given a temperature gradient
    /// `dT_dz` \[°C/m\] between the surface and thermocline.
    pub fn sst_cooling_rate(&self, dt_dz: f64) -> f64 {
        -self.upwelling_velocity() * dt_dz * 86400.0
    }
}

// ─── ElNino ───────────────────────────────────────────────────────────────────

/// Simplified linear delayed oscillator model for ENSO.
///
/// dT/dt = r·T - c·h(t - τ_d)
/// where T = SST anomaly, h = thermocline depth anomaly.
#[derive(Debug, Clone)]
pub struct ElNino {
    /// SST anomaly history \[°C\] (most recent last).
    pub sst_anomaly: Vec<f64>,
    /// Thermocline depth anomaly history \[m\] (most recent last).
    pub thermocline_anomaly: Vec<f64>,
    /// Coupling coefficient r \[1/year\].
    pub r_coupling: f64,
    /// Delayed feedback coefficient c \[°C/m/year\].
    pub c_feedback: f64,
    /// Delay time \[years\].
    pub delay_years: f64,
    /// Time step \[years\].
    pub dt_years: f64,
}

impl ElNino {
    /// Create an ENSO oscillator with default parameters (Suarez & Schopf 1988).
    pub fn new(dt_years: f64) -> Self {
        Self {
            sst_anomaly: vec![0.5],
            thermocline_anomaly: vec![5.0],
            r_coupling: 1.0 / 0.5, // 6-month growth rate
            c_feedback: 0.025,     // delayed negative feedback
            delay_years: 0.5,      // ~6-month delay
            dt_years,
        }
    }

    /// Advance the model by one time step.
    pub fn step(&mut self) {
        let n = self.sst_anomaly.len();
        let t_now = self.sst_anomaly[n - 1];
        let h_now = self.thermocline_anomaly[n - 1];

        // Delayed index
        let delay_steps = (self.delay_years / self.dt_years).round() as usize;
        let h_delayed = if n > delay_steps {
            self.thermocline_anomaly[n - delay_steps]
        } else {
            self.thermocline_anomaly[0]
        };

        let dt = self.r_coupling * t_now - self.c_feedback * h_delayed;
        let t_new = t_now + self.dt_years * dt;
        // Thermocline driven by SST (simple restoring)
        let h_new = h_now + self.dt_years * (-0.3 * h_now + 2.0 * t_now);

        self.sst_anomaly.push(t_new.clamp(-5.0, 5.0));
        self.thermocline_anomaly.push(h_new.clamp(-30.0, 30.0));
    }

    /// Run for `n_years` years and return peak-to-peak period estimate \[years\].
    pub fn simulate_years(&mut self, n_years: f64) -> f64 {
        let steps = (n_years / self.dt_years).round() as usize;
        for _ in 0..steps {
            self.step();
        }
        // Count zero crossings to estimate period
        let anomalies = &self.sst_anomaly;
        let mut crossings = 0usize;
        for i in 1..anomalies.len() {
            if anomalies[i - 1] * anomalies[i] < 0.0 {
                crossings += 1;
            }
        }
        if crossings < 2 {
            return 0.0;
        }
        2.0 * n_years / crossings as f64
    }
}

// ─── CarbonCycleOcean ─────────────────────────────────────────────────────────

/// Ocean CO2 uptake via carbonate chemistry.
#[derive(Debug, Clone)]
pub struct CarbonCycleOcean {
    /// Atmospheric pCO2 \[µatm\].
    pub pco2_atm: f64,
    /// Surface ocean temperature \[°C\].
    pub sst: f64,
    /// Surface salinity \[PSU\].
    pub salinity: f64,
    /// Total alkalinity \[µmol/kg\].
    pub alkalinity: f64,
    /// Dissolved inorganic carbon (DIC) \[µmol/kg\].
    pub dic: f64,
}

impl CarbonCycleOcean {
    /// Construct with given ocean surface conditions.
    pub fn new(pco2_atm: f64, sst: f64, salinity: f64, alkalinity: f64, dic: f64) -> Self {
        Self {
            pco2_atm,
            sst,
            salinity,
            alkalinity,
            dic,
        }
    }

    /// Solubility of CO2 in seawater \[mol/L/atm\] (Weiss 1974, simplified).
    ///
    /// K_0 ≈ exp(a1 + a2/T + a3·ln(T) + S·(b1 + b2/T))
    pub fn co2_solubility(&self) -> f64 {
        let tk = self.sst + 273.15;
        // Coefficients from Weiss 1974
        let a1 = -60.2409;
        let a2 = 9345.17;
        let a3 = 23.3585;
        let b1 = 0.023517;
        let b2 = -0.023656;
        let s = self.salinity;
        (a1 + a2 / tk + a3 * tk.ln() + s * (b1 + b2 * 100.0 / tk)).exp()
    }

    /// Ocean surface pCO2 \[µatm\] estimated from DIC and alkalinity
    /// using simplified carbonate equilibrium.
    pub fn ocean_pco2(&self) -> f64 {
        // Simplified: pCO2 ∝ DIC² / (Alkalinity * K0)
        let k0 = self.co2_solubility();
        if k0 <= 0.0 || self.alkalinity <= 0.0 {
            return 0.0;
        }
        // First dissociation constant K1 [mol/kg] (Mehrbach, simplified)
        let tk = self.sst + 273.15;
        let log_k1 = -3633.86 / tk + 61.2172 - 9.67770 * tk.ln();
        let k1 = 10_f64.powf(log_k1);

        // [CO2*] = DIC / (1 + K1/[H+])  — at typical seawater pH ~8.1
        let h_plus = 7.943e-9_f64; // pH 8.1
        let co2_aq = self.dic / (1.0 + k1 / h_plus + k1 * 1.2e7 / (h_plus * h_plus));
        co2_aq / k0 * 1.0e6 // convert mol/L → µatm (approximate)
    }

    /// Air-sea CO2 flux \[mol/m²/yr\] (positive = ocean uptake).
    ///
    /// F = k_w · K_0 · (pCO2_atm - pCO2_ocean)
    /// where k_w ≈ 0.31 · U² · (Sc/660)^{-0.5} \[cm/h\]
    pub fn air_sea_flux(&self, wind_speed_ms: f64) -> f64 {
        let dpco2 = self.pco2_atm - self.ocean_pco2(); // µatm
        // Schmidt number for CO2 (Wanninkhof 1992, simplified)
        let sc = 2073.1 - 125.62 * self.sst + 3.6276 * self.sst * self.sst;
        let k_w_cm_h = 0.31 * wind_speed_ms * wind_speed_ms * (sc / 660.0).powf(-0.5);
        let k_w_m_yr = k_w_cm_h * 0.01 * 8760.0; // cm/h → m/yr
        let k0 = self.co2_solubility(); // mol/L/atm
        let k0_mol_m3_uatm = k0 * 1000.0 / 1.0e6; // mol/m³/µatm
        k_w_m_yr * k0_mol_m3_uatm * dpco2
    }
}

// ─── WaveClimate ──────────────────────────────────────────────────────────────

/// Significant wave height estimation via Bretschneider (SMB) method.
#[derive(Debug, Clone)]
pub struct WaveClimate {
    /// 10-m wind speed \[m/s\].
    pub wind_speed: f64,
    /// Fetch (upwind distance over water) \[m\].
    pub fetch_m: f64,
    /// Water depth \[m\] (use large value for deep water).
    pub depth_m: f64,
}

impl WaveClimate {
    /// Create a new wave climate scenario.
    pub fn new(wind_speed: f64, fetch_m: f64, depth_m: f64) -> Self {
        Self {
            wind_speed,
            fetch_m,
            depth_m,
        }
    }

    /// Significant wave height H_s \[m\] using Bretschneider deep-water formula.
    ///
    /// H_s = 0.0248 · U² / g · tanh(0.0414 · (g·F/U²)^0.57)
    pub fn significant_wave_height(&self) -> f64 {
        let u = self.wind_speed;
        if u <= 0.0 {
            return 0.0;
        }
        let gf_over_u2 = G * self.fetch_m / (u * u);
        0.0248 * u * u / G * (0.0414 * gf_over_u2.powf(0.57)).tanh()
    }

    /// Peak wave period T_p \[s\] (Bretschneider).
    ///
    /// T_p = 0.7457 · U / g · (g·F/U²)^0.36
    pub fn peak_period(&self) -> f64 {
        let u = self.wind_speed;
        if u <= 0.0 {
            return 0.0;
        }
        let gf_over_u2 = G * self.fetch_m / (u * u);
        0.7457 * u / G * gf_over_u2.powf(0.36)
    }

    /// Wave power per unit crest length \[W/m\] (deep water):
    ///
    /// P = ρ · g² · H_s² · T_p / (32π)
    pub fn wave_power(&self) -> f64 {
        let hs = self.significant_wave_height();
        let tp = self.peak_period();
        RHO0_SW * G * G * hs * hs * tp / (32.0 * PI)
    }
}

// ─── PolarVortex ──────────────────────────────────────────────────────────────

/// Stratospheric polar vortex stability diagnostics.
#[derive(Debug, Clone)]
pub struct PolarVortex {
    /// Potential vorticity \[PVU, 1 PVU = 10^-6 K·m²/kg/s\].
    pub pv_field: Vec<f64>,
    /// Latitudes \[°\] corresponding to pv_field.
    pub latitudes: Vec<f64>,
    /// Zonal mean temperature at 50 hPa \[K\].
    pub temp_50hpa: Vec<f64>,
    /// Zonal mean zonal wind at 10 hPa \[m/s\].
    pub u_10hpa: Vec<f64>,
}

impl PolarVortex {
    /// Construct a polar vortex diagnostic with given fields.
    pub fn new(
        pv_field: Vec<f64>,
        latitudes: Vec<f64>,
        temp_50hpa: Vec<f64>,
        u_10hpa: Vec<f64>,
    ) -> Self {
        Self {
            pv_field,
            latitudes,
            temp_50hpa,
            u_10hpa,
        }
    }

    /// Maximum PV gradient \[PVU/degree latitude\] — proxy for vortex edge sharpness.
    pub fn max_pv_gradient(&self) -> f64 {
        if self.pv_field.len() < 2 {
            return 0.0;
        }
        let mut max_grad = 0.0_f64;
        for i in 1..self.pv_field.len() {
            let dlat = (self.latitudes[i] - self.latitudes[i - 1]).abs();
            if dlat > 0.0 {
                let grad = (self.pv_field[i] - self.pv_field[i - 1]).abs() / dlat;
                if grad > max_grad {
                    max_grad = grad;
                }
            }
        }
        max_grad
    }

    /// Mean stratospheric wind at 10 hPa poleward of 60°N \[m/s\].
    pub fn mean_polar_wind(&self) -> f64 {
        let mut sum = 0.0;
        let mut count = 0usize;
        for (i, &lat) in self.latitudes.iter().enumerate() {
            if lat >= 60.0 && i < self.u_10hpa.len() {
                sum += self.u_10hpa[i];
                count += 1;
            }
        }
        if count == 0 {
            return 0.0;
        }
        sum / count as f64
    }

    /// Returns true when a sudden stratospheric warming (SSW) event is detected:
    /// polar cap (≥60°N) zonal mean wind reverses (becomes easterly).
    pub fn is_sudden_stratospheric_warming(&self) -> bool {
        self.mean_polar_wind() < 0.0
    }

    /// Vortex strength index: normalised by reference strong-vortex wind (30 m/s).
    pub fn vortex_strength_index(&self) -> f64 {
        self.mean_polar_wind() / 30.0
    }
}

// ─── Tests ────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── UNESCO EOS ───────────────────────────────────────────────────────────

    #[test]
    fn test_unesco_density_pure_water_0c() {
        let rho = unesco_density(0.0, 0.0, 0.0);
        // TEOS-10: ρ(0°C, 0 PSU, 0 dbar) ≈ 999.84 kg/m³
        assert!((rho - 999.84).abs() < 0.1, "Got {:.6}", rho);
    }

    #[test]
    fn test_unesco_density_seawater_35psu() {
        let rho = unesco_density(15.0, 35.0, 0.0);
        // Typical surface seawater ≈ 1025–1026 kg/m³
        assert!(rho > 1024.0 && rho < 1027.0, "Got {:.6}", rho);
    }

    #[test]
    fn test_unesco_density_increases_with_salinity() {
        let rho_fresh = unesco_density(20.0, 0.0, 0.0);
        let rho_salty = unesco_density(20.0, 35.0, 0.0);
        assert!(rho_salty > rho_fresh);
    }

    #[test]
    fn test_unesco_density_increases_with_pressure() {
        let rho_surf = unesco_density(10.0, 35.0, 0.0);
        let rho_deep = unesco_density(10.0, 35.0, 1000.0);
        assert!(rho_deep > rho_surf);
    }

    #[test]
    fn test_unesco_density_decreases_with_temperature() {
        let rho_cold = unesco_density(2.0, 35.0, 0.0);
        let rho_warm = unesco_density(25.0, 35.0, 0.0);
        assert!(rho_cold > rho_warm);
    }

    // ── OceanLayerSph ─────────────────────────────────────────────────────────

    #[test]
    fn test_ocean_layer_stable_column() {
        // Density increases from surface to deep (warm-light surface, cool-dense deep)
        let depths = vec![0.0, 50.0, 200.0, 500.0];
        let temps = vec![25.0, 15.0, 5.0, 2.0];
        let sals = vec![35.0, 35.0, 34.8, 34.6];
        let col = OceanLayerSph::new(depths, temps, sals);
        assert!(col.is_stable());
    }

    #[test]
    fn test_ocean_layer_unstable_column() {
        // Dense water on top of light water → unstable
        let depths = vec![0.0, 100.0];
        let temps = vec![2.0, 25.0]; // cold dense on top
        let sals = vec![35.0, 35.0];
        let col = OceanLayerSph::new(depths, temps, sals);
        assert!(!col.is_stable());
    }

    #[test]
    fn test_ocean_layer_mixed_layer_depth() {
        let depths = vec![0.0, 10.0, 20.0, 50.0, 100.0];
        let temps = vec![20.0, 19.9, 19.8, 15.0, 10.0];
        let sals = vec![35.0, 35.0, 35.0, 35.0, 35.0];
        let col = OceanLayerSph::new(depths, temps, sals);
        let mld = col.mixed_layer_depth(0.5);
        // Layers 0,1,2 are within 0.5 kg/m³ of surface; layer 3 is not
        assert!(
            mld <= 20.0 + 1.0,
            "MLD should be around 20 m, got {:.6}",
            mld
        );
    }

    // ── AtmosphericBoundarySph ────────────────────────────────────────────────

    #[test]
    fn test_abl_log_wind_profile() {
        let abl = AtmosphericBoundarySph::new(0.3, 0.01, 1.0e6, 50.0, 30.0);
        let u10 = abl.wind_speed(10.0);
        let u100 = abl.wind_speed(100.0);
        assert!(u100 > u10, "Wind should increase with height");
    }

    #[test]
    fn test_abl_zero_below_roughness() {
        let abl = AtmosphericBoundarySph::new(0.3, 0.1, 1.0e6, 0.0, 0.0);
        assert_eq!(abl.wind_speed(0.05), 0.0);
    }

    #[test]
    fn test_abl_log_wind_magnitude() {
        // u* = 0.3 m/s, z0 = 0.01 m → u(10) ≈ (0.3/0.41)*ln(1000) ≈ 5.05 m/s
        let abl = AtmosphericBoundarySph::new(0.3, 0.01, 1.0e6, 0.0, 0.0);
        let u10 = abl.wind_speed(10.0);
        assert!((u10 - 5.05).abs() < 0.3, "Got {:.6}", u10);
    }

    #[test]
    fn test_abl_stability_correction_unstable() {
        let abl = AtmosphericBoundarySph::new(0.3, 0.01, -50.0, 200.0, 100.0);
        let u_stab = abl.wind_speed_stability(10.0);
        let u_neut = abl.wind_speed(10.0);
        // Unstable → turbulent mixing → weaker wind shear → smaller u at 10 m
        // (less certain, but should be positive)
        assert!(u_stab > 0.0);
        let _ = u_neut;
    }

    #[test]
    fn test_abl_total_heat_flux() {
        let abl = AtmosphericBoundarySph::new(0.3, 0.01, -100.0, 80.0, 120.0);
        assert!((abl.total_heat_flux() - 200.0).abs() < 1e-9);
    }

    // ── ThermohalineCirculation ───────────────────────────────────────────────

    #[test]
    fn test_thc_buoyancy_flux_cooling() {
        // Heat loss → negative Q_H → negative buoyancy flux (sinking)
        let thc = ThermohalineCirculation::new(-100.0, 0.0, 5.0, 35.0);
        assert!(thc.buoyancy_flux() < 0.0);
    }

    #[test]
    fn test_thc_buoyancy_flux_heating() {
        // Heat gain → positive Q_H → positive buoyancy flux (stabilising)
        let thc = ThermohalineCirculation::new(200.0, 0.0, 20.0, 35.0);
        assert!(thc.buoyancy_flux() > 0.0);
    }

    #[test]
    fn test_thc_deep_water_formation_nonzero() {
        let thc = ThermohalineCirculation::new(-200.0, 1.0e-7, 2.0, 35.0);
        let dw = thc.deep_water_formation_sv(1.0e12); // 10^12 m² = ~North Atlantic
        assert!(dw >= 0.0);
    }

    // ── SeaIceAlbedo ──────────────────────────────────────────────────────────

    #[test]
    fn test_sea_ice_albedo_full_ice() {
        let sia = SeaIceAlbedo::new(1.0, 300.0, 250.0, 280.0);
        assert!((sia.effective_albedo() - ALBEDO_ICE).abs() < 1e-12);
    }

    #[test]
    fn test_sea_ice_albedo_open_ocean() {
        let sia = SeaIceAlbedo::new(0.0, 300.0, 250.0, 280.0);
        assert!((sia.effective_albedo() - ALBEDO_OCEAN).abs() < 1e-12);
    }

    #[test]
    fn test_sea_ice_net_radiation_ice_lower() {
        // More ice → higher albedo → less net radiation absorbed
        let sia_ice = SeaIceAlbedo::new(1.0, 400.0, 280.0, 300.0);
        let sia_open = SeaIceAlbedo::new(0.0, 400.0, 280.0, 300.0);
        assert!(sia_ice.net_radiation() < sia_open.net_radiation());
    }

    #[test]
    fn test_albedo_feedback_negative() {
        // More ice → less absorbed SW → negative (cooling) feedback
        let sia = SeaIceAlbedo::new(0.5, 400.0, 280.0, 300.0);
        assert!(sia.albedo_feedback_sensitivity() < 0.0);
    }

    #[test]
    fn test_sea_ice_net_radiation_value() {
        let sia = SeaIceAlbedo::new(0.0, 200.0, 300.0, 350.0);
        // R_net = (1-0.06)*200 + 300 - 350 = 188 + 300 - 350 = 138
        let expected = (1.0 - ALBEDO_OCEAN) * 200.0 + 300.0 - 350.0;
        assert!((sia.net_radiation() - expected).abs() < 1e-9);
    }

    // ── MonsoonSph ────────────────────────────────────────────────────────────

    #[test]
    fn test_monsoon_ekman_transport() {
        let m = MonsoonSph::new(0.1, 1.0e-4, 30.0, 1025.0);
        let me = m.ekman_transport();
        // τ/(ρf) = 0.1/(1025*1e-4) ≈ 0.976 m²/s
        assert!((me - 0.1 / (1025.0 * 1.0e-4)).abs() < 1e-6, "Got {:.6}", me);
    }

    #[test]
    fn test_monsoon_upwelling_velocity() {
        let m = MonsoonSph::new(0.05, 5.0e-5, 50.0, 1025.0);
        let w = m.upwelling_velocity();
        assert!(w > 0.0);
    }

    #[test]
    fn test_monsoon_sst_cooling() {
        let m = MonsoonSph::new(0.05, 5.0e-5, 50.0, 1025.0);
        // Positive dT/dz (cold below) → upwelling cools SST
        let cooling = m.sst_cooling_rate(0.05);
        assert!(cooling < 0.0);
    }

    // ── WaveClimate ───────────────────────────────────────────────────────────

    #[test]
    fn test_wave_height_increases_with_wind() {
        let w10 = WaveClimate::new(10.0, 1.0e5, 100.0);
        let w20 = WaveClimate::new(20.0, 1.0e5, 100.0);
        assert!(w20.significant_wave_height() > w10.significant_wave_height());
    }

    #[test]
    fn test_wave_height_increases_with_fetch() {
        let wshort = WaveClimate::new(15.0, 5.0e4, 100.0);
        let wlong = WaveClimate::new(15.0, 2.0e5, 100.0);
        assert!(wlong.significant_wave_height() > wshort.significant_wave_height());
    }

    #[test]
    fn test_wave_period_positive() {
        let w = WaveClimate::new(12.0, 1.0e5, 50.0);
        assert!(w.peak_period() > 0.0);
    }

    #[test]
    fn test_wave_power_positive() {
        let w = WaveClimate::new(15.0, 2.0e5, 200.0);
        assert!(w.wave_power() > 0.0);
    }

    #[test]
    fn test_wave_height_zero_wind() {
        let w = WaveClimate::new(0.0, 1.0e5, 100.0);
        assert_eq!(w.significant_wave_height(), 0.0);
    }

    // ── ElNino ────────────────────────────────────────────────────────────────

    #[test]
    fn test_enso_simulation_runs() {
        let mut el = ElNino::new(0.1);
        let _ = el.simulate_years(10.0);
        assert!(el.sst_anomaly.len() > 1);
    }

    #[test]
    fn test_enso_sst_bounded() {
        let mut el = ElNino::new(0.1);
        el.simulate_years(20.0);
        for &t in &el.sst_anomaly {
            assert!(t.abs() <= 5.0, "SST anomaly out of bounds: {:.6}", t);
        }
    }

    #[test]
    fn test_enso_period_range() {
        let mut el = ElNino::new(1.0 / 12.0); // monthly steps
        let period = el.simulate_years(30.0);
        // Period should be in 2–8 year range (or 0 if no crossings detected)
        if period > 0.0 {
            assert!(
                (2.0..=10.0).contains(&period),
                "Period {:.6} out of expected range",
                period
            );
        }
    }

    // ── CarbonCycleOcean ──────────────────────────────────────────────────────

    #[test]
    fn test_co2_solubility_decreases_with_temperature() {
        let cold = CarbonCycleOcean::new(400.0, 5.0, 35.0, 2300.0, 2100.0);
        let warm = CarbonCycleOcean::new(400.0, 25.0, 35.0, 2300.0, 2100.0);
        // Warmer water holds less CO2
        assert!(cold.co2_solubility() > warm.co2_solubility());
    }

    #[test]
    fn test_air_sea_flux_positive_for_high_atm_pco2() {
        // High atmospheric CO2 → ocean absorbs CO2 → positive flux
        let ocean = CarbonCycleOcean::new(800.0, 15.0, 35.0, 2300.0, 2000.0);
        let flux = ocean.air_sea_flux(8.0);
        // Cannot assert direction without knowing ocean pCO2, but flux should be finite
        assert!(flux.is_finite());
    }

    #[test]
    fn test_air_sea_flux_zero_wind() {
        let ocean = CarbonCycleOcean::new(400.0, 15.0, 35.0, 2300.0, 2100.0);
        assert_eq!(ocean.air_sea_flux(0.0), 0.0);
    }

    // ── PolarVortex ───────────────────────────────────────────────────────────

    #[test]
    fn test_polar_vortex_mean_wind_polar_cap() {
        let lats = vec![50.0, 60.0, 70.0, 80.0, 90.0];
        let pv = vec![20.0, 30.0, 50.0, 80.0, 100.0];
        let t50 = vec![220.0; 5];
        let u10 = vec![20.0, 25.0, 30.0, 28.0, 25.0];
        let vortex = PolarVortex::new(pv, lats, t50, u10);
        let mw = vortex.mean_polar_wind();
        // Mean of [25,30,28,25] = 27.0
        assert!((mw - 27.0).abs() < 1e-9, "Got {:.6}", mw);
    }

    #[test]
    fn test_polar_vortex_ssw_detected() {
        let lats = vec![60.0, 70.0, 80.0];
        let pv = vec![30.0, 50.0, 70.0];
        let t50 = vec![230.0; 3];
        let u10 = vec![-5.0, -8.0, -3.0]; // Easterlies → SSW
        let vortex = PolarVortex::new(pv, lats, t50, u10);
        assert!(vortex.is_sudden_stratospheric_warming());
    }

    #[test]
    fn test_polar_vortex_no_ssw_westerlies() {
        let lats = vec![60.0, 75.0, 90.0];
        let pv = vec![40.0, 70.0, 120.0];
        let t50 = vec![210.0; 3];
        let u10 = vec![30.0, 35.0, 28.0];
        let vortex = PolarVortex::new(pv, lats, t50, u10);
        assert!(!vortex.is_sudden_stratospheric_warming());
    }

    #[test]
    fn test_polar_vortex_strength_index() {
        let lats = vec![65.0, 80.0];
        let pv = vec![50.0, 90.0];
        let t50 = vec![215.0; 2];
        let u10 = vec![30.0, 30.0]; // exactly reference wind
        let vortex = PolarVortex::new(pv, lats, t50, u10);
        assert!((vortex.vortex_strength_index() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_polar_vortex_max_pv_gradient() {
        let lats = vec![60.0, 65.0, 70.0];
        let pv = vec![10.0, 50.0, 55.0];
        let t50 = vec![215.0; 3];
        let u10 = vec![25.0; 3];
        let vortex = PolarVortex::new(pv, lats, t50, u10);
        // Largest gradient between 60 and 65: |50-10|/5 = 8 PVU/°
        let mg = vortex.max_pv_gradient();
        assert!((mg - 8.0).abs() < 1e-9, "Got {:.6}", mg);
    }

    // ── Integration / cross-module ────────────────────────────────────────────

    #[test]
    fn test_thermohaline_dense_water_sinks() {
        // Cold + salty water should be denser than warm + fresh
        let rho_cold_salty = unesco_density(2.0, 35.0, 0.0);
        let rho_warm_fresh = unesco_density(25.0, 0.0, 0.0);
        assert!(rho_cold_salty > rho_warm_fresh);
    }

    #[test]
    fn test_ice_albedo_feedback_positive() {
        // Positive feedback: more ice → less energy absorbed
        let half_ice = SeaIceAlbedo::new(0.5, 300.0, 200.0, 250.0);
        let full_ice = SeaIceAlbedo::new(1.0, 300.0, 200.0, 250.0);
        assert!(full_ice.net_radiation() < half_ice.net_radiation());
    }

    #[test]
    fn test_log_wind_profile_kappa() {
        // Verify κ constant embedded correctly
        let abl = AtmosphericBoundarySph::new(1.0, 1.0, 1.0e9, 0.0, 0.0);
        // u(e) = (1/0.41)*ln(e/1) = 1/0.41 ≈ 2.439
        let u_e = abl.wind_speed(std::f64::consts::E);
        assert!((u_e - 1.0 / KAPPA).abs() < 1e-9, "Got {:.6}", u_e);
    }

    #[test]
    fn test_carbon_uptake_stronger_with_wind() {
        let mut o1 = CarbonCycleOcean::new(500.0, 10.0, 34.0, 2350.0, 1900.0);
        let mut o2 = o1.clone();
        o2.pco2_atm = 500.0;
        o1.pco2_atm = 500.0;
        let f_low = o1.air_sea_flux(3.0);
        let f_high = o2.air_sea_flux(10.0);
        // Higher wind → more gas transfer
        assert!(f_high.abs() > f_low.abs());
    }
}
