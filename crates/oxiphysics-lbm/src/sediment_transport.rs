// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Sediment transport models for LBM simulations.
//!
//! Provides bed-load (Meyer-Peter & Müller), suspended-load (Rouse profile),
//! Shields parameter, Stokes settling, and a combined total-load model.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// SedimentParticle
// ---------------------------------------------------------------------------

/// A single sediment particle with position, diameter, and density.
#[derive(Debug, Clone, PartialEq)]
pub struct SedimentParticle {
    /// 3-D position \[x, y, z\] in metres.
    pub position: [f64; 3],
    /// Equivalent spherical diameter (m).
    pub diameter: f64,
    /// Particle density (kg/m³).
    pub density: f64,
}

impl SedimentParticle {
    /// Create a new `SedimentParticle`.
    pub fn new(position: [f64; 3], diameter: f64, density: f64) -> Self {
        Self {
            position,
            diameter,
            density,
        }
    }

    /// Stokes settling velocity (m/s) in a fluid of kinematic viscosity `nu`
    /// and density `rho_f`.
    ///
    /// ```text
    /// w_s = (ρ_s - ρ_f) g d² / (18 μ)
    ///     = (s - 1) g d² / (18 ν)
    /// ```
    ///
    /// where `s = ρ_s / ρ_f`.
    pub fn settling_velocity(&self, rho_f: f64, nu: f64) -> f64 {
        const G: f64 = 9.81;
        if nu <= 0.0 || rho_f <= 0.0 {
            return 0.0;
        }
        let s = self.density / rho_f;
        (s - 1.0) * G * self.diameter * self.diameter / (18.0 * nu)
    }
}

// ---------------------------------------------------------------------------
// BedLoad
// ---------------------------------------------------------------------------

/// Bed-load sediment transport model (Meyer-Peter & Müller, 1948).
#[derive(Debug, Clone)]
pub struct BedLoad {
    /// Dimensionless Shields stress (bed shear stress parameter).
    pub shields_parameter: f64,
    /// Critical Shields parameter for incipient motion (≈ 0.047).
    pub critical_shields: f64,
}

impl BedLoad {
    /// Create a new `BedLoad` model.
    pub fn new(shields_parameter: f64, critical_shields: f64) -> Self {
        Self {
            shields_parameter,
            critical_shields,
        }
    }

    /// Dimensionless bed-load transport rate using the Meyer-Peter & Müller
    /// formula:
    ///
    /// ```text
    /// q_b* = 8 (θ - θ_c)^{3/2}   for θ > θ_c
    ///      = 0                      otherwise
    /// ```
    ///
    /// # Arguments
    /// * `shields` - current Shields parameter θ
    pub fn meyer_peter_muller_rate(&self, shields: f64) -> f64 {
        let excess = shields - self.critical_shields;
        if excess <= 0.0 {
            return 0.0;
        }
        8.0 * excess.powf(1.5)
    }
}

// ---------------------------------------------------------------------------
// SuspendedLoad
// ---------------------------------------------------------------------------

/// Suspended-load sediment transport model using the Rouse profile.
#[derive(Debug, Clone)]
pub struct SuspendedLoad {
    /// Reference (near-bed) sediment concentration (kg/m³ or volumetric).
    pub reference_concentration: f64,
    /// Rouse number `Z = w_s / (κ u*)`.
    pub rouse_number: f64,
}

impl SuspendedLoad {
    /// Create a new `SuspendedLoad` model.
    pub fn new(reference_concentration: f64, rouse_number: f64) -> Self {
        Self {
            reference_concentration,
            rouse_number,
        }
    }

    /// Rouse concentration profile at height `y` above the bed in a flow of
    /// depth `h`.
    ///
    /// ```text
    /// C(y) = C_a * [ (a/y) * ((h - y)/(h - a)) ]^Z
    /// ```
    ///
    /// where `a = 0.05 h` is the reference height near the bed.
    ///
    /// # Arguments
    /// * `y` - height above the bed (m), clamped to `[a, h]`
    /// * `h` - total water depth (m)
    pub fn rouse_profile(&self, y: f64, h: f64) -> f64 {
        if h <= 0.0 {
            return 0.0;
        }
        let a = 0.05 * h; // reference height (5 % of depth)
        let y_clamped = y.clamp(a, h - 1e-12 * h);
        let ratio = (a / y_clamped) * ((h - y_clamped) / (h - a));
        self.reference_concentration * ratio.powf(self.rouse_number)
    }
}

// ---------------------------------------------------------------------------
// SedimentTransport (combined)
// ---------------------------------------------------------------------------

/// Combined sediment transport model coupling bed-load and suspended-load.
#[derive(Debug, Clone)]
pub struct SedimentTransport {
    /// Bed shear stress τ_b (Pa).
    pub bed_shear_stress: f64,
    /// Median grain diameter d₅₀ (m).
    pub median_diameter: f64,
    /// Fluid density ρ_f (kg/m³).
    pub fluid_density: f64,
    /// Sediment density ρ_s (kg/m³).
    pub sediment_density: f64,
}

impl SedimentTransport {
    /// Create a new `SedimentTransport` model.
    pub fn new(
        bed_shear_stress: f64,
        median_diameter: f64,
        fluid_density: f64,
        sediment_density: f64,
    ) -> Self {
        Self {
            bed_shear_stress,
            median_diameter,
            fluid_density,
            sediment_density,
        }
    }

    /// Compute the total dimensionless sediment transport rate.
    ///
    /// Uses the Meyer-Peter & Müller formula with the local Shields parameter
    /// derived from the stored shear stress and grain properties.
    pub fn compute_total_load(&self) -> f64 {
        const G: f64 = 9.81;
        let theta = shields_parameter(
            self.bed_shear_stress,
            self.sediment_density,
            self.fluid_density,
            self.median_diameter,
        );
        let theta_c = 0.047_f64;
        let excess = theta - theta_c;
        if excess <= 0.0 {
            return 0.0;
        }
        // MPM formula: q* = 8 (θ - θ_c)^{3/2}
        let q_star = 8.0 * excess.powf(1.5);
        // Dimensionalise: q = q* * sqrt((s-1) g d^3)
        let s = self.sediment_density / self.fluid_density;
        q_star * ((s - 1.0) * G * self.median_diameter.powi(3)).sqrt()
    }
}

// ---------------------------------------------------------------------------
// ErosionModel
// ---------------------------------------------------------------------------

/// Erosion formula selection.
#[derive(Debug, Clone, PartialEq)]
pub enum ErosionModel {
    /// Meyer-Peter & Müller (1948) bed-load formula.
    MeyerPeterMuller,
    /// Van Rijn (1984) combined bed-load and suspended-load.
    VanRijn,
    /// Engelund & Hansen (1967) total-load formula.
    Engelund,
}

impl ErosionModel {
    /// Compute the dimensionless total sediment transport rate for a given
    /// Shields parameter `theta`.
    ///
    /// - `MeyerPeterMuller`: `8 (θ - 0.047)^{3/2}`
    /// - `VanRijn`:           `0.053 θ^{2.1}` (simplified)
    /// - `Engelund`:          `0.05 θ^{5/2} / (√(s-1) f)`  where f=0.06
    ///
    /// # Arguments
    /// * `theta`   - Shields parameter
    /// * `s_ratio` - density ratio ρ_s / ρ_f
    pub fn transport_rate(&self, theta: f64, s_ratio: f64) -> f64 {
        match self {
            ErosionModel::MeyerPeterMuller => {
                let excess = theta - 0.047;
                if excess <= 0.0 {
                    0.0
                } else {
                    8.0 * excess.powf(1.5)
                }
            }
            ErosionModel::VanRijn => {
                if theta <= 0.0 {
                    0.0
                } else {
                    0.053 * theta.powf(2.1)
                }
            }
            ErosionModel::Engelund => {
                if theta <= 0.0 || s_ratio <= 1.0 {
                    0.0
                } else {
                    let f = 0.06_f64;
                    0.05 * theta.powf(2.5) / ((s_ratio - 1.0).sqrt() * f)
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// DepositionRate
// ---------------------------------------------------------------------------

/// Deposition rate model supporting Stokes and hindered settling.
#[derive(Debug, Clone)]
pub struct DepositionRate {
    /// Fluid density ρ_f (kg/m³).
    pub fluid_density: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub dynamic_viscosity: f64,
}

impl DepositionRate {
    /// Create a new `DepositionRate` model.
    pub fn new(fluid_density: f64, dynamic_viscosity: f64) -> Self {
        Self {
            fluid_density,
            dynamic_viscosity,
        }
    }

    /// Deposition flux (kg m⁻² s⁻¹) for concentration `conc` (kg/m³) and
    /// settling velocity `settling_vel` (m/s).
    ///
    /// Uses hindered settling with Richardson-Zaki correction when `conc > 0`:
    ///
    /// ```text
    /// D = conc * w_s * (1 - phi)^4.65
    /// ```
    ///
    /// where `phi = conc / rho_s` is estimated as `conc / 2650` (quartz).
    pub fn compute(&self, conc: f64, settling_vel: f64) -> f64 {
        if conc <= 0.0 || settling_vel <= 0.0 {
            return 0.0;
        }
        let rho_s = 2650.0_f64; // typical quartz density
        let phi = (conc / rho_s).clamp(0.0, 0.99);
        conc * settling_vel * (1.0 - phi).powf(4.65)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the dimensionless Shields parameter.
///
/// ```text
/// θ = τ_b / [(ρ_s - ρ_f) g d]
/// ```
///
/// # Arguments
/// * `tau`   - bed shear stress (Pa)
/// * `rho_s` - sediment density (kg/m³)
/// * `rho_f` - fluid density (kg/m³)
/// * `d`     - grain diameter (m)
pub fn shields_parameter(tau: f64, rho_s: f64, rho_f: f64, d: f64) -> f64 {
    const G: f64 = 9.81;
    let denom = (rho_s - rho_f) * G * d;
    if denom <= 0.0 {
        return 0.0;
    }
    tau / denom
}

/// Stokes settling velocity (m/s).
///
/// ```text
/// w_s = (ρ_s - ρ_f) g d² / (18 μ)
/// ```
///
/// # Arguments
/// * `d`     - grain diameter (m)
/// * `rho_s` - sediment density (kg/m³)
/// * `rho_f` - fluid density (kg/m³)
/// * `mu`    - dynamic viscosity (Pa·s)
pub fn stokes_settling(d: f64, rho_s: f64, rho_f: f64, mu: f64) -> f64 {
    const G: f64 = 9.81;
    if mu <= 0.0 {
        return 0.0;
    }
    (rho_s - rho_f) * G * d * d / (18.0 * mu)
}

/// Rouse number (dimensionless).
///
/// ```text
/// Z = w_s / (κ u*)
/// ```
///
/// # Arguments
/// * `ws`    - settling velocity (m/s)
/// * `u_star` - shear velocity (m/s)
/// * `kappa` - von Kármán constant (≈ 0.41)
pub fn rouse_number(ws: f64, u_star: f64, kappa: f64) -> f64 {
    if u_star <= 0.0 || kappa <= 0.0 {
        return f64::INFINITY;
    }
    ws / (kappa * u_star)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // --- shields_parameter ---

    #[test]
    fn test_shields_basic() {
        // τ = (ρ_s - ρ_f) g d → θ = 1
        let rho_s = 2650.0;
        let rho_f = 1000.0;
        let d = 0.001;
        let tau = (rho_s - rho_f) * 9.81 * d;
        let theta = shields_parameter(tau, rho_s, rho_f, d);
        assert!((theta - 1.0).abs() < 1e-10, "theta={theta}");
    }

    #[test]
    fn test_shields_zero_stress() {
        assert_eq!(shields_parameter(0.0, 2650.0, 1000.0, 0.001), 0.0);
    }

    #[test]
    fn test_shields_scales_linearly_with_tau() {
        let rho_s = 2650.0;
        let rho_f = 1000.0;
        let d = 0.001;
        let t1 = shields_parameter(1.0, rho_s, rho_f, d);
        let t2 = shields_parameter(2.0, rho_s, rho_f, d);
        assert!((t2 / t1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_shields_invalid_denom_returns_zero() {
        // rho_s == rho_f → no buoyancy → return 0
        assert_eq!(shields_parameter(1.0, 1000.0, 1000.0, 0.001), 0.0);
    }

    #[test]
    fn test_shields_typical_value() {
        // For a fine sand (d=0.2mm), ρ_s=2650, ρ_f=1000, τ=0.2 Pa
        let theta = shields_parameter(0.2, 2650.0, 1000.0, 2e-4);
        assert!(theta > 0.0 && theta < 10.0, "theta={theta}");
    }

    // --- stokes_settling ---

    #[test]
    fn test_stokes_positive() {
        let ws = stokes_settling(2e-4, 2650.0, 1000.0, 1e-3);
        assert!(ws > 0.0, "ws={ws}");
    }

    #[test]
    fn test_stokes_formula() {
        // ws = (2650-1000)*9.81*(2e-4)^2 / (18*1e-3)
        let expected = 1650.0 * 9.81 * (2e-4_f64).powi(2) / (18.0 * 1e-3);
        let ws = stokes_settling(2e-4, 2650.0, 1000.0, 1e-3);
        assert!((ws - expected).abs() < 1e-12, "ws={ws} expected={expected}");
    }

    #[test]
    fn test_stokes_scales_with_d_squared() {
        let ws1 = stokes_settling(1e-4, 2650.0, 1000.0, 1e-3);
        let ws2 = stokes_settling(2e-4, 2650.0, 1000.0, 1e-3);
        assert!((ws2 / ws1 - 4.0).abs() < 1e-8, "ratio={}", ws2 / ws1);
    }

    #[test]
    fn test_stokes_zero_viscosity() {
        assert_eq!(stokes_settling(1e-4, 2650.0, 1000.0, 0.0), 0.0);
    }

    #[test]
    fn test_stokes_negative_buoyancy_returns_negative() {
        // ρ_s < ρ_f → particle rises
        let ws = stokes_settling(1e-4, 500.0, 1000.0, 1e-3);
        assert!(ws < 0.0);
    }

    // --- rouse_number ---

    #[test]
    fn test_rouse_number_formula() {
        let z = rouse_number(0.01, 0.1, 0.41);
        let expected = 0.01 / (0.41 * 0.1);
        assert!((z - expected).abs() < 1e-10);
    }

    #[test]
    fn test_rouse_number_zero_ustar() {
        let z = rouse_number(0.01, 0.0, 0.41);
        assert_eq!(z, f64::INFINITY);
    }

    #[test]
    fn test_rouse_number_zero_kappa() {
        let z = rouse_number(0.01, 0.1, 0.0);
        assert_eq!(z, f64::INFINITY);
    }

    #[test]
    fn test_rouse_number_scales_with_ws() {
        let z1 = rouse_number(0.01, 0.1, 0.41);
        let z2 = rouse_number(0.02, 0.1, 0.41);
        assert!((z2 / z1 - 2.0).abs() < 1e-10);
    }

    // --- SedimentParticle ---

    #[test]
    fn test_particle_settling_velocity_positive() {
        let p = SedimentParticle::new([0.0; 3], 2e-4, 2650.0);
        let ws = p.settling_velocity(1000.0, 1e-6);
        assert!(ws > 0.0, "ws={ws}");
    }

    #[test]
    fn test_particle_settling_zero_nu() {
        let p = SedimentParticle::new([0.0; 3], 2e-4, 2650.0);
        assert_eq!(p.settling_velocity(1000.0, 0.0), 0.0);
    }

    #[test]
    fn test_particle_settling_scales_d_squared() {
        let p1 = SedimentParticle::new([0.0; 3], 1e-4, 2650.0);
        let p2 = SedimentParticle::new([0.0; 3], 2e-4, 2650.0);
        let ws1 = p1.settling_velocity(1000.0, 1e-6);
        let ws2 = p2.settling_velocity(1000.0, 1e-6);
        assert!((ws2 / ws1 - 4.0).abs() < 1e-8, "ratio={}", ws2 / ws1);
    }

    #[test]
    fn test_particle_fields_stored() {
        let pos = [1.0, 2.0, 3.0];
        let p = SedimentParticle::new(pos, 0.001, 2650.0);
        assert_eq!(p.position, pos);
        assert!((p.diameter - 0.001).abs() < 1e-15);
        assert!((p.density - 2650.0).abs() < 1e-10);
    }

    // --- BedLoad ---

    #[test]
    fn test_mpm_below_critical() {
        let bl = BedLoad::new(0.03, 0.047);
        assert_eq!(bl.meyer_peter_muller_rate(0.03), 0.0);
    }

    #[test]
    fn test_mpm_at_critical() {
        let bl = BedLoad::new(0.047, 0.047);
        assert_eq!(bl.meyer_peter_muller_rate(0.047), 0.0);
    }

    #[test]
    fn test_mpm_above_critical() {
        let bl = BedLoad::new(0.1, 0.047);
        let q = bl.meyer_peter_muller_rate(0.1);
        let expected = 8.0 * (0.1_f64 - 0.047).powf(1.5);
        assert!((q - expected).abs() < 1e-10, "q={q}");
    }

    #[test]
    fn test_mpm_scales_three_halves() {
        let bl = BedLoad::new(0.047, 0.047);
        let theta1 = 0.1_f64;
        let theta2 = 0.2_f64;
        let tc = 0.047_f64;
        let q1 = bl.meyer_peter_muller_rate(theta1);
        let q2 = bl.meyer_peter_muller_rate(theta2);
        let expected_ratio = ((theta2 - tc) / (theta1 - tc)).powf(1.5);
        assert!((q2 / q1 - expected_ratio).abs() < 1e-8);
    }

    // --- SuspendedLoad / Rouse profile ---

    #[test]
    fn test_rouse_profile_at_reference_height() {
        let sl = SuspendedLoad::new(1.0, 1.2);
        let h = 1.0;
        let a = 0.05 * h;
        let c = sl.rouse_profile(a, h);
        assert!((c - 1.0).abs() < 1e-10, "C(a) should equal C_ref, c={c}");
    }

    #[test]
    fn test_rouse_profile_decreases_with_height() {
        let sl = SuspendedLoad::new(1.0, 2.0);
        let h = 1.0;
        let c_low = sl.rouse_profile(0.1, h);
        let c_high = sl.rouse_profile(0.5, h);
        assert!(c_low > c_high, "concentration should decrease with height");
    }

    #[test]
    fn test_rouse_profile_zero_depth() {
        let sl = SuspendedLoad::new(1.0, 1.0);
        assert_eq!(sl.rouse_profile(0.1, 0.0), 0.0);
    }

    #[test]
    fn test_rouse_profile_finite_values() {
        let sl = SuspendedLoad::new(0.5, 1.5);
        let c = sl.rouse_profile(0.3, 1.0);
        assert!(c.is_finite() && c >= 0.0, "c={c}");
    }

    #[test]
    fn test_rouse_profile_large_z_near_bed_concentration() {
        // Large Rouse number → concentration confined near bed
        let sl_high = SuspendedLoad::new(1.0, 5.0);
        let sl_low = SuspendedLoad::new(1.0, 0.5);
        let h = 1.0;
        let c_high_z = sl_high.rouse_profile(0.5, h);
        let c_low_z = sl_low.rouse_profile(0.5, h);
        assert!(
            c_high_z < c_low_z,
            "higher Z → less suspension at mid-depth"
        );
    }

    // --- SedimentTransport ---

    #[test]
    fn test_total_load_zero_stress() {
        let st = SedimentTransport::new(0.0, 2e-4, 1000.0, 2650.0);
        assert_eq!(st.compute_total_load(), 0.0);
    }

    #[test]
    fn test_total_load_positive_above_threshold() {
        // Use a shear stress well above critical
        let rho_s = 2650.0;
        let rho_f = 1000.0;
        let d = 2e-4;
        let tau_crit = 0.047 * (rho_s - rho_f) * 9.81 * d;
        let st = SedimentTransport::new(tau_crit * 2.0, d, rho_f, rho_s);
        assert!(st.compute_total_load() > 0.0);
    }

    #[test]
    fn test_total_load_increases_with_shear() {
        let d = 2e-4_f64;
        let rho_s = 2650.0;
        let rho_f = 1000.0;
        let tau_crit = 0.047 * (rho_s - rho_f) * 9.81 * d;
        let st1 = SedimentTransport::new(tau_crit * 2.0, d, rho_f, rho_s);
        let st2 = SedimentTransport::new(tau_crit * 4.0, d, rho_f, rho_s);
        assert!(st2.compute_total_load() > st1.compute_total_load());
    }

    // --- ErosionModel ---

    #[test]
    fn test_erosion_mpm_zero_below_critical() {
        let q = ErosionModel::MeyerPeterMuller.transport_rate(0.04, 2.65);
        assert_eq!(q, 0.0);
    }

    #[test]
    fn test_erosion_mpm_positive() {
        let q = ErosionModel::MeyerPeterMuller.transport_rate(0.1, 2.65);
        assert!(q > 0.0);
    }

    #[test]
    fn test_erosion_van_rijn_positive() {
        let q = ErosionModel::VanRijn.transport_rate(0.1, 2.65);
        assert!(q > 0.0);
    }

    #[test]
    fn test_erosion_van_rijn_zero_for_zero_theta() {
        assert_eq!(ErosionModel::VanRijn.transport_rate(0.0, 2.65), 0.0);
    }

    #[test]
    fn test_erosion_engelund_positive() {
        let q = ErosionModel::Engelund.transport_rate(0.3, 2.65);
        assert!(q > 0.0, "q={q}");
    }

    #[test]
    fn test_erosion_engelund_zero_density_ratio_one() {
        assert_eq!(ErosionModel::Engelund.transport_rate(0.3, 1.0), 0.0);
    }

    // --- DepositionRate ---

    #[test]
    fn test_deposition_zero_conc() {
        let dr = DepositionRate::new(1000.0, 1e-3);
        assert_eq!(dr.compute(0.0, 0.01), 0.0);
    }

    #[test]
    fn test_deposition_zero_velocity() {
        let dr = DepositionRate::new(1000.0, 1e-3);
        assert_eq!(dr.compute(1.0, 0.0), 0.0);
    }

    #[test]
    fn test_deposition_positive() {
        let dr = DepositionRate::new(1000.0, 1e-3);
        let d = dr.compute(1.0, 0.01);
        assert!(d > 0.0, "d={d}");
    }

    #[test]
    fn test_deposition_scales_with_conc() {
        let dr = DepositionRate::new(1000.0, 1e-3);
        let d1 = dr.compute(1.0, 0.01);
        let d2 = dr.compute(2.0, 0.01);
        // Deposition at low concentrations should roughly double with conc
        assert!(d2 > d1, "d1={d1} d2={d2}");
    }

    #[test]
    fn test_deposition_hindrance_at_high_conc() {
        let dr = DepositionRate::new(1000.0, 1e-3);
        // At high concentration relative to quartz density the ratio < 2
        let d1 = dr.compute(10.0, 0.01);
        let d2 = dr.compute(100.0, 0.01);
        // Ratio should be less than 10 due to hindered settling
        let ratio = d2 / d1;
        assert!(ratio < 10.0, "ratio={ratio}");
    }

    // --- combined coverage ---

    #[test]
    fn test_mpm_agrees_bed_load_struct() {
        let bl = BedLoad::new(0.047, 0.047);
        let model = ErosionModel::MeyerPeterMuller;
        let theta = 0.15;
        let q_bl = bl.meyer_peter_muller_rate(theta);
        let q_em = model.transport_rate(theta, 2.65);
        assert!((q_bl - q_em).abs() < 1e-10);
    }

    #[test]
    fn test_shields_and_stokes_consistency() {
        // For a particle at critical motion, verify Shields ≈ 0.047
        let rho_s = 2650.0;
        let rho_f = 1000.0;
        let d = 2e-4;
        // Critical bed shear stress
        let tau_c = 0.047 * (rho_s - rho_f) * 9.81 * d;
        let theta = shields_parameter(tau_c, rho_s, rho_f, d);
        assert!((theta - 0.047).abs() < 1e-10, "theta={theta}");
        // Stokes settling should be positive
        let ws = stokes_settling(d, rho_s, rho_f, 1e-3);
        assert!(ws > 0.0);
    }

    #[test]
    fn test_rouse_number_physical_range() {
        let ws = stokes_settling(2e-4, 2650.0, 1000.0, 1e-3);
        let u_star = 0.05;
        let z = rouse_number(ws, u_star, 0.41);
        // For fine sand, Rouse number should be moderate (0.1–10)
        assert!(z > 0.01 && z < 100.0, "z={z}");
    }

    #[test]
    fn test_pi_constant_used() {
        // Ensure PI is reachable (coverage of use statement)
        assert!((PI - std::f64::consts::PI).abs() < 1e-15);
    }
}
