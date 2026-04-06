#![allow(clippy::manual_range_contains)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Microfluidic channel simulation tools.
//!
//! Covers Hagen-Poiseuille flow, slip flow (Knudsen regime), electroosmotic
//! flow, passive micromixers and droplet microfluidics.  All quantities are
//! in SI units unless stated otherwise.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Reynolds number for microfluidic channels:
/// Re = ρ · u · D_h / μ
pub fn reynolds_number_micro(rho: f64, u: f64, dh: f64, mu: f64) -> f64 {
    rho * u * dh / mu
}

/// Capillary number Ca = μ · u / γ.
pub fn capillary_number(mu: f64, u: f64, gamma: f64) -> f64 {
    mu * u / gamma
}

/// Knudsen number Kn = λ / D_h.
pub fn knudsen_number(lambda: f64, dh: f64) -> f64 {
    lambda / dh
}

// ---------------------------------------------------------------------------
// MicroChannel
// ---------------------------------------------------------------------------

/// Rectangular microfluidic channel geometry.
#[derive(Debug, Clone)]
pub struct MicroChannel {
    /// Channel width (m).
    pub width: f64,
    /// Channel height (m).
    pub height: f64,
    /// Channel length (m).
    pub length: f64,
}

impl MicroChannel {
    /// Create a new micro-channel.
    pub fn new(width: f64, height: f64, length: f64) -> Self {
        Self {
            width,
            height,
            length,
        }
    }

    /// Hydraulic diameter D_h = 4A / P for a rectangular cross-section.
    ///
    /// D_h = 2·w·h / (w + h)
    pub fn hydraulic_diameter(&self) -> f64 {
        4.0 * self.width * self.height / (2.0 * (self.width + self.height))
    }

    /// Aspect ratio α = h / w  (height / width, ≤ 1 by convention).
    pub fn aspect_ratio(&self) -> f64 {
        if self.width >= self.height {
            self.height / self.width
        } else {
            self.width / self.height
        }
    }
}

// ---------------------------------------------------------------------------
// HagenPoiseuilleFlow
// ---------------------------------------------------------------------------

/// Hagen-Poiseuille laminar flow in a circular capillary.
#[derive(Debug, Clone)]
pub struct HagenPoiseuilleFlow {
    /// Pressure drop ΔP along the capillary (Pa).
    pub pressure_drop: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Capillary radius R (m).
    pub radius: f64,
    /// Capillary length L (m).
    pub length: f64,
}

impl HagenPoiseuilleFlow {
    /// Create a new Hagen-Poiseuille flow configuration.
    pub fn new(pressure_drop: f64, viscosity: f64, radius: f64, length: f64) -> Self {
        Self {
            pressure_drop,
            viscosity,
            radius,
            length,
        }
    }

    /// Volumetric flow rate Q = π R⁴ ΔP / (8 μ L).
    pub fn flow_rate(&self) -> f64 {
        PI * self.radius.powi(4) * self.pressure_drop / (8.0 * self.viscosity * self.length)
    }

    /// Maximum (centreline) velocity u_max = R² ΔP / (4 μ L).
    pub fn max_velocity(&self) -> f64 {
        self.radius.powi(2) * self.pressure_drop / (4.0 * self.viscosity * self.length)
    }

    /// Parabolic velocity profile u(r) = u_max · (1 − (r/R)²).
    ///
    /// `r` is the radial position (m); returns 0 for r ≥ R.
    pub fn velocity_profile(&self, r: f64) -> f64 {
        if r >= self.radius {
            return 0.0;
        }
        self.max_velocity() * (1.0 - (r / self.radius).powi(2))
    }
}

// ---------------------------------------------------------------------------
// SlipFlow
// ---------------------------------------------------------------------------

/// Gas flow in a microchannel with velocity slip at the walls (Navier slip).
#[derive(Debug, Clone)]
pub struct SlipFlow {
    /// Knudsen number Kn = λ / D_h (dimensionless).
    pub knudsen_number: f64,
    /// Navier slip length b (m).
    pub slip_length: f64,
}

impl SlipFlow {
    /// Create a new slip flow configuration.
    pub fn new(knudsen_number: f64, slip_length: f64) -> Self {
        Self {
            knudsen_number,
            slip_length,
        }
    }

    /// Wall slip velocity u_slip = b · (∂u/∂n)|_wall = b · τ_w / μ.
    ///
    /// `wall_shear` is τ_w / μ, the wall shear rate (1/s).
    pub fn velocity_slip(&self, wall_shear: f64) -> f64 {
        self.slip_length * wall_shear
    }

    /// First-order Navier-slip correction factor for flow rate:
    /// C = 1 + 8 Kn  (from Maxwell's slip boundary condition).
    pub fn navier_slip_correction(&self) -> f64 {
        1.0 + 8.0 * self.knudsen_number
    }
}

// ---------------------------------------------------------------------------
// ElectroosmoticFlow
// ---------------------------------------------------------------------------

/// Electroosmotic (EOF) plug flow in a microchannel.
#[derive(Debug, Clone)]
pub struct ElectroosmoticFlow {
    /// Zeta potential ζ (V).
    pub zeta_potential: f64,
    /// Applied electric field E (V/m).
    pub electric_field: f64,
    /// Debye length λ_D (m).
    pub debye_length: f64,
    /// Permittivity of the medium ε = ε_r · ε_0 (F/m).
    pub permittivity: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
}

impl ElectroosmoticFlow {
    /// Create a new electroosmotic flow configuration.
    pub fn new(
        zeta_potential: f64,
        electric_field: f64,
        debye_length: f64,
        permittivity: f64,
        viscosity: f64,
    ) -> Self {
        Self {
            zeta_potential,
            electric_field,
            debye_length,
            permittivity,
            viscosity,
        }
    }

    /// Helmholtz-Smoluchowski EOF velocity:
    /// u_EOF = −ε ζ E / μ
    pub fn eof_velocity(&self) -> f64 {
        -self.permittivity * self.zeta_potential * self.electric_field / self.viscosity
    }

    /// Electroosmotic mobility μ_eo = −ε ζ / μ (m²/(V·s)).
    pub fn eof_mobility(&self) -> f64 {
        -self.permittivity * self.zeta_potential / self.viscosity
    }
}

// ---------------------------------------------------------------------------
// MicroMixer
// ---------------------------------------------------------------------------

/// Herringbone-groove passive micromixer.
#[derive(Debug, Clone)]
pub struct MicroMixer {
    /// Channel width w (m).
    pub channel_width: f64,
    /// Groove depth d (m).
    pub groove_depth: f64,
    /// Number of groove cycles along the channel.
    pub number_of_grooves: usize,
}

impl MicroMixer {
    /// Create a new micro-mixer geometry.
    pub fn new(channel_width: f64, groove_depth: f64, number_of_grooves: usize) -> Self {
        Self {
            channel_width,
            groove_depth,
            number_of_grooves,
        }
    }

    /// Empirical mixing efficiency η(Re) ∈ \[0, 1\].
    ///
    /// Uses the Stroock et al. model:
    /// η = 1 − exp(−α · N · (d/w)^0.5)
    /// where α = min(1, Re/100) captures the Re dependence, N is the
    /// number of grooves, d is groove depth, and w is channel width.
    pub fn mixing_efficiency(&self, re: f64) -> f64 {
        let alpha = (re / 100.0).min(1.0_f64);
        let dw = self.groove_depth / self.channel_width;
        1.0 - (-alpha * self.number_of_grooves as f64 * dw.sqrt()).exp()
    }

    /// Dean number De = Re · (w / 2R_c)^0.5 characterising secondary flow
    /// in a curved channel.
    ///
    /// `re` is the Reynolds number; `curvature_ratio` = w / R_c (channel
    /// width divided by centreline radius of curvature).
    pub fn dean_number(&self, re: f64, curvature_ratio: f64) -> f64 {
        re * (curvature_ratio / 2.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// DropletMicrofluidics
// ---------------------------------------------------------------------------

/// Droplet microfluidics — T-junction or flow-focusing droplet generation.
///
/// Note: the field `ca_number` stores the capillary number Ca = μ u / γ.
#[derive(Debug, Clone)]
pub struct DropletMicrofluidics {
    /// Droplet generation frequency f (Hz).
    pub generation_frequency: f64,
    /// Mean droplet diameter d (m).
    pub droplet_size: f64,
    /// Capillary number Ca = μ u / γ (dimensionless).
    pub ca_number: f64,
}

impl DropletMicrofluidics {
    /// Create a new droplet microfluidics configuration.
    pub fn new(generation_frequency: f64, droplet_size: f64, ca_number: f64) -> Self {
        Self {
            generation_frequency,
            droplet_size,
            ca_number,
        }
    }

    /// Droplet volume assuming spherical shape V = (4/3) π (d/2)³.
    pub fn volume(&self) -> f64 {
        let r = self.droplet_size / 2.0;
        (4.0 / 3.0) * PI * r.powi(3)
    }

    /// Volumetric throughput Q = f · V (m³/s).
    pub fn throughput(&self) -> f64 {
        self.generation_frequency * self.volume()
    }

    /// Squeezing-to-dripping transition criterion.
    ///
    /// Returns `true` when Ca < 0.1, indicating the squeezing (plug) regime.
    pub fn is_squeezing_regime(&self) -> bool {
        self.ca_number < 0.1
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── Free functions ───────────────────────────────────────────────────────

    #[test]
    fn test_reynolds_number_micro_basic() {
        // Re = 1000 * 0.01 * 100e-6 / 1e-3 = 1.0
        let re = reynolds_number_micro(1000.0, 0.01, 100e-6, 1e-3);
        assert!((re - 1.0).abs() < 1e-10, "Re = {re}");
    }

    #[test]
    fn test_reynolds_number_micro_water() {
        // Typical micro-channel: Re ≈ 1
        let re = reynolds_number_micro(997.0, 1e-3, 100e-6, 1e-3);
        assert!(re > 0.0 && re < 1.0);
    }

    #[test]
    fn test_capillary_number() {
        let ca = capillary_number(1e-3, 0.01, 0.072);
        let expected = 1e-3 * 0.01 / 0.072;
        assert!((ca - expected).abs() < 1e-15);
    }

    #[test]
    fn test_capillary_number_zero_velocity() {
        assert_eq!(capillary_number(1e-3, 0.0, 0.072), 0.0);
    }

    #[test]
    fn test_knudsen_number() {
        // λ = 68 nm, Dh = 1 µm → Kn = 0.068
        let kn = knudsen_number(68e-9, 1e-6);
        assert!((kn - 0.068).abs() < 1e-12);
    }

    #[test]
    fn test_knudsen_number_continuum() {
        // λ → 0 → Kn → 0
        assert_eq!(knudsen_number(0.0, 1e-6), 0.0);
    }

    // ── MicroChannel ─────────────────────────────────────────────────────────

    #[test]
    fn test_hydraulic_diameter_square() {
        let ch = MicroChannel::new(100e-6, 100e-6, 1e-3);
        // Square: Dh = w = h
        assert!((ch.hydraulic_diameter() - 100e-6).abs() < 1e-15);
    }

    #[test]
    fn test_hydraulic_diameter_rectangle() {
        let ch = MicroChannel::new(200e-6, 100e-6, 1e-3);
        let expected = 2.0 * 200e-6 * 100e-6 / (200e-6 + 100e-6);
        assert!((ch.hydraulic_diameter() - expected).abs() < 1e-20);
    }

    #[test]
    fn test_aspect_ratio_square() {
        let ch = MicroChannel::new(100e-6, 100e-6, 1e-3);
        assert!((ch.aspect_ratio() - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_aspect_ratio_tall() {
        let ch = MicroChannel::new(50e-6, 200e-6, 1e-3);
        // aspect_ratio = min(w,h)/max(w,h) = 0.25
        assert!((ch.aspect_ratio() - 0.25).abs() < 1e-15);
    }

    #[test]
    fn test_aspect_ratio_wide() {
        let ch = MicroChannel::new(400e-6, 100e-6, 1e-3);
        assert!((ch.aspect_ratio() - 0.25).abs() < 1e-15);
    }

    // ── HagenPoiseuilleFlow ──────────────────────────────────────────────────

    #[test]
    fn test_hp_max_velocity() {
        // R = 50 µm, ΔP = 1000 Pa, μ = 1e-3 Pa·s, L = 1 cm
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        let expected = (50e-6_f64).powi(2) * 1000.0 / (4.0 * 1e-3 * 1e-2);
        assert!((hp.max_velocity() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_hp_flow_rate() {
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        let expected = PI * (50e-6_f64).powi(4) * 1000.0 / (8.0 * 1e-3 * 1e-2);
        assert!((hp.flow_rate() - expected).abs() < 1e-25);
    }

    #[test]
    fn test_hp_velocity_profile_centre() {
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        // r = 0 → u(0) = u_max
        assert!((hp.velocity_profile(0.0) - hp.max_velocity()).abs() < 1e-20);
    }

    #[test]
    fn test_hp_velocity_profile_wall() {
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        // r = R → u = 0
        assert_eq!(hp.velocity_profile(50e-6), 0.0);
    }

    #[test]
    fn test_hp_velocity_profile_outside() {
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        assert_eq!(hp.velocity_profile(100e-6), 0.0);
    }

    #[test]
    fn test_hp_velocity_profile_midpoint() {
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        let r_half = 50e-6 / 2.0;
        let u_half = hp.velocity_profile(r_half);
        let u_max = hp.max_velocity();
        // u(R/2) = u_max * (1 - 1/4) = 0.75 u_max
        assert!((u_half - 0.75 * u_max).abs() < 1e-20);
    }

    #[test]
    fn test_hp_flow_rate_proportional_to_r4() {
        let hp1 = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        let hp2 = HagenPoiseuilleFlow::new(1000.0, 1e-3, 100e-6, 1e-2);
        // Doubling R → flow rate × 16
        let ratio = hp2.flow_rate() / hp1.flow_rate();
        assert!((ratio - 16.0).abs() < 1e-8);
    }

    // ── SlipFlow ─────────────────────────────────────────────────────────────

    #[test]
    fn test_slip_velocity() {
        let sf = SlipFlow::new(0.05, 1e-6);
        // u_slip = b * shear = 1e-6 * 1e6 = 1.0
        assert!((sf.velocity_slip(1e6) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_slip_velocity_zero_shear() {
        let sf = SlipFlow::new(0.05, 1e-6);
        assert_eq!(sf.velocity_slip(0.0), 0.0);
    }

    #[test]
    fn test_navier_slip_correction_continuum() {
        // Kn → 0 → C → 1
        let sf = SlipFlow::new(0.0, 0.0);
        assert_eq!(sf.navier_slip_correction(), 1.0);
    }

    #[test]
    fn test_navier_slip_correction_kn01() {
        let sf = SlipFlow::new(0.1, 1e-6);
        // C = 1 + 8 * 0.1 = 1.8
        assert!((sf.navier_slip_correction() - 1.8).abs() < 1e-14);
    }

    // ── ElectroosmoticFlow ────────────────────────────────────────────────────

    #[test]
    fn test_eof_velocity_sign() {
        // Negative zeta potential, positive E → positive EOF velocity
        let eof = ElectroosmoticFlow::new(-50e-3, 1e4, 10e-9, 7.1e-10, 1e-3);
        assert!(eof.eof_velocity() > 0.0);
    }

    #[test]
    fn test_eof_velocity_magnitude() {
        // ε = 7.1e-10, ζ = -50 mV, E = 10 kV/m, μ = 1e-3 Pa·s
        // u_EOF = -7.1e-10 * (-0.05) * 1e4 / 1e-3 = 3.55e-3 m/s
        let eof = ElectroosmoticFlow::new(-50e-3, 1e4, 10e-9, 7.1e-10, 1e-3);
        let expected = -7.1e-10 * (-50e-3) * 1e4 / 1e-3;
        assert!((eof.eof_velocity() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_eof_zero_zeta() {
        let eof = ElectroosmoticFlow::new(0.0, 1e4, 10e-9, 7.1e-10, 1e-3);
        assert_eq!(eof.eof_velocity(), 0.0);
    }

    #[test]
    fn test_eof_zero_field() {
        let eof = ElectroosmoticFlow::new(-50e-3, 0.0, 10e-9, 7.1e-10, 1e-3);
        assert_eq!(eof.eof_velocity(), 0.0);
    }

    #[test]
    fn test_eof_mobility() {
        let eof = ElectroosmoticFlow::new(-50e-3, 1e4, 10e-9, 7.1e-10, 1e-3);
        let expected = -7.1e-10 * (-50e-3) / 1e-3;
        assert!((eof.eof_mobility() - expected).abs() < 1e-20);
    }

    #[test]
    fn test_eof_velocity_proportional_to_field() {
        let eof1 = ElectroosmoticFlow::new(-50e-3, 1e4, 10e-9, 7.1e-10, 1e-3);
        let eof2 = ElectroosmoticFlow::new(-50e-3, 2e4, 10e-9, 7.1e-10, 1e-3);
        let ratio = eof2.eof_velocity() / eof1.eof_velocity();
        assert!((ratio - 2.0).abs() < 1e-10);
    }

    // ── MicroMixer ────────────────────────────────────────────────────────────

    #[test]
    fn test_mixing_efficiency_zero_grooves() {
        let mm = MicroMixer::new(100e-6, 50e-6, 0);
        // N = 0 → η = 1 − exp(0) = 0
        assert_eq!(mm.mixing_efficiency(10.0), 0.0);
    }

    #[test]
    fn test_mixing_efficiency_between_zero_and_one() {
        let mm = MicroMixer::new(100e-6, 20e-6, 10);
        let eta = mm.mixing_efficiency(50.0);
        assert!(eta >= 0.0 && eta <= 1.0, "η = {eta}");
    }

    #[test]
    fn test_mixing_efficiency_increases_with_grooves() {
        let mm1 = MicroMixer::new(100e-6, 20e-6, 5);
        let mm2 = MicroMixer::new(100e-6, 20e-6, 20);
        let re = 10.0;
        assert!(mm2.mixing_efficiency(re) > mm1.mixing_efficiency(re));
    }

    #[test]
    fn test_mixing_efficiency_alpha_clamp() {
        // Re > 100 → α clamped to 1
        let mm = MicroMixer::new(100e-6, 20e-6, 10);
        let eta_high = mm.mixing_efficiency(200.0);
        let eta_clamped = mm.mixing_efficiency(100.0);
        assert!((eta_high - eta_clamped).abs() < 1e-14);
    }

    #[test]
    fn test_dean_number_basic() {
        let mm = MicroMixer::new(100e-6, 20e-6, 10);
        // Re = 10, curvature = 0.1 → De = 10 * sqrt(0.05)
        let de = mm.dean_number(10.0, 0.1);
        let expected = 10.0 * (0.05_f64).sqrt();
        assert!((de - expected).abs() < 1e-12);
    }

    #[test]
    fn test_dean_number_zero_curvature() {
        let mm = MicroMixer::new(100e-6, 20e-6, 10);
        assert_eq!(mm.dean_number(10.0, 0.0), 0.0);
    }

    #[test]
    fn test_dean_number_proportional_to_re() {
        let mm = MicroMixer::new(100e-6, 20e-6, 10);
        let de1 = mm.dean_number(10.0, 0.1);
        let de2 = mm.dean_number(20.0, 0.1);
        assert!((de2 / de1 - 2.0).abs() < 1e-12);
    }

    // ── DropletMicrofluidics ──────────────────────────────────────────────────

    #[test]
    fn test_droplet_volume() {
        let dm = DropletMicrofluidics::new(100.0, 100e-6, 0.01);
        let expected = (4.0 / 3.0) * PI * (50e-6_f64).powi(3);
        assert!((dm.volume() - expected).abs() < 1e-30);
    }

    #[test]
    fn test_droplet_throughput() {
        let dm = DropletMicrofluidics::new(100.0, 100e-6, 0.01);
        let expected = 100.0 * (4.0 / 3.0) * PI * (50e-6_f64).powi(3);
        let rel_err = (dm.throughput() - expected).abs() / expected.abs();
        assert!(rel_err < 1e-12, "relative error = {rel_err}");
    }

    #[test]
    fn test_droplet_squeezing_regime_true() {
        let dm = DropletMicrofluidics::new(100.0, 100e-6, 0.05);
        assert!(dm.is_squeezing_regime());
    }

    #[test]
    fn test_droplet_squeezing_regime_false() {
        let dm = DropletMicrofluidics::new(100.0, 100e-6, 0.2);
        assert!(!dm.is_squeezing_regime());
    }

    #[test]
    fn test_droplet_volume_scales_as_r3() {
        let dm1 = DropletMicrofluidics::new(100.0, 100e-6, 0.01);
        let dm2 = DropletMicrofluidics::new(100.0, 200e-6, 0.01);
        let ratio = dm2.volume() / dm1.volume();
        assert!((ratio - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_hp_mean_velocity() {
        // Mean velocity = Q / (π R²) = u_max / 2
        let hp = HagenPoiseuilleFlow::new(1000.0, 1e-3, 50e-6, 1e-2);
        let q = hp.flow_rate();
        let a = PI * (50e-6_f64).powi(2);
        let u_mean = q / a;
        assert!((u_mean - hp.max_velocity() / 2.0).abs() < 1e-18);
    }

    #[test]
    fn test_micro_channel_new() {
        let ch = MicroChannel::new(50e-6, 30e-6, 5e-3);
        assert_eq!(ch.width, 50e-6);
        assert_eq!(ch.height, 30e-6);
        assert_eq!(ch.length, 5e-3);
    }

    #[test]
    fn test_slip_flow_new() {
        let sf = SlipFlow::new(0.01, 50e-9);
        assert_eq!(sf.knudsen_number, 0.01);
        assert_eq!(sf.slip_length, 50e-9);
    }

    #[test]
    fn test_eof_new_fields() {
        let eof = ElectroosmoticFlow::new(-25e-3, 5e3, 5e-9, 7e-10, 1e-3);
        assert_eq!(eof.zeta_potential, -25e-3);
        assert_eq!(eof.electric_field, 5e3);
    }

    #[test]
    fn test_capillary_number_values() {
        // Ca typical microfluidics: ~1e-5 to 1e-1
        let ca = capillary_number(1e-3, 0.001, 0.05);
        assert!(ca > 0.0);
        assert!(ca < 1.0);
    }

    #[test]
    fn test_mixing_efficiency_low_re() {
        let mm = MicroMixer::new(100e-6, 20e-6, 10);
        // At very low Re, alpha → 0 → efficiency → 0
        let eta = mm.mixing_efficiency(0.0);
        assert_eq!(eta, 0.0);
    }
}
