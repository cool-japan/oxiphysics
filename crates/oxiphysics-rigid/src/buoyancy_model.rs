// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Buoyancy and hydrostatics model for rigid body simulation.
//!
//! Implements Archimedes' principle, equilibrium draft, metacentric stability,
//! hull form parameters, Froude number, simplified hull resistance (Holtrop–Mennen),
//! shallow-water squat, water-entry slamming, and added mass for a sphere.
//!
//! # Overview
//!
//! - [`FluidEnvironment`] — fluid properties (density, gravity, viscosity, surface height).
//! - [`SubmergedBody`] — geometric and hydrostatic properties of a floating body.
//! - [`HullForm`] — ship hull form parameters (L, B, T, Cb).
//! - [`WaterEntry`] — water-entry (slamming) model based on von Kármán theory.
//! - [`archimedes_principle`] — `F_b = ρ·V·g`.
//! - [`equilibrium_draft`] — draft for static equilibrium.
//! - [`stability_criterion`] — textual stability classification.
//! - [`froude_number`] — dimensionless speed.
//! - [`hull_resistance_estimate`] — simplified total resistance.
//! - [`sinkage_squat`] — shallow-water squat formula.
//! - [`added_mass_sphere`] — fluid added mass of a sphere.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// FluidEnvironment
// ─────────────────────────────────────────────────────────────────────────────

/// Properties of a fluid environment (liquid or gas).
#[derive(Debug, Clone)]
pub struct FluidEnvironment {
    /// Fluid density \[kg/m³\].
    pub density: f64,
    /// Gravitational acceleration \[m/s²\] (positive downward).
    pub gravity: f64,
    /// Height of the free surface in the world frame \[m\].
    pub surface_height: f64,
    /// Dynamic viscosity \[Pa·s\].
    pub viscosity: f64,
}

impl FluidEnvironment {
    /// Create a new fluid environment with given density and gravity.
    ///
    /// `surface_height` defaults to `0.0` and `viscosity` to `1.0e-3 Pa·s`
    /// (water at 20 °C).
    pub fn new(density: f64, gravity: f64) -> Self {
        Self {
            density,
            gravity,
            surface_height: 0.0,
            viscosity: 1.0e-3,
        }
    }

    /// Hydrostatic pressure at a given depth below the free surface \[Pa\].
    ///
    /// `depth` is the positive distance below the surface.
    pub fn hydrostatic_pressure(&self, depth: f64) -> f64 {
        self.density * self.gravity * depth
    }

    /// Standard atmosphere at sea level (ρ ≈ 1.225 kg/m³, g = 9.81 m/s²).
    pub fn atmospheric() -> Self {
        Self {
            density: 1.225,
            gravity: 9.81,
            surface_height: 0.0,
            viscosity: 1.81e-5,
        }
    }

    /// Fresh water at 20 °C (ρ = 998 kg/m³, g = 9.81 m/s²).
    pub fn water() -> Self {
        Self::new(998.0, 9.81)
    }

    /// Salt water at 15 °C (ρ = 1025 kg/m³, g = 9.81 m/s²).
    pub fn seawater() -> Self {
        Self {
            density: 1025.0,
            gravity: 9.81,
            surface_height: 0.0,
            viscosity: 1.07e-3,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// SubmergedBody
// ─────────────────────────────────────────────────────────────────────────────

/// Geometric and hydrostatic description of a (partially) submerged body.
#[derive(Debug, Clone)]
pub struct SubmergedBody {
    /// Displaced volume of fluid \[m³\].
    pub volume: f64,
    /// Draft (depth of the keel below the waterline) \[m\].
    pub draft: f64,
    /// Area of the waterplane (intersection with free surface) \[m²\].
    pub waterplane_area: f64,
    /// Centre of buoyancy in body frame `[x, y, z]` \[m\].
    pub center_of_buoyancy: [f64; 3],
    /// Centre of gravity in body frame `[x, y, z]` \[m\].
    pub center_of_gravity: [f64; 3],
}

impl SubmergedBody {
    /// Buoyant force magnitude \[N\].
    ///
    /// `F_b = ρ · V · g`
    pub fn buoyant_force(&self, fluid: &FluidEnvironment) -> f64 {
        archimedes_principle(self.volume, fluid.density, fluid.gravity)
    }

    /// Net vertical force on the body \[N\].
    ///
    /// Positive means upward resultant (body tends to rise).
    /// `F_net = F_b − m · g`
    pub fn net_vertical_force(&self, mass: f64, fluid: &FluidEnvironment) -> f64 {
        self.buoyant_force(fluid) - mass * fluid.gravity
    }

    /// Metacentric height GM \[m\].
    ///
    /// `GM = BM − BG`
    ///
    /// where `BM = I_wl / V` (second moment of waterplane area over displaced volume)
    /// and `BG` is the signed vertical distance from B to G.
    ///
    /// `second_moment_waterplane` — second moment of the waterplane area about
    /// the longitudinal axis \[m⁴\].
    pub fn metacentric_height(&self, second_moment_waterplane: f64) -> f64 {
        let bm = second_moment_waterplane / self.volume.max(1e-14);
        // BG = z_G − z_B  (positive when G is above B)
        let bg = self.center_of_gravity[2] - self.center_of_buoyancy[2];
        bm - bg
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free functions — basic hydrostatics
// ─────────────────────────────────────────────────────────────────────────────

/// Buoyant force according to Archimedes' principle \[N\].
///
/// `F_b = ρ · V · g`
///
/// * `volume` — displaced fluid volume \[m³\]
/// * `fluid_density` — fluid density \[kg/m³\]
/// * `g` — gravitational acceleration \[m/s²\]
pub fn archimedes_principle(volume: f64, fluid_density: f64, g: f64) -> f64 {
    fluid_density * volume * g
}

/// Equilibrium draft of a freely floating body \[m\].
///
/// Derived from `ρ · A_wp · T = m`  →  `T = m / (ρ · A_wp)`.
///
/// * `mass` — mass of the body \[kg\]
/// * `waterplane_area` — waterplane area \[m²\]
/// * `fluid_density` — fluid density \[kg/m³\]
pub fn equilibrium_draft(mass: f64, waterplane_area: f64, fluid_density: f64) -> f64 {
    mass / (fluid_density * waterplane_area.max(1e-14))
}

/// Classify hydrostatic stability based on metacentric height GM \[m\].
///
/// Returns:
/// * `"stable"`  — GM > 0
/// * `"neutral"` — GM ≈ 0 (|GM| < 1e-6)
/// * `"unstable"` — GM < 0
pub fn stability_criterion(gm: f64) -> &'static str {
    if gm > 1e-6 {
        "stable"
    } else if gm < -1e-6 {
        "unstable"
    } else {
        "neutral"
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// HullForm
// ─────────────────────────────────────────────────────────────────────────────

/// Ship hull form parameters.
#[derive(Debug, Clone)]
pub struct HullForm {
    /// Length between perpendiculars \[m\].
    pub length: f64,
    /// Beam (breadth moulded) \[m\].
    pub beam: f64,
    /// Design draft \[m\].
    pub draft: f64,
    /// Block coefficient `Cb = ∇ / (L · B · T)` (dimensionless, ∈ (0, 1]).
    pub block_coefficient: f64,
}

impl HullForm {
    /// Displaced volume \[m³\].
    ///
    /// `∇ = Cb · L · B · T`
    pub fn displaced_volume(&self) -> f64 {
        self.block_coefficient * self.length * self.beam * self.draft
    }

    /// Second moment of the waterplane area about the longitudinal axis \[m⁴\].
    ///
    /// Approximated as `I_L ≈ (1/12) · L · B³` (rectangular waterplane).
    pub fn waterplane_inertia(&self) -> f64 {
        self.length * self.beam.powi(3) / 12.0
    }

    /// Prismatic coefficient `Cp`.
    ///
    /// `Cp = Cb / Cm` where `Cm` (midship section coefficient) is approximated
    /// as `0.98` for typical commercial hulls.
    pub fn prismatic_coefficient(&self) -> f64 {
        let cm = 0.98_f64;
        self.block_coefficient / cm
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Free functions — ship resistance & seakeeping
// ─────────────────────────────────────────────────────────────────────────────

/// Froude number (dimensionless speed) \[-\].
///
/// `Fr = v / √(g · L)`
///
/// * `speed` — ship speed \[m/s\]
/// * `length` — ship length \[m\]
pub fn froude_number(speed: f64, length: f64) -> f64 {
    let g = 9.81_f64;
    speed / (g * length.max(1e-14)).sqrt()
}

/// Simplified total hull resistance \[N\] (Holtrop–Mennen spirit).
///
/// Uses a simplified formulation:
/// `R_T = 0.5 · ρ · v² · S_wet · Cf`
///
/// where the wetted surface is estimated as `S ≈ L · (2T + B)` and
/// the frictional coefficient `Cf = 0.075 / (log10(Rn) − 2)²` (ITTC 1957).
///
/// A wave-making resistance term `R_w = 0.03 · Fr² · Δ · g` is added,
/// where `Δ = ρ · ∇` is the displacement in \[kg\].
pub fn hull_resistance_estimate(hull: &HullForm, speed: f64, fluid: &FluidEnvironment) -> f64 {
    if speed < 1e-9 {
        return 0.0;
    }
    // Wetted surface (approximation)
    let s_wet = hull.length * (2.0 * hull.draft + hull.beam);
    // Reynolds number
    let rn = fluid.density * speed * hull.length / fluid.viscosity.max(1e-14);
    let log_rn = rn.log10().max(2.001);
    let cf = 0.075 / (log_rn - 2.0).powi(2);
    // Frictional resistance
    let rf = 0.5 * fluid.density * speed * speed * s_wet * cf;
    // Wave-making resistance (Froude-based, very simplified)
    let fr = froude_number(speed, hull.length);
    let displacement_kg = fluid.density * hull.displaced_volume();
    let rw = 0.03 * fr * fr * displacement_kg * fluid.gravity;
    rf + rw
}

/// Shallow-water squat (Barras simplified formula) \[m\].
///
/// `S_q = Cb · v² / (20 · √h)` where `h` is water depth \[m\].
///
/// * `speed` — ship speed \[m/s\]
/// * `length` — ship length (used as scale reference, not used in Barras basic) \[m\]
/// * `block_coeff` — block coefficient Cb \[-\]
/// * `water_depth` — water depth \[m\]
pub fn sinkage_squat(speed: f64, _length: f64, block_coeff: f64, water_depth: f64) -> f64 {
    block_coeff * speed * speed / (20.0 * water_depth.max(1e-3).sqrt())
}

// ─────────────────────────────────────────────────────────────────────────────
// WaterEntry
// ─────────────────────────────────────────────────────────────────────────────

/// Water-entry (slamming) body parameters.
#[derive(Debug, Clone)]
pub struct WaterEntry {
    /// Body mass \[kg\].
    pub mass: f64,
    /// Entry velocity (positive downward) \[m/s\].
    pub velocity: f64,
    /// Deadrise angle \[deg\] (0° = flat plate, 90° = knife edge).
    pub deadrise_angle: f64,
}

impl WaterEntry {
    /// Von Kármán impact force \[N\] at time `t` after first contact.
    ///
    /// Approximated as: `F = ρ · π · (v · t · cot β)² · v`
    ///
    /// where `β` is the deadrise angle (clamped away from 0°).
    pub fn von_karman_impact_force(&self, fluid_density: f64, t: f64) -> f64 {
        let beta_deg = self.deadrise_angle.max(1.0);
        let beta_rad = beta_deg.to_radians();
        let cot_beta = (beta_rad.cos() / beta_rad.sin()).abs();
        // expanding wetted half-beam: c(t) = v * t * cot(β)
        let c = self.velocity * t * cot_beta;
        fluid_density * PI * c * c * self.velocity
    }

    /// Peak impact pressure \[Pa\] using the simplified Wagner formula.
    ///
    /// `p_peak = 0.5 · ρ · v² · (π / (2 · tan β))²`
    pub fn peak_pressure(&self, fluid_density: f64) -> f64 {
        let beta_deg = self.deadrise_angle.max(1.0);
        let beta_rad = beta_deg.to_radians();
        let tan_beta = beta_rad.tan().abs().max(1e-6);
        let factor = PI / (2.0 * tan_beta);
        0.5 * fluid_density * self.velocity * self.velocity * factor * factor
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Added mass
// ─────────────────────────────────────────────────────────────────────────────

/// Fluid added mass of a sphere \[kg\].
///
/// `m_added = 0.5 · ρ · V_sphere` where `V = (4/3)·π·r³`.
///
/// * `radius` — sphere radius \[m\]
/// * `fluid_density` — fluid density \[kg/m³\]
pub fn added_mass_sphere(radius: f64, fluid_density: f64) -> f64 {
    let volume = (4.0 / 3.0) * PI * radius.powi(3);
    0.5 * fluid_density * volume
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // ── archimedes_principle ─────────────────────────────────────────────

    #[test]
    fn archimedes_positive_volume() {
        let fb = archimedes_principle(1.0, 1000.0, 9.81);
        assert!((fb - 9810.0).abs() < 1e-6);
    }

    #[test]
    fn archimedes_zero_volume_is_zero() {
        assert_eq!(archimedes_principle(0.0, 1000.0, 9.81), 0.0);
    }

    #[test]
    fn archimedes_scales_linearly_with_volume() {
        let fb1 = archimedes_principle(1.0, 1000.0, 9.81);
        let fb2 = archimedes_principle(2.0, 1000.0, 9.81);
        assert!((fb2 - 2.0 * fb1).abs() < EPS);
    }

    #[test]
    fn archimedes_scales_linearly_with_density() {
        let fb1 = archimedes_principle(1.0, 1000.0, 9.81);
        let fb2 = archimedes_principle(1.0, 2000.0, 9.81);
        assert!((fb2 - 2.0 * fb1).abs() < EPS);
    }

    #[test]
    fn archimedes_scales_with_gravity() {
        let fb1 = archimedes_principle(1.0, 1000.0, 9.81);
        let fb2 = archimedes_principle(1.0, 1000.0, 19.62);
        assert!((fb2 - 2.0 * fb1).abs() < EPS);
    }

    // ── equilibrium_draft ────────────────────────────────────────────────

    #[test]
    fn equilibrium_draft_basic() {
        // T = m / (rho * A_wp) = 1000 / (1000 * 1.0) = 1.0 m
        let t = equilibrium_draft(1000.0, 1.0, 1000.0);
        assert!((t - 1.0).abs() < EPS);
    }

    #[test]
    fn equilibrium_draft_positive() {
        let t = equilibrium_draft(500.0, 2.0, 998.0);
        assert!(t > 0.0);
    }

    #[test]
    fn equilibrium_draft_scales_with_mass() {
        let t1 = equilibrium_draft(1000.0, 2.0, 1000.0);
        let t2 = equilibrium_draft(2000.0, 2.0, 1000.0);
        assert!((t2 - 2.0 * t1).abs() < EPS);
    }

    #[test]
    fn equilibrium_draft_large_waterplane_gives_shallow_draft() {
        let t = equilibrium_draft(1000.0, 1000.0, 1000.0);
        assert!(t < 0.01);
    }

    // ── stability_criterion ──────────────────────────────────────────────

    #[test]
    fn stability_positive_gm_is_stable() {
        assert_eq!(stability_criterion(0.5), "stable");
    }

    #[test]
    fn stability_negative_gm_is_unstable() {
        assert_eq!(stability_criterion(-0.5), "unstable");
    }

    #[test]
    fn stability_zero_gm_is_neutral() {
        assert_eq!(stability_criterion(0.0), "neutral");
    }

    #[test]
    fn stability_tiny_positive_is_stable() {
        assert_eq!(stability_criterion(1e-4), "stable");
    }

    #[test]
    fn stability_tiny_negative_is_unstable() {
        assert_eq!(stability_criterion(-1e-4), "unstable");
    }

    // ── FluidEnvironment ─────────────────────────────────────────────────

    #[test]
    fn water_density_approx_998() {
        let fluid = FluidEnvironment::water();
        assert!((fluid.density - 998.0).abs() < 1.0);
    }

    #[test]
    fn atmospheric_density_approx_1p225() {
        let fluid = FluidEnvironment::atmospheric();
        assert!((fluid.density - 1.225).abs() < 0.01);
    }

    #[test]
    fn hydrostatic_pressure_at_10m() {
        let fluid = FluidEnvironment::water();
        let p = fluid.hydrostatic_pressure(10.0);
        // approx 998 * 9.81 * 10 = 97938 Pa
        assert!(p > 90_000.0 && p < 110_000.0);
    }

    #[test]
    fn hydrostatic_pressure_zero_at_surface() {
        let fluid = FluidEnvironment::water();
        assert!((fluid.hydrostatic_pressure(0.0)).abs() < EPS);
    }

    #[test]
    fn seawater_denser_than_fresh() {
        let sw = FluidEnvironment::seawater();
        let fw = FluidEnvironment::water();
        assert!(sw.density > fw.density);
    }

    // ── SubmergedBody ────────────────────────────────────────────────────

    #[test]
    fn buoyant_force_greater_than_weight_body_floats() {
        let body = SubmergedBody {
            volume: 10.0,
            draft: 0.5,
            waterplane_area: 5.0,
            center_of_buoyancy: [0.0, 0.0, -0.25],
            center_of_gravity: [0.0, 0.0, 0.0],
        };
        let fluid = FluidEnvironment::water();
        let mass = 5000.0; // lighter than displaced water (998 * 10 = 9980 kg)
        assert!(body.net_vertical_force(mass, &fluid) > 0.0);
    }

    #[test]
    fn net_vertical_force_negative_means_sinks() {
        let body = SubmergedBody {
            volume: 1.0,
            draft: 1.0,
            waterplane_area: 1.0,
            center_of_buoyancy: [0.0, 0.0, -0.5],
            center_of_gravity: [0.0, 0.0, 0.0],
        };
        let fluid = FluidEnvironment::water();
        let mass = 5000.0; // heavier than displaced water (~998 kg)
        assert!(body.net_vertical_force(mass, &fluid) < 0.0);
    }

    #[test]
    fn buoyant_force_equals_archimedes() {
        let body = SubmergedBody {
            volume: 2.0,
            draft: 1.0,
            waterplane_area: 2.0,
            center_of_buoyancy: [0.0, 0.0, -0.5],
            center_of_gravity: [0.0, 0.0, 0.0],
        };
        let fluid = FluidEnvironment::water();
        let expected = archimedes_principle(2.0, fluid.density, fluid.gravity);
        assert!((body.buoyant_force(&fluid) - expected).abs() < EPS);
    }

    #[test]
    fn metacentric_height_positive_is_stable() {
        let body = SubmergedBody {
            volume: 10.0,
            draft: 1.0,
            waterplane_area: 5.0,
            // B below G → BG positive
            center_of_buoyancy: [0.0, 0.0, -0.5],
            center_of_gravity: [0.0, 0.0, 0.2],
        };
        // I_wl large → BM > BG
        let gm = body.metacentric_height(100.0);
        assert!(gm > 0.0);
    }

    // ── HullForm ─────────────────────────────────────────────────────────

    #[test]
    fn block_coefficient_in_valid_range() {
        let hull = HullForm {
            length: 100.0,
            beam: 20.0,
            draft: 5.0,
            block_coefficient: 0.75,
        };
        assert!(hull.block_coefficient > 0.0 && hull.block_coefficient <= 1.0);
    }

    #[test]
    fn displaced_volume_formula() {
        let hull = HullForm {
            length: 100.0,
            beam: 20.0,
            draft: 5.0,
            block_coefficient: 0.75,
        };
        let expected = 0.75 * 100.0 * 20.0 * 5.0;
        assert!((hull.displaced_volume() - expected).abs() < EPS);
    }

    #[test]
    fn waterplane_inertia_positive() {
        let hull = HullForm {
            length: 80.0,
            beam: 15.0,
            draft: 4.0,
            block_coefficient: 0.7,
        };
        assert!(hull.waterplane_inertia() > 0.0);
    }

    #[test]
    fn prismatic_coefficient_reasonable() {
        let hull = HullForm {
            length: 100.0,
            beam: 20.0,
            draft: 5.0,
            block_coefficient: 0.7,
        };
        let cp = hull.prismatic_coefficient();
        assert!(cp > 0.5 && cp < 1.0);
    }

    // ── froude_number ────────────────────────────────────────────────────

    #[test]
    fn froude_number_zero_at_rest() {
        assert_eq!(froude_number(0.0, 100.0), 0.0);
    }

    #[test]
    fn froude_number_positive() {
        assert!(froude_number(5.0, 100.0) > 0.0);
    }

    #[test]
    fn froude_number_increases_with_speed() {
        let fr1 = froude_number(5.0, 100.0);
        let fr2 = froude_number(10.0, 100.0);
        assert!(fr2 > fr1);
    }

    #[test]
    fn froude_number_unit_gravity_check() {
        // Fr = v / sqrt(g * L) = 1 / sqrt(9.81 * 9.81) = 1/9.81
        let fr = froude_number(1.0, 9.81);
        let expected = 1.0 / (9.81_f64 * 9.81_f64).sqrt();
        assert!((fr - expected).abs() < 1e-9);
    }

    // ── hull_resistance_estimate ─────────────────────────────────────────

    #[test]
    fn hull_resistance_zero_at_rest() {
        let hull = HullForm {
            length: 100.0,
            beam: 20.0,
            draft: 5.0,
            block_coefficient: 0.75,
        };
        let fluid = FluidEnvironment::water();
        assert_eq!(hull_resistance_estimate(&hull, 0.0, &fluid), 0.0);
    }

    #[test]
    fn hull_resistance_positive_at_speed() {
        let hull = HullForm {
            length: 100.0,
            beam: 20.0,
            draft: 5.0,
            block_coefficient: 0.75,
        };
        let fluid = FluidEnvironment::water();
        assert!(hull_resistance_estimate(&hull, 5.0, &fluid) > 0.0);
    }

    // ── sinkage_squat ────────────────────────────────────────────────────

    #[test]
    fn sinkage_squat_zero_at_rest() {
        assert_eq!(sinkage_squat(0.0, 100.0, 0.75, 15.0), 0.0);
    }

    #[test]
    fn sinkage_squat_positive_at_speed() {
        let sq = sinkage_squat(5.0, 100.0, 0.75, 15.0);
        assert!(sq > 0.0);
    }

    #[test]
    fn sinkage_squat_increases_with_speed() {
        let sq1 = sinkage_squat(3.0, 100.0, 0.75, 15.0);
        let sq2 = sinkage_squat(6.0, 100.0, 0.75, 15.0);
        assert!(sq2 > sq1);
    }

    // ── added_mass_sphere ────────────────────────────────────────────────

    #[test]
    fn added_mass_sphere_formula() {
        let r = 1.0_f64;
        let rho = 1000.0_f64;
        let expected = 0.5 * rho * (4.0 / 3.0) * PI * r.powi(3);
        assert!((added_mass_sphere(r, rho) - expected).abs() < 1e-9);
    }

    #[test]
    fn added_mass_sphere_positive() {
        assert!(added_mass_sphere(0.5, 1000.0) > 0.0);
    }

    #[test]
    fn added_mass_sphere_scales_with_cube_of_radius() {
        let m1 = added_mass_sphere(1.0, 1000.0);
        let m2 = added_mass_sphere(2.0, 1000.0);
        assert!((m2 - 8.0 * m1).abs() < 1e-6);
    }

    // ── WaterEntry ───────────────────────────────────────────────────────

    #[test]
    fn peak_pressure_positive() {
        let we = WaterEntry {
            mass: 10.0,
            velocity: 5.0,
            deadrise_angle: 10.0,
        };
        assert!(we.peak_pressure(1000.0) > 0.0);
    }

    #[test]
    fn von_karman_force_positive_at_nonzero_t() {
        let we = WaterEntry {
            mass: 10.0,
            velocity: 5.0,
            deadrise_angle: 15.0,
        };
        assert!(we.von_karman_impact_force(1000.0, 0.01) > 0.0);
    }

    #[test]
    fn von_karman_force_zero_at_t_zero() {
        let we = WaterEntry {
            mass: 10.0,
            velocity: 5.0,
            deadrise_angle: 15.0,
        };
        assert!(we.von_karman_impact_force(1000.0, 0.0).abs() < EPS);
    }
}
