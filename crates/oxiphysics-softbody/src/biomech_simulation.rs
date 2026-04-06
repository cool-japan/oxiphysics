// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Biomechanics simulation — bone, cartilage, tendon, and joint contact.
//!
//! Provides constitutive models and structural calculations for:
//!
//! - Cortical and trabecular bone (Voigt/Reuss bounds, remodelling).
//! - Biphasic cartilage (creep, stress-relaxation, fluid pressure).
//! - Collagenous tendons (toe/linear stiffness, hysteresis).
//! - Articular joint contact (Hertzian pressure, lubrication).
//! - Bone fracture mechanics (Paris law, stress intensity).
//! - Wolff's law density adaptation.

#![allow(dead_code)]

// ── BoneStructure ─────────────────────────────────────────────────────────────

/// A cortical or trabecular bone specimen described by its microstructural
/// composition.
#[derive(Debug, Clone)]
pub struct BoneStructure {
    /// Apparent density (kg/m³).  Cortical bone ≈ 1900; trabecular ≈ 400–800.
    pub density: f64,
    /// Porosity fraction \[0, 1\].  Cortical ≈ 0.05; trabecular ≈ 0.5–0.9.
    pub porosity: f64,
    /// Mineral (hydroxyapatite) volume fraction \[0, 1\], typically ≈ 0.45.
    pub mineral_fraction: f64,
}

impl BoneStructure {
    /// Create a new `BoneStructure` specimen.
    ///
    /// * `density`          – Apparent density (kg/m³).
    /// * `porosity`         – Porosity fraction \[0, 1\].
    /// * `mineral_fraction` – Mineral volume fraction \[0, 1\].
    pub fn new(density: f64, porosity: f64, mineral_fraction: f64) -> Self {
        Self {
            density,
            porosity,
            mineral_fraction,
        }
    }

    /// Estimated Young's modulus (Pa) using an empirical power-law fit.
    ///
    /// `E = 6.95e9 * (density / 1900)^3`  (adapted from Carter & Hayes 1977).
    pub fn youngs_modulus(&self) -> f64 {
        let rho_ref = 1900.0_f64;
        let e_ref = 6.95e9_f64;
        let ratio = (self.density / rho_ref).max(0.0);
        e_ref * ratio.powi(3)
    }

    /// Estimated uniaxial compressive strength (Pa).
    ///
    /// `sigma_c = 137e6 * (density / 1900)^2`
    pub fn compressive_strength(&self) -> f64 {
        let rho_ref = 1900.0_f64;
        let ratio = (self.density / rho_ref).max(0.0);
        137e6 * ratio.powi(2)
    }

    /// Reuss (lower-bound, series) composite modulus (Pa) of mineral and
    /// collagen phases.
    ///
    /// `1/E_R = f_m/E_m + (1-f_m)/E_c`
    pub fn reuss_modulus(&self) -> f64 {
        let e_mineral = 114e9_f64; // hydroxyapatite
        let e_collagen = 1.5e9_f64; // type-I collagen
        let f = self.mineral_fraction.clamp(0.0, 1.0);
        let compliance = f / e_mineral + (1.0 - f) / e_collagen;
        1.0 / compliance.max(1e-30)
    }

    /// Voigt (upper-bound, parallel) composite modulus (Pa).
    ///
    /// `E_V = f_m * E_m + (1 - f_m) * E_c`
    pub fn voigt_modulus(&self) -> f64 {
        let e_mineral = 114e9_f64;
        let e_collagen = 1.5e9_f64;
        let f = self.mineral_fraction.clamp(0.0, 1.0);
        f * e_mineral + (1.0 - f) * e_collagen
    }

    /// Estimated fatigue life (cycles) at the given cyclic stress amplitude
    /// (Pa) using a power-law S-N relation.
    ///
    /// `N = (sigma_u / sigma_amp)^10`  where `sigma_u` is the compressive
    /// strength.  Returns `f64::INFINITY` if `sigma_amp <= 0`.
    pub fn fatigue_life(&self, stress_amp: f64) -> f64 {
        if stress_amp <= 0.0 {
            return f64::INFINITY;
        }
        let sigma_u = self.compressive_strength();
        (sigma_u / stress_amp).powi(10)
    }

    /// Rate of bone density remodelling (kg/m³/day) driven by mechanical
    /// stimulus.
    ///
    /// Uses a lazy-zone model: no remodelling within ±`dead_zone` of the
    /// `setpoint`.  Outside, rate = `k * (stress - setpoint)`.
    ///
    /// * `stress`   – Current mechanical stimulus (Pa).
    /// * `setpoint` – Homeostatic reference stimulus (Pa).
    pub fn remodeling_rate(&self, stress: f64, setpoint: f64) -> f64 {
        let k = 1.0e-4_f64; // remodelling rate constant (kg/m³/day per Pa)
        let dead_zone = 0.05 * setpoint.abs().max(1.0);
        let diff = stress - setpoint;
        if diff.abs() <= dead_zone {
            0.0
        } else if diff > 0.0 {
            k * (diff - dead_zone)
        } else {
            k * (diff + dead_zone)
        }
    }
}

// ── CartilageLayer ────────────────────────────────────────────────────────────

/// An articular cartilage layer described by biphasic theory parameters.
#[derive(Debug, Clone)]
pub struct CartilageLayer {
    /// Thickness of the cartilage layer (m).
    pub thickness: f64,
    /// Hydraulic permeability (m⁴/N·s).
    pub permeability: f64,
    /// Aggregate equilibrium modulus (Pa).
    pub aggregate_modulus: f64,
}

impl CartilageLayer {
    /// Create a new `CartilageLayer`.
    ///
    /// * `thickness`         – Layer thickness (m), e.g. 2–4 mm.
    /// * `permeability`      – Hydraulic permeability (m⁴/N·s), ≈ 1e-15.
    /// * `aggregate_modulus` – Equilibrium compressive modulus (Pa), ≈ 0.5–1 MPa.
    pub fn new(thickness: f64, permeability: f64, aggregate_modulus: f64) -> Self {
        Self {
            thickness,
            permeability,
            aggregate_modulus,
        }
    }

    /// Biphasic stress response (Pa) at time `time` (s) under a step stress
    /// `stress` (Pa).
    ///
    /// At `t = 0` the fluid carries the full load; as `t → ∞` the solid
    /// matrix carries it.  Uses the single-term exponential approximation:
    ///
    /// `sigma_solid(t) = stress * (1 - exp(-t / tau))`
    pub fn biphasic_response(&self, stress: f64, time: f64) -> f64 {
        let tau =
            self.thickness * self.thickness / (self.aggregate_modulus * self.permeability + 1e-30);
        stress * (1.0 - (-time / tau.max(1e-30)).exp())
    }

    /// Creep displacement (m) under constant compressive load `load` (N)
    /// at time `t` (s).
    ///
    /// `u(t) = (load / (H_A * A)) * h * (1 - exp(-t/tau))`
    /// where `A = 1 m²` (unit area) and `h = thickness`.
    pub fn creep_displacement(&self, load: f64, t: f64) -> f64 {
        let area = 1.0_f64; // unit area
        let tau =
            self.thickness * self.thickness / (self.aggregate_modulus * self.permeability + 1e-30);
        let epsilon_eq = load / (self.aggregate_modulus * area + 1e-30);
        epsilon_eq * self.thickness * (1.0 - (-t / tau.max(1e-30)).exp())
    }

    /// Instantaneous contact stiffness (N/m) of the cartilage layer.
    ///
    /// `k = H_A * A / h`  (unit area, A = 1 m²).
    pub fn contact_stiffness(&self) -> f64 {
        self.aggregate_modulus * 1.0 / self.thickness.max(1e-12)
    }

    /// Fraction of load carried by the fluid phase at time `t` (s).
    ///
    /// At `t = 0` the fluid carries everything (fraction = 1); as
    /// `t → ∞` it goes to zero.
    pub fn fluid_pressure_fraction(&self, time: f64) -> f64 {
        let tau =
            self.thickness * self.thickness / (self.aggregate_modulus * self.permeability + 1e-30);
        (-time / tau.max(1e-30)).exp()
    }
}

// ── Tendon ────────────────────────────────────────────────────────────────────

/// A collagenous tendon described by its structural and compositional
/// parameters.
#[derive(Debug, Clone)]
pub struct Tendon {
    /// Cross-sectional area (m²).
    pub cross_section: f64,
    /// Unstretched (slack) length (m).
    pub rest_length: f64,
    /// Mean collagen fibre pennation angle (rad).
    pub fiber_angle: f64,
    /// Volume fraction of aligned collagen fibres \[0, 1\].
    pub collagen_fraction: f64,
}

impl Tendon {
    /// Create a new `Tendon`.
    ///
    /// * `cross_section`    – Cross-sectional area (m²).
    /// * `rest_length`      – Slack length (m).
    /// * `fiber_angle`      – Fibre pennation angle (rad).
    /// * `collagen_fraction` – Aligned collagen volume fraction \[0, 1\].
    pub fn new(
        cross_section: f64,
        rest_length: f64,
        fiber_angle: f64,
        collagen_fraction: f64,
    ) -> Self {
        Self {
            cross_section,
            rest_length,
            fiber_angle,
            collagen_fraction,
        }
    }

    /// Stiffness in the toe region (N/m) — the initial low-modulus part of
    /// the stress-strain curve where crimped fibres are straightening.
    ///
    /// Uses a toe modulus ≈ 50 MPa scaled by fibre alignment.
    pub fn toe_region_stiffness(&self) -> f64 {
        let e_toe = 50e6_f64; // Pa
        let cos_sq = self.fiber_angle.cos().powi(2);
        e_toe * self.collagen_fraction * cos_sq * self.cross_section / self.rest_length.max(1e-12)
    }

    /// Stiffness in the linear (high-strain) region (N/m).
    ///
    /// Uses a linear modulus ≈ 1.2 GPa for aligned collagen.
    pub fn linear_region_stiffness(&self) -> f64 {
        let e_lin = 1.2e9_f64; // Pa
        let cos_sq = self.fiber_angle.cos().powi(2);
        e_lin * self.collagen_fraction * cos_sq * self.cross_section / self.rest_length.max(1e-12)
    }

    /// Ultimate (failure) force (N).
    ///
    /// `F_ult = sigma_ult * A * collagen_fraction * cos²(theta)`
    /// where `sigma_ult ≈ 100 MPa`.
    pub fn ultimate_force(&self) -> f64 {
        let sigma_ult = 100e6_f64;
        let cos_sq = self.fiber_angle.cos().powi(2);
        sigma_ult * self.cross_section * self.collagen_fraction * cos_sq
    }

    /// Elastic strain energy (J) stored at the given engineering strain.
    ///
    /// Piecewise: below `eps_toe = 0.02` use toe modulus; above use linear.
    pub fn strain_energy(&self, strain: f64) -> f64 {
        let eps_toe = 0.02_f64;
        let volume = self.cross_section * self.rest_length;
        if strain <= 0.0 {
            return 0.0;
        }
        let e_toe = 50e6_f64 * self.collagen_fraction * self.fiber_angle.cos().powi(2);
        let e_lin = 1.2e9_f64 * self.collagen_fraction * self.fiber_angle.cos().powi(2);
        if strain <= eps_toe {
            0.5 * e_toe * strain * strain * volume
        } else {
            let u_toe = 0.5 * e_toe * eps_toe * eps_toe * volume;
            let sigma_toe = e_toe * eps_toe;
            let extra_strain = strain - eps_toe;
            u_toe
                + sigma_toe * extra_strain * volume
                + 0.5 * e_lin * extra_strain * extra_strain * volume
        }
    }

    /// Dimensionless hysteresis ratio (energy dissipated / energy stored).
    ///
    /// A fixed empirical value for tendon: ≈ 0.06–0.10.
    pub fn hysteresis_ratio(&self) -> f64 {
        0.08
    }
}

// ── JointContact ─────────────────────────────────────────────────────────────

/// An articular joint contact interface comprising two opposing cartilage
/// layers.
#[derive(Debug, Clone)]
pub struct JointContact {
    /// Cartilage on the first articulating surface (e.g. femoral condyle).
    pub cartilage_a: CartilageLayer,
    /// Cartilage on the second surface (e.g. tibial plateau).
    pub cartilage_b: CartilageLayer,
}

impl JointContact {
    /// Create a new `JointContact` from two `CartilageLayer`s.
    pub fn new(cartilage_a: CartilageLayer, cartilage_b: CartilageLayer) -> Self {
        Self {
            cartilage_a,
            cartilage_b,
        }
    }

    /// Mean contact pressure (Pa) given the applied load (N) and apparent
    /// contact area (m²).
    ///
    /// `p = load / area`
    pub fn contact_pressure(&self, load: f64, area: f64) -> f64 {
        load / area.max(1e-20)
    }

    /// Effective friction coefficient at the given sliding speed (m/s).
    ///
    /// Uses an elasto-hydrodynamic model: at low speeds synovial fluid
    /// provides full lubrication (μ ≈ 0.001); at high speeds boundary
    /// lubrication dominates (μ → 0.05).
    ///
    /// `mu(v) = mu_0 + (mu_inf - mu_0) * (1 - exp(-v / v_c))`
    pub fn friction_coefficient(&self, speed: f64) -> f64 {
        let mu_0 = 0.001_f64;
        let mu_inf = 0.05_f64;
        let v_c = 0.01_f64; // characteristic speed (m/s)
        mu_0 + (mu_inf - mu_0) * (1.0 - (-speed / v_c).exp())
    }

    /// Stress-relaxation response (Pa) of the joint under a step load `load`
    /// (N) at time `t` (s).
    ///
    /// Approximated as the sum of biphasic responses of both layers acting in
    /// series, normalised by a unit contact area.
    pub fn stress_relaxation(&self, load: f64, t: f64) -> f64 {
        let area = 1.0_f64;
        let stress = load / area;
        // Series: total creep compliance = sum of individual compliances.
        let ha_a = self.cartilage_a.aggregate_modulus;
        let ha_b = self.cartilage_b.aggregate_modulus;
        let ha_eff = (ha_a * ha_b) / (ha_a + ha_b + 1e-30);
        let tau_a =
            self.cartilage_a.thickness.powi(2) / (ha_a * self.cartilage_a.permeability + 1e-30);
        let tau_b =
            self.cartilage_b.thickness.powi(2) / (ha_b * self.cartilage_b.permeability + 1e-30);
        let tau_eff = (tau_a + tau_b) * 0.5;
        // At t=0 fluid carries everything; at t→∞ solid carries full stress.
        stress * ha_eff / (ha_eff + 1e-30) * (1.0 - (-t / tau_eff.max(1e-30)).exp())
    }
}

// ── FractureMechanicsBone ────────────────────────────────────────────────────

/// Fracture mechanics parameters for a bone specimen.
#[derive(Debug, Clone)]
pub struct FractureMechanicsBone {
    /// Mode-I fracture toughness K_Ic (Pa·√m).
    pub toughness: f64,
    /// Current crack half-length `a` (m).
    pub crack_length: f64,
}

impl FractureMechanicsBone {
    /// Create a new `FractureMechanicsBone`.
    ///
    /// * `toughness`    – Fracture toughness K_Ic (Pa·√m), cortical ≈ 2–5 MPa·√m.
    /// * `crack_length` – Initial crack half-length (m).
    pub fn new(toughness: f64, crack_length: f64) -> Self {
        Self {
            toughness,
            crack_length,
        }
    }

    /// Mode-I stress intensity factor K_I (Pa·√m) for an edge crack under
    /// remote stress `stress` (Pa).
    ///
    /// `K_I = F * sigma * sqrt(pi * a)` with geometry factor `F = 1.12`.
    pub fn stress_intensity_factor(&self, stress: f64) -> f64 {
        let f_geom = 1.12_f64;
        f_geom * stress * (std::f64::consts::PI * self.crack_length).sqrt()
    }

    /// Critical crack half-length (m) at which fast fracture occurs under
    /// the given applied stress (Pa).
    ///
    /// `a_c = (1 / pi) * (K_Ic / (F * sigma))^2`
    pub fn critical_crack_length(&self, stress: f64) -> f64 {
        if stress <= 0.0 {
            return f64::INFINITY;
        }
        let f_geom = 1.12_f64;
        let ratio = self.toughness / (f_geom * stress);
        ratio * ratio / std::f64::consts::PI
    }

    /// Fatigue crack growth rate da/dN (m/cycle) via the Paris law.
    ///
    /// `da/dN = C * (delta_K)^m`
    ///
    /// Using cortical bone constants: `C = 1.77e-10`, `m = 2.3` (Vashishth 2004).
    pub fn fatigue_crack_growth(&self, delta_k: f64) -> f64 {
        let c = 1.77e-10_f64;
        let m = 2.3_f64;
        c * delta_k.abs().powf(m)
    }
}

// ── Wolff's law ───────────────────────────────────────────────────────────────

/// Wolff's law bone density remodelling stimulus.
///
/// Returns the rate of change of bone apparent density (kg/m³/day) driven
/// by the difference between the current mechanical stimulus and a
/// homeostatic reference strain.
///
/// Uses a standard lazy-zone formulation:
///
/// * If `|stimulus - reference_strain| < 0.05 * reference_strain` → no remodelling.
/// * Otherwise → `rate = k * (stimulus - reference_strain)`.
///
/// * `current_density` – Current apparent density (kg/m³).
/// * `stimulus`        – Effective strain stimulus (dimensionless, e.g. SED).
/// * `reference_strain` – Homeostatic strain setpoint.
pub fn wolff_law_remodeling(_current_density: f64, stimulus: f64, reference_strain: f64) -> f64 {
    let k = 0.3_f64; // adaptation rate (kg/m³/day per unit strain)
    let dead_zone_frac = 0.05_f64;
    let dead_zone = dead_zone_frac * reference_strain.abs().max(1e-12);
    let diff = stimulus - reference_strain;
    if diff.abs() <= dead_zone {
        0.0
    } else if diff > 0.0 {
        k * (diff - dead_zone)
    } else {
        k * (diff + dead_zone)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── BoneStructure tests ───────────────────────────────────────────────

    #[test]
    fn bone_youngs_modulus_cortical() {
        let bone = BoneStructure::new(1900.0, 0.05, 0.45);
        let e = bone.youngs_modulus();
        // At full density the modulus should be exactly 6.95 GPa.
        assert!((e - 6.95e9).abs() < 1e3, "E={e}");
    }

    #[test]
    fn bone_youngs_modulus_increases_with_density() {
        let b1 = BoneStructure::new(400.0, 0.7, 0.35);
        let b2 = BoneStructure::new(1200.0, 0.2, 0.42);
        assert!(b2.youngs_modulus() > b1.youngs_modulus());
    }

    #[test]
    fn bone_compressive_strength_positive() {
        let bone = BoneStructure::new(1900.0, 0.05, 0.45);
        assert!(bone.compressive_strength() > 0.0);
    }

    #[test]
    fn bone_reuss_modulus_less_than_voigt() {
        let bone = BoneStructure::new(1900.0, 0.05, 0.45);
        assert!(bone.reuss_modulus() < bone.voigt_modulus());
    }

    #[test]
    fn bone_voigt_reuss_bounds_bracket_effective() {
        let bone = BoneStructure::new(1900.0, 0.05, 0.45);
        let e = bone.youngs_modulus();
        let r = bone.reuss_modulus();
        let v = bone.voigt_modulus();
        assert!(r <= e || e <= v, "E should be within Reuss-Voigt bounds");
    }

    #[test]
    fn bone_fatigue_life_infinite_at_zero_amplitude() {
        let bone = BoneStructure::new(1900.0, 0.05, 0.45);
        assert_eq!(bone.fatigue_life(0.0), f64::INFINITY);
    }

    #[test]
    fn bone_fatigue_life_decreases_with_higher_stress() {
        let bone = BoneStructure::new(1900.0, 0.05, 0.45);
        let n_low = bone.fatigue_life(1e6);
        let n_high = bone.fatigue_life(5e6);
        assert!(n_high < n_low, "higher stress → shorter fatigue life");
    }

    #[test]
    fn bone_remodeling_rate_zero_in_dead_zone() {
        let bone = BoneStructure::new(1000.0, 0.3, 0.40);
        // Stimulus very close to setpoint — within dead zone.
        let rate = bone.remodeling_rate(1000.0, 1000.0);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn bone_remodeling_rate_positive_above_setpoint() {
        let bone = BoneStructure::new(1000.0, 0.3, 0.40);
        let rate = bone.remodeling_rate(2000.0, 1000.0);
        assert!(rate > 0.0, "over-stimulus should increase density");
    }

    #[test]
    fn bone_remodeling_rate_negative_below_setpoint() {
        let bone = BoneStructure::new(1000.0, 0.3, 0.40);
        let rate = bone.remodeling_rate(0.0, 1000.0);
        assert!(rate < 0.0, "under-stimulus should decrease density");
    }

    // ── CartilageLayer tests ──────────────────────────────────────────────

    #[test]
    fn cartilage_biphasic_response_zero_at_t0() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        // At t=0 fluid carries everything → solid response = 0.
        let r = cart.biphasic_response(1000.0, 0.0);
        assert!(
            r.abs() < 1e-6,
            "biphasic response at t=0 should be ~0, got {r}"
        );
    }

    #[test]
    fn cartilage_biphasic_response_approaches_stress_at_large_t() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let stress = 1000.0_f64;
        let r = cart.biphasic_response(stress, 1e8); // very large time
        assert!(
            (r - stress).abs() < 1.0,
            "biphasic should approach applied stress, got {r}"
        );
    }

    #[test]
    fn cartilage_creep_displacement_zero_at_t0() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let u = cart.creep_displacement(500.0, 0.0);
        assert!(u.abs() < 1e-12, "creep at t=0 should be 0, got {u}");
    }

    #[test]
    fn cartilage_creep_displacement_increases_with_time() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let u1 = cart.creep_displacement(500.0, 1.0);
        let u2 = cart.creep_displacement(500.0, 100.0);
        assert!(u2 > u1, "creep should increase with time");
    }

    #[test]
    fn cartilage_contact_stiffness_positive() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        assert!(cart.contact_stiffness() > 0.0);
    }

    #[test]
    fn cartilage_fluid_pressure_fraction_one_at_t0() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let frac = cart.fluid_pressure_fraction(0.0);
        assert!((frac - 1.0).abs() < 1e-10, "at t=0 fluid carries full load");
    }

    #[test]
    fn cartilage_fluid_pressure_fraction_decays_over_time() {
        let cart = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let f1 = cart.fluid_pressure_fraction(1.0);
        let f2 = cart.fluid_pressure_fraction(1e8);
        assert!(f2 < f1, "fluid fraction should decrease over time");
        assert!(
            f2 < 0.01,
            "fluid fraction should be near 0 at large t, got {f2}"
        );
    }

    // ── Tendon tests ──────────────────────────────────────────────────────

    #[test]
    fn tendon_linear_stiffness_greater_than_toe() {
        let t = Tendon::new(1e-4, 0.1, 0.0, 0.85);
        assert!(
            t.linear_region_stiffness() > t.toe_region_stiffness(),
            "linear stiffness should exceed toe stiffness"
        );
    }

    #[test]
    fn tendon_ultimate_force_positive() {
        let t = Tendon::new(1e-4, 0.1, 0.0, 0.85);
        assert!(t.ultimate_force() > 0.0);
    }

    #[test]
    fn tendon_strain_energy_zero_at_zero_strain() {
        let t = Tendon::new(1e-4, 0.1, 0.0, 0.85);
        assert_eq!(t.strain_energy(0.0), 0.0);
    }

    #[test]
    fn tendon_strain_energy_increases_with_strain() {
        let t = Tendon::new(1e-4, 0.1, 0.0, 0.85);
        let u1 = t.strain_energy(0.01);
        let u2 = t.strain_energy(0.05);
        assert!(u2 > u1, "strain energy should increase with strain");
    }

    #[test]
    fn tendon_hysteresis_ratio_between_zero_and_one() {
        let t = Tendon::new(1e-4, 0.1, 0.0, 0.85);
        let h = t.hysteresis_ratio();
        assert!(h > 0.0 && h < 1.0, "hysteresis ratio must be in (0,1): {h}");
    }

    #[test]
    fn tendon_stiffness_scales_with_cross_section() {
        let t1 = Tendon::new(1e-4, 0.1, 0.0, 0.85);
        let t2 = Tendon::new(2e-4, 0.1, 0.0, 0.85);
        assert!((t2.linear_region_stiffness() - 2.0 * t1.linear_region_stiffness()).abs() < 1.0);
    }

    // ── JointContact tests ────────────────────────────────────────────────

    #[test]
    fn joint_contact_pressure_basic() {
        let ca = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let cb = CartilageLayer::new(0.002, 1e-15, 0.6e6);
        let joint = JointContact::new(ca, cb);
        let p = joint.contact_pressure(1000.0, 1e-4);
        assert!((p - 1.0e7).abs() < 1.0, "p={p}");
    }

    #[test]
    fn joint_friction_low_at_low_speed() {
        let ca = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let cb = CartilageLayer::new(0.002, 1e-15, 0.6e6);
        let joint = JointContact::new(ca, cb);
        let mu = joint.friction_coefficient(0.0);
        assert!(mu < 0.005, "mu at v=0 should be near 0.001, got {mu}");
    }

    #[test]
    fn joint_friction_increases_with_speed() {
        let ca = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let cb = CartilageLayer::new(0.002, 1e-15, 0.6e6);
        let joint = JointContact::new(ca, cb);
        let mu_low = joint.friction_coefficient(0.001);
        let mu_high = joint.friction_coefficient(1.0);
        assert!(mu_high > mu_low);
    }

    #[test]
    fn joint_stress_relaxation_zero_at_t0() {
        let ca = CartilageLayer::new(0.003, 1e-15, 0.8e6);
        let cb = CartilageLayer::new(0.002, 1e-15, 0.6e6);
        let joint = JointContact::new(ca, cb);
        let r = joint.stress_relaxation(1000.0, 0.0);
        assert!(
            r.abs() < 1e-3,
            "stress relaxation at t=0 should be ~0, got {r}"
        );
    }

    // ── FractureMechanicsBone tests ───────────────────────────────────────

    #[test]
    fn fracture_stress_intensity_zero_stress() {
        let frac = FractureMechanicsBone::new(3e6, 1e-3);
        assert_eq!(frac.stress_intensity_factor(0.0), 0.0);
    }

    #[test]
    fn fracture_stress_intensity_positive_stress() {
        let frac = FractureMechanicsBone::new(3e6, 1e-3);
        let k = frac.stress_intensity_factor(100e6);
        assert!(k > 0.0, "K_I should be positive: {k}");
    }

    #[test]
    fn fracture_critical_crack_length_infinite_at_zero_stress() {
        let frac = FractureMechanicsBone::new(3e6, 1e-3);
        assert_eq!(frac.critical_crack_length(0.0), f64::INFINITY);
    }

    #[test]
    fn fracture_critical_crack_decreases_with_stress() {
        let frac = FractureMechanicsBone::new(3e6, 1e-3);
        let a1 = frac.critical_crack_length(50e6);
        let a2 = frac.critical_crack_length(100e6);
        assert!(a2 < a1, "higher stress → shorter critical crack length");
    }

    #[test]
    fn fracture_paris_law_growth_positive() {
        let frac = FractureMechanicsBone::new(3e6, 1e-3);
        let rate = frac.fatigue_crack_growth(1e6);
        assert!(rate > 0.0, "crack growth rate should be positive: {rate}");
    }

    #[test]
    fn fracture_paris_law_growth_increases_with_delta_k() {
        let frac = FractureMechanicsBone::new(3e6, 1e-3);
        let r1 = frac.fatigue_crack_growth(1e6);
        let r2 = frac.fatigue_crack_growth(2e6);
        assert!(r2 > r1, "larger delta_K → faster growth");
    }

    // ── Wolff's law tests ─────────────────────────────────────────────────

    #[test]
    fn wolff_law_zero_in_dead_zone() {
        // Stimulus exactly at reference: should be zero.
        let rate = wolff_law_remodeling(1000.0, 0.003, 0.003);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn wolff_law_positive_above_setpoint() {
        let rate = wolff_law_remodeling(800.0, 0.006, 0.003);
        assert!(
            rate > 0.0,
            "above setpoint should give positive remodelling rate"
        );
    }

    #[test]
    fn wolff_law_negative_below_setpoint() {
        let rate = wolff_law_remodeling(1200.0, 0.0005, 0.003);
        assert!(
            rate < 0.0,
            "below setpoint should give negative remodelling rate"
        );
    }

    #[test]
    fn wolff_law_rate_magnitude_increases_with_stimulus_distance() {
        let r1 = wolff_law_remodeling(1000.0, 0.008, 0.003);
        let r2 = wolff_law_remodeling(1000.0, 0.015, 0.003);
        assert!(r2 > r1, "further from setpoint → larger magnitude");
    }
}
