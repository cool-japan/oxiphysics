// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Extended tribology — friction, wear, and lubrication.
//!
//! Provides Hertz contact mechanics, Stribeck lubrication curves, Archard
//! wear model, Sommerfeld number, elasto-hydrodynamic film thickness,
//! flash temperature, and hydrodynamic lift for tribological analysis.

use std::f64::consts::PI;

// ─────────────────────────────────────────────────────────────────────────────
// ContactInterface
// ─────────────────────────────────────────────────────────────────────────────

/// Tribological contact interface describing surface and material properties.
#[derive(Debug, Clone, Copy)]
pub struct ContactInterface {
    /// Root-mean-square surface roughness \[m\].
    pub roughness_rms: f64,
    /// Young's modulus of the contact body \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio of the contact body (dimensionless).
    pub poisson_ratio: f64,
    /// Vickers hardness \[Pa\].
    pub hardness: f64,
}

impl ContactInterface {
    /// Create a new [`ContactInterface`].
    pub fn new(roughness_rms: f64, young_modulus: f64, poisson_ratio: f64, hardness: f64) -> Self {
        Self {
            roughness_rms,
            young_modulus,
            poisson_ratio,
            hardness,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Reduced modulus
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the reduced (combined) elastic modulus E* for two contacting bodies.
///
/// E* = 1 / ((1−ν₁²)/E₁ + (1−ν₂²)/E₂)
///
/// # Arguments
/// * `e1` – Young's modulus of body 1 \[Pa\].
/// * `nu1` – Poisson's ratio of body 1.
/// * `e2` – Young's modulus of body 2 \[Pa\].
/// * `nu2` – Poisson's ratio of body 2.
pub fn reduced_modulus(e1: f64, nu1: f64, e2: f64, nu2: f64) -> f64 {
    1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2)
}

// ─────────────────────────────────────────────────────────────────────────────
// Hertz contact radius
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Hertz contact radius for two spheres in contact.
///
/// a = (3·F·R_eff / (4·E*))^(1/3)
///
/// where R_eff = R₁·R₂/(R₁+R₂) is the effective radius.
///
/// # Arguments
/// * `r1` – Radius of sphere 1 \[m\].
/// * `r2` – Radius of sphere 2 \[m\] (use a very large value for flat).
/// * `e_star` – Reduced modulus E* \[Pa\].
/// * `load` – Normal force F \[N\].
pub fn hertz_contact_radius(r1: f64, r2: f64, e_star: f64, load: f64) -> f64 {
    let r_eff = r1 * r2 / (r1 + r2);
    (3.0 * load * r_eff / (4.0 * e_star)).cbrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Hertz maximum pressure
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Hertz maximum contact pressure (peak Hertzian pressure).
///
/// p₀ = 3F / (2π·a²)
///
/// # Arguments
/// * `load` – Normal force F \[N\].
/// * `contact_radius` – Contact radius a \[m\].
pub fn hertz_max_pressure(load: f64, contact_radius: f64) -> f64 {
    3.0 * load / (2.0 * PI * contact_radius * contact_radius)
}

// ─────────────────────────────────────────────────────────────────────────────
// Stribeck curve
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Stribeck number (dimensionless lubrication parameter).
///
/// S = η·v / (P·R)
///
/// The Stribeck number governs the lubrication regime transition from boundary
/// to mixed to hydrodynamic lubrication.
///
/// # Arguments
/// * `speed` – Sliding speed v \[m/s\].
/// * `viscosity` – Dynamic viscosity η \[Pa·s\].
/// * `pressure` – Mean contact pressure P \[Pa\].
/// * `r` – Characteristic contact radius R \[m\].
pub fn stribeck_curve(speed: f64, viscosity: f64, pressure: f64, r: f64) -> f64 {
    viscosity * speed / (pressure * r)
}

// ─────────────────────────────────────────────────────────────────────────────
// Hydrodynamic lift
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the hydrodynamic lift force from the Reynolds lubrication equation.
///
/// F ≈ η·v·L·W / h²  (simplified wedge-film approximation)
///
/// # Arguments
/// * `viscosity` – Dynamic viscosity η \[Pa·s\].
/// * `speed` – Surface speed v \[m/s\].
/// * `gap` – Minimum film thickness h \[m\].
/// * `length` – Bearing length L \[m\].
/// * `width` – Bearing width W \[m\].
pub fn hydrodynamic_lift(viscosity: f64, speed: f64, gap: f64, length: f64, width: f64) -> f64 {
    viscosity * speed * length * width / (gap * gap)
}

// ─────────────────────────────────────────────────────────────────────────────
// ArchardWear
// ─────────────────────────────────────────────────────────────────────────────

/// Archard adhesive wear model.
///
/// The Archard wear equation: V = k·W·s / H  where V is wear volume,
/// k is the dimensionless wear coefficient, W is normal load, s is sliding
/// distance, and H is hardness.
#[derive(Debug, Clone, Copy)]
pub struct ArchardWear {
    /// Dimensionless Archard wear coefficient k.
    pub k: f64,
    /// Material hardness H \[Pa\].
    pub hardness: f64,
}

impl ArchardWear {
    /// Create a new [`ArchardWear`] model.
    pub fn new(k: f64, hardness: f64) -> Self {
        Self { k, hardness }
    }

    /// Compute wear volume \[m³\] for a given normal load and sliding distance.
    ///
    /// V = k·W·s / H
    ///
    /// # Arguments
    /// * `load` – Normal force W \[N\].
    /// * `sliding_dist` – Total sliding distance s \[m\].
    pub fn wear_volume(&self, load: f64, sliding_dist: f64) -> f64 {
        self.k * load * sliding_dist / self.hardness
    }

    /// Compute wear depth \[m\] for a given normal load, sliding distance, and
    /// nominal contact area.
    ///
    /// d = V / A = k·W·s / (H·A)
    ///
    /// # Arguments
    /// * `load` – Normal force W \[N\].
    /// * `sliding_dist` – Total sliding distance s \[m\].
    /// * `area` – Nominal contact area A \[m²\].
    pub fn wear_depth(&self, load: f64, sliding_dist: f64, area: f64) -> f64 {
        self.wear_volume(load, sliding_dist) / area
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Sommerfeld number
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Sommerfeld number for a journal bearing.
///
/// S = η·ω·R·L / (W·(c/R)²)
///
/// The Sommerfeld number characterises the operating conditions of a
/// hydrodynamic journal bearing.
///
/// # Arguments
/// * `viscosity` – Dynamic viscosity η \[Pa·s\].
/// * `speed` – Rotational speed ω \[rad/s\].
/// * `load` – Bearing load per unit projected area W \[N\].
/// * `r` – Journal radius R \[m\].
/// * `l` – Bearing length L \[m\].
/// * `c` – Radial clearance c \[m\].
pub fn sommerfeld_number(viscosity: f64, speed: f64, load: f64, r: f64, l: f64, c: f64) -> f64 {
    let c_over_r = c / r;
    viscosity * speed * r * l / (load * c_over_r * c_over_r)
}

// ─────────────────────────────────────────────────────────────────────────────
// Elasto-hydrodynamic film thickness
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate the minimum EHL (elasto-hydrodynamic lubrication) film thickness
/// using a simplified Hamrock-Dowson formula.
///
/// h_min ≈ 2.69 · (η₀·u)^0.67 · R'^0.464 / (E'^0.073 · W^0.067)
///
/// where η₀ is the dynamic viscosity, u is the entraining velocity, R' is
/// the effective radius, E' = 2·E* is the combined modulus, and W is the
/// dimensionless load parameter W = load/(E'·R'^2).
///
/// # Arguments
/// * `r_eff` – Effective radius R' \[m\].
/// * `e_star` – Reduced modulus E* \[Pa\].
/// * `u` – Entraining velocity u \[m/s\].
/// * `eta` – Dynamic viscosity η₀ \[Pa·s\].
/// * `load` – Normal force \[N\].
pub fn elasto_hydrodynamic_film_thickness(
    r_eff: f64,
    e_star: f64,
    u: f64,
    eta: f64,
    load: f64,
) -> f64 {
    // E' = 2 * E*
    let e_prime = 2.0 * e_star;
    // Dimensionless load W* = load / (E' * R'^2)
    let w_star = load / (e_prime * r_eff * r_eff);
    // Dimensionless speed U* = eta * u / (E' * R')
    let u_star = eta * u / (e_prime * r_eff);
    // Hamrock-Dowson: h_min/R' = 2.69 * U*^0.67 * W*^(-0.067) * (point contact)
    r_eff * 2.69 * u_star.powf(0.67) * w_star.powf(-0.067)
}

// ─────────────────────────────────────────────────────────────────────────────
// Friction power dissipation
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the frictional power dissipation rate.
///
/// P = F_friction · v_sliding
///
/// # Arguments
/// * `friction_force` – Friction force F \[N\].
/// * `sliding_velocity` – Relative sliding velocity v \[m/s\].
pub fn friction_power_dissipation(friction_force: f64, sliding_velocity: f64) -> f64 {
    friction_force * sliding_velocity
}

// ─────────────────────────────────────────────────────────────────────────────
// Flash temperature
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Blok flash temperature rise at the contact interface.
///
/// Δθ_flash ≈ q·r_contact / (kappa·sqrt(π·kappa_diff / (v·r_contact)))
///
/// Simplified form: Δθ ≈ 0.31 · q / (kappa · sqrt(r_c · v / κ))
/// where κ = kappa/(rho·cp) is thermal diffusivity.
///
/// # Arguments
/// * `q_friction` – Frictional heat flux Q \[W/m²\].
/// * `r_contact` – Contact radius r_c \[m\].
/// * `kappa` – Thermal conductivity κ \[W/(m·K)\].
/// * `rho_cp` – Volumetric heat capacity ρ·cp \[J/(m³·K)\].
/// * `v` – Sliding velocity v \[m/s\].
pub fn flash_temperature(q_friction: f64, r_contact: f64, kappa: f64, rho_cp: f64, v: f64) -> f64 {
    // Thermal diffusivity κ_diff = kappa / (rho * cp)
    let kappa_diff = kappa / rho_cp;
    // Blok flash temperature: ΔT = q * sqrt(π * r_c / (v * κ_diff)) / (2 * kappa)
    let peclet_term = (PI * r_contact / (v * kappa_diff)).sqrt();
    q_friction * peclet_term / (2.0 * kappa)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. reduced_modulus is positive for valid inputs.
    #[test]
    fn test_reduced_modulus_positive() {
        let e_star = reduced_modulus(200e9, 0.3, 200e9, 0.3);
        assert!(
            e_star > 0.0,
            "Reduced modulus must be positive, got {e_star}"
        );
    }

    // 2. reduced_modulus for equal materials equals E/(1-nu^2).
    #[test]
    fn test_reduced_modulus_equal_materials() {
        let e = 200e9_f64;
        let nu = 0.3_f64;
        let e_star = reduced_modulus(e, nu, e, nu);
        let expected = e / (2.0 * (1.0 - nu * nu));
        assert!(
            (e_star - expected).abs() / expected < 1e-9,
            "Equal materials: E* = {e_star} vs expected {expected}"
        );
    }

    // 3. reduced_modulus increases when one modulus increases.
    #[test]
    fn test_reduced_modulus_increases_with_e() {
        let e1 = reduced_modulus(100e9, 0.3, 200e9, 0.3);
        let e2 = reduced_modulus(200e9, 0.3, 200e9, 0.3);
        assert!(e2 > e1, "Larger E should give larger E*");
    }

    // 4. hertz_contact_radius is positive.
    #[test]
    fn test_hertz_contact_radius_positive() {
        let a = hertz_contact_radius(0.01, 0.01, 100e9, 100.0);
        assert!(a > 0.0, "Contact radius must be positive, got {a}");
    }

    // 5. hertz_contact_radius increases with load.
    #[test]
    fn test_hertz_contact_radius_increases_with_load() {
        let a1 = hertz_contact_radius(0.01, 0.01, 100e9, 100.0);
        let a2 = hertz_contact_radius(0.01, 0.01, 100e9, 1000.0);
        assert!(
            a2 > a1,
            "Contact radius should increase with load: a1={a1}, a2={a2}"
        );
    }

    // 6. hertz_contact_radius scales as load^(1/3).
    #[test]
    fn test_hertz_contact_radius_scaling() {
        let r1 = 0.01_f64;
        let e_star = 100e9_f64;
        let a1 = hertz_contact_radius(r1, r1, e_star, 100.0);
        let a8 = hertz_contact_radius(r1, r1, e_star, 800.0);
        // a scales as F^(1/3): a(8F)/a(F) = 2
        let ratio = a8 / a1;
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "Contact radius scales as F^1/3, got ratio={ratio}"
        );
    }

    // 7. hertz_max_pressure formula: p0 = 3F/(2*pi*a^2).
    #[test]
    fn test_hertz_max_pressure_formula() {
        let load = 100.0_f64;
        let a = 1e-3_f64;
        let p0 = hertz_max_pressure(load, a);
        let expected = 3.0 * load / (2.0 * PI * a * a);
        assert!(
            (p0 - expected).abs() < EPS,
            "Max pressure formula mismatch: {p0} vs {expected}"
        );
    }

    // 8. hertz_max_pressure is positive.
    #[test]
    fn test_hertz_max_pressure_positive() {
        let p0 = hertz_max_pressure(50.0, 5e-4);
        assert!(p0 > 0.0, "Max pressure should be positive, got {p0}");
    }

    // 9. hertz_max_pressure increases with load.
    #[test]
    fn test_hertz_max_pressure_increases_with_load() {
        let a = 1e-3_f64;
        let p1 = hertz_max_pressure(100.0, a);
        let p2 = hertz_max_pressure(200.0, a);
        assert!(p2 > p1, "Higher load → higher max pressure");
    }

    // 10. stribeck_curve is positive.
    #[test]
    fn test_stribeck_positive() {
        let s = stribeck_curve(1.0, 0.001, 1e6, 1e-3);
        assert!(s > 0.0, "Stribeck number must be positive, got {s}");
    }

    // 11. stribeck_curve increases with speed.
    #[test]
    fn test_stribeck_increases_with_speed() {
        let s1 = stribeck_curve(1.0, 0.001, 1e6, 1e-3);
        let s2 = stribeck_curve(10.0, 0.001, 1e6, 1e-3);
        assert!(s2 > s1, "Stribeck number should increase with speed");
    }

    // 12. stribeck_curve scales linearly with viscosity.
    #[test]
    fn test_stribeck_scales_with_viscosity() {
        let s1 = stribeck_curve(1.0, 0.001, 1e6, 1e-3);
        let s2 = stribeck_curve(1.0, 0.002, 1e6, 1e-3);
        assert!(
            (s2 / s1 - 2.0).abs() < EPS,
            "Stribeck scales with viscosity: {s2} vs 2*{s1}"
        );
    }

    // 13. hydrodynamic_lift is positive.
    #[test]
    fn test_hydro_lift_positive() {
        let f = hydrodynamic_lift(0.001, 5.0, 1e-5, 0.1, 0.05);
        assert!(f > 0.0, "Hydrodynamic lift must be positive, got {f}");
    }

    // 14. hydrodynamic_lift increases with speed.
    #[test]
    fn test_hydro_lift_increases_with_speed() {
        let f1 = hydrodynamic_lift(0.001, 1.0, 1e-5, 0.1, 0.05);
        let f2 = hydrodynamic_lift(0.001, 2.0, 1e-5, 0.1, 0.05);
        assert!(f2 > f1, "Lift should increase with speed");
    }

    // 15. hydrodynamic_lift scales inversely with gap^2.
    #[test]
    fn test_hydro_lift_gap_scaling() {
        let f1 = hydrodynamic_lift(0.001, 1.0, 1e-5, 0.1, 0.05);
        let f2 = hydrodynamic_lift(0.001, 1.0, 2e-5, 0.1, 0.05);
        // F ∝ 1/h^2, so doubling h gives 1/4 force
        let ratio = f1 / f2;
        assert!(
            (ratio - 4.0).abs() < 1e-6,
            "Lift scales as 1/h^2: got ratio={ratio}"
        );
    }

    // 16. ArchardWear::wear_volume scales with load.
    #[test]
    fn test_archard_wear_scales_with_load() {
        let wear = ArchardWear::new(1e-4, 1e9);
        let v1 = wear.wear_volume(100.0, 1.0);
        let v2 = wear.wear_volume(200.0, 1.0);
        assert!(
            (v2 / v1 - 2.0).abs() < EPS,
            "Wear volume should scale with load: v1={v1}, v2={v2}"
        );
    }

    // 17. ArchardWear::wear_volume scales with sliding distance.
    #[test]
    fn test_archard_wear_scales_with_distance() {
        let wear = ArchardWear::new(1e-4, 1e9);
        let v1 = wear.wear_volume(100.0, 1.0);
        let v2 = wear.wear_volume(100.0, 5.0);
        assert!(
            (v2 / v1 - 5.0).abs() < EPS,
            "Wear volume should scale with sliding distance: {v2} vs 5*{v1}"
        );
    }

    // 18. ArchardWear::wear_volume is positive.
    #[test]
    fn test_archard_wear_positive() {
        let wear = ArchardWear::new(1e-5, 5e9);
        let v = wear.wear_volume(500.0, 10.0);
        assert!(v > 0.0, "Wear volume must be positive, got {v}");
    }

    // 19. ArchardWear::wear_depth is consistent with wear_volume.
    #[test]
    fn test_archard_wear_depth_consistent() {
        let wear = ArchardWear::new(1e-4, 1e9);
        let area = 1e-4_f64; // 1 cm^2
        let vol = wear.wear_volume(100.0, 0.5);
        let depth = wear.wear_depth(100.0, 0.5, area);
        assert!(
            (depth - vol / area).abs() < EPS,
            "Wear depth must be vol/area: {depth} vs {}",
            vol / area
        );
    }

    // 20. ArchardWear: harder material gives less wear.
    #[test]
    fn test_archard_harder_less_wear() {
        let soft = ArchardWear::new(1e-4, 1e9);
        let hard = ArchardWear::new(1e-4, 5e9);
        let v_soft = soft.wear_volume(100.0, 1.0);
        let v_hard = hard.wear_volume(100.0, 1.0);
        assert!(v_hard < v_soft, "Harder material should wear less");
    }

    // 21. sommerfeld_number is positive.
    #[test]
    fn test_sommerfeld_positive() {
        let s = sommerfeld_number(0.01, 100.0, 1000.0, 0.05, 0.05, 0.0001);
        assert!(s > 0.0, "Sommerfeld number must be positive, got {s}");
    }

    // 22. sommerfeld_number increases with viscosity.
    #[test]
    fn test_sommerfeld_increases_with_viscosity() {
        let s1 = sommerfeld_number(0.01, 100.0, 1000.0, 0.05, 0.05, 0.0001);
        let s2 = sommerfeld_number(0.02, 100.0, 1000.0, 0.05, 0.05, 0.0001);
        assert!(s2 > s1, "Sommerfeld increases with viscosity");
    }

    // 23. sommerfeld_number scales linearly with speed.
    #[test]
    fn test_sommerfeld_scales_with_speed() {
        let s1 = sommerfeld_number(0.01, 100.0, 1000.0, 0.05, 0.05, 0.0001);
        let s2 = sommerfeld_number(0.01, 200.0, 1000.0, 0.05, 0.05, 0.0001);
        assert!(
            (s2 / s1 - 2.0).abs() < EPS,
            "Sommerfeld scales with speed: got ratio {}",
            s2 / s1
        );
    }

    // 24. elasto_hydrodynamic_film_thickness is positive.
    #[test]
    fn test_ehl_film_positive() {
        let h = elasto_hydrodynamic_film_thickness(0.01, 100e9, 1.0, 0.01, 1000.0);
        assert!(h > 0.0, "EHL film thickness must be positive, got {h}");
    }

    // 25. elasto_hydrodynamic_film_thickness increases with speed.
    #[test]
    fn test_ehl_film_increases_with_speed() {
        let h1 = elasto_hydrodynamic_film_thickness(0.01, 100e9, 0.5, 0.01, 1000.0);
        let h2 = elasto_hydrodynamic_film_thickness(0.01, 100e9, 1.0, 0.01, 1000.0);
        assert!(h2 > h1, "EHL film increases with speed: h1={h1}, h2={h2}");
    }

    // 26. friction_power_dissipation is positive for positive inputs.
    #[test]
    fn test_friction_power_positive() {
        let p = friction_power_dissipation(50.0, 2.0);
        assert!(p > 0.0, "Friction power must be positive, got {p}");
    }

    // 27. friction_power_dissipation formula P = F * v.
    #[test]
    fn test_friction_power_formula() {
        let f = 75.0_f64;
        let v = 3.5_f64;
        let p = friction_power_dissipation(f, v);
        assert!(
            (p - f * v).abs() < EPS,
            "P = F*v: got {p}, expected {}",
            f * v
        );
    }

    // 28. friction_power_dissipation scales linearly with force.
    #[test]
    fn test_friction_power_scales_with_force() {
        let p1 = friction_power_dissipation(10.0, 2.0);
        let p2 = friction_power_dissipation(20.0, 2.0);
        assert!(
            (p2 / p1 - 2.0).abs() < EPS,
            "Friction power scales with force"
        );
    }

    // 29. flash_temperature is positive.
    #[test]
    fn test_flash_temperature_positive() {
        let dt = flash_temperature(1e6, 1e-4, 50.0, 3e6, 1.0);
        assert!(dt > 0.0, "Flash temperature must be positive, got {dt}");
    }

    // 30. flash_temperature increases with friction heat flux.
    #[test]
    fn test_flash_temperature_increases_with_flux() {
        let dt1 = flash_temperature(1e6, 1e-4, 50.0, 3e6, 1.0);
        let dt2 = flash_temperature(2e6, 1e-4, 50.0, 3e6, 1.0);
        assert!(dt2 > dt1, "Higher heat flux → higher flash temperature");
    }

    // 31. flash_temperature decreases with higher thermal conductivity.
    #[test]
    fn test_flash_temperature_decreases_with_conductivity() {
        let dt1 = flash_temperature(1e6, 1e-4, 50.0, 3e6, 1.0);
        let dt2 = flash_temperature(1e6, 1e-4, 100.0, 3e6, 1.0);
        assert!(dt2 < dt1, "Higher conductivity → lower flash temperature");
    }

    // 32. ContactInterface stores fields correctly.
    #[test]
    fn test_contact_interface_fields() {
        let ci = ContactInterface::new(1e-6, 200e9, 0.3, 2e9);
        assert!((ci.roughness_rms - 1e-6).abs() < EPS);
        assert!((ci.young_modulus - 200e9).abs() < 1.0);
        assert!((ci.poisson_ratio - 0.3).abs() < EPS);
        assert!((ci.hardness - 2e9).abs() < 1.0);
    }

    // 33. hertz_contact_radius increases with larger radius.
    #[test]
    fn test_hertz_contact_radius_larger_sphere() {
        let a1 = hertz_contact_radius(0.01, 0.01, 100e9, 100.0);
        let a2 = hertz_contact_radius(0.02, 0.02, 100e9, 100.0);
        assert!(a2 > a1, "Larger spheres → larger contact radius");
    }

    // 34. ArchardWear::new stores k and hardness.
    #[test]
    fn test_archard_new_fields() {
        let w = ArchardWear::new(2e-4, 3e9);
        assert!((w.k - 2e-4).abs() < EPS);
        assert!((w.hardness - 3e9).abs() < 1.0);
    }

    // 35. stribeck_curve is zero when speed is zero.
    #[test]
    fn test_stribeck_zero_speed() {
        let s = stribeck_curve(0.0, 0.001, 1e6, 1e-3);
        assert_eq!(s, 0.0, "Stribeck = 0 at zero speed");
    }

    // 36. hertz_max_pressure decreases as contact radius grows.
    #[test]
    fn test_hertz_pressure_decreases_with_radius() {
        let p1 = hertz_max_pressure(100.0, 1e-3);
        let p2 = hertz_max_pressure(100.0, 2e-3);
        assert!(
            p2 < p1,
            "Larger contact radius → lower pressure for same load"
        );
    }

    // 37. sommerfeld_number decreases with higher load.
    #[test]
    fn test_sommerfeld_decreases_with_load() {
        let s1 = sommerfeld_number(0.01, 100.0, 500.0, 0.05, 0.05, 0.0001);
        let s2 = sommerfeld_number(0.01, 100.0, 1000.0, 0.05, 0.05, 0.0001);
        assert!(s2 < s1, "Higher load → lower Sommerfeld number");
    }

    // 38. ehl_film_thickness increases with viscosity.
    #[test]
    fn test_ehl_film_increases_with_viscosity() {
        let h1 = elasto_hydrodynamic_film_thickness(0.01, 100e9, 1.0, 0.01, 1000.0);
        let h2 = elasto_hydrodynamic_film_thickness(0.01, 100e9, 1.0, 0.02, 1000.0);
        assert!(h2 > h1, "Higher viscosity → thicker EHL film");
    }
}
