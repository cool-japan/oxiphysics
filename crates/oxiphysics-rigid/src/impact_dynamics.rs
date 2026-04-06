#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Impact and collision dynamics.
//!
//! Provides models for Hertz elastic impact, ballistic penetration, stress-wave
//! propagation, crash simulation, and common impact coefficients.
//!
//! # Example
//!
//! ```no_run
//! use oxiphysics_rigid::impact_dynamics::{HertzImpact, coefficient_of_restitution, impulse};
//!
//! let hertz = HertzImpact::new(1.0e9, 0.01);
//! let force = hertz.compute_force(0.001);
//! assert!(force > 0.0);
//!
//! let cor = coefficient_of_restitution(5.0, 3.0);
//! assert!((cor - 0.6).abs() < 1e-10);
//!
//! let imp = impulse(2.0, [3.0, 0.0, 0.0]);
//! assert!((imp[0] - 6.0).abs() < 1e-10);
//! ```

#![allow(dead_code)]

use std::f64::consts::PI;

// ── vector helpers ────────────────────────────────────────────────────────────

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

// ── ImpactEvent ───────────────────────────────────────────────────────────────

/// Describes a collision event between two bodies.
#[derive(Debug, Clone)]
pub struct ImpactEvent {
    /// Identifier for the first body (index or handle).
    pub body_a: usize,
    /// Identifier for the second body.
    pub body_b: usize,
    /// World-space contact point \[m\].
    pub contact_point: [f64; 3],
    /// Outward contact normal (unit vector from A to B).
    pub contact_normal: [f64; 3],
    /// Relative velocity at the contact point \[m/s\].
    pub relative_velocity: [f64; 3],
    /// Coefficient of restitution (0 = fully plastic, 1 = fully elastic).
    pub restitution: f64,
    /// Coulomb friction coefficient.
    pub friction: f64,
}

impl ImpactEvent {
    /// Create a new impact event with the given parameters.
    pub fn new(
        body_a: usize,
        body_b: usize,
        contact_point: [f64; 3],
        contact_normal: [f64; 3],
        relative_velocity: [f64; 3],
        restitution: f64,
        friction: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            contact_point,
            contact_normal,
            relative_velocity,
            restitution,
            friction,
        }
    }

    /// Compute the normal component of the relative velocity (closing speed, m/s).
    pub fn closing_speed(&self) -> f64 {
        let dot = self.relative_velocity[0] * self.contact_normal[0]
            + self.relative_velocity[1] * self.contact_normal[1]
            + self.relative_velocity[2] * self.contact_normal[2];
        -dot // positive means approaching
    }
}

// ── HertzImpact ───────────────────────────────────────────────────────────────

/// Hertz elastic contact model for sphere-on-sphere or sphere-on-flat impacts.
///
/// The force–penetration relationship is `F = (4/3) * E* * sqrt(R*) * delta^(3/2)`.
#[derive(Debug, Clone)]
pub struct HertzImpact {
    /// Effective elastic modulus `E*` \[Pa\].
    pub elastic_modulus_eff: f64,
    /// Effective contact radius `R*` \[m\].
    pub radius_eff: f64,
}

impl HertzImpact {
    /// Create a new `HertzImpact` model.
    ///
    /// * `elastic_modulus_eff` – Combined elastic modulus E* \[Pa\].
    /// * `radius_eff` – Effective radius R* \[m\].
    pub fn new(elastic_modulus_eff: f64, radius_eff: f64) -> Self {
        Self {
            elastic_modulus_eff,
            radius_eff,
        }
    }

    /// Compute contact force \[N\] for a given penetration depth \[m\].
    ///
    /// Returns zero for non-positive penetration.
    pub fn compute_force(&self, penetration: f64) -> f64 {
        if penetration <= 0.0 {
            return 0.0;
        }
        (4.0 / 3.0) * self.elastic_modulus_eff * self.radius_eff.sqrt() * penetration.powf(1.5)
    }

    /// Estimate the Hertz contact duration \[s\] for a sphere impact.
    ///
    /// Uses the classical formula: `t_c = 2.87 * (m^2 / (R* E*^2 * v0))^(1/5)`.
    ///
    /// * `mass` – Combined/reduced mass \[kg\].
    /// * `initial_velocity` – Approach velocity magnitude \[m/s\].
    pub fn compute_contact_duration(&self, mass: f64, initial_velocity: f64) -> f64 {
        if initial_velocity <= 0.0 || mass <= 0.0 {
            return 0.0;
        }
        let e_star = self.elastic_modulus_eff;
        let r_star = self.radius_eff;
        // t_c = 2.87 * (m^2 / (R* * E*^2 * v0))^(1/5)
        let arg = mass * mass / (r_star * e_star * e_star * initial_velocity);
        2.87 * arg.powf(0.2)
    }
}

// ── BallisticImpact ───────────────────────────────────────────────────────────

/// Ballistic impact model for projectile penetration.
///
/// Uses a simplified energy-balance penetration model based on the
/// projectile's kinetic energy and the target's resistance.
#[derive(Debug, Clone)]
pub struct BallisticImpact {
    /// Projectile mass \[kg\].
    pub projectile_mass: f64,
    /// Projectile velocity vector \[m/s\].
    pub velocity: [f64; 3],
    /// Target hardness / resistance stress \[Pa\].
    pub target_hardness: f64,
}

impl BallisticImpact {
    /// Create a new `BallisticImpact` instance.
    pub fn new(projectile_mass: f64, velocity: [f64; 3], target_hardness: f64) -> Self {
        Self {
            projectile_mass,
            velocity,
            target_hardness,
        }
    }

    /// Compute the speed magnitude \[m/s\] of the projectile.
    pub fn speed(&self) -> f64 {
        vec3_norm(self.velocity)
    }

    /// Estimate the penetration depth \[m\] for a target of given thickness \[m\].
    ///
    /// Uses the energy-balance formula: `p = KE / (sigma_y * A)` where A is a
    /// representative area derived from the projectile's kinetic energy and the
    /// target hardness.  Returns the minimum of the computed depth and the
    /// target thickness.
    ///
    /// This is a simplified model; real penetration depends on geometry.
    pub fn penetration_depth(&self, thickness: f64) -> f64 {
        let v = self.speed();
        let ke = 0.5 * self.projectile_mass * v * v;
        if self.target_hardness <= 0.0 {
            return thickness;
        }
        // Simplified: depth = KE / (H * cross_section_area)
        // Assume cross_section_area ~ 1 cm^2 = 1e-4 m^2 reference
        let area_ref = 1.0e-4_f64;
        let depth = ke / (self.target_hardness * area_ref);
        depth.min(thickness)
    }

    /// Compute the residual velocity \[m/s\] after perforation of a plate.
    ///
    /// Returns 0 if the projectile is stopped (depth < thickness).
    pub fn residual_velocity(&self, thickness: f64) -> f64 {
        let v = self.speed();
        let ke = 0.5 * self.projectile_mass * v * v;
        let area_ref = 1.0e-4_f64;
        let energy_absorbed = self.target_hardness * area_ref * thickness;
        let ke_residual = ke - energy_absorbed;
        if ke_residual <= 0.0 {
            return 0.0;
        }
        (2.0 * ke_residual / self.projectile_mass).sqrt()
    }
}

// ── WaveImpact ────────────────────────────────────────────────────────────────

/// Stress-wave propagation model for impact-induced waves in materials.
#[derive(Debug, Clone)]
pub struct WaveImpact {
    /// Longitudinal stress-wave speed in the medium \[m/s\].
    pub stress_wave_speed: f64,
    /// Specific acoustic impedance of the primary medium Z = rho * c \[Pa·s/m\].
    pub impedance: f64,
}

impl WaveImpact {
    /// Create a new `WaveImpact` model.
    pub fn new(stress_wave_speed: f64, impedance: f64) -> Self {
        Self {
            stress_wave_speed,
            impedance,
        }
    }

    /// Compute the reflection coefficient at an interface between media with
    /// impedances `z1` and `z2`.
    ///
    /// `R = (Z2 - Z1) / (Z2 + Z1)`
    ///
    /// Range: −1 (free surface) to +1 (rigid wall).
    #[allow(non_snake_case)]
    pub fn reflection_coeff(z1: f64, z2: f64) -> f64 {
        if (z1 + z2).abs() < 1e-30 {
            return 0.0;
        }
        (z2 - z1) / (z2 + z1)
    }

    /// Compute the transmission coefficient at the same interface.
    ///
    /// `T = 2 * Z1 / (Z1 + Z2)`
    #[allow(non_snake_case)]
    pub fn transmission_coeff(z1: f64, z2: f64) -> f64 {
        if (z1 + z2).abs() < 1e-30 {
            return 0.0;
        }
        2.0 * z1 / (z1 + z2)
    }

    /// Compute the wave arrival time \[s\] at distance `x` \[m\].
    pub fn arrival_time(&self, x: f64) -> f64 {
        if self.stress_wave_speed <= 0.0 {
            return f64::INFINITY;
        }
        x / self.stress_wave_speed
    }

    /// Compute peak stress \[Pa\] using the Hugoniot relation for a velocity
    /// step `delta_v` \[m/s\]: `sigma = Z * delta_v`.
    pub fn peak_stress(&self, delta_v: f64) -> f64 {
        self.impedance * delta_v
    }
}

// ── CrashSimulation ───────────────────────────────────────────────────────────

/// Simplified crash energy absorption model for vehicle structural analysis.
#[derive(Debug, Clone)]
pub struct CrashSimulation {
    /// Total collision (kinetic) energy at impact \[J\].
    pub collision_energy: f64,
    /// Energy absorbed by structural deformation \[J\].
    pub deformation_energy: f64,
    /// Crush/deformation distance \[m\].
    pub crush_distance: f64,
    /// Mean crush force \[N\].
    pub mean_force: f64,
    /// Structural mass involved \[kg\].
    pub structural_mass: f64,
}

impl CrashSimulation {
    /// Create a new `CrashSimulation`.
    pub fn new(
        collision_energy: f64,
        deformation_energy: f64,
        crush_distance: f64,
        mean_force: f64,
        structural_mass: f64,
    ) -> Self {
        Self {
            collision_energy,
            deformation_energy,
            crush_distance,
            mean_force,
            structural_mass,
        }
    }

    /// Compute the specific energy absorption (SEA) \[J/kg\].
    ///
    /// SEA = deformation_energy / structural_mass.
    pub fn specific_energy_absorption(&self) -> f64 {
        if self.structural_mass <= 0.0 {
            return 0.0;
        }
        self.deformation_energy / self.structural_mass
    }

    /// Fraction of collision energy absorbed (0–1).
    pub fn energy_absorption_ratio(&self) -> f64 {
        if self.collision_energy <= 0.0 {
            return 0.0;
        }
        (self.deformation_energy / self.collision_energy).min(1.0)
    }

    /// Estimate the mean deceleration \[m/s²\] given a vehicle mass \[kg\] and
    /// a closing speed \[m/s\].
    pub fn mean_deceleration(&self, vehicle_mass: f64, closing_speed: f64) -> f64 {
        if vehicle_mass <= 0.0 || closing_speed <= 0.0 {
            return 0.0;
        }
        // Using work-energy theorem: F * d = 0.5 * m * v^2
        let f = 0.5 * vehicle_mass * closing_speed * closing_speed / self.crush_distance.max(1e-9);
        f / vehicle_mass
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Compute the coefficient of restitution from speeds before and after impact.
///
/// COR = |v_after| / |v_before|.  Clamps to \[0, 1\].
pub fn coefficient_of_restitution(v_before: f64, v_after: f64) -> f64 {
    if v_before.abs() < 1e-30 {
        return 0.0;
    }
    (v_after.abs() / v_before.abs()).min(1.0)
}

/// Compute the impulse vector \[N·s\] = mass × Δv.
pub fn impulse(mass: f64, delta_v: [f64; 3]) -> [f64; 3] {
    vec3_scale(delta_v, mass)
}

/// Compute Johnson's impact number (dimensionless), which characterises
/// whether an impact is elastic or plastic.
///
/// `J = rho * v^2 / sigma_y`
///
/// * `rho`     – material density \[kg/m³\].
/// * `v`       – impact velocity \[m/s\].
/// * `sigma_y` – yield stress \[Pa\].
///
/// Values > 1 indicate plastic deformation.
pub fn johnson_impact_number(rho: f64, v: f64, sigma_y: f64) -> f64 {
    if sigma_y <= 0.0 {
        return f64::INFINITY;
    }
    rho * v * v / sigma_y
}

/// Compute the maximum penetration depth in Hertz contact \[m\].
///
/// `delta_max = (15 m v^2 / (16 E* sqrt(R*)))^(2/5)`
pub fn hertz_max_penetration(elastic_modulus_eff: f64, radius_eff: f64, mass: f64, v: f64) -> f64 {
    let numerator = 15.0 * mass * v * v;
    let denominator = 16.0 * elastic_modulus_eff * radius_eff.sqrt();
    if denominator <= 0.0 {
        return 0.0;
    }
    (numerator / denominator).powf(0.4)
}

/// Compute the post-impact velocities for a 1-D central collision.
///
/// Returns `(v_a_after, v_b_after)` using the conservation of momentum and
/// COR equations.
///
/// * `ma`, `mb` – masses of bodies A and B.
/// * `va`, `vb` – pre-impact velocities (1-D, signed).
/// * `e`        – coefficient of restitution \[0, 1\].
pub fn central_impact_velocities(ma: f64, mb: f64, va: f64, vb: f64, e: f64) -> (f64, f64) {
    let total_mass = ma + mb;
    if total_mass < 1e-30 {
        return (va, vb);
    }
    let va_new = (ma * va + mb * vb + mb * e * (vb - va)) / total_mass;
    let vb_new = (ma * va + mb * vb + ma * e * (va - vb)) / total_mass;
    (va_new, vb_new)
}

/// Compute the contact area radius \[m\] in Hertz contact.
///
/// `a = (R* * delta)^(1/2)`
pub fn hertz_contact_radius(radius_eff: f64, penetration: f64) -> f64 {
    if radius_eff <= 0.0 || penetration <= 0.0 {
        return 0.0;
    }
    (radius_eff * penetration).sqrt()
}

/// Compute the maximum contact pressure \[Pa\] in Hertz contact.
///
/// `p_max = (3/2) * F / (pi * a^2)`
pub fn hertz_max_pressure(force: f64, contact_radius: f64) -> f64 {
    let area = PI * contact_radius * contact_radius;
    if area < 1e-30 {
        return 0.0;
    }
    1.5 * force / area
}

// ── New impact dynamics functions (extended API) ──────────────────────────────

/// Estimate the coefficient of restitution for a named material.
///
/// Known materials: `"steel"` → 0.7, `"rubber"` → 0.8, `"clay"` → 0.1,
/// `"glass"` → 0.65. Returns 0.5 for unknown materials.
pub fn coefficient_of_restitution_estimate(material: &str) -> f64 {
    match material {
        "steel" => 0.7,
        "rubber" => 0.8,
        "clay" => 0.1,
        "glass" => 0.65,
        _ => 0.5,
    }
}

/// Compute Hertz contact force \[N\] using `F = k * delta^exponent`.
///
/// For sphere-sphere contact `exponent = 1.5` (3/2 power law).
///
/// * `k_hertz`  – contact stiffness coefficient.
/// * `delta`    – penetration/overlap \[m\].
/// * `exponent` – power law exponent (typically 1.5).
pub fn hertz_impact_force(k_hertz: f64, delta: f64, exponent: f64) -> f64 {
    if delta <= 0.0 {
        return 0.0;
    }
    k_hertz * delta.powf(exponent)
}

/// Compute Hertz contact stiffness `k = (4/3) * E* * sqrt(R*)`.
///
/// * `e1`, `nu1` – Young's modulus \[Pa\] and Poisson ratio of body 1.
/// * `e2`, `nu2` – Young's modulus \[Pa\] and Poisson ratio of body 2.
/// * `r1`, `r2`  – Radii of curvature \[m\] of the two bodies.
///
/// Returns `k` \[N/m^(3/2)\].
pub fn hertz_contact_stiffness(e1: f64, nu1: f64, e2: f64, nu2: f64, r1: f64, r2: f64) -> f64 {
    // Effective modulus: 1/E* = (1-nu1^2)/E1 + (1-nu2^2)/E2
    let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
    // Effective radius: 1/R* = 1/R1 + 1/R2
    let r_star = 1.0 / (1.0 / r1 + 1.0 / r2);
    (4.0 / 3.0) * e_star * r_star.sqrt()
}

/// Estimate impact duration for Hertz sphere impact \[s\].
///
/// Approximate formula: `t = 2.94 * sqrt(m / (k^(2/5) * v0^(1/5)))`.
///
/// * `mass` – reduced mass \[kg\].
/// * `k`    – Hertz contact stiffness.
/// * `v0`   – impact velocity \[m/s\].
pub fn impact_duration_hertz(mass: f64, k: f64, v0: f64) -> f64 {
    if mass <= 0.0 || k <= 0.0 || v0 <= 0.0 {
        return 0.0;
    }
    2.94 * (mass / (k.powf(0.4) * v0.powf(0.2))).sqrt()
}

/// Estimate peak Hertz impact force \[N\].
///
/// `F_peak ≈ k * (m * v0^2 / k)^0.6`
///
/// * `mass` – reduced/effective mass \[kg\].
/// * `v0`   – impact velocity \[m/s\].
/// * `k`    – Hertz stiffness.
pub fn peak_impact_force(mass: f64, v0: f64, k: f64) -> f64 {
    if k <= 0.0 || mass <= 0.0 || v0 <= 0.0 {
        return 0.0;
    }
    k * (mass * v0 * v0 / k).powf(0.6)
}

/// Compute the scalar impulse from Newton's restitution law for 1-D collision.
///
/// `j = -(1+e) * (v1 - v2) / (1/m1 + 1/m2)`
///
/// Positive `j` means body 1 receives a positive impulse.
pub fn impulse_from_restitution(m1: f64, m2: f64, v1: f64, v2: f64, e: f64) -> f64 {
    let denom = 1.0 / m1 + 1.0 / m2;
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    -(1.0 + e) * (v1 - v2) / denom
}

/// Compute post-impact velocities for a 1-D elastic-inelastic collision.
///
/// Returns `(v1_after, v2_after)`.
///
/// * `m1`, `m2` – masses \[kg\].
/// * `v1`, `v2` – pre-impact velocities \[m/s\].
/// * `e`        – coefficient of restitution \[0, 1\].
pub fn velocities_after_impact(m1: f64, m2: f64, v1: f64, v2: f64, e: f64) -> (f64, f64) {
    let total = m1 + m2;
    if total < 1e-30 {
        return (v1, v2);
    }
    let v1_new = (m1 * v1 + m2 * v2 + m2 * e * (v2 - v1)) / total;
    let v2_new = (m1 * v1 + m2 * v2 + m1 * e * (v1 - v2)) / total;
    (v1_new, v2_new)
}

/// Crash box simulation: vehicle-to-barrier collision model.
///
/// Uses a linear spring model for the crush zone.
#[derive(Debug, Clone)]
pub struct CrashBox {
    /// Vehicle mass \[kg\].
    pub vehicle_mass: f64,
    /// Crush zone stiffness (barrier stiffness) \[N/m\].
    pub barrier_stiffness: f64,
    /// Change in velocity (delta-V) \[m/s\].
    pub delta_v: f64,
    /// Maximum crush distance \[m\].
    pub crush_distance: f64,
}

impl CrashBox {
    /// Create a new `CrashBox` from mass, stiffness, and delta-V.
    ///
    /// The crush distance is computed as `delta_v * sqrt(m/k)`.
    pub fn new(mass: f64, k: f64, delta_v: f64) -> Self {
        let crush = if k > 0.0 {
            delta_v * (mass / k).sqrt()
        } else {
            0.0
        };
        Self {
            vehicle_mass: mass,
            barrier_stiffness: k,
            delta_v,
            crush_distance: crush,
        }
    }

    /// Estimate stopping time \[s\] = pi * sqrt(m/k) / 2.
    pub fn stopping_time(&self) -> f64 {
        if self.barrier_stiffness <= 0.0 {
            return 0.0;
        }
        PI / 2.0 * (self.vehicle_mass / self.barrier_stiffness).sqrt()
    }

    /// Mean deceleration \[m/s²\] = delta_v / stopping_time.
    pub fn mean_deceleration(&self) -> f64 {
        let t = self.stopping_time();
        if t <= 0.0 {
            return 0.0;
        }
        self.delta_v / t
    }

    /// Peak impact force \[N\] = k * crush_distance.
    pub fn peak_force(&self) -> f64 {
        self.barrier_stiffness * self.crush_distance
    }

    /// Energy absorbed by the crash box \[J\] = 0.5 * m * delta_v^2.
    pub fn energy_absorbed(&self) -> f64 {
        0.5 * self.vehicle_mass * self.delta_v * self.delta_v
    }
}

/// Check the HHC (Human Head/Chest) deflection criterion.
///
/// Returns `true` (safe) if `chest_deflection_mm` < 76 mm.
pub fn hhc_criterion(chest_deflection_mm: f64) -> bool {
    chest_deflection_mm < 76.0
}

/// Compute the Head Injury Criterion (HIC) from a time-series of accelerations.
///
/// HIC = max over all (t1, t2) pairs of `(t2 - t1) * (mean_a)^2.5`
/// where `mean_a` is the mean absolute acceleration over the interval.
///
/// * `a_t` – slice of acceleration values \[g or m/s²\] sampled at interval `dt`.
/// * `dt`  – time step \[s\].
///
/// Returns 0 if the slice has fewer than 2 samples.
pub fn head_injury_criterion(a_t: &[f64], dt: f64) -> f64 {
    let n = a_t.len();
    if n < 2 {
        return 0.0;
    }
    let mut hic_max = 0.0_f64;
    for i in 0..n {
        let mut sum = 0.0_f64;
        for j in (i + 1)..n {
            sum += a_t[j];
            let count = (j - i) as f64;
            let t_interval = count * dt;
            let mean_a = (sum / count).abs();
            let hic = t_interval * mean_a.powf(2.5);
            if hic > hic_max {
                hic_max = hic;
            }
        }
    }
    hic_max
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ImpactEvent ──────────────────────────────────────────────────────────

    #[test]
    fn impact_event_new_fields_set() {
        let ev = ImpactEvent::new(
            0,
            1,
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -2.0, 0.0],
            0.8,
            0.3,
        );
        assert_eq!(ev.body_a, 0);
        assert_eq!(ev.body_b, 1);
        assert!((ev.restitution - 0.8).abs() < 1e-12);
        assert!((ev.friction - 0.3).abs() < 1e-12);
    }

    #[test]
    fn impact_event_closing_speed_approaching() {
        // normal = [0,1,0], relative_velocity = [0,-3,0]  →  closing = 3
        let ev = ImpactEvent::new(
            0,
            1,
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -3.0, 0.0],
            1.0,
            0.0,
        );
        assert!((ev.closing_speed() - 3.0).abs() < 1e-12);
    }

    #[test]
    fn impact_event_closing_speed_separating() {
        let ev = ImpactEvent::new(
            0,
            1,
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 3.0, 0.0],
            1.0,
            0.0,
        );
        // separating → negative closing speed
        assert!(ev.closing_speed() < 0.0);
    }

    // ── HertzImpact ──────────────────────────────────────────────────────────

    #[test]
    fn hertz_force_zero_for_no_penetration() {
        let h = HertzImpact::new(1.0e9, 0.01);
        assert_eq!(h.compute_force(0.0), 0.0);
        assert_eq!(h.compute_force(-0.001), 0.0);
    }

    #[test]
    fn hertz_force_positive_for_penetration() {
        let h = HertzImpact::new(1.0e9, 0.01);
        let f = h.compute_force(0.001);
        assert!(f > 0.0, "force={f}");
    }

    #[test]
    fn hertz_force_increases_with_penetration() {
        let h = HertzImpact::new(1.0e9, 0.01);
        let f1 = h.compute_force(0.001);
        let f2 = h.compute_force(0.002);
        assert!(f2 > f1, "force should increase with penetration");
    }

    #[test]
    fn hertz_force_formula_check() {
        // F = (4/3) * E* * sqrt(R*) * delta^1.5
        let e_star = 1.0e9_f64;
        let r_star = 0.01_f64;
        let delta = 0.001_f64;
        let expected = (4.0 / 3.0) * e_star * r_star.sqrt() * delta.powf(1.5);
        let h = HertzImpact::new(e_star, r_star);
        let f = h.compute_force(delta);
        assert!(
            (f - expected).abs() / expected < 1e-12,
            "f={f}, expected={expected}"
        );
    }

    #[test]
    fn hertz_contact_duration_positive() {
        let h = HertzImpact::new(1.0e9, 0.01);
        let t = h.compute_contact_duration(0.1, 1.0);
        assert!(t > 0.0, "contact duration={t}");
    }

    #[test]
    fn hertz_contact_duration_zero_for_zero_velocity() {
        let h = HertzImpact::new(1.0e9, 0.01);
        assert_eq!(h.compute_contact_duration(1.0, 0.0), 0.0);
    }

    #[test]
    fn hertz_contact_duration_decreases_with_velocity() {
        let h = HertzImpact::new(1.0e9, 0.01);
        let t1 = h.compute_contact_duration(0.1, 1.0);
        let t2 = h.compute_contact_duration(0.1, 10.0);
        assert!(t2 < t1, "faster impact → shorter contact time");
    }

    // ── BallisticImpact ──────────────────────────────────────────────────────

    #[test]
    fn ballistic_speed_computed_correctly() {
        let b = BallisticImpact::new(0.01, [3.0, 4.0, 0.0], 1.0e9);
        assert!((b.speed() - 5.0).abs() < 1e-10, "speed={}", b.speed());
    }

    #[test]
    fn ballistic_penetration_capped_by_thickness() {
        let b = BallisticImpact::new(1.0, [1000.0, 0.0, 0.0], 1.0e6);
        let depth = b.penetration_depth(0.01);
        assert!(depth <= 0.01 + 1e-12, "depth={depth}");
    }

    #[test]
    fn ballistic_penetration_zero_velocity_is_zero() {
        let b = BallisticImpact::new(1.0, [0.0, 0.0, 0.0], 1.0e9);
        assert_eq!(b.penetration_depth(0.1), 0.0);
    }

    #[test]
    fn ballistic_residual_velocity_stopped() {
        // Very hard target, slow projectile → stopped
        let b = BallisticImpact::new(0.001, [10.0, 0.0, 0.0], 1.0e12);
        assert_eq!(b.residual_velocity(1.0), 0.0);
    }

    #[test]
    fn ballistic_residual_velocity_positive_for_thin_plate() {
        // High velocity, thin plate, low hardness → perforates
        let b = BallisticImpact::new(0.01, [1000.0, 0.0, 0.0], 1.0e6);
        let v_res = b.residual_velocity(0.001);
        assert!(v_res > 0.0, "v_res={v_res}");
    }

    // ── WaveImpact ───────────────────────────────────────────────────────────

    #[test]
    fn wave_reflection_free_surface() {
        // Z2 = 0 → reflection = -1
        let r = WaveImpact::reflection_coeff(1.0e6, 0.0);
        assert!((r + 1.0).abs() < 1e-10, "r={r}");
    }

    #[test]
    fn wave_reflection_rigid_wall() {
        // Z2 >> Z1 → reflection ≈ +1
        let r = WaveImpact::reflection_coeff(1.0e6, 1.0e12);
        assert!(r > 0.99, "r={r}");
    }

    #[test]
    fn wave_reflection_equal_impedance() {
        let r = WaveImpact::reflection_coeff(1.0e6, 1.0e6);
        assert!(r.abs() < 1e-10, "r={r}");
    }

    #[test]
    fn wave_transmission_energy_conservation() {
        let z1 = 2.0e6_f64;
        let z2 = 8.0e6_f64;
        let r = WaveImpact::reflection_coeff(z1, z2);
        let t = WaveImpact::transmission_coeff(z1, z2);
        // Energy conservation (power): R^2 + T^2 * (Z2/Z1) = 1
        let check = r * r + t * t * (z2 / z1) - 1.0;
        assert!(check.abs() < 1e-10, "check={check}, R={r}, T={t}");
    }

    #[test]
    fn wave_transmission_equal_impedance() {
        let t = WaveImpact::transmission_coeff(1.0e6, 1.0e6);
        assert!((t - 1.0).abs() < 1e-10, "t={t}");
    }

    #[test]
    fn wave_arrival_time_correct() {
        let w = WaveImpact::new(5000.0, 1.0e7);
        let t = w.arrival_time(10.0);
        assert!((t - 0.002).abs() < 1e-10, "t={t}");
    }

    #[test]
    fn wave_peak_stress_proportional_to_velocity() {
        let w = WaveImpact::new(5000.0, 1.0e7);
        let s1 = w.peak_stress(1.0);
        let s2 = w.peak_stress(2.0);
        assert!((s2 - 2.0 * s1).abs() < 1e-6, "s1={s1} s2={s2}");
    }

    // ── CrashSimulation ──────────────────────────────────────────────────────

    #[test]
    fn crash_sea_computed_correctly() {
        let crash = CrashSimulation::new(50_000.0, 40_000.0, 0.3, 133_333.0, 20.0);
        assert!((crash.specific_energy_absorption() - 2000.0).abs() < 1e-6);
    }

    #[test]
    fn crash_sea_zero_mass_returns_zero() {
        let crash = CrashSimulation::new(50_000.0, 40_000.0, 0.3, 133_333.0, 0.0);
        assert_eq!(crash.specific_energy_absorption(), 0.0);
    }

    #[test]
    fn crash_energy_absorption_ratio_clamped() {
        let crash = CrashSimulation::new(10_000.0, 15_000.0, 0.2, 50_000.0, 5.0);
        assert!(crash.energy_absorption_ratio() <= 1.0);
    }

    #[test]
    fn crash_energy_absorption_ratio_correct() {
        let crash = CrashSimulation::new(100_000.0, 80_000.0, 0.3, 266_666.0, 50.0);
        assert!((crash.energy_absorption_ratio() - 0.8).abs() < 1e-10);
    }

    #[test]
    fn crash_mean_deceleration_positive() {
        let crash = CrashSimulation::new(100_000.0, 90_000.0, 0.5, 180_000.0, 50.0);
        let decel = crash.mean_deceleration(1500.0, 10.0);
        assert!(decel > 0.0, "decel={decel}");
    }

    // ── coefficient_of_restitution ────────────────────────────────────────────

    #[test]
    fn cor_elastic_collision() {
        let cor = coefficient_of_restitution(5.0, 5.0);
        assert!((cor - 1.0).abs() < 1e-12, "cor={cor}");
    }

    #[test]
    fn cor_plastic_collision() {
        let cor = coefficient_of_restitution(5.0, 0.0);
        assert!(cor.abs() < 1e-12, "cor={cor}");
    }

    #[test]
    fn cor_partial_restitution() {
        let cor = coefficient_of_restitution(5.0, 3.0);
        assert!((cor - 0.6).abs() < 1e-12, "cor={cor}");
    }

    #[test]
    fn cor_zero_before_is_zero() {
        let cor = coefficient_of_restitution(0.0, 1.0);
        assert_eq!(cor, 0.0);
    }

    // ── impulse ──────────────────────────────────────────────────────────────

    #[test]
    fn impulse_basic() {
        let imp = impulse(2.0, [3.0, 0.0, 0.0]);
        assert!((imp[0] - 6.0).abs() < 1e-12);
        assert!(imp[1].abs() < 1e-12);
        assert!(imp[2].abs() < 1e-12);
    }

    #[test]
    fn impulse_zero_mass() {
        let imp = impulse(0.0, [10.0, 5.0, 2.0]);
        assert!(imp[0].abs() < 1e-12);
    }

    #[test]
    fn impulse_vector_proportional_to_mass() {
        let imp1 = impulse(1.0, [1.0, 2.0, 3.0]);
        let imp2 = impulse(3.0, [1.0, 2.0, 3.0]);
        for i in 0..3 {
            assert!((imp2[i] - 3.0 * imp1[i]).abs() < 1e-12);
        }
    }

    // ── johnson_impact_number ─────────────────────────────────────────────────

    #[test]
    fn johnson_impact_elastic_regime() {
        // Low velocity → J << 1 (elastic)
        let j = johnson_impact_number(7800.0, 1.0, 250.0e6);
        assert!(j < 1.0, "j={j}");
    }

    #[test]
    fn johnson_impact_plastic_regime() {
        // High velocity → J >> 1 (plastic)
        let j = johnson_impact_number(7800.0, 1000.0, 250.0e6);
        assert!(j > 1.0, "j={j}");
    }

    #[test]
    fn johnson_impact_zero_sigma() {
        let j = johnson_impact_number(7800.0, 100.0, 0.0);
        assert!(j.is_infinite());
    }

    // ── central_impact_velocities ────────────────────────────────────────────

    #[test]
    fn central_impact_elastic_equal_mass_exchange() {
        // Equal mass, elastic → velocities exchange
        let (va, vb) = central_impact_velocities(1.0, 1.0, 3.0, -1.0, 1.0);
        assert!((va + 1.0).abs() < 1e-10, "va={va}");
        assert!((vb - 3.0).abs() < 1e-10, "vb={vb}");
    }

    #[test]
    fn central_impact_plastic_equal_mass() {
        // Equal mass, perfectly plastic → common velocity
        let (va, vb) = central_impact_velocities(1.0, 1.0, 4.0, 0.0, 0.0);
        assert!((va - 2.0).abs() < 1e-10, "va={va}");
        assert!((vb - 2.0).abs() < 1e-10, "vb={vb}");
    }

    #[test]
    fn central_impact_momentum_conserved() {
        let ma = 3.0_f64;
        let mb = 1.0_f64;
        let va_init = 2.0_f64;
        let vb_init = -1.0_f64;
        let e = 0.5_f64;
        let (va_new, vb_new) = central_impact_velocities(ma, mb, va_init, vb_init, e);
        let p_before = ma * va_init + mb * vb_init;
        let p_after = ma * va_new + mb * vb_new;
        assert!((p_before - p_after).abs() < 1e-10, "momentum not conserved");
    }

    // ── hertz_max_penetration ────────────────────────────────────────────────

    #[test]
    fn hertz_max_penetration_positive() {
        let delta = hertz_max_penetration(1.0e9, 0.01, 0.1, 1.0);
        assert!(delta > 0.0, "delta={delta}");
    }

    #[test]
    fn hertz_max_penetration_increases_with_velocity() {
        let d1 = hertz_max_penetration(1.0e9, 0.01, 0.1, 1.0);
        let d2 = hertz_max_penetration(1.0e9, 0.01, 0.1, 10.0);
        assert!(d2 > d1, "d1={d1} d2={d2}");
    }

    // ── hertz contact radius / pressure ──────────────────────────────────────

    #[test]
    fn hertz_contact_radius_zero_penetration() {
        assert_eq!(hertz_contact_radius(0.01, 0.0), 0.0);
    }

    #[test]
    fn hertz_contact_radius_positive() {
        let a = hertz_contact_radius(0.01, 0.001);
        assert!(a > 0.0, "a={a}");
    }

    #[test]
    fn hertz_max_pressure_zero_area() {
        assert_eq!(hertz_max_pressure(100.0, 0.0), 0.0);
    }

    #[test]
    fn hertz_max_pressure_formula() {
        let force = 1000.0_f64;
        let a = 0.01_f64;
        let expected = 1.5 * force / (PI * a * a);
        let p = hertz_max_pressure(force, a);
        assert!((p - expected).abs() / expected < 1e-12, "p={p}");
    }

    // ── coefficient_of_restitution_estimate ──────────────────────────────────

    #[test]
    fn cor_estimate_steel() {
        assert!((coefficient_of_restitution_estimate("steel") - 0.7).abs() < 1e-12);
    }

    #[test]
    fn cor_estimate_rubber() {
        assert!((coefficient_of_restitution_estimate("rubber") - 0.8).abs() < 1e-12);
    }

    #[test]
    fn cor_estimate_clay() {
        assert!((coefficient_of_restitution_estimate("clay") - 0.1).abs() < 1e-12);
    }

    #[test]
    fn cor_estimate_glass() {
        assert!((coefficient_of_restitution_estimate("glass") - 0.65).abs() < 1e-12);
    }

    #[test]
    fn cor_estimate_in_range() {
        for mat in &["steel", "rubber", "clay", "glass"] {
            let e = coefficient_of_restitution_estimate(mat);
            assert!((0.0..=1.0).contains(&e), "e={e} out of range for {mat}");
        }
    }

    // ── hertz_impact_force ────────────────────────────────────────────────────

    #[test]
    fn hertz_impact_force_zero_delta() {
        assert_eq!(hertz_impact_force(1.0e9, 0.0, 1.5), 0.0);
    }

    #[test]
    fn hertz_impact_force_positive_delta() {
        let f = hertz_impact_force(1.0e9, 0.001, 1.5);
        assert!(f > 0.0, "f={f}");
    }

    #[test]
    fn hertz_impact_force_increases_with_delta() {
        let f1 = hertz_impact_force(1.0e9, 0.001, 1.5);
        let f2 = hertz_impact_force(1.0e9, 0.002, 1.5);
        assert!(f2 > f1);
    }

    #[test]
    fn hertz_impact_force_formula_check() {
        let k = 1.0e8_f64;
        let delta = 0.005_f64;
        let exp = 1.5_f64;
        let expected = k * delta.powf(exp);
        assert!((hertz_impact_force(k, delta, exp) - expected).abs() < 1e-6 * expected);
    }

    // ── hertz_contact_stiffness ───────────────────────────────────────────────

    #[test]
    fn hertz_contact_stiffness_positive() {
        let k = hertz_contact_stiffness(2.0e11, 0.3, 2.0e11, 0.3, 0.05, 0.05);
        assert!(k > 0.0, "k={k}");
    }

    #[test]
    fn hertz_contact_stiffness_increases_with_modulus() {
        let k1 = hertz_contact_stiffness(1.0e11, 0.3, 1.0e11, 0.3, 0.05, 0.05);
        let k2 = hertz_contact_stiffness(2.0e11, 0.3, 2.0e11, 0.3, 0.05, 0.05);
        assert!(k2 > k1);
    }

    // ── impact_duration_hertz ─────────────────────────────────────────────────

    #[test]
    fn impact_duration_hertz_positive() {
        let t = impact_duration_hertz(0.1, 1.0e6, 1.0);
        assert!(t > 0.0, "t={t}");
    }

    #[test]
    fn impact_duration_hertz_zero_velocity() {
        assert_eq!(impact_duration_hertz(0.1, 1.0e6, 0.0), 0.0);
    }

    #[test]
    fn impact_duration_hertz_decreases_with_velocity() {
        let t1 = impact_duration_hertz(0.1, 1.0e6, 1.0);
        let t2 = impact_duration_hertz(0.1, 1.0e6, 10.0);
        assert!(t2 < t1, "t1={t1} t2={t2}");
    }

    // ── peak_impact_force ─────────────────────────────────────────────────────

    #[test]
    fn peak_impact_force_positive() {
        let f = peak_impact_force(1.0, 1.0, 1.0e6);
        assert!(f > 0.0, "f={f}");
    }

    #[test]
    fn peak_impact_force_increases_with_velocity() {
        let f1 = peak_impact_force(1.0, 1.0, 1.0e6);
        let f2 = peak_impact_force(1.0, 10.0, 1.0e6);
        assert!(f2 > f1);
    }

    // ── impulse_from_restitution ──────────────────────────────────────────────

    #[test]
    fn impulse_from_restitution_elastic() {
        // Equal masses, v1=1, v2=0, e=1: j = -2*(1-0)/(1+1) = -1
        let j = impulse_from_restitution(1.0, 1.0, 1.0, 0.0, 1.0);
        assert!((j - (-1.0)).abs() < 1e-10, "j={j}");
    }

    #[test]
    fn impulse_from_restitution_plastic() {
        // e=0: j = -(v1-v2)/(1/m1+1/m2)
        let j = impulse_from_restitution(1.0, 1.0, 4.0, 0.0, 0.0);
        assert!((j - (-2.0)).abs() < 1e-10, "j={j}");
    }

    // ── velocities_after_impact ───────────────────────────────────────────────

    #[test]
    fn velocities_after_impact_elastic_equal_mass() {
        let (v1, v2) = velocities_after_impact(1.0, 1.0, 3.0, 0.0, 1.0);
        // Velocities exchange for equal mass elastic
        assert!((v1 - 0.0).abs() < 1e-10, "v1={v1}");
        assert!((v2 - 3.0).abs() < 1e-10, "v2={v2}");
    }

    #[test]
    fn velocities_after_impact_plastic_equal_mass() {
        let (v1, v2) = velocities_after_impact(1.0, 1.0, 4.0, 0.0, 0.0);
        assert!((v1 - 2.0).abs() < 1e-10, "v1={v1}");
        assert!((v2 - 2.0).abs() < 1e-10, "v2={v2}");
    }

    #[test]
    fn velocities_after_impact_conserves_momentum() {
        let m1 = 3.0_f64;
        let m2 = 1.0_f64;
        let v1i = 2.0_f64;
        let v2i = -1.0_f64;
        let e = 0.6_f64;
        let (v1f, v2f) = velocities_after_impact(m1, m2, v1i, v2i, e);
        let p_before = m1 * v1i + m2 * v2i;
        let p_after = m1 * v1f + m2 * v2f;
        assert!((p_before - p_after).abs() < 1e-9, "momentum not conserved");
    }

    #[test]
    fn velocities_after_impact_elastic_conserves_energy() {
        let m1 = 2.0_f64;
        let m2 = 3.0_f64;
        let v1i = 5.0_f64;
        let v2i = -1.0_f64;
        let (v1f, v2f) = velocities_after_impact(m1, m2, v1i, v2i, 1.0);
        let ke_before = 0.5 * m1 * v1i * v1i + 0.5 * m2 * v2i * v2i;
        let ke_after = 0.5 * m1 * v1f * v1f + 0.5 * m2 * v2f * v2f;
        assert!((ke_before - ke_after).abs() < 1e-6, "energy not conserved");
    }

    // ── CrashBox ──────────────────────────────────────────────────────────────

    #[test]
    fn crash_box_energy_absorbed_formula() {
        let c = CrashBox::new(1500.0, 5.0e5, 10.0);
        let expected = 0.5 * 1500.0 * 100.0;
        assert!((c.energy_absorbed() - expected).abs() < 1e-6 * expected);
    }

    #[test]
    fn crash_box_stopping_time_positive() {
        let c = CrashBox::new(1500.0, 5.0e5, 10.0);
        assert!(c.stopping_time() > 0.0);
    }

    #[test]
    fn crash_box_mean_deceleration_positive() {
        let c = CrashBox::new(1500.0, 5.0e5, 10.0);
        assert!(c.mean_deceleration() > 0.0);
    }

    #[test]
    fn crash_box_peak_force_positive() {
        let c = CrashBox::new(1500.0, 5.0e5, 10.0);
        assert!(c.peak_force() > 0.0);
    }

    // ── hhc_criterion ─────────────────────────────────────────────────────────

    #[test]
    fn hhc_safe_below_76mm() {
        assert!(hhc_criterion(50.0));
    }

    #[test]
    fn hhc_unsafe_at_76mm() {
        assert!(!hhc_criterion(76.0));
    }

    #[test]
    fn hhc_unsafe_above_76mm() {
        assert!(!hhc_criterion(100.0));
    }

    // ── head_injury_criterion ─────────────────────────────────────────────────

    #[test]
    fn hic_empty_returns_zero() {
        assert_eq!(head_injury_criterion(&[], 0.001), 0.0);
    }

    #[test]
    fn hic_single_sample_returns_zero() {
        assert_eq!(head_injury_criterion(&[100.0], 0.001), 0.0);
    }

    #[test]
    fn hic_uniform_acceleration_positive() {
        let a_t: Vec<f64> = vec![100.0; 100];
        let hic = head_injury_criterion(&a_t, 0.001);
        assert!(hic > 0.0, "hic={hic}");
    }

    #[test]
    fn hic_zero_acceleration_is_zero() {
        let a_t = vec![0.0_f64; 50];
        let hic = head_injury_criterion(&a_t, 0.001);
        assert_eq!(hic, 0.0);
    }
}
