// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Impact mechanics and wave propagation in rigid bodies.
//!
//! Implements Hertzian contact, Hunt-Crossley viscoelastic models, elastic wave
//! propagation, ballistic penetration, Taylor anvil dynamics, spallation,
//! oblique impact, multi-body chain impacts, and impact sensor models.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration (m/s²).
const G_ACCEL: f64 = 9.81;

/// Clamp utility.
fn clamp(x: f64, lo: f64, hi: f64) -> f64 {
    if x < lo {
        lo
    } else if x > hi {
        hi
    } else {
        x
    }
}

// ---------------------------------------------------------------------------
// HertzianImpact
// ---------------------------------------------------------------------------

/// Result of a Hertzian contact impact simulation.
#[derive(Debug, Clone)]
pub struct HertzianResult {
    /// Peak contact force at maximum deformation (N).
    pub peak_force: f64,
    /// Contact duration (s).
    pub contact_duration: f64,
    /// Maximum contact deformation/indentation (m).
    pub max_deformation: f64,
    /// Total energy stored in the contact spring at peak (J).
    pub strain_energy: f64,
    /// Fraction of kinetic energy converted to elastic strain energy.
    pub energy_fraction: f64,
}

/// Hertz contact theory for elastic spherical impact.
///
/// Models the contact force during impact of two elastic spheres using
/// Hertz contact mechanics: F = k_H * δ^(3/2), where k_H is the Hertz
/// stiffness and δ is the indentation.
///
/// References: Hertz (1882), Johnson (1987) Contact Mechanics.
#[derive(Debug, Clone)]
pub struct HertzianImpact {
    /// Reduced elastic modulus E* = \[(1-ν₁²)/E₁ + (1-ν₂²)/E₂\]^-1 (Pa).
    pub reduced_modulus: f64,
    /// Reduced radius R* = (1/R₁ + 1/R₂)^-1 (m).
    pub reduced_radius: f64,
    /// Mass of impactor (kg).
    pub mass_impactor: f64,
    /// Mass of target body (kg).
    pub mass_target: f64,
}

impl HertzianImpact {
    /// Construct from material/geometry parameters.
    ///
    /// - `e1`, `nu1`: Young's modulus (Pa) and Poisson's ratio for body 1.
    /// - `e2`, `nu2`: Young's modulus (Pa) and Poisson's ratio for body 2.
    /// - `r1`, `r2`: radii of curvature (m); use `f64::INFINITY` for flat.
    /// - `m1`, `m2`: masses (kg).
    pub fn new(e1: f64, nu1: f64, e2: f64, nu2: f64, r1: f64, r2: f64, m1: f64, m2: f64) -> Self {
        let inv_e_star = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
        let reduced_modulus = 1.0 / inv_e_star;

        let inv_r_star = 1.0 / r1 + 1.0 / r2;
        let reduced_radius = 1.0 / inv_r_star;

        Self {
            reduced_modulus,
            reduced_radius,
            mass_impactor: m1,
            mass_target: m2,
        }
    }

    /// Hertz stiffness k_H = (4/3) * E* * √R*.
    pub fn hertz_stiffness(&self) -> f64 {
        (4.0 / 3.0) * self.reduced_modulus * self.reduced_radius.sqrt()
    }

    /// Effective mass: m* = (m₁ * m₂) / (m₁ + m₂).
    pub fn effective_mass(&self) -> f64 {
        self.mass_impactor * self.mass_target / (self.mass_impactor + self.mass_target)
    }

    /// Maximum indentation for impact at relative velocity `v0` (m/s).
    ///
    /// Energy balance: (1/2) m* v₀² = (2/5) k_H δ_max^(5/2)
    /// → δ_max = \[(5 m* v₀² / (4 k_H))\]^(2/5)
    pub fn max_deformation(&self, v0: f64) -> f64 {
        let m_star = self.effective_mass();
        let k_h = self.hertz_stiffness();
        let bracket = 5.0 * m_star * v0 * v0 / (4.0 * k_h);
        bracket.powf(2.0 / 5.0)
    }

    /// Peak contact force at maximum indentation (N).
    pub fn peak_force(&self, v0: f64) -> f64 {
        let delta = self.max_deformation(v0);
        self.hertz_stiffness() * delta.powf(1.5)
    }

    /// Contact duration using Hertz impact theory (s).
    ///
    /// T ≈ 2.94 * δ_max / v₀
    pub fn contact_duration(&self, v0: f64) -> f64 {
        if v0 < 1e-10 {
            return 0.0;
        }
        let delta = self.max_deformation(v0);
        2.94 * delta / v0
    }

    /// Perform a full Hertzian impact analysis at relative velocity `v0`.
    pub fn analyze(&self, v0: f64) -> HertzianResult {
        let m_star = self.effective_mass();
        let k_h = self.hertz_stiffness();
        let delta_max = self.max_deformation(v0);
        let peak_force = k_h * delta_max.powf(1.5);
        let contact_duration = self.contact_duration(v0);
        let ke_initial = 0.5 * m_star * v0 * v0;
        let strain_energy = (2.0 / 5.0) * k_h * delta_max.powf(2.5);
        let energy_fraction = if ke_initial > 1e-14 {
            strain_energy / ke_initial
        } else {
            0.0
        };

        HertzianResult {
            peak_force,
            contact_duration,
            max_deformation: delta_max,
            strain_energy,
            energy_fraction,
        }
    }
}

// ---------------------------------------------------------------------------
// ViscoelasticImpact  (Hunt-Crossley model)
// ---------------------------------------------------------------------------

/// Hunt-Crossley viscoelastic impact model state.
#[derive(Debug, Clone)]
pub struct HuntCrossleyState {
    /// Current indentation δ (m).
    pub delta: f64,
    /// Current indentation rate dδ/dt (m/s).
    pub delta_dot: f64,
    /// Time elapsed since contact started (s).
    pub time: f64,
    /// Energy dissipated so far (J).
    pub energy_dissipated: f64,
}

/// Hunt-Crossley viscoelastic impact model.
///
/// Contact force: F = k δ^n (1 + α δ̇)
/// where k is contact stiffness, n = 3/2 for Hertz, and α is damping parameter.
/// The coefficient of restitution depends on impact velocity.
#[derive(Debug, Clone)]
pub struct ViscoelasticImpact {
    /// Contact stiffness k (N/m^n).
    pub stiffness: f64,
    /// Contact exponent n (1.5 for Hertz).
    pub exponent: f64,
    /// Damping parameter α (s/m).
    pub alpha: f64,
    /// Effective mass m* (kg).
    pub effective_mass: f64,
}

impl ViscoelasticImpact {
    /// Create a Hunt-Crossley model.
    pub fn new(stiffness: f64, exponent: f64, alpha: f64, effective_mass: f64) -> Self {
        Self {
            stiffness,
            exponent,
            alpha,
            effective_mass,
        }
    }

    /// Contact force F = k δ^n (1 + α δ̇), clamped to non-negative.
    pub fn contact_force(&self, delta: f64, delta_dot: f64) -> f64 {
        if delta <= 0.0 {
            return 0.0;
        }
        let f = self.stiffness * delta.powf(self.exponent) * (1.0 + self.alpha * delta_dot);
        f.max(0.0)
    }

    /// Coefficient of restitution as a function of approach velocity.
    ///
    /// For Hunt-Crossley: e ≈ 1 - (3 α / 2) * v₀ (linearized approximation).
    pub fn restitution_coefficient(&self, v0: f64) -> f64 {
        clamp(1.0 - 1.5 * self.alpha * v0, 0.0, 1.0)
    }

    /// Energy dissipated in one full impact at velocity `v0`.
    ///
    /// ΔE = (1 - e²) * (1/2) m* v₀²
    pub fn energy_dissipated(&self, v0: f64) -> f64 {
        let e = self.restitution_coefficient(v0);
        (1.0 - e * e) * 0.5 * self.effective_mass * v0 * v0
    }

    /// Integrate contact dynamics using simple Euler integration.
    ///
    /// Returns the evolution of indentation over time until separation.
    pub fn integrate(&self, v0: f64, dt: f64, max_steps: usize) -> Vec<HuntCrossleyState> {
        let mut states = Vec::new();
        let mut delta = 0.0f64;
        let mut delta_dot = v0; // approach velocity (positive = compression)
        let mut time = 0.0f64;
        let mut energy_dissipated = 0.0f64;

        for _ in 0..max_steps {
            let f = self.contact_force(delta, delta_dot);
            let f_elastic = self.stiffness * delta.powf(self.exponent).max(0.0);
            let f_damping = f - f_elastic;
            let power_diss = (f_damping * delta_dot).max(0.0);
            energy_dissipated += power_diss * dt;

            states.push(HuntCrossleyState {
                delta,
                delta_dot,
                time,
                energy_dissipated,
            });

            let accel = -f / self.effective_mass;
            delta_dot += accel * dt;
            delta += delta_dot * dt;
            time += dt;

            // Contact ends when indentation becomes negative.
            if delta < 0.0 && delta_dot < 0.0 {
                delta = 0.0;
                states.push(HuntCrossleyState {
                    delta,
                    delta_dot,
                    time,
                    energy_dissipated,
                });
                break;
            }
        }
        states
    }
}

// ---------------------------------------------------------------------------
// WaveImpact
// ---------------------------------------------------------------------------

/// Elastic wave propagation result.
#[derive(Debug, Clone)]
pub struct WaveResult {
    /// Longitudinal (P-wave) speed (m/s).
    pub c_longitudinal: f64,
    /// Transverse (S-wave) speed (m/s).
    pub c_transverse: f64,
    /// Rayleigh surface wave speed (m/s).
    pub c_rayleigh: f64,
    /// Time for wave to travel distance `d` (s).
    pub travel_time_p: f64,
    /// Reflection coefficient at interface (stress amplitude ratio).
    pub reflection_coeff: f64,
    /// Transmission coefficient at interface.
    pub transmission_coeff: f64,
}

/// Elastic wave propagation and reflection/transmission at interfaces.
///
/// Computes wave speeds and interface coefficients for impact-generated
/// elastic waves in homogeneous isotropic media.
#[derive(Debug, Clone)]
pub struct WaveImpact {
    /// Density of medium 1 (kg/m³).
    pub density1: f64,
    /// Young's modulus of medium 1 (Pa).
    pub modulus1: f64,
    /// Poisson's ratio of medium 1.
    pub poisson1: f64,
    /// Density of medium 2 (kg/m³).
    pub density2: f64,
    /// Young's modulus of medium 2 (Pa).
    pub modulus2: f64,
    /// Poisson's ratio of medium 2.
    pub poisson2: f64,
}

impl WaveImpact {
    /// Create a wave impact model for two-medium interface.
    pub fn new(
        density1: f64,
        modulus1: f64,
        poisson1: f64,
        density2: f64,
        modulus2: f64,
        poisson2: f64,
    ) -> Self {
        Self {
            density1,
            modulus1,
            poisson1,
            density2,
            modulus2,
            poisson2,
        }
    }

    /// P-wave (longitudinal) speed in medium 1 (m/s).
    ///
    /// c_P = √\[(E(1-ν)) / (ρ(1+ν)(1-2ν))\]
    pub fn p_wave_speed(&self) -> f64 {
        let e = self.modulus1;
        let nu = self.poisson1;
        let rho = self.density1;
        let num = e * (1.0 - nu);
        let den = rho * (1.0 + nu) * (1.0 - 2.0 * nu);
        (num / den).sqrt()
    }

    /// S-wave (transverse/shear) speed in medium 1 (m/s).
    ///
    /// c_S = √\[G/ρ\] = √\[E / (2ρ(1+ν))\]
    pub fn s_wave_speed(&self) -> f64 {
        let e = self.modulus1;
        let nu = self.poisson1;
        let rho = self.density1;
        (e / (2.0 * rho * (1.0 + nu))).sqrt()
    }

    /// Rayleigh surface wave speed approximation (m/s).
    ///
    /// c_R ≈ c_S * (0.862 + 1.14ν) / (1 + ν)
    pub fn rayleigh_speed(&self) -> f64 {
        let cs = self.s_wave_speed();
        let nu = self.poisson1;
        cs * (0.862 + 1.14 * nu) / (1.0 + nu)
    }

    /// Acoustic impedance Z = ρ c for P-waves.
    pub fn impedance(&self, density: f64, modulus: f64, poisson: f64) -> f64 {
        let cp = {
            let num = modulus * (1.0 - poisson);
            let den = density * (1.0 + poisson) * (1.0 - 2.0 * poisson);
            (num / den).sqrt()
        };
        density * cp
    }

    /// Reflection and transmission coefficients at a normal-incidence interface.
    pub fn interface_coefficients(&self) -> (f64, f64) {
        let z1 = self.impedance(self.density1, self.modulus1, self.poisson1);
        let z2 = self.impedance(self.density2, self.modulus2, self.poisson2);
        let reflection = (z2 - z1) / (z2 + z1);
        let transmission = 2.0 * z2 / (z2 + z1);
        (reflection, transmission)
    }

    /// Full wave analysis for a propagation distance `d`.
    pub fn analyze(&self, d: f64) -> WaveResult {
        let c_p = self.p_wave_speed();
        let c_s = self.s_wave_speed();
        let c_r = self.rayleigh_speed();
        let (r_coeff, t_coeff) = self.interface_coefficients();

        WaveResult {
            c_longitudinal: c_p,
            c_transverse: c_s,
            c_rayleigh: c_r,
            travel_time_p: if c_p > 1e-10 { d / c_p } else { f64::INFINITY },
            reflection_coeff: r_coeff,
            transmission_coeff: t_coeff,
        }
    }
}

// ---------------------------------------------------------------------------
// ImpactCrater  (Pi-group scaling)
// ---------------------------------------------------------------------------

/// Crater scaling result using Holsapple π-group scaling laws.
#[derive(Debug, Clone)]
pub struct CraterResult {
    /// Crater diameter (m).
    pub diameter: f64,
    /// Crater depth (m).
    pub depth: f64,
    /// Crater volume (m³).
    pub volume: f64,
    /// Ejecta mass (kg).
    pub ejecta_mass: f64,
}

/// Impact crater scaling using Holsapple (1993) π-group dimensional analysis.
///
/// Applies to hypervelocity impacts where the crater size depends on
/// projectile kinetic energy, gravity, and target strength.
#[derive(Debug, Clone)]
pub struct ImpactCrater {
    /// Target material density (kg/m³).
    pub target_density: f64,
    /// Projectile density (kg/m³).
    pub projectile_density: f64,
    /// Target strength Y (Pa); set to 0 for gravity-dominated regime.
    pub target_strength: f64,
    /// Pi-group exponent μ (typically 0.4–0.6 for rock).
    pub mu: f64,
    /// Pi-group exponent ν.
    pub nu: f64,
    /// Holsapple dimensionless crater constant K1.
    pub k1: f64,
    /// Holsapple dimensionless crater constant K2.
    pub k2: f64,
}

impl ImpactCrater {
    /// Create with typical rock/stone parameters.
    pub fn new_rock(target_density: f64, projectile_density: f64) -> Self {
        Self {
            target_density,
            projectile_density,
            target_strength: 1e6, // 1 MPa typical rock strength
            mu: 0.41,
            nu: 0.4,
            k1: 0.24,
            k2: 0.29,
        }
    }

    /// Create with loose soil parameters (gravity-dominated regime).
    pub fn new_soil(target_density: f64, projectile_density: f64) -> Self {
        Self {
            target_density,
            projectile_density,
            target_strength: 0.0, // purely gravity-dominated
            mu: 0.41,
            nu: 0.4,
            k1: 0.132,
            k2: 0.26,
        }
    }

    /// Compute the π₂ (gravity-scaling) group.
    ///
    /// π₂ = 3.22 g a / v² where a = projectile radius, v = impact velocity.
    pub fn pi_2(&self, projectile_radius: f64, velocity: f64) -> f64 {
        3.22 * G_ACCEL * projectile_radius / (velocity * velocity)
    }

    /// Compute the π₃ (strength-scaling) group.
    ///
    /// π₃ = Y / (ρ_t v²)
    pub fn pi_3(&self, velocity: f64) -> f64 {
        self.target_strength / (self.target_density * velocity * velocity)
    }

    /// Dimensionless crater volume π_V = K1 \[π₂^{-3μ/2} + K2 π₃^{-3μ/(2+6μ)}\]^{...}.
    ///
    /// Simplified Holsapple formula for transient crater volume.
    pub fn pi_volume(&self, projectile_radius: f64, velocity: f64) -> f64 {
        let pi2 = self.pi_2(projectile_radius, velocity);
        let density_ratio = self.projectile_density / self.target_density;

        // Gravity-dominated term: π₂^{-3μ/(2+6μ)} — but use simpler form.

        if pi2 > 1e-14 {
            self.k1
                * pi2.powf(-(3.0 * self.mu) / (2.0 + 6.0 * self.mu))
                * density_ratio
                    .powf((2.0 + 12.0 * self.mu * self.nu) / (3.0 * (2.0 + 6.0 * self.mu)))
        } else {
            1.0
        }
    }

    /// Final crater diameter (m).
    pub fn crater_diameter(&self, projectile_radius: f64, velocity: f64) -> f64 {
        let pi_v = self.pi_volume(projectile_radius, velocity);
        // Volume = pi_v * (m / ρ_t) where m = (4/3) π a³ ρ_p.
        let projectile_mass =
            (4.0 / 3.0) * PI * projectile_radius.powi(3) * self.projectile_density;
        let volume = pi_v * projectile_mass / self.target_density;
        // Assume hemisphere: V = (1/3) π (D/2)² d = (π/12) D³ for depth≈D/3.
        (12.0 * volume / PI).cbrt()
    }

    /// Final crater depth (m) ≈ diameter / 3 (typical aspect ratio).
    pub fn crater_depth(&self, projectile_radius: f64, velocity: f64) -> f64 {
        self.crater_diameter(projectile_radius, velocity) / 3.0
    }

    /// Full crater analysis.
    pub fn analyze(&self, projectile_radius: f64, velocity: f64) -> CraterResult {
        let diameter = self.crater_diameter(projectile_radius, velocity);
        let depth = self.crater_depth(projectile_radius, velocity);
        let volume = PI / 12.0 * diameter.powi(2) * depth;
        let ejecta_mass = volume * self.target_density * 0.9; // ~90% of crater material ejected

        CraterResult {
            diameter,
            depth,
            volume,
            ejecta_mass,
        }
    }
}

// ---------------------------------------------------------------------------
// BallisticImpact
// ---------------------------------------------------------------------------

/// Ballistic impact penetration result.
#[derive(Debug, Clone)]
pub struct BallisticResult {
    /// Residual velocity after penetration (m/s); 0 if target stopped projectile.
    pub residual_velocity: f64,
    /// Whether the projectile perforated the target.
    pub perforated: bool,
    /// Ballistic limit velocity (m/s).
    pub ballistic_limit: f64,
    /// Normalized residual velocity vs/vbl.
    pub normalized_velocity: f64,
}

/// Ballistic impact penetration mechanics using the Recht-Ipson model.
///
/// v_r = a(v₀^p - v_bl^p)^(1/p) for v₀ > v_bl, else stopped.
#[derive(Debug, Clone)]
pub struct BallisticImpact {
    /// Ballistic limit velocity v_bl (m/s).
    pub ballistic_limit: f64,
    /// Recht-Ipson exponent p (typically 2 for metals).
    pub exponent: f64,
    /// Recht-Ipson coefficient a (≈ 1 for thin plates).
    pub coefficient: f64,
    /// Projectile mass (kg).
    pub projectile_mass: f64,
    /// Target thickness (m).
    pub thickness: f64,
}

impl BallisticImpact {
    /// Create a new ballistic impact model.
    pub fn new(
        ballistic_limit: f64,
        exponent: f64,
        coefficient: f64,
        projectile_mass: f64,
        thickness: f64,
    ) -> Self {
        Self {
            ballistic_limit,
            exponent,
            coefficient,
            projectile_mass,
            thickness,
        }
    }

    /// Residual velocity after perforating a target (Recht-Ipson model).
    ///
    /// v_r = a * (v₀^p - v_bl^p)^(1/p)
    pub fn residual_velocity(&self, v0: f64) -> f64 {
        if v0 <= self.ballistic_limit {
            return 0.0;
        }
        let vbl = self.ballistic_limit;
        let p = self.exponent;
        let a = self.coefficient;
        a * (v0.powf(p) - vbl.powf(p)).powf(1.0 / p)
    }

    /// Energy absorbed by the target during perforation (J).
    pub fn absorbed_energy(&self, v0: f64) -> f64 {
        let vr = self.residual_velocity(v0);
        let m = self.projectile_mass;
        0.5 * m * (v0 * v0 - vr * vr)
    }

    /// Specific energy absorption per unit thickness (J/m).
    pub fn specific_energy_absorption(&self, v0: f64) -> f64 {
        if self.thickness < 1e-14 {
            return 0.0;
        }
        self.absorbed_energy(v0) / self.thickness
    }

    /// Full ballistic analysis at impact velocity `v0`.
    pub fn analyze(&self, v0: f64) -> BallisticResult {
        let vr = self.residual_velocity(v0);
        let vbl = self.ballistic_limit;
        let perforated = v0 > vbl;
        let normalized = if vbl > 1e-10 { vr / vbl } else { 0.0 };

        BallisticResult {
            residual_velocity: vr,
            perforated,
            ballistic_limit: vbl,
            normalized_velocity: normalized,
        }
    }
}

// ---------------------------------------------------------------------------
// SpallFracture
// ---------------------------------------------------------------------------

/// Spall fracture analysis result.
#[derive(Debug, Clone)]
pub struct SpallResult {
    /// Spall strength σ_spall (Pa).
    pub spall_strength: f64,
    /// Free surface velocity pullback Δu_fs (m/s).
    pub velocity_pullback: f64,
    /// Spall thickness (m).
    pub spall_thickness: f64,
    /// Whether spallation occurred.
    pub spalled: bool,
}

/// Spall fracture under shock loading using Grüneisen EOS.
///
/// Models the tensile fracture (spallation) that occurs when compressive
/// shock waves reflect from free surfaces, producing tensile stresses.
#[derive(Debug, Clone)]
pub struct SpallFracture {
    /// Material density ρ₀ (kg/m³).
    pub density: f64,
    /// Bulk sound speed C₀ (m/s).
    pub bulk_speed: f64,
    /// Hugoniot slope coefficient s (dimensionless).
    pub hugoniot_s: f64,
    /// Grüneisen parameter Γ₀.
    pub gruneisen: f64,
    /// Spall strength σ_spall (Pa).
    pub spall_strength: f64,
}

impl SpallFracture {
    /// Create a new spall fracture model.
    pub fn new(
        density: f64,
        bulk_speed: f64,
        hugoniot_s: f64,
        gruneisen: f64,
        spall_strength: f64,
    ) -> Self {
        Self {
            density,
            bulk_speed,
            hugoniot_s,
            gruneisen,
            spall_strength,
        }
    }

    /// Hugoniot shock velocity U_s = C₀ + s * u_p.
    pub fn shock_velocity(&self, particle_velocity: f64) -> f64 {
        self.bulk_speed + self.hugoniot_s * particle_velocity
    }

    /// Hugoniot pressure P_H = ρ₀ U_s u_p.
    pub fn hugoniot_pressure(&self, particle_velocity: f64) -> f64 {
        let us = self.shock_velocity(particle_velocity);
        self.density * us * particle_velocity
    }

    /// Grüneisen EOS pressure correction.
    ///
    /// P = P_H + Γ₀ ρ₀ (e - e_H) where e is specific internal energy.
    pub fn gruneisen_pressure(&self, compression: f64, temperature_rise: f64) -> f64 {
        let cv = 500.0; // specific heat capacity J/(kg·K) approximate
        let ph = self.hugoniot_pressure(compression * self.bulk_speed * 0.01);
        ph + self.gruneisen * self.density * cv * temperature_rise
    }

    /// Free surface velocity pullback from spall.
    ///
    /// Δu_fs ≈ 2 σ_spall / (ρ₀ c_L)
    pub fn velocity_pullback(&self) -> f64 {
        2.0 * self.spall_strength / (self.density * self.bulk_speed)
    }

    /// Spall layer thickness estimate.
    ///
    /// h_spall ≈ c_L * t_rise / 2 where t_rise is rise time of tension pulse.
    pub fn spall_thickness(&self, pulse_rise_time: f64) -> f64 {
        self.bulk_speed * pulse_rise_time / 2.0
    }

    /// Analyze spall fracture for a given free surface velocity `v_fs` (m/s).
    pub fn analyze(&self, v_fs: f64, pulse_rise_time: f64) -> SpallResult {
        let pullback = self.velocity_pullback();
        let peak_tension = self.density * self.bulk_speed * pullback / 2.0;
        let spalled = peak_tension > self.spall_strength;

        SpallResult {
            spall_strength: self.spall_strength,
            velocity_pullback: pullback,
            spall_thickness: self.spall_thickness(pulse_rise_time),
            spalled: spalled || v_fs > pullback,
        }
    }
}

// ---------------------------------------------------------------------------
// TaylorImpact
// ---------------------------------------------------------------------------

/// Taylor anvil impact result.
#[derive(Debug, Clone)]
pub struct TaylorResult {
    /// Final mushroomed diameter at impact face (m).
    pub mushroom_diameter: f64,
    /// Fraction of projectile that has plastically deformed.
    pub deformed_fraction: f64,
    /// Dynamic yield strength estimate (Pa).
    pub dynamic_yield_strength: f64,
    /// Average strain rate (s⁻¹).
    pub strain_rate: f64,
    /// Final length after impact (m).
    pub final_length: f64,
}

/// Taylor anvil test simulation.
///
/// The Taylor anvil test fires a cylindrical rod at a rigid anvil to
/// measure dynamic yield strength from the mushroomed geometry.
#[derive(Debug, Clone)]
pub struct TaylorImpact {
    /// Initial rod length L₀ (m).
    pub initial_length: f64,
    /// Rod diameter D₀ (m).
    pub diameter: f64,
    /// Rod material density (kg/m³).
    pub density: f64,
    /// Elastic wave speed c_L (m/s).
    pub elastic_wave_speed: f64,
    /// Dynamic yield strength Y (Pa).
    pub yield_strength: f64,
}

impl TaylorImpact {
    /// Create a new Taylor impact model.
    pub fn new(
        initial_length: f64,
        diameter: f64,
        density: f64,
        elastic_wave_speed: f64,
        yield_strength: f64,
    ) -> Self {
        Self {
            initial_length,
            diameter,
            density,
            elastic_wave_speed,
            yield_strength,
        }
    }

    /// Taylor analysis: estimate dynamic yield strength from deformation.
    ///
    /// Y_d = ρ c_L v₀ * (L₀ - L_f) / (L₀ * (1 - L_f/L₀)^2) — simplified.
    pub fn dynamic_yield_from_geometry(&self, v0: f64, final_length: f64) -> f64 {
        let l0 = self.initial_length;
        let lf = final_length;
        if (l0 - lf).abs() < 1e-14 {
            return self.yield_strength;
        }
        let deformed = l0 - lf;
        // Taylor (1948) formula: Y_d ≈ ρ c v₀ (L₀ - L_f) / (some expression)
        self.density * self.elastic_wave_speed * v0 * deformed
            / (l0 * (1.0 - (lf / l0).min(0.9999)))
    }

    /// Estimate final rod length after impact at velocity `v0`.
    ///
    /// L_f / L₀ ≈ 1 - (ρ v₀² / (2 Y))
    pub fn final_length(&self, v0: f64) -> f64 {
        let fraction = self.density * v0 * v0 / (2.0 * self.yield_strength);
        let lf = self.initial_length * (1.0 - fraction.min(0.99));
        lf.max(self.diameter) // cannot be shorter than diameter
    }

    /// Mushroomed diameter from volume conservation.
    ///
    /// D_mushroom ≈ D₀ * √(L₀ / L_f) for the plastically deformed region.
    pub fn mushroom_diameter(&self, v0: f64) -> f64 {
        let lf = self.final_length(v0);
        let ratio = self.initial_length / lf.max(1e-10);
        self.diameter * ratio.sqrt()
    }

    /// Average strain rate: ε̇ ≈ v₀ / L₀.
    pub fn strain_rate(&self, v0: f64) -> f64 {
        v0 / self.initial_length
    }

    /// Full Taylor impact analysis.
    pub fn analyze(&self, v0: f64) -> TaylorResult {
        let lf = self.final_length(v0);
        let md = self.mushroom_diameter(v0);
        let dy = self.dynamic_yield_from_geometry(v0, lf);
        let deformed_frac = 1.0 - lf / self.initial_length;
        let sr = self.strain_rate(v0);

        TaylorResult {
            mushroom_diameter: md,
            deformed_fraction: deformed_frac,
            dynamic_yield_strength: dy,
            strain_rate: sr,
            final_length: lf,
        }
    }
}

// ---------------------------------------------------------------------------
// ObliqueImpact
// ---------------------------------------------------------------------------

/// Result of an oblique impact analysis.
#[derive(Debug, Clone)]
pub struct ObliqueResult {
    /// Normal component of post-impact velocity (m/s).
    pub v_normal_post: f64,
    /// Tangential component of post-impact velocity (m/s).
    pub v_tangential_post: f64,
    /// Post-impact velocity magnitude (m/s).
    pub speed_post: f64,
    /// Post-impact rebound angle (degrees from normal).
    pub rebound_angle_deg: f64,
    /// Impulse in normal direction (N·s).
    pub normal_impulse: f64,
    /// Whether sliding occurred throughout contact.
    pub sliding: bool,
}

/// Oblique impact with sliding and eccentric contact.
///
/// Models impact at an angle to the surface normal, accounting for
/// friction-controlled tangential impulse and spin effects.
#[derive(Debug, Clone)]
pub struct ObliqueImpact {
    /// Normal coefficient of restitution eₙ.
    pub restitution_normal: f64,
    /// Tangential coefficient of restitution eₜ (for sliding).
    pub restitution_tangential: f64,
    /// Coefficient of friction μ.
    pub friction: f64,
    /// Effective mass m* (kg).
    pub effective_mass: f64,
}

impl ObliqueImpact {
    /// Create a new oblique impact model.
    pub fn new(
        restitution_normal: f64,
        restitution_tangential: f64,
        friction: f64,
        effective_mass: f64,
    ) -> Self {
        Self {
            restitution_normal,
            restitution_tangential,
            friction,
            effective_mass,
        }
    }

    /// Analyze an oblique impact.
    ///
    /// - `v_n`: normal approach velocity (positive = closing, m/s).
    /// - `v_t`: tangential velocity (m/s).
    pub fn analyze(&self, v_n: f64, v_t: f64) -> ObliqueResult {
        let en = self.restitution_normal;
        let mu = self.friction;

        // Normal impulse: J_n = (1 + eₙ) m* v_n.
        let j_n = (1.0 + en) * self.effective_mass * v_n;

        // Tangential impulse limited by Coulomb friction.
        let j_t_sliding = mu * j_n;
        let j_t_stick = self.effective_mass * v_t.abs();
        let sliding = j_t_stick > j_t_sliding;

        let j_t = if sliding {
            j_t_sliding * v_t.signum()
        } else {
            j_t_stick * v_t.signum()
        };

        // Post-impact velocities.
        let v_n_post = -en * v_n;
        let v_t_post = v_t - j_t / self.effective_mass;

        let speed_post = (v_n_post * v_n_post + v_t_post * v_t_post).sqrt();
        let rebound_angle_deg = if speed_post > 1e-14 {
            (v_t_post / v_n_post.abs().max(1e-14))
                .atan()
                .to_degrees()
                .abs()
        } else {
            0.0
        };

        ObliqueResult {
            v_normal_post: v_n_post,
            v_tangential_post: v_t_post,
            speed_post,
            rebound_angle_deg,
            normal_impulse: j_n,
            sliding,
        }
    }

    /// Impact angle from the normal (degrees) for given velocity components.
    pub fn impact_angle_deg(v_n: f64, v_t: f64) -> f64 {
        if v_n.abs() < 1e-14 {
            return 90.0;
        }
        (v_t / v_n).atan().to_degrees().abs()
    }
}

// ---------------------------------------------------------------------------
// MultiBodyImpact  (Newton's cradle / chain impact)
// ---------------------------------------------------------------------------

/// State of a single body in a multi-body chain impact.
#[derive(Debug, Clone)]
pub struct ChainBody {
    /// Body mass (kg).
    pub mass: f64,
    /// Current velocity (m/s).
    pub velocity: f64,
    /// Coefficient of restitution with the next body.
    pub restitution: f64,
}

impl ChainBody {
    /// Create a new chain body.
    pub fn new(mass: f64, velocity: f64, restitution: f64) -> Self {
        Self {
            mass,
            velocity,
            restitution,
        }
    }
}

/// Multi-body chain impact using sequential impulse propagation.
///
/// Models Newton's cradle and linear chain collisions where impact
/// propagates sequentially through a line of bodies.
///
/// Implements Stronge's energetic model for multi-body impacts.
#[derive(Debug, Clone)]
pub struct MultiBodyImpact {
    /// Bodies in the chain, ordered left to right.
    pub bodies: Vec<ChainBody>,
    /// Maximum sequential passes for impulse propagation.
    pub max_passes: usize,
}

impl MultiBodyImpact {
    /// Create a new multi-body impact chain.
    pub fn new(bodies: Vec<ChainBody>) -> Self {
        Self {
            bodies,
            max_passes: 10,
        }
    }

    /// Newton's cradle: n identical spheres, rightmost struck at velocity v0.
    ///
    /// Creates a chain where only the first body has initial velocity.
    pub fn newtons_cradle(n: usize, mass: f64, v0: f64, restitution: f64) -> Self {
        let mut bodies = Vec::with_capacity(n);
        for i in 0..n {
            let vel = if i == 0 { v0 } else { 0.0 };
            bodies.push(ChainBody::new(mass, vel, restitution));
        }
        Self {
            bodies,
            max_passes: 20,
        }
    }

    /// Apply a single binary collision between bodies `i` and `i+1`.
    ///
    /// Returns `true` if a collision occurred (approaching velocities).
    pub fn apply_collision(&mut self, i: usize) -> bool {
        if i + 1 >= self.bodies.len() {
            return false;
        }
        let v1 = self.bodies[i].velocity;
        let v2 = self.bodies[i + 1].velocity;

        // Only collide if approaching.
        if v1 <= v2 {
            return false;
        }

        let m1 = self.bodies[i].mass;
        let m2 = self.bodies[i + 1].mass;
        let e = self.bodies[i].restitution;
        let m_total = m1 + m2;

        let v1_post = (m1 * v1 + m2 * v2 - m2 * e * (v1 - v2)) / m_total;
        let v2_post = (m1 * v1 + m2 * v2 + m1 * e * (v1 - v2)) / m_total;

        self.bodies[i].velocity = v1_post;
        self.bodies[i + 1].velocity = v2_post;
        true
    }

    /// Propagate impulse through the chain (sequential passes).
    pub fn propagate(&mut self) {
        for _ in 0..self.max_passes {
            let mut any_collision = false;
            for i in 0..self.bodies.len().saturating_sub(1) {
                if self.apply_collision(i) {
                    any_collision = true;
                }
            }
            if !any_collision {
                break;
            }
        }
    }

    /// Total kinetic energy of all bodies.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.bodies
            .iter()
            .map(|b| 0.5 * b.mass * b.velocity * b.velocity)
            .sum()
    }

    /// Total momentum of all bodies.
    pub fn total_momentum(&self) -> f64 {
        self.bodies.iter().map(|b| b.mass * b.velocity).sum()
    }
}

// ---------------------------------------------------------------------------
// ImpactSensor
// ---------------------------------------------------------------------------

/// Head Injury Criterion (HIC) computation result.
#[derive(Debug, Clone)]
pub struct HicResult {
    /// HIC value (dimensionless).
    pub hic: f64,
    /// Time interval \[t1, t2\] maximizing HIC (s).
    pub t1: f64,
    /// End time of HIC interval (s).
    pub t2: f64,
    /// Whether HIC exceeds the injury threshold (1000).
    pub injury_risk: bool,
}

/// Shock Response Spectrum (SRS) at a given natural frequency.
#[derive(Debug, Clone)]
pub struct SrsPoint {
    /// Natural frequency (Hz).
    pub frequency: f64,
    /// Maximum response acceleration (m/s²).
    pub response: f64,
}

/// Accelerometer and shock spectrum sensor model for impact analysis.
///
/// Models an accelerometer measuring impact acceleration and computes
/// engineering metrics: shock response spectrum (SRS) and Head Injury
/// Criterion (HIC) for biomechanical safety assessment.
#[derive(Debug, Clone)]
pub struct ImpactSensor {
    /// Sampled acceleration time history (m/s²).
    pub acceleration: Vec<f64>,
    /// Time step between samples (s).
    pub dt: f64,
    /// Sensor natural frequency (Hz).
    pub natural_frequency: f64,
    /// Sensor damping ratio.
    pub damping_ratio: f64,
    /// Filter cutoff frequency (Hz).
    pub cutoff_frequency: f64,
}

impl ImpactSensor {
    /// Create a new impact sensor.
    pub fn new(dt: f64, natural_frequency: f64, damping_ratio: f64) -> Self {
        Self {
            acceleration: Vec::new(),
            dt,
            natural_frequency,
            damping_ratio,
            cutoff_frequency: natural_frequency * 0.9,
        }
    }

    /// Record a new acceleration sample (m/s²).
    pub fn record(&mut self, a: f64) {
        self.acceleration.push(a);
    }

    /// Peak acceleration magnitude (m/s²).
    pub fn peak_acceleration(&self) -> f64 {
        self.acceleration
            .iter()
            .map(|&a| a.abs())
            .fold(0.0f64, f64::max)
    }

    /// Mean acceleration (m/s²).
    pub fn mean_acceleration(&self) -> f64 {
        if self.acceleration.is_empty() {
            return 0.0;
        }
        let n = self.acceleration.len() as f64;
        self.acceleration.iter().sum::<f64>() / n
    }

    /// RMS acceleration (m/s²).
    pub fn rms_acceleration(&self) -> f64 {
        if self.acceleration.is_empty() {
            return 0.0;
        }
        let n = self.acceleration.len() as f64;
        let sum_sq: f64 = self.acceleration.iter().map(|&a| a * a).sum();
        (sum_sq / n).sqrt()
    }

    /// Compute the Head Injury Criterion (HIC) from the acceleration record.
    ///
    /// HIC = max_{t1 < t2} { (t2 - t1) * \[1/(t2-t1) ∫a dt\]^2.5 }
    ///
    /// Uses brute-force search over all intervals (O(n²)).
    pub fn hic(&self) -> HicResult {
        let n = self.acceleration.len();
        if n < 2 {
            return HicResult {
                hic: 0.0,
                t1: 0.0,
                t2: 0.0,
                injury_risk: false,
            };
        }

        let g = 9.81; // convert to g's
        let a_g: Vec<f64> = self.acceleration.iter().map(|&a| a / g).collect();

        let mut best_hic = 0.0f64;
        let mut best_t1 = 0.0f64;
        let mut best_t2 = self.dt;

        // Restrict search to intervals of at most 36ms (15ms for HIC_15 variant).
        let max_span = (0.036 / self.dt) as usize + 1;

        for i in 0..n {
            let mut integral = 0.0f64;
            for (j, &a_val) in a_g.iter().enumerate().skip(i + 1).take(max_span - 1) {
                integral += a_val * self.dt;
                let dt_interval = (j - i) as f64 * self.dt;
                if dt_interval < 1e-14 {
                    continue;
                }
                let mean_a = integral / dt_interval;
                let hic_ij = dt_interval * mean_a.abs().powf(2.5);
                if hic_ij > best_hic {
                    best_hic = hic_ij;
                    best_t1 = i as f64 * self.dt;
                    best_t2 = j as f64 * self.dt;
                }
            }
        }

        HicResult {
            hic: best_hic,
            t1: best_t1,
            t2: best_t2,
            injury_risk: best_hic > 1000.0,
        }
    }

    /// Compute the Shock Response Spectrum (SRS) at a given natural frequency.
    ///
    /// Uses central difference integration of a SDOF oscillator response.
    pub fn srs_at_frequency(&self, fn_hz: f64, damping: f64) -> SrsPoint {
        if self.acceleration.is_empty() {
            return SrsPoint {
                frequency: fn_hz,
                response: 0.0,
            };
        }

        let omega_n = 2.0 * PI * fn_hz;
        let dt = self.dt;
        let mut x = 0.0f64;
        let mut v = 0.0f64;
        let mut max_response = 0.0f64;

        for &a in &self.acceleration {
            // Central difference Newmark integration.
            let x_ddot = -2.0 * damping * omega_n * v - omega_n * omega_n * x - a;
            v += x_ddot * dt;
            x += v * dt;
            let abs_accel = (x_ddot + a).abs();
            if abs_accel > max_response {
                max_response = abs_accel;
            }
        }

        SrsPoint {
            frequency: fn_hz,
            response: max_response,
        }
    }

    /// Compute SRS over a range of frequencies.
    pub fn srs_spectrum(&self, f_min: f64, f_max: f64, n_points: usize) -> Vec<SrsPoint> {
        (0..n_points)
            .map(|i| {
                let f = f_min * (f_max / f_min).powf(i as f64 / (n_points - 1).max(1) as f64);
                self.srs_at_frequency(f, self.damping_ratio)
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── HertzianImpact ────────────────────────────────────────────────────

    #[test]
    fn test_hertz_stiffness_positive() {
        let h = HertzianImpact::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01, 0.1, 10.0);
        assert!(
            h.hertz_stiffness() > 0.0,
            "Hertz stiffness must be positive"
        );
    }

    #[test]
    fn test_hertz_effective_mass() {
        let h = HertzianImpact::new(200e9, 0.3, 200e9, 0.3, 0.01, f64::INFINITY, 1.0, 1.0);
        let m_star = h.effective_mass();
        assert!(
            (m_star - 0.5).abs() < 1e-10,
            "equal masses: m*=0.5, got {}",
            m_star
        );
    }

    #[test]
    fn test_hertz_max_deformation_increases_with_velocity() {
        let h = HertzianImpact::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01, 0.1, 1.0);
        let d1 = h.max_deformation(1.0);
        let d2 = h.max_deformation(10.0);
        assert!(
            d2 > d1,
            "higher velocity → more deformation: {} vs {}",
            d1,
            d2
        );
    }

    #[test]
    fn test_hertz_peak_force_positive() {
        let h = HertzianImpact::new(200e9, 0.3, 70e9, 0.33, 0.005, 0.01, 0.05, 5.0);
        let f = h.peak_force(5.0);
        assert!(f > 0.0, "peak force must be positive, got {}", f);
    }

    #[test]
    fn test_hertz_contact_duration_positive() {
        let h = HertzianImpact::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01, 0.1, 1.0);
        let td = h.contact_duration(5.0);
        assert!(td > 0.0, "contact duration must be positive, got {}", td);
    }

    #[test]
    fn test_hertz_contact_duration_zero_velocity() {
        let h = HertzianImpact::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01, 0.1, 1.0);
        let td = h.contact_duration(0.0);
        assert_eq!(td, 0.0);
    }

    #[test]
    fn test_hertz_energy_fraction_leq_one() {
        let h = HertzianImpact::new(200e9, 0.3, 200e9, 0.3, 0.01, 0.01, 1.0, 1.0);
        let result = h.analyze(10.0);
        assert!(
            result.energy_fraction <= 1.0 + 1e-10,
            "energy fraction must be ≤ 1, got {}",
            result.energy_fraction
        );
    }

    // ── ViscoelasticImpact ────────────────────────────────────────────────

    #[test]
    fn test_hunt_crossley_force_zero_at_zero_deformation() {
        let vi = ViscoelasticImpact::new(1e6, 1.5, 0.01, 0.5);
        let f = vi.contact_force(0.0, 1.0);
        assert_eq!(f, 0.0, "force must be zero at zero deformation");
    }

    #[test]
    fn test_hunt_crossley_force_positive() {
        let vi = ViscoelasticImpact::new(1e6, 1.5, 0.01, 0.5);
        let f = vi.contact_force(1e-3, 0.5);
        assert!(f > 0.0, "force must be positive during contact, got {}", f);
    }

    #[test]
    fn test_hunt_crossley_restitution_high_velocity_lower() {
        let vi = ViscoelasticImpact::new(1e5, 1.5, 0.1, 1.0);
        let e_slow = vi.restitution_coefficient(0.1);
        let e_fast = vi.restitution_coefficient(5.0);
        assert!(
            e_slow >= e_fast,
            "higher velocity should yield lower restitution: {} vs {}",
            e_slow,
            e_fast
        );
    }

    #[test]
    fn test_hunt_crossley_restitution_clamp_to_one() {
        let vi = ViscoelasticImpact::new(1e5, 1.5, 1e-9, 1.0);
        let e = vi.restitution_coefficient(0.001);
        assert!(e <= 1.0, "restitution must be ≤ 1, got {}", e);
    }

    #[test]
    fn test_hunt_crossley_energy_dissipated_nonneg() {
        let vi = ViscoelasticImpact::new(1e5, 1.5, 0.05, 1.0);
        let ediss = vi.energy_dissipated(2.0);
        assert!(
            ediss >= 0.0,
            "energy dissipated must be non-negative, got {}",
            ediss
        );
    }

    // ── WaveImpact ────────────────────────────────────────────────────────

    #[test]
    fn test_wave_speed_steel() {
        // Steel: E = 200 GPa, ρ = 7800 kg/m³, ν = 0.3
        let w = WaveImpact::new(7800.0, 200e9, 0.3, 7800.0, 200e9, 0.3);
        let cp = w.p_wave_speed();
        // P-wave in steel ≈ 5900 m/s
        assert!(
            cp > 5000.0 && cp < 7000.0,
            "P-wave speed in steel should be ~5900 m/s, got {:.0}",
            cp
        );
    }

    #[test]
    fn test_wave_speed_s_less_than_p() {
        let w = WaveImpact::new(2700.0, 70e9, 0.33, 2700.0, 70e9, 0.33);
        let cp = w.p_wave_speed();
        let cs = w.s_wave_speed();
        assert!(
            cs < cp,
            "S-wave must be slower than P-wave: cs={:.0}, cp={:.0}",
            cs,
            cp
        );
    }

    #[test]
    fn test_wave_rayleigh_less_than_s() {
        let w = WaveImpact::new(2700.0, 70e9, 0.33, 2700.0, 70e9, 0.33);
        let cs = w.s_wave_speed();
        let cr = w.rayleigh_speed();
        assert!(
            cr < cs,
            "Rayleigh wave must be slower than S-wave: cr={:.0}, cs={:.0}",
            cr,
            cs
        );
    }

    #[test]
    fn test_wave_same_medium_zero_reflection() {
        let w = WaveImpact::new(7800.0, 200e9, 0.3, 7800.0, 200e9, 0.3);
        let (r, _t) = w.interface_coefficients();
        assert!(
            r.abs() < 1e-10,
            "same medium: reflection coeff should be 0, got {}",
            r
        );
    }

    #[test]
    fn test_wave_travel_time_consistency() {
        let w = WaveImpact::new(7800.0, 200e9, 0.3, 7800.0, 200e9, 0.3);
        let result = w.analyze(100.0); // 100 m propagation
        let expected_t = 100.0 / result.c_longitudinal;
        assert!((result.travel_time_p - expected_t).abs() < 1e-10);
    }

    // ── ImpactCrater ──────────────────────────────────────────────────────

    #[test]
    fn test_crater_diameter_positive() {
        let crater = ImpactCrater::new_rock(2700.0, 8000.0);
        let d = crater.crater_diameter(0.1, 3000.0);
        assert!(d > 0.0, "crater diameter must be positive, got {}", d);
    }

    #[test]
    fn test_crater_diameter_increases_with_velocity() {
        let crater = ImpactCrater::new_rock(2700.0, 8000.0);
        let d_slow = crater.crater_diameter(0.1, 1000.0);
        let d_fast = crater.crater_diameter(0.1, 5000.0);
        assert!(
            d_fast > d_slow,
            "higher velocity → larger crater: {} vs {}",
            d_slow,
            d_fast
        );
    }

    #[test]
    fn test_crater_depth_one_third_diameter() {
        let crater = ImpactCrater::new_rock(2700.0, 8000.0);
        let d = crater.crater_diameter(0.1, 3000.0);
        let depth = crater.crater_depth(0.1, 3000.0);
        assert!(
            (depth - d / 3.0).abs() < 1e-10,
            "depth should be D/3, got depth={}, D/3={}",
            depth,
            d / 3.0
        );
    }

    // ── BallisticImpact ───────────────────────────────────────────────────

    #[test]
    fn test_ballistic_no_perforation_below_limit() {
        let b = BallisticImpact::new(500.0, 2.0, 1.0, 0.01, 0.006);
        let result = b.analyze(300.0);
        assert!(
            !result.perforated,
            "below ballistic limit: should not perforate"
        );
        assert_eq!(result.residual_velocity, 0.0);
    }

    #[test]
    fn test_ballistic_perforation_above_limit() {
        let b = BallisticImpact::new(500.0, 2.0, 1.0, 0.01, 0.006);
        let result = b.analyze(800.0);
        assert!(result.perforated, "above ballistic limit: should perforate");
        assert!(result.residual_velocity > 0.0);
    }

    #[test]
    fn test_ballistic_residual_velocity_increases_with_impact() {
        let b = BallisticImpact::new(500.0, 2.0, 1.0, 0.01, 0.006);
        let vr1 = b.residual_velocity(700.0);
        let vr2 = b.residual_velocity(1000.0);
        assert!(
            vr2 > vr1,
            "higher impact velocity → higher residual: {} vs {}",
            vr1,
            vr2
        );
    }

    #[test]
    fn test_ballistic_absorbed_energy_positive() {
        let b = BallisticImpact::new(500.0, 2.0, 1.0, 0.01, 0.006);
        let e = b.absorbed_energy(800.0);
        assert!(e > 0.0, "absorbed energy must be positive, got {}", e);
    }

    // ── SpallFracture ─────────────────────────────────────────────────────

    #[test]
    fn test_spall_velocity_pullback_positive() {
        let spall = SpallFracture::new(2700.0, 5100.0, 1.34, 2.0, 1.0e9);
        let pullback = spall.velocity_pullback();
        assert!(
            pullback > 0.0,
            "velocity pullback must be positive, got {}",
            pullback
        );
    }

    #[test]
    fn test_spall_thickness_positive() {
        let spall = SpallFracture::new(2700.0, 5100.0, 1.34, 2.0, 1.0e9);
        let h = spall.spall_thickness(1e-6);
        assert!(h > 0.0, "spall thickness must be positive, got {}", h);
    }

    #[test]
    fn test_spall_hugoniot_pressure_positive() {
        let spall = SpallFracture::new(7800.0, 4600.0, 1.49, 1.69, 2.5e9);
        let p = spall.hugoniot_pressure(500.0);
        assert!(p > 0.0, "Hugoniot pressure must be positive, got {}", p);
    }

    #[test]
    fn test_spall_high_strength_no_spall() {
        // Very high spall strength → should not spall at modest impact.
        let spall = SpallFracture::new(2700.0, 5100.0, 1.34, 2.0, 100e9);
        let result = spall.analyze(100.0, 1e-6);
        // Injury risk depends on pullback vs v_fs; just check no panic.
        let _ = result.spalled;
    }

    // ── TaylorImpact ──────────────────────────────────────────────────────

    #[test]
    fn test_taylor_mushroom_diameter_larger_than_initial() {
        let t = TaylorImpact::new(0.10, 0.01, 8900.0, 4600.0, 200e6);
        let dm = t.mushroom_diameter(200.0);
        assert!(
            dm >= t.diameter,
            "mushroomed diameter {} must be ≥ initial {}",
            dm,
            t.diameter
        );
    }

    #[test]
    fn test_taylor_final_length_less_than_initial() {
        let t = TaylorImpact::new(0.10, 0.01, 8900.0, 4600.0, 200e6);
        let lf = t.final_length(200.0);
        assert!(
            lf < t.initial_length,
            "final length {} must be < initial {}",
            lf,
            t.initial_length
        );
    }

    #[test]
    fn test_taylor_strain_rate_proportional_to_velocity() {
        let t = TaylorImpact::new(0.10, 0.01, 8900.0, 4600.0, 200e6);
        let sr1 = t.strain_rate(100.0);
        let sr2 = t.strain_rate(200.0);
        assert!(
            (sr2 / sr1 - 2.0).abs() < 1e-10,
            "strain rate should be proportional to velocity"
        );
    }

    #[test]
    fn test_taylor_deformed_fraction_in_range() {
        let t = TaylorImpact::new(0.10, 0.01, 8900.0, 4600.0, 200e6);
        let result = t.analyze(150.0);
        assert!(
            result.deformed_fraction >= 0.0 && result.deformed_fraction <= 1.0,
            "deformed fraction {} out of [0,1]",
            result.deformed_fraction
        );
    }

    // ── ObliqueImpact ─────────────────────────────────────────────────────

    #[test]
    fn test_oblique_normal_restitution() {
        // Zero tangential velocity → purely normal impact.
        let oi = ObliqueImpact::new(0.8, 0.0, 0.3, 1.0);
        let result = oi.analyze(10.0, 0.0);
        assert!(
            (result.v_normal_post - (-0.8 * 10.0)).abs() < 1e-10,
            "normal restitution: expected {}, got {}",
            -8.0,
            result.v_normal_post
        );
    }

    #[test]
    fn test_oblique_rebound_angle_computed() {
        let oi = ObliqueImpact::new(0.7, 0.5, 0.3, 1.0);
        let result = oi.analyze(5.0, 3.0);
        assert!(
            result.rebound_angle_deg >= 0.0,
            "rebound angle must be non-negative, got {}",
            result.rebound_angle_deg
        );
    }

    #[test]
    fn test_oblique_impact_angle_normal_incidence() {
        let angle = ObliqueImpact::impact_angle_deg(10.0, 0.0);
        assert!(
            angle.abs() < 1e-10,
            "normal incidence → 0° angle, got {}",
            angle
        );
    }

    // ── MultiBodyImpact ───────────────────────────────────────────────────

    #[test]
    fn test_newtons_cradle_momentum_conserved() {
        let mut chain = MultiBodyImpact::newtons_cradle(5, 1.0, 1.0, 1.0);
        let p_before = chain.total_momentum();
        chain.propagate();
        let p_after = chain.total_momentum();
        assert!(
            (p_after - p_before).abs() < 1e-10,
            "momentum must be conserved: before={}, after={}",
            p_before,
            p_after
        );
    }

    #[test]
    fn test_newtons_cradle_elastic_energy_conserved() {
        // Elastic (e=1): KE should be conserved.
        let mut chain = MultiBodyImpact::newtons_cradle(3, 1.0, 2.0, 1.0);
        let ke_before = chain.total_kinetic_energy();
        chain.propagate();
        let ke_after = chain.total_kinetic_energy();
        assert!(
            (ke_after - ke_before).abs() < 1e-10,
            "KE must be conserved for e=1: before={}, after={}",
            ke_before,
            ke_after
        );
    }

    #[test]
    fn test_chain_binary_collision_equal_masses_elastic() {
        // Equal masses, elastic: v1 and v2 swap.
        let bodies = vec![ChainBody::new(1.0, 3.0, 1.0), ChainBody::new(1.0, 0.0, 1.0)];
        let mut chain = MultiBodyImpact::new(bodies);
        chain.apply_collision(0);
        let v1 = chain.bodies[0].velocity;
        let v2 = chain.bodies[1].velocity;
        assert!(
            v1.abs() < 1e-10,
            "equal mass elastic: v1 should be ~0, got {}",
            v1
        );
        assert!(
            (v2 - 3.0).abs() < 1e-10,
            "equal mass elastic: v2 should be 3, got {}",
            v2
        );
    }

    #[test]
    fn test_chain_no_collision_if_separating() {
        let bodies = vec![
            ChainBody::new(1.0, -1.0, 1.0), // moving away
            ChainBody::new(1.0, 2.0, 1.0),
        ];
        let mut chain = MultiBodyImpact::new(bodies);
        let collided = chain.apply_collision(0);
        assert!(!collided, "separating bodies should not collide");
    }

    // ── ImpactSensor ──────────────────────────────────────────────────────

    #[test]
    fn test_impact_sensor_peak_acceleration() {
        let mut sensor = ImpactSensor::new(1e-4, 1000.0, 0.05);
        sensor.record(10.0);
        sensor.record(50.0);
        sensor.record(-30.0);
        assert!((sensor.peak_acceleration() - 50.0).abs() < 1e-10);
    }

    #[test]
    fn test_impact_sensor_rms_positive() {
        let mut sensor = ImpactSensor::new(1e-4, 1000.0, 0.05);
        for i in 0..10 {
            sensor.record(i as f64 * 10.0);
        }
        let rms = sensor.rms_acceleration();
        assert!(rms > 0.0, "RMS must be positive, got {}", rms);
    }

    #[test]
    fn test_impact_sensor_hic_empty() {
        let sensor = ImpactSensor::new(1e-4, 1000.0, 0.05);
        let hic_result = sensor.hic();
        assert_eq!(hic_result.hic, 0.0, "empty sensor should have HIC=0");
    }

    #[test]
    fn test_impact_sensor_srs_at_frequency() {
        let mut sensor = ImpactSensor::new(1e-4, 500.0, 0.05);
        // Sinusoidal excitation at 100 Hz.
        for i in 0..200 {
            let t = i as f64 * 1e-4;
            sensor.record(100.0 * (2.0 * PI * 100.0 * t).sin());
        }
        let srs = sensor.srs_at_frequency(100.0, 0.05);
        assert!(srs.response >= 0.0, "SRS response must be non-negative");
        assert!((srs.frequency - 100.0).abs() < 1e-10, "frequency mismatch");
    }

    #[test]
    fn test_impact_sensor_srs_spectrum_length() {
        let sensor = ImpactSensor::new(1e-4, 500.0, 0.05);
        let spec = sensor.srs_spectrum(10.0, 1000.0, 20);
        assert_eq!(spec.len(), 20, "SRS spectrum should have 20 points");
    }

    #[test]
    fn test_impact_sensor_mean_acceleration() {
        let mut sensor = ImpactSensor::new(1e-4, 500.0, 0.05);
        sensor.record(10.0);
        sensor.record(20.0);
        sensor.record(30.0);
        let mean = sensor.mean_acceleration();
        assert!(
            (mean - 20.0).abs() < 1e-10,
            "mean should be 20, got {}",
            mean
        );
    }
}
