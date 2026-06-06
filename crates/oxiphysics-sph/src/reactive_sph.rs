// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reactive SPH for combustion and detonation simulations.
//!
//! Implements Arrhenius reaction kinetics, species transport, heat release,
//! laminar flame speed, Chapman-Jouguet detonation velocity, and a simple
//! reactive SPH time-step integrator.

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ---------------------------------------------------------------------------
// SPH kernel (cubic spline, 3-D)
// ---------------------------------------------------------------------------

fn cubic_kernel_grad(r_ij: [f64; 3], h: f64) -> [f64; 3] {
    let r = len3(r_ij);
    if r < 1e-12 || h < 1e-12 {
        return [0.0; 3];
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    let dw_dr = if q < 1.0 {
        alpha * (-3.0 * q + 2.25 * q * q) / h
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * (-0.75 * t * t) / h
    } else {
        0.0
    };
    [
        r_ij[0] / r * dw_dr,
        r_ij[1] / r * dw_dr,
        r_ij[2] / r * dw_dr,
    ]
}

// ---------------------------------------------------------------------------
// ReactiveParticle
// ---------------------------------------------------------------------------

/// SPH particle carrying reactive-flow state.
#[derive(Clone, Debug)]
pub struct ReactiveParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Mass fractions of each species, sum = 1.
    pub species_fractions: Vec<f64>,
    /// Specific enthalpy h (J/kg).
    pub enthalpy: f64,
}

impl ReactiveParticle {
    /// Create a new [`ReactiveParticle`] with uniform species distribution.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        density: f64,
        temperature: f64,
        n_species: usize,
    ) -> Self {
        let frac = if n_species > 0 {
            1.0 / n_species as f64
        } else {
            0.0
        };
        Self {
            position,
            velocity,
            mass,
            density,
            pressure: 0.0,
            temperature,
            species_fractions: vec![frac; n_species],
            enthalpy: 0.0,
        }
    }

    /// Set specific species fraction; clamps to \[0, 1\].
    pub fn set_species(&mut self, idx: usize, frac: f64) {
        if idx < self.species_fractions.len() {
            self.species_fractions[idx] = frac.clamp(0.0, 1.0);
        }
    }

    /// First species fraction (fuel), or 0 if no species.
    pub fn fuel_fraction(&self) -> f64 {
        *self.species_fractions.first().unwrap_or(&0.0)
    }
}

// ---------------------------------------------------------------------------
// ReactionParams
// ---------------------------------------------------------------------------

/// Parameters for a single-step Arrhenius reaction.
#[derive(Clone, Debug)]
pub struct ReactionParams {
    /// Activation energy Eₐ (J/mol).
    pub activation_energy: f64,
    /// Pre-exponential factor A (1/s or appropriate units).
    pub pre_exponential: f64,
    /// Heat of combustion Q (J/kg of fuel).
    pub heat_of_combustion: f64,
    /// Number of species tracked.
    pub n_species: usize,
}

impl ReactionParams {
    /// Create a new [`ReactionParams`].
    pub fn new(
        activation_energy: f64,
        pre_exponential: f64,
        heat_of_combustion: f64,
        n_species: usize,
    ) -> Self {
        Self {
            activation_energy,
            pre_exponential,
            heat_of_combustion,
            n_species,
        }
    }
}

// ---------------------------------------------------------------------------
// Reaction kinetics
// ---------------------------------------------------------------------------

/// Universal gas constant R (J/mol/K).
const R_GAS: f64 = 8.314;

/// Arrhenius reaction rate ω = A · exp(−Eₐ / (R T)).
///
/// Returns 0 if `temp` ≤ 0 to avoid NaN.
pub fn arrhenius_rate(temp: f64, params: &ReactionParams) -> f64 {
    if temp <= 0.0 {
        return 0.0;
    }
    params.pre_exponential * (-(params.activation_energy / (R_GAS * temp))).exp()
}

/// Species source terms dY_k/dt for a single-step fuel-oxidiser reaction.
///
/// - Species 0 is fuel: dY_fuel/dt = −ω
/// - Species 1 is products: dY_prod/dt = +ω
/// - All others remain unchanged (dY_k/dt = 0).
pub fn species_source_term(particle: &ReactiveParticle, params: &ReactionParams) -> Vec<f64> {
    let rate = arrhenius_rate(particle.temperature, params);
    let fuel = particle.fuel_fraction();
    let omega = rate * fuel; // consumption rate (1/s)

    let n = params.n_species;
    let mut source = vec![0.0f64; n];
    if n > 0 {
        source[0] = -omega;
    } // fuel consumed
    if n > 1 {
        source[1] = omega;
    } // products formed
    source
}

/// Volumetric heat release rate Q̇ (W/kg) for a particle.
///
/// Q̇ = Q · ω · Y_fuel
pub fn heat_release_rate(particle: &ReactiveParticle, params: &ReactionParams) -> f64 {
    let rate = arrhenius_rate(particle.temperature, params);
    let fuel = particle.fuel_fraction();
    params.heat_of_combustion * rate * fuel
}

// ---------------------------------------------------------------------------
// Flame and detonation
// ---------------------------------------------------------------------------

/// Laminar flame speed estimate (m/s) using a thin-flame diffusion model.
///
/// S_L ≈ √(2 D ω_b)  where ω_b is the bulk Arrhenius rate evaluated at T_burnt.
///
/// - `temp_unburnt` – unburnt mixture temperature (K).
/// - `temp_burnt`   – adiabatic flame temperature (K).
/// - `diffusivity`  – thermal diffusivity α (m²/s).
/// - `rate`         – pre-computed Arrhenius rate at T_burnt (1/s).
pub fn flame_speed_laminar(
    _temp_unburnt: f64,
    _temp_burnt: f64,
    diffusivity: f64,
    rate: f64,
) -> f64 {
    (2.0 * diffusivity * rate).max(0.0).sqrt()
}

/// Chapman-Jouguet detonation velocity D_CJ (m/s).
///
/// D_CJ = c₀ · √(1 + (γ + 1) Q / (2 c₀²)) + √(Q (γ² − 1) / 2)
///
/// Simplified to the common form:
/// D_CJ ≈ c₀ · (1 + (γ − 1) √(Q / c₀²))
///
/// A robust approximate formula:
/// D_CJ = √(2 (γ² − 1) Q) + c₀
///
/// - `q`           – specific energy release (J/kg).
/// - `gamma`       – ratio of specific heats of products.
/// - `speed_sound` – speed of sound in unburnt mixture c₀ (m/s).
pub fn detonation_velocity_cj(q: f64, gamma: f64, speed_sound: f64) -> f64 {
    let term = 2.0 * (gamma * gamma - 1.0) * q;
    speed_sound + term.max(0.0).sqrt()
}

// ---------------------------------------------------------------------------
// Time stepping
// ---------------------------------------------------------------------------

/// Advance all reactive particles by one time step `dt`.
///
/// For each particle:
/// 1. Update species fractions from reaction source terms.
/// 2. Update temperature from heat release (assuming constant c_p = 1000 J/kg/K).
/// 3. Update pressure via ideal gas: p = ρ R_spec T (R_spec = 287 J/kg/K).
/// 4. Advance position by velocity.
pub fn reactive_sph_step(particles: &mut [ReactiveParticle], params: &ReactionParams, dt: f64) {
    const C_P: f64 = 1000.0; // J/(kg K) — constant-pressure specific heat
    const R_SPEC: f64 = 287.0; // J/(kg K) — specific gas constant (air-like)

    for p in particles.iter_mut() {
        // Reaction source
        let src = species_source_term(p, params);
        for (frac, &ds) in p.species_fractions.iter_mut().zip(src.iter()) {
            *frac = (*frac + ds * dt).clamp(0.0, 1.0);
        }
        // Re-normalise species fractions
        let total: f64 = p.species_fractions.iter().sum();
        if total > 1e-300 {
            for frac in p.species_fractions.iter_mut() {
                *frac /= total;
            }
        }

        // Heat release → temperature rise
        let q_dot = heat_release_rate(p, params);
        p.temperature += q_dot * dt / C_P;
        p.temperature = p.temperature.max(0.0);

        // Pressure (ideal gas)
        if p.density > 1e-300 {
            p.pressure = p.density * R_SPEC * p.temperature;
        }

        // Position advection
        p.position = add3(p.position, scale3(p.velocity, dt));
    }
}

// ---------------------------------------------------------------------------
// Ignition detection
// ---------------------------------------------------------------------------

/// Return indices of particles whose temperature exceeds `temp_threshold`.
pub fn detect_ignition(particles: &[ReactiveParticle], temp_threshold: f64) -> Vec<usize> {
    particles
        .iter()
        .enumerate()
        .filter_map(|(i, p)| {
            if p.temperature >= temp_threshold {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Supplementary functions (SPH heat/species diffusion helpers)
// ---------------------------------------------------------------------------

/// SPH diffusion contribution to species transport for particle `i`.
///
/// Returns the scalar source Σ_j m_j/ρ_j * 2 (Y_i - Y_j) / |r_ij|² * ∇W · r_ij
/// multiplied by the diffusivity D.
pub fn sph_species_diffusion(
    particles: &[ReactiveParticle],
    i: usize,
    h: f64,
    species_idx: usize,
    diffusivity: f64,
) -> f64 {
    let pi = &particles[i];
    let yi = *pi.species_fractions.get(species_idx).unwrap_or(&0.0);
    let mut result = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let r_ij = sub3(pi.position, pj.position);
        let r2 = dot3(r_ij, r_ij);
        if r2 < 1e-300 || pj.density < 1e-300 {
            continue;
        }
        let yj = *pj.species_fractions.get(species_idx).unwrap_or(&0.0);
        let grad_w = cubic_kernel_grad(r_ij, h);
        let dot_grad = dot3(grad_w, r_ij);
        result += pj.mass / pj.density * 2.0 * (yi - yj) / r2 * dot_grad;
    }
    diffusivity * result
}

/// SPH thermal diffusion contribution dT/dt for particle `i`.
///
/// Uses the Brookshaw operator with thermal conductivity assumed uniform λ = 1.
pub fn sph_thermal_diffusion_reactive(
    particles: &[ReactiveParticle],
    i: usize,
    h: f64,
    thermal_diffusivity: f64,
) -> f64 {
    let pi = &particles[i];
    let mut result = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        if j == i {
            continue;
        }
        let r_ij = sub3(pi.position, pj.position);
        let r2 = dot3(r_ij, r_ij);
        if r2 < 1e-300 || pj.density < 1e-300 {
            continue;
        }
        let dt = pi.temperature - pj.temperature;
        let grad_w = cubic_kernel_grad(r_ij, h);
        let dot_grad = dot3(grad_w, r_ij);
        result += pj.mass / pj.density * 2.0 * dt / r2 * dot_grad;
    }
    thermal_diffusivity * result
}

/// Approximate adiabatic flame temperature.
///
/// T_ad = T_unburnt + Q / c_p
pub fn adiabatic_flame_temperature(t_unburnt: f64, q: f64, c_p: f64) -> f64 {
    if c_p < 1e-300 {
        return t_unburnt;
    }
    t_unburnt + q / c_p
}

/// Normalise species fractions so they sum to 1.
pub fn normalise_species(fracs: &mut [f64]) {
    let total: f64 = fracs.iter().sum();
    if total > 1e-300 {
        for f in fracs.iter_mut() {
            *f /= total;
        }
    }
}

/// Compute the mixture molecular weight given species molecular weights and fractions.
pub fn mixture_molecular_weight(fracs: &[f64], mol_weights: &[f64]) -> f64 {
    // 1 / M_mix = Σ Y_k / M_k
    let mut inv = 0.0;
    for (&y, &m) in fracs.iter().zip(mol_weights.iter()) {
        if m > 1e-300 {
            inv += y / m;
        }
    }
    if inv < 1e-300 {
        return 0.0;
    }
    1.0 / inv
}

/// Compute the mole fractions from mass fractions and molecular weights.
pub fn mass_to_mole_fractions(fracs: &[f64], mol_weights: &[f64]) -> Vec<f64> {
    let m_mix = mixture_molecular_weight(fracs, mol_weights);
    if m_mix < 1e-300 {
        return vec![0.0; fracs.len()];
    }
    fracs
        .iter()
        .zip(mol_weights.iter())
        .map(|(&y, &m)| if m > 1e-300 { y * m_mix / m } else { 0.0 })
        .collect()
}

/// Estimate density from ideal-gas law: ρ = p M / (R_universal T).
pub fn ideal_gas_density(pressure: f64, mol_weight: f64, temperature: f64) -> f64 {
    if temperature < 1e-300 {
        return 0.0;
    }
    pressure * mol_weight / (R_GAS * temperature * 1000.0)
    // mol_weight in g/mol → kg/mol by /1000
}

/// Compute the specific heat at constant pressure for a mixture.
///
/// c_p_mix = Σ Y_k c_p_k
pub fn mixture_cp(fracs: &[f64], cp_values: &[f64]) -> f64 {
    fracs
        .iter()
        .zip(cp_values.iter())
        .map(|(&y, &cp)| y * cp)
        .sum()
}

/// Count the number of particles with non-negligible fuel (Y_fuel > threshold).
pub fn count_burning_particles(particles: &[ReactiveParticle], threshold: f64) -> usize {
    particles
        .iter()
        .filter(|p| p.fuel_fraction() > threshold)
        .count()
}

/// Compute the average temperature of all particles.
pub fn average_temperature(particles: &[ReactiveParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let sum: f64 = particles.iter().map(|p| p.temperature).sum();
    sum / particles.len() as f64
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- helpers ---

    fn make_params() -> ReactionParams {
        ReactionParams::new(
            50_000.0, // Eₐ = 50 kJ/mol
            1.0e6,    // A = 1e6 /s
            43.0e6,   // Q = 43 MJ/kg (propane-like)
            2,        // fuel + products
        )
    }

    fn make_hot_particle() -> ReactiveParticle {
        let mut p = ReactiveParticle::new([0.0; 3], [0.0; 3], 1e-6, 1.2, 1200.0, 2);
        p.set_species(0, 0.9); // 90% fuel
        p.set_species(1, 0.1);
        p
    }

    fn make_cold_particle() -> ReactiveParticle {
        let mut p = ReactiveParticle::new([0.0; 3], [0.0; 3], 1e-6, 1.2, 300.0, 2);
        p.set_species(0, 0.9);
        p.set_species(1, 0.1);
        p
    }

    // --- ReactiveParticle ---

    #[test]
    fn test_particle_new_uniform_species() {
        let p = ReactiveParticle::new([0.0; 3], [0.0; 3], 1e-6, 1.0, 300.0, 4);
        let total: f64 = p.species_fractions.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "Uniform fracs sum to 1: {total}"
        );
    }

    #[test]
    fn test_particle_new_zero_species() {
        let p = ReactiveParticle::new([0.0; 3], [0.0; 3], 1e-6, 1.0, 300.0, 0);
        assert!(p.species_fractions.is_empty());
        assert_eq!(p.fuel_fraction(), 0.0);
    }

    #[test]
    fn test_particle_set_species_clamp() {
        let mut p = ReactiveParticle::new([0.0; 3], [0.0; 3], 1e-6, 1.0, 300.0, 2);
        p.set_species(0, 1.5);
        assert_eq!(p.species_fractions[0], 1.0, "Should clamp to 1");
        p.set_species(1, -0.5);
        assert_eq!(p.species_fractions[1], 0.0, "Should clamp to 0");
    }

    #[test]
    fn test_particle_fuel_fraction() {
        let mut p = ReactiveParticle::new([0.0; 3], [0.0; 3], 1e-6, 1.0, 300.0, 3);
        p.set_species(0, 0.7);
        assert!((p.fuel_fraction() - 0.7).abs() < 1e-12);
    }

    // --- ReactionParams ---

    #[test]
    fn test_params_new() {
        let p = make_params();
        assert_eq!(p.n_species, 2);
        assert!(p.activation_energy > 0.0);
        assert!(p.pre_exponential > 0.0);
        assert!(p.heat_of_combustion > 0.0);
    }

    // --- arrhenius_rate ---

    #[test]
    fn test_arrhenius_rate_zero_temp() {
        let params = make_params();
        assert_eq!(arrhenius_rate(0.0, &params), 0.0);
    }

    #[test]
    fn test_arrhenius_rate_negative_temp() {
        let params = make_params();
        assert_eq!(arrhenius_rate(-100.0, &params), 0.0);
    }

    #[test]
    fn test_arrhenius_rate_positive_at_high_temp() {
        let params = make_params();
        let rate = arrhenius_rate(2000.0, &params);
        assert!(rate > 0.0, "Rate must be positive: {rate}");
    }

    #[test]
    fn test_arrhenius_rate_increases_with_temp() {
        let params = make_params();
        let r1 = arrhenius_rate(500.0, &params);
        let r2 = arrhenius_rate(1000.0, &params);
        assert!(r2 > r1, "Rate must increase with temperature");
    }

    // --- species_source_term ---

    #[test]
    fn test_species_source_fuel_decreases() {
        let p = make_hot_particle();
        let params = make_params();
        let src = species_source_term(&p, &params);
        assert!(src[0] <= 0.0, "Fuel source must be ≤ 0");
    }

    #[test]
    fn test_species_source_products_increase() {
        let p = make_hot_particle();
        let params = make_params();
        let src = species_source_term(&p, &params);
        assert!(src[1] >= 0.0, "Products source must be ≥ 0");
    }

    #[test]
    fn test_species_source_cold_nearly_zero() {
        let p = make_cold_particle();
        let params = make_params();
        let src = species_source_term(&p, &params);
        // At 300 K, rate is negligible for high Eₐ
        assert!(src[0].abs() < 1.0, "Cold particle: source ≈ 0 : {}", src[0]);
    }

    #[test]
    fn test_species_source_zero_fuel() {
        let mut p = make_hot_particle();
        p.set_species(0, 0.0);
        p.set_species(1, 1.0);
        let params = make_params();
        let src = species_source_term(&p, &params);
        assert_eq!(src[0], 0.0, "No fuel → no consumption");
    }

    // --- heat_release_rate ---

    #[test]
    fn test_heat_release_positive_for_hot_particle() {
        let p = make_hot_particle();
        let params = make_params();
        let q = heat_release_rate(&p, &params);
        assert!(q >= 0.0, "Heat release must be ≥ 0: {q}");
    }

    #[test]
    fn test_heat_release_zero_for_no_fuel() {
        let mut p = make_hot_particle();
        p.set_species(0, 0.0);
        p.set_species(1, 1.0);
        let params = make_params();
        let q = heat_release_rate(&p, &params);
        assert_eq!(q, 0.0, "No fuel → no heat release");
    }

    // --- flame_speed_laminar ---

    #[test]
    fn test_flame_speed_positive() {
        let sl = flame_speed_laminar(300.0, 2200.0, 2e-5, 1000.0);
        assert!(sl >= 0.0, "Flame speed must be ≥ 0: {sl}");
    }

    #[test]
    fn test_flame_speed_zero_diffusivity() {
        let sl = flame_speed_laminar(300.0, 2200.0, 0.0, 1000.0);
        assert_eq!(sl, 0.0);
    }

    #[test]
    fn test_flame_speed_zero_rate() {
        let sl = flame_speed_laminar(300.0, 2200.0, 2e-5, 0.0);
        assert_eq!(sl, 0.0);
    }

    // --- detonation_velocity_cj ---

    #[test]
    fn test_detonation_velocity_positive() {
        let d_cj = detonation_velocity_cj(2.0e6, 1.4, 350.0);
        assert!(d_cj > 350.0, "D_CJ must exceed speed of sound: {d_cj}");
    }

    #[test]
    fn test_detonation_velocity_increases_with_q() {
        let d1 = detonation_velocity_cj(1.0e6, 1.4, 350.0);
        let d2 = detonation_velocity_cj(4.0e6, 1.4, 350.0);
        assert!(d2 > d1, "Higher Q → higher D_CJ");
    }

    #[test]
    fn test_detonation_velocity_zero_q() {
        // With Q=0, D_CJ = c₀
        let c0 = 350.0;
        let d = detonation_velocity_cj(0.0, 1.4, c0);
        assert!((d - c0).abs() < 1e-10, "D_CJ = c₀ when Q=0: {d}");
    }

    // --- reactive_sph_step ---

    #[test]
    fn test_step_temperature_increases_for_hot_fuel() {
        let mut particles = vec![make_hot_particle()];
        let params = make_params();
        let t0 = particles[0].temperature;
        reactive_sph_step(&mut particles, &params, 1e-4);
        assert!(particles[0].temperature >= t0, "Temp should not decrease");
    }

    #[test]
    fn test_step_fuel_decreases() {
        let mut particles = vec![make_hot_particle()];
        let params = make_params();
        let f0 = particles[0].fuel_fraction();
        reactive_sph_step(&mut particles, &params, 1e-3);
        assert!(
            particles[0].fuel_fraction() <= f0,
            "Fuel should not increase"
        );
    }

    #[test]
    fn test_step_species_sum_remains_one() {
        let mut particles = vec![make_hot_particle()];
        let params = make_params();
        reactive_sph_step(&mut particles, &params, 1e-3);
        let total: f64 = particles[0].species_fractions.iter().sum();
        assert!((total - 1.0).abs() < 1e-10, "Species sum = {total}");
    }

    #[test]
    fn test_step_position_advances() {
        let mut p = ReactiveParticle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1e-6, 1.0, 300.0, 2);
        p.set_species(0, 0.5);
        p.set_species(1, 0.5);
        let params = make_params();
        let dt = 0.1;
        reactive_sph_step(
            &mut vec![p.clone()].into_iter().collect::<Vec<_>>(),
            &params,
            dt,
        );
        // Re-create to test position
        let mut particles = vec![p];
        reactive_sph_step(&mut particles, &params, dt);
        assert!(
            (particles[0].position[0] - 0.1).abs() < 1e-10,
            "x = {}",
            particles[0].position[0]
        );
    }

    #[test]
    fn test_step_empty_particles() {
        let mut particles: Vec<ReactiveParticle> = Vec::new();
        let params = make_params();
        reactive_sph_step(&mut particles, &params, 0.01); // should not panic
    }

    // --- detect_ignition ---

    #[test]
    fn test_detect_ignition_hot_particle() {
        let particles = vec![make_cold_particle(), make_hot_particle()];
        let ignited = detect_ignition(&particles, 1000.0);
        assert_eq!(ignited, vec![1], "Only hot particle ignites");
    }

    #[test]
    fn test_detect_ignition_none() {
        let particles = vec![make_cold_particle()];
        let ignited = detect_ignition(&particles, 1000.0);
        assert!(ignited.is_empty(), "Cold particle below threshold");
    }

    #[test]
    fn test_detect_ignition_all() {
        let particles = vec![make_hot_particle(), make_hot_particle()];
        let ignited = detect_ignition(&particles, 500.0);
        assert_eq!(ignited.len(), 2, "Both particles above 500 K");
    }

    #[test]
    fn test_detect_ignition_empty() {
        let ignited = detect_ignition(&[], 1000.0);
        assert!(ignited.is_empty());
    }

    // --- auxiliary functions ---

    #[test]
    fn test_adiabatic_flame_temperature() {
        let t_ad = adiabatic_flame_temperature(300.0, 43.0e6, 1000.0);
        assert!(t_ad > 300.0, "T_ad must exceed T_unburnt");
        assert!((t_ad - 43300.0).abs() < 1.0, "T_ad = {t_ad}");
    }

    #[test]
    fn test_normalise_species() {
        let mut fracs = vec![2.0, 2.0, 4.0];
        normalise_species(&mut fracs);
        let total: f64 = fracs.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "Normalised sum = {total}");
    }

    #[test]
    fn test_normalise_species_all_zero() {
        let mut fracs = vec![0.0, 0.0, 0.0];
        normalise_species(&mut fracs); // should not panic
        let total: f64 = fracs.iter().sum();
        assert_eq!(total, 0.0);
    }

    #[test]
    fn test_mixture_molecular_weight() {
        // Equal fracs, MW = [2, 4] → 1/M = 0.5/2 + 0.5/4 = 0.375 → M = 2.667
        let fracs = vec![0.5, 0.5];
        let mw = vec![2.0, 4.0];
        let m = mixture_molecular_weight(&fracs, &mw);
        assert!((m - 8.0 / 3.0).abs() < 1e-10, "M_mix = {m}");
    }

    #[test]
    fn test_mass_to_mole_fractions_pure() {
        // Pure species 0 (MW=2): mole fraction = 1
        let fracs = vec![1.0, 0.0];
        let mw = vec![2.0, 4.0];
        let x = mass_to_mole_fractions(&fracs, &mw);
        assert!(
            (x[0] - 1.0).abs() < 1e-10,
            "Pure species mole frac = {}",
            x[0]
        );
    }

    #[test]
    fn test_count_burning_particles() {
        let mut p1 = make_cold_particle();
        p1.set_species(0, 0.01); // barely any fuel
        p1.set_species(1, 0.99);
        let p2 = make_hot_particle();
        let count = count_burning_particles(&[p1, p2], 0.05);
        assert_eq!(count, 1, "Only hot particle has fuel > 0.05");
    }

    #[test]
    fn test_average_temperature_single() {
        let p = make_hot_particle();
        let avg = average_temperature(std::slice::from_ref(&p));
        assert!((avg - p.temperature).abs() < 1e-10);
    }

    #[test]
    fn test_average_temperature_empty() {
        assert_eq!(average_temperature(&[]), 0.0);
    }

    #[test]
    fn test_mixture_cp_single_species() {
        let fracs = vec![1.0];
        let cps = vec![1005.0];
        assert!((mixture_cp(&fracs, &cps) - 1005.0).abs() < 1e-10);
    }

    #[test]
    fn test_sph_species_diffusion_single_particle() {
        let p = make_hot_particle();
        // Only one particle: no neighbours → diffusion = 0
        let result = sph_species_diffusion(&[p], 0, 0.01, 0, 1.0);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_sph_thermal_diffusion_reactive_single() {
        let p = make_hot_particle();
        let result = sph_thermal_diffusion_reactive(&[p], 0, 0.01, 1.0);
        assert_eq!(result, 0.0);
    }

    #[test]
    fn test_ideal_gas_density_positive() {
        let rho = ideal_gas_density(101325.0, 28.97, 300.0);
        assert!(rho > 0.0, "Density must be positive: {rho}");
    }

    #[test]
    fn test_flame_speed_increases_with_rate() {
        let s1 = flame_speed_laminar(300.0, 2200.0, 2e-5, 100.0);
        let s2 = flame_speed_laminar(300.0, 2200.0, 2e-5, 10000.0);
        assert!(s2 > s1, "Higher rate → higher flame speed");
    }

    #[test]
    fn test_step_pressure_nonzero_after_step() {
        let mut particles = vec![make_hot_particle()];
        particles[0].density = 1.2;
        let params = make_params();
        reactive_sph_step(&mut particles, &params, 1e-3);
        assert!(
            particles[0].pressure > 0.0,
            "Pressure should be > 0 after step"
        );
    }
}
