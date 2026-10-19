// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermal SPH for heat transfer and phase change.
//!
//! Implements:
//! - [`ThermalParticle`] struct with temperature, thermal conductivity, specific heat
//! - [`compute_heat_flux`] – SPH discretisation of the heat equation (Cleary–Monaghan /
//!   Brookshaw conduction operator)
//! - [`phase_change`] – synchronise the liquid fraction and the specific enthalpy with a
//!   particle's temperature (it does not change the temperature itself)
//! - [`thermal_conduction_step`] – explicit Euler time integration of the sensible heat
//! - [`thermal_conduction_step_with_phase`] / [`integrate_enthalpy`] – explicit conduction
//!   step with the enthalpy method: the conducted heat is added to the specific enthalpy
//!   `e(T) = c_p T + f L`, which is then inverted for the temperature and the liquid
//!   fraction `f`, so latent heat holds the temperature inside the mushy zone
//! - [`nucleation`] – nucleation-site probability model
//! - [`solidification_front`] – Stefan problem front velocity

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

// ---------------------------------------------------------------------------
// Cubic-spline SPH kernel
// ---------------------------------------------------------------------------

/// Cubic-spline SPH kernel W(r, h).
///
/// Returns the kernel value at distance `r` for smoothing length `h`.
pub fn cubic_kernel(r: f64, h: f64) -> f64 {
    if h < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        alpha * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Derivative of the cubic-spline kernel dW/dr.
///
/// Returns dW/dr (scalar) for distance `r` and smoothing length `h`.
pub fn cubic_kernel_deriv(r: f64, h: f64) -> f64 {
    if h < 1e-300 || r < 1e-300 {
        return 0.0;
    }
    let q = r / h;
    let alpha = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        alpha / h * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        -alpha / h * 0.75 * t * t
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// ThermalParticle
// ---------------------------------------------------------------------------

/// A single thermal SPH particle carrying thermodynamic state.
pub struct ThermalParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Thermal conductivity λ (W/(m·K)).
    pub thermal_conductivity: f64,
    /// Specific heat at constant pressure c_p (J/(kg·K)).
    pub specific_heat: f64,
    /// Rate of temperature change dT/dt (K/s) caused by conduction, i.e. the
    /// sensible-heat rate `c_p dT/dt = (1/ρ) ∇·(λ∇T)` — updated by
    /// [`compute_heat_flux`] each step.
    pub d_temperature: f64,
    /// Phase indicator (liquid fraction): 0.0 = fully solid, 1.0 = fully liquid.
    pub phase: f64,
    /// Specific enthalpy per unit mass (J/kg).  [`ThermalParticle::new`] stores
    /// the sensible part `c_p T`; [`phase_change`] and
    /// [`thermal_conduction_step_with_phase`] store `c_p T + f L`.
    pub internal_energy: f64,
}

impl ThermalParticle {
    /// Create a new thermal particle with default dynamic fields set to zero.
    ///
    /// # Arguments
    /// * `position`            – initial position (m)
    /// * `mass`                – particle mass (kg)
    /// * `density`             – initial density (kg/m³)
    /// * `h`                   – smoothing length (m)
    /// * `temperature`         – initial temperature (K)
    /// * `thermal_conductivity` – λ (W/(m·K))
    /// * `specific_heat`       – c_p (J/(kg·K))
    pub fn new(
        position: [f64; 3],
        mass: f64,
        density: f64,
        h: f64,
        temperature: f64,
        thermal_conductivity: f64,
        specific_heat: f64,
    ) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density,
            h,
            temperature,
            thermal_conductivity,
            specific_heat,
            d_temperature: 0.0,
            phase: if temperature >= 273.15 { 1.0 } else { 0.0 },
            internal_energy: specific_heat * temperature,
        }
    }
}

// ---------------------------------------------------------------------------
// Heat flux computation (Brookshaw operator)
// ---------------------------------------------------------------------------

/// Compute the SPH heat-equation source term dT/dt for all particles.
///
/// Uses the Cleary–Monaghan (1999) conduction operator, which is the Brookshaw
/// (1985) discretisation of `∇·(λ∇T)` with the harmonic-mean pair conductivity:
///
/// ```text
/// ρᵢ c_pᵢ dTᵢ/dt = Σⱼ (mⱼ/ρⱼ) · 4 λᵢ λⱼ/(λᵢ+λⱼ) · (Tᵢ−Tⱼ) (rᵢⱼ·∇ᵢWᵢⱼ) / (|rᵢⱼ|²+ε)
/// ```
///
/// with `rᵢⱼ = xᵢ − xⱼ`, `rᵢⱼ·∇ᵢWᵢⱼ = |rᵢⱼ| ∂W/∂r` (negative for a decreasing
/// kernel, so heat flows from hot to cold), pair smoothing length
/// `hᵢⱼ = (hᵢ+hⱼ)/2` and `ε = 10⁻⁴ hᵢⱼ²`.  For a uniform conductivity the
/// harmonic mean `4λ²/(2λ)` equals `2λ`, so the sum is exactly the Brookshaw
/// form `λ ∇²T`; the pair terms are antisymmetric, so `Σ mᵢ c_pᵢ dTᵢ/dt = 0`.
///
/// Results are written into `particle.d_temperature` for each particle.
///
/// # Arguments
/// * `particles` – mutable slice of all thermal particles
pub fn compute_heat_flux(particles: &mut [ThermalParticle]) {
    let n = particles.len();
    // Collect positional/scalar state to avoid borrow conflicts.
    let pos: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
    let temp: Vec<f64> = particles.iter().map(|p| p.temperature).collect();
    let mass: Vec<f64> = particles.iter().map(|p| p.mass).collect();
    let rho: Vec<f64> = particles.iter().map(|p| p.density).collect();
    let lam: Vec<f64> = particles.iter().map(|p| p.thermal_conductivity).collect();
    let cp: Vec<f64> = particles.iter().map(|p| p.specific_heat).collect();
    let h_sph: Vec<f64> = particles.iter().map(|p| p.h).collect();

    let mut d_temp = vec![0.0f64; n];

    for i in 0..n {
        let mut sum = 0.0f64;
        for j in 0..n {
            if i == j || rho[j] <= 0.0 {
                continue;
            }
            let rij = sub3(pos[i], pos[j]);
            let r = len3(rij);
            let h_avg = 0.5 * (h_sph[i] + h_sph[j]);
            let dw_dr = cubic_kernel_deriv(r, h_avg);
            if dw_dr == 0.0 {
                continue;
            }
            let lam_ij = if lam[i] + lam[j] > 0.0 {
                4.0 * lam[i] * lam[j] / (lam[i] + lam[j])
            } else {
                0.0
            };
            let eps = 1e-4 * h_avg * h_avg;
            // r_ij · ∇_i W_ij = r ∂W/∂r
            let contrib =
                mass[j] / rho[j] * lam_ij * (temp[i] - temp[j]) * r * dw_dr / (r * r + eps);
            sum += contrib;
        }
        // c_p ρ dT/dt = Σ ...
        d_temp[i] = if cp[i] * rho[i] > 0.0 {
            sum / (cp[i] * rho[i])
        } else {
            0.0
        };
    }

    for (i, p) in particles.iter_mut().enumerate() {
        p.d_temperature = d_temp[i];
    }
}

// ---------------------------------------------------------------------------
// Phase change – enthalpy method
// ---------------------------------------------------------------------------

/// Thermal parameters for a phase-change material.
pub struct PhaseChangeMaterial {
    /// Melting temperature T_melt (K).
    pub melting_temperature: f64,
    /// Latent heat of fusion L (J/kg).
    pub latent_heat: f64,
    /// Solidus temperature (K) — start of mushy zone.
    pub solidus_temperature: f64,
    /// Liquidus temperature (K) — end of mushy zone.
    pub liquidus_temperature: f64,
}

impl PhaseChangeMaterial {
    /// Water/ice properties (approximate).
    pub fn water_ice() -> Self {
        Self {
            melting_temperature: 273.15,
            latent_heat: 334_000.0,
            solidus_temperature: 272.0,
            liquidus_temperature: 274.0,
        }
    }

    /// Aluminium properties (approximate).
    pub fn aluminium() -> Self {
        Self {
            melting_temperature: 933.5,
            latent_heat: 397_000.0,
            solidus_temperature: 930.0,
            liquidus_temperature: 937.0,
        }
    }

    /// Width `T_liquidus − T_solidus` of the mushy zone (K); zero for a material
    /// that melts at a single temperature (or a malformed, inverted zone).
    pub fn mushy_zone_width(&self) -> f64 {
        (self.liquidus_temperature - self.solidus_temperature).max(0.0)
    }

    /// Liquid fraction `f(T)` of the material: 0 below the solidus, 1 above the
    /// liquidus and linear in between.
    ///
    /// A material with a zero-width mushy zone melts isothermally: `f` is a step
    /// at the solidus temperature, and exactly at that temperature the fraction
    /// is not determined by `T` alone (`0.5` is returned here; the particle
    /// routines use the particle's own liquid fraction instead).
    pub fn liquid_fraction(&self, temperature: f64) -> f64 {
        self.liquid_fraction_with_state(temperature, 0.5)
    }

    /// Liquid fraction at `temperature`; `stored` is only consulted for the
    /// isothermal case at exactly the melting temperature, where the fraction is
    /// a free state variable of the particle.
    fn liquid_fraction_with_state(&self, temperature: f64, stored: f64) -> f64 {
        let ts = self.solidus_temperature;
        let width = self.mushy_zone_width();
        if width > 0.0 {
            if temperature <= ts {
                0.0
            } else if temperature >= self.liquidus_temperature {
                1.0
            } else {
                (temperature - ts) / width
            }
        } else if temperature < ts {
            0.0
        } else if temperature > ts {
            1.0
        } else {
            stored.clamp(0.0, 1.0)
        }
    }

    /// Specific enthalpy `e = c_p T + f L` (J/kg) of a particle with specific heat
    /// `specific_heat`, temperature `temperature` and liquid fraction `liquid_fraction`.
    pub fn specific_enthalpy(
        &self,
        specific_heat: f64,
        temperature: f64,
        liquid_fraction: f64,
    ) -> f64 {
        specific_heat * temperature + liquid_fraction * self.latent_heat
    }

    /// Invert the enthalpy relation: the temperature (K) and liquid fraction that
    /// belong to the specific enthalpy `enthalpy` (J/kg) of a material with
    /// specific heat `specific_heat` (J/(kg·K)).
    ///
    /// With `ΔT = T_liquidus − T_solidus`:
    /// * `e ≤ c_p T_s` (solid): `T = e/c_p`, `f = 0`;
    /// * `e ≥ c_p T_l + L` (liquid): `T = (e − L)/c_p`, `f = 1`;
    /// * otherwise (mushy zone): `T = (e + L T_s/ΔT)/(c_p + L/ΔT)` and
    ///   `f = (T − T_s)/ΔT`.  For `ΔT = 0` the temperature stays at `T_s` and
    ///   `f = (e − c_p T_s)/L`.
    ///
    /// The map `T ↦ e` is strictly increasing for `c_p > 0`, so this is its exact inverse.
    pub fn temperature_and_fraction(&self, specific_heat: f64, enthalpy: f64) -> (f64, f64) {
        let ts = self.solidus_temperature;
        let latent = self.latent_heat;
        let width = self.mushy_zone_width();
        let e_solid = specific_heat * ts;
        let e_liquid = specific_heat * (ts + width) + latent;
        if enthalpy <= e_solid {
            (enthalpy / specific_heat, 0.0)
        } else if enthalpy >= e_liquid {
            ((enthalpy - latent) / specific_heat, 1.0)
        } else if width > 0.0 {
            let t = (enthalpy + latent * ts / width) / (specific_heat + latent / width);
            (t, ((t - ts) / width).clamp(0.0, 1.0))
        } else {
            (ts, ((enthalpy - e_solid) / latent).clamp(0.0, 1.0))
        }
    }
}

/// Synchronise a particle's liquid fraction and specific enthalpy with its temperature.
///
/// Sets `particle.phase` to the liquid fraction `f(T)` of the mushy-zone model
/// (see [`PhaseChangeMaterial::liquid_fraction`]) and `particle.internal_energy`
/// to the specific enthalpy `c_p T + f L`.  The temperature itself is **not**
/// modified: this call is a diagnostic/initialisation step that puts a particle
/// into a state consistent with the enthalpy method.  Latent heat only enters the
/// temperature evolution through [`thermal_conduction_step_with_phase`] (or
/// [`integrate_enthalpy`]), which advance the enthalpy and invert it for `T`.
///
/// For a material that melts at a single temperature the liquid fraction of a
/// particle sitting exactly at that temperature is a free state variable and is
/// kept (clamped to `[0, 1]`).
///
/// # Arguments
/// * `particle` – mutable reference to the thermal particle
/// * `material` – phase-change material parameters
pub fn phase_change(particle: &mut ThermalParticle, material: &PhaseChangeMaterial) {
    particle.phase = material.liquid_fraction_with_state(particle.temperature, particle.phase);
    particle.internal_energy =
        material.specific_enthalpy(particle.specific_heat, particle.temperature, particle.phase);
}

/// Advance the specific enthalpy of every particle by its current conduction
/// rate and recover the temperature and liquid fraction (enthalpy method).
///
/// For each particle the enthalpy is rebuilt from its present temperature and
/// liquid fraction, `e = c_p T + f L`, advanced with the conduction rate stored in
/// `d_temperature` (`e ← e + c_p dT/dt · dt`, the rate is the *sensible-heat* rate
/// that [`compute_heat_flux`] writes) and inverted with
/// [`PhaseChangeMaterial::temperature_and_fraction`].  While a particle is inside
/// the mushy zone the supplied heat therefore melts (or, when removed, freezes)
/// material instead of changing the temperature: the plateau lasts
/// `(c_p ΔT + L)/(c_p |dT/dt|)` seconds.
///
/// `particle.internal_energy` is updated with the new enthalpy.  Particles with a
/// non-positive specific heat are left untouched, and the temperature is clamped
/// at absolute zero as in [`thermal_conduction_step`].
///
/// # Arguments
/// * `particles` – mutable slice of all thermal particles
/// * `dt`        – time step (s)
/// * `material`  – phase-change material parameters
pub fn integrate_enthalpy(
    particles: &mut [ThermalParticle],
    dt: f64,
    material: &PhaseChangeMaterial,
) {
    for p in particles.iter_mut() {
        let cp = p.specific_heat;
        if !cp.is_finite() || cp <= 0.0 {
            continue;
        }
        let fraction = material.liquid_fraction_with_state(p.temperature, p.phase);
        let enthalpy =
            material.specific_enthalpy(cp, p.temperature, fraction) + cp * p.d_temperature * dt;
        if !enthalpy.is_finite() {
            continue;
        }
        let (mut temperature, fraction) = material.temperature_and_fraction(cp, enthalpy);
        let mut enthalpy = enthalpy;
        if temperature < 0.0 {
            temperature = 0.0;
            enthalpy = material.specific_enthalpy(cp, 0.0, fraction);
        }
        p.temperature = temperature;
        p.phase = fraction;
        p.internal_energy = enthalpy;
    }
}

/// Advance all particles by one explicit conduction step **with latent heat**.
///
/// Calls [`compute_heat_flux`] and then [`integrate_enthalpy`]: the conducted heat
/// is added to the specific enthalpy `e = c_p T + f L` and the temperature and
/// liquid fraction are recovered from it, so the temperature of a melting or
/// freezing particle stays inside `[T_solidus, T_liquidus]` until the latent heat
/// has been absorbed or released.  The total enthalpy `Σ mᵢ eᵢ` is conserved.
///
/// The particles should carry a state consistent with the material (use
/// [`phase_change`] once after setting the temperatures); the enthalpy is always
/// rebuilt from the temperature, so a stale `internal_energy` cannot corrupt the step.
///
/// # Arguments
/// * `particles` – mutable slice of all thermal particles
/// * `dt`        – time step (s)
/// * `material`  – phase-change material parameters
pub fn thermal_conduction_step_with_phase(
    particles: &mut [ThermalParticle],
    dt: f64,
    material: &PhaseChangeMaterial,
) {
    compute_heat_flux(particles);
    integrate_enthalpy(particles, dt, material);
}

// ---------------------------------------------------------------------------
// Explicit thermal conduction step
// ---------------------------------------------------------------------------

/// Advance temperatures of all particles by one explicit Euler step.
///
/// Calls [`compute_heat_flux`] first, then integrates T += dT/dt · dt.  Only the
/// sensible heat is modelled; use [`thermal_conduction_step_with_phase`] when
/// melting or solidification has to absorb or release latent heat.
///
/// # Arguments
/// * `particles` – mutable slice of all thermal particles
/// * `dt`        – time step (s)
pub fn thermal_conduction_step(particles: &mut [ThermalParticle], dt: f64) {
    compute_heat_flux(particles);
    for p in particles.iter_mut() {
        p.temperature += p.d_temperature * dt;
        // Clamp to physical minimum
        if p.temperature < 0.0 {
            p.temperature = 0.0;
        }
    }
}

// ---------------------------------------------------------------------------
// Nucleation
// ---------------------------------------------------------------------------

/// Classical nucleation theory: nucleation rate J (m⁻³ s⁻¹).
///
/// Uses the Volmer–Weber exponential model:
///
/// ```text
/// J = A · exp(−ΔG* / (k_B · T))
/// ```
///
/// where ΔG* = 16π γ³ / (3 ΔGᵥ²) and ΔGᵥ is the volumetric free energy
/// difference approximated as L_f · ρ_solid · ΔT / T_melt.
///
/// # Arguments
/// * `temperature`         – current temperature (K)
/// * `melting_temperature` – equilibrium melting point (K)
/// * `latent_heat`         – latent heat of fusion (J/kg)
/// * `density_solid`       – density of the solid phase (kg/m³)
/// * `surface_energy`      – solid–liquid interface energy γ (J/m²)
/// * `prefactor`           – kinetic prefactor A (m⁻³ s⁻¹)
pub fn nucleation(
    temperature: f64,
    melting_temperature: f64,
    latent_heat: f64,
    density_solid: f64,
    surface_energy: f64,
    prefactor: f64,
) -> f64 {
    if temperature >= melting_temperature {
        return 0.0;
    }
    let delta_t = melting_temperature - temperature;
    let k_b = 1.380_649e-23_f64; // J/K
    // Volumetric free energy difference approximation
    let dg_v = latent_heat * density_solid * delta_t / melting_temperature;
    if dg_v.abs() < 1e-300 {
        return 0.0;
    }
    let dg_star = 16.0 * PI * surface_energy.powi(3) / (3.0 * dg_v * dg_v);
    prefactor * (-dg_star / (k_b * temperature.max(1e-6))).exp()
}

// ---------------------------------------------------------------------------
// Stefan problem – solidification front velocity
// ---------------------------------------------------------------------------

/// Stefan problem: front velocity from energy balance at the solid–liquid interface.
///
/// The Stefan condition reads:
///
/// ```text
/// v_n = (λ_s · ∂T_s/∂n − λ_l · ∂T_l/∂n) / (ρ_s · L_f)
/// ```
///
/// Approximated here from discrete temperature gradients on each side.
///
/// # Arguments
/// * `heat_flux_solid`  – λ_s · ∂T/∂n on the solid side (W/m²)
/// * `heat_flux_liquid` – λ_l · ∂T/∂n on the liquid side (W/m²)
/// * `density_solid`    – ρ_s (kg/m³)
/// * `latent_heat`      – L_f (J/kg)
///
/// Returns the normal front velocity v_n (m/s).
pub fn solidification_front(
    heat_flux_solid: f64,
    heat_flux_liquid: f64,
    density_solid: f64,
    latent_heat: f64,
) -> f64 {
    let denom = density_solid * latent_heat;
    if denom.abs() < 1e-300 {
        return 0.0;
    }
    (heat_flux_solid - heat_flux_liquid) / denom
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod conduction_tests;

#[cfg(test)]
mod tests {
    use super::*;

    // -- ThermalParticle::new ------------------------------------------------

    #[test]
    fn test_thermal_particle_new_temperature() {
        let p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, 300.0, 0.6, 4186.0);
        assert_eq!(p.temperature, 300.0);
    }

    #[test]
    fn test_thermal_particle_new_velocity_zero() {
        let p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, 300.0, 0.6, 4186.0);
        assert_eq!(p.velocity, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_thermal_particle_phase_liquid_above_melt() {
        let p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, 300.0, 0.6, 4186.0);
        assert_eq!(p.phase, 1.0);
    }

    #[test]
    fn test_thermal_particle_phase_solid_below_melt() {
        let p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, 200.0, 0.6, 4186.0);
        assert_eq!(p.phase, 0.0);
    }

    #[test]
    fn test_thermal_particle_internal_energy() {
        let cp = 4186.0;
        let t = 300.0;
        let p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, t, 0.6, cp);
        assert!((p.internal_energy - cp * t).abs() < 1e-6);
    }

    // -- cubic_kernel --------------------------------------------------------

    #[test]
    fn test_cubic_kernel_zero_beyond_support() {
        let w = cubic_kernel(2.1, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_cubic_kernel_positive_within_support() {
        let w = cubic_kernel(0.5, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_cubic_kernel_peak_at_zero() {
        let w0 = cubic_kernel(0.0, 1.0);
        let w1 = cubic_kernel(0.5, 1.0);
        assert!(w0 > w1);
    }

    #[test]
    fn test_cubic_kernel_zero_h() {
        let w = cubic_kernel(0.5, 0.0);
        assert_eq!(w, 0.0);
    }

    // -- cubic_kernel_deriv --------------------------------------------------

    #[test]
    fn test_kernel_deriv_zero_beyond_support() {
        let dw = cubic_kernel_deriv(2.1, 1.0);
        assert_eq!(dw, 0.0);
    }

    #[test]
    fn test_kernel_deriv_negative_within_support() {
        let dw = cubic_kernel_deriv(0.5, 1.0);
        assert!(dw < 0.0, "kernel derivative should be negative: {dw}");
    }

    // -- compute_heat_flux ---------------------------------------------------

    #[test]
    fn test_heat_flux_two_equal_particles_zero() {
        // Two particles at the same temperature → zero heat flux
        let mut particles = vec![
            ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.05, 300.0, 0.6, 4186.0),
            ThermalParticle::new([0.04, 0.0, 0.0], 1e-3, 1000.0, 0.05, 300.0, 0.6, 4186.0),
        ];
        compute_heat_flux(&mut particles);
        assert!(particles[0].d_temperature.abs() < 1e-10);
        assert!(particles[1].d_temperature.abs() < 1e-10);
    }

    #[test]
    fn test_heat_flux_direction() {
        // Hot particle cools, cold particle warms
        let mut particles = vec![
            ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.05, 400.0, 0.6, 4186.0),
            ThermalParticle::new([0.04, 0.0, 0.0], 1e-3, 1000.0, 0.05, 200.0, 0.6, 4186.0),
        ];
        compute_heat_flux(&mut particles);
        assert!(particles[0].d_temperature < 0.0, "hot should cool");
        assert!(particles[1].d_temperature > 0.0, "cold should warm");
    }

    #[test]
    fn test_heat_flux_energy_conservation() {
        // Sum of m*c_p*dT/dt should be ~0 (energy conserved)
        let mut particles = vec![
            ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.05, 400.0, 0.6, 4186.0),
            ThermalParticle::new([0.04, 0.0, 0.0], 1e-3, 1000.0, 0.05, 200.0, 0.6, 4186.0),
        ];
        compute_heat_flux(&mut particles);
        let total: f64 = particles
            .iter()
            .map(|p| p.mass * p.specific_heat * p.d_temperature)
            .sum();
        // Should be near zero (anti-symmetric fluxes)
        assert!(total.abs() < 1e-8, "energy should be conserved: {total}");
    }

    #[test]
    fn test_heat_flux_single_particle_zero() {
        let mut particles = vec![ThermalParticle::new(
            [0.0, 0.0, 0.0],
            1e-3,
            1000.0,
            0.05,
            300.0,
            0.6,
            4186.0,
        )];
        compute_heat_flux(&mut particles);
        assert_eq!(particles[0].d_temperature, 0.0);
    }

    // -- thermal_conduction_step --------------------------------------------

    #[test]
    fn test_conduction_step_updates_temperature() {
        let mut particles = vec![
            ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.05, 400.0, 0.6, 4186.0),
            ThermalParticle::new([0.04, 0.0, 0.0], 1e-3, 1000.0, 0.05, 200.0, 0.6, 4186.0),
        ];
        let t0 = particles[0].temperature;
        let t1 = particles[1].temperature;
        thermal_conduction_step(&mut particles, 1.0);
        // Temperatures should change
        assert!((particles[0].temperature - t0).abs() > 1e-15);
        assert!((particles[1].temperature - t1).abs() > 1e-15);
    }

    #[test]
    fn test_conduction_step_temperature_non_negative() {
        let mut particles = vec![ThermalParticle::new(
            [0.0, 0.0, 0.0],
            1e-3,
            1000.0,
            0.05,
            0.001,
            0.6,
            4186.0,
        )];
        thermal_conduction_step(&mut particles, 1000.0);
        assert!(particles[0].temperature >= 0.0);
    }

    // -- phase_change --------------------------------------------------------

    #[test]
    fn test_phase_change_fully_solid() {
        let material = PhaseChangeMaterial::water_ice();
        let mut p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 900.0, 0.01, 260.0, 2.2, 2090.0);
        phase_change(&mut p, &material);
        assert_eq!(p.phase, 0.0);
    }

    #[test]
    fn test_phase_change_fully_liquid() {
        let material = PhaseChangeMaterial::water_ice();
        let mut p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, 280.0, 0.6, 4186.0);
        phase_change(&mut p, &material);
        assert_eq!(p.phase, 1.0);
    }

    #[test]
    fn test_phase_change_mushy_zone() {
        let material = PhaseChangeMaterial::water_ice();
        // Temperature exactly halfway between solidus and liquidus
        let t_mid = 0.5 * (material.solidus_temperature + material.liquidus_temperature);
        let mut p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 1000.0, 0.01, t_mid, 0.6, 4186.0);
        phase_change(&mut p, &material);
        assert!((p.phase - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_phase_change_aluminium_above_liquidus() {
        let material = PhaseChangeMaterial::aluminium();
        let mut p = ThermalParticle::new([0.0, 0.0, 0.0], 1e-3, 2700.0, 0.01, 950.0, 200.0, 900.0);
        phase_change(&mut p, &material);
        assert_eq!(p.phase, 1.0);
    }

    // -- nucleation ----------------------------------------------------------

    #[test]
    fn test_nucleation_zero_above_melt() {
        let j = nucleation(280.0, 273.15, 334_000.0, 900.0, 0.033, 1e30);
        assert_eq!(j, 0.0);
    }

    #[test]
    fn test_nucleation_increases_with_supercooling() {
        // More supercooling → higher nucleation rate
        let j1 = nucleation(270.0, 273.15, 334_000.0, 900.0, 0.033, 1e30);
        let j2 = nucleation(265.0, 273.15, 334_000.0, 900.0, 0.033, 1e30);
        assert!(j2 >= j1, "nucleation should increase with supercooling");
    }

    #[test]
    fn test_nucleation_finite() {
        let j = nucleation(260.0, 273.15, 334_000.0, 900.0, 0.033, 1e10);
        assert!(j.is_finite());
        assert!(j >= 0.0);
    }

    #[test]
    fn test_nucleation_at_melt_zero() {
        let j = nucleation(273.15, 273.15, 334_000.0, 900.0, 0.033, 1e30);
        assert_eq!(j, 0.0);
    }

    // -- solidification_front ------------------------------------------------

    #[test]
    fn test_solidification_front_positive_flux_diff() {
        let v = solidification_front(1000.0, 500.0, 900.0, 334_000.0);
        assert!(v > 0.0);
    }

    #[test]
    fn test_solidification_front_zero_when_equal() {
        let v = solidification_front(1000.0, 1000.0, 900.0, 334_000.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_solidification_front_zero_latent_heat() {
        let v = solidification_front(1000.0, 500.0, 900.0, 0.0);
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_solidification_front_scaling() {
        // Doubling density halves the front velocity
        let v1 = solidification_front(1000.0, 0.0, 900.0, 334_000.0);
        let v2 = solidification_front(1000.0, 0.0, 1800.0, 334_000.0);
        assert!((v1 / v2 - 2.0).abs() < 1e-10);
    }
}
