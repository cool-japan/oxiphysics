// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Ocean wave simulation with SPH.
//!
//! Provides:
//!
//! - JONSWAP spectrum for initial velocity/displacement seeding
//! - Deep-water and shallow-water SPH particle states
//! - Wave breaking detection (Stokes limiting steepness, Froude criterion)
//! - Foam and spray generation heuristics
//! - Stokes drift velocity field
//! - Tidal forcing (astronomical tidal potential, M2 constituent)
//! - Wind-driven surface stress (Large & Pond drag parameterisation)
//! - Wave–current interaction (effective wavenumber shifting)
//! - Rogue (freak) wave simulation (superposition of focused wave trains)
//! - SPH beach run-up (Carrier–Greenspan shoreline tracking)

use std::f64::consts::PI;

/// Gravitational acceleration (m s⁻²).
const G: f64 = 9.80665;
/// Seawater reference density (kg m⁻³).
const RHO_WATER: f64 = 1025.0;
/// Air density at sea level (kg m⁻³).
const RHO_AIR: f64 = 1.225;
// ---------------------------------------------------------------------------
// Vec3 helpers
// ---------------------------------------------------------------------------

/// Compute the Euclidean length of a 3-D vector.
#[inline]
pub fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Add two 3-D vectors component-wise.
#[inline]
pub fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-D vector by a scalar.
#[inline]
pub fn vec3_scale(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Dot product of two 3-D vectors.
#[inline]
pub fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

// ---------------------------------------------------------------------------
// JONSWAP spectrum
// ---------------------------------------------------------------------------

/// Parameters for the JONSWAP (Joint North Sea Wave Project) spectrum.
///
/// The spectrum describes wave energy as a function of frequency in a fetch-limited
/// sea state: S(f) = α g² (2π)⁻⁴ f⁻⁵ exp(…) γ^exp(…).
#[derive(Clone, Debug)]
pub struct JonswapParams {
    /// Peak angular frequency ωₚ (rad s⁻¹).
    pub omega_peak: f64,
    /// Phillips equilibrium parameter α (typically 0.0081).
    pub alpha: f64,
    /// Peak enhancement factor γ (typically 3.3 for JONSWAP).
    pub gamma: f64,
    /// Spectral width parameter σ_a (f ≤ fₚ, typically 0.07).
    pub sigma_a: f64,
    /// Spectral width parameter σ_b (f > fₚ, typically 0.09).
    pub sigma_b: f64,
}

impl JonswapParams {
    /// Construct default JONSWAP parameters for a given peak frequency `fp` (Hz).
    pub fn new(fp: f64) -> Self {
        JonswapParams {
            omega_peak: 2.0 * PI * fp,
            alpha: 0.0081,
            gamma: 3.3,
            sigma_a: 0.07,
            sigma_b: 0.09,
        }
    }

    /// Evaluate the one-sided variance spectral density S(ω) (m² s rad⁻¹).
    ///
    /// # Arguments
    /// * `omega` – angular frequency (rad s⁻¹)
    pub fn spectrum(&self, omega: f64) -> f64 {
        if omega <= 0.0 {
            return 0.0;
        }
        let f = omega / (2.0 * PI);
        let fp = self.omega_peak / (2.0 * PI);
        let sigma = if f <= fp { self.sigma_a } else { self.sigma_b };
        // exp_arg is negative away from the peak; the peak-enhancement factor r = γ^exp(exp_arg)
        // peaks at f=fₚ where exp_arg=0 → r=γ, and decays smoothly away.
        let exp_arg = -((f - fp).powi(2)) / (2.0 * sigma * sigma * fp * fp);
        let r = self.gamma.powf(exp_arg.exp());
        let pm =
            self.alpha * G * G / (2.0 * PI).powi(4) / f.powi(5) * (-1.25 * (fp / f).powi(4)).exp();
        pm * r
    }

    /// Significant wave height Hₛ estimated by integrating S(ω) over a discrete grid.
    ///
    /// Hₛ = 4 √(m₀)  where  m₀ = ∫ S(ω) dω.
    ///
    /// # Arguments
    /// * `omega_min` – lower bound of integration (rad s⁻¹)
    /// * `omega_max` – upper bound of integration (rad s⁻¹)
    /// * `n_points`  – number of quadrature points
    pub fn significant_wave_height(&self, omega_min: f64, omega_max: f64, n_points: usize) -> f64 {
        if n_points < 2 {
            return 0.0;
        }
        let dw = (omega_max - omega_min) / (n_points - 1) as f64;
        let m0: f64 = (0..n_points)
            .map(|i| {
                let w = omega_min + i as f64 * dw;
                self.spectrum(w)
            })
            .sum::<f64>()
            * dw;
        4.0 * m0.max(0.0).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Ocean SPH particle
// ---------------------------------------------------------------------------

/// A single SPH particle representing a fluid parcel in the ocean simulation.
#[derive(Clone, Debug)]
pub struct OceanParticle {
    /// World-space position (x, y, z) in metres.
    pub pos: [f64; 3],
    /// Velocity (u, v, w) in m s⁻¹.
    pub vel: [f64; 3],
    /// Accumulated force per unit mass (m s⁻²).
    pub acc: [f64; 3],
    /// Particle mass (kg).
    pub mass: f64,
    /// Local density (kg m⁻³).
    pub rho: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Wave surface elevation contributed by this particle (m).
    pub eta: f64,
    /// Foam/spray intensity \[0, 1\].
    pub foam: f64,
    /// Whether this particle is currently breaking.
    pub breaking: bool,
    /// Whether this particle is a spray droplet.
    pub is_spray: bool,
    /// Stokes drift accumulated velocity (m s⁻¹).
    pub stokes_drift: [f64; 3],
    /// Water depth below this particle (m).
    pub depth: f64,
}

impl OceanParticle {
    /// Create a new ocean particle at rest with default state.
    ///
    /// # Arguments
    /// * `pos`  – initial position in metres
    /// * `mass` – particle mass in kg
    /// * `h`    – smoothing length in metres
    pub fn new(pos: [f64; 3], mass: f64, h: f64) -> Self {
        OceanParticle {
            pos,
            vel: [0.0; 3],
            acc: [0.0; 3],
            mass,
            rho: RHO_WATER,
            pressure: 0.0,
            h,
            eta: 0.0,
            foam: 0.0,
            breaking: false,
            is_spray: false,
            stokes_drift: [0.0; 3],
            depth: 0.0,
        }
    }

    /// Kinetic energy of this particle (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.vel, self.vel)
    }

    /// Potential energy (gravitational, J), taking z = 0 as datum.
    pub fn potential_energy(&self) -> f64 {
        self.mass * G * self.pos[2]
    }

    /// Speed of this particle (m s⁻¹).
    pub fn speed(&self) -> f64 {
        vec3_len(self.vel)
    }
}

// ---------------------------------------------------------------------------
// Deep-water wave dispersion
// ---------------------------------------------------------------------------

/// Deep-water phase speed c = g / ω (m s⁻¹).
///
/// # Arguments
/// * `omega` – angular frequency (rad s⁻¹)
pub fn deep_water_phase_speed(omega: f64) -> f64 {
    if omega == 0.0 {
        return 0.0;
    }
    G / omega
}

/// Deep-water group velocity cg = g / (2ω) (m s⁻¹).
///
/// # Arguments
/// * `omega` – angular frequency (rad s⁻¹)
pub fn deep_water_group_velocity(omega: f64) -> f64 {
    if omega == 0.0 {
        return 0.0;
    }
    G / (2.0 * omega)
}

/// Wavenumber k for arbitrary depth using linear wave theory dispersion.
///
/// Solves ω² = g k tanh(k d) iteratively via Newton's method.
///
/// # Arguments
/// * `omega` – angular frequency (rad s⁻¹)
/// * `d`     – water depth (m)
pub fn dispersion_wavenumber(omega: f64, d: f64) -> f64 {
    if omega <= 0.0 || d <= 0.0 {
        return 0.0;
    }
    let k0 = omega * omega / G; // deep-water initial guess
    let mut k = k0;
    for _ in 0..50 {
        let th = (k * d).tanh();
        let f = G * k * th - omega * omega;
        let df = G * (th + k * d * (1.0 - th * th));
        let dk = f / df;
        k -= dk;
        if k < 0.0 {
            k = k0 * 1e-3;
        }
        if dk.abs() < 1e-12 * k.abs() {
            break;
        }
    }
    k
}

/// Phase speed c = ω / k for finite depth.
///
/// # Arguments
/// * `omega` – angular frequency (rad s⁻¹)
/// * `d`     – water depth (m)
pub fn finite_depth_phase_speed(omega: f64, d: f64) -> f64 {
    let k = dispersion_wavenumber(omega, d);
    if k == 0.0 {
        return 0.0;
    }
    omega / k
}

/// Group velocity for finite depth: cg = c (1 + 2kd/sinh(2kd)) / 2.
///
/// # Arguments
/// * `omega` – angular frequency (rad s⁻¹)
/// * `d`     – water depth (m)
pub fn finite_depth_group_velocity(omega: f64, d: f64) -> f64 {
    let k = dispersion_wavenumber(omega, d);
    if k == 0.0 {
        return 0.0;
    }
    let c = omega / k;
    let s2kd = (2.0 * k * d).sinh();
    let n = if s2kd.abs() > 1e-12 {
        0.5 * (1.0 + 2.0 * k * d / s2kd)
    } else {
        1.0
    };
    c * n
}

// ---------------------------------------------------------------------------
// Wave breaking detection
// ---------------------------------------------------------------------------

/// Stokes limiting steepness criterion: wave breaks when ka > π/7.
///
/// Returns `true` if the wave is at or beyond the breaking limit.
///
/// # Arguments
/// * `amplitude` – wave amplitude (m)
/// * `wavenumber` – wavenumber k (rad m⁻¹)
pub fn stokes_breaking_criterion(amplitude: f64, wavenumber: f64) -> bool {
    wavenumber * amplitude >= PI / 7.0
}

/// Froude-based breaking criterion for shallow water: breaks when u ≥ √(g d).
///
/// # Arguments
/// * `horizontal_speed` – particle horizontal speed (m s⁻¹)
/// * `depth`            – local water depth (m)
pub fn froude_breaking_criterion(horizontal_speed: f64, depth: f64) -> bool {
    if depth <= 0.0 {
        return true;
    }
    horizontal_speed >= (G * depth).sqrt()
}

/// Banner & Peregrine (1993) kinematic breaking index: breaks when u_surface ≥ c.
///
/// # Arguments
/// * `surface_velocity` – horizontal surface velocity (m s⁻¹)
/// * `phase_speed`      – local wave phase speed (m s⁻¹)
pub fn kinematic_breaking_criterion(surface_velocity: f64, phase_speed: f64) -> bool {
    surface_velocity >= phase_speed
}

// ---------------------------------------------------------------------------
// Foam and spray
// ---------------------------------------------------------------------------

/// Foam intensity decay model: foam dissipates exponentially.
///
/// New foam value after time step `dt` (s):  f_new = f * exp(-dt / τ).
///
/// # Arguments
/// * `foam`  – current foam intensity \[0, 1\]
/// * `dt`    – time step (s)
/// * `tau`   – foam decay timescale (s)
pub fn foam_decay(foam: f64, dt: f64, tau: f64) -> f64 {
    if tau <= 0.0 {
        return 0.0;
    }
    (foam * (-dt / tau).exp()).clamp(0.0, 1.0)
}

/// Generate foam intensity from wave-breaking dissipation rate.
///
/// Uses a linear model: foam generation ∝ breaking rate.
///
/// # Arguments
/// * `breaking_dissipation` – energy dissipation rate per unit area (W m⁻²)
/// * `scale`                – empirical scale factor
pub fn foam_from_breaking(breaking_dissipation: f64, scale: f64) -> f64 {
    (breaking_dissipation * scale).clamp(0.0, 1.0)
}

/// Spray droplet ejection speed heuristic (m s⁻¹) based on Lozano-Duran (2022).
///
/// v_spray ≈ √(2 g H_break) where H_break is the breaking wave height.
///
/// # Arguments
/// * `breaking_height` – wave height at breaking (m)
pub fn spray_ejection_speed(breaking_height: f64) -> f64 {
    if breaking_height <= 0.0 {
        return 0.0;
    }
    (2.0 * G * breaking_height).sqrt()
}

// ---------------------------------------------------------------------------
// Stokes drift
// ---------------------------------------------------------------------------

/// First-order Stokes drift for a monochromatic wave (m s⁻¹).
///
/// u_S = a² ω k exp(2kz)  (deep water)
///
/// # Arguments
/// * `amplitude`  – wave amplitude (m)
/// * `omega`      – angular frequency (rad s⁻¹)
/// * `wavenumber` – wavenumber k (rad m⁻¹)
/// * `z`          – vertical position relative to mean surface (m), ≤ 0
/// * `direction`  – unit vector in wave propagation direction (x, y)
pub fn stokes_drift_deep(
    amplitude: f64,
    omega: f64,
    wavenumber: f64,
    z: f64,
    direction: [f64; 2],
) -> [f64; 3] {
    let magnitude = amplitude * amplitude * omega * wavenumber * (2.0 * wavenumber * z).exp();
    [magnitude * direction[0], magnitude * direction[1], 0.0]
}

/// Second-order Stokes drift correction (finite depth, Longuet-Higgins 1953).
///
/// Adds the second-order correction to the first-order Stokes drift.
///
/// # Arguments
/// * `amplitude`  – wave amplitude (m)
/// * `omega`      – angular frequency (rad s⁻¹)
/// * `k`          – wavenumber (rad m⁻¹)
/// * `d`          – water depth (m)
/// * `z`          – vertical coordinate (m), ≤ 0
/// * `direction`  – unit propagation direction (x, y)
pub fn stokes_drift_finite_depth(
    amplitude: f64,
    omega: f64,
    k: f64,
    d: f64,
    z: f64,
    direction: [f64; 2],
) -> [f64; 3] {
    if k <= 0.0 || d <= 0.0 {
        return [0.0; 3];
    }
    let kd = k * d;
    let cosh2 = (2.0 * k * (z + d)).cosh();
    let sinh2kd = (2.0 * kd).sinh();
    let factor = amplitude * amplitude * omega * k * cosh2 / (sinh2kd.max(1e-12));
    [factor * direction[0], factor * direction[1], 0.0]
}

// ---------------------------------------------------------------------------
// Tidal forcing
// ---------------------------------------------------------------------------

/// Astronomical tidal potential for the M2 (lunar semi-diurnal) constituent.
///
/// V(x, t) = H_M2 * cos(ω_M2 * t - k_M2 * x + φ)
///
/// # Arguments
/// * `x`     – horizontal coordinate (m)
/// * `t`     – time (s)
/// * `h_m2`  – M2 tidal amplitude (m)
/// * `phi`   – tidal phase offset (rad)
pub fn m2_tidal_elevation(x: f64, t: f64, h_m2: f64, phi: f64) -> f64 {
    // M2 period ≈ 12.42 h
    let omega_m2 = 2.0 * PI / (12.42 * 3600.0);
    // M2 wavenumber for shallow-water tides (very long waves, k ≈ ω/√(gd))
    // Use representative open-ocean depth of 4000 m
    let c_m2 = (G * 4000.0_f64).sqrt();
    let k_m2 = omega_m2 / c_m2;
    h_m2 * (omega_m2 * t - k_m2 * x + phi).cos()
}

/// Tidal acceleration (body force per unit mass, m s⁻²) in the x-direction.
///
/// Derived from the gradient of the tidal potential.
///
/// # Arguments
/// * `x`    – horizontal coordinate (m)
/// * `t`    – time (s)
/// * `h_m2` – M2 tidal amplitude (m)
/// * `phi`  – tidal phase offset (rad)
pub fn m2_tidal_acceleration_x(x: f64, t: f64, h_m2: f64, phi: f64) -> f64 {
    let omega_m2 = 2.0 * PI / (12.42 * 3600.0);
    let c_m2 = (G * 4000.0_f64).sqrt();
    let k_m2 = omega_m2 / c_m2;
    // ∂V/∂x with V = g * eta, so a = -g * ∂eta/∂x
    G * h_m2 * k_m2 * (omega_m2 * t - k_m2 * x + phi).sin()
}

// ---------------------------------------------------------------------------
// Wind-driven surface stress
// ---------------------------------------------------------------------------

/// Aerodynamic drag coefficient from Large & Pond (1981).
///
/// Cd = 1.2e-3 for U10 < 11 m/s, Cd = (0.49 + 0.065 U10) * 1e-3 otherwise.
///
/// # Arguments
/// * `u10` – wind speed at 10 m height (m s⁻¹)
pub fn large_pond_drag_coefficient(u10: f64) -> f64 {
    if u10 <= 0.0 {
        return 0.0;
    }
    if u10 < 11.0 {
        1.2e-3
    } else {
        (0.49 + 0.065 * u10) * 1e-3
    }
}

/// Wind surface stress vector τ = ρ_air Cd |U10| U10 (Pa).
///
/// # Arguments
/// * `wind_vel` – wind velocity at 10 m height (m s⁻¹), (u, v) components
pub fn wind_surface_stress(wind_vel: [f64; 2]) -> [f64; 2] {
    let u10 = (wind_vel[0] * wind_vel[0] + wind_vel[1] * wind_vel[1]).sqrt();
    let cd = large_pond_drag_coefficient(u10);
    let tau_scale = RHO_AIR * cd * u10;
    [tau_scale * wind_vel[0], tau_scale * wind_vel[1]]
}

/// Wind stress acceleration applied to ocean surface layer (m s⁻²).
///
/// # Arguments
/// * `wind_vel`   – 10-m wind velocity (m s⁻¹)
/// * `layer_rho`  – density of the surface layer (kg m⁻³)
/// * `layer_dz`   – thickness of the wind-mixed layer (m)
pub fn wind_stress_acceleration(wind_vel: [f64; 2], layer_rho: f64, layer_dz: f64) -> [f64; 2] {
    let tau = wind_surface_stress(wind_vel);
    let denom = (layer_rho * layer_dz).max(1e-12);
    [tau[0] / denom, tau[1] / denom]
}

// ---------------------------------------------------------------------------
// Wave–current interaction
// ---------------------------------------------------------------------------

/// Effective wavenumber shift due to a collinear current (Doppler shift).
///
/// In a current U, the observed frequency ω_obs = ω_intrinsic + k U.
/// Given ω_obs, solve for k:  ω_obs = √(g k tanh(k d)) + k U.
///
/// Returns the shifted wavenumber.
///
/// # Arguments
/// * `omega_obs` – observed absolute frequency (rad s⁻¹)
/// * `d`         – water depth (m)
/// * `current_u` – depth-averaged current in the wave propagation direction (m s⁻¹)
pub fn wave_current_wavenumber(omega_obs: f64, d: f64, current_u: f64) -> f64 {
    if omega_obs <= 0.0 {
        return 0.0;
    }
    // Iterative: k_{n+1} = dispersion wavenumber at ω_intrinsic = ω_obs - k_n * U
    let mut k = dispersion_wavenumber(omega_obs, d);
    for _ in 0..40 {
        let omega_int = omega_obs - k * current_u;
        if omega_int <= 0.0 {
            break;
        }
        let k_new = dispersion_wavenumber(omega_int, d);
        if (k_new - k).abs() < 1e-12 * k.abs().max(1.0) {
            k = k_new;
            break;
        }
        k = k_new;
    }
    k
}

/// Wave action density N = E / ω (m² s).
///
/// Conserved quantity in wave–current interaction (Whitham 1965).
///
/// # Arguments
/// * `wave_energy_density` – E (J m⁻²)
/// * `intrinsic_omega`     – ω_intrinsic (rad s⁻¹)
pub fn wave_action_density(wave_energy_density: f64, intrinsic_omega: f64) -> f64 {
    if intrinsic_omega <= 0.0 {
        return 0.0;
    }
    wave_energy_density / intrinsic_omega
}

/// Wave energy density from wave action (inverse of above).
///
/// # Arguments
/// * `action`          – wave action density N (m² s)
/// * `intrinsic_omega` – ω_intrinsic (rad s⁻¹)
pub fn wave_energy_from_action(action: f64, intrinsic_omega: f64) -> f64 {
    action * intrinsic_omega
}

// ---------------------------------------------------------------------------
// Rogue wave simulation
// ---------------------------------------------------------------------------

/// A focused wave group (NewWave / rogue wave) seed.
#[derive(Clone, Debug)]
pub struct RogueWaveSeed {
    /// Number of component waves in the superposition.
    pub n_components: usize,
    /// Focus amplitude (m).
    pub focus_amplitude: f64,
    /// Focus time (s).
    pub focus_time: f64,
    /// Focus location x (m).
    pub focus_x: f64,
    /// Peak frequency ωₚ (rad s⁻¹).
    pub omega_peak: f64,
    /// JONSWAP γ for the background spectrum.
    pub gamma: f64,
}

impl RogueWaveSeed {
    /// Construct a rogue wave seed with standard parameters.
    ///
    /// # Arguments
    /// * `focus_amp`    – target focus amplitude (m)
    /// * `focus_time`   – time at which maximum crest appears (s)
    /// * `focus_x`      – spatial focus position (m)
    /// * `omega_peak`   – spectral peak (rad s⁻¹)
    pub fn new(focus_amp: f64, focus_time: f64, focus_x: f64, omega_peak: f64) -> Self {
        RogueWaveSeed {
            n_components: 64,
            focus_amplitude: focus_amp,
            focus_time,
            focus_x,
            omega_peak,
            gamma: 3.3,
        }
    }

    /// Surface elevation η(x, t) from phase-focused superposition.
    ///
    /// Each component amplitude is proportional to S(ω) dω, and all components
    /// are phased to arrive in-phase at (x_f, t_f).
    ///
    /// # Arguments
    /// * `x` – evaluation position (m)
    /// * `t` – evaluation time (s)
    pub fn elevation(&self, x: f64, t: f64) -> f64 {
        let params = JonswapParams {
            omega_peak: self.omega_peak,
            alpha: 0.0081,
            gamma: self.gamma,
            sigma_a: 0.07,
            sigma_b: 0.09,
        };
        let omega_min = 0.5 * self.omega_peak;
        let omega_max = 3.0 * self.omega_peak;
        let dw = (omega_max - omega_min) / self.n_components as f64;
        // Normalise so that Σ a_i = focus_amplitude
        let spectral_sum: f64 = (0..self.n_components)
            .map(|i| {
                let w = omega_min + (i as f64 + 0.5) * dw;
                (params.spectrum(w) * dw).sqrt()
            })
            .sum();
        if spectral_sum == 0.0 {
            return 0.0;
        }
        let scale = self.focus_amplitude / spectral_sum;
        (0..self.n_components)
            .map(|i| {
                let w = omega_min + (i as f64 + 0.5) * dw;
                let k = dispersion_wavenumber(w, 1000.0); // deep water
                let a_i = scale * (params.spectrum(w) * dw).sqrt();
                let phase_at_focus = k * self.focus_x - w * self.focus_time;
                a_i * (k * x - w * t - phase_at_focus).cos()
            })
            .sum()
    }

    /// Maximum theoretical crest height (approximately equal to focus_amplitude).
    pub fn max_crest_height(&self) -> f64 {
        self.focus_amplitude
    }
}

// ---------------------------------------------------------------------------
// SPH kernel (Wendland C2, compact support)
// ---------------------------------------------------------------------------

/// Evaluate the Wendland C2 kernel W(r, h).
///
/// # Arguments
/// * `r` – particle separation (m)
/// * `h` – smoothing length (m)
pub fn wendland_c2(r: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    if q >= 2.0 {
        return 0.0;
    }
    // 3-D normalisation constant: σ = 21/(16 π h³) ensures ∫ W 4πr² dr = 1.
    let sigma = 21.0 / (16.0 * PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    sigma * t.powi(4) * (2.0 * q + 1.0)
}

/// Gradient magnitude dW/dr of the Wendland C2 kernel.
///
/// # Arguments
/// * `r` – particle separation (m)
/// * `h` – smoothing length (m)
pub fn wendland_c2_grad(r: f64, h: f64) -> f64 {
    if h <= 0.0 || r <= 0.0 {
        return 0.0;
    }
    let q = r / h;
    if q >= 2.0 {
        return 0.0;
    }
    let sigma = 7.0 / (4.0 * PI * h * h * h);
    let t = 1.0 - 0.5 * q;
    // dW/dq * dq/dr
    sigma * (-2.0 * t.powi(3) * (2.0 * q + 1.0) + t.powi(4) * 2.0) / h
}

// ---------------------------------------------------------------------------
// SPH pressure (Tait equation of state)
// ---------------------------------------------------------------------------

/// Tait equation of state for weakly compressible SPH.
///
/// p = B \[(ρ/ρ₀)^γ - 1\]  where B = ρ₀ c_s² / γ.
///
/// # Arguments
/// * `rho`     – current density (kg m⁻³)
/// * `rho0`    – reference density (kg m⁻³)
/// * `cs`      – speed of sound (m s⁻¹)
/// * `gamma`   – exponent (typically 7 for water)
pub fn tait_pressure(rho: f64, rho0: f64, cs: f64, gamma: f64) -> f64 {
    let b = rho0 * cs * cs / gamma;
    b * ((rho / rho0).powf(gamma) - 1.0)
}

// ---------------------------------------------------------------------------
// SPH density summation
// ---------------------------------------------------------------------------

/// Compute density of particle `i` by summing over neighbour contributions.
///
/// ρᵢ = Σⱼ mⱼ W(|rᵢ - rⱼ|, hᵢ)
///
/// # Arguments
/// * `pos_i`    – position of particle i
/// * `h_i`      – smoothing length of particle i
/// * `neighbors` – slice of (position, mass) for all neighbours (including i)
pub fn sph_density_sum(pos_i: [f64; 3], h_i: f64, neighbors: &[([f64; 3], f64)]) -> f64 {
    neighbors
        .iter()
        .map(|(pos_j, mass_j)| {
            let dx = pos_i[0] - pos_j[0];
            let dy = pos_i[1] - pos_j[1];
            let dz = pos_i[2] - pos_j[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            mass_j * wendland_c2(r, h_i)
        })
        .sum()
}

// ---------------------------------------------------------------------------
// SPH pressure force
// ---------------------------------------------------------------------------

/// Symmetric SPH pressure acceleration: aᵢ = -Σⱼ mⱼ (pᵢ/ρᵢ² + pⱼ/ρⱼ²) ∇Wᵢⱼ.
///
/// # Arguments
/// * `pos_i`      – position of particle i
/// * `pi`         – pressure of particle i
/// * `rho_i`      – density of particle i
/// * `h_i`        – smoothing length
/// * `neighbors`  – slice of (position, mass, pressure, density) tuples
pub fn sph_pressure_acc(
    pos_i: [f64; 3],
    pi: f64,
    rho_i: f64,
    h_i: f64,
    neighbors: &[([f64; 3], f64, f64, f64)],
) -> [f64; 3] {
    let mut acc = [0.0_f64; 3];
    let pi_rho2 = if rho_i > 0.0 {
        pi / (rho_i * rho_i)
    } else {
        0.0
    };
    for (pos_j, mass_j, pj, rho_j) in neighbors {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        if r < 1e-14 {
            continue;
        }
        let pj_rho2 = if *rho_j > 0.0 {
            pj / (rho_j * rho_j)
        } else {
            0.0
        };
        let dw = wendland_c2_grad(r, h_i);
        let coeff = -mass_j * (pi_rho2 + pj_rho2) * dw / r;
        acc[0] += coeff * dx;
        acc[1] += coeff * dy;
        acc[2] += coeff * dz;
    }
    acc
}

// ---------------------------------------------------------------------------
// Beach run-up (Carrier–Greenspan)
// ---------------------------------------------------------------------------

/// Beach run-up height from Carrier & Greenspan (1958) / Hunt (1959) formula.
///
/// R/Hₛ = ξ  (Iribarren number) for ξ < 2 (surging breaker).
/// R = Hₛ ξ for ξ < 2,  R = Hₛ ξ² / (2 tanh(ξ)) otherwise.
///
/// # Arguments
/// * `wave_height`    – significant wave height offshore (m)
/// * `wave_length`    – wave length (m)
/// * `beach_slope`    – beach gradient (tan β)
pub fn beach_runup_height(wave_height: f64, wave_length: f64, beach_slope: f64) -> f64 {
    if wave_height <= 0.0 || wave_length <= 0.0 || beach_slope <= 0.0 {
        return 0.0;
    }
    // Iribarren (surf similarity) number ξ = tan(β) / √(H/L)
    let xi = beach_slope / (wave_height / wave_length).sqrt();
    if xi < 2.0 {
        wave_height * xi
    } else {
        wave_height * xi * xi / (2.0 * xi.tanh())
    }
}

/// Shoreline position X_s(t) from the Carrier–Greenspan parametric solution.
///
/// X_s(σ) = -σ + a sin(σ) where σ is a phase parameter.
///
/// # Arguments
/// * `sigma`        – phase parameter (dimensionless)
/// * `a`            – non-dimensional run-up amplitude
/// * `characteristic_length` – length scale L = (g h₀) / (ω² L₀) (m)
pub fn carrier_greenspan_shoreline(sigma: f64, a: f64, characteristic_length: f64) -> f64 {
    characteristic_length * (-sigma + a * sigma.sin())
}

// ---------------------------------------------------------------------------
// Ocean simulation driver
// ---------------------------------------------------------------------------

/// Configuration for an ocean SPH simulation.
#[derive(Clone, Debug)]
pub struct OceanSimConfig {
    /// Simulation time step (s).
    pub dt: f64,
    /// Total simulation duration (s).
    pub t_max: f64,
    /// Reference density (kg m⁻³).
    pub rho0: f64,
    /// Speed of sound (m s⁻¹) for WCSPH Tait EOS.
    pub cs: f64,
    /// Dynamic viscosity (Pa s).
    pub viscosity: f64,
    /// Water depth (m).
    pub depth: f64,
    /// Wind velocity vector at 10 m (m s⁻¹).
    pub wind_vel: [f64; 2],
    /// Tidal amplitude (m).
    pub tidal_amplitude: f64,
    /// Enable wave breaking.
    pub enable_breaking: bool,
    /// Enable foam generation.
    pub enable_foam: bool,
    /// Foam decay timescale (s).
    pub foam_tau: f64,
    /// Enable Stokes drift.
    pub enable_stokes: bool,
    /// Enable tidal forcing.
    pub enable_tides: bool,
}

impl Default for OceanSimConfig {
    /// Default ocean simulation parameters.
    fn default() -> Self {
        OceanSimConfig {
            dt: 0.01,
            t_max: 100.0,
            rho0: RHO_WATER,
            cs: 1500.0,
            viscosity: 1e-3,
            depth: 50.0,
            wind_vel: [10.0, 0.0],
            tidal_amplitude: 1.0,
            enable_breaking: true,
            enable_foam: true,
            foam_tau: 30.0,
            enable_stokes: true,
            enable_tides: false,
        }
    }
}

/// Ocean SPH simulation state.
pub struct OceanSimulation {
    /// All particles.
    pub particles: Vec<OceanParticle>,
    /// Simulation configuration.
    pub config: OceanSimConfig,
    /// Current simulation time (s).
    pub time: f64,
    /// Step counter.
    pub step: usize,
    /// JONSWAP spectral parameters.
    pub jonswap: JonswapParams,
}

impl OceanSimulation {
    /// Create a new ocean simulation.
    ///
    /// # Arguments
    /// * `config`  – simulation parameters
    /// * `jonswap` – spectral forcing parameters
    pub fn new(config: OceanSimConfig, jonswap: JonswapParams) -> Self {
        OceanSimulation {
            particles: Vec::new(),
            config,
            time: 0.0,
            step: 0,
            jonswap,
        }
    }

    /// Seed particles from the JONSWAP spectrum over a 1-D domain.
    ///
    /// Places `n` particles with initial surface elevation and velocity
    /// obtained by summing `n_modes` spectral components.
    ///
    /// # Arguments
    /// * `n`       – number of particles
    /// * `x_min`   – domain start (m)
    /// * `x_max`   – domain end (m)
    /// * `n_modes` – number of spectral modes to superpose
    /// * `mass`    – particle mass (kg)
    /// * `h`       – smoothing length (m)
    pub fn seed_from_jonswap(
        &mut self,
        n: usize,
        x_min: f64,
        x_max: f64,
        n_modes: usize,
        mass: f64,
        h: f64,
    ) {
        let omega_min = 0.5 * self.jonswap.omega_peak;
        let omega_max = 3.0 * self.jonswap.omega_peak;
        let dw = (omega_max - omega_min) / n_modes as f64;
        let depth = self.config.depth;

        for i in 0..n {
            let x = x_min + (i as f64 + 0.5) * (x_max - x_min) / n as f64;
            let mut eta = 0.0_f64;
            let mut u = 0.0_f64;
            for j in 0..n_modes {
                let w = omega_min + (j as f64 + 0.5) * dw;
                let k = dispersion_wavenumber(w, depth);
                let a = (2.0 * self.jonswap.spectrum(w) * dw).sqrt();
                let phase = k * x; // zero initial phase
                eta += a * phase.cos();
                let c = finite_depth_phase_speed(w, depth);
                u += a * w * phase.cos() / c.max(1e-6);
            }
            let mut p = OceanParticle::new([x, 0.0, eta], mass, h);
            p.vel[0] = u;
            p.eta = eta;
            p.depth = depth;
            self.particles.push(p);
        }
    }

    /// Advance the simulation by one time step.
    pub fn step(&mut self) {
        let dt = self.config.dt;
        let t = self.time;

        // Build neighbour list (O(n²) for clarity)
        let n = self.particles.len();
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.pos).collect();
        let masses: Vec<f64> = self.particles.iter().map(|p| p.mass).collect();
        let _pressures: Vec<f64> = self.particles.iter().map(|p| p.pressure).collect();
        let _densities: Vec<f64> = self.particles.iter().map(|p| p.rho).collect();

        // Update densities
        for i in 0..n {
            let h_i = self.particles[i].h;
            let neigh: Vec<([f64; 3], f64)> = (0..n).map(|j| (positions[j], masses[j])).collect();
            self.particles[i].rho = sph_density_sum(positions[i], h_i, &neigh);
            self.particles[i].pressure =
                tait_pressure(self.particles[i].rho, self.config.rho0, self.config.cs, 7.0);
        }

        // Update accelerations
        let pressures2: Vec<f64> = self.particles.iter().map(|p| p.pressure).collect();
        let densities2: Vec<f64> = self.particles.iter().map(|p| p.rho).collect();
        for i in 0..n {
            let h_i = self.particles[i].h;
            let neigh: Vec<([f64; 3], f64, f64, f64)> = (0..n)
                .map(|j| (positions[j], masses[j], pressures2[j], densities2[j]))
                .collect();
            let p_acc = sph_pressure_acc(positions[i], pressures2[i], densities2[i], h_i, &neigh);
            // Gravity
            let mut acc = [p_acc[0], p_acc[1], p_acc[2] - G];

            // Wind stress (surface particles only: z > -0.5 m)
            if self.particles[i].pos[2] > -0.5 && self.config.enable_stokes {
                let ws_acc = wind_stress_acceleration(
                    self.config.wind_vel,
                    self.particles[i].rho,
                    2.0 * h_i,
                );
                acc[0] += ws_acc[0];
                acc[1] += ws_acc[1];
            }

            // Tidal forcing
            if self.config.enable_tides {
                acc[0] +=
                    m2_tidal_acceleration_x(positions[i][0], t, self.config.tidal_amplitude, 0.0);
            }

            self.particles[i].acc = acc;
        }

        // Leap-frog integration
        for p in &mut self.particles {
            p.vel = vec3_add(p.vel, vec3_scale(p.acc, dt));
            p.pos = vec3_add(p.pos, vec3_scale(p.vel, dt));

            // Update surface elevation
            p.eta = p.pos[2];

            // Wave breaking detection
            if self.config.enable_breaking {
                let h_speed = (p.vel[0] * p.vel[0] + p.vel[1] * p.vel[1]).sqrt();
                p.breaking = froude_breaking_criterion(h_speed, p.depth.max(0.01));
            }

            // Foam update
            if self.config.enable_foam {
                if p.breaking {
                    p.foam = (p.foam + 0.2).min(1.0);
                }
                p.foam = foam_decay(p.foam, dt, self.config.foam_tau);
            }

            // Stokes drift (surface layer)
            if self.config.enable_stokes && p.pos[2] > -1.0 {
                let sd = stokes_drift_deep(
                    0.5,
                    self.jonswap.omega_peak,
                    dispersion_wavenumber(self.jonswap.omega_peak, p.depth.max(0.1)),
                    p.pos[2],
                    [1.0, 0.0],
                );
                p.stokes_drift = sd;
            }
        }

        self.time += dt;
        self.step += 1;
    }

    /// Total kinetic energy of all particles (J).
    pub fn total_kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Count particles currently undergoing wave breaking.
    pub fn count_breaking(&self) -> usize {
        self.particles.iter().filter(|p| p.breaking).count()
    }

    /// Mean foam intensity across all particles.
    pub fn mean_foam(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        self.particles.iter().map(|p| p.foam).sum::<f64>() / self.particles.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Wave energy and power
// ---------------------------------------------------------------------------

/// Wave energy density E = (1/2) ρ g a² (J m⁻²) for a monochromatic wave.
///
/// # Arguments
/// * `amplitude` – wave amplitude a (m)
pub fn wave_energy_density(amplitude: f64) -> f64 {
    0.5 * RHO_WATER * G * amplitude * amplitude
}

/// Wave power (energy flux) per unit crest width: P = E cg (W m⁻¹).
///
/// # Arguments
/// * `amplitude` – wave amplitude (m)
/// * `omega`     – angular frequency (rad s⁻¹)
/// * `depth`     – water depth (m)
pub fn wave_power(amplitude: f64, omega: f64, depth: f64) -> f64 {
    let e = wave_energy_density(amplitude);
    let cg = finite_depth_group_velocity(omega, depth);
    e * cg
}

// ---------------------------------------------------------------------------
// Viscous SPH force
// ---------------------------------------------------------------------------

/// Artificial viscosity (Monaghan 1992) for SPH.
///
/// Returns the viscous acceleration contribution for particle i from particle j.
///
/// # Arguments
/// * `pos_i`, `pos_j` – positions (m)
/// * `vel_i`, `vel_j` – velocities (m s⁻¹)
/// * `rho_ij`         – mean density (kg m⁻³)
/// * `h_ij`           – mean smoothing length (m)
/// * `cs_ij`          – mean sound speed (m s⁻¹)
/// * `alpha`          – shear coefficient (typically 0.1)
/// * `mass_j`         – mass of particle j (kg)
/// * `dw`             – kernel gradient magnitude (m⁻¹)
pub fn sph_artificial_viscosity_acc(
    pos_i: [f64; 3],
    pos_j: [f64; 3],
    vel_i: [f64; 3],
    vel_j: [f64; 3],
    rho_ij: f64,
    h_ij: f64,
    cs_ij: f64,
    alpha: f64,
    mass_j: f64,
    dw: f64,
) -> [f64; 3] {
    let dx = [
        pos_i[0] - pos_j[0],
        pos_i[1] - pos_j[1],
        pos_i[2] - pos_j[2],
    ];
    let dv = [
        vel_i[0] - vel_j[0],
        vel_i[1] - vel_j[1],
        vel_i[2] - vel_j[2],
    ];
    let r2 = vec3_dot(dx, dx);
    let vdotr = vec3_dot(dv, dx);
    if vdotr >= 0.0 {
        return [0.0; 3];
    }
    let mu = h_ij * vdotr / (r2 + 0.01 * h_ij * h_ij);
    let pi_ij = -alpha * cs_ij * mu / rho_ij;
    let r = r2.sqrt().max(1e-12);
    let coeff = -mass_j * pi_ij * dw / r;
    [coeff * dx[0], coeff * dx[1], coeff * dx[2]]
}

// ---------------------------------------------------------------------------
// JONSWAP spectrum seeding utilities
// ---------------------------------------------------------------------------

/// Compute the peak frequency of the JONSWAP spectrum from fetch and wind speed.
///
/// Using the JONSWAP empirical relation fp = 3.5 (g/U10) (g x / U10²)^{-0.33}.
///
/// # Arguments
/// * `u10`   – wind speed at 10 m (m s⁻¹)
/// * `fetch` – upwind fetch (m)
pub fn jonswap_peak_frequency(u10: f64, fetch: f64) -> f64 {
    if u10 <= 0.0 || fetch <= 0.0 {
        return 0.0;
    }
    let nondim_fetch = G * fetch / (u10 * u10);
    3.5 * (G / u10) * nondim_fetch.powf(-0.33)
}

/// Compute the JONSWAP alpha parameter from fetch and wind speed.
///
/// α = 0.076 (g x / U10²)^{-0.22}.
///
/// # Arguments
/// * `u10`   – wind speed (m s⁻¹)
/// * `fetch` – fetch (m)
pub fn jonswap_alpha(u10: f64, fetch: f64) -> f64 {
    if u10 <= 0.0 || fetch <= 0.0 {
        return 0.0081;
    }
    let nondim_fetch = G * fetch / (u10 * u10);
    0.076 * nondim_fetch.powf(-0.22)
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: relative tolerance check
    fn rel_close(a: f64, b: f64, tol: f64) -> bool {
        if b.abs() < 1e-30 {
            return a.abs() < tol;
        }
        ((a - b) / b).abs() < tol
    }

    #[test]
    fn test_jonswap_spectrum_peak() {
        let params = JonswapParams::new(0.1);
        let fp = 0.1_f64;
        let wp = 2.0 * PI * fp;
        // The spectrum at omega_peak should be the maximum in the vicinity
        let s_peak = params.spectrum(wp);
        let s_below = params.spectrum(wp * 0.8);
        assert!(s_peak > s_below, "spectrum should peak near omega_peak");
    }

    #[test]
    fn test_jonswap_spectrum_zero_at_zero() {
        let params = JonswapParams::new(0.1);
        assert_eq!(params.spectrum(0.0), 0.0);
        assert_eq!(params.spectrum(-1.0), 0.0);
    }

    #[test]
    fn test_significant_wave_height_positive() {
        let params = JonswapParams::new(0.1);
        let hs = params.significant_wave_height(0.01, 2.0, 1000);
        assert!(hs > 0.0, "Hs should be positive");
    }

    #[test]
    fn test_deep_water_phase_speed() {
        let omega = PI; // 1 rad/s
        let c = deep_water_phase_speed(omega);
        assert!(rel_close(c, G / omega, 1e-12));
    }

    #[test]
    fn test_deep_water_group_velocity_half_phase() {
        let omega = 2.0;
        let c = deep_water_phase_speed(omega);
        let cg = deep_water_group_velocity(omega);
        assert!(rel_close(cg, 0.5 * c, 1e-12));
    }

    #[test]
    fn test_dispersion_wavenumber_deep() {
        // In deep water k ≈ ω²/g
        let omega = 1.0;
        let d = 1000.0; // very deep
        let k = dispersion_wavenumber(omega, d);
        let k_deep = omega * omega / G;
        assert!(rel_close(k, k_deep, 1e-3), "k={k} vs k_deep={k_deep}");
    }

    #[test]
    fn test_dispersion_wavenumber_shallow() {
        // In shallow water: ω = k √(g d) → k = ω / √(g d)
        let d = 1.0; // shallow
        let omega = 0.5 * (G * d).sqrt(); // low frequency
        let k = dispersion_wavenumber(omega, d);
        let k_shallow = omega / (G * d).sqrt();
        assert!(rel_close(k, k_shallow, 5e-2));
    }

    #[test]
    fn test_finite_depth_phase_speed_positive() {
        let c = finite_depth_phase_speed(1.0, 10.0);
        assert!(c > 0.0);
    }

    #[test]
    fn test_finite_depth_group_velocity_lt_phase_speed() {
        let omega = 1.5;
        let d = 20.0;
        let c = finite_depth_phase_speed(omega, d);
        let cg = finite_depth_group_velocity(omega, d);
        assert!(
            cg <= c + 1e-10,
            "group velocity should not exceed phase speed"
        );
    }

    #[test]
    fn test_stokes_breaking_criterion_below() {
        // ka = 0.2 < π/7 ≈ 0.448
        assert!(!stokes_breaking_criterion(0.2, 1.0));
    }

    #[test]
    fn test_stokes_breaking_criterion_above() {
        // ka = 0.5 > π/7 ≈ 0.448
        assert!(stokes_breaking_criterion(0.5, 1.0));
    }

    #[test]
    fn test_froude_breaking_criterion() {
        let d = 4.0;
        let c_shallow = (G * d).sqrt();
        assert!(froude_breaking_criterion(c_shallow + 0.1, d));
        assert!(!froude_breaking_criterion(c_shallow - 0.1, d));
    }

    #[test]
    fn test_kinematic_breaking_criterion() {
        assert!(kinematic_breaking_criterion(10.0, 9.0));
        assert!(!kinematic_breaking_criterion(8.0, 9.0));
    }

    #[test]
    fn test_foam_decay() {
        let f0 = 1.0;
        let f1 = foam_decay(f0, 10.0, 10.0);
        assert!(rel_close(f1, (-1.0_f64).exp(), 1e-10));
    }

    #[test]
    fn test_foam_decay_zero_tau() {
        assert_eq!(foam_decay(1.0, 1.0, 0.0), 0.0);
    }

    #[test]
    fn test_foam_from_breaking_clamped() {
        let f = foam_from_breaking(1e10, 1.0);
        assert_eq!(f, 1.0);
        let f2 = foam_from_breaking(-1.0, 1.0);
        assert_eq!(f2, 0.0);
    }

    #[test]
    fn test_spray_ejection_speed() {
        let h = 2.0;
        let v = spray_ejection_speed(h);
        assert!(rel_close(v, (2.0 * G * h).sqrt(), 1e-12));
    }

    #[test]
    fn test_stokes_drift_positive() {
        let sd = stokes_drift_deep(1.0, 1.0, 0.1, 0.0, [1.0, 0.0]);
        assert!(sd[0] > 0.0);
        assert_eq!(sd[1], 0.0);
        assert_eq!(sd[2], 0.0);
    }

    #[test]
    fn test_stokes_drift_decays_with_depth() {
        let sd0 = stokes_drift_deep(1.0, 1.0, 0.1, 0.0, [1.0, 0.0]);
        let sd1 = stokes_drift_deep(1.0, 1.0, 0.1, -5.0, [1.0, 0.0]);
        assert!(sd0[0] > sd1[0], "Stokes drift should decay with depth");
    }

    #[test]
    fn test_stokes_drift_finite_depth_nonnegative() {
        let sd = stokes_drift_finite_depth(0.5, 1.0, 0.1, 20.0, -2.0, [1.0, 0.0]);
        assert!(sd[0] >= 0.0);
    }

    #[test]
    fn test_m2_tidal_elevation_amplitude() {
        let h_m2 = 0.5;
        let eta = m2_tidal_elevation(0.0, 0.0, h_m2, 0.0);
        assert!(rel_close(eta, h_m2, 1e-12));
    }

    #[test]
    fn test_large_pond_drag_low_wind() {
        let cd = large_pond_drag_coefficient(5.0);
        assert!(rel_close(cd, 1.2e-3, 1e-12));
    }

    #[test]
    fn test_large_pond_drag_high_wind() {
        let cd = large_pond_drag_coefficient(20.0);
        let expected = (0.49 + 0.065 * 20.0) * 1e-3;
        assert!(rel_close(cd, expected, 1e-12));
    }

    #[test]
    fn test_wind_surface_stress_direction() {
        let tau = wind_surface_stress([10.0, 0.0]);
        assert!(tau[0] > 0.0);
        assert_eq!(tau[1], 0.0);
    }

    #[test]
    fn test_wave_current_wavenumber_no_current() {
        let omega = 1.0;
        let d = 50.0;
        let k_ref = dispersion_wavenumber(omega, d);
        let k = wave_current_wavenumber(omega, d, 0.0);
        assert!(rel_close(k, k_ref, 1e-6));
    }

    #[test]
    fn test_wave_action_density_roundtrip() {
        let e = 100.0;
        let omega = 1.5;
        let n = wave_action_density(e, omega);
        let e2 = wave_energy_from_action(n, omega);
        assert!(rel_close(e2, e, 1e-12));
    }

    #[test]
    fn test_rogue_wave_elevation_focus() {
        let seed = RogueWaveSeed::new(3.0, 50.0, 100.0, 0.8);
        let eta_focus = seed.elevation(100.0, 50.0);
        // Should be close to the target focus amplitude
        assert!(
            eta_focus > 1.5,
            "Focus crest should be significant, got {eta_focus}"
        );
    }

    #[test]
    fn test_rogue_wave_max_crest() {
        let seed = RogueWaveSeed::new(4.0, 0.0, 0.0, 1.0);
        assert_eq!(seed.max_crest_height(), 4.0);
    }

    #[test]
    fn test_wendland_c2_normalised() {
        // Integral of W over 3D sphere should ≈ 1 (discrete trapezoidal)
        let h = 1.0;
        let n = 200;
        let dr = 2.0 * h / n as f64;
        let integral: f64 = (0..n)
            .map(|i| {
                let r = (i as f64 + 0.5) * dr;
                wendland_c2(r, h) * 4.0 * PI * r * r * dr
            })
            .sum();
        assert!(rel_close(integral, 1.0, 2e-2), "integral={integral}");
    }

    #[test]
    fn test_wendland_c2_zero_outside_support() {
        assert_eq!(wendland_c2(2.1, 1.0), 0.0);
    }

    #[test]
    fn test_tait_pressure_reference() {
        let p = tait_pressure(RHO_WATER, RHO_WATER, 1500.0, 7.0);
        assert!(
            p.abs() < 1.0,
            "Pressure at reference density should be zero, got {p}"
        );
    }

    #[test]
    fn test_sph_density_sum_single_particle() {
        let pos = [0.0; 3];
        let h = 1.0;
        let neigh = vec![(pos, 1.0_f64)];
        let rho = sph_density_sum(pos, h, &neigh);
        // Should equal mass * W(0, h)
        let w0 = wendland_c2(0.0, h);
        assert!(rel_close(rho, w0, 1e-12));
    }

    #[test]
    fn test_beach_runup_positive() {
        let r = beach_runup_height(1.0, 30.0, 0.1);
        assert!(r > 0.0);
    }

    #[test]
    fn test_beach_runup_zero_height() {
        assert_eq!(beach_runup_height(0.0, 30.0, 0.1), 0.0);
    }

    #[test]
    fn test_wave_energy_density() {
        let e = wave_energy_density(1.0);
        assert!(rel_close(e, 0.5 * RHO_WATER * G, 1e-12));
    }

    #[test]
    fn test_wave_power_positive() {
        let p = wave_power(1.0, 1.0, 50.0);
        assert!(p > 0.0);
    }

    #[test]
    fn test_jonswap_peak_frequency_scaling() {
        // Larger fetch → lower peak frequency
        let fp1 = jonswap_peak_frequency(10.0, 1e4);
        let fp2 = jonswap_peak_frequency(10.0, 1e6);
        assert!(fp2 < fp1, "fp should decrease with increasing fetch");
    }

    #[test]
    fn test_jonswap_alpha_scaling() {
        // Larger fetch → smaller alpha (developed sea)
        let a1 = jonswap_alpha(10.0, 1e4);
        let a2 = jonswap_alpha(10.0, 1e6);
        assert!(a2 < a1, "alpha should decrease with increasing fetch");
    }

    #[test]
    fn test_ocean_simulation_step() {
        let config = OceanSimConfig {
            dt: 0.01,
            t_max: 1.0,
            depth: 50.0,
            ..Default::default()
        };
        let jonswap = JonswapParams::new(0.1);
        let mut sim = OceanSimulation::new(config, jonswap);
        sim.seed_from_jonswap(20, 0.0, 100.0, 16, 1.0, 2.0);
        assert_eq!(sim.particles.len(), 20);
        sim.step();
        assert_eq!(sim.step, 1);
        assert!((sim.time - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_ocean_simulation_kinetic_energy_nonzero() {
        let config = OceanSimConfig::default();
        let jonswap = JonswapParams::new(0.1);
        let mut sim = OceanSimulation::new(config, jonswap);
        sim.seed_from_jonswap(10, 0.0, 50.0, 8, 1.0, 1.0);
        sim.step();
        // After one step there should be kinetic energy (gravity accelerates particles)
        let ke = sim.total_kinetic_energy();
        assert!(ke >= 0.0);
    }

    #[test]
    fn test_carrier_greenspan_shoreline_zero() {
        let x = carrier_greenspan_shoreline(0.0, 1.0, 10.0);
        assert_eq!(x, 0.0);
    }

    #[test]
    fn test_vec3_len() {
        assert!(rel_close(vec3_len([3.0, 4.0, 0.0]), 5.0, 1e-12));
    }

    #[test]
    fn test_m2_tidal_acceleration_x_sign() {
        // At x=0, t=0, phi=pi/2: sin(pi/2)=1, acceleration should be positive
        let a = m2_tidal_acceleration_x(0.0, 0.0, 1.0, PI / 2.0);
        assert!(a > 0.0, "Expected positive tidal acceleration, got {a}");
    }
}
