// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Volcanic flow simulation using SPH.
//!
//! Covers:
//! - Lava rheology with Bingham plastic model and temperature-dependent yield stress
//! - Pyroclastic density current (PDC) dynamics
//! - Tephra fallout with ballistic trajectories and drag
//! - Volcanic bomb aerodynamics
//! - Ash dispersion via advection-diffusion
//! - Caldera collapse dynamics
//! - Magma chamber pressurization
//! - Magma fragmentation threshold
//! - Volatile degassing and exsolution
//! - Lava tube formation and conduit flow

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Gravitational acceleration (m s⁻²).
const G: f64 = 9.81;

/// Universal gas constant (J mol⁻¹ K⁻¹).
const R_GAS: f64 = 8.314_462_618;

/// Stefan-Boltzmann constant (W m⁻² K⁻⁴).
const SIGMA_SB: f64 = 5.670_374_419e-8;

/// Atmospheric pressure (Pa).
const P_ATM: f64 = 101_325.0;

/// Reference temperature for rheology (K).
const T_REF: f64 = 1_273.15; // 1000 °C

// ---------------------------------------------------------------------------
// Vector math helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm of a 3-vector.
#[inline]
pub fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Add two 3-vectors.
#[inline]
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-vectors.
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector.
#[inline]
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Cross product of two 3-vectors.
#[inline]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// Lava rheology – Bingham plastic
// ---------------------------------------------------------------------------

/// Parameters describing the Bingham plastic rheology of lava.
pub struct BinghamLavaParams {
    /// Reference viscosity at `T_REF` (Pa s).
    pub eta_ref: f64,
    /// Activation energy for viscous flow (J mol⁻¹).
    pub activation_energy: f64,
    /// Yield stress at `T_REF` (Pa).
    pub yield_stress_ref: f64,
    /// Temperature exponent for yield stress decay.
    pub yield_temp_exponent: f64,
    /// Crystal volume fraction (0–1) that increases effective viscosity.
    pub crystal_fraction: f64,
}

impl BinghamLavaParams {
    /// Typical basaltic lava at 1200 K.
    pub fn basaltic() -> Self {
        Self {
            eta_ref: 1.0e3,
            activation_energy: 1.2e5,
            yield_stress_ref: 500.0,
            yield_temp_exponent: 2.5,
            crystal_fraction: 0.1,
        }
    }

    /// Rhyolitic lava – highly viscous.
    pub fn rhyolitic() -> Self {
        Self {
            eta_ref: 1.0e10,
            activation_energy: 2.5e5,
            yield_stress_ref: 1.0e5,
            yield_temp_exponent: 3.5,
            crystal_fraction: 0.3,
        }
    }
}

/// Compute temperature-dependent dynamic viscosity using Arrhenius law and
/// crystal fraction correction (Einstein-Roscoe model, Einstein 1906).
///
/// η(T) = η_ref · exp(E_a/R · (1/T − 1/T_ref)) · (1 − 1.35·φ)^{−2.5}
pub fn lava_viscosity(params: &BinghamLavaParams, temperature_k: f64) -> f64 {
    let arrhenius =
        (params.activation_energy / R_GAS * (1.0 / temperature_k.max(300.0) - 1.0 / T_REF)).exp();
    let crystal_corr = (1.0 - 1.35 * params.crystal_fraction.min(0.74)).powf(-2.5);
    params.eta_ref * arrhenius * crystal_corr
}

/// Compute temperature-dependent yield stress.
///
/// τ_y(T) = τ_{y,ref} · (T_ref / T)^n
pub fn lava_yield_stress(params: &BinghamLavaParams, temperature_k: f64) -> f64 {
    params.yield_stress_ref * (T_REF / temperature_k.max(300.0)).powf(params.yield_temp_exponent)
}

/// Evaluate the Bingham effective viscosity for a given strain rate magnitude.
///
/// η_eff = η + τ_y / (γ̇ + ε)   (regularised)
pub fn bingham_effective_viscosity(
    params: &BinghamLavaParams,
    temperature_k: f64,
    strain_rate: f64,
) -> f64 {
    let eta = lava_viscosity(params, temperature_k);
    let tau_y = lava_yield_stress(params, temperature_k);
    let eps = 1.0e-12;
    eta + tau_y / (strain_rate.abs() + eps)
}

/// Check whether lava is flowing (strain rate exceeds yielding criterion).
pub fn is_flowing(params: &BinghamLavaParams, temperature_k: f64, shear_stress: f64) -> bool {
    shear_stress > lava_yield_stress(params, temperature_k)
}

// ---------------------------------------------------------------------------
// SPH kernel – cubic spline (Monaghan 1992)
// ---------------------------------------------------------------------------

/// Compute the cubic-spline SPH kernel value.
///
/// `q = r / h` where `h` is the smoothing length.
pub fn cubic_kernel(q: f64, h: f64) -> f64 {
    let sigma = 10.0 / (7.0 * PI * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q * (1.0 - 0.5 * q))
    } else if q < 2.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}

/// Compute the gradient of the cubic-spline kernel in 1-D (dW/dr).
pub fn cubic_kernel_grad(q: f64, h: f64) -> f64 {
    let sigma = 10.0 / (7.0 * PI * h * h);
    let dw_dq = if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q)
    } else if q < 2.0 {
        sigma * (-0.75 * (2.0 - q).powi(2))
    } else {
        0.0
    };
    dw_dq / h
}

// ---------------------------------------------------------------------------
// SPH lava particle
// ---------------------------------------------------------------------------

/// State of a single SPH lava particle.
pub struct LavaParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Acceleration (m s⁻²).
    pub acceleration: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Density (kg m⁻³).
    pub density: f64,
    /// Pressure (Pa).
    pub pressure: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Crystal volume fraction.
    pub crystal_fraction: f64,
    /// Dissolved volatile mass fraction (H₂O, SO₂ …).
    pub volatile_content: f64,
    /// Whether this particle is in a lava tube.
    pub in_tube: bool,
}

impl LavaParticle {
    /// Create a particle with default basaltic properties.
    pub fn new_basaltic(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            acceleration: [0.0; 3],
            mass: 10.0,
            density: 2_700.0,
            pressure: P_ATM,
            temperature: 1_473.15, // 1200 °C
            h: 0.5,
            crystal_fraction: 0.05,
            volatile_content: 0.03,
            in_tube: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Weakly-compressible EOS for lava (Tait equation)
// ---------------------------------------------------------------------------

/// Parameters for the Tait equation of state for lava.
pub struct LavaEos {
    /// Reference density (kg m⁻³).
    pub rho0: f64,
    /// Reference pressure (Pa).
    pub p0: f64,
    /// Isentropic bulk modulus at reference state (Pa).
    pub kappa: f64,
    /// Tait exponent γ.
    pub gamma: f64,
}

impl LavaEos {
    /// Default basalt EOS.
    pub fn basalt() -> Self {
        Self {
            rho0: 2_700.0,
            p0: P_ATM,
            kappa: 1.0e10,
            gamma: 7.0,
        }
    }

    /// Compute pressure from density using Tait EOS.
    pub fn pressure(&self, rho: f64) -> f64 {
        self.p0 + self.kappa / self.gamma * ((rho / self.rho0).powf(self.gamma) - 1.0)
    }

    /// Compute speed of sound (m s⁻¹).
    pub fn speed_of_sound(&self, rho: f64) -> f64 {
        (self.kappa * (rho / self.rho0).powf(self.gamma - 1.0) / rho).sqrt()
    }
}

// ---------------------------------------------------------------------------
// Pyroclastic density current
// ---------------------------------------------------------------------------

/// Particle representing material in a pyroclastic density current.
pub struct PdcParticle {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Particle mass (kg).
    pub mass: f64,
    /// Bulk density including gas (kg m⁻³).
    pub density: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Gas volume fraction.
    pub gas_fraction: f64,
    /// Pyroclast grain size (m).
    pub grain_size: f64,
}

impl PdcParticle {
    /// Create a typical pyroclastic density current particle.
    pub fn new(position: [f64; 3], temperature_k: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass: 5.0,
            density: 50.0,
            temperature: temperature_k,
            gas_fraction: 0.8,
            grain_size: 1.0e-3,
        }
    }
}

/// Parameters controlling PDC dynamics.
pub struct PdcParams {
    /// Ambient air density (kg m⁻³).
    pub rho_air: f64,
    /// Drag coefficient for PDC cloud.
    pub cd: f64,
    /// Heat capacity of mixture (J kg⁻¹ K⁻¹).
    pub cp: f64,
    /// Radiative cooling emissivity.
    pub emissivity: f64,
    /// Ambient temperature (K).
    pub t_ambient: f64,
}

impl Default for PdcParams {
    fn default() -> Self {
        Self {
            rho_air: 1.2,
            cd: 0.44,
            cp: 1_200.0,
            emissivity: 0.95,
            t_ambient: 293.15,
        }
    }
}

/// Compute drag force on a PDC particle (opposing velocity).
pub fn pdc_drag_force(particle: &PdcParticle, params: &PdcParams) -> [f64; 3] {
    let r = particle.grain_size * 0.5;
    let area = PI * r * r;
    let v = particle.velocity;
    let v_mag = norm3(v);
    if v_mag < 1.0e-12 {
        return [0.0; 3];
    }
    let f_drag = 0.5 * params.rho_air * params.cd * area * v_mag * v_mag;
    scale3(scale3(v, 1.0 / v_mag), -f_drag)
}

/// Compute radiative cooling rate (W) for a PDC particle.
pub fn pdc_radiative_cooling(particle: &PdcParticle, params: &PdcParams) -> f64 {
    let r = particle.grain_size * 0.5;
    let area = 4.0 * PI * r * r;
    params.emissivity * SIGMA_SB * area * (particle.temperature.powi(4) - params.t_ambient.powi(4))
}

/// Integrate one PDC particle step with drag, gravity, and radiative cooling.
pub fn pdc_step(particle: &mut PdcParticle, params: &PdcParams, dt: f64) {
    let drag = pdc_drag_force(particle, params);
    let gravity = [0.0, 0.0, -G];

    // Acceleration = (gravity + drag/m)
    let acc = add3(gravity, scale3(drag, 1.0 / particle.mass.max(1.0e-30)));

    particle.velocity = add3(particle.velocity, scale3(acc, dt));
    particle.position = add3(particle.position, scale3(particle.velocity, dt));

    // Thermal cooling
    let q_rad = pdc_radiative_cooling(particle, params);
    let cp_mass = params.cp * particle.mass;
    particle.temperature -= q_rad * dt / cp_mass.max(1.0e-30);
    particle.temperature = particle.temperature.max(params.t_ambient);
}

// ---------------------------------------------------------------------------
// Tephra ballistics
// ---------------------------------------------------------------------------

/// A single tephra clast following a ballistic trajectory.
pub struct TephraClast {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Diameter (m).
    pub diameter: f64,
    /// Density of clast material (kg m⁻³).
    pub rho_clast: f64,
    /// Whether the clast has landed.
    pub landed: bool,
}

impl TephraClast {
    /// Create a new tephra clast with given properties.
    pub fn new(position: [f64; 3], velocity: [f64; 3], diameter: f64, rho_clast: f64) -> Self {
        let vol = PI / 6.0 * diameter.powi(3);
        Self {
            position,
            velocity,
            mass: vol * rho_clast,
            diameter,
            rho_clast,
            landed: false,
        }
    }
}

/// Parameters for tephra fallout environment.
pub struct TephraEnv {
    /// Air density (kg m⁻³).
    pub rho_air: f64,
    /// Wind vector (m s⁻¹).
    pub wind: [f64; 3],
    /// Drag coefficient (spherical default ~0.47).
    pub cd: f64,
}

impl Default for TephraEnv {
    fn default() -> Self {
        Self {
            rho_air: 1.2,
            wind: [5.0, 0.0, 0.0],
            cd: 0.47,
        }
    }
}

/// Compute drag force on a tephra clast relative to local wind.
pub fn tephra_drag(clast: &TephraClast, env: &TephraEnv) -> [f64; 3] {
    let rel_vel = sub3(clast.velocity, env.wind);
    let v_mag = norm3(rel_vel);
    if v_mag < 1.0e-12 {
        return [0.0; 3];
    }
    let area = PI * (clast.diameter * 0.5).powi(2);
    let f = 0.5 * env.rho_air * env.cd * area * v_mag * v_mag;
    scale3(scale3(rel_vel, 1.0 / v_mag), -f)
}

/// Ballistic range for a projectile launched at angle θ (no drag).
///
/// R = v² sin(2θ) / g
pub fn ballistic_range_no_drag(v0: f64, theta_rad: f64) -> f64 {
    v0 * v0 * (2.0 * theta_rad).sin() / G
}

/// Integrate one tephra clast step.
pub fn tephra_step(clast: &mut TephraClast, env: &TephraEnv, dt: f64) {
    if clast.landed {
        return;
    }
    let drag = tephra_drag(clast, env);
    let gravity = [0.0, 0.0, -G * clast.mass];
    let total_force = add3(gravity, drag);
    let acc = scale3(total_force, 1.0 / clast.mass.max(1.0e-30));

    clast.velocity = add3(clast.velocity, scale3(acc, dt));
    clast.position = add3(clast.position, scale3(clast.velocity, dt));

    if clast.position[2] <= 0.0 {
        clast.position[2] = 0.0;
        clast.velocity = [0.0; 3];
        clast.landed = true;
    }
}

// ---------------------------------------------------------------------------
// Volcanic bomb drag
// ---------------------------------------------------------------------------

/// A volcanic bomb (large ballistic fragment) in flight.
pub struct VolcanicBomb {
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Effective cross-sectional area (m²).
    pub area: f64,
    /// Drag coefficient (spinning/fusiform ~0.3).
    pub cd: f64,
    /// Surface temperature of the bomb (K).
    pub surface_temp: f64,
}

impl VolcanicBomb {
    /// Create a spherical bomb from diameter and density.
    pub fn spherical(diameter: f64, rho: f64, surface_temp: f64) -> Self {
        let r = diameter * 0.5;
        Self {
            position: [0.0; 3],
            velocity: [0.0; 3],
            mass: rho * 4.0 / 3.0 * PI * r.powi(3),
            area: PI * r * r,
            cd: 0.47,
            surface_temp,
        }
    }

    /// Compute aerodynamic drag force vector.
    pub fn drag_force(&self, rho_air: f64, wind: [f64; 3]) -> [f64; 3] {
        let rel = sub3(self.velocity, wind);
        let v_mag = norm3(rel);
        if v_mag < 1.0e-12 {
            return [0.0; 3];
        }
        let f = 0.5 * rho_air * self.cd * self.area * v_mag * v_mag;
        scale3(scale3(rel, 1.0 / v_mag), -f)
    }

    /// Integrate bomb trajectory one step.
    pub fn step(&mut self, rho_air: f64, wind: [f64; 3], dt: f64) {
        let drag = self.drag_force(rho_air, wind);
        let gravity_force = [0.0, 0.0, -G * self.mass];
        let total = add3(gravity_force, drag);
        let acc = scale3(total, 1.0 / self.mass.max(1.0e-30));
        self.velocity = add3(self.velocity, scale3(acc, dt));
        self.position = add3(self.position, scale3(self.velocity, dt));
    }
}

// ---------------------------------------------------------------------------
// Ash dispersion – advection-diffusion
// ---------------------------------------------------------------------------

/// State of an ash concentration field on a regular 2D grid.
pub struct AshField {
    /// Grid dimension in x.
    pub nx: usize,
    /// Grid dimension in y.
    pub ny: usize,
    /// Grid cell size (m).
    pub dx: f64,
    /// Ash concentration (kg m⁻³) stored row-major.
    pub concentration: Vec<f64>,
    /// Atmospheric diffusivity (m² s⁻¹).
    pub diffusivity: f64,
    /// Wind velocity (m s⁻¹) in x-direction.
    pub wind_x: f64,
    /// Wind velocity (m s⁻¹) in y-direction.
    pub wind_y: f64,
    /// Gravitational settling velocity (m s⁻¹).
    pub settling_velocity: f64,
}

impl AshField {
    /// Create a zero-concentration ash field.
    pub fn new(nx: usize, ny: usize, dx: f64, diffusivity: f64, wind_x: f64, wind_y: f64) -> Self {
        Self {
            nx,
            ny,
            dx,
            concentration: vec![0.0; nx * ny],
            diffusivity,
            wind_x,
            wind_y,
            settling_velocity: 0.01,
        }
    }

    /// Index helper.
    pub fn idx(&self, ix: usize, iy: usize) -> usize {
        iy * self.nx + ix
    }

    /// Add a point source of ash at grid cell (ix, iy).
    pub fn add_source(&mut self, ix: usize, iy: usize, amount: f64) {
        if ix < self.nx && iy < self.ny {
            let i = self.idx(ix, iy);
            self.concentration[i] += amount;
        }
    }

    /// Advance ash concentration one time step using FTCS advection-diffusion.
    pub fn step(&mut self, dt: f64) {
        let nx = self.nx;
        let ny = self.ny;
        let dx = self.dx;
        let d = self.diffusivity;
        let ux = self.wind_x;
        let uy = self.wind_y;
        let vs = self.settling_velocity;

        let mut next = self.concentration.clone();

        for iy in 1..ny.saturating_sub(1) {
            for ix in 1..nx.saturating_sub(1) {
                let c = self.concentration[self.idx(ix, iy)];
                let cx_p = self.concentration[self.idx(ix + 1, iy)];
                let cx_m = self.concentration[self.idx(ix - 1, iy)];
                let cy_p = self.concentration[self.idx(ix, iy + 1)];
                let cy_m = self.concentration[self.idx(ix, iy - 1)];

                // Diffusion (central difference)
                let diff = d * (cx_p - 2.0 * c + cx_m + cy_p - 2.0 * c + cy_m) / (dx * dx);

                // Advection (upwind)
                let adv_x = if ux >= 0.0 {
                    ux * (c - cx_m) / dx
                } else {
                    ux * (cx_p - c) / dx
                };
                let adv_y = if uy >= 0.0 {
                    uy * (c - cy_m) / dx
                } else {
                    uy * (cy_p - c) / dx
                };

                let settling_loss = vs * c / dx;

                let i = self.idx(ix, iy);
                next[i] = c + dt * (diff - adv_x - adv_y - settling_loss);
                next[i] = next[i].max(0.0);
            }
        }

        self.concentration = next;
    }

    /// Total ash mass in the domain (kg m⁻¹ in 2D).
    pub fn total_mass(&self) -> f64 {
        self.concentration.iter().sum::<f64>() * self.dx * self.dx
    }
}

// ---------------------------------------------------------------------------
// Caldera collapse dynamics
// ---------------------------------------------------------------------------

/// Parameters describing a collapsing caldera system.
pub struct CalderaParams {
    /// Caldera radius (m).
    pub radius: f64,
    /// Initial caldera depth (m).
    pub depth: f64,
    /// Chamber roof thickness (m).
    pub roof_thickness: f64,
    /// Rock density (kg m⁻³).
    pub rho_rock: f64,
    /// Young's modulus of roof rock (Pa).
    pub youngs_modulus: f64,
    /// Fracture toughness (Pa m^0.5).
    pub fracture_toughness: f64,
    /// Chamber overpressure (Pa).
    pub overpressure: f64,
}

impl CalderaParams {
    /// Typical caldera collapse parameters.
    pub fn typical() -> Self {
        Self {
            radius: 5_000.0,
            depth: 2_000.0,
            roof_thickness: 1_500.0,
            rho_rock: 2_800.0,
            youngs_modulus: 5.0e10,
            fracture_toughness: 2.0e6,
            overpressure: 5.0e7,
        }
    }
}

/// Compute the stress intensity factor K_I for circular caldera roof failure.
///
/// K_I ≈ σ · √(π a)  where σ is the bending stress and a is crack half-length.
pub fn caldera_stress_intensity(params: &CalderaParams) -> f64 {
    // Simple circular plate bending: σ ~ 3(3+ν)qR²/(8h²), ν≈0.25
    let nu = 0.25_f64;
    let q = params.overpressure;
    let r = params.radius;
    let h = params.roof_thickness;
    let sigma = 3.0 * (3.0 + nu) * q * r * r / (8.0 * h * h);
    let crack_half = h * 0.1; // initial crack = 10% of thickness
    sigma * (PI * crack_half).sqrt()
}

/// Check whether caldera collapse is imminent.
pub fn caldera_collapse_imminent(params: &CalderaParams) -> bool {
    caldera_stress_intensity(params) >= params.fracture_toughness
}

/// Estimate collapse duration using free-fall approximation (s).
pub fn caldera_collapse_time(params: &CalderaParams) -> f64 {
    (2.0 * params.depth / G).sqrt()
}

/// Estimate the volume flux during caldera collapse (m³ s⁻¹).
pub fn caldera_volume_flux(params: &CalderaParams) -> f64 {
    let area = PI * params.radius * params.radius;
    let depth = params.depth;
    let t_collapse = caldera_collapse_time(params);
    area * depth / t_collapse
}

// ---------------------------------------------------------------------------
// Magma chamber pressurization
// ---------------------------------------------------------------------------

/// State of a magma chamber.
pub struct MagmaChamber {
    /// Chamber volume (m³).
    pub volume: f64,
    /// Current magma pressure (Pa).
    pub pressure: f64,
    /// Temperature of magma (K).
    pub temperature: f64,
    /// Magma bulk modulus (Pa).
    pub bulk_modulus: f64,
    /// Country rock shear modulus (Pa).
    pub rock_shear_modulus: f64,
    /// Mass flux of new magma entering (kg s⁻¹).
    pub input_flux: f64,
    /// Magma density (kg m⁻³).
    pub magma_density: f64,
}

impl MagmaChamber {
    /// Create a chamber with typical arc volcano parameters.
    pub fn arc_volcano() -> Self {
        Self {
            volume: 1.0e12,
            pressure: 2.0e8,
            temperature: 1_173.15,
            bulk_modulus: 1.0e10,
            rock_shear_modulus: 3.0e10,
            input_flux: 1_000.0,
            magma_density: 2_600.0,
        }
    }

    /// Effective compressibility of chamber + country rock system.
    ///
    /// 1/β_eff = 1/β_magma + 1/β_rock  (Mogi sphere approximation)
    pub fn effective_compressibility(&self) -> f64 {
        let beta_magma = 1.0 / self.bulk_modulus;
        let beta_rock = 1.0 / (4.0 * self.rock_shear_modulus / 3.0);
        1.0 / (1.0 / beta_magma + 1.0 / beta_rock)
    }

    /// Pressurization rate dP/dt (Pa s⁻¹) from magma influx.
    pub fn pressurization_rate(&self) -> f64 {
        let beta_eff = self.effective_compressibility();
        // dP/dt = Q / (β_eff · V)   where Q is volumetric flux
        let vol_flux = self.input_flux / self.magma_density;
        vol_flux / (beta_eff * self.volume)
    }

    /// Advance chamber pressure one time step.
    pub fn step(&mut self, dt: f64) {
        self.pressure += self.pressurization_rate() * dt;
    }

    /// Compute surface uplift (m) using Mogi model.
    ///
    /// u_z = (ΔP · V) / (π G_rock r³)  for r >> chamber depth.
    pub fn mogi_uplift(&self, surface_distance: f64) -> f64 {
        let delta_p = self.pressure - P_ATM;
        let r = surface_distance.max(1.0);
        delta_p * self.volume / (PI * self.rock_shear_modulus * r * r * r)
    }
}

// ---------------------------------------------------------------------------
// Magma fragmentation
// ---------------------------------------------------------------------------

/// Fragmentation criteria for silicic magma (Dingwell 1996 / Papale 1999).
pub struct FragmentationParams {
    /// Strain rate threshold for fragmentation (s⁻¹).
    pub strain_rate_threshold: f64,
    /// Vesicularity threshold (gas volume fraction).
    pub vesicularity_threshold: f64,
    /// Tensile strength of melt (Pa).
    pub tensile_strength: f64,
    /// Relaxation time of melt (s).
    pub relaxation_time: f64,
}

impl FragmentationParams {
    /// Parameters for rhyolitic magma fragmentation.
    pub fn rhyolite() -> Self {
        Self {
            strain_rate_threshold: 0.01,
            vesicularity_threshold: 0.75,
            tensile_strength: 3.0e6,
            relaxation_time: 1.0e-3,
        }
    }
}

/// Determine whether fragmentation occurs by strain rate criterion.
///
/// Fragmentation when: γ̇ · τ_relax > 1 (brittle transition, Dingwell 1996).
pub fn fragmentation_by_strain_rate(params: &FragmentationParams, strain_rate: f64) -> bool {
    strain_rate * params.relaxation_time > 1.0
}

/// Determine whether fragmentation occurs by vesicularity threshold.
pub fn fragmentation_by_vesicularity(params: &FragmentationParams, phi_gas: f64) -> bool {
    phi_gas > params.vesicularity_threshold
}

/// Compute fragmentation pressure from bubble overpressure model.
///
/// ΔP_frag = 4 T / d_b  (Laplace pressure with T = surface tension)
pub fn bubble_fragmentation_pressure(surface_tension: f64, bubble_diameter: f64) -> f64 {
    4.0 * surface_tension / bubble_diameter.max(1.0e-12)
}

// ---------------------------------------------------------------------------
// Degassing – volatile exsolution
// ---------------------------------------------------------------------------

/// Parameters for H₂O exsolution in silicate melt (Liu et al. 2005 simplified).
pub struct DegassingParams {
    /// Solubility constant (Pa^{-0.5} for water).
    pub solubility_coeff: f64,
    /// Solubility exponent (0.5 for Henry's law with √P).
    pub solubility_exponent: f64,
    /// Molar mass of dissolved volatile (kg mol⁻¹).
    pub molar_mass: f64,
    /// Initial dissolved content (mass fraction).
    pub initial_content: f64,
}

impl DegassingParams {
    /// Typical water degassing parameters for rhyolite.
    pub fn water_rhyolite() -> Self {
        Self {
            solubility_coeff: 3.44e-6,
            solubility_exponent: 0.5,
            molar_mass: 0.018,
            initial_content: 0.04,
        }
    }
}

/// Equilibrium water solubility at given pressure (mass fraction).
pub fn water_solubility(params: &DegassingParams, pressure_pa: f64) -> f64 {
    params.solubility_coeff * pressure_pa.max(0.0).powf(params.solubility_exponent)
}

/// Exsolved gas fraction at given pressure.
pub fn exsolved_gas_fraction(params: &DegassingParams, pressure_pa: f64) -> f64 {
    let sol = water_solubility(params, pressure_pa);
    (params.initial_content - sol).max(0.0)
}

/// Compute bubble nucleation rate (m⁻³ s⁻¹) using classical nucleation theory.
///
/// J = A · exp(−16πγ³ / (3 k_B T (ΔP)²))
pub fn nucleation_rate(
    surface_tension: f64,
    temperature_k: f64,
    delta_p: f64,
    pre_factor: f64,
) -> f64 {
    let kb = 1.380_649e-23_f64;
    if delta_p.abs() < 1.0 {
        return 0.0;
    }
    let exponent =
        -16.0 * PI * surface_tension.powi(3) / (3.0 * kb * temperature_k * delta_p * delta_p);
    pre_factor * exponent.exp()
}

// ---------------------------------------------------------------------------
// Lava tube formation
// ---------------------------------------------------------------------------

/// Represents a lava tube conduit.
pub struct LavaTube {
    /// Start position of tube (m).
    pub start: [f64; 3],
    /// End position of tube (m).
    pub end: [f64; 3],
    /// Tube radius (m).
    pub radius: f64,
    /// Crust thickness insulating the tube (m).
    pub crust_thickness: f64,
    /// Internal lava temperature (K).
    pub lava_temp: f64,
    /// Whether the tube is active (flowing).
    pub active: bool,
}

impl LavaTube {
    /// Create a new lava tube segment.
    pub fn new(start: [f64; 3], end: [f64; 3], radius: f64) -> Self {
        Self {
            start,
            end,
            radius,
            crust_thickness: 0.1,
            lava_temp: 1_400.0,
            active: true,
        }
    }

    /// Length of the tube segment (m).
    pub fn length(&self) -> f64 {
        norm3(sub3(self.end, self.start))
    }

    /// Compute volume flux through tube using Hagen-Poiseuille (m³ s⁻¹).
    ///
    /// Q = π r⁴ ΔP / (8 η L)
    pub fn hagen_poiseuille_flux(&self, viscosity: f64, pressure_drop: f64) -> f64 {
        let len = self.length().max(1.0e-3);
        PI * self.radius.powi(4) * pressure_drop / (8.0 * viscosity * len)
    }

    /// Thermal conduction heat loss through crust (W).
    ///
    /// Q = k A (T_lava − T_surface) / d_crust
    pub fn heat_loss_conduction(&self, thermal_conductivity: f64, t_surface: f64) -> f64 {
        let area = 2.0 * PI * self.radius * self.length();
        thermal_conductivity * area * (self.lava_temp - t_surface)
            / self.crust_thickness.max(1.0e-6)
    }

    /// Update crust thickness growth over time step (Stefan solidification).
    ///
    /// d(crust) ≈ k_crust (T_lava − T_solidus) / (L_f ρ_crust) · dt / crust_thickness
    pub fn grow_crust(&mut self, k: f64, t_solidus: f64, latent_heat: f64, rho: f64, dt: f64) {
        let delta = k * (self.lava_temp - t_solidus).max(0.0)
            / (latent_heat * rho * self.crust_thickness.max(1.0e-6));
        self.crust_thickness += delta * dt;
    }
}

/// Determine whether a lava tube will form given surface cooling conditions.
///
/// Returns `true` when the crust solidification Peclet number Pe = v·L/κ > 1
/// indicating that advection dominates over diffusion (tube-forming regime).
pub fn lava_tube_formation_criterion(
    flow_velocity: f64,
    flow_length: f64,
    thermal_diffusivity: f64,
) -> bool {
    let pe = flow_velocity * flow_length / thermal_diffusivity.max(1.0e-12);
    pe > 1.0
}

// ---------------------------------------------------------------------------
// Volcanic eruption column (Plinian)
// ---------------------------------------------------------------------------

/// A 1D Plinian eruption column model.
pub struct PlumeSolver {
    /// Heights at which state is stored (m).
    pub heights: Vec<f64>,
    /// Vertical velocity at each level (m s⁻¹).
    pub velocity: Vec<f64>,
    /// Temperature at each level (K).
    pub temperature: Vec<f64>,
    /// Density at each level (kg m⁻³).
    pub density: Vec<f64>,
    /// Column radius at each level (m).
    pub radius: Vec<f64>,
    /// Entrainment coefficient.
    pub k_entrain: f64,
}

impl PlumeSolver {
    /// Create a plume solver with uniform initial conditions.
    pub fn new(n_levels: usize, dz: f64, v0: f64, t0: f64, rho0: f64, r0: f64) -> Self {
        let mut heights = Vec::with_capacity(n_levels);
        for i in 0..n_levels {
            heights.push(i as f64 * dz);
        }
        Self {
            heights,
            velocity: vec![v0; n_levels],
            temperature: vec![t0; n_levels],
            density: vec![rho0; n_levels],
            radius: vec![r0; n_levels],
            k_entrain: 0.09,
        }
    }

    /// Advance the plume one vertical integration step (simple upwind).
    pub fn integrate_step(&mut self, dz: f64, rho_atm: f64, t_atm: f64) {
        let n = self.heights.len();
        for i in 1..n {
            let v = self.velocity[i - 1];
            let t = self.temperature[i - 1];
            let rho = self.density[i - 1];
            let r = self.radius[i - 1];

            // Entrainment: Q_entr = 2 π r k v
            let mass_flux = rho * PI * r * r * v;
            let entr_flux = 2.0 * PI * r * self.k_entrain * v * rho_atm;
            let new_mass_flux = mass_flux + entr_flux * dz;

            // New radius from mass conservation
            let new_r = (new_mass_flux / (rho * PI * v.max(0.1))).sqrt();

            // Buoyancy term
            let buoyancy = G * (rho_atm - rho) / rho_atm;

            // Momentum
            let dv = (-G + buoyancy) / v.max(0.1);
            let new_v = (v + dv * dz).max(0.1);

            // Temperature mixing with entrained air
            let new_t = (t * mass_flux + t_atm * entr_flux * dz) / new_mass_flux.max(1.0e-30);

            self.velocity[i] = new_v;
            self.temperature[i] = new_t;
            self.radius[i] = new_r;
            self.density[i] = rho * v / new_v;
        }
    }
}

// ---------------------------------------------------------------------------
// SPH force computation for volcanic flow
// ---------------------------------------------------------------------------

/// Compute SPH pressure gradient force on particle i from particle j.
pub fn sph_pressure_force(
    pos_i: [f64; 3],
    pos_j: [f64; 3],
    p_i: f64,
    p_j: f64,
    rho_i: f64,
    rho_j: f64,
    mass_j: f64,
    h: f64,
) -> [f64; 3] {
    let r_ij = sub3(pos_i, pos_j);
    let r_mag = norm3(r_ij).max(1.0e-12);
    let q = r_mag / h;
    let dw = cubic_kernel_grad(q, h);
    let grad_w = scale3(r_ij, dw / r_mag);
    let pressure_term = p_i / (rho_i * rho_i) + p_j / (rho_j * rho_j);
    scale3(grad_w, -mass_j * pressure_term)
}

/// Compute SPH viscous force on particle i from particle j (Monaghan 1992).
pub fn sph_viscous_force(
    pos_i: [f64; 3],
    pos_j: [f64; 3],
    vel_i: [f64; 3],
    vel_j: [f64; 3],
    rho_i: f64,
    rho_j: f64,
    mass_j: f64,
    viscosity: f64,
    h: f64,
) -> [f64; 3] {
    let r_ij = sub3(pos_i, pos_j);
    let r_mag = norm3(r_ij).max(1.0e-12);
    let v_ij = sub3(vel_i, vel_j);
    let q = r_mag / h;
    let dw = cubic_kernel_grad(q, h);
    let pi_ij = 2.0 * viscosity / ((rho_i + rho_j) * 0.5 * r_mag * r_mag + 1.0e-6);
    let pi_term = dot3(v_ij, r_ij) * pi_ij;
    let grad_w = scale3(r_ij, dw / r_mag);
    scale3(grad_w, mass_j * pi_term)
}

/// Compute SPH density sum for a set of neighbor distances.
pub fn sph_density_sum(masses: &[f64], neighbor_q: &[f64], h: f64) -> f64 {
    masses
        .iter()
        .zip(neighbor_q.iter())
        .map(|(m, &q)| m * cubic_kernel(q, h))
        .sum()
}

// ---------------------------------------------------------------------------
// Thermal model for lava cooling
// ---------------------------------------------------------------------------

/// Compute conductive heat flux between two lava particles.
pub fn lava_conductive_heat_flux(t_i: f64, t_j: f64, r_ij: f64, thermal_conductivity: f64) -> f64 {
    thermal_conductivity * (t_j - t_i) / r_ij.max(1.0e-12)
}

/// Compute radiative heat loss from lava surface.
pub fn lava_radiative_cooling(temperature_k: f64, emissivity: f64, area: f64) -> f64 {
    emissivity * SIGMA_SB * area * temperature_k.powi(4)
}

/// Compute crystal fraction growth rate at given temperature.
///
/// dφ/dt ≈ A exp(−E_a / (R T)) for T < T_liquidus
pub fn crystal_growth_rate(temperature_k: f64, t_liquidus: f64, a: f64, ea: f64) -> f64 {
    if temperature_k >= t_liquidus {
        return 0.0;
    }
    a * (-ea / (R_GAS * temperature_k)).exp()
}

// ---------------------------------------------------------------------------
// Hazard utilities
// ---------------------------------------------------------------------------

/// Tephra accumulation at ground level from Gaussian plume model (kg m⁻²).
///
/// Σ = M / (2π σ_x σ_y) · exp(−x²/(2σ_x²) − y²/(2σ_y²))
pub fn tephra_accumulation_gaussian(
    total_mass: f64,
    x: f64,
    y: f64,
    sigma_x: f64,
    sigma_y: f64,
) -> f64 {
    let norm = 2.0 * PI * sigma_x * sigma_y;
    total_mass / norm * (-0.5 * (x * x / (sigma_x * sigma_x) + y * y / (sigma_y * sigma_y))).exp()
}

/// Estimate dynamic pressure of a PDC.
///
/// P_dyn = 0.5 ρ v²  (Pa)
pub fn pdc_dynamic_pressure(rho: f64, velocity_mag: f64) -> f64 {
    0.5 * rho * velocity_mag * velocity_mag
}

/// Convert PDC dynamic pressure to building damage level (qualitative).
///
/// Returns damage level 0 (none) through 5 (total destruction).
pub fn pdc_damage_level(dynamic_pressure_pa: f64) -> u8 {
    match dynamic_pressure_pa as u64 {
        0..=999 => 0,
        1_000..=4_999 => 1,
        5_000..=9_999 => 2,
        10_000..=24_999 => 3,
        25_000..=99_999 => 4,
        _ => 5,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Rheology tests ---

    #[test]
    fn test_lava_viscosity_arrhenius_increases_with_lower_temp() {
        let p = BinghamLavaParams::basaltic();
        let eta_hot = lava_viscosity(&p, 1500.0);
        let eta_cold = lava_viscosity(&p, 1000.0);
        assert!(eta_cold > eta_hot, "viscosity must increase as T drops");
    }

    #[test]
    fn test_lava_yield_stress_increases_with_lower_temp() {
        let p = BinghamLavaParams::basaltic();
        let tau_hot = lava_yield_stress(&p, 1500.0);
        let tau_cold = lava_yield_stress(&p, 900.0);
        assert!(tau_cold > tau_hot);
    }

    #[test]
    fn test_bingham_effective_viscosity_zero_strain_is_large() {
        let p = BinghamLavaParams::basaltic();
        let eta = bingham_effective_viscosity(&p, 1_200.0, 0.0);
        assert!(eta > 1.0e6);
    }

    #[test]
    fn test_bingham_effective_viscosity_high_strain_approaches_plastic() {
        let p = BinghamLavaParams::basaltic();
        let eta_low = bingham_effective_viscosity(&p, 1_200.0, 1.0e6);
        let eta_high = bingham_effective_viscosity(&p, 1_200.0, 1.0e-3);
        assert!(eta_low < eta_high);
    }

    #[test]
    fn test_is_flowing_above_yield() {
        let p = BinghamLavaParams::basaltic();
        let tau_y = lava_yield_stress(&p, 1_200.0);
        assert!(is_flowing(&p, 1_200.0, tau_y * 2.0));
        assert!(!is_flowing(&p, 1_200.0, tau_y * 0.5));
    }

    // --- Kernel tests ---

    #[test]
    fn test_cubic_kernel_positive_within_support() {
        let w = cubic_kernel(0.5, 1.0);
        assert!(w > 0.0);
    }

    #[test]
    fn test_cubic_kernel_zero_outside_support() {
        let w = cubic_kernel(2.5, 1.0);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_cubic_kernel_grad_negative_for_positive_q() {
        let dw = cubic_kernel_grad(0.5, 1.0);
        assert!(
            dw < 0.0,
            "kernel gradient should be negative for q in (0,1)"
        );
    }

    // --- EOS tests ---

    #[test]
    fn test_tait_eos_reference_pressure() {
        let eos = LavaEos::basalt();
        let p = eos.pressure(eos.rho0);
        assert!((p - eos.p0).abs() < 1.0e-6);
    }

    #[test]
    fn test_tait_eos_pressure_increases_with_density() {
        let eos = LavaEos::basalt();
        let p1 = eos.pressure(eos.rho0 * 1.01);
        let p2 = eos.pressure(eos.rho0);
        assert!(p1 > p2);
    }

    #[test]
    fn test_tait_eos_speed_of_sound_positive() {
        let eos = LavaEos::basalt();
        assert!(eos.speed_of_sound(eos.rho0) > 0.0);
    }

    // --- Tephra ballistics ---

    #[test]
    fn test_ballistic_range_45_degrees() {
        let r = ballistic_range_no_drag(100.0, PI / 4.0);
        let expected = 100.0_f64.powi(2) / G;
        assert!((r - expected).abs() < 1.0e-6);
    }

    #[test]
    fn test_tephra_step_lands() {
        let mut clast = TephraClast::new([0.0, 0.0, 10.0], [0.0, 0.0, -5.0], 0.01, 2_500.0);
        let env = TephraEnv::default();
        for _ in 0..200 {
            tephra_step(&mut clast, &env, 0.1);
        }
        assert!(clast.landed);
    }

    #[test]
    fn test_tephra_drag_opposes_motion() {
        let clast = TephraClast::new([0.0; 3], [10.0, 0.0, 0.0], 0.1, 2_500.0);
        let env = TephraEnv {
            wind: [0.0; 3],
            ..Default::default()
        };
        let f = tephra_drag(&clast, &env);
        assert!(f[0] < 0.0, "drag must oppose x-velocity");
    }

    // --- Volcanic bomb ---

    #[test]
    fn test_bomb_drag_zero_at_rest() {
        let bomb = VolcanicBomb::spherical(0.5, 2_800.0, 1_200.0);
        let f = bomb.drag_force(1.2, [0.0; 3]);
        assert_eq!(f, [0.0; 3]);
    }

    #[test]
    fn test_bomb_step_position_changes() {
        let mut bomb = VolcanicBomb::spherical(0.5, 2_800.0, 1_200.0);
        bomb.velocity = [50.0, 0.0, 50.0];
        let x0 = bomb.position[0];
        bomb.step(1.2, [0.0; 3], 0.1);
        assert!(bomb.position[0] > x0);
    }

    // --- PDC ---

    #[test]
    fn test_pdc_drag_force_opposes_velocity() {
        let mut p = PdcParticle::new([0.0; 3], 500.0);
        p.velocity = [20.0, 0.0, 0.0];
        let env = PdcParams::default();
        let f = pdc_drag_force(&p, &env);
        assert!(f[0] < 0.0);
    }

    #[test]
    fn test_pdc_radiative_cooling_positive() {
        let p = PdcParticle::new([0.0; 3], 1_000.0);
        let env = PdcParams::default();
        let q = pdc_radiative_cooling(&p, &env);
        assert!(q > 0.0);
    }

    #[test]
    fn test_pdc_step_temperature_decreases() {
        let mut p = PdcParticle::new([0.0; 3], 1_000.0);
        p.mass = 10.0;
        let env = PdcParams::default();
        let t0 = p.temperature;
        pdc_step(&mut p, &env, 1.0);
        assert!(p.temperature < t0);
    }

    // --- Ash dispersion ---

    #[test]
    fn test_ash_field_source_increases_concentration() {
        let mut field = AshField::new(10, 10, 1.0, 0.1, 0.0, 0.0);
        field.add_source(5, 5, 1.0);
        assert!(field.concentration[field.idx(5, 5)] > 0.0);
    }

    #[test]
    fn test_ash_field_total_mass_decreases_due_to_settling() {
        let mut field = AshField::new(20, 20, 1.0, 0.01, 0.0, 0.0);
        field.settling_velocity = 0.1;
        for ix in 5..15 {
            for iy in 5..15 {
                field.add_source(ix, iy, 1.0);
            }
        }
        let m0 = field.total_mass();
        for _ in 0..10 {
            field.step(0.1);
        }
        assert!(field.total_mass() < m0);
    }

    #[test]
    fn test_ash_field_concentration_nonnegative() {
        let mut field = AshField::new(10, 10, 1.0, 1.0, 5.0, 5.0);
        field.add_source(5, 5, 10.0);
        for _ in 0..50 {
            field.step(0.01);
        }
        assert!(field.concentration.iter().all(|&c| c >= 0.0));
    }

    // --- Caldera ---

    #[test]
    fn test_caldera_collapse_imminent_high_overpressure() {
        let p = CalderaParams::typical();
        assert!(caldera_collapse_imminent(&p));
    }

    #[test]
    fn test_caldera_stress_intensity_positive() {
        let p = CalderaParams::typical();
        assert!(caldera_stress_intensity(&p) > 0.0);
    }

    #[test]
    fn test_caldera_collapse_time_positive() {
        let p = CalderaParams::typical();
        assert!(caldera_collapse_time(&p) > 0.0);
    }

    #[test]
    fn test_caldera_volume_flux_positive() {
        let p = CalderaParams::typical();
        assert!(caldera_volume_flux(&p) > 0.0);
    }

    // --- Magma chamber ---

    #[test]
    fn test_chamber_pressurization_rate_positive() {
        let ch = MagmaChamber::arc_volcano();
        assert!(ch.pressurization_rate() > 0.0);
    }

    #[test]
    fn test_chamber_pressure_increases_after_step() {
        let mut ch = MagmaChamber::arc_volcano();
        let p0 = ch.pressure;
        ch.step(86_400.0); // 1 day
        assert!(ch.pressure > p0);
    }

    #[test]
    fn test_mogi_uplift_decreases_with_distance() {
        let ch = MagmaChamber::arc_volcano();
        let u1 = ch.mogi_uplift(5_000.0);
        let u2 = ch.mogi_uplift(10_000.0);
        assert!(u1 > u2);
    }

    // --- Fragmentation ---

    #[test]
    fn test_fragmentation_strain_rate_criterion() {
        let p = FragmentationParams::rhyolite();
        // relaxation_time = 1e-3 s, so strain_rate > 1000 s⁻¹ satisfies γ̇ · τ > 1
        assert!(fragmentation_by_strain_rate(&p, 2_000.0));
        assert!(!fragmentation_by_strain_rate(&p, 0.001));
    }

    #[test]
    fn test_fragmentation_vesicularity_criterion() {
        let p = FragmentationParams::rhyolite();
        assert!(fragmentation_by_vesicularity(&p, 0.8));
        assert!(!fragmentation_by_vesicularity(&p, 0.5));
    }

    #[test]
    fn test_bubble_fragmentation_pressure_positive() {
        let fp = bubble_fragmentation_pressure(0.1, 1.0e-3);
        assert!(fp > 0.0);
    }

    // --- Degassing ---

    #[test]
    fn test_water_solubility_increases_with_pressure() {
        let p = DegassingParams::water_rhyolite();
        let s1 = water_solubility(&p, 1.0e8);
        let s2 = water_solubility(&p, 1.0e7);
        assert!(s1 > s2);
    }

    #[test]
    fn test_exsolved_fraction_positive_below_saturation_pressure() {
        let p = DegassingParams::water_rhyolite();
        let frac = exsolved_gas_fraction(&p, 1.0e5);
        assert!(frac > 0.0);
    }

    #[test]
    fn test_exsolved_fraction_zero_at_high_pressure() {
        let p = DegassingParams::water_rhyolite();
        let frac = exsolved_gas_fraction(&p, 1.0e12);
        assert_eq!(frac, 0.0);
    }

    // --- Lava tube ---

    #[test]
    fn test_lava_tube_length_correct() {
        let tube = LavaTube::new([0.0; 3], [3.0, 4.0, 0.0], 2.0);
        assert!((tube.length() - 5.0).abs() < 1.0e-10);
    }

    #[test]
    fn test_hagen_poiseuille_flux_positive() {
        let tube = LavaTube::new([0.0; 3], [100.0, 0.0, 0.0], 1.0);
        let q = tube.hagen_poiseuille_flux(1.0e3, 1.0e4);
        assert!(q > 0.0);
    }

    #[test]
    fn test_heat_loss_conduction_positive() {
        let tube = LavaTube::new([0.0; 3], [100.0, 0.0, 0.0], 1.0);
        let q = tube.heat_loss_conduction(2.0, 300.0);
        assert!(q > 0.0);
    }

    #[test]
    fn test_lava_tube_formation_criterion_high_pe() {
        assert!(lava_tube_formation_criterion(10.0, 1_000.0, 1.0e-6));
    }

    // --- Hazard utilities ---

    #[test]
    fn test_tephra_accumulation_gaussian_peak_at_origin() {
        let c0 = tephra_accumulation_gaussian(1.0e9, 0.0, 0.0, 1_000.0, 1_000.0);
        let c1 = tephra_accumulation_gaussian(1.0e9, 2_000.0, 0.0, 1_000.0, 1_000.0);
        assert!(c0 > c1);
    }

    #[test]
    fn test_pdc_dynamic_pressure() {
        let p = pdc_dynamic_pressure(50.0, 100.0);
        assert!((p - 250_000.0).abs() < 1.0);
    }

    #[test]
    fn test_damage_level_total_destruction() {
        assert_eq!(pdc_damage_level(200_000.0), 5);
    }

    #[test]
    fn test_damage_level_none() {
        assert_eq!(pdc_damage_level(500.0), 0);
    }

    // --- SPH force ---

    #[test]
    fn test_sph_pressure_force_antisymmetric() {
        let pi = [0.0, 0.0, 0.0];
        let pj = [0.5, 0.0, 0.0];
        let f_ij = sph_pressure_force(pi, pj, 1.0e5, 1.0e5, 2_700.0, 2_700.0, 10.0, 1.0);
        let f_ji = sph_pressure_force(pj, pi, 1.0e5, 1.0e5, 2_700.0, 2_700.0, 10.0, 1.0);
        // Forces should be equal and opposite in x
        assert!((f_ij[0] + f_ji[0]).abs() < 1.0e-6 * f_ij[0].abs() + 1.0e-6);
    }

    #[test]
    fn test_sph_density_sum_positive() {
        let masses = vec![10.0_f64; 5];
        let qs = vec![0.3, 0.5, 0.8, 1.2, 1.8];
        let rho = sph_density_sum(&masses, &qs, 1.0);
        assert!(rho > 0.0);
    }

    // --- Crystal growth ---

    #[test]
    fn test_crystal_growth_rate_zero_above_liquidus() {
        let rate = crystal_growth_rate(1_600.0, 1_500.0, 1.0e10, 1.0e5);
        assert_eq!(rate, 0.0);
    }

    #[test]
    fn test_crystal_growth_rate_positive_below_liquidus() {
        let rate = crystal_growth_rate(1_200.0, 1_400.0, 1.0e10, 1.0e5);
        assert!(rate > 0.0);
    }

    // --- Plume solver ---

    #[test]
    fn test_plume_solver_velocity_array_length() {
        let solver = PlumeSolver::new(20, 100.0, 200.0, 1_273.0, 100.0, 50.0);
        assert_eq!(solver.velocity.len(), 20);
    }

    #[test]
    fn test_plume_solver_integrate_does_not_panic() {
        let mut solver = PlumeSolver::new(20, 100.0, 200.0, 1_273.0, 100.0, 50.0);
        solver.integrate_step(100.0, 1.2, 293.15);
    }
}
