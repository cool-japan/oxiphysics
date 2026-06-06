// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thin fluid film simulation for soft body physics.
//!
//! Provides shallow water equations on surfaces, film thickness evolution,
//! Marangoni flow, surface tension gradients, droplet spreading (Tanner's law),
//! thin film rupture, lubrication theory, rivulet flow, and coating uniformity.
//!
//! All quantities are in SI units (Pa, m, N, kg, s) unless otherwise noted.
//! Uses `[f64; 3]` arrays for vectors (no nalgebra dependency).

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector helpers (plain [f64; 3], no nalgebra)
// ---------------------------------------------------------------------------

/// Add two 3-vectors.
#[inline]
fn vec_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract b from a.
#[inline]
fn vec_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a vector.
#[inline]
fn vec_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
#[inline]
fn vec_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean length.
#[inline]
fn vec_len(a: [f64; 3]) -> f64 {
    vec_dot(a, a).sqrt()
}

/// Normalize, returning zero if degenerate.
#[inline]
fn vec_normalize(a: [f64; 3]) -> [f64; 3] {
    let l = vec_len(a);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        vec_scale(a, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Standard gravity (m/s²).
const G_ACCEL: f64 = 9.80665;
/// Water density at 20°C (kg/m³).
const RHO_WATER: f64 = 998.2;
/// Water dynamic viscosity at 20°C (Pa·s).
const MU_WATER: f64 = 1.002e-3;
/// Water surface tension at 20°C (N/m).
const GAMMA_WATER: f64 = 0.0728;

// ---------------------------------------------------------------------------
// Fluid film properties
// ---------------------------------------------------------------------------

/// Properties of a thin fluid film.
#[derive(Clone, Debug)]
pub struct FluidFilmProperties {
    /// Fluid density ρ (kg/m³).
    pub density: f64,
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Surface tension γ (N/m).
    pub surface_tension: f64,
    /// Temperature (K).
    pub temperature: f64,
    /// Contact angle θ (radians).
    pub contact_angle: f64,
}

impl FluidFilmProperties {
    /// Create a new fluid film properties set.
    pub fn new(
        density: f64,
        viscosity: f64,
        surface_tension: f64,
        temperature: f64,
        contact_angle: f64,
    ) -> Self {
        Self {
            density,
            viscosity,
            surface_tension,
            temperature,
            contact_angle,
        }
    }

    /// Water at 20°C with given contact angle.
    pub fn water(contact_angle: f64) -> Self {
        Self::new(RHO_WATER, MU_WATER, GAMMA_WATER, 293.15, contact_angle)
    }

    /// Silicone oil (polydimethylsiloxane) at 25°C, ~10 cSt.
    pub fn silicone_oil() -> Self {
        Self::new(930.0, 9.3e-3, 0.021, 298.15, 0.0)
    }

    /// Kinematic viscosity ν = μ/ρ (m²/s).
    pub fn kinematic_viscosity(&self) -> f64 {
        self.viscosity / self.density
    }

    /// Capillary length l_c = √(γ / (ρ g)) (m).
    pub fn capillary_length(&self) -> f64 {
        (self.surface_tension / (self.density * G_ACCEL)).sqrt()
    }

    /// Capillary number Ca = μ U / γ.
    pub fn capillary_number(&self, velocity: f64) -> f64 {
        self.viscosity * velocity / self.surface_tension
    }

    /// Weber number We = ρ U² L / γ.
    pub fn weber_number(&self, velocity: f64, length: f64) -> f64 {
        self.density * velocity * velocity * length / self.surface_tension
    }

    /// Reynolds number Re = ρ U L / μ.
    pub fn reynolds_number(&self, velocity: f64, length: f64) -> f64 {
        self.density * velocity * length / self.viscosity
    }

    /// Bond number Bo = ρ g L² / γ.
    pub fn bond_number(&self, length: f64) -> f64 {
        self.density * G_ACCEL * length * length / self.surface_tension
    }

    /// Ohnesorge number Oh = μ / √(ρ γ L).
    pub fn ohnesorge_number(&self, length: f64) -> f64 {
        self.viscosity / (self.density * self.surface_tension * length).sqrt()
    }

    /// Spreading coefficient S = γ_sv - γ_sl - γ_lv.
    ///
    /// Positive S means complete wetting (contact angle = 0).
    pub fn spreading_coefficient(gamma_sv: f64, gamma_sl: f64, gamma_lv: f64) -> f64 {
        gamma_sv - gamma_sl - gamma_lv
    }
}

// ---------------------------------------------------------------------------
// Shallow water equations on surfaces (1D)
// ---------------------------------------------------------------------------

/// State of a 1D shallow water cell.
#[derive(Clone, Debug)]
pub struct ShallowWaterCell1D {
    /// Film thickness h (m).
    pub h: f64,
    /// Depth-averaged velocity u (m/s).
    pub u: f64,
}

/// 1D shallow water solver on a surface.
#[derive(Clone, Debug)]
pub struct ShallowWaterSolver1D {
    /// Grid cells.
    pub cells: Vec<ShallowWaterCell1D>,
    /// Cell width Δx (m).
    pub dx: f64,
    /// Gravity component along the surface (m/s²).
    pub gravity_along: f64,
    /// Bottom friction coefficient.
    pub friction: f64,
}

impl ShallowWaterSolver1D {
    /// Create a new solver with `n` cells of width `dx`.
    pub fn new(n: usize, dx: f64, h0: f64, gravity_along: f64, friction: f64) -> Self {
        let cells = vec![ShallowWaterCell1D { h: h0, u: 0.0 }; n];
        Self {
            cells,
            dx,
            gravity_along,
            friction,
        }
    }

    /// Number of cells.
    pub fn num_cells(&self) -> usize {
        self.cells.len()
    }

    /// Total fluid volume per unit width (m²).
    pub fn total_volume(&self) -> f64 {
        self.cells.iter().map(|c| c.h * self.dx).sum()
    }

    /// Advance one time step using Lax-Friedrichs scheme.
    pub fn step(&mut self, dt: f64) {
        let n = self.cells.len();
        if n < 3 {
            return;
        }

        let mut new_h = vec![0.0; n];
        let mut new_hu = vec![0.0; n];

        // Interior cells: Lax-Friedrichs
        for i in 1..n - 1 {
            let hl = self.cells[i - 1].h;
            let hc = self.cells[i].h;
            let hr = self.cells[i + 1].h;
            let ul = self.cells[i - 1].u;
            let uc = self.cells[i].u;
            let ur = self.cells[i + 1].u;

            let hul = hl * ul;
            let _huc = hc * uc;
            let hur = hr * ur;

            // Flux for h: F_h = h u
            let fh_l = hul;
            let fh_r = hur;

            // Flux for hu: F_hu = h u² + 0.5 g h²
            let g = G_ACCEL;
            let fhu_l = hl * ul * ul + 0.5 * g * hl * hl;
            let fhu_r = hr * ur * ur + 0.5 * g * hr * hr;

            // Lax-Friedrichs
            new_h[i] = 0.5 * (hl + hr) - 0.5 * dt / self.dx * (fh_r - fh_l);
            new_hu[i] = 0.5 * (hul + hur) - 0.5 * dt / self.dx * (fhu_r - fhu_l);

            // Source terms
            new_hu[i] += dt * new_h[i] * self.gravity_along;
            // Bottom friction
            if new_h[i] > 1e-12 {
                let u_new = new_hu[i] / new_h[i];
                new_hu[i] -= dt * self.friction * u_new;
            }
        }

        // Boundary: reflective
        new_h[0] = new_h[1];
        new_hu[0] = -new_hu[1];
        new_h[n - 1] = new_h[n - 2];
        new_hu[n - 1] = -new_hu[n - 2];

        // Enforce non-negative thickness
        for i in 0..n {
            new_h[i] = new_h[i].max(0.0);
            self.cells[i].h = new_h[i];
            self.cells[i].u = if new_h[i] > 1e-12 {
                new_hu[i] / new_h[i]
            } else {
                0.0
            };
        }
    }

    /// CFL-stable time step estimate.
    pub fn cfl_dt(&self, cfl: f64) -> f64 {
        let mut max_speed = 1e-15;
        for c in &self.cells {
            let speed = c.u.abs() + (G_ACCEL * c.h).sqrt();
            if speed > max_speed {
                max_speed = speed;
            }
        }
        cfl * self.dx / max_speed
    }
}

// ---------------------------------------------------------------------------
// Film thickness evolution (diffusion equation)
// ---------------------------------------------------------------------------

/// Thin film thickness evolution under gravity and surface tension.
///
/// ∂h/∂t + ∂/∂x \[ (ρ g sin(α) h³ / 3μ) - (γ h³ / 3μ) ∂³h/∂x³ \] = 0
///
/// This is a 1D solver on a flat inclined surface.
#[derive(Clone, Debug)]
pub struct FilmThicknessEvolution1D {
    /// Film thickness at each grid point (m).
    pub h: Vec<f64>,
    /// Grid spacing (m).
    pub dx: f64,
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Surface tension (N/m).
    pub surface_tension: f64,
    /// Inclination angle (radians from horizontal).
    pub inclination: f64,
}

impl FilmThicknessEvolution1D {
    /// Create a new solver.
    pub fn new(
        n: usize,
        dx: f64,
        h0: f64,
        density: f64,
        viscosity: f64,
        surface_tension: f64,
        inclination: f64,
    ) -> Self {
        Self {
            h: vec![h0; n],
            dx,
            density,
            viscosity,
            surface_tension,
            inclination,
        }
    }

    /// Number of grid points.
    pub fn num_points(&self) -> usize {
        self.h.len()
    }

    /// Total film volume per unit width.
    pub fn total_volume(&self) -> f64 {
        self.h.iter().map(|hi| hi * self.dx).sum()
    }

    /// Gravity-driven flux coefficient: ρ g sin(α) / (3μ).
    pub fn gravity_flux_coeff(&self) -> f64 {
        self.density * G_ACCEL * self.inclination.sin() / (3.0 * self.viscosity)
    }

    /// Capillary flux coefficient: γ / (3μ).
    pub fn capillary_flux_coeff(&self) -> f64 {
        self.surface_tension / (3.0 * self.viscosity)
    }

    /// Advance one explicit time step.
    ///
    /// Uses finite differences for the thin film equation.
    pub fn step(&mut self, dt: f64) {
        let n = self.h.len();
        if n < 5 {
            return;
        }

        let grav_coeff = self.gravity_flux_coeff();
        let cap_coeff = self.capillary_flux_coeff();
        let dx = self.dx;
        let dx2 = dx * dx;
        let dx4 = dx2 * dx2;

        let old_h = self.h.clone();
        // Interior points
        for i in 2..n - 2 {
            let hi = old_h[i];
            if hi < 1e-15 {
                continue;
            }
            let h3 = hi * hi * hi;

            // Gravity-driven flux: ∂/∂x (h³)
            let dh3_dx = (old_h[i + 1].powi(3) - old_h[i - 1].powi(3)) / (2.0 * dx);

            // Capillary flux: ∂/∂x (h³ ∂³h/∂x³)
            // Fourth derivative approximation
            let d4h = if i >= 2 && i + 2 < n {
                (old_h[i + 2] - 4.0 * old_h[i + 1] + 6.0 * old_h[i] - 4.0 * old_h[i - 1]
                    + old_h[i - 2])
                    / dx4
            } else {
                0.0
            };

            self.h[i] = hi - dt * grav_coeff * dh3_dx - dt * cap_coeff * h3 * d4h;
            self.h[i] = self.h[i].max(0.0);
        }
    }

    /// Stable time step estimate (explicit scheme).
    pub fn stable_dt(&self) -> f64 {
        let h_max = self.h.iter().cloned().fold(0.0_f64, f64::max);
        if h_max < 1e-15 {
            return 1e-6;
        }
        let cap_coeff = self.capillary_flux_coeff();
        let dx4 = self.dx.powi(4);
        // Stability: dt < dx⁴ / (cap_coeff * h_max³ * some_factor)
        let h3 = h_max * h_max * h_max;
        dx4 / (16.0 * cap_coeff * h3 + 1e-30)
    }
}

// ---------------------------------------------------------------------------
// Marangoni flow
// ---------------------------------------------------------------------------

/// Marangoni flow model: surface tension gradient-driven flow.
///
/// The Marangoni stress τ_M = dγ/dT * dT/dx drives a shear at the free surface.
#[derive(Clone, Debug)]
pub struct MarangoniFlow {
    /// Surface tension temperature coefficient dγ/dT (N/m·K), typically negative.
    pub dgamma_dt: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Film thickness (m).
    pub film_thickness: f64,
}

impl MarangoniFlow {
    /// Create a new Marangoni flow model.
    pub fn new(dgamma_dt: f64, viscosity: f64, film_thickness: f64) -> Self {
        Self {
            dgamma_dt,
            viscosity,
            film_thickness,
        }
    }

    /// Marangoni stress for a given temperature gradient (Pa).
    ///
    /// τ_M = (dγ/dT) (dT/dx)
    pub fn marangoni_stress(&self, dt_dx: f64) -> f64 {
        self.dgamma_dt * dt_dx
    }

    /// Marangoni number Ma = -(dγ/dT) ΔT L / (μ α).
    ///
    /// Here `alpha` is thermal diffusivity (m²/s).
    pub fn marangoni_number(&self, delta_t: f64, length: f64, thermal_diff: f64) -> f64 {
        -(self.dgamma_dt) * delta_t * length / (self.viscosity * thermal_diff)
    }

    /// Free surface velocity driven by Marangoni stress (m/s).
    ///
    /// u_surface = τ_M h / (2μ) for linear velocity profile.
    pub fn surface_velocity(&self, dt_dx: f64) -> f64 {
        let tau = self.marangoni_stress(dt_dx);
        tau * self.film_thickness / (2.0 * self.viscosity)
    }

    /// Depth-averaged velocity (m/s).
    ///
    /// ū = τ_M h / (3μ) for parabolic profile with no-slip at substrate.
    pub fn average_velocity(&self, dt_dx: f64) -> f64 {
        let tau = self.marangoni_stress(dt_dx);
        tau * self.film_thickness / (3.0 * self.viscosity)
    }

    /// Volume flow rate per unit width Q = τ_M h² / (3μ) (m²/s).
    pub fn flow_rate(&self, dt_dx: f64) -> f64 {
        let tau = self.marangoni_stress(dt_dx);
        tau * self.film_thickness * self.film_thickness / (3.0 * self.viscosity)
    }

    /// Solutal Marangoni stress: τ = (dγ/dc) (dc/dx).
    pub fn solutal_stress(dgamma_dc: f64, dc_dx: f64) -> f64 {
        dgamma_dc * dc_dx
    }

    /// Thermocapillary velocity scale U_tc = |dγ/dT| ΔT / μ.
    pub fn thermocapillary_velocity(&self, delta_t: f64) -> f64 {
        self.dgamma_dt.abs() * delta_t / self.viscosity
    }
}

// ---------------------------------------------------------------------------
// Surface tension gradient model
// ---------------------------------------------------------------------------

/// Surface tension as a function of temperature: γ(T) = γ_0 + (dγ/dT)(T - T_0).
#[derive(Clone, Debug)]
pub struct SurfaceTensionModel {
    /// Reference surface tension γ_0 (N/m).
    pub gamma0: f64,
    /// Temperature coefficient dγ/dT (N/m·K).
    pub dgamma_dt: f64,
    /// Reference temperature T_0 (K).
    pub t_ref: f64,
}

impl SurfaceTensionModel {
    /// Create a new model.
    pub fn new(gamma0: f64, dgamma_dt: f64, t_ref: f64) -> Self {
        Self {
            gamma0,
            dgamma_dt,
            t_ref,
        }
    }

    /// Water surface tension model.
    pub fn water() -> Self {
        Self::new(GAMMA_WATER, -1.5e-4, 293.15)
    }

    /// Evaluate surface tension at temperature T (N/m).
    pub fn gamma(&self, t: f64) -> f64 {
        (self.gamma0 + self.dgamma_dt * (t - self.t_ref)).max(0.0)
    }

    /// Surface tension gradient at temperature T (N/m·K).
    pub fn gradient(&self) -> f64 {
        self.dgamma_dt
    }

    /// Young-Laplace pressure for a spherical droplet: ΔP = 2γ/R.
    pub fn young_laplace_sphere(&self, t: f64, radius: f64) -> f64 {
        2.0 * self.gamma(t) / radius
    }

    /// Young-Laplace pressure for a cylindrical film: ΔP = γ/R.
    pub fn young_laplace_cylinder(&self, t: f64, radius: f64) -> f64 {
        self.gamma(t) / radius
    }

    /// Capillary pressure for a thin film with curvature κ: ΔP = γ κ.
    pub fn capillary_pressure(&self, t: f64, curvature: f64) -> f64 {
        self.gamma(t) * curvature
    }
}

// ---------------------------------------------------------------------------
// Droplet spreading: Tanner's law
// ---------------------------------------------------------------------------

/// Tanner's law for droplet spreading on a flat surface.
///
/// R(t) ∝ t^(1/10)  for small contact angles.
/// θ(t) ∝ t^(-3/10)
#[derive(Clone, Debug)]
pub struct TannerSpreading {
    /// Droplet volume V (m³).
    pub volume: f64,
    /// Surface tension (N/m).
    pub surface_tension: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Initial contact radius (m).
    pub r0: f64,
}

impl TannerSpreading {
    /// Create a new Tanner spreading model.
    pub fn new(volume: f64, surface_tension: f64, viscosity: f64, r0: f64) -> Self {
        Self {
            volume,
            surface_tension,
            viscosity,
            r0,
        }
    }

    /// Characteristic spreading time scale τ = μ R0 / γ (s).
    pub fn time_scale(&self) -> f64 {
        self.viscosity * self.r0 / self.surface_tension
    }

    /// Contact radius at time t using Tanner's law: R(t) = R0 (t/τ)^(1/10).
    pub fn radius(&self, t: f64) -> f64 {
        let tau = self.time_scale();
        if tau < 1e-30 {
            return self.r0;
        }
        self.r0 * (t / tau).powf(0.1)
    }

    /// Contact angle at time t: θ(t) ∝ (t/τ)^(-3/10).
    ///
    /// Uses spherical cap approximation: θ ≈ (4V / (π R³))^(1/3) if small.
    pub fn contact_angle(&self, t: f64) -> f64 {
        let r = self.radius(t);
        if r < 1e-15 {
            return PI / 2.0;
        }
        // Spherical cap: V ≈ π R³ θ / 4 for small θ
        let theta = (4.0 * self.volume / (PI * r * r * r)).powf(1.0 / 3.0);
        theta.min(PI / 2.0)
    }

    /// Spreading velocity dR/dt (m/s).
    pub fn spreading_velocity(&self, t: f64) -> f64 {
        let tau = self.time_scale();
        if tau < 1e-30 || t < 1e-30 {
            return 0.0;
        }
        0.1 * self.r0 / tau * (t / tau).powf(-0.9)
    }

    /// Capillary-viscous spreading parameter (Tanner's coefficient).
    ///
    /// V_cl = γ θ³ / (6 μ ln(R/ε))  where ε is a microscopic cutoff.
    pub fn tanner_velocity(&self, theta: f64, radius: f64, epsilon: f64) -> f64 {
        if epsilon < 1e-30 || radius <= epsilon {
            return 0.0;
        }
        self.surface_tension * theta.powi(3) / (6.0 * self.viscosity * (radius / epsilon).ln())
    }
}

// ---------------------------------------------------------------------------
// Thin film rupture
// ---------------------------------------------------------------------------

/// Thin film rupture model based on van der Waals disjoining pressure.
///
/// Π(h) = -A / (6π h³)  where A is the Hamaker constant.
#[derive(Clone, Debug)]
pub struct ThinFilmRupture {
    /// Hamaker constant A (J).
    pub hamaker: f64,
    /// Surface tension (N/m).
    pub surface_tension: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Initial film thickness (m).
    pub h0: f64,
}

impl ThinFilmRupture {
    /// Create a new thin film rupture model.
    pub fn new(hamaker: f64, surface_tension: f64, viscosity: f64, h0: f64) -> Self {
        Self {
            hamaker,
            surface_tension,
            viscosity,
            h0,
        }
    }

    /// Disjoining pressure Π(h) = -A / (6π h³) (Pa).
    pub fn disjoining_pressure(&self, h: f64) -> f64 {
        if h < 1e-15 {
            return 0.0;
        }
        -self.hamaker / (6.0 * PI * h * h * h)
    }

    /// Critical rupture thickness h_c = (A² / (2πγ))^(1/4) × correction.
    ///
    /// Simplified: h_c ≈ (A / (2πγ))^(1/2) for spinodal dewetting threshold.
    pub fn critical_thickness(&self) -> f64 {
        (self.hamaker / (2.0 * PI * self.surface_tension)).sqrt()
    }

    /// Rupture time scale τ_r ≈ 12 π² μ γ h₀⁵ / A² for the most unstable mode.
    pub fn rupture_time(&self) -> f64 {
        12.0 * PI * PI * self.viscosity * self.surface_tension * self.h0.powi(5)
            / (self.hamaker * self.hamaker)
    }

    /// Most unstable wavelength λ_m = 2π h₀² √(2πγ / A).
    pub fn most_unstable_wavelength(&self) -> f64 {
        2.0 * PI * self.h0 * self.h0 * (2.0 * PI * self.surface_tension / self.hamaker).sqrt()
    }

    /// Growth rate of the most unstable mode (1/s).
    pub fn max_growth_rate(&self) -> f64 {
        let tau = self.rupture_time();
        if tau < 1e-30 {
            return 0.0;
        }
        1.0 / tau
    }

    /// Film stability: stable if h > h_c.
    pub fn is_stable(&self) -> bool {
        self.h0 > self.critical_thickness()
    }

    /// Spinodal dewetting condition: ∂Π/∂h > 0 → unstable.
    pub fn spinodal_parameter(&self, h: f64) -> f64 {
        // dΠ/dh = 3A / (6π h⁴) = A / (2π h⁴)
        if h < 1e-15 {
            return f64::INFINITY;
        }
        self.hamaker / (2.0 * PI * h.powi(4))
    }
}

// ---------------------------------------------------------------------------
// Lubrication theory
// ---------------------------------------------------------------------------

/// Lubrication (Reynolds) theory for thin film flow between surfaces.
#[derive(Clone, Debug)]
pub struct LubricationModel {
    /// Dynamic viscosity μ (Pa·s).
    pub viscosity: f64,
    /// Film thickness h (m).
    pub gap: f64,
    /// Length of the bearing/film region L (m).
    pub length: f64,
    /// Width of the bearing W (m).
    pub width: f64,
}

impl LubricationModel {
    /// Create a new lubrication model.
    pub fn new(viscosity: f64, gap: f64, length: f64, width: f64) -> Self {
        Self {
            viscosity,
            gap,
            length,
            width,
        }
    }

    /// Poiseuille flow rate per unit width Q = -h³/(12μ) dP/dx (m²/s).
    pub fn poiseuille_flow_rate(&self, dp_dx: f64) -> f64 {
        -self.gap.powi(3) / (12.0 * self.viscosity) * dp_dx
    }

    /// Couette flow rate per unit width Q = U h / 2 (m²/s).
    pub fn couette_flow_rate(&self, surface_velocity: f64) -> f64 {
        surface_velocity * self.gap / 2.0
    }

    /// Combined Couette-Poiseuille flow rate (m²/s).
    pub fn combined_flow_rate(&self, surface_velocity: f64, dp_dx: f64) -> f64 {
        self.couette_flow_rate(surface_velocity) + self.poiseuille_flow_rate(dp_dx)
    }

    /// Squeeze film force F = 3μ W L³ V / (2 h³) (N).
    ///
    /// `v_approach` is the approach velocity of the surfaces (m/s).
    pub fn squeeze_film_force(&self, v_approach: f64) -> f64 {
        3.0 * self.viscosity * self.width * self.length.powi(3) * v_approach
            / (2.0 * self.gap.powi(3))
    }

    /// Squeeze film damping coefficient c = 3μ W L³ / (2 h³) (N·s/m).
    pub fn squeeze_film_damping(&self) -> f64 {
        3.0 * self.viscosity * self.width * self.length.powi(3) / (2.0 * self.gap.powi(3))
    }

    /// Pressure in a slider bearing at position x.
    ///
    /// For a linear wedge h(x) = h1 + (h2-h1)x/L with surface velocity U:
    /// Returns the mid-plane pressure estimate (Pa).
    pub fn slider_bearing_pressure(&self, h1: f64, h2: f64, u_surface: f64) -> f64 {
        // Peak pressure for a linearly tapering gap
        let h_avg = (h1 + h2) / 2.0;
        6.0 * self.viscosity * u_surface * self.length * (h1 - h2) / (h_avg * h_avg * (h1 + h2))
    }

    /// Load capacity of a 1D slider bearing (N/m width).
    ///
    /// W = 6μ U L² / (h2² (K-1)²) \[ln(K) - 2(K-1)/(K+1)\]
    /// where K = h1/h2.
    pub fn slider_bearing_load(&self, h1: f64, h2: f64, u_surface: f64) -> f64 {
        if h2 < 1e-15 {
            return 0.0;
        }
        let k = h1 / h2;
        if (k - 1.0).abs() < 1e-10 {
            return 0.0; // Parallel surfaces, no load
        }
        let km1 = k - 1.0;
        let bracket = (k).ln() - 2.0 * km1 / (k + 1.0);
        6.0 * self.viscosity * u_surface * self.length * self.length * bracket
            / (h2 * h2 * km1 * km1)
    }

    /// Reynolds equation residual for a given pressure distribution.
    ///
    /// ∂/∂x(h³ ∂P/∂x) = 6μ U ∂h/∂x + 12μ ∂h/∂t
    pub fn reynolds_rhs(viscosity: f64, u_surface: f64, dh_dx: f64, dh_dt: f64) -> f64 {
        6.0 * viscosity * u_surface * dh_dx + 12.0 * viscosity * dh_dt
    }

    /// Sommerfeld number S = μ N L D / (W (c/R)²).
    ///
    /// For journal bearings. Returns dimensionless bearing parameter.
    pub fn sommerfeld_number(
        viscosity: f64,
        speed_rps: f64,
        load_per_area: f64,
        clearance_ratio: f64,
    ) -> f64 {
        viscosity * speed_rps / (load_per_area * clearance_ratio * clearance_ratio)
    }
}

// ---------------------------------------------------------------------------
// Rivulet flow
// ---------------------------------------------------------------------------

/// Rivulet: a thin stream flowing down an inclined surface.
#[derive(Clone, Debug)]
pub struct RivuletFlow {
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Surface tension (N/m).
    pub surface_tension: f64,
    /// Inclination angle (radians).
    pub inclination: f64,
    /// Contact angle (radians).
    pub contact_angle: f64,
    /// Volume flow rate Q (m³/s).
    pub flow_rate: f64,
}

impl RivuletFlow {
    /// Create a new rivulet flow model.
    pub fn new(
        density: f64,
        viscosity: f64,
        surface_tension: f64,
        inclination: f64,
        contact_angle: f64,
        flow_rate: f64,
    ) -> Self {
        Self {
            density,
            viscosity,
            surface_tension,
            inclination,
            contact_angle,
            flow_rate,
        }
    }

    /// Capillary length (m).
    pub fn capillary_length(&self) -> f64 {
        (self.surface_tension / (self.density * G_ACCEL)).sqrt()
    }

    /// Maximum rivulet width (m) from static equilibrium.
    ///
    /// w_max = 2 l_c sin(θ/2) where l_c is capillary length.
    pub fn max_width(&self) -> f64 {
        let lc = self.capillary_length();
        2.0 * lc * (self.contact_angle / 2.0).sin()
    }

    /// Maximum rivulet height (center) for a given width w (m).
    ///
    /// h_max ≈ w tan(θ) / 4  for small angles.
    pub fn max_height(&self, width: f64) -> f64 {
        width * self.contact_angle.tan() / 4.0
    }

    /// Mean velocity in the rivulet (m/s).
    ///
    /// ū ≈ ρ g sin(α) h² / (3μ) using film thickness ~ h_max.
    pub fn mean_velocity(&self) -> f64 {
        let w = self.max_width();
        let h = self.max_height(w);
        self.density * G_ACCEL * self.inclination.sin() * h * h / (3.0 * self.viscosity)
    }

    /// Cross-sectional area A ≈ w h / 2 (triangular approximation) (m²).
    pub fn cross_section_area(&self) -> f64 {
        let w = self.max_width();
        let h = self.max_height(w);
        w * h / 2.0
    }

    /// Rivulet Reynolds number Re = ρ Q / (μ w).
    pub fn reynolds_number(&self) -> f64 {
        let w = self.max_width();
        if w < 1e-15 {
            return 0.0;
        }
        self.density * self.flow_rate / (self.viscosity * w)
    }

    /// Meandering instability wavelength λ ≈ 2π w √(We).
    pub fn meander_wavelength(&self) -> f64 {
        let w = self.max_width();
        let u = self.mean_velocity();
        let we = self.density * u * u * w / self.surface_tension;
        2.0 * PI * w * we.sqrt()
    }
}

// ---------------------------------------------------------------------------
// Coating uniformity
// ---------------------------------------------------------------------------

/// Dip-coating model (Landau-Levich-Derjaguin).
///
/// Film thickness: h = 0.94 (μ U)^(2/3) / (γ^(1/6) (ρ g)^(1/2))
#[derive(Clone, Debug)]
pub struct DipCoating {
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Surface tension (N/m).
    pub surface_tension: f64,
    /// Withdrawal velocity (m/s).
    pub velocity: f64,
}

impl DipCoating {
    /// Create a new dip-coating model.
    pub fn new(density: f64, viscosity: f64, surface_tension: f64, velocity: f64) -> Self {
        Self {
            density,
            viscosity,
            surface_tension,
            velocity,
        }
    }

    /// Landau-Levich film thickness (m).
    pub fn film_thickness(&self) -> f64 {
        let ca = self.capillary_number();
        let lc = (self.surface_tension / (self.density * G_ACCEL)).sqrt();
        0.94 * lc * ca.powf(2.0 / 3.0)
    }

    /// Capillary number Ca = μ U / γ.
    pub fn capillary_number(&self) -> f64 {
        self.viscosity * self.velocity / self.surface_tension
    }

    /// Film thickness for large Ca (draining regime): h ≈ √(μ U / (ρ g)).
    pub fn draining_thickness(&self) -> f64 {
        (self.viscosity * self.velocity / (self.density * G_ACCEL)).sqrt()
    }

    /// Gravitational drainage: film thinning with time.
    ///
    /// h(t) = h0 / (1 + ρ g h0² t / (3μ))^(1/2) for vertical film.
    pub fn draining_film_thickness(&self, h0: f64, t: f64) -> f64 {
        let factor = 1.0 + self.density * G_ACCEL * h0 * h0 * t / (3.0 * self.viscosity);
        h0 / factor.sqrt()
    }

    /// Coating uniformity metric: fractional variation across a sample.
    ///
    /// For a gravity-draining film at positions x (heights from bottom),
    /// returns the maximum relative variation in thickness.
    pub fn thickness_variation(&self, h0: f64, t: f64, x_positions: &[f64]) -> f64 {
        if x_positions.is_empty() {
            return 0.0;
        }
        let thicknesses: Vec<f64> = x_positions
            .iter()
            .map(|&_x| self.draining_film_thickness(h0, t))
            .collect();
        let h_max = thicknesses.iter().cloned().fold(0.0_f64, f64::max);
        let h_min = thicknesses.iter().cloned().fold(f64::INFINITY, f64::min);
        if h_max < 1e-15 {
            return 0.0;
        }
        (h_max - h_min) / h_max
    }
}

/// Spin-coating model.
///
/// Film thickness: h(t) = h0 / √(1 + 4 ρ ω² h0² t / (3μ))
#[derive(Clone, Debug)]
pub struct SpinCoating {
    /// Fluid density (kg/m³).
    pub density: f64,
    /// Dynamic viscosity (Pa·s).
    pub viscosity: f64,
    /// Angular velocity ω (rad/s).
    pub omega: f64,
    /// Initial film thickness (m).
    pub h0: f64,
}

impl SpinCoating {
    /// Create a new spin-coating model.
    pub fn new(density: f64, viscosity: f64, omega: f64, h0: f64) -> Self {
        Self {
            density,
            viscosity,
            omega,
            h0,
        }
    }

    /// Film thickness at time t (m).
    pub fn thickness(&self, t: f64) -> f64 {
        let factor = 1.0
            + 4.0 * self.density * self.omega * self.omega * self.h0 * self.h0 * t
                / (3.0 * self.viscosity);
        self.h0 / factor.sqrt()
    }

    /// Final (equilibrium) thickness for evaporation-limited spin coating.
    ///
    /// h_f ≈ (3μ e / (2ρ ω²))^(1/3) where `e` is evaporation rate (m/s).
    pub fn final_thickness(&self, evaporation_rate: f64) -> f64 {
        (3.0 * self.viscosity * evaporation_rate / (2.0 * self.density * self.omega * self.omega))
            .powf(1.0 / 3.0)
    }

    /// Thinning rate dh/dt (m/s) at time t.
    pub fn thinning_rate(&self, t: f64) -> f64 {
        let h = self.thickness(t);
        -2.0 * self.density * self.omega * self.omega * h * h * h / (3.0 * self.viscosity)
    }

    /// Uniformity: radial thickness variation.
    ///
    /// In ideal spin coating, the thickness is uniform. Deviations come from
    /// edge effects. Returns the ratio of edge bead height to bulk.
    pub fn edge_bead_ratio(&self, radius: f64) -> f64 {
        let _h_bulk = self.thickness(1.0); // after 1 second
        let ca = self.viscosity * self.omega * radius / 0.072; // approximate γ
        // Edge bead is amplified at high Ca
        1.0 + ca.min(2.0) * 0.5
    }
}

// ---------------------------------------------------------------------------
// 2D shallow water on triangulated surface
// ---------------------------------------------------------------------------

/// A node in a 2D surface fluid film simulation.
#[derive(Clone, Debug)]
pub struct FilmNode {
    /// Position on the surface.
    pub position: [f64; 3],
    /// Surface normal at this node.
    pub normal: [f64; 3],
    /// Film thickness (m).
    pub thickness: f64,
    /// Depth-averaged velocity (m/s), tangent to surface.
    pub velocity: [f64; 3],
    /// Source/sink rate (m/s): positive = adding fluid.
    pub source_rate: f64,
}

impl FilmNode {
    /// Create a new film node.
    pub fn new(position: [f64; 3], normal: [f64; 3], thickness: f64) -> Self {
        Self {
            position,
            normal: vec_normalize(normal),
            thickness,
            velocity: [0.0; 3],
            source_rate: 0.0,
        }
    }

    /// Project a vector onto the surface tangent plane.
    pub fn project_tangent(&self, v: [f64; 3]) -> [f64; 3] {
        let d = vec_dot(v, self.normal);
        vec_sub(v, vec_scale(self.normal, d))
    }

    /// Gravity component tangent to the surface.
    pub fn tangent_gravity(&self, gravity: [f64; 3]) -> [f64; 3] {
        self.project_tangent(gravity)
    }
}

/// Simple 2D fluid film simulation on a mesh.
#[derive(Clone, Debug)]
pub struct FluidFilmSimulation {
    /// Nodes with film data.
    pub nodes: Vec<FilmNode>,
    /// Edges as pairs of node indices.
    pub edges: Vec<(usize, usize)>,
    /// Fluid properties.
    pub properties: FluidFilmProperties,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Minimum film thickness (numerical floor).
    pub h_min: f64,
}

impl FluidFilmSimulation {
    /// Create a new fluid film simulation.
    pub fn new(
        nodes: Vec<FilmNode>,
        edges: Vec<(usize, usize)>,
        properties: FluidFilmProperties,
        gravity: [f64; 3],
    ) -> Self {
        Self {
            nodes,
            edges,
            properties,
            gravity,
            h_min: 1e-10,
        }
    }

    /// Total fluid volume (m³), assuming unit area per node.
    pub fn total_volume(&self) -> f64 {
        self.nodes.iter().map(|n| n.thickness).sum()
    }

    /// Average film thickness (m).
    pub fn average_thickness(&self) -> f64 {
        let n = self.nodes.len() as f64;
        if n < 1.0 {
            return 0.0;
        }
        self.total_volume() / n
    }

    /// Maximum film thickness (m).
    pub fn max_thickness(&self) -> f64 {
        self.nodes
            .iter()
            .map(|n| n.thickness)
            .fold(0.0_f64, f64::max)
    }

    /// Advance the simulation by one time step.
    ///
    /// Uses a simple diffusion-like scheme along edges.
    pub fn step(&mut self, dt: f64) {
        let n = self.nodes.len();
        let mut dh = vec![0.0_f64; n];

        let rho = self.properties.density;
        let mu = self.properties.viscosity;
        let _gamma = self.properties.surface_tension;

        for &(i, j) in &self.edges {
            let dx = vec_sub(self.nodes[j].position, self.nodes[i].position);
            let dist = vec_len(dx);
            if dist < 1e-15 {
                continue;
            }

            let hi = self.nodes[i].thickness;
            let hj = self.nodes[j].thickness;
            let h_avg = 0.5 * (hi + hj);

            // Gravity-driven flux along edge
            let g_tang_i = self.nodes[i].tangent_gravity(self.gravity);
            let g_tang_j = self.nodes[j].tangent_gravity(self.gravity);
            let g_avg = vec_scale(vec_add(g_tang_i, g_tang_j), 0.5);
            let g_along = vec_dot(g_avg, vec_normalize(dx));

            // Film equation flux: Q = ρ g h³ / (3μ) per unit width
            let flux_grav = rho * g_along * h_avg.powi(3) / (3.0 * mu);

            // Thickness-gradient driven flux (leveling)
            let dh_dx = (hj - hi) / dist;
            let flux_level = rho * G_ACCEL * h_avg.powi(3) * dh_dx / (3.0 * mu);

            let total_flux = (flux_grav + flux_level) * dt / dist;

            dh[i] += total_flux;
            dh[j] -= total_flux;
        }

        // Apply source terms and update
        for (i, &dhi) in dh.iter().enumerate().take(n) {
            self.nodes[i].thickness += dhi + self.nodes[i].source_rate * dt;
            self.nodes[i].thickness = self.nodes[i].thickness.max(self.h_min);
        }
    }

    /// Uniformity metric: coefficient of variation of thickness.
    pub fn uniformity_cv(&self) -> f64 {
        let avg = self.average_thickness();
        if avg < 1e-15 {
            return 0.0;
        }
        let variance: f64 = self
            .nodes
            .iter()
            .map(|n| {
                let d = n.thickness - avg;
                d * d
            })
            .sum::<f64>()
            / self.nodes.len() as f64;
        variance.sqrt() / avg
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    // 1. Fluid properties: capillary length of water
    #[test]
    fn test_water_capillary_length() {
        let fp = FluidFilmProperties::water(PI / 6.0);
        let lc = fp.capillary_length();
        // Water capillary length ≈ 2.7 mm
        assert!(
            (lc - 2.7e-3).abs() < 0.5e-3,
            "Water capillary length should be ~2.7 mm, got {:.6e}",
            lc
        );
    }

    // 2. Kinematic viscosity
    #[test]
    fn test_kinematic_viscosity() {
        let fp = FluidFilmProperties::water(0.0);
        let nu = fp.kinematic_viscosity();
        assert!(
            (nu - 1.0e-6).abs() < 0.1e-6,
            "Water kinematic viscosity ~1e-6 m²/s, got {:.6e}",
            nu
        );
    }

    // 3. Capillary number
    #[test]
    fn test_capillary_number() {
        let fp = FluidFilmProperties::water(0.0);
        let ca = fp.capillary_number(0.01);
        assert!(ca > 0.0 && ca < 1.0, "Ca should be small for slow flow");
    }

    // 4. Shallow water solver conserves volume
    #[test]
    fn test_shallow_water_conservation() {
        let mut sw = ShallowWaterSolver1D::new(50, 0.01, 0.001, 0.0, 0.0);
        // Add a bump
        sw.cells[25].h = 0.005;
        let v0 = sw.total_volume();
        let dt = sw.cfl_dt(0.5);
        for _ in 0..10 {
            sw.step(dt);
        }
        let v1 = sw.total_volume();
        // Volume should be approximately conserved (reflective BC)
        assert!(
            (v1 - v0).abs() / v0 < 0.1,
            "Volume change too large: {:.6e} -> {:.6e}",
            v0,
            v1
        );
    }

    // 5. CFL dt is positive
    #[test]
    fn test_cfl_dt_positive() {
        let sw = ShallowWaterSolver1D::new(50, 0.01, 0.001, 0.1, 0.0);
        let dt = sw.cfl_dt(0.5);
        assert!(dt > 0.0, "CFL dt must be positive");
    }

    // 6. Film thickness evolution: stable dt positive
    #[test]
    fn test_film_evolution_stable_dt() {
        let fte =
            FilmThicknessEvolution1D::new(100, 0.001, 1e-4, RHO_WATER, MU_WATER, GAMMA_WATER, 0.1);
        let dt = fte.stable_dt();
        assert!(dt > 0.0, "Stable dt must be positive");
    }

    // 7. Film thickness evolution conserves volume (flat film, no inclination)
    #[test]
    fn test_film_evolution_flat_conservation() {
        let mut fte =
            FilmThicknessEvolution1D::new(50, 0.001, 1e-4, RHO_WATER, MU_WATER, GAMMA_WATER, 0.0);
        let v0 = fte.total_volume();
        let dt = fte.stable_dt().min(1e-6);
        for _ in 0..5 {
            fte.step(dt);
        }
        let v1 = fte.total_volume();
        assert!(
            (v1 - v0).abs() / v0 < 0.01,
            "Volume should be conserved for flat film"
        );
    }

    // 8. Marangoni stress proportional to temperature gradient
    #[test]
    fn test_marangoni_stress_linear() {
        let mg = MarangoniFlow::new(-1.5e-4, 1e-3, 1e-4);
        let tau1 = mg.marangoni_stress(100.0);
        let tau2 = mg.marangoni_stress(200.0);
        assert!(
            (tau2 - 2.0 * tau1).abs() < TOL,
            "Marangoni stress should be linear in dT/dx"
        );
    }

    // 9. Marangoni number positive for negative dγ/dT and positive ΔT
    #[test]
    fn test_marangoni_number_sign() {
        let mg = MarangoniFlow::new(-1.5e-4, 1e-3, 1e-4);
        let ma = mg.marangoni_number(10.0, 0.01, 1e-7);
        assert!(ma > 0.0, "Ma should be positive, got {:.6}", ma);
    }

    // 10. Surface tension decreases with temperature (water)
    #[test]
    fn test_surface_tension_temperature() {
        let stm = SurfaceTensionModel::water();
        let g1 = stm.gamma(293.15);
        let g2 = stm.gamma(373.15);
        assert!(g2 < g1, "Surface tension should decrease with temperature");
    }

    // 11. Young-Laplace pressure positive for sphere
    #[test]
    fn test_young_laplace() {
        let stm = SurfaceTensionModel::water();
        let dp = stm.young_laplace_sphere(293.15, 1e-3);
        assert!(dp > 0.0, "Young-Laplace pressure must be positive");
        // ΔP = 2γ/R ≈ 2*0.0728/0.001 ≈ 145.6 Pa
        assert!(
            (dp - 145.6).abs() < 1.0,
            "ΔP should be ~146 Pa, got {:.6}",
            dp
        );
    }

    // 12. Tanner's law: radius increases with time
    #[test]
    fn test_tanner_radius_increases() {
        let ts = TannerSpreading::new(1e-9, 0.072, 1e-3, 1e-3);
        let r1 = ts.radius(1.0);
        let r2 = ts.radius(10.0);
        assert!(r2 > r1, "Droplet should spread with time");
    }

    // 13. Tanner spreading velocity positive at early times
    #[test]
    fn test_tanner_velocity_positive() {
        let ts = TannerSpreading::new(1e-9, 0.072, 1e-3, 1e-3);
        let v = ts.spreading_velocity(1.0);
        assert!(v > 0.0, "Spreading velocity must be positive");
    }

    // 14. Disjoining pressure sign
    #[test]
    fn test_disjoining_pressure() {
        let tfr = ThinFilmRupture::new(1e-20, 0.072, 1e-3, 100e-9);
        let pi_h = tfr.disjoining_pressure(100e-9);
        // Hamaker > 0, so Π = -A/(6πh³) < 0 (attractive → destabilizing)
        assert!(pi_h < 0.0, "Disjoining pressure should be negative for A>0");
    }

    // 15. Rupture time positive
    #[test]
    fn test_rupture_time_positive() {
        let tfr = ThinFilmRupture::new(1e-20, 0.072, 1e-3, 100e-9);
        let tau = tfr.rupture_time();
        assert!(tau > 0.0, "Rupture time must be positive");
    }

    // 16. Most unstable wavelength positive
    #[test]
    fn test_unstable_wavelength() {
        let tfr = ThinFilmRupture::new(1e-20, 0.072, 1e-3, 100e-9);
        let lam = tfr.most_unstable_wavelength();
        assert!(lam > 0.0, "Wavelength must be positive");
    }

    // 17. Couette flow rate = U h / 2
    #[test]
    fn test_couette_flow() {
        let lub = LubricationModel::new(1e-3, 1e-4, 0.01, 0.01);
        let q = lub.couette_flow_rate(1.0);
        let expected = 1.0 * 1e-4 / 2.0;
        assert!(
            (q - expected).abs() < TOL,
            "Couette flow mismatch: {:.6e} vs {:.6e}",
            q,
            expected
        );
    }

    // 18. Poiseuille flow opposes positive pressure gradient
    #[test]
    fn test_poiseuille_flow_sign() {
        let lub = LubricationModel::new(1e-3, 1e-4, 0.01, 0.01);
        let q = lub.poiseuille_flow_rate(1e6);
        assert!(q < 0.0, "Poiseuille flow opposes positive dP/dx");
    }

    // 19. Squeeze film force opposes approach
    #[test]
    fn test_squeeze_film() {
        let lub = LubricationModel::new(1e-3, 1e-4, 0.01, 0.01);
        let f = lub.squeeze_film_force(0.001);
        assert!(f > 0.0, "Squeeze film force should resist approach");
    }

    // 20. Slider bearing: zero load for parallel surfaces
    #[test]
    fn test_slider_parallel_zero_load() {
        let lub = LubricationModel::new(1e-3, 1e-4, 0.01, 0.01);
        let load = lub.slider_bearing_load(1e-4, 1e-4, 1.0);
        assert!(
            load.abs() < TOL,
            "Parallel surfaces should give zero load, got {:.6e}",
            load
        );
    }

    // 21. Rivulet width depends on contact angle
    #[test]
    fn test_rivulet_width() {
        let r1 = RivuletFlow::new(RHO_WATER, MU_WATER, GAMMA_WATER, 0.5, 0.3, 1e-6);
        let r2 = RivuletFlow::new(RHO_WATER, MU_WATER, GAMMA_WATER, 0.5, 0.6, 1e-6);
        assert!(
            r2.max_width() > r1.max_width(),
            "Larger contact angle → wider rivulet"
        );
    }

    // 22. Dip coating: Landau-Levich film thickness positive
    #[test]
    fn test_dip_coating_thickness() {
        let dc = DipCoating::new(RHO_WATER, MU_WATER, GAMMA_WATER, 0.01);
        let h = dc.film_thickness();
        assert!(h > 0.0, "Film thickness must be positive");
        // For water at 0.01 m/s, expect ~10-100 μm
        assert!(
            h < 1e-3,
            "Film thickness should be sub-mm for dip coating, got {:.6e}",
            h
        );
    }

    // 23. Spin coating: thickness decreases with time
    #[test]
    fn test_spin_coating_thinning() {
        let sc = SpinCoating::new(RHO_WATER, MU_WATER, 200.0 * 2.0 * PI / 60.0, 1e-3);
        let h1 = sc.thickness(0.0);
        let h2 = sc.thickness(1.0);
        let h3 = sc.thickness(10.0);
        assert!((h1 - 1e-3).abs() < TOL, "t=0 thickness should be h0");
        assert!(h2 < h1, "Film should thin with time");
        assert!(h3 < h2, "Film should continue thinning");
    }

    // 24. Spin coating: thinning rate negative
    #[test]
    fn test_spin_coating_rate() {
        let sc = SpinCoating::new(RHO_WATER, MU_WATER, 200.0 * 2.0 * PI / 60.0, 1e-3);
        let rate = sc.thinning_rate(1.0);
        assert!(rate < 0.0, "Thinning rate should be negative");
    }

    // 25. Film simulation: nodes maintain positive thickness
    #[test]
    fn test_film_simulation_positive_thickness() {
        let nodes = vec![
            FilmNode::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
            FilmNode::new([0.01, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
            FilmNode::new([0.02, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
        ];
        let edges = vec![(0, 1), (1, 2)];
        let props = FluidFilmProperties::water(0.5);
        let mut sim = FluidFilmSimulation::new(nodes, edges, props, [0.0, -G_ACCEL, 0.0]);
        for _ in 0..10 {
            sim.step(1e-4);
        }
        for (i, n) in sim.nodes.iter().enumerate() {
            assert!(
                n.thickness > 0.0,
                "Node {} thickness must be positive, got {:.6e}",
                i,
                n.thickness
            );
        }
    }

    // 26. Film simulation: total volume approximately conserved (no source)
    #[test]
    fn test_film_simulation_volume_conservation() {
        let nodes = vec![
            FilmNode::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2e-4),
            FilmNode::new([0.01, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
            FilmNode::new([0.02, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
        ];
        let edges = vec![(0, 1), (1, 2)];
        let props = FluidFilmProperties::water(0.5);
        let mut sim = FluidFilmSimulation::new(nodes, edges, props, [0.0, 0.0, 0.0]);
        let v0 = sim.total_volume();
        for _ in 0..100 {
            sim.step(1e-5);
        }
        let v1 = sim.total_volume();
        assert!(
            (v1 - v0).abs() / v0 < 0.01,
            "Volume should be conserved: {:.6e} -> {:.6e}",
            v0,
            v1
        );
    }

    // 27. Uniformity CV = 0 for uniform film
    #[test]
    fn test_uniformity_cv_uniform() {
        let nodes = vec![
            FilmNode::new([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
            FilmNode::new([0.01, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
            FilmNode::new([0.02, 0.0, 0.0], [0.0, 1.0, 0.0], 1e-4),
        ];
        let edges = vec![(0, 1), (1, 2)];
        let props = FluidFilmProperties::water(0.5);
        let sim = FluidFilmSimulation::new(nodes, edges, props, [0.0, 0.0, 0.0]);
        let cv = sim.uniformity_cv();
        assert!(
            cv.abs() < TOL,
            "CV should be 0 for uniform film, got {:.6e}",
            cv
        );
    }

    // 28. Bond number dimensionless and positive
    #[test]
    fn test_bond_number() {
        let fp = FluidFilmProperties::water(0.0);
        let bo = fp.bond_number(0.001);
        assert!(bo > 0.0, "Bond number must be positive");
    }

    // 29. Marangoni average velocity proportional to film thickness
    #[test]
    fn test_marangoni_velocity_thickness() {
        let mg1 = MarangoniFlow::new(-1.5e-4, 1e-3, 1e-4);
        let mg2 = MarangoniFlow::new(-1.5e-4, 1e-3, 2e-4);
        let v1 = mg1.average_velocity(100.0);
        let v2 = mg2.average_velocity(100.0);
        assert!(
            (v2 / v1 - 2.0).abs() < 1e-6,
            "Velocity should scale linearly with h: ratio = {:.6}",
            v2 / v1
        );
    }

    // 30. Thin film stability check
    #[test]
    fn test_thin_film_stability() {
        let tfr_thick = ThinFilmRupture::new(1e-20, 0.072, 1e-3, 1e-3);
        let tfr_thin = ThinFilmRupture::new(1e-20, 0.072, 1e-3, 1e-12);
        assert!(tfr_thick.is_stable(), "Thick film should be stable");
        assert!(!tfr_thin.is_stable(), "Very thin film should be unstable");
    }
}
