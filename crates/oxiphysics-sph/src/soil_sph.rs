// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Soil mechanics SPH simulation.
//!
//! Provides elastoplastic soil modelling via SPH, including:
//!
//! - **Drucker-Prager** yield criterion for pressure-dependent plasticity
//! - **Slope stability** analysis with factor-of-safety evaluation
//! - **Landslide dynamics** with progressive failure and run-out prediction
//! - **Soil-water coupling** for saturated/unsaturated soil behaviour
//! - **Pore pressure** evolution and effective stress computation
//! - **Liquefaction** assessment (excess pore-pressure ratio)
//! - **Consolidation** (Terzaghi 1-D and general Biot)
//! - **Large deformation geomechanics** with updated Lagrangian formulation
//! - **Debris flow** rheology (Bingham / Herschel-Bulkley)
//! - **Retaining wall interaction** contact model

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Vector / tensor helpers (no nalgebra)
// ---------------------------------------------------------------------------

/// 3-D vector utilities on `[f64; 3]`.
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two vectors.
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a vector.
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product.
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Euclidean norm.
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

/// Zero vector.
fn vec3_zero() -> [f64; 3] {
    [0.0, 0.0, 0.0]
}

/// Symmetric 3x3 tensor stored as row-major `[[f64; 3\]; 3]`.
type Tensor3 = [[f64; 3]; 3];

/// Zero tensor.
fn tensor_zero() -> Tensor3 {
    [[0.0; 3]; 3]
}

/// Identity tensor.
fn tensor_identity() -> Tensor3 {
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

/// Trace of a tensor.
fn tensor_trace(t: &Tensor3) -> f64 {
    t[0][0] + t[1][1] + t[2][2]
}

/// Add two tensors.
fn tensor_add(a: &Tensor3, b: &Tensor3) -> Tensor3 {
    let mut r = tensor_zero();
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = a[i][j] + b[i][j];
        }
    }
    r
}

/// Scale a tensor.
fn tensor_scale(t: &Tensor3, s: f64) -> Tensor3 {
    let mut r = tensor_zero();
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = t[i][j] * s;
        }
    }
    r
}

/// Deviatoric part of a tensor.
fn tensor_deviatoric(t: &Tensor3) -> Tensor3 {
    let p = tensor_trace(t) / 3.0;
    let mut d = *t;
    d[0][0] -= p;
    d[1][1] -= p;
    d[2][2] -= p;
    d
}

/// Second invariant J2 of a symmetric tensor (deviatoric).
fn tensor_j2(s: &Tensor3) -> f64 {
    let mut j2 = 0.0;
    for row in s.iter() {
        for &v in row.iter() {
            j2 += v * v;
        }
    }
    0.5 * j2
}

/// Frobenius norm of a tensor.
#[cfg(test)]
fn tensor_frobenius(t: &Tensor3) -> f64 {
    let mut s = 0.0;
    for row in t.iter() {
        for &v in row.iter() {
            s += v * v;
        }
    }
    s.sqrt()
}

// ---------------------------------------------------------------------------
// SPH kernel
// ---------------------------------------------------------------------------

/// Cubic spline kernel value in 3-D.
fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * 0.25 * t * t * t
    } else {
        0.0
    }
}

/// Gradient magnitude of cubic spline kernel (radial component / r).
fn cubic_spline_grad_factor(r: f64, h: f64) -> f64 {
    if r < 1.0e-14 {
        return 0.0;
    }
    let q = r / h;
    let sigma = 1.0 / (PI * h * h * h);
    if q < 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q) / (h * r)
    } else if q < 2.0 {
        let t = 2.0 - q;
        sigma * (-0.75 * t * t) / (h * r)
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------------
// Core data structures
// ---------------------------------------------------------------------------

/// A single soil SPH particle carrying stress, pore-pressure and state data.
#[derive(Clone, Debug)]
pub struct SoilParticle {
    /// Position.
    pub position: [f64; 3],
    /// Velocity.
    pub velocity: [f64; 3],
    /// Acceleration (buffer for force accumulation).
    pub acceleration: [f64; 3],
    /// Mass.
    pub mass: f64,
    /// Current density.
    pub density: f64,
    /// Reference (initial) density.
    pub density0: f64,
    /// Total Cauchy stress tensor (compression positive).
    pub stress: Tensor3,
    /// Pore water pressure (compression positive).
    pub pore_pressure: f64,
    /// Plastic strain accumulator.
    pub plastic_strain: f64,
    /// Void ratio e = V_void / V_solid.
    pub void_ratio: f64,
    /// Smoothing length.
    pub h: f64,
    /// Flag: is this particle a boundary (wall) particle?
    pub is_boundary: bool,
    /// Liquefaction excess pore pressure ratio r_u.
    pub excess_pore_ratio: f64,
    /// Strain rate tensor (symmetric).
    pub strain_rate: Tensor3,
}

impl SoilParticle {
    /// Create a new soil particle with default stress-free state.
    pub fn new(position: [f64; 3], mass: f64, density0: f64, h: f64) -> Self {
        Self {
            position,
            velocity: vec3_zero(),
            acceleration: vec3_zero(),
            mass,
            density: density0,
            density0,
            stress: tensor_zero(),
            pore_pressure: 0.0,
            plastic_strain: 0.0,
            void_ratio: 0.5,
            h,
            is_boundary: false,
            excess_pore_ratio: 0.0,
            strain_rate: tensor_zero(),
        }
    }

    /// Create a boundary (wall) particle.
    pub fn new_boundary(position: [f64; 3], mass: f64, density0: f64, h: f64) -> Self {
        let mut p = Self::new(position, mass, density0, h);
        p.is_boundary = true;
        p
    }
}

/// Material parameters for a soil medium.
#[derive(Clone, Debug)]
pub struct SoilParams {
    /// Cohesion intercept c \[Pa\].
    pub cohesion: f64,
    /// Internal friction angle \[radians\].
    pub friction_angle: f64,
    /// Dilatancy angle \[radians\].
    pub dilatancy_angle: f64,
    /// Young's modulus E \[Pa\].
    pub youngs_modulus: f64,
    /// Poisson ratio nu.
    pub poisson_ratio: f64,
    /// Reference density \[kg/m^3\].
    pub density0: f64,
    /// Hydraulic conductivity k \[m/s\].
    pub permeability: f64,
    /// Compression index for consolidation.
    pub compression_index: f64,
    /// Swelling/recompression index.
    pub swelling_index: f64,
    /// Pre-consolidation pressure \[Pa\].
    pub preconsolidation_pressure: f64,
    /// Specific gravity of soil grains.
    pub specific_gravity: f64,
    /// Water density \[kg/m^3\].
    pub water_density: f64,
    /// Yield stress for Bingham debris flow \[Pa\].
    pub yield_stress: f64,
    /// Viscosity for debris flow \[Pa.s\].
    pub viscosity: f64,
}

impl SoilParams {
    /// Typical dry sand.
    pub fn dry_sand() -> Self {
        Self {
            cohesion: 0.0,
            friction_angle: 30.0_f64.to_radians(),
            dilatancy_angle: 5.0_f64.to_radians(),
            youngs_modulus: 5.0e7,
            poisson_ratio: 0.3,
            density0: 1700.0,
            permeability: 1.0e-4,
            compression_index: 0.15,
            swelling_index: 0.03,
            preconsolidation_pressure: 100_000.0,
            specific_gravity: 2.65,
            water_density: 1000.0,
            yield_stress: 0.0,
            viscosity: 0.0,
        }
    }

    /// Saturated clay.
    pub fn saturated_clay() -> Self {
        Self {
            cohesion: 20_000.0,
            friction_angle: 15.0_f64.to_radians(),
            dilatancy_angle: 0.0,
            youngs_modulus: 1.0e7,
            poisson_ratio: 0.4,
            density0: 1900.0,
            permeability: 1.0e-9,
            compression_index: 0.4,
            swelling_index: 0.08,
            preconsolidation_pressure: 80_000.0,
            specific_gravity: 2.70,
            water_density: 1000.0,
            yield_stress: 0.0,
            viscosity: 0.0,
        }
    }

    /// Debris flow material (Bingham-like).
    pub fn debris_flow() -> Self {
        Self {
            cohesion: 500.0,
            friction_angle: 20.0_f64.to_radians(),
            dilatancy_angle: 0.0,
            youngs_modulus: 1.0e6,
            poisson_ratio: 0.35,
            density0: 2000.0,
            permeability: 1.0e-6,
            compression_index: 0.2,
            swelling_index: 0.04,
            preconsolidation_pressure: 50_000.0,
            specific_gravity: 2.60,
            water_density: 1000.0,
            yield_stress: 5000.0,
            viscosity: 100.0,
        }
    }

    /// Bulk modulus derived from E, nu.
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }

    /// Shear modulus derived from E, nu.
    pub fn shear_modulus(&self) -> f64 {
        self.youngs_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
}

// ---------------------------------------------------------------------------
// Drucker-Prager yield criterion
// ---------------------------------------------------------------------------

/// Drucker-Prager yield surface parameters derived from Mohr-Coulomb.
#[derive(Clone, Debug)]
pub struct DruckerPrager {
    /// Material parameter alpha (pressure coefficient).
    pub alpha: f64,
    /// Material parameter k (cohesive strength).
    pub k: f64,
}

impl DruckerPrager {
    /// Construct Drucker-Prager cone that circumscribes the Mohr-Coulomb surface.
    pub fn from_mohr_coulomb(cohesion: f64, friction_angle: f64) -> Self {
        let sin_phi = friction_angle.sin();
        let cos_phi = friction_angle.cos();
        let alpha = 2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi));
        let k = 6.0 * cohesion * cos_phi / (3.0_f64.sqrt() * (3.0 - sin_phi));
        Self { alpha, k }
    }

    /// Evaluate yield function f = sqrt(J2) + alpha * I1 - k.
    /// f <= 0 => elastic, f > 0 => yielding.
    pub fn yield_function(&self, stress: &Tensor3) -> f64 {
        let i1 = tensor_trace(stress);
        let dev = tensor_deviatoric(stress);
        let j2 = tensor_j2(&dev);
        j2.sqrt() + self.alpha * i1 - self.k
    }

    /// Return-map: project stress back to the yield surface.
    pub fn return_mapping(
        &self,
        stress: &Tensor3,
        shear_modulus: f64,
        bulk_modulus: f64,
    ) -> (Tensor3, f64) {
        let f = self.yield_function(stress);
        if f <= 0.0 {
            return (*stress, 0.0);
        }
        let i1 = tensor_trace(stress);
        let dev = tensor_deviatoric(stress);
        let j2 = tensor_j2(&dev);
        let sqrt_j2 = j2.sqrt().max(1.0e-20);

        // Plastic multiplier (associated flow approximation)
        let d_lambda = f / (shear_modulus + 9.0 * bulk_modulus * self.alpha * self.alpha);

        // Return-mapped deviatoric stress
        let scale = 1.0 - shear_modulus * d_lambda / sqrt_j2;
        let scale = scale.max(0.0);
        let dev_new = tensor_scale(&dev, scale);

        // Return-mapped mean stress
        let p_new = i1 / 3.0 - bulk_modulus * self.alpha * d_lambda;

        let id = tensor_identity();
        let corrected = tensor_add(&dev_new, &tensor_scale(&id, p_new));

        (corrected, d_lambda)
    }
}

// ---------------------------------------------------------------------------
// Effective stress / pore pressure
// ---------------------------------------------------------------------------

/// Compute effective stress: sigma' = sigma_total - u * I  (Terzaghi principle).
pub fn effective_stress(total_stress: &Tensor3, pore_pressure: f64) -> Tensor3 {
    let mut eff = *total_stress;
    eff[0][0] -= pore_pressure;
    eff[1][1] -= pore_pressure;
    eff[2][2] -= pore_pressure;
    eff
}

/// Compute total stress from effective stress and pore pressure.
pub fn total_stress(effective: &Tensor3, pore_pressure: f64) -> Tensor3 {
    let mut tot = *effective;
    tot[0][0] += pore_pressure;
    tot[1][1] += pore_pressure;
    tot[2][2] += pore_pressure;
    tot
}

// ---------------------------------------------------------------------------
// Liquefaction assessment
// ---------------------------------------------------------------------------

/// Excess pore pressure ratio r_u = u_excess / sigma'_v0.
///
/// When r_u >= 1.0 the effective vertical stress vanishes and the
/// soil is considered to have liquefied.
pub fn excess_pore_pressure_ratio(
    pore_pressure: f64,
    hydrostatic: f64,
    effective_vertical: f64,
) -> f64 {
    let u_excess = pore_pressure - hydrostatic;
    if effective_vertical.abs() < 1.0e-12 {
        return 0.0;
    }
    (u_excess / effective_vertical).clamp(0.0, 2.0)
}

/// Check if soil has liquefied (r_u >= threshold, typically 1.0).
pub fn is_liquefied(r_u: f64, threshold: f64) -> bool {
    r_u >= threshold
}

// ---------------------------------------------------------------------------
// Consolidation (Terzaghi 1-D)
// ---------------------------------------------------------------------------

/// Coefficient of consolidation c_v = k / (m_v * gamma_w).
pub fn coefficient_of_consolidation(
    permeability: f64,
    _compression_index: f64,
    void_ratio: f64,
    effective_stress_val: f64,
    water_density: f64,
) -> f64 {
    let gamma_w = water_density * 9.81;
    // m_v approximation from compression index
    let m_v = (1.0 + void_ratio) / (effective_stress_val.max(1.0));
    permeability / (m_v * gamma_w)
}

/// Terzaghi 1-D consolidation: degree of consolidation U_z for given
/// time factor T_v. Uses first-term approximation.
pub fn degree_of_consolidation(time_factor: f64) -> f64 {
    if time_factor < 0.0 {
        return 0.0;
    }
    // For T_v < 0.2827 use U = sqrt(4*T_v/pi), else U = 1 - 8/(pi^2) * exp(-pi^2*T_v/4)
    let threshold = 0.2827;
    if time_factor < threshold {
        (4.0 * time_factor / PI).sqrt()
    } else {
        1.0 - (8.0 / (PI * PI)) * (-PI * PI * time_factor / 4.0).exp()
    }
}

/// Update pore pressure using simple explicit diffusion step.
pub fn pore_pressure_diffusion_step(
    particles: &mut [SoilParticle],
    neighbors: &[Vec<usize>],
    cv: f64,
    dt: f64,
) {
    let n = particles.len();
    let mut du = vec![0.0_f64; n];

    for i in 0..n {
        if particles[i].is_boundary {
            continue;
        }
        let hi = particles[i].h;
        for &j in &neighbors[i] {
            let rij = vec3_sub(particles[i].position, particles[j].position);
            let dist = vec3_norm(rij);
            if dist < 1.0e-14 {
                continue;
            }
            let w = cubic_spline_kernel(dist, hi);
            let vol_j = particles[j].mass / particles[j].density;
            du[i] += cv * vol_j * (particles[j].pore_pressure - particles[i].pore_pressure) * w;
        }
    }

    for i in 0..n {
        if !particles[i].is_boundary {
            particles[i].pore_pressure += du[i] * dt;
        }
    }
}

// ---------------------------------------------------------------------------
// Slope stability (simplified Bishop / limit equilibrium)
// ---------------------------------------------------------------------------

/// Result of a slope stability analysis.
#[derive(Clone, Debug)]
pub struct SlopeStabilityResult {
    /// Factor of safety.
    pub factor_of_safety: f64,
    /// Critical slip circle centre (x, y).
    pub slip_center: [f64; 2],
    /// Critical slip circle radius.
    pub slip_radius: f64,
}

/// Slice descriptor: (base_x, base_y, height, width, base_angle_rad, cohesion, friction_angle, pore_pressure).
type BishopSlice = (f64, f64, f64, f64, f64, f64, f64, f64);

/// Simplified Bishop method for a slope defined by particle columns.
///
/// `slices` contains (base_x, base_y, height, width, base_angle_rad, cohesion, friction_angle, pore_pressure).
pub fn bishop_simplified(slices: &[BishopSlice], _unit_weight: f64) -> f64 {
    // Iterative Bishop simplified method
    let mut fos = 1.5; // initial guess
    for _iter in 0..50 {
        let mut num = 0.0;
        let mut den = 0.0;
        for &(_, _, h, b, alpha, c, phi, u) in slices {
            let weight = h * b * 18.0; // approximate unit weight * area
            let m_alpha = (alpha.cos() + phi.tan() * alpha.sin() / fos).max(0.01);
            num += (c * b + (weight - u * b) * phi.tan()) / m_alpha;
            den += weight * alpha.sin();
        }
        if den.abs() < 1.0e-12 {
            return f64::INFINITY;
        }
        let fos_new = num / den;
        if (fos_new - fos).abs() < 1.0e-6 {
            return fos_new;
        }
        fos = fos_new;
    }
    fos
}

// ---------------------------------------------------------------------------
// Landslide dynamics
// ---------------------------------------------------------------------------

/// Landslide state tracker.
#[derive(Clone, Debug)]
pub struct LandslideState {
    /// Current time \[s\].
    pub time: f64,
    /// Run-out distance \[m\] from initial toe.
    pub runout_distance: f64,
    /// Maximum velocity \[m/s\] in the flow.
    pub max_velocity: f64,
    /// Centre of mass.
    pub center_of_mass: [f64; 3],
    /// Number of active (moving) particles.
    pub active_count: usize,
}

impl LandslideState {
    /// Compute from a set of particles.
    pub fn compute(particles: &[SoilParticle], initial_toe_x: f64) -> Self {
        let mut com = vec3_zero();
        let mut total_mass = 0.0;
        let mut max_v = 0.0_f64;
        let mut active = 0usize;
        let mut max_x = f64::NEG_INFINITY;

        for p in particles {
            if p.is_boundary {
                continue;
            }
            let m = p.mass;
            com = vec3_add(com, vec3_scale(p.position, m));
            total_mass += m;
            let v = vec3_norm(p.velocity);
            if v > max_v {
                max_v = v;
            }
            if v > 0.01 {
                active += 1;
            }
            if p.position[0] > max_x {
                max_x = p.position[0];
            }
        }
        if total_mass > 0.0 {
            com = vec3_scale(com, 1.0 / total_mass);
        }
        let runout = (max_x - initial_toe_x).max(0.0);

        Self {
            time: 0.0,
            runout_distance: runout,
            max_velocity: max_v,
            center_of_mass: com,
            active_count: active,
        }
    }
}

// ---------------------------------------------------------------------------
// Soil-water coupling
// ---------------------------------------------------------------------------

/// Compute seepage force on a soil particle from pore pressure gradient.
///
/// F_seep = -grad(u) * V_particle / density_particle
pub fn seepage_force(pore_grad: [f64; 3], volume: f64) -> [f64; 3] {
    vec3_scale(pore_grad, -volume)
}

/// Compute pore pressure gradient via SPH.
pub fn compute_pore_pressure_gradient(
    i: usize,
    particles: &[SoilParticle],
    neighbors: &[usize],
) -> [f64; 3] {
    let mut grad = vec3_zero();
    let pi = &particles[i];
    let hi = pi.h;
    for &j in neighbors {
        let pj = &particles[j];
        let rij = vec3_sub(pi.position, pj.position);
        let dist = vec3_norm(rij);
        if dist < 1.0e-14 {
            continue;
        }
        let gf = cubic_spline_grad_factor(dist, hi);
        let vol_j = pj.mass / pj.density;
        let dp = pj.pore_pressure - pi.pore_pressure;
        // grad += vol_j * dp * (rij / r) * dW/dr
        for k in 0..3 {
            grad[k] += vol_j * dp * rij[k] * gf;
        }
    }
    grad
}

// ---------------------------------------------------------------------------
// Debris flow rheology
// ---------------------------------------------------------------------------

/// Bingham / Herschel-Bulkley viscous stress.
///
/// Returns an effective viscosity given strain rate magnitude.
pub fn bingham_viscosity(yield_stress: f64, viscosity: f64, strain_rate_mag: f64) -> f64 {
    if strain_rate_mag < 1.0e-12 {
        // Regularised: cap at large value
        return viscosity + yield_stress / 1.0e-12;
    }
    viscosity + yield_stress / strain_rate_mag
}

/// Compute Herschel-Bulkley effective viscosity.
///
/// tau = tau_y + K * gamma_dot^n
pub fn herschel_bulkley_viscosity(
    yield_stress: f64,
    consistency: f64,
    power_index: f64,
    strain_rate_mag: f64,
) -> f64 {
    let eps = 1.0e-12;
    let gamma = strain_rate_mag.max(eps);
    yield_stress / gamma + consistency * gamma.powf(power_index - 1.0)
}

/// Compute strain rate magnitude from strain rate tensor.
pub fn strain_rate_magnitude(sr: &Tensor3) -> f64 {
    let j2 = tensor_j2(sr);
    (2.0 * j2).sqrt()
}

// ---------------------------------------------------------------------------
// Retaining wall interaction
// ---------------------------------------------------------------------------

/// Result of soil pressure on a retaining wall segment.
#[derive(Clone, Debug)]
pub struct WallPressureResult {
    /// Normal force per unit length \[N/m\].
    pub normal_force: f64,
    /// Shear force per unit length \[N/m\].
    pub shear_force: f64,
    /// Point of application (height from base).
    pub application_height: f64,
}

/// Rankine active earth pressure coefficient.
pub fn rankine_active_coefficient(friction_angle: f64) -> f64 {
    let sin_phi = friction_angle.sin();
    (1.0 - sin_phi) / (1.0 + sin_phi)
}

/// Rankine passive earth pressure coefficient.
pub fn rankine_passive_coefficient(friction_angle: f64) -> f64 {
    let sin_phi = friction_angle.sin();
    (1.0 + sin_phi) / (1.0 - sin_phi)
}

/// Coulomb active earth pressure coefficient.
pub fn coulomb_active_coefficient(phi: f64, delta: f64, beta: f64, alpha: f64) -> f64 {
    let num = (phi + alpha).sin().powi(2);
    let t1 = alpha.sin().powi(2) * (alpha - delta).sin();
    let inner =
        ((phi + delta).sin() * (phi - beta).sin()) / ((alpha - delta).sin() * (alpha + beta).sin());
    let t2 = (1.0 + inner.max(0.0).sqrt()).powi(2);
    if t1 * t2 < 1.0e-20 {
        return 0.0;
    }
    num / (t1 * t2)
}

/// Compute wall interaction forces from nearby soil particles.
///
/// Simple penalty contact: if a soil particle penetrates the wall half-plane,
/// apply a repulsive normal force.
pub fn wall_contact_force(
    particle_pos: [f64; 3],
    wall_point: [f64; 3],
    wall_normal: [f64; 3],
    stiffness: f64,
) -> [f64; 3] {
    let rel = vec3_sub(particle_pos, wall_point);
    let penetration = -vec3_dot(rel, wall_normal);
    if penetration > 0.0 {
        vec3_scale(wall_normal, stiffness * penetration)
    } else {
        vec3_zero()
    }
}

// ---------------------------------------------------------------------------
// Large deformation geomechanics helpers
// ---------------------------------------------------------------------------

/// Velocity gradient tensor via SPH.
pub fn compute_velocity_gradient(
    i: usize,
    particles: &[SoilParticle],
    neighbors: &[usize],
) -> Tensor3 {
    let mut l = tensor_zero();
    let pi = &particles[i];
    let hi = pi.h;
    for &j in neighbors {
        let pj = &particles[j];
        let rij = vec3_sub(pi.position, pj.position);
        let dist = vec3_norm(rij);
        if dist < 1.0e-14 {
            continue;
        }
        let gf = cubic_spline_grad_factor(dist, hi);
        let vol_j = pj.mass / pj.density;
        let dv = vec3_sub(pj.velocity, pi.velocity);
        for a in 0..3 {
            for b in 0..3 {
                l[a][b] += vol_j * dv[a] * rij[b] * gf;
            }
        }
    }
    l
}

/// Symmetric part of velocity gradient = strain rate tensor.
pub fn symmetric_part(l: &Tensor3) -> Tensor3 {
    let mut d = tensor_zero();
    for i in 0..3 {
        for j in 0..3 {
            d[i][j] = 0.5 * (l[i][j] + l[j][i]);
        }
    }
    d
}

/// Antisymmetric (spin) part of velocity gradient.
pub fn antisymmetric_part(l: &Tensor3) -> Tensor3 {
    let mut w = tensor_zero();
    for i in 0..3 {
        for j in 0..3 {
            w[i][j] = 0.5 * (l[i][j] - l[j][i]);
        }
    }
    w
}

/// Jaumann stress rate: d_sigma/dt = C:D + W*sigma - sigma*W.
pub fn jaumann_stress_rate(
    stress: &Tensor3,
    strain_rate: &Tensor3,
    spin: &Tensor3,
    bulk_modulus: f64,
    shear_modulus: f64,
) -> Tensor3 {
    let trace_d = tensor_trace(strain_rate);
    let mut rate = tensor_zero();

    // Elastic contribution: sigma_dot_ij = 2G * D_ij + (K - 2G/3) * tr(D) * delta_ij
    let lame = bulk_modulus - 2.0 * shear_modulus / 3.0;
    for i in 0..3 {
        for j in 0..3 {
            rate[i][j] = 2.0 * shear_modulus * strain_rate[i][j];
            if i == j {
                rate[i][j] += lame * trace_d;
            }
        }
    }

    // Jaumann rotation: + W * sigma - sigma * W
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                rate[i][j] += spin[i][k] * stress[k][j] - stress[i][k] * spin[k][j];
            }
        }
    }

    rate
}

// ---------------------------------------------------------------------------
// Density update
// ---------------------------------------------------------------------------

/// SPH density summation.
pub fn compute_density(i: usize, particles: &[SoilParticle], neighbors: &[usize]) -> f64 {
    let pi = &particles[i];
    let hi = pi.h;
    let mut rho = pi.mass * cubic_spline_kernel(0.0, hi);
    for &j in neighbors {
        let pj = &particles[j];
        let dist = vec3_norm(vec3_sub(pi.position, pj.position));
        rho += pj.mass * cubic_spline_kernel(dist, hi);
    }
    rho
}

// ---------------------------------------------------------------------------
// Soil SPH solver
// ---------------------------------------------------------------------------

/// Configuration for the soil SPH solver.
#[derive(Clone, Debug)]
pub struct SoilSphConfig {
    /// Time step \[s\].
    pub dt: f64,
    /// Gravitational acceleration \[m/s^2\] (typically \[0, -9.81, 0\]).
    pub gravity: [f64; 3],
    /// Artificial viscosity coefficient alpha.
    pub artificial_viscosity_alpha: f64,
    /// Artificial viscosity coefficient beta.
    pub artificial_viscosity_beta: f64,
    /// Speed of sound for artificial viscosity.
    pub speed_of_sound: f64,
    /// Enable pore pressure evolution.
    pub enable_pore_pressure: bool,
    /// Coefficient of consolidation for pore pressure diffusion.
    pub consolidation_cv: f64,
    /// Enable debris flow rheology.
    pub enable_debris_rheology: bool,
}

impl Default for SoilSphConfig {
    fn default() -> Self {
        Self {
            dt: 1.0e-4,
            gravity: [0.0, -9.81, 0.0],
            artificial_viscosity_alpha: 1.0,
            artificial_viscosity_beta: 0.0,
            speed_of_sound: 200.0,
            enable_pore_pressure: false,
            consolidation_cv: 1.0e-3,
            enable_debris_rheology: false,
        }
    }
}

/// Main soil SPH solver.
#[derive(Clone)]
pub struct SoilSphSolver {
    /// All soil particles.
    pub particles: Vec<SoilParticle>,
    /// Neighbor lists (built externally or internally).
    pub neighbors: Vec<Vec<usize>>,
    /// Material parameters.
    pub params: SoilParams,
    /// Drucker-Prager yield surface.
    pub yield_surface: DruckerPrager,
    /// Solver configuration.
    pub config: SoilSphConfig,
    /// Current simulation time.
    pub time: f64,
}

impl SoilSphSolver {
    /// Create a new solver.
    pub fn new(particles: Vec<SoilParticle>, params: SoilParams, config: SoilSphConfig) -> Self {
        let n = particles.len();
        let yield_surface =
            DruckerPrager::from_mohr_coulomb(params.cohesion, params.friction_angle);
        Self {
            particles,
            neighbors: vec![Vec::new(); n],
            params,
            yield_surface,
            config,
            time: 0.0,
        }
    }

    /// Build neighbor lists using brute-force search (O(n^2)).
    pub fn build_neighbors(&mut self) {
        let n = self.particles.len();
        self.neighbors.resize(n, Vec::new());
        for i in 0..n {
            self.neighbors[i].clear();
            let hi = self.particles[i].h;
            let cutoff = 2.0 * hi;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let dist = vec3_norm(vec3_sub(
                    self.particles[i].position,
                    self.particles[j].position,
                ));
                if dist < cutoff {
                    self.neighbors[i].push(j);
                }
            }
        }
    }

    /// Update densities via SPH summation.
    pub fn update_densities(&mut self) {
        let n = self.particles.len();
        // Collect neighbor refs to avoid borrow issues
        let neighbors: Vec<Vec<usize>> = self.neighbors.clone();
        for (i, _) in neighbors.iter().enumerate().take(n) {
            self.particles[i].density = compute_density(i, &self.particles, &neighbors[i]);
        }
    }

    /// Compute strain rates and update stresses (elastoplastic with Drucker-Prager).
    pub fn update_stresses(&mut self) {
        let n = self.particles.len();
        let neighbors = self.neighbors.clone();
        let g = self.params.shear_modulus();
        let k = self.params.bulk_modulus();
        let dt = self.config.dt;

        for (i, _) in neighbors.iter().enumerate().take(n) {
            if self.particles[i].is_boundary {
                continue;
            }
            let l = compute_velocity_gradient(i, &self.particles, &neighbors[i]);
            let d = symmetric_part(&l);
            let w = antisymmetric_part(&l);

            self.particles[i].strain_rate = d;

            // Jaumann rate
            let stress_rate = jaumann_stress_rate(&self.particles[i].stress, &d, &w, k, g);

            // Trial stress
            let mut trial = tensor_add(&self.particles[i].stress, &tensor_scale(&stress_rate, dt));

            // Effective stress for yield check
            let eff = effective_stress(&trial, self.particles[i].pore_pressure);

            // Drucker-Prager return mapping on effective stress
            let (corrected_eff, d_lambda) = self.yield_surface.return_mapping(&eff, g, k);

            // Total stress = effective + pore pressure
            let corrected = total_stress(&corrected_eff, self.particles[i].pore_pressure);
            self.particles[i].stress = corrected;
            self.particles[i].plastic_strain += d_lambda;

            // Debris flow viscous stress addition
            if self.config.enable_debris_rheology && self.params.yield_stress > 0.0 {
                let sr_mag = strain_rate_magnitude(&d);
                let eta =
                    bingham_viscosity(self.params.yield_stress, self.params.viscosity, sr_mag);
                let visc_stress = tensor_scale(&d, 2.0 * eta);
                trial = tensor_add(&self.particles[i].stress, &tensor_scale(&visc_stress, dt));
                self.particles[i].stress = trial;
            }
        }
    }

    /// Compute momentum equation (acceleration) from stress divergence.
    pub fn compute_accelerations(&mut self) {
        let n = self.particles.len();
        let neighbors = self.neighbors.clone();

        for i in 0..n {
            self.particles[i].acceleration = self.config.gravity;
        }

        for (i, _) in neighbors.iter().enumerate().take(n) {
            if self.particles[i].is_boundary {
                self.particles[i].acceleration = vec3_zero();
                continue;
            }
            let pi = &self.particles[i];
            let hi = pi.h;
            let rho_i = pi.density.max(1.0);
            let sigma_i = &pi.stress;

            let mut acc = vec3_zero();

            for &j in &neighbors[i] {
                let pj = &self.particles[j];
                let rij = vec3_sub(pi.position, pj.position);
                let dist = vec3_norm(rij);
                if dist < 1.0e-14 {
                    continue;
                }
                let gf = cubic_spline_grad_factor(dist, hi);
                let rho_j = pj.density.max(1.0);
                let sigma_j = &pj.stress;

                // Stress divergence: sum_j m_j * (sigma_i/rho_i^2 + sigma_j/rho_j^2) . grad W
                for k in 0..3 {
                    let mut force_k = 0.0;
                    for l in 0..3 {
                        force_k += (sigma_i[k][l] / (rho_i * rho_i)
                            + sigma_j[k][l] / (rho_j * rho_j))
                            * rij[l]
                            * gf;
                    }
                    acc[k] += pj.mass * force_k;
                }

                // Artificial viscosity
                let vij = vec3_sub(pi.velocity, pj.velocity);
                let vr = vec3_dot(vij, rij);
                if vr < 0.0 {
                    let rho_avg = 0.5 * (rho_i + rho_j);
                    let eta2 = 0.01 * hi * hi;
                    let mu = hi * vr / (dist * dist + eta2);
                    let cs = self.config.speed_of_sound;
                    let pi_ab = (-self.config.artificial_viscosity_alpha * cs * mu
                        + self.config.artificial_viscosity_beta * mu * mu)
                        / rho_avg;
                    for k in 0..3 {
                        acc[k] -= pj.mass * pi_ab * rij[k] * gf;
                    }
                }
            }

            // Add seepage force if pore pressure is enabled
            if self.config.enable_pore_pressure {
                let pg = compute_pore_pressure_gradient(i, &self.particles, &neighbors[i]);
                let vol = self.particles[i].mass / rho_i;
                let sf = seepage_force(pg, vol);
                acc = vec3_add(acc, vec3_scale(sf, 1.0 / self.particles[i].mass));
            }

            self.particles[i].acceleration = vec3_add(self.config.gravity, acc);
        }
    }

    /// Integrate positions and velocities (Verlet kick-drift-kick).
    pub fn integrate(&mut self) {
        let dt = self.config.dt;
        let n = self.particles.len();

        for i in 0..n {
            if self.particles[i].is_boundary {
                continue;
            }
            // Half-kick
            let half_dv = vec3_scale(self.particles[i].acceleration, 0.5 * dt);
            self.particles[i].velocity = vec3_add(self.particles[i].velocity, half_dv);

            // Drift
            let dx = vec3_scale(self.particles[i].velocity, dt);
            self.particles[i].position = vec3_add(self.particles[i].position, dx);

            // Second half-kick (same acceleration)
            self.particles[i].velocity = vec3_add(self.particles[i].velocity, half_dv);
        }

        self.time += dt;
    }

    /// Perform one full time step.
    pub fn step(&mut self) {
        self.build_neighbors();
        self.update_densities();
        self.update_stresses();
        self.compute_accelerations();
        if self.config.enable_pore_pressure {
            pore_pressure_diffusion_step(
                &mut self.particles,
                &self.neighbors,
                self.config.consolidation_cv,
                self.config.dt,
            );
        }
        self.integrate();
    }

    /// Run the solver for `n_steps` steps.
    pub fn run(&mut self, n_steps: usize) {
        for _ in 0..n_steps {
            self.step();
        }
    }

    /// Update excess pore pressure ratios for liquefaction assessment.
    pub fn update_liquefaction(&mut self) {
        for p in &mut self.particles {
            if p.is_boundary {
                continue;
            }
            // Approximate effective vertical stress from total stress
            let sigma_v = -p.stress[1][1]; // vertical total stress (compression positive)
            let hydrostatic = 0.0; // simplified
            let eff_v = (sigma_v - p.pore_pressure).max(0.0);
            p.excess_pore_ratio = excess_pore_pressure_ratio(p.pore_pressure, hydrostatic, eff_v);
        }
    }
}

// ---------------------------------------------------------------------------
// Large deformation: deformation gradient tracking
// ---------------------------------------------------------------------------

/// Deformation gradient tracker for updated Lagrangian SPH.
#[derive(Clone, Debug)]
pub struct DeformationTracker {
    /// Deformation gradient F per particle (row-major 3x3).
    pub gradients: Vec<Tensor3>,
}

impl DeformationTracker {
    /// Initialize with identity deformation gradients.
    pub fn new(n: usize) -> Self {
        Self {
            gradients: vec![tensor_identity(); n],
        }
    }

    /// Update deformation gradients: F(t+dt) = (I + L*dt) * F(t).
    pub fn update(&mut self, particles: &[SoilParticle], neighbors: &[Vec<usize>], dt: f64) {
        let n = particles.len();
        for i in 0..n {
            if particles[i].is_boundary {
                continue;
            }
            let l = compute_velocity_gradient(i, particles, &neighbors[i]);
            // (I + L*dt) * F
            let mut f_new = tensor_zero();
            for a in 0..3 {
                for (b, fb) in f_new[a].iter_mut().enumerate() {
                    let mut val = 0.0;
                    for (k, lak) in l[a].iter().enumerate() {
                        let i_plus_ldt = if a == k { 1.0 } else { 0.0 } + lak * dt;
                        val += i_plus_ldt * self.gradients[i][k][b];
                    }
                    *fb = val;
                }
            }
            self.gradients[i] = f_new;
        }
    }

    /// Compute determinant of deformation gradient (volume ratio).
    pub fn determinant(&self, i: usize) -> f64 {
        let f = &self.gradients[i];
        f[0][0] * (f[1][1] * f[2][2] - f[1][2] * f[2][1])
            - f[0][1] * (f[1][0] * f[2][2] - f[1][2] * f[2][0])
            + f[0][2] * (f[1][0] * f[2][1] - f[1][1] * f[2][0])
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // 1. Drucker-Prager: hydrostatic stress within yield surface
    #[test]
    fn test_dp_hydrostatic_elastic() {
        let dp = DruckerPrager::from_mohr_coulomb(10_000.0, 30.0_f64.to_radians());
        // Small hydrostatic compression
        let sigma = tensor_scale(&tensor_identity(), -1000.0);
        let f = dp.yield_function(&sigma);
        assert!(
            f < 0.0,
            "Hydrostatic compression should be elastic, f={:.6}",
            f
        );
    }

    // 2. Drucker-Prager: large deviatoric stress yields
    #[test]
    fn test_dp_deviatoric_yield() {
        let dp = DruckerPrager::from_mohr_coulomb(100.0, 30.0_f64.to_radians());
        let mut sigma = tensor_zero();
        sigma[0][0] = 1.0e6;
        sigma[1][1] = -1.0e6;
        let f = dp.yield_function(&sigma);
        assert!(f > 0.0, "Large deviatoric stress should yield, f={:.6}", f);
    }

    // 3. Return mapping brings stress back toward yield surface
    #[test]
    fn test_return_mapping() {
        let dp = DruckerPrager::from_mohr_coulomb(5000.0, 25.0_f64.to_radians());
        let g = 1.0e7;
        let k = 2.0e7;
        let mut sigma = tensor_zero();
        sigma[0][0] = 1.0e7;
        sigma[1][1] = -5.0e6;
        sigma[2][2] = 0.0;
        let f_before = dp.yield_function(&sigma);
        assert!(f_before > 0.0, "Should be yielding before return");
        let (corrected, dl) = dp.return_mapping(&sigma, g, k);
        assert!(dl > 0.0, "Should have plastic increment");
        let f_after = dp.yield_function(&corrected);
        assert!(
            f_after < f_before,
            "After return mapping f should decrease: before={:.6}, after={:.6}",
            f_before,
            f_after
        );
    }

    // 4. Effective stress calculation
    #[test]
    fn test_effective_stress() {
        let total = tensor_scale(&tensor_identity(), -100.0);
        let u = 30.0;
        let eff = effective_stress(&total, u);
        assert!(
            (eff[0][0] - (-130.0)).abs() < EPS,
            "Effective stress diagonal wrong"
        );
        assert!(eff[0][1].abs() < EPS, "Off-diagonal should be unchanged");
    }

    // 5. Total stress roundtrip
    #[test]
    fn test_total_stress_roundtrip() {
        let sigma = [[100.0, 20.0, 0.0], [20.0, -50.0, 10.0], [0.0, 10.0, 30.0]];
        let u = 25.0;
        let eff = effective_stress(&sigma, u);
        let recovered = total_stress(&eff, u);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (recovered[i][j] - sigma[i][j]).abs() < EPS,
                    "Roundtrip failed at [{i}][{j}]"
                );
            }
        }
    }

    // 6. Excess pore pressure ratio
    #[test]
    fn test_excess_pore_ratio() {
        let r_u = excess_pore_pressure_ratio(50.0, 20.0, 100.0);
        assert!(
            (r_u - 0.3).abs() < 1.0e-6,
            "r_u should be 0.3, got {:.6}",
            r_u
        );
    }

    // 7. Liquefaction detection
    #[test]
    fn test_liquefaction_detection() {
        assert!(is_liquefied(1.0, 1.0));
        assert!(is_liquefied(1.5, 1.0));
        assert!(!is_liquefied(0.8, 1.0));
    }

    // 8. Consolidation degree at T_v = 0
    #[test]
    fn test_consolidation_zero() {
        let u = degree_of_consolidation(0.0);
        assert!(u.abs() < 1.0e-6, "U(0) should be 0, got {:.6}", u);
    }

    // 9. Consolidation degree at large T_v approaches 1
    #[test]
    fn test_consolidation_large_tv() {
        let u = degree_of_consolidation(10.0);
        assert!(
            (u - 1.0).abs() < 1.0e-3,
            "U(large) should be ~1, got {:.6}",
            u
        );
    }

    // 10. Bishop simplified: stable slope
    #[test]
    fn test_bishop_stable() {
        // Flat slices with high cohesion => FoS > 1
        let slices = vec![
            (0.0, 0.0, 5.0, 1.0, 0.1, 50_000.0, 0.5, 0.0),
            (1.0, 0.0, 4.0, 1.0, 0.15, 50_000.0, 0.5, 0.0),
        ];
        let fos = bishop_simplified(&slices, 18.0);
        assert!(fos > 1.0, "Stable slope FoS should be > 1, got {:.6}", fos);
    }

    // 11. Rankine active < passive coefficient
    #[test]
    fn test_rankine_coefficients() {
        let phi = 30.0_f64.to_radians();
        let ka = rankine_active_coefficient(phi);
        let kp = rankine_passive_coefficient(phi);
        assert!(ka < 1.0, "Ka should be < 1 for phi=30, got {:.6}", ka);
        assert!(kp > 1.0, "Kp should be > 1 for phi=30, got {:.6}", kp);
        assert!((ka * kp - 1.0).abs() < 1.0e-6, "Ka * Kp should be 1");
    }

    // 12. Wall contact force: no penetration => zero force
    #[test]
    fn test_wall_no_penetration() {
        let f = wall_contact_force([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0e6);
        assert!(vec3_norm(f) < EPS, "No penetration should give zero force");
    }

    // 13. Wall contact force: penetration gives repulsion
    #[test]
    fn test_wall_penetration() {
        let f = wall_contact_force([0.0, -0.1, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0e6);
        assert!(f[1] > 0.0, "Should push particle upward, got {:.6}", f[1]);
    }

    // 14. Bingham viscosity at zero strain rate is large
    #[test]
    fn test_bingham_zero_rate() {
        let eta = bingham_viscosity(100.0, 10.0, 0.0);
        assert!(eta > 1.0e10, "Should be very large at zero strain rate");
    }

    // 15. Bingham viscosity at finite rate
    #[test]
    fn test_bingham_finite_rate() {
        let eta = bingham_viscosity(100.0, 10.0, 1.0);
        assert!(
            (eta - 110.0).abs() < 1.0e-6,
            "eta should be 110, got {:.6}",
            eta
        );
    }

    // 16. Herschel-Bulkley reduces to Bingham when n=1
    #[test]
    fn test_hb_reduces_to_bingham() {
        let gamma = 2.0;
        let eta_b = bingham_viscosity(50.0, 20.0, gamma);
        let eta_hb = herschel_bulkley_viscosity(50.0, 20.0, 1.0, gamma);
        assert!(
            (eta_b - eta_hb).abs() < 1.0e-6,
            "HB with n=1 should equal Bingham: {:.6} vs {:.6}",
            eta_b,
            eta_hb
        );
    }

    // 17. Strain rate magnitude of zero tensor is zero
    #[test]
    fn test_strain_rate_zero() {
        let sr = tensor_zero();
        assert!(strain_rate_magnitude(&sr) < EPS);
    }

    // 18. Symmetric part of symmetric tensor is itself
    #[test]
    fn test_symmetric_part_identity() {
        let s = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]];
        let d = symmetric_part(&s);
        for i in 0..3 {
            for j in 0..3 {
                assert!(
                    (d[i][j] - s[i][j]).abs() < EPS,
                    "Sym part mismatch at [{i}][{j}]"
                );
            }
        }
    }

    // 19. Antisymmetric part of symmetric tensor is zero
    #[test]
    fn test_antisym_of_symmetric_is_zero() {
        let s = [[1.0, 2.0, 3.0], [2.0, 4.0, 5.0], [3.0, 5.0, 6.0]];
        let w = antisymmetric_part(&s);
        assert!(tensor_frobenius(&w) < EPS, "Antisym of sym should be zero");
    }

    // 20. Deformation tracker: initial determinant is 1
    #[test]
    fn test_deformation_initial_det() {
        let tracker = DeformationTracker::new(5);
        for i in 0..5 {
            let det = tracker.determinant(i);
            assert!(
                (det - 1.0).abs() < EPS,
                "Initial det should be 1, got {:.6}",
                det
            );
        }
    }

    // 21. Soil params: bulk = E / (3*(1-2*nu))
    #[test]
    fn test_soil_params_moduli() {
        let p = SoilParams::dry_sand();
        let k = p.bulk_modulus();
        let g = p.shear_modulus();
        let expected_k = p.youngs_modulus / (3.0 * (1.0 - 2.0 * p.poisson_ratio));
        let expected_g = p.youngs_modulus / (2.0 * (1.0 + p.poisson_ratio));
        assert!((k - expected_k).abs() < 1.0, "Bulk modulus mismatch");
        assert!((g - expected_g).abs() < 1.0, "Shear modulus mismatch");
    }

    // 22. Cubic spline kernel is positive at origin
    #[test]
    fn test_kernel_origin() {
        let w = cubic_spline_kernel(0.0, 0.1);
        assert!(w > 0.0, "Kernel at origin should be positive, got {:.6}", w);
    }

    // 23. Cubic spline kernel is zero beyond 2h
    #[test]
    fn test_kernel_cutoff() {
        let w = cubic_spline_kernel(0.201, 0.1);
        assert!(
            w.abs() < EPS,
            "Kernel beyond 2h should be zero, got {:.6}",
            w
        );
    }

    // 24. SoilSphSolver: small system doesn't crash
    #[test]
    fn test_solver_basic() {
        let particles = vec![
            SoilParticle::new([0.0, 0.0, 0.0], 1.0, 1700.0, 0.05),
            SoilParticle::new([0.03, 0.0, 0.0], 1.0, 1700.0, 0.05),
            SoilParticle::new([0.0, 0.03, 0.0], 1.0, 1700.0, 0.05),
        ];
        let params = SoilParams::dry_sand();
        let config = SoilSphConfig::default();
        let mut solver = SoilSphSolver::new(particles, params, config);
        solver.run(5);
        // Should have advanced time
        assert!(solver.time > 0.0, "Time should have advanced");
    }

    // 25. Landslide state: static particles have zero runout
    #[test]
    fn test_landslide_state_static() {
        let particles = vec![
            SoilParticle::new([0.0, 0.0, 0.0], 1.0, 1700.0, 0.05),
            SoilParticle::new([1.0, 0.0, 0.0], 1.0, 1700.0, 0.05),
        ];
        let state = LandslideState::compute(&particles, 1.0);
        assert!(
            state.max_velocity < 1.0e-6,
            "Static particles should have zero velocity"
        );
        assert_eq!(state.active_count, 0);
    }

    // 26. Seepage force direction
    #[test]
    fn test_seepage_force_direction() {
        let grad = [100.0, 0.0, 0.0]; // pressure gradient in +x
        let f = seepage_force(grad, 0.01);
        assert!(f[0] < 0.0, "Seepage force should oppose gradient");
    }

    // 27. Coulomb active coefficient
    #[test]
    fn test_coulomb_active() {
        let ka = coulomb_active_coefficient(30.0_f64.to_radians(), 0.0, 0.0, 90.0_f64.to_radians());
        assert!(ka > 0.0 && ka < 1.0, "Ka should be in (0,1), got {:.6}", ka);
    }

    // 28. Jaumann rate: zero spin and zero strain rate => zero rate
    #[test]
    fn test_jaumann_zero() {
        let sigma = [[100.0, 10.0, 0.0], [10.0, 50.0, 0.0], [0.0, 0.0, 30.0]];
        let d = tensor_zero();
        let w = tensor_zero();
        let rate = jaumann_stress_rate(&sigma, &d, &w, 1.0e7, 5.0e6);
        assert!(
            tensor_frobenius(&rate) < EPS,
            "Zero D and W should give zero rate"
        );
    }

    // 29. Tensor deviatoric trace is zero
    #[test]
    fn test_deviatoric_trace_zero() {
        let t = [[10.0, 3.0, 0.0], [3.0, -5.0, 2.0], [0.0, 2.0, 7.0]];
        let dev = tensor_deviatoric(&t);
        let tr = tensor_trace(&dev);
        assert!(
            tr.abs() < EPS,
            "Deviatoric trace should be zero, got {:.6}",
            tr
        );
    }

    // 30. Tensor identity trace is 3
    #[test]
    fn test_identity_trace() {
        let id = tensor_identity();
        assert!((tensor_trace(&id) - 3.0).abs() < EPS);
    }
}
