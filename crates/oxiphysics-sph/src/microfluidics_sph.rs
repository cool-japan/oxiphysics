#![allow(clippy::ptr_arg, clippy::type_complexity)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Microfluidic SPH module.
//!
//! Covers low-Reynolds-number, capillary-dominated flows in microfluidic devices:
//!
//! - Droplet formation in T-junction and flow-focusing geometries
//! - Electro-osmotic flow (EOF) via the Debye-Hückel slip boundary
//! - Dielectrophoresis (DEP) particle trapping
//! - Digital microfluidics: droplet actuation by electrowetting-on-dielectric (EWOD)
//! - Lab-on-chip passive mixing (chaotic advection, staggered herringbone)
//! - Surface-tension-dominated flow; Laplace pressure at curved interfaces
//! - Wettability: Young contact angle embedded in SPH boundary conditions
//! - Droplet coalescence and breakup governed by critical capillary number
//! - Color-gradient SPH for immiscible fluid-interface tension
//! - Marangoni flow driven by surface-tension gradients

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use rand::RngExt;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Boltzmann constant (J K⁻¹).
const K_B: f64 = 1.380_649e-23;

/// Permittivity of free space (F m⁻¹).
const EPS_0: f64 = 8.854_187_817e-12;

/// Elementary charge (C).
const E_CHARGE: f64 = 1.602_176_634e-19;

/// Avogadro's number (mol⁻¹).
const N_AV: f64 = 6.022_140_76e23;

/// Viscosity of water at 25 °C (Pa s).
const ETA_WATER: f64 = 8.9e-4;

/// Surface tension of water–air at 25 °C (N m⁻¹).
const SIGMA_WATER_AIR: f64 = 0.0720;

/// Dielectric constant of water (dimensionless).
const EPS_WATER: f64 = 78.5;

/// Reference temperature (K).
const T_REF: f64 = 298.15;

// ---------------------------------------------------------------------------
// Vector helpers (plain [f64;3])
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

/// Subtract two 3-vectors (a − b).
#[inline]
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scale a 3-vector by scalar s.
#[inline]
pub fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
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

/// Normalise a 3-vector; returns zero vector if norm < eps.
#[inline]
pub fn normalise3(a: [f64; 3]) -> [f64; 3] {
    let n = norm3(a);
    if n < 1.0e-30 {
        [0.0; 3]
    } else {
        scale3(a, 1.0 / n)
    }
}

// ---------------------------------------------------------------------------
// SPH kernel functions (Wendland C2, smoothing length h)
// ---------------------------------------------------------------------------

/// Evaluate the Wendland C2 kernel W(r, h) in 3-D.
///
/// Returns the kernel value for inter-particle distance `r` and smoothing
/// length `h`.
pub fn wendland_c2(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 {
        return 0.0;
    }
    let alpha = 21.0 / (16.0 * PI * h.powi(3));
    let t = 1.0 - 0.5 * q;
    alpha * t.powi(4) * (2.0 * q + 1.0)
}

/// Gradient magnitude of the Wendland C2 kernel with respect to `r`.
///
/// Returns dW/dr (scalar; multiply by rhat to get the vector gradient).
pub fn wendland_c2_grad(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q >= 2.0 || r < 1.0e-30 {
        return 0.0;
    }
    let alpha = 21.0 / (16.0 * PI * h.powi(3));
    let t = 1.0 - 0.5 * q;
    // d/dr [ t^4 (2q+1) ] = (1/h) d/dq [ t^4 (2q+1) ]
    // = (1/h) [ 4 t^3 (-0.5)(2q+1) + t^4 * 2 ]
    // = (1/h) t^3 [ -2(2q+1) + 2t ] * wait, simpler:
    // d/dq [t^4 (2q+1)] = -2 t^3 (2q+1) + 2 t^4 = 2 t^3 (t - (2q+1))
    // t - (2q+1) = 1 - 0.5q - 2q - 1 = -2.5q
    alpha * (1.0 / h) * 2.0 * t.powi(3) * (-2.5 * q)
}

// ---------------------------------------------------------------------------
// Microfluidic particle
// ---------------------------------------------------------------------------

/// A single SPH particle representing a fluid parcel in a microfluidic device.
///
/// The `phase` field distinguishes two immiscible phases (0 = fluid A, 1 = fluid B).
#[derive(Debug, Clone)]
pub struct MicroParticle {
    /// Position (m).
    pub pos: [f64; 3],
    /// Velocity (m s⁻¹).
    pub vel: [f64; 3],
    /// Acceleration / force per unit mass (m s⁻²).
    pub acc: [f64; 3],
    /// Pressure (Pa).
    pub pressure: f64,
    /// Density (kg m⁻³).
    pub density: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Phase index (0 or 1).
    pub phase: u8,
    /// Color gradient for color-gradient SPH (vector, 3-D).
    pub color_grad: [f64; 3],
    /// Interface curvature (m⁻¹).
    pub curvature: f64,
    /// Electric potential (V).
    pub phi_e: f64,
    /// Net charge density at particle site (C m⁻³).
    pub rho_e: f64,
}

impl MicroParticle {
    /// Create a new `MicroParticle` with default zero fields.
    pub fn new(pos: [f64; 3], mass: f64, h: f64, phase: u8) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            acc: [0.0; 3],
            pressure: 0.0,
            density: 1000.0,
            mass,
            h,
            phase,
            color_grad: [0.0; 3],
            curvature: 0.0,
            phi_e: 0.0,
            rho_e: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Equation of state (weakly compressible, Tait)
// ---------------------------------------------------------------------------

/// Compute weakly-compressible Tait pressure for a particle.
///
/// p = B \[(ρ/ρ₀)^γ − 1\], where B = ρ₀ c_s² / γ.
pub fn tait_pressure(rho: f64, rho0: f64, c_s: f64, gamma: f64) -> f64 {
    let b = rho0 * c_s * c_s / gamma;
    b * ((rho / rho0).powf(gamma) - 1.0)
}

// ---------------------------------------------------------------------------
// Color-gradient SPH (Lafaurie / Brackbill CSF)
// ---------------------------------------------------------------------------

/// Assign color values (0.0 for phase-0, 1.0 for phase-1) and accumulate
/// the color gradient at each particle using SPH summation.
///
/// `particles` is mutated; the `color_grad` field is updated.
pub fn compute_color_gradient(particles: &mut Vec<MicroParticle>) {
    let n = particles.len();
    // Snapshot positions and phases to avoid borrow conflict
    let snap: Vec<([f64; 3], f64, f64, u8)> = particles
        .iter()
        .map(|p| (p.pos, p.mass, p.density, p.phase))
        .collect();

    for i in 0..n {
        let mut grad = [0.0f64; 3];
        let ci = if snap[i].3 == 0 { 0.0 } else { 1.0 };
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_ij = sub3(snap[i].0, snap[j].0);
            let r = norm3(r_ij);
            let h_avg = (particles[i].h + particles[j].h) * 0.5;
            let dw = wendland_c2_grad(r, h_avg);
            let rhat = if r > 1.0e-30 {
                scale3(r_ij, 1.0 / r)
            } else {
                [0.0; 3]
            };
            let cj = if snap[j].3 == 0 { 0.0 } else { 1.0 };
            let vol_j = snap[j].1 / snap[j].2; // m_j / rho_j
            let coeff = (cj - ci) * vol_j * dw;
            grad = add3(grad, scale3(rhat, coeff));
        }
        particles[i].color_grad = grad;
    }
}

/// Compute interface curvature from the divergence of the normalised color
/// gradient: κ = −∇ · n̂, where n̂ = ∇C / |∇C|.
///
/// `particles` is mutated; `curvature` field is updated.
pub fn compute_curvature(particles: &mut Vec<MicroParticle>) {
    let n = particles.len();
    let snap: Vec<([f64; 3], [f64; 3], f64, f64)> = particles
        .iter()
        .map(|p| (p.pos, p.color_grad, p.mass, p.density))
        .collect();

    for i in 0..n {
        let ni_hat = normalise3(snap[i].1);
        let mut div_n = 0.0f64;
        for j in 0..n {
            if i == j {
                continue;
            }
            let r_ij = sub3(snap[i].0, snap[j].0);
            let r = norm3(r_ij);
            let h_avg = (particles[i].h + particles[j].h) * 0.5;
            let dw = wendland_c2_grad(r, h_avg);
            let rhat = if r > 1.0e-30 {
                scale3(r_ij, 1.0 / r)
            } else {
                [0.0; 3]
            };
            let nj_hat = normalise3(snap[j].1);
            let dn = sub3(nj_hat, ni_hat);
            let vol_j = snap[j].2 / snap[j].3;
            div_n += dot3(dn, rhat) * dw * vol_j;
        }
        particles[i].curvature = -div_n;
    }
}

/// Apply the Continuum Surface Force (CSF) body force to interface particles.
///
/// Adds σ κ ∇C / ρ to the acceleration of each particle, where σ is the
/// surface tension coefficient (N m⁻¹).
pub fn apply_csf_force(particles: &mut Vec<MicroParticle>, sigma: f64) {
    for p in particles.iter_mut() {
        let grad_mag = norm3(p.color_grad);
        if grad_mag < 1.0e-20 {
            continue;
        }
        let coeff = sigma * p.curvature / p.density;
        let f_csf = scale3(p.color_grad, coeff);
        p.acc = add3(p.acc, f_csf);
    }
}

// ---------------------------------------------------------------------------
// Laplace pressure
// ---------------------------------------------------------------------------

/// Laplace pressure across a spherical interface: ΔP = 2σ / R.
///
/// Returns the pressure jump (Pa) for surface tension `sigma` (N m⁻¹) and
/// droplet radius `r` (m).
pub fn laplace_pressure_sphere(sigma: f64, r: f64) -> f64 {
    2.0 * sigma / r
}

/// Laplace pressure across a cylindrical interface: ΔP = σ / R.
///
/// Returns the pressure jump (Pa).
pub fn laplace_pressure_cylinder(sigma: f64, r: f64) -> f64 {
    sigma / r
}

// ---------------------------------------------------------------------------
// Wettability: Young contact angle
// ---------------------------------------------------------------------------

/// Young–Dupré relation: cos θ = (σ_SG − σ_SL) / σ_LG.
///
/// Returns the equilibrium contact angle θ (radians) given the three
/// interfacial energies (N m⁻¹).
pub fn young_contact_angle(sigma_sg: f64, sigma_sl: f64, sigma_lg: f64) -> f64 {
    let cos_theta = (sigma_sg - sigma_sl) / sigma_lg;
    cos_theta.clamp(-1.0, 1.0).acos()
}

/// Apply a wettability boundary correction to the color gradient of particles
/// near a solid wall.
///
/// The solid wall is at y = 0 (normal in +y direction). Particles within
/// `wall_dist` are corrected so that the interface satisfies `contact_angle`
/// (radians).
pub fn apply_wettability_bc(
    particles: &mut Vec<MicroParticle>,
    contact_angle: f64,
    wall_dist: f64,
) {
    for p in particles.iter_mut() {
        if p.pos[1] > wall_dist {
            continue;
        }
        // Wall normal is [0, 1, 0]; tangential component is in xz-plane
        let cg_mag = norm3(p.color_grad);
        if cg_mag < 1.0e-20 {
            continue;
        }
        // Desired surface normal direction
        let wall_normal = [0.0f64, 1.0, 0.0];
        let tan = normalise3(sub3(
            p.color_grad,
            scale3(wall_normal, dot3(p.color_grad, wall_normal)),
        ));
        let n_wall = add3(
            scale3(wall_normal, contact_angle.cos()),
            scale3(tan, contact_angle.sin()),
        );
        p.color_grad = scale3(n_wall, cg_mag);
    }
}

// ---------------------------------------------------------------------------
// T-junction droplet formation
// ---------------------------------------------------------------------------

/// Geometry parameters for a T-junction microfluidic device.
#[derive(Debug, Clone)]
pub struct TJunctionGeometry {
    /// Width of the main channel (m).
    pub main_width: f64,
    /// Width of the side channel (m).
    pub side_width: f64,
    /// Height of the channels (m).
    pub height: f64,
    /// Junction x-coordinate (m).
    pub junction_x: f64,
    /// Junction y-coordinate (m).
    pub junction_y: f64,
}

impl TJunctionGeometry {
    /// Create a default T-junction with typical microfluidic dimensions (100 µm scale).
    pub fn default_100um() -> Self {
        Self {
            main_width: 100.0e-6,
            side_width: 50.0e-6,
            height: 50.0e-6,
            junction_x: 200.0e-6,
            junction_y: 50.0e-6,
        }
    }
}

/// Estimate the droplet length formed in a T-junction using the Garstecki
/// scaling law:
///
/// L/w = 1 + α (Q_d / Q_c)
///
/// where w is the main-channel width, Q_d and Q_c are the dispersed- and
/// continuous-phase flow rates, and α ≈ 1.
pub fn t_junction_droplet_length(
    main_width: f64,
    flow_dispersed: f64,
    flow_continuous: f64,
    alpha: f64,
) -> f64 {
    main_width * (1.0 + alpha * flow_dispersed / flow_continuous)
}

/// Capillary number Ca = η V / σ.
///
/// Dimensionless ratio of viscous to capillary forces.
pub fn capillary_number(viscosity: f64, velocity: f64, sigma: f64) -> f64 {
    viscosity * velocity / sigma
}

/// Critical capillary number for droplet breakup.
///
/// Uses the simplified Taylor criterion Ca_crit ≈ 0.4 for a drop in a
/// uniform shear flow (viscosity ratio λ = 1). Returns `true` if the drop
/// will break.
pub fn droplet_breaks(ca: f64) -> bool {
    ca > 0.4
}

// ---------------------------------------------------------------------------
// Electro-osmotic flow (EOF)
// ---------------------------------------------------------------------------

/// Debye length λ_D = sqrt(ε ε₀ k_B T / (2 n₀ z² e²)).
///
/// `n0` is the bulk ion number density (m⁻³), `z` the ion valence,
/// `eps_r` the relative permittivity of the solvent, `temp` temperature (K).
pub fn debye_length(n0: f64, z: f64, eps_r: f64, temp: f64) -> f64 {
    let numerator = eps_r * EPS_0 * K_B * temp;
    let denominator = 2.0 * n0 * (z * E_CHARGE).powi(2);
    (numerator / denominator).sqrt()
}

/// Helmholtz-Smoluchowski slip velocity for electro-osmotic flow.
///
/// v_EOF = −(ε ε₀ ζ E) / η
///
/// `zeta` is the zeta potential (V), `E` the electric field magnitude (V m⁻¹),
/// `eps_r` the relative permittivity, `eta` the dynamic viscosity (Pa s).
pub fn eof_slip_velocity(zeta: f64, e_field: f64, eps_r: f64, eta: f64) -> f64 {
    -(eps_r * EPS_0 * zeta * e_field) / eta
}

/// Apply an electro-osmotic body force to all fluid particles in the channel.
///
/// Adds f_eof = (ρ_e / ρ) × E_field to particle accelerations.
/// `e_field_vec` is the applied electric field (V m⁻¹).
pub fn apply_eof_body_force(particles: &mut Vec<MicroParticle>, e_field_vec: [f64; 3]) {
    for p in particles.iter_mut() {
        if p.rho_e.abs() < 1.0e-30 {
            continue;
        }
        let f = scale3(e_field_vec, p.rho_e / p.density);
        p.acc = add3(p.acc, f);
    }
}

/// Assign Boltzmann-distributed charge density to particles.
///
/// For a symmetric z:z electrolyte: ρ_e = −2 n₀ z e sinh(z e φ / k_B T).
/// `phi` is the electric potential at the particle site (V).
pub fn boltzmann_charge_density(n0: f64, z: f64, phi: f64, temp: f64) -> f64 {
    let beta = E_CHARGE / (K_B * temp);
    -2.0 * n0 * z * E_CHARGE * (z * beta * phi).sinh()
}

// ---------------------------------------------------------------------------
// Dielectrophoresis (DEP)
// ---------------------------------------------------------------------------

/// Clausius-Mossotti factor (real part) for a spherical particle.
///
/// Re\[K\] = (ε_p − ε_m) / (ε_p + 2 ε_m)
///
/// `eps_p` particle permittivity, `eps_m` medium permittivity (both relative).
pub fn clausius_mossotti_real(eps_p: f64, eps_m: f64) -> f64 {
    (eps_p - eps_m) / (eps_p + 2.0 * eps_m)
}

/// DEP force on a spherical particle in a non-uniform electric field.
///
/// F_DEP = 2π ε_m ε₀ r³ Re\[K\] ∇|E|²
///
/// `r_particle` particle radius (m), `grad_e2` gradient of |E|² (V² m⁻³),
/// Returns force vector (N).
pub fn dep_force(r_particle: f64, eps_m: f64, k_cm: f64, grad_e2: [f64; 3]) -> [f64; 3] {
    let coeff = 2.0 * PI * eps_m * EPS_0 * r_particle.powi(3) * k_cm;
    scale3(grad_e2, coeff)
}

/// Determine if a DEP particle is trapped in a local field minimum/maximum.
///
/// Positive DEP (pDEP): particle moves to high-field regions (K > 0).
/// Negative DEP (nDEP): particle moves to low-field regions (K < 0).
/// `threshold_force` is the minimum force required to overcome Brownian motion.
pub fn dep_is_trapped(dep_force_mag: f64, threshold_force: f64, k_cm: f64) -> bool {
    dep_force_mag > threshold_force && k_cm.abs() > 0.01
}

// ---------------------------------------------------------------------------
// Digital microfluidics: EWOD actuation
// ---------------------------------------------------------------------------

/// Electrowetting-on-dielectric (EWOD) contact angle modification.
///
/// Lippmann-Young equation: cos θ_V = cos θ₀ + (C V²) / (2 σ_LG)
///
/// `theta0` is the intrinsic contact angle (rad), `capacitance` is the
/// dielectric capacitance per unit area (F m⁻²), `voltage` applied voltage (V),
/// `sigma_lg` liquid-gas surface tension (N m⁻¹).
pub fn ewod_contact_angle(theta0: f64, capacitance: f64, voltage: f64, sigma_lg: f64) -> f64 {
    let cos_v = theta0.cos() + (capacitance * voltage * voltage) / (2.0 * sigma_lg);
    cos_v.clamp(-1.0, 1.0).acos()
}

/// EWOD driving pressure for droplet actuation.
///
/// ΔP_EWOD = σ_LG (cos θ_on − cos θ_off) / d
///
/// where d is the gap between the electrodes (m).
pub fn ewod_driving_pressure(sigma_lg: f64, theta_on: f64, theta_off: f64, gap: f64) -> f64 {
    sigma_lg * (theta_on.cos() - theta_off.cos()) / gap
}

/// EWOD droplet velocity estimate (Darcy flow through gap).
///
/// v ≈ ΔP d² / (12 η L)
///
/// `dp` driving pressure (Pa), `gap` electrode gap (m), `eta` viscosity (Pa s),
/// `length` droplet length along motion direction (m).
pub fn ewod_droplet_velocity(dp: f64, gap: f64, eta: f64, length: f64) -> f64 {
    dp * gap * gap / (12.0 * eta * length)
}

// ---------------------------------------------------------------------------
// Lab-on-chip mixing
// ---------------------------------------------------------------------------

/// Péclet number Pe = V L / D.
///
/// Ratio of advective to diffusive transport.
pub fn peclet_number(velocity: f64, length: f64, diffusivity: f64) -> f64 {
    velocity * length / diffusivity
}

/// Mixing length in a staggered herringbone micromixer.
///
/// Approximate exponential mixing model: C(x) = exp(−x / L_mix)
/// L_mix ≈ d_h / (2 ln(2) × n_grooves × mixing_efficiency)
///
/// `d_h` hydraulic diameter (m), `n_grooves` number of grooves per cycle,
/// `eff` mixing efficiency coefficient.
pub fn herringbone_mixing_length(d_h: f64, n_grooves: f64, eff: f64) -> f64 {
    d_h / (2.0 * 2_f64.ln() * n_grooves * eff)
}

/// Estimate mixing time for chaotic advection: t_mix = L_mix / V.
pub fn mixing_time(l_mix: f64, velocity: f64) -> f64 {
    l_mix / velocity
}

/// Apply a body force to simulate chaotic advection mixing.
///
/// Adds a sinusoidal transverse perturbation force to each particle based on
/// position, mimicking a herringbone groove geometry.
///
/// `amplitude` (m s⁻²), `wavenumber` (m⁻¹), `phase_offset` (rad).
pub fn apply_herringbone_mixing_force(
    particles: &mut Vec<MicroParticle>,
    amplitude: f64,
    wavenumber: f64,
    phase_offset: f64,
) {
    for p in particles.iter_mut() {
        let f_y = amplitude * (wavenumber * p.pos[0] + phase_offset).sin();
        p.acc[1] += f_y;
    }
}

// ---------------------------------------------------------------------------
// Droplet coalescence and breakup
// ---------------------------------------------------------------------------

/// Check if two approaching droplets will coalesce based on capillary number.
///
/// Coalescence occurs when Ca < Ca_crit_coalesce (film drainage completes).
/// Returns `true` if coalescence is expected.
pub fn droplets_coalesce(
    ca: f64,
    ca_crit_coalesce: f64,
    surface_fraction: f64,
    surface_fraction_threshold: f64,
) -> bool {
    ca < ca_crit_coalesce && surface_fraction < surface_fraction_threshold
}

/// Compute the Ohnesorge number: Oh = η / sqrt(ρ σ R).
///
/// Ratio of viscous to inertio-capillary forces; relevant for breakup.
pub fn ohnesorge_number(eta: f64, rho: f64, sigma: f64, r: f64) -> f64 {
    eta / (rho * sigma * r).sqrt()
}

/// Weber number: We = ρ V² R / σ.
pub fn weber_number(rho: f64, vel: f64, r: f64, sigma: f64) -> f64 {
    rho * vel * vel * r / sigma
}

/// Breakup criterion: a droplet breaks when We / Oh > threshold (~10–20).
pub fn droplet_breakup_criterion(we: f64, oh: f64, threshold: f64) -> bool {
    (we / oh) > threshold
}

/// Split a droplet particle cloud into two child clouds along the x-axis.
///
/// Returns two `Vec<MicroParticle>` with the original mass distributed evenly
/// and the phase unchanged. Children are displaced by ±`offset` in x.
pub fn split_droplet_cloud(
    particles: &[MicroParticle],
    offset: f64,
) -> (Vec<MicroParticle>, Vec<MicroParticle>) {
    let mut left = Vec::new();
    let mut right = Vec::new();
    for (i, p) in particles.iter().enumerate() {
        let mut child = p.clone();
        if i % 2 == 0 {
            child.pos[0] -= offset;
            left.push(child);
        } else {
            child.pos[0] += offset;
            right.push(child);
        }
    }
    (left, right)
}

// ---------------------------------------------------------------------------
// Marangoni flow
// ---------------------------------------------------------------------------

/// Marangoni body force: f = (dσ/dT) ∇T / ρ (or (dσ/dC) ∇C for solutal).
///
/// `dsigma_dt` is dσ/dT (N m⁻¹ K⁻¹), `grad_t` is the temperature gradient
/// (K m⁻¹) vector. Returns acceleration (m s⁻²).
pub fn marangoni_force(dsigma_dt: f64, grad_t: [f64; 3], rho: f64) -> [f64; 3] {
    scale3(grad_t, dsigma_dt / rho)
}

/// Apply Marangoni acceleration to all particles given a uniform temperature
/// gradient field `grad_t` (K m⁻¹).
pub fn apply_marangoni_force(particles: &mut Vec<MicroParticle>, dsigma_dt: f64, grad_t: [f64; 3]) {
    for p in particles.iter_mut() {
        let f_m = marangoni_force(dsigma_dt, grad_t, p.density);
        p.acc = add3(p.acc, f_m);
    }
}

// ---------------------------------------------------------------------------
// SPH density summation and pressure forces
// ---------------------------------------------------------------------------

/// Recompute SPH densities via direct summation W(r_ij, h).
pub fn compute_densities(particles: &mut Vec<MicroParticle>) {
    let n = particles.len();
    let snap: Vec<([f64; 3], f64, f64)> = particles.iter().map(|p| (p.pos, p.mass, p.h)).collect();
    for i in 0..n {
        let mut rho = 0.0f64;
        for j in 0..n {
            let r_ij = sub3(snap[i].0, snap[j].0);
            let r = norm3(r_ij);
            let h_avg = (snap[i].2 + snap[j].2) * 0.5;
            rho += snap[j].1 * wendland_c2(r, h_avg);
        }
        particles[i].density = rho;
    }
}

/// Compute Tait pressures for all particles given reference density `rho0`,
/// reference speed of sound `c_s`, and adiabatic exponent `gamma`.
pub fn compute_pressures(particles: &mut Vec<MicroParticle>, rho0: f64, c_s: f64, gamma: f64) {
    for p in particles.iter_mut() {
        p.pressure = tait_pressure(p.density, rho0, c_s, gamma);
    }
}

/// Apply symmetric SPH pressure gradient acceleration (Morris et al. form).
pub fn apply_pressure_force(particles: &mut Vec<MicroParticle>) {
    let n = particles.len();
    let snap: Vec<([f64; 3], f64, f64, f64, f64)> = particles
        .iter()
        .map(|p| (p.pos, p.pressure, p.density, p.mass, p.h))
        .collect();
    let mut accs: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    for i in 0..n {
        let (pi, rhoi) = (snap[i].1, snap[i].2);
        for j in 0..n {
            if i == j {
                continue;
            }
            let (pj, rhoj, mj) = (snap[j].1, snap[j].2, snap[j].3);
            let r_ij = sub3(snap[i].0, snap[j].0);
            let r = norm3(r_ij);
            let h_avg = (snap[i].4 + snap[j].4) * 0.5;
            let dw = wendland_c2_grad(r, h_avg);
            if r < 1.0e-30 {
                continue;
            }
            let rhat = scale3(r_ij, 1.0 / r);
            let coeff = -mj * (pi / (rhoi * rhoi) + pj / (rhoj * rhoj)) * dw;
            let dv = scale3(rhat, coeff);
            accs[i] = add3(accs[i], dv);
        }
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.acc = add3(p.acc, accs[i]);
    }
}

/// Apply laminar viscosity force (Morris viscosity model).
pub fn apply_viscosity_force(particles: &mut Vec<MicroParticle>, eta: f64) {
    let n = particles.len();
    let snap: Vec<([f64; 3], [f64; 3], f64, f64, f64)> = particles
        .iter()
        .map(|p| (p.pos, p.vel, p.density, p.mass, p.h))
        .collect();
    let mut accs: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    for i in 0..n {
        let (rhoi, hi) = (snap[i].2, snap[i].4);
        for j in 0..n {
            if i == j {
                continue;
            }
            let (rhoj, mj, hj) = (snap[j].2, snap[j].3, snap[j].4);
            let r_ij = sub3(snap[i].0, snap[j].0);
            let r = norm3(r_ij);
            if r < 1.0e-30 {
                continue;
            }
            let h_avg = (hi + hj) * 0.5;
            let dw = wendland_c2_grad(r, h_avg);
            let v_ij = sub3(snap[i].1, snap[j].1);
            let rij_dot_vij = dot3(r_ij, v_ij);
            let coeff =
                2.0 * eta * mj / (rhoi * rhoj) * rij_dot_vij / (r * r + 0.01 * h_avg * h_avg) * dw;
            let dv = scale3(scale3(r_ij, 1.0 / r), coeff);
            accs[i] = add3(accs[i], dv);
        }
    }
    for (i, p) in particles.iter_mut().enumerate() {
        p.acc = add3(p.acc, accs[i]);
    }
}

// ---------------------------------------------------------------------------
// Time integration (Velocity Verlet for microfluidics)
// ---------------------------------------------------------------------------

/// Advance positions and velocities by one Verlet half-step (kick-drift-kick).
///
/// Call order: `verlet_half_kick` → update forces → `verlet_half_kick`.
pub fn verlet_half_kick(particles: &mut Vec<MicroParticle>, dt: f64) {
    for p in particles.iter_mut() {
        p.vel = add3(p.vel, scale3(p.acc, 0.5 * dt));
    }
}

/// Advance positions by a full step using current velocities.
pub fn verlet_drift(particles: &mut Vec<MicroParticle>, dt: f64) {
    for p in particles.iter_mut() {
        p.pos = add3(p.pos, scale3(p.vel, dt));
    }
}

/// Zero all accelerations before re-computing forces each step.
pub fn zero_accelerations(particles: &mut Vec<MicroParticle>) {
    for p in particles.iter_mut() {
        p.acc = [0.0; 3];
    }
}

// ---------------------------------------------------------------------------
// CFL time-step constraint
// ---------------------------------------------------------------------------

/// Compute the CFL-limited time step for an SPH microfluidics simulation.
///
/// dt_CFL = CFL_factor × h_min / (c_s + v_max)
pub fn cfl_timestep(h_min: f64, c_s: f64, v_max: f64, cfl_factor: f64) -> f64 {
    cfl_factor * h_min / (c_s + v_max)
}

/// Compute the viscous time-step constraint.
///
/// dt_visc = CFL_factor × h_min² × rho / eta
pub fn viscous_timestep(h_min: f64, rho: f64, eta: f64, cfl_factor: f64) -> f64 {
    cfl_factor * h_min * h_min * rho / eta
}

// ---------------------------------------------------------------------------
// High-level microfluidics simulation driver
// ---------------------------------------------------------------------------

/// Parameters for a microfluidics SPH simulation.
#[derive(Debug, Clone)]
pub struct MicrofluidicsParams {
    /// Reference density of the continuous phase (kg m⁻³).
    pub rho0: f64,
    /// Speed of sound for Tait EOS (m s⁻¹).
    pub c_s: f64,
    /// Tait adiabatic exponent.
    pub gamma: f64,
    /// Dynamic viscosity of the continuous phase (Pa s).
    pub eta: f64,
    /// Surface tension coefficient at the fluid-fluid interface (N m⁻¹).
    pub sigma: f64,
    /// Equilibrium contact angle with the solid wall (rad).
    pub contact_angle: f64,
    /// CFL safety factor for time stepping.
    pub cfl: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Applied electric field (V m⁻¹) for EOF simulations.
    pub e_field: [f64; 3],
    /// Bulk ion number density for EOF (m⁻³).
    pub n_ion: f64,
    /// Ion valence for EOF.
    pub z_ion: f64,
    /// Relative permittivity of the solvent.
    pub eps_r: f64,
    /// Temperature (K).
    pub temperature: f64,
}

impl MicrofluidicsParams {
    /// Default parameters suitable for a water-droplet-in-oil T-junction.
    pub fn default_water_in_oil() -> Self {
        Self {
            rho0: 1000.0,
            c_s: 1500.0,
            gamma: 7.0,
            eta: ETA_WATER,
            sigma: SIGMA_WATER_AIR * 0.5, // water-in-oil ~0.036 N/m
            contact_angle: PI / 3.0,      // 60°
            cfl: 0.2,
            h: 5.0e-6,
            e_field: [0.0; 3],
            n_ion: 1.0e23,
            z_ion: 1.0,
            eps_r: EPS_WATER,
            temperature: T_REF,
        }
    }
}

/// Run one full SPH micro-step for the microfluidics simulation.
///
/// The sequence is:
/// 1. Half-velocity kick
/// 2. Drift positions
/// 3. Zero accelerations
/// 4. Density summation
/// 5. Pressure computation
/// 6. Apply pressure + viscosity forces
/// 7. Apply color-gradient CSF surface tension
/// 8. Apply EOF body force (if electric field non-zero)
/// 9. Apply wettability boundary condition
/// 10. Half-velocity kick (complete Verlet step)
pub fn micro_step(particles: &mut Vec<MicroParticle>, params: &MicrofluidicsParams, dt: f64) {
    verlet_half_kick(particles, dt);
    verlet_drift(particles, dt);
    zero_accelerations(particles);
    compute_densities(particles);
    compute_pressures(particles, params.rho0, params.c_s, params.gamma);
    apply_pressure_force(particles);
    apply_viscosity_force(particles, params.eta);
    // Interface forces
    compute_color_gradient(particles);
    compute_curvature(particles);
    apply_csf_force(particles, params.sigma);
    // Electro-osmotic body force
    let ef = params.e_field;
    if norm3(ef) > 1.0e-20 {
        for p in particles.iter_mut() {
            p.rho_e =
                boltzmann_charge_density(params.n_ion, params.z_ion, p.phi_e, params.temperature);
        }
        apply_eof_body_force(particles, ef);
    }
    // Wettability
    apply_wettability_bc(particles, params.contact_angle, params.h * 2.0);
    verlet_half_kick(particles, dt);
}

// ---------------------------------------------------------------------------
// Utility: initialise a particle lattice in a rectangular box
// ---------------------------------------------------------------------------

/// Create a uniform lattice of SPH particles filling a box \[0,Lx\]×\[0,Ly\]×\[0,Lz\].
///
/// `nx`, `ny`, `nz` are the number of particles along each axis.
/// `mass` is the particle mass (kg), `h` the smoothing length (m).
/// `phase` is the phase index (0 or 1).
pub fn lattice_particles(
    nx: usize,
    ny: usize,
    nz: usize,
    lx: f64,
    ly: f64,
    lz: f64,
    mass: f64,
    h: f64,
    phase: u8,
) -> Vec<MicroParticle> {
    let mut particles = Vec::with_capacity(nx * ny * nz);
    let dx = lx / nx as f64;
    let dy = ly / ny as f64;
    let dz = lz / nz as f64;
    for ix in 0..nx {
        for iy in 0..ny {
            for iz in 0..nz {
                let pos = [
                    (ix as f64 + 0.5) * dx,
                    (iy as f64 + 0.5) * dy,
                    (iz as f64 + 0.5) * dz,
                ];
                particles.push(MicroParticle::new(pos, mass, h, phase));
            }
        }
    }
    particles
}

/// Add random thermal noise to particle velocities.
///
/// Draws velocity perturbations from a uniform distribution in \[−amp, amp\]
/// for each component.
pub fn add_thermal_noise(particles: &mut Vec<MicroParticle>, amp: f64) {
    let mut rng = rand::rng();
    for p in particles.iter_mut() {
        for v in p.vel.iter_mut() {
            *v += rng.random_range(-amp..amp);
        }
    }
}

// ---------------------------------------------------------------------------
// Flow rate and Reynolds number helpers
// ---------------------------------------------------------------------------

/// Volumetric flow rate from a 2D Poiseuille profile in a rectangular channel.
///
/// Q = (w h³ ΔP) / (12 η L)  (rectangular channel, w >> h approximation).
pub fn poiseuille_flow_rate(width: f64, height: f64, dp: f64, eta: f64, length: f64) -> f64 {
    width * height.powi(3) * dp / (12.0 * eta * length)
}

/// Reynolds number Re = ρ V d_h / η.
pub fn reynolds_number(rho: f64, vel: f64, d_h: f64, eta: f64) -> f64 {
    rho * vel * d_h / eta
}

/// Bond number Bo = ρ g R² / σ (gravitational vs. capillary forces).
pub fn bond_number(rho: f64, g: f64, r: f64, sigma: f64) -> f64 {
    rho * g * r * r / sigma
}

// ---------------------------------------------------------------------------
// Droplet statistics
// ---------------------------------------------------------------------------

/// Compute the centre of mass of a group of particles.
pub fn centre_of_mass(particles: &[MicroParticle]) -> [f64; 3] {
    let mut total_mass = 0.0f64;
    let mut com = [0.0f64; 3];
    for p in particles.iter() {
        total_mass += p.mass;
        com = add3(com, scale3(p.pos, p.mass));
    }
    if total_mass > 0.0 {
        scale3(com, 1.0 / total_mass)
    } else {
        [0.0; 3]
    }
}

/// Compute the mean-square radius of gyration of a droplet particle set.
pub fn radius_of_gyration_sq(particles: &[MicroParticle]) -> f64 {
    let com = centre_of_mass(particles);
    let mut rg2 = 0.0f64;
    let mut total_mass = 0.0f64;
    for p in particles.iter() {
        let dr = sub3(p.pos, com);
        rg2 += p.mass * dot3(dr, dr);
        total_mass += p.mass;
    }
    if total_mass > 0.0 {
        rg2 / total_mass
    } else {
        0.0
    }
}

/// Estimate droplet volume from particle count and reference volume per particle.
pub fn droplet_volume(n_particles: usize, vol_per_particle: f64) -> f64 {
    n_particles as f64 * vol_per_particle
}

/// Estimate equivalent spherical droplet radius from volume.
pub fn equivalent_radius(volume: f64) -> f64 {
    (3.0 * volume / (4.0 * PI)).cbrt()
}

// ---------------------------------------------------------------------------
// Color-gradient interface-tension SPH (Shan–Chen style)
// ---------------------------------------------------------------------------

/// Shan-Chen interaction force between two fluid phases.
///
/// F_SC = −G ρ_A(x) Σ_j ρ_B(x_j) ∇W(r_ij, h) (simplified scalar version).
/// Returns the force per unit volume (Pa m⁻¹).
pub fn shan_chen_force_density(
    rho_a: f64,
    rho_b_sum: f64,
    grad_w_sum: [f64; 3],
    g_sc: f64,
) -> [f64; 3] {
    scale3(grad_w_sum, -g_sc * rho_a * rho_b_sum)
}

/// Effective surface tension from Shan-Chen coupling constant G.
///
/// σ_SC ≈ G ρ₀² h / 6  (approximate for symmetric mixture).
pub fn shan_chen_surface_tension(g_sc: f64, rho0: f64, h: f64) -> f64 {
    g_sc * rho0 * rho0 * h / 6.0
}

// ---------------------------------------------------------------------------
// Electrostatic potential solver (simple iterative Poisson, 1-D)
// ---------------------------------------------------------------------------

/// Solve the linearised Poisson-Boltzmann equation in 1-D using a finite-
/// difference Gauss-Seidel iteration.
///
/// The domain has `n` nodes from 0 to L, with Dirichlet BCs `phi_left` and
/// `phi_right`. `kappa2` = 1/λ_D² (m⁻²). Returns the potential array (V).
pub fn solve_poisson_boltzmann_1d(
    n: usize,
    length: f64,
    phi_left: f64,
    phi_right: f64,
    kappa2: f64,
    n_iter: usize,
) -> Vec<f64> {
    let dx = length / (n as f64 - 1.0);
    let mut phi = vec![0.0f64; n];
    phi[0] = phi_left;
    phi[n - 1] = phi_right;
    for _ in 0..n_iter {
        for i in 1..n - 1 {
            phi[i] = (phi[i - 1] + phi[i + 1]) / (2.0 + kappa2 * dx * dx);
        }
    }
    phi
}

// ---------------------------------------------------------------------------
// Capillary wave spectrum
// ---------------------------------------------------------------------------

/// Capillary wave dispersion relation: ω² = σ k³ / ρ  (deep-fluid limit).
///
/// Returns angular frequency (rad s⁻¹) for wavenumber `k` (m⁻¹).
pub fn capillary_wave_dispersion(sigma: f64, rho: f64, k: f64) -> f64 {
    (sigma * k.powi(3) / rho).sqrt()
}

/// Thermal roughness of a liquid interface: ξ² = k_B T / (σ k_max² · A)
///
/// `k_max` upper wavenumber cutoff (m⁻¹), `area` interface area (m²),
/// Returns mean-square roughness (m²).
pub fn thermal_interface_roughness(temp: f64, sigma: f64, k_max: f64, area: f64) -> f64 {
    K_B * temp / (sigma * k_max * k_max * area)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wendland_c2_zero_distance() {
        let w = wendland_c2(0.0, 1.0e-6);
        assert!(w > 0.0);
    }

    #[test]
    fn test_wendland_c2_beyond_support() {
        let w = wendland_c2(3.0e-6, 1.0e-6);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_wendland_c2_grad_zero_distance() {
        let dw = wendland_c2_grad(0.0, 1.0e-6);
        assert_eq!(dw, 0.0);
    }

    #[test]
    fn test_wendland_c2_grad_beyond_support() {
        let dw = wendland_c2_grad(3.0e-6, 1.0e-6);
        assert_eq!(dw, 0.0);
    }

    #[test]
    fn test_tait_pressure_at_reference_density() {
        let p = tait_pressure(1000.0, 1000.0, 1500.0, 7.0);
        assert!(p.abs() < 1.0e-6);
    }

    #[test]
    fn test_tait_pressure_increases_with_density() {
        let p1 = tait_pressure(1010.0, 1000.0, 1500.0, 7.0);
        let p2 = tait_pressure(1000.0, 1000.0, 1500.0, 7.0);
        assert!(p1 > p2);
    }

    #[test]
    fn test_laplace_pressure_sphere() {
        let dp = laplace_pressure_sphere(SIGMA_WATER_AIR, 1.0e-6);
        assert!((dp - 2.0 * SIGMA_WATER_AIR / 1.0e-6).abs() < 1.0e-6);
    }

    #[test]
    fn test_laplace_pressure_cylinder_half_sphere() {
        // Cylinder should give half the sphere pressure for same R
        let sphere = laplace_pressure_sphere(0.05, 50.0e-6);
        let cyl = laplace_pressure_cylinder(0.05, 50.0e-6);
        assert!((sphere - 2.0 * cyl).abs() < 1.0e-10);
    }

    #[test]
    fn test_young_contact_angle_hydrophilic() {
        // σ_SG > σ_SL → cos θ > 0 → θ < 90°
        let theta = young_contact_angle(0.07, 0.03, 0.072);
        assert!(theta < PI / 2.0);
    }

    #[test]
    fn test_young_contact_angle_hydrophobic() {
        // σ_SL > σ_SG → cos θ < 0 → θ > 90°
        let theta = young_contact_angle(0.02, 0.07, 0.072);
        assert!(theta > PI / 2.0);
    }

    #[test]
    fn test_capillary_number_water() {
        let ca = capillary_number(ETA_WATER, 1.0e-3, SIGMA_WATER_AIR);
        assert!(ca > 0.0 && ca < 1.0e-3);
    }

    #[test]
    fn test_droplet_breaks_at_high_ca() {
        assert!(droplet_breaks(0.5));
    }

    #[test]
    fn test_droplet_does_not_break_at_low_ca() {
        assert!(!droplet_breaks(0.1));
    }

    #[test]
    fn test_debye_length_order_of_magnitude() {
        // ~10 nm for 0.1 M monovalent salt in water at 298 K
        let lam = debye_length(0.1 * 1000.0 * N_AV, 1.0, EPS_WATER, T_REF);
        assert!(lam > 1.0e-10 && lam < 1.0e-7);
    }

    #[test]
    fn test_eof_slip_velocity_sign() {
        // Negative zeta (most surfaces), forward E → forward flow
        let v = eof_slip_velocity(-50.0e-3, 1000.0, EPS_WATER, ETA_WATER);
        assert!(v > 0.0);
    }

    #[test]
    fn test_boltzmann_charge_density_zero_at_bulk() {
        let rho_e = boltzmann_charge_density(1.0e23, 1.0, 0.0, T_REF);
        assert_eq!(rho_e, 0.0);
    }

    #[test]
    fn test_clausius_mossotti_positive_dielectric_contrast() {
        let k = clausius_mossotti_real(80.0, 2.5);
        assert!(k > 0.0);
    }

    #[test]
    fn test_clausius_mossotti_negative_dielectric_contrast() {
        let k = clausius_mossotti_real(1.5, 80.0);
        assert!(k < 0.0);
    }

    #[test]
    fn test_dep_force_scales_with_radius_cubed() {
        let f1 = dep_force(1.0e-6, EPS_WATER, 0.5, [1.0e12, 0.0, 0.0]);
        let f2 = dep_force(2.0e-6, EPS_WATER, 0.5, [1.0e12, 0.0, 0.0]);
        // Volume ∝ r³ → force ratio should be 8×
        assert!((f2[0] / f1[0] - 8.0).abs() < 1.0e-9);
    }

    #[test]
    fn test_ewod_contact_angle_reduces_with_voltage() {
        let theta0 = PI / 2.0;
        let theta_v = ewod_contact_angle(theta0, 1.0e-4, 50.0, SIGMA_WATER_AIR);
        assert!(theta_v < theta0);
    }

    #[test]
    fn test_ewod_driving_pressure_positive() {
        let dp = ewod_driving_pressure(SIGMA_WATER_AIR, PI / 4.0, PI / 2.0, 100.0e-6);
        assert!(dp > 0.0);
    }

    #[test]
    fn test_peclet_number_magnitude() {
        // Water flow at 1 mm/s, 100 µm channel, D_dye ~ 1e-9 m²/s
        let pe = peclet_number(1.0e-3, 100.0e-6, 1.0e-9);
        assert!((pe - 100.0).abs() < 1.0e-9);
    }

    #[test]
    fn test_ohnesorge_number_water_droplet() {
        let oh = ohnesorge_number(ETA_WATER, 1000.0, SIGMA_WATER_AIR, 50.0e-6);
        assert!(oh > 0.0 && oh < 0.1);
    }

    #[test]
    fn test_weber_number_scales_with_velocity_squared() {
        let we1 = weber_number(1000.0, 1.0, 50.0e-6, SIGMA_WATER_AIR);
        let we2 = weber_number(1000.0, 2.0, 50.0e-6, SIGMA_WATER_AIR);
        assert!((we2 / we1 - 4.0).abs() < 1.0e-9);
    }

    #[test]
    fn test_lattice_particles_count() {
        let parts = lattice_particles(4, 4, 4, 100.0e-6, 100.0e-6, 100.0e-6, 1.0e-12, 5.0e-6, 0);
        assert_eq!(parts.len(), 64);
    }

    #[test]
    fn test_lattice_particles_position_range() {
        let parts = lattice_particles(3, 3, 3, 30.0e-6, 30.0e-6, 30.0e-6, 1.0e-12, 5.0e-6, 0);
        for p in &parts {
            assert!(p.pos[0] >= 0.0 && p.pos[0] <= 30.0e-6);
            assert!(p.pos[1] >= 0.0 && p.pos[1] <= 30.0e-6);
            assert!(p.pos[2] >= 0.0 && p.pos[2] <= 30.0e-6);
        }
    }

    #[test]
    fn test_centre_of_mass_uniform() {
        let parts = lattice_particles(2, 2, 2, 20.0e-6, 20.0e-6, 20.0e-6, 1.0e-12, 5.0e-6, 0);
        let com = centre_of_mass(&parts);
        // Symmetric lattice → com should be at (10,10,10) µm
        for &c in com.iter() {
            assert!((c - 10.0e-6).abs() < 1.0e-9);
        }
    }

    #[test]
    fn test_capillary_wave_dispersion_positive() {
        let omega = capillary_wave_dispersion(SIGMA_WATER_AIR, 1000.0, 1.0e6);
        assert!(omega > 0.0);
    }

    #[test]
    fn test_thermal_interface_roughness_positive() {
        let xi2 = thermal_interface_roughness(T_REF, SIGMA_WATER_AIR, 1.0e9, 1.0e-12);
        assert!(xi2 > 0.0);
    }

    #[test]
    fn test_shan_chen_surface_tension_positive() {
        let s = shan_chen_surface_tension(0.1, 1000.0, 5.0e-6);
        assert!(s > 0.0);
    }

    #[test]
    fn test_poiseuille_flow_rate_increases_with_dp() {
        let q1 = poiseuille_flow_rate(100.0e-6, 50.0e-6, 100.0, ETA_WATER, 1.0e-3);
        let q2 = poiseuille_flow_rate(100.0e-6, 50.0e-6, 200.0, ETA_WATER, 1.0e-3);
        assert!(q2 > q1);
    }

    #[test]
    fn test_reynolds_number_microfluidics_low() {
        // Typical micro Re << 1
        let re = reynolds_number(1000.0, 1.0e-3, 100.0e-6, ETA_WATER);
        assert!(re < 1.0);
    }

    #[test]
    fn test_bond_number_small_droplet() {
        // Bo << 1 for micro-droplet
        let bo = bond_number(1000.0, 9.81, 50.0e-6, SIGMA_WATER_AIR);
        assert!(bo < 0.01);
    }

    #[test]
    fn test_cfl_timestep_positive() {
        let dt = cfl_timestep(5.0e-6, 1500.0, 0.1, 0.2);
        assert!(dt > 0.0);
    }

    #[test]
    fn test_compute_densities_nonzero() {
        let mut parts = lattice_particles(3, 3, 3, 30.0e-6, 30.0e-6, 30.0e-6, 1.0e-12, 5.0e-6, 0);
        compute_densities(&mut parts);
        let rho_sum: f64 = parts.iter().map(|p| p.density).sum();
        assert!(rho_sum > 0.0);
    }

    #[test]
    fn test_split_droplet_cloud_count() {
        let parts = lattice_particles(2, 2, 2, 20.0e-6, 20.0e-6, 20.0e-6, 1.0e-12, 5.0e-6, 1);
        let (left, right) = split_droplet_cloud(&parts, 5.0e-6);
        assert_eq!(left.len() + right.len(), parts.len());
    }

    #[test]
    fn test_radius_of_gyration_sq_positive() {
        let parts = lattice_particles(3, 3, 3, 30.0e-6, 30.0e-6, 30.0e-6, 1.0e-12, 5.0e-6, 0);
        let rg2 = radius_of_gyration_sq(&parts);
        assert!(rg2 > 0.0);
    }

    #[test]
    fn test_droplet_volume_and_radius_consistent() {
        let vol = droplet_volume(1000, 1.0e-18);
        let r = equivalent_radius(vol);
        // Back-compute volume from r
        let vol2 = 4.0 / 3.0 * PI * r.powi(3);
        assert!((vol - vol2).abs() / vol < 1.0e-10);
    }

    #[test]
    fn test_marangoni_force_direction() {
        // dσ/dT < 0 (typical), ∇T in +x → force in −x
        let f = marangoni_force(-1.0e-4, [1.0e3, 0.0, 0.0], 1000.0);
        assert!(f[0] < 0.0);
    }

    #[test]
    fn test_t_junction_droplet_length_increases_with_flow_ratio() {
        let l1 = t_junction_droplet_length(100.0e-6, 0.5, 1.0, 1.0);
        let l2 = t_junction_droplet_length(100.0e-6, 1.0, 1.0, 1.0);
        assert!(l2 > l1);
    }

    #[test]
    fn test_herringbone_mixing_length_positive() {
        let l = herringbone_mixing_length(100.0e-6, 8.0, 0.3);
        assert!(l > 0.0);
    }

    #[test]
    fn test_add_thermal_noise_changes_velocities() {
        let mut parts = lattice_particles(2, 2, 1, 20.0e-6, 20.0e-6, 10.0e-6, 1.0e-12, 5.0e-6, 0);
        add_thermal_noise(&mut parts, 1.0e-3);
        let total_v: f64 = parts.iter().map(|p| norm3(p.vel)).sum();
        assert!(total_v > 0.0);
    }
}
