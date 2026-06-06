// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular rigid body simulation.
//!
//! This module implements a comprehensive Discrete Element Method (DEM) solver
//! for granular assemblies of rigid particles. Topics covered:
//!
//! - Hertz-Mindlin contact mechanics (normal + tangential + rolling)
//! - Particle angularity via shape factor
//! - Grain packing (random close packing, FCC/HCP lattice placement)
//! - Avalanche dynamics and angle of repose measurement
//! - Shear band formation detection via vorticity
//! - Compaction simulation under cyclic loading
//! - Grain crushing model (energy-based split criterion)
//! - Cohesion forces (capillary bridges, van der Waals)
//! - Size-based segregation in granular flow
//! - Force chain network extraction and analysis
//! - Granular temperature (velocity fluctuation measure)
//! - Janssen arching effect (wall friction pressure saturation)

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Inline vector arithmetic helpers
// ---------------------------------------------------------------------------

#[inline]
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = norm(a);
    if n < 1e-15 {
        [0.0; 3]
    } else {
        scale(a, 1.0 / n)
    }
}

#[inline]
fn norm_sq(a: [f64; 3]) -> f64 {
    dot(a, a)
}

// ---------------------------------------------------------------------------
// GrainShape — angularity model
// ---------------------------------------------------------------------------

/// Shape descriptor for a granular particle.
///
/// A shape factor `f_shape ∈ [1, ∞)` modifies effective contact stiffness
/// and rolling resistance. A sphere has `f_shape = 1.0`.
#[derive(Debug, Clone)]
pub struct GrainShape {
    /// Shape factor (1.0 = perfect sphere, >1 = angular).
    pub factor: f64,
    /// Aspect ratio (major/minor semi-axis length).
    pub aspect_ratio: f64,
    /// Roundness index in \[0, 1\] (1 = perfectly round corners).
    pub roundness: f64,
}

impl GrainShape {
    /// Perfectly spherical grain.
    pub fn sphere() -> Self {
        Self {
            factor: 1.0,
            aspect_ratio: 1.0,
            roundness: 1.0,
        }
    }

    /// Angular grain with the given parameters.
    pub fn angular(factor: f64, aspect_ratio: f64, roundness: f64) -> Self {
        Self {
            factor: factor.max(1.0),
            aspect_ratio: aspect_ratio.max(1.0),
            roundness: roundness.clamp(0.0, 1.0),
        }
    }

    /// Effective rolling resistance multiplier from shape.
    pub fn rolling_resistance_multiplier(&self) -> f64 {
        self.factor * (2.0 - self.roundness)
    }

    /// Effective contact area multiplier (larger for angular particles).
    pub fn contact_area_multiplier(&self) -> f64 {
        1.0 / self.roundness.max(0.01)
    }
}

// ---------------------------------------------------------------------------
// GranularParticle
// ---------------------------------------------------------------------------

/// A single granular particle with rigid-body kinematics.
#[derive(Debug, Clone)]
pub struct GranularParticle {
    /// Centre-of-mass position \[m\].
    pub position: [f64; 3],
    /// Linear velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Accumulated force \[N\].
    pub force: [f64; 3],
    /// Accumulated torque \[N·m\].
    pub torque: [f64; 3],
    /// Nominal radius \[m\] (used for contact detection).
    pub radius: f64,
    /// Mass \[kg\].
    pub mass: f64,
    /// Moment of inertia (solid sphere: 2/5 m r²).
    pub inertia: f64,
    /// Young's modulus \[Pa\].
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Coulomb friction coefficient.
    pub friction: f64,
    /// Rolling resistance coefficient.
    pub rolling_coeff: f64,
    /// Shape descriptor.
    pub shape: GrainShape,
    /// Whether this grain has been crushed.
    pub crushed: bool,
    /// Accumulated strain energy \[J\] (for crushing criterion).
    pub strain_energy: f64,
    /// Particle index.
    pub id: usize,
}

impl GranularParticle {
    /// Create a new spherical grain with density-derived mass.
    pub fn new(id: usize, radius: f64, density: f64) -> Self {
        let mass = (4.0 / 3.0) * PI * radius.powi(3) * density;
        let inertia = (2.0 / 5.0) * mass * radius * radius;
        Self {
            position: [0.0; 3],
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            force: [0.0; 3],
            torque: [0.0; 3],
            radius,
            mass,
            inertia,
            young: 1.0e8,
            poisson: 0.3,
            restitution: 0.7,
            friction: 0.4,
            rolling_coeff: 0.01,
            shape: GrainShape::sphere(),
            crushed: false,
            strain_energy: 0.0,
            id,
        }
    }

    /// Reset force and torque accumulators to zero.
    pub fn reset_accumulators(&mut self) {
        self.force = [0.0; 3];
        self.torque = [0.0; 3];
    }

    /// Kinetic energy of this particle (linear only).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * norm_sq(self.velocity)
    }

    /// Rotational kinetic energy.
    pub fn rotational_energy(&self) -> f64 {
        0.5 * self.inertia * norm_sq(self.angular_velocity)
    }
}

// ---------------------------------------------------------------------------
// ContactPair — Hertz-Mindlin state
// ---------------------------------------------------------------------------

/// Persistent contact state between two particles.
#[derive(Debug, Clone)]
pub struct ContactPair {
    /// Index of particle i.
    pub i: usize,
    /// Index of particle j.
    pub j: usize,
    /// Tangential spring displacement (Mindlin model) \[m\].
    pub tangential_spring: [f64; 3],
    /// Current normal overlap δ \[m\].
    pub overlap: f64,
    /// Normal contact force magnitude \[N\].
    pub fn_mag: f64,
}

impl ContactPair {
    /// Create a new contact pair.
    pub fn new(i: usize, j: usize) -> Self {
        Self {
            i,
            j,
            tangential_spring: [0.0; 3],
            overlap: 0.0,
            fn_mag: 0.0,
        }
    }
}

// ---------------------------------------------------------------------------
// Hertz normal contact force
// ---------------------------------------------------------------------------

/// Compute Hertz normal contact force magnitude.
///
/// F_n = (4/3) E* √R* · δ^(3/2)
///
/// # Arguments
/// - `overlap` – penetration depth δ \[m\], must be ≥ 0.
/// - `e_star`  – effective Young's modulus E* \[Pa\].
/// - `r_star`  – effective radius R* \[m\].
pub fn hertz_normal_force(overlap: f64, e_star: f64, r_star: f64) -> f64 {
    if overlap <= 0.0 {
        return 0.0;
    }
    (4.0 / 3.0) * e_star * r_star.sqrt() * overlap.powf(1.5)
}

// ---------------------------------------------------------------------------
// Mindlin tangential force
// ---------------------------------------------------------------------------

/// Compute the scalar Mindlin tangential spring force, capped by Coulomb friction.
///
/// k_t = 8 G* √R* |δ_t|^(1/2)  (incremental Mindlin).
/// Force = k_t · |δ_t|, direction is opposing slip.
///
/// # Arguments
/// - `ts_mag`      – magnitude of tangential spring displacement |δ_t| \[m\].
/// - `g_star`      – effective shear modulus G* \[Pa\].
/// - `r_star`      – effective radius R* \[m\].
/// - `fn_mag`      – normal force magnitude \[N\].
/// - `mu`          – friction coefficient.
pub fn mindlin_tangential_force(
    ts_mag: f64,
    g_star: f64,
    r_star: f64,
    fn_mag: f64,
    mu: f64,
) -> f64 {
    if ts_mag < 1e-15 {
        return 0.0;
    }
    let kt = 8.0 * g_star * r_star.sqrt() * ts_mag.sqrt().max(1e-15);
    let ft = kt * ts_mag;
    ft.min(mu * fn_mag)
}

// ---------------------------------------------------------------------------
// Rolling resistance torque
// ---------------------------------------------------------------------------

/// Compute rolling resistance torque magnitude.
///
/// τ_r = μ_r · R_eff · F_n · κ_shape
///
/// where κ_shape is the shape rolling resistance multiplier.
///
/// # Arguments
/// - `mu_r`   – rolling resistance coefficient.
/// - `r_eff`  – effective radius \[m\].
/// - `fn_mag` – normal force magnitude \[N\].
/// - `kappa`  – shape multiplier (1.0 for sphere).
pub fn rolling_resistance_torque(mu_r: f64, r_eff: f64, fn_mag: f64, kappa: f64) -> f64 {
    mu_r * r_eff * fn_mag * kappa
}

// ---------------------------------------------------------------------------
// Capillary cohesion force
// ---------------------------------------------------------------------------

/// Capillary bridge force between two particles in a wet granular assembly.
///
/// Uses the Willett (2000) toroidal approximation:
/// F_cap = π γ R* (1 + tanh(1.9 s/V^(1/3)))^(-1)
///
/// For simplicity this function returns a linear approximation
/// F_cap ≈ -2π γ R* / (1 + s / l_cap) where l_cap = (V/π)^(1/3).
///
/// # Arguments
/// - `surface_tension` – liquid-air surface tension γ \[N/m\].
/// - `r_star`          – effective radius R* \[m\].
/// - `separation`      – surface-to-surface gap s \[m\] (0 = touching).
/// - `bridge_volume`   – liquid bridge volume V \[m³\].
pub fn capillary_force(
    surface_tension: f64,
    r_star: f64,
    separation: f64,
    bridge_volume: f64,
) -> f64 {
    if bridge_volume < 1e-30 {
        return 0.0;
    }
    let l_cap = (bridge_volume / PI).cbrt();
    let denom = 1.0 + separation / l_cap.max(1e-15);
    2.0 * PI * surface_tension * r_star / denom
}

// ---------------------------------------------------------------------------
// Van der Waals force
// ---------------------------------------------------------------------------

/// Van der Waals adhesion force between two spheres (Hamaker model).
///
/// F_vdW = -A_H R* / (6 h²)
///
/// where A_H is the Hamaker constant and h is the surface gap.
///
/// # Arguments
/// - `hamaker`    – Hamaker constant A_H \[J\].
/// - `r_star`     – effective radius R* \[m\].
/// - `gap`        – surface-to-surface separation h \[m\] (clamped to min_gap).
/// - `min_gap`    – minimum gap to prevent divergence \[m\].
pub fn van_der_waals_force(hamaker: f64, r_star: f64, gap: f64, min_gap: f64) -> f64 {
    let h = gap.max(min_gap);
    hamaker * r_star / (6.0 * h * h)
}

// ---------------------------------------------------------------------------
// Full Hertz-Mindlin contact resolution
// ---------------------------------------------------------------------------

/// Output of a single Hertz-Mindlin contact resolution.
#[derive(Debug, Clone)]
pub struct ContactResult {
    /// Force on particle i \[N\].
    pub force_i: [f64; 3],
    /// Force on particle j \[N\].
    pub force_j: [f64; 3],
    /// Torque on particle i \[N·m\].
    pub torque_i: [f64; 3],
    /// Torque on particle j \[N·m\].
    pub torque_j: [f64; 3],
    /// Updated tangential spring displacement \[m\].
    pub new_tangential_spring: [f64; 3],
    /// Normal force magnitude \[N\].
    pub fn_mag: f64,
    /// Strain energy increment \[J\] (for crushing model).
    pub strain_energy: f64,
}

/// Resolve a Hertz-Mindlin contact between two particles.
///
/// Computes normal, tangential, rolling resistance, and optional cohesion
/// forces/torques. The tangential spring is updated with the current time step.
pub fn resolve_hertz_mindlin(
    pi: &GranularParticle,
    pj: &GranularParticle,
    tangential_spring: [f64; 3],
    dt: f64,
    capillary_params: Option<(f64, f64)>, // (surface_tension, bridge_volume)
    hamaker: Option<f64>,
) -> Option<ContactResult> {
    let rij = sub(pj.position, pi.position);
    let d = norm(rij);
    let overlap = pi.radius + pj.radius - d;

    if overlap <= 0.0 {
        // Check cohesion range even when separated
        let cohesive = match (capillary_params, hamaker) {
            (Some((gamma, v)), _) => {
                let gap = -overlap;
                let r_star = pi.radius * pj.radius / (pi.radius + pj.radius);
                let fc = capillary_force(gamma, r_star, gap, v);
                if fc < 1e-15 {
                    return None;
                }
                fc
            }
            (_, Some(ah)) => {
                let gap = (-overlap).max(1e-10);
                let r_star = pi.radius * pj.radius / (pi.radius + pj.radius);
                van_der_waals_force(ah, r_star, gap, 1e-10)
            }
            _ => return None,
        };
        let n = normalize(rij);
        let fi = scale(n, cohesive);
        let fj = scale(n, -cohesive);
        return Some(ContactResult {
            force_i: fi,
            force_j: fj,
            torque_i: [0.0; 3],
            torque_j: [0.0; 3],
            new_tangential_spring: tangential_spring,
            fn_mag: cohesive,
            strain_energy: 0.0,
        });
    }

    let n = normalize(rij);

    // Effective modulus and radius
    let ei = pi.young;
    let ej = pj.young;
    let vi = pi.poisson;
    let vj = pj.poisson;
    let e_star = 1.0 / ((1.0 - vi * vi) / ei + (1.0 - vj * vj) / ej);
    let r_star = pi.radius * pj.radius / (pi.radius + pj.radius);
    let g_star = 1.0 / (2.0 * (2.0 - vi) * (1.0 + vi) / ei + 2.0 * (2.0 - vj) * (1.0 + vj) / ej);

    // Normal force
    let fn_elastic = hertz_normal_force(overlap, e_star, r_star);

    // Viscous normal damping
    let rel_v = sub(pi.velocity, pj.velocity);
    let vn = dot(rel_v, n);
    let beta_n = 1.0 - pi.restitution.min(pj.restitution);
    let fn_damp = beta_n * fn_elastic.sqrt() * vn * 0.1;
    let fn_total = (fn_elastic - fn_damp).max(0.0);

    let fn_vec = scale(n, -fn_total);

    // Relative tangential velocity at contact point
    let omega_sum = add(
        cross(pi.angular_velocity, scale(n, pi.radius)),
        cross(pj.angular_velocity, scale(n, -pj.radius)),
    );
    let v_rel = sub(rel_v, omega_sum);
    let vn_vec = scale(n, dot(v_rel, n));
    let vt_vec = sub(v_rel, vn_vec);

    // Update tangential spring
    let new_ts = add(tangential_spring, scale(vt_vec, dt));
    let ts_mag = norm(new_ts);
    let ft_mag = mindlin_tangential_force(
        ts_mag,
        g_star,
        r_star,
        fn_total,
        pi.friction.min(pj.friction),
    );
    let ft_dir = if ts_mag > 1e-15 {
        scale(new_ts, -1.0 / ts_mag)
    } else {
        [0.0; 3]
    };
    let ft_vec = scale(ft_dir, ft_mag);

    // Cap tangential spring to Coulomb limit
    let mu = pi.friction.min(pj.friction);
    let ts_limit = mu * fn_total / (8.0 * g_star * r_star.sqrt() * ts_mag.sqrt().max(1e-15));
    let new_ts_capped = if ts_mag > ts_limit {
        scale(new_ts, ts_limit / ts_mag.max(1e-15))
    } else {
        new_ts
    };

    // Torques: r × F at contact point
    let ri = scale(n, pi.radius);
    let rj = scale(n, -pj.radius);
    let torque_i_ft = cross(ri, ft_vec);
    let torque_j_ft = cross(rj, scale(ft_vec, -1.0));

    // Rolling resistance
    let kappa_i = pi.shape.rolling_resistance_multiplier();
    let kappa_j = pj.shape.rolling_resistance_multiplier();
    let kappa = (kappa_i + kappa_j) * 0.5;
    let rr_mag = rolling_resistance_torque(
        pi.rolling_coeff.min(pj.rolling_coeff),
        r_star,
        fn_total,
        kappa,
    );
    let omega_rel = sub(pi.angular_velocity, pj.angular_velocity);
    let om_norm = norm(omega_rel);
    let rr_torque_i = if om_norm > 1e-15 {
        scale(omega_rel, -rr_mag / om_norm)
    } else {
        [0.0; 3]
    };

    // Cohesion addition
    let fc_add = match (capillary_params, hamaker) {
        (Some((gamma, v)), _) => capillary_force(gamma, r_star, 0.0, v),
        (_, Some(ah)) => van_der_waals_force(ah, r_star, 1e-10, 1e-10),
        _ => 0.0,
    };
    let fn_cohesive = scale(n, fc_add);

    let fi_total = add(add(fn_vec, ft_vec), fn_cohesive);
    let fj_total = [
        -fi_total[0] - 2.0 * ft_vec[0],
        -fi_total[1] - 2.0 * ft_vec[1],
        -fi_total[2] - 2.0 * ft_vec[2],
    ];

    let ti = add(torque_i_ft, rr_torque_i);
    let tj = torque_j_ft;

    let strain_energy = 0.4 * fn_total * overlap;

    Some(ContactResult {
        force_i: fi_total,
        force_j: fj_total,
        torque_i: ti,
        torque_j: tj,
        new_tangential_spring: new_ts_capped,
        fn_mag: fn_total,
        strain_energy,
    })
}

// ---------------------------------------------------------------------------
// Grain crushing model
// ---------------------------------------------------------------------------

/// Grain crushing criterion (energy-based).
///
/// A grain crushes when its accumulated strain energy exceeds a critical
/// threshold: E_crush = σ_t² V / (2 E) where σ_t is tensile strength.
///
/// Returns `true` if the particle should be split.
///
/// # Arguments
/// - `strain_energy`  – accumulated contact strain energy \[J\].
/// - `tensile_strength` – material tensile strength \[Pa\].
/// - `young`          – Young's modulus \[Pa\].
/// - `radius`         – particle radius \[m\].
pub fn crushing_criterion(
    strain_energy: f64,
    tensile_strength: f64,
    young: f64,
    radius: f64,
) -> bool {
    let volume = (4.0 / 3.0) * PI * radius.powi(3);
    let e_crit = tensile_strength * tensile_strength * volume / (2.0 * young);
    strain_energy >= e_crit
}

/// Split a grain into two daughter grains of half-volume (radius reduced by 2^(-1/3)).
///
/// The daughters are placed symmetrically about the parent centre along `split_dir`.
///
/// Returns `(daughter_a, daughter_b)`.
pub fn split_grain(
    parent: &GranularParticle,
    split_dir: [f64; 3],
) -> (GranularParticle, GranularParticle) {
    let new_radius = parent.radius / 2.0_f64.cbrt();
    let density = parent.mass / ((4.0 / 3.0) * PI * parent.radius.powi(3));
    let offset = normalize(split_dir);
    let d = new_radius * 0.95;

    let mut da = GranularParticle::new(parent.id * 2, new_radius, density);
    da.position = add(parent.position, scale(offset, d));
    da.velocity = parent.velocity;
    da.young = parent.young;
    da.poisson = parent.poisson;
    da.restitution = parent.restitution;
    da.friction = parent.friction;
    da.rolling_coeff = parent.rolling_coeff;
    da.shape = parent.shape.clone();

    let mut db = GranularParticle::new(parent.id * 2 + 1, new_radius, density);
    db.position = sub(parent.position, scale(offset, d));
    db.velocity = parent.velocity;
    db.young = parent.young;
    db.poisson = parent.poisson;
    db.restitution = parent.restitution;
    db.friction = parent.friction;
    db.rolling_coeff = parent.rolling_coeff;
    db.shape = parent.shape.clone();

    (da, db)
}

// ---------------------------------------------------------------------------
// Grain packing: FCC / HCP / random
// ---------------------------------------------------------------------------

/// Generate a face-centred cubic (FCC) packing within a bounding box.
///
/// The packing fills the box `[0, Lx] × [0, Ly] × [0, Lz]` with particles
/// of the given `radius`. Returns particle positions.
pub fn fcc_packing(radius: f64, lx: f64, ly: f64, lz: f64) -> Vec<[f64; 3]> {
    let a = radius * 2.0_f64.sqrt() * 2.0; // FCC lattice constant
    let mut positions = Vec::new();
    let nx = (lx / a).ceil() as usize + 1;
    let ny = (ly / a).ceil() as usize + 1;
    let nz = (lz / a).ceil() as usize + 1;
    // FCC basis vectors
    let basis: [[f64; 3]; 4] = [
        [0.0, 0.0, 0.0],
        [0.5 * a, 0.5 * a, 0.0],
        [0.5 * a, 0.0, 0.5 * a],
        [0.0, 0.5 * a, 0.5 * a],
    ];
    for iz in 0..nz {
        for iy in 0..ny {
            for ix in 0..nx {
                for &b in &basis {
                    let pos = [
                        ix as f64 * a + b[0] + radius,
                        iy as f64 * a + b[1] + radius,
                        iz as f64 * a + b[2] + radius,
                    ];
                    if pos[0] < lx && pos[1] < ly && pos[2] < lz {
                        positions.push(pos);
                    }
                }
            }
        }
    }
    positions
}

/// Generate a hexagonal close-packed (HCP) column packing.
///
/// Returns a list of `(x, y, z)` particle centres for 2-D layers stacked in Z.
pub fn hcp_packing(radius: f64, lx: f64, ly: f64, lz: f64) -> Vec<[f64; 3]> {
    let d = 2.0 * radius;
    let dx_row = d;
    let dy_layer = d * (3.0_f64.sqrt() / 2.0);
    let dz_stack = d * (2.0_f64 / 3.0).sqrt() * 2.0;
    let mut positions = Vec::new();
    let mut iz = 0usize;
    let mut z = radius;
    while z < lz {
        let mut iy = 0usize;
        let mut y = radius;
        while y < ly {
            let x_offset = if (iy + iz).is_multiple_of(2) {
                0.0
            } else {
                radius
            };
            let mut x = radius + x_offset;
            while x < lx {
                positions.push([x, y, z]);
                x += dx_row;
            }
            y += dy_layer;
            iy += 1;
        }
        z += dz_stack;
        iz += 1;
    }
    positions
}

/// Compute the packing fraction from a list of particle radii and domain volume.
pub fn packing_fraction(radii: &[f64], domain_volume: f64) -> f64 {
    if domain_volume < 1e-30 {
        return 0.0;
    }
    let vp: f64 = radii.iter().map(|&r| (4.0 / 3.0) * PI * r.powi(3)).sum();
    (vp / domain_volume).min(1.0)
}

// ---------------------------------------------------------------------------
// Granular temperature
// ---------------------------------------------------------------------------

/// Compute the granular temperature of an assembly.
///
/// T_g = (1 / 3N) Σ m_i (v_i - `v`)² / m_i = (1/3) <(δv)²>
///
/// where δv = v_i − `v` is the fluctuation velocity.
///
/// # Arguments
/// - `particles` – slice of particles.
///
/// Returns T_g in \[m²/s²\].
pub fn granular_temperature(particles: &[GranularParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let n = particles.len() as f64;
    let mean_v: [f64; 3] = {
        let sum = particles
            .iter()
            .fold([0.0_f64; 3], |acc, p| add(acc, p.velocity));
        scale(sum, 1.0 / n)
    };
    let var: f64 = particles
        .iter()
        .map(|p| norm_sq(sub(p.velocity, mean_v)))
        .sum::<f64>();
    var / (3.0 * n)
}

// ---------------------------------------------------------------------------
// Angle of repose
// ---------------------------------------------------------------------------

/// Estimate the static angle of repose from the fabric tensor's principal directions.
///
/// A simple geometric estimate: θ_repose ≈ arctan(μ + μ_r * κ_shape)
/// for a single grain resting on a slope.
///
/// # Arguments
/// - `friction`          – Coulomb friction coefficient μ.
/// - `rolling_coeff`     – Rolling resistance coefficient μ_r.
/// - `shape_factor`      – Grain angularity multiplier κ.
///
/// Returns angle in radians.
pub fn angle_of_repose_estimate(friction: f64, rolling_coeff: f64, shape_factor: f64) -> f64 {
    let effective_mu = friction + rolling_coeff * shape_factor;
    effective_mu.atan()
}

/// Measure the actual angle of repose from a heap of particles.
///
/// Fits a plane to the outer surface of the heap and returns the angle
/// between the surface normal and the vertical.
///
/// This simplified version finds the maximum slope angle from the base centre.
pub fn measure_angle_of_repose(particles: &[GranularParticle], base_y: f64) -> f64 {
    if particles.len() < 3 {
        return 0.0;
    }
    // Find extent in x and maximum height
    let (x_min, x_max) = particles
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(mn, mx), p| {
            (mn.min(p.position[0]), mx.max(p.position[0]))
        });
    let cx = (x_min + x_max) * 0.5;
    let mut max_slope: f64 = 0.0;
    for p in particles {
        let h = (p.position[1] - base_y).max(0.0);
        let r = (p.position[0] - cx).abs();
        if r > 1e-14 {
            let slope = h / r;
            if slope > max_slope {
                max_slope = slope;
            }
        }
    }
    max_slope.atan()
}

// ---------------------------------------------------------------------------
// Shear band detection
// ---------------------------------------------------------------------------

/// Estimate local vorticity (spin tensor) at a particle location.
///
/// Uses the angular velocity of the particle as a proxy for vorticity.
/// High vorticity concentrations indicate shear band formation.
///
/// Returns the vorticity magnitude ω = |angular_velocity|.
pub fn local_vorticity(particle: &GranularParticle) -> f64 {
    norm(particle.angular_velocity)
}

/// Detect shear band particles.
///
/// A particle is flagged as being in a shear band if its vorticity exceeds
/// `threshold` times the mean vorticity of the assembly.
///
/// Returns a vector of particle indices inside the shear band.
pub fn detect_shear_band(particles: &[GranularParticle], threshold_factor: f64) -> Vec<usize> {
    if particles.is_empty() {
        return Vec::new();
    }
    let mean_vort = particles.iter().map(local_vorticity).sum::<f64>() / particles.len() as f64;
    let limit = threshold_factor * mean_vort.max(1e-15);
    particles
        .iter()
        .enumerate()
        .filter(|(_, p)| local_vorticity(p) > limit)
        .map(|(i, _)| i)
        .collect()
}

// ---------------------------------------------------------------------------
// Force chain network
// ---------------------------------------------------------------------------

/// A link in the force chain network.
#[derive(Debug, Clone)]
pub struct ForceChainLink {
    /// Index of particle i.
    pub i: usize,
    /// Index of particle j.
    pub j: usize,
    /// Normal force magnitude \[N\].
    pub fn_mag: f64,
}

/// Extract the force chain network from a set of active contacts.
///
/// A force chain is a sequence of contacts carrying normal forces above
/// a given threshold (typically the mean contact force).
///
/// Returns a list of links above the threshold.
pub fn extract_force_chains(contacts: &[ContactPair], threshold: f64) -> Vec<ForceChainLink> {
    contacts
        .iter()
        .filter(|c| c.fn_mag >= threshold)
        .map(|c| ForceChainLink {
            i: c.i,
            j: c.j,
            fn_mag: c.fn_mag,
        })
        .collect()
}

/// Compute the mean normal force in a contact network.
pub fn mean_contact_force(contacts: &[ContactPair]) -> f64 {
    if contacts.is_empty() {
        return 0.0;
    }
    contacts.iter().map(|c| c.fn_mag).sum::<f64>() / contacts.len() as f64
}

// ---------------------------------------------------------------------------
// Coordination number and fabric tensor
// ---------------------------------------------------------------------------

/// Mean coordination number (contacts per particle).
pub fn coordination_number(n_particles: usize, n_contacts: usize) -> f64 {
    if n_particles == 0 {
        return 0.0;
    }
    2.0 * n_contacts as f64 / n_particles as f64
}

/// Compute the fabric tensor from the contact normal directions.
///
/// F_ij = (1/N_c) Σ n_i n_j
///
/// Returns a 3×3 symmetric matrix as `[[f64; 3]; 3]`.
pub fn fabric_tensor(particles: &[GranularParticle], contacts: &[ContactPair]) -> [[f64; 3]; 3] {
    let nc = contacts.len();
    if nc == 0 {
        return [[0.0; 3]; 3];
    }
    let mut f = [[0.0_f64; 3]; 3];
    for c in contacts {
        if c.i >= particles.len() || c.j >= particles.len() {
            continue;
        }
        let pi = &particles[c.i];
        let pj = &particles[c.j];
        let n = normalize(sub(pj.position, pi.position));
        for row in 0..3 {
            for col in 0..3 {
                f[row][col] += n[row] * n[col];
            }
        }
    }
    let s = 1.0 / nc as f64;
    for f_row in f.iter_mut() {
        for f_ij in f_row.iter_mut() {
            *f_ij *= s;
        }
    }
    f
}

// ---------------------------------------------------------------------------
// Segregation index
// ---------------------------------------------------------------------------

/// Compute a simple size-based segregation index.
///
/// The index is defined as the ratio of mean height of large particles
/// to mean height of small particles. Values > 1 indicate Brazil-nut effect.
///
/// # Arguments
/// - `particles`      – all particles.
/// - `size_threshold` – radius threshold separating "small" from "large".
pub fn segregation_index(particles: &[GranularParticle], size_threshold: f64) -> f64 {
    let small: Vec<f64> = particles
        .iter()
        .filter(|p| p.radius < size_threshold)
        .map(|p| p.position[1])
        .collect();
    let large: Vec<f64> = particles
        .iter()
        .filter(|p| p.radius >= size_threshold)
        .map(|p| p.position[1])
        .collect();
    if small.is_empty() || large.is_empty() {
        return 1.0;
    }
    let mean_small = small.iter().sum::<f64>() / small.len() as f64;
    let mean_large = large.iter().sum::<f64>() / large.len() as f64;
    mean_large / mean_small.max(1e-15)
}

// ---------------------------------------------------------------------------
// Janssen arching effect
// ---------------------------------------------------------------------------

/// Compute the vertical stress profile in a silo using the Janssen model.
///
/// σ_v(z) = (ρ g A / (k μ_w P)) (1 − e^(−k μ_w P z / A))
///
/// where A = cross-section area, P = perimeter, k = ratio of horizontal to
/// vertical stress (Rankine coefficient), μ_w = wall friction.
///
/// # Arguments
/// - `density`     – bulk density \[kg/m³\].
/// - `g`           – gravitational acceleration \[m/s²\].
/// - `radius`      – silo radius \[m\] (circular cross-section).
/// - `mu_wall`     – wall friction coefficient.
/// - `k_rankine`   – horizontal-to-vertical stress ratio (≈0.4 for sand).
/// - `z_values`    – depths at which to evaluate \[m\].
///
/// Returns a vector of vertical stresses \[Pa\] at each `z`.
pub fn janssen_stress_profile(
    density: f64,
    g: f64,
    radius: f64,
    mu_wall: f64,
    k_rankine: f64,
    z_values: &[f64],
) -> Vec<f64> {
    let area = PI * radius * radius;
    let perimeter = 2.0 * PI * radius;
    let lambda = k_rankine * mu_wall * perimeter / area;
    let sigma_inf = density * g / lambda;
    z_values
        .iter()
        .map(|&z| sigma_inf * (1.0 - (-lambda * z).exp()))
        .collect()
}

// ---------------------------------------------------------------------------
// Compaction simulation helper
// ---------------------------------------------------------------------------

/// State of a compaction cycle.
#[derive(Debug, Clone, Default)]
pub struct CompactionState {
    /// Current applied pressure \[Pa\].
    pub pressure: f64,
    /// Void ratio (V_voids / V_solid).
    pub void_ratio: f64,
    /// Number of load cycles completed.
    pub cycles: usize,
    /// Strain in the loading direction (engineering strain).
    pub strain: f64,
}

impl CompactionState {
    /// Create initial compaction state.
    pub fn new(initial_void_ratio: f64) -> Self {
        Self {
            pressure: 0.0,
            void_ratio: initial_void_ratio,
            cycles: 0,
            strain: 0.0,
        }
    }

    /// Apply a pressure increment following a semi-log compressibility model.
    ///
    /// Δe = −C_c Δlog₁₀(σ)
    ///
    /// # Arguments
    /// - `new_pressure` – new applied pressure \[Pa\].
    /// - `cc`           – compression index (slope of e–log σ line).
    pub fn apply_pressure(&mut self, new_pressure: f64, cc: f64) {
        if new_pressure <= self.pressure.max(1e-6) {
            return;
        }
        let delta_log = (new_pressure / self.pressure.max(1e-6)).log10();
        self.void_ratio -= cc * delta_log;
        self.void_ratio = self.void_ratio.max(0.01);
        let old_strain = self.strain;
        self.strain = 1.0 - 1.0 / (1.0 + self.void_ratio).max(1e-6);
        let _ = old_strain; // strain change is implicit
        self.pressure = new_pressure;
    }

    /// Increment cycle count (one load-unload cycle).
    pub fn increment_cycle(&mut self) {
        self.cycles += 1;
    }

    /// Packing fraction derived from void ratio: φ = 1 / (1 + e).
    pub fn packing_fraction(&self) -> f64 {
        1.0 / (1.0 + self.void_ratio).max(1e-6)
    }
}

// ---------------------------------------------------------------------------
// Avalanche dynamics
// ---------------------------------------------------------------------------

/// Check if the current slope angle exceeds the dynamic angle of repose.
///
/// Returns `true` (avalanche triggered) when θ > θ_dynamic.
///
/// # Arguments
/// - `slope_angle`     – current heap surface slope \[rad\].
/// - `dynamic_repose`  – dynamic angle of repose \[rad\] (typically ~5° < static).
pub fn avalanche_triggered(slope_angle: f64, dynamic_repose: f64) -> bool {
    slope_angle > dynamic_repose
}

/// Estimate the runout distance of an avalanche.
///
/// Uses the Heim (1932) fahrböschung (H/L) ratio:
/// L_runout = H / tan(φ_mob)  where φ_mob is the mobilized friction angle.
///
/// # Arguments
/// - `fall_height`   – total drop height H \[m\].
/// - `mob_friction`  – mobilised friction angle φ_mob \[rad\].
pub fn avalanche_runout(fall_height: f64, mob_friction: f64) -> f64 {
    fall_height / mob_friction.tan().max(1e-14)
}

// ---------------------------------------------------------------------------
// GranularSystem — full simulation
// ---------------------------------------------------------------------------

/// A complete granular simulation system.
#[derive(Debug, Clone)]
pub struct GranularSystem {
    /// All particles in the assembly.
    pub particles: Vec<GranularParticle>,
    /// Active contact pairs with persistent tangential spring state.
    pub contacts: Vec<ContactPair>,
    /// Gravitational acceleration \[m/s²\].
    pub gravity: [f64; 3],
    /// Domain bounding box half-extents \[m\] (for wall contacts).
    pub domain: [[f64; 2]; 3],
    /// Time step \[s\].
    pub dt: f64,
    /// Optional capillary bridge parameters: (surface_tension, bridge_volume).
    pub capillary_params: Option<(f64, f64)>,
    /// Optional Hamaker constant for van der Waals.
    pub hamaker: Option<f64>,
    /// Tensile strength for crushing criterion \[Pa\].
    pub tensile_strength: f64,
    /// Accumulated simulation time \[s\].
    pub time: f64,
}

impl GranularSystem {
    /// Create a new granular system with the given parameters.
    pub fn new(gravity: [f64; 3], domain: [[f64; 2]; 3], dt: f64) -> Self {
        Self {
            particles: Vec::new(),
            contacts: Vec::new(),
            gravity,
            domain,
            dt,
            capillary_params: None,
            hamaker: None,
            tensile_strength: 1e7,
            time: 0.0,
        }
    }

    /// Add a particle and return its index.
    pub fn add_particle(&mut self, mut p: GranularParticle) -> usize {
        let idx = self.particles.len();
        p.id = idx;
        self.particles.push(p);
        idx
    }

    /// Detect pairwise contacts (O(N²) brute force).
    pub fn detect_contacts(&mut self) {
        self.contacts.clear();
        let n = self.particles.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = &self.particles[i];
                let pj = &self.particles[j];
                let d = norm(sub(pj.position, pi.position));
                let overlap = pi.radius + pj.radius - d;
                if overlap > 0.0 {
                    let mut cp = ContactPair::new(i, j);
                    cp.overlap = overlap;
                    self.contacts.push(cp);
                }
            }
        }
    }

    /// Apply gravity body force to all particles.
    pub fn apply_gravity(&mut self) {
        for p in &mut self.particles {
            p.force = add(p.force, scale(self.gravity, p.mass));
        }
    }

    /// Apply wall contact forces (axis-aligned domain walls).
    pub fn apply_wall_contacts(&mut self) {
        for p in &mut self.particles {
            for axis in 0..3 {
                let lo = self.domain[axis][0];
                let hi = self.domain[axis][1];
                // Low wall
                let pen_lo = p.radius - (p.position[axis] - lo);
                if pen_lo > 0.0 {
                    let r_star = p.radius;
                    let e_star = p.young / (2.0 * (1.0 - p.poisson * p.poisson));
                    let fn_mag = hertz_normal_force(pen_lo, e_star, r_star);
                    p.force[axis] += fn_mag;
                }
                // High wall
                let pen_hi = p.radius - (hi - p.position[axis]);
                if pen_hi > 0.0 {
                    let r_star = p.radius;
                    let e_star = p.young / (2.0 * (1.0 - p.poisson * p.poisson));
                    let fn_mag = hertz_normal_force(pen_hi, e_star, r_star);
                    p.force[axis] -= fn_mag;
                }
            }
        }
    }

    /// Resolve all pairwise contacts and accumulate forces.
    pub fn resolve_contacts(&mut self) {
        // We need temporary contact spring storage
        let cap = self.capillary_params;
        let hk = self.hamaker;
        let dt = self.dt;

        let mut updated_springs: Vec<[f64; 3]> =
            self.contacts.iter().map(|c| c.tangential_spring).collect();
        let mut updated_fn: Vec<f64> = vec![0.0; self.contacts.len()];

        for (ci, c) in self.contacts.iter().enumerate() {
            let pi = &self.particles[c.i];
            let pj = &self.particles[c.j];
            if let Some(res) = resolve_hertz_mindlin(pi, pj, c.tangential_spring, dt, cap, hk) {
                updated_springs[ci] = res.new_tangential_spring;
                updated_fn[ci] = res.fn_mag;

                // accumulate forces/torques (borrow checker: use indices)
                self.particles[c.i].force = add(self.particles[c.i].force, res.force_i);
                self.particles[c.j].force = add(self.particles[c.j].force, res.force_j);
                self.particles[c.i].torque = add(self.particles[c.i].torque, res.torque_i);
                self.particles[c.j].torque = add(self.particles[c.j].torque, res.torque_j);
                self.particles[c.i].strain_energy += res.strain_energy;
                self.particles[c.j].strain_energy += res.strain_energy;
            }
        }

        // Write back spring and fn_mag
        for (ci, c) in self.contacts.iter_mut().enumerate() {
            c.tangential_spring = updated_springs[ci];
            c.fn_mag = updated_fn[ci];
        }
    }

    /// Integrate all particles (semi-implicit Euler).
    pub fn integrate(&mut self) {
        let dt = self.dt;
        for p in &mut self.particles {
            if p.crushed {
                continue;
            }
            // Linear
            let a = scale(p.force, 1.0 / p.mass);
            p.velocity = add(p.velocity, scale(a, dt));
            p.position = add(p.position, scale(p.velocity, dt));
            // Angular
            let alpha = scale(p.torque, 1.0 / p.inertia);
            p.angular_velocity = add(p.angular_velocity, scale(alpha, dt));
            p.reset_accumulators();
        }
        self.time += dt;
    }

    /// Check and apply grain crushing.
    pub fn apply_crushing(&mut self) {
        let ts = self.tensile_strength;
        let to_crush: Vec<usize> = self
            .particles
            .iter()
            .enumerate()
            .filter(|(_, p)| {
                !p.crushed && crushing_criterion(p.strain_energy, ts, p.young, p.radius)
            })
            .map(|(i, _)| i)
            .collect();
        for &idx in &to_crush {
            self.particles[idx].crushed = true;
        }
    }

    /// Run one full simulation step.
    pub fn step(&mut self) {
        for p in &mut self.particles {
            p.reset_accumulators();
        }
        self.apply_gravity();
        self.detect_contacts();
        self.resolve_contacts();
        self.apply_wall_contacts();
        self.integrate();
        self.apply_crushing();
    }

    /// Compute packing fraction from current particle positions within domain.
    pub fn packing_fraction(&self) -> f64 {
        let vol: f64 = self
            .particles
            .iter()
            .map(|p| (4.0 / 3.0) * PI * p.radius.powi(3))
            .sum();
        let dv: f64 = (0..3)
            .map(|k| self.domain[k][1] - self.domain[k][0])
            .product();
        if dv < 1e-30 {
            return 0.0;
        }
        (vol / dv).min(1.0)
    }

    /// Return the granular temperature.
    pub fn granular_temperature(&self) -> f64 {
        granular_temperature(&self.particles)
    }

    /// Return force chain links above mean force.
    pub fn force_chains(&self) -> Vec<ForceChainLink> {
        let mean = mean_contact_force(&self.contacts);
        extract_force_chains(&self.contacts, mean)
    }

    /// Return the fabric tensor.
    pub fn fabric_tensor(&self) -> [[f64; 3]; 3] {
        fabric_tensor(&self.particles, &self.contacts)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(id: usize, radius: f64) -> GranularParticle {
        GranularParticle::new(id, radius, 2500.0)
    }

    // -- hertz_normal_force ---------------------------------------------------

    #[test]
    fn hertz_normal_force_zero_overlap() {
        assert_eq!(hertz_normal_force(0.0, 1e9, 0.01), 0.0);
    }

    #[test]
    fn hertz_normal_force_negative_overlap() {
        assert_eq!(hertz_normal_force(-0.001, 1e9, 0.01), 0.0);
    }

    #[test]
    fn hertz_normal_force_positive() {
        let f = hertz_normal_force(1e-4, 1e9, 0.01);
        assert!(f > 0.0, "f={f}");
    }

    #[test]
    fn hertz_normal_force_monotone() {
        let f1 = hertz_normal_force(1e-4, 1e9, 0.01);
        let f2 = hertz_normal_force(2e-4, 1e9, 0.01);
        assert!(f2 > f1);
    }

    #[test]
    fn hertz_normal_force_formula() {
        let e = 1e9_f64;
        let r = 0.01_f64;
        let d = 1e-4_f64;
        let expected = (4.0 / 3.0) * e * r.sqrt() * d.powf(1.5);
        let result = hertz_normal_force(d, e, r);
        assert!(
            (result - expected).abs() < 1e-6,
            "result={result} expected={expected}"
        );
    }

    // -- mindlin_tangential_force ---------------------------------------------

    #[test]
    fn mindlin_tangential_force_zero_displacement() {
        assert_eq!(mindlin_tangential_force(0.0, 1e8, 0.01, 100.0, 0.4), 0.0);
    }

    #[test]
    fn mindlin_tangential_force_capped_by_coulomb() {
        let ft = mindlin_tangential_force(100.0, 1e8, 0.01, 100.0, 0.4);
        assert!((ft - 40.0).abs() < 5.0, "ft={ft}");
    }

    // -- rolling_resistance_torque -------------------------------------------

    #[test]
    fn rolling_resistance_torque_zero_force() {
        assert_eq!(rolling_resistance_torque(0.01, 0.05, 0.0, 1.0), 0.0);
    }

    #[test]
    fn rolling_resistance_torque_formula() {
        let t = rolling_resistance_torque(0.02, 0.03, 50.0, 1.0);
        assert!((t - 0.02 * 0.03 * 50.0).abs() < 1e-12);
    }

    #[test]
    fn rolling_resistance_torque_shape_multiplier() {
        let t1 = rolling_resistance_torque(0.01, 0.05, 100.0, 1.0);
        let t2 = rolling_resistance_torque(0.01, 0.05, 100.0, 2.0);
        assert!((t2 - 2.0 * t1).abs() < 1e-12);
    }

    // -- capillary_force ------------------------------------------------------

    #[test]
    fn capillary_force_zero_volume() {
        assert_eq!(capillary_force(0.072, 0.01, 0.0, 0.0), 0.0);
    }

    #[test]
    fn capillary_force_decreases_with_separation() {
        let f0 = capillary_force(0.072, 0.01, 0.0, 1e-9);
        let f1 = capillary_force(0.072, 0.01, 1e-6, 1e-9);
        assert!(f1 < f0, "f1={f1} f0={f0}");
    }

    // -- van_der_waals_force --------------------------------------------------

    #[test]
    fn van_der_waals_force_increases_as_gap_decreases() {
        let f1 = van_der_waals_force(1e-20, 1e-6, 1e-9, 1e-11);
        let f2 = van_der_waals_force(1e-20, 1e-6, 5e-10, 1e-11);
        assert!(f2 > f1, "f2={f2} f1={f1}");
    }

    // -- crushing_criterion ---------------------------------------------------

    #[test]
    fn crushing_criterion_no_crush_at_zero_energy() {
        assert!(!crushing_criterion(0.0, 1e6, 1e9, 0.01));
    }

    #[test]
    fn crushing_criterion_triggers_at_high_energy() {
        // Large strain energy → should trigger
        assert!(crushing_criterion(1e6, 1e6, 1e9, 0.01));
    }

    // -- split_grain ----------------------------------------------------------

    #[test]
    fn split_grain_conserves_mass_approximately() {
        let parent = make_particle(0, 0.1);
        let dir = [1.0, 0.0, 0.0];
        let (da, db) = split_grain(&parent, dir);
        let total_mass = da.mass + db.mass;
        assert!(
            (total_mass - parent.mass).abs() / parent.mass < 0.01,
            "mass conservation: {} vs {}",
            total_mass,
            parent.mass
        );
    }

    #[test]
    fn split_grain_daughters_smaller() {
        let parent = make_particle(0, 0.1);
        let (da, _db) = split_grain(&parent, [0.0, 1.0, 0.0]);
        assert!(da.radius < parent.radius);
    }

    // -- fcc_packing ----------------------------------------------------------

    #[test]
    fn fcc_packing_non_empty() {
        let positions = fcc_packing(0.05, 1.0, 1.0, 1.0);
        assert!(!positions.is_empty(), "FCC packing should have particles");
    }

    #[test]
    fn fcc_packing_no_overlap() {
        let r = 0.1;
        let positions = fcc_packing(r, 1.0, 1.0, 1.0);
        for i in 0..positions.len() {
            for j in (i + 1)..positions.len() {
                let d = norm(sub(positions[j], positions[i]));
                assert!(d > 1.99 * r, "overlap at ({i},{j}): d={d}");
            }
        }
    }

    // -- hcp_packing ----------------------------------------------------------

    #[test]
    fn hcp_packing_non_empty() {
        let positions = hcp_packing(0.05, 1.0, 1.0, 1.0);
        assert!(!positions.is_empty());
    }

    // -- packing_fraction -----------------------------------------------------

    #[test]
    fn packing_fraction_empty() {
        assert_eq!(packing_fraction(&[], 1.0), 0.0);
    }

    #[test]
    fn packing_fraction_single_sphere_in_large_domain() {
        let phi = packing_fraction(&[0.01], 1.0);
        assert!(phi > 0.0 && phi < 0.01);
    }

    // -- granular_temperature -------------------------------------------------

    #[test]
    fn granular_temperature_zero_for_uniform_motion() {
        // All particles moving at the same velocity → zero fluctuation
        let mut particles = vec![make_particle(0, 0.05), make_particle(1, 0.05)];
        particles[0].velocity = [1.0, 0.0, 0.0];
        particles[1].velocity = [1.0, 0.0, 0.0];
        let t = granular_temperature(&particles);
        assert!(t.abs() < 1e-20, "T={t}");
    }

    #[test]
    fn granular_temperature_positive_for_fluctuations() {
        let mut particles = vec![make_particle(0, 0.05), make_particle(1, 0.05)];
        particles[0].velocity = [1.0, 0.0, 0.0];
        particles[1].velocity = [-1.0, 0.0, 0.0];
        let t = granular_temperature(&particles);
        assert!(t > 0.0, "T={t}");
    }

    #[test]
    fn granular_temperature_empty_assembly() {
        assert_eq!(granular_temperature(&[]), 0.0);
    }

    // -- angle_of_repose ------------------------------------------------------

    #[test]
    fn angle_of_repose_estimate_positive() {
        let theta = angle_of_repose_estimate(0.4, 0.01, 1.0);
        assert!(theta > 0.0);
    }

    #[test]
    fn angle_of_repose_estimate_sphere_lower_than_angular() {
        let theta_sphere = angle_of_repose_estimate(0.4, 0.01, 1.0);
        let theta_angular = angle_of_repose_estimate(0.4, 0.01, 2.0);
        assert!(theta_angular > theta_sphere);
    }

    // -- shear band detection -------------------------------------------------

    #[test]
    fn detect_shear_band_empty() {
        assert!(detect_shear_band(&[], 2.0).is_empty());
    }

    #[test]
    fn detect_shear_band_finds_high_vorticity() {
        let mut particles = vec![make_particle(0, 0.05), make_particle(1, 0.05)];
        particles[0].angular_velocity = [0.0, 0.0, 100.0]; // high vorticity
        particles[1].angular_velocity = [0.0, 0.0, 0.0]; // low vorticity
        let band = detect_shear_band(&particles, 1.5);
        assert!(band.contains(&0));
    }

    // -- force chain ----------------------------------------------------------

    #[test]
    fn extract_force_chains_above_threshold() {
        let mut c1 = ContactPair::new(0, 1);
        c1.fn_mag = 100.0;
        let mut c2 = ContactPair::new(1, 2);
        c2.fn_mag = 10.0;
        let chains = extract_force_chains(&[c1, c2], 50.0);
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].i, 0);
    }

    #[test]
    fn mean_contact_force_correct() {
        let mut c1 = ContactPair::new(0, 1);
        c1.fn_mag = 100.0;
        let mut c2 = ContactPair::new(1, 2);
        c2.fn_mag = 200.0;
        let mean = mean_contact_force(&[c1, c2]);
        assert!((mean - 150.0).abs() < 1e-10);
    }

    // -- coordination number --------------------------------------------------

    #[test]
    fn coordination_number_one_contact_two_particles() {
        let cn = coordination_number(2, 1);
        assert!((cn - 1.0).abs() < 1e-12);
    }

    #[test]
    fn coordination_number_zero_particles() {
        assert_eq!(coordination_number(0, 0), 0.0);
    }

    // -- fabric tensor --------------------------------------------------------

    #[test]
    fn fabric_tensor_empty_contacts() {
        let f = fabric_tensor(&[], &[]);
        for row in f.iter() {
            for &v in row.iter() {
                assert_eq!(v, 0.0);
            }
        }
    }

    #[test]
    fn fabric_tensor_trace_is_one() {
        let mut pi = make_particle(0, 0.05);
        let mut pj = make_particle(1, 0.05);
        pi.position = [0.0; 3];
        pj.position = [0.08, 0.0, 0.0];
        let mut c = ContactPair::new(0, 1);
        c.fn_mag = 100.0;
        let f = fabric_tensor(&[pi, pj], &[c]);
        let trace = f[0][0] + f[1][1] + f[2][2];
        assert!((trace - 1.0).abs() < 1e-10, "trace={trace}");
    }

    // -- segregation index ----------------------------------------------------

    #[test]
    fn segregation_index_no_particles_returns_one() {
        assert_eq!(segregation_index(&[], 0.05), 1.0);
    }

    #[test]
    fn segregation_index_brazil_nut() {
        // Large particles higher than small ones → index > 1
        let mut small = make_particle(0, 0.02);
        small.position = [0.0, 0.1, 0.0];
        let mut large = make_particle(1, 0.08);
        large.position = [0.0, 1.0, 0.0];
        let idx = segregation_index(&[small, large], 0.05);
        assert!(idx > 1.0, "idx={idx}");
    }

    // -- Janssen model --------------------------------------------------------

    #[test]
    fn janssen_stress_increases_with_depth() {
        let z = vec![0.0, 0.5, 1.0, 2.0, 5.0];
        let sigma = janssen_stress_profile(1500.0, 9.81, 0.5, 0.4, 0.4, &z);
        for i in 0..sigma.len() - 1 {
            assert!(
                sigma[i + 1] >= sigma[i],
                "stress should increase: {}",
                sigma[i]
            );
        }
    }

    #[test]
    fn janssen_stress_zero_at_surface() {
        let z = vec![0.0];
        let sigma = janssen_stress_profile(1500.0, 9.81, 0.5, 0.4, 0.4, &z);
        assert!(sigma[0].abs() < 1e-10);
    }

    #[test]
    fn janssen_stress_saturates_at_large_depth() {
        let z = vec![100.0, 200.0, 1000.0];
        let sigma = janssen_stress_profile(1500.0, 9.81, 0.5, 0.4, 0.4, &z);
        // Large depth → σ approaches σ_inf → difference should be small
        let diff = (sigma[2] - sigma[1]).abs();
        assert!(diff < sigma[2] * 0.01, "should saturate: diff={diff}");
    }

    // -- CompactionState ------------------------------------------------------

    #[test]
    fn compaction_state_initial_packing() {
        let state = CompactionState::new(0.6);
        let phi = state.packing_fraction();
        assert!((phi - 1.0 / 1.6).abs() < 1e-12);
    }

    #[test]
    fn compaction_apply_pressure_reduces_void_ratio() {
        let mut state = CompactionState::new(0.8);
        state.pressure = 100.0;
        state.apply_pressure(1000.0, 0.3);
        assert!(state.void_ratio < 0.8, "void ratio should decrease");
    }

    #[test]
    fn compaction_cycle_counter() {
        let mut state = CompactionState::new(0.6);
        state.increment_cycle();
        state.increment_cycle();
        assert_eq!(state.cycles, 2);
    }

    // -- Avalanche dynamics ---------------------------------------------------

    #[test]
    fn avalanche_triggered_above_repose() {
        let theta_repose = 0.5_f64; // ~29°
        assert!(avalanche_triggered(0.6, theta_repose));
        assert!(!avalanche_triggered(0.4, theta_repose));
    }

    #[test]
    fn avalanche_runout_positive() {
        let l = avalanche_runout(10.0, 0.4);
        assert!(l > 0.0, "l={l}");
    }

    // -- GranularSystem -------------------------------------------------------

    #[test]
    fn system_add_and_count_particles() {
        let mut sys = GranularSystem::new(
            [0.0, -9.81, 0.0],
            [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]],
            1e-4,
        );
        sys.add_particle(make_particle(0, 0.05));
        sys.add_particle(make_particle(1, 0.05));
        assert_eq!(sys.particles.len(), 2);
    }

    #[test]
    fn system_step_does_not_panic() {
        let mut sys = GranularSystem::new(
            [0.0, -9.81, 0.0],
            [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
            1e-4,
        );
        let mut p = make_particle(0, 0.05);
        p.position = [0.0, 0.5, 0.0];
        sys.add_particle(p);
        sys.step();
    }

    #[test]
    fn system_contact_detection_overlap() {
        let mut sys = GranularSystem::new(
            [0.0, -9.81, 0.0],
            [[-5.0, 5.0], [-5.0, 5.0], [-5.0, 5.0]],
            1e-4,
        );
        let mut p0 = make_particle(0, 0.05);
        let mut p1 = make_particle(1, 0.05);
        p0.position = [0.0; 3];
        p1.position = [0.08, 0.0, 0.0]; // overlap
        sys.add_particle(p0);
        sys.add_particle(p1);
        sys.detect_contacts();
        assert_eq!(sys.contacts.len(), 1);
    }

    #[test]
    fn system_packing_fraction_positive() {
        let mut sys = GranularSystem::new(
            [0.0, -9.81, 0.0],
            [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]],
            1e-4,
        );
        sys.add_particle(make_particle(0, 0.1));
        assert!(sys.packing_fraction() > 0.0);
    }

    #[test]
    fn system_granular_temperature_zero_at_rest() {
        let mut sys = GranularSystem::new(
            [0.0, -9.81, 0.0],
            [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]],
            1e-4,
        );
        sys.add_particle(make_particle(0, 0.05));
        sys.add_particle(make_particle(1, 0.05));
        assert!((sys.granular_temperature()).abs() < 1e-20);
    }

    // -- GrainShape -----------------------------------------------------------

    #[test]
    fn grain_shape_sphere_factor_one() {
        let s = GrainShape::sphere();
        assert!((s.factor - 1.0).abs() < 1e-12);
    }

    #[test]
    fn grain_shape_angular_rolling_resistance_greater() {
        let sphere = GrainShape::sphere();
        let angular = GrainShape::angular(2.0, 1.5, 0.5);
        assert!(angular.rolling_resistance_multiplier() > sphere.rolling_resistance_multiplier());
    }

    // -- measure_angle_of_repose ----------------------------------------------

    #[test]
    fn measure_angle_of_repose_returns_nonnegative() {
        let mut particles = vec![
            make_particle(0, 0.05),
            make_particle(1, 0.05),
            make_particle(2, 0.05),
        ];
        particles[0].position = [0.0, 0.1, 0.0];
        particles[1].position = [0.5, 0.3, 0.0];
        particles[2].position = [1.0, 0.1, 0.0];
        let theta = measure_angle_of_repose(&particles, 0.0);
        assert!(theta >= 0.0);
    }
}
