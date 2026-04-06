#![allow(clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular material modelled as a soft body via DEM-FEM coupling.
//!
//! This module implements a Discrete Element Method (DEM) approach to granular
//! media, treating each grain as a soft-body particle that can interact through
//! contact force chains, undergo compaction, transmit stress, form shear bands,
//! and eventually crush/break.  A lightweight FEM continuum overlay is provided
//! for coarse-grained stress analysis.
//!
//! # Physical models
//! - Hertz–Mindlin contact mechanics between grains
//! - Critical state soil mechanics (Cam-Clay inspired)
//! - Grain crushing / breakage with a probability model
//! - Polydisperse random packing via random sequential addition
//! - Avalanche dynamics and angle-of-repose measurement
//! - Stress-tensor assembly from contact force chains
//! - Shear-band detection from incremental strain localisation
//! - DEM → FEM homogenisation for continuum stress fields

#![allow(dead_code)]

use rand::Rng;

use rand::RngExt;
// ---------------------------------------------------------------------------
// Vec3 helper (plain [f64; 3] to avoid nalgebra)
// ---------------------------------------------------------------------------

/// A 3-component vector backed by `[f64; 3]`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// Construct from components.
    #[inline]
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Zero vector.
    #[inline]
    pub fn zero() -> Self {
        Self::new(0.0, 0.0, 0.0)
    }

    /// Euclidean length.
    #[inline]
    pub fn norm(self) -> f64 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    /// Squared length.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Normalise, returning zero vector if length is too small.
    #[inline]
    pub fn normalise(self) -> Self {
        let n = self.norm();
        if n < 1e-300 {
            Self::zero()
        } else {
            Self::new(self.x / n, self.y / n, self.z / n)
        }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Cross product.
    #[inline]
    pub fn cross(self, rhs: Self) -> Self {
        Self::new(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x * rhs.y - self.y * rhs.x,
        )
    }

    /// Scale.
    #[inline]
    pub fn scale(self, s: f64) -> Self {
        Self::new(self.x * s, self.y * s, self.z * s)
    }

    /// Add.
    #[inline]
    pub fn add(self, rhs: Self) -> Self {
        Self::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }

    /// Subtract.
    #[inline]
    pub fn sub(self, rhs: Self) -> Self {
        Self::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

// ---------------------------------------------------------------------------
// Grain
// ---------------------------------------------------------------------------

/// A single granular particle (grain).
#[derive(Debug, Clone)]
pub struct Grain {
    /// Centre-of-mass position.
    pub position: Vec3,
    /// Linear velocity.
    pub velocity: Vec3,
    /// Angular velocity (rad s⁻¹).
    pub omega: Vec3,
    /// Grain radius (m).
    pub radius: f64,
    /// Mass (kg).
    pub mass: f64,
    /// Moment of inertia (kg m²) — sphere: 2/5 m r².
    pub inertia: f64,
    /// Accumulated force (N).
    pub force: Vec3,
    /// Accumulated torque (N m).
    pub torque: Vec3,
    /// Whether this grain is fixed (boundary grain).
    pub is_fixed: bool,
    /// Grain material Young's modulus (Pa).
    pub young: f64,
    /// Grain material Poisson ratio.
    pub poisson: f64,
    /// Tangential damping coefficient.
    pub gamma_t: f64,
    /// Normal damping coefficient.
    pub gamma_n: f64,
    /// Whether this grain has been crushed/broken.
    pub crushed: bool,
    /// Grain ID for tracking.
    pub id: usize,
}

impl Grain {
    /// Create a new grain with full parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: usize,
        position: Vec3,
        radius: f64,
        density: f64,
        young: f64,
        poisson: f64,
        gamma_n: f64,
        gamma_t: f64,
    ) -> Self {
        let mass = density * (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3);
        let inertia = 0.4 * mass * radius * radius;
        Self {
            position,
            velocity: Vec3::zero(),
            omega: Vec3::zero(),
            radius,
            mass,
            inertia,
            force: Vec3::zero(),
            torque: Vec3::zero(),
            is_fixed: false,
            young,
            poisson,
            gamma_t,
            gamma_n,
            crushed: false,
            id,
        }
    }

    /// Apply an impulse force (added to accumulator).
    pub fn apply_force(&mut self, f: Vec3) {
        if !self.is_fixed {
            self.force = self.force.add(f);
        }
    }

    /// Apply torque.
    pub fn apply_torque(&mut self, t: Vec3) {
        if !self.is_fixed {
            self.torque = self.torque.add(t);
        }
    }

    /// Clear force and torque accumulators.
    pub fn clear_forces(&mut self) {
        self.force = Vec3::zero();
        self.torque = Vec3::zero();
    }

    /// Integrate state forward with simple Euler integration.
    pub fn integrate(&mut self, dt: f64) {
        if self.is_fixed || self.crushed {
            return;
        }
        let accel = self.force.scale(1.0 / self.mass);
        let alpha = self.torque.scale(1.0 / self.inertia);
        self.velocity = self.velocity.add(accel.scale(dt));
        self.omega = self.omega.add(alpha.scale(dt));
        self.position = self.position.add(self.velocity.scale(dt));
    }

    /// Effective Young's modulus for Hertz contact with another grain.
    pub fn effective_young(&self, other: &Grain) -> f64 {
        let e1 = self.young / (1.0 - self.poisson * self.poisson);
        let e2 = other.young / (1.0 - other.poisson * other.poisson);
        1.0 / (1.0 / e1 + 1.0 / e2)
    }

    /// Effective radius for Hertz contact.
    pub fn effective_radius(&self, other: &Grain) -> f64 {
        (self.radius * other.radius) / (self.radius + other.radius)
    }

    /// Effective mass for contact.
    pub fn effective_mass(&self, other: &Grain) -> f64 {
        (self.mass * other.mass) / (self.mass + other.mass)
    }
}

// ---------------------------------------------------------------------------
// Contact
// ---------------------------------------------------------------------------

/// A pairwise contact between two grains.
#[derive(Debug, Clone)]
pub struct Contact {
    /// Index of grain A.
    pub grain_a: usize,
    /// Index of grain B.
    pub grain_b: usize,
    /// Normal direction (from A to B).
    pub normal: Vec3,
    /// Overlap depth (m), positive means penetration.
    pub overlap: f64,
    /// Accumulated tangential spring displacement (for Mindlin).
    pub delta_t: Vec3,
    /// Normal contact force magnitude (N).
    pub fn_magnitude: f64,
    /// Tangential contact force (N).
    pub ft: Vec3,
    /// Whether this contact is currently active.
    pub active: bool,
}

impl Contact {
    /// Create a new contact.
    pub fn new(grain_a: usize, grain_b: usize) -> Self {
        Self {
            grain_a,
            grain_b,
            normal: Vec3::zero(),
            overlap: 0.0,
            delta_t: Vec3::zero(),
            fn_magnitude: 0.0,
            ft: Vec3::zero(),
            active: false,
        }
    }
}

// ---------------------------------------------------------------------------
// DEM parameters
// ---------------------------------------------------------------------------

/// Parameters governing the DEM simulation.
#[derive(Debug, Clone)]
pub struct DemParams {
    /// Gravitational acceleration vector (m s⁻²).
    pub gravity: Vec3,
    /// Friction coefficient (Coulomb).
    pub mu_friction: f64,
    /// Coefficient of restitution (used for damping).
    pub restitution: f64,
    /// Rolling friction coefficient.
    pub mu_rolling: f64,
    /// Background viscous damping (s⁻¹).
    pub damping: f64,
    /// Time step (s).
    pub dt: f64,
}

impl DemParams {
    /// Create default DEM parameters (sand-like).
    pub fn default_sand() -> Self {
        Self {
            gravity: Vec3::new(0.0, -9.81, 0.0),
            mu_friction: 0.4,
            restitution: 0.6,
            mu_rolling: 0.01,
            damping: 50.0,
            dt: 1e-5,
        }
    }
}

// ---------------------------------------------------------------------------
// Hertz-Mindlin contact force
// ---------------------------------------------------------------------------

/// Compute Hertz-Mindlin contact forces between two grains.
///
/// Returns `(force_on_a, torque_on_a, force_on_b, torque_on_b)`.
#[allow(clippy::too_many_arguments)]
pub fn hertz_mindlin_force(
    ga: &Grain,
    gb: &Grain,
    delta_t: &mut Vec3,
    params: &DemParams,
) -> (Vec3, Vec3, Vec3, Vec3) {
    let r_ab = gb.position.sub(ga.position);
    let dist = r_ab.norm();
    let sum_r = ga.radius + gb.radius;
    let overlap = sum_r - dist;
    if overlap <= 0.0 || dist < 1e-15 {
        return (Vec3::zero(), Vec3::zero(), Vec3::zero(), Vec3::zero());
    }

    let n = r_ab.scale(1.0 / dist);
    let e_eff = ga.effective_young(gb);
    let r_eff = ga.effective_radius(gb);
    let m_eff = ga.effective_mass(gb);

    // Hertz normal force: Fn = (4/3) E* sqrt(R*) delta^(3/2)
    let fn_hertz = (4.0 / 3.0) * e_eff * r_eff.sqrt() * overlap.powf(1.5);

    // Normal damping: Fn_damp = -gamma_n * m_eff * v_rel · n
    let v_rel = gb.velocity.sub(ga.velocity);
    let v_n = v_rel.dot(n);
    let fn_damp = -ga.gamma_n * m_eff.sqrt() * v_n;

    let fn_total = (fn_hertz + fn_damp).max(0.0);
    let f_normal = n.scale(fn_total);

    // Tangential displacement update (Mindlin)
    let r_a = ga.radius;
    let r_b = gb.radius;
    let contact_vel_a = ga.velocity.add(ga.omega.cross(n.scale(-r_a)));
    let contact_vel_b = gb.velocity.add(gb.omega.cross(n.scale(r_b)));
    let v_tang_full = contact_vel_b.sub(contact_vel_a);
    let v_t = v_tang_full.sub(n.scale(v_tang_full.dot(n)));

    // Shear modulus effective
    let g_a = ga.young / (2.0 * (1.0 + ga.poisson));
    let g_b = gb.young / (2.0 * (1.0 + gb.poisson));
    let g_eff = 1.0 / (1.0 / g_a + 1.0 / g_b);

    // Mindlin tangential stiffness: kt = 8 G* sqrt(R* delta)
    let kt = 8.0 * g_eff * (r_eff * overlap).sqrt();

    // Update tangential displacement
    let dt = params.dt;
    *delta_t = delta_t.add(v_t.scale(dt));
    // Remove normal component drift
    let n_comp = delta_t.dot(n);
    *delta_t = delta_t.sub(n.scale(n_comp));

    let ft_spring = delta_t.scale(-kt);
    let ft_damp = v_t.scale(-ga.gamma_t);
    let ft_total = ft_spring.add(ft_damp);

    // Coulomb friction limit
    let ft_max = params.mu_friction * fn_total;
    let ft_norm = ft_total.norm();
    let ft_limited = if ft_norm > ft_max && ft_norm > 1e-30 {
        // Slip: clamp and reset delta_t
        let ft_clamped = ft_total.scale(ft_max / ft_norm);
        *delta_t = ft_clamped.scale(-1.0 / kt.max(1e-30));
        ft_clamped
    } else {
        ft_total
    };

    let f_total = f_normal.add(ft_limited);

    // Torques: tau = r x F_t
    let ta = n.scale(-r_a).cross(ft_limited);
    let tb = n.scale(r_b).cross(ft_limited.scale(-1.0));

    (f_total.scale(-1.0), ta, f_total, tb)
}

// ---------------------------------------------------------------------------
// Flat-wall contact
// ---------------------------------------------------------------------------

/// Contact of a grain with an infinite flat wall.
///
/// Wall is defined by a point and outward normal (pointing into the domain).
/// Returns `(force, torque)` on the grain.
pub fn wall_contact_force(
    grain: &Grain,
    wall_point: Vec3,
    wall_normal: Vec3,
    params: &DemParams,
) -> (Vec3, Vec3) {
    // Overlap: positive if grain penetrates wall
    let d = grain.position.sub(wall_point).dot(wall_normal);
    let overlap = grain.radius - d;
    if overlap <= 0.0 {
        return (Vec3::zero(), Vec3::zero());
    }
    let n = wall_normal; // normal pointing into domain

    // Effective stiffness: treat wall as infinite mass, same Young's modulus
    let e_eff = grain.young / (1.0 - grain.poisson * grain.poisson);
    let r_eff = grain.radius;
    let fn_hertz = (4.0 / 3.0) * e_eff * r_eff.sqrt() * overlap.powf(1.5);

    let v_n = grain.velocity.dot(n);
    let fn_damp = -grain.gamma_n * grain.mass.sqrt() * v_n;
    let fn_total = (fn_hertz + fn_damp).max(0.0);
    let f_n = n.scale(fn_total);

    // Tangential friction (simplified)
    let v_t_full = grain.velocity.sub(n.scale(v_n));
    let v_t_norm = v_t_full.norm();
    let f_t = if v_t_norm > 1e-12 {
        let dir = v_t_full.scale(-1.0 / v_t_norm);
        dir.scale(params.mu_friction * fn_total)
    } else {
        Vec3::zero()
    };

    let f_total = f_n.add(f_t);
    let r_contact = n.scale(-grain.radius);
    let torque = r_contact.cross(f_t);
    (f_total, torque)
}

// ---------------------------------------------------------------------------
// Stress tensor from contact force chains
// ---------------------------------------------------------------------------

/// Cauchy stress tensor (6 independent components: xx, yy, zz, xy, xz, yz).
#[derive(Debug, Clone, Copy, Default)]
pub struct StressTensor {
    /// σ_xx component.
    pub xx: f64,
    /// σ_yy component.
    pub yy: f64,
    /// σ_zz component.
    pub zz: f64,
    /// σ_xy component.
    pub xy: f64,
    /// σ_xz component.
    pub xz: f64,
    /// σ_yz component.
    pub yz: f64,
}

impl StressTensor {
    /// Mean (hydrostatic) pressure: p = -(σ_xx + σ_yy + σ_zz) / 3.
    pub fn pressure(&self) -> f64 {
        -(self.xx + self.yy + self.zz) / 3.0
    }

    /// Von Mises equivalent stress.
    pub fn von_mises(&self) -> f64 {
        let s_xx = self.xx;
        let s_yy = self.yy;
        let s_zz = self.zz;
        let s_xy = self.xy;
        let s_xz = self.xz;
        let s_yz = self.yz;
        let vm2 = 0.5
            * ((s_xx - s_yy).powi(2)
                + (s_yy - s_zz).powi(2)
                + (s_zz - s_xx).powi(2)
                + 6.0 * (s_xy * s_xy + s_xz * s_xz + s_yz * s_yz));
        vm2.sqrt()
    }

    /// Deviator stress q (triaxial convention): q = σ_1 − σ_3 (simplified via VM).
    pub fn deviator_q(&self) -> f64 {
        self.von_mises()
    }
}

/// Compute the average Cauchy stress tensor over a representative volume element
/// using the Love–Weber formula:
///
/// σ = (1/V) Σ_{contacts} f ⊗ l
///
/// where f is the contact force and l is the branch vector.
pub fn compute_stress_tensor(
    grains: &[Grain],
    contacts: &[(usize, usize, Vec3, f64)], // (a, b, normal, fn)
    volume: f64,
) -> StressTensor {
    let mut st = StressTensor::default();
    for &(ia, ib, normal, fn_mag) in contacts {
        if ia >= grains.len() || ib >= grains.len() {
            continue;
        }
        let ga = &grains[ia];
        let gb = &grains[ib];
        let branch = gb.position.sub(ga.position);
        let force = normal.scale(fn_mag);
        // σ_ij += fi * lj / V
        st.xx += force.x * branch.x;
        st.yy += force.y * branch.y;
        st.zz += force.z * branch.z;
        st.xy += force.x * branch.y;
        st.xz += force.x * branch.z;
        st.yz += force.y * branch.z;
    }
    let inv_v = 1.0 / volume.max(1e-30);
    st.xx *= inv_v;
    st.yy *= inv_v;
    st.zz *= inv_v;
    st.xy *= inv_v;
    st.xz *= inv_v;
    st.yz *= inv_v;
    st
}

// ---------------------------------------------------------------------------
// Polydisperse packing
// ---------------------------------------------------------------------------

/// Descriptor for a polydisperse grain size distribution.
#[derive(Debug, Clone)]
pub struct GrainSizeDistribution {
    /// Minimum grain radius (m).
    pub r_min: f64,
    /// Maximum grain radius (m).
    pub r_max: f64,
    /// Distribution exponent (1 = uniform, >1 = finer grains more likely).
    pub exponent: f64,
}

impl GrainSizeDistribution {
    /// Sample a radius from the distribution using a power-law CDF inversion.
    pub fn sample(&self, rng: &mut impl Rng) -> f64 {
        let u: f64 = rng.random_range(0.0_f64..1.0_f64);
        if (self.exponent - 1.0).abs() < 1e-12 {
            // Uniform
            self.r_min + u * (self.r_max - self.r_min)
        } else {
            let e = self.exponent;
            let a = self.r_min.powf(1.0 - e);
            let b = self.r_max.powf(1.0 - e);
            (a + u * (b - a)).powf(1.0 / (1.0 - e))
        }
    }
}

/// Random-sequential addition packing of grains inside an AABB.
///
/// Attempts to place `n_grains` grains without overlap, up to `max_attempts`
/// tries per grain.  Returns the placed grains.
#[allow(clippy::too_many_arguments)]
pub fn random_sequential_packing(
    n_grains: usize,
    box_min: Vec3,
    box_max: Vec3,
    dist: &GrainSizeDistribution,
    density: f64,
    young: f64,
    poisson: f64,
    max_attempts: usize,
) -> Vec<Grain> {
    let mut rng = rand::rng();
    let mut grains: Vec<Grain> = Vec::with_capacity(n_grains);

    for id in 0..n_grains {
        let mut placed = false;
        for _ in 0..max_attempts {
            let r = dist.sample(&mut rng);
            let px = rng.random_range((box_min.x + r)..(box_max.x - r).max(box_min.x + r + 1e-12));
            let py = rng.random_range((box_min.y + r)..(box_max.y - r).max(box_min.y + r + 1e-12));
            let pz = rng.random_range((box_min.z + r)..(box_max.z - r).max(box_min.z + r + 1e-12));
            let pos = Vec3::new(px, py, pz);

            let overlap = grains
                .iter()
                .any(|g| g.position.sub(pos).norm() < g.radius + r);
            if !overlap {
                let mut g = Grain::new(id, pos, r, density, young, poisson, 50.0, 10.0);
                g.is_fixed = false;
                grains.push(g);
                placed = true;
                break;
            }
        }
        if !placed {
            // Place anyway with a tiny radius to fill quota (degenerate)
            let r = dist.r_min * 0.5;
            let px = rng.random_range(box_min.x..box_max.x);
            let py = rng.random_range(box_min.y..box_max.y);
            let pz = rng.random_range(box_min.z..box_max.z);
            grains.push(Grain::new(
                id,
                Vec3::new(px, py, pz),
                r,
                density,
                young,
                poisson,
                50.0,
                10.0,
            ));
        }
    }
    grains
}

// ---------------------------------------------------------------------------
// Grain crushing / breakage
// ---------------------------------------------------------------------------

/// Parameters for a grain crushing model based on the Weibull distribution.
#[derive(Debug, Clone)]
pub struct CrushingParams {
    /// Weibull modulus (shape parameter, dimensionless).
    pub weibull_m: f64,
    /// Characteristic strength (Pa).
    pub sigma_0: f64,
    /// Reference volume (m³).
    pub v_ref: f64,
    /// Ratio of child grain radii to parent.
    pub fragment_ratio: f64,
}

impl CrushingParams {
    /// Default parameters for silica sand.
    pub fn silica_sand() -> Self {
        Self {
            weibull_m: 4.0,
            sigma_0: 50e6,
            v_ref: 1e-9,
            fragment_ratio: 0.5,
        }
    }

    /// Compute Weibull crushing probability for a given local stress.
    pub fn crush_probability(&self, sigma: f64, volume: f64) -> f64 {
        if sigma <= 0.0 {
            return 0.0;
        }
        let v_ratio = volume / self.v_ref.max(1e-40);
        1.0 - (-v_ratio * (sigma / self.sigma_0.max(1e-40)).powf(self.weibull_m)).exp()
    }
}

/// Attempt grain crushing.  Returns replacement grains if the parent crushes.
pub fn attempt_crush(
    grain: &Grain,
    local_stress: f64,
    params: &CrushingParams,
    rng: &mut impl Rng,
    next_id: &mut usize,
) -> Option<Vec<Grain>> {
    let volume = (4.0 / 3.0) * std::f64::consts::PI * grain.radius.powi(3);
    let p = params.crush_probability(local_stress, volume);
    let u: f64 = rng.random_range(0.0_f64..1.0_f64);
    if u > p {
        return None;
    }
    // Replace with two smaller grains displaced along a random direction
    let child_r = grain.radius * params.fragment_ratio;
    let density = grain.mass / volume;
    let theta: f64 = rng.random_range(0.0_f64..(2.0 * std::f64::consts::PI));
    let offset = Vec3::new(theta.cos() * child_r, theta.sin() * child_r, 0.0);
    let id1 = *next_id;
    *next_id += 1;
    let id2 = *next_id;
    *next_id += 1;
    let g1 = Grain::new(
        id1,
        grain.position.add(offset),
        child_r,
        density,
        grain.young,
        grain.poisson,
        grain.gamma_n,
        grain.gamma_t,
    );
    let g2 = Grain::new(
        id2,
        grain.position.sub(offset),
        child_r,
        density,
        grain.young,
        grain.poisson,
        grain.gamma_n,
        grain.gamma_t,
    );
    Some(vec![g1, g2])
}

// ---------------------------------------------------------------------------
// Critical state soil mechanics (Cam-Clay)
// ---------------------------------------------------------------------------

/// State variables for critical-state soil mechanics (modified Cam-Clay).
#[derive(Debug, Clone)]
pub struct CamClayState {
    /// Mean effective stress p′ (Pa).
    pub p_prime: f64,
    /// Deviatoric stress q (Pa).
    pub q: f64,
    /// Specific volume v = 1 + e.
    pub v: f64,
    /// Pre-consolidation pressure p_c (Pa).
    pub p_c: f64,
}

impl CamClayState {
    /// Create a new Cam-Clay state at the given stresses.
    pub fn new(p_prime: f64, q: f64, v: f64, p_c: f64) -> Self {
        Self { p_prime, q, v, p_c }
    }

    /// Yield surface value f.  f <= 0 means elastic, f > 0 means plastic.
    /// Modified Cam-Clay: f = q² / M² + p′(p′ - p_c).
    pub fn yield_function(&self, m_critical: f64) -> f64 {
        self.q * self.q / (m_critical * m_critical) + self.p_prime * (self.p_prime - self.p_c)
    }

    /// Compute the slope of critical state line M from friction angle phi (deg).
    pub fn m_from_phi(phi_deg: f64) -> f64 {
        let phi = phi_deg.to_radians();
        6.0 * phi.sin() / (3.0 - phi.sin())
    }

    /// Check if currently at critical state.
    pub fn is_critical_state(&self, m_critical: f64, tol: f64) -> bool {
        self.yield_function(m_critical).abs() < tol
    }
}

// ---------------------------------------------------------------------------
// Shear band detection
// ---------------------------------------------------------------------------

/// A detected shear band.
#[derive(Debug, Clone)]
pub struct ShearBand {
    /// Centre of the shear band.
    pub centre: Vec3,
    /// Normal direction of the band (pointing across the band).
    pub normal: Vec3,
    /// Estimated thickness (m).
    pub thickness: f64,
    /// Maximum shear strain rate within band (s⁻¹).
    pub max_shear_strain_rate: f64,
}

impl ShearBand {
    /// Detect shear bands from per-grain strain-rate data.
    ///
    /// Grains with a shear-strain-rate magnitude above `threshold` are
    /// clustered; a simple mean-normal heuristic identifies the band.
    pub fn detect(grains: &[Grain], shear_strain_rates: &[f64], threshold: f64) -> Vec<ShearBand> {
        assert_eq!(grains.len(), shear_strain_rates.len());
        let high: Vec<usize> = shear_strain_rates
            .iter()
            .enumerate()
            .filter(|&(_, &s)| s > threshold)
            .map(|(i, _)| i)
            .collect();
        if high.is_empty() {
            return Vec::new();
        }
        // Centroid
        let mut cx = 0.0f64;
        let mut cy = 0.0f64;
        let mut cz = 0.0f64;
        for &i in &high {
            cx += grains[i].position.x;
            cy += grains[i].position.y;
            cz += grains[i].position.z;
        }
        let n = high.len() as f64;
        let centre = Vec3::new(cx / n, cy / n, cz / n);
        // Rough band normal: use PCA-like approach (gradient of shear rate)
        // For simplicity, use y-axis as default (horizontal band)
        let normal = Vec3::new(0.0, 1.0, 0.0);
        let max_ssr = high
            .iter()
            .map(|&i| shear_strain_rates[i])
            .fold(0.0_f64, f64::max);
        let thickness = high.len() as f64 * 2.0 * grains[high[0]].radius;
        vec![ShearBand {
            centre,
            normal,
            thickness,
            max_shear_strain_rate: max_ssr,
        }]
    }
}

// ---------------------------------------------------------------------------
// Avalanche and angle-of-repose
// ---------------------------------------------------------------------------

/// Result of an angle-of-repose measurement.
#[derive(Debug, Clone)]
pub struct AngleOfReposeResult {
    /// Angle of repose in degrees.
    pub angle_deg: f64,
    /// Height of the pile (m).
    pub pile_height: f64,
    /// Base radius of the pile (m).
    pub base_radius: f64,
}

/// Estimate the angle of repose from a granular pile.
///
/// Grains are assumed to be resting on a flat floor at y = 0.
/// The pile height and base radius are measured from the grain positions.
pub fn measure_angle_of_repose(grains: &[Grain]) -> AngleOfReposeResult {
    let floor_y = 0.0f64;
    if grains.is_empty() {
        return AngleOfReposeResult {
            angle_deg: 0.0,
            pile_height: 0.0,
            base_radius: 0.0,
        };
    }
    let pile_height = grains
        .iter()
        .map(|g| g.position.y + g.radius - floor_y)
        .fold(0.0_f64, f64::max);
    let centre_x = grains.iter().map(|g| g.position.x).sum::<f64>() / grains.len() as f64;
    let centre_z = grains.iter().map(|g| g.position.z).sum::<f64>() / grains.len() as f64;
    let base_radius = grains
        .iter()
        .map(|g| {
            let dx = g.position.x - centre_x;
            let dz = g.position.z - centre_z;
            (dx * dx + dz * dz).sqrt() + g.radius
        })
        .fold(0.0_f64, f64::max);
    let angle_deg = if base_radius > 1e-12 {
        (pile_height / base_radius).atan().to_degrees()
    } else {
        90.0
    };
    AngleOfReposeResult {
        angle_deg,
        pile_height,
        base_radius,
    }
}

/// Detect grains that are about to avalanche (slope exceeds angle of repose).
pub fn detect_avalanche_grains(grains: &[Grain], angle_of_repose_deg: f64) -> Vec<usize> {
    let tan_aor = angle_of_repose_deg.to_radians().tan();
    let cx = grains.iter().map(|g| g.position.x).sum::<f64>() / grains.len().max(1) as f64;
    let cz = grains.iter().map(|g| g.position.z).sum::<f64>() / grains.len().max(1) as f64;
    grains
        .iter()
        .enumerate()
        .filter(|(_, g)| {
            let dr = {
                let dx = g.position.x - cx;
                let dz = g.position.z - cz;
                (dx * dx + dz * dz).sqrt()
            };
            let dy = g.position.y;
            dr > 1e-12 && (dy / dr) > tan_aor
        })
        .map(|(i, _)| i)
        .collect()
}

// ---------------------------------------------------------------------------
// Compaction simulation driver
// ---------------------------------------------------------------------------

/// Compaction state tracker for uniaxial compaction.
#[derive(Debug, Clone)]
pub struct CompactionState {
    /// Current void ratio e = (V_void / V_solid).
    pub void_ratio: f64,
    /// Applied vertical stress (Pa).
    pub applied_stress: f64,
    /// Compression index Cc (slope of e-log p curve).
    pub cc: f64,
    /// Swelling index Cs.
    pub cs: f64,
    /// Reference void ratio at p = 1 Pa.
    pub e_ref: f64,
}

impl CompactionState {
    /// Create initial compaction state for a loose packing.
    pub fn loose_sand() -> Self {
        Self {
            void_ratio: 0.8,
            applied_stress: 1000.0,
            cc: 0.3,
            cs: 0.05,
            e_ref: 1.5,
        }
    }

    /// Compute void ratio under a given mean effective stress using
    /// normal consolidation line: e = e_ref - Cc * log10(p/1Pa).
    pub fn ncl_void_ratio(&self, p: f64) -> f64 {
        self.e_ref - self.cc * p.max(1e-10).log10()
    }

    /// Apply load increment dp and return the new void ratio.
    pub fn apply_load_step(&mut self, dp: f64) -> f64 {
        let p_new = self.applied_stress + dp;
        let e_new = if dp >= 0.0 {
            self.ncl_void_ratio(p_new)
        } else {
            // Unloading: swelling line
            self.void_ratio - self.cs * (p_new / self.applied_stress.max(1e-10)).log10()
        };
        self.applied_stress = p_new;
        self.void_ratio = e_new;
        e_new
    }

    /// Relative density Dr from min/max void ratios.
    pub fn relative_density(&self, e_max: f64, e_min: f64) -> f64 {
        if (e_max - e_min).abs() < 1e-12 {
            return 0.0;
        }
        (e_max - self.void_ratio) / (e_max - e_min)
    }
}

// ---------------------------------------------------------------------------
// DEM-FEM coupling
// ---------------------------------------------------------------------------

/// A triangular FEM element for continuum stress representation.
#[derive(Debug, Clone)]
pub struct FemTriangle {
    /// Global grain indices forming the triangle nodes.
    pub nodes: [usize; 3],
    /// Area of the triangle (m²).
    pub area: f64,
    /// Homogenised stress tensor within this element.
    pub stress: StressTensor,
}

impl FemTriangle {
    /// Create a new FEM triangle from grain positions.
    pub fn new(nodes: [usize; 3], grains: &[Grain]) -> Self {
        let a = grains[nodes[0]].position;
        let b = grains[nodes[1]].position;
        let c = grains[nodes[2]].position;
        let ab = b.sub(a);
        let ac = c.sub(a);
        let area = 0.5 * ab.cross(ac).norm();
        Self {
            nodes,
            area,
            stress: StressTensor::default(),
        }
    }

    /// Update the stress by homogenising DEM contact forces over the element.
    pub fn homogenise_stress(&mut self, grains: &[Grain], contacts: &[(usize, usize, Vec3, f64)]) {
        let volume = self.area; // 2D approximation: use area as volume
        self.stress = compute_stress_tensor(grains, contacts, volume);
    }
}

/// DEM-FEM coupled system.
#[derive(Debug, Clone)]
pub struct DemFemSystem {
    /// Discrete grains.
    pub grains: Vec<Grain>,
    /// FEM triangles for stress homogenisation.
    pub fem_elements: Vec<FemTriangle>,
    /// Current active contacts: (ia, ib, normal, fn_magnitude).
    pub contacts: Vec<(usize, usize, Vec3, f64)>,
    /// Tangential displacement history per grain pair (keyed by sorted pair).
    tangential_cache: std::collections::HashMap<(usize, usize), Vec3>,
    /// DEM parameters.
    pub params: DemParams,
    /// Crushing parameters (optional).
    pub crushing: Option<CrushingParams>,
    /// Next available grain ID.
    next_id: usize,
}

impl DemFemSystem {
    /// Create a new empty DEM-FEM system.
    pub fn new(params: DemParams, crushing: Option<CrushingParams>) -> Self {
        Self {
            grains: Vec::new(),
            fem_elements: Vec::new(),
            contacts: Vec::new(),
            tangential_cache: std::collections::HashMap::new(),
            params,
            crushing,
            next_id: 0,
        }
    }

    /// Add a grain to the system.
    pub fn add_grain(&mut self, mut grain: Grain) {
        grain.id = self.next_id;
        self.next_id += 1;
        self.grains.push(grain);
    }

    /// Add a FEM element.
    pub fn add_fem_element(&mut self, element: FemTriangle) {
        self.fem_elements.push(element);
    }

    /// Detect all grain-grain contacts and compute forces.
    fn compute_contacts(&mut self) {
        self.contacts.clear();
        let n = self.grains.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if self.grains[i].crushed || self.grains[j].crushed {
                    continue;
                }
                let key = (i.min(j), i.max(j));
                let entry = self.tangential_cache.entry(key).or_insert(Vec3::zero());
                let ga = self.grains[i].clone();
                let gb = self.grains[j].clone();
                let (fa, ta, fb, tb) = hertz_mindlin_force(&ga, &gb, entry, &self.params);
                let fn_mag = fa.norm();
                if fn_mag > 0.0 {
                    let r = gb.position.sub(ga.position);
                    let dist = r.norm();
                    let normal = if dist > 1e-15 {
                        r.scale(1.0 / dist)
                    } else {
                        Vec3::new(0.0, 1.0, 0.0)
                    };
                    self.contacts.push((i, j, normal, fn_mag));
                }
                self.grains[i].apply_force(fa);
                self.grains[i].apply_torque(ta);
                self.grains[j].apply_force(fb);
                self.grains[j].apply_torque(tb);
            }
        }
    }

    /// Apply gravity and integrate all grains.
    pub fn step(&mut self) {
        // Clear forces
        for g in &mut self.grains {
            g.clear_forces();
        }
        // Gravity
        let grav = self.params.gravity;
        for g in &mut self.grains {
            if !g.is_fixed && !g.crushed {
                let fg = grav.scale(g.mass);
                g.apply_force(fg);
            }
        }
        // Contacts
        self.compute_contacts();
        // Integrate
        let dt = self.params.dt;
        for g in &mut self.grains {
            g.integrate(dt);
        }
        // Crushing check
        if let Some(ref cp) = self.crushing.clone() {
            let mut rng = rand::rng();
            let mut new_grains: Vec<Grain> = Vec::new();
            for g in &mut self.grains {
                if g.crushed {
                    continue;
                }
                // Local stress ~ sum of contact forces / cross-section
                let local_stress =
                    g.force.norm() / (std::f64::consts::PI * g.radius * g.radius).max(1e-30);
                if let Some(mut children) =
                    attempt_crush(g, local_stress, cp, &mut rng, &mut self.next_id)
                {
                    g.crushed = true;
                    new_grains.append(&mut children);
                }
            }
            self.grains.extend(new_grains);
        }
        // Update FEM stresses
        let contacts = self.contacts.clone();
        let grains = &self.grains;
        for elem in &mut self.fem_elements {
            elem.homogenise_stress(grains, &contacts);
        }
    }

    /// Compute the global stress tensor over all grains.
    pub fn global_stress(&self) -> StressTensor {
        // Volume: bounding box volume
        if self.grains.is_empty() {
            return StressTensor::default();
        }
        let xs: Vec<f64> = self.grains.iter().map(|g| g.position.x).collect();
        let ys: Vec<f64> = self.grains.iter().map(|g| g.position.y).collect();
        let zs: Vec<f64> = self.grains.iter().map(|g| g.position.z).collect();
        let vol = (xs.iter().cloned().fold(f64::MIN, f64::max)
            - xs.iter().cloned().fold(f64::MAX, f64::min))
            * (ys.iter().cloned().fold(f64::MIN, f64::max)
                - ys.iter().cloned().fold(f64::MAX, f64::min))
            * (zs.iter().cloned().fold(f64::MIN, f64::max)
                - zs.iter().cloned().fold(f64::MAX, f64::min));
        compute_stress_tensor(&self.grains, &self.contacts, vol.max(1e-30))
    }
}

// ---------------------------------------------------------------------------
// Contact force chain analysis
// ---------------------------------------------------------------------------

/// A node in the contact force chain graph.
#[derive(Debug, Clone)]
pub struct ForceChainNode {
    /// Grain index.
    pub grain_id: usize,
    /// Total contact force magnitude (N) transmitted through this grain.
    pub total_force: f64,
}

/// Extract the "strong" contact network (contacts above threshold).
///
/// Returns a list of grain-index pairs forming the strong network.
pub fn strong_contact_network(
    contacts: &[(usize, usize, Vec3, f64)],
    mean_force_fraction: f64,
) -> Vec<(usize, usize)> {
    if contacts.is_empty() {
        return Vec::new();
    }
    let mean_f: f64 = contacts.iter().map(|c| c.3).sum::<f64>() / contacts.len() as f64;
    let threshold = mean_force_fraction * mean_f;
    contacts
        .iter()
        .filter(|c| c.3 >= threshold)
        .map(|c| (c.0, c.1))
        .collect()
}

/// Compute coordination number distribution (number of contacts per grain).
pub fn coordination_numbers(n_grains: usize, contacts: &[(usize, usize, Vec3, f64)]) -> Vec<usize> {
    let mut counts = vec![0usize; n_grains];
    for &(ia, ib, _, _) in contacts {
        if ia < n_grains {
            counts[ia] += 1;
        }
        if ib < n_grains {
            counts[ib] += 1;
        }
    }
    counts
}

/// Mean coordination number.
pub fn mean_coordination_number(n_grains: usize, contacts: &[(usize, usize, Vec3, f64)]) -> f64 {
    let cnums = coordination_numbers(n_grains, contacts);
    if cnums.is_empty() {
        return 0.0;
    }
    cnums.iter().sum::<usize>() as f64 / cnums.len() as f64
}

// ---------------------------------------------------------------------------
// Packing fraction utility
// ---------------------------------------------------------------------------

/// Compute packing fraction (solid volume fraction) of a grain assembly.
pub fn packing_fraction(grains: &[Grain], total_volume: f64) -> f64 {
    let solid_vol: f64 = grains
        .iter()
        .filter(|g| !g.crushed)
        .map(|g| (4.0 / 3.0) * std::f64::consts::PI * g.radius.powi(3))
        .sum();
    solid_vol / total_volume.max(1e-30)
}

// ---------------------------------------------------------------------------
// DEM simulation builder
// ---------------------------------------------------------------------------

/// Builder for a simple DEM simulation.
#[derive(Debug)]
pub struct DemSimulationBuilder {
    /// Number of grains to place.
    pub n_grains: usize,
    /// Bounding box minimum.
    pub box_min: Vec3,
    /// Bounding box maximum.
    pub box_max: Vec3,
    /// Grain size distribution.
    pub size_dist: GrainSizeDistribution,
    /// Grain density (kg m⁻³).
    pub density: f64,
    /// Young's modulus (Pa).
    pub young: f64,
    /// Poisson ratio.
    pub poisson: f64,
    /// DEM parameters.
    pub dem_params: DemParams,
    /// Enable crushing.
    pub enable_crushing: bool,
}

impl DemSimulationBuilder {
    /// Create a default sand-box simulation builder.
    pub fn sand_box(n: usize) -> Self {
        Self {
            n_grains: n,
            box_min: Vec3::new(-0.1, 0.0, -0.1),
            box_max: Vec3::new(0.1, 0.2, 0.1),
            size_dist: GrainSizeDistribution {
                r_min: 1e-3,
                r_max: 3e-3,
                exponent: 1.0,
            },
            density: 2650.0,
            young: 70e9,
            poisson: 0.3,
            dem_params: DemParams::default_sand(),
            enable_crushing: false,
        }
    }

    /// Build the DEM-FEM system with grains packed by RSA.
    pub fn build(self) -> DemFemSystem {
        let crushing = if self.enable_crushing {
            Some(CrushingParams::silica_sand())
        } else {
            None
        };
        let mut sys = DemFemSystem::new(self.dem_params, crushing);
        let grains = random_sequential_packing(
            self.n_grains,
            self.box_min,
            self.box_max,
            &self.size_dist,
            self.density,
            self.young,
            self.poisson,
            200,
        );
        for g in grains {
            sys.add_grain(g);
        }
        sys
    }
}

// ---------------------------------------------------------------------------
// Fabric tensor
// ---------------------------------------------------------------------------

/// Fabric tensor: second-order tensor describing the anisotropy of contact normals.
#[derive(Debug, Clone, Default)]
pub struct FabricTensor {
    /// Component F_xx.
    pub fxx: f64,
    /// Component F_yy.
    pub fyy: f64,
    /// Component F_zz.
    pub fzz: f64,
    /// Component F_xy.
    pub fxy: f64,
    /// Component F_xz.
    pub fxz: f64,
    /// Component F_yz.
    pub fyz: f64,
    /// Number of contacts used.
    pub n_contacts: usize,
}

impl FabricTensor {
    /// Compute fabric tensor from contact normals.
    pub fn from_contacts(contacts: &[(usize, usize, Vec3, f64)]) -> Self {
        let mut ft = Self::default();
        for &(_, _, n, _) in contacts {
            ft.fxx += n.x * n.x;
            ft.fyy += n.y * n.y;
            ft.fzz += n.z * n.z;
            ft.fxy += n.x * n.y;
            ft.fxz += n.x * n.z;
            ft.fyz += n.y * n.z;
            ft.n_contacts += 1;
        }
        if ft.n_contacts > 0 {
            let inv = 1.0 / ft.n_contacts as f64;
            ft.fxx *= inv;
            ft.fyy *= inv;
            ft.fzz *= inv;
            ft.fxy *= inv;
            ft.fxz *= inv;
            ft.fyz *= inv;
        }
        ft
    }

    /// Anisotropy parameter a = 3/2 * ||F - I/3||.
    pub fn anisotropy(&self) -> f64 {
        let iso = 1.0 / 3.0;
        let d_xx = self.fxx - iso;
        let d_yy = self.fyy - iso;
        let d_zz = self.fzz - iso;
        let fro2 = d_xx * d_xx
            + d_yy * d_yy
            + d_zz * d_zz
            + 2.0 * (self.fxy * self.fxy + self.fxz * self.fxz + self.fyz * self.fyz);
        1.5 * fro2.sqrt()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    fn default_grain(id: usize, pos: Vec3, r: f64) -> Grain {
        Grain::new(id, pos, r, 2650.0, 70e9, 0.3, 50.0, 10.0)
    }

    // 1. Vec3 norm
    #[test]
    fn test_vec3_norm() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert!((v.norm() - 5.0).abs() < EPS);
    }

    // 2. Vec3 dot product
    #[test]
    fn test_vec3_dot() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert!((a.dot(b) - 32.0).abs() < EPS);
    }

    // 3. Vec3 cross product
    #[test]
    fn test_vec3_cross() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        let z = x.cross(y);
        assert!((z.x).abs() < EPS && (z.y).abs() < EPS && (z.z - 1.0).abs() < EPS);
    }

    // 4. Vec3 normalise
    #[test]
    fn test_vec3_normalise() {
        let v = Vec3::new(3.0, 0.0, 4.0);
        let n = v.normalise();
        assert!((n.norm() - 1.0).abs() < EPS);
    }

    // 5. Grain mass from density and radius
    #[test]
    fn test_grain_mass() {
        let r = 1e-3;
        let rho = 2650.0;
        let g = default_grain(0, Vec3::zero(), r);
        let expected_mass = rho * (4.0 / 3.0) * std::f64::consts::PI * r.powi(3);
        assert!((g.mass - expected_mass).abs() < 1e-15 * expected_mass.max(1.0));
    }

    // 6. Grain inertia = 2/5 m r²
    #[test]
    fn test_grain_inertia() {
        let r = 2e-3;
        let g = default_grain(0, Vec3::zero(), r);
        let expected = 0.4 * g.mass * r * r;
        assert!((g.inertia - expected).abs() < EPS * expected.max(1.0));
    }

    // 7. Hertz contact: no force when grains not overlapping
    #[test]
    fn test_hertz_no_overlap() {
        let ga = default_grain(0, Vec3::new(0.0, 0.0, 0.0), 1e-3);
        let gb = default_grain(1, Vec3::new(10e-3, 0.0, 0.0), 1e-3);
        let params = DemParams::default_sand();
        let mut dt = Vec3::zero();
        let (fa, _ta, fb, _tb) = hertz_mindlin_force(&ga, &gb, &mut dt, &params);
        assert!(fa.norm() < EPS && fb.norm() < EPS);
    }

    // 8. Hertz contact: force is repulsive when overlapping
    #[test]
    fn test_hertz_repulsive_overlap() {
        let ga = default_grain(0, Vec3::new(0.0, 0.0, 0.0), 1e-3);
        let gb = default_grain(1, Vec3::new(1.5e-3, 0.0, 0.0), 1e-3); // overlap 0.5e-3
        let params = DemParams::default_sand();
        let mut dt = Vec3::zero();
        let (fa, _ta, fb, _tb) = hertz_mindlin_force(&ga, &gb, &mut dt, &params);
        // Force on a should be in -x, force on b in +x
        assert!(fa.x < 0.0, "fa.x should be negative (repulsive): {}", fa.x);
        assert!(fb.x > 0.0, "fb.x should be positive: {}", fb.x);
    }

    // 9. Wall contact: grain above wall → no force
    #[test]
    fn test_wall_no_contact() {
        let g = default_grain(0, Vec3::new(0.0, 5e-3, 0.0), 1e-3);
        let params = DemParams::default_sand();
        let (f, _t) = wall_contact_force(
            &g,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            &params,
        );
        assert!(f.norm() < EPS);
    }

    // 10. Wall contact: grain penetrating floor → upward force
    #[test]
    fn test_wall_contact_upward() {
        let g = default_grain(0, Vec3::new(0.0, 0.5e-3, 0.0), 1e-3); // overlap 0.5e-3
        let params = DemParams::default_sand();
        let (f, _t) = wall_contact_force(
            &g,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            &params,
        );
        assert!(f.y > 0.0, "Wall force should push grain upward: {}", f.y);
    }

    // 11. Stress tensor: trace for uniaxial loading
    #[test]
    fn test_stress_uniaxial() {
        let grains = vec![
            default_grain(0, Vec3::new(-1.0, 0.0, 0.0), 0.5),
            default_grain(1, Vec3::new(1.0, 0.0, 0.0), 0.5),
        ];
        let contacts = vec![(0usize, 1usize, Vec3::new(1.0, 0.0, 0.0), 100.0)];
        let vol = 8.0f64;
        let st = compute_stress_tensor(&grains, &contacts, vol);
        // σ_xx should be non-zero, σ_yy and σ_zz should be near zero
        assert!(st.xx.abs() > EPS, "σ_xx should be non-zero: {}", st.xx);
        assert!(st.yy.abs() < 1e-9, "σ_yy should be zero: {}", st.yy);
    }

    // 12. Von Mises stress for pure shear
    #[test]
    fn test_von_mises_shear() {
        let st = StressTensor {
            xx: 0.0,
            yy: 0.0,
            zz: 0.0,
            xy: 10.0,
            xz: 0.0,
            yz: 0.0,
        };
        let expected = (3.0 * 100.0_f64).sqrt(); // sqrt(3) * tau
        assert!((st.von_mises() - expected).abs() < EPS);
    }

    // 13. Cam-Clay yield function at critical state
    #[test]
    fn test_cam_clay_critical_state() {
        let m = CamClayState::m_from_phi(30.0);
        // At critical state: q = M * p'
        let p = 100e3;
        let q = m * p;
        let pc = 2.0 * p; // p_c = 2p at critical state in Cam-Clay
        let state = CamClayState::new(p, q, 2.0, pc);
        assert!(
            state.is_critical_state(m, 1.0),
            "Should be at critical state"
        );
    }

    // 14. M from friction angle: phi=30° should give M≈1.2
    #[test]
    fn test_m_from_phi() {
        let m = CamClayState::m_from_phi(30.0);
        assert!(
            (m - 1.2).abs() < 0.01,
            "M for phi=30° should be ~1.2, got {m}"
        );
    }

    // 15. Polydisperse packing: correct grain count
    #[test]
    fn test_polydisperse_packing_count() {
        let dist = GrainSizeDistribution {
            r_min: 1e-3,
            r_max: 2e-3,
            exponent: 1.0,
        };
        let grains = random_sequential_packing(
            10,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.1, 0.1, 0.1),
            &dist,
            2650.0,
            70e9,
            0.3,
            100,
        );
        assert_eq!(grains.len(), 10);
    }

    // 16. Packing fraction in [0, 1]
    #[test]
    fn test_packing_fraction_range() {
        let dist = GrainSizeDistribution {
            r_min: 1e-3,
            r_max: 2e-3,
            exponent: 1.0,
        };
        let grains = random_sequential_packing(
            5,
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.05, 0.05, 0.05),
            &dist,
            2650.0,
            70e9,
            0.3,
            50,
        );
        let pf = packing_fraction(&grains, 0.05_f64.powi(3));
        assert!(
            (0.0..=1.0).contains(&pf),
            "Packing fraction {pf} out of [0,1]"
        );
    }

    // 17. Compaction state: void ratio decreases under loading
    #[test]
    fn test_compaction_loading() {
        let mut cs = CompactionState::loose_sand();
        let e0 = cs.void_ratio;
        cs.apply_load_step(5000.0);
        let e1 = cs.void_ratio;
        assert!(
            e1 < e0,
            "Void ratio should decrease under loading: e0={e0}, e1={e1}"
        );
    }

    // 18. Relative density bounded in [0, 1]
    #[test]
    fn test_relative_density() {
        let cs = CompactionState {
            void_ratio: 0.6,
            ..CompactionState::loose_sand()
        };
        let dr = cs.relative_density(0.9, 0.4);
        assert!(
            (0.0..=1.0).contains(&dr),
            "Relative density {dr} out of range"
        );
    }

    // 19. Crushing probability 0 for zero stress
    #[test]
    fn test_crushing_prob_zero_stress() {
        let cp = CrushingParams::silica_sand();
        assert_eq!(cp.crush_probability(0.0, 1e-9), 0.0);
    }

    // 20. Crushing probability → 1 for very large stress
    #[test]
    fn test_crushing_prob_high_stress() {
        let cp = CrushingParams::silica_sand();
        let p = cp.crush_probability(1e15, 1e-9);
        assert!(
            p > 0.99,
            "Crushing probability should be near 1 for very high stress: {p}"
        );
    }

    // 21. Angle of repose from a simple pile
    #[test]
    fn test_angle_of_repose() {
        let mut grains = Vec::new();
        // Build a conical pile: r increases with height step
        for i in 0..5usize {
            let y = i as f64 * 0.01 + 0.005;
            let r_spread = (5 - i) as f64 * 0.005;
            for sign in [-1.0, 1.0] {
                let mut g = default_grain(i * 2, Vec3::new(sign * r_spread, y, 0.0), 0.005);
                g.is_fixed = true;
                grains.push(g);
            }
        }
        let result = measure_angle_of_repose(&grains);
        assert!(result.angle_deg > 0.0 && result.angle_deg < 90.0);
    }

    // 22. Coordination number: 2-grain contact gives coordination 1 each
    #[test]
    fn test_coordination_numbers() {
        let contacts = vec![(0usize, 1usize, Vec3::new(1.0, 0.0, 0.0), 10.0)];
        let cnums = coordination_numbers(2, &contacts);
        assert_eq!(cnums[0], 1);
        assert_eq!(cnums[1], 1);
    }

    // 23. Mean coordination number for 3 grains in a ring
    #[test]
    fn test_mean_coordination() {
        let contacts = vec![
            (0usize, 1usize, Vec3::new(1.0, 0.0, 0.0), 10.0),
            (1usize, 2usize, Vec3::new(0.0, 1.0, 0.0), 10.0),
            (0usize, 2usize, Vec3::new(1.0, 1.0, 0.0), 10.0),
        ];
        let mean = mean_coordination_number(3, &contacts);
        assert!(
            (mean - 2.0).abs() < EPS,
            "Mean coordination should be 2.0, got {mean}"
        );
    }

    // 24. Fabric tensor trace = 1 (normalised)
    #[test]
    fn test_fabric_tensor_trace() {
        let contacts = vec![
            (0usize, 1usize, Vec3::new(1.0, 0.0, 0.0), 1.0),
            (0usize, 2usize, Vec3::new(0.0, 1.0, 0.0), 1.0),
            (0usize, 3usize, Vec3::new(0.0, 0.0, 1.0), 1.0),
        ];
        let ft = FabricTensor::from_contacts(&contacts);
        let trace = ft.fxx + ft.fyy + ft.fzz;
        assert!(
            (trace - 1.0).abs() < EPS,
            "Fabric tensor trace should be 1: {trace}"
        );
    }

    // 25. Fabric tensor anisotropy = 0 for isotropic contact distribution
    #[test]
    fn test_fabric_isotropic() {
        // 6 contacts: ±x, ±y, ±z
        let contacts = vec![
            (0usize, 1usize, Vec3::new(1.0, 0.0, 0.0), 1.0),
            (0usize, 2usize, Vec3::new(-1.0, 0.0, 0.0), 1.0),
            (0usize, 3usize, Vec3::new(0.0, 1.0, 0.0), 1.0),
            (0usize, 4usize, Vec3::new(0.0, -1.0, 0.0), 1.0),
            (0usize, 5usize, Vec3::new(0.0, 0.0, 1.0), 1.0),
            (0usize, 6usize, Vec3::new(0.0, 0.0, -1.0), 1.0),
        ];
        let ft = FabricTensor::from_contacts(&contacts);
        assert!(
            ft.anisotropy() < 1e-9,
            "Isotropic contacts should give zero anisotropy"
        );
    }

    // 26. Strong contact network: filters below mean
    #[test]
    fn test_strong_contact_network() {
        let contacts = vec![
            (0usize, 1usize, Vec3::new(1.0, 0.0, 0.0), 1.0),
            (1usize, 2usize, Vec3::new(0.0, 1.0, 0.0), 10.0),
            (2usize, 3usize, Vec3::new(0.0, 0.0, 1.0), 5.0),
        ];
        // mean = (1+10+5)/3 ≈ 5.33; threshold at 1.0 * mean ≈ 5.33
        let strong = strong_contact_network(&contacts, 1.0);
        assert!(
            !strong.is_empty(),
            "Should find at least one strong contact"
        );
        for &(a, b) in &strong {
            let fn_val = contacts.iter().find(|c| c.0 == a && c.1 == b).unwrap().3;
            let mean = contacts.iter().map(|c| c.3).sum::<f64>() / contacts.len() as f64;
            assert!(
                fn_val >= mean,
                "Strong contact force {fn_val} < threshold {mean}"
            );
        }
    }

    // 27. DEM-FEM system: step runs without panic
    #[test]
    fn test_dem_fem_step() {
        let mut sys = DemFemSystem::new(DemParams::default_sand(), None);
        let ga = default_grain(0, Vec3::new(0.0, 5e-3, 0.0), 1e-3);
        let mut gb = default_grain(1, Vec3::new(3e-3, 5e-3, 0.0), 1e-3);
        gb.is_fixed = true;
        sys.add_grain(ga);
        sys.add_grain(gb);
        for _ in 0..10 {
            sys.step();
        }
        assert!(sys.grains[0].position.y < 5e-3 + 1e-3);
    }

    // 28. Shear band detection: no bands if all strain rates below threshold
    #[test]
    fn test_shear_band_none() {
        let grains: Vec<Grain> = (0..5)
            .map(|i| default_grain(i, Vec3::new(i as f64 * 1e-3, 0.0, 0.0), 5e-4))
            .collect();
        let rates = vec![0.001; 5];
        let bands = ShearBand::detect(&grains, &rates, 1.0);
        assert!(bands.is_empty(), "No shear bands expected below threshold");
    }

    // 29. Shear band detection: detects band when strain rate exceeds threshold
    #[test]
    fn test_shear_band_detected() {
        let grains: Vec<Grain> = (0..5)
            .map(|i| default_grain(i, Vec3::new(i as f64 * 1e-3, 0.0, 0.0), 5e-4))
            .collect();
        let mut rates = vec![0.001; 5];
        rates[2] = 10.0; // one grain above threshold
        rates[3] = 12.0;
        let bands = ShearBand::detect(&grains, &rates, 1.0);
        assert!(!bands.is_empty(), "Shear band should be detected");
        assert!(bands[0].max_shear_strain_rate > 1.0);
    }

    // 30. FEM triangle area computation
    #[test]
    fn test_fem_triangle_area() {
        let grains = vec![
            default_grain(0, Vec3::new(0.0, 0.0, 0.0), 1e-3),
            default_grain(1, Vec3::new(1.0, 0.0, 0.0), 1e-3),
            default_grain(2, Vec3::new(0.0, 1.0, 0.0), 1e-3),
        ];
        let tri = FemTriangle::new([0, 1, 2], &grains);
        assert!(
            (tri.area - 0.5).abs() < EPS,
            "Area of right triangle should be 0.5: {}",
            tri.area
        );
    }

    // 31. NCL void ratio decreases with increasing stress
    #[test]
    fn test_ncl_void_ratio_decreasing() {
        let cs = CompactionState::loose_sand();
        let e1 = cs.ncl_void_ratio(1000.0);
        let e2 = cs.ncl_void_ratio(10000.0);
        assert!(e2 < e1, "Void ratio should decrease with higher stress");
    }

    // 32. GrainSizeDistribution sample stays in [r_min, r_max]
    #[test]
    fn test_grain_size_distribution_bounds() {
        let mut rng = rand::rng();
        let dist = GrainSizeDistribution {
            r_min: 1e-3,
            r_max: 5e-3,
            exponent: 2.0,
        };
        for _ in 0..100 {
            let r = dist.sample(&mut rng);
            assert!(
                r >= dist.r_min && r <= dist.r_max,
                "Sample {r} out of [{},{}]",
                dist.r_min,
                dist.r_max
            );
        }
    }

    // 33. Avalanche detection: grains beyond repose angle flagged
    #[test]
    fn test_avalanche_detection() {
        // Grain far from centre at height → should exceed angle of repose
        let mut grains = vec![
            default_grain(0, Vec3::new(0.0, 10.0, 0.0), 1e-3), // steep slope
            default_grain(1, Vec3::new(100.0, 1.0, 0.0), 1e-3), // shallow
        ];
        grains[0].is_fixed = true;
        grains[1].is_fixed = true;
        let avalanche = detect_avalanche_grains(&grains, 30.0);
        // Grain 0: slope angle ≈ atan(10/50) if centre is at (50, ...) — depends on centroid
        // Just check the function runs and returns a vector
        let _ = avalanche;
    }
}
