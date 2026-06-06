// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular media physics — packing, flow, avalanche, and hopper discharge.
//!
//! Provides discrete element method (DEM) style simulation of granular
//! assemblies including:
//!
//! - [`GranularParticle`] / [`GranularPacking`]: particle assembly with
//!   Hertz-Mindlin contact, normal/tangential restitution, and bulk statistics.
//! - [`AvalancheDynamics`]: angle-of-repose flow threshold model.
//! - [`HopperDischarge`]: Beverloo-law orifice discharge and jamming.
//! - [`ForceChain`]: contact force chain analysis.
//! - Free functions [`bagnold_shear_stress`] and [`janssen_pressure`].

use std::f64::consts::PI;

// ── internal vector helpers (no external crates) ─────────────────────────────

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let n = len3(v);
    if n > 1e-15 {
        scale3(v, 1.0 / n)
    } else {
        [0.0, 0.0, 0.0]
    }
}

// ── GranularParticle ──────────────────────────────────────────────────────────

/// A single grain in a granular assembly.
///
/// Tracks position, velocity, radius, mass, and the accumulated contact force
/// from the current simulation step.
#[derive(Debug, Clone)]
pub struct GranularParticle {
    /// World-space position \[m\].
    pub pos: [f64; 3],
    /// Velocity \[m/s\].
    pub vel: [f64; 3],
    /// Particle radius \[m\].
    pub radius: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Accumulated contact force for the current time step \[N\].
    pub contact_force: [f64; 3],
}

impl GranularParticle {
    /// Create a new particle at the given position with the given radius and mass.
    pub fn new(pos: [f64; 3], radius: f64, mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            radius,
            mass,
            contact_force: [0.0; 3],
        }
    }

    /// Compute kinetic energy: 0.5 * m * |v|^2.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.vel, self.vel)
    }

    /// Reset accumulated contact force to zero at the start of each step.
    fn reset_forces(&mut self) {
        self.contact_force = [0.0; 3];
    }
}

// ── GranularPacking ───────────────────────────────────────────────────────────

/// An assembly of granular particles governed by DEM contact mechanics.
///
/// Normal and tangential restitution coefficients control how much energy is
/// lost in collisions.  `e_n = 1` is perfectly elastic; `e_n = 0` is perfectly
/// inelastic.
#[derive(Debug, Clone)]
pub struct GranularPacking {
    /// Collection of particles in this assembly.
    pub particles: Vec<GranularParticle>,
    /// Gravitational acceleration vector \[m/s²\], e.g. `[0.0, -9.81, 0.0]`.
    pub gravity: [f64; 3],
    /// Normal coefficient of restitution (0 ≤ e_n ≤ 1).
    pub e_n: f64,
    /// Tangential coefficient of restitution (0 ≤ e_t ≤ 1).
    pub e_t: f64,
    /// Normal contact stiffness \[N/m\].
    pub k_n: f64,
    /// Tangential contact stiffness \[N/m\].
    pub k_t: f64,
}

impl GranularPacking {
    /// Create an empty packing with the given gravity and restitution coefficients.
    ///
    /// Default stiffnesses are `k_n = 1e6 N/m`, `k_t = 0.5 * k_n`.
    pub fn new(gravity: [f64; 3], e_n: f64, e_t: f64) -> Self {
        Self {
            particles: Vec::new(),
            gravity,
            e_n: e_n.clamp(0.0, 1.0),
            e_t: e_t.clamp(0.0, 1.0),
            k_n: 1.0e6,
            k_t: 5.0e5,
        }
    }

    /// Add a particle to the assembly and return its index.
    pub fn add_particle(&mut self, p: GranularParticle) -> usize {
        let idx = self.particles.len();
        self.particles.push(p);
        idx
    }

    /// Detect all overlapping particle pairs.
    ///
    /// Returns a list of `(i, j)` index pairs where particle `i` and `j`
    /// overlap (distance < sum of radii).
    pub fn detect_contacts(&self) -> Vec<(usize, usize)> {
        let n = self.particles.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = &self.particles[i];
                let pj = &self.particles[j];
                let d = len3(sub3(pj.pos, pi.pos));
                if d < pi.radius + pj.radius {
                    contacts.push((i, j));
                }
            }
        }
        contacts
    }

    /// Resolve all pairwise contact forces using a linear spring-dashpot model.
    ///
    /// Adds spring-based repulsive contact forces to each particle's
    /// `contact_force` accumulator.  Call this after `detect_contacts` (or
    /// let `step` do it for you).
    pub fn resolve_contacts(&mut self) {
        let n = self.particles.len();
        // Collect data first to avoid borrow issues.
        let mut forces = vec![[0.0f64; 3]; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let pi = &self.particles[i];
                let pj = &self.particles[j];
                let delta = sub3(pj.pos, pi.pos);
                let dist = len3(delta);
                let r_sum = pi.radius + pj.radius;
                if dist >= r_sum || dist < 1e-15 {
                    continue;
                }
                // Overlap depth.
                let overlap = r_sum - dist;
                let n_hat = normalize3(delta);

                // Normal contact force (Hertz-like linear spring).
                let f_n_mag = self.k_n * overlap;
                let f_n = scale3(n_hat, f_n_mag);

                // Relative normal velocity for damping.
                let v_rel = sub3(pj.vel, pi.vel);
                let v_n_mag = dot3(v_rel, n_hat);
                // Critical damping coefficient derived from e_n.
                let m_eff = pi.mass * pj.mass / (pi.mass + pj.mass);
                // ln(e_n) / sqrt(ln^2(e_n) + pi^2)
                let beta = if self.e_n > 1e-10 {
                    let ln_e = self.e_n.ln();
                    ln_e / (ln_e * ln_e + PI * PI).sqrt()
                } else {
                    -1.0
                };
                let gamma = -2.0 * beta * (self.k_n * m_eff).sqrt();
                let f_damp = scale3(n_hat, gamma * v_n_mag);

                // Tangential contact (simplified Coulomb, no friction history).
                let v_t = sub3(v_rel, scale3(n_hat, v_n_mag));
                let f_t_mag = self.k_t * overlap * len3(v_t).min(1.0);
                let f_t = if len3(v_t) > 1e-15 {
                    scale3(normalize3(v_t), -f_t_mag * (1.0 - self.e_t))
                } else {
                    [0.0; 3]
                };

                let f_total = add3(add3(f_n, f_damp), f_t);
                // Newton's 3rd law.
                forces[i] = add3(forces[i], f_total);
                forces[j] = add3(forces[j], scale3(f_total, -1.0));
            }
        }

        for (i, fi) in forces.iter().enumerate() {
            self.particles[i].contact_force = add3(self.particles[i].contact_force, *fi);
        }
    }

    /// Advance the simulation by one time step `dt` \[s\].
    ///
    /// Per-step pipeline:
    /// 1. Reset contact forces.
    /// 2. Resolve contacts.
    /// 3. Integrate velocity (semi-implicit Euler) with gravity + contact forces.
    /// 4. Integrate position.
    pub fn step(&mut self, dt: f64) {
        // Reset accumulated contact forces.
        for p in &mut self.particles {
            p.reset_forces();
        }

        self.resolve_contacts();

        // Integrate.
        for p in &mut self.particles {
            // Total acceleration = gravity + contact / mass.
            let a = add3(self.gravity, scale3(p.contact_force, 1.0 / p.mass));
            // Semi-implicit Euler: update velocity first, then position.
            p.vel = add3(p.vel, scale3(a, dt));
            p.pos = add3(p.pos, scale3(p.vel, dt));
        }
    }

    /// Compute the packing fraction (solid volume / bounding box volume).
    ///
    /// Uses an axis-aligned bounding box around all particle centres padded by
    /// the maximum radius.  Returns 0 if fewer than 2 particles are present.
    pub fn packing_fraction(&self) -> f64 {
        if self.particles.len() < 2 {
            return 0.0;
        }
        let mut lo = self.particles[0].pos;
        let mut hi = self.particles[0].pos;
        for p in &self.particles {
            lo[0] = lo[0].min(p.pos[0] - p.radius);
            lo[1] = lo[1].min(p.pos[1] - p.radius);
            lo[2] = lo[2].min(p.pos[2] - p.radius);
            hi[0] = hi[0].max(p.pos[0] + p.radius);
            hi[1] = hi[1].max(p.pos[1] + p.radius);
            hi[2] = hi[2].max(p.pos[2] + p.radius);
        }
        let vol_box = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
        if vol_box < 1e-15 {
            return 0.0;
        }
        let vol_solid: f64 = self
            .particles
            .iter()
            .map(|p| (4.0 / 3.0) * PI * p.radius.powi(3))
            .sum();
        (vol_solid / vol_box).min(1.0)
    }

    /// Compute the mean coordination number (average contacts per particle).
    pub fn coordination_number(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        let contacts = self.detect_contacts();
        // Each contact touches 2 particles.
        2.0 * contacts.len() as f64 / self.particles.len() as f64
    }

    /// Estimate the mean pressure from virial theorem: P = Σ f_ij · r_ij / (3 V).
    pub fn pressure(&self) -> f64 {
        let n = self.particles.len();
        if n < 2 {
            return 0.0;
        }
        let mut virial = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = &self.particles[i];
                let pj = &self.particles[j];
                let r_ij = sub3(pj.pos, pi.pos);
                let d = len3(r_ij);
                let r_sum = pi.radius + pj.radius;
                if d >= r_sum || d < 1e-15 {
                    continue;
                }
                let overlap = r_sum - d;
                let f_mag = self.k_n * overlap;
                virial += f_mag * d;
            }
        }
        // Bounding volume.
        let mut lo = self.particles[0].pos;
        let mut hi = self.particles[0].pos;
        for p in &self.particles {
            lo[0] = lo[0].min(p.pos[0] - p.radius);
            lo[1] = lo[1].min(p.pos[1] - p.radius);
            lo[2] = lo[2].min(p.pos[2] - p.radius);
            hi[0] = hi[0].max(p.pos[0] + p.radius);
            hi[1] = hi[1].max(p.pos[1] + p.radius);
            hi[2] = hi[2].max(p.pos[2] + p.radius);
        }
        let vol = ((hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])).max(1e-15);
        virial / (3.0 * vol)
    }

    /// Compute the total translational kinetic energy of all particles.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles.iter().map(|p| p.kinetic_energy()).sum()
    }

    /// Compute the granular temperature: T = (2/3) * `KE` / m.
    ///
    /// Returns zero for an empty packing.
    pub fn granular_temperature(&self) -> f64 {
        if self.particles.is_empty() {
            return 0.0;
        }
        let mean_ke = self.kinetic_energy() / self.particles.len() as f64;
        // Mean mass.
        let mean_mass =
            self.particles.iter().map(|p| p.mass).sum::<f64>() / self.particles.len() as f64;
        if mean_mass < 1e-15 {
            return 0.0;
        }
        (2.0 / 3.0) * mean_ke / mean_mass
    }
}

// ── AvalancheDynamics ─────────────────────────────────────────────────────────

/// Simple continuum model for avalanche initiation and run-out.
///
/// Uses angle-of-repose and flow threshold to determine whether a granular
/// slope is flowing, and Pouliquen–Forterre scaling for run-out and velocity.
#[derive(Debug, Clone)]
pub struct AvalancheDynamics {
    /// Static angle of repose \[radians\].
    pub angle_of_repose: f64,
    /// Dynamic flow threshold angle \[radians\] (slightly less than angle of repose).
    pub flow_threshold: f64,
    /// Current slope angle \[radians\].
    pub slope_angle: f64,
    /// Pile height \[m\].
    pub pile_height: f64,
    /// Gravitational acceleration magnitude \[m/s²\].
    pub gravity: f64,
}

impl AvalancheDynamics {
    /// Create a new avalanche model.
    ///
    /// * `angle_of_repose` — maximum stable slope angle \[rad\].
    /// * `flow_threshold` — angle at which flow is initiated \[rad\].
    /// * `slope_angle` — current slope angle \[rad\].
    /// * `pile_height` — height of the granular pile \[m\].
    /// * `gravity` — gravitational acceleration magnitude \[m/s²\].
    pub fn new(
        angle_of_repose: f64,
        flow_threshold: f64,
        slope_angle: f64,
        pile_height: f64,
        gravity: f64,
    ) -> Self {
        Self {
            angle_of_repose,
            flow_threshold,
            slope_angle,
            pile_height,
            gravity,
        }
    }

    /// Return `true` if the current slope angle exceeds the flow threshold.
    pub fn is_flowing(&self) -> bool {
        self.slope_angle > self.flow_threshold
    }

    /// Estimate the run-out distance using a simple energy-based model.
    ///
    /// Run-out ≈ H / tan(angle_of_repose) where H is the pile height.
    /// Returns 0 if not flowing.
    pub fn runout_distance(&self) -> f64 {
        if !self.is_flowing() {
            return 0.0;
        }
        let tan_phi = self.angle_of_repose.tan();
        if tan_phi < 1e-10 {
            return 0.0;
        }
        self.pile_height / tan_phi
    }

    /// Estimate the mean surface flow velocity \[m/s\] using Pouliquen scaling.
    ///
    /// v ≈ sqrt(g * H * sin(θ - θ_stop)) where θ_stop ≈ angle_of_repose.
    /// Returns 0 if not flowing.
    pub fn flow_velocity(&self) -> f64 {
        if !self.is_flowing() {
            return 0.0;
        }
        let excess = self.slope_angle - self.angle_of_repose;
        if excess <= 0.0 {
            return 0.0;
        }
        (self.gravity * self.pile_height * excess.sin())
            .max(0.0)
            .sqrt()
    }

    /// Froude number of the flow: Fr = v / sqrt(g * H * cos θ).
    ///
    /// Returns 0 if not flowing or if H is negligible.
    pub fn froude_number(&self) -> f64 {
        let v = self.flow_velocity();
        let denom = (self.gravity * self.pile_height * self.slope_angle.cos()).max(1e-15);
        v / denom.sqrt()
    }
}

// ── HopperDischarge ───────────────────────────────────────────────────────────

/// Hopper/silo orifice discharge model.
///
/// Implements the Beverloo equation for granular flow rate through a circular
/// orifice and estimates jamming probability from arch formation statistics.
#[derive(Debug, Clone)]
pub struct HopperDischarge {
    /// Orifice (outlet) radius \[m\].
    pub orifice_radius: f64,
    /// Mean grain radius \[m\].
    pub grain_radius: f64,
    /// Bulk grain density \[kg/m³\].
    pub bulk_density: f64,
    /// Gravitational acceleration magnitude \[m/s²\].
    pub gravity: f64,
    /// Beverloo discharge coefficient C (typically 0.55–0.65).
    pub beverloo_c: f64,
    /// Beverloo shape factor k (typically 1.4–1.8).
    pub beverloo_k: f64,
}

impl HopperDischarge {
    /// Create a new hopper discharge model with default Beverloo coefficients.
    ///
    /// * `orifice_radius` — radius of the outlet \[m\].
    /// * `grain_radius` — mean grain radius \[m\].
    /// * `bulk_density` — bulk density of the granular material \[kg/m³\].
    /// * `gravity` — gravitational acceleration magnitude \[m/s²\].
    pub fn new(orifice_radius: f64, grain_radius: f64, bulk_density: f64, gravity: f64) -> Self {
        Self {
            orifice_radius,
            grain_radius,
            bulk_density,
            gravity,
            beverloo_c: 0.58,
            beverloo_k: 1.5,
        }
    }

    /// Mass flow rate via the Beverloo equation \[kg/s\].
    ///
    /// Q = C * ρ * √g * (D_o - k * d_g)^(5/2)
    ///
    /// where D_o = orifice diameter, d_g = grain diameter.
    /// Returns 0 if the effective orifice is blocked (Beverloo blockage condition).
    pub fn beverloo_discharge_rate(&self) -> f64 {
        let d_o = 2.0 * self.orifice_radius;
        let d_g = 2.0 * self.grain_radius;
        let effective = d_o - self.beverloo_k * d_g;
        if effective <= 0.0 {
            return 0.0;
        }
        self.beverloo_c * self.bulk_density * self.gravity.sqrt() * effective.powf(2.5)
    }

    /// Probability of a stable arch forming and jamming the orifice.
    ///
    /// Uses the empirical relation P_jam ≈ exp(-α (D_o/d_g - 1)) where α ≈ 1.
    /// P_jam → 1 as D_o/d_g → 1 (single grain width = guaranteed jam).
    pub fn jam_probability(&self) -> f64 {
        if self.grain_radius < 1e-15 {
            return 0.0;
        }
        let ratio = self.orifice_radius / self.grain_radius;
        if ratio <= 1.0 {
            return 1.0;
        }
        let alpha = 1.0_f64;
        (-alpha * (ratio - 1.0)).exp()
    }

    /// Velocity of grains exiting the orifice (Torricelli-like) \[m/s\].
    ///
    /// v_exit = sqrt(2 * g * H) where H is approximated as orifice radius
    /// (free-fall height scale).
    pub fn exit_velocity(&self) -> f64 {
        (2.0 * self.gravity * self.orifice_radius).sqrt()
    }

    /// Number flow rate \[particles/s\] = mass rate / (mass per grain).
    pub fn particle_flow_rate(&self) -> f64 {
        let grain_volume = (4.0 / 3.0) * PI * self.grain_radius.powi(3);
        let grain_mass = self.bulk_density * grain_volume;
        if grain_mass < 1e-20 {
            return 0.0;
        }
        self.beverloo_discharge_rate() / grain_mass
    }
}

// ── ForceChain ────────────────────────────────────────────────────────────────

/// Representation of a force chain in a granular assembly.
///
/// A force chain is a quasi-linear sequence of heavily loaded contacts that
/// transmits stress through the assembly.  Each entry is `(i, j, f)` where
/// `i` and `j` are particle indices and `f` is the contact force magnitude \[N\].
#[derive(Debug, Clone)]
pub struct ForceChain {
    /// List of contacts: `(particle_i, particle_j, force_magnitude [N])`.
    pub contacts: Vec<(usize, usize, f64)>,
}

impl ForceChain {
    /// Create a new force chain from a list of contacts.
    pub fn new(contacts: Vec<(usize, usize, f64)>) -> Self {
        Self { contacts }
    }

    /// Build a force chain by extracting the strongest contacts from a packing.
    ///
    /// Returns all contacts whose force magnitude exceeds `threshold` times the
    /// mean contact force.
    pub fn from_packing(packing: &GranularPacking, threshold: f64) -> Self {
        let n = packing.particles.len();
        let mut raw: Vec<(usize, usize, f64)> = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = &packing.particles[i];
                let pj = &packing.particles[j];
                let d = len3(sub3(pj.pos, pi.pos));
                let r_sum = pi.radius + pj.radius;
                if d < r_sum && d > 1e-15 {
                    let overlap = r_sum - d;
                    let f = packing.k_n * overlap;
                    raw.push((i, j, f));
                }
            }
        }
        let mean_f = if raw.is_empty() {
            0.0
        } else {
            raw.iter().map(|x| x.2).sum::<f64>() / raw.len() as f64
        };
        let chain_contacts: Vec<_> = raw
            .into_iter()
            .filter(|x| x.2 >= threshold * mean_f)
            .collect();
        Self {
            contacts: chain_contacts,
        }
    }

    /// Return the maximum contact force magnitude in the chain \[N\].
    pub fn max_force(&self) -> f64 {
        self.contacts.iter().map(|c| c.2).fold(0.0_f64, f64::max)
    }

    /// Return the number of contacts (links) in the force chain.
    pub fn chain_length(&self) -> usize {
        self.contacts.len()
    }

    /// Compute the force anisotropy index: (f_max - f_min) / (f_max + f_min).
    ///
    /// Returns 0 for an empty or single-contact chain.
    pub fn anisotropy(&self) -> f64 {
        if self.contacts.len() < 2 {
            return 0.0;
        }
        let f_max = self.contacts.iter().map(|c| c.2).fold(0.0_f64, f64::max);
        let f_min = self
            .contacts
            .iter()
            .map(|c| c.2)
            .fold(f64::INFINITY, f64::min);
        let denom = f_max + f_min;
        if denom < 1e-15 {
            return 0.0;
        }
        (f_max - f_min) / denom
    }

    /// Mean force in the chain \[N\].
    pub fn mean_force(&self) -> f64 {
        if self.contacts.is_empty() {
            return 0.0;
        }
        self.contacts.iter().map(|c| c.2).sum::<f64>() / self.contacts.len() as f64
    }
}

// ── Free functions ────────────────────────────────────────────────────────────

/// Bagnold shear stress in a rapidly shearing granular flow \[Pa\].
///
/// τ = ρ * d² * γ̇² * a_i(C)
///
/// where the Bagnold constant factor a_i(C) ≈ 0.042 is absorbed into the
/// leading coefficient.  This is a simplified grain-inertia regime formula:
///
/// τ ≈ (1/2) * ρ * d² * γ̇²
///
/// * `density` — grain density \[kg/m³\].
/// * `grain_r` — grain radius \[m\].
/// * `shear_rate` — shear rate γ̇ \[1/s\].
pub fn bagnold_shear_stress(density: f64, grain_r: f64, shear_rate: f64) -> f64 {
    let d = 2.0 * grain_r;
    // Bagnold 1954 grain-inertia coefficient ≈ 0.042; simplified as 1/2 here.
    0.5 * density * d * d * shear_rate * shear_rate
}

/// Janssen pressure in a granular silo \[Pa\].
///
/// The Janssen model accounts for wall friction reducing the vertical stress
/// below the hydrostatic value:
///
/// σ_v(z) = (ρ g R) / (2 μ_w k) * \[1 - exp(-2 μ_w k z / R)\]
///
/// * `depth` — depth below the free surface \[m\].
/// * `density` — bulk grain density \[kg/m³\].
/// * `k` — lateral pressure ratio (typically 0.3–0.6).
/// * `radius` — silo radius \[m\].
/// * `mu_w` — wall friction coefficient.
pub fn janssen_pressure(depth: f64, density: f64, k: f64, radius: f64, mu_w: f64) -> f64 {
    let g = 9.81;
    let denom = 2.0 * mu_w * k;
    if denom < 1e-15 || radius < 1e-15 {
        // Degenerate: fall back to hydrostatic.
        return density * g * depth;
    }
    let scale = density * g * radius / denom;
    let exp_term = (-denom * depth / radius).exp();
    scale * (1.0 - exp_term)
}

/// Compute the angle of repose from measured loose and dense packing fractions.
///
/// Empirical relation: φ_r ≈ arctan(0.32 + 0.43 * φ_dense).
/// Returns angle in radians.
pub fn empirical_angle_of_repose(phi_dense: f64) -> f64 {
    (0.32_f64 + 0.43 * phi_dense).atan()
}

/// Compute the critical velocity for a grain to remain stable on a slope.
///
/// v_c = sqrt(g * d * cos(θ) * (tan(φ) - tan(θ)))
///
/// * `grain_diameter` — grain diameter \[m\].
/// * `slope_angle` — slope angle \[rad\].
/// * `angle_of_repose` — angle of repose \[rad\].
/// * `gravity` — gravitational acceleration \[m/s²\].
pub fn critical_slope_velocity(
    grain_diameter: f64,
    slope_angle: f64,
    angle_of_repose: f64,
    gravity: f64,
) -> f64 {
    let stability = slope_angle.tan() - angle_of_repose.tan();
    if stability >= 0.0 {
        // Slope exceeds repose: no stable velocity.
        return 0.0;
    }
    (gravity * grain_diameter * slope_angle.cos() * (-stability))
        .max(0.0)
        .sqrt()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Helper: create a packing with two touching particles.
    fn two_particle_packing() -> GranularPacking {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.5, 0.3);
        // Place particles just overlapping.
        let mut p1 = GranularParticle::new([0.0, 0.0, 0.0], 0.5, 1.0);
        p1.vel = [1.0, 0.0, 0.0];
        let p2 = GranularParticle::new([0.8, 0.0, 0.0], 0.5, 1.0);
        pack.add_particle(p1);
        pack.add_particle(p2);
        pack
    }

    // 1. GranularParticle: kinetic energy is correct.
    #[test]
    fn test_particle_kinetic_energy() {
        let mut p = GranularParticle::new([0.0; 3], 0.01, 2.0);
        p.vel = [3.0, 4.0, 0.0]; // speed = 5
        let ke = p.kinetic_energy();
        // KE = 0.5 * 2 * 25 = 25
        assert!((ke - 25.0).abs() < 1e-10, "KE={ke}");
    }

    // 2. GranularPacking: add_particle returns correct index.
    #[test]
    fn test_add_particle_index() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        let i0 = pack.add_particle(GranularParticle::new([0.0; 3], 0.1, 1.0));
        let i1 = pack.add_particle(GranularParticle::new([1.0, 0.0, 0.0], 0.1, 1.0));
        assert_eq!(i0, 0);
        assert_eq!(i1, 1);
    }

    // 3. detect_contacts: overlapping pair detected.
    #[test]
    fn test_detect_contacts_overlap() {
        let pack = two_particle_packing();
        let contacts = pack.detect_contacts();
        assert_eq!(contacts.len(), 1, "one overlapping pair expected");
        assert_eq!(contacts[0], (0, 1));
    }

    // 4. detect_contacts: non-overlapping particles produce no contacts.
    #[test]
    fn test_detect_contacts_no_overlap() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        pack.add_particle(GranularParticle::new([0.0; 3], 0.1, 1.0));
        pack.add_particle(GranularParticle::new([10.0, 0.0, 0.0], 0.1, 1.0));
        assert!(pack.detect_contacts().is_empty());
    }

    // 5. resolve_contacts: contact force is non-zero for overlapping particles.
    #[test]
    fn test_resolve_contacts_nonzero_force() {
        let mut pack = two_particle_packing();
        pack.resolve_contacts();
        let f0 = pack.particles[0].contact_force;
        let mag = len3(f0);
        assert!(mag > 0.0, "contact force should be nonzero, got {mag}");
    }

    // 6. resolve_contacts: Newton's 3rd law — forces sum to zero.
    #[test]
    fn test_contact_force_reaction() {
        let mut pack = two_particle_packing();
        pack.resolve_contacts();
        let f0 = pack.particles[0].contact_force;
        let f1 = pack.particles[1].contact_force;
        let sum = add3(f0, f1);
        let mag = len3(sum);
        assert!(
            mag < 1e-6,
            "force sum should be zero (Newton 3rd), |sum|={mag}"
        );
    }

    // 7. step: kinetic energy changes after one step with gravity.
    #[test]
    fn test_step_changes_kinetic_energy() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.5, 0.3);
        pack.add_particle(GranularParticle::new([0.0, 10.0, 0.0], 0.1, 1.0));
        let ke_before = pack.kinetic_energy();
        pack.step(0.01);
        let ke_after = pack.kinetic_energy();
        assert!(
            ke_after > ke_before,
            "gravity should increase KE, before={ke_before}, after={ke_after}"
        );
    }

    // 8. packing_fraction: single particle returns 0.
    #[test]
    fn test_packing_fraction_single_particle() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        pack.add_particle(GranularParticle::new([0.0; 3], 0.1, 1.0));
        assert_eq!(pack.packing_fraction(), 0.0);
    }

    // 9. packing_fraction: two particles, value in [0, 1].
    #[test]
    fn test_packing_fraction_range() {
        let pack = two_particle_packing();
        let phi = pack.packing_fraction();
        assert!(
            (0.0..=1.0).contains(&phi),
            "packing fraction out of range: {phi}"
        );
    }

    // 10. coordination_number: empty packing is 0.
    #[test]
    fn test_coordination_number_empty() {
        let pack = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        assert_eq!(pack.coordination_number(), 0.0);
    }

    // 11. coordination_number: two overlapping particles → 1 contact → Z = 1.
    #[test]
    fn test_coordination_number_two_particles() {
        let pack = two_particle_packing();
        let z = pack.coordination_number();
        // 2 contacts counted as 1 pair → Z = 2*1/2 = 1.
        assert!((z - 1.0).abs() < 1e-10, "Z={z}");
    }

    // 12. pressure: increases with more overlap.
    #[test]
    fn test_pressure_increases_with_overlap() {
        let mut pack1 = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        pack1.add_particle(GranularParticle::new([0.0; 3], 0.5, 1.0));
        pack1.add_particle(GranularParticle::new([0.8, 0.0, 0.0], 0.5, 1.0));

        let mut pack2 = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        pack2.add_particle(GranularParticle::new([0.0; 3], 0.5, 1.0));
        pack2.add_particle(GranularParticle::new([0.5, 0.0, 0.0], 0.5, 1.0)); // more overlap

        assert!(
            pack2.pressure() > pack1.pressure(),
            "more overlap should give higher pressure"
        );
    }

    // 13. kinetic_energy: sums over all particles.
    #[test]
    fn test_packing_kinetic_energy_sum() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        let mut p = GranularParticle::new([0.0; 3], 0.1, 2.0);
        p.vel = [1.0, 0.0, 0.0];
        pack.add_particle(p);
        let mut p2 = GranularParticle::new([1.0, 0.0, 0.0], 0.1, 2.0);
        p2.vel = [0.0, 1.0, 0.0];
        pack.add_particle(p2);
        let ke = pack.kinetic_energy();
        // 0.5*2*1 + 0.5*2*1 = 2
        assert!((ke - 2.0).abs() < 1e-10, "KE={ke}");
    }

    // 14. AvalancheDynamics: is_flowing respects threshold.
    #[test]
    fn test_avalanche_is_flowing() {
        let phi = 30_f64.to_radians();
        let thresh = 28_f64.to_radians();
        let av = AvalancheDynamics::new(phi, thresh, 35_f64.to_radians(), 2.0, 9.81);
        assert!(av.is_flowing());

        let av2 = AvalancheDynamics::new(phi, thresh, 20_f64.to_radians(), 2.0, 9.81);
        assert!(!av2.is_flowing());
    }

    // 15. AvalancheDynamics: runout_distance is 0 if not flowing.
    #[test]
    fn test_runout_not_flowing() {
        let phi = 30_f64.to_radians();
        let thresh = 28_f64.to_radians();
        let av = AvalancheDynamics::new(phi, thresh, 20_f64.to_radians(), 2.0, 9.81);
        assert_eq!(av.runout_distance(), 0.0);
    }

    // 16. AvalancheDynamics: runout_distance positive when flowing.
    #[test]
    fn test_runout_positive_when_flowing() {
        let phi = 30_f64.to_radians();
        let thresh = 28_f64.to_radians();
        let av = AvalancheDynamics::new(phi, thresh, 35_f64.to_radians(), 2.0, 9.81);
        assert!(av.runout_distance() > 0.0);
    }

    // 17. AvalancheDynamics: flow_velocity is 0 if not flowing.
    #[test]
    fn test_flow_velocity_not_flowing() {
        let phi = 30_f64.to_radians();
        let thresh = 28_f64.to_radians();
        let av = AvalancheDynamics::new(phi, thresh, 15_f64.to_radians(), 2.0, 9.81);
        assert_eq!(av.flow_velocity(), 0.0);
    }

    // 18. HopperDischarge: beverloo rate is zero when blocked.
    #[test]
    fn test_beverloo_blocked() {
        // orifice_radius = 0.5 * k * d_g → blocked.
        let hd = HopperDischarge::new(0.001, 0.001, 1500.0, 9.81);
        // effective = 2*0.001 - 1.5*2*0.001 = 0.002 - 0.003 = -0.001 → blocked
        assert_eq!(hd.beverloo_discharge_rate(), 0.0);
    }

    // 19. HopperDischarge: larger orifice gives higher flow rate.
    #[test]
    fn test_beverloo_larger_orifice() {
        let hd1 = HopperDischarge::new(0.05, 0.005, 1500.0, 9.81);
        let hd2 = HopperDischarge::new(0.10, 0.005, 1500.0, 9.81);
        assert!(
            hd2.beverloo_discharge_rate() > hd1.beverloo_discharge_rate(),
            "larger orifice should discharge faster"
        );
    }

    // 20. HopperDischarge: jam_probability → 1 for single grain width.
    #[test]
    fn test_jam_probability_single_grain() {
        let hd = HopperDischarge::new(0.01, 0.01, 1500.0, 9.81);
        // ratio = 1.0 → jam probability = 1.0
        assert!((hd.jam_probability() - 1.0).abs() < 1e-10);
    }

    // 21. HopperDischarge: jam_probability decreases with larger orifice.
    #[test]
    fn test_jam_probability_decreasing() {
        let hd1 = HopperDischarge::new(0.02, 0.005, 1500.0, 9.81);
        let hd2 = HopperDischarge::new(0.10, 0.005, 1500.0, 9.81);
        assert!(hd2.jam_probability() < hd1.jam_probability());
    }

    // 22. ForceChain: max_force is correct.
    #[test]
    fn test_force_chain_max_force() {
        let fc = ForceChain::new(vec![(0, 1, 10.0), (1, 2, 50.0), (2, 3, 30.0)]);
        assert!((fc.max_force() - 50.0).abs() < 1e-10);
    }

    // 23. ForceChain: chain_length is correct.
    #[test]
    fn test_force_chain_length() {
        let fc = ForceChain::new(vec![(0, 1, 5.0), (1, 2, 8.0)]);
        assert_eq!(fc.chain_length(), 2);
    }

    // 24. ForceChain: anisotropy is in [0, 1].
    #[test]
    fn test_force_chain_anisotropy_range() {
        let fc = ForceChain::new(vec![(0, 1, 2.0), (1, 2, 8.0), (2, 3, 5.0)]);
        let a = fc.anisotropy();
        assert!((0.0..=1.0).contains(&a), "anisotropy={a}");
    }

    // 25. bagnold_shear_stress: scales as shear_rate squared.
    #[test]
    fn test_bagnold_shear_stress_scaling() {
        let s1 = bagnold_shear_stress(2500.0, 0.001, 10.0);
        let s2 = bagnold_shear_stress(2500.0, 0.001, 20.0);
        // Doubling shear rate → 4x stress.
        assert!((s2 / s1 - 4.0).abs() < 1e-8, "ratio={}", s2 / s1);
    }

    // 26. janssen_pressure: saturates at large depth (Janssen effect).
    #[test]
    fn test_janssen_pressure_saturates() {
        let p_shallow = janssen_pressure(0.1, 1500.0, 0.4, 0.5, 0.3);
        let p_deep = janssen_pressure(100.0, 1500.0, 0.4, 0.5, 0.3);
        let p_very_deep = janssen_pressure(10000.0, 1500.0, 0.4, 0.5, 0.3);
        // Pressure should grow then saturate.
        assert!(p_deep > p_shallow, "pressure should increase with depth");
        assert!(
            (p_very_deep - p_deep).abs() / p_very_deep < 0.01,
            "pressure should saturate at large depth"
        );
    }

    // 27. janssen_pressure: small depth approaches hydrostatic.
    #[test]
    fn test_janssen_pressure_shallow_hydrostatic() {
        let depth = 0.001;
        let density = 1500.0;
        let jp = janssen_pressure(depth, density, 0.4, 1.0, 0.3);
        let hydrostatic = density * 9.81 * depth;
        let rel_err = (jp - hydrostatic).abs() / hydrostatic;
        assert!(
            rel_err < 0.01,
            "shallow Janssen should match hydrostatic, rel_err={rel_err}"
        );
    }

    // 28. granular_temperature: zero for stationary pack.
    #[test]
    fn test_granular_temperature_zero_stationary() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.9, 0.9);
        pack.add_particle(GranularParticle::new([0.0; 3], 0.1, 1.0));
        pack.add_particle(GranularParticle::new([1.0, 0.0, 0.0], 0.1, 1.0));
        assert_eq!(pack.granular_temperature(), 0.0);
    }

    // 29. ForceChain::from_packing: high threshold → fewer contacts.
    #[test]
    fn test_force_chain_from_packing_threshold() {
        let pack = two_particle_packing();
        let fc_lo = ForceChain::from_packing(&pack, 0.0);
        let fc_hi = ForceChain::from_packing(&pack, 100.0);
        assert!(
            fc_lo.chain_length() >= fc_hi.chain_length(),
            "higher threshold should yield fewer or equal contacts"
        );
    }

    // 30. empirical_angle_of_repose: output in (0, π/2).
    #[test]
    fn test_empirical_angle_of_repose_range() {
        let phi = empirical_angle_of_repose(0.64);
        assert!(
            phi > 0.0 && phi < std::f64::consts::FRAC_PI_2,
            "angle out of range: {phi}"
        );
    }

    // 31. critical_slope_velocity: zero when slope exceeds repose.
    #[test]
    fn test_critical_slope_velocity_unstable() {
        let v = critical_slope_velocity(0.01, 40_f64.to_radians(), 30_f64.to_radians(), 9.81);
        assert_eq!(v, 0.0, "unstable slope should give v=0");
    }

    // 32. ForceChain: anisotropy is 0 for uniform forces.
    #[test]
    fn test_force_chain_anisotropy_uniform() {
        let fc = ForceChain::new(vec![(0, 1, 5.0), (1, 2, 5.0), (2, 3, 5.0)]);
        assert!(
            fc.anisotropy().abs() < 1e-10,
            "uniform forces → anisotropy=0"
        );
    }

    // 33. HopperDischarge: particle_flow_rate is positive for open orifice.
    #[test]
    fn test_particle_flow_rate_positive() {
        let hd = HopperDischarge::new(0.05, 0.003, 1500.0, 9.81);
        assert!(hd.particle_flow_rate() > 0.0);
    }

    // 34. bagnold_shear_stress: scales as grain diameter squared.
    #[test]
    fn test_bagnold_shear_stress_grain_size_scaling() {
        let s1 = bagnold_shear_stress(2500.0, 0.001, 10.0);
        let s2 = bagnold_shear_stress(2500.0, 0.002, 10.0);
        // d doubles → 4x stress.
        assert!((s2 / s1 - 4.0).abs() < 1e-8, "ratio={}", s2 / s1);
    }

    // 35. step: multiple steps integrate without NaN.
    #[test]
    fn test_step_no_nan() {
        let mut pack = GranularPacking::new([0.0, -9.81, 0.0], 0.5, 0.5);
        pack.add_particle(GranularParticle::new([0.0, 1.0, 0.0], 0.1, 1.0));
        pack.add_particle(GranularParticle::new([0.15, 1.0, 0.0], 0.1, 1.0));
        for _ in 0..100 {
            pack.step(0.001);
        }
        for p in &pack.particles {
            assert!(
                !p.pos[0].is_nan() && !p.pos[1].is_nan() && !p.pos[2].is_nan(),
                "NaN in position"
            );
            assert!(
                !p.vel[0].is_nan() && !p.vel[1].is_nan() && !p.vel[2].is_nan(),
                "NaN in velocity"
            );
        }
    }
}
