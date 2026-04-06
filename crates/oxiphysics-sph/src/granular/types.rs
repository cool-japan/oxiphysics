//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
#[allow(unused_imports)]
use super::functions_2::*;
use std::f64::consts::PI;

#[allow(unused_imports)]
use super::functions::*;

/// Extended granular simulation that also tracks angular velocities.
#[allow(dead_code)]
pub struct GranularSimExtended {
    /// Underlying DEM simulator.
    pub sim: GranularSim,
    /// Angular velocities for each particle.
    pub angular_velocities: Vec<[f64; 3]>,
    /// Torques accumulated this step.
    pub torques: Vec<[f64; 3]>,
    /// Rolling resistance model.
    pub rolling_model: RollingResistanceModel,
    /// Rolling friction coefficient.
    pub mu_r: f64,
}
#[allow(dead_code)]
impl GranularSimExtended {
    /// Create a new extended simulation.
    pub fn new(
        contact_model: ContactModel,
        mu: f64,
        gravity: [f64; 3],
        rolling_model: RollingResistanceModel,
        mu_r: f64,
    ) -> Self {
        Self {
            sim: GranularSim::new(contact_model, mu, gravity),
            angular_velocities: Vec::new(),
            torques: Vec::new(),
            rolling_model,
            mu_r,
        }
    }
    /// Add a particle.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], radius: f64, mass: f64) {
        self.sim.add_particle(pos, vel, radius, mass);
        self.angular_velocities.push([0.0; 3]);
        self.torques.push([0.0; 3]);
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.sim.len()
    }
    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.sim.is_empty()
    }
    /// Advance one step with rolling resistance.
    pub fn step(&mut self, dt: f64) {
        let contacts = self.sim.detect_contacts();
        self.sim.resolve_contacts(&contacts);
        let n = self.len();
        for t in self.torques.iter_mut().take(n) {
            *t = [0.0; 3];
        }
        for c in &contacts {
            let r_i = self.sim.radii[c.i];
            let r_j = self.sim.radii[c.j];
            let fn_mag = self.sim.forces[c.i][0] * c.normal[0]
                + self.sim.forces[c.i][1] * c.normal[1]
                + self.sim.forces[c.i][2] * c.normal[2];
            let omega_rel = [
                self.angular_velocities[c.i][0] - self.angular_velocities[c.j][0],
                self.angular_velocities[c.i][1] - self.angular_velocities[c.j][1],
                self.angular_velocities[c.i][2] - self.angular_velocities[c.j][2],
            ];
            let torque = rolling_resistance_torque(
                r_i,
                r_j,
                fn_mag.abs(),
                omega_rel,
                self.mu_r,
                self.rolling_model,
            );
            for k in 0..3 {
                self.torques[c.i][k] += torque[k];
                self.torques[c.j][k] -= torque[k];
            }
        }
        for i in 0..n {
            let inv_m = 1.0 / self.sim.masses[i].max(1e-30);
            for k in 0..3 {
                self.sim.velocities[i][k] += self.sim.forces[i][k] * inv_m * dt;
                self.sim.positions[i][k] += self.sim.velocities[i][k] * dt;
            }
            let r = self.sim.radii[i];
            let inertia = 0.4 * self.sim.masses[i] * r * r;
            let inv_i = 1.0 / inertia.max(1e-30);
            for k in 0..3 {
                self.angular_velocities[i][k] += self.torques[i][k] * inv_i * dt;
            }
        }
    }
    /// Compute total rotational kinetic energy.
    pub fn rotational_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.len() {
            let r = self.sim.radii[i];
            let inertia = 0.4 * self.sim.masses[i] * r * r;
            let w2 = self.angular_velocities[i][0].powi(2)
                + self.angular_velocities[i][1].powi(2)
                + self.angular_velocities[i][2].powi(2);
            ke += 0.5 * inertia * w2;
        }
        ke
    }
}
/// A detected contact between two particles.
#[derive(Debug, Clone)]
pub struct Contact {
    /// Index of the first particle.
    pub i: usize,
    /// Index of the second particle.
    pub j: usize,
    /// Normal overlap (positive = overlapping).
    pub overlap: f64,
    /// Contact normal (from j to i, unit vector).
    pub normal: [f64; 3],
}
/// Forces and geometry produced by a DEM contact calculation.
#[derive(Debug, Clone)]
pub struct DemContactForce {
    /// Normal force vector (points from j to i, i.e. repulsive on i) \[N\].
    pub normal_force: [f64; 3],
    /// Tangential friction force vector (acts on particle i) \[N\].
    pub tangential_force: [f64; 3],
    /// Contact point (midpoint of overlap region) \[m\].
    pub contact_point: [f64; 3],
}
/// A discrete-element (DEM) granular particle with position, velocity, spin,
/// radius, and mass.
pub struct DemGranularParticle {
    /// Position in 3-D space.
    pub position: [f64; 3],
    /// Translational velocity.
    pub velocity: [f64; 3],
    /// Angular velocity (spin).
    pub angular_velocity: [f64; 3],
    /// Particle radius \[m\].
    pub radius: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
}
impl DemGranularParticle {
    /// Create a new spherical granular particle.
    pub fn new(position: [f64; 3], radius: f64, density: f64) -> Self {
        let mass = (4.0 / 3.0) * PI * radius.powi(3) * density;
        Self {
            position,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            radius,
            mass,
        }
    }
    /// Moment of inertia for a solid sphere: `I = 2/5 m r²`.
    pub fn moment_of_inertia(&self) -> f64 {
        0.4 * self.mass * self.radius * self.radius
    }
    /// Kinetic energy (translational + rotational).
    pub fn kinetic_energy(&self) -> f64 {
        let v2 = self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2);
        let w2 = self.angular_velocity[0].powi(2)
            + self.angular_velocity[1].powi(2)
            + self.angular_velocity[2].powi(2);
        0.5 * self.mass * v2 + 0.5 * self.moment_of_inertia() * w2
    }
}
/// Simple granular (DEM) simulation with contact detection and resolution.
pub struct GranularSim {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle radii.
    pub radii: Vec<f64>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Forces accumulated this step.
    pub forces: Vec<[f64; 3]>,
    /// Contact model.
    pub contact_model: ContactModel,
    /// Coulomb friction coefficient.
    pub mu: f64,
    /// Gravity vector.
    pub gravity: [f64; 3],
}
impl GranularSim {
    /// Create a new simulation with the given contact model.
    pub fn new(contact_model: ContactModel, mu: f64, gravity: [f64; 3]) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            radii: Vec::new(),
            masses: Vec::new(),
            forces: Vec::new(),
            contact_model,
            mu,
            gravity,
        }
    }
    /// Add a particle.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], radius: f64, mass: f64) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.radii.push(radius);
        self.masses.push(mass);
        self.forces.push([0.0; 3]);
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Whether empty.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Detect all pairwise contacts. O(n²) brute force.
    pub fn detect_contacts(&self) -> Vec<Contact> {
        let n = self.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.positions[i][0] - self.positions[j][0];
                let dy = self.positions[i][1] - self.positions[j][1];
                let dz = self.positions[i][2] - self.positions[j][2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let overlap = self.radii[i] + self.radii[j] - dist;
                if overlap > 0.0 && dist > 1e-14 {
                    let inv = 1.0 / dist;
                    contacts.push(Contact {
                        i,
                        j,
                        overlap,
                        normal: [dx * inv, dy * inv, dz * inv],
                    });
                }
            }
        }
        contacts
    }
    /// Resolve detected contacts by computing and applying contact forces.
    pub fn resolve_contacts(&mut self, contacts: &[Contact]) {
        for i in 0..self.len() {
            self.forces[i] = [
                self.masses[i] * self.gravity[0],
                self.masses[i] * self.gravity[1],
                self.masses[i] * self.gravity[2],
            ];
        }
        for c in contacts {
            let rel_vel = [
                self.velocities[c.i][0] - self.velocities[c.j][0],
                self.velocities[c.i][1] - self.velocities[c.j][1],
                self.velocities[c.i][2] - self.velocities[c.j][2],
            ];
            let vn = rel_vel[0] * c.normal[0] + rel_vel[1] * c.normal[1] + rel_vel[2] * c.normal[2];
            let fn_mag = match &self.contact_model {
                ContactModel::SpringDashpot { kn, gamma_n, .. } => {
                    spring_dashpot_normal(c.overlap, vn, *kn, *gamma_n)
                }
                ContactModel::Hertzian { e_star, r_star } => {
                    hertz_normal_force(c.overlap, *e_star, *r_star)
                }
            };
            for k in 0..3 {
                let f = fn_mag * c.normal[k];
                self.forces[c.i][k] += f;
                self.forces[c.j][k] -= f;
            }
        }
    }
    /// Advance one step: detect contacts, resolve, integrate.
    pub fn step(&mut self, dt: f64) {
        let contacts = self.detect_contacts();
        self.resolve_contacts(&contacts);
        for i in 0..self.len() {
            let inv_m = 1.0 / self.masses[i].max(1e-30);
            for k in 0..3 {
                self.velocities[i][k] += self.forces[i][k] * inv_m * dt;
                self.positions[i][k] += self.velocities[i][k] * dt;
            }
        }
    }
}
/// Pairwise contact data between two DEM spheres.
pub struct DemContact {
    /// Index of particle i.
    pub i: usize,
    /// Index of particle j.
    pub j: usize,
    /// Normal overlap δ_n \[m\] (positive = interpenetration).
    pub overlap: f64,
    /// Unit contact normal pointing from i toward j.
    pub normal: [f64; 3],
    /// Effective Young's modulus E* \[Pa\].
    pub e_eff: f64,
    /// Effective radius R* \[m\].
    pub r_eff: f64,
    /// Coefficient of restitution.
    pub restitution: f64,
    /// Rolling friction coefficient μ_r \[-\].
    pub mu_rolling: f64,
}
impl DemContact {
    /// Construct a new `DemContact`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        i: usize,
        j: usize,
        overlap: f64,
        normal: [f64; 3],
        e_eff: f64,
        r_eff: f64,
        restitution: f64,
        mu_rolling: f64,
    ) -> Self {
        Self {
            i,
            j,
            overlap,
            normal,
            e_eff,
            r_eff,
            restitution,
            mu_rolling,
        }
    }
    /// Hertz elastic normal contact force vector acting on particle i.
    ///
    /// Uses the standard Hertz model:
    ///   `F_n = (4/3) E* sqrt(R*) δ_n^{3/2}`
    ///
    /// The force vector acts along `self.normal` (repulsive on i, so opposite
    /// to the normal direction from i toward j).
    ///
    /// Returns `[0; 3]` if overlap ≤ 0.
    pub fn compute_hertz_contact_force(&self) -> [f64; 3] {
        if self.overlap <= 0.0 {
            return [0.0; 3];
        }
        let fn_mag = (4.0 / 3.0) * self.e_eff * self.r_eff.sqrt() * self.overlap.powf(1.5);
        [
            -fn_mag * self.normal[0],
            -fn_mag * self.normal[1],
            -fn_mag * self.normal[2],
        ]
    }
    /// Rolling resistance torque on particle i.
    ///
    /// Uses the constant rolling-resistance model (Ai et al. 2011 Model A):
    ///   `τ_r = -μ_r * R_eff * |F_n| * ω̂_rel`
    ///
    /// where `ω̂_rel` is the unit vector of relative angular velocity
    /// `ω_i - ω_j`.  If `|ω_rel| ≈ 0` the torque is zero.
    ///
    /// # Arguments
    /// * `omega_i` – Angular velocity of particle i \[rad/s\].
    /// * `omega_j` – Angular velocity of particle j \[rad/s\].
    pub fn compute_rolling_resistance(&self, omega_i: [f64; 3], omega_j: [f64; 3]) -> [f64; 3] {
        if self.overlap <= 0.0 {
            return [0.0; 3];
        }
        let fn_mag = (4.0 / 3.0) * self.e_eff * self.r_eff.sqrt() * self.overlap.powf(1.5);
        let dw = [
            omega_i[0] - omega_j[0],
            omega_i[1] - omega_j[1],
            omega_i[2] - omega_j[2],
        ];
        let dw_mag = (dw[0] * dw[0] + dw[1] * dw[1] + dw[2] * dw[2]).sqrt();
        if dw_mag < 1e-30 {
            return [0.0; 3];
        }
        let scale = -self.mu_rolling * self.r_eff * fn_mag / dw_mag;
        [scale * dw[0], scale * dw[1], scale * dw[2]]
    }
}
/// A single granular particle with full stress-state tracking.
pub struct GranularParticle {
    /// Position in 3-D space.
    pub position: [f64; 3],
    /// Velocity vector.
    pub velocity: [f64; 3],
    /// Mass of the particle.
    pub mass: f64,
    /// Current density.
    pub density: f64,
    /// Hydrostatic pressure (positive = compression).
    pub pressure: f64,
    /// Full Cauchy stress tensor (row-major 3×3).
    pub stress_tensor: [[f64; 3]; 3],
    /// Local void fraction (0 = fully packed, 1 = fully void).
    pub void_fraction: f64,
}
/// Hertz-Mindlin contact mechanics for a sphere-sphere pair.
///
/// Stores the effective elastic and shear moduli for the contact and provides
/// methods to compute normal (Hertz) and tangential (Mindlin) contact forces.
#[derive(Debug, Clone)]
#[allow(dead_code)]
#[allow(non_snake_case)]
pub struct HertzContact {
    /// Combined elastic modulus E* \[Pa\]:
    /// `1/E* = (1 - ν₁²)/E₁ + (1 - ν₂²)/E₂`.
    pub E_star: f64,
    /// Combined shear modulus G* \[Pa\]:
    /// `1/G* = (2 - ν₁)/G₁ + (2 - ν₂)/G₂`.
    pub G_star: f64,
}
#[allow(dead_code)]
#[allow(non_snake_case)]
impl HertzContact {
    /// Create a `HertzContact` from two materials.
    ///
    /// # Arguments
    /// * `E1, nu1` – Young's modulus and Poisson ratio of particle 1.
    /// * `E2, nu2` – Young's modulus and Poisson ratio of particle 2.
    pub fn from_materials(E1: f64, nu1: f64, E2: f64, nu2: f64) -> Self {
        let E_star = 1.0 / ((1.0 - nu1 * nu1) / E1 + (1.0 - nu2 * nu2) / E2);
        let G1 = E1 / (2.0 * (1.0 + nu1));
        let G2 = E2 / (2.0 * (1.0 + nu2));
        let G_star = 1.0 / ((2.0 - nu1) / G1 + (2.0 - nu2) / G2);
        Self { E_star, G_star }
    }
    /// Create a `HertzContact` directly from effective moduli.
    pub fn new(E_star: f64, G_star: f64) -> Self {
        Self { E_star, G_star }
    }
    /// Hertz normal contact force magnitude.
    ///
    /// `F_n = (4/3) * E* * sqrt(R*) * δ_n^{3/2}`
    ///
    /// Returns 0 for non-positive overlap.
    pub fn normal_force(&self, r_star: f64, delta_n: f64) -> f64 {
        if delta_n <= 0.0 {
            return 0.0;
        }
        (4.0 / 3.0) * self.E_star * r_star.sqrt() * delta_n.powf(1.5)
    }
    /// Mindlin tangential (shear) contact force vector (no-slip regime).
    ///
    /// `F_t = 8 * G* * sqrt(R* * δ_n) * δ_t`
    ///
    /// The force is clamped to the Coulomb cone `|F_t| ≤ μ * |F_n|`.
    ///
    /// # Arguments
    /// * `r_star`   – Effective contact radius R* \[m\].
    /// * `delta_n`  – Normal overlap δ_n \[m\] (must be positive).
    /// * `delta_t`  – Tangential displacement vector \[m\].
    /// * `mu`       – Coulomb friction coefficient.
    /// * `fn_mag`   – Normal force magnitude \[N\] (for cone clamp).
    pub fn tangential_force(
        &self,
        r_star: f64,
        delta_n: f64,
        delta_t: [f64; 3],
        mu: f64,
        fn_mag: f64,
    ) -> [f64; 3] {
        if delta_n <= 0.0 {
            return [0.0; 3];
        }
        let kt = 8.0 * self.G_star * (r_star * delta_n).sqrt();
        let ft_raw = [-kt * delta_t[0], -kt * delta_t[1], -kt * delta_t[2]];
        let ft_mag = (ft_raw[0] * ft_raw[0] + ft_raw[1] * ft_raw[1] + ft_raw[2] * ft_raw[2]).sqrt();
        let limit = mu * fn_mag.abs();
        if ft_mag > limit && ft_mag > 1e-30 {
            let scale = limit / ft_mag;
            [ft_raw[0] * scale, ft_raw[1] * scale, ft_raw[2] * scale]
        } else {
            ft_raw
        }
    }
    /// Compute the Hertz contact stiffness `k_n = 2 * E* * sqrt(R* * δ_n)`.
    pub fn stiffness(&self, r_star: f64, delta_n: f64) -> f64 {
        if delta_n <= 0.0 {
            return 0.0;
        }
        2.0 * self.E_star * (r_star * delta_n).sqrt()
    }
}
/// Selects the normal contact law used in DEM calculations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GranularContactModel {
    /// Hertz nonlinear elastic contact (F ∝ δ^{3/2}).
    Hertz,
    /// Linear spring model (F = k_n * δ).
    LinearSpring,
    /// Mindlin-Deresiewicz tangential contact (requires Hertz normal force).
    Mindlin,
}
/// Rolling resistance model type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum RollingResistanceModel {
    /// No rolling resistance.
    None,
    /// Constant torque model: T_r = -mu_r * r * F_n * omega_hat.
    Constant,
    /// Viscous rolling resistance: T_r = -mu_r * r * F_n * omega_rel.
    Viscous,
    /// Elastic-plastic rolling resistance.
    ElasticPlastic,
}
/// Hertz-Mindlin contact parameters for grain-grain collisions.
pub struct GranularCollision {
    /// Normal stiffness coefficient \[N/m^1.5\].
    pub k_n: f64,
    /// Normal damping coefficient \[N·s/m^1.25\].
    pub gamma_n: f64,
    /// Tangential stiffness coefficient \[N/m\].
    pub k_t: f64,
    /// Coulomb friction coefficient.
    pub mu: f64,
}
impl GranularCollision {
    /// Construct contact parameters from elastic moduli.
    ///
    /// Uses simplified Hertz contact theory:
    /// - k_n = (4/3) · E* (here E* ≈ E / (2(1-ν²)) for identical spheres)
    /// - γ_n derived from restitution ≈ 0.3 (light damping)
    /// - k_t = 2/7 · k_n (typical approximation)
    ///
    /// # Arguments
    /// * `e_modulus` – Young's modulus E \[Pa\]
    /// * `poisson`   – Poisson's ratio ν
    /// * `mu`        – Coulomb friction coefficient
    pub fn new(e_modulus: f64, poisson: f64, mu: f64) -> Self {
        let e_star = e_modulus / (2.0 * (1.0 - poisson * poisson));
        let k_n = (4.0 / 3.0) * e_star;
        let gamma_n = 0.1 * (2.0 * k_n).sqrt();
        let k_t = (2.0 / 7.0) * k_n;
        Self {
            k_n,
            gamma_n,
            k_t,
            mu,
        }
    }
    /// Hertz normal contact force.
    ///
    /// F_n = k_n · δ^1.5 − γ_n · v_rel_n · δ^0.25
    ///
    /// The result is clamped to ≥ 0 (no tensile adhesion in Hertz contact).
    ///
    /// # Arguments
    /// * `overlap`  – normal overlap δ \[m\] (≥ 0)
    /// * `v_rel_n`  – relative normal velocity (positive = approaching) \[m/s\]
    pub fn normal_force(&self, overlap: f64, v_rel_n: f64) -> f64 {
        if overlap <= 0.0 {
            return 0.0;
        }
        let f = self.k_n * overlap.powf(1.5) - self.gamma_n * v_rel_n * overlap.powf(0.25);
        f.max(0.0)
    }
    /// Tangential (sliding) force, clamped by Coulomb friction.
    ///
    /// F_t = k_t · δ_t, clamped to μ · F_n
    ///
    /// # Arguments
    /// * `delta_t` – tangential displacement increment δ_t \[m\]
    /// * `f_n`     – current normal force magnitude \[N\]
    pub fn tangential_force(&self, delta_t: f64, f_n: f64) -> f64 {
        let f_t = self.k_t * delta_t;
        let limit = self.mu * f_n.abs();
        f_t.clamp(-limit, limit)
    }
}
/// Material parameters for a granular medium.
pub struct GranularParams {
    /// Cohesion intercept c \[Pa\].
    pub cohesion: f64,
    /// Internal friction angle φ \[degrees\].
    pub friction_angle_deg: f64,
    /// Bulk modulus K \[Pa\].
    pub bulk_modulus: f64,
    /// Shear modulus G \[Pa\].
    pub shear_modulus: f64,
    /// Coefficient of restitution (0 = perfectly plastic, 1 = elastic).
    pub restitution: f64,
}
impl GranularParams {
    /// Typical dry sand: friction angle 30°, negligible cohesion.
    pub fn sand() -> Self {
        Self {
            cohesion: 1.0e2,
            friction_angle_deg: 30.0,
            bulk_modulus: 5.0e7,
            shear_modulus: 2.5e7,
            restitution: 0.5,
        }
    }
    /// Wet sand: friction angle 25°, moderate cohesion from capillary forces.
    pub fn wet_sand() -> Self {
        Self {
            cohesion: 5.0e3,
            friction_angle_deg: 25.0,
            bulk_modulus: 6.0e7,
            shear_modulus: 3.0e7,
            restitution: 0.4,
        }
    }
    /// Fine cohesive powder: high cohesion, moderate friction.
    pub fn cohesive_powder() -> Self {
        Self {
            cohesion: 5.0e4,
            friction_angle_deg: 35.0,
            bulk_modulus: 1.0e7,
            shear_modulus: 5.0e6,
            restitution: 0.3,
        }
    }
}
/// A full Discrete Element Method (DEM) simulation with:
/// - Hertz normal contact forces
/// - Mindlin tangential forces with Coulomb cone
/// - Angular velocity / torque integration
/// - Gravity
///
/// Particles are [`DemParticle`] instances and the contact model is driven by
/// a [`HertzContact`] instance that is shared across all particle pairs.
#[derive(Debug)]
#[allow(dead_code)]
pub struct GranularSimulation {
    /// All rigid-sphere particles.
    pub particles: Vec<DemParticle>,
    /// Gravity vector \[m/s²\].
    pub gravity: [f64; 3],
    /// Hertz-Mindlin contact model.
    pub contact: HertzContact,
    /// Coulomb friction coefficient.
    pub mu: f64,
    /// Normal damping coefficient (viscous dashpot) \[N·s/m\].
    pub gamma_n: f64,
    /// Integration time step \[s\].
    pub dt: f64,
}
#[allow(dead_code)]
#[allow(non_snake_case)]
impl GranularSimulation {
    /// Create a new `GranularSimulation`.
    ///
    /// # Arguments
    /// * `contact`  – Hertz-Mindlin contact model.
    /// * `gravity`  – Gravitational acceleration \[m/s²\].
    /// * `mu`       – Coulomb friction coefficient.
    /// * `gamma_n`  – Normal viscous damping \[N·s/m\].
    /// * `dt`       – Integration time step \[s\].
    pub fn new(contact: HertzContact, gravity: [f64; 3], mu: f64, gamma_n: f64, dt: f64) -> Self {
        Self {
            particles: Vec::new(),
            gravity,
            contact,
            mu,
            gamma_n,
            dt,
        }
    }
    /// Add a sphere particle.
    pub fn add_particle(&mut self, position: [f64; 3], velocity: [f64; 3], radius: f64, mass: f64) {
        self.particles
            .push(DemParticle::new_sphere(position, velocity, radius, mass));
    }
    /// Number of particles.
    pub fn len(&self) -> usize {
        self.particles.len()
    }
    /// True if no particles.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }
    /// Apply gravitational acceleration to a force accumulator.
    pub fn apply_gravity(&self, forces: &mut Vec<[f64; 3]>) {
        for (i, f) in forces.iter_mut().enumerate() {
            let m = self.particles[i].mass;
            f[0] += m * self.gravity[0];
            f[1] += m * self.gravity[1];
            f[2] += m * self.gravity[2];
        }
    }
    /// Detect all pairwise sphere-sphere contacts (O(N²)).
    ///
    /// Returns a list of `(i, j, normal, overlap, r_star)` tuples for
    /// all overlapping pairs.
    pub fn detect_contacts(&self) -> Vec<(usize, usize, [f64; 3], f64, f64)> {
        let n = self.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.particles[j].position[0] - self.particles[i].position[0];
                let dy = self.particles[j].position[1] - self.particles[i].position[1];
                let dz = self.particles[j].position[2] - self.particles[i].position[2];
                let dist2 = dx * dx + dy * dy + dz * dz;
                let dist = dist2.sqrt();
                let sum_r = self.particles[i].radius + self.particles[j].radius;
                if dist < sum_r && dist > 1e-20 {
                    let overlap = sum_r - dist;
                    let inv_d = 1.0 / dist;
                    let normal = [dx * inv_d, dy * inv_d, dz * inv_d];
                    let r_star = self.particles[i].effective_radius(&self.particles[j]);
                    contacts.push((i, j, normal, overlap, r_star));
                }
            }
        }
        contacts
    }
    /// Advance by one time step `dt` using velocity Verlet (symplectic Euler).
    ///
    /// Pipeline:
    /// 1. Clear force/torque accumulators.
    /// 2. Apply gravity.
    /// 3. Detect contacts → accumulate Hertz + Mindlin forces and torques.
    /// 4. Integrate: `v += (F/m)*dt`, `x += v*dt`, `ω += (τ/I)*dt`.
    pub fn step(&mut self, dt: f64) {
        let n = self.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut torques = vec![[0.0_f64; 3]; n];
        self.apply_gravity(&mut forces);
        let contacts = self.detect_contacts();
        for (i, j, normal, delta_n, r_star) in contacts {
            let rel_v = [
                self.particles[j].velocity[0] - self.particles[i].velocity[0],
                self.particles[j].velocity[1] - self.particles[i].velocity[1],
                self.particles[j].velocity[2] - self.particles[i].velocity[2],
            ];
            let vn = rel_v[0] * normal[0] + rel_v[1] * normal[1] + rel_v[2] * normal[2];
            let fn_hertz = self.contact.normal_force(r_star, delta_n);
            let fn_damp = -self.gamma_n * vn;
            let fn_total = (fn_hertz + fn_damp).max(0.0);
            let vt = [
                rel_v[0] - vn * normal[0],
                rel_v[1] - vn * normal[1],
                rel_v[2] - vn * normal[2],
            ];
            let delta_t = [vt[0] * dt, vt[1] * dt, vt[2] * dt];
            let ft = self
                .contact
                .tangential_force(r_star, delta_n, delta_t, self.mu, fn_total);
            for k in 0..3 {
                forces[i][k] -= fn_total * normal[k];
                forces[j][k] += fn_total * normal[k];
                forces[i][k] += ft[k];
                forces[j][k] -= ft[k];
            }
            let ri = self.particles[i].radius;
            let rj = self.particles[j].radius;
            torques[i][0] += ri * (normal[1] * ft[2] - normal[2] * ft[1]);
            torques[i][1] += ri * (normal[2] * ft[0] - normal[0] * ft[2]);
            torques[i][2] += ri * (normal[0] * ft[1] - normal[1] * ft[0]);
            torques[j][0] -= rj * (normal[1] * ft[2] - normal[2] * ft[1]);
            torques[j][1] -= rj * (normal[2] * ft[0] - normal[0] * ft[2]);
            torques[j][2] -= rj * (normal[0] * ft[1] - normal[1] * ft[0]);
        }
        for i in 0..n {
            let inv_m = 1.0 / self.particles[i].mass;
            let inv_I = 1.0 / self.particles[i].moment_of_inertia.max(1e-30);
            for k in 0..3 {
                self.particles[i].velocity[k] += forces[i][k] * inv_m * dt;
                self.particles[i].position[k] += self.particles[i].velocity[k] * dt;
                self.particles[i].angular_velocity[k] += torques[i][k] * inv_I * dt;
            }
        }
    }
    /// Total translational kinetic energy.
    pub fn kinetic_energy_translational(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.velocity[0] * p.velocity[0]
                    + p.velocity[1] * p.velocity[1]
                    + p.velocity[2] * p.velocity[2];
                0.5 * p.mass * v2
            })
            .sum()
    }
    /// Total rotational kinetic energy.
    pub fn kinetic_energy_rotational(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let w2 = p.angular_velocity[0] * p.angular_velocity[0]
                    + p.angular_velocity[1] * p.angular_velocity[1]
                    + p.angular_velocity[2] * p.angular_velocity[2];
                0.5 * p.moment_of_inertia * w2
            })
            .sum()
    }
    /// Total kinetic energy (translational + rotational).
    pub fn kinetic_energy(&self) -> f64 {
        self.kinetic_energy_translational() + self.kinetic_energy_rotational()
    }
    /// Total linear momentum vector `[px, py, pz]`.
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut mom = [0.0_f64; 3];
        for p in &self.particles {
            mom[0] += p.mass * p.velocity[0];
            mom[1] += p.mass * p.velocity[1];
            mom[2] += p.mass * p.velocity[2];
        }
        mom
    }
}
/// A simple spatial hash cell for contact detection.
#[allow(dead_code)]
pub struct SpatialHashDem {
    pub(super) cell_size: f64,
    pub(super) buckets: std::collections::HashMap<(i64, i64, i64), Vec<usize>>,
}
#[allow(dead_code)]
impl SpatialHashDem {
    /// Build a spatial hash over particle positions with the given cell size.
    pub fn new(positions: &[[f64; 3]], radii: &[f64], max_radius: f64) -> Self {
        let cell_size = 2.0 * max_radius;
        let mut buckets: std::collections::HashMap<(i64, i64, i64), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, pos) in positions.iter().enumerate() {
            let _ = radii;
            let key = Self::cell_key(pos, cell_size);
            buckets.entry(key).or_default().push(i);
        }
        Self { cell_size, buckets }
    }
    fn cell_key(pos: &[f64; 3], cell_size: f64) -> (i64, i64, i64) {
        (
            (pos[0] / cell_size).floor() as i64,
            (pos[1] / cell_size).floor() as i64,
            (pos[2] / cell_size).floor() as i64,
        )
    }
    /// Find all candidate contact pairs using the spatial hash.
    ///
    /// Returns pairs `(i, j)` where i < j and both particles could be overlapping.
    pub fn candidate_pairs(&self, positions: &[[f64; 3]], radii: &[f64]) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        let cs = self.cell_size;
        for (&(cx, cy, cz), particles_in_cell) in &self.buckets {
            for dx in -1_i64..=1 {
                for dy in -1_i64..=1 {
                    for dz in -1_i64..=1 {
                        let nkey = (cx + dx, cy + dy, cz + dz);
                        if let Some(neighbors) = self.buckets.get(&nkey) {
                            for &i in particles_in_cell {
                                for &j in neighbors {
                                    if j <= i {
                                        continue;
                                    }
                                    let px = positions[i][0] - positions[j][0];
                                    let py = positions[i][1] - positions[j][1];
                                    let pz = positions[i][2] - positions[j][2];
                                    let dist2 = px * px + py * py + pz * pz;
                                    let contact_dist = radii[i] + radii[j] + cs * 0.01;
                                    if dist2 < contact_dist * contact_dist {
                                        pairs.push((i, j));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        pairs
    }
}
/// A rigid-sphere particle used in the Discrete Element Method (DEM).
#[derive(Debug, Clone)]
pub struct DemParticle {
    /// Centre position \[m\].
    pub position: [f64; 3],
    /// Translational velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Particle radius \[m\].
    pub radius: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Moment of inertia \[kg·m²\] (sphere: 2/5 * m * r²).
    pub moment_of_inertia: f64,
}
impl DemParticle {
    /// Create a sphere with given radius, mass, and position.
    ///
    /// Moment of inertia is automatically set to `2/5 * mass * radius²`.
    pub fn new_sphere(position: [f64; 3], velocity: [f64; 3], radius: f64, mass: f64) -> Self {
        let moi = 0.4 * mass * radius * radius;
        Self {
            position,
            velocity,
            angular_velocity: [0.0; 3],
            radius,
            mass,
            moment_of_inertia: moi,
        }
    }
    /// Compute the overlap (positive when spheres interpenetrate) with another particle.
    pub fn overlap_with(&self, other: &DemParticle) -> f64 {
        let dx = other.position[0] - self.position[0];
        let dy = other.position[1] - self.position[1];
        let dz = other.position[2] - self.position[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let gap = dist - (self.radius + other.radius);
        -gap
    }
    /// Effective radius for a contact with `other`: R* = (R_i * R_j) / (R_i + R_j).
    pub fn effective_radius(&self, other: &DemParticle) -> f64 {
        (self.radius * other.radius) / (self.radius + other.radius)
    }
    /// Effective mass: m* = (m_i * m_j) / (m_i + m_j).
    pub fn effective_mass(&self, other: &DemParticle) -> f64 {
        (self.mass * other.mass) / (self.mass + other.mass)
    }
}
/// Simple DEM simulation: sphere packing with pairwise contact resolution.
pub struct DemSimulation {
    /// All particles in the simulation.
    pub particles: Vec<DemParticle>,
    /// Integration time step \[s\].
    pub dt: f64,
    /// Gravitational acceleration \[m/s²\].
    pub gravity: [f64; 3],
    /// Normal spring stiffness \[N/m\].
    pub kn: f64,
    /// Tangential spring stiffness \[N/m\].
    pub kt: f64,
    /// Normal damping coefficient \[N·s/m\].
    pub gamma_n: f64,
    /// Coulomb friction coefficient.
    pub mu: f64,
    /// Contact model to use.
    pub contact_model: GranularContactModel,
    /// Effective modulus for Hertz model \[Pa\].
    pub e_eff: f64,
}
impl DemSimulation {
    /// Create a new `DemSimulation`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        dt: f64,
        gravity: [f64; 3],
        kn: f64,
        kt: f64,
        gamma_n: f64,
        mu: f64,
        contact_model: GranularContactModel,
        e_eff: f64,
    ) -> Self {
        Self {
            particles: Vec::new(),
            dt,
            gravity,
            kn,
            kt,
            gamma_n,
            mu,
            contact_model,
            e_eff,
        }
    }
    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: DemParticle) {
        self.particles.push(p);
    }
    /// Advance the simulation by one time step.
    ///
    /// 1. Clears forces/torques.
    /// 2. Applies gravity.
    /// 3. Detects contacts and accumulates contact forces.
    /// 4. Integrates positions and velocities (forward Euler).
    pub fn step(&mut self) {
        let n = self.particles.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut torques = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            for k in 0..3 {
                forces[i][k] += self.particles[i].mass * self.gravity[k];
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.particles[j].position[0] - self.particles[i].position[0];
                let dy = self.particles[j].position[1] - self.particles[i].position[1];
                let dz = self.particles[j].position[2] - self.particles[i].position[2];
                let dist2 = dx * dx + dy * dy + dz * dz;
                let dist = dist2.sqrt();
                let sum_r = self.particles[i].radius + self.particles[j].radius;
                if dist >= sum_r || dist < 1e-20 {
                    continue;
                }
                let nx = dx / dist;
                let ny = dy / dist;
                let nz = dz / dist;
                let delta_n = sum_r - dist;
                let fn_mag = match self.contact_model {
                    GranularContactModel::Hertz | GranularContactModel::Mindlin => {
                        let r_eff = self.particles[i].effective_radius(&self.particles[j]);
                        dem_hertz_normal_force(self.e_eff, r_eff, delta_n)
                    }
                    GranularContactModel::LinearSpring => linear_spring_normal(self.kn, delta_n),
                };
                let dvx = self.particles[j].velocity[0] - self.particles[i].velocity[0];
                let dvy = self.particles[j].velocity[1] - self.particles[i].velocity[1];
                let dvz = self.particles[j].velocity[2] - self.particles[i].velocity[2];
                let vn = dvx * nx + dvy * ny + dvz * nz;
                let fn_damp = self.gamma_n * vn;
                let fn_total = (fn_mag - fn_damp).max(0.0);
                forces[i][0] -= fn_total * nx;
                forces[i][1] -= fn_total * ny;
                forces[i][2] -= fn_total * nz;
                forces[j][0] += fn_total * nx;
                forces[j][1] += fn_total * ny;
                forces[j][2] += fn_total * nz;
                let vt_x = dvx - vn * nx;
                let vt_y = dvy - vn * ny;
                let vt_z = dvz - vn * nz;
                let delta_t = [vt_x * self.dt, vt_y * self.dt, vt_z * self.dt];
                let ft = tangential_friction_force(self.kt, delta_t, self.mu, fn_total);
                forces[i][0] += ft[0];
                forces[i][1] += ft[1];
                forces[i][2] += ft[2];
                forces[j][0] -= ft[0];
                forces[j][1] -= ft[1];
                forces[j][2] -= ft[2];
                let ri = self.particles[i].radius;
                let rj = self.particles[j].radius;
                torques[i][0] += ri * (ny * ft[2] - nz * ft[1]);
                torques[i][1] += ri * (nz * ft[0] - nx * ft[2]);
                torques[i][2] += ri * (nx * ft[1] - ny * ft[0]);
                torques[j][0] -= rj * (ny * ft[2] - nz * ft[1]);
                torques[j][1] -= rj * (nz * ft[0] - nx * ft[2]);
                torques[j][2] -= rj * (nx * ft[1] - ny * ft[0]);
            }
        }
        let dt = self.dt;
        for i in 0..n {
            let inv_m = 1.0 / self.particles[i].mass;
            let inv_i = 1.0 / self.particles[i].moment_of_inertia.max(1e-30);
            for k in 0..3 {
                self.particles[i].velocity[k] += forces[i][k] * inv_m * dt;
                self.particles[i].position[k] += self.particles[i].velocity[k] * dt;
                self.particles[i].angular_velocity[k] += torques[i][k] * inv_i * dt;
            }
        }
    }
    /// Total translational kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.velocity[0] * p.velocity[0]
                    + p.velocity[1] * p.velocity[1]
                    + p.velocity[2] * p.velocity[2];
                0.5 * p.mass * v2
            })
            .sum()
    }
    /// Number of particles.
    pub fn num_particles(&self) -> usize {
        self.particles.len()
    }
}
impl DemSimulation {
    /// Compute the average number of contacts per particle (coordination number).
    ///
    /// The coordination number Z for particle `i` is the number of distinct
    /// contacts involving `i`.  The *average* coordination number is
    ///   `Z̄ = 2 * N_contacts / N_particles`
    ///
    /// Returns 0.0 for empty simulations.
    ///
    /// # Arguments
    /// * `contact_threshold` – Minimum overlap (in meters) to count a contact.
    ///   Use 0.0 to count all overlapping pairs.
    pub fn compute_coordination_number(&self, contact_threshold: f64) -> f64 {
        let n = self.particles.len();
        if n == 0 {
            return 0.0;
        }
        let mut n_contacts: usize = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let dx = self.particles[j].position[0] - self.particles[i].position[0];
                let dy = self.particles[j].position[1] - self.particles[i].position[1];
                let dz = self.particles[j].position[2] - self.particles[i].position[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                let sum_r = self.particles[i].radius + self.particles[j].radius;
                let overlap = sum_r - dist;
                if overlap > contact_threshold {
                    n_contacts += 1;
                }
            }
        }
        2.0 * n_contacts as f64 / n as f64
    }
}
/// Contact model for grain-grain interaction.
#[derive(Debug, Clone)]
pub enum ContactModel {
    /// Linear spring-dashpot (Cundall & Strack style).
    SpringDashpot {
        /// Normal spring stiffness \[N/m\].
        kn: f64,
        /// Tangential spring stiffness \[N/m\].
        kt: f64,
        /// Normal damping coefficient \[N·s/m\].
        gamma_n: f64,
        /// Tangential damping coefficient \[N·s/m\].
        gamma_t: f64,
    },
    /// Hertzian contact (nonlinear).
    Hertzian {
        /// Effective Young's modulus E* \[Pa\].
        e_star: f64,
        /// Effective radius R* \[m\].
        r_star: f64,
    },
}
