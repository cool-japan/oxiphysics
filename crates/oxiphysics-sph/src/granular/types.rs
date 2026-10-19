//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::BTreeMap;
use std::f64::consts::PI;

use super::functions::*;

/// Extended granular simulation that also tracks angular velocities.
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
    /// Advance one step with tangential friction and rolling resistance.
    ///
    /// The contacts are resolved with the particles' current spin, so tangential
    /// slip drives Coulomb-limited friction forces whose torques spin the particles
    /// up.  The rolling-resistance torque then acts on the relative spin of each
    /// contact using that contact's own normal force.
    pub fn step(&mut self, dt: f64) {
        let contacts = self.sim.detect_contacts();
        let loads = self
            .sim
            .resolve_contacts_with_spin(&contacts, dt, &self.angular_velocities);
        let n = self.len();
        for t in self.torques.iter_mut().take(n) {
            *t = [0.0; 3];
        }
        for (c, load) in contacts.iter().zip(loads.iter()) {
            let r_i = self.sim.radii[c.i];
            let r_j = self.sim.radii[c.j];
            // Friction acts at the midpoint of the overlap region, a distance
            // r − δ/2 from each centre.  `c.normal` points from j to i, so the
            // contact point lies at −lever·n from i and at +lever·n from j and the
            // friction force F_t on i (−F_t on j) produces
            // τᵢ = −leverᵢ n × F_t and τⱼ = −leverⱼ n × F_t.
            let lever_i = r_i - 0.5 * c.overlap;
            let lever_j = r_j - 0.5 * c.overlap;
            let n_cross_ft = cross3(c.normal, load.tangential);
            add_scaled(&mut self.torques[c.i], -lever_i, n_cross_ft);
            add_scaled(&mut self.torques[c.j], -lever_j, n_cross_ft);
            let omega_rel = [
                self.angular_velocities[c.i][0] - self.angular_velocities[c.j][0],
                self.angular_velocities[c.i][1] - self.angular_velocities[c.j][1],
                self.angular_velocities[c.i][2] - self.angular_velocities[c.j][2],
            ];
            let torque = rolling_resistance_torque(
                r_i,
                r_j,
                load.normal,
                omega_rel,
                self.mu_r,
                self.rolling_model,
            );
            add_scaled(&mut self.torques[c.i], 1.0, torque);
            add_scaled(&mut self.torques[c.j], -1.0, torque);
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
/// Force produced by one contact of [`GranularSim::resolve_contacts`], listed in
/// the order of the [`Contact`] slice that was resolved.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactForce {
    /// Index of the first particle (`Contact::i`).
    pub i: usize,
    /// Index of the second particle (`Contact::j`).
    pub j: usize,
    /// Normal force magnitude `F_n ≥ 0` \[N\], repulsive along `Contact::normal`
    /// (from j to i).
    pub normal: f64,
    /// Tangential (friction) force on particle `i` \[N\]; particle `j` receives the
    /// reaction `−tangential`.
    pub tangential: [f64; 3],
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
///
/// Every contact carries a normal force (spring-dashpot or Hertz) and a
/// Coulomb-limited tangential friction force with a per-contact spring history
/// ([`tangential_springs`](Self::tangential_springs)).  The spheres themselves are
/// translational only; [`GranularSimExtended`] adds their spin.
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
    /// Tangential spring displacement `δ_t` \[m\] of every persisting contact, keyed
    /// by the particle pair `(i, j)` with `i < j`.
    ///
    /// Each [`resolve_contacts`](Self::resolve_contacts) call advances the entry of
    /// every contact it is given and drops the pairs that are no longer in contact.
    /// The keys are particle indices: call `tangential_springs.clear()` after
    /// reordering or removing particles.
    pub tangential_springs: BTreeMap<(usize, usize), [f64; 3]>,
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
            tangential_springs: BTreeMap::new(),
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
    /// Resolve detected contacts by computing and applying contact forces, for
    /// spheres without spin.
    ///
    /// Equivalent to [`resolve_contacts_with_spin`](Self::resolve_contacts_with_spin)
    /// with zero angular velocities.  `dt` advances the tangential springs.
    pub fn resolve_contacts(&mut self, contacts: &[Contact], dt: f64) -> Vec<ContactForce> {
        self.resolve_contacts_with_spin(contacts, dt, &[])
    }
    /// Resolve detected contacts by computing and applying contact forces.
    ///
    /// The forces are reset to gravity first.  For every contact:
    /// - the normal force `F_n` comes from the contact model;
    /// - the sliding velocity is the tangential part of the contact-point velocity
    ///   of `i` relative to `j`, including the spin
    ///   (see [`contact_point_relative_velocity`]; missing entries of
    ///   `angular_velocities` count as zero);
    /// - the friction force `F_t = −k_t δ_t − γ_t v_t`, limited to `μ F_n`, is
    ///   applied to `i` and its reaction to `j`.  `k_t` and `γ_t` are the
    ///   `SpringDashpot` parameters; the `Hertzian` model uses the Mindlin
    ///   stiffness `k_t = 8 G* √(R* δ_n)` without damping.
    ///
    /// Returns one [`ContactForce`] per contact, in the order of `contacts`.
    pub fn resolve_contacts_with_spin(
        &mut self,
        contacts: &[Contact],
        dt: f64,
        angular_velocities: &[[f64; 3]],
    ) -> Vec<ContactForce> {
        for i in 0..self.len() {
            self.forces[i] = [
                self.masses[i] * self.gravity[0],
                self.masses[i] * self.gravity[1],
                self.masses[i] * self.gravity[2],
            ];
        }
        let mut previous = std::mem::take(&mut self.tangential_springs);
        let mut current = BTreeMap::new();
        let mut loads = Vec::with_capacity(contacts.len());
        for c in contacts {
            let omega_i = angular_velocities.get(c.i).copied().unwrap_or([0.0; 3]);
            let omega_j = angular_velocities.get(c.j).copied().unwrap_or([0.0; 3]);
            let lever_i = self.radii[c.i] - 0.5 * c.overlap;
            let lever_j = self.radii[c.j] - 0.5 * c.overlap;
            // `c.normal` points from j to i; the contact-point velocity is
            // expressed with the normal from i to j.
            let normal_ij = [-c.normal[0], -c.normal[1], -c.normal[2]];
            let slip = contact_point_relative_velocity(
                self.velocities[c.i],
                self.velocities[c.j],
                omega_i,
                omega_j,
                lever_i,
                lever_j,
                normal_ij,
            );
            // Spin is perpendicular to the normal, so this is the translational
            // normal velocity (positive = separating).
            let vn = dot3(slip, c.normal);
            let fn_mag = match &self.contact_model {
                ContactModel::SpringDashpot { kn, gamma_n, .. } => {
                    spring_dashpot_normal(c.overlap, vn, *kn, *gamma_n)
                }
                ContactModel::Hertzian { e_star, r_star, .. } => {
                    hertz_normal_force(c.overlap, *e_star, *r_star)
                }
            };
            let law = match &self.contact_model {
                ContactModel::SpringDashpot { kt, gamma_t, .. } => TangentialContactLaw {
                    kt: *kt,
                    gamma_t: *gamma_t,
                    mu: self.mu,
                },
                ContactModel::Hertzian { g_star, r_star, .. } => TangentialContactLaw {
                    kt: 8.0 * g_star * (r_star * c.overlap).sqrt(),
                    gamma_t: 0.0,
                    mu: self.mu,
                },
            };
            let mut spring = previous.remove(&(c.i, c.j)).unwrap_or([0.0; 3]);
            let ft = law.force(
                &mut spring,
                c.normal,
                tangential_component(slip, c.normal),
                fn_mag,
                dt,
            );
            current.insert((c.i, c.j), spring);
            let mut force_on_i = ft;
            add_scaled(&mut force_on_i, fn_mag, c.normal);
            add_scaled(&mut self.forces[c.i], 1.0, force_on_i);
            add_scaled(&mut self.forces[c.j], -1.0, force_on_i);
            loads.push(ContactForce {
                i: c.i,
                j: c.j,
                normal: fn_mag,
                tangential: ft,
            });
        }
        self.tangential_springs = current;
        loads
    }
    /// Advance one step: detect contacts, resolve, integrate.
    pub fn step(&mut self, dt: f64) {
        let contacts = self.detect_contacts();
        self.resolve_contacts(&contacts, dt);
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
pub struct HertzContact {
    /// Combined elastic modulus E* \[Pa\]:
    /// `1/E* = (1 - ν₁²)/E₁ + (1 - ν₂²)/E₂`.
    pub e_star: f64,
    /// Combined shear modulus G* \[Pa\]:
    /// `1/G* = (2 - ν₁)/G₁ + (2 - ν₂)/G₂`.
    pub g_star: f64,
}
impl HertzContact {
    /// Create a `HertzContact` from two materials.
    ///
    /// # Arguments
    /// * `E1, nu1` – Young's modulus and Poisson ratio of particle 1.
    /// * `E2, nu2` – Young's modulus and Poisson ratio of particle 2.
    pub fn from_materials(e1: f64, nu1: f64, e2: f64, nu2: f64) -> Self {
        let e_star = 1.0 / ((1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2);
        let g1 = e1 / (2.0 * (1.0 + nu1));
        let g2 = e2 / (2.0 * (1.0 + nu2));
        let g_star = 1.0 / ((2.0 - nu1) / g1 + (2.0 - nu2) / g2);
        Self { e_star, g_star }
    }
    /// Create a `HertzContact` directly from effective moduli.
    pub fn new(e_star: f64, g_star: f64) -> Self {
        Self { e_star, g_star }
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
        (4.0 / 3.0) * self.e_star * r_star.sqrt() * delta_n.powf(1.5)
    }
    /// Mindlin tangential (shear) contact force vector (no-slip regime).
    ///
    /// `F_t = −k_t δ_t` with `k_t = 8 * G* * sqrt(R* * δ_n)`
    /// (see [`tangential_stiffness`](Self::tangential_stiffness)).
    ///
    /// The force is clamped to the Coulomb cone `|F_t| ≤ μ * |F_n|`.  It acts on the
    /// particle whose contact point has been displaced by `delta_t` relative to its
    /// partner and opposes that displacement.
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
        let kt = self.tangential_stiffness(r_star, delta_n);
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
    /// Mindlin no-slip tangential stiffness `k_t = 8 G* √(R* δ_n)` \[N/m\].
    ///
    /// Returns 0 for non-positive overlap.
    pub fn tangential_stiffness(&self, r_star: f64, delta_n: f64) -> f64 {
        if delta_n <= 0.0 {
            return 0.0;
        }
        8.0 * self.g_star * (r_star * delta_n).sqrt()
    }
    /// Compute the Hertz contact stiffness `k_n = 2 * E* * sqrt(R* * δ_n)`.
    pub fn stiffness(&self, r_star: f64, delta_n: f64) -> f64 {
        if delta_n <= 0.0 {
            return 0.0;
        }
        2.0 * self.e_star * (r_star * delta_n).sqrt()
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
    /// Tangential spring displacement `δ_t` \[m\] of every persisting contact,
    /// keyed by the particle pair `(i, j)` with `i < j`.
    ///
    /// Advanced by [`step`](Self::step), which drops the pairs that separated.
    /// Call `tangential_springs.clear()` after reordering or removing particles.
    pub tangential_springs: BTreeMap<(usize, usize), [f64; 3]>,
}
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
            tangential_springs: BTreeMap::new(),
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
    pub fn apply_gravity(&self, forces: &mut [[f64; 3]]) {
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
    /// Advance by one time step `dt` using symplectic Euler.
    ///
    /// Pipeline:
    /// 1. Clear force/torque accumulators.
    /// 2. Apply gravity.
    /// 3. Detect contacts → accumulate the Hertz normal force and the Mindlin
    ///    tangential force.  The tangential spring integrates the sliding velocity
    ///    of the contact point (translation plus spin, see
    ///    [`contact_point_relative_velocity`]) per contact pair, so friction opposes
    ///    the slip, is capped by `μ F_n`, and its torques spin the particles.
    /// 4. Integrate: `v += (F/m)*dt`, `x += v*dt`, `ω += (τ/I)*dt`.
    pub fn step(&mut self, dt: f64) {
        let n = self.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut torques = vec![[0.0_f64; 3]; n];
        self.apply_gravity(&mut forces);
        let contacts = self.detect_contacts();
        let mut previous = std::mem::take(&mut self.tangential_springs);
        let mut current = BTreeMap::new();
        for (i, j, normal, delta_n, r_star) in contacts {
            let (particle_i, particle_j) = (&self.particles[i], &self.particles[j]);
            // `normal` points from i to j.
            let rel_v = [
                particle_j.velocity[0] - particle_i.velocity[0],
                particle_j.velocity[1] - particle_i.velocity[1],
                particle_j.velocity[2] - particle_i.velocity[2],
            ];
            let vn = dot3(rel_v, normal);
            let fn_hertz = self.contact.normal_force(r_star, delta_n);
            let fn_damp = -self.gamma_n * vn;
            let fn_total = (fn_hertz + fn_damp).max(0.0);
            let lever_i = particle_i.radius - 0.5 * delta_n;
            let lever_j = particle_j.radius - 0.5 * delta_n;
            let slip = contact_point_relative_velocity(
                particle_i.velocity,
                particle_j.velocity,
                particle_i.angular_velocity,
                particle_j.angular_velocity,
                lever_i,
                lever_j,
                normal,
            );
            let law = TangentialContactLaw {
                kt: self.contact.tangential_stiffness(r_star, delta_n),
                gamma_t: 0.0,
                mu: self.mu,
            };
            let mut spring = previous.remove(&(i, j)).unwrap_or([0.0; 3]);
            let ft = law.force(
                &mut spring,
                normal,
                tangential_component(slip, normal),
                fn_total,
                dt,
            );
            current.insert((i, j), spring);
            add_scaled(&mut forces[i], -fn_total, normal);
            add_scaled(&mut forces[j], fn_total, normal);
            add_scaled(&mut forces[i], 1.0, ft);
            add_scaled(&mut forces[j], -1.0, ft);
            // Contact point at +lever_i·n from i and −lever_j·n from j, force ft on
            // i and −ft on j: both torques are lever·(n × ft).
            let n_cross_ft = cross3(normal, ft);
            add_scaled(&mut torques[i], lever_i, n_cross_ft);
            add_scaled(&mut torques[j], lever_j, n_cross_ft);
        }
        self.tangential_springs = current;
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
pub struct SpatialHashDem {
    pub(super) cell_size: f64,
    pub(super) buckets: std::collections::HashMap<(i64, i64, i64), Vec<usize>>,
}
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
    /// Tangential spring displacement `δ_t` \[m\] of every persisting contact,
    /// keyed by the particle pair `(i, j)` with `i < j`.
    ///
    /// Advanced by [`step`](Self::step), which drops the pairs that separated.
    /// Call `tangential_springs.clear()` after reordering or removing particles.
    pub tangential_springs: BTreeMap<(usize, usize), [f64; 3]>,
}
impl DemSimulation {
    /// Create a new `DemSimulation`.
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
            tangential_springs: BTreeMap::new(),
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
    /// 3. Detects contacts and accumulates the normal forces and the
    ///    Coulomb-limited tangential forces.  The tangential spring `kt`
    ///    integrates the sliding velocity of the contact point (translation plus
    ///    spin, see [`contact_point_relative_velocity`]) per contact pair, so the
    ///    friction force opposes the slip and its torques spin the particles.
    /// 4. Integrates positions, velocities and spins (semi-implicit Euler).
    pub fn step(&mut self) {
        let n = self.particles.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut torques = vec![[0.0_f64; 3]; n];
        for (i, f) in forces.iter_mut().enumerate() {
            for (k, fk) in f.iter_mut().enumerate() {
                *fk += self.particles[i].mass * self.gravity[k];
            }
        }
        let mut previous = std::mem::take(&mut self.tangential_springs);
        let mut current = BTreeMap::new();
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
                // Unit normal from i to j.
                let normal = [dx / dist, dy / dist, dz / dist];
                let delta_n = sum_r - dist;
                let fn_mag = match self.contact_model {
                    GranularContactModel::Hertz | GranularContactModel::Mindlin => {
                        let r_eff = self.particles[i].effective_radius(&self.particles[j]);
                        dem_hertz_normal_force(self.e_eff, r_eff, delta_n)
                    }
                    GranularContactModel::LinearSpring => linear_spring_normal(self.kn, delta_n),
                };
                let (particle_i, particle_j) = (&self.particles[i], &self.particles[j]);
                let rel_v = [
                    particle_j.velocity[0] - particle_i.velocity[0],
                    particle_j.velocity[1] - particle_i.velocity[1],
                    particle_j.velocity[2] - particle_i.velocity[2],
                ];
                let vn = dot3(rel_v, normal);
                let fn_damp = self.gamma_n * vn;
                let fn_total = (fn_mag - fn_damp).max(0.0);
                let lever_i = particle_i.radius - 0.5 * delta_n;
                let lever_j = particle_j.radius - 0.5 * delta_n;
                let slip = contact_point_relative_velocity(
                    particle_i.velocity,
                    particle_j.velocity,
                    particle_i.angular_velocity,
                    particle_j.angular_velocity,
                    lever_i,
                    lever_j,
                    normal,
                );
                let law = TangentialContactLaw {
                    kt: self.kt,
                    gamma_t: 0.0,
                    mu: self.mu,
                };
                let mut spring = previous.remove(&(i, j)).unwrap_or([0.0; 3]);
                let ft = law.force(
                    &mut spring,
                    normal,
                    tangential_component(slip, normal),
                    fn_total,
                    self.dt,
                );
                current.insert((i, j), spring);
                add_scaled(&mut forces[i], -fn_total, normal);
                add_scaled(&mut forces[j], fn_total, normal);
                add_scaled(&mut forces[i], 1.0, ft);
                add_scaled(&mut forces[j], -1.0, ft);
                // Contact point at +lever_i·n from i and −lever_j·n from j, force
                // ft on i and −ft on j: both torques are lever·(n × ft).
                let n_cross_ft = cross3(normal, ft);
                add_scaled(&mut torques[i], lever_i, n_cross_ft);
                add_scaled(&mut torques[j], lever_j, n_cross_ft);
            }
        }
        self.tangential_springs = current;
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
    /// Hertz–Mindlin contact (nonlinear): normal force `(4/3) E* √R* δ^{3/2}`,
    /// tangential stiffness `k_t = 8 G* √(R* δ)`.
    Hertzian {
        /// Effective Young's modulus E* \[Pa\].
        e_star: f64,
        /// Effective shear modulus G* \[Pa\]:
        /// `1/G* = (2 − ν₁)/G₁ + (2 − ν₂)/G₂`.
        g_star: f64,
        /// Effective radius R* \[m\].
        r_star: f64,
    },
}
/// Coulomb-limited tangential contact law with a history spring.
///
/// The friction force on particle `i` of a contact is
/// `F_t = −k_t δ_t − γ_t v_t` with the sliding velocity `v_t` of `i` relative to
/// `j` (see [`contact_point_relative_velocity`]) and the spring displacement
/// `δ_t`, which accumulates `v_t Δt` for as long as the contact lasts.  `|F_t|` is
/// capped at `μ F_n`; during slip the spring is rescaled so that the stored force
/// stays on the Coulomb cone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TangentialContactLaw {
    /// Tangential spring stiffness `k_t` \[N/m\].
    pub kt: f64,
    /// Tangential damping coefficient `γ_t` \[N·s/m\].
    pub gamma_t: f64,
    /// Coulomb friction coefficient `μ`.
    pub mu: f64,
}
impl TangentialContactLaw {
    /// Advance the spring displacement of one contact by `v_t Δt` and return the
    /// friction force on particle `i`.
    ///
    /// # Arguments
    /// * `spring` – stored tangential displacement `δ_t` (in/out) \[m\]; it is first
    ///   carried into the tangent plane of the current `normal` (keeping its
    ///   length), is left on the Coulomb cone after a slip and is cleared when the
    ///   contact cannot transmit a tangential force (`μ F_n = 0`).
    /// * `normal` – unit contact normal (its sign is irrelevant).
    /// * `v_t`    – tangential sliding velocity of `i` relative to `j` \[m/s\].
    /// * `fn_mag` – normal force magnitude `F_n ≥ 0` \[N\].
    /// * `dt`     – time step \[s\].
    pub fn force(
        &self,
        spring: &mut [f64; 3],
        normal: [f64; 3],
        v_t: [f64; 3],
        fn_mag: f64,
        dt: f64,
    ) -> [f64; 3] {
        let kt = self.kt.max(0.0);
        let gamma_t = self.gamma_t.max(0.0);
        let limit = self.mu.max(0.0) * fn_mag.max(0.0);
        if limit <= 0.0 {
            // Nothing can be transmitted without a normal load (or without
            // friction): no force, and no history is kept, damper or not.
            *spring = [0.0; 3];
            return [0.0; 3];
        }

        // The contact normal turns during the contact: carry the spring into the
        // current tangent plane without changing its length.
        let old_length = dot3(*spring, *spring).sqrt();
        let projected = tangential_component(*spring, normal);
        let projected_length = dot3(projected, projected).sqrt();
        let mut delta = [0.0_f64; 3];
        if projected_length > 1e-300 {
            add_scaled(&mut delta, old_length / projected_length, projected);
        }
        add_scaled(&mut delta, dt, v_t);

        let mut force = [0.0_f64; 3];
        add_scaled(&mut force, -kt, delta);
        add_scaled(&mut force, -gamma_t, v_t);
        let force_length = dot3(force, force).sqrt();
        if force_length > limit {
            // Slip: scale the force onto the cone and move the spring so that
            // k_t δ_t + γ_t v_t carries the same scale factor.
            let scale = limit / force_length;
            let mut slipped = [0.0_f64; 3];
            if kt > 0.0 {
                let damping_offset = gamma_t / kt;
                add_scaled(&mut slipped, scale, delta);
                add_scaled(&mut slipped, (scale - 1.0) * damping_offset, v_t);
            }
            *spring = slipped;
            return [force[0] * scale, force[1] * scale, force[2] * scale];
        }
        // Without a spring there is no history to keep (a pure dashpot).
        *spring = if kt > 0.0 { delta } else { [0.0; 3] };
        force
    }
}
