#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! DEM-SPH coupling: discrete element method coupled with smoothed-particle
//! hydrodynamics.
//!
//! Provides Hertz contact mechanics, rigid-body DEM integration, and
//! fluid-particle coupling via drag and buoyancy forces.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Helper: 3-D vector arithmetic
// ---------------------------------------------------------------------------

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// DEM particle
// ---------------------------------------------------------------------------

/// A rigid sphere in a DEM simulation.
#[derive(Debug, Clone)]
pub struct DemParticle {
    /// Position of the particle centre \[m\].
    pub position: [f64; 3],
    /// Translational velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Angular (rotational) velocity \[rad/s\].
    pub angular_velocity: [f64; 3],
    /// Particle radius \[m\].
    pub radius: f64,
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Moment of inertia about any axis (sphere: 2/5 * m * r²) \[kg·m²\].
    pub moment_of_inertia: f64,
}

impl DemParticle {
    /// Create a uniform-density sphere.
    ///
    /// `density` is the particle density \[kg/m³\].
    pub fn sphere(position: [f64; 3], radius: f64, density: f64) -> Self {
        let volume = 4.0 / 3.0 * PI * radius.powi(3);
        let mass = density * volume;
        let moi = 2.0 / 5.0 * mass * radius * radius;
        Self {
            position,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            radius,
            mass,
            moment_of_inertia: moi,
        }
    }

    /// Kinetic energy (translational only) \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * vec3_dot(self.velocity, self.velocity)
    }

    /// Rotational kinetic energy \[J\].
    pub fn rotational_energy(&self) -> f64 {
        0.5 * self.moment_of_inertia * vec3_dot(self.angular_velocity, self.angular_velocity)
    }
}

// ---------------------------------------------------------------------------
// Contact forces
// ---------------------------------------------------------------------------

/// Contact force decomposition between two particles.
#[derive(Debug, Clone, Default)]
pub struct ContactForce {
    /// Normal contact force \[N\] (points from particle 1 to particle 2).
    pub normal_force: [f64; 3],
    /// Tangential (friction) force \[N\].
    pub tangential_force: [f64; 3],
    /// Rolling resistance torque \[N·m\].
    pub rolling_torque: [f64; 3],
}

impl ContactForce {
    /// Total force on particle 1 \[N\] (reaction to particle 2 force).
    pub fn total_on_1(&self) -> [f64; 3] {
        vec3_add(
            vec3_scale(self.normal_force, -1.0),
            vec3_scale(self.tangential_force, -1.0),
        )
    }

    /// Total force on particle 2 \[N\].
    pub fn total_on_2(&self) -> [f64; 3] {
        vec3_add(self.normal_force, self.tangential_force)
    }
}

// ---------------------------------------------------------------------------
// Hertz contact
// ---------------------------------------------------------------------------

/// Hertz elastic contact model for spherical particles.
#[derive(Debug, Clone)]
pub struct HertzContact {
    /// Young's modulus of the contact material \[Pa\].
    pub elastic_modulus: f64,
    /// Poisson's ratio of the contact material.
    pub poisson_ratio: f64,
}

impl HertzContact {
    /// Create a new `HertzContact`.
    pub fn new(elastic_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            elastic_modulus,
            poisson_ratio,
        }
    }

    /// Compute the Hertz normal force for a given overlap between two spheres.
    ///
    /// F_n = (4/3) E* √R* δ^(3/2)
    ///
    /// where E* is the effective modulus and R* the effective radius.
    pub fn compute_normal_force(&self, overlap: f64, radius1: f64, radius2: f64) -> f64 {
        if overlap <= 0.0 {
            return 0.0;
        }
        let e_star = self.effective_modulus();
        let r_star = (radius1 * radius2) / (radius1 + radius2);
        (4.0 / 3.0) * e_star * r_star.sqrt() * overlap.powf(1.5)
    }

    /// Effective modulus E* = E / (2(1-ν²)) for identical material.
    pub fn effective_modulus(&self) -> f64 {
        self.elastic_modulus / (2.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }
}

// ---------------------------------------------------------------------------
// DEM simulation
// ---------------------------------------------------------------------------

/// Simple DEM simulation container.
#[derive(Debug, Clone)]
pub struct DemSimulation {
    /// All particles in the simulation.
    pub particles: Vec<DemParticle>,
    /// Contact model.
    pub contact_model: HertzContact,
    /// Coefficient of friction (tangential / normal force ratio limit).
    pub friction_coefficient: f64,
    /// Accumulated forces on each particle \[N\] (length == particles.len()).
    forces: Vec<[f64; 3]>,
    /// Accumulated torques on each particle \[N·m\].
    torques: Vec<[f64; 3]>,
}

impl DemSimulation {
    /// Create a new `DemSimulation`.
    pub fn new(contact_model: HertzContact, friction_coefficient: f64) -> Self {
        Self {
            particles: Vec::new(),
            contact_model,
            friction_coefficient,
            forces: Vec::new(),
            torques: Vec::new(),
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: DemParticle) {
        self.forces.push([0.0; 3]);
        self.torques.push([0.0; 3]);
        self.particles.push(p);
    }

    /// Detect all overlapping particle pairs and return (i, j, overlap, normal).
    pub fn detect_contacts(&self) -> Vec<(usize, usize, f64, [f64; 3])> {
        let n = self.particles.len();
        let mut contacts = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let pi = &self.particles[i];
                let pj = &self.particles[j];
                let delta = vec3_sub(pj.position, pi.position);
                let dist = vec3_norm(delta);
                let sum_r = pi.radius + pj.radius;
                if dist < sum_r && dist > 1.0e-15 {
                    let overlap = sum_r - dist;
                    let normal = vec3_scale(delta, 1.0 / dist);
                    contacts.push((i, j, overlap, normal));
                }
            }
        }
        contacts
    }

    /// Compute and accumulate contact forces on each particle.
    pub fn compute_forces(&mut self) {
        // Reset
        for f in self.forces.iter_mut() {
            *f = [0.0; 3];
        }
        for t in self.torques.iter_mut() {
            *t = [0.0; 3];
        }

        let contacts = self.detect_contacts();
        for (i, j, overlap, normal) in contacts {
            let ri = self.particles[i].radius;
            let rj = self.particles[j].radius;
            let fn_mag = self.contact_model.compute_normal_force(overlap, ri, rj);
            let fn_vec = vec3_scale(normal, fn_mag);

            // Apply Newton's third law
            for k in 0..3 {
                self.forces[i][k] -= fn_vec[k];
                self.forces[j][k] += fn_vec[k];
            }
        }
    }

    /// Integrate one time step with symplectic Euler (velocity Verlet-like).
    pub fn integrate(&mut self, dt: f64) {
        self.compute_forces();
        let n = self.particles.len();
        for idx in 0..n {
            let m = self.particles[idx].mass;
            let i_moi = self.particles[idx].moment_of_inertia;
            // Translational
            for k in 0..3 {
                self.particles[idx].velocity[k] += self.forces[idx][k] / m * dt;
                self.particles[idx].position[k] += self.particles[idx].velocity[k] * dt;
            }
            // Rotational
            for k in 0..3 {
                self.particles[idx].angular_velocity[k] += self.torques[idx][k] / i_moi * dt;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// SPH-DEM coupling
// ---------------------------------------------------------------------------

/// Fluid–particle coupling model using drag and buoyancy.
#[derive(Debug, Clone)]
pub struct SphDemCoupling {
    /// Drag coefficient C_d (dimensionless, sphere ≈ 0.44 at high Re).
    pub drag_coefficient: f64,
    /// Fluid density for buoyancy computation \[kg/m³\].
    pub buoyancy: f64,
}

impl SphDemCoupling {
    /// Create a new `SphDemCoupling`.
    pub fn new(drag_coefficient: f64, buoyancy: f64) -> Self {
        Self {
            drag_coefficient,
            buoyancy,
        }
    }

    /// Compute the total fluid force (drag + buoyancy) on a DEM particle.
    ///
    /// Drag: F_drag = (1/2) C_d ρ A |v_rel| v_rel
    ///
    /// where A is the cross-sectional area and v_rel = v_fluid – v_particle.
    /// Buoyancy: F_b = ρ_fluid × V_particle × g (upward, z-direction).
    pub fn compute_fluid_force(
        &self,
        particle: &DemParticle,
        fluid_vel: [f64; 3],
        fluid_density: f64,
    ) -> [f64; 3] {
        let v_rel = vec3_sub(fluid_vel, particle.velocity);
        let v_rel_mag = vec3_norm(v_rel);
        let cross_area = PI * particle.radius * particle.radius;
        let drag_mag =
            0.5 * self.drag_coefficient * fluid_density * cross_area * v_rel_mag * v_rel_mag;
        let drag = if v_rel_mag > 1.0e-15 {
            vec3_scale(v_rel, drag_mag / v_rel_mag)
        } else {
            [0.0; 3]
        };
        // Buoyancy in z-direction
        const G: f64 = 9.81;
        let vol = 4.0 / 3.0 * PI * particle.radius.powi(3);
        let buoyancy_force = fluid_density * vol * G;
        [drag[0], drag[1], drag[2] + buoyancy_force]
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Effective Hertz contact stiffness for two spheres of different materials.
///
/// 1/E* = (1 − ν₁²)/E₁ + (1 − ν₂²)/E₂
pub fn hertz_stiffness(e1: f64, e2: f64, nu1: f64, nu2: f64) -> f64 {
    let inv_e_star = (1.0 - nu1 * nu1) / e1 + (1.0 - nu2 * nu2) / e2;
    1.0 / inv_e_star
}

/// CFL-based stable time step for DEM simulation.
///
/// Δt_CFL = min_i { r_i √(m_i / k) } × safety_factor
///
/// where k is the contact stiffness.  Safety factor = 0.2.
pub fn dem_cfl_timestep(particles: &[DemParticle], stiffness: f64) -> f64 {
    if particles.is_empty() || stiffness <= 0.0 {
        return 1.0e-6;
    }
    let safety = 0.2_f64;
    let dt_min = particles
        .iter()
        .map(|p| p.radius * (p.mass / stiffness).sqrt())
        .fold(f64::INFINITY, f64::min);
    safety * dt_min
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(pos: [f64; 3], radius: f64) -> DemParticle {
        DemParticle::sphere(pos, radius, 2500.0) // quartz density
    }

    // --- DemParticle ---

    #[test]
    fn test_sphere_mass() {
        let p = make_particle([0.0; 3], 0.01);
        let expected = 2500.0 * 4.0 / 3.0 * PI * 0.01_f64.powi(3);
        assert!((p.mass - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_sphere_moi() {
        let p = make_particle([0.0; 3], 0.01);
        let expected = 2.0 / 5.0 * p.mass * 0.01_f64 * 0.01;
        assert!((p.moment_of_inertia - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_kinetic_energy_zero_at_rest() {
        let p = make_particle([0.0; 3], 0.01);
        assert_eq!(p.kinetic_energy(), 0.0);
    }

    #[test]
    fn test_kinetic_energy_nonzero() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.velocity = [1.0, 0.0, 0.0];
        let expected = 0.5 * p.mass;
        assert!((p.kinetic_energy() - expected).abs() < 1e-12);
    }

    #[test]
    fn test_rotational_energy_nonzero() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.angular_velocity = [0.0, 0.0, 10.0];
        let expected = 0.5 * p.moment_of_inertia * 100.0;
        assert!((p.rotational_energy() - expected).abs() < 1e-12);
    }

    // --- HertzContact ---

    #[test]
    fn test_hertz_zero_overlap() {
        let h = HertzContact::new(70e9, 0.33);
        assert_eq!(h.compute_normal_force(0.0, 0.01, 0.01), 0.0);
    }

    #[test]
    fn test_hertz_negative_overlap() {
        let h = HertzContact::new(70e9, 0.33);
        assert_eq!(h.compute_normal_force(-1e-4, 0.01, 0.01), 0.0);
    }

    #[test]
    fn test_hertz_force_positive() {
        let h = HertzContact::new(70e9, 0.33);
        let f = h.compute_normal_force(1e-4, 0.01, 0.01);
        assert!(f > 0.0);
    }

    #[test]
    fn test_hertz_force_scales_with_overlap() {
        let h = HertzContact::new(70e9, 0.33);
        let f1 = h.compute_normal_force(1e-4, 0.01, 0.01);
        let f2 = h.compute_normal_force(2e-4, 0.01, 0.01);
        // F ∝ δ^(3/2): doubling δ should give 2^1.5 ≈ 2.83× force
        assert!(f2 > f1 * 2.0);
    }

    #[test]
    fn test_hertz_effective_modulus() {
        let h = HertzContact::new(70e9, 0.33);
        let e_star = h.effective_modulus();
        let expected = 70e9 / (2.0 * (1.0 - 0.33 * 0.33));
        assert!((e_star - expected).abs() / expected < 1e-10);
    }

    #[test]
    fn test_hertz_different_radii() {
        let h = HertzContact::new(70e9, 0.33);
        let f_equal = h.compute_normal_force(1e-4, 0.01, 0.01);
        let f_diff = h.compute_normal_force(1e-4, 0.005, 0.02);
        // Different R* → different force
        assert!((f_equal - f_diff).abs() > 1.0);
    }

    // --- hertz_stiffness ---

    #[test]
    fn test_hertz_stiffness_identical_materials() {
        let e_star = hertz_stiffness(70e9, 70e9, 0.33, 0.33);
        let expected = 1.0 / (2.0 * (1.0 - 0.33 * 0.33) / 70e9);
        assert!((e_star - expected).abs() / expected < 1e-9);
    }

    #[test]
    fn test_hertz_stiffness_positive() {
        let e_star = hertz_stiffness(210e9, 70e9, 0.3, 0.33);
        assert!(e_star > 0.0);
    }

    #[test]
    fn test_hertz_stiffness_hard_vs_soft() {
        let hard = hertz_stiffness(210e9, 210e9, 0.3, 0.3);
        let soft = hertz_stiffness(1e9, 1e9, 0.3, 0.3);
        assert!(hard > soft);
    }

    // --- dem_cfl_timestep ---

    #[test]
    fn test_cfl_empty_particles() {
        let dt = dem_cfl_timestep(&[], 1e6);
        assert_eq!(dt, 1e-6);
    }

    #[test]
    fn test_cfl_zero_stiffness() {
        let p = make_particle([0.0; 3], 0.01);
        let dt = dem_cfl_timestep(&[p], 0.0);
        assert_eq!(dt, 1e-6);
    }

    #[test]
    fn test_cfl_single_particle() {
        let p = make_particle([0.0; 3], 0.01);
        let k = 1e6;
        let dt = dem_cfl_timestep(std::slice::from_ref(&p), k);
        let expected = 0.2 * p.radius * (p.mass / k).sqrt();
        assert!((dt - expected).abs() / expected < 1e-9);
    }

    #[test]
    fn test_cfl_uses_minimum() {
        let p1 = make_particle([0.0; 3], 0.01);
        let p2 = make_particle([0.1; 3], 0.001); // smaller particle → smaller dt
        let k = 1e6;
        let dt_single = dem_cfl_timestep(std::slice::from_ref(&p2), k);
        let dt_both = dem_cfl_timestep(&[p1, p2], k);
        assert!((dt_both - dt_single).abs() / dt_single < 1e-6);
    }

    #[test]
    fn test_cfl_positive() {
        let p = make_particle([0.0; 3], 0.01);
        let dt = dem_cfl_timestep(&[p], 1e8);
        assert!(dt > 0.0);
    }

    // --- DemSimulation contact detection ---

    #[test]
    fn test_no_contact_far_apart() {
        let p1 = make_particle([0.0, 0.0, 0.0], 0.01);
        let p2 = make_particle([1.0, 0.0, 0.0], 0.01);
        let contact_model = HertzContact::new(70e9, 0.33);
        let mut sim = DemSimulation::new(contact_model, 0.3);
        sim.add_particle(p1);
        sim.add_particle(p2);
        let contacts = sim.detect_contacts();
        assert!(contacts.is_empty());
    }

    #[test]
    fn test_contact_overlapping() {
        let p1 = make_particle([0.0, 0.0, 0.0], 0.01);
        let p2 = make_particle([0.015, 0.0, 0.0], 0.01); // overlap = 0.005
        let contact_model = HertzContact::new(70e9, 0.33);
        let mut sim = DemSimulation::new(contact_model, 0.3);
        sim.add_particle(p1);
        sim.add_particle(p2);
        let contacts = sim.detect_contacts();
        assert_eq!(contacts.len(), 1);
        let (_, _, overlap, _) = contacts[0];
        assert!((overlap - 0.005).abs() < 1e-10);
    }

    #[test]
    fn test_contact_normal_direction() {
        let p1 = make_particle([0.0, 0.0, 0.0], 0.01);
        let p2 = make_particle([0.015, 0.0, 0.0], 0.01);
        let contact_model = HertzContact::new(70e9, 0.33);
        let mut sim = DemSimulation::new(contact_model, 0.3);
        sim.add_particle(p1);
        sim.add_particle(p2);
        let contacts = sim.detect_contacts();
        let (_, _, _, normal) = contacts[0];
        assert!((normal[0] - 1.0).abs() < 1e-10); // points from p1 to p2
    }

    #[test]
    fn test_multiple_contacts() {
        let p0 = make_particle([0.0, 0.0, 0.0], 0.01);
        let p1 = make_particle([0.015, 0.0, 0.0], 0.01);
        let p2 = make_particle([0.0, 0.015, 0.0], 0.01);
        let contact_model = HertzContact::new(70e9, 0.33);
        let mut sim = DemSimulation::new(contact_model, 0.3);
        sim.add_particle(p0);
        sim.add_particle(p1);
        sim.add_particle(p2);
        let contacts = sim.detect_contacts();
        assert_eq!(contacts.len(), 2);
    }

    // --- DemSimulation integration ---

    #[test]
    fn test_integrate_no_contact_linear_motion() {
        let mut p = make_particle([0.0; 3], 0.01);
        p.velocity = [1.0, 0.0, 0.0];
        let contact_model = HertzContact::new(70e9, 0.33);
        let mut sim = DemSimulation::new(contact_model, 0.3);
        sim.add_particle(p);
        sim.integrate(0.01);
        // No contacts: position should advance by v*dt = 0.01
        assert!((sim.particles[0].position[0] - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_integrate_contact_repulsion() {
        let p1 = make_particle([0.0, 0.0, 0.0], 0.01);
        let p2 = make_particle([0.015, 0.0, 0.0], 0.01); // overlap
        let contact_model = HertzContact::new(70e9, 0.33);
        let mut sim = DemSimulation::new(contact_model, 0.3);
        sim.add_particle(p1);
        sim.add_particle(p2);
        sim.integrate(1e-6);
        // p1 should have negative velocity (pushed left), p2 positive (pushed right)
        assert!(sim.particles[0].velocity[0] < 0.0);
        assert!(sim.particles[1].velocity[0] > 0.0);
    }

    // --- SphDemCoupling ---

    #[test]
    fn test_drag_zero_relative_velocity() {
        let coupling = SphDemCoupling::new(0.44, 1000.0);
        let p = make_particle([0.0; 3], 0.01);
        let f = coupling.compute_fluid_force(&p, [0.0; 3], 1000.0);
        // Drag should be zero; buoyancy acts in +z
        assert_eq!(f[0], 0.0);
        assert_eq!(f[1], 0.0);
        assert!(f[2] > 0.0);
    }

    #[test]
    fn test_drag_in_x_direction() {
        let coupling = SphDemCoupling::new(0.44, 1000.0);
        let p = make_particle([0.0; 3], 0.01);
        let f = coupling.compute_fluid_force(&p, [1.0, 0.0, 0.0], 1000.0);
        assert!(f[0] > 0.0); // drag pushes particle in flow direction
    }

    #[test]
    fn test_buoyancy_upward() {
        let coupling = SphDemCoupling::new(0.44, 1000.0);
        let p = make_particle([0.0; 3], 0.01);
        let f = coupling.compute_fluid_force(&p, [0.0; 3], 1000.0);
        let vol = 4.0 / 3.0 * PI * 0.01_f64.powi(3);
        let expected_b = 1000.0 * vol * 9.81;
        assert!((f[2] - expected_b).abs() / expected_b < 1e-9);
    }

    #[test]
    fn test_drag_scales_with_velocity_squared() {
        let coupling = SphDemCoupling::new(0.44, 1000.0);
        let p1 = make_particle([0.0; 3], 0.01);
        let p2 = make_particle([0.0; 3], 0.01);
        let f1 = coupling.compute_fluid_force(&p1, [1.0, 0.0, 0.0], 1000.0);
        let f2 = coupling.compute_fluid_force(&p2, [2.0, 0.0, 0.0], 1000.0);
        // At v_rel = 2 vs 1, drag should be 4× larger (v² law)
        let drag1 = f1[0]; // buoyancy is only in z
        let drag2 = f2[0];
        assert!((drag2 / drag1 - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_drag_zero_fluid_density() {
        let coupling = SphDemCoupling::new(0.44, 1000.0);
        let p = make_particle([0.0; 3], 0.01);
        let f = coupling.compute_fluid_force(&p, [5.0, 0.0, 0.0], 0.0);
        // Zero fluid density → zero drag, zero buoyancy
        assert_eq!(f[0], 0.0);
        assert_eq!(f[2], 0.0);
    }

    // --- ContactForce ---

    #[test]
    fn test_contact_force_newton_third_law() {
        let cf = ContactForce {
            normal_force: [10.0, 0.0, 0.0],
            tangential_force: [0.0, 2.0, 0.0],
            rolling_torque: [0.0; 3],
        };
        let f1 = cf.total_on_1();
        let f2 = cf.total_on_2();
        for k in 0..3 {
            assert!((f1[k] + f2[k]).abs() < 1e-12);
        }
    }

    // --- vec3 helpers ---

    #[test]
    fn test_vec3_cross_product() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = vec3_cross(a, b);
        assert!((c[0]).abs() < 1e-15);
        assert!((c[1]).abs() < 1e-15);
        assert!((c[2] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn test_vec3_norm() {
        let v = [3.0, 4.0, 0.0];
        assert!((vec3_norm(v) - 5.0).abs() < 1e-12);
    }
}
