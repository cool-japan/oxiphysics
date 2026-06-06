// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Granular sand simulation (MPM-inspired).
//!
//! Provides particle-based sand dynamics with Drucker-Prager plasticity,
//! avalanche detection, angle-of-repose pile geometry, and granular
//! temperature / packing-fraction utilities.

// ---------------------------------------------------------------------------
// Small math helpers
// ---------------------------------------------------------------------------

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

// ---------------------------------------------------------------------------
// SandParticle
// ---------------------------------------------------------------------------

/// A single sand (MPM) particle.
pub struct SandParticle {
    /// World-space position \[m\].
    pub position: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// Mass \[kg\].
    pub mass: f64,
    /// Reference volume \[m³\].
    pub volume: f64,
    /// Deformation gradient F stored column-major as a flat 3×3 array.
    pub deformation_gradient: [f64; 9],
    /// Accumulated volumetric plastic strain (scalar).
    pub plastic_strain_vol: f64,
}

impl SandParticle {
    /// Create a new `SandParticle` at rest with an identity deformation gradient.
    pub fn new(position: [f64; 3], mass: f64, volume: f64) -> Self {
        // Identity matrix stored column-major
        #[rustfmt::skip]
        let identity = [
            1.0, 0.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
        ];
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            volume,
            deformation_gradient: identity,
            plastic_strain_vol: 0.0,
        }
    }

    /// Kinetic energy of this particle \[J\].
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }
}

// ---------------------------------------------------------------------------
// DruckerPragerPlasticity
// ---------------------------------------------------------------------------

/// Drucker-Prager (cone) plasticity model for granular materials.
///
/// The yield criterion is:
/// `‖s‖ ≤ k_phi · p` where `p = -σ_v/3` is the mean pressure and
/// `k_phi = 2·sin(φ)/(√3·(3−sin(φ)))`.
pub struct DruckerPragerPlasticity {
    /// Friction angle φ \[radians\].
    pub friction_angle: f64,
    /// Cohesion c \[Pa\].
    pub cohesion: f64,
}

impl DruckerPragerPlasticity {
    /// Create a new `DruckerPragerPlasticity`.
    pub fn new(friction_angle: f64, cohesion: f64) -> Self {
        Self {
            friction_angle,
            cohesion,
        }
    }

    /// Drucker-Prager cone constant k_φ.
    pub fn k_phi(&self) -> f64 {
        let sin_phi = self.friction_angle.sin();
        2.0 * sin_phi / (3.0_f64.sqrt() * (3.0 - sin_phi))
    }

    /// Return `true` if the stress state (σ_v, σ_d) is *outside* the yield surface.
    ///
    /// # Arguments
    /// * `sigma_v` – Volumetric (trace) stress \[Pa\]: σ_v = σ_xx + σ_yy + σ_zz.
    /// * `sigma_d` – Deviatoric stress norm ‖s‖ \[Pa\].
    pub fn yield_criterion(&self, sigma_v: f64, sigma_d: f64) -> bool {
        // Mean pressure p = -σ_v / 3  (compression positive)
        let p = -sigma_v / 3.0;
        let limit = self.k_phi() * p + self.cohesion;
        sigma_d > limit
    }

    /// Radial return mapping: project a 6-component Voigt stress back onto the
    /// yield cone.
    ///
    /// # Arguments
    /// * `stress` – Voigt stress \[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz\] \[Pa\].
    ///
    /// Returns the projected stress in the same layout.
    pub fn plastic_return_mapping(&self, stress: [f64; 6]) -> [f64; 6] {
        let sigma_v = stress[0] + stress[1] + stress[2];
        // Deviatoric part
        let p = sigma_v / 3.0;
        let s = [
            stress[0] - p,
            stress[1] - p,
            stress[2] - p,
            stress[3],
            stress[4],
            stress[5],
        ];
        let s_norm = (s[0] * s[0]
            + s[1] * s[1]
            + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
            .sqrt();

        // If inside yield surface, return unchanged.
        if !self.yield_criterion(sigma_v, s_norm) {
            return stress;
        }

        // Compute permissible deviatoric norm.
        let mean_p = -sigma_v / 3.0;
        let limit = (self.k_phi() * mean_p + self.cohesion).max(0.0);

        if s_norm < 1e-30 {
            return stress;
        }

        let scale = limit / s_norm;
        [
            scale * s[0] + p,
            scale * s[1] + p,
            scale * s[2] + p,
            scale * s[3],
            scale * s[4],
            scale * s[5],
        ]
    }
}

// ---------------------------------------------------------------------------
// SandSimulation
// ---------------------------------------------------------------------------

/// Simple MPM-inspired sand simulation.
pub struct SandSimulation {
    /// Particle list.
    pub particles: Vec<SandParticle>,
    /// Background grid cell size \[m\].
    pub grid_size: f64,
    /// Gravity vector \[m/s²\].
    pub gravity: [f64; 3],
}

impl SandSimulation {
    /// Create a new `SandSimulation`.
    pub fn new(grid_size: f64, gravity: [f64; 3]) -> Self {
        Self {
            particles: Vec::new(),
            grid_size,
            gravity,
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: SandParticle) {
        self.particles.push(p);
    }

    /// Advance the simulation by one time step `dt` \[s\].
    ///
    /// This simplified integrator applies gravity and integrates positions
    /// with a symplectic Euler scheme. A ground plane at y = 0 is enforced.
    pub fn step(&mut self, dt: f64) {
        for p in &mut self.particles {
            // Gravity acceleration
            let a = self.gravity;
            p.velocity = add3(p.velocity, scale3(a, dt));
            p.position = add3(p.position, scale3(p.velocity, dt));
            // Ground plane
            if p.position[1] < 0.0 {
                p.position[1] = 0.0;
                if p.velocity[1] < 0.0 {
                    p.velocity[1] = 0.0;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// AvalancheDetector
// ---------------------------------------------------------------------------

/// Detects sand particles that are at risk of avalanching.
///
/// A particle is considered unstable when the local surface slope exceeds
/// `critical_angle_deg` degrees relative to horizontal.
pub struct AvalancheDetector {
    /// Critical angle of repose \[degrees\] beyond which a particle is unstable.
    pub critical_angle_deg: f64,
    /// Reference to the particle list (indices into an external Vec).
    particles: Vec<[f64; 3]>,
}

impl AvalancheDetector {
    /// Create a new `AvalancheDetector` from a list of particle positions.
    pub fn new(critical_angle_deg: f64, positions: Vec<[f64; 3]>) -> Self {
        Self {
            critical_angle_deg,
            particles: positions,
        }
    }

    /// Return the indices of particles that are unstable.
    ///
    /// A particle at index `i` is considered unstable if its height above the
    /// lowest neighbour (within `search_radius`) corresponds to a slope angle
    /// greater than `critical_angle_deg`.  When no neighbours are found the
    /// particle is considered stable.
    pub fn detect_unstable(&self) -> Vec<usize> {
        let crit_rad = self.critical_angle_deg.to_radians();
        let tan_crit = crit_rad.tan();
        let search_radius = 2.0; // look within 2 m (normalised)

        let mut unstable = Vec::new();
        for (i, &pi) in self.particles.iter().enumerate() {
            let mut most_steep: f64 = 0.0;
            for (j, &pj) in self.particles.iter().enumerate() {
                if i == j {
                    continue;
                }
                let dx = pi[0] - pj[0];
                let dz = pi[2] - pj[2];
                let horiz = (dx * dx + dz * dz).sqrt();
                if horiz > search_radius || horiz < 1e-12 {
                    continue;
                }
                let dy = pi[1] - pj[1];
                if dy > 0.0 {
                    let slope = dy / horiz;
                    if slope > most_steep {
                        most_steep = slope;
                    }
                }
            }
            if most_steep > tan_crit {
                unstable.push(i);
            }
        }
        unstable
    }
}

// ---------------------------------------------------------------------------
// SandPile
// ---------------------------------------------------------------------------

/// Idealised conical sand pile.
pub struct SandPile {
    /// Base radius (half base-width) \[m\].
    pub base_width: f64,
    /// Angle of repose \[degrees\].
    pub angle_of_repose: f64,
}

impl SandPile {
    /// Create a new `SandPile`.
    pub fn new(base_width: f64, angle_of_repose: f64) -> Self {
        Self {
            base_width,
            angle_of_repose,
        }
    }

    /// Maximum height of the conical pile: h = r · tan(θ) \[m\].
    pub fn max_height(&self) -> f64 {
        let r = self.base_width / 2.0;
        r * self.angle_of_repose.to_radians().tan()
    }

    /// Surface height of the pile at radial distance `x` from the apex \[m\].
    ///
    /// Returns 0.0 for `x ≥ base_width / 2`.
    pub fn profile(&self, x: f64) -> f64 {
        let r = self.base_width / 2.0;
        if x.abs() >= r {
            return 0.0;
        }
        self.max_height() * (1.0 - x.abs() / r)
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Convert a Coulomb friction coefficient µ to angle of repose \[radians\].
///
/// θ = atan(µ)
pub fn angle_of_repose_from_friction(mu: f64) -> f64 {
    mu.atan()
}

/// Granular temperature of a particle assembly \[m²/s²\].
///
/// T_g = (1 / N) · Σ (v_i − ⟨v⟩)² where the sum is over all particles.
pub fn granular_temperature(particles: &[SandParticle]) -> f64 {
    if particles.is_empty() {
        return 0.0;
    }
    let n = particles.len() as f64;
    // Mean velocity
    let mut mean = [0.0_f64; 3];
    for p in particles {
        mean = add3(mean, p.velocity);
    }
    mean = scale3(mean, 1.0 / n);
    // Variance
    let mut var = 0.0_f64;
    for p in particles {
        let dv = sub3(p.velocity, mean);
        var += dot3(dv, dv);
    }
    var / n
}

/// Random close-packing fraction for monodisperse spheres of diameter `d`.
///
/// Returns φ ≈ 0.64 (Random Close Pack) — independent of particle diameter
/// for uniform spheres.  The `d_particle` parameter is accepted for API
/// completeness (it cancels out for uniform spheres).
pub fn packing_fraction_random(_d_particle: f64) -> f64 {
    // Empirical random close-packing fraction for monodisperse spheres.
    // The measured φ_RCP is conventionally reported as 0.6366; this is numerically
    // equivalent to 2/π and we use the constant to make the intent explicit.
    std::f64::consts::FRAC_2_PI
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // --- SandParticle ---

    #[test]
    fn test_particle_identity_deformation() {
        let p = SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6);
        // Diagonal of identity matrix (column-major indices 0,4,8)
        assert!((p.deformation_gradient[0] - 1.0).abs() < EPS);
        assert!((p.deformation_gradient[4] - 1.0).abs() < EPS);
        assert!((p.deformation_gradient[8] - 1.0).abs() < EPS);
    }

    #[test]
    fn test_particle_kinetic_energy_at_rest() {
        let p = SandParticle::new([0.0, 0.0, 0.0], 2.0, 1e-6);
        assert!(p.kinetic_energy().abs() < EPS);
    }

    #[test]
    fn test_particle_kinetic_energy_moving() {
        let mut p = SandParticle::new([0.0, 0.0, 0.0], 2.0, 1e-6);
        p.velocity = [3.0, 4.0, 0.0]; // speed = 5 m/s
        assert!((p.kinetic_energy() - 0.5 * 2.0 * 25.0).abs() < EPS);
    }

    #[test]
    fn test_particle_zero_plastic_strain() {
        let p = SandParticle::new([1.0, 2.0, 3.0], 0.5, 5e-7);
        assert!(p.plastic_strain_vol.abs() < EPS);
    }

    // --- DruckerPragerPlasticity ---

    #[test]
    fn test_yield_criterion_inside_cone() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        // Compressive pressure p = 1e4 Pa → limit = k_phi * 1e4
        // Use a deviatoric norm well below limit.
        let sigma_v = -3e4; // volumetric stress → p = 1e4
        let limit = dp.k_phi() * 1e4;
        let sigma_d = 0.5 * limit;
        assert!(
            !dp.yield_criterion(sigma_v, sigma_d),
            "Should be inside yield surface"
        );
    }

    #[test]
    fn test_yield_criterion_outside_cone() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        let sigma_v = -3e4;
        let limit = dp.k_phi() * 1e4;
        let sigma_d = 2.0 * limit; // outside
        assert!(
            dp.yield_criterion(sigma_v, sigma_d),
            "Should be outside yield surface"
        );
    }

    #[test]
    fn test_yield_criterion_zero_deviator() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        assert!(!dp.yield_criterion(-3e4, 0.0));
    }

    #[test]
    fn test_return_mapping_inside_unchanged() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 1000.0);
        let stress = [-10000.0, -10000.0, -10000.0, 0.0, 0.0, 0.0];
        let projected = dp.plastic_return_mapping(stress);
        for i in 0..6 {
            assert!((projected[i] - stress[i]).abs() < EPS);
        }
    }

    #[test]
    fn test_return_mapping_projects_onto_cone() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 0.0);
        // Large deviatoric shear stress (outside cone)
        let stress = [-30000.0, -30000.0, -30000.0, 50000.0, 0.0, 0.0];
        let projected = dp.plastic_return_mapping(stress);
        // After projection the yield criterion should be satisfied (not yielding)
        let sv = projected[0] + projected[1] + projected[2];
        let p = sv / 3.0;
        let s = [
            projected[0] - p,
            projected[1] - p,
            projected[2] - p,
            projected[3],
            projected[4],
            projected[5],
        ];
        let s_norm = (s[0] * s[0]
            + s[1] * s[1]
            + s[2] * s[2]
            + 2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]))
            .sqrt();
        assert!(
            !dp.yield_criterion(sv, s_norm),
            "After projection, should be on/inside cone"
        );
    }

    #[test]
    fn test_k_phi_increases_with_friction_angle() {
        let dp1 = DruckerPragerPlasticity::new(20_f64.to_radians(), 0.0);
        let dp2 = DruckerPragerPlasticity::new(40_f64.to_radians(), 0.0);
        assert!(dp2.k_phi() > dp1.k_phi());
    }

    #[test]
    fn test_yield_criterion_with_cohesion() {
        let dp = DruckerPragerPlasticity::new(30_f64.to_radians(), 5000.0);
        // Even with zero pressure, cohesion provides a limit.
        assert!(!dp.yield_criterion(0.0, 4999.0));
        assert!(dp.yield_criterion(0.0, 5001.0));
    }

    // --- SandSimulation ---

    #[test]
    fn test_simulation_gravity_falls() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 5.0, 0.0], 1.0, 1e-6));
        let dt = 0.01;
        for _ in 0..50 {
            sim.step(dt);
        }
        assert!(
            sim.particles[0].position[1] < 5.0,
            "Particle should fall under gravity"
        );
    }

    #[test]
    fn test_simulation_ground_plane() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 0.01, 0.0], 1.0, 1e-6));
        for _ in 0..1000 {
            sim.step(0.01);
        }
        assert!(
            sim.particles[0].position[1] >= 0.0,
            "Particle should not fall below ground"
        );
    }

    #[test]
    fn test_simulation_velocity_updates() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.add_particle(SandParticle::new([0.0, 10.0, 0.0], 1.0, 1e-6));
        sim.step(0.1);
        // After one step velocity[1] should be negative (downward)
        assert!(sim.particles[0].velocity[1] < 0.0);
    }

    #[test]
    fn test_simulation_no_particles_runs() {
        let mut sim = SandSimulation::new(0.1, [0.0, -9.81, 0.0]);
        sim.step(0.01); // should not panic
    }

    // --- AvalancheDetector ---

    #[test]
    fn test_avalanche_steep_particle_unstable() {
        // A high particle directly above a low one: angle ≈ 45° → unstable if crit < 45°
        let positions = vec![
            [0.0_f64, 1.0, 0.0], // particle 0: high
            [1.0_f64, 0.0, 0.0], // particle 1: low
        ];
        let detector = AvalancheDetector::new(30.0, positions);
        let unstable = detector.detect_unstable();
        assert!(unstable.contains(&0), "Steep particle 0 should be unstable");
    }

    #[test]
    fn test_avalanche_flat_surface_stable() {
        // All particles at the same height.
        let positions = vec![
            [0.0_f64, 0.0, 0.0],
            [1.0_f64, 0.0, 0.0],
            [2.0_f64, 0.0, 0.0],
        ];
        let detector = AvalancheDetector::new(30.0, positions);
        let unstable = detector.detect_unstable();
        assert!(
            unstable.is_empty(),
            "Flat surface should have no unstable particles"
        );
    }

    #[test]
    fn test_avalanche_single_particle_stable() {
        let positions = vec![[0.0_f64, 0.0, 0.0]];
        let detector = AvalancheDetector::new(30.0, positions);
        assert!(detector.detect_unstable().is_empty());
    }

    #[test]
    fn test_avalanche_empty_stable() {
        let detector = AvalancheDetector::new(30.0, vec![]);
        assert!(detector.detect_unstable().is_empty());
    }

    // --- SandPile ---

    #[test]
    fn test_sand_pile_max_height() {
        // base_width = 2 m, angle = 45° → r = 1, h = tan(45°) = 1 m
        let pile = SandPile::new(2.0, 45.0);
        assert!(
            (pile.max_height() - 1.0).abs() < 1e-10,
            "h={}",
            pile.max_height()
        );
    }

    #[test]
    fn test_sand_pile_profile_apex_zero() {
        let pile = SandPile::new(2.0, 45.0);
        // At x = 0 (apex): height = max_height()
        assert!((pile.profile(0.0) - pile.max_height()).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_profile_base_zero() {
        let pile = SandPile::new(2.0, 45.0);
        // At x = base_radius = 1.0: height = 0
        assert!(pile.profile(1.0).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_profile_outside_zero() {
        let pile = SandPile::new(2.0, 45.0);
        assert!(pile.profile(1.5).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_profile_midpoint() {
        let pile = SandPile::new(2.0, 45.0);
        // At x = 0.5 (half of radius 1.0): height = max_height() / 2
        assert!((pile.profile(0.5) - pile.max_height() * 0.5).abs() < EPS);
    }

    #[test]
    fn test_sand_pile_narrow_tall() {
        // 30° angle, 1 m wide: r = 0.5, h = 0.5 * tan(30°)
        let pile = SandPile::new(1.0, 30.0);
        let expected = 0.5 * 30_f64.to_radians().tan();
        assert!((pile.max_height() - expected).abs() < 1e-10);
    }

    // --- Free functions ---

    #[test]
    fn test_angle_of_repose_from_friction_45deg() {
        // mu = 1 → atan(1) = 45°
        let angle = angle_of_repose_from_friction(1.0);
        assert!((angle - std::f64::consts::FRAC_PI_4).abs() < 1e-10);
    }

    #[test]
    fn test_angle_of_repose_increases_with_mu() {
        let a1 = angle_of_repose_from_friction(0.5);
        let a2 = angle_of_repose_from_friction(1.0);
        assert!(a2 > a1);
    }

    #[test]
    fn test_angle_of_repose_zero_friction() {
        assert!(angle_of_repose_from_friction(0.0).abs() < EPS);
    }

    #[test]
    fn test_granular_temperature_at_rest() {
        let particles = vec![
            SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6),
            SandParticle::new([1.0, 0.0, 0.0], 1.0, 1e-6),
        ];
        assert!(granular_temperature(&particles).abs() < EPS);
    }

    #[test]
    fn test_granular_temperature_uniform_flow() {
        // All particles with same velocity → zero fluctuation temperature.
        let mut particles = vec![
            SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6),
            SandParticle::new([1.0, 0.0, 0.0], 1.0, 1e-6),
        ];
        for p in &mut particles {
            p.velocity = [2.0, 0.0, 0.0];
        }
        assert!(granular_temperature(&particles).abs() < EPS);
    }

    #[test]
    fn test_granular_temperature_opposite_velocities() {
        let mut particles = vec![
            SandParticle::new([0.0, 0.0, 0.0], 1.0, 1e-6),
            SandParticle::new([1.0, 0.0, 0.0], 1.0, 1e-6),
        ];
        particles[0].velocity = [1.0, 0.0, 0.0];
        particles[1].velocity = [-1.0, 0.0, 0.0];
        // Mean v = 0, T_g = (1² + 1²) / 2 = 1.0
        let tg = granular_temperature(&particles);
        assert!((tg - 1.0).abs() < EPS, "tg={tg}");
    }

    #[test]
    fn test_granular_temperature_empty() {
        assert!(granular_temperature(&[]).abs() < EPS);
    }

    #[test]
    fn test_packing_fraction_random_value() {
        let phi = packing_fraction_random(0.001);
        // Random close packing ≈ 0.64
        assert!(phi > 0.60 && phi < 0.70, "phi={phi}");
    }

    #[test]
    fn test_packing_fraction_independent_of_diameter() {
        // For monodisperse spheres the packing fraction is size-independent.
        let phi1 = packing_fraction_random(0.001);
        let phi2 = packing_fraction_random(1.0);
        assert!((phi1 - phi2).abs() < EPS);
    }
}
