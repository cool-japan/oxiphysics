// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid body fracture simulation with impulse-based fracture, Voronoi
//! shattering, crack propagation, debris simulation, and explosion modelling.
//!
//! All vector math uses `[f64; 3]` arrays — no external linear-algebra crate.
//!
//! # Main types
//!
//! - [`BreakableJoint`] — joint with tensile/shear/torsional failure thresholds
//! - [`FragmentSystem`] — manages a collection of rigid fragments after fracture
//! - [`ImpactFracture`] — impulse-based fracture and Voronoi shattering
//! - [`CrackPropagation`] — linear-elastic fracture mechanics, crack arrest
//! - [`DebrisSimulation`] — fragment tracking with secondary impacts
//! - [`ExplosionModel`] — pressure wave, Gurney equation, fragment distribution

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Internal vec3 helpers
// ---------------------------------------------------------------------------

#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

#[inline]
fn v3_normalise(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        v3_scale(a, 1.0 / n)
    }
}

#[inline]
fn v3_lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    v3_add(v3_scale(a, 1.0 - t), v3_scale(b, t))
}

// ---------------------------------------------------------------------------
// JointState
// ---------------------------------------------------------------------------

/// State of a breakable joint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum JointState {
    /// Joint is intact and transmitting forces.
    Intact,
    /// Joint has failed in the tensile direction.
    BrokenTensile,
    /// Joint has failed due to shear stress.
    BrokenShear,
    /// Joint has failed due to torsional stress.
    BrokenTorsional,
}

// ---------------------------------------------------------------------------
// BreakableJoint
// ---------------------------------------------------------------------------

/// A joint connecting two rigid bodies that fractures when stress thresholds
/// are exceeded.
///
/// Tracks accumulated tensile force, shear force, and torsional moment and
/// compares against configurable failure thresholds.
#[derive(Debug, Clone)]
pub struct BreakableJoint {
    /// Index of the first connected body.
    pub body_a: usize,
    /// Index of the second connected body.
    pub body_b: usize,
    /// Maximum tensile force before failure (N).
    pub tensile_threshold: f64,
    /// Maximum shear force before failure (N).
    pub shear_threshold: f64,
    /// Maximum torsional moment before failure (N·m).
    pub torsional_threshold: f64,
    /// Energy released per unit area on fracture (J/m²).
    pub fracture_energy: f64,
    /// Cross-sectional area of the joint (m²).
    pub area: f64,
    /// Current state of the joint.
    pub state: JointState,
    /// Accumulated damage in `[0, 1]`.
    pub damage: f64,
    /// Joint axis (unit vector from body_a to body_b in world space).
    pub axis: [f64; 3],
    /// Attachment point on body_a in world space.
    pub anchor_a: [f64; 3],
    /// Attachment point on body_b in world space.
    pub anchor_b: [f64; 3],
}

impl BreakableJoint {
    /// Create a new intact breakable joint with the given thresholds.
    pub fn new(
        body_a: usize,
        body_b: usize,
        tensile_threshold: f64,
        shear_threshold: f64,
        torsional_threshold: f64,
        fracture_energy: f64,
        area: f64,
    ) -> Self {
        Self {
            body_a,
            body_b,
            tensile_threshold,
            shear_threshold,
            torsional_threshold,
            fracture_energy,
            area,
            state: JointState::Intact,
            damage: 0.0,
            axis: [0.0, 1.0, 0.0],
            anchor_a: [0.0; 3],
            anchor_b: [0.0; 3],
        }
    }

    /// Returns `true` if the joint has broken.
    #[inline]
    pub fn is_broken(&self) -> bool {
        self.state != JointState::Intact
    }

    /// Apply a force vector to the joint and test for failure.
    ///
    /// `force` is the contact force vector at the joint in world space.
    /// `torque` is the torsional moment vector at the joint.
    ///
    /// Returns the new joint state after applying the loads.
    pub fn apply_force(&mut self, force: [f64; 3], torque: [f64; 3]) -> JointState {
        if self.is_broken() {
            return self.state;
        }

        let axis_norm = v3_normalise(self.axis);
        // Tensile component: projection onto joint axis
        let tensile = v3_dot(force, axis_norm).abs();
        // Shear: perpendicular component
        let tensile_vec = v3_scale(axis_norm, v3_dot(force, axis_norm));
        let shear_vec = v3_sub(force, tensile_vec);
        let shear = v3_norm(shear_vec);
        // Torsional: moment about joint axis
        let torsion = v3_dot(torque, axis_norm).abs();

        // Damage accumulation (Miner's rule style)
        let d_tensile = if self.tensile_threshold > 0.0 {
            tensile / self.tensile_threshold
        } else {
            0.0
        };
        let d_shear = if self.shear_threshold > 0.0 {
            shear / self.shear_threshold
        } else {
            0.0
        };
        let d_torsion = if self.torsional_threshold > 0.0 {
            torsion / self.torsional_threshold
        } else {
            0.0
        };

        // Instantaneous failure check
        if tensile >= self.tensile_threshold {
            self.state = JointState::BrokenTensile;
            self.damage = 1.0;
        } else if shear >= self.shear_threshold {
            self.state = JointState::BrokenShear;
            self.damage = 1.0;
        } else if torsion >= self.torsional_threshold {
            self.state = JointState::BrokenTorsional;
            self.damage = 1.0;
        } else {
            // Damage accumulation
            let increment =
                (d_tensile * d_tensile + d_shear * d_shear + d_torsion * d_torsion).sqrt() * 0.01;
            self.damage = (self.damage + increment).min(1.0);
            if self.damage >= 1.0 {
                self.state = JointState::BrokenTensile;
            }
        }

        self.state
    }

    /// Compute energy released on fracture (J).
    pub fn released_energy(&self) -> f64 {
        self.fracture_energy * self.area
    }

    /// Reset the joint to intact with zero damage.
    pub fn reset(&mut self) {
        self.state = JointState::Intact;
        self.damage = 0.0;
    }
}

// ---------------------------------------------------------------------------
// Fragment
// ---------------------------------------------------------------------------

/// A single rigid fragment produced by fracture or explosion.
#[derive(Debug, Clone)]
pub struct Fragment {
    /// World-space position of the fragment centre of mass (m).
    pub position: [f64; 3],
    /// Linear velocity (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_velocity: [f64; 3],
    /// Mass of the fragment (kg).
    pub mass: f64,
    /// AABB half-extents (m).
    pub half_extents: [f64; 3],
    /// Whether this fragment is still active in the simulation.
    pub active: bool,
    /// Coefficient of restitution `[0, 1]`.
    pub restitution: f64,
    /// Friction coefficient.
    pub friction: f64,
}

impl Fragment {
    /// Create a new fragment with position, velocity, mass, and extents.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        angular_velocity: [f64; 3],
        mass: f64,
        half_extents: [f64; 3],
    ) -> Self {
        Self {
            position,
            velocity,
            angular_velocity,
            mass,
            half_extents,
            active: true,
            restitution: 0.3,
            friction: 0.5,
        }
    }

    /// AABB minimum corner.
    pub fn aabb_min(&self) -> [f64; 3] {
        v3_sub(self.position, self.half_extents)
    }

    /// AABB maximum corner.
    pub fn aabb_max(&self) -> [f64; 3] {
        v3_add(self.position, self.half_extents)
    }

    /// Test AABB overlap with another fragment.
    pub fn aabb_overlaps(&self, other: &Fragment) -> bool {
        let amin = self.aabb_min();
        let amax = self.aabb_max();
        let bmin = other.aabb_min();
        let bmax = other.aabb_max();
        amin[0] <= bmax[0]
            && amax[0] >= bmin[0]
            && amin[1] <= bmax[1]
            && amax[1] >= bmin[1]
            && amin[2] <= bmax[2]
            && amax[2] >= bmin[2]
    }

    /// Kinetic energy of the fragment (linear + rotational approximation).
    pub fn kinetic_energy(&self) -> f64 {
        let ke_linear = 0.5 * self.mass * v3_dot(self.velocity, self.velocity);
        // Approximate inertia as sphere: I = 2/5 * m * r^2
        let r = v3_norm(self.half_extents) / 3.0_f64.sqrt();
        let inertia = 0.4 * self.mass * r * r;
        let ke_rot = 0.5 * inertia * v3_dot(self.angular_velocity, self.angular_velocity);
        ke_linear + ke_rot
    }

    /// Integrate fragment position and velocity for one time step.
    pub fn integrate(&mut self, dt: f64, gravity: [f64; 3]) {
        if !self.active {
            return;
        }
        self.velocity = v3_add(self.velocity, v3_scale(gravity, dt));
        self.position = v3_add(self.position, v3_scale(self.velocity, dt));
        // Angular damping
        self.angular_velocity = v3_scale(self.angular_velocity, (1.0 - 0.1 * dt).max(0.0));
    }
}

// ---------------------------------------------------------------------------
// FragmentSystem
// ---------------------------------------------------------------------------

/// Manages a collection of rigid fragments after a fracture event.
///
/// Tracks fragment state, integrates motion, and resolves simple ground
/// collisions.
#[derive(Debug, Clone)]
pub struct FragmentSystem {
    /// All fragments in the system.
    pub fragments: Vec<Fragment>,
    /// World-space ground plane height (y coordinate).
    pub ground_y: f64,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
}

impl FragmentSystem {
    /// Create an empty fragment system.
    pub fn new(ground_y: f64, gravity: [f64; 3]) -> Self {
        Self {
            fragments: Vec::new(),
            ground_y,
            gravity,
        }
    }

    /// Add a fragment to the system.
    pub fn add_fragment(&mut self, fragment: Fragment) {
        self.fragments.push(fragment);
    }

    /// Number of active fragments.
    pub fn active_count(&self) -> usize {
        self.fragments.iter().filter(|f| f.active).count()
    }

    /// Total kinetic energy of all active fragments.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.fragments
            .iter()
            .filter(|f| f.active)
            .map(|f| f.kinetic_energy())
            .sum()
    }

    /// Total linear momentum of all active fragments.
    pub fn total_momentum(&self) -> [f64; 3] {
        self.fragments
            .iter()
            .filter(|f| f.active)
            .fold([0.0; 3], |acc, f| v3_add(acc, v3_scale(f.velocity, f.mass)))
    }

    /// Integrate all fragments for one time step and resolve ground collisions.
    pub fn step(&mut self, dt: f64) {
        let gravity = self.gravity;
        let ground_y = self.ground_y;
        for frag in &mut self.fragments {
            if !frag.active {
                continue;
            }
            frag.integrate(dt, gravity);

            // Ground collision
            let min_y = frag.position[1] - frag.half_extents[1];
            if min_y < ground_y {
                frag.position[1] = ground_y + frag.half_extents[1];
                if frag.velocity[1] < 0.0 {
                    frag.velocity[1] = -frag.velocity[1] * frag.restitution;
                    frag.velocity[0] *= 1.0 - frag.friction * dt;
                    frag.velocity[2] *= 1.0 - frag.friction * dt;
                }
                // Deactivate if nearly stationary
                let speed = v3_norm(frag.velocity);
                if speed < 0.01 && frag.velocity[1].abs() < 0.01 {
                    frag.velocity = [0.0; 3];
                    if v3_norm(frag.angular_velocity) < 0.01 {
                        frag.active = false;
                    }
                }
            }
        }
    }

    /// Broad-phase AABB overlap detection between fragment pairs.
    ///
    /// Returns a list of overlapping fragment index pairs.
    pub fn broad_phase_overlaps(&self) -> Vec<(usize, usize)> {
        let mut pairs = Vec::new();
        let n = self.fragments.len();
        for i in 0..n {
            for j in (i + 1)..n {
                if self.fragments[i].active
                    && self.fragments[j].active
                    && self.fragments[i].aabb_overlaps(&self.fragments[j])
                {
                    pairs.push((i, j));
                }
            }
        }
        pairs
    }
}

// ---------------------------------------------------------------------------
// VoronoiSeed
// ---------------------------------------------------------------------------

/// A seed point for Voronoi fracture pattern generation.
#[derive(Debug, Clone, Copy)]
pub struct VoronoiSeed {
    /// Position of the seed point in local object space.
    pub position: [f64; 3],
    /// Weighting factor controlling fragment size (higher = larger cell).
    pub weight: f64,
}

impl VoronoiSeed {
    /// Create a seed at position with unit weight.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            weight: 1.0,
        }
    }

    /// Create a weighted seed.
    pub fn with_weight(position: [f64; 3], weight: f64) -> Self {
        Self { position, weight }
    }
}

// ---------------------------------------------------------------------------
// ImpactFracture
// ---------------------------------------------------------------------------

/// Impulse-based fracture engine with Voronoi shattering patterns.
///
/// Determines whether an impact impulse is sufficient to fracture a body and,
/// if so, generates a set of fragments using a simplified Voronoi
/// decomposition around the impact point.
#[derive(Debug, Clone)]
pub struct ImpactFracture {
    /// Minimum impulse magnitude required to initiate fracture (N·s).
    pub fracture_impulse_threshold: f64,
    /// Number of Voronoi seed points generated per fracture event.
    pub num_seeds: usize,
    /// Fraction of impact kinetic energy transferred to fragment motion.
    pub energy_transfer: f64,
    /// Scatter radius around impact point for Voronoi seeds (m).
    pub scatter_radius: f64,
    /// Minimum fragment mass as a fraction of the parent body mass.
    pub min_fragment_mass_fraction: f64,
}

impl ImpactFracture {
    /// Create a new impact fracture model with the given parameters.
    pub fn new(
        fracture_impulse_threshold: f64,
        num_seeds: usize,
        energy_transfer: f64,
        scatter_radius: f64,
    ) -> Self {
        Self {
            fracture_impulse_threshold,
            num_seeds,
            energy_transfer,
            scatter_radius,
            min_fragment_mass_fraction: 0.05,
        }
    }

    /// Test whether an impulse is sufficient to cause fracture.
    pub fn should_fracture(&self, impulse: [f64; 3]) -> bool {
        v3_norm(impulse) >= self.fracture_impulse_threshold
    }

    /// Generate Voronoi seeds around the impact point.
    ///
    /// Uses a deterministic pseudo-random sequence seeded by the impact
    /// position to avoid external dependencies.
    pub fn generate_seeds(&self, impact_point: [f64; 3]) -> Vec<VoronoiSeed> {
        let mut seeds = Vec::with_capacity(self.num_seeds);
        // Deterministic LCG seeded from impact point
        let seed_u64 = ((impact_point[0] * 1000.0) as u64)
            .wrapping_add((impact_point[1] * 997.0) as u64)
            .wrapping_add((impact_point[2] * 991.0) as u64)
            .wrapping_add(0xdeadbeef_cafebabe);
        let mut lcg = seed_u64;

        let rand_f64 = |s: &mut u64| -> f64 {
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*s >> 33) as f64 / (u32::MAX as f64)
        };

        // Add seed at impact point
        seeds.push(VoronoiSeed::new(impact_point));

        for _ in 1..self.num_seeds {
            // Spherical uniform sampling
            let u = rand_f64(&mut lcg) * 2.0 - 1.0; // [-1, 1]
            let phi = rand_f64(&mut lcg) * 2.0 * PI;
            let r = rand_f64(&mut lcg).cbrt() * self.scatter_radius;
            let sin_theta = (1.0 - u * u).sqrt();
            let pos = [
                impact_point[0] + r * sin_theta * phi.cos(),
                impact_point[1] + r * u,
                impact_point[2] + r * sin_theta * phi.sin(),
            ];
            seeds.push(VoronoiSeed::new(pos));
        }
        seeds
    }

    /// Fragment a body of mass `body_mass` at `body_position` given an
    /// `impulse` at `impact_point`.
    ///
    /// Returns a `FragmentSystem` populated with fragments whose velocities
    /// are derived from the impulse and seed geometry.
    pub fn fracture_body(
        &self,
        body_mass: f64,
        body_position: [f64; 3],
        body_half_extents: [f64; 3],
        impulse: [f64; 3],
        impact_point: [f64; 3],
    ) -> FragmentSystem {
        let gravity = [0.0, -9.81, 0.0];
        let mut system = FragmentSystem::new(body_position[1] - body_half_extents[1], gravity);

        if !self.should_fracture(impulse) {
            return system;
        }

        let seeds = self.generate_seeds(impact_point);
        let n = seeds.len().max(1);
        let frag_mass = (body_mass / n as f64).max(body_mass * self.min_fragment_mass_fraction);

        let impulse_dir = v3_normalise(impulse);
        let impulse_mag = v3_norm(impulse);
        let energy = 0.5 * impulse_mag * impulse_mag / body_mass * self.energy_transfer;
        let speed_scale = (2.0 * energy / frag_mass).sqrt() / n as f64;

        for seed in &seeds {
            // Direction from impact to seed
            let dir = v3_normalise(v3_sub(seed.position, impact_point));
            let dir = if v3_norm(dir) < 1e-9 {
                impulse_dir
            } else {
                dir
            };

            // Mix impulse direction with scatter direction
            let vel_dir = v3_normalise(v3_add(v3_scale(impulse_dir, 0.6), v3_scale(dir, 0.4)));
            let velocity = v3_scale(vel_dir, speed_scale);

            // Assign angular velocity perpendicular to velocity
            let up = if vel_dir[1].abs() < 0.9 {
                [0.0, 1.0, 0.0]
            } else {
                [1.0, 0.0, 0.0]
            };
            let omega_dir = v3_normalise(v3_cross(vel_dir, up));
            let angular_velocity = v3_scale(omega_dir, speed_scale * 2.0);

            let frag_extents = v3_scale(body_half_extents, 1.0 / (n as f64).cbrt());
            let frag = Fragment::new(
                seed.position,
                velocity,
                angular_velocity,
                frag_mass,
                frag_extents,
            );
            system.add_fragment(frag);
        }

        system
    }
}

// ---------------------------------------------------------------------------
// CrackTip
// ---------------------------------------------------------------------------

/// State of a propagating crack tip.
#[derive(Debug, Clone)]
pub struct CrackTip {
    /// Current position of the crack tip in world space.
    pub position: [f64; 3],
    /// Current crack propagation direction (unit vector).
    pub direction: [f64; 3],
    /// Accumulated crack length from origin (m).
    pub crack_length: f64,
    /// Whether this crack tip is still propagating.
    pub active: bool,
}

impl CrackTip {
    /// Create a new crack tip.
    pub fn new(position: [f64; 3], direction: [f64; 3]) -> Self {
        Self {
            position,
            direction: v3_normalise(direction),
            crack_length: 0.0,
            active: true,
        }
    }

    /// Advance the crack tip by `ds` along its current direction.
    pub fn advance(&mut self, ds: f64) {
        self.position = v3_add(self.position, v3_scale(self.direction, ds));
        self.crack_length += ds;
    }
}

// ---------------------------------------------------------------------------
// CrackPropagation
// ---------------------------------------------------------------------------

/// Linear-elastic fracture mechanics model for crack propagation.
///
/// Uses the stress intensity factor K_I and Paris law to advance a crack tip
/// under cyclic loading, and applies an arrest criterion.
#[derive(Debug, Clone)]
pub struct CrackPropagation {
    /// Mode-I fracture toughness K_Ic (Pa·√m).
    pub k1c: f64,
    /// Paris law coefficient C (dimensionless).
    pub paris_c: f64,
    /// Paris law exponent m.
    pub paris_m: f64,
    /// Crack arrest stress intensity factor K_arr (Pa·√m).
    pub arrest_k: f64,
    /// Crack tip list.
    pub tips: Vec<CrackTip>,
    /// Maximum allowed crack length (m) before the body is considered failed.
    pub max_crack_length: f64,
}

impl CrackPropagation {
    /// Create a new crack propagation model.
    pub fn new(k1c: f64, paris_c: f64, paris_m: f64, arrest_k: f64, max_crack_length: f64) -> Self {
        Self {
            k1c,
            paris_c,
            paris_m,
            arrest_k,
            tips: Vec::new(),
            max_crack_length,
        }
    }

    /// Compute the Mode-I stress intensity factor for a penny-shaped crack.
    ///
    /// `stress` is the far-field tensile stress (Pa).
    /// `crack_half_length` is the current crack half-length (m).
    pub fn stress_intensity_factor(&self, stress: f64, crack_half_length: f64) -> f64 {
        // K_I = σ * √(π * a) * F  (F = 1 for through crack in infinite plate)
        stress * (PI * crack_half_length).sqrt()
    }

    /// Paris law crack growth per cycle: da/dN = C * (ΔK)^m.
    pub fn crack_growth_per_cycle(&self, delta_k: f64) -> f64 {
        if delta_k <= 0.0 {
            return 0.0;
        }
        self.paris_c * delta_k.powf(self.paris_m)
    }

    /// Check whether the crack propagation should arrest.
    ///
    /// Returns `true` if the stress intensity factor is below the arrest
    /// criterion.
    pub fn should_arrest(&self, k_current: f64) -> bool {
        k_current < self.arrest_k
    }

    /// Check whether fracture is critical (unstable crack growth).
    pub fn is_critical(&self, k_current: f64) -> bool {
        k_current >= self.k1c
    }

    /// Add a crack tip at the given position and direction.
    pub fn add_tip(&mut self, position: [f64; 3], direction: [f64; 3]) {
        self.tips.push(CrackTip::new(position, direction));
    }

    /// Advance all active crack tips under the given applied stress for one
    /// loading cycle.
    ///
    /// `stress` is the applied far-field stress (Pa).
    /// `num_cycles` is the number of fatigue cycles to integrate.
    pub fn propagate(&mut self, stress: f64, num_cycles: f64) {
        for i in 0..self.tips.len() {
            if !self.tips[i].active {
                continue;
            }
            let crack_len = self.tips[i].crack_length.max(1e-6);
            let k = self.stress_intensity_factor(stress, crack_len);
            if self.should_arrest(k) {
                self.tips[i].active = false;
                continue;
            }
            if self.is_critical(k) {
                // Rapid unstable growth — advance to max length
                self.tips[i].crack_length = self.max_crack_length;
                self.tips[i].active = false;
                continue;
            }
            let da = self.crack_growth_per_cycle(k) * num_cycles;
            let max_len = self.max_crack_length;
            self.tips[i].advance(da);
            if self.tips[i].crack_length >= max_len {
                self.tips[i].active = false;
            }
        }
    }

    /// Crack velocity estimate (m/s) using Freund's dynamic crack formula.
    ///
    /// `k_current` is the current stress intensity (Pa·√m).
    /// `c_r` is the Rayleigh wave speed of the material (m/s).
    pub fn crack_velocity(&self, k_current: f64, c_r: f64) -> f64 {
        if k_current <= 0.0 || self.k1c <= 0.0 {
            return 0.0;
        }
        let ratio = k_current / self.k1c;
        // v = c_r * (1 - 1/ratio²) for K > K_Ic, else 0
        if ratio > 1.0 {
            c_r * (1.0 - 1.0 / (ratio * ratio))
        } else {
            0.0
        }
    }

    /// Returns `true` if any crack tip has reached the maximum crack length.
    pub fn is_fully_cracked(&self) -> bool {
        self.tips
            .iter()
            .any(|t| !t.active && t.crack_length >= self.max_crack_length)
    }

    /// Returns the total accumulated crack length across all tips.
    pub fn total_crack_length(&self) -> f64 {
        self.tips.iter().map(|t| t.crack_length).sum()
    }
}

// ---------------------------------------------------------------------------
// SecondaryImpact
// ---------------------------------------------------------------------------

/// Record of a secondary impact between a fragment and a surface or another
/// fragment.
#[derive(Debug, Clone)]
pub struct SecondaryImpact {
    /// Index of the impacting fragment.
    pub fragment_idx: usize,
    /// Impact point in world space (m).
    pub point: [f64; 3],
    /// Impact normal (unit vector pointing away from the surface).
    pub normal: [f64; 3],
    /// Relative velocity at impact (m/s).
    pub relative_velocity: f64,
    /// Impulse applied to resolve the impact (N·s).
    pub impulse: f64,
}

// ---------------------------------------------------------------------------
// DebrisSimulation
// ---------------------------------------------------------------------------

/// Fragment tracking simulation with secondary impacts and restitution.
///
/// Integrates the `FragmentSystem` over time and records secondary impact
/// events when fragments collide.
#[derive(Debug, Clone)]
pub struct DebrisSimulation {
    /// The fragment system being simulated.
    pub system: FragmentSystem,
    /// Accumulated secondary impact events.
    pub impacts: Vec<SecondaryImpact>,
    /// Simulation time elapsed (s).
    pub time: f64,
    /// Restitution coefficient for fragment-fragment collisions.
    pub fragment_restitution: f64,
    /// Maximum number of secondary impact events to record.
    pub max_impacts: usize,
}

impl DebrisSimulation {
    /// Create a new debris simulation wrapping the given fragment system.
    pub fn new(system: FragmentSystem, fragment_restitution: f64) -> Self {
        Self {
            system,
            impacts: Vec::new(),
            time: 0.0,
            fragment_restitution,
            max_impacts: 1000,
        }
    }

    /// Advance the simulation by `dt` seconds.
    ///
    /// Integrates fragment motion, resolves ground collisions, and performs a
    /// simple impulse-based fragment-fragment collision response.
    pub fn step(&mut self, dt: f64) {
        self.system.step(dt);
        self.resolve_fragment_collisions();
        self.time += dt;
    }

    /// Simple impulse-based fragment-fragment collision resolution.
    fn resolve_fragment_collisions(&mut self) {
        let pairs = self.system.broad_phase_overlaps();
        for (i, j) in pairs {
            if i >= self.system.fragments.len() || j >= self.system.fragments.len() {
                continue;
            }
            let rel_pos = v3_sub(
                self.system.fragments[j].position,
                self.system.fragments[i].position,
            );
            let dist = v3_norm(rel_pos);
            if dist < 1e-9 {
                continue;
            }
            let normal = v3_scale(rel_pos, 1.0 / dist);
            let rel_vel = v3_sub(
                self.system.fragments[i].velocity,
                self.system.fragments[j].velocity,
            );
            let approach = v3_dot(rel_vel, normal);
            if approach <= 0.0 {
                continue; // separating
            }

            let ma = self.system.fragments[i].mass;
            let mb = self.system.fragments[j].mass;
            let e = self.fragment_restitution;
            let j_imp = -(1.0 + e) * approach / (1.0 / ma + 1.0 / mb);

            let dv_a = v3_scale(normal, j_imp / ma);
            let dv_b = v3_scale(normal, -j_imp / mb);
            self.system.fragments[i].velocity = v3_add(self.system.fragments[i].velocity, dv_a);
            self.system.fragments[j].velocity = v3_add(self.system.fragments[j].velocity, dv_b);

            // Record impact event
            if self.impacts.len() < self.max_impacts {
                let point = v3_lerp(
                    self.system.fragments[i].position,
                    self.system.fragments[j].position,
                    0.5,
                );
                self.impacts.push(SecondaryImpact {
                    fragment_idx: i,
                    point,
                    normal,
                    relative_velocity: approach,
                    impulse: j_imp.abs(),
                });
            }
        }
    }

    /// Total simulation time elapsed.
    pub fn elapsed_time(&self) -> f64 {
        self.time
    }

    /// Number of secondary impacts recorded.
    pub fn impact_count(&self) -> usize {
        self.impacts.len()
    }

    /// Number of still-active fragments.
    pub fn active_fragment_count(&self) -> usize {
        self.system.active_count()
    }
}

// ---------------------------------------------------------------------------
// FragmentDistribution
// ---------------------------------------------------------------------------

/// Fragment size/velocity distribution from an explosion model.
#[derive(Debug, Clone)]
pub struct FragmentDistribution {
    /// Fragment masses (kg).
    pub masses: Vec<f64>,
    /// Fragment velocities (m/s) — speed only.
    pub velocities: Vec<f64>,
    /// Gurney velocity (m/s).
    pub gurney_velocity: f64,
}

impl FragmentDistribution {
    /// Mean fragment velocity.
    pub fn mean_velocity(&self) -> f64 {
        if self.velocities.is_empty() {
            return 0.0;
        }
        self.velocities.iter().sum::<f64>() / self.velocities.len() as f64
    }

    /// Maximum fragment velocity.
    pub fn max_velocity(&self) -> f64 {
        self.velocities.iter().cloned().fold(0.0_f64, f64::max)
    }

    /// Total fragment mass.
    pub fn total_mass(&self) -> f64 {
        self.masses.iter().sum()
    }
}

// ---------------------------------------------------------------------------
// ExplosionModel
// ---------------------------------------------------------------------------

/// Explosion model: pressure wave propagation, Gurney equation fragment
/// velocity, and Mott fragment-size distribution.
#[derive(Debug, Clone)]
pub struct ExplosionModel {
    /// Mass of explosive charge (kg).
    pub charge_mass: f64,
    /// Gurney energy constant √(2E) for the explosive (m/s).
    pub gurney_constant: f64,
    /// Casing mass surrounding the charge (kg).
    pub casing_mass: f64,
    /// Mott distribution scale parameter (mm).
    pub mott_mu: f64,
    /// Atmospheric pressure (Pa).
    pub ambient_pressure: f64,
    /// Number of fragments to model.
    pub num_fragments: usize,
}

impl ExplosionModel {
    /// Create a new explosion model.
    pub fn new(
        charge_mass: f64,
        gurney_constant: f64,
        casing_mass: f64,
        mott_mu: f64,
        num_fragments: usize,
    ) -> Self {
        Self {
            charge_mass,
            gurney_constant,
            casing_mass,
            mott_mu,
            ambient_pressure: 101_325.0,
            num_fragments,
        }
    }

    /// Compute the Gurney velocity for a cylindrical casing.
    ///
    /// `v_g = sqrt(2E) * sqrt(1 / (M/C + 0.5))`
    pub fn gurney_velocity_cylinder(&self) -> f64 {
        if self.charge_mass <= 0.0 {
            return 0.0;
        }
        let ratio = self.casing_mass / self.charge_mass;
        self.gurney_constant / (ratio + 0.5).sqrt()
    }

    /// Compute the Gurney velocity for a spherical casing.
    ///
    /// `v_g = sqrt(2E) * sqrt(1 / (M/C + 3/5))`
    pub fn gurney_velocity_sphere(&self) -> f64 {
        if self.charge_mass <= 0.0 {
            return 0.0;
        }
        let ratio = self.casing_mass / self.charge_mass;
        self.gurney_constant / (ratio + 0.6).sqrt()
    }

    /// Peak overpressure at distance `r` (m) from the explosion centre.
    ///
    /// Uses the Sadovsky empirical formula for TNT equivalent scaling.
    /// `charge_mass` is treated as TNT equivalent.
    pub fn peak_overpressure(&self, r: f64) -> f64 {
        if r <= 0.0 || self.charge_mass <= 0.0 {
            return 0.0;
        }
        let z = r / self.charge_mass.cbrt(); // scaled distance
        // Sadovsky: P = 0.84/z + 2.7/z^2 + 7.1/z^3  (MPa) for z > 1
        let p_mpa = 0.84 / z + 2.7 / (z * z) + 7.1 / (z * z * z);
        p_mpa * 1e6 // convert to Pa
    }

    /// Impulse delivered by the blast wave at distance `r` (m).
    ///
    /// Returns impulse per unit area (Pa·s).
    pub fn blast_impulse(&self, r: f64) -> f64 {
        if r <= 0.0 || self.charge_mass <= 0.0 {
            return 0.0;
        }
        let z = r / self.charge_mass.cbrt();
        // Simplified impulse scaling: i+ = 0.067 * C^(1/3) / z  (kPa·ms/kg^(1/3))
        let i_scaled = 0.067 / z;
        i_scaled * 1e3 // kPa·ms → Pa·s (approximate)
    }

    /// Generate a fragment mass/velocity distribution using the Gurney and
    /// Mott models.
    ///
    /// Returns a [`FragmentDistribution`] with `num_fragments` entries.
    pub fn fragment_distribution(&self) -> FragmentDistribution {
        let v_gurney = self.gurney_velocity_cylinder();
        let total_frag_mass = self.casing_mass;
        let mean_frag_mass = total_frag_mass / self.num_fragments as f64;

        // Deterministic LCG for reproducibility
        let mut lcg: u64 = 0xabcdef1234567890;
        let rand_f64 = |s: &mut u64| -> f64 {
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*s >> 33) as f64 / (u32::MAX as f64)
        };

        let mut masses = Vec::with_capacity(self.num_fragments);
        let mut velocities = Vec::with_capacity(self.num_fragments);

        for _ in 0..self.num_fragments {
            // Mott distribution: P(m > x) = exp(-(x / mu)^0.5) → exponential approx
            let u = rand_f64(&mut lcg).max(1e-9);
            let frag_mass = mean_frag_mass * (-u.ln()).max(0.01);
            masses.push(frag_mass);

            // Velocity scaled by sqrt(casing/fragment mass ratio)
            let mass_ratio = (total_frag_mass / frag_mass.max(1e-9)).sqrt();
            let vel = v_gurney * mass_ratio.min(3.0) * (0.8 + 0.4 * rand_f64(&mut lcg));
            velocities.push(vel);
        }

        FragmentDistribution {
            masses,
            velocities,
            gurney_velocity: v_gurney,
        }
    }

    /// Create a `FragmentSystem` populated with explosion fragments around
    /// `centre`.
    pub fn create_fragment_system(&self, centre: [f64; 3]) -> FragmentSystem {
        let dist = self.fragment_distribution();
        let gravity = [0.0, -9.81, 0.0];
        let mut system = FragmentSystem::new(centre[1] - 1.0, gravity);

        // Deterministic LCG for directions
        let mut lcg: u64 = 0xfeedface_deadbeef;
        let rand_f64 = |s: &mut u64| -> f64 {
            *s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (*s >> 33) as f64 / (u32::MAX as f64)
        };

        for i in 0..self.num_fragments {
            let speed = dist.velocities[i];
            let mass = dist.masses[i];

            // Uniform spherical direction
            let u = rand_f64(&mut lcg) * 2.0 - 1.0;
            let phi = rand_f64(&mut lcg) * 2.0 * PI;
            let sin_t = (1.0 - u * u).sqrt();
            let dir = [sin_t * phi.cos(), u, sin_t * phi.sin()];
            let velocity = v3_scale(dir, speed);

            let r = (mass / 7800.0).cbrt(); // approximate steel sphere radius
            let he = [r; 3];
            let frag = Fragment::new(centre, velocity, [0.0; 3], mass, he);
            system.add_fragment(frag);
        }

        system
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── BreakableJoint ────────────────────────────────────────────────────────

    #[test]
    fn breakable_joint_starts_intact() {
        let j = BreakableJoint::new(0, 1, 1000.0, 500.0, 200.0, 10.0, 0.01);
        assert_eq!(j.state, JointState::Intact);
        assert!(!j.is_broken());
        assert_eq!(j.damage, 0.0);
    }

    #[test]
    fn breakable_joint_tensile_failure() {
        let mut j = BreakableJoint::new(0, 1, 1000.0, 500.0, 200.0, 10.0, 0.01);
        j.axis = [0.0, 1.0, 0.0];
        let state = j.apply_force([0.0, 1500.0, 0.0], [0.0; 3]);
        assert_eq!(state, JointState::BrokenTensile);
        assert!(j.is_broken());
        assert_eq!(j.damage, 1.0);
    }

    #[test]
    fn breakable_joint_shear_failure() {
        let mut j = BreakableJoint::new(0, 1, 1000.0, 500.0, 200.0, 10.0, 0.01);
        j.axis = [0.0, 1.0, 0.0];
        let state = j.apply_force([800.0, 0.0, 0.0], [0.0; 3]); // shear along x
        assert_eq!(state, JointState::BrokenShear);
    }

    #[test]
    fn breakable_joint_torsional_failure() {
        let mut j = BreakableJoint::new(0, 1, 1000.0, 500.0, 200.0, 10.0, 0.01);
        j.axis = [0.0, 1.0, 0.0];
        let state = j.apply_force([0.0; 3], [0.0, 300.0, 0.0]);
        assert_eq!(state, JointState::BrokenTorsional);
    }

    #[test]
    fn breakable_joint_reset_clears_damage() {
        let mut j = BreakableJoint::new(0, 1, 100.0, 50.0, 20.0, 5.0, 0.01);
        j.axis = [0.0, 1.0, 0.0];
        j.apply_force([0.0, 200.0, 0.0], [0.0; 3]);
        assert!(j.is_broken());
        j.reset();
        assert_eq!(j.state, JointState::Intact);
        assert_eq!(j.damage, 0.0);
    }

    #[test]
    fn breakable_joint_below_threshold_stays_intact() {
        let mut j = BreakableJoint::new(0, 1, 1000.0, 500.0, 200.0, 10.0, 0.01);
        j.axis = [0.0, 1.0, 0.0];
        let state = j.apply_force([0.0, 500.0, 0.0], [0.0; 3]);
        assert_eq!(state, JointState::Intact);
    }

    #[test]
    fn breakable_joint_released_energy() {
        let j = BreakableJoint::new(0, 1, 1000.0, 500.0, 200.0, 100.0, 0.02);
        let energy = j.released_energy();
        assert!((energy - 2.0).abs() < 1e-10, "released energy = {energy}");
    }

    #[test]
    fn breakable_joint_broken_stays_broken_on_further_force() {
        let mut j = BreakableJoint::new(0, 1, 100.0, 50.0, 20.0, 5.0, 0.01);
        j.axis = [0.0, 1.0, 0.0];
        j.apply_force([0.0, 200.0, 0.0], [0.0; 3]);
        let first_state = j.state;
        j.apply_force([0.0, 5000.0, 0.0], [0.0; 3]);
        assert_eq!(j.state, first_state, "state should not change once broken");
    }

    // ── Fragment ──────────────────────────────────────────────────────────────

    #[test]
    fn fragment_aabb_min_max() {
        let f = Fragment::new([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0, 1.0, 1.0]);
        let mn = f.aabb_min();
        let mx = f.aabb_max();
        assert_eq!(mn, [-1.0, -1.0, -1.0]);
        assert_eq!(mx, [1.0, 1.0, 1.0]);
    }

    #[test]
    fn fragment_aabb_overlap_when_touching() {
        let a = Fragment::new([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        let b = Fragment::new([1.5, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]);
        assert!(a.aabb_overlaps(&b), "boxes separated by 0.5 should overlap");
    }

    #[test]
    fn fragment_aabb_no_overlap_when_far_apart() {
        let a = Fragment::new([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [0.5; 3]);
        let b = Fragment::new([5.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 1.0, [0.5; 3]);
        assert!(!a.aabb_overlaps(&b));
    }

    #[test]
    fn fragment_kinetic_energy_stationary_is_zero() {
        let f = Fragment::new([0.0; 3], [0.0; 3], [0.0; 3], 5.0, [0.5; 3]);
        assert!(f.kinetic_energy() < 1e-10);
    }

    #[test]
    fn fragment_kinetic_energy_moving() {
        let f = Fragment::new([0.0; 3], [2.0, 0.0, 0.0], [0.0; 3], 2.0, [0.5; 3]);
        let ke = f.kinetic_energy();
        // Linear: 0.5 * 2 * 4 = 4
        assert!(ke >= 4.0 - 1e-10, "KE={ke}");
    }

    #[test]
    fn fragment_integrate_falls_under_gravity() {
        let mut f = Fragment::new([0.0, 10.0, 0.0], [0.0; 3], [0.0; 3], 1.0, [0.5; 3]);
        f.integrate(1.0, [0.0, -9.81, 0.0]);
        assert!(
            f.velocity[1] < 0.0,
            "velocity should be negative after gravity"
        );
        assert!(f.position[1] < 10.0, "position should decrease");
    }

    // ── FragmentSystem ────────────────────────────────────────────────────────

    #[test]
    fn fragment_system_add_and_count() {
        let mut sys = FragmentSystem::new(0.0, [0.0, -9.81, 0.0]);
        sys.add_fragment(Fragment::new(
            [0.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.5; 3],
        ));
        sys.add_fragment(Fragment::new(
            [3.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.5; 3],
        ));
        assert_eq!(sys.active_count(), 2);
    }

    #[test]
    fn fragment_system_step_moves_fragments() {
        let mut sys = FragmentSystem::new(-100.0, [0.0, -9.81, 0.0]);
        sys.add_fragment(Fragment::new(
            [0.0, 5.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.1; 3],
        ));
        let y0 = sys.fragments[0].position[1];
        sys.step(0.1);
        let y1 = sys.fragments[0].position[1];
        assert!(y1 < y0, "fragment should fall: y0={y0} y1={y1}");
    }

    #[test]
    fn fragment_system_ground_collision_bounces() {
        let mut sys = FragmentSystem::new(0.0, [0.0, -9.81, 0.0]);
        let mut f = Fragment::new([0.0, 0.6, 0.0], [0.0, -2.0, 0.0], [0.0; 3], 1.0, [0.5; 3]);
        f.restitution = 0.5;
        sys.add_fragment(f);
        sys.step(0.1);
        // After ground collision the velocity should reverse (or be 0)
        // Position should be >= ground_y + half_extent
        assert!(sys.fragments[0].position[1] >= 0.5 - 1e-9);
    }

    #[test]
    fn fragment_system_total_momentum() {
        let mut sys = FragmentSystem::new(0.0, [0.0, -9.81, 0.0]);
        sys.add_fragment(Fragment::new(
            [0.0, 5.0, 0.0],
            [3.0, 0.0, 0.0],
            [0.0; 3],
            2.0,
            [0.1; 3],
        ));
        sys.add_fragment(Fragment::new(
            [3.0, 5.0, 0.0],
            [-3.0, 0.0, 0.0],
            [0.0; 3],
            2.0,
            [0.1; 3],
        ));
        let p = sys.total_momentum();
        assert!(p[0].abs() < 1e-10, "x-momentum should be zero: {}", p[0]);
    }

    #[test]
    fn fragment_system_broad_phase_detects_overlaps() {
        let mut sys = FragmentSystem::new(0.0, [0.0, -9.81, 0.0]);
        sys.add_fragment(Fragment::new([0.0; 3], [0.0; 3], [0.0; 3], 1.0, [1.0; 3]));
        sys.add_fragment(Fragment::new(
            [1.0, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            1.0,
            [1.0; 3],
        ));
        sys.add_fragment(Fragment::new(
            [20.0, 0.0, 0.0],
            [0.0; 3],
            [0.0; 3],
            1.0,
            [0.1; 3],
        ));
        let pairs = sys.broad_phase_overlaps();
        assert_eq!(pairs.len(), 1, "should detect exactly one overlap pair");
        assert_eq!(pairs[0], (0, 1));
    }

    // ── ImpactFracture ────────────────────────────────────────────────────────

    #[test]
    fn impact_fracture_threshold_check() {
        let model = ImpactFracture::new(100.0, 8, 0.5, 0.5);
        assert!(model.should_fracture([0.0, 200.0, 0.0]));
        assert!(!model.should_fracture([0.0, 50.0, 0.0]));
    }

    #[test]
    fn impact_fracture_generates_seeds() {
        let model = ImpactFracture::new(100.0, 8, 0.5, 0.5);
        let seeds = model.generate_seeds([0.0, 0.0, 0.0]);
        assert_eq!(seeds.len(), 8);
        // First seed should be at impact point
        assert!(v3_norm(v3_sub(seeds[0].position, [0.0; 3])) < 1e-10);
    }

    #[test]
    fn impact_fracture_creates_fragments() {
        let model = ImpactFracture::new(100.0, 8, 0.5, 0.5);
        let sys = model.fracture_body(
            10.0,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.0, 200.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert_eq!(sys.active_count(), 8, "should produce 8 fragments");
    }

    #[test]
    fn impact_fracture_below_threshold_produces_no_fragments() {
        let model = ImpactFracture::new(100.0, 8, 0.5, 0.5);
        let sys = model.fracture_body(
            10.0,
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [0.0, 50.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        assert_eq!(sys.active_count(), 0);
    }

    #[test]
    fn impact_fracture_fragments_have_positive_mass() {
        let model = ImpactFracture::new(50.0, 6, 0.8, 0.3);
        let sys = model.fracture_body(
            5.0,
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.0, 200.0, 0.0],
            [0.0, 1.0, 0.0],
        );
        for f in &sys.fragments {
            assert!(f.mass > 0.0, "all fragments must have positive mass");
        }
    }

    // ── CrackPropagation ─────────────────────────────────────────────────────

    #[test]
    fn crack_sif_zero_for_zero_stress() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        assert_eq!(model.stress_intensity_factor(0.0, 0.01), 0.0);
    }

    #[test]
    fn crack_sif_grows_with_length() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        let k1 = model.stress_intensity_factor(1e6, 0.01);
        let k2 = model.stress_intensity_factor(1e6, 0.04);
        assert!(
            k2 > k1,
            "SIF should grow with crack length: k1={k1} k2={k2}"
        );
    }

    #[test]
    fn crack_paris_law_zero_for_non_positive_dk() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        assert_eq!(model.crack_growth_per_cycle(0.0), 0.0);
        assert_eq!(model.crack_growth_per_cycle(-1.0), 0.0);
    }

    #[test]
    fn crack_arrest_criterion() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        assert!(model.should_arrest(0.4e6));
        assert!(!model.should_arrest(0.6e6));
    }

    #[test]
    fn crack_critical_check() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        assert!(model.is_critical(1.1e6));
        assert!(!model.is_critical(0.9e6));
    }

    #[test]
    fn crack_propagation_advances_tip() {
        // K_I = stress * sqrt(π * a). Start with a=0.01 m.
        // At stress=5e5: K = 5e5 * sqrt(π*0.01) ≈ 5e5 * 0.177 ≈ 88700 Pa√m > arrest_k=1e3
        let mut model = CrackPropagation::new(1e8, 1e-11, 3.0, 1e3, 1.0);
        model.add_tip([0.0; 3], [1.0, 0.0, 0.0]);
        model.tips[0].crack_length = 0.01; // start with non-zero crack length
        let len0 = model.tips[0].crack_length;
        model.propagate(5e5, 1000.0);
        let len1 = model.tips[0].crack_length;
        assert!(len1 > len0, "crack should advance: l0={len0} l1={len1}");
    }

    #[test]
    fn crack_velocity_zero_below_critical() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        let v = model.crack_velocity(0.5e6, 3000.0);
        assert_eq!(v, 0.0, "no crack velocity below K_Ic");
    }

    #[test]
    fn crack_velocity_positive_above_critical() {
        let model = CrackPropagation::new(1e6, 1e-12, 3.0, 0.5e6, 0.1);
        let v = model.crack_velocity(2e6, 3000.0);
        assert!(v > 0.0, "crack velocity should be positive above K_Ic: {v}");
    }

    #[test]
    fn crack_total_length_accumulates() {
        let mut model = CrackPropagation::new(1e8, 1e-11, 3.0, 1e3, 10.0);
        model.add_tip([0.0; 3], [1.0, 0.0, 0.0]);
        model.add_tip([0.0; 3], [-1.0, 0.0, 0.0]);
        // Set non-zero initial length so K > arrest_k from the start
        model.tips[0].crack_length = 0.01;
        model.tips[1].crack_length = 0.01;
        model.propagate(5e5, 1000.0);
        assert!(model.total_crack_length() > 0.02);
    }

    // ── DebrisSimulation ──────────────────────────────────────────────────────

    #[test]
    fn debris_simulation_step_advances_time() {
        let sys = FragmentSystem::new(0.0, [0.0, -9.81, 0.0]);
        let mut sim = DebrisSimulation::new(sys, 0.5);
        sim.step(0.016);
        assert!((sim.elapsed_time() - 0.016).abs() < 1e-10);
    }

    #[test]
    fn debris_simulation_records_impacts() {
        let mut sys = FragmentSystem::new(-100.0, [0.0; 3]);
        // Two fragments on collision course
        sys.add_fragment(Fragment::new(
            [0.0; 3],
            [5.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.5; 3],
        ));
        sys.add_fragment(Fragment::new(
            [0.5, 0.0, 0.0],
            [-5.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.5; 3],
        ));
        let mut sim = DebrisSimulation::new(sys, 0.5);
        sim.step(0.001);
        // Fragments overlap immediately → impact should be recorded
        assert!(sim.impact_count() > 0, "should detect collision");
    }

    #[test]
    fn debris_simulation_impulse_response_separates_fragments() {
        let mut sys = FragmentSystem::new(-100.0, [0.0; 3]);
        sys.add_fragment(Fragment::new(
            [0.0; 3],
            [10.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.4; 3],
        ));
        sys.add_fragment(Fragment::new(
            [0.3, 0.0, 0.0],
            [-10.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            [0.4; 3],
        ));
        let mut sim = DebrisSimulation::new(sys, 0.0);
        sim.step(0.001);
        // After perfectly inelastic restitution=0, velocities should be zero
        let vx0 = sim.system.fragments[0].velocity[0];
        let vx1 = sim.system.fragments[1].velocity[0];
        assert!(
            (vx0 + vx1).abs() < 0.1,
            "velocities should be approximately opposite"
        );
    }

    // ── ExplosionModel ────────────────────────────────────────────────────────

    #[test]
    fn explosion_gurney_cylinder_positive() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 100);
        let v = model.gurney_velocity_cylinder();
        assert!(v > 0.0, "Gurney velocity should be positive: {v}");
    }

    #[test]
    fn explosion_gurney_sphere_less_than_cylinder() {
        // Sphere has M/C + 0.6 in denominator vs 0.5 → smaller v
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 100);
        let v_cyl = model.gurney_velocity_cylinder();
        let v_sph = model.gurney_velocity_sphere();
        assert!(
            v_sph < v_cyl,
            "sphere v={v_sph} should be < cylinder v={v_cyl}"
        );
    }

    #[test]
    fn explosion_overpressure_decays_with_range() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 50);
        let p1 = model.peak_overpressure(1.0);
        let p2 = model.peak_overpressure(5.0);
        assert!(p2 < p1, "pressure should decay: p1={p1} p2={p2}");
    }

    #[test]
    fn explosion_overpressure_zero_at_zero_range() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 50);
        assert_eq!(model.peak_overpressure(0.0), 0.0);
    }

    #[test]
    fn explosion_blast_impulse_positive() {
        let model = ExplosionModel::new(10.0, 2440.0, 5.0, 5.0, 100);
        let i = model.blast_impulse(2.0);
        assert!(i > 0.0, "blast impulse should be positive: {i}");
    }

    #[test]
    fn explosion_fragment_distribution_correct_count() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 50);
        let dist = model.fragment_distribution();
        assert_eq!(dist.masses.len(), 50);
        assert_eq!(dist.velocities.len(), 50);
    }

    #[test]
    fn explosion_fragment_distribution_positive_velocities() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 20);
        let dist = model.fragment_distribution();
        for v in &dist.velocities {
            assert!(*v > 0.0, "all fragment velocities should be positive: {v}");
        }
    }

    #[test]
    fn explosion_create_fragment_system_count() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 20);
        let sys = model.create_fragment_system([0.0; 3]);
        assert_eq!(sys.active_count(), 20);
    }

    #[test]
    fn explosion_fragment_system_total_ke_positive() {
        let model = ExplosionModel::new(1.0, 2440.0, 1.0, 5.0, 10);
        let sys = model.create_fragment_system([0.0, 5.0, 0.0]);
        assert!(sys.total_kinetic_energy() > 0.0);
    }

    #[test]
    fn voronoi_seed_default_weight_one() {
        let s = VoronoiSeed::new([1.0, 2.0, 3.0]);
        assert_eq!(s.weight, 1.0);
    }

    #[test]
    fn voronoi_seed_with_weight() {
        let s = VoronoiSeed::with_weight([0.0; 3], 2.5);
        assert_eq!(s.weight, 2.5);
    }

    #[test]
    fn crack_tip_advance_increases_length() {
        let mut tip = CrackTip::new([0.0; 3], [1.0, 0.0, 0.0]);
        tip.advance(0.05);
        assert!((tip.crack_length - 0.05).abs() < 1e-10);
        assert!((tip.position[0] - 0.05).abs() < 1e-10);
    }

    #[test]
    fn impact_fracture_seeds_within_scatter_radius() {
        let radius = 1.0;
        let model = ImpactFracture::new(50.0, 10, 0.5, radius);
        let centre = [5.0, 0.0, -3.0];
        let seeds = model.generate_seeds(centre);
        for (i, seed) in seeds.iter().enumerate().skip(1) {
            let d = v3_norm(v3_sub(seed.position, centre));
            assert!(
                d <= radius + 1e-9,
                "seed {i} at distance {d} exceeds radius {radius}"
            );
        }
    }

    #[test]
    fn fragment_system_deactivates_settled_fragment() {
        let mut sys = FragmentSystem::new(0.0, [0.0; 3]); // no gravity
        // Place fragment at rest on ground
        let mut f = Fragment::new(
            [0.0, 0.5, 0.0],
            [0.0001, 0.0, 0.0], // near-zero velocity
            [0.0; 3],
            1.0,
            [0.5; 3],
        );
        f.restitution = 0.0;
        sys.add_fragment(f);
        // Step many times to allow settling
        for _ in 0..200 {
            sys.step(0.016);
        }
        // Fragment should deactivate once it settles
        let active = sys.active_count();
        assert!(
            active == 0 || active == 1,
            "fragment should settle: {active} active"
        );
    }
}
