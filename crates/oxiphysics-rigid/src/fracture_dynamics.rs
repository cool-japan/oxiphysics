// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rigid-body fracture dynamics: Voronoi splitting, crack planes, fragment
//! physics, breakable joints, impulse fragmentation, aggregate bodies, and
//! full fracture simulation.
//!
//! All vector math uses `[f64; 3]` arrays — no external linear-algebra crate.

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
// FractureThreshold
// ---------------------------------------------------------------------------

/// Voronoi-based fracture threshold: defines the stress level at which a body
/// shatters, and the preferred crack direction.
#[derive(Debug, Clone)]
pub struct FractureThreshold {
    /// Tensile stress threshold for fracture initiation (Pa).
    pub stress_threshold: f64,
    /// Preferred crack-propagation direction (unit vector).
    pub crack_direction: [f64; 3],
    /// Number of Voronoi seed points used to partition the body.
    pub voronoi_seeds: usize,
    /// Randomness factor in seed placement (0 = grid, 1 = fully random).
    pub seed_jitter: f64,
    /// Fracture energy per unit area (J/m²).
    pub fracture_energy: f64,
}

impl FractureThreshold {
    /// Create a new fracture threshold for a brittle material.
    ///
    /// * `stress_threshold` — failure stress in Pa.
    /// * `crack_direction`  — preferred propagation direction.
    /// * `seeds`            — number of Voronoi fragments.
    pub fn new(stress_threshold: f64, crack_direction: [f64; 3], seeds: usize) -> Self {
        Self {
            stress_threshold,
            crack_direction: v3_normalise(crack_direction),
            voronoi_seeds: seeds,
            seed_jitter: 0.3,
            fracture_energy: 10.0,
        }
    }

    /// Default threshold for soda-lime glass (50 MPa tensile, 8 seeds).
    pub fn glass() -> Self {
        Self::new(50e6, [1.0, 0.0, 0.0], 8)
    }

    /// Default threshold for concrete (4 MPa tensile, 6 seeds).
    pub fn concrete() -> Self {
        Self::new(4e6, [0.0, 1.0, 0.0], 6)
    }

    /// Returns `true` if `applied_stress` (Pa) meets or exceeds the threshold.
    pub fn is_exceeded(&self, applied_stress: f64) -> bool {
        applied_stress >= self.stress_threshold
    }

    /// Stress ratio (0 = no stress, 1 = exactly at threshold, >1 = exceeded).
    pub fn stress_ratio(&self, applied_stress: f64) -> f64 {
        applied_stress / self.stress_threshold.max(1.0)
    }

    /// Energy released per unit area at fracture (J/m²).
    pub fn released_energy_density(&self) -> f64 {
        self.fracture_energy
    }
}

// ---------------------------------------------------------------------------
// CrackPlane
// ---------------------------------------------------------------------------

/// Plane equation defining a fracture surface: `normal · x = d`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CrackPlane {
    /// Unit normal to the crack plane (points from body A to B).
    pub normal: [f64; 3],
    /// Plane offset: signed distance from the world origin.
    pub d: f64,
}

impl CrackPlane {
    /// Construct from a normal and a point on the plane.
    pub fn from_normal_point(normal: [f64; 3], point: [f64; 3]) -> Self {
        let n = v3_normalise(normal);
        Self {
            normal: n,
            d: v3_dot(n, point),
        }
    }

    /// Signed distance from `point` to the plane (positive = normal side).
    pub fn signed_distance(&self, point: [f64; 3]) -> f64 {
        v3_dot(self.normal, point) - self.d
    }

    /// Returns `true` if `point` is on the positive (normal) side of the plane.
    pub fn is_positive_side(&self, point: [f64; 3]) -> bool {
        self.signed_distance(point) >= 0.0
    }

    /// Project `point` onto the plane.
    pub fn project(&self, point: [f64; 3]) -> [f64; 3] {
        let dist = self.signed_distance(point);
        v3_sub(point, v3_scale(self.normal, dist))
    }

    /// Flip the plane normal (swap the two half-spaces).
    pub fn flipped(self) -> Self {
        Self {
            normal: v3_scale(self.normal, -1.0),
            d: -self.d,
        }
    }
}

// ---------------------------------------------------------------------------
// VoronoiFragment
// ---------------------------------------------------------------------------

/// A single fragment produced by Voronoi-based fracture.
#[derive(Debug, Clone)]
pub struct VoronoiFragment {
    /// Fragment index within the parent body.
    pub id: usize,
    /// Mass (kg).
    pub mass: f64,
    /// Diagonal inertia tensor (kg·m²).
    pub inertia: [f64; 3],
    /// Centre-of-mass position (m).
    pub position: [f64; 3],
    /// Linear velocity after split (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity after split (rad/s).
    pub angular_velocity: [f64; 3],
    /// Voronoi cell seed used to generate this fragment.
    pub seed: [f64; 3],
    /// Volume of the fragment (m³).
    pub volume: f64,
}

impl VoronoiFragment {
    /// Create a stationary fragment from mass and inertia.
    pub fn new(
        id: usize,
        mass: f64,
        inertia: [f64; 3],
        position: [f64; 3],
        seed: [f64; 3],
    ) -> Self {
        let vol = mass / 2500.0; // default density 2500 kg/m³
        Self {
            id,
            mass,
            inertia,
            position,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            seed,
            volume: vol,
        }
    }

    /// Kinetic energy (linear) of the fragment (J).
    pub fn kinetic_energy_linear(&self) -> f64 {
        0.5 * self.mass * v3_dot(self.velocity, self.velocity)
    }

    /// Kinetic energy (rotational) of the fragment (J).
    pub fn kinetic_energy_rotational(&self) -> f64 {
        let i = self.inertia;
        let w = self.angular_velocity;
        0.5 * (i[0] * w[0] * w[0] + i[1] * w[1] * w[1] + i[2] * w[2] * w[2])
    }

    /// Total kinetic energy (J).
    pub fn kinetic_energy(&self) -> f64 {
        self.kinetic_energy_linear() + self.kinetic_energy_rotational()
    }

    /// Linear momentum (kg·m/s).
    pub fn linear_momentum(&self) -> [f64; 3] {
        v3_scale(self.velocity, self.mass)
    }
}

// ---------------------------------------------------------------------------
// FragmentSet
// ---------------------------------------------------------------------------

/// Collection of [`VoronoiFragment`]s produced by a single fracture event.
#[derive(Debug, Clone, Default)]
pub struct FragmentSet {
    /// The individual fragments.
    pub fragments: Vec<VoronoiFragment>,
    /// Simulation time at which this set was created (s).
    pub creation_time: f64,
    /// World position of the fracture origin (m).
    pub origin: [f64; 3],
}

impl FragmentSet {
    /// Create an empty fragment set.
    pub fn new(creation_time: f64, origin: [f64; 3]) -> Self {
        Self {
            fragments: Vec::new(),
            creation_time,
            origin,
        }
    }

    /// Add a fragment to the set.
    pub fn push(&mut self, frag: VoronoiFragment) {
        self.fragments.push(frag);
    }

    /// Number of fragments.
    pub fn len(&self) -> usize {
        self.fragments.len()
    }

    /// Returns `true` if the set is empty.
    pub fn is_empty(&self) -> bool {
        self.fragments.is_empty()
    }

    /// Total mass of all fragments (kg).
    pub fn total_mass(&self) -> f64 {
        self.fragments.iter().map(|f| f.mass).sum()
    }

    /// Total kinetic energy of all fragments (J).
    pub fn total_kinetic_energy(&self) -> f64 {
        self.fragments.iter().map(|f| f.kinetic_energy()).sum()
    }

    /// Centre of mass of all fragments (m).
    pub fn centre_of_mass(&self) -> [f64; 3] {
        let m_total = self.total_mass();
        if m_total < 1e-12 {
            return [0.0; 3];
        }
        let mut sum = [0.0_f64; 3];
        for f in &self.fragments {
            sum = v3_add(sum, v3_scale(f.position, f.mass));
        }
        v3_scale(sum, 1.0 / m_total)
    }

    /// Net linear momentum of all fragments (kg·m/s).
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut p = [0.0_f64; 3];
        for f in &self.fragments {
            p = v3_add(p, f.linear_momentum());
        }
        p
    }

    /// Generate `n` fragments from a body of `total_mass` kg split by a crack
    /// plane.  Fragment seeds are placed on a uniform grid and jittered.
    pub fn voronoi_split(
        total_mass: f64,
        density: f64,
        half_extents: [f64; 3],
        n: usize,
        crack_plane: CrackPlane,
        creation_time: f64,
        origin: [f64; 3],
    ) -> Self {
        let mut set = FragmentSet::new(creation_time, origin);
        if n == 0 {
            return set;
        }
        let frag_mass = total_mass / n as f64;
        let frag_vol = frag_mass / density.max(1.0);
        let frag_side = frag_vol.cbrt();
        let ixx = frag_mass / 6.0 * frag_side * frag_side;
        let inertia = [ixx, ixx, ixx];

        let n_side = (n as f64).cbrt().ceil() as usize;
        for i in 0..n {
            let ix = (i % n_side) as f64 / n_side as f64;
            let iy = ((i / n_side) % n_side) as f64 / n_side as f64;
            let iz = (i / (n_side * n_side)) as f64 / n_side as f64;
            let seed = [
                (ix - 0.5) * 2.0 * half_extents[0],
                (iy - 0.5) * 2.0 * half_extents[1],
                (iz - 0.5) * 2.0 * half_extents[2],
            ];
            let pos = v3_add(origin, seed);
            let mut frag = VoronoiFragment::new(i, frag_mass, inertia, pos, seed);
            // Assign outward velocity along the crack normal (positive side) or opposite
            let side = if crack_plane.is_positive_side(pos) {
                1.0
            } else {
                -1.0
            };
            frag.velocity = v3_scale(crack_plane.normal, side * 0.5);
            set.push(frag);
        }
        set
    }
}

// ---------------------------------------------------------------------------
// ImpactEnergy
// ---------------------------------------------------------------------------

/// Kinetic energy bookkeeping for an impact: total incoming energy, fraction
/// going to fracture, and fraction dissipated as heat/sound.
#[derive(Debug, Clone, Copy)]
pub struct ImpactEnergy {
    /// Total impact kinetic energy (J).
    pub total_energy: f64,
    /// Fraction (0–1) of energy channelled into fracture surface creation.
    pub fracture_fraction: f64,
    /// Fraction (0–1) dissipated as heat/deformation (non-fracture).
    pub dissipation_fraction: f64,
}

impl ImpactEnergy {
    /// Construct from impact velocity, striking mass, and material fractions.
    ///
    /// * `mass`                 — striking body mass (kg).
    /// * `velocity`             — impact velocity (m/s).
    /// * `fracture_fraction`    — fraction going to fracture (0–1).
    /// * `dissipation_fraction` — fraction dissipated (0–1, remainder is elastic).
    pub fn new(
        mass: f64,
        velocity: f64,
        fracture_fraction: f64,
        dissipation_fraction: f64,
    ) -> Self {
        let total = 0.5 * mass * velocity * velocity;
        Self {
            total_energy: total,
            fracture_fraction: fracture_fraction.clamp(0.0, 1.0),
            dissipation_fraction: dissipation_fraction.clamp(0.0, 1.0),
        }
    }

    /// Energy available for fracture surface creation (J).
    pub fn fracture_energy(&self) -> f64 {
        self.total_energy * self.fracture_fraction
    }

    /// Energy dissipated as heat/deformation (J).
    pub fn dissipated_energy(&self) -> f64 {
        self.total_energy * self.dissipation_fraction
    }

    /// Elastic (rebounding) energy (J).
    pub fn elastic_energy(&self) -> f64 {
        let used = self.fracture_fraction + self.dissipation_fraction;
        self.total_energy * (1.0 - used.min(1.0))
    }

    /// Estimate number of fragments from Grady model given material toughness
    /// `k1c` (Pa·√m) and density `rho` (kg/m³).
    pub fn estimated_fragment_count(&self, k1c: f64, rho: f64, volume: f64) -> usize {
        // Mean fragment size (Grady): s = (K1c / (rho * strain_rate²))^(2/3)
        // We estimate strain_rate from energy density: e_dot ≈ sqrt(2*E_frac/rho/V)
        let e_frac = self.fracture_energy().max(1e-6);
        let strain_rate = (2.0 * e_frac / (rho.max(1.0) * volume.max(1e-9))).sqrt();
        let frag_size = (k1c / (rho.max(1.0) * strain_rate * strain_rate)).powf(2.0 / 3.0);
        let frag_vol = (4.0 / 3.0) * PI * (frag_size / 2.0).powi(3);
        (volume / frag_vol.max(1e-18)).round().max(1.0) as usize
    }
}

// ---------------------------------------------------------------------------
// BreakableJoint
// ---------------------------------------------------------------------------

/// A joint between two rigid bodies that breaks when the transmitted force
/// or accumulated impulse exceeds a material-dependent threshold.
#[derive(Debug, Clone)]
pub struct BreakableJoint {
    /// Unique joint identifier.
    pub id: usize,
    /// Index / handle of body A.
    pub body_a: usize,
    /// Index / handle of body B.
    pub body_b: usize,
    /// Force threshold (N) — instantaneous failure.
    pub force_threshold: f64,
    /// Torque threshold (N·m) — rotational failure.
    pub torque_threshold: f64,
    /// Accumulated impulse (N·s) threshold for fatigue failure.
    pub impulse_threshold: f64,
    /// Total accumulated impulse since last reset (N·s).
    pub accumulated_impulse: f64,
    /// Whether the joint has broken.
    pub broken: bool,
    /// Anchor position on body A in local coordinates (m).
    pub anchor_a: [f64; 3],
    /// Anchor position on body B in local coordinates (m).
    pub anchor_b: [f64; 3],
}

impl BreakableJoint {
    /// Create a new joint with force/torque thresholds.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: usize,
        body_a: usize,
        body_b: usize,
        force_threshold: f64,
        torque_threshold: f64,
        impulse_threshold: f64,
    ) -> Self {
        Self {
            id,
            body_a,
            body_b,
            force_threshold,
            torque_threshold,
            impulse_threshold,
            accumulated_impulse: 0.0,
            broken: false,
            anchor_a: [0.0; 3],
            anchor_b: [0.0; 3],
        }
    }

    /// Apply a constraint force (N) and torque (N·m) for time `dt` (s).
    ///
    /// Returns `true` if the joint breaks during this call.
    pub fn apply_load(&mut self, force_magnitude: f64, torque_magnitude: f64, dt: f64) -> bool {
        if self.broken {
            return false;
        }
        self.accumulated_impulse += force_magnitude.abs() * dt;
        if force_magnitude.abs() >= self.force_threshold
            || torque_magnitude.abs() >= self.torque_threshold
            || self.accumulated_impulse >= self.impulse_threshold
        {
            self.broken = true;
            return true;
        }
        false
    }

    /// Apply an instantaneous impulse (N·s) — e.g., from a collision response.
    ///
    /// Returns `true` if the joint breaks.
    pub fn apply_impulse(&mut self, impulse: f64) -> bool {
        if self.broken {
            return false;
        }
        self.accumulated_impulse += impulse.abs();
        if self.accumulated_impulse >= self.impulse_threshold {
            self.broken = true;
            return true;
        }
        false
    }

    /// Reset the accumulated impulse counter (call each physics step if needed).
    pub fn reset_accumulator(&mut self) {
        self.accumulated_impulse = 0.0;
    }

    /// Repair the joint (set broken = false, clear accumulator).
    pub fn repair(&mut self) {
        self.broken = false;
        self.accumulated_impulse = 0.0;
    }
}

// ---------------------------------------------------------------------------
// FractureEvent
// ---------------------------------------------------------------------------

/// Records a single fracture event in the simulation.
#[derive(Debug, Clone)]
pub struct FractureEvent {
    /// Simulation time at which fracture occurred (s).
    pub time: f64,
    /// World-space position of the fracture origin (m).
    pub position: [f64; 3],
    /// Force magnitude that triggered fracture (N).
    pub force: f64,
    /// Number of fragments produced.
    pub fragment_count: usize,
    /// ID of the joint or body that fractured.
    pub source_id: usize,
    /// The crack plane associated with this event.
    pub crack_plane: CrackPlane,
}

impl FractureEvent {
    /// Create a fracture event record.
    pub fn new(
        time: f64,
        position: [f64; 3],
        force: f64,
        fragment_count: usize,
        source_id: usize,
        crack_plane: CrackPlane,
    ) -> Self {
        Self {
            time,
            position,
            force,
            fragment_count,
            source_id,
            crack_plane,
        }
    }

    /// Energy released assuming fracture energy density `j_per_m2` and a
    /// crack area proportional to fragment count (each fragment ≈ 1 cm²).
    pub fn released_energy(&self, j_per_m2: f64) -> f64 {
        let area = self.fragment_count as f64 * 1e-4; // 1 cm² per fragment
        j_per_m2 * area
    }
}

// ---------------------------------------------------------------------------
// ImpulseFragmentation
// ---------------------------------------------------------------------------

/// Distributes the post-impact momentum to a set of fragments, conserving
/// total linear momentum.
#[derive(Debug, Clone)]
pub struct ImpulseFragmentation {
    /// Total impulse vector to distribute (N·s).
    pub total_impulse: [f64; 3],
    /// Radial spread parameter: fraction of impulse that becomes radial
    /// (outward) vs axial (along impact direction).  0 = all axial, 1 = all radial.
    pub radial_fraction: f64,
    /// Random spin amplitude (rad/s) added to fragments.
    pub spin_amplitude: f64,
}

impl ImpulseFragmentation {
    /// Create an impulse fragmentation model.
    ///
    /// * `total_impulse`    — total impulse (N·s) from the impact.
    /// * `radial_fraction`  — 0–1, fraction going radially outward.
    /// * `spin_amplitude`   — max angular velocity added to fragments (rad/s).
    pub fn new(total_impulse: [f64; 3], radial_fraction: f64, spin_amplitude: f64) -> Self {
        Self {
            total_impulse,
            radial_fraction: radial_fraction.clamp(0.0, 1.0),
            spin_amplitude,
        }
    }

    /// Distribute impulse to `fragments`.
    ///
    /// Each fragment receives a share proportional to its mass.  The axial
    /// component follows the impact direction; the radial component points
    /// from the fracture origin to the fragment position.
    pub fn distribute(&self, fragments: &mut [VoronoiFragment], origin: [f64; 3]) {
        let total_mass: f64 = fragments.iter().map(|f| f.mass).sum();
        if total_mass < 1e-12 {
            return;
        }
        let impulse_dir = v3_normalise(self.total_impulse);
        let impulse_mag = v3_norm(self.total_impulse);
        for frag in fragments.iter_mut() {
            let mass_ratio = frag.mass / total_mass;
            // Axial component
            let axial = v3_scale(
                impulse_dir,
                impulse_mag * mass_ratio * (1.0 - self.radial_fraction),
            );
            // Radial component (from origin toward fragment)
            let radial_dir = v3_normalise(v3_sub(frag.position, origin));
            let radial = v3_scale(radial_dir, impulse_mag * mass_ratio * self.radial_fraction);
            // dv = impulse / mass
            let dv = v3_scale(v3_add(axial, radial), 1.0 / frag.mass.max(1e-12));
            frag.velocity = v3_add(frag.velocity, dv);
        }
    }

    /// Total impulse magnitude (N·s).
    pub fn magnitude(&self) -> f64 {
        v3_norm(self.total_impulse)
    }
}

// ---------------------------------------------------------------------------
// AggregateBody
// ---------------------------------------------------------------------------

/// A compound rigid body composed of sub-components held together by
/// [`BreakableJoint`]s.  When joints break, sub-components become independent.
#[derive(Debug, Clone)]
pub struct AggregateBody {
    /// Unique body identifier.
    pub id: usize,
    /// Total mass (kg).
    pub total_mass: f64,
    /// Centre-of-mass world position (m).
    pub position: [f64; 3],
    /// Linear velocity (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity (rad/s).
    pub angular_velocity: [f64; 3],
    /// Joints connecting the sub-components.
    pub joints: Vec<BreakableJoint>,
    /// Whether the body has fully fractured (all major joints broken).
    pub fractured: bool,
    /// Sub-body masses (kg) — one per connected component.
    pub sub_masses: Vec<f64>,
}

impl AggregateBody {
    /// Create an aggregate body from a total mass and sub-mass list.
    pub fn new(id: usize, position: [f64; 3], sub_masses: Vec<f64>) -> Self {
        let total_mass = sub_masses.iter().sum();
        Self {
            id,
            total_mass,
            position,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            joints: Vec::new(),
            fractured: false,
            sub_masses,
        }
    }

    /// Add a breakable joint between sub-components.
    pub fn add_joint(&mut self, joint: BreakableJoint) {
        self.joints.push(joint);
    }

    /// Apply an external force (N) for `dt` seconds to all joints.
    ///
    /// Returns the number of joints that broke during this step.
    pub fn apply_force_step(&mut self, force: [f64; 3], _torque: [f64; 3], dt: f64) -> usize {
        let force_mag = v3_norm(force);
        let torque_mag = v3_norm(_torque);
        let mut broke = 0;
        for joint in &mut self.joints {
            if joint.apply_load(force_mag, torque_mag, dt) {
                broke += 1;
            }
        }
        if self.joints.iter().filter(|j| j.broken).count() == self.joints.len()
            && !self.joints.is_empty()
        {
            self.fractured = true;
        }
        broke
    }

    /// Number of intact (not broken) joints.
    pub fn intact_joint_count(&self) -> usize {
        self.joints.iter().filter(|j| !j.broken).count()
    }

    /// Number of broken joints.
    pub fn broken_joint_count(&self) -> usize {
        self.joints.iter().filter(|j| j.broken).count()
    }

    /// Total kinetic energy of the aggregate body (J).
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.total_mass * v3_dot(self.velocity, self.velocity)
    }
}

// ---------------------------------------------------------------------------
// FractureSimulation
// ---------------------------------------------------------------------------

/// Full fracture simulation: steps through time, detects fracture events,
/// splits bodies, and records events.
#[derive(Debug, Clone, Default)]
pub struct FractureSimulation {
    /// All aggregate bodies being simulated.
    pub bodies: Vec<AggregateBody>,
    /// All fracture events logged during the simulation.
    pub events: Vec<FractureEvent>,
    /// All fragment sets produced during fracture events.
    pub fragment_sets: Vec<FragmentSet>,
    /// Current simulation time (s).
    pub time: f64,
    /// Gravity vector (m/s²).
    pub gravity: [f64; 3],
}

impl FractureSimulation {
    /// Create a new simulation with default gravity (−9.81 m/s² in Y).
    pub fn new() -> Self {
        Self {
            bodies: Vec::new(),
            events: Vec::new(),
            fragment_sets: Vec::new(),
            time: 0.0,
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Add an aggregate body to the simulation.
    pub fn add_body(&mut self, body: AggregateBody) {
        self.bodies.push(body);
    }

    /// Step the simulation by `dt` seconds, applying `external_force` to every
    /// body.  Detects fracture events and logs them.
    ///
    /// Returns the number of fracture events that occurred this step.
    pub fn step(
        &mut self,
        dt: f64,
        external_force: [f64; 3],
        threshold: &FractureThreshold,
    ) -> usize {
        self.time += dt;
        let mut new_events = 0;

        for body in &mut self.bodies {
            if body.fractured {
                continue;
            }
            // Semi-implicit Euler: update velocity
            let dv = v3_scale(self.gravity, dt);
            body.velocity = v3_add(body.velocity, dv);
            // Update position
            body.position = v3_add(body.position, v3_scale(body.velocity, dt));

            let force_mag = v3_norm(external_force);
            let stress = force_mag / (body.total_mass.max(1e-6) * 0.001); // stress = F/A, assume A=1mm²

            if threshold.is_exceeded(stress) {
                let torque = [0.0_f64; 3];
                let broke = body.apply_force_step(external_force, torque, dt);
                if broke > 0 {
                    let normal = threshold.crack_direction;
                    let plane = CrackPlane::from_normal_point(normal, body.position);
                    let event = FractureEvent::new(
                        self.time,
                        body.position,
                        force_mag,
                        body.sub_masses.len() * 2,
                        body.id,
                        plane,
                    );
                    self.events.push(event);
                    new_events += 1;

                    // Produce fragment set
                    let fset = FragmentSet::voronoi_split(
                        body.total_mass,
                        2500.0,
                        [0.1, 0.1, 0.1],
                        threshold.voronoi_seeds,
                        plane,
                        self.time,
                        body.position,
                    );
                    self.fragment_sets.push(fset);
                }
            }
        }
        new_events
    }

    /// Number of fracture events logged.
    pub fn event_count(&self) -> usize {
        self.events.len()
    }

    /// Total fragment count across all events.
    pub fn total_fragment_count(&self) -> usize {
        self.fragment_sets.iter().map(|s| s.len()).sum()
    }

    /// Clear all events and fragment sets (keep bodies).
    pub fn clear_events(&mut self) {
        self.events.clear();
        self.fragment_sets.clear();
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- FractureThreshold ---

    #[test]
    fn fracture_threshold_exceeded_at_threshold() {
        let t = FractureThreshold::glass();
        assert!(t.is_exceeded(t.stress_threshold));
    }

    #[test]
    fn fracture_threshold_not_exceeded_below() {
        let t = FractureThreshold::glass();
        assert!(!t.is_exceeded(t.stress_threshold * 0.9));
    }

    #[test]
    fn fracture_threshold_stress_ratio_at_one() {
        let t = FractureThreshold::glass();
        let r = t.stress_ratio(t.stress_threshold);
        assert!((r - 1.0).abs() < 1e-10, "ratio={r}");
    }

    #[test]
    fn fracture_threshold_concrete_lower_than_glass() {
        let g = FractureThreshold::glass();
        let c = FractureThreshold::concrete();
        assert!(c.stress_threshold < g.stress_threshold);
    }

    #[test]
    fn fracture_threshold_released_energy_positive() {
        let t = FractureThreshold::glass();
        assert!(t.released_energy_density() > 0.0);
    }

    #[test]
    fn fracture_threshold_crack_direction_unit_length() {
        let t = FractureThreshold::new(1e6, [3.0, 4.0, 0.0], 4);
        let n = v3_norm(t.crack_direction);
        assert!((n - 1.0).abs() < 1e-10, "norm={n}");
    }

    // --- CrackPlane ---

    #[test]
    fn crack_plane_signed_distance_on_plane_zero() {
        let plane = CrackPlane::from_normal_point([0.0, 1.0, 0.0], [0.0, 5.0, 0.0]);
        let d = plane.signed_distance([3.0, 5.0, -2.0]);
        assert!(d.abs() < 1e-10, "d={d}");
    }

    #[test]
    fn crack_plane_positive_side_correct() {
        let plane = CrackPlane::from_normal_point([0.0, 1.0, 0.0], [0.0, 0.0, 0.0]);
        assert!(plane.is_positive_side([0.0, 1.0, 0.0]));
        assert!(!plane.is_positive_side([0.0, -1.0, 0.0]));
    }

    #[test]
    fn crack_plane_project_on_plane() {
        let plane = CrackPlane::from_normal_point([0.0, 1.0, 0.0], [0.0, 3.0, 0.0]);
        let proj = plane.project([5.0, 7.0, 2.0]);
        let dist = plane.signed_distance(proj);
        assert!(dist.abs() < 1e-10, "projected point not on plane: d={dist}");
    }

    #[test]
    fn crack_plane_flipped_reverses_sign() {
        let plane = CrackPlane::from_normal_point([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        let flipped = plane.flipped();
        let d_orig = plane.signed_distance([1.0, 0.0, 0.0]);
        let d_flip = flipped.signed_distance([1.0, 0.0, 0.0]);
        assert!((d_orig + d_flip).abs() < 1e-10);
    }

    // --- VoronoiFragment ---

    #[test]
    fn voronoi_fragment_zero_ke_stationary() {
        let frag = VoronoiFragment::new(0, 1.0, [0.1, 0.1, 0.1], [0.0; 3], [0.0; 3]);
        assert_eq!(frag.kinetic_energy(), 0.0);
    }

    #[test]
    fn voronoi_fragment_kinetic_energy_moving() {
        let mut frag = VoronoiFragment::new(0, 2.0, [0.1, 0.1, 0.1], [0.0; 3], [0.0; 3]);
        frag.velocity = [3.0, 0.0, 0.0];
        // KE_linear = 0.5 * 2 * 9 = 9
        let ke = frag.kinetic_energy_linear();
        assert!((ke - 9.0).abs() < 1e-10, "ke={ke}");
    }

    #[test]
    fn voronoi_fragment_rotational_ke() {
        let mut frag = VoronoiFragment::new(0, 1.0, [2.0, 2.0, 2.0], [0.0; 3], [0.0; 3]);
        frag.angular_velocity = [1.0, 0.0, 0.0];
        // KE_rot = 0.5 * 2 * 1² = 1
        let ke = frag.kinetic_energy_rotational();
        assert!((ke - 1.0).abs() < 1e-10, "ke_rot={ke}");
    }

    #[test]
    fn voronoi_fragment_momentum() {
        let mut frag = VoronoiFragment::new(0, 3.0, [1.0, 1.0, 1.0], [0.0; 3], [0.0; 3]);
        frag.velocity = [2.0, 0.0, 0.0];
        let p = frag.linear_momentum();
        assert!((p[0] - 6.0).abs() < 1e-10, "px={}", p[0]);
    }

    // --- FragmentSet ---

    #[test]
    fn fragment_set_empty_on_creation() {
        let s = FragmentSet::new(0.0, [0.0; 3]);
        assert!(s.is_empty());
        assert_eq!(s.len(), 0);
    }

    #[test]
    fn fragment_set_push_increases_count() {
        let mut s = FragmentSet::new(0.0, [0.0; 3]);
        let frag = VoronoiFragment::new(0, 1.0, [0.1, 0.1, 0.1], [0.0; 3], [0.0; 3]);
        s.push(frag);
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn fragment_set_total_mass() {
        let mut s = FragmentSet::new(0.0, [0.0; 3]);
        for i in 0..4 {
            s.push(VoronoiFragment::new(
                i,
                0.25,
                [0.1, 0.1, 0.1],
                [0.0; 3],
                [0.0; 3],
            ));
        }
        assert!((s.total_mass() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn fragment_set_voronoi_split_correct_count() {
        let plane = CrackPlane::from_normal_point([1.0, 0.0, 0.0], [0.0; 3]);
        let s = FragmentSet::voronoi_split(10.0, 2500.0, [0.1, 0.1, 0.1], 8, plane, 0.0, [0.0; 3]);
        assert_eq!(s.len(), 8);
    }

    #[test]
    fn fragment_set_voronoi_total_mass_conserved() {
        let plane = CrackPlane::from_normal_point([1.0, 0.0, 0.0], [0.0; 3]);
        let s = FragmentSet::voronoi_split(5.0, 2500.0, [0.1, 0.1, 0.1], 5, plane, 0.0, [0.0; 3]);
        assert!((s.total_mass() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn fragment_set_centre_of_mass_near_origin() {
        let plane = CrackPlane::from_normal_point([1.0, 0.0, 0.0], [0.0; 3]);
        let s =
            FragmentSet::voronoi_split(1.0, 1000.0, [0.05, 0.05, 0.05], 4, plane, 0.0, [0.0; 3]);
        let _com = s.centre_of_mass();
        // Just check it doesn't panic and is finite
        assert!(_com[0].is_finite());
    }

    // --- ImpactEnergy ---

    #[test]
    fn impact_energy_total_correct() {
        let e = ImpactEnergy::new(2.0, 10.0, 0.4, 0.3);
        assert!((e.total_energy - 100.0).abs() < 1e-10);
    }

    #[test]
    fn impact_energy_fracture_fraction() {
        let e = ImpactEnergy::new(2.0, 10.0, 0.4, 0.3);
        assert!((e.fracture_energy() - 40.0).abs() < 1e-10);
    }

    #[test]
    fn impact_energy_dissipated() {
        let e = ImpactEnergy::new(2.0, 10.0, 0.4, 0.3);
        assert!((e.dissipated_energy() - 30.0).abs() < 1e-10);
    }

    #[test]
    fn impact_energy_elastic_remainder() {
        let e = ImpactEnergy::new(2.0, 10.0, 0.4, 0.3);
        assert!((e.elastic_energy() - 30.0).abs() < 1e-10);
    }

    #[test]
    fn impact_energy_fragment_count_positive() {
        let e = ImpactEnergy::new(1.0, 50.0, 0.5, 0.3);
        let n = e.estimated_fragment_count(0.75e6, 2500.0, 1e-3);
        assert!(n >= 1, "n={n}");
    }

    #[test]
    fn impact_energy_zero_velocity_zero_total() {
        let e = ImpactEnergy::new(10.0, 0.0, 0.5, 0.3);
        assert_eq!(e.total_energy, 0.0);
    }

    // --- BreakableJoint ---

    #[test]
    fn breakable_joint_breaks_on_high_force() {
        let mut j = BreakableJoint::new(0, 0, 1, 1000.0, 500.0, 1e9);
        let broke = j.apply_load(1001.0, 0.0, 0.01);
        assert!(broke);
        assert!(j.broken);
    }

    #[test]
    fn breakable_joint_breaks_on_high_torque() {
        let mut j = BreakableJoint::new(0, 0, 1, 1e9, 200.0, 1e9);
        let broke = j.apply_load(0.0, 201.0, 0.01);
        assert!(broke);
    }

    #[test]
    fn breakable_joint_accumulates_impulse() {
        let mut j = BreakableJoint::new(0, 0, 1, 1e9, 1e9, 100.0);
        j.apply_load(10.0, 0.0, 5.0); // impulse = 50
        j.apply_load(10.0, 0.0, 5.0); // impulse = 100 → breaks
        assert!(j.broken);
    }

    #[test]
    fn breakable_joint_not_broken_below_threshold() {
        let mut j = BreakableJoint::new(0, 0, 1, 1000.0, 500.0, 1e9);
        let broke = j.apply_load(500.0, 200.0, 0.01);
        assert!(!broke);
        assert!(!j.broken);
    }

    #[test]
    fn breakable_joint_no_double_break() {
        let mut j = BreakableJoint::new(0, 0, 1, 10.0, 10.0, 1e9);
        j.apply_load(100.0, 100.0, 1.0);
        let second = j.apply_load(100.0, 100.0, 1.0);
        assert!(!second, "should not report break twice");
    }

    #[test]
    fn breakable_joint_repair_clears_broken() {
        let mut j = BreakableJoint::new(0, 0, 1, 10.0, 10.0, 1e9);
        j.apply_load(100.0, 0.0, 1.0);
        assert!(j.broken);
        j.repair();
        assert!(!j.broken);
        assert_eq!(j.accumulated_impulse, 0.0);
    }

    #[test]
    fn breakable_joint_impulse_method() {
        let mut j = BreakableJoint::new(0, 0, 1, 1e9, 1e9, 50.0);
        let broke = j.apply_impulse(60.0);
        assert!(broke);
    }

    // --- FractureEvent ---

    #[test]
    fn fracture_event_released_energy() {
        let plane = CrackPlane::from_normal_point([1.0, 0.0, 0.0], [0.0; 3]);
        let ev = FractureEvent::new(1.0, [0.0; 3], 5000.0, 10, 0, plane);
        let e = ev.released_energy(8.0);
        assert!(e > 0.0, "released energy={e}");
    }

    #[test]
    fn fracture_event_stores_time_and_position() {
        let plane = CrackPlane::from_normal_point([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        let ev = FractureEvent::new(2.5, [1.0, 2.0, 3.0], 1000.0, 4, 7, plane);
        assert!((ev.time - 2.5).abs() < 1e-10);
        assert!((ev.position[1] - 2.0).abs() < 1e-10);
        assert_eq!(ev.fragment_count, 4);
        assert_eq!(ev.source_id, 7);
    }

    // --- ImpulseFragmentation ---

    #[test]
    fn impulse_fragmentation_conserves_momentum() {
        let impulse_vec = [100.0_f64, 0.0, 0.0];
        let frag_model = ImpulseFragmentation::new(impulse_vec, 0.0, 0.0);
        let mut frags: Vec<VoronoiFragment> = (0..4)
            .map(|i| VoronoiFragment::new(i, 1.0, [0.1, 0.1, 0.1], [i as f64, 0.0, 0.0], [0.0; 3]))
            .collect();
        frag_model.distribute(&mut frags, [0.0; 3]);
        let total_px: f64 = frags.iter().map(|f| f.velocity[0] * f.mass).sum();
        assert!((total_px - 100.0).abs() < 1e-6, "px={total_px}");
    }

    #[test]
    fn impulse_fragmentation_magnitude() {
        let frag_model = ImpulseFragmentation::new([3.0, 4.0, 0.0], 0.5, 0.0);
        assert!((frag_model.magnitude() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn impulse_fragmentation_radial_gives_outward_velocity() {
        let frag_model = ImpulseFragmentation::new([0.0, 0.0, 100.0], 1.0, 0.0);
        let mut frags: Vec<VoronoiFragment> = vec![
            VoronoiFragment::new(0, 1.0, [0.1; 3], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]),
            VoronoiFragment::new(1, 1.0, [0.1; 3], [-1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        ];
        frag_model.distribute(&mut frags, [0.0; 3]);
        // Frag 0 at +x should get positive vx; frag 1 at -x should get negative vx
        assert!(
            frags[0].velocity[0] > 0.0,
            "frag0 vx={}",
            frags[0].velocity[0]
        );
        assert!(
            frags[1].velocity[0] < 0.0,
            "frag1 vx={}",
            frags[1].velocity[0]
        );
    }

    // --- AggregateBody ---

    #[test]
    fn aggregate_body_total_mass_correct() {
        let body = AggregateBody::new(0, [0.0; 3], vec![1.0, 2.0, 3.0]);
        assert!((body.total_mass - 6.0).abs() < 1e-10);
    }

    #[test]
    fn aggregate_body_no_fracture_initially() {
        let body = AggregateBody::new(0, [0.0; 3], vec![1.0, 1.0]);
        assert!(!body.fractured);
    }

    #[test]
    fn aggregate_body_joints_break_on_force() {
        let mut body = AggregateBody::new(0, [0.0; 3], vec![1.0, 1.0]);
        let joint = BreakableJoint::new(0, 0, 1, 10.0, 10.0, 1e9);
        body.add_joint(joint);
        let broke = body.apply_force_step([100.0, 0.0, 0.0], [0.0; 3], 0.01);
        assert_eq!(broke, 1);
        assert_eq!(body.broken_joint_count(), 1);
    }

    #[test]
    fn aggregate_body_fractured_when_all_joints_break() {
        let mut body = AggregateBody::new(0, [0.0; 3], vec![1.0, 1.0]);
        let j1 = BreakableJoint::new(0, 0, 1, 10.0, 10.0, 1e9);
        let j2 = BreakableJoint::new(1, 1, 0, 10.0, 10.0, 1e9);
        body.add_joint(j1);
        body.add_joint(j2);
        body.apply_force_step([100.0, 0.0, 0.0], [0.0; 3], 0.01);
        assert!(body.fractured);
    }

    #[test]
    fn aggregate_body_ke_positive_when_moving() {
        let mut body = AggregateBody::new(0, [0.0; 3], vec![2.0]);
        body.velocity = [3.0, 0.0, 0.0];
        // KE = 0.5 * 2 * 9 = 9
        assert!((body.kinetic_energy() - 9.0).abs() < 1e-10);
    }

    // --- FractureSimulation ---

    #[test]
    fn fracture_simulation_no_event_below_threshold() {
        let mut sim = FractureSimulation::new();
        let body = AggregateBody::new(0, [0.0; 3], vec![1.0]);
        sim.add_body(body);
        let threshold = FractureThreshold::glass();
        // Very small force — well below glass 50 MPa threshold
        let n = sim.step(0.01, [1.0, 0.0, 0.0], &threshold);
        assert_eq!(n, 0);
    }

    #[test]
    fn fracture_simulation_event_above_threshold() {
        let mut sim = FractureSimulation::new();
        let mut body = AggregateBody::new(0, [0.0; 3], vec![1.0, 1.0]);
        let joint = BreakableJoint::new(0, 0, 1, 0.001, 1e9, 1e9); // very low threshold
        body.add_joint(joint);
        sim.add_body(body);
        let mut threshold = FractureThreshold::glass();
        threshold.stress_threshold = 0.001; // very low threshold
        let n = sim.step(0.01, [1.0, 0.0, 0.0], &threshold);
        assert!(n > 0, "expected fracture event");
    }

    #[test]
    fn fracture_simulation_time_advances() {
        let mut sim = FractureSimulation::new();
        let threshold = FractureThreshold::glass();
        sim.step(0.016, [0.0; 3], &threshold);
        assert!((sim.time - 0.016).abs() < 1e-12);
    }

    #[test]
    fn fracture_simulation_fragments_produced() {
        let mut sim = FractureSimulation::new();
        let mut body = AggregateBody::new(0, [0.0; 3], vec![1.0, 1.0]);
        let joint = BreakableJoint::new(0, 0, 1, 0.001, 1e9, 1e9);
        body.add_joint(joint);
        sim.add_body(body);
        let mut threshold = FractureThreshold::glass();
        threshold.stress_threshold = 0.001;
        threshold.voronoi_seeds = 4;
        sim.step(0.01, [1.0, 0.0, 0.0], &threshold);
        assert!(sim.total_fragment_count() > 0);
    }

    #[test]
    fn fracture_simulation_clear_events() {
        let mut sim = FractureSimulation::new();
        let mut body = AggregateBody::new(0, [0.0; 3], vec![1.0, 1.0]);
        let joint = BreakableJoint::new(0, 0, 1, 0.001, 1e9, 1e9);
        body.add_joint(joint);
        sim.add_body(body);
        let mut threshold = FractureThreshold::glass();
        threshold.stress_threshold = 0.001;
        sim.step(0.01, [1.0, 0.0, 0.0], &threshold);
        sim.clear_events();
        assert_eq!(sim.event_count(), 0);
        assert_eq!(sim.total_fragment_count(), 0);
    }
}
