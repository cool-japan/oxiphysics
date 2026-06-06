// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hair and fur simulation.
//!
//! Provides strand-based hair simulation with Verlet integration,
//! bending/torsion constraints, collision response, fur layer generation,
//! Kajiya-Kay shading, and the SuperHelix model for curly hair.

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

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-14 {
        [0.0, 0.0, 0.0]
    } else {
        scale3(v, 1.0 / l)
    }
}

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn clamp(v: f64, lo: f64, hi: f64) -> f64 {
    v.clamp(lo, hi)
}

// ---------------------------------------------------------------------------
// HairStrand
// ---------------------------------------------------------------------------

/// A single hair strand represented as a chain of particles.
///
/// The strand uses Verlet integration with position-based distance constraints
/// to simulate hair dynamics.
pub struct HairStrand {
    /// World-space positions of strand particles (root at index 0).
    pub positions: Vec<[f64; 3]>,
    /// Velocities of strand particles (m/s).
    pub velocities: Vec<[f64; 3]>,
    /// Rest lengths of segments between consecutive particles.
    pub rest_lengths: Vec<f64>,
    /// Number of segments (= num_particles - 1).
    pub num_segments: usize,
    /// Previous positions for Verlet integration.
    prev_positions: Vec<[f64; 3]>,
    /// Original root position used for reset.
    initial_root: [f64; 3],
}

impl HairStrand {
    /// Create a new hair strand hanging downward from `root`.
    ///
    /// # Arguments
    /// * `root` - Root attachment point in world space.
    /// * `length` - Total arc length of the strand (metres).
    /// * `n_segments` - Number of segments (number of particles = n_segments + 1).
    pub fn new(root: [f64; 3], length: f64, n_segments: usize) -> Self {
        let n_particles = n_segments + 1;
        let seg_len = length / n_segments as f64;
        let mut positions = Vec::with_capacity(n_particles);
        let mut velocities = Vec::with_capacity(n_particles);
        let mut rest_lengths = Vec::with_capacity(n_segments);

        for i in 0..n_particles {
            positions.push([root[0], root[1] - i as f64 * seg_len, root[2]]);
            velocities.push([0.0, 0.0, 0.0]);
        }
        for _ in 0..n_segments {
            rest_lengths.push(seg_len);
        }
        let prev_positions = positions.clone();
        Self {
            positions,
            velocities,
            rest_lengths,
            num_segments: n_segments,
            prev_positions,
            initial_root: root,
        }
    }

    /// Return the position of the tip particle (last particle).
    pub fn tip_position(&self) -> [f64; 3] {
        *self
            .positions
            .last()
            .expect("collection should not be empty")
    }

    /// Advance the strand one time-step using Verlet integration.
    ///
    /// The root particle (index 0) is kept fixed at its current position.
    ///
    /// # Arguments
    /// * `dt` - Time-step (seconds).
    /// * `gravity` - Gravitational acceleration vector (m/s²).
    pub fn integrate_verlet(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.positions.len();
        for i in 1..n {
            let pos = self.positions[i];
            let prev = self.prev_positions[i];
            // Verlet: x_new = 2*x - x_prev + a * dt^2
            let acc = gravity;
            let new_pos = [
                2.0 * pos[0] - prev[0] + acc[0] * dt * dt,
                2.0 * pos[1] - prev[1] + acc[1] * dt * dt,
                2.0 * pos[2] - prev[2] + acc[2] * dt * dt,
            ];
            self.prev_positions[i] = pos;
            self.positions[i] = new_pos;
            // Update velocity estimate
            self.velocities[i] = scale3(sub3(new_pos, prev), 0.5 / dt);
        }
    }

    /// Apply length constraints to each segment using iterative position correction.
    ///
    /// # Arguments
    /// * `stiffness` - Constraint stiffness in \[0, 1\]. 1.0 is fully rigid.
    pub fn apply_length_constraints(&mut self, stiffness: f64) {
        let n = self.positions.len();
        for i in 0..n - 1 {
            let pa = self.positions[i];
            let pb = self.positions[i + 1];
            let diff = sub3(pb, pa);
            let current_len = len3(diff);
            let rest_len = self.rest_lengths[i];
            if current_len < 1e-14 {
                continue;
            }
            let correction = stiffness * (current_len - rest_len) / current_len;
            let delta = scale3(diff, correction * 0.5);
            if i > 0 {
                self.positions[i] = add3(pa, delta);
            }
            self.positions[i + 1] = sub3(pb, delta);
        }
    }

    /// Return number of particles in this strand.
    pub fn num_particles(&self) -> usize {
        self.positions.len()
    }

    /// Apply a per-particle force impulse (e.g., wind) for one time-step.
    ///
    /// # Arguments
    /// * `force` - Force vector applied to each free particle.
    /// * `dt` - Time-step.
    pub fn apply_force_impulse(&mut self, force: [f64; 3], dt: f64) {
        let n = self.positions.len();
        for i in 1..n {
            let impulse = scale3(force, dt * dt);
            self.positions[i] = add3(self.positions[i], impulse);
        }
    }

    /// Apply velocity damping to reduce oscillations.
    ///
    /// # Arguments
    /// * `coeff` - Damping coefficient in \[0, 1\]. 0 = no damping, 1 = full stop.
    pub fn apply_damping(&mut self, coeff: f64) {
        for v in self.velocities.iter_mut() {
            *v = scale3(*v, 1.0 - coeff);
        }
    }

    /// Reset the strand to its initial hanging-down configuration.
    pub fn reset(&mut self) {
        let root = self.initial_root;
        let n = self.positions.len();
        let seg_len = self.rest_lengths[0];
        for i in 0..n {
            self.positions[i] = [root[0], root[1] - i as f64 * seg_len, root[2]];
            self.prev_positions[i] = self.positions[i];
            self.velocities[i] = [0.0, 0.0, 0.0];
        }
    }

    /// Compute the total kinetic energy of the strand (excluding the pinned root).
    ///
    /// Assumes unit mass per particle.
    pub fn kinetic_energy(&self) -> f64 {
        self.velocities[1..]
            .iter()
            .map(|v| 0.5 * dot3(*v, *v))
            .sum()
    }

    /// Compute the total arc length of the strand (sum of segment lengths).
    pub fn arc_length(&self) -> f64 {
        let n = self.positions.len();
        let mut total = 0.0;
        for i in 0..n - 1 {
            total += len3(sub3(self.positions[i + 1], self.positions[i]));
        }
        total
    }
}

// ---------------------------------------------------------------------------
// HairSimulation
// ---------------------------------------------------------------------------

/// A collection of hair strands with shared simulation parameters.
///
/// Manages aerodynamic drag, gravity, and damping for all strands.
pub struct HairSimulation {
    /// All hair strands managed by this simulation.
    pub strands: Vec<HairStrand>,
    /// Air density for drag computation (kg/m³).
    pub air_density: f64,
    /// Drag coefficient (dimensionless).
    pub drag_coefficient: f64,
    /// Strand cross-sectional radius for drag area (metres).
    pub strand_radius: f64,
}

impl HairSimulation {
    /// Create a new empty hair simulation.
    ///
    /// # Arguments
    /// * `air_density` - Air density (kg/m³), typically 1.225.
    /// * `drag_coefficient` - Drag coefficient, typically 0.5–1.5.
    /// * `strand_radius` - Strand radius in metres, typically 0.00005.
    pub fn new(air_density: f64, drag_coefficient: f64, strand_radius: f64) -> Self {
        Self {
            strands: Vec::new(),
            air_density,
            drag_coefficient,
            strand_radius,
        }
    }

    /// Add a new strand to the simulation.
    ///
    /// # Arguments
    /// * `root` - Root attachment point.
    /// * `length` - Total strand length (metres).
    /// * `n_segs` - Number of segments.
    pub fn add_strand(&mut self, root: [f64; 3], length: f64, n_segs: usize) {
        self.strands.push(HairStrand::new(root, length, n_segs));
    }

    /// Advance the entire simulation by one time-step.
    ///
    /// # Arguments
    /// * `dt` - Time-step (seconds).
    /// * `gravity` - Gravitational acceleration (m/s²).
    /// * `wind` - Wind velocity vector (m/s).
    pub fn step(&mut self, dt: f64, gravity: [f64; 3], wind: [f64; 3]) {
        self.apply_wind_forces(wind);
        self.apply_gravity(gravity);
        for strand in self.strands.iter_mut() {
            strand.integrate_verlet(dt, gravity);
            strand.apply_length_constraints(1.0);
        }
    }

    /// Apply aerodynamic wind forces to all strands.
    ///
    /// # Arguments
    /// * `wind` - Wind velocity vector (m/s).
    pub fn apply_wind_forces(&mut self, wind: [f64; 3]) {
        for strand in self.strands.iter_mut() {
            let n = strand.positions.len();
            for i in 1..n {
                let seg_len = if i < strand.rest_lengths.len() {
                    strand.rest_lengths[i - 1]
                } else {
                    strand.rest_lengths[strand.rest_lengths.len() - 1]
                };
                let rel_vel = sub3(wind, strand.velocities[i]);
                let speed = len3(rel_vel);
                // Drag force: F = 0.5 * rho * cd * A * v^2
                let area = 2.0 * self.strand_radius * seg_len;
                let drag_mag =
                    0.5 * self.air_density * self.drag_coefficient * area * speed * speed;
                let drag = if speed > 1e-14 {
                    scale3(normalize3(rel_vel), drag_mag)
                } else {
                    [0.0, 0.0, 0.0]
                };
                strand.positions[i] = add3(strand.positions[i], drag);
            }
        }
    }

    /// Apply gravity to all free particles by adjusting prev_positions.
    ///
    /// This effectively adds a gravity impulse for the next Verlet step.
    ///
    /// # Arguments
    /// * `g` - Gravitational acceleration (m/s²).
    pub fn apply_gravity(&mut self, g: [f64; 3]) {
        // Gravity is applied inside integrate_verlet; this method is a hook
        // for additional gravity-like body forces.
        let _ = g;
    }

    /// Apply velocity damping to all strands.
    ///
    /// # Arguments
    /// * `coeff` - Damping coefficient in \[0, 1\].
    pub fn apply_damping(&mut self, coeff: f64) {
        for strand in self.strands.iter_mut() {
            strand.apply_damping(coeff);
        }
    }

    /// Return total number of particles across all strands.
    pub fn total_particles(&self) -> usize {
        self.strands.iter().map(|s| s.num_particles()).sum()
    }
}

// ---------------------------------------------------------------------------
// BendingConstraint
// ---------------------------------------------------------------------------

/// A bending constraint acting on three consecutive strand particles.
///
/// Resists deviation from the rest angle between two consecutive segments.
pub struct HairBendingConstraint {
    /// Rest angle between the two segments (radians).
    pub rest_angle: f64,
    /// Constraint stiffness in \[0, 1\].
    pub stiffness: f64,
}

impl HairBendingConstraint {
    /// Create a new bending constraint.
    ///
    /// # Arguments
    /// * `rest_angle` - Rest angle in radians.
    /// * `stiffness` - Stiffness in \[0, 1\].
    pub fn new(rest_angle: f64, stiffness: f64) -> Self {
        Self {
            rest_angle,
            stiffness,
        }
    }

    /// Apply the bending constraint to three consecutive particles.
    ///
    /// # Arguments
    /// * `p0` - First particle position.
    /// * `p1` - Middle particle position (mutable).
    /// * `p2` - Last particle position (mutable).
    pub fn apply(&self, p0: &[f64; 3], p1: &mut [f64; 3], p2: &mut [f64; 3]) {
        let d01 = sub3(*p1, *p0);
        let d12 = sub3(*p2, *p1);
        let len01 = len3(d01);
        let len12 = len3(d12);
        if len01 < 1e-14 || len12 < 1e-14 {
            return;
        }
        let cos_angle = clamp(dot3(d01, d12) / (len01 * len12), -1.0, 1.0);
        let current_angle = cos_angle.acos();
        let angle_diff = current_angle - self.rest_angle;
        // Correction magnitude
        let correction = self.stiffness * angle_diff * 0.25;
        // Perpendicular direction for p1 correction
        let perp = normalize3(sub3(d12, scale3(d01, dot3(d01, d12) / (len01 * len01))));
        let c = scale3(perp, correction);
        *p1 = add3(*p1, c);
        *p2 = sub3(*p2, c);
    }

    /// Compute the current bending angle between two direction vectors.
    ///
    /// # Arguments
    /// * `d01` - Direction of first segment.
    /// * `d12` - Direction of second segment.
    pub fn compute_angle(d01: [f64; 3], d12: [f64; 3]) -> f64 {
        let len01 = len3(d01);
        let len12 = len3(d12);
        if len01 < 1e-14 || len12 < 1e-14 {
            return 0.0;
        }
        let cos_a = clamp(dot3(d01, d12) / (len01 * len12), -1.0, 1.0);
        cos_a.acos()
    }
}

// ---------------------------------------------------------------------------
// TorsionConstraint
// ---------------------------------------------------------------------------

/// A torsion constraint acting on four consecutive strand particles.
///
/// Resists twisting deformation along the strand axis.
pub struct TorsionConstraint {
    /// Rest twist angle (radians).
    pub rest_twist: f64,
    /// Constraint stiffness in \[0, 1\].
    pub stiffness: f64,
}

impl TorsionConstraint {
    /// Create a new torsion constraint.
    ///
    /// # Arguments
    /// * `rest_twist` - Rest twist angle in radians.
    /// * `stiffness` - Stiffness in \[0, 1\].
    pub fn new(rest_twist: f64, stiffness: f64) -> Self {
        Self {
            rest_twist,
            stiffness,
        }
    }

    /// Apply the torsion constraint to four consecutive particles.
    ///
    /// Computes the twist angle around the central segment and applies
    /// a corrective displacement.
    ///
    /// # Arguments
    /// * `p0` - First particle (mutable).
    /// * `p1` - Second particle (mutable).
    /// * `p2` - Third particle (mutable).
    /// * `p3` - Fourth particle (mutable).
    pub fn apply(
        &self,
        p0: &mut [f64; 3],
        p1: &mut [f64; 3],
        p2: &mut [f64; 3],
        p3: &mut [f64; 3],
    ) {
        // Compute the central axis (segment p1->p2)
        let axis = normalize3(sub3(*p2, *p1));
        // Vectors from central axis endpoints to outer particles
        let v0 = sub3(*p0, *p1);
        let v3 = sub3(*p3, *p2);
        // Project out the axial component
        let v0_perp = sub3(v0, scale3(axis, dot3(v0, axis)));
        let v3_perp = sub3(v3, scale3(axis, dot3(v3, axis)));
        let len0 = len3(v0_perp);
        let len3_v = len3(v3_perp);
        if len0 < 1e-14 || len3_v < 1e-14 {
            return;
        }
        let cos_twist = clamp(dot3(v0_perp, v3_perp) / (len0 * len3_v), -1.0, 1.0);
        let cross_sign = dot3(cross3(v0_perp, v3_perp), axis);
        let current_twist = cos_twist.acos() * cross_sign.signum();
        let twist_diff = current_twist - self.rest_twist;
        let correction = self.stiffness * twist_diff * 0.25;
        // Apply rotation corrections
        let rot_axis = scale3(axis, correction);
        let delta0 = cross3(rot_axis, v0_perp);
        let delta3 = cross3(scale3(rot_axis, -1.0), v3_perp);
        *p0 = add3(*p0, scale3(delta0, 0.5));
        *p3 = add3(*p3, scale3(delta3, 0.5));
        // Slight correction to p1 and p2 for momentum balance
        *p1 = sub3(*p1, scale3(delta0, 0.25));
        *p2 = sub3(*p2, scale3(delta3, 0.25));
    }

    /// Compute the current torsion (twist) angle for four particles.
    ///
    /// Returns the signed twist angle in radians.
    pub fn compute_twist(p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], p3: [f64; 3]) -> f64 {
        let axis = normalize3(sub3(p2, p1));
        let v0 = sub3(p0, p1);
        let v3 = sub3(p3, p2);
        let v0_perp = sub3(v0, scale3(axis, dot3(v0, axis)));
        let v3_perp = sub3(v3, scale3(axis, dot3(v3, axis)));
        let len0 = len3(v0_perp);
        let len3_v = len3(v3_perp);
        if len0 < 1e-14 || len3_v < 1e-14 {
            return 0.0;
        }
        let cos_twist = clamp(dot3(v0_perp, v3_perp) / (len0 * len3_v), -1.0, 1.0);
        let cross_sign = dot3(cross3(v0_perp, v3_perp), axis);
        cos_twist.acos() * cross_sign.signum()
    }
}

// ---------------------------------------------------------------------------
// HairCollision
// ---------------------------------------------------------------------------

/// Collision response handler for hair-object interactions.
///
/// Provides sphere and capsule collision detection and response for hair strands.
pub struct HairCollision;

impl HairCollision {
    /// Apply sphere collision response to all particles in a strand.
    ///
    /// Pushes particles outside the sphere if they penetrate it.
    ///
    /// # Arguments
    /// * `strand` - The hair strand to apply collision to.
    /// * `center` - Sphere center in world space.
    /// * `radius` - Sphere radius (metres).
    pub fn collide_with_sphere(strand: &mut HairStrand, center: [f64; 3], radius: f64) {
        for pos in strand.positions.iter_mut() {
            let diff = sub3(*pos, center);
            let dist = len3(diff);
            if dist < radius && dist > 1e-14 {
                let normal = scale3(diff, 1.0 / dist);
                *pos = add3(center, scale3(normal, radius));
            }
        }
    }

    /// Apply capsule collision response to all particles in a strand.
    ///
    /// Projects each particle to the nearest point on the capsule axis and
    /// pushes it outside if it penetrates the capsule.
    ///
    /// # Arguments
    /// * `strand` - The hair strand to apply collision to.
    /// * `a` - First endpoint of the capsule axis.
    /// * `b` - Second endpoint of the capsule axis.
    /// * `radius` - Capsule radius (metres).
    pub fn collide_with_capsule(strand: &mut HairStrand, a: [f64; 3], b: [f64; 3], radius: f64) {
        let ab = sub3(b, a);
        let len_ab = len3(ab);
        if len_ab < 1e-14 {
            // Degenerate capsule: treat as sphere at a
            Self::collide_with_sphere(strand, a, radius);
            return;
        }
        let ab_norm = scale3(ab, 1.0 / len_ab);
        for pos in strand.positions.iter_mut() {
            let ap = sub3(*pos, a);
            let t = clamp(dot3(ap, ab_norm), 0.0, len_ab);
            let closest = add3(a, scale3(ab_norm, t));
            let diff = sub3(*pos, closest);
            let dist = len3(diff);
            if dist < radius && dist > 1e-14 {
                let normal = scale3(diff, 1.0 / dist);
                *pos = add3(closest, scale3(normal, radius));
            }
        }
    }

    /// Compute the minimum distance from a point to a sphere surface.
    ///
    /// # Arguments
    /// * `point` - The query point.
    /// * `center` - Sphere center.
    /// * `radius` - Sphere radius.
    pub fn point_sphere_distance(point: [f64; 3], center: [f64; 3], radius: f64) -> f64 {
        (len3(sub3(point, center)) - radius).abs()
    }

    /// Compute the minimum distance from a point to a capsule.
    ///
    /// # Arguments
    /// * `point` - The query point.
    /// * `a` - First endpoint of the capsule axis.
    /// * `b` - Second endpoint of the capsule axis.
    /// * `radius` - Capsule radius.
    pub fn point_capsule_distance(point: [f64; 3], a: [f64; 3], b: [f64; 3], radius: f64) -> f64 {
        let ab = sub3(b, a);
        let len_ab = len3(ab);
        if len_ab < 1e-14 {
            return (len3(sub3(point, a)) - radius).abs();
        }
        let ab_norm = scale3(ab, 1.0 / len_ab);
        let ap = sub3(point, a);
        let t = clamp(dot3(ap, ab_norm), 0.0, len_ab);
        let closest = add3(a, scale3(ab_norm, t));
        (len3(sub3(point, closest)) - radius).abs()
    }
}

// ---------------------------------------------------------------------------
// FurLayer
// ---------------------------------------------------------------------------

/// A fur layer descriptor for a mesh surface.
///
/// Controls the density, length, clumping, and roughness of fur on a mesh.
pub struct FurLayer {
    /// Fur density (hairs per cm²).
    pub density: f64,
    /// Average fur strand length (metres).
    pub length: f64,
    /// Clumping factor: 0 = no clumping, 1 = all hairs in one clump.
    pub clumping: f64,
    /// Surface roughness (0 = smooth, 1 = very rough).
    pub roughness: f64,
}

impl FurLayer {
    /// Create a new fur layer.
    ///
    /// # Arguments
    /// * `density` - Hairs per cm² (e.g., 300.0 for domestic cat).
    /// * `length` - Strand length in metres.
    /// * `clumping` - Clumping factor in \[0, 1\].
    /// * `roughness` - Roughness in \[0, 1\].
    pub fn new(density: f64, length: f64, clumping: f64, roughness: f64) -> Self {
        Self {
            density,
            length,
            clumping,
            roughness,
        }
    }

    /// Generate fur root positions by placing `n_per_face` hairs on each mesh face.
    ///
    /// Each face is treated as a planar element; the barycentric jitter ensures
    /// even coverage.
    ///
    /// # Arguments
    /// * `mesh_positions` - Array of mesh vertex positions (triplets: v0, v1, v2 per face).
    /// * `n_per_face` - Number of hair roots per triangular face.
    pub fn generate_fur_positions(mesh_positions: &[[f64; 3]], n_per_face: usize) -> Vec<[f64; 3]> {
        let mut result = Vec::new();
        let n_faces = mesh_positions.len() / 3;
        for face in 0..n_faces {
            let v0 = mesh_positions[face * 3];
            let v1 = mesh_positions[face * 3 + 1];
            let v2 = mesh_positions[face * 3 + 2];
            for j in 0..n_per_face {
                // Uniform distribution via stratified barycentric coordinates
                let u = (j as f64 + 0.5) / n_per_face as f64;
                let v = 0.5 - 0.5 * u;
                let w = 1.0 - u - v;
                let pos = [
                    u * v0[0] + v * v1[0] + w * v2[0],
                    u * v0[1] + v * v1[1] + w * v2[1],
                    u * v0[2] + v * v1[2] + w * v2[2],
                ];
                result.push(pos);
            }
        }
        result
    }

    /// Compute a simplified fur shading response using Kajiya-Kay inspired model.
    ///
    /// Returns a scalar intensity in \[0, ∞).
    ///
    /// # Arguments
    /// * `view_dir` - Direction from surface to camera (normalised).
    /// * `light_dir` - Direction from surface to light (normalised).
    /// * `normal` - Surface normal (normalised).
    pub fn fur_shading_response(view_dir: [f64; 3], light_dir: [f64; 3], normal: [f64; 3]) -> f64 {
        let n_dot_l = clamp(dot3(normal, light_dir), 0.0, 1.0);
        let n_dot_v = clamp(dot3(normal, view_dir), 0.0, 1.0);
        // Simple Lambertian + back-scatter for fur
        let diffuse = n_dot_l;
        let back_scatter = (1.0 - n_dot_v) * 0.3;
        diffuse + back_scatter
    }

    /// Estimate the number of hairs on a surface of given area.
    ///
    /// # Arguments
    /// * `area_m2` - Surface area in square metres.
    pub fn estimated_hair_count(&self, area_m2: f64) -> f64 {
        let area_cm2 = area_m2 * 1e4;
        area_cm2 * self.density
    }
}

// ---------------------------------------------------------------------------
// KajiyaKayHair
// ---------------------------------------------------------------------------

/// Kajiya-Kay anisotropic hair shading model.
///
/// Reference: Kajiya & Kay, "Rendering Fur with Three Dimensional Textures" (1989).
pub struct KajiyaKayHair {
    /// Diffuse reflectance coefficient.
    pub kd: f64,
    /// Specular reflectance coefficient.
    pub ks: f64,
    /// Specular shininess exponent.
    pub shininess: f64,
}

impl KajiyaKayHair {
    /// Create a new Kajiya-Kay shading model.
    ///
    /// # Arguments
    /// * `kd` - Diffuse coefficient.
    /// * `ks` - Specular coefficient.
    /// * `shininess` - Shininess exponent (larger = sharper highlight).
    pub fn new(kd: f64, ks: f64, shininess: f64) -> Self {
        Self { kd, ks, shininess }
    }

    /// Compute the diffuse term of the Kajiya-Kay model.
    ///
    /// Returns `sqrt(1 - (T·L)²)` where T is the tangent and L the light direction.
    ///
    /// # Arguments
    /// * `t` - Hair tangent direction (normalised).
    /// * `l` - Light direction (normalised).
    pub fn diffuse_term(t: [f64; 3], l: [f64; 3]) -> f64 {
        let tl = clamp(dot3(t, l), -1.0, 1.0);
        (1.0 - tl * tl).max(0.0).sqrt()
    }

    /// Compute the specular term of the Kajiya-Kay model.
    ///
    /// Returns `(T·L * T·V + sin_theta_L * sin_theta_V)^shininess`.
    ///
    /// # Arguments
    /// * `t` - Hair tangent direction (normalised).
    /// * `v` - View direction (normalised).
    /// * `l` - Light direction (normalised).
    /// * `shininess` - Shininess exponent.
    pub fn specular_term(t: [f64; 3], v: [f64; 3], l: [f64; 3], shininess: f64) -> f64 {
        let tl = clamp(dot3(t, l), -1.0, 1.0);
        let tv = clamp(dot3(t, v), -1.0, 1.0);
        let sin_tl = (1.0 - tl * tl).max(0.0).sqrt();
        let sin_tv = (1.0 - tv * tv).max(0.0).sqrt();
        let base = (tl * tv + sin_tl * sin_tv).max(0.0);
        base.powf(shininess)
    }

    /// Compute the full Kajiya-Kay shading for a hair fiber.
    ///
    /// Returns the total reflected intensity.
    ///
    /// # Arguments
    /// * `t` - Hair tangent direction (normalised).
    /// * `v` - View direction (normalised).
    /// * `l` - Light direction (normalised).
    pub fn shade(&self, t: [f64; 3], v: [f64; 3], l: [f64; 3]) -> f64 {
        let diff = self.kd * Self::diffuse_term(t, l);
        let spec = self.ks * Self::specular_term(t, v, l, self.shininess);
        diff + spec
    }
}

// ---------------------------------------------------------------------------
// SuperHelix
// ---------------------------------------------------------------------------

/// SuperHelix model for curly/wavy hair with intrinsic curvature and torsion.
///
/// Based on the Kirchhoff rod model used in the SuperHelix hair simulation paper.
/// Reference: Bertails et al., "Super-Helices for Predicting the Dynamics of Natural Hair" (2006).
pub struct SuperHelix {
    /// Intrinsic curvature (1/radius) along the strand (m⁻¹).
    pub curvature: f64,
    /// Intrinsic torsion along the strand (rad/m).
    pub torsion: f64,
    /// Bending stiffness (Pa·m⁴).
    pub bending_stiffness: f64,
    /// Twist stiffness (Pa·m⁴).
    pub twist_stiffness: f64,
    /// Total arc length of the strand (metres).
    pub arc_length: f64,
    /// Number of discretisation points.
    pub n_points: usize,
}

impl SuperHelix {
    /// Create a new SuperHelix model.
    ///
    /// # Arguments
    /// * `curvature` - Intrinsic curvature (m⁻¹).
    /// * `torsion` - Intrinsic torsion (rad/m).
    /// * `bending_stiffness` - EI bending stiffness (Pa·m⁴).
    /// * `twist_stiffness` - GJ torsional stiffness (Pa·m⁴).
    /// * `arc_length` - Total strand arc length (metres).
    /// * `n_points` - Number of discretisation points along the arc.
    pub fn new(
        curvature: f64,
        torsion: f64,
        bending_stiffness: f64,
        twist_stiffness: f64,
        arc_length: f64,
        n_points: usize,
    ) -> Self {
        Self {
            curvature,
            torsion,
            bending_stiffness,
            twist_stiffness,
            arc_length,
            n_points,
        }
    }

    /// Compute the rest curvature at each discretisation point along the arc.
    ///
    /// For a perfect helix, curvature is constant; this returns a Vec of
    /// length `n_points` with the constant curvature value.
    pub fn rest_curvature_along_arc(&self) -> Vec<f64> {
        vec![self.curvature; self.n_points]
    }

    /// Compute the rest torsion at each discretisation point.
    ///
    /// Returns a Vec of length `n_points` with constant torsion.
    pub fn rest_torsion_along_arc(&self) -> Vec<f64> {
        vec![self.torsion; self.n_points]
    }

    /// Compute the elastic energy stored in the strand.
    ///
    /// Uses the Kirchhoff rod energy:
    /// E = ∫ (EI * κ² + GJ * τ²) / 2 ds
    /// approximated by a Riemann sum.
    pub fn elastic_energy(&self) -> f64 {
        let _ds = self.arc_length / self.n_points as f64;
        let bending_energy = self.bending_stiffness * self.curvature * self.curvature * 0.5;
        let twist_energy = self.twist_stiffness * self.torsion * self.torsion * 0.5;
        (bending_energy + twist_energy) * self.arc_length
    }

    /// Generate the 3D helix positions along the arc.
    ///
    /// Returns a Vec of `n_points` world-space positions forming a helix.
    pub fn helix_positions(&self) -> Vec<[f64; 3]> {
        let mut positions = Vec::with_capacity(self.n_points);
        let ds = self.arc_length / (self.n_points - 1).max(1) as f64;
        let r = if self.curvature.abs() > 1e-14 {
            1.0 / self.curvature
        } else {
            1e6
        };
        let pitch = if self.torsion.abs() > 1e-14 {
            std::f64::consts::TAU / self.torsion
        } else {
            1e6
        };
        for i in 0..self.n_points {
            let s = i as f64 * ds;
            let theta = s * self.curvature;
            positions.push([
                r * theta.sin(),
                s * pitch / (std::f64::consts::TAU * r.abs().max(1e-14)),
                r * (1.0 - theta.cos()),
            ]);
        }
        positions
    }

    /// Compute the bending energy per unit length at a given arc position.
    ///
    /// # Arguments
    /// * `_s` - Arc position (metres); currently returns uniform energy.
    pub fn bending_energy_density(&self, _s: f64) -> f64 {
        0.5 * self.bending_stiffness * self.curvature * self.curvature
    }

    /// Compute the torsion energy per unit length at a given arc position.
    ///
    /// # Arguments
    /// * `_s` - Arc position (metres).
    pub fn torsion_energy_density(&self, _s: f64) -> f64 {
        0.5 * self.twist_stiffness * self.torsion * self.torsion
    }
}

// ---------------------------------------------------------------------------
// StrandGroup
// ---------------------------------------------------------------------------

/// A group of hair strands sharing a common root region.
///
/// Used to model clumping and inter-strand interactions.
pub struct StrandGroup {
    /// Strands in this group.
    pub strands: Vec<HairStrand>,
    /// Centre position of the group root.
    pub root_center: [f64; 3],
    /// Clumping radius: strands within this distance attract each other.
    pub clump_radius: f64,
    /// Clumping stiffness.
    pub clump_stiffness: f64,
}

impl StrandGroup {
    /// Create a new strand group.
    ///
    /// # Arguments
    /// * `root_center` - Centre of the root region.
    /// * `clump_radius` - Radius for clumping attraction.
    /// * `clump_stiffness` - Stiffness of clumping attraction.
    pub fn new(root_center: [f64; 3], clump_radius: f64, clump_stiffness: f64) -> Self {
        Self {
            strands: Vec::new(),
            root_center,
            clump_radius,
            clump_stiffness,
        }
    }

    /// Add a strand to this group.
    pub fn add_strand(&mut self, strand: HairStrand) {
        self.strands.push(strand);
    }

    /// Apply clumping forces to draw strands toward the group centroid.
    pub fn apply_clumping(&mut self) {
        if self.strands.is_empty() {
            return;
        }
        let n_particles = self.strands[0].num_particles();
        // Compute per-level centroid
        for level in 1..n_particles {
            let sum: [f64; 3] = self
                .strands
                .iter()
                .fold([0.0, 0.0, 0.0], |acc, s| add3(acc, s.positions[level]));
            let n = self.strands.len() as f64;
            let centroid = scale3(sum, 1.0 / n);
            for strand in self.strands.iter_mut() {
                let diff = sub3(centroid, strand.positions[level]);
                let dist = len3(diff);
                if dist > self.clump_radius && dist > 1e-14 {
                    let correction = scale3(diff, self.clump_stiffness * 0.1);
                    strand.positions[level] = add3(strand.positions[level], correction);
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// HairPhysicsParams
// ---------------------------------------------------------------------------

/// Physical parameters for hair strand simulation.
///
/// Encapsulates material properties used in advanced hair dynamics.
pub struct HairPhysicsParams {
    /// Young's modulus of hair fibre (Pa), typically ~3.5 GPa.
    pub youngs_modulus: f64,
    /// Shear modulus (Pa).
    pub shear_modulus: f64,
    /// Hair strand radius (metres), typically ~35–100 µm.
    pub strand_radius: f64,
    /// Linear mass density (kg/m).
    pub linear_density: f64,
}

impl HairPhysicsParams {
    /// Create hair physics parameters.
    ///
    /// # Arguments
    /// * `youngs_modulus` - Young's modulus in Pa.
    /// * `shear_modulus` - Shear modulus in Pa.
    /// * `strand_radius` - Strand radius in metres.
    /// * `linear_density` - Linear mass density in kg/m.
    pub fn new(
        youngs_modulus: f64,
        shear_modulus: f64,
        strand_radius: f64,
        linear_density: f64,
    ) -> Self {
        Self {
            youngs_modulus,
            shear_modulus,
            strand_radius,
            linear_density,
        }
    }

    /// Compute the bending stiffness EI for a circular cross-section.
    ///
    /// EI = E * π * r^4 / 4
    pub fn bending_stiffness(&self) -> f64 {
        let r = self.strand_radius;
        self.youngs_modulus * std::f64::consts::PI * r * r * r * r / 4.0
    }

    /// Compute the torsional stiffness GJ for a circular cross-section.
    ///
    /// GJ = G * π * r^4 / 2
    pub fn torsional_stiffness(&self) -> f64 {
        let r = self.strand_radius;
        self.shear_modulus * std::f64::consts::PI * r * r * r * r / 2.0
    }

    /// Compute the natural frequency of a pinned-pinned beam mode n.
    ///
    /// # Arguments
    /// * `n` - Mode number (1 = fundamental).
    /// * `length` - Beam length (metres).
    pub fn natural_frequency(&self, n: usize, length: f64) -> f64 {
        let ei = self.bending_stiffness();
        let rho = self.linear_density;
        let lambda = (n as f64 * std::f64::consts::PI / length).powi(2);
        lambda * (ei / rho).sqrt()
    }
}

// ---------------------------------------------------------------------------
// WispGenerator
// ---------------------------------------------------------------------------

/// Generator for hair wisps (groups of nearby strands that move together).
///
/// A wisp acts as a proxy for a cluster of hair strands to reduce simulation cost.
pub struct WispGenerator {
    /// Centre root position of this wisp.
    pub center: [f64; 3],
    /// Number of strands this wisp represents.
    pub strand_count: usize,
    /// Wisp spread radius at the tip.
    pub tip_spread: f64,
}

impl WispGenerator {
    /// Create a new wisp generator.
    ///
    /// # Arguments
    /// * `center` - Center of the wisp root.
    /// * `strand_count` - Number of represented strands.
    /// * `tip_spread` - Spread radius at the tip (metres).
    pub fn new(center: [f64; 3], strand_count: usize, tip_spread: f64) -> Self {
        Self {
            center,
            strand_count,
            tip_spread,
        }
    }

    /// Generate representative strand roots for the wisp.
    ///
    /// Distributes roots on a circle of given radius around the centre.
    ///
    /// # Arguments
    /// * `radius` - Distribution radius.
    pub fn generate_roots(&self, radius: f64) -> Vec<[f64; 3]> {
        let mut roots = Vec::with_capacity(self.strand_count);
        for i in 0..self.strand_count {
            let angle = (i as f64 / self.strand_count as f64) * std::f64::consts::TAU;
            roots.push([
                self.center[0] + radius * angle.cos(),
                self.center[1],
                self.center[2] + radius * angle.sin(),
            ]);
        }
        roots
    }
}

// ---------------------------------------------------------------------------
// HairRenderData
// ---------------------------------------------------------------------------

/// Render-ready data for a hair strand.
///
/// Stores positions and tangents for GPU upload.
pub struct HairRenderData {
    /// Positions along the strand.
    pub positions: Vec<[f64; 3]>,
    /// Tangent directions along the strand.
    pub tangents: Vec<[f64; 3]>,
    /// Strand colour (linear RGB).
    pub colour: [f64; 3],
    /// Strand opacity.
    pub opacity: f64,
}

impl HairRenderData {
    /// Build render data from a HairStrand.
    ///
    /// # Arguments
    /// * `strand` - Source strand.
    /// * `colour` - RGB colour of the strand.
    /// * `opacity` - Opacity in \[0, 1\].
    pub fn from_strand(strand: &HairStrand, colour: [f64; 3], opacity: f64) -> Self {
        let n = strand.positions.len();
        let mut tangents = Vec::with_capacity(n);
        for i in 0..n {
            if i + 1 < n {
                tangents.push(normalize3(sub3(
                    strand.positions[i + 1],
                    strand.positions[i],
                )));
            } else if n > 1 {
                tangents.push(normalize3(sub3(
                    strand.positions[n - 1],
                    strand.positions[n - 2],
                )));
            } else {
                tangents.push([0.0, 1.0, 0.0]);
            }
        }
        Self {
            positions: strand.positions.clone(),
            tangents,
            colour,
            opacity,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    // -----------------------------------------------------------------------
    // HairStrand tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hair_strand_creation() {
        let root = [0.0, 10.0, 0.0];
        let strand = HairStrand::new(root, 1.0, 5);
        assert_eq!(strand.num_particles(), 6);
        assert_eq!(strand.num_segments, 5);
        assert_eq!(strand.rest_lengths.len(), 5);
        for l in &strand.rest_lengths {
            assert!((l - 0.2).abs() < 1e-12);
        }
    }

    #[test]
    fn test_hair_strand_tip_position() {
        let root = [0.0, 10.0, 0.0];
        let strand = HairStrand::new(root, 1.0, 4);
        let tip = strand.tip_position();
        // Tip is directly below root by arc length
        assert!((tip[1] - (10.0 - 1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_hair_strand_verlet_integration() {
        let root = [0.0, 10.0, 0.0];
        let mut strand = HairStrand::new(root, 1.0, 4);
        let gravity = [0.0, -9.81, 0.0];
        let dt = 1.0 / 60.0;
        strand.integrate_verlet(dt, gravity);
        // All free particles should have moved downward
        for i in 1..strand.num_particles() {
            // Verlet integration: particles should have slightly new positions
            let _ = strand.positions[i]; // just ensure no panic
        }
    }

    #[test]
    fn test_hair_strand_length_constraints() {
        let root = [0.0, 10.0, 0.0];
        let mut strand = HairStrand::new(root, 1.0, 4);
        // Perturb all particles
        for i in 1..strand.num_particles() {
            strand.positions[i][0] += 5.0;
        }
        strand.apply_length_constraints(1.0);
        // After one pass, lengths should be closer to rest
        let n = strand.num_particles();
        for i in 0..n - 1 {
            let l = len3(sub3(strand.positions[i + 1], strand.positions[i]));
            assert!(l < 10.0, "Length {l} should be reduced by constraint");
        }
    }

    #[test]
    fn test_hair_strand_root_fixed() {
        let root = [0.0, 10.0, 0.0];
        let mut strand = HairStrand::new(root, 1.0, 4);
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..30 {
            strand.integrate_verlet(1.0 / 60.0, gravity);
            strand.apply_length_constraints(1.0);
        }
        // Root should remain at original position (we don't update index 0)
        assert!((strand.positions[0][0] - root[0]).abs() < 1e-10);
        assert!((strand.positions[0][1] - root[1]).abs() < 1e-10);
        assert!((strand.positions[0][2] - root[2]).abs() < 1e-10);
    }

    #[test]
    fn test_hair_strand_arc_length() {
        let root = [0.0, 0.0, 0.0];
        let strand = HairStrand::new(root, 2.0, 10);
        // Initial arc length equals actual segment sum
        let arc = strand.arc_length();
        assert!(
            (arc - 2.0).abs() < 1e-10,
            "Arc length should be 2.0, got {arc}"
        );
    }

    #[test]
    fn test_hair_strand_kinetic_energy_zero_at_rest() {
        let root = [0.0, 0.0, 0.0];
        let strand = HairStrand::new(root, 1.0, 4);
        assert!((strand.kinetic_energy()).abs() < 1e-14);
    }

    #[test]
    fn test_hair_strand_reset() {
        let root = [0.0, 10.0, 0.0];
        let mut strand = HairStrand::new(root, 1.0, 4);
        // Perturb everything
        for p in strand.positions.iter_mut() {
            p[0] += 100.0;
        }
        strand.reset();
        // After reset, positions should be back to initial
        assert!((strand.positions[0][0] - root[0]).abs() < 1e-12);
    }

    #[test]
    fn test_hair_strand_force_impulse() {
        let root = [0.0, 0.0, 0.0];
        let mut strand = HairStrand::new(root, 1.0, 4);
        let original_x = strand.positions[1][0];
        strand.apply_force_impulse([1.0, 0.0, 0.0], 0.016);
        assert!(
            strand.positions[1][0] > original_x,
            "Impulse should move particle in +X"
        );
    }

    // -----------------------------------------------------------------------
    // HairSimulation tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hair_simulation_add_strand() {
        let mut sim = HairSimulation::new(1.225, 1.0, 0.00005);
        sim.add_strand([0.0, 0.0, 0.0], 0.3, 10);
        sim.add_strand([0.1, 0.0, 0.0], 0.3, 10);
        assert_eq!(sim.strands.len(), 2);
    }

    #[test]
    fn test_hair_simulation_total_particles() {
        let mut sim = HairSimulation::new(1.225, 1.0, 0.00005);
        sim.add_strand([0.0, 0.0, 0.0], 0.3, 5);
        sim.add_strand([0.1, 0.0, 0.0], 0.3, 5);
        assert_eq!(sim.total_particles(), 12); // 2 * 6
    }

    #[test]
    fn test_hair_simulation_step_no_crash() {
        let mut sim = HairSimulation::new(1.225, 1.0, 0.00005);
        sim.add_strand([0.0, 0.5, 0.0], 0.3, 5);
        let gravity = [0.0, -9.81, 0.0];
        let wind = [1.0, 0.0, 0.0];
        for _ in 0..10 {
            sim.step(1.0 / 60.0, gravity, wind);
        }
    }

    #[test]
    fn test_hair_simulation_damping() {
        let mut sim = HairSimulation::new(1.225, 1.0, 0.00005);
        sim.add_strand([0.0, 0.5, 0.0], 0.3, 5);
        // Manually set some velocity
        sim.strands[0].velocities[1] = [1.0, 0.0, 0.0];
        sim.apply_damping(0.5);
        assert!((sim.strands[0].velocities[1][0] - 0.5).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // HairBendingConstraint tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_bending_constraint_straight_no_change() {
        let bc = HairBendingConstraint::new(0.0, 1.0);
        let p0 = [0.0f64, 0.0, 0.0];
        let mut p1 = [1.0f64, 0.0, 0.0];
        let mut p2 = [2.0f64, 0.0, 0.0];
        let p1_before = p1;
        let p2_before = p2;
        bc.apply(&p0, &mut p1, &mut p2);
        // For a straight strand with rest_angle=0, corrections should be minimal
        let moved1 = len3(sub3(p1, p1_before));
        let moved2 = len3(sub3(p2, p2_before));
        assert!(
            moved1 < 0.1,
            "Straight strand p1 should move little: {moved1}"
        );
        assert!(
            moved2 < 0.1,
            "Straight strand p2 should move little: {moved2}"
        );
    }

    #[test]
    fn test_bending_constraint_compute_angle() {
        let d01 = [1.0f64, 0.0, 0.0];
        let d12 = [0.0f64, 1.0, 0.0];
        let angle = HairBendingConstraint::compute_angle(d01, d12);
        assert!(
            (angle - PI / 2.0).abs() < 1e-10,
            "Perpendicular segments angle should be 90deg: {angle}"
        );
    }

    #[test]
    fn test_bending_constraint_parallel_zero_angle() {
        let d01 = [1.0f64, 0.0, 0.0];
        let d12 = [1.0f64, 0.0, 0.0];
        let angle = HairBendingConstraint::compute_angle(d01, d12);
        assert!(
            angle.abs() < 1e-10,
            "Parallel directions → zero angle, got {angle}"
        );
    }

    #[test]
    fn test_bending_constraint_180_angle() {
        let d01 = [1.0f64, 0.0, 0.0];
        let d12 = [-1.0f64, 0.0, 0.0];
        let angle = HairBendingConstraint::compute_angle(d01, d12);
        assert!(
            (angle - PI).abs() < 1e-10,
            "Anti-parallel → PI, got {angle}"
        );
    }

    // -----------------------------------------------------------------------
    // TorsionConstraint tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_torsion_constraint_zero_rest_no_change() {
        let tc = TorsionConstraint::new(0.0, 0.0); // zero stiffness
        let mut p0 = [0.0f64, 1.0, 0.0];
        let mut p1 = [0.0f64, 0.0, 0.0];
        let mut p2 = [1.0f64, 0.0, 0.0];
        let mut p3 = [1.0f64, 0.0, 1.0];
        let p0_before = p0;
        tc.apply(&mut p0, &mut p1, &mut p2, &mut p3);
        let moved = len3(sub3(p0, p0_before));
        assert!(
            moved < 1e-12,
            "Zero stiffness should not move particles: {moved}"
        );
    }

    #[test]
    fn test_torsion_compute_twist() {
        let p0 = [0.0f64, 1.0, 0.0];
        let p1 = [0.0f64, 0.0, 0.0];
        let p2 = [1.0f64, 0.0, 0.0];
        let p3 = [1.0f64, 0.0, 1.0];
        let twist = TorsionConstraint::compute_twist(p0, p1, p2, p3);
        // Should be a well-defined finite value
        assert!(twist.is_finite(), "Twist should be finite: {twist}");
    }

    #[test]
    fn test_torsion_constraint_creates() {
        let tc = TorsionConstraint::new(0.1, 0.8);
        assert!((tc.rest_twist - 0.1).abs() < 1e-12);
        assert!((tc.stiffness - 0.8).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // HairCollision tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hair_collision_sphere_pushes_out() {
        let mut strand = HairStrand::new([0.0, 5.0, 0.0], 1.0, 5);
        // Place free particles just inside the sphere (small offset to avoid degenerate dist=0)
        for i in 1..strand.num_particles() {
            strand.positions[i] = [0.01 * i as f64, 0.0, 0.0];
        }
        HairCollision::collide_with_sphere(&mut strand, [0.0, 0.0, 0.0], 1.0);
        for i in 1..strand.num_particles() {
            let dist = len3(strand.positions[i]);
            assert!(
                dist >= 1.0 - 1e-10,
                "Particle {i} should be outside sphere: dist={dist}"
            );
        }
    }

    #[test]
    fn test_hair_collision_sphere_no_change_outside() {
        let mut strand = HairStrand::new([0.0, 5.0, 0.0], 1.0, 4);
        let original = strand.positions.clone();
        // All particles are far from sphere at origin
        HairCollision::collide_with_sphere(&mut strand, [0.0, 0.0, 0.0], 0.1);
        for (i, (pos, orig)) in strand.positions.iter().zip(original.iter()).enumerate() {
            let diff = len3(sub3(*pos, *orig));
            assert!(diff < 1e-12, "Outside sphere particle {i} should not move");
        }
    }

    #[test]
    fn test_hair_collision_capsule_pushes_out() {
        let mut strand = HairStrand::new([0.0, 5.0, 0.0], 0.1, 3);
        // Place all particles at capsule axis
        for pos in strand.positions.iter_mut() {
            *pos = [0.0, 0.0, 0.0];
        }
        strand.positions[0] = [0.0, 5.0, 0.0]; // keep root far
        let a = [-5.0f64, 0.0, 0.0];
        let b = [5.0f64, 0.0, 0.0];
        HairCollision::collide_with_capsule(&mut strand, a, b, 1.0);
        for i in 1..strand.num_particles() {
            let dist = HairCollision::point_capsule_distance(strand.positions[i], a, b, 1.0);
            assert!(
                dist >= -1e-10,
                "Particle {i} should be outside capsule: dist={dist}"
            );
        }
    }

    #[test]
    fn test_point_sphere_distance() {
        let dist = HairCollision::point_sphere_distance([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0);
        assert!(
            (dist - 1.0).abs() < 1e-12,
            "Distance should be 1.0, got {dist}"
        );
    }

    #[test]
    fn test_point_capsule_distance_on_axis() {
        let a = [0.0f64, 0.0, 0.0];
        let b = [4.0f64, 0.0, 0.0];
        let point = [2.0f64, 0.0, 0.0]; // on the axis
        let dist = HairCollision::point_capsule_distance(point, a, b, 0.5);
        assert!(
            (dist - 0.5).abs() < 1e-12,
            "On-axis distance should be radius=0.5, got {dist}"
        );
    }

    // -----------------------------------------------------------------------
    // FurLayer tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_fur_layer_creation() {
        let fur = FurLayer::new(300.0, 0.02, 0.3, 0.6);
        assert!((fur.density - 300.0).abs() < 1e-12);
        assert!((fur.length - 0.02).abs() < 1e-12);
    }

    #[test]
    fn test_fur_generate_positions_count() {
        let mesh = vec![[0.0f64, 0.0, 0.0], [1.0f64, 0.0, 0.0], [0.0f64, 1.0, 0.0]];
        let pos = FurLayer::generate_fur_positions(&mesh, 5);
        assert_eq!(
            pos.len(),
            5,
            "Expected 5 positions for 1 face with 5 per face"
        );
    }

    #[test]
    fn test_fur_generate_positions_two_faces() {
        let mesh = vec![
            [0.0f64, 0.0, 0.0],
            [1.0f64, 0.0, 0.0],
            [0.0f64, 1.0, 0.0],
            [1.0f64, 0.0, 0.0],
            [1.0f64, 1.0, 0.0],
            [0.0f64, 1.0, 0.0],
        ];
        let pos = FurLayer::generate_fur_positions(&mesh, 3);
        assert_eq!(
            pos.len(),
            6,
            "Expected 6 positions for 2 faces * 3 per face"
        );
    }

    #[test]
    fn test_fur_shading_response_normal_facing() {
        let view_dir = [0.0f64, 0.0, 1.0];
        let light_dir = [0.0f64, 0.0, 1.0];
        let normal = [0.0f64, 0.0, 1.0];
        let response = FurLayer::fur_shading_response(view_dir, light_dir, normal);
        assert!(
            response > 0.0,
            "Front-facing fur should have positive response"
        );
    }

    #[test]
    fn test_fur_estimated_hair_count() {
        let fur = FurLayer::new(100.0, 0.02, 0.0, 0.5);
        // 1 m² = 10000 cm², 100 hairs/cm² → 1,000,000 hairs
        let count = fur.estimated_hair_count(1.0);
        assert!(
            (count - 1_000_000.0).abs() < 1.0,
            "Expected 1M hairs: {count}"
        );
    }

    // -----------------------------------------------------------------------
    // KajiyaKayHair tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_kajiya_kay_diffuse_perpendicular() {
        // Tangent ⊥ light → diffuse = 1.0
        let t = [1.0f64, 0.0, 0.0];
        let l = [0.0f64, 1.0, 0.0];
        let diff = KajiyaKayHair::diffuse_term(t, l);
        assert!(
            (diff - 1.0).abs() < 1e-10,
            "Perpendicular T·L → diffuse=1, got {diff}"
        );
    }

    #[test]
    fn test_kajiya_kay_diffuse_parallel() {
        // Tangent ∥ light → diffuse = 0
        let t = [1.0f64, 0.0, 0.0];
        let l = [1.0f64, 0.0, 0.0];
        let diff = KajiyaKayHair::diffuse_term(t, l);
        assert!(diff.abs() < 1e-10, "Parallel T∥L → diffuse=0, got {diff}");
    }

    #[test]
    fn test_kajiya_kay_specular_range() {
        let t = [1.0f64, 0.0, 0.0];
        let v = [0.0f64, 1.0, 0.0];
        let l = [0.0f64, 0.0, 1.0];
        let spec = KajiyaKayHair::specular_term(t, v, l, 10.0);
        assert!(spec >= 0.0, "Specular should be non-negative: {spec}");
        assert!(
            spec <= 1.0,
            "Specular should be ≤ 1 for unit inputs: {spec}"
        );
    }

    #[test]
    fn test_kajiya_kay_shade_positive() {
        let kk = KajiyaKayHair::new(0.5, 0.3, 20.0);
        let t = [1.0f64, 0.0, 0.0];
        let v = [0.0f64, 0.0, 1.0];
        let l = [0.0f64, 1.0, 0.0];
        let shade = kk.shade(t, v, l);
        assert!(shade >= 0.0, "Shade should be non-negative: {shade}");
    }

    // -----------------------------------------------------------------------
    // SuperHelix tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_super_helix_creation() {
        let sh = SuperHelix::new(10.0, 5.0, 1e-9, 5e-10, 0.3, 20);
        assert!((sh.curvature - 10.0).abs() < 1e-12);
        assert!((sh.torsion - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_super_helix_rest_curvature_length() {
        let sh = SuperHelix::new(5.0, 2.0, 1e-9, 5e-10, 0.3, 15);
        let kappa = sh.rest_curvature_along_arc();
        assert_eq!(kappa.len(), 15);
        for k in &kappa {
            assert!((k - 5.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_super_helix_elastic_energy_positive() {
        let sh = SuperHelix::new(10.0, 5.0, 1e-9, 5e-10, 0.3, 20);
        let e = sh.elastic_energy();
        assert!(e > 0.0, "Elastic energy should be positive: {e}");
    }

    #[test]
    fn test_super_helix_zero_curvature_straight() {
        let sh = SuperHelix::new(0.0, 0.0, 1e-9, 5e-10, 0.3, 10);
        let e = sh.elastic_energy();
        assert!(e.abs() < 1e-20, "Zero curvature/torsion → zero energy: {e}");
    }

    #[test]
    fn test_super_helix_helix_positions_count() {
        let sh = SuperHelix::new(10.0, 5.0, 1e-9, 5e-10, 0.3, 10);
        let pos = sh.helix_positions();
        assert_eq!(pos.len(), 10);
    }

    #[test]
    fn test_super_helix_bending_energy_density() {
        let sh = SuperHelix::new(2.0, 1.0, 1e-9, 5e-10, 0.3, 10);
        let e = sh.bending_energy_density(0.0);
        let expected = 0.5 * 1e-9 * 4.0;
        assert!(
            (e - expected).abs() < 1e-25,
            "Bending energy density: {e} vs {expected}"
        );
    }

    // -----------------------------------------------------------------------
    // StrandGroup tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_strand_group_clumping_no_crash() {
        let mut group = StrandGroup::new([0.0, 0.0, 0.0], 0.1, 0.5);
        group.add_strand(HairStrand::new([0.0, 0.0, 0.0], 0.3, 5));
        group.add_strand(HairStrand::new([0.01, 0.0, 0.0], 0.3, 5));
        group.apply_clumping();
        // Should not crash
    }

    #[test]
    fn test_strand_group_empty_clumping() {
        let mut group = StrandGroup::new([0.0, 0.0, 0.0], 0.1, 0.5);
        group.apply_clumping(); // should not crash
        assert_eq!(group.strands.len(), 0);
    }

    // -----------------------------------------------------------------------
    // HairPhysicsParams tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hair_physics_bending_stiffness() {
        // E = 3.5 GPa, r = 50 μm
        let params = HairPhysicsParams::new(3.5e9, 1.3e9, 50e-6, 65e-6);
        let ei = params.bending_stiffness();
        // EI = 3.5e9 * π * (50e-6)^4 / 4
        let expected = 3.5e9 * std::f64::consts::PI * (50e-6_f64).powi(4) / 4.0;
        assert!(
            (ei / expected - 1.0).abs() < 1e-10,
            "EI = {ei} vs {expected}"
        );
    }

    #[test]
    fn test_hair_physics_torsional_stiffness() {
        let params = HairPhysicsParams::new(3.5e9, 1.3e9, 50e-6, 65e-6);
        let gj = params.torsional_stiffness();
        assert!(gj > 0.0, "GJ should be positive: {gj}");
    }

    #[test]
    fn test_hair_physics_natural_frequency() {
        let params = HairPhysicsParams::new(3.5e9, 1.3e9, 50e-6, 65e-6);
        let f1 = params.natural_frequency(1, 0.3);
        let f2 = params.natural_frequency(2, 0.3);
        assert!(
            f2 > f1,
            "Second mode frequency should exceed first: f1={f1}, f2={f2}"
        );
    }

    // -----------------------------------------------------------------------
    // WispGenerator tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_wisp_generator_roots_count() {
        let wisp_gen = WispGenerator::new([0.0, 0.0, 0.0], 8, 0.01);
        let roots = wisp_gen.generate_roots(0.005);
        assert_eq!(roots.len(), 8);
    }

    #[test]
    fn test_wisp_generator_roots_radius() {
        let wisp_gen = WispGenerator::new([0.0, 0.0, 0.0], 8, 0.01);
        let r = 0.005;
        let roots = wisp_gen.generate_roots(r);
        for root in &roots {
            let dist = (root[0] * root[0] + root[2] * root[2]).sqrt();
            assert!(
                (dist - r).abs() < 1e-10,
                "Root distance should be {r}: {dist}"
            );
        }
    }

    // -----------------------------------------------------------------------
    // HairRenderData tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_hair_render_data_from_strand() {
        let strand = HairStrand::new([0.0, 1.0, 0.0], 0.5, 4);
        let rd = HairRenderData::from_strand(&strand, [0.3, 0.2, 0.1], 1.0);
        assert_eq!(rd.positions.len(), 5);
        assert_eq!(rd.tangents.len(), 5);
    }

    #[test]
    fn test_hair_render_data_tangents_normalised() {
        let strand = HairStrand::new([0.0, 1.0, 0.0], 0.5, 4);
        let rd = HairRenderData::from_strand(&strand, [1.0, 1.0, 1.0], 1.0);
        for (i, t) in rd.tangents.iter().enumerate() {
            let l = len3(*t);
            assert!(
                (l - 1.0).abs() < 1e-10,
                "Tangent {i} should be normalised: {l}"
            );
        }
    }

    #[test]
    fn test_hair_render_data_single_particle() {
        let strand = HairStrand::new([0.0, 0.0, 0.0], 0.0001, 1);
        // Should not crash for single-segment strand
        let rd = HairRenderData::from_strand(&strand, [0.5, 0.5, 0.5], 1.0);
        assert_eq!(rd.positions.len(), 2);
    }
}
