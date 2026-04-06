#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hair simulation using discrete elastic rods and Cosserat rod theory.
//!
//! Provides a full hair simulation pipeline:
//!
//! - [`DiscreteElasticRod`]: stretch, bend, and twist energy for a single strand.
//! - [`CosseratHairStrand`]: Cosserat rod model with material frames.
//! - [`HairHairInteraction`]: penalty-based collision between nearby strands.
//! - [`StylingForce`]: target-shape forces for artistic control.
//! - [`GravityDraper`]: initial drape solver under gravity.
//! - [`WindResponseModel`]: aerodynamic drag and flutter on hair strands.
//! - [`HairTangentFrame`]: rendering data (tangent, normal, binormal per vertex).
//! - [`GuideHairInterpolator`]: guide-hair + interpolated-strand scheme.
//! - [`HairLod`]: level-of-detail system for hair bundles.
//!
//! All vectors use `[f64; 3]` arrays (no nalgebra dependency).

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Small vector helpers
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
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-12 {
        [0.0; 3]
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
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    len3(sub3(a, b))
}

#[inline]
fn lerp3(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

#[inline]
fn neg3(v: [f64; 3]) -> [f64; 3] {
    [-v[0], -v[1], -v[2]]
}

#[inline]
fn zero3() -> [f64; 3] {
    [0.0; 3]
}

#[inline]
fn clamp_f64(v: f64, lo: f64, hi: f64) -> f64 {
    if v < lo {
        lo
    } else if v > hi {
        hi
    } else {
        v
    }
}

// ===========================================================================
// Discrete Elastic Rod (DER)
// ===========================================================================

/// A single hair strand modelled as a Discrete Elastic Rod.
///
/// Stores vertex positions, velocities, rest lengths, and curvatures.
/// Supports stretch, bend, and twist energy computation and force application.
#[derive(Debug, Clone)]
pub struct DiscreteElasticRod {
    /// World-space vertex positions (root at index 0).
    pub positions: Vec<[f64; 3]>,
    /// World-space vertex velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Inverse mass per vertex (0 = fixed/root).
    pub inv_mass: Vec<f64>,
    /// Rest edge lengths.
    pub rest_lengths: Vec<f64>,
    /// Rest curvature (kappa) per interior vertex.
    pub rest_curvatures: Vec<f64>,
    /// Stretching stiffness \[N\].
    pub stretch_stiffness: f64,
    /// Bending stiffness \[N m^2\].
    pub bend_stiffness: f64,
    /// Twist stiffness \[N m^2\].
    pub twist_stiffness: f64,
    /// Damping coefficient.
    pub damping: f64,
    /// Rest twist angles per edge.
    pub rest_twists: Vec<f64>,
}

impl DiscreteElasticRod {
    /// Create a straight rod from root to tip with `n_segments` segments.
    pub fn new_straight(
        root: [f64; 3],
        tip: [f64; 3],
        n_segments: usize,
        mass_per_vertex: f64,
        stretch_stiffness: f64,
        bend_stiffness: f64,
        twist_stiffness: f64,
        damping: f64,
    ) -> Self {
        let n_verts = n_segments + 1;
        let mut positions = Vec::with_capacity(n_verts);
        for i in 0..n_verts {
            let t = i as f64 / n_segments as f64;
            positions.push(lerp3(root, tip, t));
        }
        let mut inv_mass = vec![1.0 / mass_per_vertex; n_verts];
        inv_mass[0] = 0.0; // root is fixed

        let mut rest_lengths = Vec::with_capacity(n_segments);
        for i in 0..n_segments {
            rest_lengths.push(dist3(positions[i], positions[i + 1]));
        }

        let rest_curvatures = vec![0.0; n_verts.saturating_sub(2)];
        let rest_twists = vec![0.0; n_segments];
        let velocities = vec![zero3(); n_verts];

        Self {
            positions,
            velocities,
            inv_mass,
            rest_lengths,
            rest_curvatures,
            stretch_stiffness,
            bend_stiffness,
            twist_stiffness,
            damping,
            rest_twists,
        }
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.positions.len()
    }

    /// Number of edges (segments).
    pub fn num_edges(&self) -> usize {
        self.positions.len().saturating_sub(1)
    }

    /// Edge vector from vertex i to i+1.
    pub fn edge(&self, i: usize) -> [f64; 3] {
        sub3(self.positions[i + 1], self.positions[i])
    }

    /// Edge tangent (unit vector) for edge i.
    pub fn tangent(&self, i: usize) -> [f64; 3] {
        norm3(self.edge(i))
    }

    /// Current curvature at interior vertex i (1..n_verts-1).
    pub fn curvature(&self, i: usize) -> f64 {
        if i == 0 || i >= self.positions.len() - 1 {
            return 0.0;
        }
        let e0 = norm3(self.edge(i - 1));
        let e1 = norm3(self.edge(i));
        let cross = cross3(e0, e1);
        let sin_theta = len3(cross);
        let cos_theta = dot3(e0, e1);
        sin_theta.atan2(cos_theta)
    }

    /// Total stretch energy: sum of 0.5 * k_s * ((|e_i| - L_i) / L_i)^2 * L_i.
    pub fn stretch_energy(&self) -> f64 {
        let mut energy = 0.0;
        for i in 0..self.num_edges() {
            let l = len3(self.edge(i));
            let l0 = self.rest_lengths[i];
            if l0 > 1e-15 {
                let strain = (l - l0) / l0;
                energy += 0.5 * self.stretch_stiffness * strain * strain * l0;
            }
        }
        energy
    }

    /// Total bend energy: sum of 0.5 * k_b * (kappa_i - kappa_i0)^2 / l_bar.
    pub fn bend_energy(&self) -> f64 {
        let mut energy = 0.0;
        let n = self.positions.len();
        for i in 1..n.saturating_sub(1) {
            let kappa = self.curvature(i);
            let kappa0 = self.rest_curvatures[i - 1];
            let l_bar = (self.rest_lengths[i - 1] + self.rest_lengths[i]) / 2.0;
            if l_bar > 1e-15 {
                energy += 0.5 * self.bend_stiffness * (kappa - kappa0).powi(2) / l_bar;
            }
        }
        energy
    }

    /// Compute stretch forces and accumulate into the given force array.
    pub fn accumulate_stretch_forces(&self, forces: &mut [[f64; 3]]) {
        for i in 0..self.num_edges() {
            let e = self.edge(i);
            let l = len3(e);
            let l0 = self.rest_lengths[i];
            if l < 1e-15 || l0 < 1e-15 {
                continue;
            }
            let strain = (l - l0) / l0;
            let magnitude = self.stretch_stiffness * strain / l0;
            let dir = scale3(e, magnitude / l);
            forces[i] = add3(forces[i], dir);
            forces[i + 1] = sub3(forces[i + 1], dir);
        }
    }

    /// Compute bend forces and accumulate into the given force array.
    pub fn accumulate_bend_forces(&self, forces: &mut [[f64; 3]]) {
        let n = self.positions.len();
        for i in 1..n.saturating_sub(1) {
            let kappa = self.curvature(i);
            let kappa0 = if i - 1 < self.rest_curvatures.len() {
                self.rest_curvatures[i - 1]
            } else {
                0.0
            };
            let l_bar = (self.rest_lengths[i - 1] + self.rest_lengths[i]) / 2.0;
            if l_bar < 1e-15 {
                continue;
            }

            let dk = kappa - kappa0;
            let magnitude = self.bend_stiffness * dk / l_bar;

            // Approximate gradient: curvature binormal direction
            let e0 = self.edge(i - 1);
            let e1 = self.edge(i);
            let binormal = norm3(cross3(e0, e1));
            // Force pushes vertex i towards reducing curvature
            let f = scale3(binormal, -magnitude);
            forces[i] = add3(forces[i], f);
            // Distribute reaction to neighbors
            forces[i - 1] = add3(forces[i - 1], scale3(f, -0.5));
            forces[i + 1] = add3(forces[i + 1], scale3(f, -0.5));
        }
    }

    /// Apply damping to velocities.
    pub fn apply_damping(&mut self, dt: f64) {
        let factor = (-self.damping * dt).exp();
        for v in &mut self.velocities {
            *v = scale3(*v, factor);
        }
    }

    /// Integrate positions using symplectic Euler with given external force.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.num_vertices();
        let mut forces = vec![zero3(); n];

        // Gravity
        for i in 0..n {
            if self.inv_mass[i] > 0.0 {
                let m = 1.0 / self.inv_mass[i];
                forces[i] = scale3(gravity, m);
            }
        }

        self.accumulate_stretch_forces(&mut forces);
        self.accumulate_bend_forces(&mut forces);

        // Update velocities and positions
        for i in 0..n {
            if self.inv_mass[i] <= 0.0 {
                continue;
            }
            let acc = scale3(forces[i], self.inv_mass[i]);
            self.velocities[i] = add3(self.velocities[i], scale3(acc, dt));
        }

        self.apply_damping(dt);

        for i in 0..n {
            if self.inv_mass[i] <= 0.0 {
                continue;
            }
            self.positions[i] = add3(self.positions[i], scale3(self.velocities[i], dt));
        }
    }

    /// Total strand length (sum of edge lengths).
    pub fn total_length(&self) -> f64 {
        (0..self.num_edges()).map(|i| len3(self.edge(i))).sum()
    }

    /// Strand tip position.
    pub fn tip(&self) -> [f64; 3] {
        *self.positions.last().unwrap_or(&zero3())
    }

    /// Strand root position.
    pub fn root(&self) -> [f64; 3] {
        self.positions[0]
    }
}

// ===========================================================================
// Cosserat Rod for Hair
// ===========================================================================

/// Material frame for a Cosserat rod segment.
#[derive(Debug, Clone, Copy)]
pub struct MaterialFrame {
    /// Tangent (along the strand).
    pub tangent: [f64; 3],
    /// Normal (material direction 1).
    pub normal: [f64; 3],
    /// Binormal (material direction 2).
    pub binormal: [f64; 3],
}

impl MaterialFrame {
    /// Construct from tangent and a reference "up" direction.
    pub fn from_tangent_and_up(tangent: [f64; 3], up: [f64; 3]) -> Self {
        let t = norm3(tangent);
        let mut binormal = norm3(cross3(t, up));
        if len3(binormal) < 1e-10 {
            // tangent parallel to up; pick arbitrary perpendicular
            let arb = if t[0].abs() < 0.9 {
                [1.0, 0.0, 0.0]
            } else {
                [0.0, 1.0, 0.0]
            };
            binormal = norm3(cross3(t, arb));
        }
        let normal = norm3(cross3(binormal, t));
        Self {
            tangent: t,
            normal,
            binormal,
        }
    }

    /// Identity frame (tangent along +Y).
    pub fn identity() -> Self {
        Self {
            tangent: [0.0, -1.0, 0.0],
            normal: [1.0, 0.0, 0.0],
            binormal: [0.0, 0.0, 1.0],
        }
    }
}

/// Cosserat rod hair strand with material frames and twist.
///
/// Extends the DER model with explicit material frames for twist tracking,
/// suitable for curly/coiled hair.
#[derive(Debug, Clone)]
pub struct CosseratHairStrand {
    /// Underlying discrete elastic rod.
    pub rod: DiscreteElasticRod,
    /// Material frame per edge.
    pub frames: Vec<MaterialFrame>,
    /// Twist angle per edge \[radians\].
    pub twist_angles: Vec<f64>,
    /// Cross-section radius \[m\].
    pub radius: f64,
    /// Shear modulus \[Pa\].
    pub shear_modulus: f64,
}

impl CosseratHairStrand {
    /// Build a Cosserat strand from an existing DER.
    pub fn from_rod(rod: DiscreteElasticRod, radius: f64, shear_modulus: f64) -> Self {
        let ne = rod.num_edges();
        let mut frames = Vec::with_capacity(ne);
        let up = [0.0, 0.0, 1.0];
        for i in 0..ne {
            let t = rod.tangent(i);
            frames.push(MaterialFrame::from_tangent_and_up(t, up));
        }
        let twist_angles = vec![0.0; ne];
        Self {
            rod,
            frames,
            twist_angles,
            radius,
            shear_modulus,
        }
    }

    /// Create a straight Cosserat strand.
    pub fn new_straight(
        root: [f64; 3],
        tip: [f64; 3],
        n_segments: usize,
        mass_per_vertex: f64,
        stretch_stiffness: f64,
        bend_stiffness: f64,
        twist_stiffness: f64,
        damping: f64,
        radius: f64,
        shear_modulus: f64,
    ) -> Self {
        let rod = DiscreteElasticRod::new_straight(
            root,
            tip,
            n_segments,
            mass_per_vertex,
            stretch_stiffness,
            bend_stiffness,
            twist_stiffness,
            damping,
        );
        Self::from_rod(rod, radius, shear_modulus)
    }

    /// Update material frames by parallel-transporting along the strand.
    pub fn update_frames(&mut self) {
        if self.rod.num_edges() == 0 {
            return;
        }
        // First frame from tangent
        let t0 = self.rod.tangent(0);
        let up = if self.frames.is_empty() {
            [0.0, 0.0, 1.0]
        } else {
            self.frames[0].normal
        };
        self.frames[0] = MaterialFrame::from_tangent_and_up(t0, up);

        // Parallel transport subsequent frames
        for i in 1..self.rod.num_edges() {
            let t_prev = self.rod.tangent(i - 1);
            let t_curr = self.rod.tangent(i);
            let b = cross3(t_prev, t_curr);
            let b_len = len3(b);

            let prev_normal = self.frames[i - 1].normal;
            let new_normal = if b_len < 1e-12 {
                prev_normal
            } else {
                let axis = norm3(b);
                let cos_a = clamp_f64(dot3(t_prev, t_curr), -1.0, 1.0);
                let angle = cos_a.acos();
                rotate_around_axis(prev_normal, axis, angle)
            };
            self.frames[i] = MaterialFrame::from_tangent_and_up(t_curr, new_normal);
        }
    }

    /// Twist energy: sum of 0.5 * G * J * (theta_i - theta_i0)^2 / L_i.
    pub fn twist_energy(&self) -> f64 {
        let j = PI * self.radius.powi(4) / 2.0; // polar moment
        let mut energy = 0.0;
        for i in 0..self.rod.num_edges() {
            let l = self.rod.rest_lengths[i];
            if l < 1e-15 {
                continue;
            }
            let dtheta = self.twist_angles[i] - self.rod.rest_twists[i];
            energy += 0.5 * self.shear_modulus * j * dtheta * dtheta / l;
        }
        energy
    }

    /// Total energy (stretch + bend + twist).
    pub fn total_energy(&self) -> f64 {
        self.rod.stretch_energy() + self.rod.bend_energy() + self.twist_energy()
    }

    /// Step the strand forward in time.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        self.rod.step(dt, gravity);
        self.update_frames();
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.rod.num_vertices()
    }
}

/// Rotate vector `v` around unit axis `axis` by `angle` radians (Rodrigues).
fn rotate_around_axis(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let cos_a = angle.cos();
    let sin_a = angle.sin();
    let d = dot3(v, axis);
    let c = cross3(axis, v);
    [
        v[0] * cos_a + c[0] * sin_a + axis[0] * d * (1.0 - cos_a),
        v[1] * cos_a + c[1] * sin_a + axis[1] * d * (1.0 - cos_a),
        v[2] * cos_a + c[2] * sin_a + axis[2] * d * (1.0 - cos_a),
    ]
}

// ===========================================================================
// Hair-Hair Interaction (Penalty)
// ===========================================================================

/// Penalty-based hair-hair collision and repulsion.
///
/// Applies repulsive forces when two strand segments come within a
/// threshold distance.
#[derive(Debug, Clone)]
pub struct HairHairInteraction {
    /// Collision radius (sum of two strand radii).
    pub collision_radius: f64,
    /// Penalty stiffness \[N/m\].
    pub penalty_stiffness: f64,
    /// Friction coefficient for tangential damping.
    pub friction: f64,
}

impl HairHairInteraction {
    /// Create a new interaction model.
    pub fn new(collision_radius: f64, penalty_stiffness: f64, friction: f64) -> Self {
        Self {
            collision_radius,
            penalty_stiffness,
            friction,
        }
    }

    /// Compute repulsive force between two points if closer than threshold.
    ///
    /// Returns (force_on_a, force_on_b). Forces are equal and opposite.
    pub fn point_point_force(&self, a: [f64; 3], b: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let d = sub3(a, b);
        let dist = len3(d);
        if dist >= self.collision_radius || dist < 1e-15 {
            return (zero3(), zero3());
        }
        let overlap = self.collision_radius - dist;
        let magnitude = self.penalty_stiffness * overlap;
        let dir = scale3(d, magnitude / dist);
        (dir, neg3(dir))
    }

    /// Check if two segments (a0-a1, b0-b1) are within collision distance.
    pub fn segments_collide(&self, a0: [f64; 3], a1: [f64; 3], b0: [f64; 3], b1: [f64; 3]) -> bool {
        // Use midpoint approximation
        let ma = lerp3(a0, a1, 0.5);
        let mb = lerp3(b0, b1, 0.5);
        dist3(ma, mb) < self.collision_radius
    }

    /// Apply hair-hair forces between two strands.
    ///
    /// For each pair of vertices at the same index, applies penalty if close.
    pub fn apply_forces(
        &self,
        strand_a: &mut [[f64; 3]],
        strand_b: &mut [[f64; 3]],
        forces_a: &mut [[f64; 3]],
        forces_b: &mut [[f64; 3]],
    ) {
        let n = strand_a.len().min(strand_b.len());
        for i in 0..n {
            let (fa, fb) = self.point_point_force(strand_a[i], strand_b[i]);
            forces_a[i] = add3(forces_a[i], fa);
            forces_b[i] = add3(forces_b[i], fb);
        }
    }

    /// Number of collisions detected between two position arrays.
    pub fn count_collisions(&self, positions_a: &[[f64; 3]], positions_b: &[[f64; 3]]) -> usize {
        let n = positions_a.len().min(positions_b.len());
        let mut count = 0;
        for i in 0..n {
            if dist3(positions_a[i], positions_b[i]) < self.collision_radius {
                count += 1;
            }
        }
        count
    }
}

// ===========================================================================
// Styling Forces
// ===========================================================================

/// Target-shape styling force for artistic hair control.
///
/// Applies spring forces that pull each vertex towards a target "styled"
/// position, blended by a per-vertex weight.
#[derive(Debug, Clone)]
pub struct StylingForce {
    /// Target positions (styled shape).
    pub targets: Vec<[f64; 3]>,
    /// Per-vertex stiffness \[N/m\].
    pub stiffness: f64,
    /// Per-vertex blend weight \[0..1\].
    pub weights: Vec<f64>,
}

impl StylingForce {
    /// Create a styling force from target positions and uniform stiffness.
    pub fn new(targets: Vec<[f64; 3]>, stiffness: f64) -> Self {
        let n = targets.len();
        Self {
            targets,
            stiffness,
            weights: vec![1.0; n],
        }
    }

    /// Create with per-vertex weights.
    pub fn with_weights(targets: Vec<[f64; 3]>, stiffness: f64, weights: Vec<f64>) -> Self {
        Self {
            targets,
            stiffness,
            weights,
        }
    }

    /// Compute styling force for vertex i.
    pub fn force_at(&self, i: usize, current: [f64; 3]) -> [f64; 3] {
        if i >= self.targets.len() {
            return zero3();
        }
        let w = if i < self.weights.len() {
            self.weights[i]
        } else {
            1.0
        };
        let diff = sub3(self.targets[i], current);
        scale3(diff, self.stiffness * w)
    }

    /// Accumulate styling forces into a force array.
    pub fn accumulate(&self, positions: &[[f64; 3]], forces: &mut [[f64; 3]]) {
        let n = positions.len().min(self.targets.len()).min(forces.len());
        for i in 0..n {
            forces[i] = add3(forces[i], self.force_at(i, positions[i]));
        }
    }

    /// Total potential energy of the styling springs.
    pub fn energy(&self, positions: &[[f64; 3]]) -> f64 {
        let n = positions.len().min(self.targets.len());
        let mut e = 0.0;
        for i in 0..n {
            let w = if i < self.weights.len() {
                self.weights[i]
            } else {
                1.0
            };
            let d = dist3(positions[i], self.targets[i]);
            e += 0.5 * self.stiffness * w * d * d;
        }
        e
    }
}

// ===========================================================================
// Gravity Draping
// ===========================================================================

/// Solver for draping hair under gravity to find the rest pose.
///
/// Iteratively integrates with heavy damping until convergence.
#[derive(Debug, Clone)]
pub struct GravityDraper {
    /// Gravity vector \[m/s^2\].
    pub gravity: [f64; 3],
    /// Integration time step \[s\].
    pub dt: f64,
    /// Damping factor (0..1, high = more damping).
    pub damping: f64,
    /// Maximum iterations.
    pub max_iterations: usize,
    /// Convergence tolerance (max velocity magnitude).
    pub tolerance: f64,
}

impl Default for GravityDraper {
    fn default() -> Self {
        Self::new()
    }
}

impl GravityDraper {
    /// Create a default gravity draper.
    pub fn new() -> Self {
        Self {
            gravity: [0.0, -9.81, 0.0],
            dt: 0.01,
            damping: 0.95,
            max_iterations: 1000,
            tolerance: 1e-6,
        }
    }

    /// Create with custom gravity direction and parameters.
    pub fn with_params(
        gravity: [f64; 3],
        dt: f64,
        damping: f64,
        max_iterations: usize,
        tolerance: f64,
    ) -> Self {
        Self {
            gravity,
            dt,
            damping,
            max_iterations,
            tolerance,
        }
    }

    /// Drape a DER strand under gravity until convergence.
    /// Returns the number of iterations used.
    pub fn drape(&self, rod: &mut DiscreteElasticRod) -> usize {
        let original_damping = rod.damping;
        rod.damping = self.damping;

        let mut iter = 0;
        for i in 0..self.max_iterations {
            rod.step(self.dt, self.gravity);

            // Check convergence
            let max_vel = rod
                .velocities
                .iter()
                .map(|v| len3(*v))
                .fold(0.0_f64, f64::max);
            iter = i + 1;
            if max_vel < self.tolerance {
                break;
            }
        }

        rod.damping = original_damping;
        // Zero velocities after draping
        for v in &mut rod.velocities {
            *v = zero3();
        }
        iter
    }

    /// Drape a Cosserat strand.
    pub fn drape_cosserat(&self, strand: &mut CosseratHairStrand) -> usize {
        let iters = self.drape(&mut strand.rod);
        strand.update_frames();
        iters
    }

    /// Check if a strand is approximately in equilibrium.
    pub fn is_converged(&self, rod: &DiscreteElasticRod) -> bool {
        rod.velocities.iter().all(|v| len3(*v) < self.tolerance)
    }
}

// ===========================================================================
// Wind Response
// ===========================================================================

/// Aerodynamic model for wind interacting with hair strands.
///
/// Applies drag and lift forces based on relative wind velocity and
/// strand tangent direction.
#[derive(Debug, Clone)]
pub struct WindResponseModel {
    /// Air density \[kg/m^3\].
    pub air_density: f64,
    /// Drag coefficient (normal to strand).
    pub cd_normal: f64,
    /// Drag coefficient (tangential to strand).
    pub cd_tangent: f64,
    /// Strand diameter \[m\].
    pub diameter: f64,
    /// Flutter amplitude (random perturbation scale).
    pub flutter_amplitude: f64,
}

impl WindResponseModel {
    /// Create a new wind response model.
    pub fn new(
        air_density: f64,
        cd_normal: f64,
        cd_tangent: f64,
        diameter: f64,
        flutter_amplitude: f64,
    ) -> Self {
        Self {
            air_density,
            cd_normal,
            cd_tangent,
            diameter,
            flutter_amplitude,
        }
    }

    /// Default wind model for typical hair.
    pub fn default_hair() -> Self {
        Self::new(1.225, 1.2, 0.1, 0.0001, 0.01)
    }

    /// Compute aerodynamic force on a segment given wind velocity and tangent.
    ///
    /// Decomposes relative velocity into normal and tangential components
    /// and applies corresponding drag coefficients.
    pub fn segment_force(
        &self,
        wind: [f64; 3],
        velocity: [f64; 3],
        tangent: [f64; 3],
        segment_length: f64,
    ) -> [f64; 3] {
        let v_rel = sub3(wind, velocity);
        let v_rel_mag = len3(v_rel);
        if v_rel_mag < 1e-15 {
            return zero3();
        }

        // Decompose into tangential and normal
        let v_tang_mag = dot3(v_rel, tangent);
        let v_tang = scale3(tangent, v_tang_mag);
        let v_norm = sub3(v_rel, v_tang);
        let v_norm_mag = len3(v_norm);

        // Drag = 0.5 * rho * Cd * A * |v|^2 * dir
        let area_normal = self.diameter * segment_length;
        let area_tangent = PI * self.diameter * self.diameter / 4.0;

        let f_normal = if v_norm_mag > 1e-15 {
            let mag = 0.5 * self.air_density * self.cd_normal * area_normal * v_norm_mag;
            scale3(v_norm, mag / v_norm_mag * v_norm_mag)
        } else {
            zero3()
        };

        let f_tangent = if v_tang_mag.abs() > 1e-15 {
            let mag = 0.5 * self.air_density * self.cd_tangent * area_tangent * v_tang_mag.abs();
            scale3(tangent, mag * v_tang_mag.signum())
        } else {
            zero3()
        };

        add3(f_normal, f_tangent)
    }

    /// Accumulate wind forces on a DER strand.
    pub fn apply_to_strand(
        &self,
        rod: &DiscreteElasticRod,
        wind: [f64; 3],
        forces: &mut [[f64; 3]],
    ) {
        for i in 0..rod.num_edges() {
            let t = rod.tangent(i);
            let vel_avg = scale3(add3(rod.velocities[i], rod.velocities[i + 1]), 0.5);
            let seg_len = rod.rest_lengths[i];
            let f = self.segment_force(wind, vel_avg, t, seg_len);
            let half_f = scale3(f, 0.5);
            forces[i] = add3(forces[i], half_f);
            forces[i + 1] = add3(forces[i + 1], half_f);
        }
    }

    /// Total drag force magnitude on a strand.
    pub fn total_drag(&self, rod: &DiscreteElasticRod, wind: [f64; 3]) -> f64 {
        let mut total = zero3();
        for i in 0..rod.num_edges() {
            let t = rod.tangent(i);
            let vel_avg = scale3(add3(rod.velocities[i], rod.velocities[i + 1]), 0.5);
            let seg_len = rod.rest_lengths[i];
            let f = self.segment_force(wind, vel_avg, t, seg_len);
            total = add3(total, f);
        }
        len3(total)
    }
}

// ===========================================================================
// Hair Rendering Data (Tangent Frames)
// ===========================================================================

/// Per-vertex rendering data for hair (tangent, normal, binormal).
#[derive(Debug, Clone)]
pub struct HairTangentFrame {
    /// Per-vertex tangent vectors.
    pub tangents: Vec<[f64; 3]>,
    /// Per-vertex normal vectors.
    pub normals: Vec<[f64; 3]>,
    /// Per-vertex binormal vectors.
    pub binormals: Vec<[f64; 3]>,
}

impl HairTangentFrame {
    /// Build tangent frames from a list of positions.
    pub fn from_positions(positions: &[[f64; 3]]) -> Self {
        let n = positions.len();
        let mut tangents = vec![zero3(); n];
        let mut normals = vec![zero3(); n];
        let mut binormals = vec![zero3(); n];

        if n < 2 {
            return Self {
                tangents,
                normals,
                binormals,
            };
        }

        // Compute tangents
        for i in 0..n {
            tangents[i] = if i == 0 {
                norm3(sub3(positions[1], positions[0]))
            } else if i == n - 1 {
                norm3(sub3(positions[n - 1], positions[n - 2]))
            } else {
                norm3(sub3(positions[i + 1], positions[i - 1]))
            };
        }

        // Build orthonormal frames via parallel transport
        let up = [0.0, 0.0, 1.0];
        normals[0] = {
            let b = norm3(cross3(tangents[0], up));
            if len3(b) < 1e-10 {
                let arb = if tangents[0][0].abs() < 0.9 {
                    [1.0, 0.0, 0.0]
                } else {
                    [0.0, 1.0, 0.0]
                };
                norm3(cross3(norm3(cross3(tangents[0], arb)), tangents[0]))
            } else {
                norm3(cross3(b, tangents[0]))
            }
        };
        binormals[0] = norm3(cross3(tangents[0], normals[0]));

        for i in 1..n {
            let b = cross3(tangents[i - 1], tangents[i]);
            let b_len = len3(b);
            if b_len < 1e-12 {
                normals[i] = normals[i - 1];
            } else {
                let axis = norm3(b);
                let cos_a = clamp_f64(dot3(tangents[i - 1], tangents[i]), -1.0, 1.0);
                let angle = cos_a.acos();
                normals[i] = norm3(rotate_around_axis(normals[i - 1], axis, angle));
            }
            binormals[i] = norm3(cross3(tangents[i], normals[i]));
        }

        Self {
            tangents,
            normals,
            binormals,
        }
    }

    /// Number of vertices with frame data.
    pub fn len(&self) -> usize {
        self.tangents.len()
    }

    /// Check if empty.
    pub fn is_empty(&self) -> bool {
        self.tangents.is_empty()
    }

    /// Get the frame at vertex i as (tangent, normal, binormal).
    pub fn frame_at(&self, i: usize) -> ([f64; 3], [f64; 3], [f64; 3]) {
        (self.tangents[i], self.normals[i], self.binormals[i])
    }
}

// ===========================================================================
// Guide Hair + Interpolation
// ===========================================================================

/// Guide-hair interpolation system.
///
/// A small number of "guide" strands are simulated, and many more
/// "interpolated" strands are generated by blending nearby guides.
#[derive(Debug, Clone)]
pub struct GuideHairInterpolator {
    /// Root positions of guide strands.
    pub guide_roots: Vec<[f64; 3]>,
    /// Number of vertices per guide strand.
    pub verts_per_strand: usize,
    /// Interpolated strand root positions.
    pub interp_roots: Vec<[f64; 3]>,
    /// For each interpolated strand: (guide_index_0, guide_index_1, weight).
    pub interp_weights: Vec<(usize, usize, f64)>,
}

impl GuideHairInterpolator {
    /// Create an interpolator from guide roots and interpolated roots.
    ///
    /// Automatically computes nearest-two-guide blending weights.
    pub fn new(
        guide_roots: Vec<[f64; 3]>,
        verts_per_strand: usize,
        interp_roots: Vec<[f64; 3]>,
    ) -> Self {
        let mut interp_weights = Vec::with_capacity(interp_roots.len());
        for root in &interp_roots {
            let (i0, i1, w) = Self::find_nearest_two(&guide_roots, *root);
            interp_weights.push((i0, i1, w));
        }
        Self {
            guide_roots,
            verts_per_strand,
            interp_roots,
            interp_weights,
        }
    }

    /// Find the two nearest guides and compute blend weight.
    fn find_nearest_two(guides: &[[f64; 3]], point: [f64; 3]) -> (usize, usize, f64) {
        if guides.len() < 2 {
            return (0, 0, 0.5);
        }
        let mut dists: Vec<(usize, f64)> = guides
            .iter()
            .enumerate()
            .map(|(i, g)| (i, dist3(*g, point)))
            .collect();
        dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let i0 = dists[0].0;
        let i1 = dists[1].0;
        let d0 = dists[0].1;
        let d1 = dists[1].1;
        let total = d0 + d1;
        let w = if total < 1e-15 { 0.5 } else { d1 / total };
        (i0, i1, w)
    }

    /// Interpolate positions for one interpolated strand from guide data.
    ///
    /// `guide_positions` is a flat array: guide_positions\[guide_idx * verts_per_strand + vert_idx\].
    pub fn interpolate_strand(
        &self,
        interp_idx: usize,
        guide_positions: &[[f64; 3]],
    ) -> Vec<[f64; 3]> {
        let (i0, i1, w) = self.interp_weights[interp_idx];
        let base0 = i0 * self.verts_per_strand;
        let base1 = i1 * self.verts_per_strand;
        let mut result = Vec::with_capacity(self.verts_per_strand);

        // Offset from guide root to interpolated root
        let root_offset = sub3(
            self.interp_roots[interp_idx],
            lerp3(self.guide_roots[i0], self.guide_roots[i1], 1.0 - w),
        );

        for v in 0..self.verts_per_strand {
            let p0 = guide_positions[base0 + v];
            let p1 = guide_positions[base1 + v];
            let blended = lerp3(p0, p1, 1.0 - w);
            result.push(add3(blended, root_offset));
        }
        result
    }

    /// Number of guide strands.
    pub fn num_guides(&self) -> usize {
        self.guide_roots.len()
    }

    /// Number of interpolated strands.
    pub fn num_interpolated(&self) -> usize {
        self.interp_roots.len()
    }

    /// Total number of strands (guides + interpolated).
    pub fn total_strands(&self) -> usize {
        self.num_guides() + self.num_interpolated()
    }
}

// ===========================================================================
// LOD for Hair Bundles
// ===========================================================================

/// Level of detail for hair bundles.
///
/// Manages multiple LOD levels that reduce strand count and vertex count
/// at increasing distances from the camera.
#[derive(Debug, Clone)]
pub struct HairLod {
    /// LOD levels: (max_distance, strand_fraction, vertex_skip).
    pub levels: Vec<LodLevel>,
}

/// A single LOD level specification.
#[derive(Debug, Clone, Copy)]
pub struct LodLevel {
    /// Maximum camera distance for this level.
    pub max_distance: f64,
    /// Fraction of strands to render \[0..1\].
    pub strand_fraction: f64,
    /// Vertex decimation: keep every N-th vertex (1 = full detail).
    pub vertex_stride: usize,
}

impl LodLevel {
    /// Create a new LOD level.
    pub fn new(max_distance: f64, strand_fraction: f64, vertex_stride: usize) -> Self {
        Self {
            max_distance,
            strand_fraction: clamp_f64(strand_fraction, 0.0, 1.0),
            vertex_stride: vertex_stride.max(1),
        }
    }
}

impl HairLod {
    /// Create a default 4-level LOD system.
    pub fn default_lod() -> Self {
        Self {
            levels: vec![
                LodLevel::new(5.0, 1.0, 1),      // Full detail
                LodLevel::new(15.0, 0.5, 2),     // Half strands, half verts
                LodLevel::new(30.0, 0.25, 4),    // Quarter strands
                LodLevel::new(f64::MAX, 0.1, 8), // Distant
            ],
        }
    }

    /// Create with custom levels.
    pub fn new(levels: Vec<LodLevel>) -> Self {
        Self { levels }
    }

    /// Select the appropriate LOD level for a given camera distance.
    pub fn select_level(&self, distance: f64) -> &LodLevel {
        for level in &self.levels {
            if distance <= level.max_distance {
                return level;
            }
        }
        self.levels.last().unwrap_or(&LodLevel {
            max_distance: f64::MAX,
            strand_fraction: 0.1,
            vertex_stride: 8,
        })
    }

    /// Compute how many strands to render from a total count at a given distance.
    pub fn strand_count(&self, total_strands: usize, distance: f64) -> usize {
        let level = self.select_level(distance);
        ((total_strands as f64) * level.strand_fraction).ceil() as usize
    }

    /// Compute how many vertices per strand at a given distance,
    /// given the original vertex count.
    pub fn vertex_count(&self, original_verts: usize, distance: f64) -> usize {
        let level = self.select_level(distance);
        let count = original_verts.div_ceil(level.vertex_stride);
        count.max(2) // at least 2 vertices (root + tip)
    }

    /// Decimate a position array by the stride of the selected LOD level.
    pub fn decimate_positions(&self, positions: &[[f64; 3]], distance: f64) -> Vec<[f64; 3]> {
        let level = self.select_level(distance);
        let mut result = Vec::new();
        for (i, p) in positions.iter().enumerate() {
            if i % level.vertex_stride == 0 {
                result.push(*p);
            }
        }
        // Always include last vertex (tip)
        if positions.len() > 1 {
            let last = *positions.last().expect("collection should not be empty");
            if result.last() != Some(&last) {
                result.push(last);
            }
        }
        result
    }

    /// Number of LOD levels.
    pub fn num_levels(&self) -> usize {
        self.levels.len()
    }
}

// ===========================================================================
// Hair Bundle (top-level container)
// ===========================================================================

/// A bundle of hair strands with simulation and rendering support.
#[derive(Debug, Clone)]
pub struct HairBundle {
    /// Simulated guide strands.
    pub guides: Vec<DiscreteElasticRod>,
    /// Guide-to-interpolated mapping.
    pub interpolator: Option<GuideHairInterpolator>,
    /// LOD manager.
    pub lod: HairLod,
    /// Wind model.
    pub wind: WindResponseModel,
    /// Global gravity.
    pub gravity: [f64; 3],
}

impl HairBundle {
    /// Create a new hair bundle with given guide strands.
    pub fn new(guides: Vec<DiscreteElasticRod>) -> Self {
        Self {
            guides,
            interpolator: None,
            lod: HairLod::default_lod(),
            wind: WindResponseModel::default_hair(),
            gravity: [0.0, -9.81, 0.0],
        }
    }

    /// Step all guide strands forward.
    pub fn step(&mut self, dt: f64, wind_velocity: [f64; 3]) {
        for guide in &mut self.guides {
            let n = guide.num_vertices();
            let mut forces = vec![zero3(); n];
            self.wind.apply_to_strand(guide, wind_velocity, &mut forces);

            // Add forces to velocities before stepping
            for i in 0..n {
                if guide.inv_mass[i] > 0.0 {
                    let acc = scale3(forces[i], guide.inv_mass[i]);
                    guide.velocities[i] = add3(guide.velocities[i], scale3(acc, dt));
                }
            }

            guide.step(dt, self.gravity);
        }
    }

    /// Number of guide strands.
    pub fn num_guides(&self) -> usize {
        self.guides.len()
    }

    /// Total vertex count across all guides.
    pub fn total_guide_vertices(&self) -> usize {
        self.guides.iter().map(|g| g.num_vertices()).sum()
    }

    /// Get all guide positions as a flat array.
    pub fn flat_guide_positions(&self) -> Vec<[f64; 3]> {
        let mut result = Vec::new();
        for guide in &self.guides {
            result.extend_from_slice(&guide.positions);
        }
        result
    }
}

// ===========================================================================
// Tests
// ===========================================================================

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-6;

    // ---- Vector helpers ----

    #[test]
    fn test_vec3_ops() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let s = add3(a, b);
        assert!((s[0] - 5.0).abs() < TOL);
        assert!((s[1] - 7.0).abs() < TOL);
        assert!((s[2] - 9.0).abs() < TOL);
    }

    #[test]
    fn test_cross_product() {
        let x = [1.0, 0.0, 0.0];
        let y = [0.0, 1.0, 0.0];
        let z = cross3(x, y);
        assert!((z[0]).abs() < TOL);
        assert!((z[1]).abs() < TOL);
        assert!((z[2] - 1.0).abs() < TOL);
    }

    // ---- Discrete Elastic Rod ----

    #[test]
    fn test_der_creation() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            10,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        assert_eq!(rod.num_vertices(), 11);
        assert_eq!(rod.num_edges(), 10);
    }

    #[test]
    fn test_der_root_fixed() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        assert!((rod.inv_mass[0]).abs() < TOL, "Root should be fixed");
        assert!(rod.inv_mass[1] > 0.0, "Non-root should be dynamic");
    }

    #[test]
    fn test_der_rest_lengths() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -2.0, 0.0],
            4,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        for l in &rod.rest_lengths {
            assert!((*l - 0.5).abs() < TOL, "Rest length should be 0.5, got {l}");
        }
    }

    #[test]
    fn test_der_stretch_energy_at_rest() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        let e = rod.stretch_energy();
        assert!(
            e.abs() < TOL,
            "Stretch energy at rest should be ~0, got {e}"
        );
    }

    #[test]
    fn test_der_bend_energy_straight() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        let e = rod.bend_energy();
        assert!(
            e.abs() < TOL,
            "Bend energy of straight rod should be ~0, got {e}"
        );
    }

    #[test]
    fn test_der_step_gravity() {
        let mut rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.1,
            10.0,
            0.01,
            0.01,
            5.0,
        );
        let initial_tip_y = rod.tip()[1];
        let dt = 1.0 / 120.0;
        for _ in 0..20 {
            rod.step(dt, [0.0, -9.81, 0.0]);
        }
        assert!(
            rod.tip()[1] < initial_tip_y,
            "Tip should fall under gravity: initial={initial_tip_y}, now={}",
            rod.tip()[1]
        );
    }

    #[test]
    fn test_der_root_stays_fixed() {
        let mut rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        for _ in 0..20 {
            rod.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        }
        assert!(
            dist3(rod.root(), [0.0, 0.0, 0.0]) < TOL,
            "Root should stay at origin"
        );
    }

    #[test]
    fn test_der_total_length() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -2.0, 0.0],
            10,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        let total = rod.total_length();
        assert!(
            (total - 2.0).abs() < TOL,
            "Total length should be 2.0, got {total}"
        );
    }

    // ---- Cosserat Strand ----

    #[test]
    fn test_cosserat_creation() {
        let strand = CosseratHairStrand::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            8,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
            0.00005,
            80e9,
        );
        assert_eq!(strand.num_vertices(), 9);
        assert_eq!(strand.frames.len(), 8);
    }

    #[test]
    fn test_cosserat_twist_energy_zero() {
        let strand = CosseratHairStrand::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
            0.00005,
            80e9,
        );
        let te = strand.twist_energy();
        assert!(te.abs() < TOL, "Twist energy should be 0 at rest, got {te}");
    }

    #[test]
    fn test_cosserat_frame_update() {
        let mut strand = CosseratHairStrand::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
            0.00005,
            80e9,
        );
        strand.update_frames();
        // All tangents should point downward
        for f in &strand.frames {
            assert!(
                f.tangent[1] < -0.9,
                "Tangent should point -Y, got {:?}",
                f.tangent
            );
        }
    }

    #[test]
    fn test_material_frame_orthogonal() {
        let frame = MaterialFrame::from_tangent_and_up([0.0, -1.0, 0.0], [0.0, 0.0, 1.0]);
        let d_tn = dot3(frame.tangent, frame.normal).abs();
        let d_tb = dot3(frame.tangent, frame.binormal).abs();
        let d_nb = dot3(frame.normal, frame.binormal).abs();
        assert!(d_tn < TOL, "T.N should be 0, got {d_tn}");
        assert!(d_tb < TOL, "T.B should be 0, got {d_tb}");
        assert!(d_nb < TOL, "N.B should be 0, got {d_nb}");
    }

    // ---- Hair-Hair Interaction ----

    #[test]
    fn test_hair_interaction_repulsion() {
        let hhi = HairHairInteraction::new(0.01, 1000.0, 0.3);
        let a = [0.0, 0.0, 0.0];
        let b = [0.005, 0.0, 0.0]; // within collision radius
        let (fa, fb) = hhi.point_point_force(a, b);
        // fa should push a away from b (negative x)
        assert!(
            fa[0] < 0.0,
            "Force on a should be negative x, got {:.6}",
            fa[0]
        );
        // fb should push b away from a (positive x)
        assert!(
            fb[0] > 0.0,
            "Force on b should be positive x, got {:.6}",
            fb[0]
        );
    }

    #[test]
    fn test_hair_interaction_no_force_far() {
        let hhi = HairHairInteraction::new(0.01, 1000.0, 0.3);
        let a = [0.0, 0.0, 0.0];
        let b = [1.0, 0.0, 0.0]; // far apart
        let (fa, fb) = hhi.point_point_force(a, b);
        assert!(len3(fa) < TOL, "No force when far apart");
        assert!(len3(fb) < TOL, "No force when far apart");
    }

    #[test]
    fn test_hair_collision_count() {
        let hhi = HairHairInteraction::new(0.1, 1000.0, 0.3);
        let a = vec![[0.0, 0.0, 0.0], [0.0, -0.5, 0.0], [0.0, -1.0, 0.0]];
        let b = vec![[0.05, 0.0, 0.0], [0.05, -0.5, 0.0], [2.0, -1.0, 0.0]];
        let count = hhi.count_collisions(&a, &b);
        assert_eq!(count, 2, "Expected 2 collisions, got {count}");
    }

    // ---- Styling Force ----

    #[test]
    fn test_styling_force_pulls_to_target() {
        let targets = vec![[1.0, 0.0, 0.0], [1.0, -0.5, 0.0]];
        let sf = StylingForce::new(targets.clone(), 10.0);
        let f = sf.force_at(0, [0.0, 0.0, 0.0]);
        assert!(f[0] > 0.0, "Force should pull toward target x=1");
    }

    #[test]
    fn test_styling_energy_zero_at_target() {
        let targets = vec![[1.0, 2.0, 3.0]];
        let sf = StylingForce::new(targets.clone(), 10.0);
        let e = sf.energy(&targets);
        assert!(e.abs() < TOL, "Energy should be 0 at target, got {e}");
    }

    // ---- Gravity Draping ----

    #[test]
    fn test_draper_converges() {
        let mut rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -0.5, 0.0],
            5,
            0.01,
            1000.0,
            1.0,
            0.01,
            5.0,
        );
        let draper = GravityDraper::with_params([0.0, -9.81, 0.0], 0.005, 5.0, 500, 1e-4);
        let iters = draper.drape(&mut rod);
        assert!(iters <= 500, "Draping should converge within 500 iters");
        // All velocities should be zeroed after draping
        for v in &rod.velocities {
            assert!(len3(*v) < TOL, "Velocities should be zeroed after drape");
        }
    }

    // ---- Wind Response ----

    #[test]
    fn test_wind_zero_no_force() {
        let wind_model = WindResponseModel::default_hair();
        let f = wind_model.segment_force(zero3(), zero3(), [0.0, -1.0, 0.0], 0.1);
        assert!(len3(f) < TOL, "Zero wind should give zero force");
    }

    #[test]
    fn test_wind_produces_force() {
        let wind_model = WindResponseModel::default_hair();
        let wind = [5.0, 0.0, 0.0]; // wind blowing in +X
        let tangent = [0.0, -1.0, 0.0]; // strand hangs down
        let f = wind_model.segment_force(wind, zero3(), tangent, 0.1);
        assert!(
            f[0] > 0.0,
            "Wind in +X should produce +X force, got {:.6}",
            f[0]
        );
    }

    // ---- Tangent Frame ----

    #[test]
    fn test_tangent_frame_straight() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.0, -0.5, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, -1.5, 0.0],
        ];
        let frames = HairTangentFrame::from_positions(&positions);
        assert_eq!(frames.len(), 4);
        for t in &frames.tangents {
            // Should all point in -Y direction
            assert!(t[1] < -0.9, "Tangent should point -Y, got {:?}", t);
        }
    }

    #[test]
    fn test_tangent_frame_orthogonality() {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [0.0, -0.5, 0.1],
            [0.0, -1.0, 0.0],
            [0.1, -1.5, 0.0],
        ];
        let frames = HairTangentFrame::from_positions(&positions);
        for i in 0..frames.len() {
            let (t, n, b) = frames.frame_at(i);
            let d_tn = dot3(t, n).abs();
            let d_tb = dot3(t, b).abs();
            let d_nb = dot3(n, b).abs();
            assert!(d_tn < 0.01, "T.N should be ~0 at vertex {i}, got {d_tn}");
            assert!(d_tb < 0.01, "T.B should be ~0 at vertex {i}, got {d_tb}");
            assert!(d_nb < 0.01, "N.B should be ~0 at vertex {i}, got {d_nb}");
        }
    }

    // ---- Guide Hair Interpolation ----

    #[test]
    fn test_guide_interpolation_exact() {
        let guides = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let guide_positions = vec![
            // Guide 0: 3 verts
            [0.0, 0.0, 0.0],
            [0.0, -0.5, 0.0],
            [0.0, -1.0, 0.0],
            // Guide 1: 3 verts
            [1.0, 0.0, 0.0],
            [1.0, -0.5, 0.0],
            [1.0, -1.0, 0.0],
        ];
        let interp_roots = vec![[0.5, 0.0, 0.0]]; // midpoint
        let interp = GuideHairInterpolator::new(guides, 3, interp_roots);
        let result = interp.interpolate_strand(0, &guide_positions);
        assert_eq!(result.len(), 3);
        // Midpoint interpolation: x should be ~0.5
        assert!(
            (result[1][0] - 0.5).abs() < 0.1,
            "Interpolated x should be ~0.5, got {:.6}",
            result[1][0]
        );
    }

    #[test]
    fn test_guide_count() {
        let guides = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let interp_roots = vec![[0.5, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let interp = GuideHairInterpolator::new(guides, 5, interp_roots);
        assert_eq!(interp.num_guides(), 3);
        assert_eq!(interp.num_interpolated(), 2);
        assert_eq!(interp.total_strands(), 5);
    }

    // ---- LOD ----

    #[test]
    fn test_lod_close() {
        let lod = HairLod::default_lod();
        let count = lod.strand_count(100, 3.0);
        assert_eq!(count, 100, "Full strands at close distance");
    }

    #[test]
    fn test_lod_far() {
        let lod = HairLod::default_lod();
        let count = lod.strand_count(100, 20.0);
        assert!(count < 100, "Fewer strands at medium distance: got {count}");
    }

    #[test]
    fn test_lod_vertex_decimation() {
        let lod = HairLod::default_lod();
        let verts = lod.vertex_count(20, 20.0);
        assert!(
            verts < 20,
            "Should decimate vertices at distance: got {verts}"
        );
        assert!(verts >= 2, "At least 2 vertices");
    }

    #[test]
    fn test_lod_decimate_positions() {
        let lod = HairLod::default_lod();
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [0.0, -(i as f64) * 0.1, 0.0]).collect();
        let decimated = lod.decimate_positions(&positions, 20.0);
        assert!(
            decimated.len() < positions.len(),
            "Decimated should have fewer vertices: {} vs {}",
            decimated.len(),
            positions.len()
        );
    }

    // ---- Hair Bundle ----

    #[test]
    fn test_hair_bundle_creation() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        let bundle = HairBundle::new(vec![rod]);
        assert_eq!(bundle.num_guides(), 1);
        assert_eq!(bundle.total_guide_vertices(), 6);
    }

    #[test]
    fn test_hair_bundle_step() {
        let rod = DiscreteElasticRod::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            5,
            0.01,
            100.0,
            0.1,
            0.01,
            0.5,
        );
        let mut bundle = HairBundle::new(vec![rod]);
        let initial_tip = bundle.guides[0].tip();
        bundle.step(1.0 / 60.0, [0.0, 0.0, 0.0]);
        let new_tip = bundle.guides[0].tip();
        assert!(
            dist3(initial_tip, new_tip) > 1e-10,
            "Tip should move after stepping"
        );
    }

    // ---- Rodrigues rotation ----

    #[test]
    fn test_rodrigues_identity() {
        let v = [1.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let rotated = rotate_around_axis(v, axis, 0.0);
        assert!(
            (rotated[0] - 1.0).abs() < TOL,
            "No rotation should preserve vector"
        );
        assert!(rotated[1].abs() < TOL);
    }

    #[test]
    fn test_rodrigues_90_deg() {
        let v = [1.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let rotated = rotate_around_axis(v, axis, PI / 2.0);
        assert!(
            rotated[0].abs() < TOL,
            "x should be ~0 after 90 deg rotation"
        );
        assert!(
            (rotated[1] - 1.0).abs() < TOL,
            "y should be ~1 after 90 deg rotation"
        );
    }
}
