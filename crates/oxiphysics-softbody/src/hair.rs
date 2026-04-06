// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hair simulation.
//!
//! Provides strand-based hair physics including stretch/bend/twist constraints,
//! collision response, guide interpolation, shading models, clustering, grooming,
//! and aerodynamics.

#[allow(unused_imports)]
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Small vector helpers (no nalgebra)
// ---------------------------------------------------------------------------

#[allow(dead_code)]
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[allow(dead_code)]
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[allow(dead_code)]
#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[allow(dead_code)]
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[allow(dead_code)]
#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[allow(dead_code)]
#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-12 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

#[allow(dead_code)]
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ---------------------------------------------------------------------------
// SimHairStrand — a single hair strand (segments ordered tip-to-root)
// ---------------------------------------------------------------------------

/// A single hair strand composed of linked segments.
///
/// Segments are stored in root-to-tip order.
#[allow(dead_code)]
pub struct SimHairStrand {
    /// World-space positions of each particle along the strand.
    pub positions: Vec<[f64; 3]>,
    /// World-space velocities of each particle.
    pub velocities: Vec<[f64; 3]>,
    /// Rest length between consecutive particles (metres).
    pub rest_length: f64,
    /// Rest direction (unit vector) for bending constraints.
    pub rest_direction: [f64; 3],
    /// Stiffness coefficient \[0, 1\].
    pub stiffness: f64,
    /// Whether the root particle is pinned (fixed).
    pub root_pinned: bool,
}

impl SimHairStrand {
    /// Create a strand with `segments` particles along `growth_dir`.
    pub fn new(
        root: [f64; 3],
        segment_count: usize,
        rest_length: f64,
        growth_dir: [f64; 3],
    ) -> Self {
        let dir = norm3(growth_dir);
        let mut positions = Vec::with_capacity(segment_count);
        for i in 0..segment_count {
            positions.push(add3(root, scale3(dir, i as f64 * rest_length)));
        }
        let velocities = vec![[0.0; 3]; segment_count];
        Self {
            positions,
            velocities,
            rest_length,
            rest_direction: dir,
            stiffness: 0.9,
            root_pinned: true,
        }
    }

    /// Number of particles in the strand.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// Returns `true` if the strand has no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Tip position (last particle).
    pub fn tip(&self) -> [f64; 3] {
        *self.positions.last().unwrap_or(&[0.0; 3])
    }

    /// Root position (first particle).
    pub fn root(&self) -> [f64; 3] {
        *self.positions.first().unwrap_or(&[0.0; 3])
    }
}

// ---------------------------------------------------------------------------
// HairFollicle
// ---------------------------------------------------------------------------

/// Root attachment point for a hair strand on the scalp mesh.
#[allow(dead_code)]
pub struct HairFollicle {
    /// World-space position of the follicle on the scalp.
    pub position: [f64; 3],
    /// Local frame: tangent, bitangent, normal.
    pub frame: [[f64; 3]; 3],
    /// Growth direction (local normal direction).
    pub growth_dir: [f64; 3],
    /// Index of the strand this follicle spawns.
    pub strand_index: usize,
}

impl HairFollicle {
    /// Create a follicle at `position` with `normal` as the growth direction.
    pub fn new(position: [f64; 3], normal: [f64; 3], strand_index: usize) -> Self {
        let n = norm3(normal);
        // Build tangent frame
        let t = if n[0].abs() < 0.9 {
            norm3(cross3(n, [1.0, 0.0, 0.0]))
        } else {
            norm3(cross3(n, [0.0, 1.0, 0.0]))
        };
        let bt = cross3(n, t);
        Self {
            position,
            frame: [t, bt, n],
            growth_dir: n,
            strand_index,
        }
    }
}

// ---------------------------------------------------------------------------
// HairConstraint
// ---------------------------------------------------------------------------

/// Hair constraint types.
#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HairConstraintKind {
    /// Distance (stretch) constraint between two particles.
    Stretch,
    /// Angular (bend) constraint between three consecutive particles.
    Bend,
    /// Twist constraint around the strand axis.
    Twist,
}

/// A single hair constraint.
#[allow(dead_code)]
pub struct HairConstraint {
    /// Constraint type.
    pub kind: HairConstraintKind,
    /// Index of first particle.
    pub i: usize,
    /// Index of second particle.
    pub j: usize,
    /// Rest value (length or angle).
    pub rest: f64,
    /// Compliance (XPBD).
    pub compliance: f64,
}

impl HairConstraint {
    /// Create a stretch constraint.
    pub fn stretch(i: usize, j: usize, rest_length: f64) -> Self {
        Self {
            kind: HairConstraintKind::Stretch,
            i,
            j,
            rest: rest_length,
            compliance: 1e-8,
        }
    }

    /// Create a bend constraint.
    pub fn bend(i: usize, j: usize, rest_angle: f64) -> Self {
        Self {
            kind: HairConstraintKind::Bend,
            i,
            j,
            rest: rest_angle,
            compliance: 1e-4,
        }
    }

    /// Create a twist constraint.
    pub fn twist(i: usize, j: usize) -> Self {
        Self {
            kind: HairConstraintKind::Twist,
            i,
            j,
            rest: 0.0,
            compliance: 1e-5,
        }
    }
}

// ---------------------------------------------------------------------------
// HairCollision
// ---------------------------------------------------------------------------

/// Capsule-vs-strand collision detector.
#[allow(dead_code)]
pub struct HairCollision {
    /// Capsule axis start.
    pub capsule_a: [f64; 3],
    /// Capsule axis end.
    pub capsule_b: [f64; 3],
    /// Capsule radius (metres).
    pub capsule_radius: f64,
    /// Restitution coefficient for strand bounce.
    pub restitution: f64,
}

impl HairCollision {
    /// Create a hair collision capsule.
    pub fn new(a: [f64; 3], b: [f64; 3], radius: f64) -> Self {
        Self {
            capsule_a: a,
            capsule_b: b,
            capsule_radius: radius,
            restitution: 0.1,
        }
    }

    /// Closest point on capsule axis to point `p`.
    pub fn closest_on_axis(&self, p: [f64; 3]) -> [f64; 3] {
        let ab = sub3(self.capsule_b, self.capsule_a);
        let ap = sub3(p, self.capsule_a);
        let t = (dot3(ap, ab) / dot3(ab, ab).max(1e-12)).clamp(0.0, 1.0);
        add3(self.capsule_a, scale3(ab, t))
    }

    /// Compute impulse correction for a strand particle at `pos` with `vel`.
    ///
    /// Returns displacement vector to resolve collision (zero if no collision).
    pub fn resolve_particle(&self, pos: [f64; 3]) -> [f64; 3] {
        let closest = self.closest_on_axis(pos);
        let diff = sub3(pos, closest);
        let dist = len3(diff);
        if dist < self.capsule_radius && dist > 1e-12 {
            let penetration = self.capsule_radius - dist;
            scale3(norm3(diff), penetration)
        } else {
            [0.0; 3]
        }
    }

    /// Apply collision response to all particles in a strand.
    pub fn apply_to_strand(&self, strand: &mut SimHairStrand) {
        let start = if strand.root_pinned { 1 } else { 0 };
        for i in start..strand.positions.len() {
            let correction = self.resolve_particle(strand.positions[i]);
            if len3(correction) > 1e-12 {
                strand.positions[i] = add3(strand.positions[i], correction);
                // Damp velocity component into capsule
                let n = norm3(correction);
                let vn = dot3(strand.velocities[i], n);
                if vn < 0.0 {
                    strand.velocities[i] = add3(
                        strand.velocities[i],
                        scale3(n, -vn * (1.0 + self.restitution)),
                    );
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// HairSimulation
// ---------------------------------------------------------------------------

/// Full hair simulation: all strands, gravity, wind, damping.
#[allow(dead_code)]
pub struct HairSimulation {
    /// All hair strands.
    pub strands: Vec<SimHairStrand>,
    /// Gravity vector.
    pub gravity: [f64; 3],
    /// Wind force vector.
    pub wind: [f64; 3],
    /// Velocity damping coefficient \[0, 1\].
    pub damping: f64,
    /// Number of constraint solver iterations per step.
    pub solver_iterations: u32,
    /// Whether self-collision is enabled.
    pub self_collision: bool,
}

impl HairSimulation {
    /// Create a hair simulation with default parameters.
    pub fn new() -> Self {
        Self {
            strands: Vec::new(),
            gravity: [0.0, -9.81, 0.0],
            wind: [0.0; 3],
            damping: 0.98,
            solver_iterations: 4,
            self_collision: false,
        }
    }

    /// Add a strand.
    pub fn add_strand(&mut self, strand: SimHairStrand) {
        self.strands.push(strand);
    }

    /// Integrate all strands forward by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let gravity = self.gravity;
        let wind = self.wind;
        let damping = self.damping;
        let iters = self.solver_iterations;
        for strand in &mut self.strands {
            let start = if strand.root_pinned { 1 } else { 0 };
            // Apply external forces
            for i in start..strand.positions.len() {
                let f = add3(gravity, wind);
                strand.velocities[i] = add3(strand.velocities[i], scale3(f, dt));
                strand.velocities[i] = scale3(strand.velocities[i], damping);
                strand.positions[i] = add3(strand.positions[i], scale3(strand.velocities[i], dt));
            }
            // Solve stretch constraints
            for _ in 0..iters {
                for seg in 0..strand.positions.len().saturating_sub(1) {
                    let diff = sub3(strand.positions[seg + 1], strand.positions[seg]);
                    let dist = len3(diff);
                    if dist < 1e-12 {
                        continue;
                    }
                    let correction = scale3(diff, (dist - strand.rest_length) / dist * 0.5);
                    if seg >= start {
                        strand.positions[seg] = add3(strand.positions[seg], correction);
                    }
                    if seg + 1 >= start {
                        strand.positions[seg + 1] = sub3(strand.positions[seg + 1], correction);
                    }
                }
            }
        }
    }
}

impl Default for HairSimulation {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// HairGuideInterp
// ---------------------------------------------------------------------------

/// Guide-hair interpolation system.
///
/// Simulates a smaller set of guide hairs and interpolates the full
/// resolution from them using barycentric blending.
#[allow(dead_code)]
pub struct HairGuideInterp {
    /// Guide strands (simulated at full fidelity).
    pub guides: Vec<SimHairStrand>,
    /// Weights for each full-res strand: list of (guide_index, weight) pairs.
    pub strand_weights: Vec<Vec<(usize, f64)>>,
    /// Number of full-res strands.
    pub full_strand_count: usize,
}

impl HairGuideInterp {
    /// Create a guide interpolation system.
    pub fn new(full_strand_count: usize) -> Self {
        Self {
            guides: Vec::new(),
            strand_weights: vec![Vec::new(); full_strand_count],
            full_strand_count,
        }
    }

    /// Add a guide strand.
    pub fn add_guide(&mut self, guide: SimHairStrand) {
        self.guides.push(guide);
    }

    /// Assign interpolation weights for full-res strand `strand_idx`.
    pub fn set_weights(&mut self, strand_idx: usize, weights: Vec<(usize, f64)>) {
        if strand_idx < self.full_strand_count {
            self.strand_weights[strand_idx] = weights;
        }
    }

    /// Interpolate position of particle `particle_idx` on full-res strand `strand_idx`.
    pub fn interpolated_position(&self, strand_idx: usize, particle_idx: usize) -> [f64; 3] {
        if strand_idx >= self.full_strand_count {
            return [0.0; 3];
        }
        let mut result = [0.0_f64; 3];
        let mut total_w = 0.0;
        for &(guide_idx, w) in &self.strand_weights[strand_idx] {
            if let Some(guide) = self.guides.get(guide_idx)
                && let Some(&pos) = guide.positions.get(particle_idx)
            {
                result = add3(result, scale3(pos, w));
                total_w += w;
            }
        }
        if total_w > 1e-12 {
            scale3(result, 1.0 / total_w)
        } else {
            result
        }
    }
}

// ---------------------------------------------------------------------------
// HairShading — Kajiya-Kay model
// ---------------------------------------------------------------------------

/// Kajiya-Kay shading model parameters for hair.
#[allow(dead_code)]
pub struct HairShading {
    /// Diffuse scatter coefficient.
    pub diffuse: f64,
    /// Specular scatter coefficient.
    pub specular: f64,
    /// Specular exponent (shininess).
    pub shininess: f64,
    /// Transmission fraction.
    pub transmission: f64,
    /// Anisotropy angle shift (radians).
    pub anisotropy_shift: f64,
}

impl HairShading {
    /// Create Kajiya-Kay parameters.
    pub fn new(diffuse: f64, specular: f64, shininess: f64) -> Self {
        Self {
            diffuse,
            specular,
            shininess,
            transmission: 0.3,
            anisotropy_shift: -0.1,
        }
    }

    /// Compute the Kajiya-Kay diffuse term.
    ///
    /// `light_dir` and `tangent` are unit vectors.
    pub fn diffuse_term(&self, light_dir: [f64; 3], tangent: [f64; 3]) -> f64 {
        let cos_theta = dot3(light_dir, tangent).clamp(-1.0, 1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        self.diffuse * sin_theta
    }

    /// Compute the Kajiya-Kay specular highlight.
    ///
    /// `light_dir`, `view_dir`, `tangent` are unit vectors.
    pub fn specular_term(&self, light_dir: [f64; 3], view_dir: [f64; 3], tangent: [f64; 3]) -> f64 {
        kajiya_kay_highlight(
            light_dir,
            view_dir,
            tangent,
            self.shininess,
            self.anisotropy_shift,
        ) * self.specular
    }

    /// Combined shading.
    pub fn shade(&self, light_dir: [f64; 3], view_dir: [f64; 3], tangent: [f64; 3]) -> f64 {
        self.diffuse_term(light_dir, tangent) + self.specular_term(light_dir, view_dir, tangent)
    }
}

/// Compute Kajiya-Kay specular highlight.
#[allow(dead_code)]
pub fn kajiya_kay_highlight(
    light_dir: [f64; 3],
    view_dir: [f64; 3],
    tangent: [f64; 3],
    shininess: f64,
    shift: f64,
) -> f64 {
    let t_shifted = norm3(add3(tangent, scale3([0.0, 1.0, 0.0], shift)));
    let dot_lt = dot3(light_dir, t_shifted).clamp(-1.0, 1.0);
    let dot_vt = dot3(view_dir, t_shifted).clamp(-1.0, 1.0);
    let sin_lt = (1.0 - dot_lt * dot_lt).sqrt();
    let sin_vt = (1.0 - dot_vt * dot_vt).sqrt();
    let h = sin_lt * sin_vt - dot_lt * dot_vt;
    h.max(0.0).powf(shininess)
}

// ---------------------------------------------------------------------------
// HairClustering
// ---------------------------------------------------------------------------

/// Hair clustering for LOD and wind stiffness.
#[allow(dead_code)]
pub struct HairClustering {
    /// Cluster centres (world-space XZ position of cluster root).
    pub cluster_centres: Vec<[f64; 2]>,
    /// Strand-to-cluster assignment.
    pub assignments: Vec<usize>,
    /// Wind stiffness per cluster.
    pub wind_stiffness: Vec<f64>,
    /// LOD distance thresholds (metres).
    pub lod_distances: Vec<f64>,
}

impl HairClustering {
    /// Create an empty clustering.
    pub fn new() -> Self {
        Self {
            cluster_centres: Vec::new(),
            assignments: Vec::new(),
            wind_stiffness: Vec::new(),
            lod_distances: vec![5.0, 15.0, 30.0],
        }
    }

    /// Compute k-means clustering of strand roots.
    ///
    /// `roots` are XZ positions of strand roots; `k` is the number of clusters.
    pub fn cluster(&mut self, roots: &[[f64; 2]], k: usize) {
        if roots.is_empty() || k == 0 {
            return;
        }
        let k = k.min(roots.len());
        // Initialise centres from first k roots
        self.cluster_centres = roots[..k].to_vec();
        self.assignments = vec![0; roots.len()];
        self.wind_stiffness = vec![1.0; k];

        for _ in 0..10 {
            // Assign
            for (i, root) in roots.iter().enumerate() {
                let nearest = self
                    .cluster_centres
                    .iter()
                    .enumerate()
                    .min_by(|(_, a), (_, b)| {
                        let da = dist2(*root, **a);
                        let db = dist2(*root, **b);
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .map(|(idx, _)| idx)
                    .unwrap_or(0);
                self.assignments[i] = nearest;
            }
            // Update centres
            let mut sums = vec![[0.0_f64; 2]; k];
            let mut counts = vec![0usize; k];
            for (i, root) in roots.iter().enumerate() {
                let c = self.assignments[i];
                sums[c] = [sums[c][0] + root[0], sums[c][1] + root[1]];
                counts[c] += 1;
            }
            for c in 0..k {
                if counts[c] > 0 {
                    self.cluster_centres[c] =
                        [sums[c][0] / counts[c] as f64, sums[c][1] / counts[c] as f64];
                }
            }
        }
    }

    /// Determine LOD level for a given camera distance.
    pub fn lod_level(&self, distance: f64) -> usize {
        for (i, &thresh) in self.lod_distances.iter().enumerate() {
            if distance < thresh {
                return i;
            }
        }
        self.lod_distances.len()
    }
}

impl Default for HairClustering {
    fn default() -> Self {
        Self::new()
    }
}

fn dist2(a: [f64; 2], b: [f64; 2]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    dx * dx + dy * dy
}

// ---------------------------------------------------------------------------
// HairGroom
// ---------------------------------------------------------------------------

/// Grooming operator for hair strands.
#[allow(dead_code)]
pub struct HairGroom {
    /// Comb direction field (one direction per cluster or control point).
    pub comb_directions: Vec<[f64; 3]>,
    /// Curl radius (metres); 0 = no curl.
    pub curl_radius: f64,
    /// Curl frequency (turns per metre).
    pub curl_frequency: f64,
    /// Cut length (metres; 0 = no cut).
    pub cut_length: f64,
}

impl HairGroom {
    /// Create a groom operator.
    pub fn new() -> Self {
        Self {
            comb_directions: Vec::new(),
            curl_radius: 0.0,
            curl_frequency: 0.0,
            cut_length: 0.0,
        }
    }

    /// Apply combing direction to a strand.
    pub fn comb(&self, strand: &mut SimHairStrand, direction: [f64; 3]) {
        let d = norm3(direction);
        let n = strand.positions.len();
        for i in 1..n {
            let from = strand.positions[i - 1];
            let desired = add3(from, scale3(d, strand.rest_length));
            // Blend current position toward combed position
            let w = 0.5_f64.powi(i as i32);
            strand.positions[i] = add3(scale3(strand.positions[i], 1.0 - w), scale3(desired, w));
        }
    }

    /// Apply curl to a strand.
    pub fn curl(&self, strand: &mut SimHairStrand) {
        if self.curl_radius < 1e-9 || self.curl_frequency < 1e-9 {
            return;
        }
        let r = self.curl_radius;
        let f = self.curl_frequency;
        let root = strand.positions[0];
        for (i, pos) in strand.positions.iter_mut().enumerate().skip(1) {
            let t = i as f64 * strand.rest_length;
            let angle = TAU * f * t;
            *pos = [root[0] + r * angle.cos(), pos[1], root[2] + r * angle.sin()];
        }
    }

    /// Cut a strand to `cut_length`.
    pub fn cut(&self, strand: &mut SimHairStrand) {
        if self.cut_length < 1e-9 {
            return;
        }
        let max_particles = (self.cut_length / strand.rest_length).ceil() as usize + 1;
        strand.positions.truncate(max_particles);
        strand.velocities.truncate(max_particles);
    }
}

impl Default for HairGroom {
    fn default() -> Self {
        Self::new()
    }
}

use std::f64::consts::TAU;

// ---------------------------------------------------------------------------
// HairAerodynamics
// ---------------------------------------------------------------------------

/// Per-strand aerodynamics using drag and lift on each segment.
#[allow(dead_code)]
pub struct HairAerodynamics {
    /// Drag coefficient.
    pub drag_coeff: f64,
    /// Lift coefficient.
    pub lift_coeff: f64,
    /// Air density (kg/m³).
    pub air_density: f64,
    /// Segment cross-sectional area (m²).
    pub segment_area: f64,
}

impl HairAerodynamics {
    /// Create an aerodynamics model.
    pub fn new(drag_coeff: f64, lift_coeff: f64) -> Self {
        Self {
            drag_coeff,
            lift_coeff,
            air_density: 1.225,
            segment_area: 1e-6,
        }
    }

    /// Compute aerodynamic force on a strand segment.
    ///
    /// `wind_vel` is the wind velocity; `tangent` is the segment unit tangent.
    /// Returns a force vector (N) to apply to the midpoint of the segment.
    pub fn segment_force(&self, wind_vel: [f64; 3], tangent: [f64; 3]) -> [f64; 3] {
        // Relative wind velocity projected perpendicular to strand
        let dot_wt = dot3(wind_vel, tangent);
        let wind_perp = sub3(wind_vel, scale3(tangent, dot_wt));
        let v_perp_len = len3(wind_perp);
        if v_perp_len < 1e-9 {
            return [0.0; 3];
        }
        let q = 0.5 * self.air_density * v_perp_len * v_perp_len * self.segment_area;
        let drag = scale3(norm3(wind_perp), self.drag_coeff * q);

        // Lift: perpendicular to both wind_perp and tangent
        let lift_dir = norm3(cross3(wind_perp, tangent));
        let lift = scale3(lift_dir, self.lift_coeff * q);

        add3(drag, lift)
    }

    /// Apply aerodynamic forces to all non-root particles of a strand.
    pub fn apply_to_strand(
        &self,
        strand: &mut SimHairStrand,
        wind_vel: [f64; 3],
        dt: f64,
        inv_mass: f64,
    ) {
        let n = strand.positions.len();
        let start = if strand.root_pinned { 1 } else { 0 };
        for i in start..n {
            let tangent = if i > 0 {
                let seg = sub3(strand.positions[i], strand.positions[i - 1]);
                norm3(seg)
            } else if i + 1 < n {
                norm3(sub3(strand.positions[i + 1], strand.positions[i]))
            } else {
                strand.rest_direction
            };
            let force = self.segment_force(wind_vel, tangent);
            let accel = scale3(force, inv_mass);
            strand.velocities[i] = add3(strand.velocities[i], scale3(accel, dt));
        }
    }
}

// ---------------------------------------------------------------------------
// Helper functions
// ---------------------------------------------------------------------------

/// Compute the stretch (potential) energy of a strand.
#[allow(dead_code)]
pub fn strand_stretch_energy(strand: &SimHairStrand, stiffness: f64) -> f64 {
    let mut energy = 0.0;
    for i in 0..strand.positions.len().saturating_sub(1) {
        let dist = len3(sub3(strand.positions[i + 1], strand.positions[i]));
        let extension = dist - strand.rest_length;
        energy += 0.5 * stiffness * extension * extension;
    }
    energy
}

/// Compute the bend angle between three consecutive particles.
#[allow(dead_code)]
pub fn bend_angle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let v1 = norm3(sub3(a, b));
    let v2 = norm3(sub3(c, b));
    dot3(v1, v2).clamp(-1.0, 1.0).acos()
}

/// Build a local frame for a hair segment given `tangent`.
///
/// Returns `[tangent, normal, binormal]`.
#[allow(dead_code)]
pub fn local_frame_hair(tangent: [f64; 3]) -> [[f64; 3]; 3] {
    let t = norm3(tangent);
    let up = if t[1].abs() < 0.9 {
        [0.0, 1.0, 0.0_f64]
    } else {
        [1.0, 0.0, 0.0_f64]
    };
    let n = norm3(cross3(t, up));
    let b = cross3(t, n);
    [t, n, b]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // 1. SimHairStrand: correct number of particles
    #[test]
    fn test_strand_particle_count() {
        let s = SimHairStrand::new([0.0; 3], 10, 0.05, [0.0, 1.0, 0.0]);
        assert_eq!(s.len(), 10);
    }

    // 2. SimHairStrand: root is pinned
    #[test]
    fn test_strand_root_pinned() {
        let s = SimHairStrand::new([0.0; 3], 5, 0.05, [0.0, 1.0, 0.0]);
        assert!(s.root_pinned);
    }

    // 3. SimHairStrand: tip is at root + (n-1)*rest*dir
    #[test]
    fn test_strand_tip_position() {
        let s = SimHairStrand::new([0.0; 3], 5, 0.1, [0.0, 1.0, 0.0]);
        let tip = s.tip();
        assert!((tip[1] - 0.4).abs() < 1e-9);
    }

    // 4. SimHairStrand: is_empty false
    #[test]
    fn test_strand_not_empty() {
        let s = SimHairStrand::new([0.0; 3], 3, 0.1, [1.0, 0.0, 0.0]);
        assert!(!s.is_empty());
    }

    // 5. HairFollicle: frame is orthonormal
    #[test]
    fn test_follicle_frame_orthogonal() {
        let f = HairFollicle::new([0.0; 3], [0.0, 1.0, 0.0], 0);
        let t = f.frame[0];
        let n = f.frame[2];
        assert!(dot3(t, n).abs() < 1e-9);
    }

    // 6. HairConstraint: stretch compliance lower than bend
    #[test]
    fn test_constraint_compliance_ordering() {
        let s = HairConstraint::stretch(0, 1, 0.05);
        let b = HairConstraint::bend(0, 2, 0.0);
        assert!(s.compliance < b.compliance);
    }

    // 7. HairCollision: no correction when far away
    #[test]
    fn test_collision_no_correction() {
        let hc = HairCollision::new([0.0; 3], [0.0, 1.0, 0.0], 0.1);
        let corr = hc.resolve_particle([1.0, 0.5, 0.0]);
        assert!(len3(corr) < 1e-9);
    }

    // 8. HairCollision: correction when overlapping
    #[test]
    fn test_collision_correction() {
        let hc = HairCollision::new([0.0; 3], [0.0, 1.0, 0.0], 0.2);
        let corr = hc.resolve_particle([0.1, 0.5, 0.0]);
        assert!(len3(corr) > 0.0);
    }

    // 9. HairSimulation: step advances positions
    #[test]
    fn test_simulation_step() {
        let mut sim = HairSimulation::new();
        let s = SimHairStrand::new([0.0; 3], 5, 0.05, [0.0, 1.0, 0.0]);
        sim.add_strand(s);
        let before = sim.strands[0].positions[1];
        sim.step(0.016);
        let after = sim.strands[0].positions[1];
        assert!(before != after);
    }

    // 10. HairSimulation: root stays fixed
    #[test]
    fn test_simulation_root_fixed() {
        let mut sim = HairSimulation::new();
        let s = SimHairStrand::new([0.0; 3], 5, 0.05, [0.0, 1.0, 0.0]);
        sim.add_strand(s);
        let root_before = sim.strands[0].positions[0];
        sim.step(0.1);
        let root_after = sim.strands[0].positions[0];
        assert_eq!(root_before, root_after);
    }

    // 11. HairGuideInterp: interpolated position blends guides
    #[test]
    fn test_guide_interp() {
        let mut hgi = HairGuideInterp::new(1);
        let mut g1 = SimHairStrand::new([0.0; 3], 3, 0.1, [0.0, 1.0, 0.0]);
        g1.positions[1] = [0.0, 0.1, 0.0];
        let mut g2 = SimHairStrand::new([0.0; 3], 3, 0.1, [0.0, 1.0, 0.0]);
        g2.positions[1] = [0.0, 0.3, 0.0];
        hgi.add_guide(g1);
        hgi.add_guide(g2);
        hgi.set_weights(0, vec![(0, 0.5), (1, 0.5)]);
        let p = hgi.interpolated_position(0, 1);
        assert!((p[1] - 0.2).abs() < 1e-9);
    }

    // 12. HairShading: diffuse zero along tangent
    #[test]
    fn test_shading_diffuse_zero() {
        let hs = HairShading::new(1.0, 0.5, 20.0);
        // Light parallel to tangent → sin=0
        let d = hs.diffuse_term([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-9);
    }

    // 13. HairShading: diffuse max perpendicular to tangent
    #[test]
    fn test_shading_diffuse_max() {
        let hs = HairShading::new(1.0, 0.5, 20.0);
        let d = hs.diffuse_term([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 1.0).abs() < 1e-9);
    }

    // 14. kajiya_kay_highlight: non-negative
    #[test]
    fn test_kajiya_kay_non_negative() {
        let h = kajiya_kay_highlight(
            [0.577, 0.577, 0.577],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            20.0,
            0.0,
        );
        assert!(h >= 0.0);
    }

    // 15. HairClustering: assigns all strands
    #[test]
    fn test_clustering_assigns_all() {
        let mut hc = HairClustering::new();
        let roots: Vec<[f64; 2]> = (0..10).map(|i| [i as f64 * 0.1, 0.0]).collect();
        hc.cluster(&roots, 3);
        assert_eq!(hc.assignments.len(), 10);
    }

    // 16. HairClustering: LOD level 0 close up
    #[test]
    fn test_clustering_lod_close() {
        let hc = HairClustering::new();
        assert_eq!(hc.lod_level(1.0), 0);
    }

    // 17. HairClustering: LOD level 3 far away
    #[test]
    fn test_clustering_lod_far() {
        let hc = HairClustering::new();
        assert_eq!(hc.lod_level(100.0), 3);
    }

    // 18. HairGroom: comb modifies non-root positions
    #[test]
    fn test_groom_comb() {
        let groom = HairGroom::new();
        let mut strand = SimHairStrand::new([0.0; 3], 5, 0.1, [0.0, 1.0, 0.0]);
        let before = strand.positions[2];
        groom.comb(&mut strand, [1.0, 0.0, 0.0]);
        assert!(strand.positions[2] != before);
    }

    // 19. HairGroom: cut shortens strand
    #[test]
    fn test_groom_cut() {
        let mut groom = HairGroom::new();
        groom.cut_length = 0.25;
        let mut strand = SimHairStrand::new([0.0; 3], 10, 0.1, [0.0, 1.0, 0.0]);
        groom.cut(&mut strand);
        assert!(strand.positions.len() <= 4);
    }

    // 20. HairAerodynamics: zero force with zero wind
    #[test]
    fn test_aero_zero_wind() {
        let ha = HairAerodynamics::new(1.0, 0.5);
        let f = ha.segment_force([0.0; 3], [0.0, 1.0, 0.0]);
        assert_eq!(f, [0.0; 3]);
    }

    // 21. HairAerodynamics: force non-zero with cross-wind
    #[test]
    fn test_aero_cross_wind() {
        let ha = HairAerodynamics::new(1.0, 0.5);
        let f = ha.segment_force([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(len3(f) > 0.0);
    }

    // 22. strand_stretch_energy: zero for unstretched strand
    #[test]
    fn test_stretch_energy_zero() {
        let s = SimHairStrand::new([0.0; 3], 5, 0.1, [0.0, 1.0, 0.0]);
        let e = strand_stretch_energy(&s, 100.0);
        assert!(e < 1e-9);
    }

    // 23. bend_angle: straight strand has angle PI
    #[test]
    fn test_bend_angle_straight() {
        let a = [0.0, 0.0, 0.0_f64];
        let b = [0.0, 1.0, 0.0_f64];
        let c = [0.0, 2.0, 0.0_f64];
        let angle = bend_angle(a, b, c);
        assert!((angle - PI).abs() < 1e-9);
    }

    // 24. local_frame_hair: tangent preserved
    #[test]
    fn test_local_frame_tangent() {
        let frame = local_frame_hair([0.0, 1.0, 0.0]);
        let t = frame[0];
        assert!((t[1] - 1.0).abs() < 1e-9);
    }

    // 25. local_frame_hair: orthogonal
    #[test]
    fn test_local_frame_ortho() {
        let frame = local_frame_hair([1.0, 0.0, 0.0]);
        assert!(dot3(frame[0], frame[1]).abs() < 1e-9);
        assert!(dot3(frame[0], frame[2]).abs() < 1e-9);
    }

    // 26. HairCollision: apply_to_strand resolves penetrations for non-root particles
    #[test]
    fn test_collision_apply_to_strand() {
        let hc = HairCollision::new([0.0; 3], [0.0, 1.0, 0.0], 0.3);
        let mut strand = SimHairStrand::new([0.0; 3], 5, 0.05, [0.1, 0.5, 0.0]);
        hc.apply_to_strand(&mut strand);
        // Non-root particles (index >= 1) should be outside or on the capsule surface
        for pos in strand.positions.iter().skip(1) {
            let closest = hc.closest_on_axis(*pos);
            let d = len3(sub3(*pos, closest));
            assert!(d >= hc.capsule_radius - 1e-9);
        }
    }

    // 27. HairSimulation: default gravity is -9.81 in y
    #[test]
    fn test_simulation_default_gravity() {
        let sim = HairSimulation::new();
        assert!((sim.gravity[1] + 9.81).abs() < 1e-9);
    }

    // 28. HairClustering: cluster count matches k
    #[test]
    fn test_clustering_cluster_count() {
        let mut hc = HairClustering::new();
        let roots: Vec<[f64; 2]> = (0..10).map(|i| [i as f64 * 0.1, 0.0]).collect();
        hc.cluster(&roots, 4);
        assert_eq!(hc.cluster_centres.len(), 4);
    }
}
