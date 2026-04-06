//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use crate::constraint::{DistanceConstraint, SoftConstraint};
use crate::particle::{SoftBody, SoftParticle};
use crate::solver::XpbdSolver;
use oxiphysics_core::math::{Real, Vec3};

/// A rope / hair chain simulated with Verlet integration and PBD constraints.
#[derive(Debug, Clone)]
pub struct Rope {
    /// Nodes that make up the chain.
    pub nodes: Vec<RopeNode>,
    /// Rest length of each segment between adjacent nodes.
    pub segment_length: f64,
    /// Stretch stiffness in \[0, 1\]; higher values resist elongation more.
    pub stiffness: f64,
    /// Bending stiffness in \[0, 1\]; higher values resist bending more.
    pub bending_stiffness: f64,
    /// Velocity damping coefficient.
    pub damping: f64,
    /// Visual radius of the rope (not used in physics).
    pub radius: f64,
}
impl Rope {
    /// Create a straight rope with `n_nodes` nodes starting at `start` and
    /// extending in `direction` (which will be normalised) for `total_length`.
    ///
    /// All nodes begin with mass `mass_per_node`; none are fixed initially.
    pub fn new_straight(
        n_nodes: usize,
        start: [f64; 3],
        direction: [f64; 3],
        total_length: f64,
        mass_per_node: f64,
    ) -> Self {
        assert!(n_nodes >= 2, "A rope needs at least 2 nodes");
        let segment_length = total_length / (n_nodes - 1) as f64;
        let dir = normalize3(direction);
        let mut nodes = Vec::with_capacity(n_nodes);
        for i in 0..n_nodes {
            let pos = add3(start, scale3(dir, segment_length * i as f64));
            nodes.push(RopeNode::new(pos, mass_per_node));
        }
        Self {
            nodes,
            segment_length,
            stiffness: 1.0,
            bending_stiffness: 0.3,
            damping: 0.01,
            radius: 0.01,
        }
    }
    /// Fix one end of the rope.
    ///
    /// * `end == 0` fixes node 0 (the start).
    /// * `end == 1` fixes the last node.
    pub fn fix_end(&mut self, end: usize) {
        match end {
            0 => self.nodes[0].fixed = true,
            1 => {
                let last = self.nodes.len() - 1;
                self.nodes[last].fixed = true;
            }
            _ => panic!("end must be 0 (start) or 1 (end)"),
        }
    }
    /// Apply gravitational acceleration to every free node by directly
    /// displacing `prev_position` so that the Verlet step picks it up.
    ///
    /// Internally we shift `prev_position` backward by `a * dt^2` which adds
    /// `a * dt^2` to the Verlet displacement `pos - prev_pos`.
    pub fn apply_gravity(&mut self, dt: f64, gravity: [f64; 3]) {
        let accel_disp = scale3(gravity, dt * dt);
        for node in &mut self.nodes {
            if node.fixed {
                continue;
            }
            node.prev_position = sub3(node.prev_position, accel_disp);
        }
    }
    /// Verlet integration step followed by constraint solving.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        self.apply_gravity(dt, gravity);
        let damping_factor = 1.0 - self.damping;
        for node in &mut self.nodes {
            if node.fixed {
                continue;
            }
            let vel = scale3(sub3(node.position, node.prev_position), 1.0 / dt);
            let vel_damped = scale3(vel, damping_factor);
            let new_pos = add3(node.position, scale3(vel_damped, dt));
            node.prev_position = node.position;
            node.position = new_pos;
        }
        self.solve_distance_constraints();
        self.solve_bending_constraints();
        for node in &mut self.nodes {
            if node.fixed {
                continue;
            }
            node.velocity = scale3(sub3(node.position, node.prev_position), 1.0 / dt);
        }
    }
    /// Project distance (stretch) constraints using PBD.
    pub fn solve_distance_constraints(&mut self) {
        let n = self.nodes.len();
        for i in 0..n - 1 {
            let pi = self.nodes[i].position;
            let pj = self.nodes[i + 1].position;
            let delta = sub3(pj, pi);
            let dist = len3(delta);
            if dist < 1e-15 {
                continue;
            }
            let error = dist - self.segment_length;
            let correction = scale3(normalize3(delta), error * self.stiffness * 0.5);
            let wi = if self.nodes[i].fixed { 0.0 } else { 1.0 };
            let wj = if self.nodes[i + 1].fixed { 0.0 } else { 1.0 };
            let total_w = wi + wj;
            if total_w < 1e-15 {
                continue;
            }
            if !self.nodes[i].fixed {
                self.nodes[i].position =
                    add3(self.nodes[i].position, scale3(correction, wi / total_w));
            }
            if !self.nodes[i + 1].fixed {
                self.nodes[i + 1].position =
                    sub3(self.nodes[i + 1].position, scale3(correction, wj / total_w));
            }
        }
    }
    /// Project bending constraints using consecutive triplets (i-1, i, i+1).
    pub fn solve_bending_constraints(&mut self) {
        let n = self.nodes.len();
        if n < 3 {
            return;
        }
        for i in 1..n - 1 {
            let pa = self.nodes[i - 1].position;
            let pb = self.nodes[i].position;
            let pc = self.nodes[i + 1].position;
            let mid = scale3(add3(pa, pc), 0.5);
            let delta = sub3(pb, mid);
            let correction = scale3(delta, self.bending_stiffness * 0.5);
            if !self.nodes[i - 1].fixed {
                self.nodes[i - 1].position =
                    add3(self.nodes[i - 1].position, scale3(correction, 0.25));
            }
            if !self.nodes[i].fixed {
                self.nodes[i].position = sub3(self.nodes[i].position, scale3(correction, 0.5));
            }
            if !self.nodes[i + 1].fixed {
                self.nodes[i + 1].position =
                    add3(self.nodes[i + 1].position, scale3(correction, 0.25));
            }
        }
    }
    /// Sum of inter-node distances (current total length of the rope).
    pub fn total_length(&self) -> f64 {
        self.nodes
            .windows(2)
            .map(|w| len3(sub3(w[1].position, w[0].position)))
            .sum()
    }
    /// Discrete curvature at node `i` (must satisfy `1 <= i <= n-2`).
    ///
    /// Computed as `|t2 - t1| / segment_length` where `t1` and `t2` are the
    /// unit tangents of the two adjacent segments.
    ///
    /// Returns 0.0 for boundary nodes or degenerate segments.
    pub fn curvature_at(&self, i: usize) -> f64 {
        let n = self.nodes.len();
        if i == 0 || i >= n - 1 {
            return 0.0;
        }
        let t1 = normalize3(sub3(self.nodes[i].position, self.nodes[i - 1].position));
        let t2 = normalize3(sub3(self.nodes[i + 1].position, self.nodes[i].position));
        let dt = sub3(t2, t1);
        if self.segment_length < 1e-15 {
            0.0
        } else {
            len3(dt) / self.segment_length
        }
    }
}
impl Rope {
    /// Resolve collisions between rope nodes and a sphere centred at `center`
    /// with `radius`. Nodes inside the sphere are pushed to the surface.
    ///
    /// Returns the list of collisions that were resolved.
    #[allow(dead_code)]
    pub fn collide_sphere(&mut self, center: [f64; 3], radius: f64) -> Vec<RopeCollision> {
        let mut collisions = Vec::new();
        for i in 0..self.nodes.len() {
            if self.nodes[i].fixed {
                continue;
            }
            let diff = sub3(self.nodes[i].position, center);
            let dist = len3(diff);
            if dist < radius && dist > 1e-15 {
                let normal = normalize3(diff);
                let penetration = radius - dist;
                self.nodes[i].position = add3(center, scale3(normal, radius));
                collisions.push(RopeCollision {
                    node_index: i,
                    penetration,
                    normal,
                });
            }
        }
        collisions
    }
    /// Resolve collisions between rope nodes and an infinite plane defined by
    /// `plane_point` and `plane_normal` (must be unit length).
    ///
    /// Nodes below the plane are pushed onto it.
    #[allow(dead_code)]
    pub fn collide_plane(
        &mut self,
        plane_point: [f64; 3],
        plane_normal: [f64; 3],
    ) -> Vec<RopeCollision> {
        let mut collisions = Vec::new();
        for i in 0..self.nodes.len() {
            if self.nodes[i].fixed {
                continue;
            }
            let to_node = sub3(self.nodes[i].position, plane_point);
            let dist = dot3(to_node, plane_normal);
            if dist < 0.0 {
                self.nodes[i].position = sub3(self.nodes[i].position, scale3(plane_normal, dist));
                collisions.push(RopeCollision {
                    node_index: i,
                    penetration: -dist,
                    normal: plane_normal,
                });
            }
        }
        collisions
    }
    /// Self-collision detection and response.
    ///
    /// Checks all non-adjacent node pairs separated by more than `skip_adjacent`
    /// segments. If two nodes are closer than `min_distance`, they are pushed
    /// apart symmetrically.
    #[allow(dead_code)]
    pub fn resolve_self_collision(&mut self, min_distance: f64, skip_adjacent: usize) {
        let n = self.nodes.len();
        for i in 0..n {
            if self.nodes[i].fixed {
                continue;
            }
            for j in (i + skip_adjacent + 1)..n {
                if self.nodes[j].fixed {
                    continue;
                }
                let diff = sub3(self.nodes[j].position, self.nodes[i].position);
                let dist = len3(diff);
                if dist < min_distance && dist > 1e-15 {
                    let normal = normalize3(diff);
                    let overlap = min_distance - dist;
                    let half = overlap * 0.5;
                    self.nodes[i].position = sub3(self.nodes[i].position, scale3(normal, half));
                    self.nodes[j].position = add3(self.nodes[j].position, scale3(normal, half));
                }
            }
        }
    }
    /// Rope–surface interaction: slide rope along a cylinder surface.
    ///
    /// The cylinder axis runs from `axis_start` to `axis_end` with `radius`.
    /// Nodes inside the cylinder are projected onto its surface.
    #[allow(dead_code)]
    #[allow(clippy::too_many_arguments)]
    pub fn collide_cylinder(
        &mut self,
        axis_start: [f64; 3],
        axis_end: [f64; 3],
        radius: f64,
    ) -> Vec<RopeCollision> {
        let axis = sub3(axis_end, axis_start);
        let axis_len = len3(axis);
        if axis_len < 1e-15 {
            return Vec::new();
        }
        let axis_dir = normalize3(axis);
        let mut collisions = Vec::new();
        for i in 0..self.nodes.len() {
            if self.nodes[i].fixed {
                continue;
            }
            let to_node = sub3(self.nodes[i].position, axis_start);
            let proj = dot3(to_node, axis_dir);
            if proj < 0.0 || proj > axis_len {
                continue;
            }
            let closest_on_axis = add3(axis_start, scale3(axis_dir, proj));
            let radial = sub3(self.nodes[i].position, closest_on_axis);
            let radial_dist = len3(radial);
            if radial_dist < radius && radial_dist > 1e-15 {
                let normal = normalize3(radial);
                let penetration = radius - radial_dist;
                self.nodes[i].position = add3(closest_on_axis, scale3(normal, radius));
                collisions.push(RopeCollision {
                    node_index: i,
                    penetration,
                    normal,
                });
            }
        }
        collisions
    }
    /// Compute the catenary shape for a rope of length `rope_length` suspended
    /// between two endpoints under gravity.
    ///
    /// Returns the positions of `n_nodes` evenly spaced nodes along the
    /// catenary curve in the XY plane (gravity in -Y). The first and last
    /// positions are `start` and `end` respectively.
    ///
    /// This is an approximate catenary using iterative Newton's method on the
    /// catenary parameter `a`.
    #[allow(dead_code)]
    pub fn catenary_positions(
        start: [f64; 3],
        end: [f64; 3],
        rope_length: f64,
        n_nodes: usize,
    ) -> Vec<[f64; 3]> {
        assert!(n_nodes >= 2, "Need at least 2 nodes");
        let dx = end[0] - start[0];
        let dy = end[1] - start[1];
        let dz = end[2] - start[2];
        let horiz_dist = (dx * dx + dz * dz).sqrt();
        let vert_diff = dy;
        let span = (dx * dx + dy * dy + dz * dz).sqrt();
        if rope_length <= span + 1e-12 {
            let mut positions = Vec::with_capacity(n_nodes);
            for i in 0..n_nodes {
                let t = i as f64 / (n_nodes - 1) as f64;
                positions.push([start[0] + dx * t, start[1] + dy * t, start[2] + dz * t]);
            }
            return positions;
        }
        let mut a = 1.0_f64;
        for _ in 0..50 {
            let rhs = (rope_length * rope_length - vert_diff * vert_diff).sqrt();
            let arg = horiz_dist / (2.0 * a);
            let f = 2.0 * a * arg.sinh() - rhs;
            let df = 2.0 * arg.sinh() - horiz_dist * arg.cosh() / a;
            if df.abs() < 1e-15 {
                break;
            }
            a -= f / df;
            if a < 1e-6 {
                a = 1e-6;
            }
        }
        let x0 = horiz_dist * 0.5 - a * ((vert_diff / rope_length).atanh()).clamp(-10.0, 10.0);
        let y0 = start[1] - a * ((0.0 - x0) / a).cosh();
        let horiz_dir = if horiz_dist > 1e-15 {
            [dx / horiz_dist, 0.0, dz / horiz_dist]
        } else {
            [1.0, 0.0, 0.0]
        };
        let mut positions = Vec::with_capacity(n_nodes);
        for i in 0..n_nodes {
            let t = i as f64 / (n_nodes - 1) as f64;
            let x_local = horiz_dist * t;
            let y_cat = a * ((x_local - x0) / a).cosh() + y0;
            positions.push([
                start[0] + horiz_dir[0] * x_local,
                y_cat,
                start[2] + horiz_dir[2] * x_local,
            ]);
        }
        positions
    }
    /// Detect whether the rope is taut (stretched to near its rest length).
    ///
    /// Returns `true` if `total_length() >= rest_length * (1.0 - tolerance)`.
    #[allow(dead_code)]
    pub fn is_taut(&self, tolerance: f64) -> bool {
        let rest = self.segment_length * (self.nodes.len() - 1) as f64;
        self.total_length() >= rest * (1.0 - tolerance)
    }
    /// Apply a constraint that limits the maximum stretch per segment.
    ///
    /// If any segment exceeds `max_stretch_ratio * segment_length`, the nodes
    /// are pulled back. This is a refinement pass that can be called after
    /// the standard distance constraints.
    #[allow(dead_code)]
    pub fn enforce_max_stretch(&mut self, max_stretch_ratio: f64) {
        let max_len = self.segment_length * max_stretch_ratio;
        let n = self.nodes.len();
        let max_iters = n * 10;
        for _iter in 0..max_iters {
            let mut converged = true;
            for i in 0..n - 1 {
                let diff = sub3(self.nodes[i + 1].position, self.nodes[i].position);
                let dist = len3(diff);
                if dist > max_len && dist > 1e-15 {
                    converged = false;
                    let correction_dir = normalize3(diff);
                    let excess = dist - max_len;
                    let half = excess * 0.5;
                    if !self.nodes[i].fixed {
                        self.nodes[i].position =
                            add3(self.nodes[i].position, scale3(correction_dir, half));
                    }
                    if !self.nodes[i + 1].fixed {
                        self.nodes[i + 1].position =
                            sub3(self.nodes[i + 1].position, scale3(correction_dir, half));
                    }
                }
            }
            if converged {
                break;
            }
        }
    }
    /// Compute the tension at each segment as the force magnitude.
    ///
    /// Returns a vector of length `n_nodes - 1` with the estimated tension
    /// in each segment, based on the displacement from rest length and stiffness.
    #[allow(dead_code)]
    pub fn segment_tensions(&self) -> Vec<f64> {
        let n = self.nodes.len();
        let mut tensions = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let diff = sub3(self.nodes[i + 1].position, self.nodes[i].position);
            let dist = len3(diff);
            let stretch = (dist - self.segment_length).max(0.0);
            tensions.push(self.stiffness * stretch / self.segment_length.max(1e-15));
        }
        tensions
    }
    /// Maximum tension across all segments.
    #[allow(dead_code)]
    pub fn max_tension(&self) -> f64 {
        self.segment_tensions().into_iter().fold(0.0_f64, f64::max)
    }
    /// Compute kinetic energy of the rope: 0.5 * sum(m * |v|^2).
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> f64 {
        self.nodes
            .iter()
            .filter(|n| !n.fixed)
            .map(|n| 0.5 * n.mass * dot3(n.velocity, n.velocity))
            .sum()
    }
    /// Compute gravitational potential energy relative to `ref_y`.
    #[allow(dead_code)]
    pub fn potential_energy(&self, gravity_y: f64, ref_y: f64) -> f64 {
        self.nodes
            .iter()
            .filter(|n| !n.fixed)
            .map(|n| n.mass * gravity_y.abs() * (n.position[1] - ref_y))
            .sum()
    }
    /// Compute the catenary tension distribution along the rope.
    ///
    /// For a rope hanging under gravity, the tension at each node is given by
    /// the catenary equation.  Requires the first node to be fixed (the anchor).
    ///
    /// `T(s) = w * a`  where `a` is the catenary parameter, `w` is the weight
    /// per unit length, and `s` is the arc-length from the lowest point.
    ///
    /// In this discrete approximation, the horizontal tension component
    /// `H = w * a` is assumed constant and computed from the total weight
    /// and the horizontal span.
    ///
    /// # Arguments
    /// * `gravity_y` – gravitational acceleration magnitude (m/s²).
    ///
    /// Returns a vector of tension magnitudes at each node (N).
    #[allow(dead_code)]
    pub fn compute_catenary_tension(&self, gravity_y: f64) -> Vec<f64> {
        let n = self.nodes.len();
        if n < 2 {
            return vec![0.0; n];
        }
        let m_avg: f64 = self.nodes.iter().map(|nd| nd.mass).sum::<f64>() / n as f64;
        let w = m_avg * gravity_y.abs() / self.segment_length.max(1e-15);
        let x_start = self.nodes[0].position[0];
        let x_end = self.nodes[n - 1].position[0];
        let span = (x_end - x_start).abs();
        let y_start = self.nodes[0].position[1];
        let y_end = self.nodes[n - 1].position[1];
        let y_min = self
            .nodes
            .iter()
            .map(|nd| nd.position[1])
            .fold(f64::INFINITY, f64::min);
        let sag = ((y_start + y_end) / 2.0 - y_min).max(1e-6);
        let a = span * span / (8.0 * sag).max(1e-15);
        let h = w * a;
        let x_mid = (x_start + x_end) / 2.0;
        self.nodes
            .iter()
            .map(|nd| {
                let s = (nd.position[0] - x_mid).abs();
                (h * h + (w * s) * (w * s)).sqrt()
            })
            .collect()
    }
    /// Compute the dynamic bending stiffness of the rope.
    ///
    /// The effective dynamic bending stiffness accounts for both material
    /// bending resistance and the tension-induced geometric stiffness:
    /// ```text
    /// EI_eff = EI_material + T * L²
    /// ```
    /// where `T` is the average tension, `L` is the segment length, and
    /// `EI_material = E * I` is the cross-section bending rigidity.
    ///
    /// # Arguments
    /// * `youngs_modulus` – E (Pa).
    /// * `radius`         – cross-section radius (m) for a circular rod.
    /// * `gravity_y`      – gravitational acceleration magnitude (m/s²).
    ///
    /// Returns the effective bending stiffness (N·m²).
    #[allow(dead_code)]
    pub fn compute_dynamic_stiffness(
        &self,
        youngs_modulus: f64,
        radius: f64,
        gravity_y: f64,
    ) -> f64 {
        let area_moment = std::f64::consts::PI * radius.powi(4) / 4.0;
        let ei_material = youngs_modulus * area_moment;
        let tensions = self.compute_catenary_tension(gravity_y);
        let avg_tension = if tensions.is_empty() {
            0.0
        } else {
            tensions.iter().sum::<f64>() / tensions.len() as f64
        };
        let l = self.segment_length;
        ei_material + avg_tension * l * l
    }
    /// Apply a twist constraint along the rope to limit relative twist between segments.
    ///
    /// Each pair of adjacent nodes has an associated twist angle (the rotation about
    /// the segment tangent).  This constraint projects the positions to bring the
    /// twist angle within `[-max_twist, max_twist]`.
    ///
    /// In the discrete approximation, twist is estimated from the rotation of an
    /// arbitrary perpendicular frame.  When the estimated twist exceeds the limit,
    /// a corrective displacement is applied orthogonal to the tangent.
    ///
    /// # Arguments
    /// * `max_twist` – maximum allowed twist per segment (radians).
    /// * `compliance` – XPBD compliance for the twist correction (m²/N).  0 = rigid.
    /// * `dt`         – current time step (s).
    #[allow(dead_code)]
    pub fn apply_twist_constraint(&mut self, max_twist: f64, compliance: f64, dt: f64) {
        let n = self.nodes.len();
        if n < 3 {
            return;
        }
        for i in 0..n - 2 {
            let t0 = normalize3(sub3(self.nodes[i + 1].position, self.nodes[i].position));
            let t1 = normalize3(sub3(self.nodes[i + 2].position, self.nodes[i + 1].position));
            let cross = [
                t0[1] * t1[2] - t0[2] * t1[1],
                t0[2] * t1[0] - t0[0] * t1[2],
                t0[0] * t1[1] - t0[1] * t1[0],
            ];
            let sin_angle = len3(cross);
            let cos_angle = dot3(t0, t1).clamp(-1.0, 1.0);
            let angle = sin_angle.atan2(cos_angle);
            if angle.abs() <= max_twist {
                continue;
            }
            let excess = angle - angle.signum() * max_twist;
            let w_sum = if !self.nodes[i + 1].fixed && !self.nodes[i + 2].fixed {
                2.0
            } else if !self.nodes[i + 1].fixed || !self.nodes[i + 2].fixed {
                1.0
            } else {
                continue;
            };
            let alpha_tilde = compliance / (dt * dt).max(1e-30);
            let delta = excess / (w_sum + alpha_tilde);
            let perp = normalize3(cross);
            let correction = scale3(perp, delta * 0.5);
            if !self.nodes[i + 1].fixed {
                self.nodes[i + 1].position = add3(self.nodes[i + 1].position, correction);
            }
            if !self.nodes[i + 2].fixed {
                self.nodes[i + 2].position = sub3(self.nodes[i + 2].position, correction);
            }
        }
    }
}
/// A single elastic segment connecting two particles.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RopeSegment {
    /// Position of the first endpoint.
    pub p0: [f64; 3],
    /// Position of the second endpoint.
    pub p1: [f64; 3],
    /// Natural (rest) length of the segment.
    pub rest_length: f64,
    /// Mass of the segment (equally split to each endpoint for accumulation).
    pub mass: f64,
    /// Stretch stiffness (spring constant, N/m).
    pub k_stretch: f64,
}
#[allow(dead_code)]
impl RopeSegment {
    /// Create a new segment.
    pub fn new(p0: [f64; 3], p1: [f64; 3], rest_length: f64, mass: f64, k_stretch: f64) -> Self {
        Self {
            p0,
            p1,
            rest_length,
            mass,
            k_stretch,
        }
    }
    /// Spring force on p0 (first return value) and p1 (second return value).
    ///
    /// Uses Hooke's law: F = k * (d - rest) * dir
    /// where dir points from p0 to p1.
    pub fn elastic_force(&self) -> ([f64; 3], [f64; 3]) {
        let diff = sub3(self.p1, self.p0);
        let d = len3(diff);
        if d < 1e-15 {
            return ([0.0; 3], [0.0; 3]);
        }
        let dir = normalize3(diff);
        let stretch = d - self.rest_length;
        let f_mag = self.k_stretch * stretch;
        let force = scale3(dir, f_mag);
        (force, scale3(force, -1.0))
    }
}
/// Kelvin–Voigt viscoelastic segment connecting two particles.
///
/// The constitutive model is `σ = E ε + η dε/dt` where:
/// - `E` is the elastic (Young's) modulus (spring constant).
/// - `η` is the viscosity coefficient (damper constant).
/// - `ε` is the engineering strain `(l - l₀) / l₀`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ViscoelasticSegment {
    /// Index of the first node.
    pub node_a: usize,
    /// Index of the second node.
    pub node_b: usize,
    /// Natural (rest) length l₀ (m).
    pub rest_length: f64,
    /// Elastic modulus times cross-section area E·A (N).
    pub ea: f64,
    /// Viscous damping coefficient η·A (N·s).
    pub eta_a: f64,
    /// Previous strain (for finite-difference strain rate).
    pub prev_strain: f64,
}
#[allow(dead_code)]
impl ViscoelasticSegment {
    /// Create a new segment between nodes `a` and `b`.
    pub fn new(node_a: usize, node_b: usize, rest_length: f64, ea: f64, eta_a: f64) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            ea,
            eta_a,
            prev_strain: 0.0,
        }
    }
    /// Compute the current engineering strain.
    pub fn strain(&self, positions: &[[f64; 3]]) -> f64 {
        let pa = positions[self.node_a];
        let pb = positions[self.node_b];
        let l = len3(sub3(pb, pa));
        (l - self.rest_length) / self.rest_length.max(1e-15)
    }
    /// Compute the Kelvin–Voigt force on node_a (and equal-and-opposite on node_b).
    ///
    /// `dt` is the time step (used to estimate strain rate).
    /// Returns `(force_on_a, force_on_b)`.
    pub fn kelvin_voigt_force(&mut self, positions: &[[f64; 3]], dt: f64) -> ([f64; 3], [f64; 3]) {
        let pa = positions[self.node_a];
        let pb = positions[self.node_b];
        let diff = sub3(pb, pa);
        let l = len3(diff);
        if l < 1e-15 {
            return ([0.0; 3], [0.0; 3]);
        }
        let dir = normalize3(diff);
        let eps = (l - self.rest_length) / self.rest_length.max(1e-15);
        let eps_dot = (eps - self.prev_strain) / dt.max(1e-15);
        self.prev_strain = eps;
        let sigma = self.ea * eps + self.eta_a * eps_dot;
        let fa = scale3(dir, sigma);
        (fa, scale3(fa, -1.0))
    }
    /// Elastic potential energy stored in this segment.
    pub fn elastic_energy(&self, positions: &[[f64; 3]]) -> f64 {
        let eps = self.strain(positions);
        0.5 * self.ea * eps * eps * self.rest_length
    }
}
/// A single node in a [`Rope`] chain.
#[derive(Debug, Clone)]
pub struct RopeNode {
    /// Current world-space position.
    pub position: [f64; 3],
    /// Position at the previous time step (used by Verlet integration).
    pub prev_position: [f64; 3],
    /// Current velocity (derived from Verlet; also used to apply impulses).
    pub velocity: [f64; 3],
    /// Mass of the node in kg.
    pub mass: f64,
    /// When `true` the node is kinematically fixed and will not move.
    pub fixed: bool,
}
impl RopeNode {
    fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            prev_position: position,
            velocity: [0.0; 3],
            mass,
            fixed: false,
        }
    }
}
/// Result of a rope–sphere collision check.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RopeCollision {
    /// Index of the colliding node.
    pub node_index: usize,
    /// Penetration depth (positive = overlapping).
    pub penetration: f64,
    /// Surface normal at the contact point (points away from the sphere centre).
    pub normal: [f64; 3],
}
/// A rope built from [`ViscoelasticSegment`]s.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct ViscoelasticRope {
    /// Node positions.
    pub positions: Vec<[f64; 3]>,
    /// Node velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Node masses.
    pub masses: Vec<f64>,
    /// Fixed flags.
    pub fixed: Vec<bool>,
    /// Viscoelastic segments.
    pub segments: Vec<ViscoelasticSegment>,
}
#[allow(dead_code)]
impl ViscoelasticRope {
    /// Build a straight viscoelastic rope with `n` nodes.
    #[allow(clippy::too_many_arguments)]
    pub fn new_straight(
        n: usize,
        start: [f64; 3],
        direction: [f64; 3],
        seg_len: f64,
        mass_per_node: f64,
        ea: f64,
        eta_a: f64,
    ) -> Self {
        assert!(n >= 2);
        let dir = normalize3(direction);
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|i| add3(start, scale3(dir, seg_len * i as f64)))
            .collect();
        let velocities = vec![[0.0; 3]; n];
        let masses = vec![mass_per_node; n];
        let fixed = vec![false; n];
        let segments: Vec<ViscoelasticSegment> = (0..n - 1)
            .map(|i| ViscoelasticSegment::new(i, i + 1, seg_len, ea, eta_a))
            .collect();
        Self {
            positions,
            velocities,
            masses,
            fixed,
            segments,
        }
    }
    /// One semi-implicit Euler step under gravity.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let positions_snap = self.positions.clone();
        for seg in &mut self.segments {
            let (fa, fb) = seg.kelvin_voigt_force(&positions_snap, dt);
            forces[seg.node_a] = add3(forces[seg.node_a], fa);
            forces[seg.node_b] = add3(forces[seg.node_b], fb);
        }
        for i in 0..n {
            if self.fixed[i] {
                continue;
            }
            let inv_m = 1.0 / self.masses[i].max(1e-30);
            let a = [
                forces[i][0] * inv_m + gravity[0],
                forces[i][1] * inv_m + gravity[1],
                forces[i][2] * inv_m + gravity[2],
            ];
            self.velocities[i] = add3(self.velocities[i], scale3(a, dt));
            self.positions[i] = add3(self.positions[i], scale3(self.velocities[i], dt));
        }
    }
    /// Total elastic energy of all segments.
    pub fn total_elastic_energy(&self) -> f64 {
        self.segments
            .iter()
            .map(|s| s.elastic_energy(&self.positions))
            .sum()
    }
    /// Total kinetic energy.
    pub fn total_kinetic_energy(&self) -> f64 {
        self.positions
            .iter()
            .enumerate()
            .filter(|(i, _)| !self.fixed[*i])
            .map(|(i, _)| 0.5 * self.masses[i] * dot3(self.velocities[i], self.velocities[i]))
            .sum()
    }
}
/// A rope represented as a chain of particles connected by XPBD distance
/// constraints.
///
/// This is the original implementation.  For Verlet-based hair simulation
/// see [`Rope`] and [`HairSystem`].
#[derive(Debug, Clone)]
pub struct XpbdRope {
    /// The underlying soft body.
    pub body: SoftBody,
    /// Distance constraints along the chain.
    pub constraints: Vec<DistanceConstraint>,
}
impl XpbdRope {
    /// Create a straight rope between `start` and `end`.
    ///
    /// * `num_segments` – number of segments (one fewer than particles).
    /// * `mass_per_segment` – mass of each particle (except the first, which
    ///   is fixed by default).
    /// * `compliance` – XPBD compliance for distance constraints.
    pub fn new(
        start: Vec3,
        end: Vec3,
        num_segments: usize,
        mass_per_segment: Real,
        compliance: Real,
    ) -> Self {
        assert!(num_segments >= 1, "Rope needs at least one segment");
        let num_particles = num_segments + 1;
        let mut particles = Vec::with_capacity(num_particles);
        for i in 0..num_particles {
            let t = i as Real / num_segments as Real;
            let pos = start + (end - start) * t;
            if i == 0 {
                particles.push(SoftParticle::new_static(pos));
            } else {
                particles.push(SoftParticle::new(pos, mass_per_segment));
            }
        }
        let body = SoftBody::from_particles(particles);
        let mut constraints = Vec::with_capacity(num_segments);
        for i in 0..num_segments {
            constraints.push(DistanceConstraint::from_particles(
                i,
                i + 1,
                &body.particles,
                compliance,
            ));
        }
        Self { body, constraints }
    }
    /// Number of particles in the rope.
    pub fn num_particles(&self) -> usize {
        self.body.particles.len()
    }
    /// Number of segments (distance constraints).
    pub fn num_segments(&self) -> usize {
        self.constraints.len()
    }
}
/// A single hair strand: a [`Rope`] with an associated display colour.
#[derive(Debug, Clone)]
pub struct HairStrand {
    /// The underlying rope simulation.
    pub rope: Rope,
    /// RGB colour in \[0, 1\] range.
    pub color: [f64; 3],
}
impl HairStrand {
    /// Create a hair strand hanging straight down from `root`.
    ///
    /// The strand has `n_nodes` nodes and a rest length of `length`.
    pub fn new(root: [f64; 3], length: f64, n_nodes: usize) -> Self {
        let direction = [0.0, -1.0, 0.0];
        let mut rope = Rope::new_straight(n_nodes, root, direction, length, 1.0 / n_nodes as f64);
        rope.fix_end(0);
        Self {
            rope,
            color: [0.6, 0.4, 0.2],
        }
    }
}
/// A Verlet-integrated rope with per-segment rest lengths.
///
/// This is the flat-array representation (positions / velocities / masses /
/// rest_lengths) requested for the new API.  The original [`Rope`] (node-based)
/// is kept for backward compatibility.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct RopeVerlet {
    /// Particle positions.
    pub particles: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses (kg).
    pub masses: Vec<f64>,
    /// Rest length of each segment (length = particles.len() - 1).
    pub rest_lengths: Vec<f64>,
    /// Whether each particle is fixed (pinned).
    pub fixed: Vec<bool>,
}
#[allow(dead_code)]
impl RopeVerlet {
    /// Create a straight rope with `n` particles, starting at `start` and
    /// extending in `direction` (normalised internally).
    ///
    /// All particles have mass `mass_per_segment` and none are fixed initially.
    pub fn new(
        n: usize,
        start: [f64; 3],
        direction: [f64; 3],
        segment_length: f64,
        mass_per_segment: f64,
    ) -> Self {
        assert!(n >= 2, "Rope needs at least 2 particles");
        let dir = normalize3(direction);
        let mut particles = Vec::with_capacity(n);
        for i in 0..n {
            particles.push(add3(start, scale3(dir, segment_length * i as f64)));
        }
        let velocities = vec![[0.0_f64; 3]; n];
        let masses = vec![mass_per_segment; n];
        let rest_lengths = vec![segment_length; n - 1];
        let fixed = vec![false; n];
        Self {
            particles,
            velocities,
            masses,
            rest_lengths,
            fixed,
        }
    }
    /// Total arc length (sum of current inter-particle distances).
    pub fn total_length(&self) -> f64 {
        self.particles
            .windows(2)
            .map(|w| len3(sub3(w[1], w[0])))
            .sum()
    }
    /// One Verlet + PBD step under gravity with velocity damping.
    ///
    /// `k_stretch` is a PBD stiffness in \[0, 1\];
    /// `damping` is a velocity damping coefficient per step.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3], k_stretch: f64, damping: f64) {
        let n = self.particles.len();
        for i in 0..n {
            if self.fixed[i] {
                continue;
            }
            self.velocities[i] = add3(self.velocities[i], scale3(gravity, dt));
            self.velocities[i] = scale3(self.velocities[i], 1.0 - damping);
            self.particles[i] = add3(self.particles[i], scale3(self.velocities[i], dt));
        }
        self.apply_constraint_iterations(k_stretch, 1);
    }
    /// Run `n_iters` PBD distance constraint iterations with stiffness `k`.
    pub fn apply_constraint_iterations(&mut self, k: f64, n_iters: usize) {
        let n = self.particles.len();
        for _ in 0..n_iters {
            for seg in 0..n - 1 {
                let pi = self.particles[seg];
                let pj = self.particles[seg + 1];
                let diff = sub3(pj, pi);
                let dist = len3(diff);
                if dist < 1e-15 {
                    continue;
                }
                let rest = self.rest_lengths[seg];
                let error = dist - rest;
                let correction = scale3(normalize3(diff), error * k * 0.5);
                let wi = if self.fixed[seg] { 0.0 } else { 1.0 };
                let wj = if self.fixed[seg + 1] { 0.0 } else { 1.0 };
                let total_w = wi + wj;
                if total_w < 1e-15 {
                    continue;
                }
                if !self.fixed[seg] {
                    self.particles[seg] = add3(pi, scale3(correction, wi / total_w));
                }
                if !self.fixed[seg + 1] {
                    self.particles[seg + 1] = sub3(pj, scale3(correction, wj / total_w));
                }
            }
        }
    }
    /// Compute the catenary SAG parameter `a` for the rope in its current shape.
    ///
    /// Fits to y = a * cosh(x/a) + c in the vertical plane.  Returns `a`
    /// (positive) or 0.0 if the rope is straight.
    ///
    /// This is an approximate Newton solve on the horizontal span and total length.
    pub fn catenary_shape(&self) -> Vec<f64> {
        let n = self.particles.len();
        if n < 2 {
            return vec![0.0];
        }
        let start = self.particles[0];
        let end = self.particles[n - 1];
        let dx = end[0] - start[0];
        let dy = end[1] - start[1];
        let dz = end[2] - start[2];
        let horiz_dist = (dx * dx + dz * dz).sqrt();
        let vert_diff = dy;
        let rope_length = self.total_length();
        let span = (dx * dx + dy * dy + dz * dz).sqrt();
        if rope_length <= span + 1e-12 || horiz_dist < 1e-12 {
            return vec![0.0];
        }
        let mut a = horiz_dist.max(1e-6);
        for _ in 0..50 {
            let rhs = (rope_length * rope_length - vert_diff * vert_diff).sqrt();
            let arg = horiz_dist / (2.0 * a);
            let f = 2.0 * a * arg.sinh() - rhs;
            let df = 2.0 * arg.sinh() - horiz_dist * arg.cosh() / a;
            if df.abs() < 1e-15 {
                break;
            }
            a -= f / df;
            if a < 1e-6 {
                a = 1e-6;
            }
        }
        vec![a]
    }
}
/// A collection of [`HairStrand`]s that are stepped together.
#[derive(Debug, Clone)]
pub struct HairSystem {
    /// All strands in the system.
    pub strands: Vec<HairStrand>,
}
impl HairSystem {
    /// Create an empty hair system.
    pub fn new() -> Self {
        Self {
            strands: Vec::new(),
        }
    }
    /// Add a strand to the system.
    pub fn add_strand(&mut self, strand: HairStrand) {
        self.strands.push(strand);
    }
    /// Step all strands forward by `dt` under `gravity`.
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        for strand in &mut self.strands {
            strand.rope.step(dt, gravity);
        }
    }
    /// Number of strands currently in the system.
    pub fn strand_count(&self) -> usize {
        self.strands.len()
    }
}
/// Material frame associated with a single rod segment.
///
/// The Cosserat (special) Kirchhoff rod model attaches a right-handed
/// orthonormal frame (d1, d2, d3) to every segment:
/// - `d3` is the tangent (unit) to the centreline.
/// - `d1`, `d2` are two directors that track twist around the tangent.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy)]
pub struct MaterialFrame {
    /// First director (spans the cross-section).
    pub d1: [f64; 3],
    /// Second director (spans the cross-section, perpendicular to d1).
    pub d2: [f64; 3],
    /// Tangent director (along the centreline).
    pub d3: [f64; 3],
}
#[allow(dead_code)]
impl MaterialFrame {
    /// Construct a frame from the segment tangent.
    ///
    /// `d3` is set to the normalised tangent; `d1` and `d2` are built
    /// with a stable Gram-Schmidt complement.
    pub fn from_tangent(tangent: [f64; 3]) -> Self {
        let d3 = normalize3(tangent);
        let d1 = Self::perpendicular_to(d3);
        let d2 = cross3(d3, d1);
        Self { d1, d2, d3 }
    }
    /// Rotate the frame by angle `phi` around `d3` (twist).
    pub fn twist(&self, phi: f64) -> Self {
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let d1_new = [
            cos_phi * self.d1[0] + sin_phi * self.d2[0],
            cos_phi * self.d1[1] + sin_phi * self.d2[1],
            cos_phi * self.d1[2] + sin_phi * self.d2[2],
        ];
        let d2_new = cross3(self.d3, d1_new);
        Self {
            d1: d1_new,
            d2: d2_new,
            d3: self.d3,
        }
    }
    /// Build a stable unit vector perpendicular to `n`.
    fn perpendicular_to(n: [f64; 3]) -> [f64; 3] {
        let ax: [f64; 3] = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
            [1.0, 0.0, 0.0]
        } else if n[1].abs() <= n[2].abs() {
            [0.0, 1.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        let d = dot3(ax, n);
        normalize3([ax[0] - d * n[0], ax[1] - d * n[1], ax[2] - d * n[2]])
    }
}
/// Discrete Kirchhoff (Cosserat) rod with per-segment material frames.
///
/// The rod stores a sequence of positions and one material frame per *segment*
/// (one fewer than nodes).  Bending and twisting energies are computed from
/// the discrete Darboux vector at each interior node.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct KirchhoffRod {
    /// Node positions along the centreline.
    pub positions: Vec<[f64; 3]>,
    /// Material frame for each segment (length = positions.len() - 1).
    pub frames: Vec<MaterialFrame>,
    /// Rest segment length (uniform).
    pub rest_length: f64,
    /// Bending stiffness EI (N·m²).
    pub bending_stiffness: f64,
    /// Twisting stiffness GJ (N·m²).
    pub twisting_stiffness: f64,
}
#[allow(dead_code)]
impl KirchhoffRod {
    /// Build a straight rod with `n` nodes and uniform rest length `seg_len`.
    ///
    /// All frames are initialised from the tangent of the straight line.
    pub fn new_straight(
        n: usize,
        start: [f64; 3],
        direction: [f64; 3],
        seg_len: f64,
        bending_stiffness: f64,
        twisting_stiffness: f64,
    ) -> Self {
        assert!(n >= 2, "rod needs at least 2 nodes");
        let dir = normalize3(direction);
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|i| add3(start, scale3(dir, seg_len * i as f64)))
            .collect();
        let frames: Vec<MaterialFrame> = (0..n - 1)
            .map(|_| MaterialFrame::from_tangent(dir))
            .collect();
        Self {
            positions,
            frames,
            rest_length: seg_len,
            bending_stiffness,
            twisting_stiffness,
        }
    }
    /// Discrete Darboux (curvature-twist) vector at interior node `i`.
    ///
    /// `i` must satisfy `1 <= i <= n-2`.  Returns `[kappa1, kappa2, tau]`
    /// in the material frame of segment `i-1`.
    pub fn darboux_vector(&self, i: usize) -> [f64; 3] {
        let n = self.positions.len();
        if i == 0 || i >= n - 1 {
            return [0.0; 3];
        }
        let f_prev = &self.frames[i - 1];
        let f_next = &self.frames[i];
        let inv_l = 1.0 / self.rest_length.max(1e-15);
        let kappa1 = dot3(f_prev.d2, f_next.d3) * inv_l;
        let kappa2 = dot3(f_prev.d3, f_next.d1) * inv_l;
        let tau = dot3(f_prev.d1, f_next.d2) * inv_l;
        [kappa1, kappa2, tau]
    }
    /// Bending energy for the entire rod (sum over interior nodes).
    ///
    /// `E_bend = 0.5 * EI * Σ (κ₁² + κ₂²) * l`
    pub fn bending_energy(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 1..n - 1 {
            let dv = self.darboux_vector(i);
            e += dv[0] * dv[0] + dv[1] * dv[1];
        }
        0.5 * self.bending_stiffness * e * self.rest_length
    }
    /// Twisting energy for the entire rod.
    ///
    /// `E_twist = 0.5 * GJ * Σ τ² * l`
    pub fn twisting_energy(&self) -> f64 {
        let n = self.positions.len();
        let mut e = 0.0;
        for i in 1..n - 1 {
            let dv = self.darboux_vector(i);
            e += dv[2] * dv[2];
        }
        0.5 * self.twisting_stiffness * e * self.rest_length
    }
    /// Total elastic energy (bending + twisting).
    pub fn total_elastic_energy(&self) -> f64 {
        self.bending_energy() + self.twisting_energy()
    }
    /// Update material frames from current node positions.
    ///
    /// The tangent of each segment is recomputed; `d1` and `d2` are
    /// transported parallel from the previous frame (approximate parallel
    /// transport by projection and Gram-Schmidt re-orthogonalisation).
    pub fn update_frames(&mut self) {
        let n = self.positions.len();
        for s in 0..n - 1 {
            let tangent = normalize3(sub3(self.positions[s + 1], self.positions[s]));
            let d1_old = self.frames[s].d1;
            let proj = dot3(d1_old, tangent);
            let d1_proj = sub3(d1_old, scale3(tangent, proj));
            let d1_new = if len3(d1_proj) > 1e-14 {
                normalize3(d1_proj)
            } else {
                MaterialFrame::perpendicular_to(tangent)
            };
            let d2_new = cross3(tangent, d1_new);
            self.frames[s] = MaterialFrame {
                d1: d1_new,
                d2: d2_new,
                d3: tangent,
            };
        }
    }
}
/// Convenience solver wrapper for [`XpbdRope`] simulation.
#[derive(Debug)]
pub struct RopeSolver {
    /// Inner XPBD solver.
    pub solver: XpbdSolver,
}
impl RopeSolver {
    /// Create a rope solver with the given number of substeps.
    pub fn new(num_substeps: usize) -> Self {
        Self {
            solver: XpbdSolver::new(num_substeps),
        }
    }
    /// Step the rope simulation forward by `dt`.
    pub fn step(&mut self, rope: &mut XpbdRope, dt: Real) {
        let mut boxed: Vec<Box<dyn SoftConstraint>> = rope
            .constraints
            .iter()
            .cloned()
            .map(|c| Box::new(c) as Box<dyn SoftConstraint>)
            .collect();
        self.solver.solve(&mut rope.body, &mut boxed, dt);
    }
}
