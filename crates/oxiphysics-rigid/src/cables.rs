#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cable and catenary dynamics.
//!
//! Provides analytical catenary models, lumped-mass cable simulation,
//! mooring line dynamics, and utility functions for cable engineering.
//!
//! # Overview
//!
//! - [`CatenarySegment`] — analytical catenary curve y = a·cosh((x−x₀)/a) + y₀
//! - [`CableNode`] — single node in a lumped-mass cable model
//! - [`CableSimulation`] — spring-damper cable simulation with gravity
//! - [`MooringLine`] — seabed-contact mooring cable for offshore applications
//! - [`wire_sag_midpoint`] — quick parabolic sag formula
//! - [`cable_eigenfrequency`] — natural frequencies of taut cables

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Physical constants
// ---------------------------------------------------------------------------

/// Standard gravitational acceleration in m/s².
const G: f64 = 9.80665;

// ---------------------------------------------------------------------------
// Helper – 3-D vector arithmetic on plain arrays
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
fn vec3_len(a: [f64; 3]) -> f64 {
    vec3_dot(a, a).sqrt()
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> [f64; 3] {
    let len = vec3_len(a);
    if len < 1e-300 {
        [0.0; 3]
    } else {
        vec3_scale(a, 1.0 / len)
    }
}

// ---------------------------------------------------------------------------
// CatenarySegment
// ---------------------------------------------------------------------------

/// Analytical catenary segment defined by y = a·cosh((x − x₀)/a) + y₀.
///
/// The parameter `a` (catenary constant) equals the ratio of horizontal
/// tension to weight per unit length.
pub struct CatenarySegment {
    /// Catenary constant a = H / (ρ_l · g), where H is the horizontal tension.
    pub a: f64,
    /// Horizontal offset x₀ of the catenary vertex.
    pub x_offset: f64,
    /// Vertical offset y₀ (height of the catenary vertex).
    pub y_offset: f64,
}

impl CatenarySegment {
    /// Construct a catenary segment from two support endpoints and the cable
    /// length between them.
    ///
    /// The catenary constant `a` is found by Newton-Raphson iteration on the
    /// arc-length equation:
    ///
    /// ```text
    /// L = a * (sinh(x2'/a) − sinh(x1'/a))
    /// ```
    ///
    /// where x1' = x1 − x₀ and x2' = x2 − x₀.
    ///
    /// Returns a default (a = 1, zero offsets) when no finite solution is
    /// found within the iteration limit.
    pub fn from_endpoints(x1: f64, y1: f64, x2: f64, y2: f64, length: f64) -> Self {
        let span = x2 - x1;
        let rise = y2 - y1;
        let chord = (span * span + rise * rise).sqrt();

        // Degenerate: endpoints coincide or cable shorter than chord
        if chord < 1e-12 || length <= chord {
            return Self {
                a: 1.0,
                x_offset: 0.0,
                y_offset: 0.0,
            };
        }

        // Newton-Raphson iteration to solve for a
        // Residual: f(a) = 2a*sinh(span/(2a)) - L
        // where we work in the symmetric half-span frame first
        let half_span = span * 0.5;
        let mut a = chord * 0.5; // initial guess
        for _ in 0..100 {
            let ratio = half_span / a;
            let sinh_r = ratio.sinh();
            let cosh_r = ratio.cosh();
            let f = 2.0 * a * sinh_r - length;
            let df = 2.0 * (sinh_r - ratio * cosh_r);
            if df.abs() < 1e-14 {
                break;
            }
            let da = -f / df;
            a += da;
            if a < 1e-9 {
                a = 1e-9;
            }
            if da.abs() < 1e-10 {
                break;
            }
        }

        // x₀ is such that the vertex is equidistant in arc-length from both ends.
        // For level ends the vertex x₀ = (x1 + x2) / 2.
        // For unequal heights we use: x₀ = (x1 + x2)/2 − (a/2)*ln((y2−y1+L_half)/(y2−y1−L_half+2*L))
        // A simpler closed-form from hyperbolic inversion:
        let x0 = (x1 + x2) * 0.5
            - a * ((rise) / (length)).asinh() * 0.0  // correction term
            // correct offset for asymmetric case
            + a * (rise / (2.0 * a)).asinh();

        // y₀ follows from y1 = a*cosh((x1 - x0)/a) + y0
        let y0 = y1 - a * ((x1 - x0) / a).cosh();

        Self {
            a,
            x_offset: x0,
            y_offset: y0,
        }
    }

    /// Height of the catenary at horizontal position `x`.
    ///
    /// y = a · cosh((x − x₀) / a) + y₀
    pub fn height_at(&self, x: f64) -> f64 {
        self.a * ((x - self.x_offset) / self.a).cosh() + self.y_offset
    }

    /// Tension magnitude at horizontal position `x`.
    ///
    /// The tension along the cable is T = ρ_l · g · a · cosh((x − x₀)/a),
    /// where ρ_l is the linear mass density (kg/m).
    ///
    /// # Parameters
    /// - `x`     – horizontal position
    /// - `rho_l` – linear mass density of the cable (kg/m)
    pub fn tension_at(&self, x: f64, rho_l: f64) -> f64 {
        rho_l * G * self.a * ((x - self.x_offset) / self.a).cosh()
    }

    /// Maximum sag (downward displacement from the chord) over the span
    /// \[x1, x2\].
    ///
    /// Evaluated by finding the vertex (lowest point) and comparing to the
    /// linear chord interpolation.
    pub fn max_sag(&self, x1: f64, x2: f64) -> f64 {
        // Sample the catenary at many points and subtract the chord
        let n = 200;
        let mut max_sag: f64 = 0.0;
        let y_chord_at = |x: f64| {
            let t = (x - x1) / (x2 - x1).max(1e-12);
            self.height_at(x1) * (1.0 - t) + self.height_at(x2) * t
        };
        for i in 0..=n {
            let x = x1 + (x2 - x1) * i as f64 / n as f64;
            let sag = y_chord_at(x) - self.height_at(x);
            if sag > max_sag {
                max_sag = sag;
            }
        }
        max_sag
    }

    /// Arc length of the catenary between x1 and x2.
    ///
    /// L = a · (sinh((x2 − x₀)/a) − sinh((x1 − x₀)/a))
    pub fn arc_length(&self, x1: f64, x2: f64) -> f64 {
        let s1 = ((x1 - self.x_offset) / self.a).sinh();
        let s2 = ((x2 - self.x_offset) / self.a).sinh();
        self.a * (s2 - s1)
    }
}

// ---------------------------------------------------------------------------
// CableNode
// ---------------------------------------------------------------------------

/// A single node in a lumped-mass cable model.
pub struct CableNode {
    /// Position in world space \[m\].
    pub pos: [f64; 3],
    /// Velocity in world space \[m/s\].
    pub vel: [f64; 3],
    /// Node mass \[kg\].
    pub mass: f64,
    /// If `true` the node is pinned (position and velocity are not updated).
    pub fixed: bool,
}

impl CableNode {
    /// Create a new free node at the given position.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            fixed: false,
        }
    }

    /// Create a fixed (pinned) node at the given position.
    pub fn new_fixed(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            vel: [0.0; 3],
            mass,
            fixed: true,
        }
    }
}

// ---------------------------------------------------------------------------
// CableSimulation
// ---------------------------------------------------------------------------

/// Lumped-mass cable simulation using spring-damper segments.
///
/// The cable is discretised into `n_nodes` point masses connected by
/// Hookean springs.  Each segment has rest length `L₀ = total_length /
/// (n_nodes − 1)` and spring constant `stiffness` (N/m).
pub struct CableSimulation {
    /// Node list.
    pub nodes: Vec<CableNode>,
    /// Rest lengths of segments between consecutive nodes \[m\].
    pub rest_lengths: Vec<f64>,
    /// Spring stiffness \[N/m\].
    pub stiffness: f64,
    /// Viscous damping coefficient \[N·s/m\].
    pub damping: f64,
}

impl CableSimulation {
    /// Create a new cable hanging vertically from the origin.
    ///
    /// Nodes are evenly spaced along the negative-Z axis.
    /// The first node is pinned (fixed).
    ///
    /// # Parameters
    /// - `n_nodes`          – number of nodes (≥ 2)
    /// - `total_length`     – undeformed cable length \[m\]
    /// - `mass_per_length`  – linear mass density \[kg/m\]
    /// - `stiffness`        – axial spring stiffness \[N/m\]
    /// - `damping`          – viscous damping coefficient \[N·s/m\]
    pub fn new(
        n_nodes: usize,
        total_length: f64,
        mass_per_length: f64,
        stiffness: f64,
        damping: f64,
    ) -> Self {
        assert!(n_nodes >= 2, "cable must have at least 2 nodes");
        let seg_len = total_length / (n_nodes - 1) as f64;
        let node_mass = mass_per_length * total_length / n_nodes as f64;

        let mut nodes = Vec::with_capacity(n_nodes);
        for i in 0..n_nodes {
            let z = -(i as f64) * seg_len;
            let mut node = CableNode::new([0.0, 0.0, z], node_mass);
            if i == 0 {
                node.fixed = true;
            }
            nodes.push(node);
        }

        let rest_lengths = vec![seg_len; n_nodes - 1];

        Self {
            nodes,
            rest_lengths,
            stiffness,
            damping,
        }
    }

    /// Compute the net elastic (spring + damper) force on node `i` from its
    /// two adjacent segments (if they exist).
    ///
    /// Returns `[fx, fy, fz]` in world space.
    pub fn elastic_force(&self, i: usize) -> [f64; 3] {
        let mut force = [0.0_f64; 3];
        let n = self.nodes.len();

        // Contribution from segment (i-1, i)
        if i > 0 {
            let j = i - 1;
            let f = self.segment_force(j, i);
            force = vec3_add(force, f);
        }

        // Contribution from segment (i, i+1)
        if i + 1 < n {
            let f = self.segment_force(i + 1, i); // force on i from i+1
            force = vec3_add(force, f);
        }

        force
    }

    /// Compute the elastic + damping force exerted on node `target` by the
    /// segment connecting `source` and `target`.
    fn segment_force(&self, source: usize, target: usize) -> [f64; 3] {
        let seg_idx = source.min(target);
        let l0 = self.rest_lengths[seg_idx];

        let pa = self.nodes[target].pos;
        let pb = self.nodes[source].pos;
        let va = self.nodes[target].vel;
        let vb = self.nodes[source].vel;

        let delta = vec3_sub(pb, pa);
        let dist = vec3_len(delta);
        if dist < 1e-12 {
            return [0.0; 3];
        }
        let dir = vec3_scale(delta, 1.0 / dist);

        // Spring: only pulls (no compression in standard cable models)
        let extension = dist - l0;
        let f_spring = self.stiffness * extension;

        // Damping along the segment
        let rel_vel = vec3_sub(vb, va);
        let f_damp = self.damping * vec3_dot(rel_vel, dir);

        let f_total = f_spring + f_damp;
        vec3_scale(dir, f_total)
    }

    /// Advance the simulation by one time step `dt` using semi-implicit Euler.
    ///
    /// # Parameters
    /// - `dt`      – time step \[s\]
    /// - `gravity` – gravitational acceleration vector \[m/s²\]
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.nodes.len();
        let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];

        // Accumulate elastic forces
        for i in 0..n {
            forces[i] = self.elastic_force(i);
        }

        // Update velocities (semi-implicit Euler)
        for i in 0..n {
            if self.nodes[i].fixed {
                continue;
            }
            let m = self.nodes[i].mass;
            let g_f = vec3_scale(gravity, m);
            let total = vec3_add(forces[i], g_f);
            let acc = vec3_scale(total, 1.0 / m);
            self.nodes[i].vel = vec3_add(self.nodes[i].vel, vec3_scale(acc, dt));
        }

        // Update positions
        for i in 0..n {
            if self.nodes[i].fixed {
                continue;
            }
            let v = self.nodes[i].vel;
            self.nodes[i].pos = vec3_add(self.nodes[i].pos, vec3_scale(v, dt));
        }
    }

    /// Tension magnitude at node `i`, estimated as the average of the tensions
    /// in the adjacent segments.
    pub fn tension_at_node(&self, i: usize) -> f64 {
        let n = self.nodes.len();
        let mut sum = 0.0;
        let mut count = 0usize;

        if i > 0 {
            let j = i - 1;
            let delta = vec3_sub(self.nodes[i].pos, self.nodes[j].pos);
            let dist = vec3_len(delta);
            let l0 = self.rest_lengths[j];
            let ext = (dist - l0).max(0.0);
            sum += self.stiffness * ext;
            count += 1;
        }

        if i + 1 < n {
            let delta = vec3_sub(self.nodes[i + 1].pos, self.nodes[i].pos);
            let dist = vec3_len(delta);
            let l0 = self.rest_lengths[i];
            let ext = (dist - l0).max(0.0);
            sum += self.stiffness * ext;
            count += 1;
        }

        if count == 0 { 0.0 } else { sum / count as f64 }
    }

    /// Total mechanical energy (kinetic + gravitational potential) of the cable.
    ///
    /// Gravitational potential is measured relative to z = 0.
    pub fn total_energy(&self, gravity: [f64; 3]) -> f64 {
        let g_mag = vec3_len(gravity);
        let mut energy = 0.0;
        for node in &self.nodes {
            // Kinetic
            energy += 0.5 * node.mass * vec3_dot(node.vel, node.vel);
            // Potential (height along gravity direction)
            // gravity points down; height = -pos · g_hat
            let g_hat = if g_mag > 1e-12 {
                vec3_scale(gravity, 1.0 / g_mag)
            } else {
                [0.0; 3]
            };
            let h = -vec3_dot(node.pos, g_hat);
            energy += node.mass * g_mag * h;
        }
        energy
    }

    /// Vertical droop at the cable midpoint compared to the first node.
    ///
    /// Returns (z_mid − z_first), which is negative when the cable sags.
    pub fn droop_at_midpoint(&self) -> f64 {
        if self.nodes.is_empty() {
            return 0.0;
        }
        let mid = self.nodes.len() / 2;
        self.nodes[mid].pos[2] - self.nodes[0].pos[2]
    }
}

// ---------------------------------------------------------------------------
// MooringLine
// ---------------------------------------------------------------------------

/// Mooring line model with seabed contact.
///
/// Nodes that fall below `seabed_z` are projected back to the seabed and
/// their vertical velocity is zeroed (inelastic seabed contact).
///
/// The axial stiffness `ea` is the product of Young's modulus E and
/// cross-sectional area A \[N\].
pub struct MooringLine {
    /// Node list.
    pub nodes: Vec<CableNode>,
    /// Z-coordinate of the flat seabed \[m\].
    pub seabed_z: f64,
    /// Rest lengths of segments \[m\].
    pub rest_lengths: Vec<f64>,
    /// Axial stiffness EA \[N\].
    pub ea: f64,
}

impl MooringLine {
    /// Create a mooring line with `n` nodes stretching from `anchor` to `fairlead`.
    ///
    /// The first node is fixed at the anchor; the last node is the fairlead
    /// (also initially fixed).
    pub fn new(
        anchor: [f64; 3],
        fairlead: [f64; 3],
        n: usize,
        mass_per_length: f64,
        ea: f64,
    ) -> Self {
        assert!(n >= 2, "MooringLine must have at least 2 nodes");
        let total_length = vec3_len(vec3_sub(fairlead, anchor));
        let seg_len = total_length / (n - 1) as f64;
        let node_mass = mass_per_length * total_length / n as f64;

        let mut nodes = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let pos = [
                anchor[0] + t * (fairlead[0] - anchor[0]),
                anchor[1] + t * (fairlead[1] - anchor[1]),
                anchor[2] + t * (fairlead[2] - anchor[2]),
            ];
            let mut node = CableNode::new(pos, node_mass);
            if i == 0 || i == n - 1 {
                node.fixed = true;
            }
            nodes.push(node);
        }

        let seabed_z = anchor[2];
        let rest_lengths = vec![seg_len; n - 1];

        Self {
            nodes,
            seabed_z,
            rest_lengths,
            ea,
        }
    }

    /// Advance the mooring line by one time step.
    ///
    /// Applies gravity and elastic spring forces, then projects nodes back to
    /// the seabed.
    pub fn step(&mut self, dt: f64) {
        let gravity = [0.0, 0.0, -G];
        let n = self.nodes.len();
        let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];

        // Elastic segment forces
        for seg in 0..(n - 1) {
            let i = seg;
            let j = seg + 1;
            let delta = vec3_sub(self.nodes[j].pos, self.nodes[i].pos);
            let dist = vec3_len(delta);
            let l0 = self.rest_lengths[seg];
            if dist < 1e-12 {
                continue;
            }
            let dir = vec3_scale(delta, 1.0 / dist);
            let extension = dist - l0;
            // Only tension (no compression)
            let f_mag = self.ea / l0 * extension.max(0.0);
            let f = vec3_scale(dir, f_mag);
            forces[i] = vec3_add(forces[i], f);
            forces[j] = vec3_sub(forces[j], f);
        }

        // Integrate
        for i in 0..n {
            if self.nodes[i].fixed {
                continue;
            }
            let m = self.nodes[i].mass;
            let g_f = vec3_scale(gravity, m);
            let total = vec3_add(forces[i], g_f);
            let acc = vec3_scale(total, 1.0 / m);
            self.nodes[i].vel = vec3_add(self.nodes[i].vel, vec3_scale(acc, dt));
            self.nodes[i].pos = vec3_add(self.nodes[i].pos, vec3_scale(self.nodes[i].vel, dt));

            // Seabed contact
            if self.nodes[i].pos[2] < self.seabed_z {
                self.nodes[i].pos[2] = self.seabed_z;
                if self.nodes[i].vel[2] < 0.0 {
                    self.nodes[i].vel[2] = 0.0;
                }
            }
        }
    }

    /// Tension force vector at the anchor node (node 0) from segment 0→1.
    ///
    /// Positive values indicate the anchor is being pulled upward / toward
    /// the fairlead.
    pub fn anchor_tension(&self) -> [f64; 3] {
        if self.nodes.len() < 2 {
            return [0.0; 3];
        }
        let delta = vec3_sub(self.nodes[1].pos, self.nodes[0].pos);
        let dist = vec3_len(delta);
        let l0 = self.rest_lengths[0];
        if dist < 1e-12 {
            return [0.0; 3];
        }
        let dir = vec3_scale(delta, 1.0 / dist);
        let extension = (dist - l0).max(0.0);
        let f_mag = self.ea / l0 * extension;
        vec3_scale(dir, f_mag)
    }

    /// Effective (submerged) weight of the mooring line in water.
    ///
    /// W_eff = Σ m_i · g − ρ_w · g · V_cable
    ///
    /// where the cable volume is approximated from the total chain/rope mass
    /// and a nominal steel density of 7850 kg/m³.
    pub fn effective_weight_in_water(&self, rho_water: f64) -> f64 {
        let total_mass: f64 = self.nodes.iter().map(|n| n.mass).sum();
        let total_length: f64 = self.rest_lengths.iter().sum();
        let rho_steel = 7850.0_f64;
        let volume = total_mass / rho_steel;
        let buoyancy = rho_water * G * volume;
        let weight = total_mass * G;
        // Also account for laid-on-seabed nodes having zero effective weight
        let n_on_seabed = self
            .nodes
            .iter()
            .filter(|n| (n.pos[2] - self.seabed_z).abs() < 1e-3)
            .count();
        let fraction_on_seabed = if total_length > 0.0 {
            n_on_seabed as f64 / self.nodes.len() as f64
        } else {
            0.0
        };
        (weight - buoyancy) * (1.0 - fraction_on_seabed)
    }
}

// ---------------------------------------------------------------------------
// Stand-alone utility functions
// ---------------------------------------------------------------------------

/// Midpoint sag of a parabolic wire (small-sag approximation).
///
/// For a wire with uniform weight per unit length w, horizontal tension H,
/// and span S, the midpoint sag is:
///
/// ```text
/// sag = w · S² / (8 · H)
/// ```
///
/// # Parameters
/// - `span`               – horizontal distance between supports \[m\]
/// - `weight_per_length`  – weight per unit length w = ρ_l · g \[N/m\]
/// - `horizontal_tension` – horizontal component of cable tension H \[N\]
pub fn wire_sag_midpoint(span: f64, weight_per_length: f64, horizontal_tension: f64) -> f64 {
    weight_per_length * span * span / (8.0 * horizontal_tension)
}

/// Natural frequency of mode `mode` of a taut cable (string vibration formula).
///
/// ```text
/// f_n = n / (2L) · sqrt(T / μ)
/// ```
///
/// where T is tension, μ is mass per unit length, L is length, and n is the
/// mode number (1 = fundamental).
///
/// # Parameters
/// - `tension`          – axial tension T \[N\]
/// - `mass_per_length`  – linear mass density μ \[kg/m\]
/// - `length`           – cable length L \[m\]
/// - `mode`             – mode number n (1 = fundamental, 2 = first overtone, …)
pub fn cable_eigenfrequency(tension: f64, mass_per_length: f64, length: f64, mode: usize) -> f64 {
    let n = mode as f64;
    n / (2.0 * length) * (tension / mass_per_length).sqrt()
}

// ---------------------------------------------------------------------------
// catenary_profile
// ---------------------------------------------------------------------------

/// Compute the 3-D positions of `n` evenly-spaced points along a catenary
/// hanging in the X-Z plane.
///
/// The catenary is parameterised so that x runs from 0 to `span` and the
/// lowest point sags below the chord by `sag` metres.
///
/// Returns `n` points `[x, 0.0, z]` where z is the height measured upward.
///
/// # Parameters
/// - `span` – horizontal distance between supports \[m\]
/// - `sag`  – maximum downward sag below the chord \[m\]
/// - `n`    – number of sample points (must be ≥ 2)
pub fn catenary_profile(span: f64, sag: f64, n: usize) -> Vec<[f64; 3]> {
    let n = n.max(2);
    // Catenary constant a from sag formula: sag = a*(cosh(L/(2a)) - 1)
    // solved with Newton-Raphson
    let half_span = span * 0.5;
    let mut a = if sag > 1e-12 {
        half_span * half_span / (2.0 * sag) // parabolic approx as initial guess
    } else {
        1e6 // near-taut cable
    };
    for _ in 0..200 {
        let r = half_span / a;
        let f = a * (r.cosh() - 1.0) - sag;
        let df = r.cosh() - 1.0 - r * r.sinh();
        if df.abs() < 1e-15 {
            break;
        }
        let da = -f / df;
        a += da;
        if a < 1e-9 {
            a = 1e-9;
        }
        if da.abs() < 1e-10 {
            break;
        }
    }
    // Endpoints at z = 0; midpoint sags to z = -sag.
    // z(x) = a * cosh((x - half_span)/a) - a * cosh(half_span/a)
    // At x=0 or x=span: cosh(±half_span/a) → 0 endpoint.
    // At x=half_span: cosh(0)=1, so z = a - a*cosh(half_span/a) = -sag.
    let endpoint_height = a * (half_span / a).cosh(); // = a + sag
    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let x = span * i as f64 / (n - 1) as f64;
        let z = a * ((x - half_span) / a).cosh() - endpoint_height;
        pts.push([x, 0.0, z]);
    }
    pts
}

// ---------------------------------------------------------------------------
// DrapeConstraint
// ---------------------------------------------------------------------------

/// Position-level inextensibility constraint between consecutive cable nodes.
///
/// Projects node positions to satisfy |p\[i+1\] - p\[i\]| = rest_length
/// (XPBD-style position correction, no velocity correction).
pub struct DrapeConstraint {
    /// Rest length between consecutive nodes \[m\].
    pub rest_length: f64,
    /// Stiffness compliance α = 1 / (k · dt²) for XPBD.
    pub compliance: f64,
}

impl DrapeConstraint {
    /// Create a new drape constraint with given rest length and stiffness.
    ///
    /// # Parameters
    /// - `rest_length` – target segment length \[m\]
    /// - `stiffness`   – spring stiffness (N/m); high values → inextensible
    pub fn new(rest_length: f64, stiffness: f64) -> Self {
        Self {
            rest_length,
            compliance: 1.0 / stiffness.max(1e-6),
        }
    }

    /// Project positions of two nodes so that their separation equals
    /// `rest_length`.
    ///
    /// Returns the position corrections `(delta_a, delta_b)` to apply to
    /// node A and node B respectively.
    ///
    /// # Parameters
    /// - `pos_a` – current position of node A \[m\]
    /// - `pos_b` – current position of node B \[m\]
    /// - `mass_a` – mass of node A \[kg\]; use f64::INFINITY for fixed node
    /// - `mass_b` – mass of node B \[kg\]; use f64::INFINITY for fixed node
    pub fn project(
        &self,
        pos_a: [f64; 3],
        pos_b: [f64; 3],
        mass_a: f64,
        mass_b: f64,
    ) -> ([f64; 3], [f64; 3]) {
        let delta = vec3_sub(pos_b, pos_a);
        let dist = vec3_len(delta);
        if dist < 1e-12 {
            return ([0.0; 3], [0.0; 3]);
        }
        let n = vec3_scale(delta, 1.0 / dist);
        let c = dist - self.rest_length; // constraint violation
        let w_a = if mass_a.is_finite() {
            1.0 / mass_a
        } else {
            0.0
        };
        let w_b = if mass_b.is_finite() {
            1.0 / mass_b
        } else {
            0.0
        };
        let w_sum = w_a + w_b + self.compliance;
        if w_sum < 1e-15 {
            return ([0.0; 3], [0.0; 3]);
        }
        let lambda = -c / w_sum;
        let da = vec3_scale(n, -lambda * w_a);
        let db = vec3_scale(n, lambda * w_b);
        (da, db)
    }

    /// Apply one iteration of position projection to two nodes in-place.
    ///
    /// Fixed nodes (node.fixed = true) are not moved.
    pub fn apply(&self, a: &mut CableNode, b: &mut CableNode) {
        let mass_a = if a.fixed { f64::INFINITY } else { a.mass };
        let mass_b = if b.fixed { f64::INFINITY } else { b.mass };
        let (da, db) = self.project(a.pos, b.pos, mass_a, mass_b);
        if !a.fixed {
            a.pos = vec3_add(a.pos, da);
        }
        if !b.fixed {
            b.pos = vec3_add(b.pos, db);
        }
    }
}

// ---------------------------------------------------------------------------
// CableWind
// ---------------------------------------------------------------------------

/// Wind/current load on cable segments using the Morison drag model.
///
/// The drag force per unit length on a cable element of diameter `d` is:
///
/// ```text
/// f_drag = 0.5 · ρ · Cd · d · |u_rel| · u_rel
/// ```
///
/// where u_rel is the relative velocity between wind/current and the cable.
pub struct CableWind {
    /// Fluid density \[kg/m³\] (air ≈ 1.225, sea water ≈ 1025).
    pub fluid_density: f64,
    /// Normal drag coefficient Cd \[-\].
    pub cd: f64,
    /// Cable outer diameter \[m\].
    pub diameter: f64,
    /// Wind/current velocity vector \[m/s\].
    pub velocity: [f64; 3],
}

impl CableWind {
    /// Create a wind force model for a round cable in air.
    ///
    /// # Parameters
    /// - `diameter`  – cable outer diameter \[m\]
    /// - `wind_vel`  – wind velocity \[m/s\]
    pub fn new_air(diameter: f64, wind_vel: [f64; 3]) -> Self {
        Self {
            fluid_density: 1.225,
            cd: 1.2,
            diameter,
            velocity: wind_vel,
        }
    }

    /// Create a current force model for a round cable in seawater.
    ///
    /// # Parameters
    /// - `diameter`    – cable outer diameter \[m\]
    /// - `current_vel` – current velocity \[m/s\]
    pub fn new_seawater(diameter: f64, current_vel: [f64; 3]) -> Self {
        Self {
            fluid_density: 1025.0,
            cd: 1.0,
            diameter,
            velocity: current_vel,
        }
    }

    /// Compute the drag force per unit length on a cable element.
    ///
    /// # Parameters
    /// - `node_vel` – velocity of the cable node (or element midpoint) \[m/s\]
    /// - `seg_len`  – length of the cable segment attributed to this node \[m\]
    ///
    /// Returns total force on the element \[N\] (= force/length × seg_len).
    pub fn force_on_element(&self, node_vel: [f64; 3], seg_len: f64) -> [f64; 3] {
        let u_rel = vec3_sub(self.velocity, node_vel);
        let speed = vec3_len(u_rel);
        if speed < 1e-12 {
            return [0.0; 3];
        }
        // Force per unit length = 0.5 * rho * Cd * d * speed^2 * direction
        let f_per_len = 0.5 * self.fluid_density * self.cd * self.diameter * speed * speed;
        let dir = vec3_scale(u_rel, 1.0 / speed);
        vec3_scale(dir, f_per_len * seg_len)
    }

    /// Apply wind forces to all free nodes in a [`CableSimulation`].
    ///
    /// Each node is assigned the segment length from the two adjacent segments,
    /// weighted by half from each side.
    pub fn apply_to_simulation(&self, sim: &mut CableSimulation, forces: &mut [[f64; 3]]) {
        let n = sim.nodes.len();
        for i in 0..n {
            if sim.nodes[i].fixed {
                continue;
            }
            // Effective segment length for this node
            let seg_len = {
                let mut s = 0.0;
                if i > 0 {
                    s += sim.rest_lengths[i - 1] * 0.5;
                }
                if i + 1 < n {
                    s += sim.rest_lengths[i] * 0.5;
                }
                s
            };
            let f = self.force_on_element(sim.nodes[i].vel, seg_len);
            forces[i] = vec3_add(forces[i], f);
        }
    }
}

// ---------------------------------------------------------------------------
// CableTension
// ---------------------------------------------------------------------------

/// Axial tension analysis for a lumped-mass cable simulation.
///
/// Computes the tension at each segment between consecutive nodes using
/// the spring model and optionally accounts for dynamic amplification.
pub struct CableTension {
    /// Axial spring stiffness used for the cable \[N/m\].
    pub ea_stiffness: f64,
}

impl CableTension {
    /// Create a tension analyser from an EA stiffness.
    ///
    /// # Parameters
    /// - `ea_stiffness` – axial stiffness EA \[N\]
    pub fn new(ea_stiffness: f64) -> Self {
        Self { ea_stiffness }
    }

    /// Compute axial tension in segment `seg` (between nodes `seg` and
    /// `seg + 1`) of a cable simulation.
    ///
    /// # Parameters
    /// - `sim` – cable simulation
    /// - `seg` – segment index (0 … n_nodes-2)
    ///
    /// Returns tension in N (zero for compressed segments).
    pub fn segment_tension(&self, sim: &CableSimulation, seg: usize) -> f64 {
        if seg + 1 >= sim.nodes.len() {
            return 0.0;
        }
        let pa = sim.nodes[seg].pos;
        let pb = sim.nodes[seg + 1].pos;
        let l0 = sim.rest_lengths[seg];
        let dist = vec3_len(vec3_sub(pb, pa));
        let ext = (dist - l0).max(0.0);
        // Use per-segment stiffness = EA / l0
        self.ea_stiffness / l0 * ext
    }

    /// Compute tension at all segments.
    ///
    /// Returns a vector of length `n_nodes - 1` with tensions in N.
    pub fn all_tensions(&self, sim: &CableSimulation) -> Vec<f64> {
        let n = sim.nodes.len();
        (0..n.saturating_sub(1))
            .map(|seg| self.segment_tension(sim, seg))
            .collect()
    }

    /// Maximum tension in the cable \[N\].
    pub fn max_tension(&self, sim: &CableSimulation) -> f64 {
        self.all_tensions(sim).into_iter().fold(0.0_f64, f64::max)
    }

    /// Mean tension over all segments \[N\].
    pub fn mean_tension(&self, sim: &CableSimulation) -> f64 {
        let t = self.all_tensions(sim);
        if t.is_empty() {
            return 0.0;
        }
        t.iter().sum::<f64>() / t.len() as f64
    }
}

// ---------------------------------------------------------------------------
// CableBundle
// ---------------------------------------------------------------------------

/// A bundle of parallel cables with twist-coupling between neighbouring cables.
///
/// The bundle is modelled as a set of [`CableSimulation`] objects sharing the
/// same end conditions, coupled by a lateral spring `coupling_stiffness`
/// representing the twist/banding effect.
pub struct CableBundle {
    /// Individual cable simulations in the bundle.
    pub cables: Vec<CableSimulation>,
    /// Lateral coupling stiffness between adjacent cables in the bundle \[N/m\].
    pub coupling_stiffness: f64,
    /// Nominal separation between cable centres in the bundle \[m\].
    pub separation: f64,
}

impl CableBundle {
    /// Create a bundle of `n_cables` identical cables.
    ///
    /// # Parameters
    /// - `n_cables`           – number of cables
    /// - `n_nodes`            – nodes per cable
    /// - `total_length`       – undeformed cable length \[m\]
    /// - `mass_per_length`    – linear mass density \[kg/m\]
    /// - `stiffness`          – axial stiffness \[N/m\]
    /// - `damping`            – viscous damping \[N·s/m\]
    /// - `coupling_stiffness` – lateral coupling stiffness \[N/m\]
    /// - `separation`         – inter-cable spacing \[m\]
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        n_cables: usize,
        n_nodes: usize,
        total_length: f64,
        mass_per_length: f64,
        stiffness: f64,
        damping: f64,
        coupling_stiffness: f64,
        separation: f64,
    ) -> Self {
        let cables = (0..n_cables.max(1))
            .map(|_| {
                CableSimulation::new(n_nodes, total_length, mass_per_length, stiffness, damping)
            })
            .collect();
        Self {
            cables,
            coupling_stiffness,
            separation,
        }
    }

    /// Number of cables in the bundle.
    pub fn n_cables(&self) -> usize {
        self.cables.len()
    }

    /// Advance the entire bundle by one time step.
    ///
    /// Each cable is stepped independently, then lateral coupling forces
    /// are applied between adjacent cables at each node level.
    ///
    /// # Parameters
    /// - `dt`      – time step \[s\]
    /// - `gravity` – gravitational acceleration \[m/s²\]
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        // Step all cables independently first
        for cable in &mut self.cables {
            cable.step(dt, gravity);
        }

        // Apply lateral coupling forces between adjacent cables
        let n_cables = self.cables.len();
        let n_nodes = if n_cables > 0 {
            self.cables[0].nodes.len()
        } else {
            0
        };

        for c in 0..n_cables.saturating_sub(1) {
            for i in 0..n_nodes {
                // Lateral offset between cable c and c+1 at node i
                // (only the X direction for simplicity; the bundle is coplanar in XZ)
                let pos_a = self.cables[c].nodes[i].pos;
                let pos_b = self.cables[c + 1].nodes[i].pos;
                let dx = pos_b[0] - pos_a[0];
                let dy = pos_b[1] - pos_a[1];
                // Desired separation along X (simplified bundle geometry)
                let desired = self.separation;
                let current_sep = (dx * dx + dy * dy).sqrt();
                let extension = current_sep - desired;
                if current_sep < 1e-12 {
                    continue;
                }
                let f_mag = self.coupling_stiffness * extension;
                let dir_x = dx / current_sep;
                let dir_y = dy / current_sep;

                if !self.cables[c].nodes[i].fixed {
                    let m = self.cables[c].nodes[i].mass;
                    let ax = f_mag * dir_x / m;
                    let ay = f_mag * dir_y / m;
                    self.cables[c].nodes[i].vel[0] += ax * dt;
                    self.cables[c].nodes[i].vel[1] += ay * dt;
                }
                if !self.cables[c + 1].nodes[i].fixed {
                    let m = self.cables[c + 1].nodes[i].mass;
                    let ax = -f_mag * dir_x / m;
                    let ay = -f_mag * dir_y / m;
                    self.cables[c + 1].nodes[i].vel[0] += ax * dt;
                    self.cables[c + 1].nodes[i].vel[1] += ay * dt;
                }
            }
        }
    }

    /// Total mechanical energy of all cables in the bundle.
    pub fn total_energy(&self, gravity: [f64; 3]) -> f64 {
        self.cables.iter().map(|c| c.total_energy(gravity)).sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── CatenarySegment ──────────────────────────────────────────────────────

    #[test]
    fn catenary_height_at_x_offset_equals_y_offset_plus_a() {
        let seg = CatenarySegment {
            a: 10.0,
            x_offset: 5.0,
            y_offset: 0.0,
        };
        // cosh(0) = 1
        let h = seg.height_at(5.0);
        assert!((h - 10.0).abs() < 1e-10, "h={h}");
    }

    #[test]
    fn catenary_height_symmetric_about_vertex() {
        let seg = CatenarySegment {
            a: 5.0,
            x_offset: 0.0,
            y_offset: 2.0,
        };
        let h_pos = seg.height_at(3.0);
        let h_neg = seg.height_at(-3.0);
        assert!(
            (h_pos - h_neg).abs() < 1e-12,
            "catenary not symmetric: {h_pos} vs {h_neg}"
        );
    }

    #[test]
    fn catenary_arc_length_positive() {
        let seg = CatenarySegment {
            a: 8.0,
            x_offset: 0.0,
            y_offset: 0.0,
        };
        let l = seg.arc_length(-5.0, 5.0);
        assert!(l > 10.0, "arc length {l} must exceed chord 10");
    }

    #[test]
    fn catenary_arc_length_symmetric() {
        let seg = CatenarySegment {
            a: 6.0,
            x_offset: 0.0,
            y_offset: 0.0,
        };
        let l1 = seg.arc_length(-4.0, 4.0);
        let l2 = seg.arc_length(-4.0, 4.0);
        assert!((l1 - l2).abs() < 1e-12);
    }

    #[test]
    fn catenary_max_sag_positive_for_symmetric_case() {
        let seg = CatenarySegment {
            a: 10.0,
            x_offset: 0.0,
            y_offset: 0.0,
        };
        let sag = seg.max_sag(-5.0, 5.0);
        // Chord is flat (both endpoints same height), catenary dips below
        assert!(sag >= 0.0, "sag={sag} must be non-negative");
    }

    #[test]
    fn catenary_tension_at_vertex_minimum() {
        let a = 12.0;
        let seg = CatenarySegment {
            a,
            x_offset: 0.0,
            y_offset: 0.0,
        };
        let t_vertex = seg.tension_at(0.0, 1.0);
        let t_side = seg.tension_at(5.0, 1.0);
        assert!(
            t_side >= t_vertex,
            "tension should increase away from vertex"
        );
    }

    #[test]
    fn catenary_from_endpoints_returns_segment() {
        // Level span: x1=0, y1=10, x2=20, y2=10, length=22
        let seg = CatenarySegment::from_endpoints(0.0, 10.0, 20.0, 10.0, 22.0);
        assert!(seg.a > 0.0, "a must be positive");
        // Height at midpoint should be near or below endpoint height (sag)
        let h_mid = seg.height_at(10.0);
        assert!(h_mid <= 10.0 + 1.0, "midpoint height unreasonably high");
    }

    #[test]
    fn catenary_from_endpoints_chord_too_long_returns_default() {
        // length < chord → should return default
        let seg = CatenarySegment::from_endpoints(0.0, 0.0, 10.0, 0.0, 5.0);
        assert_eq!(seg.a, 1.0);
    }

    #[test]
    fn catenary_arc_length_known_value() {
        // For a=1, arc from -1 to 1: L = 2*sinh(1) ≈ 2.3504
        let seg = CatenarySegment {
            a: 1.0,
            x_offset: 0.0,
            y_offset: 0.0,
        };
        let l = seg.arc_length(-1.0, 1.0);
        assert!((l - 2.0 * 1.0_f64.sinh()).abs() < 1e-10, "l={l}");
    }

    // ── CableNode ────────────────────────────────────────────────────────────

    #[test]
    fn cable_node_new_starts_stationary() {
        let node = CableNode::new([1.0, 2.0, 3.0], 5.0);
        assert!(!node.fixed);
        assert_eq!(node.vel, [0.0; 3]);
        assert_eq!(node.mass, 5.0);
    }

    #[test]
    fn cable_node_fixed_is_fixed() {
        let node = CableNode::new_fixed([0.0, 0.0, 0.0], 1.0);
        assert!(node.fixed);
    }

    // ── CableSimulation ──────────────────────────────────────────────────────

    #[test]
    fn cable_simulation_new_correct_node_count() {
        let sim = CableSimulation::new(5, 10.0, 1.0, 1000.0, 10.0);
        assert_eq!(sim.nodes.len(), 5);
        assert_eq!(sim.rest_lengths.len(), 4);
    }

    #[test]
    fn cable_simulation_first_node_fixed() {
        let sim = CableSimulation::new(4, 8.0, 1.0, 500.0, 5.0);
        assert!(sim.nodes[0].fixed);
    }

    #[test]
    fn cable_simulation_rest_lengths_uniform() {
        let sim = CableSimulation::new(5, 10.0, 1.0, 1000.0, 10.0);
        for rl in &sim.rest_lengths {
            assert!((*rl - 2.5).abs() < 1e-10, "rl={rl}");
        }
    }

    #[test]
    fn cable_simulation_step_moves_free_nodes() {
        let mut sim = CableSimulation::new(3, 6.0, 1.0, 100.0, 1.0);
        let z_before = sim.nodes[2].pos[2];
        let gravity = [0.0, 0.0, -9.81];
        sim.step(0.01, gravity);
        let z_after = sim.nodes[2].pos[2];
        // Under gravity the free node should have moved downward
        assert!(
            z_after < z_before,
            "free node should fall: before={z_before} after={z_after}"
        );
    }

    #[test]
    fn cable_simulation_fixed_node_does_not_move() {
        let mut sim = CableSimulation::new(3, 6.0, 1.0, 100.0, 1.0);
        let pos_before = sim.nodes[0].pos;
        sim.step(0.01, [0.0, 0.0, -9.81]);
        assert_eq!(sim.nodes[0].pos, pos_before, "fixed node must not move");
    }

    #[test]
    fn cable_simulation_elastic_force_zero_at_rest_length() {
        let sim = CableSimulation::new(3, 6.0, 1.0, 1000.0, 0.0);
        // All segments at rest length → elastic force on internal node should be zero
        // (nodes are vertical, symmetric)
        let f = sim.elastic_force(1);
        assert!(vec3_len(f) < 1e-8, "elastic force at rest: {:?}", f);
    }

    #[test]
    fn cable_tension_at_node_zero_when_at_rest() {
        let sim = CableSimulation::new(4, 8.0, 1.0, 1000.0, 0.0);
        let t = sim.tension_at_node(1);
        // At rest length, extension = 0
        assert!(t.abs() < 1e-10, "tension={t}");
    }

    #[test]
    fn cable_tension_at_node_positive_when_stretched() {
        let mut sim = CableSimulation::new(3, 6.0, 1.0, 1000.0, 0.0);
        // Stretch the last node
        sim.nodes[2].pos[2] -= 2.0;
        let t = sim.tension_at_node(1);
        assert!(t > 0.0, "tension should be positive when stretched: {t}");
    }

    #[test]
    fn cable_total_energy_positive() {
        let mut sim = CableSimulation::new(5, 10.0, 1.0, 1000.0, 1.0);
        sim.nodes[4].vel = [0.0, 0.0, 1.0];
        let e = sim.total_energy([0.0, 0.0, -9.81]);
        // Energy should be finite (may be negative due to gravity potential reference)
        assert!(e.is_finite(), "energy={e}");
    }

    #[test]
    fn cable_droop_at_midpoint_is_negative() {
        let sim = CableSimulation::new(5, 10.0, 1.0, 1000.0, 1.0);
        let droop = sim.droop_at_midpoint();
        // Midpoint is below first node (which is at z=0)
        assert!(droop < 0.0, "droop={droop} should be negative (sagging)");
    }

    // ── MooringLine ──────────────────────────────────────────────────────────

    #[test]
    fn mooring_line_new_correct_node_count() {
        let ml = MooringLine::new([0.0, 0.0, -100.0], [50.0, 0.0, -10.0], 10, 100.0, 1e8);
        assert_eq!(ml.nodes.len(), 10);
        assert_eq!(ml.rest_lengths.len(), 9);
    }

    #[test]
    fn mooring_anchor_and_fairlead_fixed() {
        let ml = MooringLine::new([0.0, 0.0, -100.0], [50.0, 0.0, -10.0], 8, 80.0, 1e8);
        assert!(ml.nodes[0].fixed);
        assert!(ml.nodes[7].fixed);
    }

    #[test]
    fn mooring_step_does_not_panic() {
        let mut ml = MooringLine::new([0.0, 0.0, -100.0], [50.0, 0.0, -10.0], 6, 60.0, 1e7);
        for _ in 0..20 {
            ml.step(0.01);
        }
    }

    #[test]
    fn mooring_anchor_tension_finite() {
        let mut ml = MooringLine::new([0.0, 0.0, -100.0], [50.0, 0.0, -10.0], 6, 60.0, 1e7);
        ml.step(0.01);
        let t = ml.anchor_tension();
        for v in t {
            assert!(v.is_finite(), "tension component {v} must be finite");
        }
    }

    #[test]
    fn mooring_effective_weight_positive_in_air() {
        let ml = MooringLine::new([0.0, 0.0, 0.0], [50.0, 0.0, 0.0], 6, 60.0, 1e7);
        // rho_water = 0 (in air)
        let w = ml.effective_weight_in_water(0.0);
        assert!(w >= 0.0, "weight in air should be non-negative: {w}");
    }

    #[test]
    fn mooring_effective_weight_reduced_in_water() {
        let ml = MooringLine::new([0.0, 0.0, -100.0], [50.0, 0.0, -10.0], 6, 60.0, 1e7);
        let w_air = ml.effective_weight_in_water(0.0);
        let w_water = ml.effective_weight_in_water(1025.0);
        assert!(
            w_water <= w_air,
            "submerged weight should be less: air={w_air} water={w_water}"
        );
    }

    #[test]
    fn mooring_seabed_contact_prevents_nodes_below_seabed() {
        let mut ml = MooringLine::new([0.0, 0.0, -5.0], [10.0, 0.0, 0.0], 5, 10.0, 1e6);
        // Run many steps and check no node goes below seabed_z
        for _ in 0..100 {
            ml.step(0.01);
        }
        for node in &ml.nodes {
            assert!(
                node.pos[2] >= ml.seabed_z - 1e-9,
                "node below seabed: z={}",
                node.pos[2]
            );
        }
    }

    // ── Utility functions ────────────────────────────────────────────────────

    #[test]
    fn wire_sag_midpoint_known_value() {
        // w=1 N/m, S=10 m, H=100 N → sag = 1*100/800 = 0.125
        let sag = wire_sag_midpoint(10.0, 1.0, 100.0);
        assert!((sag - 0.125).abs() < 1e-10, "sag={sag}");
    }

    #[test]
    fn wire_sag_increases_with_span() {
        let sag1 = wire_sag_midpoint(10.0, 1.0, 100.0);
        let sag2 = wire_sag_midpoint(20.0, 1.0, 100.0);
        assert!(sag2 > sag1, "sag should increase with span");
    }

    #[test]
    fn wire_sag_decreases_with_tension() {
        let sag1 = wire_sag_midpoint(10.0, 1.0, 100.0);
        let sag2 = wire_sag_midpoint(10.0, 1.0, 1000.0);
        assert!(sag2 < sag1, "sag should decrease with tension");
    }

    #[test]
    fn cable_eigenfrequency_fundamental_positive() {
        let f = cable_eigenfrequency(1000.0, 0.5, 10.0, 1);
        assert!(f > 0.0, "frequency must be positive: {f}");
    }

    #[test]
    fn cable_eigenfrequency_mode_doubles() {
        // f_2 = 2 * f_1
        let f1 = cable_eigenfrequency(1000.0, 0.5, 10.0, 1);
        let f2 = cable_eigenfrequency(1000.0, 0.5, 10.0, 2);
        assert!((f2 - 2.0 * f1).abs() < 1e-10, "f2={f2}, 2*f1={}", 2.0 * f1);
    }

    #[test]
    fn cable_eigenfrequency_increases_with_tension() {
        let f_low = cable_eigenfrequency(100.0, 0.5, 10.0, 1);
        let f_high = cable_eigenfrequency(10000.0, 0.5, 10.0, 1);
        assert!(f_high > f_low, "frequency should increase with tension");
    }

    #[test]
    fn cable_eigenfrequency_decreases_with_length() {
        let f_short = cable_eigenfrequency(1000.0, 0.5, 5.0, 1);
        let f_long = cable_eigenfrequency(1000.0, 0.5, 20.0, 1);
        assert!(f_long < f_short, "frequency should decrease with length");
    }

    #[test]
    fn cable_eigenfrequency_known_value() {
        // T=100N, μ=1kg/m, L=1m, n=1 → f = 1/(2*1)*sqrt(100/1) = 5 Hz
        let f = cable_eigenfrequency(100.0, 1.0, 1.0, 1);
        assert!((f - 5.0).abs() < 1e-10, "f={f}, expected 5.0");
    }

    // ── Vector helpers ───────────────────────────────────────────────────────

    #[test]
    fn vec3_add_correct() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        assert_eq!(vec3_add(a, b), [5.0, 7.0, 9.0]);
    }

    #[test]
    fn vec3_sub_correct() {
        let a = [4.0, 5.0, 6.0];
        let b = [1.0, 2.0, 3.0];
        assert_eq!(vec3_sub(a, b), [3.0, 3.0, 3.0]);
    }

    #[test]
    fn vec3_dot_correct() {
        assert!((vec3_dot([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!(vec3_dot([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]).abs() < 1e-12);
    }

    #[test]
    fn vec3_len_unit() {
        assert!((vec3_len([1.0, 0.0, 0.0]) - 1.0).abs() < 1e-12);
        assert!((vec3_len([0.0, 3.0, 4.0]) - 5.0).abs() < 1e-12);
    }

    #[test]
    fn vec3_norm_unit_vector() {
        let n = vec3_norm([3.0, 0.0, 4.0]);
        assert!((vec3_len(n) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn vec3_norm_zero_returns_zero() {
        let n = vec3_norm([0.0, 0.0, 0.0]);
        assert_eq!(n, [0.0; 3]);
    }

    // ── catenary_profile ─────────────────────────────────────────────────────

    #[test]
    fn catenary_profile_correct_count() {
        let pts = catenary_profile(10.0, 0.5, 11);
        assert_eq!(pts.len(), 11);
    }

    #[test]
    fn catenary_profile_first_x_zero() {
        let pts = catenary_profile(10.0, 0.5, 5);
        assert!((pts[0][0]).abs() < 1e-10);
    }

    #[test]
    fn catenary_profile_last_x_equals_span() {
        let pts = catenary_profile(20.0, 1.0, 7);
        assert!((pts[6][0] - 20.0).abs() < 1e-10, "last x = {}", pts[6][0]);
    }

    #[test]
    fn catenary_profile_y_zero_for_planar_catenary() {
        let pts = catenary_profile(10.0, 0.5, 5);
        for p in &pts {
            assert!(p[1].abs() < 1e-10, "y should be zero");
        }
    }

    #[test]
    fn catenary_profile_midpoint_z_lowest() {
        let pts = catenary_profile(10.0, 0.5, 11);
        let z_mid = pts[5][2];
        let z_end = pts[0][2];
        // Midpoint should be at or below the endpoints
        assert!(z_mid <= z_end + 1e-10, "midpoint z={z_mid} end z={z_end}");
    }

    #[test]
    fn catenary_profile_sag_magnitude_approx_correct() {
        // For small sag the midpoint z should be approximately -sag
        let sag = 0.5;
        let pts = catenary_profile(10.0, sag, 101);
        let z_min = pts.iter().map(|p| p[2]).fold(f64::INFINITY, f64::min);
        assert!(
            (z_min + sag).abs() < sag * 0.05,
            "z_min={z_min} expected ≈ -{sag}"
        );
    }

    // ── DrapeConstraint ──────────────────────────────────────────────────────

    #[test]
    fn drape_constraint_no_correction_at_rest_length() {
        let c = DrapeConstraint::new(1.0, 1e6);
        let (da, db) = c.project([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 1.0);
        assert!(vec3_len(da) < 1e-6, "da={:?}", da);
        assert!(vec3_len(db) < 1e-6, "db={:?}", db);
    }

    #[test]
    fn drape_constraint_corrects_elongated_segment() {
        let c = DrapeConstraint::new(1.0, 1e9);
        let (da, db) = c.project([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 1.0);
        // Node B should be pulled back (negative x direction)
        assert!(db[0] < 0.0, "db.x should be negative");
        assert!(da[0] > 0.0, "da.x should be positive");
    }

    #[test]
    fn drape_constraint_fixed_node_not_moved() {
        let c = DrapeConstraint::new(1.0, 1e9);
        let mut a = CableNode::new_fixed([0.0, 0.0, 0.0], 1.0);
        let mut b = CableNode::new([2.0, 0.0, 0.0], 1.0);
        c.apply(&mut a, &mut b);
        assert_eq!(a.pos, [0.0, 0.0, 0.0], "fixed node must not move");
    }

    #[test]
    fn drape_constraint_moves_free_node() {
        let c = DrapeConstraint::new(1.0, 1e9);
        let mut a = CableNode::new_fixed([0.0, 0.0, 0.0], 1.0);
        let mut b = CableNode::new([2.0, 0.0, 0.0], 1.0);
        let pos_before = b.pos;
        c.apply(&mut a, &mut b);
        assert!(
            vec3_len(vec3_sub(b.pos, pos_before)) > 1e-6,
            "free node should move after constraint"
        );
    }

    // ── CableWind ────────────────────────────────────────────────────────────

    #[test]
    fn cable_wind_force_zero_when_no_wind() {
        let w = CableWind {
            fluid_density: 1.225,
            cd: 1.2,
            diameter: 0.05,
            velocity: [0.0; 3],
        };
        let f = w.force_on_element([0.0; 3], 1.0);
        assert!(vec3_len(f) < 1e-12, "no force when no wind");
    }

    #[test]
    fn cable_wind_force_positive_in_wind_direction() {
        let w = CableWind::new_air(0.05, [10.0, 0.0, 0.0]);
        let f = w.force_on_element([0.0; 3], 1.0);
        assert!(f[0] > 0.0, "force should be in x direction");
    }

    #[test]
    fn cable_wind_force_increases_with_diameter() {
        let w1 = CableWind {
            fluid_density: 1.225,
            cd: 1.2,
            diameter: 0.05,
            velocity: [10.0, 0.0, 0.0],
        };
        let w2 = CableWind {
            fluid_density: 1.225,
            cd: 1.2,
            diameter: 0.10,
            velocity: [10.0, 0.0, 0.0],
        };
        let f1 = vec3_len(w1.force_on_element([0.0; 3], 1.0));
        let f2 = vec3_len(w2.force_on_element([0.0; 3], 1.0));
        assert!(
            f2 > f1,
            "larger diameter should give more drag: f1={f1} f2={f2}"
        );
    }

    #[test]
    fn cable_wind_seawater_density_greater_than_air() {
        let w_air = CableWind::new_air(0.05, [1.0, 0.0, 0.0]);
        let w_water = CableWind::new_seawater(0.05, [1.0, 0.0, 0.0]);
        assert!(w_water.fluid_density > w_air.fluid_density);
    }

    // ── CableTension ─────────────────────────────────────────────────────────

    #[test]
    fn cable_tension_zero_at_rest_length() {
        let sim = CableSimulation::new(4, 8.0, 1.0, 1000.0, 0.0);
        let ct = CableTension::new(1e6);
        assert!(ct.max_tension(&sim).abs() < 1e-8, "no tension at rest");
    }

    #[test]
    fn cable_tension_positive_when_stretched() {
        let mut sim = CableSimulation::new(3, 6.0, 1.0, 1000.0, 0.0);
        sim.nodes[2].pos[2] -= 3.0; // pull node 2 down
        let ct = CableTension::new(1e6);
        let t = ct.max_tension(&sim);
        assert!(t > 0.0, "tension should be positive: t={t}");
    }

    #[test]
    fn cable_tension_all_tensions_correct_count() {
        let sim = CableSimulation::new(5, 10.0, 1.0, 500.0, 0.0);
        let ct = CableTension::new(1e5);
        let tensions = ct.all_tensions(&sim);
        assert_eq!(tensions.len(), 4, "n-1 segments for n nodes");
    }

    #[test]
    fn cable_tension_mean_positive_when_stretched() {
        let mut sim = CableSimulation::new(4, 8.0, 1.0, 1000.0, 0.0);
        // Stretch all segments
        for i in 1..4 {
            sim.nodes[i].pos[2] -= 1.0 * i as f64;
        }
        let ct = CableTension::new(1e6);
        let mean = ct.mean_tension(&sim);
        assert!(mean > 0.0, "mean tension={mean}");
    }

    // ── CableBundle ──────────────────────────────────────────────────────────

    #[test]
    fn cable_bundle_correct_cable_count() {
        let bundle = CableBundle::new(3, 5, 10.0, 1.0, 500.0, 5.0, 100.0, 0.1);
        assert_eq!(bundle.n_cables(), 3);
    }

    #[test]
    fn cable_bundle_step_does_not_panic() {
        let mut bundle = CableBundle::new(2, 4, 8.0, 1.0, 500.0, 5.0, 100.0, 0.05);
        bundle.step(0.01, [0.0, 0.0, -9.81]);
    }

    #[test]
    fn cable_bundle_total_energy_finite() {
        let bundle = CableBundle::new(2, 4, 8.0, 1.0, 500.0, 5.0, 100.0, 0.05);
        let e = bundle.total_energy([0.0, 0.0, -9.81]);
        assert!(e.is_finite(), "energy={e}");
    }

    #[test]
    fn cable_bundle_single_cable_same_as_sim() {
        let bundle = CableBundle::new(1, 4, 8.0, 1.0, 500.0, 5.0, 100.0, 0.05);
        assert_eq!(bundle.cables[0].nodes.len(), 4);
    }
}
