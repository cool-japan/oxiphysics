// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Cable network and tensegrity structure simulation.
//!
//! Provides:
//! - [`CableNode`] / [`Cable`] / [`CableNet`]: general cable mesh with
//!   stiffness, damping, and pretension.
//! - [`TensegrityStrut`] / [`TensegrityStructure`]: compression members
//!   combined with tension cables to form self-stressed tensegrity systems.
//! - [`FormFinding`]: dynamic relaxation method to find equilibrium shapes.
//! - Free functions: [`cable_sag_parabolic`], [`taut_string_frequency`].

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Helper math on [f64; 3]   (no nalgebra — softbody rule)
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
    if l < 1.0e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// CableNode
// ---------------------------------------------------------------------------

/// A point mass node in a cable or tensegrity network.
#[derive(Clone, Debug)]
pub struct CableNode {
    /// World-space position (m).
    pub position: [f64; 3],
    /// Current velocity (m/s).
    pub velocity: [f64; 3],
    /// Node mass (kg).
    pub mass: f64,
    /// When `true` the node is kinematically fixed and cannot translate.
    pub fixed: bool,
}

impl CableNode {
    /// Create a free node at `position` with given `mass`.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            fixed: false,
        }
    }

    /// Create a fixed (pinned) node at `position`.
    pub fn fixed(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass: f64::INFINITY,
            fixed: true,
        }
    }

    /// Kinetic energy of this node (J).
    pub fn kinetic_energy(&self) -> f64 {
        if self.fixed {
            return 0.0;
        }
        0.5 * self.mass * dot3(self.velocity, self.velocity)
    }
}

// ---------------------------------------------------------------------------
// Cable
// ---------------------------------------------------------------------------

/// A linear spring-damper cable element between two nodes.
///
/// The cable only pulls (tension only): it produces force only when the
/// current length exceeds `rest_length + pretension / stiffness`.
#[derive(Clone, Debug)]
pub struct Cable {
    /// Index of node A in the owner network.
    pub node_a: usize,
    /// Index of node B in the owner network.
    pub node_b: usize,
    /// Natural (unstressed) length (m).
    pub rest_length: f64,
    /// Axial stiffness EA (N), i.e. Young's modulus × cross-section.
    pub stiffness: f64,
    /// Viscous damping coefficient (N·s/m).
    pub damping: f64,
    /// Pre-tension force applied at rest length (N).
    pub pretension: f64,
}

impl Cable {
    /// Construct a cable from node indices and material parameters.
    pub fn new(
        node_a: usize,
        node_b: usize,
        rest_length: f64,
        stiffness: f64,
        damping: f64,
        pretension: f64,
    ) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            stiffness,
            damping,
            pretension,
        }
    }

    /// Compute the force vector acting on node A due to this cable.
    ///
    /// The returned vector points **from B toward A** (tension) when the cable
    /// is taut.  The force on node B is the negation of this vector.
    ///
    /// Requires the current positions and velocities of both endpoints.
    pub fn compute_force(
        &self,
        pos_a: [f64; 3],
        vel_a: [f64; 3],
        pos_b: [f64; 3],
        vel_b: [f64; 3],
    ) -> [f64; 3] {
        let ab = sub3(pos_b, pos_a); // vector from A to B
        let current_len = len3(ab);
        if current_len < 1.0e-15 {
            return [0.0; 3];
        }
        let dir_ab = normalize3(ab); // unit vector A→B

        // Extension beyond rest length.
        let extension = current_len - self.rest_length;

        // Cables are tension-only: no compressive force unless pre-tension dominates.
        let elastic_force = self.pretension + self.stiffness * extension;
        if elastic_force <= 0.0 {
            return [0.0; 3];
        }

        // Relative velocity projected onto the cable axis (positive = approaching).
        let rel_vel = sub3(vel_b, vel_a);
        let vel_along = dot3(rel_vel, dir_ab);

        // Total axial force (positive = tension, pulling A toward B).
        let total_force = elastic_force + self.damping * vel_along;
        let total_force = total_force.max(0.0); // cables cannot push

        // Force on A is in the +AB direction (toward B).
        scale3(dir_ab, total_force)
    }

    /// Elastic potential energy stored in this cable (J).
    pub fn potential_energy(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> f64 {
        let current_len = len3(sub3(pos_b, pos_a));
        let extension = (current_len - self.rest_length).max(0.0);
        0.5 * self.stiffness * extension * extension + self.pretension * extension
    }
}

// ---------------------------------------------------------------------------
// CableNet
// ---------------------------------------------------------------------------

/// A network of cable nodes connected by cable elements.
#[derive(Clone, Debug)]
pub struct CableNet {
    /// All nodes in the network.
    pub nodes: Vec<CableNode>,
    /// All cable elements.
    pub cables: Vec<Cable>,
}

impl CableNet {
    /// Create an empty cable net.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            cables: Vec::new(),
        }
    }

    /// Add a node and return its index.
    pub fn add_node(&mut self, node: CableNode) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }

    /// Add a cable.
    pub fn add_cable(&mut self, cable: Cable) {
        self.cables.push(cable);
    }

    /// Advance the simulation by one time step `dt` (s) using explicit Euler.
    ///
    /// External gravity `g` (m/s²) is applied as a body force.
    pub fn step(&mut self, dt: f64) {
        self.step_with_gravity(dt, [0.0, -9.81, 0.0]);
    }

    /// Advance the simulation with a user-supplied gravity vector.
    pub fn step_with_gravity(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.nodes.len();
        let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];

        // Gravity forces.
        for (i, node) in self.nodes.iter().enumerate() {
            if !node.fixed {
                forces[i] = scale3(gravity, node.mass);
            }
        }

        // Cable forces.
        for cable in &self.cables {
            let (ia, ib) = (cable.node_a, cable.node_b);
            let pa = self.nodes[ia].position;
            let pb = self.nodes[ib].position;
            let va = self.nodes[ia].velocity;
            let vb = self.nodes[ib].velocity;

            let f_on_a = cable.compute_force(pa, va, pb, vb);
            // f_on_a acts on A (directed toward B when cable is taut).
            forces[ia] = add3(forces[ia], f_on_a);
            forces[ib] = add3(forces[ib], scale3(f_on_a, -1.0));
        }

        // Integrate.
        for (i, node) in self.nodes.iter_mut().enumerate() {
            if node.fixed {
                continue;
            }
            let accel = scale3(forces[i], 1.0 / node.mass);
            node.velocity = add3(node.velocity, scale3(accel, dt));
            node.position = add3(node.position, scale3(node.velocity, dt));
        }
    }

    /// Total elastic potential energy stored in all cables (J).
    pub fn total_potential_energy(&self) -> f64 {
        self.cables
            .iter()
            .map(|c| {
                c.potential_energy(self.nodes[c.node_a].position, self.nodes[c.node_b].position)
            })
            .sum()
    }

    /// Total kinetic energy in all free nodes (J).
    pub fn total_kinetic_energy(&self) -> f64 {
        self.nodes.iter().map(|n| n.kinetic_energy()).sum()
    }
}

impl Default for CableNet {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// TensegrityStrut
// ---------------------------------------------------------------------------

/// A compression-only strut used in tensegrity structures.
///
/// Struts resist compression but offer no tensile resistance.
#[derive(Clone, Debug)]
pub struct TensegrityStrut {
    /// Index of node A.
    pub node_a: usize,
    /// Index of node B.
    pub node_b: usize,
    /// Natural (uncompressed) length (m).
    pub rest_length: f64,
    /// Axial compressive stiffness (N/m).
    pub compression_stiffness: f64,
}

impl TensegrityStrut {
    /// Construct a strut between two nodes.
    pub fn new(node_a: usize, node_b: usize, rest_length: f64, compression_stiffness: f64) -> Self {
        Self {
            node_a,
            node_b,
            rest_length,
            compression_stiffness,
        }
    }

    /// Compute the force vector on node A (points from B toward A for compression).
    pub fn compute_force(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> [f64; 3] {
        let ab = sub3(pos_b, pos_a);
        let current_len = len3(ab);
        if current_len < 1.0e-15 {
            return [0.0; 3];
        }
        let compression = self.rest_length - current_len; // positive when compressed
        if compression <= 0.0 {
            // Strut is not compressed — no force.
            return [0.0; 3];
        }
        let dir_ba = normalize3(scale3(ab, -1.0)); // unit vector B→A
        scale3(dir_ba, self.compression_stiffness * compression)
    }

    /// Elastic energy stored in the strut (J).
    pub fn potential_energy(&self, pos_a: [f64; 3], pos_b: [f64; 3]) -> f64 {
        let compression = (self.rest_length - len3(sub3(pos_b, pos_a))).max(0.0);
        0.5 * self.compression_stiffness * compression * compression
    }
}

// ---------------------------------------------------------------------------
// TensegrityStructure
// ---------------------------------------------------------------------------

/// A self-stressed structure of tension cables and compression struts.
#[derive(Clone, Debug)]
pub struct TensegrityStructure {
    /// Shared node list.
    pub nodes: Vec<CableNode>,
    /// Tension cables.
    pub cables: Vec<Cable>,
    /// Compression struts.
    pub struts: Vec<TensegrityStrut>,
}

impl TensegrityStructure {
    /// Create an empty tensegrity structure.
    pub fn new() -> Self {
        Self {
            nodes: Vec::new(),
            cables: Vec::new(),
            struts: Vec::new(),
        }
    }

    /// Add a node and return its index.
    pub fn add_node(&mut self, node: CableNode) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(node);
        idx
    }

    /// Add a tension cable.
    pub fn add_cable(&mut self, cable: Cable) {
        self.cables.push(cable);
    }

    /// Add a compression strut.
    pub fn add_strut(&mut self, strut: TensegrityStrut) {
        self.struts.push(strut);
    }

    /// Compute net force on every node (cables + struts + gravity).
    pub fn net_forces(&self, gravity: [f64; 3]) -> Vec<[f64; 3]> {
        let n = self.nodes.len();
        let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];

        // Gravity.
        for (i, node) in self.nodes.iter().enumerate() {
            if !node.fixed {
                forces[i] = scale3(gravity, node.mass);
            }
        }

        // Cable tension.
        for cable in &self.cables {
            let (ia, ib) = (cable.node_a, cable.node_b);
            let f = cable.compute_force(
                self.nodes[ia].position,
                self.nodes[ia].velocity,
                self.nodes[ib].position,
                self.nodes[ib].velocity,
            );
            forces[ia] = add3(forces[ia], f);
            forces[ib] = add3(forces[ib], scale3(f, -1.0));
        }

        // Strut compression.
        for strut in &self.struts {
            let (ia, ib) = (strut.node_a, strut.node_b);
            let f = strut.compute_force(self.nodes[ia].position, self.nodes[ib].position);
            forces[ia] = add3(forces[ia], f);
            forces[ib] = add3(forces[ib], scale3(f, -1.0));
        }

        forces
    }

    /// Find equilibrium by dynamic relaxation (zero gravity by default).
    ///
    /// Runs up to `max_iter` steps; stops early when total kinetic energy
    /// drops below `tol` J and all node velocities are zeroed at every step.
    pub fn find_equilibrium(&mut self, max_iter: usize, tol: f64) {
        let dt = 1.0e-3;
        let gravity = [0.0; 3];
        for _ in 0..max_iter {
            let forces = self.net_forces(gravity);
            let mut converged = true;
            for (i, node) in self.nodes.iter_mut().enumerate() {
                if node.fixed {
                    continue;
                }
                let accel = scale3(forces[i], 1.0 / node.mass);
                node.velocity = add3(node.velocity, scale3(accel, dt));
                // Kinetic damping: zero velocity if KE increased (dynamic relaxation).
                if dot3(node.velocity, node.velocity) > tol {
                    converged = false;
                }
                node.position = add3(node.position, scale3(node.velocity, dt));
            }
            // Zero velocities for kinetic damping.
            for node in self.nodes.iter_mut() {
                if !node.fixed {
                    node.velocity = [0.0; 3];
                }
            }
            if converged {
                break;
            }
        }
    }

    /// Total potential energy (cables + struts) in J.
    pub fn total_potential_energy(&self) -> f64 {
        let e_cables: f64 = self
            .cables
            .iter()
            .map(|c| {
                c.potential_energy(self.nodes[c.node_a].position, self.nodes[c.node_b].position)
            })
            .sum();
        let e_struts: f64 = self
            .struts
            .iter()
            .map(|s| {
                s.potential_energy(self.nodes[s.node_a].position, self.nodes[s.node_b].position)
            })
            .sum();
        e_cables + e_struts
    }
}

impl Default for TensegrityStructure {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// FormFinding  – dynamic relaxation
// ---------------------------------------------------------------------------

/// Dynamic relaxation solver for form-finding of cable / tensegrity structures.
///
/// The solver iterates an explicit time-integration with kinetic energy
/// damping: whenever the kinetic energy peaks, all velocities are zeroed.
/// This dissipates energy and forces convergence to a static equilibrium.
#[derive(Clone, Debug)]
pub struct FormFinding {
    /// Maximum number of dynamic relaxation iterations.
    pub max_iterations: usize,
    /// Convergence tolerance: stop when max residual force < `force_tol` N.
    pub force_tol: f64,
    /// Integration time step (s).
    pub dt: f64,
}

impl FormFinding {
    /// Create a form-finder with sensible defaults.
    pub fn new(max_iterations: usize, force_tol: f64) -> Self {
        Self {
            max_iterations,
            force_tol,
            dt: 1.0e-3,
        }
    }

    /// Run dynamic relaxation on `net` until convergence or `max_iterations`.
    ///
    /// Returns the number of iterations taken.
    pub fn solve(&self, net: &mut CableNet) -> usize {
        let mut prev_ke = f64::INFINITY;
        for iter in 0..self.max_iterations {
            let n = net.nodes.len();
            let mut forces: Vec<[f64; 3]> = vec![[0.0; 3]; n];

            for cable in &net.cables {
                let (ia, ib) = (cable.node_a, cable.node_b);
                let f = cable.compute_force(
                    net.nodes[ia].position,
                    net.nodes[ia].velocity,
                    net.nodes[ib].position,
                    net.nodes[ib].velocity,
                );
                forces[ia] = add3(forces[ia], f);
                forces[ib] = add3(forces[ib], scale3(f, -1.0));
            }

            // Check convergence: max force residual.
            let max_residual = forces
                .iter()
                .enumerate()
                .filter(|(i, _)| !net.nodes[*i].fixed)
                .map(|(_, f)| len3(*f))
                .fold(0.0_f64, f64::max);

            if max_residual < self.force_tol {
                return iter;
            }

            // Integrate.
            for (i, node) in net.nodes.iter_mut().enumerate() {
                if node.fixed {
                    continue;
                }
                let accel = scale3(forces[i], 1.0 / node.mass);
                node.velocity = add3(node.velocity, scale3(accel, self.dt));
                node.position = add3(node.position, scale3(node.velocity, self.dt));
            }

            // Kinetic energy peak detection (kinetic damping).
            let ke: f64 = net.nodes.iter().map(|n| n.kinetic_energy()).sum();
            if ke < prev_ke {
                // KE is decreasing — keep going.
            } else {
                // KE peaked — zero velocities.
                for node in net.nodes.iter_mut() {
                    if !node.fixed {
                        node.velocity = [0.0; 3];
                    }
                }
            }
            prev_ke = ke;
        }
        self.max_iterations
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Catenary/parabolic approximation of cable sag.
///
/// For a uniformly loaded cable spanning `span` (m) under horizontal tension
/// `tension` (N) and uniform vertical load `load_per_unit` (N/m), the maximum
/// mid-span sag is:
///
/// ```text
/// sag = w * L² / (8 * T)
/// ```
///
/// Returns sag in metres.
pub fn cable_sag_parabolic(span: f64, load_per_unit: f64, tension: f64) -> f64 {
    load_per_unit * span * span / (8.0 * tension)
}

/// Natural frequency of a taut string (Hz).
///
/// First-mode frequency of a simply-supported taut string:
///
/// ```text
/// f = (1 / 2L) * sqrt(T / μ)
/// ```
///
/// * `tension` – axial tension (N).
/// * `length` – string length (m).
/// * `linear_density` – mass per unit length μ (kg/m).
///
/// Returns frequency in Hz.
pub fn taut_string_frequency(tension: f64, length: f64, linear_density: f64) -> f64 {
    (1.0 / (2.0 * length)) * (tension / linear_density).sqrt()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    // -----------------------------------------------------------------------
    // Helper math
    // -----------------------------------------------------------------------

    #[test]
    fn test_len3_unit_x() {
        assert!(approx_eq(len3([1.0, 0.0, 0.0]), 1.0, EPS));
    }

    #[test]
    fn test_normalize3_unit() {
        let v = normalize3([3.0, 0.0, 0.0]);
        assert!(approx_eq(v[0], 1.0, EPS));
        assert!(approx_eq(v[1], 0.0, EPS));
        assert!(approx_eq(v[2], 0.0, EPS));
    }

    #[test]
    fn test_normalize3_zero_returns_zero() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert!(approx_eq(len3(v), 0.0, EPS));
    }

    #[test]
    fn test_dot3() {
        assert!(approx_eq(dot3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]), 32.0, EPS));
    }

    // -----------------------------------------------------------------------
    // CableNode
    // -----------------------------------------------------------------------

    #[test]
    fn test_cable_node_kinetic_energy() {
        let mut n = CableNode::new([0.0; 3], 2.0);
        n.velocity = [3.0, 0.0, 0.0];
        // KE = 0.5 * 2 * 9 = 9 J
        assert!(approx_eq(n.kinetic_energy(), 9.0, EPS));
    }

    #[test]
    fn test_fixed_node_zero_kinetic_energy() {
        let mut n = CableNode::fixed([0.0; 3]);
        n.velocity = [10.0, 10.0, 10.0];
        assert!(approx_eq(n.kinetic_energy(), 0.0, EPS));
    }

    // -----------------------------------------------------------------------
    // Cable force
    // -----------------------------------------------------------------------

    #[test]
    fn test_cable_force_at_rest_length_is_pretension() {
        let cable = Cable::new(0, 1, 1.0, 1000.0, 0.0, 10.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [1.0, 0.0, 0.0]; // exactly at rest length
        let f = cable.compute_force(pa, [0.0; 3], pb, [0.0; 3]);
        // Force magnitude should equal pretension.
        assert!(
            approx_eq(len3(f), 10.0, 1e-8),
            "Force at rest length should equal pretension, got {f:?}"
        );
    }

    #[test]
    fn test_cable_force_zero_when_slack() {
        // Cable with zero pretension, nodes closer than rest length.
        let cable = Cable::new(0, 1, 2.0, 1000.0, 0.0, 0.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [0.5, 0.0, 0.0]; // compressed → slack
        let f = cable.compute_force(pa, [0.0; 3], pb, [0.0; 3]);
        assert!(
            approx_eq(len3(f), 0.0, EPS),
            "Slack cable should produce zero force, got {f:?}"
        );
    }

    #[test]
    fn test_cable_force_increases_with_extension() {
        let cable = Cable::new(0, 1, 1.0, 1000.0, 0.0, 0.0);
        let pa = [0.0, 0.0, 0.0];
        let f1 = cable.compute_force(pa, [0.0; 3], [1.5, 0.0, 0.0], [0.0; 3]);
        let f2 = cable.compute_force(pa, [0.0; 3], [2.0, 0.0, 0.0], [0.0; 3]);
        assert!(len3(f2) > len3(f1), "Force must increase with extension");
    }

    #[test]
    fn test_cable_force_direction() {
        // Cable stretched in +Y direction: force on A must point in +Y.
        let cable = Cable::new(0, 1, 1.0, 1000.0, 0.0, 0.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [0.0, 2.0, 0.0];
        let f = cable.compute_force(pa, [0.0; 3], pb, [0.0; 3]);
        assert!(
            f[1] > 0.0,
            "Force on A should point toward B (+Y), got {f:?}"
        );
    }

    #[test]
    fn test_cable_damping_modifies_force() {
        // With positive damping and nodes moving apart (A away from B),
        // the relative velocity along the cable is negative → damping reduces force.
        // With nodes moving together, damping increases force.
        let cable_no_damp = Cable::new(0, 1, 0.5, 1000.0, 0.0, 0.0);
        let cable_damp = Cable::new(0, 1, 0.5, 1000.0, 50.0, 0.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [1.0, 0.0, 0.0];
        // A is moving in -X (away from B in terms of relative velocity).
        let va = [-1.0, 0.0, 0.0];
        let vb = [0.0, 0.0, 0.0];
        let f_no = len3(cable_no_damp.compute_force(pa, va, pb, vb));
        let f_d = len3(cable_damp.compute_force(pa, va, pb, vb));
        // rel_vel = vb - va = (1, 0, 0); projected onto dir_ab = (1,0,0) → +1 → increases force.
        assert!(
            f_d >= f_no,
            "Positive rel_vel component should increase total damped force"
        );
    }

    #[test]
    fn test_cable_potential_energy_zero_at_rest() {
        let cable = Cable::new(0, 1, 1.0, 1000.0, 0.0, 0.0);
        let pa = [0.0, 0.0, 0.0];
        let pb = [1.0, 0.0, 0.0];
        assert!(approx_eq(cable.potential_energy(pa, pb), 0.0, EPS));
    }

    #[test]
    fn test_cable_potential_energy_quadratic() {
        let k = 1000.0;
        let cable = Cable::new(0, 1, 1.0, k, 0.0, 0.0);
        let ext = 0.1;
        let pa = [0.0; 3];
        let pb = [1.0 + ext, 0.0, 0.0];
        let pe = cable.potential_energy(pa, pb);
        // PE = 0.5 * k * ext² = 0.5 * 1000 * 0.01 = 5 J
        assert!(approx_eq(pe, 5.0, 1e-8));
    }

    // -----------------------------------------------------------------------
    // CableNet
    // -----------------------------------------------------------------------

    #[test]
    fn test_cable_net_single_cable_step() {
        let mut net = CableNet::new();
        let a = net.add_node(CableNode::fixed([0.0, 0.0, 0.0]));
        let b = net.add_node(CableNode::new([2.0, 0.0, 0.0], 1.0));
        net.add_cable(Cable::new(a, b, 1.0, 1000.0, 0.0, 0.0));
        // Node B should move toward A under cable tension.
        let x0 = net.nodes[b].position[0];
        net.step_with_gravity(0.01, [0.0; 3]);
        let x1 = net.nodes[b].position[0];
        assert!(x1 < x0, "Free node should be pulled toward fixed node");
    }

    #[test]
    fn test_cable_net_potential_energy_decreases() {
        let mut net = CableNet::new();
        let a = net.add_node(CableNode::fixed([0.0, 0.0, 0.0]));
        let b = net.add_node(CableNode::new([3.0, 0.0, 0.0], 1.0));
        net.add_cable(Cable::new(a, b, 1.0, 500.0, 10.0, 0.0));
        let pe0 = net.total_potential_energy();
        for _ in 0..50 {
            net.step_with_gravity(0.005, [0.0; 3]);
        }
        let pe1 = net.total_potential_energy();
        assert!(
            pe1 < pe0,
            "Potential energy should decrease as cable contracts"
        );
    }

    #[test]
    fn test_cable_net_gravity_falls() {
        let mut net = CableNet::new();
        let b = net.add_node(CableNode::new([0.0, 10.0, 0.0], 1.0));
        // No cable — free fall under gravity for 2 seconds.
        let dt = 0.01;
        for _ in 0..200 {
            net.step_with_gravity(dt, [0.0, -9.81, 0.0]);
        }
        // After ~2 s: y ≈ 10 - 0.5*9.81*4 ≈ -9.6 (explicit Euler slightly larger).
        assert!(
            net.nodes[b].position[1] < 0.0,
            "Node should have fallen well below y=0 under 2 s of gravity, got y={}",
            net.nodes[b].position[1]
        );
    }

    #[test]
    fn test_cable_net_fixed_node_does_not_move() {
        let mut net = CableNet::new();
        let a = net.add_node(CableNode::fixed([5.0, 5.0, 5.0]));
        net.add_cable(Cable::new(a, a, 0.0, 1000.0, 0.0, 100.0));
        for _ in 0..100 {
            net.step_with_gravity(0.01, [0.0, -9.81, 0.0]);
        }
        let p = net.nodes[a].position;
        assert!(
            approx_eq(p[0], 5.0, EPS) && approx_eq(p[1], 5.0, EPS) && approx_eq(p[2], 5.0, EPS),
            "Fixed node should not move"
        );
    }

    // -----------------------------------------------------------------------
    // TensegrityStrut
    // -----------------------------------------------------------------------

    #[test]
    fn test_strut_force_zero_when_stretched() {
        // Strut is compression-only: no force when extended beyond rest length.
        let strut = TensegrityStrut::new(0, 1, 1.0, 5000.0);
        let f = strut.compute_force([0.0; 3], [2.0, 0.0, 0.0]);
        assert!(
            approx_eq(len3(f), 0.0, EPS),
            "Strut should produce zero force in tension"
        );
    }

    #[test]
    fn test_strut_force_nonzero_when_compressed() {
        let strut = TensegrityStrut::new(0, 1, 2.0, 5000.0);
        let f = strut.compute_force([0.0; 3], [1.0, 0.0, 0.0]); // length=1, rest=2
        assert!(len3(f) > 0.0, "Compressed strut should push nodes apart");
    }

    #[test]
    fn test_strut_force_direction_pushes_apart() {
        let strut = TensegrityStrut::new(0, 1, 2.0, 5000.0);
        // A at origin, B at (1, 0, 0).  Compressed → force on A should point in -X.
        let f = strut.compute_force([0.0; 3], [1.0, 0.0, 0.0]);
        assert!(f[0] < 0.0, "Strut force on A should push in -X direction");
    }

    #[test]
    fn test_strut_potential_energy_zero_when_extended() {
        let strut = TensegrityStrut::new(0, 1, 1.0, 5000.0);
        let pe = strut.potential_energy([0.0; 3], [2.0, 0.0, 0.0]);
        assert!(approx_eq(pe, 0.0, EPS));
    }

    #[test]
    fn test_strut_potential_energy_quadratic() {
        let k = 5000.0;
        let strut = TensegrityStrut::new(0, 1, 2.0, k);
        let _comp = 0.5; // compression = 0.5 m
        let pe = strut.potential_energy([0.0; 3], [1.5, 0.0, 0.0]);
        // PE = 0.5 * k * _comp² = 0.5 * 5000 * 0.25 = 625 J
        assert!(approx_eq(pe, 625.0, 1e-6));
    }

    // -----------------------------------------------------------------------
    // TensegrityStructure
    // -----------------------------------------------------------------------

    #[test]
    fn test_tensegrity_structure_find_equilibrium_runs() {
        // Simple two-node tensegrity: one cable + one strut opposing each other.
        let mut ts = TensegrityStructure::new();
        let a = ts.add_node(CableNode::fixed([0.0, 0.0, 0.0]));
        let b = ts.add_node(CableNode::new([1.0, 0.0, 0.0], 1.0));
        // Cable wants rest length 0.5 (pull B toward A).
        ts.add_cable(Cable::new(a, b, 0.5, 1000.0, 0.0, 0.0));
        // Strut wants rest length 1.5 (push B away from A).
        ts.add_strut(TensegrityStrut::new(a, b, 1.5, 1000.0));
        // Equilibrium at L = 1.0 (between 0.5 and 1.5 based on equal stiffness).
        ts.find_equilibrium(5000, 1e-8);
        // After equilibrium, net forces should be small.
        let forces = ts.net_forces([0.0; 3]);
        let max_f = forces.iter().map(|f| len3(*f)).fold(0.0_f64, f64::max);
        assert!(
            max_f < 0.5,
            "Equilibrium residual force {max_f:.3e} too large"
        );
    }

    #[test]
    fn test_tensegrity_potential_energy_nonneg() {
        let mut ts = TensegrityStructure::new();
        let a = ts.add_node(CableNode::fixed([0.0; 3]));
        let b = ts.add_node(CableNode::new([1.5, 0.0, 0.0], 1.0));
        ts.add_cable(Cable::new(a, b, 1.0, 500.0, 0.0, 0.0));
        ts.add_strut(TensegrityStrut::new(a, b, 2.0, 500.0));
        assert!(ts.total_potential_energy() >= 0.0);
    }

    // -----------------------------------------------------------------------
    // FormFinding
    // -----------------------------------------------------------------------

    #[test]
    fn test_form_finding_converges() {
        // A single free node connected by cable to a fixed node.
        // The cable will pull the free node to rest_length.
        let mut net = CableNet::new();
        let a = net.add_node(CableNode::fixed([0.0, 0.0, 0.0]));
        let b = net.add_node(CableNode::new([5.0, 0.0, 0.0], 1.0));
        net.add_cable(Cable::new(a, b, 1.0, 1000.0, 5.0, 0.0));

        let ff = FormFinding::new(10_000, 1e-3);
        let iters = ff.solve(&mut net);
        let final_len = len3(sub3(net.nodes[b].position, net.nodes[a].position));
        assert!(
            (final_len - 1.0).abs() < 0.05,
            "Form-finding did not converge to rest length; got {final_len:.4} in {iters} iters"
        );
    }

    #[test]
    fn test_form_finding_at_rest_returns_early() {
        // A cable already at rest length — should converge immediately.
        let mut net = CableNet::new();
        let a = net.add_node(CableNode::fixed([0.0, 0.0, 0.0]));
        let b = net.add_node(CableNode::new([1.0, 0.0, 0.0], 1.0));
        net.add_cable(Cable::new(a, b, 1.0, 1000.0, 0.0, 0.0));

        let ff = FormFinding::new(1000, 1e-6);
        let iters = ff.solve(&mut net);
        assert!(iters < 1000, "Should converge early when already at rest");
    }

    // -----------------------------------------------------------------------
    // cable_sag_parabolic
    // -----------------------------------------------------------------------

    #[test]
    fn test_cable_sag_basic() {
        // 10 m span, 1 kN/m load, 50 kN tension → sag = 1000*100/(8*50000) = 0.25 m
        let sag = cable_sag_parabolic(10.0, 1000.0, 50_000.0);
        assert!(
            approx_eq(sag, 0.25, 1e-10),
            "Sag should be 0.25 m, got {sag}"
        );
    }

    #[test]
    fn test_cable_sag_increases_with_load() {
        let sag1 = cable_sag_parabolic(10.0, 500.0, 10_000.0);
        let sag2 = cable_sag_parabolic(10.0, 1000.0, 10_000.0);
        assert!(sag2 > sag1, "Heavier load should produce larger sag");
    }

    #[test]
    fn test_cable_sag_decreases_with_tension() {
        let sag1 = cable_sag_parabolic(10.0, 1000.0, 10_000.0);
        let sag2 = cable_sag_parabolic(10.0, 1000.0, 20_000.0);
        assert!(sag2 < sag1, "Higher tension should reduce sag");
    }

    #[test]
    fn test_cable_sag_scales_with_span_squared() {
        let w = 1000.0;
        let t = 50_000.0;
        let s1 = cable_sag_parabolic(10.0, w, t);
        let s2 = cable_sag_parabolic(20.0, w, t);
        // sag ∝ L² → doubling span quadruples sag.
        assert!(approx_eq(s2, 4.0 * s1, 1e-8));
    }

    // -----------------------------------------------------------------------
    // taut_string_frequency
    // -----------------------------------------------------------------------

    #[test]
    fn test_taut_string_frequency_basic() {
        // Guitar string: T = 70 N, L = 0.65 m, μ = 0.00038 kg/m
        // f = (1 / 1.3) * sqrt(70 / 0.00038) ≈ 328 Hz  (roughly E4)
        let f = taut_string_frequency(70.0, 0.65, 3.8e-4);
        assert!(f > 200.0 && f < 500.0, "Expected ~328 Hz, got {f:.1} Hz");
    }

    #[test]
    fn test_taut_string_frequency_increases_with_tension() {
        let f1 = taut_string_frequency(100.0, 1.0, 0.01);
        let f2 = taut_string_frequency(400.0, 1.0, 0.01);
        // Doubling tension increases f by sqrt(4) = 2.
        assert!(approx_eq(f2 / f1, 2.0, 1e-8));
    }

    #[test]
    fn test_taut_string_frequency_halves_with_double_length() {
        let f1 = taut_string_frequency(100.0, 1.0, 0.01);
        let f2 = taut_string_frequency(100.0, 2.0, 0.01);
        assert!(approx_eq(f2, 0.5 * f1, 1e-10));
    }

    #[test]
    fn test_taut_string_frequency_positive() {
        let f = taut_string_frequency(50.0, 0.5, 0.005);
        assert!(f > 0.0, "Frequency must be positive");
    }

    #[test]
    fn test_taut_string_frequency_decreases_with_density() {
        let f1 = taut_string_frequency(100.0, 1.0, 0.01);
        let f2 = taut_string_frequency(100.0, 1.0, 0.04); // 4x heavier
        assert!(approx_eq(f2, 0.5 * f1, 1e-10));
    }
}
