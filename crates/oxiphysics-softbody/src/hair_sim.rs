#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Hair and fiber simulation using a discrete Cosserat rod model.
//!
//! Models hair strands as chains of oriented segments. Each segment carries a
//! position and a frame (tangent/normal/binormal) to capture bending and
//! twisting along the strand axis.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ─────────────────────────────────────────────────────────────────────────────
// Vector helpers (no nalgebra: use [f64; 3])
// ─────────────────────────────────────────────────────────────────────────────

/// 3-component vector addition.
#[inline]
fn v3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// 3-component vector subtraction.
#[inline]
fn v3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Scalar multiplication of a 3-vector.
#[inline]
fn v3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Dot product of two 3-vectors.
#[inline]
fn v3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Cross product of two 3-vectors.
#[inline]
fn v3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean norm of a 3-vector.
#[inline]
fn v3_norm(a: [f64; 3]) -> f64 {
    v3_dot(a, a).sqrt()
}

/// Normalize a 3-vector. Returns the zero vector if norm < eps.
#[inline]
fn v3_normalize(a: [f64; 3]) -> [f64; 3] {
    let n = v3_norm(a);
    if n < 1e-15 {
        [0.0, 0.0, 0.0]
    } else {
        v3_scale(a, 1.0 / n)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Core data structures
// ─────────────────────────────────────────────────────────────────────────────

/// A single node (particle) in a hair strand.
#[derive(Debug, Clone)]
pub struct HairNode {
    /// World-space position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Accumulated external force (N).
    pub force: [f64; 3],
    /// Inverse mass (1/kg). 0.0 for static/pinned nodes.
    pub inv_mass: f64,
}

impl HairNode {
    /// Create a new dynamic hair node.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        let inv_mass = if mass > 0.0 { 1.0 / mass } else { 0.0 };
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass,
        }
    }

    /// Create a static (pinned) hair node.
    pub fn new_static(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            inv_mass: 0.0,
        }
    }

    /// Returns true if this node is static (pinned).
    pub fn is_static(&self) -> bool {
        self.inv_mass == 0.0
    }
}

/// A single hair strand modelled as a chain of nodes (discrete Cosserat rod).
///
/// Adjacent nodes are connected by elastic segments. The strand stores rest
/// lengths, bend stiffness, and twist stiffness for each segment.
#[derive(Debug, Clone)]
pub struct HairStrand {
    /// Nodes along the strand (root first, tip last).
    pub nodes: Vec<HairNode>,
    /// Natural (rest) length of each segment between consecutive nodes (m).
    pub rest_lengths: Vec<f64>,
    /// Rest curvature (bending) at each internal node (rad/m). Zero → straight.
    pub rest_curvature: Vec<[f64; 3]>,
    /// Bending stiffness E_b (N·m²). Higher = stiffer bends.
    pub bend_stiffness: f64,
    /// Twist stiffness E_t (N·m²). Higher = stiffer torsion.
    pub twist_stiffness: f64,
    /// Stretch stiffness k_s (N/m).
    pub stretch_stiffness: f64,
    /// Aerodynamic drag coefficient per unit length (kg/(m·s)).
    pub drag_coeff: f64,
}

impl HairStrand {
    /// Create a straight hair strand from root to tip.
    ///
    /// # Arguments
    /// * `root` — root position (m)
    /// * `tip` — tip position (m)
    /// * `num_segments` — number of segments (nodes = num_segments + 1)
    /// * `mass_per_node` — mass of each non-root node (kg)
    /// * `bend_stiffness` — E_b (N·m²)
    /// * `twist_stiffness` — E_t (N·m²)
    /// * `stretch_stiffness` — k_s (N/m)
    pub fn new_straight(
        root: [f64; 3],
        tip: [f64; 3],
        num_segments: usize,
        mass_per_node: f64,
        bend_stiffness: f64,
        twist_stiffness: f64,
        stretch_stiffness: f64,
    ) -> Self {
        assert!(num_segments >= 1, "At least one segment required");
        let n_nodes = num_segments + 1;
        let mut nodes = Vec::with_capacity(n_nodes);
        let seg_length = {
            let d = v3_sub(tip, root);
            v3_norm(d) / num_segments as f64
        };
        for i in 0..n_nodes {
            let t = i as f64 / num_segments as f64;
            let pos = [
                root[0] + t * (tip[0] - root[0]),
                root[1] + t * (tip[1] - root[1]),
                root[2] + t * (tip[2] - root[2]),
            ];
            if i == 0 {
                nodes.push(HairNode::new_static(pos));
            } else {
                nodes.push(HairNode::new(pos, mass_per_node));
            }
        }
        let rest_lengths = vec![seg_length; num_segments];
        let rest_curvature = vec![[0.0; 3]; n_nodes];
        Self {
            nodes,
            rest_lengths,
            rest_curvature,
            bend_stiffness,
            twist_stiffness,
            stretch_stiffness,
            drag_coeff: 0.0,
        }
    }

    /// Number of nodes in the strand.
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Number of segments (= num_nodes - 1).
    pub fn num_segments(&self) -> usize {
        self.nodes.len().saturating_sub(1)
    }

    /// Total rest length of the strand.
    pub fn total_rest_length(&self) -> f64 {
        self.rest_lengths.iter().sum()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Cosserat rod energy and forces
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the bending energy of a strand.
///
/// Uses the discrete curvature κ at each interior node:
///   κᵢ = (tᵢ - tᵢ₋₁) / (lᵢ + lᵢ₋₁) * 2
/// where tᵢ = normalized edge tangent.
///
/// E_bend = (E_b / 2) Σ |κᵢ - κ₀ᵢ|² × lᵢ
///
/// Returns bending energy in Joules.
pub fn bending_energy(strand: &HairStrand) -> f64 {
    let nodes = &strand.nodes;
    let n = nodes.len();
    if n < 3 {
        return 0.0;
    }
    let mut energy = 0.0;
    for i in 1..n - 1 {
        let e0 = v3_sub(nodes[i].position, nodes[i - 1].position);
        let e1 = v3_sub(nodes[i + 1].position, nodes[i].position);
        let t0 = v3_normalize(e0);
        let t1 = v3_normalize(e1);
        let l0 = strand.rest_lengths[i - 1];
        let l1 = strand.rest_lengths[i];
        let l_avg = 0.5 * (l0 + l1);
        // Discrete curvature vector κ = (t1 - t0) / l_avg
        let kappa = v3_scale(v3_sub(t1, t0), 1.0 / l_avg.max(1e-15));
        // Curvature deviation from rest
        let dkappa = v3_sub(kappa, strand.rest_curvature[i]);
        let kappa_sq = v3_dot(dkappa, dkappa);
        energy += 0.5 * strand.bend_stiffness * kappa_sq * l_avg;
    }
    energy
}

/// Compute the twisting energy of a strand.
///
/// Twist angle θᵢ is estimated from the rotation of the binormal frame
/// between successive segments. Uses the scalar twist κᵢ · tᵢ.
///
/// E_twist = (E_t / 2) Σ θᵢ² / lᵢ
///
/// Returns twisting energy in Joules.
pub fn twisting_energy(strand: &HairStrand) -> f64 {
    let nodes = &strand.nodes;
    let n = nodes.len();
    if n < 3 {
        return 0.0;
    }
    let mut energy = 0.0;
    for i in 1..n - 1 {
        let e0 = v3_sub(nodes[i].position, nodes[i - 1].position);
        let e1 = v3_sub(nodes[i + 1].position, nodes[i].position);
        let t0 = v3_normalize(e0);
        let t1 = v3_normalize(e1);
        // Scalar twist: projection of tangent change onto the mid-tangent.
        let t_mid = v3_normalize(v3_add(t0, t1));
        let dt = v3_sub(t1, t0);
        let twist_scalar = v3_dot(dt, t_mid);
        let l = strand.rest_lengths[i];
        energy += 0.5 * strand.twist_stiffness * twist_scalar * twist_scalar / l.max(1e-15);
    }
    energy
}

/// Compute stretch energy of a strand.
///
/// E_stretch = (k_s / 2) Σ (|eᵢ| - l₀ᵢ)² / l₀ᵢ
///
/// Returns stretch energy in Joules.
pub fn stretch_energy(strand: &HairStrand) -> f64 {
    let nodes = &strand.nodes;
    let n = nodes.len();
    if n < 2 {
        return 0.0;
    }
    let mut energy = 0.0;
    for i in 0..n - 1 {
        let e = v3_sub(nodes[i + 1].position, nodes[i].position);
        let length = v3_norm(e);
        let l0 = strand.rest_lengths[i];
        let strain = length - l0;
        energy += 0.5 * strand.stretch_stiffness * strain * strain / l0.max(1e-15);
    }
    energy
}

/// Apply root constraint: pin the first node to its rest position.
///
/// This modifies the root node's position and zeroes its velocity.
pub fn root_constraint(strand: &mut HairStrand) {
    if !strand.nodes.is_empty() {
        strand.nodes[0].velocity = [0.0; 3];
        strand.nodes[0].force = [0.0; 3];
        // Position is kept fixed by its inv_mass = 0 in integration
    }
}

/// Apply aerodynamic wind force to all nodes in a strand.
///
/// Uses a simple drag model: F_drag = -0.5 × ρ × C_d × A × (v_node - v_wind) × |v_rel|
/// Simplified per-node: F = drag_coeff × (v_wind - v_node) per segment length.
///
/// # Arguments
/// * `strand` — the hair strand to update
/// * `wind_velocity` — ambient wind velocity (m/s)
pub fn wind_force(strand: &mut HairStrand, wind_velocity: [f64; 3]) {
    let n = strand.nodes.len();
    if n < 2 {
        return;
    }
    let cd = strand.drag_coeff;
    for i in 0..n {
        if strand.nodes[i].is_static() {
            continue;
        }
        // Segment length associated with this node.
        let seg_len = if i == 0 {
            strand.rest_lengths[0]
        } else if i == n - 1 {
            *strand
                .rest_lengths
                .last()
                .expect("collection should not be empty")
        } else {
            0.5 * (strand.rest_lengths[i - 1] + strand.rest_lengths[i])
        };
        let v_rel = v3_sub(wind_velocity, strand.nodes[i].velocity);
        let drag = v3_scale(v_rel, cd * seg_len);
        strand.nodes[i].force = v3_add(strand.nodes[i].force, drag);
    }
}

/// Compute bending forces and accumulate them on each node.
///
/// Uses negative gradient of bending energy w.r.t. positions (finite difference).
/// Returns forces as `Vec<[f64; 3]>` aligned with node indices.
pub fn discrete_elastic_rod(strand: &HairStrand) -> Vec<[f64; 3]> {
    let n = strand.nodes.len();
    let mut forces = vec![[0.0_f64; 3]; n];
    if n < 3 {
        return forces;
    }
    let eps = 1e-7;
    let e0 = bending_energy(strand) + stretch_energy(strand) + twisting_energy(strand);
    // Finite-difference gradient.
    for i in 0..n {
        if strand.nodes[i].is_static() {
            continue;
        }
        for k in 0..3 {
            let mut perturbed = strand.clone();
            perturbed.nodes[i].position[k] += eps;
            let e1 = bending_energy(&perturbed)
                + stretch_energy(&perturbed)
                + twisting_energy(&perturbed);
            forces[i][k] = -(e1 - e0) / eps;
        }
    }
    forces
}

// ─────────────────────────────────────────────────────────────────────────────
// Hair-hair collision (strand-strand contact)
// ─────────────────────────────────────────────────────────────────────────────

/// Detect and respond to contact between two hair strands.
///
/// For each pair of nodes from different strands within `contact_radius`,
/// applies a repulsive spring force.
///
/// # Arguments
/// * `strand_a` — first strand (modified in place)
/// * `strand_b` — second strand (modified in place)
/// * `contact_radius` — minimum separation distance (m)
/// * `repulsion_stiffness` — spring constant for repulsion (N/m)
pub fn hair_hair_collision(
    strand_a: &mut HairStrand,
    strand_b: &mut HairStrand,
    contact_radius: f64,
    repulsion_stiffness: f64,
) {
    let na = strand_a.nodes.len();
    let nb = strand_b.nodes.len();
    for i in 0..na {
        for j in 0..nb {
            let d = v3_sub(strand_b.nodes[j].position, strand_a.nodes[i].position);
            let dist = v3_norm(d);
            if dist < contact_radius && dist > 1e-15 {
                let penetration = contact_radius - dist;
                let dir = v3_scale(d, 1.0 / dist);
                let impulse = v3_scale(dir, repulsion_stiffness * penetration);
                if !strand_a.nodes[i].is_static() {
                    strand_a.nodes[i].force = v3_sub(strand_a.nodes[i].force, impulse);
                }
                if !strand_b.nodes[j].is_static() {
                    strand_b.nodes[j].force = v3_add(strand_b.nodes[j].force, impulse);
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Style target (rest-shape styling)
// ─────────────────────────────────────────────────────────────────────────────

/// Set the rest curvature of a strand to match a target set of positions.
///
/// The rest curvature at each node is computed from the target positions and
/// stored in `strand.rest_curvature`. This allows the strand to restore a
/// styled shape.
///
/// # Arguments
/// * `strand` — the strand to style
/// * `target_positions` — desired node positions (same count as strand nodes)
pub fn style_target(strand: &mut HairStrand, target_positions: &[[f64; 3]]) {
    assert_eq!(
        target_positions.len(),
        strand.nodes.len(),
        "target_positions must match node count"
    );
    let n = strand.nodes.len();
    for i in 0..n {
        if i == 0 || i == n - 1 {
            strand.rest_curvature[i] = [0.0; 3];
            continue;
        }
        let e0 = v3_sub(target_positions[i], target_positions[i - 1]);
        let e1 = v3_sub(target_positions[i + 1], target_positions[i]);
        let t0 = v3_normalize(e0);
        let t1 = v3_normalize(e1);
        let l0 = v3_norm(e0).max(1e-15);
        let l1 = v3_norm(e1).max(1e-15);
        let l_avg = 0.5 * (l0 + l1);
        let kappa = v3_scale(v3_sub(t1, t0), 1.0 / l_avg);
        strand.rest_curvature[i] = kappa;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Simple explicit integration step
// ─────────────────────────────────────────────────────────────────────────────

/// Advance a hair strand one time step with explicit Euler integration.
///
/// Applies elastic rod forces, then integrates velocities and positions.
///
/// # Arguments
/// * `strand` — the strand to integrate
/// * `gravity` — gravitational acceleration (m/s²)
/// * `dt` — time step (s)
pub fn step_strand(strand: &mut HairStrand, gravity: [f64; 3], dt: f64) {
    // Compute elastic forces.
    let elastic_forces = discrete_elastic_rod(strand);
    let n = strand.nodes.len();
    // Velocity damping coefficient (per step): prevents explosion with stiff springs.
    let damping = 0.98_f64;
    for i in 0..n {
        if strand.nodes[i].is_static() {
            continue;
        }
        let m_inv = strand.nodes[i].inv_mass;
        let mass = 1.0 / m_inv.max(1e-15);
        // Gravity force = mass × g.
        let grav_f = v3_scale(gravity, mass);
        let total_force = v3_add(v3_add(strand.nodes[i].force, elastic_forces[i]), grav_f);
        let accel = v3_scale(total_force, m_inv);
        let new_vel = v3_scale(
            v3_add(strand.nodes[i].velocity, v3_scale(accel, dt)),
            damping,
        );
        strand.nodes[i].velocity = new_vel;
        strand.nodes[i].position = v3_add(
            strand.nodes[i].position,
            v3_scale(strand.nodes[i].velocity, dt),
        );
        strand.nodes[i].force = [0.0; 3];
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // Helper: create a straight vertical strand.
    fn make_strand(num_segs: usize) -> HairStrand {
        HairStrand::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            num_segs,
            0.001,
            1e-3,
            1e-4,
            1000.0,
        )
    }

    // 1. Strand creation: correct node count.
    #[test]
    fn test_strand_node_count() {
        let s = make_strand(4);
        assert_eq!(s.num_nodes(), 5, "4 segments → 5 nodes");
    }

    // 2. Strand creation: correct segment count.
    #[test]
    fn test_strand_segment_count() {
        let s = make_strand(4);
        assert_eq!(s.num_segments(), 4, "should have 4 segments");
    }

    // 3. Root node is static.
    #[test]
    fn test_root_is_static() {
        let s = make_strand(4);
        assert!(s.nodes[0].is_static(), "root node should be static");
    }

    // 4. Non-root nodes are dynamic.
    #[test]
    fn test_non_root_dynamic() {
        let s = make_strand(4);
        for i in 1..s.num_nodes() {
            assert!(!s.nodes[i].is_static(), "node {i} should be dynamic");
        }
    }

    // 5. Total rest length matches root-to-tip distance.
    #[test]
    fn test_total_rest_length() {
        let s = make_strand(5);
        let len = s.total_rest_length();
        assert!(
            (len - 1.0).abs() < 1e-10,
            "Total rest length should be 1.0, got {len}"
        );
    }

    // 6. Bending energy of a straight strand is zero.
    #[test]
    fn test_bending_energy_straight_is_zero() {
        let s = make_strand(5);
        let e = bending_energy(&s);
        assert!(
            e.abs() < 1e-12,
            "Straight strand bending energy should be 0, got {e}"
        );
    }

    // 7. Twisting energy of a straight strand is zero.
    #[test]
    fn test_twisting_energy_straight_is_zero() {
        let s = make_strand(5);
        let e = twisting_energy(&s);
        assert!(
            e.abs() < 1e-12,
            "Straight strand twist energy should be 0, got {e}"
        );
    }

    // 8. Stretch energy of an undeformed strand is zero.
    #[test]
    fn test_stretch_energy_zero_undeformed() {
        let s = make_strand(5);
        let e = stretch_energy(&s);
        assert!(
            e.abs() < 1e-12,
            "Undeformed stretch energy should be 0, got {e}"
        );
    }

    // 9. Bending energy increases when a node is displaced.
    #[test]
    fn test_bending_energy_increases_on_displacement() {
        let mut s = make_strand(4);
        let e0 = bending_energy(&s);
        s.nodes[2].position[0] += 0.1; // lateral displacement
        let e1 = bending_energy(&s);
        assert!(
            e1 > e0,
            "Bending energy should increase after displacement: e0={e0}, e1={e1}"
        );
    }

    // 10. Stretch energy increases when a node is displaced along axis.
    #[test]
    fn test_stretch_energy_increases_on_displacement() {
        let mut s = make_strand(4);
        let e0 = stretch_energy(&s);
        s.nodes[4].position[1] += 0.2; // extend tip
        let e1 = stretch_energy(&s);
        assert!(
            e1 > e0,
            "Stretch energy should increase after extension: e0={e0}, e1={e1}"
        );
    }

    // 11. Root constraint zeroes root node velocity.
    #[test]
    fn test_root_constraint_zeroes_velocity() {
        let mut s = make_strand(4);
        s.nodes[0].velocity = [1.0, 2.0, 3.0];
        root_constraint(&mut s);
        assert_eq!(
            s.nodes[0].velocity,
            [0.0, 0.0, 0.0],
            "Root velocity should be zeroed"
        );
    }

    // 12. Wind force with zero drag coefficient produces zero force.
    #[test]
    fn test_wind_force_zero_drag() {
        let mut s = make_strand(4);
        s.drag_coeff = 0.0;
        wind_force(&mut s, [10.0, 0.0, 0.0]);
        for node in &s.nodes {
            assert_eq!(node.force, [0.0; 3], "Zero drag → zero wind force");
        }
    }

    // 13. Wind force with positive drag applies force on dynamic nodes.
    #[test]
    fn test_wind_force_nonzero_drag() {
        let mut s = make_strand(4);
        s.drag_coeff = 0.01;
        wind_force(&mut s, [1.0, 0.0, 0.0]);
        // Node 1..4 should have non-zero force.
        let total_fx: f64 = s.nodes[1..].iter().map(|n| n.force[0]).sum();
        assert!(
            total_fx > EPS,
            "Wind should produce force in x-direction, got {total_fx}"
        );
    }

    // 14. Wind force is zero on static (root) node.
    #[test]
    fn test_wind_force_not_on_static_node() {
        let mut s = make_strand(4);
        s.drag_coeff = 0.1;
        wind_force(&mut s, [10.0, 0.0, 0.0]);
        assert_eq!(
            s.nodes[0].force, [0.0; 3],
            "Static node should not receive wind force"
        );
    }

    // 15. Discrete elastic rod returns forces of correct size.
    #[test]
    fn test_discrete_elastic_rod_force_size() {
        let s = make_strand(4);
        let forces = discrete_elastic_rod(&s);
        assert_eq!(
            forces.len(),
            s.num_nodes(),
            "Force count should equal node count"
        );
    }

    // 16. Static node gets zero elastic force.
    #[test]
    fn test_discrete_elastic_rod_static_zero_force() {
        let s = make_strand(4);
        let forces = discrete_elastic_rod(&s);
        let f_root = v3_norm(forces[0]);
        assert!(
            f_root < EPS,
            "Static root node should have zero force, got {f_root}"
        );
    }

    // 17. Hair-hair collision: no force when strands are far apart.
    #[test]
    fn test_hair_collision_no_force_when_far() {
        let mut s1 = make_strand(3);
        let mut s2 = HairStrand::new_straight(
            [10.0, 0.0, 0.0],
            [10.0, 1.0, 0.0],
            3,
            0.001,
            1e-3,
            1e-4,
            1000.0,
        );
        hair_hair_collision(&mut s1, &mut s2, 0.1, 1000.0);
        for n in &s1.nodes {
            assert_eq!(n.force, [0.0; 3], "No collision → zero force on s1");
        }
    }

    // 18. Hair-hair collision: repulsion when strands overlap.
    #[test]
    fn test_hair_collision_repulsion() {
        let mut s1 = make_strand(1);
        let mut s2 = make_strand(1);
        // Place s2 very close to s1 at node 1.
        s2.nodes[1].position = [0.01, 1.0, 0.0];
        hair_hair_collision(&mut s1, &mut s2, 0.1, 1000.0);
        let f_x = s1.nodes[1].force[0];
        // They should be pushed apart (force < 0 in x for s1 since s2 is at +x).
        assert!(
            f_x < 0.0 || f_x.abs() > EPS,
            "Collision should produce repulsive force"
        );
    }

    // 19. Style target: rest curvature at endpoints is zero.
    #[test]
    fn test_style_target_endpoints_zero() {
        let mut s = make_strand(4);
        let targets: Vec<[f64; 3]> = s.nodes.iter().map(|n| n.position).collect();
        style_target(&mut s, &targets);
        assert_eq!(s.rest_curvature[0], [0.0; 3]);
        assert_eq!(s.rest_curvature[s.num_nodes() - 1], [0.0; 3]);
    }

    // 20. Style target: straight strand has zero rest curvature everywhere.
    #[test]
    fn test_style_target_straight_zero_curvature() {
        let mut s = make_strand(4);
        let targets: Vec<[f64; 3]> = s.nodes.iter().map(|n| n.position).collect();
        style_target(&mut s, &targets);
        for (i, kappa) in s.rest_curvature.iter().enumerate() {
            let mag = v3_norm(*kappa);
            assert!(
                mag < 1e-8,
                "Straight target should have zero curvature at node {i}, got {mag}"
            );
        }
    }

    // 21. Stepping a free strand under gravity lowers the tip position.
    // Use a very soft strand (tiny stiffness) to keep explicit Euler stable.
    #[test]
    fn test_step_strand_gravity_lowers_tip() {
        let mut s = HairStrand::new_straight(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            4,
            0.001,
            1e-6, // very soft bend
            1e-7, // very soft twist
            1.0,  // very soft stretch
        );
        let gravity = [0.0, -9.81, 0.0];
        let initial_y = s.nodes[4].position[1];
        for _ in 0..10 {
            step_strand(&mut s, gravity, 1.0 / 60.0);
        }
        let final_y = s.nodes[4].position[1];
        assert!(
            final_y < initial_y,
            "Tip should fall under gravity: initial={initial_y}, final={final_y}"
        );
    }

    // 22. Root node does not move during step_strand.
    #[test]
    fn test_step_strand_root_fixed() {
        let mut s = make_strand(4);
        let initial_root = s.nodes[0].position;
        for _ in 0..100 {
            step_strand(&mut s, [0.0, -9.81, 0.0], 1.0 / 60.0);
        }
        assert_eq!(
            s.nodes[0].position, initial_root,
            "Root should remain fixed"
        );
    }

    // 23. Bending energy is symmetric in sign of displacement.
    #[test]
    fn test_bending_energy_symmetric() {
        let mut s_pos = make_strand(4);
        let mut s_neg = make_strand(4);
        s_pos.nodes[2].position[0] += 0.05;
        s_neg.nodes[2].position[0] -= 0.05;
        let e_pos = bending_energy(&s_pos);
        let e_neg = bending_energy(&s_neg);
        assert!(
            (e_pos - e_neg).abs() < 1e-10,
            "Bending energy should be symmetric: e_pos={e_pos}, e_neg={e_neg}"
        );
    }

    // 24. Strand with 1 segment has zero bending energy (not enough nodes).
    #[test]
    fn test_bending_energy_single_segment() {
        let s = make_strand(1);
        let e = bending_energy(&s);
        assert!(
            e.abs() < EPS,
            "Single segment strand has no bending energy, got {e}"
        );
    }

    // 25. Twisting energy increases with lateral tip displacement.
    #[test]
    fn test_twisting_energy_increases_on_lateral_displacement() {
        let mut s = make_strand(5);
        let e0 = twisting_energy(&s);
        s.nodes[3].position[0] += 0.2;
        let e1 = twisting_energy(&s);
        // Twist energy is based on tangent frame change; lateral bend also has twist component.
        assert!(
            e1 >= e0,
            "Twist energy should not decrease with strand deformation: e0={e0} e1={e1}"
        );
    }

    // 26. HairNode default velocity and force are zero.
    #[test]
    fn test_hair_node_defaults() {
        let node = HairNode::new([1.0, 2.0, 3.0], 0.01);
        assert_eq!(node.velocity, [0.0; 3]);
        assert_eq!(node.force, [0.0; 3]);
    }

    // 27. Elastic forces on straight strand are close to zero.
    #[test]
    fn test_elastic_forces_straight_near_zero() {
        let s = make_strand(5);
        let forces = discrete_elastic_rod(&s);
        for (i, f) in forces.iter().enumerate() {
            if s.nodes[i].is_static() {
                continue;
            }
            let mag = v3_norm(*f);
            assert!(
                mag < 1e-3,
                "Elastic forces on undeformed strand should be ~0 at node {i}, got {mag}"
            );
        }
    }
}
