// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rope and cable simulation using Position-Based Dynamics (PBD).
//!
//! Provides a standalone rope model with catenary sag, natural frequency,
//! wind drag, and energy computations using plain `[f64; 3]` arrays.

// ---------------------------------------------------------------------------
// Vector helpers
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

// ---------------------------------------------------------------------------
// RopeNode
// ---------------------------------------------------------------------------

/// A single mass node in a rope chain.
#[derive(Debug, Clone)]
pub struct RopeNode {
    /// World-space position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Mass of the node (kg).
    pub mass: f64,
    /// When `true` the node is kinematically fixed.
    pub pinned: bool,
}

impl RopeNode {
    /// Create a new rope node at `position` with given `mass`.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            pinned: false,
        }
    }
}

// ---------------------------------------------------------------------------
// RopeSegment
// ---------------------------------------------------------------------------

/// A constraint segment connecting two adjacent rope nodes.
#[derive(Debug, Clone)]
pub struct RopeSegment {
    /// Index of the first node.
    pub i: usize,
    /// Index of the second node.
    pub j: usize,
    /// Rest (natural) length of the segment (m).
    pub rest_length: f64,
    /// Stretch stiffness coefficient (0–1).
    pub stiffness: f64,
    /// Velocity damping coefficient (s⁻¹).
    pub damping: f64,
}

// ---------------------------------------------------------------------------
// Rope
// ---------------------------------------------------------------------------

/// A chain of mass nodes connected by distance constraints.
#[derive(Debug, Clone)]
pub struct Rope {
    /// All nodes in the rope.
    pub nodes: Vec<RopeNode>,
    /// All segment constraints.
    pub segments: Vec<RopeSegment>,
    /// Global bending stiffness applied between non-adjacent nodes.
    pub bend_stiffness: f64,
}

// ---------------------------------------------------------------------------
// create_rope
// ---------------------------------------------------------------------------

/// Create a straight rope hanging vertically from y=0 downward.
///
/// * `n` – number of nodes (segments = n − 1).
/// * `length` – total rope length (m).
/// * `mass_per_node` – mass assigned to every node (kg).
/// * `stiffness` – stretch constraint stiffness (0–1).
///
/// The first node (index 0) is pinned at the origin.
pub fn create_rope(n: usize, length: f64, mass_per_node: f64, stiffness: f64) -> Rope {
    assert!(n >= 2, "A rope needs at least 2 nodes");
    let seg_length = length / (n - 1) as f64;
    let mut nodes = Vec::with_capacity(n);
    for i in 0..n {
        let y = -(i as f64) * seg_length;
        nodes.push(RopeNode::new([0.0, y, 0.0], mass_per_node));
    }
    nodes[0].pinned = true;

    let mut segments = Vec::with_capacity(n - 1);
    for i in 0..n - 1 {
        segments.push(RopeSegment {
            i,
            j: i + 1,
            rest_length: seg_length,
            stiffness,
            damping: 0.0,
        });
    }

    Rope {
        nodes,
        segments,
        bend_stiffness: 0.0,
    }
}

// ---------------------------------------------------------------------------
// rope_catenary_sag
// ---------------------------------------------------------------------------

/// Maximum sag of a catenary rope at mid-span (m).
///
/// ```text
/// sag ≈ w·L² / (8·T)      (parabolic approximation)
/// ```
///
/// where `w = linear_density · g` (N/m), `L` = span (m), `T` = tension (N).
pub fn rope_catenary_sag(length: f64, span: f64, linear_density: f64, tension: f64) -> f64 {
    let g = 9.81;
    let w = linear_density * g;
    if tension < 1e-30 {
        return f64::INFINITY;
    }
    // Parabolic catenary approximation valid for sag << span.
    let _length = length; // catenary parameter — used for validation reference
    w * span * span / (8.0 * tension)
}

// ---------------------------------------------------------------------------
// rope_natural_frequency
// ---------------------------------------------------------------------------

/// Fundamental natural frequency of a taut rope (Hz).
///
/// ```text
/// f₁ = (1 / 2L) · √(T / μ)
/// ```
///
/// where `L` = length (m), `T` = tension (N), `μ` = linear density (kg/m).
pub fn rope_natural_frequency(length: f64, tension: f64, linear_density: f64) -> f64 {
    if length < 1e-30 || linear_density < 1e-30 {
        return 0.0;
    }
    (tension / linear_density).sqrt() / (2.0 * length)
}

// ---------------------------------------------------------------------------
// apply_rope_wind
// ---------------------------------------------------------------------------

/// Apply aerodynamic drag from a uniform wind to each rope node.
///
/// The drag force on each node is:
/// ```text
/// F_drag = 0.5 · ρ_air · C_d · A · |v_rel|² · v̂_rel
/// ```
/// simplified to `F_drag = drag_coeff · |v_rel|² · v̂_rel` where `drag_coeff`
/// absorbs density, area, and 0.5.
///
/// The impulse `F_drag · dt` is added to each un-pinned node's velocity.
pub fn apply_rope_wind(rope: &mut Rope, wind: [f64; 3], drag_coeff: f64, dt: f64) {
    for node in rope.nodes.iter_mut() {
        if node.pinned {
            continue;
        }
        let v_rel = sub3(wind, node.velocity);
        let speed_sq = dot3(v_rel, v_rel);
        if speed_sq < 1e-30 {
            continue;
        }
        let speed = speed_sq.sqrt();
        let f_drag = scale3(v_rel, drag_coeff * speed);
        let impulse = scale3(f_drag, dt / node.mass.max(1e-30));
        node.velocity = add3(node.velocity, impulse);
    }
}

// ---------------------------------------------------------------------------
// resolve_rope_distance_constraint
// ---------------------------------------------------------------------------

/// Project a single rope distance constraint (one PBD iteration).
///
/// Adjusts positions of nodes `i` and `j` so their distance equals
/// `rest_length`, weighted by inverse mass.
pub fn resolve_rope_distance_constraint(rope: &mut Rope, seg_idx: usize) {
    let seg = rope.segments[seg_idx].clone();
    let pi = rope.nodes[seg.i].position;
    let pj = rope.nodes[seg.j].position;
    let delta = sub3(pj, pi);
    let dist = len3(delta);
    if dist < 1e-15 {
        return;
    }
    let correction = (dist - seg.rest_length) / dist;
    let wi = if rope.nodes[seg.i].pinned {
        0.0
    } else {
        1.0 / rope.nodes[seg.i].mass.max(1e-30)
    };
    let wj = if rope.nodes[seg.j].pinned {
        0.0
    } else {
        1.0 / rope.nodes[seg.j].mass.max(1e-30)
    };
    let sum_w = wi + wj;
    if sum_w < 1e-30 {
        return;
    }
    let corr = scale3(delta, seg.stiffness * correction / sum_w);
    if !rope.nodes[seg.i].pinned {
        rope.nodes[seg.i].position = add3(pi, scale3(corr, wi));
    }
    if !rope.nodes[seg.j].pinned {
        rope.nodes[seg.j].position = sub3(pj, scale3(corr, wj));
    }
}

// ---------------------------------------------------------------------------
// rope_step_pbd
// ---------------------------------------------------------------------------

/// Advance the rope by one PBD time step.
///
/// 1. Apply gravity to velocities of un-pinned nodes.
/// 2. Predict new positions.
/// 3. Resolve all distance constraints for `n_iters` iterations.
/// 4. Update velocities from position change.
pub fn rope_step_pbd(rope: &mut Rope, gravity: [f64; 3], dt: f64, n_iters: usize) {
    let old_positions: Vec<[f64; 3]> = rope.nodes.iter().map(|n| n.position).collect();

    for node in rope.nodes.iter_mut() {
        if node.pinned {
            continue;
        }
        node.velocity = add3(node.velocity, scale3(gravity, dt));
        node.position = add3(node.position, scale3(node.velocity, dt));
    }

    for _ in 0..n_iters {
        for seg_idx in 0..rope.segments.len() {
            resolve_rope_distance_constraint(rope, seg_idx);
        }
    }

    let inv_dt = if dt > 1e-30 { 1.0 / dt } else { 0.0 };
    for (idx, node) in rope.nodes.iter_mut().enumerate() {
        if !node.pinned {
            node.velocity = scale3(sub3(node.position, old_positions[idx]), inv_dt);
        }
    }
}

// ---------------------------------------------------------------------------
// rope_end_reaction_force
// ---------------------------------------------------------------------------

/// Compute the reaction force at the pinned end of the rope (N).
///
/// Estimates the tension at the top pinned node by summing gravitational
/// loads from all nodes below it.  Returns `[0, 0, 0]` if no pinned node
/// exists.
pub fn rope_end_reaction_force(rope: &Rope) -> [f64; 3] {
    let g = 9.81f64;
    // Find the first pinned node.
    let _pin_idx = match rope.nodes.iter().position(|n| n.pinned) {
        Some(i) => i,
        None => return [0.0; 3],
    };
    // Sum gravitational weight of all nodes (the pinned node must support them).
    let total_mass: f64 = rope.nodes.iter().map(|n| n.mass).sum();
    [0.0, total_mass * g, 0.0]
}

// ---------------------------------------------------------------------------
// rope_total_energy
// ---------------------------------------------------------------------------

/// Total mechanical energy (kinetic + gravitational potential) of the rope (J).
///
/// Gravitational PE is measured relative to y = 0 (positive upward).
pub fn rope_total_energy(rope: &Rope, gravity: [f64; 3]) -> f64 {
    let g_magnitude = len3(gravity);
    let mut energy = 0.0;
    for node in &rope.nodes {
        let ke = 0.5 * node.mass * dot3(node.velocity, node.velocity);
        // PE: m * g * h, where h is the y-component of position.
        let pe = node.mass * g_magnitude * node.position[1];
        energy += ke + pe;
    }
    energy
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    // 1. create_rope: correct node count.
    #[test]
    fn test_create_rope_node_count() {
        let rope = create_rope(10, 5.0, 0.1, 1.0);
        assert_eq!(rope.nodes.len(), 10);
    }

    // 2. create_rope: correct segment count.
    #[test]
    fn test_create_rope_segment_count() {
        let rope = create_rope(10, 5.0, 0.1, 1.0);
        assert_eq!(rope.segments.len(), 9);
    }

    // 3. create_rope: first node is pinned.
    #[test]
    fn test_create_rope_first_pinned() {
        let rope = create_rope(5, 2.0, 0.5, 1.0);
        assert!(rope.nodes[0].pinned);
    }

    // 4. create_rope: subsequent nodes are not pinned.
    #[test]
    fn test_create_rope_others_not_pinned() {
        let rope = create_rope(5, 2.0, 0.5, 1.0);
        assert!(rope.nodes[1..].iter().all(|n| !n.pinned));
    }

    // 5. create_rope: rest lengths are equal.
    #[test]
    fn test_create_rope_rest_lengths() {
        let n = 6usize;
        let len = 3.0f64;
        let rope = create_rope(n, len, 0.1, 1.0);
        let seg_len = len / (n - 1) as f64;
        for seg in &rope.segments {
            assert!((seg.rest_length - seg_len).abs() < EPS);
        }
    }

    // 6. create_rope: nodes are spaced correctly.
    #[test]
    fn test_create_rope_node_spacing() {
        let rope = create_rope(5, 4.0, 0.1, 1.0);
        let seg_len = 4.0 / 4.0;
        for i in 0..rope.nodes.len() - 1 {
            let d = len3(sub3(rope.nodes[i + 1].position, rope.nodes[i].position));
            assert!((d - seg_len).abs() < EPS);
        }
    }

    // 7. rope_catenary_sag: proportional to w*L²/8T.
    #[test]
    fn test_catenary_sag_formula() {
        let sag = rope_catenary_sag(5.0, 4.0, 0.5, 100.0);
        let g = 9.81;
        let expected = 0.5 * g * 4.0 * 4.0 / (8.0 * 100.0);
        assert!(
            (sag - expected).abs() < EPS * expected,
            "sag={sag} expected={expected}"
        );
    }

    // 8. rope_catenary_sag: zero tension → infinity.
    #[test]
    fn test_catenary_sag_zero_tension() {
        assert!(rope_catenary_sag(5.0, 4.0, 0.5, 0.0).is_infinite());
    }

    // 9. rope_catenary_sag: larger tension → smaller sag.
    #[test]
    fn test_catenary_sag_larger_tension() {
        let s1 = rope_catenary_sag(5.0, 4.0, 0.5, 100.0);
        let s2 = rope_catenary_sag(5.0, 4.0, 0.5, 400.0);
        assert!(s1 > s2, "Higher tension should reduce sag");
    }

    // 10. rope_natural_frequency: positive for positive inputs.
    #[test]
    fn test_natural_frequency_positive() {
        let f = rope_natural_frequency(1.0, 100.0, 0.5);
        assert!(f > 0.0);
    }

    // 11. rope_natural_frequency: zero length → 0.
    #[test]
    fn test_natural_frequency_zero_length() {
        assert_eq!(rope_natural_frequency(0.0, 100.0, 0.5), 0.0);
    }

    // 12. rope_natural_frequency: formula check.
    #[test]
    fn test_natural_frequency_formula() {
        let l = 2.0f64;
        let t = 100.0f64;
        let mu = 0.5f64;
        let f = rope_natural_frequency(l, t, mu);
        let expected = (t / mu).sqrt() / (2.0 * l);
        assert!((f - expected).abs() < EPS * expected);
    }

    // 13. apply_rope_wind: zero wind → no velocity change.
    #[test]
    fn test_wind_zero() {
        let mut rope = create_rope(5, 2.0, 0.1, 1.0);
        apply_rope_wind(&mut rope, [0.0, 0.0, 0.0], 1.0, 0.01);
        assert!(rope.nodes.iter().all(|n| n.velocity == [0.0; 3]));
    }

    // 14. apply_rope_wind: non-zero wind imparts velocity.
    #[test]
    fn test_wind_nonzero() {
        let mut rope = create_rope(5, 2.0, 0.1, 1.0);
        apply_rope_wind(&mut rope, [10.0, 0.0, 0.0], 0.1, 0.01);
        let any_moved = rope.nodes.iter().any(|n| n.velocity[0].abs() > 1e-15);
        assert!(any_moved);
    }

    // 15. apply_rope_wind: pinned nodes are unaffected.
    #[test]
    fn test_wind_pinned_unaffected() {
        let mut rope = create_rope(5, 2.0, 0.1, 1.0);
        apply_rope_wind(&mut rope, [100.0, 0.0, 0.0], 1.0, 0.1);
        assert_eq!(rope.nodes[0].velocity, [0.0; 3]);
    }

    // 16. resolve_rope_distance_constraint: constraint converges.
    #[test]
    fn test_resolve_constraint_convergence() {
        let mut rope = create_rope(2, 1.0, 1.0, 1.0);
        // Stretch second node to 3.0 m.
        rope.nodes[1].position = [0.0, -3.0, 0.0];
        for _ in 0..200 {
            resolve_rope_distance_constraint(&mut rope, 0);
        }
        let dist = len3(sub3(rope.nodes[1].position, rope.nodes[0].position));
        let rest = rope.segments[0].rest_length;
        assert!(
            (dist - rest).abs() < 0.05,
            "dist={dist} should ≈ rest={rest}"
        );
    }

    // 17. resolve_rope_distance_constraint: pinned–pinned → no movement.
    #[test]
    fn test_resolve_constraint_pinned() {
        let mut rope = create_rope(2, 1.0, 1.0, 1.0);
        rope.nodes[1].pinned = true;
        rope.nodes[1].position = [0.0, -5.0, 0.0];
        let before = rope.nodes[1].position;
        resolve_rope_distance_constraint(&mut rope, 0);
        assert_eq!(rope.nodes[1].position, before);
    }

    // 18. rope_step_pbd: gravity causes nodes to fall.
    #[test]
    fn test_pbd_step_gravity() {
        let mut rope = create_rope(5, 2.0, 0.1, 1.0);
        let dt = 1.0 / 60.0;
        for _ in 0..60 {
            rope_step_pbd(&mut rope, [0.0, -9.81, 0.0], dt, 5);
        }
        // Last (free) node should be below its initial position.
        let last_y = rope.nodes.last().unwrap().position[1];
        assert!(last_y < -2.0, "Last node should fall, y={last_y}");
    }

    // 19. rope_step_pbd: pinned node stays at origin.
    #[test]
    fn test_pbd_step_pinned_fixed() {
        let mut rope = create_rope(5, 2.0, 0.1, 1.0);
        let orig = rope.nodes[0].position;
        rope_step_pbd(&mut rope, [0.0, -9.81, 0.0], 0.01, 5);
        assert_eq!(rope.nodes[0].position, orig);
    }

    // 20. rope_end_reaction_force: returns non-zero for a rope with weight.
    #[test]
    fn test_reaction_force_nonzero() {
        let rope = create_rope(5, 2.0, 1.0, 1.0);
        let f = rope_end_reaction_force(&rope);
        assert!(
            f[1] > 0.0,
            "Reaction force should support weight, got {:?}",
            f
        );
    }

    // 21. rope_end_reaction_force: scales with total mass.
    #[test]
    fn test_reaction_force_scales_mass() {
        let r1 = create_rope(5, 2.0, 1.0, 1.0);
        let r2 = create_rope(5, 2.0, 2.0, 1.0);
        let f1 = rope_end_reaction_force(&r1);
        let f2 = rope_end_reaction_force(&r2);
        assert!((f2[1] - 2.0 * f1[1]).abs() < EPS);
    }

    // 22. rope_end_reaction_force: no pinned node → [0,0,0].
    #[test]
    fn test_reaction_no_pin() {
        let mut rope = create_rope(3, 1.0, 1.0, 1.0);
        rope.nodes[0].pinned = false;
        let f = rope_end_reaction_force(&rope);
        // All nodes unpinned → still returns weight (sum of masses).
        // Actually our impl always finds pinned index; with no pin returns zero.
        // Since we un-pinned all, expect [0,0,0].
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    // 23. rope_total_energy: KE is zero for stationary rope.
    #[test]
    fn test_total_energy_zero_velocity() {
        let rope = create_rope(5, 2.0, 1.0, 1.0);
        let g = 9.81f64;
        let e = rope_total_energy(&rope, [0.0, -g, 0.0]);
        // All velocities are zero; PE should be negative (nodes below y=0).
        let total_mass: f64 = rope.nodes.iter().map(|n| n.mass).sum();
        // PE can be negative when nodes are below origin.
        assert!(
            e < total_mass * g,
            "Energy should reflect downward positions"
        );
    }

    // 24. rope_total_energy: adding velocity increases total energy.
    #[test]
    fn test_total_energy_kinetic() {
        let rope_still = create_rope(2, 1.0, 1.0, 1.0);
        let mut rope_moving = create_rope(2, 1.0, 1.0, 1.0);
        // Give the last node a velocity.
        rope_moving.nodes[1].velocity = [3.0, 0.0, 0.0];
        let e_still = rope_total_energy(&rope_still, [0.0, -9.81, 0.0]);
        let e_moving = rope_total_energy(&rope_moving, [0.0, -9.81, 0.0]);
        let ke = 0.5 * 1.0 * 9.0; // 0.5 * m * v²
        // Moving rope should have exactly ke more energy than still rope.
        assert!(
            (e_moving - e_still - ke).abs() < EPS,
            "Energy with velocity should be higher by KE: e_moving={e_moving}, e_still={e_still}, ke={ke}"
        );
    }

    // 25. rope_natural_frequency: larger length → lower frequency.
    #[test]
    fn test_natural_frequency_length_effect() {
        let f1 = rope_natural_frequency(1.0, 100.0, 0.5);
        let f2 = rope_natural_frequency(2.0, 100.0, 0.5);
        assert!(f1 > f2, "Longer rope should have lower frequency");
    }

    // 26. rope_catenary_sag: proportional to L².
    #[test]
    fn test_catenary_sag_span_squared() {
        let s1 = rope_catenary_sag(10.0, 2.0, 1.0, 500.0);
        let s2 = rope_catenary_sag(10.0, 4.0, 1.0, 500.0);
        assert!((s2 - 4.0 * s1).abs() < EPS * s1, "Sag scales as L²");
    }

    // 27. apply_rope_wind: drag scales with drag_coeff.
    #[test]
    fn test_wind_scales_with_drag_coeff() {
        let mut r1 = create_rope(3, 1.0, 1.0, 1.0);
        let mut r2 = create_rope(3, 1.0, 1.0, 1.0);
        apply_rope_wind(&mut r1, [10.0, 0.0, 0.0], 0.1, 0.01);
        apply_rope_wind(&mut r2, [10.0, 0.0, 0.0], 0.2, 0.01);
        let v1 = r1.nodes[1].velocity[0];
        let v2 = r2.nodes[1].velocity[0];
        assert!(
            (v2 - 2.0 * v1).abs() < EPS,
            "Velocity should scale with drag_coeff"
        );
    }

    // 28. RopeNode: default velocity is zero.
    #[test]
    fn test_rope_node_default_velocity() {
        let n = RopeNode::new([1.0, 2.0, 3.0], 0.5);
        assert_eq!(n.velocity, [0.0; 3]);
    }

    // 29. rope_step_pbd: velocity is updated after step.
    #[test]
    fn test_pbd_velocity_updated() {
        let mut rope = create_rope(3, 1.0, 1.0, 1.0);
        rope_step_pbd(&mut rope, [0.0, -9.81, 0.0], 0.01, 1);
        // At least one non-pinned node should have non-zero velocity.
        assert!(rope.nodes[1..].iter().any(|n| n.velocity[1].abs() > 1e-10));
    }

    // 30. resolve_rope_distance_constraint: both nodes move toward each other.
    #[test]
    fn test_resolve_moves_both_nodes() {
        let mut rope = create_rope(2, 1.0, 1.0, 1.0);
        rope.nodes[0].pinned = false;
        rope.nodes[0].position = [0.0, 0.0, 0.0];
        rope.nodes[1].position = [0.0, -5.0, 0.0];
        let p0_before = rope.nodes[0].position;
        let p1_before = rope.nodes[1].position;
        resolve_rope_distance_constraint(&mut rope, 0);
        let p0_after = rope.nodes[0].position;
        let p1_after = rope.nodes[1].position;
        assert_ne!(p0_after, p0_before, "Node 0 should have moved");
        assert_ne!(p1_after, p1_before, "Node 1 should have moved");
    }
}
