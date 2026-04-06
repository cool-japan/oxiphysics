//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::Rope;

/// Dot product of two 3-vectors.
pub fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
/// Euclidean length of a 3-vector.
pub fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}
/// Normalise a 3-vector.  Returns the zero vector when the input is
/// near-zero to avoid NaN propagation.
pub fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-15 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}
/// Component-wise addition.
pub fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
/// Component-wise subtraction.
pub fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
/// Scalar multiplication.
pub fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}
/// Cross product of two 3-vectors.
#[allow(dead_code)]
pub fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::HairStrand;
    use crate::HairSystem;
    use crate::Rope;
    // taut_string_frequency from super::*

    use crate::material_point::len3;
    use crate::material_point::sub3;

    pub(super) const EPS: f64 = 1e-9;
    #[test]
    fn test_new_straight_node_count_and_spacing() {
        let n = 6;
        let rope = Rope::new_straight(n, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 5.0, 1.0);
        assert_eq!(rope.nodes.len(), n, "Should have {n} nodes");
        let expected_seg = 5.0 / (n - 1) as f64;
        for i in 0..n - 1 {
            let d = len3(sub3(rope.nodes[i + 1].position, rope.nodes[i].position));
            assert!(
                (d - expected_seg).abs() < 1e-10,
                "Segment {i} length {d} ≠ {expected_seg}"
            );
        }
    }
    #[test]
    fn test_fixed_end_does_not_move() {
        let mut rope = Rope::new_straight(5, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 4.0, 0.2);
        rope.fix_end(0);
        let fixed_pos = rope.nodes[0].position;
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..60 {
            rope.step(1.0 / 60.0, gravity);
        }
        let diff = len3(sub3(rope.nodes[0].position, fixed_pos));
        assert!(diff < EPS, "Fixed node moved by {diff}");
    }
    #[test]
    fn test_total_length() {
        let rope = Rope::new_straight(4, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0, 1.0);
        let tl = rope.total_length();
        assert!(
            (tl - 3.0).abs() < 1e-10,
            "total_length should be 3.0, got {tl}"
        );
    }
    #[test]
    fn test_hair_system_strand_count() {
        let mut sys = HairSystem::new();
        assert_eq!(sys.strand_count(), 0);
        sys.add_strand(HairStrand::new([0.0, 1.0, 0.0], 1.0, 5));
        sys.add_strand(HairStrand::new([0.1, 1.0, 0.0], 1.0, 5));
        assert_eq!(sys.strand_count(), 2);
    }
    #[test]
    fn test_rope_end_falls_under_gravity() {
        let mut rope = Rope::new_straight(5, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 1.0, 0.2);
        rope.fix_end(0);
        let initial_y = rope.nodes[rope.nodes.len() - 1].position[1];
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..120 {
            rope.step(1.0 / 60.0, gravity);
        }
        let final_y = rope.nodes[rope.nodes.len() - 1].position[1];
        assert!(
            final_y < initial_y,
            "Rope end should fall: initial_y={initial_y}, final_y={final_y}"
        );
    }
    #[test]
    fn test_collide_sphere_pushes_nodes() {
        let mut rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        let center = [2.0, 0.0, 0.0];
        let radius = 1.5;
        let collisions = rope.collide_sphere(center, radius);
        for c in &collisions {
            let pos = rope.nodes[c.node_index].position;
            let dist = len3(sub3(pos, center));
            assert!(
                dist >= radius - 1e-6,
                "Node {} should be at or outside sphere, dist={dist}",
                c.node_index
            );
        }
    }
    #[test]
    fn test_collide_plane_pushes_nodes_up() {
        let mut rope = Rope::new_straight(5, [0.0, -1.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        let collisions = rope.collide_plane([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(!collisions.is_empty(), "Nodes below plane should collide");
        for node in &rope.nodes {
            assert!(
                node.position[1] >= -1e-10,
                "All nodes should be at or above y=0 after collision"
            );
        }
    }
    #[test]
    fn test_self_collision_separates_nodes() {
        let mut rope = Rope::new_straight(10, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 9.0, 1.0);
        rope.nodes[8].position = [0.05, 0.0, 0.0];
        let min_dist = 0.5;
        rope.resolve_self_collision(min_dist, 2);
        let dist = len3(sub3(rope.nodes[8].position, rope.nodes[0].position));
        assert!(
            dist >= min_dist - 1e-6,
            "Nodes should be at least {min_dist} apart, got {dist}"
        );
    }
    #[test]
    fn test_collide_cylinder_pushes_nodes() {
        let mut rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        let collisions = rope.collide_cylinder([0.0, 0.0, -5.0], [0.0, 0.0, 5.0], 2.0);
        for c in &collisions {
            assert!(c.penetration > 0.0);
        }
    }
    #[test]
    fn test_catenary_positions_count() {
        let positions = Rope::catenary_positions([0.0, 0.0, 0.0], [10.0, 0.0, 0.0], 15.0, 20);
        assert_eq!(positions.len(), 20);
    }
    #[test]
    fn test_catenary_sags() {
        let start = [0.0, 5.0, 0.0];
        let end = [10.0, 5.0, 0.0];
        let positions = Rope::catenary_positions(start, end, 15.0, 21);
        let mid_y = positions[10][1];
        assert!(
            mid_y < 5.0,
            "Catenary midpoint should sag below 5.0, got {mid_y}"
        );
    }
    #[test]
    fn test_catenary_straight_when_short() {
        let start = [0.0, 0.0, 0.0];
        let end = [10.0, 0.0, 0.0];
        let positions = Rope::catenary_positions(start, end, 9.0, 5);
        for (i, pos) in positions.iter().enumerate() {
            let t = i as f64 / 4.0;
            assert!((pos[0] - 10.0 * t).abs() < 1e-6, "Expected straight line x");
        }
    }
    #[test]
    fn test_is_taut() {
        let rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        assert!(rope.is_taut(0.01), "Fresh rope should be taut");
    }
    #[test]
    fn test_enforce_max_stretch() {
        let mut rope = Rope::new_straight(3, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 1.0);
        rope.nodes[2].position = [10.0, 0.0, 0.0];
        rope.enforce_max_stretch(1.5);
        let max_len = rope.segment_length * 1.5;
        for i in 0..rope.nodes.len() - 1 {
            let dist = len3(sub3(rope.nodes[i + 1].position, rope.nodes[i].position));
            assert!(
                dist <= max_len + 1e-6,
                "Segment {i} length {dist} exceeds max {max_len}"
            );
        }
    }
    #[test]
    fn test_segment_tensions_count() {
        let rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        let tensions = rope.segment_tensions();
        assert_eq!(tensions.len(), 4);
    }
    #[test]
    fn test_max_tension_at_rest() {
        let rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        assert!(
            rope.max_tension() < 1e-10,
            "Unstretched rope should have zero tension"
        );
    }
    #[test]
    fn test_kinetic_energy_zero_at_rest() {
        let rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        assert!(
            rope.kinetic_energy() < 1e-15,
            "Stationary rope should have zero KE"
        );
    }
    #[test]
    fn test_potential_energy_positive_above_ref() {
        let rope = Rope::new_straight(3, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 2.0, 1.0);
        let pe = rope.potential_energy(9.81, 0.0);
        assert!(
            pe > 0.0,
            "Rope above ref_y=0 should have positive PE, got {pe}"
        );
    }
    #[test]
    fn test_cross3() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(a, b);
        assert!((c[0]).abs() < EPS);
        assert!((c[1]).abs() < EPS);
        assert!((c[2] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_collide_sphere_no_collision() {
        let mut rope = Rope::new_straight(3, [5.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 1.0);
        let collisions = rope.collide_sphere([0.0, 0.0, 0.0], 1.0);
        assert!(
            collisions.is_empty(),
            "No nodes should collide with distant sphere"
        );
    }
    #[test]
    fn test_collide_plane_no_collision_above() {
        let mut rope = Rope::new_straight(3, [0.0, 5.0, 0.0], [1.0, 0.0, 0.0], 2.0, 1.0);
        let collisions = rope.collide_plane([0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(
            collisions.is_empty(),
            "Nodes above plane should not collide"
        );
    }
    #[test]
    fn test_fixed_nodes_not_moved_by_sphere() {
        let mut rope = Rope::new_straight(3, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 1.0);
        rope.fix_end(0);
        let original_pos = rope.nodes[0].position;
        let _collisions = rope.collide_sphere([0.0, 0.0, 0.0], 5.0);
        assert!(
            len3(sub3(rope.nodes[0].position, original_pos)) < EPS,
            "Fixed node should not be moved"
        );
    }
    #[test]
    fn test_hair_strand_root_fixed() {
        let strand = HairStrand::new([0.0, 2.0, 0.0], 1.0, 5);
        assert!(strand.rope.nodes[0].fixed, "Root node should be fixed");
    }
    #[test]
    fn test_curvature_at_boundary() {
        let rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        assert!(rope.curvature_at(0).abs() < EPS);
        assert!(rope.curvature_at(4).abs() < EPS);
    }
    #[test]
    fn test_curvature_at_interior_straight() {
        let rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        for i in 1..4 {
            assert!(
                rope.curvature_at(i).abs() < 1e-6,
                "Straight rope should have zero curvature at node {i}"
            );
        }
    }
}
#[cfg(test)]
mod rope_verlet_tests {

    // taut_string_frequency from super::*

    use crate::material_point::len3;
    use crate::material_point::sub3;

    use crate::rope::RopeSegment;
    use crate::rope::RopeVerlet;

    #[test]
    fn test_rope_verlet_new_particle_count() {
        let n = 8;
        let r = RopeVerlet::new(n, [0.0; 3], [1.0, 0.0, 0.0], 0.5, 1.0);
        assert_eq!(r.particles.len(), n);
        assert_eq!(r.velocities.len(), n);
        assert_eq!(r.masses.len(), n);
        assert_eq!(r.rest_lengths.len(), n - 1);
    }
    #[test]
    fn test_rope_verlet_step_changes_positions() {
        let mut r = RopeVerlet::new(5, [0.0, 5.0, 0.0], [0.0, -1.0, 0.0], 0.5, 0.2);
        r.fixed[0] = true;
        let pos_before: Vec<_> = r.particles.clone();
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..30 {
            r.step(1.0 / 60.0, gravity, 1.0, 0.01);
        }
        let moved = (1..r.particles.len()).any(|i| {
            let d = len3(sub3(r.particles[i], pos_before[i]));
            d > 1e-6
        });
        assert!(moved, "Particles should move under gravity");
    }
    #[test]
    fn test_rope_verlet_elastic_force_restores() {
        let seg = RopeSegment {
            p0: [0.0; 3],
            p1: [2.0, 0.0, 0.0],
            rest_length: 1.0,
            mass: 1.0,
            k_stretch: 100.0,
        };
        let (f0, f1) = seg.elastic_force();
        assert!(
            f0[0] > 0.0,
            "Force on p0 should be in +x direction, got {:?}",
            f0
        );
        assert!(
            f1[0] < 0.0,
            "Force on p1 should be in -x direction, got {:?}",
            f1
        );
    }
    #[test]
    fn test_rope_verlet_total_length() {
        let n = 5;
        let seg_len = 1.0;
        let r = RopeVerlet::new(n, [0.0; 3], [1.0, 0.0, 0.0], seg_len, 1.0);
        let total = r.total_length();
        let expected = seg_len * (n - 1) as f64;
        assert!(
            (total - expected).abs() < 1e-10,
            "total_length should be {expected}, got {total}"
        );
    }
}
/// Compute the `n`-th natural frequency (Hz) of a taut string.
///
/// For a string of length `L` (m), tension `T` (N) and linear mass density
/// `mu` (kg/m), the `n`-th mode frequency is:
///
/// ```text
/// f_n = n / (2 L) * sqrt(T / μ)
/// ```
///
/// `n` must be >= 1.
#[allow(dead_code)]
pub fn taut_string_frequency(n: usize, length: f64, tension: f64, linear_density: f64) -> f64 {
    assert!(n >= 1, "mode number must be >= 1");
    let wave_speed = (tension / linear_density.max(1e-30)).sqrt();
    (n as f64) * wave_speed / (2.0 * length.max(1e-30))
}
/// Compute the first `n_modes` standing-wave frequencies for a taut string.
///
/// Returns a `Vec` of length `n_modes` with frequencies in ascending order.
#[allow(dead_code)]
pub fn taut_string_modes(
    n_modes: usize,
    length: f64,
    tension: f64,
    linear_density: f64,
) -> Vec<f64> {
    (1..=n_modes)
        .map(|n| taut_string_frequency(n, length, tension, linear_density))
        .collect()
}
/// Transverse displacement of the `n`-th standing-wave mode at position `x`
/// and time `t`:
///
/// `y(x, t) = amplitude * sin(n π x / L) * cos(2π f_n t)`
#[allow(dead_code)]
pub fn standing_wave_displacement(
    x: f64,
    t: f64,
    n: usize,
    length: f64,
    tension: f64,
    linear_density: f64,
    amplitude: f64,
) -> f64 {
    let freq = taut_string_frequency(n, length, tension, linear_density);
    let k = (n as f64) * std::f64::consts::PI / length.max(1e-30);
    amplitude * (k * x).sin() * (2.0 * std::f64::consts::PI * freq * t).cos()
}
/// Detect whether the rope is coiling.
///
/// Coiling occurs when the cumulative curvature over the lower portion of the
/// rope exceeds a threshold — indicating that the rope cannot bear its own
/// weight under compression.
///
/// Returns the index of the first node where coiling is detected, or `None`.
#[allow(dead_code)]
pub fn detect_coiling(rope: &Rope, curvature_threshold: f64, window: usize) -> Option<usize> {
    let n = rope.nodes.len();
    if n < 3 || window == 0 {
        return None;
    }
    let mut acc = 0.0;
    let start = n.saturating_sub(window + 1);
    for i in (start.max(1)..n.saturating_sub(1)).rev() {
        acc += rope.curvature_at(i);
        if acc > curvature_threshold {
            return Some(i);
        }
    }
    None
}
/// Maximum tension threshold for coiling: if any segment tension exceeds
/// `max_tension`, the rope is considered taut enough that coiling cannot occur.
#[allow(dead_code)]
pub fn coil_tension_exceeded(rope: &Rope, max_tension: f64) -> bool {
    rope.max_tension() > max_tension
}
#[cfg(test)]
mod cosserat_tests {
    use super::*;

    use crate::Rope;
    // taut_string_frequency from super::*
    use crate::material_point::dot3;
    use crate::material_point::len3;
    use crate::material_point::sub3;
    use crate::rope::KirchhoffRod;
    use crate::rope::MaterialFrame;

    use crate::rope::ViscoelasticRope;
    use crate::rope::ViscoelasticSegment;
    use crate::rope::coil_tension_exceeded;
    use crate::rope::detect_coiling;
    use crate::rope::standing_wave_displacement;
    use crate::rope::taut_string_modes;
    #[test]
    fn test_material_frame_from_tangent_orthonormal() {
        let tangent = [0.0, 0.0, 1.0];
        let frame = MaterialFrame::from_tangent(tangent);
        assert!((frame.d3[2] - 1.0).abs() < 1e-10);
        assert!((len3(frame.d1) - 1.0).abs() < 1e-10);
        assert!((len3(frame.d2) - 1.0).abs() < 1e-10);
        assert!((len3(frame.d3) - 1.0).abs() < 1e-10);
        assert!(dot3(frame.d1, frame.d2).abs() < 1e-10);
        assert!(dot3(frame.d1, frame.d3).abs() < 1e-10);
        assert!(dot3(frame.d2, frame.d3).abs() < 1e-10);
    }
    #[test]
    fn test_material_frame_twist_preserves_d3() {
        let tangent = [1.0, 0.0, 0.0];
        let frame = MaterialFrame::from_tangent(tangent);
        let twisted = frame.twist(std::f64::consts::FRAC_PI_4);
        for k in 0..3 {
            assert!(
                (twisted.d3[k] - frame.d3[k]).abs() < 1e-10,
                "d3 component {k} changed after twist"
            );
        }
        assert!((len3(twisted.d1) - 1.0).abs() < 1e-10);
        assert!((len3(twisted.d2) - 1.0).abs() < 1e-10);
        assert!(dot3(twisted.d1, twisted.d2).abs() < 1e-10);
    }
    #[test]
    fn test_kirchhoff_rod_straight_zero_bending_energy() {
        let rod = KirchhoffRod::new_straight(5, [0.0; 3], [1.0, 0.0, 0.0], 0.1, 1.0, 0.5);
        let e_bend = rod.bending_energy();
        assert!(
            e_bend.abs() < 1e-12,
            "Straight rod should have zero bending energy, got {e_bend}"
        );
    }
    #[test]
    fn test_kirchhoff_rod_straight_zero_twisting_energy() {
        let rod = KirchhoffRod::new_straight(5, [0.0; 3], [1.0, 0.0, 0.0], 0.1, 1.0, 0.5);
        let e_twist = rod.twisting_energy();
        assert!(
            e_twist.abs() < 1e-12,
            "Straight rod should have zero twisting energy, got {e_twist}"
        );
    }
    #[test]
    fn test_kirchhoff_rod_update_frames_preserves_orthogonality() {
        let mut rod = KirchhoffRod::new_straight(4, [0.0; 3], [1.0, 0.0, 0.0], 0.25, 1.0, 0.5);
        rod.positions[2][1] = 0.1;
        rod.positions[3][1] = 0.2;
        rod.update_frames();
        for (s, frame) in rod.frames.iter().enumerate() {
            assert!(
                (len3(frame.d1) - 1.0).abs() < 1e-10,
                "frame {s} d1 not unit"
            );
            assert!(
                (len3(frame.d2) - 1.0).abs() < 1e-10,
                "frame {s} d2 not unit"
            );
            assert!(
                (len3(frame.d3) - 1.0).abs() < 1e-10,
                "frame {s} d3 not unit"
            );
            assert!(
                dot3(frame.d1, frame.d2).abs() < 1e-10,
                "frame {s} d1·d2 not zero"
            );
            assert!(
                dot3(frame.d1, frame.d3).abs() < 1e-10,
                "frame {s} d1·d3 not zero"
            );
        }
    }
    #[test]
    fn test_taut_string_frequency_fundamental() {
        let f1 = taut_string_frequency(1, 1.0, 100.0, 0.01);
        assert!((f1 - 50.0).abs() < 1e-6, "f1 should be 50 Hz, got {f1}");
    }
    #[test]
    fn test_taut_string_modes_ordered() {
        let modes = taut_string_modes(5, 1.0, 100.0, 0.01);
        assert_eq!(modes.len(), 5);
        for i in 1..5 {
            assert!(
                modes[i] > modes[i - 1],
                "modes should be in ascending order"
            );
        }
    }
    #[test]
    fn test_taut_string_harmonics_integer_ratios() {
        let f1 = taut_string_frequency(1, 1.0, 100.0, 0.01);
        let f2 = taut_string_frequency(2, 1.0, 100.0, 0.01);
        let f3 = taut_string_frequency(3, 1.0, 100.0, 0.01);
        assert!(
            (f2 / f1 - 2.0).abs() < 1e-10,
            "2nd harmonic should be 2× fundamental"
        );
        assert!(
            (f3 / f1 - 3.0).abs() < 1e-10,
            "3rd harmonic should be 3× fundamental"
        );
    }
    #[test]
    fn test_standing_wave_nodes_zero() {
        let _y0 = standing_wave_displacement(0.0, 0.0, 1, 1.0, 100.0, 0.01, 1.0);
        let _y_l = standing_wave_displacement(1.0, 0.0, 1, 1.0, 100.0, 0.01, 1.0);
        assert!(_y0.abs() < 1e-12, "displacement at x=0 should be zero");
        assert!(_y_l.abs() < 1e-12, "displacement at x=L should be zero");
    }
    #[test]
    fn test_standing_wave_antinode_maximum() {
        let amp = 0.05;
        let y = standing_wave_displacement(0.5, 0.0, 1, 1.0, 100.0, 0.01, amp);
        assert!(
            (y - amp).abs() < 1e-10,
            "antinode displacement should equal amplitude, got {y}"
        );
    }
    #[test]
    fn test_detect_coiling_straight_rope_no_coiling() {
        let rope = Rope::new_straight(8, [0.0; 3], [1.0, 0.0, 0.0], 7.0, 1.0);
        let result = detect_coiling(&rope, 0.5, 4);
        assert!(result.is_none(), "straight rope should not coil");
    }
    #[test]
    fn test_detect_coiling_high_curvature() {
        let mut rope = Rope::new_straight(6, [0.0; 3], [1.0, 0.0, 0.0], 5.0, 1.0);
        rope.nodes[2].position = [1.0, 2.0, 0.0];
        rope.nodes[3].position = [2.0, 0.0, 0.0];
        rope.nodes[4].position = [3.0, 2.0, 0.0];
        let result = detect_coiling(&rope, 0.01, 5);
        assert!(result.is_some(), "bent rope should detect coiling");
    }
    #[test]
    fn test_coil_tension_exceeded() {
        let mut rope = Rope::new_straight(3, [0.0; 3], [1.0, 0.0, 0.0], 2.0, 1.0);
        rope.nodes[2].position = [10.0, 0.0, 0.0];
        assert!(
            coil_tension_exceeded(&rope, 0.001),
            "stretched rope should exceed tension threshold"
        );
        assert!(
            !coil_tension_exceeded(&rope, 1e9),
            "very high threshold should not be exceeded"
        );
    }
    #[test]
    fn test_viscoelastic_segment_zero_force_at_rest() {
        let mut seg = ViscoelasticSegment::new(0, 1, 1.0, 100.0, 10.0);
        let positions = [[0.0; 3], [1.0, 0.0, 0.0]];
        let (fa, fb) = seg.kelvin_voigt_force(&positions, 0.01);
        assert!(len3(fa) < 1e-10, "no force at rest, got {:?}", fa);
        assert!(len3(fb) < 1e-10, "no force at rest, got {:?}", fb);
    }
    #[test]
    fn test_viscoelastic_segment_stretched_force_direction() {
        let mut seg = ViscoelasticSegment::new(0, 1, 1.0, 100.0, 0.0);
        let positions = [[0.0; 3], [1.5, 0.0, 0.0]];
        let (fa, _fb) = seg.kelvin_voigt_force(&positions, 1.0);
        assert!(
            fa[0] > 0.0,
            "stretched segment: force on a should be in +x, got {:?}",
            fa
        );
    }
    #[test]
    fn test_viscoelastic_rope_creation() {
        let rope =
            ViscoelasticRope::new_straight(5, [0.0; 3], [1.0, 0.0, 0.0], 0.25, 1.0, 1000.0, 10.0);
        assert_eq!(rope.positions.len(), 5);
        assert_eq!(rope.segments.len(), 4);
    }
    #[test]
    fn test_viscoelastic_rope_step_moves_free_nodes() {
        let mut rope = ViscoelasticRope::new_straight(
            4,
            [0.0, 5.0, 0.0],
            [0.0, -1.0, 0.0],
            0.5,
            1.0,
            500.0,
            5.0,
        );
        rope.fixed[0] = true;
        let init_y = rope.positions[3][1];
        let gravity = [0.0, -9.81, 0.0];
        for _ in 0..30 {
            rope.step(1.0 / 60.0, gravity);
        }
        let final_y = rope.positions[3][1];
        assert!(final_y < init_y, "free node should fall under gravity");
    }
    #[test]
    fn test_viscoelastic_rope_elastic_energy_zero_at_rest() {
        let rope =
            ViscoelasticRope::new_straight(4, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 1.0, 100.0, 5.0);
        let e = rope.total_elastic_energy();
        assert!(
            e.abs() < 1e-10,
            "unstretched rope should have zero elastic energy, got {e}"
        );
    }
    #[test]
    fn test_viscoelastic_rope_kinetic_energy_zero_at_rest() {
        let rope =
            ViscoelasticRope::new_straight(4, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 1.0, 100.0, 5.0);
        let ke = rope.total_kinetic_energy();
        assert!(
            ke.abs() < 1e-15,
            "stationary rope should have zero KE, got {ke}"
        );
    }
    #[test]
    fn test_kirchhoff_rod_darboux_vector_boundary_returns_zero() {
        let rod = KirchhoffRod::new_straight(4, [0.0; 3], [1.0, 0.0, 0.0], 0.1, 1.0, 0.5);
        let dv_start = rod.darboux_vector(0);
        let dv_end = rod.darboux_vector(3);
        assert_eq!(dv_start, [0.0; 3]);
        assert_eq!(dv_end, [0.0; 3]);
    }
    #[test]
    fn test_catenary_tension_node_count() {
        let rope = Rope::new_straight(8, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 5.0, 1.0);
        let tensions = rope.compute_catenary_tension(9.81);
        assert_eq!(tensions.len(), 8, "Should have one tension per node");
    }
    #[test]
    fn test_catenary_tension_nonnegative() {
        let rope = Rope::new_straight(6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0, 1.0);
        let tensions = rope.compute_catenary_tension(9.81);
        for (i, &t) in tensions.iter().enumerate() {
            assert!(
                t >= 0.0,
                "Tension at node {i} should be non-negative, got {t}"
            );
        }
    }
    #[test]
    fn test_catenary_tension_higher_at_endpoints() {
        let n = 9;
        let mut rope = Rope::new_straight(n, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 1.0);
        rope.nodes[n / 2].position[1] = -0.3;
        rope.nodes[0].fixed = true;
        rope.nodes[n - 1].fixed = true;
        let tensions = rope.compute_catenary_tension(9.81);
        let t_mid = tensions[n / 2];
        let t_end = tensions[0];
        assert!(
            t_end >= t_mid,
            "Catenary tension should be at least as high at endpoints: t_end={t_end}, t_mid={t_mid}"
        );
    }
    #[test]
    fn test_dynamic_stiffness_positive() {
        let rope = Rope::new_straight(6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0, 1.0);
        let ks = rope.compute_dynamic_stiffness(70e9, 0.005, 9.81);
        assert!(ks > 0.0, "Dynamic stiffness should be positive, got {ks}");
    }
    #[test]
    fn test_dynamic_stiffness_scales_with_modulus() {
        let rope = Rope::new_straight(6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0, 1.0);
        let ks_soft = rope.compute_dynamic_stiffness(1e6, 0.005, 9.81);
        let ks_stiff = rope.compute_dynamic_stiffness(200e9, 0.005, 9.81);
        assert!(
            ks_stiff > ks_soft,
            "Stiffer material → larger dynamic stiffness: soft={ks_soft}, stiff={ks_stiff}"
        );
    }
    #[test]
    fn test_dynamic_stiffness_scales_with_radius() {
        let rope = Rope::new_straight(6, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 3.0, 1.0);
        let ks_thin = rope.compute_dynamic_stiffness(70e9, 0.001, 9.81);
        let ks_thick = rope.compute_dynamic_stiffness(70e9, 0.01, 9.81);
        assert!(
            ks_thick > ks_thin,
            "Thicker rope → larger bending stiffness: thin={ks_thin}, thick={ks_thick}"
        );
    }
    #[test]
    fn test_twist_constraint_straight_rope_unchanged() {
        let mut rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        let positions_before: Vec<_> = rope.nodes.iter().map(|n| n.position).collect();
        rope.apply_twist_constraint(0.1, 0.0, 1.0 / 60.0);
        for (i, (before, after)) in positions_before
            .iter()
            .zip(rope.nodes.iter().map(|n| n.position))
            .enumerate()
        {
            let dist = len3(sub3(*before, after));
            assert!(
                dist < 1e-10,
                "Straight rope should not be modified by twist constraint at node {i}: dist={dist}"
            );
        }
    }
    #[test]
    fn test_twist_constraint_does_not_move_fixed_nodes() {
        let mut rope = Rope::new_straight(5, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        rope.nodes[0].fixed = true;
        rope.nodes[4].fixed = true;
        rope.nodes[2].position[2] = 0.5;
        let pos_first = rope.nodes[0].position;
        let pos_last = rope.nodes[4].position;
        rope.apply_twist_constraint(0.05, 0.0, 1.0 / 60.0);
        let dist_first = len3(sub3(rope.nodes[0].position, pos_first));
        let dist_last = len3(sub3(rope.nodes[4].position, pos_last));
        assert!(
            dist_first < 1e-12,
            "Fixed node 0 should not move: dist={dist_first}"
        );
        assert!(
            dist_last < 1e-12,
            "Fixed node 4 should not move: dist={dist_last}"
        );
    }
    #[test]
    fn test_catenary_tension_scales_with_gravity() {
        let mut rope = Rope::new_straight(7, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 4.0, 1.0);
        rope.nodes[3].position[1] = -0.2;
        let tensions_g1 = rope.compute_catenary_tension(1.0);
        let tensions_g2 = rope.compute_catenary_tension(2.0);
        let avg_t1: f64 = tensions_g1.iter().sum::<f64>() / tensions_g1.len() as f64;
        let avg_t2: f64 = tensions_g2.iter().sum::<f64>() / tensions_g2.len() as f64;
        assert!(
            avg_t2 > avg_t1,
            "Higher gravity should produce higher tension: g=1→{avg_t1}, g=2→{avg_t2}"
        );
    }
    #[test]
    fn test_twist_constraint_short_rope_no_panic() {
        let mut rope = Rope::new_straight(2, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 1.0);
        rope.apply_twist_constraint(0.1, 0.0, 1.0 / 60.0);
        assert_eq!(rope.nodes.len(), 2);
    }
}
