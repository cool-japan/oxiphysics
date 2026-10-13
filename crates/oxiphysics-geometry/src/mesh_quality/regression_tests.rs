// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Exact-value regression tests for the tetrahedron quality metrics.
//!
//! The dihedral angles and the scaled Jacobian of well-known elements (the regular
//! tetrahedron and the unit right-corner tetrahedron) are known in closed form, so
//! these tests pin the metrics to those values for both vertex windings.

use super::*;

/// Dihedral angle of the regular tetrahedron in degrees, `acos(1/3)`.
fn regular_dihedral_deg() -> f64 {
    (1.0_f64 / 3.0).acos().to_degrees()
}

/// A regular tetrahedron with vertices on alternate corners of the cube
/// `[-1, 1]^3` in the winding with a *negative* edge-matrix determinant.
fn regular_tet_negative() -> [[f64; 3]; 4] {
    [
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ]
}

/// The same regular tetrahedron with `v1` and `v2` swapped (positive winding).
fn regular_tet_positive() -> [[f64; 3]; 4] {
    let [v0, v1, v2, v3] = regular_tet_negative();
    [v0, v2, v1, v3]
}

/// The unit right-corner tetrahedron.
fn right_corner_tet() -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

/// Deterministic pseudo-random numbers in `[-1, 1)` (64-bit LCG).
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let mantissa = (self.0 >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * mantissa - 1.0
    }

    fn point(&mut self) -> [f64; 3] {
        [self.next_unit(), self.next_unit(), self.next_unit()]
    }
}

#[test]
fn regular_tet_has_six_equal_dihedral_angles_for_both_windings() {
    let expected = regular_dihedral_deg();
    assert!(
        (expected - 70.528_779_365_5).abs() < 1e-9,
        "acos(1/3) = {expected}"
    );
    for tet in [regular_tet_negative(), regular_tet_positive()] {
        let q = compute_tet_quality(tet[0], tet[1], tet[2], tet[3]);
        assert!(
            (q.min_angle_deg - expected).abs() < 1e-3,
            "min dihedral {} != {expected}",
            q.min_angle_deg
        );
        assert!(
            (q.max_angle_deg - expected).abs() < 1e-3,
            "max dihedral {} != {expected}",
            q.max_angle_deg
        );
        assert!(
            q.skewness < 1e-4,
            "regular tet skewness {} must be ~0",
            q.skewness
        );
        assert!(q.is_acceptable(), "a regular tet must be acceptable: {q:?}");
    }
}

#[test]
fn right_corner_tet_dihedral_angles_are_exact() {
    // Three right dihedral angles along the edges at the right-angle corner and
    // three angles acos(1/sqrt(3)) along the edges of the slanted face.
    let tet = right_corner_tet();
    let slanted = (1.0_f64 / 3.0_f64.sqrt()).acos().to_degrees();
    let q = compute_tet_quality(tet[0], tet[1], tet[2], tet[3]);
    assert!(
        (q.min_angle_deg - slanted).abs() < 1e-9,
        "min dihedral {} != {slanted}",
        q.min_angle_deg
    );
    assert!(
        (q.max_angle_deg - 90.0).abs() < 1e-9,
        "max dihedral {} != 90",
        q.max_angle_deg
    );
    let (min_d, max_d) = tet_dihedral_angle_extremes(tet[0], tet[1], tet[2], tet[3]);
    assert!((min_d - slanted).abs() < 1e-9 && (max_d - 90.0).abs() < 1e-9);
}

#[test]
fn dihedral_angles_do_not_depend_on_winding_or_vertex_order_of_the_base() {
    let tet = right_corner_tet();
    let a = compute_tet_quality(tet[0], tet[1], tet[2], tet[3]);
    // Swapping two vertices flips the winding.
    let b = compute_tet_quality(tet[0], tet[2], tet[1], tet[3]);
    assert!((a.min_angle_deg - b.min_angle_deg).abs() < 1e-9);
    assert!((a.max_angle_deg - b.max_angle_deg).abs() < 1e-9);
    assert!((a.jacobian - b.jacobian).abs() < 1e-12);
}

#[test]
fn regular_tet_scaled_jacobian_is_one() {
    let pos = regular_tet_positive();
    let neg = regular_tet_negative();
    let single = tet_scaled_jacobian(pos[0], pos[1], pos[2], pos[3]);
    assert!((single - 1.0).abs() < 1e-12, "single-corner value {single}");
    let min_corner = tet_min_scaled_jacobian(pos[0], pos[1], pos[2], pos[3]);
    assert!(
        (min_corner - 1.0).abs() < 1e-12,
        "min over corners {min_corner}"
    );
    let q = compute_tet_quality(pos[0], pos[1], pos[2], pos[3]);
    assert!((q.jacobian - 1.0).abs() < 1e-12, "jacobian {}", q.jacobian);

    // The inverted winding is the mirror image: -1 for the signed metrics and the
    // magnitude 1 for the winding-independent element quality.
    let single_neg = tet_scaled_jacobian(neg[0], neg[1], neg[2], neg[3]);
    assert!(
        (single_neg + 1.0).abs() < 1e-12,
        "inverted value {single_neg}"
    );
    let min_neg = tet_min_scaled_jacobian(neg[0], neg[1], neg[2], neg[3]);
    assert!((min_neg + 1.0).abs() < 1e-12, "inverted min {min_neg}");
    let q_neg = compute_tet_quality(neg[0], neg[1], neg[2], neg[3]);
    assert!((q_neg.jacobian - 1.0).abs() < 1e-12);
}

#[test]
fn right_corner_tet_matches_the_verdict_scaled_jacobian() {
    let tet = right_corner_tet();
    // Verdict: sqrt(2) * det / max(l_i l_j l_k) = sqrt(2) * 1 / 2.
    let expected = std::f64::consts::FRAC_1_SQRT_2;
    let min_corner = tet_min_scaled_jacobian(tet[0], tet[1], tet[2], tet[3]);
    assert!(
        (min_corner - expected).abs() < 1e-12,
        "min scaled Jacobian {min_corner} != {expected}"
    );
    // The orthogonal corner alone would score sqrt(2); the documented range is
    // [-1, 1], so it saturates at exactly 1.
    let corner0 = tet_scaled_jacobian(tet[0], tet[1], tet[2], tet[3]);
    assert!((corner0 - 1.0).abs() < 1e-12, "corner-0 value {corner0}");
    // The element quality compares |det| with the regular tetrahedron of the same
    // mean edge length: sqrt(2) / mean_edge^3.
    let mean_edge = (3.0 + 3.0 * std::f64::consts::SQRT_2) / 6.0;
    let q = compute_tet_quality(tet[0], tet[1], tet[2], tet[3]);
    let expected_quality = std::f64::consts::SQRT_2 / mean_edge.powi(3);
    assert!(
        (q.jacobian - expected_quality).abs() < 1e-12,
        "jacobian {} != {expected_quality}",
        q.jacobian
    );
}

#[test]
fn tet_quality_metrics_are_scale_invariant_and_bounded() {
    let mut rng = Lcg(0x5EED_1234_ABCD_0001);
    let mut non_degenerate = 0usize;
    for _ in 0..500 {
        let tet = [rng.point(), rng.point(), rng.point(), rng.point()];
        let q = compute_tet_quality(tet[0], tet[1], tet[2], tet[3]);
        let min_corner = tet_min_scaled_jacobian(tet[0], tet[1], tet[2], tet[3]);
        let single = tet_scaled_jacobian(tet[0], tet[1], tet[2], tet[3]);
        assert!(
            (0.0..=1.0).contains(&q.jacobian),
            "jacobian {} outside [0, 1]",
            q.jacobian
        );
        assert!(
            (-1.0..=1.0).contains(&min_corner) && (-1.0..=1.0).contains(&single),
            "scaled Jacobian outside [-1, 1]: {single}, {min_corner}"
        );
        assert!(
            min_corner <= single + 1e-12,
            "min over corners exceeds corner 0"
        );

        let scale = 37.5;
        let big = tet.map(|v| [v[0] * scale, v[1] * scale, v[2] * scale]);
        let qb = compute_tet_quality(big[0], big[1], big[2], big[3]);
        assert!((q.jacobian - qb.jacobian).abs() < 1e-9);
        assert!((q.min_angle_deg - qb.min_angle_deg).abs() < 1e-7);
        assert!((q.max_angle_deg - qb.max_angle_deg).abs() < 1e-7);
        if q.jacobian > 1e-3 {
            non_degenerate += 1;
            // A non-degenerate tet has dihedral angles strictly inside (0, 180).
            assert!(q.min_angle_deg > 0.0 && q.max_angle_deg < 180.0);
        }
    }
    assert!(non_degenerate > 300, "only {non_degenerate} usable samples");
}

#[test]
fn quality_guided_smoothing_never_lowers_the_worst_tet() {
    // A unit cube split into six Kuhn tetrahedra (every one of them has a
    // right-angle corner, which saturates the single-corner scaled Jacobian) whose
    // vertices are perturbed.  The guard compares the minimum over all four
    // corners of every affected tetrahedron, so the global minimum can only grow.
    let mut rng = Lcg(0xC0FF_EE00_1234_5678);
    let corners: [[f64; 3]; 8] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
    ];
    let tets: Vec<[usize; 4]> = vec![
        [0, 1, 3, 7],
        [0, 1, 5, 7],
        [0, 2, 3, 7],
        [0, 2, 6, 7],
        [0, 4, 5, 7],
        [0, 4, 6, 7],
    ];
    let worst = |positions: &[[f64; 3]]| -> f64 {
        tets.iter()
            .map(|t| {
                tet_min_scaled_jacobian(
                    positions[t[0]],
                    positions[t[1]],
                    positions[t[2]],
                    positions[t[3]],
                )
                .abs()
            })
            .fold(f64::INFINITY, f64::min)
    };
    for _ in 0..25 {
        let mut positions: Vec<[f64; 3]> = corners
            .iter()
            .map(|c| {
                [
                    c[0] + 0.08 * rng.next_unit(),
                    c[1] + 0.08 * rng.next_unit(),
                    c[2] + 0.08 * rng.next_unit(),
                ]
            })
            .collect();
        let before = worst(&positions);
        quality_guided_smoothing(&mut positions, &tets, 3);
        let after = worst(&positions);
        assert!(
            after >= before - 1e-12,
            "worst tet degraded: {before} -> {after}"
        );
        assert!(positions.iter().flatten().all(|c| c.is_finite()));
    }
}
