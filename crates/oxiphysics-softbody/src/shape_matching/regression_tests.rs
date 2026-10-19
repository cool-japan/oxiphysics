// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the plastic flow of the rest shape
//! ([`ShapeMatchingBody`], [`PlasticityImprovedBody`]) and for the least-squares
//! quadratic fit of [`QuadraticDeformationBody`].

use super::*;

const TIGHT: f64 = 1.0e-9;

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    v3_norm(v3_sub(a, b))
}

/// Rotation by `angle` about the (not necessarily unit) `axis` (Rodrigues).
fn rotation_about(axis: [f64; 3], angle: f64) -> Mat3x3 {
    let k = v3_scale(axis, 1.0 / v3_norm(axis));
    let (s, c) = angle.sin_cos();
    let cross = [[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]];
    let mut r = mat3_identity();
    for (i, row) in r.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = c * *entry + s * cross[i][j] + (1.0 - c) * k[i] * k[j];
        }
    }
    r
}

fn rotate_z_quarter() -> Mat3x3 {
    rotation_about([0.0, 0.0, 1.0], std::f64::consts::FRAC_PI_2)
}

/// Six particles at `±x`, `±y`, `±z`.
fn octahedron() -> Vec<[f64; 3]> {
    vec![
        [1.0, 0.0, 0.0],
        [-1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, -1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ]
}

/// The octahedron stretched by 1.5 along `x`, then rotated by `rotation`.
fn stretched_octahedron(rotation: Mat3x3) -> Vec<[f64; 3]> {
    octahedron()
        .iter()
        .map(|p| mat3_mul_vec(rotation, [1.5 * p[0], p[1], p[2]]))
        .collect()
}

/// 3×3×3 grid with coordinates in `{-1, 0, 1}`: centred, and rich enough to
/// determine all nine quadratic modes.
fn grid27() -> Vec<[f64; 3]> {
    let mut points = Vec::new();
    for x in [-1.0, 0.0, 1.0] {
        for y in [-1.0, 0.0, 1.0] {
            for z in [-1.0, 0.0, 1.0] {
                points.push([x, y, z]);
            }
        }
    }
    points
}

/// Zero-mean quadratic deformation `x' = x + 0.2 (y² − z²)`: it is exactly
/// representable by the nine quadratic features and keeps the centre of mass.
fn bent(points: &[[f64; 3]]) -> Vec<[f64; 3]> {
    points
        .iter()
        .map(|p| [p[0] + 0.2 * (p[1] * p[1] - p[2] * p[2]), p[1], p[2]])
        .collect()
}

fn rigidly_moved(points: &[[f64; 3]], rotation: Mat3x3, offset: [f64; 3]) -> Vec<[f64; 3]> {
    points
        .iter()
        .map(|p| v3_add(mat3_mul_vec(rotation, *p), offset))
        .collect()
}

// ---------------------------------------------------------------------------
// ShapeMatchingBody plasticity
// ---------------------------------------------------------------------------

/// Sum over the particles of the distance to their goal position.
fn goal_error(body: &mut ShapeMatchingBody) -> f64 {
    let goals = body.goal_positions();
    goals
        .iter()
        .zip(body.positions.iter())
        .map(|(g, p)| dist3(*g, *p))
        .sum()
}

fn plastic_body(rotation: Mat3x3, threshold: f64, rate: f64) -> ShapeMatchingBody {
    let mut body = ShapeMatchingBody::new(octahedron(), vec![1.0; 6], 0.0);
    body.positions = stretched_octahedron(rotation);
    body.plasticity_threshold = threshold;
    body.plasticity_rate = rate;
    body
}

fn assert_error_decreases_monotonically(rotation: Mat3x3) {
    let mut body = plastic_body(rotation, 1.0e-3, 0.3);
    let initial = goal_error(&mut body);
    assert!(
        (initial - 1.0).abs() < TIGHT,
        "two particles are 0.5 off their goals, total {initial}"
    );
    let mut previous = initial;
    for call in 1..=20 {
        body.apply_shape_matching_force(0.0, 0.016);
        let error = goal_error(&mut body);
        assert!(
            error < previous - 1.0e-6,
            "call {call}: the goal error must shrink, {previous} -> {error}"
        );
        previous = error;
    }
    // Geometric decay towards the yield threshold: 2 * (1e-3 + 0.499 * 0.7^20).
    assert!(previous < 0.01, "after 20 calls the error is {previous}");
    // The positions were held (stiffness 0), only the rest shape moved.
    let held = stretched_octahedron(rotation);
    for (now, before) in body.positions.iter().zip(held.iter()) {
        assert!(dist3(*now, *before) < TIGHT);
    }
}

#[test]
fn plasticity_moves_the_goal_towards_the_deformation() {
    assert_error_decreases_monotonically(mat3_identity());
}

#[test]
fn plasticity_works_in_the_rest_frame_of_a_rotated_body() {
    // Rotated by 90 degrees about z the deviation of the +x particle is along
    // +y in the world frame. Adding it to the rest position as it is would
    // move the goal sideways instead of closer; Rᵀ must be applied.
    assert_error_decreases_monotonically(rotate_z_quarter());
}

#[test]
fn plasticity_works_for_a_general_rotation() {
    assert_error_decreases_monotonically(rotation_about([1.0, 2.0, 3.0], 0.9));
}

#[test]
fn plasticity_absorbs_exactly_the_excess_error_at_unit_rate() {
    let threshold = 0.05;
    let mut body = plastic_body(mat3_identity(), threshold, 1.0);
    body.apply_shape_matching_force(0.0, 0.016);
    let goals = body.goal_positions();
    // The two stretched particles end up exactly `threshold` away from their
    // goals, the others never were off.
    assert!((dist3(goals[0], body.positions[0]) - threshold).abs() < TIGHT);
    assert!((dist3(goals[1], body.positions[1]) - threshold).abs() < TIGHT);
    for (i, goal) in goals.iter().enumerate().skip(2) {
        assert!(dist3(*goal, body.positions[i]) < TIGHT, "particle {i}");
    }
    assert!((body.accumulated_plastic_deformation - 2.0 * (0.5 - threshold)).abs() < TIGHT);
}

#[test]
fn plasticity_rate_above_one_does_not_overshoot() {
    let mut unit = plastic_body(mat3_identity(), 0.05, 1.0);
    let mut excessive = plastic_body(mat3_identity(), 0.05, 5.0);
    unit.apply_shape_matching_force(0.0, 0.016);
    excessive.apply_shape_matching_force(0.0, 0.016);
    for (a, b) in unit
        .rest_positions
        .iter()
        .zip(excessive.rest_positions.iter())
    {
        assert!(
            dist3(*a, *b) < TIGHT,
            "rate is clamped to 1: {a:?} vs {b:?}"
        );
    }
    assert_eq!(
        unit.accumulated_plastic_deformation,
        excessive.accumulated_plastic_deformation
    );
}

#[test]
fn plasticity_ignores_non_positive_and_nan_rates() {
    for rate in [0.0, -0.5, f64::NAN] {
        let mut body = plastic_body(mat3_identity(), 0.05, rate);
        let before = body.rest_positions.clone();
        body.apply_shape_matching_force(0.0, 0.016);
        assert_eq!(body.rest_positions, before, "rate {rate}");
        assert_eq!(body.accumulated_plastic_deformation, 0.0);
    }
}

#[test]
fn no_plastic_flow_below_the_threshold() {
    let mut body = plastic_body(mat3_identity(), 0.6, 1.0);
    let before = body.rest_positions.clone();
    body.apply_shape_matching_force(0.0, 0.016);
    assert_eq!(body.rest_positions, before, "errors of 0.5 stay elastic");
    assert_eq!(body.accumulated_plastic_deformation, 0.0);
}

#[test]
fn plastic_flow_keeps_the_rest_shape_centred() {
    // A single stretched particle: the individual rest-position moves are not
    // balanced, the re-centring must restore Σ m q = 0 after every step.
    let masses = vec![1.0, 2.0, 1.0, 3.0, 1.0, 0.5];
    let mut body = ShapeMatchingBody::new(octahedron(), masses.clone(), 0.0);
    body.plasticity_threshold = 0.01;
    body.plasticity_rate = 0.4;
    body.positions[0] = [3.0, 0.0, 0.0];
    let mut previous_goal_distance = f64::MAX;
    for step in 0..15 {
        let goals = body.goal_positions();
        let distance = dist3(goals[0], body.positions[0]);
        assert!(
            distance < previous_goal_distance,
            "step {step}: the stretched particle's goal must approach it, {previous_goal_distance} -> {distance}"
        );
        previous_goal_distance = distance;
        body.apply_shape_matching_force(0.0, 0.016);
        let moment = body
            .rest_positions
            .iter()
            .zip(masses.iter())
            .fold([0.0_f64; 3], |acc, (q, &m)| v3_add(acc, v3_scale(*q, m)));
        assert!(
            v3_norm(moment) < 1.0e-12,
            "step {step}: Σ m q must vanish, got {moment:?}"
        );
    }
}

#[test]
fn elastic_pull_and_plastic_flow_settle_on_the_same_shape() {
    let mut body = plastic_body(rotation_about([0.3, -1.0, 0.5], 1.2), 1.0e-3, 0.2);
    body.stiffness = 0.3;
    for _ in 0..300 {
        body.apply_shape_matching_force(0.3, 0.016);
    }
    let error = goal_error(&mut body);
    assert!(error < 6.0 * 1.0e-3, "residual goal error {error}");
    assert!(body.accumulated_plastic_deformation > 0.1);
    // The relaxed shape has been absorbed: it is no longer the octahedron.
    let stretched_extent = body.rest_positions[0][0]
        .abs()
        .max(body.rest_positions[0][1].abs());
    assert!(
        stretched_extent > 1.05,
        "the rest shape must have taken part of the stretch, extent {stretched_extent}"
    );
}

// ---------------------------------------------------------------------------
// PlasticityImprovedBody
// ---------------------------------------------------------------------------

fn improved_goal_error(body: &PlasticityImprovedBody) -> f64 {
    let xcm = body.centre_of_mass();
    let mut a = mat3_zero();
    for ((pos, q), &mi) in body
        .positions
        .iter()
        .zip(body.plastic_rest.iter())
        .zip(body.masses.iter())
    {
        a = mat3_add_scaled(a, outer_product_f64(v3_sub(*pos, xcm), *q), mi);
    }
    let r = polar_rotation_f64(a);
    body.positions
        .iter()
        .zip(body.plastic_rest.iter())
        .map(|(p, q)| dist3(v3_add(mat3_mul_vec(r, *q), xcm), *p))
        .sum()
}

fn improved_body(rotation: Mat3x3) -> PlasticityImprovedBody {
    let mut body = PlasticityImprovedBody::new(octahedron(), vec![1.0; 6], 0.0);
    body.positions = stretched_octahedron(rotation);
    body.yield_stress = 1.0e-3;
    body.plasticity_rate = 0.3;
    body
}

#[test]
fn improved_plasticity_moves_the_goal_towards_the_deformation() {
    for rotation in [
        mat3_identity(),
        rotate_z_quarter(),
        rotation_about([1.0, 2.0, 3.0], 0.9),
    ] {
        let mut body = improved_body(rotation);
        let mut previous = improved_goal_error(&body);
        assert!((previous - 1.0).abs() < TIGHT);
        for call in 1..=20 {
            body.step(0.016);
            let error = improved_goal_error(&body);
            assert!(
                error < previous - 1.0e-6,
                "call {call}: the goal error must shrink, {previous} -> {error}"
            );
            previous = error;
        }
        assert!(previous < 0.01, "after 20 calls the error is {previous}");
        assert!(body.total_plastic_strain() > 0.9);
    }
}

#[test]
fn improved_plasticity_keeps_the_rest_shape_centred() {
    let masses = vec![1.0, 2.0, 1.0, 3.0, 1.0, 0.5];
    let mut body = PlasticityImprovedBody::new(octahedron(), masses.clone(), 0.0);
    body.yield_stress = 0.01;
    body.plasticity_rate = 0.4;
    body.positions[0] = [3.0, 0.0, 0.0];
    for step in 0..15 {
        body.step(0.016);
        let moment = body
            .plastic_rest
            .iter()
            .zip(masses.iter())
            .fold([0.0_f64; 3], |acc, (q, &m)| v3_add(acc, v3_scale(*q, m)));
        assert!(
            v3_norm(moment) < 1.0e-12,
            "step {step}: Σ m q must vanish, got {moment:?}"
        );
    }
}

#[test]
fn improved_plasticity_hardening_raises_the_yield_and_stops_the_flow_earlier() {
    let mut soft = improved_body(mat3_identity());
    let mut hard = improved_body(mat3_identity());
    hard.hardening_modulus = 1.0;
    for _ in 0..40 {
        soft.step(0.016);
        hard.step(0.016);
    }
    let soft_error = improved_goal_error(&soft);
    let hard_error = improved_goal_error(&hard);
    assert!(
        hard_error > 4.0 * soft_error,
        "hardening leaves a larger elastic residual: {hard_error} vs {soft_error}"
    );
    assert!(hard.total_plastic_strain() < soft.total_plastic_strain());
}

#[test]
fn improved_plasticity_below_yield_changes_nothing() {
    let mut body = improved_body(mat3_identity());
    body.yield_stress = 0.6;
    let before = body.plastic_rest.clone();
    body.step(0.016);
    assert_eq!(body.plastic_rest, before);
    assert_eq!(body.total_plastic_strain(), 0.0);
}

// ---------------------------------------------------------------------------
// Pseudo-inverse of the 9x9 moment matrix
// ---------------------------------------------------------------------------

/// Small deterministic generator of values in `[-1, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_value(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    }
}

/// `B Bᵀ` for a random 9×`rank` matrix `B`: symmetric PSD of rank `rank`.
fn random_psd(rank: usize, seed: u64) -> Mat9 {
    let mut rng = Lcg(seed);
    let columns: Vec<[f64; 9]> = (0..rank)
        .map(|_| {
            let mut column = [0.0; 9];
            for entry in column.iter_mut() {
                *entry = rng.next_value();
            }
            column
        })
        .collect();
    let mut m = [[0.0; 9]; 9];
    for column in &columns {
        for (i, row) in m.iter_mut().enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                *entry += column[i] * column[j];
            }
        }
    }
    m
}

fn mat9_mul(a: &Mat9, b: &Mat9) -> Mat9 {
    let mut c = [[0.0; 9]; 9];
    for (i, row) in c.iter_mut().enumerate() {
        for (j, entry) in row.iter_mut().enumerate() {
            *entry = (0..9).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

fn mat9_max_diff(a: &Mat9, b: &Mat9) -> f64 {
    a.iter()
        .flatten()
        .zip(b.iter().flatten())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

#[test]
fn pseudo_inverse_of_a_full_rank_matrix_is_the_inverse() {
    for seed in [1, 2, 3] {
        let m = random_psd(12, seed);
        let inverse = symmetric_pseudo_inverse(&m);
        let mut identity = [[0.0; 9]; 9];
        for (i, row) in identity.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        let product = mat9_mul(&m, &inverse);
        assert!(
            mat9_max_diff(&product, &identity) < 1.0e-9,
            "seed {seed}: M M⁺ must be the identity"
        );
    }
}

#[test]
fn pseudo_inverse_of_a_rank_deficient_matrix_satisfies_the_penrose_conditions() {
    for (rank, seed) in [(1, 11), (4, 12), (8, 13)] {
        let m = random_psd(rank, seed);
        let p = symmetric_pseudo_inverse(&m);
        let mpm = mat9_mul(&mat9_mul(&m, &p), &m);
        assert!(mat9_max_diff(&mpm, &m) < 1.0e-9, "rank {rank}: M M⁺ M = M");
        let pmp = mat9_mul(&mat9_mul(&p, &m), &p);
        assert!(
            mat9_max_diff(&pmp, &p) < 1.0e-8,
            "rank {rank}: M⁺ M M⁺ = M⁺"
        );
        for (i, row) in p.iter().enumerate() {
            for (j, &value) in row.iter().enumerate() {
                assert!(
                    (value - p[j][i]).abs() < 1.0e-12,
                    "rank {rank}: M⁺ must be symmetric"
                );
            }
        }
    }
}

#[test]
fn pseudo_inverse_of_the_zero_matrix_is_zero() {
    let zero = [[0.0; 9]; 9];
    assert_eq!(symmetric_pseudo_inverse(&zero), zero);
}

// ---------------------------------------------------------------------------
// QuadraticDeformationBody
// ---------------------------------------------------------------------------

fn quadratic_body(rest: &[[f64; 3]], beta: f64) -> QuadraticDeformationBody {
    let mut body = QuadraticDeformationBody::new(rest.to_vec(), vec![1.0; rest.len()], 1.0);
    body.beta = beta;
    body
}

fn assert_goals_equal_positions(body: &QuadraticDeformationBody, tolerance: f64, label: &str) {
    let goals = body.goal_positions();
    assert_eq!(goals.len(), body.positions.len());
    for (i, (goal, position)) in goals.iter().zip(body.positions.iter()).enumerate() {
        assert!(
            dist3(*goal, *position) < tolerance,
            "{label}: particle {i}: goal {goal:?} vs position {position:?}"
        );
    }
}

/// Rest shapes of different rank in the quadratic feature space: the five
/// particles of the old test, a box corner set (`qx² = qy² = qz²`), a planar
/// sheet and the full grid.
fn rest_shapes() -> Vec<(&'static str, Vec<[f64; 3]>)> {
    let corners = vec![
        [1.0, 1.0, 1.0],
        [1.0, 1.0, -1.0],
        [1.0, -1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, 1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [-1.0, -1.0, -1.0],
    ];
    let sheet = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 1.0, 0.0],
        [0.0, 2.0, 0.0],
        [1.0, 2.0, 0.0],
        [2.5, 2.0, 0.0],
    ];
    vec![
        (
            "five particles",
            vec![
                [1.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, -1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
        ),
        ("box corners", corners),
        ("planar sheet", sheet),
        ("grid", grid27()),
    ]
}

#[test]
fn undeformed_quadratic_body_sits_at_its_goals_for_every_beta() {
    for (name, rest) in rest_shapes() {
        for beta in [0.0, 0.25, 0.5, 1.0] {
            let body = quadratic_body(&rest, beta);
            assert_goals_equal_positions(&body, TIGHT, &format!("{name}, beta {beta}"));
        }
    }
}

#[test]
fn rigidly_moved_quadratic_body_sits_at_its_goals_for_every_beta() {
    let rotation = rotation_about([1.0, -2.0, 0.5], 1.1);
    let offset = [0.3, -2.0, 1.1];
    for (name, rest) in rest_shapes() {
        for beta in [0.0, 0.5, 1.0] {
            let mut body = quadratic_body(&rest, beta);
            body.positions = rigidly_moved(&rest, rotation, offset);
            assert_goals_equal_positions(&body, TIGHT, &format!("{name}, beta {beta}"));
        }
    }
}

#[test]
fn quadratic_goals_do_not_depend_on_mass_or_length_units() {
    // Raw moments scale with mass * length²: the old goals did too. The fit
    // is invariant to both, and equivariant to a uniform scaling.
    let rest = grid27();
    let deformed = bent(&rest);
    let reference = {
        let mut body = quadratic_body(&rest, 0.5);
        body.positions = deformed.clone();
        body.goal_positions()
    };

    let mut heavy = QuadraticDeformationBody::new(rest.clone(), vec![1000.0; 27], 1.0);
    heavy.beta = 0.5;
    heavy.positions = deformed.clone();
    for (a, b) in heavy.goal_positions().iter().zip(reference.iter()) {
        assert!(dist3(*a, *b) < TIGHT, "mass independence: {a:?} vs {b:?}");
    }

    for scale in [1.0e-4, 1.0e3] {
        let scaled_rest: Vec<[f64; 3]> = rest.iter().map(|p| v3_scale(*p, scale)).collect();
        let mut body = quadratic_body(&scaled_rest, 0.5);
        body.positions = deformed.iter().map(|p| v3_scale(*p, scale)).collect();
        for (goal, expected) in body.goal_positions().iter().zip(reference.iter()) {
            let want = v3_scale(*expected, scale);
            assert!(
                dist3(*goal, want) < TIGHT * scale,
                "scale {scale}: {goal:?} vs {want:?}"
            );
        }
        // An undeformed body of any size is at its goals.
        let at_rest = quadratic_body(&scaled_rest, 1.0);
        assert_goals_equal_positions(&at_rest, TIGHT * scale, &format!("scale {scale}"));
    }
}

#[test]
fn quadratic_fit_reproduces_an_exactly_quadratic_deformation() {
    let rest = grid27();
    let mut body = quadratic_body(&rest, 1.0);
    body.positions = bent(&rest);
    assert_goals_equal_positions(&body, TIGHT, "bent grid, beta 1");

    // The deformation matrix itself: p = Ã · [qx, qy, qz, qx², qy², qz², qxqy, qyqz, qzqx].
    let a = body.compute_a39();
    let expected = [
        [1.0, 0.0, 0.0, 0.0, 0.2, -0.2, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    for (row, (got, want)) in a.iter().zip(expected.iter()).enumerate() {
        for (col, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
            assert!((g - w).abs() < TIGHT, "Ã[{row}][{col}] = {g}, expected {w}");
        }
    }
}

#[test]
fn quadratic_deformation_matrix_of_the_undeformed_grid_is_the_linear_identity() {
    // Full-rank rest shape: Ã = [I 0 0] exactly.
    let a = quadratic_body(&grid27(), 1.0).compute_a39();
    for (row, values) in a.iter().enumerate() {
        for (col, &value) in values.iter().enumerate() {
            let expected = if row == col { 1.0 } else { 0.0 };
            assert!(
                (value - expected).abs() < TIGHT,
                "Ã[{row}][{col}] = {value}, expected {expected}"
            );
        }
    }
}

#[test]
fn quadratic_deformation_matrix_maps_every_rest_feature_onto_its_rest_offset() {
    // Also for rank-deficient rest shapes, where Ã is only determined on the
    // span of the features that occur.
    for (name, rest) in rest_shapes() {
        let body = quadratic_body(&rest, 1.0);
        let a = body.compute_a39();
        let centre = rest.iter().fold([0.0; 3], |acc, p| {
            v3_add(acc, v3_scale(*p, 1.0 / rest.len() as f64))
        });
        for p in &rest {
            let offset = v3_sub(*p, centre);
            let mapped = mul_3x9_vec(&a, &QuadraticDeformationBody::feature(offset));
            assert!(
                dist3(mapped, offset) < 1.0e-7,
                "{name}: Ã · feature(q) = {mapped:?} must reproduce q = {offset:?}"
            );
        }
    }
}

#[test]
fn quadratic_beta_zero_is_rigid_shape_matching() {
    let rest = grid27();
    let deformed = bent(&rest);
    let mut quadratic = quadratic_body(&rest, 0.0);
    quadratic.positions = deformed.clone();
    let mut rigid = ShapeMatchingBody::new(rest.clone(), vec![1.0; rest.len()], 1.0);
    rigid.positions = deformed;
    let rigid_goals = rigid.goal_positions();
    let quadratic_goals = quadratic.goal_positions();
    let mut max_residual = 0.0_f64;
    for (q, r) in quadratic_goals.iter().zip(rigid_goals.iter()) {
        assert!(dist3(*q, *r) < TIGHT, "{q:?} vs {r:?}");
    }
    for (g, p) in quadratic_goals.iter().zip(quadratic.positions.iter()) {
        max_residual = max_residual.max(dist3(*g, *p));
    }
    assert!(
        max_residual > 0.15,
        "a rigid goal cannot follow the bend, residual {max_residual}"
    );
}

#[test]
fn quadratic_goals_blend_linearly_between_rigid_and_quadratic() {
    let rest = grid27();
    let deformed = bent(&rest);
    let goals_at = |beta: f64| {
        let mut body = quadratic_body(&rest, beta);
        body.positions = deformed.clone();
        body.goal_positions()
    };
    let rigid = goals_at(0.0);
    let full = goals_at(1.0);
    let mid = goals_at(0.3);
    for ((r, f), m) in rigid.iter().zip(full.iter()).zip(mid.iter()) {
        let expected = v3_add(v3_scale(*f, 0.3), v3_scale(*r, 0.7));
        assert!(dist3(*m, expected) < TIGHT);
    }
    // beta = 1 follows the bend, beta = 0 does not.
    let body = quadratic_body(&rest, 1.0);
    assert_eq!(body.beta, 1.0);
    for (g, p) in full.iter().zip(deformed.iter()) {
        assert!(dist3(*g, *p) < TIGHT);
    }
}

#[test]
fn quadratic_body_defaults_to_an_even_blend() {
    let body = QuadraticDeformationBody::new(grid27(), vec![1.0; 27], 0.7);
    assert_eq!(body.beta, 0.5);
    assert_eq!(body.stiffness, 0.7);
}

#[test]
fn quadratic_correction_relaxes_a_bent_body_into_its_rigid_goal() {
    let rest = grid27();
    let mut body = quadratic_body(&rest, 0.0);
    body.positions = bent(&rest);
    body.apply_correction(1.0, 0.016);
    // With beta = 0 and stiffness 1 the body jumps to the rigid image of the
    // rest shape, where it then sits at its goals.
    assert_goals_equal_positions(&body, TIGHT, "after the correction");
    let mut rigid = ShapeMatchingBody::new(rest.clone(), vec![1.0; rest.len()], 1.0);
    rigid.positions = bent(&rest);
    let expected = rigid.goal_positions();
    for (p, e) in body.positions.iter().zip(expected.iter()) {
        assert!(dist3(*p, *e) < TIGHT);
    }
}

#[test]
fn quadratic_body_with_coincident_particles_stays_finite() {
    let body = QuadraticDeformationBody::new(vec![[1.0, 2.0, 3.0]; 4], vec![1.0; 4], 1.0);
    for goal in body.goal_positions() {
        for component in goal {
            assert!(component.is_finite());
        }
        assert!(dist3(goal, [1.0, 2.0, 3.0]) < TIGHT);
    }
}
