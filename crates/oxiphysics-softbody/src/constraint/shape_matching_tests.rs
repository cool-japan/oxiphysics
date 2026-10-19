// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for [`ShapeMatchingConstraint`].
//!
//! The old rotation orthogonalised the columns of the cross-covariance matrix
//! `A = sum m p q^T` (Gram-Schmidt, i.e. the QR factor of `A`) instead of
//! extracting the rotation of its polar decomposition.  At rest `A` is the
//! symmetric rest covariance matrix, whose QR factor is a rotation unless the
//! matrix is diagonal, so resting clusters were pulled off their rest shape; and
//! for a planar group `A` has rank 2, the orthogonalisation broke down and the
//! identity was returned, snapping the goals back to the rest orientation.
//!
//! The expected values follow from the definition of the polar factor: for
//! `A = R S` with `S` symmetric positive semi-definite the rotation `R` is
//! recovered exactly, and it is the proper rotation closest to `A`.

use super::functions::SoftConstraint;
use super::types::{ShapeMatchingConstraint, polar_rotation};
use crate::particle::SoftParticle;
use oxiphysics_core::math::{Mat3, Quat, Real, Unit, Vec3};
use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

const DT: Real = 1.0 / 60.0;

/// Deterministic 64-bit linear congruential generator for randomised tests.
struct Lcg(u64);

impl Lcg {
    /// Uniform sample in `[0, 1)`.
    fn next_unit(&mut self) -> Real {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as Real / (1_u64 << 53) as Real
    }

    /// Uniform sample in `[-1, 1)`.
    fn next_signed(&mut self) -> Real {
        2.0 * self.next_unit() - 1.0
    }

    /// Unit vector with a random direction.
    fn next_axis(&mut self) -> Vec3 {
        loop {
            let candidate = Vec3::new(self.next_signed(), self.next_signed(), self.next_signed());
            if candidate.norm() > 0.2 {
                return candidate.normalize();
            }
        }
    }

    /// Rotation about a random axis by a random angle.
    fn next_rotation(&mut self) -> Mat3 {
        let axis = Unit::new_normalize(self.next_axis());
        let angle = PI * self.next_signed();
        Quat::from_axis_angle(&axis, angle)
            .to_rotation_matrix()
            .into_inner()
    }

    /// Matrix with independent entries in `[-1, 1)`.
    fn next_matrix(&mut self) -> Mat3 {
        Mat3::from_fn(|_, _| self.next_signed())
    }
}

/// Rotation by `angle` (radians) about `axis`.
fn rotation(axis: [Real; 3], angle: Real) -> Mat3 {
    let axis = Unit::new_normalize(Vec3::new(axis[0], axis[1], axis[2]));
    Quat::from_axis_angle(&axis, angle)
        .to_rotation_matrix()
        .into_inner()
}

/// Rotations that exercise the extraction: about the normal of the `z = 0` plane,
/// about axes inside that plane (these turn a planar group over), half turns and a
/// generic one.
fn rotations() -> Vec<Mat3> {
    vec![
        rotation([0.0, 0.0, 1.0], FRAC_PI_2),
        rotation([1.0, 1.0, 0.0], 70.0_f64.to_radians()),
        rotation([1.0, 0.0, 0.0], PI),
        rotation([1.0, 1.0, 0.0], PI),
        rotation([1.0, 2.0, 3.0], 2.5),
        rotation([0.3, -1.0, 0.5], 1.0e-3),
    ]
}

fn points_of(rows: &[[Real; 3]]) -> Vec<Vec3> {
    rows.iter().map(|r| Vec3::new(r[0], r[1], r[2])).collect()
}

/// A named point set with one mass per point.
struct Cluster {
    name: &'static str,
    points: Vec<Vec3>,
    masses: Vec<Real>,
}

impl Cluster {
    fn rest_particles(&self) -> Vec<SoftParticle> {
        self.mapped(&Mat3::identity(), Vec3::zeros())
    }

    /// The cluster after the affine map `x -> map * x + shift`.
    fn mapped(&self, map: &Mat3, shift: Vec3) -> Vec<SoftParticle> {
        self.points
            .iter()
            .zip(&self.masses)
            .map(|(p, &m)| SoftParticle::new(map * p + shift, m))
            .collect()
    }

    /// A constraint built from the rest configuration.
    fn constraint(&self, stiffness: Real) -> ShapeMatchingConstraint {
        ShapeMatchingConstraint::from_particles(
            (0..self.points.len()).collect(),
            &self.rest_particles(),
            stiffness,
        )
    }
}

/// Corners of the cube `[-1, 1]^3`: an isotropic cluster, `sum q q^T = 8 I`.
fn cube() -> Cluster {
    let mut points = Vec::new();
    for x in [-1.0, 1.0] {
        for y in [-1.0, 1.0] {
            for z in [-1.0, 1.0] {
                points.push(Vec3::new(x, y, z));
            }
        }
    }
    Cluster {
        name: "cube",
        masses: vec![1.0; points.len()],
        points,
    }
}

/// Antipodal pairs `+-L e_j` with `L L^T = [[2, 1, 0], [1, 2, 0], [0, 0, 1]]`, so
/// that the rest covariance is `2 L L^T = [[4, 2, 0], [2, 4, 0], [0, 0, 2]]`: not
/// diagonal, and its QR factor is a rotation by about 26.6 degrees.
fn sheared() -> Cluster {
    let a = 2.0_f64.sqrt();
    let b = 1.5_f64.sqrt();
    let half = [[a, 1.0 / a, 0.0], [0.0, b, 0.0], [0.0, 0.0, 1.0]];
    let mut points = Vec::new();
    for h in half {
        points.push(Vec3::new(h[0], h[1], h[2]));
        points.push(Vec3::new(-h[0], -h[1], -h[2]));
    }
    Cluster {
        name: "sheared",
        masses: vec![1.0; points.len()],
        points,
    }
}

/// Scalene solid cluster with distinct masses.
fn solid() -> Cluster {
    Cluster {
        name: "solid",
        points: points_of(&[
            [0.0, 0.0, 0.0],
            [1.2, 0.1, -0.3],
            [-0.4, 1.1, 0.2],
            [0.3, -0.5, 1.4],
            [-1.0, -0.8, -0.6],
            [0.6, 0.9, 0.8],
            [-0.7, 0.2, 1.1],
        ]),
        masses: vec![1.0, 2.0, 1.5, 0.7, 1.2, 0.9, 1.8],
    }
}

/// Scalene planar cluster in the plane `z = 0` with distinct masses.
fn planar() -> Cluster {
    Cluster {
        name: "planar",
        points: points_of(&[
            [0.0, 0.0, 0.0],
            [1.5, 0.2, 0.0],
            [-0.6, 1.1, 0.0],
            [0.4, -1.2, 0.0],
            [-1.1, -0.5, 0.0],
        ]),
        masses: vec![1.0, 2.0, 1.0, 1.5, 0.8],
    }
}

/// Four corners of the rectangle `[-2, 2] x [-1, 1]` in the plane `z = 0`.
fn rectangle() -> Cluster {
    Cluster {
        name: "rectangle",
        points: points_of(&[
            [-2.0, -1.0, 0.0],
            [2.0, -1.0, 0.0],
            [2.0, 1.0, 0.0],
            [-2.0, 1.0, 0.0],
        ]),
        masses: vec![1.0; 4],
    }
}

/// Collinear cluster along `(1, 2, 3)`.
fn collinear() -> Cluster {
    let direction = Vec3::new(1.0, 2.0, 3.0) / 14.0_f64.sqrt();
    Cluster {
        name: "collinear",
        points: [-1.0, -0.2, 0.7, 1.5]
            .iter()
            .map(|t| direction * *t)
            .collect(),
        masses: vec![1.0, 2.0, 1.0, 0.5],
    }
}

fn solid_clusters() -> Vec<Cluster> {
    vec![cube(), sheared(), solid()]
}

fn planar_clusters() -> Vec<Cluster> {
    vec![planar(), rectangle()]
}

fn positions(particles: &[SoftParticle]) -> Vec<Vec3> {
    particles.iter().map(|p| p.position).collect()
}

fn max_displacement(particles: &[SoftParticle], reference: &[Vec3]) -> Real {
    particles
        .iter()
        .zip(reference)
        .map(|(p, r)| (p.position - r).norm())
        .fold(0.0, Real::max)
}

fn max_entry_difference(a: &Mat3, b: &Mat3) -> Real {
    (a - b).iter().fold(0.0, |largest, e| largest.max(e.abs()))
}

fn mass_weighted_centroid(particles: &[SoftParticle]) -> Vec3 {
    let mut total = 0.0;
    let mut sum = Vec3::zeros();
    for p in particles {
        let mass = 1.0 / p.inverse_mass;
        total += mass;
        sum += p.position * mass;
    }
    sum / total
}

/// Rotation reported by the constraint for the particles in `current`.
fn extracted_rotation(constraint: &ShapeMatchingConstraint, current: &[SoftParticle]) -> Mat3 {
    let com = constraint.current_com(current);
    constraint.compute_rotation(current, &com)
}

/// Rigidly moves `cluster` by every rotation of [`rotations`] and a shift and
/// checks that the rotation is recovered and that nothing needs correcting.
fn assert_rigid_motions_are_recovered(cluster: &Cluster) {
    for r0 in rotations() {
        for shift in [Vec3::zeros(), Vec3::new(-4.0, 2.5, 1.0)] {
            let mut current = cluster.mapped(&r0, shift);
            let before = positions(&current);
            let mut constraint = cluster.constraint(1.0);
            let r = extracted_rotation(&constraint, &current);
            assert!(
                max_entry_difference(&r, &r0) < 1e-9,
                "{}: rotation not recovered, R = {r}, expected {r0}",
                cluster.name
            );
            constraint.project(&mut current, DT);
            assert!(
                max_displacement(&current, &before) < 1e-9,
                "{}: a rigidly moved group must be its own goal, moved by {}",
                cluster.name,
                max_displacement(&current, &before)
            );
        }
    }
}

#[test]
fn resting_clusters_keep_the_identity_rotation_and_are_not_corrected() {
    // The rotation of a collinear group is fixed only up to a twist about its
    // line, so it is checked separately; its goals are checked here as well.
    let mut clusters = solid_clusters();
    clusters.extend(planar_clusters());
    clusters.push(collinear());
    for cluster in clusters {
        for shift in [Vec3::zeros(), Vec3::new(3.0, -2.0, 7.5)] {
            let mut current = cluster.mapped(&Mat3::identity(), shift);
            let before = positions(&current);
            let mut constraint = cluster.constraint(1.0);
            let r = extracted_rotation(&constraint, &current);
            if cluster.name != "collinear" {
                assert!(
                    max_entry_difference(&r, &Mat3::identity()) < 1e-12,
                    "{}: a group at rest must not be rotated, R = {r}",
                    cluster.name
                );
            }
            constraint.project(&mut current, DT);
            assert!(
                max_displacement(&current, &before) < 1e-12,
                "{}: a group at rest must not be corrected",
                cluster.name
            );
        }
    }
}

#[test]
fn rigid_rotations_of_solid_clusters_are_recovered() {
    for cluster in solid_clusters() {
        assert_rigid_motions_are_recovered(&cluster);
    }
}

#[test]
fn rigid_rotations_of_planar_clusters_are_recovered() {
    // The cross-covariance matrix of a planar group has rank 2.  The rotations
    // about axes inside the plane turn the group over, which only the sign choice
    // for the vanishing singular value gets right.
    for cluster in planar_clusters() {
        assert_rigid_motions_are_recovered(&cluster);
    }
}

#[test]
fn rigid_rotations_of_collinear_clusters_leave_the_goals_unchanged() {
    // The rotation of a rank-1 group is only fixed up to a twist about the line,
    // which does not move any goal position.
    let cluster = collinear();
    let mut motions = rotations();
    motions.push(Mat3::identity());
    for r0 in motions {
        let mut current = cluster.mapped(&r0, Vec3::new(1.0, -1.0, 2.0));
        let before = positions(&current);
        let mut constraint = cluster.constraint(1.0);
        let r = extracted_rotation(&constraint, &current);
        let axis = cluster.points[3] - cluster.points[0];
        assert!(
            (r * axis - r0 * axis).norm() < 1e-9,
            "the rotation must carry the line onto its image"
        );
        assert!(
            (r.transpose() * r - Mat3::identity()).norm() < 1e-12
                && (r.determinant() - 1.0).abs() < 1e-12,
            "the rotation must be proper"
        );
        constraint.project(&mut current, DT);
        assert!(max_displacement(&current, &before) < 1e-9);
    }
}

#[test]
fn quarter_turn_of_a_planar_group_about_its_normal_is_not_undone() {
    // The case that used to return the identity: a cloth patch rotated by 90
    // degrees about z must keep its orientation instead of snapping back.
    let cluster = planar();
    let r0 = rotation([0.0, 0.0, 1.0], FRAC_PI_2);
    let mut current = cluster.mapped(&r0, Vec3::new(2.0, 0.0, -1.0));
    let before = positions(&current);
    let mut constraint = cluster.constraint(1.0);
    constraint.project(&mut current, DT);
    assert!(
        max_displacement(&current, &before) < 1e-9,
        "the patch was moved back towards its rest orientation by {}",
        max_displacement(&current, &before)
    );
}

#[test]
fn rotation_with_anisotropic_stretch_recovers_the_exact_rotation() {
    // For an isotropic rest shape (the cube) A = 8 R0 F, so the polar factor is R0
    // for every symmetric positive definite stretch F, however anisotropic and
    // however it is oriented with respect to the cube.
    let cube = cube();
    let tilt = rotation([2.0, -1.0, 1.0], 0.7);
    let stretch = tilt * Mat3::from_diagonal(&Vec3::new(1.4, 0.8, 1.1)) * tilt.transpose();
    let r0 = rotation([1.0, 2.0, 3.0], 2.5);
    let shift = Vec3::new(0.5, -3.0, 2.0);
    let mut current = cube.mapped(&(r0 * stretch), shift);
    let mut constraint = cube.constraint(1.0);
    let r = extracted_rotation(&constraint, &current);
    assert!(
        max_entry_difference(&r, &r0) < 1e-9,
        "R = {r}, expected {r0}"
    );

    // Full stiffness replaces the stretched cube by the rotated rigid cube.
    constraint.project(&mut current, DT);
    let rigid = cube.mapped(&r0, shift);
    assert!(max_displacement(&current, &positions(&rigid)) < 1e-9);

    // The same holds for a planar group stretched along its own axes: the
    // rectangle has the diagonal covariance diag(16, 4, 0), so F S is diagonal
    // too and the rotation is again R0 (and the group is turned over by it).
    let rectangle = rectangle();
    let r0 = rotation([1.0, 1.0, 0.0], 70.0_f64.to_radians());
    let stretch = Mat3::from_diagonal(&Vec3::new(1.5, 0.6, 1.0));
    let mut current = rectangle.mapped(&(r0 * stretch), shift);
    let mut constraint = rectangle.constraint(1.0);
    let r = extracted_rotation(&constraint, &current);
    assert!(
        max_entry_difference(&r, &r0) < 1e-9,
        "R = {r}, expected {r0}"
    );
    constraint.project(&mut current, DT);
    let rigid = rectangle.mapped(&r0, shift);
    assert!(max_displacement(&current, &positions(&rigid)) < 1e-9);
}

#[test]
fn rotation_is_found_for_tiny_and_huge_masses_and_extents() {
    // The cross-covariance matrix of a milligram patch of a millimetre is about
    // 1e-15; an absolute tolerance on it used to discard the orientation.
    let r0 = rotation([1.0, 1.0, 0.0], 70.0_f64.to_radians());
    for (length, mass) in [(1.0e-3, 1.0e-9), (1.0, 1.0), (1.0e3, 1.0e6)] {
        let reference = planar();
        let cluster = Cluster {
            name: "scaled planar",
            points: reference.points.iter().map(|p| p * length).collect(),
            masses: reference.masses.iter().map(|m| m * mass).collect(),
        };
        let shift = Vec3::new(1.0, -2.0, 0.5) * length;
        let current = cluster.mapped(&r0, shift);
        let constraint = cluster.constraint(1.0);
        let r = extracted_rotation(&constraint, &current);
        assert!(
            max_entry_difference(&r, &r0) < 1e-9,
            "length {length}, mass {mass}: R = {r}"
        );
    }
}

#[test]
fn mirrored_planar_group_is_matched_by_turning_it_over() {
    // Mirroring a planar group inside its plane is the same shape seen from the
    // other side: the half turn about x reproduces it exactly.
    let cluster = planar();
    let mirror = Mat3::from_diagonal(&Vec3::new(1.0, -1.0, 1.0));
    let mut current = cluster.mapped(&mirror, Vec3::new(0.0, 1.0, 0.0));
    let before = positions(&current);
    let mut constraint = cluster.constraint(1.0);
    let r = extracted_rotation(&constraint, &current);
    let half_turn = Mat3::from_diagonal(&Vec3::new(1.0, -1.0, -1.0));
    assert!(max_entry_difference(&r, &half_turn) < 1e-9, "R = {r}");
    constraint.project(&mut current, DT);
    assert!(max_displacement(&current, &before) < 1e-9);
}

#[test]
fn mirrored_solid_group_gets_the_closest_proper_rotation() {
    // A mirrored solid group (det A < 0) has no rotation that matches it; the
    // rotation must still be proper and the best fit among all rotations.
    for cluster in [solid(), sheared()] {
        let mirror = Mat3::from_diagonal(&Vec3::new(1.0, 1.0, -1.0));
        let current = cluster.mapped(&mirror, Vec3::zeros());
        let constraint = cluster.constraint(1.0);
        let com = constraint.current_com(&current);
        let r = constraint.compute_rotation(&current, &com);
        assert!(
            (r.transpose() * r - Mat3::identity()).norm() < 1e-12,
            "{}: R is not orthogonal",
            cluster.name
        );
        assert!(
            (r.determinant() - 1.0).abs() < 1e-12,
            "{}: R must be a rotation, det = {}",
            cluster.name,
            r.determinant()
        );

        let mut a = Mat3::zeros();
        for (p, q) in current.iter().zip(&constraint.rest_relative) {
            a += (p.position - com) * q.transpose() * (1.0 / p.inverse_mass);
        }
        assert!(a.determinant() < 0.0, "the scenario must mirror the group");
        let fit = (r.transpose() * a).trace();
        let mut rng = Lcg(0xFEED_0001);
        for _ in 0..500 {
            let rival = rng.next_rotation();
            assert!(
                fit >= (rival.transpose() * a).trace() - 1e-12,
                "{}: a rotation fits better than the extracted one",
                cluster.name
            );
        }
    }
}

#[test]
fn polar_rotation_is_the_closest_proper_rotation_of_random_matrices() {
    let mut rng = Lcg(0x5EED_1234_ABCD_EF01);
    let rivals: Vec<Mat3> = (0..400).map(|_| rng.next_rotation()).collect();
    let mut seen = [0_usize; 2];
    for _ in 0..300 {
        let a = rng.next_matrix();
        let r = polar_rotation(&a);
        assert!((r.transpose() * r - Mat3::identity()).norm() < 1e-12);
        assert!((r.determinant() - 1.0).abs() < 1e-12);

        let fit = (r.transpose() * a).trace();
        for rival in &rivals {
            assert!(fit >= (rival.transpose() * a).trace() - 1e-12);
        }
        // The optimum is s0 + s1 + sign(det A) s2 (Kabsch).
        let mut s: Vec<Real> = a.singular_values().iter().copied().collect();
        s.sort_by(|x, y| y.total_cmp(x));
        let sign = if a.determinant() < 0.0 { -1.0 } else { 1.0 };
        let optimum = s[0] + s[1] + sign * s[2];
        assert!(
            (fit - optimum).abs() < 1e-9 * (1.0 + optimum.abs()),
            "tr(R^T A) = {fit}, optimum {optimum}"
        );
        seen[usize::from(a.determinant() < 0.0)] += 1;

        // Positive scaling does not change the rotation.
        for factor in [1.0e-12, 1.0e12] {
            assert!(max_entry_difference(&polar_rotation(&(a * factor)), &r) < 1e-9);
        }
    }
    assert!(
        seen[0] > 20 && seen[1] > 20,
        "both orientations must be covered: {seen:?}"
    );
}

#[test]
fn polar_rotation_inverts_a_rotation_times_a_symmetric_positive_definite_matrix() {
    let mut rng = Lcg(0xC0FFEE);
    for _ in 0..200 {
        let r0 = rng.next_rotation();
        let b = rng.next_matrix();
        let s = b.transpose() * b + Mat3::identity() * 0.1;
        let r = polar_rotation(&(r0 * s));
        assert!(max_entry_difference(&r, &r0) < 1e-9, "R = {r}, R0 = {r0}");
    }
}

#[test]
fn polar_rotation_of_rank_deficient_matrices_agrees_with_them_on_their_range() {
    let mut rng = Lcg(0xBAD5EED);
    for _ in 0..200 {
        // Rank 2: A = s1 u1 v1^T + s2 u2 v2^T.  The unique proper rotation that
        // agrees with A is u1 v1^T + u2 v2^T + (u1 x u2)(v1 x v2)^T.
        let u = rng.next_rotation();
        let v = rng.next_rotation();
        let (u1, u2) = (u.column(0).into_owned(), u.column(1).into_owned());
        let (v1, v2) = (v.column(0).into_owned(), v.column(1).into_owned());
        let a = u1 * v1.transpose() * 2.0 + u2 * v2.transpose() * 0.5;
        let expected =
            u1 * v1.transpose() + u2 * v2.transpose() + u1.cross(&u2) * v1.cross(&v2).transpose();
        let r = polar_rotation(&a);
        assert!(
            max_entry_difference(&r, &expected) < 1e-9,
            "rank 2: R = {r}"
        );

        // Rank 1: A = s u v^T.  Any proper rotation with R v = u is a polar factor.
        let a = u1 * v1.transpose() * 3.0;
        let r = polar_rotation(&a);
        assert!((r.transpose() * r - Mat3::identity()).norm() < 1e-12);
        assert!((r.determinant() - 1.0).abs() < 1e-12);
        assert!((r * v1 - u1).norm() < 1e-9, "rank 1: R v = {}", r * v1);
    }
}

#[test]
fn polar_rotation_of_unusable_matrices_is_the_identity() {
    assert_eq!(polar_rotation(&Mat3::zeros()), Mat3::identity());
    let mut broken = Mat3::identity();
    broken[(1, 2)] = Real::NAN;
    assert_eq!(polar_rotation(&broken), Mat3::identity());
    broken[(1, 2)] = Real::INFINITY;
    assert_eq!(polar_rotation(&broken), Mat3::identity());
}

#[test]
fn projection_makes_the_group_rigid_and_conserves_the_centroid() {
    let mut clusters = solid_clusters();
    clusters.extend(planar_clusters());
    for cluster in clusters {
        let mut rng = Lcg(0x1234_5678);
        let r0 = rng.next_rotation();
        let shift = Vec3::new(0.3, -0.7, 1.1);
        let mut current = cluster.mapped(&r0, shift);
        for p in &mut current {
            p.position += Vec3::new(rng.next_signed(), rng.next_signed(), rng.next_signed()) * 0.2;
        }
        let centroid_before = mass_weighted_centroid(&current);
        let mut constraint = cluster.constraint(1.0);
        constraint.project(&mut current, DT);

        let centroid_after = mass_weighted_centroid(&current);
        assert!(
            (centroid_after - centroid_before).norm() < 1e-12,
            "{}: the centroid moved by {}",
            cluster.name,
            (centroid_after - centroid_before).norm()
        );
        // Full stiffness puts every particle on its goal R q + com: the group is
        // a rigid copy of the rest shape.
        for (i, first) in current.iter().enumerate() {
            for (j, second) in current.iter().enumerate().take(i) {
                let rest = (cluster.points[i] - cluster.points[j]).norm();
                let now = (first.position - second.position).norm();
                assert!(
                    (rest - now).abs() < 1e-12,
                    "{}: distance {i}-{j} is {now}, rest {rest}",
                    cluster.name
                );
            }
        }
    }
}

#[test]
fn stiffness_scales_the_step_towards_the_goal() {
    let cluster = solid();
    let mut rng = Lcg(99);
    let r0 = rng.next_rotation();
    let mut start = cluster.mapped(&r0, Vec3::new(1.0, 2.0, 3.0));
    for p in &mut start {
        p.position += Vec3::new(rng.next_signed(), rng.next_signed(), rng.next_signed()) * 0.15;
    }
    let initial = positions(&start);
    let mut rigid = start.clone();
    cluster.constraint(1.0).project(&mut rigid, DT);
    for stiffness in [0.0, 0.25, 0.5, 1.0] {
        let mut stepped = start.clone();
        cluster.constraint(stiffness).project(&mut stepped, DT);
        for (k, ((moved, goal), origin)) in stepped.iter().zip(&rigid).zip(&initial).enumerate() {
            let expected = origin + (goal.position - origin) * stiffness;
            assert!(
                (moved.position - expected).norm() < 1e-12,
                "stiffness {stiffness}, particle {k}"
            );
        }
    }
}

#[test]
fn half_stiffness_conserves_the_unweighted_centroid_of_equal_masses() {
    // With equal masses and no pinned particles sum(q) = 0, so the goals have
    // the current centroid as their mean and the step does not move it.
    let cluster = cube();
    let mut rng = Lcg(2024);
    let r0 = rng.next_rotation();
    let mut current = cluster.mapped(&r0, Vec3::new(-2.0, 0.5, 4.0));
    for p in &mut current {
        p.position += Vec3::new(rng.next_signed(), rng.next_signed(), rng.next_signed()) * 0.3;
    }
    let mean = |particles: &[SoftParticle]| {
        particles
            .iter()
            .fold(Vec3::zeros(), |acc, p| acc + p.position)
            / particles.len() as Real
    };
    let before = mean(&current);
    let mut constraint = cluster.constraint(0.5);
    constraint.project(&mut current, DT);
    assert!(
        (mean(&current) - before).norm() < 1e-13,
        "the centroid moved by {}",
        (mean(&current) - before).norm()
    );
}

#[test]
fn pinned_particles_stay_put_and_the_rest_follows_their_orientation() {
    // Pin one particle of the group.  A rigid motion of the whole group, pinned
    // particle included, is still recovered...
    let mut cluster = solid();
    cluster.masses[1] = 0.0;
    let r0 = rotation([1.0, -1.0, 2.0], 1.3);
    let mut current = cluster.mapped(&r0, Vec3::new(1.0, 0.0, -2.0));
    let before = positions(&current);
    let mut constraint = cluster.constraint(1.0);
    assert!(
        max_entry_difference(&extracted_rotation(&constraint, &current), &r0) < 1e-9,
        "the rotation of a rigid motion must be found with a pinned particle"
    );
    constraint.project(&mut current, DT);
    assert!(max_displacement(&current, &before) < 1e-9);

    // ...and when the group is deformed the pinned particle is not moved at all
    // while the others are corrected.
    let mut rng = Lcg(5);
    for p in &mut current {
        p.position += Vec3::new(rng.next_signed(), rng.next_signed(), rng.next_signed()) * 0.1;
    }
    let deformed = positions(&current);
    constraint.project(&mut current, DT);
    assert_eq!(current[1].position, deformed[1], "pinned particle moved");
    let moved: Real = current
        .iter()
        .zip(&deformed)
        .map(|(p, d)| (p.position - d).norm())
        .sum();
    assert!(moved > 1e-3, "the free particles must be corrected");
}

#[test]
fn pinned_planar_patch_swings_down_under_gravity_like_a_rigid_pendulum() {
    use crate::particle::SoftBody;
    use crate::solver::XpbdSolver;

    // A 3 x 3 patch in the plane z = 0 hangs from its corner (0, 2, 0); gravity
    // pulls along -y.  The patch is a planar group, so its cross-covariance matrix
    // has rank 2: with the old rotation (the identity for planar groups) the goals
    // never turned and the patch stayed in its rest orientation.  With the polar
    // factor it swings as a rigid body about the pivot.
    let mut particles = Vec::new();
    for y in 0..3 {
        for x in 0..3 {
            let position = Vec3::new(x as Real, y as Real, 0.0);
            if x == 0 && y == 2 {
                particles.push(SoftParticle::new_static(position));
            } else {
                let mut particle = SoftParticle::new(position, 1.0);
                particle.external_force = Vec3::new(0.0, -9.81, 0.0);
                particles.push(particle);
            }
        }
    }
    let rest = positions(&particles);
    let pivot = Vec3::new(0.0, 2.0, 0.0);
    let group: Vec<usize> = (0..particles.len()).collect();
    let mut constraints: [Box<dyn SoftConstraint>; 1] = [Box::new(
        ShapeMatchingConstraint::from_particles(group, &particles, 1.0),
    )];
    let mut body = SoftBody::from_particles(particles);
    body.damping = 0.02;
    let mut solver = XpbdSolver::new(4);
    // Keep the solver awake: the test checks the settled state, not sleeping.
    solver.sleep_threshold = 0.0;
    for _ in 0..800 {
        solver.solve(&mut body, &mut constraints, DT);
    }

    // At rest the centroid of the free particles, (1.125, 0.875, 0), hangs straight
    // below the pivot: the patch has turned by -45 degrees about z.
    let swing = rotation([0.0, 0.0, 1.0], -FRAC_PI_4);
    for (k, (particle, rest_position)) in body.particles.iter().zip(&rest).enumerate() {
        let expected = pivot + swing * (rest_position - pivot);
        assert!(
            (particle.position - expected).norm() < 1e-3,
            "particle {k} is at {:?}, expected {expected:?}",
            particle.position
        );
    }
    assert_eq!(body.particles[6].position, pivot, "the pivot must not move");
    for (i, first) in body.particles.iter().enumerate() {
        for (j, second) in body.particles.iter().enumerate().take(i) {
            let rest_distance = (rest[i] - rest[j]).norm();
            let distance = (first.position - second.position).norm();
            assert!(
                (distance - rest_distance).abs() < 1e-3,
                "the patch must stay rigid: distance {i}-{j} is {distance}, rest {rest_distance}"
            );
        }
    }
}

#[test]
fn degenerate_groups_are_left_alone_without_panicking() {
    // Empty group, single particle, coincident particles.
    let particles = [
        SoftParticle::new(Vec3::new(1.0, 2.0, 3.0), 1.0),
        SoftParticle::new(Vec3::new(1.0, 2.0, 3.0), 2.0),
        SoftParticle::new(Vec3::new(-1.0, 0.0, 4.0), 1.0),
    ];
    let before = positions(&particles);
    for indices in [vec![], vec![2], vec![0, 1]] {
        let mut constraint = ShapeMatchingConstraint::from_particles(indices, &particles, 1.0);
        let mut current = particles.clone();
        constraint.project(&mut current, DT);
        assert!(max_displacement(&current, &before) < 1e-12);
    }
}

#[test]
fn invalid_input_is_ignored() {
    let cluster = solid();
    let mut constraint = cluster.constraint(1.0);

    // An index outside the particle array.
    let mut too_short = cluster.rest_particles();
    too_short.truncate(3);
    let before = positions(&too_short);
    constraint.project(&mut too_short, DT);
    assert!(max_displacement(&too_short, &before) < 1e-15);

    // A non-finite particle must not poison the rest of the group.
    let mut broken = cluster.mapped(&rotation([0.0, 1.0, 0.0], 0.4), Vec3::zeros());
    broken[2].position.y = Real::NAN;
    let snapshot = positions(&broken);
    constraint.project(&mut broken, DT);
    for (k, (p, s)) in broken.iter().zip(&snapshot).enumerate() {
        if k != 2 {
            assert_eq!(p.position, *s, "particle {k} was disturbed");
        }
    }
}
