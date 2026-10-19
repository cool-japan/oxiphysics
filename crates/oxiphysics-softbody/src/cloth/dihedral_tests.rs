// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the dihedral bending constraint of [`PbdClothMesh`].
//!
//! The old bend arm of `PbdClothMesh::step` pushed the two wing particles along
//! the line from the hinge midpoint to the wing itself.  That direction lies in
//! the wing's own half plane, so the dihedral angle never changed (the wing
//! only slid towards the hinge), and the angle helper it relied on used
//! unnormalised face normals and an unsigned `acos`.  The tests below check the
//! replacement against independent constructions:
//!
//! * hinges folded by a *known* angle (built from the rotation of one half
//!   plane about the edge), also with oblique wings,
//! * central finite differences of the angle for the analytic gradient,
//! * the conservation laws of an internal constraint (linear and angular
//!   momentum), and
//! * whole-cloth behaviour: a creased patch unfolds, a cantilever droops less.

use super::dihedral::{
    bending_compliance, hinge_state, interior_edge_hinges, project_dihedral_xpbd, wrap_angle,
};
use super::functions::{add, compute_dihedral_angle, cross, dot, length, scale, sub};
use super::types::{ClothConstraint, PbdClothMesh};
use std::f64::consts::{FRAC_PI_2, PI};

/// Deterministic 64-bit linear congruential generator for test geometry.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn next_in(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    fn point(&mut self) -> [f64; 3] {
        [
            self.next_in(-1.0, 1.0),
            self.next_in(-1.0, 1.0),
            self.next_in(-1.0, 1.0),
        ]
    }
}

/// Hinge `[p0, p1, pk, pl]` with the edge along `+x`, folded by the signed
/// dihedral angle `fold`.
///
/// Face A lies in the plane `z = 0` on the `+y` side of the edge; face B is the
/// `-y` half plane rotated about the edge by `fold`.  `along_*` is the position
/// of a wing's foot point along the edge and `height_*` its distance from the
/// edge, so oblique wings (`along != 0` or `along != edge_len`) are allowed.
fn folded_hinge(
    fold: f64,
    edge_len: f64,
    (along_k, height_k): (f64, f64),
    (along_l, height_l): (f64, f64),
) -> [[f64; 3]; 4] {
    [
        [0.0, 0.0, 0.0],
        [edge_len, 0.0, 0.0],
        [along_k, height_k, 0.0],
        [along_l, -height_l * fold.cos(), -height_l * fold.sin()],
    ]
}

/// Rotation of `v` about the unit vector `axis` by `angle` (Rodrigues).
fn rotate(v: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let (sin, cos) = angle.sin_cos();
    add(
        add(scale(v, cos), scale(cross(axis, v), sin)),
        scale(axis, dot(axis, v) * (1.0 - cos)),
    )
}

fn angle_of(hinge: &[[f64; 3]; 4]) -> f64 {
    hinge_state(hinge[0], hinge[1], hinge[2], hinge[3])
        .expect("a regular hinge has an angle")
        .angle
}

/// The helper `compute_dihedral_angle` had before it was made signed: the
/// angle between `e x k_hat` and `e x l_hat`, which are not unit vectors unless
/// the wings are perpendicular to the edge, followed by an unsigned `acos`.
fn legacy_dihedral_angle(p0: [f64; 3], p1: [f64; 3], pk: [f64; 3], pl: [f64; 3]) -> f64 {
    let normalize = |a: [f64; 3]| scale(a, 1.0 / length(a));
    let e_hat = normalize(sub(p1, p0));
    let nk = cross(e_hat, normalize(sub(pk, p0)));
    let nl = cross(e_hat, normalize(sub(pl, p0)));
    dot(nk, nl).clamp(-1.0, 1.0).acos()
}

#[test]
fn signed_angle_matches_the_constructed_fold_for_perpendicular_and_oblique_wings() {
    let wing_layouts = [
        ((0.5, 1.0), (0.5, 1.0)),
        ((0.0, 0.7), (1.0, 1.3)),
        ((0.1, 0.7), (0.9, 1.3)),
        ((-0.4, 2.0), (1.7, 0.3)),
        ((2.5, 0.05), (-1.5, 4.0)),
    ];
    for degrees in [-170.0, -90.0, -30.0, -1.0, 0.0, 1.0, 30.0, 90.0, 170.0] {
        let fold = f64::to_radians(degrees);
        for &(wing_k, wing_l) in &wing_layouts {
            let hinge = folded_hinge(fold, 1.3, wing_k, wing_l);
            let angle = angle_of(&hinge);
            assert!(
                (angle - fold).abs() < 1e-12,
                "fold {degrees} deg, wings {wing_k:?} {wing_l:?}: got {angle}"
            );
            let helper = compute_dihedral_angle(hinge[0], hinge[1], hinge[2], hinge[3]);
            assert_eq!(helper, angle, "the helper must report the same angle");
        }
    }
}

#[test]
fn sign_convention_folds_towards_the_normals_negative() {
    // Face normals of the flat hinge point along +z (see `folded_hinge`).
    let raised = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.1],
        [0.5, -1.0, 0.1],
    ];
    let lowered = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, -0.1],
        [0.5, -1.0, -0.1],
    ];
    let expected = 2.0 * 0.1_f64.atan();
    assert!((angle_of(&raised) + expected).abs() < 1e-12);
    assert!((angle_of(&lowered) - expected).abs() < 1e-12);
}

#[test]
fn flat_hinge_has_zero_angle_whatever_the_wing_positions() {
    // The old helper compared two vectors of length sin(angle(edge, wing)),
    // which only have unit length for perpendicular wings, so it reported a
    // fold for a perfectly flat hinge with oblique wings.
    let flat = folded_hinge(0.0, 1.0, (0.2, 0.8), (0.9, 1.4));
    assert!(angle_of(&flat).abs() < 1e-14);
    let legacy = legacy_dihedral_angle(flat[0], flat[1], flat[2], flat[3]);
    assert!(
        (legacy - PI).abs() > 0.1,
        "the legacy helper was expected to be wrong here, got {legacy}"
    );
    // Perpendicular wings: the legacy value was pi for "flat"; the signed angle
    // is 0 for the same hinge.
    let perpendicular = folded_hinge(0.0, 1.0, (0.0, 1.0), (0.0, 1.0));
    let legacy = legacy_dihedral_angle(
        perpendicular[0],
        perpendicular[1],
        perpendicular[2],
        perpendicular[3],
    );
    assert!((legacy - PI).abs() < 1e-12);
    assert!(angle_of(&perpendicular).abs() < 1e-14);
}

#[test]
fn angle_is_invariant_under_rigid_motions() {
    let mut rng = Lcg(17);
    for _ in 0..20 {
        let fold = rng.next_in(-3.0, 3.0);
        let hinge = folded_hinge(
            fold,
            rng.next_in(0.5, 2.0),
            (rng.next_in(-1.0, 3.0), rng.next_in(0.2, 2.0)),
            (rng.next_in(-1.0, 3.0), rng.next_in(0.2, 2.0)),
        );
        let axis_raw = rng.point();
        let axis = scale(axis_raw, 1.0 / length(axis_raw));
        let turn = rng.next_in(-PI, PI);
        let shift = [
            rng.next_in(-50.0, 50.0),
            rng.next_in(-50.0, 50.0),
            rng.next_in(-50.0, 50.0),
        ];
        let moved = hinge.map(|p| add(rotate(p, axis, turn), shift));
        assert!(
            (angle_of(&moved) - fold).abs() < 1e-9,
            "fold {fold} changed under a rigid motion to {}",
            angle_of(&moved)
        );
    }
}

#[test]
fn degenerate_hinges_have_no_angle() {
    let flat = folded_hinge(0.3, 1.0, (0.5, 1.0), (0.5, 1.0));
    // Zero-length edge.
    assert!(hinge_state(flat[0], flat[0], flat[2], flat[3]).is_none());
    // A wing on the hinge line leaves its face without an orientation.
    assert!(hinge_state(flat[0], flat[1], [0.4, 0.0, 0.0], flat[3]).is_none());
    assert!(hinge_state(flat[0], flat[1], flat[2], [2.0, 0.0, 0.0]).is_none());
    // A wing that coincides with a hinge end.
    assert!(hinge_state(flat[0], flat[1], flat[0], flat[3]).is_none());
    // Non-finite coordinates.
    assert!(hinge_state(flat[0], flat[1], [f64::NAN, 1.0, 0.0], flat[3]).is_none());
    assert!(hinge_state(flat[0], [f64::INFINITY, 0.0, 0.0], flat[2], flat[3]).is_none());
    // The helper falls back to a flat angle instead of returning garbage.
    assert_eq!(
        compute_dihedral_angle(flat[0], flat[0], flat[2], flat[3]),
        0.0
    );
    // Nothing to project.
    let positions = [flat[0], flat[0], flat[2], flat[3]];
    assert!(project_dihedral_xpbd(positions, [1.0; 4], 0.0, 0.0, 0.0).is_none());
}

#[test]
fn analytic_gradient_matches_central_finite_differences() {
    let mut rng = Lcg(2024);
    let step = 1e-6;
    let mut worst = 0.0_f64;
    for _ in 0..60 {
        let hinge = [rng.point(), rng.point(), rng.point(), rng.point()];
        let state = hinge_state(hinge[0], hinge[1], hinge[2], hinge[3])
            .expect("random points form a regular hinge");
        for particle in 0..4 {
            for axis in 0..3 {
                let mut plus = hinge;
                let mut minus = hinge;
                plus[particle][axis] += step;
                minus[particle][axis] -= step;
                let difference = wrap_angle(angle_of(&plus) - angle_of(&minus));
                let numeric = difference / (2.0 * step);
                let analytic = state.gradients[particle][axis];
                let error = (numeric - analytic).abs() / (1.0 + analytic.abs());
                worst = worst.max(error);
                assert!(
                    error < 1e-6,
                    "d angle / d x[{particle}][{axis}]: analytic {analytic}, numeric {numeric}"
                );
            }
        }
    }
    assert!(worst < 1e-6, "worst relative gradient error {worst}");
}

#[test]
fn gradients_are_translation_and_rotation_invariant() {
    let mut rng = Lcg(99);
    for _ in 0..40 {
        let hinge = [rng.point(), rng.point(), rng.point(), rng.point()];
        let state = hinge_state(hinge[0], hinge[1], hinge[2], hinge[3]).expect("regular hinge");
        let magnitude: f64 = state.gradients.iter().map(|g| length(*g)).sum();
        // Translation invariance: the gradients sum to zero.
        let net = state.gradients.iter().fold([0.0; 3], |acc, g| add(acc, *g));
        assert!(length(net) < 1e-12 * (1.0 + magnitude), "net {net:?}");
        // Rotation invariance: no net torque about the origin.
        let torque = hinge
            .iter()
            .zip(&state.gradients)
            .fold([0.0; 3], |acc, (x, g)| add(acc, cross(*x, *g)));
        assert!(
            length(torque) < 1e-11 * (1.0 + magnitude),
            "torque {torque:?}"
        );
        // Sliding a hinge end along the edge does not change the angle.
        let edge_hat = {
            let e = sub(hinge[1], hinge[0]);
            scale(e, 1.0 / length(e))
        };
        assert!(dot(state.gradients[0], edge_hat).abs() < 1e-12 * (1.0 + magnitude));
        assert!(dot(state.gradients[1], edge_hat).abs() < 1e-12 * (1.0 + magnitude));
    }
}

/// Mass-weighted displacement `sum m_i dx_i` and angular momentum change
/// `sum x_i x (m_i dx_i)` of one projection.
fn momentum_changes(
    positions: &[[f64; 3]; 4],
    inv_masses: &[f64; 4],
    corrections: &[[f64; 3]; 4],
) -> ([f64; 3], [f64; 3]) {
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for ((x, w), dx) in positions.iter().zip(inv_masses).zip(corrections) {
        let momentum = scale(*dx, 1.0 / w);
        linear = add(linear, momentum);
        angular = add(angular, cross(*x, momentum));
    }
    (linear, angular)
}

#[test]
fn projection_restores_a_hinge_bent_by_thirty_degrees() {
    let inv_masses = [1.0, 0.5, 0.25, 2.0];
    let mut hinge = folded_hinge(30.0_f64.to_radians(), 1.0, (0.3, 0.8), (0.7, 1.1));
    assert!((angle_of(&hinge) - 30.0_f64.to_radians()).abs() < 1e-12);
    let mut previous = angle_of(&hinge).abs();
    for iteration in 0..12 {
        let projection = project_dihedral_xpbd(hinge, inv_masses, 0.0, 0.0, 0.0)
            .expect("a bent hinge can be projected");
        let (linear, angular) = momentum_changes(&hinge, &inv_masses, &projection.corrections);
        assert!(length(linear) < 1e-12, "iteration {iteration}: {linear:?}");
        assert!(
            length(angular) < 1e-12,
            "iteration {iteration}: {angular:?}"
        );
        for (position, correction) in hinge.iter_mut().zip(projection.corrections) {
            *position = add(*position, correction);
        }
        let now = angle_of(&hinge).abs();
        assert!(
            now < previous || now < 1e-12,
            "iteration {iteration}: the violation must shrink ({previous} -> {now})"
        );
        previous = now;
    }
    assert!(previous < 1e-10, "residual angle {previous}");
}

#[test]
fn projection_drives_the_angle_to_a_nonzero_rest_angle_with_either_sign() {
    for (start_deg, rest_deg) in [(30.0, -20.0), (-45.0, 25.0), (10.0, 80.0), (-10.0, -80.0)] {
        let mut hinge = folded_hinge(f64::to_radians(start_deg), 1.0, (0.2, 0.9), (0.6, 1.2));
        let rest = f64::to_radians(rest_deg);
        for _ in 0..30 {
            if let Some(projection) = project_dihedral_xpbd(hinge, [1.0; 4], rest, 0.0, 0.0) {
                for (position, correction) in hinge.iter_mut().zip(projection.corrections) {
                    *position = add(*position, correction);
                }
            }
        }
        assert!(
            (angle_of(&hinge) - rest).abs() < 1e-9,
            "start {start_deg} deg, rest {rest_deg} deg, got {} deg",
            angle_of(&hinge).to_degrees()
        );
    }
}

#[test]
fn projection_takes_the_short_way_round_the_branch_cut() {
    // Angle -170 deg, rest +170 deg: the constraint violation is +20 deg
    // measured across the -pi / pi cut, not -340 deg through zero.
    let mut hinge = folded_hinge(f64::to_radians(-170.0), 1.0, (0.5, 1.0), (0.5, 1.0));
    let rest = f64::to_radians(170.0);
    let projection = project_dihedral_xpbd(hinge, [1.0; 4], rest, 0.0, 0.0).expect("regular hinge");
    for (position, correction) in hinge.iter_mut().zip(projection.corrections) {
        *position = add(*position, correction);
    }
    // After the first step the hinge is closer to the rest angle along the cut.
    let violation = wrap_angle(angle_of(&hinge) - rest);
    assert!(
        violation.abs() < f64::to_radians(20.0),
        "violation after one step: {} deg",
        violation.to_degrees()
    );
    for _ in 0..30 {
        if let Some(projection) = project_dihedral_xpbd(hinge, [1.0; 4], rest, 0.0, 0.0) {
            for (position, correction) in hinge.iter_mut().zip(projection.corrections) {
                *position = add(*position, correction);
            }
        }
    }
    assert!(wrap_angle(angle_of(&hinge) - rest).abs() < 1e-9);
}

#[test]
fn non_finite_inputs_produce_no_correction() {
    let hinge = folded_hinge(0.4, 1.0, (0.3, 0.8), (0.7, 1.1));
    let bad_inputs = [
        (f64::NAN, 0.0, 0.0),
        (f64::INFINITY, 0.0, 0.0),
        (0.0, f64::NAN, 0.0),
        (0.0, f64::INFINITY, 0.0),
        (0.0, 0.0, f64::NAN),
        (0.0, 0.0, f64::NEG_INFINITY),
    ];
    for (rest, alpha_tilde, lambda) in bad_inputs {
        assert!(
            project_dihedral_xpbd(hinge, [1.0; 4], rest, alpha_tilde, lambda).is_none(),
            "rest {rest}, alpha_tilde {alpha_tilde}, lambda {lambda}"
        );
    }
    let mut broken = hinge;
    broken[3][1] = f64::NAN;
    assert!(project_dihedral_xpbd(broken, [1.0; 4], 0.0, 0.0, 0.0).is_none());
}

#[test]
fn immovable_particles_stay_put_and_the_free_ones_do_the_work() {
    let hinge = folded_hinge(0.5, 1.0, (0.3, 0.8), (0.7, 1.1));
    let inv_masses = [0.0, 1.0, 0.0, 1.0];
    let projection =
        project_dihedral_xpbd(hinge, inv_masses, 0.0, 0.0, 0.0).expect("two free particles");
    assert_eq!(projection.corrections[0], [0.0; 3]);
    assert_eq!(projection.corrections[2], [0.0; 3]);
    assert!(length(projection.corrections[1]) > 0.0);
    assert!(length(projection.corrections[3]) > 0.0);
    // All four immovable: nothing to do.
    assert!(project_dihedral_xpbd(hinge, [0.0; 4], 0.0, 0.0, 0.0).is_none());
}

#[test]
fn compliance_makes_the_constraint_soft_and_satisfies_the_xpbd_fixed_point() {
    let start = folded_hinge(0.6, 1.0, (0.4, 0.9), (0.5, 1.0));
    let initial = angle_of(&start).abs();
    let mut residuals = Vec::new();
    for alpha_tilde in [1e-3, 1.0, 1e3] {
        let mut hinge = start;
        let mut lambda = 0.0;
        for _ in 0..400 {
            let projection = project_dihedral_xpbd(hinge, [1.0; 4], 0.0, alpha_tilde, lambda)
                .expect("regular hinge");
            lambda += projection.delta_lambda;
            for (position, correction) in hinge.iter_mut().zip(projection.corrections) {
                *position = add(*position, correction);
            }
        }
        let violation = angle_of(&hinge);
        // The XPBD stationarity condition C + alpha_tilde * lambda = 0.
        assert!(
            (violation + alpha_tilde * lambda).abs() < 1e-9,
            "alpha_tilde {alpha_tilde}: C = {violation}, lambda = {lambda}"
        );
        residuals.push(violation.abs());
    }
    assert!(residuals[0] < residuals[1] && residuals[1] < residuals[2]);
    assert!(residuals[0] < 1e-3 * initial, "{residuals:?}");
    assert!(
        residuals[2] > 0.5 * initial,
        "a very soft constraint barely moves: {residuals:?}"
    );
}

#[test]
fn compliance_follows_k_dt_squared() {
    assert_eq!(bending_compliance(2.0, 0.1), Some(1.0 / (2.0 * 0.1 * 0.1)));
    let large = bending_compliance(1.0, 1.0 / 60.0).expect("positive stiffness");
    assert!((large - 3600.0).abs() < 1e-9);
    for stiffness in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(bending_compliance(stiffness, 0.01), None, "{stiffness}");
    }
    for dt in [0.0, -0.01, f64::NAN, f64::INFINITY] {
        assert_eq!(bending_compliance(1.0, dt), None, "{dt}");
    }
}

#[test]
fn wrap_angle_maps_to_minus_pi_exclusive_pi_inclusive() {
    for angle in [0.0, 1e-17, -1e-17, 0.3, -0.3, FRAC_PI_2, PI, -FRAC_PI_2] {
        assert_eq!(wrap_angle(angle), angle, "{angle} is already in range");
    }
    assert_eq!(wrap_angle(-PI), PI);
    assert!((wrap_angle(3.0 * PI) - PI).abs() < 1e-12);
    assert!((wrap_angle(-1.5 * PI) - FRAC_PI_2).abs() < 1e-12);
    assert!((wrap_angle(1.5 * PI) + FRAC_PI_2).abs() < 1e-12);
    assert!((wrap_angle(2.0 * PI + 0.25) - 0.25).abs() < 1e-12);
}

#[test]
fn interior_edge_hinges_pair_the_triangles_that_share_an_edge() {
    // Two triangles sharing the edge (1, 2); the second is wound the other way
    // round on that edge.
    let hinges = interior_edge_hinges(&[[0, 1, 2], [1, 3, 2]]);
    assert_eq!(hinges, vec![[1, 2, 0, 3]]);
    // A boundary-only triangle has no hinge.
    assert!(interior_edge_hinges(&[[0, 1, 2]]).is_empty());
    // Three triangles on one edge are non-manifold: no hinge for that edge.
    let fan = interior_edge_hinges(&[[0, 1, 2], [1, 0, 3], [0, 1, 4]]);
    assert!(fan.is_empty(), "{fan:?}");
    // Hinges are ordered by their sorted edge; the edge keeps the direction of
    // the first triangle that holds it (here `3 -> 2` in the second triangle).
    let strip = interior_edge_hinges(&[[0, 1, 2], [1, 3, 2], [2, 3, 4]]);
    assert_eq!(strip, vec![[1, 2, 0, 3], [3, 2, 1, 4]]);
}

/// Largest absolute signed dihedral angle over all interior edges of the
/// triangulation of `cloth` (independent of which constraints the cloth has).
fn max_abs_dihedral(cloth: &PbdClothMesh) -> f64 {
    interior_edge_hinges(&cloth.triangles)
        .iter()
        .map(|&[i, j, k, l]| {
            compute_dihedral_angle(
                cloth.positions[i],
                cloth.positions[j],
                cloth.positions[k],
                cloth.positions[l],
            )
            .abs()
        })
        .fold(0.0, f64::max)
}

fn bend_constraints(cloth: &PbdClothMesh) -> Vec<(usize, usize, usize, usize, f64)> {
    cloth
        .constraints
        .iter()
        .filter_map(|c| match c {
            ClothConstraint::Bend {
                i,
                j,
                k,
                l,
                rest_angle,
                ..
            } => Some((*i, *j, *k, *l, *rest_angle)),
            ClothConstraint::Stretch { .. } => None,
        })
        .collect()
}

fn set_bend_stiffness(cloth: &mut PbdClothMesh, value: f64) {
    for constraint in cloth.constraints.iter_mut() {
        if let ClothConstraint::Bend { stiffness, .. } = constraint {
            *stiffness = value;
        }
    }
}

fn without_bending(cloth: &PbdClothMesh) -> PbdClothMesh {
    let mut copy = cloth.clone();
    copy.constraints
        .retain(|c| matches!(c, ClothConstraint::Stretch { .. }));
    copy
}

#[test]
fn grid_has_one_regular_flat_hinge_per_interior_edge() {
    let (nx, ny) = (5, 4);
    let cloth = PbdClothMesh::grid(nx, ny, 0.25, 1.0);
    let bends = bend_constraints(&cloth);
    // Horizontal, vertical and diagonal interior edges of the triangulation.
    let analytic = (nx - 1) * (ny - 2) + (nx - 2) * (ny - 1) + (nx - 1) * (ny - 1);
    assert_eq!(bends.len(), analytic);
    for (i, j, k, l, rest_angle) in bends {
        assert!(
            rest_angle.abs() < 1e-12,
            "flat grid, rest angle {rest_angle}"
        );
        let mut distinct = [i, j, k, l];
        distinct.sort_unstable();
        assert!(distinct.windows(2).all(|w| w[0] != w[1]), "{distinct:?}");
        // The wings are the apexes of the two triangles that share the edge.
        let mut apexes: Vec<usize> = cloth
            .triangles
            .iter()
            .filter(|tri| tri.contains(&i) && tri.contains(&j))
            .flat_map(|tri| tri.iter().copied().filter(|&v| v != i && v != j))
            .collect();
        apexes.sort_unstable();
        let mut wings = [k, l];
        wings.sort_unstable();
        assert_eq!(apexes, wings, "hinge ({i}, {j}) wings");
        let state = hinge_state(
            cloth.positions[i],
            cloth.positions[j],
            cloth.positions[k],
            cloth.positions[l],
        )
        .expect("the flat grid has no degenerate hinge");
        assert!(state.angle.abs() < 1e-12);
    }
}

/// Fold the patch about the grid row `axis_row`: every particle in a row below
/// it is rotated rigidly about that row's line by `fold` radians.  All edge
/// lengths of the grid are preserved by this motion.
fn crease_about_row(cloth: &mut PbdClothMesh, axis_row: usize, spacing: f64, fold: f64) {
    let axis_z = axis_row as f64 * spacing;
    let nx = cloth.nx;
    for (index, position) in cloth.positions.iter_mut().enumerate() {
        if index / nx < axis_row {
            let lever = position[2] - axis_z;
            *position = [
                position[0],
                -lever * fold.sin(),
                axis_z + lever * fold.cos(),
            ];
        }
    }
}

fn centre_of_mass(cloth: &PbdClothMesh) -> [f64; 3] {
    let total: f64 = cloth.masses.iter().sum();
    cloth
        .positions
        .iter()
        .zip(&cloth.masses)
        .fold([0.0; 3], |acc, (p, m)| add(acc, scale(*p, *m / total)))
}

fn total_momentum(cloth: &PbdClothMesh) -> [f64; 3] {
    cloth
        .velocities
        .iter()
        .zip(&cloth.masses)
        .fold([0.0; 3], |acc, (v, m)| add(acc, scale(*v, *m)))
}

fn creased_patch(fold_degrees: f64) -> PbdClothMesh {
    let spacing = 0.1;
    let mut cloth = PbdClothMesh::grid(7, 7, spacing, 0.02);
    crease_about_row(&mut cloth, 3, spacing, fold_degrees.to_radians());
    cloth
}

#[test]
fn creasing_the_patch_keeps_every_edge_length() {
    let creased = creased_patch(40.0);
    for constraint in &creased.constraints {
        if let ClothConstraint::Stretch { i, j, rest, .. } = constraint {
            let d = length(sub(creased.positions[*j], creased.positions[*i]));
            assert!((d - rest).abs() < 1e-12, "edge ({i}, {j}): {d} vs {rest}");
        }
    }
    assert!((max_abs_dihedral(&creased) - 40.0_f64.to_radians()).abs() < 1e-9);
}

#[test]
fn bending_constraints_unfold_a_creased_patch() {
    let fold = 40.0_f64.to_radians();
    let mut cloth = creased_patch(40.0);
    set_bend_stiffness(&mut cloth, 50.0);
    let mut control = without_bending(&cloth);
    let centre_before = centre_of_mass(&cloth);
    let dt = 1.0 / 60.0;
    for _ in 0..180 {
        cloth.step(dt, [0.0; 3], 10);
        control.step(dt, [0.0; 3], 10);
    }
    // Without bend constraints nothing acts on the creased patch: stretch is
    // satisfied, gravity is off, so the crease stays exactly where it was.
    assert!(
        (max_abs_dihedral(&control) - fold).abs() < 1e-6,
        "control crease {}",
        max_abs_dihedral(&control)
    );
    // The bend constraints flatten it.
    let remaining = max_abs_dihedral(&cloth);
    assert!(
        remaining < 0.02 * fold,
        "crease {} deg of {} deg left after 3 s",
        remaining.to_degrees(),
        fold.to_degrees()
    );
    // Internal forces do not move the centre of mass.
    let centre_after = centre_of_mass(&cloth);
    assert!(
        length(sub(centre_after, centre_before)) < 1e-9,
        "centre of mass drifted from {centre_before:?} to {centre_after:?}"
    );
}

#[test]
fn a_single_bend_step_already_moves_the_crease_towards_flat() {
    // The old arm slid the wings along their own half planes, so the dihedral
    // angle did not change at all within a step.
    for (stiffness, bound) in [(5.0, 0.5), (0.3, 0.75)] {
        let mut cloth = creased_patch(60.0);
        set_bend_stiffness(&mut cloth, stiffness);
        let initial = max_abs_dihedral(&cloth);
        cloth.step(1.0 / 60.0, [0.0; 3], 8);
        let after = max_abs_dihedral(&cloth);
        assert!(
            after < bound * initial,
            "stiffness {stiffness}: {} deg -> {} deg in one step",
            initial.to_degrees(),
            after.to_degrees()
        );
    }
}

#[test]
fn the_default_grid_bend_stiffness_unfolds_a_crease_of_either_sign() {
    for degrees in [40.0, -40.0] {
        let mut cloth = creased_patch(degrees);
        let initial = max_abs_dihedral(&cloth);
        for _ in 0..150 {
            cloth.step(1.0 / 60.0, [0.0; 3], 10);
        }
        let remaining = max_abs_dihedral(&cloth);
        assert!(
            remaining < 0.1 * initial,
            "fold {degrees} deg: {} deg left",
            remaining.to_degrees()
        );
    }
}

#[test]
fn a_bend_step_conserves_linear_momentum_and_the_centre_of_mass() {
    let mut cloth = creased_patch(50.0);
    set_bend_stiffness(&mut cloth, 20.0);
    let centre_before = centre_of_mass(&cloth);
    for _ in 0..20 {
        cloth.step(1.0 / 60.0, [0.0; 3], 6);
        let momentum = total_momentum(&cloth);
        assert!(
            length(momentum) < 1e-9,
            "internal constraints must not create momentum: {momentum:?}"
        );
    }
    assert!(length(sub(centre_of_mass(&cloth), centre_before)) < 1e-9);
}

#[test]
fn pinned_particles_do_not_move_under_bending() {
    let mut cloth = creased_patch(50.0);
    set_bend_stiffness(&mut cloth, 20.0);
    // Pin the crease row and one corner of the flat half.
    for index in 3 * cloth.nx..4 * cloth.nx {
        cloth.pin(index);
    }
    cloth.pin(6 * cloth.nx);
    let pinned_before: Vec<[f64; 3]> = cloth.pinned.iter().map(|&p| cloth.positions[p]).collect();
    let initial = max_abs_dihedral(&cloth);
    for _ in 0..30 {
        cloth.step(1.0 / 60.0, [0.0; 3], 8);
    }
    for (&p, before) in cloth.pinned.iter().zip(&pinned_before) {
        assert_eq!(cloth.positions[p], *before, "pinned particle {p} moved");
        assert_eq!(cloth.velocities[p], [0.0; 3]);
    }
    assert!(
        max_abs_dihedral(&cloth) < 0.5 * initial,
        "the free half must still unfold around the pinned row"
    );
}

#[test]
fn non_positive_bending_stiffness_leaves_the_cloth_untouched_by_bending() {
    let mut zero = creased_patch(35.0);
    set_bend_stiffness(&mut zero, 0.0);
    let mut negative = creased_patch(35.0);
    set_bend_stiffness(&mut negative, -3.0);
    let mut reference = without_bending(&zero);
    for _ in 0..15 {
        zero.step(1.0 / 60.0, [0.0, -9.81, 0.0], 5);
        negative.step(1.0 / 60.0, [0.0, -9.81, 0.0], 5);
        reference.step(1.0 / 60.0, [0.0, -9.81, 0.0], 5);
    }
    assert_eq!(zero.positions, reference.positions);
    assert_eq!(negative.positions, reference.positions);
}

/// Mean height of the free edge of a clamped cantilever: an `8 x 8` cloth lying
/// in the XZ plane, loaded by gravity along `-y`.
///
/// Rows 0 and 1 are pinned.  Pinning a single row would only hinge the sheet at
/// that row (a bend constraint needs a wing on each side of its hinge edge, and
/// the pinned boundary row has none), so the sheet would swing down like a
/// curtain whatever its bending stiffness; the second row fixes the slope at the
/// root.
fn cantilever_free_edge_height(bend_stiffness: Option<f64>) -> f64 {
    let n = 8;
    let mut cloth = PbdClothMesh::grid(n, n, 0.1, 0.05);
    match bend_stiffness {
        Some(value) => set_bend_stiffness(&mut cloth, value),
        None => cloth = without_bending(&cloth),
    }
    for index in 0..2 * n {
        cloth.pin(index);
    }
    for _ in 0..240 {
        cloth.step(1.0 / 60.0, [0.0, -9.81, 0.0], 10);
    }
    let free_edge = &cloth.positions[(n - 1) * n..];
    free_edge.iter().map(|p| p[1]).sum::<f64>() / n as f64
}

#[test]
fn a_clamped_cantilever_droops_less_the_stiffer_its_bending_constraints() {
    let free = cantilever_free_edge_height(None);
    let soft = cantilever_free_edge_height(Some(2.0));
    let stiff = cantilever_free_edge_height(Some(200.0));
    // Without bending resistance the free part of the sheet hangs down.
    assert!(
        free < -0.5,
        "an unsupported sheet droops, free edge at {free}"
    );
    // Every increase of the bending stiffness keeps the free edge higher.
    assert!(
        soft > free + 0.05,
        "soft hinges must hold the edge up: {soft} (soft) vs {free} (none)"
    );
    assert!(
        stiff > soft + 0.1,
        "stiff hinges must hold it higher still: {stiff} (stiff) vs {soft} (soft)"
    );
    assert!(
        stiff > -0.25,
        "a stiff sheet stays nearly horizontal, free edge at {stiff}"
    );
}
