// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the exact force-closure test and of the Ferrari-Canny Q1 / epsilon
//! grasp metrics, including independent cross-checks (closed forms, brute-force
//! facet enumeration and support-function sampling).

use super::inscribed_ball::{epsilon_ball, q1_ball};
use super::polar_cone::{Budget, ConeEnumerator};
use super::*;

// ── deterministic helpers ────────────────────────────────────────────────────

/// Small xorshift64* generator so that the tests are reproducible.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Uniform sample from `(0, 1)`.
    fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64 + 0.5) / (1_u64 << 53) as f64
    }

    /// Standard normal sample (Box-Muller).
    fn gaussian(&mut self) -> f64 {
        let u = self.uniform();
        let v = self.uniform();
        (-2.0 * u.ln()).sqrt() * (2.0 * std::f64::consts::PI * v).cos()
    }

    /// Uniformly distributed unit vector of `R^6`.
    fn unit6(&mut self) -> [f64; 6] {
        loop {
            let mut v = [0.0; 6];
            for c in v.iter_mut() {
                *c = self.gaussian();
            }
            let norm = vec6_len(&v);
            if norm > 1e-3 {
                return [
                    v[0] / norm,
                    v[1] / norm,
                    v[2] / norm,
                    v[3] / norm,
                    v[4] / norm,
                    v[5] / norm,
                ];
            }
        }
    }
}

/// Support function `h(u) = max_i u . w_i` for a unit direction.
fn support(points: &[Wrench], u: &[f64; 6]) -> f64 {
    points
        .iter()
        .map(|w| dot6(w, u))
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Every wrench scaled to unit length.
fn unit_wrenches(wrenches: &[Wrench]) -> Vec<Wrench> {
    wrenches
        .iter()
        .map(|w| {
            let l = vec6_len(w);
            [w[0] / l, w[1] / l, w[2] / l, w[3] / l, w[4] / l, w[5] / l]
        })
        .collect()
}

/// The wrenches `+-scale[k] * e_k` (a cross-polytope).
fn axis_wrenches(scale: [f64; 6]) -> Vec<Wrench> {
    let mut out = Vec::new();
    for (k, s) in scale.iter().enumerate() {
        for sign in [1.0, -1.0] {
            let mut w = [0.0; 6];
            w[k] = sign * s;
            out.push(w);
        }
    }
    out
}

fn contact(position: [f64; 3], normal: [f64; 3], mu: f64) -> ContactPoint {
    ContactPoint {
        position,
        normal,
        mu,
    }
}

/// Six contacts on the coordinate axes at distance one, normals pointing at
/// the origin (a sphere grasp).
fn sphere_grasp(mu: f64) -> Vec<ContactPoint> {
    vec![
        contact([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], mu),
        contact([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], mu),
        contact([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], mu),
        contact([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], mu),
        contact([0.0, 0.0, 1.0], [0.0, 0.0, -1.0], mu),
        contact([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], mu),
    ]
}

/// All index sets `i_0 < .. < i_5` of `{0, .., n-1}`.
fn six_subsets(n: usize) -> Vec<[usize; 6]> {
    let mut out = Vec::new();
    let mut idx = [0_usize, 1, 2, 3, 4, 5];
    if n < 6 {
        return out;
    }
    loop {
        out.push(idx);
        let mut pos = 6;
        while pos > 0 {
            pos -= 1;
            if idx[pos] < n - 6 + pos {
                idx[pos] += 1;
                for next in pos + 1..6 {
                    idx[next] = idx[next - 1] + 1;
                }
                break;
            }
            if pos == 0 {
                return out;
            }
        }
    }
}

/// Solve `A u = (1, .., 1)` by Gaussian elimination with partial pivoting.
fn solve_ones(mut a: [[f64; 6]; 6]) -> Option<[f64; 6]> {
    let mut b = [1.0_f64; 6];
    for col in 0..6 {
        let mut pivot = col;
        for r in col + 1..6 {
            if a[r][col].abs() > a[pivot][col].abs() {
                pivot = r;
            }
        }
        if a[pivot][col].abs() < 1e-10 {
            return None;
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let pivot_row = a[col];
        let pivot_rhs = b[col];
        for r in col + 1..6 {
            let factor = a[r][col] / pivot_row[col];
            for (target, p) in a[r].iter_mut().zip(pivot_row.iter()).skip(col) {
                *target -= factor * p;
            }
            b[r] -= factor * pivot_rhs;
        }
    }
    let mut u = [0.0; 6];
    for row in (0..6).rev() {
        let mut acc = b[row];
        for c in row + 1..6 {
            acc -= a[row][c] * u[c];
        }
        u[row] = acc / a[row][row];
    }
    Some(u)
}

/// Inradius about the origin by brute-force facet enumeration: every
/// supporting hyperplane through six affinely independent points is a facet
/// hyperplane `u . x = 1`; the answer is the smallest `1 / |u|`.  Only valid
/// when the origin is an interior point of the hull.
fn brute_force_inradius(points: &[Wrench]) -> f64 {
    let mut best = f64::INFINITY;
    for subset in six_subsets(points.len()) {
        let mut matrix = [[0.0; 6]; 6];
        for (row, &index) in matrix.iter_mut().zip(subset.iter()) {
            *row = points[index];
        }
        let Some(u) = solve_ones(matrix) else {
            continue;
        };
        if points.iter().all(|w| dot6(w, &u) <= 1.0 + 1e-9) {
            best = best.min(1.0 / vec6_len(&u));
        }
    }
    best
}

// ── the finding's regression: per-axis sign test is not force closure ────────

/// Eleven wrenches from the audit: the first two coordinates carry three
/// vectors that all lie in an open half-plane, the other four axes are
/// covered by `+-e`.
fn half_plane_wrenches(third: [f64; 2]) -> Vec<Wrench> {
    let mut w: Vec<Wrench> = vec![
        [1.0, 0.1, 0.0, 0.0, 0.0, 0.0],
        [-1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        [third[0], third[1], 0.0, 0.0, 0.0, 0.0],
    ];
    for axis in 2..6 {
        for sign in [1.0, -1.0] {
            let mut e = [0.0; 6];
            e[axis] = sign;
            w.push(e);
        }
    }
    w
}

#[test]
fn per_axis_test_accepts_but_exact_test_rejects_a_half_plane_set() {
    let wrenches = half_plane_wrenches([1.0, -1.0]);
    assert_eq!(wrenches.len(), 11);
    // Every axis shows both signs ...
    assert!(force_closure_test(&wrenches));
    // ... yet the origin lies outside the hull of the first three wrenches, so
    // nothing in the xy-plane can be balanced.
    assert!(!force_closure_exact(&wrenches));
    assert!(!force_closure_test_full(&wrenches));
    assert_eq!(grasp_quality_q1(&wrenches), 0.0);
    assert_eq!(grasp_quality_epsilon(&wrenches), 0.0);
}

#[test]
fn repairing_the_half_plane_set_restores_force_closure() {
    // (-0.2, 1) closes the xy-plane: the three directions positively span it.
    let wrenches = half_plane_wrenches([-0.2, 1.0]);
    assert!(force_closure_test(&wrenches));
    assert!(force_closure_exact(&wrenches));
    assert!(grasp_quality_q1(&wrenches) > 0.0);
    assert!(grasp_quality_epsilon(&wrenches) > 0.0);
}

/// The wrenches `+-e_2 .. +-e_5` that complete a planar set to eleven or more.
fn upper_axis_wrenches() -> Vec<Wrench> {
    let mut out = Vec::new();
    for axis in 2..6 {
        for sign in [1.0, -1.0] {
            let mut e = [0.0; 6];
            e[axis] = sign;
            out.push(e);
        }
    }
    out
}

#[test]
fn origin_on_the_boundary_is_not_force_closure() {
    // (1, 1) and (-1, -1) cancel, but (1, -1) cannot be balanced by any
    // non-negative combination: the origin lies on the boundary of the hull of
    // the planar part.  Both signs occur on every axis, so the per-axis
    // prefilter accepts the set, and no strictly positive combination
    // vanishes.
    let mut wrenches: Vec<Wrench> = vec![
        [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
        [-1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    ];
    wrenches.extend(upper_axis_wrenches());
    assert_eq!(wrenches.len(), 11);
    assert!(force_closure_test(&wrenches));
    assert!(!force_closure_exact(&wrenches));
    assert_eq!(grasp_quality_q1(&wrenches), 0.0);
    assert_eq!(grasp_quality_epsilon(&wrenches), 0.0);

    // The missing planar direction (-1, 1) puts the origin in the interior.
    wrenches.push([-1.0, 1.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(force_closure_exact(&wrenches));
}

#[test]
fn free_sum_of_a_square_and_a_cross_polytope_has_a_closed_form_quality() {
    // conv(square (+-1, +-1) in the xy-plane  U  +-e_2 .. +-e_5) is the free sum
    // of a square and a 4-D cross-polytope, so its polar body is the product of
    // their polars: the diamond |y_0| + |y_1| <= 1 (farthest vertex at distance
    // 1) and the cube [-1, 1]^4 (farthest vertex at distance 2).  Hence
    // Q1 = 1 / sqrt(1 + 4).
    let mut wrenches: Vec<Wrench> = Vec::new();
    for sx in [1.0, -1.0] {
        for sy in [1.0, -1.0] {
            wrenches.push([sx, sy, 0.0, 0.0, 0.0, 0.0]);
        }
    }
    wrenches.extend(upper_axis_wrenches());
    assert!(force_closure_exact(&wrenches));
    let q1 = grasp_quality_q1(&wrenches);
    assert!((q1 - 1.0 / 5.0_f64.sqrt()).abs() < 1e-9, "q1 = {q1}");
    // Normalised, the square has vertices at (+-1, +-1) / sqrt(2): its polar
    // diamond reaches sqrt(2), so epsilon = 1 / sqrt(2 + 4).
    let eps = grasp_quality_epsilon(&wrenches);
    assert!((eps - 1.0 / 6.0_f64.sqrt()).abs() < 1e-9, "epsilon = {eps}");
}

// ── wrench-set clean-up ──────────────────────────────────────────────────────

#[test]
fn wrench_clean_up_removes_zero_and_repeated_wrenches() {
    use super::wrench_set::{distinct_nonzero, distinct_unit};
    let wrenches = vec![
        [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0; 6],
        [2.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 3.0, 0.0, 0.0, 0.0, 0.0],
        [1e-14, 0.0, 0.0, 0.0, 0.0, 0.0],
        [4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ];
    // Zero and tiny wrenches vanish, exact repeats collapse.
    let raw = distinct_nonzero(&wrenches).expect("finite input");
    assert_eq!(raw.len(), 3);
    // After normalisation the parallel wrenches coincide as well.
    let unit = distinct_unit(&wrenches).expect("finite input");
    assert_eq!(unit.len(), 2);
    for w in &unit {
        assert!((vec6_len(w) - 1.0).abs() < 1e-15, "{w:?}");
    }
    let mut with_nan = wrenches.clone();
    with_nan.push([f64::NAN; 6]);
    assert!(distinct_nonzero(&with_nan).is_none());
    assert!(distinct_unit(&with_nan).is_none());
}

// ── exact force-closure test ─────────────────────────────────────────────────

#[test]
fn cross_polytope_has_force_closure() {
    assert!(force_closure_exact(&axis_wrenches([1.0; 6])));
}

#[test]
fn six_wrenches_can_never_close_a_six_dimensional_space() {
    // Positive spanning of R^6 needs at least seven wrenches.
    let wrenches = &axis_wrenches([1.0; 6])[..6];
    assert!(!force_closure_exact(wrenches));
    assert_eq!(grasp_quality_q1(wrenches), 0.0);
}

#[test]
fn rank_deficient_wrenches_do_not_close() {
    // Wrenches confined to the hyperplane x_4 = x_5 of wrench space: both signs
    // occur on every axis, so the per-axis test passes, yet the wrenches span
    // only five dimensions and cannot resist a wrench along (0,0,0,0,1,-1).
    let mut wrenches = Vec::new();
    for k in 0..4 {
        for sign in [1.0, -1.0] {
            let mut w = [0.0; 6];
            w[k] = sign;
            wrenches.push(w);
        }
    }
    wrenches.push([0.0, 0.0, 0.0, 0.0, 1.0, 1.0]);
    wrenches.push([0.0, 0.0, 0.0, 0.0, -1.0, -1.0]);
    assert!(force_closure_test(&wrenches));
    assert!(!force_closure_exact(&wrenches));
    assert_eq!(grasp_quality_q1(&wrenches), 0.0);
    assert_eq!(grasp_quality_epsilon(&wrenches), 0.0);
    // The missing direction restores full rank and force closure.
    wrenches.push([0.0, 0.0, 0.0, 0.0, 1.0, -1.0]);
    wrenches.push([0.0, 0.0, 0.0, 0.0, -1.0, 1.0]);
    assert!(force_closure_exact(&wrenches));
    assert!(grasp_quality_q1(&wrenches) > 0.0);
}

#[test]
fn non_finite_and_empty_input_is_rejected() {
    assert!(!force_closure_exact(&[]));
    let mut wrenches = axis_wrenches([1.0; 6]);
    wrenches.push([f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(!force_closure_exact(&wrenches));
    assert_eq!(grasp_quality_q1(&wrenches), 0.0);
    assert_eq!(grasp_quality_epsilon(&wrenches), 0.0);
    wrenches.pop();
    wrenches.push([f64::INFINITY, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(!force_closure_exact(&wrenches));
    assert_eq!(grasp_quality_q1(&[]), 0.0);
}

#[test]
fn duplicates_and_zero_wrenches_do_not_matter() {
    let base = axis_wrenches([1.0, 2.0, 1.0, 1.0, 0.5, 1.0]);
    let mut padded = base.clone();
    padded.extend(base.iter().copied());
    padded.push([0.0; 6]);
    padded.push([1e-15, 0.0, 0.0, 0.0, 0.0, 0.0]);
    assert!(force_closure_exact(&padded));
    let q_base = grasp_quality_q1(&base);
    let q_padded = grasp_quality_q1(&padded);
    assert!((q_base - q_padded).abs() < 1e-12, "{q_base} vs {q_padded}");
    assert!((grasp_quality_epsilon(&base) - grasp_quality_epsilon(&padded)).abs() < 1e-12);
}

/// Scale the wrenches by different powers of ten (positive scaling keeps the
/// positive span).
fn scaled_unevenly(wrenches: &[Wrench]) -> Vec<Wrench> {
    wrenches
        .iter()
        .enumerate()
        .map(|(i, w)| {
            let factor = 10.0_f64.powi((i % 5) as i32 - 2);
            [
                w[0] * factor,
                w[1] * factor,
                w[2] * factor,
                w[3] * factor,
                w[4] * factor,
                w[5] * factor,
            ]
        })
        .collect()
}

#[test]
fn closure_is_independent_of_wrench_magnitudes() {
    let closed = half_plane_wrenches([-0.2, 1.0]);
    assert!(force_closure_exact(&closed));
    assert!(force_closure_exact(&scaled_unevenly(&closed)));
    let open = half_plane_wrenches([1.0, -1.0]);
    assert!(!force_closure_exact(&open));
    assert!(!force_closure_exact(&scaled_unevenly(&open)));
}

#[test]
fn sphere_grasp_with_friction_has_force_closure() {
    let wrenches = multi_contact_wrench_set(&sphere_grasp(0.8), 8);
    assert_eq!(wrenches.len(), 48);
    assert!(force_closure_exact(&wrenches));
    assert!(force_closure_test_full(&wrenches));
    assert!(ContactWrenchCone::new(&sphere_grasp(0.8), 8).is_force_closure());
}

#[test]
fn frictionless_sphere_grasp_cannot_resist_torques() {
    // Six normal forces through the centre produce no torque at all.
    let wrenches = multi_contact_wrench_set(&sphere_grasp(0.0), 8);
    assert!(!force_closure_exact(&wrenches));
    assert_eq!(grasp_quality_q1(&wrenches), 0.0);
}

#[test]
fn two_antipodal_contacts_never_close() {
    for mu in [0.0, 0.5, 2.0] {
        let contacts = [
            contact([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], mu),
            contact([-1.0, 0.0, 0.0], [1.0, 0.0, 0.0], mu),
        ];
        let wrenches = multi_contact_wrench_set(&contacts, 8);
        assert!(!force_closure_exact(&wrenches), "mu = {mu}");
        assert_eq!(grasp_quality_q1(&wrenches), 0.0, "mu = {mu}");
        assert_eq!(grasp_quality_epsilon(&wrenches), 0.0, "mu = {mu}");
        assert!(!ContactWrenchCone::new(&contacts, 8).is_force_closure());
    }
}

// ── Q1 and epsilon: closed forms ─────────────────────────────────────────────

#[test]
fn cross_polytope_inradius_is_one_over_sqrt_six() {
    let wrenches = axis_wrenches([1.0; 6]);
    let expected = 1.0 / 6.0_f64.sqrt();
    let q1 = grasp_quality_q1(&wrenches);
    let eps = grasp_quality_epsilon(&wrenches);
    assert!((q1 - expected).abs() < 1e-9, "q1 = {q1}");
    assert!((eps - expected).abs() < 1e-9, "epsilon = {eps}");
}

#[test]
fn scaled_cross_polytope_matches_the_closed_form() {
    // conv(+-a_k e_k) has the facets sum(+-x_k / a_k) = 1, at distance
    // 1 / sqrt(sum 1 / a_k^2) from the origin.
    let cases = [
        [1.0, 2.0, 3.0, 0.5, 1.5, 1.0],
        [4.0, 4.0, 4.0, 4.0, 4.0, 4.0],
        [0.1, 1.0, 1.0, 1.0, 1.0, 10.0],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1e-3],
    ];
    for scale in cases {
        let expected = 1.0 / scale.iter().map(|a| 1.0 / (a * a)).sum::<f64>().sqrt();
        let q1 = grasp_quality_q1(&axis_wrenches(scale));
        assert!(
            (q1 - expected).abs() <= 1e-9 * expected.max(1.0),
            "scale {scale:?}: q1 = {q1}, expected {expected}"
        );
        // Epsilon normalises the wrench lengths away.
        let eps = grasp_quality_epsilon(&axis_wrenches(scale));
        assert!((eps - 1.0 / 6.0_f64.sqrt()).abs() < 1e-9, "epsilon = {eps}");
    }
}

/// Vertices of a regular simplex centred at the origin, with circumradius one.
fn regular_simplex_6d() -> Vec<Wrench> {
    // Project e_i - mean(e) of R^7 onto an orthonormal (Helmert) basis of the
    // hyperplane sum(x) = 0.
    let helmert = |k: usize, v: &[f64; 7]| -> f64 {
        let kf = k as f64;
        let mut acc = 0.0;
        for item in v.iter().take(k) {
            acc += item;
        }
        acc -= kf * v[k];
        acc / (kf * (kf + 1.0)).sqrt()
    };
    let mut out = Vec::new();
    for i in 0..7 {
        let mut v = [-1.0 / 7.0; 7];
        v[i] += 1.0;
        let mut w = [0.0; 6];
        for (k, c) in w.iter_mut().enumerate() {
            *c = helmert(k + 1, &v);
        }
        let norm = vec6_len(&w);
        for c in w.iter_mut() {
            *c /= norm;
        }
        out.push(w);
    }
    out
}

#[test]
fn regular_simplex_inradius_is_one_sixth_of_the_circumradius() {
    let simplex = regular_simplex_6d();
    assert_eq!(simplex.len(), 7);
    assert!(force_closure_exact(&simplex));
    let q1 = grasp_quality_q1(&simplex);
    assert!((q1 - 1.0 / 6.0).abs() < 1e-9, "q1 = {q1}");
    assert!((grasp_quality_epsilon(&simplex) - 1.0 / 6.0).abs() < 1e-9);
}

// ── Q1 and epsilon: behaviour on contact grasps ──────────────────────────────

#[test]
fn epsilon_is_not_constant_and_responds_to_the_geometry() {
    // The old epsilon normalised every wrench and then returned the minimum
    // wrench norm, i.e. identically 1.
    let reference = multi_contact_wrench_set(&sphere_grasp(0.8), 8);
    let eps_ref = grasp_quality_epsilon(&reference);
    let q1_ref = grasp_quality_q1(&reference);
    assert!(eps_ref > 0.05 && eps_ref < 0.9, "epsilon = {eps_ref}");
    assert!(q1_ref > 0.05 && q1_ref < 2.0, "q1 = {q1_ref}");

    // More friction widens the cones and improves the grasp.
    let slippery = multi_contact_wrench_set(&sphere_grasp(0.2), 8);
    let eps_slippery = grasp_quality_epsilon(&slippery);
    assert!(
        eps_slippery > 0.0 && eps_slippery < eps_ref,
        "epsilon(mu=0.2) = {eps_slippery}, epsilon(mu=0.8) = {eps_ref}"
    );

    // Moving one contact off the axis changes both metrics.
    let mut moved_contacts = sphere_grasp(0.8);
    moved_contacts[4] = contact([0.4, 0.2, 1.0], [0.0, 0.0, -1.0], 0.8);
    let moved = multi_contact_wrench_set(&moved_contacts, 8);
    assert!(force_closure_exact(&moved));
    let eps_moved = grasp_quality_epsilon(&moved);
    let q1_moved = grasp_quality_q1(&moved);
    assert!(
        (eps_moved - eps_ref).abs() > 1e-3,
        "{eps_moved} vs {eps_ref}"
    );
    assert!((q1_moved - q1_ref).abs() > 1e-3, "{q1_moved} vs {q1_ref}");
}

#[test]
fn quality_is_positive_exactly_when_force_closure_holds() {
    let mut rng = Rng::new(7);
    let mut closed = 0;
    let mut open = 0;
    for trial in 0..60 {
        let n = 7 + (trial % 8);
        let wrenches: Vec<Wrench> = (0..n).map(|_| rng.unit6()).collect();
        let exact = force_closure_exact(&wrenches);
        let q1 = grasp_quality_q1(&wrenches);
        assert_eq!(exact, q1 > 1e-9, "trial {trial}: exact {exact}, q1 {q1}");
        if exact {
            closed += 1;
        } else {
            open += 1;
        }
    }
    assert!(
        closed >= 8 && open >= 8,
        "random sets must exercise both outcomes: {closed} closed, {open} open"
    );
}

#[test]
fn q1_matches_brute_force_facet_enumeration() {
    let mut rng = Rng::new(2024);
    let mut compared = 0;
    for trial in 0..40 {
        let n = 10 + (trial % 5);
        let wrenches: Vec<Wrench> = (0..n)
            .map(|_| {
                let u = rng.unit6();
                let scale = 0.5 + rng.uniform();
                [
                    u[0] * scale,
                    u[1] * scale,
                    u[2] * scale,
                    u[3] * scale,
                    u[4] * scale,
                    u[5] * scale,
                ]
            })
            .collect();
        let q1 = grasp_quality_q1(&wrenches);
        if q1 <= 0.0 {
            assert!(!force_closure_exact(&wrenches));
            continue;
        }
        let reference = brute_force_inradius(&wrenches);
        assert!(
            (q1 - reference).abs() < 1e-8,
            "trial {trial}: q1 = {q1}, brute force = {reference}"
        );
        compared += 1;
    }
    assert!(compared >= 15, "only {compared} closed sets compared");
}

#[test]
fn q1_matches_brute_force_on_a_degenerate_contact_grasp() {
    // Six contacts with square (4-sided) friction pyramids: 24 wrenches with
    // many coplanar subsets.
    let contacts = sphere_grasp(0.7);
    let wrenches = multi_contact_wrench_set(&contacts, 4);
    assert_eq!(wrenches.len(), 24);
    let q1 = grasp_quality_q1(&wrenches);
    let reference = brute_force_inradius(&wrenches);
    assert!(q1 > 0.0);
    assert!((q1 - reference).abs() < 1e-8, "{q1} vs {reference}");
    let unit = unit_wrenches(&wrenches);
    let eps = grasp_quality_epsilon(&wrenches);
    let eps_reference = brute_force_inradius(&unit);
    assert!(
        (eps - eps_reference).abs() < 1e-8,
        "{eps} vs {eps_reference}"
    );
}

#[test]
fn support_function_never_drops_below_q1_and_touches_it() {
    let mut sets: Vec<Vec<Wrench>> = vec![
        axis_wrenches([1.0, 2.0, 3.0, 0.5, 1.5, 1.0]),
        regular_simplex_6d(),
        multi_contact_wrench_set(&sphere_grasp(0.8), 8),
        multi_contact_wrench_set(&sphere_grasp(0.3), 6),
    ];
    let mut rng = Rng::new(99);
    for _ in 0..4 {
        sets.push((0..40).map(|_| rng.unit6()).collect());
    }
    for (index, set) in sets.iter().enumerate() {
        let ball = q1_ball(set, Budget::STANDARD);
        assert!(ball.exact);
        assert!(ball.radius > 0.0, "set {index} has no closure");
        // The minimiser u* attains the radius ...
        let at_minimiser = support(set, &ball.direction);
        assert!(
            (at_minimiser - ball.radius).abs() < 1e-9,
            "set {index}: h(u*) = {at_minimiser}, Q1 = {}",
            ball.radius
        );
        // ... and no direction does better.
        let mut best_sample = f64::INFINITY;
        let mut best_direction = [0.0; 6];
        for _ in 0..20_000 {
            let u = rng.unit6();
            let h = support(set, &u);
            assert!(
                h >= ball.radius - 1e-9,
                "set {index}: h(u) = {h} < Q1 = {}",
                ball.radius
            );
            if h < best_sample {
                best_sample = h;
                best_direction = u;
            }
        }
        // Randomised local descent from the best sample cannot undercut Q1.
        let mut step = 0.2;
        let mut current = best_direction;
        let mut current_value = best_sample;
        for _ in 0..4_000 {
            let mut trial = current;
            for c in trial.iter_mut() {
                *c += step * rng.gaussian();
            }
            let norm = vec6_len(&trial);
            for c in trial.iter_mut() {
                *c /= norm;
            }
            let value = support(set, &trial);
            if value < current_value {
                current = trial;
                current_value = value;
            } else {
                step = (step * 0.995).max(1e-6);
            }
            assert!(
                current_value >= ball.radius - 1e-9,
                "set {index}: descent reached {current_value} < Q1 = {}",
                ball.radius
            );
        }
    }
}

#[test]
fn detailed_quality_reports_the_weakest_direction() {
    let wrenches = axis_wrenches([1.0; 6]);
    let detail = grasp_quality_q1_detailed(&wrenches);
    assert!(detail.exact);
    assert!((detail.q1 - grasp_quality_q1(&wrenches)).abs() < 1e-15);
    // The nearest facets of the cross-polytope have normals (+-1, .., +-1)/sqrt(6).
    for c in detail.weakest_direction {
        assert!((c.abs() - 1.0 / 6.0_f64.sqrt()).abs() < 1e-9, "{c}");
    }
    assert!((support(&wrenches, &detail.weakest_direction) - detail.q1).abs() < 1e-9);

    let open = grasp_quality_q1_detailed(&half_plane_wrenches([1.0, -1.0]));
    assert_eq!(open.q1, 0.0);
    assert_eq!(open.weakest_direction, [0.0; 6]);
    assert!(open.exact);
}

#[test]
fn epsilon_is_q1_of_the_normalised_wrenches() {
    let wrenches = multi_contact_wrench_set(&sphere_grasp(0.5), 6);
    let unit = unit_wrenches(&wrenches);
    let direct = q1_ball(&unit, Budget::STANDARD).radius;
    let via_epsilon = epsilon_ball(&wrenches, Budget::STANDARD).radius;
    assert!((direct - via_epsilon).abs() < 1e-12);
    assert!((grasp_quality_epsilon(&wrenches) - via_epsilon).abs() < 1e-15);
}

#[test]
fn exhausted_budget_returns_a_lower_bound() {
    let wrenches = multi_contact_wrench_set(&sphere_grasp(0.8), 8);
    let exact = q1_ball(&wrenches, Budget::STANDARD);
    assert!(exact.exact && exact.radius > 0.0);
    // Allow a handful of rays only: the enumeration stops early.
    let tight = Budget {
        max_rays: 40,
        max_work: Budget::STANDARD.max_work,
    };
    let partial = q1_ball(&wrenches, tight);
    assert!(!partial.exact);
    assert!(
        partial.radius <= exact.radius + 1e-12,
        "partial {} exceeds exact {}",
        partial.radius,
        exact.radius
    );
    let starved = Budget {
        max_rays: Budget::STANDARD.max_rays,
        max_work: 2_000,
    };
    let starved_ball = q1_ball(&wrenches, starved);
    assert!(!starved_ball.exact);
    assert!(starved_ball.radius <= exact.radius + 1e-12);
}

// ── double-description enumerator ────────────────────────────────────────────

/// All extreme rays of `{ x : a . x <= 0 }`: the enumerator is started from
/// the best-conditioned independent rows and every remaining row is inserted.
fn all_extreme_rays(rows: &[Vec<f64>], dim: usize) -> Option<Vec<Vec<f64>>> {
    let mut cone = ConeEnumerator::new(rows, dim, &[], Budget::STANDARD)?;
    for index in 0..cone.row_count() {
        cone.insert(index).ok()?;
    }
    Some(cone.rays().map(|ray| ray.to_vec()).collect())
}

/// Rows `(+-e_k, -1)` of the cone over the cube `[-1, 1]^3` (dimension 4).
fn cube_cone_rows() -> Vec<Vec<f64>> {
    let mut rows = Vec::new();
    for k in 0..3 {
        for sign in [1.0, -1.0] {
            let mut row = vec![0.0; 4];
            row[k] = sign;
            row[3] = -1.0;
            rows.push(row);
        }
    }
    rows
}

#[test]
fn cube_cone_has_eight_extreme_rays() {
    let rays = all_extreme_rays(&cube_cone_rows(), 4).expect("pointed cone");
    assert_eq!(rays.len(), 8);
    for ray in &rays {
        // (+-1, +-1, +-1, 1) / 2
        assert!((ray[3] - 0.5).abs() < 1e-12, "{ray:?}");
        for c in &ray[..3] {
            assert!((c.abs() - 0.5).abs() < 1e-12, "{ray:?}");
        }
    }
}

#[test]
fn redundant_and_repeated_rows_do_not_change_the_rays() {
    let mut rows = cube_cone_rows();
    let copy = rows.clone();
    rows.extend(copy);
    // Loose copies of the faces and the redundant half-space s >= 0.
    for k in 0..3 {
        for sign in [1.0, -1.0] {
            let mut row = vec![0.0; 4];
            row[k] = sign;
            row[3] = -3.0;
            rows.push(row);
        }
    }
    rows.push(vec![0.0, 0.0, 0.0, -1.0]);
    rows.push(vec![0.0; 4]);
    let rays = all_extreme_rays(&rows, 4).expect("pointed cone");
    assert_eq!(rays.len(), 8);
}

#[test]
fn degenerate_vertices_are_enumerated_once() {
    // Polar of the octahedron conv(+-e_k) in R^3: the facets sum(+-y_k) <= 1
    // describe a cube whose eight rows meet in the six vertices four at a time
    // (a degenerate arrangement).
    let mut rows = Vec::new();
    for s0 in [1.0, -1.0] {
        for s1 in [1.0, -1.0] {
            for s2 in [1.0, -1.0] {
                rows.push(vec![s0, s1, s2, -1.0]);
            }
        }
    }
    let rays = all_extreme_rays(&rows, 4).expect("pointed cone");
    assert_eq!(rays.len(), 6);
    for ray in &rays {
        let nonzero = ray[..3].iter().filter(|c| c.abs() > 1e-9).count();
        assert_eq!(nonzero, 1, "{ray:?}");
        assert!((ray[3] - (1.0 / 2.0_f64.sqrt())).abs() < 1e-12);
    }
}

#[test]
fn insertion_order_does_not_change_the_rays() {
    let rows = cube_cone_rows();
    let forward = all_extreme_rays(&rows, 4).expect("pointed cone");
    let mut reversed_rows = rows.clone();
    reversed_rows.reverse();
    let reversed = all_extreme_rays(&reversed_rows, 4).expect("pointed cone");
    assert_eq!(forward.len(), reversed.len());
    for ray in &forward {
        assert!(
            reversed
                .iter()
                .any(|other| ray.iter().zip(other).all(|(a, b)| (a - b).abs() < 1e-9)),
            "ray {ray:?} missing after reordering"
        );
    }
}

#[test]
fn cones_with_a_line_have_no_extreme_rays() {
    // Only constrains the first coordinate of R^3.
    assert!(ConeEnumerator::new(&[vec![1.0, 0.0, 0.0]], 3, &[], Budget::STANDARD).is_none());
    assert!(ConeEnumerator::new(&[], 2, &[], Budget::STANDARD).is_none());
}

#[test]
fn malformed_rows_are_rejected() {
    let budget = Budget::STANDARD;
    assert!(ConeEnumerator::new(&[vec![1.0, 0.0], vec![0.0, 1.0, 0.0]], 2, &[], budget).is_none());
    assert!(ConeEnumerator::new(&[vec![f64::NAN, 0.0], vec![0.0, 1.0]], 2, &[], budget).is_none());
    assert!(ConeEnumerator::new(&[vec![1.0]], 0, &[], budget).is_none());
    // Seed rows must exist and be independent.
    let rows = [vec![-1.0, 0.0], vec![0.0, -1.0], vec![-2.0, 0.0]];
    assert!(ConeEnumerator::new(&rows, 2, &[7], budget).is_none());
    assert!(ConeEnumerator::new(&rows, 2, &[0, 2], budget).is_none());
    assert!(ConeEnumerator::new(&rows, 2, &[0, 1, 2], budget).is_none());
    assert!(ConeEnumerator::new(&rows, 2, &[0, 1], budget).is_some());
}

#[test]
fn planar_cone_rays_are_its_two_boundary_directions() {
    // {x >= 0, y >= 0} written as -x <= 0, -y <= 0.
    let rays = all_extreme_rays(&[vec![-1.0, 0.0], vec![0.0, -1.0]], 2).expect("pointed cone");
    assert_eq!(rays.len(), 2);
    let mut firsts: Vec<f64> = rays.iter().map(|r| r[0]).collect();
    firsts.sort_by(f64::total_cmp);
    assert!(firsts[0].abs() < 1e-12 && (firsts[1] - 1.0).abs() < 1e-12);
}

#[test]
fn cutting_the_planar_cone_keeps_two_rays() {
    // x >= 0, y >= 0, then x <= y.
    let rows = [vec![-1.0, 0.0], vec![0.0, -1.0], vec![1.0, -1.0]];
    let rays = all_extreme_rays(&rows, 2).expect("pointed cone");
    assert_eq!(rays.len(), 2);
    let diagonal = 1.0 / 2.0_f64.sqrt();
    assert!(
        rays.iter()
            .any(|r| (r[0] - diagonal).abs() < 1e-12 && (r[1] - diagonal).abs() < 1e-12)
    );
    assert!(
        rays.iter()
            .any(|r| r[0].abs() < 1e-12 && (r[1] - 1.0).abs() < 1e-12)
    );
}

#[test]
fn refused_insertion_leaves_the_enumerator_unchanged() {
    let starved = Budget {
        max_rays: Budget::STANDARD.max_rays,
        max_work: 10,
    };
    let rows = cube_cone_rows();
    let mut cone = ConeEnumerator::new(&rows, 4, &[], starved).expect("pointed cone");
    // The starting cone is simplicial: four rays in dimension four.
    let before: Vec<Vec<f64>> = cone.rays().map(|ray| ray.to_vec()).collect();
    assert_eq!(before.len(), 4);
    let mut refused = 0;
    for index in 0..cone.row_count() {
        if !cone.is_processed(index) && cone.insert(index).is_err() {
            refused += 1;
            assert!(!cone.is_processed(index));
        }
    }
    assert!(refused > 0, "the tiny budget must refuse an insertion");
    let after: Vec<Vec<f64>> = cone.rays().map(|ray| ray.to_vec()).collect();
    assert_eq!(
        before, after,
        "a refused insertion must not change the rays"
    );
    // With the standard budget the same insertions go through.
    let mut cone = ConeEnumerator::new(&rows, 4, &[], Budget::STANDARD).expect("pointed cone");
    for index in 0..cone.row_count() {
        cone.insert(index).expect("within the standard budget");
        assert!(cone.is_processed(index));
    }
    assert_eq!(cone.rays().count(), 8);
}

/// All `k`-element index subsets of `{0, .., n-1}` in lexicographic order.
fn subsets(n: usize, k: usize) -> Vec<Vec<usize>> {
    fn extend(
        start: usize,
        n: usize,
        k: usize,
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == k {
            out.push(current.clone());
            return;
        }
        for next in start..n {
            current.push(next);
            extend(next + 1, n, k, current, out);
            current.pop();
        }
    }
    let mut out = Vec::new();
    extend(0, n, k, &mut Vec::new(), &mut out);
    out
}

/// Unit vector orthogonal to all `dim - 1` rows, if they are independent.
fn null_vector(rows: &[&Vec<f64>], dim: usize) -> Option<Vec<f64>> {
    let mut m: Vec<Vec<f64>> = rows.iter().map(|r| (*r).clone()).collect();
    let mut pivot_columns = Vec::new();
    let mut row = 0;
    for col in 0..dim {
        if row == m.len() {
            break;
        }
        let mut best = row;
        for r in row + 1..m.len() {
            if m[r][col].abs() > m[best][col].abs() {
                best = r;
            }
        }
        if m[best][col].abs() < 1e-9 {
            continue;
        }
        m.swap(row, best);
        let pivot = m[row][col];
        for value in m[row].iter_mut() {
            *value /= pivot;
        }
        let pivot_row = m[row].clone();
        for (r, other) in m.iter_mut().enumerate() {
            if r != row {
                let factor = other[col];
                for (value, p) in other.iter_mut().zip(pivot_row.iter()) {
                    *value -= factor * p;
                }
            }
        }
        pivot_columns.push(col);
        row += 1;
    }
    if pivot_columns.len() != dim - 1 {
        return None;
    }
    let free = (0..dim).find(|c| !pivot_columns.contains(c))?;
    let mut x = vec![0.0; dim];
    x[free] = 1.0;
    for (i, &p) in pivot_columns.iter().enumerate() {
        x[p] = -m[i][free];
    }
    let norm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    Some(x.iter().map(|v| v / norm).collect())
}

/// Rank of a set of rows (Gaussian elimination with partial pivoting).
fn rank_of(rows: &[Vec<f64>], dim: usize) -> usize {
    let mut m: Vec<Vec<f64>> = rows.to_vec();
    let mut rank = 0;
    for col in 0..dim {
        if rank == m.len() {
            break;
        }
        let mut best = rank;
        for r in rank + 1..m.len() {
            if m[r][col].abs() > m[best][col].abs() {
                best = r;
            }
        }
        if m[best][col].abs() < 1e-9 {
            continue;
        }
        m.swap(rank, best);
        let pivot_row = m[rank].clone();
        for other in m.iter_mut().skip(rank + 1) {
            let factor = other[col] / pivot_row[col];
            for (value, p) in other.iter_mut().zip(pivot_row.iter()) {
                *value -= factor * p;
            }
        }
        rank += 1;
    }
    rank
}

/// Extreme rays of `{ x : a . x <= 0 }` by exhaustive enumeration: a feasible
/// unit vector orthogonal to `dim - 1` independent rows.
fn brute_force_rays(rows: &[Vec<f64>], dim: usize) -> Vec<Vec<f64>> {
    let unit: Vec<Vec<f64>> = rows
        .iter()
        .filter_map(|r| {
            let norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
            (norm > 0.0).then(|| r.iter().map(|v| v / norm).collect())
        })
        .collect();
    let mut found: Vec<Vec<f64>> = Vec::new();
    for subset in subsets(unit.len(), dim - 1) {
        let chosen: Vec<&Vec<f64>> = subset.iter().map(|&i| &unit[i]).collect();
        let Some(x) = null_vector(&chosen, dim) else {
            continue;
        };
        for sign in [1.0, -1.0] {
            let candidate: Vec<f64> = x.iter().map(|v| sign * v).collect();
            let feasible = unit.iter().all(|a| {
                a.iter()
                    .zip(candidate.iter())
                    .map(|(p, q)| p * q)
                    .sum::<f64>()
                    <= 1e-9
            });
            let known = found.iter().any(|r| {
                r.iter()
                    .zip(candidate.iter())
                    .map(|(p, q)| p * q)
                    .sum::<f64>()
                    > 1.0 - 1e-9
            });
            if feasible && !known {
                found.push(candidate);
            }
        }
    }
    found
}

fn assert_same_rays(actual: &[Vec<f64>], expected: &[Vec<f64>], label: &str) {
    assert_eq!(actual.len(), expected.len(), "{label}: ray count");
    for ray in expected {
        let hit = actual
            .iter()
            .any(|a| a.iter().zip(ray.iter()).map(|(p, q)| p * q).sum::<f64>() > 1.0 - 1e-7);
        assert!(hit, "{label}: missing ray {ray:?}");
    }
}

#[test]
fn enumerator_matches_brute_force_on_random_cones() {
    let mut rng = Rng::new(31);
    let mut nontrivial = 0;
    for trial in 0..120 {
        let dim = 3 + trial % 3;
        let count = dim + 2 + (rng.next_u64() % 6) as usize;
        // Rows with a common strictly negative component along (1, .., 1), so
        // the cone has a non-empty interior.
        let rows: Vec<Vec<f64>> = (0..count)
            .map(|_| {
                let mut row: Vec<f64> = (0..dim).map(|_| rng.gaussian()).collect();
                let mean = row.iter().sum::<f64>() / dim as f64;
                for v in row.iter_mut() {
                    *v -= mean + 0.3;
                }
                row
            })
            .collect();
        if rank_of(&rows, dim) < dim {
            // The cone contains a line and has no extreme rays.
            assert!(all_extreme_rays(&rows, dim).is_none(), "trial {trial}");
            continue;
        }
        let expected = brute_force_rays(&rows, dim);
        let actual = all_extreme_rays(&rows, dim).expect("pointed cone");
        assert_same_rays(&actual, &expected, &format!("trial {trial}"));
        if actual.len() > dim {
            nontrivial += 1;
        }
    }
    assert!(
        nontrivial >= 40,
        "only {nontrivial} cones had more than dim rays"
    );
}

#[test]
fn enumerator_matches_brute_force_on_degenerate_cones() {
    // Rows on the coarse grid {-1, 0, 1}^dim: many rows share rays, several
    // facets meet in one ray and some rows are redundant or repeated.
    let mut rng = Rng::new(77);
    let mut degenerate = 0;
    for trial in 0..80 {
        let dim = 3 + trial % 2;
        let count = dim + 3 + (rng.next_u64() % 6) as usize;
        let mut rows: Vec<Vec<f64>> = (0..count)
            .map(|_| {
                (0..dim)
                    .map(|_| ((rng.next_u64() % 3) as f64) - 1.0)
                    .collect::<Vec<f64>>()
            })
            .collect();
        // Make (1, .., 1) an interior direction of the cone.
        for row in rows.iter_mut() {
            if row.iter().sum::<f64>() >= 0.0 {
                for v in row.iter_mut() {
                    *v = -*v;
                }
            }
        }
        rows.retain(|r| r.iter().sum::<f64>() < 0.0);
        if rows.len() < dim {
            continue;
        }
        if rank_of(&rows, dim) < dim {
            assert!(all_extreme_rays(&rows, dim).is_none(), "trial {trial}");
            continue;
        }
        let expected = brute_force_rays(&rows, dim);
        let actual = all_extreme_rays(&rows, dim).expect("pointed cone");
        assert_same_rays(&actual, &expected, &format!("trial {trial}"));
        // Degenerate: some ray lies on more than dim - 1 of the hyperplanes.
        let on_planes = |ray: &Vec<f64>| {
            rows.iter()
                .filter(|r| {
                    let norm = r.iter().map(|v| v * v).sum::<f64>().sqrt();
                    let along: f64 = r.iter().zip(ray.iter()).map(|(a, b)| a * b).sum();
                    along.abs() / norm < 1e-9
                })
                .count()
        };
        if actual.iter().any(|ray| on_planes(ray) > dim - 1) {
            degenerate += 1;
        }
    }
    assert!(
        degenerate >= 10,
        "only {degenerate} degenerate cones exercised"
    );
}

// ── consumers ────────────────────────────────────────────────────────────────

#[test]
fn contact_stability_report_uses_the_exact_metrics() {
    let contacts = sphere_grasp(0.8);
    let report = analyse_contact_stability(&contacts, [0.0; 3], [0.0, 0.0, -9.81], 1.0, 8);
    assert!(report.force_closure);
    let wrenches = multi_contact_wrench_set(&contacts, 8);
    assert!((report.q1 - grasp_quality_q1(&wrenches)).abs() < 1e-12);
    assert!((report.epsilon - grasp_quality_epsilon(&wrenches)).abs() < 1e-12);
    assert!(report.epsilon > 0.05 && report.epsilon < 0.9);
    // The epsilon term of the score is no longer the constant 0.3.
    let expected_epsilon_term = report.epsilon * 0.3;
    let closure_term = 0.4;
    assert!(report.score >= closure_term + expected_epsilon_term - 1e-12);
    assert!(report.score <= 1.0 + 1e-12);
    let direct = multi_contact_stability_score(&contacts, [0.0; 3], 8);
    assert!((report.score - direct).abs() < 1e-12);
}

#[test]
fn open_grasps_score_without_the_closure_bonus() {
    let contacts = [contact([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.5)];
    let report = analyse_contact_stability(&contacts, [0.0, 0.0, 1.0], [0.0, 0.0, -9.81], 1.0, 8);
    assert!(!report.force_closure);
    assert_eq!(report.q1, 0.0);
    assert_eq!(report.epsilon, 0.0);
    assert!(report.score < 0.4, "score = {}", report.score);
    assert!(multi_contact_stability_score(&contacts, [0.0, 0.0, 1.0], 8) < 0.4);
    assert_eq!(multi_contact_stability_score(&[], [0.0; 3], 8), 0.0);
}

#[test]
fn friction_cone_quality_grows_with_the_friction_coefficient() {
    let mut previous = 0.0;
    for mu in [0.15, 0.3, 0.6, 1.2] {
        let wrenches = multi_contact_wrench_set(&sphere_grasp(mu), 8);
        let q1 = grasp_quality_q1(&wrenches);
        assert!(q1 > previous, "mu = {mu}: q1 = {q1} !> {previous}");
        previous = q1;
    }
}
