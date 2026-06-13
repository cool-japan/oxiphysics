// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the sparse `LᵀDL` mass-matrix factorization.

use oxiphysics_articulated::body::RigidBody;
use oxiphysics_articulated::crba::compute_mass_matrix_crba;
use oxiphysics_articulated::joint::RevoluteJoint;
use oxiphysics_articulated::ltdl::*;
use oxiphysics_articulated::model::ArticulatedModel;
use oxiphysics_articulated::spatial::{SpatialInertia, SpatialTransform};

fn link(mass: f64, com: [f64; 3], parent: Option<usize>) -> RigidBody {
    let inertia = SpatialInertia::from_com(mass, com, [[0.0; 3]; 3]);
    RigidBody::new("link", inertia, parent, SpatialTransform::IDENTITY)
}

/// torso(0) → L1(1) → L2(2) ; torso(0) → R1(3) → R2(4). 5 revolute DOFs, real branch at the torso.
fn build_branched_humanoid() -> ArticulatedModel {
    let mut m = ArticulatedModel::new([0.0, 0.0, -9.81]);
    let torso = m.add_body(
        link(2.0, [0.0, 0.0, 0.2], None),
        Box::new(RevoluteJoint::new([0.0, 0.0, 1.0])),
    );
    let l1 = m.add_body(
        link(1.0, [0.3, 0.0, 0.0], Some(torso)),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );
    let _l2 = m.add_body(
        link(0.8, [0.3, 0.0, 0.0], Some(l1)),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );
    let r1 = m.add_body(
        link(1.0, [-0.3, 0.0, 0.0], Some(torso)),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );
    let _r2 = m.add_body(
        link(0.8, [-0.3, 0.0, 0.0], Some(r1)),
        Box::new(RevoluteJoint::new([0.0, 1.0, 0.0])),
    );
    m
}

/// Build a deep + branched tree of revolute links for the performance test.
///
/// A central spine of `spine_len` links, with `n_limbs` limbs (each
/// `limb_len` links long) hanging off the *top* spine link. Offset COMs and
/// nonzero masses keep the joint-space inertia strictly SPD. The returned
/// model has `spine_len + n_limbs * limb_len` revolute DOFs and a genuine
/// branching structure (so the sparse pattern differs materially from dense).
fn build_perf_tree(spine_len: usize, n_limbs: usize, limb_len: usize) -> ArticulatedModel {
    let mut m = ArticulatedModel::new([0.0, 0.0, -9.81]);
    // Spine: a vertical chain rotating about mixed axes for richness.
    let mut prev: Option<usize> = None;
    let mut spine_top = 0usize;
    for s in 0..spine_len {
        let axis = if s % 2 == 0 {
            [0.0, 1.0, 0.0]
        } else {
            [1.0, 0.0, 0.0]
        };
        let com = [0.05, 0.02, 0.2 + 0.01 * s as f64];
        let id = m.add_body(
            link(1.5 + 0.05 * s as f64, com, prev),
            Box::new(RevoluteJoint::new(axis)),
        );
        prev = Some(id);
        spine_top = id;
    }
    // Limbs branch off the top of the spine.
    for limb in 0..n_limbs {
        let mut parent = Some(spine_top);
        for j in 0..limb_len {
            let sign = if limb % 2 == 0 { 1.0 } else { -1.0 };
            let axis = if j % 2 == 0 {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, 1.0]
            };
            let com = [sign * (0.3 + 0.02 * j as f64), 0.03, 0.01];
            let id = m.add_body(
                link(0.8 + 0.03 * j as f64, com, parent),
                Box::new(RevoluteJoint::new(axis)),
            );
            parent = Some(id);
        }
    }
    m
}

// ---------------------------------------------------------------------------
// Local dense Cholesky (copied from src/aba_derivatives.rs, lines 91-137) so
// the perf test compares against a real dense baseline with no extra deps.
// These are private test helpers; on a non-positive pivot they just `.expect`
// because the performance matrix is known to be SPD.
// ---------------------------------------------------------------------------

/// In-place lower-triangular Cholesky factor `L` of an SPD matrix (`a = L·Lᵀ`).
fn dense_cholesky_factor(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for (k, l_jk) in l[j].iter().enumerate().take(j) {
                sum -= l[i][k] * l_jk;
            }
            if i == j {
                assert!(sum > 0.0, "dense Cholesky: non-positive pivot at {i}");
                l[i][i] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    l
}

/// Solve `M·x = b` from a precomputed dense Cholesky factor `l`.
fn dense_cholesky_solve_vec(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut y = vec![0.0f64; n];
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i];
    }
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut sum = y[i];
        for k in (i + 1)..n {
            sum -= l[k][i] * x[k];
        }
        x[i] = sum / l[i][i];
    }
    x
}

/// Dense matrix-vector product `H·x` using the full symmetric `h`.
fn spmv(h: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
    let n = h.len();
    let mut y = vec![0.0f64; n];
    for (i, yi) in y.iter_mut().enumerate() {
        let mut sum = 0.0;
        for (j, xj) in x.iter().enumerate() {
            sum += h[i][j] * xj;
        }
        *yi = sum;
    }
    y
}

/// Enumerate the (1-based) ancestor set of DOF `k` by walking the `λ` chain.
fn ancestors_of(lambda: &[usize], k: usize) -> Vec<usize> {
    let mut set = Vec::new();
    let mut i = lambda[k];
    while i != 0 {
        set.push(i);
        i = lambda[i];
    }
    set
}

#[test]
fn test_lambda_correctness() {
    let model = build_branched_humanoid();
    let n = model.total_dof();
    let lambda = build_lambda(&model);

    // Invariant: parent DOF index strictly precedes the child for non-roots.
    for (k, &parent) in lambda.iter().enumerate().take(n + 1).skip(1) {
        assert!(parent < k, "lambda[{k}] = {parent} must be < {k}");
    }

    // Torso = body 0; its first (and only) DOF is a root ⇒ λ == 0.
    let torso_first = model.dof_start(0) + 1;
    assert_eq!(lambda[torso_first], 0, "root DOF must have λ == 0");

    // Torso's last 1-based DOF index.
    let torso_last = model.dof_start(0) + model.dof_count(0);

    // L1 = body 1, R1 = body 3: both branch directly off the torso, so the
    // first DOF of each must point at the torso's last DOF.
    let l1_first = model.dof_start(1) + 1;
    let r1_first = model.dof_start(3) + 1;
    assert_eq!(
        lambda[l1_first], torso_last,
        "L1 first DOF must map to torso last DOF"
    );
    assert_eq!(
        lambda[r1_first], torso_last,
        "R1 first DOF must map to torso last DOF"
    );

    // L2 = body 2 chains off L1; R2 = body 4 chains off R1.
    let l2_first = model.dof_start(2) + 1;
    let r2_first = model.dof_start(4) + 1;
    assert_eq!(lambda[l2_first], model.dof_start(1) + model.dof_count(1));
    assert_eq!(lambda[r2_first], model.dof_start(3) + model.dof_count(3));
}

#[test]
fn test_ltdl_factor_solve_vs_dense_cholesky() {
    // HARD GATE: sparse LᵀDL solve must match the original H·x = b to 1e-10.
    let mut model = build_branched_humanoid();
    let n = model.total_dof();
    let q = vec![0.0; n];
    let h = compute_mass_matrix_crba(&mut model, &q);
    let lambda = build_lambda(&model);

    let fact = SparseMassFactorization::factor(h.clone(), lambda, n)
        .expect("branched humanoid mass matrix must be SPD");

    let b: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) * 0.1 + 0.37).collect();
    let x = fact.solve(&b).expect("solve must succeed");

    // Residual check via dense SpMV using the ORIGINAL full-symmetric h.
    for i in 0..n {
        let mut row = 0.0;
        for j in 0..n {
            row += h[i][j] * x[j];
        }
        assert!(
            (row - b[i]).abs() < 1e-10,
            "row {i}: H·x = {row}, b = {}, |diff| = {}",
            b[i],
            (row - b[i]).abs()
        );
    }
}

#[test]
fn test_ltdl_reconstruction() {
    let mut model = build_branched_humanoid();
    let n = model.total_dof();
    let q = vec![0.0; n];
    let h = compute_mass_matrix_crba(&mut model, &q);
    let lambda = build_lambda(&model);
    let fact = SparseMassFactorization::factor(h.clone(), lambda, n).expect("SPD");

    // Solve H⁻¹ e_k for each k, SpMV back and check ≈ δ_ik.
    for k in 0..n {
        let mut e = vec![0.0; n];
        e[k] = 1.0;
        let col = fact.solve(&e).expect("solve");
        let back = spmv(&h, &col);
        for (i, &bi) in back.iter().enumerate() {
            let expected = if i == k { 1.0 } else { 0.0 };
            assert!(
                (bi - expected).abs() < 1e-9,
                "reconstruction[{i}][{k}] = {bi}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_sparsity_no_fill_in() {
    let mut model = build_branched_humanoid();
    let n = model.total_dof();
    let q = vec![0.0; n];
    let h = compute_mass_matrix_crba(&mut model, &q);
    let lambda = build_lambda(&model);
    let fact = SparseMassFactorization::factor(h, lambda.clone(), n).expect("SPD");

    // For each row, any strictly-lower column that is NOT an ancestor must be
    // (and remain) exactly zero — i.e. zero fill-in.
    for k in 1..=n {
        let anc = ancestors_of(&lambda, k);
        let r = k - 1;
        for c in 0..r {
            let one_based = c + 1;
            if !anc.contains(&one_based) {
                let v = fact.factored_entry(r, c);
                assert!(
                    v.abs() < 1e-12,
                    "fill-in at ({r},{c}) = {v} (DOF {one_based} not an ancestor of {k})"
                );
            }
        }
    }
}

#[test]
fn test_factor_crba_convenience() {
    let mut model = build_branched_humanoid();
    let n = model.total_dof();
    let q = vec![0.0; n];
    let fact = factor_crba(&mut model, &q).expect("factor_crba must succeed");
    assert_eq!(fact.n_dof(), n);

    let b = vec![1.0; n];
    let x = fact.solve(&b).expect("solve");
    assert_eq!(x.len(), n);
    assert!(x.iter().all(|v| v.is_finite()), "all solution entries finite");
}

#[test]
fn test_solve_columns_matches_individual() {
    let mut model = build_branched_humanoid();
    let n = model.total_dof();
    let q = vec![0.0; n];
    let fact = factor_crba(&mut model, &q).expect("factor_crba");

    let cols: Vec<Vec<f64>> = vec![
        (0..n).map(|i| i as f64 + 1.0).collect(),
        (0..n).map(|i| (n - i) as f64 * 0.5).collect(),
        vec![1.0; n],
    ];

    let batched = fact.solve_columns(&cols).expect("solve_columns");
    assert_eq!(batched.len(), cols.len());
    for (c, col) in cols.iter().enumerate() {
        let individual = fact.solve(col).expect("solve");
        for (i, (&a, &b)) in batched[c].iter().zip(individual.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-15,
                "column {c} entry {i}: batched {a} vs individual {b}"
            );
        }
    }
}

#[test]
fn test_dimension_mismatch_errors() {
    // Non-square h ⇒ DimensionMismatch.
    let bad_h = vec![vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]];
    let res = SparseMassFactorization::factor(bad_h, vec![0, 0, 0], 2);
    assert!(
        matches!(res, Err(LtdlError::DimensionMismatch { .. })),
        "non-square h must yield DimensionMismatch, got {res:?}"
    );

    // Wrong-length b ⇒ DimensionMismatch on solve.
    let good_h = vec![vec![2.0, 0.0], vec![0.0, 3.0]];
    let fact = SparseMassFactorization::factor(good_h, vec![0, 0, 0], 2).expect("SPD 2x2");
    let res = fact.solve(&[1.0, 2.0, 3.0]);
    assert!(
        matches!(res, Err(LtdlError::DimensionMismatch { .. })),
        "wrong-length b must yield DimensionMismatch, got {res:?}"
    );
}

#[test]
fn test_not_positive_definite_errors() {
    let res = SparseMassFactorization::factor(vec![vec![-1.0, 0.0], vec![0.0, 1.0]], vec![0, 0, 0], 2);
    let err = res.expect_err("factor of indefinite matrix must error");
    assert_eq!(err, LtdlError::NotPositiveDefinite { index: 0 });
}

#[test]
fn test_sparse_beats_dense_cholesky_2x() {
    // PERF GATE: sparse LᵀDL must be ≥ 2x faster than dense Cholesky on a real
    // deep + branched tree. Spine(6) + 5 limbs × 6 = 36 revolute DOFs.
    let mut model = build_perf_tree(6, 5, 6);
    let n = model.total_dof();
    let q = vec![0.0; n];
    let h = compute_mass_matrix_crba(&mut model, &q);
    let lambda = build_lambda(&model);
    let b: Vec<f64> = (0..n).map(|i| (i as f64 + 1.0) * 0.13 + 0.41).collect();

    // Warm up both paths once (also validates correctness / SPD).
    {
        let l = dense_cholesky_factor(&h);
        let xd = dense_cholesky_solve_vec(&l, &b);
        let fact = SparseMassFactorization::factor(h.clone(), lambda.clone(), n)
            .expect("perf tree must be SPD");
        let xs = fact.solve(&b).expect("solve");
        // Sanity: both solvers must agree on this SPD system.
        for i in 0..n {
            assert!(
                (xd[i] - xs[i]).abs() < 1e-8,
                "warmup mismatch at {i}: dense {} vs sparse {}",
                xd[i],
                xs[i]
            );
        }
        std::hint::black_box((xd, xs));
    }

    let iters = 2000usize;
    let mut best_ratio = 0.0f64;
    let mut checksum = 0.0f64;

    for _attempt in 0..4 {
        // Time dense: clone → factor → solve.
        let t0 = std::time::Instant::now();
        for _ in 0..iters {
            let l = dense_cholesky_factor(&h);
            let x = dense_cholesky_solve_vec(&l, &b);
            checksum += std::hint::black_box(x[0]);
        }
        let dense_elapsed = t0.elapsed().as_secs_f64();

        // Time sparse: factor (owned clone) → solve.
        let t1 = std::time::Instant::now();
        for _ in 0..iters {
            let fact = SparseMassFactorization::factor(h.clone(), lambda.clone(), n)
                .expect("perf tree must be SPD");
            let x = fact.solve(&b).expect("solve");
            checksum += std::hint::black_box(x[0]);
        }
        let sparse_elapsed = t1.elapsed().as_secs_f64();

        let ratio = dense_elapsed / sparse_elapsed;
        if ratio > best_ratio {
            best_ratio = ratio;
        }
        if best_ratio >= 2.0 {
            break;
        }
    }

    // Read the checksum after timing so the optimizer cannot elide the work.
    std::hint::black_box(checksum);

    eprintln!("dense/sparse perf ratio = {best_ratio} (n = {n})");
    assert!(
        best_ratio >= 2.0,
        "sparse LᵀDL must be ≥ 2x faster than dense Cholesky (n = {n}); best ratio = {best_ratio}"
    );
}
