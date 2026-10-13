// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the single-sweep Dirichlet elimination used by
//! `boundary::apply_dirichlet`, `assembly::apply_dirichlet_bc` and
//! `assembly::apply_dirichlet_bc_values`.
//!
//! The pre-fix implementations rescanned the whole CSR matrix for every
//! constrained DOF (`O(n_bc * nnz)`) and panicked on out-of-range DOFs.  The
//! old sequential algorithm is kept below as a test oracle: the new code must
//! reproduce its output exactly (up to floating-point summation order) on
//! well-posed input, and must report errors instead of panicking otherwise.

use oxiphysics_core::math::Vec3;
use oxiphysics_fem::Error;
use oxiphysics_fem::assembly::{
    LinearElasticAssembler, apply_dirichlet_bc, apply_dirichlet_bc_values, assemble_force_vector,
    assemble_global_stiffness,
};
use oxiphysics_fem::boundary::{
    DirichletBc, TimeDependentBC, TimeFunction, apply_dirichlet, apply_time_dependent_sparse,
};
use oxiphysics_fem::constitutive::LinearElasticMaterial;
use oxiphysics_fem::mesh::{cantilever_fixed_dofs, cantilever_tip_nodes, generate_cantilever_mesh};
use oxiphysics_fem::solvers::{jacobi_preconditioner, pcg_solve};
use oxiphysics_fem::sparse::CsrMatrix;

// ── helpers ──────────────────────────────────────────────────────────────────

/// Deterministic pseudo-random number in `[0, 1)` (64-bit LCG).
fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// Symmetric, strictly diagonally dominant (hence SPD) banded matrix.
fn spd_matrix(n: usize, seed: u64) -> CsrMatrix {
    let mut state = seed;
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    let mut off_diag_abs_sum = vec![0.0; n];
    for i in 0..n {
        for j in (i + 1)..n.min(i + 6) {
            let v = -(0.2 + lcg(&mut state));
            triplets.push((i, j, v));
            triplets.push((j, i, v));
            off_diag_abs_sum[i] += v.abs();
            off_diag_abs_sum[j] += v.abs();
        }
    }
    for (i, s) in off_diag_abs_sum.iter().enumerate() {
        triplets.push((i, i, s + 1.0 + (i % 3) as f64));
    }
    CsrMatrix::from_triplets(n, n, &triplets)
}

/// The pre-fix sequential `O(n_fixed * nnz)` algorithm, kept as an oracle.
fn reference_apply(k: &mut CsrMatrix, f: &mut [f64], fixed: &[(usize, f64)]) {
    for &(dof, val) in fixed {
        if val.abs() > 0.0 {
            for (row, f_row) in f.iter_mut().enumerate().take(k.nrows) {
                if row == dof {
                    continue;
                }
                for idx in k.row_ptr[row]..k.row_ptr[row + 1] {
                    if k.col_indices[idx] == dof {
                        *f_row -= k.values[idx] * val;
                    }
                }
            }
        }
        for idx in k.row_ptr[dof]..k.row_ptr[dof + 1] {
            k.values[idx] = if k.col_indices[idx] == dof { 1.0 } else { 0.0 };
        }
        for row in 0..k.nrows {
            if row == dof {
                continue;
            }
            for idx in k.row_ptr[row]..k.row_ptr[row + 1] {
                if k.col_indices[idx] == dof {
                    k.values[idx] = 0.0;
                }
            }
        }
        f[dof] = val;
    }
}

fn assert_close(a: f64, b: f64, what: &str) {
    assert!(
        (a - b).abs() <= 1e-10 * (1.0 + a.abs().max(b.abs())),
        "{what}: {a} vs {b}"
    );
}

fn assert_same_system(k: &CsrMatrix, f: &[f64], k_ref: &CsrMatrix, f_ref: &[f64]) {
    assert_eq!(k.row_ptr, k_ref.row_ptr);
    assert_eq!(k.col_indices, k_ref.col_indices);
    for (idx, (a, b)) in k.values.iter().zip(&k_ref.values).enumerate() {
        assert_close(*a, *b, &format!("K value #{idx}"));
    }
    for (i, (a, b)) in f.iter().zip(f_ref).enumerate() {
        assert_close(*a, *b, &format!("f[{i}]"));
    }
}

fn ten_bcs() -> Vec<DirichletBc> {
    vec![
        DirichletBc::new(0, 0, 0.0),
        DirichletBc::new(0, 1, 0.3),
        DirichletBc::new(2, 2, -0.7),
        DirichletBc::new(3, 0, 1.1),
        DirichletBc::new(3, 1, 0.0),
        DirichletBc::new(7, 2, 0.25),
        DirichletBc::new(8, 0, -0.4),
        DirichletBc::new(8, 1, 0.9),
        DirichletBc::new(12, 2, 0.6),
        DirichletBc::new(16, 1, -1.3),
    ]
}

fn cantilever_system() -> (CsrMatrix, Vec<f64>, Vec<usize>) {
    let (length, width, height) = (1.0_f64, 0.2_f64, 0.2_f64);
    let mesh = generate_cantilever_mesh(length, width, height, 4, 2, 2);
    let d = LinearElasticMaterial::new(200.0e9, 0.3).constitutive_matrix();
    let k = assemble_global_stiffness(&mesh.nodes, &mesh.elements, &d);
    let tip = cantilever_tip_nodes(&mesh, length, 1e-9);
    let n_tip = tip.len() as f64;
    let forces: Vec<(usize, Vec3)> = tip
        .iter()
        .map(|&i| (i, Vec3::new(0.0, -1000.0 / n_tip, 0.0)))
        .collect();
    let f = assemble_force_vector(mesh.n_nodes(), &forces);
    let fixed = cantilever_fixed_dofs(&mesh, 1e-9);
    (k, f, fixed)
}

// ── boundary::apply_dirichlet ────────────────────────────────────────────────

#[test]
fn matches_sequential_reference_on_51_dof_spd_matrix_with_ten_bcs() {
    let n = 51;
    let k0 = spd_matrix(n, 7);
    let f0: Vec<f64> = (0..n).map(|i| 1.0 + 0.1 * i as f64).collect();
    let bcs = ten_bcs();

    let mut k_new = k0.clone();
    let mut f_new = f0.clone();
    apply_dirichlet(&mut k_new, &mut f_new, &bcs).expect("valid BCs");

    let mut k_ref = k0.clone();
    let mut f_ref = f0.clone();
    let prescribed: Vec<(usize, f64)> = bcs.iter().map(|b| (b.global_dof(), b.value)).collect();
    reference_apply(&mut k_ref, &mut f_ref, &prescribed);

    assert_same_system(&k_new, &f_new, &k_ref, &f_ref);
}

#[test]
fn eliminated_system_reproduces_prescribed_values_and_free_equilibrium() {
    let n = 51;
    let k0 = spd_matrix(n, 11);
    let f0: Vec<f64> = (0..n).map(|i| ((i * 7 % 5) as f64) - 2.0).collect();
    let bcs = ten_bcs();

    let mut k = k0.clone();
    let mut f = f0.clone();
    apply_dirichlet(&mut k, &mut f, &bcs).expect("valid BCs");
    let precond = jacobi_preconditioner(&k);
    let (u, stats) = pcg_solve(&k, &f, &precond, 2000, 1e-14);
    assert!(stats.converged, "PCG did not converge: {stats:?}");

    let fixed: Vec<usize> = bcs.iter().map(DirichletBc::global_dof).collect();
    for bc in &bcs {
        assert_close(u[bc.global_dof()], bc.value, "prescribed value");
    }
    // The *original* equations must hold on every free DOF.
    let residual = k0.mul_vec(&u);
    for i in (0..n).filter(|i| !fixed.contains(i)) {
        assert!(
            (residual[i] - f0[i]).abs() < 1e-9,
            "free equilibrium violated at DOF {i}: {} vs {}",
            residual[i],
            f0[i]
        );
    }
}

#[test]
fn out_of_range_dof_is_an_error_and_leaves_the_system_untouched() {
    let k0 = spd_matrix(6, 3);
    let f0 = vec![1.0; 6];
    let mut k = k0.clone();
    let mut f = f0.clone();
    // Valid BC first, invalid (global DOF 3 * 2 + 0 = 6) second: nothing may change.
    let bcs = [DirichletBc::new(0, 0, 0.5), DirichletBc::new(2, 0, 0.0)];
    let result = apply_dirichlet(&mut k, &mut f, &bcs);
    assert!(
        matches!(result, Err(Error::DofOutOfRange { dof: 6, n: 6 })),
        "unexpected result: {result:?}"
    );
    assert_eq!(k.values, k0.values);
    assert_eq!(f, f0);
}

#[test]
fn local_dof_above_two_is_rejected_instead_of_aliasing_the_next_node() {
    let k0 = spd_matrix(6, 5);
    let mut k = k0.clone();
    let mut f = vec![0.0; 6];
    let result = apply_dirichlet(&mut k, &mut f, &[DirichletBc::new(0, 4, 1.0)]);
    assert!(
        matches!(result, Err(Error::InvalidInput(_))),
        "unexpected result: {result:?}"
    );
    assert_eq!(k.values, k0.values);
}

#[test]
fn huge_node_index_does_not_overflow() {
    let mut k = spd_matrix(3, 1);
    let mut f = vec![0.0; 3];
    let result = apply_dirichlet(&mut k, &mut f, &[DirichletBc::new(usize::MAX, 0, 0.0)]);
    assert!(
        matches!(result, Err(Error::DofOutOfRange { .. })),
        "unexpected result: {result:?}"
    );
}

#[test]
fn identical_duplicates_are_fine_conflicting_duplicates_are_rejected() {
    let k0 = spd_matrix(6, 9);
    let f0 = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

    let mut k_once = k0.clone();
    let mut f_once = f0.clone();
    apply_dirichlet(&mut k_once, &mut f_once, &[DirichletBc::new(0, 1, 0.75)]).expect("ok");
    let mut k_twice = k0.clone();
    let mut f_twice = f0.clone();
    let twice = [DirichletBc::new(0, 1, 0.75), DirichletBc::new(0, 1, 0.75)];
    apply_dirichlet(&mut k_twice, &mut f_twice, &twice).expect("identical duplicates are legal");
    assert_same_system(&k_twice, &f_twice, &k_once, &f_once);

    let mut k = k0.clone();
    let mut f = f0.clone();
    let conflicting = [DirichletBc::new(0, 1, 0.75), DirichletBc::new(0, 1, -0.75)];
    let result = apply_dirichlet(&mut k, &mut f, &conflicting);
    assert!(
        matches!(result, Err(Error::InvalidInput(_))),
        "unexpected result: {result:?}"
    );
    assert_eq!(k.values, k0.values);
    assert_eq!(f, f0);
}

#[test]
fn non_finite_prescribed_values_are_rejected() {
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut k = spd_matrix(3, 2);
        let mut f = vec![0.0; 3];
        let result = apply_dirichlet(&mut k, &mut f, &[DirichletBc::new(0, 0, bad)]);
        assert!(
            matches!(result, Err(Error::InvalidInput(_))),
            "value {bad} gave {result:?}"
        );
    }
}

#[test]
fn size_mismatches_are_rejected() {
    let mut k = spd_matrix(6, 4);
    let mut short_f = vec![0.0; 5];
    let result = apply_dirichlet(&mut k, &mut short_f, &[DirichletBc::new(0, 0, 0.0)]);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");

    let mut rect = CsrMatrix::from_triplets(3, 4, &[(0, 0, 1.0), (1, 1, 1.0), (2, 2, 1.0)]);
    let mut f = vec![0.0; 3];
    let result = apply_dirichlet(&mut rect, &mut f, &[DirichletBc::new(0, 0, 0.0)]);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");

    // Inconsistent CSR arrays are reported instead of indexing out of bounds.
    let mut broken = spd_matrix(3, 6);
    broken.row_ptr.pop();
    let mut f = vec![0.0; 3];
    let result = apply_dirichlet(&mut broken, &mut f, &[DirichletBc::new(0, 0, 0.0)]);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");
    let mut truncated = spd_matrix(3, 6);
    truncated.values.pop();
    let result = apply_dirichlet(&mut truncated, &mut f, &[DirichletBc::new(0, 0, 0.0)]);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");

    // A row pointer that jumps past the stored entries and comes back would
    // index out of bounds mid-sweep; it is rejected up front.
    let mut jumping = spd_matrix(3, 6);
    jumping.row_ptr[1] = jumping.values.len() + 7;
    let before = jumping.values.clone();
    let result = apply_dirichlet(&mut jumping, &mut f, &[DirichletBc::new(0, 0, 0.0)]);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");
    assert_eq!(jumping.values, before);
    let mut shifted = spd_matrix(3, 6);
    shifted.row_ptr[0] = 1;
    let result = apply_dirichlet(&mut shifted, &mut f, &[DirichletBc::new(0, 0, 0.0)]);
    assert!(matches!(result, Err(Error::InvalidInput(_))), "{result:?}");
}

#[test]
fn fixed_dof_without_stored_diagonal_still_becomes_an_identity_row() {
    // DOF 1 is an orphan (no element touches it): its CSR row is empty.
    let triplets = [(0, 0, 2.0), (0, 2, -1.0), (2, 0, -1.0), (2, 2, 2.0)];
    let mut k = CsrMatrix::from_triplets(3, 3, &triplets);
    assert_eq!(k.nnz(), 4);
    let mut f = vec![1.0, 5.0, 1.0];
    apply_dirichlet(&mut k, &mut f, &[DirichletBc::new(0, 1, 0.5)]).expect("valid BC");
    assert_eq!(k.nnz(), 5, "the missing diagonal must be inserted");
    assert_eq!(k.get(1, 1), 1.0);
    assert_eq!(f[1], 0.5);
    assert_eq!(k.get(0, 0), 2.0);
    assert_eq!(k.get(2, 2), 2.0);
    assert_eq!(f[0], 1.0);
    assert_eq!(f[2], 1.0);
}

#[test]
fn empty_bc_list_is_a_no_op() {
    let k0 = spd_matrix(6, 8);
    let mut k = k0.clone();
    let mut f = vec![1.0; 6];
    apply_dirichlet(&mut k, &mut f, &[]).expect("no BCs");
    assert_eq!(k.values, k0.values);
    assert_eq!(f, vec![1.0; 6]);
}

#[test]
fn time_dependent_sparse_applies_the_value_at_time_t_and_propagates_errors() {
    let k0 = spd_matrix(6, 13);
    let f0 = vec![1.0; 6];
    let bcs = [TimeDependentBC::new(
        0,
        1,
        2.0,
        TimeFunction::Ramp { ramp_time: 4.0 },
    )];
    let mut k = k0.clone();
    let mut f = f0.clone();
    apply_time_dependent_sparse(&mut k, &mut f, &bcs, 1.0).expect("valid BC");
    assert_close(f[1], 0.5, "ramp(1/4) * amplitude 2");
    assert_eq!(k.get(1, 1), 1.0);

    let bad = [TimeDependentBC::new(5, 0, 1.0, TimeFunction::Constant(1.0))];
    let mut k = k0.clone();
    let mut f = f0.clone();
    let result = apply_time_dependent_sparse(&mut k, &mut f, &bad, 0.0);
    assert!(
        matches!(result, Err(Error::DofOutOfRange { .. })),
        "{result:?}"
    );
    assert_eq!(k.values, k0.values);
}

// ── assembly::apply_dirichlet_bc / apply_dirichlet_bc_values ─────────────────

#[test]
fn apply_dirichlet_bc_matches_full_sweep_reference_on_cantilever_mesh() {
    let (k0, f0, fixed) = cantilever_system();
    assert!(fixed.len() > 20, "the clamped face must fix many DOFs");

    let mut k = k0.clone();
    let mut f = f0.clone();
    apply_dirichlet_bc(&mut k, &mut f, &fixed).expect("valid fixed DOFs");

    let mut k_ref = k0.clone();
    let mut f_ref = f0.clone();
    let zero_values: Vec<(usize, f64)> = fixed.iter().map(|&d| (d, 0.0)).collect();
    reference_apply(&mut k_ref, &mut f_ref, &zero_values);

    assert_same_system(&k, &f, &k_ref, &f_ref);
}

#[test]
fn apply_dirichlet_bc_values_matches_sequential_reference_with_nonzero_values() {
    let (k0, f0, fixed) = cantilever_system();
    let prescribed: Vec<(usize, f64)> = fixed
        .iter()
        .enumerate()
        .map(|(i, &d)| (d, 1.0e-4 * ((i % 7) as f64 - 3.0)))
        .collect();

    let mut k = k0.clone();
    let mut f = f0.clone();
    apply_dirichlet_bc_values(&mut k, &mut f, &prescribed).expect("valid DOFs");

    let mut k_ref = k0.clone();
    let mut f_ref = f0.clone();
    reference_apply(&mut k_ref, &mut f_ref, &prescribed);

    assert_same_system(&k, &f, &k_ref, &f_ref);
}

#[test]
fn assembly_out_of_range_dof_is_an_error_not_a_panic_and_leaves_the_system_untouched() {
    let (k0, f0, fixed) = cantilever_system();
    let n = k0.nrows;

    let mut dofs = fixed.clone();
    dofs.push(n); // one past the end, listed after valid DOFs
    let mut k = k0.clone();
    let mut f = f0.clone();
    let result = apply_dirichlet_bc(&mut k, &mut f, &dofs);
    assert!(
        matches!(result, Err(Error::DofOutOfRange { dof, n: size }) if dof == n && size == n),
        "unexpected result: {result:?}"
    );
    assert_eq!(k.values, k0.values);
    assert_eq!(f, f0);

    let result = apply_dirichlet_bc_values(&mut k, &mut f, &[(0, 1.0), (n + 5, 2.0)]);
    assert!(
        matches!(result, Err(Error::DofOutOfRange { .. })),
        "unexpected result: {result:?}"
    );
    assert_eq!(k.values, k0.values);
    assert_eq!(f, f0);
}

#[test]
fn assembly_invalid_values_and_shapes_are_rejected() {
    let (k0, f0, _) = cantilever_system();
    let mut k = k0.clone();
    let mut f = f0.clone();
    let nan = apply_dirichlet_bc_values(&mut k, &mut f, &[(0, f64::NAN)]);
    assert!(matches!(nan, Err(Error::InvalidInput(_))), "{nan:?}");
    let conflict = apply_dirichlet_bc_values(&mut k, &mut f, &[(3, 1.0), (3, 2.0)]);
    assert!(
        matches!(conflict, Err(Error::InvalidInput(_))),
        "{conflict:?}"
    );
    let mut short = vec![0.0; k0.nrows - 1];
    let mismatch = apply_dirichlet_bc(&mut k, &mut short, &[0]);
    assert!(
        matches!(mismatch, Err(Error::InvalidInput(_))),
        "{mismatch:?}"
    );
    assert_eq!(k.values, k0.values);
    assert_eq!(f, f0);
}

#[test]
fn repeating_a_fixed_dof_is_harmless() {
    let (k0, f0, fixed) = cantilever_system();
    let mut doubled = fixed.clone();
    doubled.extend_from_slice(&fixed);

    let mut k_once = k0.clone();
    let mut f_once = f0.clone();
    apply_dirichlet_bc(&mut k_once, &mut f_once, &fixed).expect("valid");
    let mut k_twice = k0.clone();
    let mut f_twice = f0.clone();
    apply_dirichlet_bc(&mut k_twice, &mut f_twice, &doubled).expect("valid");
    assert_eq!(k_once.values, k_twice.values);
    assert_eq!(f_once, f_twice);
}

#[test]
fn assembler_wrapper_applies_and_propagates_the_error() {
    let (k0, f0, fixed) = cantilever_system();
    let nodes: Vec<Vec3> = Vec::new();
    let elements: Vec<[usize; 4]> = Vec::new();
    let assembler = LinearElasticAssembler::new(&nodes, &elements, [[0.0; 6]; 6]);

    let mut k = k0.clone();
    let mut f = f0.clone();
    assembler
        .apply_dirichlet(&mut k, &mut f, &fixed)
        .expect("valid fixed DOFs");
    for &dof in &fixed {
        assert_eq!(k.get(dof, dof), 1.0);
        assert_eq!(f[dof], 0.0);
    }

    let mut k = k0.clone();
    let mut f = f0.clone();
    let result = assembler.apply_dirichlet(&mut k, &mut f, &[k0.nrows]);
    assert!(
        matches!(result, Err(Error::DofOutOfRange { .. })),
        "unexpected result: {result:?}"
    );
}

// ── complexity ───────────────────────────────────────────────────────────────

/// Regression for the `O(n_bc * nnz)` sweep: with 100 000 constrained DOFs of a
/// 200 000-DOF system the old algorithm needs about 6e10 inner steps (minutes);
/// the single sweep needs a few milliseconds.  The generous 20 s bound only
/// trips for the quadratic behaviour.
#[test]
fn cost_is_linear_in_nnz_not_nnz_times_constraints() {
    let n = 200_000;
    let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity(3 * n);
    for i in 0..n {
        triplets.push((i, i, 4.0));
        if i + 1 < n {
            triplets.push((i, i + 1, -1.0));
            triplets.push((i + 1, i, -1.0));
        }
    }
    let mut k = CsrMatrix::from_triplets(n, n, &triplets);
    let mut f = vec![1.0; n];
    // Every even DOF is constrained to 0.5.
    let prescribed: Vec<(usize, f64)> = (0..n).step_by(2).map(|d| (d, 0.5)).collect();
    assert_eq!(prescribed.len(), 100_000);

    let start = std::time::Instant::now();
    apply_dirichlet_bc_values(&mut k, &mut f, &prescribed).expect("valid DOFs");
    let elapsed = start.elapsed();
    assert!(
        elapsed < std::time::Duration::from_secs(20),
        "elimination took {elapsed:?}"
    );

    // Free (odd) DOFs lose both neighbours: f = 1 - (-1 * 0.5) * 2 = 2 for
    // interior ones; the last DOF (n - 1 is odd) only has one neighbour.
    assert_close(f[1], 2.0, "interior free DOF");
    assert_close(f[n - 1], 1.5, "last free DOF with one fixed neighbour");
    assert_eq!(f[0], 0.5);
    assert_eq!(k.get(1, 0), 0.0);
    assert_eq!(k.get(0, 1), 0.0);
    assert_eq!(k.get(1, 1), 4.0);
}
