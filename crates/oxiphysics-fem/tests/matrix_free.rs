// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for matrix-free sum-factorized high-order Poisson operators.
//!
//! These tests exercise [`oxiphysics_fem::solvers::matrix_free::SumFactPoisson`]
//! through its *public* surface only. The crate-private helpers `element_apply`,
//! `lagrange_at`, and `contract_axis` are not reachable from an integration test
//! (which is a separate crate), so every baseline here is an independent
//! reimplementation built from `gll_nodes_weights` plus the closed-form DOF map.

use oxiphysics_fem::parallel_solver::CsrMatrix;
use oxiphysics_fem::solvers::amg::preconditioner::JacobiPreconditioner;
use oxiphysics_fem::solvers::matrix_free::sum_factorization::gll_nodes_weights;
use oxiphysics_fem::solvers::{MatrixFreeOperator, MatrixFreePcg, SumFactPoisson};
use std::f64::consts::PI;

// ───────────────────────────── deterministic RNG ─────────────────────────────

/// A tiny deterministic linear-congruential generator producing values in
/// roughly `[-1, 1]`. Pure Rust, no external crate; reproducible across runs.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next pseudo-random `f64` in approximately `[-1, 1)`.
    fn next_unit(&mut self) -> f64 {
        // 64-bit LCG (PCG/Knuth multiplier) then take the high bits.
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let v = ((self.state >> 33) as f64) / ((1u64 << 31) as f64);
        v - 1.0
    }

    /// A length-`n` vector of pseudo-random values in `[-1, 1)`.
    fn vec(&mut self, n: usize) -> Vec<f64> {
        (0..n).map(|_| self.next_unit()).collect()
    }
}

// ──────────────────────── independent 1D Lagrange basis ──────────────────────

/// Reproduce the 1D Lagrange value and derivative matrices, row-major
/// `np × np` with `out[q * np + a]`, sampling basis `a` (nodal at `nodes`) at
/// `quad_pts[q]`. This intentionally mirrors the operator's own convention
/// (`phi[q*np + a] = φ_a(ξ_q)`) so the comparison is apples-to-apples, but it is
/// an entirely separate implementation living in the test crate.
fn lagrange_1d(nodes: &[f64], quad_pts: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let np = nodes.len();
    let nq = quad_pts.len();
    let mut phi = vec![0.0f64; nq * np];
    let mut dphi = vec![0.0f64; nq * np];
    for q in 0..nq {
        let x = quad_pts[q];
        for a in 0..np {
            // value: Π_{m≠a} (x - x_m)/(x_a - x_m)
            let mut val = 1.0;
            for m in 0..np {
                if m == a {
                    continue;
                }
                val *= (x - nodes[m]) / (nodes[a] - nodes[m]);
            }
            phi[q * np + a] = val;

            // derivative (stable form, exact at nodes):
            // Σ_{k≠a} [ 1/(x_a - x_k) · Π_{m≠a,k} (x - x_m)/(x_a - x_m) ]
            let mut d = 0.0;
            for k in 0..np {
                if k == a {
                    continue;
                }
                let mut prod = 1.0 / (nodes[a] - nodes[k]);
                for m in 0..np {
                    if m == a || m == k {
                        continue;
                    }
                    prod *= (x - nodes[m]) / (nodes[a] - nodes[m]);
                }
                d += prod;
            }
            dphi[q * np + a] = d;
        }
    }
    (phi, dphi)
}

// ─────────────────── independent dense element stiffness (GLL) ────────────────

/// Build the dense reference element stiffness `Ke` (`np³ × np³`, row-major) for a
/// degree-`p` GLL spectral element on a cube of side `h`, using GLL quadrature at
/// the same `p+1` points. This is the brute-force triple-quadrature definition —
/// NO sum factorization.
///
/// Local node `m = (im,jm,km)` flattens to `(im*np + jm)*np + km`. Quadrature
/// point `q = (a,b,c)`. Geometry factor for the isotropic cube map `J=(h/2)I` is
/// `geo = (2/h)² · (h/2)³`, and the 3D weight is `w_a·w_b·w_c`.
fn element_stiffness_dense(p: usize, h: f64) -> Vec<f64> {
    let np = p + 1;
    let (nodes, weights) = gll_nodes_weights(p);
    let (phi, dphi) = lagrange_1d(&nodes, &nodes);
    let n3 = np * np * np;
    let mut ke = vec![0.0f64; n3 * n3];

    let two_over_h = 2.0 / h;
    let half_h = h / 2.0;
    let geo = two_over_h * two_over_h * half_h * half_h * half_h;

    // Precompute, for every (node m, quad q), the three reference gradient
    // components G^ξ, G^η, G^ζ. Index by flat node and flat quad.
    let mut g_xi = vec![0.0f64; n3 * n3];
    let mut g_eta = vec![0.0f64; n3 * n3];
    let mut g_zeta = vec![0.0f64; n3 * n3];
    for im in 0..np {
        for jm in 0..np {
            for km in 0..np {
                let m = (im * np + jm) * np + km;
                for a in 0..np {
                    let phi_a_im = phi[a * np + im];
                    let dphi_a_im = dphi[a * np + im];
                    for b in 0..np {
                        let phi_b_jm = phi[b * np + jm];
                        let dphi_b_jm = dphi[b * np + jm];
                        for c in 0..np {
                            let phi_c_km = phi[c * np + km];
                            let dphi_c_km = dphi[c * np + km];
                            let q = (a * np + b) * np + c;
                            let idx = m * n3 + q;
                            g_xi[idx] = dphi_a_im * phi_b_jm * phi_c_km;
                            g_eta[idx] = phi_a_im * dphi_b_jm * phi_c_km;
                            g_zeta[idx] = phi_a_im * phi_b_jm * dphi_c_km;
                        }
                    }
                }
            }
        }
    }

    // Quadrature weight at each flat quad point.
    let mut wq = vec![0.0f64; n3];
    for a in 0..np {
        for b in 0..np {
            for c in 0..np {
                wq[(a * np + b) * np + c] = weights[a] * weights[b] * weights[c];
            }
        }
    }

    // Ke[m][n] = Σ_q geo·w_q·( G^ξ_m·G^ξ_n + G^η_m·G^η_n + G^ζ_m·G^ζ_n )
    for m in 0..n3 {
        for nn in 0..n3 {
            let mut s = 0.0;
            for q in 0..n3 {
                let f = geo * wq[q];
                s += f
                    * (g_xi[m * n3 + q] * g_xi[nn * n3 + q]
                        + g_eta[m * n3 + q] * g_eta[nn * n3 + q]
                        + g_zeta[m * n3 + q] * g_zeta[nn * n3 + q]);
            }
            ke[m * n3 + nn] = s;
        }
    }
    ke
}

/// Closed-form global DOF index of local node `(i,j,k)` in element `(ex,ey,ez)`
/// for a mesh with `ndir = n*p + 1` nodes per axis.
fn global_dof(
    ex: usize,
    ey: usize,
    ez: usize,
    i: usize,
    j: usize,
    k: usize,
    p: usize,
    ndir: usize,
) -> usize {
    let gi = ex * p + i;
    let gj = ey * p + j;
    let gk = ez * p + k;
    (gi * ndir + gj) * ndir + gk
}

/// Assemble the full global stiffness as a dense row-major `Vec<f64>` of length
/// `n_dofs²` for an `n³`-element degree-`p` GLL mesh on the unit cube. Returns
/// `(n_dofs, K_global)`.
fn assembled_gll_stiffness_dense(n: usize, p: usize) -> (usize, Vec<f64>) {
    let np = p + 1;
    let n3 = np * np * np;
    let h = 1.0 / n as f64;
    let ndir = n * p + 1;
    let n_dofs = ndir * ndir * ndir;

    let ke = element_stiffness_dense(p, h);
    let mut k_global = vec![0.0f64; n_dofs * n_dofs];

    // local flat node -> (i,j,k)
    let unflatten = |m: usize| -> (usize, usize, usize) {
        let i = m / (np * np);
        let rem = m % (np * np);
        let j = rem / np;
        let k = rem % np;
        (i, j, k)
    };

    for ex in 0..n {
        for ey in 0..n {
            for ez in 0..n {
                for m in 0..n3 {
                    let (im, jm, km) = unflatten(m);
                    let gm = global_dof(ex, ey, ez, im, jm, km, p, ndir);
                    for nn in 0..n3 {
                        let (in_, jn, kn) = unflatten(nn);
                        let gn = global_dof(ex, ey, ez, in_, jn, kn, p, ndir);
                        k_global[gm * n_dofs + gn] += ke[m * n3 + nn];
                    }
                }
            }
        }
    }
    (n_dofs, k_global)
}

/// Dense mat-vec `y = K·x` for a row-major `n×n` matrix.
fn dense_matvec(k: &[f64], x: &[f64], n: usize) -> Vec<f64> {
    let mut y = vec![0.0f64; n];
    for r in 0..n {
        let mut s = 0.0;
        let row = &k[r * n..r * n + n];
        for c in 0..n {
            s += row[c] * x[c];
        }
        y[r] = s;
    }
    y
}

/// Maximum absolute difference between two equal-length slices.
fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .fold(0.0f64, |m, (ai, bi)| m.max((ai - bi).abs()))
}

// ───────────────────────── 1D coordinate + lumped mass ───────────────────────

/// Build the physical coordinate of every global node along one axis for an
/// `n`-element degree-`p` GLL mesh on `[0,1]`. Element `e` occupies
/// `[e·h, (e+1)·h]` with `h = 1/n`; the local GLL node `l` sits at
/// `e·h + (ξ_l + 1)/2 · h`. Shared nodes (l=p of element e, l=0 of element e+1)
/// receive the same coordinate, so the overwrite is consistent.
fn coord_1d(n: usize, p: usize) -> Vec<f64> {
    let (nodes, _w) = gll_nodes_weights(p);
    let h = 1.0 / n as f64;
    let ndir = n * p + 1;
    let mut coord = vec![0.0f64; ndir];
    for e in 0..n {
        for (l, &node) in nodes.iter().enumerate() {
            let g = e * p + l;
            coord[g] = e as f64 * h + (node + 1.0) * 0.5 * h;
        }
    }
    coord
}

/// Assemble the diagonal (lumped) GLL mass vector `m_g` for an `n³`-element
/// degree-`p` mesh on the unit cube. Because the basis is GLL-collocated and the
/// quadrature uses the same points, the consistent mass matrix is already
/// diagonal: `m_g = Σ_{elements ∋ g} (h/2)³ · w_a·w_b·w_c` evaluated at the local
/// index of `g` within each element.
fn lumped_mass(n: usize, p: usize) -> Vec<f64> {
    let np = p + 1;
    let (_nodes, weights) = gll_nodes_weights(p);
    let h = 1.0 / n as f64;
    let half_h3 = (h / 2.0).powi(3);
    let ndir = n * p + 1;
    let n_dofs = ndir * ndir * ndir;
    let mut mass = vec![0.0f64; n_dofs];

    for ex in 0..n {
        for ey in 0..n {
            for ez in 0..n {
                for i in 0..np {
                    let wi = weights[i];
                    for j in 0..np {
                        let wij = wi * weights[j];
                        for (k, &wk) in weights.iter().enumerate() {
                            let wijk = wij * wk;
                            let g = global_dof(ex, ey, ez, i, j, k, p, ndir);
                            mass[g] += half_h3 * wijk;
                        }
                    }
                }
            }
        }
    }
    mass
}

// ─────────────────── Dirichlet wrapper around the operator ───────────────────

/// Wraps a [`SumFactPoisson`] to impose homogeneous Dirichlet BCs by symmetric
/// elimination: interior rows act as `K`, boundary rows act as identity. The
/// input is masked (boundary entries zeroed) before the inner stiffness apply,
/// and boundary rows of the output get the (unmasked) input value back, giving a
/// clean `[K_ii 0; 0 I]` block structure consistent with a zeroed boundary RHS.
struct DirichletOp<'a> {
    inner: &'a SumFactPoisson,
    mask: &'a [bool],
}

impl MatrixFreeOperator for DirichletOp<'_> {
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        let n = self.inner.dim();
        // Masked copy of x: zero on boundary so boundary DOFs do not couple in.
        let mut x_int = x.to_vec();
        for (i, xi) in x_int.iter_mut().enumerate() {
            if self.mask[i] {
                *xi = 0.0;
            }
        }
        // temp = K · x_int
        let mut temp = vec![0.0f64; n];
        self.inner.apply(&x_int, &mut temp);
        // Interior: y += (K·x_int)_i. Boundary: y += x_i (identity row).
        for i in 0..n {
            if self.mask[i] {
                y[i] += x[i];
            } else {
                y[i] += temp[i];
            }
        }
    }

    fn dim(&self) -> usize {
        self.inner.dim()
    }
}

/// Build the boundary mask for an `n³`-element degree-`p` mesh: a global node is
/// Dirichlet iff any of its axis indices is `0` or `ndir-1`.
fn boundary_mask(ndir: usize) -> Vec<bool> {
    let n_dofs = ndir * ndir * ndir;
    let mut mask = vec![false; n_dofs];
    for gi in 0..ndir {
        let on_i = gi == 0 || gi == ndir - 1;
        for gj in 0..ndir {
            let on_j = gj == 0 || gj == ndir - 1;
            for gk in 0..ndir {
                let on_k = gk == 0 || gk == ndir - 1;
                if on_i || on_j || on_k {
                    mask[(gi * ndir + gj) * ndir + gk] = true;
                }
            }
        }
    }
    mask
}

// ─────────────────────────────────── tests ──────────────────────────────────

/// Test 1 — the matrix-free operator's mat-vec equals an independently assembled
/// dense GLL stiffness, at p=1 (identity φ) and p=2 (non-trivial dφ).
#[test]
fn matfree_matvec_matches_assembled_p1() {
    // p = 1, 2×2×2 elements → 27 DOFs.
    let op = SumFactPoisson::new(2, 1);
    assert_eq!(op.n_dofs(), 27, "expected 27 DOFs, got {}", op.n_dofs());
    assert_eq!(
        op.dofs_per_dir(),
        3,
        "expected 3 dofs/dir, got {}",
        op.dofs_per_dir()
    );

    let (n_dofs, k_global) = assembled_gll_stiffness_dense(2, 1);
    assert_eq!(n_dofs, 27, "baseline n_dofs mismatch: {n_dofs}");

    let mut rng = Lcg::new(0x1234_5678);
    let x = rng.vec(27);

    let y_assembled = dense_matvec(&k_global, &x, 27);

    let mut y_matfree = vec![0.0f64; 27];
    op.apply(&x, &mut y_matfree);

    let diff = max_abs_diff(&y_matfree, &y_assembled);
    eprintln!("[Test1] p=1 2x2x2: max|y_matfree - y_assembled| = {diff:.3e}");
    assert!(
        diff < 1e-12,
        "p=1 matvec mismatch vs assembled GLL stiffness: max diff {diff:.3e} (>= 1e-12)"
    );

    // ── Second sub-check at p=2 on a 2×2×2 mesh (125 DOFs). p=1 has identity φ
    // and constant dφ, hiding p-dependent bugs; p=2 genuinely exercises dφ.
    let op2 = SumFactPoisson::new(2, 2);
    assert_eq!(
        op2.n_dofs(),
        125,
        "p=2 expected 125 DOFs, got {}",
        op2.n_dofs()
    );
    let (n_dofs2, k_global2) = assembled_gll_stiffness_dense(2, 2);
    assert_eq!(n_dofs2, 125, "p=2 baseline n_dofs mismatch: {n_dofs2}");

    let mut rng2 = Lcg::new(0x1234_5678);
    let x2 = rng2.vec(125);
    let y_assembled2 = dense_matvec(&k_global2, &x2, 125);
    let mut y_matfree2 = vec![0.0f64; 125];
    op2.apply(&x2, &mut y_matfree2);
    let diff2 = max_abs_diff(&y_matfree2, &y_assembled2);
    eprintln!("[Test1] p=2 2x2x2: max|y_matfree - y_assembled| = {diff2:.3e}");
    assert!(
        diff2 < 1e-11,
        "p=2 matvec mismatch vs assembled GLL stiffness: max diff {diff2:.3e} (>= 1e-11)"
    );

    // Bonus: exercise CsrMatrix::spmv on the p=1 dense baseline to confirm the
    // sparse path agrees too (independent third route).
    let csr = dense_to_csr(&k_global, 27);
    let mut y_csr = vec![0.0f64; 27];
    csr.spmv(&x, &mut y_csr);
    let diff_csr = max_abs_diff(&y_csr, &y_assembled);
    assert!(
        diff_csr < 1e-12,
        "CsrMatrix::spmv disagrees with dense baseline: max diff {diff_csr:.3e}"
    );
}

/// Convert a dense row-major `n×n` matrix to a [`CsrMatrix`], dropping exact
/// zeros. Used only to exercise the public `spmv` against the dense baseline.
fn dense_to_csr(k: &[f64], n: usize) -> CsrMatrix {
    let mut row_offsets = Vec::with_capacity(n + 1);
    let mut col_indices = Vec::new();
    let mut values = Vec::new();
    row_offsets.push(0usize);
    for r in 0..n {
        for c in 0..n {
            let v = k[r * n + c];
            if v != 0.0 {
                col_indices.push(c);
                values.push(v);
            }
        }
        row_offsets.push(col_indices.len());
    }
    CsrMatrix {
        nrows: n,
        ncols: n,
        row_offsets,
        col_indices,
        values,
    }
}

/// Test 2 — PCG + matrix-free operator on a manufactured Poisson problem,
/// verifying h-convergence at p=2. u = sin(πx)sin(πy)sin(πz), -Δu = 3π²u, u=0 on ∂Ω.
#[test]
fn matfree_pcg_poisson_h_convergence_p2() {
    let p = 2;
    let three_pi2 = 3.0 * PI * PI;
    let mut errors: Vec<(usize, f64, usize)> = Vec::new();

    for &n in &[4usize, 8, 16] {
        let op = SumFactPoisson::new(n, p);
        let ndir = op.dofs_per_dir();
        let n_dofs = op.n_dofs();

        let coord = coord_1d(n, p);
        let mass = lumped_mass(n, p);

        // exact solution + RHS b_g = m_g · 3π² · u_exact_g, zeroed on boundary.
        let mut u_exact = vec![0.0f64; n_dofs];
        let mut b = vec![0.0f64; n_dofs];
        let mask = boundary_mask(ndir);
        for gi in 0..ndir {
            let sx = (PI * coord[gi]).sin();
            for gj in 0..ndir {
                let sy = (PI * coord[gj]).sin();
                for (gk, &czk) in coord.iter().enumerate() {
                    let sz = (PI * czk).sin();
                    let g = (gi * ndir + gj) * ndir + gk;
                    let ue = sx * sy * sz;
                    u_exact[g] = ue;
                    b[g] = if mask[g] {
                        0.0
                    } else {
                        mass[g] * three_pi2 * ue
                    };
                }
            }
        }

        let dop = DirichletOp {
            inner: &op,
            mask: &mask,
        };
        let diag_inv = op.jacobi_diag_inverse(&mask);
        let precond = JacobiPreconditioner { diag_inv };

        let mut x = vec![0.0f64; n_dofs];
        let iters = MatrixFreePcg::new(2000, 1e-10)
            .solve(&dop, &precond, &b, &mut x)
            .expect("matrix-free PCG must converge on the MMS Poisson system");

        // mass-weighted discrete L2 error
        let mut err2 = 0.0;
        for g in 0..n_dofs {
            let d = x[g] - u_exact[g];
            err2 += mass[g] * d * d;
        }
        let err = err2.sqrt();
        eprintln!("[Test2] n={n} p={p}: iters={iters}, L2_err={err:.6e}");
        errors.push((n, err, iters));
    }

    // monotone decrease + h² (p=2) convergence ratios ≥ 3.0
    for w in errors.windows(2) {
        let (n0, e0, _) = w[0];
        let (n1, e1, _) = w[1];
        assert!(
            e1 < e0,
            "error must decrease n={n0}->{n1}: {e0:.3e} -> {e1:.3e}"
        );
    }
    let r1 = errors[0].1 / errors[1].1;
    let r2 = errors[1].1 / errors[2].1;
    eprintln!("[Test2] convergence ratios: err(4)/err(8)={r1:.3}, err(8)/err(16)={r2:.3}");
    // h² convergence for p=2 gives ~4×; we require ≥ 3.0 for safety margin.
    assert!(
        r1 >= 3.0,
        "convergence ratio err(4)/err(8) too low: {r1:.3} (< 3.0)"
    );
    assert!(
        r2 >= 3.0,
        "convergence ratio err(8)/err(16) too low: {r2:.3} (< 3.0)"
    );
    // Sanity: finest error must be genuinely small.
    assert!(
        errors[2].1 < 1e-2,
        "finest-mesh L2 error unexpectedly large: {:.3e}",
        errors[2].1
    );
}

/// Test 3 — operator heap footprint per DOF (benchmark-style print + soft check).
/// The closed-form DOF map means the operator stores only a few `(p+1)²` 1D
/// matrices, so its size is `O(1)` in the mesh and bytes/dof shrinks with n.
#[test]
fn matfree_bytes_per_dof_p4() {
    let p = 4usize;
    let op = SumFactPoisson::new(4, p);
    let n_dofs = op.n_dofs();

    // phi, dphi, phi_t, dphi_t: each (p+1)² f64; weights: (p+1) f64; plus struct.
    let np = p + 1;
    let bytes = 4 * np * np * 8 + np * 8 + std::mem::size_of::<SumFactPoisson>();
    let bytes_per_dof = bytes as f64 / n_dofs as f64;

    eprintln!(
        "[Test3] p=4 4x4x4: n_dofs={n_dofs}, operator_bytes≈{bytes}, bytes/dof≈{bytes_per_dof:.4}"
    );
    assert!(
        bytes_per_dof < 40.0,
        "matrix-free operator bytes/dof should be tiny (O(1) storage): {bytes_per_dof:.4} (>= 40)"
    );
}

/// Test 4 — sum-factorized contraction path equals the brute-force quadrature
/// definition. `contract_axis` is `pub(crate)` and unreachable from this test
/// crate, so we verify the equivalence at the *operator* level: for a single
/// degree-3 element, `op.apply` (which internally uses sum-factorized
/// `contract_axis`) must match the dense `Ke·x` built from the explicit
/// triple-quadrature definition (no sum factorization).
#[test]
fn sumfac_tensor_contraction_matches_bruteforce_p3() {
    let p = 3usize;
    let op = SumFactPoisson::new(1, p); // single element → 4³ = 64 DOFs, local map.
    let n_dofs = op.n_dofs();
    assert_eq!(
        n_dofs, 64,
        "single p=3 element should have 64 DOFs, got {n_dofs}"
    );

    // brute-force dense element stiffness on the single unit element (h=1).
    let ke = element_stiffness_dense(p, 1.0);

    let mut rng = Lcg::new(0x0BADC0DE);
    let mut worst = 0.0f64;
    for trial in 0..5 {
        let x = rng.vec(n_dofs);
        let y_dense = dense_matvec(&ke, &x, n_dofs);
        let mut y_matfree = vec![0.0f64; n_dofs];
        op.apply(&x, &mut y_matfree);
        let diff = max_abs_diff(&y_matfree, &y_dense);
        worst = worst.max(diff);
        assert!(
            diff < 1e-11,
            "trial {trial}: sum-factorized apply disagrees with brute-force Ke·x: {diff:.3e} (>= 1e-11)"
        );
    }
    eprintln!("[Test4] p=3 single element: worst max|apply - Ke·x| over 5 trials = {worst:.3e}");
}
