// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Analytic partial derivatives of forward dynamics (ABA) with respect to the
//! joint state.
//!
//! Given forward dynamics `q̈ = ABA(q, q̇, τ)` and the joint-space equation of
//! motion `M(q)·q̈ + C(q, q̇, q̇) = τ` (where `RNEA(q, q̇, q̈) = M·q̈ + drift`),
//! the sensitivities of the joint accelerations follow from differentiating the
//! implicit relation `RNEA(q, q̇, q̈) = τ` (Carpentier & Mansard 2018):
//!
//! ```text
//! ∂q̈/∂τ  =  M(q)^{-1}
//! ∂q̈/∂q  = −M(q)^{-1} · ∂RNEA/∂q   evaluated at (q, q̇, q̈ = ABA(q, q̇, τ))
//! ∂q̈/∂q̇ = −M(q)^{-1} · ∂RNEA/∂q̇  evaluated at the same (q, q̇, q̈)
//! ```
//!
//! The crucial subtlety is that the RNEA derivatives must be evaluated at the
//! acceleration `q̈ = ABA(q, q̇, τ)` produced by the operating point, because
//! `∂RNEA/∂q` itself depends on `q̈`. This module therefore first runs ABA, then
//! evaluates [`crate::rnea_derivatives::rnea_derivatives`] at that `q̈`, and
//! finally reuses a single Cholesky factorisation of `M(q)` for all the linear
//! solves.
//!
//! Reference: Justin Carpentier and Nicolas Mansard, "Analytical Derivatives of
//! Rigid Body Dynamics Algorithms", Robotics: Science and Systems (RSS), 2018.

use crate::{
    aba::aba, crba::compute_mass_matrix_crba, model::ArticulatedModel,
    rnea_derivatives::rnea_derivatives,
};

/// Tolerance below which a Cholesky pivot is treated as non-positive-definite.
const CHOLESKY_EPS: f64 = 1e-12;

/// Error returned when the mass matrix cannot be Cholesky-factorised
/// (not symmetric-positive-definite within tolerance).
#[derive(Debug, Clone, PartialEq)]
pub enum AbaDerivativeError {
    /// The mass matrix was not positive-definite at the given pivot index.
    NotPositiveDefinite {
        /// Zero-based pivot index where factorisation failed.
        index: usize,
    },
    /// Input slice length did not match `model.total_dof()`.
    DimensionMismatch {
        /// Expected length (`model.total_dof()`).
        expected: usize,
        /// Actual length supplied.
        got: usize,
    },
}

impl std::fmt::Display for AbaDerivativeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NotPositiveDefinite { index } => write!(
                f,
                "mass matrix is not positive-definite at pivot index {index}"
            ),
            Self::DimensionMismatch { expected, got } => {
                write!(f, "dimension mismatch: expected {expected}, got {got}")
            }
        }
    }
}

impl std::error::Error for AbaDerivativeError {}

/// Partial derivatives of forward dynamics q̈ = ABA(q, q̇, τ).
///
/// All three are `n×n` row-major (`n = model.total_dof()`):
/// `dqdd_dq[i][j] = ∂q̈_i/∂q_j`, `dqdd_dqd[i][j] = ∂q̈_i/∂q̇_j`,
/// `dqdd_dtau[i][j] = ∂q̈_i/∂τ_j` (the latter equals `M(q)^{-1}`).
#[derive(Debug, Clone)]
pub struct AbaDerivatives {
    /// ∂q̈/∂q  (n×n).
    pub dqdd_dq: Vec<Vec<f64>>,
    /// ∂q̈/∂q̇ (n×n).
    pub dqdd_dqd: Vec<Vec<f64>>,
    /// ∂q̈/∂τ = M(q)^{-1} (n×n, symmetric).
    pub dqdd_dtau: Vec<Vec<f64>>,
}

/// In-place lower-triangular Cholesky factor `L` of an SPD matrix `a` (n×n
/// row-major `Vec<Vec<f64>>`), such that `a = L·Lᵀ`. Returns
/// `Err(NotPositiveDefinite)` if a pivot is ≤ `eps`.
///
/// Uses the standard Cholesky–Banachiewicz recursion; the upper-triangular
/// entries of the returned factor are zero.
fn cholesky_factor(a: &[Vec<f64>], eps: f64) -> Result<Vec<Vec<f64>>, AbaDerivativeError> {
    let n = a.len();
    let mut l = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a[i][j];
            for (k, l_jk) in l[j].iter().enumerate().take(j) {
                sum -= l[i][k] * l_jk;
            }
            if i == j {
                if sum <= eps {
                    return Err(AbaDerivativeError::NotPositiveDefinite { index: i });
                }
                l[i][i] = sum.sqrt();
            } else {
                l[i][j] = sum / l[j][j];
            }
        }
    }
    Ok(l)
}

/// Solve `M·x = b` for a single right-hand side using the precomputed Cholesky
/// factor `l` (lower-triangular). Performs a forward substitution `L·y = b`
/// followed by a back substitution `Lᵀ·x = y`. `b` length = `n`.
fn cholesky_solve_vec(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    // Forward solve: L·y = b.
    let mut y = vec![0.0f64; n];
    for i in 0..n {
        let mut sum = b[i];
        for k in 0..i {
            sum -= l[i][k] * y[k];
        }
        y[i] = sum / l[i][i];
    }
    // Back solve: Lᵀ·x = y.
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

/// Compute the inverse `M(q)^{-1}` of the joint-space mass matrix via a single
/// Cholesky factorisation.
///
/// Returns an `n×n` symmetric matrix (`n = model.total_dof()`). Each column is
/// obtained by solving `M·x = e_j` for the `j`-th unit vector.
///
/// `model` is taken by `&mut` because [`compute_mass_matrix_crba`] requires
/// mutable access to the model.
///
/// # Errors
/// Returns [`AbaDerivativeError::DimensionMismatch`] if `q.len()` differs from
/// `model.total_dof()`, or [`AbaDerivativeError::NotPositiveDefinite`] if the
/// mass matrix is not SPD within tolerance.
pub fn mass_matrix_inverse(
    model: &mut ArticulatedModel,
    q: &[f64],
) -> Result<Vec<Vec<f64>>, AbaDerivativeError> {
    let n = model.total_dof();
    if q.len() != n {
        return Err(AbaDerivativeError::DimensionMismatch {
            expected: n,
            got: q.len(),
        });
    }
    let m = compute_mass_matrix_crba(model, q);
    let l = cholesky_factor(&m, CHOLESKY_EPS)?;
    Ok(invert_from_factor(&l, n))
}

/// Assemble `M^{-1}` column-by-column from a Cholesky factor `l` of an `n×n`
/// matrix, solving `M·x = e_j` for each unit column `j`.
fn invert_from_factor(l: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let mut inv = vec![vec![0.0f64; n]; n];
    let mut e = vec![0.0f64; n];
    for j in 0..n {
        e[j] = 1.0;
        let col = cholesky_solve_vec(l, &e);
        for (i, inv_row) in inv.iter_mut().enumerate() {
            inv_row[j] = col[i];
        }
        e[j] = 0.0;
    }
    inv
}

/// Compute ∂q̈/∂q, ∂q̈/∂q̇, and ∂q̈/∂τ for forward dynamics via the identity
/// ∂q̈/∂x = −M(q)^{-1}·∂RNEA/∂x (x∈{q,q̇}) and ∂q̈/∂τ = M(q)^{-1}.
///
/// Internally: q̈ = aba(q,q̇,τ); M = CRBA(q); ∂RNEA evaluated at (q,q̇,q̈);
/// then a Cholesky factorisation of M is reused for all 2n+ n solves.
///
/// `model` is taken by `&mut` because [`compute_mass_matrix_crba`] requires
/// mutable access to the model.
///
/// # Errors
/// Returns [`AbaDerivativeError`] if dimensions mismatch or `M` is not SPD.
pub fn aba_derivatives(
    model: &mut ArticulatedModel,
    q: &[f64],
    q_dot: &[f64],
    tau: &[f64],
) -> Result<AbaDerivatives, AbaDerivativeError> {
    let n = model.total_dof();
    if q.len() != n {
        return Err(AbaDerivativeError::DimensionMismatch {
            expected: n,
            got: q.len(),
        });
    }
    if q_dot.len() != n {
        return Err(AbaDerivativeError::DimensionMismatch {
            expected: n,
            got: q_dot.len(),
        });
    }
    if tau.len() != n {
        return Err(AbaDerivativeError::DimensionMismatch {
            expected: n,
            got: tau.len(),
        });
    }

    // Forward dynamics at the operating point.
    let q_ddot = aba(model, q, q_dot, tau);

    // Joint-space mass matrix and its single Cholesky factorisation.
    let m = compute_mass_matrix_crba(model, q);
    let l = cholesky_factor(&m, CHOLESKY_EPS)?;

    // RNEA derivatives evaluated AT the operating-point acceleration q̈.
    let derivs = rnea_derivatives(model, q, q_dot, &q_ddot);

    // ∂q̈/∂τ = M^{-1}, assembled column-by-column.
    let dqdd_dtau = invert_from_factor(&l, n);

    // ∂q̈/∂q = −M^{-1}·∂RNEA/∂q and ∂q̈/∂q̇ = −M^{-1}·∂RNEA/∂q̇.
    // Each column j is a single Cholesky solve against the j-th column of the
    // corresponding RNEA-derivative matrix.
    let dqdd_dq = neg_minv_times(&l, &derivs.dtau_dq, n);
    let dqdd_dqd = neg_minv_times(&l, &derivs.dtau_dqd, n);

    Ok(AbaDerivatives {
        dqdd_dq,
        dqdd_dqd,
        dqdd_dtau,
    })
}

/// Compute `−M^{-1}·B` column-by-column, given the Cholesky factor `l` of `M`
/// and an `n×n` right-hand-side matrix `b` (row-major). For each column `j` the
/// right-hand side is the `j`-th column of `b`; the solved column is negated and
/// written back as column `j` of the result.
fn neg_minv_times(l: &[Vec<f64>], b: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0f64; n]; n];
    let mut rhs = vec![0.0f64; n];
    for j in 0..n {
        for (i, rhs_i) in rhs.iter_mut().enumerate() {
            *rhs_i = b[i][j];
        }
        let col = cholesky_solve_vec(l, &rhs);
        for (i, out_row) in out.iter_mut().enumerate() {
            out_row[j] = -col[i];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        body::RigidBody,
        joint::RevoluteJoint,
        model::ArticulatedModel,
        spatial::{SpatialInertia, SpatialTransform},
    };

    // ── Numerical helpers ─────────────────────────────────────────────────────

    fn assert_close(a: f64, b: f64, eps: f64, label: &str) {
        assert!(
            (a - b).abs() < eps,
            "{label}: expected {b:.12}, got {a:.12} (|Δ|={:.3e})",
            (a - b).abs()
        );
    }

    /// Multiply an `n×n` row-major matrix by an `n`-vector.
    fn matvec(m: &[Vec<f64>], x: &[f64]) -> Vec<f64> {
        m.iter()
            .map(|row| row.iter().zip(x.iter()).map(|(a, b)| a * b).sum())
            .collect()
    }

    /// Compute `−A·B` for `n×n` row-major matrices.
    fn neg_matmul(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = a.len();
        let mut out = vec![vec![0.0f64; n]; n];
        for i in 0..n {
            for j in 0..n {
                let mut s = 0.0;
                for k in 0..n {
                    s += a[i][k] * b[k][j];
                }
                out[i][j] = -s;
            }
        }
        out
    }

    // ── Model builders (mirroring rnea_derivatives.rs) ─────────────────────────

    fn build_pendulum(m: f64, l: f64, g: f64) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -g, 0.0]);
        let joint = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let com = [0.0, -l, 0.0];
        let inertia = SpatialInertia::from_com(m, com, [[0.0; 3]; 3]);
        let body = RigidBody::new("link1", inertia, None, SpatialTransform::IDENTITY);
        model.add_body(body, joint);
        model
    }

    fn build_two_link(m1: f64, m2: f64, l1: f64, l2: f64, g: f64) -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -g, 0.0]);

        let j1 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let i1 = SpatialInertia::from_com(m1, [l1 * 0.5, 0.0, 0.0], [[0.0; 3]; 3]);
        let b1 = RigidBody::new("link1", i1, None, SpatialTransform::IDENTITY);
        model.add_body(b1, j1);

        let j2 = Box::new(RevoluteJoint::new([0.0, 0.0, 1.0]));
        let i2 = SpatialInertia::from_com(m2, [l2 * 0.5, 0.0, 0.0], [[0.0; 3]; 3]);
        let b2 = RigidBody::new(
            "link2",
            i2,
            Some(0),
            SpatialTransform::from_translation([l1, 0.0, 0.0]),
        );
        model.add_body(b2, j2);

        model
    }

    fn build_seven_dof() -> ArticulatedModel {
        let mut model = ArticulatedModel::new([0.0, -9.81, 0.0]);
        let axes = [
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let rot_com = [[0.01, 0.0, 0.0], [0.0, 0.01, 0.0], [0.0, 0.0, 0.01]];
        for (idx, axis) in axes.iter().enumerate() {
            let joint = Box::new(RevoluteJoint::new(*axis));
            let mass = 1.0 + 0.1 * idx as f64;
            let inertia = SpatialInertia::from_com(mass, [0.15, 0.0, 0.0], rot_com);
            let (parent_id, parent_transform) = if idx == 0 {
                (None, SpatialTransform::IDENTITY)
            } else {
                (
                    Some(idx - 1),
                    SpatialTransform::from_translation([0.3, 0.0, 0.0]),
                )
            };
            let body = RigidBody::new(format!("link{idx}"), inertia, parent_id, parent_transform);
            model.add_body(body, joint);
        }
        model
    }

    /// Central finite-difference of `aba` with respect to `tau`, step `h`.
    fn fd_aba_dtau(
        model: &ArticulatedModel,
        q: &[f64],
        q_dot: &[f64],
        tau: &[f64],
        h: f64,
    ) -> Vec<Vec<f64>> {
        let n = model.total_dof();
        let mut out = vec![vec![0.0f64; n]; n];
        let inv_2h = 1.0 / (2.0 * h);
        let mut t = tau.to_vec();
        for j in 0..n {
            let orig = t[j];
            t[j] = orig + h;
            let plus = aba(model, q, q_dot, &t);
            t[j] = orig - h;
            let minus = aba(model, q, q_dot, &t);
            t[j] = orig;
            for i in 0..n {
                out[i][j] = (plus[i] - minus[i]) * inv_2h;
            }
        }
        out
    }

    /// Central finite-difference of `aba` with respect to `q`, step `h`.
    fn fd_aba_dq(
        model: &ArticulatedModel,
        q: &[f64],
        q_dot: &[f64],
        tau: &[f64],
        h: f64,
    ) -> Vec<Vec<f64>> {
        let n = model.total_dof();
        let mut out = vec![vec![0.0f64; n]; n];
        let inv_2h = 1.0 / (2.0 * h);
        let mut x = q.to_vec();
        for j in 0..n {
            let orig = x[j];
            x[j] = orig + h;
            let plus = aba(model, &x, q_dot, tau);
            x[j] = orig - h;
            let minus = aba(model, &x, q_dot, tau);
            x[j] = orig;
            for i in 0..n {
                out[i][j] = (plus[i] - minus[i]) * inv_2h;
            }
        }
        out
    }

    /// Central finite-difference of `aba` with respect to `q̇`, step `h`.
    fn fd_aba_dqd(
        model: &ArticulatedModel,
        q: &[f64],
        q_dot: &[f64],
        tau: &[f64],
        h: f64,
    ) -> Vec<Vec<f64>> {
        let n = model.total_dof();
        let mut out = vec![vec![0.0f64; n]; n];
        let inv_2h = 1.0 / (2.0 * h);
        let mut x = q_dot.to_vec();
        for j in 0..n {
            let orig = x[j];
            x[j] = orig + h;
            let plus = aba(model, q, &x, tau);
            x[j] = orig - h;
            let minus = aba(model, q, &x, tau);
            x[j] = orig;
            for i in 0..n {
                out[i][j] = (plus[i] - minus[i]) * inv_2h;
            }
        }
        out
    }

    /// Maximum absolute entry-wise difference between two `n×n` matrices.
    fn max_abs_diff(a: &[Vec<f64>], b: &[Vec<f64>]) -> f64 {
        let mut worst = 0.0f64;
        for (ar, br) in a.iter().zip(b.iter()) {
            for (av, bv) in ar.iter().zip(br.iter()) {
                worst = worst.max((av - bv).abs());
            }
        }
        worst
    }

    // ── Cholesky unit tests ───────────────────────────────────────────────────

    #[test]
    fn test_cholesky_identity() {
        let a = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let l = cholesky_factor(&a, CHOLESKY_EPS).expect("identity is SPD");
        assert_close(l[0][0], 1.0, 1e-12, "L[0][0]");
        assert_close(l[0][1], 0.0, 1e-12, "L[0][1]");
        assert_close(l[1][0], 0.0, 1e-12, "L[1][0]");
        assert_close(l[1][1], 1.0, 1e-12, "L[1][1]");

        let b = [3.0, -5.0];
        let x = cholesky_solve_vec(&l, &b);
        assert_close(x[0], 3.0, 1e-12, "x[0]");
        assert_close(x[1], -5.0, 1e-12, "x[1]");
    }

    #[test]
    fn test_cholesky_spd_solve() {
        let a = vec![vec![4.0, 1.0], vec![1.0, 3.0]];
        let l = cholesky_factor(&a, CHOLESKY_EPS).expect("matrix is SPD");
        let b = [1.0, 2.0];
        let x = cholesky_solve_vec(&l, &b);
        // Multiply back: M·x should recover b.
        let recovered = matvec(&a, &x);
        assert_close(recovered[0], b[0], 1e-10, "M·x[0]");
        assert_close(recovered[1], b[1], 1e-10, "M·x[1]");
    }

    #[test]
    fn test_cholesky_rejects_non_spd() {
        // Indefinite matrix [[1,2],[2,1]] (eigenvalues 3 and -1).
        let a = vec![vec![1.0, 2.0], vec![2.0, 1.0]];
        let result = cholesky_factor(&a, CHOLESKY_EPS);
        assert!(
            matches!(result, Err(AbaDerivativeError::NotPositiveDefinite { .. })),
            "expected NotPositiveDefinite, got {result:?}"
        );
    }

    // ── ∂q̈/∂τ tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_dqdd_dtau_equals_minv_pendulum() {
        let mut model = build_pendulum(1.5, 0.8, 9.81);
        let q = [0.6];
        let qd = [0.4];
        let tau = [0.2];

        let m = compute_mass_matrix_crba(&mut model, &q);
        let derivs = aba_derivatives(&mut model, &q, &qd, &tau).expect("SPD pendulum");
        // 1×1: ∂q̈/∂τ == 1/M[0][0].
        assert_close(derivs.dqdd_dtau[0][0], 1.0 / m[0][0], 1e-10, "1/M[0][0]");
    }

    #[test]
    fn test_dqdd_dtau_matches_fd_two_link() {
        let mut model = build_two_link(1.2, 0.9, 0.5, 0.4, 9.81);
        let q = [0.3, -0.5];
        let qd = [0.2, 0.7];
        let tau = [0.15, -0.25];

        let derivs = aba_derivatives(&mut model, &q, &qd, &tau).expect("SPD two-link");
        let fd = fd_aba_dtau(&model, &q, &qd, &tau, 1e-6);
        let worst = max_abs_diff(&derivs.dqdd_dtau, &fd);
        assert!(worst < 1e-6, "∂q̈/∂τ vs FD: worst |Δ|={worst:.3e}");
    }

    #[test]
    fn test_dqdd_dtau_symmetry_two_link() {
        let mut model = build_two_link(1.2, 0.9, 0.5, 0.4, 9.81);
        let q = [0.3, -0.5];
        let qd = [0.2, 0.7];
        let tau = [0.15, -0.25];

        let derivs = aba_derivatives(&mut model, &q, &qd, &tau).expect("SPD two-link");
        // M^{-1} is symmetric.
        assert_close(
            derivs.dqdd_dtau[0][1],
            derivs.dqdd_dtau[1][0],
            1e-10,
            "M^{-1} symmetry",
        );
    }

    // ── ∂q̈/∂q tests ───────────────────────────────────────────────────────────

    #[test]
    fn test_dqdd_dq_equals_minv_dtaudq_two_link() {
        let mut model = build_two_link(1.2, 0.9, 0.5, 0.4, 9.81);
        let q = [0.3, -0.5];
        let qd = [0.2, 0.7];
        let tau = [0.15, -0.25];

        let derivs = aba_derivatives(&mut model, &q, &qd, &tau).expect("SPD two-link");

        // Independently: q̈ = aba; M^{-1}; ∂τ/∂q at q̈; then −M^{-1}·∂τ/∂q.
        let q_ddot = aba(&model, &q, &qd, &tau);
        let minv = mass_matrix_inverse(&mut model, &q).expect("SPD inverse");
        let rnea_d = rnea_derivatives(&model, &q, &qd, &q_ddot);
        let expected = neg_matmul(&minv, &rnea_d.dtau_dq);

        let worst = max_abs_diff(&derivs.dqdd_dq, &expected);
        assert!(worst < 1e-8, "∂q̈/∂q vs −M⁻¹·∂τ/∂q: worst |Δ|={worst:.3e}");
    }

    #[test]
    fn test_dqdd_dq_matches_fd_two_link() {
        let mut model = build_two_link(1.2, 0.9, 0.5, 0.4, 9.81);
        let q = [0.3, -0.5];
        let qd = [0.2, 0.7];
        let tau = [0.15, -0.25];

        let derivs = aba_derivatives(&mut model, &q, &qd, &tau).expect("SPD two-link");

        // Try 1e-6 first; the FD of forward dynamics passes at this tolerance.
        let fd_q = fd_aba_dq(&model, &q, &qd, &tau, 1e-6);
        let fd_qd = fd_aba_dqd(&model, &q, &qd, &tau, 1e-6);

        let worst_q = max_abs_diff(&derivs.dqdd_dq, &fd_q);
        let worst_qd = max_abs_diff(&derivs.dqdd_dqd, &fd_qd);
        assert!(worst_q < 1e-6, "∂q̈/∂q vs FD: worst |Δ|={worst_q:.3e}");
        assert!(worst_qd < 1e-6, "∂q̈/∂q̇ vs FD: worst |Δ|={worst_qd:.3e}");
    }

    // ── Full-pipeline scale test ──────────────────────────────────────────────

    #[test]
    fn test_seven_dof_dqdd_dtau_vs_fd() {
        let mut model = build_seven_dof();
        let q = [0.2, -0.4, 0.6, -0.1, 0.3, -0.5, 0.15];
        let qd = [0.1, 0.2, -0.3, 0.4, -0.2, 0.1, 0.25];
        let tau = [0.05, -0.1, 0.2, -0.15, 0.1, -0.05, 0.12];

        let derivs = aba_derivatives(&mut model, &q, &qd, &tau).expect("SPD 7-DoF");
        let fd = fd_aba_dtau(&model, &q, &qd, &tau, 1e-6);
        let worst = max_abs_diff(&derivs.dqdd_dtau, &fd);
        assert!(worst < 1e-5, "7-DoF ∂q̈/∂τ vs FD: worst |Δ|={worst:.3e}");
    }
}
