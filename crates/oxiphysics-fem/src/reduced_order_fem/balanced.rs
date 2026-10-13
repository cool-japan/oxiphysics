// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Balanced truncation of stable state-space systems for the reduced-order FEM
//! toolkit: the projectors returned by the square-root method of
//! [`crate::reduced_order::balancing_transform`] packaged with the helpers that
//! apply them to full-order matrices and vectors.

use super::{mat_mul, mat_vec};
use crate::analysis::FemError;
use crate::reduced_order::balancing_transform;

/// Result of balanced truncation of a stable state-space system `(A, B, C)`.
///
/// `T_l` (`r × N`) and `T_r` (`N × r`) are the biorthogonal balancing
/// projectors of the square-root method, `T_l T_r = I_r`.  The reduced system
/// is `(T_l A T_r, T_l B, C T_r)`.
#[derive(Clone, Debug)]
pub struct BalancedTruncation {
    /// Hankel singular values `σ₁ ≥ … ≥ σ_N` of the **full** system (length
    /// `N`); the first `reduced_dim` of them belong to the retained states.
    pub hankel_singular_values: Vec<f64>,
    /// Left transformation matrix T_l (r × N).
    pub transform_left: Vec<Vec<f64>>,
    /// Right transformation matrix T_r (N × r).
    pub transform_right: Vec<Vec<f64>>,
    /// Reduced dimension.
    pub reduced_dim: usize,
}

impl BalancedTruncation {
    /// Project a full-order state matrix A → T_l A T_r  (r × r).
    pub fn project_matrix(&self, a: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if self.reduced_dim == 0 {
            return Vec::new();
        }
        mat_mul(&self.transform_left, &mat_mul(a, &self.transform_right))
    }

    /// Project a full-order input matrix B (N × m) → T_l B  (r × m).
    pub fn project_input(&self, b: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if self.reduced_dim == 0 {
            return Vec::new();
        }
        mat_mul(&self.transform_left, b)
    }

    /// Project a full-order output matrix C (p × N) → C T_r  (p × r).
    pub fn project_output(&self, c: &[Vec<f64>]) -> Vec<Vec<f64>> {
        if self.reduced_dim == 0 {
            return vec![Vec::new(); c.len()];
        }
        mat_mul(c, &self.transform_right)
    }

    /// Project a full-order vector b → T_l b.
    pub fn project_vector(&self, b: &[f64]) -> Vec<f64> {
        mat_vec(&self.transform_left, b)
    }

    /// Lift a reduced state z (length r) back to the full space: T_r z.
    pub fn lift_state(&self, z: &[f64]) -> Vec<f64> {
        mat_vec(&self.transform_right, z)
    }

    /// A-priori bound `2 (σ_{r+1} + … + σ_N)` on the `H∞` norm of the
    /// difference between the full and the reduced transfer function.
    pub fn error_bound(&self) -> f64 {
        2.0 * self
            .hankel_singular_values
            .iter()
            .skip(self.reduced_dim)
            .sum::<f64>()
    }
}

/// Balanced truncation of the stable linear system `ẋ = A x + B u`, `y = C x`
/// by the square-root method.
///
/// The controllability and observability Gramians are the solutions of the
/// Lyapunov equations `A Wc + Wc Aᵀ + B Bᵀ = 0` and `Aᵀ Wo + Wo A + Cᵀ C = 0`
/// (see [`crate::reduced_order::solve_lyapunov`]); their Cholesky-like factors
/// are combined through an SVD to obtain the Hankel singular values and the
/// biorthogonal projectors `T_l`, `T_r`.  This is the same numerical core as
/// [`crate::reduced_order::balanced_truncation`].
///
/// - `a`          – system matrix A (N × N), asymptotically stable.
/// - `b_mat`      – input matrix B (N × n_inputs).
/// - `c_mat`      – output matrix C (n_outputs × N).
/// - `max_modes`  – number of balanced modes to retain (fewer are kept if the
///   system has fewer states that are both controllable and observable).
///
/// # Errors
/// [`FemError::InvalidInput`] for inconsistent shapes or an unstable `A`,
/// [`FemError::SingularSystem`] if `A` has eigenvalues on the imaginary axis,
/// [`FemError::SolverNotConverged`] if the Lyapunov solver fails.
pub fn balanced_truncation(
    a: &[Vec<f64>],
    b_mat: &[Vec<f64>],
    c_mat: &[Vec<f64>],
    max_modes: usize,
) -> Result<BalancedTruncation, FemError> {
    let transform = balancing_transform(a, b_mat, c_mat, max_modes)?;
    let reduced_dim = transform.reduced_dim();
    Ok(BalancedTruncation {
        hankel_singular_values: transform.hankel_singular_values,
        transform_left: transform.t_left,
        transform_right: transform.t_right,
        reduced_dim,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// State-space triple `(A, B, C)`.
    type StateSpace = (Vec<Vec<f64>>, Vec<Vec<f64>>, Vec<Vec<f64>>);

    /// Stable test system: the diagonal pole pair -1, -2 with unit input/output.
    fn two_pole_system() -> StateSpace {
        (
            vec![vec![-1.0, 0.0], vec![0.0, -2.0]],
            vec![vec![1.0], vec![1.0]],
            vec![vec![1.0, 1.0]],
        )
    }

    #[test]
    fn test_balanced_truncation_hankel_singular_values_of_two_pole_system() {
        // Wc = Wo = [[1/2, 1/3], [1/3, 1/4]], so sigma = eig(Wc).
        let (a, b_mat, c_mat) = two_pole_system();
        let bt = balanced_truncation(&a, &b_mat, &c_mat, 2).expect("stable system");
        let half_trace: f64 = 0.375;
        let root = (half_trace * half_trace - 1.0 / 72.0).sqrt();
        assert_eq!(bt.hankel_singular_values.len(), 2);
        assert!((bt.hankel_singular_values[0] - (half_trace + root)).abs() < 1e-12);
        assert!((bt.hankel_singular_values[1] - (half_trace - root)).abs() < 1e-12);
        assert!((bt.hankel_singular_values[0] - 0.7310).abs() < 1e-4);
        assert!((bt.hankel_singular_values[1] - 0.0190).abs() < 1e-4);
        assert_eq!(bt.reduced_dim, 2);
    }

    #[test]
    fn test_balanced_truncation_returns_biorthogonal_projectors() {
        let (a, b_mat, c_mat) = two_pole_system();
        let bt = balanced_truncation(&a, &b_mat, &c_mat, 1).expect("stable system");
        assert_eq!(bt.reduced_dim, 1);
        assert_eq!(bt.transform_left.len(), 1);
        assert_eq!(bt.transform_left[0].len(), 2);
        assert_eq!(bt.transform_right.len(), 2);
        assert_eq!(bt.transform_right[0].len(), 1);
        // T_l T_r = I_1
        let product = mat_mul(&bt.transform_left, &bt.transform_right);
        assert!(
            (product[0][0] - 1.0).abs() < 1e-12,
            "T_l T_r = {}",
            product[0][0]
        );
        // Bound: 2 * sigma_2.
        assert!((bt.error_bound() - 2.0 * bt.hankel_singular_values[1]).abs() < 1e-14);
    }

    #[test]
    fn test_balanced_truncation_project_matrix_uses_the_right_transform() {
        // A_r = T_l A T_r must be a 1x1 matrix for r = 1 (the previous code
        // multiplied by T_lᵀ and returned a mis-shaped N x r product).
        let (a, b_mat, c_mat) = two_pole_system();
        let bt = balanced_truncation(&a, &b_mat, &c_mat, 1).expect("stable system");
        let a_r = bt.project_matrix(&a);
        assert_eq!(a_r.len(), 1);
        assert_eq!(a_r[0].len(), 1);
        let expected = mat_mul(&bt.transform_left, &mat_mul(&a, &bt.transform_right));
        assert!((a_r[0][0] - expected[0][0]).abs() < 1e-14);
        assert!(
            a_r[0][0] < 0.0,
            "reduced model must stay stable: {}",
            a_r[0][0]
        );
        let b_r = bt.project_input(&b_mat);
        let c_r = bt.project_output(&c_mat);
        assert_eq!((b_r.len(), b_r[0].len()), (1, 1));
        assert_eq!((c_r.len(), c_r[0].len()), (1, 1));
        // DC gain of the reduced model approximates 1/1 + 1/2 = 1.5 within the bound.
        let dc_full = 1.0 + 0.5;
        let dc_reduced = -c_r[0][0] * b_r[0][0] / a_r[0][0];
        assert!(
            (dc_full - dc_reduced).abs() <= bt.error_bound() + 1e-12,
            "DC gain {dc_reduced} vs {dc_full} exceeds bound {}",
            bt.error_bound()
        );
    }

    #[test]
    fn test_balanced_truncation_project_vector() {
        let (a, b_mat, c_mat) = two_pole_system();
        let bt = balanced_truncation(&a, &b_mat, &c_mat, 1).expect("stable system");
        let proj = bt.project_vector(&[1.0, 0.0]);
        assert_eq!(proj.len(), bt.reduced_dim);
        let lifted = bt.lift_state(&proj);
        assert_eq!(lifted.len(), 2);
        // T_r (T_l v) is the oblique projection of v onto the retained subspace:
        // projecting again must not change it.
        let again = bt.lift_state(&bt.project_vector(&lifted));
        for (x, y) in lifted.iter().zip(again.iter()) {
            assert!((x - y).abs() < 1e-12);
        }
    }

    #[test]
    fn test_balanced_truncation_rejects_unstable_systems() {
        // The previous implementation happily "reduced" the identity matrix.
        let a = vec![
            vec![1.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0],
            vec![0.0, 0.0, 1.0, 0.0],
            vec![0.0, 0.0, 0.0, 1.0],
        ];
        let b_mat = vec![vec![1.0], vec![0.5], vec![0.0], vec![0.0]];
        let c_mat = vec![vec![1.0, 0.0, 0.5, 0.0]];
        let result = balanced_truncation(&a, &b_mat, &c_mat, 2);
        assert!(
            matches!(result, Err(FemError::InvalidInput(_))),
            "unstable system must be rejected: {result:?}"
        );
    }

    #[test]
    fn test_balanced_truncation_marginally_stable_system_is_singular() {
        let a = vec![vec![0.0, 0.0], vec![0.0, -1.0]];
        let b_mat = vec![vec![1.0], vec![1.0]];
        let c_mat = vec![vec![1.0, 1.0]];
        assert_eq!(
            balanced_truncation(&a, &b_mat, &c_mat, 1).err(),
            Some(FemError::SingularSystem)
        );
    }

    #[test]
    fn test_balanced_truncation_rejects_inconsistent_shapes() {
        let (a, b_mat, c_mat) = two_pole_system();
        let short_b = vec![vec![1.0]];
        assert!(matches!(
            balanced_truncation(&a, &short_b, &c_mat, 1),
            Err(FemError::InvalidInput(_))
        ));
        let wide_c = vec![vec![1.0, 1.0, 1.0]];
        assert!(matches!(
            balanced_truncation(&a, &b_mat, &wide_c, 1),
            Err(FemError::InvalidInput(_))
        ));
    }
}
