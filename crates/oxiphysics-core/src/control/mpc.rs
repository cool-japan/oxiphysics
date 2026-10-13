// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Linear model predictive control solved as a condensed quadratic program.
//!
//! For the plant `x[k+1] = A x[k] + B u[k]` the controller minimises over the
//! input sequence `U = (u_0, ..., u_{N-1})`
//!
//! ```text
//! J(U) = sum_{k=0}^{N-1} ( q |x_{k+1} - r_{k+1}|^2 + r |u_k|^2 ),   x_0 = x,
//! ```
//!
//! optionally subject to `|u_k[i]| <= u_max`.  Stacking the predicted states
//! `X = (x_1, ..., x_N)` gives `X = Phi x + Gamma U` with
//! `Phi = [A; A^2; ...; A^N]` and the block lower-triangular
//! `Gamma[i][j] = A^(i-j) B` (`i >= j`), so that
//!
//! ```text
//! J(U) = U' H U + 2 g' U + const,   H = q Gamma' Gamma + r I,
//!                                   g = q Gamma' (Phi x - X_ref).
//! ```
//!
//! The Hessian `H` is positive semi-definite (positive definite for `r > 0`).
//! Without input bounds the minimiser solves `H U = -g` (Cholesky, with an
//! eigendecomposition pseudo-inverse if `H` is singular).  With bounds the
//! box-constrained QP is solved by an accelerated projected-gradient method
//! with step `1 / lambda_max(H)`, adaptive momentum restart and a
//! warm start at the clamped unconstrained solution.

use nalgebra::{DMatrix, DVector, SymmetricEigen};

/// Maximum number of projected-gradient iterations for the bounded problem.
const MAX_PROJECTED_GRADIENT_ITERATIONS: usize = 50_000;

/// Convergence tolerance on the max-norm of the iterate change.
const PROJECTED_GRADIENT_TOLERANCE: f64 = 1e-10;

/// Relative eigenvalue cut-off of the pseudo-inverse used for a singular Hessian.
const PSEUDO_INVERSE_RELATIVE_CUTOFF: f64 = 1e-12;

/// Linear MPC with prediction horizon `N`, tracking and input penalties and an
/// optional elementwise input bound.
///
/// At every step the controller minimises, over the inputs `u_0..u_{N-1}` of the
/// plant `x[k+1] = A x[k] + B u[k]`,
///
/// ```text
/// J = sum_{k=0}^{N-1} ( q_weight |x_{k+1} - ref_{k+1}|^2 + r_weight |u_k|^2 )
/// ```
///
/// subject to `|u_k[i]| <= u_max`, and returns the first input
/// ([`compute`](Self::compute)) or the whole sequence
/// ([`optimal_sequence`](Self::optimal_sequence)).  The problem is condensed
/// into a quadratic program in the stacked inputs and solved exactly:
/// Cholesky for the unconstrained case, an accelerated projected-gradient
/// method for the box-constrained case.
#[derive(Debug, Clone)]
pub struct ModelPredictiveControl {
    /// System A matrix (n×n).
    pub a: Vec<Vec<f64>>,
    /// System B matrix (n×m).
    pub b: Vec<Vec<f64>>,
    /// Prediction horizon.
    pub horizon: usize,
    /// State tracking weight (diagonal Q).
    pub q_weight: f64,
    /// Input weight (diagonal R).
    pub r_weight: f64,
    /// Input constraint: max magnitude of every input component
    /// (`f64::MAX` or larger means unconstrained).
    pub u_max: f64,
}

/// The condensed QP `min U'HU + 2 g'U` of one control step.
struct CondensedProblem {
    /// Hessian `H = q Gamma' Gamma + r I`.
    hessian: DMatrix<f64>,
    /// Linear term `g = q Gamma' (Phi x - X_ref)`.
    gradient: DVector<f64>,
    /// Number of inputs per step.
    input_dim: usize,
}

impl CondensedProblem {
    /// Minimiser of the unconstrained problem (minimum-norm if `H` is singular).
    fn solve_unconstrained(&self) -> Option<DVector<f64>> {
        let rhs = -&self.gradient;
        if let Some(cholesky) = self.hessian.clone().cholesky() {
            let solution = cholesky.solve(&rhs);
            if solution.iter().all(|v| v.is_finite()) {
                return Some(solution);
            }
        }
        // H is only positive semi-definite (r = 0 and a rank-deficient Gamma):
        // the system H U = -g is still consistent, take the minimum-norm solution.
        let eigen = SymmetricEigen::try_new(self.hessian.clone(), f64::EPSILON, 0)?;
        let largest = eigen
            .eigenvalues
            .iter()
            .fold(0.0_f64, |acc, &lambda| acc.max(lambda.abs()));
        let cutoff = PSEUDO_INVERSE_RELATIVE_CUTOFF * largest;
        let coefficients = eigen.eigenvectors.transpose() * &rhs;
        let scaled = DVector::from_fn(coefficients.len(), |i, _| {
            let lambda = eigen.eigenvalues[i];
            if lambda > cutoff {
                coefficients[i] / lambda
            } else {
                0.0
            }
        });
        let solution = &eigen.eigenvectors * scaled;
        solution.iter().all(|v| v.is_finite()).then_some(solution)
    }

    /// Upper bound of the largest eigenvalue of `H` (exact when the
    /// eigendecomposition converges, Gershgorin otherwise).
    fn lipschitz_constant(&self) -> f64 {
        let exact = SymmetricEigen::try_new(self.hessian.clone(), f64::EPSILON, 0)
            .map(|eigen| eigen.eigenvalues.iter().fold(0.0_f64, |acc, &l| acc.max(l)));
        match exact {
            Some(value) if value.is_finite() && value > 0.0 => value,
            _ => (0..self.hessian.nrows())
                .map(|i| self.hessian.row(i).iter().map(|v| v.abs()).sum::<f64>())
                .fold(0.0_f64, f64::max),
        }
    }

    /// Minimiser of the problem restricted to the box `|u_i| <= bound`.
    ///
    /// Accelerated projected gradient with step `1 / lambda_max(H)` and gradient
    /// restart (O'Donoghue and Candes), started from the clamped
    /// unconstrained solution.
    fn solve_box_constrained(&self, unconstrained: DVector<f64>, bound: f64) -> DVector<f64> {
        let clamp = |v: f64| v.clamp(-bound, bound);
        let mut current = unconstrained.map(clamp);
        let lipschitz = self.lipschitz_constant();
        if lipschitz <= 0.0 || !lipschitz.is_finite() {
            return current;
        }
        let step = 1.0 / lipschitz;
        let mut extrapolated = current.clone();
        let mut momentum = 1.0_f64;
        for _ in 0..MAX_PROJECTED_GRADIENT_ITERATIONS {
            let gradient = &self.hessian * &extrapolated + &self.gradient;
            let next = (&extrapolated - gradient * step).map(clamp);
            let from_extrapolated = max_abs(&(&next - &extrapolated));
            let from_current = max_abs(&(&next - &current));
            let restart = (&extrapolated - &next).dot(&(&next - &current)) > 0.0;
            if restart {
                momentum = 1.0;
                extrapolated = next.clone();
            } else {
                let next_momentum = 0.5 * (1.0 + (1.0 + 4.0 * momentum * momentum).sqrt());
                let beta = (momentum - 1.0) / next_momentum;
                extrapolated = &next + (&next - &current) * beta;
                momentum = next_momentum;
            }
            current = next;
            if from_extrapolated < PROJECTED_GRADIENT_TOLERANCE
                && from_current < PROJECTED_GRADIENT_TOLERANCE
            {
                break;
            }
        }
        current
    }
}

/// Largest absolute component of a vector.
fn max_abs(v: &DVector<f64>) -> f64 {
    v.iter().fold(0.0_f64, |acc, x| acc.max(x.abs()))
}

impl ModelPredictiveControl {
    /// Create a new MPC controller.
    pub fn new(
        a: Vec<Vec<f64>>,
        b: Vec<Vec<f64>>,
        horizon: usize,
        q_weight: f64,
        r_weight: f64,
    ) -> Self {
        Self {
            a,
            b,
            horizon,
            q_weight,
            r_weight,
            u_max: f64::MAX,
        }
    }

    /// Set input constraint.
    pub fn set_input_constraint(&mut self, u_max: f64) {
        self.u_max = u_max;
    }

    /// Compute the first action of the optimal input sequence.
    ///
    /// `x` is the current state and `x_ref` the reference for the predicted
    /// states `x_1, ..., x_N`; see [`optimal_sequence`](Self::optimal_sequence)
    /// for the accepted layouts.  Returns the `m` components of `u_0`, or an
    /// empty vector when the problem is ill-posed (inconsistent dimensions,
    /// non-finite data, negative weights or bound, zero horizon).
    pub fn compute(&self, x: &[f64], x_ref: &[f64]) -> Vec<f64> {
        self.optimal_sequence(x, x_ref)
            .and_then(|sequence| sequence.into_iter().next())
            .unwrap_or_default()
    }

    /// Solve the finite-horizon problem and return the whole optimal input
    /// sequence `(u_0, ..., u_{N-1})`, each of length `m`.
    ///
    /// The reference `x_ref` gives the target of the predicted states: it is
    /// either
    /// * empty (regulation to the origin),
    /// * one state of length `n`, held constant over the horizon, or
    /// * a trajectory of `k * n` entries whose block `j` is the target of
    ///   `x_{j+1}`; the last block is held if `k < N`, surplus blocks beyond
    ///   the horizon are ignored.
    ///
    /// Without a bound (`u_max >= f64::MAX`) the exact minimiser of the
    /// quadratic cost is returned.  With a bound it is the minimiser over
    /// `|u_k[i]| <= u_max`: the unconstrained solution if it already satisfies
    /// the bound, otherwise the limit of an accelerated projected-gradient
    /// iteration that stops once both the extrapolated and the previous iterate
    /// are within `1e-10` (max-norm) of the next one, or after 50 000
    /// iterations.  Every returned iterate satisfies the bound.
    ///
    /// Returns `None` when the problem is ill-posed: zero horizon or no
    /// state/input, `A` not `n x n`, `B` not `n x m`, `x_ref.len()` not a
    /// multiple of `n`, any non-finite entry, a negative or non-finite weight,
    /// or a negative/NaN input bound.
    pub fn optimal_sequence(&self, x: &[f64], x_ref: &[f64]) -> Option<Vec<Vec<f64>>> {
        let problem = self.condense(x, x_ref)?;
        let unconstrained = problem.solve_unconstrained()?;
        let solution = if self.u_max >= f64::MAX {
            unconstrained
        } else {
            let inside_box = unconstrained.iter().all(|v| v.abs() <= self.u_max);
            if inside_box {
                unconstrained
            } else {
                problem.solve_box_constrained(unconstrained, self.u_max)
            }
        };
        Some(
            solution
                .as_slice()
                .chunks(problem.input_dim)
                .map(<[f64]>::to_vec)
                .collect(),
        )
    }

    /// Validate the data and build the condensed QP of one control step.
    fn condense(&self, x: &[f64], x_ref: &[f64]) -> Option<CondensedProblem> {
        let n = x.len();
        let steps = self.horizon;
        let m = self.b.first().map(Vec::len)?;
        if n == 0 || m == 0 || steps == 0 {
            return None;
        }
        if self.a.len() != n || self.a.iter().any(|row| row.len() != n) {
            return None;
        }
        if self.b.len() != n || self.b.iter().any(|row| row.len() != m) {
            return None;
        }
        if !x_ref.len().is_multiple_of(n) {
            return None;
        }
        let weights_valid = |w: f64| w.is_finite() && w >= 0.0;
        if !weights_valid(self.q_weight) || !weights_valid(self.r_weight) {
            return None;
        }
        if self.u_max.is_nan() || self.u_max < 0.0 {
            return None;
        }
        let finite = |rows: &[Vec<f64>]| rows.iter().flatten().all(|v| v.is_finite());
        if !finite(&self.a)
            || !finite(&self.b)
            || x.iter().chain(x_ref.iter()).any(|v| !v.is_finite())
        {
            return None;
        }

        let a = DMatrix::from_fn(n, n, |i, j| self.a[i][j]);
        let b = DMatrix::from_fn(n, m, |i, j| self.b[i][j]);
        let initial = DVector::from_column_slice(x);

        // Phi x (free response) and the blocks A^k B of Gamma.
        let mut a_power = DMatrix::<f64>::identity(n, n);
        let mut a_power_b: Vec<DMatrix<f64>> = Vec::with_capacity(steps);
        let mut free_response = DVector::<f64>::zeros(steps * n);
        for k in 0..steps {
            a_power_b.push(&a_power * &b);
            a_power = &a * &a_power;
            free_response
                .rows_mut(k * n, n)
                .copy_from(&(&a_power * &initial));
        }
        let mut gamma = DMatrix::<f64>::zeros(steps * n, steps * m);
        for i in 0..steps {
            for j in 0..=i {
                gamma
                    .view_mut((i * n, j * m), (n, m))
                    .copy_from(&a_power_b[i - j]);
            }
        }

        // Reference trajectory X_ref for x_1..x_N.
        let blocks = x_ref.len() / n;
        let reference = DVector::from_fn(steps * n, |idx, _| {
            if blocks == 0 {
                0.0
            } else {
                x_ref[(idx / n).min(blocks - 1) * n + idx % n]
            }
        });

        let residual = free_response - reference;
        let hessian = gamma.tr_mul(&gamma) * self.q_weight
            + DMatrix::<f64>::identity(steps * m, steps * m) * self.r_weight;
        let gradient = gamma.tr_mul(&residual) * self.q_weight;
        Some(CondensedProblem {
            hessian,
            gradient,
            input_dim: m,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn double_integrator() -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        (
            vec![vec![1.0, 0.1], vec![0.0, 1.0]],
            vec![vec![0.005], vec![0.1]],
        )
    }

    fn to_matrix(rows: &[Vec<f64>]) -> DMatrix<f64> {
        DMatrix::from_fn(rows.len(), rows[0].len(), |i, j| rows[i][j])
    }

    /// Finite-horizon LQR gains `K_0..K_{N-1}` for the cost
    /// `sum_k q |x_{k+1}|^2 + r |u_k|^2`: backward Riccati recursion started
    /// from `P_N = q I`.
    fn riccati_gains(
        a: &DMatrix<f64>,
        b: &DMatrix<f64>,
        q: f64,
        r: f64,
        horizon: usize,
    ) -> Vec<DMatrix<f64>> {
        let n = a.nrows();
        let m = b.ncols();
        let mut p = DMatrix::<f64>::identity(n, n) * q;
        let mut gains = Vec::with_capacity(horizon);
        for _ in 0..horizon {
            let s = DMatrix::<f64>::identity(m, m) * r + b.transpose() * &p * b;
            let s_inv = s.try_inverse().expect("R + B'PB must be invertible");
            let k = s_inv * b.transpose() * &p * a;
            p = DMatrix::<f64>::identity(n, n) * q + a.transpose() * &p * a
                - a.transpose() * &p * b * &k;
            gains.push(k);
        }
        gains.reverse();
        gains
    }

    /// Cost of the input sequence by forward simulation; `references[k]` is
    /// the target of `x_{k+1}`.  Independent of the condensed formulation.
    fn simulated_cost(
        mpc: &ModelPredictiveControl,
        x0: &[f64],
        references: &[Vec<f64>],
        inputs: &[f64],
    ) -> f64 {
        let m = mpc.b[0].len();
        let mut x = x0.to_vec();
        let mut cost = 0.0;
        for (k, reference) in references.iter().enumerate() {
            let u = &inputs[k * m..(k + 1) * m];
            let mut next = vec![0.0; x.len()];
            for (i, slot) in next.iter_mut().enumerate() {
                *slot = mpc.a[i].iter().zip(&x).map(|(a, xj)| a * xj).sum::<f64>()
                    + mpc.b[i].iter().zip(u).map(|(b, uj)| b * uj).sum::<f64>();
            }
            x = next;
            cost += mpc.q_weight
                * x.iter()
                    .zip(reference)
                    .map(|(xi, ri)| (xi - ri).powi(2))
                    .sum::<f64>();
            cost += mpc.r_weight * u.iter().map(|v| v * v).sum::<f64>();
        }
        cost
    }

    /// Central finite-difference gradient of the simulated cost.
    fn numeric_gradient(
        mpc: &ModelPredictiveControl,
        x0: &[f64],
        references: &[Vec<f64>],
        inputs: &[f64],
    ) -> Vec<f64> {
        let h = 1e-5;
        (0..inputs.len())
            .map(|i| {
                let mut plus = inputs.to_vec();
                let mut minus = inputs.to_vec();
                plus[i] += h;
                minus[i] -= h;
                (simulated_cost(mpc, x0, references, &plus)
                    - simulated_cost(mpc, x0, references, &minus))
                    / (2.0 * h)
            })
            .collect()
    }

    fn flatten(sequence: &[Vec<f64>]) -> Vec<f64> {
        sequence.iter().flatten().copied().collect()
    }

    /// Asserts the KKT conditions of a minimiser over the box `|u_i| <= bound`
    /// (free components have a vanishing gradient, components at a bound have
    /// a gradient pointing out of the box) and returns the number of active
    /// bounds.
    fn assert_box_kkt(inputs: &[f64], gradient: &[f64], bound: f64) -> usize {
        let mut active = 0;
        for (i, (&u, &g)) in inputs.iter().zip(gradient.iter()).enumerate() {
            assert!(u.abs() <= bound + 1e-12, "u[{i}] = {u} violates the bound");
            if u >= bound - 1e-9 {
                active += 1;
                assert!(g <= 1e-6, "u[{i}] at +bound needs dJ/du <= 0, got {g}");
            } else if u <= -bound + 1e-9 {
                active += 1;
                assert!(g >= -1e-6, "u[{i}] at -bound needs dJ/du >= 0, got {g}");
            } else {
                assert!(g.abs() < 1e-6, "free u[{i}] needs dJ/du = 0, got {g}");
            }
        }
        active
    }

    #[test]
    fn first_input_matches_finite_horizon_riccati_lqr() {
        let (a, b) = double_integrator();
        let mpc = ModelPredictiveControl::new(a.clone(), b.clone(), 10, 1.0, 0.1);
        let x0 = [1.0, 0.5];
        let u = mpc.compute(&x0, &[0.0, 0.0]);
        assert_eq!(u.len(), 1);
        let gains = riccati_gains(&to_matrix(&a), &to_matrix(&b), 1.0, 0.1, 10);
        let x = DVector::from_column_slice(&x0);
        let expected = -(&gains[0] * x);
        assert!(
            (u[0] - expected[0]).abs() < 1e-8,
            "mpc {} vs riccati {}",
            u[0],
            expected[0]
        );
        // The truncated-gradient predecessor ended far from this optimum.
        assert!(u[0].abs() > 1e-3);
    }

    #[test]
    fn whole_sequence_matches_time_varying_riccati_feedback() {
        // Three states, two inputs, weights q = 2, r = 0.3, horizon 6.
        let a = vec![
            vec![1.0, 0.1, 0.0],
            vec![0.0, 0.95, 0.1],
            vec![0.05, 0.0, 0.9],
        ];
        let b = vec![vec![0.0, 0.1], vec![0.1, 0.0], vec![0.02, 0.05]];
        let (q, r, horizon) = (2.0, 0.3, 6);
        let mpc = ModelPredictiveControl::new(a.clone(), b.clone(), horizon, q, r);
        let x0 = [0.8, -0.4, 1.3];
        let sequence = mpc
            .optimal_sequence(&x0, &[0.0; 3])
            .expect("well-posed problem");
        assert_eq!(sequence.len(), horizon);
        let a_m = to_matrix(&a);
        let b_m = to_matrix(&b);
        let gains = riccati_gains(&a_m, &b_m, q, r, horizon);
        let mut x = DVector::from_column_slice(&x0);
        for (k, u_mpc) in sequence.iter().enumerate() {
            assert_eq!(u_mpc.len(), 2);
            let u_lqr = -(&gains[k] * &x);
            for (i, (got, want)) in u_mpc.iter().zip(u_lqr.iter()).enumerate() {
                assert!(
                    (got - want).abs() < 1e-8,
                    "step {k} input {i}: {got} vs {want}"
                );
            }
            x = &a_m * &x + &b_m * u_lqr;
        }
    }

    #[test]
    fn trajectory_reference_gradient_vanishes_at_the_optimum() {
        // A time-varying reference (ramp on position, constant velocity) over
        // the whole horizon; the stationarity check uses an independent
        // forward-simulation cost, so it also validates the block ordering.
        let (a, b) = double_integrator();
        let horizon = 8;
        let mpc = ModelPredictiveControl::new(a, b, horizon, 1.5, 0.2);
        let x0 = [0.2, -0.1];
        let references: Vec<Vec<f64>> = (0..horizon)
            .map(|k| vec![0.1 * (k + 1) as f64, 0.5])
            .collect();
        let x_ref = flatten(&references);
        let sequence = mpc
            .optimal_sequence(&x0, &x_ref)
            .expect("well-posed problem");
        let inputs = flatten(&sequence);
        let gradient = numeric_gradient(&mpc, &x0, &references, &inputs);
        for (i, g) in gradient.iter().enumerate() {
            assert!(g.abs() < 1e-6, "dJ/du[{i}] = {g}");
        }
        // Any perturbation increases the cost.
        let optimum = simulated_cost(&mpc, &x0, &references, &inputs);
        for i in 0..inputs.len() {
            let mut perturbed = inputs.clone();
            perturbed[i] += 0.05;
            assert!(simulated_cost(&mpc, &x0, &references, &perturbed) > optimum);
        }
    }

    #[test]
    fn single_state_reference_is_held_over_the_horizon() {
        let (a, b) = double_integrator();
        let horizon = 5;
        let mpc = ModelPredictiveControl::new(a, b, horizon, 1.0, 0.1);
        let x0 = [0.0, 0.0];
        let target = [1.0, 0.0];
        let held = mpc.compute(&x0, &target);
        let repeated: Vec<f64> = (0..horizon).flat_map(|_| target).collect();
        let explicit = mpc.compute(&x0, &repeated);
        assert_eq!(held.len(), 1);
        assert!((held[0] - explicit[0]).abs() < 1e-12);
        assert!(held[0] > 0.0, "must push towards the positive target");
    }

    #[test]
    fn short_trajectory_holds_its_last_block() {
        let (a, b) = double_integrator();
        let mpc = ModelPredictiveControl::new(a, b, 6, 1.0, 0.1);
        let x0 = [0.3, 0.0];
        let short = [0.0, 0.0, 1.0, 0.0];
        let padded = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0];
        let a_seq = mpc.optimal_sequence(&x0, &short).expect("well-posed");
        let b_seq = mpc.optimal_sequence(&x0, &padded).expect("well-posed");
        for (u, v) in flatten(&a_seq).iter().zip(flatten(&b_seq).iter()) {
            assert!((u - v).abs() < 1e-12);
        }
        // Surplus blocks beyond the horizon are ignored.
        let long: Vec<f64> = padded.iter().copied().chain([9.0, 9.0, 9.0, 9.0]).collect();
        let c_seq = mpc.optimal_sequence(&x0, &long).expect("well-posed");
        for (u, v) in flatten(&b_seq).iter().zip(flatten(&c_seq).iter()) {
            assert!((u - v).abs() < 1e-12);
        }
    }

    #[test]
    fn empty_reference_regulates_to_the_origin() {
        let (a, b) = double_integrator();
        let mpc = ModelPredictiveControl::new(a, b, 7, 1.0, 0.1);
        // Positive position and positive velocity: the regulator must brake.
        let x0 = [1.0, 0.5];
        let regulated = mpc.compute(&x0, &[]);
        let zero = mpc.compute(&x0, &[0.0, 0.0]);
        assert_eq!(regulated, zero);
        assert!(regulated[0] < -1e-3, "u0 = {}", regulated[0]);
    }

    #[test]
    fn closed_loop_regulates_the_double_integrator() {
        let (a, b) = double_integrator();
        let mpc = ModelPredictiveControl::new(a.clone(), b.clone(), 12, 1.0, 0.1);
        let a_m = to_matrix(&a);
        let b_m = to_matrix(&b);
        let mut x = DVector::from_column_slice(&[1.0, 0.0]);
        for _ in 0..120 {
            let u = mpc.compute(x.as_slice(), &[0.0, 0.0]);
            assert_eq!(u.len(), 1);
            x = &a_m * &x + &b_m * DVector::from_column_slice(&u);
        }
        assert!(x.norm() < 1e-3, "final state {x}");
    }

    #[test]
    fn inactive_bound_returns_the_unconstrained_solution() {
        let (a, b) = double_integrator();
        let free = ModelPredictiveControl::new(a.clone(), b.clone(), 6, 1.0, 0.1);
        let mut bounded = free.clone();
        bounded.set_input_constraint(1e3);
        let x0 = [0.4, 0.1];
        let u_free = flatten(&free.optimal_sequence(&x0, &[0.0, 0.0]).expect("ok"));
        let u_bounded = flatten(&bounded.optimal_sequence(&x0, &[0.0, 0.0]).expect("ok"));
        assert_eq!(u_free, u_bounded);
    }

    #[test]
    fn single_step_bound_equals_the_clamped_unconstrained_input() {
        // N = 1, scalar plant: J(u) = q (a x + b u - ref)^2 + r u^2 is a 1-D
        // convex quadratic, so the bounded minimiser is exactly the clamp.
        let (a, b, q, r) = (0.9, 0.5, 2.0, 0.1);
        let mut mpc = ModelPredictiveControl::new(vec![vec![a]], vec![vec![b]], 1, q, r);
        let (x, reference) = (3.0, 1.0);
        let free = -q * b * (a * x - reference) / (q * b * b + r);
        let u_free = mpc.compute(&[x], &[reference]);
        assert!((u_free[0] - free).abs() < 1e-12);
        for bound in [0.2, 0.5, 5.0] {
            mpc.set_input_constraint(bound);
            let u = mpc.compute(&[x], &[reference]);
            let expected = free.clamp(-bound, bound);
            assert!(
                (u[0] - expected).abs() < 1e-9,
                "bound {bound}: {} vs {expected}",
                u[0]
            );
        }
    }

    #[test]
    fn bounded_solution_satisfies_the_kkt_conditions() {
        let (a, b) = double_integrator();
        let horizon = 10;
        let mut mpc = ModelPredictiveControl::new(a, b, horizon, 1.0, 0.05);
        let bound = 0.4;
        mpc.set_input_constraint(bound);
        let x0 = [2.0, 0.0];
        let references: Vec<Vec<f64>> = vec![vec![0.0, 0.0]; horizon];
        let sequence = mpc
            .optimal_sequence(&x0, &[0.0, 0.0])
            .expect("well-posed problem");
        let inputs = flatten(&sequence);
        let gradient = numeric_gradient(&mpc, &x0, &references, &inputs);
        let active = assert_box_kkt(&inputs, &gradient, bound);
        assert!(active > 0, "the bound must bite in this scenario");
        assert!(active < inputs.len(), "some inputs must stay free");
    }

    #[test]
    fn bounded_multi_input_solution_satisfies_the_kkt_conditions() {
        // Three states, two inputs, constant non-zero target: the second input
        // saturates at every step while the first stays free.
        let a = vec![
            vec![1.0, 0.1, 0.0],
            vec![0.0, 0.95, 0.1],
            vec![0.05, 0.0, 0.9],
        ];
        let b = vec![vec![0.0, 0.1], vec![0.1, 0.0], vec![0.02, 0.05]];
        let horizon = 6;
        let mut mpc = ModelPredictiveControl::new(a, b, horizon, 1.0, 0.1);
        let bound = 0.8;
        mpc.set_input_constraint(bound);
        let x0 = [1.5, -0.5, 1.0];
        let target = [0.2, 0.0, -0.1];
        let references: Vec<Vec<f64>> = vec![target.to_vec(); horizon];
        let sequence = mpc
            .optimal_sequence(&x0, &target)
            .expect("well-posed problem");
        assert_eq!(sequence.len(), horizon);
        assert!(sequence.iter().all(|u| u.len() == 2));
        let inputs = flatten(&sequence);
        let gradient = numeric_gradient(&mpc, &x0, &references, &inputs);
        let active = assert_box_kkt(&inputs, &gradient, bound);
        assert!(
            active >= horizon,
            "the second input must saturate: {active}"
        );
        assert!(active < inputs.len(), "the first input must stay free");
    }

    #[test]
    fn bounded_cost_lies_between_the_free_optimum_and_the_clamped_free_solution() {
        let (a, b) = double_integrator();
        let horizon = 10;
        let mut mpc = ModelPredictiveControl::new(a, b, horizon, 1.0, 0.05);
        let x0 = [2.0, 0.0];
        let references: Vec<Vec<f64>> = vec![vec![0.0, 0.0]; horizon];
        let free = flatten(&mpc.optimal_sequence(&x0, &[0.0, 0.0]).expect("ok"));
        let bound = 0.4;
        mpc.set_input_constraint(bound);
        let bounded = flatten(&mpc.optimal_sequence(&x0, &[0.0, 0.0]).expect("ok"));
        let clamped: Vec<f64> = free.iter().map(|v| v.clamp(-bound, bound)).collect();
        let j_free = simulated_cost(&mpc, &x0, &references, &free);
        let j_bounded = simulated_cost(&mpc, &x0, &references, &bounded);
        let j_clamped = simulated_cost(&mpc, &x0, &references, &clamped);
        assert!(j_free <= j_bounded + 1e-12);
        assert!(
            j_bounded < j_clamped,
            "projected gradient must beat naive clamping: {j_bounded} vs {j_clamped}"
        );
    }

    #[test]
    fn zero_bound_forces_zero_input() {
        let (a, b) = double_integrator();
        let mut mpc = ModelPredictiveControl::new(a, b, 4, 1.0, 0.1);
        mpc.set_input_constraint(0.0);
        let sequence = mpc.optimal_sequence(&[1.0, 1.0], &[0.0, 0.0]).expect("ok");
        assert!(flatten(&sequence).iter().all(|&v| v == 0.0));
    }

    #[test]
    fn singular_hessian_falls_back_to_the_minimum_norm_solution() {
        // r = 0 and a second input that does not act on the plant: H is
        // singular, the useless input must be exactly zero and the useful one
        // must equal the single-input solution.
        let a = vec![vec![1.0, 0.1], vec![0.0, 1.0]];
        let two_inputs = ModelPredictiveControl::new(
            a.clone(),
            vec![vec![0.005, 0.0], vec![0.1, 0.0]],
            3,
            1.0,
            0.0,
        );
        let one_input = ModelPredictiveControl::new(a, vec![vec![0.005], vec![0.1]], 3, 1.0, 0.0);
        let x0 = [1.0, 0.2];
        let wide = two_inputs.optimal_sequence(&x0, &[0.0, 0.0]).expect("ok");
        let narrow = one_input.optimal_sequence(&x0, &[0.0, 0.0]).expect("ok");
        for (w, n) in wide.iter().zip(narrow.iter()) {
            assert!((w[0] - n[0]).abs() < 1e-8, "{} vs {}", w[0], n[0]);
            assert!(w[1].abs() < 1e-12, "dead input = {}", w[1]);
        }
    }

    #[test]
    fn ill_posed_problems_yield_nothing() {
        let (a, b) = double_integrator();
        let ok = ModelPredictiveControl::new(a.clone(), b.clone(), 5, 1.0, 0.1);
        let x0 = [1.0, 0.0];
        assert!(ok.optimal_sequence(&x0, &[0.0, 0.0]).is_some());
        // Zero horizon: the old code indexed an empty sequence and panicked.
        let mut zero_horizon = ok.clone();
        zero_horizon.horizon = 0;
        assert!(zero_horizon.compute(&x0, &[0.0, 0.0]).is_empty());
        // No inputs.
        let no_inputs = ModelPredictiveControl::new(a.clone(), vec![], 5, 1.0, 0.1);
        assert!(no_inputs.compute(&x0, &[0.0, 0.0]).is_empty());
        let empty_columns =
            ModelPredictiveControl::new(a.clone(), vec![vec![], vec![]], 5, 1.0, 0.1);
        assert!(empty_columns.compute(&x0, &[0.0, 0.0]).is_empty());
        // State of the wrong size / reference not a multiple of n.
        assert!(ok.compute(&[1.0], &[0.0, 0.0]).is_empty());
        assert!(ok.compute(&x0, &[0.0, 0.0, 0.0]).is_empty());
        // Ragged A and B.
        let ragged_a =
            ModelPredictiveControl::new(vec![vec![1.0, 0.1], vec![0.0]], b.clone(), 5, 1.0, 0.1);
        assert!(ragged_a.compute(&x0, &[0.0, 0.0]).is_empty());
        let wrong_b = ModelPredictiveControl::new(a.clone(), vec![vec![0.1]], 5, 1.0, 0.1);
        assert!(wrong_b.compute(&x0, &[0.0, 0.0]).is_empty());
        // Non-finite data and invalid weights / bound.
        assert!(ok.compute(&[f64::NAN, 0.0], &[0.0, 0.0]).is_empty());
        assert!(ok.compute(&x0, &[f64::INFINITY, 0.0]).is_empty());
        let negative_weight = ModelPredictiveControl::new(a.clone(), b.clone(), 5, -1.0, 0.1);
        assert!(negative_weight.compute(&x0, &[0.0, 0.0]).is_empty());
        let nan_weight = ModelPredictiveControl::new(a, b, 5, 1.0, f64::NAN);
        assert!(nan_weight.compute(&x0, &[0.0, 0.0]).is_empty());
        let mut negative_bound = ok.clone();
        negative_bound.set_input_constraint(-1.0);
        assert!(negative_bound.compute(&x0, &[0.0, 0.0]).is_empty());
        let mut nan_bound = ok;
        nan_bound.set_input_constraint(f64::NAN);
        assert!(nan_bound.compute(&x0, &[0.0, 0.0]).is_empty());
    }

    #[test]
    fn infinite_bound_means_unconstrained() {
        let (a, b) = double_integrator();
        let free = ModelPredictiveControl::new(a, b, 5, 1.0, 0.1);
        let mut unbounded = free.clone();
        unbounded.set_input_constraint(f64::INFINITY);
        let x0 = [1.0, 0.3];
        assert_eq!(
            free.compute(&x0, &[0.0, 0.0]),
            unbounded.compute(&x0, &[0.0, 0.0])
        );
    }

    #[test]
    fn output_length_is_the_input_dimension() {
        let (a, b) = double_integrator();
        let mpc = ModelPredictiveControl::new(a, b, 5, 1.0, 0.1);
        let u = mpc.compute(&[1.0, 0.0], &[0.0, 0.0]);
        assert_eq!(u.len(), 1);
    }
}
