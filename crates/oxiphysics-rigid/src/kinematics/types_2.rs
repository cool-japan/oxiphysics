//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Pfaffian constraint: A(q) · q̇ = 0.
///
/// For a wheeled robot with no-slip condition:
/// ẋ·sin θ − ẏ·cos θ = 0
#[derive(Debug, Clone)]
pub struct PfaffianConstraint {
    /// Constraint matrix A(q) (m × n)
    pub a_matrix: Vec<Vec<f64>>,
    /// Number of constraints m
    pub n_constraints: usize,
    /// Number of degrees of freedom n
    pub n_dof: usize,
}
impl PfaffianConstraint {
    /// Create a unicycle non-holonomic constraint.
    ///
    /// Coordinates: (x, y, θ). Constraint: ẋ·sin θ − ẏ·cos θ = 0.
    pub fn unicycle(theta: f64) -> Self {
        Self {
            a_matrix: vec![vec![theta.sin(), -theta.cos(), 0.0]],
            n_constraints: 1,
            n_dof: 3,
        }
    }
    /// Create a car-like constraint (bicycle model).
    ///
    /// Coordinates: (x, y, θ, φ). Constraints:
    /// 1. ẋ·sin θ − ẏ·cos θ = 0  (rear wheel no-slip)
    /// 2. ẋ·sin(θ+φ) − ẏ·cos(θ+φ) − L·θ̇·cos φ = 0  (front wheel no-slip)
    pub fn bicycle(theta: f64, phi: f64, wheelbase: f64) -> Self {
        Self {
            a_matrix: vec![
                vec![theta.sin(), -theta.cos(), 0.0, 0.0],
                vec![
                    (theta + phi).sin(),
                    -(theta + phi).cos(),
                    -wheelbase * phi.cos(),
                    0.0,
                ],
            ],
            n_constraints: 2,
            n_dof: 4,
        }
    }
    /// Check if the constraint is satisfied: ‖A · q̇‖ < tol.
    pub fn is_satisfied(&self, q_dot: &[f64], tol: f64) -> bool {
        for row in &self.a_matrix {
            let val: f64 = row.iter().zip(q_dot.iter()).map(|(a, q)| a * q).sum();
            if val.abs() > tol {
                return false;
            }
        }
        true
    }
    /// Project velocity q̇ onto the constraint manifold (null space of A).
    ///
    /// q̇_proj = (I − Aᵀ(AAᵀ)⁻¹A) · q̇, the orthogonal projection of q̇ onto
    /// the null space of A, for any number of constraints m.
    ///
    /// The projector is built from an orthonormal basis of the row space of
    /// A (Gram–Schmidt with re-orthogonalisation) instead of an explicit
    /// inverse of AAᵀ.  This is numerically better conditioned and handles
    /// rank-deficient constraint sets (duplicated or dependent rows) by
    /// dropping the dependent rows, which yields the pseudo-inverse form
    /// `I − Aᵀ(AAᵀ)⁺A`.
    ///
    /// The constraint rows are zero-padded or truncated to `q_dot.len()`, and
    /// the returned vector has the same length as `q_dot`.
    pub fn project_velocity(&self, q_dot: &[f64]) -> Vec<f64> {
        let basis = orthonormal_row_space(&self.a_matrix, q_dot.len());
        let mut projected = q_dot.to_vec();
        for _ in 0..2 {
            for direction in &basis {
                let coefficient = dot_slices(direction, &projected);
                for (value, component) in projected.iter_mut().zip(direction) {
                    *value -= coefficient * component;
                }
            }
        }
        projected
    }
    /// Degree of freedom = n − m (number of independent velocities).
    pub fn degrees_of_freedom(&self) -> usize {
        self.n_dof.saturating_sub(self.n_constraints)
    }
    /// Frobenius test: is the Pfaffian system `A(q)·q̇ = 0` integrable
    /// (holonomic) at configuration `q`?
    ///
    /// `constraint_fn(q)` must return the `m × n` constraint matrix `A(q)` as
    /// rows (shorter rows are zero-padded to `n = q.len()`).  The system is
    /// integrable iff the distribution `Δ = ker A(q)` is involutive, i.e. for
    /// every constraint covector `ωᵏ` the two-form `dωᵏ` vanishes on `Δ`:
    ///
    /// `Σᵢⱼ (∂ᵢ aᵏⱼ − ∂ⱼ aᵏᵢ) Xᵢ Yⱼ = 0` for all `X, Y ∈ Δ`.
    ///
    /// The partial derivatives `∂ᵢ aᵏⱼ` are central finite differences with
    /// step `step`; `tol` is a relative tolerance on the two-form (scaled by
    /// the largest derivative magnitude).  Note that this is strictly weaker
    /// than requiring `∂ᵢ aⱼ = ∂ⱼ aᵢ` (closedness of `ω`), which is sufficient
    /// but not necessary: a system may be integrable through an integrating
    /// factor even when `ω` itself is not closed.
    ///
    /// The unicycle and bicycle constraints are non-integrable (non-holonomic);
    /// constraints with constant coefficients are always integrable.
    pub fn is_integrable_at<F>(constraint_fn: F, q: &[f64], step: f64, tol: f64) -> bool
    where
        F: Fn(&[f64]) -> Vec<Vec<f64>>,
    {
        let n = q.len();
        let a0 = constraint_fn(q);
        let null_basis = orthonormal_null_space(&a0, n);
        if null_basis.len() < 2 {
            // A distribution of dimension 0 or 1 is always involutive.
            return true;
        }
        let h = step.abs().max(f64::EPSILON);
        // per_coordinate[i][k][j] = d a^k_j / d q_i (central differences).
        let mut per_coordinate: Vec<Vec<Vec<f64>>> = Vec::with_capacity(n);
        let mut scale = 1.0_f64;
        for i in 0..n {
            let mut forward = q.to_vec();
            let mut backward = q.to_vec();
            forward[i] += h;
            backward[i] -= h;
            let a_plus = constraint_fn(&forward);
            let a_minus = constraint_fn(&backward);
            let mut per_row = Vec::with_capacity(a_plus.len().min(a_minus.len()));
            for (plus_row, minus_row) in a_plus.iter().zip(a_minus.iter()) {
                let d: Vec<f64> = (0..n)
                    .map(|j| {
                        let plus = plus_row.get(j).copied().unwrap_or(0.0);
                        let minus = minus_row.get(j).copied().unwrap_or(0.0);
                        (plus - minus) / (2.0 * h)
                    })
                    .collect();
                scale = d.iter().fold(scale, |acc, v| acc.max(v.abs()));
                per_row.push(d);
            }
            per_coordinate.push(per_row);
        }
        let threshold = tol.abs() * scale;
        let rows = per_coordinate.iter().map(Vec::len).min().unwrap_or(0);
        // by_row[k][i][j] = d a^k_j / d q_i
        let by_row: Vec<Vec<&Vec<f64>>> = (0..rows)
            .map(|k| per_coordinate.iter().map(|rows_i| &rows_i[k]).collect())
            .collect();
        for derivative in &by_row {
            for (index_x, x) in null_basis.iter().enumerate() {
                for y in null_basis.iter().skip(index_x + 1) {
                    let mut value = 0.0;
                    for (i, (d_i, x_i)) in derivative.iter().zip(x.iter()).enumerate() {
                        for (j, (d_j, y_j)) in derivative.iter().zip(y.iter()).enumerate() {
                            // (d_i a_j - d_j a_i) x_i y_j
                            value += (d_i[j] - d_j[i]) * x_i * y_j;
                        }
                    }
                    if value.abs() > threshold {
                        return false;
                    }
                }
            }
        }
        true
    }
}
/// Dot product of two equally long slices (extra entries are ignored).
fn dot_slices(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
/// Relative rank tolerance used when orthogonalising constraint rows.
const RANK_TOLERANCE: f64 = 1e-10;
/// Orthonormalise `vectors` (modified Gram–Schmidt with a second pass) and
/// drop every vector that is numerically dependent on the previous ones.
fn orthonormalize(vectors: impl IntoIterator<Item = Vec<f64>>) -> Vec<Vec<f64>> {
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for mut v in vectors {
        let original_norm = dot_slices(&v, &v).sqrt();
        if !original_norm.is_finite() || original_norm <= 0.0 {
            continue;
        }
        for _ in 0..2 {
            for u in &basis {
                let coefficient = dot_slices(u, &v);
                for (value, component) in v.iter_mut().zip(u) {
                    *value -= coefficient * component;
                }
            }
        }
        let norm = dot_slices(&v, &v).sqrt();
        if norm > RANK_TOLERANCE * original_norm {
            for value in v.iter_mut() {
                *value /= norm;
            }
            basis.push(v);
        }
    }
    basis
}
/// Rows of `a_matrix`, zero-padded or truncated to length `n`.
fn padded_rows(a_matrix: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    a_matrix
        .iter()
        .map(|row| (0..n).map(|j| row.get(j).copied().unwrap_or(0.0)).collect())
        .collect()
}
/// Orthonormal basis of the row space of `a_matrix` (as vectors of length `n`).
fn orthonormal_row_space(a_matrix: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    orthonormalize(padded_rows(a_matrix, n))
}
/// Orthonormal basis of the null space of `a_matrix` (as vectors of length `n`).
fn orthonormal_null_space(a_matrix: &[Vec<f64>], n: usize) -> Vec<Vec<f64>> {
    let row_space = orthonormal_row_space(a_matrix, n);
    let rank = row_space.len();
    // Complete the row-space basis with the standard basis; the vectors that
    // survive after the row-space ones form an orthonormal basis of the
    // orthogonal complement (the null space).
    let completion = (0..n).map(|j| {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        e
    });
    let full = orthonormalize(row_space.into_iter().chain(completion));
    full.into_iter().skip(rank).collect()
}
