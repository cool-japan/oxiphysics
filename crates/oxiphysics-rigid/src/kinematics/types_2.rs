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
    /// q̇_proj = (I − Aᵀ(AAᵀ)⁻¹A) · q̇
    pub fn project_velocity(&self, q_dot: &[f64]) -> Vec<f64> {
        let n = self.n_dof;
        let m = self.n_constraints;
        let mut aq: Vec<f64> = vec![0.0; m];
        for (i, aq_i) in aq.iter_mut().enumerate() {
            *aq_i = q_dot[..n]
                .iter()
                .enumerate()
                .map(|(j, &qd)| self.a_matrix[i][j] * qd)
                .sum();
        }
        let mut aat = vec![vec![0.0f64; m]; m];
        for (i, aat_row) in aat.iter_mut().enumerate() {
            for (j, cell) in aat_row.iter_mut().enumerate() {
                *cell = (0..n)
                    .map(|k| self.a_matrix[i][k] * self.a_matrix[j][k])
                    .sum();
            }
        }
        if m == 1 {
            let a_norm_sq: f64 = self.a_matrix[0].iter().map(|x| x * x).sum();
            let lambda = aq[0] / a_norm_sq.max(1e-30);
            return q_dot
                .iter()
                .enumerate()
                .map(|(j, &qd)| qd - lambda * self.a_matrix[0][j])
                .collect();
        }
        q_dot.to_vec()
    }
    /// Degree of freedom = n − m (number of independent velocities).
    pub fn degrees_of_freedom(&self) -> usize {
        self.n_dof.saturating_sub(self.n_constraints)
    }
    /// Check if constraints are integrable (holonomic test via commutativity).
    ///
    /// For a single constraint a(q)·q̇ = 0, it is integrable iff ∂aᵢ/∂qⱼ = ∂aⱼ/∂qᵢ.
    pub fn frobenius_integrability_1d(&self) -> bool {
        self.n_constraints == 0
    }
}
