#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Meshless / meshfree methods module.
//!
//! Implements:
//! - **Moving Least Squares (MLS)** shape function approximation
//! - **Element-Free Galerkin (EFG)** method for boundary value problems
//! - **Reproducing Kernel Particle Method (RKPM)** with consistency correction
//!
//! All computations use plain `f64` arrays — no external linear-algebra crate.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Small linear-algebra helpers (2×2 and 3×3, row-major)
// ---------------------------------------------------------------------------

/// Solve 2×2 linear system A x = b.  Returns None if singular.
fn solve2(a: [[f64; 2]; 2], b: [f64; 2]) -> Option<[f64; 2]> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < 1e-14 {
        return None;
    }
    Some([
        (b[0] * a[1][1] - b[1] * a[0][1]) / det,
        (a[0][0] * b[1] - a[1][0] * b[0]) / det,
    ])
}

/// Solve 3×3 linear system A x = b.  Returns None if singular.
fn solve3(a: [[f64; 3]; 3], b: [f64; 3]) -> Option<[f64; 3]> {
    let det = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if det.abs() < 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    let x0 = inv_det
        * (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            + b[1] * (a[0][2] * a[2][1] - a[0][1] * a[2][2])
            + b[2] * (a[0][1] * a[1][2] - a[0][2] * a[1][1]));
    let x1 = inv_det
        * (b[0] * (a[1][2] * a[2][0] - a[1][0] * a[2][2])
            + b[1] * (a[0][0] * a[2][2] - a[0][2] * a[2][0])
            + b[2] * (a[0][2] * a[1][0] - a[0][0] * a[1][2]));
    let x2 = inv_det
        * (b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
            + b[1] * (a[0][1] * a[2][0] - a[0][0] * a[2][1])
            + b[2] * (a[0][0] * a[1][1] - a[0][1] * a[1][0]));
    Some([x0, x1, x2])
}

// ---------------------------------------------------------------------------
// Weight functions
// ---------------------------------------------------------------------------

/// Weight function type for MLS and EFG.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum WeightType {
    /// Cubic spline weight function.
    CubicSpline,
    /// Gaussian weight function.
    Gaussian,
    /// Quartic (4th-order) spline weight.
    QuarticSpline,
}

/// Evaluate a scalar weight function at normalised distance `r_norm = |x - xI| / dI`.
///
/// * `r_norm`  — normalised distance (non-negative)
/// * `wtype`   — kernel type to use
///
/// Returns the weight value w(r_norm) ≥ 0.
pub fn weight_function(r_norm: f64, wtype: WeightType) -> f64 {
    let s = r_norm.abs();
    match wtype {
        WeightType::CubicSpline => {
            if s <= 0.5 {
                2.0 / 3.0 - 4.0 * s * s + 4.0 * s * s * s
            } else if s <= 1.0 {
                4.0 / 3.0 - 4.0 * s + 4.0 * s * s - 4.0 / 3.0 * s * s * s
            } else {
                0.0
            }
        }
        WeightType::Gaussian => {
            if s <= 1.0 {
                (-s * s / 0.3).exp()
            } else {
                0.0
            }
        }
        WeightType::QuarticSpline => {
            if s <= 1.0 {
                let t = 1.0 - s;
                t * t * t * t
            } else {
                0.0
            }
        }
    }
}

/// Derivative of the weight function with respect to the normalised distance.
///
/// Returns dw/ds evaluated at `r_norm`.
pub fn weight_function_deriv(r_norm: f64, wtype: WeightType) -> f64 {
    let s = r_norm.abs();
    match wtype {
        WeightType::CubicSpline => {
            if s <= 0.5 {
                -8.0 * s + 12.0 * s * s
            } else if s <= 1.0 {
                -4.0 + 8.0 * s - 4.0 * s * s
            } else {
                0.0
            }
        }
        WeightType::Gaussian => {
            if s <= 1.0 {
                -2.0 * s / 0.3 * (-s * s / 0.3).exp()
            } else {
                0.0
            }
        }
        WeightType::QuarticSpline => {
            if s <= 1.0 {
                let t = 1.0 - s;
                -4.0 * t * t * t
            } else {
                0.0
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Moving Least Squares (MLS) approximation
// ---------------------------------------------------------------------------

/// Moving Least Squares approximation data.
///
/// Stores the list of nodes (support centres), their nodal values, and
/// the influence radius used to build weight functions.
#[derive(Debug, Clone)]
pub struct MlsApproximation {
    /// Node x-coordinates.
    pub node_x: Vec<f64>,
    /// Node y-coordinates (set to 0.0 for 1-D problems).
    pub node_y: Vec<f64>,
    /// Nodal values associated with the nodes.
    pub nodal_values: Vec<f64>,
    /// Influence radius for weight functions.
    pub influence_radius: f64,
    /// Weight function type.
    pub weight_type: WeightType,
}

impl MlsApproximation {
    /// Construct a new MLS approximation with the given nodes and parameters.
    pub fn new(
        node_x: Vec<f64>,
        node_y: Vec<f64>,
        nodal_values: Vec<f64>,
        influence_radius: f64,
        weight_type: WeightType,
    ) -> Self {
        Self {
            node_x,
            node_y,
            nodal_values,
            influence_radius,
            weight_type,
        }
    }

    /// Evaluate the MLS approximant at point `(x, y)`.
    ///
    /// Returns the approximated scalar field value.
    pub fn evaluate(&self, x: f64, y: f64) -> f64 {
        let phi = mls_shape_functions(
            x,
            y,
            &self.node_x,
            &self.node_y,
            self.influence_radius,
            self.weight_type,
        );
        phi.iter()
            .zip(self.nodal_values.iter())
            .map(|(p, v)| p * v)
            .sum()
    }
}

/// Compute MLS shape functions Φ_I(x) for a set of 2-D nodes.
///
/// Uses a linear (complete first-order) polynomial basis p(x) = \[1, x, y\].
/// The shape functions satisfy the partition-of-unity property: Σ Φ_I = 1.
///
/// * `x`, `y`         — evaluation point
/// * `node_x`, `node_y` — node coordinates (length N)
/// * `h`              — support radius
/// * `wtype`          — weight function type
///
/// Returns a vector of length N with the shape-function values.
pub fn mls_shape_functions(
    x: f64,
    y: f64,
    node_x: &[f64],
    node_y: &[f64],
    h: f64,
    wtype: WeightType,
) -> Vec<f64> {
    let n = node_x.len();
    assert_eq!(node_y.len(), n);

    // p(x) = [1, x, y]  (linear basis, m=3)
    let p_eval = [1.0_f64, x, y];

    // Moment matrix A = Σ w_I p_I p_I^T  (3×3)
    let mut a = [[0.0_f64; 3]; 3];
    // B matrix: columns are w_I * p_I  (3 × N)
    let mut weights = vec![0.0_f64; n];

    for i in 0..n {
        let dx = x - node_x[i];
        let dy = y - node_y[i];
        let dist = (dx * dx + dy * dy).sqrt();
        let r_norm = dist / h;
        let w = weight_function(r_norm, wtype);
        weights[i] = w;
        let p_i = [1.0_f64, node_x[i], node_y[i]];
        for row in 0..3 {
            for col in 0..3 {
                a[row][col] += w * p_i[row] * p_i[col];
            }
        }
    }

    // Solve A c = p_eval  for each node contribution
    // Shape function: Φ_I(x) = p^T A^{-1} (w_I p_I)
    let mut phi = vec![0.0_f64; n];

    // Compute A^{-1} p_eval  via Cramer's rule / solve3
    if let Some(coeff) = solve3(a, p_eval) {
        for i in 0..n {
            let p_i = [1.0_f64, node_x[i], node_y[i]];
            let dot: f64 = coeff.iter().zip(p_i.iter()).map(|(c, p)| c * p).sum();
            phi[i] = weights[i] * dot;
        }
    }
    phi
}

// ---------------------------------------------------------------------------
// Element-Free Galerkin (EFG) method
// ---------------------------------------------------------------------------

/// Element-Free Galerkin method data structure.
///
/// Stores the nodal cloud, background integration grid, and material
/// parameters for a 2-D linear-elasticity EFG problem.
#[derive(Debug, Clone)]
pub struct EfgMethod {
    /// Node x-coordinates.
    pub node_x: Vec<f64>,
    /// Node y-coordinates.
    pub node_y: Vec<f64>,
    /// Number of background cells in x direction.
    pub n_cells_x: usize,
    /// Number of background cells in y direction.
    pub n_cells_y: usize,
    /// Domain extent: \[x_min, x_max, y_min, y_max\].
    pub domain: [f64; 4],
    /// Young's modulus.
    pub young: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Influence radius for MLS.
    pub influence_radius: f64,
    /// Weight function type.
    pub weight_type: WeightType,
}

impl EfgMethod {
    /// Construct a new EFG method instance.
    pub fn new(
        node_x: Vec<f64>,
        node_y: Vec<f64>,
        n_cells_x: usize,
        n_cells_y: usize,
        domain: [f64; 4],
        young: f64,
        poisson: f64,
        influence_radius: f64,
        weight_type: WeightType,
    ) -> Self {
        Self {
            node_x,
            node_y,
            n_cells_x,
            n_cells_y,
            domain,
            young,
            poisson,
            influence_radius,
            weight_type,
        }
    }

    /// Number of nodes.
    pub fn num_nodes(&self) -> usize {
        self.node_x.len()
    }
}

/// Compute the plane-stress constitutive matrix D (3×3, Voigt notation).
///
/// Returns `[[D\]]` for the given Young's modulus and Poisson's ratio.
pub fn plane_stress_d_matrix(young: f64, poisson: f64) -> [[f64; 3]; 3] {
    let factor = young / (1.0 - poisson * poisson);
    [
        [factor, factor * poisson, 0.0],
        [factor * poisson, factor, 0.0],
        [0.0, 0.0, factor * (1.0 - poisson) / 2.0],
    ]
}

/// Assemble the EFG global stiffness matrix using background Gauss quadrature.
///
/// Uses a 2×2 Gauss rule in each background cell.  Returns the (2N × 2N)
/// stiffness matrix as a flat row-major vector.
pub fn efg_stiffness(efg: &EfgMethod) -> Vec<f64> {
    let n = efg.num_nodes();
    let dof = 2 * n;
    let mut k = vec![0.0_f64; dof * dof];

    let d = plane_stress_d_matrix(efg.young, efg.poisson);

    let x_min = efg.domain[0];
    let x_max = efg.domain[1];
    let y_min = efg.domain[2];
    let y_max = efg.domain[3];
    let hx = (x_max - x_min) / efg.n_cells_x as f64;
    let hy = (y_max - y_min) / efg.n_cells_y as f64;

    // 2-point Gauss rule: points ±1/√3, weight 1
    let gp = [-1.0_f64 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
    let gw = [1.0_f64, 1.0];

    for ci in 0..efg.n_cells_x {
        for cj in 0..efg.n_cells_y {
            let x_c = x_min + (ci as f64 + 0.5) * hx;
            let y_c = y_min + (cj as f64 + 0.5) * hy;

            for (xi, wi_x) in gp.iter().zip(gw.iter()) {
                let (xi, wi_x) = (*xi, *wi_x);
                for (eta, wi_y) in gp.iter().zip(gw.iter()) {
                    let (eta, wi_y) = (*eta, *wi_y);
                    let xq = x_c + 0.5 * hx * xi;
                    let yq = y_c + 0.5 * hy * eta;
                    let jac = 0.25 * hx * hy;
                    let wq = wi_x * wi_y * jac;

                    let phi = mls_shape_functions(
                        xq,
                        yq,
                        &efg.node_x,
                        &efg.node_y,
                        efg.influence_radius,
                        efg.weight_type,
                    );

                    // Approximate gradient via finite difference (step δ)
                    let delta = efg.influence_radius * 1e-5;
                    let phi_px = mls_shape_functions(
                        xq + delta,
                        yq,
                        &efg.node_x,
                        &efg.node_y,
                        efg.influence_radius,
                        efg.weight_type,
                    );
                    let phi_py = mls_shape_functions(
                        xq,
                        yq + delta,
                        &efg.node_x,
                        &efg.node_y,
                        efg.influence_radius,
                        efg.weight_type,
                    );

                    // B matrix column for node I: [dΦ/dx, 0; 0, dΦ/dy; dΦ/dy, dΦ/dx]
                    for i in 0..n {
                        let dphi_i_dx = (phi_px[i] - phi[i]) / delta;
                        let dphi_i_dy = (phi_py[i] - phi[i]) / delta;
                        let bi = [[dphi_i_dx, 0.0], [0.0, dphi_i_dy], [dphi_i_dy, dphi_i_dx]];

                        for j in 0..n {
                            let dphi_j_dx = (phi_px[j] - phi[j]) / delta;
                            let dphi_j_dy = (phi_py[j] - phi[j]) / delta;
                            let bj = [[dphi_j_dx, 0.0], [0.0, dphi_j_dy], [dphi_j_dy, dphi_j_dx]];

                            // k_ij += w * B_i^T D B_j  (2×2 block)
                            for a in 0..2 {
                                for b in 0..2 {
                                    let mut val = 0.0;
                                    for p in 0..3 {
                                        for q in 0..3 {
                                            val += bi[p][a] * d[p][q] * bj[q][b];
                                        }
                                    }
                                    k[(2 * i + a) * dof + (2 * j + b)] += wq * val;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    k
}

// ---------------------------------------------------------------------------
// Reproducing Kernel Particle Method (RKPM)
// ---------------------------------------------------------------------------

/// RKPM basis data for a nodal cloud.
///
/// Stores node positions and the correction factors needed to enforce
/// reproducing-kernel consistency up to first order.
#[derive(Debug, Clone)]
pub struct RkpmBasis {
    /// Node x-coordinates.
    pub node_x: Vec<f64>,
    /// Node y-coordinates (0.0 for 1-D).
    pub node_y: Vec<f64>,
    /// Smoothing length (kernel support radius).
    pub h: f64,
    /// Weight function type.
    pub weight_type: WeightType,
}

impl RkpmBasis {
    /// Construct a new RKPM basis.
    pub fn new(node_x: Vec<f64>, node_y: Vec<f64>, h: f64, weight_type: WeightType) -> Self {
        Self {
            node_x,
            node_y,
            h,
            weight_type,
        }
    }

    /// Number of particles.
    pub fn num_particles(&self) -> usize {
        self.node_x.len()
    }
}

/// Compute RKPM corrected kernel values Ψ_I(x) at point `(x, y)`.
///
/// Applies zeroth- and first-order consistency corrections so that
/// Σ Ψ_I = 1  and  Σ Ψ_I x_I = x  (linear completeness).
///
/// Returns a vector of length N with the corrected kernel values.
pub fn rkpm_correction(x: f64, y: f64, basis: &RkpmBasis) -> Vec<f64> {
    let n = basis.num_particles();
    let mut raw_w = vec![0.0_f64; n];

    for i in 0..n {
        let dx = x - basis.node_x[i];
        let dy = y - basis.node_y[i];
        let r = (dx * dx + dy * dy).sqrt() / basis.h;
        raw_w[i] = weight_function(r, basis.weight_type);
    }

    // Moment conditions: sum w_I = C0, sum w_I x_I = C1 x, sum w_I y_I = C2 y
    // Solve for correction coefficients c = [c0, c1, c2] from:
    //   M c = e_1  where M_kl = sum_I w_I p_k(x_I) p_l(x_I)
    //   and e_1 = [1, x, y]
    let mut m = [[0.0_f64; 3]; 3];
    for i in 0..n {
        let pi = [1.0_f64, basis.node_x[i], basis.node_y[i]];
        let w = raw_w[i];
        for row in 0..3 {
            for col in 0..3 {
                m[row][col] += w * pi[row] * pi[col];
            }
        }
    }

    let e1 = [1.0_f64, x, y];
    let c = solve3(m, e1).unwrap_or([0.0, 0.0, 0.0]);

    let mut psi = vec![0.0_f64; n];
    for i in 0..n {
        let pi = [1.0_f64, basis.node_x[i], basis.node_y[i]];
        let dot: f64 = c.iter().zip(pi.iter()).map(|(ci, pi_k)| ci * pi_k).sum();
        psi[i] = raw_w[i] * dot;
    }
    psi
}

/// Evaluate a scalar field approximated by RKPM at point `(x, y)`.
///
/// * `nodal_values` — scalar value at each particle (length N)
pub fn rkpm_evaluate(x: f64, y: f64, basis: &RkpmBasis, nodal_values: &[f64]) -> f64 {
    let psi = rkpm_correction(x, y, basis);
    psi.iter()
        .zip(nodal_values.iter())
        .map(|(p, v)| p * v)
        .sum()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- weight functions ----------

    #[test]
    fn test_cubic_weight_at_zero() {
        let w = weight_function(0.0, WeightType::CubicSpline);
        assert!((w - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_cubic_weight_at_boundary() {
        let w = weight_function(1.0, WeightType::CubicSpline);
        assert!(w.abs() < 1e-12);
    }

    #[test]
    fn test_cubic_weight_outside_support() {
        let w = weight_function(1.5, WeightType::CubicSpline);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_gaussian_weight_at_zero() {
        let w = weight_function(0.0, WeightType::Gaussian);
        assert!((w - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_gaussian_weight_outside() {
        let w = weight_function(1.1, WeightType::Gaussian);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_quartic_weight_at_zero() {
        let w = weight_function(0.0, WeightType::QuarticSpline);
        assert!((w - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_quartic_weight_at_boundary() {
        let w = weight_function(1.0, WeightType::QuarticSpline);
        assert!(w.abs() < 1e-12);
    }

    #[test]
    fn test_quartic_weight_outside() {
        let w = weight_function(2.0, WeightType::QuarticSpline);
        assert_eq!(w, 0.0);
    }

    #[test]
    fn test_weight_non_negative() {
        for wtype in [
            WeightType::CubicSpline,
            WeightType::Gaussian,
            WeightType::QuarticSpline,
        ] {
            for i in 0..20 {
                let s = i as f64 / 20.0 * 1.2;
                assert!(weight_function(s, wtype) >= 0.0);
            }
        }
    }

    #[test]
    fn test_weight_deriv_cubic_at_zero() {
        let dw = weight_function_deriv(0.0, WeightType::CubicSpline);
        assert!(dw.abs() < 1e-12, "dw at 0 should be 0, got {dw}");
    }

    // ---------- solve helpers ----------

    #[test]
    fn test_solve2_identity() {
        let a = [[1.0, 0.0], [0.0, 1.0]];
        let b = [3.0, 7.0];
        let x = solve2(a, b).unwrap();
        assert!((x[0] - 3.0).abs() < 1e-12);
        assert!((x[1] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_solve2_singular() {
        let a = [[1.0, 1.0], [1.0, 1.0]];
        let b = [1.0, 2.0];
        assert!(solve2(a, b).is_none());
    }

    #[test]
    fn test_solve3_identity() {
        let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let b = [1.0, 2.0, 3.0];
        let x = solve3(a, b).unwrap();
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_solve3_known() {
        // 2x + y = 5, x + 3y + z = 10, x + y + 2z = 8
        let a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [1.0, 1.0, 2.0]];
        let b = [5.0, 10.0, 8.0];
        let x = solve3(a, b).unwrap();
        // Check A x == b
        for row in 0..3 {
            let lhs: f64 = (0..3).map(|col| a[row][col] * x[col]).sum();
            assert!(
                (lhs - b[row]).abs() < 1e-10,
                "row {row}: {lhs} != {}",
                b[row]
            );
        }
    }

    // ---------- MLS shape functions ----------

    fn uniform_nodes_1d(n: usize, x0: f64, x1: f64) -> (Vec<f64>, Vec<f64>) {
        let node_x: Vec<f64> = (0..n)
            .map(|i| x0 + i as f64 * (x1 - x0) / (n - 1) as f64)
            .collect();
        let node_y = vec![0.0_f64; n];
        (node_x, node_y)
    }

    fn uniform_nodes_2d(nx: usize, ny: usize) -> (Vec<f64>, Vec<f64>) {
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                xs.push(i as f64 / (nx - 1) as f64);
                ys.push(j as f64 / (ny - 1) as f64);
            }
        }
        (xs, ys)
    }

    #[test]
    fn test_mls_partition_of_unity_1d() {
        // Use a 2D grid (required: basis [1,x,y] needs nodes spread in both x and y)
        let (nx, ny) = uniform_nodes_2d(5, 4);
        let h = 0.5;
        let phi = mls_shape_functions(0.5, 0.5, &nx, &ny, h, WeightType::CubicSpline);
        let sum: f64 = phi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5, "PoU failed: sum = {sum}");
    }

    #[test]
    fn test_mls_partition_of_unity_2d() {
        let (nx, ny) = uniform_nodes_2d(4, 4);
        let h = 0.5;
        let phi = mls_shape_functions(0.4, 0.6, &nx, &ny, h, WeightType::CubicSpline);
        let sum: f64 = phi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-5, "PoU failed: sum = {sum}");
    }

    #[test]
    fn test_mls_linear_consistency_x() {
        // Σ Φ_I * x_I ≈ x
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let h = 0.4;
        let x_eval = 0.3;
        let y_eval = 0.5;
        let phi = mls_shape_functions(x_eval, y_eval, &nx, &ny, h, WeightType::Gaussian);
        let repr: f64 = phi.iter().zip(nx.iter()).map(|(p, xi)| p * xi).sum();
        assert!(
            (repr - x_eval).abs() < 1e-5,
            "linear consistency x: {repr} != {x_eval}"
        );
    }

    #[test]
    fn test_mls_linear_consistency_y() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let h = 0.4;
        let x_eval = 0.5;
        let y_eval = 0.7;
        let phi = mls_shape_functions(x_eval, y_eval, &nx, &ny, h, WeightType::Gaussian);
        let repr: f64 = phi.iter().zip(ny.iter()).map(|(p, yi)| p * yi).sum();
        assert!(
            (repr - y_eval).abs() < 1e-5,
            "linear consistency y: {repr} != {y_eval}"
        );
    }

    #[test]
    fn test_mls_approx_constant_field() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let nodal_vals = vec![3.125_f64; nx.len()];
        let mls = MlsApproximation::new(nx, ny, nodal_vals, 0.4, WeightType::CubicSpline);
        let val = mls.evaluate(0.5, 0.5);
        assert!((val - 3.125).abs() < 1e-5, "constant field: {val}");
    }

    #[test]
    fn test_mls_approx_linear_field() {
        let (nx, ny) = uniform_nodes_2d(6, 6);
        // f(x,y) = 2x + 3y
        let nodal_vals: Vec<f64> = nx
            .iter()
            .zip(ny.iter())
            .map(|(xi, yi)| 2.0 * xi + 3.0 * yi)
            .collect();
        let mls = MlsApproximation::new(nx, ny, nodal_vals, 0.4, WeightType::CubicSpline);
        let (xq, yq) = (0.4, 0.6);
        let val = mls.evaluate(xq, yq);
        let expected = 2.0 * xq + 3.0 * yq;
        assert!(
            (val - expected).abs() < 1e-4,
            "linear field: {val} != {expected}"
        );
    }

    // ---------- EFG ----------

    #[test]
    fn test_efg_stiffness_size() {
        let (nx, ny) = uniform_nodes_2d(3, 3);
        let n = nx.len();
        let efg = EfgMethod::new(
            nx,
            ny,
            2,
            2,
            [0.0, 1.0, 0.0, 1.0],
            1e6,
            0.3,
            0.4,
            WeightType::CubicSpline,
        );
        let k = efg_stiffness(&efg);
        assert_eq!(k.len(), (2 * n) * (2 * n));
    }

    #[test]
    fn test_efg_stiffness_symmetry() {
        let (nx, ny) = uniform_nodes_2d(3, 3);
        let n = nx.len();
        let efg = EfgMethod::new(
            nx,
            ny,
            2,
            2,
            [0.0, 1.0, 0.0, 1.0],
            1e6,
            0.3,
            0.5,
            WeightType::CubicSpline,
        );
        let k = efg_stiffness(&efg);
        let dof = 2 * n;
        for i in 0..dof {
            for j in i + 1..dof {
                let diff = (k[i * dof + j] - k[j * dof + i]).abs();
                assert!(diff < 1e-4, "asymmetry at ({i},{j}): {diff}");
            }
        }
    }

    #[test]
    fn test_plane_stress_d_matrix() {
        let d = plane_stress_d_matrix(1e6, 0.0);
        // With ν=0: D = diag(E, E, E/2)
        assert!((d[0][0] - 1e6).abs() < 1.0);
        assert!(d[0][1].abs() < 1.0);
        assert!((d[2][2] - 0.5e6).abs() < 1.0);
    }

    // ---------- RKPM ----------

    #[test]
    fn test_rkpm_partition_of_unity() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let basis = RkpmBasis::new(nx, ny, 0.4, WeightType::CubicSpline);
        let psi = rkpm_correction(0.5, 0.5, &basis);
        let sum: f64 = psi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-6, "RKPM PoU failed: {sum}");
    }

    #[test]
    fn test_rkpm_linear_consistency_x() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let x_eval = 0.3;
        let y_eval = 0.6;
        let basis = RkpmBasis::new(nx.clone(), ny, 0.4, WeightType::CubicSpline);
        let psi = rkpm_correction(x_eval, y_eval, &basis);
        let repr: f64 = psi.iter().zip(nx.iter()).map(|(p, xi)| p * xi).sum();
        assert!((repr - x_eval).abs() < 1e-5, "RKPM linear x: {repr}");
    }

    #[test]
    fn test_rkpm_linear_consistency_y() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let x_eval = 0.6;
        let y_eval = 0.4;
        let basis = RkpmBasis::new(nx, ny.clone(), 0.4, WeightType::CubicSpline);
        let psi = rkpm_correction(x_eval, y_eval, &basis);
        let repr: f64 = psi.iter().zip(ny.iter()).map(|(p, yi)| p * yi).sum();
        assert!((repr - y_eval).abs() < 1e-5, "RKPM linear y: {repr}");
    }

    #[test]
    fn test_rkpm_evaluate_constant() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let vals = vec![2.71_f64; nx.len()];
        let basis = RkpmBasis::new(nx, ny, 0.4, WeightType::Gaussian);
        let val = rkpm_evaluate(0.5, 0.5, &basis, &vals);
        assert!((val - 2.71).abs() < 1e-5, "RKPM constant field: {val}");
    }

    #[test]
    fn test_rkpm_num_particles() {
        let (nx, ny) = uniform_nodes_2d(4, 4);
        let n = nx.len();
        let basis = RkpmBasis::new(nx, ny, 0.3, WeightType::QuarticSpline);
        assert_eq!(basis.num_particles(), n);
    }
}
