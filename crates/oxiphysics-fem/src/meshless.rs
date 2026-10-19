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
//!
//! The shape functions are built in the basis `[1, (x_I − x)/h, (y_I − y)/h]`
//! shifted to the evaluation point and scaled by the support radius, so the
//! moment matrix is well conditioned at any length scale, 1-D clouds (all `y`
//! equal) are handled with the reduced basis, and a singular moment matrix is
//! reported as `None` instead of being replaced by zero shape functions.

use crate::error::{Error, Result};

// ---------------------------------------------------------------------------
// Small dense solver (Gaussian elimination with partial pivoting)
// ---------------------------------------------------------------------------

/// Relative pivot tolerance of the moment-matrix solve.
///
/// A pivot that is not larger than this fraction of the largest entry of the
/// matrix marks the matrix as singular.  The moment matrices are assembled in
/// coordinates scaled by the support radius, so their entries are `O(1)` at every
/// length scale and a relative test is meaningful.
const PIVOT_RELATIVE_TOL: f64 = 1.0e-12;

/// Solve the small dense system `A x = b` by Gaussian elimination with partial
/// pivoting.
///
/// Returns `None` if `A` is singular, i.e. if a pivot is not larger than
/// `rel_tol · max|A_ij|` (a zero matrix is singular), or if `A` or `b` hold a
/// non-finite value.
fn solve_pivoted<const N: usize>(
    mut a: [[f64; N]; N],
    mut b: [f64; N],
    rel_tol: f64,
) -> Option<[f64; N]> {
    if a.iter().flatten().any(|v| !v.is_finite()) || b.iter().any(|v| !v.is_finite()) {
        return None;
    }
    let scale = a.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()));
    if scale == 0.0 {
        return None;
    }
    let tolerance = rel_tol * scale;
    for col in 0..N {
        let mut pivot_row = col;
        for row in (col + 1)..N {
            if a[row][col].abs() > a[pivot_row][col].abs() {
                pivot_row = row;
            }
        }
        if a[pivot_row][col].abs() <= tolerance {
            return None;
        }
        a.swap(col, pivot_row);
        b.swap(col, pivot_row);
        let upper = a[col];
        for row in (col + 1)..N {
            let factor = a[row][col] / upper[col];
            for (entry, u) in a[row].iter_mut().zip(upper.iter()).skip(col) {
                *entry -= factor * u;
            }
            b[row] -= factor * b[col];
        }
    }
    let mut x = [0.0_f64; N];
    for row in (0..N).rev() {
        let tail: f64 = ((row + 1)..N).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - tail) / a[row][row];
    }
    Some(x)
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
    /// Returns the approximated scalar field value, or `None` if the MLS shape
    /// functions do not exist at the point (see [`mls_shape_functions`]) or if
    /// `nodal_values` does not hold exactly one value per node (a shorter vector
    /// would silently drop nodes from the sum).
    pub fn evaluate(&self, x: f64, y: f64) -> Option<f64> {
        if self.nodal_values.len() != self.node_x.len() {
            return None;
        }
        let phi = mls_shape_functions(
            x,
            y,
            &self.node_x,
            &self.node_y,
            self.influence_radius,
            self.weight_type,
        )?;
        Some(
            phi.iter()
                .zip(self.nodal_values.iter())
                .map(|(p, v)| p * v)
                .sum(),
        )
    }
}

/// Polynomial basis used for the moment matrix of a node cloud.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BasisLayout {
    /// Nodes spread in both directions: linear basis `[1, x, y]`.
    Plane,
    /// All nodes share one `y` (a 1-D cloud along `x`): basis `[1, x]`.
    AlongX,
    /// All nodes share one `x` (a 1-D cloud along `y`): basis `[1, y]`.
    AlongY,
}

/// Spread of the node coordinates, relative to the support radius, below which
/// the coordinates of a cloud count as equal.
const LAYOUT_RELATIVE_TOL: f64 = 1.0e-12;

/// Decide the basis of a cloud: a 1-D cloud (all `y`, or all `x`, equal) cannot
/// support the full 2-D linear basis because its moment matrix is exactly singular.
fn basis_layout(node_x: &[f64], node_y: &[f64], h: f64) -> BasisLayout {
    let spread = |values: &[f64]| -> f64 {
        let (lo, hi) = values
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
        hi - lo
    };
    let tolerance = LAYOUT_RELATIVE_TOL * h;
    if spread(node_y) <= tolerance {
        BasisLayout::AlongX
    } else if spread(node_x) <= tolerance {
        BasisLayout::AlongY
    } else {
        BasisLayout::Plane
    }
}

/// Shape values `Φ_I = w_I · p_Iᵀ A⁻¹ e₁` for the moment matrix `A = Σ w_I p_I p_Iᵀ`.
///
/// `basis[I]` is the polynomial basis `p_I` evaluated at node `I` in coordinates
/// shifted to the evaluation point, so the basis at the evaluation point itself is
/// `e₁ = [1, 0, …]`.  `None` if `A` is singular.
fn moment_shape_values<const N: usize>(weights: &[f64], basis: &[[f64; N]]) -> Option<Vec<f64>> {
    let mut moment = [[0.0_f64; N]; N];
    for (&w, p) in weights.iter().zip(basis.iter()) {
        if w != 0.0 {
            for (row, moment_row) in moment.iter_mut().enumerate() {
                for (col, entry) in moment_row.iter_mut().enumerate() {
                    *entry += w * p[row] * p[col];
                }
            }
        }
    }
    let mut e1 = [0.0_f64; N];
    e1[0] = 1.0;
    let coefficients = solve_pivoted(moment, e1, PIVOT_RELATIVE_TOL)?;
    Some(
        weights
            .iter()
            .zip(basis.iter())
            .map(|(&w, p)| {
                let polynomial: f64 = coefficients
                    .iter()
                    .zip(p.iter())
                    .map(|(c, pk)| c * pk)
                    .sum();
                w * polynomial
            })
            .collect(),
    )
}

/// Reproducing-kernel shape values shared by MLS and RKPM.
///
/// Both approximations build `Φ_I(x) = w_I(x) · c(x)ᵀ p_I` with the correction
/// vector `c` solving `A(x) c = e₁`, `A = Σ w_I p_I p_Iᵀ`, in the basis
/// `p_I = [1, (x_I − x)/h, (y_I − y)/h]` shifted to the evaluation point and
/// scaled by the support radius `h`.  Scaling keeps the moment matrix
/// dimensionless and well conditioned for any node spacing; shifting makes the
/// right-hand side the constant vector `e₁`.
fn kernel_shape_values(
    x: f64,
    y: f64,
    node_x: &[f64],
    node_y: &[f64],
    h: f64,
    wtype: WeightType,
) -> Option<Vec<f64>> {
    if node_y.len() != node_x.len()
        || !x.is_finite()
        || !y.is_finite()
        || !h.is_finite()
        || h <= 0.0
    {
        return None;
    }
    let local: Vec<[f64; 2]> = node_x
        .iter()
        .zip(node_y.iter())
        .map(|(&px, &py)| [(px - x) / h, (py - y) / h])
        .collect();
    let weights: Vec<f64> = local
        .iter()
        .map(|p| weight_function(p[0].hypot(p[1]), wtype))
        .collect();
    match basis_layout(node_x, node_y, h) {
        BasisLayout::Plane => {
            let basis: Vec<[f64; 3]> = local.iter().map(|p| [1.0, p[0], p[1]]).collect();
            moment_shape_values(&weights, &basis)
        }
        BasisLayout::AlongX => {
            let basis: Vec<[f64; 2]> = local.iter().map(|p| [1.0, p[0]]).collect();
            moment_shape_values(&weights, &basis)
        }
        BasisLayout::AlongY => {
            let basis: Vec<[f64; 2]> = local.iter().map(|p| [1.0, p[1]]).collect();
            moment_shape_values(&weights, &basis)
        }
    }
}

/// Compute MLS shape functions Φ_I(x) for a set of 2-D nodes.
///
/// Uses the complete linear basis p(x) = \[1, x, y\], written in coordinates
/// shifted to the evaluation point and scaled by the support radius,
/// `p_I = [1, (x_I − x)/h, (y_I − y)/h]`.  The moment matrix
/// `A = Σ w_I p_I p_Iᵀ` is then dimensionless and well conditioned at any length
/// scale, and `Φ_I = w_I · p_Iᵀ A⁻¹ e₁` with `e₁ = [1, 0, 0]ᵀ` (the basis at the
/// shifted origin).  The shape functions satisfy the partition-of-unity property
/// Σ Φ_I = 1 and linear reproduction Σ Φ_I x_I = x, Σ Φ_I y_I = y.
///
/// When all nodes share one `y` coordinate (a 1-D problem, `node_y` = 0) the `y`
/// term is dropped and the basis is `[1, (x_I − x)/h]`; likewise for a common `x`.
///
/// * `x`, `y`         — evaluation point
/// * `node_x`, `node_y` — node coordinates (length N)
/// * `h`              — support radius
/// * `wtype`          — weight function type
///
/// Returns a vector of length N with the shape-function values, or `None` if they
/// do not exist at `(x, y)`: the moment matrix is singular (a 2-D cloud needs at
/// least three non-collinear nodes in range, a 1-D cloud two distinct ones), `h`
/// is not positive and finite, the point is not finite, or `node_x` and `node_y`
/// differ in length.
pub fn mls_shape_functions(
    x: f64,
    y: f64,
    node_x: &[f64],
    node_y: &[f64],
    h: f64,
    wtype: WeightType,
) -> Option<Vec<f64>> {
    kernel_shape_values(x, y, node_x, node_y, h, wtype)
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
///
/// # Errors
///
/// [`Error::InvalidInput`] if the node arrays differ in length, a background grid
/// direction has no cells, the influence radius is not positive and finite, or the
/// MLS shape functions do not exist at a quadrature point (the influence radius is
/// too small for the node distribution: every quadrature point needs at least
/// three non-collinear nodes in range).
pub fn efg_stiffness(efg: &EfgMethod) -> Result<Vec<f64>> {
    if efg.node_x.len() != efg.node_y.len() {
        return Err(Error::InvalidInput(format!(
            "node_x has {} entries but node_y has {}",
            efg.node_x.len(),
            efg.node_y.len()
        )));
    }
    if efg.n_cells_x == 0 || efg.n_cells_y == 0 {
        return Err(Error::InvalidInput(
            "the background grid needs at least one cell in each direction".to_string(),
        ));
    }
    if !efg.influence_radius.is_finite() || efg.influence_radius <= 0.0 {
        return Err(Error::InvalidInput(format!(
            "the influence radius must be positive and finite, got {}",
            efg.influence_radius
        )));
    }
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

    let shape_at = |px: f64, py: f64| -> Result<Vec<f64>> {
        mls_shape_functions(
            px,
            py,
            &efg.node_x,
            &efg.node_y,
            efg.influence_radius,
            efg.weight_type,
        )
        .ok_or_else(|| {
            Error::InvalidInput(format!(
                "the MLS shape functions are undefined at the quadrature point ({px}, {py}): \
                 the influence radius {} is too small for the node distribution",
                efg.influence_radius
            ))
        })
    };

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

                    let phi = shape_at(xq, yq)?;

                    // Approximate gradient via finite difference (step δ)
                    let delta = efg.influence_radius * 1e-5;
                    let phi_px = shape_at(xq + delta, yq)?;
                    let phi_py = shape_at(xq, yq + delta)?;

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
    Ok(k)
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
/// The raw kernel `w_I(x)` is multiplied by the polynomial correction
/// `c(x)ᵀ p_I` whose coefficients solve the moment equations `M(x) c = e₁`, so
/// that the zeroth- and first-order consistency conditions hold:
/// Σ Ψ_I = 1  and  Σ Ψ_I x_I = x  (linear completeness).  Like
/// [`mls_shape_functions`] the moment matrix is assembled in the basis
/// `p_I = [1, (x_I − x)/h, (y_I − y)/h]` shifted to the evaluation point and scaled
/// by the smoothing length, which keeps it well conditioned at any length scale,
/// and a 1-D cloud (all `y`, or all `x`, equal) uses the reduced basis.
///
/// Returns a vector of length N with the corrected kernel values, or `None` if the
/// moment matrix is singular at `(x, y)` (too few, or collinear, particles in
/// range of a 2-D cloud), the smoothing length is not positive and finite, the
/// point is not finite, or `node_x` and `node_y` differ in length.
pub fn rkpm_correction(x: f64, y: f64, basis: &RkpmBasis) -> Option<Vec<f64>> {
    kernel_shape_values(
        x,
        y,
        &basis.node_x,
        &basis.node_y,
        basis.h,
        basis.weight_type,
    )
}

/// Evaluate a scalar field approximated by RKPM at point `(x, y)`.
///
/// * `nodal_values` — scalar value at each particle (length N)
///
/// Returns `None` where the corrected kernel does not exist (see
/// [`rkpm_correction`]) or when `nodal_values` does not hold exactly one value per
/// particle (a shorter slice would silently drop particles from the sum).
pub fn rkpm_evaluate(x: f64, y: f64, basis: &RkpmBasis, nodal_values: &[f64]) -> Option<f64> {
    if nodal_values.len() != basis.num_particles() {
        return None;
    }
    let psi = rkpm_correction(x, y, basis)?;
    Some(
        psi.iter()
            .zip(nodal_values.iter())
            .map(|(p, v)| p * v)
            .sum(),
    )
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

    // ---------- pivoted solver ----------

    #[test]
    fn test_solve_pivoted_identity() {
        let a = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let b = [1.0, 2.0, 3.0];
        let x = solve_pivoted(a, b, PIVOT_RELATIVE_TOL).expect("identity is regular");
        for i in 0..3 {
            assert!((x[i] - b[i]).abs() < 1e-12);
        }
        let x2 = solve_pivoted([[1.0, 0.0], [0.0, 1.0]], [3.0, 7.0], PIVOT_RELATIVE_TOL)
            .expect("2x2 identity is regular");
        assert!((x2[0] - 3.0).abs() < 1e-12 && (x2[1] - 7.0).abs() < 1e-12);
    }

    #[test]
    fn test_solve_pivoted_known_3x3() {
        // 2x + y = 5, x + 3y + z = 10, x + y + 2z = 8
        let a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [1.0, 1.0, 2.0]];
        let b = [5.0, 10.0, 8.0];
        let x = solve_pivoted(a, b, PIVOT_RELATIVE_TOL).expect("regular system");
        for row in 0..3 {
            let lhs: f64 = (0..3).map(|col| a[row][col] * x[col]).sum();
            assert!(
                (lhs - b[row]).abs() < 1e-10,
                "row {row}: {lhs} != {}",
                b[row]
            );
        }
    }

    #[test]
    fn test_solve_pivoted_needs_a_row_exchange() {
        // The leading entry is zero: elimination without pivoting divides by zero.
        let x = solve_pivoted([[0.0, 1.0], [1.0, 0.0]], [2.0, 3.0], PIVOT_RELATIVE_TOL)
            .expect("regular after pivoting");
        assert!((x[0] - 3.0).abs() < 1e-14 && (x[1] - 2.0).abs() < 1e-14);
        // A tiny leading pivot next to a large entry (loses all digits without pivoting).
        let y = solve_pivoted([[1e-14, 1.0], [1.0, 1.0]], [1.0, 2.0], 1e-16)
            .expect("regular after pivoting");
        assert!(
            (y[0] - 1.0).abs() < 1e-12 && (y[1] - 1.0).abs() < 1e-12,
            "{y:?}"
        );
    }

    #[test]
    fn test_solve_pivoted_rejects_singular_and_non_finite_input() {
        // Rank-deficient 2x2.
        assert!(solve_pivoted([[1.0, 1.0], [1.0, 1.0]], [1.0, 2.0], PIVOT_RELATIVE_TOL).is_none());
        // Rank-deficient 3x3 (third row = first + second).
        let a = [[1.0, 2.0, 3.0], [0.0, 1.0, 4.0], [1.0, 3.0, 7.0]];
        assert!(solve_pivoted(a, [1.0, 1.0, 1.0], PIVOT_RELATIVE_TOL).is_none());
        // Zero matrix and non-finite data.
        assert!(solve_pivoted([[0.0; 3]; 3], [1.0; 3], PIVOT_RELATIVE_TOL).is_none());
        assert!(
            solve_pivoted(
                [[f64::NAN, 0.0], [0.0, 1.0]],
                [1.0, 1.0],
                PIVOT_RELATIVE_TOL
            )
            .is_none()
        );
        assert!(
            solve_pivoted(
                [[1.0, 0.0], [0.0, 1.0]],
                [f64::INFINITY, 1.0],
                PIVOT_RELATIVE_TOL
            )
            .is_none()
        );
    }

    #[test]
    fn test_solve_pivoted_tolerance_is_relative_to_the_matrix_scale() {
        // The same regular system at wildly different scales gives the same solution; the
        // previous absolute `det < 1e-14` test rejected the small-scale one.
        let a = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]];
        let b = [1.0, 2.0, 3.0];
        let reference = solve_pivoted(a, b, PIVOT_RELATIVE_TOL).expect("regular");
        for &scale in &[1.0e-30, 1.0e-12, 1.0e12, 1.0e30] {
            let scaled = a.map(|row| row.map(|v| v * scale));
            let x = solve_pivoted(scaled, b.map(|v| v * scale), PIVOT_RELATIVE_TOL)
                .expect("regular at every scale");
            for i in 0..3 {
                assert!((x[i] - reference[i]).abs() < 1e-12 * (1.0 + reference[i].abs()));
            }
        }
    }

    // ---------- MLS shape functions ----------

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
    fn test_mls_partition_of_unity_grid() {
        let (nx, ny) = uniform_nodes_2d(5, 4);
        let h = 0.5;
        let phi = mls_shape_functions(0.5, 0.5, &nx, &ny, h, WeightType::CubicSpline)
            .expect("enough nodes in range");
        let sum: f64 = phi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "PoU failed: sum = {sum}");
    }

    #[test]
    fn test_mls_partition_of_unity_2d() {
        let (nx, ny) = uniform_nodes_2d(4, 4);
        let h = 0.5;
        let phi = mls_shape_functions(0.4, 0.6, &nx, &ny, h, WeightType::CubicSpline)
            .expect("enough nodes in range");
        let sum: f64 = phi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "PoU failed: sum = {sum}");
    }

    #[test]
    fn test_mls_linear_consistency_x() {
        // Σ Φ_I * x_I ≈ x
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let h = 0.4;
        let x_eval = 0.3;
        let y_eval = 0.5;
        let phi = mls_shape_functions(x_eval, y_eval, &nx, &ny, h, WeightType::Gaussian)
            .expect("enough nodes in range");
        let repr: f64 = phi.iter().zip(nx.iter()).map(|(p, xi)| p * xi).sum();
        assert!(
            (repr - x_eval).abs() < 1e-12,
            "linear consistency x: {repr} != {x_eval}"
        );
    }

    #[test]
    fn test_mls_linear_consistency_y() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let h = 0.4;
        let x_eval = 0.5;
        let y_eval = 0.7;
        let phi = mls_shape_functions(x_eval, y_eval, &nx, &ny, h, WeightType::Gaussian)
            .expect("enough nodes in range");
        let repr: f64 = phi.iter().zip(ny.iter()).map(|(p, yi)| p * yi).sum();
        assert!(
            (repr - y_eval).abs() < 1e-12,
            "linear consistency y: {repr} != {y_eval}"
        );
    }

    #[test]
    fn test_mls_approx_constant_field() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let nodal_vals = vec![3.125_f64; nx.len()];
        let mls = MlsApproximation::new(nx, ny, nodal_vals, 0.4, WeightType::CubicSpline);
        let val = mls.evaluate(0.5, 0.5).expect("enough nodes in range");
        assert!((val - 3.125).abs() < 1e-12, "constant field: {val}");
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
        let val = mls.evaluate(xq, yq).expect("enough nodes in range");
        let expected = 2.0 * xq + 3.0 * yq;
        assert!(
            (val - expected).abs() < 1e-12,
            "linear field: {val} != {expected}"
        );
    }

    /// Regular `n x n` grid with the given spacing.
    fn grid(n: usize, spacing: f64) -> (Vec<f64>, Vec<f64>) {
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for i in 0..n {
            for j in 0..n {
                xs.push(i as f64 * spacing);
                ys.push(j as f64 * spacing);
            }
        }
        (xs, ys)
    }

    #[test]
    fn test_mls_reproduces_linear_fields_at_micrometre_spacing() {
        // 1 µm node spacing: the unscaled moment matrix has entries ~1e-12 and its
        // determinant fell below the old absolute 1e-14 threshold, so every shape
        // function silently came back as 0.
        let spacing = 1.0e-6;
        let (xs, ys) = grid(5, spacing);
        let h = 2.5 * spacing;
        for wtype in [
            WeightType::CubicSpline,
            WeightType::Gaussian,
            WeightType::QuarticSpline,
        ] {
            for &(xe, ye) in &[(2.3e-6, 1.7e-6), (0.4e-6, 3.1e-6), (2.0e-6, 2.0e-6)] {
                let phi = mls_shape_functions(xe, ye, &xs, &ys, h, wtype)
                    .expect("shape functions exist at micrometre scale");
                assert!(phi.iter().any(|p| p.abs() > 0.05), "{wtype:?}: all zero");
                let sum: f64 = phi.iter().sum();
                let rx: f64 = phi.iter().zip(xs.iter()).map(|(p, x)| p * x).sum();
                let ry: f64 = phi.iter().zip(ys.iter()).map(|(p, y)| p * y).sum();
                assert!((sum - 1.0).abs() < 1e-12, "{wtype:?}: sum {sum}");
                assert!(
                    (rx - xe).abs() < 1e-12 * spacing,
                    "{wtype:?}: Σ Φ x = {rx}, x = {xe}"
                );
                assert!(
                    (ry - ye).abs() < 1e-12 * spacing,
                    "{wtype:?}: Σ Φ y = {ry}, y = {ye}"
                );
            }
        }
    }

    #[test]
    fn test_mls_shape_functions_do_not_depend_on_the_length_scale() {
        // The same node pattern scaled by a factor must give the same Φ_I.
        let (xs1, ys1) = grid(5, 1.0);
        let reference = mls_shape_functions(2.3, 1.7, &xs1, &ys1, 2.5, WeightType::CubicSpline)
            .expect("regular");
        for &scale in &[1.0e-9, 1.0e-6, 1.0e3, 1.0e6] {
            let xs: Vec<f64> = xs1.iter().map(|v| v * scale).collect();
            let ys: Vec<f64> = ys1.iter().map(|v| v * scale).collect();
            let phi = mls_shape_functions(
                2.3 * scale,
                1.7 * scale,
                &xs,
                &ys,
                2.5 * scale,
                WeightType::CubicSpline,
            )
            .expect("regular at every scale");
            for (a, b) in phi.iter().zip(reference.iter()) {
                assert!((a - b).abs() < 1e-11, "scale {scale}: {a} vs {b}");
            }
        }
    }

    #[test]
    fn test_mls_one_dimensional_cloud_with_zero_y() {
        // The documented 1-D usage: node_y = 0.  The 3x3 basis [1, x, y] is exactly
        // singular here, so the shape functions used to come back as zeros.
        let xs: Vec<f64> = (0..6).map(|i| 0.2 * i as f64).collect();
        let ys = vec![0.0; xs.len()];
        let h = 0.45;
        let phi = mls_shape_functions(0.33, 0.0, &xs, &ys, h, WeightType::CubicSpline)
            .expect("1-D cloud uses the reduced basis");
        let sum: f64 = phi.iter().sum();
        let repr: f64 = phi.iter().zip(xs.iter()).map(|(p, x)| p * x).sum();
        assert!((sum - 1.0).abs() < 1e-13, "PoU {sum}");
        assert!((repr - 0.33).abs() < 1e-13, "linear reproduction {repr}");
        // Only the nodes within the support radius carry weight.
        for (i, p) in phi.iter().enumerate() {
            if (xs[i] - 0.33).abs() >= h {
                assert_eq!(*p, 0.0, "node {i} is outside the support");
            }
        }
        // A linear nodal field is reproduced exactly.
        let values: Vec<f64> = xs.iter().map(|x| 3.0 * x + 1.0).collect();
        let mls = MlsApproximation::new(xs, ys, values, h, WeightType::CubicSpline);
        let v = mls.evaluate(0.33, 0.0).expect("1-D cloud");
        assert!((v - (3.0 * 0.33 + 1.0)).abs() < 1e-13, "{v}");
    }

    #[test]
    fn test_mls_one_dimensional_cloud_along_y() {
        let ys: Vec<f64> = (0..6).map(|i| 0.2 * i as f64).collect();
        let xs = vec![0.0; ys.len()];
        let phi = mls_shape_functions(0.0, 0.57, &xs, &ys, 0.45, WeightType::QuarticSpline)
            .expect("1-D cloud along y uses the reduced basis");
        let sum: f64 = phi.iter().sum();
        let repr: f64 = phi.iter().zip(ys.iter()).map(|(p, y)| p * y).sum();
        assert!((sum - 1.0).abs() < 1e-13);
        assert!((repr - 0.57).abs() < 1e-13);
    }

    #[test]
    fn test_mls_returns_none_instead_of_zero_shape_functions() {
        let xs = vec![0.0, 1.0, 2.0, 0.0];
        let ys = vec![0.0, 0.0, 0.0, 10.0];
        // Collinear support in a genuinely 2-D cloud: the moment matrix has rank 2.
        assert!(mls_shape_functions(1.0, 0.2, &xs, &ys, 1.5, WeightType::CubicSpline).is_none());
        // Two nodes in range of a 2-D cloud (needs three).
        assert!(mls_shape_functions(0.5, 0.0, &xs, &ys, 0.75, WeightType::CubicSpline).is_none());
        // No node in range.
        assert!(mls_shape_functions(50.0, 50.0, &xs, &ys, 1.0, WeightType::CubicSpline).is_none());
        // One node in range of a 1-D cloud (needs two).
        let line_x = vec![0.0, 1.0, 2.0];
        let line_y = vec![0.0; 3];
        assert!(
            mls_shape_functions(0.1, 0.0, &line_x, &line_y, 0.5, WeightType::CubicSpline).is_none()
        );
        // Invalid support radius and non-finite evaluation point.
        let (gx, gy) = uniform_nodes_2d(4, 4);
        for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(mls_shape_functions(0.5, 0.5, &gx, &gy, h, WeightType::CubicSpline).is_none());
        }
        assert!(
            mls_shape_functions(f64::NAN, 0.5, &gx, &gy, 0.6, WeightType::CubicSpline).is_none()
        );
        // Empty cloud.
        assert!(mls_shape_functions(0.5, 0.5, &[], &[], 0.6, WeightType::CubicSpline).is_none());
        // And the approximant propagates the failure.
        let mls = MlsApproximation::new(xs, ys, vec![1.0; 4], 1.5, WeightType::CubicSpline);
        assert!(mls.evaluate(1.0, 0.2).is_none());
    }

    #[test]
    fn test_inconsistent_array_lengths_give_none_instead_of_a_wrong_sum() {
        let (xs, ys) = uniform_nodes_2d(5, 5);
        let n = xs.len();
        // Coordinate arrays of different length: no shape functions, and no panic.
        let short_y = &ys[..n - 1];
        assert!(
            mls_shape_functions(0.5, 0.5, &xs, short_y, 0.4, WeightType::CubicSpline).is_none()
        );
        let ragged = RkpmBasis::new(xs.clone(), short_y.to_vec(), 0.4, WeightType::CubicSpline);
        assert!(rkpm_correction(0.5, 0.5, &ragged).is_none());

        // One nodal value per node is required; a short vector used to be zipped with the
        // shape functions and silently dropped the remaining nodes from the sum.
        let values = vec![2.0; n];
        let complete = MlsApproximation::new(
            xs.clone(),
            ys.clone(),
            values.clone(),
            0.4,
            WeightType::CubicSpline,
        );
        let v = complete.evaluate(0.5, 0.5).expect("complete data");
        assert!((v - 2.0).abs() < 1e-12, "{v}");
        let truncated = MlsApproximation::new(
            xs.clone(),
            ys.clone(),
            values[..n - 1].to_vec(),
            0.4,
            WeightType::CubicSpline,
        );
        assert!(truncated.evaluate(0.5, 0.5).is_none());
        let padded = MlsApproximation::new(
            xs.clone(),
            ys.clone(),
            vec![2.0; n + 1],
            0.4,
            WeightType::CubicSpline,
        );
        assert!(padded.evaluate(0.5, 0.5).is_none());

        let basis = RkpmBasis::new(xs, ys, 0.4, WeightType::CubicSpline);
        let v = rkpm_evaluate(0.5, 0.5, &basis, &values).expect("complete data");
        assert!((v - 2.0).abs() < 1e-12, "{v}");
        assert!(rkpm_evaluate(0.5, 0.5, &basis, &values[..n - 1]).is_none());
        assert!(rkpm_evaluate(0.5, 0.5, &basis, &[]).is_none());
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
            0.6,
            WeightType::CubicSpline,
        );
        let k = efg_stiffness(&efg).expect("support covers every quadrature point");
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
        let k = efg_stiffness(&efg).expect("support covers every quadrature point");
        let dof = 2 * n;
        for i in 0..dof {
            for j in i + 1..dof {
                let diff = (k[i * dof + j] - k[j * dof + i]).abs();
                assert!(diff < 1e-4, "asymmetry at ({i},{j}): {diff}");
            }
        }
    }

    #[test]
    fn test_efg_stiffness_annihilates_rigid_translations() {
        let (nx, ny) = uniform_nodes_2d(4, 4);
        let n = nx.len();
        let efg = EfgMethod::new(
            nx,
            ny,
            3,
            3,
            [0.0, 1.0, 0.0, 1.0],
            1e6,
            0.3,
            0.7,
            WeightType::CubicSpline,
        );
        let k = efg_stiffness(&efg).expect("support covers every quadrature point");
        let dof = 2 * n;
        let k_max = k.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        assert!(k_max > 1.0);
        for direction in 0..2 {
            let u: Vec<f64> = (0..dof)
                .map(|i| if i % 2 == direction { 1.0 } else { 0.0 })
                .collect();
            for row in 0..dof {
                let f: f64 = (0..dof).map(|col| k[row * dof + col] * u[col]).sum();
                assert!(
                    f.abs() < 1e-6 * k_max,
                    "row {row}: rigid translation force {f}"
                );
            }
        }
    }

    #[test]
    fn test_efg_stiffness_reports_an_undersized_influence_radius() {
        // With radius 0.4 the corner quadrature points of the 3x3 cloud see a single node,
        // so the shape functions do not exist; they used to be silently replaced by zeros.
        let (nx, ny) = uniform_nodes_2d(3, 3);
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
        match efg_stiffness(&efg) {
            Err(Error::InvalidInput(message)) => {
                assert!(message.contains("quadrature point"), "{message}");
            }
            other => panic!("expected InvalidInput, got {other:?}"),
        }
    }

    #[test]
    fn test_efg_stiffness_rejects_malformed_methods() {
        let (nx, ny) = uniform_nodes_2d(3, 3);
        let base = EfgMethod::new(
            nx.clone(),
            ny.clone(),
            2,
            2,
            [0.0, 1.0, 0.0, 1.0],
            1e6,
            0.3,
            0.6,
            WeightType::CubicSpline,
        );
        let mut no_cells = base.clone();
        no_cells.n_cells_x = 0;
        assert!(matches!(
            efg_stiffness(&no_cells),
            Err(Error::InvalidInput(_))
        ));
        let mut bad_radius = base.clone();
        bad_radius.influence_radius = f64::NAN;
        assert!(matches!(
            efg_stiffness(&bad_radius),
            Err(Error::InvalidInput(_))
        ));
        let mut short_y = base;
        short_y.node_y.pop();
        assert!(matches!(
            efg_stiffness(&short_y),
            Err(Error::InvalidInput(_))
        ));
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
        let psi = rkpm_correction(0.5, 0.5, &basis).expect("enough particles in range");
        let sum: f64 = psi.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12, "RKPM PoU failed: {sum}");
    }

    #[test]
    fn test_rkpm_linear_consistency_x() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let x_eval = 0.3;
        let y_eval = 0.6;
        let basis = RkpmBasis::new(nx.clone(), ny, 0.4, WeightType::CubicSpline);
        let psi = rkpm_correction(x_eval, y_eval, &basis).expect("enough particles in range");
        let repr: f64 = psi.iter().zip(nx.iter()).map(|(p, xi)| p * xi).sum();
        assert!((repr - x_eval).abs() < 1e-12, "RKPM linear x: {repr}");
    }

    #[test]
    fn test_rkpm_linear_consistency_y() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let x_eval = 0.6;
        let y_eval = 0.4;
        let basis = RkpmBasis::new(nx, ny.clone(), 0.4, WeightType::CubicSpline);
        let psi = rkpm_correction(x_eval, y_eval, &basis).expect("enough particles in range");
        let repr: f64 = psi.iter().zip(ny.iter()).map(|(p, yi)| p * yi).sum();
        assert!((repr - y_eval).abs() < 1e-12, "RKPM linear y: {repr}");
    }

    #[test]
    fn test_rkpm_evaluate_constant() {
        let (nx, ny) = uniform_nodes_2d(5, 5);
        let vals = vec![2.71_f64; nx.len()];
        let basis = RkpmBasis::new(nx, ny, 0.4, WeightType::Gaussian);
        let val = rkpm_evaluate(0.5, 0.5, &basis, &vals).expect("enough particles in range");
        assert!((val - 2.71).abs() < 1e-12, "RKPM constant field: {val}");
    }

    #[test]
    fn test_rkpm_num_particles() {
        let (nx, ny) = uniform_nodes_2d(4, 4);
        let n = nx.len();
        let basis = RkpmBasis::new(nx, ny, 0.3, WeightType::QuarticSpline);
        assert_eq!(basis.num_particles(), n);
    }

    #[test]
    fn test_rkpm_reproduces_linear_fields_at_micrometre_spacing() {
        let spacing = 1.0e-6;
        let (xs, ys) = grid(5, spacing);
        let values: Vec<f64> = xs
            .iter()
            .zip(ys.iter())
            .map(|(x, y)| 4.0e6 * x - 2.0e6 * y + 0.5)
            .collect();
        let basis = RkpmBasis::new(xs, ys, 2.5 * spacing, WeightType::CubicSpline);
        for &(xe, ye) in &[(2.3e-6, 1.7e-6), (0.4e-6, 3.1e-6)] {
            let psi = rkpm_correction(xe, ye, &basis).expect("exists at micrometre scale");
            let sum: f64 = psi.iter().sum();
            assert!((sum - 1.0).abs() < 1e-12, "sum {sum}");
            let v = rkpm_evaluate(xe, ye, &basis, &values).expect("exists");
            let expected = 4.0e6 * xe - 2.0e6 * ye + 0.5;
            assert!((v - expected).abs() < 1e-11, "{v} vs {expected}");
        }
    }

    #[test]
    fn test_rkpm_one_dimensional_particles_with_zero_y() {
        let xs: Vec<f64> = (0..8).map(|i| 0.125 * i as f64).collect();
        let ys = vec![0.0; xs.len()];
        let values: Vec<f64> = xs.iter().map(|x| 2.0 - 5.0 * x).collect();
        let basis = RkpmBasis::new(xs, ys, 0.3, WeightType::CubicSpline);
        let psi = rkpm_correction(0.41, 0.0, &basis).expect("1-D cloud uses the reduced basis");
        assert!((psi.iter().sum::<f64>() - 1.0).abs() < 1e-13);
        let v = rkpm_evaluate(0.41, 0.0, &basis, &values).expect("1-D cloud");
        assert!((v - (2.0 - 5.0 * 0.41)).abs() < 1e-13, "{v}");
    }

    #[test]
    fn test_rkpm_returns_none_where_the_moment_matrix_is_singular() {
        // Collinear particles in range of a 2-D cloud, then an empty neighbourhood.
        let basis = RkpmBasis::new(
            vec![0.0, 1.0, 2.0, 0.0],
            vec![0.0, 0.0, 0.0, 10.0],
            1.5,
            WeightType::CubicSpline,
        );
        assert!(rkpm_correction(1.0, 0.2, &basis).is_none());
        assert!(rkpm_evaluate(1.0, 0.2, &basis, &[1.0; 4]).is_none());
        assert!(rkpm_correction(100.0, 100.0, &basis).is_none());
        let bad_h = RkpmBasis::new(
            vec![0.0, 1.0],
            vec![0.0, 1.0],
            -1.0,
            WeightType::CubicSpline,
        );
        assert!(rkpm_correction(0.5, 0.5, &bad_h).is_none());
    }
}
