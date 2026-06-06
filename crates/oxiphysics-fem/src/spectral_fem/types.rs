//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Chebyshev spectral collocation differentiation matrix.
///
/// The N×N matrix D_{ij} = l'_j(x_i) at Chebyshev-Gauss-Lobatto nodes.
pub struct ChebyshevDiffMatrix {
    /// Number of CGL nodes N+1.
    pub n_nodes: usize,
    /// CGL nodes x_k = cos(kπ/N).
    pub nodes: Vec<f64>,
    /// Differentiation matrix.
    pub d_matrix: Vec<Vec<f64>>,
}
impl ChebyshevDiffMatrix {
    /// Construct the Chebyshev differentiation matrix of order `n`.
    pub fn new(n: usize) -> Self {
        let nodes: Vec<f64> = (0..=n)
            .map(|k| (k as f64 * std::f64::consts::PI / n as f64).cos())
            .collect();
        let np1 = n + 1;
        let c: Vec<f64> = (0..=n)
            .map(|k| if k == 0 || k == n { 2.0 } else { 1.0 })
            .collect();
        let mut d = vec![vec![0.0; np1]; np1];
        for i in 0..np1 {
            for j in 0..np1 {
                if i == j {
                    continue;
                }
                let ci = c[i];
                let cj = c[j];
                d[i][j] =
                    (ci / cj) * (if (i + j) % 2 == 0 { 1.0 } else { -1.0 }) / (nodes[i] - nodes[j]);
            }
        }
        for (i, row) in d.iter_mut().enumerate() {
            let sum: f64 = row
                .iter()
                .enumerate()
                .filter(|&(j, _)| j != i)
                .map(|(_, &v)| v)
                .sum();
            row[i] = -sum;
        }
        d[0][0] = (2.0 * (n as f64) * (n as f64) + 1.0) / 6.0;
        d[n][n] = -(2.0 * (n as f64) * (n as f64) + 1.0) / 6.0;
        ChebyshevDiffMatrix {
            n_nodes: np1,
            nodes,
            d_matrix: d,
        }
    }
    /// Apply the differentiation matrix to a vector of nodal values.
    pub fn differentiate(&self, u: &[f64]) -> Vec<f64> {
        let n = self.n_nodes;
        let mut du = vec![0.0; n];
        for (i, du_i) in du.iter_mut().enumerate() {
            for (j, &u_j) in u.iter().enumerate() {
                *du_i += self.d_matrix[i][j] * u_j;
            }
        }
        du
    }
    /// Row sum (should be zero for constant field).
    pub fn row_sum(&self, i: usize) -> f64 {
        self.d_matrix[i].iter().sum()
    }
}
/// 2-D spectral element on a quadrilateral using tensor-product GLL nodes.
pub struct SpectralElement2D {
    /// Polynomial degree in each direction (square element).
    pub degree: usize,
    /// GLL nodes in 1-D.
    pub nodes_1d: Vec<f64>,
    /// GLL weights in 1-D.
    pub weights_1d: Vec<f64>,
    /// 1-D derivative matrix.
    pub d1: Vec<Vec<f64>>,
    /// Physical dimensions \[lx, ly\].
    pub dims: [f64; 2],
}
impl SpectralElement2D {
    /// Construct a 2-D spectral element.
    pub fn new(degree: usize, dims: [f64; 2]) -> Self {
        let (nodes_1d, weights_1d) = gll_nodes_weights(degree);
        let d1 = spectral_derivative_matrix(&nodes_1d);
        SpectralElement2D {
            degree,
            nodes_1d,
            weights_1d,
            d1,
            dims,
        }
    }
    /// Number of nodes per direction.
    pub fn n_nodes_1d(&self) -> usize {
        self.nodes_1d.len()
    }
    /// Total number of nodes.
    pub fn n_nodes(&self) -> usize {
        let n = self.n_nodes_1d();
        n * n
    }
    /// Jacobian determinant (constant for affine quad): J = (lx/2)*(ly/2).
    pub fn jacobian(&self) -> f64 {
        (self.dims[0] / 2.0) * (self.dims[1] / 2.0)
    }
    /// Diagonal mass matrix using tensor product GLL quadrature.
    pub fn mass_matrix_diagonal(&self) -> Vec<f64> {
        let n = self.n_nodes_1d();
        let j = self.jacobian();
        let mut m = Vec::with_capacity(n * n);
        for i in 0..n {
            for j_idx in 0..n {
                m.push(self.weights_1d[i] * self.weights_1d[j_idx] * j);
            }
        }
        m
    }
    /// Assemble element stiffness matrix (2-D, Laplace operator).
    ///
    /// Returns a dense `(n²) × (n²)` matrix.
    pub fn stiffness_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_nodes_1d();
        let nn = n * n;
        let jac = self.jacobian();
        let jx = self.dims[0] / 2.0;
        let jy = self.dims[1] / 2.0;
        let mut k = vec![vec![0.0; nn]; nn];
        let mut kx = vec![vec![0.0; n]; n];
        let mut ky = vec![vec![0.0; n]; n];
        for a in 0..n {
            for b in 0..n {
                for q in 0..n {
                    kx[a][b] += self.weights_1d[q] * self.d1[q][a] * self.d1[q][b];
                    ky[a][b] += self.weights_1d[q] * self.d1[q][a] * self.d1[q][b];
                }
            }
        }
        for (i, kx_row) in kx.iter().enumerate() {
            for (j_idx, &w_jidx) in self.weights_1d.iter().enumerate() {
                let row = i * n + j_idx;
                for (kk, &kx_ikk) in kx_row.iter().enumerate() {
                    for (l, &ky_jl) in ky[j_idx].iter().enumerate() {
                        let col = kk * n + l;
                        let term_x =
                            (jy / jx) * kx_ikk * w_jidx * (if j_idx == l { 1.0 } else { 0.0 });
                        let term_y = (jx / jy)
                            * self.weights_1d[i]
                            * (if i == kk { 1.0 } else { 0.0 })
                            * ky_jl;
                        k[row][col] = (term_x + term_y) * jac / (jx * jy);
                    }
                }
            }
        }
        k
    }
}
/// Spectral convergence rate estimator.
///
/// Measures the L2 error between a spectral solution and a reference function
/// on a sequence of polynomial degrees to estimate exponential convergence.
pub struct SpectralConvergence {
    /// GLL nodes used for the current degree.
    pub nodes: Vec<f64>,
    /// GLL weights.
    pub weights: Vec<f64>,
}
impl SpectralConvergence {
    /// Construct for degree `p`.
    pub fn new(degree: usize) -> Self {
        let (nodes, weights) = gll_nodes_weights(degree);
        SpectralConvergence { nodes, weights }
    }
    /// L2 error of the polynomial interpolant of `f` on GLL nodes.
    ///
    /// Evaluates on a fine equidistant grid of 4*(degree+1) points and returns
    /// √(mean((u_h - f)²)) for a coarser estimate of off-node interpolation error.
    pub fn l2_interp_error<F: Fn(f64) -> f64>(&self, f: F) -> f64 {
        let nodal_values: Vec<f64> = self.nodes.iter().map(|&x| f(x)).collect();
        let interp = SpectralInterpolation::new(self.nodes.clone());
        let n_fine = 4 * (self.nodes.len() + 1);
        let mut err_sq = 0.0;
        for k in 0..n_fine {
            let xi = -1.0 + (2.0 * k as f64 + 1.0) / n_fine as f64;
            let u_h = interp.interpolate(&nodal_values, xi);
            let err = u_h - f(xi);
            err_sq += err * err;
        }
        (err_sq / n_fine as f64).sqrt()
    }
    /// Estimate spectral convergence rate between degree `p1` and `p2` for function `f`.
    ///
    /// Returns the slope in log(error) vs log(p) (p-convergence rate).
    pub fn convergence_rate<F: Fn(f64) -> f64 + Copy>(f: F, p1: usize, p2: usize) -> f64 {
        let sc1 = SpectralConvergence::new(p1);
        let sc2 = SpectralConvergence::new(p2);
        let e1 = sc1.l2_interp_error(f).max(1e-300);
        let e2 = sc2.l2_interp_error(f).max(1e-300);
        (e2.ln() - e1.ln()) / ((p2 as f64).ln() - (p1 as f64).ln())
    }
}
/// 1-D spectral element convection operator with upwind flux.
///
/// Solves: ∂u/∂t + c ∂u/∂x = 0 using the spectral derivative matrix.
pub struct SpectralConvection {
    /// Spectral element.
    pub element: SpectralElement1D,
    /// Convection velocity c.
    pub velocity: f64,
}
impl SpectralConvection {
    /// Construct a spectral convection operator.
    pub fn new(degree: usize, length: f64, velocity: f64) -> Self {
        SpectralConvection {
            element: SpectralElement1D::new(degree, length),
            velocity,
        }
    }
    /// Compute the convective flux divergence c * du/dx at GLL nodes.
    pub fn flux_divergence(&self, u: &[f64]) -> Vec<f64> {
        let du = self.element.differentiate(u);
        du.iter().map(|&d| self.velocity * d).collect()
    }
    /// Explicit Euler time step: u^{n+1} = u^n - dt * c * du/dx.
    pub fn step_euler(&self, u: &[f64], dt: f64) -> Vec<f64> {
        let rhs = self.flux_divergence(u);
        u.iter()
            .zip(rhs.iter())
            .map(|(&un, &r)| un - dt * r)
            .collect()
    }
    /// Upwind numerical flux at the left boundary (x = −1).
    pub fn upwind_flux_left(&self, u_left: f64, u_right: f64) -> f64 {
        if self.velocity >= 0.0 {
            self.velocity * u_left
        } else {
            self.velocity * u_right
        }
    }
}
/// Chebyshev spectral method using discrete cosine transform (DCT-II / DCT-III).
///
/// Represents a function on Chebyshev-Gauss-Lobatto nodes using Chebyshev
/// polynomial coefficients.
pub struct ChebyshevSpectral {
    /// Number of Chebyshev modes N (number of nodes = N+1).
    pub n_modes: usize,
    /// CGL nodes: x_k = cos(k π / N), k = 0..N.
    pub nodes: Vec<f64>,
}
impl ChebyshevSpectral {
    /// Construct a Chebyshev spectral method with `n` modes (n+1 nodes).
    pub fn new(n_modes: usize) -> Self {
        let n = n_modes;
        let nodes: Vec<f64> = (0..=n)
            .map(|k| (k as f64 * std::f64::consts::PI / n as f64).cos())
            .collect();
        ChebyshevSpectral { n_modes, nodes }
    }
    /// Forward DCT: convert nodal values to Chebyshev coefficients.
    ///
    /// Uses the Chebyshev expansion on CGL nodes:
    /// â_k = (2/N) Σ_{j=0}^{N} '' f_j cos(kπj/N)
    /// where '' means the j=0 and j=N terms are halved.
    pub fn to_coefficients(&self, values: &[f64]) -> Vec<f64> {
        let n = self.n_modes;
        let np1 = n + 1;
        let mut coeff = vec![0.0; np1];
        for (k, coeff_k) in coeff.iter_mut().enumerate() {
            let sum: f64 = values
                .iter()
                .enumerate()
                .map(|(j, &vj)| {
                    let wj = if j == 0 || j == n { 0.5 } else { 1.0 };
                    wj * vj * (k as f64 * std::f64::consts::PI * j as f64 / n as f64).cos()
                })
                .sum();
            *coeff_k = 2.0 / n as f64 * sum;
        }
        coeff
    }
    /// Inverse DCT: convert Chebyshev coefficients to nodal values.
    ///
    /// f_j = Σ_{k=0}^{N} '' â_k cos(kπj/N)
    /// where '' means the k=0 and k=N coefficients are halved.
    pub fn to_nodal(&self, coeff: &[f64]) -> Vec<f64> {
        let n = self.n_modes;
        let np1 = n + 1;
        let mut values = vec![0.0; np1];
        for (j, val_j) in values.iter_mut().enumerate() {
            *val_j = coeff
                .iter()
                .enumerate()
                .map(|(k, &ck)| {
                    let wk = if k == 0 || k == n { 0.5 } else { 1.0 };
                    wk * ck * (k as f64 * std::f64::consts::PI * j as f64 / n as f64).cos()
                })
                .sum();
        }
        values
    }
    /// Differentiate nodal values via Chebyshev spectral differentiation.
    ///
    /// Converts to coefficients, applies the Chebyshev derivative recurrence,
    /// and converts back to nodal values.
    pub fn differentiate(&self, values: &[f64]) -> Vec<f64> {
        let n = self.n_modes;
        let np1 = n + 1;
        let coeff = self.to_coefficients(values);
        let mut a_std = coeff.clone();
        a_std[0] /= 2.0;
        a_std[n] /= 2.0;
        let mut da_std = vec![0.0; np1];
        for k in (0..n).rev() {
            let next2 = if k + 2 <= n { da_std[k + 2] } else { 0.0 };
            let next1 = if k < n { a_std[k + 1] } else { 0.0 };
            let val = 2.0 * (k + 1) as f64 * next1 + next2;
            da_std[k] = if k == 0 { val / 2.0 } else { val };
        }
        let mut da = da_std;
        da[0] *= 2.0;
        da[n] *= 2.0;
        self.to_nodal(&da)
    }
    /// Aliasing dealiasing filter: zero-out modes above a threshold fraction.
    ///
    /// Common choice: retain the lower 2/3 of modes (2/3 rule).
    pub fn dealias_filter(&self, coeff: &mut [f64], retain_fraction: f64) {
        let keep = ((self.n_modes as f64 * retain_fraction) as usize).min(self.n_modes);
        for c in coeff.iter_mut().skip(keep) {
            *c = 0.0;
        }
    }
}
/// Spectral element solver for the Helmholtz equation: Δu + k²u = f.
///
/// Assembles the system matrix A = K − k² M and a simple diagonal solver.
pub struct SpectralHelmholtz {
    /// Spectral element (1-D).
    pub element: SpectralElement1D,
    /// Wave number k.
    pub k: f64,
}
impl SpectralHelmholtz {
    /// Construct a spectral Helmholtz solver.
    pub fn new(degree: usize, length: f64, k: f64) -> Self {
        SpectralHelmholtz {
            element: SpectralElement1D::new(degree, length),
            k,
        }
    }
    /// Assemble the Helmholtz system matrix A = K − k² diag(M).
    pub fn system_matrix(&self) -> Vec<Vec<f64>> {
        let kmat = self.element.stiffness_matrix();
        let mdiag = self.element.mass_matrix_diagonal();
        let n = kmat.len();
        let mut a = kmat;
        for i in 0..n {
            a[i][i] -= self.k * self.k * mdiag[i];
        }
        a
    }
    /// Apply Dirichlet boundary conditions by zeroing rows/columns at endpoints.
    pub fn apply_dirichlet(&self, a: &mut [Vec<f64>], rhs: &mut [f64]) {
        let n = a.len();
        if n == 0 {
            return;
        }
        for v in a[0].iter_mut() {
            *v = 0.0;
        }
        a[0][0] = 1.0;
        rhs[0] = 0.0;
        for v in a[n - 1].iter_mut() {
            *v = 0.0;
        }
        a[n - 1][n - 1] = 1.0;
        rhs[n - 1] = 0.0;
    }
}
/// SEM solver for the 1-D elliptic PDE: −d²u/dx² = f on \[a, b\]
/// with Dirichlet boundary conditions u(a) = u(b) = 0.
///
/// Uses direct Gaussian elimination on the assembled global stiffness matrix.
pub struct SemEllipticSolver {
    /// hp mesh.
    pub mesh: HpMesh1D,
    /// Global-local mapping.
    pub mapping: GlobalLocalMapping,
    /// Assembled global stiffness matrix (sparse stored as dense for small N).
    pub k_global: Vec<Vec<f64>>,
    /// Right-hand-side vector.
    pub rhs: Vec<f64>,
}
impl SemEllipticSolver {
    /// Construct and assemble the global stiffness matrix.
    pub fn new(mesh: HpMesh1D) -> Self {
        let mapping = GlobalLocalMapping::from_hp_mesh(&mesh);
        let ng = mapping.n_global;
        let mut k_global = vec![vec![0.0; ng]; ng];
        let mut rhs = vec![0.0; ng];
        for e in 0..mesh.n_elem() {
            let elem = SpectralElement1D::new(mesh.degrees[e], mesh.lengths[e]);
            let k_loc = elem.stiffness_matrix();
            for (i, k_loc_row) in k_loc.iter().enumerate() {
                let gi = mapping.conn[e][i];
                for (j, &k_loc_ij) in k_loc_row.iter().enumerate() {
                    let gj = mapping.conn[e][j];
                    k_global[gi][gj] += k_loc_ij;
                }
            }
        }
        Self::apply_dirichlet_bc(&mut k_global, &mut rhs, 0);
        Self::apply_dirichlet_bc(&mut k_global, &mut rhs, ng - 1);
        SemEllipticSolver {
            mesh,
            mapping,
            k_global,
            rhs,
        }
    }
    /// Apply Dirichlet BC at DOF `idx` (zero value).
    fn apply_dirichlet_bc(k: &mut [Vec<f64>], rhs: &mut [f64], idx: usize) {
        for v in k[idx].iter_mut() {
            *v = 0.0;
        }
        k[idx][idx] = 1.0;
        rhs[idx] = 0.0;
        let n = k.len();
        for i in 0..n {
            if i != idx {
                rhs[i] -= k[i][idx] * rhs[idx];
                k[i][idx] = 0.0;
            }
        }
    }
    /// Load the RHS from a source function f(x) via GLL quadrature.
    pub fn load_rhs<F: Fn(f64) -> f64>(&mut self, f: F) {
        let ng = self.mapping.n_global;
        self.rhs = vec![0.0; ng];
        let mut x_offset = 0.0;
        for e in 0..self.mesh.n_elem() {
            let elem = SpectralElement1D::new(self.mesh.degrees[e], self.mesh.lengths[e]);
            let j = elem.jacobian();
            let _n_loc = elem.nodes.len();
            let m_diag = elem.mass_matrix_diagonal();
            for (i, (&node, &mdi)) in elem.nodes.iter().zip(m_diag.iter()).enumerate() {
                let x_phys = x_offset + (node + 1.0) / 2.0 * self.mesh.lengths[e];
                let gi = self.mapping.conn[e][i];
                self.rhs[gi] += mdi / j * f(x_phys);
            }
            x_offset += self.mesh.lengths[e];
        }
        self.rhs[0] = 0.0;
        let last = self.mapping.n_global - 1;
        self.rhs[last] = 0.0;
    }
    /// Solve the assembled system K u = rhs using Gaussian elimination.
    pub fn solve(&self) -> Vec<f64> {
        let n = self.k_global.len();
        if n == 0 {
            return vec![];
        }
        let mut a = self.k_global.clone();
        let mut b = self.rhs.clone();
        for col in 0..n {
            let mut pivot_row = col;
            let mut max_val = a[col][col].abs();
            for (row, a_row) in a.iter().enumerate().skip(col + 1) {
                if a_row[col].abs() > max_val {
                    max_val = a_row[col].abs();
                    pivot_row = row;
                }
            }
            a.swap(col, pivot_row);
            b.swap(col, pivot_row);
            let diag = a[col][col];
            if diag.abs() < 1e-14 {
                continue;
            }
            let col_slice: Vec<f64> = a[col][col..].to_vec();
            for row in col + 1..n {
                let factor = a[row][col] / diag;
                for (off, &cv) in col_slice.iter().enumerate() {
                    a[row][col + off] -= factor * cv;
                }
                b[row] -= factor * b[col];
            }
        }
        let mut x = vec![0.0; n];
        for i in (0..n).rev() {
            let mut sum = b[i];
            for j in i + 1..n {
                sum -= a[i][j] * x[j];
            }
            x[i] = if a[i][i].abs() > 1e-14 {
                sum / a[i][i]
            } else {
                0.0
            };
        }
        x
    }
}
/// Diagonal mass matrix for spectral elements (GLL mass lumping).
pub struct SpectralMassMatrix {
    /// Diagonal entries.
    pub diagonal: Vec<f64>,
}
impl SpectralMassMatrix {
    /// Construct from a 1-D GLL rule and element length.
    pub fn from_1d_element(degree: usize, length: f64) -> Self {
        let elem = SpectralElement1D::new(degree, length);
        SpectralMassMatrix {
            diagonal: elem.mass_matrix_diagonal(),
        }
    }
    /// Construct from a 2-D GLL element.
    pub fn from_2d_element(degree: usize, dims: [f64; 2]) -> Self {
        let elem = SpectralElement2D::new(degree, dims);
        SpectralMassMatrix {
            diagonal: elem.mass_matrix_diagonal(),
        }
    }
    /// Scale the mass matrix by factor `s`.
    pub fn scale(&mut self, s: f64) {
        for d in &mut self.diagonal {
            *d *= s;
        }
    }
    /// Invert the diagonal mass matrix.
    pub fn inverse(&self) -> Vec<f64> {
        self.diagonal
            .iter()
            .map(|&m| if m.abs() > 1e-300 { 1.0 / m } else { 0.0 })
            .collect()
    }
}
/// Stiffness matrix assembly for 1-D spectral elements.
pub struct SpectralStiffness {
    /// Dense stiffness matrix.
    pub matrix: Vec<Vec<f64>>,
    /// Polynomial degree.
    pub degree: usize,
}
impl SpectralStiffness {
    /// Assemble stiffness for a 1-D element.
    pub fn from_1d_element(degree: usize, length: f64) -> Self {
        let elem = SpectralElement1D::new(degree, length);
        SpectralStiffness {
            matrix: elem.stiffness_matrix(),
            degree,
        }
    }
    /// Retrieve entry K\[i\]\[j\].
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.matrix[i][j]
    }
    /// Row sum (should be zero for pure Laplace).
    pub fn row_sum(&self, i: usize) -> f64 {
        self.matrix[i].iter().sum()
    }
}
/// Descriptor for an hp-refined spectral element mesh.
///
/// Each element is described by its physical length and polynomial degree,
/// enabling independent h- and p-refinement per element.
pub struct HpMesh1D {
    /// Physical length of each element (m).
    pub lengths: Vec<f64>,
    /// Polynomial degree of each element.
    pub degrees: Vec<usize>,
}
impl HpMesh1D {
    /// Construct a uniform hp mesh of `n_elem` elements with given length and degree.
    pub fn uniform(n_elem: usize, total_length: f64, degree: usize) -> Self {
        let h = total_length / n_elem as f64;
        HpMesh1D {
            lengths: vec![h; n_elem],
            degrees: vec![degree; n_elem],
        }
    }
    /// Number of elements.
    pub fn n_elem(&self) -> usize {
        self.lengths.len()
    }
    /// Total number of unique DOFs (global, assuming C^0 continuity).
    pub fn n_dofs(&self) -> usize {
        if self.lengths.is_empty() {
            return 0;
        }
        self.degrees.iter().sum::<usize>() + 1
    }
    /// Apply p-refinement: increase degree of element `i` by `delta`.
    pub fn p_refine(&mut self, elem_idx: usize, delta: usize) {
        if elem_idx < self.degrees.len() {
            self.degrees[elem_idx] += delta;
        }
    }
    /// Apply h-refinement: split element `i` into two equal halves.
    pub fn h_refine(&mut self, elem_idx: usize) {
        if elem_idx < self.lengths.len() {
            let h = self.lengths[elem_idx] / 2.0;
            let deg = self.degrees[elem_idx];
            self.lengths.remove(elem_idx);
            self.degrees.remove(elem_idx);
            self.lengths.insert(elem_idx, h);
            self.lengths.insert(elem_idx + 1, h);
            self.degrees.insert(elem_idx, deg);
            self.degrees.insert(elem_idx + 1, deg);
        }
    }
}
/// Catalogue of Legendre polynomial utilities.
pub struct LegendrePolynomial {
    /// Polynomial degree.
    pub n: usize,
}
impl LegendrePolynomial {
    /// Create a Legendre polynomial of degree `n`.
    pub fn new(n: usize) -> Self {
        LegendrePolynomial { n }
    }
    /// Evaluate P_n(x).
    pub fn evaluate(&self, x: f64) -> f64 {
        legendre_pn(self.n, x)
    }
    /// Evaluate derivative P'_n(x).
    pub fn derivative(&self, x: f64) -> f64 {
        legendre_evaluate(self.n, x).1
    }
    /// Compute the n zeros of P_n.
    pub fn zeros(&self) -> Vec<f64> {
        legendre_zeros(self.n)
    }
}
/// Spectral boundary element collocation for 2-D Laplace.
pub struct SpectralBoundaryIntegral {
    /// Boundary node positions.
    pub nodes: Vec<[f64; 2]>,
    /// Outward unit normals at each node.
    pub normals: Vec<[f64; 2]>,
}
impl SpectralBoundaryIntegral {
    /// Construct from boundary nodes and normals.
    pub fn new(nodes: Vec<[f64; 2]>, normals: Vec<[f64; 2]>) -> Self {
        SpectralBoundaryIntegral { nodes, normals }
    }
    /// Assemble the BEM influence matrix H_{ij} = ∂G/∂n(x_i, y_j).
    pub fn h_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let mut h = vec![vec![0.0; n]; n];
        for (i, h_row) in h.iter_mut().enumerate() {
            for (j, h_ij) in h_row.iter_mut().enumerate() {
                if i == j {
                    continue;
                }
                *h_ij = green_normal_derivative_2d(self.nodes[i], self.nodes[j], self.normals[j]);
            }
        }
        h
    }
    /// Assemble the BEM G matrix G_{ij} = G(x_i, y_j).
    pub fn g_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let mut g = vec![vec![0.0; n]; n];
        for (i, g_row) in g.iter_mut().enumerate() {
            for (j, g_ij) in g_row.iter_mut().enumerate() {
                if i == j {
                    continue;
                }
                *g_ij = green_function_2d(self.nodes[i], self.nodes[j]);
            }
        }
        g
    }
}
/// Gauss–Lobatto–Legendre quadrature rule.
pub struct GaussLobattoLegendre {
    /// Polynomial degree n (number of points = n+1).
    pub degree: usize,
    /// GLL nodes on \[−1, 1\].
    pub nodes: Vec<f64>,
    /// GLL weights.
    pub weights: Vec<f64>,
}
impl GaussLobattoLegendre {
    /// Construct a GLL rule of degree `n`.
    pub fn new(degree: usize) -> Self {
        let (nodes, weights) = gll_nodes_weights(degree);
        GaussLobattoLegendre {
            degree,
            nodes,
            weights,
        }
    }
    /// Integrate `f` over \[−1, 1\] using this GLL rule.
    pub fn integrate<F: Fn(f64) -> f64>(&self, f: F) -> f64 {
        self.nodes
            .iter()
            .zip(self.weights.iter())
            .map(|(&x, &w)| w * f(x))
            .sum()
    }
    /// Number of quadrature points.
    pub fn n_points(&self) -> usize {
        self.nodes.len()
    }
}
/// Barycentric spectral interpolation at GLL nodes.
pub struct SpectralInterpolation {
    /// GLL nodes.
    pub nodes: Vec<f64>,
    /// Barycentric weights.
    pub bary_weights: Vec<f64>,
}
impl SpectralInterpolation {
    /// Construct from GLL nodes.
    pub fn new(nodes: Vec<f64>) -> Self {
        let n = nodes.len();
        let mut bary_weights = vec![1.0f64; n];
        for j in 0..n {
            for k in 0..n {
                if k == j {
                    continue;
                }
                let diff = nodes[j] - nodes[k];
                if diff.abs() > 1e-300 {
                    bary_weights[j] /= diff;
                }
            }
        }
        SpectralInterpolation {
            nodes,
            bary_weights,
        }
    }
    /// Interpolate `values` at point `x`.
    pub fn interpolate(&self, values: &[f64], x: f64) -> f64 {
        chebyshev_interp(&self.nodes, values, x)
    }
    /// Interpolate at multiple points.
    pub fn interpolate_many(&self, values: &[f64], points: &[f64]) -> Vec<f64> {
        points
            .iter()
            .map(|&x| self.interpolate(values, x))
            .collect()
    }
}
/// Global-to-local DOF mapping for a 1-D SEM mesh.
///
/// Maps global node indices to local element-node pairs, handling shared
/// interface nodes with C^0 continuity.
pub struct GlobalLocalMapping {
    /// Local-to-global connectivity: `conn[e][i]` = global DOF index.
    pub conn: Vec<Vec<usize>>,
    /// Total number of global DOFs.
    pub n_global: usize,
}
impl GlobalLocalMapping {
    /// Build the connectivity for an hp mesh with C^0 continuity.
    pub fn from_hp_mesh(mesh: &HpMesh1D) -> Self {
        let n_elem = mesh.n_elem();
        let mut conn: Vec<Vec<usize>> = Vec::with_capacity(n_elem);
        let mut global_idx = 0usize;
        for e in 0..n_elem {
            let n_local = mesh.degrees[e] + 1;
            let mut local_conn = Vec::with_capacity(n_local);
            if e == 0 {
                for i in 0..n_local {
                    local_conn.push(global_idx + i);
                }
                global_idx += n_local;
            } else {
                local_conn.push(
                    *conn[e - 1]
                        .last()
                        .expect("previous element connectivity is non-empty"),
                );
                for i in 1..n_local {
                    local_conn.push(global_idx + i - 1);
                }
                global_idx += n_local - 1;
            }
            conn.push(local_conn);
        }
        GlobalLocalMapping {
            conn,
            n_global: global_idx,
        }
    }
    /// Scatter local element vector `u_loc` into global vector `u_glob`.
    pub fn scatter(&self, elem: usize, u_loc: &[f64], u_glob: &mut [f64]) {
        for (i, &g) in self.conn[elem].iter().enumerate() {
            if g < u_glob.len() {
                u_glob[g] += u_loc[i];
            }
        }
    }
    /// Gather global vector into local element vector.
    pub fn gather(&self, elem: usize, u_glob: &[f64]) -> Vec<f64> {
        self.conn[elem]
            .iter()
            .map(|&g| if g < u_glob.len() { u_glob[g] } else { 0.0 })
            .collect()
    }
}
/// Modal-to-nodal and nodal-to-modal basis transformations for SEM.
///
/// In the modal basis, coefficients represent amplitudes of Legendre modes P_k(ξ).
/// In the nodal basis, values are given at GLL quadrature nodes.
pub struct ModalNodalTransform {
    /// Polynomial degree N.
    pub degree: usize,
    /// GLL nodes.
    pub nodes: Vec<f64>,
    /// GLL weights.
    pub weights: Vec<f64>,
    /// Vandermonde matrix V\[i\]\[k\] = P_k(ξ_i) / norm_k.
    pub vandermonde: Vec<Vec<f64>>,
}
impl ModalNodalTransform {
    /// Construct the modal-nodal transform for degree `n`.
    pub fn new(degree: usize) -> Self {
        let (nodes, weights) = gll_nodes_weights(degree);
        let n = nodes.len();
        let mut vandermonde = vec![vec![0.0; n]; n];
        for (v_row, &node_i) in vandermonde.iter_mut().zip(nodes.iter()) {
            for (k, v_ik) in v_row.iter_mut().enumerate() {
                let pk = legendre_pn(k, node_i);
                let norm = ((2 * k + 1) as f64 / 2.0).sqrt();
                *v_ik = pk * norm;
            }
        }
        ModalNodalTransform {
            degree,
            nodes,
            weights,
            vandermonde,
        }
    }
    /// Nodal to modal: modal coefficients = V^{-T} u_nodal.
    ///
    /// Exploits GLL quadrature exactness: â_k = Σ_i w_i P_k(ξ_i) u_i / norm_k^2 * norm_k.
    pub fn nodal_to_modal(&self, u_nodal: &[f64]) -> Vec<f64> {
        let n = self.nodes.len();
        let mut modal = vec![0.0; n];
        for (k, modal_k) in modal.iter_mut().enumerate() {
            let norm_sq = 2.0 / (2 * k + 1) as f64;
            let sum: f64 = self
                .weights
                .iter()
                .zip(self.nodes.iter())
                .zip(u_nodal.iter())
                .map(|((&w, &node), &u)| w * legendre_pn(k, node) * u)
                .sum();
            *modal_k = sum / norm_sq;
        }
        modal
    }
    /// Modal to nodal: u_nodal\[i\] = Σ_k â_k P_k(ξ_i).
    pub fn modal_to_nodal(&self, modal: &[f64]) -> Vec<f64> {
        let n = self.nodes.len();
        let mut nodal = vec![0.0; n];
        for (nodal_i, &node_i) in nodal.iter_mut().zip(self.nodes.iter()) {
            *nodal_i = modal
                .iter()
                .enumerate()
                .map(|(k, &m)| m * legendre_pn(k, node_i))
                .sum();
        }
        nodal
    }
}
/// 1-D spectral element with Gauss–Lobatto–Legendre nodes.
pub struct SpectralElement1D {
    /// Polynomial degree (number of nodes = degree + 1).
    pub degree: usize,
    /// GLL nodes on \[−1, 1\].
    pub nodes: Vec<f64>,
    /// GLL weights.
    pub weights: Vec<f64>,
    /// Derivative matrix D_{ij}.
    pub d_matrix: Vec<Vec<f64>>,
    /// Element length (Jacobian = length/2).
    pub length: f64,
}
impl SpectralElement1D {
    /// Construct a 1-D spectral element of given degree and physical length.
    pub fn new(degree: usize, length: f64) -> Self {
        let (nodes, weights) = gll_nodes_weights(degree);
        let d_matrix = spectral_derivative_matrix(&nodes);
        SpectralElement1D {
            degree,
            nodes,
            weights,
            d_matrix,
            length,
        }
    }
    /// Jacobian dx/dξ = length/2.
    pub fn jacobian(&self) -> f64 {
        self.length / 2.0
    }
    /// Mass matrix entry M_{ij} = ∫_{-1}^{1} l_i(ξ) l_j(ξ) J dξ.
    ///
    /// For GLL quadrature, the mass matrix is diagonal: M_{ii} = w_i * J.
    pub fn mass_matrix_diagonal(&self) -> Vec<f64> {
        let j = self.jacobian();
        self.weights.iter().map(|&w| w * j).collect()
    }
    /// Stiffness matrix entry K_{ij} = ∫_{-1}^{1} l'_i l'_j / J dξ.
    pub fn stiffness_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.nodes.len();
        let j = self.jacobian();
        let mut k = vec![vec![0.0; n]; n];
        for (i, k_row) in k.iter_mut().enumerate() {
            for (j_idx, k_ij) in k_row.iter_mut().enumerate() {
                let kij: f64 = self
                    .weights
                    .iter()
                    .zip(self.d_matrix.iter())
                    .map(|(&w, d_q)| w * d_q[i] * d_q[j_idx])
                    .sum();
                *k_ij = kij / j;
            }
        }
        k
    }
    /// Differentiate nodal values `u` at the GLL nodes.
    pub fn differentiate(&self, u: &[f64]) -> Vec<f64> {
        let j = self.jacobian();
        let mut du = vec![0.0; self.nodes.len()];
        for (du_i, d_row) in du.iter_mut().zip(self.d_matrix.iter()) {
            *du_i = d_row
                .iter()
                .zip(u.iter())
                .map(|(&d, &u_j)| d * u_j)
                .sum::<f64>()
                / j;
        }
        du
    }
}
