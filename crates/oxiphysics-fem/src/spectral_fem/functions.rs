//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::f64::consts::PI;

use super::types::{GlobalLocalMapping, HpMesh1D, SemEllipticSolver, SpectralElement1D};

/// Evaluates the Legendre polynomial P_n(x) and its derivative P'_n(x).
///
/// Uses the three-term recurrence relation.
pub fn legendre_evaluate(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    if n == 1 {
        return (x, 1.0);
    }
    let mut p_prev = 1.0_f64;
    let mut p_curr = x;
    for k in 2..=(n as u64) {
        let kf = k as f64;
        let p_next = ((2.0 * kf - 1.0) * x * p_curr - (kf - 1.0) * p_prev) / kf;
        p_prev = p_curr;
        p_curr = p_next;
    }
    let dp = if (x * x - 1.0).abs() < 1e-14 {
        let sign = if x > 0.0 || n.is_multiple_of(2) {
            1.0
        } else {
            -1.0
        };
        sign * (n as f64) * (n as f64 + 1.0) / 2.0
    } else {
        let nf = n as f64;
        nf * (x * p_curr - p_prev) / (x * x - 1.0)
    };
    (p_curr, dp)
}
/// Evaluates the n-th Legendre polynomial P_n(x).
pub fn legendre_pn(n: usize, x: f64) -> f64 {
    legendre_evaluate(n, x).0
}
/// Computes the n zeros of P_n(x) using Newton's method starting from Chebyshev nodes.
///
/// Returns the zeros sorted in ascending order on \[−1, 1\].
pub fn legendre_zeros(n: usize) -> Vec<f64> {
    if n == 0 {
        return vec![];
    }
    let mut zeros = Vec::with_capacity(n);
    let half = n.div_ceil(2);
    for k in 1..=half {
        let mut x = -((2 * k - 1) as f64 * PI / (2 * n) as f64).cos();
        for _ in 0..50 {
            let (p, dp) = legendre_evaluate(n, x);
            if dp.abs() < 1e-300 {
                break;
            }
            let dx = -p / dp;
            x += dx;
            if dx.abs() < 1e-14 {
                break;
            }
        }
        zeros.push(x);
    }
    let mut all_zeros = Vec::with_capacity(n);
    if n % 2 == 1 {
        all_zeros.push(0.0);
    }
    for &z in zeros.iter().rev() {
        if z.abs() > 1e-14 {
            all_zeros.insert(0, -z);
        }
    }
    for &z in zeros.iter() {
        if z.abs() > 1e-14 {
            all_zeros.push(z);
        }
    }

    all_zeros.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    all_zeros.truncate(n);
    all_zeros
}
/// Compute GLL nodes and weights for a given polynomial degree `n` (producing n+1 points).
///
/// The GLL nodes are the zeros of (1-x²) P'_n(x) and include ±1.
/// Returns `(nodes, weights)`.
pub fn gll_nodes_weights(n: usize) -> (Vec<f64>, Vec<f64>) {
    let np1 = n + 1;
    if n == 0 {
        return (vec![0.0], vec![2.0]);
    }
    if n == 1 {
        return (vec![-1.0, 1.0], vec![1.0, 1.0]);
    }
    let mut nodes = vec![0.0; np1];
    nodes[0] = -1.0;
    nodes[n] = 1.0;
    let half_interior = (n - 1).div_ceil(2);
    let n_interior = n - 1;
    let mut interior = Vec::with_capacity(n_interior);
    for k in 1..=half_interior {
        let mut x = -((k as f64 * PI) / n as f64).cos();
        for _ in 0..50 {
            let (p, dp) = legendre_evaluate(n, x);
            if dp.abs() < 1e-300 {
                break;
            }
            let d2p = if (1.0 - x * x).abs() < 1e-14 {
                0.0
            } else {
                (2.0 * x * dp - n as f64 * (n as f64 + 1.0) * p) / (1.0 - x * x)
            };
            if d2p.abs() < 1e-300 {
                break;
            }
            let dx = -dp / d2p;
            x += dx;
            if dx.abs() < 1e-14 {
                break;
            }
        }
        interior.push(x);
    }
    let mut all_interior = Vec::with_capacity(n_interior);
    if n_interior % 2 == 1 {
        all_interior.push(0.0);
    }
    for &z in interior.iter().rev() {
        if z.abs() > 1e-14 {
            all_interior.insert(0, -z);
        }
    }
    for &z in interior.iter() {
        if z.abs() > 1e-14 {
            all_interior.push(z);
        }
    }
    all_interior.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    all_interior.truncate(n_interior);
    for (i, &x) in all_interior.iter().enumerate() {
        nodes[i + 1] = x;
    }
    let mut weights = vec![0.0; np1];
    let nf = n as f64;
    let denom = nf * (nf + 1.0);
    for i in 0..np1 {
        let p = legendre_pn(n, nodes[i]);
        weights[i] = 2.0 / (denom * p * p);
    }
    (nodes, weights)
}
/// Compute the spectral derivative matrix D_{ij} = l'_j(x_i) at GLL nodes.
///
/// The entries are the derivatives of the Lagrange basis polynomials at the
/// GLL nodes.
pub fn spectral_derivative_matrix(nodes: &[f64]) -> Vec<Vec<f64>> {
    let n = nodes.len();
    let mut d = vec![vec![0.0; n]; n];
    let deg = n - 1;
    let p_nodes: Vec<f64> = nodes.iter().map(|&x| legendre_pn(deg, x)).collect();
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let diff = nodes[i] - nodes[j];
            if diff.abs() < 1e-14 {
                continue;
            }
            d[i][j] = (p_nodes[i] / p_nodes[j]) / diff;
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
    d
}
/// Barycentric interpolation at GLL nodes.
///
/// Given nodal values `f_j = f(x_j)` at GLL nodes, interpolates to any point `x`.
pub fn chebyshev_interp(nodes: &[f64], values: &[f64], x: f64) -> f64 {
    for (&xi, &fi) in nodes.iter().zip(values.iter()) {
        if (x - xi).abs() < 1e-14 {
            return fi;
        }
    }
    let n = nodes.len();
    let mut weights = vec![1.0f64; n];
    for j in 0..n {
        for k in 0..n {
            if k == j {
                continue;
            }
            weights[j] *= nodes[j] - nodes[k];
        }
        if weights[j].abs() > 1e-300 {
            weights[j] = 1.0 / weights[j];
        }
    }
    let num: f64 = (0..n)
        .map(|j| weights[j] * values[j] / (x - nodes[j]))
        .sum();
    let den: f64 = (0..n).map(|j| weights[j] / (x - nodes[j])).sum();
    if den.abs() < 1e-300 {
        return values[0];
    }
    num / den
}
/// Boundary element method: 2-D Green's function for Laplace equation.
///
/// G(x, y) = -1/(2π) * ln|x - y|
pub fn green_function_2d(x: [f64; 2], y: [f64; 2]) -> f64 {
    let r = ((x[0] - y[0]).powi(2) + (x[1] - y[1]).powi(2)).sqrt();
    if r < 1e-300 {
        return 0.0;
    }
    -0.5 / PI * r.ln()
}
/// Normal derivative of 2-D Laplace Green's function.
///
/// ∂G/∂n_y = -1/(2π) * (x-y)·n / |x-y|²
pub fn green_normal_derivative_2d(x: [f64; 2], y: [f64; 2], n_y: [f64; 2]) -> f64 {
    let dx = x[0] - y[0];
    let dy = x[1] - y[1];
    let r2 = dx * dx + dy * dy;
    if r2 < 1e-300 {
        return 0.0;
    }
    let dot = dx * n_y[0] + dy * n_y[1];
    -dot / (2.0 * PI * r2)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::spectral_fem::*;
    pub(super) const EPS: f64 = 1e-10;
    #[test]
    fn test_legendre_p0() {
        let (p, dp) = legendre_evaluate(0, 0.5);
        assert!((p - 1.0).abs() < EPS);
        assert!((dp - 0.0).abs() < EPS);
    }
    #[test]
    fn test_legendre_p1() {
        let x = 0.7;
        let (p, dp) = legendre_evaluate(1, x);
        assert!((p - x).abs() < EPS);
        assert!((dp - 1.0).abs() < EPS);
    }
    #[test]
    fn test_legendre_p2() {
        let x = 0.5;
        let (p, _) = legendre_evaluate(2, x);
        let expected = (3.0 * x * x - 1.0) / 2.0;
        assert!((p - expected).abs() < EPS);
    }
    #[test]
    fn test_legendre_p3_zeros() {
        let lp = LegendrePolynomial::new(3);
        let zeros = lp.zeros();
        assert_eq!(zeros.len(), 3);
        for &z in &zeros {
            let val = lp.evaluate(z);
            assert!(val.abs() < 1e-8, "P_3({})={}", z, val);
        }
    }
    #[test]
    fn test_gll_n1_nodes() {
        let (nodes, _weights) = gll_nodes_weights(1);
        assert_eq!(nodes.len(), 2);
        assert!((nodes[0] - (-1.0)).abs() < EPS);
        assert!((nodes[1] - 1.0).abs() < EPS);
    }
    #[test]
    fn test_gll_n2_weights_sum() {
        let (_, weights) = gll_nodes_weights(2);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-8);
    }
    #[test]
    fn test_gll_n3_weights_sum() {
        let (_, weights) = gll_nodes_weights(3);
        let sum: f64 = weights.iter().sum();
        assert!((sum - 2.0).abs() < 1e-8);
    }
    #[test]
    fn test_gll_integrates_constants() {
        let gll = GaussLobattoLegendre::new(4);
        let result = gll.integrate(|_x| 1.0);
        assert!((result - 2.0).abs() < 1e-8);
    }
    #[test]
    fn test_gll_integrates_linear() {
        let gll = GaussLobattoLegendre::new(4);
        let result = gll.integrate(|x| x);
        assert!(result.abs() < 1e-12);
    }
    #[test]
    fn test_gll_integrates_quadratic() {
        let gll = GaussLobattoLegendre::new(4);
        let result = gll.integrate(|x| x * x);
        assert!((result - 2.0 / 3.0).abs() < 1e-10, "result={}", result);
    }
    #[test]
    fn test_spectral_derivative_matrix_size() {
        let (nodes, _) = gll_nodes_weights(4);
        let d = spectral_derivative_matrix(&nodes);
        assert_eq!(d.len(), 5);
        assert_eq!(d[0].len(), 5);
    }
    #[test]
    fn test_spectral_derivative_matrix_differentiates_linear() {
        let deg = 4;
        let (nodes, _) = gll_nodes_weights(deg);
        let d = spectral_derivative_matrix(&nodes);
        let n = nodes.len();
        let mut du = vec![0.0; n];
        for i in 0..n {
            for j in 0..n {
                du[i] += d[i][j] * nodes[j];
            }
        }
        for val in &du {
            assert!((val - 1.0).abs() < 1e-8, "du={}", val);
        }
    }
    #[test]
    fn test_spectral_derivative_rowsum_zero() {
        let deg = 3;
        let (nodes, _) = gll_nodes_weights(deg);
        let d = spectral_derivative_matrix(&nodes);
        for row in &d {
            let sum: f64 = row.iter().sum();
            assert!(sum.abs() < 1e-10, "row_sum={}", sum);
        }
    }
    #[test]
    fn test_spectral_interpolation_at_node() {
        let gll = GaussLobattoLegendre::new(3);
        let interp = SpectralInterpolation::new(gll.nodes.clone());
        let values: Vec<f64> = gll.nodes.iter().map(|&x| x * x).collect();
        let x = gll.nodes[2];
        let result = interp.interpolate(&values, x);
        assert!((result - x * x).abs() < 1e-10, "result={}", result);
    }
    #[test]
    fn test_spectral_interpolation_midpoint() {
        let gll = GaussLobattoLegendre::new(4);
        let interp = SpectralInterpolation::new(gll.nodes.clone());
        let values: Vec<f64> = gll.nodes.iter().map(|&x| x * x).collect();
        let result = interp.interpolate(&values, 0.0);
        assert!(result.abs() < 1e-8, "result={}", result);
    }
    #[test]
    fn test_spectral_element_1d_mass_positive() {
        let elem = SpectralElement1D::new(4, 1.0);
        let m = elem.mass_matrix_diagonal();
        for &mi in &m {
            assert!(mi > 0.0);
        }
    }
    #[test]
    fn test_spectral_element_1d_mass_sum() {
        let length = 2.5;
        let elem = SpectralElement1D::new(4, length);
        let m = elem.mass_matrix_diagonal();
        let sum: f64 = m.iter().sum();
        assert!(
            (sum - length).abs() < 1e-10,
            "sum={}, length={}",
            sum,
            length
        );
    }
    #[test]
    fn test_spectral_element_1d_stiffness_rowsum() {
        let elem = SpectralElement1D::new(3, 1.0);
        let k = elem.stiffness_matrix();
        for row in &k {
            let sum: f64 = row.iter().sum();
            assert!(sum.abs() < 1e-8, "row_sum={}", sum);
        }
    }
    #[test]
    fn test_spectral_element_2d_mass_positive() {
        let elem = SpectralElement2D::new(3, [1.0, 1.0]);
        let m = elem.mass_matrix_diagonal();
        for &mi in &m {
            assert!(mi > 0.0);
        }
    }
    #[test]
    fn test_spectral_element_2d_mass_sum() {
        let lx = 2.0;
        let ly = 3.0;
        let elem = SpectralElement2D::new(3, [lx, ly]);
        let m = elem.mass_matrix_diagonal();
        let sum: f64 = m.iter().sum();
        assert!((sum - lx * ly).abs() < 1e-8, "sum={}", sum);
    }
    #[test]
    fn test_spectral_mass_matrix_inverse() {
        let sm = SpectralMassMatrix::from_1d_element(3, 1.0);
        let inv = sm.inverse();
        for (&m, &mi) in sm.diagonal.iter().zip(inv.iter()) {
            assert!((m * mi - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn test_spectral_stiffness_1d_rowsum() {
        let ss = SpectralStiffness::from_1d_element(3, 1.0);
        for i in 0..ss.matrix.len() {
            assert!(ss.row_sum(i).abs() < 1e-8);
        }
    }
    #[test]
    fn test_spectral_convection_constant_field() {
        let conv = SpectralConvection::new(4, 1.0, 1.0);
        let n = conv.element.nodes.len();
        let u = vec![1.0; n];
        let rhs = conv.flux_divergence(&u);
        for r in rhs {
            assert!(r.abs() < 1e-8, "r={}", r);
        }
    }
    #[test]
    fn test_spectral_helmholtz_system_matrix_size() {
        let helm = SpectralHelmholtz::new(4, 1.0, 1.0);
        let a = helm.system_matrix();
        assert_eq!(a.len(), 5);
    }
    #[test]
    fn test_green_function_symmetry() {
        let x = [1.0, 0.0];
        let y = [0.0, 1.0];
        let g1 = green_function_2d(x, y);
        let g2 = green_function_2d(y, x);
        assert!((g1 - g2).abs() < EPS);
    }
    #[test]
    fn test_green_function_nonzero() {
        let x = [2.0, 0.0];
        let y = [0.0, 0.0];
        let g = green_function_2d(x, y);
        assert!(g.abs() > 0.01);
    }
    #[test]
    fn test_bem_h_matrix_size() {
        let nodes = vec![[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]];
        let normals = vec![[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]];
        let bem = SpectralBoundaryIntegral::new(nodes, normals);
        let h = bem.h_matrix();
        assert_eq!(h.len(), 3);
        assert_eq!(h[0].len(), 3);
    }
    #[test]
    fn test_bem_g_matrix_diagonal_zero() {
        let nodes = vec![[1.0, 0.0], [0.0, 1.0]];
        let normals = vec![[1.0, 0.0], [0.0, 1.0]];
        let bem = SpectralBoundaryIntegral::new(nodes, normals);
        let g = bem.g_matrix();
        assert_eq!(g[0][0], 0.0);
        assert_eq!(g[1][1], 0.0);
    }
    #[test]
    fn test_chebyshev_interp_at_nodes() {
        let nodes = vec![-1.0, 0.0, 1.0];
        let values = vec![1.0, 0.0, 1.0];
        let result = chebyshev_interp(&nodes, &values, 0.0);
        assert!(result.abs() < 1e-10);
    }
    #[test]
    fn test_legendre_derivative() {
        let lp = LegendrePolynomial::new(2);
        let dp = lp.derivative(0.5);
        assert!((dp - 1.5).abs() < EPS, "dp={}", dp);
    }
    #[test]
    fn test_spectral_element_1d_differentiate_linear() {
        let deg = 4;
        let elem = SpectralElement1D::new(deg, 2.0);
        let u: Vec<f64> = elem.nodes.to_vec();
        let du = elem.differentiate(&u);
        for &d in &du {
            assert!((d - 1.0).abs() < 1e-8, "d={}", d);
        }
    }
}
/// Assemble the global stiffness matrix for a uniform 1-D SEM mesh.
///
/// Returns the (N_dof × N_dof) dense stiffness matrix where continuity
/// is enforced by adding shared node contributions.
pub fn assemble_sem_stiffness_1d(n_elem: usize, total_length: f64, degree: usize) -> Vec<Vec<f64>> {
    let mesh = HpMesh1D::uniform(n_elem, total_length, degree);
    let solver = SemEllipticSolver::new(mesh);
    solver.k_global
}
/// Assemble the global mass vector (diagonal) for a uniform 1-D SEM mesh.
pub fn assemble_sem_mass_diagonal_1d(n_elem: usize, total_length: f64, degree: usize) -> Vec<f64> {
    let mesh = HpMesh1D::uniform(n_elem, total_length, degree);
    let mapping = GlobalLocalMapping::from_hp_mesh(&mesh);
    let ng = mapping.n_global;
    let mut m_glob = vec![0.0; ng];
    for e in 0..n_elem {
        let elem = SpectralElement1D::new(mesh.degrees[e], mesh.lengths[e]);
        let m_loc = elem.mass_matrix_diagonal();
        for (i, &m) in m_loc.iter().enumerate() {
            let gi = mapping.conn[e][i];
            m_glob[gi] += m;
        }
    }
    m_glob
}
#[cfg(test)]
mod tests_part2 {
    use super::*;
    use crate::spectral_fem::*;
    #[test]
    fn test_hp_mesh_uniform_n_dofs() {
        let mesh = HpMesh1D::uniform(3, 3.0, 2);
        assert_eq!(mesh.n_dofs(), 7);
    }
    #[test]
    fn test_hp_mesh_p_refine() {
        let mut mesh = HpMesh1D::uniform(3, 3.0, 2);
        mesh.p_refine(1, 2);
        assert_eq!(mesh.degrees[1], 4);
    }
    #[test]
    fn test_hp_mesh_h_refine_count() {
        let mut mesh = HpMesh1D::uniform(2, 2.0, 3);
        assert_eq!(mesh.n_elem(), 2);
        mesh.h_refine(0);
        assert_eq!(mesh.n_elem(), 3);
    }
    #[test]
    fn test_hp_mesh_h_refine_lengths() {
        let mut mesh = HpMesh1D::uniform(1, 2.0, 3);
        mesh.h_refine(0);
        assert!((mesh.lengths[0] - 1.0).abs() < 1e-14);
        assert!((mesh.lengths[1] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn test_spectral_convergence_l2_error_zero_for_polynomial() {
        let sc = SpectralConvergence::new(4);
        let err = sc.l2_interp_error(|x| x * x * x * x);
        assert!(err < 1e-10, "err={}", err);
    }
    #[test]
    fn test_spectral_convergence_error_decreases_with_degree() {
        let sc_lo = SpectralConvergence::new(4);
        let sc_hi = SpectralConvergence::new(8);
        let f = |x: f64| (std::f64::consts::PI * x).sin();
        let e_lo = sc_lo.l2_interp_error(f);
        let e_hi = sc_hi.l2_interp_error(f);
        assert!(e_hi <= e_lo, "e_hi={}, e_lo={}", e_hi, e_lo);
    }
    #[test]
    fn test_chebyshev_spectral_roundtrip() {
        let cs = ChebyshevSpectral::new(5);
        let vals: Vec<f64> = cs.nodes.iter().map(|&x| x * x + 1.0).collect();
        let coeffs = cs.to_coefficients(&vals);
        let reconstructed = cs.to_nodal(&coeffs);
        for (v, r) in vals.iter().zip(reconstructed.iter()) {
            assert!((v - r).abs() < 1e-8, "v={}, r={}", v, r);
        }
    }
    #[test]
    fn test_chebyshev_dealias_zeroes_high_modes() {
        let cs = ChebyshevSpectral::new(8);
        let vals: Vec<f64> = cs.nodes.iter().map(|&x| x.sin()).collect();
        let mut coeffs = cs.to_coefficients(&vals);
        cs.dealias_filter(&mut coeffs, 2.0 / 3.0);
        for (k, &c) in coeffs.iter().enumerate().skip(6) {
            assert!(c.abs() < 1e-14, "coeff[{}]={}", k, c);
        }
    }
    #[test]
    fn test_modal_nodal_roundtrip() {
        let mnt = ModalNodalTransform::new(4);
        let u_nodal: Vec<f64> = mnt.nodes.iter().map(|&x| x * x).collect();
        let modal = mnt.nodal_to_modal(&u_nodal);
        let reconstructed = mnt.modal_to_nodal(&modal);
        for (a, b) in u_nodal.iter().zip(reconstructed.iter()) {
            assert!((a - b).abs() < 1e-8, "a={}, b={}", a, b);
        }
    }
    #[test]
    fn test_modal_nodal_constant_field() {
        let mnt = ModalNodalTransform::new(3);
        let u_nodal = vec![1.0; mnt.nodes.len()];
        let modal = mnt.nodal_to_modal(&u_nodal);
        assert!((modal[0] - 1.0).abs() < 1e-8, "modal[0]={}", modal[0]);
        for (k, &m) in modal.iter().enumerate().skip(1) {
            assert!(m.abs() < 1e-8, "modal[{}]={}", k, m);
        }
    }
    #[test]
    fn test_global_local_mapping_n_global() {
        let mesh = HpMesh1D::uniform(3, 3.0, 2);
        let map = GlobalLocalMapping::from_hp_mesh(&mesh);
        assert_eq!(map.n_global, 7);
    }
    #[test]
    fn test_global_local_mapping_shared_node() {
        let mesh = HpMesh1D::uniform(2, 2.0, 2);
        let map = GlobalLocalMapping::from_hp_mesh(&mesh);
        let last_e0 = *map.conn[0].last().unwrap();
        let first_e1 = map.conn[1][0];
        assert_eq!(last_e0, first_e1);
    }
    #[test]
    fn test_global_local_mapping_scatter_gather_roundtrip() {
        let mesh = HpMesh1D::uniform(2, 2.0, 2);
        let map = GlobalLocalMapping::from_hp_mesh(&mesh);
        let u_glob: Vec<f64> = (0..map.n_global).map(|i| i as f64).collect();
        let u_loc = map.gather(0, &u_glob);
        let mut u_glob_out = vec![0.0; map.n_global];
        map.scatter(0, &u_loc, &mut u_glob_out);
        for &g in &map.conn[0] {
            assert!((u_glob_out[g] - u_glob[g]).abs() < 1e-14);
        }
    }
    #[test]
    fn test_sem_elliptic_solver_global_stiffness_positive_diagonal() {
        let mesh = HpMesh1D::uniform(2, 2.0, 2);
        let solver = SemEllipticSolver::new(mesh);
        let ng = solver.mapping.n_global;
        for i in 0..ng {
            assert!(
                solver.k_global[i][i] > 0.0,
                "diag[{}]={}",
                i,
                solver.k_global[i][i]
            );
        }
    }
    #[test]
    fn test_sem_elliptic_solver_solve_trivial_f_zero() {
        let mesh = HpMesh1D::uniform(3, 1.0, 3);
        let mut solver = SemEllipticSolver::new(mesh);
        solver.load_rhs(|_x| 0.0);
        let u = solver.solve();
        for &ui in &u {
            assert!(ui.abs() < 1e-10, "u={}", ui);
        }
    }
    #[test]
    fn test_assemble_sem_stiffness_1d_size() {
        let k = assemble_sem_stiffness_1d(3, 3.0, 2);
        assert_eq!(k.len(), 7);
        assert_eq!(k[0].len(), 7);
    }
    #[test]
    fn test_assemble_sem_mass_diagonal_sum() {
        let total = 3.0;
        let m = assemble_sem_mass_diagonal_1d(3, total, 3);
        let sum: f64 = m.iter().sum();
        assert!((sum - total).abs() < 1e-8, "sum={}", sum);
    }
    #[test]
    fn test_chebyshev_diff_matrix_constant() {
        let cdn = ChebyshevDiffMatrix::new(4);
        let u = vec![1.0; cdn.n_nodes];
        let du = cdn.differentiate(&u);
        for val in du {
            assert!(val.abs() < 1e-8, "du={}", val);
        }
    }
    #[test]
    fn test_chebyshev_diff_matrix_row_sum_zero() {
        let cdn = ChebyshevDiffMatrix::new(5);
        for i in 0..cdn.n_nodes {
            let s = cdn.row_sum(i);
            assert!(s.abs() < 1e-8, "row {} sum={}", i, s);
        }
    }
    #[test]
    fn test_chebyshev_diff_matrix_differentiates_x() {
        let cdn = ChebyshevDiffMatrix::new(4);
        let u: Vec<f64> = cdn.nodes.clone();
        let du = cdn.differentiate(&u);
        for val in &du {
            assert!((val - 1.0).abs() < 1e-8, "du={}", val);
        }
    }
    #[test]
    fn test_chebyshev_spectral_derivative_of_x() {
        let cs = ChebyshevSpectral::new(4);
        let u: Vec<f64> = cs.nodes.clone();
        let du = cs.differentiate(&u);
        let du_len = du.len();
        for (i, &dv) in du.iter().enumerate().skip(1).take(du_len.saturating_sub(2)) {
            assert!((dv - 1.0).abs() < 1e-6, "du[{}]={}", i, dv);
        }
    }
    #[test]
    fn test_spectral_convergence_rate_negative() {
        let f = |x: f64| x.powi(6) + x.powi(5) + x.powi(4);
        let sc_lo = SpectralConvergence::new(3);
        let sc_hi = SpectralConvergence::new(6);
        let e_lo = sc_lo.l2_interp_error(f);
        let e_hi = sc_hi.l2_interp_error(f);
        assert!(e_hi < e_lo, "e_hi={} e_lo={}", e_hi, e_lo);
    }
    #[test]
    fn test_gll_nodes_include_endpoints() {
        for deg in [2, 3, 4, 5, 6] {
            let (nodes, _) = gll_nodes_weights(deg);
            assert!((nodes[0] - (-1.0)).abs() < 1e-12, "deg={}", deg);
            assert!((nodes[deg] - 1.0).abs() < 1e-12, "deg={}", deg);
        }
    }
    #[test]
    fn test_gll_weights_sum_to_2_for_various_degrees() {
        for deg in 1..=7 {
            let (_, weights) = gll_nodes_weights(deg);
            let sum: f64 = weights.iter().sum();
            assert!((sum - 2.0).abs() < 1e-8, "deg={} sum={}", deg, sum);
        }
    }
    #[test]
    fn test_spectral_convergence_exact_for_degree_polynomial() {
        for p in 1..=4usize {
            let sc = SpectralConvergence::new(p);
            let err = sc.l2_interp_error(|x| x.powi(p as i32));
            assert!(err < 1e-8, "p={} err={}", p, err);
        }
    }
}
