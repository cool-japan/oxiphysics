// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermal FEM for 3-D tetrahedral elements.
//!
//! Provides steady-state and transient (explicit Euler) thermal analysis with
//! Dirichlet boundary conditions, heat-flux computation, and thermo-elastic
//! strain coupling.

// ─────────────────────────────────────────────────────────────────────────────
// § 1  Data structures
// ─────────────────────────────────────────────────────────────────────────────

/// A node in a thermal finite-element mesh.
#[derive(Debug, Clone)]
pub struct ThermalNode {
    /// Node coordinates \[x, y, z\] \[m\].
    pub position: [f64; 3],
    /// Current temperature \[K or °C\].
    pub temperature: f64,
    /// Whether a Dirichlet (prescribed) temperature is applied.
    pub is_dirichlet: bool,
    /// Prescribed temperature value (used when `is_dirichlet` is true).
    pub prescribed_temp: f64,
}

impl ThermalNode {
    /// Create a free node at `position` with initial temperature.
    pub fn free(position: [f64; 3], temperature: f64) -> Self {
        Self {
            position,
            temperature,
            is_dirichlet: false,
            prescribed_temp: 0.0,
        }
    }

    /// Create a Dirichlet node with a fixed temperature.
    pub fn dirichlet(position: [f64; 3], prescribed_temp: f64) -> Self {
        Self {
            position,
            temperature: prescribed_temp,
            is_dirichlet: true,
            prescribed_temp,
        }
    }
}

/// A linear tetrahedral thermal element with anisotropic conductivity.
#[derive(Debug, Clone)]
pub struct ThermalElement {
    /// Global node indices (4-node linear tetrahedron).
    pub nodes: [usize; 4],
    /// Thermal conductivity components \[kx, ky, kz\] \[W/(m·K)\].
    pub conductivity: [f64; 3],
}

impl ThermalElement {
    /// Create a new `ThermalElement`.
    pub fn new(nodes: [usize; 4], conductivity: [f64; 3]) -> Self {
        Self {
            nodes,
            conductivity,
        }
    }

    /// Isotropic convenience constructor.
    pub fn isotropic(nodes: [usize; 4], k: f64) -> Self {
        Self {
            nodes,
            conductivity: [k, k, k],
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 2  Fourier heat flux
// ─────────────────────────────────────────────────────────────────────────────

/// Fourier heat flux in one direction: `q = -k dT/dx`.
///
/// # Arguments
/// * `k`      – thermal conductivity \[W/(m·K)\]
/// * `grad_t` – temperature gradient dT/dx \[K/m\]
pub fn fourier_heat_flux(k: f64, grad_t: f64) -> f64 {
    -k * grad_t
}

// ─────────────────────────────────────────────────────────────────────────────
// § 3  Element stiffness (conductivity) matrix
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the 4×4 element thermal conductivity matrix for a linear tetrahedron.
///
/// Uses the constant-gradient (linear) tet formulation:
/// `K_e = V B^T D B` where B encodes the shape-function gradients and D is
/// the anisotropic conductivity tensor (diagonal \[kx, ky, kz\]).
///
/// Shape-function gradients are computed via the volume coordinates of the
/// four vertices.
///
/// # Arguments
/// * `nodes` – four corner positions \[\[x,y,z\\]; 4]
/// * `k`     – conductivity components \[kx, ky, kz\]
pub fn tet_thermal_conductivity_matrix(nodes: &[[f64; 3]; 4], k: [f64; 3]) -> [[f64; 4]; 4] {
    // Jacobian matrix [3×3]: columns are (p1-p0), (p2-p0), (p3-p0)
    let p0 = nodes[0];
    let p1 = nodes[1];
    let p2 = nodes[2];
    let p3 = nodes[3];

    let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let c = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];

    // det(J) = a · (b × c)
    let bxc = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    let det_j = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];
    let vol = det_j.abs() / 6.0;

    if vol < 1e-30 {
        return [[0.0; 4]; 4];
    }

    // Inverse of J (3×3, column-major from a,b,c):
    // J = [a | b | c] as column vectors → J^{-T} rows give grad(N)
    let inv_det = 1.0 / det_j;

    // Cofactor matrix rows (transposed inverse scaled by det)
    // Row i of J^{-1} = cofactor col / det
    let j_inv = [
        [
            (b[1] * c[2] - b[2] * c[1]) * inv_det,
            -(a[1] * c[2] - a[2] * c[1]) * inv_det,
            (a[1] * b[2] - a[2] * b[1]) * inv_det,
        ],
        [
            -(b[0] * c[2] - b[2] * c[0]) * inv_det,
            (a[0] * c[2] - a[2] * c[0]) * inv_det,
            -(a[0] * b[2] - a[2] * b[0]) * inv_det,
        ],
        [
            (b[0] * c[1] - b[1] * c[0]) * inv_det,
            -(a[0] * c[1] - a[1] * c[0]) * inv_det,
            (a[0] * b[1] - a[1] * b[0]) * inv_det,
        ],
    ];

    // Shape-function gradients in global coords (B matrix, 3×4)
    // grad N0 = -sum(grad Ni for i=1..3)
    // grad N1 = J^{-T} e1, grad N2 = J^{-T} e2, grad N3 = J^{-T} e3
    // where ei is the i-th column of J^{-T} = rows of J^{-1}
    let dn: [[f64; 3]; 4] = [
        [
            -(j_inv[0][0] + j_inv[1][0] + j_inv[2][0]),
            -(j_inv[0][1] + j_inv[1][1] + j_inv[2][1]),
            -(j_inv[0][2] + j_inv[1][2] + j_inv[2][2]),
        ],
        [j_inv[0][0], j_inv[0][1], j_inv[0][2]],
        [j_inv[1][0], j_inv[1][1], j_inv[1][2]],
        [j_inv[2][0], j_inv[2][1], j_inv[2][2]],
    ];

    // K_e = V * B^T D B  where D = diag(kx,ky,kz)
    let mut ke = [[0.0_f64; 4]; 4];
    for i in 0..4 {
        for j in 0..4 {
            ke[i][j] = vol
                * (k[0] * dn[i][0] * dn[j][0]
                    + k[1] * dn[i][1] * dn[j][1]
                    + k[2] * dn[i][2] * dn[j][2]);
        }
    }
    ke
}

// ─────────────────────────────────────────────────────────────────────────────
// § 4  Global assembly
// ─────────────────────────────────────────────────────────────────────────────

/// Assemble the global thermal stiffness matrix K (n_nodes × n_nodes).
///
/// Each element contributes via scatter-add from its 4×4 local matrix.
pub fn assemble_thermal_stiffness(
    elements: &[ThermalElement],
    nodes: &[ThermalNode],
) -> Vec<Vec<f64>> {
    let n = nodes.len();
    let mut k_global = vec![vec![0.0_f64; n]; n];

    for elem in elements {
        let pos: [[f64; 3]; 4] = [
            nodes[elem.nodes[0]].position,
            nodes[elem.nodes[1]].position,
            nodes[elem.nodes[2]].position,
            nodes[elem.nodes[3]].position,
        ];
        let ke = tet_thermal_conductivity_matrix(&pos, elem.conductivity);
        for a in 0..4 {
            for b in 0..4 {
                k_global[elem.nodes[a]][elem.nodes[b]] += ke[a][b];
            }
        }
    }
    k_global
}

// ─────────────────────────────────────────────────────────────────────────────
// § 5  Dirichlet boundary conditions
// ─────────────────────────────────────────────────────────────────────────────

/// Apply Dirichlet (prescribed-temperature) boundary conditions in place.
///
/// For each Dirichlet node i: zero row i and column i, set K\[i\]\[i\]=1, f\[i\]=T_prescribed.
pub fn apply_thermal_dirichlet(k: &mut [Vec<f64>], f: &mut [f64], nodes: &[ThermalNode]) {
    let n = nodes.len();
    for (i, node) in nodes.iter().enumerate() {
        if node.is_dirichlet {
            let t_p = node.prescribed_temp;
            // Modify RHS to account for off-diagonal terms already zeroed
            for j in 0..n {
                if j != i {
                    f[j] -= k[j][i] * t_p;
                    k[j][i] = 0.0;
                    k[i][j] = 0.0;
                }
            }
            k[i][i] = 1.0;
            f[i] = t_p;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 6  Steady-state solver (Gaussian elimination)
// ─────────────────────────────────────────────────────────────────────────────

/// Solve the linear system K T = f using Gaussian elimination with
/// partial pivoting.
///
/// Returns the nodal temperature vector.
pub fn solve_thermal_steady(k: &[Vec<f64>], f: &[f64]) -> Vec<f64> {
    let n = f.len();
    // Build augmented matrix [K | f]
    let mut aug: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut row = k[i].clone();
            row.push(f[i]);
            row
        })
        .collect();

    for col in 0..n {
        // Partial pivot
        let mut max_val = aug[col][col].abs();
        let mut max_row = col;
        for (row, aug_row) in aug.iter().enumerate().skip(col + 1) {
            if aug_row[col].abs() > max_val {
                max_val = aug_row[col].abs();
                max_row = row;
            }
        }
        aug.swap(col, max_row);

        let pivot = aug[col][col];
        if pivot.abs() < 1e-30 {
            continue;
        }
        let col_slice: Vec<f64> = aug[col][col..=n].to_vec();
        for aug_row in aug[col + 1..].iter_mut() {
            let factor = aug_row[col] / pivot;
            for (off, &cv) in col_slice.iter().enumerate() {
                aug_row[col + off] -= cv * factor;
            }
        }
    }

    // Back substitution
    let mut x = vec![0.0_f64; n];
    for i in (0..n).rev() {
        let mut rhs = aug[i][n];
        for j in i + 1..n {
            rhs -= aug[i][j] * x[j];
        }
        let diag = aug[i][i];
        x[i] = if diag.abs() < 1e-30 { 0.0 } else { rhs / diag };
    }
    x
}

// ─────────────────────────────────────────────────────────────────────────────
// § 7  Transient solver (explicit Euler)
// ─────────────────────────────────────────────────────────────────────────────

/// Advance the temperature field one explicit Euler time step.
///
/// `T_{n+1} = T_n + dt * C^{-1} (f - K T_n)`
///
/// The mass (capacitance) matrix C is assumed diagonal (lumped).
///
/// # Arguments
/// * `temps` – nodal temperatures at time n
/// * `k`     – global stiffness (conductivity) matrix
/// * `c`     – global capacitance matrix (full, but only diagonal is used here)
/// * `f`     – global load vector
/// * `dt`    – time step \[s\]
pub fn thermal_transient_euler(
    temps: &[f64],
    k: &[Vec<f64>],
    c: &[Vec<f64>],
    f: &[f64],
    dt: f64,
) -> Vec<f64> {
    let n = temps.len();
    let mut t_new = vec![0.0_f64; n];
    for i in 0..n {
        // residual r_i = f_i - (K T)_i
        let kt_i: f64 = k[i]
            .iter()
            .zip(temps.iter())
            .map(|(kij, tj)| kij * tj)
            .sum();
        let r_i = f[i] - kt_i;
        let c_ii = c[i][i].max(1e-30);
        t_new[i] = temps[i] + dt * r_i / c_ii;
    }
    t_new
}

// ─────────────────────────────────────────────────────────────────────────────
// § 8  Element heat flux
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the heat-flux vector for a single tet element.
///
/// `q = -D B T_e`  where D = diag(kx,ky,kz), B = shape-function gradient matrix.
///
/// # Arguments
/// * `temps` – global nodal temperatures (indexed by element node indices)
/// * `nodes` – four corner positions \[\[x,y,z\\]; 4]
/// * `k`     – conductivity components \[kx, ky, kz\]
pub fn heat_flux_element(temps: &[f64], nodes: &[[f64; 3]; 4], k: [f64; 3]) -> [f64; 3] {
    // We need the shape-function gradients.  Recompute via the same approach
    // used in tet_thermal_conductivity_matrix.
    let p0 = nodes[0];
    let p1 = nodes[1];
    let p2 = nodes[2];
    let p3 = nodes[3];

    let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let c_vec = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];

    let bxc = [
        b[1] * c_vec[2] - b[2] * c_vec[1],
        b[2] * c_vec[0] - b[0] * c_vec[2],
        b[0] * c_vec[1] - b[1] * c_vec[0],
    ];
    let det_j = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];
    if det_j.abs() < 1e-30 {
        return [0.0; 3];
    }
    let inv_det = 1.0 / det_j;

    let j_inv = [
        [
            (b[1] * c_vec[2] - b[2] * c_vec[1]) * inv_det,
            -(a[1] * c_vec[2] - a[2] * c_vec[1]) * inv_det,
            (a[1] * b[2] - a[2] * b[1]) * inv_det,
        ],
        [
            -(b[0] * c_vec[2] - b[2] * c_vec[0]) * inv_det,
            (a[0] * c_vec[2] - a[2] * c_vec[0]) * inv_det,
            -(a[0] * b[2] - a[2] * b[0]) * inv_det,
        ],
        [
            (b[0] * c_vec[1] - b[1] * c_vec[0]) * inv_det,
            -(a[0] * c_vec[1] - a[1] * c_vec[0]) * inv_det,
            (a[0] * b[1] - a[1] * b[0]) * inv_det,
        ],
    ];

    let dn: [[f64; 3]; 4] = [
        [
            -(j_inv[0][0] + j_inv[1][0] + j_inv[2][0]),
            -(j_inv[0][1] + j_inv[1][1] + j_inv[2][1]),
            -(j_inv[0][2] + j_inv[1][2] + j_inv[2][2]),
        ],
        [j_inv[0][0], j_inv[0][1], j_inv[0][2]],
        [j_inv[1][0], j_inv[1][1], j_inv[1][2]],
        [j_inv[2][0], j_inv[2][1], j_inv[2][2]],
    ];

    // grad T = B T_e
    let mut grad_t = [0.0_f64; 3];
    for m in 0..4 {
        let t_m = temps[m];
        grad_t[0] += dn[m][0] * t_m;
        grad_t[1] += dn[m][1] * t_m;
        grad_t[2] += dn[m][2] * t_m;
    }

    [-k[0] * grad_t[0], -k[1] * grad_t[1], -k[2] * grad_t[2]]
}

// ─────────────────────────────────────────────────────────────────────────────
// § 9  Thermal stress coupling
// ─────────────────────────────────────────────────────────────────────────────

/// Compute thermal stress components due to a temperature change.
///
/// For an isotropic material the thermal strain is `ε_th = α ΔT`,
/// and the induced stress (with the body fully constrained) is
/// `σ_th = -E/(1-2ν) α ΔT` (volumetric component only).
///
/// Returns a vector of length `n` with one value per node.
///
/// # Arguments
/// * `delta_t` – nodal temperature increment ΔT \[K\]
/// * `alpha`   – thermal expansion coefficient \[1/K\]
/// * `e`       – Young's modulus \[Pa\]
/// * `nu`      – Poisson's ratio
pub fn thermal_stress_coupling(delta_t: &[f64], alpha: f64, e: f64, nu: f64) -> Vec<f64> {
    let factor = -e / (1.0 - 2.0 * nu) * alpha;
    delta_t.iter().map(|&dt| factor * dt).collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// § 10  Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── ThermalNode ──────────────────────────────────────────────────────────

    #[test]
    fn test_thermal_node_free() {
        let n = ThermalNode::free([0.0, 0.0, 0.0], 300.0);
        assert!(!n.is_dirichlet);
        assert!((n.temperature - 300.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_node_dirichlet() {
        let n = ThermalNode::dirichlet([1.0, 2.0, 3.0], 500.0);
        assert!(n.is_dirichlet);
        assert!((n.prescribed_temp - 500.0).abs() < 1e-12);
    }

    // ── ThermalElement ───────────────────────────────────────────────────────

    #[test]
    fn test_thermal_element_new() {
        let e = ThermalElement::new([0, 1, 2, 3], [1.0, 2.0, 3.0]);
        assert_eq!(e.nodes, [0, 1, 2, 3]);
        assert!((e.conductivity[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_element_isotropic() {
        let e = ThermalElement::isotropic([0, 1, 2, 3], 5.0);
        assert!((e.conductivity[0] - 5.0).abs() < 1e-12);
        assert!((e.conductivity[1] - 5.0).abs() < 1e-12);
        assert!((e.conductivity[2] - 5.0).abs() < 1e-12);
    }

    // ── fourier_heat_flux ────────────────────────────────────────────────────

    #[test]
    fn test_fourier_heat_flux_sign() {
        // Positive gradient → negative flux (heat flows downhill)
        assert!(fourier_heat_flux(1.0, 1.0) < 0.0);
    }

    #[test]
    fn test_fourier_heat_flux_magnitude() {
        let q = fourier_heat_flux(2.0, 3.0);
        assert!((q + 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_fourier_heat_flux_zero_gradient() {
        assert_eq!(fourier_heat_flux(5.0, 0.0), 0.0);
    }

    // ── tet_thermal_conductivity_matrix ─────────────────────────────────────

    fn unit_tet_nodes() -> [[f64; 3]; 4] {
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    }

    #[test]
    fn test_tet_ke_symmetric() {
        let nodes = unit_tet_nodes();
        let ke = tet_thermal_conductivity_matrix(&nodes, [1.0, 1.0, 1.0]);
        for (i, row) in ke.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - ke[j][i]).abs() < 1e-12, "asymmetry at ({i},{j})");
            }
        }
    }

    #[test]
    fn test_tet_ke_row_sum_zero() {
        // Rows of K sum to zero (constant temperature gives zero heat flux)
        let nodes = unit_tet_nodes();
        let ke = tet_thermal_conductivity_matrix(&nodes, [1.0, 1.0, 1.0]);
        for (i, row) in ke.iter().enumerate() {
            let row_sum: f64 = row.iter().sum();
            assert!(row_sum.abs() < 1e-10, "row {i} sum = {row_sum}");
        }
    }

    #[test]
    fn test_tet_ke_degenerate() {
        // All nodes coincident → should return zeros
        let nodes = [[0.0_f64; 3]; 4];
        let ke = tet_thermal_conductivity_matrix(&nodes, [1.0, 1.0, 1.0]);
        for row in &ke {
            for &v in row {
                assert_eq!(v, 0.0);
            }
        }
    }

    #[test]
    fn test_tet_ke_scales_with_conductivity() {
        let nodes = unit_tet_nodes();
        let ke1 = tet_thermal_conductivity_matrix(&nodes, [1.0, 1.0, 1.0]);
        let ke2 = tet_thermal_conductivity_matrix(&nodes, [2.0, 2.0, 2.0]);
        for (i, row2) in ke2.iter().enumerate() {
            for (j, &v2) in row2.iter().enumerate() {
                assert!((v2 - 2.0 * ke1[i][j]).abs() < 1e-12);
            }
        }
    }

    // ── assemble_thermal_stiffness ───────────────────────────────────────────

    #[test]
    fn test_assemble_single_element() {
        let nodes = vec![
            ThermalNode::free([0.0, 0.0, 0.0], 0.0),
            ThermalNode::free([1.0, 0.0, 0.0], 0.0),
            ThermalNode::free([0.0, 1.0, 0.0], 0.0),
            ThermalNode::free([0.0, 0.0, 1.0], 0.0),
        ];
        let elements = vec![ThermalElement::isotropic([0, 1, 2, 3], 1.0)];
        let k = assemble_thermal_stiffness(&elements, &nodes);
        assert_eq!(k.len(), 4);
        // Symmetry check
        for (i, row) in k.iter().enumerate() {
            for (j, &val) in row.iter().enumerate() {
                assert!((val - k[j][i]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_assemble_size() {
        let nodes: Vec<ThermalNode> = (0..5)
            .map(|i| ThermalNode::free([i as f64, 0.0, 0.0], 0.0))
            .collect();
        let elements: Vec<ThermalElement> = vec![];
        let k = assemble_thermal_stiffness(&elements, &nodes);
        assert_eq!(k.len(), 5);
        assert_eq!(k[0].len(), 5);
    }

    // ── apply_thermal_dirichlet ──────────────────────────────────────────────

    #[test]
    fn test_dirichlet_diagonal_one() {
        let nodes = vec![
            ThermalNode::dirichlet([0.0, 0.0, 0.0], 100.0),
            ThermalNode::free([1.0, 0.0, 0.0], 0.0),
        ];
        let mut k = vec![vec![2.0, -1.0], vec![-1.0, 2.0]];
        let mut f = vec![0.0, 50.0];
        apply_thermal_dirichlet(&mut k, &mut f, &nodes);
        assert!((k[0][0] - 1.0).abs() < 1e-12);
        assert!((k[0][1]).abs() < 1e-12);
        assert!((f[0] - 100.0).abs() < 1e-12);
    }

    // ── solve_thermal_steady ─────────────────────────────────────────────────

    #[test]
    fn test_solve_thermal_steady_identity() {
        let k = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let f = vec![3.0, 7.0];
        let t = solve_thermal_steady(&k, &f);
        assert!((t[0] - 3.0).abs() < 1e-10);
        assert!((t[1] - 7.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_thermal_steady_two_node() {
        // Simple 1-D rod: K = [[1,-1],[-1,1]] with Dirichlet applied
        let k = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let f = vec![100.0, 200.0];
        let t = solve_thermal_steady(&k, &f);
        assert!((t[0] - 100.0).abs() < 1e-10);
        assert!((t[1] - 200.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_thermal_steady_residual() {
        // Verify K T ≈ f
        let k = vec![vec![4.0, -1.0], vec![-1.0, 3.0]];
        let f = vec![9.0, 8.0];
        let t = solve_thermal_steady(&k, &f);
        let r0 = 4.0 * t[0] - t[1] - 9.0;
        let r1 = -t[0] + 3.0 * t[1] - 8.0;
        assert!(r0.abs() < 1e-10);
        assert!(r1.abs() < 1e-10);
    }

    // ── thermal_transient_euler ──────────────────────────────────────────────

    #[test]
    fn test_transient_euler_no_forcing() {
        // With f=0 and K=I, temperature should decay: T_new = T - dt/C * K * T
        let k = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let c = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let f = vec![0.0, 0.0];
        let t = vec![10.0, 20.0];
        let t_new = thermal_transient_euler(&t, &k, &c, &f, 0.1);
        assert!((t_new[0] - 9.0).abs() < 1e-10);
        assert!((t_new[1] - 18.0).abs() < 1e-10);
    }

    #[test]
    fn test_transient_euler_length() {
        let k = vec![vec![1.0]];
        let c = vec![vec![1.0]];
        let f = vec![0.0];
        let t = vec![5.0];
        let t_new = thermal_transient_euler(&t, &k, &c, &f, 0.01);
        assert_eq!(t_new.len(), 1);
    }

    #[test]
    fn test_transient_euler_steady_state() {
        // If T is the steady-state solution, K T = f, so residual = 0 → T unchanged
        let k = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let c = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let f = vec![5.0, 10.0];
        let t = vec![5.0, 10.0]; // steady state: K*T = f
        let t_new = thermal_transient_euler(&t, &k, &c, &f, 1.0);
        assert!((t_new[0] - t[0]).abs() < 1e-10);
        assert!((t_new[1] - t[1]).abs() < 1e-10);
    }

    // ── heat_flux_element ────────────────────────────────────────────────────

    #[test]
    fn test_heat_flux_uniform_temperature() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let temps = [300.0, 300.0, 300.0, 300.0];
        let q = heat_flux_element(&temps, &nodes, [1.0, 1.0, 1.0]);
        assert!(q[0].abs() < 1e-10);
        assert!(q[1].abs() < 1e-10);
        assert!(q[2].abs() < 1e-10);
    }

    #[test]
    fn test_heat_flux_direction() {
        // Temperature gradient purely in x direction
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        // T(x) = 100x → T values at nodes
        let temps = [0.0, 100.0, 0.0, 0.0];
        let q = heat_flux_element(&temps, &nodes, [1.0, 1.0, 1.0]);
        // grad T in x = 100 → q_x = -k_x * 100 = -100
        assert!((q[0] + 100.0).abs() < 1e-8);
    }

    // ── thermal_stress_coupling ──────────────────────────────────────────────

    #[test]
    fn test_thermal_stress_sign() {
        // Positive ΔT with constrained body → compressive (negative) stress
        let sigma = thermal_stress_coupling(&[10.0], 1e-5, 2e11, 0.3);
        assert!(sigma[0] < 0.0);
    }

    #[test]
    fn test_thermal_stress_zero_delta_t() {
        let sigma = thermal_stress_coupling(&[0.0, 0.0], 1e-5, 2e11, 0.3);
        assert_eq!(sigma[0], 0.0);
        assert_eq!(sigma[1], 0.0);
    }

    #[test]
    fn test_thermal_stress_scales_linearly() {
        let s1 = thermal_stress_coupling(&[1.0], 1e-5, 2e11, 0.3);
        let s2 = thermal_stress_coupling(&[2.0], 1e-5, 2e11, 0.3);
        assert!((s2[0] - 2.0 * s1[0]).abs() < 1e-4);
    }

    #[test]
    fn test_thermal_stress_length() {
        let dt = vec![1.0; 10];
        let sigma = thermal_stress_coupling(&dt, 1e-5, 2e11, 0.3);
        assert_eq!(sigma.len(), 10);
    }

    #[test]
    fn test_tet_ke_positive_diagonal() {
        let nodes = unit_tet_nodes();
        let ke = tet_thermal_conductivity_matrix(&nodes, [1.0, 1.0, 1.0]);
        for (i, row) in ke.iter().enumerate() {
            assert!(row[i] >= 0.0, "diagonal [{i}] is negative");
        }
    }

    #[test]
    fn test_solve_3x3() {
        // 3×3 system with known solution
        let k = vec![
            vec![2.0, 1.0, 0.0],
            vec![1.0, 3.0, 1.0],
            vec![0.0, 1.0, 2.0],
        ];
        let f = vec![5.0, 10.0, 7.0];
        let t = solve_thermal_steady(&k, &f);
        let r0 = 2.0 * t[0] + t[1] - 5.0;
        let r1 = t[0] + 3.0 * t[1] + t[2] - 10.0;
        let r2 = t[1] + 2.0 * t[2] - 7.0;
        assert!(r0.abs() < 1e-9);
        assert!(r1.abs() < 1e-9);
        assert!(r2.abs() < 1e-9);
    }

    #[test]
    fn test_fourier_heat_flux_negative_gradient() {
        // Negative gradient → positive flux
        let q = fourier_heat_flux(1.0, -5.0);
        assert!((q - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_thermal_transient_euler_heating() {
        // With f > K*T the temperature should increase
        let k = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let c = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let t = vec![0.0, 0.0];
        let f = vec![10.0, 5.0];
        let t_new = thermal_transient_euler(&t, &k, &c, &f, 0.1);
        assert!(t_new[0] > t[0]);
        assert!(t_new[1] > t[1]);
    }

    #[test]
    fn test_heat_flux_z_gradient() {
        // Temperature gradient purely in z direction
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let temps = [0.0, 0.0, 0.0, 50.0];
        let q = heat_flux_element(&temps, &nodes, [1.0, 1.0, 1.0]);
        // grad T in z = 50 → q_z = -50
        assert!((q[2] + 50.0).abs() < 1e-8);
    }

    #[test]
    fn test_assemble_two_elements_accumulate() {
        // Two elements sharing nodes — stiffness should accumulate at shared nodes
        let nodes = vec![
            ThermalNode::free([0.0, 0.0, 0.0], 0.0),
            ThermalNode::free([1.0, 0.0, 0.0], 0.0),
            ThermalNode::free([0.0, 1.0, 0.0], 0.0),
            ThermalNode::free([0.0, 0.0, 1.0], 0.0),
            ThermalNode::free([1.0, 1.0, 0.0], 0.0),
        ];
        let elements = vec![
            ThermalElement::isotropic([0, 1, 2, 3], 1.0),
            ThermalElement::isotropic([1, 2, 4, 3], 1.0),
        ];
        let k = assemble_thermal_stiffness(&elements, &nodes);
        // Node 1 and 2 are shared; their diagonal entries should be > single-element contribution
        let k_single_nodes = {
            let elems1 = vec![ThermalElement::isotropic([0, 1, 2, 3], 1.0)];
            assemble_thermal_stiffness(&elems1, &nodes)
        };
        assert!(k[1][1] > k_single_nodes[1][1]);
    }

    #[test]
    fn test_thermal_stress_coupling_formula() {
        // σ = -E/(1-2ν) * α * ΔT
        let e = 1.0;
        let nu = 0.25;
        let alpha = 1.0;
        let dt = 2.0;
        let sigma = thermal_stress_coupling(&[dt], alpha, e, nu);
        let expected = -e / (1.0 - 2.0 * nu) * alpha * dt;
        assert!((sigma[0] - expected).abs() < 1e-12);
    }

    #[test]
    fn test_solve_thermal_steady_single() {
        // 1×1 system
        let k = vec![vec![5.0]];
        let f = vec![15.0];
        let t = solve_thermal_steady(&k, &f);
        assert!((t[0] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_tet_ke_unit_volume() {
        // Unit tet has volume 1/6; check the sum of all entries is zero
        let nodes = unit_tet_nodes();
        let ke = tet_thermal_conductivity_matrix(&nodes, [1.0, 1.0, 1.0]);
        let total: f64 = ke.iter().flat_map(|r| r.iter()).sum();
        assert!(total.abs() < 1e-10, "total KE sum = {total}");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 11  ThermalFemParams
// ─────────────────────────────────────────────────────────────────────────────

/// Parameters for a thermal FEM analysis.
#[derive(Debug, Clone)]
pub struct ThermalFemParams {
    /// Thermal conductivity tensor diagonal \[kx, ky, kz\] \[W/(m·K)\].
    pub conductivity: [f64; 3],
    /// Specific heat capacity c_p \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Mass density ρ \[kg/m³\].
    pub density: f64,
    /// Convection coefficient h \[W/(m²·K)\].
    pub convection_coeff: f64,
    /// Ambient / reference temperature for convection \[K\].
    pub t_ambient: f64,
}

impl ThermalFemParams {
    /// Create isotropic thermal parameters.
    pub fn isotropic(k: f64, specific_heat: f64, density: f64) -> Self {
        Self {
            conductivity: [k, k, k],
            specific_heat,
            density,
            convection_coeff: 0.0,
            t_ambient: 293.15_f64,
        }
    }

    /// Volumetric heat capacity ρ c_p \[J/(m³·K)\].
    pub fn volumetric_heat_capacity(&self) -> f64 {
        self.density * self.specific_heat
    }

    /// Thermal diffusivity α = k / (ρ c_p) \[m²/s\].
    pub fn thermal_diffusivity(&self) -> f64 {
        let k_avg = (self.conductivity[0] + self.conductivity[1] + self.conductivity[2]) / 3.0_f64;
        let rho_cp = self.volumetric_heat_capacity();
        if rho_cp < 1e-30 { 0.0 } else { k_avg / rho_cp }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 12  ThermalElementMatrix — consistent heat capacity + conductivity
// ─────────────────────────────────────────────────────────────────────────────

/// Consistent heat-capacity matrix for a linear tetrahedron.
///
/// `M_e = ρ c_p V / 20 * [2 1 1 1; 1 2 1 1; 1 1 2 1; 1 1 1 2]`
pub fn tet_consistent_heat_capacity(
    params: &ThermalFemParams,
    nodes: &[[f64; 3]; 4],
) -> [[f64; 4]; 4] {
    let p0 = nodes[0];
    let p1 = nodes[1];
    let p2 = nodes[2];
    let p3 = nodes[3];
    let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
    let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
    let c = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
    let bxc = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    let det_j = a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2];
    let vol = det_j.abs() / 6.0_f64;
    let factor = params.volumetric_heat_capacity() * vol / 20.0_f64;
    let mut m = [[0.0_f64; 4]; 4];
    for (i, row) in m.iter_mut().enumerate() {
        for (j, val) in row.iter_mut().enumerate() {
            *val = factor * if i == j { 2.0_f64 } else { 1.0_f64 };
        }
    }
    m
}

/// Combined element matrices: conductivity (K_e) and heat capacity (M_e).
pub struct ThermalElementMatrix {
    /// 4×4 element conductivity matrix.
    pub conductivity: [[f64; 4]; 4],
    /// 4×4 consistent heat capacity matrix.
    pub capacity: [[f64; 4]; 4],
}

impl ThermalElementMatrix {
    /// Compute both matrices for a linear tet element.
    pub fn compute(params: &ThermalFemParams, nodes: &[[f64; 3]; 4]) -> Self {
        let conductivity = tet_thermal_conductivity_matrix(nodes, params.conductivity);
        let capacity = tet_consistent_heat_capacity(params, nodes);
        Self {
            conductivity,
            capacity,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 13  ThermalBoundaryConditions
// ─────────────────────────────────────────────────────────────────────────────

/// Type of thermal boundary condition.
#[derive(Debug, Clone, PartialEq)]
pub enum ThermalBcKind {
    /// Dirichlet — prescribed temperature.
    Dirichlet,
    /// Neumann — prescribed heat flux \[W/m²\].
    Neumann,
    /// Robin / convection — h (T - T_inf).
    Robin,
}

/// A thermal boundary condition applied to a node.
#[derive(Debug, Clone)]
pub struct ThermalBc {
    /// Global node index.
    pub node_index: usize,
    /// Kind of BC.
    pub kind: ThermalBcKind,
    /// Value: temperature \[K\] for Dirichlet, flux \[W/m²\] for Neumann, h*T_inf for Robin.
    pub value: f64,
    /// Convection coefficient \[W/(m²·K)\] (used only for Robin BC).
    pub h_conv: f64,
}

impl ThermalBc {
    /// Create a Dirichlet BC: fixed temperature.
    pub fn dirichlet(node_index: usize, temperature: f64) -> Self {
        Self {
            node_index,
            kind: ThermalBcKind::Dirichlet,
            value: temperature,
            h_conv: 0.0_f64,
        }
    }

    /// Create a Neumann BC: prescribed inward heat flux \[W/m²\].
    pub fn neumann(node_index: usize, flux: f64) -> Self {
        Self {
            node_index,
            kind: ThermalBcKind::Neumann,
            value: flux,
            h_conv: 0.0_f64,
        }
    }

    /// Create a Robin (convection) BC: h*(T_inf - T).
    pub fn robin(node_index: usize, h_conv: f64, t_inf: f64) -> Self {
        Self {
            node_index,
            kind: ThermalBcKind::Robin,
            value: t_inf,
            h_conv,
        }
    }
}

/// Apply a list of `ThermalBc`s to the global system (K, f).
///
/// - Dirichlet: penalty / elimination method.
/// - Neumann: add flux to load vector.
/// - Robin: add h_conv to diagonal of K, add h_conv * T_inf to f.
pub fn apply_thermal_bcs(k: &mut [Vec<f64>], f: &mut [f64], bcs: &[ThermalBc]) {
    let n = f.len();
    for bc in bcs {
        let i = bc.node_index;
        if i >= n {
            continue;
        }
        match bc.kind {
            ThermalBcKind::Dirichlet => {
                let t_p = bc.value;
                for j in 0..n {
                    if j != i {
                        f[j] -= k[j][i] * t_p;
                        k[j][i] = 0.0_f64;
                        k[i][j] = 0.0_f64;
                    }
                }
                k[i][i] = 1.0_f64;
                f[i] = t_p;
            }
            ThermalBcKind::Neumann => {
                f[i] += bc.value;
            }
            ThermalBcKind::Robin => {
                k[i][i] += bc.h_conv;
                f[i] += bc.h_conv * bc.value;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 14  TransientThermalSolver — theta-method (Euler / Crank–Nicolson)
// ─────────────────────────────────────────────────────────────────────────────

/// Transient solver using the generalised θ-method.
///
/// θ = 0  → explicit Euler.
/// θ = 0.5 → Crank–Nicolson (second-order accurate, unconditionally stable).
/// θ = 1  → implicit Euler.
pub struct TransientThermalSolver {
    /// θ parameter in \[0, 1\].
    pub theta: f64,
    /// Time step \[s\].
    pub dt: f64,
}

impl TransientThermalSolver {
    /// Create a solver with explicit Euler (θ = 0).
    pub fn explicit_euler(dt: f64) -> Self {
        Self { theta: 0.0_f64, dt }
    }

    /// Create a Crank–Nicolson solver (θ = 0.5).
    pub fn crank_nicolson(dt: f64) -> Self {
        Self { theta: 0.5_f64, dt }
    }

    /// Create an implicit Euler solver (θ = 1).
    pub fn implicit_euler(dt: f64) -> Self {
        Self { theta: 1.0_f64, dt }
    }

    /// Advance one time step.
    ///
    /// Solves `(M + θ dt K) T_{n+1} = (M - (1-θ) dt K) T_n + dt f`.
    ///
    /// Uses the lumped mass matrix (diagonal of `m_full`) for simplicity.
    pub fn step(&self, temps: &[f64], k: &[Vec<f64>], m: &[Vec<f64>], f: &[f64]) -> Vec<f64> {
        let n = temps.len();
        let dt = self.dt;
        let theta = self.theta;

        // Lumped capacitance: sum of each row of M
        let c_lumped: Vec<f64> = (0..n)
            .map(|i| m[i].iter().sum::<f64>().max(1e-30))
            .collect();

        // RHS: c_lumped * T_n - (1-theta)*dt*(K T_n) + dt*f
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            let kt_i: f64 = k[i]
                .iter()
                .zip(temps.iter())
                .map(|(kij, tj)| kij * tj)
                .sum();
            rhs[i] = c_lumped[i] * temps[i] - (1.0_f64 - theta) * dt * kt_i + dt * f[i];
        }

        if theta == 0.0_f64 {
            // Explicit: T_{n+1} = rhs / c_lumped
            rhs.iter()
                .zip(c_lumped.iter())
                .map(|(r, c)| r / c)
                .collect()
        } else {
            // Implicit: (c_lumped + theta*dt*K_diag)*T_{n+1} ≈ rhs  (diagonal approx)
            // For a full solve we use the diagonal approximation to avoid linear algebra
            (0..n)
                .map(|i| {
                    let lhs_diag = c_lumped[i] + theta * dt * k[i][i];
                    if lhs_diag.abs() < 1e-30 {
                        0.0_f64
                    } else {
                        rhs[i] / lhs_diag
                    }
                })
                .collect()
        }
    }

    /// Run `n_steps` time integration steps, returning temperature history.
    pub fn integrate(
        &self,
        t0: &[f64],
        k: &[Vec<f64>],
        m: &[Vec<f64>],
        f: &[f64],
        n_steps: usize,
    ) -> Vec<Vec<f64>> {
        let mut history = Vec::with_capacity(n_steps + 1);
        let mut temps = t0.to_vec();
        history.push(temps.clone());
        for _ in 0..n_steps {
            temps = self.step(&temps, k, m, f);
            history.push(temps.clone());
        }
        history
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 15  RadiationBoundary — Stefan–Boltzmann
// ─────────────────────────────────────────────────────────────────────────────

/// Stefan–Boltzmann radiation boundary condition.
///
/// Heat flux: `q_rad = ε σ (T_surf⁴ - T_surr⁴)` \[W/m²\].
pub struct RadiationBoundary {
    /// Surface emissivity (0 ≤ ε ≤ 1).
    pub emissivity: f64,
    /// Surrounding temperature T_surr \[K\].
    pub t_surrounding: f64,
}

/// Stefan–Boltzmann constant σ \[W/(m²·K⁴)\].
pub const STEFAN_BOLTZMANN: f64 = 5.670374419e-8_f64;

impl RadiationBoundary {
    /// Create a radiation boundary.
    pub fn new(emissivity: f64, t_surrounding: f64) -> Self {
        Self {
            emissivity,
            t_surrounding,
        }
    }

    /// Net radiation heat flux from surface at temperature `t_surf` \[W/m²\].
    ///
    /// Positive value = heat leaving the surface.
    pub fn heat_flux(&self, t_surf: f64) -> f64 {
        self.emissivity * STEFAN_BOLTZMANN * (t_surf.powi(4) - self.t_surrounding.powi(4))
    }

    /// Linearised radiation conductance at temperature `t_surf` \[W/(m²·K)\].
    ///
    /// Used to assemble a linearised radiation term into the FEM matrix:
    /// `h_rad = 4 ε σ T_surf³`.
    pub fn linearised_conductance(&self, t_surf: f64) -> f64 {
        4.0_f64 * self.emissivity * STEFAN_BOLTZMANN * t_surf.powi(3)
    }

    /// Apply linearised radiation BC to node `i` of the global system.
    pub fn apply_to_node(&self, k: &mut [Vec<f64>], f: &mut [f64], node: usize, t_surf: f64) {
        if node >= f.len() {
            return;
        }
        let h_rad = self.linearised_conductance(t_surf);
        k[node][node] += h_rad;
        f[node] += h_rad * self.t_surrounding + self.heat_flux(t_surf) - h_rad * t_surf;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 16  ThermoElasticCoupling — thermal strains in structural elements
// ─────────────────────────────────────────────────────────────────────────────

/// Thermo-elastic coupling: compute the thermal strain vector and
/// equivalent nodal forces for a structural element.
pub struct ThermoElasticCoupling {
    /// Thermal expansion coefficient α \[1/K\].
    pub alpha: f64,
    /// Young's modulus E \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
}

impl ThermoElasticCoupling {
    /// Create a thermo-elastic coupling for an isotropic material.
    pub fn new(alpha: f64, young_modulus: f64, poisson_ratio: f64) -> Self {
        Self {
            alpha,
            young_modulus,
            poisson_ratio,
        }
    }

    /// Volumetric thermal strain for a temperature increment ΔT.
    pub fn volumetric_strain(&self, delta_t: f64) -> f64 {
        3.0_f64 * self.alpha * delta_t
    }

    /// Engineering thermal strain vector `[εx, εy, εz, γxy, γyz, γxz]` for ΔT.
    ///
    /// For an isotropic material, `εx = εy = εz = α ΔT`, shear strains = 0.
    pub fn strain_vector(&self, delta_t: f64) -> [f64; 6] {
        let eps = self.alpha * delta_t;
        [eps, eps, eps, 0.0_f64, 0.0_f64, 0.0_f64]
    }

    /// Thermal stress vector `σ_th = -D ε_th` for a fully-constrained body \[Pa\].
    ///
    /// Returns `[σx, σy, σz, τxy, τyz, τxz]`.
    pub fn stress_vector(&self, delta_t: f64) -> [f64; 6] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let factor = -e * self.alpha * delta_t / (1.0_f64 - 2.0_f64 * nu);
        [factor, factor, factor, 0.0_f64, 0.0_f64, 0.0_f64]
    }

    /// Equivalent nodal thermal load for a tet element.
    ///
    /// `f_th = ∫ B^T D ε_th dV = V B^T D ε_th` (constant strain tet).
    ///
    /// Returns a 12-component vector (3 DOF × 4 nodes).
    pub fn nodal_thermal_load(&self, nodes: &[[f64; 3]; 4], delta_t_avg: f64) -> Vec<f64> {
        let stress = self.stress_vector(delta_t_avg);
        // Volume of tet
        let p0 = nodes[0];
        let p1 = nodes[1];
        let p2 = nodes[2];
        let p3 = nodes[3];
        let a = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
        let b = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
        let cv = [p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]];
        let bxc = [
            b[1] * cv[2] - b[2] * cv[1],
            b[2] * cv[0] - b[0] * cv[2],
            b[0] * cv[1] - b[1] * cv[0],
        ];
        let vol = (a[0] * bxc[0] + a[1] * bxc[1] + a[2] * bxc[2]).abs() / 6.0_f64;
        // For a tet, B^T σ_th gives nodal forces; simplified: distribute equally
        let nodal_force = vol / 4.0_f64;
        let mut load = vec![0.0_f64; 12];
        for nd in 0..4 {
            load[nd * 3] = nodal_force * stress[0];
            load[nd * 3 + 1] = nodal_force * stress[1];
            load[nd * 3 + 2] = nodal_force * stress[2];
        }
        load
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 17  PhaseThermal — latent heat model
// ─────────────────────────────────────────────────────────────────────────────

/// Phase-change thermal model with latent heat.
///
/// Uses the apparent heat capacity method: the effective c_p is increased
/// in the mushy zone to account for latent heat release/absorption.
#[derive(Debug, Clone)]
pub struct PhaseThermal {
    /// Solidus temperature T_s \[K\].
    pub t_solidus: f64,
    /// Liquidus temperature T_l \[K\].
    pub t_liquidus: f64,
    /// Latent heat L \[J/kg\].
    pub latent_heat: f64,
    /// Baseline specific heat c_p \[J/(kg·K)\].
    pub specific_heat: f64,
    /// Mass density ρ \[kg/m³\].
    pub density: f64,
}

impl PhaseThermal {
    /// Create a phase-change model.
    pub fn new(
        t_solidus: f64,
        t_liquidus: f64,
        latent_heat: f64,
        specific_heat: f64,
        density: f64,
    ) -> Self {
        Self {
            t_solidus,
            t_liquidus,
            latent_heat,
            specific_heat,
            density,
        }
    }

    /// Liquid fraction at temperature T (0 = solid, 1 = liquid).
    pub fn liquid_fraction(&self, temperature: f64) -> f64 {
        if temperature <= self.t_solidus {
            0.0_f64
        } else if temperature >= self.t_liquidus {
            1.0_f64
        } else {
            (temperature - self.t_solidus) / (self.t_liquidus - self.t_solidus)
        }
    }

    /// Apparent heat capacity (including latent heat contribution) \[J/(kg·K)\].
    pub fn apparent_heat_capacity(&self, temperature: f64) -> f64 {
        if temperature > self.t_solidus && temperature < self.t_liquidus {
            let delta_t = self.t_liquidus - self.t_solidus;
            self.specific_heat + self.latent_heat / delta_t
        } else {
            self.specific_heat
        }
    }

    /// Volumetric apparent heat capacity ρ c_eff \[J/(m³·K)\].
    pub fn volumetric_apparent_capacity(&self, temperature: f64) -> f64 {
        self.density * self.apparent_heat_capacity(temperature)
    }

    /// Enthalpy at temperature T \[J/kg\] (relative to T_solidus).
    pub fn enthalpy(&self, temperature: f64) -> f64 {
        if temperature <= self.t_solidus {
            self.specific_heat * (temperature - self.t_solidus)
        } else if temperature >= self.t_liquidus {
            self.latent_heat + self.specific_heat * (temperature - self.t_solidus)
        } else {
            let fl = self.liquid_fraction(temperature);
            fl * self.latent_heat + self.specific_heat * (temperature - self.t_solidus)
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 18  1-D rod steady-state helper
// ─────────────────────────────────────────────────────────────────────────────

/// Solve a 1-D rod steady-state thermal problem with two Dirichlet BCs.
///
/// The rod has `n_elem` equal elements of length `L / n_elem`, conductivity `k`.
/// Returns the nodal temperature vector.
pub fn solve_1d_rod_steady(
    n_elem: usize,
    length: f64,
    k: f64,
    t_left: f64,
    t_right: f64,
) -> Vec<f64> {
    let n_nodes = n_elem + 1;
    let le = length / n_elem as f64;
    let ke_val = k / le;

    // Build global stiffness
    let mut k_global = vec![vec![0.0_f64; n_nodes]; n_nodes];
    for e in 0..n_elem {
        let i0 = e;
        let i1 = e + 1;
        k_global[i0][i0] += ke_val;
        k_global[i0][i1] -= ke_val;
        k_global[i1][i0] -= ke_val;
        k_global[i1][i1] += ke_val;
    }

    let mut f = vec![0.0_f64; n_nodes];
    // Apply Dirichlet BCs
    for j in 0..n_nodes {
        if j != 0 {
            f[j] -= k_global[j][0] * t_left;
            k_global[j][0] = 0.0_f64;
            k_global[0][j] = 0.0_f64;
        }
    }
    k_global[0][0] = 1.0_f64;
    f[0] = t_left;

    let last = n_nodes - 1;
    for j in 0..n_nodes {
        if j != last {
            f[j] -= k_global[j][last] * t_right;
            k_global[j][last] = 0.0_f64;
            k_global[last][j] = 0.0_f64;
        }
    }
    k_global[last][last] = 1.0_f64;
    f[last] = t_right;

    solve_thermal_steady(&k_global, &f)
}

// ─────────────────────────────────────────────────────────────────────────────
// § 19  Lumped-mass heat capacity assembly
// ─────────────────────────────────────────────────────────────────────────────

/// Assemble the global lumped heat capacity matrix for a mesh.
///
/// Lumping strategy: sum the row of each element's consistent capacity matrix
/// and put it on the diagonal.
pub fn assemble_lumped_capacity(
    elements: &[ThermalElement],
    nodes: &[ThermalNode],
    params: &ThermalFemParams,
) -> Vec<Vec<f64>> {
    let n = nodes.len();
    let mut c_global = vec![vec![0.0_f64; n]; n];

    for elem in elements {
        let pos: [[f64; 3]; 4] = [
            nodes[elem.nodes[0]].position,
            nodes[elem.nodes[1]].position,
            nodes[elem.nodes[2]].position,
            nodes[elem.nodes[3]].position,
        ];
        let me = tet_consistent_heat_capacity(params, &pos);
        for a in 0..4 {
            // Lump: add row sum to diagonal
            let row_sum: f64 = me[a].iter().sum();
            c_global[elem.nodes[a]][elem.nodes[a]] += row_sum;
        }
    }
    c_global
}

// ─────────────────────────────────────────────────────────────────────────────
// § 20  Additional tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests_extended {
    use super::*;

    // ── ThermalFemParams ──────────────────────────────────────────────────────

    #[test]
    fn test_params_isotropic_conductivity() {
        let p = ThermalFemParams::isotropic(50.0_f64, 500.0_f64, 7800.0_f64);
        assert!((p.conductivity[0] - 50.0_f64).abs() < 1e-12);
        assert!((p.conductivity[1] - 50.0_f64).abs() < 1e-12);
        assert!((p.conductivity[2] - 50.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_params_volumetric_heat_capacity() {
        let p = ThermalFemParams::isotropic(1.0_f64, 500.0_f64, 8000.0_f64);
        let expected = 500.0_f64 * 8000.0_f64;
        assert!((p.volumetric_heat_capacity() - expected).abs() < 1e-6);
    }

    #[test]
    fn test_params_thermal_diffusivity_positive() {
        let p = ThermalFemParams::isotropic(50.0_f64, 500.0_f64, 7800.0_f64);
        assert!(p.thermal_diffusivity() > 0.0_f64);
    }

    #[test]
    fn test_params_thermal_diffusivity_formula() {
        let p = ThermalFemParams::isotropic(10.0_f64, 100.0_f64, 1000.0_f64);
        // k_avg = 10, rho*cp = 100_000  → alpha = 10 / 100_000 = 1e-4
        let expected = 10.0_f64 / (100.0_f64 * 1000.0_f64);
        assert!((p.thermal_diffusivity() - expected).abs() < 1e-15);
    }

    // ── ThermalElementMatrix ──────────────────────────────────────────────────

    #[test]
    fn test_element_matrix_capacity_positive_diag() {
        let params = ThermalFemParams::isotropic(1.0_f64, 500.0_f64, 7800.0_f64);
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let em = ThermalElementMatrix::compute(&params, &nodes);
        for i in 0..4 {
            assert!(
                em.capacity[i][i] > 0.0_f64,
                "capacity diagonal must be positive"
            );
        }
    }

    #[test]
    fn test_element_matrix_capacity_symmetric() {
        let params = ThermalFemParams::isotropic(1.0_f64, 500.0_f64, 7800.0_f64);
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let em = ThermalElementMatrix::compute(&params, &nodes);
        for i in 0..4 {
            for j in 0..4 {
                assert!((em.capacity[i][j] - em.capacity[j][i]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_element_matrix_conductivity_row_sum_zero() {
        let params = ThermalFemParams::isotropic(1.0_f64, 500.0_f64, 7800.0_f64);
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let em = ThermalElementMatrix::compute(&params, &nodes);
        for i in 0..4 {
            let rs: f64 = em.conductivity[i].iter().sum();
            assert!(rs.abs() < 1e-10, "conductivity row {i} sum = {rs}");
        }
    }

    // ── ThermalBc ─────────────────────────────────────────────────────────────

    #[test]
    fn test_bc_dirichlet_kind() {
        let bc = ThermalBc::dirichlet(0, 300.0_f64);
        assert_eq!(bc.kind, ThermalBcKind::Dirichlet);
        assert!((bc.value - 300.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_bc_neumann_kind() {
        let bc = ThermalBc::neumann(1, 100.0_f64);
        assert_eq!(bc.kind, ThermalBcKind::Neumann);
        assert!((bc.value - 100.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_bc_robin_kind() {
        let bc = ThermalBc::robin(2, 25.0_f64, 293.15_f64);
        assert_eq!(bc.kind, ThermalBcKind::Robin);
        assert!((bc.h_conv - 25.0_f64).abs() < 1e-12);
    }

    #[test]
    fn test_apply_bcs_neumann_adds_to_rhs() {
        let mut k = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let mut f = vec![0.0_f64, 0.0_f64];
        let bcs = vec![ThermalBc::neumann(0, 50.0_f64)];
        apply_thermal_bcs(&mut k, &mut f, &bcs);
        assert!(
            (f[0] - 50.0_f64).abs() < 1e-12,
            "Neumann should add to f[0]"
        );
    }

    #[test]
    fn test_apply_bcs_robin_modifies_k_and_f() {
        let mut k = vec![vec![10.0_f64, 0.0_f64], vec![0.0_f64, 10.0_f64]];
        let mut f = vec![0.0_f64, 0.0_f64];
        let h = 5.0_f64;
        let t_inf = 300.0_f64;
        let bcs = vec![ThermalBc::robin(1, h, t_inf)];
        apply_thermal_bcs(&mut k, &mut f, &bcs);
        assert!(
            (k[1][1] - 15.0_f64).abs() < 1e-12,
            "Robin should add h to K diagonal"
        );
        assert!(
            (f[1] - h * t_inf).abs() < 1e-12,
            "Robin should add h*T_inf to f"
        );
    }

    // ── TransientThermalSolver ────────────────────────────────────────────────

    #[test]
    fn test_transient_solver_explicit_no_forcing_decay() {
        let solver = TransientThermalSolver::explicit_euler(0.01_f64);
        let k = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let m = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let f = vec![0.0_f64, 0.0_f64];
        let t0 = vec![100.0_f64, 50.0_f64];
        let t1 = solver.step(&t0, &k, &m, &f);
        // With decay, temperatures should decrease
        assert!(t1[0] < t0[0]);
        assert!(t1[1] < t0[1]);
    }

    #[test]
    fn test_transient_solver_integrate_returns_history() {
        let solver = TransientThermalSolver::crank_nicolson(0.1_f64);
        let k = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let m = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let f = vec![0.0_f64, 0.0_f64];
        let t0 = vec![100.0_f64, 50.0_f64];
        let history = solver.integrate(&t0, &k, &m, &f, 5);
        assert_eq!(history.len(), 6, "should have initial + 5 steps");
    }

    #[test]
    fn test_transient_solver_implicit_stable() {
        // Implicit Euler should remain finite even for large dt
        let solver = TransientThermalSolver::implicit_euler(10.0_f64);
        let k = vec![vec![100.0_f64, 0.0_f64], vec![0.0_f64, 100.0_f64]];
        let m = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let f = vec![0.0_f64, 0.0_f64];
        let t0 = vec![100.0_f64, 50.0_f64];
        let t1 = solver.step(&t0, &k, &m, &f);
        assert!(t1[0].is_finite());
        assert!(t1[1].is_finite());
    }

    // ── RadiationBoundary ─────────────────────────────────────────────────────

    #[test]
    fn test_radiation_heat_flux_positive_for_hot_surface() {
        let rb = RadiationBoundary::new(0.8_f64, 300.0_f64);
        let q = rb.heat_flux(1000.0_f64);
        assert!(q > 0.0_f64, "hot surface should radiate away heat");
    }

    #[test]
    fn test_radiation_heat_flux_negative_for_cold_surface() {
        let rb = RadiationBoundary::new(0.8_f64, 1000.0_f64);
        let q = rb.heat_flux(300.0_f64);
        assert!(q < 0.0_f64, "cold surface should absorb heat");
    }

    #[test]
    fn test_radiation_equilibrium_zero_flux() {
        let rb = RadiationBoundary::new(1.0_f64, 500.0_f64);
        let q = rb.heat_flux(500.0_f64);
        assert!(q.abs() < 1e-6, "at equilibrium flux should be zero: {q}");
    }

    #[test]
    fn test_radiation_linearised_conductance_positive() {
        let rb = RadiationBoundary::new(0.9_f64, 300.0_f64);
        let h = rb.linearised_conductance(500.0_f64);
        assert!(h > 0.0_f64, "linearised conductance must be positive");
    }

    #[test]
    fn test_radiation_stefan_boltzmann_const() {
        assert!((STEFAN_BOLTZMANN - 5.670374419e-8_f64).abs() < 1e-20);
    }

    // ── ThermoElasticCoupling ─────────────────────────────────────────────────

    #[test]
    fn test_thermo_elastic_strain_isotropy() {
        let te = ThermoElasticCoupling::new(1e-5_f64, 2e11_f64, 0.3_f64);
        let eps = te.strain_vector(100.0_f64);
        // All three normal strains should be equal
        assert!((eps[0] - eps[1]).abs() < 1e-20);
        assert!((eps[1] - eps[2]).abs() < 1e-20);
        // Shear strains should be zero
        assert_eq!(eps[3], 0.0_f64);
        assert_eq!(eps[4], 0.0_f64);
        assert_eq!(eps[5], 0.0_f64);
    }

    #[test]
    fn test_thermo_elastic_volumetric_strain() {
        let te = ThermoElasticCoupling::new(1e-5_f64, 2e11_f64, 0.3_f64);
        let vol = te.volumetric_strain(10.0_f64);
        // 3 * alpha * delta_t
        let expected = 3.0_f64 * 1e-5_f64 * 10.0_f64;
        assert!((vol - expected).abs() < 1e-20);
    }

    #[test]
    fn test_thermo_elastic_stress_sign() {
        // Positive heating in constrained body → compressive (negative) stress
        let te = ThermoElasticCoupling::new(1e-5_f64, 2e11_f64, 0.3_f64);
        let sigma = te.stress_vector(100.0_f64);
        assert!(sigma[0] < 0.0_f64, "thermal stress should be compressive");
    }

    #[test]
    fn test_thermo_elastic_nodal_load_length() {
        let te = ThermoElasticCoupling::new(1e-5_f64, 2e11_f64, 0.3_f64);
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let load = te.nodal_thermal_load(&nodes, 50.0_f64);
        assert_eq!(load.len(), 12, "12 dof for 4-node tet");
    }

    // ── PhaseThermal ──────────────────────────────────────────────────────────

    #[test]
    fn test_phase_liquid_fraction_below_solidus() {
        let pt = PhaseThermal::new(500.0_f64, 600.0_f64, 250_000.0_f64, 500.0_f64, 7000.0_f64);
        assert_eq!(pt.liquid_fraction(400.0_f64), 0.0_f64);
    }

    #[test]
    fn test_phase_liquid_fraction_above_liquidus() {
        let pt = PhaseThermal::new(500.0_f64, 600.0_f64, 250_000.0_f64, 500.0_f64, 7000.0_f64);
        assert_eq!(pt.liquid_fraction(700.0_f64), 1.0_f64);
    }

    #[test]
    fn test_phase_liquid_fraction_in_mushy_zone() {
        let pt = PhaseThermal::new(500.0_f64, 600.0_f64, 250_000.0_f64, 500.0_f64, 7000.0_f64);
        let fl = pt.liquid_fraction(550.0_f64);
        assert!(
            (fl - 0.5_f64).abs() < 1e-12,
            "fl at midpoint should be 0.5: {fl}"
        );
    }

    #[test]
    fn test_phase_apparent_heat_capacity_in_mushy() {
        let pt = PhaseThermal::new(500.0_f64, 600.0_f64, 250_000.0_f64, 500.0_f64, 7000.0_f64);
        let cp_eff = pt.apparent_heat_capacity(550.0_f64);
        // c_p + L/ΔT = 500 + 250000/100 = 3000
        let expected = 500.0_f64 + 250_000.0_f64 / 100.0_f64;
        assert!((cp_eff - expected).abs() < 1e-8);
    }

    #[test]
    fn test_phase_enthalpy_increases_with_temperature() {
        let pt = PhaseThermal::new(500.0_f64, 600.0_f64, 250_000.0_f64, 500.0_f64, 7000.0_f64);
        let h1 = pt.enthalpy(550.0_f64);
        let h2 = pt.enthalpy(700.0_f64);
        assert!(h2 > h1, "enthalpy should increase with temperature");
    }

    #[test]
    fn test_phase_latent_heat_jump() {
        let l = 250_000.0_f64;
        let pt = PhaseThermal::new(500.0_f64, 600.0_f64, l, 500.0_f64, 7000.0_f64);
        let h_below = pt.enthalpy(499.0_f64);
        let h_above = pt.enthalpy(601.0_f64);
        let sensible = 500.0_f64 * (601.0_f64 - 499.0_f64);
        let total = h_above - h_below;
        // Total should exceed sensible by approximately L
        assert!(total > sensible, "enthalpy jump should include latent heat");
    }

    // ── 1-D rod steady-state ──────────────────────────────────────────────────

    #[test]
    fn test_1d_rod_steady_linear_distribution() {
        // Linear temp distribution expected
        let temps = solve_1d_rod_steady(4, 1.0_f64, 10.0_f64, 100.0_f64, 200.0_f64);
        assert_eq!(temps.len(), 5);
        assert!((temps[0] - 100.0_f64).abs() < 1e-10);
        assert!((temps[4] - 200.0_f64).abs() < 1e-10);
        // Interior should be linear: 125, 150, 175
        assert!((temps[1] - 125.0_f64).abs() < 1e-8, "t[1]={}", temps[1]);
        assert!((temps[2] - 150.0_f64).abs() < 1e-8, "t[2]={}", temps[2]);
        assert!((temps[3] - 175.0_f64).abs() < 1e-8, "t[3]={}", temps[3]);
    }

    #[test]
    fn test_1d_rod_steady_equal_temps() {
        // If both ends are at the same temperature, all nodes should be equal
        let temps = solve_1d_rod_steady(3, 1.0_f64, 5.0_f64, 300.0_f64, 300.0_f64);
        for &t in &temps {
            assert!((t - 300.0_f64).abs() < 1e-8);
        }
    }

    #[test]
    fn test_1d_rod_steady_single_element() {
        let temps = solve_1d_rod_steady(1, 1.0_f64, 1.0_f64, 0.0_f64, 100.0_f64);
        assert_eq!(temps.len(), 2);
        assert!((temps[0] - 0.0_f64).abs() < 1e-10);
        assert!((temps[1] - 100.0_f64).abs() < 1e-10);
    }

    // ── Lumped capacity assembly ──────────────────────────────────────────────

    #[test]
    fn test_lumped_capacity_positive_diagonal() {
        let params = ThermalFemParams::isotropic(1.0_f64, 500.0_f64, 7800.0_f64);
        let nodes = vec![
            ThermalNode::free([0.0, 0.0, 0.0], 0.0_f64),
            ThermalNode::free([1.0, 0.0, 0.0], 0.0_f64),
            ThermalNode::free([0.0, 1.0, 0.0], 0.0_f64),
            ThermalNode::free([0.0, 0.0, 1.0], 0.0_f64),
        ];
        let elements = vec![ThermalElement::isotropic([0, 1, 2, 3], 1.0_f64)];
        let c = assemble_lumped_capacity(&elements, &nodes, &params);
        for (i, row) in c.iter().enumerate() {
            assert!(
                row[i] > 0.0_f64,
                "lumped capacity diagonal must be positive"
            );
        }
    }

    #[test]
    fn test_lumped_capacity_size() {
        let params = ThermalFemParams::isotropic(1.0_f64, 500.0_f64, 7800.0_f64);
        let nodes: Vec<ThermalNode> = (0..6)
            .map(|i| ThermalNode::free([i as f64, 0.0_f64, 0.0_f64], 0.0_f64))
            .collect();
        let elements: Vec<ThermalElement> = vec![];
        let c = assemble_lumped_capacity(&elements, &nodes, &params);
        assert_eq!(c.len(), 6);
        assert_eq!(c[0].len(), 6);
    }

    #[test]
    fn test_transient_solver_steady_state_preserved() {
        // At steady state K*T = f → solver should leave T unchanged (explicit)
        let solver = TransientThermalSolver::explicit_euler(0.001_f64);
        let k = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let m = vec![vec![1.0_f64, 0.0_f64], vec![0.0_f64, 1.0_f64]];
        let f = vec![5.0_f64, 10.0_f64];
        let t0 = vec![5.0_f64, 10.0_f64];
        let t1 = solver.step(&t0, &k, &m, &f);
        assert!((t1[0] - t0[0]).abs() < 1e-10, "should be at steady state");
        assert!((t1[1] - t0[1]).abs() < 1e-10, "should be at steady state");
    }
}
