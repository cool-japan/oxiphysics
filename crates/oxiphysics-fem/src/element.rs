// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Finite element formulations: linear tetrahedron, trilinear hexahedron,
//! serendipity elements, Jacobian computation, element quality metrics,
//! and stress recovery (nodal averaging).

use oxiphysics_core::math::Vec3;

/// A linear tetrahedral element (4-node, constant strain).
///
/// This is the simplest 3D solid element, producing constant stress
/// and strain within each element.
pub struct LinearTetrahedron;

/// Type alias for [`LinearTetrahedron`].
///
/// Provided as a more descriptive name for the 4-node tetrahedral element
/// used throughout the FEM pipeline.
pub type TetrahedralElement = LinearTetrahedron;

impl LinearTetrahedron {
    /// Compute the volume of a tetrahedron from its four node coordinates.
    ///
    /// Uses the formula `V = |det(J)| / 6` where J is formed from edge vectors.
    pub fn volume(nodes: &[Vec3; 4]) -> f64 {
        let x10 = nodes[1] - nodes[0];
        let x20 = nodes[2] - nodes[0];
        let x30 = nodes[3] - nodes[0];
        // Volume = |det([x10, x20, x30])| / 6
        let det = x10.x * (x20.y * x30.z - x20.z * x30.y) - x10.y * (x20.x * x30.z - x20.z * x30.x)
            + x10.z * (x20.x * x30.y - x20.y * x30.x);
        det.abs() / 6.0
    }

    /// Compute the B-matrix (strain-displacement) for a linear tetrahedron.
    ///
    /// Returns a 6x12 matrix in row-major order. The B-matrix relates nodal
    /// displacements to element strains in Voigt notation.
    ///
    /// Also returns the signed volume (can be negative if nodes are ordered
    /// such that det < 0).
    pub fn b_matrix(nodes: &[Vec3; 4]) -> ([f64; 72], f64) {
        // Edge vectors
        let x10 = nodes[1] - nodes[0];
        let x20 = nodes[2] - nodes[0];
        let x30 = nodes[3] - nodes[0];

        // Jacobian determinant (= 6 * signed_volume)
        let det_j = x10.x * (x20.y * x30.z - x20.z * x30.y)
            - x10.y * (x20.x * x30.z - x20.z * x30.x)
            + x10.z * (x20.x * x30.y - x20.y * x30.x);

        let vol = det_j / 6.0;
        let inv_6v = 1.0 / det_j; // = 1 / (6V)

        // Cofactors of J (columns = [x10, x20, x30]).
        //   a_ij = cof(J)[i-1, j-1] = (-1)^(i+j) * det( minor removing row i-1,
        //   col j-1 ).
        // Note the natural-number 1-based indexing: `a11` is the (1,1) entry.
        let a11 = x20.y * x30.z - x20.z * x30.y;
        let a12 = -(x20.x * x30.z - x20.z * x30.x);
        let a13 = x20.x * x30.y - x20.y * x30.x;

        let a21 = -(x10.y * x30.z - x10.z * x30.y);
        let a22 = x10.x * x30.z - x10.z * x30.x;
        let a23 = -(x10.x * x30.y - x10.y * x30.x);

        let a31 = x10.y * x20.z - x10.z * x20.y;
        let a32 = -(x10.x * x20.z - x10.z * x20.x);
        let a33 = x10.x * x20.y - x10.y * x20.x;

        // Shape-function gradients: dN_k/dx_l = J^{-1}[k-1, l-1]
        //   J^{-1} = adj(J) / det(J) = cof(J)^T / det(J),
        // so dN_k/dx_l = a_{l,k} / det(J). The `a_ij` variable naming above
        // follows matrix convention cof(J)[i-1, j-1]; transposing gives
        // J^{-1}[k-1, l-1] = a_{l,k}/det, i.e. the ROW of J^{-1} for node k
        // is `(a_{1,k}, a_{2,k}, a_{3,k}) / det`. For a symmetric Jacobian
        // (e.g. the axis-aligned unit tet) adj(J) = cof(J) and the
        // transpose is silent, which is why this bug stayed hidden on
        // axis-aligned test cases until diagnosed on a sheared tet.
        let dn1dx = a11 * inv_6v;
        let dn1dy = a12 * inv_6v;
        let dn1dz = a13 * inv_6v;

        let dn2dx = a21 * inv_6v;
        let dn2dy = a22 * inv_6v;
        let dn2dz = a23 * inv_6v;

        let dn3dx = a31 * inv_6v;
        let dn3dy = a32 * inv_6v;
        let dn3dz = a33 * inv_6v;

        // dN0/dx = -(dN1/dx + dN2/dx + dN3/dx)
        let dn0dx = -(dn1dx + dn2dx + dn3dx);
        let dn0dy = -(dn1dy + dn2dy + dn3dy);
        let dn0dz = -(dn1dz + dn2dz + dn3dz);

        let dnx = [dn0dx, dn1dx, dn2dx, dn3dx];
        let dny = [dn0dy, dn1dy, dn2dy, dn3dy];
        let dnz = [dn0dz, dn1dz, dn2dz, dn3dz];

        // Build B matrix (6 x 12), stored in row-major order
        // Voigt order: [eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]
        let mut b = [0.0f64; 72]; // 6 rows x 12 cols

        for i in 0..4 {
            let col = i * 3;
            // Row 0: eps_xx = dN/dx * u_x
            b[col] = dnx[i];
            // Row 1: eps_yy = dN/dy * u_y
            b[12 + col + 1] = dny[i];
            // Row 2: eps_zz = dN/dz * u_z
            b[2 * 12 + col + 2] = dnz[i];
            // Row 3: gamma_xy = dN/dy * u_x + dN/dx * u_y
            b[3 * 12 + col] = dny[i];
            b[3 * 12 + col + 1] = dnx[i];
            // Row 4: gamma_yz = dN/dz * u_y + dN/dy * u_z
            b[4 * 12 + col + 1] = dnz[i];
            b[4 * 12 + col + 2] = dny[i];
            // Row 5: gamma_xz = dN/dz * u_x + dN/dx * u_z
            b[5 * 12 + col] = dnz[i];
            b[5 * 12 + col + 2] = dnx[i];
        }

        (b, vol)
    }

    /// Compute the 12x12 element stiffness matrix for a linear tetrahedron.
    ///
    /// The stiffness matrix is `K_e = V * B^T D B`, where `V` is the element
    /// volume, `B` is the strain-displacement matrix, and `D` is the 6x6
    /// constitutive matrix.
    ///
    /// # Arguments
    ///
    /// * `nodes` - The four node coordinates
    /// * `d_matrix` - The 6x6 constitutive (elasticity) matrix in Voigt notation
    pub fn element_stiffness(nodes: &[Vec3; 4], d_matrix: &[[f64; 6]; 6]) -> [[f64; 12]; 12] {
        let (b, vol) = Self::b_matrix(nodes);
        let vol_abs = vol.abs();

        // Compute DB = D * B (6x12)
        let mut db = [0.0f64; 72];
        for i in 0..6 {
            for j in 0..12 {
                let mut sum = 0.0;
                for k in 0..6 {
                    sum += d_matrix[i][k] * b[k * 12 + j];
                }
                db[i * 12 + j] = sum;
            }
        }

        // Compute K = V * B^T * DB (12x12)
        let mut ke = [[0.0f64; 12]; 12];
        for i in 0..12 {
            for j in 0..12 {
                let mut sum = 0.0;
                for k in 0..6 {
                    sum += b[k * 12 + i] * db[k * 12 + j];
                }
                ke[i][j] = vol_abs * sum;
            }
        }

        ke
    }

    /// Compute the shape function values at a point given by barycentric
    /// coordinates `(xi, eta, zeta)` with `N_0 = 1 - xi - eta - zeta`.
    ///
    /// Returns `[N0, N1, N2, N3]`.
    pub fn shape_functions(xi: f64, eta: f64, zeta: f64) -> [f64; 4] {
        [1.0 - xi - eta - zeta, xi, eta, zeta]
    }

    /// Assemble the 12-component body-force load vector for a single element.
    ///
    /// Uses the consistent formulation:
    /// `f_e = V / 4 * [b_x, b_y, b_z, b_x, b_y, b_z, ...]`
    /// (lumped: each of the 4 nodes receives `V/4` of the total body force).
    ///
    /// # Arguments
    ///
    /// * `nodes`       – The four node coordinates.
    /// * `body_force`  – Body force per unit volume `[b_x, b_y, b_z]`.
    pub fn load_vector(nodes: &[Vec3; 4], body_force: [f64; 3]) -> [f64; 12] {
        let vol = Self::volume(nodes);
        let f_node = vol / 4.0;
        let mut f = [0.0f64; 12];
        for i in 0..4 {
            f[i * 3] = f_node * body_force[0];
            f[i * 3 + 1] = f_node * body_force[1];
            f[i * 3 + 2] = f_node * body_force[2];
        }
        f
    }
}

// ── Trilinear hexahedral element (8-node brick) ─────────────────────────────

/// 8-node trilinear hexahedral (brick) element.
///
/// Natural coordinates `(xi, eta, zeta)` are in `[-1, 1]^3`.
/// Node ordering follows the standard convention:
/// ```text
///   7----6
///  /|   /|
/// 4----5 |
/// | 3--|-2
/// |/   |/
/// 0----1
/// ```
#[allow(dead_code)]
pub struct HexahedralElement;

#[allow(dead_code)]
impl HexahedralElement {
    /// Evaluate the 8 shape functions at natural coordinates `(xi, eta, zeta)`.
    ///
    /// Each shape function is `N_i = 1/8 * (1 + xi_i*xi) * (1 + eta_i*eta) * (1 + zeta_i*zeta)`.
    pub fn shape_functions(xi: f64, eta: f64, zeta: f64) -> [f64; 8] {
        // Node natural coordinates: (xi_i, eta_i, zeta_i)
        let nodes = [
            (-1.0, -1.0, -1.0), // 0
            (1.0, -1.0, -1.0),  // 1
            (1.0, 1.0, -1.0),   // 2
            (-1.0, 1.0, -1.0),  // 3
            (-1.0, -1.0, 1.0),  // 4
            (1.0, -1.0, 1.0),   // 5
            (1.0, 1.0, 1.0),    // 6
            (-1.0, 1.0, 1.0),   // 7
        ];
        let mut n = [0.0; 8];
        for (i, &(xi_i, eta_i, zeta_i)) in nodes.iter().enumerate() {
            n[i] = 0.125 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta);
        }
        n
    }

    /// Compute shape function derivatives w.r.t. natural coordinates.
    ///
    /// Returns `(dN_dxi, dN_deta, dN_dzeta)` each of length 8.
    pub fn shape_function_derivatives(
        xi: f64,
        eta: f64,
        zeta: f64,
    ) -> ([f64; 8], [f64; 8], [f64; 8]) {
        let nodes = [
            (-1.0, -1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, 1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0),
        ];
        let mut dn_dxi = [0.0; 8];
        let mut dn_deta = [0.0; 8];
        let mut dn_dzeta = [0.0; 8];
        for (i, &(xi_i, eta_i, zeta_i)) in nodes.iter().enumerate() {
            dn_dxi[i] = 0.125 * xi_i * (1.0 + eta_i * eta) * (1.0 + zeta_i * zeta);
            dn_deta[i] = 0.125 * (1.0 + xi_i * xi) * eta_i * (1.0 + zeta_i * zeta);
            dn_dzeta[i] = 0.125 * (1.0 + xi_i * xi) * (1.0 + eta_i * eta) * zeta_i;
        }
        (dn_dxi, dn_deta, dn_dzeta)
    }

    /// Compute the 3x3 Jacobian matrix at a given point in natural coordinates.
    ///
    /// `nodes` are the 8 physical node coordinates as `[[x, y, z\]; 8]`.
    /// Returns the 3x3 Jacobian `J[i][j] = dx_i / d_xi_j` in row-major.
    pub fn jacobian(nodes: &[[f64; 3]; 8], xi: f64, eta: f64, zeta: f64) -> [[f64; 3]; 3] {
        let (dn_dxi, dn_deta, dn_dzeta) = Self::shape_function_derivatives(xi, eta, zeta);
        let mut j = [[0.0; 3]; 3];
        for n in 0..8 {
            for d in 0..3 {
                j[0][d] += dn_dxi[n] * nodes[n][d];
                j[1][d] += dn_deta[n] * nodes[n][d];
                j[2][d] += dn_dzeta[n] * nodes[n][d];
            }
        }
        j
    }

    /// Compute the determinant of a 3x3 matrix.
    pub fn det3x3(m: &[[f64; 3]; 3]) -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    /// Compute the inverse of a 3x3 matrix.
    ///
    /// Returns `None` if the determinant is near zero.
    pub fn inv3x3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
        let det = Self::det3x3(m);
        if det.abs() < 1e-30 {
            return None;
        }
        let inv_det = 1.0 / det;
        Some([
            [
                (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det,
                (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det,
                (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det,
            ],
            [
                (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det,
                (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det,
                (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det,
            ],
            [
                (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det,
                (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det,
                (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det,
            ],
        ])
    }

    /// Compute the volume of a hexahedral element using 2x2x2 Gauss quadrature.
    pub fn volume(nodes: &[[f64; 3]; 8]) -> f64 {
        let gp = [-1.0 / 3.0_f64.sqrt(), 1.0 / 3.0_f64.sqrt()];
        let gw = [1.0, 1.0];
        let mut vol = 0.0;
        for (i, &xi) in gp.iter().enumerate() {
            for (j, &eta) in gp.iter().enumerate() {
                for (k, &zeta) in gp.iter().enumerate() {
                    let jac = Self::jacobian(nodes, xi, eta, zeta);
                    let det_j = Self::det3x3(&jac);
                    vol += gw[i] * gw[j] * gw[k] * det_j.abs();
                }
            }
        }
        vol
    }
}

// ── Serendipity element (20-node hexahedron) ────────────────────────────────

/// 20-node serendipity hexahedral element.
///
/// This element has 8 corner nodes and 12 mid-edge nodes, providing
/// quadratic interpolation along the edges without interior nodes.
#[allow(dead_code)]
pub struct SerendipityElement;

#[allow(dead_code)]
impl SerendipityElement {
    /// Evaluate the 20 shape functions at natural coordinates `(xi, eta, zeta)`.
    ///
    /// The 8 corner nodes are numbered 0-7 (same as `HexahedralElement`),
    /// followed by 12 mid-edge nodes 8-19.
    pub fn shape_functions(xi: f64, eta: f64, zeta: f64) -> [f64; 20] {
        // Corner node natural coordinates
        let corners: [(f64, f64, f64); 8] = [
            (-1.0, -1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, 1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0),
        ];

        // Mid-edge node natural coordinates
        let midedges: [(f64, f64, f64); 12] = [
            (0.0, -1.0, -1.0), // 8:  edge 0-1
            (1.0, 0.0, -1.0),  // 9:  edge 1-2
            (0.0, 1.0, -1.0),  // 10: edge 2-3
            (-1.0, 0.0, -1.0), // 11: edge 3-0
            (0.0, -1.0, 1.0),  // 12: edge 4-5
            (1.0, 0.0, 1.0),   // 13: edge 5-6
            (0.0, 1.0, 1.0),   // 14: edge 6-7
            (-1.0, 0.0, 1.0),  // 15: edge 7-4
            (-1.0, -1.0, 0.0), // 16: edge 0-4
            (1.0, -1.0, 0.0),  // 17: edge 1-5
            (1.0, 1.0, 0.0),   // 18: edge 2-6
            (-1.0, 1.0, 0.0),  // 19: edge 3-7
        ];

        let mut n = [0.0; 20];

        // Corner nodes: N_i = 1/8 * (1+xi_i*xi)(1+eta_i*eta)(1+zeta_i*zeta)(xi_i*xi+eta_i*eta+zeta_i*zeta-2)
        for (i, &(xi_i, eta_i, zeta_i)) in corners.iter().enumerate() {
            let f1 = 1.0 + xi_i * xi;
            let f2 = 1.0 + eta_i * eta;
            let f3 = 1.0 + zeta_i * zeta;
            n[i] = 0.125 * f1 * f2 * f3 * (xi_i * xi + eta_i * eta + zeta_i * zeta - 2.0);
        }

        // Mid-edge nodes
        for (idx, &(xi_m, eta_m, zeta_m)) in midedges.iter().enumerate() {
            let i = idx + 8;
            if xi_m.abs() < 1e-14 {
                // xi = 0 mid-edge node
                n[i] = 0.25 * (1.0 - xi * xi) * (1.0 + eta_m * eta) * (1.0 + zeta_m * zeta);
            } else if eta_m.abs() < 1e-14 {
                // eta = 0 mid-edge node
                n[i] = 0.25 * (1.0 + xi_m * xi) * (1.0 - eta * eta) * (1.0 + zeta_m * zeta);
            } else {
                // zeta = 0 mid-edge node
                n[i] = 0.25 * (1.0 + xi_m * xi) * (1.0 + eta_m * eta) * (1.0 - zeta * zeta);
            }
        }

        n
    }

    /// Check partition of unity: shape functions must sum to 1.
    pub fn verify_partition_of_unity(xi: f64, eta: f64, zeta: f64) -> f64 {
        let n = Self::shape_functions(xi, eta, zeta);
        let sum: f64 = n.iter().sum();
        (sum - 1.0).abs()
    }
}

// ── Jacobian computation utilities ──────────────────────────────────────────

/// Compute the condition number of a 3x3 Jacobian matrix
/// (ratio of largest to smallest singular value, approximated via
/// the Frobenius norm and determinant).
///
/// A condition number close to 1 indicates a well-shaped element.
/// Large condition numbers indicate poor element quality (distorted).
#[allow(dead_code)]
pub fn jacobian_condition_number(j: &[[f64; 3]; 3]) -> f64 {
    // Frobenius norm
    let mut frob2 = 0.0;
    for row in j {
        for &val in row {
            frob2 += val * val;
        }
    }
    let frob = frob2.sqrt();

    let det = HexahedralElement::det3x3(j);
    if det.abs() < 1e-30 {
        return f64::MAX;
    }

    // Approximate condition number: ||J||_F * ||J^{-1}||_F
    // Using ||J^{-1}||_F ≈ ||adj(J)||_F / |det(J)|
    let inv = match HexahedralElement::inv3x3(j) {
        Some(inv) => inv,
        None => return f64::MAX,
    };
    let mut inv_frob2 = 0.0;
    for row in &inv {
        for &val in row {
            inv_frob2 += val * val;
        }
    }
    let inv_frob = inv_frob2.sqrt();

    frob * inv_frob / 3.0 // Normalize by dimension
}

// ── Element quality metrics ─────────────────────────────────────────────────

/// Compute the aspect ratio of a tetrahedral element.
///
/// The aspect ratio is defined as the ratio of the longest edge
/// to the shortest edge. A value of 1.0 indicates a regular (equilateral)
/// tetrahedron. Higher values indicate more elongated elements.
#[allow(dead_code)]
pub fn tet_aspect_ratio(nodes: &[[f64; 3]; 4]) -> f64 {
    let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

    let mut min_len = f64::MAX;
    let mut max_len = 0.0_f64;

    for &(i, j) in &edges {
        let dx = nodes[j][0] - nodes[i][0];
        let dy = nodes[j][1] - nodes[i][1];
        let dz = nodes[j][2] - nodes[i][2];
        let len = (dx * dx + dy * dy + dz * dz).sqrt();
        if len < min_len {
            min_len = len;
        }
        if len > max_len {
            max_len = len;
        }
    }

    if min_len < 1e-30 {
        return f64::MAX;
    }

    max_len / min_len
}

/// Compute the skewness of a tetrahedral element.
///
/// Skewness is defined as `1 - V / V_ideal` where `V_ideal` is the
/// volume of an equilateral tetrahedron with the same average edge length.
/// A value of 0 is ideal; values approaching 1 indicate degenerate elements.
#[allow(dead_code)]
pub fn tet_skewness(nodes: &[[f64; 3]; 4]) -> f64 {
    let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];

    let mut sum_len = 0.0;
    for &(i, j) in &edges {
        let dx = nodes[j][0] - nodes[i][0];
        let dy = nodes[j][1] - nodes[i][1];
        let dz = nodes[j][2] - nodes[i][2];
        sum_len += (dx * dx + dy * dy + dz * dz).sqrt();
    }
    let avg_len = sum_len / 6.0;

    // Volume of an equilateral tetrahedron with edge length `a` is
    // V_eq = a^3 / (6 * sqrt(2))
    let v_ideal = avg_len * avg_len * avg_len / (6.0 * 2.0_f64.sqrt());

    // Compute actual volume
    let x10 = [
        nodes[1][0] - nodes[0][0],
        nodes[1][1] - nodes[0][1],
        nodes[1][2] - nodes[0][2],
    ];
    let x20 = [
        nodes[2][0] - nodes[0][0],
        nodes[2][1] - nodes[0][1],
        nodes[2][2] - nodes[0][2],
    ];
    let x30 = [
        nodes[3][0] - nodes[0][0],
        nodes[3][1] - nodes[0][1],
        nodes[3][2] - nodes[0][2],
    ];
    let det = x10[0] * (x20[1] * x30[2] - x20[2] * x30[1])
        - x10[1] * (x20[0] * x30[2] - x20[2] * x30[0])
        + x10[2] * (x20[0] * x30[1] - x20[1] * x30[0]);
    let vol = det.abs() / 6.0;

    if v_ideal < 1e-30 {
        return 1.0;
    }

    (1.0 - vol / v_ideal).clamp(0.0, 1.0)
}

/// Compute the radius ratio quality metric for a tetrahedron.
///
/// This is defined as `3 * r_in / r_circ` where `r_in` is the inradius
/// and `r_circ` is the circumradius. The value is 1 for a regular
/// tetrahedron and approaches 0 for degenerate elements.
#[allow(dead_code)]
pub fn tet_radius_ratio(nodes: &[[f64; 3]; 4]) -> f64 {
    // Compute volume
    let x10 = [
        nodes[1][0] - nodes[0][0],
        nodes[1][1] - nodes[0][1],
        nodes[1][2] - nodes[0][2],
    ];
    let x20 = [
        nodes[2][0] - nodes[0][0],
        nodes[2][1] - nodes[0][1],
        nodes[2][2] - nodes[0][2],
    ];
    let x30 = [
        nodes[3][0] - nodes[0][0],
        nodes[3][1] - nodes[0][1],
        nodes[3][2] - nodes[0][2],
    ];
    let det = x10[0] * (x20[1] * x30[2] - x20[2] * x30[1])
        - x10[1] * (x20[0] * x30[2] - x20[2] * x30[0])
        + x10[2] * (x20[0] * x30[1] - x20[1] * x30[0]);
    let vol = det.abs() / 6.0;

    if vol < 1e-30 {
        return 0.0;
    }

    // Compute face areas (4 faces)
    let face_area = |a: &[f64; 3], b: &[f64; 3], c: &[f64; 3]| -> f64 {
        let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cross = [
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        ];
        0.5 * (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt()
    };

    let a0 = face_area(&nodes[0], &nodes[1], &nodes[2]);
    let a1 = face_area(&nodes[0], &nodes[1], &nodes[3]);
    let a2 = face_area(&nodes[0], &nodes[2], &nodes[3]);
    let a3 = face_area(&nodes[1], &nodes[2], &nodes[3]);

    let total_area = a0 + a1 + a2 + a3;

    // Inradius: r_in = 3 * V / total_surface_area
    let r_in = 3.0 * vol / total_area;

    // Compute edge lengths
    let edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
    let mut edge_lens = [0.0; 6];
    for (idx, &(i, j)) in edges.iter().enumerate() {
        let dx = nodes[j][0] - nodes[i][0];
        let dy = nodes[j][1] - nodes[i][1];
        let dz = nodes[j][2] - nodes[i][2];
        edge_lens[idx] = (dx * dx + dy * dy + dz * dz).sqrt();
    }

    // Circumradius formula: R = (abc) / (8 * V * area_factor)
    // Use simplified approximation: product of edge lengths around a face
    // divided by (8 * V) for one face, then take the maximum
    // Simplified: use the longest edge as estimate of 2*R
    let max_edge = edge_lens.iter().cloned().fold(0.0_f64, f64::max);
    let r_circ_approx = max_edge / 2.0;

    if r_circ_approx < 1e-30 {
        return 0.0;
    }

    (3.0 * r_in / r_circ_approx).min(1.0)
}

// ── Stress recovery (nodal averaging) ───────────────────────────────────────

/// Perform nodal stress averaging (stress recovery).
///
/// Given element stresses (one stress tensor per element) and an element
/// connectivity table, compute averaged nodal stresses by distributing
/// each element's stress to its nodes and averaging.
///
/// # Arguments
///
/// * `num_nodes` - Total number of nodes in the mesh
/// * `elements` - Element connectivity: each element is `[n0, n1, n2, n3]`
/// * `element_stresses` - Stress tensor for each element in Voigt notation
///   `[sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_xz]`
///
/// # Returns
///
/// Nodal stresses: `Vec` of length `num_nodes`, each entry is `[f64; 6]`.
#[allow(dead_code)]
pub fn nodal_stress_averaging(
    num_nodes: usize,
    elements: &[[usize; 4]],
    element_stresses: &[[f64; 6]],
) -> Vec<[f64; 6]> {
    let mut nodal_stress = vec![[0.0_f64; 6]; num_nodes];
    let mut nodal_count = vec![0usize; num_nodes];

    for (elem, stress) in elements.iter().zip(element_stresses.iter()) {
        for &node in elem {
            for c in 0..6 {
                nodal_stress[node][c] += stress[c];
            }
            nodal_count[node] += 1;
        }
    }

    for i in 0..num_nodes {
        if nodal_count[i] > 0 {
            let inv = 1.0 / nodal_count[i] as f64;
            for s in &mut nodal_stress[i] {
                *s *= inv;
            }
        }
    }

    nodal_stress
}

/// Compute the von Mises equivalent stress from a Voigt stress tensor.
///
/// `sigma = [sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_xz]`
///
/// Formula: `sigma_vm = sqrt(0.5 * ((s1-s2)^2 + (s2-s3)^2 + (s3-s1)^2) + 3*(t12^2+t23^2+t13^2))`
#[allow(dead_code)]
pub fn von_mises_stress(sigma: &[f64; 6]) -> f64 {
    let s = sigma;
    let vm2 = 0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2))
        + 3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
    vm2.sqrt()
}

/// Compute the hydrostatic (mean) stress from a Voigt stress tensor.
#[allow(dead_code)]
pub fn hydrostatic_stress(sigma: &[f64; 6]) -> f64 {
    (sigma[0] + sigma[1] + sigma[2]) / 3.0
}

/// Compute deviatoric stress tensor from a Voigt stress tensor.
///
/// Returns the deviatoric part: `s_ij = sigma_ij - sigma_h * delta_ij`
#[allow(dead_code)]
pub fn deviatoric_stress(sigma: &[f64; 6]) -> [f64; 6] {
    let p = hydrostatic_stress(sigma);
    [
        sigma[0] - p,
        sigma[1] - p,
        sigma[2] - p,
        sigma[3],
        sigma[4],
        sigma[5],
    ]
}

/// Extrapolate integration point values to nodes for a tetrahedral element
/// using the shape function extrapolation matrix.
///
/// For a linear tetrahedron with a single integration point (centroid),
/// the extrapolated nodal values are all equal to the centroid value.
#[allow(dead_code)]
pub fn extrapolate_to_nodes_tet(centroid_value: &[f64; 6]) -> [[f64; 6]; 4] {
    [*centroid_value; 4]
}

/// Superconvergent Patch Recovery (SPR) for a patch of elements sharing a node.
///
/// Given a set of element centroid stresses and their centroid coordinates,
/// fits a linear polynomial to the stress field and evaluates it at the
/// target node position.
///
/// # Arguments
///
/// * `target` - Physical coordinates of the target node
/// * `centroids` - Centroid coordinates of surrounding elements
/// * `stresses` - Centroid stresses `[sig_xx, sig_yy, sig_zz, tau_xy, tau_yz, tau_xz]`
#[allow(dead_code)]
pub fn spr_recovery(target: [f64; 3], centroids: &[[f64; 3]], stresses: &[[f64; 6]]) -> [f64; 6] {
    let n = centroids.len();
    if n == 0 {
        return [0.0; 6];
    }
    if n == 1 {
        return stresses[0];
    }

    // For each stress component, fit a linear polynomial:
    // sigma(x,y,z) = a0 + a1*x + a2*y + a3*z
    // using least squares.
    let mut result = [0.0; 6];
    for comp in 0..6 {
        // Build normal equations: A^T A * coeffs = A^T b
        // A is n x 4 matrix: [1, x, y, z]
        let mut ata = [[0.0; 4]; 4];
        let mut atb = [0.0; 4];

        for i in 0..n {
            let row = [1.0, centroids[i][0], centroids[i][1], centroids[i][2]];
            let b_val = stresses[i][comp];
            for r in 0..4 {
                for c in 0..4 {
                    ata[r][c] += row[r] * row[c];
                }
                atb[r] += row[r] * b_val;
            }
        }

        // Solve 4x4 system using Gaussian elimination with partial pivoting
        let coeffs = solve_4x4(&ata, &atb);

        // Evaluate at target
        result[comp] =
            coeffs[0] + coeffs[1] * target[0] + coeffs[2] * target[1] + coeffs[3] * target[2];
    }

    result
}

/// Solve a 4x4 linear system using Gaussian elimination with partial pivoting.
#[allow(dead_code)]
#[allow(clippy::needless_range_loop)]
fn solve_4x4(a: &[[f64; 4]; 4], b: &[f64; 4]) -> [f64; 4] {
    let mut m = [[0.0; 5]; 4]; // augmented matrix
    for i in 0..4 {
        for j in 0..4 {
            m[i][j] = a[i][j];
        }
        m[i][4] = b[i];
    }

    // Forward elimination with partial pivoting
    for col in 0..4 {
        // Find pivot
        let mut max_val = m[col][col].abs();
        let mut max_row = col;
        for row in (col + 1)..4 {
            if m[row][col].abs() > max_val {
                max_val = m[row][col].abs();
                max_row = row;
            }
        }
        if max_val < 1e-30 {
            return [0.0; 4]; // Singular
        }
        // Swap rows
        if max_row != col {
            m.swap(col, max_row);
        }
        // Eliminate
        let pivot = m[col][col];
        for row in (col + 1)..4 {
            let factor = m[row][col] / pivot;
            for c in col..5 {
                m[row][c] -= factor * m[col][c];
            }
        }
    }

    // Back substitution
    let mut x = [0.0; 4];
    for i in (0..4).rev() {
        let mut sum = m[i][4];
        for j in (i + 1)..4 {
            sum -= m[i][j] * x[j];
        }
        x[i] = sum / m[i][i];
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_element_volume() {
        // Unit tetrahedron with vertices at origin and unit vectors
        let nodes = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let vol = LinearTetrahedron::volume(&nodes);
        assert!((vol - 1.0 / 6.0).abs() < 1e-12, "volume = {vol}");
    }

    #[test]
    #[allow(clippy::needless_range_loop)]
    fn test_element_stiffness_symmetry() {
        let nodes = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        // Simple isotropic D for testing
        let e = 200.0e9;
        let nu = 0.3;
        let d = crate::constitutive::LinearElasticMaterial::new(e, nu).constitutive_matrix();

        let ke = LinearTetrahedron::element_stiffness(&nodes, &d);

        // Check symmetry
        for i in 0..12 {
            for j in 0..12 {
                let diff = (ke[i][j] - ke[j][i]).abs();
                let scale = ke[i][j].abs().max(ke[j][i].abs()).max(1.0);
                assert!(
                    diff / scale < 1e-10,
                    "Stiffness not symmetric at ({i},{j}): {} vs {}",
                    ke[i][j],
                    ke[j][i]
                );
            }
        }
    }

    /// Shape functions must sum to 1 at any interior point.
    #[test]
    fn test_shape_functions_partition_of_unity() {
        let test_points = [
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.25, 0.25, 0.25),
            (0.1, 0.2, 0.3),
        ];
        for (xi, eta, zeta) in test_points {
            let n = LinearTetrahedron::shape_functions(xi, eta, zeta);
            let sum: f64 = n.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-14,
                "Shape functions do not sum to 1 at ({xi},{eta},{zeta}): sum = {sum}"
            );
        }
    }

    /// Shape function N_i evaluates to 1 at node i and 0 at other nodes.
    #[test]
    fn test_shape_functions_nodal_values() {
        let n0 = LinearTetrahedron::shape_functions(0.0, 0.0, 0.0);
        assert!((n0[0] - 1.0).abs() < 1e-14);
        assert!(n0[1].abs() < 1e-14);
        assert!(n0[2].abs() < 1e-14);
        assert!(n0[3].abs() < 1e-14);

        let n1 = LinearTetrahedron::shape_functions(1.0, 0.0, 0.0);
        assert!(n1[0].abs() < 1e-14);
        assert!((n1[1] - 1.0).abs() < 1e-14);

        let n3 = LinearTetrahedron::shape_functions(0.0, 0.0, 1.0);
        assert!(n3[0].abs() < 1e-14);
        assert!((n3[3] - 1.0).abs() < 1e-14);
    }

    /// The load vector must equal the body force times the element volume in total.
    #[test]
    fn test_load_vector_body_force() {
        let nodes = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let body_force = [0.0, -9.81, 0.0]; // gravity-like
        let f = LinearTetrahedron::load_vector(&nodes, body_force);

        let vol = LinearTetrahedron::volume(&nodes); // = 1/6
        let total_fy: f64 = f.iter().skip(1).step_by(3).sum();
        let expected_fy = body_force[1] * vol;
        assert!(
            (total_fy - expected_fy).abs() < 1e-12,
            "total Fy = {total_fy}, expected {expected_fy}"
        );
        // All x and z components should be zero
        let total_fx: f64 = f.iter().step_by(3).sum();
        assert!(
            total_fx.abs() < 1e-12,
            "total Fx should be zero, got {total_fx}"
        );
    }

    /// Patch test: uniform uniaxial strain should be reproduced exactly.
    #[test]
    fn test_single_element_patch_test_uniform_strain() {
        let eps = 0.001;

        let nodes = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];

        let u_e = [
            0.0 * eps,
            0.0,
            0.0,
            1.0 * eps,
            0.0,
            0.0,
            0.0 * eps,
            0.0,
            0.0,
            0.0 * eps,
            0.0,
            0.0,
        ];

        let (b, _vol) = LinearTetrahedron::b_matrix(&nodes);

        let mut strain = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..12 {
                strain[i] += b[i * 12 + j] * u_e[j];
            }
        }

        assert!(
            (strain[0] - eps).abs() < 1e-12,
            "eps_xx = {}, expected {eps}",
            strain[0]
        );
        assert!(
            strain[1].abs() < 1e-12,
            "eps_yy should be 0, got {}",
            strain[1]
        );
        assert!(
            strain[2].abs() < 1e-12,
            "eps_zz should be 0, got {}",
            strain[2]
        );
        assert!(
            strain[3].abs() < 1e-12,
            "gamma_xy should be 0, got {}",
            strain[3]
        );
        assert!(
            strain[4].abs() < 1e-12,
            "gamma_yz should be 0, got {}",
            strain[4]
        );
        assert!(
            strain[5].abs() < 1e-12,
            "gamma_xz should be 0, got {}",
            strain[5]
        );
    }

    /// Rigid body mode: zero internal forces under uniform translation.
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn test_element_stiffness_rigid_body_zero_force() {
        let nodes = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.0, 0.0, 2.0),
        ];
        let d = crate::constitutive::LinearElasticMaterial::new(1.0e6, 0.25).constitutive_matrix();
        let ke = LinearTetrahedron::element_stiffness(&nodes, &d);

        let u_rb = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
        for i in 0..12 {
            let mut fi = 0.0;
            for j in 0..12 {
                fi += ke[i][j] * u_rb[j];
            }
            assert!(fi.abs() < 1e-6, "K*u_rb[{i}] = {fi:.3e}, should be ~0");
        }
    }

    /// `TetrahedralElement` type alias must work identically to `LinearTetrahedron`.
    #[test]
    fn test_tetrahedral_element_alias() {
        let nodes = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        let vol_alias = TetrahedralElement::volume(&nodes);
        let vol_direct = LinearTetrahedron::volume(&nodes);
        assert!((vol_alias - vol_direct).abs() < 1e-14);
    }

    // ── Hexahedral element tests ────────────────────────────────────────────

    /// Hex shape functions must sum to 1 at any point (partition of unity).
    #[test]
    fn test_hex_shape_functions_partition_of_unity() {
        let test_points = [
            (0.0, 0.0, 0.0),
            (1.0, -1.0, 0.5),
            (-0.5, 0.3, -0.7),
            (-1.0, -1.0, -1.0),
            (1.0, 1.0, 1.0),
        ];
        for (xi, eta, zeta) in test_points {
            let n = HexahedralElement::shape_functions(xi, eta, zeta);
            let sum: f64 = n.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-14,
                "Hex shape functions do not sum to 1 at ({xi},{eta},{zeta}): sum = {sum}"
            );
        }
    }

    /// At each corner node, only one shape function should be 1.
    #[test]
    fn test_hex_shape_functions_nodal_values() {
        let corners = [
            (-1.0, -1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, 1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0),
        ];
        for (node_idx, &(xi, eta, zeta)) in corners.iter().enumerate() {
            let n = HexahedralElement::shape_functions(xi, eta, zeta);
            for (i, &val) in n.iter().enumerate() {
                if i == node_idx {
                    assert!(
                        (val - 1.0).abs() < 1e-14,
                        "N_{i} at node {node_idx} should be 1, got {val}"
                    );
                } else {
                    assert!(
                        val.abs() < 1e-14,
                        "N_{i} at node {node_idx} should be 0, got {val}"
                    );
                }
            }
        }
    }

    /// Jacobian of a unit cube should be diagonal with determinant = 1/8.
    #[test]
    fn test_hex_jacobian_unit_cube() {
        // Unit cube [0,1]^3 mapped from [-1,1]^3
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let jac = HexahedralElement::jacobian(&nodes, 0.0, 0.0, 0.0);
        let det = HexahedralElement::det3x3(&jac);
        // Mapping from [-1,1]^3 to [0,1]^3: dx/dxi = 0.5
        // So det(J) = (0.5)^3 = 0.125
        assert!(
            (det - 0.125).abs() < 1e-12,
            "det(J) = {det}, expected 0.125"
        );
    }

    /// Volume of a unit cube should be 1.
    #[test]
    fn test_hex_volume_unit_cube() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ];
        let vol = HexahedralElement::volume(&nodes);
        assert!(
            (vol - 1.0).abs() < 1e-12,
            "unit cube volume = {vol}, expected 1.0"
        );
    }

    /// Volume of a 2x3x4 box should be 24.
    #[test]
    fn test_hex_volume_stretched_box() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [2.0, 3.0, 0.0],
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
            [2.0, 0.0, 4.0],
            [2.0, 3.0, 4.0],
            [0.0, 3.0, 4.0],
        ];
        let vol = HexahedralElement::volume(&nodes);
        assert!(
            (vol - 24.0).abs() < 1e-10,
            "box volume = {vol}, expected 24.0"
        );
    }

    /// Inverse of a 3x3 identity should be identity.
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn test_inv3x3_identity() {
        let id = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let inv = HexahedralElement::inv3x3(&id).unwrap();
        for i in 0..3 {
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!(
                    (inv[i][j] - expected).abs() < 1e-14,
                    "inv[{i}][{j}] = {}, expected {expected}",
                    inv[i][j]
                );
            }
        }
    }

    // ── Serendipity element tests ───────────────────────────────────────────

    /// Serendipity shape functions must sum to 1 at any point.
    #[test]
    fn test_serendipity_partition_of_unity() {
        let pts = [
            (0.0, 0.0, 0.0),
            (0.5, -0.3, 0.7),
            (-1.0, -1.0, -1.0),
            (1.0, 1.0, 1.0),
            (0.0, -1.0, -1.0), // mid-edge node
        ];
        for (xi, eta, zeta) in pts {
            let err = SerendipityElement::verify_partition_of_unity(xi, eta, zeta);
            assert!(
                err < 1e-12,
                "Serendipity partition of unity error = {err} at ({xi},{eta},{zeta})"
            );
        }
    }

    /// At corner nodes, only the corresponding corner shape function is 1.
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn test_serendipity_corner_nodes() {
        let corners = [
            (-1.0, -1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, 1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0),
        ];
        for (node_idx, &(xi, eta, zeta)) in corners.iter().enumerate() {
            let n = SerendipityElement::shape_functions(xi, eta, zeta);
            assert!(
                (n[node_idx] - 1.0).abs() < 1e-12,
                "Corner node {node_idx}: N = {}, expected 1.0",
                n[node_idx]
            );
            // All mid-edge functions should be 0 at corner nodes
            for i in 8..20 {
                assert!(
                    n[i].abs() < 1e-12,
                    "Mid-edge N_{i} at corner {node_idx} = {}, expected 0",
                    n[i]
                );
            }
        }
    }

    // ── Element quality tests ───────────────────────────────────────────────

    /// A regular tetrahedron should have aspect ratio close to 1.
    #[test]
    fn test_tet_aspect_ratio_regular() {
        // Regular tetrahedron with unit edge length
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
        ];
        let ar = tet_aspect_ratio(&nodes);
        assert!(
            (ar - 1.0).abs() < 0.05,
            "regular tet aspect ratio = {ar}, expected ~1.0"
        );
    }

    /// An elongated tetrahedron should have high aspect ratio.
    #[test]
    fn test_tet_aspect_ratio_elongated() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, 0.0, 0.1],
        ];
        let ar = tet_aspect_ratio(&nodes);
        assert!(
            ar > 5.0,
            "elongated tet should have high aspect ratio, got {ar}"
        );
    }

    /// A regular tetrahedron should have low skewness.
    #[test]
    fn test_tet_skewness_regular() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
        ];
        let sk = tet_skewness(&nodes);
        assert!(sk < 0.1, "regular tet skewness should be near 0, got {sk}");
    }

    /// Radius ratio for a regular tetrahedron should be close to 1.
    #[test]
    fn test_tet_radius_ratio_regular() {
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 2.0, 0.0],
            [0.5, 3.0_f64.sqrt() / 6.0, (2.0_f64 / 3.0).sqrt()],
        ];
        let rr = tet_radius_ratio(&nodes);
        assert!(
            rr > 0.5,
            "regular tet radius ratio should be high, got {rr}"
        );
    }

    // ── Jacobian condition number tests ─────────────────────────────────────

    /// Identity Jacobian should have condition number of 1.
    #[test]
    fn test_jacobian_condition_number_identity() {
        let jac = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let cn = jacobian_condition_number(&jac);
        assert!(
            (cn - 1.0).abs() < 1e-10,
            "identity Jacobian condition number = {cn}, expected 1.0"
        );
    }

    /// Scaled identity should still have condition number of 1.
    #[test]
    fn test_jacobian_condition_number_scaled() {
        let jac = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]];
        let cn = jacobian_condition_number(&jac);
        assert!(
            (cn - 1.0).abs() < 1e-10,
            "scaled identity condition number = {cn}, expected 1.0"
        );
    }

    // ── Stress recovery tests ───────────────────────────────────────────────

    /// Nodal averaging for a single element should return the element stress.
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn test_nodal_stress_averaging_single_element() {
        let elements = vec![[0, 1, 2, 3]];
        let stresses = vec![[100.0, 50.0, 30.0, 10.0, 5.0, 3.0]];
        let nodal = nodal_stress_averaging(4, &elements, &stresses);
        for node in 0..4 {
            for c in 0..6 {
                assert!(
                    (nodal[node][c] - stresses[0][c]).abs() < 1e-12,
                    "node {node} comp {c}: {} vs {}",
                    nodal[node][c],
                    stresses[0][c]
                );
            }
        }
    }

    /// Two elements sharing a node should average the stress at that node.
    #[test]
    fn test_nodal_stress_averaging_two_elements() {
        let elements = vec![[0, 1, 2, 3], [1, 4, 5, 6]];
        let stresses = vec![
            [100.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [200.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let nodal = nodal_stress_averaging(7, &elements, &stresses);
        // Node 1 is shared: average of 100 and 200 = 150
        assert!(
            (nodal[1][0] - 150.0).abs() < 1e-12,
            "node 1 sig_xx = {}, expected 150",
            nodal[1][0]
        );
        // Node 0 belongs only to element 0: stress = 100
        assert!(
            (nodal[0][0] - 100.0).abs() < 1e-12,
            "node 0 sig_xx = {}, expected 100",
            nodal[0][0]
        );
    }

    /// Von Mises stress for uniaxial stress state.
    #[test]
    fn test_von_mises_uniaxial() {
        let sigma = [100.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = von_mises_stress(&sigma);
        assert!(
            (vm - 100.0).abs() < 1e-10,
            "uniaxial von Mises = {vm}, expected 100"
        );
    }

    /// Von Mises stress for pure shear.
    #[test]
    fn test_von_mises_pure_shear() {
        let sigma = [0.0, 0.0, 0.0, 50.0, 0.0, 0.0];
        let vm = von_mises_stress(&sigma);
        let expected = (3.0_f64 * 50.0 * 50.0).sqrt();
        assert!(
            (vm - expected).abs() < 1e-10,
            "pure shear von Mises = {vm}, expected {expected}"
        );
    }

    /// Hydrostatic stress for equal triaxial tension.
    #[test]
    fn test_hydrostatic_stress() {
        let sigma = [100.0, 200.0, 300.0, 0.0, 0.0, 0.0];
        let p = hydrostatic_stress(&sigma);
        assert!((p - 200.0).abs() < 1e-12, "hydrostatic = {p}, expected 200");
    }

    /// Deviatoric stress should have zero trace.
    #[test]
    fn test_deviatoric_zero_trace() {
        let sigma = [100.0, 200.0, 300.0, 10.0, 20.0, 30.0];
        let dev = deviatoric_stress(&sigma);
        let trace = dev[0] + dev[1] + dev[2];
        assert!(
            trace.abs() < 1e-12,
            "deviatoric trace = {trace}, should be 0"
        );
        // Shear components unchanged
        assert!((dev[3] - 10.0).abs() < 1e-12);
        assert!((dev[4] - 20.0).abs() < 1e-12);
        assert!((dev[5] - 30.0).abs() < 1e-12);
    }

    /// SPR recovery with single element returns the element stress.
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn test_spr_recovery_single() {
        let target = [0.0, 0.0, 0.0];
        let centroids = vec![[0.25, 0.25, 0.25]];
        let stresses = vec![[100.0, 50.0, 30.0, 10.0, 5.0, 3.0]];
        let recovered = spr_recovery(target, &centroids, &stresses);
        // With single point, any polynomial fit is exact at that point
        // but extrapolated to target; result depends on linear fit
        // Just verify it returns a finite result
        for c in 0..6 {
            assert!(recovered[c].is_finite(), "comp {c} is not finite");
        }
    }

    /// SPR recovery with constant stress field should recover the exact stress everywhere.
    #[test]
    fn test_spr_recovery_constant_field() {
        let target = [1.0, 2.0, 3.0];
        let centroids = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
        ];
        let constant = [42.0, 21.0, 10.0, 5.0, 2.0, 1.0];
        let stresses = vec![constant; 5];
        let recovered = spr_recovery(target, &centroids, &stresses);
        for c in 0..6 {
            assert!(
                (recovered[c] - constant[c]).abs() < 1e-8,
                "comp {c}: recovered {} vs expected {}",
                recovered[c],
                constant[c]
            );
        }
    }

    /// Extrapolation to nodes for tet should give uniform value.
    #[allow(clippy::needless_range_loop)]
    #[test]
    fn test_extrapolate_to_nodes_tet() {
        let val = [100.0, 50.0, 30.0, 10.0, 5.0, 3.0];
        let nodal = extrapolate_to_nodes_tet(&val);
        for i in 0..4 {
            for c in 0..6 {
                assert!((nodal[i][c] - val[c]).abs() < 1e-14);
            }
        }
    }
}
