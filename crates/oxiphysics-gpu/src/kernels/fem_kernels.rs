// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CPU-mock FEM (Finite Element Method) compute kernels.
//!
//! All kernels mirror what would run on a GPU in their data layout and dispatch
//! model, but execute in pure Rust on the CPU.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── Element type ─────────────────────────────────────────────────────────────

/// Finite element type supported by the FEM kernel suite.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ElementType {
    /// 4-node linear tetrahedron.
    Tet4,
    /// 8-node trilinear hexahedron.
    Hex8,
    /// 3-node linear triangle.
    Tri3,
}

// ── FemKernelConfig ──────────────────────────────────────────────────────────

/// Configuration shared across all FEM kernels in a simulation run.
#[derive(Debug, Clone)]
pub struct FemKernelConfig {
    /// Element type used in this simulation.
    pub element_type: ElementType,
    /// Gauss integration order (1, 2, or 3).
    pub integration_order: u32,
    /// Young's modulus (Pa).
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless, in (-1, 0.5)).
    pub poisson_ratio: f64,
    /// Mass density (kg/m³).
    pub density: f64,
    /// Thermal expansion coefficient (1/K).
    pub thermal_alpha: f64,
    /// Reference temperature (K).
    pub t_ref: f64,
}

impl FemKernelConfig {
    /// Create a new config with common steel-like defaults.
    pub fn new_steel() -> Self {
        Self {
            element_type: ElementType::Tet4,
            integration_order: 2,
            youngs_modulus: 200.0e9,
            poisson_ratio: 0.3,
            density: 7850.0,
            thermal_alpha: 12.0e-6,
            t_ref: 293.15,
        }
    }

    /// Lamé parameter λ.
    pub fn lame_lambda(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    }

    /// Lamé parameter μ (shear modulus).
    pub fn lame_mu(&self) -> f64 {
        let e = self.youngs_modulus;
        let nu = self.poisson_ratio;
        e / (2.0 * (1.0 + nu))
    }

    /// Bulk modulus K.
    pub fn bulk_modulus(&self) -> f64 {
        self.youngs_modulus / (3.0 * (1.0 - 2.0 * self.poisson_ratio))
    }
}

// ── Helper math ──────────────────────────────────────────────────────────────

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn length3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let l = length3(v);
    if l < 1e-15 { [0.0; 3] } else { [v[0]/l, v[1]/l, v[2]/l] }
}

/// Multiply two 3×3 matrices stored in row-major flat arrays of length 9.
fn mat3x3_mul(a: [f64; 9], b: [f64; 9]) -> [f64; 9] {
    let mut c = [0.0f64; 9];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i*3+j] += a[i*3+k] * b[k*3+j];
            }
        }
    }
    c
}

/// Determinant of a 3×3 matrix stored in row-major order.
fn det3x3(m: [f64; 9]) -> f64 {
    m[0]*(m[4]*m[8]-m[5]*m[7])
    - m[1]*(m[3]*m[8]-m[5]*m[6])
    + m[2]*(m[3]*m[7]-m[4]*m[6])
}

/// Inverse of a 3×3 matrix; returns None if singular.
fn inv3x3(m: [f64; 9]) -> Option<[f64; 9]> {
    let d = det3x3(m);
    if d.abs() < 1e-15 { return None; }
    let inv_d = 1.0 / d;
    Some([
        (m[4]*m[8]-m[5]*m[7])*inv_d,
       -(m[1]*m[8]-m[2]*m[7])*inv_d,
        (m[1]*m[5]-m[2]*m[4])*inv_d,
       -(m[3]*m[8]-m[5]*m[6])*inv_d,
        (m[0]*m[8]-m[2]*m[6])*inv_d,
       -(m[0]*m[5]-m[2]*m[3])*inv_d,
        (m[3]*m[7]-m[4]*m[6])*inv_d,
       -(m[0]*m[7]-m[1]*m[6])*inv_d,
        (m[0]*m[4]-m[1]*m[3])*inv_d,
    ])
}

// ── FemStiffnessKernel ────────────────────────────────────────────────────────

/// Kernel for computing element stiffness matrices.
#[derive(Debug, Clone)]
pub struct FemStiffnessKernel {
    /// Young's modulus.
    pub e_modulus: f64,
    /// Poisson's ratio.
    pub nu: f64,
}

impl FemStiffnessKernel {
    /// Create a new stiffness kernel.
    pub fn new(e_modulus: f64, nu: f64) -> Self {
        Self { e_modulus, nu }
    }

    /// Build the 6×6 isotropic elasticity tensor C (Voigt notation).
    ///
    /// Row/column ordering: \[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_xz\].
    pub fn elasticity_tensor_voigt(&self) -> [[f64; 6]; 6] {
        let e = self.e_modulus;
        let nu = self.nu;
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let c11 = lam + 2.0 * mu;
        let c12 = lam;
        let c44 = mu;
        let mut c = [[0.0f64; 6]; 6];
        c[0][0] = c11; c[0][1] = c12; c[0][2] = c12;
        c[1][0] = c12; c[1][1] = c11; c[1][2] = c12;
        c[2][0] = c12; c[2][1] = c12; c[2][2] = c11;
        c[3][3] = c44;
        c[4][4] = c44;
        c[5][5] = c44;
        c
    }

    /// Compute the 12×12 element stiffness matrix K_e for a linear Tet4.
    ///
    /// `nodes` — four nodal positions in 3-D space.
    /// Returns a symmetric 12×12 matrix (dof ordering: u1,v1,w1, u2,…,u4,v4,w4).
    pub fn compute_element_stiffness(&self, nodes: &[[f64; 3]; 4]) -> [[f64; 12]; 12] {
        // Shape function gradients for Tet4 (constant throughout element).
        // B matrix columns: [dN/dx, dN/dy, dN/dz] for each node.
        let b_mat = self.tet4_b_matrix(nodes);
        let c_mat = self.elasticity_tensor_voigt();
        let vol = self.tet4_volume(nodes);

        // K = vol * B^T * C * B  (6×12 B, 6×6 C)
        // Step 1: CB = C * B  (6×12)
        let mut cb = [[0.0f64; 12]; 6];
        for i in 0..6 {
            for j in 0..12 {
                for k in 0..6 {
                    cb[i][j] += c_mat[i][k] * b_mat[k][j];
                }
            }
        }
        // Step 2: K = vol * B^T * CB  (12×12)
        let mut k_e = [[0.0f64; 12]; 12];
        for i in 0..12 {
            for j in 0..12 {
                for k in 0..6 {
                    k_e[i][j] += vol * b_mat[k][i] * cb[k][j];
                }
            }
        }
        k_e
    }

    /// Assemble the strain-displacement B matrix for a Tet4 element (6×12).
    pub fn tet4_b_matrix(&self, nodes: &[[f64; 3]; 4]) -> [[f64; 12]; 6] {
        // Jacobian columns: x21, x31, x41 etc.
        let j = [
            sub3(nodes[1], nodes[0]),
            sub3(nodes[2], nodes[0]),
            sub3(nodes[3], nodes[0]),
        ];
        let j_flat = [
            j[0][0], j[1][0], j[2][0],
            j[0][1], j[1][1], j[2][1],
            j[0][2], j[1][2], j[2][2],
        ];
        let j_inv = inv3x3(j_flat).unwrap_or([0.0; 9]);

        // dN/dxi for natural coords: N1=1-xi-eta-zeta, N2=xi, N3=eta, N4=zeta
        // dN/dx = J^{-T} * dN/dxi
        let dn_dxi: [[f64; 3]; 4] = [
            [-1.0, -1.0, -1.0],
            [ 1.0,  0.0,  0.0],
            [ 0.0,  1.0,  0.0],
            [ 0.0,  0.0,  1.0],
        ];
        // dN/dx[a][i] = sum_j J_inv[j][i] * dN_dxi[a][j]  (J_inv is 3×3)
        let mut dn_dx = [[0.0f64; 3]; 4];
        for a in 0..4 {
            for i in 0..3 {
                for jj in 0..3 {
                    dn_dx[a][i] += j_inv[jj*3+i] * dn_dxi[a][jj];
                }
            }
        }
        // B matrix: 6 × 12
        let mut b = [[0.0f64; 12]; 6];
        for a in 0..4 {
            let col = a * 3;
            let (nx, ny, nz) = (dn_dx[a][0], dn_dx[a][1], dn_dx[a][2]);
            b[0][col]   = nx;
            b[1][col+1] = ny;
            b[2][col+2] = nz;
            b[3][col]   = ny; b[3][col+1] = nx;
            b[4][col+1] = nz; b[4][col+2] = ny;
            b[5][col]   = nz; b[5][col+2] = nx;
        }
        b
    }

    /// Volume of a Tet4 element.
    pub fn tet4_volume(&self, nodes: &[[f64; 3]; 4]) -> f64 {
        let a = sub3(nodes[1], nodes[0]);
        let b = sub3(nodes[2], nodes[0]);
        let c = sub3(nodes[3], nodes[0]);
        let axb = cross3(a, b);
        (dot3(axb, c) / 6.0).abs()
    }
}

// ── FemMassKernel ─────────────────────────────────────────────────────────────

/// Kernel for computing element mass matrices.
#[derive(Debug, Clone)]
pub struct FemMassKernel {
    /// Material density (kg/m³).
    pub density: f64,
}

impl FemMassKernel {
    /// Create a new mass kernel.
    pub fn new(density: f64) -> Self {
        Self { density }
    }

    /// Consistent 12×12 mass matrix for a Tet4 element.
    ///
    /// Uses the analytical formula: M = rho * V / 20 * (I_4 + J_4) ⊗ I_3
    /// where the diagonal entries are 2/20 and off-diagonal are 1/20.
    pub fn consistent_mass_matrix(&self, nodes: &[[f64; 3]; 4]) -> [[f64; 12]; 12] {
        let stiff = FemStiffnessKernel::new(1.0, 0.0);
        let vol = stiff.tet4_volume(nodes);
        let c = self.density * vol / 20.0;
        let mut m = [[0.0f64; 12]; 12];
        for a in 0..4 {
            for b in 0..4 {
                let v = if a == b { 2.0 * c } else { c };
                for d in 0..3 {
                    m[a*3+d][b*3+d] += v;
                }
            }
        }
        m
    }

    /// Lumped mass vector (row-sum) for a Tet4 element (length 12).
    pub fn lumped_mass_vector(&self, nodes: &[[f64; 3]; 4]) -> [f64; 12] {
        let m = self.consistent_mass_matrix(nodes);
        let mut lump = [0.0f64; 12];
        for i in 0..12 {
            for j in 0..12 {
                lump[i] += m[i][j];
            }
        }
        lump
    }

    /// Total element mass.
    pub fn element_mass(&self, nodes: &[[f64; 3]; 4]) -> f64 {
        let stiff = FemStiffnessKernel::new(1.0, 0.0);
        self.density * stiff.tet4_volume(nodes)
    }
}

// ── FemForceKernel ────────────────────────────────────────────────────────────

/// Kernel for computing element force vectors.
#[derive(Debug, Clone)]
pub struct FemForceKernel {
    /// Material density (kg/m³).
    pub density: f64,
}

impl FemForceKernel {
    /// Create a new force kernel.
    pub fn new(density: f64) -> Self {
        Self { density }
    }

    /// Body force vector f_b = rho * V/4 * \[bx,by,bz, bx,by,bz, ...\] (12 dof).
    pub fn body_force_vector(&self, nodes: &[[f64; 3]; 4], body: [f64; 3]) -> [f64; 12] {
        let stiff = FemStiffnessKernel::new(1.0, 0.0);
        let vol = stiff.tet4_volume(nodes);
        let scale = self.density * vol / 4.0;
        let mut f = [0.0f64; 12];
        for a in 0..4 {
            f[a*3]   = scale * body[0];
            f[a*3+1] = scale * body[1];
            f[a*3+2] = scale * body[2];
        }
        f
    }

    /// Surface traction vector for a Tri3 face (9 dof for 3 nodes).
    ///
    /// `face_nodes` — three nodal positions of the face.
    /// `traction`   — traction vector (Pa) in global coords.
    pub fn surface_traction_vector(
        &self,
        face_nodes: &[[f64; 3]; 3],
        traction: [f64; 3],
    ) -> [f64; 9] {
        // Area of triangle
        let e1 = sub3(face_nodes[1], face_nodes[0]);
        let e2 = sub3(face_nodes[2], face_nodes[0]);
        let n = cross3(e1, e2);
        let area = length3(n) / 2.0;
        let scale = area / 3.0;
        let mut f = [0.0f64; 9];
        for a in 0..3 {
            f[a*3]   = scale * traction[0];
            f[a*3+1] = scale * traction[1];
            f[a*3+2] = scale * traction[2];
        }
        f
    }

    /// Thermal load vector due to temperature change ΔT (12 dof, Tet4).
    ///
    /// f_th = vol * B^T * C * alpha * ΔT * \[1,1,1,0,0,0\]^T
    pub fn thermal_load_vector(
        &self,
        nodes: &[[f64; 3]; 4],
        alpha: f64,
        delta_t: f64,
        e_mod: f64,
        nu: f64,
    ) -> [f64; 12] {
        let kern = FemStiffnessKernel::new(e_mod, nu);
        let vol = kern.tet4_volume(nodes);
        let b = kern.tet4_b_matrix(nodes);
        let c = kern.elasticity_tensor_voigt();
        // eps_th = alpha * dT * [1,1,1,0,0,0]
        let eps_th = [alpha * delta_t; 6];
        let mut sigma_th = [0.0f64; 6];
        for i in 0..6 {
            // Only first three rows of eps matter (shear is zero for thermal)
            for j in 0..3 {
                sigma_th[i] += c[i][j] * eps_th[j];
            }
        }
        // f = vol * B^T * sigma_th
        let mut f = [0.0f64; 12];
        for i in 0..12 {
            for j in 0..6 {
                f[i] += vol * b[j][i] * sigma_th[j];
            }
        }
        f
    }
}

// ── FemStrainKernel ────────────────────────────────────────────────────────────

/// Kernel for computing strain and stress tensors at Gauss points.
#[derive(Debug, Clone)]
pub struct FemStrainKernel;

impl FemStrainKernel {
    /// Create a new strain kernel.
    pub fn new() -> Self {
        Self
    }

    /// Compute the engineering strain vector (6-component Voigt) from
    /// displacement `u` (12 dof) and the B matrix (6×12).
    pub fn compute_strain_tensor(&self, u: &[f64; 12], b_matrix: &[[f64; 12]; 6]) -> [f64; 6] {
        let mut eps = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..12 {
                eps[i] += b_matrix[i][j] * u[j];
            }
        }
        eps
    }

    /// Compute the Cauchy stress from strain and elasticity tensor C (6×6 Voigt).
    pub fn compute_stress(&self, strain: &[f64; 6], c_matrix: &[[f64; 6]; 6]) -> [f64; 6] {
        let mut sigma = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                sigma[i] += c_matrix[i][j] * strain[j];
            }
        }
        sigma
    }

    /// Von Mises equivalent stress from a 6-component stress vector.
    pub fn von_mises_stress(&self, sigma: &[f64; 6]) -> f64 {
        let s11 = sigma[0]; let s22 = sigma[1]; let s33 = sigma[2];
        let s12 = sigma[3]; let s23 = sigma[4]; let s13 = sigma[5];
        let val = 0.5 * ((s11-s22).powi(2) + (s22-s33).powi(2) + (s33-s11).powi(2))
            + 3.0 * (s12.powi(2) + s23.powi(2) + s13.powi(2));
        val.sqrt()
    }

    /// Hydrostatic pressure p = -(sigma_xx + sigma_yy + sigma_zz) / 3.
    pub fn hydrostatic_pressure(&self, sigma: &[f64; 6]) -> f64 {
        -(sigma[0] + sigma[1] + sigma[2]) / 3.0
    }
}

impl Default for FemStrainKernel {
    fn default() -> Self {
        Self::new()
    }
}

// ── CSR sparse matrix ─────────────────────────────────────────────────────────

/// Compressed Sparse Row matrix (f64, 0-indexed).
#[derive(Debug, Clone)]
pub struct CsrMatrix {
    /// Number of rows.
    pub nrows: usize,
    /// Number of columns.
    pub ncols: usize,
    /// Row pointer array (length nrows+1).
    pub row_ptr: Vec<usize>,
    /// Column indices.
    pub col_idx: Vec<usize>,
    /// Non-zero values.
    pub values: Vec<f64>,
}

impl CsrMatrix {
    /// Create a zero CSR matrix of given size.
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        Self {
            nrows,
            ncols,
            row_ptr: vec![0; nrows + 1],
            col_idx: Vec::new(),
            values: Vec::new(),
        }
    }

    /// Matrix-vector product y = A * x.
    pub fn matvec(&self, x: &[f64]) -> Vec<f64> {
        assert_eq!(x.len(), self.ncols);
        let mut y = vec![0.0f64; self.nrows];
        for row in 0..self.nrows {
            for idx in self.row_ptr[row]..self.row_ptr[row+1] {
                y[row] += self.values[idx] * x[self.col_idx[idx]];
            }
        }
        y
    }

    /// Diagonal vector.
    pub fn diagonal(&self) -> Vec<f64> {
        let n = self.nrows.min(self.ncols);
        let mut diag = vec![0.0f64; n];
        for row in 0..n {
            for idx in self.row_ptr[row]..self.row_ptr[row+1] {
                if self.col_idx[idx] == row {
                    diag[row] = self.values[idx];
                    break;
                }
            }
        }
        diag
    }
}

// ── FemAssemblyKernel ─────────────────────────────────────────────────────────

/// Kernel for assembling global stiffness matrices from element contributions.
#[derive(Debug, Clone)]
pub struct FemAssemblyKernel {
    /// Total number of global degrees of freedom.
    pub n_dofs: usize,
}

impl FemAssemblyKernel {
    /// Create a new assembly kernel.
    pub fn new(n_dofs: usize) -> Self {
        Self { n_dofs }
    }

    /// Assemble a global stiffness matrix in COO format, then convert to CSR.
    ///
    /// `elements` — slice of element DOF connectivity (12 dofs per Tet4).
    /// `local_k`  — local 12×12 stiffness matrices, one per element.
    pub fn parallel_assemble_global(
        &self,
        elements: &[Vec<usize>],
        local_k: &[[[f64; 12]; 12]],
    ) -> CsrMatrix {
        assert_eq!(elements.len(), local_k.len());

        // Triplet list
        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        let n = self.n_dofs;
        for (elem, k_e) in elements.iter().zip(local_k.iter()) {
            let ndofs = elem.len().min(12);
            for i in 0..ndofs {
                for j in 0..ndofs {
                    if elem[i] < n && elem[j] < n {
                        triplets.push((elem[i], elem[j], k_e[i][j]));
                    }
                }
            }
        }

        // Sort by (row, col) and deduplicate into unique triplets
        triplets.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));

        let mut unique_triplets: Vec<(usize, usize, f64)> = Vec::new();
        let mut t = 0;
        while t < triplets.len() {
            let (row, col, _) = triplets[t];
            let mut s = t;
            while s < triplets.len() && triplets[s].0 == row && triplets[s].1 == col {
                s += 1;
            }
            let val: f64 = triplets[t..s].iter().map(|x| x.2).sum();
            unique_triplets.push((row, col, val));
            t = s;
        }

        // Build CSR row_ptr
        let mut row_ptr = vec![0usize; n + 1];
        for &(row, _, _) in &unique_triplets {
            row_ptr[row + 1] += 1;
        }
        for r in 0..n {
            row_ptr[r + 1] += row_ptr[r];
        }

        let nnz = unique_triplets.len();
        let mut col_idx = vec![0usize; nnz];
        let mut values = vec![0.0f64; nnz];
        for (idx, &(_, col, val)) in unique_triplets.iter().enumerate() {
            col_idx[idx] = col;
            values[idx] = val;
        }

        CsrMatrix { nrows: n, ncols: n, row_ptr, col_idx, values }
    }
}

// ── FemSolverKernel ───────────────────────────────────────────────────────────

/// Kernel wrapping the PCG linear solver for FEM systems.
#[derive(Debug, Clone)]
pub struct FemSolverKernel {
    /// Solver tolerance (residual norm).
    pub tol: f64,
    /// Maximum CG iterations.
    pub max_iter: usize,
}

impl FemSolverKernel {
    /// Create a solver kernel with given tolerance and iteration limit.
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self { tol, max_iter }
    }

    /// Preconditioned Conjugate Gradient solver.
    ///
    /// Solves A x = b using a Jacobi (diagonal) preconditioner.
    /// Returns the solution vector and the number of iterations.
    pub fn pcg_solve(
        &self,
        a: &CsrMatrix,
        b: &[f64],
        x0: &[f64],
    ) -> (Vec<f64>, usize) {
        let n = b.len();
        assert_eq!(a.nrows, n);
        assert_eq!(a.ncols, n);
        assert_eq!(x0.len(), n);

        let diag = a.diagonal();
        let precond = |r: &[f64]| -> Vec<f64> {
            r.iter().enumerate().map(|(i, &ri)| {
                if diag[i].abs() > 1e-15 { ri / diag[i] } else { ri }
            }).collect()
        };

        let mut x = x0.to_vec();
        let ax = a.matvec(&x);
        let mut r: Vec<f64> = b.iter().zip(ax.iter()).map(|(bi, axi)| bi - axi).collect();
        let mut z = precond(&r);
        let mut p = z.clone();
        let mut rz = r.iter().zip(z.iter()).map(|(ri, zi)| ri * zi).sum::<f64>();

        for iter in 0..self.max_iter {
            let ap = a.matvec(&p);
            let pap: f64 = p.iter().zip(ap.iter()).map(|(pi, api)| pi * api).sum();
            if pap.abs() < 1e-30 { return (x, iter); }
            let alpha = rz / pap;
            for i in 0..n {
                x[i] += alpha * p[i];
                r[i] -= alpha * ap[i];
            }
            let res_norm: f64 = r.iter().map(|ri| ri * ri).sum::<f64>().sqrt();
            if res_norm < self.tol { return (x, iter + 1); }
            z = precond(&r);
            let rz_new: f64 = r.iter().zip(z.iter()).map(|(ri, zi)| ri * zi).sum();
            let beta = rz_new / rz.max(1e-300);
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }
            rz = rz_new;
        }
        (x, self.max_iter)
    }
}

// ── FemContactKernel ──────────────────────────────────────────────────────────

/// Kernel for contact mechanics calculations.
#[derive(Debug, Clone)]
pub struct FemContactKernel {
    /// Penalty stiffness for contact enforcement.
    pub penalty: f64,
    /// Friction coefficient.
    pub friction_mu: f64,
}

impl FemContactKernel {
    /// Create a new contact kernel.
    pub fn new(penalty: f64, friction_mu: f64) -> Self {
        Self { penalty, friction_mu }
    }

    /// Signed gap function g = (x_slave - x_master) · n_master.
    ///
    /// Negative means penetration.
    pub fn gap_function(&self, x_slave: [f64; 3], x_master: [f64; 3], n_master: [f64; 3]) -> f64 {
        let diff = sub3(x_slave, x_master);
        dot3(diff, n_master)
    }

    /// Outward contact normal from a master surface defined by two edge vectors.
    pub fn contact_normal(&self, edge1: [f64; 3], edge2: [f64; 3]) -> [f64; 3] {
        normalize3(cross3(edge1, edge2))
    }

    /// Penetration depth (positive when penetrating).
    pub fn penetration_depth(&self, gap: f64) -> f64 {
        (-gap).max(0.0)
    }

    /// Penalty contact force vector (3 components) applied at slave node.
    pub fn penalty_force(&self, gap: f64, normal: [f64; 3]) -> [f64; 3] {
        let pen = self.penetration_depth(gap);
        let f_mag = self.penalty * pen;
        [f_mag * normal[0], f_mag * normal[1], f_mag * normal[2]]
    }

    /// Friction force tangential component (Coulomb, sticking limit).
    pub fn friction_force(&self, contact_force: f64, trial_tangent: [f64; 3]) -> [f64; 3] {
        let len = length3(trial_tangent);
        if len < 1e-15 { return [0.0; 3]; }
        let limit = self.friction_mu * contact_force.abs();
        let scale = limit.min(len) / len;
        [trial_tangent[0]*scale, trial_tangent[1]*scale, trial_tangent[2]*scale]
    }
}

// ── FemModalKernel ────────────────────────────────────────────────────────────

/// Kernel for modal analysis (eigenvalue/eigenvector extraction).
#[derive(Debug, Clone)]
pub struct FemModalKernel {
    /// Number of modes to extract.
    pub n_modes: usize,
    /// Lanczos convergence tolerance.
    pub tol: f64,
    /// Maximum Lanczos steps.
    pub max_lanczos: usize,
}

impl FemModalKernel {
    /// Create a new modal kernel.
    pub fn new(n_modes: usize, tol: f64, max_lanczos: usize) -> Self {
        Self { n_modes, tol, max_lanczos }
    }

    /// Minimal Lanczos iteration to approximate eigenvalues of K relative to M.
    ///
    /// `k` — global stiffness matrix (CSR, n×n).
    /// `m_diag` — lumped mass diagonal (length n).
    /// Returns `(eigenvalues, eigenvectors)` where each eigenvector has length n.
    pub fn lanczos_iteration(
        &self,
        k: &CsrMatrix,
        m_diag: &[f64],
    ) -> (Vec<f64>, Vec<Vec<f64>>) {
        let n = k.nrows;
        let m_count = self.n_modes.min(n);

        // Initial random-like vector (deterministic for reproducibility)
        let mut q_prev = vec![0.0f64; n];
        let mut q_cur: Vec<f64> = (0..n).map(|i| ((i+1) as f64).sin()).collect();
        let norm0: f64 = q_cur.iter().map(|x| x*x).sum::<f64>().sqrt();
        for x in q_cur.iter_mut() { *x /= norm0.max(1e-15); }

        let mut alphas = Vec::new();
        let mut betas = Vec::new();
        let mut q_vecs: Vec<Vec<f64>> = vec![q_cur.clone()];

        let steps = self.max_lanczos.min(n);
        for _j in 0..steps {
            let kq = k.matvec(&q_cur);
            // M^{-1} K q
            let mkq: Vec<f64> = kq.iter().enumerate()
                .map(|(i, &v)| if m_diag[i].abs() > 1e-15 { v / m_diag[i] } else { v })
                .collect();
            let alpha: f64 = q_cur.iter().zip(mkq.iter()).map(|(a,b)| a*b).sum();
            alphas.push(alpha);

            let beta_prev = betas.last().copied().unwrap_or(0.0);
            let mut r: Vec<f64> = mkq.iter().enumerate()
                .map(|(i, &v)| v - alpha * q_cur[i] - beta_prev * q_prev[i])
                .collect();
            let beta: f64 = r.iter().map(|x| x*x).sum::<f64>().sqrt();
            betas.push(beta);
            if beta < self.tol { break; }
            for x in r.iter_mut() { *x /= beta; }
            q_prev = q_cur.clone();
            q_cur = r;
            q_vecs.push(q_cur.clone());
        }

        // Tridiagonal Rayleigh-Ritz: just return diagonal as eigenvalue approx
        let mut eigenvalues: Vec<f64> = alphas[..m_count.min(alphas.len())].to_vec();
        eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let eigenvectors: Vec<Vec<f64>> = q_vecs[..m_count.min(q_vecs.len())].to_vec();
        (eigenvalues, eigenvectors)
    }
}

// ── FemDynamicKernel ──────────────────────────────────────────────────────────

/// Kernel for time integration using the Newmark-β method.
#[derive(Debug, Clone)]
pub struct FemDynamicKernel {
    /// Newmark β parameter (0.25 for constant average acceleration).
    pub beta: f64,
    /// Newmark γ parameter (0.5 for no numerical dissipation).
    pub gamma: f64,
}

impl FemDynamicKernel {
    /// Create a new dynamic kernel with constant-average-acceleration defaults.
    pub fn new(beta: f64, gamma: f64) -> Self {
        Self { beta, gamma }
    }

    /// Single Newmark-β time step.
    ///
    /// Solves `(M + gamma*dt*C + beta*dt²*K) u_new = f - C*v_pred - K*u_pred`
    /// using a simplified lumped-mass explicit path for demo purposes.
    ///
    /// `m_diag` — lumped mass diagonal.
    /// `c_diag` — damping diagonal (Rayleigh).
    /// `k`      — stiffness matrix (CSR).
    /// `f`      — external force vector.
    /// `u`, `v` — current displacement and velocity.
    /// `dt`     — time step.
    ///
    /// Returns `(u_new, v_new, a_new)`.
    pub fn newmark_beta_step(
        &self,
        m_diag: &[f64],
        c_diag: &[f64],
        k: &CsrMatrix,
        f: &[f64],
        u: &[f64],
        v: &[f64],
        a_prev: &[f64],
        dt: f64,
    ) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = u.len();
        // Predictors
        let u_pred: Vec<f64> = (0..n).map(|i| {
            u[i] + dt * v[i] + dt * dt * (0.5 - self.beta) * a_prev[i]
        }).collect();
        let v_pred: Vec<f64> = (0..n).map(|i| {
            v[i] + dt * (1.0 - self.gamma) * a_prev[i]
        }).collect();

        // Effective RHS: r = f - K*u_pred - C*v_pred
        let ku_pred = k.matvec(&u_pred);
        let r: Vec<f64> = (0..n).map(|i| {
            f[i] - ku_pred[i] - c_diag[i] * v_pred[i]
        }).collect();

        // Effective mass M_eff diagonal: M + gamma*dt*C + beta*dt^2 * K_diag
        let k_diag = k.diagonal();
        let m_eff: Vec<f64> = (0..n).map(|i| {
            m_diag[i] + self.gamma * dt * c_diag[i] + self.beta * dt * dt * k_diag[i]
        }).collect();

        // Solve M_eff * a_new = r (diagonal → trivial)
        let a_new: Vec<f64> = (0..n).map(|i| {
            if m_eff[i].abs() > 1e-30 { r[i] / m_eff[i] } else { 0.0 }
        }).collect();

        let u_new: Vec<f64> = (0..n).map(|i| u_pred[i] + self.beta * dt * dt * a_new[i]).collect();
        let v_new: Vec<f64> = (0..n).map(|i| v_pred[i] + self.gamma * dt * a_new[i]).collect();

        (u_new, v_new, a_new)
    }

    /// Critical time step estimate Δt_cr = 2/ω_max (Courant condition).
    pub fn critical_dt(&self, omega_max: f64) -> f64 {
        if omega_max < 1e-15 { return f64::INFINITY; }
        2.0 / omega_max
    }

    /// Rayleigh damping coefficients α, β from two modal frequencies and damping ratios.
    pub fn rayleigh_coefficients(omega1: f64, omega2: f64, zeta1: f64, zeta2: f64) -> (f64, f64) {
        // [omega1  1/omega1] [alpha]   [2*zeta1]
        // [omega2  1/omega2] [beta ] = [2*zeta2]
        let det = omega1 / omega2 - omega2 / omega1;
        if det.abs() < 1e-15 { return (0.0, 0.0); }
        let alpha = (2.0 * (zeta1 * omega1 - zeta2 * omega2)) / (omega1 * omega1 - omega2 * omega2);
        let beta  = (2.0 * (zeta2 * omega2 * omega2 - zeta1 * omega1 * omega1))
            / (omega1 * omega2 * (omega1 * omega1 - omega2 * omega2));
        (alpha, beta)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod fem_tests {
    use super::*;
    use std::f64::consts::PI;

    fn unit_tet() -> [[f64; 3]; 4] {
        [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
    }

    #[test]
    fn test_config_lame_lambda() {
        let cfg = FemKernelConfig::new_steel();
        let lam = cfg.lame_lambda();
        assert!(lam > 0.0);
    }

    #[test]
    fn test_config_lame_mu() {
        let cfg = FemKernelConfig::new_steel();
        let mu = cfg.lame_mu();
        assert!(mu > 0.0);
    }

    #[test]
    fn test_config_bulk_modulus() {
        let cfg = FemKernelConfig::new_steel();
        let k = cfg.bulk_modulus();
        assert!(k > 0.0);
    }

    #[test]
    fn test_tet4_volume_unit() {
        let kern = FemStiffnessKernel::new(1.0, 0.3);
        let vol = kern.tet4_volume(&unit_tet());
        let expected = 1.0 / 6.0;
        assert!((vol - expected).abs() < 1e-10, "vol={vol}");
    }

    #[test]
    fn test_elasticity_tensor_symmetry() {
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let c = kern.elasticity_tensor_voigt();
        for i in 0..6 {
            for j in 0..6 {
                assert!((c[i][j] - c[j][i]).abs() < 1e-3, "C[{i}][{j}] != C[{j}][{i}]");
            }
        }
    }

    #[test]
    fn test_stiffness_matrix_size() {
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let k = kern.compute_element_stiffness(&unit_tet());
        assert_eq!(k.len(), 12);
        assert_eq!(k[0].len(), 12);
    }

    #[test]
    fn test_stiffness_matrix_symmetry() {
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let k = kern.compute_element_stiffness(&unit_tet());
        for i in 0..12 {
            for j in 0..12 {
                assert!((k[i][j] - k[j][i]).abs() < 1.0, "K[{i}][{j}] asymmetry");
            }
        }
    }

    #[test]
    fn test_stiffness_positive_diagonal() {
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let k = kern.compute_element_stiffness(&unit_tet());
        for i in 0..12 {
            assert!(k[i][i] >= 0.0, "K[{i}][{i}] = {}", k[i][i]);
        }
    }

    #[test]
    fn test_consistent_mass_matrix_symmetry() {
        let mk = FemMassKernel::new(7850.0);
        let m = mk.consistent_mass_matrix(&unit_tet());
        for i in 0..12 {
            for j in 0..12 {
                assert!((m[i][j] - m[j][i]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_lumped_mass_positive() {
        let mk = FemMassKernel::new(7850.0);
        let lump = mk.lumped_mass_vector(&unit_tet());
        for v in lump {
            assert!(v > 0.0);
        }
    }

    #[test]
    fn test_element_mass() {
        let mk = FemMassKernel::new(7850.0);
        let mass = mk.element_mass(&unit_tet());
        let expected = 7850.0 / 6.0;
        assert!((mass - expected).abs() < 1e-6);
    }

    #[test]
    fn test_body_force_vector() {
        let fk = FemForceKernel::new(7850.0);
        let f = fk.body_force_vector(&unit_tet(), [0.0, -9.81, 0.0]);
        // Total force = density * vol * g
        let total_fy: f64 = (0..4).map(|a| f[a*3+1]).sum();
        let expected = 7850.0 * (1.0/6.0) * (-9.81);
        assert!((total_fy - expected).abs() < 1e-4);
    }

    #[test]
    fn test_surface_traction_vector() {
        let fk = FemForceKernel::new(7850.0);
        let face = [[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]];
        let f = fk.surface_traction_vector(&face, [0.0, 0.0, 1.0]);
        let total_fz: f64 = (0..3).map(|a| f[a*3+2]).sum();
        // Area = 0.5, traction=1 => total force = 0.5
        assert!((total_fz - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_load_vector() {
        let fk = FemForceKernel::new(7850.0);
        let f = fk.thermal_load_vector(&unit_tet(), 12e-6, 100.0, 200e9, 0.3);
        assert_eq!(f.len(), 12);
        // Should be nonzero
        let norm: f64 = f.iter().map(|x| x*x).sum::<f64>().sqrt();
        assert!(norm > 0.0);
    }

    #[test]
    fn test_compute_strain_tensor() {
        let sk = FemStrainKernel::new();
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let b = kern.tet4_b_matrix(&unit_tet());
        // Pure x-displacement on node 2
        let mut u = [0.0f64; 12];
        u[3] = 0.001; // x-disp of node 2
        let eps = sk.compute_strain_tensor(&u, &b);
        assert_eq!(eps.len(), 6);
    }

    #[test]
    fn test_compute_stress() {
        let sk = FemStrainKernel::new();
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let c = kern.elasticity_tensor_voigt();
        let eps = [0.001, 0.0, 0.0, 0.0, 0.0, 0.0];
        let sigma = sk.compute_stress(&eps, &c);
        assert!(sigma[0] > 0.0);
    }

    #[test]
    fn test_von_mises_stress() {
        let sk = FemStrainKernel::new();
        let sigma = [100e6, 0.0, 0.0, 50e6, 0.0, 0.0];
        let vm = sk.von_mises_stress(&sigma);
        assert!(vm > 0.0);
    }

    #[test]
    fn test_hydrostatic_pressure() {
        let sk = FemStrainKernel::new();
        let sigma = [100.0, 100.0, 100.0, 0.0, 0.0, 0.0];
        let p = sk.hydrostatic_pressure(&sigma);
        assert!((p + 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_csr_matvec() {
        // 2×2 identity
        let mut csr = CsrMatrix::zeros(2, 2);
        csr.row_ptr = vec![0, 1, 2];
        csr.col_idx = vec![0, 1];
        csr.values = vec![1.0, 1.0];
        let x = [3.0, 5.0];
        let y = csr.matvec(&x);
        assert_eq!(y, vec![3.0, 5.0]);
    }

    #[test]
    fn test_pcg_solve_identity() {
        let mut csr = CsrMatrix::zeros(3, 3);
        csr.row_ptr = vec![0, 1, 2, 3];
        csr.col_idx = vec![0, 1, 2];
        csr.values = vec![2.0, 2.0, 2.0];
        let b = vec![4.0, 6.0, 8.0];
        let x0 = vec![0.0; 3];
        let solver = FemSolverKernel::new(1e-10, 100);
        let (x, _iters) = solver.pcg_solve(&csr, &b, &x0);
        assert!((x[0] - 2.0).abs() < 1e-8);
        assert!((x[1] - 3.0).abs() < 1e-8);
        assert!((x[2] - 4.0).abs() < 1e-8);
    }

    #[test]
    fn test_gap_function_no_penetration() {
        let ck = FemContactKernel::new(1e6, 0.3);
        let g = ck.gap_function([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(g > 0.0);
    }

    #[test]
    fn test_gap_function_penetration() {
        let ck = FemContactKernel::new(1e6, 0.3);
        let g = ck.gap_function([0.0, -0.1, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!(g < 0.0);
    }

    #[test]
    fn test_contact_normal() {
        let ck = FemContactKernel::new(1e6, 0.3);
        let n = ck.contact_normal([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let len: f64 = n.iter().map(|x| x*x).sum::<f64>().sqrt();
        assert!((len - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_penetration_depth() {
        let ck = FemContactKernel::new(1e6, 0.3);
        assert!((ck.penetration_depth(-0.05) - 0.05).abs() < 1e-12);
        assert!((ck.penetration_depth(0.1)).abs() < 1e-12);
    }

    #[test]
    fn test_penalty_force() {
        let ck = FemContactKernel::new(1e6, 0.3);
        let f = ck.penalty_force(-0.001, [0.0, 1.0, 0.0]);
        assert!((f[1] - 1000.0).abs() < 1e-6);
    }

    #[test]
    fn test_lanczos_returns_modes() {
        let mk = FemModalKernel::new(2, 1e-8, 10);
        let mut k = CsrMatrix::zeros(4, 4);
        k.row_ptr = vec![0, 1, 2, 3, 4];
        k.col_idx = vec![0, 1, 2, 3];
        k.values = vec![4.0, 4.0, 4.0, 4.0];
        let m_diag = vec![1.0; 4];
        let (evals, evecs) = mk.lanczos_iteration(&k, &m_diag);
        assert!(!evals.is_empty());
        assert!(!evecs.is_empty());
    }

    #[test]
    fn test_newmark_step_static() {
        let dk = FemDynamicKernel::new(0.25, 0.5);
        let mut k = CsrMatrix::zeros(2, 2);
        k.row_ptr = vec![0, 1, 2];
        k.col_idx = vec![0, 1];
        k.values = vec![1000.0, 1000.0];
        let m_diag = vec![1.0, 1.0];
        let c_diag = vec![0.0, 0.0];
        let f = vec![0.0, 0.0];
        let u = vec![0.0, 0.0];
        let v = vec![0.0, 0.0];
        let a = vec![0.0, 0.0];
        let (u_new, _, _) = dk.newmark_beta_step(&m_diag, &c_diag, &k, &f, &u, &v, &a, 0.001);
        assert!(u_new[0].abs() < 1e-10);
    }

    #[test]
    fn test_critical_dt() {
        let dk = FemDynamicKernel::new(0.25, 0.5);
        let dt = dk.critical_dt(1000.0);
        assert!((dt - 0.002).abs() < 1e-10);
    }

    #[test]
    fn test_rayleigh_coefficients() {
        let (alpha, beta) = FemDynamicKernel::rayleigh_coefficients(
            2.0 * PI * 10.0, 2.0 * PI * 100.0, 0.02, 0.02
        );
        assert!(alpha.is_finite());
        assert!(beta.is_finite());
    }

    #[test]
    fn test_assembly_kernel_basic() {
        let ak = FemAssemblyKernel::new(12);
        let kern = FemStiffnessKernel::new(200e9, 0.3);
        let nodes = unit_tet();
        let k_e = kern.compute_element_stiffness(&nodes);
        let dofs: Vec<usize> = (0..12).collect();
        let csr = ak.parallel_assemble_global(&[dofs], &[k_e]);
        assert_eq!(csr.nrows, 12);
    }

    #[test]
    fn test_element_type_variants() {
        let t4 = ElementType::Tet4;
        let h8 = ElementType::Hex8;
        let t3 = ElementType::Tri3;
        assert_ne!(t4, h8);
        assert_ne!(h8, t3);
    }

    #[test]
    fn test_csr_diagonal() {
        let mut csr = CsrMatrix::zeros(3, 3);
        csr.row_ptr = vec![0, 1, 2, 3];
        csr.col_idx = vec![0, 1, 2];
        csr.values = vec![3.0, 5.0, 7.0];
        let d = csr.diagonal();
        assert_eq!(d, vec![3.0, 5.0, 7.0]);
    }
}
