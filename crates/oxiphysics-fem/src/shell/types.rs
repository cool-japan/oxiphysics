//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Vibration analysis of thin circular cylindrical shells using the
/// Donnell simplified equations for flexural modes.
#[derive(Debug, Clone)]
pub struct ShellVibrationSolver {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Mass density \[kg/m³\].
    pub density: f64,
    /// Shell thickness \[m\].
    pub thickness: f64,
    /// Shell radius \[m\].
    pub radius: f64,
    /// Shell length \[m\].
    pub length: f64,
}
impl ShellVibrationSolver {
    /// Create a new shell vibration solver.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        thickness: f64,
        radius: f64,
        length: f64,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            density,
            thickness,
            radius,
            length,
        }
    }
    /// Ring (breathing) frequency for circumferential mode n (n ≥ 2).
    ///
    /// ω_n = (c_L / R) * √(n² (n²−1)² h² / (12 R² (n²+1))) + n²
    /// simplified Donnell approximation for long shells.
    pub fn ring_frequency_mode(&self, n: u32) -> f64 {
        let n = n as f64;
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let rho = self.density;
        let r = self.radius;
        let h = self.thickness;
        let c = e / (rho * (1.0 - nu * nu));
        let omega_sq = c * (h * h) / (12.0 * r * r * r * r) * n * n * (n * n - 1.0) * (n * n - 1.0)
            / (n * n + 1.0);
        omega_sq.abs().sqrt()
    }
    /// Breathing (axisymmetric) frequency, n=0 mode.
    ///
    /// ω_0 = c_L / R  (ring frequency)
    pub fn breathing_frequency(&self) -> f64 {
        let c_l = (self.young_modulus
            / (self.density * (1.0 - self.poisson_ratio * self.poisson_ratio)))
            .sqrt();
        c_l / self.radius
    }
    /// Longitudinal natural frequency for axial mode m.
    ///
    /// ω_m = m * π / L * c_L
    pub fn longitudinal_frequency(&self, m: u32) -> f64 {
        let c_l = (self.young_modulus / self.density).sqrt();
        m as f64 * std::f64::consts::PI / self.length * c_l
    }
}
/// Doubly-curved shell geometry (principal radii of curvature).
pub struct ShellGeometry {
    /// First principal radius of curvature.
    pub r1: f64,
    /// Second principal radius of curvature.
    pub r2: f64,
    /// Shell thickness.
    pub thickness: f64,
}
impl ShellGeometry {
    /// Create a shell geometry.
    pub fn new(r1: f64, r2: f64, thickness: f64) -> Self {
        Self { r1, r2, thickness }
    }
}
/// Classical Lamination Theory (CLT) stacking sequence.
///
/// Computes the A (extensional), B (coupling), D (bending) sub-matrices
/// of the full 6×6 ABD stiffness matrix.
#[derive(Debug, Clone)]
pub struct CompositeLaminate {
    /// Ply layers, from bottom to top.
    pub plies: Vec<PlyLayer>,
    /// Fibre angles (radians) for each ply.
    pub angles: Vec<f64>,
}
impl CompositeLaminate {
    /// Create a laminate from a stack of plies and their fibre angles.
    pub fn new(plies: Vec<PlyLayer>, angles: Vec<f64>) -> Self {
        assert_eq!(plies.len(), angles.len());
        Self { plies, angles }
    }
    /// Total laminate thickness.
    pub fn total_thickness(&self) -> f64 {
        self.plies.iter().map(|p| p.thickness).sum()
    }
    /// Compute the full 6×6 ABD stiffness matrix.
    ///
    /// Rows/cols 0–2: extensional (A), 3–5: coupling (B) / bending (D).
    pub fn abd_matrix(&self) -> [[f64; 6]; 6] {
        let n = self.plies.len();
        let total_h = self.total_thickness();
        let mut z = vec![0.0_f64; n + 1];
        z[0] = -total_h / 2.0;
        for k in 0..n {
            z[k + 1] = z[k] + self.plies[k].thickness;
        }
        let mut a = [[0.0_f64; 3]; 3];
        let mut b = [[0.0_f64; 3]; 3];
        let mut d = [[0.0_f64; 3]; 3];
        for k in 0..n {
            let q = self.plies[k].q_matrix_global(self.angles[k]);
            let zk = z[k];
            let zk1 = z[k + 1];
            let dz1 = zk1 - zk;
            let dz2 = (zk1 * zk1 - zk * zk) / 2.0;
            let dz3 = (zk1 * zk1 * zk1 - zk * zk * zk) / 3.0;
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] += q[i][j] * dz1;
                    b[i][j] += q[i][j] * dz2;
                    d[i][j] += q[i][j] * dz3;
                }
            }
        }
        let mut abd = [[0.0_f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                abd[i][j] = a[i][j];
                abd[i][j + 3] = b[i][j];
                abd[i + 3][j] = b[i][j];
                abd[i + 3][j + 3] = d[i][j];
            }
        }
        abd
    }
}
/// Thin shell element with in-plane (membrane) and buckling capabilities.
pub struct ShellElement {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Thickness \[m\].
    pub thickness: f64,
}
impl ShellElement {
    /// Create a new shell element.
    pub fn new(young_modulus: f64, poisson_ratio: f64, thickness: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            thickness,
        }
    }
    /// Compute in-plane membrane stress resultants N = \[Nx, Ny, Nxy\] \[N/m\]
    /// for a CST triangle using plane-stress constitutive relation.
    ///
    /// # Arguments
    /// * `nodes` – 3 × \[x, y\] nodal coordinates
    /// * `disp`  – 6-element displacement vector \[u0,v0, u1,v1, u2,v2\]
    ///
    /// Returns `[Nx, Ny, Nxy]`.
    pub fn compute_membrane_forces(&self, nodes: &[[f64; 2]; 3], disp: &[f64; 6]) -> [f64; 3] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let t = self.thickness;
        let (x0, y0) = (nodes[0][0], nodes[0][1]);
        let (x1, y1) = (nodes[1][0], nodes[1][1]);
        let (x2, y2) = (nodes[2][0], nodes[2][1]);
        let two_a = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        if two_a.abs() < 1e-30 {
            return [0.0; 3];
        }
        let b = [
            [
                (y1 - y2) / two_a,
                0.0,
                (y2 - y0) / two_a,
                0.0,
                (y0 - y1) / two_a,
                0.0,
            ],
            [
                0.0,
                (x2 - x1) / two_a,
                0.0,
                (x0 - x2) / two_a,
                0.0,
                (x1 - x0) / two_a,
            ],
            [
                (x2 - x1) / two_a,
                (y1 - y2) / two_a,
                (x0 - x2) / two_a,
                (y2 - y0) / two_a,
                (x1 - x0) / two_a,
                (y0 - y1) / two_a,
            ],
        ];
        let mut eps = [0.0_f64; 3];
        for i in 0..3 {
            for j in 0..6 {
                eps[i] += b[i][j] * disp[j];
            }
        }
        let c_fac = e / (1.0 - nu * nu);
        let c = [
            [c_fac, c_fac * nu, 0.0],
            [c_fac * nu, c_fac, 0.0],
            [0.0, 0.0, c_fac * (1.0 - nu) / 2.0],
        ];
        let mut sigma = [0.0_f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                sigma[i] += c[i][j] * eps[j];
            }
        }
        [sigma[0] * t, sigma[1] * t, sigma[2] * t]
    }
    /// Compute the linear buckling critical load parameter N_cr for the element \[N/m\].
    ///
    /// Uses the classical plate buckling formula for a simply supported
    /// rectangular equivalent:  N_cr = k * π² D / b²
    ///
    /// For a triangle with legs `a` and `b`, an equivalent rectangle is formed
    /// with sides equal to the two leg lengths.  k = 4.0 (four-sided buckling).
    ///
    /// # Arguments
    /// * `nodes` – 3 × \[x, y\] nodal coordinates
    pub fn compute_buckling_load(&self, nodes: &[[f64; 2]; 3]) -> f64 {
        let t = self.thickness;
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let d = e * t * t * t / (12.0 * (1.0 - nu * nu));
        let sides = [
            ((nodes[1][0] - nodes[0][0]).powi(2) + (nodes[1][1] - nodes[0][1]).powi(2)).sqrt(),
            ((nodes[2][0] - nodes[1][0]).powi(2) + (nodes[2][1] - nodes[1][1]).powi(2)).sqrt(),
            ((nodes[0][0] - nodes[2][0]).powi(2) + (nodes[0][1] - nodes[2][1]).powi(2)).sqrt(),
        ];
        let b = sides
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
            .max(1e-30);
        4.0 * std::f64::consts::PI * std::f64::consts::PI * d / (b * b)
    }
    /// Compute the von Karman large-deflection correction vector (9-element, one per DOF).
    ///
    /// The von Karman correction to the bending-stiffness residual is:
    ///   F_nl = ∫ B_b^T * (N_m · ∇w ∇w^T) dA
    ///
    /// For a CST triangle this is evaluated by the product of average membrane
    /// resultants and the gradient of the transverse displacement field:
    ///   ∇w = \[∂w/∂x, ∂w/∂y\]
    ///
    /// # Arguments
    /// * `nodes` – 3 × \[x, y\] nodal positions
    /// * `disp`  – 9-element displacement: \[w0,θx0,θy0, w1,θx1,θy1, w2,θx2,θy2\]
    ///
    /// Returns a 9-element correction force vector.
    pub fn compute_large_deflection_correction(
        &self,
        nodes: &[[f64; 2]; 3],
        disp: &[f64; 9],
    ) -> [f64; 9] {
        let (x0, y0) = (nodes[0][0], nodes[0][1]);
        let (x1, y1) = (nodes[1][0], nodes[1][1]);
        let (x2, y2) = (nodes[2][0], nodes[2][1]);
        let two_a = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        if two_a.abs() < 1e-30 {
            return [0.0; 9];
        }
        let area = two_a.abs() / 2.0;
        let dndx = [(y1 - y2) / two_a, (y2 - y0) / two_a, (y0 - y1) / two_a];
        let dndy = [(x2 - x1) / two_a, (x0 - x2) / two_a, (x1 - x0) / two_a];
        let w = [disp[0], disp[3], disp[6]];
        let dw_dx: f64 = (0..3).map(|i| dndx[i] * w[i]).sum();
        let dw_dy: f64 = (0..3).map(|i| dndy[i] * w[i]).sum();
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let t = self.thickness;
        let c_fac = e * t / (1.0 - nu * nu);
        let eps_xx_nl = 0.5 * dw_dx * dw_dx;
        let eps_yy_nl = 0.5 * dw_dy * dw_dy;
        let eps_xy_nl = dw_dx * dw_dy;
        let nx = c_fac * (eps_xx_nl + nu * eps_yy_nl);
        let ny = c_fac * (nu * eps_xx_nl + eps_yy_nl);
        let nxy = c_fac * (1.0 - nu) / 2.0 * eps_xy_nl;
        let mut f_corr = [0.0_f64; 9];
        for i in 0..3 {
            let fi = area
                * (nx * dndx[i] * dw_dx
                    + ny * dndy[i] * dw_dy
                    + nxy * (dndx[i] * dw_dy + dndy[i] * dw_dx));
            f_corr[i * 3] = fi;
        }
        f_corr
    }
}
/// 4-node rectangular plate element (Kirchhoff, simplified diagonal stiffness).
///
/// DOF per node: (w, θx, θy) → 12 DOF total.
pub struct RectPlateElement {
    /// Plate material parameters.
    pub params: PlateParams,
    /// Half-length in x direction.
    pub a: f64,
    /// Half-width in y direction.
    pub b: f64,
    /// Node positions (centered at origin): (±a, ±b).
    pub node_positions: [[f64; 2]; 4],
}
impl RectPlateElement {
    /// Create element centered at origin with half-dimensions `a` × `b`.
    pub fn new(params: PlateParams, a: f64, b: f64) -> Self {
        let node_positions = [[-a, -b], [a, -b], [a, b], [-a, b]];
        Self {
            params,
            a,
            b,
            node_positions,
        }
    }
    /// Approximate 12×12 stiffness matrix using classical plate theory (diagonal blocks).
    ///
    /// Uses D·(1/a²+1/b²) scaling for bending DOFs – a closed-form diagonal approximation.
    pub fn stiffness_matrix_12x12(&self) -> Vec<Vec<f64>> {
        let d = self.params.flexural_rigidity();
        let a = self.a;
        let b = self.b;
        let a2 = a * a;
        let b2 = b * b;
        let kw = d * (2.0 * b / (a2 * a) + 2.0 * a / (b2 * b));
        let ktx = d * (2.0 * b / a);
        let kty = d * (2.0 * a / b);
        let mut k = vec![vec![0.0f64; 12]; 12];
        for node in 0..4 {
            let base = node * 3;
            k[base][base] = kw;
            k[base + 1][base + 1] = ktx;
            k[base + 2][base + 2] = kty;
        }
        k
    }
    /// Consistent 12×12 mass matrix: M = ρ·h·∫N_I·N_J dA (diagonal lumping approximation).
    ///
    /// Each translational DOF gets ρ·h·A/4; rotational DOFs get ρ·h·A/4·(a²+b²)/12.
    pub fn consistent_mass_matrix(&self, rho: f64) -> Vec<Vec<f64>> {
        let h = self.params.thickness;
        let a = self.a;
        let b = self.b;
        let area = (2.0 * a) * (2.0 * b);
        let mw = rho * h * area / 4.0;
        let mr = rho * h * area / 4.0 * (a * a + b * b) / 12.0;
        let mut m = vec![vec![0.0f64; 12]; 12];
        for node in 0..4 {
            let base = node * 3;
            m[base][base] = mw;
            m[base + 1][base + 1] = mr;
            m[base + 2][base + 2] = mr;
        }
        m
    }
}
/// Flat 3-node triangular shell element combining:
/// - Membrane: CST (constant strain triangle, 2 DOF/node in-plane)
/// - Bending:  DKT-style (3 DOF/node: w, θ_x, θ_y)
///
/// Each node has 6 DOF: (u, v, w, θ_x, θ_y, θ_z).
/// Returns 18×18 stiffness matrix.
pub struct FlatShellElement {
    /// Young's modulus.
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Density (kg/m³).
    pub density: f64,
    /// Shell thickness.
    pub thickness: f64,
}
impl FlatShellElement {
    /// Create a new flat shell element.
    pub fn new(young_modulus: f64, poisson_ratio: f64, density: f64, thickness: f64) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            density,
            thickness,
        }
    }
    /// Compute 18×18 stiffness matrix by assembling membrane (6×6) and bending (9×9) parts.
    ///
    /// `nodes`: 3 planar (x,y) node positions.
    ///
    /// DOF ordering per node i: (u_i, v_i, w_i, θx_i, θy_i, θz_i).
    /// θ_z is the drilling DOF (assigned a small penalty stiffness).
    pub fn stiffness_triangle(&self, nodes: &[[f64; 2]; 3]) -> Vec<Vec<f64>> {
        let k_mem = self.membrane_stiffness(nodes);
        let k_bend = self.bending_stiffness(nodes);
        let mut k = vec![vec![0.0_f64; 18]; 18];
        let mem_map = [0, 1, 6, 7, 12, 13];
        for i in 0..6 {
            for j in 0..6 {
                k[mem_map[i]][mem_map[j]] += k_mem[i][j];
            }
        }
        let bend_map = [2, 3, 4, 8, 9, 10, 14, 15, 16];
        for i in 0..9 {
            for j in 0..9 {
                k[bend_map[i]][bend_map[j]] += k_bend[i][j];
            }
        }
        let area = {
            let (x0, y0) = (nodes[0][0], nodes[0][1]);
            let (x1, y1) = (nodes[1][0], nodes[1][1]);
            let (x2, y2) = (nodes[2][0], nodes[2][1]);
            0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)).abs()
        };
        let drill_stiff = self.young_modulus * self.thickness * area * 1e-6;
        k[5][5] += drill_stiff;
        k[11][11] += drill_stiff;
        k[17][17] += drill_stiff;
        k
    }
    fn membrane_stiffness(&self, nodes: &[[f64; 2]; 3]) -> [[f64; 6]; 6] {
        let elem = MembraneTriangle::new(self.thickness, self.young_modulus, self.poisson_ratio);
        elem.stiffness_2d(nodes)
    }
    fn bending_stiffness(&self, nodes: &[[f64; 2]; 3]) -> [[f64; 9]; 9] {
        let plate = KirchhoffPlate::new(
            self.thickness,
            self.young_modulus,
            self.poisson_ratio,
            self.density,
        );
        let d_mat = plate.constitutive_matrix();
        let (x0, y0) = (nodes[0][0], nodes[0][1]);
        let (x1, y1) = (nodes[1][0], nodes[1][1]);
        let (x2, y2) = (nodes[2][0], nodes[2][1]);
        let area = 0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)).abs();
        if area < 1e-15 {
            return [[0.0; 9]; 9];
        }
        let two_a = 2.0 * area;
        let b1x = (y1 - y2) / two_a;
        let b2x = (y2 - y0) / two_a;
        let b3x = (y0 - y1) / two_a;
        let b1y = (x2 - x1) / two_a;
        let b2y = (x0 - x2) / two_a;
        let b3y = (x1 - x0) / two_a;
        let b_bend = [
            [0.0, 0.0, b1x, 0.0, 0.0, b2x, 0.0, 0.0, b3x],
            [0.0, -b1y, 0.0, 0.0, -b2y, 0.0, 0.0, -b3y, 0.0],
            [0.0, -b1x, b1y, 0.0, -b2x, b2y, 0.0, -b3x, b3y],
        ];
        let mut db = [[0.0_f64; 9]; 3];
        for i in 0..3 {
            for j in 0..9 {
                for l in 0..3 {
                    db[i][j] += d_mat[i][l] * b_bend[l][j];
                }
            }
        }
        let mut k_bend = [[0.0_f64; 9]; 9];
        for i in 0..9 {
            for j in 0..9 {
                let mut val = 0.0;
                for l in 0..3 {
                    val += b_bend[l][i] * db[l][j];
                }
                k_bend[i][j] = area * val;
            }
        }
        k_bend
    }
}
/// Membrane element (2D plane-stress triangle, 3 nodes, 2 DOF/node = 6×6 stiffness).
///
/// This is the CST (Constant Strain Triangle) element with uniform strain field.
pub struct MembraneTriangle {
    /// Element thickness.
    pub thickness: f64,
    /// Young's modulus.
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
}
impl MembraneTriangle {
    /// Create a new membrane triangle element.
    pub fn new(thickness: f64, e: f64, nu: f64) -> Self {
        Self {
            thickness,
            young_modulus: e,
            poisson_ratio: nu,
        }
    }
    /// 3×3 plane stress constitutive matrix.
    ///
    /// D = E/(1-ν²) * \[\[1, ν, 0\\], \[ν, 1, 0\], \[0, 0, (1-ν)/2\]]
    pub fn constitutive_matrix(&self) -> [[f64; 3]; 3] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let c = e / (1.0 - nu * nu);
        [
            [c, c * nu, 0.0],
            [c * nu, c, 0.0],
            [0.0, 0.0, c * (1.0 - nu) / 2.0],
        ]
    }
    /// Compute 6×6 stiffness matrix for CST (constant strain triangle).
    ///
    /// `nodes`: 3 planar positions (x,y) for the triangle vertices.
    ///
    /// The B-matrix for CST is 3×6, constant within the element:
    /// K = t * A * B^T * D * B
    pub fn stiffness_2d(&self, nodes: &[[f64; 2]; 3]) -> [[f64; 6]; 6] {
        let area = Self::area(nodes);
        if area.abs() < 1e-15 {
            return [[0.0; 6]; 6];
        }
        let x0 = nodes[0][0];
        let y0 = nodes[0][1];
        let x1 = nodes[1][0];
        let y1 = nodes[1][1];
        let x2 = nodes[2][0];
        let y2 = nodes[2][1];
        let two_a = 2.0 * area;
        let dn0dx = (y1 - y2) / two_a;
        let dn0dy = (x2 - x1) / two_a;
        let dn1dx = (y2 - y0) / two_a;
        let dn1dy = (x0 - x2) / two_a;
        let dn2dx = (y0 - y1) / two_a;
        let dn2dy = (x1 - x0) / two_a;
        let bmat = [
            [dn0dx, 0.0, dn1dx, 0.0, dn2dx, 0.0],
            [0.0, dn0dy, 0.0, dn1dy, 0.0, dn2dy],
            [dn0dy, dn0dx, dn1dy, dn1dx, dn2dy, dn2dx],
        ];
        let d = self.constitutive_matrix();
        let mut db = [[0.0f64; 6]; 3];
        for i in 0..3 {
            for j in 0..6 {
                for l in 0..3 {
                    db[i][j] += d[i][l] * bmat[l][j];
                }
            }
        }
        let mut k = [[0.0f64; 6]; 6];
        let ta = self.thickness * area;
        for i in 0..6 {
            for j in 0..6 {
                let mut sum = 0.0;
                for l in 0..3 {
                    sum += bmat[l][i] * db[l][j];
                }
                k[i][j] = ta * sum;
            }
        }
        k
    }
    /// Triangle area from 2D coordinates.
    ///
    /// Returns 0.5 * |det(\[x1-x0, y1-y0; x2-x0, y2-y0\])|
    pub fn area(nodes: &[[f64; 2]; 3]) -> f64 {
        let (x0, y0) = (nodes[0][0], nodes[0][1]);
        let (x1, y1) = (nodes[1][0], nodes[1][1]);
        let (x2, y2) = (nodes[2][0], nodes[2][1]);
        0.5 * ((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)).abs()
    }
}
/// Reissner-Mindlin thick plate/shell (accounts for transverse shear).
pub struct MindlinShell {
    /// Plate/shell thickness.
    pub thickness: f64,
    /// Young's modulus.
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Shear correction factor (κ = 5/6 for homogeneous plates).
    pub shear_correction: f64,
}
impl MindlinShell {
    /// Create a new Mindlin shell element with default shear correction κ = 5/6.
    pub fn new(thickness: f64, e: f64, nu: f64) -> Self {
        Self {
            thickness,
            young_modulus: e,
            poisson_ratio: nu,
            shear_correction: 5.0 / 6.0,
        }
    }
    /// Shear modulus G = E / (2*(1+ν)).
    pub fn shear_modulus(&self) -> f64 {
        self.young_modulus / (2.0 * (1.0 + self.poisson_ratio))
    }
    /// Transverse shear stiffness coefficient: Ks = κ * G * t.
    pub fn shear_stiffness(&self) -> f64 {
        self.shear_correction * self.shear_modulus() * self.thickness
    }
    /// Bending stiffness matrix (3×3) — same form as Kirchhoff.
    ///
    /// Relates bending moments {M} to curvatures {κ}.
    pub fn bending_matrix(&self) -> [[f64; 3]; 3] {
        let t = self.thickness;
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let d = e * t * t * t / (12.0 * (1.0 - nu * nu));
        [
            [d, d * nu, 0.0],
            [d * nu, d, 0.0],
            [0.0, 0.0, d * (1.0 - nu) / 2.0],
        ]
    }
    /// Shear stiffness matrix (2×2) for \[Q_xz, Q_yz\].
    ///
    /// Returns Ks * I₂ where Ks = κ*G*t.
    pub fn shear_matrix(&self) -> [[f64; 2]; 2] {
        let ks = self.shear_stiffness();
        [[ks, 0.0], [0.0, ks]]
    }
    /// Thickness-to-span ratio check: thin plate if t/L < 0.05.
    pub fn is_thin(&self, span: f64) -> bool {
        self.thickness / span < 0.05
    }
}
/// Kirchhoff-Love thin plate element (4-node bilinear, reduced integration).
///
/// Assumes small deformations and thin plate theory (transverse shear ignored).
pub struct KirchhoffPlate {
    /// Plate thickness.
    pub thickness: f64,
    /// Young's modulus.
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Mass density.
    pub density: f64,
}
impl KirchhoffPlate {
    /// Create a new Kirchhoff plate element.
    pub fn new(thickness: f64, e: f64, nu: f64, density: f64) -> Self {
        Self {
            thickness,
            young_modulus: e,
            poisson_ratio: nu,
            density,
        }
    }
    /// Flexural rigidity: D = E*t³ / (12*(1-ν²))
    pub fn flexural_rigidity(&self) -> f64 {
        let t = self.thickness;
        self.young_modulus * t * t * t / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio))
    }
    /// 3×3 plate constitutive matrix \[D_plate\].
    ///
    /// Relates bending moments to curvatures: {M} = \[D_plate\] * {κ}
    ///
    /// D_plate = D * \[\[1, ν, 0\\], \[ν, 1, 0\], \[0, 0, (1-ν)/2\]]
    pub fn constitutive_matrix(&self) -> [[f64; 3]; 3] {
        let d = self.flexural_rigidity();
        let nu = self.poisson_ratio;
        [
            [d, d * nu, 0.0],
            [d * nu, d, 0.0],
            [0.0, 0.0, d * (1.0 - nu) / 2.0],
        ]
    }
    /// Stiffness matrix for a rectangular plate element a × b.
    ///
    /// Returns 12×12 matrix (3 DOF per node: w, θx, θy).
    /// Uses Gauss quadrature over the element of the curvature energy:
    /// K\[I,J\] = ∫ κ_I^T D_plate κ_J dA, where κ contains second derivatives
    /// of the Hermite shape functions.
    ///
    /// DOF ordering: node 0=(0,0), node 1=(a,0), node 2=(a,b), node 3=(0,b).
    /// Per node: \[w, θx (=∂w/∂y), θy (=∂w/∂x)\].
    pub fn rectangular_stiffness(&self, a: f64, b: f64) -> [[f64; 12]; 12] {
        let d = self.flexural_rigidity();
        let nu = self.poisson_ratio;
        let a2 = a * a;
        let b2 = b * b;
        let gp = [-f64::sqrt(3.0 / 5.0), 0.0, f64::sqrt(3.0 / 5.0)];
        let gw = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        let hermite_val = |xi: f64, which: usize| -> f64 {
            match which {
                0 => 0.25 * (2.0 - 3.0 * xi + xi * xi * xi),
                1 => 0.25 * (1.0 - xi - xi * xi + xi * xi * xi),
                2 => 0.25 * (2.0 + 3.0 * xi - xi * xi * xi),
                3 => 0.25 * (-1.0 - xi + xi * xi + xi * xi * xi),
                _ => 0.0,
            }
        };
        let hermite_d1 = |xi: f64, which: usize| -> f64 {
            match which {
                0 => 0.25 * (-3.0 + 3.0 * xi * xi),
                1 => 0.25 * (-1.0 - 2.0 * xi + 3.0 * xi * xi),
                2 => 0.25 * (3.0 - 3.0 * xi * xi),
                3 => 0.25 * (-1.0 + 2.0 * xi + 3.0 * xi * xi),
                _ => 0.0,
            }
        };
        let hermite_d2 = |xi: f64, which: usize| -> f64 {
            match which {
                0 => 3.0 * xi / 2.0,
                1 => (3.0 * xi - 1.0) / 2.0,
                2 => -3.0 * xi / 2.0,
                3 => (3.0 * xi + 1.0) / 2.0,
                _ => 0.0,
            }
        };
        let dof_desc: [(usize, f64, usize, f64); 12] = [
            (0, 1.0, 0, 1.0),
            (0, 1.0, 1, b),
            (1, a, 0, 1.0),
            (2, 1.0, 0, 1.0),
            (2, 1.0, 1, b),
            (3, a, 0, 1.0),
            (2, 1.0, 2, 1.0),
            (2, 1.0, 3, b),
            (3, a, 2, 1.0),
            (0, 1.0, 2, 1.0),
            (0, 1.0, 3, b),
            (1, a, 2, 1.0),
        ];
        let mut k = [[0.0f64; 12]; 12];
        for i in 0..12 {
            let (xi_b_i, sx_i, eta_b_i, sy_i) = dof_desc[i];
            for j in 0..12 {
                let (xi_b_j, sx_j, eta_b_j, sy_j) = dof_desc[j];
                let mut val = 0.0;
                for (ig, &xi) in gp.iter().enumerate() {
                    for (jg, &eta) in gp.iter().enumerate() {
                        let w = gw[ig] * gw[jg];
                        let d2ni_dx2 = sx_i
                            * (4.0 / a2)
                            * hermite_d2(xi, xi_b_i)
                            * sy_i
                            * hermite_val(eta, eta_b_i);
                        let d2ni_dy2 = sx_i
                            * hermite_val(xi, xi_b_i)
                            * sy_i
                            * (4.0 / b2)
                            * hermite_d2(eta, eta_b_i);
                        let d2ni_dxy = sx_i
                            * (2.0 / a)
                            * hermite_d1(xi, xi_b_i)
                            * sy_i
                            * (2.0 / b)
                            * hermite_d1(eta, eta_b_i);
                        let d2nj_dx2 = sx_j
                            * (4.0 / a2)
                            * hermite_d2(xi, xi_b_j)
                            * sy_j
                            * hermite_val(eta, eta_b_j);
                        let d2nj_dy2 = sx_j
                            * hermite_val(xi, xi_b_j)
                            * sy_j
                            * (4.0 / b2)
                            * hermite_d2(eta, eta_b_j);
                        let d2nj_dxy = sx_j
                            * (2.0 / a)
                            * hermite_d1(xi, xi_b_j)
                            * sy_j
                            * (2.0 / b)
                            * hermite_d1(eta, eta_b_j);
                        let integrand = d
                            * (d2ni_dx2 * d2nj_dx2
                                + d2ni_dy2 * d2nj_dy2
                                + nu * (d2ni_dx2 * d2nj_dy2 + d2ni_dy2 * d2nj_dx2)
                                + 2.0 * (1.0 - nu) * d2ni_dxy * d2nj_dxy);
                        val += w * integrand * (a * b / 4.0);
                    }
                }
                k[i][j] = val;
            }
        }
        k
    }
    /// Consistent mass matrix for rectangular element a × b.
    ///
    /// Uses Gauss quadrature of ρ*t*N_I*N_J over the element domain.
    /// Returns a 12×12 consistent mass matrix.
    pub fn rectangular_mass(&self, a: f64, b: f64) -> [[f64; 12]; 12] {
        let rho_t = self.density * self.thickness;
        let gp = [-f64::sqrt(3.0 / 5.0), 0.0, f64::sqrt(3.0 / 5.0)];
        let gw = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0];
        let hermite_val = |xi: f64, which: usize| -> f64 {
            match which {
                0 => 0.25 * (2.0 - 3.0 * xi + xi * xi * xi),
                1 => 0.25 * (1.0 - xi - xi * xi + xi * xi * xi),
                2 => 0.25 * (2.0 + 3.0 * xi - xi * xi * xi),
                3 => 0.25 * (-1.0 - xi + xi * xi + xi * xi * xi),
                _ => 0.0,
            }
        };
        let dof_desc: [(usize, f64, usize, f64); 12] = [
            (0, 1.0, 0, 1.0),
            (0, 1.0, 1, b),
            (1, a, 0, 1.0),
            (2, 1.0, 0, 1.0),
            (2, 1.0, 1, b),
            (3, a, 0, 1.0),
            (2, 1.0, 2, 1.0),
            (2, 1.0, 3, b),
            (3, a, 2, 1.0),
            (0, 1.0, 2, 1.0),
            (0, 1.0, 3, b),
            (1, a, 2, 1.0),
        ];
        let mut m = [[0.0f64; 12]; 12];
        for i in 0..12 {
            let (xi_b_i, sx_i, eta_b_i, sy_i) = dof_desc[i];
            for j in 0..12 {
                let (xi_b_j, sx_j, eta_b_j, sy_j) = dof_desc[j];
                let mut val = 0.0;
                for (ig, &xi) in gp.iter().enumerate() {
                    for (jg, &eta) in gp.iter().enumerate() {
                        let w = gw[ig] * gw[jg];
                        let ni = sx_i * hermite_val(xi, xi_b_i) * sy_i * hermite_val(eta, eta_b_i);
                        let nj = sx_j * hermite_val(xi, xi_b_j) * sy_j * hermite_val(eta, eta_b_j);
                        val += w * rho_t * ni * nj * (a * b / 4.0);
                    }
                }
                m[i][j] = val;
            }
        }
        m
    }
    /// Natural frequencies of a simply-supported rectangular plate.
    ///
    /// ω_mn = π² * sqrt(D/(ρ*t)) * ((m/a)² + (n/b)²)
    pub fn natural_frequency_ss(&self, m: u32, n: u32, a: f64, b: f64) -> f64 {
        let d = self.flexural_rigidity();
        let mf = m as f64;
        let nf = n as f64;
        std::f64::consts::PI
            * std::f64::consts::PI
            * f64::sqrt(d / (self.density * self.thickness))
            * ((mf / a).powi(2) + (nf / b).powi(2))
    }
}
/// Mindlin-Reissner 4-node plate element (size: full dimensions).
pub struct MindlinElement {
    /// Mindlin plate parameters.
    pub params: MindlinParams,
    /// Element dimensions \[width_x, width_y\].
    pub size: [f64; 2],
}
impl MindlinElement {
    /// Create a Mindlin element with given parameters and size.
    pub fn new(params: MindlinParams, size: [f64; 2]) -> Self {
        Self { params, size }
    }
}
/// Sandwich plate model with stiff face sheets and a flexible core.
///
/// Assumes equal face-sheet thickness on both sides and a homogeneous core.
#[derive(Debug, Clone)]
pub struct SandwichPlate {
    /// Face-sheet Young's modulus \[Pa\].
    pub face_e: f64,
    /// Face-sheet Poisson's ratio.
    pub face_nu: f64,
    /// Face-sheet thickness (each face) \[m\].
    pub face_thickness: f64,
    /// Core shear modulus \[Pa\].
    pub core_g: f64,
    /// Core Poisson's ratio (rarely used, included for completeness).
    pub core_nu: f64,
    /// Core thickness \[m\].
    pub core_thickness: f64,
}
impl SandwichPlate {
    /// Create a new sandwich plate.
    pub fn new(
        face_e: f64,
        face_nu: f64,
        face_thickness: f64,
        core_g: f64,
        core_nu: f64,
        core_thickness: f64,
    ) -> Self {
        Self {
            face_e,
            face_nu,
            face_thickness,
            core_g,
            core_nu,
            core_thickness,
        }
    }
    /// Total plate thickness: 2 * face + core.
    pub fn total_thickness(&self) -> f64 {
        2.0 * self.face_thickness + self.core_thickness
    }
    /// Effective flexural rigidity D based on face sheets only (parallel axis theorem).
    ///
    /// D ≈ E_f * t_f * d² / (1 - ν_f²) where d = (core_thickness + face_thickness) / 2.
    pub fn effective_flexural_rigidity(&self) -> f64 {
        let d = (self.core_thickness + self.face_thickness) / 2.0;
        self.face_e * self.face_thickness * d * d / (1.0 - self.face_nu * self.face_nu)
    }
    /// Transverse shear stiffness of the core: Ks = G_core * core_thickness.
    pub fn transverse_shear_stiffness(&self) -> f64 {
        self.core_g * self.core_thickness
    }
    /// Face-sheet membrane stiffness: A = E_f * t_f / (1 - ν_f²).
    pub fn face_membrane_stiffness(&self) -> f64 {
        self.face_e * self.face_thickness / (1.0 - self.face_nu * self.face_nu)
    }
    /// Critical wrinkling stress (face-sheet wrinkling under compression).
    ///
    /// σ_cr = 0.5 (E_f E_c G_c)^{1/3}   (approximate formula)
    pub fn wrinkling_stress(&self, core_e: f64) -> f64 {
        0.5 * (self.face_e * core_e * self.core_g).powf(1.0 / 3.0)
    }
}
/// Material and geometry parameters for a Mindlin-Reissner (shear-deformable) plate.
pub struct MindlinParams {
    /// Young's modulus.
    pub e: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Plate thickness.
    pub thickness: f64,
    /// Shear correction factor (κ = 5/6 for homogeneous plates).
    pub shear_correction_factor: f64,
}
impl MindlinParams {
    /// Create Mindlin plate parameters with given shear correction factor.
    pub fn new(e: f64, nu: f64, thickness: f64, shear_correction_factor: f64) -> Self {
        Self {
            e,
            nu,
            thickness,
            shear_correction_factor,
        }
    }
    /// Transverse shear stiffness: Ks = κ·G·h where G = E / (2·(1+ν)).
    pub fn shear_stiffness(&self) -> f64 {
        let g = self.e / (2.0 * (1.0 + self.nu));
        self.shear_correction_factor * g * self.thickness
    }
}
/// Thin circular cylindrical shell element (Donnell-Mushtari-Vlasov theory).
#[derive(Debug, Clone)]
pub struct CylindricalShellElement {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Mass density \[kg/m³\].
    pub density: f64,
    /// Shell thickness \[m\].
    pub thickness: f64,
    /// Shell length \[m\].
    pub length: f64,
    /// Shell radius \[m\].
    pub radius: f64,
}
impl CylindricalShellElement {
    /// Create a new cylindrical shell element.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        thickness: f64,
        length: f64,
        radius: f64,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            density,
            thickness,
            length,
            radius,
        }
    }
    /// Ring frequency (n=1 breathing mode): ω_ring = c_L / R where c_L = √(E/(ρ(1-ν²))).
    pub fn ring_frequency(&self) -> f64 {
        let c_l = (self.young_modulus / (self.density * (1.0 - self.poisson_ratio.powi(2)))).sqrt();
        c_l / self.radius
    }
    /// Classical critical external pressure (Timoshenko formula, long shell).
    ///
    /// p_cr = 2E/(1−ν²) * (t/2R)³
    pub fn critical_pressure_external(&self) -> f64 {
        let t = self.thickness;
        let r = self.radius;
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        2.0 * e / (1.0 - nu * nu) * (t / (2.0 * r)).powi(3)
    }
    /// Hoop (circumferential) stress under internal pressure p: σ_θ = p*R/t.
    pub fn hoop_stress_internal_pressure(&self, p: f64) -> f64 {
        p * self.radius / self.thickness
    }
    /// Axial stress for closed-end vessel under internal pressure: σ_z = p*R/(2t).
    pub fn axial_stress_closed_end(&self, p: f64) -> f64 {
        p * self.radius / (2.0 * self.thickness)
    }
    /// Classical axial buckling load (Euler formula for cylindrical shell):
    ///
    /// N_cr = 2π²EI/L² = E*h / (R * sqrt(3*(1-ν²)))
    pub fn axial_buckling_load(&self) -> f64 {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let t = self.thickness;
        let r = self.radius;
        e * t / (r * (3.0 * (1.0 - nu * nu)).sqrt())
    }
    /// Longitudinal-to-radius ratio (slenderness parameter).
    pub fn slenderness(&self) -> f64 {
        self.length / self.radius
    }
}
/// Orthotropic thin plate using classical Kirchhoff theory.
///
/// Material has different stiffnesses in x and y directions.
pub struct OrthotropicPlate {
    /// Young's modulus in x direction \[Pa\].
    pub e1: f64,
    /// Young's modulus in y direction \[Pa\].
    pub e2: f64,
    /// In-plane shear modulus \[Pa\].
    pub g12: f64,
    /// Poisson's ratio ν₁₂.
    pub nu12: f64,
    /// Plate thickness \[m\].
    pub thickness: f64,
}
impl OrthotropicPlate {
    /// Create a new orthotropic plate.
    pub fn new(e1: f64, e2: f64, g12: f64, nu12: f64, thickness: f64) -> Self {
        Self {
            e1,
            e2,
            g12,
            nu12,
            thickness,
        }
    }
    /// Compute the four orthotropic flexural rigidities (Dx, Dy, Dxy, Dν).
    ///
    /// Returns (Dx, Dy, Dxy, Dν) where:
    /// - Dx = E1 h³ / (12(1-ν12 ν21))
    /// - Dy = E2 h³ / (12(1-ν12 ν21))
    /// - Dxy = G12 h³ / 12  (twisting)
    /// - Dν = ν12 Dy
    pub fn rigidities(&self) -> (f64, f64, f64, f64) {
        let h = self.thickness;
        let nu21 = self.nu12 * self.e2 / self.e1;
        let denom = 1.0 - self.nu12 * nu21;
        let h3_12 = h * h * h / 12.0;
        let dx = self.e1 * h3_12 / denom;
        let dy = self.e2 * h3_12 / denom;
        let dxy = self.g12 * h3_12;
        let dnu = self.nu12 * dy;
        (dx, dy, dxy, dnu)
    }
    /// Natural frequency of a simply-supported orthotropic plate, mode (m, n).
    ///
    /// ω_mn = π² \[(m/a)⁴ Dx + 2(m/a)²(n/b)² (Dν+2Dxy) + (n/b)⁴ Dy\]^{1/2} / √(ρ·h)
    pub fn natural_frequency_ss(&self, m: u32, n: u32, a: f64, b: f64, rho: f64) -> f64 {
        let (dx, dy, dxy, dnu) = self.rigidities();
        let mf = m as f64;
        let nf = n as f64;
        let pi2 = std::f64::consts::PI * std::f64::consts::PI;
        let alpha = mf / a;
        let beta = nf / b;
        let num = dx * alpha.powi(4)
            + 2.0 * (dnu + 2.0 * dxy) * alpha * alpha * beta * beta
            + dy * beta.powi(4);
        let rho_h = rho * self.thickness;
        pi2 * (num / rho_h).sqrt()
    }
}
/// Material and geometry parameters for a thin plate (Kirchhoff-Love theory).
pub struct PlateParams {
    /// Young's modulus.
    pub e: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Plate thickness.
    pub thickness: f64,
}
impl PlateParams {
    /// Create new plate parameters.
    pub fn new(e: f64, nu: f64, thickness: f64) -> Self {
        Self { e, nu, thickness }
    }
    /// Flexural rigidity: D = E·h³ / (12·(1−ν²)).
    pub fn flexural_rigidity(&self) -> f64 {
        let h = self.thickness;
        self.e * h * h * h / (12.0 * (1.0 - self.nu * self.nu))
    }
    /// 3×3 bending stiffness matrix (Voigt notation): Db = D · \[\[1,ν,0\\],\[ν,1,0\],\[0,0,(1-ν)/2\]].
    pub fn bending_stiffness_matrix(&self) -> [[f64; 3]; 3] {
        let d = self.flexural_rigidity();
        let nu = self.nu;
        [
            [d, d * nu, 0.0],
            [d * nu, d, 0.0],
            [0.0, 0.0, d * (1.0 - nu) / 2.0],
        ]
    }
}
/// A single orthotropic ply layer.
///
/// Properties are defined in the principal material axes:
/// 1 = fibre direction, 2 = transverse direction.
#[derive(Debug, Clone)]
pub struct PlyLayer {
    /// Young's modulus along fibre direction \[Pa\].
    pub e1: f64,
    /// Young's modulus transverse to fibres \[Pa\].
    pub e2: f64,
    /// In-plane shear modulus \[Pa\].
    pub g12: f64,
    /// Major Poisson's ratio ν₁₂.
    pub nu12: f64,
    /// Ply thickness \[m\].
    pub thickness: f64,
}
impl PlyLayer {
    /// Create a new ply layer.
    pub fn new(e1: f64, e2: f64, g12: f64, nu12: f64, thickness: f64) -> Self {
        Self {
            e1,
            e2,
            g12,
            nu12,
            thickness,
        }
    }
    /// Reduced stiffness matrix Q (3×3) in principal material axes.
    ///
    /// Q = \[ Q11 Q12  0  \]
    ///     \[ Q12 Q22  0  \]
    ///     \[  0   0  Q66 \]
    pub fn q_matrix(&self) -> [[f64; 3]; 3] {
        let nu21 = self.nu12 * self.e2 / self.e1;
        let denom = 1.0 - self.nu12 * nu21;
        let q11 = self.e1 / denom;
        let q22 = self.e2 / denom;
        let q12 = self.nu12 * self.e2 / denom;
        let q66 = self.g12;
        [[q11, q12, 0.0], [q12, q22, 0.0], [0.0, 0.0, q66]]
    }
    /// Transformed reduced stiffness matrix Q̄ (3×3) at fibre angle θ (radians).
    ///
    /// Applies the standard tensor rotation to bring Q into global (x,y) axes.
    pub fn q_matrix_global(&self, theta: f64) -> [[f64; 3]; 3] {
        let q = self.q_matrix();
        let c = theta.cos();
        let s = theta.sin();
        let c2 = c * c;
        let s2 = s * s;
        let c4 = c2 * c2;
        let s4 = s2 * s2;
        let c2s2 = c2 * s2;
        let q11 = q[0][0];
        let q22 = q[1][1];
        let q12 = q[0][1];
        let q66 = q[2][2];
        [
            [
                q11 * c4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * s4,
                (q11 + q22 - 4.0 * q66) * c2s2 + q12 * (c4 + s4),
                (q11 - q12 - 2.0 * q66) * s * c.powi(3) - (q22 - q12 - 2.0 * q66) * c * s.powi(3),
            ],
            [
                (q11 + q22 - 4.0 * q66) * c2s2 + q12 * (c4 + s4),
                q11 * s4 + 2.0 * (q12 + 2.0 * q66) * c2s2 + q22 * c4,
                (q11 - q12 - 2.0 * q66) * c * s.powi(3) - (q22 - q12 - 2.0 * q66) * s * c.powi(3),
            ],
            [
                (q11 - q12 - 2.0 * q66) * s * c.powi(3) - (q22 - q12 - 2.0 * q66) * c * s.powi(3),
                (q11 - q12 - 2.0 * q66) * c * s.powi(3) - (q22 - q12 - 2.0 * q66) * s * c.powi(3),
                (q11 + q22 - 2.0 * q12 - 2.0 * q66) * c2s2 + q66 * (c4 + s4),
            ],
        ]
    }
}
/// Thin conical shell (apex at origin, generator at half-apex angle α).
#[derive(Debug, Clone)]
pub struct ConicalShell {
    /// Young's modulus \[Pa\].
    pub young_modulus: f64,
    /// Poisson's ratio.
    pub poisson_ratio: f64,
    /// Mass density \[kg/m³\].
    pub density: f64,
    /// Shell thickness \[m\].
    pub thickness: f64,
    /// Half-apex angle \[radians\].
    pub half_angle: f64,
    /// Base radius at the small end \[m\].
    pub r0: f64,
}
impl ConicalShell {
    /// Create a new conical shell.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        density: f64,
        thickness: f64,
        half_angle: f64,
        r0: f64,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            density,
            thickness,
            half_angle,
            r0,
        }
    }
    /// Slant length between two meridional positions s1 and s2 (measured along the slant).
    pub fn slant_length(&self, s1: f64, s2: f64) -> f64 {
        (s2 - s1).abs()
    }
    /// Radius at slant position s from the small end: R(s) = r0 + s * sin(α).
    pub fn mean_radius_at(&self, s: f64) -> f64 {
        self.r0 + s * self.half_angle.sin()
    }
    /// Meridional membrane stress due to internal pressure at slant position s.
    ///
    /// σ_s = p * r(s) / (2 * t * cos(α))
    pub fn meridional_stress_pressure(&self, p: f64, s: f64) -> f64 {
        let r = self.mean_radius_at(s);
        p * r / (2.0 * self.thickness * self.half_angle.cos())
    }
    /// Hoop membrane stress at slant position s.
    ///
    /// σ_θ = p * r(s) / (t * cos(α))
    pub fn hoop_stress_pressure(&self, p: f64, s: f64) -> f64 {
        let r = self.mean_radius_at(s);
        p * r / (self.thickness * self.half_angle.cos())
    }
}
