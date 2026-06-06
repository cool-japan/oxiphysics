//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

/// Extended piezoelectric constitutive relations.
///
/// In matrix form:
/// ```text
/// {T} = [C^E] {S} - [e]^T {E}
/// {D} = [e] {S} + [ε^S] {E}
/// ```
/// where C^E is stiffness at constant E, ε^S is permittivity at constant S.
pub struct PiezoCoupling {
    /// Coupling constant e33 (C/m²).
    pub e33: f64,
    /// Coupling constant e31 (C/m²).
    pub e31: f64,
    /// Coupling constant e15 (C/m²).
    pub e15: f64,
    /// Young's modulus (Pa).
    pub e_mod: f64,
    /// Poisson ratio.
    pub nu: f64,
}
impl PiezoCoupling {
    /// Create from material parameters.
    pub fn new(e_mod: f64, nu: f64, e33: f64, e31: f64, e15: f64) -> Self {
        Self {
            e33,
            e31,
            e15,
            e_mod,
            nu,
        }
    }
    /// Compute the d33 piezoelectric strain coefficient: d33 = e33 / C33.
    pub fn d33(&self) -> f64 {
        let c33 = self.e_mod * (1.0 - self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu));
        if c33.abs() > 1e-30 {
            self.e33 / c33
        } else {
            0.0
        }
    }
    /// Compute the d31 piezoelectric strain coefficient.
    pub fn d31(&self) -> f64 {
        let lam = self.e_mod * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu));
        let c33 = lam + 2.0 * self.e_mod / (2.0 * (1.0 + self.nu));
        if c33.abs() > 1e-30 {
            self.e31 / c33
        } else {
            0.0
        }
    }
    /// Compute the coupling factor k33 = e33 / sqrt(C33 * eps33).
    pub fn k33(&self, eps33: f64) -> f64 {
        let c33 = self.e_mod * (1.0 - self.nu) / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu));
        let denom = (c33 * eps33).max(0.0).sqrt();
        if denom < 1e-30 {
            return 0.0;
        }
        self.e33 / denom
    }
    /// Compute the induced strain for a given electric field E3 (V/m).
    ///
    /// ε33 = d33 * E3
    pub fn induced_strain_33(&self, e_field: f64) -> f64 {
        self.d33() * e_field
    }
    /// Compute the open-circuit voltage for a given applied stress T33 (Pa).
    ///
    /// Voc = -d33 * T33 * L / eps33
    pub fn open_circuit_voltage(&self, stress: f64, length: f64, eps33: f64) -> f64 {
        -self.d33() * stress * length / eps33.max(1e-30)
    }
}
/// FEM node carrying both mechanical and electrical degrees of freedom.
pub struct ElectromechanicalNode {
    /// X coordinate.
    pub x: f64,
    /// Y coordinate.
    pub y: f64,
    /// Z coordinate.
    pub z: f64,
    /// Mechanical DOF – displacement in X.
    pub u: f64,
    /// Mechanical DOF – displacement in Y.
    pub v: f64,
    /// Mechanical DOF – displacement in Z.
    pub w: f64,
    /// Electric potential DOF.
    pub phi: f64,
}
impl ElectromechanicalNode {
    /// Create a new node at position `(x, y, z)` with all DOFs zeroed.
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        ElectromechanicalNode {
            x,
            y,
            z,
            u: 0.0,
            v: 0.0,
            w: 0.0,
            phi: 0.0,
        }
    }
}
/// Simple 2-D electrostatic solver on a uniform Cartesian grid using Gauss–
/// Seidel iterations to solve the Laplace (or Poisson) equation
/// `∇²φ = -ρ/ε₀`.
pub struct ElectrostaticSolver {
    /// Number of grid points in x.
    pub nx: usize,
    /// Number of grid points in y.
    pub ny: usize,
    /// Nodal potential values (row-major: index = iy * nx + ix).
    pub phi: Vec<f64>,
    /// Volumetric charge density ρ at each node.
    pub charge_density: Vec<f64>,
    /// Grid spacing (uniform, square grid).
    pub dx: f64,
}
impl ElectrostaticSolver {
    /// Create a new solver on an `nx × ny` grid with spacing `dx`.
    ///
    /// All potentials and charge densities are initialised to zero.
    pub fn new(nx: usize, ny: usize, dx: f64) -> Self {
        let n = nx * ny;
        ElectrostaticSolver {
            nx,
            ny,
            phi: vec![0.0; n],
            charge_density: vec![0.0; n],
            dx,
        }
    }
    /// Run Gauss–Seidel iterations until convergence or `max_iter` is reached.
    ///
    /// Boundary nodes (the outermost ring) are treated as Dirichlet (held fixed
    /// at their initial value).
    ///
    /// Returns the number of iterations performed.
    pub fn solve_laplace_gauss_seidel(&mut self, max_iter: usize, tol: f64) -> usize {
        let nx = self.nx;
        let ny = self.ny;
        let dx2 = self.dx * self.dx;
        const EPS0: f64 = 8.854_187_817e-12;
        for iter in 0..max_iter {
            let mut max_change = 0.0_f64;
            for iy in 1..ny - 1 {
                for ix in 1..nx - 1 {
                    let idx = iy * nx + ix;
                    let rhs = self.charge_density[idx] / EPS0 * dx2;
                    let new_phi = 0.25
                        * (self.phi[idx - 1]
                            + self.phi[idx + 1]
                            + self.phi[idx - nx]
                            + self.phi[idx + nx]
                            + rhs);
                    let change = (new_phi - self.phi[idx]).abs();
                    if change > max_change {
                        max_change = change;
                    }
                    self.phi[idx] = new_phi;
                }
            }
            if max_change < tol {
                return iter + 1;
            }
        }
        max_iter
    }
    /// Compute the electric field magnitude `|E| = |−∇φ|` at each interior
    /// node using central differences.  Boundary nodes are set to zero.
    pub fn electric_field_magnitude(&self) -> Vec<f64> {
        let nx = self.nx;
        let ny = self.ny;
        let inv2dx = 1.0 / (2.0 * self.dx);
        let mut mag = vec![0.0_f64; nx * ny];
        for iy in 1..ny - 1 {
            for ix in 1..nx - 1 {
                let idx = iy * nx + ix;
                let ex = -(self.phi[idx + 1] - self.phi[idx - 1]) * inv2dx;
                let ey = -(self.phi[idx + nx] - self.phi[idx - nx]) * inv2dx;
                mag[idx] = (ex * ex + ey * ey).sqrt();
            }
        }
        mag
    }
}
/// Piezoelectric sensor: converts mechanical strain to voltage.
pub struct PiezoSensor {
    /// Material properties.
    pub material: PiezoCoupling,
    /// Electrode area (m²).
    pub area: f64,
    /// Thickness of sensing layer (m).
    pub thickness: f64,
    /// Permittivity at constant stress ε^T (F/m).
    pub permittivity: f64,
    /// Load resistance (Ohm).
    pub load_resistance: f64,
}
impl PiezoSensor {
    /// Create a new piezoelectric sensor.
    pub fn new(
        e_mod: f64,
        nu: f64,
        e33: f64,
        e31: f64,
        e15: f64,
        area: f64,
        thickness: f64,
        eps_r: f64,
        r_load: f64,
    ) -> Self {
        const EPS0: f64 = 8.854_187_817e-12;
        Self {
            material: PiezoCoupling::new(e_mod, nu, e33, e31, e15),
            area,
            thickness,
            permittivity: eps_r * EPS0,
            load_resistance: r_load,
        }
    }
    /// Short-circuit charge for a given applied stress T33 (Pa).
    ///
    /// Q = d33 * T33 * A
    pub fn short_circuit_charge(&self, stress_33: f64) -> f64 {
        self.material.d33() * stress_33 * self.area
    }
    /// Open-circuit voltage for stress T33.
    ///
    /// V_oc = Q / C = d33 * T33 * t / ε^T
    pub fn open_circuit_voltage(&self, stress_33: f64) -> f64 {
        self.material.d33() * stress_33 * self.thickness / self.permittivity.max(1e-60)
    }
    /// Capacitance: C = ε^T * A / t.
    pub fn capacitance(&self) -> f64 {
        self.permittivity * self.area / self.thickness.max(1e-30)
    }
    /// Power output for sinusoidal stress at frequency f_hz.
    ///
    /// P = V²_rms / R = (V_oc / sqrt(2))² / R
    pub fn power_output(&self, stress_33_amplitude: f64, _f_hz: f64) -> f64 {
        let v_oc = self.open_circuit_voltage(stress_33_amplitude);
        let v_rms = v_oc / 2.0_f64.sqrt();
        v_rms * v_rms / self.load_resistance.max(1e-30)
    }
}
/// Piezoelectric constitutive law in d-form (compliance form).
///
/// ```text
/// {S} = [s^E]{T} + [d]^T{E}
/// {D} = [d]{T}   + [ε^T]{E}
/// ```
/// where s^E is compliance at constant E and ε^T is permittivity at constant T.
pub struct PiezoD {
    /// Compliance matrix s^E (6×6).
    pub compliance: [[f64; 6]; 6],
    /// Piezoelectric d-coefficients (3×6): rows = E direction, cols = Voigt stress.
    pub d_matrix: [[f64; 6]; 3],
    /// Permittivity at constant stress ε^T (3×3).
    pub permittivity_t: [[f64; 3]; 3],
}
impl PiezoD {
    /// Create from d33, d31, d15 coefficients and material constants.
    pub fn from_pzt(e_mod: f64, nu: f64, eps_r: f64, d33: f64, d31: f64, d15: f64) -> Self {
        let s11 = 1.0 / e_mod;
        let s12 = -nu / e_mod;
        let s44 = 2.0 * (1.0 + nu) / e_mod;
        let mut s = [[0.0_f64; 6]; 6];
        for (i, row) in s.iter_mut().enumerate().take(3) {
            for (j, val) in row.iter_mut().enumerate().take(3) {
                *val = if i == j { s11 } else { s12 };
            }
        }
        s[3][3] = s44;
        s[4][4] = s44;
        s[5][5] = s44;
        let mut d = [[0.0_f64; 6]; 3];
        d[2][0] = d31;
        d[2][1] = d31;
        d[2][2] = d33;
        d[1][3] = d15;
        d[0][4] = d15;
        const EPS0: f64 = 8.854_187_817e-12;
        let eps_val = eps_r * EPS0;
        let mut eps = [[0.0_f64; 3]; 3];
        eps[0][0] = eps_val;
        eps[1][1] = eps_val;
        eps[2][2] = eps_val;
        PiezoD {
            compliance: s,
            d_matrix: d,
            permittivity_t: eps,
        }
    }
    /// Compute strain from stress {T} and electric field {E}.
    /// S = s^E * T + d^T * E  (6-vector)
    pub fn strain(&self, stress: &[f64; 6], e_field: &[f64; 3]) -> [f64; 6] {
        let mut s = [0.0_f64; 6];
        for (i, si) in s.iter_mut().enumerate() {
            let s_comp: f64 = (0..6).map(|j| self.compliance[i][j] * stress[j]).sum();
            let d_comp: f64 = (0..3).map(|k| self.d_matrix[k][i] * e_field[k]).sum();
            *si = s_comp + d_comp;
        }
        s
    }
    /// Compute electric displacement from stress and electric field.
    /// D = d * T + ε^T * E  (3-vector)
    pub fn electric_displacement(&self, stress: &[f64; 6], e_field: &[f64; 3]) -> [f64; 3] {
        let mut d = [0.0_f64; 3];
        for (i, di) in d.iter_mut().enumerate() {
            let dt: f64 = (0..6).map(|j| self.d_matrix[i][j] * stress[j]).sum();
            let de: f64 = (0..3).map(|j| self.permittivity_t[i][j] * e_field[j]).sum();
            *di = dt + de;
        }
        d
    }
}
/// Linear tetrahedral element with 4 nodes and an associated material.
pub struct ElectromechanicalElement {
    /// Global node indices (must reference nodes in the assembler's node list).
    pub nodes: [usize; 4],
    /// Index into the assembler's material list.
    pub material_id: usize,
}
/// Electrostriction tensor (simplified isotropic model).
///
/// The electrostrictive strain is:
/// ε_ij = M_ijkl E_k E_l
///
/// For an isotropic material with two electrostriction coefficients M1, M2:
/// ε_ij = M1 * E_i * E_j + M2 * δ_ij * |E|²
pub struct ElectrostrictiveMaterial {
    /// Electrostriction coefficient M1 (m²/V²).
    pub m1: f64,
    /// Electrostriction coefficient M2 (m²/V²).
    pub m2: f64,
}
impl ElectrostrictiveMaterial {
    /// Create a new electrostrictive material.
    pub fn new(m1: f64, m2: f64) -> Self {
        Self { m1, m2 }
    }
    /// Compute the electrostrictive strain tensor (3×3) for a given E field.
    pub fn strain_tensor(&self, e_field: &[f64; 3]) -> [[f64; 3]; 3] {
        let e_sq: f64 = e_field.iter().map(|e| e * e).sum();
        let mut eps = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                eps[i][j] = self.m1 * e_field[i] * e_field[j] + self.m2 * delta * e_sq;
            }
        }
        eps
    }
    /// Electrostrictive energy density: U_e = (M1/2)|E|⁴ (simplified).
    pub fn energy_density(&self, e_field: &[f64; 3]) -> f64 {
        let e_sq: f64 = e_field.iter().map(|e| e * e).sum();
        0.5 * (self.m1 + 3.0 * self.m2) * e_sq * e_sq
    }
}
/// Dielectric elastomer actuator (DEA) model.
///
/// A soft dielectric membrane that deforms under electric field.
/// Uses the neo-Hookean + electrostatic energy formulation.
pub struct DielectricElastomer {
    /// Shear modulus (Pa).
    pub mu: f64,
    /// Bulk modulus (Pa).
    pub kappa: f64,
    /// Permittivity (F/m).
    pub permittivity: f64,
    /// Thickness at rest (m).
    pub thickness: f64,
}
impl DielectricElastomer {
    /// Create a new DEA.
    pub fn new(mu: f64, kappa: f64, eps_r: f64, thickness: f64) -> Self {
        const EPS0: f64 = 8.854_187_817e-12;
        Self {
            mu,
            kappa,
            permittivity: eps_r * EPS0,
            thickness,
        }
    }
    /// Maxwell stress tensor for a given electric field E (3-vector).
    ///
    /// T_ij^{Maxwell} = ε (E_i E_j - δ_ij |E|²/2)
    pub fn maxwell_stress_tensor(&self, e_field: &[f64; 3]) -> [[f64; 3]; 3] {
        let e_sq: f64 = e_field.iter().map(|e| e * e).sum();
        let mut t = [[0.0_f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let delta = if i == j { 1.0 } else { 0.0 };
                t[i][j] = self.permittivity * (e_field[i] * e_field[j] - 0.5 * delta * e_sq);
            }
        }
        t
    }
    /// Equivalent pressure from the Maxwell stress in the z-direction.
    ///
    /// p_Maxwell = ε E_z² (compressive, drives thinning)
    pub fn maxwell_pressure(&self, voltage: f64) -> f64 {
        let e_z = voltage / self.thickness.max(1e-30);
        self.permittivity * e_z * e_z
    }
    /// Actuation strain in the thickness direction under Maxwell pressure.
    ///
    /// ε_z = -p_Maxwell / (2 * mu + kappa * ...)  (linearised)
    /// Simplified: ε_z ≈ -ε E² / (4 mu)
    pub fn actuation_strain(&self, voltage: f64) -> f64 {
        let p = self.maxwell_pressure(voltage);
        -p / (4.0 * self.mu.max(1e-30))
    }
    /// Electrostatic energy stored in the DEA per unit volume.
    ///
    /// U_e = ε E² / 2
    pub fn electrostatic_energy_density(&self, voltage: f64) -> f64 {
        let e = voltage / self.thickness.max(1e-30);
        0.5 * self.permittivity * e * e
    }
}
/// Piezoelectric material constants in Voigt notation.
///
/// The constitutive relation is:
///
/// ```text
/// {T} = [C]{S} - [e]^T {E}
/// {D} = [e]{S} + [eps_s]{E}
/// ```
pub struct PiezoMaterial {
    /// Elastic stiffness matrix C (6×6, Voigt notation).
    pub elastic_stiffness: [[f64; 6]; 6],
    /// Piezoelectric coupling matrix e (6×3): relates strain to electric
    /// displacement.  (Row index = Voigt strain component, col = field dir.)
    pub piezo_coupling: [[f64; 3]; 6],
    /// Permittivity at constant strain ε_s (3×3).
    pub permittivity: [[f64; 3]; 3],
    /// Mass density (kg/m³).
    pub density: f64,
}
impl PiezoMaterial {
    /// Build a transversely-isotropic piezoelectric material from the minimal
    /// parameters `E` (Young's modulus), `nu` (Poisson ratio), `eps_r`
    /// (relative permittivity), and the piezoelectric stress coefficients
    /// `e33`, `e31`, `e15` (C/m²).
    pub fn isotropic_with_piezo(
        e_mod: f64,
        nu: f64,
        eps_r: f64,
        e33: f64,
        e31: f64,
        e15: f64,
    ) -> Self {
        let lam = e_mod * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e_mod / (2.0 * (1.0 + nu));
        let mut c = [[0.0_f64; 6]; 6];
        for (i, row) in c.iter_mut().enumerate().take(3) {
            for val in row.iter_mut().take(3) {
                *val = lam;
            }
            row[i] = lam + 2.0 * mu;
        }
        c[3][3] = mu;
        c[4][4] = mu;
        c[5][5] = mu;
        let mut e = [[0.0_f64; 3]; 6];
        e[0][2] = e31;
        e[1][2] = e31;
        e[2][2] = e33;
        e[3][1] = e15;
        e[4][0] = e15;
        const EPS0: f64 = 8.854_187_817e-12;
        let eps_val = eps_r * EPS0;
        let mut eps = [[0.0_f64; 3]; 3];
        eps[0][0] = eps_val;
        eps[1][1] = eps_val;
        eps[2][2] = eps_val;
        PiezoMaterial {
            elastic_stiffness: c,
            piezo_coupling: e,
            permittivity: eps,
            density: 7750.0,
        }
    }
}
/// Simplified 1-D piezoelectric actuator model (poled along the z-axis).
pub struct PiezoActuator {
    /// Total length of the actuator in the poling direction (m).
    pub length: f64,
    /// Cross-sectional area (m²).
    pub area: f64,
    /// Material properties.
    pub material: PiezoMaterial,
    /// Applied voltage (V).
    pub voltage: f64,
}
impl PiezoActuator {
    /// Blocked force (force when tip displacement is zero).
    ///
    /// `F_blocked = e33 * (V / L) * A`
    pub fn blocked_force(&self) -> f64 {
        let e33 = self.material.piezo_coupling[2][2];
        let field = self.voltage / self.length;
        e33 * field * self.area
    }
    /// Free strain (strain when no mechanical load is applied).
    ///
    /// `d33 = e33 / C33`, `ε_free = d33 * V / L`
    pub fn free_strain(&self) -> f64 {
        let e33 = self.material.piezo_coupling[2][2];
        let c33 = self.material.elastic_stiffness[2][2];
        let d33 = if c33.abs() > 1e-30 { e33 / c33 } else { 0.0 };
        d33 * self.voltage / self.length
    }
    /// Free tip displacement `δ = ε_free * L`.
    pub fn tip_displacement(&self) -> f64 {
        self.free_strain() * self.length
    }
}
/// Assembles the coupled electromechanical stiffness matrix for a mesh of
/// linear tetrahedral piezoelectric elements.
pub struct PiezoFemAssembler {
    /// All nodes in the mesh.
    pub nodes: Vec<ElectromechanicalNode>,
    /// All elements.
    pub elements: Vec<ElectromechanicalElement>,
    /// Material database.
    pub materials: Vec<PiezoMaterial>,
}
impl PiezoFemAssembler {
    /// Create an empty assembler.
    pub fn new() -> Self {
        PiezoFemAssembler {
            nodes: Vec::new(),
            elements: Vec::new(),
            materials: Vec::new(),
        }
    }
    /// Compute the strain-displacement matrix **B** (6×12) for a linear
    /// tetrahedral element.
    ///
    /// `nodes` – coordinates of the four element nodes: `[[x,y,z\]; 4]`.
    ///
    /// The 12 mechanical DOFs are ordered as `[u1,v1,w1, u2,v2,w2, u3,v3,w3,
    /// u4,v4,w4]`.
    pub fn compute_b_matrix(nodes: &[[f64; 3]; 4]) -> [[f64; 12]; 6] {
        let [p0, p1, p2, p3] = nodes;
        let j = [
            [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
            [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
            [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
        ];
        let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
        let inv_det = 1.0 / det;
        let ji = [
            [
                (j[1][1] * j[2][2] - j[1][2] * j[2][1]) * inv_det,
                (j[0][2] * j[2][1] - j[0][1] * j[2][2]) * inv_det,
                (j[0][1] * j[1][2] - j[0][2] * j[1][1]) * inv_det,
            ],
            [
                (j[1][2] * j[2][0] - j[1][0] * j[2][2]) * inv_det,
                (j[0][0] * j[2][2] - j[0][2] * j[2][0]) * inv_det,
                (j[0][2] * j[1][0] - j[0][0] * j[1][2]) * inv_det,
            ],
            [
                (j[1][0] * j[2][1] - j[1][1] * j[2][0]) * inv_det,
                (j[0][1] * j[2][0] - j[0][0] * j[2][1]) * inv_det,
                (j[0][0] * j[1][1] - j[0][1] * j[1][0]) * inv_det,
            ],
        ];
        let grad_ref: [[f64; 3]; 4] = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut dn = [[0.0_f64; 3]; 4];
        for i in 0..4 {
            for row in 0..3 {
                let mut s = 0.0;
                for col in 0..3 {
                    s += ji[col][row] * grad_ref[i][col];
                }
                dn[i][row] = s;
            }
        }
        let mut b = [[0.0_f64; 12]; 6];
        for (i, &[dndx, dndy, dndz]) in dn.iter().enumerate() {
            let col_u = 3 * i;
            let col_v = 3 * i + 1;
            let col_w = 3 * i + 2;
            b[0][col_u] = dndx;
            b[1][col_v] = dndy;
            b[2][col_w] = dndz;
            b[3][col_v] = dndz;
            b[3][col_w] = dndy;
            b[4][col_u] = dndz;
            b[4][col_w] = dndx;
            b[5][col_u] = dndy;
            b[5][col_v] = dndx;
        }
        b
    }
    /// Compute the volume of a linear tetrahedron: `|det(J)| / 6`.
    pub fn compute_element_volume(nodes: &[[f64; 3]; 4]) -> f64 {
        let [p0, p1, p2, p3] = nodes;
        let j = [
            [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
            [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
            [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
        ];
        let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
        det.abs() / 6.0
    }
    /// Compute the combined 16×16 element stiffness matrix for element
    /// `elem_idx`.
    ///
    /// The DOF ordering is `[u1,v1,w1,φ1, u2,v2,w2,φ2, …, u4,v4,w4,φ4]`
    /// (4 nodes × 4 DOFs per node = 16).
    ///
    /// The sub-blocks are:
    /// - **K_uu** (12×12): mechanical stiffness `∫ B^T C B dV`
    /// - **K_up** (12×4): mechanical-electrical coupling `∫ B^T e B_φ dV`
    /// - **K_pp** (4×4): electric stiffness `∫ B_φ^T ε_s B_φ dV`
    pub fn compute_element_stiffness(&self, elem_idx: usize) -> Vec<Vec<f64>> {
        let elem = &self.elements[elem_idx];
        let mat = &self.materials[elem.material_id];
        let mut coords = [[0.0_f64; 3]; 4];
        for (i, &ni) in elem.nodes.iter().enumerate() {
            let n = &self.nodes[ni];
            coords[i] = [n.x, n.y, n.z];
        }
        let vol = Self::compute_element_volume(&coords);
        let b = Self::compute_b_matrix(&coords);
        let [p0, p1, p2, p3] = &coords;
        let jac = [
            [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
            [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
            [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
        ];
        let det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1])
            - jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0])
            + jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);
        let inv_det = 1.0 / det;
        let ji = [
            [
                (jac[1][1] * jac[2][2] - jac[1][2] * jac[2][1]) * inv_det,
                (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) * inv_det,
                (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) * inv_det,
            ],
            [
                (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) * inv_det,
                (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) * inv_det,
                (jac[0][2] * jac[1][0] - jac[0][0] * jac[1][2]) * inv_det,
            ],
            [
                (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]) * inv_det,
                (jac[0][1] * jac[2][0] - jac[0][0] * jac[2][1]) * inv_det,
                (jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0]) * inv_det,
            ],
        ];
        let grad_ref: [[f64; 3]; 4] = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut bp = [[0.0_f64; 4]; 3];
        for i in 0..4 {
            for row in 0..3 {
                let mut s = 0.0;
                for col in 0..3 {
                    s += ji[col][row] * grad_ref[i][col];
                }
                bp[row][i] = s;
            }
        }
        let mut cb = [[0.0_f64; 12]; 6];
        for (r, _) in (0..6usize).enumerate() {
            for c in 0..12 {
                let mut s = 0.0;
                for (k, _) in (0..6usize).enumerate() {
                    s += mat.elastic_stiffness[r][k] * b[k][c];
                }
                cb[r][c] = s;
            }
        }
        let mut kuu = vec![vec![0.0_f64; 12]; 12];
        for r in 0..12 {
            for c in 0..12 {
                let mut s = 0.0;
                for k in 0..6 {
                    s += b[k][r] * cb[k][c];
                }
                kuu[r][c] = s * vol;
            }
        }
        let mut ebp = [[0.0_f64; 4]; 6];
        for (r, _) in (0..6usize).enumerate() {
            for c in 0..4 {
                let mut s = 0.0;
                for (k, _) in (0..3usize).enumerate() {
                    s += mat.piezo_coupling[r][k] * bp[k][c];
                }
                ebp[r][c] = s;
            }
        }
        let mut kup = vec![vec![0.0_f64; 4]; 12];
        for r in 0..12 {
            for c in 0..4 {
                let mut s = 0.0;
                for k in 0..6 {
                    s += b[k][r] * ebp[k][c];
                }
                kup[r][c] = s * vol;
            }
        }
        let mut eps_bp = [[0.0_f64; 4]; 3];
        for (r, _) in (0..3usize).enumerate() {
            for c in 0..4 {
                let mut s = 0.0;
                for (k, _) in (0..3usize).enumerate() {
                    s += mat.permittivity[r][k] * bp[k][c];
                }
                eps_bp[r][c] = s;
            }
        }
        let mut kpp = vec![vec![0.0_f64; 4]; 4];
        for r in 0..4 {
            for c in 0..4 {
                let mut s = 0.0;
                for k in 0..3 {
                    s += bp[k][r] * eps_bp[k][c];
                }
                kpp[r][c] = s * vol;
            }
        }
        let mut ke = vec![vec![0.0_f64; 16]; 16];
        for r in 0..12 {
            for c in 0..12 {
                ke[r][c] = kuu[r][c];
            }
        }
        let mech_to_global = |i: usize| -> usize { (i / 3) * 4 + (i % 3) };
        let phi_to_global = |j: usize| -> usize { j * 4 + 3 };
        let mut ke16 = vec![vec![0.0_f64; 16]; 16];
        for r in 0..12 {
            for c in 0..12 {
                ke16[mech_to_global(r)][mech_to_global(c)] = kuu[r][c];
            }
        }
        for r in 0..12 {
            for c in 0..4 {
                ke16[mech_to_global(r)][phi_to_global(c)] = kup[r][c];
                ke16[phi_to_global(c)][mech_to_global(r)] = kup[r][c];
            }
        }
        for r in 0..4 {
            for c in 0..4 {
                ke16[phi_to_global(r)][phi_to_global(c)] = -kpp[r][c];
            }
        }
        let _ = ke;
        ke16
    }
}

impl Default for PiezoFemAssembler {
    fn default() -> Self {
        Self::new()
    }
}

/// A 3D tetrahedral element for electrostatic analysis.
///
/// Each element solves ∇·(ε ∇φ) = -ρ using standard Galerkin FEM.
pub struct ElectrostaticElement {
    /// Node coordinates: \[\[x,y,z\\]; 4].
    pub coords: [[f64; 3]; 4],
    /// Permittivity of the material (F/m).
    pub permittivity: f64,
}
impl ElectrostaticElement {
    /// Compute the 4×4 element stiffness matrix for Laplace/Poisson FEM.
    ///
    /// K_e = ε * ∫ B_φ^T B_φ dV = ε * vol * B_φ^T B_φ
    pub fn element_stiffness(&self) -> [[f64; 4]; 4] {
        let vol = PiezoFemAssembler::compute_element_volume(&self.coords);
        let bp = self.b_phi_matrix();
        let mut k = [[0.0f64; 4]; 4];
        for r in 0..4 {
            for c in 0..4 {
                let mut s = 0.0;
                for (k_idx, bprow) in bp.iter().enumerate().take(3) {
                    s += bprow[r] * bprow[c];
                    let _ = k_idx;
                }
                k[r][c] = self.permittivity * s * vol;
            }
        }
        k
    }
    /// Compute the electric field vector \[Ex, Ey, Ez\] from nodal potentials.
    ///
    /// E = -∇φ = -B_φ · {φ}
    pub fn electric_field(&self, nodal_phi: &[f64; 4]) -> [f64; 3] {
        let bp = self.b_phi_matrix();
        let mut grad = [0.0f64; 3];
        for row in 0..3 {
            for col in 0..4 {
                grad[row] += bp[row][col] * nodal_phi[col];
            }
        }
        [-grad[0], -grad[1], -grad[2]]
    }
    /// Compute B_φ: the 3×4 gradient-of-potential matrix for this element.
    fn b_phi_matrix(&self) -> [[f64; 4]; 3] {
        let [p0, p1, p2, p3] = &self.coords;
        let j = [
            [p1[0] - p0[0], p2[0] - p0[0], p3[0] - p0[0]],
            [p1[1] - p0[1], p2[1] - p0[1], p3[1] - p0[1]],
            [p1[2] - p0[2], p2[2] - p0[2], p3[2] - p0[2]],
        ];
        let det = j[0][0] * (j[1][1] * j[2][2] - j[1][2] * j[2][1])
            - j[0][1] * (j[1][0] * j[2][2] - j[1][2] * j[2][0])
            + j[0][2] * (j[1][0] * j[2][1] - j[1][1] * j[2][0]);
        let inv_det = 1.0 / det;
        let ji = [
            [
                (j[1][1] * j[2][2] - j[1][2] * j[2][1]) * inv_det,
                (j[0][2] * j[2][1] - j[0][1] * j[2][2]) * inv_det,
                (j[0][1] * j[1][2] - j[0][2] * j[1][1]) * inv_det,
            ],
            [
                (j[1][2] * j[2][0] - j[1][0] * j[2][2]) * inv_det,
                (j[0][0] * j[2][2] - j[0][2] * j[2][0]) * inv_det,
                (j[0][2] * j[1][0] - j[0][0] * j[1][2]) * inv_det,
            ],
            [
                (j[1][0] * j[2][1] - j[1][1] * j[2][0]) * inv_det,
                (j[0][1] * j[2][0] - j[0][0] * j[2][1]) * inv_det,
                (j[0][0] * j[1][1] - j[0][1] * j[1][0]) * inv_det,
            ],
        ];
        let grad_ref: [[f64; 3]; 4] = [
            [-1.0, -1.0, -1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let mut bp = [[0.0f64; 4]; 3];
        for i in 0..4 {
            for row in 0..3 {
                let mut s = 0.0;
                for col in 0..3 {
                    s += ji[col][row] * grad_ref[i][col];
                }
                bp[row][i] = s;
            }
        }
        bp
    }
}
/// Electric field components evaluated per element.
pub struct ElectricField {
    /// X-component of E per element.
    pub ex: Vec<f64>,
    /// Y-component of E per element.
    pub ey: Vec<f64>,
    /// Z-component of E per element.
    pub ez: Vec<f64>,
}
impl ElectricField {
    /// Create a zero electric field for `n` elements.
    pub fn zeros(n: usize) -> Self {
        ElectricField {
            ex: vec![0.0; n],
            ey: vec![0.0; n],
            ez: vec![0.0; n],
        }
    }
}
/// Dielectric breakdown criterion.
///
/// A material undergoes dielectric breakdown when |E| > E_breakdown.
#[derive(Debug, Clone)]
pub struct DielectricBreakdown {
    /// Breakdown electric field strength (V/m).
    pub e_breakdown: f64,
    /// Material relative permittivity.
    pub eps_r: f64,
}
impl DielectricBreakdown {
    /// Create a new breakdown model.
    pub fn new(e_breakdown: f64, eps_r: f64) -> Self {
        Self { e_breakdown, eps_r }
    }
    /// Check if the electric field magnitude exceeds the breakdown threshold.
    pub fn is_broken_down(&self, e_magnitude: f64) -> bool {
        e_magnitude >= self.e_breakdown
    }
    /// Safety margin: ratio of breakdown field to actual field.
    /// Values < 1.0 indicate breakdown.
    pub fn safety_margin(&self, e_magnitude: f64) -> f64 {
        if e_magnitude < 1e-30 {
            return f64::INFINITY;
        }
        self.e_breakdown / e_magnitude
    }
    /// Estimate energy density at breakdown: U = ε₀ εᵣ E_bd² / 2 (J/m³).
    pub fn energy_density_at_breakdown(&self) -> f64 {
        const EPS0: f64 = 8.854_187_817e-12;
        0.5 * EPS0 * self.eps_r * self.e_breakdown * self.e_breakdown
    }
}
/// MEMS electrostatic actuator model (parallel-plate approximation).
///
/// Models a capacitive actuator: two parallel plates of area A separated
/// by gap g.  An applied voltage V creates an electrostatic force.
#[derive(Debug, Clone)]
pub struct MemsActuator {
    /// Plate area (m²).
    pub area: f64,
    /// Initial (rest) gap (m).
    pub rest_gap: f64,
    /// Mechanical spring stiffness (N/m).
    pub spring_stiffness: f64,
    /// Permittivity of the gap medium (F/m).
    pub permittivity: f64,
}
impl MemsActuator {
    /// Create a new MEMS actuator.
    pub fn new(area: f64, rest_gap: f64, spring_stiffness: f64, eps_r: f64) -> Self {
        const EPS0: f64 = 8.854_187_817e-12;
        Self {
            area,
            rest_gap,
            spring_stiffness,
            permittivity: eps_r * EPS0,
        }
    }
    /// Compute the electrostatic force: F = ε A V² / (2 g²).
    pub fn electrostatic_force(&self, voltage: f64, gap: f64) -> f64 {
        let g = gap.max(1e-20);
        self.permittivity * self.area * voltage * voltage / (2.0 * g * g)
    }
    /// Find equilibrium gap for a given applied voltage.
    ///
    /// Equilibrium: F_elec(g) = k * (g0 - g)
    ///
    /// The stable equilibrium (if it exists) lies in (g0/3, g0).
    /// Bisection searches for the larger-gap root.
    pub fn equilibrium_gap(&self, voltage: f64) -> Option<f64> {
        let g0 = self.rest_gap;
        let h = |g: f64| self.spring_stiffness * (g0 - g) - self.electrostatic_force(voltage, g);
        let g_lo = g0 * 0.334;
        let g_hi = g0 * (1.0 - 1e-9);
        let h_lo = h(g_lo);
        let h_hi = h(g_hi);
        if h_lo * h_hi > 0.0 {
            return None;
        }
        let (mut lo, mut hi) = if h_lo > 0.0 {
            (g_lo, g_hi)
        } else {
            (g_hi, g_lo)
        };
        for _ in 0..80 {
            let mid = 0.5 * (lo + hi);
            if h(mid) > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let g_eq = 0.5 * (lo + hi);
        if g_eq > 0.0 && g_eq < g0 {
            Some(g_eq)
        } else {
            None
        }
    }
    /// Pull-in voltage: V_pi = sqrt(8 k g0³ / (27 ε A)).
    pub fn pull_in_voltage(&self) -> f64 {
        let g0 = self.rest_gap;
        let num = 8.0 * self.spring_stiffness * g0.powi(3);
        let den = 27.0 * self.permittivity * self.area;
        if den.abs() < 1e-60 {
            return f64::INFINITY;
        }
        (num / den).sqrt()
    }
    /// Capacitance at a given gap: C = ε A / g.
    pub fn capacitance(&self, gap: f64) -> f64 {
        self.permittivity * self.area / gap.max(1e-20)
    }
}
/// Global assembler for the coupled piezoelectric FEM problem.
///
/// Assembles the full global stiffness matrix from a list of elements.
/// DOF ordering: `[u1,v1,w1,φ1, u2,v2,w2,φ2, …]` (4 DOFs per node).
pub struct GlobalPiezoAssembler {
    /// Inner element-level assembler.
    pub local: PiezoFemAssembler,
    /// Global stiffness matrix (n_dof × n_dof, dense).
    pub k_global: Vec<Vec<f64>>,
    /// Total number of DOFs (4 per node).
    pub n_dof: usize,
}
impl GlobalPiezoAssembler {
    /// Create a new global assembler from an existing local assembler.
    pub fn new(local: PiezoFemAssembler) -> Self {
        let n_nodes = local.nodes.len();
        let n_dof = n_nodes * 4;
        let k_global = vec![vec![0.0; n_dof]; n_dof];
        Self {
            local,
            k_global,
            n_dof,
        }
    }
    /// Assemble the global stiffness matrix from all elements.
    pub fn assemble(&mut self) {
        let n_elems = self.local.elements.len();
        for ei in 0..n_elems {
            let ke = self.local.compute_element_stiffness(ei);
            let nodes = self.local.elements[ei].nodes;
            for (a, &na) in nodes.iter().enumerate() {
                for da in 0..4 {
                    let global_a = na * 4 + da;
                    let local_a = a * 4 + da;
                    for (b, &nb) in nodes.iter().enumerate() {
                        for db in 0..4 {
                            let global_b = nb * 4 + db;
                            let local_b = b * 4 + db;
                            if global_a < self.n_dof && global_b < self.n_dof {
                                self.k_global[global_a][global_b] += ke[local_a][local_b];
                            }
                        }
                    }
                }
            }
        }
    }
    /// Apply a Dirichlet boundary condition: set DOF `idx` to value `val`.
    pub fn apply_dirichlet(&mut self, idx: usize, _val: f64) {
        if idx >= self.n_dof {
            return;
        }
        for j in 0..self.n_dof {
            self.k_global[idx][j] = 0.0;
            self.k_global[j][idx] = 0.0;
        }
        self.k_global[idx][idx] = 1.0;
    }
    /// Extract the mechanical sub-block K_uu from the global matrix.
    ///
    /// Returns K_uu as a dense matrix of size (3*n_nodes) × (3*n_nodes).
    pub fn extract_k_uu(&self) -> Vec<Vec<f64>> {
        let n_nodes = self.local.nodes.len();
        let n_mech = 3 * n_nodes;
        let mut kuu = vec![vec![0.0; n_mech]; n_mech];
        for na in 0..n_nodes {
            for da in 0..3 {
                let row_global = na * 4 + da;
                let row_mech = na * 3 + da;
                for nb in 0..n_nodes {
                    for db in 0..3 {
                        let col_global = nb * 4 + db;
                        let col_mech = nb * 3 + db;
                        kuu[row_mech][col_mech] = self.k_global[row_global][col_global];
                    }
                }
            }
        }
        kuu
    }
    /// Extract the electrical sub-block K_pp from the global matrix.
    ///
    /// Returns K_pp as a dense matrix of size n_nodes × n_nodes.
    pub fn extract_k_pp(&self) -> Vec<Vec<f64>> {
        let n_nodes = self.local.nodes.len();
        let mut kpp = vec![vec![0.0; n_nodes]; n_nodes];
        for (na, kpp_row) in kpp.iter_mut().enumerate().take(n_nodes) {
            let row_global = na * 4 + 3;
            for (nb, kpp_val) in kpp_row.iter_mut().enumerate().take(n_nodes) {
                let col_global = nb * 4 + 3;
                *kpp_val = self.k_global[row_global][col_global];
            }
        }
        kpp
    }
}
