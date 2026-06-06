//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::GAS_CONSTANT;

/// Acoustic pressure loading on a structural surface.
///
/// The acoustic pressure acts as a distributed traction on the wet surface:
/// f_struct = ∫ p · n dA
#[derive(Debug, Clone)]
pub struct AcousticPressureLoad {
    /// Node indices on the wet surface.
    pub node_indices: Vec<usize>,
    /// Outward structural normal at each node.
    pub normals: Vec<[f64; 3]>,
    /// Acoustic pressure at each node \[Pa\].
    pub pressure: Vec<f64>,
    /// Area weight (tributary area) at each node \[m²\].
    pub areas: Vec<f64>,
}
impl AcousticPressureLoad {
    /// Create a new pressure load with zero initial pressure.
    pub fn new(node_indices: Vec<usize>, normals: Vec<[f64; 3]>, areas: Vec<f64>) -> Self {
        let n = node_indices.len();
        Self {
            node_indices,
            normals,
            pressure: vec![0.0; n],
            areas,
        }
    }
    /// Set acoustic pressure at each interface node.
    pub fn set_pressure(&mut self, pressure: Vec<f64>) {
        assert_eq!(pressure.len(), self.node_indices.len());
        self.pressure = pressure;
    }
    /// Compute the equivalent nodal force vector \[N\] on the structure.
    ///
    /// f_i = p_i · n_i · A_i
    pub fn nodal_forces(&self) -> Vec<[f64; 3]> {
        self.node_indices
            .iter()
            .enumerate()
            .map(|(k, _)| {
                let p = self.pressure[k];
                let n = self.normals[k];
                let a = self.areas[k];
                [p * n[0] * a, p * n[1] * a, p * n[2] * a]
            })
            .collect()
    }
    /// Total force magnitude \[N\].
    pub fn total_force_magnitude(&self) -> f64 {
        let forces = self.nodal_forces();
        let fx: f64 = forces.iter().map(|f| f[0]).sum();
        let fy: f64 = forces.iter().map(|f| f[1]).sum();
        let fz: f64 = forces.iter().map(|f| f[2]).sum();
        (fx * fx + fy * fy + fz * fz).sqrt()
    }
}
/// Radiation damping matrix for an acoustic-structural interface.
///
/// Models the energy loss from the structure into the infinite fluid domain
/// using the plane-wave approximation: c_rad = ρ_f · c_f · A (per unit area).
#[derive(Debug, Clone)]
pub struct RadiationDamping {
    /// Fluid density \[kg/m³\].
    pub rho_fluid: f64,
    /// Speed of sound in fluid \[m/s\].
    pub c_fluid: f64,
    /// Interface area \[m²\].
    pub area: f64,
}
impl RadiationDamping {
    /// Create a new radiation damping model.
    pub fn new(rho_fluid: f64, c_fluid: f64, area: f64) -> Self {
        Self {
            rho_fluid,
            c_fluid,
            area,
        }
    }
    /// Radiation damping coefficient \[N·s/m\].
    pub fn damping_coeff(&self) -> f64 {
        self.rho_fluid * self.c_fluid * self.area
    }
    /// Damping force on structure: f_rad = −c_rad · v_n (normal velocity).
    pub fn damping_force(&self, normal_velocity: f64) -> f64 {
        -self.damping_coeff() * normal_velocity
    }
    /// Power radiated into the fluid \[W\].
    pub fn radiated_power(&self, normal_velocity: f64) -> f64 {
        self.damping_coeff() * normal_velocity.powi(2)
    }
}
/// 1-D moisture diffusion solver coupled with elastic body.
#[derive(Debug, Clone)]
pub struct HygroDiffusionSolver {
    /// Number of nodes.
    pub n_nodes: usize,
    /// Bar length \[m\].
    pub length: f64,
    /// Moisture concentration at each node \[kg/m³\].
    pub concentration: Vec<f64>,
    /// Material.
    pub material: HygroMechanicalMaterial,
    /// Mechanical displacement at each node \[m\].
    pub displacement: Vec<f64>,
    /// Hygrostress at each node \[Pa\].
    pub hygrostress: Vec<f64>,
}
impl HygroDiffusionSolver {
    /// Create a solver with uniform initial concentration.
    pub fn new(n_nodes: usize, length: f64, c0: f64, material: HygroMechanicalMaterial) -> Self {
        Self {
            n_nodes,
            length,
            concentration: vec![c0; n_nodes],
            material,
            displacement: vec![0.0; n_nodes],
            hygrostress: vec![0.0; n_nodes],
        }
    }
    /// Advance moisture diffusion by one time step using explicit Euler.
    ///
    /// ∂c/∂t = D · ∂²c/∂x²  (Fickian diffusion, 1-D, explicit)
    pub fn diffusion_step(&mut self, dt: f64, c_left: f64, c_right: f64) {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let d = self.material.diffusivity;
        let r = d * dt / (dx * dx);
        let c_old = self.concentration.clone();
        self.concentration[0] = c_left;
        self.concentration[self.n_nodes - 1] = c_right;
        for i in 1..self.n_nodes - 1 {
            self.concentration[i] = c_old[i] + r * (c_old[i + 1] - 2.0 * c_old[i] + c_old[i - 1]);
        }
    }
    /// Update hygrostress from the current concentration field.
    pub fn update_hygrostress(&mut self) {
        for i in 0..self.n_nodes {
            self.hygrostress[i] = self
                .material
                .hydrostatic_stress_from_concentration(0.0, self.concentration[i]);
        }
    }
    /// Total moisture content \[kg/m²\] by trapezoidal integration.
    pub fn total_moisture(&self) -> f64 {
        let dx = self.length / (self.n_nodes - 1) as f64;
        let mut total = 0.0;
        for i in 0..self.n_nodes - 1 {
            total += 0.5 * (self.concentration[i] + self.concentration[i + 1]) * dx;
        }
        total
    }
}
/// Electro-mechanical FEM node.
#[derive(Debug, Clone)]
pub struct ElectroMechNode {
    /// Position \[m\]
    pub pos: [f64; 2],
    /// Displacement \[m\]
    pub disp: [f64; 2],
    /// Electric potential \[V\]
    pub phi: f64,
}
impl ElectroMechNode {
    /// Create a new node.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            pos: [x, y],
            disp: [0.0; 2],
            phi: 0.0,
        }
    }
    /// Electric field from potential gradient (finite difference approximation).
    pub fn e_field_from_gradient(phi_plus: f64, phi_minus: f64, h: f64) -> f64 {
        -(phi_plus - phi_minus) / (2.0 * h)
    }
}
/// Block 2×2 monolithic system for coupled problems.
///
/// ┌ K_uu  K_up ┐ ┌ u ┐   ┌ f_u ┐
/// │ K_pu  K_pp │ │ p │ = │ f_p │
/// └            ┘ └   ┘   └     ┘
#[derive(Debug, Clone)]
pub struct MonolithicBlock2x2 {
    /// Block K_uu (n_u × n_u)
    pub k_uu: Vec<Vec<f64>>,
    /// Block K_up (n_u × n_p) — coupling
    pub k_up: Vec<Vec<f64>>,
    /// Block K_pu (n_p × n_u) — coupling transpose
    pub k_pu: Vec<Vec<f64>>,
    /// Block K_pp (n_p × n_p)
    pub k_pp: Vec<Vec<f64>>,
    /// RHS for mechanics
    pub f_u: Vec<f64>,
    /// RHS for second field
    pub f_p: Vec<f64>,
    /// Displacement solution
    pub u: Vec<f64>,
    /// Second field solution (pressure / potential / etc.)
    pub p: Vec<f64>,
    /// Mechanical DOF count
    pub n_u: usize,
    /// Second-field DOF count
    pub n_p: usize,
}
impl MonolithicBlock2x2 {
    /// Create a zero-initialized block system.
    pub fn new(n_u: usize, n_p: usize) -> Self {
        Self {
            k_uu: vec![vec![0.0; n_u]; n_u],
            k_up: vec![vec![0.0; n_p]; n_u],
            k_pu: vec![vec![0.0; n_u]; n_p],
            k_pp: vec![vec![0.0; n_p]; n_p],
            f_u: vec![0.0; n_u],
            f_p: vec![0.0; n_p],
            u: vec![0.0; n_u],
            p: vec![0.0; n_p],
            n_u,
            n_p,
        }
    }
    /// Assemble the full (n_u+n_p) × (n_u+n_p) matrix.
    pub fn full_matrix(&self) -> Vec<Vec<f64>> {
        let n = self.n_u + self.n_p;
        let mut mat = vec![vec![0.0; n]; n];
        for (i, mat_i) in mat.iter_mut().enumerate().take(self.n_u) {
            for (j, mat_ij) in mat_i.iter_mut().enumerate().take(self.n_u) {
                *mat_ij = self.k_uu[i][j];
            }
            for (j, mat_ij) in mat_i[self.n_u..].iter_mut().enumerate().take(self.n_p) {
                *mat_ij = self.k_up[i][j];
            }
        }
        for (i, mat_row) in mat[self.n_u..].iter_mut().enumerate().take(self.n_p) {
            for (j, mat_rj) in mat_row.iter_mut().enumerate().take(self.n_u) {
                *mat_rj = self.k_pu[i][j];
            }
            for (j, mat_rj) in mat_row[self.n_u..].iter_mut().enumerate().take(self.n_p) {
                *mat_rj = self.k_pp[i][j];
            }
        }
        mat
    }
    /// Full RHS vector \[f_u; f_p\].
    pub fn full_rhs(&self) -> Vec<f64> {
        let mut rhs = self.f_u.clone();
        rhs.extend_from_slice(&self.f_p);
        rhs
    }
    /// Solve using block Gauss-Seidel (one iteration for demonstration).
    ///
    /// u = K_uu⁻¹ · (f_u − K_up · p)
    /// p = K_pp⁻¹ · (f_p − K_pu · u)
    pub fn block_gauss_seidel_step(&mut self) {
        for i in 0..self.n_u {
            let mut rhs_i = self.f_u[i];
            for j in 0..self.n_p {
                rhs_i -= self.k_up[i][j] * self.p[j];
            }
            let diag = self.k_uu[i][i].max(1e-30);
            self.u[i] = rhs_i / diag;
        }
        for i in 0..self.n_p {
            let mut rhs_i = self.f_p[i];
            for j in 0..self.n_u {
                rhs_i -= self.k_pu[i][j] * self.u[j];
            }
            let diag = self.k_pp[i][i].max(1e-30);
            self.p[i] = rhs_i / diag;
        }
    }
    /// Residual norm of the full coupled system.
    pub fn residual_norm(&self) -> f64 {
        let mat = self.full_matrix();
        let rhs = self.full_rhs();
        let n = self.n_u + self.n_p;
        let sol: Vec<f64> = self.u.iter().chain(self.p.iter()).cloned().collect();
        let mut res = 0.0f64;
        for i in 0..n {
            let ax: f64 = (0..n).map(|j| mat[i][j] * sol[j]).sum();
            res += (ax - rhs[i]).powi(2);
        }
        res.sqrt()
    }
}
/// Temperature-dependent 3-D isotropic stiffness tensor (Voigt 6×6).
///
/// E and ν are evaluated at the current temperature.
#[derive(Debug, Clone)]
pub struct TempDependentStiffness {
    /// Elastic modulus function: E(T).
    pub e0: f64,
    /// Thermal degradation coefficient \[1/K\].
    pub beta: f64,
    /// Reference temperature \[K\].
    pub t_ref: f64,
    /// Poisson's ratio (assumed constant).
    pub nu: f64,
}
impl TempDependentStiffness {
    /// Create a new temperature-dependent stiffness.
    pub fn new(e0: f64, beta: f64, t_ref: f64, nu: f64) -> Self {
        Self {
            e0,
            beta,
            t_ref,
            nu,
        }
    }
    /// Elastic modulus at temperature T.
    pub fn modulus(&self, temp: f64) -> f64 {
        (self.e0 * (1.0 - self.beta * (temp - self.t_ref))).max(0.0)
    }
    /// Build the 6×6 isotropic stiffness matrix \[Pa\] at temperature T.
    ///
    /// Uses Voigt ordering: \[σ_xx, σ_yy, σ_zz, τ_xy, τ_yz, τ_zx\].
    pub fn stiffness_6x6(&self, temp: f64) -> [[f64; 6]; 6] {
        let e = self.modulus(temp);
        let nu = self.nu;
        let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let c11 = c * (1.0 - nu);
        let c12 = c * nu;
        let c44 = e / (2.0 * (1.0 + nu));
        let mut d = [[0.0f64; 6]; 6];
        for (i, d_row) in d.iter_mut().enumerate().take(3) {
            d_row[i] = c11;
        }
        for (i, d_row) in d.iter_mut().enumerate().take(3) {
            for (j, d_ij) in d_row.iter_mut().enumerate().take(3) {
                if i != j {
                    *d_ij = c12;
                }
            }
        }
        for (i, d_row) in d.iter_mut().enumerate().skip(3).take(3) {
            d_row[i] = c44;
        }
        d
    }
    /// Compute stress from total strain and temperature.
    ///
    /// σ = C(T)·(ε_total − ε_th(T))
    pub fn stress_from_total_strain(
        &self,
        total_strain: [f64; 6],
        temp: f64,
        alpha: f64,
    ) -> [f64; 6] {
        let mut state = ThermalExpansionState::new(self.t_ref, alpha);
        state.update(temp);
        let el = state.elastic_strain(total_strain);
        let c = self.stiffness_6x6(temp);
        let mut sigma = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                sigma[i] += c[i][j] * el[j];
            }
        }
        sigma
    }
}
/// Arbitrary Lagrangian-Eulerian (ALE) mesh node.
#[derive(Debug, Clone)]
pub struct AleNode {
    /// Material (reference) coordinates X
    pub material_pos: [f64; 2],
    /// Spatial (current) coordinates x
    pub spatial_pos: [f64; 2],
    /// Mesh velocity w (ALE frame velocity)
    pub mesh_velocity: [f64; 2],
    /// Fluid velocity v at this node
    pub fluid_velocity: [f64; 2],
    /// Fluid pressure p
    pub pressure: f64,
    /// Structural displacement u
    pub displacement: [f64; 2],
}
impl AleNode {
    /// Create a new ALE node at the given material coordinate.
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            material_pos: [x, y],
            spatial_pos: [x, y],
            mesh_velocity: [0.0; 2],
            fluid_velocity: [0.0; 2],
            pressure: 0.0,
            displacement: [0.0; 2],
        }
    }
    /// Convective velocity (fluid velocity relative to mesh motion).
    #[inline]
    pub fn convective_velocity(&self) -> [f64; 2] {
        [
            self.fluid_velocity[0] - self.mesh_velocity[0],
            self.fluid_velocity[1] - self.mesh_velocity[1],
        ]
    }
    /// Update spatial position from displacement.
    pub fn update_position(&mut self) {
        self.spatial_pos[0] = self.material_pos[0] + self.displacement[0];
        self.spatial_pos[1] = self.material_pos[1] + self.displacement[1];
    }
}
/// Eddy-current induced body force model for a conducting 1-D bar.
///
/// Models induction heating and body force due to a time-varying magnetic field.
#[derive(Debug, Clone)]
pub struct EddyCurrentModel {
    /// Number of elements in the bar.
    pub n_elem: usize,
    /// Bar length \[m\].
    pub length: f64,
    /// Electrical conductivity \[S/m\].
    pub sigma_e: f64,
    /// Magnetic permeability \[H/m\].
    pub mu_r: f64,
    /// Eddy current density at each element \[A/m²\].
    pub eddy_current: Vec<f64>,
    /// Joule heat source at each element \[W/m³\].
    pub heat_source: Vec<f64>,
}
impl EddyCurrentModel {
    /// Create a new eddy current model.
    pub fn new(n_elem: usize, length: f64, sigma_e: f64, mu_r: f64) -> Self {
        Self {
            n_elem,
            length,
            sigma_e,
            mu_r,
            eddy_current: vec![0.0; n_elem],
            heat_source: vec![0.0; n_elem],
        }
    }
    /// Compute eddy current from the time derivative of the applied B field.
    ///
    /// J_eddy ≈ σ · (dB/dt) · L / (2π)  (simplified circular cross-section).
    pub fn compute_eddy_currents(&mut self, db_dt: f64) {
        let dx = self.length / self.n_elem as f64;
        for i in 0..self.n_elem {
            let x_mid = (i as f64 + 0.5) * dx;
            self.eddy_current[i] = self.sigma_e * db_dt * x_mid;
            self.heat_source[i] = self.eddy_current[i].powi(2) / self.sigma_e.max(1e-30);
        }
    }
    /// Return the total Joule heating power \[W/m²\].
    pub fn total_joule_power(&self) -> f64 {
        let dx = self.length / self.n_elem as f64;
        self.heat_source.iter().sum::<f64>() * dx
    }
    /// Skin depth δ = √(2 / (ω·μ·σ)) \[m\].
    pub fn skin_depth(&self, angular_frequency: f64) -> f64 {
        let mu0 = 4.0 * std::f64::consts::PI * 1e-7;
        let mu = self.mu_r * mu0;
        (2.0 / (angular_frequency * mu * self.sigma_e).max(1e-60)).sqrt()
    }
}
/// Mortar coupling for non-conforming interface meshes.
///
/// Enforces displacement continuity across a non-matching interface using
/// Lagrange multipliers.
#[derive(Debug, Clone)]
pub struct MortarInterface {
    /// Mortar (master) side node positions \[m\]
    pub master_nodes: Vec<f64>,
    /// Non-mortar (slave) side node positions \[m\]
    pub slave_nodes: Vec<f64>,
    /// Mortar coupling matrix D (slave × slave)
    pub d_mat: Vec<Vec<f64>>,
    /// Mortar coupling matrix M (slave × master)
    pub m_mat: Vec<Vec<f64>>,
    /// Lagrange multipliers (traction at interface)
    pub lambdas: Vec<f64>,
}
impl MortarInterface {
    /// Create a 1D mortar interface.
    pub fn new_1d(master_nodes: Vec<f64>, slave_nodes: Vec<f64>) -> Self {
        let n_slave = slave_nodes.len();
        let n_master = master_nodes.len();
        let mut iface = Self {
            master_nodes,
            slave_nodes,
            d_mat: vec![vec![0.0; n_slave]; n_slave],
            m_mat: vec![vec![0.0; n_master]; n_slave],
            lambdas: vec![0.0; n_slave],
        };
        iface.compute_mortar_matrices();
        iface
    }
    /// Compute the D and M mortar matrices using Lobatto quadrature.
    fn compute_mortar_matrices(&mut self) {
        let n_slave = self.slave_nodes.len();
        let n_master = self.master_nodes.len();
        for i in 0..n_slave {
            self.d_mat[i][i] = if i == 0 || i == n_slave - 1 { 0.5 } else { 1.0 };
            if i > 0 {
                self.d_mat[i][i] *= self.slave_nodes[i] - self.slave_nodes[i - 1];
            }
        }
        for i in 0..n_slave {
            for j in 0..n_master {
                let dist = (self.slave_nodes[i] - self.master_nodes[j]).abs();
                let h = 1.0 / n_master as f64;
                self.m_mat[i][j] = (1.0 - dist / h).max(0.0) * h;
            }
        }
    }
    /// Compute the gap function g = u_slave − M · u_master.
    pub fn gap_function(&self, u_slave: &[f64], u_master: &[f64]) -> Vec<f64> {
        (0..self.slave_nodes.len())
            .map(|i| {
                let mu: f64 = (0..self.master_nodes.len())
                    .map(|j| self.m_mat[i][j] * u_master[j])
                    .sum();
                u_slave[i] - mu
            })
            .collect()
    }
    /// Enforce compatibility using penalty method: f_penalty = k_p · g.
    pub fn penalty_force(&self, u_slave: &[f64], u_master: &[f64], k_p: f64) -> Vec<f64> {
        let gap = self.gap_function(u_slave, u_master);
        gap.iter().map(|&g| -k_p * g).collect()
    }
    /// L² norm of the gap (compatibility error).
    pub fn compatibility_error(&self, u_slave: &[f64], u_master: &[f64]) -> f64 {
        let gap = self.gap_function(u_slave, u_master);
        gap.iter().map(|g| g * g).sum::<f64>().sqrt()
    }
}
/// Simplified ALE-FSI partitioned solver using Gauss-Seidel iteration.
#[derive(Debug, Clone)]
pub struct AleFsiSolver {
    /// All ALE nodes
    pub nodes: Vec<AleNode>,
    /// FSI interface
    pub interface: AleInterface,
    /// Structural mass density \[kg/m³\]
    pub struct_rho: f64,
    /// Structural thickness \[m\]
    pub thickness: f64,
    /// Fluid dynamic viscosity \[Pa·s\]
    pub mu: f64,
    /// Fluid density \[kg/m³\]
    pub rho_f: f64,
    /// Under-relaxation factor
    pub omega: f64,
}
impl AleFsiSolver {
    /// Create a new ALE-FSI solver.
    pub fn new(nodes: Vec<AleNode>, interface: AleInterface) -> Self {
        Self {
            nodes,
            interface,
            struct_rho: 7850.0,
            thickness: 0.01,
            mu: 1e-3,
            rho_f: 1000.0,
            omega: 0.5,
        }
    }
    /// One partitioned coupling step (loosely coupled, one sub-iteration).
    pub fn step(&mut self, dt: f64, inlet_velocity: f64) {
        for node in &mut self.nodes {
            node.fluid_velocity[0] = inlet_velocity;
            node.pressure = self.rho_f * 9.81 * (1.0 - node.spatial_pos[1]);
            node.pressure = node.pressure.max(0.0);
        }
        let segment_lengths =
            vec![1.0 / self.interface.node_indices.len() as f64; self.interface.node_indices.len()];
        let tractions = self
            .interface
            .compute_traction(&self.nodes, &segment_lengths);
        for (k, &idx) in self.interface.node_indices.iter().enumerate() {
            let f_y = tractions[k][1];
            let a_y = f_y / (self.struct_rho * self.thickness);
            self.nodes[idx].displacement[1] += 0.5 * a_y * dt * dt;
            self.nodes[idx].mesh_velocity[1] += a_y * dt;
        }
        for node in &mut self.nodes {
            node.update_position();
        }
    }
    /// Compute the total fluid force on the structure.
    pub fn total_fluid_force(&self) -> [f64; 2] {
        let segment_lengths =
            vec![1.0 / self.interface.node_indices.len() as f64; self.interface.node_indices.len()];
        let tractions = self
            .interface
            .compute_traction(&self.nodes, &segment_lengths);
        tractions
            .iter()
            .fold([0.0; 2], |acc, t| [acc[0] + t[0], acc[1] + t[1]])
    }
    /// Compute maximum structural deflection.
    pub fn max_deflection(&self) -> f64 {
        self.interface
            .node_indices
            .iter()
            .map(|&idx| self.nodes[idx].displacement[1].abs())
            .fold(0.0_f64, f64::max)
    }
}
/// Aitken dynamic under-relaxation for partitioned FSI coupling.
///
/// Computes the relaxation factor ω using the Aitken ΔΔ method:
/// ω_{k+1} = −ω_k · (rᵢ · (rᵢ − rᵢ₋₁)) / ‖rᵢ − rᵢ₋₁‖²
#[derive(Debug, Clone)]
pub struct AitkenRelaxation {
    /// Current relaxation factor
    pub omega: f64,
    /// Previous residual vector
    pub residual_prev: Vec<f64>,
    /// Current residual vector
    pub residual_curr: Vec<f64>,
    /// Previous solution
    pub solution_prev: Vec<f64>,
    /// Minimum allowed relaxation factor
    pub omega_min: f64,
    /// Maximum allowed relaxation factor
    pub omega_max: f64,
}
impl AitkenRelaxation {
    /// Create a new Aitken relaxation controller.
    pub fn new(n_dof: usize, omega0: f64) -> Self {
        Self {
            omega: omega0,
            residual_prev: vec![0.0; n_dof],
            residual_curr: vec![0.0; n_dof],
            solution_prev: vec![0.0; n_dof],
            omega_min: 0.001,
            omega_max: 1.0,
        }
    }
    /// Update the Aitken factor given the current and previous residuals.
    pub fn update_omega(&mut self) {
        let n = self.residual_curr.len();
        let mut num = 0.0f64;
        let mut denom = 0.0f64;
        for i in 0..n {
            let dr = self.residual_curr[i] - self.residual_prev[i];
            num += self.residual_prev[i] * dr;
            denom += dr * dr;
        }
        if denom.abs() > 1e-30 {
            self.omega = (-self.omega * num / denom).clamp(self.omega_min, self.omega_max);
        }
    }
    /// Apply relaxation to the interface solution update.
    ///
    /// x_{k+1} = x_k + ω · (x_tilde − x_k)
    pub fn relax(&self, x_old: &[f64], x_new_pred: &[f64]) -> Vec<f64> {
        x_old
            .iter()
            .zip(x_new_pred.iter())
            .map(|(&xo, &xn)| xo + self.omega * (xn - xo))
            .collect()
    }
    /// Check convergence: ‖r‖ / ‖x‖ < tol.
    pub fn is_converged(&self, tol: f64) -> bool {
        let r_norm: f64 = self.residual_curr.iter().map(|r| r * r).sum::<f64>().sqrt();
        let x_norm: f64 = self.solution_prev.iter().map(|x| x * x).sum::<f64>().sqrt();
        r_norm < tol * (1.0 + x_norm)
    }
    /// Number of degrees of freedom.
    pub fn n_dof(&self) -> usize {
        self.residual_curr.len()
    }
}
/// Lorentz force density on a conducting material carrying a current.
///
/// f = J × B  \[N/m³\]  where J is current density \[A/m²\] and B is magnetic flux \[T\].
#[derive(Debug, Clone)]
pub struct LorentzForceDensity {
    /// Applied magnetic flux density \[Bx, By, Bz\] \[T\].
    pub b_field: [f64; 3],
    /// Electrical conductivity of the material \[S/m\].
    pub conductivity: f64,
}
impl LorentzForceDensity {
    /// Create a new Lorentz force model.
    pub fn new(b_field: [f64; 3], conductivity: f64) -> Self {
        Self {
            b_field,
            conductivity,
        }
    }
    /// Compute the force density f = J × B given current density J \[A/m²\].
    pub fn force(&self, j_current: [f64; 3]) -> [f64; 3] {
        let b = self.b_field;
        let j = j_current;
        [
            j[1] * b[2] - j[2] * b[1],
            j[2] * b[0] - j[0] * b[2],
            j[0] * b[1] - j[1] * b[0],
        ]
    }
    /// Compute current density from electric field: J = σ · E.
    pub fn current_from_field(&self, e_field: [f64; 3]) -> [f64; 3] {
        [
            self.conductivity * e_field[0],
            self.conductivity * e_field[1],
            self.conductivity * e_field[2],
        ]
    }
    /// Ohmic power dissipation density: q = J · E = |J|² / σ \[W/m³\].
    pub fn ohmic_dissipation(&self, j_current: [f64; 3]) -> f64 {
        let j2 = j_current[0].powi(2) + j_current[1].powi(2) + j_current[2].powi(2);
        if self.conductivity > 0.0 {
            j2 / self.conductivity
        } else {
            0.0
        }
    }
    /// Compute magnetic pressure: p_mag = |B|² / (2·μ₀) \[Pa\].
    pub fn magnetic_pressure(&self) -> f64 {
        let mu0 = 4.0 * std::f64::consts::PI * 1e-7;
        let b2 = self.b_field[0].powi(2) + self.b_field[1].powi(2) + self.b_field[2].powi(2);
        b2 / (2.0 * mu0)
    }
}
/// A staggered coupling scheme that alternates between field solvers.
///
/// At each step:
///  1. Freeze field B, update field A.
///  2. Freeze field A, update field B.
///  3. Check coupling convergence.
#[derive(Debug, Clone)]
pub struct StaggeredCouplingScheme {
    /// Maximum number of coupling iterations per time step.
    pub max_coupling_iter: usize,
    /// Convergence tolerance on the coupling residual.
    pub coupling_tol: f64,
    /// Under-relaxation factor for field A (e.g., structural displacement).
    pub omega_a: f64,
    /// Under-relaxation factor for field B (e.g., temperature).
    pub omega_b: f64,
    /// History of coupling residuals.
    pub residual_history: Vec<f64>,
}
impl StaggeredCouplingScheme {
    /// Create a new staggered scheme with given parameters.
    pub fn new(max_iter: usize, tol: f64) -> Self {
        Self {
            max_coupling_iter: max_iter,
            coupling_tol: tol,
            omega_a: 0.8,
            omega_b: 0.8,
            residual_history: Vec::new(),
        }
    }
    /// Apply under-relaxation: x_new = x_old + ω·(x_iter − x_old).
    pub fn relax(&self, x_old: &[f64], x_iter: &[f64], omega: f64) -> Vec<f64> {
        x_old
            .iter()
            .zip(x_iter.iter())
            .map(|(xo, xi)| xo + omega * (xi - xo))
            .collect()
    }
    /// Check coupling convergence given old and new field values.
    ///
    /// Returns (converged, relative_change).
    pub fn check_convergence(&mut self, x_old: &[f64], x_new: &[f64]) -> (bool, f64) {
        let norm_old: f64 = x_old.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
        let diff: f64 = x_old
            .iter()
            .zip(x_new.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let rel = diff / norm_old;
        self.residual_history.push(rel);
        (rel < self.coupling_tol, rel)
    }
    /// Returns true if the scheme is converging (last residual < tol).
    pub fn is_converged(&self) -> bool {
        self.residual_history
            .last()
            .map(|&r| r < self.coupling_tol)
            .unwrap_or(false)
    }
    /// Number of coupling iterations taken.
    pub fn iter_count(&self) -> usize {
        self.residual_history.len()
    }
}
/// Piezoelectric coupling tensor (e-form, 3×6 in full 3D Voigt notation).
///
/// For simplicity we store the 2D plane-stress version: e\[3×3\].
#[derive(Debug, Clone)]
pub struct PiezoelectricCoupling {
    /// Piezoelectric coupling matrix e \[C/m²\] in reduced Voigt form
    pub e_mat: [[f64; 3]; 2],
    /// Permittivity tensor \[F/m\]
    pub epsilon: [f64; 2],
    /// Mechanical stiffness matrix at constant E-field \[Pa\]
    pub c_mat: [[f64; 3]; 3],
}
impl PiezoelectricCoupling {
    /// Create a standard PZT-5H-like coupling in plane stress.
    pub fn pzt5h() -> Self {
        Self {
            e_mat: [[0.0, 0.0, 17.0], [-6.5, 23.3, 0.0]],
            epsilon: [1.503e-8, 1.301e-8],
            c_mat: [
                [12.6e10, 8.41e10, 0.0],
                [8.41e10, 11.7e10, 0.0],
                [0.0, 0.0, 2.33e10],
            ],
        }
    }
    /// Compute mechanical stress from strain and electric field.
    ///
    /// σ = C · ε − eᵀ · E
    pub fn compute_stress(&self, strain: [f64; 3], e_field: [f64; 2]) -> [f64; 3] {
        let mut sigma = [0.0f64; 3];
        for (i, sigma_i) in sigma.iter_mut().enumerate() {
            for (j, &s_j) in strain.iter().enumerate() {
                *sigma_i += self.c_mat[i][j] * s_j;
            }
            for (k, &e_k) in e_field.iter().enumerate() {
                *sigma_i -= self.e_mat[k][i] * e_k;
            }
        }
        sigma
    }
    /// Compute electric displacement from strain and electric field.
    ///
    /// D = e · ε + ε_r · E
    pub fn compute_electric_displacement(&self, strain: [f64; 3], e_field: [f64; 2]) -> [f64; 2] {
        let mut d = [0.0f64; 2];
        for (k, d_k) in d.iter_mut().enumerate() {
            for (j, &s_j) in strain.iter().enumerate() {
                *d_k += self.e_mat[k][j] * s_j;
            }
            *d_k += self.epsilon[k] * e_field[k];
        }
        d
    }
    /// Electromechanical coupling coefficient k_eff.
    pub fn coupling_coefficient(&self) -> f64 {
        let e33 = self.e_mat[1][1];
        let eps = self.epsilon[1];
        let c33 = self.c_mat[1][1];
        if eps * c33 > 0.0 {
            (e33 * e33 / (eps * c33)).sqrt()
        } else {
            0.0
        }
    }
}
/// Thermo-mechanical staggered coupling driver.
///
/// Solves heat equation first, then uses resulting temperatures as input to the
/// mechanical problem, iterating until convergence.
#[derive(Debug, Clone)]
pub struct ThermoMechanicalCoupling {
    /// Thermal field
    pub thermal: ThermalField,
    /// Mechanical displacement field \[u_x, u_y\] per node
    pub displacements: Vec<[f64; 2]>,
    /// Material model
    pub material: ThermoElasticMaterial,
    /// Convergence tolerance for the staggered iteration
    pub tol: f64,
    /// Maximum iterations
    pub max_iter: usize,
}
impl ThermoMechanicalCoupling {
    /// Create a new thermo-mechanical coupling problem on an nx × ny mesh.
    pub fn new(nx: usize, ny: usize, material: ThermoElasticMaterial) -> Self {
        let n_nodes = nx * ny;
        let thermal = ThermalField::uniform(nx, ny, material.t_ref);
        let displacements = vec![[0.0; 2]; n_nodes];
        Self {
            thermal,
            displacements,
            material,
            tol: 1e-6,
            max_iter: 50,
        }
    }
    /// One staggered step: solve heat first, then mechanics.
    ///
    /// Returns the change in displacement norm (convergence indicator).
    pub fn staggered_step(&mut self, dt: f64, heat_source: f64) -> f64 {
        let nx = self.thermal.nx;
        let ny = self.thermal.ny;
        let alpha_diff = self.material.kappa / (self.material.rho * self.material.cp);
        let h = 1.0 / (nx.max(ny) as f64);
        let fo = alpha_diff * dt / (h * h);
        let mut new_temps = self.thermal.temperatures.clone();
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let c = self.thermal.temperatures[self.thermal.idx(i, j)];
                let xp = self.thermal.temperatures[self.thermal.idx(i + 1, j)];
                let xm = self.thermal.temperatures[self.thermal.idx(i - 1, j)];
                let yp = self.thermal.temperatures[self.thermal.idx(i, j + 1)];
                let ym = self.thermal.temperatures[self.thermal.idx(i, j - 1)];
                new_temps[self.thermal.idx(i, j)] = c
                    + fo * (xp + xm + yp + ym - 4.0 * c)
                    + heat_source * dt / (self.material.rho * self.material.cp);
            }
        }
        self.thermal.temperatures = new_temps;
        let mut disp_change = 0.0f64;
        for node in 0..self.displacements.len() {
            let temp = self.thermal.temperatures[node];
            let dt_t = temp - self.material.t_ref;
            let du = self.material.alpha * dt_t * h;
            disp_change += du * du;
            self.displacements[node][0] += du * 0.5;
            self.displacements[node][1] += du * 0.5;
        }
        disp_change.sqrt()
    }
    /// Run the staggered iteration until convergence.
    pub fn solve(&mut self, dt: f64, heat_source: f64) -> usize {
        for iter in 0..self.max_iter {
            let change = self.staggered_step(dt, heat_source);
            if change < self.tol {
                return iter + 1;
            }
        }
        self.max_iter
    }
    /// Compute the von Mises stress at a given node (plane-stress approximation).
    pub fn von_mises_stress(&self, node: usize) -> f64 {
        let temp = self.thermal.temperatures[node];
        let e = self.material.modulus_at(temp);
        let nu = self.material.nu;
        let ux = self.displacements[node][0];
        let uy = self.displacements[node][1];
        let eps_x = ux;
        let eps_y = uy;
        let sigma_x = e / (1.0 - nu * nu) * (eps_x + nu * eps_y);
        let sigma_y = e / (1.0 - nu * nu) * (nu * eps_x + eps_y);
        let tau_xy = 0.0_f64;
        (sigma_x * sigma_x - sigma_x * sigma_y + sigma_y * sigma_y + 3.0 * tau_xy * tau_xy).sqrt()
    }
}
/// Node temperatures for a 2-D quadrilateral mesh patch.
#[derive(Debug, Clone)]
pub struct ThermalField {
    /// Temperature at each node \[K\]
    pub temperatures: Vec<f64>,
    /// Number of nodes in x-direction
    pub nx: usize,
    /// Number of nodes in y-direction
    pub ny: usize,
}
impl ThermalField {
    /// Create a uniform temperature field.
    pub fn uniform(nx: usize, ny: usize, t: f64) -> Self {
        Self {
            temperatures: vec![t; nx * ny],
            nx,
            ny,
        }
    }
    /// Node index (row-major).
    #[inline]
    pub fn idx(&self, i: usize, j: usize) -> usize {
        i * self.ny + j
    }
    /// Bilinear interpolation of temperature at (ξ, η) ∈ \[0,1\]² within element (i,j).
    pub fn interpolate(&self, i: usize, j: usize, xi: f64, eta: f64) -> f64 {
        let t00 = self.temperatures[self.idx(i, j)];
        let t10 = self.temperatures[self.idx(i + 1, j)];
        let t01 = self.temperatures[self.idx(i, j + 1)];
        let t11 = self.temperatures[self.idx(i + 1, j + 1)];
        (1.0 - xi) * (1.0 - eta) * t00
            + xi * (1.0 - eta) * t10
            + (1.0 - xi) * eta * t01
            + xi * eta * t11
    }
    /// Apply a Dirichlet temperature boundary condition on the left face.
    pub fn set_left_bc(&mut self, t_wall: f64) {
        for j in 0..self.ny {
            let idx = self.idx(0, j);
            self.temperatures[idx] = t_wall;
        }
    }
    /// Apply a Dirichlet temperature boundary condition on the right face.
    pub fn set_right_bc(&mut self, t_wall: f64) {
        let nx = self.nx;
        for j in 0..self.ny {
            let idx = self.idx(nx - 1, j);
            self.temperatures[idx] = t_wall;
        }
    }
    /// Compute the average temperature.
    pub fn average(&self) -> f64 {
        self.temperatures.iter().sum::<f64>() / self.temperatures.len() as f64
    }
    /// Compute the maximum temperature.
    pub fn max_temp(&self) -> f64 {
        self.temperatures
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
    /// Compute the minimum temperature.
    pub fn min_temp(&self) -> f64 {
        self.temperatures
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min)
    }
}
/// Quadrature-free matching interpolation for conforming interface DOF transfer.
#[derive(Debug, Clone)]
pub struct ConformingTransfer {
    /// Interpolation weights: weights\[i\]\[j\] = weight from source DOF j to target DOF i
    pub weights: Vec<Vec<f64>>,
    /// Source positions
    pub src_positions: Vec<f64>,
    /// Target positions
    pub tgt_positions: Vec<f64>,
}
impl ConformingTransfer {
    /// Build linear interpolation transfer from source to target.
    pub fn linear_interpolation(src_pos: Vec<f64>, tgt_pos: Vec<f64>) -> Self {
        let n_src = src_pos.len();
        let n_tgt = tgt_pos.len();
        let mut weights = vec![vec![0.0; n_src]; n_tgt];
        for (i, &xt) in tgt_pos.iter().enumerate() {
            let j = src_pos.partition_point(|&xs| xs <= xt);
            let j = j.saturating_sub(1).min(n_src - 2);
            let x0 = src_pos[j];
            let x1 = src_pos[j + 1];
            let xi = if (x1 - x0).abs() > 1e-30 {
                (xt - x0) / (x1 - x0)
            } else {
                0.0
            };
            weights[i][j] = 1.0 - xi;
            weights[i][j + 1] = xi;
        }
        Self {
            weights,
            src_positions: src_pos,
            tgt_positions: tgt_pos,
        }
    }
    /// Transfer field values from source to target.
    pub fn transfer(&self, src_values: &[f64]) -> Vec<f64> {
        self.weights
            .iter()
            .map(|row| row.iter().zip(src_values.iter()).map(|(w, v)| w * v).sum())
            .collect()
    }
    /// Conservative transfer (preserves integral).
    pub fn conservative_transfer(&self, src_values: &[f64], src_areas: &[f64]) -> Vec<f64> {
        let raw = self.transfer(src_values);
        let src_total: f64 = src_values
            .iter()
            .zip(src_areas.iter())
            .map(|(v, a)| v * a)
            .sum();
        let tgt_total: f64 = raw.iter().sum::<f64>();
        if tgt_total.abs() > 1e-30 {
            raw.iter().map(|v| v * src_total / tgt_total).collect()
        } else {
            raw
        }
    }
}
/// ALE formulation fluid-structure interface.
///
/// Enforces kinematic compatibility: structural velocity = fluid velocity at interface.
#[derive(Debug, Clone)]
pub struct AleInterface {
    /// Interface node indices in the global mesh
    pub node_indices: Vec<usize>,
    /// Fluid-side normal (outward from fluid domain)
    pub normals: Vec<[f64; 2]>,
}
impl AleInterface {
    /// Create a new ALE interface.
    pub fn new(node_indices: Vec<usize>, normals: Vec<[f64; 2]>) -> Self {
        assert_eq!(node_indices.len(), normals.len());
        Self {
            node_indices,
            normals,
        }
    }
    /// Enforce no-penetration condition: set mesh velocity to match structural velocity.
    pub fn enforce_kinematic_bc(&self, nodes: &mut [AleNode], struct_vel: &[[f64; 2]]) {
        for (k, &idx) in self.node_indices.iter().enumerate() {
            nodes[idx].mesh_velocity = struct_vel[k];
            nodes[idx].fluid_velocity = struct_vel[k];
        }
    }
    /// Compute the traction force on the structure from fluid pressure.
    ///
    /// f = -p · n · area (per unit depth in 2D)
    pub fn compute_traction(&self, nodes: &[AleNode], segment_lengths: &[f64]) -> Vec<[f64; 2]> {
        self.node_indices
            .iter()
            .zip(self.normals.iter())
            .zip(segment_lengths.iter())
            .map(|((&idx, n), &ds)| {
                let p = nodes[idx].pressure;
                [-p * n[0] * ds, -p * n[1] * ds]
            })
            .collect()
    }
}
/// Thermal expansion state for a single integration point.
///
/// Stores the effective thermal strain vector in Voigt notation (3-D: εxx, εyy, εzz).
#[derive(Debug, Clone)]
pub struct ThermalExpansionState {
    /// Temperature at this point \[K\].
    pub temperature: f64,
    /// Thermal strain vector \[ε_xx, ε_yy, ε_zz, γ_xy, γ_yz, γ_zx\].
    pub thermal_strain: [f64; 6],
    /// Reference temperature used when the strain was computed \[K\].
    pub t_ref: f64,
    /// Isotropic thermal expansion coefficient \[1/K\].
    pub alpha: f64,
}
impl ThermalExpansionState {
    /// Create a new expansion state at the reference temperature.
    pub fn new(t_ref: f64, alpha: f64) -> Self {
        Self {
            temperature: t_ref,
            thermal_strain: [0.0; 6],
            t_ref,
            alpha,
        }
    }
    /// Update thermal strain from the current temperature.
    ///
    /// ε_th = α·(T − T_ref)·\[1,1,1,0,0,0\]
    pub fn update(&mut self, temperature: f64) {
        self.temperature = temperature;
        let eps = self.alpha * (temperature - self.t_ref);
        self.thermal_strain = [eps, eps, eps, 0.0, 0.0, 0.0];
    }
    /// Compute the elastic strain given total strain and thermal strain.
    ///
    /// ε_el = ε_total − ε_th
    pub fn elastic_strain(&self, total_strain: [f64; 6]) -> [f64; 6] {
        let mut el = [0.0f64; 6];
        for k in 0..6 {
            el[k] = total_strain[k] - self.thermal_strain[k];
        }
        el
    }
    /// Effective volumetric thermal strain.
    pub fn volumetric_thermal_strain(&self) -> f64 {
        self.thermal_strain[0] + self.thermal_strain[1] + self.thermal_strain[2]
    }
}
/// ALE mesh with topology for mesh-velocity update and quality control.
#[derive(Debug, Clone)]
pub struct AleMesh {
    /// Node positions in reference configuration.
    pub ref_positions: Vec<[f64; 2]>,
    /// Node positions in current spatial configuration.
    pub positions: Vec<[f64; 2]>,
    /// Mesh velocities \[m/s\].
    pub mesh_velocity: Vec<[f64; 2]>,
    /// Connectivity: pairs of node indices forming edges.
    pub edges: Vec<[usize; 2]>,
    /// Smoothing strategy.
    pub smoothing: MeshSmoothing,
}
impl AleMesh {
    /// Create a 1-D chain mesh from a list of reference positions.
    pub fn chain(ref_positions: Vec<[f64; 2]>) -> Self {
        let n = ref_positions.len();
        let edges = (0..n.saturating_sub(1)).map(|i| [i, i + 1]).collect();
        let positions = ref_positions.clone();
        let mesh_velocity = vec![[0.0; 2]; n];
        Self {
            ref_positions,
            positions,
            mesh_velocity,
            edges,
            smoothing: MeshSmoothing::Laplacian,
        }
    }
    /// Apply Laplacian smoothing to interior nodes for one sweep.
    ///
    /// Interior nodes are moved to the arithmetic mean of their neighbours.
    pub fn smooth_laplacian(&mut self) {
        let n = self.positions.len();
        if n < 3 {
            return;
        }
        let mut new_pos = self.positions.clone();
        let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &[a, b] in &self.edges {
            neighbours[a].push(b);
            neighbours[b].push(a);
        }
        for i in 1..n - 1 {
            let nb = &neighbours[i];
            if nb.is_empty() {
                continue;
            }
            let mut sum = [0.0f64; 2];
            for &j in nb {
                sum[0] += self.positions[j][0];
                sum[1] += self.positions[j][1];
            }
            let cnt = nb.len() as f64;
            new_pos[i] = [sum[0] / cnt, sum[1] / cnt];
        }
        self.positions = new_pos;
    }
    /// Compute mesh velocity from the change in position over `dt`.
    pub fn update_mesh_velocity(&mut self, old_positions: &[[f64; 2]], dt: f64) {
        let dt = dt.max(1e-30);
        for (i, mv) in self.mesh_velocity.iter_mut().enumerate() {
            if i < old_positions.len() {
                mv[0] = (self.positions[i][0] - old_positions[i][0]) / dt;
                mv[1] = (self.positions[i][1] - old_positions[i][1]) / dt;
            }
        }
    }
    /// Compute the minimum edge length (mesh quality metric).
    pub fn min_edge_length(&self) -> f64 {
        self.edges.iter().fold(f64::INFINITY, |acc, &[a, b]| {
            let dx = self.positions[a][0] - self.positions[b][0];
            let dy = self.positions[a][1] - self.positions[b][1];
            acc.min((dx * dx + dy * dy).sqrt())
        })
    }
    /// Compute a mesh distortion metric (ratio of current to reference edge lengths).
    ///
    /// Returns 1.0 for undistorted mesh, > 1.0 for stretched meshes.
    pub fn distortion_metric(&self) -> f64 {
        if self.edges.is_empty() {
            return 1.0;
        }
        let mut max_ratio = 0.0f64;
        for &[a, b] in &self.edges {
            let dx_ref = self.ref_positions[a][0] - self.ref_positions[b][0];
            let dy_ref = self.ref_positions[a][1] - self.ref_positions[b][1];
            let l_ref = (dx_ref * dx_ref + dy_ref * dy_ref).sqrt().max(1e-30);
            let dx_cur = self.positions[a][0] - self.positions[b][0];
            let dy_cur = self.positions[a][1] - self.positions[b][1];
            let l_cur = (dx_cur * dx_cur + dy_cur * dy_cur).sqrt();
            let ratio = (l_cur / l_ref - 1.0).abs();
            if ratio > max_ratio {
                max_ratio = ratio;
            }
        }
        1.0 + max_ratio
    }
}
/// Chemical species diffusion coupled with mechanical stress.
///
/// Concentration field c coupled to volumetric strain:
/// ∂c/∂t = ∇·(D·∇c) + D·Ω·c / (R·T) · ∇²σ_h
#[derive(Debug, Clone)]
pub struct ChemoMechanicalProblem {
    /// Concentration field \[mol/m³\]
    pub concentration: Vec<f64>,
    /// Stress field (hydrostatic) \[Pa\]
    pub hydrostatic_stress: Vec<f64>,
    /// Number of nodes
    pub n_nodes: usize,
    /// Element size h \[m\]
    pub h: f64,
    /// Diffusion coefficient D \[m²/s\]
    pub diff_coeff: f64,
    /// Partial molar volume Ω \[m³/mol\]
    pub partial_molar_vol: f64,
    /// Elastic modulus \[Pa\]
    pub young_mod: f64,
    /// Poisson ratio
    pub nu: f64,
    /// Temperature \[K\]
    pub temperature: f64,
}
impl ChemoMechanicalProblem {
    /// Create a 1D chemo-mechanical problem.
    pub fn new_1d(n_nodes: usize, length: f64, c0: f64) -> Self {
        let h = length / (n_nodes as f64 - 1.0);
        Self {
            concentration: vec![c0; n_nodes],
            hydrostatic_stress: vec![0.0; n_nodes],
            n_nodes,
            h,
            diff_coeff: 1e-14,
            partial_molar_vol: 3.49e-6,
            young_mod: 100e9,
            nu: 0.3,
            temperature: 298.15,
        }
    }
    /// Compute chemical diffusion potential: μ = μ₀ + R·T·ln(c) − Ω·σ_h
    pub fn chemical_potential(&self, node: usize) -> f64 {
        let c = self.concentration[node].max(1e-15);
        let sigma_h = self.hydrostatic_stress[node];
        GAS_CONSTANT * self.temperature * c.ln() - self.partial_molar_vol * sigma_h
    }
    /// Compute hydrostatic stress from concentration using isotropic expansion.
    ///
    /// σ_h = -E·Ω·(c − c0) / (3·(1 − 2·ν))
    pub fn update_stress_from_concentration(&mut self, c0: f64) {
        let factor = -self.young_mod * self.partial_molar_vol / (3.0 * (1.0 - 2.0 * self.nu));
        for i in 0..self.n_nodes {
            self.hydrostatic_stress[i] = factor * (self.concentration[i] - c0);
        }
    }
    /// Explicit time step for the diffusion equation.
    pub fn diffusion_step(&mut self, dt: f64, c_left: f64, c_right: f64) {
        let fo = self.diff_coeff * dt / (self.h * self.h);
        assert!(fo <= 0.5, "Fourier number {fo:.3} > 0.5: reduce dt");
        let mut new_c = self.concentration.clone();
        new_c[0] = c_left;
        new_c[self.n_nodes - 1] = c_right;
        for (i, new_c_i) in new_c.iter_mut().enumerate().skip(1).take(self.n_nodes - 2) {
            *new_c_i = self.concentration[i]
                + fo * (self.concentration[i + 1] - 2.0 * self.concentration[i]
                    + self.concentration[i - 1]);
        }
        self.concentration = new_c;
    }
    /// Staggered chemo-mechanical step.
    pub fn coupled_step(&mut self, dt: f64, c_left: f64, c_right: f64, c0: f64) {
        self.diffusion_step(dt, c_left, c_right);
        self.update_stress_from_concentration(c0);
    }
    /// Total amount of species in the domain.
    pub fn total_species(&self) -> f64 {
        self.concentration.iter().sum::<f64>() * self.h
    }
}
/// Coupling strategy enumeration.
#[derive(Debug, Clone, PartialEq)]
pub enum CouplingStrategy {
    /// Staggered (partitioned, explicit transfer)
    Staggered,
    /// Iterative partitioned with relaxation
    IterativePartitioned {
        /// Maximum number of sub-iterations per time step.
        max_iter: usize,
        /// Convergence tolerance for residual.
        tol: f64,
    },
    /// Monolithic (fully coupled single system matrix)
    Monolithic,
    /// Semi-monolithic (block-triangular approximation)
    SemiMonolithic,
}
/// Global residual vector from a multi-physics coupled problem.
///
/// Assembles contributions from structural, thermal, and other field residuals.
#[derive(Debug, Clone)]
pub struct MultiPhysicsResidual {
    /// Structural residual R_u = K_uu · u − f_u.
    pub r_struct: Vec<f64>,
    /// Thermal residual R_T = K_TT · T − f_T.
    pub r_thermal: Vec<f64>,
    /// Optional third-field residual (e.g., moisture, electric potential).
    pub r_aux: Vec<f64>,
}
impl MultiPhysicsResidual {
    /// Create a zero residual for given DOF counts.
    pub fn new(n_struct: usize, n_thermal: usize, n_aux: usize) -> Self {
        Self {
            r_struct: vec![0.0; n_struct],
            r_thermal: vec![0.0; n_thermal],
            r_aux: vec![0.0; n_aux],
        }
    }
    /// L2 norm of the structural residual.
    pub fn struct_norm(&self) -> f64 {
        self.r_struct.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
    /// L2 norm of the thermal residual.
    pub fn thermal_norm(&self) -> f64 {
        self.r_thermal.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
    /// L2 norm of the auxiliary residual.
    pub fn aux_norm(&self) -> f64 {
        self.r_aux.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
    /// Combined (stacked) residual norm.
    pub fn total_norm(&self) -> f64 {
        let s = self.struct_norm().powi(2) + self.thermal_norm().powi(2) + self.aux_norm().powi(2);
        s.sqrt()
    }
    /// Assemble structural residual: R_u = K · u − f.
    pub fn assemble_structural(&mut self, k: &[Vec<f64>], u: &[f64], f: &[f64]) {
        let n = self.r_struct.len().min(u.len()).min(f.len());
        for i in 0..n {
            let ku: f64 = (0..k[i].len().min(u.len())).map(|j| k[i][j] * u[j]).sum();
            self.r_struct[i] = ku - f[i];
        }
    }
    /// Assemble thermal residual: R_T = K_T · T − f_T.
    pub fn assemble_thermal(&mut self, k_t: &[Vec<f64>], temps: &[f64], f_t: &[f64]) {
        let n = self.r_thermal.len().min(temps.len()).min(f_t.len());
        for i in 0..n {
            let kt: f64 = (0..k_t[i].len().min(temps.len()))
                .map(|j| k_t[i][j] * temps[j])
                .sum();
            self.r_thermal[i] = kt - f_t[i];
        }
    }
}
/// Proper Orthogonal Decomposition (POD) reduced-order model for coupling.
///
/// Projects the full-order system onto a low-dimensional subspace.
#[derive(Debug, Clone)]
pub struct PodRom {
    /// POD basis vectors (columns), shape: n_full × n_modes
    pub basis: Vec<Vec<f64>>,
    /// Reduced-order stiffness matrix (n_modes × n_modes)
    pub k_rom: Vec<Vec<f64>>,
    /// Reduced-order mass matrix (n_modes × n_modes)
    pub m_rom: Vec<Vec<f64>>,
    /// Reduced-order force vector
    pub f_rom: Vec<f64>,
    /// Reduced coordinates
    pub q: Vec<f64>,
    /// Full-order DOFs
    pub n_full: usize,
    /// Number of POD modes retained
    pub n_modes: usize,
}
impl PodRom {
    /// Create a POD-ROM from snapshot data.
    ///
    /// Snapshots: n_full × n_snapshots matrix.
    pub fn from_snapshots(snapshots: &[Vec<f64>], n_modes: usize) -> Self {
        let n_snap = snapshots.len();
        let n_full = if n_snap > 0 { snapshots[0].len() } else { 0 };
        let n_modes = n_modes.min(n_snap).min(n_full);
        let mut basis = vec![vec![0.0; n_full]; n_modes];
        for (m, snap) in snapshots.iter().take(n_modes).enumerate() {
            let norm: f64 = snap.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
            for i in 0..n_full {
                basis[m][i] = snap[i] / norm;
            }
        }
        for m in 1..n_modes {
            for k in 0..m {
                let dot: f64 = (0..n_full).map(|i| basis[m][i] * basis[k][i]).sum();
                let bk = basis[k].clone();
                for i in 0..n_full {
                    basis[m][i] -= dot * bk[i];
                }
            }
            let norm: f64 = basis[m]
                .iter()
                .map(|x| x * x)
                .sum::<f64>()
                .sqrt()
                .max(1e-30);
            for val in basis[m].iter_mut() {
                *val /= norm;
            }
        }
        let k_rom = vec![vec![0.0; n_modes]; n_modes];
        let m_rom = vec![vec![0.0; n_modes]; n_modes];
        let f_rom = vec![0.0; n_modes];
        let q = vec![0.0; n_modes];
        Self {
            basis,
            k_rom,
            m_rom,
            f_rom,
            q,
            n_full,
            n_modes,
        }
    }
    /// Project a full-order vector onto the reduced basis.
    pub fn project(&self, v_full: &[f64]) -> Vec<f64> {
        (0..self.n_modes)
            .map(|m| (0..self.n_full).map(|i| self.basis[m][i] * v_full[i]).sum())
            .collect()
    }
    /// Reconstruct the full-order vector from reduced coordinates.
    pub fn reconstruct(&self, q: &[f64]) -> Vec<f64> {
        let mut v = vec![0.0; self.n_full];
        for (m, &q_m) in q.iter().enumerate().take(self.n_modes) {
            for (i, v_i) in v.iter_mut().enumerate().take(self.n_full) {
                *v_i += self.basis[m][i] * q_m;
            }
        }
        v
    }
    /// Project a full-order stiffness matrix into the ROM space.
    ///
    /// K_r = Φᵀ K Φ
    pub fn project_matrix(&self, k_full: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = self.n_modes;
        let mut k_r = vec![vec![0.0; n]; n];
        for (i, k_r_row) in k_r.iter_mut().enumerate() {
            for (j, k_r_ij) in k_r_row.iter_mut().enumerate() {
                let mut val = 0.0f64;
                for (row, k_full_row) in k_full.iter().enumerate().take(self.n_full) {
                    let kphi_j: f64 = (0..self.n_full)
                        .map(|col| k_full_row[col] * self.basis[j][col])
                        .sum();
                    val += self.basis[i][row] * kphi_j;
                }
                *k_r_ij = val;
            }
        }
        k_r
    }
    /// Projection error ‖v − Φ Φᵀ v‖ / ‖v‖.
    pub fn projection_error(&self, v_full: &[f64]) -> f64 {
        let q = self.project(v_full);
        let v_approx = self.reconstruct(&q);
        let err: f64 = v_full
            .iter()
            .zip(v_approx.iter())
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt();
        let norm: f64 = v_full.iter().map(|x| x * x).sum::<f64>().sqrt().max(1e-30);
        err / norm
    }
}
/// Hygro-mechanical material: moisture-induced swelling coupled with diffusion.
///
/// Moisture strain: ε_hm = β_h · (c − c_ref)  \[analogous to thermal strain\]
#[derive(Debug, Clone)]
pub struct HygroMechanicalMaterial {
    /// Elastic modulus \[Pa\].
    pub e_mod: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Moisture expansion coefficient \[1/(kg/m³) or 1/wt%\].
    pub beta_h: f64,
    /// Reference moisture concentration c_ref \[kg/m³\].
    pub c_ref: f64,
    /// Moisture diffusivity D_m \[m²/s\].
    pub diffusivity: f64,
}
impl HygroMechanicalMaterial {
    /// Create a new hygro-mechanical material.
    pub fn new(e_mod: f64, nu: f64, beta_h: f64, c_ref: f64, diffusivity: f64) -> Self {
        Self {
            e_mod,
            nu,
            beta_h,
            c_ref,
            diffusivity,
        }
    }
    /// Moisture-induced volumetric strain at concentration c.
    pub fn moisture_strain(&self, concentration: f64) -> f64 {
        self.beta_h * (concentration - self.c_ref)
    }
    /// Hygrostress from total strain and moisture concentration.
    ///
    /// σ = E/(1−2ν) · (ε_vol − ε_hm)  (simplified volumetric form)
    pub fn hydrostatic_stress_from_concentration(
        &self,
        vol_strain: f64,
        concentration: f64,
    ) -> f64 {
        let eps_hm = self.moisture_strain(concentration);
        let k_bulk = self.e_mod / (3.0 * (1.0 - 2.0 * self.nu).max(1e-10));
        k_bulk * (vol_strain - eps_hm)
    }
}
/// Temperature-dependent elastic modulus model.
///
/// E(T) = E₀ · (1 − β · (T − T_ref))
#[derive(Debug, Clone)]
pub struct ThermoElasticMaterial {
    /// Reference Young's modulus \[Pa\] at T_ref
    pub e0: f64,
    /// Thermal degradation coefficient \[1/K\]
    pub beta: f64,
    /// Reference temperature \[K\]
    pub t_ref: f64,
    /// Poisson's ratio (temperature independent)
    pub nu: f64,
    /// Thermal expansion coefficient \[1/K\]
    pub alpha: f64,
    /// Thermal conductivity \[W/(m·K)\]
    pub kappa: f64,
    /// Specific heat capacity \[J/(kg·K)\]
    pub cp: f64,
    /// Density \[kg/m³\]
    pub rho: f64,
}
impl ThermoElasticMaterial {
    /// Create a new thermo-elastic material.
    pub fn new(e0: f64, nu: f64, alpha: f64, kappa: f64, cp: f64, rho: f64) -> Self {
        Self {
            e0,
            beta: 0.0,
            t_ref: 293.15,
            nu,
            alpha,
            kappa,
            cp,
            rho,
        }
    }
    /// Evaluate Young's modulus at temperature T.
    pub fn modulus_at(&self, temp: f64) -> f64 {
        self.e0 * (1.0 - self.beta * (temp - self.t_ref)).max(0.01)
    }
    /// Evaluate the plane-stress constitutive matrix D(T) \[3×3 Voigt\].
    pub fn constitutive_2d(&self, temp: f64) -> [[f64; 3]; 3] {
        let e = self.modulus_at(temp);
        let nu = self.nu;
        let factor = e / (1.0 - nu * nu);
        [
            [factor, factor * nu, 0.0],
            [factor * nu, factor, 0.0],
            [0.0, 0.0, factor * (1.0 - nu) * 0.5],
        ]
    }
    /// Thermal load vector contribution for a single element.
    ///
    /// f_th = D · α · ΔT · {1, 1, 0}ᵀ  (plane stress, 2D)
    pub fn thermal_load_vector_2d(&self, temp: f64, delta_t: f64) -> [f64; 3] {
        let d = self.constitutive_2d(temp);
        let eps_th = self.alpha * delta_t;
        let strain_th = [eps_th, eps_th, 0.0];
        let mut f = [0.0f64; 3];
        for i in 0..3 {
            for j in 0..3 {
                f[i] += d[i][j] * strain_th[j];
            }
        }
        f
    }
}
/// Moving Least Squares (MLS) interpolation for data transfer between
/// non-conforming interface meshes.
///
/// Given source data at irregularly spaced points, evaluates the MLS
/// approximant at arbitrary target points.
#[derive(Debug, Clone)]
pub struct MlsInterpolator {
    /// Source point positions \[m\].
    pub source_positions: Vec<f64>,
    /// Weight function radius (support radius) \[m\].
    pub support_radius: f64,
}
impl MlsInterpolator {
    /// Create a new MLS interpolator.
    pub fn new(source_positions: Vec<f64>, support_radius: f64) -> Self {
        Self {
            source_positions,
            support_radius,
        }
    }
    /// Gaussian weight function: w(r) = exp(−(r/h)²).
    fn weight(&self, x_src: f64, x_tgt: f64) -> f64 {
        let r = (x_src - x_tgt).abs();
        let h = self.support_radius;
        if r > h { 0.0 } else { (-(r / h).powi(2)).exp() }
    }
    /// MLS interpolation using linear basis {1, x}.
    ///
    /// Returns the interpolated value at `x_target`.
    pub fn interpolate(&self, source_values: &[f64], x_target: f64) -> f64 {
        let mut a00 = 0.0f64;
        let mut a01 = 0.0f64;
        let mut a11 = 0.0f64;
        let mut b0 = 0.0f64;
        let mut b1 = 0.0f64;
        for (j, &xs) in self.source_positions.iter().enumerate() {
            let fj = if j < source_values.len() {
                source_values[j]
            } else {
                0.0
            };
            let w = self.weight(xs, x_target);
            if w < 1e-30 {
                continue;
            }
            a00 += w;
            a01 += w * xs;
            a11 += w * xs * xs;
            b0 += w * fj;
            b1 += w * xs * fj;
        }
        let det = a00 * a11 - a01 * a01;
        if det.abs() < 1e-30 {
            return self
                .source_positions
                .iter()
                .enumerate()
                .min_by(|&(_, a), &(_, b)| {
                    (a - x_target)
                        .abs()
                        .partial_cmp(&(b - x_target).abs())
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .and_then(|(i, _)| source_values.get(i).copied())
                .unwrap_or(0.0);
        }
        let c0 = (a11 * b0 - a01 * b1) / det;
        let c1 = (a00 * b1 - a01 * b0) / det;
        c0 + c1 * x_target
    }
    /// Transfer all source values to a set of target positions.
    pub fn transfer(&self, source_values: &[f64], target_positions: &[f64]) -> Vec<f64> {
        target_positions
            .iter()
            .map(|&xt| self.interpolate(source_values, xt))
            .collect()
    }
    /// Compute the maximum interpolation error at source points (leave-one-out).
    pub fn self_test_error(&self, source_values: &[f64]) -> f64 {
        self.source_positions
            .iter()
            .enumerate()
            .map(|(i, &xs)| {
                let mut pos_excl = self.source_positions.clone();
                pos_excl.remove(i);
                let mut val_excl: Vec<f64> = source_values.to_vec();
                if i < val_excl.len() {
                    val_excl.remove(i);
                }
                let interp =
                    MlsInterpolator::new(pos_excl, self.support_radius).interpolate(&val_excl, xs);
                let fi = if i < source_values.len() {
                    source_values[i]
                } else {
                    0.0
                };
                (interp - fi).abs()
            })
            .fold(0.0_f64, f64::max)
    }
}
/// Mesh-smoothing strategy for ALE moving meshes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshSmoothing {
    /// Laplacian (Jacobi) smoothing.
    Laplacian,
    /// Algebraic smoothing (spring analogy).
    SpringAnalogy,
    /// No smoothing (rigid body mesh motion only).
    Rigid,
}
/// Acoustic fluid element (rectangular cavity, 1D model).
#[derive(Debug, Clone)]
pub struct AcousticCavity1D {
    /// Acoustic pressure at each node \[Pa\]
    pub pressure: Vec<f64>,
    /// Particle velocity at each node \[m/s\]
    pub velocity: Vec<f64>,
    /// Speed of sound \[m/s\]
    pub c_sound: f64,
    /// Fluid density \[kg/m³\]
    pub rho_f: f64,
    /// Element length \[m\]
    pub h: f64,
    /// Number of nodes
    pub n_nodes: usize,
}
impl AcousticCavity1D {
    /// Create a 1D acoustic cavity.
    pub fn new(n_nodes: usize, length: f64, c_sound: f64, rho_f: f64) -> Self {
        let h = length / (n_nodes as f64 - 1.0);
        Self {
            pressure: vec![0.0; n_nodes],
            velocity: vec![0.0; n_nodes],
            c_sound,
            rho_f,
            h,
            n_nodes,
        }
    }
    /// Acoustic impedance Z = ρ·c.
    #[inline]
    pub fn impedance(&self) -> f64 {
        self.rho_f * self.c_sound
    }
    /// Apply a harmonic pressure excitation at the left boundary.
    pub fn apply_pressure_source(&mut self, p_amp: f64, omega: f64, t: f64) {
        self.pressure[0] = p_amp * (omega * t).sin();
    }
    /// Apply a wall boundary condition (zero velocity) at the right boundary.
    pub fn apply_wall_bc(&mut self) {
        self.velocity[self.n_nodes - 1] = 0.0;
    }
    /// One-step Lax-Wendroff update (linearized acoustic equations).
    pub fn step(&mut self, dt: f64) {
        let r = self.c_sound * dt / self.h;
        assert!(r <= 1.0, "CFL {r:.3} > 1: reduce dt");
        let mut new_p = self.pressure.clone();
        let mut new_v = self.velocity.clone();
        for i in 1..self.n_nodes - 1 {
            new_p[i] = self.pressure[i]
                - 0.5
                    * r
                    * self.rho_f
                    * self.c_sound
                    * (self.velocity[i + 1] - self.velocity[i - 1]);
            new_v[i] = self.velocity[i]
                - 0.5 * r / (self.rho_f * self.c_sound)
                    * (self.pressure[i + 1] - self.pressure[i - 1]);
        }
        self.pressure = new_p;
        self.velocity = new_v;
    }
    /// Structural response force on the left boundary from acoustic pressure.
    ///
    /// F = p · A (here A = 1.0 m² unit area)
    pub fn boundary_force(&self) -> f64 {
        self.pressure[0]
    }
    /// Compute total acoustic energy density.
    pub fn total_energy(&self) -> f64 {
        let mut e = 0.0;
        for i in 0..self.n_nodes {
            let e_p = self.pressure[i] * self.pressure[i]
                / (2.0 * self.rho_f * self.c_sound * self.c_sound);
            let e_v = 0.5 * self.rho_f * self.velocity[i] * self.velocity[i];
            e += (e_p + e_v) * self.h;
        }
        e
    }
}
/// Magnetostrictive material model (Terfenol-D like).
#[derive(Debug, Clone)]
pub struct MagnetostrictiveMaterial {
    /// Saturation magnetostiction strain λₛ \[dimensionless\]
    pub lambda_s: f64,
    /// Saturation magnetization Mₛ \[A/m\]
    pub m_sat: f64,
    /// Elastic modulus at zero field \[Pa\]
    pub e_h0: f64,
    /// Delta-E effect coefficient
    pub delta_e: f64,
    /// Piezomagnetic coefficient d₃₃ \[m/A\]
    pub d33: f64,
}
impl MagnetostrictiveMaterial {
    /// Terfenol-D approximate properties.
    pub fn terfenol_d() -> Self {
        Self {
            lambda_s: 1600e-6,
            m_sat: 7.65e5,
            e_h0: 30e9,
            delta_e: 0.25,
            d33: 7e-9,
        }
    }
    /// Magnetostrictive strain at applied field H.
    ///
    /// ε_ms = λₛ · (M(H)/Mₛ)² using Langevin model approximation.
    pub fn magnetostrictive_strain(&self, h_field: f64) -> f64 {
        let m_ratio = self.langevin_magnetization(h_field) / self.m_sat;
        self.lambda_s * m_ratio * m_ratio
    }
    /// Langevin magnetization model: M = Mₛ · (coth(x) - 1/x), x = μ₀·H·m / (kB·T).
    ///
    /// Simplified version: M/Mₛ = tanh(H / H_c) where H_c is the coercive-like field.
    pub fn langevin_magnetization(&self, h_field: f64) -> f64 {
        let h_c = self.m_sat / 10.0;
        self.m_sat * (h_field / h_c).tanh()
    }
    /// Effective elastic modulus including delta-E effect.
    pub fn effective_modulus(&self, h_field: f64) -> f64 {
        let m_ratio = self.langevin_magnetization(h_field).abs() / self.m_sat;
        self.e_h0 * (1.0 + self.delta_e * m_ratio * m_ratio)
    }
    /// Compute total stress: σ = E_eff · (ε − ε_ms)
    pub fn stress(&self, strain: f64, h_field: f64) -> f64 {
        let eps_ms = self.magnetostrictive_strain(h_field);
        let e_eff = self.effective_modulus(h_field);
        e_eff * (strain - eps_ms)
    }
}
/// Partitioned vs monolithic coupling comparison result.
#[derive(Debug, Clone)]
pub struct CouplingComparison {
    /// Number of iterations for staggered (partitioned) scheme.
    pub staggered_iters: usize,
    /// Whether partitioned scheme converged.
    pub staggered_converged: bool,
    /// Residual for monolithic (direct) solve.
    pub monolithic_residual: f64,
}
impl CouplingComparison {
    /// Compare staggered vs monolithic coupling for a 2×2 block system.
    ///
    /// The monolithic solve uses exact Gaussian elimination.
    pub fn compare_2x2(
        k_aa: f64,
        k_ab: f64,
        k_ba: f64,
        k_bb: f64,
        f_a: f64,
        f_b: f64,
        max_iter: usize,
        tol: f64,
    ) -> Self {
        let det = k_aa * k_bb - k_ab * k_ba;
        let monolithic_residual = if det.abs() > 1e-30 {
            let _ua = (f_a * k_bb - k_ab * f_b) / det;
            let _ub = (k_aa * f_b - f_a * k_ba) / det;
            0.0
        } else {
            f64::INFINITY
        };
        let mut ua = 0.0f64;
        let mut ub = 0.0f64;
        let mut scheme = StaggeredCouplingScheme::new(max_iter, tol);
        let mut iters = 0;
        for _ in 0..max_iter {
            let ua_new = if k_aa.abs() > 1e-30 {
                (f_a - k_ab * ub) / k_aa
            } else {
                0.0
            };
            let ub_new = if k_bb.abs() > 1e-30 {
                (f_b - k_ba * ua_new) / k_bb
            } else {
                0.0
            };
            let (conv, _rel) = scheme.check_convergence(&[ua, ub], &[ua_new, ub_new]);
            ua = ua_new;
            ub = ub_new;
            iters += 1;
            if conv {
                break;
            }
        }
        CouplingComparison {
            staggered_iters: iters,
            staggered_converged: scheme.is_converged(),
            monolithic_residual,
        }
    }
}
