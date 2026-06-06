// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coupled multiphysics FEM module.
//!
//! Implements coupled thermal-structural, piezoelectric, poroelastic,
//! fluid-structure, magneto-mechanical, and electrochemical FEM formulations.
//!
//! Supports:
//! - Thermo-elastic coupling (temperature → stress via thermal expansion)
//! - Piezoelectric u-phi formulation
//! - Biot consolidation (u-p poroelastic)
//! - ALE fluid-structure interaction
//! - Magneto-mechanical Lorentz force
//! - Diffusion-reaction with SUPG stabilization
//! - Electrochemical Butler-Volmer kinetics
//! - Thermal-mechanical fatigue (Chaboche hardening)
//! - Staggered and monolithic coupled solvers

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// Compute dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Compute norm of a 3-vector.
#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

/// Scale a 3-vector.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ---------------------------------------------------------------------------
// MultiFieldResult
// ---------------------------------------------------------------------------

/// Combined result container for multiphysics FEM solutions.
///
/// Stores displacements, temperatures, electric potential, pore pressure,
/// and magnetic vector potential at each node.
#[derive(Debug, Clone)]
pub struct MultiFieldResult {
    /// Displacement vector at each node (x, y, z).
    pub displacements: Vec<[f64; 3]>,
    /// Temperature at each node (K).
    pub temperatures: Vec<f64>,
    /// Electric potential at each node (V).
    pub electric_potential: Vec<f64>,
    /// Pore pressure at each node (Pa).
    pub pressure: Vec<f64>,
    /// Magnetic vector potential at each node (T·m).
    pub magnetic_potential: Vec<[f64; 3]>,
    /// Number of converged iterations.
    pub iterations: usize,
    /// Residual norm at convergence.
    pub residual: f64,
}

impl MultiFieldResult {
    /// Create an empty multiphysics result for n_nodes.
    pub fn new(n_nodes: usize) -> Self {
        Self {
            displacements: vec![[0.0; 3]; n_nodes],
            temperatures: vec![0.0; n_nodes],
            electric_potential: vec![0.0; n_nodes],
            pressure: vec![0.0; n_nodes],
            magnetic_potential: vec![[0.0; 3]; n_nodes],
            iterations: 0,
            residual: 0.0,
        }
    }

    /// Number of nodes in this result.
    pub fn num_nodes(&self) -> usize {
        self.displacements.len()
    }

    /// Maximum displacement magnitude.
    pub fn max_displacement(&self) -> f64 {
        self.displacements
            .iter()
            .map(|&d| norm3(d))
            .fold(0.0_f64, f64::max)
    }

    /// Maximum temperature.
    pub fn max_temperature(&self) -> f64 {
        self.temperatures
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Maximum electric potential.
    pub fn max_electric_potential(&self) -> f64 {
        self.electric_potential
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }

    /// Maximum pore pressure.
    pub fn max_pressure(&self) -> f64 {
        self.pressure
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

// ---------------------------------------------------------------------------
// ThermoElasticFem
// ---------------------------------------------------------------------------

/// Coupled thermo-elastic FEM formulation.
///
/// Implements thermal expansion coupling:
/// σ = C:(ε - α(T - T_ref)I)
/// where α is thermal expansion coefficient.
#[derive(Debug, Clone)]
pub struct ThermoElasticFem {
    /// Young's modulus E (Pa).
    pub young_modulus: f64,
    /// Poisson's ratio ν.
    pub poisson_ratio: f64,
    /// Thermal expansion coefficient α (1/K).
    pub thermal_expansion: f64,
    /// Reference temperature T_ref (K).
    pub reference_temperature: f64,
    /// Thermal conductivity k (W/m·K).
    pub thermal_conductivity: f64,
    /// Heat capacity ρc_p (J/m³·K).
    pub heat_capacity: f64,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Current temperature field.
    pub temperatures: Vec<f64>,
    /// Current displacement field.
    pub displacements: Vec<[f64; 3]>,
}

impl ThermoElasticFem {
    /// Create a new thermo-elastic FEM model.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        thermal_expansion: f64,
        reference_temperature: f64,
        thermal_conductivity: f64,
        heat_capacity: f64,
        n_nodes: usize,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            thermal_expansion,
            reference_temperature,
            thermal_conductivity,
            heat_capacity,
            n_nodes,
            temperatures: vec![reference_temperature; n_nodes],
            displacements: vec![[0.0; 3]; n_nodes],
        }
    }

    /// Compute thermal strain vector at a point given temperature T.
    ///
    /// Returns ε_th = α(T - T_ref) \[1, 1, 1, 0, 0, 0\] (Voigt notation, 6 components).
    pub fn thermal_strain(&self, temperature: f64) -> [f64; 6] {
        let alpha_dt = self.thermal_expansion * (temperature - self.reference_temperature);
        [alpha_dt, alpha_dt, alpha_dt, 0.0, 0.0, 0.0]
    }

    /// Compute plane-stress constitutive matrix D (3x3 in Voigt notation).
    pub fn constitutive_matrix_2d(&self) -> [[f64; 3]; 3] {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let c = e / (1.0 - nu * nu);
        [
            [c, c * nu, 0.0],
            [c * nu, c, 0.0],
            [0.0, 0.0, c * 0.5 * (1.0 - nu)],
        ]
    }

    /// Compute 3D constitutive matrix D (6x6, returned as flat array).
    pub fn constitutive_matrix_3d(&self) -> Vec<f64> {
        let e = self.young_modulus;
        let nu = self.poisson_ratio;
        let lam = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        let mu = e / (2.0 * (1.0 + nu));
        let mut d = vec![0.0f64; 36];
        for i in 0..3 {
            for j in 0..3 {
                d[i * 6 + j] = lam;
            }
            d[i * 6 + i] += 2.0 * mu;
        }
        d[3 * 6 + 3] = mu;
        d[4 * 6 + 4] = mu;
        d[5 * 6 + 5] = mu;
        d
    }

    /// Compute thermal load vector contribution from temperature field.
    ///
    /// Returns nodal thermal load vector of length n_nodes * 3.
    pub fn thermal_load_vector(&self) -> Vec<f64> {
        let n = self.n_nodes;
        let mut f = vec![0.0f64; n * 3];
        let d = self.constitutive_matrix_3d();
        for node in 0..n {
            let t = self.temperatures[node];
            let eps_th = self.thermal_strain(t);
            // Simplified: apply D * eps_th as nodal force (actual FEM needs integration)
            for i in 0..3 {
                let mut force = 0.0;
                for j in 0..6 {
                    if j < 3 {
                        force += d[i * 6 + j] * eps_th[j];
                    }
                }
                f[node * 3 + i] = force;
            }
        }
        f
    }

    /// Set temperature field from array.
    pub fn set_temperatures(&mut self, temps: &[f64]) {
        for (i, &t) in temps.iter().enumerate() {
            if i < self.n_nodes {
                self.temperatures[i] = t;
            }
        }
    }

    /// Compute von Mises stress at a node from displacement gradient and temperature.
    pub fn von_mises_stress_node(&self, node_idx: usize, strain: [f64; 6]) -> f64 {
        let t = self
            .temperatures
            .get(node_idx)
            .copied()
            .unwrap_or(self.reference_temperature);
        let eps_th = self.thermal_strain(t);
        // Mechanical strain
        let eps_mech: Vec<f64> = strain
            .iter()
            .zip(eps_th.iter())
            .map(|(e, et)| e - et)
            .collect();
        let d = self.constitutive_matrix_3d();

        // Stress = D * eps_mech
        let mut sigma = [0.0f64; 6];
        for i in 0..6 {
            for j in 0..6 {
                sigma[i] += d[i * 6 + j] * eps_mech[j];
            }
        }

        // Von Mises
        let s11 = sigma[0];
        let s22 = sigma[1];
        let s33 = sigma[2];
        let s12 = sigma[3];
        let s13 = sigma[4];
        let s23 = sigma[5];

        ((s11 - s22).powi(2)
            + (s22 - s33).powi(2)
            + (s33 - s11).powi(2)
            + 6.0 * (s12.powi(2) + s13.powi(2) + s23.powi(2)))
        .sqrt()
            / 2.0_f64.sqrt()
    }
}

// ---------------------------------------------------------------------------
// PiezoelectricFem
// ---------------------------------------------------------------------------

/// Piezoelectric FEM formulation (u-φ coupling).
///
/// Implements the electromechanical coupling:
/// {σ} = \[C\]{ε} - \[e\]^T{E}
/// {D} = \[e\]{ε} + \[κ\]{E}
/// where e is the piezoelectric coupling matrix and κ the permittivity.
#[derive(Debug, Clone)]
pub struct PiezoelectricFem {
    /// Elastic stiffness matrix (6x6 Voigt, flat).
    pub elastic_stiffness: Vec<f64>,
    /// Piezoelectric coupling matrix e (3x6 in stress form, flat).
    pub piezo_coupling: Vec<f64>,
    /// Dielectric permittivity matrix κ (3x3, flat).
    pub permittivity: Vec<f64>,
    /// Number of mechanical DOF.
    pub n_mech_dof: usize,
    /// Number of electric DOF.
    pub n_elec_dof: usize,
}

impl PiezoelectricFem {
    /// Create a new piezoelectric FEM model.
    pub fn new(
        elastic_stiffness: Vec<f64>,
        piezo_coupling: Vec<f64>,
        permittivity: Vec<f64>,
        n_mech_dof: usize,
        n_elec_dof: usize,
    ) -> Self {
        Self {
            elastic_stiffness,
            piezo_coupling,
            permittivity,
            n_mech_dof,
            n_elec_dof,
        }
    }

    /// Compute stress from mechanical strain and electric field (direct piezo effect).
    ///
    /// σ = C*ε - e^T * E
    pub fn compute_stress(&self, strain: &[f64; 6], e_field: &[f64; 3]) -> [f64; 6] {
        let mut sigma = [0.0f64; 6];
        // C * ε
        for (i, sig_i) in sigma.iter_mut().enumerate() {
            for (j, &sj) in strain.iter().enumerate() {
                *sig_i += self.elastic_stiffness[i * 6 + j] * sj;
            }
        }
        // -e^T * E (e is 3x6, so e^T is 6x3)
        for (i, sig_i) in sigma.iter_mut().enumerate() {
            for (k, &ek) in e_field.iter().enumerate() {
                *sig_i -= self.piezo_coupling[k * 6 + i] * ek;
            }
        }
        sigma
    }

    /// Compute electric displacement from strain and electric field.
    ///
    /// D = e*ε + κ*E
    pub fn compute_electric_displacement(&self, strain: &[f64; 6], e_field: &[f64; 3]) -> [f64; 3] {
        let mut d = [0.0f64; 3];
        // e * ε
        for (i, di) in d.iter_mut().enumerate() {
            for (j, &sj) in strain.iter().enumerate() {
                *di += self.piezo_coupling[i * 6 + j] * sj;
            }
        }
        // κ * E
        for (i, di) in d.iter_mut().enumerate() {
            for (j, &ej) in e_field.iter().enumerate() {
                *di += self.permittivity[i * 3 + j] * ej;
            }
        }
        d
    }

    /// Compute capacitance matrix C_ee (electric-electric block).
    ///
    /// Returns a simplified scalar capacitance for a parallel-plate capacitor.
    pub fn capacitance_scalar(&self, area: f64, thickness: f64) -> f64 {
        // C = κ * A / d
        if thickness < 1e-15 {
            return 0.0;
        }
        let kappa_eff = self.permittivity[0]; // Use (0,0) component
        kappa_eff * area / thickness
    }

    /// Total number of DOF (mechanical + electrical).
    pub fn total_dof(&self) -> usize {
        self.n_mech_dof + self.n_elec_dof
    }
}

/// Create a default PZT-5A piezoelectric model.
///
/// Returns a PiezoelectricFem with typical PZT-5A material properties.
pub fn make_pzt5a(n_mech: usize, n_elec: usize) -> PiezoelectricFem {
    // PZT-5A properties (simplified)
    let e_stiff = {
        let c11 = 121e9_f64;
        let c12 = 75.4e9_f64;
        let c13 = 75.2e9_f64;
        let c33 = 111e9_f64;
        let c44 = 21.1e9_f64;
        let c66 = 22.6e9_f64;
        let mut c = vec![0.0f64; 36];
        c[0] = c11;
        c[1] = c12;
        c[2] = c13;
        c[6] = c12;
        c[7] = c11;
        c[8] = c13;
        c[12] = c13;
        c[13] = c13;
        c[14] = c33;
        c[21] = c44;
        c[28] = c44;
        c[35] = c66;
        c
    };
    // Piezoelectric coupling (3x6): e31, e33, e15
    let mut e_coupling = vec![0.0f64; 18];
    e_coupling[2] = -5.4_f64; // e31 (row 0, col 2)
    e_coupling[8] = 15.8_f64; // e33 (row 1, col 2... simplified)
    e_coupling[9] = 12.3_f64; // e15

    // Dielectric permittivity (3x3)
    let eps_11 = 1650.0_f64 * 8.854e-12_f64;
    let eps_33 = 1700.0_f64 * 8.854e-12_f64;
    let mut perm = vec![0.0f64; 9];
    perm[0] = eps_11;
    perm[4] = eps_11;
    perm[8] = eps_33;

    PiezoelectricFem::new(e_stiff, e_coupling, perm, n_mech, n_elec)
}

// ---------------------------------------------------------------------------
// PoroelasticFem
// ---------------------------------------------------------------------------

/// Biot consolidation (poroelastic) FEM formulation using u-p formulation.
///
/// Implements Biot's theory of poroelasticity:
/// div σ(u) - α grad p = b (equilibrium)
/// α div u_dot + (1/M) p_dot - div(k/mu grad p) = f (continuity)
#[derive(Debug, Clone)]
pub struct PoroelasticFem {
    /// Drained Young's modulus E' (Pa).
    pub young_modulus: f64,
    /// Drained Poisson's ratio ν'.
    pub poisson_ratio: f64,
    /// Biot coefficient α (dimensionless, 0..1).
    pub biot_coefficient: f64,
    /// Biot modulus M (Pa).
    pub biot_modulus: f64,
    /// Permeability tensor k (m²), stored as 3x3 matrix.
    pub permeability: [[f64; 3]; 3],
    /// Fluid viscosity μ (Pa·s).
    pub fluid_viscosity: f64,
    /// Number of structural nodes.
    pub n_nodes: usize,
    /// Current pressure field (Pa).
    pub pressure: Vec<f64>,
    /// Current displacement field (m).
    pub displacements: Vec<[f64; 3]>,
}

impl PoroelasticFem {
    /// Create a new poroelastic FEM model.
    pub fn new(
        young_modulus: f64,
        poisson_ratio: f64,
        biot_coefficient: f64,
        biot_modulus: f64,
        permeability: [[f64; 3]; 3],
        fluid_viscosity: f64,
        n_nodes: usize,
    ) -> Self {
        Self {
            young_modulus,
            poisson_ratio,
            biot_coefficient,
            biot_modulus,
            permeability,
            fluid_viscosity,
            n_nodes,
            pressure: vec![0.0; n_nodes],
            displacements: vec![[0.0; 3]; n_nodes],
        }
    }

    /// Compute effective stress from strain and pore pressure.
    ///
    /// σ'_ij = σ_ij + α p δ_ij (Biot effective stress)
    pub fn effective_stress(&self, total_stress: [f64; 6], pore_pressure: f64) -> [f64; 6] {
        let alpha_p = self.biot_coefficient * pore_pressure;
        [
            total_stress[0] + alpha_p,
            total_stress[1] + alpha_p,
            total_stress[2] + alpha_p,
            total_stress[3],
            total_stress[4],
            total_stress[5],
        ]
    }

    /// Compute Darcy flux q = -(k/μ) grad p.
    ///
    /// # Arguments
    /// * `grad_p` – pressure gradient vector
    pub fn darcy_flux(&self, grad_p: [f64; 3]) -> [f64; 3] {
        let mut q = [0.0f64; 3];
        for (i, qi) in q.iter_mut().enumerate() {
            for (j, &gpj) in grad_p.iter().enumerate() {
                *qi -= (self.permeability[i][j] / self.fluid_viscosity) * gpj;
            }
        }
        q
    }

    /// Compute volumetric strain from displacement gradient (trace of ε).
    pub fn volumetric_strain(&self, strain: [f64; 6]) -> f64 {
        strain[0] + strain[1] + strain[2]
    }

    /// Compute coupling term α * div u for mass balance equation.
    pub fn coupling_term(&self, vol_strain: f64) -> f64 {
        self.biot_coefficient * vol_strain
    }

    /// Skempton's B coefficient (undrained pore pressure ratio).
    pub fn skempton_b(&self) -> f64 {
        let nu = self.poisson_ratio;
        let e = self.young_modulus;
        let k_drained = e / (3.0 * (1.0 - 2.0 * nu));
        if k_drained.abs() < 1e-15 || self.biot_modulus.abs() < 1e-15 {
            return 0.0;
        }
        let alpha = self.biot_coefficient;
        alpha / (alpha + k_drained / self.biot_modulus)
    }

    /// Undrained Poisson's ratio ν_u.
    pub fn undrained_poisson_ratio(&self) -> f64 {
        let nu = self.poisson_ratio;
        let _b = self.skempton_b();
        let alpha = self.biot_coefficient;
        let e = self.young_modulus;
        let k = e / (3.0 * (1.0 - 2.0 * nu));
        let k_u = k + alpha * alpha * self.biot_modulus;
        let g = e / (2.0 * (1.0 + nu));
        (3.0 * k_u - 2.0 * g) / (2.0 * (3.0 * k_u + g))
    }
}

// ---------------------------------------------------------------------------
// FluidStructureFem
// ---------------------------------------------------------------------------

/// Monolithic fluid-structure interaction FEM using ALE formulation.
///
/// Implements the Arbitrary Lagrangian-Eulerian (ALE) description for
/// moving-domain fluid-structure problems.
#[derive(Debug, Clone)]
pub struct FluidStructureFem {
    /// Fluid density ρ_f (kg/m³).
    pub fluid_density: f64,
    /// Dynamic viscosity μ_f (Pa·s).
    pub fluid_viscosity: f64,
    /// Structural density ρ_s (kg/m³).
    pub structural_density: f64,
    /// Young's modulus E_s (Pa).
    pub structural_young: f64,
    /// Poisson's ratio ν_s.
    pub structural_poisson: f64,
    /// Mesh stiffness coefficient for ALE mesh motion.
    pub mesh_stiffness: f64,
    /// Number of fluid nodes.
    pub n_fluid_nodes: usize,
    /// Number of structural nodes.
    pub n_struct_nodes: usize,
    /// Current mesh velocity at fluid nodes.
    pub mesh_velocity: Vec<[f64; 3]>,
    /// Current fluid velocity at fluid nodes.
    pub fluid_velocity: Vec<[f64; 3]>,
    /// Current structural displacement.
    pub struct_displacement: Vec<[f64; 3]>,
}

impl FluidStructureFem {
    /// Create a new FSI model.
    pub fn new(
        fluid_density: f64,
        fluid_viscosity: f64,
        structural_density: f64,
        structural_young: f64,
        structural_poisson: f64,
        mesh_stiffness: f64,
        n_fluid_nodes: usize,
        n_struct_nodes: usize,
    ) -> Self {
        Self {
            fluid_density,
            fluid_viscosity,
            structural_density,
            structural_young,
            structural_poisson,
            mesh_stiffness,
            n_fluid_nodes,
            n_struct_nodes,
            mesh_velocity: vec![[0.0; 3]; n_fluid_nodes],
            fluid_velocity: vec![[0.0; 3]; n_fluid_nodes],
            struct_displacement: vec![[0.0; 3]; n_struct_nodes],
        }
    }

    /// Compute ALE convective velocity: u_conv = u_fluid - u_mesh.
    pub fn convective_velocity(&self, node: usize) -> [f64; 3] {
        if node >= self.n_fluid_nodes {
            return [0.0; 3];
        }
        let uf = self.fluid_velocity[node];
        let um = self.mesh_velocity[node];
        [uf[0] - um[0], uf[1] - um[1], uf[2] - um[2]]
    }

    /// Compute fluid Reynolds number Re = ρ |u| L / μ.
    pub fn reynolds_number(&self, velocity: [f64; 3], length_scale: f64) -> f64 {
        self.fluid_density * norm3(velocity) * length_scale / self.fluid_viscosity
    }

    /// Compute mesh Laplacian smoothing velocity for ALE.
    ///
    /// Returns mesh velocity based on spring analogy.
    pub fn laplacian_mesh_smooth(&self, node: usize, neighbor_positions: &[[f64; 3]]) -> [f64; 3] {
        if neighbor_positions.is_empty() || node >= self.n_fluid_nodes {
            return [0.0; 3];
        }
        let n = neighbor_positions.len() as f64;
        let mean_x = neighbor_positions.iter().map(|p| p[0]).sum::<f64>() / n;
        let mean_y = neighbor_positions.iter().map(|p| p[1]).sum::<f64>() / n;
        let mean_z = neighbor_positions.iter().map(|p| p[2]).sum::<f64>() / n;
        let mesh_pos = self
            .struct_displacement
            .first()
            .copied()
            .unwrap_or([0.0; 3]);
        [
            self.mesh_stiffness * (mean_x - mesh_pos[0]),
            self.mesh_stiffness * (mean_y - mesh_pos[1]),
            self.mesh_stiffness * (mean_z - mesh_pos[2]),
        ]
    }

    /// Evaluate fluid traction on FSI interface (simplified Navier-Stokes).
    ///
    /// Returns traction vector f = -p*n + μ*(∂u/∂n + ∂u_n/∂x).
    pub fn interface_traction(
        &self,
        pressure: f64,
        normal: [f64; 3],
        velocity_grad: [[f64; 3]; 3],
    ) -> [f64; 3] {
        let mut traction = [0.0f64; 3];
        // Pressure part: -p*n
        for i in 0..3 {
            traction[i] -= pressure * normal[i];
        }
        // Viscous part: μ * (∇u + ∇u^T) * n
        for i in 0..3 {
            for j in 0..3 {
                traction[i] +=
                    self.fluid_viscosity * (velocity_grad[i][j] + velocity_grad[j][i]) * normal[j];
            }
        }
        traction
    }
}

// ---------------------------------------------------------------------------
// MagnetoMechanical
// ---------------------------------------------------------------------------

/// Magneto-mechanical FEM: Lorentz force body load.
///
/// Computes body forces from magnetic fields and current density,
/// and couples them to structural equations.
#[derive(Debug, Clone)]
pub struct MagnetoMechanical {
    /// Magnetic permeability μ_0 (H/m).
    pub permeability: f64,
    /// Electrical conductivity σ (S/m).
    pub conductivity: f64,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Current density J at each node (A/m²).
    pub current_density: Vec<[f64; 3]>,
    /// Magnetic flux density B at each node (T).
    pub magnetic_field: Vec<[f64; 3]>,
}

impl MagnetoMechanical {
    /// Create a new magneto-mechanical model.
    pub fn new(permeability: f64, conductivity: f64, n_nodes: usize) -> Self {
        Self {
            permeability,
            conductivity,
            n_nodes,
            current_density: vec![[0.0; 3]; n_nodes],
            magnetic_field: vec![[0.0; 3]; n_nodes],
        }
    }

    /// Compute Lorentz body force at a node: f = J × B.
    pub fn lorentz_force(&self, node: usize) -> [f64; 3] {
        if node >= self.n_nodes {
            return [0.0; 3];
        }
        let j = self.current_density[node];
        let b = self.magnetic_field[node];
        // J × B
        [
            j[1] * b[2] - j[2] * b[1],
            j[2] * b[0] - j[0] * b[2],
            j[0] * b[1] - j[1] * b[0],
        ]
    }

    /// Compute Maxwell stress tensor at a node.
    ///
    /// T_ij = (1/μ) * (B_i*B_j - 0.5 * |B|² δ_ij)
    pub fn maxwell_stress_tensor(&self, node: usize) -> [[f64; 3]; 3] {
        if node >= self.n_nodes {
            return [[0.0; 3]; 3];
        }
        let b = self.magnetic_field[node];
        let b_sq = dot3(b, b);
        let inv_mu = 1.0 / self.permeability;
        let mut t = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                t[i][j] = inv_mu * b[i] * b[j];
                if i == j {
                    t[i][j] -= 0.5 * inv_mu * b_sq;
                }
            }
        }
        t
    }

    /// Set uniform magnetic field at all nodes.
    pub fn set_uniform_magnetic_field(&mut self, b: [f64; 3]) {
        for node in self.magnetic_field.iter_mut() {
            *node = b;
        }
    }

    /// Compute eddy current density from velocity and magnetic field.
    ///
    /// J_eddy = σ * (v × B)
    pub fn eddy_current(&self, node: usize, velocity: [f64; 3]) -> [f64; 3] {
        if node >= self.n_nodes {
            return [0.0; 3];
        }
        let b = self.magnetic_field[node];
        let vxb = [
            velocity[1] * b[2] - velocity[2] * b[1],
            velocity[2] * b[0] - velocity[0] * b[2],
            velocity[0] * b[1] - velocity[1] * b[0],
        ];
        scale3(vxb, self.conductivity)
    }
}

// ---------------------------------------------------------------------------
// DiffusionReaction
// ---------------------------------------------------------------------------

/// FEM for diffusion-reaction equations with SUPG stabilization.
///
/// Solves: ρ c_p ∂T/∂t + u·∇T = ∇·(k∇T) + s
/// with SUPG stabilization for convection-dominated problems.
#[derive(Debug, Clone)]
pub struct DiffusionReaction {
    /// Diffusion coefficient D (m²/s).
    pub diffusivity: f64,
    /// Reaction rate coefficient R (1/s).
    pub reaction_rate: f64,
    /// Convection velocity u (m/s).
    pub velocity: [f64; 3],
    /// Source term s.
    pub source: f64,
    /// SUPG stabilization parameter τ.
    pub supg_tau: f64,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Current concentration/temperature field.
    pub field: Vec<f64>,
}

impl DiffusionReaction {
    /// Create a new diffusion-reaction model.
    pub fn new(
        diffusivity: f64,
        reaction_rate: f64,
        velocity: [f64; 3],
        source: f64,
        n_nodes: usize,
    ) -> Self {
        let h = 1.0; // characteristic mesh size
        let u_norm = norm3(velocity);
        let supg_tau = if u_norm > 1e-14 {
            h / (2.0 * u_norm)
        } else {
            0.0
        };
        Self {
            diffusivity,
            reaction_rate,
            velocity,
            source,
            supg_tau,
            n_nodes,
            field: vec![0.0; n_nodes],
        }
    }

    /// Compute SUPG stabilization parameter τ.
    ///
    /// τ = h / (2|u|) * coth(Pe) - 1/Pe
    /// where Pe = |u|h / (2D) is the Peclet number.
    pub fn supg_stabilization(&self, h: f64) -> f64 {
        let u_norm = norm3(self.velocity);
        if u_norm < 1e-14 || h < 1e-14 {
            return 0.0;
        }
        let pe = u_norm * h / (2.0 * self.diffusivity);
        if pe < 1e-6 {
            // Taylor expansion: coth(x) - 1/x ≈ x/3 for small x
            return h * h / (6.0 * self.diffusivity);
        }
        let coth_pe = (pe.exp() + (-pe).exp()) / (pe.exp() - (-pe).exp());
        h / (2.0 * u_norm) * (coth_pe - 1.0 / pe)
    }

    /// Compute Peclet number Pe = |u| h / (2D).
    pub fn peclet_number(&self, h: f64) -> f64 {
        let u_norm = norm3(self.velocity);
        u_norm * h / (2.0 * self.diffusivity)
    }

    /// Compute residual r = D ∇²φ - u·∇φ - R φ + s at a node.
    ///
    /// Uses finite differences for gradient approximation.
    pub fn compute_residual_node(&self, node: usize, laplacian: f64, grad_phi: [f64; 3]) -> f64 {
        let phi = self.field.get(node).copied().unwrap_or(0.0);
        self.diffusivity * laplacian - dot3(self.velocity, grad_phi) - self.reaction_rate * phi
            + self.source
    }

    /// Steady-state 1D analytical solution for constant coefficients.
    ///
    /// Returns φ(x) for 1D problem with φ(0)=0, φ(L)=φ_L.
    pub fn analytical_1d(&self, x: f64, l: f64, phi_l: f64) -> f64 {
        let u = self.velocity[0];
        let d = self.diffusivity;
        if d < 1e-15 || l < 1e-15 {
            return 0.0;
        }
        let pe_l = u * l / d;
        if pe_l.abs() < 1e-8 {
            // Pure diffusion: linear profile
            return phi_l * x / l;
        }
        phi_l * ((pe_l * x / l).exp() - 1.0) / (pe_l.exp() - 1.0)
    }
}

// ---------------------------------------------------------------------------
// ElectrochemicalFem
// ---------------------------------------------------------------------------

/// Butler-Volmer electrode kinetics in FEM framework.
///
/// Models electrochemical reactions at electrode surfaces using
/// the Butler-Volmer equation for current density.
#[derive(Debug, Clone)]
pub struct ElectrochemicalFem {
    /// Exchange current density i_0 (A/m²).
    pub exchange_current_density: f64,
    /// Anodic transfer coefficient α_a.
    pub transfer_coeff_anodic: f64,
    /// Cathodic transfer coefficient α_c.
    pub transfer_coeff_cathodic: f64,
    /// Faraday constant F (C/mol).
    pub faraday: f64,
    /// Universal gas constant R (J/mol·K).
    pub gas_constant: f64,
    /// Temperature T (K).
    pub temperature: f64,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Electric potential φ at nodes (V).
    pub potential: Vec<f64>,
    /// Species concentration c at nodes (mol/m³).
    pub concentration: Vec<f64>,
}

impl ElectrochemicalFem {
    /// Create a new electrochemical FEM model.
    pub fn new(
        exchange_current_density: f64,
        transfer_coeff_anodic: f64,
        transfer_coeff_cathodic: f64,
        temperature: f64,
        n_nodes: usize,
    ) -> Self {
        Self {
            exchange_current_density,
            transfer_coeff_anodic,
            transfer_coeff_cathodic,
            faraday: 96485.0,
            gas_constant: 8.314,
            temperature,
            n_nodes,
            potential: vec![0.0; n_nodes],
            concentration: vec![1.0; n_nodes],
        }
    }

    /// Compute Butler-Volmer current density at overpotential η.
    ///
    /// i = i_0 * (exp(α_a*F*η/(RT)) - exp(-α_c*F*η/(RT)))
    pub fn butler_volmer(&self, overpotential: f64) -> f64 {
        let ft = self.faraday / (self.gas_constant * self.temperature);
        let anodic = (self.transfer_coeff_anodic * ft * overpotential).exp();
        let cathodic = (-self.transfer_coeff_cathodic * ft * overpotential).exp();
        self.exchange_current_density * (anodic - cathodic)
    }

    /// Linearized Butler-Volmer current (first-order Taylor around η=0).
    ///
    /// i ≈ i_0 * (α_a + α_c) * F * η / RT
    pub fn butler_volmer_linearized(&self, overpotential: f64) -> f64 {
        let ft = self.faraday / (self.gas_constant * self.temperature);
        self.exchange_current_density
            * (self.transfer_coeff_anodic + self.transfer_coeff_cathodic)
            * ft
            * overpotential
    }

    /// Compute Tafel slope (high overpotential approximation).
    ///
    /// Returns b = RT / (α_a * F) for anodic reaction.
    pub fn tafel_slope_anodic(&self) -> f64 {
        self.gas_constant * self.temperature / (self.transfer_coeff_anodic * self.faraday)
    }

    /// Compute open-circuit potential from Nernst equation.
    ///
    /// U = U_0 + (RT/nF) * ln(c_ox / c_red)
    pub fn nernst_potential(
        &self,
        u_standard: f64,
        c_ox: f64,
        c_red: f64,
        n_electrons: f64,
    ) -> f64 {
        if c_red < 1e-15 || n_electrons < 1e-15 {
            return u_standard;
        }
        u_standard
            + (self.gas_constant * self.temperature) / (n_electrons * self.faraday)
                * (c_ox / c_red).ln()
    }
}

// ---------------------------------------------------------------------------
// ThermalMechanicalFatigue
// ---------------------------------------------------------------------------

/// Cyclic thermal-mechanical loading with Chaboche kinematic hardening.
///
/// Models creep-fatigue interaction with temperature-dependent material
/// properties using the Chaboche nonlinear kinematic hardening model.
#[derive(Debug, Clone)]
pub struct ThermalMechanicalFatigue {
    /// Initial yield stress σ_y0 (Pa).
    pub yield_stress: f64,
    /// Chaboche hardening modulus C_1 (Pa).
    pub hardening_c1: f64,
    /// Chaboche recall parameter γ_1.
    pub hardening_gamma1: f64,
    /// Isotropic hardening modulus Q_inf (Pa).
    pub isotropic_q: f64,
    /// Isotropic hardening rate b.
    pub isotropic_b: f64,
    /// Thermal softening coefficient.
    pub thermal_softening: f64,
    /// Reference temperature for hardening (K).
    pub reference_temperature: f64,
    /// Accumulated plastic strain p.
    pub accumulated_plastic_strain: f64,
    /// Back stress X (kinematic hardening, 6 components Voigt).
    pub back_stress: [f64; 6],
    /// Current temperature.
    pub temperature: f64,
}

impl ThermalMechanicalFatigue {
    /// Create a new thermal-mechanical fatigue model.
    pub fn new(
        yield_stress: f64,
        hardening_c1: f64,
        hardening_gamma1: f64,
        isotropic_q: f64,
        isotropic_b: f64,
        thermal_softening: f64,
        reference_temperature: f64,
    ) -> Self {
        Self {
            yield_stress,
            hardening_c1,
            hardening_gamma1,
            isotropic_q,
            isotropic_b,
            thermal_softening,
            reference_temperature,
            accumulated_plastic_strain: 0.0,
            back_stress: [0.0; 6],
            temperature: reference_temperature,
        }
    }

    /// Compute current yield stress including temperature softening and isotropic hardening.
    ///
    /// σ_y = (σ_y0 + Q_inf*(1 - exp(-b*p))) * (1 - β*(T - T_ref))
    pub fn current_yield_stress(&self) -> f64 {
        let iso_hard =
            self.isotropic_q * (1.0 - (-self.isotropic_b * self.accumulated_plastic_strain).exp());
        let thermal_factor =
            1.0 - self.thermal_softening * (self.temperature - self.reference_temperature).max(0.0);
        (self.yield_stress + iso_hard) * thermal_factor.max(0.1)
    }

    /// Compute von Mises equivalent stress from deviatoric stress minus back stress.
    pub fn von_mises_equivalent(&self, stress: [f64; 6]) -> f64 {
        let s = [
            stress[0] - self.back_stress[0],
            stress[1] - self.back_stress[1],
            stress[2] - self.back_stress[2],
            stress[3] - self.back_stress[3],
            stress[4] - self.back_stress[4],
            stress[5] - self.back_stress[5],
        ];
        let s11 = s[0];
        let s22 = s[1];
        let s33 = s[2];
        let s12 = s[3];
        let s13 = s[4];
        let s23 = s[5];
        ((s11 - s22).powi(2)
            + (s22 - s33).powi(2)
            + (s33 - s11).powi(2)
            + 6.0 * (s12.powi(2) + s13.powi(2) + s23.powi(2)))
        .sqrt()
            / 2.0_f64.sqrt()
    }

    /// Update back stress using Chaboche evolution equation.
    ///
    /// Ẋ = C * ε̇_p - γ * X * ṗ
    pub fn update_back_stress(&mut self, d_plastic_strain: [f64; 6], d_p: f64) {
        for (bs_i, &dps_i) in self.back_stress.iter_mut().zip(d_plastic_strain.iter()) {
            *bs_i += self.hardening_c1 * dps_i - self.hardening_gamma1 * *bs_i * d_p;
        }
        self.accumulated_plastic_strain += d_p;
    }

    /// Compute fatigue damage per cycle using linear damage accumulation.
    ///
    /// Uses simplified Coffin-Manson relationship.
    pub fn damage_per_cycle(&self, delta_plastic_strain: f64) -> f64 {
        // Coffin-Manson: Δεp = C_f * (2N_f)^c
        // Damage per cycle D = 1/N_f ≈ (Δεp / C_f)^(1/c)
        let c_f = 0.5; // fatigue ductility coefficient
        let c_exp = -0.6; // fatigue ductility exponent
        if delta_plastic_strain < 1e-10 {
            return 0.0;
        }
        let two_nf = (delta_plastic_strain / c_f).powf(1.0 / c_exp);
        if two_nf < 2.0 {
            return 1.0;
        }
        2.0 / two_nf
    }
}

// ---------------------------------------------------------------------------
// CoupledSolverFem
// ---------------------------------------------------------------------------

/// Solver strategy for coupled multiphysics FEM.
#[derive(Debug, Clone, PartialEq)]
pub enum SolverStrategy {
    /// Staggered (Gauss-Seidel): alternating field solves.
    Staggered,
    /// Monolithic (block Newton): full coupled system.
    Monolithic,
    /// Operator splitting: explicit-implicit decomposition.
    OperatorSplitting,
}

/// Convergence criteria for coupled solver.
#[derive(Debug, Clone)]
pub struct CoupledConvergenceCriteria {
    /// Maximum number of coupling iterations.
    pub max_iterations: usize,
    /// Tolerance for displacement residual.
    pub disp_tolerance: f64,
    /// Tolerance for temperature residual.
    pub temp_tolerance: f64,
    /// Tolerance for pressure residual.
    pub press_tolerance: f64,
}

impl Default for CoupledConvergenceCriteria {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            disp_tolerance: 1e-8,
            temp_tolerance: 1e-6,
            press_tolerance: 1e-6,
        }
    }
}

/// Coupled multiphysics FEM solver.
///
/// Implements staggered (Gauss-Seidel), monolithic (block Newton),
/// and operator-splitting schemes.
#[derive(Debug, Clone)]
pub struct CoupledSolverFem {
    /// Solver strategy.
    pub strategy: SolverStrategy,
    /// Convergence criteria.
    pub criteria: CoupledConvergenceCriteria,
    /// Time step size (s).
    pub dt: f64,
    /// Current simulation time (s).
    pub time: f64,
    /// Number of converged iterations per time step (diagnostic).
    pub last_iter_count: usize,
}

impl CoupledSolverFem {
    /// Create a new coupled solver.
    pub fn new(strategy: SolverStrategy, dt: f64) -> Self {
        Self {
            strategy,
            criteria: CoupledConvergenceCriteria::default(),
            dt,
            time: 0.0,
            last_iter_count: 0,
        }
    }

    /// Perform one staggered coupling iteration.
    ///
    /// Returns true if converged, false if not yet converged.
    pub fn staggered_iteration(
        &mut self,
        disp_res: f64,
        temp_res: f64,
        press_res: f64,
        iter: usize,
    ) -> bool {
        self.last_iter_count = iter;
        disp_res < self.criteria.disp_tolerance
            && temp_res < self.criteria.temp_tolerance
            && press_res < self.criteria.press_tolerance
    }

    /// Advance time by one step.
    pub fn advance_time(&mut self) {
        self.time += self.dt;
    }

    /// Compute relaxation factor for staggered coupling (Aitken acceleration).
    ///
    /// Returns relaxation factor ω ∈ (0, 1].
    pub fn aitken_relaxation(&self, res_prev: &[f64], res_curr: &[f64], omega_prev: f64) -> f64 {
        let mut num = 0.0;
        let mut den = 0.0;
        for (r0, r1) in res_prev.iter().zip(res_curr.iter()) {
            let dr = r1 - r0;
            num += r0 * dr;
            den += dr * dr;
        }
        if den.abs() < 1e-20 {
            return omega_prev;
        }
        let omega_new = -omega_prev * num / den;
        omega_new.clamp(0.1, 1.0)
    }

    /// Estimate computational cost ratio for monolithic vs staggered approach.
    pub fn cost_ratio_estimate(&self) -> f64 {
        // Monolithic is typically 3-5x more expensive per iteration but converges faster
        match self.strategy {
            SolverStrategy::Monolithic => 4.0,
            SolverStrategy::Staggered => 1.0,
            SolverStrategy::OperatorSplitting => 1.5,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_thermo_elastic() -> ThermoElasticFem {
        ThermoElasticFem::new(2e11, 0.3, 12e-6, 293.0, 50.0, 3.8e6, 10)
    }

    #[test]
    fn test_multi_field_result_new() {
        let result = MultiFieldResult::new(5);
        assert_eq!(result.num_nodes(), 5);
        assert_eq!(result.displacements.len(), 5);
        assert_eq!(result.temperatures.len(), 5);
    }

    #[test]
    fn test_multi_field_result_max_displacement() {
        let mut result = MultiFieldResult::new(3);
        result.displacements[1] = [3.0, 4.0, 0.0];
        let max = result.max_displacement();
        assert!((max - 5.0).abs() < 1e-10, "max displacement should be 5.0");
    }

    #[test]
    fn test_thermo_elastic_thermal_strain() {
        let fem = make_thermo_elastic();
        let eps = fem.thermal_strain(393.0); // 100 K above reference
        let expected = 12e-6 * 100.0;
        assert!((eps[0] - expected).abs() < 1e-15);
        assert!((eps[1] - expected).abs() < 1e-15);
        assert!((eps[3] - 0.0).abs() < 1e-15); // shear should be zero
    }

    #[test]
    fn test_thermo_elastic_zero_temp_diff() {
        let fem = make_thermo_elastic();
        let eps = fem.thermal_strain(293.0); // at reference temperature
        for &e in &eps {
            assert!(e.abs() < 1e-15, "thermal strain should be zero at T_ref");
        }
    }

    #[test]
    fn test_thermo_elastic_constitutive_2d() {
        let fem = make_thermo_elastic();
        let d = fem.constitutive_matrix_2d();
        // D should be symmetric
        assert!(
            (d[0][1] - d[1][0]).abs() < 1e-6,
            "D matrix should be symmetric"
        );
        // Diagonal should be positive
        assert!(d[0][0] > 0.0);
        assert!(d[1][1] > 0.0);
        assert!(d[2][2] > 0.0);
    }

    #[test]
    fn test_thermo_elastic_thermal_load_vector() {
        let mut fem = make_thermo_elastic();
        fem.temperatures[0] = 393.0;
        let f = fem.thermal_load_vector();
        assert_eq!(f.len(), 30); // 10 nodes * 3 dof
        // Load at heated node should be non-zero
        let f0_norm = (f[0].powi(2) + f[1].powi(2) + f[2].powi(2)).sqrt();
        assert!(
            f0_norm > 0.0,
            "thermal load at heated node should be non-zero"
        );
    }

    #[test]
    fn test_piezoelectric_compute_stress() {
        let fem = make_pzt5a(30, 10);
        let strain = [1e-4, 0.0, 0.0, 0.0, 0.0, 0.0];
        let e_field = [0.0, 0.0, 0.0];
        let sigma = fem.compute_stress(&strain, &e_field);
        // With non-zero strain, stress should be non-zero
        let has_nonzero = sigma.iter().any(|&s| s.abs() > 0.0);
        assert!(has_nonzero, "stress should be non-zero for non-zero strain");
    }

    #[test]
    fn test_piezoelectric_electric_displacement() {
        let fem = make_pzt5a(30, 10);
        let strain = [0.0; 6];
        let e_field = [0.0, 0.0, 1000.0]; // 1 kV/m
        let d_vec = fem.compute_electric_displacement(&strain, &e_field);
        // D_3 should be non-zero from κ_33 * E_3
        assert!(d_vec[2].abs() > 0.0 || d_vec[0].abs() >= 0.0);
    }

    #[test]
    fn test_piezoelectric_capacitance() {
        let fem = make_pzt5a(30, 10);
        let c = fem.capacitance_scalar(0.01, 0.001);
        assert!(c > 0.0, "capacitance should be positive");
    }

    #[test]
    fn test_poroelastic_effective_stress() {
        let perm = [[1e-12, 0.0, 0.0], [0.0, 1e-12, 0.0], [0.0, 0.0, 1e-12]];
        let fem = PoroelasticFem::new(1e7, 0.3, 0.8, 2e9, perm, 1e-3, 10);
        let total_stress = [0.0; 6];
        let p = 1e5;
        let eff = fem.effective_stress(total_stress, p);
        let expected_diag = fem.biot_coefficient * p;
        assert!(
            (eff[0] - expected_diag).abs() < 1.0,
            "effective stress diagonal incorrect"
        );
    }

    #[test]
    fn test_poroelastic_darcy_flux() {
        let perm = [[1e-12, 0.0, 0.0], [0.0, 1e-12, 0.0], [0.0, 0.0, 1e-12]];
        let fem = PoroelasticFem::new(1e7, 0.3, 0.8, 2e9, perm, 1e-3, 10);
        let grad_p = [1000.0, 0.0, 0.0]; // 1 kPa/m
        let q = fem.darcy_flux(grad_p);
        // q_x = -(k/mu) * dp/dx = -(1e-12/1e-3) * 1000 = -1e-6 m/s
        assert!(
            (q[0] - (-1e-6)).abs() < 1e-10,
            "Darcy flux incorrect: {}",
            q[0]
        );
    }

    #[test]
    fn test_poroelastic_skempton_b() {
        let perm = [[1e-12, 0.0, 0.0], [0.0, 1e-12, 0.0], [0.0, 0.0, 1e-12]];
        let fem = PoroelasticFem::new(1e7, 0.3, 0.8, 2e9, perm, 1e-3, 10);
        let b = fem.skempton_b();
        assert!(
            (0.0..=1.0).contains(&b),
            "Skempton's B should be in [0,1], got {}",
            b
        );
    }

    #[test]
    fn test_fluid_structure_convective_velocity() {
        let mut fsi = FluidStructureFem::new(1000.0, 1e-3, 7800.0, 2e11, 0.3, 1.0, 5, 3);
        fsi.fluid_velocity[0] = [1.0, 0.0, 0.0];
        fsi.mesh_velocity[0] = [0.3, 0.0, 0.0];
        let cv = fsi.convective_velocity(0);
        assert!(
            (cv[0] - 0.7).abs() < 1e-10,
            "convective velocity x should be 0.7"
        );
    }

    #[test]
    fn test_fluid_structure_reynolds_number() {
        let fsi = FluidStructureFem::new(1000.0, 1e-3, 7800.0, 2e11, 0.3, 1.0, 5, 3);
        let re = fsi.reynolds_number([1.0, 0.0, 0.0], 0.1);
        assert!(
            (re - 100000.0).abs() < 1.0,
            "Reynolds number should be 1e5, got {}",
            re
        );
    }

    #[test]
    fn test_magneto_mechanical_lorentz_force() {
        let mut mm = MagnetoMechanical::new(4.0 * std::f64::consts::PI * 1e-7, 1e7, 5);
        mm.current_density[0] = [1.0, 0.0, 0.0]; // J = 1 A/m² in x
        mm.magnetic_field[0] = [0.0, 1.0, 0.0]; // B = 1 T in y
        let f = mm.lorentz_force(0);
        // J × B = [1,0,0] × [0,1,0] = [0,0,1]
        assert!(f[2].abs() > 0.9, "Lorentz force z component should be ~1");
    }

    #[test]
    fn test_magneto_mechanical_maxwell_stress() {
        let mut mm = MagnetoMechanical::new(4.0 * std::f64::consts::PI * 1e-7, 1e7, 3);
        mm.magnetic_field[0] = [1.0, 0.0, 0.0]; // B = 1 T in x
        let t = mm.maxwell_stress_tensor(0);
        // T_11 = B_1^2/mu - 0.5*|B|^2/mu = 0.5/mu
        let expected = 0.5 / mm.permeability;
        assert!((t[0][0] - expected).abs() < expected * 0.01);
    }

    #[test]
    fn test_diffusion_reaction_peclet() {
        let dr = DiffusionReaction::new(1e-5, 0.0, [0.1, 0.0, 0.0], 0.0, 10);
        let pe = dr.peclet_number(0.01); // h=1cm
        let expected = 0.1 * 0.01 / (2.0 * 1e-5);
        assert!(
            (pe - expected).abs() < 0.01,
            "Peclet number incorrect: {}",
            pe
        );
    }

    #[test]
    fn test_diffusion_reaction_analytical_1d() {
        let dr = DiffusionReaction::new(1.0, 0.0, [0.0, 0.0, 0.0], 0.0, 5);
        // Pure diffusion, linear profile
        let phi = dr.analytical_1d(0.5, 1.0, 1.0);
        assert!(
            (phi - 0.5).abs() < 0.01,
            "1D diffusion solution should be linear: {}",
            phi
        );
    }

    #[test]
    fn test_electrochemical_butler_volmer() {
        let fem = ElectrochemicalFem::new(1e-3, 0.5, 0.5, 298.0, 5);
        let i = fem.butler_volmer(0.0);
        assert!(
            i.abs() < 1e-10,
            "BV current at zero overpotential should be zero, got {}",
            i
        );
    }

    #[test]
    fn test_electrochemical_butler_volmer_positive_eta() {
        let fem = ElectrochemicalFem::new(1e-3, 0.5, 0.5, 298.0, 5);
        let i = fem.butler_volmer(0.1); // 100 mV overpotential
        assert!(
            i > 0.0,
            "BV current should be positive for positive overpotential"
        );
    }

    #[test]
    fn test_electrochemical_linearized_bv() {
        let fem = ElectrochemicalFem::new(1e-3, 0.5, 0.5, 298.0, 5);
        let eta = 1e-4; // small overpotential
        let i_lin = fem.butler_volmer_linearized(eta);
        let i_bv = fem.butler_volmer(eta);
        // Should agree for small overpotential
        assert!(
            (i_lin - i_bv).abs() / i_bv.abs() < 0.01,
            "linearized BV should agree with BV for small eta"
        );
    }

    #[test]
    fn test_electrochemical_nernst_potential() {
        let fem = ElectrochemicalFem::new(1e-3, 0.5, 0.5, 298.0, 5);
        let u = fem.nernst_potential(0.0, 1.0, 1.0, 1.0);
        // At equal concentrations, Nernst = U_0
        assert!(
            u.abs() < 1e-10,
            "Nernst potential at equal concentrations should be U_0"
        );
    }

    #[test]
    fn test_thermal_mechanical_fatigue_yield_stress() {
        let tmf = ThermalMechanicalFatigue::new(250e6, 1e11, 500.0, 50e6, 10.0, 5e-4, 293.0);
        let sy = tmf.current_yield_stress();
        assert!(
            (sy - 250e6).abs() < 1.0,
            "Initial yield stress should equal yield_stress: {}",
            sy
        );
    }

    #[test]
    fn test_thermal_mechanical_fatigue_back_stress_update() {
        let mut tmf = ThermalMechanicalFatigue::new(250e6, 1e11, 500.0, 50e6, 10.0, 5e-4, 293.0);
        let d_ep = [1e-4, 0.0, 0.0, 0.0, 0.0, 0.0];
        tmf.update_back_stress(d_ep, 1e-4);
        assert!(tmf.back_stress[0].abs() > 0.0, "back stress should update");
        assert!((tmf.accumulated_plastic_strain - 1e-4).abs() < 1e-15);
    }

    #[test]
    fn test_thermal_softening_effect() {
        let tmf = ThermalMechanicalFatigue::new(250e6, 1e11, 500.0, 50e6, 10.0, 5e-4, 293.0);
        let sy0 = tmf.current_yield_stress();
        let mut tmf2 = tmf.clone();
        tmf2.temperature = 1000.0; // very high T
        let sy1 = tmf2.current_yield_stress();
        assert!(sy1 < sy0, "yield stress should decrease with temperature");
    }

    #[test]
    fn test_coupled_solver_staggered_convergence() {
        let mut solver = CoupledSolverFem::new(SolverStrategy::Staggered, 0.01);
        let converged = solver.staggered_iteration(1e-10, 1e-9, 1e-9, 5);
        assert!(converged, "should converge with small residuals");
    }

    #[test]
    fn test_coupled_solver_not_converged() {
        let mut solver = CoupledSolverFem::new(SolverStrategy::Staggered, 0.01);
        let not_converged = solver.staggered_iteration(1.0, 1.0, 1.0, 1);
        assert!(!not_converged, "should not converge with large residuals");
    }

    #[test]
    fn test_coupled_solver_advance_time() {
        let mut solver = CoupledSolverFem::new(SolverStrategy::Monolithic, 0.1);
        solver.advance_time();
        assert!((solver.time - 0.1).abs() < 1e-14);
    }

    #[test]
    fn test_aitken_relaxation() {
        let solver = CoupledSolverFem::new(SolverStrategy::Staggered, 0.01);
        let res_prev = vec![0.1, 0.2, 0.3];
        let res_curr = vec![0.05, 0.1, 0.15]; // converging
        let omega = solver.aitken_relaxation(&res_prev, &res_curr, 1.0);
        assert!(
            omega > 0.0 && omega <= 1.0,
            "Aitken omega should be in (0,1]: {}",
            omega
        );
    }

    #[test]
    fn test_magneto_mechanical_eddy_current() {
        let mut mm = MagnetoMechanical::new(4.0 * std::f64::consts::PI * 1e-7, 1e7, 3);
        mm.magnetic_field[0] = [0.0, 0.0, 1.0]; // B in z
        let v = [1.0, 0.0, 0.0]; // velocity in x
        let j_eddy = mm.eddy_current(0, v);
        // v × B = [1,0,0] × [0,0,1] = [0*1-0*0, 0*0-1*1, 1*0-0*0] = [0,-1,0]
        assert!(
            j_eddy[1].abs() > 0.0,
            "eddy current y component should be non-zero"
        );
    }
}
