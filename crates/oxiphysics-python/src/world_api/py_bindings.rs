// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! PyO3-style binding wrappers for MD, LBM, FEM, SPH simulations
//! and material query helpers.

#![allow(missing_docs)]

use super::lbm::{PyLbmConfig, PyLbmGrid};
use super::sph::{PySphConfig, PySphSim};

// ===========================================================================
// MD Simulation Bindings (PyO3-style API wrappers)
// ===========================================================================

use crate::md_api::{PyMdConfig, PyMdSimulation};

/// PyO3-style binding handle for an MD simulation.
///
/// Wraps `PyMdSimulation` and exposes flat-array query helpers that mirror the
/// `#[pymethods]` surface that would be used with the real PyO3 macro
/// expansion.  All methods take and return primitive Rust types so that
/// future `#[pyfunction]` / `#[pymethods]` annotations require no changes.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyMdBinding {
    sim: PyMdSimulation,
}

impl PyMdBinding {
    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// Create a new MD binding with the default (argon-reduced) configuration.
    pub fn new_default() -> Self {
        Self {
            sim: PyMdSimulation::new(PyMdConfig::default()),
        }
    }

    /// Create from explicit configuration.
    pub fn new(config: PyMdConfig) -> Self {
        Self {
            sim: PyMdSimulation::new(config),
        }
    }

    // ------------------------------------------------------------------
    // Particle management
    // ------------------------------------------------------------------

    /// Add an atom at `[x, y, z]` with optional velocity `[vx, vy, vz]`.
    ///
    /// Returns the index of the newly added atom.
    #[allow(clippy::too_many_arguments)]
    pub fn add_atom(
        &mut self,
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        atom_type: u32,
    ) -> usize {
        let idx = self.sim.add_atom([x, y, z], atom_type);
        self.sim.set_velocity(idx, [vx, vy, vz]);
        idx
    }

    /// Return the number of atoms.
    pub fn atom_count(&self) -> usize {
        self.sim.atom_count()
    }

    // ------------------------------------------------------------------
    // Stepping
    // ------------------------------------------------------------------

    /// Advance the simulation by `dt` (reduced units or SI, consistent with config).
    pub fn step(&mut self, dt: f64) {
        self.sim.step(dt);
    }

    /// Advance the simulation by `n_steps` steps, each of size `dt`.
    pub fn run_steps(&mut self, n_steps: u32, dt: f64) {
        for _ in 0..n_steps {
            self.sim.step(dt);
        }
    }

    // ------------------------------------------------------------------
    // State queries
    // ------------------------------------------------------------------

    /// Return positions as a flat `[x0,y0,z0, x1,y1,z1, …]` array.
    pub fn get_positions_flat(&self) -> Vec<f64> {
        self.sim.all_positions()
    }

    /// Return velocities as a flat array.
    pub fn get_velocities_flat(&self) -> Vec<f64> {
        self.sim.all_velocities()
    }

    /// Return `[kinetic, potential, total]` energy.
    pub fn get_energies(&self) -> [f64; 3] {
        let ke = self.sim.kinetic_energy();
        let pe = self.sim.potential_energy();
        [ke, pe, ke + pe]
    }

    /// Instantaneous temperature in reduced units (from equipartition theorem).
    pub fn get_temperature(&self) -> f64 {
        self.sim.temperature()
    }

    /// Total simulation time elapsed.
    pub fn elapsed_time(&self) -> f64 {
        self.sim.time()
    }

    /// Number of completed time steps.
    pub fn step_count(&self) -> u64 {
        self.sim.step_count()
    }

    // ------------------------------------------------------------------
    // Thermostat control
    // ------------------------------------------------------------------

    /// Enable or disable the velocity-rescaling thermostat at runtime.
    pub fn set_thermostat_enabled(&mut self, enabled: bool) {
        self.sim.set_thermostat(enabled);
    }

    /// Set the target temperature (reduced units).
    pub fn set_target_temperature(&mut self, temp: f64) {
        self.sim.set_target_temperature(temp);
    }
}

// ===========================================================================
// LBM Simulation Bindings (PyO3-style API wrappers)
// ===========================================================================

/// PyO3-style binding handle for a 2-D LBM simulation.
///
/// Uses the `PyLbmGrid` type (D2Q9 BGK) defined earlier in this file.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyLbmBinding {
    grid: PyLbmGrid,
}

impl PyLbmBinding {
    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// Create a new LBM binding with the given grid dimensions and viscosity.
    pub fn new(width: usize, height: usize, viscosity: f64) -> Self {
        let config = PyLbmConfig::new(width, height, viscosity);
        Self {
            grid: PyLbmGrid::new(&config),
        }
    }

    /// Create from an explicit `PyLbmConfig`.
    pub fn from_config(config: PyLbmConfig) -> Self {
        Self {
            grid: PyLbmGrid::new(&config),
        }
    }

    // ------------------------------------------------------------------
    // Grid dimensions
    // ------------------------------------------------------------------

    /// Grid width (number of cells along X).
    pub fn grid_width(&self) -> usize {
        self.grid.width()
    }

    /// Grid height (number of cells along Y).
    pub fn grid_height(&self) -> usize {
        self.grid.height()
    }

    /// Total number of cells.
    pub fn cell_count(&self) -> usize {
        self.grid.width() * self.grid.height()
    }

    // ------------------------------------------------------------------
    // Stepping
    // ------------------------------------------------------------------

    /// Advance the simulation by one LBM step.
    pub fn step(&mut self) {
        self.grid.step();
    }

    /// Advance the simulation by `n_steps` steps.
    pub fn run_steps(&mut self, n_steps: u32) {
        for _ in 0..n_steps {
            self.grid.step();
        }
    }

    // ------------------------------------------------------------------
    // Field queries
    // ------------------------------------------------------------------

    /// Return the density field as a flat row-major array.
    pub fn get_density_flat(&self) -> Vec<f64> {
        self.grid.density_field()
    }

    /// Return velocity field as interleaved `[ux0,uy0, ux1,uy1, …]`.
    pub fn get_velocity_flat(&self) -> Vec<f64> {
        self.grid.velocity_field()
    }

    /// Return velocity magnitude (speed) for each cell as a flat array.
    pub fn get_speed_flat(&self) -> Vec<f64> {
        let vel = self.grid.velocity_field();
        vel.chunks(2)
            .map(|c| (c[0] * c[0] + c[1] * c[1]).sqrt())
            .collect()
    }

    /// Mean density across all cells.
    pub fn mean_density(&self) -> f64 {
        let rho = self.get_density_flat();
        if rho.is_empty() {
            return 0.0;
        }
        rho.iter().sum::<f64>() / rho.len() as f64
    }

    /// Maximum speed across all cells.
    pub fn max_speed(&self) -> f64 {
        self.get_speed_flat()
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
            .max(0.0)
    }

    /// Number of completed steps.
    pub fn step_count(&self) -> u64 {
        self.grid.step_count()
    }
}

// ===========================================================================
// FEM Simulation Bindings (PyO3-style API wrappers)
// ===========================================================================

use crate::fem_api::{PyFemDirichletBC, PyFemMaterial, PyFemMesh, PyFemSolveResult, PyFemSolver};

/// PyO3-style binding handle for a 2-D FEM solver.
///
/// Wraps `PyFemSolver` and `PyFemMesh` and exposes helpers suitable for
/// future `#[pymethods]` annotation.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyFemBinding {
    mesh: PyFemMesh,
    solver: PyFemSolver,
    last_result: Option<PyFemSolveResult>,
}

impl PyFemBinding {
    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// Create a new FEM binding with an empty mesh.
    pub fn new() -> Self {
        Self {
            mesh: PyFemMesh::new(),
            solver: PyFemSolver::new(),
            last_result: None,
        }
    }

    // ------------------------------------------------------------------
    // Mesh construction helpers
    // ------------------------------------------------------------------

    /// Add a 2-D node at position `(x, y)`.  Returns the new node index.
    pub fn add_node(&mut self, x: f64, y: f64) -> usize {
        self.mesh.add_node(x, y)
    }

    /// Add a CST triangular element given three node indices.
    ///
    /// Returns the new element index.
    pub fn add_tri_element(&mut self, n0: usize, n1: usize, n2: usize) -> usize {
        self.mesh.add_element(n0, n1, n2)
    }

    /// Add a CST triangular element with explicit material and thickness.
    #[allow(clippy::too_many_arguments)]
    pub fn add_tri_element_with_material(
        &mut self,
        n0: usize,
        n1: usize,
        n2: usize,
        material_id: usize,
        thickness: f64,
    ) -> usize {
        self.mesh
            .add_element_with_material(n0, n1, n2, material_id, thickness)
    }

    /// Add a material to the mesh library. Returns the material index.
    pub fn add_material(&mut self, material: PyFemMaterial) -> usize {
        self.mesh.add_material(material)
    }

    /// Return the number of nodes.
    pub fn node_count(&self) -> usize {
        self.mesh.nodes.len()
    }

    /// Return the number of elements.
    pub fn element_count(&self) -> usize {
        self.mesh.elements.len()
    }

    // ------------------------------------------------------------------
    // Boundary conditions and loads
    // ------------------------------------------------------------------

    /// Pin a node (fix both x and y DOFs to zero).
    pub fn pin_node(&mut self, node: usize) {
        self.mesh.pin_node(node);
    }

    /// Apply a Dirichlet (fixed DOF) boundary condition.
    pub fn add_dirichlet_bc(&mut self, bc: PyFemDirichletBC) {
        self.mesh.add_dirichlet_bc(bc);
    }

    /// Apply a nodal force `[fx, fy]` at `node`.
    pub fn apply_nodal_force(&mut self, node: usize, fx: f64, fy: f64) {
        self.mesh.apply_nodal_force(node, fx, fy);
    }

    // ------------------------------------------------------------------
    // Solve
    // ------------------------------------------------------------------

    /// Solve the linear-elastic static problem.
    ///
    /// Returns `true` on success.  Query results via `get_displacement_flat`.
    pub fn solve(&mut self) -> bool {
        match self.solver.solve(&self.mesh) {
            Some(result) => {
                self.last_result = Some(result);
                true
            }
            None => false,
        }
    }

    /// Return nodal displacements as a flat `[ux0,uy0, ux1,uy1, …]`
    /// array, or an empty `Vec` if no solve has been run.
    pub fn get_displacement_flat(&self) -> Vec<f64> {
        match &self.last_result {
            Some(r) => r.displacements.clone(),
            None => Vec::new(),
        }
    }

    /// Return element von-Mises stresses as a flat array.
    pub fn get_stress_flat(&self) -> Vec<f64> {
        match &self.last_result {
            Some(r) => r.von_mises_stress.clone(),
            None => Vec::new(),
        }
    }

    /// Return the maximum von-Mises stress, or `0.0` if no result exists.
    pub fn max_von_mises_stress(&self) -> f64 {
        match &self.last_result {
            Some(r) if !r.von_mises_stress.is_empty() => r.max_von_mises().max(0.0),
            _ => 0.0,
        }
    }

    /// Return whether the last solve converged.
    pub fn last_solve_converged(&self) -> bool {
        self.last_result
            .as_ref()
            .map(|r| r.converged)
            .unwrap_or(false)
    }
}

impl Default for PyFemBinding {
    fn default() -> Self {
        Self::new()
    }
}

// ===========================================================================
// SPH Simulation Bindings (PyO3-style API wrappers)
// ===========================================================================

/// PyO3-style binding handle for a 3-D SPH simulation.
///
/// Uses the `PySphSim` type (WCSPH) defined earlier in this file.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PySphBinding {
    sim: PySphSim,
}

impl PySphBinding {
    // ------------------------------------------------------------------
    // Construction
    // ------------------------------------------------------------------

    /// Create a new SPH binding with the water-like default configuration.
    pub fn new_water() -> Self {
        Self {
            sim: PySphSim::new(PySphConfig::water()),
        }
    }

    /// Create from an explicit `PySphConfig`.
    pub fn from_config(config: PySphConfig) -> Self {
        Self {
            sim: PySphSim::new(config),
        }
    }

    // ------------------------------------------------------------------
    // Particle management
    // ------------------------------------------------------------------

    /// Add a fluid particle at `[x, y, z]` with zero initial velocity.
    ///
    /// Returns the new particle index.
    pub fn add_particle_at(&mut self, x: f64, y: f64, z: f64) -> usize {
        let idx = self.sim.particle_count();
        self.sim.add_particle([x, y, z]);
        idx
    }

    /// Add a fluid particle at `[x, y, z]` with the given initial velocity
    /// (`vx`, `vy`, `vz`).
    pub fn add_particle(&mut self, x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64) -> usize {
        let idx = self.sim.particle_count();
        self.sim.add_particle([x, y, z]);
        // Set velocity via the position-then-velocity pattern
        if let Some(vel) = self.sim.velocities_mut().get_mut(idx) {
            *vel = [vx, vy, vz];
        }
        idx
    }

    /// Return the number of particles.
    pub fn particle_count(&self) -> usize {
        self.sim.particle_count()
    }

    // ------------------------------------------------------------------
    // Stepping
    // ------------------------------------------------------------------

    /// Advance the simulation by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        self.sim.step(dt);
    }

    /// Advance by `n_steps` steps.
    pub fn run_steps(&mut self, n_steps: u32, dt: f64) {
        for _ in 0..n_steps {
            self.sim.step(dt);
        }
    }

    // ------------------------------------------------------------------
    // State queries
    // ------------------------------------------------------------------

    /// Return positions as a flat `[x0,y0,z0, x1,y1,z1, …]` array.
    pub fn get_positions_flat(&self) -> Vec<f64> {
        self.sim.all_positions()
    }

    /// Return velocities as a flat array.
    pub fn get_velocities_flat(&self) -> Vec<f64> {
        self.sim.all_velocities()
    }

    /// Return densities as a flat array (one value per particle).
    pub fn get_densities_flat(&self) -> Vec<f64> {
        (0..self.sim.particle_count())
            .filter_map(|i| self.sim.density(i))
            .collect()
    }

    /// Return pressures as a flat array.
    pub fn get_pressures_flat(&self) -> Vec<f64> {
        (0..self.sim.particle_count())
            .filter_map(|i| self.sim.pressure(i))
            .collect()
    }

    /// Return speeds (velocity magnitudes) as a flat array.
    pub fn get_speeds_flat(&self) -> Vec<f64> {
        let vels = self.sim.all_velocities();
        vels.chunks(3)
            .map(|c| (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt())
            .collect()
    }

    /// Mean density of all particles.
    pub fn mean_density(&self) -> f64 {
        let n = self.sim.particle_count();
        if n == 0 {
            return 0.0;
        }
        let sum: f64 = (0..n).filter_map(|i| self.sim.density(i)).sum();
        sum / n as f64
    }

    /// Maximum speed across all particles.
    pub fn max_speed(&self) -> f64 {
        self.get_speeds_flat()
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
            .max(0.0)
    }

    /// Total kinetic energy of all particles.
    pub fn total_kinetic_energy(&self) -> f64 {
        let mass = self.sim.particle_mass();
        let vels = self.sim.all_velocities();
        vels.chunks(3)
            .map(|c| 0.5 * mass * (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]))
            .sum()
    }

    // ------------------------------------------------------------------
    // Configuration queries
    // ------------------------------------------------------------------

    /// Smoothing length h.
    pub fn smoothing_length(&self) -> f64 {
        self.sim.smoothing_length()
    }

    /// Rest density ρ₀.
    pub fn rest_density(&self) -> f64 {
        self.sim.rest_density()
    }
}

// ===========================================================================
// Material query helpers (PyO3-style standalone functions)
// ===========================================================================

/// Material property record returned by query functions.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub struct MaterialProperties {
    /// Material name.
    pub name: String,
    /// Density in kg/m³.
    pub density: f64,
    /// Young's modulus in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub poisson_ratio: f64,
    /// Ultimate tensile strength in Pa.
    pub tensile_strength: f64,
    /// Dynamic viscosity in Pa·s (0 for solids).
    pub viscosity: f64,
    /// Yield strength in Pa (0 if not applicable).
    pub yield_strength: f64,
}

impl MaterialProperties {
    /// Steel (structural grade).
    pub fn steel() -> Self {
        Self {
            name: "steel".into(),
            density: 7850.0,
            youngs_modulus: 200e9,
            poisson_ratio: 0.3,
            tensile_strength: 400e6,
            viscosity: 0.0,
            yield_strength: 250e6,
        }
    }

    /// Aluminium alloy 6061-T6.
    pub fn aluminium() -> Self {
        Self {
            name: "aluminium".into(),
            density: 2700.0,
            youngs_modulus: 69e9,
            poisson_ratio: 0.33,
            tensile_strength: 310e6,
            viscosity: 0.0,
            yield_strength: 276e6,
        }
    }

    /// Concrete (typical structural).
    pub fn concrete() -> Self {
        Self {
            name: "concrete".into(),
            density: 2400.0,
            youngs_modulus: 30e9,
            poisson_ratio: 0.2,
            tensile_strength: 3e6,
            viscosity: 0.0,
            yield_strength: 0.0,
        }
    }

    /// Natural rubber.
    pub fn rubber() -> Self {
        Self {
            name: "rubber".into(),
            density: 1100.0,
            youngs_modulus: 0.01e9,
            poisson_ratio: 0.49,
            tensile_strength: 20e6,
            viscosity: 0.0,
            yield_strength: 0.0,
        }
    }

    /// Water at 20 °C.
    pub fn water() -> Self {
        Self {
            name: "water".into(),
            density: 998.0,
            youngs_modulus: 0.0,
            poisson_ratio: 0.0,
            tensile_strength: 0.0,
            viscosity: 1.002e-3,
            yield_strength: 0.0,
        }
    }

    /// Air at 20 °C, 1 atm.
    pub fn air() -> Self {
        Self {
            name: "air".into(),
            density: 1.204,
            youngs_modulus: 0.0,
            poisson_ratio: 0.0,
            tensile_strength: 0.0,
            viscosity: 1.825e-5,
            yield_strength: 0.0,
        }
    }

    /// Titanium Ti-6Al-4V.
    pub fn titanium() -> Self {
        Self {
            name: "titanium".into(),
            density: 4430.0,
            youngs_modulus: 114e9,
            poisson_ratio: 0.34,
            tensile_strength: 950e6,
            viscosity: 0.0,
            yield_strength: 880e6,
        }
    }

    /// Lookup a material by name (case-insensitive).
    ///
    /// Returns `None` if the name is unknown.
    pub fn lookup(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "steel" => Some(Self::steel()),
            "aluminium" | "aluminum" => Some(Self::aluminium()),
            "concrete" => Some(Self::concrete()),
            "rubber" => Some(Self::rubber()),
            "water" => Some(Self::water()),
            "air" => Some(Self::air()),
            "titanium" => Some(Self::titanium()),
            _ => None,
        }
    }

    /// Compute the speed of longitudinal acoustic wave (P-wave) in m/s.
    ///
    /// Returns `None` for fluids where `youngs_modulus == 0`.
    pub fn p_wave_speed(&self) -> Option<f64> {
        if self.youngs_modulus < 1.0 || self.density < 1e-6 {
            return None;
        }
        let nu = self.poisson_ratio;
        let e = self.youngs_modulus;
        let rho = self.density;
        let c = ((e * (1.0 - nu)) / (rho * (1.0 + nu) * (1.0 - 2.0 * nu))).sqrt();
        Some(c)
    }

    /// Compute the shear wave speed (S-wave) in m/s.
    ///
    /// Returns `None` for fluids.
    pub fn s_wave_speed(&self) -> Option<f64> {
        if self.youngs_modulus < 1.0 || self.density < 1e-6 {
            return None;
        }
        let nu = self.poisson_ratio;
        let e = self.youngs_modulus;
        let rho = self.density;
        let g = e / (2.0 * (1.0 + nu));
        Some((g / rho).sqrt())
    }
}

// PyO3 standalone function mirrors (named with `py_` prefix to avoid future
// name conflicts when the real `#[pyfunction]` wrappers are generated).

/// PyO3-callable: look up a material and return `[density, E, nu, UTS, visc, yield]`
/// or an empty vec if the name is unknown.
#[allow(dead_code)]
pub fn py_query_material(name: &str) -> Vec<f64> {
    match MaterialProperties::lookup(name) {
        Some(m) => vec![
            m.density,
            m.youngs_modulus,
            m.poisson_ratio,
            m.tensile_strength,
            m.viscosity,
            m.yield_strength,
        ],
        None => Vec::new(),
    }
}

/// PyO3-callable: return P-wave speed for a named material, or `0.0`.
#[allow(dead_code)]
pub fn py_p_wave_speed(name: &str) -> f64 {
    MaterialProperties::lookup(name)
        .and_then(|m| m.p_wave_speed())
        .unwrap_or(0.0)
}

/// PyO3-callable: return S-wave speed for a named material, or `0.0`.
#[allow(dead_code)]
pub fn py_s_wave_speed(name: &str) -> f64 {
    MaterialProperties::lookup(name)
        .and_then(|m| m.s_wave_speed())
        .unwrap_or(0.0)
}
