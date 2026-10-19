// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Thermoelastic stress analysis.
//!
//! Provides thermoelastic coupling for FEM analysis including:
//! - Thermal strain computation from temperature change
//! - Plane-stress and plane-strain thermal stress vectors
//! - 1-D and 2-D thermoelastic element stiffness and load vectors
//! - One-way and two-way thermal–mechanical coupling
//! - Von Mises and principal stress post-processing
//!
//! # Quick start
//!
//! ```no_run
//! use oxiphysics_fem::thermal_stress::{compute_thermal_strain, plane_stress_thermal};
//!
//! let eps_th = compute_thermal_strain(50.0, 12e-6);
//! assert!((eps_th - 6e-4).abs() < 1e-12);
//!
//! let sig = plane_stress_thermal(200e9, 0.3, 12e-6, 50.0);
//! assert!(sig[0] < 0.0); // compressive on heating
//! ```

use crate::analysis::FemError;

// ─────────────────────────────────────────────────────────────────────────────
// § 1  FREE FUNCTIONS
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the thermal (free) strain for a given temperature change.
///
/// `ε_th = α · ΔT`
///
/// # Arguments
/// * `delta_t` – temperature change from the reference state \[K or °C\]
/// * `alpha`   – linear thermal expansion coefficient \[1/K\]
///
/// # Returns
/// Thermal strain (dimensionless).
pub fn compute_thermal_strain(delta_t: f64, alpha: f64) -> f64 {
    alpha * delta_t
}

/// Compute the plane-stress thermal stress vector `[σ_xx, σ_yy, τ_xy]`.
///
/// Under plane-stress (σ_zz = 0) the thermal stress is:
///
/// ```text
/// σ_xx = σ_yy = -E α ΔT / (1 − ν)
/// τ_xy = 0
/// ```
///
/// A *positive* ΔT (heating) produces *compressive* stresses (negative values).
///
/// # Arguments
/// * `e`     – Young's modulus \[Pa\]
/// * `nu`    – Poisson's ratio (dimensionless, 0 < ν < 0.5)
/// * `alpha` – linear thermal expansion coefficient \[1/K\]
/// * `d_t`   – temperature change ΔT \[K\]
///
/// # Returns
/// `[σ_xx, σ_yy, τ_xy]` in Pa.
pub fn plane_stress_thermal(e: f64, nu: f64, alpha: f64, d_t: f64) -> [f64; 3] {
    let s = -e * alpha * d_t / (1.0 - nu);
    [s, s, 0.0]
}

/// Compute the plane-strain thermal stress vector `[σ_xx, σ_yy, τ_xy]`.
///
/// Under plane-strain (ε_zz = 0) the thermal stress is:
///
/// ```text
/// σ_xx = σ_yy = -E α ΔT (1 + ν) / (1 − 2ν)  × 1/(1+ν)
///             = -E α ΔT / (1 − 2ν)
/// τ_xy = 0
/// ```
///
/// # Arguments
/// * `e`     – Young's modulus \[Pa\]
/// * `nu`    – Poisson's ratio
/// * `alpha` – linear thermal expansion coefficient \[1/K\]
/// * `d_t`   – temperature change ΔT \[K\]
///
/// # Returns
/// `[σ_xx, σ_yy, τ_xy]` in Pa.
pub fn plane_strain_thermal(e: f64, nu: f64, alpha: f64, d_t: f64) -> [f64; 3] {
    let s = -e * alpha * d_t / (1.0 - 2.0 * nu);
    [s, s, 0.0]
}

/// Compute the 3-D isotropic thermal stress `σ_th`.
///
/// For a fully constrained isotropic body:
/// `σ_th = -E α ΔT / (1 − 2ν)`
pub fn thermal_stress_3d(e: f64, nu: f64, alpha: f64, d_t: f64) -> f64 {
    -e * alpha * d_t / (1.0 - 2.0 * nu)
}

/// Compute the von Mises stress from the 2-D stress state.
///
/// `σ_vm = sqrt(σ_xx² − σ_xx σ_yy + σ_yy² + 3 τ_xy²)`
pub fn von_mises_2d(sxx: f64, syy: f64, txy: f64) -> f64 {
    (sxx * sxx - sxx * syy + syy * syy + 3.0 * txy * txy).sqrt()
}

/// Compute the principal stresses in 2-D.
///
/// Returns `[σ_1, σ_2, 0.0]` where σ_1 ≥ σ_2.
pub fn principal_stresses_2d(sxx: f64, syy: f64, txy: f64) -> [f64; 3] {
    let avg = 0.5 * (sxx + syy);
    let r = ((0.5 * (sxx - syy)).powi(2) + txy * txy).sqrt();
    [avg + r, avg - r, 0.0]
}

// ─────────────────────────────────────────────────────────────────────────────
// § 2  DATA STRUCTURES
// ─────────────────────────────────────────────────────────────────────────────

/// Thermal load applied to a thermoelastic analysis.
#[derive(Debug, Clone)]
pub struct ThermalLoad {
    /// Nodal temperature field (one entry per node) \[K or °C\].
    pub temperature_field: Vec<f64>,
    /// Reference (stress-free) temperature \[K or °C\].
    pub reference_temp: f64,
    /// Linear thermal expansion coefficient \[1/K\].
    pub alpha: f64,
}

impl ThermalLoad {
    /// Create a new `ThermalLoad`.
    ///
    /// # Arguments
    /// * `temperature_field` – nodal temperatures
    /// * `reference_temp`    – reference (stress-free) temperature
    /// * `alpha`             – thermal expansion coefficient
    pub fn new(temperature_field: Vec<f64>, reference_temp: f64, alpha: f64) -> Self {
        Self {
            temperature_field,
            reference_temp,
            alpha,
        }
    }

    /// Temperature change (ΔT) at node `i`.
    pub fn delta_t(&self, i: usize) -> f64 {
        self.temperature_field[i] - self.reference_temp
    }

    /// Thermal strain at node `i`.
    pub fn thermal_strain_at(&self, i: usize) -> f64 {
        compute_thermal_strain(self.delta_t(i), self.alpha)
    }

    /// Mean temperature change across all nodes.
    pub fn mean_delta_t(&self) -> f64 {
        if self.temperature_field.is_empty() {
            return 0.0;
        }
        let sum: f64 = self
            .temperature_field
            .iter()
            .map(|&t| t - self.reference_temp)
            .sum();
        sum / self.temperature_field.len() as f64
    }

    /// Maximum temperature change across all nodes.
    pub fn max_delta_t(&self) -> f64 {
        self.temperature_field
            .iter()
            .map(|&t| (t - self.reference_temp).abs())
            .fold(0.0_f64, f64::max)
    }
}

/// Coupling mode for thermoelastic analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CouplingMode {
    /// Thermal problem solved first; result used as load for mechanics (no feedback).
    OneWay,
    /// Fully coupled: mechanical deformation feeds back into thermal solution.
    ///
    /// The *steady-state* 1-D solver of [`ThermalStressAnalysis`] has no feedback
    /// path: the thermoelastic term of the heat equation is proportional to the
    /// strain rate and vanishes at steady state, so the temperature field cannot
    /// depend on the deformation.  [`ThermalStressAnalysis::solve`] therefore
    /// rejects this mode with [`FemError::InvalidInput`] instead of silently
    /// returning a one-way result.
    TwoWay,
}

/// Thermoelastic finite element for 1-D rod / 2-D element.
///
/// Carries node connectivity, Young's modulus, Poisson's ratio, and the
/// thermal expansion coefficient.  Call [`ThermoelasticElement::k_mechanical`]
/// for the elastic stiffness and [`ThermoelasticElement::k_thermal`] for the
/// thermal-conductivity sub-matrix.
#[derive(Debug, Clone)]
pub struct ThermoelasticElement {
    /// Global node IDs for this element (length = number of nodes).
    pub node_ids: Vec<usize>,
    /// Temperature DOF indices (one per node).
    pub temperature_dofs: Vec<usize>,
    /// Mechanical DOF indices (one per node × spatial dimension).
    pub mechanical_dofs: Vec<usize>,
    /// Young's modulus \[Pa\].
    pub e_modulus: f64,
    /// Poisson's ratio (dimensionless).
    pub nu: f64,
    /// Linear thermal expansion coefficient \[1/K\].
    pub alpha: f64,
    /// Thermal conductivity \[W/(m·K)\].
    pub conductivity: f64,
    /// Element length (for 1-D) or characteristic size \[m\].
    pub length: f64,
    /// Cross-sectional area (for 1-D rod) \[m²\].
    pub area: f64,
}

impl ThermoelasticElement {
    /// Create a 1-D thermoelastic rod element.
    ///
    /// Node IDs are `[n0, n1]`; temperature DOFs are `[n0, n1]`;
    /// mechanical DOFs are `[2*n0, 2*n0+1, 2*n1, 2*n1+1]` (u, v per node).
    pub fn new_1d_rod(
        n0: usize,
        n1: usize,
        e: f64,
        nu: f64,
        alpha: f64,
        k: f64,
        length: f64,
        area: f64,
    ) -> Self {
        Self {
            node_ids: vec![n0, n1],
            temperature_dofs: vec![n0, n1],
            mechanical_dofs: vec![2 * n0, 2 * n0 + 1, 2 * n1, 2 * n1 + 1],
            e_modulus: e,
            nu,
            alpha,
            conductivity: k,
            length,
            area,
        }
    }

    /// 1-D thermal conductance matrix (2×2) `K_th = (k A / L) [[1,-1\],[-1,1]]`.
    ///
    /// Returns a flat row-major 2×2 array: `[k00, k01, k10, k11]`.
    pub fn k_thermal(&self) -> [f64; 4] {
        let g = self.conductivity * self.area / self.length;
        [g, -g, -g, g]
    }

    /// 1-D mechanical (axial) stiffness matrix (2×2) `K_m = (E A / L) [[1,-1\],[-1,1]]`.
    ///
    /// Returns a flat row-major 2×2 array: `[k00, k01, k10, k11]`.
    pub fn k_mechanical(&self) -> [f64; 4] {
        let g = self.e_modulus * self.area / self.length;
        [g, -g, -g, g]
    }

    /// Thermal load vector for 1-D rod: `f_th = E A α ΔT_mean * [-1, 1]`.
    ///
    /// `delta_t_mean` is the mean temperature change along the element.
    pub fn f_thermal(&self, delta_t_mean: f64) -> [f64; 2] {
        let v = self.e_modulus * self.area * self.alpha * delta_t_mean;
        [-v, v]
    }

    /// Axial thermal stress `σ = -E α ΔT` (negative = compressive on heating).
    pub fn axial_thermal_stress(&self, delta_t: f64) -> f64 {
        -self.e_modulus * self.alpha * delta_t
    }
}

/// Stress result at a single point in a thermoelastic analysis.
#[derive(Debug, Clone)]
pub struct ThermalStress {
    /// Von Mises equivalent stress \[Pa\].
    pub von_mises: f64,
    /// Principal stresses `[σ_1, σ_2, σ_3]` \[Pa\].
    pub principal: [f64; 3],
    /// Temperature at the evaluation point \[K or °C\].
    pub temperature: f64,
    /// Direct stress components `[σ_xx, σ_yy, σ_zz]` \[Pa\].
    pub direct: [f64; 3],
    /// Shear stress components `[τ_xy, τ_yz, τ_zx]` \[Pa\].
    pub shear: [f64; 3],
}

impl ThermalStress {
    /// Construct a `ThermalStress` from 2-D plane-stress components.
    pub fn from_plane_stress(sxx: f64, syy: f64, txy: f64, temperature: f64) -> Self {
        let vm = von_mises_2d(sxx, syy, txy);
        let pr = principal_stresses_2d(sxx, syy, txy);
        Self {
            von_mises: vm,
            principal: pr,
            temperature,
            direct: [sxx, syy, 0.0],
            shear: [txy, 0.0, 0.0],
        }
    }

    /// Whether the stress state is compressive (von Mises dominated by compressive σ).
    pub fn is_compressive(&self) -> bool {
        self.direct[0] < 0.0 && self.direct[1] < 0.0
    }
}

/// Boundary condition type for thermoelastic analysis.
#[derive(Debug, Clone)]
pub struct BoundaryCondition {
    /// DOF index.
    pub dof: usize,
    /// Prescribed value (temperature \[K\] or displacement \[m\]).
    pub value: f64,
}

impl BoundaryCondition {
    /// Create a new boundary condition.
    pub fn new(dof: usize, value: f64) -> Self {
        Self { dof, value }
    }
}

/// Simple 1-D mesh for thermoelastic analysis.
#[derive(Debug, Clone)]
pub struct ThermoelasticMesh {
    /// Node positions along the rod \[m\].
    pub node_positions: Vec<f64>,
    /// Element connectivity: each entry is `[n0, n1]`.
    pub connectivity: Vec<[usize; 2]>,
    /// Cross-sectional area (uniform) \[m²\].
    pub area: f64,
}

impl ThermoelasticMesh {
    /// Create a uniform 1-D mesh with `n_elems` elements spanning `[0, length]`.
    pub fn uniform_1d(n_elems: usize, length: f64, area: f64) -> Self {
        let n_nodes = n_elems + 1;
        let dx = length / n_elems as f64;
        let node_positions: Vec<f64> = (0..n_nodes).map(|i| i as f64 * dx).collect();
        let connectivity: Vec<[usize; 2]> = (0..n_elems).map(|i| [i, i + 1]).collect();
        Self {
            node_positions,
            connectivity,
            area,
        }
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.node_positions.len()
    }

    /// Number of elements.
    pub fn n_elements(&self) -> usize {
        self.connectivity.len()
    }

    /// Element length for element `i`.
    pub fn element_length(&self, i: usize) -> f64 {
        let [n0, n1] = self.connectivity[i];
        (self.node_positions[n1] - self.node_positions[n0]).abs()
    }
}

/// Result of a thermoelastic solve.
#[derive(Debug, Clone)]
pub struct ThermoelasticResult {
    /// Nodal temperatures \[K or °C\].
    pub temperatures: Vec<f64>,
    /// Nodal displacements \[m\].
    pub displacements: Vec<f64>,
    /// Element-level thermal stresses \[Pa\].
    pub stresses: Vec<ThermalStress>,
}

/// High-level thermoelastic analysis driver.
///
/// Performs a one-way (thermal → mechanical) analysis on a 1-D rod mesh.  The
/// steady heat-conduction problem is solved first with a direct method; the
/// resulting temperatures, measured against the stress-free
/// [`reference_temp`](Self::reference_temp), load the axial mechanical problem
/// through the thermal strain `α (T − T_ref)`.  Element stresses are recovered
/// from the *total* strain of the displacement solution,
/// `σ = E (du/dx − α ΔT)`, so a freely expanding rod carries no stress while a
/// fully constrained one carries `−E α ΔT`.
#[derive(Debug, Clone)]
pub struct ThermalStressAnalysis {
    /// Discretised mesh.
    pub mesh: ThermoelasticMesh,
    /// Thermal boundary conditions (prescribed temperatures).
    pub thermal_bc: Vec<BoundaryCondition>,
    /// Mechanical boundary conditions (prescribed displacements).
    pub mechanical_bc: Vec<BoundaryCondition>,
    /// Young's modulus \[Pa\].
    pub e_modulus: f64,
    /// Poisson's ratio.
    pub nu: f64,
    /// Thermal expansion coefficient \[1/K\].
    pub alpha: f64,
    /// Thermal conductivity \[W/(m·K)\].
    pub conductivity: f64,
    /// Uniform heat generation rate \[W/m³\] (zero for pure-conduction).
    pub heat_source: f64,
    /// Stress-free reference temperature \[K or °C\]: the thermal strain is
    /// `α (T − reference_temp)`.
    pub reference_temp: f64,
    /// Coupling mode.  Only [`CouplingMode::OneWay`] is supported by the
    /// steady-state solver (see [`CouplingMode::TwoWay`]).
    pub coupling_mode: CouplingMode,
}

impl ThermalStressAnalysis {
    /// Create a `ThermalStressAnalysis` with default one-way coupling.
    ///
    /// `reference_temp` is the stress-free reference temperature; the thermal
    /// strain of an element is `α (T_mean − reference_temp)`.
    pub fn new(
        mesh: ThermoelasticMesh,
        thermal_bc: Vec<BoundaryCondition>,
        mechanical_bc: Vec<BoundaryCondition>,
        e_modulus: f64,
        nu: f64,
        alpha: f64,
        conductivity: f64,
        heat_source: f64,
        reference_temp: f64,
    ) -> Self {
        Self {
            mesh,
            thermal_bc,
            mechanical_bc,
            e_modulus,
            nu,
            alpha,
            conductivity,
            heat_source,
            reference_temp,
            coupling_mode: CouplingMode::OneWay,
        }
    }

    /// Solve the thermoelastic problem.
    ///
    /// Returns a [`ThermoelasticResult`] containing nodal temperatures,
    /// nodal displacements, and element stresses.
    ///
    /// # Errors
    /// [`FemError::InvalidInput`] when
    /// * `coupling_mode` is [`CouplingMode::TwoWay`] (no feedback path exists in
    ///   the steady-state solver),
    /// * the mesh or the material data are unusable (fewer than one element,
    ///   non-positive or non-finite area / modulus / conductivity / element
    ///   length, connectivity outside the node list, non-finite parameters),
    /// * a boundary condition addresses a node outside the mesh, has a
    ///   non-finite value, or conflicts with another one on the same node, or
    /// * the thermal or the mechanical problem has no boundary condition at all
    ///   (the steady problem would be singular).
    ///
    /// [`FemError::SingularSystem`] if a direct solve meets a vanishing pivot.
    pub fn solve(&self) -> Result<ThermoelasticResult, FemError> {
        self.validate()?;
        let temperatures = self.solve_thermal()?;
        let displacements = self.solve_mechanical(&temperatures)?;
        let stresses = self.compute_stresses(&temperatures, &displacements);
        Ok(ThermoelasticResult {
            temperatures,
            displacements,
            stresses,
        })
    }

    // ── Internal: input validation ─────────────────────────────────────────────

    fn validate(&self) -> Result<(), FemError> {
        let invalid = |message: String| Err(FemError::InvalidInput(message));
        if self.coupling_mode == CouplingMode::TwoWay {
            return invalid(
                "CouplingMode::TwoWay is not supported: the steady-state 1-D solver has no \
                 feedback path from the deformation to the temperature field"
                    .to_string(),
            );
        }
        let n_nodes = self.mesh.n_nodes();
        if self.mesh.n_elements() == 0 {
            return invalid("the mesh must contain at least one element".to_string());
        }
        for (name, value) in [
            ("cross-sectional area", self.mesh.area),
            ("Young's modulus", self.e_modulus),
            ("thermal conductivity", self.conductivity),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return invalid(format!("{name} must be positive and finite, got {value}"));
            }
        }
        for (name, value) in [
            ("thermal expansion coefficient", self.alpha),
            ("heat source", self.heat_source),
            ("reference temperature", self.reference_temp),
        ] {
            if !value.is_finite() {
                return invalid(format!("{name} must be finite, got {value}"));
            }
        }
        for (i, &[n0, n1]) in self.mesh.connectivity.iter().enumerate() {
            if n0 >= n_nodes || n1 >= n_nodes {
                return invalid(format!(
                    "element {i} connects nodes ({n0}, {n1}) but the mesh has {n_nodes} nodes"
                ));
            }
            let len = self.mesh.element_length(i);
            if !len.is_finite() || len <= 0.0 {
                return invalid(format!(
                    "element {i} must have a positive finite length, got {len}"
                ));
            }
        }
        Ok(())
    }

    // ── Internal: thermal solve ────────────────────────────────────────────────

    fn solve_thermal(&self) -> Result<Vec<f64>, FemError> {
        let n = self.mesh.n_nodes();
        // Assemble global conductivity matrix K (diagonal + off-diagonal).
        let mut k_global = vec![vec![0.0_f64; n]; n];
        let mut f_global = vec![0.0_f64; n];

        for i in 0..self.mesh.n_elements() {
            let [n0, n1] = self.mesh.connectivity[i];
            let len = self.mesh.element_length(i);
            let g = self.conductivity * self.mesh.area / len;
            k_global[n0][n0] += g;
            k_global[n0][n1] -= g;
            k_global[n1][n0] -= g;
            k_global[n1][n1] += g;
            // Distributed heat source (lumped equally to both nodes).
            let q = self.heat_source * self.mesh.area * len * 0.5;
            f_global[n0] += q;
            f_global[n1] += q;
        }

        Self::apply_dirichlet(&mut k_global, &mut f_global, &self.thermal_bc, "thermal")?;

        // Solve K T = f using Gaussian elimination (small systems only).
        Self::gauss_solve(&k_global, &f_global)
    }

    // ── Internal: mechanical solve ─────────────────────────────────────────────

    fn solve_mechanical(&self, temperatures: &[f64]) -> Result<Vec<f64>, FemError> {
        let n = self.mesh.n_nodes();
        let mut k_global = vec![vec![0.0_f64; n]; n];
        let mut f_global = vec![0.0_f64; n];

        for i in 0..self.mesh.n_elements() {
            let [n0, n1] = self.mesh.connectivity[i];
            let dx = self.mesh.node_positions[n1] - self.mesh.node_positions[n0];
            let len = dx.abs();
            let g = self.e_modulus * self.mesh.area / len;
            k_global[n0][n0] += g;
            k_global[n0][n1] -= g;
            k_global[n1][n0] -= g;
            k_global[n1][n1] += g;
            // Thermal load  f_th = ∫ Bᵀ E α ΔT A dx = E A α ΔT · sign(dx) · [-1, +1],
            // with ΔT measured from the stress-free reference temperature.
            let d_t_mean = 0.5 * (temperatures[n0] + temperatures[n1]) - self.reference_temp;
            let fth = dx.signum() * self.e_modulus * self.mesh.area * self.alpha * d_t_mean;
            f_global[n0] -= fth;
            f_global[n1] += fth;
        }

        // Dirichlet BCs for mechanics.
        Self::apply_dirichlet(
            &mut k_global,
            &mut f_global,
            &self.mechanical_bc,
            "mechanical",
        )?;

        Self::gauss_solve(&k_global, &f_global)
    }

    // ── Internal: stress recovery ──────────────────────────────────────────────

    /// Axial stress of every element from the total strain of the displacement
    /// field: `σ = E (du/dx − α (T_mean − T_ref))`.
    fn compute_stresses(&self, temperatures: &[f64], displacements: &[f64]) -> Vec<ThermalStress> {
        let mut stresses = Vec::with_capacity(self.mesh.n_elements());
        for i in 0..self.mesh.n_elements() {
            let [n0, n1] = self.mesh.connectivity[i];
            let dx = self.mesh.node_positions[n1] - self.mesh.node_positions[n0];
            let t_mean = 0.5 * (temperatures[n0] + temperatures[n1]);
            let thermal_strain = compute_thermal_strain(t_mean - self.reference_temp, self.alpha);
            let total_strain = (displacements[n1] - displacements[n0]) / dx;
            let sxx = self.e_modulus * (total_strain - thermal_strain);
            stresses.push(ThermalStress::from_plane_stress(sxx, 0.0, 0.0, t_mean));
        }
        stresses
    }

    // ── Internal: Dirichlet boundary conditions ────────────────────────────────

    /// Impose prescribed values `bcs` on the system `K x = f` by exact
    /// elimination: the known column is moved to the right-hand side of every
    /// free equation and the constrained row and column become the identity.
    ///
    /// Fails with [`FemError::InvalidInput`] for an out-of-range node, a
    /// non-finite value, two different values on one node, or an empty list
    /// (the steady problem is singular without a boundary condition).
    fn apply_dirichlet(
        k: &mut [Vec<f64>],
        f: &mut [f64],
        bcs: &[BoundaryCondition],
        kind: &str,
    ) -> Result<(), FemError> {
        let n = f.len();
        let mut prescribed: Vec<Option<f64>> = vec![None; n];
        for bc in bcs {
            if bc.dof >= n {
                return Err(FemError::InvalidInput(format!(
                    "{kind} boundary condition addresses node {} but the mesh has {n} nodes",
                    bc.dof
                )));
            }
            if !bc.value.is_finite() {
                return Err(FemError::InvalidInput(format!(
                    "{kind} boundary condition on node {} has the non-finite value {}",
                    bc.dof, bc.value
                )));
            }
            match prescribed[bc.dof] {
                Some(existing) if existing != bc.value => {
                    return Err(FemError::InvalidInput(format!(
                        "conflicting {kind} boundary conditions on node {} ({existing} and {})",
                        bc.dof, bc.value
                    )));
                }
                _ => prescribed[bc.dof] = Some(bc.value),
            }
        }
        if prescribed.iter().all(Option::is_none) {
            return Err(FemError::InvalidInput(format!(
                "at least one {kind} boundary condition is required: the steady problem is \
                 singular without it"
            )));
        }
        for (d, value) in prescribed.iter().enumerate() {
            if let Some(v) = *value {
                for (i, f_i) in f.iter_mut().enumerate() {
                    if prescribed[i].is_none() {
                        *f_i -= k[i][d] * v;
                    }
                }
            }
        }
        for (d, value) in prescribed.iter().enumerate() {
            if let Some(v) = *value {
                for j in 0..n {
                    k[d][j] = 0.0;
                    k[j][d] = 0.0;
                }
                k[d][d] = 1.0;
                f[d] = v;
            }
        }
        Ok(())
    }

    // ── Internal: Gaussian elimination (dense, n×n) ────────────────────────────

    fn gauss_solve(a: &[Vec<f64>], b: &[f64]) -> Result<Vec<f64>, FemError> {
        let n = b.len();
        let mut mat: Vec<Vec<f64>> = a.to_vec();
        let mut rhs: Vec<f64> = b.to_vec();

        for col in 0..n {
            // Find pivot.
            let mut pivot_row = col;
            let mut max_val = mat[col][col].abs();
            for (row, mat_row) in mat.iter().enumerate().skip(col + 1) {
                if mat_row[col].abs() > max_val {
                    max_val = mat_row[col].abs();
                    pivot_row = row;
                }
            }
            mat.swap(col, pivot_row);
            rhs.swap(col, pivot_row);

            let diag = mat[col][col];
            if !diag.is_finite() || diag.abs() < 1e-300 {
                return Err(FemError::SingularSystem);
            }
            let col_slice: Vec<f64> = mat[col][col..].to_vec();
            for row in (col + 1)..n {
                let factor = mat[row][col] / diag;
                rhs[row] -= factor * rhs[col];
                for (off, &cv) in col_slice.iter().enumerate() {
                    mat[row][col + off] -= factor * cv;
                }
            }
        }

        // Back substitution.
        let mut x = vec![0.0_f64; n];
        for row in (0..n).rev() {
            let mut sum = rhs[row];
            for c in (row + 1)..n {
                sum -= mat[row][c] * x[c];
            }
            x[row] = sum / mat[row][row];
        }
        Ok(x)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// § 3  UNIT TESTS
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod regression_tests;

#[cfg(test)]
mod tests {
    use super::*;

    // ── § 3.1  Free functions ──────────────────────────────────────────────────

    #[test]
    fn thermal_strain_zero_at_zero_dt() {
        assert_eq!(compute_thermal_strain(0.0, 12e-6), 0.0);
    }

    #[test]
    fn thermal_strain_proportional_to_alpha() {
        let eps1 = compute_thermal_strain(100.0, 12e-6);
        let eps2 = compute_thermal_strain(100.0, 24e-6);
        assert!((eps2 - 2.0 * eps1).abs() < 1e-15);
    }

    #[test]
    fn thermal_strain_proportional_to_delta_t() {
        let eps1 = compute_thermal_strain(50.0, 12e-6);
        let eps2 = compute_thermal_strain(100.0, 12e-6);
        assert!((eps2 - 2.0 * eps1).abs() < 1e-15);
    }

    #[test]
    fn thermal_strain_steel_50k() {
        // α = 12e-6, ΔT = 50 → ε = 6e-4
        let eps = compute_thermal_strain(50.0, 12e-6);
        assert!((eps - 6e-4).abs() < 1e-12);
    }

    #[test]
    fn plane_stress_thermal_compressive_on_heating() {
        let sig = plane_stress_thermal(200e9, 0.3, 12e-6, 100.0);
        assert!(sig[0] < 0.0, "σ_xx should be compressive");
        assert!(sig[1] < 0.0, "σ_yy should be compressive");
        assert_eq!(sig[2], 0.0, "τ_xy should be zero");
    }

    #[test]
    fn plane_stress_thermal_tensile_on_cooling() {
        let sig = plane_stress_thermal(200e9, 0.3, 12e-6, -50.0);
        assert!(sig[0] > 0.0, "σ_xx should be tensile on cooling");
    }

    #[test]
    fn plane_stress_thermal_zero_dt() {
        let sig = plane_stress_thermal(200e9, 0.3, 12e-6, 0.0);
        assert_eq!(sig, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn plane_stress_thermal_symmetry() {
        let sig = plane_stress_thermal(200e9, 0.3, 12e-6, 100.0);
        assert!(
            (sig[0] - sig[1]).abs() < 1e-6,
            "σ_xx == σ_yy in biaxial case"
        );
    }

    #[test]
    fn plane_strain_thermal_magnitude_greater_than_plane_stress() {
        // plane-strain denominator (1-2ν) < plane-stress denominator (1-ν)
        // so |σ_plane_strain| > |σ_plane_stress|
        let sig_ps = plane_stress_thermal(200e9, 0.3, 12e-6, 100.0);
        let sig_pe = plane_strain_thermal(200e9, 0.3, 12e-6, 100.0);
        assert!(sig_pe[0].abs() > sig_ps[0].abs());
    }

    #[test]
    fn plane_strain_thermal_zero_shear() {
        let sig = plane_strain_thermal(200e9, 0.3, 12e-6, 100.0);
        assert_eq!(sig[2], 0.0);
    }

    #[test]
    fn thermal_stress_3d_sign_positive_on_cooling() {
        let s = thermal_stress_3d(200e9, 0.3, 12e-6, -100.0);
        assert!(s > 0.0);
    }

    #[test]
    fn von_mises_biaxial_equal_principal() {
        // σ_xx = σ_yy = σ, τ = 0  →  σ_vm = σ (equal biaxial case: √(σ²-σ²+σ²) = σ)
        let vm = von_mises_2d(100.0, 100.0, 0.0);
        assert!((vm - 100.0).abs() < 1e-8);
    }

    #[test]
    fn von_mises_uniaxial() {
        // σ_xx = σ, σ_yy = τ = 0  →  σ_vm = σ
        let vm = von_mises_2d(200.0, 0.0, 0.0);
        assert!((vm - 200.0).abs() < 1e-8);
    }

    #[test]
    fn von_mises_pure_shear() {
        // σ_xx = σ_yy = 0, τ = τ  →  σ_vm = √3 τ
        let tau = 100.0_f64;
        let vm = von_mises_2d(0.0, 0.0, tau);
        assert!((vm - 3.0_f64.sqrt() * tau).abs() < 1e-8);
    }

    #[test]
    fn principal_stresses_uniaxial() {
        let [s1, s2, _] = principal_stresses_2d(100.0, 0.0, 0.0);
        assert!((s1 - 100.0).abs() < 1e-8);
        assert!((s2).abs() < 1e-8);
    }

    #[test]
    fn principal_stresses_biaxial_equal() {
        let [s1, s2, _] = principal_stresses_2d(100.0, 100.0, 0.0);
        assert!((s1 - 100.0).abs() < 1e-8);
        assert!((s2 - 100.0).abs() < 1e-8);
    }

    #[test]
    fn principal_stresses_pure_shear() {
        // σ_xx = σ_yy = 0, τ = τ  →  σ_1 = τ, σ_2 = -τ
        let tau = 50.0_f64;
        let [s1, s2, _] = principal_stresses_2d(0.0, 0.0, tau);
        assert!((s1 - tau).abs() < 1e-8);
        assert!((s2 + tau).abs() < 1e-8);
    }

    // ── § 3.2  ThermalLoad ─────────────────────────────────────────────────────

    #[test]
    fn thermal_load_delta_t_at_node() {
        let load = ThermalLoad::new(vec![350.0, 400.0, 450.0], 300.0, 12e-6);
        assert!((load.delta_t(0) - 50.0).abs() < 1e-12);
        assert!((load.delta_t(1) - 100.0).abs() < 1e-12);
        assert!((load.delta_t(2) - 150.0).abs() < 1e-12);
    }

    #[test]
    fn thermal_load_strain_at_node() {
        let load = ThermalLoad::new(vec![350.0], 300.0, 12e-6);
        let eps = load.thermal_strain_at(0);
        assert!((eps - 6e-4).abs() < 1e-12);
    }

    #[test]
    fn thermal_load_mean_delta_t() {
        let load = ThermalLoad::new(vec![310.0, 320.0, 330.0], 300.0, 12e-6);
        // ΔT values: 10, 20, 30 → mean = 20
        assert!((load.mean_delta_t() - 20.0).abs() < 1e-10);
    }

    #[test]
    fn thermal_load_mean_delta_t_empty() {
        let load = ThermalLoad::new(vec![], 300.0, 12e-6);
        assert_eq!(load.mean_delta_t(), 0.0);
    }

    #[test]
    fn thermal_load_max_delta_t() {
        let load = ThermalLoad::new(vec![350.0, 280.0, 400.0], 300.0, 12e-6);
        // |ΔT|: 50, 20, 100  → max = 100
        assert!((load.max_delta_t() - 100.0).abs() < 1e-10);
    }

    // ── § 3.3  ThermoelasticElement ────────────────────────────────────────────

    #[test]
    fn thermoelastic_element_k_thermal_symmetry() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        let k = elem.k_thermal();
        // k[0][0] = k[1][1] and k[0][1] = k[1][0]
        assert!((k[0] - k[3]).abs() < 1e-10);
        assert!((k[1] - k[2]).abs() < 1e-10);
    }

    #[test]
    fn thermoelastic_element_k_thermal_row_sum_zero() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        let k = elem.k_thermal();
        assert!((k[0] + k[1]).abs() < 1e-10);
        assert!((k[2] + k[3]).abs() < 1e-10);
    }

    #[test]
    fn thermoelastic_element_k_mechanical_symmetry() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        let k = elem.k_mechanical();
        assert!((k[0] - k[3]).abs() < 1e-10);
        assert!((k[1] - k[2]).abs() < 1e-10);
    }

    #[test]
    fn thermoelastic_element_k_mechanical_row_sum_zero() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        let k = elem.k_mechanical();
        assert!((k[0] + k[1]).abs() < 1e-10);
        assert!((k[2] + k[3]).abs() < 1e-10);
    }

    #[test]
    fn thermoelastic_element_f_thermal_sign() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        let f = elem.f_thermal(50.0); // heating → node 0 pulls left, node 1 pushes right
        assert!(f[0] < 0.0);
        assert!(f[1] > 0.0);
    }

    #[test]
    fn thermoelastic_element_axial_stress_compressive_heating() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        let s = elem.axial_thermal_stress(100.0);
        assert!(s < 0.0);
    }

    #[test]
    fn thermoelastic_element_axial_stress_zero_dt() {
        let elem = ThermoelasticElement::new_1d_rod(0, 1, 200e9, 0.3, 12e-6, 50.0, 1.0, 0.01);
        assert_eq!(elem.axial_thermal_stress(0.0), 0.0);
    }

    // ── § 3.4  ThermalStress ───────────────────────────────────────────────────

    #[test]
    fn thermal_stress_from_plane_stress_compressive() {
        let ts = ThermalStress::from_plane_stress(-100.0e6, -100.0e6, 0.0, 400.0);
        assert!(ts.is_compressive());
        assert!((ts.von_mises - 100.0e6).abs() < 1.0);
    }

    #[test]
    fn thermal_stress_principal_order() {
        // σ_1 ≥ σ_2
        let ts = ThermalStress::from_plane_stress(150.0, -50.0, 0.0, 300.0);
        assert!(ts.principal[0] >= ts.principal[1]);
    }

    #[test]
    fn thermal_stress_temperature_stored() {
        let ts = ThermalStress::from_plane_stress(0.0, 0.0, 0.0, 500.0);
        assert_eq!(ts.temperature, 500.0);
    }

    // ── § 3.5  ThermoelasticMesh ───────────────────────────────────────────────

    #[test]
    fn mesh_uniform_1d_node_count() {
        let mesh = ThermoelasticMesh::uniform_1d(4, 1.0, 0.01);
        assert_eq!(mesh.n_nodes(), 5);
    }

    #[test]
    fn mesh_uniform_1d_element_count() {
        let mesh = ThermoelasticMesh::uniform_1d(4, 1.0, 0.01);
        assert_eq!(mesh.n_elements(), 4);
    }

    #[test]
    fn mesh_uniform_1d_element_length() {
        let mesh = ThermoelasticMesh::uniform_1d(4, 1.0, 0.01);
        for i in 0..4 {
            assert!((mesh.element_length(i) - 0.25).abs() < 1e-12);
        }
    }

    #[test]
    fn mesh_node_positions_start_end() {
        let mesh = ThermoelasticMesh::uniform_1d(5, 2.0, 0.01);
        assert!((mesh.node_positions[0]).abs() < 1e-12);
        assert!((mesh.node_positions[5] - 2.0).abs() < 1e-12);
    }

    // ── § 3.6  ThermalStressAnalysis (integration) ────────────────────────────

    #[test]
    fn analysis_linear_temperature_no_panic() {
        let mesh = ThermoelasticMesh::uniform_1d(5, 1.0, 1e-4);
        let thermal_bc = vec![
            BoundaryCondition::new(0, 0.0),   // T(0) = 0
            BoundaryCondition::new(5, 100.0), // T(L) = 100
        ];
        let mech_bc = vec![BoundaryCondition::new(0, 0.0)]; // fixed at x=0
        let analysis = ThermalStressAnalysis::new(
            mesh, thermal_bc, mech_bc, 200e9, 0.3, 12e-6, 50.0, 0.0, 0.0,
        );
        let result = analysis.solve().expect("well-posed analysis");
        assert_eq!(result.temperatures.len(), 6);
        assert_eq!(result.displacements.len(), 6);
        assert_eq!(result.stresses.len(), 5);
    }

    #[test]
    fn analysis_uniform_temperature_gives_uniform_field() {
        let mesh = ThermoelasticMesh::uniform_1d(3, 1.0, 1e-4);
        let thermal_bc = vec![
            BoundaryCondition::new(0, 50.0),
            BoundaryCondition::new(3, 50.0),
        ];
        let mech_bc = vec![BoundaryCondition::new(0, 0.0)];
        let analysis = ThermalStressAnalysis::new(
            mesh, thermal_bc, mech_bc, 200e9, 0.3, 12e-6, 50.0, 0.0, 0.0,
        );
        let result = analysis.solve().expect("well-posed analysis");
        // All temperatures should be ~50 °C
        for &t in &result.temperatures {
            assert!((t - 50.0).abs() < 1.0);
        }
    }

    #[test]
    fn analysis_coupling_mode_default_one_way() {
        let mesh = ThermoelasticMesh::uniform_1d(2, 1.0, 1e-4);
        let analysis = ThermalStressAnalysis::new(
            mesh,
            vec![
                BoundaryCondition::new(0, 0.0),
                BoundaryCondition::new(2, 100.0),
            ],
            vec![BoundaryCondition::new(0, 0.0)],
            200e9,
            0.3,
            12e-6,
            50.0,
            0.0,
            0.0,
        );
        assert_eq!(analysis.coupling_mode, CouplingMode::OneWay);
    }

    #[test]
    fn coupling_mode_equality() {
        assert_eq!(CouplingMode::OneWay, CouplingMode::OneWay);
        assert_ne!(CouplingMode::OneWay, CouplingMode::TwoWay);
    }
}
