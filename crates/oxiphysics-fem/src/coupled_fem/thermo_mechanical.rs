// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Coupled thermo-elastic analysis of a rectangular plate (plane stress).
//!
//! The plate `[0, lx] × [0, ly]` (unit thickness) is discretised by a structured
//! grid of bilinear Q4 elements shared by both fields:
//!
//! * **Heat**: ρc_p Ṫ = ∇·(κ∇T) + Q − T₀ β tr ε̇ is integrated with backward
//!   Euler (unconditionally stable, lumped capacity).  `β = E(T)α/(1−ν)` is the
//!   plane-stress thermo-elastic coupling modulus and `T₀` the absolute
//!   reference temperature.
//! * **Mechanics**: K(T) u = f_th + f_ext with the plane-stress constitutive
//!   matrix `D(T)` and the thermal load `f_th = ∫ Bᵀ D α ΔT {1, 1, 0}ᵀ dA`.
//!
//! A time step couples the two fields with a staggered (block Gauss-Seidel)
//! iteration until the displacement update falls below `tol`.

use super::types::{ThermalField, ThermoElasticMaterial};
use crate::error::{Error, Result};
use crate::solvers::pcg_solve;
use crate::sparse::CsrMatrix;

/// Natural coordinates of the Q4 corner nodes (counter-clockwise).
const Q4_CORNERS: [[f64; 2]; 4] = [[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]];
/// 2-point Gauss abscissa (weights are 1).
const GAUSS_2: [f64; 2] = [-0.577_350_269_189_625_8, 0.577_350_269_189_625_8];
/// Relative residual tolerance of the PCG solves.
const PCG_RTOL: f64 = 1e-12;

/// (current displacement iterate, displacements at the start of the step).
type DisplacementPair<'a> = (&'a [[f64; 2]], &'a [[f64; 2]]);

/// Convergence report of one staggered thermo-mechanical time step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StaggeredReport {
    /// Number of heat/mechanics passes performed.
    pub iterations: usize,
    /// Relative displacement change ‖uᵏ⁺¹ − uᵏ‖ / ‖uᵏ⁺¹‖ of the last pass.
    pub displacement_change: f64,
    /// Whether the change dropped below the tolerance.
    pub converged: bool,
}

/// Thermo-mechanical staggered coupling driver (plane stress, Q4 elements).
#[derive(Debug, Clone)]
pub struct ThermoMechanicalCoupling {
    /// Nodal temperature field \[K\] on the `nx × ny` grid.
    pub thermal: ThermalField,
    /// Nodal displacement field \[u_x, u_y\] \[m\].
    pub displacements: Vec<[f64; 2]>,
    /// Material model.
    pub material: ThermoElasticMaterial,
    /// Relative displacement tolerance of the staggered iteration.
    pub tol: f64,
    /// Maximum staggered iterations per time step.
    pub max_iter: usize,
    /// Plate length in x \[m\].
    pub lx: f64,
    /// Plate length in y \[m\].
    pub ly: f64,
    /// Constrained displacement DOFs (`2·node + component`, held at zero).
    ///
    /// Defaults to the minimal statically determinate support: node (0, 0)
    /// fixed in x and y, node (nx−1, 0) fixed in y.  This removes the rigid-body
    /// modes without restraining free thermal expansion.
    pub fixed_dofs: Vec<usize>,
    /// Prescribed nodal temperatures `(node, T)`; all other boundary nodes are
    /// insulated (zero heat flux).
    pub thermal_dirichlet: Vec<(usize, f64)>,
    /// External nodal forces \[N per unit thickness\].
    pub external_forces: Vec<[f64; 2]>,
    /// Include the thermo-elastic heat term −T₀β tr ε̇ (two-way coupling).
    pub thermoelastic_coupling: bool,
}

impl ThermoMechanicalCoupling {
    /// Create a coupled problem on the unit square with an `nx × ny` node grid.
    ///
    /// # Errors
    /// See [`Self::with_domain`].
    pub fn new(nx: usize, ny: usize, material: ThermoElasticMaterial) -> Result<Self> {
        Self::with_domain(nx, ny, 1.0, 1.0, material)
    }

    /// Create a coupled problem on the plate `[0, lx] × [0, ly]`.
    ///
    /// The initial temperature is the material reference temperature, so the
    /// plate starts stress free.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if `nx < 2`, `ny < 2`, the plate size is not
    /// positive, or the material constants are not physical.
    pub fn with_domain(
        nx: usize,
        ny: usize,
        lx: f64,
        ly: f64,
        material: ThermoElasticMaterial,
    ) -> Result<Self> {
        if nx < 2 || ny < 2 {
            return Err(Error::InvalidInput(format!(
                "thermo-mechanical grid needs at least 2×2 nodes, got {nx}×{ny}"
            )));
        }
        if !(lx.is_finite() && ly.is_finite() && lx > 0.0 && ly > 0.0) {
            return Err(Error::InvalidInput(format!(
                "plate dimensions must be positive, got {lx}×{ly}"
            )));
        }
        let m = &material;
        let physical = m.e0 > 0.0
            && m.nu > -1.0
            && m.nu < 0.5
            && m.rho > 0.0
            && m.cp > 0.0
            && m.kappa >= 0.0
            && m.alpha.is_finite()
            && m.t_ref.is_finite();
        if !physical {
            return Err(Error::InvalidInput(
                "thermo-elastic material constants are not physical".into(),
            ));
        }
        let n_nodes = nx * ny;
        let thermal = ThermalField::uniform(nx, ny, material.t_ref);
        let corner = thermal.idx(nx - 1, 0);
        Ok(Self {
            thermal,
            displacements: vec![[0.0; 2]; n_nodes],
            material,
            tol: 1e-8,
            max_iter: 50,
            lx,
            ly,
            fixed_dofs: vec![0, 1, 2 * corner + 1],
            thermal_dirichlet: Vec::new(),
            external_forces: vec![[0.0; 2]; n_nodes],
            thermoelastic_coupling: true,
        })
    }

    /// Grid spacing in x.
    pub fn hx(&self) -> f64 {
        self.lx / (self.thermal.nx - 1) as f64
    }

    /// Grid spacing in y.
    pub fn hy(&self) -> f64 {
        self.ly / (self.thermal.ny - 1) as f64
    }

    /// Number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.thermal.nx * self.thermal.ny
    }

    /// Prescribe temperature `t_wall` on the left edge (x = 0).
    pub fn set_left_temperature(&mut self, t_wall: f64) {
        for j in 0..self.thermal.ny {
            let node = self.thermal.idx(0, j);
            self.prescribe_temperature(node, t_wall);
        }
    }

    /// Prescribe temperature `t_wall` on the right edge (x = lx).
    pub fn set_right_temperature(&mut self, t_wall: f64) {
        let nx = self.thermal.nx;
        for j in 0..self.thermal.ny {
            let node = self.thermal.idx(nx - 1, j);
            self.prescribe_temperature(node, t_wall);
        }
    }

    fn prescribe_temperature(&mut self, node: usize, value: f64) {
        if let Some(entry) = self.thermal_dirichlet.iter_mut().find(|(n, _)| *n == node) {
            entry.1 = value;
        } else {
            self.thermal_dirichlet.push((node, value));
        }
        if let Some(t) = self.thermal.temperatures.get_mut(node) {
            *t = value;
        }
    }

    /// Constrain every displacement DOF on the plate boundary (fully clamped edges).
    pub fn clamp_all_edges(&mut self) {
        let (nx, ny) = (self.thermal.nx, self.thermal.ny);
        let mut dofs = Vec::new();
        for i in 0..nx {
            for j in 0..ny {
                if i == 0 || j == 0 || i == nx - 1 || j == ny - 1 {
                    let node = self.thermal.idx(i, j);
                    dofs.push(2 * node);
                    dofs.push(2 * node + 1);
                }
            }
        }
        self.fixed_dofs = dofs;
    }

    // -----------------------------------------------------------------------
    // Element kinematics
    // -----------------------------------------------------------------------

    /// Global node indices of element (ei, ej) in counter-clockwise order.
    fn element_nodes(&self, ei: usize, ej: usize) -> [usize; 4] {
        let t = &self.thermal;
        [
            t.idx(ei, ej),
            t.idx(ei + 1, ej),
            t.idx(ei + 1, ej + 1),
            t.idx(ei, ej + 1),
        ]
    }

    /// Shape functions and physical gradients at (ξ, η).
    fn shape(&self, xi: f64, eta: f64) -> ([f64; 4], [[f64; 2]; 4]) {
        let (hx, hy) = (self.hx(), self.hy());
        let mut n = [0.0; 4];
        let mut dn = [[0.0; 2]; 4];
        for (a, c) in Q4_CORNERS.iter().enumerate() {
            n[a] = 0.25 * (1.0 + c[0] * xi) * (1.0 + c[1] * eta);
            dn[a][0] = 0.25 * c[0] * (1.0 + c[1] * eta) * 2.0 / hx;
            dn[a][1] = 0.25 * c[1] * (1.0 + c[0] * xi) * 2.0 / hy;
        }
        (n, dn)
    }

    /// Strain [εxx, εyy, γxy] from element displacements at gradient `dn`.
    fn strain(&self, nodes: &[usize; 4], dn: &[[f64; 2]; 4], u: &[[f64; 2]]) -> [f64; 3] {
        let mut eps = [0.0; 3];
        for (a, &node) in nodes.iter().enumerate() {
            let [ux, uy] = u[node];
            eps[0] += dn[a][0] * ux;
            eps[1] += dn[a][1] * uy;
            eps[2] += dn[a][1] * ux + dn[a][0] * uy;
        }
        eps
    }

    fn interpolate_temperature(&self, nodes: &[usize; 4], n: &[f64; 4]) -> f64 {
        nodes
            .iter()
            .zip(n.iter())
            .map(|(&node, &na)| na * self.thermal.temperatures[node])
            .sum()
    }

    // -----------------------------------------------------------------------
    // Heat
    // -----------------------------------------------------------------------

    /// Advance the temperature field by one backward-Euler step of length `dt`
    /// with volumetric heat source `heat_source` \[W/m³\], without mechanical
    /// feedback (one-way coupling).  Call [`Self::solve_mechanics`] afterwards
    /// to update the displacements.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] for a non-positive `dt` or an invalid prescribed
    /// node; [`Error::General`] if the linear solve fails to converge.
    pub fn advance_heat(&mut self, dt: f64, heat_source: f64) -> Result<()> {
        let t_old = self.thermal.temperatures.clone();
        let t_new = self.heat_solve(&t_old, dt, heat_source, None)?;
        self.thermal.temperatures = t_new;
        Ok(())
    }

    /// Backward-Euler heat solve starting from `t_old`.  `coupling` holds the
    /// current displacement iterate and the displacements at the start of the
    /// step for the thermo-elastic term.
    fn heat_solve(
        &self,
        t_old: &[f64],
        dt: f64,
        heat_source: f64,
        coupling: Option<DisplacementPair<'_>>,
    ) -> Result<Vec<f64>> {
        if !(dt.is_finite() && dt > 0.0) {
            return Err(Error::InvalidInput(format!(
                "time step must be positive, got {dt}"
            )));
        }
        let n = self.n_nodes();
        let (nx, ny) = (self.thermal.nx, self.thermal.ny);
        let mat = &self.material;
        let det_j = 0.25 * self.hx() * self.hy();
        let rho_cp = mat.rho * mat.cp;
        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity((nx - 1) * (ny - 1) * 16);
        let mut capacity = vec![0.0; n];
        let mut rhs = vec![0.0; n];
        for ei in 0..nx - 1 {
            for ej in 0..ny - 1 {
                let nodes = self.element_nodes(ei, ej);
                let mut ke = [[0.0; 4]; 4];
                for &xi in &GAUSS_2 {
                    for &eta in &GAUSS_2 {
                        let (shape_n, dn) = self.shape(xi, eta);
                        for a in 0..4 {
                            for b in 0..4 {
                                ke[a][b] +=
                                    mat.kappa * (dn[a][0] * dn[b][0] + dn[a][1] * dn[b][1]) * det_j;
                            }
                            // Lumped capacity (row sum of the consistent matrix).
                            capacity[nodes[a]] += rho_cp * shape_n[a] * det_j;
                            rhs[nodes[a]] += shape_n[a] * heat_source * det_j;
                        }
                        if let Some((u_iter, u_start)) = coupling {
                            let t_gp = self.interpolate_temperature(&nodes, &shape_n);
                            let beta = mat.modulus_at(t_gp) * mat.alpha / (1.0 - mat.nu);
                            let e_new = self.strain(&nodes, &dn, u_iter);
                            let e_old = self.strain(&nodes, &dn, u_start);
                            let tr_rate = (e_new[0] + e_new[1] - e_old[0] - e_old[1]) / dt;
                            let q_c = mat.t_ref * beta * tr_rate;
                            for a in 0..4 {
                                rhs[nodes[a]] -= shape_n[a] * q_c * det_j;
                            }
                        }
                    }
                }
                for a in 0..4 {
                    for b in 0..4 {
                        triplets.push((nodes[a], nodes[b], ke[a][b]));
                    }
                }
            }
        }
        for i in 0..n {
            triplets.push((i, i, capacity[i] / dt));
            rhs[i] += capacity[i] / dt * t_old[i];
        }
        let mut prescribed: Vec<Option<f64>> = vec![None; n];
        for &(node, value) in &self.thermal_dirichlet {
            if node >= n {
                return Err(Error::DofOutOfRange { dof: node, n });
            }
            prescribed[node] = Some(value);
        }
        solve_with_dirichlet(n, &triplets, rhs, &prescribed)
    }

    // -----------------------------------------------------------------------
    // Mechanics
    // -----------------------------------------------------------------------

    /// Assemble and solve the plane-stress equilibrium K(T) u = f_th + f_ext for
    /// the current temperature field, replacing [`Self::displacements`].
    ///
    /// Returns the norm of the displacement change ‖u_new − u_old‖.
    ///
    /// # Errors
    /// [`Error::InvalidInput`] if no displacement DOF is constrained,
    /// [`Error::DofOutOfRange`] for an invalid constrained DOF and
    /// [`Error::General`] if the linear solve fails to converge.
    pub fn solve_mechanics(&mut self) -> Result<f64> {
        let u_new = self.mechanics_solve()?;
        let change = u_new
            .iter()
            .zip(self.displacements.iter())
            .map(|(a, b)| (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2))
            .sum::<f64>()
            .sqrt();
        self.displacements = u_new;
        Ok(change)
    }

    fn mechanics_solve(&self) -> Result<Vec<[f64; 2]>> {
        let n_nodes = self.n_nodes();
        let n_dof = 2 * n_nodes;
        if self.fixed_dofs.is_empty() {
            return Err(Error::InvalidInput(
                "no displacement constraints: the plate stiffness is singular (rigid-body modes)"
                    .into(),
            ));
        }
        let (nx, ny) = (self.thermal.nx, self.thermal.ny);
        let mat = &self.material;
        let det_j = 0.25 * self.hx() * self.hy();
        let mut triplets: Vec<(usize, usize, f64)> = Vec::with_capacity((nx - 1) * (ny - 1) * 64);
        let mut rhs = vec![0.0; n_dof];
        for (node, f) in self.external_forces.iter().enumerate().take(n_nodes) {
            rhs[2 * node] += f[0];
            rhs[2 * node + 1] += f[1];
        }
        for ei in 0..nx - 1 {
            for ej in 0..ny - 1 {
                let nodes = self.element_nodes(ei, ej);
                let mut ke = [[0.0; 8]; 8];
                let mut fe = [0.0; 8];
                for &xi in &GAUSS_2 {
                    for &eta in &GAUSS_2 {
                        let (shape_n, dn) = self.shape(xi, eta);
                        let t_gp = self.interpolate_temperature(&nodes, &shape_n);
                        let d = mat.constitutive_2d(t_gp);
                        let sigma_th = mat.thermal_load_vector_2d(t_gp, t_gp - mat.t_ref);
                        let b = b_matrix(&dn);
                        // DB (3×8)
                        let mut db = [[0.0; 8]; 3];
                        for r in 0..3 {
                            for c in 0..8 {
                                db[r][c] = (0..3).map(|k| d[r][k] * b[k][c]).sum();
                            }
                        }
                        for p in 0..8 {
                            for q in 0..8 {
                                ke[p][q] += (0..3).map(|k| b[k][p] * db[k][q]).sum::<f64>() * det_j;
                            }
                            fe[p] += (0..3).map(|k| b[k][p] * sigma_th[k]).sum::<f64>() * det_j;
                        }
                    }
                }
                let dofs = element_dofs(&nodes);
                for p in 0..8 {
                    rhs[dofs[p]] += fe[p];
                    for q in 0..8 {
                        triplets.push((dofs[p], dofs[q], ke[p][q]));
                    }
                }
            }
        }
        let mut prescribed: Vec<Option<f64>> = vec![None; n_dof];
        for &dof in &self.fixed_dofs {
            if dof >= n_dof {
                return Err(Error::DofOutOfRange { dof, n: n_dof });
            }
            prescribed[dof] = Some(0.0);
        }
        let u = solve_with_dirichlet(n_dof, &triplets, rhs, &prescribed)?;
        Ok((0..n_nodes).map(|i| [u[2 * i], u[2 * i + 1]]).collect())
    }

    // -----------------------------------------------------------------------
    // Coupled time stepping
    // -----------------------------------------------------------------------

    /// Advance the coupled problem by one time step `dt` with volumetric heat
    /// source `heat_source` \[W/m³\].
    ///
    /// Each staggered pass solves the heat equation (including the
    /// thermo-elastic term from the current displacement iterate when
    /// [`Self::thermoelastic_coupling`] is set) and then the mechanical
    /// equilibrium for the new temperatures, until the relative displacement
    /// update is below [`Self::tol`].  Without thermo-elastic coupling the
    /// mechanics does not feed back into the heat equation and a single pass is
    /// exact.
    ///
    /// The displacement field on entry is taken as the converged state at the
    /// start of the step; call [`Self::solve_mechanics`] first if the
    /// temperatures were modified externally.
    ///
    /// # Errors
    /// Propagates the errors of the heat and mechanics solves.
    pub fn staggered_step(&mut self, dt: f64, heat_source: f64) -> Result<StaggeredReport> {
        let t_start = self.thermal.temperatures.clone();
        let u_start = self.displacements.clone();
        let passes = if self.thermoelastic_coupling {
            self.max_iter.max(1)
        } else {
            1
        };
        let mut report = StaggeredReport {
            iterations: 0,
            displacement_change: f64::INFINITY,
            converged: false,
        };
        for pass in 0..passes {
            let coupling = if self.thermoelastic_coupling {
                Some((self.displacements.as_slice(), u_start.as_slice()))
            } else {
                None
            };
            let t_new = self.heat_solve(&t_start, dt, heat_source, coupling)?;
            self.thermal.temperatures = t_new;
            let change = self.solve_mechanics()?;
            let norm = self
                .displacements
                .iter()
                .map(|u| u[0] * u[0] + u[1] * u[1])
                .sum::<f64>()
                .sqrt();
            let rel = if norm > 0.0 { change / norm } else { change };
            report = StaggeredReport {
                iterations: pass + 1,
                displacement_change: rel,
                converged: rel < self.tol || !self.thermoelastic_coupling,
            };
            if report.converged {
                break;
            }
        }
        Ok(report)
    }

    /// Run `n_steps` coupled time steps; returns the per-step reports.
    ///
    /// # Errors
    /// Propagates the errors of [`Self::staggered_step`].
    pub fn solve(
        &mut self,
        dt: f64,
        heat_source: f64,
        n_steps: usize,
    ) -> Result<Vec<StaggeredReport>> {
        (0..n_steps)
            .map(|_| self.staggered_step(dt, heat_source))
            .collect()
    }

    // -----------------------------------------------------------------------
    // Post-processing
    // -----------------------------------------------------------------------

    /// Stress [σxx, σyy, τxy] in element (ei, ej) at natural coordinates (ξ, η):
    /// σ = D(T) (B u − α ΔT {1, 1, 0}).
    ///
    /// Returns `None` for an element index outside the grid.
    pub fn element_stress(&self, ei: usize, ej: usize, xi: f64, eta: f64) -> Option<[f64; 3]> {
        if ei + 1 >= self.thermal.nx || ej + 1 >= self.thermal.ny {
            return None;
        }
        let nodes = self.element_nodes(ei, ej);
        let (shape_n, dn) = self.shape(xi, eta);
        let t = self.interpolate_temperature(&nodes, &shape_n);
        let eps = self.strain(&nodes, &dn, &self.displacements);
        let eps_th = self.material.alpha * (t - self.material.t_ref);
        let mech = [eps[0] - eps_th, eps[1] - eps_th, eps[2]];
        let d = self.material.constitutive_2d(t);
        let mut sigma = [0.0; 3];
        for (r, s) in sigma.iter_mut().enumerate() {
            *s = (0..3).map(|k| d[r][k] * mech[k]).sum();
        }
        Some(sigma)
    }

    /// Nodal stress [σxx, σyy, τxy], averaged over the elements sharing the node
    /// (each evaluated at that corner).  Returns `None` for an invalid node.
    pub fn nodal_stress(&self, node: usize) -> Option<[f64; 3]> {
        if node >= self.n_nodes() {
            return None;
        }
        let ny = self.thermal.ny;
        let (i, j) = (node / ny, node % ny);
        let mut acc = [0.0; 3];
        let mut count = 0usize;
        for (di, dj) in [(0usize, 0usize), (1, 0), (0, 1), (1, 1)] {
            let (Some(ei), Some(ej)) = (i.checked_sub(di), j.checked_sub(dj)) else {
                continue;
            };
            // Corner of element (ei, ej) at which the node sits.
            let xi = if di == 0 { -1.0 } else { 1.0 };
            let eta = if dj == 0 { -1.0 } else { 1.0 };
            if let Some(s) = self.element_stress(ei, ej, xi, eta) {
                for (a, v) in acc.iter_mut().zip(s.iter()) {
                    *a += v;
                }
                count += 1;
            }
        }
        if count == 0 {
            return None;
        }
        Some(acc.map(|v| v / count as f64))
    }

    /// Von Mises stress (plane stress) at a node from the averaged nodal stress:
    /// σ_vm = √(σx² − σxσy + σy² + 3τ²).  Returns `None` for an invalid node.
    pub fn von_mises_stress(&self, node: usize) -> Option<f64> {
        self.nodal_stress(node).map(|[sx, sy, txy]| {
            (sx * sx - sx * sy + sy * sy + 3.0 * txy * txy)
                .max(0.0)
                .sqrt()
        })
    }
}

/// Strain-displacement matrix B (3 × 8) for engineering strains [εxx, εyy, γxy].
fn b_matrix(dn: &[[f64; 2]; 4]) -> [[f64; 8]; 3] {
    let mut b = [[0.0; 8]; 3];
    for (a, g) in dn.iter().enumerate() {
        b[0][2 * a] = g[0];
        b[1][2 * a + 1] = g[1];
        b[2][2 * a] = g[1];
        b[2][2 * a + 1] = g[0];
    }
    b
}

fn element_dofs(nodes: &[usize; 4]) -> [usize; 8] {
    let mut d = [0usize; 8];
    for (a, &n) in nodes.iter().enumerate() {
        d[2 * a] = 2 * n;
        d[2 * a + 1] = 2 * n + 1;
    }
    d
}

/// Solve the SPD system given in triplet form with Dirichlet values eliminated
/// symmetrically (known columns moved to the right-hand side, constrained rows
/// replaced by identity rows), using Jacobi-preconditioned CG.
fn solve_with_dirichlet(
    n: usize,
    triplets: &[(usize, usize, f64)],
    mut rhs: Vec<f64>,
    prescribed: &[Option<f64>],
) -> Result<Vec<f64>> {
    let mut reduced: Vec<(usize, usize, f64)> = Vec::with_capacity(triplets.len() + n);
    for &(r, c, v) in triplets {
        match (prescribed[r], prescribed[c]) {
            (None, None) => reduced.push((r, c, v)),
            (None, Some(g)) => rhs[r] -= v * g,
            _ => {}
        }
    }
    for (i, p) in prescribed.iter().enumerate() {
        if let Some(g) = p {
            reduced.push((i, i, 1.0));
            rhs[i] = *g;
        }
    }
    let a = CsrMatrix::from_triplets(n, n, &reduced);
    let mut precond = vec![0.0; n];
    for (i, m) in precond.iter_mut().enumerate() {
        let d = a.diagonal(i);
        if !(d.is_finite() && d > 0.0) {
            return Err(Error::SingularSystem);
        }
        *m = 1.0 / d;
    }
    let (x, stats) = pcg_solve(&a, &rhs, &precond, 20 * n + 100, PCG_RTOL);
    if !stats.converged || x.iter().any(|v| !v.is_finite()) {
        return Err(Error::General(format!(
            "thermo-mechanical PCG solve did not converge (residual {:e} after {} iterations)",
            stats.residual, stats.iterations
        )));
    }
    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn steel() -> ThermoElasticMaterial {
        ThermoElasticMaterial::new(200e9, 0.3, 12e-6, 50.0, 500.0, 7850.0)
    }

    #[test]
    fn rejects_degenerate_grid() {
        assert!(ThermoMechanicalCoupling::new(1, 5, steel()).is_err());
        assert!(ThermoMechanicalCoupling::new(5, 0, steel()).is_err());
        assert!(ThermoMechanicalCoupling::with_domain(3, 3, 0.0, 1.0, steel()).is_err());
    }

    #[test]
    fn free_plate_uniform_heating_is_stress_free() {
        let mat = steel();
        let mut c = ThermoMechanicalCoupling::with_domain(5, 4, 2.0, 1.0, mat.clone())
            .expect("valid problem");
        let delta_t = 100.0;
        for t in &mut c.thermal.temperatures {
            *t = mat.t_ref + delta_t;
        }
        c.solve_mechanics().expect("mechanics");
        let eps = mat.alpha * delta_t;
        let (hx, hy) = (c.hx(), c.hy());
        for i in 0..5 {
            for j in 0..4 {
                let node = c.thermal.idx(i, j);
                let [ux, uy] = c.displacements[node];
                assert!((ux - eps * i as f64 * hx).abs() < 1e-12, "ux at ({i},{j})");
                assert!((uy - eps * j as f64 * hy).abs() < 1e-12, "uy at ({i},{j})");
            }
        }
        let scale = mat.e0 * eps;
        for node in 0..c.n_nodes() {
            let s = c.nodal_stress(node).expect("stress");
            assert!(s.iter().all(|v| v.abs() < 1e-7 * scale), "σ = {s:?}");
            assert!(c.von_mises_stress(node).expect("vm") < 1e-7 * scale);
        }
    }

    #[test]
    fn clamped_plate_uniform_heating_biaxial_stress() {
        let mut mat = steel();
        mat.beta = 2e-4; // temperature-dependent modulus
        let mut c = ThermoMechanicalCoupling::new(4, 4, mat.clone()).expect("valid problem");
        c.clamp_all_edges();
        let delta_t = 80.0;
        let t = mat.t_ref + delta_t;
        for v in &mut c.thermal.temperatures {
            *v = t;
        }
        c.solve_mechanics().expect("mechanics");
        let expected = -mat.modulus_at(t) * mat.alpha * delta_t / (1.0 - mat.nu);
        for node in 0..c.n_nodes() {
            let [sx, sy, txy] = c.nodal_stress(node).expect("stress");
            assert!((sx - expected).abs() < 1e-9 * expected.abs());
            assert!((sy - expected).abs() < 1e-9 * expected.abs());
            assert!(txy.abs() < 1e-9 * expected.abs());
            let vm = c.von_mises_stress(node).expect("vm");
            assert!((vm - expected.abs()).abs() < 1e-8 * expected.abs());
        }
    }

    #[test]
    fn mechanics_requires_constraints() {
        let mut c = ThermoMechanicalCoupling::new(3, 3, steel()).expect("valid problem");
        c.fixed_dofs.clear();
        assert!(matches!(c.solve_mechanics(), Err(Error::InvalidInput(_))));
        c.fixed_dofs = vec![1000];
        assert!(matches!(
            c.solve_mechanics(),
            Err(Error::DofOutOfRange { dof: 1000, n: 18 })
        ));
    }

    #[test]
    fn implicit_heat_reaches_linear_steady_state() {
        let mut c = ThermoMechanicalCoupling::new(9, 4, steel()).expect("valid problem");
        c.set_left_temperature(600.0);
        c.set_right_temperature(300.0);
        // Steps far beyond the explicit Fourier limit stay stable.
        for _ in 0..50 {
            c.advance_heat(1e6, 0.0).expect("heat");
        }
        let hx = c.hx();
        for i in 0..9 {
            for j in 0..4 {
                let x = i as f64 * hx;
                let t = c.thermal.temperatures[c.thermal.idx(i, j)];
                assert!((t - (600.0 - 300.0 * x)).abs() < 1e-6, "T({i},{j}) = {t}");
            }
        }
    }

    #[test]
    fn implicit_heat_large_step_respects_maximum_principle() {
        let mut c = ThermoMechanicalCoupling::new(7, 7, steel()).expect("valid problem");
        for (k, t) in c.thermal.temperatures.iter_mut().enumerate() {
            *t = if k % 3 == 0 { 600.0 } else { 300.0 };
        }
        c.set_left_temperature(450.0);
        for _ in 0..5 {
            c.advance_heat(10.0, 0.0).expect("heat");
            assert!(
                c.thermal
                    .temperatures
                    .iter()
                    .all(|&t| (300.0 - 1e-9..=600.0 + 1e-9).contains(&t))
            );
        }
    }

    #[test]
    fn insulated_plate_one_way_energy_balance() {
        let mat = steel();
        let mut c = ThermoMechanicalCoupling::new(4, 4, mat.clone()).expect("valid problem");
        c.thermoelastic_coupling = false;
        let (dt, q) = (2.0, 1e6);
        let report = c.staggered_step(dt, q).expect("step");
        assert_eq!(report.iterations, 1);
        let expected = mat.t_ref + q * dt / (mat.rho * mat.cp);
        assert!(
            c.thermal
                .temperatures
                .iter()
                .all(|&t| (t - expected).abs() < 1e-9)
        );
        // Mechanics follows the new temperature: free expansion.
        let eps = mat.alpha * (expected - mat.t_ref);
        let corner = c.thermal.idx(3, 3);
        assert!((c.displacements[corner][0] - eps).abs() < 1e-12);
        assert!((c.displacements[corner][1] - eps).abs() < 1e-12);
    }

    #[test]
    fn two_way_coupled_heating_matches_thermoelastic_balance() {
        // Uniform heating of a free, insulated plate: ρc ΔT + T₀β·2αΔT = Q dt.
        let mat = steel();
        let mut c = ThermoMechanicalCoupling::new(4, 4, mat.clone()).expect("valid problem");
        c.tol = 1e-12;
        let (dt, q) = (1.0, 5e6);
        let beta = mat.e0 * mat.alpha / (1.0 - mat.nu);
        let rho_cp = mat.rho * mat.cp;
        let delta = q * dt / (rho_cp + 2.0 * mat.t_ref * beta * mat.alpha);
        let report = c.staggered_step(dt, q).expect("step");
        assert!(report.converged, "{report:?}");
        assert!(report.iterations > 1 && report.iterations < 15);
        for &t in &c.thermal.temperatures {
            assert!((t - (mat.t_ref + delta)).abs() < 1e-9 * delta);
        }
        // Thermo-elastic cooling makes the coupled rise smaller than the uncoupled one.
        assert!(delta < q * dt / rho_cp);
        let reports = c.solve(dt, q, 3).expect("steps");
        assert!(reports.iter().all(|r| r.converged));
        for &t in &c.thermal.temperatures {
            assert!((t - (mat.t_ref + 4.0 * delta)).abs() < 1e-8 * delta);
        }
    }

    #[test]
    fn heated_edge_produces_thermal_stress_in_clamped_plate() {
        let mut c = ThermoMechanicalCoupling::new(6, 6, steel()).expect("valid problem");
        c.clamp_all_edges();
        c.set_left_temperature(500.0);
        let report = c.staggered_step(5.0, 0.0).expect("step");
        assert!(report.converged);
        // The heated, clamped edge carries ≈ EαΔT/(1−ν) compressive stress.
        let node = c.thermal.idx(0, 3);
        let vm = c.von_mises_stress(node).expect("vm");
        assert!(vm > 1e8, "heated region must be stressed, got {vm}");
        // The thermal stress decays away from the heated edge.
        let far = c.thermal.idx(5, 3);
        let vm_far = c.von_mises_stress(far).expect("vm");
        assert!(vm_far < 0.2 * vm, "far field {vm_far} vs heated edge {vm}");
        assert!(c.von_mises_stress(10_000).is_none());
        assert!(c.element_stress(5, 0, 0.0, 0.0).is_none());
    }
}
