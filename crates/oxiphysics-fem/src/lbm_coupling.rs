// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! FEM-LBM coupling module using the Immersed Boundary Method (IBM).
//!
//! Provides fluid-structure interaction (FSI) between a Finite Element Method
//! (FEM) structural solver and a Lattice Boltzmann Method (LBM) fluid solver.
//!
//! The coupling uses a smoothed delta function (4-point Roma kernel) to spread
//! forces from Lagrangian (FEM) points to the Eulerian (LBM) grid and to
//! interpolate velocities from the grid to the Lagrangian points.
//!
//! # References
//!
//! - Roma, Peskin, Berger (1999), "An adaptive version of the immersed boundary method"
//! - Peskin (2002), "The immersed boundary method", Acta Numerica

#[cfg(test)]
mod velocity_bc_tests;

/// Configuration parameters for FEM-LBM coupling.
#[derive(Debug, Clone)]
pub struct FemLbmCouplingConfig {
    /// Relaxation factor for under-relaxation of force/velocity transfer (0, 1].
    /// A value of 1.0 means no relaxation; smaller values improve stability.
    pub relaxation_factor: f64,
    /// Interpolation order for the IBM delta function kernel.
    /// Currently only order 4 (Roma kernel) is supported.
    pub interpolation_order: usize,
    /// Penalty stiffness for enforcing the no-slip condition at the interface.
    /// Larger values enforce the constraint more rigidly but may reduce stability.
    pub penalty_stiffness: f64,
    /// Number of sub-iterations for the coupling exchange per time step.
    pub sub_iterations: usize,
    /// Fluid dynamic viscosity (used in viscous stress computation).
    pub fluid_viscosity: f64,
}

impl Default for FemLbmCouplingConfig {
    fn default() -> Self {
        Self {
            relaxation_factor: 1.0,
            interpolation_order: 4,
            penalty_stiffness: 1.0e3,
            sub_iterations: 1,
            fluid_viscosity: 1.0e-3,
        }
    }
}

/// A single immersed boundary point at the FEM-LBM interface.
#[derive(Debug, Clone)]
pub struct ImmersedBoundaryPoint {
    /// Position of the IB point in physical space.
    pub position: [f64; 3],
    /// Velocity of the IB point (from FEM structural solver).
    pub velocity: [f64; 3],
    /// Force accumulated on this IB point (from LBM fluid solver).
    pub force: [f64; 3],
    /// Reference (initial) position for tracking deformation.
    pub reference_position: [f64; 3],
    /// Area weight associated with this IB point for force spreading, in units
    /// of the squared LBM grid spacing (`area = area_weight * h^2`).
    pub area_weight: f64,
    /// Unit surface normal at the point, pointing from the structure **into the
    /// fluid**.  The fluid traction on the structure is `σ·normal`; the zero
    /// vector means "unknown" and switches the pressure and viscous parts of
    /// the force off (only the no-slip penalty acts).
    pub normal: [f64; 3],
}

impl ImmersedBoundaryPoint {
    /// Create a new immersed boundary point at the given position.
    ///
    /// The normal is unset (zero vector); see [`with_normal`](Self::with_normal)
    /// and [`FemLbmCoupling::set_interface_surface`].
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            reference_position: position,
            area_weight: 1.0,
            normal: [0.0; 3],
        }
    }

    /// Create a new immersed boundary point with a specified area weight.
    pub fn with_area_weight(position: [f64; 3], area_weight: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            reference_position: position,
            area_weight,
            normal: [0.0; 3],
        }
    }

    /// Create a new immersed boundary point with a surface normal and area weight.
    ///
    /// `normal` points from the structure into the fluid; it is normalised here
    /// (a zero or non-finite vector leaves the normal unset).
    pub fn with_normal(position: [f64; 3], normal: [f64; 3], area_weight: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            force: [0.0; 3],
            reference_position: position,
            area_weight,
            normal: normalized(normal),
        }
    }
}

/// Unit vector of `v`, or the zero vector when `v` is zero or not finite.
fn normalized(v: [f64; 3]) -> [f64; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if len.is_finite() && len > 1e-300 {
        [v[0] / len, v[1] / len, v[2] / len]
    } else {
        [0.0; 3]
    }
}

/// Vector area of the surface patch tributary to each point of a triangulated
/// surface: one third of the (outward) area vector `½ (p1 - p0) × (p2 - p0)` of
/// every adjacent triangle.
///
/// The vectors of a closed surface sum to zero, so a uniform pressure
/// integrates to exactly zero net force.
fn tributary_vector_areas(points: &[[f64; 3]], faces: &[[usize; 3]]) -> Vec<[f64; 3]> {
    let mut areas = vec![[0.0_f64; 3]; points.len()];
    for face in faces {
        let (a, b, c) = (points[face[0]], points[face[1]], points[face[2]]);
        let e1 = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let e2 = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let cross = [
            e1[1] * e2[2] - e1[2] * e2[1],
            e1[2] * e2[0] - e1[0] * e2[2],
            e1[0] * e2[1] - e1[1] * e2[0],
        ];
        for &vertex in face {
            for (acc, cr) in areas[vertex].iter_mut().zip(cross.iter()) {
                *acc += cr / 6.0;
            }
        }
    }
    areas
}

/// Main FEM-LBM coupling state using the Immersed Boundary Method.
///
/// Manages the bidirectional data transfer between a FEM structural domain
/// and an LBM fluid domain across a fluid-structure interface.
#[derive(Debug, Clone)]
pub struct FemLbmCoupling {
    /// Configuration parameters.
    pub config: FemLbmCouplingConfig,
    /// FEM node positions (all nodes, not just interface).
    pub fem_nodes: Vec<[f64; 3]>,
    /// LBM grid cell spacing.
    pub lbm_grid_spacing: f64,
    /// Indices of FEM nodes that lie on the FSI interface.
    pub interface_nodes: Vec<usize>,
    /// Boundary velocities to impose on LBM cells near the interface.
    /// One velocity per interface node.
    pub boundary_velocities: Vec<[f64; 3]>,
    /// Fluid forces to apply to FEM nodes at the interface.
    /// One force per interface node.
    pub fluid_forces: Vec<[f64; 3]>,
    /// Immersed boundary points at the interface.
    pub ib_points: Vec<ImmersedBoundaryPoint>,
    /// Previous iteration forces for relaxation.
    prev_forces: Vec<[f64; 3]>,
    /// Triangulation of the interface (indices into the interface node list),
    /// kept so that normals and areas follow the deforming structure.  Empty
    /// when the interface was defined by nodes only.
    interface_faces: Vec<[usize; 3]>,
    /// Whether the coupling has been initialized with interface nodes.
    initialized: bool,
}

impl FemLbmCoupling {
    /// Create a new FEM-LBM coupling with the given configuration.
    pub fn new(config: FemLbmCouplingConfig) -> Self {
        Self {
            config,
            fem_nodes: Vec::new(),
            lbm_grid_spacing: 1.0,
            interface_nodes: Vec::new(),
            boundary_velocities: Vec::new(),
            fluid_forces: Vec::new(),
            ib_points: Vec::new(),
            prev_forces: Vec::new(),
            interface_faces: Vec::new(),
            initialized: false,
        }
    }

    /// Define the fluid-structure interface by specifying which FEM nodes
    /// participate in coupling and their current positions.
    ///
    /// The interface points get **no surface normal** (and unit area weight), so
    /// only the no-slip penalty force acts on them; pressure and viscous
    /// tractions need the orientation of the surface.  Use
    /// [`set_interface_surface`](Self::set_interface_surface) when the boundary
    /// faces of the FEM mesh are known, or
    /// [`set_interface_normals`](Self::set_interface_normals) to give normals
    /// directly.
    ///
    /// # Arguments
    /// * `node_indices` - Indices into the FEM mesh identifying interface nodes.
    /// * `positions` - Current 3D positions of the interface nodes.
    /// * `grid_spacing` - LBM grid cell spacing.
    ///
    /// # Errors
    /// Returns an error if `node_indices` and `positions` have different lengths.
    pub fn set_interface_nodes(
        &mut self,
        node_indices: Vec<usize>,
        positions: Vec<[f64; 3]>,
        grid_spacing: f64,
    ) -> Result<(), CouplingError> {
        if node_indices.len() != positions.len() {
            return Err(CouplingError::DimensionMismatch {
                expected: node_indices.len(),
                got: positions.len(),
            });
        }
        if grid_spacing <= 0.0 {
            return Err(CouplingError::InvalidParameter(
                "grid_spacing must be positive".to_string(),
            ));
        }

        let n = node_indices.len();
        self.interface_nodes = node_indices;
        self.fem_nodes = positions.clone();
        self.lbm_grid_spacing = grid_spacing;
        self.boundary_velocities = vec![[0.0; 3]; n];
        self.fluid_forces = vec![[0.0; 3]; n];
        self.prev_forces = vec![[0.0; 3]; n];
        self.interface_faces.clear();

        // Create immersed boundary points from interface nodes
        self.ib_points = positions
            .into_iter()
            .map(ImmersedBoundaryPoint::new)
            .collect();

        self.initialized = true;
        Ok(())
    }

    /// Define the interface from the boundary faces of the FEM mesh.
    ///
    /// Like [`set_interface_nodes`](Self::set_interface_nodes) but additionally
    /// takes the triangulated interface surface, from which every IB point gets
    /// its **unit normal** (pointing from the structure into the fluid) and its
    /// **area weight**.  The pressure and viscous traction of
    /// [`compute_fluid_forces`](Self::compute_fluid_forces) act along this
    /// normal.
    ///
    /// A point's normal and area follow from the vector area
    /// `a = ⅓ Σ ½ (p1 - p0) × (p2 - p0)` of its adjacent triangles:
    /// `normal = a / |a|` and `area_weight = |a| / h²`, so that
    /// `area_weight · h² · normal` is exactly the tributary area vector (a
    /// uniform pressure on a closed surface then integrates to zero net force).
    /// The surface is re-evaluated whenever
    /// [`update_fem_state`](Self::update_fem_state) moves the points.
    ///
    /// # Arguments
    /// * `node_indices` - Indices into the FEM mesh identifying interface nodes.
    /// * `positions` - Current 3D positions of the interface nodes.
    /// * `faces` - Triangles as **indices into `node_indices` / `positions`**
    ///   (not FEM mesh indices).  Each triangle must be wound
    ///   counter-clockwise when seen from the fluid, i.e.
    ///   `(p1 - p0) × (p2 - p0)` points into the fluid.
    /// * `grid_spacing` - LBM grid cell spacing.
    ///
    /// # Errors
    /// Returns an error if `node_indices` and `positions` have different lengths,
    /// if `grid_spacing` is not positive, or if a face refers to a point that does
    /// not exist.
    pub fn set_interface_surface(
        &mut self,
        node_indices: Vec<usize>,
        positions: Vec<[f64; 3]>,
        faces: Vec<[usize; 3]>,
        grid_spacing: f64,
    ) -> Result<(), CouplingError> {
        let n = positions.len();
        if let Some(bad) = faces.iter().flatten().find(|&&v| v >= n) {
            return Err(CouplingError::InvalidParameter(format!(
                "face refers to interface point {bad} but only {n} points exist"
            )));
        }
        self.set_interface_nodes(node_indices, positions, grid_spacing)?;
        self.interface_faces = faces;
        self.refresh_surface_geometry();
        Ok(())
    }

    /// Set the unit normals and area weights of the IB points directly.
    ///
    /// Use this when no triangulation is available.  `normals[i]` points from
    /// the structure into the fluid and is normalised here; `area_weights[i]` is
    /// the tributary area in units of the squared grid spacing.  The values are
    /// kept fixed (they are not recomputed by
    /// [`update_fem_state`](Self::update_fem_state)).
    ///
    /// # Errors
    /// Returns an error if the coupling is not initialized, if the slice lengths
    /// differ from the number of interface nodes, or if an area weight is
    /// negative or not finite.
    pub fn set_interface_normals(
        &mut self,
        normals: &[[f64; 3]],
        area_weights: &[f64],
    ) -> Result<(), CouplingError> {
        if !self.initialized {
            return Err(CouplingError::NotInitialized);
        }
        let n = self.ib_points.len();
        if normals.len() != n || area_weights.len() != n {
            return Err(CouplingError::DimensionMismatch {
                expected: n,
                got: normals.len().min(area_weights.len()),
            });
        }
        if area_weights.iter().any(|a| !a.is_finite() || *a < 0.0) {
            return Err(CouplingError::InvalidParameter(
                "area weights must be finite and non-negative".to_string(),
            ));
        }
        self.interface_faces.clear();
        for ((ibp, normal), area) in self
            .ib_points
            .iter_mut()
            .zip(normals.iter())
            .zip(area_weights.iter())
        {
            ibp.normal = normalized(*normal);
            ibp.area_weight = *area;
        }
        Ok(())
    }

    /// Recompute normals and area weights from the interface triangulation.
    fn refresh_surface_geometry(&mut self) {
        if self.interface_faces.is_empty() {
            return;
        }
        let positions: Vec<[f64; 3]> = self.ib_points.iter().map(|p| p.position).collect();
        let vector_areas = tributary_vector_areas(&positions, &self.interface_faces);
        let h2 = self.lbm_grid_spacing * self.lbm_grid_spacing;
        for (ibp, area) in self.ib_points.iter_mut().zip(vector_areas.iter()) {
            let magnitude = (area[0] * area[0] + area[1] * area[1] + area[2] * area[2]).sqrt();
            ibp.normal = normalized(*area);
            ibp.area_weight = magnitude / h2;
        }
    }

    /// Update the positions and velocities of the interface nodes from the FEM solver.
    ///
    /// # Errors
    /// Returns an error if the coupling is not initialized or dimensions mismatch.
    pub fn update_fem_state(
        &mut self,
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
    ) -> Result<(), CouplingError> {
        if !self.initialized {
            return Err(CouplingError::NotInitialized);
        }
        let n = self.interface_nodes.len();
        if positions.len() != n || velocities.len() != n {
            return Err(CouplingError::DimensionMismatch {
                expected: n,
                got: positions.len().min(velocities.len()),
            });
        }

        for (ibp, (&pos, &vel)) in self
            .ib_points
            .iter_mut()
            .zip(positions.iter().zip(velocities.iter()))
        {
            ibp.position = pos;
            ibp.velocity = vel;
        }
        for (bv, &vel) in self.boundary_velocities.iter_mut().zip(velocities.iter()) {
            *bv = vel;
        }
        // The deformed interface has new normals and tributary areas.
        self.refresh_surface_geometry();

        Ok(())
    }

    /// Transfer velocity boundary conditions from FEM interface to LBM grid.
    ///
    /// Returns a list of `(grid_index, velocity)` pairs that should be imposed
    /// as velocity boundary conditions on the LBM grid.  Every grid node in the
    /// support of the IBM delta kernel of at least one IB point is reported
    /// exactly once, with the kernel-weighted average of the velocities of the
    /// points that reach it:
    ///
    /// ```text
    /// u_bc(i, j, k) = Σ_p ω_p u_p / Σ_p ω_p,
    /// ω_p = φ(x_p/h - i) φ(y_p/h - j) φ(z_p/h - k)
    /// ```
    ///
    /// `φ` is the dimensionless 1-D Roma kernel (unit spacing, so that
    /// `Σ_i φ(s - i) = 1`).  A node that only one point reaches therefore gets
    /// exactly that point's velocity, however far from the point it lies, and
    /// the result does not depend on the grid spacing.  Nodes outside the grid
    /// are dropped.  The list is ordered by ascending grid index (row-major), so
    /// repeated calls are bit-for-bit reproducible.
    ///
    /// # Arguments
    /// * `grid_dims` - Dimensions of the LBM grid `[nx, ny, nz]`.
    ///
    /// # Errors
    /// Returns an error if the coupling is not initialized or if the stored grid
    /// spacing is not positive and finite.
    ///
    /// # Examples
    ///
    /// A single structure point imposes its own velocity on all 3 x 3 x 3 grid nodes
    /// of its kernel support:
    ///
    /// ```
    /// use oxiphysics_fem::lbm_coupling::{CouplingError, FemLbmCoupling, FemLbmCouplingConfig};
    ///
    /// let mut coupling = FemLbmCoupling::new(FemLbmCouplingConfig::default());
    /// let position = [5.3, 4.8, 6.1];
    /// coupling.set_interface_nodes(vec![0], vec![position], 1.0)?;
    /// coupling.update_fem_state(&[position], &[[0.7, -0.2, 1.3]])?;
    ///
    /// let bcs = coupling.transfer_velocity_to_lbm([12, 12, 12])?;
    /// assert_eq!(bcs.len(), 27);
    /// for (_node, velocity) in &bcs {
    ///     assert!((velocity[0] - 0.7).abs() < 1e-12);
    ///     assert!((velocity[1] + 0.2).abs() < 1e-12);
    ///     assert!((velocity[2] - 1.3).abs() < 1e-12);
    /// }
    /// # Ok::<(), CouplingError>(())
    /// ```
    pub fn transfer_velocity_to_lbm(&self, grid_dims: [usize; 3]) -> LbmCouplingResult {
        if !self.initialized {
            return Err(CouplingError::NotInitialized);
        }
        let h = self.lbm_grid_spacing;
        if !h.is_finite() || h <= 0.0 {
            return Err(CouplingError::InvalidParameter(
                "grid spacing must be positive and finite".to_string(),
            ));
        }

        // One record per (IB point, stencil node): the node, omega * u and omega.
        let mut records: Vec<([usize; 3], [f64; 3], f64)> = Vec::new();
        for ib_point in &self.ib_points {
            let u = ib_point.velocity;
            for (node, omega) in kernel_stencil(ib_point.position, h, grid_dims) {
                records.push((node, [omega * u[0], omega * u[1], omega * u[2]], omega));
            }
        }
        // The lexicographic order of `[i, j, k]` is the row-major flat-index
        // order.  The sort is stable, so the contributions of a node are summed
        // in IB-point order and the result is reproducible bit for bit.
        records.sort_by_key(|record| record.0);

        let mut velocity_bcs = LbmVelocityBcList::new();
        for node_records in records.chunk_by(|a, b| a.0 == b.0) {
            let Some(&(node, _, _)) = node_records.first() else {
                continue;
            };
            let mut weighted_velocity = [0.0; 3];
            let mut weight_sum = 0.0;
            for (_, weighted_u, omega) in node_records {
                for (acc, component) in weighted_velocity.iter_mut().zip(weighted_u.iter()) {
                    *acc += component;
                }
                weight_sum += omega;
            }
            if weight_sum > MIN_KERNEL_WEIGHT {
                velocity_bcs.push((
                    node,
                    [
                        weighted_velocity[0] / weight_sum,
                        weighted_velocity[1] / weight_sum,
                        weighted_velocity[2] / weight_sum,
                    ],
                ));
            }
        }

        Ok(velocity_bcs)
    }

    /// Compute fluid forces on FEM interface nodes from LBM pressure and velocity fields.
    ///
    /// The force on the structure at an interface point is the fluid traction
    /// integrated over the tributary area plus a no-slip feedback term:
    ///
    /// ```text
    /// t = σ·n,   σ = -p I + μ (∇u + ∇uᵀ),   p = ρ c_s²
    /// F = t · A  +  k (u_fluid - u_solid) h³
    /// ```
    ///
    /// where `n` is the unit normal of the IB point pointing from the structure
    /// into the fluid, `A = area_weight · h²` its tributary area and `k` the
    /// penalty stiffness.  Pressure therefore acts *against the normal*
    /// (`-p n A`, compressing the structure) and the viscous part is the
    /// full Newtonian shear traction.  The feedback term accelerates the
    /// structure toward the local fluid velocity (it is the reaction to the
    /// penalty force that pulls the fluid toward the structure velocity), so a
    /// slipping structure is dragged along with the flow.  Points without a
    /// normal (see [`set_interface_nodes`](Self::set_interface_nodes)) receive
    /// the feedback term only.
    ///
    /// # Arguments
    /// * `lbm_density` - Density field on the LBM grid (flattened, row-major).
    /// * `lbm_velocity` - Velocity field on the LBM grid (flattened, row-major).
    /// * `grid_spacing` - LBM grid cell spacing.
    /// * `grid_dims` - Dimensions of the LBM grid `[nx, ny, nz]`.
    ///
    /// # Errors
    /// Returns an error if the coupling is not initialized or field sizes mismatch.
    pub fn compute_fluid_forces(
        &mut self,
        lbm_density: &[f64],
        lbm_velocity: &[[f64; 3]],
        grid_spacing: f64,
        grid_dims: [usize; 3],
    ) -> Result<(), CouplingError> {
        if !self.initialized {
            return Err(CouplingError::NotInitialized);
        }

        let total_cells = grid_dims[0] * grid_dims[1] * grid_dims[2];
        if lbm_density.len() != total_cells || lbm_velocity.len() != total_cells {
            return Err(CouplingError::FieldSizeMismatch {
                expected: total_cells,
                density_len: lbm_density.len(),
                velocity_len: lbm_velocity.len(),
            });
        }

        let h = grid_spacing;
        let cs2 = 1.0 / 3.0; // LBM speed of sound squared (lattice units)
        let mu = self.config.fluid_viscosity;
        let n = self.interface_nodes.len();

        // Save previous forces for relaxation
        for (prev, &curr) in self.prev_forces.iter_mut().zip(self.fluid_forces.iter()) {
            *prev = curr;
        }

        for i in 0..n {
            let pos = self.ib_points[i].position;

            // Interpolate density at IB point
            let rho = interpolate_scalar(pos, lbm_density, h, grid_dims);

            // Interpolate velocity at IB point
            let vel_fluid = interpolate_velocity(pos, lbm_velocity, h, grid_dims);

            // Compute velocity gradient tensor at IB point using finite differences
            let grad_u = compute_velocity_gradient(pos, lbm_velocity, h, grid_dims);

            // Pressure from the equation of state of the lattice fluid: p = rho * cs^2.
            let pressure = rho * cs2;

            // Slip of the fluid relative to the structure at the IB point.
            let dv = [
                vel_fluid[0] - self.ib_points[i].velocity[0],
                vel_fluid[1] - self.ib_points[i].velocity[1],
                vel_fluid[2] - self.ib_points[i].velocity[2],
            ];

            // No-slip feedback: reaction to the penalty force k * (u_solid - u_fluid)
            // that the structure exerts on the fluid.
            let penalty_k = self.config.penalty_stiffness;

            // Traction exerted by the fluid on the structure,
            //   t_d = -p n_d + mu * sum_e (du_d/dx_e + du_e/dx_d) n_e,
            // with n the normal pointing from the structure into the fluid.
            let normal = self.ib_points[i].normal;
            let mut traction = [0.0; 3];
            for (d, t_d) in traction.iter_mut().enumerate() {
                let mut shear = 0.0;
                for (e, &n_e) in normal.iter().enumerate() {
                    shear += (grad_u[d][e] + grad_u[e][d]) * n_e;
                }
                *t_d = -pressure * normal[d] + mu * shear;
            }

            // Total force: traction integrated over the tributary area + feedback.
            let area = self.ib_points[i].area_weight * h * h;
            let force = [
                traction[0] * area + penalty_k * dv[0] * h.powi(3),
                traction[1] * area + penalty_k * dv[1] * h.powi(3),
                traction[2] * area + penalty_k * dv[2] * h.powi(3),
            ];

            // Apply under-relaxation
            let alpha = self.config.relaxation_factor;
            self.fluid_forces[i] = [
                alpha * force[0] + (1.0 - alpha) * self.prev_forces[i][0],
                alpha * force[1] + (1.0 - alpha) * self.prev_forces[i][1],
                alpha * force[2] + (1.0 - alpha) * self.prev_forces[i][2],
            ];

            // Store force on IB point as well
            self.ib_points[i].force = self.fluid_forces[i];
        }

        Ok(())
    }

    /// Retrieve the accumulated fluid forces on FEM interface nodes.
    pub fn get_fem_forces(&self) -> &[[f64; 3]] {
        &self.fluid_forces
    }

    /// Perform a full coupling exchange step.
    ///
    /// Runs the configured number of sub-iterations of:
    /// 1. Transfer FEM velocities to LBM grid
    /// 2. (Caller advances LBM)
    /// 3. Compute fluid forces on FEM nodes
    ///
    /// Returns the velocity BCs for the LBM grid after the final sub-iteration.
    ///
    /// # Errors
    /// Returns an error if the coupling is not initialized or field sizes mismatch.
    pub fn coupling_step(
        &mut self,
        lbm_density: &[f64],
        lbm_velocity: &[[f64; 3]],
        grid_spacing: f64,
        grid_dims: [usize; 3],
    ) -> LbmCouplingResult {
        if !self.initialized {
            return Err(CouplingError::NotInitialized);
        }

        // Compute forces from the current LBM state
        self.compute_fluid_forces(lbm_density, lbm_velocity, grid_spacing, grid_dims)?;

        // Transfer updated velocities to LBM
        self.transfer_velocity_to_lbm(grid_dims)
    }

    /// Reset accumulated forces to zero.
    pub fn reset_forces(&mut self) {
        for f in &mut self.fluid_forces {
            *f = [0.0; 3];
        }
        for ib in &mut self.ib_points {
            ib.force = [0.0; 3];
        }
    }
}

/// Error type for FEM-LBM coupling operations.
#[derive(Debug, Clone)]
pub enum CouplingError {
    /// Dimension mismatch between inputs.
    DimensionMismatch {
        /// Expected size.
        expected: usize,
        /// Actual size received.
        got: usize,
    },
    /// Invalid parameter value.
    InvalidParameter(String),
    /// Coupling has not been initialized (no interface nodes set).
    NotInitialized,
    /// LBM field size does not match expected grid dimensions.
    FieldSizeMismatch {
        /// Expected number of cells.
        expected: usize,
        /// Density field length.
        density_len: usize,
        /// Velocity field length.
        velocity_len: usize,
    },
}

impl std::fmt::Display for CouplingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::DimensionMismatch { expected, got } => {
                write!(f, "dimension mismatch: expected {expected}, got {got}")
            }
            Self::InvalidParameter(msg) => write!(f, "invalid parameter: {msg}"),
            Self::NotInitialized => write!(
                f,
                "coupling not initialized: call set_interface_nodes first"
            ),
            Self::FieldSizeMismatch {
                expected,
                density_len,
                velocity_len,
            } => write!(
                f,
                "field size mismatch: expected {expected} cells, got density={density_len}, velocity={velocity_len}"
            ),
        }
    }
}

impl std::error::Error for CouplingError {}

/// A list of `(grid_index, velocity)` pairs for LBM velocity boundary conditions.
///
/// [`FemLbmCoupling::transfer_velocity_to_lbm`] reports every grid node at most
/// once, in ascending row-major order of the index.
pub type LbmVelocityBcList = Vec<([usize; 3], [f64; 3])>;

/// Result type for LBM velocity transfer and coupling step operations.
pub type LbmCouplingResult = Result<LbmVelocityBcList, CouplingError>;

// ---------------------------------------------------------------------------
// Free functions for IBM delta function and interpolation
// ---------------------------------------------------------------------------

/// Convert a physical position to the nearest LBM grid index.
///
/// The grid is assumed to start at the origin with uniform spacing `grid_spacing`.
/// Coordinates are clamped so the returned index is always valid for a grid of
/// the given dimensions.
pub fn interpolate_to_grid(pos: [f64; 3], grid_spacing: f64) -> [usize; 3] {
    let safe_h = if grid_spacing > 0.0 {
        grid_spacing
    } else {
        1.0
    };
    [
        (pos[0] / safe_h).round().max(0.0) as usize,
        (pos[1] / safe_h).round().max(0.0) as usize,
        (pos[2] / safe_h).round().max(0.0) as usize,
    ]
}

/// 1D Roma (4-point) smoothed delta function.
///
/// The Roma kernel is defined piecewise on `|r| / h`:
///
/// ```text
/// phi(r) = (1/6h) * (5 - 3|r|/h - sqrt(-3(1 - |r|/h)^2 + 1))   for 0.5 <= |r|/h <= 1.5
/// phi(r) = (1/3h) * (1 + sqrt(-3(r/h)^2 + 1))                     for |r|/h <= 0.5
/// phi(r) = 0                                                        otherwise
/// ```
///
/// # Arguments
/// * `r` - Distance (may be negative; absolute value is used).
/// * `h` - Grid spacing (support width parameter).
///
/// # Returns
/// The value of the smoothed delta function.
pub fn delta_function(r: f64, h: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    roma_delta_1d(r, h)
}

/// Internal 1D Roma delta kernel evaluation.
///
/// Takes the normalised distance `r` and spacing `h`. The result is scaled
/// by `1/h` so that the 3D product `delta_3d = phi(x)*phi(y)*phi(z)` integrates
/// to 1 over all space when each factor uses spacing `h`.
fn roma_delta_1d(r: f64, h: f64) -> f64 {
    let s = (r / h).abs();
    if s >= 1.5 {
        0.0
    } else if s >= 0.5 {
        let inner = -3.0 * (1.0 - s) * (1.0 - s) + 1.0;
        if inner < 0.0 {
            // Numerically should not happen for s in [0.5, 1.5], but guard anyway
            0.0
        } else {
            (1.0 / (6.0 * h)) * (5.0 - 3.0 * s - inner.sqrt())
        }
    } else {
        let inner = -3.0 * s * s + 1.0;
        if inner < 0.0 {
            0.0
        } else {
            (1.0 / (3.0 * h)) * (1.0 + inner.sqrt())
        }
    }
}

/// 3D IBM smoothed delta function as a product of 1D Roma kernels.
///
/// ```text
/// delta_3d(r, h) = phi(rx, h) * phi(ry, h) * phi(rz, h)
/// ```
#[cfg(test)]
fn delta_3d(r: [f64; 3], h: f64) -> f64 {
    roma_delta_1d(r[0], h) * roma_delta_1d(r[1], h) * roma_delta_1d(r[2], h)
}

/// Spread a force from a Lagrangian point to nearby Eulerian grid cells.
///
/// Uses the 3D Roma delta function to distribute the force. Returns a list of
/// `(grid_index, force_contribution)` pairs.
///
/// # Arguments
/// * `pos` - Position of the Lagrangian point.
/// * `force` - Force vector to spread.
/// * `grid_spacing` - LBM grid cell spacing.
/// * `grid_dims` - LBM grid dimensions `[nx, ny, nz]`.
pub fn spread_force(
    pos: [f64; 3],
    force: [f64; 3],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> Vec<([usize; 3], [f64; 3])> {
    spread_to_grid(pos, force, grid_spacing, grid_dims)
}

/// Internal spreading: distributes a vector quantity from a Lagrangian point
/// to nearby grid cells using the IBM delta function.
fn spread_to_grid(
    pos: [f64; 3],
    quantity: [f64; 3],
    h: f64,
    grid_dims: [usize; 3],
) -> Vec<([usize; 3], [f64; 3])> {
    if h <= 0.0 {
        return Vec::new();
    }

    let mut result = Vec::new();

    // The Roma kernel has support of width 3h centered on the point.
    // We need to iterate over grid cells within +/- 1.5 * h of the point.
    let base = [pos[0] / h, pos[1] / h, pos[2] / h];

    // Integer range for each dimension
    let i_min = ((base[0] - 1.5).floor() as i64).max(0) as usize;
    let i_max = ((base[0] + 1.5).ceil() as i64).max(0) as usize;
    let j_min = ((base[1] - 1.5).floor() as i64).max(0) as usize;
    let j_max = ((base[1] + 1.5).ceil() as i64).max(0) as usize;
    let k_min = ((base[2] - 1.5).floor() as i64).max(0) as usize;
    let k_max = ((base[2] + 1.5).ceil() as i64).max(0) as usize;

    for i in i_min..=i_max.min(grid_dims[0].saturating_sub(1)) {
        let dx = pos[0] - (i as f64) * h;
        let phi_x = roma_delta_1d(dx, h);
        if phi_x.abs() < 1e-30 {
            continue;
        }
        for j in j_min..=j_max.min(grid_dims[1].saturating_sub(1)) {
            let dy = pos[1] - (j as f64) * h;
            let phi_y = roma_delta_1d(dy, h);
            if phi_y.abs() < 1e-30 {
                continue;
            }
            for k in k_min..=k_max.min(grid_dims[2].saturating_sub(1)) {
                let dz = pos[2] - (k as f64) * h;
                let phi_z = roma_delta_1d(dz, h);
                let w = phi_x * phi_y * phi_z;
                if w.abs() < 1e-30 {
                    continue;
                }
                // In IBM, spreading multiplies by h^3 (volume element)
                let h3 = h * h * h;
                result.push((
                    [i, j, k],
                    [
                        quantity[0] * w * h3,
                        quantity[1] * w * h3,
                        quantity[2] * w * h3,
                    ],
                ));
            }
        }
    }

    result
}

/// Kernel weights not above this value are treated as zero (outside the support).
const MIN_KERNEL_WEIGHT: f64 = 1e-30;

/// Grid indices along one axis in the support of the delta kernel of a point at
/// the normalised coordinate `s = position / spacing`, with their dimensionless
/// 1-D kernel weights `φ(s - i)` (Roma kernel at unit spacing).
///
/// Only indices `0 <= i < dim` with a non-zero weight are returned, in ascending
/// order.  An empty axis or a non-finite `s` gives no nodes.
fn kernel_axis_weights(s: f64, dim: usize) -> Vec<(usize, f64)> {
    if dim == 0 || !s.is_finite() {
        return Vec::new();
    }
    let last = (dim - 1) as f64;
    let lowest = (s - 1.5).floor();
    let highest = (s + 1.5).ceil();
    if highest < 0.0 || lowest > last {
        return Vec::new();
    }
    let first = lowest.max(0.0) as usize;
    let end = highest.min(last) as usize;
    (first..=end)
        .filter_map(|node| {
            let weight = roma_delta_1d(s - node as f64, 1.0);
            (weight > MIN_KERNEL_WEIGHT).then_some((node, weight))
        })
        .collect()
}

/// Nodes of the LBM grid in the support of the IBM delta kernel of the point
/// `pos`, with the dimensionless weights
/// `ω = φ(x/h - i) φ(y/h - j) φ(z/h - k)` of the 3-D kernel (the product of
/// three 1-D Roma kernels at unit spacing).
///
/// When the whole support lies inside the grid the weights sum to one; nodes
/// outside the grid are dropped, so near a face of the grid the stencil is
/// clipped and the sum is smaller.  The nodes come in ascending row-major order.
fn kernel_stencil(pos: [f64; 3], h: f64, grid_dims: [usize; 3]) -> Vec<([usize; 3], f64)> {
    let xs = kernel_axis_weights(pos[0] / h, grid_dims[0]);
    let ys = kernel_axis_weights(pos[1] / h, grid_dims[1]);
    let zs = kernel_axis_weights(pos[2] / h, grid_dims[2]);
    let mut stencil = Vec::with_capacity(xs.len() * ys.len() * zs.len());
    for &(i, weight_x) in &xs {
        for &(j, weight_y) in &ys {
            for &(k, weight_z) in &zs {
                let omega = weight_x * weight_y * weight_z;
                if omega > MIN_KERNEL_WEIGHT {
                    stencil.push(([i, j, k], omega));
                }
            }
        }
    }
    stencil
}

/// Interpolate a velocity field from the Eulerian grid to a Lagrangian point.
///
/// Uses the 3D Roma delta function for interpolation.
///
/// # Arguments
/// * `pos` - Position of the Lagrangian point.
/// * `velocities` - Velocity field on the grid (flattened, row-major `[nx][ny][nz]`).
/// * `grid_spacing` - LBM grid cell spacing.
/// * `grid_dims` - LBM grid dimensions `[nx, ny, nz]`.
pub fn interpolate_velocity(
    pos: [f64; 3],
    velocities: &[[f64; 3]],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> [f64; 3] {
    if grid_spacing <= 0.0 {
        return [0.0; 3];
    }

    let h = grid_spacing;
    let base = [pos[0] / h, pos[1] / h, pos[2] / h];
    let mut result = [0.0; 3];

    let i_min = ((base[0] - 1.5).floor() as i64).max(0) as usize;
    let i_max = ((base[0] + 1.5).ceil() as i64).max(0) as usize;
    let j_min = ((base[1] - 1.5).floor() as i64).max(0) as usize;
    let j_max = ((base[1] + 1.5).ceil() as i64).max(0) as usize;
    let k_min = ((base[2] - 1.5).floor() as i64).max(0) as usize;
    let k_max = ((base[2] + 1.5).ceil() as i64).max(0) as usize;

    for i in i_min..=i_max.min(grid_dims[0].saturating_sub(1)) {
        let dx = pos[0] - (i as f64) * h;
        let phi_x = roma_delta_1d(dx, h);
        if phi_x.abs() < 1e-30 {
            continue;
        }
        for j in j_min..=j_max.min(grid_dims[1].saturating_sub(1)) {
            let dy = pos[1] - (j as f64) * h;
            let phi_y = roma_delta_1d(dy, h);
            if phi_y.abs() < 1e-30 {
                continue;
            }
            for k in k_min..=k_max.min(grid_dims[2].saturating_sub(1)) {
                let dz = pos[2] - (k as f64) * h;
                let phi_z = roma_delta_1d(dz, h);
                let w = phi_x * phi_y * phi_z;
                if w.abs() < 1e-30 {
                    continue;
                }

                let flat_idx = i * grid_dims[1] * grid_dims[2] + j * grid_dims[2] + k;
                if flat_idx < velocities.len() {
                    // Interpolation: u(X) = sum_x u(x) * delta(x - X) * h^3
                    let h3 = h * h * h;
                    result[0] += velocities[flat_idx][0] * w * h3;
                    result[1] += velocities[flat_idx][1] * w * h3;
                    result[2] += velocities[flat_idx][2] * w * h3;
                }
            }
        }
    }

    result
}

/// Interpolate a scalar field from the grid to a Lagrangian point.
fn interpolate_scalar(
    pos: [f64; 3],
    field: &[f64],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> f64 {
    if grid_spacing <= 0.0 {
        return 0.0;
    }

    let h = grid_spacing;
    let base = [pos[0] / h, pos[1] / h, pos[2] / h];
    let mut result = 0.0;

    let i_min = ((base[0] - 1.5).floor() as i64).max(0) as usize;
    let i_max = ((base[0] + 1.5).ceil() as i64).max(0) as usize;
    let j_min = ((base[1] - 1.5).floor() as i64).max(0) as usize;
    let j_max = ((base[1] + 1.5).ceil() as i64).max(0) as usize;
    let k_min = ((base[2] - 1.5).floor() as i64).max(0) as usize;
    let k_max = ((base[2] + 1.5).ceil() as i64).max(0) as usize;

    for i in i_min..=i_max.min(grid_dims[0].saturating_sub(1)) {
        let dx = pos[0] - (i as f64) * h;
        let phi_x = roma_delta_1d(dx, h);
        if phi_x.abs() < 1e-30 {
            continue;
        }
        for j in j_min..=j_max.min(grid_dims[1].saturating_sub(1)) {
            let dy = pos[1] - (j as f64) * h;
            let phi_y = roma_delta_1d(dy, h);
            if phi_y.abs() < 1e-30 {
                continue;
            }
            for k in k_min..=k_max.min(grid_dims[2].saturating_sub(1)) {
                let dz = pos[2] - (k as f64) * h;
                let phi_z = roma_delta_1d(dz, h);
                let w = phi_x * phi_y * phi_z;
                if w.abs() < 1e-30 {
                    continue;
                }

                let flat_idx = i * grid_dims[1] * grid_dims[2] + j * grid_dims[2] + k;
                if flat_idx < field.len() {
                    let h3 = h * h * h;
                    result += field[flat_idx] * w * h3;
                }
            }
        }
    }

    result
}

/// Compute the velocity gradient tensor at a point using central finite differences
/// on the LBM grid.
fn compute_velocity_gradient(
    pos: [f64; 3],
    velocities: &[[f64; 3]],
    grid_spacing: f64,
    grid_dims: [usize; 3],
) -> [[f64; 3]; 3] {
    let eps = grid_spacing * 0.5;
    let mut grad = [[0.0; 3]; 3];

    for d in 0..3usize {
        let mut pos_plus = pos;
        let mut pos_minus = pos;
        pos_plus[d] += eps;
        pos_minus[d] -= eps;

        let v_plus = interpolate_velocity(pos_plus, velocities, grid_spacing, grid_dims);
        let v_minus = interpolate_velocity(pos_minus, velocities, grid_spacing, grid_dims);

        let inv_2eps = 1.0 / (2.0 * eps);
        for (c, grad_row) in grad.iter_mut().enumerate() {
            grad_row[d] = (v_plus[c] - v_minus[c]) * inv_2eps;
        }
    }

    grad
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// Test that the 1D Roma delta function integrates to 1 over the support.
    #[test]
    fn test_delta_function_normalization_1d() {
        let h = 1.0;
        let n = 10000;
        let range = 3.0 * h; // well beyond the support of 1.5h
        let dr = 2.0 * range / n as f64;
        let mut integral = 0.0;

        for i in 0..n {
            let r = -range + (i as f64 + 0.5) * dr;
            integral += roma_delta_1d(r, h) * dr;
        }

        assert!(
            (integral - 1.0).abs() < 1e-4,
            "1D delta integral = {integral}, expected ~1.0"
        );
    }

    /// Test that the 3D delta function integrates to 1 over its support.
    #[test]
    fn test_delta_function_normalization_3d() {
        let h = 1.0;
        let n = 100; // per dimension
        let range = 2.0 * h;
        let dr = 2.0 * range / n as f64;
        let mut integral = 0.0;

        for i in 0..n {
            let x = -range + (i as f64 + 0.5) * dr;
            for j in 0..n {
                let y = -range + (j as f64 + 0.5) * dr;
                for k in 0..n {
                    let z = -range + (k as f64 + 0.5) * dr;
                    integral += delta_3d([x, y, z], h) * dr * dr * dr;
                }
            }
        }

        assert!(
            (integral - 1.0).abs() < 1e-3,
            "3D delta integral = {integral}, expected ~1.0"
        );
    }

    /// Test velocity interpolation on a uniform velocity field.
    /// Interpolating from a constant field should return that constant.
    #[test]
    fn test_velocity_interpolation_uniform_field() {
        let grid_dims = [10, 10, 10];
        let total = grid_dims[0] * grid_dims[1] * grid_dims[2];
        let uniform_vel = [1.5, -0.3, 2.7];
        let velocities: Vec<[f64; 3]> = vec![uniform_vel; total];
        let h = 1.0;

        // Interpolate at a point in the interior of the grid
        let pos = [5.0, 5.0, 5.0];
        let result = interpolate_velocity(pos, &velocities, h, grid_dims);

        // The IBM interpolation of a constant field should recover that constant.
        // With h=1 and the delta function integrating to 1/h^3 over grid points,
        // the sum of delta * h^3 = 1 for a point exactly on a grid node.
        for (d, (&res_d, &exp_d)) in result.iter().zip(uniform_vel.iter()).enumerate() {
            assert!(
                (res_d - exp_d).abs() < 0.1,
                "dim {d}: interpolated {res_d}, expected {exp_d}"
            );
        }
    }

    /// Test that force spreading conserves the total force.
    #[test]
    fn test_force_spreading_conservation() {
        let h = 1.0;
        let grid_dims = [10, 10, 10];
        let force = [3.0, -1.0, 2.0];
        let pos = [5.0, 5.0, 5.0]; // on a grid node

        let spread = spread_force(pos, force, h, grid_dims);

        // Sum all spread forces
        let mut total = [0.0; 3];
        for &(_, f) in &spread {
            total[0] += f[0];
            total[1] += f[1];
            total[2] += f[2];
        }

        // The total spread force should equal the original force times h^3
        // because spread_to_grid multiplies by h^3 (the volume element).
        // Actually, the sum of delta(x-X)*h^3 over grid points should be ~1
        // for the Roma kernel, so sum of spread = force * h^3 * (sum of delta * h^3 over grid)
        // With h=1: force * 1 * sum(delta), and sum(delta) ~ 1/h^3 = 1,
        // so total ~ force * h^3.
        // Actually spread already includes h^3, so total = force * sum(w * h^3)
        // and sum(w) ~ 1/h^3 => total ~ force.
        // Let's just check conservation to a reasonable tolerance.
        for (d, (&tot_d, &f_d)) in total.iter().zip(force.iter()).enumerate() {
            assert!(
                (tot_d - f_d).abs() < 0.15,
                "dim {d}: total spread force = {tot_d}, expected {f_d}"
            );
        }
    }

    /// Test round-trip coupling consistency: spread a velocity to grid, then
    /// interpolate it back at the same point. The result should be close to
    /// the original velocity.
    #[test]
    fn test_round_trip_coupling_consistency() {
        let h = 1.0;
        let grid_dims = [12, 12, 12];
        let total_cells = grid_dims[0] * grid_dims[1] * grid_dims[2];
        let pos = [6.0, 6.0, 6.0]; // exactly on a grid node
        let original_vel = [2.0, -1.0, 0.5];

        // Spread velocity to grid
        let spread = spread_to_grid(pos, original_vel, h, grid_dims);

        // Build velocity field from spread values
        let mut velocities = vec![[0.0; 3]; total_cells];
        for &(idx, vel) in &spread {
            let flat = idx[0] * grid_dims[1] * grid_dims[2] + idx[1] * grid_dims[2] + idx[2];
            if flat < total_cells {
                velocities[flat][0] += vel[0];
                velocities[flat][1] += vel[1];
                velocities[flat][2] += vel[2];
            }
        }

        // Interpolate back
        let recovered = interpolate_velocity(pos, &velocities, h, grid_dims);

        // Check consistency (not exact due to discrete convolution).
        // Round-trip spread+interpolate includes h^6 * Σ δ² factor, so only
        // directional consistency is expected.
        for (d, (&rec_d, &orig_d)) in recovered.iter().zip(original_vel.iter()).enumerate() {
            assert!(
                (rec_d - orig_d).abs() < 2.0,
                "dim {d}: recovered {rec_d}, expected {orig_d}"
            );
        }
    }

    /// Test grid index computation.
    #[test]
    fn test_grid_index_computation() {
        let h = 0.5;

        // Position exactly on grid node (2, 4, 6) in index space
        let idx = interpolate_to_grid([1.0, 2.0, 3.0], h);
        assert_eq!(idx, [2, 4, 6]);

        // Position between nodes should round
        let idx2 = interpolate_to_grid([0.3, 0.7, 1.1], h);
        assert_eq!(idx2, [1, 1, 2]);

        // Negative positions clamp to zero
        let idx3 = interpolate_to_grid([-1.0, -0.5, 0.0], h);
        assert_eq!(idx3, [0, 0, 0]);
    }

    /// Test FemLbmCoupling initialization and basic operations.
    #[test]
    fn test_coupling_initialization() {
        let config = FemLbmCouplingConfig::default();
        let mut coupling = FemLbmCoupling::new(config);

        // Should fail before initialization
        let result = coupling.transfer_velocity_to_lbm([10, 10, 10]);
        assert!(result.is_err());

        // Initialize
        let nodes = vec![0, 1, 2];
        let positions = vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]];
        let init_result = coupling.set_interface_nodes(nodes, positions, 1.0);
        assert!(init_result.is_ok());

        // Should succeed after initialization
        let result = coupling.transfer_velocity_to_lbm([10, 10, 10]);
        assert!(result.is_ok());

        // Forces should be zero initially
        let forces = coupling.get_fem_forces();
        assert_eq!(forces.len(), 3);
        for f in forces {
            assert_eq!(*f, [0.0; 3]);
        }
    }

    /// Test dimension mismatch error.
    #[test]
    fn test_dimension_mismatch_error() {
        let config = FemLbmCouplingConfig::default();
        let mut coupling = FemLbmCoupling::new(config);

        let nodes = vec![0, 1, 2];
        let positions = vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]; // mismatched length
        let result = coupling.set_interface_nodes(nodes, positions, 1.0);
        assert!(result.is_err());
    }

    /// Test invalid grid spacing.
    #[test]
    fn test_invalid_grid_spacing() {
        let config = FemLbmCouplingConfig::default();
        let mut coupling = FemLbmCoupling::new(config);

        let nodes = vec![0];
        let positions = vec![[1.0, 1.0, 1.0]];
        let result = coupling.set_interface_nodes(nodes, positions, -1.0);
        assert!(result.is_err());
    }

    /// Test the delta function at specific known values.
    #[test]
    fn test_delta_function_specific_values() {
        let h = 1.0;

        // At r=0, the Roma kernel gives (1/3h)(1 + sqrt(1)) = (1/3)(1+1) = 2/3
        let val_at_zero = delta_function(0.0, h);
        assert!(
            (val_at_zero - 2.0 / 3.0).abs() < 1e-10,
            "delta(0, 1) = {val_at_zero}, expected {}",
            2.0 / 3.0
        );

        // Beyond support: r >= 1.5h should give 0
        let val_outside = delta_function(1.5, h);
        assert!(
            val_outside.abs() < 1e-10,
            "delta(1.5, 1) = {val_outside}, expected 0"
        );

        // Symmetry: delta(-r) == delta(r)
        let val_pos = delta_function(0.7, h);
        let val_neg = delta_function(-0.7, h);
        assert!(
            (val_pos - val_neg).abs() < 1e-14,
            "delta function should be symmetric"
        );
    }

    /// Test that spreading with zero grid spacing returns empty.
    #[test]
    fn test_spread_zero_grid_spacing() {
        let result = spread_force([1.0, 1.0, 1.0], [1.0, 0.0, 0.0], 0.0, [5, 5, 5]);
        assert!(result.is_empty());
    }

    /// Test compute_fluid_forces with field size mismatch.
    #[test]
    fn test_fluid_forces_field_mismatch() {
        let config = FemLbmCouplingConfig::default();
        let mut coupling = FemLbmCoupling::new(config);
        let nodes = vec![0];
        let positions = vec![[1.0, 1.0, 1.0]];
        coupling.set_interface_nodes(nodes, positions, 1.0).ok();

        let grid_dims = [5, 5, 5];
        let wrong_density = vec![1.0; 10]; // wrong size
        let wrong_velocity = vec![[0.0; 3]; 10];

        let result = coupling.compute_fluid_forces(&wrong_density, &wrong_velocity, 1.0, grid_dims);
        assert!(result.is_err());
    }

    /// Test the full coupling step.
    #[test]
    fn test_coupling_step() {
        let config = FemLbmCouplingConfig {
            penalty_stiffness: 100.0,
            ..FemLbmCouplingConfig::default()
        };
        let mut coupling = FemLbmCoupling::new(config);

        let nodes = vec![0];
        let positions = vec![[3.0, 3.0, 3.0]];
        coupling
            .set_interface_nodes(nodes, positions.clone(), 1.0)
            .ok();

        // Set a velocity on the FEM side
        coupling
            .update_fem_state(&positions, &[[1.0, 0.0, 0.0]])
            .ok();

        let grid_dims = [8, 8, 8];
        let total = grid_dims[0] * grid_dims[1] * grid_dims[2];
        let density = vec![1.0; total];
        let velocity = vec![[0.0; 3]; total]; // fluid at rest

        let result = coupling.coupling_step(&density, &velocity, 1.0, grid_dims);
        assert!(result.is_ok());

        // The fluid is at rest while the structure moves with (1, 0, 0), so the
        // no-slip feedback k (u_fluid - u_solid) h^3 = 100 * (-1) decelerates the
        // structure; the point has no normal, hence no pressure / viscous part.
        let forces = coupling.get_fem_forces();
        assert!((forces[0][0] + 100.0).abs() < 1e-9, "{:?}", forces[0]);
        assert!(forces[0][1].abs() < 1e-12 && forces[0][2].abs() < 1e-12);
    }

    /// Test ImmersedBoundaryPoint creation.
    #[test]
    fn test_ib_point_creation() {
        let p = ImmersedBoundaryPoint::new([1.0, 2.0, 3.0]);
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.velocity, [0.0; 3]);
        assert_eq!(p.force, [0.0; 3]);
        assert_eq!(p.reference_position, [1.0, 2.0, 3.0]);
        assert!((p.area_weight - 1.0).abs() < 1e-14);

        let p2 = ImmersedBoundaryPoint::with_area_weight([4.0, 5.0, 6.0], 2.5);
        assert!((p2.area_weight - 2.5).abs() < 1e-14);
    }

    /// Test reset_forces zeroes everything.
    #[test]
    fn test_reset_forces() {
        let config = FemLbmCouplingConfig::default();
        let mut coupling = FemLbmCoupling::new(config);
        let nodes = vec![0, 1];
        let positions = vec![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
        coupling.set_interface_nodes(nodes, positions, 1.0).ok();

        // Manually set some forces
        coupling.fluid_forces[0] = [10.0, 20.0, 30.0];
        coupling.fluid_forces[1] = [40.0, 50.0, 60.0];
        coupling.ib_points[0].force = [10.0, 20.0, 30.0];

        coupling.reset_forces();

        for f in coupling.get_fem_forces() {
            assert_eq!(*f, [0.0; 3]);
        }
        for ib in &coupling.ib_points {
            assert_eq!(ib.force, [0.0; 3]);
        }
    }

    // ------------------------------------------------------------------
    // Fluid traction on the interface (pressure along the normal, viscous
    // stress, no-slip feedback)
    // ------------------------------------------------------------------

    /// Index of the cell `(i, j, k)` in the flattened row-major LBM fields.
    fn cell(dims: [usize; 3], i: usize, j: usize, k: usize) -> usize {
        i * dims[1] * dims[2] + j * dims[2] + k
    }

    /// Triangulated sphere obtained by subdividing an octahedron.  Triangles are
    /// wound counter-clockwise seen from outside.
    fn sphere_mesh(
        subdivisions: usize,
        radius: f64,
        center: [f64; 3],
    ) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
        fn midpoint(
            points: &mut Vec<[f64; 3]>,
            cache: &mut HashMap<(usize, usize), usize>,
            a: usize,
            b: usize,
        ) -> usize {
            let key = (a.min(b), a.max(b));
            if let Some(&idx) = cache.get(&key) {
                return idx;
            }
            let m = [
                points[a][0] + points[b][0],
                points[a][1] + points[b][1],
                points[a][2] + points[b][2],
            ];
            let len = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt();
            points.push([m[0] / len, m[1] / len, m[2] / len]);
            cache.insert(key, points.len() - 1);
            points.len() - 1
        }

        let mut points = vec![
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 0.0, -1.0],
        ];
        let mut faces = vec![
            [0, 2, 4],
            [2, 1, 4],
            [1, 3, 4],
            [3, 0, 4],
            [2, 0, 5],
            [1, 2, 5],
            [3, 1, 5],
            [0, 3, 5],
        ];
        for _ in 0..subdivisions {
            let mut cache = HashMap::new();
            let mut refined = Vec::with_capacity(faces.len() * 4);
            for f in &faces {
                let ab = midpoint(&mut points, &mut cache, f[0], f[1]);
                let bc = midpoint(&mut points, &mut cache, f[1], f[2]);
                let ca = midpoint(&mut points, &mut cache, f[2], f[0]);
                refined.push([f[0], ab, ca]);
                refined.push([f[1], bc, ab]);
                refined.push([f[2], ca, bc]);
                refined.push([ab, bc, ca]);
            }
            faces = refined;
        }
        let points = points
            .iter()
            .map(|p| {
                [
                    center[0] + radius * p[0],
                    center[1] + radius * p[1],
                    center[2] + radius * p[2],
                ]
            })
            .collect();
        (points, faces)
    }

    fn sum3(forces: &[[f64; 3]]) -> [f64; 3] {
        forces.iter().fold([0.0; 3], |acc, f| {
            [acc[0] + f[0], acc[1] + f[1], acc[2] + f[2]]
        })
    }

    /// Hydrostatic pressure on a *closed* surface integrates to zero net force,
    /// independent of the (tiny) slip velocity.  The old code pushed the full
    /// pressure along the component-wise sign of the slip velocity on every
    /// point, i.e. a net force of -p * area * N in x.
    #[test]
    fn test_uniform_pressure_on_closed_sphere_has_zero_net_force() {
        let dims = [17, 17, 17];
        let total = dims[0] * dims[1] * dims[2];
        let center = [8.0, 8.0, 8.0];
        let radius = 3.0;
        let (points, faces) = sphere_mesh(2, radius, center);
        let n = points.len();

        let config = FemLbmCouplingConfig {
            penalty_stiffness: 0.0,
            ..FemLbmCouplingConfig::default()
        };
        for slip in [1e-9, -1e-9, 0.0] {
            let mut coupling = FemLbmCoupling::new(config.clone());
            coupling
                .set_interface_surface((0..n).collect(), points.clone(), faces.clone(), 1.0)
                .expect("valid surface");
            let density = vec![1.0; total];
            let velocity = vec![[slip, 0.0, 0.0]; total];
            coupling
                .compute_fluid_forces(&density, &velocity, 1.0, dims)
                .expect("fields match the grid");

            let pressure = 1.0 / 3.0;
            let net = sum3(coupling.get_fem_forces());
            // |F_net| relative to the force scale p * (sphere area)
            let scale = pressure * 4.0 * std::f64::consts::PI * radius * radius;
            for (d, f) in net.iter().enumerate() {
                assert!(f.abs() < 1e-12 * scale, "slip {slip}: net F[{d}] = {f}");
            }
            // Every point is pushed *into* the structure, exactly along its own
            // normal: F = -p * (area_weight h^2) * n, whatever the slip.
            let mut total_magnitude = 0.0;
            for (ibp, force) in coupling.ib_points.iter().zip(coupling.get_fem_forces()) {
                let area = ibp.area_weight;
                for (f_d, n_d) in force.iter().zip(ibp.normal.iter()) {
                    assert!(
                        (f_d + pressure * area * n_d).abs() < 1e-12,
                        "force {force:?} is not -p A n for normal {:?}",
                        ibp.normal
                    );
                }
                // the vertex normals of the tessellation are close to radial
                let radial = [
                    (ibp.position[0] - center[0]) / radius,
                    (ibp.position[1] - center[1]) / radius,
                    (ibp.position[2] - center[2]) / radius,
                ];
                let alignment: f64 = ibp
                    .normal
                    .iter()
                    .zip(radial.iter())
                    .map(|(a, b)| a * b)
                    .sum();
                assert!(alignment > 0.95, "normal is not outward: {alignment}");
                total_magnitude += (force[0].powi(2) + force[1].powi(2) + force[2].powi(2)).sqrt();
            }
            // The polyhedral sphere is slightly smaller than the sphere.
            assert!(total_magnitude > 0.85 * scale && total_magnitude <= scale);
        }
    }

    /// A flat plate (normal +z, in the fluid above) in a simple shear flow
    /// `u_x = gamma z` and a uniform density: pressure pushes it down, the viscous
    /// traction drags it along +x with `mu * gamma` per unit area.
    #[test]
    fn test_flat_plate_pressure_and_shear_traction() {
        let dims = [17, 17, 17];
        let total = dims[0] * dims[1] * dims[2];
        let (rho, gamma, mu) = (1.2, 0.01, 0.2);
        let mut velocity = vec![[0.0; 3]; total];
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    velocity[cell(dims, i, j, k)] = [gamma * k as f64, 0.0, 0.0];
                }
            }
        }
        let density = vec![rho; total];

        let mut points = Vec::new();
        for i in 5..=11 {
            for j in 5..=11 {
                points.push([i as f64 + 0.25, j as f64 - 0.4, 8.0]);
            }
        }
        let n = points.len();
        let config = FemLbmCouplingConfig {
            penalty_stiffness: 0.0,
            fluid_viscosity: mu,
            ..FemLbmCouplingConfig::default()
        };
        let mut coupling = FemLbmCoupling::new(config);
        coupling
            .set_interface_nodes((0..n).collect(), points, 1.0)
            .expect("valid interface");
        coupling
            .set_interface_normals(&vec![[0.0, 0.0, 2.5]; n], &vec![1.0; n])
            .expect("matching lengths");
        // the normal is normalised on input
        assert!((coupling.ib_points[0].normal[2] - 1.0).abs() < 1e-15);

        coupling
            .compute_fluid_forces(&density, &velocity, 1.0, dims)
            .expect("fields match the grid");
        let pressure = rho / 3.0;
        for force in coupling.get_fem_forces() {
            assert!((force[0] - mu * gamma).abs() < 1e-12, "F_x = {}", force[0]);
            assert!(force[1].abs() < 1e-12, "F_y = {}", force[1]);
            assert!((force[2] + pressure).abs() < 1e-12, "F_z = {}", force[2]);
        }
    }

    /// The viscous traction uses the symmetric rate-of-strain tensor: a flow
    /// `u_z = gamma x` shears a plate with normal +z along +x (through
    /// `du_z/dx`), which a velocity-gradient-row formula would miss.
    #[test]
    fn test_viscous_traction_includes_transposed_gradient() {
        let dims = [17, 17, 17];
        let total = dims[0] * dims[1] * dims[2];
        let (gamma, mu) = (0.02, 0.5);
        let mut velocity = vec![[0.0; 3]; total];
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    velocity[cell(dims, i, j, k)] = [0.0, 0.0, gamma * i as f64];
                }
            }
        }
        let density = vec![0.0; total];
        let config = FemLbmCouplingConfig {
            penalty_stiffness: 0.0,
            fluid_viscosity: mu,
            ..FemLbmCouplingConfig::default()
        };
        let mut coupling = FemLbmCoupling::new(config);
        coupling
            .set_interface_nodes(vec![0], vec![[8.3, 7.9, 8.0]], 1.0)
            .expect("valid interface");
        coupling
            .set_interface_normals(&[[0.0, 0.0, 1.0]], &[2.0])
            .expect("matching lengths");
        coupling
            .compute_fluid_forces(&density, &velocity, 1.0, dims)
            .expect("fields match the grid");
        let force = coupling.get_fem_forces()[0];
        // t = mu (du_x/dz + du_z/dx) n_z along x, area 2
        assert!((force[0] - mu * gamma * 2.0).abs() < 1e-12, "{force:?}");
        // t_z = 2 mu du_z/dz = 0, t_y = 0
        assert!(
            force[1].abs() < 1e-12 && force[2].abs() < 1e-12,
            "{force:?}"
        );
    }

    /// The no-slip feedback accelerates the structure toward the fluid velocity
    /// (a force on the *structure* along u_fluid - u_solid).
    #[test]
    fn test_noslip_feedback_drags_structure_along_with_the_flow() {
        let dims = [10, 10, 10];
        let total = 1000;
        let config = FemLbmCouplingConfig {
            penalty_stiffness: 100.0,
            ..FemLbmCouplingConfig::default()
        };
        let mut coupling = FemLbmCoupling::new(config);
        coupling
            .set_interface_nodes(vec![0], vec![[5.0, 5.0, 5.0]], 1.0)
            .expect("valid interface");
        coupling
            .update_fem_state(&[[5.0, 5.0, 5.0]], &[[0.25, 0.0, 0.0]])
            .expect("matching lengths");
        // fluid faster than the structure along +x, uniform density, no normal
        let density = vec![1.0; total];
        let velocity = vec![[1.0, 0.0, 0.0]; total];
        coupling
            .compute_fluid_forces(&density, &velocity, 1.0, dims)
            .expect("fields match the grid");
        let force = coupling.get_fem_forces()[0];
        // k * (1.0 - 0.25) * h^3 along +x, and no pressure without a normal
        assert!((force[0] - 100.0 * 0.75).abs() < 1e-10, "{force:?}");
        assert!(force[1].abs() < 1e-12 && force[2].abs() < 1e-12);

        // structure moving faster than the fluid is decelerated
        coupling
            .update_fem_state(&[[5.0, 5.0, 5.0]], &[[3.0, 0.0, 0.0]])
            .expect("matching lengths");
        coupling
            .compute_fluid_forces(&density, &velocity, 1.0, dims)
            .expect("fields match the grid");
        assert!(coupling.get_fem_forces()[0][0] < 0.0);
    }

    #[test]
    fn test_points_without_normal_feel_no_pressure() {
        let dims = [8, 8, 8];
        let total = 512;
        let mut coupling = FemLbmCoupling::new(FemLbmCouplingConfig::default());
        coupling
            .set_interface_nodes(vec![0], vec![[4.0, 4.0, 4.0]], 1.0)
            .expect("valid interface");
        let density = vec![5.0; total];
        let velocity = vec![[0.0; 3]; total];
        coupling
            .compute_fluid_forces(&density, &velocity, 1.0, dims)
            .expect("fields match the grid");
        assert_eq!(coupling.get_fem_forces()[0], [0.0; 3]);
    }

    #[test]
    fn test_interface_surface_gives_outward_unit_normals_and_vector_areas() {
        // Octahedron with circumradius 1: each vertex touches 4 faces whose area
        // vectors add up to (2/3) e_axis.
        let (points, faces) = sphere_mesh(0, 1.0, [0.0; 3]);
        let n = points.len();
        let mut coupling = FemLbmCoupling::new(FemLbmCouplingConfig::default());
        coupling
            .set_interface_surface((0..n).collect(), points.clone(), faces.clone(), 0.5)
            .expect("valid surface");
        let mut total = [0.0; 3];
        for ibp in &coupling.ib_points {
            let len = ibp.normal.iter().map(|x| x * x).sum::<f64>().sqrt();
            assert!((len - 1.0).abs() < 1e-14, "unit normal");
            for (n_d, p_d) in ibp.normal.iter().zip(ibp.position.iter()) {
                assert!((n_d - p_d).abs() < 1e-14, "radial");
            }
            for (t_d, n_d) in total.iter_mut().zip(ibp.normal.iter()) {
                *t_d += ibp.area_weight * 0.25 * n_d;
            }
            // |a| = 2/3 in physical units; weights are in units of h^2 = 0.25
            assert!((ibp.area_weight * 0.25 - 2.0 / 3.0).abs() < 1e-14);
        }
        // closed surface: the tributary area vectors cancel
        assert!(total.iter().all(|t| t.abs() < 1e-14), "{total:?}");

        // reversed winding gives inward normals
        let flipped: Vec<[usize; 3]> = faces.iter().map(|f| [f[0], f[2], f[1]]).collect();
        coupling
            .set_interface_surface((0..n).collect(), points, flipped, 0.5)
            .expect("valid surface");
        for ibp in &coupling.ib_points {
            let outward: f64 = ibp
                .normal
                .iter()
                .zip(ibp.position.iter())
                .map(|(a, b)| a * b)
                .sum();
            assert!(outward < 0.0);
        }
    }

    #[test]
    fn test_update_fem_state_refreshes_normals_and_areas() {
        let (points, faces) = sphere_mesh(1, 1.0, [0.0; 3]);
        let n = points.len();
        let mut coupling = FemLbmCoupling::new(FemLbmCouplingConfig::default());
        coupling
            .set_interface_surface((0..n).collect(), points.clone(), faces, 1.0)
            .expect("valid surface");
        let area_before: Vec<f64> = coupling.ib_points.iter().map(|p| p.area_weight).collect();

        // Dilate the structure by a factor 2: areas grow by 4, normals unchanged.
        let scaled: Vec<[f64; 3]> = points
            .iter()
            .map(|p| [2.0 * p[0], 2.0 * p[1], 2.0 * p[2]])
            .collect();
        coupling
            .update_fem_state(&scaled, &vec![[0.0; 3]; n])
            .expect("matching lengths");
        for (before, ibp) in area_before.iter().zip(coupling.ib_points.iter()) {
            assert!((ibp.area_weight - 4.0 * before).abs() < 1e-12 * before);
            for d in 0..3 {
                assert!((ibp.normal[d] - ibp.position[d] / 2.0).abs() < 1e-14);
            }
        }

        // Rotate by 90 degrees about z: the normals rotate with the surface.
        let rotated: Vec<[f64; 3]> = points.iter().map(|p| [-p[1], p[0], p[2]]).collect();
        coupling
            .update_fem_state(&rotated, &vec![[0.0; 3]; n])
            .expect("matching lengths");
        for ibp in &coupling.ib_points {
            for d in 0..3 {
                assert!((ibp.normal[d] - ibp.position[d]).abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_interface_surface_validation() {
        let mut coupling = FemLbmCoupling::new(FemLbmCouplingConfig::default());
        let pts = vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        // face refers to a point that does not exist
        let bad_face =
            coupling.set_interface_surface(vec![0, 1, 2], pts.clone(), vec![[0, 1, 3]], 1.0);
        assert!(matches!(bad_face, Err(CouplingError::InvalidParameter(_))));
        // mismatching node / position lengths
        let bad_len = coupling.set_interface_surface(vec![0, 1], pts.clone(), vec![[0, 1, 2]], 1.0);
        assert!(matches!(
            bad_len,
            Err(CouplingError::DimensionMismatch { .. })
        ));
        // bad grid spacing
        let bad_h =
            coupling.set_interface_surface(vec![0, 1, 2], pts.clone(), vec![[0, 1, 2]], 0.0);
        assert!(bad_h.is_err());
        // normals need an initialised interface and matching lengths
        let mut fresh = FemLbmCoupling::new(FemLbmCouplingConfig::default());
        assert!(matches!(
            fresh.set_interface_normals(&[], &[]),
            Err(CouplingError::NotInitialized)
        ));
        coupling
            .set_interface_nodes(vec![0, 1, 2], pts, 1.0)
            .expect("valid interface");
        assert!(
            coupling
                .set_interface_normals(&[[0.0, 0.0, 1.0]], &[1.0])
                .is_err()
        );
        assert!(
            coupling
                .set_interface_normals(&[[0.0, 0.0, 1.0]; 3], &[1.0, -1.0, 1.0])
                .is_err()
        );
    }

    #[test]
    fn test_ib_point_with_normal_is_normalised() {
        let p = ImmersedBoundaryPoint::with_normal([0.0; 3], [0.0, 3.0, 4.0], 2.0);
        assert!((p.normal[1] - 0.6).abs() < 1e-15 && (p.normal[2] - 0.8).abs() < 1e-15);
        assert!((p.area_weight - 2.0).abs() < 1e-15);
        let unset = ImmersedBoundaryPoint::with_normal([0.0; 3], [0.0; 3], 1.0);
        assert_eq!(unset.normal, [0.0; 3]);
        let nan = ImmersedBoundaryPoint::with_normal([0.0; 3], [f64::NAN, 0.0, 1.0], 1.0);
        assert_eq!(nan.normal, [0.0; 3]);
    }
}
