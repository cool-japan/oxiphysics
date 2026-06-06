//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;
/// MPM particle material type.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaterialType {
    /// Elastic neo-Hookean material.
    NeoHookean,
    /// Elastic-perfectly-plastic material with von Mises yield criterion.
    VonMises,
    /// Snow material with hardening (Stomakhin et al. 2013).
    Snow,
    /// Sand material with Drucker-Prager yield criterion (Klar et al. 2016).
    Sand,
    /// Fluid (weakly compressible or pressure projection).
    Fluid,
    /// Rigid body (kinematic, not updated by MPM).
    Rigid,
}
/// A uniform background Eulerian grid for MPM.
#[derive(Debug, Clone)]
pub struct MpmGrid {
    /// Number of nodes in each dimension: (nx, ny, nz).
    pub dims: [usize; 3],
    /// Grid spacing (cell size).
    pub h: f64,
    /// Grid origin (min corner position).
    pub origin: [f64; 3],
    /// Flat array of grid nodes.
    pub nodes: Vec<GridNode>,
}
impl MpmGrid {
    /// Create a new uniform grid.
    pub fn new(dims: [usize; 3], h: f64, origin: [f64; 3]) -> Self {
        let total = dims[0] * dims[1] * dims[2];
        Self {
            dims,
            h,
            origin,
            nodes: vec![GridNode::default(); total],
        }
    }
    /// Flat index from (ix, iy, iz).
    #[inline]
    pub fn idx(&self, ix: usize, iy: usize, iz: usize) -> usize {
        ix * self.dims[1] * self.dims[2] + iy * self.dims[2] + iz
    }
    /// World-space position of node (ix, iy, iz).
    #[inline]
    pub fn node_pos(&self, ix: usize, iy: usize, iz: usize) -> [f64; 3] {
        [
            self.origin[0] + ix as f64 * self.h,
            self.origin[1] + iy as f64 * self.h,
            self.origin[2] + iz as f64 * self.h,
        ]
    }
    /// Base grid index and fractional offset for a given world position.
    #[inline]
    pub fn base_and_frac(&self, pos: [f64; 3]) -> ([i64; 3], [f64; 3]) {
        let mut base = [0i64; 3];
        let mut frac = [0.0f64; 3];
        for d in 0..3 {
            let x = (pos[d] - self.origin[d]) / self.h;
            base[d] = x.floor() as i64;
            frac[d] = x - base[d] as f64;
        }
        (base, frac)
    }
    /// Reset all grid nodes.
    pub fn reset(&mut self) {
        for node in &mut self.nodes {
            node.reset();
        }
    }
    /// Total number of nodes.
    #[inline]
    pub fn num_nodes(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }
    /// Check whether a grid index is in bounds.
    #[inline]
    pub fn in_bounds(&self, ix: i64, iy: i64, iz: i64) -> bool {
        ix >= 0
            && iy >= 0
            && iz >= 0
            && (ix as usize) < self.dims[0]
            && (iy as usize) < self.dims[1]
            && (iz as usize) < self.dims[2]
    }
    /// Apply Dirichlet (no-slip) boundary conditions on the boundary layer.
    pub fn apply_boundary_conditions(&mut self, thickness: usize) {
        let [nx, ny, nz] = self.dims;
        let t = thickness;
        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let boundary =
                        ix < t || ix >= nx - t || iy < t || iy >= ny - t || iz < t || iz >= nz - t;
                    if boundary {
                        let idx = self.idx(ix, iy, iz);
                        self.nodes[idx].velocity = [0.0; 3];
                        self.nodes[idx].is_boundary = true;
                    }
                }
            }
        }
    }
    /// Sticky floor boundary: zero velocity at bottom layer.
    pub fn apply_floor_bc(&mut self, floor_layer: usize) {
        let [nx, ny, nz] = self.dims;
        for ix in 0..nx {
            for iy in 0..floor_layer.min(ny) {
                for iz in 0..nz {
                    let idx = self.idx(ix, iy, iz);
                    self.nodes[idx].velocity = [0.0; 3];
                }
            }
        }
    }
    /// Slip floor: zero normal (y) velocity at bottom layer.
    pub fn apply_slip_floor_bc(&mut self, floor_layer: usize) {
        let [nx, ny, nz] = self.dims;
        for ix in 0..nx {
            for iy in 0..floor_layer.min(ny) {
                for iz in 0..nz {
                    let idx = self.idx(ix, iy, iz);
                    self.nodes[idx].velocity[1] = self.nodes[idx].velocity[1].max(0.0);
                }
            }
        }
    }
}
/// Snow simulation parameters.
#[derive(Debug, Clone)]
pub struct SnowParams {
    /// Initial Lamé μ.
    pub mu0: f64,
    /// Initial Lamé λ.
    pub lambda0: f64,
    /// Hardening coefficient.
    pub hardening: f64,
    /// Critical compression.
    pub theta_c: f64,
    /// Critical stretch.
    pub theta_s: f64,
    /// Initial density.
    pub rho0: f64,
}
/// A single background Eulerian grid node.
#[derive(Debug, Clone, Default)]
pub struct GridNode {
    /// Accumulated mass.
    pub mass: f64,
    /// Accumulated momentum (pre-normalized to velocity).
    pub momentum: [f64; 3],
    /// Grid velocity (after dividing by mass).
    pub velocity: [f64; 3],
    /// Force accumulated from particles.
    pub force: [f64; 3],
    /// Whether this node is on a fixed (Dirichlet) boundary.
    pub is_boundary: bool,
    /// Material ID (for multi-material).
    pub material_id: usize,
}
impl GridNode {
    /// Normalize momentum to velocity: v = p / m.
    #[inline]
    pub fn normalize_velocity(&mut self) {
        if self.mass > 1e-15 {
            self.velocity = scale3(self.momentum, 1.0 / self.mass);
        } else {
            self.velocity = [0.0; 3];
        }
    }
    /// Apply force impulse: v += f * dt / m.
    #[inline]
    pub fn apply_force(&mut self, dt: f64) {
        if self.mass > 1e-15 {
            let accel = scale3(self.force, 1.0 / self.mass);
            self.velocity = add3(self.velocity, scale3(accel, dt));
        }
    }
    /// Reset all accumulated quantities for a new step.
    #[inline]
    pub fn reset(&mut self) {
        self.mass = 0.0;
        self.momentum = [0.0; 3];
        self.velocity = [0.0; 3];
        self.force = [0.0; 3];
    }
}
/// A single MPM particle carrying position, velocity, mass,
/// deformation gradient, and APIC affine velocity field.
#[derive(Debug, Clone)]
pub struct MpmParticle {
    /// World-space position.
    pub position: [f64; 3],
    /// World-space velocity.
    pub velocity: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// Initial volume (J0 = det(F0) = 1 at creation).
    pub volume0: f64,
    /// Deformation gradient F.
    pub deformation_gradient: [[f64; 3]; 3],
    /// Elastic part of the deformation gradient (multiplicative decomp).
    pub fe: [[f64; 3]; 3],
    /// Plastic part of the deformation gradient.
    pub fp: [[f64; 3]; 3],
    /// APIC affine velocity matrix C (B-matrix in MLS-MPM).
    pub affine_c: [[f64; 3]; 3],
    /// Material type.
    pub material: MaterialType,
    /// Hardening coefficient (for snow/sand).
    pub hardening: f64,
    /// Lamé parameter μ₀.
    pub mu0: f64,
    /// Lamé parameter λ₀.
    pub lambda0: f64,
    /// Yield stress (von Mises) or friction angle (Drucker-Prager).
    pub yield_param: f64,
    /// Level-set signed distance value at this particle.
    pub level_set_phi: f64,
    /// Material ID for multi-material MPM.
    pub material_id: usize,
    /// Whether this particle is active (not culled).
    pub active: bool,
}
impl MpmParticle {
    /// Create a new MPM particle with given position, velocity, and mass.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        volume0: f64,
        material: MaterialType,
        mu0: f64,
        lambda0: f64,
    ) -> Self {
        Self {
            position,
            velocity,
            mass,
            volume0,
            deformation_gradient: mat3_identity(),
            fe: mat3_identity(),
            fp: mat3_identity(),
            affine_c: mat3_zero(),
            material,
            hardening: 10.0,
            mu0,
            lambda0,
            yield_param: 1e4,
            level_set_phi: f64::MAX,
            material_id: 0,
            active: true,
        }
    }
    /// Compute the current volume: V = V₀ * det(F).
    #[inline]
    pub fn volume(&self) -> f64 {
        self.volume0 * mat3_det(self.deformation_gradient).abs()
    }
    /// Compute the current density: ρ = m / V.
    #[inline]
    pub fn density(&self) -> f64 {
        let v = self.volume();
        if v < 1e-15 { 0.0 } else { self.mass / v }
    }
    /// Compute the Jacobian J = det(F).
    #[inline]
    pub fn jacobian(&self) -> f64 {
        mat3_det(self.deformation_gradient)
    }
}
/// High-level MPM simulation context.
pub struct MpmSimulation {
    /// The background grid.
    pub grid: MpmGrid,
    /// All particles.
    pub particles: Vec<MpmParticle>,
    /// Gravity acceleration.
    pub gravity: [f64; 3],
    /// Snow-specific parameters (used if snow particles exist).
    pub snow_params: SnowParams,
    /// Sand-specific parameters.
    pub sand_params: SandParams,
    /// Current simulation time.
    pub time: f64,
    /// Step counter.
    pub step: u64,
    /// CFL safety factor.
    pub cfl: f64,
    /// Fixed time step (if > 0), otherwise adaptive.
    pub fixed_dt: f64,
    /// Boundary thickness (number of grid layers for Dirichlet BC).
    pub boundary_thickness: usize,
}
impl MpmSimulation {
    /// Create a new MPM simulation.
    pub fn new(grid_dims: [usize; 3], h: f64, grid_origin: [f64; 3], gravity: [f64; 3]) -> Self {
        Self {
            grid: MpmGrid::new(grid_dims, h, grid_origin),
            particles: Vec::new(),
            gravity,
            snow_params: SnowParams::default(),
            sand_params: SandParams::default(),
            time: 0.0,
            step: 0,
            cfl: 0.4,
            fixed_dt: 0.0,
            boundary_thickness: 2,
        }
    }
    /// Add particles to the simulation.
    pub fn add_particles(&mut self, particles: Vec<MpmParticle>) {
        self.particles.extend(particles);
    }
    /// Get the number of active particles.
    pub fn num_active_particles(&self) -> usize {
        self.particles.iter().filter(|p| p.active).count()
    }
    /// Advance the simulation by one step.
    pub fn step(&mut self, dt: f64) {
        let theta_c = self.snow_params.theta_c;
        let theta_s = self.snow_params.theta_s;
        let bt = self.boundary_thickness;
        mpm_step(
            &mut self.grid,
            &mut self.particles,
            dt,
            self.gravity,
            theta_c,
            theta_s,
            bt,
        );
        self.time += dt;
        self.step += 1;
    }
    /// Compute total kinetic energy.
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| p.active)
            .map(|p| 0.5 * p.mass * dot3(p.velocity, p.velocity))
            .sum()
    }
    /// Compute total momentum.
    pub fn total_momentum(&self) -> [f64; 3] {
        self.particles
            .iter()
            .filter(|p| p.active)
            .fold([0.0f64; 3], |acc, p| add3(acc, scale3(p.velocity, p.mass)))
    }
    /// Compute total mass.
    pub fn total_mass(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| p.active)
            .map(|p| p.mass)
            .sum()
    }
    /// Remove particles that have left the grid domain.
    pub fn cull_out_of_bounds_particles(&mut self) {
        let origin = self.grid.origin;
        let dims = self.grid.dims;
        let h = self.grid.h;
        let max_x = origin[0] + dims[0] as f64 * h;
        let max_y = origin[1] + dims[1] as f64 * h;
        let max_z = origin[2] + dims[2] as f64 * h;
        for p in &mut self.particles {
            if p.position[0] < origin[0]
                || p.position[0] > max_x
                || p.position[1] < origin[1]
                || p.position[1] > max_y
                || p.position[2] < origin[2]
                || p.position[2] > max_z
            {
                p.active = false;
            }
        }
    }
}
/// Sand simulation parameters.
#[derive(Debug, Clone)]
pub struct SandParams {
    /// Initial Lamé μ.
    pub mu0: f64,
    /// Initial Lamé λ.
    pub lambda0: f64,
    /// Friction angle in degrees.
    pub friction_angle: f64,
    /// Initial density.
    pub rho0: f64,
}
/// Multi-material MPM grid: one velocity field per material.
#[derive(Debug, Clone)]
pub struct MultiMaterialGrid {
    /// Base grid geometry.
    pub grid: MpmGrid,
    /// Number of materials.
    pub num_materials: usize,
    /// Per-material velocity fields (flat: material_id * num_nodes + node_idx).
    pub velocities: Vec<[f64; 3]>,
    /// Per-material mass fields.
    pub masses: Vec<f64>,
}
impl MultiMaterialGrid {
    /// Create a new multi-material grid.
    pub fn new(dims: [usize; 3], h: f64, origin: [f64; 3], num_materials: usize) -> Self {
        let n = dims[0] * dims[1] * dims[2];
        Self {
            grid: MpmGrid::new(dims, h, origin),
            num_materials,
            velocities: vec![[0.0; 3]; num_materials * n],
            masses: vec![0.0; num_materials * n],
        }
    }
    /// Reset all per-material fields.
    pub fn reset(&mut self) {
        for v in &mut self.velocities {
            *v = [0.0; 3];
        }
        for m in &mut self.masses {
            *m = 0.0;
        }
        self.grid.reset();
    }
    /// Get per-material velocity for node `node_idx` and material `mat_id`.
    #[inline]
    pub fn material_velocity(&self, mat_id: usize, node_idx: usize) -> [f64; 3] {
        let n = self.grid.num_nodes();
        self.velocities[mat_id * n + node_idx]
    }
    /// Get per-material mass for node `node_idx` and material `mat_id`.
    #[inline]
    pub fn material_mass(&self, mat_id: usize, node_idx: usize) -> f64 {
        let n = self.grid.num_nodes();
        self.masses[mat_id * n + node_idx]
    }
    /// Apply no-slip multi-material coupling at shared nodes.
    pub fn apply_contact_coupling(&mut self, friction_coeff: f64) {
        let n = self.grid.num_nodes();
        for node_idx in 0..n {
            let mut total_mass = 0.0;
            let mut total_momentum = [0.0f64; 3];
            for mat_id in 0..self.num_materials {
                let m = self.material_mass(mat_id, node_idx);
                let v = self.material_velocity(mat_id, node_idx);
                total_mass += m;
                total_momentum = add3(total_momentum, scale3(v, m));
            }
            if total_mass < 1e-15 {
                continue;
            }
            let v_cm = scale3(total_momentum, 1.0 / total_mass);
            for mat_id in 0..self.num_materials {
                let m = self.material_mass(mat_id, node_idx);
                if m < 1e-15 {
                    continue;
                }
                let v_mat = self.material_velocity(mat_id, node_idx);
                let v_rel = sub3(v_mat, v_cm);
                let v_coupled = add3(v_cm, scale3(v_rel, (1.0 - friction_coeff).max(0.0)));
                let idx = mat_id * n + node_idx;
                self.velocities[idx] = v_coupled;
            }
        }
    }
}
