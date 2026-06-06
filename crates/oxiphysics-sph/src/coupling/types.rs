//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// Interpolation weight for SPH-FEM data transfer.
#[derive(Debug, Clone)]
pub struct SphFemWeight {
    /// SPH particle index.
    pub sph_idx: usize,
    /// FEM node index.
    pub fem_idx: usize,
    /// Interpolation weight.
    pub weight: f64,
}
/// Region type for domain decomposition in atomistic-continuum coupling.
#[derive(Debug, Clone, PartialEq)]
pub enum RegionType {
    /// Pure atomistic/particle region.
    Atomistic,
    /// Pure continuum region.
    Continuum,
    /// Handshake/overlap region bridging atomistic and continuum.
    Handshake,
}
/// FEM-SPH interface coupling: transfers fields across a shared interface
/// between an FEM structural domain and an SPH fluid domain.
///
/// Provides methods for computing interface normals and transferring
/// stress tensors between SPH particles and FEM nodes at their boundary.
#[derive(Debug, Clone)]
pub struct FemSphCoupling {
    /// Positions of the FEM interface nodes.
    pub fem_node_positions: Vec<[f64; 3]>,
    /// Outward normals at each FEM interface node (pointing into SPH domain).
    pub fem_node_normals: Vec<[f64; 3]>,
    /// Smoothing length for SPH kernel evaluations.
    pub smoothing_length: f64,
    /// Stress tensor at each FEM node (Voigt notation: \[σxx,σyy,σzz,τxy,τyz,τxz\]).
    pub nodal_stress: Vec<[f64; 6]>,
    /// Force accumulated at each FEM node from stress transfer.
    pub nodal_forces: Vec<[f64; 3]>,
}
impl FemSphCoupling {
    /// Create a new FEM-SPH interface coupling.
    ///
    /// - `fem_node_positions` : positions of FEM interface nodes.
    /// - `smoothing_length`   : SPH kernel smoothing length h.
    pub fn new(fem_node_positions: Vec<[f64; 3]>, smoothing_length: f64) -> Self {
        let n = fem_node_positions.len();
        Self {
            fem_node_positions,
            fem_node_normals: vec![[0.0; 3]; n],
            smoothing_length,
            nodal_stress: vec![[0.0; 6]; n],
            nodal_forces: vec![[0.0; 3]; n],
        }
    }
    /// Compute outward interface normals at each FEM node using SPH particle
    /// positions as the neighbouring fluid domain.
    ///
    /// The interface normal at node k is estimated as the weighted average of
    /// the unit vectors pointing from each nearby SPH particle toward the node:
    ///
    /// ```text
    /// n_k = -( Σ_j W(r_kj, h) * r_hat_kj )  ,  then normalised.
    /// ```
    ///
    /// The sign convention is such that `n_k` points **from** the SPH fluid
    /// domain **into** the FEM solid domain.
    pub fn compute_interface_normal(&mut self, sph_positions: &[[f64; 3]]) {
        let h = self.smoothing_length;
        let n_nodes = self.fem_node_positions.len();
        for k in 0..n_nodes {
            let nk = self.fem_node_positions[k];
            let mut acc = [0.0_f64; 3];
            let mut w_sum = 0.0_f64;
            for &sp in sph_positions.iter() {
                let dx = [nk[0] - sp[0], nk[1] - sp[1], nk[2] - sp[2]];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let w = cubic_kernel(r, h);
                let inv_r = 1.0 / r;
                acc[0] -= w * dx[0] * inv_r;
                acc[1] -= w * dx[1] * inv_r;
                acc[2] -= w * dx[2] * inv_r;
                w_sum += w;
            }
            if w_sum > 1e-30 {
                let mag = (acc[0] * acc[0] + acc[1] * acc[1] + acc[2] * acc[2]).sqrt();
                if mag > 1e-14 {
                    let inv_mag = 1.0 / mag;
                    self.fem_node_normals[k] =
                        [acc[0] * inv_mag, acc[1] * inv_mag, acc[2] * inv_mag];
                } else {
                    self.fem_node_normals[k] = [0.0; 3];
                }
            } else {
                self.fem_node_normals[k] = [0.0; 3];
            }
        }
    }
    /// Transfer SPH stress tensors to FEM interface nodal forces.
    ///
    /// For each FEM node k the traction `t_k = σ_sph · n_k` is accumulated
    /// from nearby SPH particles using a kernel-weighted average:
    ///
    /// ```text
    /// f_k = Σ_j  W(r_kj, h) * (σ_j · n_k)
    /// ```
    ///
    /// The SPH stress tensor is supplied in Voigt notation
    /// `[σxx, σyy, σzz, τxy, τyz, τxz]`.
    ///
    /// The result is stored in `self.nodal_forces`.
    pub fn transfer_stress(&mut self, sph_positions: &[[f64; 3]], sph_stress: &[[f64; 6]]) {
        let h = self.smoothing_length;
        let n_nodes = self.fem_node_positions.len();
        for f in &mut self.nodal_forces {
            *f = [0.0; 3];
        }
        for k in 0..n_nodes {
            let nk = self.fem_node_positions[k];
            let nhat = self.fem_node_normals[k];
            let mut force = [0.0_f64; 3];
            let mut w_sum = 0.0_f64;
            for (j, &sp) in sph_positions.iter().enumerate() {
                let dx = [nk[0] - sp[0], nk[1] - sp[1], nk[2] - sp[2]];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let w = cubic_kernel(r, h);
                let s = sph_stress[j];
                let tx = s[0] * nhat[0] + s[3] * nhat[1] + s[5] * nhat[2];
                let ty = s[3] * nhat[0] + s[1] * nhat[1] + s[4] * nhat[2];
                let tz = s[5] * nhat[0] + s[4] * nhat[1] + s[2] * nhat[2];
                force[0] += w * tx;
                force[1] += w * ty;
                force[2] += w * tz;
                w_sum += w;
            }
            if w_sum > 1e-30 {
                let inv_w = 1.0 / w_sum;
                self.nodal_forces[k] = [force[0] * inv_w, force[1] * inv_w, force[2] * inv_w];
            } else {
                self.nodal_forces[k] = [0.0; 3];
            }
        }
    }
}
/// A lightweight fluid-structure interaction (FSI) helper that aggregates
/// buoyancy, drag, added-mass, and pressure-gradient forces on a rigid body
/// from SPH particle data.
#[derive(Debug, Clone)]
pub struct FluidStructureInteraction {
    /// Fluid density (kg m⁻³).
    pub fluid_density: f64,
    /// Gravitational acceleration (m s⁻²), positive downward convention.
    pub gravity: f64,
    /// Drag coefficient C_d.
    pub cd: f64,
    /// Added-mass coefficient C_m.
    pub cm: f64,
}
impl FluidStructureInteraction {
    /// Create a new FSI object.
    pub fn new(fluid_density: f64, gravity: f64, cd: f64, cm: f64) -> Self {
        Self {
            fluid_density,
            gravity,
            cd,
            cm,
        }
    }
    /// Total hydrodynamic force on a rigid body at position `body_pos`
    /// moving at velocity `body_vel` with volume `body_vol`.
    ///
    /// Includes:
    /// 1. Archimedes buoyancy: F_b = ρ_f g V_body ẑ
    /// 2. SPH drag from nearby fluid particles (two-way coupling)
    /// 3. Added-mass correction when `a_body` is provided.
    ///
    /// # Arguments
    /// - `body_pos`       : body centroid.
    /// - `body_vel`       : body velocity.
    /// - `body_vol`       : body volume (m³).
    /// - `a_body`         : body acceleration (used for added-mass term).
    /// - `sph_positions`  : fluid particle positions.
    /// - `sph_velocities` : fluid particle velocities.
    /// - `sph_volumes`    : fluid particle volumes.
    /// - `support_radius` : neighbourhood radius for drag computation.
    ///
    /// # Returns
    /// Total force vector (N).
    pub fn total_force(
        &self,
        body_pos: [f64; 3],
        body_vel: [f64; 3],
        body_vol: f64,
        a_body: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_velocities: &[[f64; 3]],
        sph_volumes: &[f64],
        support_radius: f64,
    ) -> [f64; 3] {
        let fb_z = self.fluid_density * self.gravity * body_vol;
        let f_buoy = [0.0, 0.0, fb_z];
        let coupling = TwoWayCouplingForce::new(self.fluid_density, support_radius, self.cd);
        let f_drag = coupling.body_force_from_fluid(
            body_pos,
            body_vel,
            sph_positions,
            sph_velocities,
            sph_volumes,
            support_radius,
        );
        let m_added = self.cm * self.fluid_density * body_vol;
        let f_added = [
            -m_added * a_body[0],
            -m_added * a_body[1],
            -m_added * a_body[2],
        ];
        [
            f_buoy[0] + f_drag[0] + f_added[0],
            f_buoy[1] + f_drag[1] + f_added[1],
            f_buoy[2] + f_drag[2] + f_added[2],
        ]
    }
    /// Compute pressure force on a rigid body from SPH pressure field.
    ///
    /// F_p = -∫ p n dA  ≈ -Σ_k m_k/ρ_k * p_k * ∇W(x_body - x_k, h)
    ///
    /// Each SPH particle contributes to the body's pressure force via the
    /// kernel gradient evaluated at the body surface.
    pub fn pressure_force(
        &self,
        body_pos: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_pressures: &[f64],
        sph_masses: &[f64],
        sph_densities: &[f64],
        h: f64,
    ) -> [f64; 3] {
        let mut force = [0.0_f64; 3];
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r_vec = sub3(body_pos, pos);
            let r = norm3(r_vec);
            if r < 1e-30 || r > 2.0 * h {
                continue;
            }
            let r_hat = scale3(1.0 / r, r_vec);
            let dw_dr = cubic_kernel_grad(r, h);
            let grad_w = scale3(dw_dr, r_hat);
            let weight = sph_masses[k] / sph_densities[k].max(1e-30);
            let contrib = scale3(-weight * sph_pressures[k], grad_w);
            force = add3(force, contrib);
        }
        force
    }
}
/// DEM-SPH coupling engine.
///
/// Computes hydrodynamic forces (buoyancy + drag + pressure) acting on each
/// DEM particle from the surrounding SPH fluid, and the equal/opposite
/// reaction forces on the SPH particles.
///
/// Algorithm follows:
/// - Lagrangian-based DEM-SPH coupling (Tsuji et al. 1992 adapted for SPH).
/// - Each DEM particle is treated as a moving boundary.
/// - Reaction forces on fluid use the conservative momentum exchange form.
#[derive(Debug, Clone)]
pub struct DemSphCoupling {
    /// Fluid density (kg m⁻³).
    pub fluid_density: f64,
    /// Gravitational acceleration vector (m s⁻²).
    pub gravity: [f64; 3],
    /// Drag coefficient.
    pub cd: f64,
    /// Added-mass coefficient.
    pub cm: f64,
    /// SPH smoothing length (m).
    pub h: f64,
}
impl DemSphCoupling {
    /// Create a new DEM-SPH coupling object.
    pub fn new(fluid_density: f64, gravity: [f64; 3], cd: f64, cm: f64, h: f64) -> Self {
        Self {
            fluid_density,
            gravity,
            cd,
            cm,
            h,
        }
    }
    /// Compute hydrodynamic force on a single DEM particle.
    ///
    /// Combines:
    /// 1. Buoyancy: F_b = -ρ_f g V
    /// 2. Drag:     F_d = -cd ρ_f V (v_dem - v_fluid_avg)
    /// 3. Pressure: F_p from SPH pressure field
    ///
    /// where `v_fluid_avg` is the kernel-weighted average fluid velocity
    /// at the DEM particle position.
    ///
    /// # Arguments
    /// - `dem`            : the DEM particle.
    /// - `dem_accel`      : acceleration of the DEM particle (for added mass).
    /// - `sph_positions`  : SPH particle positions.
    /// - `sph_velocities` : SPH particle velocities.
    /// - `sph_pressures`  : SPH particle pressures.
    /// - `sph_masses`     : SPH particle masses.
    /// - `sph_densities`  : SPH particle densities.
    ///
    /// # Returns
    /// Hydrodynamic force on the DEM particle (N).
    pub fn hydro_force_on_dem(
        &self,
        dem: &DemSphParticle,
        dem_accel: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_velocities: &[[f64; 3]],
        sph_pressures: &[f64],
        sph_masses: &[f64],
        sph_densities: &[f64],
    ) -> [f64; 3] {
        let support = 2.0 * self.h;
        let vol = dem.volume();
        let f_buoy = [
            -self.fluid_density * self.gravity[0] * vol,
            -self.fluid_density * self.gravity[1] * vol,
            -self.fluid_density * self.gravity[2] * vol,
        ];
        let mut v_fluid = [0.0_f64; 3];
        let mut weight_sum = 0.0_f64;
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r_vec = sub3(pos, dem.position);
            let r = norm3(r_vec);
            if r > support {
                continue;
            }
            let w = cubic_kernel(r, self.h);
            let vol_k = sph_masses[k] / sph_densities[k].max(1e-30);
            let wv = w * vol_k;
            v_fluid[0] += wv * sph_velocities[k][0];
            v_fluid[1] += wv * sph_velocities[k][1];
            v_fluid[2] += wv * sph_velocities[k][2];
            weight_sum += wv;
        }
        if weight_sum > 1e-30 {
            v_fluid[0] /= weight_sum;
            v_fluid[1] /= weight_sum;
            v_fluid[2] /= weight_sum;
        }
        let dv = sub3(dem.velocity, v_fluid);
        let f_drag = scale3(-self.cd * self.fluid_density * vol, dv);
        let mut f_press = [0.0_f64; 3];
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r_vec = sub3(dem.position, pos);
            let r = norm3(r_vec);
            if r < 1e-30 || r > support {
                continue;
            }
            let r_hat = scale3(1.0 / r, r_vec);
            let dw_dr = cubic_kernel_grad(r, self.h);
            let grad_w = scale3(dw_dr, r_hat);
            let vol_k = sph_masses[k] / sph_densities[k].max(1e-30);
            f_press = add3(f_press, scale3(-vol * vol_k * sph_pressures[k], grad_w));
        }
        let m_added = self.cm * self.fluid_density * vol;
        let f_added = [
            -m_added * dem_accel[0],
            -m_added * dem_accel[1],
            -m_added * dem_accel[2],
        ];
        [
            f_buoy[0] + f_drag[0] + f_press[0] + f_added[0],
            f_buoy[1] + f_drag[1] + f_press[1] + f_added[1],
            f_buoy[2] + f_drag[2] + f_press[2] + f_added[2],
        ]
    }
    /// Compute reaction forces on SPH particles from a DEM particle.
    ///
    /// Each SPH particle within the support radius receives an equal and
    /// opposite force contribution weighted by the SPH kernel.
    ///
    /// # Returns
    /// Per-particle force corrections (same length as `sph_positions`).
    pub fn reaction_forces_on_sph(
        &self,
        dem: &DemSphParticle,
        dem_force: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_masses: &[f64],
        sph_densities: &[f64],
    ) -> Vec<[f64; 3]> {
        let support = 2.0 * self.h;
        let n = sph_positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        let mut weight_sum = 0.0_f64;
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r = norm3(sub3(pos, dem.position));
            if r > support {
                continue;
            }
            let w = cubic_kernel(r, self.h);
            let vol_k = sph_masses[k] / sph_densities[k].max(1e-30);
            weight_sum += w * vol_k;
        }
        if weight_sum < 1e-30 {
            return forces;
        }
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r = norm3(sub3(pos, dem.position));
            if r > support {
                continue;
            }
            let w = cubic_kernel(r, self.h);
            let vol_k = sph_masses[k] / sph_densities[k].max(1e-30);
            let frac = w * vol_k / weight_sum;
            forces[k] = scale3(-frac, dem_force);
        }
        forces
    }
}
/// Sub-stepping strategy for coupled systems.
#[derive(Debug, Clone, PartialEq)]
pub enum SubstepStrategy {
    /// Both systems share the same time step.
    Monolithic,
    /// SPH uses a finer sub-step.
    SphSubStep {
        /// Number of SPH sub-steps per rigid-body step.
        sph_steps_per_rigid: usize,
    },
    /// Rigid body uses a finer sub-step.
    RigidSubStep {
        /// Number of rigid-body sub-steps per SPH step.
        rigid_steps_per_sph: usize,
    },
}
/// An axis-aligned bounding box representing a coupling region.
#[derive(Debug, Clone)]
pub struct CouplingRegion {
    /// Minimum corner of the region.
    pub min: [f64; 3],
    /// Maximum corner of the region.
    pub max: [f64; 3],
    /// Type of this region.
    pub region_type: RegionType,
}
impl CouplingRegion {
    /// Returns `true` if the given position lies inside this region (inclusive).
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        pos[0] >= self.min[0]
            && pos[0] <= self.max[0]
            && pos[1] >= self.min[1]
            && pos[1] <= self.max[1]
            && pos[2] >= self.min[2]
            && pos[2] <= self.max[2]
    }
    /// Returns the volume of this region.
    pub fn volume(&self) -> f64 {
        let dx = (self.max[0] - self.min[0]).max(0.0);
        let dy = (self.max[1] - self.min[1]).max(0.0);
        let dz = (self.max[2] - self.min[2]).max(0.0);
        dx * dy * dz
    }
    /// Returns the center of this region.
    pub fn center(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }
}
/// Represents a rigid body for coupling with SPH particles.
#[derive(Debug, Clone)]
pub struct RigidBody {
    /// Center of mass position.
    pub position: [f64; 3],
    /// Velocity of center of mass.
    pub velocity: [f64; 3],
    /// Angular velocity.
    pub angular_velocity: [f64; 3],
    /// Mass (kg).
    pub mass: f64,
    /// Moment of inertia (simplified scalar for sphere).
    pub inertia: f64,
    /// Representative radius (m).
    pub radius: f64,
}
impl RigidBody {
    /// Create a new spherical rigid body.
    pub fn sphere(position: [f64; 3], radius: f64, density: f64) -> Self {
        let mass = (4.0 / 3.0) * std::f64::consts::PI * radius.powi(3) * density;
        let inertia = 0.4 * mass * radius * radius;
        Self {
            position,
            velocity: [0.0; 3],
            angular_velocity: [0.0; 3],
            mass,
            inertia,
            radius,
        }
    }
    /// Velocity at a surface point (includes angular velocity contribution).
    pub fn surface_velocity(&self, surface_point: [f64; 3]) -> [f64; 3] {
        let r = [
            surface_point[0] - self.position[0],
            surface_point[1] - self.position[1],
            surface_point[2] - self.position[2],
        ];
        let omega_cross_r = cross3(self.angular_velocity, r);
        [
            self.velocity[0] + omega_cross_r[0],
            self.velocity[1] + omega_cross_r[1],
            self.velocity[2] + omega_cross_r[2],
        ]
    }
    /// Signed distance from a point to the sphere surface (negative = inside).
    pub fn signed_distance(&self, pos: [f64; 3]) -> f64 {
        let dx = pos[0] - self.position[0];
        let dy = pos[1] - self.position[1];
        let dz = pos[2] - self.position[2];
        (dx * dx + dy * dy + dz * dz).sqrt() - self.radius
    }
    /// Apply force and torque to update the rigid body.
    pub fn apply_force(&mut self, force: [f64; 3], torque: [f64; 3], dt: f64) {
        let inv_m = 1.0 / self.mass;
        self.velocity[0] += force[0] * inv_m * dt;
        self.velocity[1] += force[1] * inv_m * dt;
        self.velocity[2] += force[2] * inv_m * dt;
        if self.inertia > 1e-30 {
            let inv_i = 1.0 / self.inertia;
            self.angular_velocity[0] += torque[0] * inv_i * dt;
            self.angular_velocity[1] += torque[1] * inv_i * dt;
            self.angular_velocity[2] += torque[2] * inv_i * dt;
        }
        self.position[0] += self.velocity[0] * dt;
        self.position[1] += self.velocity[1] * dt;
        self.position[2] += self.velocity[2] * dt;
    }
}
/// FEM element type used in weak coupling.
#[derive(Debug, Clone, PartialEq)]
pub enum FemElementType {
    /// Linear tetrahedral element (4 nodes).
    Tet4,
    /// Quadratic tetrahedral element (10 nodes).
    Tet10,
    /// Linear hexahedral element (8 nodes).
    Hex8,
}
/// Computes **two-way coupling** forces between a rigid body and the
/// surrounding SPH fluid.
///
/// The fluid exerts a drag + pressure force on the rigid body, and by
/// Newton's third law the rigid body exerts an equal and opposite reaction
/// force on each fluid particle.
///
/// Reference: Akinci et al. (2012), "Versatile Rigid-Fluid Coupling for
/// Incompressible SPH".
#[derive(Debug, Clone)]
pub struct TwoWayCouplingForce {
    /// Fluid density (kg m⁻³).
    pub fluid_density: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Drag coefficient (dimensionless).
    pub cd: f64,
}
impl TwoWayCouplingForce {
    /// Create a new two-way coupling force calculator.
    pub fn new(fluid_density: f64, h: f64, cd: f64) -> Self {
        Self {
            fluid_density,
            h,
            cd,
        }
    }
    /// Compute force on the rigid body from nearby SPH particles.
    ///
    /// Each boundary particle `k` at position `x_k` and velocity `v_k`
    /// contributes a drag force proportional to the velocity difference
    /// relative to the rigid body velocity `v_body`:
    ///
    /// ```text
    /// dF = -cd * rho_f * V_k * (v_k - v_body)
    /// ```
    ///
    /// # Arguments
    /// - `body_pos`       : rigid body reference position.
    /// - `body_vel`       : rigid body velocity.
    /// - `sph_positions`  : SPH particle positions.
    /// - `sph_velocities` : SPH particle velocities.
    /// - `sph_volumes`    : SPH particle volumes.
    /// - `support_radius` : interaction radius (particles farther away are ignored).
    ///
    /// # Returns
    /// Force vector on the rigid body (N).
    pub fn body_force_from_fluid(
        &self,
        body_pos: [f64; 3],
        body_vel: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_velocities: &[[f64; 3]],
        sph_volumes: &[f64],
        support_radius: f64,
    ) -> [f64; 3] {
        let mut force = [0.0_f64; 3];
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r = sub3(pos, body_pos);
            if norm3(r) > support_radius {
                continue;
            }
            let dv = sub3(sph_velocities[k], body_vel);
            let contrib = scale3(-self.cd * self.fluid_density * sph_volumes[k], dv);
            force = add3(force, contrib);
        }
        force
    }
    /// Compute reaction forces on each SPH fluid particle (Newton's 3rd law).
    ///
    /// For each SPH particle close to the body the reaction force is
    /// the negative of its individual contribution to the body force.
    ///
    /// # Returns
    /// Vector of per-particle force corrections (same length as `sph_positions`).
    pub fn fluid_reaction_forces(
        &self,
        body_pos: [f64; 3],
        body_vel: [f64; 3],
        sph_positions: &[[f64; 3]],
        sph_velocities: &[[f64; 3]],
        sph_volumes: &[f64],
        support_radius: f64,
    ) -> Vec<[f64; 3]> {
        let n = sph_positions.len();
        let mut forces = vec![[0.0_f64; 3]; n];
        for (k, &pos) in sph_positions.iter().enumerate() {
            let r = sub3(pos, body_pos);
            if norm3(r) > support_radius {
                continue;
            }
            let dv = sub3(sph_velocities[k], body_vel);
            forces[k] = scale3(self.cd * self.fluid_density * sph_volumes[k], dv);
        }
        forces
    }
}
/// Simple interface tracking using a color function.
///
/// Each particle carries a color (0 or 1) indicating which phase it belongs to.
/// The interface is located where the color function transitions.
pub struct InterfaceTracker {
    /// Color values for each particle (0.0 = phase A, 1.0 = phase B).
    pub colors: Vec<f64>,
    /// Threshold for interface detection.
    pub threshold: f64,
}
impl InterfaceTracker {
    /// Create a new interface tracker.
    pub fn new(n_particles: usize) -> Self {
        Self {
            colors: vec![0.0; n_particles],
            threshold: 0.5,
        }
    }
    /// Identify particles near the interface (color between threshold bounds).
    ///
    /// A particle is "near the interface" if its smoothed color value is
    /// between `threshold - margin` and `threshold + margin`.
    pub fn interface_particles(&self, margin: f64) -> Vec<usize> {
        let lo = self.threshold - margin;
        let hi = self.threshold + margin;
        self.colors
            .iter()
            .enumerate()
            .filter(|&(_, &c)| c > lo && c < hi)
            .map(|(i, _)| i)
            .collect()
    }
    /// Set color for a particle.
    pub fn set_color(&mut self, idx: usize, color: f64) {
        if idx < self.colors.len() {
            self.colors[idx] = color;
        }
    }
    /// Smooth the color field using SPH kernel interpolation.
    pub fn smooth_colors(&mut self, positions: &[[f64; 3]], h: f64) {
        let n = positions.len();
        let old_colors = self.colors.clone();
        for i in 0..n {
            let mut weighted_c = 0.0_f64;
            let mut total_w = 0.0_f64;
            for j in 0..n {
                let r = dist(positions[i], positions[j]);
                let w = cubic_kernel(r, h);
                weighted_c += w * old_colors[j];
                total_w += w;
            }
            if total_w > 0.0 {
                self.colors[i] = weighted_c / total_w;
            }
        }
    }
    /// Compute the interface normal at a particle position.
    ///
    /// n = -∇C / |∇C| where C is the color field.
    pub fn interface_normal(&self, idx: usize, positions: &[[f64; 3]], h: f64) -> [f64; 3] {
        let pos_i = positions[idx];
        let mut grad_c = [0.0_f64; 3];
        for (j, pos_j) in positions.iter().enumerate() {
            let rij = [
                pos_i[0] - pos_j[0],
                pos_i[1] - pos_j[1],
                pos_i[2] - pos_j[2],
            ];
            let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
            if r < 1e-14 || r > 2.0 * h {
                continue;
            }
            let dw_dr = cubic_kernel_grad(r, h);
            let fac = dw_dr / r * (self.colors[j] - self.colors[idx]);
            for (gc, &rij_d) in grad_c.iter_mut().zip(rij.iter()) {
                *gc += fac * rij_d;
            }
        }
        let mag = (grad_c[0] * grad_c[0] + grad_c[1] * grad_c[1] + grad_c[2] * grad_c[2]).sqrt();
        if mag < 1e-14 {
            return [0.0; 3];
        }
        [-grad_c[0] / mag, -grad_c[1] / mag, -grad_c[2] / mag]
    }
}
/// A boundary (ghost/dummy) particle used to enforce no-slip conditions
/// at rigid walls in the SPH framework.
///
/// The boundary particle mirrors the velocity of the wall and applies a
/// repulsion force to prevent fluid particles from penetrating.
#[derive(Debug, Clone)]
pub struct SphBoundaryParticle {
    /// Position of the dummy particle.
    pub position: [f64; 3],
    /// Wall velocity (no-slip: equals wall velocity at that point).
    pub velocity: [f64; 3],
    /// Pressure assigned to the dummy particle.
    pub pressure: f64,
    /// Density assigned to the dummy particle.
    pub density: f64,
    /// Volume of the dummy particle.
    pub volume: f64,
}
impl SphBoundaryParticle {
    /// Create a new boundary particle.
    pub fn new(
        position: [f64; 3],
        wall_velocity: [f64; 3],
        pressure: f64,
        density: f64,
        volume: f64,
    ) -> Self {
        Self {
            position,
            velocity: wall_velocity,
            pressure,
            density,
            volume,
        }
    }
    /// Generalized wall boundary condition (Adami et al. 2012):
    /// set the dummy particle pressure to enforce zero normal acceleration
    /// at the wall.
    ///
    /// p_wall = (Σ_f p_f W_wf + Σ_f ρ_f (g - a_wall) · r_wf W_wf) / Σ_f W_wf
    ///
    /// `fluid_positions`, `fluid_pressures`, `fluid_densities` are for nearby
    /// fluid particles; `a_wall` is the wall acceleration; `h` is the kernel
    /// smoothing length.
    pub fn update_pressure_adami(
        &mut self,
        fluid_positions: &[[f64; 3]],
        fluid_pressures: &[f64],
        fluid_densities: &[f64],
        a_wall: [f64; 3],
        gravity: [f64; 3],
        h: f64,
    ) {
        let mut num = 0.0_f64;
        let mut denom = 0.0_f64;
        for (k, &pos_f) in fluid_positions.iter().enumerate() {
            let r_vec = sub3(self.position, pos_f);
            let r = norm3(r_vec);
            if r > 2.0 * h {
                continue;
            }
            let w = cubic_kernel(r, h);
            let g_minus_a = sub3(gravity, a_wall);
            let body_term = fluid_densities[k] * dot3(g_minus_a, r_vec);
            num += (fluid_pressures[k] + body_term) * w;
            denom += w;
        }
        self.pressure = if denom > 1e-30 { num / denom } else { 0.0 };
    }
}
/// SPH-to-FEM coupling: maps fluid pressure fields onto FEM nodal forces.
#[derive(Debug, Clone)]
pub struct FemCoupling {
    /// FEM node indices this coupling applies to.
    pub node_ids: Vec<usize>,
    /// FEM node positions.
    pub node_positions: Vec<[f64; 3]>,
    /// Accumulated nodal forces (updated by `transfer_sph_pressure_to_nodes`).
    pub nodal_forces: Vec<[f64; 3]>,
}
impl FemCoupling {
    /// Create a new FEM coupling for the given nodes.
    pub fn new(node_ids: Vec<usize>, node_positions: Vec<[f64; 3]>) -> Self {
        let n = node_positions.len();
        Self {
            node_ids,
            node_positions,
            nodal_forces: vec![[0.0; 3]; n],
        }
    }
    /// Transfer SPH pressure field to FEM nodal forces.
    ///
    /// Each SPH particle contributes a pressure-weighted kernel value to
    /// each FEM node within the smoothing radius h:
    ///   f_node += W(r, h) · p_sph · r_hat
    pub fn transfer_sph_pressure_to_nodes(
        &mut self,
        sph_pos: &[[f64; 3]],
        sph_pressure: &[f64],
        h: f64,
    ) {
        for f in &mut self.nodal_forces {
            *f = [0.0; 3];
        }
        for (k, &np) in self.node_positions.iter().enumerate() {
            for (j, &sp) in sph_pos.iter().enumerate() {
                let rij = [np[0] - sp[0], np[1] - sp[1], np[2] - sp[2]];
                let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
                if r < 1e-14 || r > 2.0 * h {
                    continue;
                }
                let w = cubic_kernel(r, h);
                if w <= 0.0 {
                    continue;
                }
                let p = sph_pressure[j];
                let fac = w * p / r;
                self.nodal_forces[k][0] += fac * rij[0];
                self.nodal_forces[k][1] += fac * rij[1];
                self.nodal_forces[k][2] += fac * rij[2];
            }
        }
    }
}
/// Coupling descriptor for a rigid body embedded in an SPH fluid.
///
/// Stores the body's kinematic state and provides helpers to compute
/// SPH-based fluid forces on the body.
#[derive(Debug, Clone)]
pub struct RigidBodyCoupling {
    /// Unique identifier for this rigid body.
    pub body_id: usize,
    /// Body mass (kg).
    pub mass: f64,
    /// Center-of-mass position (m).
    pub position: [f64; 3],
    /// Center-of-mass velocity (m/s).
    pub velocity: [f64; 3],
}
impl RigidBodyCoupling {
    /// Create a new rigid body coupling descriptor.
    pub fn new(body_id: usize, mass: f64, position: [f64; 3], velocity: [f64; 3]) -> Self {
        Self {
            body_id,
            mass,
            position,
            velocity,
        }
    }
    /// Compute the SPH coupling force on the body from the surrounding fluid.
    ///
    /// Pressure-weighted interpolation:
    ///   F = Σ_j  (m_j / ρ_j) · p_j · ∇W(|r - r_j|, h)
    ///
    /// where the gradient points from the body centre toward particle j and
    /// is approximated using the cubic-spline kernel derivative.
    pub fn sph_coupling_force(
        &self,
        fluid_particles: &[[f64; 3]],
        fluid_pressures: &[f64],
        fluid_rho: &[f64],
        h: f64,
    ) -> [f64; 3] {
        let mut force = [0.0_f64; 3];
        let pos = self.position;
        for (k, &fp) in fluid_particles.iter().enumerate() {
            let rij = [fp[0] - pos[0], fp[1] - pos[1], fp[2] - pos[2]];
            let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
            if r < 1e-14 || r > 2.0 * h {
                continue;
            }
            let dw_dr = cubic_kernel_grad(r, h);
            let rho_j = fluid_rho[k].max(1e-14);
            let weight = fluid_pressures[k] / rho_j;
            let fac = weight * dw_dr / r;
            force[0] += fac * rij[0];
            force[1] += fac * rij[1];
            force[2] += fac * rij[2];
        }
        force
    }
}
/// Simplified Bridging Domain Method (BDM) coupling structure.
///
/// Holds the three domain regions and provides helpers for overlap detection
/// and coupling weight queries.
#[derive(Debug, Clone)]
pub struct BridgingDomainMethod {
    /// Pure atomistic region.
    pub atomistic_region: CouplingRegion,
    /// Pure continuum region.
    pub continuum_region: CouplingRegion,
    /// Overlap (handshake) region shared between the two descriptions.
    pub overlap_region: CouplingRegion,
}
impl BridgingDomainMethod {
    /// Creates a new `BridgingDomainMethod` from three pre-constructed regions.
    pub fn new(
        atomistic: CouplingRegion,
        continuum: CouplingRegion,
        overlap: CouplingRegion,
    ) -> Self {
        Self {
            atomistic_region: atomistic,
            continuum_region: continuum,
            overlap_region: overlap,
        }
    }
    /// Returns `true` if `pos` lies inside the overlap region.
    pub fn is_in_overlap(&self, pos: [f64; 3]) -> bool {
        self.overlap_region.contains(pos)
    }
    /// Computes the coupling weight for a position in the overlap region.
    ///
    /// The weight is based on the distance along the first axis from the
    /// overlap region's minimum to its maximum, using the quintic transition
    /// function (1 near the atomistic side, 0 near the continuum side).
    pub fn coupling_weight(&self, pos: [f64; 3]) -> f64 {
        let d = pos[0];
        let d_min = self.overlap_region.min[0];
        let d_max = self.overlap_region.max[0];
        weight_function_transition(d, d_min, d_max)
    }
}
/// Hertz contact force model for DEM-DEM and DEM-wall collisions.
///
/// Normal force: F_n = k_n * δ^(3/2)  (Hertz)
/// Damping:      F_d = -γ_n * ṙ_n   (linear dashpot)
/// Tangential:   F_t = min(μ F_n, k_t * δ_t)  (Coulomb friction limit)
#[derive(Debug, Clone)]
pub struct DemContactForce {
    /// Normal stiffness coefficient (N m⁻³/²).
    pub k_n: f64,
    /// Normal damping coefficient (N s m⁻¹).
    pub gamma_n: f64,
    /// Tangential stiffness (N m⁻¹).
    pub k_t: f64,
    /// Friction coefficient μ.
    pub mu: f64,
}
impl DemContactForce {
    /// Create a new Hertz contact force model.
    pub fn new(k_n: f64, gamma_n: f64, k_t: f64, mu: f64) -> Self {
        Self {
            k_n,
            gamma_n,
            k_t,
            mu,
        }
    }
    /// Compute the normal contact force between two DEM spheres.
    ///
    /// # Arguments
    /// - `pos_i`, `pos_j`   : positions of the two spheres.
    /// - `vel_i`, `vel_j`   : velocities of the two spheres.
    /// - `radius_i`, `radius_j` : radii of the two spheres.
    ///
    /// # Returns
    /// Force on particle i from particle j (N). Force on j is the negative.
    /// Returns `[0;3]` if particles are not in contact.
    pub fn normal_contact_force(
        &self,
        pos_i: [f64; 3],
        pos_j: [f64; 3],
        vel_i: [f64; 3],
        vel_j: [f64; 3],
        radius_i: f64,
        radius_j: f64,
    ) -> [f64; 3] {
        let r_vec = sub3(pos_i, pos_j);
        let r = norm3(r_vec);
        let r_contact = radius_i + radius_j;
        if r >= r_contact || r < 1e-30 {
            return [0.0; 3];
        }
        let n_hat = scale3(1.0 / r, r_vec);
        let delta = r_contact - r;
        let f_hertz = self.k_n * delta.powf(1.5);
        let dv = sub3(vel_i, vel_j);
        let v_n = dot3(dv, n_hat);
        let f_damp = -self.gamma_n * v_n;
        let f_n = (f_hertz + f_damp).max(0.0);
        scale3(f_n, n_hat)
    }
    /// Compute Coulomb-limited tangential friction force.
    ///
    /// Requires the accumulated tangential displacement `delta_t` (a 3-vector
    /// maintained by the caller across time steps).
    ///
    /// Returns the tangential force on particle i.
    pub fn tangential_force(
        &self,
        pos_i: [f64; 3],
        pos_j: [f64; 3],
        vel_i: [f64; 3],
        vel_j: [f64; 3],
        radius_i: f64,
        radius_j: f64,
        delta_t: [f64; 3],
    ) -> [f64; 3] {
        let r_vec = sub3(pos_i, pos_j);
        let r = norm3(r_vec);
        let r_contact = radius_i + radius_j;
        if r >= r_contact || r < 1e-30 {
            return [0.0; 3];
        }
        let n_hat = scale3(1.0 / r, r_vec);
        let dv = sub3(vel_i, vel_j);
        let v_n = dot3(dv, n_hat);
        let v_t = sub3(dv, scale3(v_n, n_hat));
        let _ = v_t;
        let f_t_spring = scale3(-self.k_t, delta_t);
        let f_t_mag = norm3(f_t_spring);
        let delta = r_contact - r;
        let f_n_mag = self.k_n * delta.powf(1.5);
        if f_t_mag > self.mu * f_n_mag {
            let dt_hat = if f_t_mag > 1e-30 {
                scale3(1.0 / f_t_mag, f_t_spring)
            } else {
                [0.0; 3]
            };
            scale3(self.mu * f_n_mag, dt_hat)
        } else {
            f_t_spring
        }
    }
}
/// Track coupled energy exchange between SPH fluid and rigid body.
///
/// Keeps a running balance of:
/// - Work done by fluid on body (from SPH boundary forces).
/// - Work done by body on fluid (reaction).
#[derive(Debug, Clone, Default)]
pub struct CoupledEnergyBalance {
    /// Total energy transferred from fluid to rigid body.
    pub fluid_to_body: f64,
    /// Total energy transferred from rigid body to fluid.
    pub body_to_fluid: f64,
    /// Accumulated kinetic energy of rigid body over time.
    pub rigid_ke_integral: f64,
}
impl CoupledEnergyBalance {
    /// Create a new energy balance tracker.
    pub fn new() -> Self {
        Self::default()
    }
    /// Update the energy balance for one time step.
    ///
    /// - `force_on_body`  : total SPH force on the rigid body.
    /// - `body_velocity`  : velocity of rigid body CoM.
    /// - `body_ke`        : kinetic energy of the rigid body.
    /// - `dt`             : time step size.
    pub fn update(
        &mut self,
        force_on_body: [f64; 3],
        body_velocity: [f64; 3],
        body_ke: f64,
        dt: f64,
    ) {
        let power: f64 = force_on_body
            .iter()
            .zip(body_velocity.iter())
            .map(|(f, v)| f * v)
            .sum();
        if power >= 0.0 {
            self.fluid_to_body += power * dt;
        } else {
            self.body_to_fluid += (-power) * dt;
        }
        self.rigid_ke_integral += body_ke * dt;
    }
    /// Net energy transferred from fluid to body (positive = body gained energy).
    pub fn net_transfer(&self) -> f64 {
        self.fluid_to_body - self.body_to_fluid
    }
}
/// A simple FEM element (tetrahedron) for coupling.
#[derive(Debug, Clone)]
pub struct FemElement {
    /// Node indices.
    pub nodes: [usize; 4],
    /// Element volume.
    pub volume: f64,
}
/// A simple FEM node for coupling purposes.
#[derive(Debug, Clone)]
pub struct FemNode {
    /// Position.
    pub position: [f64; 3],
    /// Velocity.
    pub velocity: [f64; 3],
    /// Force from SPH coupling.
    pub coupling_force: [f64; 3],
}
/// SPH rigid immersed boundary: a collection of boundary particles
/// with outward normals.
///
/// Used to compute reaction forces exerted by the fluid on the boundary.
#[derive(Debug, Clone)]
pub struct SphRigidImmersedBoundary {
    /// Positions of boundary (marker) particles.
    pub boundary_particles: Vec<[f64; 3]>,
    /// Outward unit normals at each boundary particle.
    pub normals: Vec<[f64; 3]>,
}
impl SphRigidImmersedBoundary {
    /// Create a new immersed boundary.
    pub fn new(boundary_particles: Vec<[f64; 3]>, normals: Vec<[f64; 3]>) -> Self {
        Self {
            boundary_particles,
            normals,
        }
    }
    /// Compute the reaction force on each boundary particle from the fluid.
    ///
    /// For each boundary marker, sum the pressure contributions from nearby
    /// fluid particles:
    ///   f_k = Σ_j  W(r_kj, h) · p_j · n_k
    pub fn compute_reaction_forces(
        &self,
        fluid_particles: &[[f64; 3]],
        pressures: &[f64],
        h: f64,
    ) -> Vec<[f64; 3]> {
        let nb = self.boundary_particles.len();
        let mut forces = vec![[0.0_f64; 3]; nb];
        for (force, (bp, nrm)) in forces
            .iter_mut()
            .zip(self.boundary_particles.iter().zip(self.normals.iter()))
        {
            let mut p_sum = 0.0_f64;
            for (&fp, &pj) in fluid_particles.iter().zip(pressures.iter()) {
                let rij = [bp[0] - fp[0], bp[1] - fp[1], bp[2] - fp[2]];
                let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
                if r > 2.0 * h {
                    continue;
                }
                let w = cubic_kernel(r, h);
                p_sum += w * pj;
            }
            *force = [p_sum * nrm[0], p_sum * nrm[1], p_sum * nrm[2]];
        }
        forces
    }
}
/// State of a single DEM (Discrete Element Method) particle immersed in an
/// SPH fluid domain.
#[derive(Debug, Clone)]
pub struct DemSphParticle {
    /// Particle position (m).
    pub position: [f64; 3],
    /// Particle velocity (m s⁻¹).
    pub velocity: [f64; 3],
    /// Particle angular velocity (rad s⁻¹).
    pub omega: [f64; 3],
    /// Particle radius (m).
    pub radius: f64,
    /// Particle mass (kg).
    pub mass: f64,
    /// Moment of inertia scalar (kg m²) — assumes a sphere.
    pub inertia: f64,
    /// Material density (kg m⁻³).
    pub density: f64,
}
impl DemSphParticle {
    /// Create a new DEM-SPH particle (spherical, mass and inertia computed
    /// automatically from radius and density).
    pub fn new_sphere(position: [f64; 3], velocity: [f64; 3], radius: f64, density: f64) -> Self {
        let vol = (4.0 / 3.0) * std::f64::consts::PI * radius * radius * radius;
        let mass = density * vol;
        let inertia = 0.4 * mass * radius * radius;
        Self {
            position,
            velocity,
            omega: [0.0; 3],
            radius,
            mass,
            inertia,
            density,
        }
    }
    /// Volume of the particle (m³).
    pub fn volume(&self) -> f64 {
        (4.0 / 3.0) * std::f64::consts::PI * self.radius * self.radius * self.radius
    }
    /// Integrate position and velocity using a symplectic Euler step.
    ///
    /// `f_total` is the total force (hydrodynamic + contact + gravity) acting
    /// on this particle.  `torque` is the total torque.
    pub fn integrate(&mut self, f_total: [f64; 3], torque: [f64; 3], dt: f64) {
        self.velocity[0] += dt * f_total[0] / self.mass;
        self.velocity[1] += dt * f_total[1] / self.mass;
        self.velocity[2] += dt * f_total[2] / self.mass;
        self.omega[0] += dt * torque[0] / self.inertia;
        self.omega[1] += dt * torque[1] / self.inertia;
        self.omega[2] += dt * torque[2] / self.inertia;
        self.position[0] += dt * self.velocity[0];
        self.position[1] += dt * self.velocity[1];
        self.position[2] += dt * self.velocity[2];
    }
}
/// SPH-FEM coupling interface.
///
/// Maps forces between SPH particles and FEM nodes using kernel interpolation.
pub struct SphFemCoupling {
    /// FEM nodes.
    pub fem_nodes: Vec<FemNode>,
    /// FEM elements.
    pub fem_elements: Vec<FemElement>,
    /// Smoothing length for interpolation.
    pub h: f64,
}
impl SphFemCoupling {
    /// Create a new coupling interface.
    pub fn new(h: f64) -> Self {
        Self {
            fem_nodes: Vec::new(),
            fem_elements: Vec::new(),
            h,
        }
    }
    /// Add a FEM node. Returns the node index.
    pub fn add_node(&mut self, pos: [f64; 3]) -> usize {
        let idx = self.fem_nodes.len();
        self.fem_nodes.push(FemNode {
            position: pos,
            velocity: [0.0; 3],
            coupling_force: [0.0; 3],
        });
        idx
    }
    /// Transfer forces from SPH particles to FEM nodes.
    ///
    /// Each SPH particle distributes its force to nearby FEM nodes
    /// weighted by the kernel function.
    pub fn transfer_forces_to_fem(&mut self, sph_positions: &[[f64; 3]], sph_forces: &[[f64; 3]]) {
        for node in &mut self.fem_nodes {
            node.coupling_force = [0.0; 3];
        }
        let h = self.h;
        for (i, &pos) in sph_positions.iter().enumerate() {
            let force = sph_forces[i];
            for node in &mut self.fem_nodes {
                let r2 = dist_sq(pos, node.position);
                let r = r2.sqrt();
                let w = cubic_kernel(r, h);
                if w > 0.0 {
                    node.coupling_force[0] += w * force[0];
                    node.coupling_force[1] += w * force[1];
                    node.coupling_force[2] += w * force[2];
                }
            }
        }
    }
    /// Number of FEM nodes.
    pub fn node_count(&self) -> usize {
        self.fem_nodes.len()
    }
}
/// Full 6-DOF (6 degrees of freedom) rigid body state.
#[derive(Debug, Clone)]
pub struct RigidBody6Dof {
    /// Centre of mass position.
    pub position: [f64; 3],
    /// Linear velocity of the centre of mass.
    pub velocity: [f64; 3],
    /// Orientation quaternion (w, x, y, z).
    pub quaternion: [f64; 4],
    /// Angular velocity in world frame (rad/s).
    pub angular_velocity: [f64; 3],
    /// Total mass (kg).
    pub mass: f64,
    /// Moment of inertia in body frame (diagonal, kg·m²).
    pub inertia: [f64; 3],
}
impl RigidBody6Dof {
    /// Create a new rigid body at rest.
    pub fn new(mass: f64, inertia: [f64; 3], position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            quaternion: [1.0, 0.0, 0.0, 0.0],
            angular_velocity: [0.0; 3],
            mass,
            inertia,
        }
    }
    /// Integrate linear motion: v += F/m * dt,  x += v * dt.
    pub fn integrate_linear(&mut self, force: [f64; 3], dt: f64) {
        let inv_m = 1.0 / self.mass.max(1e-30);
        for ((v, p), &f) in self
            .velocity
            .iter_mut()
            .zip(self.position.iter_mut())
            .zip(force.iter())
        {
            *v += f * inv_m * dt;
            *p += *v * dt;
        }
    }
    /// Integrate angular motion: ω += I⁻¹ τ dt,  q += ½ Ω q dt.
    pub fn integrate_angular(&mut self, torque: [f64; 3], dt: f64) {
        let alpha = [
            torque[0] / self.inertia[0].max(1e-30),
            torque[1] / self.inertia[1].max(1e-30),
            torque[2] / self.inertia[2].max(1e-30),
        ];
        for (av, &al) in self.angular_velocity.iter_mut().zip(alpha.iter()) {
            *av += al * dt;
        }
        let [w, qx, qy, qz] = self.quaternion;
        let [wx, wy, wz] = self.angular_velocity;
        let dw = 0.5 * (-wx * qx - wy * qy - wz * qz);
        let dx = 0.5 * (wx * w + wy * qz - wz * qy);
        let dy = 0.5 * (-wx * qz + wy * w + wz * qx);
        let dz = 0.5 * (wx * qy - wy * qx + wz * w);
        self.quaternion = [w + dw * dt, qx + dx * dt, qy + dy * dt, qz + dz * dt];
        self.normalize_quaternion();
    }
    /// Normalise the orientation quaternion to unit length.
    fn normalize_quaternion(&mut self) {
        let [w, x, y, z] = self.quaternion;
        let len = (w * w + x * x + y * y + z * z).sqrt().max(1e-30);
        self.quaternion = [w / len, x / len, y / len, z / len];
    }
    /// Kinetic energy: ½ m v² + ½ ω·I·ω.
    pub fn kinetic_energy(&self) -> f64 {
        let ke_lin = 0.5
            * self.mass
            * (self.velocity[0].powi(2) + self.velocity[1].powi(2) + self.velocity[2].powi(2));
        let ke_rot = 0.5
            * (self.inertia[0] * self.angular_velocity[0].powi(2)
                + self.inertia[1] * self.angular_velocity[1].powi(2)
                + self.inertia[2] * self.angular_velocity[2].powi(2));
        ke_lin + ke_rot
    }
    /// Velocity of a body point given its offset from CoM in world frame.
    pub fn point_velocity(&self, offset: [f64; 3]) -> [f64; 3] {
        let [wx, wy, wz] = self.angular_velocity;
        let [rx, ry, rz] = offset;
        [
            self.velocity[0] + wy * rz - wz * ry,
            self.velocity[1] + wz * rx - wx * rz,
            self.velocity[2] + wx * ry - wy * rx,
        ]
    }
}
/// Shared data buffer between atomistic and continuum regions.
///
/// Stores particle data that is exchanged at the coupling interface.
#[derive(Debug, Clone)]
pub struct CouplingBuffer {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle densities.
    pub densities: Vec<f64>,
    /// Coupling weights (e.g., from the transition function).
    pub weights: Vec<f64>,
}
impl CouplingBuffer {
    /// Creates an empty `CouplingBuffer`.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            densities: Vec::new(),
            weights: Vec::new(),
        }
    }
    /// Appends a particle to the buffer.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], rho: f64, weight: f64) {
        self.positions.push(pos);
        self.velocities.push(vel);
        self.densities.push(rho);
        self.weights.push(weight);
    }
    /// Returns the number of particles stored in the buffer.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Returns `true` if the buffer contains no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Kernel-weighted interpolation of velocity at a given query position.
    ///
    /// Uses a Gaussian-like kernel with support radius `h`.
    /// Returns the zero vector if no particles are within support.
    pub fn interpolate_velocity_at(&self, pos: [f64; 3], h: f64) -> [f64; 3] {
        let mut weighted_vel = [0.0_f64; 3];
        let mut total_weight = 0.0_f64;
        for (i, &p) in self.positions.iter().enumerate() {
            let r2 = dist_sq(pos, p);
            let w = cubic_kernel(r2.sqrt(), h) * self.weights[i];
            weighted_vel[0] += w * self.velocities[i][0];
            weighted_vel[1] += w * self.velocities[i][1];
            weighted_vel[2] += w * self.velocities[i][2];
            total_weight += w;
        }
        if total_weight > 0.0 {
            [
                weighted_vel[0] / total_weight,
                weighted_vel[1] / total_weight,
                weighted_vel[2] / total_weight,
            ]
        } else {
            [0.0, 0.0, 0.0]
        }
    }
    /// Kernel-weighted interpolation of density at a given query position.
    ///
    /// Returns `0.0` if no particles are within support.
    pub fn interpolate_density_at(&self, pos: [f64; 3], h: f64) -> f64 {
        let mut weighted_rho = 0.0_f64;
        let mut total_weight = 0.0_f64;
        for (i, &p) in self.positions.iter().enumerate() {
            let r2 = dist_sq(pos, p);
            let w = cubic_kernel(r2.sqrt(), h) * self.weights[i];
            weighted_rho += w * self.densities[i];
            total_weight += w;
        }
        if total_weight > 0.0 {
            weighted_rho / total_weight
        } else {
            0.0
        }
    }
    /// Interpolate a scalar field at a query position.
    pub fn interpolate_scalar_at(&self, pos: [f64; 3], h: f64, values: &[f64]) -> f64 {
        let mut weighted_val = 0.0_f64;
        let mut total_weight = 0.0_f64;
        for (i, &p) in self.positions.iter().enumerate() {
            let r2 = dist_sq(pos, p);
            let w = cubic_kernel(r2.sqrt(), h) * self.weights[i];
            weighted_val += w * values[i];
            total_weight += w;
        }
        if total_weight > 0.0 {
            weighted_val / total_weight
        } else {
            0.0
        }
    }
    /// Clear all data from the buffer.
    pub fn clear(&mut self) {
        self.positions.clear();
        self.velocities.clear();
        self.densities.clear();
        self.weights.clear();
    }
}
/// Coupled time-stepper state for an SPH-rigid body system.
#[derive(Debug, Clone)]
pub struct CoupledTimeStepper {
    /// Physical time elapsed.
    pub time: f64,
    /// Global time-step size.
    pub dt: f64,
    /// Accumulated SPH-to-rigid force (to be applied over the current step).
    pub accumulated_force: [f64; 3],
    /// Accumulated SPH-to-rigid torque.
    pub accumulated_torque: [f64; 3],
    /// Sub-stepping strategy.
    pub strategy: SubstepStrategy,
    /// Number of coupling cycles completed.
    pub cycles: usize,
}
impl CoupledTimeStepper {
    /// Create a new coupled time-stepper.
    pub fn new(dt: f64, strategy: SubstepStrategy) -> Self {
        Self {
            time: 0.0,
            dt,
            accumulated_force: [0.0; 3],
            accumulated_torque: [0.0; 3],
            strategy,
            cycles: 0,
        }
    }
    /// Accumulate force and torque from SPH onto the rigid body.
    pub fn accumulate_force(&mut self, force: [f64; 3], torque: [f64; 3]) {
        for i in 0..3 {
            self.accumulated_force[i] += force[i];
            self.accumulated_torque[i] += torque[i];
        }
    }
    /// Reset accumulated force/torque after applying to rigid body.
    pub fn reset_accumulation(&mut self) {
        self.accumulated_force = [0.0; 3];
        self.accumulated_torque = [0.0; 3];
    }
    /// Advance the clock by one global time step.
    pub fn advance(&mut self) {
        self.time += self.dt;
        self.cycles += 1;
        self.reset_accumulation();
    }
    /// Returns the sub-step size for SPH when using SphSubStep strategy.
    pub fn sph_substep_dt(&self) -> f64 {
        match self.strategy {
            SubstepStrategy::SphSubStep {
                sph_steps_per_rigid,
            } => self.dt / sph_steps_per_rigid.max(1) as f64,
            _ => self.dt,
        }
    }
    /// Returns the sub-step size for the rigid body when using RigidSubStep.
    pub fn rigid_substep_dt(&self) -> f64 {
        match self.strategy {
            SubstepStrategy::RigidSubStep {
                rigid_steps_per_sph,
            } => self.dt / rigid_steps_per_sph.max(1) as f64,
            _ => self.dt,
        }
    }
}
/// Added-mass effect for a rigid body accelerating through fluid.
///
/// The added mass M_a = C_m * ρ_f * V_body accounts for the inertia of the
/// surrounding fluid that must be accelerated with the body.
///
/// Reference: Lamb (1932), "Hydrodynamics", §§72–91.
#[derive(Debug, Clone)]
pub struct AddedMass {
    /// Added-mass coefficient (0.5 for a sphere).
    pub cm: f64,
    /// Fluid density (kg m⁻³).
    pub fluid_density: f64,
    /// Volume of the rigid body (m³).
    pub body_volume: f64,
}
impl AddedMass {
    /// Create a new added-mass model.
    pub fn new(cm: f64, fluid_density: f64, body_volume: f64) -> Self {
        Self {
            cm,
            fluid_density,
            body_volume,
        }
    }
    /// Effective added mass (kg).
    pub fn added_mass(&self) -> f64 {
        self.cm * self.fluid_density * self.body_volume
    }
    /// Effective total inertial mass = body_mass + added_mass.
    pub fn effective_mass(&self, body_mass: f64) -> f64 {
        body_mass + self.added_mass()
    }
    /// Compute the added-mass force F_a = -M_a * a_body.
    ///
    /// This force opposes the body acceleration `a_body`, effectively
    /// increasing the body's apparent inertia.
    pub fn force(&self, a_body: [f64; 3]) -> [f64; 3] {
        let m_a = self.added_mass();
        [-m_a * a_body[0], -m_a * a_body[1], -m_a * a_body[2]]
    }
    /// Correct acceleration for added-mass effect.
    ///
    /// Given the net external force `f_ext` and body mass `m_body`,
    /// return the actual acceleration including added-mass resistance.
    ///
    /// ```text
    /// a = f_ext / (m_body + M_a)
    /// ```
    pub fn corrected_acceleration(&self, f_ext: [f64; 3], m_body: f64) -> [f64; 3] {
        let m_eff = self.effective_mass(m_body).max(1e-30);
        [f_ext[0] / m_eff, f_ext[1] / m_eff, f_ext[2] / m_eff]
    }
}
/// Impulse-based coupling: when a rigid body hits an SPH free surface or
/// internal region, an impulse is exchanged to conserve momentum.
#[derive(Debug, Clone)]
pub struct FluidImpulseTransfer {
    /// Restitution coefficient (0 = fully inelastic, 1 = elastic).
    pub restitution: f64,
    /// Fluid density (kg m⁻³).
    pub fluid_density: f64,
}
impl FluidImpulseTransfer {
    /// Create a new impulse transfer model.
    pub fn new(restitution: f64, fluid_density: f64) -> Self {
        Self {
            restitution,
            fluid_density,
        }
    }
    /// Compute impulse on the rigid body when it impacts the fluid surface.
    ///
    /// The impulse magnitude is computed from momentum conservation:
    ///
    /// ```text
    /// J = -(1 + e) * m_body * m_fluid / (m_body + m_fluid) * Δv_n
    /// ```
    ///
    /// where Δv_n is the relative normal velocity and m_fluid is the
    /// effective fluid mass within the interaction volume.
    ///
    /// # Arguments
    /// - `m_body`         : rigid body mass (kg).
    /// - `v_body`         : rigid body velocity (m s⁻¹).
    /// - `v_fluid`        : effective fluid velocity at contact (m s⁻¹).
    /// - `n_hat`          : outward normal at contact point.
    /// - `fluid_volume`   : volume of fluid involved in impact (m³).
    ///
    /// # Returns
    /// Impulse vector on the body (N·s).
    pub fn compute_impulse(
        &self,
        m_body: f64,
        v_body: [f64; 3],
        v_fluid: [f64; 3],
        n_hat: [f64; 3],
        fluid_volume: f64,
    ) -> [f64; 3] {
        let m_fluid = self.fluid_density * fluid_volume;
        let dv = sub3(v_body, v_fluid);
        let dv_n = dot3(dv, n_hat);
        if dv_n >= 0.0 {
            return [0.0; 3];
        }
        let j_mag =
            -(1.0 + self.restitution) * dv_n * m_body * m_fluid / (m_body + m_fluid).max(1e-30);
        scale3(j_mag, n_hat)
    }
}
/// Rigid-SPH coupling that computes buoyancy forces **and torques** on a
/// floating or submerged rigid body.
///
/// The buoyancy torque arises when the centre of buoyancy (CoB) does not
/// coincide with the centre of mass (CoM): τ = r_CoB × F_buoyancy.
#[derive(Debug, Clone)]
pub struct RigidSphCoupling {
    /// Centre of mass of the rigid body in world frame.
    pub centre_of_mass: [f64; 3],
    /// Density of the fluid (kg m⁻³).
    pub fluid_density: f64,
    /// Gravitational acceleration magnitude (m s⁻²), positive downward.
    pub gravity: f64,
}
impl RigidSphCoupling {
    /// Create a new rigid-SPH coupling.
    ///
    /// - `centre_of_mass` : CoM position of the rigid body.
    /// - `fluid_density`  : ambient fluid density ρ_f (kg m⁻³).
    /// - `gravity`        : gravitational acceleration g (m s⁻²).
    pub fn new(centre_of_mass: [f64; 3], fluid_density: f64, gravity: f64) -> Self {
        Self {
            centre_of_mass,
            fluid_density,
            gravity,
        }
    }
    /// Compute the buoyancy torque on the rigid body.
    ///
    /// Each SPH particle that is **inside** the rigid body (the "displaced
    /// fluid") contributes a buoyancy force
    ///
    /// ```text
    /// dF_b = ρ_f * g * V_p * ẑ
    /// ```
    ///
    /// applied at the particle position `x_p`.  The torque about the CoM is
    ///
    /// ```text
    /// τ += (x_p - x_CoM) × dF_b
    /// ```
    ///
    /// Particles are considered "displaced" when they lie within `radius` of
    /// the rigid-body centroid (a simplified spherical body model).
    ///
    /// # Arguments
    /// - `sph_positions`  : positions of all SPH particles.
    /// - `particle_volume`: volume V_p of each SPH particle.
    /// - `radius`         : effective radius of the rigid body for inclusion test.
    ///
    /// # Returns
    /// Buoyancy torque vector τ (N·m) in the world frame.
    pub fn compute_buoyancy_torque(
        &self,
        sph_positions: &[[f64; 3]],
        particle_volume: &[f64],
        radius: f64,
    ) -> [f64; 3] {
        let mut torque = [0.0_f64; 3];
        let cm = self.centre_of_mass;
        let rho_f = self.fluid_density;
        let g = self.gravity;
        for (i, &pos) in sph_positions.iter().enumerate() {
            let r = [pos[0] - cm[0], pos[1] - cm[1], pos[2] - cm[2]];
            let dist = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
            if dist > radius {
                continue;
            }
            let d_f = rho_f * g * particle_volume[i];
            let fb = [0.0, 0.0, d_f];
            torque[0] += r[1] * fb[2] - r[2] * fb[1];
            torque[1] += r[2] * fb[0] - r[0] * fb[2];
            torque[2] += r[0] * fb[1] - r[1] * fb[0];
        }
        torque
    }
    /// Compute the total buoyancy force on the rigid body.
    ///
    /// Sums `ρ_f * g * V_p` over all displaced (interior) particles.
    ///
    /// # Returns
    /// Buoyancy force vector F_b (N) — directed along +z (upward).
    pub fn compute_buoyancy_force(
        &self,
        sph_positions: &[[f64; 3]],
        particle_volume: &[f64],
        radius: f64,
    ) -> [f64; 3] {
        let mut total_vol = 0.0_f64;
        let cm = self.centre_of_mass;
        for (i, &pos) in sph_positions.iter().enumerate() {
            let dx = pos[0] - cm[0];
            let dy = pos[1] - cm[1];
            let dz = pos[2] - cm[2];
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist <= radius {
                total_vol += particle_volume[i];
            }
        }
        let fb = self.fluid_density * self.gravity * total_vol;
        [0.0, 0.0, fb]
    }
}
/// Contact model for DEM particle interacting with an infinite planar wall.
///
/// The wall is defined by a point `p0` on the wall and an inward normal `n`.
#[derive(Debug, Clone)]
pub struct DemWallContact {
    /// A point on the wall.
    pub wall_point: [f64; 3],
    /// Inward normal (unit vector pointing into the fluid domain).
    pub wall_normal: [f64; 3],
    /// Normal stiffness (N m⁻¹).
    pub k_n: f64,
    /// Damping coefficient (N s m⁻¹).
    pub gamma_n: f64,
}
impl DemWallContact {
    /// Create a new DEM-wall contact model.
    pub fn new(wall_point: [f64; 3], wall_normal: [f64; 3], k_n: f64, gamma_n: f64) -> Self {
        let mag = norm3(wall_normal).max(1e-30);
        let n = scale3(1.0 / mag, wall_normal);
        Self {
            wall_point,
            wall_normal: n,
            k_n,
            gamma_n,
        }
    }
    /// Signed distance from the wall surface (positive = inside domain).
    pub fn signed_distance(&self, pos: [f64; 3]) -> f64 {
        dot3(sub3(pos, self.wall_point), self.wall_normal)
    }
    /// Repulsion force on a DEM particle from the wall.
    ///
    /// Uses a linear spring-dashpot model.
    pub fn contact_force(&self, dem: &DemSphParticle) -> [f64; 3] {
        let d = self.signed_distance(dem.position);
        let penetration = dem.radius - d;
        if penetration <= 0.0 {
            return [0.0; 3];
        }
        let v_n = dot3(dem.velocity, self.wall_normal);
        let f_n = self.k_n * penetration - self.gamma_n * v_n;
        scale3(f_n.max(0.0), self.wall_normal)
    }
}
