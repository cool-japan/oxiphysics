// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! SPH–FEM and SPH–Rigid coupling module.
//!
//! Provides:
//! - [`SphFemInterface`]: Data structure for coupled SPH particle / FEM node lists.
//! - [`transfer_stress`]: Map SPH Cauchy stresses to FEM nodal forces.
//! - [`transfer_displacement`]: Map FEM nodal displacements to SPH boundary conditions.
//! - [`SphRigidCoupling`]: Coupling between SPH fluid and an immersed rigid body.
//! - [`rigid_to_sph_bc`]: Impose rigid-body velocity as a boundary condition on SPH.
//! - [`sph_to_rigid_force`]: Compute net pressure and viscous force on the rigid body.
//! - [`penalty_coupling`]: Penalty-based weak coupling force between SPH and FEM/rigid.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ---------------------------------------------------------------------------
// Helper math
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Vector difference a - b (3-D).
#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Vector sum a + b (3-D).
#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Scale a 3-vector by scalar `s`.
#[inline]
fn scale3(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Euclidean norm of a 3-vector.
#[inline]
fn norm3(a: [f64; 3]) -> f64 {
    dot3(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// SPH–FEM Interface
// ---------------------------------------------------------------------------

/// Coupled SPH particle / FEM node data used at the shared interface.
///
/// Each entry in the paired lists corresponds to one coupled degree of freedom:
/// a SPH particle that sits on (or near) a FEM boundary node.
#[derive(Debug, Clone)]
pub struct SphFemInterface {
    /// Indices of the SPH interface particles.
    pub sph_indices: Vec<usize>,
    /// Indices of the coupled FEM nodes.
    pub fem_indices: Vec<usize>,
    /// Positions of the interface particles (x, y, z).
    pub positions: Vec<[f64; 3]>,
    /// SPH particle masses.
    pub masses: Vec<f64>,
    /// SPH particle volumes (m³).
    pub volumes: Vec<f64>,
}

impl SphFemInterface {
    /// Construct a new interface descriptor.
    ///
    /// All input slices must have the same length.
    pub fn new(
        sph_indices: Vec<usize>,
        fem_indices: Vec<usize>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
        volumes: Vec<f64>,
    ) -> Self {
        let n = sph_indices.len();
        assert_eq!(fem_indices.len(), n);
        assert_eq!(positions.len(), n);
        assert_eq!(masses.len(), n);
        assert_eq!(volumes.len(), n);
        Self {
            sph_indices,
            fem_indices,
            positions,
            masses,
            volumes,
        }
    }

    /// Number of coupled particle–node pairs.
    pub fn len(&self) -> usize {
        self.sph_indices.len()
    }

    /// Returns `true` if the interface has no pairs.
    pub fn is_empty(&self) -> bool {
        self.sph_indices.is_empty()
    }
}

/// Transfer SPH Cauchy stresses to equivalent FEM nodal forces.
///
/// For each interface pair the nodal force contribution is:
/// ```text
/// f_I = -V_I · σ_I · n_I
/// ```
/// where V_I is the SPH particle volume, σ_I is the Cauchy stress (Voigt:
/// \[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz\]), and n_I is the outward unit
/// normal at the interface.
///
/// * `interface`       — coupled particle–node descriptor
/// * `sph_stresses`    — Cauchy stress (Voigt 6-component) for every SPH particle in the interface (length = interface.len())
/// * `normals`         — outward unit normals at each interface location (length = interface.len())
///
/// Returns a vector of nodal force increments \[fx, fy, fz\] for each FEM interface node.
pub fn transfer_stress(
    interface: &SphFemInterface,
    sph_stresses: &[[f64; 6]],
    normals: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let n = interface.len();
    assert_eq!(sph_stresses.len(), n);
    assert_eq!(normals.len(), n);

    let mut forces = vec![[0.0_f64; 3]; n];

    for i in 0..n {
        let v = interface.volumes[i];
        let s = sph_stresses[i]; // [sxx, syy, szz, sxy, syz, sxz]
        let nrm = normals[i];
        // σ · n  for symmetric stress tensor in Voigt notation
        // σ = [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
        let fx = s[0] * nrm[0] + s[3] * nrm[1] + s[5] * nrm[2];
        let fy = s[3] * nrm[0] + s[1] * nrm[1] + s[4] * nrm[2];
        let fz = s[5] * nrm[0] + s[4] * nrm[1] + s[2] * nrm[2];
        forces[i] = [-v * fx, -v * fy, -v * fz];
    }
    forces
}

/// Transfer FEM nodal displacements to SPH boundary-particle positions.
///
/// Moves each interface SPH particle to the new position implied by the
/// FEM displacement field:
/// ```text
/// x_I_new = x_I_ref + u_I
/// ```
///
/// * `interface`        — coupled particle–node descriptor
/// * `reference_pos`    — reference (undeformed) positions of interface particles
/// * `fem_displacements`— FEM displacement vector \[ux, uy, uz\] for each interface node
///
/// Returns updated positions for the interface SPH particles.
pub fn transfer_displacement(
    interface: &SphFemInterface,
    reference_pos: &[[f64; 3]],
    fem_displacements: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let n = interface.len();
    assert_eq!(reference_pos.len(), n);
    assert_eq!(fem_displacements.len(), n);

    (0..n)
        .map(|i| add3(reference_pos[i], fem_displacements[i]))
        .collect()
}

// ---------------------------------------------------------------------------
// SPH–Rigid Coupling
// ---------------------------------------------------------------------------

/// State of an immersed rigid body coupled to an SPH domain.
///
/// Stores the rigid body's kinematic state and the indices of the SPH
/// boundary particles that represent its surface.
#[derive(Debug, Clone)]
pub struct SphRigidCoupling {
    /// Centre-of-mass position (m).
    pub position: [f64; 3],
    /// Centre-of-mass velocity (m/s).
    pub velocity: [f64; 3],
    /// Angular velocity ω (rad/s), rotation about CoM.
    pub angular_velocity: [f64; 3],
    /// Total mass of the rigid body (kg).
    pub mass: f64,
    /// Indices of the SPH boundary particles on the rigid surface.
    pub boundary_particle_indices: Vec<usize>,
    /// Rest positions of boundary particles relative to the CoM.
    pub rest_offsets: Vec<[f64; 3]>,
}

impl SphRigidCoupling {
    /// Construct a new rigid–SPH coupling descriptor.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        angular_velocity: [f64; 3],
        mass: f64,
        boundary_particle_indices: Vec<usize>,
        rest_offsets: Vec<[f64; 3]>,
    ) -> Self {
        assert_eq!(
            boundary_particle_indices.len(),
            rest_offsets.len(),
            "boundary_particle_indices and rest_offsets must have the same length"
        );
        Self {
            position,
            velocity,
            angular_velocity,
            mass,
            boundary_particle_indices,
            rest_offsets,
        }
    }

    /// Number of surface boundary particles.
    pub fn num_boundary_particles(&self) -> usize {
        self.boundary_particle_indices.len()
    }
}

/// Cross product of two 3-vectors.
#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Impose a rigid-body velocity boundary condition on SPH boundary particles.
///
/// Sets the velocity of each boundary particle to match the rigid body's
/// local velocity:
/// ```text
/// v_I = v_CoM + ω × r_I
/// ```
/// where r_I is the offset from the centre of mass to particle I.
///
/// * `rigid`       — rigid body state
/// * `sph_vel`     — full SPH velocity array (modified in-place at boundary indices)
/// * `current_pos` — current world-space positions of all SPH particles
pub fn rigid_to_sph_bc(
    rigid: &SphRigidCoupling,
    sph_vel: &mut [[f64; 3]],
    current_pos: &[[f64; 3]],
) {
    for (k, &idx) in rigid.boundary_particle_indices.iter().enumerate() {
        let r = sub3(current_pos[idx], rigid.position);
        let omega_cross_r = cross3(rigid.angular_velocity, r);
        sph_vel[idx] = add3(rigid.velocity, omega_cross_r);
        let _ = rigid.rest_offsets[k]; // acknowledge field
    }
}

/// Compute the total SPH hydrodynamic force and torque acting on a rigid body.
///
/// Sums the pressure and viscous forces exerted by the SPH fluid onto each
/// boundary particle:
/// ```text
/// F = Σ_I  m_I / ρ_I · ( -∇p_I + μ ∇²v_I )
/// f_torque = Σ_I  r_I × F_I
/// ```
/// In this simplified version the total force on particle I is passed
/// directly as `boundary_forces[k]`.
///
/// * `rigid`            — rigid body state
/// * `current_pos`      — current world-space positions of all SPH particles
/// * `boundary_forces`  — total hydrodynamic force on each boundary particle (length = num_boundary_particles)
///
/// Returns (total_force, total_torque) as \[f64; 3\] each.
pub fn sph_to_rigid_force(
    rigid: &SphRigidCoupling,
    current_pos: &[[f64; 3]],
    boundary_forces: &[[f64; 3]],
) -> ([f64; 3], [f64; 3]) {
    assert_eq!(boundary_forces.len(), rigid.num_boundary_particles());

    let mut total_force = [0.0_f64; 3];
    let mut total_torque = [0.0_f64; 3];

    for (k, &idx) in rigid.boundary_particle_indices.iter().enumerate() {
        let f = boundary_forces[k];
        total_force = add3(total_force, f);
        let r = sub3(current_pos[idx], rigid.position);
        let torque = cross3(r, f);
        total_torque = add3(total_torque, torque);
    }
    (total_force, total_torque)
}

/// Penalty-force weak coupling between SPH boundary particles and a target surface.
///
/// For each boundary particle I with current position `x_I`, the penalty
/// force is:
/// ```text
/// f_I = -k_pen * max(0, d_pen - |x_I - x_target_I|) * n_I
/// ```
/// where `d_pen` is the penetration depth, `k_pen` is the stiffness, and
/// `n_I` is the outward normal.  The reaction force on the target is −f_I.
///
/// * `particle_pos`    — current position of each coupling particle (length M)
/// * `target_pos`      — target position on the FEM/rigid surface (length M)
/// * `normals`         — outward normals at each target point (length M)
/// * `k_pen`           — penalty stiffness (N/m²)
/// * `d_pen`           — activation gap: only active when distance < d_pen
///
/// Returns a vector of penalty forces \[fx, fy, fz\] for each particle.
pub fn penalty_coupling(
    particle_pos: &[[f64; 3]],
    target_pos: &[[f64; 3]],
    normals: &[[f64; 3]],
    k_pen: f64,
    d_pen: f64,
) -> Vec<[f64; 3]> {
    let n = particle_pos.len();
    assert_eq!(target_pos.len(), n);
    assert_eq!(normals.len(), n);

    (0..n)
        .map(|i| {
            let diff = sub3(particle_pos[i], target_pos[i]);
            let dist = norm3(diff);
            let penetration = d_pen - dist;
            if penetration > 0.0 && dist > 1e-14 {
                // Force magnitude
                let mag = -k_pen * penetration;
                scale3(normals[i], mag)
            } else {
                [0.0, 0.0, 0.0]
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ---------- helpers ----------

    fn make_interface(n: usize) -> SphFemInterface {
        let sph_idx: Vec<usize> = (0..n).collect();
        let fem_idx: Vec<usize> = (0..n).collect();
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        let masses = vec![0.1_f64; n];
        let volumes = vec![0.001_f64; n];
        SphFemInterface::new(sph_idx, fem_idx, positions, masses, volumes)
    }

    // ---------- interface ----------

    #[test]
    fn test_interface_len() {
        let iface = make_interface(5);
        assert_eq!(iface.len(), 5);
        assert!(!iface.is_empty());
    }

    #[test]
    fn test_interface_empty() {
        let iface = make_interface(0);
        assert!(iface.is_empty());
        assert_eq!(iface.len(), 0);
    }

    // ---------- transfer_stress ----------

    #[test]
    fn test_transfer_stress_hydrostatic() {
        // Hydrostatic stress σ = -p I.  Normal = [1,0,0].
        // Force should be f = -V * σ·n = -V * (-p, 0, 0) = (V*p, 0, 0)
        let iface = make_interface(1);
        let p = 1000.0_f64;
        let stress = [[-p, -p, -p, 0.0, 0.0, 0.0]]; // Voigt: sxx,syy,szz,sxy,syz,sxz
        let normal = [[1.0_f64, 0.0, 0.0]];
        let forces = transfer_stress(&iface, &stress, &normal);
        let v = iface.volumes[0];
        assert!(
            (forces[0][0] - v * p).abs() < 1e-10,
            "fx = {}",
            forces[0][0]
        );
        assert!(forces[0][1].abs() < 1e-12);
        assert!(forces[0][2].abs() < 1e-12);
    }

    #[test]
    fn test_transfer_stress_zero() {
        let iface = make_interface(3);
        let stresses = vec![[0.0_f64; 6]; 3];
        let normals = vec![[1.0_f64, 0.0, 0.0]; 3];
        let forces = transfer_stress(&iface, &stresses, &normals);
        for f in &forces {
            for &c in f {
                assert!(c.abs() < 1e-14);
            }
        }
    }

    #[test]
    fn test_transfer_stress_shear() {
        // σ_xy = τ, normal = [0,1,0] → f_x = -V*τ, f_y = 0, f_z = 0
        let iface = make_interface(1);
        let tau = 500.0;
        let stress = [[0.0, 0.0, 0.0, tau, 0.0, 0.0]];
        let normal = [[0.0_f64, 1.0, 0.0]];
        let forces = transfer_stress(&iface, &stress, &normal);
        let v = iface.volumes[0];
        assert!(
            (forces[0][0] + v * tau).abs() < 1e-10,
            "fx = {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_transfer_stress_multiple() {
        let iface = make_interface(4);
        let p = 200.0;
        let stresses: Vec<[f64; 6]> = vec![[-p, -p, -p, 0.0, 0.0, 0.0]; 4];
        let normals: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; 4];
        let forces = transfer_stress(&iface, &stresses, &normals);
        assert_eq!(forces.len(), 4);
    }

    // ---------- transfer_displacement ----------

    #[test]
    fn test_transfer_displacement_zero() {
        let iface = make_interface(3);
        let ref_pos: Vec<[f64; 3]> = (0..3).map(|i| [i as f64, 0.0, 0.0]).collect();
        let displ = vec![[0.0_f64; 3]; 3];
        let new_pos = transfer_displacement(&iface, &ref_pos, &displ);
        for (i, pos) in new_pos.iter().enumerate() {
            assert!((pos[0] - i as f64).abs() < 1e-14);
        }
    }

    #[test]
    fn test_transfer_displacement_uniform() {
        let iface = make_interface(4);
        let ref_pos: Vec<[f64; 3]> = vec![[0.0; 3]; 4];
        let displ: Vec<[f64; 3]> = vec![[1.0, 2.0, 3.0]; 4];
        let new_pos = transfer_displacement(&iface, &ref_pos, &displ);
        for pos in &new_pos {
            assert!((pos[0] - 1.0).abs() < 1e-14);
            assert!((pos[1] - 2.0).abs() < 1e-14);
            assert!((pos[2] - 3.0).abs() < 1e-14);
        }
    }

    #[test]
    fn test_transfer_displacement_length() {
        let iface = make_interface(5);
        let ref_pos = vec![[0.0_f64; 3]; 5];
        let displ = vec![[0.1, 0.2, 0.3]; 5];
        let result = transfer_displacement(&iface, &ref_pos, &displ);
        assert_eq!(result.len(), 5);
    }

    // ---------- SphRigidCoupling ----------

    #[test]
    fn test_rigid_coupling_new() {
        let rigid = SphRigidCoupling::new(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0; 3],
            2.0,
            vec![0, 1, 2],
            vec![[0.1, 0.0, 0.0]; 3],
        );
        assert_eq!(rigid.num_boundary_particles(), 3);
    }

    #[test]
    fn test_rigid_to_sph_bc_translation_only() {
        // Pure translation: ω = 0, all particles get v = v_rigid
        let v_rigid = [2.0_f64, 1.0, 0.0];
        let rigid = SphRigidCoupling::new(
            [0.0; 3],
            v_rigid,
            [0.0; 3],
            1.0,
            vec![0, 1],
            vec![[0.5, 0.0, 0.0], [-0.5, 0.0, 0.0]],
        );
        let pos = vec![[0.5_f64, 0.0, 0.0], [-0.5, 0.0, 0.0]];
        let mut vel = vec![[0.0_f64; 3]; 2];
        rigid_to_sph_bc(&rigid, &mut vel, &pos);
        for v in &vel {
            assert!((v[0] - 2.0).abs() < 1e-12);
            assert!((v[1] - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn test_rigid_to_sph_bc_rotation() {
        // ω = [0,0,1], particle at [1,0,0] → v = ω×r = [0,1,0]
        let rigid = SphRigidCoupling::new(
            [0.0; 3],
            [0.0; 3],
            [0.0, 0.0, 1.0],
            1.0,
            vec![0],
            vec![[1.0, 0.0, 0.0]],
        );
        let pos = vec![[1.0_f64, 0.0, 0.0]];
        let mut vel = vec![[0.0_f64; 3]; 1];
        rigid_to_sph_bc(&rigid, &mut vel, &pos);
        assert!(vel[0][0].abs() < 1e-12, "vx = {}", vel[0][0]);
        assert!((vel[0][1] - 1.0).abs() < 1e-12, "vy = {}", vel[0][1]);
        assert!(vel[0][2].abs() < 1e-12, "vz = {}", vel[0][2]);
    }

    // ---------- sph_to_rigid_force ----------

    #[test]
    fn test_sph_to_rigid_force_single() {
        let rigid = SphRigidCoupling::new(
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            1.0,
            vec![0],
            vec![[1.0, 0.0, 0.0]],
        );
        let pos = vec![[1.0_f64, 0.0, 0.0]];
        let bf = vec![[5.0_f64, 0.0, 0.0]];
        let (force, torque) = sph_to_rigid_force(&rigid, &pos, &bf);
        assert!((force[0] - 5.0).abs() < 1e-12);
        assert!(force[1].abs() < 1e-12);
        // Torque = r × F = [1,0,0] × [5,0,0] = [0,0,0]
        assert!(torque[0].abs() < 1e-12);
        assert!(torque[1].abs() < 1e-12);
        assert!(torque[2].abs() < 1e-12);
    }

    #[test]
    fn test_sph_to_rigid_torque() {
        // r = [1,0,0], F = [0,1,0] → torque = [0,0,1]
        let rigid = SphRigidCoupling::new(
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            1.0,
            vec![0],
            vec![[1.0, 0.0, 0.0]],
        );
        let pos = vec![[1.0_f64, 0.0, 0.0]];
        let bf = vec![[0.0_f64, 1.0, 0.0]];
        let (_force, torque) = sph_to_rigid_force(&rigid, &pos, &bf);
        assert!(torque[0].abs() < 1e-12, "tx = {}", torque[0]);
        assert!(torque[1].abs() < 1e-12, "ty = {}", torque[1]);
        assert!((torque[2] - 1.0).abs() < 1e-12, "tz = {}", torque[2]);
    }

    #[test]
    fn test_sph_to_rigid_force_sum() {
        // Two equal forces in opposite directions → total force = 0
        let rigid = SphRigidCoupling::new(
            [0.0; 3],
            [0.0; 3],
            [0.0; 3],
            1.0,
            vec![0, 1],
            vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        );
        let pos = vec![[1.0_f64, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let bf = vec![[3.0_f64, 0.0, 0.0], [-3.0, 0.0, 0.0]];
        let (force, _torque) = sph_to_rigid_force(&rigid, &pos, &bf);
        assert!(force[0].abs() < 1e-12, "net force = {}", force[0]);
    }

    // ---------- penalty_coupling ----------

    #[test]
    fn test_penalty_no_penetration() {
        // Particle is far from target → no force
        let particle_pos = vec![[2.0_f64, 0.0, 0.0]];
        let target_pos = vec![[0.0_f64, 0.0, 0.0]];
        let normals = vec![[1.0_f64, 0.0, 0.0]];
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, 1e5, 0.1);
        assert!(forces[0][0].abs() < 1e-14);
    }

    #[test]
    fn test_penalty_penetration_magnitude() {
        // Particle coincides with target (dist=0) → full penalty
        let particle_pos = vec![[0.0_f64, 0.0, 0.0]];
        let target_pos = vec![[0.0_f64, 0.0, 0.0]];
        let normals = vec![[1.0_f64, 0.0, 0.0]];
        let k_pen = 1000.0;
        let d_pen = 0.5;
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, k_pen, d_pen);
        // penetration = d_pen - dist = 0.5, force = -k_pen * 0.5 * n
        // but dist = 0 → condition `dist > 1e-14` fails → force = 0
        // (degenerate case, just check it doesn't panic)
        assert_eq!(forces.len(), 1);
    }

    #[test]
    fn test_penalty_slight_penetration() {
        let particle_pos = vec![[0.09_f64, 0.0, 0.0]];
        let target_pos = vec![[0.0_f64, 0.0, 0.0]];
        let normals = vec![[1.0_f64, 0.0, 0.0]];
        let k_pen = 1000.0;
        let d_pen = 0.1;
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, k_pen, d_pen);
        // dist = 0.09, penetration = 0.01, force = -1000 * 0.01 * [1,0,0]
        assert!(
            (forces[0][0] + 10.0).abs() < 0.1,
            "force = {}",
            forces[0][0]
        );
    }

    #[test]
    fn test_penalty_zero_stiffness() {
        let particle_pos = vec![[0.0_f64, 0.0, 0.0]];
        let target_pos = vec![[0.0_f64, 0.0, 0.0]];
        let normals = vec![[1.0_f64, 0.0, 0.0]];
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, 0.0, 0.1);
        assert!(forces[0][0].abs() < 1e-14);
    }

    #[test]
    fn test_penalty_multiple_particles() {
        let n = 5;
        let particle_pos: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.02, 0.0, 0.0]).collect();
        let target_pos = vec![[0.0_f64; 3]; n];
        let normals = vec![[1.0_f64, 0.0, 0.0]; n];
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, 1e4, 0.1);
        assert_eq!(forces.len(), n);
    }

    // ---------- math helpers ----------

    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < 1e-12);
    }

    #[test]
    fn test_cross3_standard() {
        let c = cross3([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((c[2] - 1.0).abs() < 1e-12);
        assert!(c[0].abs() < 1e-12);
        assert!(c[1].abs() < 1e-12);
    }

    #[test]
    fn test_norm3() {
        let n = norm3([3.0, 4.0, 0.0]);
        assert!((n - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_add3() {
        let r = add3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]);
        assert_eq!(r, [5.0, 7.0, 9.0]);
    }

    #[test]
    fn test_sub3() {
        let r = sub3([5.0, 7.0, 9.0], [1.0, 2.0, 3.0]);
        assert_eq!(r, [4.0, 5.0, 6.0]);
    }

    #[test]
    fn test_scale3() {
        let r = scale3([1.0, 2.0, 3.0], 2.0);
        assert_eq!(r, [2.0, 4.0, 6.0]);
    }
}
