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

/// Unit vector along `a`; the zero vector when `a` is zero or not finite.
#[inline]
fn unit3(a: [f64; 3]) -> [f64; 3] {
    let len = norm3(a);
    if len.is_finite() && len > 0.0 {
        scale3(a, 1.0 / len)
    } else {
        [0.0; 3]
    }
}

// ---------------------------------------------------------------------------
// SPH–FEM Interface
// ---------------------------------------------------------------------------

/// Coupled SPH particle / FEM node data used at the shared interface.
///
/// Each entry in the paired lists corresponds to one coupled degree of freedom:
/// a SPH particle that sits on (or near) a FEM boundary node.
///
/// A stress transferred to a FEM node becomes a force only once it is
/// multiplied by the *area* of boundary surface the node (particle) is
/// responsible for, which is why the interface carries `areas` next to the
/// particle `volumes`.
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
    /// Tributary contact area of every interface particle (m²): the part of
    /// the FEM boundary surface on which the particle's traction acts.  A
    /// particle of volume `V` and spacing `dx` on a flat interface covers
    /// `A = V / dx` (see [`SphFemInterface::with_spacing`]).
    pub areas: Vec<f64>,
}

impl SphFemInterface {
    /// Construct a new interface descriptor.
    ///
    /// All input vectors must have the same length, and every area must be
    /// finite and non-negative.
    ///
    /// # Panics
    ///
    /// Panics if the vectors differ in length or if an entry of `areas` is
    /// negative or not finite.
    pub fn new(
        sph_indices: Vec<usize>,
        fem_indices: Vec<usize>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
        volumes: Vec<f64>,
        areas: Vec<f64>,
    ) -> Self {
        let n = sph_indices.len();
        assert_eq!(fem_indices.len(), n);
        assert_eq!(positions.len(), n);
        assert_eq!(masses.len(), n);
        assert_eq!(volumes.len(), n);
        assert_eq!(areas.len(), n);
        assert!(
            areas.iter().all(|a| a.is_finite() && *a >= 0.0),
            "interface areas must be finite and non-negative"
        );
        Self {
            sph_indices,
            fem_indices,
            positions,
            masses,
            volumes,
            areas,
        }
    }

    /// Construct an interface whose tributary areas follow from the particle
    /// spacing: `A_I = V_I / dx` (for a cubic particle cell of edge `dx` this
    /// is its face area `dx²`).
    ///
    /// # Panics
    ///
    /// Panics if the vectors differ in length or if `dx` is not a finite,
    /// strictly positive length.
    pub fn with_spacing(
        sph_indices: Vec<usize>,
        fem_indices: Vec<usize>,
        positions: Vec<[f64; 3]>,
        masses: Vec<f64>,
        volumes: Vec<f64>,
        dx: f64,
    ) -> Self {
        assert!(
            dx.is_finite() && dx > 0.0,
            "particle spacing dx must be finite and positive, got {dx}"
        );
        let areas = volumes.iter().map(|v| v / dx).collect();
        Self::new(sph_indices, fem_indices, positions, masses, volumes, areas)
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
/// f_I = -A_I · σ_I · n_I
/// ```
/// where A_I is the contact area (m²) the interface particle represents
/// ([`SphFemInterface::areas`]), σ_I is the Cauchy stress (Voigt:
/// \[σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz\], Pa) and n_I is the unit outward
/// normal of the fluid at the interface, i.e. pointing from the fluid into the
/// structure.  `σ_I · n_I` is the traction (Pa) the structure exerts on the
/// fluid; multiplied by the area it is a force in N, and the minus sign gives
/// the reaction force that the fluid exerts on the FEM node.  (A hydrostatic
/// fluid, σ = −p I, therefore pushes the structure along +n with `p·A`.)
///
/// * `interface`       — coupled particle–node descriptor
/// * `sph_stresses`    — Cauchy stress (Voigt 6-component) for every SPH particle in the interface (length = interface.len())
/// * `normals`         — outward normals at each interface location (length = interface.len());
///   they are normalised internally and a zero normal contributes no force
///
/// Returns a vector of nodal forces \[fx, fy, fz\] (N) for each FEM interface node.
///
/// # Panics
///
/// Panics if the slices do not all have `interface.len()` entries.
pub fn transfer_stress(
    interface: &SphFemInterface,
    sph_stresses: &[[f64; 6]],
    normals: &[[f64; 3]],
) -> Vec<[f64; 3]> {
    let n = interface.len();
    assert_eq!(sph_stresses.len(), n);
    assert_eq!(normals.len(), n);
    assert_eq!(interface.areas.len(), n);

    let mut forces = vec![[0.0_f64; 3]; n];

    for i in 0..n {
        let a = interface.areas[i];
        let s = sph_stresses[i]; // [sxx, syy, szz, sxy, syz, sxz]
        let nrm = unit3(normals[i]);
        // σ · n  for symmetric stress tensor in Voigt notation
        // σ = [[sxx, sxy, sxz], [sxy, syy, syz], [sxz, syz, szz]]
        let fx = s[0] * nrm[0] + s[3] * nrm[1] + s[5] * nrm[2];
        let fy = s[3] * nrm[0] + s[1] * nrm[1] + s[4] * nrm[2];
        let fz = s[5] * nrm[0] + s[4] * nrm[1] + s[2] * nrm[2];
        forces[i] = [-a * fx, -a * fy, -a * fz];
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
/// For each boundary particle I the signed gap to the target surface is
/// measured along the outward normal `n_I` of the target:
/// ```text
/// g_I = (x_I - x_target_I) · n_I
/// ```
/// `g_I > 0` means the particle is on the outer (fluid) side of the surface and
/// `g_I < 0` means it has crossed the surface.  Inside the activation gap the
/// penalty force pushes the particle back out along the normal:
/// ```text
/// f_I = k_pen * max(0, d_pen - g_I) * n_I
/// ```
/// where `d_pen` is the activation gap and `k_pen` is the stiffness.  The force
/// is repulsive (`f_I · n_I >= 0`), grows linearly as the particle gets closer
/// to the surface, stays finite when the particle sits exactly on it (`g_I = 0`
/// gives `k_pen * d_pen`) and keeps growing for particles that have crossed it
/// (`g_I < 0`).  The reaction force on the target is −f_I.
///
/// * `particle_pos`    — current position of each coupling particle (length M)
/// * `target_pos`      — target position on the FEM/rigid surface (length M)
/// * `normals`         — outward normals at each target point (length M); they
///   are normalised internally and a zero normal yields a zero force
/// * `k_pen`           — penalty stiffness (N/m), non-negative
/// * `d_pen`           — activation gap (m): the force acts only while the
///   signed gap is smaller than `d_pen`
///
/// Returns a vector of penalty forces \[fx, fy, fz\] (N) for each particle.
///
/// # Panics
///
/// Panics if the slices have different lengths.
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
            let nrm = unit3(normals[i]);
            let gap = dot3(sub3(particle_pos[i], target_pos[i]), nrm);
            if gap < d_pen {
                scale3(nrm, k_pen * (d_pen - gap))
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

    /// Interface of `n` particles with volume 1e-3 m³ and contact area `area`.
    fn make_interface_with_area(n: usize, area: f64) -> SphFemInterface {
        let sph_idx: Vec<usize> = (0..n).collect();
        let fem_idx: Vec<usize> = (0..n).collect();
        let positions: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        let masses = vec![0.1_f64; n];
        let volumes = vec![0.001_f64; n];
        let areas = vec![area; n];
        SphFemInterface::new(sph_idx, fem_idx, positions, masses, volumes, areas)
    }

    /// Interface of `n` particles with a 1e-4 m² (1 cm x 1 cm) contact area each.
    fn make_interface(n: usize) -> SphFemInterface {
        make_interface_with_area(n, 1.0e-4)
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

    #[test]
    fn test_interface_keeps_areas() {
        let iface = make_interface_with_area(3, 2.5e-5);
        assert_eq!(iface.areas, vec![2.5e-5; 3]);
    }

    #[test]
    fn test_with_spacing_derives_areas_from_volume_and_spacing() {
        let dx = 0.02_f64;
        let cell = dx * dx * dx;
        let iface = SphFemInterface::with_spacing(
            vec![0, 1],
            vec![0, 1],
            vec![[0.0; 3]; 2],
            vec![0.1; 2],
            vec![cell, 2.0 * cell],
            dx,
        );
        assert!(
            (iface.areas[0] - dx * dx).abs() < 1e-18,
            "{}",
            iface.areas[0]
        );
        assert!(
            (iface.areas[1] - 2.0 * dx * dx).abs() < 1e-18,
            "{}",
            iface.areas[1]
        );
        assert_eq!(iface.volumes, vec![cell, 2.0 * cell]);
    }

    #[test]
    #[should_panic(expected = "particle spacing")]
    fn test_with_spacing_rejects_zero_spacing() {
        let _ = SphFemInterface::with_spacing(
            vec![0],
            vec![0],
            vec![[0.0; 3]],
            vec![0.1],
            vec![1.0e-3],
            0.0,
        );
    }

    #[test]
    #[should_panic(expected = "particle spacing")]
    fn test_with_spacing_rejects_negative_spacing() {
        let _ = SphFemInterface::with_spacing(
            vec![0],
            vec![0],
            vec![[0.0; 3]],
            vec![0.1],
            vec![1.0e-3],
            -0.01,
        );
    }

    #[test]
    #[should_panic(expected = "particle spacing")]
    fn test_with_spacing_rejects_non_finite_spacing() {
        let _ = SphFemInterface::with_spacing(
            vec![0],
            vec![0],
            vec![[0.0; 3]],
            vec![0.1],
            vec![1.0e-3],
            f64::NAN,
        );
    }

    #[test]
    #[should_panic(expected = "non-negative")]
    fn test_new_rejects_negative_area() {
        let _ = SphFemInterface::new(
            vec![0],
            vec![0],
            vec![[0.0; 3]],
            vec![0.1],
            vec![1.0e-3],
            vec![-1.0e-4],
        );
    }

    #[test]
    #[should_panic(expected = "non-negative")]
    fn test_new_rejects_non_finite_area() {
        let _ = SphFemInterface::new(
            vec![0],
            vec![0],
            vec![[0.0; 3]],
            vec![0.1],
            vec![1.0e-3],
            vec![f64::INFINITY],
        );
    }

    #[test]
    #[should_panic]
    fn test_new_rejects_mismatched_area_length() {
        let _ = SphFemInterface::new(
            vec![0, 1],
            vec![0, 1],
            vec![[0.0; 3]; 2],
            vec![0.1; 2],
            vec![1.0e-3; 2],
            vec![1.0e-4],
        );
    }

    // ---------- transfer_stress ----------

    #[test]
    fn test_transfer_stress_hydrostatic() {
        // Hydrostatic stress σ = -p I, fluid-outward normal n = +x.
        // f = -A σ·n = A p (+x): a 1000 Pa traction on 1e-4 m² is 0.1 N.
        let iface = make_interface(1);
        let p = 1000.0_f64;
        let stress = [[-p, -p, -p, 0.0, 0.0, 0.0]]; // Voigt: sxx,syy,szz,sxy,syz,sxz
        let normal = [[1.0_f64, 0.0, 0.0]];
        let forces = transfer_stress(&iface, &stress, &normal);
        assert!(
            (forces[0][0] - 0.1).abs() < 1e-12,
            "fx = {} (expected A p = 0.1 N)",
            forces[0][0]
        );
        assert!(forces[0][1].abs() < 1e-12);
        assert!(forces[0][2].abs() < 1e-12);
    }

    #[test]
    fn test_transfer_stress_is_a_force_not_a_traction_times_volume() {
        // A traction (Pa) times a particle volume (m³) is N·m; the nodal force must
        // depend on the contact area only, never on the particle volume.
        let mut small = make_interface(1);
        let mut large = make_interface(1);
        small.volumes[0] = 1.0e-3;
        large.volumes[0] = 8.0e-3;
        let stress = [[-1000.0, -1000.0, -1000.0, 0.0, 0.0, 0.0]];
        let normal = [[1.0_f64, 0.0, 0.0]];
        let f_small = transfer_stress(&small, &stress, &normal);
        let f_large = transfer_stress(&large, &stress, &normal);
        assert_eq!(f_small, f_large);
        assert!((f_small[0][0] - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_transfer_stress_scales_linearly_with_area() {
        let stress = [[-500.0, -500.0, -500.0, 0.0, 0.0, 0.0]];
        let normal = [[0.0_f64, 0.0, 1.0]];
        let f1 = transfer_stress(&make_interface_with_area(1, 1.0e-4), &stress, &normal);
        let f3 = transfer_stress(&make_interface_with_area(1, 3.0e-4), &stress, &normal);
        assert!((f1[0][2] - 0.05).abs() < 1e-12, "fz = {}", f1[0][2]);
        assert!((f3[0][2] / f1[0][2] - 3.0).abs() < 1e-12);
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
        // σ_xy = τ, normal = [0,1,0] → f_x = -A τ = -0.05 N, f_y = 0, f_z = 0
        let iface = make_interface(1);
        let tau = 500.0;
        let stress = [[0.0, 0.0, 0.0, tau, 0.0, 0.0]];
        let normal = [[0.0_f64, 1.0, 0.0]];
        let forces = transfer_stress(&iface, &stress, &normal);
        let a = iface.areas[0];
        assert!(
            (forces[0][0] + a * tau).abs() < 1e-12,
            "fx = {}",
            forces[0][0]
        );
        assert!((forces[0][0] + 0.05).abs() < 1e-12);
        assert!(forces[0][1].abs() < 1e-12);
        assert!(forces[0][2].abs() < 1e-12);
    }

    #[test]
    fn test_transfer_stress_full_tensor_contraction() {
        // σ = [[1,2,3],[2,4,5],[3,5,6]] (Voigt: 1,4,6,2,5,3), n = z.
        // σ·n = (3,5,6)  →  f = -A (3,5,6)
        let iface = make_interface(1);
        let stress = [[1.0, 4.0, 6.0, 2.0, 5.0, 3.0]];
        let normal = [[0.0_f64, 0.0, 1.0]];
        let f = transfer_stress(&iface, &stress, &normal)[0];
        let a = iface.areas[0];
        assert!((f[0] + 3.0 * a).abs() < 1e-15, "fx = {}", f[0]);
        assert!((f[1] + 5.0 * a).abs() < 1e-15, "fy = {}", f[1]);
        assert!((f[2] + 6.0 * a).abs() < 1e-15, "fz = {}", f[2]);
    }

    #[test]
    fn test_transfer_stress_normalises_the_normal() {
        let iface = make_interface(1);
        let stress = [[-1000.0, -1000.0, -1000.0, 0.0, 0.0, 0.0]];
        let unit = transfer_stress(&iface, &stress, &[[0.6, 0.8, 0.0]]);
        let scaled = transfer_stress(&iface, &stress, &[[3.0, 4.0, 0.0]]);
        for k in 0..3 {
            assert!((unit[0][k] - scaled[0][k]).abs() < 1e-14);
        }
        assert!((unit[0][0] - 0.06).abs() < 1e-12);
        assert!((unit[0][1] - 0.08).abs() < 1e-12);
    }

    #[test]
    fn test_transfer_stress_zero_normal_gives_zero_force() {
        let iface = make_interface(1);
        let stress = [[-1000.0, -1000.0, -1000.0, 10.0, 10.0, 10.0]];
        let forces = transfer_stress(&iface, &stress, &[[0.0, 0.0, 0.0]]);
        assert_eq!(forces[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_transfer_stress_multiple() {
        let iface = make_interface(4);
        let p = 200.0;
        let stresses: Vec<[f64; 6]> = vec![[-p, -p, -p, 0.0, 0.0, 0.0]; 4];
        let normals: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0]; 4];
        let forces = transfer_stress(&iface, &stresses, &normals);
        assert_eq!(forces.len(), 4);
        for f in &forces {
            assert!((f[0] - 0.02).abs() < 1e-12, "fx = {}", f[0]);
        }
    }

    #[test]
    fn test_transfer_stress_resultant_is_pressure_times_total_area() {
        // 5 x 5 particles of edge dx = 1 cm on a flat plate loaded with p = 2e5 Pa:
        // the resultant must be p · (5 dx)² = 500 N.
        let dx = 0.01_f64;
        let n = 25;
        let iface = SphFemInterface::with_spacing(
            (0..n).collect(),
            (0..n).collect(),
            vec![[0.0; 3]; n],
            vec![1.0e-3; n],
            vec![dx * dx * dx; n],
            dx,
        );
        let p = 2.0e5_f64;
        let stresses = vec![[-p, -p, -p, 0.0, 0.0, 0.0]; n];
        let normals = vec![[0.0_f64, -1.0, 0.0]; n];
        let forces = transfer_stress(&iface, &stresses, &normals);
        let fy: f64 = forces.iter().map(|f| f[1]).sum();
        assert!((fy + 500.0).abs() < 1e-9, "resultant fy = {fy}");
        let fx: f64 = forces.iter().map(|f| f[0]).sum();
        assert!(fx.abs() < 1e-9);
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

    /// Single particle at `offset` from a target at the origin with normal +x.
    fn penalty_force_at(offset: [f64; 3], k_pen: f64, d_pen: f64) -> [f64; 3] {
        let forces = penalty_coupling(
            &[offset],
            &[[0.0, 0.0, 0.0]],
            &[[1.0, 0.0, 0.0]],
            k_pen,
            d_pen,
        );
        forces[0]
    }

    #[test]
    fn test_penalty_no_penetration() {
        // Particle is far from target → no force
        let particle_pos = vec![[2.0_f64, 0.0, 0.0]];
        let target_pos = vec![[0.0_f64, 0.0, 0.0]];
        let normals = vec![[1.0_f64, 0.0, 0.0]];
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, 1e5, 0.1);
        assert_eq!(forces[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_penalty_gap_equal_to_activation_distance_gives_zero_force() {
        let f = penalty_force_at([0.1, 0.0, 0.0], 1000.0, 0.1);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_penalty_coincident_particle_gets_full_activation_force() {
        // The particle sits exactly on the surface (g = 0): f = k_pen * d_pen * n.
        let f = penalty_force_at([0.0, 0.0, 0.0], 1000.0, 0.5);
        assert!((f[0] - 500.0).abs() < 1e-9, "fx = {}", f[0]);
        assert!(f[1].abs() < 1e-12 && f[2].abs() < 1e-12);
    }

    #[test]
    fn test_penalty_slight_penetration() {
        // g = 0.09 < d_pen = 0.1: f = k_pen * (d_pen - g) * n = +10 N along the normal.
        let f = penalty_force_at([0.09, 0.0, 0.0], 1000.0, 0.1);
        assert!((f[0] - 10.0).abs() < 1e-9, "force = {}", f[0]);
    }

    #[test]
    fn test_penalty_crossed_surface_is_pushed_back_harder() {
        // A particle that has crossed the surface (g = -0.05) must feel a larger
        // repulsive force than one that is still outside (g = +0.05).
        let outside = penalty_force_at([0.05, 0.0, 0.0], 1000.0, 0.1);
        let crossed = penalty_force_at([-0.05, 0.0, 0.0], 1000.0, 0.1);
        assert!((outside[0] - 50.0).abs() < 1e-9, "outside = {}", outside[0]);
        assert!(
            (crossed[0] - 150.0).abs() < 1e-9,
            "crossed = {}",
            crossed[0]
        );
        assert!(crossed[0] > outside[0]);
    }

    #[test]
    fn test_penalty_force_is_repulsive_and_linear_in_gap() {
        let k_pen = 2.0e3;
        let d_pen = 0.2;
        let mut previous = f64::NEG_INFINITY;
        for step in (0..=12).rev() {
            let gap = -0.1 + 0.025 * step as f64; // from +0.2 down to -0.1
            let f = penalty_force_at([gap, 0.0, 0.0], k_pen, d_pen);
            assert!(
                f[0] >= 0.0,
                "force must point along +n, got {} at g = {gap}",
                f[0]
            );
            assert!(f[0] >= previous, "force must grow as the gap closes");
            previous = f[0];
            let expected = k_pen * (d_pen - gap).max(0.0);
            assert!(
                (f[0] - expected).abs() < 1e-9,
                "g = {gap}: {} vs {expected}",
                f[0]
            );
        }
    }

    #[test]
    fn test_penalty_tangential_offset_does_not_change_the_force() {
        // Only the normal component of the separation is a gap.
        let on_axis = penalty_force_at([0.04, 0.0, 0.0], 1000.0, 0.1);
        let sheared = penalty_force_at([0.04, 3.0, -2.0], 1000.0, 0.1);
        assert_eq!(on_axis, sheared);
        assert!((on_axis[0] - 60.0).abs() < 1e-9);
    }

    #[test]
    fn test_penalty_oblique_normal() {
        // Normal (1,1,0)/√2, particle displaced by 0.05 along it: g = 0.05, so
        // f = 1000 * 0.05 = 50 N along the normal.
        let s = std::f64::consts::FRAC_1_SQRT_2;
        let n = [s, s, 0.0];
        let offset = [0.05 * s, 0.05 * s, 0.0];
        let forces = penalty_coupling(&[offset], &[[0.0; 3]], &[n], 1000.0, 0.1);
        assert!((forces[0][0] - 50.0 * s).abs() < 1e-9);
        assert!((forces[0][1] - 50.0 * s).abs() < 1e-9);
        assert!(forces[0][2].abs() < 1e-12);
    }

    #[test]
    fn test_penalty_normalises_the_normal() {
        let unit = penalty_coupling(
            &[[0.04, 0.0, 0.0]],
            &[[0.0; 3]],
            &[[1.0, 0.0, 0.0]],
            1000.0,
            0.1,
        );
        let long = penalty_coupling(
            &[[0.04, 0.0, 0.0]],
            &[[0.0; 3]],
            &[[5.0, 0.0, 0.0]],
            1000.0,
            0.1,
        );
        assert_eq!(unit, long);
        assert!((unit[0][0] - 60.0).abs() < 1e-9);
    }

    #[test]
    fn test_penalty_zero_normal_gives_zero_force() {
        let forces = penalty_coupling(&[[0.0; 3]], &[[0.0; 3]], &[[0.0; 3]], 1000.0, 0.1);
        assert_eq!(forces[0], [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_penalty_zero_stiffness() {
        let f = penalty_force_at([0.0, 0.0, 0.0], 0.0, 0.1);
        assert!(f[0].abs() < 1e-14);
    }

    #[test]
    fn test_penalty_multiple_particles() {
        let n = 5;
        let particle_pos: Vec<[f64; 3]> = (0..n).map(|i| [i as f64 * 0.02, 0.0, 0.0]).collect();
        let target_pos = vec![[0.0_f64; 3]; n];
        let normals = vec![[1.0_f64, 0.0, 0.0]; n];
        let forces = penalty_coupling(&particle_pos, &target_pos, &normals, 1e4, 0.1);
        assert_eq!(forces.len(), n);
        // gaps 0, 0.02, 0.04, 0.06, 0.08 → forces 1000, 800, 600, 400, 200 N
        for (i, f) in forces.iter().enumerate() {
            let expected = 1.0e4 * (0.1 - 0.02 * i as f64);
            assert!(
                (f[0] - expected).abs() < 1e-9,
                "particle {i}: {} vs {expected}",
                f[0]
            );
        }
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
