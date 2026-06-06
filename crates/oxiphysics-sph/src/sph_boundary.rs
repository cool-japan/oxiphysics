// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Advanced SPH boundary conditions.
//!
//! Implements:
//! - [`BoundaryParticle`] — boundary particle with type, normal, and tangent
//! - \[`mirror_ghost_particle()`\] — create a mirror ghost for free-surface tracking
//! - \[`solid_boundary_force()`\] — Lennard-Jones repulsive boundary force
//! - \[`inflow_boundary()`\] — buffer zone particle inlet
//! - \[`outflow_boundary()`\] — absorbing outflow condition (mark particles for removal)
//! - \[`periodic_boundary()`\] — copy particles across periodic faces
//! - \[`no_penetration_correction()`\] — velocity correction at solid walls

// ── Math helpers ─────────────────────────────────────────────────────────────

/// Dot product of two 3-D vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Squared Euclidean distance between two 3-D points.
fn sq_dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
}

/// Euclidean distance between two 3-D points.
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    sq_dist3(a, b).sqrt()
}

/// Normalize a 3-D vector; returns the zero vector if the norm is negligible.
fn normalize3(v: [f64; 3]) -> [f64; 3] {
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag < 1e-300 {
        [0.0, 0.0, 0.0]
    } else {
        [v[0] / mag, v[1] / mag, v[2] / mag]
    }
}

/// Scale a 3-D vector.
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

/// Add two 3-D vectors.
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// Subtract two 3-D vectors.
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

// ═══════════════════════════════════════════════════════════════════════════
// § 1  BoundaryType
// ═══════════════════════════════════════════════════════════════════════════

/// Classification of a boundary particle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryType {
    /// Rigid solid wall.
    SolidWall,
    /// Free-surface particle.
    FreeSurface,
    /// Inlet buffer zone.
    Inlet,
    /// Outlet / absorbing layer.
    Outlet,
    /// Periodic face.
    Periodic,
}

// ═══════════════════════════════════════════════════════════════════════════
// § 2  BoundaryParticle
// ═══════════════════════════════════════════════════════════════════════════

/// A particle that represents a boundary condition.
///
/// The `normal` vector points **into** the fluid domain; `tangent` lies along
/// the boundary surface.
#[derive(Debug, Clone)]
pub struct BoundaryParticle {
    /// Position \[x, y, z\].
    pub pos: [f64; 3],
    /// Outward-pointing (fluid-side) unit normal.
    pub normal: [f64; 3],
    /// Unit tangent vector along the boundary.
    pub tangent: [f64; 3],
    /// Boundary classification.
    pub kind: BoundaryType,
    /// Particle mass.
    pub mass: f64,
    /// Smoothing length.
    pub h: f64,
    /// Prescribed velocity (used for inlet/moving walls).
    pub velocity: [f64; 3],
}

impl BoundaryParticle {
    /// Create a new [`BoundaryParticle`].
    pub fn new(
        pos: [f64; 3],
        normal: [f64; 3],
        tangent: [f64; 3],
        kind: BoundaryType,
        mass: f64,
        h: f64,
        velocity: [f64; 3],
    ) -> Self {
        Self {
            pos,
            normal: normalize3(normal),
            tangent: normalize3(tangent),
            kind,
            mass,
            h,
            velocity,
        }
    }
}

// ═══════════════════════════════════════════════════════════════════════════
// § 3  mirror_ghost_particle
// ═══════════════════════════════════════════════════════════════════════════

/// Create a mirror ghost particle for a fluid particle near a free surface or
/// solid wall.
///
/// The ghost is the reflection of `fluid_pos` across the plane defined by
/// `boundary_pos` and `normal`.  Its velocity is set to `ghost_vel`
/// (caller-supplied, e.g. the reflection of the fluid velocity).
///
/// Returns `(ghost_pos, ghost_vel)`.
pub fn mirror_ghost_particle(
    fluid_pos: [f64; 3],
    fluid_vel: [f64; 3],
    boundary_pos: [f64; 3],
    normal: [f64; 3],
    reflect_vel: bool,
) -> ([f64; 3], [f64; 3]) {
    let n = normalize3(normal);
    // signed distance from the boundary plane
    let d = dot3(sub3(fluid_pos, boundary_pos), n);
    // reflected position
    let ghost_pos = sub3(fluid_pos, scale3(n, 2.0 * d));
    // velocity: optionally reflect the normal component (for no-slip/free-slip)
    let ghost_vel = if reflect_vel {
        let vn = dot3(fluid_vel, n);
        sub3(fluid_vel, scale3(n, 2.0 * vn))
    } else {
        fluid_vel
    };
    (ghost_pos, ghost_vel)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 4  solid_boundary_force
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the Lennard-Jones boundary repulsion force on a fluid particle
/// approaching a solid boundary particle.
///
/// The force magnitude follows
///
/// `f = strength * [(r0/r)^p1 − (r0/r)^p2] * (1/r²) * n`
///
/// where `r` is the distance between the fluid and boundary particles,
/// `r0` is the reference distance, and `(p1, p2) = (4, 2)` are the
/// Lennard-Jones exponents.
///
/// Returns the force vector `[fx, fy, fz]` acting on the fluid particle.
pub fn solid_boundary_force(
    fluid_pos: [f64; 3],
    boundary_pos: [f64; 3],
    r0: f64,
    strength: f64,
) -> [f64; 3] {
    let r = dist3(fluid_pos, boundary_pos);
    if r < 1e-10 || r > 2.0 * r0 {
        return [0.0, 0.0, 0.0];
    }
    let ratio = r0 / r;
    let magnitude = strength * (ratio.powi(4) - ratio.powi(2)) / (r * r);
    let dir = normalize3(sub3(fluid_pos, boundary_pos));
    scale3(dir, magnitude)
}

// ═══════════════════════════════════════════════════════════════════════════
// § 5  inflow_boundary
// ═══════════════════════════════════════════════════════════════════════════

/// Generate inlet particles in a buffer zone.
///
/// Creates up to `max_particles` new particles at positions staggered along
/// the inflow normal with spacing `spacing`.  Each new particle is given
/// `inflow_vel` and `mass`.
///
/// Returns the list of newly created particles as `(position, velocity, mass)` tuples.
pub fn inflow_boundary(
    inlet_center: [f64; 3],
    inflow_normal: [f64; 3],
    inflow_vel: [f64; 3],
    mass: f64,
    spacing: f64,
    max_particles: usize,
) -> Vec<([f64; 3], [f64; 3], f64)> {
    let n = normalize3(inflow_normal);
    let half = (max_particles as f64 - 1.0) * 0.5;
    (0..max_particles)
        .map(|i| {
            let offset = (i as f64 - half) * spacing;
            let pos = add3(inlet_center, scale3(n, offset));
            (pos, inflow_vel, mass)
        })
        .collect()
}

// ═══════════════════════════════════════════════════════════════════════════
// § 6  outflow_boundary
// ═══════════════════════════════════════════════════════════════════════════

/// Mark particles for removal (outflow) that have crossed the outlet plane.
///
/// The outlet plane is defined by `outlet_pos` and `outlet_normal` (pointing
/// outward / downstream).  Returns the indices of particles that have crossed
/// the plane (`dot(x − x_outlet, n) ≥ 0`).
pub fn outflow_boundary(
    positions: &[[f64; 3]],
    outlet_pos: [f64; 3],
    outlet_normal: [f64; 3],
) -> Vec<usize> {
    let n = normalize3(outlet_normal);
    positions
        .iter()
        .enumerate()
        .filter(|&(_, pos)| dot3(sub3(*pos, outlet_pos), n) >= 0.0)
        .map(|(i, _)| i)
        .collect()
}

// ═══════════════════════════════════════════════════════════════════════════
// § 7  periodic_boundary
// ═══════════════════════════════════════════════════════════════════════════

/// Copy particles across a periodic boundary.
///
/// For each particle within `buffer_width` of either end of the periodic
/// axis, a ghost copy is created at the opposite end.
///
/// # Arguments
/// * `positions` — particle positions (arbitrary dimension, but only axis `axis`
///   is wrapped)
/// * `velocities` — corresponding velocities
/// * `masses` — corresponding masses
/// * `axis` — axis along which periodicity applies (0, 1, or 2)
/// * `period` — domain length along `axis`
/// * `buffer_width` — thickness of the copy buffer near each face
///
/// Returns ghost `(position, velocity, mass)` tuples.
pub fn periodic_boundary(
    positions: &[[f64; 3]],
    velocities: &[[f64; 3]],
    masses: &[f64],
    axis: usize,
    period: f64,
    buffer_width: f64,
) -> Vec<([f64; 3], [f64; 3], f64)> {
    let mut ghosts = Vec::new();
    for (idx, &pos) in positions.iter().enumerate() {
        let x = pos[axis];
        // near lower face → ghost at upper
        if x < buffer_width {
            let mut ghost_pos = pos;
            ghost_pos[axis] += period;
            ghosts.push((ghost_pos, velocities[idx], masses[idx]));
        }
        // near upper face → ghost at lower
        if x > period - buffer_width {
            let mut ghost_pos = pos;
            ghost_pos[axis] -= period;
            ghosts.push((ghost_pos, velocities[idx], masses[idx]));
        }
    }
    ghosts
}

// ═══════════════════════════════════════════════════════════════════════════
// § 8  no_penetration_correction
// ═══════════════════════════════════════════════════════════════════════════

/// Apply a no-penetration velocity correction at a solid wall.
///
/// Removes the normal component of `vel` directed into the wall.  Optionally
/// applies a free-slip condition (tangential component preserved) or a
/// partial-slip / no-slip damping.
///
/// # Arguments
/// * `vel` — fluid particle velocity
/// * `wall_normal` — unit normal pointing **into** the fluid (away from wall)
/// * `slip_factor` — 0 = no-slip (full tangential damping), 1 = free-slip
///
/// Returns the corrected velocity.
pub fn no_penetration_correction(
    vel: [f64; 3],
    wall_normal: [f64; 3],
    slip_factor: f64,
) -> [f64; 3] {
    let n = normalize3(wall_normal);
    let vn = dot3(vel, n);
    // Remove penetrating component
    let vel_corrected = if vn < 0.0 {
        sub3(vel, scale3(n, vn))
    } else {
        vel
    };
    // Apply slip factor: blend tangential component
    let vn2 = dot3(vel_corrected, n);
    let v_normal = scale3(n, vn2);
    let v_tangent = sub3(vel_corrected, v_normal);
    add3(v_normal, scale3(v_tangent, slip_factor))
}

// ═══════════════════════════════════════════════════════════════════════════
// § 9  Additional helpers
// ═══════════════════════════════════════════════════════════════════════════

/// Compute the signed distance from `point` to the plane defined by
/// `plane_pos` and `plane_normal` (unit normal).
pub fn signed_distance_to_plane(
    point: [f64; 3],
    plane_pos: [f64; 3],
    plane_normal: [f64; 3],
) -> f64 {
    dot3(sub3(point, plane_pos), normalize3(plane_normal))
}

/// Return the index of the nearest boundary particle to `query_pos`, or
/// `None` if `boundaries` is empty.
pub fn nearest_boundary(query_pos: [f64; 3], boundaries: &[BoundaryParticle]) -> Option<usize> {
    boundaries
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            dist3(query_pos, a.pos)
                .partial_cmp(&dist3(query_pos, b.pos))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(i, _)| i)
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests
// ═══════════════════════════════════════════════════════════════════════════

#[cfg(test)]
mod tests {
    use super::*;

    // ── BoundaryParticle tests ────────────────────────────────────────────

    #[test]
    fn test_boundary_particle_new() {
        let bp = BoundaryParticle::new(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            BoundaryType::SolidWall,
            1.0,
            0.1,
            [0.0; 3],
        );
        assert_eq!(bp.kind, BoundaryType::SolidWall);
        assert!((bp.normal[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boundary_particle_normal_normalized() {
        let bp = BoundaryParticle::new(
            [0.0, 0.0, 0.0],
            [3.0, 4.0, 0.0], // not unit
            [1.0, 0.0, 0.0],
            BoundaryType::FreeSurface,
            1.0,
            0.1,
            [0.0; 3],
        );
        let mag = (bp.normal[0].powi(2) + bp.normal[1].powi(2) + bp.normal[2].powi(2)).sqrt();
        assert!((mag - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_boundary_type_variants() {
        assert_ne!(BoundaryType::SolidWall, BoundaryType::FreeSurface);
        assert_ne!(BoundaryType::Inlet, BoundaryType::Outlet);
        assert_ne!(BoundaryType::Periodic, BoundaryType::SolidWall);
    }

    // ── mirror_ghost_particle tests ───────────────────────────────────────

    #[test]
    fn test_mirror_ghost_reflects_position() {
        // Fluid at y=1, boundary at y=0, normal=[0,1,0]
        let (ghost_pos, _) = mirror_ghost_particle(
            [0.0, 1.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            false,
        );
        assert!(
            (ghost_pos[1] - (-1.0)).abs() < 1e-10,
            "ghost_y={}",
            ghost_pos[1]
        );
    }

    #[test]
    fn test_mirror_ghost_no_reflect_vel() {
        let fluid_vel = [1.0, -2.0, 3.0];
        let (_, ghost_vel) = mirror_ghost_particle(
            [0.0, 1.0, 0.0],
            fluid_vel,
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            false,
        );
        assert_eq!(ghost_vel, fluid_vel);
    }

    #[test]
    fn test_mirror_ghost_reflect_vel_removes_normal_component() {
        // Normal = [0,1,0], vel = [1,-2,3] → reflected = [1,2,3]
        let (_, ghost_vel) = mirror_ghost_particle(
            [0.0, 1.0, 0.0],
            [1.0, -2.0, 3.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            true,
        );
        assert!((ghost_vel[1] - 2.0).abs() < 1e-10, "vy={}", ghost_vel[1]);
        assert!((ghost_vel[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_mirror_ghost_symmetry() {
        // Reflecting twice should recover the original position
        let fluid_pos = [1.0, 2.0, 3.0];
        let normal = [0.0, 0.0, 1.0];
        let (ghost_pos, _) =
            mirror_ghost_particle(fluid_pos, [0.0; 3], [0.0, 0.0, 0.0], normal, false);
        let (recovered, _) =
            mirror_ghost_particle(ghost_pos, [0.0; 3], [0.0, 0.0, 0.0], normal, false);
        for d in 0..3 {
            assert!((recovered[d] - fluid_pos[d]).abs() < 1e-10);
        }
    }

    // ── solid_boundary_force tests ────────────────────────────────────────

    #[test]
    fn test_solid_boundary_force_far() {
        // Particle far from wall → zero force
        let f = solid_boundary_force([10.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1.0);
        assert!(f.iter().all(|&x| x.abs() < 1e-10));
    }

    #[test]
    fn test_solid_boundary_force_direction() {
        // Force should point away from wall (positive x when fluid is at +x)
        let f = solid_boundary_force([0.05, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1.0);
        assert!(
            f[0] > 0.0,
            "force should be repulsive (positive x), got fx={}",
            f[0]
        );
    }

    #[test]
    fn test_solid_boundary_force_zero_distance() {
        let f = solid_boundary_force([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1.0);
        assert!(f.iter().all(|&x| x.abs() < 1e-10));
    }

    #[test]
    fn test_solid_boundary_force_increases_closer() {
        let f_far = solid_boundary_force([0.08, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1.0);
        let f_near = solid_boundary_force([0.04, 0.0, 0.0], [0.0, 0.0, 0.0], 0.1, 1.0);
        assert!(
            f_near[0] > f_far[0],
            "near force {} should be > far force {}",
            f_near[0],
            f_far[0]
        );
    }

    // ── inflow_boundary tests ─────────────────────────────────────────────

    #[test]
    fn test_inflow_boundary_count() {
        let particles = inflow_boundary(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.1,
            0.05,
            5,
        );
        assert_eq!(particles.len(), 5);
    }

    #[test]
    fn test_inflow_boundary_correct_velocity() {
        let vel = [2.0, 0.5, -1.0];
        let particles = inflow_boundary([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], vel, 0.1, 0.1, 3);
        for (_, pv, _) in &particles {
            assert_eq!(*pv, vel);
        }
    }

    #[test]
    fn test_inflow_boundary_mass() {
        let particles = inflow_boundary(
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            0.25,
            0.1,
            4,
        );
        for (_, _, m) in &particles {
            assert!((*m - 0.25).abs() < 1e-10);
        }
    }

    #[test]
    fn test_inflow_boundary_zero_particles() {
        let particles = inflow_boundary(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            0.1,
            0.1,
            0,
        );
        assert!(particles.is_empty());
    }

    // ── outflow_boundary tests ────────────────────────────────────────────

    #[test]
    fn test_outflow_none_removed() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let removed = outflow_boundary(&positions, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(removed.is_empty());
    }

    #[test]
    fn test_outflow_some_removed() {
        let positions = vec![[0.0_f64, 0.0, 0.0], [2.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let removed = outflow_boundary(&positions, [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(removed.contains(&1));
        assert!(!removed.contains(&2));
    }

    #[test]
    fn test_outflow_all_removed() {
        let positions = vec![[5.0_f64, 0.0, 0.0], [6.0, 0.0, 0.0]];
        let removed = outflow_boundary(&positions, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(removed.len(), 2);
    }

    // ── periodic_boundary tests ───────────────────────────────────────────

    #[test]
    fn test_periodic_boundary_no_ghosts_needed() {
        let positions = vec![[0.5_f64, 0.0, 0.0]];
        let velocities = vec![[0.0_f64; 3]];
        let masses = vec![1.0_f64];
        let ghosts = periodic_boundary(&positions, &velocities, &masses, 0, 1.0, 0.1);
        assert!(ghosts.is_empty());
    }

    #[test]
    fn test_periodic_boundary_lower_face_ghost() {
        let positions = vec![[0.05_f64, 0.0, 0.0]];
        let velocities = vec![[1.0_f64, 0.0, 0.0]];
        let masses = vec![1.0_f64];
        let ghosts = periodic_boundary(&positions, &velocities, &masses, 0, 1.0, 0.1);
        assert_eq!(ghosts.len(), 1);
        assert!((ghosts[0].0[0] - 1.05).abs() < 1e-10);
    }

    #[test]
    fn test_periodic_boundary_upper_face_ghost() {
        let positions = vec![[0.95_f64, 0.0, 0.0]];
        let velocities = vec![[0.0_f64; 3]];
        let masses = vec![2.0_f64];
        let ghosts = periodic_boundary(&positions, &velocities, &masses, 0, 1.0, 0.1);
        assert_eq!(ghosts.len(), 1);
        assert!((ghosts[0].0[0] - (-0.05)).abs() < 1e-10);
    }

    #[test]
    fn test_periodic_boundary_velocity_preserved() {
        let positions = vec![[0.02_f64, 0.0, 0.0]];
        let vel = [3.0_f64, 2.0, 1.0];
        let velocities = vec![vel];
        let masses = vec![1.0_f64];
        let ghosts = periodic_boundary(&positions, &velocities, &masses, 0, 1.0, 0.1);
        assert_eq!(ghosts[0].1, vel);
    }

    #[test]
    fn test_periodic_boundary_y_axis() {
        let positions = vec![[0.5_f64, 0.05, 0.5]];
        let velocities = vec![[0.0_f64; 3]];
        let masses = vec![1.0_f64];
        let ghosts = periodic_boundary(&positions, &velocities, &masses, 1, 1.0, 0.1);
        assert_eq!(ghosts.len(), 1);
        assert!((ghosts[0].0[1] - 1.05).abs() < 1e-10);
    }

    // ── no_penetration_correction tests ──────────────────────────────────

    #[test]
    fn test_no_penetration_removes_inward_component() {
        // vel = [0, -1, 0] moving into wall with normal [0, 1, 0]
        let corrected = no_penetration_correction([0.0, -1.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        assert!(
            corrected[1].abs() < 1e-10,
            "vy should be ~0, got {}",
            corrected[1]
        );
    }

    #[test]
    fn test_no_penetration_outward_unchanged() {
        // vel = [0, 1, 0] moving away from wall → unchanged
        let corrected = no_penetration_correction([0.0, 1.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        assert!((corrected[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_no_penetration_free_slip() {
        // vel = [2, -1, 0], normal = [0, 1, 0], slip_factor = 1 → tangential preserved
        let corrected = no_penetration_correction([2.0, -1.0, 0.0], [0.0, 1.0, 0.0], 1.0);
        assert!(
            (corrected[0] - 2.0).abs() < 1e-10,
            "vx should be 2, got {}",
            corrected[0]
        );
    }

    #[test]
    fn test_no_penetration_no_slip() {
        // slip_factor = 0 → only normal component remains (zero here after correction)
        let corrected = no_penetration_correction([2.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.0);
        assert!(
            corrected[0].abs() < 1e-10,
            "vx should be 0 (no-slip), got {}",
            corrected[0]
        );
    }

    #[test]
    fn test_no_penetration_partial_slip() {
        let corrected = no_penetration_correction([1.0, -1.0, 0.0], [0.0, 1.0, 0.0], 0.5);
        // tangential component = [1, 0, 0], scaled by 0.5 → [0.5, 0, 0]
        assert!((corrected[0] - 0.5).abs() < 1e-10);
    }

    // ── signed_distance_to_plane tests ───────────────────────────────────

    #[test]
    fn test_signed_distance_on_plane() {
        let d = signed_distance_to_plane([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_signed_distance_positive() {
        let d = signed_distance_to_plane([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_signed_distance_negative() {
        let d = signed_distance_to_plane([-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert!((d + 1.0).abs() < 1e-10);
    }

    // ── nearest_boundary tests ────────────────────────────────────────────

    #[test]
    fn test_nearest_boundary_empty() {
        let result = nearest_boundary([0.0, 0.0, 0.0], &[]);
        assert!(result.is_none());
    }

    #[test]
    fn test_nearest_boundary_single() {
        let bp = BoundaryParticle::new(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            BoundaryType::SolidWall,
            1.0,
            0.1,
            [0.0; 3],
        );
        let idx = nearest_boundary([0.0, 0.0, 0.0], &[bp]).unwrap();
        assert_eq!(idx, 0);
    }

    #[test]
    fn test_nearest_boundary_picks_closest() {
        let bp0 = BoundaryParticle::new(
            [5.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            BoundaryType::SolidWall,
            1.0,
            0.1,
            [0.0; 3],
        );
        let bp1 = BoundaryParticle::new(
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            BoundaryType::SolidWall,
            1.0,
            0.1,
            [0.0; 3],
        );
        let idx = nearest_boundary([0.0, 0.0, 0.0], &[bp0, bp1]).unwrap();
        assert_eq!(idx, 1);
    }

    // ── Math helper tests ─────────────────────────────────────────────────

    #[test]
    fn test_normalize3_unit_vector() {
        let v = normalize3([3.0, 4.0, 0.0]);
        assert!((v[0] - 0.6).abs() < 1e-10);
        assert!((v[1] - 0.8).abs() < 1e-10);
    }

    #[test]
    fn test_normalize3_zero_vector() {
        let v = normalize3([0.0, 0.0, 0.0]);
        assert_eq!(v, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_dot3() {
        assert!((dot3([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) - 32.0).abs() < 1e-10);
    }

    #[test]
    fn test_dist3_known() {
        assert!((dist3([0.0, 0.0, 0.0], [3.0, 4.0, 0.0]) - 5.0).abs() < 1e-10);
    }
}
