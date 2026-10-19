// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adapter wrapping MD / continuum regions as `CouplingDomain` objects.
//!
//! `MdContinuumAdapter` is a `DomainCoupler` that links a particle-based
//! MD region (`MockMdDomain`) to a mesh-based continuum region
//! (`MockContinuumDomain`) using a penalty force analogous to the
//! Bridging Domain Method (BDM): a spring-damper that pulls the two sides of
//! every interface site together and damps their relative velocity.

use crate::coupling::{
    CouplingDomain, CouplingReport, DomainCoupler, DomainKind, InterfaceForce, InterfaceForceVec,
    InterfaceSite, InterfaceState, InterfaceStateVec,
};

// ─────────────────────────────────────────────────────────────────────────────
// MockMdDomain
// ─────────────────────────────────────────────────────────────────────────────

/// Mock molecular-dynamics domain for coupling tests.
///
/// Stores particle positions, velocities, and masses and integrates the
/// applied forces with semi-implicit Euler (`v += F * dt / m`, then
/// `x += v * dt`).  Only particles that receive a force entry are advanced, and
/// they are advanced once per received entry, so list every interface particle
/// once per step.
pub struct MockMdDomain {
    /// Particle positions in world coordinates.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses (kg).
    pub masses: Vec<f64>,
}

impl MockMdDomain {
    /// Construct a domain with `n` stationary particles at the origin, each
    /// with the given `mass`.
    pub fn new(n: usize, mass: f64) -> Self {
        Self {
            positions: vec![[0.0; 3]; n],
            velocities: vec![[0.0; 3]; n],
            masses: vec![mass; n],
        }
    }

    /// Construct a domain from explicit particle data.
    ///
    /// # Panics (debug)
    ///
    /// Panics in debug builds if the three slices differ in length.
    pub fn from_particles(
        positions: Vec<[f64; 3]>,
        velocities: Vec<[f64; 3]>,
        masses: Vec<f64>,
    ) -> Self {
        debug_assert_eq!(
            positions.len(),
            velocities.len(),
            "positions and velocities must have the same length"
        );
        debug_assert_eq!(
            positions.len(),
            masses.len(),
            "positions and masses must have the same length"
        );
        Self {
            positions,
            velocities,
            masses,
        }
    }
}

impl CouplingDomain for MockMdDomain {
    fn kind(&self) -> DomainKind {
        DomainKind::Md
    }

    /// Sample states at the requested sites.
    ///
    /// The site `id` is used as a particle index.  If the index is out of
    /// bounds, the site's own position is returned with zero velocity.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec {
        let states = sites
            .iter()
            .map(|s| {
                let particle = InterfaceSite::local_index(s.id)
                    .and_then(|idx| Some((self.positions.get(idx)?, self.velocities.get(idx)?)));
                let (position, velocity) = match particle {
                    Some((&position, &velocity)) => (position, velocity),
                    None => (s.position, [0.0; 3]),
                };
                InterfaceState {
                    id: s.id,
                    position,
                    velocity,
                }
            })
            .collect();
        InterfaceStateVec { states }
    }

    /// Integrate the force on each particle identified by `force.id`:
    /// `v += F * dt / m`, then `x += v * dt`.  Massless particles (`m <= ε`)
    /// and unknown ids are left untouched.
    fn apply_interface_force(&mut self, forces: &InterfaceForceVec, dt: f64) {
        for iforce in &forces.forces {
            let Some(idx) = InterfaceSite::local_index(iforce.id) else {
                continue;
            };
            let (Some(&mass), Some(velocity), Some(position)) = (
                self.masses.get(idx),
                self.velocities.get_mut(idx),
                self.positions.get_mut(idx),
            ) else {
                continue;
            };
            if mass > f64::EPSILON {
                for j in 0..3 {
                    velocity[j] += iforce.force[j] / mass * dt;
                    position[j] += velocity[j] * dt;
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MockContinuumDomain
// ─────────────────────────────────────────────────────────────────────────────

/// Mock continuum (FEM/FVM) domain for coupling tests.
///
/// Stores node positions, velocities, and masses and integrates the applied
/// forces with semi-implicit Euler (`v += F * dt / m`, then `x += v * dt`).
/// Only nodes that receive a force entry are advanced, and they are advanced
/// once per received entry, so list every interface node once per step.
pub struct MockContinuumDomain {
    /// Node positions in world coordinates.
    pub node_positions: Vec<[f64; 3]>,
    /// Node velocities.
    pub node_velocities: Vec<[f64; 3]>,
    /// Nodal masses (kg).
    pub node_masses: Vec<f64>,
}

impl MockContinuumDomain {
    /// Construct a domain with `n` stationary nodes at the origin, each with
    /// the given `mass`.
    pub fn new(n: usize, mass: f64) -> Self {
        Self {
            node_positions: vec![[0.0; 3]; n],
            node_velocities: vec![[0.0; 3]; n],
            node_masses: vec![mass; n],
        }
    }

    /// Construct a domain from explicit node data.
    ///
    /// # Panics (debug)
    ///
    /// Panics in debug builds if the three slices differ in length.
    pub fn from_nodes(
        node_positions: Vec<[f64; 3]>,
        node_velocities: Vec<[f64; 3]>,
        node_masses: Vec<f64>,
    ) -> Self {
        debug_assert_eq!(
            node_positions.len(),
            node_velocities.len(),
            "positions and velocities must have the same length"
        );
        debug_assert_eq!(
            node_positions.len(),
            node_masses.len(),
            "positions and masses must have the same length"
        );
        Self {
            node_positions,
            node_velocities,
            node_masses,
        }
    }
}

impl CouplingDomain for MockContinuumDomain {
    fn kind(&self) -> DomainKind {
        DomainKind::Fem
    }

    /// Sample states at the requested sites.
    ///
    /// The site `id` is used as a node index.  If the index is out of bounds,
    /// the site's own position is returned with zero velocity.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec {
        let states = sites
            .iter()
            .map(|s| {
                let node = InterfaceSite::local_index(s.id).and_then(|idx| {
                    Some((
                        self.node_positions.get(idx)?,
                        self.node_velocities.get(idx)?,
                    ))
                });
                let (position, velocity) = match node {
                    Some((&position, &velocity)) => (position, velocity),
                    None => (s.position, [0.0; 3]),
                };
                InterfaceState {
                    id: s.id,
                    position,
                    velocity,
                }
            })
            .collect();
        InterfaceStateVec { states }
    }

    /// Integrate the force on each node identified by `force.id`:
    /// `v += F * dt / m`, then `x += v * dt`.  Massless nodes (`m <= ε`) and
    /// unknown ids are left untouched.
    fn apply_interface_force(&mut self, forces: &InterfaceForceVec, dt: f64) {
        for iforce in &forces.forces {
            let Some(idx) = InterfaceSite::local_index(iforce.id) else {
                continue;
            };
            let (Some(&mass), Some(velocity), Some(position)) = (
                self.node_masses.get(idx),
                self.node_velocities.get_mut(idx),
                self.node_positions.get_mut(idx),
            ) else {
                continue;
            };
            if mass > f64::EPSILON {
                for j in 0..3 {
                    velocity[j] += iforce.force[j] / mass * dt;
                    position[j] += velocity[j] * dt;
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// MdContinuumAdapter
// ─────────────────────────────────────────────────────────────────────────────

/// Penalty coupler bridging an MD domain (A) and a continuum domain (B).
///
/// Implements a penalty form of the Bridging Domain Method: at each interface
/// site a spring-damper joins the MD particle and the continuum node, with the
/// force (newtons)
///
/// ```text
/// f = k * (x_md - x_cont) + c * (v_md - v_cont)
/// ```
///
/// `−f` is applied to the MD side and `+f` to the continuum side (Newton's
/// third law), so the MD particle is pulled towards the continuum node and
/// their relative velocity is damped.
pub struct MdContinuumAdapter {
    /// Interface sites shared by the MD and the continuum region.  Each
    /// domain reads [`InterfaceSite::id`] as its own local particle / node
    /// index.
    pub sites: Vec<InterfaceSite>,
    /// Penalty stiffness `k` (N/m), acting on the position mismatch.
    pub stiffness: f64,
    /// Penalty damping `c` (N·s/m), acting on the relative velocity.
    pub damping: f64,
}

impl MdContinuumAdapter {
    /// Create a new adapter for the given interface `sites` and penalty
    /// coefficients.
    pub fn new(sites: Vec<InterfaceSite>, stiffness: f64, damping: f64) -> Self {
        Self {
            sites,
            stiffness,
            damping,
        }
    }
}

impl DomainCoupler for MdContinuumAdapter {
    fn interface_sites(&self) -> &[InterfaceSite] {
        &self.sites
    }

    /// Compute BDM penalty forces between an MD region (A) and a continuum
    /// region (B).
    fn step(
        &mut self,
        dt: f64,
        state_a: &InterfaceStateVec,
        state_b: &InterfaceStateVec,
    ) -> (InterfaceForceVec, InterfaceForceVec, CouplingReport) {
        let n = state_a.states.len().min(state_b.states.len());

        let mut forces_a = Vec::with_capacity(n);
        let mut forces_b = Vec::with_capacity(n);
        let mut momentum_a_to_b = [0.0_f64; 3];
        let mut momentum_b_to_a = [0.0_f64; 3];
        let mut sq_pos_err_sum = 0.0_f64;

        for i in 0..n {
            let sa = &state_a.states[i];
            let sb = &state_b.states[i];

            // Spring-damper (newtons): f = k*(x_md - x_cont) + c*(v_md - v_cont).
            let f: [f64; 3] = std::array::from_fn(|j| {
                self.stiffness * (sa.position[j] - sb.position[j])
                    + self.damping * (sa.velocity[j] - sb.velocity[j])
            });

            // MD gets −f (pulled towards the continuum node), the continuum
            // gets +f (pulled towards the MD particle).
            forces_a.push(InterfaceForce {
                id: sa.id,
                force: [-f[0], -f[1], -f[2]],
            });
            forces_b.push(InterfaceForce {
                id: sb.id,
                force: f,
            });

            for j in 0..3 {
                momentum_a_to_b[j] += f[j] * dt;
                momentum_b_to_a[j] -= f[j] * dt;
                let dx = sa.position[j] - sb.position[j];
                sq_pos_err_sum += dx * dx;
            }
        }

        let rms_position_error = if n > 0 {
            (sq_pos_err_sum / (3 * n) as f64).sqrt()
        } else {
            0.0
        };

        let report = CouplingReport {
            site_count: n,
            momentum_a_to_b,
            momentum_b_to_a,
            rms_position_error,
        };

        (
            InterfaceForceVec { forces: forces_a },
            InterfaceForceVec { forces: forces_b },
            report,
        )
    }
}
