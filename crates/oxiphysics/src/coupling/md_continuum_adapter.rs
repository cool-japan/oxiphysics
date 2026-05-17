// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adapter wrapping MD / continuum regions as `CouplingDomain` objects.
//!
//! `MdContinuumAdapter` is a `DomainCoupler` that links a particle-based
//! MD region (`MockMdDomain`) to a mesh-based continuum region
//! (`MockContinuumDomain`) using a penalty force analogous to the
//! Bridging Domain Method (BDM).

use crate::coupling::{
    CouplingDomain, CouplingReport, DomainCoupler, DomainKind, InterfaceForce, InterfaceForceVec,
    InterfaceSite, InterfaceState, InterfaceStateVec,
};

// ─────────────────────────────────────────────────────────────────────────────
// MockMdDomain
// ─────────────────────────────────────────────────────────────────────────────

/// Mock molecular-dynamics domain for coupling tests.
///
/// Stores particle positions, velocities, and masses and accumulates
/// velocity impulses from applied forces.
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
                let idx = s.id as usize;
                let (position, velocity) = if idx < self.positions.len() {
                    (self.positions[idx], self.velocities[idx])
                } else {
                    (s.position, [0.0; 3])
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

    /// Apply impulse `f / mass` to each particle identified by `force.id`.
    fn apply_interface_force(&mut self, forces: &InterfaceForceVec) {
        for iforce in &forces.forces {
            let idx = iforce.id as usize;
            if idx < self.velocities.len() {
                let mass = self.masses[idx];
                if mass > f64::EPSILON {
                    for j in 0..3 {
                        self.velocities[idx][j] += iforce.force[j] / mass;
                    }
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
/// Stores node positions, velocities, and masses and accumulates velocity
/// impulses from applied forces.
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
                let idx = s.id as usize;
                let (position, velocity) = if idx < self.node_positions.len() {
                    (self.node_positions[idx], self.node_velocities[idx])
                } else {
                    (s.position, [0.0; 3])
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

    /// Apply impulse `f / mass` to each node identified by `force.id`.
    fn apply_interface_force(&mut self, forces: &InterfaceForceVec) {
        for iforce in &forces.forces {
            let idx = iforce.id as usize;
            if idx < self.node_velocities.len() {
                let mass = self.node_masses[idx];
                if mass > f64::EPSILON {
                    for j in 0..3 {
                        self.node_velocities[idx][j] += iforce.force[j] / mass;
                    }
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
/// Implements a simplified Bridging Domain Method: the penalty force at each
/// interface site is
///
/// ```text
/// f = k * (v_md - v_cont) * dt + c * (x_md - x_cont)
/// ```
///
/// with `−f` applied to the continuum side and `+f` applied to the MD side
/// (Newton's third law).
pub struct MdContinuumAdapter {
    /// Penalty stiffness coefficient (N/m).
    pub stiffness: f64,
    /// Penalty damping coefficient (N·s/m).
    pub damping: f64,
}

impl MdContinuumAdapter {
    /// Create a new adapter with the given penalty coefficients.
    pub fn new(stiffness: f64, damping: f64) -> Self {
        Self { stiffness, damping }
    }
}

impl DomainCoupler for MdContinuumAdapter {
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

            // f = k*(v_md - v_cont)*dt + c*(x_md - x_cont)
            let f: [f64; 3] = std::array::from_fn(|j| {
                self.stiffness * (sa.velocity[j] - sb.velocity[j]) * dt
                    + self.damping * (sa.position[j] - sb.position[j])
            });

            // MD gets +f, continuum gets −f
            forces_a.push(InterfaceForce {
                id: sa.id,
                force: f,
            });
            forces_b.push(InterfaceForce {
                id: sb.id,
                force: [-f[0], -f[1], -f[2]],
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
