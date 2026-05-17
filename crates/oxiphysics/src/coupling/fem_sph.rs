// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Immersed-boundary FEM↔SPH coupler.
//!
//! `FemSphCoupler` implements the penalty coupling method: at each
//! interface site the force on the SPH side is
//!
//! ```text
//! f = k * (v_fem - v_sph) * dt + c * (x_fem - x_sph)
//! ```
//!
//! where `k` is stiffness (Pa/m) and `c` is damping (N·s/m).
//! Newton's third law is enforced by applying `−f` to the FEM side.

use crate::coupling::{
    CouplingDomain, CouplingReport, DomainCoupler, DomainKind, InterfaceForce, InterfaceForceVec,
    InterfaceSite, InterfaceState, InterfaceStateVec,
};

// ─────────────────────────────────────────────────────────────────────────────
// FemSphCoupler
// ─────────────────────────────────────────────────────────────────────────────

/// Penalty-based coupler linking a FEM domain (A) to an SPH domain (B).
///
/// The coupler owns the list of [`InterfaceSite`]s that define where
/// sampling and force application occur.  On each [`DomainCoupler::step`]
/// call it:
///
/// 1. Pairs corresponding A and B states by position in the site list.
/// 2. Computes a penalty force per site: `f = k*(v_a - v_b)*dt + c*(x_a - x_b)`.
/// 3. Returns `−f` for domain A and `+f` for domain B.
pub struct FemSphCoupler {
    /// Interface sites, shared between both coupled domains.
    pub interface_sites: Vec<InterfaceSite>,
    /// Penalty stiffness coefficient (N/m).
    pub stiffness: f64,
    /// Penalty damping coefficient (N·s/m).
    pub damping: f64,
}

impl FemSphCoupler {
    /// Create a new coupler with the given sites and penalty coefficients.
    pub fn new(interface_sites: Vec<InterfaceSite>, stiffness: f64, damping: f64) -> Self {
        Self {
            interface_sites,
            stiffness,
            damping,
        }
    }
}

impl DomainCoupler for FemSphCoupler {
    /// Compute penalty forces and return step diagnostics.
    ///
    /// `state_a` — FEM interface states (length must equal `interface_sites`).
    /// `state_b` — SPH interface states (length must equal `interface_sites`).
    ///
    /// Returns `(forces_for_fem, forces_for_sph, report)`.
    fn step(
        &mut self,
        dt: f64,
        state_a: &InterfaceStateVec,
        state_b: &InterfaceStateVec,
    ) -> (InterfaceForceVec, InterfaceForceVec, CouplingReport) {
        let n = self
            .interface_sites
            .len()
            .min(state_a.states.len())
            .min(state_b.states.len());

        let mut forces_a = Vec::with_capacity(n);
        let mut forces_b = Vec::with_capacity(n);
        let mut momentum_a_to_b = [0.0_f64; 3];
        let mut momentum_b_to_a = [0.0_f64; 3];
        let mut sq_pos_err_sum = 0.0_f64;

        for i in 0..n {
            let sa: &InterfaceState = &state_a.states[i];
            let sb: &InterfaceState = &state_b.states[i];

            // Penalty coupling: f = k*(v_fem - v_sph)*dt + c*(x_fem - x_sph)
            let f: [f64; 3] = std::array::from_fn(|j| {
                self.stiffness * (sa.velocity[j] - sb.velocity[j]) * dt
                    + self.damping * (sa.position[j] - sb.position[j])
            });

            // Newton's third law: FEM gets −f, SPH gets +f
            forces_a.push(InterfaceForce {
                id: sa.id,
                force: [-f[0], -f[1], -f[2]],
            });
            forces_b.push(InterfaceForce {
                id: sb.id,
                force: f,
            });

            // Accumulate diagnostics
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

// ─────────────────────────────────────────────────────────────────────────────
// ConstantVelocityDomain  (mock domain for tests)
// ─────────────────────────────────────────────────────────────────────────────

/// A mock [`CouplingDomain`] that moves all interface sites at a constant
/// velocity and ignores applied forces.
///
/// Used in tests to supply deterministic interface states without running a
/// real physics solver.
pub struct ConstantVelocityDomain {
    /// The constant velocity vector in world coordinates.
    pub vel: [f64; 3],
    /// The domain kind reported by this mock.
    pub kind: DomainKind,
}

impl ConstantVelocityDomain {
    /// Create a new constant-velocity mock domain.
    pub fn new(vel: [f64; 3], kind: DomainKind) -> Self {
        Self { vel, kind }
    }
}

impl CouplingDomain for ConstantVelocityDomain {
    fn kind(&self) -> DomainKind {
        self.kind.clone()
    }

    /// Return one [`InterfaceState`] per input site, all with the same
    /// constant velocity and the site's own position.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec {
        let states = sites
            .iter()
            .map(|s| InterfaceState {
                id: s.id,
                position: s.position,
                velocity: self.vel,
            })
            .collect();
        InterfaceStateVec { states }
    }

    /// This mock domain ignores applied forces (kinematically driven).
    fn apply_interface_force(&mut self, _forces: &InterfaceForceVec) {
        // Kinematic domain — forces have no effect.
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImpulseDomain  (mock domain with mass for integration tests)
// ─────────────────────────────────────────────────────────────────────────────

/// A mock [`CouplingDomain`] that accumulates velocity impulses from applied
/// forces.
///
/// Useful for verifying that Newton's third law is satisfied: the change in
/// momentum of this domain should equal the impulse applied by the coupler.
pub struct ImpulseDomain {
    /// Current velocity of each interface node.
    pub velocities: Vec<[f64; 3]>,
    /// Per-node mass (kg).
    pub masses: Vec<f64>,
    /// The domain kind reported by this mock.
    pub kind: DomainKind,
}

impl ImpulseDomain {
    /// Create a new impulse domain with uniform mass.
    pub fn new(n: usize, mass: f64, kind: DomainKind) -> Self {
        Self {
            velocities: vec![[0.0; 3]; n],
            masses: vec![mass; n],
            kind,
        }
    }
}

impl CouplingDomain for ImpulseDomain {
    fn kind(&self) -> DomainKind {
        self.kind.clone()
    }

    /// Return one [`InterfaceState`] per input site, using stored velocities.
    ///
    /// If a site index (cast from `id`) exceeds the stored count, the
    /// returned state has zero velocity.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec {
        let states = sites
            .iter()
            .map(|s| {
                let idx = s.id as usize;
                let velocity = if idx < self.velocities.len() {
                    self.velocities[idx]
                } else {
                    [0.0; 3]
                };
                InterfaceState {
                    id: s.id,
                    position: s.position,
                    velocity,
                }
            })
            .collect();
        InterfaceStateVec { states }
    }

    /// Apply impulse `f/mass` to the node identified by `force.id`.
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
