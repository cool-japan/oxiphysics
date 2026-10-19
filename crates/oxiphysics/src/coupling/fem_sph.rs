// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Immersed-boundary FEM↔SPH coupler.
//!
//! `FemSphCoupler` implements the penalty coupling method: at each
//! interface site a spring-damper joins the FEM node and the SPH particle,
//! and the force on the SPH side is
//!
//! ```text
//! f = k * (x_fem - x_sph) + c * (v_fem - v_sph)
//! ```
//!
//! where `k` is the penalty stiffness (N/m) and `c` the penalty damping
//! (N·s/m), so `f` is a force in newtons.  The SPH side receives `+f` (it is
//! pulled towards the FEM node and dragged along with it) and, by Newton's
//! third law, the FEM side receives `−f`.

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
/// 2. Computes a spring-damper force per site:
///    `f = k*(x_a - x_b) + c*(v_a - v_b)` (newtons).
/// 3. Returns `−f` for domain A and `+f` for domain B.
pub struct FemSphCoupler {
    /// Interface sites, shared between both coupled domains.
    pub interface_sites: Vec<InterfaceSite>,
    /// Penalty stiffness `k` (N/m), acting on the position mismatch.
    pub stiffness: f64,
    /// Penalty damping `c` (N·s/m), acting on the relative velocity.
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
    fn interface_sites(&self) -> &[InterfaceSite] {
        &self.interface_sites
    }

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

            // Spring-damper penalty coupling (newtons):
            //   f = k*(x_fem - x_sph) + c*(v_fem - v_sph)
            // The stiffness multiplies the position mismatch and the damping
            // the velocity mismatch, so the force is restoring and dissipative.
            let f: [f64; 3] = std::array::from_fn(|j| {
                self.stiffness * (sa.position[j] - sb.position[j])
                    + self.damping * (sa.velocity[j] - sb.velocity[j])
            });

            // Newton's third law: FEM gets −f (pulled towards the SPH side),
            // SPH gets +f (pulled towards the FEM side).
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
/// The domain is kinematically driven: a site that starts at `p` is reported at
/// `p + vel * elapsed`, where `elapsed` is the sum of the step lengths received
/// through [`CouplingDomain::apply_interface_force`] (one call advances the
/// clock by one `dt`, so couple the domain through a single scheduled pair).
///
/// Used in tests to supply deterministic interface states without running a
/// real physics solver.
pub struct ConstantVelocityDomain {
    /// The constant velocity vector in world coordinates.
    pub vel: [f64; 3],
    /// The domain kind reported by this mock.
    pub kind: DomainKind,
    /// Time advanced so far (seconds): the sum of the `dt` of all
    /// `apply_interface_force` calls.
    pub elapsed: f64,
}

impl ConstantVelocityDomain {
    /// Create a new constant-velocity mock domain at time zero.
    pub fn new(vel: [f64; 3], kind: DomainKind) -> Self {
        Self {
            vel,
            kind,
            elapsed: 0.0,
        }
    }
}

impl CouplingDomain for ConstantVelocityDomain {
    fn kind(&self) -> DomainKind {
        self.kind.clone()
    }

    /// Return one [`InterfaceState`] per input site, all with the same
    /// constant velocity and the site's own position advanced by
    /// `vel * elapsed`.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec {
        let states = sites
            .iter()
            .map(|s| InterfaceState {
                id: s.id,
                position: std::array::from_fn(|j| s.position[j] + self.vel[j] * self.elapsed),
                velocity: self.vel,
            })
            .collect();
        InterfaceStateVec { states }
    }

    /// This mock domain ignores the applied forces (kinematically driven); it
    /// only advances its clock by `dt`.
    fn apply_interface_force(&mut self, _forces: &InterfaceForceVec, dt: f64) {
        self.elapsed += dt;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// ImpulseDomain  (mock domain with mass for integration tests)
// ─────────────────────────────────────────────────────────────────────────────

/// A mock [`CouplingDomain`] that integrates the applied forces with
/// semi-implicit Euler: `v += F * dt / m`, then `x += v * dt`.
///
/// Useful for verifying that Newton's third law is satisfied: the change in
/// momentum of this domain (`m * Δv`) equals the impulse `F * dt` applied by
/// the coupler.  Node positions are tracked as a displacement relative to the
/// interface site at which the node is sampled.
pub struct ImpulseDomain {
    /// Current velocity of each interface node.
    pub velocities: Vec<[f64; 3]>,
    /// Per-node mass (kg).
    pub masses: Vec<f64>,
    /// Accumulated displacement (m) of each node from its interface site.
    pub displacements: Vec<[f64; 3]>,
    /// The domain kind reported by this mock.
    pub kind: DomainKind,
}

impl ImpulseDomain {
    /// Create a new impulse domain with uniform mass.
    pub fn new(n: usize, mass: f64, kind: DomainKind) -> Self {
        Self {
            velocities: vec![[0.0; 3]; n],
            masses: vec![mass; n],
            displacements: vec![[0.0; 3]; n],
            kind,
        }
    }
}

impl CouplingDomain for ImpulseDomain {
    fn kind(&self) -> DomainKind {
        self.kind.clone()
    }

    /// Return one [`InterfaceState`] per input site: the stored velocity and
    /// the site position shifted by the node's accumulated displacement.
    ///
    /// If a site index (cast from `id`) exceeds the stored count, the
    /// returned state is the site's own position with zero velocity.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec {
        let states = sites
            .iter()
            .map(|s| {
                let node = InterfaceSite::local_index(s.id);
                let velocity = node
                    .and_then(|idx| self.velocities.get(idx))
                    .copied()
                    .unwrap_or([0.0; 3]);
                let displacement = node
                    .and_then(|idx| self.displacements.get(idx))
                    .copied()
                    .unwrap_or([0.0; 3]);
                InterfaceState {
                    id: s.id,
                    position: std::array::from_fn(|j| s.position[j] + displacement[j]),
                    velocity,
                }
            })
            .collect();
        InterfaceStateVec { states }
    }

    /// Integrate the force on the node identified by `force.id`:
    /// `v += F * dt / m`, then `d += v * dt`.  Massless nodes (`m <= ε`) and
    /// unknown ids are left untouched.
    fn apply_interface_force(&mut self, forces: &InterfaceForceVec, dt: f64) {
        for iforce in &forces.forces {
            let Some(idx) = InterfaceSite::local_index(iforce.id) else {
                continue;
            };
            let (Some(&mass), Some(velocity), Some(displacement)) = (
                self.masses.get(idx),
                self.velocities.get_mut(idx),
                self.displacements.get_mut(idx),
            ) else {
                continue;
            };
            if mass > f64::EPSILON {
                for j in 0..3 {
                    velocity[j] += iforce.force[j] / mass * dt;
                    displacement[j] += velocity[j] * dt;
                }
            }
        }
    }
}
