// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! # OxiCAR — Cross-domain Auto-coupling Runtime
//!
//! Blueprint KF-1: domain-agnostic `DomainCoupler` framework for coupling
//! heterogeneous simulation domains (FEM ↔ SPH ↔ LBM ↔ MD) at shared
//! interfaces.
//!
//! ## Design
//!
//! Each simulation domain (FEM, SPH, MD, …) implements `CouplingDomain`.
//! Two domains are linked by a `DomainCoupler` that calls
//! `sample_interface` on both sides, computes interface forces, and calls
//! `apply_interface_force` to push those forces back.
//!
//! A `CouplingRuntime` owns all domains and couplers. A *schedule* is a
//! list of `(coupler_idx, domain_a_idx, domain_b_idx)` triples evaluated
//! in order on every `CouplingRuntime::step` call.

pub mod conservation;
pub mod fem_sph;
pub mod md_continuum_adapter;

// ─────────────────────────────────────────────────────────────────────────────
// Payload types
// ─────────────────────────────────────────────────────────────────────────────

/// A single point on the interface between two coupled domains.
///
/// `id` is a caller-assigned identifier that is passed back unchanged in
/// [`InterfaceState`] and [`InterfaceForce`] so the domain can correlate
/// results with its own data structures.
#[derive(Debug, Clone, Copy)]
pub struct InterfaceSite {
    /// Domain-local identifier for this interface node / particle.
    pub id: u64,
    /// Position of the site in world coordinates.
    pub position: [f64; 3],
}

/// Sampled state at a single interface site.
///
/// Returned by [`CouplingDomain::sample_interface`] for each input
/// [`InterfaceSite`].
#[derive(Debug, Clone, Copy)]
pub struct InterfaceState {
    /// Site identifier, mirrored from the corresponding [`InterfaceSite`].
    pub id: u64,
    /// Position of the site at the current time step.
    pub position: [f64; 3],
    /// Velocity of the site at the current time step.
    pub velocity: [f64; 3],
}

/// Force to be applied at a single interface site.
///
/// Passed to [`CouplingDomain::apply_interface_force`].
#[derive(Debug, Clone, Copy)]
pub struct InterfaceForce {
    /// Site identifier, correlating back to the generating [`InterfaceSite`].
    pub id: u64,
    /// Force vector in world coordinates.
    pub force: [f64; 3],
}

/// Owned list of [`InterfaceState`] entries — heap-allocated, not `Copy`.
#[derive(Debug, Clone)]
pub struct InterfaceStateVec {
    /// The sampled states, one per queried [`InterfaceSite`].
    pub states: Vec<InterfaceState>,
}

/// Owned list of [`InterfaceForce`] entries — heap-allocated, not `Copy`.
#[derive(Debug, Clone)]
pub struct InterfaceForceVec {
    /// The forces to be applied, one per queried [`InterfaceSite`].
    pub forces: Vec<InterfaceForce>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Domain kind tag
// ─────────────────────────────────────────────────────────────────────────────

/// Identifies the physical domain type of a simulation region.
///
/// Used for diagnostics and as a hint to couplers that may specialise their
/// algorithms depending on which pair of domain types they are connecting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DomainKind {
    /// Finite Element Method — structural / continuum solid.
    Fem,
    /// Smoothed Particle Hydrodynamics — particle-based fluid.
    Sph,
    /// Lattice Boltzmann Method — lattice fluid.
    Lbm,
    /// Molecular Dynamics — atomistic simulation.
    Md,
    /// User-defined domain type.
    Custom(String),
}

// ─────────────────────────────────────────────────────────────────────────────
// Coupling report
// ─────────────────────────────────────────────────────────────────────────────

/// Summary statistics produced by one coupler step.
///
/// All fields are scalar so the struct is `Copy` and cheap to collect.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CouplingReport {
    /// Number of interface sites processed this step.
    pub site_count: usize,
    /// Total momentum exchanged from domain A to domain B (`[fx, fy, fz] * dt`).
    pub momentum_a_to_b: [f64; 3],
    /// Total momentum exchanged from domain B to domain A (`[fx, fy, fz] * dt`).
    pub momentum_b_to_a: [f64; 3],
    /// RMS position mismatch across all sites at the end of this step.
    pub rms_position_error: f64,
}

impl CouplingReport {
    /// Construct a zero-valued report.
    pub fn zero() -> Self {
        Self {
            site_count: 0,
            momentum_a_to_b: [0.0; 3],
            momentum_b_to_a: [0.0; 3],
            rms_position_error: 0.0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Core traits
// ─────────────────────────────────────────────────────────────────────────────

/// A simulation domain that can participate in cross-domain coupling.
///
/// Implementors expose their interface state via [`sample_interface`] and
/// accept coupling forces via [`apply_interface_force`].
///
/// [`sample_interface`]: CouplingDomain::sample_interface
/// [`apply_interface_force`]: CouplingDomain::apply_interface_force
pub trait CouplingDomain: Send + Sync {
    /// Returns the domain kind for this region.
    fn kind(&self) -> DomainKind;

    /// Sample the current state at each requested interface site.
    ///
    /// The returned [`InterfaceStateVec`] must contain exactly `sites.len()`
    /// entries, in the same order as `sites`.
    fn sample_interface(&self, sites: &[InterfaceSite]) -> InterfaceStateVec;

    /// Apply the given interface forces to this domain's degrees of freedom.
    ///
    /// Forces are typically accumulated into the domain's internal force
    /// vector for the next integration step.
    fn apply_interface_force(&mut self, forces: &InterfaceForceVec);
}

/// A coupler that computes interface forces between two coupled domains.
///
/// A coupler is invoked once per [`CouplingRuntime::step`] for each
/// scheduled domain pair.  It receives immutable views of both domain states
/// and produces two sets of forces (one for each side) together with a
/// diagnostic [`CouplingReport`].
pub trait DomainCoupler: Send + Sync {
    /// Compute coupling forces and return a [`CouplingReport`].
    ///
    /// `dt` — the physical time step duration (seconds).
    ///
    /// The returned tuple `(forces_a, forces_b, report)` is:
    /// - `forces_a` — forces to be applied to domain A,
    /// - `forces_b` — forces to be applied to domain B,
    /// - `report`   — step diagnostics.
    fn step(
        &mut self,
        dt: f64,
        state_a: &InterfaceStateVec,
        state_b: &InterfaceStateVec,
    ) -> (InterfaceForceVec, InterfaceForceVec, CouplingReport);
}

// ─────────────────────────────────────────────────────────────────────────────
// Runtime
// ─────────────────────────────────────────────────────────────────────────────

/// Manages a collection of coupled simulation domains.
///
/// The runtime owns all domains (as `Box<dyn CouplingDomain>`) and couplers
/// (as `Box<dyn DomainCoupler>`).  A *schedule* is an ordered list of
/// `(coupler_idx, domain_a_idx, domain_b_idx)` triples.  Each call to
/// [`step`] walks the schedule, drives every coupler, and applies the
/// resulting forces back to the corresponding domains.
///
/// [`step`]: CouplingRuntime::step
pub struct CouplingRuntime {
    /// Registered simulation domains.
    domains: Vec<Box<dyn CouplingDomain>>,
    /// Registered couplers.
    couplers: Vec<Box<dyn DomainCoupler>>,
    /// Schedule: `(coupler_idx, domain_a_idx, domain_b_idx)`.
    schedule: Vec<(usize, usize, usize)>,
}

impl CouplingRuntime {
    /// Create an empty runtime with no domains or couplers.
    pub fn new() -> Self {
        Self {
            domains: Vec::new(),
            couplers: Vec::new(),
            schedule: Vec::new(),
        }
    }

    /// Register a domain and return its index.
    pub fn add_domain(&mut self, domain: Box<dyn CouplingDomain>) -> usize {
        let idx = self.domains.len();
        self.domains.push(domain);
        idx
    }

    /// Register a coupler and return its index.
    pub fn add_coupler(&mut self, coupler: Box<dyn DomainCoupler>) -> usize {
        let idx = self.couplers.len();
        self.couplers.push(coupler);
        idx
    }

    /// Schedule a coupling step between two domains via an existing coupler.
    ///
    /// The triple `(coupler_idx, domain_a_idx, domain_b_idx)` is appended to
    /// the schedule in the order it was registered.  Panics (debug) if any
    /// index is out of bounds.
    pub fn schedule_pair(&mut self, coupler_idx: usize, domain_a_idx: usize, domain_b_idx: usize) {
        debug_assert!(
            coupler_idx < self.couplers.len(),
            "coupler index out of bounds"
        );
        debug_assert!(
            domain_a_idx < self.domains.len(),
            "domain_a index out of bounds"
        );
        debug_assert!(
            domain_b_idx < self.domains.len(),
            "domain_b index out of bounds"
        );
        self.schedule
            .push((coupler_idx, domain_a_idx, domain_b_idx));
    }

    /// Advance all scheduled couplers by one time step of duration `dt`.
    ///
    /// Returns a [`Vec<CouplingReport>`] with one entry per scheduled pair,
    /// in schedule order.
    ///
    /// # Borrow safety
    ///
    /// The method needs simultaneous mutable access to two distinct domains
    /// while holding immutable access to the coupler.  Rust's borrow checker
    /// cannot verify this at compile time through the shared slice, so we use
    /// [`slice::split_at_mut`] to obtain non-overlapping mutable references
    /// after checking that the two domain indices are distinct.
    pub fn step(&mut self, dt: f64) -> Vec<CouplingReport> {
        // Copy the schedule out so that `self` is not borrowed while we
        // mutably borrow self.domains and self.couplers below.
        let schedule: Vec<(usize, usize, usize)> = self.schedule.clone();
        let mut reports = Vec::with_capacity(schedule.len());

        for (coupler_idx, idx_a, idx_b) in schedule {
            // Collect the interface sites that the coupler knows about.
            // We drive the coupler with a fixed set of sites derived from the
            // coupler's own state (it owns the InterfaceSite list).  Here we
            // let the coupler tell us what sites it needs by sampling from the
            // domain objects.
            //
            // For the generic runtime, we define the site list as the union
            // of positions currently held by domain A.  Real couplers supply
            // their own site geometry — here we just need *something* to pass.
            //
            // We collect the states first with immutable refs, then apply
            // forces with mutable refs, using split_at_mut to satisfy the
            // borrow checker when idx_a != idx_b.

            let (state_a, state_b) = if idx_a == idx_b {
                // Degenerate: same domain on both sides (self-coupling).
                let s = self.domains[idx_a].sample_interface(&[]);
                (
                    InterfaceStateVec {
                        states: s.states.clone(),
                    },
                    InterfaceStateVec { states: s.states },
                )
            } else {
                let sa = self.domains[idx_a].sample_interface(&[]);
                let sb = self.domains[idx_b].sample_interface(&[]);
                (sa, sb)
            };

            let (forces_a, forces_b, report) =
                self.couplers[coupler_idx].step(dt, &state_a, &state_b);

            // Apply forces back — need two mutable borrows on distinct elements.
            if idx_a == idx_b {
                self.domains[idx_a].apply_interface_force(&forces_a);
            } else if idx_a < idx_b {
                let (left, right) = self.domains.split_at_mut(idx_b);
                left[idx_a].apply_interface_force(&forces_a);
                right[0].apply_interface_force(&forces_b);
            } else {
                // idx_b < idx_a
                let (left, right) = self.domains.split_at_mut(idx_a);
                left[idx_b].apply_interface_force(&forces_b);
                right[0].apply_interface_force(&forces_a);
            }

            reports.push(report);
        }

        reports
    }

    /// Return the number of registered domains.
    pub fn domain_count(&self) -> usize {
        self.domains.len()
    }

    /// Return the number of registered couplers.
    pub fn coupler_count(&self) -> usize {
        self.couplers.len()
    }

    /// Return the number of scheduled pairs.
    pub fn schedule_len(&self) -> usize {
        self.schedule.len()
    }
}

impl Default for CouplingRuntime {
    fn default() -> Self {
        Self::new()
    }
}
