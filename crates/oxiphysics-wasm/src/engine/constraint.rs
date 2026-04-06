// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmConstraintSolver` — iterative constraint solver for the WASM engine.

#![allow(missing_docs)]

use super::WasmPhysicsEngine;

/// A single bilateral constraint between two bodies (or one body and ground).
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WasmConstraint {
    /// Handle of the first body (or `u32::MAX` for world).
    pub body_a: u32,
    /// Handle of the second body (or `u32::MAX` for world).
    pub body_b: u32,
    /// Constraint Jacobian direction (unit vector along which relative motion
    /// must be zero).
    pub normal: [f64; 3],
    /// Bias (penetration depth / dt * Baumgarte factor).
    pub bias: f64,
    /// Accumulated lambda for warm-starting.
    pub lambda: f64,
}

/// Contact data for the constraint solver.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ContactData {
    /// First body handle.
    pub body_a: u32,
    /// Second body handle.
    pub body_b: u32,
    /// Contact normal (pointing from B to A).
    pub normal: [f64; 3],
    /// Penetration depth.
    pub depth: f64,
    /// Combined restitution.
    pub restitution: f64,
}

/// A simple iterative (Gauss-Seidel) constraint solver.
///
/// Works in conjunction with a `WasmPhysicsEngine` by accepting contact data,
/// building constraint rows, and iterating to convergence.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct WasmConstraintSolver {
    constraints: Vec<WasmConstraint>,
    contacts: Vec<ContactData>,
    /// Number of Gauss-Seidel iterations per solve call.
    pub(crate) iterations: u32,
    /// Baumgarte stabilisation factor (typical: 0.1–0.3).
    baumgarte: f64,
    /// Allowed penetration slop before correction (metres).
    slop: f64,
    /// Velocity corrections per body: keyed by body handle.
    velocity_corrections: std::collections::HashMap<u32, [f64; 3]>,
}

impl WasmConstraintSolver {
    /// Create a new constraint solver.
    pub fn new() -> Self {
        Self {
            constraints: Vec::new(),
            contacts: Vec::new(),
            iterations: 10,
            baumgarte: 0.2,
            slop: 0.001,
            velocity_corrections: std::collections::HashMap::new(),
        }
    }

    /// Set the number of solver iterations.
    pub fn set_iterations(&mut self, n: u32) {
        self.iterations = n.max(1);
    }

    /// Set the Baumgarte stabilisation factor.
    pub fn set_baumgarte(&mut self, b: f64) {
        self.baumgarte = b.clamp(0.0, 1.0);
    }

    /// Add a contact for the next solve call.
    #[allow(clippy::too_many_arguments)]
    pub fn add_contact(
        &mut self,
        body_a: u32,
        body_b: u32,
        nx: f64,
        ny: f64,
        nz: f64,
        depth: f64,
        restitution: f64,
    ) {
        self.contacts.push(ContactData {
            body_a,
            body_b,
            normal: [nx, ny, nz],
            depth,
            restitution,
        });
    }

    /// Clear all accumulated contacts.
    pub fn clear_contacts(&mut self) {
        self.contacts.clear();
    }

    /// Add an explicit bilateral constraint.
    pub fn add_constraint(
        &mut self,
        body_a: u32,
        body_b: u32,
        nx: f64,
        ny: f64,
        nz: f64,
        bias: f64,
    ) {
        self.constraints.push(WasmConstraint {
            body_a,
            body_b,
            normal: [nx, ny, nz],
            bias,
            lambda: 0.0,
        });
    }

    /// Clear bilateral constraints.
    pub fn clear_constraints(&mut self) {
        self.constraints.clear();
    }

    /// Solve the constraints against the provided engine.
    ///
    /// Updates body velocities in the engine to satisfy the constraints.
    /// Returns the number of constraint rows processed.
    pub fn solve(&mut self, engine: &mut WasmPhysicsEngine, dt: f64) -> u32 {
        self.velocity_corrections.clear();

        // Build constraint rows from contacts
        let mut rows: Vec<WasmConstraint> = Vec::new();
        for c in &self.contacts {
            let bias = (self.baumgarte / dt) * (c.depth - self.slop).max(0.0);
            rows.push(WasmConstraint {
                body_a: c.body_a,
                body_b: c.body_b,
                normal: c.normal,
                bias,
                lambda: 0.0,
            });
        }
        // Also include explicit constraints
        for con in &self.constraints {
            rows.push(con.clone());
        }

        let n_rows = rows.len() as u32;

        // Gauss-Seidel iterations
        for _iter in 0..self.iterations {
            for row in &mut rows {
                let va = engine.get_velocity(row.body_a);
                let vb = engine.get_velocity(row.body_b);

                let nx = row.normal[0];
                let ny = row.normal[1];
                let nz = row.normal[2];

                let rel_vn = (va[0] - vb[0]) * nx + (va[1] - vb[1]) * ny + (va[2] - vb[2]) * nz;

                let delta_lambda = -(rel_vn + row.bias);
                let lambda_old = row.lambda;
                row.lambda = (lambda_old + delta_lambda).max(0.0);
                let d_lam = row.lambda - lambda_old;

                if d_lam.abs() < 1e-12 {
                    continue;
                }

                // Apply velocity correction
                let inv_ma = if engine.body_is_dynamic(row.body_a) {
                    engine.body_inv_mass(row.body_a)
                } else {
                    0.0
                };

                let inv_mb = if engine.body_is_dynamic(row.body_b) {
                    engine.body_inv_mass(row.body_b)
                } else {
                    0.0
                };

                let _ = engine.apply_impulse(
                    row.body_a,
                    d_lam * nx * inv_ma,
                    d_lam * ny * inv_ma,
                    d_lam * nz * inv_ma,
                );
                let _ = engine.apply_impulse(
                    row.body_b,
                    -d_lam * nx * inv_mb,
                    -d_lam * ny * inv_mb,
                    -d_lam * nz * inv_mb,
                );
            }
        }

        n_rows
    }

    /// Number of contacts currently queued.
    pub fn contact_count(&self) -> u32 {
        self.contacts.len() as u32
    }

    /// Number of explicit constraints.
    pub fn constraint_count(&self) -> u32 {
        self.constraints.len() as u32
    }
}

impl Default for WasmConstraintSolver {
    fn default() -> Self {
        Self::new()
    }
}
