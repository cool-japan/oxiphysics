// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Constraint builders and constraint resolution.

#![allow(missing_docs)]

use super::PyPhysicsWorld;

// ===========================================================================
// Constraint Builders
// ===========================================================================

/// Type of constraint between two bodies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)]
pub enum ConstraintType {
    /// Distance constraint: maintain a fixed distance between two anchor points.
    Distance,
    /// Point-to-point (ball-socket): anchor points coincide.
    PointToPoint,
    /// Hinge: bodies rotate around a shared axis.
    Hinge,
    /// Slider: bodies translate along a shared axis.
    Slider,
}

/// A constraint between two rigid bodies.
///
/// Constraints are stored in the world and resolved during each step.
/// Currently implemented as soft position-level corrections.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct PyConstraint {
    /// Unique handle.
    pub handle: u32,
    /// Type of constraint.
    pub constraint_type: ConstraintType,
    /// Handle of the first body.
    pub body_a: u32,
    /// Handle of the second body.
    pub body_b: u32,
    /// Anchor point on body_a in local coordinates (or world if no body).
    pub anchor_a: [f64; 3],
    /// Anchor point on body_b in local coordinates.
    pub anchor_b: [f64; 3],
    /// Target distance (for Distance constraints).
    pub target_distance: f64,
    /// Hinge axis in world space (for Hinge constraints).
    pub axis: [f64; 3],
    /// Constraint stiffness (0..1, 1 = rigid).
    pub stiffness: f64,
    /// Whether the constraint is currently active.
    pub enabled: bool,
}

impl PyConstraint {
    /// Create a distance constraint between two bodies.
    pub fn distance(
        handle: u32,
        body_a: u32,
        body_b: u32,
        anchor_a: [f64; 3],
        anchor_b: [f64; 3],
        distance: f64,
    ) -> Self {
        Self {
            handle,
            constraint_type: ConstraintType::Distance,
            body_a,
            body_b,
            anchor_a,
            anchor_b,
            target_distance: distance,
            axis: [0.0, 1.0, 0.0],
            stiffness: 0.5,
            enabled: true,
        }
    }

    /// Create a point-to-point constraint (ball-socket joint).
    pub fn point_to_point(handle: u32, body_a: u32, body_b: u32, pivot: [f64; 3]) -> Self {
        Self {
            handle,
            constraint_type: ConstraintType::PointToPoint,
            body_a,
            body_b,
            anchor_a: pivot,
            anchor_b: pivot,
            target_distance: 0.0,
            axis: [0.0, 1.0, 0.0],
            stiffness: 0.5,
            enabled: true,
        }
    }

    /// Create a hinge constraint around `axis` at `pivot`.
    pub fn hinge(handle: u32, body_a: u32, body_b: u32, pivot: [f64; 3], axis: [f64; 3]) -> Self {
        Self {
            handle,
            constraint_type: ConstraintType::Hinge,
            body_a,
            body_b,
            anchor_a: pivot,
            anchor_b: pivot,
            target_distance: 0.0,
            axis,
            stiffness: 0.5,
            enabled: true,
        }
    }

    /// Set constraint stiffness (builder pattern).
    pub fn with_stiffness(mut self, s: f64) -> Self {
        self.stiffness = s.clamp(0.0, 1.0);
        self
    }

    /// Enable or disable the constraint.
    pub fn set_enabled(&mut self, enabled: bool) {
        self.enabled = enabled;
    }
}

/// Extension: constraint storage and resolution on `PyPhysicsWorld`.
impl PyPhysicsWorld {
    /// Add a constraint to the world. Returns the constraint handle.
    ///
    /// Note: constraint handles are auto-incremented independently of body handles.
    pub fn add_constraint(&mut self, mut c: PyConstraint) -> u32 {
        let h = self.next_constraint_handle();
        c.handle = h;
        self.constraints.push(c);
        h
    }

    /// Remove a constraint by handle. Returns `true` if found.
    pub fn remove_constraint(&mut self, handle: u32) -> bool {
        let before = self.constraints.len();
        self.constraints.retain(|c| c.handle != handle);
        self.constraints.len() < before
    }

    /// Return the number of active constraints.
    pub fn constraint_count(&self) -> usize {
        self.constraints.len()
    }

    /// Resolve all enabled constraints (soft position correction).
    ///
    /// Called internally by `step()` after contact resolution.
    pub(super) fn resolve_constraints(&mut self) {
        // Clone to avoid borrow issues
        let constraints: Vec<PyConstraint> = self.constraints.clone();
        for c in &constraints {
            if !c.enabled {
                continue;
            }
            match c.constraint_type {
                ConstraintType::Distance | ConstraintType::PointToPoint => {
                    let (pa, pb, inv_ma, inv_mb, static_a, static_b) = {
                        let ba = match self.get_body(c.body_a) {
                            Some(b) => b,
                            None => continue,
                        };
                        let bb = match self.get_body(c.body_b) {
                            Some(b) => b,
                            None => continue,
                        };
                        (
                            ba.position,
                            bb.position,
                            ba.inv_mass(),
                            bb.inv_mass(),
                            ba.is_static,
                            bb.is_static,
                        )
                    };
                    let diff = [pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]];
                    let dist = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
                    let target = if c.constraint_type == ConstraintType::PointToPoint {
                        0.0
                    } else {
                        c.target_distance
                    };
                    if dist < 1e-12 {
                        continue;
                    }
                    let error = dist - target;
                    if error.abs() < 1e-8 {
                        continue;
                    }
                    let n = [diff[0] / dist, diff[1] / dist, diff[2] / dist];
                    let total_inv_m = inv_ma + inv_mb;
                    if total_inv_m < 1e-15 {
                        continue;
                    }
                    let correction = error * c.stiffness / total_inv_m;
                    if !static_a && let Some(ba) = self.get_body_mut(c.body_a) {
                        ba.position[0] += n[0] * correction * inv_ma;
                        ba.position[1] += n[1] * correction * inv_ma;
                        ba.position[2] += n[2] * correction * inv_ma;
                    }
                    if !static_b && let Some(bb) = self.get_body_mut(c.body_b) {
                        bb.position[0] -= n[0] * correction * inv_mb;
                        bb.position[1] -= n[1] * correction * inv_mb;
                        bb.position[2] -= n[2] * correction * inv_mb;
                    }
                }
                ConstraintType::Hinge | ConstraintType::Slider => {
                    // Simplified: treat as point-to-point for now
                    // A full implementation would project motion onto the axis
                }
            }
        }
    }

    /// Compute the next available constraint handle.
    fn next_constraint_handle(&self) -> u32 {
        self.constraints
            .iter()
            .map(|c| c.handle)
            .max()
            .map(|h| h + 1)
            .unwrap_or(0)
    }
}
