// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation statistics.

use super::PyPhysicsWorld;
use pyo3::prelude::*;

// ===========================================================================
// Simulation Statistics
// ===========================================================================

/// Per-step simulation performance and state statistics.
///
/// Retrieved via `PyPhysicsWorld::stats()` after each `step()` call.
#[pyclass(from_py_object)]
#[derive(Debug, Clone, Default)]
pub struct SimStats {
    /// Total number of active (non-removed) bodies.
    #[pyo3(get, set)]
    pub body_count: usize,
    /// Number of bodies currently sleeping.
    #[pyo3(get, set)]
    pub sleeping_count: usize,
    /// Number of bodies that are awake (body_count - sleeping_count).
    #[pyo3(get, set)]
    pub awake_count: usize,
    /// Number of contacts detected in the most recent step.
    #[pyo3(get, set)]
    pub contact_count: usize,
    /// Accumulated simulation time (seconds).
    #[pyo3(get, set)]
    pub simulation_time: f64,
    /// Total kinetic energy summed over all dynamic bodies (½mv²).
    #[pyo3(get, set)]
    pub total_kinetic_energy: f64,
    /// Largest linear speed among all dynamic bodies.
    #[pyo3(get, set)]
    pub max_linear_speed: f64,
    /// Largest angular speed among all dynamic bodies.
    #[pyo3(get, set)]
    pub max_angular_speed: f64,
}

impl PyPhysicsWorld {
    /// Return simulation statistics for the current state.
    pub fn stats(&self) -> SimStats {
        let mut s = SimStats {
            body_count: self.body_count(),
            sleeping_count: self.sleeping_count(),
            contact_count: self.contacts.len(),
            simulation_time: self.time,
            ..SimStats::default()
        };
        s.awake_count = s.body_count.saturating_sub(s.sleeping_count);

        for slot in &self.slots {
            if let Some(body) = slot.body.as_ref() {
                if body.is_static || body.is_kinematic {
                    continue;
                }
                let v2 =
                    body.velocity[0].powi(2) + body.velocity[1].powi(2) + body.velocity[2].powi(2);
                let speed = v2.sqrt();
                s.total_kinetic_energy += 0.5 * body.mass * v2;
                if speed > s.max_linear_speed {
                    s.max_linear_speed = speed;
                }
                let w = (body.angular_velocity[0].powi(2)
                    + body.angular_velocity[1].powi(2)
                    + body.angular_velocity[2].powi(2))
                .sqrt();
                if w > s.max_angular_speed {
                    s.max_angular_speed = w;
                }
            }
        }
        s
    }

    /// Return the kinetic energy of a single body, or `None` if not found.
    pub fn body_kinetic_energy(&self, handle: u32) -> Option<f64> {
        let body = self.get_body(handle)?;
        let v2 = body.velocity[0].powi(2) + body.velocity[1].powi(2) + body.velocity[2].powi(2);
        Some(0.5 * body.mass * v2)
    }

    /// Return the closest active body handle to `point`, or `None` if empty.
    pub fn closest_body(&self, point: [f64; 3]) -> Option<(u32, f64)> {
        let mut best_handle: Option<u32> = None;
        let mut best_dist = f64::MAX;
        for (i, slot) in self.slots.iter().enumerate() {
            if let Some(body) = slot.body.as_ref() {
                let dx = body.position[0] - point[0];
                let dy = body.position[1] - point[1];
                let dz = body.position[2] - point[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < best_dist {
                    best_dist = dist;
                    best_handle = Some(i as u32);
                }
            }
        }
        best_handle.map(|h| (h, best_dist))
    }

    /// Return handles of all bodies within `radius` of `center`.
    pub fn bodies_in_sphere(&self, center: [f64; 3], radius: f64) -> Vec<u32> {
        let r2 = radius * radius;
        self.slots
            .iter()
            .enumerate()
            .filter_map(|(i, slot)| {
                let body = slot.body.as_ref()?;
                let dx = body.position[0] - center[0];
                let dy = body.position[1] - center[1];
                let dz = body.position[2] - center[2];
                if dx * dx + dy * dy + dz * dz <= r2 {
                    Some(i as u32)
                } else {
                    None
                }
            })
            .collect()
    }

    /// Find the pair of bodies with the smallest center-to-center distance.
    ///
    /// Returns `None` if fewer than 2 bodies are present.
    pub fn closest_pair(&self) -> Option<(u32, u32, f64)> {
        let active: Vec<(u32, [f64; 3])> = self
            .slots
            .iter()
            .enumerate()
            .filter_map(|(i, s)| s.body.as_ref().map(|b| (i as u32, b.position)))
            .collect();
        if active.len() < 2 {
            return None;
        }
        let mut best = (0u32, 0u32, f64::MAX);
        for i in 0..active.len() {
            for j in (i + 1)..active.len() {
                let (h_a, p_a) = active[i];
                let (h_b, p_b) = active[j];
                let dx = p_a[0] - p_b[0];
                let dy = p_a[1] - p_b[1];
                let dz = p_a[2] - p_b[2];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                if dist < best.2 {
                    best = (h_a, h_b, dist);
                }
            }
        }
        Some(best)
    }
}
