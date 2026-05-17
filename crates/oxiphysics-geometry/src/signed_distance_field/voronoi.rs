// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Voronoi-based signed distance field.

use super::helpers::{norm3, sub3};

// ─────────────────────────────────────────────────────────────────────────────
// Voronoi SDF
// ─────────────────────────────────────────────────────────────────────────────

/// SDF defined by a set of seed points (Voronoi diagram).
///
/// The SDF value at `p` is the distance to the nearest seed minus the distance
/// to the second-nearest seed, normalized by 2.
#[derive(Debug, Clone)]
pub struct VoronoiSdf {
    /// Seed points.
    pub seeds: Vec<[f64; 3]>,
}

impl VoronoiSdf {
    /// Construct from a list of seed points.
    pub fn new(seeds: Vec<[f64; 3]>) -> Self {
        Self { seeds }
    }

    /// Evaluate the Voronoi SDF at point `p`.
    ///
    /// Returns the distance to the nearest Voronoi cell boundary (positive outside
    /// the nearest cell, negative inside — using the half-plane formulation).
    pub fn evaluate(&self, p: [f64; 3]) -> f64 {
        if self.seeds.is_empty() {
            return f64::MAX;
        }
        let mut d1 = f64::MAX;
        let mut d2 = f64::MAX;
        for &s in &self.seeds {
            let d = norm3(sub3(p, s));
            if d < d1 {
                d2 = d1;
                d1 = d;
            } else if d < d2 {
                d2 = d;
            }
        }
        // Half-plane distance: positive if closer to nearest than second nearest
        (d2 - d1) * 0.5
    }

    /// Index of the nearest seed to point `p`.
    pub fn nearest_seed(&self, p: [f64; 3]) -> usize {
        self.seeds
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                norm3(sub3(p, **a))
                    .partial_cmp(&norm3(sub3(p, **b)))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }

    /// 3D Voronoi cell membership array on a grid.
    pub fn cell_ids(&self, nx: usize, ny: usize, nz: usize, bounds: [f64; 6]) -> Vec<usize> {
        let dx = (bounds[1] - bounds[0]) / nx as f64;
        let dy = (bounds[3] - bounds[2]) / ny as f64;
        let dz = (bounds[5] - bounds[4]) / nz as f64;
        let mut ids = Vec::with_capacity(nx * ny * nz);
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = [
                        bounds[0] + (ix as f64 + 0.5) * dx,
                        bounds[2] + (iy as f64 + 0.5) * dy,
                        bounds[4] + (iz as f64 + 0.5) * dz,
                    ];
                    ids.push(self.nearest_seed(p));
                }
            }
        }
        ids
    }
}
