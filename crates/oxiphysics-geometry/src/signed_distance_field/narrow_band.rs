// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Narrow-band signed distance field storage.

// ─────────────────────────────────────────────────────────────────────────────
// Narrow-Band SDF
// ─────────────────────────────────────────────────────────────────────────────

/// Narrow-band SDF: only stores values within ±bandwidth of the zero level set.
#[derive(Debug, Clone)]
pub struct NarrowBandSdf {
    /// Grid dimensions.
    pub nx: usize,
    /// Grid dimensions.
    pub ny: usize,
    /// Grid dimensions.
    pub nz: usize,
    /// Grid spacing.
    pub dx: f64,
    /// Origin of the grid.
    pub origin: [f64; 3],
    /// Bandwidth (only |d| < bandwidth stored explicitly).
    pub bandwidth: f64,
    /// SDF values; `f64::MAX` for far-band cells.
    pub values: Vec<f64>,
}

impl NarrowBandSdf {
    /// Build a narrow-band SDF from a function `f`.
    pub fn from_function<F>(
        f: &F,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: f64,
        origin: [f64; 3],
        bandwidth: f64,
    ) -> Self
    where
        F: Fn([f64; 3]) -> f64,
    {
        let mut values = vec![f64::MAX; nx * ny * nz];
        for iz in 0..nz {
            for iy in 0..ny {
                for ix in 0..nx {
                    let p = [
                        origin[0] + ix as f64 * dx,
                        origin[1] + iy as f64 * dx,
                        origin[2] + iz as f64 * dx,
                    ];
                    let d = f(p);
                    if d.abs() <= bandwidth {
                        values[iz * ny * nx + iy * nx + ix] = d;
                    }
                }
            }
        }
        Self {
            nx,
            ny,
            nz,
            dx,
            origin,
            bandwidth,
            values,
        }
    }

    /// Evaluate the narrow-band SDF at grid index (ix, iy, iz).
    pub fn at(&self, ix: usize, iy: usize, iz: usize) -> f64 {
        self.values[iz * self.ny * self.nx + iy * self.nx + ix]
    }

    /// Check if a cell is in the narrow band.
    pub fn in_band(&self, ix: usize, iy: usize, iz: usize) -> bool {
        self.at(ix, iy, iz).abs() < self.bandwidth
    }

    /// Count active (in-band) cells.
    pub fn active_count(&self) -> usize {
        self.values.iter().filter(|&&v| v != f64::MAX).count()
    }
}
