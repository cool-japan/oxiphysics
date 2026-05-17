// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! `WasmLbmSim` — D2Q9 Lattice-Boltzmann simulation for the WASM boundary.
//! Note: `EngineWasmLbmConfig` is the engine-internal config; `WasmLbmConfig` in
//! `fluid_bridge` is the canonical public type exposed to JavaScript.

use wasm_bindgen::prelude::*;

use super::Float64View;

/// Configuration for the engine-internal LBM simulation.
///
/// This type is distinct from `fluid_bridge::WasmLbmConfig` which is the
/// canonical JS-facing configuration type.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct EngineWasmLbmConfig {
    /// Grid width (cells in X).
    pub width: usize,
    /// Grid height (cells in Y).
    pub height: usize,
    /// Kinematic viscosity.
    pub viscosity: f64,
}

#[wasm_bindgen]
impl EngineWasmLbmConfig {
    /// Default 64×64 cavity configuration.
    pub fn default_cavity() -> Self {
        Self {
            width: 64,
            height: 64,
            viscosity: 0.01,
        }
    }

    /// Compute the BGK relaxation rate.
    pub fn omega(&self) -> f64 {
        1.0 / (3.0 * self.viscosity + 0.5)
    }

    /// Grid width as `u32` (JS-compatible).
    pub fn width_js(&self) -> u32 {
        self.width as u32
    }

    /// Grid height as `u32` (JS-compatible).
    pub fn height_js(&self) -> u32 {
        self.height as u32
    }
}

// D2Q9 weights / velocities (same as Python side)
const LBM_W: [f64; 9] = [
    4.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 9.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
    1.0 / 36.0,
];
const LBM_EX: [f64; 9] = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0];
const LBM_EY: [f64; 9] = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0];

/// JS-accessible D2Q9 Lattice-Boltzmann simulation.
///
/// Velocity and density fields can be retrieved as `Float64View` or
/// `js_sys::Float64Array` for zero-copy hand-off to JavaScript.
#[wasm_bindgen]
#[derive(Debug, Clone)]
pub struct WasmLbmSim {
    width: usize,
    height: usize,
    omega: f64,
    /// Distribution functions, indexed `[cell * 9 + q]`.
    f: Vec<f64>,
    /// Streaming buffer.
    f_tmp: Vec<f64>,
    step_count: u64,
}

#[wasm_bindgen]
impl WasmLbmSim {
    /// Create a new LBM simulation from configuration.
    pub fn new(config: &EngineWasmLbmConfig) -> Self {
        let n = config.width * config.height;
        let mut f = vec![0.0f64; n * 9];
        for i in 0..n {
            for q in 0..9 {
                f[i * 9 + q] = LBM_W[q];
            }
        }
        Self {
            width: config.width,
            height: config.height,
            omega: config.omega(),
            f: f.clone(),
            f_tmp: f,
            step_count: 0,
        }
    }

    /// Grid width as `u32` (JS-compatible).
    pub fn width(&self) -> u32 {
        self.width as u32
    }

    /// Grid height as `u32` (JS-compatible).
    pub fn height(&self) -> u32 {
        self.height as u32
    }

    /// Number of completed steps as `f64` (JS-compatible; `u64` unsupported by wasm-bindgen).
    pub fn step_count(&self) -> f64 {
        self.step_count as f64
    }

    /// Number of completed steps, low 32 bits.
    pub fn step_count_lo(&self) -> u32 {
        (self.step_count & 0xffff_ffff) as u32
    }

    /// Number of completed steps, high 32 bits.
    pub fn step_count_hi(&self) -> u32 {
        (self.step_count >> 32) as u32
    }

    /// Advance the simulation by one time step.
    pub fn step(&mut self) {
        let w = self.width;
        let h = self.height;
        let omega = self.omega;

        // Collision
        for y in 0..h {
            for x in 0..w {
                let idx = (y * w + x) * 9;
                let rho: f64 = self.f[idx..idx + 9].iter().sum();
                let mut ux = 0.0f64;
                let mut uy = 0.0f64;
                for q in 0..9 {
                    ux += LBM_EX[q] * self.f[idx + q];
                    uy += LBM_EY[q] * self.f[idx + q];
                }
                let rho_inv = if rho > 1e-15 { 1.0 / rho } else { 0.0 };
                ux *= rho_inv;
                uy *= rho_inv;
                let u2 = ux * ux + uy * uy;
                for q in 0..9 {
                    let eu = LBM_EX[q] * ux + LBM_EY[q] * uy;
                    let feq = LBM_W[q] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * u2);
                    self.f_tmp[idx + q] = self.f[idx + q] * (1.0 - omega) + feq * omega;
                }
            }
        }

        // Streaming with periodic BCs
        let src = self.f_tmp.clone();
        for y in 0..h {
            for x in 0..w {
                let dst = (y * w + x) * 9;
                for q in 0..9 {
                    let sx = ((x as isize - LBM_EX[q] as isize).rem_euclid(w as isize)) as usize;
                    let sy = ((y as isize - LBM_EY[q] as isize).rem_euclid(h as isize)) as usize;
                    self.f[dst + q] = src[(sy * w + sx) * 9 + q];
                }
            }
        }

        self.step_count += 1;
    }

    /// Get macroscopic density at `(x, y)`.
    pub fn density_at(&self, x: u32, y: u32) -> f64 {
        let (x, y) = (x as usize, y as usize);
        if x >= self.width || y >= self.height {
            return 0.0;
        }
        let idx = (y * self.width + x) * 9;
        self.f[idx..idx + 9].iter().sum()
    }

    /// Get macroscopic velocity `[ux, uy]` at `(x, y)` as `Vec<f64>` (JS-compatible).
    pub fn velocity_at_js(&self, x: u32, y: u32) -> Vec<f64> {
        let (x, y) = (x as usize, y as usize);
        if x >= self.width || y >= self.height {
            return vec![0.0, 0.0];
        }
        let idx = (y * self.width + x) * 9;
        let rho: f64 = self.f[idx..idx + 9].iter().sum();
        if rho < 1e-15 {
            return vec![0.0, 0.0];
        }
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        for q in 0..9 {
            ux += LBM_EX[q] * self.f[idx + q];
            uy += LBM_EY[q] * self.f[idx + q];
        }
        vec![ux / rho, uy / rho]
    }

    /// Return the full velocity field as a `Float64View` of `[ux, uy]` pairs.
    pub fn velocity_field(&self) -> Float64View {
        let n = self.width * self.height;
        let mut out = vec![0.0f64; n * 2];
        for y in 0..self.height {
            for x in 0..self.width {
                let cell = y * self.width + x;
                let v = self.velocity_at_inner(x, y);
                out[cell * 2] = v[0];
                out[cell * 2 + 1] = v[1];
            }
        }
        Float64View::from_vec(out)
    }

    /// Return the full density field as a `Float64View`.
    pub fn density_field(&self) -> Float64View {
        let n = self.width * self.height;
        let mut out = vec![0.0f64; n];
        for y in 0..self.height {
            for x in 0..self.width {
                out[y * self.width + x] = self.density_at(x as u32, y as u32);
            }
        }
        Float64View::from_vec(out)
    }

    /// Return the full velocity field as a `js_sys::Float64Array`.
    pub fn velocity_field_typed(&self) -> js_sys::Float64Array {
        self.velocity_field().to_typed_array()
    }

    /// Return the full density field as a `js_sys::Float64Array`.
    pub fn density_field_typed(&self) -> js_sys::Float64Array {
        self.density_field().to_typed_array()
    }
}

impl WasmLbmSim {
    /// Get macroscopic velocity `[ux, uy]` at `(x, y)` (Rust-only).
    pub fn velocity_at_inner(&self, x: usize, y: usize) -> [f64; 2] {
        if x >= self.width || y >= self.height {
            return [0.0; 2];
        }
        let idx = (y * self.width + x) * 9;
        let rho: f64 = self.f[idx..idx + 9].iter().sum();
        if rho < 1e-15 {
            return [0.0; 2];
        }
        let mut ux = 0.0f64;
        let mut uy = 0.0f64;
        for q in 0..9 {
            ux += LBM_EX[q] * self.f[idx + q];
            uy += LBM_EY[q] * self.f[idx + q];
        }
        [ux / rho, uy / rho]
    }

    /// Raw step count as `u64` (Rust-only; use `step_count()` from JS).
    pub fn step_count_u64(&self) -> u64 {
        self.step_count
    }

    /// Grid width as `usize` (Rust-only).
    pub fn width_usize(&self) -> usize {
        self.width
    }

    /// Grid height as `usize` (Rust-only).
    pub fn height_usize(&self) -> usize {
        self.height
    }
}
