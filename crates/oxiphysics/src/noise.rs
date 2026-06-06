// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Procedural 3D noise for turbulence, terrain, and force variation.
//!
//! This module implements **value noise** — a hash-lattice approach that
//! trilinearly interpolates random values at integer grid corners.  Value
//! noise is correct by construction, portable, and adequate for all physics
//! use cases (turbulent force fields, wind gusts, stochastic particle jitter,
//! heightfield terrain).
//!
//! ## Types
//!
//! - `ValueNoise3D` — single-octave 3D value noise, output in `[-1, 1]`.
//! - `FractalNoise` — fractal Brownian motion (fBm) — sums `octaves`
//!   octaves of value noise with `persistence` / `lacunarity` scaling.
//!   Also exposes `FractalNoise::turbulence` and `FractalNoise::ridged`
//!   variants.
//!
//! ## Example
//!
//! ```rust,no_run
//! use oxiphysics::noise::FractalNoise;
//!
//! let noise = FractalNoise::new(42, 4);
//! let turbulence = noise.turbulence(1.0, 2.0, 3.0);
//! assert!(turbulence >= 0.0 && turbulence <= 1.0);
//! ```

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Hash and interpolation helpers — private
// ---------------------------------------------------------------------------

/// Quintic fade: `6t⁵ − 15t⁴ + 10t³` (second derivative is zero at 0 and 1).
#[inline]
fn fade(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

#[inline]
fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}

/// Deterministic hash of integer cell coordinates to a float in `[-1, 1]`.
///
/// Uses a multiply-xorshift mix with the seed to avoid visible axis-aligned
/// banding.
#[inline]
fn hash_cell(ix: i32, iy: i32, iz: i32, seed: u32) -> f64 {
    let mut h = seed;
    h = h.wrapping_add((ix as u32).wrapping_mul(0x6b43_a9b5));
    h = h.wrapping_add((iy as u32).wrapping_mul(0x9e37_79b9));
    h = h.wrapping_add((iz as u32).wrapping_mul(0x517c_c1b7));
    h ^= h >> 16;
    h = h.wrapping_mul(0x45d9_f3b7);
    h ^= h >> 16;
    // Map [0, 2^32) → [-1, 1]
    (h as f64) / (u32::MAX as f64 / 2.0) - 1.0
}

// ---------------------------------------------------------------------------
// ValueNoise3D
// ---------------------------------------------------------------------------

/// Single-octave 3D value noise with quintic interpolation.
///
/// Outputs are in `[-1.0, 1.0]` and continuous everywhere.
///
/// ```rust,no_run
/// use oxiphysics::noise::ValueNoise3D;
///
/// let n = ValueNoise3D::new(0);
/// let v = n.sample(1.5, -2.3, 4.1);
/// assert!(v >= -1.0 && v <= 1.0);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValueNoise3D {
    /// Seed used to initialise the hash function.
    pub seed: u32,
}

impl ValueNoise3D {
    /// Create a new noise instance with the given seed.
    pub fn new(seed: u32) -> Self {
        Self { seed }
    }

    /// Sample noise at `(x, y, z)`.  Returns a value in `[-1.0, 1.0]`.
    pub fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        let ix = x.floor() as i32;
        let iy = y.floor() as i32;
        let iz = z.floor() as i32;
        let fx = x - x.floor();
        let fy = y - y.floor();
        let fz = z - z.floor();

        let ux = fade(fx);
        let uy = fade(fy);
        let uz = fade(fz);

        let seed = self.seed;
        let v000 = hash_cell(ix, iy, iz, seed);
        let v100 = hash_cell(ix + 1, iy, iz, seed);
        let v010 = hash_cell(ix, iy + 1, iz, seed);
        let v110 = hash_cell(ix + 1, iy + 1, iz, seed);
        let v001 = hash_cell(ix, iy, iz + 1, seed);
        let v101 = hash_cell(ix + 1, iy, iz + 1, seed);
        let v011 = hash_cell(ix, iy + 1, iz + 1, seed);
        let v111 = hash_cell(ix + 1, iy + 1, iz + 1, seed);

        let x00 = lerp(v000, v100, ux);
        let x10 = lerp(v010, v110, ux);
        let x01 = lerp(v001, v101, ux);
        let x11 = lerp(v011, v111, ux);

        let y0 = lerp(x00, x10, uy);
        let y1 = lerp(x01, x11, uy);

        lerp(y0, y1, uz)
    }

    /// Sample noise at a world-space position array.
    #[inline]
    pub fn sample_pos(&self, pos: [f64; 3]) -> f64 {
        self.sample(pos[0], pos[1], pos[2])
    }

    /// Sample noise and remap output to `[0.0, 1.0]`.
    #[inline]
    pub fn sample_normalized(&self, x: f64, y: f64, z: f64) -> f64 {
        self.sample(x, y, z) * 0.5 + 0.5
    }

    /// Sample a 3D noise vector — useful for turbulent displacement.
    ///
    /// Uses three independent hash seeds for the three output components.
    pub fn sample_vec3(&self, x: f64, y: f64, z: f64) -> [f64; 3] {
        [
            self.sample(x, y, z),
            ValueNoise3D::new(self.seed.wrapping_add(0x9e37_79b9)).sample(x, y, z),
            ValueNoise3D::new(self.seed.wrapping_add(0x6b43_a9b5)).sample(x, y, z),
        ]
    }
}

// ---------------------------------------------------------------------------
// FractalNoise (fBm)
// ---------------------------------------------------------------------------

/// Fractal Brownian Motion — multiple octaves of [`ValueNoise3D`].
///
/// Each successive octave has its frequency multiplied by `lacunarity` and its
/// amplitude multiplied by `persistence`, then all octaves are summed and
/// normalised to `[-1, 1]`.
///
/// ```rust,no_run
/// use oxiphysics::noise::FractalNoise;
///
/// let mut n = FractalNoise::new(7, 6);
/// n.persistence = 0.5;
/// n.lacunarity  = 2.0;
/// let v = n.sample(0.5, 1.5, 2.5);
/// assert!(v >= -1.0 && v <= 1.0);
/// ```
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FractalNoise {
    /// Number of octaves to sum.  More octaves → finer detail, higher cost.
    pub octaves: u32,
    /// Base frequency (cycles per unit).  Larger = finer base detail.
    pub frequency: f64,
    /// Amplitude multiplier per octave (0 < persistence < 1 for natural looks).
    pub persistence: f64,
    /// Frequency multiplier per octave.  Typically 2.0.
    pub lacunarity: f64,
    noise: ValueNoise3D,
}

impl FractalNoise {
    /// Create with sensible defaults (`frequency=1, persistence=0.5,
    /// lacunarity=2`).
    pub fn new(seed: u32, octaves: u32) -> Self {
        Self {
            octaves: octaves.max(1),
            frequency: 1.0,
            persistence: 0.5,
            lacunarity: 2.0,
            noise: ValueNoise3D::new(seed),
        }
    }

    /// Seed used by the underlying noise.
    pub fn seed(&self) -> u32 {
        self.noise.seed
    }

    /// Sample fBm noise at `(x, y, z)`.  Output is in `[-1, 1]`.
    pub fn sample(&self, x: f64, y: f64, z: f64) -> f64 {
        let mut value = 0.0_f64;
        let mut amplitude = 1.0_f64;
        let mut frequency = self.frequency;
        let mut max_value = 0.0_f64;

        // Use a different seed for each octave to avoid coherent layer overlap.
        let mut oct_seed = self.noise.seed;
        for _ in 0..self.octaves {
            let n = ValueNoise3D::new(oct_seed);
            value += n.sample(x * frequency, y * frequency, z * frequency) * amplitude;
            max_value += amplitude;
            amplitude *= self.persistence;
            frequency *= self.lacunarity;
            oct_seed = oct_seed.wrapping_add(0x9e37_79b9);
        }

        if max_value > 0.0 {
            value / max_value
        } else {
            0.0
        }
    }

    /// Sample at a position array.
    #[inline]
    pub fn sample_pos(&self, pos: [f64; 3]) -> f64 {
        self.sample(pos[0], pos[1], pos[2])
    }

    /// Remap output to `[0.0, 1.0]`.
    #[inline]
    pub fn sample_normalized(&self, x: f64, y: f64, z: f64) -> f64 {
        self.sample(x, y, z) * 0.5 + 0.5
    }

    /// **Turbulence**: absolute value of fBm — creates sharper, cloud-like features.
    ///
    /// Output is in `[0.0, 1.0]`.
    pub fn turbulence(&self, x: f64, y: f64, z: f64) -> f64 {
        let mut value = 0.0_f64;
        let mut amplitude = 1.0_f64;
        let mut frequency = self.frequency;
        let mut max_value = 0.0_f64;

        let mut oct_seed = self.noise.seed;
        for _ in 0..self.octaves {
            let n = ValueNoise3D::new(oct_seed);
            value += n.sample(x * frequency, y * frequency, z * frequency).abs() * amplitude;
            max_value += amplitude;
            amplitude *= self.persistence;
            frequency *= self.lacunarity;
            oct_seed = oct_seed.wrapping_add(0x9e37_79b9);
        }

        if max_value > 0.0 {
            (value / max_value).clamp(0.0, 1.0)
        } else {
            0.0
        }
    }

    /// **Ridged multifractal**: `1 − |fBm|` — creates sharp ridges.
    ///
    /// Output is in `[0.0, 1.0]`.
    pub fn ridged(&self, x: f64, y: f64, z: f64) -> f64 {
        (1.0 - self.turbulence(x, y, z)).clamp(0.0, 1.0)
    }

    /// Sample a turbulence displacement vector (all three axes).
    ///
    /// Useful for particle jitter or wind gusts.
    pub fn turbulence_vec3(&self, x: f64, y: f64, z: f64) -> [f64; 3] {
        // Sample at three offset positions for uncorrelated axes.
        [
            self.sample(x, y, z),
            self.sample(x + 1.7, y + 9.2, z + 3.4),
            self.sample(x + 8.3, y + 2.8, z + 5.1),
        ]
    }

    /// Wrap `sample` to a `[f64;3]` position for use with force-field systems.
    pub fn turbulence_force(&self, pos: [f64; 3], scale: f64) -> [f64; 3] {
        let [tx, ty, tz] = self.turbulence_vec3(pos[0], pos[1], pos[2]);
        [tx * scale, ty * scale, tz * scale]
    }
}
