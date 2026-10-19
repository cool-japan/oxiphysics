// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Helpers shared by the sand simulation tests.

use super::SandSimulation;
use super::mpm::{Mat3, Vec3, mat_mul};

/// Deterministic linear congruential generator (reproducible test data).
pub(super) struct Lcg(u64);

impl Lcg {
    pub(super) fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Uniform sample in `[0, 1)`.
    pub(super) fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// Uniform sample in `[lo, hi)`.
    pub(super) fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }

    /// Vector with every component in `[lo, hi)`.
    pub(super) fn vec3(&mut self, lo: f64, hi: f64) -> Vec3 {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }

    /// Matrix with every entry in `[lo, hi)`.
    pub(super) fn mat3(&mut self, lo: f64, hi: f64) -> Mat3 {
        [self.vec3(lo, hi), self.vec3(lo, hi), self.vec3(lo, hi)]
    }
}

/// `|a − b| ≤ rel · (1 + max(|a|, |b|))`.
pub(super) fn assert_close(a: f64, b: f64, rel: f64, what: &str) {
    let tol = rel * (1.0 + a.abs().max(b.abs()));
    assert!((a - b).abs() <= tol, "{what}: {a} vs {b} (tolerance {tol})");
}

/// Component-wise [`assert_close`].
pub(super) fn assert_vec_close(a: Vec3, b: Vec3, rel: f64, what: &str) {
    for (k, (x, y)) in a.iter().zip(&b).enumerate() {
        assert_close(*x, *y, rel, &format!("{what}[{k}]"));
    }
}

/// Largest entry of `|a − b|`.
pub(super) fn max_abs_diff(a: &Mat3, b: &Mat3) -> f64 {
    a.iter()
        .flatten()
        .zip(b.iter().flatten())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

/// Rotation matrix about the (normalised) `axis` by `angle` (Rodrigues).
pub(super) fn rotation(axis: Vec3, angle: f64) -> Mat3 {
    let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    let [x, y, z] = [axis[0] / n, axis[1] / n, axis[2] / n];
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

/// `R diag(d) Rᵀ`-style helper: `a · diag(d) · b`.
pub(super) fn sandwich(a: &Mat3, d: [f64; 3], b: &Mat3) -> Mat3 {
    let diag = [[d[0], 0.0, 0.0], [0.0, d[1], 0.0], [0.0, 0.0, d[2]]];
    mat_mul(&mat_mul(a, &diag), b)
}

/// Largest vertical coordinate of the particles.
pub(super) fn max_height(sim: &SandSimulation) -> f64 {
    sim.particles
        .iter()
        .map(|p| p.position[1])
        .fold(f64::NEG_INFINITY, f64::max)
}

/// Mean vertical coordinate of the particles.
pub(super) fn mean_height(sim: &SandSimulation) -> f64 {
    sim.particles.iter().map(|p| p.position[1]).sum::<f64>() / sim.particles.len() as f64
}

/// Largest `|x|` of the particles.
pub(super) fn half_width(sim: &SandSimulation) -> f64 {
    sim.particles
        .iter()
        .map(|p| p.position[0].abs())
        .fold(0.0, f64::max)
}

/// A pre-built triangular heap `|x| ≤ half_base`, `y ≤ (half_base − |x|) tan β`
/// in a one-cell plane-strain slab, filled with sand of density 1600 kg/m³.
pub(super) fn wedge(
    friction_angle_deg: f64,
    slope_deg: f64,
    half_base: f64,
    spacing: f64,
) -> SandSimulation {
    let h = 2.0 * spacing;
    let mut sim = SandSimulation::new(h, [0.0, -9.81, 0.0]);
    sim.friction_angle_deg = friction_angle_deg;
    sim.cfl = 0.6;
    sim.plane_strain_z = Some([0.0, h]);
    let tan_slope = slope_deg.to_radians().tan();
    let mut block = SandSimulation::new(h, [0.0, -9.81, 0.0]);
    let filled = block.fill_box(
        [-half_base, 0.0, 0.0],
        [half_base, half_base * tan_slope, h],
        spacing,
        1600.0,
    );
    assert!(filled.is_ok(), "wedge block must be fillable");
    for p in block.particles {
        if p.position[1] <= (half_base - p.position[0].abs()) * tan_slope {
            sim.add_particle(p);
        }
    }
    sim
}

/// A rectangular column of sand `[−w/2, w/2] × [0, height]` in a one-cell
/// plane-strain slab, density 1600 kg/m³.
pub(super) fn slab_column(
    friction_angle_deg: f64,
    width: f64,
    height: f64,
    spacing: f64,
) -> SandSimulation {
    let h = 2.0 * spacing;
    let mut sim = SandSimulation::new(h, [0.0, -9.81, 0.0]);
    sim.friction_angle_deg = friction_angle_deg;
    sim.cfl = 0.6;
    sim.plane_strain_z = Some([0.0, h]);
    let filled = sim.fill_box(
        [-0.5 * width, 0.0, 0.0],
        [0.5 * width, height, h],
        spacing,
        1600.0,
    );
    assert!(filled.is_ok(), "column must be fillable");
    sim
}
