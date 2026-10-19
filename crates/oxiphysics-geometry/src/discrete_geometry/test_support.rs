// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Meshes and pseudo-random numbers shared by the unit tests of the
//! discrete-geometry modules.

use super::DiscreteMesh;
use crate::test_meshes::{RawMesh, flat_grid};
use std::f64::consts::TAU;

/// Deterministic xorshift64 generator.
pub(super) struct XorShift(u64);

impl XorShift {
    /// A generator for `seed` (a zero seed is replaced by one).
    pub(super) fn new(seed: u64) -> Self {
        Self(seed.max(1))
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform in `[0, 1)`.
    pub(super) fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in `[-1, 1)`.
    pub(super) fn signed(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
}

/// Wrap raw vertices and triangles into a [`DiscreteMesh`].
pub(super) fn discrete(raw: RawMesh) -> DiscreteMesh {
    DiscreteMesh::new(raw.0, raw.1)
}

/// An open cylinder around the z axis with `nu` vertices per ring and `nv + 1`
/// rings (vertex `(i, j)` has index `j * nu + i`), outward counter-clockwise
/// triangles.  It is an annulus: two boundary loops, Euler characteristic zero.
pub(super) fn open_cylinder(nu: usize, nv: usize, radius: f64, height: f64) -> RawMesh {
    let mut vertices = Vec::with_capacity(nu * (nv + 1));
    for j in 0..=nv {
        for i in 0..nu {
            let u = TAU * i as f64 / nu as f64;
            vertices.push([
                radius * u.cos(),
                radius * u.sin(),
                height * j as f64 / nv as f64,
            ]);
        }
    }
    let mut faces = Vec::with_capacity(2 * nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let a = j * nu + i;
            let b = j * nu + (i + 1) % nu;
            let c = (j + 1) * nu + (i + 1) % nu;
            let d = (j + 1) * nu + i;
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (vertices, faces)
}

/// A torus of revolution around the z axis (`major` is the distance of the tube
/// centre from the axis, `minor` the tube radius) with `nu` vertices around the
/// axis and `nv` around the tube (vertex `(i, j)` has index `j * nu + i`),
/// outward counter-clockwise triangles.  Closed, genus one.
pub(super) fn torus(nu: usize, nv: usize, major: f64, minor: f64) -> RawMesh {
    let mut vertices = Vec::with_capacity(nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let u = TAU * i as f64 / nu as f64;
            let v = TAU * j as f64 / nv as f64;
            let ring = major + minor * v.cos();
            vertices.push([ring * u.cos(), ring * u.sin(), minor * v.sin()]);
        }
    }
    let mut faces = Vec::with_capacity(2 * nu * nv);
    for j in 0..nv {
        for i in 0..nu {
            let a = j * nu + i;
            let b = j * nu + (i + 1) % nu;
            let c = ((j + 1) % nv) * nu + (i + 1) % nu;
            let d = ((j + 1) % nv) * nu + i;
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (vertices, faces)
}

/// A unit-spacing `n x n` grid whose vertices are displaced by up to `amplitude`
/// in x and y and up to `amplitude / 2` in z, which produces obtuse triangles on a
/// non-planar surface.
pub(super) fn jittered_grid(n: usize, amplitude: f64, seed: u64) -> RawMesh {
    let (mut vertices, faces) = flat_grid(n, n, 1.0);
    let mut rng = XorShift::new(seed);
    for p in &mut vertices {
        p[0] += amplitude * rng.signed();
        p[1] += amplitude * rng.signed();
        p[2] += 0.5 * amplitude * rng.signed();
    }
    (vertices, faces)
}

/// A cap of the unit sphere around the north pole, reaching the polar angle
/// `max_polar`: vertex 0 is the pole and vertex `1 + (r - 1) * nu + i` is the `i`-th
/// of the `nu` vertices of ring `r` (`1 <= r <= rings`).  Outward counter-clockwise
/// triangles; a disc whose only boundary loop is the last ring.
pub(super) fn spherical_cap(rings: usize, nu: usize, max_polar: f64) -> RawMesh {
    let mut vertices = vec![[0.0, 0.0, 1.0]];
    for r in 1..=rings {
        let theta = max_polar * r as f64 / rings as f64;
        for i in 0..nu {
            let phi = TAU * i as f64 / nu as f64;
            vertices.push([
                theta.sin() * phi.cos(),
                theta.sin() * phi.sin(),
                theta.cos(),
            ]);
        }
    }
    let mut faces = Vec::with_capacity(nu * (2 * rings - 1));
    for i in 0..nu {
        faces.push([0, 1 + i, 1 + (i + 1) % nu]);
    }
    for r in 1..rings {
        for i in 0..nu {
            let a = 1 + (r - 1) * nu + i;
            let b = 1 + r * nu + i;
            let c = 1 + r * nu + (i + 1) % nu;
            let d = 1 + (r - 1) * nu + (i + 1) % nu;
            faces.push([a, b, c]);
            faces.push([a, c, d]);
        }
    }
    (vertices, faces)
}
