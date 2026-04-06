// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Corotational FEM element and soft body using nalgebra types.

#![allow(clippy::needless_range_loop)]

use oxiphysics_core::math::{Mat3, Real, Vec3};

use crate::particle::SoftParticle;

// ---------------------------------------------------------------------------
// Corotational element
// ---------------------------------------------------------------------------

/// A single tetrahedral corotational FEM element.
#[derive(Debug, Clone)]
pub struct CorotationalElement {
    /// Indices of the four vertices of the tetrahedron.
    pub indices: [usize; 4],
    /// Rest-shape inverse of the edge matrix (Dm^{-1}).
    pub rest_inv: Mat3,
    /// Rest volume.
    pub rest_volume: Real,
    /// Young's modulus.
    pub young: Real,
    /// Poisson's ratio.
    pub poisson: Real,
}

impl CorotationalElement {
    /// Create a new corotational element from the current particle positions.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        indices: [usize; 4],
        particles: &[SoftParticle],
        young: Real,
        poisson: Real,
    ) -> Self {
        let dm = Self::edge_matrix(
            &particles[indices[0]].position,
            &particles[indices[1]].position,
            &particles[indices[2]].position,
            &particles[indices[3]].position,
        );
        let rest_volume = dm.determinant().abs() / 6.0;
        // Safe inversion - if singular, fall back to identity.
        let rest_inv = dm.try_inverse().unwrap_or_else(Mat3::identity);

        Self {
            indices,
            rest_inv,
            rest_volume,
            young,
            poisson,
        }
    }

    /// Compute the 3x3 edge matrix \[e1 | e2 | e3\] where ei = pi - p0.
    fn edge_matrix(p0: &Vec3, p1: &Vec3, p2: &Vec3, p3: &Vec3) -> Mat3 {
        let e1 = p1 - p0;
        let e2 = p2 - p0;
        let e3 = p3 - p0;
        Mat3::new(e1.x, e2.x, e3.x, e1.y, e2.y, e3.y, e1.z, e2.z, e3.z)
    }

    /// Polar decomposition F = R * S. Returns R (the rotation part).
    ///
    /// Uses an iterative method (successive SVD-like correction).
    fn polar_rotation(f: &Mat3) -> Mat3 {
        let mut r = *f;
        // Iterative polar decomposition (a few iterations suffice).
        for _ in 0..10 {
            let r_inv_t = match r.try_inverse() {
                Some(inv) => inv.transpose(),
                None => return Mat3::identity(),
            };
            r = (r + r_inv_t) * 0.5;
        }
        // Ensure determinant is positive (proper rotation).
        if r.determinant() < 0.0 {
            r = -r;
        }
        r
    }

    /// Compute elastic forces on the four vertices and return them as an array.
    ///
    /// The corotational model computes:
    ///   deformation gradient F = Ds * Dm^{-1}
    ///   rotation R = polar(F)
    ///   strain = R^T * Ds * Dm^{-1} - I  (Green-like, linearised)
    ///   stress via Lame parameters
    ///   forces = -V * P * Dm^{-T}  distributed to vertices
    pub fn compute_forces(&self, particles: &[SoftParticle]) -> [Vec3; 4] {
        let [i0, i1, i2, i3] = self.indices;
        let p0 = particles[i0].position;
        let p1 = particles[i1].position;
        let p2 = particles[i2].position;
        let p3 = particles[i3].position;

        // Deformed edge matrix.
        let ds = Self::edge_matrix(&p0, &p1, &p2, &p3);
        let f = ds * self.rest_inv;

        let r = Self::polar_rotation(&f);
        let rt = r.transpose();

        // Strain in the rotated frame.
        let strain = rt * f - Mat3::identity();

        // Lame parameters.
        let mu = self.young / (2.0 * (1.0 + self.poisson));
        let lambda =
            self.young * self.poisson / ((1.0 + self.poisson) * (1.0 - 2.0 * self.poisson));

        // Stress (linear elastic in rotated frame).
        let trace_strain = strain[(0, 0)] + strain[(1, 1)] + strain[(2, 2)];
        let stress = strain * (2.0 * mu) + Mat3::identity() * (lambda * trace_strain);

        // First Piola-Kirchhoff in world frame: P = R * stress.
        let piola = r * stress;

        // Force on edges: H = -V * P * Dm^{-T}.
        let h = piola * self.rest_inv.transpose() * (-self.rest_volume);

        let f1 = Vec3::new(h[(0, 0)], h[(1, 0)], h[(2, 0)]);
        let f2 = Vec3::new(h[(0, 1)], h[(1, 1)], h[(2, 1)]);
        let f3 = Vec3::new(h[(0, 2)], h[(1, 2)], h[(2, 2)]);
        let f0 = -(f1 + f2 + f3);

        [f0, f1, f2, f3]
    }
}

// ---------------------------------------------------------------------------
// FEM soft body
// ---------------------------------------------------------------------------

/// A soft body simulated with corotational FEM and explicit time integration.
#[derive(Debug, Clone)]
pub struct FemSoftBody {
    /// Particles (shared representation with PBD bodies).
    pub particles: Vec<SoftParticle>,
    /// Tetrahedral elements.
    pub elements: Vec<CorotationalElement>,
    /// Global damping coefficient.
    pub damping: Real,
}

impl FemSoftBody {
    /// Create a new FEM soft body.
    pub fn new(
        particles: Vec<SoftParticle>,
        elements: Vec<CorotationalElement>,
        damping: Real,
    ) -> Self {
        Self {
            particles,
            elements,
            damping,
        }
    }

    /// Compute the total kinetic energy: sum of 0.5 * m * |v|^2.
    #[allow(dead_code)]
    pub fn kinetic_energy(&self) -> Real {
        let mut ke = 0.0;
        for p in &self.particles {
            if !p.is_static() {
                let mass = 1.0 / p.inverse_mass;
                ke += 0.5 * mass * p.velocity.norm_squared();
            }
        }
        ke
    }

    /// Compute the gravitational potential energy: sum of m * g * y.
    ///
    /// Uses the Y component of `gravity` as the gravitational acceleration
    /// magnitude (assumes gravity points in -Y).
    #[allow(dead_code)]
    pub fn potential_energy(&self, gravity: &Vec3) -> Real {
        let g_mag = gravity.norm();
        if g_mag < 1e-30 {
            return 0.0;
        }
        let g_dir = gravity / g_mag;
        let mut pe = 0.0;
        for p in &self.particles {
            if !p.is_static() {
                let mass = 1.0 / p.inverse_mass;
                // PE = -m * g . r  (so that lower = less PE)
                pe += -mass * g_mag * g_dir.dot(&p.position);
            }
        }
        pe
    }

    /// Perform one explicit time step.
    pub fn step(&mut self, dt: Real, gravity: &Vec3) {
        let n = self.particles.len();

        // Accumulate forces.
        let mut forces = vec![Vec3::zeros(); n];
        for p_idx in 0..n {
            if !self.particles[p_idx].is_static() {
                let mass = 1.0 / self.particles[p_idx].inverse_mass;
                forces[p_idx] += gravity * mass;
                forces[p_idx] += self.particles[p_idx].external_force;
            }
        }

        // Element forces.
        for elem in &self.elements {
            let f = elem.compute_forces(&self.particles);
            for k in 0..4 {
                forces[elem.indices[k]] += f[k];
            }
        }

        // Integrate (symplectic Euler).
        for i in 0..n {
            let p = &mut self.particles[i];
            if p.is_static() {
                continue;
            }
            let accel = forces[i] * p.inverse_mass;
            p.velocity += accel * dt;
            p.velocity *= 1.0 - self.damping;
            p.position += p.velocity * dt;
        }
    }
}
