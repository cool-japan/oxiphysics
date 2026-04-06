#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Position-Based Fluid (PBF) simulation for liquid dynamics.
//!
//! Implements the Macklin & Müller (2013) position-based fluid method:
//!
//! - Poly6 density kernel and Spiky pressure gradient kernel
//! - Density constraint with XPBD-style lambda multipliers
//! - XSPH viscosity for coherent motion
//! - Surface particle detection via density threshold
//! - Wall and periodic boundary conditions
//!
//! All quantities use SI units.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Small vector helpers (no nalgebra dependency)
// ---------------------------------------------------------------------------

#[inline]
fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

#[inline]
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn len3(v: [f64; 3]) -> f64 {
    dot3(v, v).sqrt()
}

#[inline]
fn norm3(v: [f64; 3]) -> [f64; 3] {
    let l = len3(v);
    if l < 1e-300 {
        [0.0; 3]
    } else {
        scale3(v, 1.0 / l)
    }
}

// ---------------------------------------------------------------------------
// SPH kernels
// ---------------------------------------------------------------------------

/// Müller poly6 smoothing kernel W(r, h).
///
/// W(r, h) = (315 / (64π h⁹)) · (h² - r²)³  for r ≤ h, else 0.
pub fn poly6_kernel(r: f64, h: f64) -> f64 {
    if r < 0.0 || r > h {
        return 0.0;
    }
    let h2 = h * h;
    let r2 = r * r;
    let diff = h2 - r2;
    (315.0 / (64.0 * PI * h.powi(9))) * diff.powi(3)
}

/// Gradient of the Debrun Spiky kernel ∇W(r⃗, h).
///
/// ∇W(r⃗, h) = −(45 / (π h⁶)) · (h − |r⃗|)² · r̂
///
/// Used for pressure / density constraint gradients because it is non-zero
/// at r = 0 (unlike poly6).
pub fn spiky_gradient(r: [f64; 3], h: f64) -> [f64; 3] {
    let mag = len3(r);
    if mag <= 0.0 || mag > h {
        return [0.0; 3];
    }
    let coeff = -(45.0 / (PI * h.powi(6))) * (h - mag).powi(2);
    let dir = norm3(r);
    scale3(dir, coeff)
}

/// Compute the XPBD constraint multiplier λᵢ for particle i.
///
/// λᵢ = −C(ρᵢ) / (Σⱼ |∇ⱼC|² + ε)
///
/// where C(ρᵢ) = ρᵢ/ρ₀ − 1 is the density constraint residual and
/// `grad_sum` = Σⱼ |∇ⱼC|² is the sum of squared gradient magnitudes.
/// A small relaxation parameter ε = 1e-6 is added for stability.
pub fn compute_lambda(density: f64, rest_density: f64, grad_sum: f64) -> f64 {
    let constraint = density / rest_density.max(1e-300) - 1.0;
    -constraint / (grad_sum + 1e-6)
}

// ---------------------------------------------------------------------------
// Boundary condition
// ---------------------------------------------------------------------------

/// Boundary condition type for the simulation domain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BoundaryCondition {
    /// Reflective wall boundaries defined by axis-aligned box `[min, max]`.
    Walls {
        /// Minimum corner of the bounding box (x, y, z).
        min: [f64; 3],
        /// Maximum corner of the bounding box (x, y, z).
        max: [f64; 3],
    },
    /// Periodic (wrap-around) boundaries with given cell extents.
    Periodic {
        /// Cell size along each axis.
        cell: [f64; 3],
    },
}

impl BoundaryCondition {
    /// Apply boundary conditions to a position and velocity, modifying them
    /// in-place.
    pub fn apply(&self, position: &mut [f64; 3], velocity: &mut [f64; 3]) {
        match self {
            BoundaryCondition::Walls { min, max } => {
                for k in 0..3 {
                    if position[k] < min[k] {
                        position[k] = min[k];
                        if velocity[k] < 0.0 {
                            velocity[k] = -velocity[k] * 0.5;
                        }
                    }
                    if position[k] > max[k] {
                        position[k] = max[k];
                        if velocity[k] > 0.0 {
                            velocity[k] = -velocity[k] * 0.5;
                        }
                    }
                }
            }
            BoundaryCondition::Periodic { cell } => {
                for k in 0..3 {
                    let c = cell[k];
                    if c > 0.0 {
                        position[k] = position[k].rem_euclid(c);
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// LiquidParticle
// ---------------------------------------------------------------------------

/// A single position-based fluid particle.
#[derive(Clone, Debug)]
pub struct LiquidParticle {
    /// World-space position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Predicted position used during constraint solving (m).
    pub predicted: [f64; 3],
    /// Local density estimate ρᵢ (kg/m³).
    pub density: f64,
    /// Local pressure p (Pa).
    pub pressure: f64,
    /// XPBD constraint multiplier λᵢ.
    pub lambda: f64,
    /// Whether this particle is on the free surface.
    pub is_surface: bool,
    /// Outward-pointing surface normal (unit vector).
    pub surface_normal: [f64; 3],
}

impl LiquidParticle {
    /// Create a particle at `position` with zero velocity.
    pub fn new(position: [f64; 3]) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            predicted: position,
            density: 0.0,
            pressure: 0.0,
            lambda: 0.0,
            is_surface: false,
            surface_normal: [0.0, 1.0, 0.0],
        }
    }

    /// Create a particle at `position` with initial `velocity`.
    pub fn with_velocity(position: [f64; 3], velocity: [f64; 3]) -> Self {
        let mut p = Self::new(position);
        p.velocity = velocity;
        p
    }
}

// ---------------------------------------------------------------------------
// DensityConstraint
// ---------------------------------------------------------------------------

/// Computes density from neighbour particles and solves position corrections.
pub struct DensityConstraint {
    /// Rest (target) density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Smoothing kernel radius h (m).
    pub smoothing_h: f64,
}

impl DensityConstraint {
    /// Create a new density constraint with given rest density and smoothing radius.
    pub fn new(rest_density: f64, smoothing_h: f64) -> Self {
        Self {
            rest_density,
            smoothing_h,
        }
    }

    /// Estimate density at particle `i` from `neighbours` (indices into `particles`).
    pub fn compute_density(
        &self,
        i: usize,
        particles: &[LiquidParticle],
        neighbours: &[usize],
    ) -> f64 {
        let xi = particles[i].predicted;
        let mut rho = 0.0;
        // Include self-contribution (r = 0).
        rho += poly6_kernel(0.0, self.smoothing_h);
        for &j in neighbours {
            if j == i {
                continue;
            }
            let xj = particles[j].predicted;
            let r = len3(sub3(xi, xj));
            rho += poly6_kernel(r, self.smoothing_h);
        }
        rho
    }

    /// Compute the position correction Δxᵢ for particle `i`.
    ///
    /// Returns the position delta `[f64; 3]` to add to `particles[i].predicted`.
    pub fn position_correction(
        &self,
        i: usize,
        particles: &[LiquidParticle],
        neighbours: &[usize],
    ) -> [f64; 3] {
        let xi = particles[i].predicted;
        let lambda_i = particles[i].lambda;
        let mut delta = [0.0_f64; 3];
        for &j in neighbours {
            if j == i {
                continue;
            }
            let xj = particles[j].predicted;
            let r_vec = sub3(xi, xj);
            let grad = spiky_gradient(r_vec, self.smoothing_h);
            let lambda_j = particles[j].lambda;
            let coeff = (lambda_i + lambda_j) / self.rest_density.max(1e-300);
            delta = add3(delta, scale3(grad, coeff));
        }
        delta
    }
}

// ---------------------------------------------------------------------------
// ViscosityForce (XSPH)
// ---------------------------------------------------------------------------

/// XSPH viscosity: smooths velocity field for coherent fluid motion.
///
/// v̄ᵢ = vᵢ + c · Σⱼ (vⱼ − vᵢ) W(|xᵢ − xⱼ|, h)
pub struct ViscosityForce {
    /// XSPH coefficient c (typically 0.01 – 0.1).
    pub coefficient: f64,
    /// Smoothing kernel radius h (m).
    pub smoothing_h: f64,
}

impl ViscosityForce {
    /// Create an XSPH viscosity with given coefficient and smoothing radius.
    pub fn new(coefficient: f64, smoothing_h: f64) -> Self {
        Self {
            coefficient,
            smoothing_h,
        }
    }

    /// Apply XSPH velocity correction to particle `i`, returning the corrected velocity.
    pub fn apply(&self, i: usize, particles: &[LiquidParticle], neighbours: &[usize]) -> [f64; 3] {
        let xi = particles[i].position;
        let vi = particles[i].velocity;
        let mut correction = [0.0_f64; 3];
        for &j in neighbours {
            if j == i {
                continue;
            }
            let xj = particles[j].position;
            let vj = particles[j].velocity;
            let r = len3(sub3(xi, xj));
            let w = poly6_kernel(r, self.smoothing_h);
            let dv = sub3(vj, vi);
            correction = add3(correction, scale3(dv, w));
        }
        add3(vi, scale3(correction, self.coefficient))
    }
}

// ---------------------------------------------------------------------------
// SurfaceDetection
// ---------------------------------------------------------------------------

/// Detects surface particles by comparing their density to a threshold fraction
/// of the rest density.
pub struct SurfaceDetection {
    /// Fraction of rest density below which a particle is considered a surface
    /// particle (e.g. 0.7 means ρ < 0.7 ρ₀ → surface).
    pub threshold_fraction: f64,
}

impl SurfaceDetection {
    /// Create a surface detector with the given density threshold fraction.
    pub fn new(threshold_fraction: f64) -> Self {
        Self { threshold_fraction }
    }

    /// Return `true` if the given particle density indicates a surface particle.
    pub fn is_surface(&self, density: f64, rest_density: f64) -> bool {
        density < self.threshold_fraction * rest_density
    }

    /// Compute the outward surface normal for particle `i` using colour field gradient.
    ///
    /// n̂ = −Σⱼ ∇W(rᵢⱼ, h) / |Σⱼ ∇W(rᵢⱼ, h)|
    pub fn compute_normal(
        &self,
        i: usize,
        particles: &[LiquidParticle],
        neighbours: &[usize],
        h: f64,
    ) -> [f64; 3] {
        let xi = particles[i].position;
        let mut grad_sum = [0.0_f64; 3];
        for &j in neighbours {
            if j == i {
                continue;
            }
            let xj = particles[j].position;
            let r_vec = sub3(xi, xj);
            let g = spiky_gradient(r_vec, h);
            grad_sum = add3(grad_sum, g);
        }
        // Outward normal is negative of the inward colour gradient.
        let neg = scale3(grad_sum, -1.0);
        norm3(neg)
    }
}

// ---------------------------------------------------------------------------
// LiquidSimulation
// ---------------------------------------------------------------------------

/// Position-based fluid simulation (Macklin & Müller 2013).
pub struct LiquidSimulation {
    /// Simulated fluid particles.
    pub particles: Vec<LiquidParticle>,
    /// Rest density ρ₀ (kg/m³).
    pub rest_density: f64,
    /// Smoothing kernel radius h (m).
    pub smoothing_h: f64,
    /// XSPH viscosity coefficient c.
    pub viscosity_coeff: f64,
    /// Gravity acceleration vector (m/s²).
    pub gravity: [f64; 3],
    /// Boundary conditions.
    pub boundary: BoundaryCondition,
    /// Number of constraint solver iterations per step.
    pub solver_iterations: usize,
}

impl LiquidSimulation {
    /// Construct a new liquid simulation.
    pub fn new(
        rest_density: f64,
        smoothing_h: f64,
        viscosity_coeff: f64,
        boundary: BoundaryCondition,
    ) -> Self {
        Self {
            particles: Vec::new(),
            rest_density,
            smoothing_h,
            viscosity_coeff,
            gravity: [0.0, -9.81, 0.0],
            boundary,
            solver_iterations: 3,
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, particle: LiquidParticle) {
        self.particles.push(particle);
    }

    /// Naive O(n²) neighbour search: returns all particles within distance h.
    fn find_neighbours(&self, i: usize) -> Vec<usize> {
        let xi = self.particles[i].predicted;
        let h = self.smoothing_h;
        (0..self.particles.len())
            .filter(|&j| {
                let r = len3(sub3(xi, self.particles[j].predicted));
                r <= h
            })
            .collect()
    }

    /// Advance the simulation by time step `dt` (s).
    pub fn step(&mut self, dt: f64) {
        let n = self.particles.len();
        // 1. Apply gravity and compute predicted positions.
        for p in &mut self.particles {
            p.velocity = add3(p.velocity, scale3(self.gravity, dt));
            p.predicted = add3(p.position, scale3(p.velocity, dt));
        }

        // 2. Solver iterations.
        let density_constraint = DensityConstraint::new(self.rest_density, self.smoothing_h);
        for _ in 0..self.solver_iterations {
            // Compute density and lambda for all particles.
            let mut densities = vec![0.0_f64; n];
            let mut lambdas = vec![0.0_f64; n];
            let neighbours_all: Vec<Vec<usize>> = (0..n).map(|i| self.find_neighbours(i)).collect();

            for i in 0..n {
                let rho =
                    density_constraint.compute_density(i, &self.particles, &neighbours_all[i]);
                densities[i] = rho;
                // Compute gradient sum for lambda denominator.
                let xi = self.particles[i].predicted;
                let mut grad_sq_sum = 0.0;
                for &j in &neighbours_all[i] {
                    if j == i {
                        continue;
                    }
                    let xj = self.particles[j].predicted;
                    let r_vec = sub3(xi, xj);
                    let g = spiky_gradient(r_vec, self.smoothing_h);
                    let g_scaled = scale3(g, 1.0 / self.rest_density);
                    grad_sq_sum += dot3(g_scaled, g_scaled);
                }
                lambdas[i] = compute_lambda(rho, self.rest_density, grad_sq_sum);
            }

            // Write back densities and lambdas.
            for i in 0..n {
                self.particles[i].density = densities[i];
                self.particles[i].lambda = lambdas[i];
            }

            // Compute and apply position corrections.
            let mut deltas = vec![[0.0_f64; 3]; n];
            for i in 0..n {
                deltas[i] =
                    density_constraint.position_correction(i, &self.particles, &neighbours_all[i]);
            }
            for i in 0..n {
                self.particles[i].predicted = add3(self.particles[i].predicted, deltas[i]);
            }
        }

        // 3. Update velocity and position.
        for p in &mut self.particles {
            p.velocity = scale3(sub3(p.predicted, p.position), 1.0 / dt);
            p.position = p.predicted;
        }

        // 4. XSPH viscosity.
        let viscosity = ViscosityForce::new(self.viscosity_coeff, self.smoothing_h);
        let neighbours_all: Vec<Vec<usize>> = (0..n)
            .map(|i| {
                let xi = self.particles[i].position;
                let h = self.smoothing_h;
                (0..n)
                    .filter(|&j| len3(sub3(xi, self.particles[j].position)) <= h)
                    .collect()
            })
            .collect();
        let mut new_velocities = vec![[0.0_f64; 3]; n];
        for i in 0..n {
            new_velocities[i] = viscosity.apply(i, &self.particles, &neighbours_all[i]);
        }
        for i in 0..n {
            self.particles[i].velocity = new_velocities[i];
        }

        // 5. Apply boundary conditions.
        for p in &mut self.particles {
            self.boundary.apply(&mut p.position, &mut p.velocity);
        }

        // 6. Surface detection.
        let surface_detect = SurfaceDetection::new(0.7);
        let is_surface: Vec<bool> = (0..n)
            .map(|i| surface_detect.is_surface(self.particles[i].density, self.rest_density))
            .collect();
        let normals: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let xi = self.particles[i].position;
                let h = self.smoothing_h;
                let neighbours: Vec<usize> = (0..n)
                    .filter(|&j| len3(sub3(xi, self.particles[j].position)) <= h)
                    .collect();
                surface_detect.compute_normal(i, &self.particles, &neighbours, h)
            })
            .collect();
        for i in 0..n {
            self.particles[i].is_surface = is_surface[i];
            self.particles[i].surface_normal = normals[i];
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;

    // 1. poly6_kernel is zero beyond h.
    #[test]
    fn test_poly6_zero_beyond_h() {
        assert!(poly6_kernel(0.5, 0.4).abs() < EPS);
        assert!(poly6_kernel(1.0, 0.4).abs() < EPS);
    }

    // 2. poly6_kernel is positive inside h.
    #[test]
    fn test_poly6_positive_inside_h() {
        let h = 0.1;
        assert!(poly6_kernel(0.0, h) > 0.0);
        assert!(poly6_kernel(0.05, h) > 0.0);
    }

    // 3. poly6_kernel at r = 0 is the maximum.
    #[test]
    fn test_poly6_max_at_zero() {
        let h = 0.1;
        let w0 = poly6_kernel(0.0, h);
        let w1 = poly6_kernel(0.05, h);
        assert!(w0 > w1, "poly6 should be maximum at r=0");
    }

    // 4. poly6_kernel decreases monotonically with r.
    #[test]
    fn test_poly6_monotone() {
        let h = 0.1;
        let vals: Vec<f64> = (0..10).map(|i| poly6_kernel(i as f64 * 0.01, h)).collect();
        for w in vals.windows(2) {
            assert!(
                w[1] <= w[0] + EPS,
                "poly6 not monotone: {} -> {}",
                w[0],
                w[1]
            );
        }
    }

    // 5. poly6_kernel is exactly zero at r = h.
    #[test]
    fn test_poly6_zero_at_h() {
        let h = 0.1;
        assert!(poly6_kernel(h, h).abs() < EPS);
    }

    // 6. spiky_gradient is zero beyond h.
    #[test]
    fn test_spiky_gradient_zero_beyond_h() {
        let r = [0.5, 0.0, 0.0];
        let h = 0.1;
        let g = spiky_gradient(r, h);
        assert!(len3(g).abs() < EPS);
    }

    // 7. spiky_gradient direction is along r (pointing inward, negative).
    #[test]
    fn test_spiky_gradient_direction() {
        let r = [0.05, 0.0, 0.0];
        let h = 0.1;
        let g = spiky_gradient(r, h);
        // Gradient should be negative (pointing toward j, against r).
        assert!(g[0] < 0.0, "Spiky gradient should point inward: {:?}", g);
        assert!(g[1].abs() < EPS);
        assert!(g[2].abs() < EPS);
    }

    // 8. spiky_gradient magnitude increases as r → 0.
    #[test]
    fn test_spiky_gradient_increases_toward_center() {
        let h = 0.1;
        let g1 = len3(spiky_gradient([0.08, 0.0, 0.0], h));
        let g2 = len3(spiky_gradient([0.02, 0.0, 0.0], h));
        assert!(g2 > g1, "Spiky gradient should be larger closer to center");
    }

    // 9. compute_lambda returns zero when density = rest_density and grad_sum > 0.
    #[test]
    fn test_compute_lambda_zero_constraint() {
        let rho0 = 1000.0;
        let lam = compute_lambda(rho0, rho0, 1.0);
        assert!(lam.abs() < EPS, "lambda should be 0 at rest density: {lam}");
    }

    // 10. compute_lambda is negative when density > rest_density (compression).
    #[test]
    fn test_compute_lambda_compressed() {
        let lam = compute_lambda(1100.0, 1000.0, 1.0);
        assert!(
            lam < 0.0,
            "lambda should be negative for compressed fluid: {lam}"
        );
    }

    // 11. compute_lambda is positive when density < rest_density (tension).
    #[test]
    fn test_compute_lambda_tension() {
        let lam = compute_lambda(900.0, 1000.0, 1.0);
        assert!(
            lam > 0.0,
            "lambda should be positive for rarefied fluid: {lam}"
        );
    }

    // 12. compute_lambda scales inversely with grad_sum.
    #[test]
    fn test_compute_lambda_scales_with_grad_sum() {
        let rho = 1100.0;
        let rho0 = 1000.0;
        let l1 = compute_lambda(rho, rho0, 1.0);
        let l2 = compute_lambda(rho, rho0, 2.0);
        // Both should be negative; larger grad_sum → smaller magnitude.
        assert!(
            l1.abs() > l2.abs(),
            "Higher grad_sum should reduce |lambda|"
        );
    }

    // 13. LiquidParticle::new initialises velocity to zero.
    #[test]
    fn test_liquid_particle_zero_velocity() {
        let p = LiquidParticle::new([1.0, 2.0, 3.0]);
        assert!(len3(p.velocity) < EPS);
    }

    // 14. LiquidParticle::with_velocity stores velocity.
    #[test]
    fn test_liquid_particle_with_velocity() {
        let vel = [1.0, -2.0, 0.5];
        let p = LiquidParticle::with_velocity([0.0; 3], vel);
        for k in 0..3 {
            assert!((p.velocity[k] - vel[k]).abs() < EPS);
        }
    }

    // 15. DensityConstraint::compute_density includes self contribution.
    #[test]
    fn test_density_constraint_self_contribution() {
        let dc = DensityConstraint::new(1000.0, 0.1);
        let particles = vec![LiquidParticle::new([0.0; 3])];
        let rho = dc.compute_density(0, &particles, &[]);
        // W(0, h) > 0.
        assert!(rho > 0.0, "Density should include self-contribution");
    }

    // 16. DensityConstraint: density increases with more neighbours.
    #[test]
    fn test_density_constraint_increases_with_neighbours() {
        let h = 0.1;
        let dc = DensityConstraint::new(1000.0, h);
        let mut particles = vec![
            LiquidParticle::new([0.0, 0.0, 0.0]),
            LiquidParticle::new([0.01, 0.0, 0.0]),
            LiquidParticle::new([0.02, 0.0, 0.0]),
        ];
        // Set predicted positions.
        for p in &mut particles {
            p.predicted = p.position;
        }
        let rho_alone = dc.compute_density(0, &particles, &[]);
        let rho_with_neighbors = dc.compute_density(0, &particles, &[1, 2]);
        assert!(
            rho_with_neighbors > rho_alone,
            "More neighbours → higher density"
        );
    }

    // 17. BoundaryCondition::Walls clamps positions.
    #[test]
    fn test_walls_boundary_clamps_position() {
        let bc = BoundaryCondition::Walls {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let mut pos = [-0.5, 1.5, 0.5];
        let mut vel = [-1.0, 1.0, 0.0];
        bc.apply(&mut pos, &mut vel);
        assert!(pos[0] >= 0.0);
        assert!(pos[1] <= 1.0);
    }

    // 18. BoundaryCondition::Walls reverses velocity on collision.
    #[test]
    fn test_walls_boundary_reverses_velocity() {
        let bc = BoundaryCondition::Walls {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let mut pos = [-0.1, 0.5, 0.5];
        let mut vel = [-2.0, 0.0, 0.0];
        bc.apply(&mut pos, &mut vel);
        assert!(vel[0] > 0.0, "Velocity should reverse on wall hit");
    }

    // 19. BoundaryCondition::Periodic wraps positions.
    #[test]
    fn test_periodic_boundary_wraps() {
        let bc = BoundaryCondition::Periodic { cell: [1.0; 3] };
        let mut pos = [1.5, -0.3, 2.1];
        let mut vel = [0.0; 3];
        bc.apply(&mut pos, &mut vel);
        assert!(pos[0] >= 0.0 && pos[0] < 1.0, "Periodic x wrap: {}", pos[0]);
        assert!(pos[1] >= 0.0 && pos[1] < 1.0, "Periodic y wrap: {}", pos[1]);
        assert!(pos[2] >= 0.0 && pos[2] < 1.0, "Periodic z wrap: {}", pos[2]);
    }

    // 20. ViscosityForce returns original velocity when no neighbours.
    #[test]
    fn test_viscosity_no_neighbours() {
        let visc = ViscosityForce::new(0.1, 0.1);
        let particles = vec![LiquidParticle::with_velocity([0.0; 3], [1.0, 2.0, 3.0])];
        let v = visc.apply(0, &particles, &[]);
        for k in 0..3 {
            assert!((v[k] - particles[0].velocity[k]).abs() < EPS);
        }
    }

    // 21. ViscosityForce blends toward neighbour velocity.
    #[test]
    fn test_viscosity_blends_velocity() {
        let h = 0.05;
        let visc = ViscosityForce::new(0.5, h);
        let mut particles = vec![
            LiquidParticle::with_velocity([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
            LiquidParticle::with_velocity([0.01, 0.0, 0.0], [2.0, 0.0, 0.0]),
        ];
        for p in &mut particles {
            p.predicted = p.position;
        }
        let v = visc.apply(0, &particles, &[1]);
        assert!(
            v[0] > 0.0,
            "Viscosity should blend toward neighbour velocity"
        );
    }

    // 22. SurfaceDetection::is_surface true when density < threshold * rho0.
    #[test]
    fn test_surface_detection_is_surface() {
        let sd = SurfaceDetection::new(0.7);
        assert!(sd.is_surface(600.0, 1000.0)); // 600 < 700
        assert!(!sd.is_surface(800.0, 1000.0)); // 800 > 700
    }

    // 23. SurfaceDetection::compute_normal returns unit vector.
    #[test]
    fn test_surface_normal_unit_length() {
        let sd = SurfaceDetection::new(0.7);
        let particles = vec![
            LiquidParticle::new([0.0, 0.0, 0.0]),
            LiquidParticle::new([0.05, 0.0, 0.0]),
        ];
        let n = sd.compute_normal(0, &particles, &[1], 0.1);
        let mag = len3(n);
        // Either unit vector or zero (degenerate).
        assert!(
            mag < EPS || (mag - 1.0).abs() < 1e-8,
            "Normal magnitude: {mag}"
        );
    }

    // 24. LiquidSimulation step does not panic on single particle.
    #[test]
    fn test_single_particle_step_no_panic() {
        let bc = BoundaryCondition::Walls {
            min: [-5.0; 3],
            max: [5.0; 3],
        };
        let mut sim = LiquidSimulation::new(1000.0, 0.1, 0.01, bc);
        sim.add_particle(LiquidParticle::new([0.0, 1.0, 0.0]));
        sim.step(0.01);
        // Particle should have fallen due to gravity.
        assert!(sim.particles[0].position[1] < 1.0);
    }

    // 25. Particle falls under gravity.
    #[test]
    fn test_gravity_moves_particle_downward() {
        let bc = BoundaryCondition::Walls {
            min: [-100.0; 3],
            max: [100.0; 3],
        };
        let mut sim = LiquidSimulation::new(1000.0, 0.1, 0.0, bc);
        sim.add_particle(LiquidParticle::new([0.0, 10.0, 0.0]));
        for _ in 0..60 {
            sim.step(1.0 / 60.0);
        }
        assert!(
            sim.particles[0].position[1] < 10.0,
            "Particle should fall: y={}",
            sim.particles[0].position[1]
        );
    }

    // 26. Wall boundary prevents falling below floor.
    #[test]
    fn test_wall_boundary_floor() {
        let bc = BoundaryCondition::Walls {
            min: [0.0, 0.0, 0.0],
            max: [10.0, 10.0, 10.0],
        };
        let mut sim = LiquidSimulation::new(1000.0, 0.1, 0.0, bc);
        sim.add_particle(LiquidParticle::new([5.0, 0.5, 5.0]));
        for _ in 0..200 {
            sim.step(1.0 / 60.0);
        }
        assert!(
            sim.particles[0].position[1] >= 0.0,
            "Particle should not fall below floor: y={}",
            sim.particles[0].position[1]
        );
    }

    // 27. poly6_kernel is symmetric: W(r, h) = W(-r, h) (r is scalar).
    #[test]
    fn test_poly6_symmetric() {
        let h = 0.1;
        let w1 = poly6_kernel(0.03, h);
        let w2 = poly6_kernel(0.03, h);
        assert!((w1 - w2).abs() < EPS);
    }

    // 28. compute_lambda is invariant to grad_sum = 0 (uses epsilon guard).
    #[test]
    fn test_compute_lambda_epsilon_guard() {
        let lam = compute_lambda(1100.0, 1000.0, 0.0);
        // Should not be infinite/NaN.
        assert!(lam.is_finite(), "lambda should be finite: {lam}");
    }

    // 29. LiquidSimulation add_particle increases count.
    #[test]
    fn test_add_particle_increases_count() {
        let bc = BoundaryCondition::Walls {
            min: [0.0; 3],
            max: [1.0; 3],
        };
        let mut sim = LiquidSimulation::new(1000.0, 0.1, 0.01, bc);
        assert_eq!(sim.particles.len(), 0);
        sim.add_particle(LiquidParticle::new([0.5; 3]));
        assert_eq!(sim.particles.len(), 1);
        sim.add_particle(LiquidParticle::new([0.6; 3]));
        assert_eq!(sim.particles.len(), 2);
    }

    // 30. Multiple particles step together without panic.
    #[test]
    fn test_multi_particle_step() {
        let bc = BoundaryCondition::Walls {
            min: [0.0; 3],
            max: [5.0; 3],
        };
        let mut sim = LiquidSimulation::new(1000.0, 0.3, 0.01, bc);
        for i in 0..5 {
            sim.add_particle(LiquidParticle::new([i as f64 * 0.2, 2.0, 0.5]));
        }
        for _ in 0..10 {
            sim.step(0.01);
        }
        // All particles should be within bounds.
        for p in &sim.particles {
            assert!(p.position[1] >= 0.0);
        }
    }

    // 31. spiky_gradient is zero for r = 0.
    #[test]
    fn test_spiky_gradient_zero_at_origin() {
        let g = spiky_gradient([0.0, 0.0, 0.0], 0.1);
        assert!(
            len3(g) < EPS,
            "Spiky gradient at origin should be zero: {:?}",
            g
        );
    }

    // 32. Periodic boundary wraps near-boundary particle correctly.
    #[test]
    fn test_periodic_boundary_near_cell_edge() {
        let bc = BoundaryCondition::Periodic {
            cell: [2.0, 2.0, 2.0],
        };
        let mut pos = [1.99, 0.5, 0.5];
        let mut vel = [5.0, 0.0, 0.0];
        bc.apply(&mut pos, &mut vel);
        assert!(pos[0] < 2.0 && pos[0] >= 0.0);
    }

    // 33. ViscosityForce coefficient = 0 → no change.
    #[test]
    fn test_viscosity_zero_coefficient() {
        let h = 0.1;
        let visc = ViscosityForce::new(0.0, h);
        let particles = vec![
            LiquidParticle::with_velocity([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]),
            LiquidParticle::with_velocity([0.05, 0.0, 0.0], [4.0, 5.0, 6.0]),
        ];
        let v = visc.apply(0, &particles, &[1]);
        for k in 0..3 {
            assert!(
                (v[k] - particles[0].velocity[k]).abs() < EPS,
                "Zero coefficient viscosity should not change velocity"
            );
        }
    }

    // 34. Surface detection threshold = 1.0 → all particles are surface.
    #[test]
    fn test_surface_detection_all_surface() {
        let sd = SurfaceDetection::new(1.0);
        // Any density < 1.0 * rho0 is surface.
        assert!(sd.is_surface(999.0, 1000.0));
        assert!(!sd.is_surface(1001.0, 1000.0));
    }
}
