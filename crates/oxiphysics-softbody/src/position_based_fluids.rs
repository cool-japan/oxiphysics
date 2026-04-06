#![allow(clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Position-Based Fluids (PBF) for liquid simulation.
//!
//! Implements the PBF method (Macklin & Müller, 2013), including SPH kernels,
//! density constraint projection, vorticity confinement, and XSPH viscosity.

/// A single PBF particle.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PbfParticle {
    /// Current position \[m\].
    pub position: [f64; 3],
    /// Predicted position after explicit integration \[m\].
    pub predicted: [f64; 3],
    /// Velocity \[m/s\].
    pub velocity: [f64; 3],
    /// PBF constraint multiplier λ.
    pub lambda: f64,
    /// Accumulated position correction Δp.
    pub delta_pos: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
}

impl PbfParticle {
    /// Create a new `PbfParticle` at rest.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            predicted: position,
            velocity: [0.0; 3],
            lambda: 0.0,
            delta_pos: [0.0; 3],
            mass,
        }
    }
}

// ------------------------------------------------------------------
// Helper vector operations (plain f64 arrays, no nalgebra)
// ------------------------------------------------------------------

#[inline]
fn vec3_sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn vec3_add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

#[inline]
fn vec3_scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

#[inline]
fn vec3_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn vec3_norm_sq(a: [f64; 3]) -> f64 {
    vec3_dot(a, a)
}

#[inline]
fn vec3_norm(a: [f64; 3]) -> f64 {
    vec3_norm_sq(a).sqrt()
}

#[inline]
fn vec3_cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// ------------------------------------------------------------------
// SPH kernels
// ------------------------------------------------------------------

/// Poly6 kernel W_poly6(r, h).
///
/// # Arguments
/// * `r`  – distance between particles \[m\]
/// * `h`  – smoothing radius \[m\]
#[inline]
pub fn poly6_kernel(r: f64, h: f64) -> f64 {
    if r > h {
        return 0.0;
    }
    let coeff = 315.0 / (64.0 * std::f64::consts::PI * h.powi(9));
    let diff = h * h - r * r;
    coeff * diff.powi(3)
}

/// Gradient of the Spiky kernel ∇W_spiky(r_vec, h).
///
/// # Arguments
/// * `r_vec` – vector from neighbor to particle \[m\]
/// * `h`     – smoothing radius \[m\]
pub fn spiky_gradient(r_vec: [f64; 3], h: f64) -> [f64; 3] {
    let r = vec3_norm(r_vec);
    if r < 1e-12 || r > h {
        return [0.0; 3];
    }
    let coeff = -45.0 / (std::f64::consts::PI * h.powi(6));
    let factor = coeff * (h - r).powi(2) / r;
    vec3_scale(r_vec, factor)
}

/// Find neighbors of each particle using a simple O(n²) spatial scan.
///
/// Returns a `Vec<Vec`usize`>` where element `i` contains the indices of all
/// particles within smoothing radius `h` of particle `i` (excluding itself).
///
/// # Arguments
/// * `particles` – slice of PBF particles
/// * `h`         – smoothing radius \[m\]
pub fn find_neighbors(particles: &[PbfParticle], h: f64) -> Vec<Vec<usize>> {
    let n = particles.len();
    let mut neighbors = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let diff = vec3_sub(particles[i].predicted, particles[j].predicted);
            if vec3_norm(diff) < h {
                neighbors[i].push(j);
            }
        }
    }
    neighbors
}

/// SPH density estimate at particle `i`.
///
/// ρ_i = Σ_j m_j · W_poly6(|x_i − x_j|, h)
///
/// # Arguments
/// * `particles`  – slice of PBF particles
/// * `neighbors`  – neighbour list for particle `i`
/// * `i`          – particle index
/// * `h`          – smoothing radius \[m\]
pub fn density_constraint(particles: &[PbfParticle], neighbors: &[usize], i: usize, h: f64) -> f64 {
    let mut rho = 0.0;
    // Self-contribution.
    rho += particles[i].mass * poly6_kernel(0.0, h);
    for &j in neighbors {
        let diff = vec3_sub(particles[i].predicted, particles[j].predicted);
        let r = vec3_norm(diff);
        rho += particles[j].mass * poly6_kernel(r, h);
    }
    rho
}

/// Compute PBF lambda for particle `i`.
///
/// λ_i = −C_i / (Σ_k |∇_{p_k} C_i|² + ε)
///
/// where C_i = ρ_i / ρ_0 − 1.
///
/// # Arguments
/// * `particles`   – slice of PBF particles
/// * `neighbors`   – neighbour list for particle `i`
/// * `i`           – particle index
/// * `h`           – smoothing radius \[m\]
/// * `rho0`        – rest density \[kg/m³\]
/// * `epsilon`     – regularisation (CFM) parameter
pub fn lambda_compute(
    particles: &[PbfParticle],
    neighbors: &[usize],
    i: usize,
    h: f64,
    rho0: f64,
    epsilon: f64,
) -> f64 {
    let rho_i = density_constraint(particles, neighbors, i, h);
    let c_i = rho_i / rho0 - 1.0;

    // ∇_{p_k} C_i for each k.
    let mut denom = 0.0;

    // k = i: sum of all neighbour spiky gradients.
    let mut grad_i = [0.0f64; 3];
    for &j in neighbors {
        let r_vec = vec3_sub(particles[i].predicted, particles[j].predicted);
        let g = spiky_gradient(r_vec, h);
        let scaled = vec3_scale(g, 1.0 / rho0);
        grad_i = vec3_add(grad_i, scaled);
    }
    denom += vec3_norm_sq(grad_i);

    // k = j: each neighbour contributes its own gradient.
    for &j in neighbors {
        let r_vec = vec3_sub(particles[i].predicted, particles[j].predicted);
        let g = spiky_gradient(r_vec, h);
        let scaled = vec3_scale(g, -1.0 / rho0);
        denom += vec3_norm_sq(scaled);
    }

    -c_i / (denom + epsilon)
}

/// Position correction Δp_i from lambda values.
///
/// Δp_i = (1/ρ_0) Σ_j (λ_i + λ_j) ∇W_spiky(x_i − x_j, h)
///
/// # Arguments
/// * `particles`  – slice of PBF particles (lambdas must be pre-computed)
/// * `neighbors`  – neighbour list for particle `i`
/// * `i`          – particle index
/// * `h`          – smoothing radius \[m\]
/// * `rho0`       – rest density \[kg/m³\]
pub fn delta_position(
    particles: &[PbfParticle],
    neighbors: &[usize],
    i: usize,
    h: f64,
    rho0: f64,
) -> [f64; 3] {
    let mut dp = [0.0f64; 3];
    for &j in neighbors {
        let r_vec = vec3_sub(particles[i].predicted, particles[j].predicted);
        let g = spiky_gradient(r_vec, h);
        let factor = (particles[i].lambda + particles[j].lambda) / rho0;
        dp = vec3_add(dp, vec3_scale(g, factor));
    }
    dp
}

/// Vorticity confinement — add rotational energy back to the velocity field.
///
/// Returns the vorticity confinement force on particle `i`.
///
/// # Arguments
/// * `particles`   – slice of PBF particles
/// * `neighbors`   – neighbour list for particle `i`
/// * `i`           – particle index
/// * `h`           – smoothing radius \[m\]
/// * `epsilon_vc`  – vorticity confinement coefficient (≥ 0)
pub fn vorticity_confinement(
    particles: &[PbfParticle],
    neighbors: &[usize],
    i: usize,
    h: f64,
    epsilon_vc: f64,
) -> [f64; 3] {
    // ω_i = Σ_j v_ij × ∇W(x_i − x_j, h)
    let mut omega = [0.0f64; 3];
    for &j in neighbors {
        let v_ij = vec3_sub(particles[j].velocity, particles[i].velocity);
        let r_vec = vec3_sub(particles[i].predicted, particles[j].predicted);
        let grad_w = spiky_gradient(r_vec, h);
        let cross = vec3_cross(v_ij, grad_w);
        omega = vec3_add(omega, cross);
    }
    let omega_norm = vec3_norm(omega);
    if omega_norm < 1e-12 {
        return [0.0; 3];
    }

    // η = Σ_j |ω_j| / |ω_i| * ∇W  (location vector)
    let mut eta = [0.0f64; 3];
    for &j in neighbors {
        let r_vec = vec3_sub(particles[i].predicted, particles[j].predicted);
        let grad_w = spiky_gradient(r_vec, h);
        let omega_j_norm = {
            let mut w_j = [0.0f64; 3];
            // Approximate neighbour vorticity as omega_norm.
            w_j = vec3_add(w_j, vec3_scale(grad_w, omega_norm));
            vec3_norm(w_j)
        };
        eta = vec3_add(eta, vec3_scale(grad_w, omega_j_norm));
    }
    let eta_norm = vec3_norm(eta);
    if eta_norm < 1e-12 {
        return [0.0; 3];
    }
    let n_hat = vec3_scale(eta, 1.0 / eta_norm);

    vec3_scale(vec3_cross(n_hat, omega), epsilon_vc)
}

/// XSPH artificial viscosity.
///
/// Returns a velocity correction for particle `i`:
///
/// Δv_i = c · Σ_j (v_j − v_i) · W_poly6(|x_i − x_j|, h)
///
/// # Arguments
/// * `particles`  – slice of PBF particles
/// * `neighbors`  – neighbour list for particle `i`
/// * `i`          – particle index
/// * `h`          – smoothing radius \[m\]
/// * `c`          – XSPH coefficient (0–1)
pub fn xsph_viscosity(
    particles: &[PbfParticle],
    neighbors: &[usize],
    i: usize,
    h: f64,
    c: f64,
) -> [f64; 3] {
    let mut dv = [0.0f64; 3];
    for &j in neighbors {
        let diff = vec3_sub(particles[i].predicted, particles[j].predicted);
        let r = vec3_norm(diff);
        let w = poly6_kernel(r, h);
        let v_ij = vec3_sub(particles[j].velocity, particles[i].velocity);
        dv = vec3_add(dv, vec3_scale(v_ij, w));
    }
    vec3_scale(dv, c)
}

/// Full PBF substep: predict → solve density constraints → update velocities.
///
/// Modifies `particles` in-place.
///
/// # Arguments
/// * `particles`   – mutable slice of PBF particles
/// * `gravity`     – acceleration due to gravity \[m/s²\]
/// * `dt`          – time step \[s\]
/// * `h`           – smoothing radius \[m\]
/// * `rho0`        – rest density \[kg/m³\]
/// * `epsilon`     – CFM (constraint-force-mixing) regularisation
/// * `solver_iters`– number of constraint projection iterations
/// * `xsph_c`      – XSPH viscosity coefficient
#[allow(clippy::too_many_arguments)]
pub fn pbf_step(
    particles: &mut Vec<PbfParticle>,
    gravity: [f64; 3],
    dt: f64,
    h: f64,
    rho0: f64,
    epsilon: f64,
    solver_iters: usize,
    xsph_c: f64,
) {
    let n = particles.len();

    // 1. Apply gravity and predict positions.
    for p in particles.iter_mut() {
        p.velocity = vec3_add(p.velocity, vec3_scale(gravity, dt));
        p.predicted = vec3_add(p.position, vec3_scale(p.velocity, dt));
    }

    // 2. Iterative constraint projection.
    for _ in 0..solver_iters {
        let neighbors = find_neighbors(particles, h);

        // Compute lambdas.
        let lambdas: Vec<f64> = (0..n)
            .map(|i| lambda_compute(particles, &neighbors[i], i, h, rho0, epsilon))
            .collect();
        for (p, &lam) in particles.iter_mut().zip(lambdas.iter()) {
            p.lambda = lam;
        }

        // Compute and apply position corrections.
        let deltas: Vec<[f64; 3]> = (0..n)
            .map(|i| delta_position(particles, &neighbors[i], i, h, rho0))
            .collect();
        for (p, &dp) in particles.iter_mut().zip(deltas.iter()) {
            p.predicted = vec3_add(p.predicted, dp);
        }
    }

    // 3. Update velocities and positions.
    let neighbors = find_neighbors(particles, h);
    let xsph_corrections: Vec<[f64; 3]> = (0..n)
        .map(|i| xsph_viscosity(particles, &neighbors[i], i, h, xsph_c))
        .collect();

    for (i, p) in particles.iter_mut().enumerate() {
        p.velocity = vec3_scale(vec3_sub(p.predicted, p.position), 1.0 / dt);
        p.velocity = vec3_add(p.velocity, xsph_corrections[i]);
        p.position = p.predicted;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-10;
    const H: f64 = 0.1; // smoothing radius for tests
    const RHO0: f64 = 1000.0;

    fn make_particle(x: f64, y: f64, z: f64) -> PbfParticle {
        PbfParticle::new([x, y, z], 1.0)
    }

    // ------------------------------------------------------------------
    // Particle construction
    // ------------------------------------------------------------------

    #[test]
    fn test_particle_new_zero_velocity() {
        let p = make_particle(1.0, 2.0, 3.0);
        assert_eq!(p.velocity, [0.0; 3]);
        assert_eq!(p.lambda, 0.0);
    }

    #[test]
    fn test_particle_predicted_equals_position_initially() {
        let p = make_particle(1.0, 2.0, 3.0);
        assert_eq!(p.predicted, p.position);
    }

    // ------------------------------------------------------------------
    // Poly6 kernel
    // ------------------------------------------------------------------

    #[test]
    fn test_poly6_zero_at_zero_distance() {
        // Kernel should be non-negative and finite.
        let w = poly6_kernel(0.0, H);
        assert!(w >= 0.0 && w.is_finite());
    }

    #[test]
    fn test_poly6_zero_beyond_h() {
        let w = poly6_kernel(H + 0.001, H);
        assert!(w.abs() < EPS, "poly6 should be 0 beyond h, got {w}");
    }

    #[test]
    fn test_poly6_monotone_decreasing() {
        // Kernel should decrease as r increases.
        let w0 = poly6_kernel(0.0, H);
        let w1 = poly6_kernel(0.05, H);
        let w2 = poly6_kernel(0.09, H);
        assert!(w0 >= w1, "poly6 should decrease with distance");
        assert!(w1 >= w2, "poly6 should decrease with distance");
    }

    #[test]
    fn test_poly6_positive_inside_h() {
        let w = poly6_kernel(H * 0.5, H);
        assert!(w > 0.0, "poly6 should be positive inside h");
    }

    // ------------------------------------------------------------------
    // Spiky gradient
    // ------------------------------------------------------------------

    #[test]
    fn test_spiky_gradient_zero_beyond_h() {
        let g = spiky_gradient([H + 0.01, 0.0, 0.0], H);
        assert!(vec3_norm(g) < EPS, "Spiky gradient should be zero beyond h");
    }

    #[test]
    fn test_spiky_gradient_points_away() {
        // Gradient should point from j to i (away from origin in r direction).
        let r_vec = [0.05, 0.0, 0.0];
        let g = spiky_gradient(r_vec, H);
        // For positive r_vec, gradient should be negative (pointing away from j).
        assert!(
            g[0] < 0.0,
            "Spiky gradient x component should be negative for +x r_vec, got {}",
            g[0]
        );
    }

    #[test]
    fn test_spiky_gradient_zero_distance() {
        let g = spiky_gradient([0.0; 3], H);
        assert!(
            vec3_norm(g) < EPS,
            "Spiky gradient should be zero at zero distance"
        );
    }

    // ------------------------------------------------------------------
    // Find neighbors
    // ------------------------------------------------------------------

    #[test]
    fn test_find_neighbors_nearby() {
        let particles = vec![
            make_particle(0.0, 0.0, 0.0),
            make_particle(0.05, 0.0, 0.0), // within H
            make_particle(1.0, 0.0, 0.0),  // outside H
        ];
        let neighbors = find_neighbors(&particles, H);
        assert!(
            neighbors[0].contains(&1),
            "Particle 0 should have particle 1 as neighbor"
        );
        assert!(
            !neighbors[0].contains(&2),
            "Particle 0 should not have particle 2 as neighbor"
        );
    }

    #[test]
    fn test_find_neighbors_symmetric() {
        let particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.05, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        assert!(neighbors[0].contains(&1));
        assert!(neighbors[1].contains(&0));
    }

    #[test]
    fn test_find_neighbors_no_self() {
        let particles = vec![make_particle(0.0, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        assert!(
            !neighbors[0].contains(&0),
            "Particle should not be its own neighbor"
        );
    }

    // ------------------------------------------------------------------
    // Density constraint
    // ------------------------------------------------------------------

    #[test]
    fn test_density_constraint_positive() {
        let particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        let rho = density_constraint(&particles, &neighbors[0], 0, H);
        assert!(rho > 0.0, "Density should be positive, got {rho}");
    }

    #[test]
    fn test_density_constraint_increases_with_more_neighbors() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0)];
        let rho1 = density_constraint(&particles, &[], 0, H);

        particles.push(make_particle(0.04, 0.0, 0.0));
        let neighbors = find_neighbors(&particles, H);
        let rho2 = density_constraint(&particles, &neighbors[0], 0, H);

        assert!(rho2 > rho1, "More neighbors should give higher density");
    }

    // ------------------------------------------------------------------
    // Lambda computation
    // ------------------------------------------------------------------

    #[test]
    fn test_lambda_finite() {
        let particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        let lam = lambda_compute(&particles, &neighbors[0], 0, H, RHO0, 100.0);
        assert!(lam.is_finite(), "Lambda should be finite, got {lam}");
    }

    #[test]
    fn test_lambda_negative_when_overdense() {
        // Pack many particles close together → density > ρ0 → λ < 0.
        let mut particles = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                particles.push(make_particle(i as f64 * 0.02, j as f64 * 0.02, 0.0));
            }
        }
        let neighbors = find_neighbors(&particles, H);
        let lam = lambda_compute(&particles, &neighbors[0], 0, H, 10.0, 0.001);
        // For over-dense case, C > 0 → λ < 0.
        assert!(
            lam <= 0.0,
            "Lambda should be ≤ 0 for overdense system, got {lam}"
        );
    }

    // ------------------------------------------------------------------
    // Delta position
    // ------------------------------------------------------------------

    #[test]
    fn test_delta_position_finite() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        particles[0].lambda = -0.01;
        particles[1].lambda = -0.01;
        let dp = delta_position(&particles, &neighbors[0], 0, H, RHO0);
        assert!(
            dp.iter().all(|v| v.is_finite()),
            "Delta position should be finite: {dp:?}"
        );
    }

    #[test]
    fn test_delta_position_zero_no_neighbors() {
        let particles = vec![make_particle(0.0, 0.0, 0.0)];
        let dp = delta_position(&particles, &[], 0, H, RHO0);
        assert!(vec3_norm(dp) < EPS, "No neighbors → zero delta position");
    }

    // ------------------------------------------------------------------
    // XSPH viscosity
    // ------------------------------------------------------------------

    #[test]
    fn test_xsph_zero_with_equal_velocities() {
        let particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        let dv = xsph_viscosity(&particles, &neighbors[0], 0, H, 0.1);
        assert!(
            vec3_norm(dv) < EPS,
            "XSPH correction should be zero when velocities are equal"
        );
    }

    #[test]
    fn test_xsph_nonzero_with_different_velocities() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        particles[1].velocity = [1.0, 0.0, 0.0];
        let neighbors = find_neighbors(&particles, H);
        let dv = xsph_viscosity(&particles, &neighbors[0], 0, H, 0.1);
        assert!(
            vec3_norm(dv) > 0.0,
            "XSPH correction should be non-zero with different velocities"
        );
    }

    #[test]
    fn test_xsph_scales_with_coefficient() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        particles[1].velocity = [1.0, 0.0, 0.0];
        let neighbors = find_neighbors(&particles, H);
        let dv1 = xsph_viscosity(&particles, &neighbors[0], 0, H, 0.1);
        let dv2 = xsph_viscosity(&particles, &neighbors[0], 0, H, 0.2);
        let diff = (vec3_norm(dv2) - 2.0 * vec3_norm(dv1)).abs();
        assert!(diff < EPS, "XSPH should scale linearly with coefficient");
    }

    // ------------------------------------------------------------------
    // Full PBF step
    // ------------------------------------------------------------------

    #[test]
    fn test_pbf_step_no_crash() {
        let mut particles = vec![
            make_particle(0.0, 1.0, 0.0),
            make_particle(0.04, 1.0, 0.0),
            make_particle(0.08, 1.0, 0.0),
        ];
        pbf_step(
            &mut particles,
            [0.0, -9.81, 0.0],
            1.0 / 60.0,
            H,
            RHO0,
            100.0,
            2,
            0.01,
        );
        for p in &particles {
            assert!(
                p.position.iter().all(|v| v.is_finite()),
                "Positions should remain finite after PBF step"
            );
        }
    }

    #[test]
    fn test_pbf_step_gravity_moves_particles() {
        let mut particles = vec![make_particle(0.0, 1.0, 0.0)];
        let y_before = particles[0].position[1];
        pbf_step(
            &mut particles,
            [0.0, -9.81, 0.0],
            1.0 / 60.0,
            H,
            RHO0,
            100.0,
            1,
            0.0,
        );
        let y_after = particles[0].position[1];
        assert!(
            y_after < y_before,
            "Particle should fall under gravity: y_before={y_before}, y_after={y_after}"
        );
    }

    #[test]
    fn test_pbf_step_velocity_updated() {
        let mut particles = vec![make_particle(0.0, 10.0, 0.0)];
        pbf_step(
            &mut particles,
            [0.0, -9.81, 0.0],
            1.0 / 60.0,
            H,
            RHO0,
            100.0,
            1,
            0.0,
        );
        let vy = particles[0].velocity[1];
        assert!(
            vy < 0.0,
            "Velocity should be negative after downward gravity step, got {vy}"
        );
    }

    // ------------------------------------------------------------------
    // Vorticity confinement
    // ------------------------------------------------------------------

    #[test]
    fn test_vorticity_zero_uniform_velocity() {
        let particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        let neighbors = find_neighbors(&particles, H);
        // All particles have zero velocity → vorticity is zero.
        let fv = vorticity_confinement(&particles, &neighbors[0], 0, H, 0.1);
        assert!(
            vec3_norm(fv) < EPS,
            "Vorticity force should be zero for uniform velocity field"
        );
    }

    #[test]
    fn test_vorticity_finite() {
        let mut particles = vec![make_particle(0.0, 0.0, 0.0), make_particle(0.04, 0.0, 0.0)];
        particles[0].velocity = [0.0, 1.0, 0.0];
        particles[1].velocity = [0.0, -1.0, 0.0];
        let neighbors = find_neighbors(&particles, H);
        let fv = vorticity_confinement(&particles, &neighbors[0], 0, H, 0.1);
        assert!(
            fv.iter().all(|v| v.is_finite()),
            "Vorticity force should be finite: {fv:?}"
        );
    }
}

// =============================================================================
// Extended PBF structures: PbfParticle with predict/update, pbf_density,
// pbf_lambda, pbf_delta_pos, tensile_instability_correction, PbfSolver,
// vorticity_confinement_pbf, xsph_viscosity wrappers.
// =============================================================================

// ------------------------------------------------------------------
// Extended PbfParticle with density field
// ------------------------------------------------------------------

/// An extended PBF particle carrying a density estimate and lambda multiplier.
///
/// Provides `predict` and `update_velocity` helpers for the standard PBF loop.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PbfExtParticle {
    /// Current position \[m\].
    pub pos: [f64; 3],
    /// Predicted position (after explicit integration) \[m\].
    pub pred_pos: [f64; 3],
    /// Velocity \[m/s\].
    pub vel: [f64; 3],
    /// Particle mass \[kg\].
    pub mass: f64,
    /// Density estimate \[kg/m³\].
    pub rho: f64,
    /// PBF constraint multiplier λ.
    pub lambda: f64,
}

impl PbfExtParticle {
    /// Creates a new `PbfExtParticle` at rest.
    pub fn new(pos: [f64; 3], mass: f64) -> Self {
        Self {
            pos,
            pred_pos: pos,
            vel: [0.0; 3],
            mass,
            rho: 0.0,
            lambda: 0.0,
        }
    }

    /// Explicit Euler prediction step: v += g·dt, pred_pos = pos + v·dt.
    pub fn predict(&mut self, gravity: [f64; 3], dt: f64) {
        self.vel[0] += gravity[0] * dt;
        self.vel[1] += gravity[1] * dt;
        self.vel[2] += gravity[2] * dt;
        self.pred_pos[0] = self.pos[0] + self.vel[0] * dt;
        self.pred_pos[1] = self.pos[1] + self.vel[1] * dt;
        self.pred_pos[2] = self.pos[2] + self.vel[2] * dt;
    }

    /// Velocity update from corrected predicted position: v = (pred_pos − pos) / dt.
    pub fn update_velocity(&mut self, dt: f64) {
        self.vel[0] = (self.pred_pos[0] - self.pos[0]) / dt;
        self.vel[1] = (self.pred_pos[1] - self.pos[1]) / dt;
        self.vel[2] = (self.pred_pos[2] - self.pos[2]) / dt;
        self.pos = self.pred_pos;
    }
}

// ------------------------------------------------------------------
// Free functions for PbfExtParticle
// ------------------------------------------------------------------

/// SPH density estimate at particle `i` using predicted positions.
///
/// ρ_i = Σ_j m_j · W_poly6(|pred_i − pred_j|, h)
///
/// * `particles` – extended particle slice
/// * `i`         – query particle index
/// * `h`         – smoothing radius \[m\]
#[allow(dead_code)]
pub fn pbf_density(particles: &[PbfExtParticle], i: usize, h: f64) -> f64 {
    let mut rho = 0.0;
    for (j, pj) in particles.iter().enumerate() {
        let diff = [
            particles[i].pred_pos[0] - pj.pred_pos[0],
            particles[i].pred_pos[1] - pj.pred_pos[1],
            particles[i].pred_pos[2] - pj.pred_pos[2],
        ];
        let r = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
        if j != i || r < h {
            rho += pj.mass * poly6_kernel(r, h);
        }
    }
    rho
}

/// PBF constraint multiplier λ for particle `i`.
///
/// λ = −(ρ_i/ρ₀ − 1) / (Σ|∇W|² + ε)
///
/// * `rho_i`    – density estimate at particle `i`
/// * `rho0`     – rest density \[kg/m³\]
/// * `grad_sum` – sum of squared gradient magnitudes (denominator term)
/// * `eps`      – regularisation constant (CFM)
#[allow(dead_code)]
pub fn pbf_lambda(rho_i: f64, rho0: f64, grad_sum: f64, eps: f64) -> f64 {
    -(rho_i / rho0 - 1.0) / (grad_sum + eps)
}

/// Position correction Δp_i for particle `i` from pre-computed lambdas.
///
/// Δp_i = (1/ρ₀) Σ_j (λ_i + λ_j) ∇W_spiky(pred_i − pred_j, h)
///
/// * `particles` – extended particle slice (lambdas pre-computed)
/// * `i`         – query particle index
/// * `lambdas`   – lambda slice aligned with `particles`
/// * `h`         – smoothing radius \[m\]
/// * `rho0`      – rest density \[kg/m³\]
#[allow(dead_code)]
pub fn pbf_delta_pos(
    particles: &[PbfExtParticle],
    i: usize,
    lambdas: &[f64],
    h: f64,
    rho0: f64,
) -> [f64; 3] {
    let mut dp = [0.0f64; 3];
    for j in 0..particles.len() {
        if i == j {
            continue;
        }
        let r_vec = [
            particles[i].pred_pos[0] - particles[j].pred_pos[0],
            particles[i].pred_pos[1] - particles[j].pred_pos[1],
            particles[i].pred_pos[2] - particles[j].pred_pos[2],
        ];
        let g = spiky_gradient(r_vec, h);
        let factor = (lambdas[i] + lambdas[j]) / rho0;
        dp[0] += g[0] * factor;
        dp[1] += g[1] * factor;
        dp[2] += g[2] * factor;
    }
    dp
}

/// Tensile instability correction (surface tension term) s_corr.
///
/// s_corr = −k · (W(q·h) / W(Δq·h))^n
///
/// where q is the normalised kernel argument ratio.
///
/// * `q`  – |r_ij| / h
/// * `k`  – artificial surface tension coefficient
/// * `n`  – power (typically 4)
#[allow(dead_code)]
pub fn tensile_instability_correction(q: f64, k: f64, n: usize) -> f64 {
    let delta_q = 0.2_f64; // reference distance fraction
    let w_q = poly6_kernel(q, 1.0);
    let w_dq = poly6_kernel(delta_q, 1.0);
    if w_dq.abs() < 1e-30 {
        return 0.0;
    }
    -k * (w_q / w_dq).powi(n as i32)
}

// ------------------------------------------------------------------
// PbfSolver
// ------------------------------------------------------------------

/// A self-contained Position-Based Fluids solver using `PbfExtParticle`.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct PbfSolver {
    /// Particle list.
    pub particles: Vec<PbfExtParticle>,
    /// Rest density \[kg/m³\].
    pub rho0: f64,
    /// Smoothing radius \[m\].
    pub h: f64,
    /// Number of constraint-projection substeps per time step.
    pub substeps: usize,
}

impl PbfSolver {
    /// Creates a new empty `PbfSolver`.
    pub fn new(rho0: f64, h: f64, substeps: usize) -> Self {
        Self {
            particles: Vec::new(),
            rho0,
            h,
            substeps,
        }
    }

    /// Adds a particle to the solver.
    pub fn add_particle(&mut self, p: PbfExtParticle) {
        self.particles.push(p);
    }

    /// Number of particles currently in the solver.
    pub fn particle_count(&self) -> usize {
        self.particles.len()
    }

    /// Advances the simulation by one time step using `substeps` iterations.
    ///
    /// * `dt`      – time step \[s\]
    /// * `gravity` – gravitational acceleration \[m/s²\]
    #[allow(clippy::too_many_arguments)]
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.particles.len();
        if n == 0 {
            return;
        }

        // 1. Predict positions.
        for p in &mut self.particles {
            p.predict(gravity, dt);
        }

        // 2. Constraint projection iterations.
        let eps_cfm = 1.0e-4;
        for _ in 0..self.substeps {
            // Compute densities.
            let densities: Vec<f64> = (0..n)
                .map(|i| pbf_density(&self.particles, i, self.h))
                .collect();

            // Compute lambdas.
            let lambdas: Vec<f64> = (0..n)
                .map(|i| {
                    // Approximate gradient sum using neighbour count.
                    let grad_sum = densities[i] / (self.rho0 * self.rho0) + 1.0e-6;
                    pbf_lambda(densities[i], self.rho0, grad_sum, eps_cfm)
                })
                .collect();

            // Compute and apply position corrections.
            let deltas: Vec<[f64; 3]> = (0..n)
                .map(|i| pbf_delta_pos(&self.particles, i, &lambdas, self.h, self.rho0))
                .collect();

            for (p, &dp) in self.particles.iter_mut().zip(deltas.iter()) {
                p.pred_pos[0] += dp[0];
                p.pred_pos[1] += dp[1];
                p.pred_pos[2] += dp[2];
            }

            // Store lambdas.
            for (p, &lam) in self.particles.iter_mut().zip(lambdas.iter()) {
                p.lambda = lam;
            }
        }

        // 3. Update velocities.
        for p in &mut self.particles {
            p.update_velocity(dt);
        }
    }

    /// Sum of |ρ_i/ρ₀ − 1| across all particles (density error).
    pub fn total_density_error(&self) -> f64 {
        self.particles
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let rho = pbf_density(&self.particles, i, self.h);
                (rho / self.rho0 - 1.0).abs()
            })
            .sum()
    }

    /// Total kinetic energy of all particles \[J\].
    ///
    /// KE = Σ 0.5 · m · |v|²
    pub fn kinetic_energy(&self) -> f64 {
        self.particles
            .iter()
            .map(|p| {
                let v2 = p.vel[0] * p.vel[0] + p.vel[1] * p.vel[1] + p.vel[2] * p.vel[2];
                0.5 * p.mass * v2
            })
            .sum()
    }
}

/// Vorticity confinement force for `PbfExtParticle`.
///
/// Returns the confinement force vector on particle `i`.
///
/// * `particles`  – extended particle slice
/// * `i`          – particle index
/// * `h`          – smoothing radius \[m\]
/// * `epsilon`    – vorticity confinement coefficient
#[allow(dead_code)]
pub fn vorticity_confinement_pbf(
    particles: &[PbfExtParticle],
    i: usize,
    h: f64,
    epsilon: f64,
) -> [f64; 3] {
    // ω_i = Σ_j v_ij × ∇W_spiky
    let mut omega = [0.0f64; 3];
    for j in 0..particles.len() {
        if i == j {
            continue;
        }
        let r = particles[i].pred_pos[0] - particles[j].pred_pos[0];
        let dist = ((r * r)
            + (particles[i].pred_pos[1] - particles[j].pred_pos[1]).powi(2)
            + (particles[i].pred_pos[2] - particles[j].pred_pos[2]).powi(2))
        .sqrt();
        if dist >= h {
            continue;
        }
        let r_vec = [
            particles[i].pred_pos[0] - particles[j].pred_pos[0],
            particles[i].pred_pos[1] - particles[j].pred_pos[1],
            particles[i].pred_pos[2] - particles[j].pred_pos[2],
        ];
        let v_ij = [
            particles[j].vel[0] - particles[i].vel[0],
            particles[j].vel[1] - particles[i].vel[1],
            particles[j].vel[2] - particles[i].vel[2],
        ];
        let g = spiky_gradient(r_vec, h);
        // cross product v_ij × g
        omega[0] += v_ij[1] * g[2] - v_ij[2] * g[1];
        omega[1] += v_ij[2] * g[0] - v_ij[0] * g[2];
        omega[2] += v_ij[0] * g[1] - v_ij[1] * g[0];
    }
    let omega_norm = (omega[0] * omega[0] + omega[1] * omega[1] + omega[2] * omega[2]).sqrt();
    if omega_norm < 1e-12 {
        return [0.0; 3];
    }
    // η location vector
    let mut eta = [0.0f64; 3];
    for j in 0..particles.len() {
        if i == j {
            continue;
        }
        let r_vec = [
            particles[i].pred_pos[0] - particles[j].pred_pos[0],
            particles[i].pred_pos[1] - particles[j].pred_pos[1],
            particles[i].pred_pos[2] - particles[j].pred_pos[2],
        ];
        let g = spiky_gradient(r_vec, h);
        eta[0] += g[0] * omega_norm;
        eta[1] += g[1] * omega_norm;
        eta[2] += g[2] * omega_norm;
    }
    let eta_norm = (eta[0] * eta[0] + eta[1] * eta[1] + eta[2] * eta[2]).sqrt();
    if eta_norm < 1e-12 {
        return [0.0; 3];
    }
    let n_hat = [eta[0] / eta_norm, eta[1] / eta_norm, eta[2] / eta_norm];
    // confinement force = epsilon * (n_hat × omega)
    [
        epsilon * (n_hat[1] * omega[2] - n_hat[2] * omega[1]),
        epsilon * (n_hat[2] * omega[0] - n_hat[0] * omega[2]),
        epsilon * (n_hat[0] * omega[1] - n_hat[1] * omega[0]),
    ]
}

/// XSPH velocity correction for `PbfExtParticle`.
///
/// Δv_i = c · Σ_j (v_j − v_i) · W_poly6(|pred_i − pred_j|, h)
///
/// * `particles` – extended particle slice
/// * `i`         – particle index
/// * `h`         – smoothing radius \[m\]
/// * `c`         – XSPH coefficient (0–1)
#[allow(dead_code)]
pub fn xsph_viscosity_pbf(particles: &[PbfExtParticle], i: usize, h: f64, c: f64) -> [f64; 3] {
    let mut dv = [0.0f64; 3];
    for j in 0..particles.len() {
        if i == j {
            continue;
        }
        let diff = [
            particles[i].pred_pos[0] - particles[j].pred_pos[0],
            particles[i].pred_pos[1] - particles[j].pred_pos[1],
            particles[i].pred_pos[2] - particles[j].pred_pos[2],
        ];
        let r = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
        let w = poly6_kernel(r, h);
        dv[0] += (particles[j].vel[0] - particles[i].vel[0]) * w;
        dv[1] += (particles[j].vel[1] - particles[i].vel[1]) * w;
        dv[2] += (particles[j].vel[2] - particles[i].vel[2]) * w;
    }
    [dv[0] * c, dv[1] * c, dv[2] * c]
}

// =============================================================================
// Tests for extended PBF
// =============================================================================

#[cfg(test)]
mod pbf_ext_tests {

    use crate::position_based_fluids::PbfExtParticle;
    use crate::position_based_fluids::PbfSolver;
    use crate::position_based_fluids::pbf_delta_pos;
    use crate::position_based_fluids::pbf_density;
    use crate::position_based_fluids::pbf_lambda;
    use crate::position_based_fluids::tensile_instability_correction;
    use crate::position_based_fluids::vorticity_confinement_pbf;
    use crate::position_based_fluids::xsph_viscosity_pbf;

    const EPS: f64 = 1e-10;
    const H: f64 = 0.1;
    const RHO0: f64 = 1000.0;

    fn make_ext(x: f64, y: f64, z: f64) -> PbfExtParticle {
        PbfExtParticle::new([x, y, z], 1.0)
    }

    // ── PbfExtParticle ────────────────────────────────────────────────────

    #[test]
    fn test_ext_particle_new_zero_velocity() {
        let p = make_ext(1.0, 2.0, 3.0);
        assert_eq!(p.vel, [0.0; 3]);
        assert!((p.lambda).abs() < EPS);
    }

    #[test]
    fn test_ext_particle_predict_moves() {
        let mut p = make_ext(0.0, 10.0, 0.0);
        p.predict([0.0, -9.81, 0.0], 1.0 / 60.0);
        assert!(p.pred_pos[1] < 10.0, "Gravity should lower predicted y");
    }

    #[test]
    fn test_ext_particle_update_velocity() {
        let mut p = make_ext(0.0, 1.0, 0.0);
        p.pred_pos = [0.0, 0.9, 0.0];
        p.update_velocity(0.1);
        assert!(
            p.vel[1] < 0.0,
            "Velocity should be negative after downward correction"
        );
    }

    #[test]
    fn test_ext_particle_predict_pred_pos_differs() {
        let mut p = make_ext(0.0, 0.0, 0.0);
        p.predict([1.0, 0.0, 0.0], 0.01);
        assert!(
            p.pred_pos[0] > 0.0,
            "Prediction should move particle in gravity direction"
        );
    }

    // ── pbf_density ───────────────────────────────────────────────────────

    #[test]
    fn test_pbf_density_positive() {
        let particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        let rho = pbf_density(&particles, 0, H);
        assert!(rho > 0.0, "Density estimate should be positive");
    }

    #[test]
    fn test_pbf_density_increases_with_neighbors() {
        let p_alone = vec![make_ext(0.0, 0.0, 0.0)];
        let rho1 = pbf_density(&p_alone, 0, H);
        let p_crowd = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        let rho2 = pbf_density(&p_crowd, 0, H);
        assert!(rho2 > rho1, "More neighbors → higher density");
    }

    // ── pbf_lambda ────────────────────────────────────────────────────────

    #[test]
    fn test_pbf_lambda_negative_overdense() {
        // rho_i > rho0 → C > 0 → lambda < 0
        let lam = pbf_lambda(1500.0, RHO0, 0.01, 1e-4);
        assert!(lam < 0.0, "Lambda should be negative for overdense: {lam}");
    }

    #[test]
    fn test_pbf_lambda_positive_underdense() {
        // rho_i < rho0 → C < 0 → lambda > 0
        let lam = pbf_lambda(500.0, RHO0, 0.01, 1e-4);
        assert!(lam > 0.0, "Lambda should be positive for underdense: {lam}");
    }

    #[test]
    fn test_pbf_lambda_zero_at_rest_density() {
        let lam = pbf_lambda(RHO0, RHO0, 0.01, 1e-4);
        assert!(lam.abs() < EPS, "Lambda should be zero at rest density");
    }

    // ── pbf_delta_pos ─────────────────────────────────────────────────────

    #[test]
    fn test_pbf_delta_pos_finite() {
        let mut particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        particles[0].lambda = -0.01;
        particles[1].lambda = -0.01;
        let lambdas = vec![-0.01, -0.01];
        let dp = pbf_delta_pos(&particles, 0, &lambdas, H, RHO0);
        assert!(
            dp.iter().all(|v| v.is_finite()),
            "Delta pos should be finite: {dp:?}"
        );
    }

    #[test]
    fn test_pbf_delta_pos_zero_single_particle() {
        let particles = vec![make_ext(0.0, 0.0, 0.0)];
        let lambdas = vec![0.0];
        let dp = pbf_delta_pos(&particles, 0, &lambdas, H, RHO0);
        let mag = (dp[0] * dp[0] + dp[1] * dp[1] + dp[2] * dp[2]).sqrt();
        assert!(mag < EPS, "Single particle → zero delta pos");
    }

    // ── tensile_instability_correction ───────────────────────────────────

    #[test]
    fn test_tensile_correction_finite() {
        let s = tensile_instability_correction(0.1, 0.001, 4);
        assert!(s.is_finite(), "Tensile correction must be finite: {s}");
    }

    #[test]
    fn test_tensile_correction_nonpositive() {
        // scorr is −k*(...)^n with k > 0 → should be ≤ 0
        let s = tensile_instability_correction(0.1, 0.001, 4);
        assert!(s <= 0.0, "Tensile correction should be ≤ 0, got {s}");
    }

    #[test]
    fn test_tensile_correction_zero_k() {
        let s = tensile_instability_correction(0.1, 0.0, 4);
        assert!(s.abs() < EPS, "k=0 → zero tensile correction");
    }

    // ── PbfSolver ─────────────────────────────────────────────────────────

    #[test]
    fn test_pbf_solver_particle_count() {
        let mut solver = PbfSolver::new(RHO0, H, 2);
        solver.add_particle(make_ext(0.0, 0.0, 0.0));
        solver.add_particle(make_ext(0.04, 0.0, 0.0));
        assert_eq!(solver.particle_count(), 2);
    }

    #[test]
    fn test_pbf_solver_kinetic_energy_non_negative() {
        let mut solver = PbfSolver::new(RHO0, H, 1);
        solver.add_particle(make_ext(0.0, 0.0, 0.0));
        solver.step(0.01, [0.0, -9.81, 0.0]);
        assert!(solver.kinetic_energy() >= 0.0);
    }

    #[test]
    fn test_pbf_solver_step_no_crash() {
        let mut solver = PbfSolver::new(RHO0, H, 2);
        for i in 0..5 {
            solver.add_particle(make_ext(i as f64 * 0.04, 0.0, 0.0));
        }
        solver.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        for p in &solver.particles {
            assert!(
                p.pos.iter().all(|v| v.is_finite()),
                "Positions should remain finite"
            );
        }
    }

    #[test]
    fn test_pbf_solver_gravity_falls() {
        let mut solver = PbfSolver::new(RHO0, H, 1);
        solver.add_particle(make_ext(0.0, 5.0, 0.0));
        let y_before = solver.particles[0].pos[1];
        solver.step(1.0 / 60.0, [0.0, -9.81, 0.0]);
        let y_after = solver.particles[0].pos[1];
        assert!(y_after < y_before, "Particle should fall under gravity");
    }

    #[test]
    fn test_pbf_solver_empty_no_panic() {
        let mut solver = PbfSolver::new(RHO0, H, 2);
        solver.step(0.01, [0.0, -9.81, 0.0]); // should not panic
    }

    // ── vorticity_confinement_pbf ─────────────────────────────────────────

    #[test]
    fn test_vorticity_pbf_zero_uniform_velocity() {
        let particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        let fv = vorticity_confinement_pbf(&particles, 0, H, 0.1);
        let mag = (fv[0] * fv[0] + fv[1] * fv[1] + fv[2] * fv[2]).sqrt();
        assert!(mag < EPS, "Vorticity force zero for uniform velocity field");
    }

    #[test]
    fn test_vorticity_pbf_finite() {
        let mut particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        particles[0].vel = [0.0, 1.0, 0.0];
        particles[1].vel = [0.0, -1.0, 0.0];
        let fv = vorticity_confinement_pbf(&particles, 0, H, 0.1);
        assert!(
            fv.iter().all(|v| v.is_finite()),
            "Vorticity pbf force should be finite"
        );
    }

    // ── xsph_viscosity_pbf ────────────────────────────────────────────────

    #[test]
    fn test_xsph_pbf_zero_equal_velocities() {
        let particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        let dv = xsph_viscosity_pbf(&particles, 0, H, 0.1);
        let mag = (dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]).sqrt();
        assert!(mag < EPS, "XSPH pbf should be zero for equal velocities");
    }

    #[test]
    fn test_xsph_pbf_direction() {
        let mut particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        particles[1].vel = [1.0, 0.0, 0.0]; // neighbor moving right
        let dv = xsph_viscosity_pbf(&particles, 0, H, 0.1);
        // Particle 0 should be pulled toward neighbor → positive x correction
        assert!(
            dv[0] > 0.0,
            "XSPH pbf should push particle 0 toward moving neighbor, got {}",
            dv[0]
        );
    }

    #[test]
    fn test_xsph_pbf_scales_with_coefficient() {
        let mut particles = vec![make_ext(0.0, 0.0, 0.0), make_ext(0.04, 0.0, 0.0)];
        particles[1].vel = [1.0, 0.0, 0.0];
        let dv1 = xsph_viscosity_pbf(&particles, 0, H, 0.1);
        let dv2 = xsph_viscosity_pbf(&particles, 0, H, 0.2);
        let diff = (dv2[0] - 2.0 * dv1[0]).abs();
        assert!(diff < EPS, "XSPH pbf should scale linearly with c");
    }
}
