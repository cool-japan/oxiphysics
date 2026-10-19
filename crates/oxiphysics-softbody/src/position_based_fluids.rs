// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Position-Based Fluids (PBF) for liquid simulation.
//!
//! Implements the PBF method (Macklin & Müller, 2013), including SPH kernels,
//! density constraint projection, vorticity confinement, and XSPH viscosity.
//!
//! Neighbour search uses a uniform spatial hash
//! ([`oxiphysics_geometry::spatial_hash::SpatialHash3D`]) so a step costs
//! O(n · k) instead of O(n²), and the neighbour lists are built **once per
//! step** (after the position prediction) exactly as in Algorithm 1 of the
//! paper. The all-pairs entry points ([`pbf_density`], [`pbf_delta_pos`], ...)
//! are kept for callers that do not have a neighbour list at hand; they share
//! their arithmetic with the neighbour-list variants (`*_nb`).
//!
//! # Mass handling
//!
//! The density constraint is `C_i = ρ_i / ρ₀ − 1` with `ρ_i = Σ_j m_j W(x_i − x_j)`.
//! Its gradients therefore carry the neighbour masses:
//!
//! * `∇_{p_i} C_i = (1/ρ₀) Σ_j m_j ∇W(x_i − x_j)`
//! * `∇_{p_j} C_i = −(m_j/ρ₀) ∇W(x_i − x_j)`
//!
//! and the position correction is `Δp_i = (1/ρ₀) Σ_j (m_j λ_i + m_i λ_j) ∇W(x_i − x_j)`.
//! For unit masses this reduces to the formulas printed in the paper.
//!
//! # Relaxation
//!
//! `λ_i = −C_i / (Σ_k |∇_{p_k} C_i|² + ε)`: the constraint-force-mixing term `ε`
//! is part of the method, not a numerical afterthought. All particles are
//! corrected simultaneously (Jacobi), so each one receives the corrections of
//! every constraint it takes part in; without `ε` of the order of the gradient
//! sum (`≈ k / h²`) the iteration overshoots and expels free-surface
//! particles. [`default_cfm_epsilon`] returns `6 / h²` (`ε = 600` at `h = 0.1`,
//! the value of the paper).

use oxiphysics_geometry::spatial_hash::SpatialHash3D;

/// A single PBF particle.
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

#[inline]
fn vec3_dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    vec3_norm_sq(vec3_sub(a, b))
}

// ------------------------------------------------------------------
// Neighbour search
// ------------------------------------------------------------------

/// Relative enlargement of the spatial-hash cell edge over `h`.
///
/// With a cell edge of exactly `h`, two particles closer than `h` could end
/// up two cells apart after floating-point rounding of the cell index; the
/// margin makes the 27-cell stencil provably sufficient.
const CELL_EDGE_MARGIN: f64 = 1.0e-9;

/// Largest extent (in cells of edge `h`) the hash path is trusted with.
///
/// The hash works on coordinates translated to the cloud's lower corner, so
/// each coordinate carries a rounding error of about `2.2e-16 · extent`. The
/// [`CELL_EDGE_MARGIN`] only absorbs that error while it stays below
/// `1e-9 · h`, i.e. for extents up to a few million cells; larger clouds (and
/// the `i32` cell coordinates of [`SpatialHash3D`] would overflow long before
/// `2e9` cells) fall back to the exhaustive scan, which has neither limit.
const MAX_HASH_CELLS_PER_AXIS: f64 = 1.0e6;

/// Particle neighbour lists in compressed-sparse-row (CSR) layout.
///
/// `neighbors(i)` yields, in ascending order, the indices of all particles
/// `j ≠ i` with `|x_i − x_j| < h`. The relation is symmetric, and particles
/// with a non-finite coordinate neither have neighbours nor appear in anyone's
/// list. Lists are produced with a uniform spatial hash, so building them
/// costs O(n · k) for `k` candidate neighbours per particle; the result is
/// identical to the exhaustive O(n²) scan.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct NeighborCsr {
    /// `offsets[i]..offsets[i + 1]` is the slice of `indices` owned by `i`.
    offsets: Vec<usize>,
    /// Concatenated, per-particle sorted neighbour indices.
    indices: Vec<usize>,
}

impl NeighborCsr {
    /// Builds the neighbour lists of `positions` for smoothing radius `h`.
    ///
    /// A non-positive or NaN `h` yields empty lists (no pair can satisfy
    /// `|r| < h`).
    pub fn from_positions(positions: &[[f64; 3]], h: f64) -> Self {
        Self::build(positions.len(), |i| positions[i], h, false)
    }

    /// Reference implementation: exhaustive O(n²) pair scan.
    ///
    /// Produces exactly the same lists as [`NeighborCsr::from_positions`]; it
    /// exists for validation and as the fallback for point clouds whose extent
    /// is too large for the spatial hash.
    pub fn from_positions_all_pairs(positions: &[[f64; 3]], h: f64) -> Self {
        Self::build(positions.len(), |i| positions[i], h, true)
    }

    /// Number of particles the lists were built for.
    pub fn len(&self) -> usize {
        self.offsets.len().saturating_sub(1)
    }

    /// `true` if the lists were built for zero particles.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Total number of directed `(i, j)` pairs stored (twice the pair count).
    pub fn pair_count(&self) -> usize {
        self.indices.len()
    }

    /// Sorted neighbour indices of particle `i`; empty for an unknown `i`.
    pub fn neighbors(&self, i: usize) -> &[usize] {
        match (self.offsets.get(i), self.offsets.get(i + 1)) {
            (Some(&begin), Some(&end)) => self.indices.get(begin..end).unwrap_or(&[]),
            _ => &[],
        }
    }

    /// Converts to one `Vec<usize>` per particle.
    pub fn to_lists(&self) -> Vec<Vec<usize>> {
        (0..self.len())
            .map(|i| self.neighbors(i).to_vec())
            .collect()
    }

    /// Shared builder: `pos(i)` supplies the position of particle `i`.
    fn build<P: Fn(usize) -> [f64; 3]>(n: usize, pos: P, h: f64, force_all_pairs: bool) -> Self {
        let mut offsets = vec![0usize; n + 1];
        let mut indices: Vec<usize> = Vec::new();
        if n < 2 || h.is_nan() || h <= 0.0 {
            return Self { offsets, indices };
        }
        let h2 = h * h;
        let points: Vec<[f64; 3]> = (0..n).map(&pos).collect();
        let finite: Vec<bool> = points
            .iter()
            .map(|p| p.iter().all(|c| c.is_finite()))
            .collect();

        // Bounding box of the finite points; translating by its lower corner
        // keeps cell indices small and non-negative.
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for (p, &ok) in points.iter().zip(&finite) {
            if ok {
                lo = std::array::from_fn(|k| lo[k].min(p[k]));
                hi = std::array::from_fn(|k| hi[k].max(p[k]));
            }
        }
        if lo[0] > hi[0] {
            // No finite point at all.
            return Self { offsets, indices };
        }
        let extent = (hi[0] - lo[0]).max(hi[1] - lo[1]).max(hi[2] - lo[2]);
        let use_hash = !force_all_pairs && h.is_finite() && extent / h < MAX_HASH_CELLS_PER_AXIS;

        if use_hash {
            let mut hash: SpatialHash3D<()> = SpatialHash3D::new(h * (1.0 + CELL_EDGE_MARGIN));
            for (i, p) in points.iter().enumerate() {
                if finite[i] {
                    hash.insert(i, vec3_sub(*p, lo), ());
                }
            }
            let mut scratch: Vec<usize> = Vec::new();
            for i in 0..n {
                if finite[i] {
                    let p = points[i];
                    scratch.clear();
                    for (j, _) in hash.query_radius(vec3_sub(p, lo), h) {
                        if j != i && vec3_dist_sq(p, points[j]) < h2 {
                            scratch.push(j);
                        }
                    }
                    scratch.sort_unstable();
                    indices.extend_from_slice(&scratch);
                }
                offsets[i + 1] = indices.len();
            }
        } else {
            for i in 0..n {
                if finite[i] {
                    for j in 0..n {
                        if j != i && finite[j] && vec3_dist_sq(points[i], points[j]) < h2 {
                            indices.push(j);
                        }
                    }
                }
                offsets[i + 1] = indices.len();
            }
        }
        Self { offsets, indices }
    }
}

/// Find the neighbours of each particle (by predicted position) with a
/// uniform spatial hash.
///
/// Returns a `Vec<Vec<usize>>` where element `i` contains, in ascending
/// order, the indices of all particles at distance `< h` from particle `i`
/// (excluding itself). The result is identical to an exhaustive O(n²) scan;
/// see [`NeighborCsr`] for the compact representation used internally.
///
/// # Arguments
/// * `particles` – slice of PBF particles
/// * `h`         – smoothing radius \[m\]
pub fn find_neighbors(particles: &[PbfParticle], h: f64) -> Vec<Vec<usize>> {
    NeighborCsr::build(particles.len(), |i| particles[i].predicted, h, false).to_lists()
}

/// Σ_k |∇_{p_k} C_i|² of the density constraint for one particle.
///
/// `pairs` yields `(x_i − x_j, m_j)` for every neighbour `j`. With
/// `g_j = (m_j / ρ₀) ∇W_spiky(x_i − x_j)` the sum is `|Σ_j g_j|² + Σ_j |g_j|²`
/// (the first term is `k = i`, the second the `k = j` contributions).
fn density_gradient_sum<I>(pairs: I, h: f64, rho0: f64) -> f64
where
    I: Iterator<Item = ([f64; 3], f64)>,
{
    let mut grad_i = [0.0f64; 3];
    let mut sum_sq = 0.0;
    for (r_vec, m_j) in pairs {
        let g = vec3_scale(spiky_gradient(r_vec, h), m_j / rho0);
        grad_i = vec3_add(grad_i, g);
        sum_sq += vec3_norm_sq(g);
    }
    vec3_norm_sq(grad_i) + sum_sq
}

/// PBF multiplier `λ = −C / (Σ_k |∇_{p_k} C|² + ε)`.
///
/// Returns `0` when the denominator vanishes (no neighbours and no
/// relaxation): without a constraint gradient the density cannot be corrected,
/// and `−C / 0` would otherwise leak `±∞` / `NaN` into the solver state.
#[inline]
fn constraint_multiplier(c: f64, grad_sum: f64, epsilon: f64) -> f64 {
    let denom = grad_sum + epsilon;
    if denom > 0.0 { -c / denom } else { 0.0 }
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
/// where C_i = ρ_i / ρ_0 − 1 and ρ_i = Σ_j m_j W_poly6. The gradients use the
/// Spiky kernel and carry the neighbour masses (see the module documentation);
/// for unit masses they reduce to the expressions of Macklin & Müller (2013).
///
/// # Arguments
/// * `particles`   – slice of PBF particles
/// * `neighbors`   – neighbour list for particle `i`
/// * `i`           – particle index
/// * `h`           – smoothing radius \[m\]
/// * `rho0`        – rest density \[kg/m³\]
/// * `epsilon`     – regularisation (CFM) parameter, see [`default_cfm_epsilon`]
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
    let denom = density_gradient_sum(
        neighbors.iter().map(|&j| {
            (
                vec3_sub(particles[i].predicted, particles[j].predicted),
                particles[j].mass,
            )
        }),
        h,
        rho0,
    );
    constraint_multiplier(c_i, denom, epsilon)
}

/// Position correction Δp_i from lambda values.
///
/// Δp_i = (1/ρ_0) Σ_j (m_j λ_i + m_i λ_j) ∇W_spiky(x_i − x_j, h)
///
/// This is `Σ_c λ_c ∇_{p_i} C_c` over the constraints that depend on `p_i`
/// (its own and its neighbours'); for unit masses it is the familiar
/// `(1/ρ_0) Σ_j (λ_i + λ_j) ∇W`.
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
        let factor = (particles[j].mass * particles[i].lambda
            + particles[i].mass * particles[j].lambda)
            / rho0;
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
/// Follows Algorithm 1 of Macklin & Müller (2013): the neighbour lists are
/// built **once**, from the predicted positions, and reused by every solver
/// iteration and by the XSPH pass (kernels vanish beyond `h`, so pairs that
/// drift apart during the iterations simply stop contributing).
///
/// The XSPH correction smooths the velocities *produced by* the constraint
/// solve, `v_i ← (x*_i − x_i)/Δt + c Σ_j (v_j − v_i) W(x*_i − x*_j)`, with all
/// corrections evaluated from the same velocity snapshot.
///
/// Modifies `particles` in-place. A non-positive or NaN `dt` leaves them
/// untouched.
///
/// # Arguments
/// * `particles`   – mutable slice of PBF particles
/// * `gravity`     – acceleration due to gravity \[m/s²\]
/// * `dt`          – time step \[s\]
/// * `h`           – smoothing radius \[m\]
/// * `rho0`        – rest density \[kg/m³\]
/// * `epsilon`     – CFM (constraint-force-mixing) relaxation \[1/m²\]; see
///   [`default_cfm_epsilon`] (a value ≪ Σ|∇C|² makes the Jacobi iteration unstable)
/// * `solver_iters`– number of constraint projection iterations
/// * `xsph_c`      – XSPH viscosity coefficient
pub fn pbf_step(
    particles: &mut [PbfParticle],
    gravity: [f64; 3],
    dt: f64,
    h: f64,
    rho0: f64,
    epsilon: f64,
    solver_iters: usize,
    xsph_c: f64,
) {
    let n = particles.len();
    if n == 0 || dt.is_nan() || dt <= 0.0 {
        return;
    }

    // 1. Apply gravity and predict positions.
    for p in particles.iter_mut() {
        p.velocity = vec3_add(p.velocity, vec3_scale(gravity, dt));
        p.predicted = vec3_add(p.position, vec3_scale(p.velocity, dt));
    }

    // 2. Neighbour search, once per step.
    let neighbors = {
        let snapshot: &[PbfParticle] = particles;
        NeighborCsr::build(n, |i| snapshot[i].predicted, h, false)
    };

    // 3. Iterative constraint projection.
    for _ in 0..solver_iters {
        let lambdas: Vec<f64> = (0..n)
            .map(|i| lambda_compute(particles, neighbors.neighbors(i), i, h, rho0, epsilon))
            .collect();
        for (p, &lam) in particles.iter_mut().zip(lambdas.iter()) {
            p.lambda = lam;
        }

        let deltas: Vec<[f64; 3]> = (0..n)
            .map(|i| delta_position(particles, neighbors.neighbors(i), i, h, rho0))
            .collect();
        for (p, &dp) in particles.iter_mut().zip(deltas.iter()) {
            p.predicted = vec3_add(p.predicted, dp);
        }
    }

    // 4. Velocities from the corrected positions, then XSPH on those velocities.
    let inv_dt = 1.0 / dt;
    for p in particles.iter_mut() {
        p.velocity = vec3_scale(vec3_sub(p.predicted, p.position), inv_dt);
    }
    let xsph_corrections: Vec<[f64; 3]> = (0..n)
        .map(|i| xsph_viscosity(particles, neighbors.neighbors(i), i, h, xsph_c))
        .collect();
    for (p, &dv) in particles.iter_mut().zip(xsph_corrections.iter()) {
        p.velocity = vec3_add(p.velocity, dv);
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

/// Shared SPH density sum over the candidate neighbours `js` of particle `i`.
///
/// Pairs with `r² ≥ h²` are rejected before the square root (the kernel is
/// zero there anyway); the self term is always included for `h > 0`.
fn density_over<I: Iterator<Item = usize>>(
    particles: &[PbfExtParticle],
    i: usize,
    h: f64,
    js: I,
) -> f64 {
    let Some(pi) = particles.get(i) else {
        return 0.0;
    };
    if h.is_nan() || h <= 0.0 {
        return 0.0;
    }
    let h2 = h * h;
    let mut rho = pi.mass * poly6_kernel(0.0, h);
    for j in js {
        if j == i {
            continue;
        }
        let Some(pj) = particles.get(j) else {
            continue;
        };
        let d2 = vec3_dist_sq(pi.pred_pos, pj.pred_pos);
        if d2 < h2 {
            rho += pj.mass * poly6_kernel(d2.sqrt(), h);
        }
    }
    rho
}

/// SPH density estimate at particle `i` using predicted positions (all-pairs).
///
/// ρ_i = Σ_j m_j · W_poly6(|pred_i − pred_j|, h)
///
/// This scans every particle; use [`pbf_density_nb`] with a neighbour list
/// (see [`NeighborCsr`]) inside solver loops. Both return identical values.
///
/// * `particles` – extended particle slice
/// * `i`         – query particle index
/// * `h`         – smoothing radius \[m\]
pub fn pbf_density(particles: &[PbfExtParticle], i: usize, h: f64) -> f64 {
    density_over(particles, i, h, 0..particles.len())
}

/// SPH density estimate at particle `i` summed over the neighbour list `nbrs`.
///
/// Identical to [`pbf_density`] provided `nbrs` contains every particle
/// within `h` of particle `i`; the particle's own contribution is always
/// included and an `i` appearing in `nbrs` is ignored.
///
/// * `particles` – extended particle slice
/// * `i`         – query particle index
/// * `nbrs`      – neighbour indices of `i` (e.g. [`NeighborCsr::neighbors`])
/// * `h`         – smoothing radius \[m\]
pub fn pbf_density_nb(particles: &[PbfExtParticle], i: usize, nbrs: &[usize], h: f64) -> f64 {
    density_over(particles, i, h, nbrs.iter().copied())
}

/// Σ_k |∇_{p_k} C_i|² for the density constraint of particle `i`, summed
/// over the candidate neighbours `js`.
fn gradient_sum_over<I: Iterator<Item = usize>>(
    particles: &[PbfExtParticle],
    i: usize,
    h: f64,
    rho0: f64,
    js: I,
) -> f64 {
    let Some(pi) = particles.get(i) else {
        return 0.0;
    };
    if h.is_nan() || h <= 0.0 {
        return 0.0;
    }
    let h2 = h * h;
    density_gradient_sum(
        js.filter(|&j| j != i).filter_map(|j| {
            let pj = particles.get(j)?;
            let r_vec = vec3_sub(pi.pred_pos, pj.pred_pos);
            (vec3_norm_sq(r_vec) < h2).then_some((r_vec, pj.mass))
        }),
        h,
        rho0,
    )
}

/// Σ_k |∇_{p_k} C_i|² of the density constraint `C_i = ρ_i/ρ₀ − 1` (all-pairs).
///
/// With `g_j = (m_j/ρ₀) ∇W_spiky(pred_i − pred_j)` this equals
/// `|Σ_j g_j|² + Σ_j |g_j|²`: the first term is the gradient with respect to
/// `p_i`, the second the gradients with respect to each neighbour. It is the
/// `grad_sum` argument of [`pbf_lambda`].
///
/// * `particles` – extended particle slice
/// * `i`         – query particle index
/// * `h`         – smoothing radius \[m\]
/// * `rho0`      – rest density \[kg/m³\]
pub fn pbf_grad_sum(particles: &[PbfExtParticle], i: usize, h: f64, rho0: f64) -> f64 {
    gradient_sum_over(particles, i, h, rho0, 0..particles.len())
}

/// Neighbour-list variant of [`pbf_grad_sum`]; identical when `nbrs` holds
/// every particle within `h` of particle `i`.
///
/// * `particles` – extended particle slice
/// * `i`         – query particle index
/// * `nbrs`      – neighbour indices of `i`
/// * `h`         – smoothing radius \[m\]
/// * `rho0`      – rest density \[kg/m³\]
pub fn pbf_grad_sum_nb(
    particles: &[PbfExtParticle],
    i: usize,
    nbrs: &[usize],
    h: f64,
    rho0: f64,
) -> f64 {
    gradient_sum_over(particles, i, h, rho0, nbrs.iter().copied())
}

/// PBF multiplier λ_i from the particle's neighbour list.
///
/// Evaluates `ρ_i` ([`pbf_density_nb`]), `Σ_k |∇_{p_k} C_i|²`
/// ([`pbf_grad_sum_nb`]) and combines them with [`pbf_lambda`]. For unit
/// masses it equals [`lambda_compute`] on the equivalent [`PbfParticle`] set.
///
/// * `particles` – extended particle slice
/// * `i`         – query particle index
/// * `nbrs`      – neighbour indices of `i`
/// * `h`         – smoothing radius \[m\]
/// * `rho0`      – rest density \[kg/m³\]
/// * `eps`       – regularisation constant (CFM)
pub fn pbf_lambda_nb(
    particles: &[PbfExtParticle],
    i: usize,
    nbrs: &[usize],
    h: f64,
    rho0: f64,
    eps: f64,
) -> f64 {
    let rho_i = pbf_density_nb(particles, i, nbrs, h);
    let grad_sum = pbf_grad_sum_nb(particles, i, nbrs, h, rho0);
    pbf_lambda(rho_i, rho0, grad_sum, eps)
}

/// PBF constraint multiplier λ for particle `i`.
///
/// λ = −(ρ_i/ρ₀ − 1) / (Σ_k |∇_{p_k} C_i|² + ε)
///
/// Returns `0` if the denominator is not positive (see [`pbf_grad_sum`] for
/// the gradient sum).
///
/// * `rho_i`    – density estimate at particle `i`
/// * `rho0`     – rest density \[kg/m³\]
/// * `grad_sum` – Σ_k |∇_{p_k} C_i|² \[1/m²\] (see [`pbf_grad_sum`])
/// * `eps`      – regularisation constant (CFM), see [`default_cfm_epsilon`]
pub fn pbf_lambda(rho_i: f64, rho0: f64, grad_sum: f64, eps: f64) -> f64 {
    constraint_multiplier(rho_i / rho0 - 1.0, grad_sum, eps)
}

/// Shared position-correction sum over the candidate neighbours `js`.
fn delta_pos_over<I: Iterator<Item = usize>>(
    particles: &[PbfExtParticle],
    i: usize,
    lambdas: &[f64],
    h: f64,
    rho0: f64,
    js: I,
) -> [f64; 3] {
    let mut dp = [0.0f64; 3];
    let (Some(pi), Some(&lambda_i)) = (particles.get(i), lambdas.get(i)) else {
        return dp;
    };
    if h.is_nan() || h <= 0.0 {
        return dp;
    }
    let h2 = h * h;
    for j in js {
        if j == i {
            continue;
        }
        let (Some(pj), Some(&lambda_j)) = (particles.get(j), lambdas.get(j)) else {
            continue;
        };
        let r_vec = vec3_sub(pi.pred_pos, pj.pred_pos);
        if vec3_norm_sq(r_vec) >= h2 {
            continue;
        }
        let g = spiky_gradient(r_vec, h);
        let factor = (pj.mass * lambda_i + pi.mass * lambda_j) / rho0;
        dp = vec3_add(dp, vec3_scale(g, factor));
    }
    dp
}

/// Position correction Δp_i for particle `i` from pre-computed lambdas
/// (all-pairs).
///
/// Δp_i = (1/ρ₀) Σ_j (m_j λ_i + m_i λ_j) ∇W_spiky(pred_i − pred_j, h)
///
/// For unit masses this is the familiar `(1/ρ₀) Σ_j (λ_i + λ_j) ∇W`. Scans
/// every particle; use [`pbf_delta_pos_nb`] inside solver loops.
///
/// * `particles` – extended particle slice (lambdas pre-computed)
/// * `i`         – query particle index
/// * `lambdas`   – lambda slice aligned with `particles`
/// * `h`         – smoothing radius \[m\]
/// * `rho0`      – rest density \[kg/m³\]
pub fn pbf_delta_pos(
    particles: &[PbfExtParticle],
    i: usize,
    lambdas: &[f64],
    h: f64,
    rho0: f64,
) -> [f64; 3] {
    delta_pos_over(particles, i, lambdas, h, rho0, 0..particles.len())
}

/// Position correction Δp_i summed over the neighbour list `nbrs`.
///
/// Identical to [`pbf_delta_pos`] provided `nbrs` contains every particle
/// within `h` of particle `i`.
///
/// * `particles` – extended particle slice (lambdas pre-computed)
/// * `i`         – query particle index
/// * `nbrs`      – neighbour indices of `i`
/// * `lambdas`   – lambda slice aligned with `particles`
/// * `h`         – smoothing radius \[m\]
/// * `rho0`      – rest density \[kg/m³\]
pub fn pbf_delta_pos_nb(
    particles: &[PbfExtParticle],
    i: usize,
    nbrs: &[usize],
    lambdas: &[f64],
    h: f64,
    rho0: f64,
) -> [f64; 3] {
    delta_pos_over(particles, i, lambdas, h, rho0, nbrs.iter().copied())
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

/// Relaxation coefficient κ of the default CFM parameter `ε = κ / h²`.
///
/// `Σ_k |∇_{p_k} C_i|²` grows like `(number of neighbours) / h²`, so a
/// relaxation proportional to `1/h²` keeps `ε / Σ|∇C|²` independent of the
/// length scale. `κ = 6` reproduces the `ε = 600` that Macklin & Müller (2013)
/// use for `h = 0.1`.
///
/// The relaxation is what keeps the Jacobi iteration stable: every particle is
/// moved by the constraints of all its ≈ 50 neighbours at once, and with
/// `ε ≪ Σ|∇C|²` (a full Newton step per constraint) free-surface particles are
/// thrown out of the fluid within a few steps.
pub const PBF_RELAXATION_COEFFICIENT: f64 = 6.0;

/// Default CFM regularisation ε \[1/m²\] for smoothing radius `h`:
/// `κ / h²` with `κ` = [`PBF_RELAXATION_COEFFICIENT`] (`0` for a non-positive
/// or non-finite `h`, for which no neighbour interaction exists anyway).
pub fn default_cfm_epsilon(h: f64) -> f64 {
    if h.is_finite() && h > 0.0 {
        PBF_RELAXATION_COEFFICIENT / (h * h)
    } else {
        0.0
    }
}

/// A self-contained Position-Based Fluids solver using `PbfExtParticle`.
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
    /// CFM regularisation ε \[1/m²\] added to Σ_k |∇C_i|² in
    /// `λ_i = −C_i / (Σ_k |∇C_i|² + ε)`.
    ///
    /// Larger values under-relax the projection; they bound the correction of
    /// sparsely populated neighbourhoods (free surface, spray) where the
    /// gradient sum is small, and keep the Jacobi iteration stable. Set to
    /// [`default_cfm_epsilon`]`(h)` (`6/h²`) by [`PbfSolver::new`]; recompute it
    /// if `h` is changed afterwards.
    pub epsilon: f64,
}

impl PbfSolver {
    /// Creates a new empty `PbfSolver`.
    pub fn new(rho0: f64, h: f64, substeps: usize) -> Self {
        Self {
            particles: Vec::new(),
            rho0,
            h,
            substeps,
            epsilon: default_cfm_epsilon(h),
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

    /// Neighbour lists of the current predicted positions (radius `h`).
    fn neighbor_lists(&self) -> NeighborCsr {
        NeighborCsr::build(
            self.particles.len(),
            |i| self.particles[i].pred_pos,
            self.h,
            false,
        )
    }

    /// Advances the simulation by one time step using `substeps` iterations.
    ///
    /// Predicts positions, builds the neighbour lists once (Algorithm 1 of
    /// Macklin & Müller 2013), runs `substeps` Jacobi iterations of the density
    /// constraint `C_i = ρ_i/ρ₀ − 1` with
    /// `λ_i = −C_i / (Σ_k |∇_{p_k} C_i|² + ε)` and
    /// `Δp_i = (1/ρ₀) Σ_j (m_j λ_i + m_i λ_j) ∇W_spiky`, and finally derives
    /// the velocities from the corrected positions. A non-positive or NaN `dt`
    /// leaves the particles untouched.
    ///
    /// * `dt`      – time step \[s\]
    /// * `gravity` – gravitational acceleration \[m/s²\]
    pub fn step(&mut self, dt: f64, gravity: [f64; 3]) {
        let n = self.particles.len();
        if n == 0 || dt.is_nan() || dt <= 0.0 {
            return;
        }

        // 1. Predict positions.
        for p in &mut self.particles {
            p.predict(gravity, dt);
        }

        // 2. Neighbour search, once per step.
        let neighbors = self.neighbor_lists();

        // 3. Constraint projection iterations.
        for _ in 0..self.substeps {
            let densities: Vec<f64> = (0..n)
                .map(|i| pbf_density_nb(&self.particles, i, neighbors.neighbors(i), self.h))
                .collect();

            let lambdas: Vec<f64> = (0..n)
                .map(|i| {
                    let grad_sum = pbf_grad_sum_nb(
                        &self.particles,
                        i,
                        neighbors.neighbors(i),
                        self.h,
                        self.rho0,
                    );
                    pbf_lambda(densities[i], self.rho0, grad_sum, self.epsilon)
                })
                .collect();

            let deltas: Vec<[f64; 3]> = (0..n)
                .map(|i| {
                    pbf_delta_pos_nb(
                        &self.particles,
                        i,
                        neighbors.neighbors(i),
                        &lambdas,
                        self.h,
                        self.rho0,
                    )
                })
                .collect();

            for ((p, &dp), (&lam, &rho)) in self
                .particles
                .iter_mut()
                .zip(deltas.iter())
                .zip(lambdas.iter().zip(densities.iter()))
            {
                p.pred_pos = vec3_add(p.pred_pos, dp);
                p.lambda = lam;
                p.rho = rho;
            }
        }

        // 4. Update velocities.
        for p in &mut self.particles {
            p.update_velocity(dt);
        }
    }

    /// Sum of |ρ_i/ρ₀ − 1| across all particles (density error).
    ///
    /// Densities are evaluated at the particles' predicted positions (equal to
    /// their positions after [`PbfSolver::step`]) with a spatial-hash
    /// neighbour search.
    pub fn total_density_error(&self) -> f64 {
        let neighbors = self.neighbor_lists();
        (0..self.particles.len())
            .map(|i| {
                let rho = pbf_density_nb(&self.particles, i, neighbors.neighbors(i), self.h);
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

#[cfg(test)]
mod nb_tests;
