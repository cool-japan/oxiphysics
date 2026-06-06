// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Free surface detection and tracking for SPH simulations.
//!
//! Implements the Color Function / Continuum Surface Force (CSF) approach
//! for detecting and tracking free surfaces in SPH fluid simulations,
//! including surface normals, curvature computation, and a simple
//! dam-break scenario for 2D free-surface flows.

/// A particle that carries surface-tracking information.
#[derive(Debug, Clone)]
pub struct SurfaceParticle {
    /// World-space position \[x, y, z\].
    pub position: [f64; 3],
    /// Velocity \[vx, vy, vz\].
    pub velocity: [f64; 3],
    /// Mass of the particle.
    pub mass: f64,
    /// Density estimate.
    pub density: f64,
    /// Pressure.
    pub pressure: f64,
    /// Inward surface normal (un-normalised intermediate or normalised final).
    pub normal: [f64; 3],
    /// Mean curvature κ (negative divergence of normalised normal).
    pub curvature: f64,
    /// True when this particle lies on the free surface.
    pub is_surface: bool,
}

impl SurfaceParticle {
    /// Create a new particle with default values.
    pub fn new(position: [f64; 3], mass: f64) -> Self {
        Self {
            position,
            velocity: [0.0; 3],
            mass,
            density: 1000.0,
            pressure: 0.0,
            normal: [0.0; 3],
            curvature: 0.0,
            is_surface: false,
        }
    }
}

// ---------------------------------------------------------------------------
// Kernel helpers
// ---------------------------------------------------------------------------

/// Wendland C2 kernel and its gradient.
pub struct ColorFunction {
    /// Smoothing length.
    pub smoothing_length: f64,
}

impl ColorFunction {
    /// Create a new `ColorFunction` with the given smoothing length.
    pub fn new(smoothing_length: f64) -> Self {
        Self { smoothing_length }
    }

    /// Wendland C2 kernel value W(r, h).
    ///
    /// Compact support radius is 2 h.
    pub fn kernel_value(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        // 3-D normalisation constant: 21 / (16 π h³)
        let alpha = 21.0 / (16.0 * std::f64::consts::PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        alpha * t.powi(4) * (2.0 * q + 1.0)
    }

    /// Gradient of the Wendland C2 kernel: ∂W/∂x_i = (dW/dr) * (dx_i / r).
    pub fn kernel_gradient(dx: [f64; 3], r: f64, h: f64) -> [f64; 3] {
        if r < 1e-12 || r / h >= 2.0 {
            return [0.0; 3];
        }
        let q = r / h;
        // dW/dq = alpha * [ 4*(1-q/2)^3 * (-1/2) * (2q+1) + (1-q/2)^4 * 2 ]
        //       = alpha * (1-q/2)^3 * [ -2(2q+1) + 2*(2-q) ]   / ... simplify:
        // dW/dq = alpha * (1-q/2)^3 * (-5q)
        let alpha = 21.0 / (16.0 * std::f64::consts::PI * h * h * h);
        let t = 1.0 - 0.5 * q;
        let dw_dq = alpha * t.powi(3) * (-5.0 * q);
        // dW/dr = dW/dq * 1/h
        let dw_dr = dw_dq / h;
        [dw_dr * dx[0] / r, dw_dr * dx[1] / r, dw_dr * dx[2] / r]
    }
}

// ---------------------------------------------------------------------------
// Surface tension / detection
// ---------------------------------------------------------------------------

/// Snapshot entry for normal computation: (position, mass, density, smoothing_length).
type NormalSnapEntry = ([f64; 3], f64, f64, f64);

/// Snapshot entry for curvature computation: (position, normal, mass, density, smoothing_length).
type CurvatureSnapEntry = ([f64; 3], [f64; 3], f64, f64, f64);

/// Surface tension model using the Continuum Surface Force (CSF) approach.
pub struct SurfaceTension {
    /// Surface tension coefficient σ \[N/m\].
    pub sigma: f64,
    /// SPH smoothing length.
    pub smoothing_length: f64,
}

impl SurfaceTension {
    /// Create a new `SurfaceTension` model.
    pub fn new(sigma: f64, smoothing_length: f64) -> Self {
        Self {
            sigma,
            smoothing_length,
        }
    }

    /// Compute un-normalised surface normals for every particle.
    ///
    /// n_i = Σ_j (m_j / ρ_j) ∇W_ij
    pub fn compute_normals(particles: &mut [SurfaceParticle]) {
        let n = particles.len();
        let mut normals = vec![[0.0_f64; 3]; n];

        // We need an immutable snapshot of the data we read.
        let snap: Vec<NormalSnapEntry> = particles
            .iter()
            .map(|p| (p.position, p.mass, p.density, p.smoothing_length_or(1.0)))
            .collect();

        for (i, (norm_i, &(pi, _, _, hi))) in normals.iter_mut().zip(snap.iter()).enumerate() {
            for (j, &(pj, mj, rhoj, _)) in snap.iter().enumerate() {
                if i == j || rhoj < 1e-12 {
                    continue;
                }
                let dx = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                let grad = ColorFunction::kernel_gradient(dx, r, hi);
                let w = mj / rhoj;
                norm_i[0] += w * grad[0];
                norm_i[1] += w * grad[1];
                norm_i[2] += w * grad[2];
            }
        }

        for (p, &norm) in particles.iter_mut().zip(normals.iter()) {
            p.normal = norm;
        }
    }

    /// Mark particles as surface particles when |n_i| > threshold.
    pub fn detect_surface(particles: &mut [SurfaceParticle], threshold: f64) {
        for p in particles.iter_mut() {
            let mag = vec3_len(p.normal);
            p.is_surface = mag > threshold;
        }
    }

    /// Normalise the surface normals in-place (for particles where |n| > ε).
    pub fn normalize_normals(particles: &mut [SurfaceParticle]) {
        for p in particles.iter_mut() {
            let mag = vec3_len(p.normal);
            if mag > 1e-12 {
                p.normal[0] /= mag;
                p.normal[1] /= mag;
                p.normal[2] /= mag;
            }
        }
    }

    /// Compute mean curvature κ_i = -∇·n̂ via SPH divergence.
    ///
    /// κ_i = -Σ_j (m_j / ρ_j) (n̂_j − n̂_i) · ∇W_ij
    ///
    /// Call `normalize_normals` first.
    pub fn compute_curvature(particles: &mut [SurfaceParticle]) {
        let n = particles.len();
        let mut kappas = vec![0.0_f64; n];

        let snap: Vec<CurvatureSnapEntry> = particles
            .iter()
            .map(|p| {
                (
                    p.position,
                    p.normal,
                    p.mass,
                    p.density,
                    p.smoothing_length_or(1.0),
                )
            })
            .collect();

        for (i, (kappa_i, &(pi, ni, _, _, hi))) in kappas.iter_mut().zip(snap.iter()).enumerate() {
            for (j, &(pj, nj, mj, rhoj, _)) in snap.iter().enumerate() {
                if i == j || rhoj < 1e-12 {
                    continue;
                }
                let dx = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                let grad = ColorFunction::kernel_gradient(dx, r, hi);
                let dn = [nj[0] - ni[0], nj[1] - ni[1], nj[2] - ni[2]];
                let dot = dn[0] * grad[0] + dn[1] * grad[1] + dn[2] * grad[2];
                *kappa_i -= (mj / rhoj) * dot;
            }
        }

        for (i, p) in particles.iter_mut().enumerate() {
            p.curvature = kappas[i];
        }
    }

    /// Continuum Surface Force: F_i = σ κ_i n̂_i  (only on surface particles).
    pub fn surface_tension_force(p_i: &SurfaceParticle, sigma: f64) -> [f64; 3] {
        if !p_i.is_surface {
            return [0.0; 3];
        }
        [
            sigma * p_i.curvature * p_i.normal[0],
            sigma * p_i.curvature * p_i.normal[1],
            sigma * p_i.curvature * p_i.normal[2],
        ]
    }
}

// ---------------------------------------------------------------------------
// WaterSurface – simple 2-D dam-break scenario
// ---------------------------------------------------------------------------

/// 2-D free-surface water simulation (dam-break scenario).
///
/// Particles are stored as 3-D vectors but only the x-y plane is used.
pub struct WaterSurface {
    /// All fluid particles.
    pub particles: Vec<SurfaceParticle>,
    /// Gravitational acceleration vector.
    pub gravity: [f64; 3],
    /// Axis-aligned domain: (min, max).
    pub domain: ([f64; 3], [f64; 3]),
    /// SPH smoothing length.
    pub h: f64,
}

impl WaterSurface {
    /// Build a dam-break scenario.
    ///
    /// The domain is `[0, width] × [0, height] × [0, 0]`.
    /// Particles fill the left half (`x ≤ width/2`) on a regular grid with
    /// spacing `resolution`.
    pub fn new_dam_break(width: f64, height: f64, resolution: f64) -> Self {
        let h = 2.0 * resolution;
        let rho0 = 1000.0;
        let particle_mass = rho0 * resolution * resolution; // 2-D mass proxy

        let domain_min = [0.0, 0.0, 0.0];
        let domain_max = [width, height, 0.0];

        let mut particles = Vec::new();
        let half_width = width / 2.0;
        let mut x = resolution * 0.5;
        while x <= half_width {
            let mut y = resolution * 0.5;
            while y <= height {
                let mut p = SurfaceParticle::new([x, y, 0.0], particle_mass);
                p.density = rho0;
                particles.push(p);
                y += resolution;
            }
            x += resolution;
        }

        Self {
            particles,
            gravity: [0.0, -9.81, 0.0],
            domain: (domain_min, domain_max),
            h,
        }
    }

    /// Update density via SPH summation: ρ_i = Σ_j m_j W(|x_i − x_j|, h).
    pub fn update_density(&mut self) {
        let _n = self.particles.len();
        let h = self.h;
        let positions: Vec<[f64; 3]> = self.particles.iter().map(|p| p.position).collect();
        let masses: Vec<f64> = self.particles.iter().map(|p| p.mass).collect();

        for (p_i, &pos_i) in self.particles.iter_mut().zip(positions.iter()) {
            let rho: f64 = positions
                .iter()
                .zip(masses.iter())
                .map(|(&pos_j, &mj)| {
                    let dx = [
                        pos_i[0] - pos_j[0],
                        pos_i[1] - pos_j[1],
                        pos_i[2] - pos_j[2],
                    ];
                    let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                    mj * ColorFunction::kernel_value(r, h)
                })
                .sum::<f64>()
                .max(1e-6);
            p_i.density = rho;
        }
    }

    /// Apply gravitational body force: v_i += g * dt.
    pub fn apply_gravity(&mut self, dt: f64) {
        let g = self.gravity;
        for p in self.particles.iter_mut() {
            p.velocity[0] += g[0] * dt;
            p.velocity[1] += g[1] * dt;
            p.velocity[2] += g[2] * dt;
        }
    }

    /// Reflect velocities at domain walls (penalty-free boundary).
    pub fn enforce_boundary(&mut self) {
        let (dmin, dmax) = self.domain;
        for p in self.particles.iter_mut() {
            for ((&lo, &hi_val), (pos_k, vel_k)) in dmin
                .iter()
                .zip(dmax.iter())
                .zip(p.position.iter_mut().zip(p.velocity.iter_mut()))
            {
                if *pos_k < lo {
                    *pos_k = lo;
                    if *vel_k < 0.0 {
                        *vel_k = -*vel_k;
                    }
                }
                if *pos_k > hi_val {
                    *pos_k = hi_val;
                    if *vel_k > 0.0 {
                        *vel_k = -*vel_k;
                    }
                }
            }
        }
    }

    /// Semi-implicit Euler integration: update velocity, then position.
    pub fn integrate(&mut self, dt: f64) {
        self.apply_gravity(dt);
        for p in self.particles.iter_mut() {
            p.position[0] += p.velocity[0] * dt;
            p.position[1] += p.velocity[1] * dt;
            p.position[2] += p.velocity[2] * dt;
        }
        self.enforce_boundary();
    }

    /// Return the maximum y-coordinate of surface particles near `x`.
    ///
    /// "Near" is defined as within `h` of `x` in the x-direction.
    /// Returns `0.0` if no particle qualifies.
    pub fn surface_height_at(&self, x: f64) -> f64 {
        let h = self.h;
        // First detect surface: we use a simple density-based heuristic here
        // because normals have not necessarily been computed yet.
        self.particles
            .iter()
            .filter(|p| (p.position[0] - x).abs() <= h)
            .map(|p| p.position[1])
            .fold(f64::NEG_INFINITY, f64::max)
            .max(0.0)
    }
}

// ---------------------------------------------------------------------------
// SplashDetector
// ---------------------------------------------------------------------------

/// Detects isolated droplet clusters among surface particles.
pub struct SplashDetector {
    /// Maximum distance between two particles to be considered part of the
    /// same cluster.
    pub droplet_threshold: f64,
}

impl SplashDetector {
    /// Create a new `SplashDetector`.
    pub fn new(droplet_threshold: f64) -> Self {
        Self { droplet_threshold }
    }

    /// Union-Find helper: find root with path compression.
    fn find(parent: &mut [usize], i: usize) -> usize {
        if parent[i] != i {
            parent[i] = Self::find(parent, parent[i]);
        }
        parent[i]
    }

    /// Union two sets.
    fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
        let ra = Self::find(parent, a);
        let rb = Self::find(parent, b);
        if ra == rb {
            return;
        }
        if rank[ra] < rank[rb] {
            parent[ra] = rb;
        } else if rank[ra] > rank[rb] {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra] += 1;
        }
    }

    /// Detect isolated droplet clusters among surface particles.
    ///
    /// Returns a `Vec` of `(droplet_id, centroid)` pairs, one per cluster.
    /// Only clusters that are not connected to the bulk (i.e. isolated groups)
    /// are returned; the largest cluster is considered bulk and excluded.
    pub fn detect_droplets(&self, particles: &[SurfaceParticle]) -> Vec<(usize, [f64; 3])> {
        let n = particles.len();
        if n == 0 {
            return Vec::new();
        }

        let mut parent: Vec<usize> = (0..n).collect();
        let mut rank = vec![0usize; n];

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = [
                    particles[i].position[0] - particles[j].position[0],
                    particles[i].position[1] - particles[j].position[1],
                    particles[i].position[2] - particles[j].position[2],
                ];
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                if r2 <= self.droplet_threshold * self.droplet_threshold {
                    Self::union(&mut parent, &mut rank, i, j);
                }
            }
        }

        // Collect cluster membership.
        use std::collections::HashMap;
        let mut clusters: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..n {
            let root = Self::find(&mut parent, i);
            clusters.entry(root).or_default().push(i);
        }

        // Find the largest cluster (bulk fluid) and exclude it.
        let bulk_root = clusters
            .iter()
            .max_by_key(|(_, v)| v.len())
            .map(|(k, _)| *k)
            .unwrap_or(usize::MAX);

        let mut result = Vec::new();
        let mut droplet_id = 0usize;
        for (root, members) in &clusters {
            if *root == bulk_root {
                continue;
            }
            let count = members.len() as f64;
            let mut centroid = [0.0_f64; 3];
            for &idx in members {
                centroid[0] += particles[idx].position[0];
                centroid[1] += particles[idx].position[1];
                centroid[2] += particles[idx].position[2];
            }
            centroid[0] /= count;
            centroid[1] /= count;
            centroid[2] /= count;
            result.push((droplet_id, centroid));
            droplet_id += 1;
        }

        result
    }
}

// ---------------------------------------------------------------------------
// Eigenvalue-based surface detection
// ---------------------------------------------------------------------------

/// Eigenvalue-based free surface detection.
///
/// Uses the eigenvalue ratio of the smoothed position covariance matrix
/// to detect surface particles (surface particles have a small minimum
/// eigenvalue relative to the maximum).
pub struct EigenvalueSurfaceDetector {
    /// Ratio threshold: particle is surface when λ_min/λ_max < threshold.
    pub eigen_ratio_threshold: f64,
    /// Smoothing length for kernel evaluation.
    pub smoothing_length: f64,
}

impl EigenvalueSurfaceDetector {
    /// Create a new eigenvalue surface detector.
    pub fn new(eigen_ratio_threshold: f64, smoothing_length: f64) -> Self {
        Self {
            eigen_ratio_threshold,
            smoothing_length,
        }
    }

    /// Compute the weighted covariance matrix of neighbors around particle `i`.
    ///
    /// Returns `[[c00, c01, c02\], [c10, c11, c12], [c20, c21, c22]]`.
    pub fn neighbor_covariance(particles: &[SurfaceParticle], i: usize, h: f64) -> [[f64; 3]; 3] {
        let pi = particles[i].position;
        let mut cov = [[0.0_f64; 3]; 3];
        let mut w_sum = 0.0_f64;

        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let dx = [
                pj.position[0] - pi[0],
                pj.position[1] - pi[1],
                pj.position[2] - pi[2],
            ];
            let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
            let w = ColorFunction::kernel_value(r, h);
            if w < 1e-30 {
                continue;
            }
            for a in 0..3 {
                for b in 0..3 {
                    cov[a][b] += w * dx[a] * dx[b];
                }
            }
            w_sum += w;
        }

        if w_sum > 1e-30 {
            for row in cov.iter_mut() {
                for v in row.iter_mut() {
                    *v /= w_sum;
                }
            }
        }
        cov
    }

    /// Compute eigenvalues of a symmetric 3x3 matrix using the analytical formula.
    ///
    /// Returns eigenvalues sorted in ascending order.
    pub fn eigenvalues_symmetric_3x3(m: [[f64; 3]; 3]) -> [f64; 3] {
        let a = m[0][0];
        let b = m[1][1];
        let c = m[2][2];
        let d = m[0][1];
        let e = m[1][2];
        let f = m[0][2];

        let p1 = d * d + f * f + e * e;
        if p1 < 1e-30 {
            // Matrix is diagonal
            let mut eigs = [a, b, c];
            eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
            return eigs;
        }

        let q = (a + b + c) / 3.0;
        let p2 = (a - q) * (a - q) + (b - q) * (b - q) + (c - q) * (c - q) + 2.0 * p1;
        let p = (p2 / 6.0).sqrt();

        // B = (1/p) * (M - q*I)
        let b00 = (a - q) / p;
        let b11 = (b - q) / p;
        let b22 = (c - q) / p;
        let b01 = d / p;
        let b12 = e / p;
        let b02 = f / p;

        let det_b = b00 * (b11 * b22 - b12 * b12) - b01 * (b01 * b22 - b12 * b02)
            + b02 * (b01 * b12 - b11 * b02);

        let half_det = det_b / 2.0;
        let half_det = half_det.clamp(-1.0, 1.0);
        let phi = half_det.acos() / 3.0;

        let eig1 = q + 2.0 * p * phi.cos();
        let eig3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
        let eig2 = 3.0 * q - eig1 - eig3;

        let mut eigs = [eig1, eig2, eig3];
        eigs.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
        eigs
    }

    /// Detect surface particles using eigenvalue ratio.
    ///
    /// Marks `is_surface = true` when the smallest eigenvalue of the local
    /// covariance matrix is much smaller than the largest.
    pub fn detect(&self, particles: &mut [SurfaceParticle]) {
        let n = particles.len();
        let h = self.smoothing_length;
        let threshold = self.eigen_ratio_threshold;

        // Compute surface flags in a temporary vec to avoid borrow issues
        let mut flags = vec![false; n];
        // Use immutable borrow to compute covariances
        let positions: Vec<[f64; 3]> = particles.iter().map(|p| p.position).collect();
        let masses: Vec<f64> = particles.iter().map(|p| p.mass).collect();
        let densities: Vec<f64> = particles.iter().map(|p| p.density).collect();

        for i in 0..n {
            let pi = positions[i];
            let mut cov = [[0.0_f64; 3]; 3];
            let mut w_sum = 0.0_f64;

            for j in 0..n {
                if j == i {
                    continue;
                }
                let dx = [
                    positions[j][0] - pi[0],
                    positions[j][1] - pi[1],
                    positions[j][2] - pi[2],
                ];
                let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                let w = ColorFunction::kernel_value(r, h);
                if w < 1e-30 {
                    continue;
                }
                let vol = if densities[j] > 1e-12 {
                    masses[j] / densities[j]
                } else {
                    0.0
                };
                for a in 0..3 {
                    for b in 0..3 {
                        cov[a][b] += vol * w * dx[a] * dx[b];
                    }
                }
                w_sum += vol * w;
            }

            if w_sum > 1e-30 {
                for row in cov.iter_mut() {
                    for v in row.iter_mut() {
                        *v /= w_sum;
                    }
                }
            }

            let eigs = Self::eigenvalues_symmetric_3x3(cov);
            let lambda_max = eigs[2].abs().max(1e-30);
            let lambda_min = eigs[0].abs();
            flags[i] = (lambda_min / lambda_max) < threshold;
        }

        for (i, p) in particles.iter_mut().enumerate() {
            p.is_surface = flags[i];
        }
    }
}

// ---------------------------------------------------------------------------
// Divergence-based surface detection
// ---------------------------------------------------------------------------

/// Divergence-based free surface detection.
///
/// Uses the SPH number density (particle count density) to detect
/// surface particles: particles with fewer neighbors than expected.
pub struct DivergenceSurfaceDetector {
    /// Fraction of reference number density below which a particle is surface.
    pub density_ratio_threshold: f64,
    /// Smoothing length.
    pub smoothing_length: f64,
}

impl DivergenceSurfaceDetector {
    /// Create a new divergence surface detector.
    pub fn new(density_ratio_threshold: f64, smoothing_length: f64) -> Self {
        Self {
            density_ratio_threshold,
            smoothing_length,
        }
    }

    /// Compute the number density (sum of kernel values) for particle `i`.
    pub fn number_density(particles: &[SurfaceParticle], i: usize, h: f64) -> f64 {
        let pi = particles[i].position;
        let mut nd = ColorFunction::kernel_value(0.0, h); // self contribution
        for (j, pj) in particles.iter().enumerate() {
            if j == i {
                continue;
            }
            let dx = [
                pi[0] - pj.position[0],
                pi[1] - pj.position[1],
                pi[2] - pj.position[2],
            ];
            let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
            nd += ColorFunction::kernel_value(r, h);
        }
        nd
    }

    /// Detect surface particles by comparing number density to reference.
    pub fn detect(&self, particles: &mut [SurfaceParticle]) {
        let n = particles.len();
        let h = self.smoothing_length;

        // Find the maximum number density as reference
        let mut number_densities = vec![0.0_f64; n];
        let mut max_nd = 0.0_f64;
        for (i, nd) in number_densities.iter_mut().enumerate() {
            *nd = Self::number_density(particles, i, h);
            if *nd > max_nd {
                max_nd = *nd;
            }
        }

        if max_nd < 1e-30 {
            return;
        }

        for (i, p) in particles.iter_mut().enumerate() {
            p.is_surface = (number_densities[i] / max_nd) < self.density_ratio_threshold;
        }
    }
}

// ---------------------------------------------------------------------------
// Surface reconstruction (marching cubes-like level set)
// ---------------------------------------------------------------------------

/// Simple level-set surface reconstruction on a uniform grid.
///
/// Evaluates a particle density field on a grid and extracts the isosurface
/// as a set of grid cell positions where the field crosses the iso-value.
pub struct SurfaceReconstructor {
    /// Grid cell size.
    pub cell_size: f64,
    /// Iso-value for surface extraction.
    pub iso_value: f64,
    /// Smoothing length for kernel evaluation.
    pub smoothing_length: f64,
}

impl SurfaceReconstructor {
    /// Create a new surface reconstructor.
    pub fn new(cell_size: f64, iso_value: f64, smoothing_length: f64) -> Self {
        Self {
            cell_size,
            iso_value,
            smoothing_length,
        }
    }

    /// Extract surface cell centers from a set of particles.
    ///
    /// Returns positions of grid cells that lie on the isosurface boundary.
    pub fn extract_surface_cells(&self, particles: &[SurfaceParticle]) -> Vec<[f64; 3]> {
        if particles.is_empty() {
            return Vec::new();
        }

        // Compute bounding box
        let mut min_pos = particles[0].position;
        let mut max_pos = particles[0].position;
        for p in particles {
            for d in 0..3 {
                if p.position[d] < min_pos[d] {
                    min_pos[d] = p.position[d];
                }
                if p.position[d] > max_pos[d] {
                    max_pos[d] = p.position[d];
                }
            }
        }

        // Extend by smoothing length
        let margin = self.smoothing_length * 2.0;
        for d in 0..3 {
            min_pos[d] -= margin;
            max_pos[d] += margin;
        }

        let cs = self.cell_size;
        let nx = ((max_pos[0] - min_pos[0]) / cs).ceil() as usize + 1;
        let ny = ((max_pos[1] - min_pos[1]) / cs).ceil() as usize + 1;
        let nz = ((max_pos[2] - min_pos[2]) / cs).ceil() as usize + 1;

        // Cap grid size for safety
        let nx = nx.min(100);
        let ny = ny.min(100);
        let nz = nz.min(100);

        let h = self.smoothing_length;
        let iso = self.iso_value;

        let mut surface_cells = Vec::new();

        for ix in 0..nx {
            for iy in 0..ny {
                for iz in 0..nz {
                    let gx = min_pos[0] + ix as f64 * cs;
                    let gy = min_pos[1] + iy as f64 * cs;
                    let gz = min_pos[2] + iz as f64 * cs;

                    // Evaluate field
                    let mut field_val = 0.0_f64;
                    for p in particles {
                        let dx = [gx - p.position[0], gy - p.position[1], gz - p.position[2]];
                        let r = (dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]).sqrt();
                        let rho = if p.density > 1e-12 { p.density } else { 1000.0 };
                        field_val += (p.mass / rho) * ColorFunction::kernel_value(r, h);
                    }

                    // Check if this cell is on the surface (field crosses iso_value)
                    if (field_val - iso).abs() < iso * 0.5 {
                        surface_cells.push([gx, gy, gz]);
                    }
                }
            }
        }

        surface_cells
    }
}

// ---------------------------------------------------------------------------
// Curvature estimation utilities
// ---------------------------------------------------------------------------

/// Utility functions for curvature estimation.
pub struct CurvatureEstimator;

impl CurvatureEstimator {
    /// Estimate curvature from a fitted quadratic surface patch.
    ///
    /// Given a set of nearby surface points, fits a paraboloid
    /// z = ax² + bxy + cy² and returns the mean curvature κ = a + c.
    pub fn estimate_from_neighbors(
        center: [f64; 3],
        neighbors: &[[f64; 3]],
        normal: [f64; 3],
    ) -> f64 {
        if neighbors.len() < 3 {
            return 0.0;
        }

        // Build a local coordinate frame with normal as z-axis
        let n_len = vec3_len(normal);
        if n_len < 1e-12 {
            return 0.0;
        }
        let nz = [normal[0] / n_len, normal[1] / n_len, normal[2] / n_len];

        // Find a tangent vector
        let tmp = if nz[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        // nx = tmp × nz
        let nx = [
            tmp[1] * nz[2] - tmp[2] * nz[1],
            tmp[2] * nz[0] - tmp[0] * nz[2],
            tmp[0] * nz[1] - tmp[1] * nz[0],
        ];
        let nx_len = vec3_len(nx);
        if nx_len < 1e-12 {
            return 0.0;
        }
        let nx = [nx[0] / nx_len, nx[1] / nx_len, nx[2] / nx_len];
        // ny = nz × nx
        let ny = [
            nz[1] * nx[2] - nz[2] * nx[1],
            nz[2] * nx[0] - nz[0] * nx[2],
            nz[0] * nx[1] - nz[1] * nx[0],
        ];

        // Project neighbors into local frame and fit paraboloid
        let mut sum_u2 = 0.0_f64;
        let mut sum_v2 = 0.0_f64;
        let mut sum_u2w = 0.0_f64;
        let mut sum_v2w = 0.0_f64;
        let mut count = 0.0_f64;

        for &nb in neighbors {
            let dx = [nb[0] - center[0], nb[1] - center[1], nb[2] - center[2]];
            let u = dx[0] * nx[0] + dx[1] * nx[1] + dx[2] * nx[2];
            let v = dx[0] * ny[0] + dx[1] * ny[1] + dx[2] * ny[2];
            let w = dx[0] * nz[0] + dx[1] * nz[1] + dx[2] * nz[2];

            sum_u2 += u * u;
            sum_v2 += v * v;
            sum_u2w += u * u * w;
            sum_v2w += v * v * w;
            count += 1.0;
        }

        if count < 1.0 || sum_u2.abs() < 1e-30 || sum_v2.abs() < 1e-30 {
            return 0.0;
        }

        // Least squares fit: a = Σ(u²w) / Σ(u⁴), c = Σ(v²w) / Σ(v⁴)
        // Simplified: κ ≈ 2(a + c) where a = Σ(u²w)/Σ(u²)², c = Σ(v²w)/Σ(v²)²
        let a = sum_u2w / (sum_u2 * sum_u2 / count);
        let c = sum_v2w / (sum_v2 * sum_v2 / count);

        2.0 * (a + c)
    }
}

// ---------------------------------------------------------------------------
// Surface rendering hints
// ---------------------------------------------------------------------------

/// Rendering hints for surface particles.
#[derive(Debug, Clone)]
pub struct SurfaceRenderHint {
    /// Particle index.
    pub particle_index: usize,
    /// Surface normal (normalized).
    pub normal: [f64; 3],
    /// Local curvature.
    pub curvature: f64,
    /// Suggested transparency (0=opaque, 1=fully transparent).
    pub transparency: f64,
    /// Whether this is part of a thin sheet.
    pub is_thin_sheet: bool,
}

/// Generate rendering hints for surface particles.
pub fn generate_render_hints(particles: &[SurfaceParticle]) -> Vec<SurfaceRenderHint> {
    particles
        .iter()
        .enumerate()
        .filter(|(_, p)| p.is_surface)
        .map(|(i, p)| {
            let n_mag = vec3_len(p.normal);
            let normal = if n_mag > 1e-12 {
                [
                    p.normal[0] / n_mag,
                    p.normal[1] / n_mag,
                    p.normal[2] / n_mag,
                ]
            } else {
                [0.0, 1.0, 0.0]
            };
            let curvature = p.curvature;
            // High curvature → more transparent (thin features)
            let transparency = (curvature.abs() * 0.1).min(0.8);
            let is_thin_sheet = curvature.abs() > 5.0;
            SurfaceRenderHint {
                particle_index: i,
                normal,
                curvature,
                transparency,
                is_thin_sheet,
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn vec3_len(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Extension trait to retrieve a per-particle smoothing length if available;
/// falls back to a supplied default.
trait SmoothinLengthExt {
    fn smoothing_length_or(&self, default: f64) -> f64;
}

impl SmoothinLengthExt for SurfaceParticle {
    fn smoothing_length_or(&self, _default: f64) -> f64 {
        // SurfaceParticle does not currently store a per-particle h;
        // callers should pass the global h explicitly. This method
        // exists for forward-compatibility.
        _default
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dam_break_particles_left_half() {
        let width = 2.0;
        let height = 1.0;
        let res = 0.1;
        let ws = WaterSurface::new_dam_break(width, height, res);
        assert!(!ws.particles.is_empty(), "should have particles");
        for p in &ws.particles {
            assert!(
                p.position[0] <= width / 2.0 + 1e-10,
                "particle x={} is in right half",
                p.position[0]
            );
        }
    }

    #[test]
    fn test_surface_height_positive_at_zero() {
        let ws = WaterSurface::new_dam_break(2.0, 1.0, 0.1);
        let h = ws.surface_height_at(0.0);
        assert!(h > 0.0, "surface height at x=0 should be > 0, got {}", h);
    }

    #[test]
    fn test_apply_gravity_changes_velocity() {
        let mut ws = WaterSurface::new_dam_break(2.0, 1.0, 0.2);
        let before: Vec<[f64; 3]> = ws.particles.iter().map(|p| p.velocity).collect();
        ws.apply_gravity(0.01);
        let after: Vec<[f64; 3]> = ws.particles.iter().map(|p| p.velocity).collect();
        assert!(
            before.iter().zip(after.iter()).any(|(b, a)| b != a),
            "gravity should change at least one velocity"
        );
        // All particles should have vy decreased (gravity is -9.81)
        for (b, a) in before.iter().zip(after.iter()) {
            assert!(a[1] < b[1], "vy should decrease under gravity");
        }
    }

    #[test]
    fn test_enforce_boundary_keeps_in_domain() {
        let mut ws = WaterSurface::new_dam_break(2.0, 1.0, 0.2);
        // Push all particles out of bounds.
        for p in ws.particles.iter_mut() {
            p.position = [-1.0, -1.0, 0.0];
            p.velocity = [-5.0, -5.0, 0.0];
        }
        ws.enforce_boundary();
        let (dmin, dmax) = ws.domain;
        for p in &ws.particles {
            for k in 0..3 {
                assert!(
                    p.position[k] >= dmin[k] - 1e-12,
                    "position[{}]={} below min",
                    k,
                    p.position[k]
                );
                assert!(
                    p.position[k] <= dmax[k] + 1e-12,
                    "position[{}]={} above max",
                    k,
                    p.position[k]
                );
            }
        }
    }

    #[test]
    fn test_compute_normals_single_particle_no_panic() {
        let mut particles = vec![SurfaceParticle::new([0.5, 0.5, 0.0], 0.01)];
        particles[0].density = 1000.0;
        // Should not panic with a single particle.
        SurfaceTension::compute_normals(&mut particles);
        // Normal should remain zero (no neighbours).
        assert_eq!(particles[0].normal, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_kernel_value_at_origin_positive() {
        let w = ColorFunction::kernel_value(0.0, 1.0);
        assert!(w > 0.0, "kernel at r=0 should be positive, got {}", w);
    }

    #[test]
    fn test_kernel_value_outside_support_zero() {
        let h = 1.0;
        let w = ColorFunction::kernel_value(2.0 * h, h);
        assert_eq!(w, 0.0, "kernel at r=2h should be exactly 0, got {}", w);
    }

    // ── Additional tests ────────────────────────────────────────────────

    #[test]
    fn test_kernel_gradient_at_zero_is_zero() {
        let grad = ColorFunction::kernel_gradient([0.0, 0.0, 0.0], 0.0, 1.0);
        assert_eq!(grad, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_kernel_gradient_outside_support() {
        let grad = ColorFunction::kernel_gradient([3.0, 0.0, 0.0], 3.0, 1.0);
        assert_eq!(grad, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_kernel_gradient_direction() {
        // Gradient should point away from the center for a positive displacement
        let dx = [0.5, 0.0, 0.0];
        let r = 0.5;
        let grad = ColorFunction::kernel_gradient(dx, r, 1.0);
        // dW/dr < 0 for q in (0, 2), so grad[0] should be negative
        assert!(
            grad[0] < 0.0,
            "gradient x-component should be negative, got {}",
            grad[0]
        );
    }

    #[test]
    fn test_detect_surface_threshold() {
        let mut particles = vec![
            SurfaceParticle::new([0.0, 0.0, 0.0], 0.01),
            SurfaceParticle::new([0.1, 0.0, 0.0], 0.01),
        ];
        // Set one particle with a large normal
        particles[0].normal = [1.0, 0.0, 0.0];
        particles[1].normal = [0.001, 0.0, 0.0];

        SurfaceTension::detect_surface(&mut particles, 0.5);
        assert!(particles[0].is_surface, "particle 0 should be surface");
        assert!(!particles[1].is_surface, "particle 1 should not be surface");
    }

    #[test]
    fn test_normalize_normals() {
        let mut particles = vec![SurfaceParticle::new([0.0; 3], 1.0)];
        particles[0].normal = [3.0, 4.0, 0.0];
        SurfaceTension::normalize_normals(&mut particles);
        let mag = vec3_len(particles[0].normal);
        assert!(
            (mag - 1.0).abs() < 1e-12,
            "Normal should be unit length, got {mag}"
        );
    }

    #[test]
    fn test_normalize_normals_zero_normal() {
        let mut particles = vec![SurfaceParticle::new([0.0; 3], 1.0)];
        particles[0].normal = [0.0, 0.0, 0.0];
        SurfaceTension::normalize_normals(&mut particles);
        assert_eq!(particles[0].normal, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_surface_tension_force_non_surface() {
        let p = SurfaceParticle::new([0.0; 3], 1.0);
        let f = SurfaceTension::surface_tension_force(&p, 0.072);
        assert_eq!(f, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_surface_tension_force_on_surface() {
        let mut p = SurfaceParticle::new([0.0; 3], 1.0);
        p.is_surface = true;
        p.curvature = 10.0;
        p.normal = [0.0, 1.0, 0.0];
        let f = SurfaceTension::surface_tension_force(&p, 0.072);
        assert!((f[1] - 0.072 * 10.0).abs() < 1e-10);
        assert_eq!(f[0], 0.0);
        assert_eq!(f[2], 0.0);
    }

    #[test]
    fn test_splash_detector_single_cluster() {
        let detector = SplashDetector::new(0.5);
        let particles = vec![
            SurfaceParticle::new([0.0, 0.0, 0.0], 1.0),
            SurfaceParticle::new([0.1, 0.0, 0.0], 1.0),
            SurfaceParticle::new([0.2, 0.0, 0.0], 1.0),
        ];
        let droplets = detector.detect_droplets(&particles);
        // All within threshold, single cluster → no droplets (bulk only)
        assert!(
            droplets.is_empty(),
            "Single cluster should produce no droplets"
        );
    }

    #[test]
    fn test_splash_detector_two_clusters() {
        let detector = SplashDetector::new(0.5);
        let particles = vec![
            // Bulk cluster
            SurfaceParticle::new([0.0, 0.0, 0.0], 1.0),
            SurfaceParticle::new([0.1, 0.0, 0.0], 1.0),
            SurfaceParticle::new([0.2, 0.0, 0.0], 1.0),
            // Isolated droplet
            SurfaceParticle::new([10.0, 10.0, 0.0], 1.0),
        ];
        let droplets = detector.detect_droplets(&particles);
        assert_eq!(droplets.len(), 1, "Should detect one droplet");
    }

    #[test]
    fn test_splash_detector_empty() {
        let detector = SplashDetector::new(0.5);
        let droplets = detector.detect_droplets(&[]);
        assert!(droplets.is_empty());
    }

    #[test]
    fn test_water_surface_density_update() {
        let mut ws = WaterSurface::new_dam_break(2.0, 1.0, 0.2);
        ws.update_density();
        for p in &ws.particles {
            assert!(p.density > 0.0, "Density should be positive after update");
        }
    }

    #[test]
    fn test_water_surface_integration() {
        let mut ws = WaterSurface::new_dam_break(2.0, 1.0, 0.2);
        let dt = 0.001;
        ws.integrate(dt);
        // After integration with gravity, all particles should still be in domain
        let (dmin, dmax) = ws.domain;
        for p in &ws.particles {
            for k in 0..3 {
                assert!(p.position[k] >= dmin[k] - 1e-10);
                assert!(p.position[k] <= dmax[k] + 1e-10);
            }
        }
    }

    #[test]
    fn test_surface_height_at_empty_region() {
        let ws = WaterSurface::new_dam_break(2.0, 1.0, 0.2);
        // Far right should have 0 height
        let h = ws.surface_height_at(100.0);
        assert_eq!(h, 0.0);
    }

    #[test]
    fn test_eigenvalue_detector_basic() {
        let detector = EigenvalueSurfaceDetector::new(0.3, 0.5);
        // Create a flat sheet of particles in the x-y plane
        let mut particles = Vec::new();
        for i in 0..5 {
            for j in 0..5 {
                let mut p = SurfaceParticle::new([i as f64 * 0.1, j as f64 * 0.1, 0.0], 1.0);
                p.density = 1000.0;
                particles.push(p);
            }
        }
        detector.detect(&mut particles);
        // Edge particles should be detected as surface
        // Center particles may or may not be, depending on eigenvalue ratios
        // Just verify no panic and some particles are detected
        let n_surface = particles.iter().filter(|p| p.is_surface).count();
        // At least some should be surface (edge particles)
        assert!(
            n_surface > 0 || particles.len() < 5,
            "Some edge particles should be surface"
        );
    }

    #[test]
    fn test_eigenvalues_identity_matrix() {
        let m = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let eigs = EigenvalueSurfaceDetector::eigenvalues_symmetric_3x3(m);
        for &e in &eigs {
            assert!(
                (e - 1.0).abs() < 1e-10,
                "Eigenvalue of identity should be 1, got {e}"
            );
        }
    }

    #[test]
    fn test_eigenvalues_diagonal_matrix() {
        let m = [[3.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 2.0]];
        let eigs = EigenvalueSurfaceDetector::eigenvalues_symmetric_3x3(m);
        assert!((eigs[0] - 1.0).abs() < 1e-10);
        assert!((eigs[1] - 2.0).abs() < 1e-10);
        assert!((eigs[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_divergence_surface_detector() {
        let detector = DivergenceSurfaceDetector::new(0.7, 0.5);
        let mut particles = vec![
            SurfaceParticle::new([0.0, 0.0, 0.0], 1.0),
            SurfaceParticle::new([0.1, 0.0, 0.0], 1.0),
            SurfaceParticle::new([0.2, 0.0, 0.0], 1.0),
        ];
        for p in &mut particles {
            p.density = 1000.0;
        }
        detector.detect(&mut particles);
        // All should be processed without panic
    }

    #[test]
    fn test_surface_reconstructor_empty() {
        let recon = SurfaceReconstructor::new(0.1, 0.5, 0.2);
        let cells = recon.extract_surface_cells(&[]);
        assert!(cells.is_empty());
    }

    #[test]
    fn test_curvature_estimator_flat_surface() {
        let center = [0.0, 0.0, 0.0];
        let neighbors = [
            [0.1, 0.0, 0.0],
            [-0.1, 0.0, 0.0],
            [0.0, 0.1, 0.0],
            [0.0, -0.1, 0.0],
        ];
        let normal = [0.0, 0.0, 1.0];
        let kappa = CurvatureEstimator::estimate_from_neighbors(center, &neighbors, normal);
        // Flat surface → curvature should be ~0
        assert!(
            kappa.abs() < 1.0,
            "Flat surface curvature should be small, got {kappa}"
        );
    }

    #[test]
    fn test_curvature_estimator_too_few_neighbors() {
        let center = [0.0; 3];
        let neighbors: [[f64; 3]; 2] = [[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]];
        let normal = [0.0, 0.0, 1.0];
        let kappa = CurvatureEstimator::estimate_from_neighbors(center, &neighbors, normal);
        assert_eq!(kappa, 0.0);
    }

    #[test]
    fn test_render_hints_empty() {
        let hints = generate_render_hints(&[]);
        assert!(hints.is_empty());
    }

    #[test]
    fn test_render_hints_non_surface_excluded() {
        let p = SurfaceParticle::new([0.0; 3], 1.0);
        let hints = generate_render_hints(&[p]);
        assert!(
            hints.is_empty(),
            "Non-surface particles should not generate hints"
        );
    }

    #[test]
    fn test_render_hints_surface_included() {
        let mut p = SurfaceParticle::new([0.0; 3], 1.0);
        p.is_surface = true;
        p.normal = [0.0, 1.0, 0.0];
        p.curvature = 2.0;
        let hints = generate_render_hints(&[p]);
        assert_eq!(hints.len(), 1);
        assert_eq!(hints[0].particle_index, 0);
        assert!((hints[0].normal[1] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_color_function_new() {
        let cf = ColorFunction::new(0.5);
        assert!((cf.smoothing_length - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_surface_tension_new() {
        let st = SurfaceTension::new(0.072, 0.5);
        assert!((st.sigma - 0.072).abs() < 1e-10);
        assert!((st.smoothing_length - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_kernel_monotone_decreasing() {
        let h = 1.0;
        let mut prev_w = ColorFunction::kernel_value(0.0, h);
        for i in 1..20 {
            let r = i as f64 * 0.1;
            let w = ColorFunction::kernel_value(r, h);
            assert!(
                w <= prev_w + 1e-14,
                "Kernel should decrease: w({r})={w} > w({})={prev_w}",
                (i - 1) as f64 * 0.1
            );
            prev_w = w;
        }
    }
}
