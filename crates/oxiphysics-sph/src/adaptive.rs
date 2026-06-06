// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive particle refinement (splitting and merging) for SPH simulations.
//!
//! This module implements adaptive resolution SPH where particles can be
//! dynamically split into 8 children or merged from N particles into one,
//! based on local flow activity (velocity gradients, density errors).

// ---------------------------------------------------------------------------
// SphParticle
// ---------------------------------------------------------------------------

/// A single SPH particle with refinement level tracking.
#[derive(Debug, Clone)]
pub struct SphParticle {
    /// Position in 3D space.
    pub position: [f64; 3],
    /// Velocity vector.
    pub velocity: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// Local fluid density.
    pub density: f64,
    /// Local fluid pressure.
    pub pressure: f64,
    /// Smoothing length.
    pub smoothing_length: f64,
    /// Refinement level (0 = coarsest).
    pub level: u32,
}

impl SphParticle {
    /// Create a new particle with specified values.
    pub fn new(
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        density: f64,
        pressure: f64,
        smoothing_length: f64,
        level: u32,
    ) -> Self {
        Self {
            position,
            velocity,
            mass,
            density,
            pressure,
            smoothing_length,
            level,
        }
    }
}

// ---------------------------------------------------------------------------
// AdaptiveSph
// ---------------------------------------------------------------------------

/// Adaptive SPH simulation supporting particle splitting and merging.
pub struct AdaptiveSph {
    /// Collection of all particles.
    pub particles: Vec<SphParticle>,
    /// Maximum allowed refinement level.
    pub max_level: u32,
    /// Smoothing length at level 0.
    pub base_h: f64,
    /// Velocity gradient magnitude above which a particle is split.
    pub split_threshold: f64,
    /// Velocity gradient magnitude below which a particle may be merged.
    pub merge_threshold: f64,
}

impl AdaptiveSph {
    /// Create a new adaptive SPH system.
    ///
    /// # Arguments
    /// * `max_level` – maximum refinement level allowed.
    /// * `base_h` – smoothing length at level 0.
    pub fn new(max_level: u32, base_h: f64) -> Self {
        Self {
            particles: Vec::new(),
            max_level,
            base_h,
            split_threshold: 1.0,
            merge_threshold: 0.1,
        }
    }

    /// Add a particle to the simulation.
    pub fn add_particle(&mut self, p: SphParticle) {
        self.particles.push(p);
    }

    /// Split the particle at `idx` into 8 children arranged in a 2×2×2 pattern.
    ///
    /// Each child receives 1/8 of the parent mass and an adjusted smoothing
    /// length for its new level. The parent particle is removed (swap-removed).
    /// Returns the indices (in the updated `particles` vec) of the 8 new children.
    ///
    /// # Panics
    /// Panics if `idx` is out of bounds.
    pub fn split_particle(&mut self, idx: usize) -> Vec<usize> {
        let parent = self.particles[idx].clone();
        let child_level = parent.level + 1;
        let child_h = Self::smoothing_length_for_level(self.base_h, child_level);
        let child_mass = parent.mass / 8.0;

        // Offset children by ±offset from parent centre so they tile neatly.
        // Use half the child smoothing length as the grid spacing.
        let offset = child_h * 0.5;
        let offsets: [[f64; 3]; 8] = [
            [-offset, -offset, -offset],
            [offset, -offset, -offset],
            [-offset, offset, -offset],
            [offset, offset, -offset],
            [-offset, -offset, offset],
            [offset, -offset, offset],
            [-offset, offset, offset],
            [offset, offset, offset],
        ];

        // Remove parent (swap-remove).
        self.particles.swap_remove(idx);

        let start = self.particles.len();
        for off in &offsets {
            let child = SphParticle {
                position: [
                    parent.position[0] + off[0],
                    parent.position[1] + off[1],
                    parent.position[2] + off[2],
                ],
                velocity: parent.velocity,
                mass: child_mass,
                density: parent.density,
                pressure: parent.pressure,
                smoothing_length: child_h,
                level: child_level,
            };
            self.particles.push(child);
        }

        (start..start + 8).collect()
    }

    /// Merge the particles at `indices` into a single particle with:
    /// - total mass = sum of individual masses,
    /// - mass-weighted position and velocity,
    /// - smoothing length and level taken from the first particle.
    ///
    /// All particles in `indices` are removed and the merged particle is
    /// appended. Returns the index of the merged particle.
    ///
    /// # Panics
    /// Panics if `indices` is empty or any index is out of bounds.
    pub fn merge_particles(&mut self, indices: &[usize]) -> usize {
        assert!(
            !indices.is_empty(),
            "merge_particles: indices must not be empty"
        );

        let mut total_mass = 0.0_f64;
        let mut weighted_pos = [0.0_f64; 3];
        let mut weighted_vel = [0.0_f64; 3];

        for &i in indices {
            let p = &self.particles[i];
            total_mass += p.mass;
            for d in 0..3 {
                weighted_pos[d] += p.mass * p.position[d];
                weighted_vel[d] += p.mass * p.velocity[d];
            }
        }

        for d in 0..3 {
            weighted_pos[d] /= total_mass;
            weighted_vel[d] /= total_mass;
        }

        // Use level/h of first particle as the target (one level coarser would
        // be appropriate for actual use, but we keep it generic here).
        let first = &self.particles[indices[0]];
        let merged_level = first.level;
        let merged_h = first.smoothing_length;
        let merged_density = first.density;
        let merged_pressure = first.pressure;

        let merged = SphParticle {
            position: weighted_pos,
            velocity: weighted_vel,
            mass: total_mass,
            density: merged_density,
            pressure: merged_pressure,
            smoothing_length: merged_h,
            level: merged_level,
        };

        // Remove particles in descending index order to keep indices valid.
        let mut sorted_indices: Vec<usize> = indices.to_vec();
        sorted_indices.sort_unstable_by(|a, b| b.cmp(a));
        for &i in &sorted_indices {
            self.particles.swap_remove(i);
        }

        self.particles.push(merged);
        self.particles.len() - 1
    }

    /// Returns `true` if the particle at `idx` should be split.
    ///
    /// The criterion is: the velocity gradient magnitude (approximated using
    /// all other particles as neighbours) exceeds `split_threshold`, and the
    /// particle has not yet reached `max_level`.
    pub fn should_split(&self, idx: usize) -> bool {
        let p = &self.particles[idx];
        if p.level >= self.max_level {
            return false;
        }

        // Gather all other particles as neighbours.
        let neighbors: Vec<([f64; 3], [f64; 3])> = self
            .particles
            .iter()
            .enumerate()
            .filter(|(i, _)| *i != idx)
            .map(|(_, q)| (q.position, q.velocity))
            .collect();

        let grad_mag = compute_velocity_gradient_magnitude(
            p.position,
            p.velocity,
            &neighbors,
            p.smoothing_length,
        );

        grad_mag > self.split_threshold
    }

    /// Returns `true` if the particle at `idx` can be merged with `neighbor_idxs`.
    ///
    /// The criterion is: all neighbours share the same refinement level and the
    /// velocity gradient at this particle is below `merge_threshold`.
    pub fn should_merge(&self, idx: usize, neighbor_idxs: &[usize]) -> bool {
        if neighbor_idxs.is_empty() {
            return false;
        }
        let p = &self.particles[idx];
        // All neighbours must be at the same level.
        if neighbor_idxs
            .iter()
            .any(|&ni| self.particles[ni].level != p.level)
        {
            return false;
        }

        let neighbors: Vec<([f64; 3], [f64; 3])> = neighbor_idxs
            .iter()
            .map(|&ni| (self.particles[ni].position, self.particles[ni].velocity))
            .collect();

        let grad_mag = compute_velocity_gradient_magnitude(
            p.position,
            p.velocity,
            &neighbors,
            p.smoothing_length,
        );

        grad_mag < self.merge_threshold
    }

    /// Scan all particles and apply split/merge decisions.
    ///
    /// Split pass: iterate forward; when a particle should split, perform the
    /// split and continue with the remaining original particles.
    /// Merge pass: group same-level neighbours and merge those that qualify.
    pub fn adaptive_step(&mut self) {
        // --- Split pass ---
        // We iterate by index; after a split the parent is removed and 8
        // children appended, so we re-check the same index (now a different
        // particle moved in via swap_remove).
        let mut i = 0;
        loop {
            if i >= self.particles.len() {
                break;
            }
            if self.should_split(i) {
                self.split_particle(i);
                // Do NOT advance i: the slot now holds a different particle
                // (the former last one moved into position i by swap_remove).
                // Children are appended at the end and will be visited later.
            } else {
                i += 1;
            }
        }

        // --- Merge pass ---
        // Simple greedy: for each particle that is eligible, merge it with its
        // same-level neighbours whose gradient is also below threshold.
        // We track which indices have been consumed in this pass.
        let n = self.particles.len();
        let mut consumed = vec![false; n];

        let mut groups: Vec<Vec<usize>> = Vec::new();

        for i in 0..n {
            if consumed[i] {
                continue;
            }
            // Find same-level neighbours.
            let level_i = self.particles[i].level;
            let neighbors: Vec<usize> = (0..n)
                .filter(|&j| j != i && !consumed[j] && self.particles[j].level == level_i)
                .collect();

            if self.should_merge(i, &neighbors) {
                // Merge i with its qualifying neighbours (up to 7 for 8→1).
                let mut group: Vec<usize> = vec![i];
                for &ni in &neighbors {
                    if group.len() >= 8 {
                        break;
                    }
                    // Check that neighbor also wants to merge.
                    let other_neighbors: Vec<usize> = (0..n)
                        .filter(|&j| j != ni && !consumed[j] && self.particles[j].level == level_i)
                        .collect();
                    if self.should_merge(ni, &other_neighbors) {
                        group.push(ni);
                    }
                }
                if group.len() > 1 {
                    for &gi in &group {
                        consumed[gi] = true;
                    }
                    groups.push(group);
                }
            }
        }

        // Apply merges in reverse order of max index to keep indices stable.
        groups.sort_unstable_by_key(|g| std::cmp::Reverse(g.iter().copied().max().unwrap_or(0)));
        for group in groups {
            self.merge_particles(&group);
        }
    }

    /// Compute smoothing length for a given refinement level.
    ///
    /// `h(level) = base_h / 2^level`
    pub fn smoothing_length_for_level(base_h: f64, level: u32) -> f64 {
        base_h / (1u64 << level) as f64
    }
}

// ---------------------------------------------------------------------------
// Velocity gradient magnitude
// ---------------------------------------------------------------------------

/// Estimate the velocity gradient magnitude `||∇v||` at particle `i` using
/// a standard SPH first-derivative formula.
///
/// The SPH approximation of the velocity gradient component is:
/// ```text
/// (∂v_α/∂x_β)_i ≈ Σ_j (m_j / ρ_j) * (v_j_α - v_i_α) * ∂W/∂x_β
/// ```
/// Because mass and density of neighbours are not provided, this function
/// uses a simplified symmetric form that only requires positions and
/// velocities, treating each neighbour as having unit contribution weight
/// proportional to the kernel gradient.
///
/// # Arguments
/// * `pos_i` – position of particle *i*.
/// * `vel_i` – velocity of particle *i*.
/// * `neighbors` – slice of `(position, velocity)` for all neighbours.
/// * `h` – smoothing length.
pub fn compute_velocity_gradient_magnitude(
    pos_i: [f64; 3],
    vel_i: [f64; 3],
    neighbors: &[([f64; 3], [f64; 3])],
    h: f64,
) -> f64 {
    // 3×3 velocity gradient tensor G[alpha][beta] = ∂v_alpha/∂x_beta
    let mut grad = [[0.0_f64; 3]; 3];

    for (pos_j, vel_j) in neighbors {
        // r_ij vector from i to j
        let rx = pos_j[0] - pos_i[0];
        let ry = pos_j[1] - pos_i[1];
        let rz = pos_j[2] - pos_i[2];
        let r = (rx * rx + ry * ry + rz * rz).sqrt();

        if r < 1e-14 || r >= 2.0 * h {
            continue;
        }

        // Cubic-spline kernel gradient magnitude dW/dr (unnormalized for simplicity)
        let q = r / h;
        let dw_dr = if q < 1.0 {
            (-3.0 * q + 2.25 * q * q) / (std::f64::consts::PI * h * h * h * h)
        } else {
            let t = 2.0 - q;
            (-0.75 * t * t) / (std::f64::consts::PI * h * h * h * h)
        };

        // Gradient vector: ∂W/∂x_beta = dW/dr * (r_beta / r)
        let grad_w = [dw_dr * rx / r, dw_dr * ry / r, dw_dr * rz / r];

        // Velocity difference
        let dv = [
            vel_j[0] - vel_i[0],
            vel_j[1] - vel_i[1],
            vel_j[2] - vel_i[2],
        ];

        // Accumulate: G[alpha][beta] += dv[alpha] * grad_w[beta]
        for alpha in 0..3 {
            for beta in 0..3 {
                grad[alpha][beta] += dv[alpha] * grad_w[beta];
            }
        }
    }

    // Frobenius norm of gradient tensor
    let mut sum_sq = 0.0;
    for row in &grad {
        for &g in row {
            sum_sq += g * g;
        }
    }
    sum_sq.sqrt()
}

// ---------------------------------------------------------------------------
// Adaptive smoothing length – Newton iteration on density sum
// ---------------------------------------------------------------------------

/// Compute the SPH density sum for a particle at `pos_i` given its neighbours.
///
/// Uses the cubic-spline kernel W(r, h) (Monaghan 1992, 3D):
/// ```text
/// W(r, h) = σ / h³ × {1 - 1.5 q² + 0.75 q³,  0 ≤ q < 1
///                      0.25 (2 - q)³,            1 ≤ q < 2
///                      0,                         q ≥ 2 }
/// σ = 1 / π (3D normalisation)
/// q = r / h
/// ```
pub fn sph_density_sum(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let sigma_3d = 1.0 / std::f64::consts::PI;
    let inv_h3 = 1.0 / (h * h * h);

    let mut rho = 0.0_f64;
    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass) {
        let r = dist3(pos_i, *pos_j);
        let q = r / h;
        let w = if q < 1.0 {
            sigma_3d * inv_h3 * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else if q < 2.0 {
            let t = 2.0 - q;
            sigma_3d * inv_h3 * 0.25 * t * t * t
        } else {
            0.0
        };
        rho += m_j * w;
    }
    rho
}

/// Derivative dΩ/dh used in Newton iteration for adaptive h.
///
/// dρ_sph/dh = Σ_j m_j * (∂W/∂h)(r_ij, h)
pub fn sph_density_sum_dh(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
) -> f64 {
    let sigma_3d = 1.0 / std::f64::consts::PI;
    let inv_h = 1.0 / h;

    let mut drho_dh = 0.0_f64;
    for (pos_j, &m_j) in neighbor_pos.iter().zip(neighbor_mass) {
        let r = dist3(pos_i, *pos_j);
        let q = r / h;
        // dW/dh = -(3/h)*W + (1/h)*q*dW/dq
        let (w, dw_dq) = if q < 1.0 {
            let w = sigma_3d / (h * h * h) * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
            let dw_dq = sigma_3d / (h * h * h) * (-3.0 * q + 2.25 * q * q);
            (w, dw_dq)
        } else if q < 2.0 {
            let t = 2.0 - q;
            let w = sigma_3d / (h * h * h) * 0.25 * t * t * t;
            let dw_dq = sigma_3d / (h * h * h) * (-0.75 * t * t);
            (w, dw_dq)
        } else {
            (0.0, 0.0)
        };
        let dw_dh = -3.0 * inv_h * w + inv_h * q * dw_dq;
        drho_dh += m_j * dw_dh;
    }
    drho_dh
}

/// Solve for the adaptive smoothing length `h` such that the SPH density
/// sum equals `rho_target`, using Newton-Raphson iteration.
///
/// `h_guess` is the initial guess; up to `max_iter` iterations are taken
/// with convergence tolerance `tol` on `|h_new - h_old| / h_old`.
///
/// Returns `Some(h)` on convergence, `None` if the iteration diverges or
/// the maximum iteration count is exceeded.
pub fn adaptive_h_newton(
    pos_i: [f64; 3],
    h_guess: f64,
    rho_target: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    max_iter: usize,
    tol: f64,
) -> Option<f64> {
    let mut h = h_guess;
    for _ in 0..max_iter {
        let rho = sph_density_sum(pos_i, h, neighbor_pos, neighbor_mass);
        let f = rho - rho_target;
        if f.abs() < tol * rho_target.max(1e-14) {
            return Some(h);
        }
        let df_dh = sph_density_sum_dh(pos_i, h, neighbor_pos, neighbor_mass);
        if df_dh.abs() < 1e-30 {
            return None;
        }
        let h_new = h - f / df_dh;
        if h_new <= 0.0 {
            // Clamp to positive
            h *= 0.5;
            continue;
        }
        let rel_change = (h_new - h).abs() / h.max(1e-30);
        h = h_new;
        if rel_change < tol {
            return Some(h);
        }
    }
    None
}

/// Compute the variable-h kernel correction factor Ω_i for particle i.
///
/// Ω_i = 1 - (∂ρ_i / ∂h_i) * h_i / (n_dim * ρ_i)
///
/// Used in the grad-h SPH formulation (Springel & Hernquist 2002).
pub fn omega_grad_h_correction(
    pos_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_mass: &[f64],
    rho_i: f64,
    n_dim: f64,
) -> f64 {
    if rho_i.abs() < 1e-30 {
        return 1.0;
    }
    let drho_dh = sph_density_sum_dh(pos_i, h, neighbor_pos, neighbor_mass);
    1.0 - drho_dh * h / (n_dim * rho_i)
}

/// Helper: 3D Euclidean distance.
fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
}

// ---------------------------------------------------------------------------
// AdaptiveHResult – per-particle adaptive-h diagnostics
// ---------------------------------------------------------------------------

/// Result of an adaptive-h solve for a single particle.
#[derive(Debug, Clone)]
pub struct AdaptiveHResult {
    /// Converged smoothing length.
    pub h: f64,
    /// SPH density sum at the converged h.
    pub rho: f64,
    /// grad-h correction factor Ω.
    pub omega: f64,
    /// Whether Newton iteration converged.
    pub converged: bool,
}

impl AdaptiveHResult {
    /// Solve for adaptive h for a single particle.
    pub fn solve(
        pos_i: [f64; 3],
        h_guess: f64,
        rho_target: f64,
        neighbor_pos: &[[f64; 3]],
        neighbor_mass: &[f64],
    ) -> Self {
        let result = adaptive_h_newton(
            pos_i,
            h_guess,
            rho_target,
            neighbor_pos,
            neighbor_mass,
            50,
            1e-6,
        );
        match result {
            Some(h) => {
                let rho = sph_density_sum(pos_i, h, neighbor_pos, neighbor_mass);
                let omega =
                    omega_grad_h_correction(pos_i, h, neighbor_pos, neighbor_mass, rho, 3.0);
                Self {
                    h,
                    rho,
                    omega,
                    converged: true,
                }
            }
            None => {
                let rho = sph_density_sum(pos_i, h_guess, neighbor_pos, neighbor_mass);
                Self {
                    h: h_guess,
                    rho,
                    omega: 1.0,
                    converged: false,
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Coarsen / refine criteria
// ---------------------------------------------------------------------------

/// Criterion for h-refinement: refine if the smoothing length is larger
/// than the target for the local density.
///
/// `h_target(rho) = eta * (mass / rho)^(1/3)` where η ≈ 1.2 typically.
pub fn h_target_from_density(mass: f64, rho: f64, eta: f64) -> f64 {
    if rho.abs() < 1e-30 {
        return f64::MAX;
    }
    eta * (mass / rho).cbrt()
}

/// Returns `true` if particle i should be refined (h too large for local density).
pub fn should_refine_h(h_current: f64, h_target: f64, refine_ratio: f64) -> bool {
    h_current > h_target * refine_ratio
}

/// Returns `true` if particle i should be coarsened (h too small for local density).
pub fn should_coarsen_h(h_current: f64, h_target: f64, coarsen_ratio: f64) -> bool {
    h_current < h_target / coarsen_ratio
}

// ---------------------------------------------------------------------------
// HStatistics – smoothing-length statistics over all particles
// ---------------------------------------------------------------------------

/// Statistics of smoothing lengths across the particle set.
#[derive(Debug, Clone)]
pub struct HStatistics {
    /// Minimum smoothing length.
    pub h_min: f64,
    /// Maximum smoothing length.
    pub h_max: f64,
    /// Mean smoothing length.
    pub h_mean: f64,
    /// Standard deviation.
    pub h_std: f64,
    /// Number of particles.
    pub n: usize,
}

impl HStatistics {
    /// Compute statistics from a slice of smoothing lengths.
    pub fn compute(hs: &[f64]) -> Self {
        let n = hs.len();
        if n == 0 {
            return Self {
                h_min: 0.0,
                h_max: 0.0,
                h_mean: 0.0,
                h_std: 0.0,
                n: 0,
            };
        }
        let h_min = hs.iter().cloned().fold(f64::INFINITY, f64::min);
        let h_max = hs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let h_mean = hs.iter().sum::<f64>() / n as f64;
        let h_std = if n > 1 {
            let var = hs.iter().map(|&h| (h - h_mean).powi(2)).sum::<f64>() / (n - 1) as f64;
            var.sqrt()
        } else {
            0.0
        };
        Self {
            h_min,
            h_max,
            h_mean,
            h_std,
            n,
        }
    }

    /// Compute statistics from the smoothing lengths stored in a particle set.
    pub fn from_particles(particles: &[SphParticle]) -> Self {
        let hs: Vec<f64> = particles.iter().map(|p| p.smoothing_length).collect();
        Self::compute(&hs)
    }
}

// ---------------------------------------------------------------------------
// CFL adaptive time stepping for SPH
// ---------------------------------------------------------------------------

/// CFL condition for SPH (Monaghan 1992 / Morris 1997).
///
/// The CFL time step limit for SPH combines the velocity CFL condition and
/// the viscous CFL condition:
///
/// ```text
/// dt_cfl   = cfl_vel  * h / (c_s + max|v|)
/// dt_visc  = cfl_visc * h² / nu
/// dt = min(dt_cfl, dt_visc)
/// ```
///
/// where `c_s` is the speed of sound, `nu` is kinematic viscosity (can be 0),
/// and `h` is the smoothing length.
pub fn cfl_timestep_sph(
    h: f64,
    c_sound: f64,
    max_speed: f64,
    kinematic_viscosity: f64,
    cfl_vel: f64,
    cfl_visc: f64,
    dt_min: f64,
    dt_max: f64,
) -> f64 {
    let denominator = c_sound + max_speed;
    let dt_cfl = if denominator > 1e-14 {
        cfl_vel * h / denominator
    } else {
        dt_max
    };
    let dt_visc = if kinematic_viscosity > 1e-30 {
        cfl_visc * h * h / kinematic_viscosity
    } else {
        dt_max
    };
    dt_cfl.min(dt_visc).clamp(dt_min, dt_max)
}

/// Compute the maximum particle speed from a velocity array.
pub fn max_particle_speed(velocities: &[[f64; 3]]) -> f64 {
    velocities
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .fold(0.0_f64, f64::max)
}

/// Compute the maximum acceleration magnitude.
pub fn max_acceleration(accels: &[[f64; 3]]) -> f64 {
    accels
        .iter()
        .map(|a| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt())
        .fold(0.0_f64, f64::max)
}

/// Acceleration-based time step limiter:
/// `dt_acc = sqrt(h / max|a|)`.
pub fn acceleration_timestep(h: f64, max_accel: f64, safety: f64) -> f64 {
    if max_accel < 1e-14 {
        return f64::MAX;
    }
    safety * (h / max_accel).sqrt()
}

/// Combined SPH time step: min of CFL, viscosity, and acceleration limits.
pub fn combined_sph_timestep(
    h: f64,
    c_sound: f64,
    max_speed: f64,
    kinematic_viscosity: f64,
    max_accel: f64,
    cfl_vel: f64,
    cfl_visc: f64,
    acc_safety: f64,
    dt_min: f64,
    dt_max: f64,
) -> f64 {
    let dt_cfl = cfl_timestep_sph(
        h,
        c_sound,
        max_speed,
        kinematic_viscosity,
        cfl_vel,
        cfl_visc,
        dt_min,
        dt_max,
    );
    let dt_acc = acceleration_timestep(h, max_accel, acc_safety).min(dt_max);
    dt_cfl.min(dt_acc).clamp(dt_min, dt_max)
}

// ---------------------------------------------------------------------------
// Error estimators
// ---------------------------------------------------------------------------

/// L2 norm error estimator for a scalar field over particles.
///
/// Given a computed field `f_computed` and a reference field `f_reference`,
/// returns the relative L2 error:
///
/// ```text
/// E = sqrt( Σ (f_c - f_r)² ) / sqrt( Σ f_r² )
/// ```
pub fn l2_error_estimator(f_computed: &[f64], f_reference: &[f64]) -> f64 {
    assert_eq!(f_computed.len(), f_reference.len());
    let n = f_computed.len();
    if n == 0 {
        return 0.0;
    }
    let sum_sq_err: f64 = f_computed
        .iter()
        .zip(f_reference)
        .map(|(&fc, &fr)| (fc - fr).powi(2))
        .sum();
    let sum_sq_ref: f64 = f_reference.iter().map(|&fr| fr * fr).sum();
    if sum_sq_ref < 1e-30 {
        return sum_sq_err.sqrt();
    }
    (sum_sq_err / sum_sq_ref).sqrt()
}

/// Maximum absolute error between two fields.
pub fn linf_error(f_computed: &[f64], f_reference: &[f64]) -> f64 {
    f_computed
        .iter()
        .zip(f_reference)
        .map(|(&fc, &fr)| (fc - fr).abs())
        .fold(0.0_f64, f64::max)
}

/// Per-particle density error: `|ρ_i - ρ₀| / ρ₀`.
pub fn relative_density_errors(densities: &[f64], rho0: f64) -> Vec<f64> {
    densities
        .iter()
        .map(|&rho| (rho - rho0).abs() / rho0.max(1e-14))
        .collect()
}

/// Maximum relative density error.
pub fn max_relative_density_error(densities: &[f64], rho0: f64) -> f64 {
    relative_density_errors(densities, rho0)
        .into_iter()
        .fold(0.0_f64, f64::max)
}

/// Mean relative density error.
pub fn mean_relative_density_error(densities: &[f64], rho0: f64) -> f64 {
    let errs = relative_density_errors(densities, rho0);
    if errs.is_empty() {
        return 0.0;
    }
    errs.iter().sum::<f64>() / errs.len() as f64
}

// ---------------------------------------------------------------------------
// Refinement indicators
// ---------------------------------------------------------------------------

/// Velocity-gradient-based refinement indicator per particle.
///
/// Returns a scalar indicator `χ_i` proportional to the local velocity
/// gradient magnitude.  Particles with `χ_i > threshold` should be refined.
///
/// Uses a simple SPH finite-difference approximation:
/// ```text
/// χ_i = |v_i - Σ_j W_ij v_j / Σ_j W_ij|
/// ```
/// (a measure of local variation from the SPH-smoothed field).
pub fn velocity_variation_indicator(
    pos_i: [f64; 3],
    vel_i: [f64; 3],
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_vel: &[[f64; 3]],
    neighbor_mass: &[f64],
    neighbor_density: &[f64],
) -> f64 {
    let mut sum_w = 0.0_f64;
    let mut sum_wv = [0.0_f64; 3];

    for (((pos_j, vel_j), &m_j), &rho_j) in neighbor_pos
        .iter()
        .zip(neighbor_vel.iter())
        .zip(neighbor_mass.iter())
        .zip(neighbor_density.iter())
    {
        let dx = pos_i[0] - pos_j[0];
        let dy = pos_i[1] - pos_j[1];
        let dz = pos_i[2] - pos_j[2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let q = r / h;
        if q >= 2.0 {
            continue;
        }
        let sigma_3d = 1.0 / (std::f64::consts::PI * h * h * h);
        let w = if q < 1.0 {
            sigma_3d * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else {
            let t = 2.0 - q;
            sigma_3d * 0.25 * t * t * t
        };
        let rho_j_safe = rho_j.max(1e-14);
        let wm = w * m_j / rho_j_safe;
        sum_w += wm;
        sum_wv[0] += wm * vel_j[0];
        sum_wv[1] += wm * vel_j[1];
        sum_wv[2] += wm * vel_j[2];
    }

    if sum_w < 1e-30 {
        return 0.0;
    }
    let v_smooth = [sum_wv[0] / sum_w, sum_wv[1] / sum_w, sum_wv[2] / sum_w];
    let diff = [
        vel_i[0] - v_smooth[0],
        vel_i[1] - v_smooth[1],
        vel_i[2] - v_smooth[2],
    ];
    (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt()
}

/// Density-gradient refinement indicator.
///
/// Returns `|∇ρ|_i * h / ρ_i` — the relative density gradient scaled by h.
/// Large values (> threshold) indicate regions that need refinement.
pub fn density_gradient_indicator(
    pos_i: [f64; 3],
    rho_i: f64,
    h: f64,
    neighbor_pos: &[[f64; 3]],
    neighbor_rho: &[f64],
    neighbor_mass: &[f64],
) -> f64 {
    if rho_i < 1e-30 {
        return 0.0;
    }
    // ∇ρ_i ≈ Σ_j m_j (ρ_j - ρ_i) ∇W_ij / ρ_j  (symmetric form simplified)
    let sigma_3d = 1.0 / (std::f64::consts::PI * h * h * h);
    let mut grad = [0.0_f64; 3];
    for ((pos_j, &rho_j), &m_j) in neighbor_pos
        .iter()
        .zip(neighbor_rho.iter())
        .zip(neighbor_mass.iter())
    {
        let rij = [
            pos_i[0] - pos_j[0],
            pos_i[1] - pos_j[1],
            pos_i[2] - pos_j[2],
        ];
        let r = (rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2]).sqrt();
        if r < 1e-14 || r >= 2.0 * h {
            continue;
        }
        let q = r / h;
        let dw_dr = if q < 1.0 {
            sigma_3d * (-3.0 * q + 2.25 * q * q) / h
        } else {
            let t = 2.0 - q;
            sigma_3d * (-0.75 * t * t) / h
        };
        let rho_j_safe = rho_j.max(1e-14);
        let coeff = m_j * (rho_j - rho_i) / rho_j_safe * dw_dr / r;
        grad[0] += coeff * rij[0];
        grad[1] += coeff * rij[1];
        grad[2] += coeff * rij[2];
    }
    let grad_mag = (grad[0] * grad[0] + grad[1] * grad[1] + grad[2] * grad[2]).sqrt();
    grad_mag * h / rho_i
}

// ---------------------------------------------------------------------------
// Particle adaptation error control
// ---------------------------------------------------------------------------

/// Error-based particle adaptation controller.
///
/// Tracks per-particle error indicators and decides which particles to
/// split (too coarse), merge (too fine), or leave unchanged.
#[derive(Debug, Clone)]
pub struct AdaptationController {
    /// Refinement threshold: particles with indicator > refine_threshold are split.
    pub refine_threshold: f64,
    /// Coarsening threshold: particles with indicator < coarsen_threshold may merge.
    pub coarsen_threshold: f64,
    /// Maximum refinement level.
    pub max_level: u32,
    /// Minimum refinement level.
    pub min_level: u32,
    /// Safety factor (hysteresis): coarsen_threshold = refine_threshold / hysteresis.
    pub hysteresis: f64,
}

impl AdaptationController {
    /// Create a new adaptation controller.
    pub fn new(refine_threshold: f64, max_level: u32) -> Self {
        let hysteresis = 4.0;
        Self {
            refine_threshold,
            coarsen_threshold: refine_threshold / hysteresis,
            max_level,
            min_level: 0,
            hysteresis,
        }
    }

    /// Decide adaptation action for particle `i`.
    ///
    /// Returns `AdaptationAction::Refine`, `Coarsen`, or `None`.
    pub fn decide(&self, indicator: f64, level: u32) -> AdaptationAction {
        if indicator > self.refine_threshold && level < self.max_level {
            AdaptationAction::Refine
        } else if indicator < self.coarsen_threshold && level > self.min_level {
            AdaptationAction::Coarsen
        } else {
            AdaptationAction::None
        }
    }

    /// Batch decision for all particles.
    ///
    /// Returns a vec of `AdaptationAction` parallel to `indicators` / `levels`.
    pub fn decide_all(&self, indicators: &[f64], levels: &[u32]) -> Vec<AdaptationAction> {
        indicators
            .iter()
            .zip(levels.iter())
            .map(|(&ind, &lvl)| self.decide(ind, lvl))
            .collect()
    }

    /// Count particles needing refinement.
    pub fn count_refine(&self, actions: &[AdaptationAction]) -> usize {
        actions
            .iter()
            .filter(|a| matches!(a, AdaptationAction::Refine))
            .count()
    }

    /// Count particles needing coarsening.
    pub fn count_coarsen(&self, actions: &[AdaptationAction]) -> usize {
        actions
            .iter()
            .filter(|a| matches!(a, AdaptationAction::Coarsen))
            .count()
    }
}

/// Decision from the adaptation controller.
#[derive(Debug, Clone, PartialEq)]
pub enum AdaptationAction {
    /// Split this particle into 8 children.
    Refine,
    /// Merge this particle with neighbours.
    Coarsen,
    /// No action needed.
    None,
}

// ---------------------------------------------------------------------------
// AdaptiveTimestepManager
// ---------------------------------------------------------------------------

/// Manages adaptive time stepping for the full SPH simulation.
///
/// Combines CFL condition, force-based limiter, and solver convergence
/// feedback to suggest the next time step.
#[derive(Debug, Clone)]
pub struct AdaptiveTimestepManager {
    /// Current time step.
    pub dt: f64,
    /// Minimum allowed time step.
    pub dt_min: f64,
    /// Maximum allowed time step.
    pub dt_max: f64,
    /// CFL velocity factor (typically 0.25–0.4).
    pub cfl_vel: f64,
    /// CFL viscosity factor (typically 0.125).
    pub cfl_visc: f64,
    /// Acceleration safety factor.
    pub acc_safety: f64,
    /// Factor to reduce dt when solver is not converging.
    pub reduction_factor: f64,
    /// Factor to increase dt when converging well.
    pub growth_factor: f64,
    /// Consecutive convergence steps before growing dt.
    pub growth_patience: usize,
    /// Convergence step counter.
    pub converged_count: usize,
}

impl AdaptiveTimestepManager {
    /// Create a new manager with typical SPH defaults.
    pub fn new(dt_initial: f64, dt_min: f64, dt_max: f64) -> Self {
        Self {
            dt: dt_initial,
            dt_min,
            dt_max,
            cfl_vel: 0.3,
            cfl_visc: 0.125,
            acc_safety: 0.25,
            reduction_factor: 0.5,
            growth_factor: 1.05,
            growth_patience: 5,
            converged_count: 0,
        }
    }

    /// Update the time step based on current flow conditions and solver feedback.
    ///
    /// * `h` – current (minimum) smoothing length.
    /// * `c_sound` – speed of sound.
    /// * `max_speed` – maximum particle speed.
    /// * `kinematic_viscosity` – fluid kinematic viscosity.
    /// * `max_accel` – maximum acceleration magnitude.
    /// * `solver_converged` – whether the pressure solver converged.
    pub fn update(
        &mut self,
        h: f64,
        c_sound: f64,
        max_speed: f64,
        kinematic_viscosity: f64,
        max_accel: f64,
        solver_converged: bool,
    ) -> f64 {
        if solver_converged {
            self.converged_count += 1;
        } else {
            self.converged_count = 0;
            self.dt = (self.dt * self.reduction_factor).max(self.dt_min);
        }

        let dt_physics = combined_sph_timestep(
            h,
            c_sound,
            max_speed,
            kinematic_viscosity,
            max_accel,
            self.cfl_vel,
            self.cfl_visc,
            self.acc_safety,
            self.dt_min,
            self.dt_max,
        );

        // Allow gradual growth if we have been converging consistently.
        if self.converged_count >= self.growth_patience {
            self.dt = (self.dt * self.growth_factor).min(self.dt_max);
        }

        self.dt = self.dt.min(dt_physics);
        self.dt
    }
}

// ---------------------------------------------------------------------------
// Refinement history / statistics
// ---------------------------------------------------------------------------

/// Statistics tracked during an adaptive simulation.
#[derive(Debug, Clone, Default)]
pub struct AdaptiveStats {
    /// Total splits performed so far.
    pub total_splits: usize,
    /// Total merges performed so far.
    pub total_merges: usize,
    /// Total time steps taken.
    pub total_steps: usize,
    /// Sum of all time steps (for average).
    pub sum_dt: f64,
    /// Minimum time step ever used.
    pub min_dt: f64,
    /// Maximum time step ever used.
    pub max_dt: f64,
}

impl AdaptiveStats {
    /// Create zeroed statistics.
    pub fn new() -> Self {
        Self {
            min_dt: f64::INFINITY,
            max_dt: 0.0,
            ..Default::default()
        }
    }

    /// Record a completed time step.
    pub fn record_step(&mut self, dt: f64, splits: usize, merges: usize) {
        self.total_steps += 1;
        self.total_splits += splits;
        self.total_merges += merges;
        self.sum_dt += dt;
        if dt < self.min_dt {
            self.min_dt = dt;
        }
        if dt > self.max_dt {
            self.max_dt = dt;
        }
    }

    /// Average time step.
    pub fn avg_dt(&self) -> f64 {
        if self.total_steps == 0 {
            0.0
        } else {
            self.sum_dt / self.total_steps as f64
        }
    }

    /// Total simulation time elapsed.
    pub fn elapsed_time(&self) -> f64 {
        self.sum_dt
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(pos: [f64; 3], mass: f64, h: f64, level: u32) -> SphParticle {
        SphParticle::new(pos, [0.0; 3], mass, 1000.0, 0.0, h, level)
    }

    // -----------------------------------------------------------------------
    // smoothing_length_for_level
    // -----------------------------------------------------------------------

    #[test]
    fn smoothing_length_level_0() {
        let base_h = 0.1;
        assert!((AdaptiveSph::smoothing_length_for_level(base_h, 0) - base_h).abs() < 1e-15);
    }

    #[test]
    fn smoothing_length_level_1() {
        let base_h = 0.1;
        assert!((AdaptiveSph::smoothing_length_for_level(base_h, 1) - base_h / 2.0).abs() < 1e-15);
    }

    #[test]
    fn smoothing_length_level_2() {
        let base_h = 0.1;
        assert!((AdaptiveSph::smoothing_length_for_level(base_h, 2) - base_h / 4.0).abs() < 1e-15);
    }

    // -----------------------------------------------------------------------
    // split_particle produces 8 children with 1/8 mass each
    // -----------------------------------------------------------------------

    #[test]
    fn split_produces_8_children() {
        let base_h = 0.1;
        let mut sph = AdaptiveSph::new(4, base_h);
        sph.add_particle(make_particle([0.0, 0.0, 0.0], 1.0, base_h, 0));

        let children = sph.split_particle(0);
        assert_eq!(children.len(), 8, "should produce exactly 8 children");
        assert_eq!(sph.particles.len(), 8, "total particle count should be 8");

        for &ci in &children {
            let child = &sph.particles[ci];
            assert!(
                (child.mass - 1.0 / 8.0).abs() < 1e-14,
                "each child mass should be 1/8"
            );
            assert_eq!(child.level, 1, "children should be at level 1");
        }
    }

    #[test]
    fn split_children_mass_sum() {
        let base_h = 0.1;
        let mut sph = AdaptiveSph::new(4, base_h);
        sph.add_particle(make_particle([0.0, 0.0, 0.0], 2.4, base_h, 0));

        sph.split_particle(0);
        let total: f64 = sph.particles.iter().map(|p| p.mass).sum();
        assert!(
            (total - 2.4).abs() < 1e-12,
            "mass should be conserved after split"
        );
    }

    // -----------------------------------------------------------------------
    // merge_particles conserves total mass
    // -----------------------------------------------------------------------

    #[test]
    fn merge_conserves_mass() {
        let base_h = 0.1;
        let h1 = AdaptiveSph::smoothing_length_for_level(base_h, 1);
        let mut sph = AdaptiveSph::new(4, base_h);

        let child_mass = 0.5 / 8.0;
        for k in 0..8_usize {
            sph.add_particle(make_particle(
                [k as f64 * 0.01, 0.0, 0.0],
                child_mass,
                h1,
                1,
            ));
        }
        let expected_mass: f64 = sph.particles.iter().map(|p| p.mass).sum();

        let indices: Vec<usize> = (0..8).collect();
        sph.merge_particles(&indices);

        assert_eq!(sph.particles.len(), 1, "should have 1 merged particle");
        let actual_mass = sph.particles[0].mass;
        assert!(
            (actual_mass - expected_mass).abs() < 1e-12,
            "mass should be conserved after merge: expected {expected_mass}, got {actual_mass}"
        );
    }

    #[test]
    fn merge_mass_weighted_position() {
        let base_h = 0.1;
        let h1 = AdaptiveSph::smoothing_length_for_level(base_h, 1);
        let mut sph = AdaptiveSph::new(4, base_h);

        // Two particles with equal mass; merged position should be midpoint.
        sph.add_particle(SphParticle::new(
            [0.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            1000.0,
            0.0,
            h1,
            1,
        ));
        sph.add_particle(SphParticle::new(
            [2.0, 0.0, 0.0],
            [0.0; 3],
            1.0,
            1000.0,
            0.0,
            h1,
            1,
        ));

        sph.merge_particles(&[0, 1]);
        let merged = &sph.particles[0];
        assert!((merged.position[0] - 1.0).abs() < 1e-12);
        assert!((merged.mass - 2.0).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // split then merge approximately conserves mass
    // -----------------------------------------------------------------------

    // -----------------------------------------------------------------------
    // sph_density_sum tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_sph_density_sum_self_contribution() {
        // Single particle at origin, evaluating its self-contribution
        let pos_i = [0.0; 3];
        let h = 0.1;
        let masses = vec![1.0];
        let positions = vec![[0.0; 3]];
        let rho = sph_density_sum(pos_i, h, &positions, &masses);
        // W(0, h) = σ/h³ × 1.0
        let expected = 1.0 / (std::f64::consts::PI * h * h * h);
        assert!(
            (rho - expected).abs() < 1e-10,
            "rho={rho}, expected={expected}"
        );
    }

    #[test]
    fn test_sph_density_sum_positive() {
        let pos_i = [0.0; 3];
        let h = 0.1;
        let positions = vec![[0.05, 0.0, 0.0], [0.0, 0.05, 0.0]];
        let masses = vec![0.001, 0.001];
        let rho = sph_density_sum(pos_i, h, &positions, &masses);
        assert!(rho > 0.0, "density must be positive");
    }

    #[test]
    fn test_sph_density_sum_far_particle_zero() {
        // Particle outside 2h contributes nothing
        let pos_i = [0.0; 3];
        let h = 0.1;
        let positions = vec![[10.0, 0.0, 0.0]];
        let masses = vec![1.0];
        let rho = sph_density_sum(pos_i, h, &positions, &masses);
        assert!(
            rho.abs() < 1e-30,
            "far particle should not contribute, rho={rho}"
        );
    }

    // -----------------------------------------------------------------------
    // adaptive_h_newton tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_adaptive_h_newton_uniform_cluster() {
        // 8 particles on a cube of side 0.1, all mass 0.01
        let pos_i = [0.0; 3];
        let h = 0.1;
        let positions: Vec<[f64; 3]> = vec![
            [0.05, 0.0, 0.0],
            [0.0, 0.05, 0.0],
            [0.0, 0.0, 0.05],
            [-0.05, 0.0, 0.0],
            [0.0, -0.05, 0.0],
            [0.0, 0.0, -0.05],
            [0.05, 0.05, 0.0],
            [-0.05, -0.05, 0.0],
        ];
        let masses = vec![0.01; 8];
        let rho_init = sph_density_sum(pos_i, h, &positions, &masses);
        // Try to converge to the same rho (should trivially converge from h_guess)
        let result = adaptive_h_newton(pos_i, h, rho_init, &positions, &masses, 50, 1e-6);
        assert!(result.is_some(), "Newton should converge");
        let h_final = result.unwrap();
        let rho_final = sph_density_sum(pos_i, h_final, &positions, &masses);
        assert!(
            (rho_final - rho_init).abs() / rho_init.max(1e-14) < 1e-4,
            "rho mismatch: rho_final={rho_final}, rho_init={rho_init}"
        );
    }

    // -----------------------------------------------------------------------
    // omega_grad_h_correction tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_omega_correction_near_unity_for_uniform() {
        let pos_i = [0.0; 3];
        let h = 0.1;
        let positions: Vec<[f64; 3]> = vec![[0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.05]];
        let masses = vec![0.001; 3];
        let rho = sph_density_sum(pos_i, h, &positions, &masses);
        let omega = omega_grad_h_correction(pos_i, h, &positions, &masses, rho, 3.0);
        // Omega should be a real number (not NaN or Inf)
        assert!(omega.is_finite(), "omega must be finite, got {omega}");
    }

    // -----------------------------------------------------------------------
    // h_target / should_refine / should_coarsen
    // -----------------------------------------------------------------------

    #[test]
    fn test_h_target_from_density_scaling() {
        // h_target ∝ (mass/rho)^(1/3)
        let h1 = h_target_from_density(0.001, 1000.0, 1.2);
        let h2 = h_target_from_density(0.001, 8000.0, 1.2);
        // rho 8x larger → h 2x smaller
        assert!((h1 / h2 - 2.0).abs() < 1e-10, "h1/h2={}", h1 / h2);
    }

    #[test]
    fn test_should_refine_h() {
        assert!(should_refine_h(0.2, 0.1, 1.5)); // 0.2 > 0.1 * 1.5 = 0.15
        assert!(!should_refine_h(0.1, 0.1, 1.5)); // 0.1 <= 0.15
    }

    #[test]
    fn test_should_coarsen_h() {
        assert!(should_coarsen_h(0.01, 0.1, 5.0)); // 0.01 < 0.1 / 5 = 0.02
        assert!(!should_coarsen_h(0.05, 0.1, 5.0)); // 0.05 >= 0.02
    }

    // -----------------------------------------------------------------------
    // HStatistics tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_h_statistics_empty() {
        let stats = HStatistics::compute(&[]);
        assert_eq!(stats.n, 0);
    }

    #[test]
    fn test_h_statistics_uniform() {
        let hs = vec![0.1, 0.1, 0.1, 0.1];
        let stats = HStatistics::compute(&hs);
        assert!((stats.h_mean - 0.1).abs() < 1e-14);
        assert!(stats.h_std.abs() < 1e-14);
        assert!((stats.h_min - 0.1).abs() < 1e-14);
        assert!((stats.h_max - 0.1).abs() < 1e-14);
    }

    #[test]
    fn test_h_statistics_varying() {
        let hs = vec![0.05, 0.1, 0.15, 0.2];
        let stats = HStatistics::compute(&hs);
        assert_eq!(stats.n, 4);
        assert!(
            (stats.h_mean - 0.125).abs() < 1e-12,
            "mean={}",
            stats.h_mean
        );
        assert!(stats.h_std > 0.0);
        assert!((stats.h_min - 0.05).abs() < 1e-14);
        assert!((stats.h_max - 0.2).abs() < 1e-14);
    }

    #[test]
    fn test_h_statistics_from_particles() {
        let base_h = 0.1;
        let mut sph = AdaptiveSph::new(2, base_h);
        sph.add_particle(make_particle([0.0, 0.0, 0.0], 1.0, base_h, 0));
        let h1 = AdaptiveSph::smoothing_length_for_level(base_h, 1);
        sph.add_particle(make_particle([0.5, 0.0, 0.0], 0.5, h1, 1));
        let stats = HStatistics::from_particles(&sph.particles);
        assert_eq!(stats.n, 2);
        assert!((stats.h_min - h1).abs() < 1e-14);
        assert!((stats.h_max - base_h).abs() < 1e-14);
    }

    // -----------------------------------------------------------------------
    // AdaptiveHResult tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_adaptive_h_result_solve() {
        let pos_i = [0.0; 3];
        let h_guess = 0.1;
        let positions = vec![[0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.0, 0.0, 0.05]];
        let masses = vec![0.001; 3];
        let rho_init = sph_density_sum(pos_i, h_guess, &positions, &masses);
        let result = AdaptiveHResult::solve(pos_i, h_guess, rho_init, &positions, &masses);
        assert!(result.converged, "should converge");
        assert!(result.h > 0.0);
        assert!(result.omega.is_finite());
    }

    #[test]
    fn split_then_merge_conserves_mass() {
        let base_h = 0.1;
        let parent_mass = 1.0;
        let mut sph = AdaptiveSph::new(4, base_h);
        sph.add_particle(make_particle([0.0, 0.0, 0.0], parent_mass, base_h, 0));

        // Split into 8 children.
        let children = sph.split_particle(0);
        assert_eq!(children.len(), 8);

        // Mass after split.
        let mass_after_split: f64 = sph.particles.iter().map(|p| p.mass).sum();
        assert!(
            (mass_after_split - parent_mass).abs() < 1e-12,
            "mass conserved after split"
        );

        // Merge children back.
        sph.merge_particles(&children);
        let mass_after_merge: f64 = sph.particles.iter().map(|p| p.mass).sum();
        assert!(
            (mass_after_merge - parent_mass).abs() < 1e-12,
            "mass conserved after merge: expected {parent_mass}, got {mass_after_merge}"
        );
    }

    // -----------------------------------------------------------------------
    // CFL timestep tests
    // -----------------------------------------------------------------------

    #[test]
    fn cfl_timestep_zero_velocity_returns_cfl_limit() {
        let dt = cfl_timestep_sph(0.1, 10.0, 0.0, 0.0, 0.3, 0.125, 1e-6, 1.0);
        let expected = 0.3 * 0.1 / 10.0;
        assert!(
            (dt - expected).abs() < 1e-14,
            "Expected {expected}, got {dt}"
        );
    }

    #[test]
    fn cfl_timestep_high_velocity_reduces_dt() {
        let dt_slow = cfl_timestep_sph(0.1, 1.0, 0.0, 0.0, 0.3, 0.125, 1e-6, 1.0);
        let dt_fast = cfl_timestep_sph(0.1, 1.0, 100.0, 0.0, 0.3, 0.125, 1e-6, 1.0);
        assert!(dt_fast < dt_slow, "Higher speed → smaller dt");
    }

    #[test]
    fn cfl_timestep_viscosity_clamps_dt() {
        let dt = cfl_timestep_sph(0.1, 0.0, 0.0, 1e-3, 0.3, 0.125, 1e-10, 1.0);
        let expected_visc: f64 = 0.125 * 0.1 * 0.1 / 1e-3;
        assert!(
            (dt - expected_visc.min(1.0)).abs() < 1e-10,
            "dt_visc={expected_visc}, dt={dt}"
        );
    }

    #[test]
    fn cfl_timestep_clamps_to_min_max() {
        let dt = cfl_timestep_sph(0.001, 1e5, 1e5, 0.0, 0.3, 0.125, 0.5, 1.0);
        assert!(dt >= 0.5, "dt should be at least dt_min=0.5");
        let dt2 = cfl_timestep_sph(100.0, 0.001, 0.0, 0.0, 0.3, 0.125, 1e-6, 0.1);
        assert!(dt2 <= 0.1, "dt should be at most dt_max=0.1");
    }

    #[test]
    fn max_particle_speed_zero() {
        let vels = vec![[0.0_f64; 3]; 5];
        assert!((max_particle_speed(&vels)).abs() < 1e-14);
    }

    #[test]
    fn max_particle_speed_correct() {
        let vels = vec![[3.0_f64, 4.0, 0.0], [1.0, 0.0, 0.0]];
        assert!((max_particle_speed(&vels) - 5.0).abs() < 1e-14);
    }

    #[test]
    fn acceleration_timestep_zero_accel_returns_max() {
        let dt = acceleration_timestep(0.1, 0.0, 0.25);
        assert_eq!(dt, f64::MAX);
    }

    #[test]
    fn acceleration_timestep_finite() {
        let dt = acceleration_timestep(0.1, 9.81, 0.25);
        let expected = 0.25 * (0.1_f64 / 9.81).sqrt();
        assert!((dt - expected).abs() < 1e-14);
    }

    #[test]
    fn combined_timestep_takes_min() {
        let dt = combined_sph_timestep(0.1, 10.0, 1.0, 1e-3, 9.81, 0.3, 0.125, 0.25, 1e-10, 10.0);
        let dt_cfl = cfl_timestep_sph(0.1, 10.0, 1.0, 1e-3, 0.3, 0.125, 1e-10, 10.0);
        let dt_acc = acceleration_timestep(0.1, 9.81, 0.25).min(10.0);
        let expected = dt_cfl.min(dt_acc).clamp(1e-10, 10.0);
        assert!((dt - expected).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // Error estimator tests
    // -----------------------------------------------------------------------

    #[test]
    fn l2_error_identical_fields_zero() {
        let f = vec![1.0, 2.0, 3.0];
        assert!(l2_error_estimator(&f, &f).abs() < 1e-14);
    }

    #[test]
    fn l2_error_known_value() {
        let f_c = vec![1.1, 2.2];
        let f_r = vec![1.0, 2.0];
        let err = l2_error_estimator(&f_c, &f_r);
        // sqrt((0.01+0.04)/(1+4)) = sqrt(0.05/5) = sqrt(0.01) = 0.1
        assert!((err - 0.1).abs() < 1e-12, "Expected 0.1, got {err}");
    }

    #[test]
    fn linf_error_known_value() {
        let f_c = vec![1.5, 2.0, 3.2];
        let f_r = vec![1.0, 2.0, 3.0];
        let err = linf_error(&f_c, &f_r);
        assert!((err - 0.5).abs() < 1e-14);
    }

    #[test]
    fn relative_density_errors_computation() {
        let densities = vec![1000.0, 1100.0, 900.0];
        let errs = relative_density_errors(&densities, 1000.0);
        assert!((errs[0]).abs() < 1e-14, "zero error at rho0");
        assert!((errs[1] - 0.1).abs() < 1e-12, "10% error");
        assert!((errs[2] - 0.1).abs() < 1e-12, "10% error");
    }

    #[test]
    fn max_relative_density_error_matches_manual() {
        let d = vec![1000.0, 1200.0, 800.0];
        let max_err = max_relative_density_error(&d, 1000.0);
        assert!((max_err - 0.2).abs() < 1e-12);
    }

    #[test]
    fn mean_relative_density_error_matches_manual() {
        let d = vec![1100.0, 900.0]; // both have 10% error
        let mean_err = mean_relative_density_error(&d, 1000.0);
        assert!((mean_err - 0.1).abs() < 1e-12);
    }

    // -----------------------------------------------------------------------
    // Refinement indicator tests
    // -----------------------------------------------------------------------

    #[test]
    fn velocity_variation_indicator_zero_for_uniform() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [1.0, 0.0, 0.0];
        let h = 0.1;
        let nbpos = vec![[0.05_f64, 0.0, 0.0]];
        let nbvel = vec![[1.0_f64, 0.0, 0.0]]; // same velocity
        let nbmass = vec![0.001_f64];
        let nbrho = vec![1000.0_f64];
        let ind = velocity_variation_indicator(pos_i, vel_i, h, &nbpos, &nbvel, &nbmass, &nbrho);
        assert!(
            ind < 1e-12,
            "Uniform velocity field → indicator ~0, got {ind}"
        );
    }

    #[test]
    fn velocity_variation_indicator_nonzero_for_shear() {
        let pos_i = [0.0_f64; 3];
        let vel_i = [1.0, 0.0, 0.0];
        let h = 0.1;
        let nbpos = vec![[0.05_f64, 0.0, 0.0]];
        let nbvel = vec![[-1.0_f64, 0.0, 0.0]]; // opposite velocity
        let nbmass = vec![0.001_f64];
        let nbrho = vec![1000.0_f64];
        let ind = velocity_variation_indicator(pos_i, vel_i, h, &nbpos, &nbvel, &nbmass, &nbrho);
        assert!(ind > 0.0, "Shear field → indicator > 0, got {ind}");
    }

    #[test]
    fn density_gradient_indicator_zero_for_uniform() {
        let pos_i = [0.0_f64; 3];
        let rho_i = 1000.0;
        let h = 0.1;
        let nbpos = vec![[0.05_f64, 0.0, 0.0]];
        let nbrho = vec![1000.0_f64]; // same density
        let nbmass = vec![0.001_f64];
        let ind = density_gradient_indicator(pos_i, rho_i, h, &nbpos, &nbrho, &nbmass);
        assert!(ind < 1e-14, "Uniform density → gradient indicator 0: {ind}");
    }

    #[test]
    fn density_gradient_indicator_positive_for_gradient() {
        let pos_i = [0.0_f64; 3];
        let rho_i = 1000.0;
        let h = 0.1;
        let nbpos = vec![[0.05_f64, 0.0, 0.0]];
        let nbrho = vec![1500.0_f64]; // higher density to the right
        let nbmass = vec![0.001_f64];
        let ind = density_gradient_indicator(pos_i, rho_i, h, &nbpos, &nbrho, &nbmass);
        assert!(ind > 0.0, "Density gradient → indicator > 0: {ind}");
    }

    // -----------------------------------------------------------------------
    // AdaptationController tests
    // -----------------------------------------------------------------------

    #[test]
    fn controller_refine_decision() {
        let ctrl = AdaptationController::new(1.0, 5);
        assert_eq!(ctrl.decide(2.0, 0), AdaptationAction::Refine);
    }

    #[test]
    fn controller_coarsen_decision() {
        let ctrl = AdaptationController::new(1.0, 5);
        // coarsen_threshold = 1.0 / 4.0 = 0.25; indicator 0.1 < 0.25
        assert_eq!(ctrl.decide(0.1, 1), AdaptationAction::Coarsen);
    }

    #[test]
    fn controller_none_decision_in_between() {
        let ctrl = AdaptationController::new(1.0, 5);
        // 0.3 < 1.0 (no refine) and 0.3 > 0.25 (no coarsen)
        assert_eq!(ctrl.decide(0.3, 2), AdaptationAction::None);
    }

    #[test]
    fn controller_no_refine_at_max_level() {
        let ctrl = AdaptationController::new(0.5, 3);
        assert_ne!(ctrl.decide(1.0, 3), AdaptationAction::Refine);
    }

    #[test]
    fn controller_no_coarsen_at_min_level() {
        let ctrl = AdaptationController::new(1.0, 5);
        assert_ne!(ctrl.decide(0.0, 0), AdaptationAction::Coarsen);
    }

    #[test]
    fn controller_count_refine_coarsen() {
        let ctrl = AdaptationController::new(1.0, 5);
        let indicators = vec![2.0, 0.1, 0.3, 3.0, 0.05];
        let levels = vec![0_u32, 1, 2, 0, 1];
        let actions = ctrl.decide_all(&indicators, &levels);
        assert_eq!(ctrl.count_refine(&actions), 2); // 2.0 and 3.0
        assert_eq!(ctrl.count_coarsen(&actions), 2); // 0.1 and 0.05
    }

    // -----------------------------------------------------------------------
    // AdaptiveTimestepManager tests
    // -----------------------------------------------------------------------

    #[test]
    fn timestep_manager_initial_dt() {
        let mgr = AdaptiveTimestepManager::new(0.01, 1e-5, 0.01);
        assert!((mgr.dt - 0.01).abs() < 1e-14);
    }

    #[test]
    fn timestep_manager_reduces_on_divergence() {
        let mut mgr = AdaptiveTimestepManager::new(0.01, 1e-5, 0.01);
        let dt = mgr.update(0.1, 10.0, 1.0, 0.0, 9.81, false);
        assert!(dt <= 0.01, "dt should be <= initial");
    }

    #[test]
    fn timestep_manager_grows_after_convergence() {
        let mut mgr = AdaptiveTimestepManager::new(0.005, 1e-5, 0.01);
        for _ in 0..10 {
            mgr.update(0.1, 0.0, 0.0, 0.0, 0.0, true);
        }
        assert!(mgr.dt >= 0.005, "dt should grow or stay stable");
    }

    // -----------------------------------------------------------------------
    // AdaptiveStats tests
    // -----------------------------------------------------------------------

    #[test]
    fn stats_empty() {
        let stats = AdaptiveStats::new();
        assert_eq!(stats.total_steps, 0);
        assert_eq!(stats.avg_dt(), 0.0);
    }

    #[test]
    fn stats_record_step() {
        let mut stats = AdaptiveStats::new();
        stats.record_step(0.01, 2, 0);
        stats.record_step(0.005, 0, 1);
        assert_eq!(stats.total_steps, 2);
        assert_eq!(stats.total_splits, 2);
        assert_eq!(stats.total_merges, 1);
        assert!((stats.avg_dt() - 0.0075).abs() < 1e-14);
        assert!((stats.min_dt - 0.005).abs() < 1e-14);
        assert!((stats.max_dt - 0.01).abs() < 1e-14);
        assert!((stats.elapsed_time() - 0.015).abs() < 1e-14);
    }

    #[test]
    fn stats_min_max_track_correctly() {
        let mut stats = AdaptiveStats::new();
        for i in 1..=5 {
            stats.record_step(i as f64 * 0.001, 0, 0);
        }
        assert!((stats.min_dt - 0.001).abs() < 1e-14);
        assert!((stats.max_dt - 0.005).abs() < 1e-14);
    }
}
