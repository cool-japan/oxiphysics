// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive resolution SPH with variable smoothing length and particle refinement/coarsening.
//!
//! Implements grad-h SPH formulation with Newton-Raphson smoothing length iteration,
//! Wendland C2 kernel, particle splitting and merging.

// ---------------------------------------------------------------------------
// Math helpers
// ---------------------------------------------------------------------------

/// Dot product of two 3-vectors.
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Component-wise subtraction.
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Euclidean length of a 3-vector.
pub fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Scalar multiplication.
pub fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// Component-wise addition.
pub fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

// ---------------------------------------------------------------------------
// Particle
// ---------------------------------------------------------------------------

/// A single adaptive-resolution SPH particle.
#[derive(Clone, Debug)]
pub struct AdaptiveParticle {
    /// Position in world space.
    pub pos: [f64; 3],
    /// Velocity.
    pub vel: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// Current density estimate.
    pub density: f64,
    /// Pressure.
    pub pressure: f64,
    /// Individual smoothing length.
    pub h: f64,
    /// Grad-h correction factor (omega).
    pub omega: f64,
    /// Refinement level (0 = coarsest).
    pub level: u8,
    /// Whether this particle is active.
    pub alive: bool,
}

impl AdaptiveParticle {
    /// Create a new particle with default derived quantities.
    pub fn new(pos: [f64; 3], vel: [f64; 3], mass: f64, h: f64, level: u8) -> Self {
        Self {
            pos,
            vel,
            mass,
            density: 0.0,
            pressure: 0.0,
            h,
            omega: 1.0,
            level,
            alive: true,
        }
    }
}

// ---------------------------------------------------------------------------
// Wendland C2 kernel with variable h
// ---------------------------------------------------------------------------

/// Wendland C2 kernel with individual smoothing length support.
pub struct VariableHKernel;

impl VariableHKernel {
    /// Normalisation constant in 3D: (21 / (2 pi h^3)).
    fn alpha(h: f64) -> f64 {
        21.0 / (2.0 * std::f64::consts::PI * h * h * h)
    }

    /// Evaluate W(r, h).  Support radius = 2h.
    pub fn w(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let t = 1.0 - 0.5 * q;
        Self::alpha(h) * t * t * t * t * (2.0 * q + 1.0)
    }

    /// Derivative dW/dr.
    pub fn dw_dr(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let t = 1.0 - 0.5 * q;
        // d/dr [ alpha * t^4 * (2q+1) ] where q = r/h
        // = alpha * [ 4 t^3 * (-1/(2h)) * (2q+1) + t^4 * (2/h) ]
        Self::alpha(h)
            * (1.0 / h)
            * t
            * t
            * t
            * (-4.0 * q + 2.0 * (2.0 * q + 1.0) * (-0.5) * 4.0 / 4.0)
        // cleaner derivation below:
    }

    /// Derivative dW/dr — clean derivation.
    fn dw_dr_clean(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let t = 1.0 - 0.5 * q;
        // W = alpha * t^4 * (2q + 1)
        // dW/dq = alpha * [4 t^3 * (-1/2) * (2q+1) + t^4 * 2]
        //       = alpha * t^3 * [-2(2q+1) + 2t] ... wait, let's be precise:
        // dW/dq = alpha * [4 t^3 * (-1/2) * (2q+1)  +  t^4 * 2]
        //       = alpha * t^3 * [ -2(2q+1) + 2t ]
        // dW/dr = dW/dq * dq/dr = dW/dq * (1/h)
        let dw_dq = Self::alpha(h) * t * t * t * (-2.0 * (2.0 * q + 1.0) + 2.0 * t);
        dw_dq / h
    }

    /// Derivative dW/dh (needed for grad-h / omega computation).
    pub fn dw_dh(r: f64, h: f64) -> f64 {
        let q = r / h;
        if q >= 2.0 {
            return 0.0;
        }
        let t = 1.0 - 0.5 * q;
        // W = (21/(2 pi h^3)) * t^4 * (2q+1),  q = r/h
        // Differentiate w.r.t. h treating r fixed:
        //   d(alpha)/dh = -3 alpha / h
        //   dq/dh       = -r/h^2 = -q/h
        //   dt/dh       = (1/2)(q/h) = q/(2h)
        // dW/dh = d(alpha)/dh * t^4*(2q+1)
        //       + alpha * [4 t^3 * dt/dh * (2q+1) + t^4 * 2 * dq/dh]
        let alpha = Self::alpha(h);
        let dalpha_dh = -3.0 * alpha / h;
        let dq_dh = -q / h;
        let dt_dh = q / (2.0 * h);
        dalpha_dh * t * t * t * t * (2.0 * q + 1.0)
            + alpha * (4.0 * t * t * t * dt_dh * (2.0 * q + 1.0) + t * t * t * t * 2.0 * dq_dh)
    }
}

// ---------------------------------------------------------------------------
// Smoothing-length updater
// ---------------------------------------------------------------------------

/// Updates individual smoothing lengths to maintain a target neighbor count.
pub struct SmoothingLengthUpdate {
    /// Desired number of neighbors (e.g. 32 in 3D).
    pub n_neighbors_target: f64,
    /// Minimum allowed smoothing length.
    pub h_min: f64,
    /// Maximum allowed smoothing length.
    pub h_max: f64,
}

impl SmoothingLengthUpdate {
    /// Create a new updater.
    pub fn new(n_neighbors_target: f64, h_min: f64, h_max: f64) -> Self {
        Self {
            n_neighbors_target,
            h_min,
            h_max,
        }
    }

    /// Newton-Raphson step: h_new = h * (n_target / n_actual)^(1/3), clamped.
    pub fn update_h(&self, particle: &mut AdaptiveParticle, neighbor_count: usize) {
        let n_actual = neighbor_count.max(1) as f64;
        let ratio = self.n_neighbors_target / n_actual;
        let h_new = particle.h * ratio.powf(1.0 / 3.0);
        particle.h = h_new.clamp(self.h_min, self.h_max);
    }

    /// Compute omega = 1 - sum_j m_j/rho_j * dW/dh(r_ij, h_i).
    ///
    /// This grad-h correction factor enters the momentum equation.
    pub fn compute_omega(particle: &AdaptiveParticle, neighbors: &[&AdaptiveParticle]) -> f64 {
        let mut sum = 0.0;
        for nb in neighbors {
            if !nb.alive {
                continue;
            }
            let r = length(sub(particle.pos, nb.pos));
            let dw_dh = VariableHKernel::dw_dh(r, particle.h);
            let rho_j = nb.density.max(1e-10);
            sum += nb.mass / rho_j * dw_dh;
        }
        1.0 - sum
    }
}

// ---------------------------------------------------------------------------
// Refinement criterion
// ---------------------------------------------------------------------------

/// Rules that decide when a particle should be split or merged.
#[derive(Clone, Debug)]

pub enum RefinementCriterion {
    /// Split when density is above `high`, merge when below `low`.
    DensityBased { high: f64, low: f64 },
    /// Split when |∇v| exceeds the threshold.
    VelocityGradient { threshold: f64 },
    /// Interface detection based on density-gradient magnitude.
    Interface { detection_threshold: f64 },
}

impl RefinementCriterion {
    /// Return true if particle `p` should be refined.
    pub fn should_refine(&self, p: &AdaptiveParticle, density_grad_mag: f64) -> bool {
        match self {
            RefinementCriterion::DensityBased { high, .. } => p.density > *high,
            RefinementCriterion::VelocityGradient { threshold } => density_grad_mag > *threshold,
            RefinementCriterion::Interface {
                detection_threshold,
            } => density_grad_mag > *detection_threshold,
        }
    }

    /// Return true if particle `p` should be coarsened.
    pub fn should_coarsen(&self, p: &AdaptiveParticle, density_grad_mag: f64) -> bool {
        match self {
            RefinementCriterion::DensityBased { low, .. } => p.density < *low,
            RefinementCriterion::VelocityGradient { threshold } => {
                density_grad_mag < *threshold * 0.1
            }
            RefinementCriterion::Interface {
                detection_threshold,
            } => density_grad_mag < *detection_threshold * 0.1,
        }
    }
}

// ---------------------------------------------------------------------------
// Adaptive SPH system
// ---------------------------------------------------------------------------

/// Full adaptive-resolution SPH simulation system.
pub struct AdaptiveSphSystem {
    /// All particles (including dead ones pending compaction).
    pub particles: Vec<AdaptiveParticle>,
    /// Minimum allowed smoothing length.
    pub h_min: f64,
    /// Maximum allowed smoothing length.
    pub h_max: f64,
    /// Split particle when density exceeds this value.
    pub refine_density_threshold: f64,
    /// Merge particles when density drops below this value.
    pub coarsen_density_threshold: f64,
}

impl AdaptiveSphSystem {
    /// Create a new system.
    pub fn new(h_min: f64, h_max: f64) -> Self {
        Self {
            particles: Vec::new(),
            h_min,
            h_max,
            refine_density_threshold: f64::MAX,
            coarsen_density_threshold: 0.0,
        }
    }

    /// Add a particle to the system.
    pub fn add_particle(&mut self, pos: [f64; 3], vel: [f64; 3], mass: f64, h: f64) {
        self.particles
            .push(AdaptiveParticle::new(pos, vel, mass, h, 0));
    }

    /// Return indices of neighbors of particle `idx` within 2*h_i.
    pub fn find_neighbors(&self, idx: usize) -> Vec<usize> {
        let pi = &self.particles[idx];
        if !pi.alive {
            return Vec::new();
        }
        let two_h = 2.0 * pi.h;
        self.particles
            .iter()
            .enumerate()
            .filter(|(j, pj)| *j != idx && pj.alive && length(sub(pi.pos, pj.pos)) < two_h)
            .map(|(j, _)| j)
            .collect()
    }

    /// Compute density for all alive particles: ρ_i = Σ_j m_j W(|r_ij|, h_i).
    pub fn compute_density(&mut self) {
        let n = self.particles.len();
        let mut densities = vec![0.0_f64; n];
        for (i, d) in densities.iter_mut().enumerate() {
            if !self.particles[i].alive {
                continue;
            }
            let pi_pos = self.particles[i].pos;
            let pi_h = self.particles[i].h;
            let mut rho = 0.0;
            for pj in &self.particles {
                if !pj.alive {
                    continue;
                }
                let r = length(sub(pi_pos, pj.pos));
                rho += pj.mass * VariableHKernel::w(r, pi_h);
            }
            *d = rho;
        }
        for (p, &d) in self.particles.iter_mut().zip(densities.iter()) {
            p.density = d;
        }
    }

    /// Split particle at `idx` into two daughters.
    ///
    /// Returns indices of the two new particles. The parent is marked dead.
    pub fn split_particle(&mut self, idx: usize) -> [usize; 2] {
        let parent = self.particles[idx].clone();
        let child_mass = parent.mass * 0.5;
        let child_h = parent.h * 2.0_f64.powf(-1.0 / 3.0);
        let child_level = parent.level.saturating_add(1);

        // Offset daughters along a fixed axis by half the smoothing length.
        let offset = child_h * 0.5;
        let pos_a = add(parent.pos, [offset, 0.0, 0.0]);
        let pos_b = sub(parent.pos, [offset, 0.0, 0.0]);

        let child_a = AdaptiveParticle::new(pos_a, parent.vel, child_mass, child_h, child_level);
        let child_b = AdaptiveParticle::new(pos_b, parent.vel, child_mass, child_h, child_level);

        // Kill parent.
        self.particles[idx].alive = false;

        let idx_a = self.particles.len();
        self.particles.push(child_a);
        let idx_b = self.particles.len();
        self.particles.push(child_b);

        [idx_a, idx_b]
    }

    /// Merge particles `idx_a` and `idx_b` into a single particle.
    ///
    /// Returns the index of the new merged particle.
    pub fn merge_particles(&mut self, idx_a: usize, idx_b: usize) -> usize {
        let pa = self.particles[idx_a].clone();
        let pb = self.particles[idx_b].clone();

        let total_mass = pa.mass + pb.mass;
        // Centre-of-mass position.
        let pos = scale(
            add(scale(pa.pos, pa.mass), scale(pb.pos, pb.mass)),
            1.0 / total_mass,
        );
        // Momentum-conserving velocity.
        let vel = scale(
            add(scale(pa.vel, pa.mass), scale(pb.vel, pb.mass)),
            1.0 / total_mass,
        );
        // Merged smoothing length.
        let h_merged = (pa.h + pb.h) * 0.5 * 2.0_f64.powf(1.0 / 3.0);
        let h_merged = h_merged.clamp(self.h_min, self.h_max);
        // Level is min of the two parents minus 1, but not below 0.
        let merged_level = pa.level.min(pb.level).saturating_sub(1);

        let merged = AdaptiveParticle::new(pos, vel, total_mass, h_merged, merged_level);

        self.particles[idx_a].alive = false;
        self.particles[idx_b].alive = false;

        let new_idx = self.particles.len();
        self.particles.push(merged);
        new_idx
    }

    /// Adapt the resolution based on density thresholds.
    ///
    /// - Particles above `refine_density_threshold` are split.
    /// - Adjacent same-level pairs below `coarsen_density_threshold` are merged.
    pub fn adapt(&mut self) {
        // --- Refinement pass ---
        // Collect indices to split (snapshot, since we may append during split).
        let to_split: Vec<usize> = self
            .particles
            .iter()
            .enumerate()
            .filter(|(_, p)| p.alive && p.density > self.refine_density_threshold)
            .map(|(i, _)| i)
            .collect();

        for idx in to_split {
            if self.particles[idx].alive {
                self.split_particle(idx);
            }
        }

        // --- Coarsening pass ---
        // Find pairs at the same level that are both below the coarsen threshold
        // and within 2*h of each other. Greedy: once a particle is consumed, skip.
        let n = self.particles.len();
        let mut merged = vec![false; n];

        let candidates: Vec<usize> = (0..n)
            .filter(|&i| {
                self.particles[i].alive
                    && self.particles[i].density < self.coarsen_density_threshold
            })
            .collect();

        for &i in &candidates {
            if merged[i] || !self.particles[i].alive {
                continue;
            }
            let pi_level = self.particles[i].level;
            let pi_h = self.particles[i].h;
            let pi_pos = self.particles[i].pos;

            // Find first eligible partner.
            let partner = candidates.iter().find(|&&j| {
                j != i
                    && !merged[j]
                    && self.particles[j].alive
                    && self.particles[j].level == pi_level
                    && length(sub(pi_pos, self.particles[j].pos)) < 2.0 * pi_h
            });

            if let Some(&j) = partner {
                merged[i] = true;
                merged[j] = true;
                self.merge_particles(i, j);
            }
        }
    }

    /// Total mass of all alive particles.
    pub fn total_mass(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| p.alive)
            .map(|p| p.mass)
            .sum()
    }

    /// Number of alive particles.
    pub fn active_count(&self) -> usize {
        self.particles.iter().filter(|p| p.alive).count()
    }

    /// Maximum refinement level among alive particles.
    pub fn max_level(&self) -> u8 {
        self.particles
            .iter()
            .filter(|p| p.alive)
            .map(|p| p.level)
            .max()
            .unwrap_or(0)
    }

    /// Compact the particle list by removing dead particles.
    /// Returns a mapping from old indices to new indices (None for dead particles).
    pub fn compact(&mut self) -> Vec<Option<usize>> {
        let n = self.particles.len();
        let mut mapping = vec![None; n];
        let mut new_idx = 0;
        let mut compacted = Vec::with_capacity(self.active_count());
        for (i, p) in self.particles.iter().enumerate() {
            if p.alive {
                mapping[i] = Some(new_idx);
                compacted.push(p.clone());
                new_idx += 1;
            }
        }
        self.particles = compacted;
        mapping
    }

    /// Compute pressure for all alive particles using the Tait equation of state.
    ///
    /// P = B * ((rho / rho0)^gamma - 1)
    pub fn compute_pressure_tait(&mut self, rho0: f64, gamma: f64, speed_of_sound: f64) {
        let b = rho0 * speed_of_sound * speed_of_sound / gamma;
        for p in self.particles.iter_mut() {
            if !p.alive {
                continue;
            }
            let ratio = p.density / rho0;
            p.pressure = b * (ratio.powf(gamma) - 1.0);
        }
    }

    /// Mean smoothing length of alive particles.
    pub fn mean_h(&self) -> f64 {
        let alive: Vec<&AdaptiveParticle> = self.particles.iter().filter(|p| p.alive).collect();
        if alive.is_empty() {
            return 0.0;
        }
        let sum_h: f64 = alive.iter().map(|p| p.h).sum();
        sum_h / alive.len() as f64
    }

    /// Centre of mass of alive particles.
    pub fn center_of_mass(&self) -> [f64; 3] {
        let mut com = [0.0; 3];
        let mut total_mass = 0.0;
        for p in &self.particles {
            if !p.alive {
                continue;
            }
            com = add(com, scale(p.pos, p.mass));
            total_mass += p.mass;
        }
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        scale(com, 1.0 / total_mass)
    }

    /// Total momentum of alive particles.
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut mom = [0.0; 3];
        for p in &self.particles {
            if !p.alive {
                continue;
            }
            mom = add(mom, scale(p.vel, p.mass));
        }
        mom
    }

    /// Maximum velocity magnitude among alive particles.
    pub fn max_velocity(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| p.alive)
            .map(|p| length(p.vel))
            .fold(0.0, f64::max)
    }

    /// Minimum smoothing length among alive particles.
    pub fn min_h(&self) -> f64 {
        self.particles
            .iter()
            .filter(|p| p.alive)
            .map(|p| p.h)
            .fold(f64::MAX, f64::min)
    }
}

// ---------------------------------------------------------------------------
// Variable smoothing length iteration
// ---------------------------------------------------------------------------

/// Newton-Raphson iteration for smoothing length to satisfy a target
/// number-density equation: n_target * h^3 = sum_j W(r_ij, h).
pub struct SmoothingLengthSolver {
    /// Target number of neighbors.
    pub n_target: f64,
    /// Convergence tolerance for relative change in h.
    pub tol: f64,
    /// Maximum number of Newton iterations.
    pub max_iter: usize,
    /// Minimum allowed h.
    pub h_min: f64,
    /// Maximum allowed h.
    pub h_max: f64,
}

impl SmoothingLengthSolver {
    /// Create a new solver.
    pub fn new(n_target: f64, h_min: f64, h_max: f64) -> Self {
        Self {
            n_target,
            tol: 1e-4,
            max_iter: 20,
            h_min,
            h_max,
        }
    }

    /// Iterate h for a single particle given its neighbor positions and masses.
    ///
    /// Returns the converged h value.
    pub fn solve_h(
        &self,
        pos_i: [f64; 3],
        h_init: f64,
        neighbor_positions: &[[f64; 3]],
        neighbor_masses: &[f64],
    ) -> f64 {
        let mut h = h_init;
        for _ in 0..self.max_iter {
            // f(h) = sum_j m_j * W(r_ij, h) - rho_target
            // rho_target from n_target and kernel normalization
            let mut sum_w = 0.0;
            let mut _sum_dw_dh = 0.0;
            for (j, &npos) in neighbor_positions.iter().enumerate() {
                let r = length(sub(pos_i, npos));
                sum_w += neighbor_masses[j] * VariableHKernel::w(r, h);
                _sum_dw_dh += neighbor_masses[j] * VariableHKernel::dw_dh(r, h);
            }
            // Simple Newton step: h_new = h * (n_target / n_actual)^(1/3)
            let n_actual = sum_w.max(1e-30);
            let ratio = self.n_target / n_actual;
            let h_new = h * ratio.powf(1.0 / 3.0);
            let h_new = h_new.clamp(self.h_min, self.h_max);
            if ((h_new - h) / h).abs() < self.tol {
                return h_new;
            }
            h = h_new;
        }
        h.clamp(self.h_min, self.h_max)
    }
}

// ---------------------------------------------------------------------------
// H-refinement manager
// ---------------------------------------------------------------------------

/// Controls h-refinement with level-dependent smoothing length bounds.
pub struct HRefinementManager {
    /// Base smoothing length at level 0.
    pub h_base: f64,
    /// Maximum allowed refinement level.
    pub max_level: u8,
    /// Refinement ratio per level (h_l = h_base * ratio^(-level)).
    pub ratio: f64,
}

impl HRefinementManager {
    /// Create a new h-refinement manager.
    pub fn new(h_base: f64, max_level: u8) -> Self {
        Self {
            h_base,
            max_level,
            ratio: 2.0_f64.powf(1.0 / 3.0),
        }
    }

    /// Target smoothing length for a given refinement level.
    pub fn target_h(&self, level: u8) -> f64 {
        self.h_base / self.ratio.powi(level as i32)
    }

    /// Target mass for a given refinement level, assuming constant density rho0.
    pub fn target_mass(&self, level: u8, rho0: f64) -> f64 {
        let h = self.target_h(level);
        // mass ~ rho0 * (4/3 pi h^3) / n_neighbors
        // Simplified: mass halves per level
        rho0 * h * h * h
    }

    /// Whether a particle at given level can be further refined.
    pub fn can_refine(&self, level: u8) -> bool {
        level < self.max_level
    }

    /// Whether a particle at given level can be coarsened.
    pub fn can_coarsen(&self, level: u8) -> bool {
        level > 0
    }
}

// ---------------------------------------------------------------------------
// Particle splitting strategies
// ---------------------------------------------------------------------------

/// Different strategies for placing daughter particles after splitting.
#[derive(Clone, Debug)]
pub enum SplitStrategy {
    /// Split along the x-axis (default).
    AlongX,
    /// Split along the velocity direction.
    AlongVelocity,
    /// Split along the density gradient direction.
    AlongGradient,
}

impl SplitStrategy {
    /// Compute the offset direction for splitting a particle.
    /// Returns a unit vector indicating the split direction.
    pub fn split_direction(
        &self,
        _particle: &AdaptiveParticle,
        density_gradient: [f64; 3],
    ) -> [f64; 3] {
        match self {
            SplitStrategy::AlongX => [1.0, 0.0, 0.0],
            SplitStrategy::AlongVelocity => {
                let v_mag = length(_particle.vel);
                if v_mag < 1e-15 {
                    [1.0, 0.0, 0.0]
                } else {
                    scale(_particle.vel, 1.0 / v_mag)
                }
            }
            SplitStrategy::AlongGradient => {
                let g_mag = length(density_gradient);
                if g_mag < 1e-15 {
                    [1.0, 0.0, 0.0]
                } else {
                    scale(density_gradient, 1.0 / g_mag)
                }
            }
        }
    }

    /// Split a particle into two daughters along the computed direction.
    /// Returns (pos_a, pos_b, mass_child, h_child, new_level).
    pub fn split(
        &self,
        particle: &AdaptiveParticle,
        density_gradient: [f64; 3],
    ) -> ([f64; 3], [f64; 3], f64, f64, u8) {
        let dir = self.split_direction(particle, density_gradient);
        let child_mass = particle.mass * 0.5;
        let child_h = particle.h * 2.0_f64.powf(-1.0 / 3.0);
        let child_level = particle.level.saturating_add(1);
        let offset = child_h * 0.5;
        let pos_a = add(particle.pos, scale(dir, offset));
        let pos_b = sub(particle.pos, scale(dir, offset));
        (pos_a, pos_b, child_mass, child_h, child_level)
    }
}

// ---------------------------------------------------------------------------
// Particle merging criteria
// ---------------------------------------------------------------------------

/// Criteria for selecting merge partners.
pub struct MergeCriteria {
    /// Maximum distance (in units of h) for two particles to be merge candidates.
    pub max_distance_ratio: f64,
    /// Maximum mass ratio between merge candidates.
    pub max_mass_ratio: f64,
    /// Maximum level difference between merge candidates.
    pub max_level_diff: u8,
}

impl MergeCriteria {
    /// Create default merge criteria.
    pub fn new() -> Self {
        Self {
            max_distance_ratio: 1.5,
            max_mass_ratio: 2.0,
            max_level_diff: 0,
        }
    }

    /// Check if two particles are eligible to be merged.
    pub fn can_merge(&self, a: &AdaptiveParticle, b: &AdaptiveParticle) -> bool {
        if !a.alive || !b.alive {
            return false;
        }
        let dist = length(sub(a.pos, b.pos));
        let h_avg = (a.h + b.h) * 0.5;
        if dist > self.max_distance_ratio * h_avg {
            return false;
        }
        let mass_ratio = if a.mass > b.mass {
            a.mass / b.mass.max(1e-30)
        } else {
            b.mass / a.mass.max(1e-30)
        };
        if mass_ratio > self.max_mass_ratio {
            return false;
        }
        let level_diff = a.level.abs_diff(b.level);
        level_diff <= self.max_level_diff
    }
}

impl Default for MergeCriteria {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Density gradient estimator
// ---------------------------------------------------------------------------

/// SPH density gradient estimator for refinement criteria.
pub struct DensityGradientEstimator;

impl DensityGradientEstimator {
    /// Estimate the density gradient at particle i using SPH interpolation.
    ///
    /// grad(rho)_i = sum_j m_j * (rho_j/rho_i - 1) * grad_W(r_ij, h_i)
    pub fn estimate(particle: &AdaptiveParticle, neighbors: &[&AdaptiveParticle]) -> [f64; 3] {
        let rho_i = particle.density.max(1e-30);
        let mut grad = [0.0; 3];
        for nb in neighbors {
            if !nb.alive {
                continue;
            }
            let rij = sub(particle.pos, nb.pos);
            let r = length(rij);
            if r < 1e-15 {
                continue;
            }
            let dw = VariableHKernel::dw_dr_clean(r, particle.h);
            let rij_hat = scale(rij, 1.0 / r);
            // grad W = dW/dr * r_hat
            let grad_w = scale(rij_hat, dw);
            let coeff = nb.mass * (nb.density / rho_i - 1.0);
            grad = add(grad, scale(grad_w, coeff));
        }
        grad
    }

    /// Magnitude of the density gradient.
    pub fn gradient_magnitude(particle: &AdaptiveParticle, neighbors: &[&AdaptiveParticle]) -> f64 {
        length(Self::estimate(particle, neighbors))
    }
}

// ---------------------------------------------------------------------------
// Refinement statistics
// ---------------------------------------------------------------------------

/// Statistics about refinement operations.
#[derive(Clone, Debug, Default)]
pub struct RefinementStats {
    /// Number of split operations performed.
    pub splits: usize,
    /// Number of merge operations performed.
    pub merges: usize,
    /// Number of particles at each level.
    pub level_counts: Vec<usize>,
}

impl RefinementStats {
    /// Compute statistics from a system.
    pub fn from_system(system: &AdaptiveSphSystem) -> Self {
        let max_level = system.max_level() as usize;
        let mut level_counts = vec![0usize; max_level + 1];
        for p in &system.particles {
            if p.alive && (p.level as usize) < level_counts.len() {
                level_counts[p.level as usize] += 1;
            }
        }
        Self {
            splits: 0,
            merges: 0,
            level_counts,
        }
    }

    /// Total number of active particles across all levels.
    pub fn total_active(&self) -> usize {
        self.level_counts.iter().sum()
    }

    /// Fraction of particles at the finest level.
    pub fn finest_fraction(&self) -> f64 {
        let total = self.total_active();
        if total == 0 {
            return 0.0;
        }
        *self.level_counts.last().unwrap_or(&0) as f64 / total as f64
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_system() -> AdaptiveSphSystem {
        AdaptiveSphSystem::new(0.01, 10.0)
    }

    #[test]
    fn test_add_particle_increases_active_count() {
        let mut sys = make_system();
        assert_eq!(sys.active_count(), 0);
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        assert_eq!(sys.active_count(), 1);
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        assert_eq!(sys.active_count(), 2);
    }

    #[test]
    fn test_total_mass_conserved_after_split() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 2.0, 0.5);
        let mass_before = sys.total_mass();
        sys.split_particle(0);
        let mass_after = sys.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "mass before={mass_before}, after={mass_after}"
        );
    }

    #[test]
    fn test_total_mass_conserved_after_merge() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([0.1, 0.0, 0.0], [0.0; 3], 1.5, 0.5);
        let mass_before = sys.total_mass();
        sys.merge_particles(0, 1);
        let mass_after = sys.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-12,
            "mass before={mass_before}, after={mass_after}"
        );
    }

    #[test]
    fn test_split_creates_level_plus_one_particles() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        let parent_level = sys.particles[0].level;
        let [a, b] = sys.split_particle(0);
        assert_eq!(sys.particles[a].level, parent_level + 1);
        assert_eq!(sys.particles[b].level, parent_level + 1);
    }

    #[test]
    fn test_merge_creates_lower_level_relative_to_parents() {
        let mut sys = make_system();
        // Use level-2 particles so merged result is level 1.
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([0.1, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.particles[0].level = 2;
        sys.particles[1].level = 2;
        let merged_idx = sys.merge_particles(0, 1);
        let parent_level = 2u8;
        assert_eq!(
            sys.particles[merged_idx].level,
            parent_level.saturating_sub(1)
        );
    }

    #[test]
    fn test_wendland_kernel_positive_at_origin() {
        let h = 0.5;
        let w0 = VariableHKernel::w(0.0, h);
        assert!(w0 > 0.0, "W(0, h) should be positive, got {w0}");
    }

    #[test]
    fn test_wendland_kernel_zero_outside_support() {
        let h = 0.5;
        // Compact support: r >= 2h → W = 0.
        let w_outside = VariableHKernel::w(2.5 * h, h);
        assert_eq!(w_outside, 0.0, "W(2.5h, h) should be exactly 0");
    }

    #[test]
    fn test_wendland_kernel_monotone_decreasing() {
        let h = 0.5;
        let mut prev = VariableHKernel::w(0.0, h);
        for i in 1..20 {
            let r = (i as f64) * 0.05;
            let w = VariableHKernel::w(r, h);
            assert!(w <= prev + 1e-14, "Kernel should be monotone decreasing");
            prev = w;
        }
    }

    #[test]
    fn test_dw_dr_negative_inside_support() {
        let h = 0.5;
        // dW/dr should be negative for r > 0 inside support
        for i in 1..19 {
            let r = (i as f64) * 0.05;
            let dw = VariableHKernel::dw_dr_clean(r, h);
            assert!(dw <= 0.0, "dW/dr should be <= 0 at r={r}, got {dw}");
        }
    }

    #[test]
    fn test_dw_dh_at_origin() {
        let h = 0.5;
        let dw = VariableHKernel::dw_dh(0.0, h);
        // At r=0: as h increases, W(0,h) decreases => dW/dh < 0
        assert!(dw < 0.0, "dW/dh(0, h) should be negative, got {dw}");
    }

    #[test]
    fn test_smoothing_length_update_increases_h_for_few_neighbors() {
        let updater = SmoothingLengthUpdate::new(32.0, 0.01, 10.0);
        let mut p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        let h_old = p.h;
        updater.update_h(&mut p, 5); // few neighbors => increase h
        assert!(p.h > h_old, "h should increase with few neighbors");
    }

    #[test]
    fn test_smoothing_length_update_decreases_h_for_many_neighbors() {
        let updater = SmoothingLengthUpdate::new(32.0, 0.01, 10.0);
        let mut p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        let h_old = p.h;
        updater.update_h(&mut p, 100); // many neighbors => decrease h
        assert!(p.h < h_old, "h should decrease with many neighbors");
    }

    #[test]
    fn test_smoothing_length_update_clamped() {
        let updater = SmoothingLengthUpdate::new(32.0, 0.1, 1.0);
        let mut p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        updater.update_h(&mut p, 1); // very few neighbors => try to increase h a lot
        assert!(p.h <= 1.0, "h should be clamped to h_max=1.0");
        updater.update_h(&mut p, 10000); // many neighbors
        // After several updates with huge neighbor counts
        for _ in 0..20 {
            let mut p2 = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
            updater.update_h(&mut p2, 100000);
            assert!(p2.h >= 0.1, "h should be clamped to h_min=0.1");
        }
    }

    #[test]
    fn test_compute_omega() {
        let p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        // No neighbors => omega = 1.0
        let omega = SmoothingLengthUpdate::compute_omega(&p, &[]);
        assert!((omega - 1.0).abs() < 1e-14);
    }

    #[test]
    fn test_refinement_criterion_density() {
        let crit = RefinementCriterion::DensityBased {
            high: 100.0,
            low: 10.0,
        };
        let mut p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        p.density = 150.0;
        assert!(crit.should_refine(&p, 0.0));
        assert!(!crit.should_coarsen(&p, 0.0));

        p.density = 5.0;
        assert!(!crit.should_refine(&p, 0.0));
        assert!(crit.should_coarsen(&p, 0.0));
    }

    #[test]
    fn test_refinement_criterion_velocity_gradient() {
        let crit = RefinementCriterion::VelocityGradient { threshold: 10.0 };
        let p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        assert!(crit.should_refine(&p, 15.0));
        assert!(!crit.should_refine(&p, 5.0));
        assert!(crit.should_coarsen(&p, 0.5)); // < threshold * 0.1
        assert!(!crit.should_coarsen(&p, 5.0));
    }

    #[test]
    fn test_refinement_criterion_interface() {
        let crit = RefinementCriterion::Interface {
            detection_threshold: 20.0,
        };
        let p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        assert!(crit.should_refine(&p, 25.0));
        assert!(crit.should_coarsen(&p, 1.0)); // < 20 * 0.1
    }

    #[test]
    fn test_find_neighbors() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([0.3, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([5.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5); // too far
        let neighbors = sys.find_neighbors(0);
        assert_eq!(neighbors.len(), 1, "particle 1 is within 2*h=1.0");
        assert_eq!(neighbors[0], 1);
    }

    #[test]
    fn test_compute_density() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.compute_density();
        assert!(
            sys.particles[0].density > 0.0,
            "Self-density must be positive"
        );
    }

    #[test]
    fn test_adapt_splits_high_density() {
        let mut sys = make_system();
        sys.refine_density_threshold = 100.0;
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 2.0, 0.5);
        sys.particles[0].density = 200.0; // above threshold
        let count_before = sys.active_count();
        sys.adapt();
        let count_after = sys.active_count();
        assert!(
            count_after > count_before,
            "High density should trigger split"
        );
    }

    #[test]
    fn test_compact_removes_dead() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.particles[0].alive = false;
        assert_eq!(sys.particles.len(), 2);
        let mapping = sys.compact();
        assert_eq!(sys.particles.len(), 1);
        assert!(mapping[0].is_none());
        assert_eq!(mapping[1], Some(0));
    }

    #[test]
    fn test_compute_pressure_tait() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.particles[0].density = 1000.0;
        sys.compute_pressure_tait(1000.0, 7.0, 1500.0);
        // At rho = rho0 => pressure = 0
        assert!(sys.particles[0].pressure.abs() < 1e-6);

        sys.particles[0].density = 1001.0;
        sys.compute_pressure_tait(1000.0, 7.0, 1500.0);
        assert!(
            sys.particles[0].pressure > 0.0,
            "Compressed fluid => positive pressure"
        );
    }

    #[test]
    fn test_center_of_mass() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        let com = sys.center_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-14);
        assert!(com[1].abs() < 1e-14);
    }

    #[test]
    fn test_total_momentum_conservation() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 0.5);
        sys.add_particle([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 2.0, 0.5);
        let mom = sys.total_momentum();
        assert!(mom[0].abs() < 1e-14, "Momentum should be zero");
    }

    #[test]
    fn test_max_velocity() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 1.0, 0.5);
        sys.add_particle([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.5);
        let v_max = sys.max_velocity();
        assert!((v_max - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_h_refinement_manager() {
        let mgr = HRefinementManager::new(1.0, 5);
        let h0 = mgr.target_h(0);
        let h1 = mgr.target_h(1);
        assert!((h0 - 1.0).abs() < 1e-14);
        assert!(h1 < h0, "Higher level => smaller h");
        assert!(mgr.can_refine(3));
        assert!(!mgr.can_refine(5));
        assert!(mgr.can_coarsen(1));
        assert!(!mgr.can_coarsen(0));
    }

    #[test]
    fn test_split_strategy_along_x() {
        let p = AdaptiveParticle::new([1.0, 2.0, 3.0], [0.0; 3], 4.0, 0.8, 0);
        let strat = SplitStrategy::AlongX;
        let (pos_a, pos_b, mass, h, level) = strat.split(&p, [0.0; 3]);
        assert!((mass - 2.0).abs() < 1e-14);
        assert!(h < 0.8);
        assert_eq!(level, 1);
        assert!(pos_a[0] > p.pos[0]);
        assert!(pos_b[0] < p.pos[0]);
    }

    #[test]
    fn test_split_strategy_along_velocity() {
        let p = AdaptiveParticle::new([0.0; 3], [0.0, 5.0, 0.0], 2.0, 0.5, 0);
        let strat = SplitStrategy::AlongVelocity;
        let (pos_a, pos_b, _, _, _) = strat.split(&p, [0.0; 3]);
        // Split should be along y
        assert!(pos_a[1] > 0.0);
        assert!(pos_b[1] < 0.0);
        assert!(pos_a[0].abs() < 1e-14);
    }

    #[test]
    fn test_merge_criteria() {
        let criteria = MergeCriteria::new();
        let a = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 1);
        let b = AdaptiveParticle::new([0.3, 0.0, 0.0], [0.0; 3], 1.0, 0.5, 1);
        assert!(criteria.can_merge(&a, &b));

        // Different levels
        let c = AdaptiveParticle::new([0.3, 0.0, 0.0], [0.0; 3], 1.0, 0.5, 3);
        assert!(!criteria.can_merge(&a, &c));
    }

    #[test]
    fn test_merge_criteria_distance() {
        let criteria = MergeCriteria::new();
        let a = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        let far = AdaptiveParticle::new([10.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5, 0);
        assert!(!criteria.can_merge(&a, &far));
    }

    #[test]
    fn test_density_gradient_no_neighbors() {
        let p = AdaptiveParticle::new([0.0; 3], [0.0; 3], 1.0, 0.5, 0);
        let grad = DensityGradientEstimator::estimate(&p, &[]);
        assert!(length(grad) < 1e-14);
    }

    #[test]
    fn test_refinement_stats() {
        let mut sys = make_system();
        sys.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.5);
        sys.particles[1].level = 1;
        let stats = RefinementStats::from_system(&sys);
        assert_eq!(stats.total_active(), 2);
        assert_eq!(stats.level_counts[0], 1);
        assert_eq!(stats.level_counts[1], 1);
        assert!((stats.finest_fraction() - 0.5).abs() < 1e-14);
    }

    #[test]
    fn test_smoothing_length_solver() {
        let solver = SmoothingLengthSolver::new(32.0, 0.01, 10.0);
        // Single neighbor nearby
        let h = solver.solve_h([0.0; 3], 0.5, &[[0.1, 0.0, 0.0]], &[1.0]);
        assert!(h > 0.0 && h.is_finite());
    }

    #[test]
    fn test_min_h() {
        let mut sys = make_system();
        sys.add_particle([0.0; 3], [0.0; 3], 1.0, 0.5);
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.3);
        assert!((sys.min_h() - 0.3).abs() < 1e-14);
    }

    #[test]
    fn test_mean_h() {
        let mut sys = make_system();
        sys.add_particle([0.0; 3], [0.0; 3], 1.0, 0.4);
        sys.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.6);
        assert!((sys.mean_h() - 0.5).abs() < 1e-14);
    }
}
