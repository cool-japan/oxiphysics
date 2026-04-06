// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Adaptive SPH refinement: particle splitting/merging, variable smoothing
//! length, CFL time stepping, and resolution management.

#![allow(dead_code)]

use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// 1. RefinementCriteria — error indicators
// ---------------------------------------------------------------------------

/// Error indicators used to decide whether a particle needs refinement.
///
/// Three physics-based indicators are provided:
/// - density gradient magnitude
/// - velocity curl magnitude (vorticity)
/// - pressure Laplacian magnitude
#[derive(Debug, Clone)]
pub struct RefinementCriteria {
    /// Threshold for density gradient |∇ρ| above which splitting is triggered.
    pub density_gradient_threshold: f64,
    /// Threshold for velocity curl |∇×v| (vorticity) for refinement.
    pub velocity_curl_threshold: f64,
    /// Threshold for |∇²p| (pressure Laplacian) for refinement.
    pub pressure_laplacian_threshold: f64,
}

impl RefinementCriteria {
    /// Create refinement criteria with the given thresholds.
    pub fn new(
        density_gradient_threshold: f64,
        velocity_curl_threshold: f64,
        pressure_laplacian_threshold: f64,
    ) -> Self {
        Self {
            density_gradient_threshold,
            velocity_curl_threshold,
            pressure_laplacian_threshold,
        }
    }

    /// Returns `true` if any indicator exceeds its threshold.
    pub fn needs_refinement(
        &self,
        density_gradient: f64,
        velocity_curl: f64,
        pressure_laplacian: f64,
    ) -> bool {
        density_gradient.abs() > self.density_gradient_threshold
            || velocity_curl.abs() > self.velocity_curl_threshold
            || pressure_laplacian.abs() > self.pressure_laplacian_threshold
    }

    /// Returns `true` if all indicators are below half the thresholds (safe to merge).
    pub fn can_merge(
        &self,
        density_gradient: f64,
        velocity_curl: f64,
        pressure_laplacian: f64,
    ) -> bool {
        density_gradient.abs() < 0.5 * self.density_gradient_threshold
            && velocity_curl.abs() < 0.5 * self.velocity_curl_threshold
            && pressure_laplacian.abs() < 0.5 * self.pressure_laplacian_threshold
    }

    /// Weighted error indicator (scalar score for ranking particles).
    pub fn error_score(
        &self,
        density_gradient: f64,
        velocity_curl: f64,
        pressure_laplacian: f64,
    ) -> f64 {
        (density_gradient / self.density_gradient_threshold)
            .abs()
            .max(
                (velocity_curl / self.velocity_curl_threshold)
                    .abs()
                    .max((pressure_laplacian / self.pressure_laplacian_threshold).abs()),
            )
    }
}

// ---------------------------------------------------------------------------
// 2. ParticleSplit — mass/volume-conserving particle splitting
// ---------------------------------------------------------------------------

/// Splits one parent particle into N children conserving mass and momentum.
///
/// Children are placed symmetrically around the parent using a regular
/// polygon layout (2-D) or tetrahedral/octahedral patterns (3-D).
#[derive(Debug, Clone)]
pub struct ParticleSplit {
    /// Number of children to produce from a single split.
    pub n_children: usize,
    /// Offset fraction of smoothing length for child placement.
    pub placement_fraction: f64,
}

impl ParticleSplit {
    /// Create a splitter producing `n_children` children.
    ///
    /// `placement_fraction` ∈ (0, 1] controls how far from the parent centre
    /// the children are placed (fraction of smoothing length h).
    pub fn new(n_children: usize, placement_fraction: f64) -> Self {
        assert!(n_children >= 2, "must split into at least 2 children");
        Self {
            n_children,
            placement_fraction: placement_fraction.clamp(1e-3, 1.0),
        }
    }

    /// Child mass (each child carries equal fraction of parent mass).
    pub fn child_mass(&self, parent_mass: f64) -> f64 {
        parent_mass / self.n_children as f64
    }

    /// Child smoothing length h_c = h_p * n^(-1/d) (3-D, d=3).
    pub fn child_smoothing_length(&self, parent_h: f64) -> f64 {
        parent_h * (self.n_children as f64).powf(-1.0 / 3.0)
    }

    /// Generate child positions in 2-D around the parent at `pos`.
    ///
    /// Children are placed on a regular polygon of radius
    /// `placement_fraction * h`.
    pub fn split_positions_2d(&self, pos: [f64; 2], h: f64) -> Vec<[f64; 2]> {
        let r = self.placement_fraction * h;
        let n = self.n_children;
        (0..n)
            .map(|k| {
                let theta = 2.0 * PI * k as f64 / n as f64;
                [pos[0] + r * theta.cos(), pos[1] + r * theta.sin()]
            })
            .collect()
    }

    /// Generate child positions in 3-D using a latitude/longitude subdivision.
    pub fn split_positions_3d(&self, pos: [f64; 3], h: f64) -> Vec<[f64; 3]> {
        let r = self.placement_fraction * h;
        let n = self.n_children;
        let mut positions = Vec::with_capacity(n);
        for k in 0..n {
            let phi = PI * (k as f64 + 0.5) / n as f64;
            let theta = 2.0 * PI * k as f64 * 1.618_033_988; // golden-angle
            positions.push([
                pos[0] + r * phi.sin() * theta.cos(),
                pos[1] + r * phi.sin() * theta.sin(),
                pos[2] + r * phi.cos(),
            ]);
        }
        positions
    }

    /// Split velocity: all children inherit parent velocity (momentum conserving).
    pub fn child_velocity(&self, parent_velocity: [f64; 3]) -> [f64; 3] {
        parent_velocity
    }

    /// Verify mass conservation: total child mass == parent mass.
    pub fn verify_mass_conservation(&self, parent_mass: f64) -> bool {
        let total = self.child_mass(parent_mass) * self.n_children as f64;
        (total - parent_mass).abs() < parent_mass * 1e-14
    }
}

// ---------------------------------------------------------------------------
// 3. ParticleMerge — mass-weighted centroid merge
// ---------------------------------------------------------------------------

/// Merges a set of nearby low-density particles into a single parent particle.
///
/// Mass, momentum, and energy are conserved through weighted averaging.
#[derive(Debug, Clone)]
pub struct ParticleMerge {
    /// Maximum number of particles that can be merged in one operation.
    pub max_merge_count: usize,
    /// Maximum allowed mass ratio between particles to be merged.
    pub mass_ratio_limit: f64,
}

impl ParticleMerge {
    /// Create a merge operation manager.
    pub fn new(max_merge_count: usize, mass_ratio_limit: f64) -> Self {
        Self {
            max_merge_count,
            mass_ratio_limit,
        }
    }

    /// Compute mass-weighted centroid position from a list of (mass, position) pairs.
    pub fn merged_position(&self, particles: &[(f64, [f64; 3])]) -> [f64; 3] {
        let total_mass: f64 = particles.iter().map(|(m, _)| m).sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut pos = [0.0f64; 3];
        for (m, p) in particles {
            pos[0] += m * p[0];
            pos[1] += m * p[1];
            pos[2] += m * p[2];
        }
        [
            pos[0] / total_mass,
            pos[1] / total_mass,
            pos[2] / total_mass,
        ]
    }

    /// Compute mass-weighted velocity from a list of (mass, velocity) pairs.
    pub fn merged_velocity(&self, particles: &[(f64, [f64; 3])]) -> [f64; 3] {
        let total_mass: f64 = particles.iter().map(|(m, _)| m).sum();
        if total_mass < 1e-30 {
            return [0.0; 3];
        }
        let mut vel = [0.0f64; 3];
        for (m, v) in particles {
            vel[0] += m * v[0];
            vel[1] += m * v[1];
            vel[2] += m * v[2];
        }
        [
            vel[0] / total_mass,
            vel[1] / total_mass,
            vel[2] / total_mass,
        ]
    }

    /// Total mass of merged particle.
    pub fn merged_mass(&self, particles: &[f64]) -> f64 {
        particles.iter().sum()
    }

    /// Check if a set of particles is eligible for merging (mass ratio within limit).
    pub fn eligible(&self, masses: &[f64]) -> bool {
        if masses.len() < 2 || masses.len() > self.max_merge_count {
            return false;
        }
        let m_min = masses.iter().cloned().fold(f64::INFINITY, f64::min);
        let m_max = masses.iter().cloned().fold(0.0f64, f64::max);
        if m_min < 1e-30 {
            return false;
        }
        m_max / m_min <= self.mass_ratio_limit
    }

    /// Merged smoothing length from volume-weighted average.
    ///
    /// h_merged = (Σ h_i³)^(1/3) (volume additivity in 3-D).
    pub fn merged_smoothing_length(&self, smoothing_lengths: &[f64]) -> f64 {
        let sum_v: f64 = smoothing_lengths.iter().map(|h| h * h * h).sum();
        sum_v.cbrt()
    }
}

// ---------------------------------------------------------------------------
// 4. SmoothingLengthUpdate — variable h: h = h₀ (ρ₀/ρ)^(1/d)
// ---------------------------------------------------------------------------

/// Updates the smoothing length based on local density for variable-resolution SPH.
///
/// Uses the relation **h = h₀ · (ρ₀ / ρ)^(1/d)** derived from constant
/// particle count within the kernel support in d dimensions.
#[derive(Debug, Clone)]
pub struct SmoothingLengthUpdate {
    /// Reference smoothing length h₀.
    pub h0: f64,
    /// Reference density ρ₀.
    pub rho0: f64,
    /// Spatial dimension d (1, 2, or 3).
    pub dimension: u32,
    /// Minimum allowed smoothing length (prevents collapse).
    pub h_min: f64,
    /// Maximum allowed smoothing length (prevents excessive coarsening).
    pub h_max: f64,
}

impl SmoothingLengthUpdate {
    /// Create a smoothing-length updater.
    pub fn new(h0: f64, rho0: f64, dimension: u32) -> Self {
        Self {
            h0,
            rho0,
            dimension,
            h_min: h0 * 0.1,
            h_max: h0 * 10.0,
        }
    }

    /// Compute updated smoothing length for a given density.
    ///
    /// h = h₀ · (ρ₀ / ρ)^(1/d), clamped to \[h_min, h_max\].
    pub fn update(&self, density: f64) -> f64 {
        if density < 1e-30 {
            return self.h_max;
        }
        let h = self.h0 * (self.rho0 / density).powf(1.0 / self.dimension as f64);
        h.clamp(self.h_min, self.h_max)
    }

    /// Iterative update (Springel & Hernquist-style): Newton solver for
    /// h such that (4π/3) h³ n = N_sph, given neighbour count `n_neighbours`.
    ///
    /// Simplified: h_new = h_old · (N_target / N_actual)^(1/d).
    pub fn iterate(&self, h_old: f64, n_actual: usize, n_target: usize) -> f64 {
        if n_actual == 0 {
            return self.h_max;
        }
        let h = h_old * (n_target as f64 / n_actual as f64).powf(1.0 / self.dimension as f64);
        h.clamp(self.h_min, self.h_max)
    }

    /// Ω correction factor for variable-h kernels: Ω = 1 - dh/dρ · ∂W/∂h.
    ///
    /// Approximate: Ω ≈ 1 + (h / (d ρ)) · ∂ρ/∂h (chain rule, first order).
    pub fn omega_correction(&self, density: f64) -> f64 {
        if density < 1e-30 {
            return 1.0;
        }
        let h = self.update(density);
        // dh/dρ = -h / (d ρ)
        let dh_drho = -h / (self.dimension as f64 * density);
        // Ω = 1 - (ρ/h) · dh/dρ (sign convention)
        1.0 - (density / h) * dh_drho
    }
}

// ---------------------------------------------------------------------------
// 5. AdaptiveParticleSet — particle set with split/merge
// ---------------------------------------------------------------------------

/// A single SPH particle with mass, position, velocity, density, smoothing length,
/// and refinement level.
#[derive(Debug, Clone)]
pub struct AdaptiveParticle {
    /// Particle mass (kg).
    pub mass: f64,
    /// Position (m).
    pub position: [f64; 3],
    /// Velocity (m/s).
    pub velocity: [f64; 3],
    /// Density (kg/m³).
    pub density: f64,
    /// Smoothing length (m).
    pub h: f64,
    /// Refinement level (0 = coarsest).
    pub level: u32,
}

impl AdaptiveParticle {
    /// Create a new particle.
    pub fn new(mass: f64, position: [f64; 3], velocity: [f64; 3], density: f64, h: f64) -> Self {
        Self {
            mass,
            position,
            velocity,
            density,
            h,
            level: 0,
        }
    }
}

/// Particle set supporting adaptive splitting and merging operations.
#[derive(Debug, Clone)]
pub struct AdaptiveParticleSet {
    /// All particles in the simulation.
    pub particles: Vec<AdaptiveParticle>,
    /// Maximum refinement level allowed.
    pub max_level: u32,
}

impl AdaptiveParticleSet {
    /// Create an empty particle set.
    pub fn new(max_level: u32) -> Self {
        Self {
            particles: Vec::new(),
            max_level,
        }
    }

    /// Add a particle to the set.
    pub fn add(&mut self, p: AdaptiveParticle) {
        self.particles.push(p);
    }

    /// Total number of particles.
    pub fn len(&self) -> usize {
        self.particles.len()
    }

    /// Returns `true` if the set is empty.
    pub fn is_empty(&self) -> bool {
        self.particles.is_empty()
    }

    /// Total mass (should be conserved through split/merge).
    pub fn total_mass(&self) -> f64 {
        self.particles.iter().map(|p| p.mass).sum()
    }

    /// Total momentum (should be conserved).
    pub fn total_momentum(&self) -> [f64; 3] {
        let mut mom = [0.0f64; 3];
        for p in &self.particles {
            mom[0] += p.mass * p.velocity[0];
            mom[1] += p.mass * p.velocity[1];
            mom[2] += p.mass * p.velocity[2];
        }
        mom
    }

    /// Split particle at index `i` into `n_children` children using a `ParticleSplit`.
    ///
    /// Returns the number of new particles added (n_children - 1).
    pub fn split_particle(&mut self, i: usize, splitter: &ParticleSplit) -> usize {
        if i >= self.particles.len() {
            return 0;
        }
        let parent = self.particles[i].clone();
        if parent.level >= self.max_level {
            return 0;
        }
        let child_mass = splitter.child_mass(parent.mass);
        let child_h = splitter.child_smoothing_length(parent.h);
        let positions = splitter.split_positions_3d(parent.position, parent.h);

        // Replace parent with first child, append the rest
        self.particles[i] = AdaptiveParticle {
            mass: child_mass,
            position: positions[0],
            velocity: parent.velocity,
            density: parent.density,
            h: child_h,
            level: parent.level + 1,
        };
        for pos in positions.iter().skip(1) {
            self.particles.push(AdaptiveParticle {
                mass: child_mass,
                position: *pos,
                velocity: parent.velocity,
                density: parent.density,
                h: child_h,
                level: parent.level + 1,
            });
        }
        splitter.n_children - 1
    }

    /// Merge particles at indices `indices` into a single particle.
    ///
    /// The merged particle replaces the first index; the rest are removed.
    pub fn merge_particles(&mut self, indices: &[usize], merger: &ParticleMerge) -> bool {
        if indices.len() < 2 {
            return false;
        }
        let masses: Vec<f64> = indices.iter().map(|&i| self.particles[i].mass).collect();
        if !merger.eligible(&masses) {
            return false;
        }
        let pos_pairs: Vec<(f64, [f64; 3])> = indices
            .iter()
            .map(|&i| (self.particles[i].mass, self.particles[i].position))
            .collect();
        let vel_pairs: Vec<(f64, [f64; 3])> = indices
            .iter()
            .map(|&i| (self.particles[i].mass, self.particles[i].velocity))
            .collect();
        let hs: Vec<f64> = indices.iter().map(|&i| self.particles[i].h).collect();

        let merged_pos = merger.merged_position(&pos_pairs);
        let merged_vel = merger.merged_velocity(&vel_pairs);
        let merged_mass = merger.merged_mass(&masses);
        let merged_h = merger.merged_smoothing_length(&hs);
        let merged_density = self.particles[indices[0]].density;
        let merged_level = indices
            .iter()
            .map(|&i| self.particles[i].level)
            .min()
            .unwrap_or(0)
            .saturating_sub(1);

        // Update first particle
        self.particles[indices[0]] = AdaptiveParticle {
            mass: merged_mass,
            position: merged_pos,
            velocity: merged_vel,
            density: merged_density,
            h: merged_h,
            level: merged_level,
        };
        // Remove the rest (in reverse order to keep indices valid)
        let mut to_remove: Vec<usize> = indices[1..].to_vec();
        to_remove.sort_unstable();
        for &idx in to_remove.iter().rev() {
            self.particles.swap_remove(idx);
        }
        true
    }
}

// ---------------------------------------------------------------------------
// 6. RefinementZone — forced fine resolution in a region
// ---------------------------------------------------------------------------

/// Shape variants for a refinement zone.
#[derive(Debug, Clone)]
pub enum ZoneShape {
    /// Axis-aligned bounding box.
    Box {
        /// Minimum corner.
        min: [f64; 3],
        /// Maximum corner.
        max: [f64; 3],
    },
    /// Sphere.
    Sphere {
        /// Centre.
        centre: [f64; 3],
        /// Radius.
        radius: f64,
    },
}

/// A region where fine particle resolution is enforced.
#[derive(Debug, Clone)]
pub struct RefinementZone {
    /// Geometric shape of the zone.
    pub shape: ZoneShape,
    /// Target smoothing length inside the zone.
    pub target_h: f64,
    /// Maximum refinement level inside the zone.
    pub max_level: u32,
}

impl RefinementZone {
    /// Create a box-shaped refinement zone.
    pub fn new_box(min: [f64; 3], max: [f64; 3], target_h: f64, max_level: u32) -> Self {
        Self {
            shape: ZoneShape::Box { min, max },
            target_h,
            max_level,
        }
    }

    /// Create a sphere-shaped refinement zone.
    pub fn new_sphere(centre: [f64; 3], radius: f64, target_h: f64, max_level: u32) -> Self {
        Self {
            shape: ZoneShape::Sphere { centre, radius },
            target_h,
            max_level,
        }
    }

    /// Check whether a point is inside the zone.
    pub fn contains(&self, point: [f64; 3]) -> bool {
        match &self.shape {
            ZoneShape::Box { min, max } => {
                point[0] >= min[0]
                    && point[0] <= max[0]
                    && point[1] >= min[1]
                    && point[1] <= max[1]
                    && point[2] >= min[2]
                    && point[2] <= max[2]
            }
            ZoneShape::Sphere { centre, radius } => {
                let dx = point[0] - centre[0];
                let dy = point[1] - centre[1];
                let dz = point[2] - centre[2];
                (dx * dx + dy * dy + dz * dz).sqrt() <= *radius
            }
        }
    }

    /// Returns `true` if a particle at `pos` with smoothing length `h` needs
    /// splitting to meet the zone's resolution requirement.
    pub fn requires_splitting(&self, pos: [f64; 3], h: f64, level: u32) -> bool {
        self.contains(pos) && (h > self.target_h || level < self.max_level)
    }
}

// ---------------------------------------------------------------------------
// 7. VariableResolutionKernel — Ω correction for variable h
// ---------------------------------------------------------------------------

/// Kernel gradient correction (Ω factor) for variable-h SPH.
///
/// The Ω correction accounts for the variation of h with position when
/// computing pressure gradient and viscous forces:
///
/// **corrected gradient = (1/Ω) · ∇W**
#[derive(Debug, Clone)]
pub struct VariableResolutionKernel {
    /// Cubic spline kernel radius multiplier (support = κ h).
    pub kappa: f64,
    /// Spatial dimension.
    pub dim: u32,
}

impl VariableResolutionKernel {
    /// Create a variable-resolution kernel helper.
    pub fn new(kappa: f64, dim: u32) -> Self {
        Self { kappa, dim }
    }

    /// Cubic spline kernel W(q) where q = r / h, normalised for `dim` dimensions.
    pub fn kernel_w(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        let sigma = self.normalization_constant(h);
        if q < 1.0 {
            sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
        } else if q < 2.0 {
            sigma * 0.25 * (2.0 - q).powi(3)
        } else {
            0.0
        }
    }

    /// Kernel gradient magnitude dW/dr (radial derivative).
    pub fn kernel_dw_dr(&self, r: f64, h: f64) -> f64 {
        let q = r / h;
        let sigma = self.normalization_constant(h);
        if r < 1e-30 {
            return 0.0;
        }
        let dw_dq = if q < 1.0 {
            sigma * (-3.0 * q + 2.25 * q * q)
        } else if q < 2.0 {
            sigma * (-0.75 * (2.0 - q).powi(2))
        } else {
            0.0
        };
        dw_dq / h
    }

    /// Ω correction factor: Ω_i = 1 − (dh_i/dρ_i) Σ_j m_j (∂W_ij/∂h_i).
    ///
    /// Simplified approximation: Ω ≈ 1 + (ρ / (d h)) · (dh/dρ · h).
    /// Given constant-neighbour-count smoothing: dh/dρ = −h/(d ρ), so Ω ≈ 1 + 1/d.
    pub fn omega_factor(&self) -> f64 {
        1.0 + 1.0 / self.dim as f64
    }

    /// Apply Ω correction to a kernel gradient component.
    pub fn corrected_gradient(&self, raw_grad: f64) -> f64 {
        raw_grad / self.omega_factor()
    }

    /// Normalisation constant σ for the cubic spline kernel.
    fn normalization_constant(&self, h: f64) -> f64 {
        match self.dim {
            1 => 2.0 / (3.0 * h),
            2 => 10.0 / (7.0 * PI * h * h),
            _ => 1.0 / (PI * h * h * h),
        }
    }
}

// ---------------------------------------------------------------------------
// 8. AdaptiveTimeStep — CFL-based adaptive dt
// ---------------------------------------------------------------------------

/// CFL-based adaptive time step controller for SPH.
///
/// Combines signal velocity (CFL), viscous, and body-force conditions.
#[derive(Debug, Clone)]
pub struct AdaptiveTimeStep {
    /// CFL safety factor α_CFL ∈ (0, 1].
    pub cfl_factor: f64,
    /// Viscous safety factor α_visc ∈ (0, 1].
    pub viscous_factor: f64,
    /// Body-force safety factor α_f ∈ (0, 1].
    pub force_factor: f64,
    /// Maximum allowed time step (s).
    pub dt_max: f64,
    /// Minimum allowed time step (s).
    pub dt_min: f64,
}

impl AdaptiveTimeStep {
    /// Create an adaptive time step controller.
    pub fn new(cfl_factor: f64, viscous_factor: f64, force_factor: f64, dt_max: f64) -> Self {
        Self {
            cfl_factor,
            viscous_factor,
            force_factor,
            dt_max,
            dt_min: dt_max * 1e-8,
        }
    }

    /// CFL condition: dt_cfl = α_CFL · h / v_sig.
    pub fn dt_cfl(&self, h: f64, v_signal: f64) -> f64 {
        if v_signal < 1e-30 {
            return self.dt_max;
        }
        self.cfl_factor * h / v_signal
    }

    /// Viscous condition: dt_visc = α_visc · h² / ν (kinematic viscosity ν).
    pub fn dt_viscous(&self, h: f64, kinematic_viscosity: f64) -> f64 {
        if kinematic_viscosity < 1e-30 {
            return self.dt_max;
        }
        self.viscous_factor * h * h / kinematic_viscosity
    }

    /// Body-force condition: dt_force = α_f · √(h / |a|).
    pub fn dt_force(&self, h: f64, acceleration_magnitude: f64) -> f64 {
        if acceleration_magnitude < 1e-30 {
            return self.dt_max;
        }
        self.force_factor * (h / acceleration_magnitude).sqrt()
    }

    /// Overall adaptive time step: minimum of all conditions.
    pub fn compute_dt(
        &self,
        h: f64,
        v_signal: f64,
        kinematic_viscosity: f64,
        acceleration_magnitude: f64,
    ) -> f64 {
        let dt = self
            .dt_cfl(h, v_signal)
            .min(self.dt_viscous(h, kinematic_viscosity))
            .min(self.dt_force(h, acceleration_magnitude))
            .min(self.dt_max);
        dt.max(self.dt_min)
    }
}

// ---------------------------------------------------------------------------
// 9. ParticleResolution — resolution management
// ---------------------------------------------------------------------------

/// Resolution management: tracks target count, current count, and refinement level.
#[derive(Debug, Clone)]
pub struct ParticleResolution {
    /// Target number of particles in the simulation.
    pub target_count: usize,
    /// Current number of particles.
    pub current_count: usize,
    /// Maximum refinement level in the current set.
    pub max_level: u32,
    /// Global refinement level counter.
    pub refinement_step: u64,
}

impl ParticleResolution {
    /// Create a resolution manager.
    pub fn new(target_count: usize) -> Self {
        Self {
            target_count,
            current_count: 0,
            max_level: 0,
            refinement_step: 0,
        }
    }

    /// Update current count and level from the particle set.
    pub fn update(&mut self, current_count: usize, max_level: u32) {
        self.current_count = current_count;
        self.max_level = max_level;
        self.refinement_step += 1;
    }

    /// Returns `true` if more particles are needed (below target).
    pub fn needs_more_particles(&self) -> bool {
        self.current_count < self.target_count
    }

    /// Returns `true` if particles should be merged (above target).
    pub fn needs_fewer_particles(&self) -> bool {
        self.current_count > self.target_count * 2
    }

    /// Ratio of current to target count.
    pub fn resolution_ratio(&self) -> f64 {
        if self.target_count == 0 {
            return 1.0;
        }
        self.current_count as f64 / self.target_count as f64
    }

    /// Returns `true` when the resolution is within ±20 % of target.
    pub fn is_converged(&self) -> bool {
        let r = self.resolution_ratio();
        (0.8..=1.2).contains(&r)
    }
}

// ---------------------------------------------------------------------------
// 10. MassTransfer — mass/momentum transfer between refined and coarse particles
// ---------------------------------------------------------------------------

/// Mass and momentum transfer operator between levels of refinement.
///
/// Ensures conservative exchange of mass, momentum, and energy when
/// particles at different refinement levels are coupled.
#[derive(Debug, Clone)]
pub struct MassTransfer {
    /// Coupling coefficient α ∈ (0, 1]: fraction of mass transferred per step.
    pub coupling: f64,
    /// True when energy transfer is also applied (in addition to mass/momentum).
    pub transfer_energy: bool,
}

impl MassTransfer {
    /// Create a mass transfer operator.
    pub fn new(coupling: f64, transfer_energy: bool) -> Self {
        Self {
            coupling: coupling.clamp(0.0, 1.0),
            transfer_energy,
        }
    }

    /// Compute the mass transferred from a coarse particle to a fine particle.
    ///
    /// δm = α · (m_c / n_fine) where n_fine is the number of fine particles
    /// overlapping the coarse particle.
    pub fn mass_transfer(&self, coarse_mass: f64, n_fine_neighbours: usize) -> f64 {
        if n_fine_neighbours == 0 {
            return 0.0;
        }
        self.coupling * coarse_mass / n_fine_neighbours as f64
    }

    /// Momentum transferred along with δm from coarse particle at velocity `v_c`.
    pub fn momentum_transfer(&self, delta_mass: f64, coarse_velocity: [f64; 3]) -> [f64; 3] {
        [
            delta_mass * coarse_velocity[0],
            delta_mass * coarse_velocity[1],
            delta_mass * coarse_velocity[2],
        ]
    }

    /// Kinetic energy transferred: δE = ½ δm |v_c|².
    pub fn energy_transfer(&self, delta_mass: f64, coarse_velocity: [f64; 3]) -> f64 {
        if !self.transfer_energy {
            return 0.0;
        }
        let v2 =
            coarse_velocity[0].powi(2) + coarse_velocity[1].powi(2) + coarse_velocity[2].powi(2);
        0.5 * delta_mass * v2
    }

    /// Apply one transfer step: returns updated (fine_mass, fine_velocity).
    pub fn apply(
        &self,
        fine_mass: f64,
        fine_velocity: [f64; 3],
        coarse_mass: f64,
        coarse_velocity: [f64; 3],
        n_fine_neighbours: usize,
    ) -> (f64, [f64; 3]) {
        let dm = self.mass_transfer(coarse_mass, n_fine_neighbours);
        let dp = self.momentum_transfer(dm, coarse_velocity);
        let new_mass = fine_mass + dm;
        if new_mass < 1e-30 {
            return (fine_mass, fine_velocity);
        }
        let new_vel = [
            (fine_mass * fine_velocity[0] + dp[0]) / new_mass,
            (fine_mass * fine_velocity[1] + dp[1]) / new_mass,
            (fine_mass * fine_velocity[2] + dp[2]) / new_mass,
        ];
        (new_mass, new_vel)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ── RefinementCriteria ────────────────────────────────────────────────

    fn make_criteria() -> RefinementCriteria {
        RefinementCriteria::new(10.0, 5.0, 100.0)
    }

    #[test]
    fn test_criteria_needs_refinement_density_grad() {
        let c = make_criteria();
        assert!(c.needs_refinement(15.0, 0.0, 0.0));
    }

    #[test]
    fn test_criteria_needs_refinement_curl() {
        let c = make_criteria();
        assert!(c.needs_refinement(0.0, 8.0, 0.0));
    }

    #[test]
    fn test_criteria_no_refinement_below_thresholds() {
        let c = make_criteria();
        assert!(!c.needs_refinement(5.0, 2.0, 50.0));
    }

    #[test]
    fn test_criteria_can_merge_below_half() {
        let c = make_criteria();
        assert!(c.can_merge(4.0, 1.0, 20.0));
    }

    #[test]
    fn test_criteria_cannot_merge_above_half() {
        let c = make_criteria();
        assert!(!c.can_merge(8.0, 1.0, 20.0));
    }

    #[test]
    fn test_criteria_error_score_above_threshold_gt_one() {
        let c = make_criteria();
        let score = c.error_score(20.0, 0.0, 0.0);
        assert!(score > 1.0);
    }

    #[test]
    fn test_criteria_error_score_below_threshold_lt_one() {
        let c = make_criteria();
        let score = c.error_score(5.0, 2.0, 50.0);
        assert!(score < 1.0);
    }

    // ── ParticleSplit ─────────────────────────────────────────────────────

    fn make_splitter() -> ParticleSplit {
        ParticleSplit::new(4, 0.25)
    }

    #[test]
    fn test_split_child_mass_conserved() {
        let s = make_splitter();
        let pm = 1.0e-3;
        assert!((s.child_mass(pm) * s.n_children as f64 - pm).abs() < 1e-20);
    }

    #[test]
    fn test_split_child_h_smaller() {
        let s = make_splitter();
        let parent_h = 0.01;
        assert!(s.child_smoothing_length(parent_h) < parent_h);
    }

    #[test]
    fn test_split_positions_2d_count() {
        let s = make_splitter();
        let positions = s.split_positions_2d([0.0, 0.0], 0.01);
        assert_eq!(positions.len(), s.n_children);
    }

    #[test]
    fn test_split_positions_3d_count() {
        let s = make_splitter();
        let positions = s.split_positions_3d([0.0; 3], 0.01);
        assert_eq!(positions.len(), s.n_children);
    }

    #[test]
    fn test_split_mass_conservation_verified() {
        let s = make_splitter();
        assert!(s.verify_mass_conservation(1.0));
    }

    // ── ParticleMerge ─────────────────────────────────────────────────────

    fn make_merger() -> ParticleMerge {
        ParticleMerge::new(4, 2.0)
    }

    #[test]
    fn test_merge_centroid_weighted() {
        let m = make_merger();
        let particles = vec![(1.0, [0.0, 0.0, 0.0]), (1.0, [2.0, 0.0, 0.0])];
        let pos = m.merged_position(&particles);
        assert!((pos[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_merge_velocity_momentum_conserving() {
        let m = make_merger();
        let particles = vec![(1.0, [1.0, 0.0, 0.0]), (3.0, [3.0, 0.0, 0.0])];
        let vel = m.merged_velocity(&particles);
        // mass-weighted: (1*1 + 3*3) / 4 = 10/4 = 2.5
        assert!((vel[0] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_merge_total_mass() {
        let m = make_merger();
        let total = m.merged_mass(&[1.0, 2.0, 3.0]);
        assert!((total - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_merge_eligible_equal_masses() {
        let m = make_merger();
        assert!(m.eligible(&[1.0, 1.0, 1.0]));
    }

    #[test]
    fn test_merge_not_eligible_large_ratio() {
        let m = make_merger();
        assert!(!m.eligible(&[1.0, 5.0]));
    }

    #[test]
    fn test_merge_smoothing_length_volume_additive() {
        let m = make_merger();
        let hs = [1.0f64, 1.0, 1.0, 1.0];
        let h_merged = m.merged_smoothing_length(&hs);
        // (4 * 1.0)^(1/3) ≈ 1.587
        assert!((h_merged - 4.0f64.cbrt()).abs() < 1e-10);
    }

    // ── SmoothingLengthUpdate ─────────────────────────────────────────────

    fn make_slu() -> SmoothingLengthUpdate {
        SmoothingLengthUpdate::new(0.01, 1000.0, 3)
    }

    #[test]
    fn test_slu_update_at_reference_density() {
        let slu = make_slu();
        let h = slu.update(1000.0);
        assert!((h - 0.01).abs() < 1e-12);
    }

    #[test]
    fn test_slu_update_denser_gives_smaller_h() {
        let slu = make_slu();
        let h = slu.update(2000.0);
        assert!(h < 0.01);
    }

    #[test]
    fn test_slu_update_clamped_to_h_max() {
        let slu = make_slu();
        let h = slu.update(1e-10);
        assert!(h <= slu.h_max);
    }

    #[test]
    fn test_slu_iterate_increases_h_when_too_few_neighbours() {
        let slu = make_slu();
        let h_new = slu.iterate(0.01, 10, 30);
        assert!(h_new > 0.01);
    }

    #[test]
    fn test_slu_omega_correction_above_one() {
        let slu = make_slu();
        let omega = slu.omega_correction(1000.0);
        assert!(omega > 1.0, "Ω should be > 1 for constant-N smoothing");
    }

    // ── AdaptiveParticleSet ───────────────────────────────────────────────

    fn make_set() -> AdaptiveParticleSet {
        let mut set = AdaptiveParticleSet::new(3);
        set.add(AdaptiveParticle::new(
            1.0e-3,
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            1000.0,
            0.01,
        ));
        set
    }

    #[test]
    fn test_set_total_mass_single_particle() {
        let set = make_set();
        assert!((set.total_mass() - 1.0e-3).abs() < 1e-20);
    }

    #[test]
    fn test_set_split_increases_count() {
        let mut set = make_set();
        let splitter = ParticleSplit::new(4, 0.25);
        let added = set.split_particle(0, &splitter);
        assert_eq!(added, 3);
        assert_eq!(set.len(), 4);
    }

    #[test]
    fn test_set_split_conserves_mass() {
        let mut set = make_set();
        let mass_before = set.total_mass();
        let splitter = ParticleSplit::new(4, 0.25);
        set.split_particle(0, &splitter);
        let mass_after = set.total_mass();
        assert!(
            (mass_before - mass_after).abs() < 1e-18,
            "mass not conserved"
        );
    }

    #[test]
    fn test_set_merge_reduces_count() {
        let mut set = AdaptiveParticleSet::new(3);
        for i in 0..3 {
            set.add(AdaptiveParticle::new(
                1.0e-3,
                [i as f64 * 0.005, 0.0, 0.0],
                [0.0; 3],
                1000.0,
                0.01,
            ));
        }
        let merger = ParticleMerge::new(4, 2.0);
        let ok = set.merge_particles(&[0, 1, 2], &merger);
        assert!(ok);
        assert_eq!(set.len(), 1);
    }

    // ── RefinementZone ────────────────────────────────────────────────────

    #[test]
    fn test_zone_box_contains_interior_point() {
        let zone = RefinementZone::new_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.005, 3);
        assert!(zone.contains([0.5, 0.5, 0.5]));
    }

    #[test]
    fn test_zone_box_does_not_contain_exterior() {
        let zone = RefinementZone::new_box([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], 0.005, 3);
        assert!(!zone.contains([1.5, 0.5, 0.5]));
    }

    #[test]
    fn test_zone_sphere_contains_centre() {
        let zone = RefinementZone::new_sphere([0.0; 3], 0.5, 0.005, 3);
        assert!(zone.contains([0.0, 0.0, 0.0]));
    }

    #[test]
    fn test_zone_sphere_does_not_contain_far_point() {
        let zone = RefinementZone::new_sphere([0.0; 3], 0.5, 0.005, 3);
        assert!(!zone.contains([1.0, 0.0, 0.0]));
    }

    #[test]
    fn test_zone_requires_splitting_large_h() {
        let zone = RefinementZone::new_box([0.0; 3], [1.0; 3], 0.005, 3);
        assert!(zone.requires_splitting([0.5, 0.5, 0.5], 0.01, 0));
    }

    // ── VariableResolutionKernel ──────────────────────────────────────────

    fn make_vrk() -> VariableResolutionKernel {
        VariableResolutionKernel::new(2.0, 3)
    }

    #[test]
    fn test_vrk_kernel_positive_inside_support() {
        let k = make_vrk();
        assert!(k.kernel_w(0.5, 1.0) > 0.0);
    }

    #[test]
    fn test_vrk_kernel_zero_outside_support() {
        let k = make_vrk();
        assert_eq!(k.kernel_w(3.0, 1.0), 0.0);
    }

    #[test]
    fn test_vrk_dw_dr_negative_for_r_gt_zero() {
        let k = make_vrk();
        // Kernel decreases away from origin → dW/dr < 0
        assert!(k.kernel_dw_dr(0.5, 1.0) < 0.0);
    }

    #[test]
    fn test_vrk_omega_factor_greater_than_one() {
        let k = make_vrk();
        assert!(k.omega_factor() > 1.0);
    }

    #[test]
    fn test_vrk_corrected_gradient_less_than_raw() {
        let k = make_vrk();
        let raw = 10.0;
        let corrected = k.corrected_gradient(raw);
        assert!(corrected < raw);
    }

    // ── AdaptiveTimeStep ──────────────────────────────────────────────────

    fn make_dt() -> AdaptiveTimeStep {
        AdaptiveTimeStep::new(0.3, 0.125, 0.25, 1.0e-3)
    }

    #[test]
    fn test_dt_cfl_decreases_with_velocity() {
        let dt = make_dt();
        let dt1 = dt.dt_cfl(0.01, 100.0);
        let dt2 = dt.dt_cfl(0.01, 200.0);
        assert!(dt1 > dt2);
    }

    #[test]
    fn test_dt_viscous_decreases_with_viscosity() {
        let dt = make_dt();
        let dt1 = dt.dt_viscous(0.01, 1.0e-6);
        let dt2 = dt.dt_viscous(0.01, 1.0e-4);
        assert!(dt1 > dt2);
    }

    #[test]
    fn test_dt_compute_bounded_by_dt_max() {
        let dt = make_dt();
        let result = dt.compute_dt(0.001, 1.0, 1.0e-6, 9.8);
        assert!(result <= dt.dt_max);
    }

    #[test]
    fn test_dt_compute_bounded_by_dt_min() {
        let dt = make_dt();
        let result = dt.compute_dt(1.0e-10, 1.0e10, 1.0, 1.0e10);
        assert!(result >= dt.dt_min);
    }

    // ── ParticleResolution ────────────────────────────────────────────────

    #[test]
    fn test_resolution_needs_more_when_below_target() {
        let mut r = ParticleResolution::new(1000);
        r.update(500, 2);
        assert!(r.needs_more_particles());
    }

    #[test]
    fn test_resolution_needs_fewer_when_above_twice_target() {
        let mut r = ParticleResolution::new(1000);
        r.update(2500, 5);
        assert!(r.needs_fewer_particles());
    }

    #[test]
    fn test_resolution_converged_near_target() {
        let mut r = ParticleResolution::new(1000);
        r.update(950, 2);
        assert!(r.is_converged());
    }

    #[test]
    fn test_resolution_ratio_correct() {
        let mut r = ParticleResolution::new(1000);
        r.update(500, 2);
        assert!((r.resolution_ratio() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_resolution_step_increments() {
        let mut r = ParticleResolution::new(1000);
        r.update(1000, 1);
        r.update(1000, 1);
        assert_eq!(r.refinement_step, 2);
    }

    // ── MassTransfer ─────────────────────────────────────────────────────

    fn make_transfer() -> MassTransfer {
        MassTransfer::new(0.1, true)
    }

    #[test]
    fn test_mass_transfer_amount_proportional_to_coupling() {
        let mt = make_transfer();
        let dm = mt.mass_transfer(1.0, 1);
        assert!((dm - 0.1).abs() < 1e-12);
    }

    #[test]
    fn test_mass_transfer_zero_when_no_neighbours() {
        let mt = make_transfer();
        assert_eq!(mt.mass_transfer(1.0, 0), 0.0);
    }

    #[test]
    fn test_momentum_transfer_proportional_to_dm() {
        let mt = make_transfer();
        let dp = mt.momentum_transfer(0.1, [2.0, 0.0, 0.0]);
        assert!((dp[0] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn test_energy_transfer_positive() {
        let mt = make_transfer();
        let e = mt.energy_transfer(0.1, [10.0, 0.0, 0.0]);
        assert!(e > 0.0);
    }

    #[test]
    fn test_energy_transfer_zero_when_disabled() {
        let mt = MassTransfer::new(0.1, false);
        let e = mt.energy_transfer(0.1, [10.0, 0.0, 0.0]);
        assert_eq!(e, 0.0);
    }

    #[test]
    fn test_apply_increases_fine_mass() {
        let mt = make_transfer();
        let (new_mass, _) = mt.apply(
            1.0e-3,
            [0.0; 3], // fine
            1.0e-2,
            [1.0, 0.0, 0.0], // coarse
            1,
        );
        assert!(new_mass > 1.0e-3);
    }

    #[test]
    fn test_apply_momentum_conserving() {
        let mt = make_transfer();
        let v_c = [3.0, 0.0, 0.0];
        let v_f = [1.0, 0.0, 0.0];
        let m_f = 1.0;
        let m_c = 1.0;
        let (new_mass, new_vel) = mt.apply(m_f, v_f, m_c, v_c, 1);
        // momentum before: m_f*v_f + dm*v_c; after: new_mass * new_vel
        let dm = mt.mass_transfer(m_c, 1);
        let p_before = m_f * v_f[0] + dm * v_c[0];
        let p_after = new_mass * new_vel[0];
        assert!(
            (p_before - p_after).abs() < 1e-12,
            "momentum not conserved: {p_before} vs {p_after}"
        );
    }
}
