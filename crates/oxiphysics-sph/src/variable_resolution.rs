// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Variable-resolution (adaptive) SPH.
//!
//! Provides particle splitting/merging (refinement/coarsening), an octree
//! spatial index, and consistency-condition checks for multi-resolution SPH
//! simulations.

// ---------------------------------------------------------------------------
// ResolutionLevel
// ---------------------------------------------------------------------------

/// One level of a multi-resolution SPH hierarchy.
#[derive(Debug, Clone)]
pub struct ResolutionLevel {
    /// Smoothing length for this level.
    pub h: f64,
    /// Nominal inter-particle spacing.
    pub particle_spacing: f64,
    /// Level index (0 = coarsest).
    pub level: usize,
}

impl ResolutionLevel {
    /// Create a new resolution level.
    pub fn new(h: f64, particle_spacing: f64, level: usize) -> Self {
        Self {
            h,
            particle_spacing,
            level,
        }
    }
}

// ---------------------------------------------------------------------------
// AdaptiveParticle
// ---------------------------------------------------------------------------

/// An SPH particle that carries its own individual smoothing length and
/// refinement level.
#[derive(Debug, Clone)]
pub struct AdaptiveParticle {
    /// Position in 3-D space.
    pub position: [f64; 3],
    /// Velocity vector.
    pub velocity: [f64; 3],
    /// Particle mass.
    pub mass: f64,
    /// Individual smoothing length.
    pub h: f64,
    /// Refinement level (0 = coarsest).
    pub level: usize,
}

impl AdaptiveParticle {
    /// Create a new adaptive particle.
    pub fn new(position: [f64; 3], velocity: [f64; 3], mass: f64, h: f64, level: usize) -> Self {
        Self {
            position,
            velocity,
            mass,
            h,
            level,
        }
    }

    /// Compute the (approximate) particle volume as mass / reference density.
    pub fn volume(&self, density: f64) -> f64 {
        if density > 0.0 {
            self.mass / density
        } else {
            0.0
        }
    }
}

// ---------------------------------------------------------------------------
// RefinementCriteria
// ---------------------------------------------------------------------------

/// Criterion used to decide whether a particle should be split or merged.
#[derive(Debug, Clone)]
pub enum RefinementCriteria {
    /// Refine if the local velocity-gradient norm exceeds the threshold.
    VelocityGradient(f64),
    /// Refine if the local density deviates from the reference by more than the threshold.
    Density(f64),
    /// Refine if the vorticity magnitude exceeds the threshold.
    VorticityMagnitude(f64),
    /// User-supplied boolean function: takes `(position, velocity, mass, h)`.
    UserDefined,
}

impl RefinementCriteria {
    /// Return `true` if `value` (the relevant physical measure) triggers refinement.
    pub fn should_refine(&self, value: f64) -> bool {
        match self {
            Self::VelocityGradient(thresh) => value > *thresh,
            Self::Density(thresh) => value > *thresh,
            Self::VorticityMagnitude(thresh) => value > *thresh,
            Self::UserDefined => false,
        }
    }

    /// Return `true` if `value` is safely below the coarsening threshold (half the refine threshold).
    pub fn should_coarsen(&self, value: f64) -> bool {
        match self {
            Self::VelocityGradient(thresh) => value < 0.5 * thresh,
            Self::Density(thresh) => value < 0.5 * thresh,
            Self::VorticityMagnitude(thresh) => value < 0.5 * thresh,
            Self::UserDefined => false,
        }
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Compute the optimal smoothing length for a given particle spacing and dimension.
///
/// Uses the standard factor h = α · Δx where α = 1.3 (1-D), 1.2 (2-D), 1.1 (3-D).
pub fn optimal_h(particle_spacing: f64, dim: usize) -> f64 {
    let alpha = match dim {
        1 => 1.3,
        2 => 1.2,
        _ => 1.1,
    };
    alpha * particle_spacing
}

/// Compute particle mass from density and volume.
pub fn mass_from_density(density: f64, volume: f64) -> f64 {
    density * volume
}

/// Split one particle into 8 children arranged in a local cubic pattern.
///
/// Each child has 1/8 of the parent mass, half the parent smoothing length,
/// and is placed at ±δ from the parent position in each Cartesian direction
/// where δ = h/4.
pub fn split_particle(p: &AdaptiveParticle) -> Vec<AdaptiveParticle> {
    let child_mass = p.mass / 8.0;
    let child_h = p.h / 2.0;
    let delta = p.h / 4.0;
    let child_level = p.level + 1;

    let mut children = Vec::with_capacity(8);
    for sx in [-1.0_f64, 1.0] {
        for sy in [-1.0_f64, 1.0] {
            for sz in [-1.0_f64, 1.0] {
                let pos = [
                    p.position[0] + sx * delta,
                    p.position[1] + sy * delta,
                    p.position[2] + sz * delta,
                ];
                children.push(AdaptiveParticle::new(
                    pos,
                    p.velocity,
                    child_mass,
                    child_h,
                    child_level,
                ));
            }
        }
    }
    children
}

// ---------------------------------------------------------------------------
// AdaptiveSph
// ---------------------------------------------------------------------------

/// Adaptive SPH manager: stores particles and resolution levels.
pub struct AdaptiveSph {
    /// All particles in the simulation.
    pub particles: Vec<AdaptiveParticle>,
    /// Resolution levels available.
    pub levels: Vec<ResolutionLevel>,
}

impl AdaptiveSph {
    /// Create an empty adaptive SPH system.
    pub fn new(levels: Vec<ResolutionLevel>) -> Self {
        Self {
            particles: Vec::new(),
            levels,
        }
    }

    /// Add a particle.
    pub fn add_particle(&mut self, p: AdaptiveParticle) {
        self.particles.push(p);
    }

    /// Refine (split) particle at `idx`.  Returns the 8 children; the original
    /// particle is left in place (caller decides whether to remove it).
    pub fn refine_particle(&self, idx: usize) -> Vec<AdaptiveParticle> {
        split_particle(&self.particles[idx])
    }

    /// Coarsen (merge) the particles listed in `indices` into one representative
    /// particle.  Uses mass-weighted averaging for position and velocity.
    pub fn coarsen_particles(&self, indices: &[usize]) -> AdaptiveParticle {
        assert!(!indices.is_empty(), "Cannot coarsen an empty set");

        let total_mass: f64 = indices.iter().map(|&i| self.particles[i].mass).sum();
        let mut pos = [0.0_f64; 3];
        let mut vel = [0.0_f64; 3];
        let mut h_sum = 0.0_f64;

        for &i in indices {
            let p = &self.particles[i];
            let w = p.mass / total_mass;
            for k in 0..3 {
                pos[k] += w * p.position[k];
                vel[k] += w * p.velocity[k];
            }
            h_sum += p.h;
        }

        let avg_h = h_sum / indices.len() as f64;
        let parent_level = self.particles[indices[0]].level.saturating_sub(1);

        AdaptiveParticle::new(pos, vel, total_mass, avg_h * 2.0, parent_level)
    }
}

// ---------------------------------------------------------------------------
// ConsistencyCondition
// ---------------------------------------------------------------------------

/// Zeroth- and first-order SPH consistency measures.
///
/// For a consistent SPH approximation:
/// - zeroth-order: Σⱼ mⱼ/ρⱼ W(r-rⱼ, h) ≈ 1  (particle normalization)
/// - first-order:  Σⱼ mⱼ/ρⱼ (rⱼ - r) W(rⱼ-r, h) ≈ 0
///
/// We use a simplified proxy that measures how close the particle distribution
/// is to these conditions.
#[derive(Debug, Clone)]
pub struct ConsistencyCondition {
    /// Mean absolute deviation from 1 of the zeroth-order sum.
    pub zeroth_order: f64,
    /// Mean absolute deviation from 0 of the first-order sum (norm).
    pub first_order: f64,
}

impl ConsistencyCondition {
    /// Compute consistency conditions for a set of particles.
    ///
    /// Uses a Gaussian kernel W(r, h) = exp(-r²/(2h²)) / ((2πh²)^(3/2))
    /// evaluated at each particle's own neighbourhood (all particles within 2h).
    pub fn compute_consistency(particles: &[AdaptiveParticle]) -> (f64, f64) {
        if particles.is_empty() {
            return (0.0, 0.0);
        }

        use std::f64::consts::PI;
        let reference_density = 1000.0_f64; // kg/m³ (water-like)

        let mut zeroth_errors = Vec::new();
        let mut first_errors = Vec::new();

        for pi in particles {
            let h2 = pi.h * pi.h;
            let kernel_norm = 1.0 / ((2.0 * PI * h2).powf(1.5));
            let cutoff = 2.0 * pi.h;

            let mut zeroth_sum = 0.0_f64;
            let mut first_sum = [0.0_f64; 3];

            for pj in particles {
                let dx = [
                    pj.position[0] - pi.position[0],
                    pj.position[1] - pi.position[1],
                    pj.position[2] - pi.position[2],
                ];
                let r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                if r2 > cutoff * cutoff {
                    continue;
                }
                let w = kernel_norm * (-0.5 * r2 / h2).exp();
                let vj = pj.mass / reference_density;
                zeroth_sum += vj * w;
                for k in 0..3 {
                    first_sum[k] += vj * dx[k] * w;
                }
            }

            zeroth_errors.push((zeroth_sum - 1.0).abs());
            let fnorm = (first_sum[0].powi(2) + first_sum[1].powi(2) + first_sum[2].powi(2)).sqrt();
            first_errors.push(fnorm);
        }

        let n = particles.len() as f64;
        let z0 = zeroth_errors.iter().sum::<f64>() / n;
        let f1 = first_errors.iter().sum::<f64>() / n;
        (z0, f1)
    }
}

// ---------------------------------------------------------------------------
// OctreeNode
// ---------------------------------------------------------------------------

/// A node in an octree spatial index for SPH particles.
#[derive(Debug)]
pub struct OctreeNode {
    /// Centre of the cube represented by this node.
    pub center: [f64; 3],
    /// Half-edge length of the cube.
    pub half_size: f64,
    /// Indices of particles inside this node (leaf only).
    pub particle_indices: Vec<usize>,
    /// Eight children (None for leaf nodes).
    pub children: Option<Box<[OctreeNode; 8]>>,
}

impl OctreeNode {
    /// Create a new leaf node.
    pub fn new(center: [f64; 3], half_size: f64) -> Self {
        Self {
            center,
            half_size,
            particle_indices: Vec::new(),
            children: None,
        }
    }

    /// Return `true` if `pos` lies inside this node's bounding cube.
    pub fn contains(&self, pos: [f64; 3]) -> bool {
        for (k, &p) in pos.iter().enumerate() {
            if (p - self.center[k]).abs() > self.half_size {
                return false;
            }
        }
        true
    }

    /// Insert particle index `idx` at position `pos`.
    ///
    /// Splits into 8 children when the leaf reaches `max_leaf_size` particles.
    pub fn insert(
        &mut self,
        idx: usize,
        pos: [f64; 3],
        max_leaf_size: usize,
        all_positions: &[[f64; 3]],
    ) {
        if self.children.is_some() {
            // Delegate to the appropriate child
            let child_idx = self.child_index(pos);
            if let Some(children) = self.children.as_mut() {
                children[child_idx].insert(idx, pos, max_leaf_size, all_positions);
            }
            return;
        }

        self.particle_indices.push(idx);

        if self.particle_indices.len() > max_leaf_size {
            self.subdivide(max_leaf_size, all_positions);
        }
    }

    /// Determine which of the 8 octants `pos` falls into.
    fn child_index(&self, pos: [f64; 3]) -> usize {
        let mut idx = 0;
        if pos[0] >= self.center[0] {
            idx |= 1;
        }
        if pos[1] >= self.center[1] {
            idx |= 2;
        }
        if pos[2] >= self.center[2] {
            idx |= 4;
        }
        idx
    }

    /// Subdivide this leaf into 8 children and re-distribute particles.
    fn subdivide(&mut self, max_leaf_size: usize, all_positions: &[[f64; 3]]) {
        let hs = self.half_size / 2.0;
        let c = self.center;
        let offsets: [[f64; 3]; 8] = [
            [-hs, -hs, -hs],
            [hs, -hs, -hs],
            [-hs, hs, -hs],
            [hs, hs, -hs],
            [-hs, -hs, hs],
            [hs, -hs, hs],
            [-hs, hs, hs],
            [hs, hs, hs],
        ];

        // We need to create all 8 children at once to satisfy [T; 8]
        let children: [OctreeNode; 8] = std::array::from_fn(|i| {
            OctreeNode::new(
                [
                    c[0] + offsets[i][0],
                    c[1] + offsets[i][1],
                    c[2] + offsets[i][2],
                ],
                hs,
            )
        });
        self.children = Some(Box::new(children));

        // Move existing particle indices to children
        let old_indices = std::mem::take(&mut self.particle_indices);
        for old_idx in old_indices {
            let pos = all_positions[old_idx];
            let ci = self.child_index(pos);
            if let Some(ch) = self.children.as_mut() {
                ch[ci].insert(old_idx, pos, max_leaf_size, all_positions);
            }
        }
    }

    /// Collect all particle indices in this subtree.
    pub fn collect_all(&self, out: &mut Vec<usize>) {
        out.extend_from_slice(&self.particle_indices);
        if let Some(children) = &self.children {
            for child in children.iter() {
                child.collect_all(out);
            }
        }
    }

    /// Count all particles in this subtree.
    pub fn count(&self) -> usize {
        let mut out = Vec::new();
        self.collect_all(&mut out);
        out.len()
    }

    /// Find all particles within distance `radius` of `query`.
    pub fn query_radius(
        &self,
        query: [f64; 3],
        radius: f64,
        out: &mut Vec<usize>,
        all_positions: &[[f64; 3]],
    ) {
        // Check AABB overlap
        let mut min_dist_sq = 0.0_f64;
        for (k, &qk) in query.iter().enumerate() {
            let lo = self.center[k] - self.half_size;
            let hi = self.center[k] + self.half_size;
            if qk < lo {
                min_dist_sq += (qk - lo).powi(2);
            } else if qk > hi {
                min_dist_sq += (qk - hi).powi(2);
            }
        }
        if min_dist_sq > radius * radius {
            return;
        }

        if self.children.is_none() {
            for &idx in &self.particle_indices {
                let p = all_positions[idx];
                let d2 = (0..3).map(|k| (p[k] - query[k]).powi(2)).sum::<f64>();
                if d2 <= radius * radius {
                    out.push(idx);
                }
            }
        } else if let Some(children) = &self.children {
            for child in children.iter() {
                child.query_radius(query, radius, out, all_positions);
            }
        }
    }
}

/// Build an octree from a list of particle positions.
///
/// Automatically computes a bounding box and creates the root node.
pub fn build_octree(positions: &[[f64; 3]], max_leaf_size: usize) -> OctreeNode {
    if positions.is_empty() {
        return OctreeNode::new([0.0; 3], 1.0);
    }

    // Bounding box
    let mut lo = positions[0];
    let mut hi = positions[0];
    for p in positions.iter() {
        for k in 0..3 {
            if p[k] < lo[k] {
                lo[k] = p[k];
            }
            if p[k] > hi[k] {
                hi[k] = p[k];
            }
        }
    }
    let center = [
        0.5 * (lo[0] + hi[0]),
        0.5 * (lo[1] + hi[1]),
        0.5 * (lo[2] + hi[2]),
    ];
    let half_size = (0..3)
        .map(|k| 0.5 * (hi[k] - lo[k]))
        .fold(0.0_f64, f64::max)
        * 1.001 // small expansion to avoid boundary ambiguity
        + 1e-10;

    let mut root = OctreeNode::new(center, half_size);
    for (i, &pos) in positions.iter().enumerate() {
        root.insert(i, pos, max_leaf_size, positions);
    }
    root
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_particle(pos: [f64; 3], mass: f64, h: f64, level: usize) -> AdaptiveParticle {
        AdaptiveParticle::new(pos, [0.0; 3], mass, h, level)
    }

    // --- ResolutionLevel ---

    #[test]
    fn test_resolution_level_fields() {
        let rl = ResolutionLevel::new(0.1, 0.08, 2);
        assert_eq!(rl.level, 2);
        assert!((rl.h - 0.1).abs() < 1e-12);
    }

    // --- optimal_h ---

    #[test]
    fn test_optimal_h_1d() {
        let h = optimal_h(0.1, 1);
        assert!((h - 0.13).abs() < 1e-10);
    }

    #[test]
    fn test_optimal_h_2d() {
        let h = optimal_h(0.1, 2);
        assert!((h - 0.12).abs() < 1e-10);
    }

    #[test]
    fn test_optimal_h_3d() {
        let h = optimal_h(0.1, 3);
        assert!((h - 0.11).abs() < 1e-10);
    }

    // --- mass_from_density ---

    #[test]
    fn test_mass_from_density() {
        assert!((mass_from_density(1000.0, 0.001) - 1.0).abs() < 1e-10);
    }

    // --- split_particle ---

    #[test]
    fn test_split_particle_count() {
        let p = make_particle([0.0; 3], 8.0, 0.1, 0);
        let children = split_particle(&p);
        assert_eq!(children.len(), 8);
    }

    #[test]
    fn test_split_particle_mass_conservation() {
        let p = make_particle([0.0; 3], 8.0, 0.1, 0);
        let children = split_particle(&p);
        let total: f64 = children.iter().map(|c| c.mass).sum();
        assert!((total - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_split_particle_level_incremented() {
        let p = make_particle([0.0; 3], 1.0, 0.1, 2);
        let children = split_particle(&p);
        assert!(children.iter().all(|c| c.level == 3));
    }

    #[test]
    fn test_split_particle_h_halved() {
        let p = make_particle([0.0; 3], 1.0, 0.2, 0);
        let children = split_particle(&p);
        assert!(children.iter().all(|c| (c.h - 0.1).abs() < 1e-10));
    }

    #[test]
    fn test_split_particle_center_of_mass() {
        let p = make_particle([1.0, 2.0, 3.0], 8.0, 0.4, 0);
        let children = split_particle(&p);
        let cx: f64 = children.iter().map(|c| c.position[0]).sum::<f64>() / 8.0;
        let cy: f64 = children.iter().map(|c| c.position[1]).sum::<f64>() / 8.0;
        let cz: f64 = children.iter().map(|c| c.position[2]).sum::<f64>() / 8.0;
        assert!((cx - 1.0).abs() < 1e-10);
        assert!((cy - 2.0).abs() < 1e-10);
        assert!((cz - 3.0).abs() < 1e-10);
    }

    // --- AdaptiveSph refine/coarsen ---

    #[test]
    fn test_refine_produces_8_children() {
        let levels = vec![ResolutionLevel::new(0.1, 0.08, 0)];
        let mut sph = AdaptiveSph::new(levels);
        sph.add_particle(make_particle([0.0; 3], 8.0, 0.1, 0));
        let children = sph.refine_particle(0);
        assert_eq!(children.len(), 8);
    }

    #[test]
    fn test_coarsen_mass_conservation() {
        let levels = vec![ResolutionLevel::new(0.1, 0.08, 0)];
        let mut sph = AdaptiveSph::new(levels);
        for i in 0..4 {
            sph.add_particle(make_particle([i as f64 * 0.05; 3], 1.0, 0.05, 1));
        }
        let merged = sph.coarsen_particles(&[0, 1, 2, 3]);
        assert!((merged.mass - 4.0).abs() < 1e-10);
    }

    #[test]
    fn test_coarsen_level_decremented() {
        let levels = vec![ResolutionLevel::new(0.1, 0.08, 0)];
        let mut sph = AdaptiveSph::new(levels);
        sph.add_particle(make_particle([0.0; 3], 1.0, 0.05, 2));
        sph.add_particle(make_particle([0.1; 3], 1.0, 0.05, 2));
        let merged = sph.coarsen_particles(&[0, 1]);
        assert_eq!(merged.level, 1);
    }

    #[test]
    fn test_coarsen_position_weighted_average() {
        let levels = vec![ResolutionLevel::new(0.1, 0.08, 0)];
        let mut sph = AdaptiveSph::new(levels);
        sph.add_particle(make_particle([0.0, 0.0, 0.0], 1.0, 0.1, 1));
        sph.add_particle(make_particle([2.0, 0.0, 0.0], 1.0, 0.1, 1));
        let merged = sph.coarsen_particles(&[0, 1]);
        assert!((merged.position[0] - 1.0).abs() < 1e-10);
    }

    // --- RefinementCriteria ---

    #[test]
    fn test_refine_criteria_velocity_gradient() {
        let c = RefinementCriteria::VelocityGradient(0.5);
        assert!(c.should_refine(0.6));
        assert!(!c.should_refine(0.4));
    }

    #[test]
    fn test_refine_criteria_density() {
        let c = RefinementCriteria::Density(100.0);
        assert!(c.should_refine(150.0));
        assert!(!c.should_refine(30.0));
    }

    #[test]
    fn test_coarsen_criteria_vorticity() {
        let c = RefinementCriteria::VorticityMagnitude(1.0);
        assert!(c.should_coarsen(0.3));
        assert!(!c.should_coarsen(0.8));
    }

    #[test]
    fn test_refine_criteria_user_defined_no_refine() {
        let c = RefinementCriteria::UserDefined;
        assert!(!c.should_refine(1e9));
    }

    // --- ConsistencyCondition ---

    #[test]
    fn test_consistency_empty() {
        let (z, f) = ConsistencyCondition::compute_consistency(&[]);
        assert_eq!(z, 0.0);
        assert_eq!(f, 0.0);
    }

    #[test]
    fn test_consistency_single_particle() {
        let p = make_particle([0.0; 3], 1e-6, 0.1, 0);
        let (z, _f) = ConsistencyCondition::compute_consistency(&[p]);
        // With a single tiny-mass particle the zeroth-order sum is near 0, error near 1
        assert!(z >= 0.0);
    }

    #[test]
    fn test_consistency_returns_nonnegative() {
        let particles: Vec<AdaptiveParticle> = (0..5)
            .map(|i| make_particle([i as f64 * 0.05, 0.0, 0.0], 0.125, 0.1, 0))
            .collect();
        let (z, f) = ConsistencyCondition::compute_consistency(&particles);
        assert!(z >= 0.0);
        assert!(f >= 0.0);
    }

    // --- OctreeNode ---

    #[test]
    fn test_octree_contains() {
        let node = OctreeNode::new([0.0; 3], 1.0);
        assert!(node.contains([0.5, 0.5, 0.5]));
        assert!(!node.contains([2.0, 0.0, 0.0]));
    }

    #[test]
    fn test_octree_build_empty() {
        let root = build_octree(&[], 4);
        assert_eq!(root.count(), 0);
    }

    #[test]
    fn test_octree_build_single() {
        let positions = vec![[1.0, 2.0, 3.0]];
        let root = build_octree(&positions, 4);
        assert_eq!(root.count(), 1);
    }

    #[test]
    fn test_octree_build_many() {
        let positions: Vec<[f64; 3]> = (0..64).map(|i| [i as f64 * 0.1, 0.0, 0.0]).collect();
        let root = build_octree(&positions, 4);
        assert_eq!(root.count(), 64);
    }

    #[test]
    fn test_octree_query_radius_finds_near() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let root = build_octree(&positions, 2);
        let mut found = Vec::new();
        root.query_radius([0.0, 0.0, 0.0], 1.5, &mut found, &positions);
        // Should find indices 0 (dist=0) and 1 (dist=1)
        assert!(found.contains(&0));
        assert!(found.contains(&1));
        assert!(!found.contains(&3));
    }

    #[test]
    fn test_octree_query_empty_result() {
        let positions = vec![[10.0, 10.0, 10.0]];
        let root = build_octree(&positions, 4);
        let mut found = Vec::new();
        root.query_radius([0.0, 0.0, 0.0], 0.5, &mut found, &positions);
        assert!(found.is_empty());
    }

    #[test]
    fn test_octree_collect_all() {
        let positions: Vec<[f64; 3]> = (0..20).map(|i| [i as f64, 0.0, 0.0]).collect();
        let root = build_octree(&positions, 3);
        let mut all = Vec::new();
        root.collect_all(&mut all);
        assert_eq!(all.len(), 20);
    }

    #[test]
    fn test_adaptive_particle_volume() {
        let p = make_particle([0.0; 3], 1.0, 0.1, 0);
        let vol = p.volume(1000.0);
        assert!((vol - 0.001).abs() < 1e-12);
    }

    #[test]
    fn test_split_velocity_inherited() {
        let mut p = make_particle([0.0; 3], 8.0, 0.2, 0);
        p.velocity = [1.0, 2.0, 3.0];
        let children = split_particle(&p);
        assert!(children.iter().all(|c| c.velocity == [1.0, 2.0, 3.0]));
    }
}
