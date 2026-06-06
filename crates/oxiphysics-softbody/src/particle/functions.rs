//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{LifeParticle, Particle, ParticleCluster};

/// Compute a 3D Morton code (Z-order curve interleaving).
///
/// Given integer grid coordinates (x, y, z) each up to 21 bits, returns a 63-bit
/// Morton code by interleaving the bits.
pub fn morton_code_3d(x: u32, y: u32, z: u32) -> u64 {
    fn spread(v: u32) -> u64 {
        let mut v64 = v as u64;
        v64 = (v64 | (v64 << 32)) & 0x1f00000000ffff;
        v64 = (v64 | (v64 << 16)) & 0x1f0000ff0000ff;
        v64 = (v64 | (v64 << 8)) & 0x100f00f00f00f00f;
        v64 = (v64 | (v64 << 4)) & 0x10c30c30c30c30c3;
        v64 = (v64 | (v64 << 2)) & 0x1249249249249249;
        v64
    }
    spread(x) | (spread(y) << 1) | (spread(z) << 2)
}
/// Sort particles by their Morton code for cache-friendly spatial access.
///
/// `cell_size` is the size of each Morton grid cell. Returns the sort permutation
/// (indices into the original array, sorted by Morton code).
pub fn sort_particles_morton(particles: &[Particle], cell_size: f64) -> Vec<usize> {
    if cell_size < 1e-15 || particles.is_empty() {
        return (0..particles.len()).collect();
    }
    let min_x = particles
        .iter()
        .map(|p| p.position[0])
        .fold(f64::INFINITY, f64::min);
    let min_y = particles
        .iter()
        .map(|p| p.position[1])
        .fold(f64::INFINITY, f64::min);
    let min_z = particles
        .iter()
        .map(|p| p.position[2])
        .fold(f64::INFINITY, f64::min);
    let mut codes: Vec<(u64, usize)> = particles
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let gx = ((p.position[0] - min_x) / cell_size).floor() as u32;
            let gy = ((p.position[1] - min_y) / cell_size).floor() as u32;
            let gz = ((p.position[2] - min_z) / cell_size).floor() as u32;
            let code = morton_code_3d(gx, gy, gz);
            (code, i)
        })
        .collect();
    codes.sort_unstable_by_key(|&(code, _)| code);
    codes.into_iter().map(|(_, idx)| idx).collect()
}
/// Cluster particles into `k` groups using Lloyd's algorithm (k-means).
///
/// Returns a list of clusters with member indices.
pub fn cluster_particles(
    particles: &[Particle],
    k: usize,
    max_iter: usize,
) -> Vec<ParticleCluster> {
    if k == 0 || particles.is_empty() {
        return Vec::new();
    }
    let k = k.min(particles.len());
    let mut clusters: Vec<ParticleCluster> = (0..k)
        .map(|i| ParticleCluster::new(particles[i].position))
        .collect();
    for _ in 0..max_iter {
        for c in &mut clusters {
            c.members.clear();
            c.total_mass = 0.0;
        }
        for (idx, p) in particles.iter().enumerate() {
            let nearest = (0..k)
                .min_by(|&a, &b| {
                    let da = {
                        let dx = p.position[0] - clusters[a].centroid[0];
                        let dy = p.position[1] - clusters[a].centroid[1];
                        let dz = p.position[2] - clusters[a].centroid[2];
                        dx * dx + dy * dy + dz * dz
                    };
                    let db = {
                        let dx = p.position[0] - clusters[b].centroid[0];
                        let dy = p.position[1] - clusters[b].centroid[1];
                        let dz = p.position[2] - clusters[b].centroid[2];
                        dx * dx + dy * dy + dz * dz
                    };
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .expect("operation should succeed");
            clusters[nearest].members.push(idx);
        }
        for c in &mut clusters {
            c.update_centroid(particles);
        }
    }
    clusters
}
/// Interpolate a per-particle scalar property at an arbitrary position.
///
/// Uses the SPH kernel: `f(x) = Σ (m_j / ρ_j) * f_j * W(|x - x_j|, h)`
///
/// If density is unknown, pass `1.0` for each `densities[j]`.
pub fn interpolate_scalar(
    query: [f64; 3],
    particles: &[Particle],
    values: &[f64],
    densities: &[f64],
    h: f64,
) -> f64 {
    debug_assert_eq!(particles.len(), values.len());
    debug_assert_eq!(particles.len(), densities.len());
    let h2 = h * h;
    let mut num = 0.0_f64;
    let mut denom = 0.0_f64;
    for (i, p) in particles.iter().enumerate() {
        let dx = query[0] - p.position[0];
        let dy = query[1] - p.position[1];
        let dz = query[2] - p.position[2];
        let r2 = dx * dx + dy * dy + dz * dz;
        if r2 >= 4.0 * h2 {
            continue;
        }
        let r = r2.sqrt();
        let q = r / h;
        let sigma = 3.0 / (2.0 * std::f64::consts::PI * h * h * h);
        let w = sigma
            * if q < 1.0 {
                2.0 / 3.0 - q * q + 0.5 * q * q * q
            } else if q < 2.0 {
                (2.0 - q).powi(3) / 6.0
            } else {
                0.0
            };
        let rho = densities[i].max(1e-15);
        num += (p.mass / rho) * values[i] * w;
        denom += (p.mass / rho) * w;
    }
    if denom > 1e-15 { num / denom } else { 0.0 }
}
/// Interpolate a per-particle 3D vector property at an arbitrary position.
pub fn interpolate_vector(
    query: [f64; 3],
    particles: &[Particle],
    vectors: &[[f64; 3]],
    densities: &[f64],
    h: f64,
) -> [f64; 3] {
    debug_assert_eq!(particles.len(), vectors.len());
    let mut result = [0.0_f64; 3];
    let components: Vec<f64> = (0..3)
        .map(|d| {
            let scalars: Vec<f64> = vectors.iter().map(|v| v[d]).collect();
            interpolate_scalar(query, particles, &scalars, densities, h)
        })
        .collect();
    result.copy_from_slice(&components);
    result
}
/// Halton low-discrepancy sequence, element `n` in base `b`.
pub(super) fn halton(mut n: usize, b: usize) -> f64 {
    let mut f = 1.0_f64;
    let mut r = 0.0_f64;
    while n > 0 {
        f /= b as f64;
        r += f * (n % b) as f64;
        n /= b;
    }
    r
}
/// Simple 4D deterministic value noise in range \[−1, +1\].
pub(super) fn value_noise_4d(x: f64, y: f64, z: f64, w: f64, seed: u64) -> f64 {
    let ix = x.floor() as i64;
    let iy = y.floor() as i64;
    let iz = z.floor() as i64;
    let iw = w.floor() as i64;
    let fx = x - x.floor();
    let fy = y - y.floor();
    let fz = z - z.floor();
    let fw = w - w.floor();
    let ux = fade(fx);
    let uy = fade(fy);
    let uz = fade(fz);
    let uw = fade(fw);
    let v0000 = lattice_value(ix, iy, iz, iw, seed);
    let v1000 = lattice_value(ix + 1, iy, iz, iw, seed);
    let v0100 = lattice_value(ix, iy + 1, iz, iw, seed);
    let v1100 = lattice_value(ix + 1, iy + 1, iz, iw, seed);
    let v0010 = lattice_value(ix, iy, iz + 1, iw, seed);
    let v1010 = lattice_value(ix + 1, iy, iz + 1, iw, seed);
    let v0110 = lattice_value(ix, iy + 1, iz + 1, iw, seed);
    let v1110 = lattice_value(ix + 1, iy + 1, iz + 1, iw, seed);
    let v0001 = lattice_value(ix, iy, iz, iw + 1, seed);
    let v1001 = lattice_value(ix + 1, iy, iz, iw + 1, seed);
    let v0101 = lattice_value(ix, iy + 1, iz, iw + 1, seed);
    let v1101 = lattice_value(ix + 1, iy + 1, iz, iw + 1, seed);
    let v0011 = lattice_value(ix, iy, iz + 1, iw + 1, seed);
    let v1011 = lattice_value(ix + 1, iy, iz + 1, iw + 1, seed);
    let v0111 = lattice_value(ix, iy + 1, iz + 1, iw + 1, seed);
    let v1111 = lattice_value(ix + 1, iy + 1, iz + 1, iw + 1, seed);
    let a = lerp(v0000, v1000, ux);
    let b = lerp(v0100, v1100, ux);
    let c = lerp(v0010, v1010, ux);
    let d = lerp(v0110, v1110, ux);
    let e = lerp(v0001, v1001, ux);
    let f = lerp(v0101, v1101, ux);
    let g = lerp(v0011, v1011, ux);
    let h = lerp(v0111, v1111, ux);
    let ab = lerp(a, b, uy);
    let cd = lerp(c, d, uy);
    let ef = lerp(e, f, uy);
    let gh = lerp(g, h, uy);
    let abcd = lerp(ab, cd, uz);
    let efgh = lerp(ef, gh, uz);
    lerp(abcd, efgh, uw)
}
/// Perlin fade function: 6t⁵ − 15t⁴ + 10t³.
#[inline]
pub(super) fn fade(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}
/// Linear interpolation.
#[inline]
pub(super) fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + t * (b - a)
}
/// Deterministic lattice value in \[−1, +1\] for integer coordinates.
pub(super) fn lattice_value(x: i64, y: i64, z: i64, w: i64, seed: u64) -> f64 {
    let mut h = (x as u64)
        .wrapping_mul(2654435761)
        .wrapping_add((y as u64).wrapping_mul(805459861))
        .wrapping_add((z as u64).wrapping_mul(3674653429))
        .wrapping_add((w as u64).wrapping_mul(1482910573))
        .wrapping_add(seed.wrapping_mul(2246822519));
    h ^= h >> 33;
    h = h.wrapping_mul(0xff51afd7ed558ccd);
    h ^= h >> 33;
    h = h.wrapping_mul(0xc4ceb9fe1a85ec53);
    h ^= h >> 33;
    (h as i64) as f64 / i64::MAX as f64
}
/// RGBA colour as four `f32` components in \[0, 1\].
pub type Rgba = [f32; 4];
/// Sort particles by their signed depth along a camera axis for correct
/// alpha-blended rendering (back-to-front painter's algorithm).
///
/// `camera_dir` is the camera's forward direction (normalized).
/// Returns the indices of particles sorted from back (most negative depth)
/// to front (most positive depth) — i.e. furthest first for back-to-front.
pub fn sort_particles_by_depth(
    particles: &[Particle],
    camera_pos: [f64; 3],
    camera_dir: [f64; 3],
) -> Vec<usize> {
    let mut depths: Vec<(f64, usize)> = particles
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let dx = p.position[0] - camera_pos[0];
            let dy = p.position[1] - camera_pos[1];
            let dz = p.position[2] - camera_pos[2];
            let depth = dx * camera_dir[0] + dy * camera_dir[1] + dz * camera_dir[2];
            (depth, i)
        })
        .collect();
    depths.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    depths.into_iter().map(|(_, i)| i).collect()
}
/// Sort `LifeParticle` entries by depth for rendering.
///
/// Returns back-to-front order (most distant particle first).
pub fn sort_life_particles_by_depth(
    particles: &[LifeParticle],
    camera_pos: [f64; 3],
    camera_dir: [f64; 3],
) -> Vec<usize> {
    let mut depths: Vec<(f64, usize)> = particles
        .iter()
        .enumerate()
        .map(|(i, lp)| {
            let p = &lp.particle;
            let dx = p.position[0] - camera_pos[0];
            let dy = p.position[1] - camera_pos[1];
            let dz = p.position[2] - camera_pos[2];
            let depth = dx * camera_dir[0] + dy * camera_dir[1] + dz * camera_dir[2];
            (depth, i)
        })
        .collect();
    depths.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    depths.into_iter().map(|(_, i)| i).collect()
}
/// Compute distance from each particle to the camera and return sorted order.
///
/// Equivalent to depth sort by distance (not signed depth).
/// Returns back-to-front order.
pub fn sort_particles_by_distance(particles: &[Particle], camera_pos: [f64; 3]) -> Vec<usize> {
    let mut dists: Vec<(f64, usize)> = particles
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let dx = p.position[0] - camera_pos[0];
            let dy = p.position[1] - camera_pos[1];
            let dz = p.position[2] - camera_pos[2];
            let dist_sq = dx * dx + dy * dy + dz * dz;
            (dist_sq, i)
        })
        .collect();
    dists.sort_unstable_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    dists.into_iter().map(|(_, i)| i).collect()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::particle::BoxEmitter;
    use crate::particle::ColorOverLifetime;
    use crate::particle::ConeEmitter;
    use crate::particle::LifeParticlePool;
    use crate::particle::LineEmitter;
    use crate::particle::ParticleCollisionHandler;
    use crate::particle::ParticleFieldSampler;
    use crate::particle::ParticleFlags;
    use crate::particle::ParticleSet;
    use crate::particle::ParticleToGrid;
    use crate::particle::PointAttractor;
    use crate::particle::PointEmitter;
    use crate::particle::SizeOverLifetime;
    use crate::particle::SphSoftBody;
    use crate::particle::SphSoftBodyParticle;
    use crate::particle::SphereEmitter;
    use crate::particle::TurbulenceField;
    #[test]
    fn test_neighbors_within_finds_close() {
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([0.0, 0.0, 0.0], 1.0, 0.1));
        ps.add_particle(Particle::new([0.3, 0.0, 0.0], 1.0, 0.1));
        ps.add_particle(Particle::new([5.0, 0.0, 0.0], 1.0, 0.1));
        let nb = ps.neighbors_within([0.0, 0.0, 0.0], 0.5);
        assert!(nb.contains(&0), "Should find self");
        assert!(nb.contains(&1), "Should find close particle");
        assert!(!nb.contains(&2), "Should not find far particle");
    }
    #[test]
    fn test_neighbors_within_empty() {
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([10.0, 0.0, 0.0], 1.0, 0.1));
        let nb = ps.neighbors_within([0.0, 0.0, 0.0], 1.0);
        assert!(nb.is_empty(), "No particles should be within range");
    }
    #[test]
    fn test_kinetic_energy() {
        let mut ps = ParticleSet::new(1.0);
        let mut p = Particle::new([0.0; 3], 2.0, 0.1);
        p.velocity = [3.0, 4.0, 0.0];
        ps.add_particle(p);
        let ke = ps.kinetic_energy();
        assert!((ke - 25.0).abs() < 1e-10, "KE should be 25, got {ke}");
    }
    #[test]
    fn test_static_particle_zero_ke_momentum() {
        let mut ps = ParticleSet::new(1.0);
        let mut p = Particle::new([0.0; 3], 5.0, 0.1);
        p.velocity = [100.0, 100.0, 100.0];
        p.flags.set(ParticleFlags::STATIC);
        ps.add_particle(p);
        assert_eq!(ps.kinetic_energy(), 0.0);
        let mom = ps.total_momentum();
        assert_eq!(mom, [0.0; 3]);
    }
    #[test]
    fn test_impulse_momentum_conservation() {
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([0.0; 3], 1.0, 0.1));
        ps.add_particle(Particle::new([1.0, 0.0, 0.0], 2.0, 0.1));
        let p0_before = ps.total_momentum();
        let j = [5.0, -3.0, 1.0];
        ps.apply_impulse(0, j);
        ps.apply_impulse(1, [-j[0], -j[1], -j[2]]);
        let p_after = ps.total_momentum();
        for i in 0..3 {
            assert!(
                (p_after[i] - p0_before[i]).abs() < 1e-12,
                "Momentum[{i}] should be conserved: before={}, after={}",
                p0_before[i],
                p_after[i]
            );
        }
    }
    #[test]
    fn test_apply_external_force() {
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([0.0; 3], 2.0, 0.1));
        ps.apply_external_force(0, [10.0, 0.0, 0.0], 1.0);
        let vx = ps.particles[0].velocity[0];
        assert!(
            (vx - 5.0).abs() < 1e-12,
            "vx should be 5.0 after force, got {vx}"
        );
    }
    #[test]
    fn test_remove_particle() {
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([0.0; 3], 1.0, 0.1));
        ps.add_particle(Particle::new([1.0, 0.0, 0.0], 1.0, 0.1));
        ps.add_particle(Particle::new([2.0, 0.0, 0.0], 1.0, 0.1));
        ps.remove_particle(1);
        assert_eq!(
            ps.particles.len(),
            2,
            "Should have 2 particles after removal"
        );
    }
    #[test]
    fn test_step_integrates_positions() {
        let mut ps = ParticleSet::new(1.0);
        let mut p = Particle::new([0.0; 3], 1.0, 0.1);
        p.velocity = [1.0, 2.0, 3.0];
        ps.add_particle(p);
        ps.step(0.5);
        let pos = ps.particles[0].position;
        assert!((pos[0] - 0.5).abs() < 1e-12, "x mismatch: {}", pos[0]);
        assert!((pos[1] - 1.0).abs() < 1e-12, "y mismatch: {}", pos[1]);
        assert!((pos[2] - 1.5).abs() < 1e-12, "z mismatch: {}", pos[2]);
    }
    #[test]
    fn test_total_momentum_sum() {
        let mut ps = ParticleSet::new(1.0);
        let mut p1 = Particle::new([0.0; 3], 1.0, 0.1);
        p1.velocity = [1.0, 0.0, 0.0];
        let mut p2 = Particle::new([1.0, 0.0, 0.0], 2.0, 0.1);
        p2.velocity = [0.0, 3.0, 0.0];
        ps.add_particle(p1);
        ps.add_particle(p2);
        let mom = ps.total_momentum();
        assert!((mom[0] - 1.0).abs() < 1e-12, "px mismatch");
        assert!((mom[1] - 6.0).abs() < 1e-12, "py mismatch");
        assert!(mom[2].abs() < 1e-12, "pz should be 0");
    }
    #[test]
    fn test_point_emitter_creates_particles() {
        let mut emitter = PointEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 5.0);
        emitter.rate = 100.0;
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 0.1);
        assert_eq!(count, 10);
        assert_eq!(ps.particles.len(), 10);
        for p in &ps.particles {
            assert!((p.velocity[1] - 5.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_point_emitter_zero_dt() {
        let mut emitter = PointEmitter::new([0.0; 3], [1.0, 0.0, 0.0], 1.0);
        emitter.rate = 10.0;
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 0.0);
        assert_eq!(count, 0);
        assert!(ps.particles.is_empty());
    }
    #[test]
    fn test_line_emitter() {
        let emitter = LineEmitter::new([0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 5);
        assert_eq!(count, 5);
        assert_eq!(ps.particles.len(), 5);
        assert!((ps.particles[0].position[0] - 0.0).abs() < 1e-10);
        assert!((ps.particles[4].position[0] - 1.0).abs() < 1e-10);
    }
    #[test]
    fn test_line_emitter_single() {
        let emitter = LineEmitter::new([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3]);
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 1);
        assert_eq!(count, 1);
        assert!((ps.particles[0].position[0] - 0.5).abs() < 1e-10);
    }
    #[test]
    fn test_line_emitter_zero() {
        let emitter = LineEmitter::new([0.0; 3], [1.0, 0.0, 0.0], [0.0; 3]);
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 0);
        assert_eq!(count, 0);
    }
    #[test]
    fn test_attractor_pulls_particles() {
        let attractor = PointAttractor::new([0.0; 3], 100.0, 10.0);
        let mut ps = ParticleSet::new(1.0);
        let mut p = Particle::new([1.0, 0.0, 0.0], 1.0, 0.1);
        p.velocity = [0.0; 3];
        ps.add_particle(p);
        attractor.apply(&mut ps, 0.01);
        assert!(
            ps.particles[0].velocity[0] < 0.0,
            "Particle should be pulled toward attractor"
        );
    }
    #[test]
    fn test_attractor_out_of_range() {
        let attractor = PointAttractor::new([0.0; 3], 100.0, 0.5);
        let mut ps = ParticleSet::new(1.0);
        let p = Particle::new([10.0, 0.0, 0.0], 1.0, 0.1);
        ps.add_particle(p);
        attractor.apply(&mut ps, 0.01);
        assert_eq!(ps.particles[0].velocity, [0.0; 3]);
    }
    #[test]
    fn test_attractor_repel() {
        let attractor = PointAttractor::new([0.0; 3], -100.0, 10.0);
        let mut ps = ParticleSet::new(1.0);
        let mut p = Particle::new([1.0, 0.0, 0.0], 1.0, 0.1);
        p.velocity = [0.0; 3];
        ps.add_particle(p);
        attractor.apply(&mut ps, 0.01);
        assert!(
            ps.particles[0].velocity[0] > 0.0,
            "Negative strength should repel"
        );
    }
    #[test]
    fn test_collision_handler_overlapping() {
        let handler = ParticleCollisionHandler::new(1.0);
        let mut ps = ParticleSet::new(1.0);
        let mut p1 = Particle::new([0.0, 0.0, 0.0], 1.0, 0.5);
        p1.velocity = [1.0, 0.0, 0.0];
        let mut p2 = Particle::new([0.5, 0.0, 0.0], 1.0, 0.5);
        p2.velocity = [-1.0, 0.0, 0.0];
        ps.add_particle(p1);
        ps.add_particle(p2);
        let count = handler.resolve_collisions(&mut ps);
        assert_eq!(count, 1, "Should detect one collision");
    }
    #[test]
    fn test_collision_handler_separated() {
        let handler = ParticleCollisionHandler::new(1.0);
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([0.0, 0.0, 0.0], 1.0, 0.1));
        ps.add_particle(Particle::new([10.0, 0.0, 0.0], 1.0, 0.1));
        let count = handler.resolve_collisions(&mut ps);
        assert_eq!(count, 0, "No collision for separated particles");
    }
    #[test]
    fn test_sample_density_at_particle() {
        let particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        let density = ParticleFieldSampler::sample_density(&particles, [0.0; 3], 0.5);
        assert!(density > 0.0, "Density at particle position should be > 0");
    }
    #[test]
    fn test_sample_density_far_away() {
        let particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        let density = ParticleFieldSampler::sample_density(&particles, [100.0, 0.0, 0.0], 0.5);
        assert!(
            density.abs() < 1e-10,
            "Density far from particle should be ~0, got {density}"
        );
    }
    #[test]
    fn test_sample_velocity_at_particle() {
        let mut particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        particles[0].velocity = [1.0, 2.0, 3.0];
        let vel = ParticleFieldSampler::sample_velocity(&particles, [0.0; 3], 0.5);
        assert!((vel[0] - 1.0).abs() < 1e-6);
        assert!((vel[1] - 2.0).abs() < 1e-6);
        assert!((vel[2] - 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_particle_to_grid_total_cells() {
        let grid = ParticleToGrid::new([0.0; 3], 1.0, [3, 4, 5]);
        assert_eq!(grid.total_cells(), 60);
    }
    #[test]
    fn test_particle_to_grid_scatter() {
        let grid = ParticleToGrid::new([0.0; 3], 1.0, [5, 5, 5]);
        let particles = vec![Particle::new([2.5, 2.5, 2.5], 1.0, 0.1)];
        let field = grid.scatter_mass(&particles);
        assert_eq!(field.len(), 125);
        let total_mass: f64 = field.iter().sum();
        assert!(
            (total_mass - 1.0).abs() < 1e-10,
            "Total scattered mass should be 1.0, got {total_mass}"
        );
    }
    #[test]
    fn test_particle_to_grid_empty() {
        let grid = ParticleToGrid::new([0.0; 3], 1.0, [3, 3, 3]);
        let field = grid.scatter_mass(&[]);
        assert_eq!(field.len(), 27);
        assert!(field.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_particle_flags_operations() {
        let mut flags = ParticleFlags::none();
        assert!(!flags.has(ParticleFlags::STATIC));
        flags.set(ParticleFlags::STATIC);
        assert!(flags.has(ParticleFlags::STATIC));
        flags.clear(ParticleFlags::STATIC);
        assert!(!flags.has(ParticleFlags::STATIC));
    }
    #[test]
    fn test_particle_flags_multiple() {
        let mut flags = ParticleFlags::none();
        flags.set(ParticleFlags::STATIC);
        flags.set(ParticleFlags::ASLEEP);
        assert!(flags.has(ParticleFlags::STATIC));
        assert!(flags.has(ParticleFlags::ASLEEP));
        assert!(!flags.has(ParticleFlags::GHOST));
    }
    #[test]
    fn test_static_particle_creation() {
        let p = Particle::new_static([1.0, 2.0, 3.0], 0.5);
        assert!(p.is_static());
        assert_eq!(p.position, [1.0, 2.0, 3.0]);
        assert_eq!(p.radius, 0.5);
    }
    #[test]
    fn test_sph_soft_body_density_positive() {
        let mut body = SphSoftBody::new(1000.0, 10.0, 0.01);
        let h = 0.1;
        for ix in 0..3_usize {
            for iy in 0..3_usize {
                body.add_particle(SphSoftBodyParticle::new(
                    [ix as f64 * h * 0.5, iy as f64 * h * 0.5, 0.0],
                    0.01,
                    h,
                ));
            }
        }
        body.compute_densities();
        for p in &body.particles {
            assert!(p.density > 0.0, "Density should be positive");
        }
    }
    #[test]
    fn test_sph_soft_body_step_no_nan() {
        let mut body = SphSoftBody::new(1000.0, 100.0, 0.1);
        let h = 0.05;
        for ix in 0..3_usize {
            body.add_particle(SphSoftBodyParticle::new(
                [ix as f64 * h * 0.8, 0.0, 0.0],
                0.001,
                h,
            ));
        }
        body.step(1e-5);
        for p in &body.particles {
            assert!(p.position[0].is_finite());
            assert!(p.velocity[0].is_finite());
        }
    }
    #[test]
    fn test_sph_soft_body_pinned_does_not_move() {
        let mut body = SphSoftBody::new(1000.0, 100.0, 0.1);
        let h = 0.1;
        let mut p = SphSoftBodyParticle::new([0.0, 0.0, 0.0], 0.01, h);
        p.pinned = true;
        body.add_particle(p);
        body.add_particle(SphSoftBodyParticle::new([0.05, 0.0, 0.0], 0.01, h));
        body.step(0.01);
        assert_eq!(
            body.particles[0].position,
            [0.0, 0.0, 0.0],
            "Pinned particle should not move"
        );
    }
    #[test]
    fn test_morton_code_origin() {
        assert_eq!(morton_code_3d(0, 0, 0), 0);
    }
    #[test]
    fn test_morton_code_unit_x() {
        let m = morton_code_3d(1, 0, 0);
        assert_eq!(m & 1, 1, "X bit should be set in bit 0");
    }
    #[test]
    fn test_morton_code_unit_y() {
        let m = morton_code_3d(0, 1, 0);
        assert_eq!((m >> 1) & 1, 1, "Y bit should be set in bit 1");
    }
    #[test]
    fn test_morton_code_unit_z() {
        let m = morton_code_3d(0, 0, 1);
        assert_eq!((m >> 2) & 1, 1, "Z bit should be set in bit 2");
    }
    #[test]
    fn test_sort_particles_morton() {
        let particles = vec![
            Particle::new([3.0, 0.0, 0.0], 1.0, 0.1),
            Particle::new([0.0, 0.0, 0.0], 1.0, 0.1),
            Particle::new([1.0, 0.0, 0.0], 1.0, 0.1),
        ];
        let order = sort_particles_morton(&particles, 1.0);
        assert_eq!(
            order[0], 1,
            "Origin particle should be first in Morton order"
        );
    }
    #[test]
    fn test_cluster_particles_k1() {
        let particles: Vec<Particle> = (0..10)
            .map(|i| Particle::new([i as f64, 0.0, 0.0], 1.0, 0.1))
            .collect();
        let clusters = cluster_particles(&particles, 1, 5);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].members.len(), 10);
    }
    #[test]
    fn test_cluster_particles_k2_separated() {
        let mut particles: Vec<Particle> = (0..5)
            .map(|i| Particle::new([i as f64 * 0.1, 0.0, 0.0], 1.0, 0.05))
            .collect();
        for i in 0..5 {
            particles.push(Particle::new([100.0 + i as f64 * 0.1, 0.0, 0.0], 1.0, 0.05));
        }
        let clusters = cluster_particles(&particles, 2, 20);
        assert_eq!(clusters.len(), 2);
        assert!(!clusters[0].members.is_empty());
        assert!(!clusters[1].members.is_empty());
    }
    #[test]
    fn test_cluster_particles_empty() {
        let clusters = cluster_particles(&[], 3, 5);
        assert!(clusters.is_empty());
    }
    #[test]
    fn test_life_particle_aging() {
        let p = Particle::new([0.0; 3], 1.0, 0.1);
        let mut lp = LifeParticle::new(p, 1.0);
        assert!(lp.is_alive());
        lp.age(0.5);
        assert!(lp.is_alive());
        assert!((lp.age_fraction() - 0.5).abs() < 1e-10);
        lp.age(0.5);
        assert!(!lp.is_alive());
    }
    #[test]
    fn test_life_particle_pool_removes_dead() {
        let mut pool = LifeParticlePool::new();
        pool.spawn(Particle::new([0.0; 3], 1.0, 0.1), 0.05);
        pool.spawn(Particle::new([1.0, 0.0, 0.0], 1.0, 0.1), 10.0);
        pool.step(0.1);
        assert_eq!(pool.alive_count(), 1);
    }
    #[test]
    fn test_life_particle_pool_step_moves() {
        let mut pool = LifeParticlePool::new();
        let mut p = Particle::new([0.0; 3], 1.0, 0.1);
        p.velocity = [1.0, 0.0, 0.0];
        pool.spawn(p, 1.0);
        pool.step(0.1);
        assert!((pool.particles[0].particle.position[0] - 0.1).abs() < 1e-10);
    }
    #[test]
    fn test_interpolate_scalar_at_particle() {
        let particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        let values = vec![42.0_f64];
        let densities = vec![1000.0_f64];
        let v = interpolate_scalar([0.0; 3], &particles, &values, &densities, 0.2);
        assert!((v - 42.0).abs() < 0.1, "Should be close to 42, got {v}");
    }
    #[test]
    fn test_interpolate_scalar_zero_far_away() {
        let particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        let values = vec![100.0_f64];
        let densities = vec![1000.0_f64];
        let v = interpolate_scalar([1000.0, 0.0, 0.0], &particles, &values, &densities, 0.2);
        assert!(v.abs() < 1e-10, "Far away should give ~0, got {v}");
    }
    #[test]
    fn test_interpolate_vector() {
        let particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        let vectors = vec![[1.0_f64, 2.0, 3.0]];
        let densities = vec![1000.0_f64];
        let v = interpolate_vector([0.0; 3], &particles, &vectors, &densities, 0.2);
        assert!(v[0] > 0.0, "X component should be > 0");
    }
    #[test]
    fn test_cone_emitter_creates_particles() {
        let emitter = ConeEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 0.3, 5.0);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 8);
        assert_eq!(ps.particles.len(), 8, "Should create 8 particles");
    }
    #[test]
    fn test_cone_emitter_zero_spread_pencil_beam() {
        let emitter = ConeEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 5.0);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 4);
        for p in &ps.particles {
            let speed = (p.velocity[0] * p.velocity[0]
                + p.velocity[1] * p.velocity[1]
                + p.velocity[2] * p.velocity[2])
                .sqrt();
            assert!(
                (speed - 5.0).abs() < 0.01,
                "Speed should be 5.0, got {speed}"
            );
            assert!(p.velocity[0].abs() < 0.01, "vx should be 0 for pencil beam");
            assert!(p.velocity[2].abs() < 0.01, "vz should be 0 for pencil beam");
        }
    }
    #[test]
    fn test_cone_emitter_rate_based() {
        let mut emitter = ConeEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 0.5, 1.0);
        emitter.rate = 100.0;
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 0.1);
        assert_eq!(count, 10);
        assert_eq!(ps.particles.len(), 10);
    }
    #[test]
    fn test_cone_emitter_solid_angle_full_sphere() {
        let emitter = ConeEmitter::new([0.0; 3], [0.0, 1.0, 0.0], std::f64::consts::FRAC_PI_2, 1.0);
        let sa = emitter.solid_angle();
        assert!(
            (sa - 2.0 * std::f64::consts::PI).abs() < 0.01,
            "Half-space solid angle should be 2π, got {sa}"
        );
    }
    #[test]
    fn test_cone_emitter_solid_angle_pencil() {
        let emitter = ConeEmitter::new([0.0; 3], [0.0, 1.0, 0.0], 0.0, 1.0);
        let sa = emitter.solid_angle();
        assert!(sa < 0.01, "Pencil beam solid angle should be ~0, got {sa}");
    }
    #[test]
    fn test_sphere_emitter_surface_creates_particles() {
        let emitter = SphereEmitter::new([0.0; 3], 2.0, 1.0);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 6);
        assert_eq!(ps.particles.len(), 6);
    }
    #[test]
    fn test_sphere_emitter_particles_on_surface() {
        let emitter = SphereEmitter::new([1.0, 2.0, 3.0], 2.0, 1.0);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 20);
        for p in &ps.particles {
            let dx = p.position[0] - 1.0;
            let dy = p.position[1] - 2.0;
            let dz = p.position[2] - 3.0;
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                (r - 2.0).abs() < 0.01,
                "Particle should be on sphere surface (r=2), got r={r}"
            );
        }
    }
    #[test]
    fn test_sphere_emitter_volume_mode() {
        let mut emitter = SphereEmitter::new([0.0; 3], 2.0, 1.0);
        emitter.volume_emit = true;
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 20);
        for p in &ps.particles {
            let r = (p.position[0] * p.position[0]
                + p.position[1] * p.position[1]
                + p.position[2] * p.position[2])
                .sqrt();
            assert!(
                r <= 2.01,
                "Volume particle should be within radius 2, got r={r}"
            );
        }
    }
    #[test]
    fn test_sphere_emitter_surface_area() {
        let emitter = SphereEmitter::new([0.0; 3], 3.0, 1.0);
        let area = emitter.surface_area();
        let expected = 4.0 * std::f64::consts::PI * 9.0;
        assert!(
            (area - expected).abs() < 0.001,
            "Surface area should be 4πr², got {area}"
        );
    }
    #[test]
    fn test_sphere_emitter_rate_based() {
        let mut emitter = SphereEmitter::new([0.0; 3], 1.0, 1.0);
        emitter.rate = 50.0;
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 0.2);
        assert_eq!(count, 10);
    }
    #[test]
    fn test_box_emitter_creates_particles() {
        let emitter = BoxEmitter::new([0.0; 3], [1.0, 1.0, 1.0], [0.0; 3]);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 5);
        assert_eq!(ps.particles.len(), 5);
    }
    #[test]
    fn test_box_emitter_particles_inside_box() {
        let min = [-1.0, -1.0, -1.0];
        let max = [1.0, 1.0, 1.0];
        let emitter = BoxEmitter::new(min, max, [0.0; 3]);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 30);
        for p in &ps.particles {
            for d in 0..3 {
                assert!(
                    p.position[d] >= min[d] - 1e-10 && p.position[d] <= max[d] + 1e-10,
                    "Particle should be inside box, pos[{d}]={}",
                    p.position[d]
                );
            }
        }
    }
    #[test]
    fn test_box_emitter_velocity_set() {
        let vel = [1.0, 2.0, 3.0];
        let emitter = BoxEmitter::new([0.0; 3], [1.0; 3], vel);
        let mut ps = ParticleSet::new(1.0);
        emitter.emit_n(&mut ps, 4);
        for p in &ps.particles {
            assert_eq!(
                p.velocity, vel,
                "Particle velocity should match emitter velocity"
            );
        }
    }
    #[test]
    fn test_box_emitter_volume() {
        let emitter = BoxEmitter::new([0.0; 3], [2.0, 3.0, 4.0], [0.0; 3]);
        assert!(
            (emitter.volume() - 24.0).abs() < 1e-10,
            "Volume should be 24"
        );
    }
    #[test]
    fn test_box_emitter_rate_based() {
        let mut emitter = BoxEmitter::new([0.0; 3], [1.0; 3], [0.0; 3]);
        emitter.rate = 20.0;
        let mut ps = ParticleSet::new(1.0);
        let count = emitter.emit(&mut ps, 0.5);
        assert_eq!(count, 10);
    }
    #[test]
    fn test_turbulence_sample_nonzero() {
        let turb = TurbulenceField::new(1.0, 1.0);
        let v = turb.sample([0.5, 0.5, 0.5], 0.0);
        let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!(mag > 1e-6, "Turbulence sample should be nonzero, got {mag}");
    }
    #[test]
    fn test_turbulence_zero_strength() {
        let turb = TurbulenceField::new(0.0, 1.0);
        let v = turb.sample([1.0, 2.0, 3.0], 0.5);
        let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!(
            mag < 1e-10,
            "Zero strength turbulence should give zero velocity"
        );
    }
    #[test]
    fn test_turbulence_different_positions() {
        let turb = TurbulenceField::new(1.0, 1.0);
        let v1 = turb.sample([0.0, 0.0, 0.0], 0.0);
        let v2 = turb.sample([10.0, 0.0, 0.0], 0.0);
        let same = v1[0] == v2[0] && v1[1] == v2[1] && v1[2] == v2[2];
        assert!(
            !same,
            "Different positions should give different turbulence"
        );
    }
    #[test]
    fn test_turbulence_apply_changes_velocity() {
        let turb = TurbulenceField::new(5.0, 2.0);
        let mut ps = ParticleSet::new(1.0);
        ps.add_particle(Particle::new([1.0, 2.0, 3.0], 1.0, 0.1));
        let v_before = ps.particles[0].velocity;
        turb.apply(&mut ps, 0.0, 0.01);
        let v_after = ps.particles[0].velocity;
        let delta = ((v_after[0] - v_before[0]).powi(2)
            + (v_after[1] - v_before[1]).powi(2)
            + (v_after[2] - v_before[2]).powi(2))
        .sqrt();
        assert!(delta > 1e-10, "Turbulence should change velocity");
    }
    #[test]
    fn test_turbulence_advance_time() {
        let mut turb = TurbulenceField::new(1.0, 1.0);
        let t0 = turb.time_offset;
        turb.advance_time(0.5);
        assert!((turb.time_offset - (t0 + 0.5)).abs() < 1e-10);
    }
    #[test]
    fn test_color_over_lifetime_at_start() {
        let col = ColorOverLifetime::white_to_transparent();
        let c = col.evaluate(0.0);
        assert!(
            (c[3] - 1.0).abs() < 1e-6,
            "Alpha at start should be 1.0, got {}",
            c[3]
        );
    }
    #[test]
    fn test_color_over_lifetime_at_end() {
        let col = ColorOverLifetime::white_to_transparent();
        let c = col.evaluate(1.0);
        assert!(c[3] < 1e-6, "Alpha at end should be 0.0, got {}", c[3]);
    }
    #[test]
    fn test_color_over_lifetime_midpoint() {
        let col = ColorOverLifetime::white_to_transparent();
        let c = col.evaluate(0.5);
        assert!(
            (c[3] - 0.5).abs() < 1e-5,
            "Midpoint alpha should be ~0.5, got {}",
            c[3]
        );
    }
    #[test]
    fn test_color_over_lifetime_clamp() {
        let col = ColorOverLifetime::white_to_transparent();
        let c_neg = col.evaluate(-1.0);
        let c_over = col.evaluate(2.0);
        assert!(
            (c_neg[3] - 1.0).abs() < 1e-6,
            "Negative t should clamp to start"
        );
        assert!(c_over[3] < 1e-6, "t>1 should clamp to end");
    }
    #[test]
    fn test_color_over_lifetime_fire() {
        let col = ColorOverLifetime::fire();
        let c_start = col.evaluate(0.0);
        let c_end = col.evaluate(1.0);
        assert!((c_start[0] - 1.0).abs() < 1e-6, "Fire start R should be 1");
        assert!(c_end[3] < 0.01, "Fire end alpha should be near 0");
    }
    #[test]
    fn test_color_alpha_at() {
        let col = ColorOverLifetime::white_to_transparent();
        let a = col.alpha_at(0.25);
        assert!(
            (a - 0.75).abs() < 0.01,
            "alpha_at(0.25) should be ~0.75, got {a}"
        );
    }
    #[test]
    fn test_size_over_lifetime_constant() {
        let sol = SizeOverLifetime::constant();
        assert!((sol.evaluate(0.0) - 1.0).abs() < 1e-6);
        assert!((sol.evaluate(0.5) - 1.0).abs() < 1e-6);
        assert!((sol.evaluate(1.0) - 1.0).abs() < 1e-6);
    }
    #[test]
    fn test_size_over_lifetime_grow_shrink() {
        let sol = SizeOverLifetime::grow_shrink();
        let s0 = sol.evaluate(0.0);
        let s_mid = sol.evaluate(0.5);
        let s1 = sol.evaluate(1.0);
        assert!(s0 < 0.01, "Start size should be 0, got {s0}");
        assert!(
            (s_mid - 1.0).abs() < 0.01,
            "Mid size should be 1, got {s_mid}"
        );
        assert!(s1 < 0.01, "End size should be 0, got {s1}");
    }
    #[test]
    fn test_size_over_lifetime_clamp() {
        let sol = SizeOverLifetime::constant();
        assert!(
            (sol.evaluate(-0.5) - 1.0).abs() < 1e-6,
            "Negative t should clamp"
        );
        assert!((sol.evaluate(1.5) - 1.0).abs() < 1e-6, "t>1 should clamp");
    }
    #[test]
    fn test_size_over_lifetime_min_max() {
        let sol = SizeOverLifetime::grow_shrink();
        assert!(sol.min_size() < 0.01, "min_size should be 0");
        assert!((sol.max_size() - 1.0).abs() < 0.01, "max_size should be 1");
    }
    #[test]
    fn test_size_over_lifetime_custom() {
        let sol = SizeOverLifetime::new(vec![(0.0, 0.5), (1.0, 2.0)]);
        assert!((sol.evaluate(0.0) - 0.5).abs() < 1e-6);
        assert!((sol.evaluate(1.0) - 2.0).abs() < 1e-6);
        assert!((sol.evaluate(0.5) - 1.25).abs() < 0.01);
    }
    #[test]
    fn test_sort_particles_by_depth_order() {
        let particles = vec![
            Particle::new([0.0, 0.0, 1.0], 1.0, 0.1),
            Particle::new([0.0, 0.0, 5.0], 1.0, 0.1),
        ];
        let camera_pos = [0.0, 0.0, 0.0];
        let camera_dir = [0.0, 0.0, 1.0];
        let order = sort_particles_by_depth(&particles, camera_pos, camera_dir);
        assert_eq!(order[0], 1, "Far particle should be first (depth=5)");
        assert_eq!(order[1], 0, "Near particle should be second (depth=1)");
    }
    #[test]
    fn test_sort_particles_by_depth_single() {
        let particles = vec![Particle::new([0.0; 3], 1.0, 0.1)];
        let order = sort_particles_by_depth(&particles, [0.0; 3], [0.0, 0.0, 1.0]);
        assert_eq!(order, vec![0]);
    }
    #[test]
    fn test_sort_particles_by_depth_empty() {
        let order = sort_particles_by_depth(&[], [0.0; 3], [0.0, 0.0, 1.0]);
        assert!(order.is_empty());
    }
    #[test]
    fn test_sort_particles_by_distance_order() {
        let particles = vec![
            Particle::new([0.0, 0.0, 1.0], 1.0, 0.1),
            Particle::new([0.0, 0.0, 10.0], 1.0, 0.1),
            Particle::new([0.0, 0.0, 3.0], 1.0, 0.1),
        ];
        let order = sort_particles_by_distance(&particles, [0.0; 3]);
        assert_eq!(order[0], 1, "Farthest should be first");
        assert_eq!(order[2], 0, "Nearest should be last");
    }
    #[test]
    fn test_sort_life_particles_by_depth() {
        let mut p_near = Particle::new([0.0, 0.0, 2.0], 1.0, 0.1);
        p_near.velocity = [0.0; 3];
        let mut p_far = Particle::new([0.0, 0.0, 8.0], 1.0, 0.1);
        p_far.velocity = [0.0; 3];
        let life_particles = vec![
            LifeParticle::new(p_near, 1.0),
            LifeParticle::new(p_far, 1.0),
        ];
        let order = sort_life_particles_by_depth(&life_particles, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0]);
        assert_eq!(order[0], 1, "Far particle (depth 8) should come first");
        assert_eq!(order[1], 0, "Near particle (depth 2) should come second");
    }
    #[test]
    fn test_halton_range() {
        for n in 1..=20 {
            let h = halton(n, 2);
            assert!(
                (0.0..1.0).contains(&h),
                "Halton value should be in [0,1), got {h}"
            );
        }
    }
    #[test]
    fn test_halton_base2_first() {
        assert!((halton(1, 2) - 0.5).abs() < 1e-10);
        assert!((halton(2, 2) - 0.25).abs() < 1e-10);
        assert!((halton(3, 2) - 0.75).abs() < 1e-10);
    }
}
