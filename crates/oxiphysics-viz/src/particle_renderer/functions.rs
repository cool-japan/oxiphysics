//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    BillboardParticle, ColorMap, ParticleBuffer, ParticleInstance, ScalarColorizer, VelocityGlyph,
};

/// Linearly interpolate between two RGBA colors.
pub(super) fn lerp_color(a: [f32; 4], b: [f32; 4], t: f32) -> [f32; 4] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
        a[3] + (b[3] - a[3]) * t,
    ]
}
/// Evaluate a piecewise-linear colormap defined by `stops` at parameter `t` ∈ \[0, 1\].
pub(super) fn piecewise_lerp(stops: &[[f32; 4]], t: f32) -> [f32; 4] {
    let t = t.clamp(0.0, 1.0);
    let n = stops.len();
    if n == 0 {
        return [0.0, 0.0, 0.0, 1.0];
    }
    if n == 1 {
        return stops[0];
    }
    let seg = (n - 1) as f32 * t;
    let idx = (seg as usize).min(n - 2);
    let local_t = seg - idx as f32;
    lerp_color(stops[idx], stops[idx + 1], local_t)
}
/// Compute a level-of-detail radius for a particle given the camera distance.
///
/// `radius = base_radius * lod_factor / camera_distance`
#[allow(dead_code)]
pub fn particle_lod_radius(camera_distance: f32, base_radius: f32, lod_factor: f32) -> f32 {
    base_radius * lod_factor / camera_distance.max(f32::EPSILON)
}
/// Generate billboard vertex data for a set of particles.
///
/// Each particle becomes a camera-facing quad with 4 vertices.
/// `camera_right` and `camera_up` are unit vectors defining the billboard
/// orientation.
///
/// Returns interleaved data: `[px, py, pz, u, v, cr, cg, cb, ca]` per vertex.
#[allow(dead_code)]
pub fn generate_billboard_vertices(
    particles: &[BillboardParticle],
    camera_right: [f32; 3],
    camera_up: [f32; 3],
) -> Vec<f32> {
    let mut out = Vec::with_capacity(particles.len() * 4 * 9);
    let offsets: [(f32, f32, f32, f32); 4] = [
        (-1.0, -1.0, 0.0, 0.0),
        (1.0, -1.0, 1.0, 0.0),
        (1.0, 1.0, 1.0, 1.0),
        (-1.0, 1.0, 0.0, 1.0),
    ];
    for p in particles {
        for &(ox, oy, u, v) in &offsets {
            let vx = p.position[0]
                + camera_right[0] * ox * p.half_size
                + camera_up[0] * oy * p.half_size;
            let vy = p.position[1]
                + camera_right[1] * ox * p.half_size
                + camera_up[1] * oy * p.half_size;
            let vz = p.position[2]
                + camera_right[2] * ox * p.half_size
                + camera_up[2] * oy * p.half_size;
            out.extend_from_slice(&[vx, vy, vz, u, v]);
            out.extend_from_slice(&p.color);
        }
    }
    out
}
/// Sort particles by distance from the camera (back-to-front for alpha blending).
///
/// Modifies `particles` in place, sorting by decreasing distance from `camera_pos`.
#[allow(dead_code)]
pub fn sort_particles_back_to_front(particles: &mut [ParticleInstance], camera_pos: [f32; 3]) {
    particles.sort_by(|a, b| {
        let da = dist_sq(a.position, camera_pos);
        let db = dist_sq(b.position, camera_pos);
        db.partial_cmp(&da).unwrap_or(std::cmp::Ordering::Equal)
    });
}
/// Sort particles front-to-back (useful for early-Z rejection in opaque passes).
#[allow(dead_code)]
pub fn sort_particles_front_to_back(particles: &mut [ParticleInstance], camera_pos: [f32; 3]) {
    particles.sort_by(|a, b| {
        let da = dist_sq(a.position, camera_pos);
        let db = dist_sq(b.position, camera_pos);
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    });
}
pub(super) fn dist_sq(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    dx * dx + dy * dy + dz * dz
}
/// Compute an attenuated particle size based on camera distance.
///
/// Uses the model: `apparent_size = base_size * reference_dist / distance`.
/// Clamped to `[min_size, max_size]`.
#[allow(dead_code)]
pub fn attenuated_size(
    base_size: f32,
    distance: f32,
    reference_dist: f32,
    min_size: f32,
    max_size: f32,
) -> f32 {
    let size = base_size * reference_dist / distance.max(f32::EPSILON);
    size.clamp(min_size, max_size)
}
/// Apply size attenuation to all particles in a buffer based on camera position.
#[allow(dead_code)]
pub fn attenuate_particle_sizes(
    buffer: &mut ParticleBuffer,
    camera_pos: [f32; 3],
    reference_dist: f32,
    min_size: f32,
    max_size: f32,
) {
    for p in &mut buffer.instances {
        let d = dist_sq(p.position, camera_pos).sqrt();
        p.radius = attenuated_size(p.radius, d, reference_dist, min_size, max_size);
    }
}
/// Color particles based on their velocity magnitude using a colorizer.
#[allow(dead_code)]
pub fn color_by_velocity(buffer: &mut ParticleBuffer, colorizer: &ScalarColorizer) {
    for p in &mut buffer.instances {
        let speed = (p.velocity[0] * p.velocity[0]
            + p.velocity[1] * p.velocity[1]
            + p.velocity[2] * p.velocity[2])
            .sqrt();
        p.color = colorizer.colorize(speed);
    }
}
/// Set velocities for particles in a buffer.
#[allow(dead_code)]
pub fn set_particle_velocities(buffer: &mut ParticleBuffer, velocities: &[[f32; 3]]) {
    for (p, v) in buffer.instances.iter_mut().zip(velocities.iter()) {
        p.velocity = *v;
    }
}
/// Color particles by an arbitrary scalar property using a colormap.
///
/// `props[i]` is the scalar property of particle `i`.
/// Values are clamped to `[vmin, vmax]` before colormap lookup.
#[allow(dead_code)]
pub fn color_by_property(
    buffer: &mut ParticleBuffer,
    props: &[f32],
    vmin: f32,
    vmax: f32,
    color_map: ColorMap,
) {
    let range = vmax - vmin;
    for (p, &prop) in buffer.instances.iter_mut().zip(props.iter()) {
        let t = if range.abs() < f32::EPSILON {
            0.5
        } else {
            ((prop - vmin) / range).clamp(0.0, 1.0)
        };
        p.color = color_map.map(t);
    }
}
/// Build velocity glyph arrows for a set of particles.
///
/// `scale` is the length multiplier applied to the velocity vector.
#[allow(dead_code)]
pub fn build_velocity_glyphs(
    positions: &[[f32; 3]],
    velocities: &[[f32; 3]],
    scale: f32,
) -> Vec<VelocityGlyph> {
    let max_speed = velocities
        .iter()
        .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
        .fold(0.0f32, f32::max);
    positions
        .iter()
        .zip(velocities.iter())
        .map(|(&base, &vel)| {
            let speed = (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]).sqrt();
            let tip = [
                base[0] + vel[0] * scale,
                base[1] + vel[1] * scale,
                base[2] + vel[2] * scale,
            ];
            let t = if max_speed > f32::EPSILON {
                speed / max_speed
            } else {
                0.0
            };
            let col = ColorMap::Hot.map(t);
            VelocityGlyph {
                base,
                tip,
                color: col,
                magnitude: speed,
            }
        })
        .collect()
}
/// Produce line vertex data from velocity glyphs:
/// `[bx, by, bz, tx, ty, tz, cr, cg, cb, ca]` per glyph (10 floats).
#[allow(dead_code)]
pub fn velocity_glyphs_to_line_data(glyphs: &[VelocityGlyph]) -> Vec<f32> {
    let mut out = Vec::with_capacity(glyphs.len() * 10);
    for g in glyphs {
        out.extend_from_slice(&g.base);
        out.extend_from_slice(&g.tip);
        out.extend_from_slice(&g.color);
    }
    out
}
/// Frustum-cull a set of spheres against six half-planes.
///
/// `frustum_planes` is an array of 6 planes, each represented as `[a, b, c, d]`
/// where `ax + by + cz + d >= 0` defines the inside half-space.
/// Returns a `Vec`bool` where `true` means the sphere is *visible* (not culled).
#[allow(dead_code)]
pub fn frustum_cull_particles(
    positions: &[[f64; 3]],
    radii: &[f64],
    frustum_planes: &[[f64; 4]; 6],
) -> Vec<bool> {
    positions
        .iter()
        .zip(radii.iter())
        .map(|(pos, &r)| {
            for plane in frustum_planes.iter() {
                let dist = plane[0] * pos[0] + plane[1] * pos[1] + plane[2] * pos[2] + plane[3];
                if dist < -r {
                    return false;
                }
            }
            true
        })
        .collect()
}
/// Return the LOD level for a particle at `distance` from the camera.
///
/// `lod_distances` is a sorted ascending slice of camera distances.
/// Level 0 = highest detail (closest).
/// If `distance` is beyond all entries the last level index is returned.
#[allow(dead_code)]
pub fn particle_lod(distance: f64, lod_distances: &[f64]) -> usize {
    for (i, &d) in lod_distances.iter().enumerate() {
        if distance < d {
            return i;
        }
    }
    lod_distances.len().saturating_sub(1)
}
/// Generate a circle of billboard vertices that approximate a sphere impostor.
///
/// Returns `n_segments` vertices arranged in a regular polygon around the
/// origin (in the XY plane).  Each vertex is `\[x, y, 0.0\]` scaled by `radius`.
#[allow(dead_code)]
pub fn generate_sphere_impostor(radius: f64, n_segments: usize) -> Vec<[f32; 3]> {
    let n = n_segments.max(3);
    let mut verts = Vec::with_capacity(n);
    for i in 0..n {
        let angle = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
        verts.push([
            (radius * angle.cos()) as f32,
            (radius * angle.sin()) as f32,
            0.0_f32,
        ]);
    }
    verts
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_viz::Metaball;
    use crate::fluid_viz::MetaballField;
    use crate::gpu_viz::LodSystem;
    use crate::particle_renderer::ArrowBuffer;
    use crate::particle_renderer::ParticleSortBuffer;
    use crate::particle_renderer::ParticleTrail;
    use crate::particle_renderer::PointSpriteConfig;
    use crate::particle_renderer::PointSpriteRenderer;
    use crate::particle_renderer::SphDensityField;
    use crate::particle_renderer::SplatRenderer;
    #[test]
    fn test_particle_buffer_new_empty() {
        let buf = ParticleBuffer::new(10);
        assert_eq!(buf.count(), 0);
        assert_eq!(buf.max_particles, 10);
    }
    #[test]
    fn test_particle_buffer_add_returns_true() {
        let mut buf = ParticleBuffer::new(5);
        let ok = buf.add([0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0, 1.0]);
        assert!(ok);
        assert_eq!(buf.count(), 1);
    }
    #[test]
    fn test_particle_buffer_full_returns_false() {
        let mut buf = ParticleBuffer::new(2);
        assert!(buf.add([0.0, 0.0, 0.0], 1.0, [1.0, 0.0, 0.0, 1.0]));
        assert!(buf.add([1.0, 0.0, 0.0], 1.0, [0.0, 1.0, 0.0, 1.0]));
        let full = buf.add([2.0, 0.0, 0.0], 1.0, [0.0, 0.0, 1.0, 1.0]);
        assert!(!full);
        assert_eq!(buf.count(), 2);
    }
    #[test]
    fn test_particle_buffer_clear() {
        let mut buf = ParticleBuffer::new(5);
        buf.add([0.0, 0.0, 0.0], 1.0, [1.0, 1.0, 1.0, 1.0]);
        buf.add([1.0, 0.0, 0.0], 1.0, [1.0, 1.0, 1.0, 1.0]);
        buf.clear();
        assert_eq!(buf.count(), 0);
    }
    #[test]
    fn test_particle_buffer_to_vertex_data_length() {
        let mut buf = ParticleBuffer::new(3);
        buf.add([0.0, 1.0, 2.0], 0.5, [1.0, 0.0, 0.0, 1.0]);
        buf.add([3.0, 4.0, 5.0], 0.5, [0.0, 1.0, 0.0, 1.0]);
        let data = buf.to_vertex_data();
        assert_eq!(data.len(), 2 * 8);
    }
    #[test]
    fn test_particle_buffer_to_vertex_data_values() {
        let mut buf = ParticleBuffer::new(1);
        buf.add([1.0, 2.0, 3.0], 0.25, [0.5, 0.6, 0.7, 1.0]);
        let data = buf.to_vertex_data();
        assert!((data[0] - 1.0).abs() < f32::EPSILON);
        assert!((data[1] - 2.0).abs() < f32::EPSILON);
        assert!((data[2] - 3.0).abs() < f32::EPSILON);
        assert!((data[3] - 0.25).abs() < f32::EPSILON);
        assert!((data[4] - 0.5).abs() < f32::EPSILON);
    }
    #[test]
    fn test_particle_buffer_from_positions() {
        let positions = [[0.0f32, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]];
        let buf = ParticleBuffer::from_positions(&positions, 0.1, [1.0, 1.0, 1.0, 1.0]);
        assert_eq!(buf.count(), 3);
        assert_eq!(buf.max_particles, 3);
    }
    #[test]
    fn test_particle_buffer_from_positions_empty() {
        let buf = ParticleBuffer::from_positions(&[], 0.1, [1.0, 1.0, 1.0, 1.0]);
        assert_eq!(buf.count(), 0);
    }
    #[test]
    fn test_colormap_grayscale_zero_is_black() {
        let c = ColorMap::Grayscale.map(0.0);
        assert!((c[0]).abs() < f32::EPSILON);
        assert!((c[1]).abs() < f32::EPSILON);
        assert!((c[2]).abs() < f32::EPSILON);
        assert!((c[3] - 1.0).abs() < f32::EPSILON);
    }
    #[test]
    fn test_colormap_grayscale_one_is_white() {
        let c = ColorMap::Grayscale.map(1.0);
        assert!((c[0] - 1.0).abs() < f32::EPSILON);
        assert!((c[1] - 1.0).abs() < f32::EPSILON);
        assert!((c[2] - 1.0).abs() < f32::EPSILON);
    }
    #[test]
    fn test_colormap_cool_midpoint() {
        let c = ColorMap::Cool.map(0.5);
        assert!((c[0] - 0.5).abs() < 1e-5);
        assert!((c[1] - 0.5).abs() < 1e-5);
        assert!((c[2] - 1.0).abs() < 1e-5);
        assert!((c[3] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_colormap_hot_zero_is_black() {
        let c = ColorMap::Hot.map(0.0);
        assert!(c[0].abs() < 1e-5);
        assert!(c[1].abs() < 1e-5);
        assert!(c[2].abs() < 1e-5);
    }
    #[test]
    fn test_colormap_hot_one_is_white() {
        let c = ColorMap::Hot.map(1.0);
        assert!((c[0] - 1.0).abs() < 1e-5);
        assert!((c[1] - 1.0).abs() < 1e-5);
        assert!((c[2] - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_colormap_viridis_endpoints_distinct() {
        let c0 = ColorMap::Viridis.map(0.0);
        let c1 = ColorMap::Viridis.map(1.0);
        let diff = (c0[0] - c1[0]).abs() + (c0[1] - c1[1]).abs() + (c0[2] - c1[2]).abs();
        assert!(diff > 0.1, "Viridis endpoints should differ");
    }
    #[test]
    fn test_colormap_plasma_midpoint_in_range() {
        let c = ColorMap::Plasma.map(0.5);
        for component in &c[..3] {
            assert!(*component >= 0.0 && *component <= 1.0);
        }
    }
    #[test]
    fn test_colormap_clamps_below_zero() {
        let c_neg = ColorMap::Cool.map(-0.5);
        let c_zero = ColorMap::Cool.map(0.0);
        assert_eq!(c_neg, c_zero);
    }
    #[test]
    fn test_colormap_clamps_above_one() {
        let c_over = ColorMap::Grayscale.map(2.0);
        let c_one = ColorMap::Grayscale.map(1.0);
        assert_eq!(c_over, c_one);
    }
    #[test]
    fn test_scalar_colorizer_normalizes() {
        let col = ScalarColorizer::new(0.0, 10.0, ColorMap::Grayscale);
        let c5 = col.colorize(5.0);
        assert!((c5[0] - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_scalar_colorizer_zero_range() {
        let col = ScalarColorizer::new(5.0, 5.0, ColorMap::Grayscale);
        let c = col.colorize(5.0);
        assert!((c[0] - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_scalar_colorizer_array() {
        let col = ScalarColorizer::new(0.0, 1.0, ColorMap::Grayscale);
        let values = [0.0f32, 0.5, 1.0];
        let colors = col.colorize_array(&values);
        assert_eq!(colors.len(), 3);
        assert!((colors[1][0] - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_arrow_buffer_add_and_count() {
        let col = ScalarColorizer::new(0.0, 1.0, ColorMap::Grayscale);
        let mut buf = ArrowBuffer::new();
        buf.add([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, &col);
        assert_eq!(buf.arrows.len(), 1);
    }
    #[test]
    fn test_arrow_buffer_to_line_data_length() {
        let col = ScalarColorizer::new(0.0, 5.0, ColorMap::Grayscale);
        let mut buf = ArrowBuffer::new();
        buf.add([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, &col);
        buf.add([1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 1.0, &col);
        let data = buf.to_line_data();
        assert_eq!(data.len(), 2 * 6);
    }
    #[test]
    fn test_arrow_buffer_max_magnitude_empty() {
        let buf = ArrowBuffer::new();
        assert!((buf.max_magnitude()).abs() < f32::EPSILON);
    }
    #[test]
    fn test_arrow_buffer_max_magnitude() {
        let col = ScalarColorizer::new(0.0, 10.0, ColorMap::Grayscale);
        let mut buf = ArrowBuffer::new();
        buf.add([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 1.0, &col);
        buf.add([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, &col);
        assert!((buf.max_magnitude() - 5.0).abs() < 1e-4);
    }
    #[test]
    fn test_arrow_zero_velocity() {
        let col = ScalarColorizer::new(0.0, 1.0, ColorMap::Grayscale);
        let mut buf = ArrowBuffer::new();
        buf.add([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, &col);
        let data = buf.to_line_data();
        assert_eq!(data.len(), 6);
        for v in &data {
            assert!(v.abs() < f32::EPSILON);
        }
    }
    #[test]
    fn test_particle_lod_radius_basic() {
        let r = particle_lod_radius(10.0, 1.0, 5.0);
        assert!((r - 0.5).abs() < 1e-5);
    }
    #[test]
    fn test_particle_lod_radius_zero_distance() {
        let r = particle_lod_radius(0.0, 1.0, 1.0);
        assert!(r.is_finite());
    }
    #[test]
    fn test_billboard_vertex_count() {
        let particles = vec![BillboardParticle {
            position: [0.0, 0.0, 0.0],
            half_size: 1.0,
            color: [1.0, 1.0, 1.0, 1.0],
            texture_id: 0,
        }];
        let data = generate_billboard_vertices(&particles, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(data.len(), 36);
    }
    #[test]
    fn test_billboard_two_particles() {
        let particles = vec![
            BillboardParticle {
                position: [0.0, 0.0, 0.0],
                half_size: 0.5,
                color: [1.0, 0.0, 0.0, 1.0],
                texture_id: 0,
            },
            BillboardParticle {
                position: [3.0, 0.0, 0.0],
                half_size: 0.5,
                color: [0.0, 1.0, 0.0, 1.0],
                texture_id: 0,
            },
        ];
        let data = generate_billboard_vertices(&particles, [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
        assert_eq!(data.len(), 72);
    }
    #[test]
    fn test_sort_back_to_front() {
        let mut particles = vec![
            ParticleInstance {
                position: [1.0, 0.0, 0.0],
                radius: 0.1,
                color: [1.0; 4],
                velocity: [0.0; 3],
            },
            ParticleInstance {
                position: [5.0, 0.0, 0.0],
                radius: 0.1,
                color: [1.0; 4],
                velocity: [0.0; 3],
            },
            ParticleInstance {
                position: [3.0, 0.0, 0.0],
                radius: 0.1,
                color: [1.0; 4],
                velocity: [0.0; 3],
            },
        ];
        sort_particles_back_to_front(&mut particles, [0.0, 0.0, 0.0]);
        assert!((particles[0].position[0] - 5.0).abs() < f32::EPSILON);
        assert!((particles[1].position[0] - 3.0).abs() < f32::EPSILON);
        assert!((particles[2].position[0] - 1.0).abs() < f32::EPSILON);
    }
    #[test]
    fn test_sort_front_to_back() {
        let mut particles = vec![
            ParticleInstance {
                position: [5.0, 0.0, 0.0],
                radius: 0.1,
                color: [1.0; 4],
                velocity: [0.0; 3],
            },
            ParticleInstance {
                position: [1.0, 0.0, 0.0],
                radius: 0.1,
                color: [1.0; 4],
                velocity: [0.0; 3],
            },
        ];
        sort_particles_front_to_back(&mut particles, [0.0, 0.0, 0.0]);
        assert!((particles[0].position[0] - 1.0).abs() < f32::EPSILON);
        assert!((particles[1].position[0] - 5.0).abs() < f32::EPSILON);
    }
    #[test]
    fn test_attenuated_size_at_reference() {
        let s = attenuated_size(1.0, 10.0, 10.0, 0.1, 5.0);
        assert!((s - 1.0).abs() < 1e-5);
    }
    #[test]
    fn test_attenuated_size_far() {
        let s = attenuated_size(1.0, 100.0, 10.0, 0.1, 5.0);
        assert!((s - 0.1).abs() < 1e-5);
    }
    #[test]
    fn test_attenuated_size_close() {
        let s = attenuated_size(1.0, 1.0, 10.0, 0.1, 5.0);
        assert!((s - 5.0).abs() < 1e-5);
    }
    #[test]
    fn test_attenuate_buffer() {
        let mut buf =
            ParticleBuffer::from_positions(&[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]], 1.0, [1.0; 4]);
        attenuate_particle_sizes(&mut buf, [0.0, 0.0, 0.0], 5.0, 0.1, 10.0);
        assert!(buf.instances[0].radius > buf.instances[1].radius);
    }
    #[test]
    fn test_color_by_velocity() {
        let mut buf = ParticleBuffer::new(2);
        buf.instances.push(ParticleInstance {
            position: [0.0; 3],
            radius: 0.1,
            color: [1.0; 4],
            velocity: [3.0, 4.0, 0.0],
        });
        buf.instances.push(ParticleInstance {
            position: [1.0, 0.0, 0.0],
            radius: 0.1,
            color: [1.0; 4],
            velocity: [0.0; 3],
        });
        let col = ScalarColorizer::new(0.0, 10.0, ColorMap::Grayscale);
        color_by_velocity(&mut buf, &col);
        assert!((buf.instances[0].color[0] - 0.5).abs() < 1e-5);
        assert!(buf.instances[1].color[0].abs() < 1e-5);
    }
    #[test]
    fn test_trail_add_point() {
        let mut trail = ParticleTrail::new(10, 5.0);
        trail.add_point([0.0, 0.0, 0.0], 0.1, [1.0; 4]);
        assert_eq!(trail.point_count(), 1);
    }
    #[test]
    fn test_trail_max_points() {
        let mut trail = ParticleTrail::new(3, 5.0);
        for i in 0..5 {
            trail.add_point([i as f32, 0.0, 0.0], 0.1, [1.0; 4]);
        }
        assert_eq!(trail.point_count(), 3);
    }
    #[test]
    fn test_trail_aging_removes_old() {
        let mut trail = ParticleTrail::new(10, 1.0);
        trail.add_point([0.0, 0.0, 0.0], 0.1, [1.0; 4]);
        trail.update(2.0);
        assert_eq!(trail.point_count(), 0);
    }
    #[test]
    fn test_trail_line_strip_data() {
        let mut trail = ParticleTrail::new(10, 5.0);
        trail.add_point([1.0, 2.0, 3.0], 0.1, [1.0, 0.0, 0.0, 1.0]);
        let data = trail.to_line_strip_data();
        assert_eq!(data.len(), 7);
        assert!((data[0] - 1.0).abs() < 1e-5);
        assert!((data[1] - 2.0).abs() < 1e-5);
    }
    #[test]
    fn test_trail_alpha_fade() {
        let mut trail = ParticleTrail::new(10, 2.0);
        trail.add_point([0.0; 3], 0.1, [1.0, 1.0, 1.0, 1.0]);
        trail.update(1.0);
        let data = trail.to_line_strip_data();
        assert!(
            (data[6] - 0.5).abs() < 1e-5,
            "expected 0.5, got {}",
            data[6]
        );
    }
    #[test]
    fn test_point_sprite_renderer_basic() {
        let mut renderer = PointSpriteRenderer::new(PointSpriteConfig {
            size_world: 0.1,
            use_velocity_stretch: false,
            stretch_factor: 1.0,
        });
        renderer.add([1.0, 2.0, 3.0], [0.0; 3], [1.0, 0.0, 0.0, 1.0]);
        renderer.add([4.0, 5.0, 6.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0, 1.0]);
        assert_eq!(renderer.count(), 2);
    }
    #[test]
    fn test_point_sprite_renderer_clear() {
        let mut renderer = PointSpriteRenderer::new(PointSpriteConfig {
            size_world: 0.1,
            use_velocity_stretch: false,
            stretch_factor: 1.0,
        });
        renderer.add([0.0; 3], [0.0; 3], [1.0; 4]);
        renderer.clear();
        assert_eq!(renderer.count(), 0);
    }
    #[test]
    fn test_point_sprite_to_vertex_data() {
        let mut renderer = PointSpriteRenderer::new(PointSpriteConfig {
            size_world: 0.5,
            use_velocity_stretch: true,
            stretch_factor: 2.0,
        });
        renderer.add([1.0, 0.0, 0.0], [3.0, 0.0, 0.0], [1.0, 1.0, 1.0, 1.0]);
        let data = renderer.to_vertex_data();
        assert_eq!(data.len(), 11);
    }
    #[test]
    fn test_point_sprite_stretch_increases_size() {
        let config = PointSpriteConfig {
            size_world: 0.1,
            use_velocity_stretch: true,
            stretch_factor: 3.0,
        };
        let mut r = PointSpriteRenderer::new(config);
        r.add([0.0; 3], [1.0, 0.0, 0.0], [1.0; 4]);
        let data = r.to_vertex_data();
        let size = data[10];
        assert!((size - 3.1).abs() < 1e-4, "stretched size = {size}");
    }
    #[test]
    fn test_metaball_field_empty() {
        let field = MetaballField::new(vec![], 1.0);
        let v = field.evaluate([0.0, 0.0, 0.0]);
        assert!((v - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_metaball_field_single_blob() {
        let field = MetaballField::new(
            vec![Metaball {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
                strength: 1.0,
            }],
            1.0,
        );
        let v_center = field.evaluate([0.0, 0.0, 0.0]);
        let v_far = field.evaluate([10.0, 0.0, 0.0]);
        assert!(v_center > v_far, "center should have higher field value");
    }
    #[test]
    fn test_metaball_isosurface_scan_crossing() {
        let field = MetaballField::new(
            vec![Metaball {
                center: [0.5, 0.5, 0.5],
                radius: 0.3,
                strength: 2.0,
            }],
            0.5,
        );
        // Field = strength/(d²+ε). Threshold crossing at d≈2.0 from center,
        // so scan range must extend beyond [0,1].
        let crossings = field.scan_x(0.5, 0.5, -3.0, 4.0, 256);
        assert!(
            !crossings.is_empty(),
            "metaball should produce scan crossings"
        );
    }
    #[test]
    fn test_metaball_two_blobs_merge() {
        let field = MetaballField::new(
            vec![
                Metaball {
                    center: [0.0, 0.0, 0.0],
                    radius: 0.5,
                    strength: 1.0,
                },
                Metaball {
                    center: [0.3, 0.0, 0.0],
                    radius: 0.5,
                    strength: 1.0,
                },
            ],
            1.0,
        );
        let between = field.evaluate([0.15, 0.0, 0.0]);
        let outside = field.evaluate([5.0, 0.0, 0.0]);
        assert!(
            between > outside,
            "between blobs should have higher density"
        );
    }
    #[test]
    fn test_sph_density_empty() {
        let sph = SphDensityField::new(0.1);
        let d = sph.density([0.0, 0.0, 0.0], &[]);
        assert!((d - 0.0).abs() < 1e-10);
    }
    #[test]
    fn test_sph_density_single_particle() {
        let sph = SphDensityField::new(1.0);
        let positions = [[0.0f32, 0.0, 0.0]];
        let masses = [1.0f32];
        let d = sph.density(
            [0.0, 0.0, 0.0],
            &positions
                .iter()
                .zip(masses.iter())
                .map(|(&p, &m)| (p, m))
                .collect::<Vec<_>>(),
        );
        assert!(d > 0.0, "density at particle center should be positive");
    }
    #[test]
    fn test_sph_density_falls_off_with_distance() {
        let sph = SphDensityField::new(1.0);
        let particles: Vec<([f32; 3], f32)> = vec![([0.0, 0.0, 0.0], 1.0)];
        let d_near = sph.density([0.1, 0.0, 0.0], &particles);
        let d_far = sph.density([2.0, 0.0, 0.0], &particles);
        assert!(d_near > d_far, "density should fall off with distance");
    }
    #[test]
    fn test_lod_level_close() {
        let mut lod = LodSystem::new();
        lod.add_level(1.0, 0);
        lod.add_level(5.0, 1);
        lod.add_level(20.0, 2);
        assert_eq!(lod.select(0.5), Some(0));
    }
    #[test]
    fn test_lod_level_mid() {
        let mut lod = LodSystem::new();
        lod.add_level(1.0, 0);
        lod.add_level(5.0, 1);
        lod.add_level(20.0, 2);
        assert_eq!(lod.select(3.0), Some(1));
    }
    #[test]
    fn test_lod_level_far() {
        let mut lod = LodSystem::new();
        lod.add_level(1.0, 0);
        lod.add_level(5.0, 1);
        lod.add_level(20.0, 2);
        assert_eq!(lod.select(10.0), Some(2));
    }
    #[test]
    fn test_lod_level_very_far() {
        let mut lod = LodSystem::new();
        lod.add_level(1.0, 0);
        lod.add_level(5.0, 1);
        lod.add_level(20.0, 2);
        let level = lod.select(100.0);
        assert_eq!(level, Some(2), "beyond all thresholds → highest LOD level");
    }
    #[test]
    fn test_lod_render_radius_decreases_with_distance() {
        let mut lod = LodSystem::new();
        lod.add_level(1.0, 0);
        lod.add_level(5.0, 1);
        lod.add_level(20.0, 2);
        let close_level = lod.select(0.5);
        let far_level = lod.select(50.0);
        // Close returns mesh 0 (highest detail), far returns mesh 2 (lowest detail)
        assert!(close_level.is_some());
        assert!(far_level.is_some());
        assert!(
            close_level.unwrap_or(0) < far_level.unwrap_or(0),
            "close particles should use lower LOD index (higher detail)"
        );
    }
    #[test]
    fn test_color_by_temperature() {
        let mut buf = ParticleBuffer::new(3);
        buf.instances.push(ParticleInstance {
            position: [0.0; 3],
            radius: 0.1,
            color: [1.0; 4],
            velocity: [0.0; 3],
        });
        buf.instances.push(ParticleInstance {
            position: [1.0, 0.0, 0.0],
            radius: 0.1,
            color: [1.0; 4],
            velocity: [0.0; 3],
        });
        let temps = vec![300.0f32, 1000.0f32];
        color_by_property(&mut buf, &temps, 300.0, 1000.0, ColorMap::Grayscale);
        assert!(
            buf.instances[0].color[0] < buf.instances[1].color[0],
            "hot particle should be brighter"
        );
    }
    #[test]
    fn test_color_by_property_clamps() {
        let mut buf = ParticleBuffer::new(1);
        buf.instances.push(ParticleInstance {
            position: [0.0; 3],
            radius: 0.1,
            color: [0.5; 4],
            velocity: [0.0; 3],
        });
        let props = vec![1e10f32];
        color_by_property(&mut buf, &props, 0.0, 100.0, ColorMap::Grayscale);
        assert!(
            (buf.instances[0].color[0] - 1.0).abs() < 1e-3,
            "out-of-range value should clamp to 1.0"
        );
    }
    #[test]
    fn test_velocity_glyph_count() {
        let positions = vec![[0.0f32; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let velocities = vec![[1.0f32, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];
        let glyphs = build_velocity_glyphs(&positions, &velocities, 0.5);
        assert_eq!(glyphs.len(), 3);
    }
    #[test]
    fn test_velocity_glyph_tip_position() {
        let positions = vec![[0.0f32, 0.0, 0.0]];
        let velocities = vec![[2.0f32, 0.0, 0.0]];
        let glyphs = build_velocity_glyphs(&positions, &velocities, 1.0);
        let g = &glyphs[0];
        assert!((g.tip[0] - 2.0).abs() < 1e-4, "tip x = {}", g.tip[0]);
        assert!(g.tip[1].abs() < 1e-4);
    }
    #[test]
    fn test_velocity_glyph_zero_velocity() {
        let positions = vec![[1.0f32, 2.0, 3.0]];
        let velocities = vec![[0.0f32; 3]];
        let glyphs = build_velocity_glyphs(&positions, &velocities, 1.0);
        let g = &glyphs[0];
        let diff = (g.tip[0] - g.base[0]).abs()
            + (g.tip[1] - g.base[1]).abs()
            + (g.tip[2] - g.base[2]).abs();
        assert!(diff < 1e-4, "zero-velocity tip should equal base");
    }
    #[test]
    fn test_sort_buffer_depth_ordering() {
        let camera = [0.0f64, 0.0, 0.0];
        let positions = vec![[1.0f64, 0.0, 0.0], [3.0f64, 0.0, 0.0], [2.0f64, 0.0, 0.0]];
        let mut buf = ParticleSortBuffer::new();
        buf.sort_back_to_front(camera, &positions);
        assert_eq!(buf.indices[0], 1, "farthest particle should come first");
        assert_eq!(buf.indices[2], 0, "nearest particle should come last");
    }
    #[test]
    fn test_sort_buffer_single_particle() {
        let mut buf = ParticleSortBuffer::new();
        buf.sort_back_to_front([0.0; 3], &[[5.0, 0.0, 0.0]]);
        assert_eq!(buf.indices.len(), 1);
        assert_eq!(buf.indices[0], 0);
    }
    #[test]
    fn test_gaussian_splat_at_center_is_one() {
        let splat = SplatRenderer {
            sigma: 1.0,
            cutoff_radius: 3.0,
        };
        let w = splat.splat_to_buffer([0.0; 3], 1.0, [0.5, 0.5], 1.0);
        assert!(
            (w - 1.0).abs() < 1e-9,
            "Gaussian at center should be 1.0, got {w}"
        );
    }
    #[test]
    fn test_gaussian_splat_falls_off_with_distance() {
        let splat = SplatRenderer {
            sigma: 0.5,
            cutoff_radius: 3.0,
        };
        let w_center = splat.splat_to_buffer([0.0; 3], 1.0, [0.5, 0.5], 1.0);
        let w_edge = splat.splat_to_buffer([0.0; 3], 1.0, [0.9, 0.5], 1.0);
        assert!(
            w_center > w_edge,
            "Gaussian weight should decrease with distance from centre"
        );
    }
    #[test]
    fn test_frustum_cull_sphere_outside_plane() {
        let planes: [[f64; 4]; 6] = [
            [0.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0, 100.0],
            [-1.0, 0.0, 0.0, 100.0],
            [0.0, 1.0, 0.0, 100.0],
            [0.0, -1.0, 0.0, 100.0],
            [0.0, 0.0, -1.0, 100.0],
        ];
        let positions = vec![[0.0f64, 0.0, -5.0]];
        let radii = vec![1.0f64];
        let visible = frustum_cull_particles(&positions, &radii, &planes);
        assert!(!visible[0], "sphere fully behind plane should be culled");
    }
    #[test]
    fn test_frustum_cull_sphere_inside_all_planes() {
        let planes: [[f64; 4]; 6] = [
            [1.0, 0.0, 0.0, 1000.0],
            [-1.0, 0.0, 0.0, 1000.0],
            [0.0, 1.0, 0.0, 1000.0],
            [0.0, -1.0, 0.0, 1000.0],
            [0.0, 0.0, 1.0, 1000.0],
            [0.0, 0.0, -1.0, 1000.0],
        ];
        let positions = vec![[0.0f64; 3]];
        let radii = vec![1.0f64];
        let visible = frustum_cull_particles(&positions, &radii, &planes);
        assert!(visible[0], "sphere inside frustum should not be culled");
    }
    #[test]
    fn test_metaball_field_add_particle_and_sample() {
        let field = MetaballField::new(
            vec![Metaball {
                center: [0.0, 0.0, 0.0],
                radius: 1.0,
                strength: 1.0,
            }],
            1.0,
        );
        let at_center = field.evaluate([0.0, 0.0, 0.0]);
        let far_away = field.evaluate([5.0, 0.0, 0.0]);
        assert!(
            at_center > far_away,
            "evaluate should return higher value at blob centre ({at_center} vs {far_away})"
        );
    }
    #[test]
    fn test_particle_lod_close_returns_zero() {
        let levels = vec![5.0f64, 15.0, 50.0];
        assert_eq!(particle_lod(1.0, &levels), 0);
    }
    #[test]
    fn test_particle_lod_far_returns_last() {
        let levels = vec![5.0f64, 15.0, 50.0];
        assert_eq!(particle_lod(100.0, &levels), 2);
    }
    #[test]
    fn test_sphere_impostor_vertex_count() {
        let verts = generate_sphere_impostor(1.0, 8);
        assert_eq!(verts.len(), 8);
    }
    #[test]
    fn test_sphere_impostor_radius() {
        let r = 2.5;
        let verts = generate_sphere_impostor(r, 12);
        for v in &verts {
            let dist = ((v[0] as f64).powi(2) + (v[1] as f64).powi(2)).sqrt();
            assert!(
                (dist - r).abs() < 1e-5,
                "impostor vertex distance from origin should equal radius, got {dist}"
            );
        }
    }
}
