//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{ParticleBuffer, ParticleRenderData, SimpleRng, SortedParticleRenderData};

/// Extract rendering data from a particle buffer.
///
/// Returns only alive particles, with color interpolated based on age.
pub fn extract_render_data(
    buffer: &ParticleBuffer,
    color_young: [f32; 4],
    color_old: [f32; 4],
    base_size: f32,
) -> Vec<ParticleRenderData> {
    let mut data = Vec::new();
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let age = buffer.ages[i];
        let initial_lifetime = age + buffer.lifetimes[i];
        let age_norm = if initial_lifetime > 0.0 {
            (age / initial_lifetime).clamp(0.0, 1.0)
        } else {
            1.0
        };
        let color = [
            color_young[0] + (color_old[0] - color_young[0]) * age_norm,
            color_young[1] + (color_old[1] - color_young[1]) * age_norm,
            color_young[2] + (color_old[2] - color_young[2]) * age_norm,
            color_young[3] + (color_old[3] - color_young[3]) * age_norm,
        ];
        let size = base_size * (1.0 - 0.5 * age_norm);
        data.push(ParticleRenderData {
            position: buffer.get_position(i),
            color,
            size,
            age_normalized: age_norm,
        });
    }
    data
}
/// Compute total momentum of all alive particles: sum(m * v).
pub fn compute_total_momentum(buffer: &ParticleBuffer) -> [f32; 3] {
    let mut px = 0.0f32;
    let mut py = 0.0f32;
    let mut pz = 0.0f32;
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let m = buffer.masses[i];
        px += m * buffer.velocities_x[i];
        py += m * buffer.velocities_y[i];
        pz += m * buffer.velocities_z[i];
    }
    [px, py, pz]
}
/// Morton code (Z-order curve) interleave for a 3D position.
///
/// Quantizes each coordinate to 10 bits, then interleaves.
pub fn morton_encode(x: u32, y: u32, z: u32) -> u32 {
    let xm = morton_part1by2(x);
    let ym = morton_part1by2(y);
    let zm = morton_part1by2(z);
    xm | (ym << 1) | (zm << 2)
}
pub(super) fn morton_part1by2(mut x: u32) -> u32 {
    x &= 0x000003ff;
    x = (x ^ (x << 16)) & 0xff0000ff;
    x = (x ^ (x << 8)) & 0x0300f00f;
    x = (x ^ (x << 4)) & 0x030c30c3;
    x = (x ^ (x << 2)) & 0x09249249;
    x
}
/// Compute Morton codes for all alive particles in a buffer.
///
/// Particles are sorted into a `[0, grid_cells)^3` grid using their positions.
pub fn compute_morton_codes(
    buffer: &ParticleBuffer,
    origin: [f32; 3],
    extent: [f32; 3],
    grid_cells: u32,
) -> Vec<(u32, usize)> {
    let mut codes: Vec<(u32, usize)> = Vec::new();
    let gc = grid_cells as f32;
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let px = ((buffer.positions_x[i] - origin[0]) / extent[0] * gc).clamp(0.0, gc - 1.0) as u32;
        let py = ((buffer.positions_y[i] - origin[1]) / extent[1] * gc).clamp(0.0, gc - 1.0) as u32;
        let pz = ((buffer.positions_z[i] - origin[2]) / extent[2] * gc).clamp(0.0, gc - 1.0) as u32;
        codes.push((morton_encode(px, py, pz), i));
    }
    codes.sort_unstable_by_key(|&(code, _)| code);
    codes
}
/// Reorder a particle buffer according to Morton-sorted indices.
///
/// Returns a new buffer with particles ordered by z-curve for better
/// cache locality.
pub fn sort_particles_morton(
    buffer: &ParticleBuffer,
    origin: [f32; 3],
    extent: [f32; 3],
    grid_cells: u32,
) -> ParticleBuffer {
    let sorted = compute_morton_codes(buffer, origin, extent, grid_cells);
    let mut new_buf = ParticleBuffer::new(buffer.count);
    for (slot, &(_, old_idx)) in sorted.iter().enumerate() {
        new_buf.positions_x[slot] = buffer.positions_x[old_idx];
        new_buf.positions_y[slot] = buffer.positions_y[old_idx];
        new_buf.positions_z[slot] = buffer.positions_z[old_idx];
        new_buf.velocities_x[slot] = buffer.velocities_x[old_idx];
        new_buf.velocities_y[slot] = buffer.velocities_y[old_idx];
        new_buf.velocities_z[slot] = buffer.velocities_z[old_idx];
        new_buf.masses[slot] = buffer.masses[old_idx];
        new_buf.lifetimes[slot] = buffer.lifetimes[old_idx];
        new_buf.ages[slot] = buffer.ages[old_idx];
    }
    new_buf
}
/// Prepare GPU particle rendering data sorted back-to-front for transparency.
///
/// `camera_forward` is the camera's forward direction (normalized).
pub fn prepare_sorted_render_data(
    buffer: &ParticleBuffer,
    color_young: [f32; 4],
    color_old: [f32; 4],
    base_size: f32,
    camera_pos: [f32; 3],
    camera_forward: [f32; 3],
) -> Vec<SortedParticleRenderData> {
    let mut result = Vec::new();
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let age = buffer.ages[i];
        let initial_lifetime = age + buffer.lifetimes[i];
        let age_norm = if initial_lifetime > 0.0 {
            (age / initial_lifetime).clamp(0.0, 1.0)
        } else {
            1.0
        };
        let color = [
            color_young[0] + (color_old[0] - color_young[0]) * age_norm,
            color_young[1] + (color_old[1] - color_young[1]) * age_norm,
            color_young[2] + (color_old[2] - color_young[2]) * age_norm,
            color_young[3] + (color_old[3] - color_young[3]) * age_norm,
        ];
        let size = base_size * (1.0 - 0.5 * age_norm);
        let dx = buffer.positions_x[i] - camera_pos[0];
        let dy = buffer.positions_y[i] - camera_pos[1];
        let dz = buffer.positions_z[i] - camera_pos[2];
        let depth = dx * camera_forward[0] + dy * camera_forward[1] + dz * camera_forward[2];
        result.push(SortedParticleRenderData {
            render_data: ParticleRenderData {
                position: [
                    buffer.positions_x[i],
                    buffer.positions_y[i],
                    buffer.positions_z[i],
                ],
                color,
                size,
                age_normalized: age_norm,
            },
            sort_key: depth,
            buffer_index: i,
        });
    }
    result.sort_unstable_by(|a, b| {
        b.sort_key
            .partial_cmp(&a.sort_key)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    result
}
/// Compute center of mass of all alive particles.
pub fn compute_center_of_mass(buffer: &ParticleBuffer) -> [f32; 3] {
    let mut total_mass = 0.0f32;
    let mut cx = 0.0f32;
    let mut cy = 0.0f32;
    let mut cz = 0.0f32;
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let m = buffer.masses[i];
        cx += m * buffer.positions_x[i];
        cy += m * buffer.positions_y[i];
        cz += m * buffer.positions_z[i];
        total_mass += m;
    }
    if total_mass > 0.0 {
        [cx / total_mass, cy / total_mass, cz / total_mass]
    } else {
        [0.0, 0.0, 0.0]
    }
}
/// Compute a velocity-speed histogram for all alive particles.
///
/// Speed bins are `[0, dv)`, `[dv, 2*dv)`, … up to `max_speed`.
/// Returns a `Vec`usize` of length `ceil(max_speed / dv)`.
pub fn compute_velocity_histogram(buffer: &ParticleBuffer, max_speed: f32, dv: f32) -> Vec<usize> {
    let dv = dv.max(1e-12);
    let n_bins = (max_speed / dv).ceil() as usize;
    let mut hist = vec![0usize; n_bins.max(1)];
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let vx = buffer.velocities_x[i];
        let vy = buffer.velocities_y[i];
        let vz = buffer.velocities_z[i];
        let speed = (vx * vx + vy * vy + vz * vz).sqrt();
        if speed < max_speed {
            let bin = (speed / dv) as usize;
            if bin < n_bins {
                hist[bin] += 1;
            }
        }
    }
    hist
}
/// Compute the total angular momentum of all alive particles about the origin.
///
/// `L = sum_i m_i * (r_i × v_i)`
///
/// Returns `\[Lx, Ly, Lz\]`.
pub fn compute_angular_momentum(buffer: &ParticleBuffer) -> [f32; 3] {
    let mut lx = 0.0f32;
    let mut ly = 0.0f32;
    let mut lz = 0.0f32;
    for i in 0..buffer.count {
        if !buffer.is_alive(i) {
            continue;
        }
        let m = buffer.masses[i];
        let rx = buffer.positions_x[i];
        let ry = buffer.positions_y[i];
        let rz = buffer.positions_z[i];
        let vx = buffer.velocities_x[i];
        let vy = buffer.velocities_y[i];
        let vz = buffer.velocities_z[i];
        lx += m * (ry * vz - rz * vy);
        ly += m * (rz * vx - rx * vz);
        lz += m * (rx * vy - ry * vx);
    }
    [lx, ly, lz]
}
/// Emit a burst of `count` particles from `origin` with the given velocity
/// and lifetime, using an internal LCG for per-particle spread.
///
/// Spread magnitude `spread` adds a random offset to velocity.
/// Returns the number of particles actually emitted (may be less than `count`
/// if the buffer is full).
pub fn emit_burst(
    buffer: &mut ParticleBuffer,
    origin: [f32; 3],
    velocity: [f32; 3],
    spread: f32,
    lifetime: f32,
    mass: f32,
    count: usize,
    seed: u64,
) -> usize {
    let mut rng = SimpleRng::new(seed);
    let mut spawned = 0usize;
    for _ in 0..count {
        let dir = rng.next_unit_sphere();
        let vel = [
            velocity[0] + dir[0] * spread,
            velocity[1] + dir[1] * spread,
            velocity[2] + dir[2] * spread,
        ];
        if buffer.add_particle(origin, vel, mass, lifetime).is_some() {
            spawned += 1;
        }
    }
    spawned
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::BoundingBoxKill;
    use crate::DragForce;
    use crate::EmitterShape;
    use crate::FloorCollision;
    use crate::GpuParticleEmitter;
    use crate::GpuParticleLayout;
    use crate::GravityForce;
    use crate::GridParticleCollision;
    use crate::ParticleEmitter;
    use crate::ParticleIntegrator;
    use crate::ParticleLifetimeManager;
    use crate::ParticleRepulsion;
    use crate::ParticleStats;
    use crate::ParticleSystem;
    use crate::ParticleSystemStats;
    use crate::RadialForceField;
    use crate::VortexForceField;
    #[test]
    fn test_particle_buffer_add_and_get_position() {
        let mut buf = ParticleBuffer::new(4);
        let idx = buf.add_particle([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], 1.0, 5.0);
        assert!(idx.is_some());
        let i = idx.unwrap();
        let pos = buf.get_position(i);
        assert!((pos[0] - 1.0).abs() < 1e-6);
        assert!((pos[1] - 2.0).abs() < 1e-6);
        assert!((pos[2] - 3.0).abs() < 1e-6);
    }
    #[test]
    fn test_particle_buffer_is_alive_after_kill() {
        let mut buf = ParticleBuffer::new(4);
        let i = buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        assert!(buf.is_alive(i));
        buf.kill(i);
        assert!(!buf.is_alive(i));
    }
    #[test]
    fn test_gravity_force_increases_downward_velocity() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([0.0; 3], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let grav = GravityForce {
            g: [0.0, -9.81, 0.0],
        };
        grav.apply(&mut buf, 1.0);
        let vel = buf.get_velocity(0);
        assert!(vel[1] < 0.0, "vy should be negative after gravity");
        assert!((vel[1] + 9.81).abs() < 1e-4);
    }
    #[test]
    fn test_integrator_moves_particles() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        ParticleIntegrator::integrate(&mut buf, 1.0);
        let pos = buf.get_position(0);
        assert!(
            (pos[0] - 1.0).abs() < 1e-6,
            "particle should move 1 unit in x"
        );
    }
    #[test]
    fn test_floor_collision_reflects_particle() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([0.0, -0.5, 0.0], [0.0, -3.0, 0.0], 1.0, 10.0)
            .unwrap();
        let floor = FloorCollision {
            y: 0.0,
            restitution: 0.8,
        };
        floor.apply(&mut buf);
        let pos = buf.get_position(0);
        let vel = buf.get_velocity(0);
        assert!(pos[1] >= 0.0, "particle should be at or above floor");
        assert!(
            vel[1] > 0.0,
            "velocity y should be positive after reflection"
        );
        assert!(
            (vel[1] - 2.4).abs() < 1e-5,
            "reflected vy = 3.0 * 0.8 = 2.4"
        );
    }
    #[test]
    fn test_bounding_box_kill_removes_out_of_bounds() {
        let mut buf = ParticleBuffer::new(3);
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 10.0)
            .unwrap();
        buf.add_particle([100.0, 0.0, 0.0], [0.0; 3], 1.0, 10.0)
            .unwrap();
        buf.add_particle([0.0, -50.0, 0.0], [0.0; 3], 1.0, 10.0)
            .unwrap();
        let kill = BoundingBoxKill {
            min: [-10.0; 3],
            max: [10.0; 3],
        };
        kill.apply(&mut buf);
        assert!(buf.is_alive(0), "particle inside box should survive");
        assert!(!buf.is_alive(1), "particle outside x should be killed");
        assert!(!buf.is_alive(2), "particle outside y should be killed");
    }
    #[test]
    fn test_particle_system_step_no_panic() {
        let mut sys = ParticleSystem::new(64);
        let emitter = ParticleEmitter::new([0.0; 3], 10.0, [0.0, 1.0, 0.0], 3.0);
        sys.add_emitter(emitter);
        sys.floor = Some(FloorCollision {
            y: -5.0,
            restitution: 0.5,
        });
        for _ in 0..30 {
            sys.step(1.0 / 60.0);
        }
    }
    #[test]
    fn test_gpu_particle_layout_round_trip() {
        let mut buf = ParticleBuffer::new(4);
        buf.add_particle([1.0, 2.0, 3.0], [4.0, 5.0, 6.0], 0.5, 7.0)
            .unwrap();
        buf.add_particle([-1.0, 0.5, 2.0], [0.1, 0.2, 0.3], 2.0, 3.5)
            .unwrap();
        let flat = GpuParticleLayout::to_f32_buffer(&buf);
        assert_eq!(flat.len(), 4 * GpuParticleLayout::stride());
        let restored = GpuParticleLayout::from_f32_buffer(&flat, 4);
        let p0 = restored.get_position(0);
        assert!((p0[0] - 1.0).abs() < 1e-6);
        assert!((p0[1] - 2.0).abs() < 1e-6);
        assert!((p0[2] - 3.0).abs() < 1e-6);
        let v1 = restored.get_velocity(1);
        assert!((v1[0] - 0.1).abs() < 1e-6);
        assert!((v1[1] - 0.2).abs() < 1e-6);
        assert!((v1[2] - 0.3).abs() < 1e-6);
        assert!(!restored.is_alive(2));
        assert!(!restored.is_alive(3));
    }
    #[test]
    fn test_radial_force_attraction() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([5.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let field = RadialForceField {
            center: [0.0, 0.0, 0.0],
            strength: 10.0,
            falloff: 1.0,
            min_distance: 0.01,
        };
        field.apply(&mut buf, 1.0);
        let vel = buf.get_velocity(0);
        assert!(
            vel[0] < 0.0,
            "should attract toward center, got vx={}",
            vel[0]
        );
    }
    #[test]
    fn test_radial_force_repulsion() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([5.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let field = RadialForceField {
            center: [0.0, 0.0, 0.0],
            strength: -10.0,
            falloff: 1.0,
            min_distance: 0.01,
        };
        field.apply(&mut buf, 1.0);
        let vel = buf.get_velocity(0);
        assert!(vel[0] > 0.0, "should repel from center, got vx={}", vel[0]);
    }
    #[test]
    fn test_radial_force_dead_particle_ignored() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([5.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        buf.add_particle([3.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        buf.kill(1);
        let field = RadialForceField {
            center: [0.0, 0.0, 0.0],
            strength: 10.0,
            falloff: 1.0,
            min_distance: 0.01,
        };
        field.apply(&mut buf, 1.0);
        let vel1 = buf.get_velocity(1);
        assert!((vel1[0]).abs() < 1e-6);
    }
    #[test]
    fn test_vortex_force_field() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let vortex = VortexForceField {
            center: [0.0, 0.0],
            angular_velocity: 10.0,
            radius: 5.0,
        };
        vortex.apply(&mut buf, 1.0);
        let vel = buf.get_velocity(0);
        assert!(
            vel[2].abs() > 0.0 || vel[0].abs() > 0.0,
            "vortex should add tangential velocity"
        );
    }
    #[test]
    fn test_vortex_outside_radius() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([100.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let vortex = VortexForceField {
            center: [0.0, 0.0],
            angular_velocity: 10.0,
            radius: 5.0,
        };
        vortex.apply(&mut buf, 1.0);
        let vel = buf.get_velocity(0);
        assert!((vel[0]).abs() < 1e-6);
        assert!((vel[2]).abs() < 1e-6);
    }
    #[test]
    fn test_particle_repulsion_two_particles() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        buf.add_particle([0.5, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let repulsion = ParticleRepulsion {
            strength: 10.0,
            radius: 1.0,
        };
        repulsion.apply(&mut buf, 1.0);
        let v0 = buf.get_velocity(0);
        let v1 = buf.get_velocity(1);
        assert!(v0[0] < 0.0, "particle 0 should move left");
        assert!(v1[0] > 0.0, "particle 1 should move right");
        let total = v0[0] + v1[0];
        assert!(total.abs() < 1e-5, "momentum not conserved: {total}");
    }
    #[test]
    fn test_particle_repulsion_outside_radius() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        buf.add_particle([10.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let repulsion = ParticleRepulsion {
            strength: 10.0,
            radius: 1.0,
        };
        repulsion.apply(&mut buf, 1.0);
        let v0 = buf.get_velocity(0);
        let v1 = buf.get_velocity(1);
        assert!((v0[0]).abs() < 1e-6);
        assert!((v1[0]).abs() < 1e-6);
    }
    #[test]
    fn test_extract_render_data_empty_buffer() {
        let buf = ParticleBuffer::new(4);
        let data = extract_render_data(&buf, [1.0; 4], [0.0; 4], 1.0);
        assert!(data.is_empty());
    }
    #[test]
    fn test_extract_render_data_alive_only() {
        let mut buf = ParticleBuffer::new(4);
        buf.add_particle([1.0, 2.0, 3.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([4.0, 5.0, 6.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let data = extract_render_data(&buf, [1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0], 2.0);
        assert_eq!(data.len(), 2);
        assert!((data[0].position[0] - 1.0).abs() < 1e-6);
        assert!(data[0].size > 0.0);
        assert!(data[0].age_normalized >= 0.0 && data[0].age_normalized <= 1.0);
    }
    #[test]
    fn test_extract_render_data_color_interpolation() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 2.0).unwrap();
        buf.ages[0] = 1.0;
        buf.lifetimes[0] = 1.0;
        let data = extract_render_data(&buf, [1.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0, 0.0], 1.0);
        assert_eq!(data.len(), 1);
        assert!((data[0].color[0] - 0.5).abs() < 1e-4);
        assert!((data[0].color[2] - 0.5).abs() < 1e-4);
    }
    #[test]
    fn test_compute_total_momentum() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0; 3], [1.0, 0.0, 0.0], 2.0, 10.0)
            .unwrap();
        buf.add_particle([0.0; 3], [-1.0, 0.0, 0.0], 3.0, 10.0)
            .unwrap();
        let p = compute_total_momentum(&buf);
        assert!((p[0] - (-1.0)).abs() < 1e-5);
    }
    #[test]
    fn test_compute_total_momentum_dead_ignored() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0; 3], [1.0, 0.0, 0.0], 2.0, 10.0)
            .unwrap();
        buf.add_particle([0.0; 3], [10.0, 0.0, 0.0], 3.0, 10.0)
            .unwrap();
        buf.kill(1);
        let p = compute_total_momentum(&buf);
        assert!((p[0] - 2.0).abs() < 1e-5);
    }
    #[test]
    fn test_compute_center_of_mass() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 10.0)
            .unwrap();
        buf.add_particle([10.0, 0.0, 0.0], [0.0; 3], 1.0, 10.0)
            .unwrap();
        let com = compute_center_of_mass(&buf);
        assert!((com[0] - 5.0).abs() < 1e-5);
    }
    #[test]
    fn test_compute_center_of_mass_weighted() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 10.0)
            .unwrap();
        buf.add_particle([10.0, 0.0, 0.0], [0.0; 3], 3.0, 10.0)
            .unwrap();
        let com = compute_center_of_mass(&buf);
        assert!((com[0] - 7.5).abs() < 1e-5);
    }
    #[test]
    fn test_compute_center_of_mass_empty() {
        let buf = ParticleBuffer::new(4);
        let com = compute_center_of_mass(&buf);
        assert!((com[0]).abs() < 1e-6);
    }
    #[test]
    fn test_particle_stats_basic() {
        let mut buf = ParticleBuffer::new(3);
        buf.add_particle([1.0, 2.0, 3.0], [1.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([4.0, 5.0, 6.0], [0.0, 2.0, 0.0], 2.0, 5.0)
            .unwrap();
        let stats = ParticleStats::compute(&buf);
        assert_eq!(stats.active, 2);
        assert!(stats.avg_speed > 0.0);
        assert!(stats.total_kinetic_energy > 0.0);
    }
    #[test]
    fn test_particle_stats_no_alive() {
        let buf = ParticleBuffer::new(4);
        let stats = ParticleStats::compute(&buf);
        assert_eq!(stats.active, 0);
        assert!((stats.avg_speed).abs() < 1e-6);
    }
    #[test]
    fn test_emitter_inactive_no_emission() {
        let mut emitter = ParticleEmitter::new([0.0; 3], 1000.0, [0.0; 3], 1.0);
        emitter.active = false;
        let mut buf = ParticleBuffer::new(100);
        let spawned = emitter.emit(&mut buf, 1.0, 42);
        assert_eq!(spawned, 0);
    }
    #[test]
    fn test_emitter_emits_particles() {
        let mut emitter = ParticleEmitter::new([0.0; 3], 100.0, [0.0, 1.0, 0.0], 5.0);
        let mut buf = ParticleBuffer::new(200);
        let spawned = emitter.emit(&mut buf, 1.0, 42);
        assert!(spawned > 0, "should emit particles");
    }
    #[test]
    fn test_emitter_with_spread() {
        let mut emitter = ParticleEmitter::new([0.0; 3], 10.0, [0.0, 1.0, 0.0], 5.0);
        emitter.velocity_spread = 0.5;
        let mut buf = ParticleBuffer::new(20);
        emitter.emit(&mut buf, 1.0, 42);
        assert!(buf.active_count() > 0);
    }
    #[test]
    fn test_drag_reduces_speed() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([0.0; 3], [10.0, 0.0, 0.0], 1.0, 10.0)
            .unwrap();
        let drag = DragForce { coefficient: 0.5 };
        drag.apply(&mut buf, 1.0);
        let vel = buf.get_velocity(0);
        assert!(vel[0] < 10.0, "drag should reduce speed");
        assert!(vel[0] > 0.0, "speed should stay positive");
    }
    #[test]
    fn test_active_count() {
        let mut buf = ParticleBuffer::new(5);
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        assert_eq!(buf.active_count(), 3);
        buf.kill(1);
        assert_eq!(buf.active_count(), 2);
    }
    #[test]
    fn test_gpu_emitter_continuous_burst_count_zero() {
        let emitter = GpuParticleEmitter::new_continuous([0.0; 3], 50.0, 2.0);
        assert_eq!(emitter.burst_count(), 0);
    }
    #[test]
    fn test_gpu_emitter_burst_count() {
        let emitter = GpuParticleEmitter::new_burst([0.0; 3], 100, 2.0);
        assert_eq!(emitter.burst_count(), 100);
    }
    #[test]
    fn test_gpu_emitter_burst_fires_once() {
        let mut emitter = GpuParticleEmitter::new_burst([0.0; 3], 5, 1.5);
        let mut buf = ParticleBuffer::new(32);
        let first = emitter.emit(&mut buf, 0.016);
        assert_eq!(first, 5);
        let second = emitter.emit(&mut buf, 0.016);
        assert_eq!(second, 0);
    }
    #[test]
    fn test_gpu_emitter_continuous_accumulates() {
        let mut emitter = GpuParticleEmitter::new_continuous([0.0; 3], 100.0, 2.0);
        let mut buf = ParticleBuffer::new(64);
        let spawned = emitter.emit(&mut buf, 0.1);
        assert_eq!(spawned, 10);
    }
    #[test]
    fn test_gpu_emitter_sphere_shape_positions_in_radius() {
        let mut emitter = GpuParticleEmitter::new_continuous([0.0; 3], 100.0, 2.0);
        emitter.shape = EmitterShape::Sphere { radius: 3.0 };
        let mut buf = ParticleBuffer::new(64);
        emitter.emit(&mut buf, 1.0);
        let count = buf.active_count();
        assert!(count > 0);
        for i in 0..count {
            let px = buf.positions_x[i];
            let py = buf.positions_y[i];
            let pz = buf.positions_z[i];
            let r = (px * px + py * py + pz * pz).sqrt();
            assert!(r <= 3.0 + 1e-3, "sphere sample outside radius: r={r}");
        }
    }
    #[test]
    fn test_gpu_emitter_box_shape_positions_in_bounds() {
        let mut emitter = GpuParticleEmitter::new_continuous([0.0; 3], 50.0, 2.0);
        emitter.shape = EmitterShape::Box {
            half_extents: [1.0, 2.0, 0.5],
        };
        let mut buf = ParticleBuffer::new(64);
        emitter.emit(&mut buf, 1.0);
        let count = buf.active_count();
        assert!(count > 0);
        for i in 0..count {
            assert!(buf.positions_x[i].abs() <= 1.0 + 1e-3);
            assert!(buf.positions_y[i].abs() <= 2.0 + 1e-3);
            assert!(buf.positions_z[i].abs() <= 0.5 + 1e-3);
        }
    }
    #[test]
    fn test_lifetime_manager_no_spawn_fraction() {
        let mgr = ParticleLifetimeManager::new();
        let buf = ParticleBuffer::new(4);
        let frac = mgr.alive_fraction(&buf);
        assert!((frac).abs() < 1e-5, "no spawns → fraction=0");
    }
    #[test]
    fn test_lifetime_manager_spawn_records_count() {
        let mut mgr = ParticleLifetimeManager::new();
        mgr.record_spawn(1.0);
        mgr.record_spawn(2.0);
        mgr.record_spawn(3.0);
        assert_eq!(mgr.total_spawned, 3);
        assert!((mgr.min_observed_lifetime - 1.0).abs() < 1e-5);
        assert!((mgr.max_observed_lifetime - 3.0).abs() < 1e-5);
    }
    #[test]
    fn test_lifetime_manager_expiration_count() {
        let mut mgr = ParticleLifetimeManager::new();
        mgr.record_spawn(5.0);
        mgr.record_spawn(5.0);
        mgr.record_expiration();
        assert_eq!(mgr.total_expired, 1);
    }
    #[test]
    fn test_lifetime_manager_alive_fraction_with_active() {
        let mut mgr = ParticleLifetimeManager::new();
        for _ in 0..4 {
            mgr.record_spawn(5.0);
        }
        let mut buf = ParticleBuffer::new(4);
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        buf.add_particle([1.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        let frac = mgr.alive_fraction(&buf);
        assert!((frac - 0.5).abs() < 1e-5, "expected 0.5 got {frac}");
    }
    #[test]
    fn test_morton_encode_zeros() {
        assert_eq!(morton_encode(0, 0, 0), 0);
    }
    #[test]
    fn test_morton_encode_axes_distinct() {
        let cx = morton_encode(1, 0, 0);
        let cy = morton_encode(0, 1, 0);
        let cz = morton_encode(0, 0, 1);
        assert_ne!(cx, 0);
        assert_ne!(cy, 0);
        assert_ne!(cz, 0);
        assert_ne!(cx, cy);
        assert_ne!(cy, cz);
        assert_ne!(cx, cz);
    }
    #[test]
    fn test_compute_morton_codes_length() {
        let mut buf = ParticleBuffer::new(4);
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([1.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 1.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 0.0, 1.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let codes = compute_morton_codes(&buf, [0.0; 3], [2.0; 3], 16);
        assert_eq!(codes.len(), 4);
    }
    #[test]
    fn test_compute_morton_codes_sorted_ascending() {
        let mut buf = ParticleBuffer::new(3);
        buf.add_particle([1.0, 1.0, 1.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.5, 0.5, 0.5], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let codes = compute_morton_codes(&buf, [0.0; 3], [2.0; 3], 16);
        for w in codes.windows(2) {
            assert!(
                w[0].0 <= w[1].0,
                "codes not sorted: {} > {}",
                w[0].0,
                w[1].0
            );
        }
    }
    #[test]
    fn test_sort_particles_morton_origin_first() {
        let mut buf = ParticleBuffer::new(3);
        buf.add_particle([1.0, 1.0, 1.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.5, 0.5, 0.5], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let sorted = sort_particles_morton(&buf, [0.0; 3], [2.0; 3], 16);
        assert!(
            (sorted.positions_x[0]).abs() < 1e-5,
            "origin should be first"
        );
    }
    #[test]
    fn test_grid_collision_no_overlap() {
        let mut buf = ParticleBuffer::new(3);
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([10.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 10.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let col = GridParticleCollision::new(2.0, 0.3, 0.8);
        col.resolve(&mut buf);
        for i in 0..3 {
            let v = buf.get_velocity(i);
            assert!((v[0]).abs() < 1e-5, "particle {i} should have zero vx");
        }
    }
    #[test]
    fn test_grid_collision_overlapping_pair_momentum_conserved() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.5, 0.0, 0.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let col = GridParticleCollision::new(1.0, 0.4, 0.5);
        col.resolve(&mut buf);
        let v0 = buf.get_velocity(0);
        let v1 = buf.get_velocity(1);
        let total_px = v0[0] * buf.masses[0] + v1[0] * buf.masses[1];
        assert!(total_px.abs() < 1e-4, "momentum not conserved: {total_px}");
    }
    #[test]
    fn test_prepare_sorted_render_data_empty() {
        let buf = ParticleBuffer::new(4);
        let data =
            prepare_sorted_render_data(&buf, [1.0; 4], [0.5; 4], 1.0, [0.0; 3], [0.0, 0.0, 1.0]);
        assert!(data.is_empty());
    }
    #[test]
    fn test_prepare_sorted_render_data_back_to_front() {
        let mut buf = ParticleBuffer::new(3);
        buf.add_particle([0.0, 0.0, 1.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 0.0, 5.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0, 0.0, 3.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let camera_fwd = [0.0_f32, 0.0, 1.0];
        let data = prepare_sorted_render_data(&buf, [1.0; 4], [0.0; 4], 1.0, [0.0; 3], camera_fwd);
        assert_eq!(data.len(), 3);
        assert!(
            data[0].sort_key >= data[1].sort_key,
            "first entry should be farthest: {} >= {}",
            data[0].sort_key,
            data[1].sort_key
        );
        assert!(
            data[1].sort_key >= data[2].sort_key,
            "second entry should be middle: {} >= {}",
            data[1].sort_key,
            data[2].sort_key
        );
    }
    #[test]
    fn test_prepare_sorted_render_data_position_preserved() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([1.0, 2.0, 3.0], [0.0; 3], 1.0, 5.0)
            .unwrap();
        let data = prepare_sorted_render_data(
            &buf,
            [1.0_f32, 0.0, 0.0, 1.0],
            [0.0_f32, 1.0, 0.0, 0.5],
            2.0,
            [0.0; 3],
            [0.0, 0.0, 1.0],
        );
        assert_eq!(data.len(), 1);
        assert!((data[0].render_data.position[0] - 1.0).abs() < 1e-5);
        assert!((data[0].render_data.position[1] - 2.0).abs() < 1e-5);
        assert!((data[0].render_data.position[2] - 3.0).abs() < 1e-5);
        assert!(data[0].render_data.size > 0.0);
        assert!(data[0].sort_key.abs() > 0.0);
    }
    #[test]
    fn test_particle_system_stats_extended_empty() {
        let buf = ParticleBuffer::new(8);
        let stats = ParticleSystemStats::compute_extended(&buf);
        assert_eq!(stats.basic.active, 0);
        assert_eq!(stats.capacity, 8);
        assert!(!stats.is_near_capacity(0.9));
    }
    #[test]
    fn test_particle_system_stats_near_capacity() {
        let mut buf = ParticleBuffer::new(4);
        buf.add_particle([0.0; 3], [1.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([1.0; 3], [0.0, 1.0, 0.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([2.0; 3], [0.0, 0.0, 1.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([3.0; 3], [1.0, 1.0, 0.0], 1.0, 5.0)
            .unwrap();
        let stats = ParticleSystemStats::compute_extended(&buf);
        assert!(
            stats.is_near_capacity(0.9),
            "4/4 active should be near capacity"
        );
    }
    #[test]
    fn test_particle_system_stats_fill_ratio() {
        let mut buf = ParticleBuffer::new(4);
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        buf.add_particle([0.0; 3], [0.0; 3], 1.0, 5.0).unwrap();
        let stats = ParticleSystemStats::compute_extended(&buf);
        assert!(
            (stats.fill_ratio - 0.5).abs() < 1e-5,
            "expected 0.5 got {}",
            stats.fill_ratio
        );
    }
    #[test]
    fn test_particle_system_stats_kinetic_energy() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0; 3], [2.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0; 3], [0.0, 4.0, 0.0], 2.0, 5.0)
            .unwrap();
        let stats = ParticleSystemStats::compute_extended(&buf);
        assert!(
            (stats.total_kinetic_energy - 18.0).abs() < 1e-3,
            "expected KE=18 got {}",
            stats.total_kinetic_energy
        );
    }
    #[test]
    fn test_particle_system_stats_mean_age_zero_at_spawn() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0; 3], [1.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0; 3], [0.0, 1.0, 0.0], 1.0, 5.0)
            .unwrap();
        let stats = ParticleSystemStats::compute_extended(&buf);
        assert!(
            (stats.mean_age).abs() < 1e-5,
            "fresh particles should have mean_age=0"
        );
    }
    #[test]
    fn test_velocity_histogram_basic_bin() {
        let mut buf = ParticleBuffer::new(2);
        buf.add_particle([0.0; 3], [1.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        buf.add_particle([0.0; 3], [3.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        let hist = compute_velocity_histogram(&buf, 5.0, 1.0);
        assert!(hist[1] >= 1, "bin 1 should contain the speed=1.0 particle");
        assert!(hist[3] >= 1, "bin 3 should contain the speed=3.0 particle");
    }
    #[test]
    fn test_velocity_histogram_empty_buffer() {
        let buf = ParticleBuffer::new(4);
        let hist = compute_velocity_histogram(&buf, 5.0, 1.0);
        let total: usize = hist.iter().sum();
        assert_eq!(total, 0, "empty buffer yields zero histogram");
    }
    #[test]
    fn test_velocity_histogram_length() {
        let buf = ParticleBuffer::new(1);
        let hist = compute_velocity_histogram(&buf, 10.0, 2.0);
        assert_eq!(hist.len(), 5, "ceil(10/2)=5 bins");
    }
    #[test]
    fn test_angular_momentum_z_axis() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([1.0, 0.0, 0.0], [0.0, 2.0, 0.0], 1.0, 5.0)
            .unwrap();
        let l = compute_angular_momentum(&buf);
        assert!((l[2] - 2.0).abs() < 1e-5, "Lz should be 2, got {}", l[2]);
        assert!(l[0].abs() < 1e-5);
        assert!(l[1].abs() < 1e-5);
    }
    #[test]
    fn test_angular_momentum_zero_for_radial_motion() {
        let mut buf = ParticleBuffer::new(1);
        buf.add_particle([1.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 5.0)
            .unwrap();
        let l = compute_angular_momentum(&buf);
        assert!(l[0].abs() < 1e-6);
        assert!(l[1].abs() < 1e-6);
        assert!(l[2].abs() < 1e-6);
    }
    #[test]
    fn test_angular_momentum_empty_buffer() {
        let buf = ParticleBuffer::new(4);
        let l = compute_angular_momentum(&buf);
        assert!(l[0].abs() < 1e-6 && l[1].abs() < 1e-6 && l[2].abs() < 1e-6);
    }
    #[test]
    fn test_emit_burst_count() {
        let mut buf = ParticleBuffer::new(10);
        let spawned = emit_burst(&mut buf, [0.0; 3], [0.0, 1.0, 0.0], 0.1, 5.0, 1.0, 5, 42);
        assert_eq!(spawned, 5, "should emit exactly 5 particles");
    }
    #[test]
    fn test_emit_burst_respects_buffer_capacity() {
        let mut buf = ParticleBuffer::new(3);
        let spawned = emit_burst(&mut buf, [0.0; 3], [0.0, 1.0, 0.0], 0.0, 5.0, 1.0, 100, 99);
        assert_eq!(spawned, 3, "cannot exceed buffer capacity");
    }
    #[test]
    fn test_emit_burst_particles_are_alive() {
        let mut buf = ParticleBuffer::new(5);
        emit_burst(&mut buf, [1.0, 2.0, 3.0], [0.0; 3], 0.0, 5.0, 1.0, 3, 7);
        let alive: usize = (0..buf.count).filter(|&i| buf.is_alive(i)).count();
        assert_eq!(alive, 3);
    }
}
