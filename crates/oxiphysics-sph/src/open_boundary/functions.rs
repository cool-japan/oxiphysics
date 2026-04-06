//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::TurbulenceInjectionConfig;

/// Uniform velocity profile: constant velocity everywhere.
pub fn uniform_velocity_profile(inflow_velocity: [f64; 3], _pos: [f64; 3]) -> [f64; 3] {
    inflow_velocity
}
/// Parabolic (Poiseuille) velocity profile for pipe flow.
///
/// `v(r) = v_max * (1 - (r/R)²)` where r is the distance from the center axis.
/// Returns the velocity in the direction of `axis`.
pub fn parabolic_velocity_profile(
    center: [f64; 3],
    axis: [f64; 3],
    radius: f64,
    v_max: f64,
    pos: [f64; 3],
) -> [f64; 3] {
    let dx = pos[0] - center[0];
    let dy = pos[1] - center[1];
    let dz = pos[2] - center[2];
    let axis_len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if axis_len < 1e-14 || radius < 1e-14 {
        return [0.0; 3];
    }
    let ax = [axis[0] / axis_len, axis[1] / axis_len, axis[2] / axis_len];
    let proj = dx * ax[0] + dy * ax[1] + dz * ax[2];
    let perp_x = dx - proj * ax[0];
    let perp_y = dy - proj * ax[1];
    let perp_z = dz - proj * ax[2];
    let r = (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt();
    let r_ratio = (r / radius).min(1.0);
    let speed = v_max * (1.0 - r_ratio * r_ratio);
    [speed * ax[0], speed * ax[1], speed * ax[2]]
}
/// 1/7th power law turbulent velocity profile.
///
/// `v(r) = v_max * (1 - r/R)^(1/7)` for turbulent pipe flow.
pub fn power_law_velocity_profile(
    center: [f64; 3],
    axis: [f64; 3],
    radius: f64,
    v_max: f64,
    pos: [f64; 3],
) -> [f64; 3] {
    let dx = pos[0] - center[0];
    let dy = pos[1] - center[1];
    let dz = pos[2] - center[2];
    let axis_len = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
    if axis_len < 1e-14 || radius < 1e-14 {
        return [0.0; 3];
    }
    let ax = [axis[0] / axis_len, axis[1] / axis_len, axis[2] / axis_len];
    let proj = dx * ax[0] + dy * ax[1] + dz * ax[2];
    let perp_x = dx - proj * ax[0];
    let perp_y = dy - proj * ax[1];
    let perp_z = dz - proj * ax[2];
    let r = (perp_x * perp_x + perp_y * perp_y + perp_z * perp_z).sqrt();
    let r_ratio = (r / radius).min(1.0);
    let speed = v_max * (1.0 - r_ratio).powf(1.0 / 7.0);
    [speed * ax[0], speed * ax[1], speed * ax[2]]
}
/// Extrapolate pressure from interior fluid particles to boundary/buffer particles.
///
/// Uses a linear extrapolation along the boundary normal:
/// p_boundary = p_interior + rho * g_n * d
///
/// where g_n is the gravity component along the normal and d is the distance.
pub fn extrapolate_pressure_linear(
    interior_pressure: f64,
    interior_density: f64,
    gravity: [f64; 3],
    normal: [f64; 3],
    distance: f64,
) -> f64 {
    let g_n = gravity[0] * normal[0] + gravity[1] * normal[1] + gravity[2] * normal[2];
    (interior_pressure + interior_density * g_n * distance).max(0.0)
}
/// Extrapolate pressure using hydrostatic condition for open boundaries.
///
/// Uses zeroth-order extrapolation (nearest neighbor pressure).
pub fn extrapolate_pressure_zeroth(interior_pressure: f64) -> f64 {
    interior_pressure.max(0.0)
}
/// Weighted pressure extrapolation from multiple interior particles.
///
/// `neighbors`: (pressure, distance) pairs.
/// Uses inverse-distance weighting.
pub fn extrapolate_pressure_weighted(neighbors: &[(f64, f64)]) -> f64 {
    if neighbors.is_empty() {
        return 0.0;
    }
    let mut total_weight = 0.0f64;
    let mut weighted_p = 0.0f64;
    for &(p, d) in neighbors {
        if d < 1e-14 {
            return p.max(0.0);
        }
        let w = 1.0 / d;
        weighted_p += w * p;
        total_weight += w;
    }
    if total_weight > 0.0 {
        (weighted_p / total_weight).max(0.0)
    } else {
        0.0
    }
}
/// Characteristic-based non-reflecting boundary condition.
///
/// Applies the characteristic velocity adjustment to reduce wave reflection
/// at an open boundary.
///
/// Returns the adjusted velocity for a particle near the boundary.
pub fn non_reflecting_bc_velocity(
    particle_vel: [f64; 3],
    reference_vel: [f64; 3],
    normal: [f64; 3],
    sound_speed: f64,
    density: f64,
    particle_pressure: f64,
    reference_pressure: f64,
) -> [f64; 3] {
    let v_n =
        particle_vel[0] * normal[0] + particle_vel[1] * normal[1] + particle_vel[2] * normal[2];
    let v_ref_n =
        reference_vel[0] * normal[0] + reference_vel[1] * normal[1] + reference_vel[2] * normal[2];
    let dp = particle_pressure - reference_pressure;
    let dv_n = v_n - v_ref_n;
    let correction = if density * sound_speed > 1e-14 {
        dp / (density * sound_speed) - dv_n
    } else {
        0.0
    };
    let damping = 0.5;
    [
        particle_vel[0] + damping * correction * normal[0],
        particle_vel[1] + damping * correction * normal[1],
        particle_vel[2] + damping * correction * normal[2],
    ]
}
/// Apply non-reflecting pressure correction at outflow boundary.
///
/// Returns the adjusted pressure.
pub fn non_reflecting_bc_pressure(
    particle_pressure: f64,
    reference_pressure: f64,
    relaxation: f64,
) -> f64 {
    let corrected = particle_pressure + relaxation * (reference_pressure - particle_pressure);
    corrected.max(0.0)
}
/// Blend the interior and exterior Riemann invariants to get the boundary velocity.
///
/// At an inflow: use exterior `W+` and interior `W-`.
/// Returns the normal velocity at the boundary.
#[allow(dead_code)]
pub fn riemann_boundary_velocity(
    interior_v_n: f64,
    exterior_v_n: f64,
    interior_c: f64,
    exterior_c: f64,
    gamma: f64,
) -> f64 {
    let denom = gamma - 1.0;
    if denom.abs() < 1e-12 {
        return exterior_v_n;
    }
    let w_plus = exterior_v_n + 2.0 * exterior_c / denom;
    let w_minus = interior_v_n - 2.0 * interior_c / denom;
    0.5 * (w_plus + w_minus)
}
/// Generate a turbulent velocity perturbation for a particle at `pos`.
///
/// Uses a simple random fluctuation scaled by `intensity * |v_mean|`.
/// In production code this would use a correlated random field (e.g. digital
/// filter method); here we use a deterministic hash for reproducibility.
#[allow(dead_code)]
pub fn turbulent_velocity_perturbation(
    pos: [f64; 3],
    v_mean: [f64; 3],
    config: &TurbulenceInjectionConfig,
) -> [f64; 3] {
    let mean_speed = (v_mean[0] * v_mean[0] + v_mean[1] * v_mean[1] + v_mean[2] * v_mean[2]).sqrt();
    let amp = config.intensity * mean_speed;
    let hash = |x: f64, salt: u64| -> f64 {
        let raw = (x * 1e6) as i64 as u64 ^ salt ^ config.seed;
        let r = (raw
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407)) as f64;
        (r / u64::MAX as f64) * 2.0 - 1.0
    };
    [
        amp * hash(pos[0] + pos[1], 1),
        amp * hash(pos[1] + pos[2], 2),
        amp * hash(pos[0] + pos[2], 3),
    ]
}
/// NSCBC (Navier-Stokes Characteristic Boundary Condition) relaxation for
/// inlet pressure.
///
/// Adjusts the local pressure toward the reference using the acoustic
/// relaxation: dp/dt = −σ · (p − p_ref) where σ = u_bc / L.
///
/// Returns the relaxed pressure after one time-like update.
#[allow(dead_code)]
pub fn nscbc_inlet_pressure(p_ref: f64, u_bc: f64, rho: f64, c: f64, l: f64) -> f64 {
    let p_local = rho * c * c;
    let sigma = u_bc.abs() / l.max(1e-14);
    let p_bc = p_local - sigma * (p_local - p_ref);
    p_bc.max(0.0)
}
/// Apply sponge-zone damping to a velocity vector.
///
/// The sponge is active for `x_start ≤ x ≤ x_end` and damps `vel` toward
/// `target` with coefficient `sigma` (1/s).  Returns the new velocity.
///
/// `vel_new[α] = vel[α] − σ * (vel[α] − target[α])` (explicit relaxation).
#[allow(dead_code)]
pub fn sponge_zone_damping(
    x: f64,
    x_start: f64,
    x_end: f64,
    sigma: f64,
    target: [f64; 3],
    vel: [f64; 3],
) -> [f64; 3] {
    if x < x_start || x > x_end {
        return vel;
    }
    let t = if (x_end - x_start).abs() < 1e-14 {
        1.0
    } else {
        ((x - x_start) / (x_end - x_start)).clamp(0.0, 1.0)
    };
    let s = sigma * t;
    [
        vel[0] - s * (vel[0] - target[0]),
        vel[1] - s * (vel[1] - target[1]),
        vel[2] - s * (vel[2] - target[2]),
    ]
}
/// Compute the number of particles to inject per time step at an inlet.
///
/// N_inject = round(Q_target · dt / V_particle)
///
/// where Q_target (m³/s) is the volumetric flow rate and V_particle (m³) is the
/// reference particle volume.
#[allow(dead_code)]
pub fn particles_to_inject(q_target: f64, dt: f64, v_particle: f64) -> usize {
    if v_particle < 1e-30 || q_target < 0.0 {
        return 0;
    }
    (q_target * dt / v_particle).round() as usize
}
/// Compute the target volumetric flow rate from a prescribed velocity and area.
#[allow(dead_code)]
pub fn target_flow_rate(velocity: f64, area: f64) -> f64 {
    velocity * area
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::open_boundary::types::*;
    fn make_inflow_zone() -> BoundaryZone {
        BoundaryZone::new(
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            1.0,
            BoundaryType::Inflow {
                velocity: [0.0, 0.0, 1.0],
                density: 1000.0,
            },
        )
    }
    fn make_outflow_zone() -> BoundaryZone {
        BoundaryZone::new(
            [0.0, 0.0, 10.0],
            [0.0, 0.0, -1.0],
            1.0,
            BoundaryType::Outflow,
        )
    }
    #[test]
    fn test_inflow_buffer_generate() {
        let zone = make_inflow_zone();
        let buf = InflowBuffer::new(zone, [0.0, 0.0, 1.0], 0.1);
        let pts = buf.generate_layer(-0.05, 0.05, 3, 3);
        assert_eq!(
            pts.len(),
            9,
            "Expected 9 particles from 3×3 grid, got {}",
            pts.len()
        );
    }
    #[test]
    fn test_outflow_removal() {
        let zone = make_outflow_zone();
        let damper = OutflowDamper::new(zone, 1.0, 10.0);
        let mut mgr = OpenBoundaryManager::new();
        mgr.add_outflow(damper);
        let mut positions: Vec<[f64; 3]> = vec![[0.0, 0.0, 9.0], [0.0, 0.0, 11.0]];
        let mut velocities: Vec<[f64; 3]> = vec![[0.0, 0.0, 1.0]; 2];
        let densities: Vec<f64> = vec![1000.0; 2];
        let (_injected, removed) = mgr.step(0.01, &mut positions, &mut velocities, &densities);
        assert!(
            removed.contains(&1),
            "Particle past outlet should be marked for removal, removed={removed:?}"
        );
        assert!(
            !removed.contains(&0),
            "Particle inside fluid should NOT be removed, removed={removed:?}"
        );
    }
    #[test]
    fn test_open_boundary_mass_conservation() {
        let zone = make_inflow_zone();
        let spacing = 0.1_f64;
        let v_in = [0.0, 0.0, 1.0];
        let mut buf = InflowBuffer::new(zone, v_in, spacing);
        buf.reservoir.push(BufferParticle {
            position: [0.0, 0.0, -spacing],
            velocity: v_in,
        });
        let out_zone = make_outflow_zone();
        let damper = OutflowDamper::new(out_zone, 1.0, 10.0);
        let mut mgr = OpenBoundaryManager::new();
        mgr.add_inflow(buf);
        mgr.add_outflow(damper);
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut velocities: Vec<[f64; 3]> = Vec::new();
        let mut densities: Vec<f64> = Vec::new();
        let dt = spacing;
        let mut total_injected = 0usize;
        let mut total_removed = 0usize;
        for _ in 0..10 {
            let (inj, rem) = mgr.step(dt, &mut positions, &mut velocities, &densities);
            total_injected += inj.len();
            total_removed += rem.len();
            for &idx in rem.iter().rev() {
                positions.swap_remove(idx);
                velocities.swap_remove(idx);
            }
            densities.resize(positions.len(), 1000.0);
        }
        assert!(
            total_injected >= 1,
            "Expected at least 1 injected particle over 10 steps, got {total_injected}"
        );
        assert!(
            total_removed <= total_injected,
            "Removed {total_removed} > injected {total_injected}; mass not conserved"
        );
    }
    #[test]
    fn test_damping_force() {
        let zone = BoundaryZone::new(
            [0.0, 0.0, 10.0],
            [0.0, 0.0, -1.0],
            1.0,
            BoundaryType::Outflow,
        );
        let damper = OutflowDamper::new(zone, 2.0, 5.0);
        let pos = [0.0, 0.0, 9.5];
        let vel = [0.0, 0.0, 1.0];
        let f = damper.apply_damping_force(pos, vel);
        assert!(
            f[2] < 0.0,
            "Damping force should oppose velocity; f[2]={:.4}",
            f[2]
        );
    }
    #[test]
    fn test_uniform_velocity_profile() {
        let v = uniform_velocity_profile([1.0, 2.0, 3.0], [5.0, 5.0, 5.0]);
        assert_eq!(v, [1.0, 2.0, 3.0]);
    }
    #[test]
    fn test_parabolic_velocity_at_center() {
        let center = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let v = parabolic_velocity_profile(center, axis, 1.0, 2.0, [0.0, 0.0, 0.0]);
        assert!(
            (v[2] - 2.0).abs() < 1e-10,
            "Parabolic at center should be v_max, got {}",
            v[2]
        );
    }
    #[test]
    fn test_parabolic_velocity_at_wall() {
        let center = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let v = parabolic_velocity_profile(center, axis, 1.0, 2.0, [1.0, 0.0, 0.0]);
        assert!(
            v[2].abs() < 1e-10,
            "Parabolic at wall should be 0, got {}",
            v[2]
        );
    }
    #[test]
    fn test_parabolic_velocity_between() {
        let center = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let v = parabolic_velocity_profile(center, axis, 1.0, 2.0, [0.5, 0.0, 0.0]);
        assert!((v[2] - 1.5).abs() < 1e-10, "Expected 1.5, got {}", v[2]);
    }
    #[test]
    fn test_power_law_velocity_at_center() {
        let center = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let v = power_law_velocity_profile(center, axis, 1.0, 2.0, [0.0, 0.0, 0.0]);
        assert!(
            (v[2] - 2.0).abs() < 1e-10,
            "Power law at center should be v_max, got {}",
            v[2]
        );
    }
    #[test]
    fn test_power_law_velocity_at_wall() {
        let center = [0.0, 0.0, 0.0];
        let axis = [0.0, 0.0, 1.0];
        let v = power_law_velocity_profile(center, axis, 1.0, 2.0, [1.0, 0.0, 0.0]);
        assert!(
            v[2].abs() < 1e-10,
            "Power law at wall should be 0, got {}",
            v[2]
        );
    }
    #[test]
    fn test_extrapolate_pressure_linear_zero_gravity() {
        let p = extrapolate_pressure_linear(100.0, 1000.0, [0.0; 3], [0.0, 0.0, 1.0], 0.1);
        assert!(
            (p - 100.0).abs() < 1e-10,
            "With zero gravity, extrapolated pressure should equal interior"
        );
    }
    #[test]
    fn test_extrapolate_pressure_linear_with_gravity() {
        let p =
            extrapolate_pressure_linear(100.0, 1000.0, [0.0, -9.81, 0.0], [0.0, -1.0, 0.0], 0.1);
        let expected = 100.0 + 1000.0 * 9.81 * 0.1;
        assert!((p - expected).abs() < 1e-6, "Expected {expected}, got {p}");
    }
    #[test]
    fn test_extrapolate_pressure_zeroth() {
        assert!((extrapolate_pressure_zeroth(50.0) - 50.0).abs() < 1e-12);
        assert!((extrapolate_pressure_zeroth(-10.0) - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_extrapolate_pressure_weighted() {
        let neighbors = vec![(100.0, 0.1), (200.0, 0.2)];
        let p = extrapolate_pressure_weighted(&neighbors);
        let expected = 2000.0 / 15.0;
        assert!((p - expected).abs() < 1e-8, "Expected {expected}, got {p}");
    }
    #[test]
    fn test_extrapolate_pressure_weighted_empty() {
        let p = extrapolate_pressure_weighted(&[]);
        assert!((p - 0.0).abs() < 1e-12);
    }
    #[test]
    fn test_non_reflecting_bc_no_change_at_reference() {
        let vel = [1.0, 0.0, 0.0];
        let ref_vel = [1.0, 0.0, 0.0];
        let normal = [1.0, 0.0, 0.0];
        let result = non_reflecting_bc_velocity(vel, ref_vel, normal, 340.0, 1.225, 0.0, 0.0);
        for i in 0..3 {
            assert!(
                (result[i] - vel[i]).abs() < 1e-8,
                "Velocity should not change at reference state"
            );
        }
    }
    #[test]
    fn test_non_reflecting_bc_pressure() {
        let p = non_reflecting_bc_pressure(120.0, 100.0, 0.5);
        assert!((p - 110.0).abs() < 1e-10, "Expected 110.0, got {p}");
    }
    #[test]
    fn test_buffer_zone_fill() {
        let zone = make_inflow_zone();
        let mut mgr = BufferZoneManager::new(zone, 2, 0.1);
        mgr.fill_buffer(3, 3);
        assert_eq!(mgr.count(), 18);
    }
    #[test]
    fn test_buffer_zone_prune() {
        let zone = make_inflow_zone();
        let mut mgr = BufferZoneManager::new(zone, 2, 0.1);
        mgr.fill_buffer(3, 3);
        let initial = mgr.count();
        mgr.prune_distant_particles(0.5);
        assert_eq!(mgr.count(), initial, "No particles should be pruned");
        mgr.prune_distant_particles(0.05);
        assert!(mgr.count() < initial, "Some particles should be pruned");
    }
    #[test]
    fn test_lifecycle_manager_can_inject() {
        let mgr = ParticleLifecycleManager::new(1000);
        assert!(mgr.can_inject(500, 100));
        assert!(!mgr.can_inject(950, 100));
    }
    #[test]
    fn test_lifecycle_manager_net_change() {
        let mut mgr = ParticleLifecycleManager::new(1000);
        mgr.record_injection(50);
        mgr.record_removal(20);
        assert_eq!(mgr.net_change(), 30);
    }
    #[test]
    fn test_lifecycle_manager_injection_rate() {
        let mut mgr = ParticleLifecycleManager::new(1000);
        mgr.record_injection(100);
        let rate = mgr.injection_rate(2.0);
        assert!((rate - 50.0).abs() < 1e-10);
    }
    #[test]
    fn test_boundary_zone_within_radius() {
        let zone = BoundaryZone::new([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1.0, BoundaryType::Outflow);
        assert!(zone.within_radius([0.5, 0.0, 0.0]));
        assert!(zone.within_radius([0.0, 0.5, 5.0]));
        assert!(!zone.within_radius([2.0, 0.0, 0.0]));
    }
    #[test]
    fn test_quadratic_damping_coefficient() {
        let zone = BoundaryZone::new(
            [0.0, 0.0, 10.0],
            [0.0, 0.0, -1.0],
            1.0,
            BoundaryType::Outflow,
        );
        let damper = OutflowDamper::new(zone, 2.0, 5.0);
        let alpha = damper.quadratic_damping_coefficient([0.0, 0.0, 9.0]);
        assert!((alpha - 1.25).abs() < 1e-10, "Expected 1.25, got {alpha}");
    }
    #[test]
    fn test_should_remove() {
        let zone = make_outflow_zone();
        let damper = OutflowDamper::new(zone, 1.0, 10.0);
        assert!(damper.should_remove([0.0, 0.0, 11.0]));
        assert!(!damper.should_remove([0.0, 0.0, 9.0]));
    }
    #[test]
    fn test_compute_damping_forces() {
        let zone = BoundaryZone::new(
            [0.0, 0.0, 10.0],
            [0.0, 0.0, -1.0],
            1.0,
            BoundaryType::Outflow,
        );
        let damper = OutflowDamper::new(zone, 2.0, 5.0);
        let mut mgr = OpenBoundaryManager::new();
        mgr.add_outflow(damper);
        let positions = vec![[0.0, 0.0, 9.5], [0.0, 0.0, 5.0]];
        let velocities = vec![[0.0, 0.0, 1.0], [0.0, 0.0, 1.0]];
        let forces = mgr.compute_damping_forces(&positions, &velocities);
        assert_eq!(forces.len(), 2);
        assert!(forces[0][2] < 0.0, "Should have damping force");
        assert!(forces[1][2].abs() < 1e-12, "No damping far from outlet");
    }
    #[test]
    fn test_total_buffer_count() {
        let zone = make_inflow_zone();
        let mut buf = InflowBuffer::new(zone, [0.0, 0.0, 1.0], 0.1);
        buf.reservoir.push(BufferParticle {
            position: [0.0, 0.0, -0.1],
            velocity: [0.0, 0.0, 1.0],
        });
        buf.reservoir.push(BufferParticle {
            position: [0.0, 0.0, -0.2],
            velocity: [0.0, 0.0, 1.0],
        });
        let mut mgr = OpenBoundaryManager::new();
        mgr.add_inflow(buf);
        assert_eq!(mgr.total_buffer_count(), 2);
    }
    #[test]
    fn test_riemann_buffer_particle_invariants() {
        let bp =
            RiemannBufferParticle::new([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1000.0, 1000.0, 343.0);
        let normal = [1.0, 0.0, 0.0];
        let gamma = 1.4;
        let w_plus = bp.riemann_invariant_plus(normal, gamma);
        let w_minus = bp.riemann_invariant_minus(normal, gamma);
        assert!(
            w_plus > w_minus,
            "W+ should be > W-, got W+={w_plus}, W-={w_minus}"
        );
    }
    #[test]
    fn test_riemann_boundary_velocity_subsonic_inflow() {
        let v_bc = riemann_boundary_velocity(1.0, 2.0, 343.0, 343.0, 1.4);
        assert!((v_bc - 1.5).abs() < 1e-8, "expected 1.5, got {v_bc}");
    }
    #[test]
    fn test_sponge_layer_damping_outside() {
        let zone = make_outflow_zone();
        let sponge = SpongeLayer::new(zone, 2.0, 10.0, 2);
        let alpha = sponge.damping_at([0.0, 0.0, 5.0]);
        assert!(alpha.abs() < 1e-12, "no damping far from outlet");
    }
    #[test]
    fn test_sponge_layer_damping_inside() {
        let zone = make_outflow_zone();
        let sponge = SpongeLayer::new(zone, 2.0, 10.0, 1);
        let alpha = sponge.damping_at([0.0, 0.0, 9.0]);
        assert!((alpha - 5.0).abs() < 1e-10, "expected 5.0, got {alpha}");
    }
    #[test]
    fn test_sponge_forcing_opposes_deviation() {
        let zone = make_outflow_zone();
        let sponge = SpongeLayer::new(zone, 2.0, 10.0, 1);
        let pos = [0.0, 0.0, 9.0];
        let vel = [0.0, 0.0, 2.0];
        let v_ref = [0.0, 0.0, 1.0];
        let f = sponge.apply_sponge_forcing(pos, vel, v_ref);
        assert!(
            f[2] < 0.0,
            "sponge should decelerate overspeed, got {}",
            f[2]
        );
    }
    #[test]
    fn test_turbulent_perturbation_zero_intensity() {
        let config = TurbulenceInjectionConfig {
            intensity: 0.0,
            ..Default::default()
        };
        let perturb = turbulent_velocity_perturbation([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], &config);
        for c in perturb {
            assert!(c.abs() < 1e-12, "zero intensity → no perturbation");
        }
    }
    #[test]
    fn test_turbulent_perturbation_nonzero() {
        let config = TurbulenceInjectionConfig::default();
        let v_mean = [10.0, 0.0, 0.0];
        let perturb = turbulent_velocity_perturbation([0.5, 0.5, 0.5], v_mean, &config);
        let amp = (perturb[0].powi(2) + perturb[1].powi(2) + perturb[2].powi(2)).sqrt();
        let max_amp = config.intensity * 10.0 * 3.0_f64.sqrt();
        assert!(amp <= max_amp + 1e-10, "perturbation too large: amp={amp}");
    }
    #[test]
    fn test_turbulent_inflow_boundary_apply() {
        let zone = make_inflow_zone();
        let config = TurbulenceInjectionConfig::default();
        let mut tib = TurbulentInflowBoundary::new(zone, [1.0, 0.0, 0.0], 0.1, config);
        tib.inflow.reservoir.push(BufferParticle {
            position: [0.1, 0.1, -0.1],
            velocity: [0.0, 0.0, 0.0],
        });
        tib.apply_turbulent_velocities();
        let bp = &tib.inflow.reservoir[0];
        assert!(bp.velocity[0].is_finite());
    }
    #[test]
    fn test_pressure_outlet_density_from_pressure() {
        let bc = PressureOutletBC::new(1000.0, 2.15e9, 7.0, 0.0);
        let rho = bc.target_density();
        assert!((rho - 1000.0).abs() < 1e-6, "zero pressure → rho = rho0");
    }
    #[test]
    fn test_pressure_outlet_tait_roundtrip() {
        let bc = PressureOutletBC::new(1000.0, 2.15e9, 7.0, 1e5);
        let rho = bc.target_density();
        let p = bc.pressure_from_density(rho);
        assert!(
            (p - 1e5).abs() < 1.0,
            "Tait roundtrip: p should recover target"
        );
    }
    #[test]
    fn test_pressure_outlet_corrected_pressure_at_boundary() {
        let bc = PressureOutletBC::new(1000.0, 2.15e9, 7.0, 1e5);
        let p = bc.corrected_pressure(1000.0, 0.0, 1.0);
        assert!((p - 1e5).abs() < 1.0, "at sd=0 should be p_outlet");
    }
    #[test]
    fn test_characteristic_state_invariants_subsonic() {
        let state = CharacteristicState {
            v_normal: 0.5,
            pressure: 1e5,
            density: 1.225,
            sound_speed: 340.0,
        };
        let (jp, jm) = state.acoustic_invariants();
        assert!(jp > jm, "J+ > J- for positive pressure");
        let v_rec = (jp + jm) / 2.0;
        assert!((v_rec - 0.5).abs() < 1e-10, "v reconstructed should match");
    }
    #[test]
    fn test_non_reflecting_outlet_matches_reference() {
        let state = CharacteristicState {
            v_normal: 0.5,
            pressure: 1e5,
            density: 1.225,
            sound_speed: 340.0,
        };
        let reference = CharacteristicState {
            v_normal: 0.5,
            pressure: 1e5,
            density: 1.225,
            sound_speed: 340.0,
        };
        let (v_bc, p_bc) = state.non_reflecting_outlet(&reference);
        assert!(v_bc.is_finite(), "v_bc should be finite");
        assert!(p_bc >= 0.0, "p_bc should be non-negative");
    }
    #[test]
    fn test_non_reflecting_outlet_wave_cancellation() {
        let reference = CharacteristicState {
            v_normal: 0.0,
            pressure: 1e5,
            density: 1.225,
            sound_speed: 340.0,
        };
        let perturbed = CharacteristicState {
            v_normal: 0.1,
            pressure: 1e5 + 340.0 * 1.225 * 0.1,
            density: 1.225,
            sound_speed: 340.0,
        };
        let (v_bc, p_bc) = perturbed.non_reflecting_outlet(&reference);
        assert!(v_bc.is_finite());
        assert!(p_bc >= 0.0);
    }
    #[test]
    fn test_inlet_buffer_grows_after_generate_layer() {
        let mut buf = InletBuffer::new([1.0, 0.0, 0.0], 1000.0, 3, 0.1);
        assert_eq!(buf.particles.len(), 0);
        buf.generate_layer(0.0);
        assert_eq!(
            buf.particles.len(),
            3,
            "One layer of n_layers=3 should add 3 particles"
        );
        buf.generate_layer(0.1);
        assert_eq!(
            buf.particles.len(),
            6,
            "Two layers should add 6 particles total"
        );
    }
    #[test]
    fn test_inlet_buffer_shift_moves_particles() {
        let vel = [2.0, 0.0, 0.0];
        let mut buf = InletBuffer::new(vel, 1000.0, 2, 0.1);
        buf.generate_layer(0.0);
        let x0 = buf.particles[0][0];
        buf.shift_buffer(0.5);
        let x1 = buf.particles[0][0];
        assert!(
            (x1 - (x0 + 2.0 * 0.5)).abs() < 1e-12,
            "Particle should have moved by v*dt"
        );
    }
    #[test]
    fn test_outlet_removes_particle_past_xmax() {
        let outlet = OutletZone::new(5.0, 0.5);
        let pos_past = [5.5, 0.0, 0.0];
        let mut vel = [1.0, 0.0, 0.0];
        let should_remove = outlet.apply_outlet_damping(&mut vel, pos_past);
        assert!(
            should_remove,
            "Particle past x_max should be flagged for removal"
        );
    }
    #[test]
    fn test_outlet_does_not_remove_particle_inside() {
        let outlet = OutletZone::new(5.0, 0.5);
        let pos_inside = [3.0, 0.0, 0.0];
        let mut vel = [1.0, 0.0, 0.0];
        let should_remove = outlet.apply_outlet_damping(&mut vel, pos_inside);
        assert!(
            !should_remove,
            "Particle inside domain should not be removed"
        );
        assert!(
            (vel[0] - 1.0).abs() < 1e-14,
            "Velocity should be unchanged inside domain"
        );
    }
    #[test]
    fn test_outlet_damping_reduces_velocity() {
        let outlet = OutletZone::new(5.0, 0.8);
        let pos_past = [6.0, 0.0, 0.0];
        let mut vel = [2.0, 1.0, 0.5];
        outlet.apply_outlet_damping(&mut vel, pos_past);
        assert!((vel[0] - 2.0 * 0.2).abs() < 1e-12);
        assert!((vel[1] - 1.0 * 0.2).abs() < 1e-12);
    }
    #[test]
    fn test_sponge_damps_toward_target() {
        let target = [0.0, 0.0, 0.0];
        let vel = [10.0, 0.0, 0.0];
        let result = sponge_zone_damping(5.0, 0.0, 5.0, 1.0, target, vel);
        assert!(
            result[0].abs() < 1e-12,
            "Full damping should reduce to target: {}",
            result[0]
        );
    }
    #[test]
    fn test_sponge_no_damping_outside_zone() {
        let target = [0.0, 0.0, 0.0];
        let vel = [5.0, 3.0, 1.0];
        let result = sponge_zone_damping(10.0, 0.0, 5.0, 1.0, target, vel);
        for i in 0..3 {
            assert!(
                (result[i] - vel[i]).abs() < 1e-14,
                "Velocity should be unchanged outside zone"
            );
        }
    }
    #[test]
    fn test_sponge_partial_damping_midzone() {
        let target = [0.0, 0.0, 0.0];
        let vel = [4.0, 0.0, 0.0];
        let result = sponge_zone_damping(2.5, 0.0, 5.0, 0.5, target, vel);
        let expected = 4.0 - 0.25 * 4.0;
        assert!(
            (result[0] - expected).abs() < 1e-12,
            "Expected {expected}, got {}",
            result[0]
        );
    }
    #[test]
    fn test_nscbc_pressure_finite() {
        let p = nscbc_inlet_pressure(1e5, 1.0, 1.225, 340.0, 0.1);
        assert!(p.is_finite(), "NSCBC pressure should be finite");
    }
    #[test]
    fn test_nscbc_pressure_nonnegative() {
        let p = nscbc_inlet_pressure(1e5, 1.0, 1.225, 340.0, 0.1);
        assert!(p >= 0.0, "NSCBC pressure should be non-negative, got {p}");
    }
    #[test]
    fn test_nscbc_pressure_relaxes_toward_ref() {
        let rho = 1.225_f64;
        let c = 340.0_f64;
        let p_ref = rho * c * c;
        let p = nscbc_inlet_pressure(p_ref, 1.0, rho, c, 1.0);
        assert!(
            (p - p_ref).abs() < 1.0,
            "At equilibrium NSCBC should stay near p_ref"
        );
    }
    #[test]
    fn test_flow_rate_controller_proportional() {
        let mut ctrl = FlowRateController::new(1.0, 1.0, 0.0, 0.0);
        let correction = ctrl.update(1.5, 0.01);
        assert!(
            (correction + 0.5).abs() < 1e-12,
            "P correction should be -0.5, got {correction}"
        );
    }
    #[test]
    fn test_flow_rate_controller_zero_error() {
        let mut ctrl = FlowRateController::new(1.0, 1.0, 0.0, 0.0);
        let correction = ctrl.update(1.0, 0.01);
        assert!(
            correction.abs() < 1e-12,
            "Zero error → zero correction, got {correction}"
        );
    }
    #[test]
    fn test_flow_rate_controller_reset() {
        let mut ctrl = FlowRateController::new(1.0, 1.0, 0.5, 0.0);
        ctrl.update(2.0, 0.01);
        ctrl.reset();
        assert!(ctrl.integral.abs() < 1e-14, "Integral should be reset");
        assert!(
            ctrl.prev_error.abs() < 1e-14,
            "Previous error should be reset"
        );
    }
    #[test]
    fn test_flow_rate_controller_clamp() {
        let mut ctrl = FlowRateController::new(0.0, 1000.0, 0.0, 0.0);
        ctrl.max_correction = 1.0;
        let correction = ctrl.update(1000.0, 0.01);
        assert!(
            correction.abs() <= 1.0 + 1e-12,
            "Correction should be clamped to max"
        );
    }
    #[test]
    fn test_estimate_flow_rate() {
        let velocities = vec![[1.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let normal = [1.0_f64, 0.0, 0.0];
        let masses = vec![0.001_f64; 2];
        let densities = vec![1000.0_f64; 2];
        let layer = 0.1_f64;
        let q =
            FlowRateController::estimate_flow_rate(&velocities, normal, &masses, &densities, layer);
        assert!(q > 0.0, "Flow rate must be positive for positive velocity");
        assert!((q - 2e-5).abs() < 1e-15);
    }
    #[test]
    fn test_riemann_invariants_r_plus_greater_r_minus() {
        let ri = RiemannInvariantBc::default();
        let c = ri.sound_speed(ri.rho_ref);
        let rp = ri.r_plus(10.0, c);
        let rm = ri.r_minus(10.0, c);
        assert!(rp > rm, "R+ must be > R- for c > 0");
    }
    #[test]
    fn test_riemann_inlet_velocity_recovery() {
        let ri = RiemannInvariantBc::new(1.4, 340.0, 1.225);
        let u = 10.0_f64;
        let c = 340.0_f64;
        let rp = ri.r_plus(u, c);
        let rm = ri.r_minus(u, c);
        let (u_rec, c_rec) = ri.inlet_velocity(rp, rm);
        assert!(
            (u_rec - u).abs() < 1e-12,
            "Recovered inlet velocity should match, got {u_rec}"
        );
        assert!(
            (c_rec - c).abs() < 1e-10,
            "Recovered sound speed should match, got {c_rec}"
        );
    }
    #[test]
    fn test_riemann_outlet_velocity_recovery() {
        let ri = RiemannInvariantBc::new(1.4, 340.0, 1.225);
        let u = 5.0_f64;
        let c = 340.0_f64;
        let rp = ri.r_plus(u, c);
        let rm = ri.r_minus(u, c);
        let (u_rec, c_rec) = ri.outlet_velocity(rm, rp);
        assert!(
            (u_rec - u).abs() < 1e-12,
            "Recovered outlet velocity should match, got {u_rec}"
        );
        assert!(
            (c_rec - c).abs() < 1e-10,
            "Recovered sound speed should match, got {c_rec}"
        );
    }
    #[test]
    fn test_riemann_sound_speed_reference_density() {
        let ri = RiemannInvariantBc::new(1.4, 340.0, 1.225);
        let c = ri.sound_speed(1.225);
        assert!(
            (c - 340.0).abs() < 1e-10,
            "At reference density c should equal c_ref"
        );
    }
    #[test]
    fn test_riemann_density_from_sound_speed() {
        let ri = RiemannInvariantBc::new(1.4, 340.0, 1.225);
        let rho = ri.density_from_sound_speed(340.0);
        assert!(
            (rho - 1.225).abs() < 1e-10,
            "At c_ref density should equal rho_ref"
        );
    }
    #[test]
    fn test_absorbing_outlet_damping() {
        let normal = [1.0_f64, 0.0, 0.0];
        let mut outlet = AbsorbingOutlet::new(normal, 340.0, 2, 1.0, 0.001);
        let vel = vec![[1.0_f64, 0.0, 0.0], [1.5, 0.0, 0.0]];
        let result = outlet.apply(&vel);
        assert!(
            result[0][0] < 1.0,
            "Absorbing outlet should reduce velocity toward zero initial"
        );
    }
    #[test]
    fn test_absorbing_outlet_no_history() {
        let normal = [1.0_f64, 0.0, 0.0];
        let mut outlet = AbsorbingOutlet::new(normal, 0.0, 1, 1.0, 0.01);
        let vel = vec![[5.0_f64, 3.0, 1.0]];
        let result = outlet.apply(&vel);
        for c in 0..3 {
            assert!(
                (result[0][c] - vel[0][c]).abs() < 1e-14,
                "Zero sound speed → no change"
            );
        }
    }
    #[test]
    fn test_absorbing_outlet_resize() {
        let normal = [0.0_f64, 0.0, 1.0];
        let mut outlet = AbsorbingOutlet::new(normal, 340.0, 3, 1.0, 0.01);
        outlet.resize(5);
        assert_eq!(outlet.vel_prev.len(), 5);
        assert_eq!(outlet.vel_prev2.len(), 5);
    }
    #[test]
    fn test_recycling_inlet_gamma_one_when_equal_thicknesses() {
        let ri = RecyclingTurbulenceInlet::new(1.0, 0.1, 0.1);
        assert!(
            (ri.gamma - 1.0).abs() < 1e-10,
            "γ should be 1 when δ_in == δ_re"
        );
    }
    #[test]
    fn test_recycling_inlet_generate_without_fluctuations() {
        let ri = RecyclingTurbulenceInlet::new(1.0, 0.1, 0.1);
        let mean = [1.0_f64, 0.0, 0.0];
        let v = ri.generate_inlet_velocity(0, mean);
        for c in 0..3 {
            assert!(
                (v[c] - mean[c]).abs() < 1e-14,
                "Without fluctuations, should return mean"
            );
        }
    }
    #[test]
    fn test_recycling_inlet_rms_fluctuation_zero() {
        let ri = RecyclingTurbulenceInlet::new(1.0, 0.1, 0.1);
        assert!(
            ri.rms_fluctuation().abs() < 1e-14,
            "Empty fluctuations → zero RMS"
        );
    }
    #[test]
    fn test_recycling_inlet_store_and_generate() {
        let mut ri = RecyclingTurbulenceInlet::new(1.0, 0.05, 0.1);
        let u_re = vec![[2.0_f64, 0.0, 0.0]];
        let u_mean = vec![[1.0_f64, 0.0, 0.0]];
        ri.store_fluctuations(&u_re, &u_mean);
        let mean_inlet = [1.0_f64, 0.0, 0.0];
        let v = ri.generate_inlet_velocity(0, mean_inlet);
        assert!(
            v[0] > 1.0,
            "Inlet velocity should be above mean with positive fluctuation"
        );
    }
    #[test]
    fn test_particles_to_inject_round_trip() {
        let n = particles_to_inject(0.01, 0.1, 0.001);
        assert_eq!(n, 1, "Should inject 1 particle");
    }
    #[test]
    fn test_particles_to_inject_zero_volume() {
        let n = particles_to_inject(1.0, 0.01, 0.0);
        assert_eq!(n, 0, "Zero particle volume → 0 injections");
    }
    #[test]
    fn test_target_flow_rate() {
        let q = target_flow_rate(2.0, 0.5);
        assert!((q - 1.0).abs() < 1e-14, "Q = v*A = 2*0.5 = 1, got {q}");
    }
    #[test]
    fn test_lodi_l1_inlet_zero_at_reference_pressure() {
        let lodi = LodiBC::default();
        let l1 = lodi.l1_inlet(101325.0, 0.1);
        assert!(
            l1.abs() < 1e-6,
            "L1 should be zero at reference pressure, got {l1}"
        );
    }
    #[test]
    fn test_lodi_pressure_update_symmetric() {
        let lodi = LodiBC::default();
        let l = 100.0_f64;
        let dp = lodi.pressure_update(l, l, 0.01);
        assert!(
            (dp + l * 0.01).abs() < 1e-12,
            "Pressure update should be -L1*dt for equal waves"
        );
    }
    #[test]
    fn test_lodi_velocity_update_zero_for_equal_waves() {
        let lodi = LodiBC::default();
        let du = lodi.velocity_update(100.0, 100.0, 0.01);
        assert!(
            du.abs() < 1e-14,
            "Velocity update should be zero for L1==L5"
        );
    }
    #[test]
    fn test_lodi_density_update_finite() {
        let lodi = LodiBC::default();
        let drho = lodi.density_update(50.0, 50.0, 0.01);
        assert!(drho.is_finite(), "Density update should be finite");
    }
    #[test]
    fn test_convective_outflow_courant_zero() {
        let mut bc = ConvectiveOutflowBC::new(0.0, 1.0, 0.01, 2);
        let vel_out = vec![[3.0_f64, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let vel_int = vec![[1.0_f64, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let result = bc.apply(&vel_out, &vel_int);
        for i in 0..2 {
            assert!(
                (result[i][0] - vel_out[i][0]).abs() < 1e-14,
                "No advection for c=0"
            );
        }
    }
    #[test]
    fn test_convective_outflow_reduces_velocity() {
        let mut bc = ConvectiveOutflowBC::new(1.0, 1.0, 0.5, 1);
        let vel_out = vec![[4.0_f64, 0.0, 0.0]];
        let vel_int = vec![[2.0_f64, 0.0, 0.0]];
        let result = bc.apply(&vel_out, &vel_int);
        assert!(
            (result[0][0] - 3.0).abs() < 1e-14,
            "Expected 3.0, got {}",
            result[0][0]
        );
    }
    #[test]
    fn test_convective_outflow_resize() {
        let mut bc = ConvectiveOutflowBC::new(1.0, 1.0, 0.01, 2);
        bc.resize(5);
        assert_eq!(bc.vel_outlet_prev.len(), 5);
    }
}
