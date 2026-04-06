//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use oxiphysics_core::math::Vec3;

#[allow(unused_imports)]
use crate::particle::types::*;

/// Public re-export of the cubic-spline gradient for use in other modules.
#[doc(hidden)]
pub fn cubic_spline_gradient_pub(r: f64, h: f64) -> f64 {
    cubic_spline_gradient(r, h)
}
/// Cubic-spline gradient magnitude: dW/dr (unnormalized, 3-D).
pub fn cubic_spline_gradient(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q > 2.0 {
        return 0.0;
    }
    let sigma = 3.0 / (2.0 * std::f64::consts::PI * h.powi(4));
    if q <= 1.0 {
        sigma * (-3.0 * q + 2.25 * q * q)
    } else {
        sigma * (-0.75 * (2.0 - q).powi(2))
    }
}
/// Cubic-spline kernel value W(r, h).
pub fn cubic_spline_kernel(r: f64, h: f64) -> f64 {
    let q = r / h;
    let sigma = 1.0 / (std::f64::consts::PI * h * h * h);
    if q <= 1.0 {
        sigma * (1.0 - 1.5 * q * q + 0.75 * q * q * q)
    } else if q <= 2.0 {
        sigma * 0.25 * (2.0 - q).powi(3)
    } else {
        0.0
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    fn make_set(n: usize) -> SphParticleSet {
        let mut ps = SphParticleSet::new();
        for i in 0..n {
            let x = i as f64;
            ps.add_raw([x, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 0.001, false);
        }
        ps
    }
    #[test]
    fn sph_particle_set_add_and_len() {
        let mut ps = SphParticleSet::new();
        assert!(ps.is_empty());
        let p = SphParticle::new(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0);
        ps.add(&p);
        assert_eq!(ps.len(), 1);
    }
    #[test]
    fn sph_particle_set_remove_swap() {
        let mut ps = make_set(3);
        ps.remove(0);
        assert_eq!(ps.len(), 2);
        assert!((ps.positions[0][0] - 2.0).abs() < 1e-12);
    }
    #[test]
    fn clear_accelerations_zeroes_all() {
        let mut ps = make_set(4);
        for a in &mut ps.accelerations {
            *a = [9.0, 8.0, 7.0];
        }
        ps.clear_accelerations();
        for a in &ps.accelerations {
            assert_eq!(*a, [0.0, 0.0, 0.0]);
        }
    }
    #[test]
    fn kinetic_energy_correctness() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [3.0, 4.0, 0.0], 2.0, 0.0, false);
        let ke = ps.kinetic_energy();
        assert!((ke - 25.0).abs() < 1e-10, "ke={ke}");
    }
    #[test]
    fn center_of_mass_two_equal_masses() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.add_raw([2.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, false);
        let com = ps.center_of_mass();
        assert!((com[0] - 1.0).abs() < 1e-12);
        assert!(com[1].abs() < 1e-12);
    }
    #[test]
    fn apply_gravity_skips_boundary() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, true);
        ps.apply_gravity([0.0, -9.81, 0.0]);
        assert!((ps.accelerations[0][1] - (-9.81)).abs() < 1e-12);
        assert_eq!(ps.accelerations[1][1], 0.0, "boundary must not get gravity");
    }
    #[test]
    fn integrate_euler_moves_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [0.0, -9.81, 0.0];
        ps.integrate(0.1);
        assert!((ps.velocities[0][0] - 1.0).abs() < 1e-12);
        assert!((ps.velocities[0][1] - (-0.981)).abs() < 1e-10);
        assert!((ps.positions[0][0] - 0.1).abs() < 1e-12);
    }
    #[test]
    fn apply_pbc_wraps_coordinates() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.05, -0.1, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.apply_pbc([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(
            (ps.positions[0][0] - 0.05).abs() < 1e-12,
            "x={}",
            ps.positions[0][0]
        );
        assert!(
            (ps.positions[0][1] - 0.9).abs() < 1e-12,
            "y={}",
            ps.positions[0][1]
        );
    }
    #[test]
    fn total_momentum_calculation() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 3.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0, false);
        let mom = ps.total_momentum();
        assert!((mom[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn add_and_remove_particle() {
        let mut ps = ParticleSet::new();
        let p = SphParticle::new(Vec3::new(1.0, 2.0, 3.0), Vec3::zeros(), 1.0);
        ps.add_particle(&p);
        assert_eq!(ps.len(), 1);
        let p2 = SphParticle::new(Vec3::new(4.0, 5.0, 6.0), Vec3::zeros(), 2.0);
        ps.add_particle(&p2);
        assert_eq!(ps.len(), 2);
        ps.remove(0);
        assert_eq!(ps.len(), 1);
        assert!((ps.positions[0].x - 4.0).abs() < 1e-10);
    }
    #[test]
    fn particle_set_is_empty() {
        let ps = ParticleSet::new();
        assert!(ps.is_empty());
    }
    #[test]
    fn integrate_euler_alias_matches_integrate() {
        let mut ps1 = make_set(3);
        let mut ps2 = ps1.clone();
        ps1.accelerations[0] = [1.0, -2.0, 3.0];
        ps2.accelerations[0] = [1.0, -2.0, 3.0];
        ps1.integrate(0.01);
        ps2.integrate_euler(0.01);
        for k in 0..3 {
            assert!((ps1.positions[0][k] - ps2.positions[0][k]).abs() < 1e-14);
            assert!((ps1.velocities[0][k] - ps2.velocities[0][k]).abs() < 1e-14);
        }
    }
    #[test]
    fn verlet_half_position_update() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [0.0, 2.0, 0.0];
        let dt = 0.1;
        ps.integrate_verlet_half(dt);
        assert!((ps.positions[0][0] - 0.1).abs() < 1e-12);
        assert!((ps.positions[0][1] - 0.01).abs() < 1e-12);
        assert!((ps.velocities[0][0] - 1.0).abs() < 1e-12);
        assert!((ps.velocities[0][1] - 0.1).abs() < 1e-12);
    }
    #[test]
    fn potential_energy_positive_height() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 10.0, 0.0], [0.0, 0.0, 0.0], 2.0, 0.0, false);
        let g = [0.0, -9.81, 0.0];
        let pe = ps.potential_energy(g);
        assert!((pe - 196.2).abs() < 1e-10);
    }
    #[test]
    fn angular_momentum_circular() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0, 0.0, false);
        let l = ps.angular_momentum();
        assert!((l[2] - 1.0).abs() < 1e-12);
        assert!(l[0].abs() < 1e-12);
        assert!(l[1].abs() < 1e-12);
    }
    #[test]
    fn bounding_box_correct() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([-1.0, 2.0, 0.5], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([3.0, -1.0, 4.0], [0.0; 3], 1.0, 0.0, false);
        let (lo, hi) = ps.bounding_box();
        assert!((lo[0] - (-1.0)).abs() < 1e-12);
        assert!((lo[1] - (-1.0)).abs() < 1e-12);
        assert!((lo[2] - 0.5).abs() < 1e-12);
        assert!((hi[0] - 3.0).abs() < 1e-12);
        assert!((hi[1] - 2.0).abs() < 1e-12);
        assert!((hi[2] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn scale_velocities_halves_speed() {
        let mut ps = make_set(3);
        let orig_speed = ps.max_speed();
        ps.scale_velocities(0.5);
        let new_speed = ps.max_speed();
        assert!((new_speed - 0.5 * orig_speed).abs() < 1e-12);
    }
    #[test]
    fn fluid_and_boundary_counts() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, true);
        ps.add_raw([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        assert_eq!(ps.fluid_count(), 2);
        assert_eq!(ps.boundary_count(), 1);
        assert_eq!(ps.total_mass(), 3.0);
    }
    #[test]
    fn clamp_velocities_enforces_limit() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [3.0, 4.0, 0.0], 1.0, 0.0, false);
        ps.clamp_velocities(2.0);
        let speed = ps.max_speed();
        assert!((speed - 2.0).abs() < 1e-10);
    }
    #[test]
    fn reflect_wall_bounces_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, -0.1, 0.0], [0.0, -1.0, 0.0], 1.0, 0.0, false);
        ps.reflect_wall(1, 0.0, 0.8);
        assert!((ps.positions[0][1] - 0.1).abs() < 1e-12);
        assert!((ps.velocities[0][1] - 0.8).abs() < 1e-12);
    }
    #[test]
    fn leapfrog_kick_drift_moves_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [2.0, 0.0, 0.0];
        let dt = 0.1;
        ps.leapfrog_kick_drift(dt);
        assert!((ps.velocities[0][0] - 1.1).abs() < 1e-12);
        assert!((ps.positions[0][0] - 0.11).abs() < 1e-12);
    }
    #[test]
    fn sph_particle_kinetic_energy() {
        let p = SphParticle::new(Vec3::zeros(), Vec3::new(3.0, 4.0, 0.0), 2.0);
        assert!((p.kinetic_energy() - 25.0).abs() < 1e-12);
    }
    #[test]
    fn average_density_calculation() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1000.0;
        ps.densities[1] = 1200.0;
        assert!((ps.average_density() - 1100.0).abs() < 1e-10);
        assert!((ps.max_density() - 1200.0).abs() < 1e-10);
    }
    #[test]
    fn mean_neighbor_distance_k1() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let dists = ps.mean_neighbor_distance(1);
        assert!((dists[0] - 1.0).abs() < 1e-12);
        assert!((dists[1] - 1.0).abs() < 1e-12);
    }
    #[test]
    fn mean_neighbor_distance_boundary_zero() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, true);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let dists = ps.mean_neighbor_distance(1);
        assert_eq!(
            dists[0], 0.0,
            "Boundary particle should have 0 mean neighbor distance"
        );
    }
    #[test]
    fn neighbor_counts_radius() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.5, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let counts = ps.neighbor_counts(1.0);
        assert_eq!(
            counts[0], 1,
            "Particle 0 should have 1 neighbor within radius 1.0"
        );
        assert_eq!(
            counts[2], 0,
            "Particle 2 should have 0 neighbors within radius 1.0"
        );
    }
    #[test]
    fn density_variance_uniform() {
        let mut ps = SphParticleSet::new();
        for _ in 0..4 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
            if let Some(d) = ps.densities.last_mut() {
                *d = 1000.0;
            }
        }
        assert!(
            ps.density_variance().abs() < 1e-12,
            "Uniform densities should have zero variance"
        );
    }
    #[test]
    fn density_std_is_sqrt_variance() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 900.0;
        ps.densities[1] = 1100.0;
        let std = ps.density_std();
        let var = ps.density_variance();
        assert!((std - var.sqrt()).abs() < 1e-12);
    }
    #[test]
    fn velocity_divergence_no_panic() {
        let mut ps = SphParticleSet::new();
        for i in 0..5 {
            ps.add_raw([i as f64 * 0.1, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, false);
            *ps.densities.last_mut().unwrap() = 1000.0;
        }
        let divv = ps.velocity_divergence(0.5);
        assert_eq!(divv.len(), ps.len());
    }
    #[test]
    fn vorticity_magnitude_no_panic() {
        let mut ps = SphParticleSet::new();
        for i in 0..4 {
            ps.add_raw([i as f64 * 0.1, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0, 0.0, false);
            *ps.densities.last_mut().unwrap() = 1000.0;
        }
        let omega = ps.vorticity_magnitude(0.5);
        assert_eq!(omega.len(), ps.len());
    }
    #[test]
    fn export_positions_csv_header() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 2.0, 3.0], [0.0; 3], 1.0, 0.0, false);
        let csv = ps.export_positions_csv();
        assert!(csv.starts_with("x,y,z"), "CSV should start with header");
        assert!(csv.contains("1"), "CSV should contain position data");
    }
    #[test]
    fn export_full_csv_has_columns() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 2.0, 0.001, false);
        let csv = ps.export_full_csv();
        assert!(
            csv.contains("density"),
            "Full CSV should have density column"
        );
        assert!(
            csv.contains("pressure"),
            "Full CSV should have pressure column"
        );
    }
    #[test]
    fn export_import_positions_flat_round_trip() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 2.0, 3.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([4.0, 5.0, 6.0], [0.0; 3], 1.0, 0.0, false);
        let flat = ps.export_positions_flat();
        let mut ps2 = ps.clone();
        ps2.positions[0] = [0.0; 3];
        ps2.import_positions_flat(&flat);
        assert!((ps2.positions[0][0] - 1.0).abs() < 1e-15);
        assert!((ps2.positions[1][2] - 6.0).abs() < 1e-15);
    }
    #[test]
    fn compute_sph_densities_positive() {
        let mut ps = SphParticleSet::new();
        for i in 0..5 {
            ps.add_raw([i as f64 * 0.05, 0.0, 0.0], [0.0; 3], 0.01, 0.0, false);
        }
        ps.compute_sph_densities(0.2);
        for &d in &ps.densities {
            assert!(d > 0.0, "SPH density should be positive, got {d}");
        }
    }
    #[test]
    fn merge_two_particles() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 2.0, 0.0, false);
        ps.add_raw([2.0, 0.0, 0.0], [3.0, 0.0, 0.0], 2.0, 0.0, false);
        let surviving = ps.merge(0, 1);
        assert_eq!(ps.len(), 1);
        assert!((ps.positions[surviving][0] - 1.0).abs() < 1e-12);
        assert!((ps.velocities[surviving][0] - 2.0).abs() < 1e-12);
        assert!((ps.masses[surviving] - 4.0).abs() < 1e-12);
    }
    #[test]
    fn split_particle_doubles_count() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 0.0, 0.0], [0.5, 0.0, 0.0], 2.0, 0.0, false);
        ps.densities[0] = 1000.0;
        ps.split(0, 0.05, 0);
        assert_eq!(ps.len(), 2);
        assert!((ps.masses[0] - 1.0).abs() < 1e-12);
        assert!((ps.masses[1] - 1.0).abs() < 1e-12);
        assert!((ps.positions[0][0] - 0.95).abs() < 1e-12);
        assert!((ps.positions[1][0] - 1.05).abs() < 1e-12);
    }
    #[test]
    fn tracker_records_trajectory() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        let mut tracker = ParticleTracker::new(0);
        tracker.record(&ps, 0.0);
        ps.positions[0][0] = 1.0;
        tracker.record(&ps, 1.0);
        ps.positions[0][0] = 2.0;
        tracker.record(&ps, 2.0);
        assert_eq!(tracker.trajectory.len(), 3);
        assert!((tracker.path_length() - 2.0).abs() < 1e-12);
        assert!((tracker.mean_speed() - 1.0).abs() < 1e-12);
    }
    #[test]
    fn tracker_reset_clears_data() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        let mut tracker = ParticleTracker::new(0);
        tracker.record(&ps, 0.0);
        tracker.record(&ps, 1.0);
        assert!(!tracker.trajectory.is_empty());
        tracker.reset();
        assert!(tracker.trajectory.is_empty());
        assert!(tracker.times.is_empty());
    }
    #[test]
    fn tracker_path_length_single_point() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        let mut tracker = ParticleTracker::new(0);
        tracker.record(&ps, 0.0);
        assert!((tracker.path_length() - 0.0).abs() < 1e-12);
    }
    #[test]
    fn speed_variance_uniform_velocity() {
        let mut ps = SphParticleSet::new();
        for _ in 0..4 {
            ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        }
        let var = ps.speed_variance();
        assert!(
            var.abs() < 1e-12,
            "Uniform speed should have zero variance, got {var}"
        );
    }
    #[test]
    fn with_capacity_creates_correct_initial_state() {
        let ps = SphParticleSet::with_capacity(100);
        assert!(ps.is_empty());
        assert_eq!(ps.len(), 0);
    }
    #[test]
    fn add_raw_then_remove_all_leaves_empty() {
        let mut ps = SphParticleSet::new();
        for i in 0..5 {
            ps.add_raw([i as f64, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        }
        while !ps.is_empty() {
            ps.remove(0);
        }
        assert!(ps.is_empty());
    }
    #[test]
    fn apply_pbc_wraps_above_box() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([2.5, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.apply_pbc([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(
            (ps.positions[0][0] - 0.5).abs() < 1e-12,
            "x after PBC={}",
            ps.positions[0][0]
        );
    }
    #[test]
    fn apply_pbc_wraps_far_below_box() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([-3.5, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.apply_pbc([0.0, 0.0, 0.0], [2.0, 2.0, 2.0]);
        assert!(
            (ps.positions[0][0] - 0.5).abs() < 1e-12,
            "x after PBC={}",
            ps.positions[0][0]
        );
    }
    #[test]
    fn apply_pbc_zero_length_dimension_is_noop() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([5.0, 3.0, -1.0], [0.0; 3], 1.0, 0.0, false);
        ps.apply_pbc([0.0, 1.0, 0.0], [0.0, 1.0, 0.0]);
        assert!((ps.positions[0][0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn velocity_verlet_half_boundary_unchanged() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [2.0, 0.0, 0.0], 1.0, 0.0, true);
        ps.accelerations[0] = [10.0, 0.0, 0.0];
        ps.integrate_verlet_half(0.1);
        assert_eq!(ps.positions[0], [0.0, 0.0, 0.0]);
        assert_eq!(ps.velocities[0], [2.0, 0.0, 0.0]);
    }
    #[test]
    fn finish_verlet_completes_velocity_update() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [2.0, 0.0, 0.0];
        let dt = 0.1;
        ps.finish_verlet(dt);
        assert!((ps.velocities[0][0] - 1.1).abs() < 1e-12);
    }
    #[test]
    fn velocity_verlet_round_trip_energy() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 10.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [0.0, -9.81, 0.0];
        let g = [0.0, -9.81, 0.0];
        let e0 = ps.total_energy(g);
        let dt = 0.001;
        ps.integrate_verlet_half(dt);
        ps.accelerations[0] = [0.0, -9.81, 0.0];
        ps.finish_verlet(dt);
        let e1 = ps.total_energy(g);
        assert!(
            (e1 - e0).abs() < 0.01,
            "Energy drift too large: delta={}",
            (e1 - e0).abs()
        );
    }
    #[test]
    fn clamp_velocities_below_limit_unchanged() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.clamp_velocities(5.0);
        assert!((ps.velocities[0][0] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn clamp_velocities_3d_direction_preserved() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [3.0, 4.0, 0.0], 1.0, 0.0, false);
        ps.clamp_velocities(2.5);
        let v = &ps.velocities[0];
        let speed = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        assert!((speed - 2.5).abs() < 1e-10);
        assert!((v[0] / v[1] - 3.0 / 4.0).abs() < 1e-10);
    }
    #[test]
    fn reflect_wall_restitution_one_perfect_elastic() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, -0.01, 0.0], [0.0, -5.0, 0.0], 1.0, 0.0, false);
        ps.reflect_wall(1, 0.0, 1.0);
        assert!((ps.velocities[0][1] - 5.0).abs() < 1e-12);
        assert!((ps.positions[0][1] - 0.01).abs() < 1e-12);
    }
    #[test]
    fn reflect_wall_restitution_zero_stops_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, -0.05, 0.0], [0.0, -3.0, 0.0], 1.0, 0.0, false);
        ps.reflect_wall(1, 0.0, 0.0);
        assert!(ps.velocities[0][1].abs() < 1e-14);
    }
    #[test]
    fn reflect_wall_skips_boundary_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, -0.1, 0.0], [0.0, -1.0, 0.0], 1.0, 0.0, true);
        ps.reflect_wall(1, 0.0, 1.0);
        assert!((ps.positions[0][1] - (-0.1)).abs() < 1e-14);
        assert!((ps.velocities[0][1] - (-1.0)).abs() < 1e-14);
    }
    #[test]
    fn reflect_wall_above_wall_unchanged() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 1.0, 0.0], [0.0, -1.0, 0.0], 1.0, 0.0, false);
        ps.reflect_wall(1, 0.0, 1.0);
        assert!((ps.positions[0][1] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn total_energy_equals_ke_plus_pe() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 2.0, 0.0], [3.0, 0.0, 0.0], 1.0, 0.0, false);
        let g = [0.0, -9.81, 0.0];
        let ke = ps.kinetic_energy();
        let pe = ps.potential_energy(g);
        let total = ps.total_energy(g);
        assert!((total - (ke + pe)).abs() < 1e-12);
    }
    #[test]
    fn angular_momentum_conserved_no_force() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0, 0.0, false);
        let l0 = ps.angular_momentum();
        ps.clear_accelerations();
        ps.integrate(0.01);
        let l1 = ps.angular_momentum();
        assert!(
            (l1[2] - l0[2]).abs() < 0.1,
            "L_z drifted: {}",
            (l1[2] - l0[2]).abs()
        );
    }
    #[test]
    fn bounding_box_single_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([3.0, 7.0, -2.0], [0.0; 3], 1.0, 0.0, false);
        let (lo, hi) = ps.bounding_box();
        assert_eq!(lo, [3.0, 7.0, -2.0]);
        assert_eq!(hi, [3.0, 7.0, -2.0]);
    }
    #[test]
    fn bounding_box_empty_returns_zero() {
        let ps = SphParticleSet::new();
        let (lo, hi) = ps.bounding_box();
        assert_eq!(lo, [0.0, 0.0, 0.0]);
        assert_eq!(hi, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn scale_velocities_to_zero() {
        let mut ps = make_set(5);
        ps.scale_velocities(0.0);
        for v in &ps.velocities {
            assert_eq!(*v, [0.0, 0.0, 0.0]);
        }
    }
    #[test]
    fn total_mass_sums_correctly() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.5, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 2.5, 0.0, false);
        assert!((ps.total_mass() - 4.0).abs() < 1e-14);
    }
    #[test]
    fn fluid_count_boundary_count_sum_to_len() {
        let mut ps = SphParticleSet::new();
        for i in 0..6 {
            ps.add_raw([i as f64, 0.0, 0.0], [0.0; 3], 1.0, 0.0, i % 2 == 0);
        }
        assert_eq!(ps.fluid_count() + ps.boundary_count(), ps.len());
    }
    #[test]
    fn sph_particle_new_boundary_flag() {
        let p = SphParticle::new_boundary(Vec3::new(1.0, 2.0, 3.0), 0.5);
        assert!(p.is_boundary);
        assert_eq!(p.mass, 0.5);
    }
    #[test]
    fn sph_particle_new_not_boundary() {
        let p = SphParticle::new(Vec3::zeros(), Vec3::new(1.0, 0.0, 0.0), 1.0);
        assert!(!p.is_boundary);
    }
    #[test]
    fn center_of_mass_weighted_three_particles() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([3.0, 0.0, 0.0], [0.0; 3], 2.0, 0.0, false);
        let com = ps.center_of_mass();
        assert!((com[0] - 2.0).abs() < 1e-14);
    }
    #[test]
    fn center_of_mass_empty_returns_zero() {
        let ps = SphParticleSet::new();
        let com = ps.center_of_mass();
        assert_eq!(com, [0.0, 0.0, 0.0]);
    }
    #[test]
    fn max_speed_empty_is_zero() {
        let ps = SphParticleSet::new();
        assert_eq!(ps.max_speed(), 0.0);
    }
    #[test]
    fn max_speed_single_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0, 3.0, 4.0], 1.0, 0.0, false);
        assert!((ps.max_speed() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn integrate_euler_boundary_not_moved() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([5.0, 5.0, 5.0], [1.0, 1.0, 1.0], 1.0, 0.0, true);
        ps.accelerations[0] = [2.0, 2.0, 2.0];
        ps.integrate(0.1);
        assert_eq!(ps.positions[0], [5.0, 5.0, 5.0]);
        assert_eq!(ps.velocities[0], [1.0, 1.0, 1.0]);
    }
    #[test]
    fn leapfrog_kick_drift_boundary_unchanged() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.0, true);
        ps.accelerations[0] = [5.0, 0.0, 0.0];
        ps.leapfrog_kick_drift(0.1);
        assert_eq!(ps.positions[0], [0.0, 0.0, 0.0]);
    }
    #[test]
    fn average_density_empty_is_zero() {
        let ps = SphParticleSet::new();
        assert_eq!(ps.average_density(), 0.0);
    }
    #[test]
    fn density_variance_empty_is_zero() {
        let ps = SphParticleSet::new();
        assert_eq!(ps.density_variance(), 0.0);
    }
    #[test]
    fn split_preserves_total_mass() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 5.0, 0.0, false);
        let mass_before = ps.total_mass();
        ps.split(0, 0.1, 2);
        assert!(
            (ps.total_mass() - mass_before).abs() < 1e-12,
            "Mass not conserved after split"
        );
    }
    #[test]
    fn merge_mass_conserved() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 3.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 7.0, 0.0, false);
        let mass_before = ps.total_mass();
        ps.merge(0, 1);
        assert!((ps.total_mass() - mass_before).abs() < 1e-12);
    }
    #[test]
    fn neighbor_counts_all_within_radius() {
        let mut ps = SphParticleSet::new();
        for i in 0..5 {
            ps.add_raw([i as f64 * 0.1, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        }
        let counts = ps.neighbor_counts(0.5);
        assert_eq!(counts[0], 4, "particle 0 should see all 4 others");
    }
    #[test]
    fn export_positions_flat_round_trip_length() {
        let mut ps = SphParticleSet::new();
        for i in 0..3 {
            ps.add_raw([i as f64; 3], [0.0; 3], 1.0, 0.0, false);
        }
        let flat = ps.export_positions_flat();
        assert_eq!(flat.len(), 9, "3 particles × 3 coords = 9 values");
    }
}
/// Encode a 3-D grid coordinate `(ix, iy, iz)` into a 21-bit-per-axis Morton
/// (Z-order) code.  Coordinates must fit in 21 bits each (< 2^21).
pub fn morton_encode(ix: u32, iy: u32, iz: u32) -> u64 {
    /// Spread the low 21 bits of `x` into every 3rd bit.
    fn expand(mut x: u32) -> u64 {
        let mut x64 = x as u64 & 0x001f_ffff;
        x64 = (x64 | (x64 << 32)) & 0x001f_0000_0000_ffff;
        x64 = (x64 | (x64 << 16)) & 0x001f_0000_ff00_00ff;
        x64 = (x64 | (x64 << 8)) & 0x100f_00f0_0f00_f00f;
        x64 = (x64 | (x64 << 4)) & 0x10c3_0c30_c30c_30c3;
        x64 = (x64 | (x64 << 2)) & 0x1249_2492_4924_9249;
        x = 0;
        let _ = x;
        x64
    }
    expand(ix) | (expand(iy) << 1) | (expand(iz) << 2)
}
#[cfg(test)]
mod tests_particle_ext {
    use super::*;
    fn make_cluster(n: usize, spacing: f64) -> SphParticleSet {
        let mut ps = SphParticleSet::new();
        let side = (n as f64).cbrt().ceil() as usize;
        let mut count = 0;
        'outer: for i in 0..side {
            for j in 0..side {
                for k in 0..side {
                    if count >= n {
                        break 'outer;
                    }
                    ps.add_raw(
                        [i as f64 * spacing, j as f64 * spacing, k as f64 * spacing],
                        [0.0; 3],
                        0.001,
                        0.0,
                        false,
                    );
                    count += 1;
                }
            }
        }
        ps
    }
    #[test]
    fn morton_encode_origin_is_zero() {
        assert_eq!(morton_encode(0, 0, 0), 0);
    }
    #[test]
    fn morton_encode_x_axis() {
        let c = morton_encode(1, 0, 0);
        assert_eq!(c, 1);
    }
    #[test]
    fn morton_encode_y_axis() {
        let c = morton_encode(0, 1, 0);
        assert_eq!(c, 2);
    }
    #[test]
    fn morton_encode_z_axis() {
        let c = morton_encode(0, 0, 1);
        assert_eq!(c, 4);
    }
    #[test]
    fn sort_by_morton_preserves_particle_count() {
        let mut ps = make_cluster(8, 0.1);
        let lo = [0.0; 3];
        let hi = [1.0; 3];
        ps.sort_by_morton(lo, hi);
        assert_eq!(ps.len(), 8);
    }
    #[test]
    fn sort_by_morton_all_positions_present() {
        let mut ps = make_cluster(8, 0.1);
        let orig_positions: Vec<[f64; 3]> = ps.positions.clone();
        let lo = [0.0; 3];
        let hi = [1.0; 3];
        ps.sort_by_morton(lo, hi);
        for orig in &orig_positions {
            assert!(
                ps.positions.iter().any(|p| {
                    (p[0] - orig[0]).abs() < 1e-14
                        && (p[1] - orig[1]).abs() < 1e-14
                        && (p[2] - orig[2]).abs() < 1e-14
                }),
                "position {orig:?} missing after Morton sort"
            );
        }
    }
    #[test]
    fn iter_fluid_positions_skips_boundary() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, true);
        ps.add_raw([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let fluid: Vec<_> = ps.iter_fluid_positions().collect();
        assert_eq!(fluid.len(), 2, "should yield 2 fluid particles");
        for (_, p) in &fluid {
            assert!(((*p)[0] - 1.0).abs() > 1e-12, "boundary must not appear");
        }
    }
    #[test]
    fn adaptive_smoothing_lengths_boundary_uses_h0() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, true);
        ps.densities[0] = 1000.0;
        let hs = ps.adaptive_smoothing_lengths(1.3, 0.05);
        assert!(
            (hs[0] - 0.05).abs() < 1e-14,
            "boundary particle should use h0"
        );
    }
    #[test]
    fn adaptive_smoothing_lengths_fluid_positive() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1000.0;
        let hs = ps.adaptive_smoothing_lengths(1.3, 0.05);
        assert!(hs[0] > 0.0, "fluid smoothing length must be positive");
    }
    #[test]
    fn leapfrog_full_constant_gravity() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0, 1.0, 0.0], [0.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [0.0, -9.81, 0.0];
        ps.leapfrog_full(0.1);
        assert!(
            ps.positions[0][1] < 1.0,
            "particle should fall under gravity"
        );
        assert!(ps.positions[0][0].abs() < 1e-14);
    }
    #[test]
    fn verlet_full_constant_acceleration() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.accelerations[0] = [2.0, 0.0, 0.0];
        let dt = 0.1;
        ps.verlet_full(dt);
        assert!((ps.positions[0][0] - 0.11).abs() < 1e-12);
        assert!((ps.velocities[0][0] - 1.2).abs() < 1e-12);
    }
    #[test]
    fn radius_of_gyration_two_symmetric_particles() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([-1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let rg = ps.radius_of_gyration();
        assert!((rg - 1.0).abs() < 1e-12, "R_g should be 1, got {rg}");
    }
    #[test]
    fn radius_of_gyration_single_particle_at_origin() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        let rg = ps.radius_of_gyration();
        assert!(rg.abs() < 1e-14, "single particle at origin → R_g = 0");
    }
    #[test]
    fn granular_temperature_stationary_zero() {
        let mut ps = SphParticleSet::new();
        for _ in 0..4 {
            ps.add_raw([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.0, false);
        }
        let t = ps.granular_temperature();
        assert!(t.abs() < 1e-12, "T_g = 0 for uniform velocity, got {t}");
    }
    #[test]
    fn granular_temperature_opposite_velocities() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 1.0, 0.0, false);
        let t = ps.granular_temperature();
        assert!((t - 1.0).abs() < 1e-12, "T_g should be 1, got {t}");
    }
    #[test]
    fn add_particle_h_returns_correct_index() {
        let mut ps = SphParticleSet::new();
        let idx0 = ps.add_particle_h([0.0; 3], [0.0; 3], 1.0, 0.1);
        let idx1 = ps.add_particle_h([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        assert_eq!(idx0, 0);
        assert_eq!(idx1, 1);
        assert_eq!(ps.len(), 2);
    }
    #[test]
    fn remove_particle_alias_works() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.remove_particle(0);
        assert_eq!(ps.len(), 1);
    }
    #[test]
    fn iter_indices_covers_all() {
        let mut ps = SphParticleSet::new();
        for _ in 0..5 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        }
        let indices: Vec<usize> = ps.iter_indices().collect();
        assert_eq!(indices, vec![0, 1, 2, 3, 4]);
    }
    #[test]
    fn wcsph_pressure_at_rest_density_near_zero() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1000.0;
        ps.update_wcsph_pressure(1000.0, 100.0, 7.0);
        assert!(
            ps.pressures[0].abs() < 1e-6,
            "P at rest density should be ~0"
        );
    }
    #[test]
    fn wcsph_pressure_compressed_positive() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1100.0;
        ps.update_wcsph_pressure(1000.0, 100.0, 7.0);
        assert!(
            ps.pressures[0] > 0.0,
            "compressed fluid must have positive pressure"
        );
    }
    #[test]
    fn compute_density_variable_h_positive() {
        let mut ps = make_cluster(8, 0.05);
        let hs = vec![0.2_f64; 8];
        ps.compute_density_variable_h(&hs);
        for &d in &ps.densities {
            assert!(d > 0.0, "density must be positive");
        }
    }
    #[test]
    fn cfl_timestep_proportional_to_h() {
        let ps = SphParticleSet::new();
        let dt1 = ps.cfl_timestep(0.1, 100.0, 0.4);
        let dt2 = ps.cfl_timestep(0.2, 100.0, 0.4);
        assert!(
            (dt2 - 2.0 * dt1).abs() < 1e-12,
            "dt must scale linearly with h"
        );
    }
    #[test]
    fn viscous_timestep_positive() {
        let ps = SphParticleSet::new();
        let dt = ps.viscous_timestep(0.1, 1e-4);
        assert!(dt > 0.0 && dt.is_finite());
    }
    #[test]
    fn adaptive_timestep_within_bounds() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        let dt = ps.adaptive_timestep(0.1, 100.0, 1e-4, 0.4, 1e-6, 0.01);
        assert!((1e-6..=0.01).contains(&dt), "dt = {dt} out of [1e-6, 0.01]");
    }
    #[test]
    fn accumulate_pressure_acceleration_no_nan() {
        let mut ps = make_cluster(8, 0.05);
        let hs = vec![0.2_f64; 8];
        ps.compute_density_variable_h(&hs);
        ps.update_wcsph_pressure(1000.0, 100.0, 7.0);
        ps.accumulate_pressure_acceleration(0.2);
        for a in &ps.accelerations {
            assert!(a[0].is_finite() && a[1].is_finite() && a[2].is_finite());
        }
    }
    #[test]
    fn accumulate_morris_viscosity_skipped_for_zero_mu() {
        let mut ps = make_cluster(4, 0.05);
        let hs = vec![0.2_f64; 4];
        ps.compute_density_variable_h(&hs);
        let acc_before: Vec<[f64; 3]> = ps.accelerations.clone();
        ps.accumulate_morris_viscosity(0.2, 0.0);
        for (a_b, a_a) in acc_before.iter().zip(ps.accelerations.iter()) {
            for k in 0..3 {
                assert!(
                    (a_b[k] - a_a[k]).abs() < 1e-30,
                    "zero mu should not change acc"
                );
            }
        }
    }
    #[test]
    fn accumulate_morris_viscosity_finite() {
        let mut ps = make_cluster(4, 0.05);
        let hs = vec![0.2_f64; 4];
        ps.compute_density_variable_h(&hs);
        ps.velocities[0] = [1.0, 0.0, 0.0];
        ps.velocities[1] = [-1.0, 0.0, 0.0];
        ps.accumulate_morris_viscosity(0.2, 1e-3);
        for a in &ps.accelerations {
            assert!(a[0].is_finite(), "viscosity acc should be finite");
        }
    }
}
#[cfg(test)]
mod tests_soa {
    use super::*;
    #[test]
    fn soa_new_is_empty() {
        let ps = ParticleSetSoA::new(10);
        assert!(ps.is_empty());
        assert_eq!(ps.len(), 0);
    }
    #[test]
    fn soa_add_particle() {
        let mut ps = ParticleSetSoA::new(4);
        ps.add_particle([1.0, 2.0, 3.0], [0.1, 0.0, 0.0], 2.5);
        assert_eq!(ps.len(), 1);
        assert!((ps.masses[0] - 2.5).abs() < 1e-15);
        assert_eq!(ps.forces[0], [0.0; 3]);
        assert_eq!(ps.densities[0], 0.0);
        assert_eq!(ps.pressures[0], 0.0);
    }
    #[test]
    fn soa_remove_particle_swap() {
        let mut ps = ParticleSetSoA::new(4);
        ps.add_particle([0.0; 3], [0.0; 3], 1.0);
        ps.add_particle([1.0, 0.0, 0.0], [0.0; 3], 2.0);
        ps.add_particle([2.0, 0.0, 0.0], [0.0; 3], 3.0);
        ps.remove_particle(0);
        assert_eq!(ps.len(), 2);
        assert!((ps.masses[0] - 3.0).abs() < 1e-15);
    }
    #[test]
    fn soa_apply_pbc_wraps_positive_overflow() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([1.05, 0.5, 0.5], [0.0; 3], 1.0);
        ps.apply_pbc([0.0; 3], [1.0; 3]);
        assert!(
            (ps.positions[0][0] - 0.05).abs() < 1e-12,
            "x={}",
            ps.positions[0][0]
        );
    }
    #[test]
    fn soa_apply_pbc_wraps_negative_overflow() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([-0.1, 0.5, 0.5], [0.0; 3], 1.0);
        ps.apply_pbc([0.0; 3], [1.0; 3]);
        assert!(
            (ps.positions[0][0] - 0.9).abs() < 1e-12,
            "x={}",
            ps.positions[0][0]
        );
    }
    #[test]
    fn soa_integrate_euler_constant_force() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([0.0; 3], [0.0; 3], 1.0);
        ps.forces[0] = [1.0, 0.0, 0.0];
        ps.integrate_euler(0.1);
        assert!(
            (ps.velocities[0][0] - 0.1).abs() < 1e-14,
            "vx={}",
            ps.velocities[0][0]
        );
        assert!(
            (ps.positions[0][0] - 0.01).abs() < 1e-14,
            "px={}",
            ps.positions[0][0]
        );
    }
    #[test]
    fn soa_leapfrog_kick_and_drift() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([0.0; 3], [1.0, 0.0, 0.0], 2.0);
        ps.forces[0] = [4.0, 0.0, 0.0];
        ps.integrate_leapfrog_kick(0.1);
        assert!((ps.velocities[0][0] - 1.2).abs() < 1e-14);
        ps.integrate_leapfrog_drift(0.1);
        assert!((ps.positions[0][0] - 0.12).abs() < 1e-14);
    }
    #[test]
    fn soa_center_of_mass_two_particles() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([0.0, 0.0, 0.0], [0.0; 3], 1.0);
        ps.add_particle([4.0, 0.0, 0.0], [0.0; 3], 3.0);
        let com = ps.center_of_mass();
        assert!((com[0] - 3.0).abs() < 1e-12, "com_x={}", com[0]);
        assert!(com[1].abs() < 1e-12);
        assert!(com[2].abs() < 1e-12);
    }
    #[test]
    fn soa_center_of_mass_empty() {
        let ps = ParticleSetSoA::new(0);
        let com = ps.center_of_mass();
        assert_eq!(com, [0.0; 3]);
    }
    #[test]
    fn soa_total_kinetic_energy_correct() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([0.0; 3], [3.0, 4.0, 0.0], 2.0);
        assert!((ps.total_kinetic_energy() - 25.0).abs() < 1e-12);
    }
    #[test]
    fn soa_total_kinetic_energy_zero_velocity() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([1.0, 2.0, 3.0], [0.0; 3], 5.0);
        assert!(ps.total_kinetic_energy().abs() < 1e-15);
    }
    #[test]
    fn soa_total_momentum() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([0.0; 3], [2.0, 0.0, 0.0], 3.0);
        ps.add_particle([0.0; 3], [-1.0, 0.0, 0.0], 1.0);
        let mom = ps.total_momentum();
        assert!((mom[0] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn soa_clear_forces() {
        let mut ps = ParticleSetSoA::new(2);
        ps.add_particle([0.0; 3], [0.0; 3], 1.0);
        ps.forces[0] = [9.0, 8.0, 7.0];
        ps.clear_forces();
        assert_eq!(ps.forces[0], [0.0; 3]);
    }
}
/// Cubic spline kernel value `W(r, h)` (3-D, normalized).
pub fn cubic_spline_w(r: f64, h: f64) -> f64 {
    let q = r / h;
    if q > 2.0 {
        return 0.0;
    }
    let sigma = 3.0 / (2.0 * std::f64::consts::PI * h * h * h);
    if q <= 1.0 {
        sigma * (2.0 / 3.0 - q * q + 0.5 * q * q * q)
    } else {
        let t = 2.0 - q;
        sigma * (t * t * t / 6.0)
    }
}
#[cfg(test)]
mod tests_particle_new {
    use super::*;
    #[test]
    fn filter_selects_fluid_particles_only() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, true);
        ps.add_raw([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let filter = ParticleFilter::new(&ps, |i| !ps.is_boundary[i]);
        assert_eq!(filter.count(), 2);
        assert!(filter.indices.contains(&0));
        assert!(filter.indices.contains(&2));
    }
    #[test]
    fn filter_empty_when_no_match() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        let filter = ParticleFilter::new(&ps, |_| false);
        assert!(filter.is_empty());
    }
    #[test]
    fn filter_kinetic_energy_matches_subset() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [3.0, 4.0, 0.0], 2.0, 0.0, false);
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        let filter = ParticleFilter::new(&ps, |i| i == 0);
        assert!((filter.kinetic_energy(&ps) - 25.0).abs() < 1e-12);
    }
    #[test]
    fn filter_center_of_mass_single_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([5.0, 3.0, 2.0], [0.0; 3], 1.0, 0.0, false);
        let filter = ParticleFilter::new(&ps, |_| true);
        let com = filter.center_of_mass(&ps);
        assert!((com[0] - 5.0).abs() < 1e-14);
        assert!((com[1] - 3.0).abs() < 1e-14);
    }
    #[test]
    fn filter_average_density_empty_is_zero() {
        let ps = SphParticleSet::new();
        let filter = ParticleFilter::new(&ps, |_| true);
        assert_eq!(filter.average_density(&ps), 0.0);
    }
    #[test]
    fn rms_speed_all_zero_velocities() {
        let mut ps = SphParticleSet::new();
        for _ in 0..5 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        }
        assert_eq!(ps.rms_speed(), 0.0);
    }
    #[test]
    fn rms_speed_single_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [3.0, 4.0, 0.0], 1.0, 0.0, false);
        assert!((ps.rms_speed() - 5.0).abs() < 1e-12);
    }
    #[test]
    fn rms_speed_empty_is_zero() {
        let ps = SphParticleSet::new();
        assert_eq!(ps.rms_speed(), 0.0);
    }
    #[test]
    fn mean_velocity_uniform_flow() {
        let mut ps = SphParticleSet::new();
        for _ in 0..4 {
            ps.add_raw([0.0; 3], [2.0, 0.0, 0.0], 1.0, 0.0, false);
        }
        let mv = ps.mean_velocity();
        assert!((mv[0] - 2.0).abs() < 1e-14);
        assert!(mv[1].abs() < 1e-14);
    }
    #[test]
    fn mean_velocity_empty_is_zero() {
        let ps = SphParticleSet::new();
        let mv = ps.mean_velocity();
        assert_eq!(mv, [0.0; 3]);
    }
    #[test]
    fn mean_pressure_uniform() {
        let mut ps = SphParticleSet::new();
        for _ in 0..4 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        }
        for p in &mut ps.pressures {
            *p = 500.0;
        }
        assert!((ps.mean_pressure() - 500.0).abs() < 1e-10);
    }
    #[test]
    fn max_pressure_returns_correct_value() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.pressures[0] = 100.0;
        ps.pressures[1] = 300.0;
        assert!((ps.max_pressure() - 300.0).abs() < 1e-14);
    }
    #[test]
    fn min_pressure_returns_correct_value() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.pressures[0] = 100.0;
        ps.pressures[1] = 300.0;
        assert!((ps.min_pressure() - 100.0).abs() < 1e-14);
    }
    #[test]
    fn pressure_std_zero_for_uniform_pressures() {
        let mut ps = SphParticleSet::new();
        for _ in 0..5 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        }
        for p in &mut ps.pressures {
            *p = 200.0;
        }
        assert!(ps.pressure_std() < 1e-12);
    }
    #[test]
    fn angular_momentum_about_origin_consistent_with_angular_momentum() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 2.0, 0.0, false);
        let l0 = ps.angular_momentum();
        let la = ps.angular_momentum_about([0.0; 3]);
        assert!((l0[0] - la[0]).abs() < 1e-14);
        assert!((l0[1] - la[1]).abs() < 1e-14);
        assert!((l0[2] - la[2]).abs() < 1e-14);
    }
    #[test]
    fn angular_momentum_about_pivot_shifts_correctly() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([2.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1.0, 0.0, false);
        let la = ps.angular_momentum_about([1.0, 0.0, 0.0]);
        assert!(
            (la[2] - 1.0).abs() < 1e-14,
            "L_z about pivot should be 1, got {}",
            la[2]
        );
    }
    #[test]
    fn translate_shifts_all_positions() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 2.0, 3.0], [0.0; 3], 1.0, 0.0, false);
        ps.translate([10.0, 20.0, 30.0]);
        assert!((ps.positions[0][0] - 11.0).abs() < 1e-14);
        assert!((ps.positions[0][1] - 22.0).abs() < 1e-14);
        assert!((ps.positions[0][2] - 33.0).abs() < 1e-14);
    }
    #[test]
    fn add_velocity_offset_moves_fluid_not_boundary() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, false);
        ps.add_raw([0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.0, true);
        ps.add_velocity_offset([5.0, 0.0, 0.0]);
        assert!(
            (ps.velocities[0][0] - 6.0).abs() < 1e-14,
            "Fluid velocity should increase"
        );
        assert!(
            (ps.velocities[1][0] - 1.0).abs() < 1e-14,
            "Boundary velocity must not change"
        );
    }
    #[test]
    fn set_uniform_density_sets_all_to_rho0() {
        let mut ps = SphParticleSet::new();
        for _ in 0..5 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        }
        ps.set_uniform_density(1000.0);
        for &d in &ps.densities {
            assert!((d - 1000.0).abs() < 1e-14);
        }
    }
    #[test]
    fn zero_pressures_clears_all() {
        let mut ps = SphParticleSet::new();
        for _ in 0..3 {
            ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        }
        for p in &mut ps.pressures {
            *p = 999.0;
        }
        ps.zero_pressures();
        for &p in &ps.pressures {
            assert_eq!(p, 0.0);
        }
    }
    #[test]
    fn update_pressures_tait_at_rest_density_near_zero() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1000.0;
        ps.update_pressures_tait(1000.0, 100.0, 7.0);
        assert!(
            ps.pressures[0].abs() < 1e-8,
            "Pressure at rest density should be ~0"
        );
    }
    #[test]
    fn update_pressures_tait_compressed_fluid_positive() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1100.0;
        ps.update_pressures_tait(1000.0, 100.0, 7.0);
        assert!(
            ps.pressures[0] > 0.0,
            "Compressed fluid should have positive pressure"
        );
    }
    #[test]
    fn box_reflection_bounces_particle_off_floor() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.5, -0.1, 0.5], [0.0, -1.0, 0.0], 1.0, 0.0, false);
        ps.apply_box_reflection([0.0; 3], [1.0; 3], 1.0);
        assert!(
            ps.positions[0][1] >= 0.0,
            "Particle should be above floor after reflection"
        );
        assert!(
            ps.velocities[0][1] >= 0.0,
            "Vertical velocity should be upward after reflection"
        );
    }
    #[test]
    fn box_reflection_particle_inside_not_moved() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.5, 0.5, 0.5], [1.0, 1.0, 1.0], 1.0, 0.0, false);
        ps.apply_box_reflection([0.0; 3], [1.0; 3], 1.0);
        assert!((ps.positions[0][0] - 0.5).abs() < 1e-14);
        assert!((ps.velocities[0][0] - 1.0).abs() < 1e-14);
    }
    #[test]
    fn box_reflection_boundary_particle_not_moved() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([-0.2, 0.5, 0.5], [-1.0, 0.0, 0.0], 1.0, 0.0, true);
        ps.apply_box_reflection([0.0; 3], [1.0; 3], 1.0);
        assert!(
            (ps.positions[0][0] - (-0.2)).abs() < 1e-14,
            "Boundary particle must not be reflected"
        );
    }
    #[test]
    fn sorted_by_distance_from_closest_first() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([10.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([5.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let sorted = ps.sorted_by_distance_from([0.0; 3]);
        assert_eq!(sorted[0].0, 1);
        assert_eq!(sorted[1].0, 2);
        assert_eq!(sorted[2].0, 0);
    }
    #[test]
    fn count_within_radius_correct() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.5, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        let count = ps.count_within_radius([0.0; 3], 1.0);
        assert_eq!(count, 2);
    }
    #[test]
    fn argsort_density_ascending_order() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities = vec![1200.0, 800.0, 1000.0];
        let idx = ps.argsort_by_field("density");
        assert_eq!(idx[0], 1);
        assert_eq!(idx[1], 2);
        assert_eq!(idx[2], 0);
    }
    #[test]
    fn argsort_unknown_field_returns_unsorted_indices() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        let idx = ps.argsort_by_field("unknown_field");
        assert_eq!(idx, vec![0, 1]);
    }
    #[test]
    fn check_invariant_consistent_set_is_ok() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        assert!(ps.check_invariant().is_ok());
    }
    #[test]
    fn check_invariant_empty_set_is_ok() {
        let ps = SphParticleSet::new();
        assert!(ps.check_invariant().is_ok());
    }
    #[test]
    fn ensemble_stats_initial_state_is_empty() {
        let stats = EnsembleStats::new(1000.0);
        assert!(stats.is_empty());
        assert_eq!(stats.count, 0);
        assert_eq!(stats.mean_ke, 0.0);
    }
    #[test]
    fn ensemble_stats_accumulate_single_snapshot() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [2.0, 0.0, 0.0], 2.0, 0.0, false);
        ps.densities[0] = 1000.0;
        let mut stats = EnsembleStats::new(1000.0);
        stats.accumulate(&ps);
        assert_eq!(stats.count, 1);
        assert!((stats.mean_ke - 4.0).abs() < 1e-12);
    }
    #[test]
    fn ensemble_stats_accumulate_two_snapshots_averages_ke() {
        let mut ps1 = SphParticleSet::new();
        ps1.add_raw([0.0; 3], [2.0, 0.0, 0.0], 2.0, 0.0, false);
        ps1.densities[0] = 1000.0;
        let mut ps2 = SphParticleSet::new();
        ps2.add_raw([0.0; 3], [4.0, 0.0, 0.0], 2.0, 0.0, false);
        ps2.densities[0] = 1000.0;
        let mut stats = EnsembleStats::new(1000.0);
        stats.accumulate(&ps1);
        stats.accumulate(&ps2);
        assert_eq!(stats.count, 2);
        assert!(
            (stats.mean_ke - 10.0).abs() < 1e-10,
            "mean_ke should be 10, got {}",
            stats.mean_ke
        );
    }
    #[test]
    fn ensemble_stats_density_error_increases_on_deviation() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1200.0;
        let mut stats = EnsembleStats::new(1000.0);
        stats.accumulate(&ps);
        assert!((stats.max_density_error - 0.2).abs() < 1e-12);
    }
    #[test]
    fn cubic_spline_w_zero_outside_support() {
        assert_eq!(cubic_spline_w(0.21, 0.1), 0.0);
    }
    #[test]
    fn cubic_spline_w_positive_at_origin() {
        assert!(cubic_spline_w(0.0, 0.1) > 0.0);
    }
    #[test]
    fn pressure_gradient_force_zero_isolated_particle() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.0; 3], 1.0, 0.0, false);
        ps.densities[0] = 1000.0;
        ps.pressures[0] = 100.0;
        let f = ps.pressure_gradient_force(0, 0.1);
        assert!(f[0].abs() < 1e-14, "Isolated particle force must be zero");
    }
    #[test]
    fn count_in_box_counts_correctly() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.5, 0.5, 0.5], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([1.5, 0.5, 0.5], [0.0; 3], 1.0, 0.0, false);
        let count = ps.count_in_box([0.0; 3], [1.0; 3]);
        assert_eq!(count, 1);
    }
    #[test]
    fn count_positive_x_half_correct() {
        let mut ps = SphParticleSet::new();
        ps.add_raw([1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        ps.add_raw([-1.0, 0.0, 0.0], [0.0; 3], 1.0, 0.0, false);
        assert_eq!(ps.count_positive_x_half(), 1);
    }
    #[test]
    fn verlet_position_free_fall_correct_displacement() {
        let dt = 0.01;
        let a = [0.0, -10.0, 0.0];
        let x0 = [0.0, 10.0, 0.0];
        let mut ps = SphParticleSet::new();
        ps.add_raw(x0, [0.0; 3], 1.0, 0.0, false);
        ps.accelerations[0] = a;
        let x_prev = [
            x0[0] - 0.0 * dt + 0.5 * a[0] * dt * dt,
            x0[1] - 0.0 * dt + 0.5 * a[1] * dt * dt,
            x0[2] - 0.0 * dt + 0.5 * a[2] * dt * dt,
        ];
        let mut prev = vec![x_prev];
        ps.integrate_verlet_position(&mut prev, dt);
        assert!(
            ps.positions[0][1] < x0[1],
            "Particle should fall under gravity"
        );
    }
}
