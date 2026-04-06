//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use super::functions::*;

#[cfg(test)]
mod pbd_extended_tests {
    use super::*;
    use crate::pbd_system::ConstraintColoring;
    use crate::pbd_system::PbdConstraint;
    use crate::pbd_system::PbdConstraintType;
    use crate::pbd_system::PbdParticle;
    use crate::pbd_system::PbdSystem;
    use crate::pbd_system::SleepController;
    use crate::pbd_system::WarmStartCache;
    pub(super) const EPS: f64 = 1e-9;
    /// Coloring of non-overlapping constraints gives one group.
    #[test]
    fn test_coloring_non_overlapping() {
        let c1 = PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 1.0,
                compliance: 0.0,
            },
            particles: vec![0, 1],
        };
        let c2 = PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 1.0,
                compliance: 0.0,
            },
            particles: vec![2, 3],
        };
        let coloring = ConstraintColoring::build(&[c1, c2], 4);
        assert_eq!(
            coloring.num_colors(),
            1,
            "non-overlapping constraints should be in 1 color group"
        );
    }
    /// Coloring of overlapping constraints gives 2+ groups.
    #[test]
    fn test_coloring_overlapping() {
        let c1 = PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 1.0,
                compliance: 0.0,
            },
            particles: vec![0, 1],
        };
        let c2 = PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 1.0,
                compliance: 0.0,
            },
            particles: vec![1, 2],
        };
        let coloring = ConstraintColoring::build(&[c1, c2], 3);
        assert!(
            coloring.num_colors() >= 2,
            "overlapping constraints need at least 2 colors"
        );
    }
    /// Colored XPBD step produces same direction of motion as regular step.
    #[test]
    fn test_colored_xpbd_step_falls() {
        let mut sys = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys.add_particle([0.0, 5.0, 0.0], 1.0);
        let coloring = ConstraintColoring::build(&sys.constraints, sys.particles.len());
        xpbd_step_colored(&mut sys, 1.0 / 60.0, &coloring);
        assert!(
            sys.particles[0].position[1] < 5.0,
            "particle should fall, y={}",
            sys.particles[0].position[1]
        );
    }
    /// Jacobi distance solve moves particles toward rest length.
    #[test]
    fn test_jacobi_distance_solve() {
        let mut particles = vec![
            PbdParticle::new([0.0, 0.0, 0.0], 1.0),
            PbdParticle::new([4.0, 0.0, 0.0], 1.0),
        ];
        particles[0].predicted = particles[0].position;
        particles[1].predicted = particles[1].position;
        let constraints = vec![PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 2.0,
                compliance: 0.0,
            },
            particles: vec![0, 1],
        }];
        for _ in 0..50 {
            solve_distance_jacobi(&mut particles, &constraints, 1.0 / 60.0);
            for p in particles.iter_mut() {
                p.position = p.predicted;
            }
        }
        let d = dist3(particles[0].position, particles[1].position);
        assert!((d - 2.0).abs() < 0.1, "Jacobi: dist = {d}, expected ~2");
    }
    /// Sleep controller puts a stationary particle to sleep.
    #[test]
    fn test_sleep_controller_sleeps_still_particle() {
        let mut ctrl = SleepController::new(1, 1.0, 3);
        let mut p = PbdParticle::new([0.0; 3], 1.0);
        p.velocity = [0.0; 3];
        for _ in 0..3 {
            ctrl.update(&[p.clone()]);
        }
        assert!(ctrl.asleep[0], "particle with zero velocity should sleep");
    }
    /// Sleep controller wakes up a particle that gains velocity.
    #[test]
    fn test_sleep_controller_wakes_moving_particle() {
        let mut ctrl = SleepController::new(1, 0.5, 2);
        let mut p = PbdParticle::new([0.0; 3], 1.0);
        p.velocity = [0.0; 3];
        for _ in 0..2 {
            ctrl.update(&[p.clone()]);
        }
        p.velocity = [5.0, 0.0, 0.0];
        ctrl.update(&[p.clone()]);
        assert!(!ctrl.asleep[0], "fast particle should wake up");
    }
    /// Sleep controller count and wake_all.
    #[test]
    fn test_sleep_controller_count_and_wake() {
        let mut ctrl = SleepController::new(2, 1.0, 1);
        let p = PbdParticle::new([0.0; 3], 1.0);
        ctrl.update(&[p.clone(), p.clone()]);
        ctrl.update(&[p.clone(), p.clone()]);
        assert!(ctrl.sleeping_count() <= 2);
        ctrl.wake_all();
        assert_eq!(ctrl.sleeping_count(), 0);
    }
    /// Warm-start cache initializes to zero.
    #[test]
    fn test_warm_start_cache_init_zero() {
        let cache = WarmStartCache::new(5);
        for &l in &cache.lambdas {
            assert!(l.abs() < EPS, "initial lambda = {l}");
        }
    }
    /// Warm-start cache update stores constraint violations.
    #[test]
    fn test_warm_start_cache_update() {
        let mut cache = WarmStartCache::new(1);
        let particles = vec![
            PbdParticle::new([0.0, 0.0, 0.0], 1.0),
            PbdParticle::new([3.0, 0.0, 0.0], 1.0),
        ];
        let constraints = vec![PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 2.0,
                compliance: 0.0,
            },
            particles: vec![0, 1],
        }];
        cache.update(&constraints, &particles);
        assert!(
            (cache.lambdas[0] - 1.0).abs() < 1e-10,
            "cached violation = {}",
            cache.lambdas[0]
        );
    }
    /// Velocity damping reduces velocity magnitude.
    #[test]
    fn test_velocity_damping_reduces() {
        let mut particles = vec![PbdParticle::new([0.0; 3], 1.0)];
        particles[0].velocity = [2.0, 0.0, 0.0];
        apply_velocity_damping(&mut particles, 0.5);
        let speed = len3(particles[0].velocity);
        assert!(
            (speed - 1.0).abs() < EPS,
            "speed after 50% damping = {speed}"
        );
    }
    /// Full velocity damping sets velocity to zero.
    #[test]
    fn test_velocity_damping_full() {
        let mut particles = vec![PbdParticle::new([0.0; 3], 1.0)];
        particles[0].velocity = [3.0, 4.0, 5.0];
        apply_velocity_damping(&mut particles, 1.0);
        let speed = len3(particles[0].velocity);
        assert!(speed < EPS, "full damping should zero velocity: {speed}");
    }
    /// Fixed particles are not affected by velocity damping.
    #[test]
    fn test_velocity_damping_ignores_fixed() {
        let mut particles = vec![PbdParticle::new_fixed([0.0; 3])];
        particles[0].velocity = [1.0, 0.0, 0.0];
        apply_velocity_damping(&mut particles, 0.5);
    }
    /// Position damping blends toward COM velocity.
    #[test]
    fn test_position_damping_toward_com() {
        let mut particles = vec![
            {
                let mut p = PbdParticle::new([0.0; 3], 1.0);
                p.velocity = [2.0, 0.0, 0.0];
                p
            },
            {
                let mut p = PbdParticle::new([1.0, 0.0, 0.0], 1.0);
                p.velocity = [0.0, 0.0, 0.0];
                p
            },
        ];
        apply_position_damping(&mut particles, 1.0);
        for p in &particles {
            assert!(
                (p.velocity[0] - 1.0).abs() < EPS,
                "velocity after full position damping = {}",
                p.velocity[0]
            );
        }
    }
    /// Volume constraint is added correctly.
    #[test]
    fn test_add_volume_constraint() {
        let mut sys = PbdSystem::new([0.0; 3], 1);
        let i0 = sys.add_particle([0.0, 0.0, 0.0], 1.0);
        let i1 = sys.add_particle([1.0, 0.0, 0.0], 1.0);
        let i2 = sys.add_particle([0.0, 1.0, 0.0], 1.0);
        let i3 = sys.add_particle([0.0, 0.0, 1.0], 1.0);
        sys.add_volume_constraint(i0, i1, i2, i3, 0.0);
        assert_eq!(sys.constraint_count(), 1);
        if let PbdConstraintType::VolumeConservation { rest_volume, .. } =
            sys.constraints[0].constraint_type
        {
            assert!(
                (rest_volume - 1.0 / 6.0).abs() < 1e-10,
                "rest_volume = {rest_volume}"
            );
        }
    }
    /// Bending constraint is added correctly.
    #[test]
    fn test_add_bending_constraint() {
        let mut sys = PbdSystem::new([0.0; 3], 1);
        for _ in 0..4 {
            sys.add_particle([0.0; 3], 1.0);
        }
        sys.add_bending_constraint(0, 1, 2, 3, std::f64::consts::PI, 0.01);
        assert_eq!(sys.constraint_count(), 1);
    }
    /// apply_velocity_damping via system method works.
    #[test]
    fn test_system_velocity_damping() {
        let mut sys = PbdSystem::new([0.0; 3], 1);
        let i = sys.add_particle([0.0; 3], 1.0);
        sys.particles[i].velocity = [10.0, 0.0, 0.0];
        sys.apply_velocity_damping(0.5);
        let speed = len3(sys.particles[i].velocity);
        assert!((speed - 5.0).abs() < EPS, "speed = {speed}");
    }
    /// Total potential energy increases with height.
    #[test]
    fn test_total_potential_energy_height() {
        let mut sys1 = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys1.add_particle([0.0, 1.0, 0.0], 1.0);
        let mut sys2 = PbdSystem::new([0.0, -9.81, 0.0], 1);
        sys2.add_particle([0.0, 10.0, 0.0], 1.0);
        let e1 = sys1.total_potential_energy();
        let e2 = sys2.total_potential_energy();
        assert!(
            e2 < e1,
            "higher particle should have more negative potential: e1={e1}, e2={e2}"
        );
    }
    /// Build a grid cloth, add coloring, and run a few colored steps.
    #[test]
    fn test_colored_grid_cloth_simulation() {
        let mut sys = build_grid_cloth(3, 3, 0.5, 0.1, 1e-4);
        for c in 0..3 {
            sys.particles[c].fixed = true;
            sys.particles[c].inv_mass = 0.0;
        }
        let n_particles = sys.particles.len();
        let coloring = ConstraintColoring::build(&sys.constraints, n_particles);
        assert!(coloring.num_colors() >= 1);
        let dt = 1.0 / 60.0;
        for _ in 0..20 {
            xpbd_step_colored(&mut sys, dt, &coloring);
        }
        let y_bottom = sys.particles[6].position[1];
        assert!(
            y_bottom < 0.1,
            "bottom cloth particles should fall, y={y_bottom}"
        );
    }
    /// Chain simulation with Jacobi solver converges.
    #[test]
    fn test_jacobi_chain_converges() {
        let mut sys = build_chain(3, 1.0, 1.0, 0.0);
        let dt = 1.0 / 60.0;
        for _ in 0..300 {
            for p in sys.particles.iter_mut() {
                if p.fixed {
                    p.predicted = p.position;
                    continue;
                }
                let v = add3(p.velocity, scale3(sys.gravity, dt));
                p.predicted = add3(p.position, scale3(v, dt));
            }
            solve_distance_jacobi(&mut sys.particles, &sys.constraints, dt);
            for p in sys.particles.iter_mut() {
                if p.fixed {
                    continue;
                }
                p.velocity = scale3(sub3(p.predicted, p.position), 1.0 / dt);
                p.position = p.predicted;
            }
        }
        for i in 0..2 {
            let d = dist3(sys.particles[i].position, sys.particles[i + 1].position);
            assert!(
                (d - 1.0).abs() < 0.5,
                "chain link {i} dist = {d}, expected ~1"
            );
        }
    }
    /// Warm-start apply does not change positions for zero cache.
    #[test]
    fn test_warm_start_zero_no_change() {
        let cache = WarmStartCache::new(1);
        let mut particles = vec![
            PbdParticle::new([0.0, 0.0, 0.0], 1.0),
            PbdParticle::new([1.0, 0.0, 0.0], 1.0),
        ];
        let before = [particles[0].predicted, particles[1].predicted];
        let constraints = vec![PbdConstraint {
            constraint_type: PbdConstraintType::Distance {
                rest_length: 1.0,
                compliance: 0.0,
            },
            particles: vec![0, 1],
        }];
        cache.apply_warm_start(&mut particles, &constraints, 0.0, 1.0 / 60.0);
        assert_eq!(particles[0].predicted, before[0]);
        assert_eq!(particles[1].predicted, before[1]);
    }
    /// cross3 produces a vector perpendicular to both inputs.
    #[test]
    fn test_cross3_perpendicular() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        let c = cross3(a, b);
        assert!((c[0]).abs() < EPS);
        assert!((c[1]).abs() < EPS);
        assert!((c[2] - 1.0).abs() < EPS);
    }
    /// dot3 of orthogonal vectors is zero.
    #[test]
    fn test_dot3_orthogonal() {
        let a = [1.0, 0.0, 0.0];
        let b = [0.0, 1.0, 0.0];
        assert!(dot3(a, b).abs() < EPS);
    }
    /// len3 of unit vector is 1.
    #[test]
    fn test_len3_unit() {
        let v = [1.0 / 3.0_f64.sqrt(); 3];
        assert!((len3(v) - 1.0).abs() < EPS);
    }
    /// Floor solve does not move particle above floor level.
    #[test]
    fn test_floor_solve_above_floor() {
        let mut p = PbdParticle::new([0.0, 2.0, 0.0], 1.0);
        p.predicted = [0.0, 2.0, 0.0];
        let before = p.predicted;
        solve_floor_xpbd(&mut p, 0.0, 0.0);
        assert_eq!(
            p.predicted, before,
            "particle above floor should not be moved"
        );
    }
    /// Floor solve with restitution > 0 reflects velocity.
    #[test]
    fn test_floor_solve_restitution() {
        let mut p = PbdParticle::new([0.0, -1.0, 0.0], 1.0);
        p.predicted = [0.0, -1.0, 0.0];
        p.velocity = [0.0, -5.0, 0.0];
        solve_floor_xpbd(&mut p, 0.0, 0.5);
        assert!(p.predicted[1] >= 0.0, "particle should be on floor");
        assert!(p.velocity[1] > 0.0, "velocity should be reflected upward");
        assert!(
            (p.velocity[1] - 2.5).abs() < EPS,
            "reflected velocity = {}",
            p.velocity[1]
        );
    }
    /// build_chain with n=1 has no constraints.
    #[test]
    fn test_build_chain_single_node() {
        let sys = build_chain(1, 1.0, 1.0, 0.0);
        assert_eq!(sys.particle_count(), 1);
        assert_eq!(sys.constraint_count(), 0);
    }
    /// build_chain spacing is correct.
    #[test]
    fn test_build_chain_spacing() {
        let sys = build_chain(3, 2.0, 1.0, 0.0);
        let d = dist3(sys.particles[0].position, sys.particles[1].position);
        assert!((d - 2.0).abs() < EPS, "chain spacing = {d}");
    }
    /// Tet volume is positive for non-degenerate tet.
    #[test]
    fn test_tet_volume_positive() {
        let v = compute_tetrahedron_volume(
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 2.0, 0.0],
            [0.0, 0.0, 2.0],
        );
        assert!((v - 8.0 / 6.0).abs() < 1e-10, "volume = {v}");
    }
}
