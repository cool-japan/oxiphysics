//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[cfg(test)]
use super::types::{DemContact, DemParticle, DemSimulation, GranularContactModel};

#[cfg(test)]
mod tests_dem_contact {
    use super::*;
    fn make_contact(overlap: f64) -> DemContact {
        DemContact::new(0, 1, overlap, [1.0, 0.0, 0.0], 1.0e9, 0.05, 0.5, 0.01)
    }
    #[test]
    fn hertz_contact_force_zero_overlap_returns_zero() {
        let c = make_contact(0.0);
        let f = c.compute_hertz_contact_force();
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn hertz_contact_force_negative_overlap_returns_zero() {
        let c = make_contact(-0.001);
        let f = c.compute_hertz_contact_force();
        assert_eq!(f, [0.0; 3]);
    }
    #[test]
    fn hertz_contact_force_positive_overlap_repulsive_on_i() {
        let c = make_contact(0.001);
        let f = c.compute_hertz_contact_force();
        assert!(
            f[0] < 0.0,
            "Force on i should be in −x (repulsive), got {}",
            f[0]
        );
        assert_eq!(f[1], 0.0);
        assert_eq!(f[2], 0.0);
    }
    #[test]
    fn hertz_contact_force_scales_as_delta_to_3_2() {
        let c1 = make_contact(0.001);
        let c2 = make_contact(0.002);
        let f1 = c1.compute_hertz_contact_force()[0].abs();
        let f2 = c2.compute_hertz_contact_force()[0].abs();
        let ratio = f2 / f1;
        let expected = 2.0_f64.powf(1.5);
        assert!(
            (ratio - expected).abs() < 1e-9,
            "ratio={ratio} expected={expected}"
        );
    }
    #[test]
    fn rolling_resistance_zero_relative_angular_velocity() {
        let c = make_contact(0.001);
        let tau = c.compute_rolling_resistance([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(tau, [0.0; 3], "Zero relative ω → zero torque");
    }
    #[test]
    fn rolling_resistance_zero_overlap_returns_zero() {
        let c = make_contact(0.0);
        let tau = c.compute_rolling_resistance([0.0, 0.0, 5.0], [0.0, 0.0, 0.0]);
        assert_eq!(tau, [0.0; 3]);
    }
    #[test]
    fn rolling_resistance_direction_opposes_relative_spin() {
        let c = make_contact(0.001);
        let tau = c.compute_rolling_resistance([0.0, 0.0, 1.0], [0.0; 3]);
        assert!(
            tau[2] < 0.0,
            "Torque should oppose relative spin, tau_z={}",
            tau[2]
        );
    }
    #[test]
    fn coordination_number_no_contact() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            50.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        let p0 = DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0);
        let p1 = DemParticle::new_sphere([10.0, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle(p0);
        sim.add_particle(p1);
        let z = sim.compute_coordination_number(0.0);
        assert_eq!(z, 0.0, "Far-apart particles have Z=0, got {z}");
    }
    #[test]
    fn coordination_number_one_contact_two_particles() {
        let mut sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            50.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        let p0 = DemParticle::new_sphere([0.0; 3], [0.0; 3], 0.1, 1.0);
        let p1 = DemParticle::new_sphere([0.1, 0.0, 0.0], [0.0; 3], 0.1, 1.0);
        sim.add_particle(p0);
        sim.add_particle(p1);
        let z = sim.compute_coordination_number(0.0);
        assert!((z - 1.0).abs() < 1e-12, "Expected Z̄=1.0, got {z}");
    }
    #[test]
    fn coordination_number_empty_simulation_is_zero() {
        let sim = DemSimulation::new(
            0.001,
            [0.0, -9.81, 0.0],
            1e5,
            5e4,
            50.0,
            0.3,
            GranularContactModel::LinearSpring,
            1e9,
        );
        let z = sim.compute_coordination_number(0.0);
        assert_eq!(z, 0.0);
    }
}
