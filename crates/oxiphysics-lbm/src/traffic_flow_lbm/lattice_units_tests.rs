//! Regression tests for the lattice-unit handling of [`TrafficLbm`]: physical
//! speeds (m/s) must enter the D1Q3 equilibrium as lattice velocities
//! `v * dt / dx`, and velocities and flows must be reported in physical units.

use super::{
    AccidentPropagation, DiagramModel, FundamentalDiagram, SpeedLimitZone, TrafficLbm,
    TrafficSignal, shock_wave_speed,
};
use crate::error::Error;

fn greenshields() -> FundamentalDiagram {
    FundamentalDiagram::new(0.1, 20.0, DiagramModel::Greenshields)
}

/// Fractional cell index at which the density profile first rises through
/// `level`, searching from cell `from` onwards (linear interpolation).
fn rising_crossing(lbm: &TrafficLbm, level: f64, from: usize) -> Option<f64> {
    (from..lbm.nx - 1).find_map(|i| {
        let (a, b) = (lbm.density(i), lbm.density(i + 1));
        if a < level && b >= level {
            Some(i as f64 + (level - a) / (b - a))
        } else {
            None
        }
    })
}

#[test]
fn equilibrium_velocity_and_flow_are_reported_in_physical_units() {
    let rho = 0.03;
    let speed = greenshields().speed(rho);
    assert!((speed - 14.0).abs() < 1e-12);
    // The same physical state on three different lattices: a uniform
    // equilibrium is stationary and must read back as 14 m/s everywhere.
    for (dx, dt) in [(1.0, 0.02), (2.0, 0.04), (0.5, 0.01)] {
        let mut lbm = TrafficLbm::new(40, dx, dt, 1.0, greenshields());
        lbm.init_density(|_| rho);
        lbm.run(25);
        for i in [0, 17, 39] {
            assert!((lbm.density(i) - rho).abs() < 1e-12, "dx={dx} dt={dt}");
            assert!(
                (lbm.velocity(i) - speed).abs() < 1e-10,
                "dx={dx} dt={dt}: velocity {} vs {speed}",
                lbm.velocity(i)
            );
            assert!(
                (lbm.flow(i) - rho * speed).abs() < 1e-12,
                "dx={dx} dt={dt}: flow {} vs {}",
                lbm.flow(i),
                rho * speed
            );
        }
        assert!((lbm.time - 25.0 * dt).abs() < 1e-12);
    }
}

#[test]
fn a_density_perturbation_travels_at_the_characteristic_speed() {
    // q(rho) = 20 rho (1 - 10 rho): the characteristic speed at rho0 = 0.02 is
    // q'(rho0) = 20 (1 - 20 rho0) = 12 m/s, independent of the lattice.
    let rho0 = 0.02;
    let amplitude = 2.0e-4;
    let characteristic_speed = 12.0;
    let travel_time = 8.0;
    for (nx, dx, dt) in [(400usize, 1.0, 0.02), (800, 0.5, 0.01)] {
        let mut lbm = TrafficLbm::new(nx, dx, dt, 1.0, greenshields());
        let centre = 100.0;
        let sigma = 6.0;
        lbm.init_density(|x| {
            rho0 + amplitude * (-(x - centre).powi(2) / (2.0 * sigma * sigma)).exp()
        });
        let centroid = |lbm: &TrafficLbm| -> f64 {
            let (mut moment, mut mass) = (0.0, 0.0);
            for i in 0..lbm.nx {
                let excess = lbm.density(i) - rho0;
                moment += (i as f64 + 0.5) * dx * excess;
                mass += excess;
            }
            moment / mass
        };
        let start = centroid(&lbm);
        lbm.run((travel_time / dt).round() as usize);
        let shift = centroid(&lbm) - start;
        let expected = characteristic_speed * travel_time;
        assert!(
            (shift - expected).abs() < 0.02 * expected,
            "dx={dx} dt={dt}: the perturbation moved {shift} m, expected {expected} m"
        );
    }
}

#[test]
fn a_jam_shock_moves_at_the_rankine_hugoniot_speed() {
    let diagram = greenshields();
    let (rho_up, rho_down) = (0.01, 0.05);
    let speed = shock_wave_speed(&diagram, rho_up, rho_down);
    assert!((speed - 8.0).abs() < 1e-12, "q(0.05) - q(0.01) over 0.04");
    let level = 0.5 * (rho_up + rho_down);
    for (nx, dx, dt, steps) in [(800usize, 1.0, 0.02, 1000usize), (1600, 0.5, 0.01, 2000)] {
        let mut lbm = TrafficLbm::new(nx, dx, dt, 1.0, diagram.clone());
        let interface = nx / 2;
        lbm.init_density(|x| {
            if x < interface as f64 * dx {
                rho_up
            } else {
                rho_down
            }
        });
        let vehicles = lbm.total_vehicles();
        let start = rising_crossing(&lbm, level, 0).expect("initial jump");
        lbm.run(steps);
        let end = rising_crossing(&lbm, level, nx / 4).expect("the shock must still exist");
        let travelled = (end - start) * dx;
        let expected = speed * steps as f64 * dt;
        assert!(
            (travelled - expected).abs() < 0.05 * expected,
            "dx={dx} dt={dt}: the shock moved {travelled} m, Rankine-Hugoniot gives {expected} m"
        );
        assert!(
            (lbm.total_vehicles() - vehicles).abs() < 1e-9 * vehicles,
            "the scheme must conserve the vehicles"
        );
    }
}

#[test]
fn the_free_flow_lattice_velocity_and_cfl_check_follow_dt_over_dx() {
    let lbm = TrafficLbm::new(10, 1.0, 0.02, 1.0, greenshields());
    assert!((lbm.free_flow_lattice_velocity() - 0.4).abs() < 1e-14);
    assert!(lbm.cfl_ok());
    let edge = TrafficLbm::new(10, 1.0, 0.028, 1.0, greenshields());
    assert!(edge.cfl_ok(), "0.56 < 1/sqrt(3)");
    let supersonic = TrafficLbm::new(10, 1.0, 0.029, 1.0, greenshields());
    assert!(!supersonic.cfl_ok(), "0.58 > 1/sqrt(3)");
    // The same dt on a coarser lattice is fine again.
    let coarse = TrafficLbm::new(10, 2.0, 0.029, 1.0, greenshields());
    assert!(coarse.cfl_ok());
    let degenerate = TrafficLbm::new(10, 0.0, 0.02, 1.0, greenshields());
    assert!(!degenerate.cfl_ok());
    assert!((TrafficLbm::max_lattice_velocity() - 1.0 / 3.0_f64.sqrt()).abs() < 1e-15);
}

#[test]
fn try_new_accepts_a_stable_lattice() {
    let lbm =
        TrafficLbm::try_new(50, 1.0, 0.02, 1.0, greenshields()).expect("u_lat = 0.4 is subsonic");
    assert_eq!(lbm.nx, 50);
    assert!(lbm.cfl_ok());
    assert!(lbm.f.iter().all(|f| f.iter().all(|&v| v == 0.0)));
}

#[test]
fn try_new_rejects_a_supersonic_lattice_velocity() {
    match TrafficLbm::try_new(50, 1.0, 0.04, 1.0, greenshields()) {
        Err(Error::ParameterOutOfBounds {
            name,
            value,
            min,
            max,
        }) => {
            assert_eq!(name, "v_free * dt / dx");
            assert!((value - 0.8).abs() < 1e-12, "value {value}");
            assert_eq!(min, 0.0);
            assert!((max - 1.0 / 3.0_f64.sqrt()).abs() < 1e-15, "max {max}");
        }
        other => panic!("expected the lattice-velocity error, got {other:?}"),
    }
}

#[test]
fn try_new_rejects_invalid_geometry_time_step_and_relaxation_time() {
    let expect_out_of_bounds = |result: Result<TrafficLbm, Error>, expected_name: &str| match result
    {
        Err(Error::ParameterOutOfBounds { name, .. }) => assert_eq!(name, expected_name),
        other => panic!("expected an out-of-bounds error for {expected_name}, got {other:?}"),
    };
    expect_out_of_bounds(
        TrafficLbm::try_new(50, 0.0, 0.02, 1.0, greenshields()),
        "dx",
    );
    expect_out_of_bounds(
        TrafficLbm::try_new(50, -1.0, 0.02, 1.0, greenshields()),
        "dx",
    );
    expect_out_of_bounds(
        TrafficLbm::try_new(50, f64::NAN, 0.02, 1.0, greenshields()),
        "dx",
    );
    expect_out_of_bounds(TrafficLbm::try_new(50, 1.0, 0.0, 1.0, greenshields()), "dt");
    expect_out_of_bounds(
        TrafficLbm::try_new(50, 1.0, f64::INFINITY, 1.0, greenshields()),
        "dt",
    );

    let mut broken = greenshields();
    broken.rho_jam = 0.0;
    expect_out_of_bounds(TrafficLbm::try_new(50, 1.0, 0.02, 1.0, broken), "rho_jam");
    let mut broken = greenshields();
    broken.v_free = -1.0;
    expect_out_of_bounds(TrafficLbm::try_new(50, 1.0, 0.02, 1.0, broken), "v_free");

    assert!(matches!(
        TrafficLbm::try_new(0, 1.0, 0.02, 1.0, greenshields()),
        Err(Error::InvalidGridDimensions { nx: 0, ny: 1 })
    ));
    assert!(matches!(
        TrafficLbm::try_new(50, 1.0, 0.02, 0.5, greenshields()),
        Err(Error::UnstableRelaxationTime { tau }) if tau == 0.5
    ));
    assert!(matches!(
        TrafficLbm::try_new(50, 1.0, 0.02, f64::NAN, greenshields()),
        Err(Error::UnstableRelaxationTime { .. })
    ));
}

#[test]
fn speed_limits_signals_and_accidents_act_on_the_physical_speed() {
    let rho = 0.03;
    let dt = 0.02;
    let mut lbm = TrafficLbm::new(30, 1.0, dt, 1.0, greenshields());
    lbm.init_density(|_| rho);
    // A zone limit of 5 m/s (below the 14 m/s of the diagram) on cells 5..=9.
    lbm.add_speed_zone(SpeedLimitZone::new(5, 9, 5.0, 1.0));
    // A red signal on cell 15 and an accident that halves the speed on cell 25.
    lbm.add_signal(TrafficSignal::new(15, 1.0, 1.0, 100.0, 2.5));
    lbm.add_accident(AccidentPropagation::new(25, 0.5, 1000.0));
    lbm.collide();

    // With tau = 1 the collision output is the equilibrium of the modified speed.
    let free = greenshields().speed(rho);
    for i in [0, 12, 20] {
        assert!((lbm.velocity(i) - free).abs() < 1e-10, "cell {i}");
    }
    for i in 5..=9 {
        assert!((lbm.velocity(i) - 5.0).abs() < 1e-10, "zone cell {i}");
        assert!((lbm.flow(i) - rho * 5.0).abs() < 1e-12, "zone cell {i}");
    }
    assert!(
        lbm.velocity(15).abs() < 1e-12,
        "red signal stops the traffic"
    );
    assert!(
        (lbm.velocity(25) - 0.5 * free).abs() < 1e-10,
        "accident cell: {}",
        lbm.velocity(25)
    );
}

#[test]
fn the_vehicle_count_is_conserved_with_physical_units() {
    let mut lbm = TrafficLbm::new(120, 2.0, 0.04, 1.0, greenshields());
    lbm.init_density(|x| 0.03 + 0.02 * (2.0 * std::f64::consts::PI * x / 240.0).sin());
    let before = lbm.total_vehicles();
    lbm.run(200);
    assert!((lbm.total_vehicles() - before).abs() < 1e-9 * before);
    assert!(lbm.density_profile().iter().all(|&rho| rho > 0.0));
}
