//! Regression tests for the sign convention of the Boussinesq buoyancy
//! coupling: `F = -rho * beta * (T - T_ref) * g_vec`, so that fluid hotter than
//! the reference temperature is pushed against the gravity vector.

use super::types::{BoussinesqCoupling, ThermalLbm};
use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

const GRAVITY_DOWN: [f64; 2] = [0.0, -1.0e-2];

#[test]
fn hot_fluid_is_pushed_against_the_gravity_vector() {
    let coupling = BoussinesqCoupling::new(GRAVITY_DOWN, 2.0e-3, 1.0);
    let force = coupling.buoyancy_force(1.2, 1.5);
    let expected = 1.2 * 2.0e-3 * 0.5 * 1.0e-2;
    assert!(force[1] > 0.0, "hot fluid must rise, F_y = {}", force[1]);
    assert!(
        (force[1] - expected).abs() < 1e-15,
        "F_y = {}, expected {expected}",
        force[1]
    );
    assert!(
        force[0].abs() < 1e-30,
        "no horizontal force, F_x = {}",
        force[0]
    );
}

#[test]
fn cold_fluid_is_pushed_along_the_gravity_vector() {
    let coupling = BoussinesqCoupling::new(GRAVITY_DOWN, 2.0e-3, 1.0);
    let force = coupling.buoyancy_force(0.9, 0.25);
    let expected = -0.9 * 2.0e-3 * 0.75 * 1.0e-2;
    assert!(force[1] < 0.0, "cold fluid must sink, F_y = {}", force[1]);
    assert!(
        (force[1] - expected).abs() < 1e-15,
        "F_y = {}, expected {expected}",
        force[1]
    );
}

#[test]
fn both_gravity_components_are_reversed_for_hot_fluid() {
    let coupling = BoussinesqCoupling::new([2.0e-2, -1.0e-2], 1.0e-3, 0.0);
    let force = coupling.buoyancy_force(1.0, 3.0);
    assert!(force[0] < 0.0 && force[1] > 0.0, "F = {force:?}");
    assert!((force[0] + 1.0 * 1.0e-3 * 3.0 * 2.0e-2).abs() < 1e-15);
    assert!((force[1] - 1.0 * 1.0e-3 * 3.0 * 1.0e-2).abs() < 1e-15);
}

#[test]
fn force_scales_linearly_with_density_expansion_and_temperature_difference() {
    let base = BoussinesqCoupling::new(GRAVITY_DOWN, 1.0e-3, 1.0).buoyancy_force(1.0, 2.0);
    let double_rho = BoussinesqCoupling::new(GRAVITY_DOWN, 1.0e-3, 1.0).buoyancy_force(2.0, 2.0);
    let double_beta = BoussinesqCoupling::new(GRAVITY_DOWN, 2.0e-3, 1.0).buoyancy_force(1.0, 2.0);
    let double_delta_t =
        BoussinesqCoupling::new(GRAVITY_DOWN, 1.0e-3, 1.0).buoyancy_force(1.0, 3.0);
    assert!(base[1] > 0.0);
    assert!((double_rho[1] / base[1] - 2.0).abs() < 1e-12);
    assert!((double_beta[1] / base[1] - 2.0).abs() < 1e-12);
    assert!((double_delta_t[1] / base[1] - 2.0).abs() < 1e-12);
}

#[test]
fn guo_buoyancy_adds_upward_momentum_to_hot_rest_fluid() {
    // Rest fluid at a node hotter than T_ref: the Guo source term carries
    // (1 - 1/(2 tau)) F of momentum and no mass.
    let tau = 0.8;
    let coupling = BoussinesqCoupling::new(GRAVITY_DOWN, 2.0e-3, 1.0);
    let force = coupling.buoyancy_force(1.0, 1.5);
    let mut f = D2Q9_WEIGHTS;
    coupling.apply_to_distribution(&mut f, 1.0, [0.0, 0.0], 1.5, tau);

    let (mut mass, mut momentum_x, mut momentum_y) = (0.0, 0.0, 0.0);
    for i in 0..9 {
        let delta = f[i] - D2Q9_WEIGHTS[i];
        mass += delta;
        momentum_x += delta * f64::from(D2Q9_VELOCITIES[i][0]);
        momentum_y += delta * f64::from(D2Q9_VELOCITIES[i][1]);
    }
    let prefactor = 1.0 - 1.0 / (2.0 * tau);
    assert!(mass.abs() < 1e-15, "the forcing must not add mass: {mass}");
    assert!(
        momentum_x.abs() < 1e-15,
        "no horizontal momentum: {momentum_x}"
    );
    assert!(momentum_y > 0.0, "hot fluid must gain upward momentum");
    assert!(
        (momentum_y - prefactor * force[1]).abs() < 1e-13,
        "momentum {momentum_y} vs (1 - 1/2tau) F_y = {}",
        prefactor * force[1]
    );
}

#[test]
fn guo_buoyancy_adds_downward_momentum_to_cold_rest_fluid() {
    let tau = 1.2;
    let coupling = BoussinesqCoupling::new(GRAVITY_DOWN, 2.0e-3, 1.0);
    let mut f = D2Q9_WEIGHTS;
    coupling.apply_to_distribution(&mut f, 1.0, [0.0, 0.0], 0.4, tau);
    let momentum_y: f64 = (0..9)
        .map(|i| (f[i] - D2Q9_WEIGHTS[i]) * f64::from(D2Q9_VELOCITIES[i][1]))
        .sum();
    assert!(momentum_y < 0.0, "cold fluid must gain downward momentum");
}

#[test]
fn thermal_lbm_buoyancy_field_points_up_in_the_hot_spot() {
    let (nx, ny) = (6, 5);
    let mut lbm = ThermalLbm::new(nx, ny, 1.0, 1.0).with_boussinesq(GRAVITY_DOWN, 2.0e-3, 1.0);
    let hot = lbm.thermal.idx(2, 3);
    let cold = lbm.thermal.idx(4, 1);
    lbm.thermal.temperature[hot] = 2.0;
    lbm.thermal.temperature[cold] = 0.0;
    let densities = vec![1.0; nx * ny];
    let forces = lbm
        .compute_buoyancy_forces(&densities)
        .expect("a Boussinesq coupling is configured");
    assert!(forces[hot][1] > 0.0, "hot spot: {:?}", forces[hot]);
    assert!(forces[cold][1] < 0.0, "cold spot: {:?}", forces[cold]);
    let neutral = lbm.thermal.idx(0, 0);
    assert!(forces[neutral][1].abs() < 1e-30);
}
