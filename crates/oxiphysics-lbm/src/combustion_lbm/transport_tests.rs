//! Regression tests for the advection-diffusion transport of [`CombustionLbm`].
//!
//! The previous solver collided and streamed only the flow populations and then diffused the
//! species and the temperature with explicit finite differences that never saw the velocity
//! (and, with the physical `dt` multiplying the lattice diffusivity, hardly diffused at all), so
//! a flame could never be convected. These tests pin the transport: uniform flows carry blobs
//! of species and heat with the flow speed, the diffusivities follow the relaxation parameters,
//! and the totals are conserved.

use super::*;

/// Kinetics that never release anything.
fn inert_kinetics() -> ArrheniusKinetics {
    ArrheniusKinetics::new(0.0, 1.0, 8.314, 1.0, 1.0, 0.0, 1.0)
}

/// Cell index helper.
fn cell(nx: usize, ix: usize, iy: usize) -> usize {
    iy * nx + ix
}

/// Gaussian profile `exp(-(x - x0)^2 / (2 sigma^2))`.
fn gaussian(x: f64, x0: f64, sigma: f64) -> f64 {
    (-(x - x0) * (x - x0) / (2.0 * sigma * sigma)).exp()
}

/// Centroid `(x, y)` of the non-negative field `phi` on an `nx`-wide grid.
fn centroid(phi: &[f64], nx: usize) -> (f64, f64) {
    let (mut mass, mut mx, mut my) = (0.0, 0.0, 0.0);
    for (k, &value) in phi.iter().enumerate() {
        mass += value;
        mx += value * (k % nx) as f64;
        my += value * (k / nx) as f64;
    }
    (mx / mass, my / mass)
}

/// Variance of the x-distribution of `phi`.
fn variance_x(phi: &[f64], nx: usize) -> f64 {
    let (mean, _) = centroid(phi, nx);
    let (mut mass, mut m2) = (0.0, 0.0);
    for (k, &value) in phi.iter().enumerate() {
        mass += value;
        m2 += value * ((k % nx) as f64 - mean).powi(2);
    }
    m2 / mass
}

#[test]
fn uniform_flow_carries_a_species_blob_with_the_flow_speed() {
    let (nx, ny) = (64, 8);
    let mut solver = CombustionLbm::new(nx, ny, 1.0, 1.2, inert_kinetics(), 300.0, 0.01);
    solver.set_uniform_flow(1.0, 0.05, 0.0);
    for iy in 0..ny {
        for ix in 0..nx {
            solver.species.fuel[cell(nx, ix, iy)] = 0.2 * gaussian(ix as f64, 16.0, 3.0);
        }
    }
    let (x0, y0) = centroid(&solver.species.fuel, nx);
    let total0: f64 = solver.species.fuel.iter().sum();
    let steps = 200;
    solver.run(steps);
    let (x1, y1) = centroid(&solver.species.fuel, nx);
    let expected = 0.05 * steps as f64;
    assert!(
        (x1 - x0 - expected).abs() < 0.05,
        "blob moved {} cells in {steps} steps, expected {expected}",
        x1 - x0
    );
    assert!((y1 - y0).abs() < 1e-9, "no cross-stream drift: {}", y1 - y0);
    let total1: f64 = solver.species.fuel.iter().sum();
    assert!(
        (total1 - total0).abs() < 1e-9 * total0,
        "fuel not conserved: {total0} -> {total1}"
    );
    // A uniform field is invariant under advection.
    for &yo in &solver.species.oxidizer {
        assert!((yo - 0.23).abs() < 1e-12, "oxidizer {yo}");
    }
}

#[test]
fn hot_spot_is_advected_by_a_diagonal_flow() {
    // The grid is large enough that the diffusing blob does not wrap around the periodic
    // boundaries, which would bias the (non-periodic) centroid.
    let (nx, ny) = (96, 96);
    let t0 = 300.0;
    let mut solver = CombustionLbm::new(nx, ny, 1.2, 1.0, inert_kinetics(), t0, 0.01);
    solver.set_uniform_flow(1.0, 0.04, 0.03);
    for iy in 0..ny {
        for ix in 0..nx {
            let r2 = (ix as f64 - 24.0).powi(2) + (iy as f64 - 30.0).powi(2);
            solver.temperature.temp[cell(nx, ix, iy)] = t0 + 500.0 * (-r2 / 18.0).exp();
        }
    }
    let excess = |solver: &CombustionLbm| -> Vec<f64> {
        solver.temperature.temp.iter().map(|&t| t - t0).collect()
    };
    let (x0, y0) = centroid(&excess(&solver), nx);
    let total0: f64 = excess(&solver).iter().sum();
    let steps = 200;
    solver.run(steps);
    let (x1, y1) = centroid(&excess(&solver), nx);
    assert!(
        (x1 - x0 - 0.04 * steps as f64).abs() < 0.1,
        "x shift {}, expected {}",
        x1 - x0,
        0.04 * steps as f64
    );
    assert!(
        (y1 - y0 - 0.03 * steps as f64).abs() < 0.1,
        "y shift {}, expected {}",
        y1 - y0,
        0.03 * steps as f64
    );
    let total1: f64 = excess(&solver).iter().sum();
    assert!(
        (total1 - total0).abs() < 1e-9 * total0,
        "heat not conserved: {total0} -> {total1}"
    );
}

#[test]
fn species_and_temperature_diffuse_with_their_own_relaxation_parameters() {
    let (nx, ny) = (96, 4);
    let (omega, omega_s) = (1.4, 1.2);
    let mut solver = CombustionLbm::new(nx, ny, omega, omega_s, inert_kinetics(), 300.0, 0.01);
    let t0 = 300.0;
    for iy in 0..ny {
        for ix in 0..nx {
            let g = gaussian(ix as f64, 48.0, 3.0);
            solver.species.fuel[cell(nx, ix, iy)] = 0.2 * g;
            solver.temperature.temp[cell(nx, ix, iy)] = t0 + 400.0 * g;
        }
    }
    let temperature_excess = |solver: &CombustionLbm| -> Vec<f64> {
        solver.temperature.temp.iter().map(|&t| t - t0).collect()
    };
    let var_fuel0 = variance_x(&solver.species.fuel, nx);
    let var_temp0 = variance_x(&temperature_excess(&solver), nx);
    let steps = 150;
    solver.run(steps);
    let d_species = CS2 * (1.0 / omega_s - 0.5);
    let d_thermal = CS2 * (1.0 / omega - 0.5);
    let growth_fuel = variance_x(&solver.species.fuel, nx) - var_fuel0;
    let growth_temp = variance_x(&temperature_excess(&solver), nx) - var_temp0;
    let expected_fuel = 2.0 * d_species * steps as f64;
    let expected_temp = 2.0 * d_thermal * steps as f64;
    assert!(
        (growth_fuel - expected_fuel).abs() < 0.02 * expected_fuel,
        "fuel variance grew by {growth_fuel}, expected 2 D n = {expected_fuel}"
    );
    assert!(
        (growth_temp - expected_temp).abs() < 0.02 * expected_temp,
        "temperature variance grew by {growth_temp}, expected 2 alpha n = {expected_temp}"
    );
    // No flow: the blobs stay where they were.
    let (cx, _) = centroid(&solver.species.fuel, nx);
    assert!((cx - 48.0).abs() < 1e-9, "centroid drifted to {cx}");
}

#[test]
fn reacting_flow_conserves_the_total_mass_fraction() {
    // Ta = Ea / R = 1000 K: a fast one-step reaction at 1500 K.
    let kinetics = ArrheniusKinetics::new(80.0, 8314.0, 8.314, 1.0, 1.0, 2.0e6, 4.0);
    let (nx, ny) = (32, 8);
    let mut solver = CombustionLbm::new(nx, ny, 1.0, 1.1, kinetics, 300.0, 1.0e-3);
    solver.set_uniform_flow(1.0, 0.05, 0.0);
    for iy in 0..ny {
        for ix in 8..14 {
            solver.temperature.temp[cell(nx, ix, iy)] = 1500.0;
        }
    }
    let mass = |solver: &CombustionLbm| -> f64 {
        solver
            .species
            .fuel
            .iter()
            .chain(solver.species.oxidizer.iter())
            .chain(solver.species.products.iter())
            .sum()
    };
    let mass0 = mass(&solver);
    let fuel0: f64 = solver.species.fuel.iter().sum();
    solver.run(120);
    let fuel1: f64 = solver.species.fuel.iter().sum();
    let products: f64 = solver.species.products.iter().sum();
    assert!(fuel1 < fuel0 - 1e-3, "no fuel consumed: {fuel0} -> {fuel1}");
    assert!(products > 1e-3, "no products formed: {products}");
    assert!(
        (mass(&solver) - mass0).abs() < 1e-9 * mass0,
        "total mass fraction not conserved: {mass0} -> {}",
        mass(&solver)
    );
}

#[test]
fn products_are_convected_downstream_of_the_reaction_zone() {
    let kinetics = ArrheniusKinetics::new(80.0, 8314.0, 8.314, 1.0, 1.0, 0.0, 4.0);
    let (nx, ny) = (48, 4);
    let mut solver = CombustionLbm::new(nx, ny, 1.0, 1.1, kinetics, 300.0, 1.0e-3);
    solver.set_uniform_flow(1.0, 0.06, 0.0);
    for iy in 0..ny {
        for ix in 6..10 {
            solver.temperature.temp[cell(nx, ix, iy)] = 1500.0;
        }
    }
    solver.run(100);
    let (cx, _) = centroid(&solver.species.products, nx);
    let zone_centre = 7.5;
    assert!(
        cx > zone_centre + 3.0,
        "products must be carried downstream of the reaction zone, centroid at {cx}"
    );
}

#[test]
fn ignition_heat_is_imported_conserved_and_advected() {
    // A symmetric, large periodic grid: the blob spreads to a standard deviation of about 6
    // cells and must not wrap around, so the centroid is meaningful.
    let (nx, ny) = (64, 41);
    let t0 = 300.0;
    let mut solver = CombustionLbm::new(nx, ny, 1.0, 1.0, inert_kinetics(), t0, 0.01);
    solver.set_uniform_flow(1.0, 0.05, 0.0);
    solver.ignite(20, 20, 800.0);
    let excess = |solver: &CombustionLbm| -> Vec<f64> {
        solver.temperature.temp.iter().map(|&t| t - t0).collect()
    };
    solver.run(100);
    let field = excess(&solver);
    let total: f64 = field.iter().sum();
    assert!(
        (total - 800.0).abs() < 1e-8,
        "the ignition heat must be conserved by the transport: {total}"
    );
    let (cx, cy) = centroid(&field, nx);
    assert!(
        (cx - 20.0 - 5.0).abs() < 0.2,
        "hot spot moved {} cells in 100 steps, expected 5",
        cx - 20.0
    );
    assert!((cy - 20.0).abs() < 1e-6, "cross-stream drift: {cy}");
}

#[test]
fn flame_front_moves_with_the_flow() {
    let (nx, ny) = (60, 4);
    let t0 = 300.0;
    let run_front = |ux: f64| -> f64 {
        let mut solver = CombustionLbm::new(nx, ny, 1.0, 1.0, inert_kinetics(), t0, 0.01);
        solver.set_uniform_flow(1.0, ux, 0.0);
        for iy in 0..ny {
            for ix in 10..16 {
                solver.temperature.temp[cell(nx, ix, iy)] = t0 + 1000.0;
            }
        }
        solver.run(100);
        solver.mean_flame_x()
    };
    let moved = run_front(0.05) - run_front(0.0);
    assert!(
        (moved - 5.0).abs() < 0.3,
        "the flame front moved {moved} cells relative to the quiescent case, expected 5"
    );
}

#[test]
fn scalar_equilibrium_carries_the_scalar_and_its_flux() {
    let (phi, ux, uy) = (0.37, 0.04, -0.03);
    let geq = scalar_equilibrium_d2q9(phi, ux, uy);
    let sum: f64 = geq.iter().sum();
    assert!((sum - phi).abs() < 1e-15, "zeroth moment {sum}");
    let (mut jx, mut jy) = (0.0, 0.0);
    for (q, &g) in geq.iter().enumerate() {
        jx += f64::from(C9[q].0) * g;
        jy += f64::from(C9[q].1) * g;
    }
    assert!((jx - phi * ux).abs() < 1e-15, "x flux {jx}");
    assert!((jy - phi * uy).abs() < 1e-15, "y flux {jy}");
}

#[test]
fn passive_scalar_relaxation_rate_follows_the_diffusivity() {
    for &d in &[0.01, 0.05, 0.1, 0.2] {
        let omega = scalar_relaxation_rate(d);
        assert!((omega - 1.0 / (3.0 * d + 0.5)).abs() < 1e-14);
        let solver = PassiveScalarLbm::uniform(4, 4, d, 1.0);
        assert!((solver.diffusivity() - d).abs() < 1e-14);
        assert!((solver.omega() - omega).abs() < 1e-15);
    }
    assert!((scalar_relaxation_rate(-0.3) - 2.0).abs() < 1e-15);
}

#[test]
fn passive_scalar_sync_distributes_changes_with_the_lattice_weights() {
    let (nx, ny) = (6, 5);
    let mut solver = PassiveScalarLbm::uniform(nx, ny, 0.1, 1.0);
    let before: Vec<[f64; 9]> = solver.populations().to_vec();
    let mut phi = vec![1.0; nx * ny];
    phi[cell(nx, 2, 3)] += 0.45;
    solver.sync_to_field(&phi);
    let k = cell(nx, 2, 3);
    for q in 0..9 {
        assert!(
            (solver.populations()[k][q] - before[k][q] - W9[q] * 0.45).abs() < 1e-15,
            "direction {q}"
        );
    }
    assert!((solver.total() - (nx * ny) as f64 - 0.45).abs() < 1e-12);
    assert!((solver.value(k).unwrap_or(f64::NAN) - 1.45).abs() < 1e-14);
    assert!(solver.value(nx * ny).is_none());
}

#[test]
fn passive_scalar_sync_with_its_own_field_is_the_identity() {
    let (nx, ny) = (8, 6);
    let mut solver = PassiveScalarLbm::uniform(nx, ny, 0.08, 0.0);
    let mut phi: Vec<f64> = (0..nx * ny).map(|k| 0.3 + 0.01 * (k % 7) as f64).collect();
    solver.sync_to_field(&phi);
    let ux = vec![0.03; nx * ny];
    let uy = vec![-0.02; nx * ny];
    for _ in 0..5 {
        solver.step(&ux, &uy);
    }
    let populations = solver.populations().to_vec();
    solver.read_field(&mut phi);
    solver.sync_to_field(&phi);
    for (cell_g, cell_prior) in solver.populations().iter().zip(populations.iter()) {
        for q in 0..9 {
            assert!(
                (cell_g[q] - cell_prior[q]).abs() < 1e-15,
                "the stored non-equilibrium part must be kept"
            );
        }
    }
}

#[test]
fn passive_scalar_from_field_rejects_a_wrong_length() {
    let result = PassiveScalarLbm::from_field(4, 3, 0.1, &[0.0; 5]);
    assert!(matches!(
        result,
        Err(crate::error::Error::ComponentCountMismatch {
            expected: 12,
            found: 5,
            ..
        })
    ));
    let ok = PassiveScalarLbm::from_field(2, 2, 0.1, &[1.0, 2.0, 3.0, 4.0]);
    assert!(ok.is_ok());
}
