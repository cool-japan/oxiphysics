//! Regression tests for the fluid–solid coupling of the conjugate heat solver.
//!
//! The previous solver ran the lattice on every node (solid ones included), let the solid
//! finite-difference stencil read the frozen initial temperature of its fluid neighbours and then
//! overwrote both interface nodes with the same interface temperature, so the interface flux
//! `k_f (T_s − T_f)/Δx` — and with it the Nusselt number and the convection coefficient — was
//! identically zero after every step. These tests pin the masked lattice, the half-cell
//! interface treatment and the steady state of a composite slab.

use super::*;

/// Fluid conductivity of the slab tests: `ρ c_p α_f` with `ρ = c_p = 1`, `τ_t = 1`.
const K_FLUID: f64 = 1.0 / 6.0;
/// Solid conductivity of the slab tests (`α_s = 0.2` with `ρ c_p = 2.5`).
const K_SOLID: f64 = 0.5;
const T_HOT: f64 = 400.0;
const T_COLD: f64 = 300.0;

/// Solid material of the slab tests: `α_s = k/(ρ c_p) = 0.2`.
fn wall_material() -> SolidRegion {
    SolidRegion::new(K_SOLID, 2.5, 1.0, "wall")
}

/// Cell index helper.
fn at(nx: usize, x: usize, y: usize) -> usize {
    y * nx + x
}

/// Steady piecewise-linear temperature of the composite slab at node column `x`.
///
/// Solid for `x ≤ 9`, fluid for `10 ≤ x ≤ 29`, solid for `x ≥ 30`; the interfaces sit at
/// `x = 9.5` and `x = 29.5`, the end nodes `x = 0` and `x = 39` are held at the hot and cold
/// temperature.
fn slab_temperature(x: usize) -> f64 {
    let flux = slab_flux();
    let t_left = T_HOT - flux * 9.5 / K_SOLID;
    let t_right = t_left - flux * 20.0 / K_FLUID;
    let x = x as f64;
    if x < 9.5 {
        T_HOT - flux * x / K_SOLID
    } else if x < 29.5 {
        t_left - flux * (x - 9.5) / K_FLUID
    } else {
        t_right - flux * (x - 29.5) / K_SOLID
    }
}

/// Steady heat flux through the slab (series of the three layers).
fn slab_flux() -> f64 {
    (T_HOT - T_COLD) / (2.0 * 9.5 / K_SOLID + 20.0 / K_FLUID)
}

/// Solid–fluid–solid slab, uniform in y: hot solid wall (`x ≤ 9`), fluid (`10 ≤ x ≤ 29`),
/// cold solid wall (`x ≥ 30`).
///
/// The outer ring of the grid is held fixed by the solver; its solid nodes are therefore
/// initialised with the steady profile (uniform in y), so that the 1-D steady state is
/// compatible with the ring. Everything else starts away from the steady state.
fn composite_slab() -> ConjugateHeatSolver {
    let (nx, ny) = (40, 3);
    let is_fluid: Vec<bool> = (0..nx * ny)
        .map(|k| (10..=29).contains(&(k % nx)))
        .collect();
    let mut solver = ConjugateHeatSolver::new(
        nx,
        ny,
        1.0,
        1.0,
        1.0,
        0.5 * (T_HOT + T_COLD),
        is_fluid,
        wall_material(),
        1.0,
        1.0,
    )
    .unwrap();
    for y in 0..ny {
        for x in 0..nx {
            let k = at(nx, x, y);
            let on_ring = x == 0 || x == nx - 1 || y == 0 || y == ny - 1;
            if !solver.is_fluid[k] {
                solver.solid_temp[k] = if on_ring {
                    slab_temperature(x)
                } else if x < 10 {
                    T_HOT
                } else {
                    T_COLD
                };
            }
        }
        solver
            .add_interface(FluidSolidInterface::new(
                at(nx, 10, y),
                at(nx, 9, y),
                [-1.0, 0.0],
                K_FLUID,
                K_SOLID,
            ))
            .unwrap();
        solver
            .add_interface(FluidSolidInterface::new(
                at(nx, 29, y),
                at(nx, 30, y),
                [1.0, 0.0],
                K_FLUID,
                K_SOLID,
            ))
            .unwrap();
    }
    solver
}

#[test]
fn composite_slab_converges_to_the_piecewise_linear_steady_state() {
    let mut solver = composite_slab();
    for _ in 0..12_000 {
        solver.step();
    }
    let nx = solver.lbm.nx;
    let field = solver.temperature_field();
    let mid_row = nx;
    // The half-way anti-bounce-back link and the half-cell ghost value reproduce the
    // piecewise-linear solution of the series resistance exactly at steady state.
    for x in 0..nx {
        let analytic = slab_temperature(x);
        assert!(
            (field[mid_row + x] - analytic).abs() < 1e-6,
            "x = {x}: T = {}, analytic {analytic}",
            field[mid_row + x]
        );
    }

    // Flux continuity k_f ∇T_f = k_s ∇T_s and a non-zero flux equal to the series-resistance one.
    let flux = slab_flux();
    let grad_fluid = (field[mid_row + 20] - field[mid_row + 15]) / 5.0;
    let grad_solid = (field[mid_row + 6] - field[mid_row + 3]) / 3.0;
    assert!(
        (K_FLUID * grad_fluid + flux).abs() < 1e-7 * flux,
        "fluid conduction flux {} vs {flux}",
        -K_FLUID * grad_fluid
    );
    assert!(
        (K_SOLID * grad_solid + flux).abs() < 1e-7 * flux,
        "solid conduction flux {} vs {flux}",
        -K_SOLID * grad_solid
    );

    // Interface fluxes are positive where heat flows from the solid into the fluid: the hot
    // wall feeds the fluid (+q) and the fluid feeds the cold wall (−q).
    let result = solver.results(T_COLD, T_HOT, K_FLUID);
    assert_eq!(result.interface_heat_flux.len(), 6);
    for (pair, q) in result.interface_heat_flux.chunks(2).enumerate() {
        assert!(
            (q[0] - flux).abs() < 1e-7 * flux,
            "row {pair}: hot-side interface flux {} vs {flux}",
            q[0]
        );
        assert!(
            (q[1] + flux).abs() < 1e-7 * flux,
            "row {pair}: cold-side interface flux {} vs {}",
            q[1],
            -flux
        );
    }
    assert!(result.convective_coefficient > 0.0);
    assert!(result.nusselt_number > 0.0);
    // The interface temperatures follow the analytic steady state.
    let t_left = T_HOT - flux * 9.5 / K_SOLID;
    let t_right = t_left - flux * 20.0 / K_FLUID;
    for (index, iface) in solver.interfaces.iter().enumerate() {
        let expected = if index % 2 == 0 { t_left } else { t_right };
        assert!(
            (iface.interface_temperature - expected).abs() < 1e-6,
            "interface {index}: T_i = {}, analytic {expected}",
            iface.interface_temperature
        );
    }
}

/// Hot fluid (x ≤ 2, 400 K) next to a cold solid wall (x ≥ 3, 300 K) on a 7 × 7 grid with a
/// single interface at the node pair in row 3.
fn hot_fluid_cold_wall() -> ConjugateHeatSolver {
    let (nx, ny) = (7, 7);
    let is_fluid: Vec<bool> = (0..nx * ny).map(|k| k % nx <= 2).collect();
    let mut solver = ConjugateHeatSolver::new(
        nx,
        ny,
        1.0,
        1.0,
        1.0,
        400.0,
        is_fluid,
        wall_material(),
        1.0,
        1.0,
    )
    .unwrap();
    for k in 0..nx * ny {
        if !solver.is_fluid[k] {
            solver.solid_temp[k] = 300.0;
        }
    }
    solver
        .add_interface(FluidSolidInterface::new(
            at(nx, 2, 3),
            at(nx, 3, 3),
            [1.0, 0.0],
            K_FLUID,
            K_SOLID,
        ))
        .unwrap();
    solver
}

#[test]
fn interface_flux_and_nusselt_number_are_non_zero_after_a_step() {
    let mut solver = hot_fluid_cold_wall();
    solver.step();
    let result = solver.results(400.0, 300.0, K_FLUID);
    assert_eq!(result.interface_heat_flux.len(), 1);
    // Heat flows from the hot fluid into the cold solid: negative flux (solid → fluid positive).
    let q = result.interface_heat_flux[0];
    assert!(q < -1.0, "interface flux {q}");
    let nx = solver.lbm.nx;
    let t_fluid = solver.lbm.temperature[at(nx, 2, 3)];
    let t_solid = solver.solid_temp[at(nx, 3, 3)];
    let expected = 2.0 * K_FLUID * K_SOLID * (t_solid - t_fluid) / ((K_FLUID + K_SOLID) * 1.0);
    assert!(
        (q - expected).abs() < 1e-12 * expected.abs(),
        "{q} vs {expected}"
    );
    assert!((result.convective_coefficient - q.abs() / 100.0).abs() < 1e-12);
    assert!(
        (result.nusselt_number - result.convective_coefficient * 7.0 / K_FLUID).abs() < 1e-9,
        "Nu = h L / k"
    );
    assert!(result.convective_coefficient > 0.0 && result.nusselt_number > 0.0);
}

#[test]
fn interface_nodes_keep_their_own_temperatures() {
    let mut solver = composite_slab();
    for _ in 0..200 {
        solver.step();
    }
    let nx = solver.lbm.nx;
    // Across the left interface the node temperatures differ by q/2 (1/k_f + 1/k_s) > 0 in a
    // developed state; the old coupling overwrote both with the same value.
    let t_solid = solver.solid_temp[at(nx, 9, 1)];
    let t_fluid = solver.lbm.temperature[at(nx, 10, 1)];
    assert!(
        (t_solid - t_fluid).abs() > 1e-3,
        "node temperatures coincide: {t_solid} vs {t_fluid}"
    );
}

#[test]
fn solid_stencil_sees_the_half_cell_interface_ghost_value() {
    // Fluid for x <= 2 at 400 K, solid for x >= 3 at 300 K; only the node pair in row 3 carries an
    // interface (k_f = 1/6, k_s = 1/2 → T_i = 325 K), the other fluid–solid links are adiabatic.
    let mut solver = hot_fluid_cold_wall();
    let nx = solver.lbm.nx;
    solver.step();
    let r = solver.solid_diffusion_number();
    assert!((r - 0.2).abs() < 1e-15);
    // T' = T + 2 r (T_i − T) through the half-cell with the interface …
    let expected = 300.0 + 2.0 * r * (325.0 - 300.0);
    assert!(
        (solver.solid_temp[at(nx, 3, 3)] - expected).abs() < 1e-12,
        "interface-adjacent solid node: {} vs {expected}",
        solver.solid_temp[at(nx, 3, 3)]
    );
    // … while a solid node whose fluid neighbour carries no interface sees an adiabatic link.
    for y in [2, 4] {
        assert!(
            (solver.solid_temp[at(nx, 3, y)] - 300.0).abs() < 1e-12,
            "adiabatic link at row {y}: {}",
            solver.solid_temp[at(nx, 3, y)]
        );
    }
    // The node temperatures were not overwritten by the interface temperature.
    assert!((solver.solid_temp[at(nx, 3, 3)] - 325.0).abs() > 1.0);
    assert!((solver.lbm.temperature[at(nx, 2, 3)] - 325.0).abs() > 1.0);
    // The merged temperature field reports the solid nodes from the finite-difference grid.
    let field = solver.temperature_field();
    assert!((field[at(nx, 3, 3)] - expected).abs() < 1e-12);
    assert!((field[at(nx, 5, 5)] - 300.0).abs() < 1e-12);
    assert!((field[at(nx, 0, 0)] - 400.0).abs() < 1e-9);
}

#[test]
fn masked_lattice_never_touches_solid_nodes_and_conserves_the_fluid_totals() {
    let (nx, ny) = (12, 8);
    let mut dd = DoubleDistribution::new(nx, ny, 0.8, 0.9, 1.0, 300.0);
    let is_fluid: Vec<bool> = (0..nx * ny)
        .map(|k| !((5..=7).contains(&(k % nx)) && (2..=4).contains(&(k / nx))))
        .collect();
    // Moving fluid with a hot spot; the solid nodes hold recognisable marker populations.
    let u0 = [0.05, 0.02];
    for (k, &fluid) in is_fluid.iter().enumerate() {
        for i in 0..9 {
            dd.f[k][i] = dd.f_equilibrium(1.0, u0, i);
        }
        for i in 0..5 {
            dd.h[k][i] = dd.h_equilibrium(300.0, u0, i);
        }
        dd.u[k] = u0;
        if !fluid {
            for i in 0..9 {
                dd.f[k][i] = 10.0 + i as f64;
            }
            for i in 0..5 {
                dd.h[k][i] = 20.0 + i as f64;
            }
        }
    }
    let hot = at(nx, 2, 3);
    for i in 0..5 {
        dd.h[hot][i] = dd.h_equilibrium(900.0, u0, i);
    }
    dd.compute_macroscopic_masked(&is_fluid);
    let solid_f: Vec<[f64; 9]> = dd.f.clone();
    let solid_h: Vec<[f64; 5]> = dd.h.clone();
    let fluid_sum = |dd: &DoubleDistribution| -> (f64, f64) {
        let mut rho = 0.0;
        let mut temperature = 0.0;
        for (k, &fluid) in is_fluid.iter().enumerate() {
            if fluid {
                rho += dd.f[k].iter().sum::<f64>();
                temperature += dd.h[k].iter().sum::<f64>();
            }
        }
        (rho, temperature)
    };
    let (rho0, temperature0) = fluid_sum(&dd);
    let no_walls = vec![[None; 5]; nx * ny];
    for _ in 0..60 {
        dd.collide_masked(&is_fluid);
        dd.stream_masked(&is_fluid, &no_walls);
        dd.compute_macroscopic_masked(&is_fluid);
    }
    for k in 0..nx * ny {
        if !is_fluid[k] {
            assert_eq!(dd.f[k], solid_f[k], "solid node {k}: f changed");
            assert_eq!(dd.h[k], solid_h[k], "solid node {k}: h changed");
        }
    }
    let (rho1, temperature1) = fluid_sum(&dd);
    assert!(
        (rho1 - rho0).abs() < 1e-9 * rho0,
        "fluid mass: {rho0} -> {rho1}"
    );
    assert!(
        (temperature1 - temperature0).abs() < 1e-9 * temperature0,
        "adiabatic walls must conserve the fluid heat: {temperature0} -> {temperature1}"
    );
    assert!(dd.f.iter().all(|cell| cell.iter().all(|v| v.is_finite())));
}

#[test]
fn bounce_back_reflects_the_flow_at_a_solid_link() {
    // A fluid node whose east neighbour is solid: the population leaving east comes back west.
    let (nx, ny) = (4, 3);
    let mut dd = DoubleDistribution::new(nx, ny, 1.0, 1.0, 1.0, 300.0);
    let is_fluid: Vec<bool> = (0..nx * ny).map(|k| k % nx != 2).collect();
    let k = at(nx, 1, 1);
    dd.f[k][1] = 0.7; // east-moving population
    dd.stream_masked(&is_fluid, &vec![[None; 5]; nx * ny]);
    assert!(
        (dd.f[k][3] - 0.7).abs() < 1e-15,
        "reflected west-moving population {}",
        dd.f[k][3]
    );
    assert!(
        (dd.f[at(nx, 2, 1)][1] - 1.0 / 9.0).abs() < 1e-15,
        "the solid node must keep its population"
    );
}

#[test]
fn dirichlet_wall_link_heats_a_closed_fluid_column_to_the_wall_temperature() {
    // Closed box: fluid x = 0..=4, solid column x = 5 (periodic wrap makes it the west
    // neighbour of x = 0 too). Only the east link of x = 4 is a Dirichlet wall at 400 K, all
    // others are adiabatic, so the box relaxes to a uniform 400 K.
    let (nx, ny) = (6, 3);
    let mut dd = DoubleDistribution::new(nx, ny, 1.0, 1.0, 1.0, 300.0);
    let is_fluid: Vec<bool> = (0..nx * ny).map(|k| k % nx != 5).collect();
    let mut walls = vec![[None; 5]; nx * ny];
    for y in 0..ny {
        walls[at(nx, 4, y)][1] = Some(400.0);
    }
    let mean = |dd: &DoubleDistribution| -> f64 {
        let (mut sum, mut count) = (0.0, 0.0);
        for (k, &fluid) in is_fluid.iter().enumerate() {
            if fluid {
                sum += dd.temperature[k];
                count += 1.0;
            }
        }
        sum / count
    };
    let mut previous = mean(&dd);
    for step in 0..6000 {
        dd.collide_masked(&is_fluid);
        dd.stream_masked(&is_fluid, &walls);
        dd.compute_macroscopic_masked(&is_fluid);
        if step % 500 == 0 {
            let now = mean(&dd);
            assert!(now >= previous - 1e-12, "step {step}: {previous} -> {now}");
            previous = now;
        }
    }
    for (k, &fluid) in is_fluid.iter().enumerate() {
        if fluid {
            assert!(
                (dd.temperature[k] - 400.0).abs() < 1e-3,
                "node {k}: T = {}",
                dd.temperature[k]
            );
        }
    }
}

#[test]
fn constructor_validates_grid_relaxation_and_diffusion_number() {
    let make = |tau: f64, tau_t: f64, mask: usize, solid: SolidRegion, dx: f64, dt: f64| {
        ConjugateHeatSolver::new(
            4,
            4,
            tau,
            tau_t,
            1.0,
            300.0,
            vec![true; mask],
            solid,
            dx,
            dt,
        )
    };
    assert!(make(1.0, 1.0, 16, wall_material(), 1.0, 1.0).is_ok());
    assert!(matches!(
        make(1.0, 1.0, 15, wall_material(), 1.0, 1.0),
        Err(Error::ComponentCountMismatch {
            expected: 16,
            found: 15,
            ..
        })
    ));
    assert!(matches!(
        make(0.5, 1.0, 16, wall_material(), 1.0, 1.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    assert!(matches!(
        make(1.0, 0.4, 16, wall_material(), 1.0, 1.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    assert!(matches!(
        make(1.0, 1.0, 16, wall_material(), 0.0, 1.0),
        Err(Error::ParameterOutOfBounds { .. })
    ));
    // α_s dt/dx² = 0.3 > 1/4: the explicit stencil is unstable.
    let unstable = SolidRegion::new(0.75, 2.5, 1.0, "hot");
    match make(1.0, 1.0, 16, unstable, 1.0, 1.0) {
        Err(Error::ParameterOutOfBounds { value, max, .. }) => {
            assert!((value - 0.3).abs() < 1e-12 && (max - 0.25).abs() < 1e-15);
        }
        other => panic!("expected a stability error, got {:?}", other.map(|_| ())),
    }
    // Exactly at the limit is accepted.
    let marginal = SolidRegion::new(0.625, 2.5, 1.0, "limit");
    assert!(make(1.0, 1.0, 16, marginal, 1.0, 1.0).is_ok());
    // Halving dt makes the unstable case stable again.
    let was_unstable = SolidRegion::new(0.75, 2.5, 1.0, "hot");
    assert!(make(1.0, 1.0, 16, was_unstable, 1.0, 0.5).is_ok());
}

#[test]
fn add_interface_rejects_invalid_links() {
    let (nx, ny) = (6, 6);
    let is_fluid: Vec<bool> = (0..nx * ny).map(|k| k % nx < 3).collect();
    let mut solver = ConjugateHeatSolver::new(
        nx,
        ny,
        1.0,
        1.0,
        1.0,
        300.0,
        is_fluid,
        wall_material(),
        1.0,
        1.0,
    )
    .unwrap();
    let valid = FluidSolidInterface::new(at(nx, 2, 2), at(nx, 3, 2), [1.0, 0.0], K_FLUID, K_SOLID);
    assert!(solver.add_interface(valid).is_ok());
    assert_eq!(solver.interfaces.len(), 1);
    let cases = [
        // out of range
        FluidSolidInterface::new(100, at(nx, 3, 2), [1.0, 0.0], K_FLUID, K_SOLID),
        // fluid index is a solid node
        FluidSolidInterface::new(at(nx, 3, 2), at(nx, 4, 2), [1.0, 0.0], K_FLUID, K_SOLID),
        // solid index is a fluid node
        FluidSolidInterface::new(at(nx, 2, 2), at(nx, 1, 2), [-1.0, 0.0], K_FLUID, K_SOLID),
        // not neighbours
        FluidSolidInterface::new(at(nx, 2, 2), at(nx, 4, 2), [1.0, 0.0], K_FLUID, K_SOLID),
        // diagonal neighbours are not axis-adjacent
        FluidSolidInterface::new(at(nx, 2, 2), at(nx, 3, 3), [1.0, 0.0], K_FLUID, K_SOLID),
        // non-positive conductivities
        FluidSolidInterface::new(at(nx, 2, 3), at(nx, 3, 3), [1.0, 0.0], 0.0, K_SOLID),
        FluidSolidInterface::new(at(nx, 2, 3), at(nx, 3, 3), [1.0, 0.0], K_FLUID, -1.0),
    ];
    for (index, iface) in cases.into_iter().enumerate() {
        assert!(
            solver.add_interface(iface).is_err(),
            "case {index} must be rejected"
        );
    }
    assert_eq!(
        solver.interfaces.len(),
        1,
        "rejected interfaces are not stored"
    );
}

#[test]
fn add_interface_checks_the_conductivity_contrast_against_the_stencil_weight() {
    // r = 0.24: the bulk limit 1/4 holds, but a highly conducting fluid (beta → −1) leaves the
    // solid node with the weight 1 − r (4 − beta) ≈ 1 − 0.24 · 5 < 0.
    let (nx, ny) = (6, 6);
    let is_fluid: Vec<bool> = (0..nx * ny).map(|k| k % nx < 3).collect();
    let solid = SolidRegion::new(0.6, 2.5, 1.0, "fast");
    let mut solver =
        ConjugateHeatSolver::new(nx, ny, 1.0, 1.0, 1.0, 300.0, is_fluid, solid, 1.0, 1.0).unwrap();
    assert!((solver.solid_diffusion_number() - 0.24).abs() < 1e-15);
    let strongly_conducting_fluid =
        FluidSolidInterface::new(at(nx, 2, 2), at(nx, 3, 2), [1.0, 0.0], 1000.0, 0.6);
    assert!(solver.add_interface(strongly_conducting_fluid).is_err());
    let comparable = FluidSolidInterface::new(at(nx, 2, 2), at(nx, 3, 2), [1.0, 0.0], 0.6, 0.6);
    assert!(solver.add_interface(comparable).is_ok());
}

#[test]
fn interface_heat_flux_matches_the_series_resistance_formula() {
    let iface = FluidSolidInterface::new(0, 1, [1.0, 0.0], 0.6, 16.0);
    let (t_fluid, t_solid, dx) = (300.0, 340.0, 0.5);
    let q = iface.heat_flux(t_fluid, t_solid, dx);
    // Series resistance of the two half-cells, each of length dx/2.
    let resistance = (dx / 2.0) / 0.6 + (dx / 2.0) / 16.0;
    assert!((q - (t_solid - t_fluid) / resistance).abs() < 1e-9 * q.abs());
    // Both half-cells carry that flux through the interface temperature.
    let mut probe = iface.clone();
    probe.compute_interface_temperature(t_fluid, t_solid);
    let from_fluid = 0.6 * (probe.interface_temperature - t_fluid) / (dx / 2.0);
    let from_solid = 16.0 * (t_solid - probe.interface_temperature) / (dx / 2.0);
    assert!((from_fluid - q).abs() < 1e-9 * q.abs());
    assert!((from_solid - q).abs() < 1e-9 * q.abs());
    assert_eq!(iface.heat_flux(310.0, 310.0, dx), 0.0);
    assert!(iface.heat_flux(400.0, 300.0, dx) < 0.0);
}

#[test]
fn link_direction_finds_axis_neighbours_with_periodic_wrap() {
    let (nx, ny) = (5, 4);
    let origin = at(nx, 0, 0);
    assert_eq!(link_direction(nx, ny, origin, at(nx, 1, 0)), Some(1));
    assert_eq!(link_direction(nx, ny, origin, at(nx, 0, 1)), Some(2));
    assert_eq!(link_direction(nx, ny, origin, at(nx, 4, 0)), Some(3));
    assert_eq!(link_direction(nx, ny, origin, at(nx, 0, 3)), Some(4));
    assert_eq!(link_direction(nx, ny, origin, at(nx, 1, 1)), None);
    assert_eq!(link_direction(nx, ny, origin, at(nx, 2, 0)), None);
    assert_eq!(link_direction(nx, ny, origin, origin), None);
    assert_eq!(link_direction(nx, ny, origin, nx * ny), None);
}
