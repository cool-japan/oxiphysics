// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the FE-driven [`MultiLoadTopOpt`] (whose old
//! implementation minimised a load proxy that *grows* with stiffness) and for
//! the casting (demoulding) filter (whose old implementation returned its
//! input unchanged).

use super::*;

/// A chain of independent springs: element `e` carries the force `f_e` of the
/// load case through a spring of stiffness `E(ρ_e)`.
///
/// Compliance is `Σ f_e² / E(ρ_e)` and the unit-stiffness element strain
/// energy `u_e² k0` with `u_e = f_e / E(ρ_e)` and `k0 = 1` is `f_e² / E(ρ_e)²`.
fn series_spring_fea(
    simp: SiMP,
) -> impl FnMut(&[f64], &LoadCase) -> PhysicsResult<(f64, Vec<f64>)> {
    move |density, load_case| {
        let mut compliance = 0.0;
        let mut strain_energy = Vec::with_capacity(density.len());
        for (&rho, &force) in density.iter().zip(load_case.loads.iter()) {
            let modulus = simp.interpolated_stiffness(rho);
            let displacement = force / modulus;
            compliance += force * displacement;
            strain_energy.push(displacement * displacement);
        }
        Ok((compliance, strain_energy))
    }
}

fn graded_loads(scale: f64) -> Vec<f64> {
    (0..16)
        .map(|e| scale * (1.0 + 0.25 * (e % 5) as f64))
        .collect()
}

fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

/// Exponent of the optimality-criteria update used for the series-spring
/// problems. Their strain energy scales like `ρ^-6` (the force through a spring
/// is fixed while its stiffness is `ρ^3`), so the update `ρ <- ρ B^η` with
/// `B ∝ ρ^-4` contracts only for `0 < η < 1/2`; at `η = 1/2` it merely flips
/// between two designs.
const SERIES_OC_DAMPING: f64 = 0.3;

fn series_problem(filter_radius: f64, load_scale: f64) -> (MultiLoadTopOpt, SiMP) {
    let grid = GridSize::new(4, 4);
    let simp = SiMP::new(grid, 0.4, 3.0, filter_radius, MaterialParams::default());
    let filter = FilteringMethods::new(filter_radius, grid);
    let mut optimiser = MultiLoadTopOpt::new(
        simp.clone(),
        filter,
        vec![LoadCase::new(graded_loads(load_scale), 1.0)],
    );
    optimiser.damping = SERIES_OC_DAMPING;
    (optimiser, simp)
}

/// Optimum of the series-spring problem: with `E = ρ^3` the stationarity
/// condition `3 f² / ρ^4 = λ` gives `ρ ∝ sqrt(f)`, scaled to the mean `volume`.
fn series_optimum(loads: &[f64], volume: f64) -> Vec<f64> {
    let roots: Vec<f64> = loads.iter().map(|f| f.sqrt()).collect();
    let scale = volume * roots.len() as f64 / roots.iter().sum::<f64>();
    roots.iter().map(|r| r * scale).collect()
}

fn invalid_field(error: &PhysicsError) -> &str {
    match error {
        PhysicsError::InvalidInput { field, .. } => field.as_str(),
        other => panic!("expected InvalidInput, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// MultiLoadTopOpt
// ---------------------------------------------------------------------------

#[test]
fn compliance_falls_as_the_structure_gets_stiffer() {
    let (optimiser, simp) = series_problem(1.0, 1.0);
    let (soft, soft_gradient) = optimiser
        .evaluate(&[0.4; 16], series_spring_fea(simp.clone()))
        .expect("the closure never fails");
    let (stiff, _) = optimiser
        .evaluate(&[0.8; 16], series_spring_fea(simp))
        .expect("the closure never fails");
    assert!(
        stiff < soft,
        "true compliance must drop with more material: {stiff} vs {soft}"
    );
    assert!(
        soft_gradient.iter().all(|&g| g < 0.0),
        "extra material can only lower the compliance: {soft_gradient:?}"
    );
}

#[test]
fn evaluate_returns_the_true_gradient_of_the_weighted_compliance() {
    let grid = GridSize::new(4, 4);
    let simp = SiMP::new(grid, 0.5, 3.0, 1.0, MaterialParams::default());
    let filter = FilteringMethods::new(1.0, grid);
    let load_cases = vec![
        LoadCase::new(graded_loads(1.0), 0.3),
        LoadCase::new(graded_loads(0.5).into_iter().rev().collect(), 0.7),
    ];
    let optimiser = MultiLoadTopOpt::new(simp.clone(), filter, load_cases);
    let density: Vec<f64> = (0..16).map(|e| 0.3 + 0.04 * e as f64).collect();

    let (compliance, gradient) = optimiser
        .evaluate(&density, series_spring_fea(simp.clone()))
        .expect("the closure never fails");
    assert!(compliance > 0.0);

    let step = 1e-6;
    for e in 0..16 {
        let mut plus = density.clone();
        plus[e] += step;
        let mut minus = density.clone();
        minus[e] -= step;
        let (c_plus, _) = optimiser
            .evaluate(&plus, series_spring_fea(simp.clone()))
            .expect("the closure never fails");
        let (c_minus, _) = optimiser
            .evaluate(&minus, series_spring_fea(simp.clone()))
            .expect("the closure never fails");
        let finite_difference = (c_plus - c_minus) / (2.0 * step);
        assert!(
            (finite_difference - gradient[e]).abs() <= 1e-5 * gradient[e].abs().max(1e-3),
            "element {e}: analytic {} vs finite difference {finite_difference}",
            gradient[e]
        );
    }
}

#[test]
fn equal_weight_split_matches_a_single_load_case() {
    let grid = GridSize::new(4, 4);
    let simp = SiMP::new(grid, 0.5, 3.0, 1.0, MaterialParams::default());
    let filter = FilteringMethods::new(1.0, grid);
    let split = MultiLoadTopOpt::new(
        simp.clone(),
        filter.clone(),
        vec![
            LoadCase::new(graded_loads(1.0), 0.5),
            LoadCase::new(graded_loads(1.0), 0.5),
        ],
    );
    let single = MultiLoadTopOpt::new(
        simp.clone(),
        filter,
        vec![LoadCase::new(graded_loads(1.0), 1.0)],
    );
    let density = vec![0.55; 16];
    let (c_split, g_split) = split
        .evaluate(&density, series_spring_fea(simp.clone()))
        .expect("the closure never fails");
    let (c_single, g_single) = single
        .evaluate(&density, series_spring_fea(simp))
        .expect("the closure never fails");
    assert!((c_split - c_single).abs() < 1e-10 * c_single.abs());
    for (a, b) in g_split.iter().zip(g_single.iter()) {
        assert!((a - b).abs() < 1e-10 * b.abs().max(1.0));
    }
}

#[test]
fn run_reduces_the_compliance_monotonically_and_meets_the_volume_target() {
    // A radius of 1 element makes the sensitivity filter the identity, so the
    // optimality-criteria update follows exact gradients.
    let (mut optimiser, simp) = series_problem(1.0, 1.0);
    let density = optimiser
        .run(80, 1e-9, series_spring_fea(simp))
        .expect("the closure never fails");

    let history = &optimiser.compliance_history;
    assert!(history.len() > 5, "history: {history:?}");
    for pair in history.windows(2) {
        assert!(
            pair[1] <= pair[0] * (1.0 + 1e-9),
            "compliance must not increase: {} -> {}",
            pair[0],
            pair[1]
        );
    }
    let (first, last) = (history[0], history[history.len() - 1]);
    assert!(
        last < 0.99 * first,
        "no real improvement: {first} -> {last}"
    );

    assert!(
        (mean(&density) - 0.4).abs() < 1e-3,
        "volume fraction {}",
        mean(&density)
    );
    assert!(density.iter().all(|&r| (0.0..=1.0).contains(&r)));
    assert_eq!(optimiser.simp.density, density);

    // It converges to the analytic optimum, rho proportional to sqrt(force).
    let optimum = series_optimum(&graded_loads(1.0), 0.4);
    for (e, (got, want)) in density.iter().zip(optimum.iter()).enumerate() {
        assert!((got - want).abs() < 1e-4, "element {e}: {got} vs {want}");
    }

    // The most heavily loaded springs get the most material.
    let loads = graded_loads(1.0);
    let heaviest = (0..16)
        .max_by(|&a, &b| loads[a].total_cmp(&loads[b]))
        .unwrap_or(0);
    let lightest = (0..16)
        .min_by(|&a, &b| loads[a].total_cmp(&loads[b]))
        .unwrap_or(0);
    assert!(
        density[heaviest] > density[lightest] + 0.05,
        "{} vs {}",
        density[heaviest],
        density[lightest]
    );
}

#[test]
fn run_is_independent_of_the_magnitude_of_the_loads() {
    // Scaling every load by 1e-6 scales the sensitivities by 1e-12; a fixed
    // linear bracket for the Lagrange multiplier cannot resolve that.
    let (mut reference, simp) = series_problem(1.0, 1.0);
    let expected = reference
        .run(30, 1e-12, series_spring_fea(simp))
        .expect("the closure never fails");
    let (mut tiny, simp) = series_problem(1.0, 1e-6);
    let density = tiny
        .run(30, 1e-12, series_spring_fea(simp))
        .expect("the closure never fails");
    assert!((mean(&density) - 0.4).abs() < 1e-3);
    for (a, b) in density.iter().zip(expected.iter()) {
        assert!((a - b).abs() < 1e-6, "{a} vs {b}");
    }
    // Both runs sit at the analytic optimum, not at a uniform fallback design.
    let optimum = series_optimum(&graded_loads(1.0), 0.4);
    for (a, b) in density.iter().zip(optimum.iter()) {
        assert!((a - b).abs() < 1e-4, "{a} vs {b}");
    }
}

#[test]
fn run_with_a_smoothing_filter_still_meets_the_volume_target() {
    let (mut optimiser, simp) = series_problem(1.5, 1.0);
    let density = optimiser
        .run(40, 1e-9, series_spring_fea(simp))
        .expect("the closure never fails");
    assert!((mean(&density) - 0.4).abs() < 1e-3);
    assert!(density.iter().all(|&r| (0.0..=1.0).contains(&r)));
}

#[test]
fn run_without_iterations_returns_the_initial_design() {
    let (mut optimiser, simp) = series_problem(1.0, 1.0);
    let density = optimiser
        .run(0, 1e-6, series_spring_fea(simp))
        .expect("nothing to evaluate");
    assert_eq!(density, vec![0.4; 16]);
    assert!(optimiser.compliance_history.is_empty());
}

#[test]
fn errors_of_the_fe_closure_are_propagated_and_leave_the_design_alone() {
    let (mut optimiser, _) = series_problem(1.0, 1.0);
    let failure = PhysicsError::NumericalDivergence {
        message: "singular stiffness matrix".to_string(),
    };
    let reported = failure.clone();
    let result = optimiser.run(10, 1e-6, move |_, _| Err(reported.clone()));
    assert_eq!(result, Err(failure));
    assert_eq!(optimiser.simp.density, vec![0.4; 16]);
}

#[test]
fn malformed_fe_results_and_setups_are_rejected() {
    let (optimiser, simp) = series_problem(1.0, 1.0);
    let density = vec![0.4; 16];

    let short = optimiser.evaluate(&density, |_, _| Ok((1.0, vec![0.0; 15])));
    assert_eq!(
        invalid_field(&short.expect_err("short energies")),
        "strain_energy"
    );

    let not_finite = optimiser.evaluate(&density, |_, _| Ok((f64::NAN, vec![0.0; 16])));
    assert_eq!(
        invalid_field(&not_finite.expect_err("NaN compliance")),
        "fea"
    );

    let infinite_energy = optimiser.evaluate(&density, |_, _| {
        let mut energy = vec![0.0; 16];
        energy[3] = f64::INFINITY;
        Ok((1.0, energy))
    });
    assert_eq!(
        invalid_field(&infinite_energy.expect_err("infinite energy")),
        "fea"
    );

    let wrong_density = optimiser.evaluate(&[0.4; 15], series_spring_fea(simp.clone()));
    assert_eq!(
        invalid_field(&wrong_density.expect_err("density length")),
        "density"
    );

    let mut no_loads = optimiser.clone();
    no_loads.load_cases.clear();
    let missing = no_loads.evaluate(&density, series_spring_fea(simp.clone()));
    assert_eq!(
        invalid_field(&missing.expect_err("no load case")),
        "load_cases"
    );

    let mut bad_weight = optimiser.clone();
    bad_weight.load_cases[0].weight = -1.0;
    let negative = bad_weight.evaluate(&density, series_spring_fea(simp.clone()));
    assert_eq!(
        invalid_field(&negative.expect_err("negative weight")),
        "load_cases[0].weight"
    );

    let mut bad_filter = optimiser.clone();
    bad_filter.filter = FilteringMethods::new(1.0, GridSize::new(3, 3));
    let mismatch = bad_filter.run(3, 1e-6, series_spring_fea(simp));
    assert_eq!(invalid_field(&mismatch.expect_err("filter grid")), "filter");
}

#[test]
fn run_records_the_compliance_reported_by_the_fe_closure() {
    // `run` must not consult any load proxy: with a closure that reports a
    // constant compliance the history is exactly that compliance.
    let (mut optimiser, _) = series_problem(1.0, 1.0);
    let _ = optimiser
        .run(3, 0.0, |_, _| Ok((7.5, vec![1.0; 16])))
        .expect("valid closure output");
    assert_eq!(optimiser.compliance_history, vec![7.5; 3]);
}

// ---------------------------------------------------------------------------
// Optimality-criteria multiplier search
// ---------------------------------------------------------------------------

#[test]
fn oc_update_meets_the_volume_target_for_tiny_sensitivities() {
    let grid = GridSize::new(4, 4);
    let simp = SiMP::new(grid, 0.4, 3.0, 1.5, MaterialParams::default());
    let filter = FilteringMethods::new(1.5, grid);
    let solver = TopOptSolver::new(simp, filter);
    for magnitude in [1e-12, 1e-9, 1e-3, 1.0, 1e6] {
        let sensitivity: Vec<f64> = (0..16).map(|e| -magnitude * (1.0 + e as f64)).collect();
        let (density, multiplier) = solver.oc_update(&[0.5; 16], &sensitivity, 0.4);
        assert!(
            (mean(&density) - 0.4).abs() < 1e-6,
            "magnitude {magnitude}: volume {}",
            mean(&density)
        );
        assert!(multiplier > 0.0 && multiplier.is_finite());
    }
}

#[test]
fn oc_update_without_any_demand_for_material_falls_to_the_lower_move_bound() {
    let grid = GridSize::new(2, 2);
    let simp = SiMP::new(grid, 0.5, 3.0, 1.0, MaterialParams::default());
    let solver = TopOptSolver::new(simp, FilteringMethods::new(1.0, grid)).with_move_limit(0.1);
    let (density, multiplier) = solver.oc_update(&[0.5; 4], &[0.0; 4], 0.5);
    assert_eq!(multiplier, 0.0);
    for rho in density {
        assert!((rho - 0.4).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// Casting filter
// ---------------------------------------------------------------------------

fn casting(nx: usize, ny: usize, direction: usize) -> ManufacturingFilters {
    ManufacturingFilters::new(GridSize::new(nx, ny), 2.0, PI / 4.0).with_build_direction(direction)
}

fn grid_of(nx: usize, ny: usize, solid: &[(usize, usize)]) -> Vec<f64> {
    let mut density = vec![0.0; nx * ny];
    for &(row, col) in solid {
        density[row * nx + col] = 1.0;
    }
    density
}

fn solid_cells(density: &[f64], nx: usize) -> Vec<(usize, usize)> {
    density
        .iter()
        .enumerate()
        .filter(|&(_, &rho)| rho > 0.5)
        .map(|(i, _)| (i / nx, i % nx))
        .collect()
}

#[test]
fn casting_fills_the_column_down_to_the_parting_plane_for_every_direction() {
    // Pull +y: parting plane at row 0, so a pixel in the top row drags its
    // whole column down.
    let out = casting(3, 3, 0).casting_filter(&grid_of(3, 3, &[(2, 1)]));
    assert_eq!(solid_cells(&out, 3), vec![(0, 1), (1, 1), (2, 1)]);

    // Pull -y: parting plane at the last row.
    let out = casting(3, 3, 1).casting_filter(&grid_of(3, 3, &[(0, 1)]));
    assert_eq!(solid_cells(&out, 3), vec![(0, 1), (1, 1), (2, 1)]);

    // Pull +x: parting plane at column 0.
    let out = casting(3, 3, 2).casting_filter(&grid_of(3, 3, &[(1, 2)]));
    assert_eq!(solid_cells(&out, 3), vec![(1, 0), (1, 1), (1, 2)]);

    // Pull -x: parting plane at the last column.
    let out = casting(3, 3, 3).casting_filter(&grid_of(3, 3, &[(1, 0)]));
    assert_eq!(solid_cells(&out, 3), vec![(1, 0), (1, 1), (1, 2)]);
}

#[test]
fn casting_leaves_material_next_to_the_parting_plane_alone() {
    // A pixel on the parting plane has nothing between itself and the plane.
    assert_eq!(
        solid_cells(
            &casting(3, 3, 0).casting_filter(&grid_of(3, 3, &[(0, 1)])),
            3
        ),
        vec![(0, 1)]
    );
    assert_eq!(
        solid_cells(
            &casting(3, 3, 1).casting_filter(&grid_of(3, 3, &[(2, 1)])),
            3
        ),
        vec![(2, 1)]
    );
    assert_eq!(
        solid_cells(
            &casting(3, 3, 2).casting_filter(&grid_of(3, 3, &[(1, 0)])),
            3
        ),
        vec![(1, 0)]
    );
    assert_eq!(
        solid_cells(
            &casting(3, 3, 3).casting_filter(&grid_of(3, 3, &[(1, 2)])),
            3
        ),
        vec![(1, 2)]
    );
}

#[test]
fn casting_leaves_an_undercut_free_design_unchanged() {
    // Pull +y: every column is solid from row 0 up to some height.
    let heights = [3_usize, 0, 2, 1, 3];
    let mut density = vec![0.0; 5 * 4];
    for (col, &height) in heights.iter().enumerate() {
        for row in 0..height {
            density[row * 5 + col] = 1.0;
        }
    }
    let out = casting(5, 4, 0).casting_filter(&density);
    assert_eq!(out, density);
}

#[test]
fn casting_removes_undercuts_so_the_result_is_monotone_idempotent_and_additive() {
    let (nx, ny) = (5, 3);
    // Deterministic pseudo-random densities in [0, 1].
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    for direction in 0..4 {
        let filter = casting(nx, ny, direction);
        let density: Vec<f64> = (0..nx * ny).map(|_| next()).collect();
        let out = filter.casting_filter(&density);
        assert_eq!(out.len(), density.len());

        for (before, after) in density.iter().zip(out.iter()) {
            assert!(
                after >= before,
                "direction {direction}: material was removed"
            );
            assert!((0.0..=1.0).contains(after));
        }
        // Along the pull axis the density never increases away from the parting plane.
        for row in 0..ny {
            for col in 0..nx {
                let here = out[row * nx + col];
                match direction {
                    0 if row + 1 < ny => assert!(here >= out[(row + 1) * nx + col]),
                    1 if row > 0 => assert!(here >= out[(row - 1) * nx + col]),
                    2 if col + 1 < nx => assert!(here >= out[row * nx + col + 1]),
                    3 if col > 0 => assert!(here >= out[row * nx + col - 1]),
                    _ => {}
                }
            }
        }
        assert_eq!(
            filter.casting_filter(&out),
            out,
            "direction {direction}: idempotence"
        );
    }
}

#[test]
fn casting_never_panics_on_degenerate_input() {
    let filter = casting(4, 4, 0);
    assert_eq!(filter.casting_filter(&[0.3; 5]), vec![0.3; 5]);
    assert!(filter.casting_filter(&[]).is_empty());
    let empty_grid = casting(0, 0, 2);
    assert!(empty_grid.casting_filter(&[]).is_empty());
    let single_row = casting(4, 1, 0);
    assert_eq!(
        single_row.casting_filter(&[0.0, 1.0, 0.5, 0.0]),
        vec![0.0, 1.0, 0.5, 0.0]
    );
}
