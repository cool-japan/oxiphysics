// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the interior-penalty diffusion terms of the DG solver.

use super::*;

const METHODS: [PenaltyMethod; 3] = [
    PenaltyMethod::Sipg,
    PenaltyMethod::Nipg,
    PenaltyMethod::Iipg,
];

/// Deterministic pseudo-random numbers in [0, 1).
fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / ((1u64 << 53) as f64)
}

/// Elements covering consecutive intervals of the given widths, starting at 0.
fn mesh_from_widths(widths: &[f64], order: usize) -> Vec<DgElement> {
    let mut x = 0.0;
    widths
        .iter()
        .map(|&w| {
            let e = DgElement::new(x, x + w, order);
            x += w;
            e
        })
        .collect()
}

fn randomize(elements: &mut [DgElement], seed: u64) {
    let mut state = seed;
    for e in elements.iter_mut() {
        for c in e.coeffs.iter_mut() {
            *c = 2.0 * lcg(&mut state) - 1.0;
        }
    }
}

const NONUNIFORM: [f64; 6] = [0.2, 0.35, 0.15, 0.3, 0.25, 0.2];

/// Pure-diffusion residual (no advection).
fn diffusion_residual(
    elements: &[DgElement],
    kappa: f64,
    method: PenaltyMethod,
    c_pen: f64,
) -> Vec<Vec<f64>> {
    dg_advection_diffusion_residual(
        elements,
        0.0,
        kappa,
        method,
        c_pen,
        NumericalFlux::LaxFriedrichs,
    )
}

#[test]
fn test_legendre_deriv_boundary_values() {
    // Regression: P_n'(-1) had the wrong sign for n >= 2.
    assert_eq!(legendre_poly_deriv(0, 1.0), 0.0);
    for n in 1..=8_usize {
        let magnitude = (n * (n + 1)) as f64 / 2.0;
        let sign_at_minus_one = if n % 2 == 1 { 1.0 } else { -1.0 };
        assert!(
            (legendre_poly_deriv(n, 1.0) - magnitude).abs() < 1e-12,
            "P_{n}'(1)"
        );
        assert!(
            (legendre_poly_deriv(n, -1.0) - sign_at_minus_one * magnitude).abs() < 1e-12,
            "P_{n}'(-1) = {}",
            legendre_poly_deriv(n, -1.0)
        );
    }
    // Explicit cases: P2' = 3x, P3' = (15x^2 - 3)/2.
    assert!((legendre_poly_deriv(2, -1.0) + 3.0).abs() < 1e-12);
    assert!((legendre_poly_deriv(3, -1.0) - 6.0).abs() < 1e-12);
    // The boundary value is the limit of the interior formula.
    for n in 2..=6_usize {
        let inside = legendre_poly_deriv(n, -1.0 + 1e-7);
        assert!((inside - legendre_poly_deriv(n, -1.0)).abs() < 1e-4 * inside.abs().max(1.0));
    }
}

#[test]
fn test_diffusion_residual_is_exact_for_smooth_continuous_data() {
    // u = x^3 is C^1 and lies in the P3 space, so on every element that has a neighbour on both
    // sides the residual equals the L2 projection of kappa * u'' = 6 kappa x for every method.
    let kappa = 0.7;
    for method in METHODS {
        let mut mesh = DgMesh1D {
            elements: mesh_from_widths(&NONUNIFORM, 3),
        };
        mesh.project_initial_condition(|x| x * x * x);
        let residual = diffusion_residual(&mesh.elements, kappa, method, 10.0);
        let mut target = DgMesh1D {
            elements: mesh_from_widths(&NONUNIFORM, 3),
        };
        target.project_initial_condition(|x| 6.0 * kappa * x);
        for (i, (res, tgt)) in residual
            .iter()
            .zip(target.elements.iter())
            .enumerate()
            .skip(1)
            .take(NONUNIFORM.len() - 2)
        {
            for (m, (&got, &want)) in res.iter().zip(tgt.coeffs.iter()).enumerate() {
                assert!(
                    (got - want).abs() < 1e-9 * (1.0 + want.abs()),
                    "{method:?}: element {i} mode {m}: {got} vs projection {want}"
                );
            }
        }
    }
}

#[test]
fn test_diffusion_decays_a_smooth_bump_instead_of_amplifying_it() {
    // Regression: the volume term had the anti-diffusive sign. At the crest of a Gaussian (u'' < 0)
    // the cell mean must decrease, in the convex flank (u'' > 0) it must increase.  The L2 projection
    // of the bump jumps across element faces, so the numerical flux differs noticeably from the exact
    // one; only the sign and the order of magnitude are asserted here (exactness is covered by the
    // polynomial test above).
    let kappa = 0.3;
    let mut mesh = DgMesh1D::uniform(0.0, 4.0, 16, 2);
    mesh.project_initial_condition(|x| (-(x - 2.0) * (x - 2.0)).exp());
    let residual = diffusion_residual(&mesh.elements, kappa, PenaltyMethod::Sipg, 10.0);
    let mean_rate = |i: usize| residual[i][0] / 2.0_f64.sqrt();
    // Elements 7 and 8 meet at the crest x = 2 (kappa u'' = -0.6 there).
    for crest in [7, 8] {
        assert!(
            mean_rate(crest) < -0.2 && mean_rate(crest) > -0.7,
            "crest element {crest}: mean rate {}",
            mean_rate(crest)
        );
    }
    // The two crest elements are mirror images of each other.
    assert!((mean_rate(7) - mean_rate(8)).abs() < 1e-9);
    assert!(mean_rate(3) > 0.05, "flank mean rate {}", mean_rate(3));
    assert!(mean_rate(12) > 0.05, "flank mean rate {}", mean_rate(12));
}

#[test]
fn test_penalty_term_acts_on_both_neighbours() {
    // p = 0: u_L = c_L/sqrt2, u_R = 0.  Regression: only the left element used to be updated.
    let kappa = 0.4;
    let c_pen = 5.0;
    for method in METHODS {
        let mut elements = vec![DgElement::new(0.0, 1.0, 0), DgElement::new(1.0, 2.0, 0)];
        elements[0].coeffs[0] = 2.0;
        let res = diffusion_residual(&elements, kappa, method, c_pen);
        // eta = c_pen / h = 5, [u] = sqrt2, flux = -kappa eta [u]; res = +-flux phi0 / J.
        let expected = 2.0 * kappa * c_pen;
        assert!(
            (res[0][0] + expected).abs() < 1e-12,
            "{method:?}: left {}",
            res[0][0]
        );
        assert!(
            (res[1][0] - expected).abs() < 1e-12,
            "{method:?}: right {}",
            res[1][0]
        );
    }
}

#[test]
fn test_diffusion_conserves_mass_for_every_method() {
    for method in METHODS {
        for order in 0..=3_usize {
            let mut elements = mesh_from_widths(&NONUNIFORM, order);
            randomize(&mut elements, 11 + order as u64);
            let residual = diffusion_residual(&elements, 0.9, method, 10.0);
            // d/dt of the integral: sum_i J_i * sqrt2 * res_i[0].
            let rate: f64 = elements
                .iter()
                .zip(residual.iter())
                .map(|(e, r)| e.jacobian() * 2.0_f64.sqrt() * r[0])
                .sum();
            let scale: f64 = residual.iter().map(|r| r[0].abs()).sum::<f64>();
            assert!(
                rate.abs() < 1e-12 * (1.0 + scale),
                "{method:?} p = {order}: d(mass)/dt = {rate}"
            );
        }
    }
}

/// The dense operator A with d(coeffs)/dt = A coeffs.
fn operator_matrix(
    widths: &[f64],
    order: usize,
    kappa: f64,
    method: PenaltyMethod,
    c_pen: f64,
) -> Vec<Vec<f64>> {
    let n_modes = order + 1;
    let n_dof = widths.len() * n_modes;
    let mut matrix = vec![vec![0.0; n_dof]; n_dof];
    for col in 0..n_dof {
        let mut elements = mesh_from_widths(widths, order);
        elements[col / n_modes].coeffs[col % n_modes] = 1.0;
        let residual = diffusion_residual(&elements, kappa, method, c_pen);
        for (i, r) in residual.iter().enumerate() {
            for (m, &v) in r.iter().enumerate() {
                matrix[i * n_modes + m][col] = v;
            }
        }
    }
    matrix
}

#[test]
fn test_sipg_operator_is_symmetric_and_the_others_are_not() {
    let order = 2;
    let n_modes = order + 1;
    let jac = |dof: usize| NONUNIFORM[dof / n_modes] / 2.0;
    for method in METHODS {
        let a = operator_matrix(&NONUNIFORM, order, 0.8, method, 10.0);
        let mut asymmetry = 0.0_f64;
        let mut scale = 0.0_f64;
        for (i, row) in a.iter().enumerate() {
            for (j, &a_ij) in row.iter().enumerate() {
                // The operator is symmetric in the J-weighted inner product (mass matrix J I).
                asymmetry = asymmetry.max((jac(i) * a_ij - jac(j) * a[j][i]).abs());
                scale = scale.max((jac(i) * a_ij).abs());
            }
        }
        match method {
            PenaltyMethod::Sipg => assert!(
                asymmetry < 1e-10 * scale,
                "SIPG must be symmetric, asymmetry {asymmetry:e} (scale {scale:e})"
            ),
            _ => assert!(
                asymmetry > 1e-3 * scale,
                "{method:?} must be non-symmetric, asymmetry {asymmetry:e} (scale {scale:e})"
            ),
        }
    }
}

#[test]
fn test_energy_identity_matches_interior_penalty_form() {
    // d/dt (1/2 ||u||^2) = sum_i J_i c_i . res_i = -(kappa sum_K int u_x^2 + sum_faces a_e(u,u)).
    let kappa = 0.6;
    let c_pen = 8.0;
    for method in METHODS {
        for order in 1..=3_usize {
            let mut elements = mesh_from_widths(&NONUNIFORM, order);
            randomize(&mut elements, 100 + order as u64);
            let residual = diffusion_residual(&elements, kappa, method, c_pen);
            let rate: f64 = elements
                .iter()
                .zip(residual.iter())
                .map(|(e, r)| {
                    e.jacobian()
                        * e.coeffs
                            .iter()
                            .zip(r.iter())
                            .map(|(c, x)| c * x)
                            .sum::<f64>()
                })
                .sum();

            let mut dissipation = 0.0;
            for e in &elements {
                let basis = DgBasis::new(e.order);
                let (qpts, qwts) = gauss_legendre_quadrature(e.order + 2);
                for (k, &xi) in qpts.iter().enumerate() {
                    let du: f64 = e
                        .coeffs
                        .iter()
                        .enumerate()
                        .map(|(m, &c)| c * basis.evaluate_deriv(m, xi))
                        .sum();
                    dissipation += kappa * qwts[k] * du * du / e.jacobian();
                }
            }
            for pair in elements.windows(2) {
                let (l, r) = (&pair[0], &pair[1]);
                let (bl, br) = (DgBasis::new(l.order), DgBasis::new(r.order));
                let grad_l = l
                    .coeffs
                    .iter()
                    .enumerate()
                    .map(|(m, &c)| c * bl.evaluate_deriv(m, 1.0))
                    .sum::<f64>()
                    / l.jacobian();
                let grad_r = r
                    .coeffs
                    .iter()
                    .enumerate()
                    .map(|(m, &c)| c * br.evaluate_deriv(m, -1.0))
                    .sum::<f64>()
                    / r.jacobian();
                let eta = penalty_parameter(l.order.max(r.order), l.width(), r.width(), c_pen);
                let (half_l, half_r) = interior_penalty_face(
                    l.eval_right(),
                    r.eval_left(),
                    grad_l,
                    grad_r,
                    kappa,
                    eta,
                    method,
                );
                dissipation += half_l + half_r;
            }
            assert!(
                (rate + dissipation).abs() < 1e-10 * (1.0 + dissipation.abs()),
                "{method:?} p = {order}: d(E)/dt = {rate}, -dissipation = {}",
                -dissipation
            );
        }
    }
}

#[test]
fn test_sipg_diffusion_never_increases_the_l2_norm() {
    for order in 1..=3_usize {
        let mut elements = mesh_from_widths(&NONUNIFORM, order);
        randomize(&mut elements, 7 * order as u64);
        let residual = diffusion_residual(&elements, 1.0, PenaltyMethod::Sipg, 10.0);
        let rate: f64 = elements
            .iter()
            .zip(residual.iter())
            .map(|(e, r)| {
                e.jacobian()
                    * e.coeffs
                        .iter()
                        .zip(r.iter())
                        .map(|(c, x)| c * x)
                        .sum::<f64>()
            })
            .sum();
        assert!(
            rate < 0.0,
            "p = {order}: d(||u||^2/2)/dt = {rate} must be negative"
        );
    }
}

#[test]
fn test_interior_penalty_face_energy_form() {
    // SIPG with continuous u: no contribution.  Pure jump: kappa * eta * jump^2 for NIPG.
    let (l, r) = interior_penalty_face(1.0, 1.0, 3.0, 3.0, 0.5, 10.0, PenaltyMethod::Sipg);
    assert_eq!((l, r), (0.0, 0.0));
    let (l, r) = interior_penalty_face(2.0, 0.0, 0.0, 0.0, 0.5, 10.0, PenaltyMethod::Nipg);
    assert!((l - 0.5 * 0.5 * 10.0 * 4.0).abs() < 1e-14 && l == r);
    // SIPG: kappa eta [u]^2 - 2 {kappa u_x} [u]; IIPG: kappa eta [u]^2 - {kappa u_x} [u].
    let (l, r) = interior_penalty_face(1.0, 0.0, 2.0, 4.0, 0.5, 10.0, PenaltyMethod::Sipg);
    let expected = 0.5 * 10.0 - 2.0 * (0.5 * 3.0);
    assert!((l + r - expected).abs() < 1e-14, "{}", l + r);
    let (l, r) = interior_penalty_face(1.0, 0.0, 2.0, 4.0, 0.5, 10.0, PenaltyMethod::Iipg);
    assert!((l + r - (0.5 * 10.0 - 1.5)).abs() < 1e-14);
}

#[test]
fn test_gaussian_diffusion_conserves_mass_and_decays_monotonically() {
    // a = 0, kappa = 0.01: Gaussian of width 0.5 spreads; mass is conserved, ||u|| and max(u) decay.
    let kappa = 0.01;
    let sigma0 = 0.5;
    let gaussian = |x: f64, t: f64| {
        let s2 = sigma0 * sigma0 + 2.0 * kappa * t;
        sigma0 / s2.sqrt() * (-x * x / (2.0 * s2)).exp()
    };
    let mut mesh = DgMesh1D::uniform(-4.0, 4.0, 40, 2);
    mesh.project_initial_condition(|x| gaussian(x, 0.0));
    let config = DgAdvDiffConfig {
        wave_speed: 0.0,
        kappa,
        flux_type: NumericalFlux::LaxFriedrichs,
        penalty_method: PenaltyMethod::Sipg,
        c_pen: 10.0,
        cfl: 1.0,
    };
    let mass0 = total_mass(&mesh);
    let mut norm_prev = mesh.l2_norm();
    let mut max_prev = solution_max(&mesh);
    let chunk = 0.25;
    for step in 1..=8 {
        let reached = solve_advection_diffusion(&mut mesh, &config, chunk);
        assert!((reached - chunk).abs() < 1e-12);
        let mass = total_mass(&mesh);
        assert!(
            (mass - mass0).abs() < 1e-11 * mass0.abs(),
            "step {step}: mass {mass} vs {mass0}"
        );
        let norm = mesh.l2_norm();
        let max = solution_max(&mesh);
        assert!(
            norm < norm_prev,
            "step {step}: ||u|| {norm} did not decay from {norm_prev}"
        );
        assert!(
            max < max_prev,
            "step {step}: max {max} did not decay from {max_prev}"
        );
        norm_prev = norm;
        max_prev = max;
    }
    let t_final = 8.0 * chunk;
    let error = mesh.l2_error(|x| gaussian(x, t_final));
    assert!(
        error < 1e-3,
        "L2 error against the exact Gaussian: {error:e}"
    );
    // The peak really decreased towards the exact value sigma0 / sqrt(sigma0^2 + 2 kappa t).
    let exact_peak = gaussian(0.0, t_final);
    assert!((solution_max(&mesh) - exact_peak).abs() < 2e-3);
}

/// L2 errors of u_t = kappa u_xx, u0 = cos(x) on [0, pi] (homogeneous Neumann) at T = 0.05.
fn cosine_decay_errors(order: usize, method: PenaltyMethod) -> Vec<f64> {
    let kappa = 1.0_f64;
    let t_final = 0.05_f64;
    let exact = |x: f64| (-kappa * t_final).exp() * x.cos();
    [4_usize, 8, 16]
        .iter()
        .map(|&n| {
            let mut mesh = DgMesh1D::uniform(0.0, std::f64::consts::PI, n, order);
            mesh.project_initial_condition(f64::cos);
            let config = DgAdvDiffConfig {
                wave_speed: 0.0,
                kappa,
                flux_type: NumericalFlux::Upwind,
                penalty_method: method,
                c_pen: 10.0,
                cfl: 1.5,
            };
            let reached = solve_advection_diffusion(&mut mesh, &config, t_final);
            assert!((reached - t_final).abs() < 1e-12);
            mesh.l2_error(exact)
        })
        .collect()
}

#[test]
fn test_diffusion_converges_at_order_p_plus_one() {
    for method in [PenaltyMethod::Sipg, PenaltyMethod::Nipg] {
        for order in 1..=2_usize {
            let errors = cosine_decay_errors(order, method);
            for (k, pair) in errors.windows(2).enumerate() {
                let rate = (pair[0] / pair[1]).log2();
                assert!(
                    rate > order as f64 + 0.8,
                    "{method:?} p = {order}: rate {rate} between levels {k} and {} (errors {errors:?})",
                    k + 1
                );
            }
        }
    }
    // Absolute accuracy of the finest level (p = 1: ~1.7e-3, p = 2: ~5e-5).
    let p1 = cosine_decay_errors(1, PenaltyMethod::Sipg);
    assert!(p1[2] < 2.5e-3, "{p1:?}");
    let p2 = cosine_decay_errors(2, PenaltyMethod::Sipg);
    assert!(p2[2] < 8e-5, "{p2:?}");
}

#[test]
fn test_advection_diffusion_travelling_gaussian() {
    // u_t + u_x = kappa u_xx: the Gaussian drifts with speed 1 and spreads diffusively.
    let (a, kappa, sigma0) = (1.0, 0.01, 0.5);
    let exact = |x: f64, t: f64| {
        let s2 = sigma0 * sigma0 + 2.0 * kappa * t;
        sigma0 / s2.sqrt() * (-(x - 3.0 - a * t).powi(2) / (2.0 * s2)).exp()
    };
    let mut mesh = DgMesh1D::uniform(0.0, 12.0, 60, 2);
    mesh.project_initial_condition(|x| exact(x, 0.0));
    let config = DgAdvDiffConfig {
        wave_speed: a,
        kappa,
        flux_type: NumericalFlux::Upwind,
        penalty_method: PenaltyMethod::Sipg,
        c_pen: 10.0,
        cfl: 1.0,
    };
    let t_final = 2.0;
    let reached = solve_advection_diffusion(&mut mesh, &config, t_final);
    assert!((reached - t_final).abs() < 1e-12);
    let error = mesh.l2_error(|x| exact(x, t_final));
    assert!(error < 2e-3, "L2 error {error:e}");
    // The pulse has been advected by a * t = 2, from x = 3 to x = 5.
    let peak_x = (0..=1200)
        .map(|k| 0.01 * k as f64)
        .max_by(|&x, &y| {
            mesh.evaluate_at(x)
                .partial_cmp(&mesh.evaluate_at(y))
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(0.0);
    assert!((peak_x - 5.0).abs() < 0.1, "peak at x = {peak_x}");
}

#[test]
fn test_advection_diffusion_mass_balance_is_exact_up_to_boundary_flux() {
    // With the extrapolation boundary condition of the advection part, d(mass)/dt equals the
    // advective boundary flux a * (u(x_left^+) - u(x_right^-)); diffusion adds nothing.
    let a = 1.3;
    for method in METHODS {
        for order in 1..=3_usize {
            let mut elements = mesh_from_widths(&NONUNIFORM, order);
            randomize(&mut elements, 500 + order as u64);
            let residual = dg_advection_diffusion_residual(
                &elements,
                a,
                0.4,
                method,
                10.0,
                NumericalFlux::Upwind,
            );
            let rate: f64 = elements
                .iter()
                .zip(residual.iter())
                .map(|(e, r)| e.jacobian() * 2.0_f64.sqrt() * r[0])
                .sum();
            let first = elements[0].eval_left();
            let last = elements[elements.len() - 1].eval_right();
            let expected = a * (first - last);
            assert!(
                (rate - expected).abs() < 1e-11 * (1.0 + expected.abs()),
                "{method:?} p = {order}: d(mass)/dt = {rate}, boundary flux balance = {expected}"
            );
        }
    }
}

#[test]
fn test_time_step_estimate_is_stable_for_rough_data() {
    // Random (non-smooth) data on a non-uniform mesh must not blow up for any method / order at the
    // advertised CFL numbers; the energy of the pure diffusion problem decays.
    for method in METHODS {
        for order in 1..=3_usize {
            let elements = mesh_from_widths(&NONUNIFORM, order);
            let mut mesh = DgMesh1D { elements };
            randomize(&mut mesh.elements, 31 + order as u64);
            let config = DgAdvDiffConfig {
                wave_speed: 0.0,
                kappa: 0.5,
                flux_type: NumericalFlux::LaxFriedrichs,
                penalty_method: method,
                c_pen: 10.0,
                cfl: 2.0,
            };
            let initial = mesh.l2_norm();
            let reached = solve_advection_diffusion(&mut mesh, &config, 0.05);
            assert!((reached - 0.05).abs() < 1e-12);
            let after = mesh.l2_norm();
            assert!(
                after.is_finite() && after < initial,
                "{method:?} p = {order}: {initial} -> {after}"
            );
        }
    }
}

#[test]
fn test_unusable_time_step_advances_nothing() {
    // cfl <= 0 or NaN used to spin forever (zero step) or take one giant step (NaN step).
    for cfl in [0.0, -1.0, f64::NAN] {
        let mut mesh = DgMesh1D::uniform(0.0, 1.0, 4, 1);
        mesh.project_initial_condition(|x| 1.0 + x);
        let before: Vec<Vec<f64>> = mesh.elements.iter().map(|e| e.coeffs.clone()).collect();
        let config = DgAdvDiffConfig {
            wave_speed: 0.5,
            kappa: 0.1,
            flux_type: NumericalFlux::Upwind,
            penalty_method: PenaltyMethod::Sipg,
            c_pen: 10.0,
            cfl,
        };
        let reached = solve_advection_diffusion(&mut mesh, &config, 1.0);
        assert_eq!(reached, 0.0, "cfl = {cfl}");
        for (element, old) in mesh.elements.iter().zip(before.iter()) {
            assert_eq!(&element.coeffs, old, "cfl = {cfl}");
        }
    }
    let mut empty = DgMesh1D { elements: vec![] };
    assert_eq!(
        solve_advection_diffusion(&mut empty, &DgAdvDiffConfig::default(), 1.0),
        0.0
    );
}
