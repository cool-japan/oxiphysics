//! Tests of the condensed-QP model predictive controller
//! ([`ModelPredictiveControl`]).
//!
//! None of the oracles uses the machinery of the solver (adjoint gradient,
//! accelerated projected gradient, augmented Lagrangian):
//!
//! * the finite-horizon LQR of this crate for the unconstrained problem;
//! * a dense condensed model `J(u) = ½ uᵀHu + fᵀu + c₀`, `x = s₀ + G u`,
//!   assembled by simulating unit inputs (exact for a linear system), solved by
//!   brute-force enumeration of the active sets: every combination of input and
//!   state bounds is tried, each equality-constrained KKT system is solved by LU
//!   and the one with non-negative multipliers and feasible inputs is the
//!   optimum;
//! * the Karush–Kuhn–Tucker conditions of larger problems, checked with the
//!   same dense model;
//! * hand-computed and numpy reference values.

use super::mpc_qp::{cost_and_gradient, curvature_estimate};
use super::types::{FiniteHorizonLqr, ModelPredictiveControl, MpcProblem, MpcSolution};
use crate::dense_linalg::{lu_solve, real_eigenvalues};

const DT: f64 = 0.1;
const INF: f64 = f64::INFINITY;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}

fn norm_inf(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |m, v| m.max(v.abs()))
}

fn diagonal(entries: &[f64]) -> Vec<f64> {
    let n = entries.len();
    let mut d = vec![0.0; n * n];
    for (i, v) in entries.iter().enumerate() {
        d[i * n + i] = *v;
    }
    d
}

/// Double integrator `x = (position, velocity)` sampled with `DT`, weights
/// `Q = diag(q)`, `R = r`, `Qf = diag(qf)` and no bounds at all.
fn double_integrator(horizon: usize, q: [f64; 2], r: f64, qf: [f64; 2]) -> MpcProblem {
    MpcProblem {
        horizon,
        n: 2,
        m: 1,
        a: vec![1.0, DT, 0.0, 1.0],
        b: vec![0.5 * DT * DT, DT],
        q: diagonal(&q),
        qf: diagonal(&qf),
        r: vec![r],
        x_lb: vec![-INF; 2],
        x_ub: vec![INF; 2],
        u_lb: vec![-INF],
        u_ub: vec![INF],
    }
}

/// The reference problem: `Q = diag(10, 1)`, `R = 0.1`, `Qf = 2 Q`.
fn reference_problem(horizon: usize) -> MpcProblem {
    double_integrator(horizon, [10.0, 1.0], 0.1, [20.0, 2.0])
}

fn solve(problem: &MpcProblem, x0: &[f64]) -> MpcSolution {
    ModelPredictiveControl::new(problem.clone()).solve(x0)
}

// ── Independent oracles ─────────────────────────────────────────────────

/// State trajectory `(N + 1) x n` of `x_{k+1} = A x_k + B u_k`.
fn simulate(p: &MpcProblem, x0: &[f64], u: &[f64]) -> Vec<f64> {
    let (n, m) = (p.n, p.m);
    let mut traj = x0.to_vec();
    for k in 0..p.horizon {
        let x = traj[k * n..(k + 1) * n].to_vec();
        for i in 0..n {
            let ax: f64 = (0..n).map(|j| p.a[i * n + j] * x[j]).sum();
            let bu: f64 = (0..m).map(|j| p.b[i * m + j] * u[k * m + j]).sum();
            traj.push(ax + bu);
        }
    }
    traj
}

/// `vᵀ M v` of a full (not necessarily symmetric) square matrix.
fn quadratic_form(matrix: &[f64], v: &[f64]) -> f64 {
    let n = v.len();
    (0..n)
        .map(|i| v[i] * (0..n).map(|j| matrix[i * n + j] * v[j]).sum::<f64>())
        .sum()
}

/// `J(u) = Σ (xₖᵀ Q xₖ + uₖᵀ R uₖ) + x_Nᵀ Qf x_N`.
fn total_cost(p: &MpcProblem, x0: &[f64], u: &[f64]) -> f64 {
    let (n, m) = (p.n, p.m);
    let traj = simulate(p, x0, u);
    let mut cost = quadratic_form(&p.qf, &traj[p.horizon * n..]);
    for k in 0..p.horizon {
        cost += quadratic_form(&p.q, &traj[k * n..(k + 1) * n]);
        cost += quadratic_form(&p.r, &u[k * m..(k + 1) * m]);
    }
    cost
}

/// Dense condensed model of the problem.
struct Condensed {
    dim: usize,
    /// Hessian of `J` (`dim x dim`).
    h: Vec<f64>,
    /// Linear term of `J`.
    f: Vec<f64>,
    /// Predicted states for zero inputs, `(N + 1) n` entries.
    s0: Vec<f64>,
    /// Response of the predicted states to the inputs, `(N + 1) n x dim`.
    g: Vec<f64>,
}

fn condense(p: &MpcProblem, x0: &[f64]) -> Condensed {
    let dim = p.horizon * p.m;
    let unit = |i: usize, value: f64| {
        let mut u = vec![0.0; dim];
        u[i] = value;
        u
    };
    let cost = |u: &[f64]| total_cost(p, x0, u);
    let c0 = cost(&vec![0.0; dim]);
    let mut h = vec![0.0; dim * dim];
    let mut f = vec![0.0; dim];
    for i in 0..dim {
        let (plus, minus) = (cost(&unit(i, 1.0)), cost(&unit(i, -1.0)));
        h[i * dim + i] = plus + minus - 2.0 * c0;
        f[i] = 0.5 * (plus - minus);
    }
    for i in 0..dim {
        for j in (i + 1)..dim {
            let mut both = unit(i, 1.0);
            both[j] = 1.0;
            let value = cost(&both) - cost(&unit(i, 1.0)) - cost(&unit(j, 1.0)) + c0;
            h[i * dim + j] = value;
            h[j * dim + i] = value;
        }
    }
    let zero_state = vec![0.0; p.n];
    let s0 = simulate(p, x0, &vec![0.0; dim]);
    let mut g = vec![0.0; s0.len() * dim];
    for j in 0..dim {
        let response = simulate(p, &zero_state, &unit(j, 1.0));
        for (r, value) in response.iter().enumerate() {
            g[r * dim + j] = *value;
        }
    }
    Condensed { dim, h, f, s0, g }
}

/// Linear inequality `aᵀ u <= b` on the input sequence.
struct Constraint {
    a: Vec<f64>,
    b: f64,
}

/// All finite input and state bounds (the states of the stages `1..=N`).
fn constraints(p: &MpcProblem, c: &Condensed) -> Vec<Constraint> {
    let dim = c.dim;
    let mut out = Vec::new();
    for k in 0..p.horizon {
        for j in 0..p.m {
            let mut e = vec![0.0; dim];
            e[k * p.m + j] = 1.0;
            if p.u_ub[j].is_finite() {
                out.push(Constraint {
                    a: e.clone(),
                    b: p.u_ub[j],
                });
            }
            if p.u_lb[j].is_finite() {
                out.push(Constraint {
                    a: e.iter().map(|v| -v).collect(),
                    b: -p.u_lb[j],
                });
            }
        }
    }
    for k in 1..=p.horizon {
        for i in 0..p.n {
            let row = k * p.n + i;
            let gradient = &c.g[row * dim..(row + 1) * dim];
            if p.x_ub[i].is_finite() {
                out.push(Constraint {
                    a: gradient.to_vec(),
                    b: p.x_ub[i] - c.s0[row],
                });
            }
            if p.x_lb[i].is_finite() {
                out.push(Constraint {
                    a: gradient.iter().map(|v| -v).collect(),
                    b: c.s0[row] - p.x_lb[i],
                });
            }
        }
    }
    out
}

/// Exact minimiser by enumeration of the active sets: the KKT system of every
/// subset of the constraints is solved and the first candidate with feasible
/// inputs and non-negative multipliers is returned.
fn enumerated_optimum(p: &MpcProblem, x0: &[f64]) -> Vec<f64> {
    let c = condense(p, x0);
    let cons = constraints(p, &c);
    assert!(
        cons.len() <= 18,
        "{} constraints are too many to enumerate",
        cons.len()
    );
    let dim = c.dim;
    for mask in 0_u32..(1_u32 << cons.len()) {
        let active: Vec<&Constraint> = (0..cons.len())
            .filter(|k| (mask >> k) & 1 == 1)
            .map(|k| &cons[k])
            .collect();
        let size = dim + active.len();
        let mut kkt = vec![0.0; size * size];
        let mut rhs = vec![0.0; size];
        for i in 0..dim {
            kkt[i * size..i * size + dim].copy_from_slice(&c.h[i * dim..(i + 1) * dim]);
            rhs[i] = -c.f[i];
        }
        for (r, con) in active.iter().enumerate() {
            for j in 0..dim {
                kkt[(dim + r) * size + j] = con.a[j];
                kkt[j * size + dim + r] = con.a[j];
            }
            rhs[dim + r] = con.b;
        }
        let Some(solution) = lu_solve(&kkt, size, &rhs, 1) else {
            continue;
        };
        let (u, lambda) = solution.split_at(dim);
        let feasible = cons
            .iter()
            .all(|con| dot(&con.a, u) <= con.b + 1e-9 * (1.0 + con.b.abs()));
        if feasible && lambda.iter().all(|&l| l >= -1e-9) {
            return u.to_vec();
        }
    }
    panic!("the enumeration found no KKT point");
}

/// Assert the Karush–Kuhn–Tucker conditions of the condensed problem at the
/// solver's inputs: feasibility, stationarity of the Lagrangian with
/// least-squares multipliers on the active constraints, and `λ >= 0`.
fn assert_kkt_point(p: &MpcProblem, x0: &[f64], sol: &MpcSolution, tol: f64) {
    let c = condense(p, x0);
    let cons = constraints(p, &c);
    let (dim, u) = (c.dim, &sol.u_seq);
    let slack: Vec<f64> = cons.iter().map(|con| con.b - dot(&con.a, u)).collect();
    for (k, s) in slack.iter().enumerate() {
        assert!(*s >= -tol, "constraint {k} is violated by {}", -s);
    }
    let gradient: Vec<f64> = (0..dim)
        .map(|i| dot(&c.h[i * dim..(i + 1) * dim], u) + c.f[i])
        .collect();
    let scale = 1.0 + norm_inf(&gradient).max(norm_inf(&c.f));
    let active: Vec<usize> = (0..cons.len()).filter(|&k| slack[k] <= 1e-7).collect();
    if active.is_empty() {
        assert!(
            norm_inf(&gradient) <= tol * scale,
            "gradient {} at an interior point",
            norm_inf(&gradient)
        );
        return;
    }
    // Least-squares multipliers: (A Aᵀ) λ = -A ∇J.
    let na = active.len();
    let mut aat = vec![0.0; na * na];
    let mut rhs = vec![0.0; na];
    for (r, &kr) in active.iter().enumerate() {
        for (s, &ks) in active.iter().enumerate() {
            aat[r * na + s] = dot(&cons[kr].a, &cons[ks].a);
        }
        rhs[r] = -dot(&cons[kr].a, &gradient);
    }
    let lambda = lu_solve(&aat, na, &rhs, 1).expect("active constraints are independent");
    let mut residual = gradient.clone();
    for (r, &kr) in active.iter().enumerate() {
        for (res, a) in residual.iter_mut().zip(cons[kr].a.iter()) {
            *res += lambda[r] * a;
        }
    }
    assert!(
        norm_inf(&residual) <= tol * scale,
        "stationarity residual {} (scale {scale})",
        norm_inf(&residual)
    );
    let worst = lambda.iter().fold(0.0_f64, |m, l| m.max(-l));
    assert!(worst <= tol * scale, "negative multiplier {}", -worst);
}

fn assert_close(actual: &[f64], expected: &[f64], tol: f64, what: &str) {
    assert_eq!(actual.len(), expected.len(), "{what}: length");
    for (k, (a, e)) in actual.iter().zip(expected.iter()).enumerate() {
        assert!(
            (a - e).abs() <= tol * e.abs().max(1.0),
            "{what}[{k}]: {a} vs {e}"
        );
    }
}

/// Finite-horizon LQR open-loop inputs and states from `x0`.
fn lqr_reference(p: &MpcProblem, x0: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let lqr = FiniteHorizonLqr::new(p.n, p.m, &p.a, &p.b, &p.q, &p.r, &p.qf, p.horizon);
    let mut x = x0.to_vec();
    let (mut inputs, mut states) = (Vec::new(), x0.to_vec());
    for k in 0..p.horizon {
        let u = lqr.control(k, &x);
        x = simulate(
            &MpcProblem {
                horizon: 1,
                ..p.clone()
            },
            &x,
            &u,
        )[p.n..]
            .to_vec();
        inputs.extend_from_slice(&u);
        states.extend_from_slice(&x);
    }
    (inputs, states)
}

// ── Unconstrained problems ──────────────────────────────────────────────

#[test]
fn unconstrained_mpc_reproduces_the_finite_horizon_lqr_solution() {
    let p = reference_problem(12);
    let x0 = [1.0, -0.5];
    let sol = solve(&p, &x0);
    assert!(sol.converged && sol.feasible);
    assert!(
        sol.iterations > 0 && sol.iterations < 500,
        "{} iterations",
        sol.iterations
    );
    let (u_lqr, x_lqr) = lqr_reference(&p, &x0);
    // The first input is the receding-horizon control law.
    assert!(
        (sol.u_seq[0] - u_lqr[0]).abs() < 1e-6,
        "first input {} vs LQR {}",
        sol.u_seq[0],
        u_lqr[0]
    );
    assert_close(&sol.u_seq, &u_lqr, 1e-6, "u");
    assert_close(&sol.x_traj, &x_lqr, 1e-6, "x");
    // numpy: -H⁻¹ f of the condensed problem.
    let numpy = [
        -5.302298250661886,
        -2.288456667084613,
        -0.369716626041461,
        0.766906575646723,
    ];
    assert_close(&sol.u_seq[..4], &numpy, 1e-7, "numpy");
    let cost = total_cost(&p, &x0, &sol.u_seq);
    assert!(
        (sol.cost - cost).abs() <= 1e-12 * cost,
        "{} vs {cost}",
        sol.cost
    );
    assert!(sol.cost < total_cost(&p, &x0, &[0.0; 12]));
    assert_eq!(sol.x_traj, simulate(&p, &x0, &sol.u_seq));
}

#[test]
fn missing_nan_and_infinite_bounds_all_mean_unbounded() {
    let x0 = [1.0, -0.5];
    let infinite = solve(&reference_problem(12), &x0);
    let missing = MpcProblem {
        x_lb: vec![],
        x_ub: vec![],
        u_lb: vec![],
        u_ub: vec![],
        ..reference_problem(12)
    };
    let not_a_number = MpcProblem {
        x_lb: vec![f64::NAN, f64::NAN],
        x_ub: vec![f64::NAN],
        u_lb: vec![f64::NAN],
        u_ub: vec![f64::NAN],
        ..reference_problem(12)
    };
    for p in [missing, not_a_number] {
        let sol = solve(&p, &x0);
        assert!(sol.converged && sol.feasible);
        assert_eq!(sol.u_seq, infinite.u_seq);
        assert_eq!(sol.cost, infinite.cost);
    }
}

#[test]
fn asymmetric_weights_act_through_their_symmetric_part() {
    // xᵀ M x only sees (M + Mᵀ)/2: the gradient needs the symmetrised matrices.
    let skew = |m: &[f64]| -> Vec<f64> { vec![m[0], m[1] + 3.0, m[2] - 3.0, m[3]] };
    let sym = reference_problem(8);
    let asym = MpcProblem {
        q: skew(&sym.q),
        qf: skew(&sym.qf),
        ..sym.clone()
    };
    let x0 = [0.7, 0.2];
    let (a, s) = (solve(&asym, &x0), solve(&sym, &x0));
    assert!(a.converged && s.converged);
    assert_close(&a.u_seq, &s.u_seq, 1e-7, "u");
    assert!((a.cost - s.cost).abs() < 1e-9 * s.cost);
    let (u_lqr, _) = lqr_reference(&sym, &x0);
    assert_close(&a.u_seq, &u_lqr, 1e-6, "lqr");
}

#[test]
fn zero_initial_state_needs_no_iterations_to_certify_the_zero_input() {
    let sol = solve(&reference_problem(10), &[0.0, 0.0]);
    assert!(sol.converged && sol.feasible);
    assert!(sol.iterations <= 1, "{} iterations", sol.iterations);
    assert!(sol.u_seq.iter().all(|u| *u == 0.0));
    assert_eq!(sol.cost, 0.0);
}

// ── Step size: the fixed step of the previous implementation ────────────

#[test]
fn stiff_problem_that_made_the_fixed_step_diverge_converges() {
    // Hessian eigenvalue 5127 > 2 / 1e-3: the old fixed step 1e-3 doubled the
    // error in the stiffest direction every iteration (J ~ 1e21 after 500).
    let p = double_integrator(12, [1e4, 1e3], 1.0, [1e4, 1e3]);
    let x0 = [1.0, -0.5];
    let curvature = curvature_estimate(&p, 1.0, None).expect("consistent problem");
    assert!(
        curvature > 5000.0 && curvature < 5200.0,
        "curvature {curvature}"
    );
    let sol = solve(&p, &x0);
    assert!(sol.converged, "{} iterations", sol.iterations);
    assert!(sol.iterations < 2000, "{} iterations", sol.iterations);
    let (u_lqr, _) = lqr_reference(&p, &x0);
    assert_close(&sol.u_seq, &u_lqr, 1e-6, "u");
    // numpy optimum 40805.73897409487.
    assert!(
        (sol.cost - 40805.73897409487).abs() < 1e-9 * 40805.74,
        "cost {}",
        sol.cost
    );
}

#[test]
fn tiny_weights_that_stalled_the_fixed_step_reach_the_same_optimum() {
    // The same problem with every weight scaled by 1e-7 / 1e-6: the optimiser
    // is scale invariant, the old gradient step (1e-3 x 1e-6) never left the
    // origin.  numpy: u0 = -1.5017374458453392, J* = 6.150418061852913e-6.
    let p = double_integrator(12, [1e-6, 1e-7], 1e-7, [1e-6, 1e-7]);
    let x0 = [1.0, -0.5];
    let sol = solve(&p, &x0);
    assert!(sol.converged);
    let (u_lqr, _) = lqr_reference(&p, &x0);
    assert_close(&sol.u_seq, &u_lqr, 1e-7, "u");
    assert!(
        (sol.u_seq[0] + 1.5017374458453392).abs() < 1e-6,
        "u0 {}",
        sol.u_seq[0]
    );
    assert!(
        (sol.cost - 6.150418061852913e-6).abs() < 1e-12,
        "cost {}",
        sol.cost
    );
}

#[test]
fn curvature_estimate_matches_the_dense_hessian() {
    for p in [
        reference_problem(12),
        double_integrator(9, [1e3, 5.0], 0.3, [1e3, 5.0]),
    ] {
        let c = condense(&p, &[1.0, -0.5]);
        let eigenvalues = real_eigenvalues(&c.h, c.dim).expect("symmetric matrix");
        let largest = eigenvalues.iter().fold(0.0_f64, |m, e| m.max(e.0));
        let estimate = curvature_estimate(&p, 1.0, None).expect("consistent problem");
        assert!(
            estimate <= largest * (1.0 + 1e-9),
            "{estimate} above {largest}"
        );
        assert!(
            estimate >= largest * (1.0 - 1e-6),
            "{estimate} far below {largest}"
        );
    }
    // State-bound curvature: Γᵀ M Γ with M selecting the bounded velocity.
    let p = MpcProblem {
        x_ub: vec![INF, 1.0],
        ..reference_problem(12)
    };
    let c = condense(&p, &[0.0, 0.0]);
    let mut gram = vec![0.0; c.dim * c.dim];
    for k in 1..=p.horizon {
        let row = &c.g[(k * 2 + 1) * c.dim..(k * 2 + 2) * c.dim];
        for i in 0..c.dim {
            for j in 0..c.dim {
                gram[i * c.dim + j] += row[i] * row[j];
            }
        }
    }
    let eigenvalues = real_eigenvalues(&gram, c.dim).expect("symmetric matrix");
    let largest = eigenvalues.iter().fold(0.0_f64, |m, e| m.max(e.0));
    let estimate = curvature_estimate(&p, 0.0, Some(1.0)).expect("consistent problem");
    assert!(
        (estimate - largest).abs() <= 1e-6 * largest,
        "{estimate} vs {largest}"
    );
}

#[test]
fn adjoint_gradient_equals_central_differences_of_the_cost() {
    // Non-symmetric weights, three states and two inputs, an arbitrary point.
    let p = MpcProblem {
        horizon: 5,
        n: 3,
        m: 2,
        a: vec![0.9, 0.2, 0.0, -0.1, 1.0, 0.3, 0.05, 0.0, 0.8],
        b: vec![0.1, 0.0, 0.0, 0.2, 0.3, -0.1],
        q: vec![2.0, 0.5, 0.0, -0.3, 1.0, 0.2, 0.1, 0.0, 3.0],
        qf: vec![4.0, 0.0, 0.7, 0.2, 5.0, 0.0, 0.0, 0.4, 2.0],
        r: vec![1.0, 0.25, -0.25, 2.0],
        x_lb: vec![-INF; 3],
        x_ub: vec![INF; 3],
        u_lb: vec![-INF; 2],
        u_ub: vec![INF; 2],
    };
    let x0 = [0.4, -1.0, 0.8];
    let u: Vec<f64> = (0..10).map(|i| ((i as f64) * 1.3).sin() * 2.0).collect();
    let (value, gradient) = cost_and_gradient(&p, &x0, &u).expect("consistent problem");
    assert!((value - total_cost(&p, &x0, &u)).abs() <= 1e-13 * value.abs());
    let h = 1e-5;
    for i in 0..u.len() {
        let (mut plus, mut minus) = (u.clone(), u.clone());
        plus[i] += h;
        minus[i] -= h;
        let numeric = (total_cost(&p, &x0, &plus) - total_cost(&p, &x0, &minus)) / (2.0 * h);
        assert!(
            (gradient[i] - numeric).abs() <= 1e-7 * (1.0 + numeric.abs()),
            "d J / d u[{i}]: adjoint {} vs central difference {numeric}",
            gradient[i]
        );
    }
    // Wrong sizes are rejected instead of panicking.
    assert!(cost_and_gradient(&p, &x0[..2], &u).is_none());
    assert!(cost_and_gradient(&p, &x0, &u[..9]).is_none());
}

// ── Input bounds ────────────────────────────────────────────────────────

#[test]
fn active_input_bound_matches_the_enumerated_optimum() {
    let p = MpcProblem {
        u_lb: vec![-2.0],
        u_ub: vec![2.0],
        ..reference_problem(4)
    };
    let x0 = [1.0, -0.5];
    let sol = solve(&p, &x0);
    assert!(sol.converged && sol.feasible);
    // The bound is hit exactly (projection), not approximately.
    assert_eq!(sol.u_seq[0], -2.0);
    assert!(sol.u_seq.iter().all(|u| (-2.0..=2.0).contains(u)));
    assert_close(
        &sol.u_seq,
        &enumerated_optimum(&p, &x0),
        1e-6,
        "enumeration",
    );
    // numpy active-set reference.
    let numpy = [-2.0, -1.727722185928, -0.013787408165, 0.881080319161];
    assert_close(&sol.u_seq, &numpy, 1e-9, "numpy");
    assert!(
        (sol.cost - 46.073953072413055).abs() < 1e-9,
        "cost {}",
        sol.cost
    );
    // The bound costs something relative to the unconstrained optimum.
    assert!(sol.cost > solve(&reference_problem(4), &x0).cost + 0.1);
}

#[test]
fn input_bounds_are_never_violated_even_by_one_ulp() {
    let p = MpcProblem {
        u_lb: vec![-0.7],
        u_ub: vec![1.3],
        ..reference_problem(20)
    };
    for x0 in [[3.0, 1.0], [-3.0, -1.0], [0.2, 0.0]] {
        let sol = solve(&p, &x0);
        assert!(sol.converged);
        assert!(sol.u_seq.iter().all(|&u| (-0.7..=1.3).contains(&u)));
        assert_kkt_point(&p, &x0, &sol, 1e-6);
    }
}

// ── State bounds ────────────────────────────────────────────────────────

fn velocity_limited_problem() -> MpcProblem {
    MpcProblem {
        x_lb: vec![-INF, -0.4],
        u_lb: vec![-3.0],
        u_ub: vec![3.0],
        ..reference_problem(4)
    }
}

#[test]
fn active_state_and_input_bounds_match_the_enumerated_optimum() {
    let p = velocity_limited_problem();
    let x0 = [1.0, 0.0];
    // The unconstrained optimum breaks the velocity bound clearly.
    let free = solve(&reference_problem(4), &x0);
    let free_velocity = simulate(&reference_problem(4), &x0, &free.u_seq);
    let fastest = (1..=4)
        .map(|k| free_velocity[2 * k + 1])
        .fold(0.0_f64, f64::min);
    assert!(fastest < -0.7, "unconstrained velocity {fastest}");

    let sol = solve(&p, &x0);
    assert!(
        sol.converged && sol.feasible,
        "{} iterations",
        sol.iterations
    );
    let enumerated = enumerated_optimum(&p, &x0);
    assert_close(&sol.u_seq, &enumerated, 1e-6, "enumeration");
    // Hand computation: u0 hits the input bound, u1 brings the velocity to the
    // bound exactly, the later inputs keep it there.
    assert_close(&sol.u_seq, &[-3.0, -1.0, 0.0, 0.0], 1e-7, "hand");
    assert_eq!(sol.u_seq[0], -3.0);
    // J = 10.9 + 9.89225 + 9.185 + 8.441 + 15.458.
    assert!((sol.cost - 53.87625).abs() < 1e-7, "cost {}", sol.cost);
    // The bound is respected to the solver's feasibility tolerance.
    for k in 1..=4 {
        let velocity = sol.x_traj[2 * k + 1];
        assert!(velocity >= -0.4 - 1e-8, "stage {k} velocity {velocity}");
    }
    assert!(sol.cost > free.cost + 1.0);
    assert_kkt_point(&p, &x0, &sol, 1e-6);
}

#[test]
fn bounds_on_several_states_match_the_enumerated_optimum() {
    // Position floor, velocity floor and an input box: 16 constraints.
    let p = MpcProblem {
        x_lb: vec![0.93, -0.35],
        u_lb: vec![-3.0],
        u_ub: vec![3.0],
        ..reference_problem(4)
    };
    for x0 in [[1.0, 0.0], [1.05, -0.1]] {
        let sol = solve(&p, &x0);
        assert!(
            sol.converged && sol.feasible,
            "{x0:?}: {} iterations",
            sol.iterations
        );
        assert_close(
            &sol.u_seq,
            &enumerated_optimum(&p, &x0),
            1e-6,
            "enumeration",
        );
        assert_kkt_point(&p, &x0, &sol, 1e-6);
    }
}

#[test]
fn upper_state_bounds_work_too() {
    // Mirror image: pushing up, the velocity ceiling and a position ceiling.
    let p = MpcProblem {
        x_ub: vec![1.1, 0.4],
        u_lb: vec![-3.0],
        u_ub: vec![3.0],
        ..reference_problem(4)
    };
    let x0 = [-1.0, 0.0];
    let sol = solve(&p, &x0);
    assert!(sol.converged && sol.feasible);
    let mirrored = MpcProblem {
        x_lb: vec![-INF, -0.4],
        x_ub: vec![INF; 2],
        ..velocity_limited_problem()
    };
    let reference = solve(&mirrored, &[1.0, 0.0]);
    let negated: Vec<f64> = reference.u_seq.iter().map(|u| -u).collect();
    assert_close(&sol.u_seq, &negated, 1e-6, "mirror");
    assert_close(
        &sol.u_seq,
        &enumerated_optimum(&p, &x0),
        1e-6,
        "enumeration",
    );
}

#[test]
fn inactive_bounds_do_not_change_the_solution() {
    let x0 = [1.0, -0.5];
    let free = solve(&reference_problem(12), &x0);
    let loose = MpcProblem {
        x_lb: vec![-100.0, -100.0],
        x_ub: vec![100.0, 100.0],
        u_lb: vec![-100.0],
        u_ub: vec![100.0],
        ..reference_problem(12)
    };
    let sol = solve(&loose, &x0);
    assert!(sol.converged && sol.feasible);
    assert_close(&sol.u_seq, &free.u_seq, 1e-7, "u");
}

#[test]
fn long_horizon_with_mixed_bounds_satisfies_the_kkt_conditions() {
    let p = MpcProblem {
        x_lb: vec![-INF, -0.45],
        x_ub: vec![INF, 0.6],
        u_lb: vec![-2.5],
        u_ub: vec![2.5],
        ..double_integrator(30, [10.0, 1.0], 0.1, [40.0, 4.0])
    };
    let x0 = [2.0, 0.0];
    let sol = solve(&p, &x0);
    assert!(
        sol.converged && sol.feasible,
        "{} iterations",
        sol.iterations
    );
    assert_kkt_point(&p, &x0, &sol, 1e-6);
    // The velocity floor and the input floor are both reached.
    let velocities: Vec<f64> = (1..=30).map(|k| sol.x_traj[2 * k + 1]).collect();
    assert!(velocities.iter().any(|v| (v + 0.45).abs() < 1e-8));
    assert!(sol.u_seq.iter().any(|&u| u == -2.5));
}

#[test]
fn mimo_problem_with_input_and_state_bounds_satisfies_the_kkt_conditions() {
    // Two double integrators coupled through the inputs: n = 4, m = 2.
    let n = 4;
    let mut a = diagonal(&[1.0; 4]);
    a[1] = DT;
    a[2 * n + 3] = DT;
    a[n + 2] = 0.02;
    let b = vec![
        0.005, 0.0, //
        0.1, 0.01, //
        0.0, 0.005, //
        0.01, 0.1,
    ];
    let p = MpcProblem {
        horizon: 25,
        n,
        m: 2,
        a,
        b,
        q: diagonal(&[5.0, 0.5, 5.0, 0.5]),
        qf: diagonal(&[20.0, 2.0, 20.0, 2.0]),
        r: vec![0.2, 0.05, 0.05, 0.4],
        x_lb: vec![-INF, -0.5, -INF, -0.5],
        x_ub: vec![INF, 0.5, INF, 0.5],
        u_lb: vec![-2.0, -1.5],
        u_ub: vec![2.0, 1.5],
    };
    let x0 = [1.5, 0.0, -1.0, 0.3];
    let sol = solve(&p, &x0);
    assert!(
        sol.converged && sol.feasible,
        "{} iterations",
        sol.iterations
    );
    assert_eq!(sol.u_seq.len(), 50);
    for pair in sol.u_seq.chunks(2) {
        assert!((-2.0..=2.0).contains(&pair[0]) && (-1.5..=1.5).contains(&pair[1]));
    }
    // Both kinds of bounds are really in play.
    let speeds = (1..=25).flat_map(|k| [sol.x_traj[4 * k + 1].abs(), sol.x_traj[4 * k + 3].abs()]);
    assert!(speeds.fold(0.0_f64, f64::max) > 0.5 - 1e-7);
    assert!(
        sol.u_seq
            .chunks(2)
            .any(|pair| pair[0].abs() == 2.0 || pair[1].abs() == 1.5)
    );
    assert_kkt_point(&p, &x0, &sol, 1e-6);
    // Without bounds the same MIMO problem is the LQR solution.
    let free = MpcProblem {
        x_lb: vec![-INF; 4],
        x_ub: vec![INF; 4],
        u_lb: vec![-INF; 2],
        u_ub: vec![INF; 2],
        ..p.clone()
    };
    let unconstrained = solve(&free, &x0);
    assert!(unconstrained.converged);
    let (u_lqr, _) = lqr_reference(&free, &x0);
    assert_close(&unconstrained.u_seq, &u_lqr, 1e-6, "lqr");
    assert!(sol.cost > unconstrained.cost);
}

#[test]
fn receding_horizon_loop_keeps_the_velocity_inside_its_bounds() {
    let p = MpcProblem {
        x_lb: vec![-INF, -0.4],
        x_ub: vec![INF, 0.4],
        u_lb: vec![-3.0],
        u_ub: vec![3.0],
        ..double_integrator(10, [10.0, 1.0], 0.1, [20.0, 2.0])
    };
    let mpc = ModelPredictiveControl::new(p.clone());
    let mut x = vec![1.0, 0.0];
    let mut slowest = 0.0_f64;
    for step in 0..80 {
        let sol = mpc.solve(&x);
        assert!(
            sol.converged && sol.feasible,
            "step {step}: {} iterations",
            sol.iterations
        );
        let one_step = MpcProblem {
            horizon: 1,
            ..p.clone()
        };
        x = simulate(&one_step, &x, &sol.u_seq[..1])[2..].to_vec();
        slowest = slowest.min(x[1]);
        assert!(x[1].abs() <= 0.4 + 1e-8, "step {step}: velocity {}", x[1]);
    }
    // The limit is reached while travelling, and the cart ends up at rest.
    assert!((slowest + 0.4).abs() < 1e-6, "slowest velocity {slowest}");
    assert!(x[0].abs() < 0.02 && x[1].abs() < 0.02, "final state {x:?}");
}

// ── Infeasible and ill-posed problems ───────────────────────────────────

#[test]
fn unreachable_state_bound_reports_infeasibility_with_the_least_violation() {
    // The position must exceed 5 although |u| <= 1 moves the cart by 0.015.
    let p = MpcProblem {
        x_lb: vec![5.0, -INF],
        u_lb: vec![-1.0],
        u_ub: vec![1.0],
        ..double_integrator(3, [1.0, 1.0], 0.1, [1.0, 1.0])
    };
    let sol = solve(&p, &[0.0, 0.0]);
    assert!(!sol.feasible);
    assert!(!sol.converged);
    assert!(sol.cost.is_finite() && sol.x_traj.iter().all(|v| v.is_finite()));
    // Pushing as hard as allowed minimises every stage's violation.
    assert!(sol.u_seq.iter().all(|&u| (-1.0..=1.0).contains(&u)));
    for u in &sol.u_seq {
        assert!(*u > 1.0 - 1e-6, "input {u} is not at its upper bound");
    }
}

#[test]
fn initial_state_outside_the_bounds_is_reported_as_infeasible() {
    // The velocity 0.5 breaks the ceiling 0.45 at stage 0, which no input can
    // change; stages 1..=N can still be (and are) brought inside.
    let p = MpcProblem {
        x_ub: vec![INF, 0.45],
        u_lb: vec![-3.0],
        u_ub: vec![3.0],
        ..reference_problem(5)
    };
    let sol = solve(&p, &[0.9, 0.5]);
    assert!(!sol.feasible);
    for k in 1..=5 {
        assert!(sol.x_traj[2 * k + 1] <= 0.45 + 1e-8, "stage {k}");
    }
    assert!(sol.cost.is_finite());
}

#[test]
fn singular_hessian_is_handled_by_the_projected_gradient_method() {
    // Only the terminal position is weighted (Q = 0, R = 0, Qf = diag(1, 0)),
    // so H = 2 Γᵀ Qf Γ has rank 1 and the problem has a whole subspace of
    // minimisers.  The cart can be steered to the origin within the horizon,
    // so the optimal cost is zero; the old fixed step of 1e-3 did not get
    // anywhere near it in 500 iterations.
    let p = double_integrator(12, [0.0, 0.0], 0.0, [1.0, 0.0]);
    let x0 = [1.0, -0.5];
    let sol = solve(&p, &x0);
    assert!(sol.converged && sol.feasible);
    assert!(sol.cost < 1e-14, "cost {}", sol.cost);
    assert!(
        sol.x_traj[24].abs() < 1e-7,
        "final position {}",
        sol.x_traj[24]
    );
    assert!(sol.u_seq.iter().all(|u| u.is_finite()));
    // The same degenerate problem with an input box is solved to optimality
    // of the box-constrained problem: the final position is as close to zero
    // as the bounded inputs allow (here they do not allow zero).
    let bounded = MpcProblem {
        u_lb: vec![-0.5],
        u_ub: vec![0.5],
        ..p
    };
    let sol = solve(&bounded, &[3.0, 0.0]);
    assert!(sol.converged);
    assert!(sol.u_seq.iter().all(|&u| (-0.5..=0.5).contains(&u)));
    assert_kkt_point(&bounded, &[3.0, 0.0], &sol, 1e-6);
    assert!(sol.u_seq.iter().all(|&u| u == -0.5), "all inputs saturate");
}

// ── Degenerate and malformed input ──────────────────────────────────────

fn assert_failed(sol: &MpcSolution, p: &MpcProblem) {
    assert!(!sol.converged && !sol.feasible);
    assert!(sol.cost.is_nan());
    assert_eq!(sol.iterations, 0);
    assert_eq!(sol.u_seq, vec![0.0; p.horizon * p.m]);
    assert_eq!(sol.x_traj.len(), (p.horizon + 1) * p.n);
    assert!(sol.x_traj.iter().all(|v| v.is_nan()));
}

#[test]
fn malformed_input_gives_a_failed_solution_instead_of_a_panic() {
    let p = reference_problem(6);
    // x0 of the wrong length or with non-finite entries.
    assert_failed(&solve(&p, &[1.0]), &p);
    assert_failed(&solve(&p, &[1.0, 2.0, 3.0]), &p);
    assert_failed(&solve(&p, &[f64::NAN, 0.0]), &p);
    assert_failed(&solve(&p, &[0.0, INF]), &p);
    // Matrices that do not match the declared dimensions.
    for broken in [
        MpcProblem {
            a: vec![1.0, 0.0, 0.0],
            ..p.clone()
        },
        MpcProblem {
            b: vec![0.1],
            ..p.clone()
        },
        MpcProblem {
            q: vec![1.0; 3],
            ..p.clone()
        },
        MpcProblem {
            qf: vec![],
            ..p.clone()
        },
        MpcProblem {
            r: vec![0.1, 0.1],
            ..p.clone()
        },
    ] {
        assert_failed(&solve(&broken, &[1.0, 0.0]), &broken);
    }
    // Non-finite matrix entries.
    for broken in [
        MpcProblem {
            a: vec![1.0, f64::NAN, 0.0, 1.0],
            ..p.clone()
        },
        MpcProblem {
            b: vec![INF, 0.1],
            ..p.clone()
        },
        MpcProblem {
            r: vec![f64::NAN],
            ..p.clone()
        },
    ] {
        assert_failed(&solve(&broken, &[1.0, 0.0]), &broken);
    }
}

#[test]
fn zero_horizon_and_input_free_systems_are_evaluated_without_optimisation() {
    // N = 0: the "trajectory" is the initial state, the cost is x0ᵀ Qf x0.
    let p = reference_problem(0);
    let sol = solve(&p, &[1.0, -0.5]);
    assert!(sol.u_seq.is_empty() && sol.converged && sol.feasible);
    assert_eq!(sol.x_traj, vec![1.0, -0.5]);
    assert!((sol.cost - (20.0 + 0.5)).abs() < 1e-14);
    // A bound violated by x0 itself still makes the answer infeasible.
    let violated = MpcProblem {
        x_ub: vec![0.5, INF],
        ..p
    };
    assert!(!solve(&violated, &[1.0, -0.5]).feasible);

    // m = 0: nothing to optimise, the autonomous rollout is returned.
    let autonomous = MpcProblem {
        horizon: 3,
        n: 2,
        m: 0,
        a: vec![1.0, DT, 0.0, 1.0],
        b: vec![],
        q: diagonal(&[1.0, 1.0]),
        qf: diagonal(&[2.0, 2.0]),
        r: vec![],
        x_lb: vec![-INF; 2],
        x_ub: vec![INF; 2],
        u_lb: vec![],
        u_ub: vec![],
    };
    let sol = solve(&autonomous, &[1.0, -1.0]);
    assert!(sol.u_seq.is_empty() && sol.converged && sol.feasible);
    assert_eq!(sol.x_traj.len(), 8);
    assert_close(&sol.x_traj[6..], &[0.7, -1.0], 1e-14, "final state");
    let expected = total_cost(&autonomous, &[1.0, -1.0], &[]);
    assert!((sol.cost - expected).abs() < 1e-14 * expected);
}
