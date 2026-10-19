//! Tests of the Guo forcing used by the porous-media driver: zeroth/first/second
//! force moments, mass conservation under a body force, Darcy and
//! Darcy-Forchheimer steady states, the half-force velocity correction and the
//! momentum-conserving equilibrium of partially porous cells.

use super::functions::{CX, CY, W, estimate_lattice_permeability, guo_forcing_term};
use super::types::*;

fn total_mass(grid: &PorousLbmGrid) -> f64 {
    grid.cells.iter().map(|c| c.f.iter().sum::<f64>()).sum()
}

fn total_momentum(grid: &PorousLbmGrid) -> [f64; 2] {
    let mut m = [0.0_f64; 2];
    for cell in &grid.cells {
        for i in 0..9 {
            m[0] += CX[i] * cell.f[i];
            m[1] += CY[i] * cell.f[i];
        }
    }
    m
}

fn rel(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1e-300)
}

// ---------------------------------------------------------------------------
// Moments of the forcing term
// ---------------------------------------------------------------------------

#[test]
fn guo_forcing_term_has_the_correct_moments_for_any_velocity() {
    let cases: [(f64, [f64; 2], [f64; 2]); 4] = [
        (0.75, [0.0, 0.0], [1.0e-5, 0.0]),
        (1.2, [0.03, -0.01], [2.0e-4, 5.0e-5]),
        (1.9, [-0.08, 0.05], [-3.0e-3, 1.0e-3]),
        (0.5, [0.2, 0.1], [0.0, -2.0e-2]),
    ];
    for &(omega, u, force) in &cases {
        let pref = 1.0 - 0.5 * omega;
        let mut m0 = 0.0;
        let mut m1 = [0.0_f64; 2];
        let mut m2 = [[0.0_f64; 2]; 2];
        for i in 0..9 {
            let fi = guo_forcing_term(i, omega, u, force);
            let c = [CX[i], CY[i]];
            m0 += fi;
            for a in 0..2 {
                m1[a] += c[a] * fi;
                for b in 0..2 {
                    m2[a][b] += c[a] * c[b] * fi;
                }
            }
        }
        let scale = force[0].abs() + force[1].abs();
        // Zeroth moment: the force creates no mass, whatever the velocity.
        assert!(
            m0.abs() < 1e-15 * scale.max(1e-300) * 10.0 + 1e-18,
            "m0 = {m0}"
        );
        // First moment: (1 - omega/2) F.
        for a in 0..2 {
            assert!(
                (m1[a] - pref * force[a]).abs() < 1e-15 + 1e-12 * scale,
                "m1[{a}]"
            );
        }
        // Second moment: (1 - omega/2)(u F + F u).
        for a in 0..2 {
            for b in 0..2 {
                let expected = pref * (u[a] * force[b] + u[b] * force[a]);
                assert!(
                    (m2[a][b] - expected).abs() < 1e-15 + 1e-12 * scale,
                    "m2[{a}][{b}] = {} expected {expected}",
                    m2[a][b]
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Mass conservation
// ---------------------------------------------------------------------------

#[test]
fn driver_conserves_total_mass_under_body_force_and_flow() {
    for &(eps, k) in &[(1.0, 2.0), (0.6, 2.0), (0.6, 1.0e-3), (0.4, 0.05)] {
        let mut drv = PorousMediaDriver::new(6, 5, 1.0 / 6.0, eps, k);
        drv.set_body_force(1.0e-5, 0.0);
        // Start from a moving fluid so that u.F is non-zero from the first step.
        for cell in drv.grid.cells.iter_mut() {
            for i in 0..9 {
                cell.f[i] = PorousLbmGrid::equilibrium(1.0, [0.01, 0.005], i);
            }
        }
        let m0 = total_mass(&drv.grid);
        drv.run(200);
        let m1 = total_mass(&drv.grid);
        assert!(
            (m1 - m0).abs() < 1e-12 * m0,
            "eps = {eps}, K = {k}: mass drifted {m0} -> {m1}"
        );
        for cell in &drv.grid.cells {
            assert!(cell.velocity[0].is_finite() && cell.velocity[1].is_finite());
        }
    }
}

// ---------------------------------------------------------------------------
// Darcy and Darcy-Forchheimer steady states
// ---------------------------------------------------------------------------

#[test]
fn darcy_velocity_equals_permeability_times_force_over_viscosity() {
    let nu = 1.0 / 6.0;
    let f = 1.0e-5;
    // The last case has nu/K = 167 per step: stable only because the drag is
    // treated implicitly.
    for &(eps, k) in &[(1.0, 2.0), (0.6, 2.0), (0.4, 0.5), (0.6, 1.0e-3)] {
        let mut drv = PorousMediaDriver::new(8, 8, nu, eps, k);
        drv.set_body_force(f, 0.0);
        drv.run(4_000);
        let u_sup = drv.mean_superficial_velocity();
        let darcy = k * f / nu;
        assert!(
            rel(u_sup[0], darcy) < 0.02,
            "eps = {eps}, K = {k}: superficial velocity {} vs K f / nu = {darcy}",
            u_sup[0]
        );
        // No transverse flow.  The y-momentum of a y-symmetric state is the
        // sum of populations of order 0.1, so double precision leaves a residue
        // of about 1e-17 whatever the (possibly tiny) Darcy velocity is; the
        // bound is that roundoff floor, not a fraction of the Darcy velocity.
        assert!(
            u_sup[1].abs() < 1e-14,
            "eps = {eps}, K = {k}: spurious transverse velocity {}",
            u_sup[1]
        );
        // The pore velocity is the superficial velocity over the porosity.
        let u_pore = drv.grid.cells[0].velocity[0];
        assert!(
            rel(eps * u_pore, darcy) < 0.02,
            "eps = {eps}, K = {k}: pore velocity {u_pore} vs K f / (nu eps)"
        );
    }
}

#[test]
fn forchheimer_drag_lowers_the_steady_flow() {
    // eps = 1, K = 1, nu = 1/6: Darcy coefficient 1/6, Forchheimer coefficient
    // C_F rho / sqrt(K) = 0.55.  Balance: 0.55 u^2 + u/6 = f.
    let nu = 1.0 / 6.0;
    let f = 5.0e-3;
    let mut with_inertia = PorousMediaDriver::new(8, 8, nu, 1.0, 1.0);
    with_inertia.set_body_force(f, 0.0);
    with_inertia.run(3_000);
    let u_f = with_inertia.mean_superficial_velocity()[0];
    let u_expected =
        (-(1.0 / 6.0) + ((1.0_f64 / 6.0).powi(2) + 4.0 * 0.55 * f).sqrt()) / (2.0 * 0.55);
    assert!(
        rel(u_f, u_expected) < 1e-3,
        "Darcy-Forchheimer velocity {u_f} vs {u_expected}"
    );

    let mut darcy_only = PorousMediaDriver::new(8, 8, nu, 1.0, 1.0);
    for model in darcy_only.dbf_models.iter_mut() {
        model.forchheimer_cf = 0.0;
    }
    darcy_only.set_body_force(f, 0.0);
    darcy_only.run(3_000);
    let u_d = darcy_only.mean_superficial_velocity()[0];
    assert!(rel(u_d, f * 6.0) < 1e-3, "Darcy velocity {u_d}");
    assert!(
        u_f < 0.95 * u_d,
        "the Forchheimer term must slow the flow: {u_f} vs {u_d}"
    );
}

// ---------------------------------------------------------------------------
// Half-force velocity correction
// ---------------------------------------------------------------------------

#[test]
fn velocity_includes_half_of_the_applied_force() {
    // Free flow (no drag): after n steps the lattice momentum is n F and the
    // physical velocity is (n F + F/2) / (rho eps).
    let (eps, f) = (0.5, 2.0e-4);
    let mut drv = PorousMediaDriver::new(4, 4, 1.0 / 6.0, eps, 1.0e30);
    drv.set_body_force(f, 0.0);
    drv.step();
    let u1 = drv.grid.cells[3].velocity[0];
    assert!(rel(u1, 1.5 * f / eps) < 1e-9, "u after one step: {u1}");
    drv.step();
    let u2 = drv.grid.cells[3].velocity[0];
    assert!(
        rel(u2 - u1, f / eps) < 1e-7,
        "velocity increment {}",
        u2 - u1
    );
    // The stored force is what the correction uses.
    assert!(rel(drv.grid.cells[3].force[0], f) < 1e-12);
}

#[test]
fn compute_macroscopic_adds_half_of_the_stored_force() {
    let mut grid = PorousLbmGrid::new(2, 2);
    let force = [4.0e-3, -2.0e-3];
    for cell in grid.cells.iter_mut() {
        cell.force = force;
    }
    grid.compute_macroscopic();
    for cell in &grid.cells {
        assert!(rel(cell.velocity[0], 0.5 * force[0]) < 1e-12);
        assert!(rel(cell.velocity[1], 0.5 * force[1]) < 1e-12);
    }
}

#[test]
fn half_step_velocity_satisfies_the_implicit_momentum_balance() {
    let dbf = DarcyBrinkmanForchheimer::new(0.02, 0.6, 1.0 / 6.0, 0.55, 1.0);
    let (rho, eps) = (1.02, 0.6);
    let momentum = [3.0e-3, -1.0e-3];
    let body_force = [2.0e-4, 5.0e-5];
    let u = dbf.half_step_velocity(momentum, body_force, rho, eps);
    let drag = dbf.total_force(u);
    for a in 0..2 {
        let lhs = rho * eps * u[a];
        let rhs = momentum[a] + 0.5 * (body_force[a] + drag[a]);
        assert!((lhs - rhs).abs() < 1e-15, "component {a}: {lhs} vs {rhs}");
    }
    // Without any force or momentum the fluid is at rest.
    let rest = dbf.half_step_velocity([0.0; 2], [0.0; 2], 1.0, 0.6);
    assert_eq!(rest, [0.0; 2]);
}

// ---------------------------------------------------------------------------
// estimate_lattice_permeability
// ---------------------------------------------------------------------------

#[test]
fn lattice_permeability_estimate_recovers_the_velocity_shift_drag_law() {
    // tau = 1 (nu = 1/6), resistance sigma = 0.05: the steady state of the
    // velocity-shift drag with Guo forcing is u = f (1 + sigma tau) / (rho sigma),
    // hence K_eff = nu u / f = nu (1 + sigma tau) / sigma = 3.5.
    let mut grid = PorousLbmGrid::new(6, 6);
    for cell in grid.cells.iter_mut() {
        cell.resistance = 0.05;
    }
    let m0 = total_mass(&grid);
    let k = estimate_lattice_permeability(&mut grid, 1.0, 1.0e-5, 800);
    let expected = (1.0 / 6.0) * 1.05 / 0.05;
    assert!(rel(k, expected) < 0.02, "K_eff = {k}, expected {expected}");
    let m1 = total_mass(&grid);
    assert!((m1 - m0).abs() < 1e-12 * m0, "mass drifted {m0} -> {m1}");
}

// ---------------------------------------------------------------------------
// Momentum-conserving equilibrium of partially porous cells
// ---------------------------------------------------------------------------

#[test]
fn collision_without_resistance_conserves_mass_and_momentum_for_partial_porosity() {
    let eps = 0.5;
    let mut grid = PorousLbmGrid::new(3, 3);
    grid.set_porosity_field(|_, _| eps);
    // Equilibrium at the superficial velocity plus a stress perturbation with
    // zero mass and momentum.
    let u_sup = [0.02, 0.01];
    for cell in grid.cells.iter_mut() {
        for i in 0..9 {
            cell.f[i] = PorousLbmGrid::equilibrium(1.0, u_sup, i)
                + 1.0e-3 * W[i] * (CX[i] * CX[i] - CY[i] * CY[i]);
        }
    }
    grid.compute_macroscopic();
    let m_before = total_momentum(&grid);
    let mass_before = total_mass(&grid);
    grid.collide(0.8);
    let m_after = total_momentum(&grid);
    let mass_after = total_mass(&grid);
    assert!((mass_after - mass_before).abs() < 1e-13);
    assert!(
        (m_after[0] - m_before[0]).abs() < 1e-13,
        "px {m_before:?} -> {m_after:?}"
    );
    assert!(
        (m_after[1] - m_before[1]).abs() < 1e-13,
        "py {m_before:?} -> {m_after:?}"
    );
    // The stored forces were reset: no external force acts in this collision.
    assert!(grid.cells.iter().all(|c| c.force == [0.0; 2]));
}
