//! Tests of the shear-rate estimate and the self-consistent relaxation time of
//! the non-Newtonian BGK collision.
//!
//! Non-equilibrium populations of a prescribed strain-rate tensor are built from
//! the first-order Chapman-Enskog expansion
//! `f_neq_i = -(w_i rho tau / cs^2) Q_i : S`, for which the shear rate
//! `sqrt(2 S:S)` is known exactly.

use super::*;
use crate::lattice::{bgk_d2q9, equilibrium_d2q9, macros_from_d2q9, shear_rate_magnitude_d2q9};
use std::cell::Cell;

/// Chapman-Enskog non-equilibrium populations of a node with density `rho` and
/// relaxation time `tau` in the strain-rate field `S = [[sxx, sxy], [sxy, syy]]`.
fn fneq_from_strain(rho: f64, tau: f64, sxx: f64, syy: f64, sxy: f64) -> [f64; 9] {
    let mut fneq = [0.0_f64; 9];
    for i in 0..9 {
        let cx = C9_NN[i][0] as f64;
        let cy = C9_NN[i][1] as f64;
        let q = (cx * cx - CS2) * sxx + (cy * cy - CS2) * syy + 2.0 * cx * cy * sxy;
        fneq[i] = -W9_NN[i] * rho * tau / CS2 * q;
    }
    fneq
}

/// Populations of a node: equilibrium at `(rho, u)` plus `fneq`.
fn node_population(rho: f64, u: [f64; 2], fneq: &[f64; 9]) -> [f64; 9] {
    let feq = equilibrium_d2q9(rho, u[0], u[1]);
    let mut f = [0.0_f64; 9];
    for i in 0..9 {
        f[i] = feq[i] + fneq[i];
    }
    f
}

/// Relaxation frequency actually applied by `apply_non_newtonian_collision`,
/// read from the contraction of the largest non-equilibrium component.
fn applied_omega(model: &dyn NonNewtonianModel, rho: f64, u: [f64; 2], fneq: &[f64; 9]) -> f64 {
    let feq = equilibrium_d2q9(rho, u[0], u[1]);
    let mut f = node_population(rho, u, fneq);
    apply_non_newtonian_collision(&mut f, rho, u, model);
    let k = (0..9)
        .max_by(|&a, &b| fneq[a].abs().total_cmp(&fneq[b].abs()))
        .unwrap_or(0);
    1.0 - (f[k] - feq[k]) / fneq[k]
}

/// Viscosity linear in the shear rate, `nu = a + b gamma`.
struct LinearViscosity {
    a: f64,
    b: f64,
}

impl NonNewtonianModel for LinearViscosity {
    fn viscosity(&self, shear_rate: f64) -> f64 {
        self.a + self.b * shear_rate.abs()
    }
}

/// Newtonian viscosity that records the last shear rate it was queried with.
struct RecordingModel {
    nu: f64,
    last_shear: Cell<f64>,
}

impl RecordingModel {
    fn new(nu: f64) -> Self {
        Self {
            nu,
            last_shear: Cell::new(f64::NAN),
        }
    }
}

impl NonNewtonianModel for RecordingModel {
    fn viscosity(&self, shear_rate: f64) -> f64 {
        self.last_shear.set(shear_rate);
        self.nu
    }
}

impl LocalViscosityModel for RecordingModel {
    fn viscosity(&self, gamma_dot: f64) -> f64 {
        self.last_shear.set(gamma_dot);
        self.nu
    }
}

fn rel(value: f64, reference: f64) -> f64 {
    (value - reference).abs() / reference.abs().max(1e-300)
}

// ---------------------------------------------------------------------------
// Shear-rate estimate
// ---------------------------------------------------------------------------

#[test]
fn chapman_enskog_populations_reproduce_the_prescribed_shear_rate() {
    // Sanity check of the test fixture against the lattice helper.
    let (rho, tau) = (1.1, 0.8);
    for &(sxx, syy, sxy) in &[
        (0.0, 0.0, 0.002),
        (0.003, -0.003, 0.0),
        (0.001, -0.001, 0.0015),
    ] {
        let fneq = fneq_from_strain(rho, tau, sxx, syy, sxy);
        let expected = (2.0 * (sxx * sxx + syy * syy + 2.0 * sxy * sxy)).sqrt();
        let gamma = shear_rate_magnitude_d2q9(&fneq, rho, tau);
        assert!(rel(gamma, expected) < 1e-12, "{gamma} vs {expected}");
    }
}

#[test]
fn collision_estimates_shear_rate_from_all_stress_components() {
    // Constant viscosity nu = 0.1 -> tau = 0.8, so the solved tau equals the
    // tau the populations were generated with and the recovered shear rate is
    // the generating one.
    let nu = 0.1;
    let tau = 0.5 + nu / CS2;
    let rho = 1.0;
    for &(sxx, syy, sxy) in &[
        (0.0, 0.0, 0.002),    // simple shear
        (0.003, -0.003, 0.0), // planar extension: Pi_xx and Pi_yy only
        (0.001, -0.001, 0.0015),
    ] {
        let fneq = fneq_from_strain(rho, tau, sxx, syy, sxy);
        let expected = (2.0 * (sxx * sxx + syy * syy + 2.0 * sxy * sxy)).sqrt();

        let model = RecordingModel::new(nu);
        let mut f = node_population(rho, [0.02, -0.01], &fneq);
        apply_non_newtonian_collision(&mut f, rho, [0.02, -0.01], &model);
        // The model is queried at the converged shear rate (last call).
        let seen = model.last_shear.get();
        assert!(
            rel(seen, expected) < 1e-9,
            "collision: shear {seen} vs {expected}"
        );

        let model = RecordingModel::new(nu);
        let mut pop = [node_population(rho, [0.02, -0.01], &fneq)];
        non_newtonian_bgk_step_2d(&mut pop, 1, 1, &model);
        let seen = model.last_shear.get();
        assert!(
            rel(seen, expected) < 1e-9,
            "bgk step: shear {seen} vs {expected}"
        );
    }
}

#[test]
fn shear_rate_does_not_depend_on_the_bulk_velocity() {
    // The same non-equilibrium state at very different translation speeds must
    // relax identically: a uniformly translating fluid is not sheared.
    let model = PowerLaw::new(0.05, 0.6);
    let rho = 1.0;
    let fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, 0.004);
    let slow = applied_omega(&model, rho, [0.001, 0.0], &fneq);
    let fast = applied_omega(&model, rho, [0.08, -0.03], &fneq);
    let still = applied_omega(&model, rho, [0.0, 0.0], &fneq);
    assert!(rel(slow, fast) < 1e-9, "omega {slow} vs {fast}");
    assert!(rel(slow, still) < 1e-9, "omega {slow} vs {still}");
}

#[test]
fn uniform_translation_keeps_the_zero_shear_relaxation_time() {
    let rho = 1.2;
    let u = [0.05, 0.0];
    // lambda = 200: a shear rate of 0.01 is deep in the shear-thinning branch.
    let carreau = Carreau::new(0.3, 0.02, 200.0, 0.5);
    let tau = non_newtonian_relaxation_time_d2q9(&[0.0; 9], rho, |g| {
        NonNewtonianModel::viscosity(&carreau, g)
    });
    assert!(
        (tau - (0.5 + 0.3 / CS2)).abs() < 1e-12,
        "tau = {tau}, zero-shear value {}",
        0.5 + 0.3 / CS2
    );
    // The populations of the translating equilibrium are left untouched.
    let mut f = node_population(rho, u, &[0.0; 9]);
    let before = f;
    apply_non_newtonian_collision(&mut f, rho, u, &carreau);
    for i in 0..9 {
        assert!((f[i] - before[i]).abs() < 1e-15);
    }
    // A pure-shear node at rest is detected as sheared (old estimate: |u| = 0).
    let fneq = fneq_from_strain(rho, 0.9, 0.0, 0.0, 0.01);
    let sheared = non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| {
        NonNewtonianModel::viscosity(&carreau, g)
    });
    assert!(
        sheared < tau - 0.05,
        "shear thinning must lower tau: {sheared} vs {tau}"
    );
}

// ---------------------------------------------------------------------------
// Self-consistent relaxation time
// ---------------------------------------------------------------------------

/// `A = sqrt(2 Pi:Pi) / (2 rho cs^2)` of a simple-shear node, i.e. `gamma tau`.
fn shear_numerator(fneq: &[f64; 9], rho: f64) -> f64 {
    shear_rate_magnitude_d2q9(fneq, rho, 1.0)
}

#[test]
fn linear_viscosity_gives_the_quadratic_closed_form() {
    // nu = a + b gamma, gamma = A / tau:
    //   tau^2 - (1/2 + a/cs^2) tau - b A / cs^2 = 0
    let model = LinearViscosity { a: 0.08, b: 3.0 };
    let rho = 1.0;
    let fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, 0.005);
    let a_num = shear_numerator(&fneq, rho);
    let c0 = 0.5 + model.a / CS2;
    let expected = 0.5 * (c0 + (c0 * c0 + 4.0 * model.b * a_num / CS2).sqrt());
    let tau =
        non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| NonNewtonianModel::viscosity(&model, g));
    assert!(
        rel(tau, expected) < 1e-10,
        "tau = {tau}, expected {expected}"
    );

    // The collision applies exactly omega = 1/tau.
    let omega = applied_omega(&model, rho, [0.0, 0.0], &fneq);
    assert!(rel(omega, 1.0 / expected) < 1e-9, "omega = {omega}");
}

#[test]
fn power_law_with_exponent_one_half_matches_the_closed_form() {
    // nu = K gamma^(-1/2), gamma = A/tau  ->  with s = sqrt(tau):
    //   s^2 - C s - 1/2 = 0,   C = K A^(-1/2) / cs^2
    let model = PowerLaw::new(0.02, 0.5);
    let rho = 1.0;
    let fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, 0.01);
    let a_num = shear_numerator(&fneq, rho);
    let c = 0.02 * a_num.powf(-0.5) / CS2;
    let s = 0.5 * (c + (c * c + 2.0).sqrt());
    let expected = s * s;
    let tau =
        non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| NonNewtonianModel::viscosity(&model, g));
    assert!(
        rel(tau, expected) < 1e-10,
        "tau = {tau}, expected {expected}"
    );
    assert!(tau > 0.5 && tau.is_finite());
}

#[test]
fn strong_shear_thickening_is_solved_where_plain_iteration_diverges() {
    // n = 3: nu = k gamma^2, T(tau) = 1/2 + D / tau^2 with D = k A^2 / cs^2.
    // For D = 4 the fixed point is tau ~ 1.77 with |T'| = 2 - 1/tau > 1, so the
    // plain iteration oscillates away from it.
    let model = PowerLaw::new(16.0 / 3.0, 3.0);
    let rho = 1.0;
    let mut fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, 1.0);
    // Scale the populations so that A = 0.5.
    let a_num = shear_numerator(&fneq, rho);
    for v in fneq.iter_mut() {
        *v *= 0.5 / a_num;
    }
    let a_num = shear_numerator(&fneq, rho);
    assert!((a_num - 0.5).abs() < 1e-12);
    let d = 16.0 / 3.0 * a_num * a_num / CS2;
    let tau =
        non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| NonNewtonianModel::viscosity(&model, g));
    // tau^3 - tau^2 / 2 = D
    let residual = tau * tau * tau - 0.5 * tau * tau - d;
    assert!(
        residual.abs() < 1e-9 * d,
        "residual {residual}, tau = {tau}"
    );
    assert!(tau > 1.0 && tau < 3.0, "tau = {tau}");
}

#[test]
fn bingham_plug_and_yielded_branch_are_distinguished() {
    // gamma_y = tau_y / mu_p = 0.1.
    let model = BinghamPlastic::new(0.01, 0.1);
    let rho = 1.0;

    // Large stress: the fluid yields, tau = 1/2 + 3 (mu_p + tau_y tau / A)
    //   =>  tau = (1/2 + 3 mu_p) / (1 - 3 tau_y / A).
    let fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, 5.0);
    let a_num = shear_numerator(&fneq, rho);
    let expected = (0.5 + model.mu_p / CS2) / (1.0 - model.tau_y / (CS2 * a_num));
    let tau =
        non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| NonNewtonianModel::viscosity(&model, g));
    assert!(
        rel(tau, expected) < 1e-9,
        "tau = {tau}, expected {expected}"
    );
    assert!(
        a_num / tau > model.tau_y / model.mu_p,
        "the converged state must be above the yield shear rate"
    );

    // Small stress: rigid plug, the relaxation time hits the rigid viscosity.
    let fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, 0.001);
    let tau =
        non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| NonNewtonianModel::viscosity(&model, g));
    assert!(tau > 1.0e7, "un-yielded plug should be frozen, tau = {tau}");
}

#[test]
fn solver_returns_fixed_points_for_every_continuous_model() {
    let rho = 1.0;
    let hb = HerschelBulkley::new(0.05, 0.7, 0.002);
    let carreau = Carreau::new(0.4, 0.02, 80.0, 0.4);
    let reg = RegularizedBingham::new(0.002, 0.03, 500.0);
    let casson = CassonFluid::new(0.001, 0.03);
    let thin = PowerLaw::new(0.03, 0.4);
    let thick = PowerLaw::new(0.02, 1.8);
    let newtonian = PowerLaw::new(1.0 / 6.0, 1.0);
    let models: [(&str, &dyn Fn(f64) -> f64); 7] = [
        ("herschel-bulkley", &|g| {
            NonNewtonianModel::viscosity(&hb, g)
        }),
        ("carreau", &|g| NonNewtonianModel::viscosity(&carreau, g)),
        ("regularized bingham", &|g| reg.viscosity(g)),
        ("casson", &|g| casson.viscosity(g)),
        ("shear thinning", &|g| {
            NonNewtonianModel::viscosity(&thin, g)
        }),
        ("shear thickening", &|g| {
            NonNewtonianModel::viscosity(&thick, g)
        }),
        ("newtonian", &|g| {
            NonNewtonianModel::viscosity(&newtonian, g)
        }),
    ];
    for (name, nu) in models {
        for exponent in -6..=1 {
            let strain = 10.0_f64.powi(exponent);
            let fneq = fneq_from_strain(rho, 1.0, 0.0, 0.0, strain);
            let tau = non_newtonian_relaxation_time_d2q9(&fneq, rho, nu);
            assert!(
                tau.is_finite() && tau >= 0.5 + 1.0e-6,
                "{name}: tau = {tau} at S = {strain}"
            );
            let gamma = shear_rate_magnitude_d2q9(&fneq, rho, tau);
            let mapped = (0.5 + nu(gamma) / CS2).max(0.5 + 1.0e-6);
            assert!(
                (mapped - tau).abs() <= 1e-7 * tau,
                "{name}: tau = {tau} is not a fixed point (T = {mapped}) at S = {strain}"
            );
        }
    }
}

#[test]
fn collision_conserves_mass_and_momentum_for_every_model() {
    let rho = 1.3;
    let u = [0.04, -0.02];
    let fneq = fneq_from_strain(rho, 0.9, 0.001, -0.001, 0.002);
    let models: [&dyn NonNewtonianModel; 4] = [
        &PowerLaw::new(0.05, 0.6),
        &Carreau::new(0.3, 0.02, 40.0, 0.5),
        &HerschelBulkley::new(0.05, 0.8, 0.0005),
        &BinghamPlastic::new(0.0005, 0.08),
    ];
    for model in models {
        let mut f = node_population(rho, u, &fneq);
        apply_non_newtonian_collision(&mut f, rho, u, model);
        let mass: f64 = f.iter().sum();
        let mut mom = [0.0_f64; 2];
        for i in 0..9 {
            mom[0] += C9_NN[i][0] as f64 * f[i];
            mom[1] += C9_NN[i][1] as f64 * f[i];
        }
        assert!((mass - rho).abs() < 1e-13, "mass {mass}");
        assert!((mom[0] - rho * u[0]).abs() < 1e-13, "px {}", mom[0]);
        assert!((mom[1] - rho * u[1]).abs() < 1e-13, "py {}", mom[1]);
    }
}

#[test]
fn newtonian_model_reproduces_plain_bgk_on_the_grid() {
    let nu = 1.0 / 6.0;
    let omega = 1.0 / (0.5 + nu / CS2);
    let model = RecordingModel::new(nu);
    let (nx, ny) = (3, 2);
    let mut pop: Vec<[f64; 9]> = (0..nx * ny)
        .map(|k| {
            let fneq = fneq_from_strain(1.0 + 0.01 * k as f64, 0.9, 0.001, -0.001, 0.002);
            node_population(1.0 + 0.01 * k as f64, [0.01 * k as f64, -0.005], &fneq)
        })
        .collect();
    let mut reference = pop.clone();
    non_newtonian_bgk_step_2d(&mut pop, nx, ny, &model);
    for cell in reference.iter_mut() {
        let rho: f64 = cell.iter().sum();
        let mut mx = 0.0;
        let mut my = 0.0;
        for i in 0..9 {
            mx += C9_NN[i][0] as f64 * cell[i];
            my += C9_NN[i][1] as f64 * cell[i];
        }
        bgk_d2q9(cell, rho, mx / rho, my / rho, omega);
    }
    for (a, b) in pop.iter().zip(reference.iter()) {
        for i in 0..9 {
            assert!((a[i] - b[i]).abs() < 1e-14, "{} vs {}", a[i], b[i]);
        }
    }
}

// ---------------------------------------------------------------------------
// End-to-end channel flows
// ---------------------------------------------------------------------------

const EX: [i64; 9] = [0, 1, 0, -1, 0, 1, -1, -1, 1];
const EY: [i64; 9] = [0, 0, 1, 0, -1, 1, 1, -1, -1];
const OPP: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];

/// Minimal channel solver: periodic in x, halfway bounce-back walls below
/// (at rest) and above (moving with `u_top`), non-Newtonian BGK collision and an
/// optional body force added as a momentum impulse before the collision.
struct Channel {
    nx: usize,
    ny: usize,
    pop: Vec<[f64; 9]>,
}

impl Channel {
    fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            pop: vec![equilibrium_d2q9(1.0, 0.0, 0.0); nx * ny],
        }
    }

    fn step<M: LocalViscosityModel>(&mut self, model: &M, force: f64, u_top: f64) {
        for cell in self.pop.iter_mut() {
            for i in 0..9 {
                cell[i] += 3.0 * W9_NN[i] * force * C9_NN[i][0] as f64;
            }
        }
        non_newtonian_bgk_step_2d(&mut self.pop, self.nx, self.ny, model);
        let old = self.pop.clone();
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                let dst = iy * self.nx + ix;
                let rho: f64 = old[dst].iter().sum();
                for i in 0..9 {
                    let sy = iy as i64 - EY[i];
                    self.pop[dst][i] = if sy < 0 {
                        old[dst][OPP[i]]
                    } else if sy >= self.ny as i64 {
                        old[dst][OPP[i]] + 2.0 * W9_NN[i] * rho * EX[i] as f64 * u_top / CS2
                    } else {
                        let sx = (ix as i64 - EX[i]).rem_euclid(self.nx as i64) as usize;
                        old[sy as usize * self.nx + sx][i]
                    };
                }
            }
        }
    }

    fn node(&self, ix: usize, iy: usize) -> (f64, f64, f64) {
        macros_from_d2q9(&self.pop[iy * self.nx + ix])
    }

    fn ux(&self, ix: usize, iy: usize) -> f64 {
        self.node(ix, iy).1
    }
}

#[test]
fn power_law_couette_flow_recovers_the_imposed_shear_rate() {
    let (nx, ny) = (4, 16);
    let u_top = 0.04;
    let h = ny as f64;
    let gamma = u_top / h;
    // nu(gamma) = 0.1 at the imposed shear rate.
    let n_exp = 0.7;
    let k = 0.1 / gamma.powf(n_exp - 1.0);
    let model = PowerLawFluid::new(n_exp, k);
    let mut channel = Channel::new(nx, ny);
    for _ in 0..8_000 {
        channel.step(&model, 0.0, u_top);
    }
    // Linear profile u(y) = u_top (y + 1/2) / H between the halfway walls.
    for iy in 0..ny {
        let expected = u_top * (iy as f64 + 0.5) / h;
        let measured = channel.ux(1, iy);
        assert!(
            (measured - expected).abs() < 0.02 * u_top,
            "row {iy}: u = {measured}, expected {expected}"
        );
    }
    // Local shear rate from the non-equilibrium populations of a mid-channel
    // node, with the self-consistent relaxation time.
    let iy = ny / 2;
    let (rho, ux, uy) = channel.node(1, iy);
    let feq = equilibrium_d2q9(rho, ux, uy);
    let mut fneq = [0.0_f64; 9];
    for i in 0..9 {
        fneq[i] = channel.pop[iy * nx + 1][i] - feq[i];
    }
    let tau = non_newtonian_relaxation_time_d2q9(&fneq, rho, |g| model.viscosity(g));
    let gamma_local = shear_rate_magnitude_d2q9(&fneq, rho, tau);
    assert!(
        rel(gamma_local, gamma) < 0.02,
        "recovered shear rate {gamma_local} vs imposed {gamma}"
    );
    // The central-difference shear rate agrees.
    let gamma_fd = 0.5 * (channel.ux(1, iy + 1) - channel.ux(1, iy - 1));
    assert!(rel(gamma_fd, gamma) < 0.02, "fd shear {gamma_fd}");
}

#[test]
fn power_law_poiseuille_flow_matches_the_analytic_profile() {
    // Pressure-gradient-like body force F (per unit mass), kinematic power-law
    // fluid nu = K |gamma|^(n-1):  K |u'|^n sign(u') = -F (y - H/2),
    //   u(y) = n/(n+1) (F/K)^(1/n) [ (H/2)^((n+1)/n) - |y - H/2|^((n+1)/n) ].
    let (nx, ny) = (4, 24);
    let n_exp = 0.8;
    let k = 0.02;
    let force = 2.0e-6;
    let h = ny as f64;
    let model = PowerLawFluid::new(n_exp, k);
    let mut channel = Channel::new(nx, ny);
    for _ in 0..20_000 {
        channel.step(&model, force, 0.0);
    }
    let p = (n_exp + 1.0) / n_exp;
    let prefactor = n_exp / (n_exp + 1.0) * (force / k).powf(1.0 / n_exp);
    let analytic = |y: f64| prefactor * ((0.5 * h).powf(p) - (y - 0.5 * h).abs().powf(p));
    let u_max = analytic(0.5 * h);
    for iy in 0..ny {
        let y = iy as f64 + 0.5;
        let expected = analytic(y);
        let measured = channel.ux(1, iy);
        assert!(
            (measured - expected).abs() < 0.05 * u_max,
            "row {iy}: u = {measured}, analytic {expected} (u_max = {u_max})"
        );
    }
    // The profile is flatter than a parabola: the shear-thinning fluid is
    // plug-like at the centre.
    let newtonian_ratio = 0.75; // u(H/4) / u_max for a parabola
    let ratio = analytic(0.25 * h) / u_max;
    assert!(ratio > newtonian_ratio);
}
