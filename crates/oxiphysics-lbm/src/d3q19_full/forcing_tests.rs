//! Regression tests for the Guo body force of [`D3q19Simulation`].
//!
//! The previous implementation modified the populations *after* a collision that used the
//! uncorrected velocity, with the `(1 - 1/(2 tau))` prefactor but without the `F/2` velocity
//! shift. It therefore injected only `(1 - 1/(2 tau)) F` of momentum per step (`F/2` for
//! `tau = 1`), so a force-driven channel developed a velocity that was too small by exactly that
//! viscosity-dependent factor. These tests pin the exact momentum injection `F` per node and
//! step for every collision operator, the corrected velocity, the Guo momentum flux and a
//! force-driven Poiseuille channel against its analytic profile.

use super::*;

/// Deterministic pseudo-random generator (64-bit LCG) returning values in `[-1, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_value(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let unit = (self.0 >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * unit - 1.0
    }
}

/// Fill every node with a positive, randomly perturbed, moving, non-equilibrium state.
fn randomise(sim: &mut D3q19Simulation, seed: u64, amplitude: f64) {
    let mut rng = Lcg(seed);
    for k in 0..sim.nx * sim.ny * sim.nz {
        for i in 0..N_DIRS {
            let feq = D3q19Simulation::equilibrium(1.05, [0.02, -0.015, 0.01], i);
            sim.f[k * N_DIRS + i] = feq * (1.0 + amplitude * rng.next_value());
        }
    }
}

/// Mass and momentum `[rho, jx, jy, jz]` of node `k`.
fn node_moments(sim: &D3q19Simulation, k: usize) -> [f64; 4] {
    let (rho, momentum) = sim.raw_moments(k);
    [rho, momentum[0], momentum[1], momentum[2]]
}

/// Second moment `Pi_ab = sum_i f_i e_ia e_ib` of node `k`.
fn node_flux(sim: &D3q19Simulation, k: usize) -> [[f64; 3]; 3] {
    let mut pi = [[0.0_f64; 3]; 3];
    for i in 0..N_DIRS {
        let e = [
            f64::from(D3Q19_EX[i]),
            f64::from(D3Q19_EY[i]),
            f64::from(D3Q19_EZ[i]),
        ];
        let fi = sim.f[k * N_DIRS + i];
        for (a, row) in pi.iter_mut().enumerate() {
            for (b, entry) in row.iter_mut().enumerate() {
                *entry += fi * e[a] * e[b];
            }
        }
    }
    pi
}

/// Collision operators exercised against the same force.
#[derive(Clone, Copy, Debug)]
enum Operator {
    Bgk,
    Trt,
    Mrt,
}

fn collide_with(sim: &mut D3q19Simulation, operator: Operator) {
    let omega = 1.0 / sim.tau;
    match operator {
        Operator::Bgk => sim.collide(),
        Operator::Trt => sim.collide_trt(omega, D3q19Simulation::optimal_omega_minus(omega)),
        Operator::Mrt => {
            let s = crate::mrt3d::build_relaxation_matrix(1.3, omega, 1.1);
            sim.collide_mrt(&s);
        }
    }
}

const OPERATORS: [Operator; 3] = [Operator::Bgk, Operator::Trt, Operator::Mrt];

#[test]
fn every_collision_operator_injects_exactly_the_body_force_per_node_and_step() {
    let force = [2.0e-4, -1.0e-4, 3.0e-4];
    for &tau in &[0.55, 0.8, 1.0, 1.6] {
        for &operator in &OPERATORS {
            let mut sim = D3q19Simulation::new(3, 2, 3, tau);
            randomise(&mut sim, 7, 0.04);
            sim.apply_body_force(force[0], force[1], force[2]);
            let before: Vec<[f64; 4]> = (0..sim.f.len() / N_DIRS)
                .map(|k| node_moments(&sim, k))
                .collect();
            collide_with(&mut sim, operator);
            for (k, prior) in before.iter().enumerate() {
                let after = node_moments(&sim, k);
                assert!(
                    (after[0] - prior[0]).abs() < 1e-13,
                    "{operator:?}, tau {tau}, node {k}: mass changed"
                );
                for a in 0..3 {
                    assert!(
                        (after[a + 1] - prior[a + 1] - force[a]).abs() < 1e-13,
                        "{operator:?}, tau {tau}, node {k}: injected {} along axis {a}, \
                         expected {}",
                        after[a + 1] - prior[a + 1],
                        force[a]
                    );
                }
            }
        }
    }
}

#[test]
fn solid_nodes_receive_no_force() {
    let mut sim = D3q19Simulation::new(3, 3, 3, 0.9);
    randomise(&mut sim, 3, 0.03);
    let solid = sim.idx(1, 1, 1);
    sim.solid[solid] = true;
    sim.apply_body_force(1.0e-3, 2.0e-3, -1.0e-3);
    sim.apply_variable_body_force(|x, _y, _z| [1.0e-3 * x as f64, 0.0, 0.0]);
    let before = sim.f[solid * N_DIRS..(solid + 1) * N_DIRS].to_vec();
    assert_eq!(sim.force_field[solid], [0.0; 3]);
    sim.collide();
    assert_eq!(
        &sim.f[solid * N_DIRS..(solid + 1) * N_DIRS],
        before.as_slice(),
        "a solid node must be left untouched by collision and force"
    );
}

#[test]
fn variable_body_force_is_injected_node_by_node() {
    let mut sim = D3q19Simulation::new(4, 3, 2, 0.85);
    randomise(&mut sim, 21, 0.05);
    let force_at = |x: usize, y: usize, z: usize| {
        [
            1.0e-4 * (1.0 + x as f64),
            -2.0e-4 * y as f64,
            5.0e-5 * (z as f64 - 0.5),
        ]
    };
    sim.apply_variable_body_force(force_at);
    sim.apply_body_force(0.0, 0.0, 1.0e-5);
    let before: Vec<[f64; 4]> = (0..sim.f.len() / N_DIRS)
        .map(|k| node_moments(&sim, k))
        .collect();
    sim.collide();
    for z in 0..sim.nz {
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                let k = sim.idx(x, y, z);
                let after = node_moments(&sim, k);
                let expected = {
                    let f = force_at(x, y, z);
                    [f[0], f[1], f[2] + 1.0e-5]
                };
                for a in 0..3 {
                    assert!(
                        (after[a + 1] - before[k][a + 1] - expected[a]).abs() < 1e-13,
                        "node ({x},{y},{z}) axis {a}: injected {}, expected {}",
                        after[a + 1] - before[k][a + 1],
                        expected[a]
                    );
                }
            }
        }
    }
}

#[test]
fn reported_velocity_is_the_guo_corrected_one() {
    let (tau, fx) = (0.9, 4.0e-4);
    let mut sim = D3q19Simulation::new(4, 3, 3, tau);
    sim.apply_body_force(fx, 0.0, 0.0);
    let (rho, u) = sim.macroscopic(0);
    assert!((rho - 1.0).abs() < 1e-14);
    assert!(
        (u[0] - 0.5 * fx).abs() < 1e-15 && u[1].abs() < 1e-15 && u[2].abs() < 1e-15,
        "the rest state with a body force reports F/(2 rho): {u:?}"
    );
    let steps = 50;
    for _ in 0..steps {
        sim.step();
    }
    let expected = (steps as f64 + 0.5) * fx;
    for (k, velocity) in sim.velocity_field().iter().enumerate() {
        assert!(
            (velocity[0] - expected).abs() < 1e-12,
            "node {k}: ux = {}, expected (n + 1/2) F = {expected}",
            velocity[0]
        );
        assert!(velocity[1].abs() < 1e-12 && velocity[2].abs() < 1e-12);
    }
    assert!((sim.max_velocity() - expected).abs() < 1e-12);
    let n = (4 * 3 * 3) as f64;
    assert!(
        (sim.total_momentum_x() - n * expected).abs() < 1e-10,
        "total momentum {}",
        sim.total_momentum_x()
    );
}

/// Expected post-collision momentum flux of the Guo scheme started from a uniform equilibrium
/// state with density `rho`, velocity `u0` and force `F`, for the symmetric-part rate `omega`:
/// `Pi = rho cs^2 I + rho u0 u0 + (u0 F + F u0) + (1 - omega/4) F F / rho`.
fn expected_guo_flux(rho: f64, u0: [f64; 3], force: [f64; 3], omega: f64) -> [[f64; 3]; 3] {
    let mut pi = [[0.0_f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let delta = if a == b { 1.0 } else { 0.0 };
            pi[a][b] = rho * CS2 * delta
                + rho * u0[a] * u0[b]
                + (u0[a] * force[b] + force[a] * u0[b])
                + (1.0 - 0.25 * omega) * force[a] * force[b] / rho;
        }
    }
    pi
}

#[test]
fn guo_forcing_reproduces_the_exact_momentum_flux() {
    let (rho, u0) = (1.15, [0.04, -0.02, 0.03]);
    let force = [0.02, 0.01, -0.015];
    // BGK, and TRT with a different anti-symmetric rate: the momentum flux is an even moment,
    // so it depends on the symmetric rate only.
    for &(omega, omega_minus) in &[(0.8, 0.8), (1.3, 1.3), (1.6, 0.9), (0.7, 1.8)] {
        let mut sim = D3q19Simulation::new(1, 1, 1, 1.0 / omega);
        sim.init_with_velocity(rho, |_, _, _| u0);
        sim.apply_body_force(force[0], force[1], force[2]);
        if (omega - omega_minus).abs() < 1e-15 {
            sim.collide();
        } else {
            sim.collide_trt(omega, omega_minus);
        }
        let momentum = node_moments(&sim, 0);
        for a in 0..3 {
            assert!(
                (momentum[a + 1] - (rho * u0[a] + force[a])).abs() < 1e-13,
                "omega {omega}, omega_minus {omega_minus}: momentum[{a}] = {}",
                momentum[a + 1]
            );
        }
        let pi = node_flux(&sim, 0);
        let expected = expected_guo_flux(rho, u0, force, omega);
        for a in 0..3 {
            for b in 0..3 {
                assert!(
                    (pi[a][b] - expected[a][b]).abs() < 1e-13,
                    "omega {omega}, omega_minus {omega_minus}: Pi[{a}][{b}] = {}, expected {}",
                    pi[a][b],
                    expected[a][b]
                );
            }
        }
    }
}

#[test]
fn mrt_with_uniform_rates_equals_bgk_with_a_body_force() {
    let tau = 0.85;
    let mut bgk = D3q19Simulation::new(3, 2, 2, tau);
    randomise(&mut bgk, 99, 0.05);
    let mut mrt = bgk.clone();
    bgk.apply_body_force(3.0e-4, -1.0e-4, 2.0e-4);
    mrt.apply_body_force(3.0e-4, -1.0e-4, 2.0e-4);
    bgk.collide();
    mrt.collide_mrt(&[1.0 / tau; 19]);
    for (slot, (a, b)) in bgk.f.iter().zip(mrt.f.iter()).enumerate() {
        assert!(
            (a - b).abs() < 1e-12,
            "population {slot}: BGK {a} vs MRT {b}"
        );
    }
}

#[test]
fn trt_with_equal_rates_equals_bgk_with_a_body_force() {
    let tau = 1.2;
    let mut bgk = D3q19Simulation::new(2, 3, 2, tau);
    randomise(&mut bgk, 5, 0.05);
    let mut trt = bgk.clone();
    bgk.apply_body_force(-2.0e-4, 1.0e-4, 4.0e-4);
    trt.apply_body_force(-2.0e-4, 1.0e-4, 4.0e-4);
    bgk.collide();
    trt.collide_trt(1.0 / tau, 1.0 / tau);
    for (slot, (a, b)) in bgk.f.iter().zip(trt.f.iter()).enumerate() {
        assert!(
            (a - b).abs() < 1e-13,
            "population {slot}: BGK {a} vs TRT {b}"
        );
    }
}

#[test]
fn entropic_rescaling_conserves_momentum_of_a_forced_node() {
    let mut sim = D3q19Simulation::new(2, 2, 2, 0.9);
    randomise(&mut sim, 17, 0.05);
    sim.apply_body_force(5.0e-4, -3.0e-4, 1.0e-4);
    sim.collide();
    let before: Vec<[f64; 4]> = (0..8).map(|k| node_moments(&sim, k)).collect();
    sim.apply_entropic_stabilisation(1.4);
    for (k, prior) in before.iter().enumerate() {
        let after = node_moments(&sim, k);
        for c in 0..4 {
            assert!(
                (after[c] - prior[c]).abs() < 1e-13,
                "node {k}, moment {c}: {} -> {}",
                prior[c],
                after[c]
            );
        }
    }
}

#[test]
fn clear_body_force_removes_both_forces() {
    let mut sim = D3q19Simulation::new(2, 2, 2, 1.0);
    sim.apply_body_force(1.0e-3, 0.0, 0.0);
    sim.apply_variable_body_force(|_, _, _| [0.0, 2.0e-3, 0.0]);
    assert_eq!(sim.body_force, [1.0e-3, 0.0, 0.0]);
    assert_eq!(sim.force_field.len(), 8);
    sim.clear_body_force();
    assert_eq!(sim.body_force, [0.0; 3]);
    assert!(sim.force_field.is_empty());
    let (_, u) = sim.macroscopic(0);
    assert_eq!(u, [0.0; 3]);
}

#[test]
fn setup_poiseuille_channel_stores_its_driving_force() {
    let mut sim = D3q19Simulation::new(3, 9, 3, 0.8);
    sim.setup_poiseuille_channel(2.5e-5);
    assert_eq!(sim.body_force, [2.5e-5, 0.0, 0.0]);
    assert!(sim.solid[sim.idx(1, 0, 1)] && sim.solid[sim.idx(1, 8, 1)]);
}

/// Force-driven channel between on-node bounce-back walls at `y = 0` and `y = ny - 1`.
///
/// The on-node scheme places the no-slip planes half-way between the wall row and the first
/// fluid row, i.e. at `y = 1/2` and `y = ny - 3/2`, so the steady profile is
/// `u(y) = F / (2 nu) (y - 1/2)(ny - 3/2 - y)` with a centreline velocity `F H^2 / (8 nu)`
/// for the gap `H = ny - 2` (see [`D3q19Simulation::poiseuille_velocity`]).
fn run_channel(
    operator: Operator,
    ny: usize,
    tau: f64,
    fx: f64,
    steps: usize,
) -> (D3q19Simulation, f64) {
    let mut sim = D3q19Simulation::new(3, ny, 3, tau);
    sim.setup_poiseuille_channel(fx);
    let omega = 1.0 / tau;
    let omega_minus = D3q19Simulation::optimal_omega_minus(omega);
    for _ in 0..steps {
        match operator {
            Operator::Bgk => sim.step(),
            Operator::Trt => sim.step_trt(omega, omega_minus),
            Operator::Mrt => sim.step_mrt(&[omega; 19]),
        }
    }
    let nu = sim.kinematic_viscosity();
    let gap = (ny - 2) as f64;
    let u_max = fx * gap * gap / (8.0 * nu);
    (sim, u_max)
}

fn assert_poiseuille_profile(
    sim: &D3q19Simulation,
    fx: f64,
    u_max: f64,
    tolerance: f64,
    label: &str,
) {
    let nu = sim.kinematic_viscosity();
    let ny = sim.ny;
    for y in 1..ny - 1 {
        let analytic = D3q19Simulation::poiseuille_velocity(y as f64, ny, fx, nu);
        let (_, u) = sim.macroscopic(sim.idx(1, y, 1));
        assert!(
            (u[0] - analytic).abs() < tolerance * u_max,
            "{label}: y = {y}: ux = {}, analytic {analytic} (centreline {u_max})",
            u[0]
        );
        assert!(
            u[1].abs() < 1e-10 && u[2].abs() < 1e-10,
            "{label}: y = {y}: {u:?}"
        );
    }
}

#[test]
fn force_driven_poiseuille_flow_reaches_the_analytic_centreline_velocity() {
    // ny = 21, tau = 0.8: the old force injected only (1 - 1/(2 tau)) = 0.375 of F and gave a
    // centreline velocity of 0.375 u_max.
    let (ny, tau, fx) = (21, 0.8, 1.0e-5);
    let (sim, u_max) = run_channel(Operator::Bgk, ny, tau, fx, 6000);
    let centre = sim.macroscopic(sim.idx(1, ny / 2, 1)).1[0];
    let analytic_centre =
        D3q19Simulation::poiseuille_velocity(10.0, ny, fx, sim.kinematic_viscosity());
    assert!(
        (analytic_centre - u_max).abs() < 1e-12 * u_max,
        "the helper's centreline {analytic_centre} must equal F H^2 / (8 nu) = {u_max}"
    );
    assert!(
        (centre - u_max).abs() < 0.01 * u_max,
        "centreline velocity {centre}, analytic {u_max} (ratio {})",
        centre / u_max
    );
    assert_poiseuille_profile(&sim, fx, u_max, 0.005, "BGK");
}

#[test]
fn force_driven_poiseuille_flow_matches_with_trt_and_mrt() {
    let (ny, tau, fx) = (21, 0.8, 1.0e-5);
    // MRT with a uniform rate is BGK.
    let (sim, u_max) = run_channel(Operator::Mrt, ny, tau, fx, 6000);
    assert_poiseuille_profile(&sim, fx, u_max, 0.005, "MRT");
    // TRT with the magic parameter 3/16 places the half-way wall exactly, so the forced
    // channel is exactly parabolic: this pins the split of the Guo source into the
    // momentum (odd, relaxed with omega_minus) and momentum-flux (even, omega_plus) parts.
    let (sim, u_max) = run_channel(Operator::Trt, ny, tau, fx, 6000);
    assert_poiseuille_profile(&sim, fx, u_max, 1.0e-5, "TRT");
}

#[test]
fn poiseuille_throughput_scales_inversely_with_the_viscosity() {
    // The injected momentum no longer depends on tau, so the centreline velocity of a
    // force-driven channel is F H^2 / (8 nu) at both a low and a high viscosity.
    let (ny, fx) = (15, 2.0e-5);
    for &tau in &[0.7, 1.3] {
        let (sim, u_max) = run_channel(Operator::Bgk, ny, tau, fx, 4000);
        let centre = sim.macroscopic(sim.idx(1, ny / 2, 1)).1[0];
        assert!(
            (centre - u_max).abs() < 0.02 * u_max,
            "tau {tau}: centreline velocity {centre}, analytic {u_max}"
        );
    }
}

#[test]
fn poiseuille_helper_places_the_no_slip_planes_half_way_between_wall_and_fluid_rows() {
    let (ny, fx, nu) = (12, 3.0e-5, 0.2);
    // Zero on and outside the wall planes at y = 1/2 and y = ny - 3/2.
    for &y in &[-1.0, 0.0, 0.5, 10.5, 11.0, 20.0] {
        assert_eq!(
            D3q19Simulation::poiseuille_velocity(y, ny, fx, nu),
            0.0,
            "y = {y}"
        );
    }
    // Parabola with the centreline value F H^2 / (8 nu) for the gap H = ny - 2.
    let gap = (ny - 2) as f64;
    let centre = D3q19Simulation::poiseuille_velocity(5.5, ny, fx, nu);
    assert!((centre - fx * gap * gap / (8.0 * nu)).abs() < 1e-15);
    // Symmetric about the channel centre and a single fluid row between two solid rows.
    let near_bottom = D3q19Simulation::poiseuille_velocity(2.0, ny, fx, nu);
    let near_top = D3q19Simulation::poiseuille_velocity(9.0, ny, fx, nu);
    assert!((near_bottom - near_top).abs() < 1e-15);
    assert!((D3q19Simulation::poiseuille_velocity(1.0, 3, fx, nu) - fx / (8.0 * nu)).abs() < 1e-15);
    // No fluid row at all.
    assert_eq!(D3q19Simulation::poiseuille_velocity(1.0, 2, fx, nu), 0.0);
    assert_eq!(D3q19Simulation::poiseuille_velocity(0.0, 0, fx, nu), 0.0);
}
