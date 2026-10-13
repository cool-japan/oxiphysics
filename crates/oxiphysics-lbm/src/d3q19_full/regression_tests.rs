//! Regression tests for the audit fixes in [`super`]: the opposite-direction table, the MRT and
//! TRT collision operators, the on-node bounce-back walls and the Zou-He boundary conditions.

use super::*;
use crate::lattice::{D3Q19_OPPOSITES, D3Q19_VELOCITIES, D3Q19_WEIGHTS as LATTICE_WEIGHTS};

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

/// Mass and momentum `[rho, jx, jy, jz]` of the 19 populations starting at `base`.
fn moments(f: &[f64], base: usize) -> [f64; 4] {
    let mut m = [0.0_f64; 4];
    for i in 0..N_DIRS {
        let fi = f[base + i];
        m[0] += fi;
        m[1] += fi * f64::from(D3Q19_EX[i]);
        m[2] += fi * f64::from(D3Q19_EY[i]);
        m[3] += fi * f64::from(D3Q19_EZ[i]);
    }
    m
}

/// Index of the direction with the reversed velocity vector, found from the velocity components
/// (independent of the `D3Q19_OPP` table under test).
fn reversed_direction(i: usize) -> usize {
    (0..N_DIRS)
        .find(|&j| {
            D3Q19_EX[j] == -D3Q19_EX[i]
                && D3Q19_EY[j] == -D3Q19_EY[i]
                && D3Q19_EZ[j] == -D3Q19_EZ[i]
        })
        .expect("every D3Q19 velocity has a reversed partner")
}

/// Flat node index for `(x, y, z)`.
fn node(sim: &D3q19Simulation, x: usize, y: usize, z: usize) -> usize {
    x + y * sim.nx + z * sim.nx * sim.ny
}

/// Fill every node with a positive, randomly perturbed set of populations.
fn randomise_all_nodes(sim: &mut D3q19Simulation, seed: u64, amplitude: f64) {
    let mut rng = Lcg(seed);
    for k in 0..sim.nx * sim.ny * sim.nz {
        for (i, &weight) in D3Q19_WEIGHTS.iter().enumerate() {
            sim.f[k * N_DIRS + i] = weight * (1.0 + amplitude * rng.next_value());
        }
    }
}

/// A single-node simulation holding a randomly perturbed, moving, non-equilibrium state.
fn perturbed_single_node(tau: f64, seed: u64) -> D3q19Simulation {
    let mut sim = D3q19Simulation::new(1, 1, 1, tau);
    let mut rng = Lcg(seed);
    for (i, f_i) in sim.f.iter_mut().enumerate() {
        let feq = D3q19Simulation::equilibrium(1.1, [0.03, -0.02, 0.01], i);
        *f_i = feq + 0.05 * D3Q19_WEIGHTS[i] * rng.next_value();
    }
    sim
}

#[test]
fn opp_table_is_the_geometric_opposite_for_every_direction() {
    for i in 0..N_DIRS {
        assert_eq!(
            [D3Q19_EX[i], D3Q19_EY[i], D3Q19_EZ[i]],
            D3Q19_VELOCITIES[i],
            "velocity ordering differs from the lattice at direction {i}"
        );
        assert_eq!(D3Q19_WEIGHTS[i], LATTICE_WEIGHTS[i], "weight {i}");
        let j = D3Q19_OPP[i];
        assert_eq!(D3Q19_EX[j], -D3Q19_EX[i], "ex of opposite of {i}");
        assert_eq!(D3Q19_EY[j], -D3Q19_EY[i], "ey of opposite of {i}");
        assert_eq!(D3Q19_EZ[j], -D3Q19_EZ[i], "ez of opposite of {i}");
        assert_eq!(D3Q19_OPP[j], i, "opposite of the opposite of {i}");
    }
    assert_eq!(D3Q19_OPP, D3Q19_OPPOSITES);
    // The edge directions are the cases the old table got wrong.
    assert_eq!(D3Q19_OPP[7], 10);
    assert_eq!(D3Q19_OPP[8], 9);
    assert_eq!(D3Q19_OPP[11], 14);
    assert_eq!(D3Q19_OPP[12], 13);
    assert_eq!(D3Q19_OPP[15], 18);
    assert_eq!(D3Q19_OPP[16], 17);
}

#[test]
fn bounce_back_reverses_the_full_velocity_vector_at_solid_nodes() {
    let mut sim = D3q19Simulation::new(3, 3, 3, 1.0);
    let k = node(&sim, 1, 1, 1);
    sim.solid[k] = true;
    for i in 0..N_DIRS {
        sim.f[k * N_DIRS + i] = 1.0 + i as f64;
    }
    let other = node(&sim, 0, 0, 0);
    let other_before: Vec<f64> = sim.f[other * N_DIRS..(other + 1) * N_DIRS].to_vec();
    sim.apply_bounce_back();
    for i in 0..N_DIRS {
        let j = reversed_direction(i);
        assert_eq!(
            sim.f[k * N_DIRS + i],
            1.0 + j as f64,
            "slot {i} must hold the population that was travelling along -e_{i}"
        );
    }
    assert_eq!(
        &sim.f[other * N_DIRS..(other + 1) * N_DIRS],
        other_before.as_slice(),
        "fluid nodes must not be touched by bounce-back"
    );
}

#[test]
fn wall_bounce_back_y_swaps_every_opposite_pair_on_both_walls() {
    let (nx, ny, nz) = (3, 5, 2);
    let mut sim = D3q19Simulation::new(nx, ny, nz, 1.0);
    for (slot, value) in sim.f.iter_mut().enumerate() {
        *value = slot as f64;
    }
    let before = sim.f.clone();
    sim.apply_wall_bounce_back_y();
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                let base = node(&sim, x, y, z) * N_DIRS;
                for i in 0..N_DIRS {
                    let expected = if y == 0 || y == ny - 1 {
                        before[base + reversed_direction(i)]
                    } else {
                        before[base + i]
                    };
                    assert_eq!(
                        sim.f[base + i],
                        expected,
                        "(x={x}, y={y}, z={z}) direction {i}"
                    );
                }
            }
        }
    }
}

#[test]
fn wall_bounce_back_y_on_a_single_row_swaps_each_pair_once() {
    let mut sim = D3q19Simulation::new(2, 1, 1, 1.0);
    for (slot, value) in sim.f.iter_mut().enumerate() {
        *value = slot as f64;
    }
    let before = sim.f.clone();
    sim.apply_wall_bounce_back_y();
    for x in 0..2 {
        let base = node(&sim, x, 0, 0) * N_DIRS;
        for i in 0..N_DIRS {
            assert_eq!(
                sim.f[base + i],
                before[base + reversed_direction(i)],
                "direction {i}"
            );
        }
    }
}

/// Add the momentum `fx` per step to every fluid node as a direct source term
/// `3 w_i e_ix fx` (its first moment is exactly `fx`).
fn add_momentum_source(sim: &mut D3q19Simulation, fx: f64) {
    for k in 0..sim.nx * sim.ny * sim.nz {
        if sim.solid[k] {
            continue;
        }
        for i in 0..N_DIRS {
            sim.f[k * N_DIRS + i] += 3.0 * D3Q19_WEIGHTS[i] * f64::from(D3Q19_EX[i]) * fx;
        }
    }
}

#[test]
fn no_slip_channel_develops_a_poiseuille_profile() {
    // Solid walls at y = 0 and y = ny - 1 (on-node bounce-back), periodic in x and z, driven by
    // a momentum source. In steady state the on-node scheme places the no-slip wall half-way
    // between the wall row and the first fluid row, so the analytic profile is
    // u(y) = F / (2 nu) (y - 1/2) (ny - 3/2 - y). With the wrong opposite table the walls leak
    // momentum between each other and the profile is far from parabolic.
    let (nx, ny, nz) = (3, 11, 3);
    let fx = 1.0e-5;
    let mut sim = D3q19Simulation::new(nx, ny, nz, 1.0);
    for z in 0..nz {
        for x in 0..nx {
            let bottom = node(&sim, x, 0, z);
            let top = node(&sim, x, ny - 1, z);
            sim.solid[bottom] = true;
            sim.solid[top] = true;
        }
    }
    for _ in 0..800 {
        sim.collide();
        add_momentum_source(&mut sim, fx);
        sim.stream();
        sim.apply_bounce_back();
    }
    let nu = sim.kinematic_viscosity();
    let (y_lo, y_hi) = (0.5, ny as f64 - 1.5);
    let analytic = |y: f64| fx / (2.0 * nu) * (y - y_lo) * (y_hi - y);
    let u_max = analytic(0.5 * (y_lo + y_hi));
    for y in 1..ny - 1 {
        let (_, u) = sim.macroscopic(node(&sim, 1, y, 1));
        assert!(
            (u[0] - analytic(y as f64)).abs() < 0.015 * u_max,
            "y={y}: ux = {} vs analytic {}",
            u[0],
            analytic(y as f64)
        );
        assert!(u[1].abs() < 1e-9 && u[2].abs() < 1e-9, "y={y}: u = {u:?}");
        let (_, mirror) = sim.macroscopic(node(&sim, 1, ny - 1 - y, 1));
        assert!(
            (u[0] - mirror[0]).abs() < 1e-12,
            "profile symmetry at y={y}"
        );
    }
}

#[test]
fn trt_collision_conserves_mass_and_momentum_for_unequal_rates() {
    let (omega_plus, omega_minus) = (1.2, 1.7);
    for seed in 1..=8u64 {
        let mut sim = perturbed_single_node(0.9, seed);
        let before = moments(&sim.f, 0);
        sim.collide_trt(omega_plus, omega_minus);
        let after = moments(&sim.f, 0);
        for (c, name) in ["rho", "jx", "jy", "jz"].iter().enumerate() {
            assert!(
                (before[c] - after[c]).abs() < 1e-13,
                "seed {seed}: {name} changed from {} to {}",
                before[c],
                after[c]
            );
        }
    }
}

#[test]
fn trt_collision_matches_the_reference_operator_and_reduces_to_bgk() {
    let mut sim = perturbed_single_node(0.9, 42);
    let cell: [f64; N_DIRS] = std::array::from_fn(|i| sim.f[i]);
    let (rho, u) = sim.macroscopic(0);
    let reference = crate::collision::collide_trt(&cell, rho, u, 1.2, 1.7);
    sim.collide_trt(1.2, 1.7);
    for (i, (&got, &want)) in sim.f.iter().zip(reference.iter()).enumerate() {
        assert!(
            (got - want).abs() < 1e-14,
            "direction {i}: {got} vs reference {want}"
        );
    }

    // omega_plus == omega_minus is plain BGK with tau = 1 / omega.
    let omega = 1.3;
    let mut trt = perturbed_single_node(1.0 / omega, 7);
    let mut bgk = perturbed_single_node(1.0 / omega, 7);
    trt.collide_trt(omega, omega);
    bgk.collide();
    for i in 0..N_DIRS {
        assert!(
            (trt.f[i] - bgk.f[i]).abs() < 1e-14,
            "direction {i}: TRT {} vs BGK {}",
            trt.f[i],
            bgk.f[i]
        );
    }
}

#[test]
fn mrt_collision_conserves_mass_and_momentum_for_non_uniform_rates() {
    let s = crate::mrt3d::build_relaxation_matrix(1.1, 1.3, 1.2);
    for seed in 1..=8u64 {
        let mut sim = perturbed_single_node(0.9, seed);
        let before = moments(&sim.f, 0);
        sim.collide_mrt(&s);
        let after = moments(&sim.f, 0);
        for (c, name) in ["rho", "jx", "jy", "jz"].iter().enumerate() {
            assert!(
                (before[c] - after[c]).abs() < 1e-13,
                "seed {seed}: {name} changed from {} to {}",
                before[c],
                after[c]
            );
        }
    }
}

#[test]
fn mrt_collision_uses_the_moment_basis_not_per_direction_relaxation() {
    let s = crate::mrt3d::build_relaxation_matrix(1.1, 1.3, 1.2);
    let mut sim = perturbed_single_node(0.9, 11);
    let cell: [f64; N_DIRS] = std::array::from_fn(|i| sim.f[i]);
    let (rho, u) = sim.macroscopic(0);
    let reference = crate::mrt3d::collide_mrt(&cell, rho, u, &s);
    sim.collide_mrt(&s);
    for (i, (&got, &want)) in sim.f.iter().zip(reference.iter()).enumerate() {
        assert!(
            (got - want).abs() < 1e-14,
            "direction {i}: {got} vs reference {want}"
        );
    }

    // The xy shear stress is the moment `pxy` (index 13): it relaxes by exactly (1 - s[13]).
    let mut shear = perturbed_single_node(0.9, 11);
    let stress_before = shear.stress_tensor(0);
    shear.collide_mrt(&s);
    let stress_after = shear.stress_tensor(0);
    for (name, c) in [("xy", 3usize), ("xz", 4), ("yz", 5)] {
        assert!(
            (stress_after[c] - (1.0 - s[13]) * stress_before[c]).abs() < 1e-13,
            "{name} stress: {} vs {}",
            stress_after[c],
            (1.0 - s[13]) * stress_before[c]
        );
    }
    // The deviatoric normal stress 2 Pi_xx - Pi_yy - Pi_zz is the moment `3pxx` (index 9).
    let dev = |pi: &[f64; 6]| 2.0 * pi[0] - pi[1] - pi[2];
    assert!(
        (dev(&stress_after) - (1.0 - s[9]) * dev(&stress_before)).abs() < 1e-13,
        "normal-stress difference: {} vs {}",
        dev(&stress_after),
        (1.0 - s[9]) * dev(&stress_before)
    );
}

#[test]
fn mrt_collision_with_uniform_rates_equals_bgk() {
    let omega = 1.15;
    let s = [omega; 19];
    let mut mrt = perturbed_single_node(1.0 / omega, 5);
    let mut bgk = perturbed_single_node(1.0 / omega, 5);
    mrt.collide_mrt(&s);
    bgk.collide();
    for i in 0..N_DIRS {
        assert!(
            (mrt.f[i] - bgk.f[i]).abs() < 1e-13,
            "direction {i}: MRT {} vs BGK {}",
            mrt.f[i],
            bgk.f[i]
        );
    }
}

#[test]
fn mrt_collision_ignores_the_rates_of_conserved_moments() {
    let mut s = crate::mrt3d::build_relaxation_matrix(1.1, 1.3, 1.2);
    let mut a = perturbed_single_node(0.9, 21);
    a.collide_mrt(&s);
    for mode in [0usize, 3, 5, 7] {
        s[mode] = 1.5;
    }
    let mut b = perturbed_single_node(0.9, 21);
    let before = moments(&b.f, 0);
    b.collide_mrt(&s);
    let after = moments(&b.f, 0);
    for i in 0..N_DIRS {
        assert_eq!(a.f[i], b.f[i], "direction {i}");
    }
    for c in 0..4 {
        assert!((before[c] - after[c]).abs() < 1e-13, "moment {c}");
    }
}

#[test]
fn mrt_full_delegates_to_the_moment_space_operator() {
    let s = crate::mrt3d::build_relaxation_matrix(1.4, 1.2, 1.6);
    let mut a = perturbed_single_node(0.8, 3);
    let mut b = perturbed_single_node(0.8, 3);
    a.collide_mrt(&s);
    b.collide_mrt_full(&s);
    assert_eq!(a.f, b.f);
}

#[test]
fn mrt_and_trt_steps_conserve_mass_and_momentum_on_a_periodic_box() {
    let omega_minus = D3q19Simulation::optimal_omega_minus(1.2);
    let s = crate::mrt3d::build_relaxation_matrix(1.1, 1.3, 1.2);
    let mut trt = D3q19Simulation::new(4, 4, 4, 1.0 / 1.2);
    let mut mrt = D3q19Simulation::new(4, 4, 4, 1.0 / 1.3);
    randomise_all_nodes(&mut trt, 99, 0.1);
    randomise_all_nodes(&mut mrt, 99, 0.1);
    let total = |sim: &D3q19Simulation| {
        let mut m = [0.0_f64; 4];
        for k in 0..sim.nx * sim.ny * sim.nz {
            let nm = moments(&sim.f, k * N_DIRS);
            for c in 0..4 {
                m[c] += nm[c];
            }
        }
        m
    };
    let trt_before = total(&trt);
    let mrt_before = total(&mrt);
    for _ in 0..25 {
        trt.step_trt(1.2, omega_minus);
        mrt.step_mrt(&s);
    }
    let trt_after = total(&trt);
    let mrt_after = total(&mrt);
    for c in 0..4 {
        assert!(
            (trt_before[c] - trt_after[c]).abs() < 1e-10,
            "TRT moment {c}: {} -> {}",
            trt_before[c],
            trt_after[c]
        );
        assert!(
            (mrt_before[c] - mrt_after[c]).abs() < 1e-10,
            "MRT moment {c}: {} -> {}",
            mrt_before[c],
            mrt_after[c]
        );
    }
}

/// Face description used by the Zou-He tests: normal axis and whether it is the max-index face.
fn face_nodes(sim: &D3q19Simulation, axis: usize, at_max: bool) -> Vec<usize> {
    let dims = [sim.nx, sim.ny, sim.nz];
    let layer = if at_max { dims[axis] - 1 } else { 0 };
    let mut nodes = Vec::new();
    for z in 0..sim.nz {
        for y in 0..sim.ny {
            for x in 0..sim.nx {
                let coords = [x, y, z];
                if coords[axis] == layer {
                    nodes.push(node(sim, x, y, z));
                }
            }
        }
    }
    nodes
}

/// Direction indices whose velocity component along `axis` equals `value`.
fn directions_with(axis: usize, value: i32) -> Vec<usize> {
    (0..N_DIRS)
        .filter(|&i| direction_component(i, axis) == value)
        .collect()
}

#[test]
fn zou_he_inlet_x_imposes_the_prescribed_velocity_with_zero_transverse_flow() {
    let mut sim = D3q19Simulation::new(6, 5, 4, 0.9);
    randomise_all_nodes(&mut sim, 17, 0.2);
    let before = sim.f.clone();
    let u_in = 0.05;
    sim.apply_zou_he_inlet_x(u_in);
    let zero = directions_with(0, 0);
    let leaving = directions_with(0, -1);
    for k in face_nodes(&sim, 0, false) {
        let base = k * N_DIRS;
        let s0: f64 = zero.iter().map(|&i| before[base + i]).sum();
        let s_out: f64 = leaving.iter().map(|&i| before[base + i]).sum();
        let rho_expected = (s0 + 2.0 * s_out) / (1.0 - u_in);
        let (rho, u) = sim.macroscopic(k);
        assert!(
            (rho - rho_expected).abs() < 1e-13,
            "rho {rho} vs {rho_expected}"
        );
        assert!((u[0] - u_in).abs() < 1e-12, "ux = {}", u[0]);
        assert!(u[1].abs() < 1e-12, "uy = {}", u[1]);
        assert!(u[2].abs() < 1e-12, "uz = {}", u[2]);
        for i in zero.iter().chain(leaving.iter()) {
            assert_eq!(sim.f[base + i], before[base + i], "known population {i}");
        }
    }
    // Nodes off the inlet face are untouched.
    for k in 0..sim.nx * sim.ny * sim.nz {
        if k % sim.nx != 0 {
            assert_eq!(
                &sim.f[k * N_DIRS..(k + 1) * N_DIRS],
                &before[k * N_DIRS..(k + 1) * N_DIRS]
            );
        }
    }
}

#[test]
fn zou_he_outlet_x_imposes_the_density_and_allows_outflow() {
    let mut sim = D3q19Simulation::new(6, 5, 4, 0.9);
    randomise_all_nodes(&mut sim, 29, 0.2);
    let before = sim.f.clone();
    let rho_out = 1.03;
    sim.apply_zou_he_outlet_x(rho_out);
    let zero = directions_with(0, 0);
    let leaving = directions_with(0, 1);
    let mut saw_outflow = false;
    for k in face_nodes(&sim, 0, true) {
        let base = k * N_DIRS;
        let s0: f64 = zero.iter().map(|&i| before[base + i]).sum();
        let s_out: f64 = leaving.iter().map(|&i| before[base + i]).sum();
        let ux_expected = -1.0 + (s0 + 2.0 * s_out) / rho_out;
        let (rho, u) = sim.macroscopic(k);
        assert!((rho - rho_out).abs() < 1e-12, "rho = {rho}");
        assert!(
            (u[0] - ux_expected).abs() < 1e-12,
            "ux {} vs {ux_expected}",
            u[0]
        );
        assert!(u[1].abs() < 1e-12, "uy = {}", u[1]);
        assert!(u[2].abs() < 1e-12, "uz = {}", u[2]);
        saw_outflow |= u[0] > 1e-6;
        for i in zero.iter().chain(leaving.iter()) {
            assert_eq!(sim.f[base + i], before[base + i], "known population {i}");
        }
    }
    assert!(saw_outflow, "the outlet must be able to carry outflow");
}

#[test]
fn zou_he_outlet_x_leaves_a_matching_equilibrium_state_unchanged() {
    // An equilibrium state moving out through the outlet with the matching density is a fixed
    // point of the boundary condition; the old clamp forced u_x <= -1e-8 and broke it.
    let mut sim = D3q19Simulation::new(5, 4, 4, 0.8);
    sim.init_with_velocity(1.0, |_, _, _| [0.05, 0.0, 0.0]);
    let before = sim.f.clone();
    sim.apply_zou_he_outlet_x(1.0);
    for (slot, (&after, &prior)) in sim.f.iter().zip(before.iter()).enumerate() {
        assert!(
            (after - prior).abs() < 1e-14,
            "population {slot} changed from {prior} to {after}"
        );
    }
    let (_, u) = sim.macroscopic(node(&sim, 4, 1, 1));
    assert!((u[0] - 0.05).abs() < 1e-12, "outlet ux = {}", u[0]);
}

#[test]
fn zou_he_inlet_x_leaves_a_matching_equilibrium_state_unchanged() {
    let mut sim = D3q19Simulation::new(5, 4, 4, 0.8);
    sim.init_with_velocity(1.0, |_, _, _| [0.04, 0.0, 0.0]);
    let before = sim.f.clone();
    sim.apply_zou_he_inlet_x(0.04);
    for (slot, (&after, &prior)) in sim.f.iter().zip(before.iter()).enumerate() {
        assert!(
            (after - prior).abs() < 1e-14,
            "population {slot} changed from {prior} to {after}"
        );
    }
}

#[test]
fn zou_he_y_faces_impose_velocity_and_density() {
    let mut sim = D3q19Simulation::new(4, 6, 5, 0.9);
    randomise_all_nodes(&mut sim, 5, 0.2);
    let before = sim.f.clone();
    let uy_in = 0.04;
    sim.apply_zou_he_inlet_y(uy_in);
    let zero = directions_with(1, 0);
    let leaving_inlet = directions_with(1, -1);
    for k in face_nodes(&sim, 1, false) {
        let base = k * N_DIRS;
        let s0: f64 = zero.iter().map(|&i| before[base + i]).sum();
        let s_out: f64 = leaving_inlet.iter().map(|&i| before[base + i]).sum();
        let rho_expected = (s0 + 2.0 * s_out) / (1.0 - uy_in);
        let (rho, u) = sim.macroscopic(k);
        assert!(
            (rho - rho_expected).abs() < 1e-13,
            "rho {rho} vs {rho_expected}"
        );
        assert!(u[0].abs() < 1e-12, "ux = {}", u[0]);
        assert!((u[1] - uy_in).abs() < 1e-12, "uy = {}", u[1]);
        assert!(u[2].abs() < 1e-12, "uz = {}", u[2]);
    }

    let rho_out = 0.97;
    let before_outlet = sim.f.clone();
    sim.apply_zou_he_outlet_y(rho_out);
    let leaving_outlet = directions_with(1, 1);
    for k in face_nodes(&sim, 1, true) {
        let base = k * N_DIRS;
        let s0: f64 = zero.iter().map(|&i| before_outlet[base + i]).sum();
        let s_out: f64 = leaving_outlet
            .iter()
            .map(|&i| before_outlet[base + i])
            .sum();
        let uy_expected = -1.0 + (s0 + 2.0 * s_out) / rho_out;
        let (rho, u) = sim.macroscopic(k);
        assert!((rho - rho_out).abs() < 1e-12, "rho = {rho}");
        assert!(u[0].abs() < 1e-12, "ux = {}", u[0]);
        assert!(
            (u[1] - uy_expected).abs() < 1e-12,
            "uy {} vs {uy_expected}",
            u[1]
        );
        assert!(u[2].abs() < 1e-12, "uz = {}", u[2]);
    }
}

#[test]
fn zou_he_y_faces_leave_matching_equilibrium_states_unchanged() {
    let mut sim = D3q19Simulation::new(4, 5, 4, 0.8);
    sim.init_with_velocity(1.0, |_, _, _| [0.0, 0.03, 0.0]);
    let before = sim.f.clone();
    sim.apply_zou_he_inlet_y(0.03);
    sim.apply_zou_he_outlet_y(1.0);
    for (slot, (&after, &prior)) in sim.f.iter().zip(before.iter()).enumerate() {
        assert!(
            (after - prior).abs() < 1e-14,
            "population {slot} changed from {prior} to {after}"
        );
    }
}

#[test]
fn zou_he_skips_solid_nodes() {
    let mut sim = D3q19Simulation::new(5, 4, 4, 0.9);
    randomise_all_nodes(&mut sim, 61, 0.2);
    let solid_inlet = node(&sim, 0, 1, 2);
    let solid_outlet = node(&sim, 4, 2, 1);
    sim.solid[solid_inlet] = true;
    sim.solid[solid_outlet] = true;
    let before = sim.f.clone();
    sim.apply_zou_he_inlet_x(0.03);
    sim.apply_zou_he_outlet_x(1.0);
    for k in [solid_inlet, solid_outlet] {
        assert_eq!(
            &sim.f[k * N_DIRS..(k + 1) * N_DIRS],
            &before[k * N_DIRS..(k + 1) * N_DIRS],
            "solid node {k} must not be modified"
        );
    }
}

#[test]
fn zou_he_inlet_and_outlet_drive_a_three_dimensional_duct() {
    // Duct with no-slip walls on four sides, velocity inlet at x = 0 and pressure outlet at
    // x = nx - 1. The boundary conditions are applied after streaming and bounce-back.
    let (nx, ny, nz) = (24, 9, 9);
    let u_in = 0.02;
    let rho_out = 1.0;
    let mut sim = D3q19Simulation::new(nx, ny, nz, 0.8);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if y == 0 || y == ny - 1 || z == 0 || z == nz - 1 {
                    let k = node(&sim, x, y, z);
                    sim.solid[k] = true;
                }
            }
        }
    }
    for _ in 0..900 {
        sim.collide();
        sim.stream();
        sim.apply_bounce_back();
        sim.apply_zou_he_inlet_x(u_in);
        sim.apply_zou_he_outlet_x(rho_out);
    }
    assert!(sim.all_distributions_finite());

    // Prescribed inlet velocity and outlet density hold on the boundary nodes.
    for z in 1..nz - 1 {
        for y in 1..ny - 1 {
            let (_, u) = sim.macroscopic(node(&sim, 0, y, z));
            assert!((u[0] - u_in).abs() < 1e-10, "inlet ux = {}", u[0]);
            assert!(u[1].abs() < 1e-12 && u[2].abs() < 1e-12, "inlet u = {u:?}");
            let (rho, u) = sim.macroscopic(node(&sim, nx - 1, y, z));
            assert!((rho - rho_out).abs() < 1e-10, "outlet rho = {rho}");
            assert!(u[1].abs() < 1e-12 && u[2].abs() < 1e-12, "outlet u = {u:?}");
            assert!(u[0] > 0.0, "the outlet must carry outflow, ux = {}", u[0]);
        }
    }

    // Mass flux through x-planes is uniform in the developed region and the profile is
    // peaked on the axis.
    let flux = |x: usize| -> f64 {
        let mut q = 0.0;
        for z in 1..nz - 1 {
            for y in 1..ny - 1 {
                let (rho, u) = sim.macroscopic(node(&sim, x, y, z));
                q += rho * u[0];
            }
        }
        q
    };
    let q_ref = flux(nx / 2);
    assert!(q_ref > 0.0, "net flow must be in +x, got {q_ref}");
    for x in 3..nx - 3 {
        let rel = (flux(x) - q_ref).abs() / q_ref;
        assert!(rel < 0.02, "flux at x={x} deviates by {rel}");
    }
    let (_, u_axis) = sim.macroscopic(node(&sim, nx / 2, ny / 2, nz / 2));
    let (_, u_near_wall) = sim.macroscopic(node(&sim, nx / 2, 1, nz / 2));
    assert!(
        u_axis[0] > u_near_wall[0],
        "axis must be faster than the wall layer"
    );
    assert!(
        u_axis[0] > u_in,
        "developed axis velocity exceeds the mean inlet velocity"
    );
}
