// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests of the D3Q19 [`LBMGpuKernels`] CPU kernels: collision, streaming and
//! bounce-back are checked against independent references (push-form
//! streaming, literal-constant equilibrium), conservation laws, and the
//! analytic viscous decay of a shear wave.

use super::*;

// ── LBMGpuKernels tests ───────────────────────────────────────────────

#[test]
fn lbm_new_equilibrium_weights() {
    let lbm = LBMGpuKernels::new(4, 4, 4, 1.0, true, vec![]);
    // sum of weights over all Q directions at one cell should equal 1
    let rho = lbm.density_at(0, 0, 0);
    assert!((rho - 1.0).abs() < 1e-10);
}

/// Deterministic pseudo-random populations: the rest equilibrium with
/// every population perturbed by up to +-15 % (xorshift64, fixed seed).
fn perturbed_lbm(
    dims: [usize; 3],
    tau: f64,
    periodic: bool,
    solid_cells: Vec<usize>,
) -> LBMGpuKernels {
    let mut lbm = LBMGpuKernels::new(dims[0], dims[1], dims[2], tau, periodic, solid_cells);
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    for v in lbm.f.front.iter_mut() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let r = (state >> 11) as f64 / (1_u64 << 53) as f64;
        *v *= 1.0 + 0.3 * (r - 0.5);
    }
    lbm.f.back.clone_from(&lbm.f.front);
    lbm
}

/// Total mass and momentum over all cells.
fn lbm_totals(lbm: &LBMGpuKernels) -> (f64, [f64; 3]) {
    let (nx, ny, nz) = (lbm.f.width, lbm.f.height, lbm.f.depth);
    let mut mass = 0.0;
    let mut momentum = [0.0; 3];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                mass += lbm.density_at(x, y, z);
                let m = lbm.momentum_at(x, y, z);
                for d in 0..3 {
                    momentum[d] += m[d];
                }
            }
        }
    }
    (mass, momentum)
}

/// Equilibrium written out with the literal D3Q19 constants (3, 4.5, 1.5),
/// independent of `d3q19_equilibrium`.
fn feq_reference(i: usize, rho: f64, u: [f64; 3]) -> f64 {
    let e = [
        f64::from(D3Q19_EX[i]),
        f64::from(D3Q19_EY[i]),
        f64::from(D3Q19_EZ[i]),
    ];
    let eu = e[0] * u[0] + e[1] * u[1] + e[2] * u[2];
    let uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    D3Q19_WEIGHTS[i] * rho * (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * uu)
}

#[test]
fn lbm_collision_conserves_mass_and_momentum_per_cell() {
    let dims = [4, 3, 5];
    let mut lbm = perturbed_lbm(dims, 0.8, false, vec![]);
    let before = lbm.clone();
    lbm.dispatch_collision();
    let mut max_change = 0.0_f64;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let (rho0, rho1) = (before.density_at(x, y, z), lbm.density_at(x, y, z));
                assert!(
                    (rho0 - rho1).abs() < 1e-12,
                    "cell ({x},{y},{z}): rho {rho0} -> {rho1}"
                );
                let (m0, m1) = (before.momentum_at(x, y, z), lbm.momentum_at(x, y, z));
                for d in 0..3 {
                    assert!(
                        (m0[d] - m1[d]).abs() < 1e-12,
                        "cell ({x},{y},{z}): momentum[{d}] {} -> {}",
                        m0[d],
                        m1[d]
                    );
                }
                for c in 0..D3Q19_Q {
                    max_change =
                        max_change.max((before.f.read(x, y, z, c) - lbm.f.read(x, y, z, c)).abs());
                }
            }
        }
    }
    // The perturbed populations are off equilibrium, so collision acts.
    assert!(max_change > 1e-3, "collision left the state unchanged");
}

#[test]
fn lbm_collision_with_unit_tau_yields_the_equilibrium() {
    let mut lbm = perturbed_lbm([3, 3, 3], 1.0, false, vec![]);
    let before = lbm.clone();
    lbm.dispatch_collision();
    for z in 0..3 {
        for y in 0..3 {
            for x in 0..3 {
                let rho = before.density_at(x, y, z);
                let m = before.momentum_at(x, y, z);
                let u = [m[0] / rho, m[1] / rho, m[2] / rho];
                for c in 0..D3Q19_Q {
                    let expected = feq_reference(c, rho, u);
                    let got = lbm.f.read(x, y, z, c);
                    assert!(
                        (got - expected).abs() < 1e-13,
                        "cell ({x},{y},{z}) direction {c}: {got} vs {expected}"
                    );
                }
            }
        }
    }
}

/// A uniform equilibrium flow is a fixed point of collision and (periodic)
/// streaming, and the velocity recovered from the populations is the
/// imposed one.
#[test]
fn lbm_uniform_flow_is_a_fixed_point() {
    let (rho, u) = (1.2, [0.05, -0.02, 0.03]);
    let mut lbm = LBMGpuKernels::new(4, 3, 5, 0.7, true, vec![]);
    for cell in 0..4 * 3 * 5 {
        for c in 0..D3Q19_Q {
            lbm.f.front[cell * D3Q19_Q + c] = feq_reference(c, rho, u);
        }
    }
    lbm.f.back.clone_from(&lbm.f.front);
    let initial = lbm.f.front.clone();
    for _ in 0..5 {
        lbm.dispatch_collision();
        lbm.dispatch_streaming();
        lbm.dispatch_boundary();
    }
    for (a, b) in lbm.f.front.iter().zip(initial.iter()) {
        assert!((a - b).abs() < 1e-13);
    }
    let m = lbm.momentum_at(1, 1, 2);
    let rho_c = lbm.density_at(1, 1, 2);
    assert!((rho_c - rho).abs() < 1e-13);
    for d in 0..3 {
        assert!((m[d] / rho_c - u[d]).abs() < 1e-13, "velocity[{d}]");
    }
}

#[test]
fn lbm_collision_skips_solid_cells() {
    let mut lbm = perturbed_lbm([3, 3, 3], 0.8, false, vec![13]);
    let before = lbm.clone();
    lbm.dispatch_collision();
    // Cell 13 is (1, 1, 1): untouched; its neighbour (0, 1, 1) relaxed.
    for c in 0..D3Q19_Q {
        assert_eq!(lbm.f.read(1, 1, 1, c), before.f.read(1, 1, 1, c));
    }
    let moved: f64 = (0..D3Q19_Q)
        .map(|c| (lbm.f.read(0, 1, 1, c) - before.f.read(0, 1, 1, c)).abs())
        .sum();
    assert!(moved > 1e-6, "fluid cells must still collide");
}

/// Push-form reference streaming: every population moves to `x + e_i`;
/// one that would leave a non-periodic domain is reflected in place into
/// direction `opp(i)`.  Independent of the pull implementation.
fn reference_stream(lbm: &LBMGpuKernels) -> Vec<f64> {
    let (nx, ny, nz) = (lbm.f.width, lbm.f.height, lbm.f.depth);
    let dims = [nx as isize, ny as isize, nz as isize];
    let mut out = vec![f64::NAN; lbm.f.front.len()];
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                for i in 0..D3Q19_Q {
                    let mut dst = [
                        x as isize + D3Q19_EX[i] as isize,
                        y as isize + D3Q19_EY[i] as isize,
                        z as isize + D3Q19_EZ[i] as isize,
                    ];
                    let mut outside = false;
                    for (d, &n) in dst.iter_mut().zip(dims.iter()) {
                        if lbm.streaming.periodic {
                            *d = d.rem_euclid(n);
                        } else if *d < 0 || *d >= n {
                            outside = true;
                        }
                    }
                    let value = lbm.f.read(x, y, z, i);
                    if outside {
                        let at = lbm.f.idx(x, y, z, D3Q19_OPP[i]);
                        assert!(out[at].is_nan(), "population written twice");
                        out[at] = value;
                    } else {
                        let at = lbm
                            .f
                            .idx(dst[0] as usize, dst[1] as usize, dst[2] as usize, i);
                        assert!(out[at].is_nan(), "population written twice");
                        out[at] = value;
                    }
                }
            }
        }
    }
    out
}

#[test]
fn lbm_streaming_matches_push_reference_periodic() {
    let mut lbm = perturbed_lbm([4, 3, 5], 1.0, true, vec![]);
    let expected = reference_stream(&lbm);
    assert!(expected.iter().all(|v| v.is_finite()));
    lbm.dispatch_streaming();
    assert_eq!(lbm.f.front, expected);
}

#[test]
fn lbm_streaming_matches_push_reference_with_walls() {
    let mut lbm = perturbed_lbm([4, 3, 5], 1.0, false, vec![]);
    let expected = reference_stream(&lbm);
    assert!(expected.iter().all(|v| v.is_finite()));
    lbm.dispatch_streaming();
    assert_eq!(lbm.f.front, expected);
}

/// Every one of the 19 populations travels along its own velocity vector
/// (with periodic wrap-around).
#[test]
fn lbm_streaming_moves_each_population_along_its_velocity() {
    let dims = [4_usize, 5, 6];
    let origin = [3_usize, 0, 2];
    let mut lbm = LBMGpuKernels::new(dims[0], dims[1], dims[2], 1.0, true, vec![]);
    lbm.f.front.fill(0.0);
    for i in 0..D3Q19_Q {
        let at = lbm.f.idx(origin[0], origin[1], origin[2], i);
        lbm.f.front[at] = 100.0 + i as f64;
    }
    lbm.dispatch_streaming();
    for i in 0..D3Q19_Q {
        let target = [
            (origin[0] as isize + D3Q19_EX[i] as isize).rem_euclid(dims[0] as isize) as usize,
            (origin[1] as isize + D3Q19_EY[i] as isize).rem_euclid(dims[1] as isize) as usize,
            (origin[2] as isize + D3Q19_EZ[i] as isize).rem_euclid(dims[2] as isize) as usize,
        ];
        assert_eq!(
            lbm.f.read(target[0], target[1], target[2], i),
            100.0 + i as f64,
            "direction {i} did not arrive at {target:?}"
        );
    }
    assert_eq!(lbm.f.front.iter().filter(|&&v| v != 0.0).count(), D3Q19_Q);
}

#[test]
fn lbm_streaming_conserves_mass_and_momentum_when_periodic() {
    let mut lbm = perturbed_lbm([4, 3, 5], 1.0, true, vec![]);
    let (mass0, mom0) = lbm_totals(&lbm);
    assert!(
        mom0.iter().any(|m| m.abs() > 1e-3),
        "need non-zero momentum"
    );
    lbm.dispatch_streaming();
    let (mass1, mom1) = lbm_totals(&lbm);
    assert!((mass0 - mass1).abs() < 1e-12);
    for d in 0..3 {
        assert!((mom0[d] - mom1[d]).abs() < 1e-12, "momentum[{d}]");
    }
}

/// Mass and momentum are conserved by a full collide + stream cycle on a
/// periodic domain without obstacles.
#[test]
fn lbm_collide_and_stream_conserve_mass_and_momentum() {
    let mut lbm = perturbed_lbm([5, 4, 3], 0.9, true, vec![]);
    let (mass0, mom0) = lbm_totals(&lbm);
    for _ in 0..25 {
        lbm.dispatch_collision();
        lbm.dispatch_streaming();
        lbm.dispatch_boundary();
    }
    let (mass1, mom1) = lbm_totals(&lbm);
    assert!(
        (mass0 - mass1).abs() < 1e-10,
        "mass drifted: {mass0} -> {mass1}"
    );
    for d in 0..3 {
        assert!(
            (mom0[d] - mom1[d]).abs() < 1e-10,
            "momentum[{d}] drifted: {} -> {}",
            mom0[d],
            mom1[d]
        );
    }
}

/// A population leaving a non-periodic domain through a face is reflected
/// back into the opposite direction of the same cell.
#[test]
fn lbm_streaming_reflects_populations_at_domain_walls() {
    let mut lbm = LBMGpuKernels::new(4, 3, 2, 1.0, false, vec![]);
    lbm.f.front.fill(0.0);
    // +x population in the last x-layer wants to leave the domain.
    let leaving = lbm.f.idx(3, 1, 0, 1);
    lbm.f.front[leaving] = 7.0;
    lbm.dispatch_streaming();
    assert_eq!(lbm.f.read(3, 1, 0, 2), 7.0, "reflected into -x");
    assert_eq!(lbm.f.front.iter().filter(|&&v| v != 0.0).count(), 1);
}

#[test]
fn lbm_boundary_reverses_populations_of_solid_cells() {
    // Cell 5 is (1, 1, 0) in a 4 x 4 x 1 grid.
    let mut lbm = LBMGpuKernels::new(4, 4, 1, 1.0, false, vec![5]);
    lbm.f.front.fill(0.0);
    let at = lbm.f.idx(1, 1, 0, 1);
    lbm.f.front[at] = 1.0;
    let other = lbm.f.idx(2, 1, 0, 1);
    lbm.f.front[other] = 3.0;
    lbm.dispatch_boundary();
    assert_eq!(lbm.f.read(1, 1, 0, 2), 1.0, "f_1 = 1 must end up in f_2");
    assert_eq!(lbm.f.read(1, 1, 0, 1), 0.0);
    assert_eq!(lbm.f.read(2, 1, 0, 1), 3.0, "fluid cells are untouched");
}

#[test]
fn lbm_boundary_reflects_every_direction_pair() {
    let mut lbm = perturbed_lbm([3, 3, 3], 1.0, false, vec![13]);
    let before = lbm.clone();
    lbm.dispatch_boundary();
    for (c, &opposite) in D3Q19_OPP.iter().enumerate() {
        assert_eq!(
            lbm.f.read(1, 1, 1, c),
            before.f.read(1, 1, 1, opposite),
            "direction {c}"
        );
    }
    for c in 0..D3Q19_Q {
        assert_eq!(lbm.f.read(0, 0, 0, c), before.f.read(0, 0, 0, c));
    }
    // Mass of the solid cell is unchanged by the permutation.
    assert!((lbm.density_at(1, 1, 1) - before.density_at(1, 1, 1)).abs() < 1e-14);
}

/// A duplicated solid index must still reflect only once, and indices
/// beyond the grid are ignored.
#[test]
fn lbm_boundary_reflects_duplicate_solid_indices_once() {
    let mut lbm = LBMGpuKernels::new(2, 2, 1, 1.0, false, vec![1, 1, 99]);
    lbm.f.front.fill(0.0);
    let at = lbm.f.idx(1, 0, 0, 3);
    lbm.f.front[at] = 2.0;
    lbm.dispatch_boundary();
    assert_eq!(lbm.f.read(1, 0, 0, 4), 2.0, "f_3 must be reversed into f_4");
    assert_eq!(lbm.f.read(1, 0, 0, 3), 0.0);
}

/// Closed box with an obstacle: a full step is a permutation plus a
/// mass-conserving relaxation, so the total mass is constant.
#[test]
fn lbm_closed_box_with_obstacle_conserves_mass() {
    let dims = [6, 5, 4];
    let cell = |x: usize, y: usize, z: usize| x + dims[0] * (y + dims[1] * z);
    let solid = vec![cell(1, 3, 1), cell(2, 3, 1), cell(1, 4, 1)];
    let mut lbm = perturbed_lbm(dims, 0.75, false, solid);
    let (mass0, _) = lbm_totals(&lbm);
    for _ in 0..100 {
        lbm.dispatch_collision();
        lbm.dispatch_streaming();
        lbm.dispatch_boundary();
    }
    let (mass1, _) = lbm_totals(&lbm);
    assert!(
        (mass0 - mass1).abs() < 1e-9 * mass0,
        "mass drifted: {mass0} -> {mass1}"
    );
    assert!(lbm.f.front.iter().all(|v| v.is_finite()));
}

#[test]
fn lbm_rest_state_stays_at_rest_in_a_closed_box() {
    let mut lbm = LBMGpuKernels::new(4, 4, 4, 0.9, false, vec![21]);
    let initial = lbm.f.front.clone();
    for _ in 0..10 {
        lbm.dispatch_collision();
        lbm.dispatch_streaming();
        lbm.dispatch_boundary();
    }
    for (a, b) in lbm.f.front.iter().zip(initial.iter()) {
        assert!((a - b).abs() < 1e-14);
    }
}

/// Physical validation: a sinusoidal shear wave `u_y(x)` decays as
/// `exp(-nu k^2 t)` with `nu = cs^2 (tau - 1/2)`.
#[test]
fn lbm_shear_wave_decays_with_the_lattice_viscosity() {
    let (nx, tau, amplitude) = (32_usize, 0.8, 0.01);
    let k = 2.0 * std::f64::consts::PI / nx as f64;
    let mut lbm = LBMGpuKernels::new(nx, 1, 1, tau, true, vec![]);
    for x in 0..nx {
        let u = [0.0, amplitude * (k * x as f64).sin(), 0.0];
        for c in 0..D3Q19_Q {
            let at = lbm.f.idx(x, 0, 0, c);
            lbm.f.front[at] = feq_reference(c, 1.0, u);
        }
    }
    lbm.f.back.clone_from(&lbm.f.front);

    let measure = |lbm: &LBMGpuKernels| -> f64 {
        (0..nx)
            .map(|x| {
                let uy = lbm.momentum_at(x, 0, 0)[1] / lbm.density_at(x, 0, 0);
                uy * (k * x as f64).sin()
            })
            .sum::<f64>()
            * 2.0
            / nx as f64
    };
    assert!((measure(&lbm) - amplitude).abs() < 1e-12);

    let steps = 100;
    for _ in 0..steps {
        lbm.dispatch_collision();
        lbm.dispatch_streaming();
        lbm.dispatch_boundary();
    }
    let nu = (tau - 0.5) / 3.0;
    let expected = amplitude * (-nu * k * k * steps as f64).exp();
    let measured = measure(&lbm);
    assert!(
        (measured - expected).abs() < 0.02 * expected,
        "shear wave amplitude {measured} vs analytic {expected}"
    );
}

#[test]
fn lbm_density_at_single_cell() {
    let lbm = LBMGpuKernels::new(2, 2, 1, 1.0, false, vec![]);
    let rho = lbm.density_at(0, 0, 0);
    assert!(rho > 0.0);
}

#[test]
fn lbm_boundary_kind_eq() {
    assert_eq!(LBMBoundaryKind::BounceBack, LBMBoundaryKind::BounceBack);
    assert_ne!(LBMBoundaryKind::BounceBack, LBMBoundaryKind::MovingWall);
}

#[test]
fn lbm_full_step_keeps_rest_equilibrium() {
    let mut lbm = LBMGpuKernels::new(4, 4, 4, 1.0, true, vec![]);
    lbm.dispatch_collision();
    lbm.dispatch_streaming();
    lbm.dispatch_boundary();
    for (c, &weight) in D3Q19_WEIGHTS.iter().enumerate() {
        assert!((lbm.f.read(2, 1, 3, c) - weight).abs() < 1e-14);
    }
}

/// The boundary spec holds no wall velocity, inlet velocity or outlet density,
/// so every `LBMBoundaryKind` reflects the solid-cell populations exactly like
/// `BounceBack` (as the enum documents); fluid cells are never touched.
#[test]
fn lbm_boundary_kinds_all_reflect_like_bounce_back() {
    let solid = vec![13, 4];
    let untouched = perturbed_lbm([3, 3, 3], 1.0, false, solid.clone());
    let mut reference = untouched.clone();
    reference.dispatch_boundary();
    assert_ne!(
        reference.f.front, untouched.f.front,
        "the reflection must change the solid cells"
    );
    for kind in [
        LBMBoundaryKind::MovingWall,
        LBMBoundaryKind::ZouHeInlet,
        LBMBoundaryKind::ZouHeOutlet,
    ] {
        let mut lbm = untouched.clone();
        lbm.boundary.kind = kind.clone();
        lbm.dispatch_boundary();
        assert_eq!(lbm.f.front, reference.f.front, "{kind:?}");
    }
    // Only the two solid cells (indices 13 and 4) differ from the input.
    for cell in 0..27_usize {
        let range = cell * D3Q19_Q..(cell + 1) * D3Q19_Q;
        let changed = reference.f.front[range.clone()] != untouched.f.front[range];
        assert_eq!(changed, solid.contains(&cell), "cell {cell}");
    }
}
