//! Regression tests for the D2Q9 multiscale building blocks: Guo forcing in
//! `TwoWayCoupling`, acoustic-scaling relaxation times with Dupuis-Chopard
//! rescaling, and vertex-centred restriction in `PatchedGrid`.

use super::functions::{
    coarse_to_fine_pop, coarse_to_fine_scale, fine_relaxation_time, fine_to_coarse_scale,
    lbm_equilibrium,
};
use super::types::{
    AdaptiveMesh, GridRefinementInterface, PatchedGrid, SpaceTimeLbm, TwoWayCoupling,
};
use crate::lattice::{D2Q9_VELOCITIES, D2Q9_WEIGHTS};

const CS2: f64 = 1.0 / 3.0;

/// Zeroth and first velocity moments `(Σ f, [Σ c_x f, Σ c_y f])`.
fn moments(f: &[f64; 9]) -> (f64, [f64; 2]) {
    let rho: f64 = f.iter().sum();
    let mut j = [0.0_f64; 2];
    for (fi, c) in f.iter().zip(D2Q9_VELOCITIES.iter()) {
        j[0] += fi * c[0] as f64;
        j[1] += fi * c[1] as f64;
    }
    (rho, j)
}

/// Second velocity moment `Σ c_a c_b f` as `[[xx, xy], [yx, yy]]`.
fn second_moment(f: &[f64; 9]) -> [[f64; 2]; 2] {
    let mut pi = [[0.0_f64; 2]; 2];
    for (fi, c) in f.iter().zip(D2Q9_VELOCITIES.iter()) {
        for a in 0..2 {
            for b in 0..2 {
                pi[a][b] += fi * c[a] as f64 * c[b] as f64;
            }
        }
    }
    pi
}

/// Population perturbation with zero mass and zero momentum (a pure shear /
/// normal-stress mode): `w_i (c_x² - c_y²) + 0.5 w_i c_x c_y`.
fn stress_mode(amplitude: f64) -> [f64; 9] {
    let mut mode = [0.0_f64; 9];
    for (i, m) in mode.iter_mut().enumerate() {
        let c = D2Q9_VELOCITIES[i];
        let (cx, cy) = (c[0] as f64, c[1] as f64);
        *m = amplitude * D2Q9_WEIGHTS[i] * (cx * cx - cy * cy + 0.5 * cx * cy);
    }
    mode
}

/// A moving, non-equilibrium node state with density `rho` and velocity `u`.
fn nonequilibrium_state(rho: f64, u: [f64; 2], amplitude: f64) -> [f64; 9] {
    let mut f = lbm_equilibrium(rho, u);
    let mode = stress_mode(amplitude);
    for (fi, m) in f.iter_mut().zip(mode.iter()) {
        *fi += m;
    }
    f
}

/// Deterministic pseudo-random populations in `[0.01, 0.2)` (64-bit LCG).
fn pseudo_random_populations(seed: u64, count: usize) -> Vec<[f64; 9]> {
    let mut state = seed;
    let mut next = move || {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        0.01 + 0.19 * ((state >> 11) as f64 / (1u64 << 53) as f64)
    };
    (0..count)
        .map(|_| {
            let mut f = [0.0_f64; 9];
            for v in f.iter_mut() {
                *v = next();
            }
            f
        })
        .collect()
}

// ---------------------------------------------------------------------------
// TwoWayCoupling: Guo forcing (F2630)
// ---------------------------------------------------------------------------

#[test]
fn two_way_coupling_force_accelerates_a_node_at_rest_by_exactly_f() {
    for &tau in &[0.6, 0.8, 1.0, 1.7] {
        let mut coupling = TwoWayCoupling::new(1, 0, 1e-3, tau);
        coupling.body_force[0] = [1e-5, 0.0];
        coupling.collide_with_force();
        let (rho, j) = moments(&coupling.f[0]);
        assert!(
            (rho - 1.0).abs() < 1e-15,
            "tau={tau}: mass changed to {rho}"
        );
        assert!(
            (j[0] - 1e-5).abs() < 1e-15,
            "tau={tau}: Σ c_x f = {} instead of +1e-5 (the force must accelerate the fluid)",
            j[0]
        );
        assert!(
            j[1].abs() < 1e-15,
            "tau={tau}: spurious y momentum {}",
            j[1]
        );
    }
}

#[test]
fn two_way_coupling_force_injects_f_per_collision_for_a_moving_nonequilibrium_node() {
    let tau = 0.9;
    let force = [3.0e-4, -1.0e-4];
    let f_start = nonequilibrium_state(1.3, [0.03, -0.02], 0.004);
    let (rho0, j0) = moments(&f_start);
    let mut coupling = TwoWayCoupling::new(1, 0, 1e-3, tau);
    coupling.f[0] = f_start;
    coupling.body_force[0] = force;
    let collisions = 25;
    for step in 1..=collisions {
        coupling.collide_with_force();
        let (rho, j) = moments(&coupling.f[0]);
        assert!(
            (rho - rho0).abs() < 1e-14,
            "step {step}: forcing must not create mass ({rho} vs {rho0})"
        );
        for axis in 0..2 {
            let expected = j0[axis] + step as f64 * force[axis];
            assert!(
                (j[axis] - expected).abs() < 1e-13,
                "step {step} axis {axis}: Σ c f = {} expected {expected}",
                j[axis]
            );
        }
    }
}

#[test]
fn two_way_coupling_force_second_moment_contains_the_guo_stress_source() {
    // The Guo source has the second moment (1 - 1/(2τ)) (u_a F_b + u_b F_a):
    // the -u·F term cancels the isotropic part, so a missing term shows up as
    // an error (1 - 1/(2τ)) (u·F) δ_ab in the normal stresses.
    let tau = 0.75;
    let force = [2.0e-4, 1.0e-4];
    let f_start = nonequilibrium_state(1.2, [0.05, 0.02], 0.003);
    let (rho, j) = moments(&f_start);
    let u_half = [(j[0] + 0.5 * force[0]) / rho, (j[1] + 0.5 * force[1]) / rho];
    let pi_pre = second_moment(&f_start);
    let guo = 1.0 - 1.0 / (2.0 * tau);
    let mut coupling = TwoWayCoupling::new(1, 0, 1e-3, tau);
    coupling.f[0] = f_start;
    coupling.body_force[0] = force;
    coupling.collide_with_force();
    let pi_post = second_moment(&coupling.f[0]);
    for a in 0..2 {
        for b in 0..2 {
            let delta = if a == b { 1.0 } else { 0.0 };
            let pi_eq = rho * (CS2 * delta + u_half[a] * u_half[b]);
            let source = guo * (u_half[a] * force[b] + u_half[b] * force[a]);
            let expected = (1.0 - 1.0 / tau) * pi_pre[a][b] + pi_eq / tau + source;
            assert!(
                (pi_post[a][b] - expected).abs() < 1e-14,
                "Π[{a}][{b}] = {} expected {expected}",
                pi_post[a][b]
            );
        }
    }
}

#[test]
fn two_way_coupling_without_force_is_plain_bgk() {
    let tau = 1.3;
    let f_start = nonequilibrium_state(0.9, [0.02, 0.01], 0.002);
    let mut coupling = TwoWayCoupling::new(1, 0, 1e-3, tau);
    coupling.f[0] = f_start;
    coupling.collide_with_force();
    let (rho, j) = moments(&f_start);
    let feq = lbm_equilibrium(rho, [j[0] / rho, j[1] / rho]);
    for i in 0..9 {
        let expected = f_start[i] - (f_start[i] - feq[i]) / tau;
        assert!((coupling.f[0][i] - expected).abs() < 1e-15);
    }
}

// ---------------------------------------------------------------------------
// PatchedGrid: restriction by injection (F2636)
// ---------------------------------------------------------------------------

#[test]
fn patched_grid_restrict_keeps_the_density_of_a_uniform_field() {
    for &(nxc, nyc, ratio) in &[(3usize, 3usize, 2usize), (4, 3, 3), (5, 2, 4)] {
        let mut grid = PatchedGrid::new(nxc, nyc, ratio);
        let rest: [f64; 9] = D2Q9_WEIGHTS.map(|w| 1.7 * w);
        for f in grid.f_fine.iter_mut() {
            *f = rest;
        }
        for f in grid.f_coarse.iter_mut() {
            *f = [0.0; 9];
        }
        grid.restrict();
        for (k, f) in grid.f_coarse.iter().enumerate() {
            let rho: f64 = f.iter().sum();
            assert!(
                (rho - 1.7).abs() < 1e-14,
                "({nxc},{nyc},{ratio}) coarse node {k}: density {rho} instead of 1.7"
            );
            assert_eq!(*f, rest);
        }
    }
}

#[test]
fn patched_grid_restrict_injects_the_coincident_fine_node() {
    let mut grid = PatchedGrid::new(4, 3, 3);
    assert_eq!((grid.nx_fine, grid.ny_fine), (10, 7));
    for (kf, f) in grid.f_fine.iter_mut().enumerate() {
        for (i, v) in f.iter_mut().enumerate() {
            *v = (10 * kf + i) as f64;
        }
    }
    grid.restrict();
    for yc in 0..3 {
        for xc in 0..4 {
            let kf = (yc * 3) * grid.nx_fine + xc * 3;
            assert_eq!(grid.f_coarse[yc * 4 + xc], grid.f_fine[kf]);
        }
    }
}

#[test]
fn patched_grid_restrict_inverts_prolongate_exactly() {
    for &(nxc, nyc, ratio) in &[(3usize, 3usize, 2usize), (5, 4, 2), (4, 5, 3)] {
        let mut grid = PatchedGrid::new(nxc, nyc, ratio);
        let original = pseudo_random_populations(7 + ratio as u64, nxc * nyc);
        grid.f_coarse = original.clone();
        grid.prolongate();
        for f in grid.f_coarse.iter_mut() {
            *f = [-1.0; 9];
        }
        grid.restrict();
        assert_eq!(
            grid.f_coarse, original,
            "({nxc},{nyc},{ratio}): restrict(prolongate(f)) must reproduce f"
        );
    }
}

#[test]
fn patched_grid_prolongate_is_exact_for_bilinear_fields() {
    // Bilinear interpolation reproduces a field that is linear in x and y, so
    // the coincident nodes and the interior nodes must all lie on the plane.
    let (nxc, nyc, ratio) = (4usize, 3usize, 2usize);
    let mut grid = PatchedGrid::new(nxc, nyc, ratio);
    for yc in 0..nyc {
        for xc in 0..nxc {
            let v = 0.5 + 0.1 * xc as f64 + 0.05 * yc as f64;
            grid.f_coarse[yc * nxc + xc] = D2Q9_WEIGHTS.map(|w| w * v);
        }
    }
    grid.prolongate();
    for yf in 0..grid.ny_fine {
        for xf in 0..grid.nx_fine {
            let v = 0.5 + 0.1 * xf as f64 / ratio as f64 + 0.05 * yf as f64 / ratio as f64;
            let rho: f64 = grid.f_fine[yf * grid.nx_fine + xf].iter().sum();
            assert!(
                (rho - v).abs() < 1e-14,
                "fine node ({xf},{yf}): {rho} vs {v}"
            );
        }
    }
    grid.restrict();
    for yc in 0..nyc {
        for xc in 0..nxc {
            let rho: f64 = grid.f_coarse[yc * nxc + xc].iter().sum();
            let v = 0.5 + 0.1 * xc as f64 + 0.05 * yc as f64;
            assert!((rho - v).abs() < 1e-14);
        }
    }
}

// ---------------------------------------------------------------------------
// Acoustic scaling and Dupuis-Chopard rescaling (F2632)
// ---------------------------------------------------------------------------

#[test]
fn fine_relaxation_time_follows_acoustic_scaling() {
    assert!((fine_relaxation_time(1.0, 2) - 1.5).abs() < 1e-15);
    assert!((fine_relaxation_time(1.0, 4) - 2.5).abs() < 1e-15);
    assert!((fine_relaxation_time(0.6, 3) - 0.8).abs() < 1e-15);
    // A refinement ratio of one is the identity.
    assert!((fine_relaxation_time(0.83, 1) - 0.83).abs() < 1e-15);
}

#[test]
fn both_refinement_constructors_conserve_the_physical_viscosity() {
    for &tau_c in &[0.55, 0.8, 1.0, 1.7] {
        for &ratio in &[2usize, 3, 4] {
            let r = ratio as f64;
            let gri = GridRefinementInterface::new(4, ratio, tau_c, 1);
            let st = SpaceTimeLbm::new(4, ratio, tau_c);
            // Acoustic scaling: dx_f = dx_c / r and dt_f = dt_c / r, so
            // ν = cs² (τ - 1/2) dx² / dt is  cs² (τ_f - 1/2) / r  on the fine
            // grid (in coarse units) and  cs² (τ_c - 1/2)  on the coarse grid.
            let nu_coarse = CS2 * (tau_c - 0.5);
            for tau_fine in [gri.tau_fine, st.tau_fine] {
                let nu_fine = CS2 * (tau_fine - 0.5) / r;
                assert!(
                    ((nu_fine - nu_coarse) / nu_coarse).abs() < 1e-13,
                    "tau_c={tau_c}, r={ratio}: nu_fine={nu_fine}, nu_coarse={nu_coarse}"
                );
                assert!((tau_fine - fine_relaxation_time(tau_c, ratio)).abs() < 1e-15);
            }
        }
    }
}

#[test]
fn dupuis_chopard_scales_are_reciprocal_and_match_the_literature_formula() {
    let (tau_c, ratio) = (1.0, 2usize);
    let tau_f = fine_relaxation_time(tau_c, ratio);
    assert!((coarse_to_fine_scale(tau_c, tau_f, ratio) - 0.75).abs() < 1e-15);
    assert!((fine_to_coarse_scale(tau_c, tau_f, ratio) - 4.0 / 3.0).abs() < 1e-15);
    for &tau_c in &[0.55, 0.8, 1.7] {
        for &ratio in &[2usize, 3, 4] {
            let tau_f = fine_relaxation_time(tau_c, ratio);
            let product = coarse_to_fine_scale(tau_c, tau_f, ratio)
                * fine_to_coarse_scale(tau_c, tau_f, ratio);
            assert!((product - 1.0).abs() < 1e-14);
        }
    }
}

#[test]
fn grid_refinement_rescale_scales_the_nonequilibrium_part_by_the_dupuis_chopard_factor() {
    let (rho, u) = (1.1, [0.04, -0.03]);
    let gri = GridRefinementInterface::new(4, 2, 0.8, 1);
    let f_coarse = nonequilibrium_state(rho, u, 0.006);
    let feq = lbm_equilibrium(rho, u);
    let f_fine = gri.coarse_to_fine_rescale(&f_coarse, rho, u);
    let scale = gri.tau_fine / (2.0 * gri.tau_coarse);
    for i in 0..9 {
        let neq_coarse = f_coarse[i] - feq[i];
        let neq_fine = f_fine[i] - feq[i];
        assert!(
            (neq_fine - scale * neq_coarse).abs() < 1e-15,
            "direction {i}: f_neq_fine = {neq_fine}, expected {}",
            scale * neq_coarse
        );
    }
    // Mass and momentum are untouched by the non-equilibrium rescaling.
    let (rho_f, j_f) = moments(&f_fine);
    assert!((rho_f - rho).abs() < 1e-14);
    assert!((j_f[0] - rho * u[0]).abs() < 1e-14 && (j_f[1] - rho * u[1]).abs() < 1e-14);
    // The reverse transfer undoes the forward transfer exactly.
    let f_back = gri.fine_to_coarse_rescale(&f_fine, rho, u);
    for i in 0..9 {
        assert!((f_back[i] - f_coarse[i]).abs() < 1e-15);
    }
}

#[test]
fn rescaling_keeps_the_physical_viscous_stress_continuous() {
    // The viscous stress felt by the fluid is (1 - 1/(2τ)) Σ c_a c_b f_neq.  With
    // acoustic scaling the Dupuis-Chopard factor makes it identical on both
    // levels; the inverted formulas of the previous implementation do not.
    let (rho, u) = (1.1, [0.04, -0.03]);
    let feq = lbm_equilibrium(rho, u);
    let neq_stress = |f: &[f64; 9]| {
        let mut neq = *f;
        for (n, e) in neq.iter_mut().zip(feq.iter()) {
            *n -= e;
        }
        second_moment(&neq)
    };
    for &(tau_c, ratio) in &[(0.6, 2usize), (0.8, 2), (1.0, 3), (1.7, 4)] {
        let tau_f = fine_relaxation_time(tau_c, ratio);
        let f_coarse = nonequilibrium_state(rho, u, 0.006);
        let f_fine = coarse_to_fine_pop(&f_coarse, rho, u, u, tau_c, tau_f, ratio);
        let pi_coarse = neq_stress(&f_coarse);
        let pi_fine = neq_stress(&f_fine);
        let mut checked = 0;
        for a in 0..2 {
            for b in 0..2 {
                let stress_coarse = (1.0 - 0.5 / tau_c) * pi_coarse[a][b];
                let stress_fine = (1.0 - 0.5 / tau_f) * pi_fine[a][b];
                assert!(
                    (stress_coarse - stress_fine).abs() < 1e-15,
                    "tau_c={tau_c}, r={ratio}, component ({a},{b}): {stress_coarse} vs {stress_fine}"
                );
                if stress_coarse.abs() > 1e-5 {
                    checked += 1;
                }
            }
        }
        assert!(checked >= 3, "the stress components must be non-trivial");
        // The previous formula (tau_f = (tau_c - 1/2)/r + 1/2 with the factor
        // tau_c/tau_f) breaks this continuity.
        let tau_old = (tau_c - 0.5) / ratio as f64 + 0.5;
        let old_scale = tau_c / tau_old;
        let stress_old = (1.0 - 0.5 / tau_old) * old_scale * pi_coarse[0][1];
        let stress_coarse = (1.0 - 0.5 / tau_c) * pi_coarse[0][1];
        assert!(((stress_old - stress_coarse) / stress_coarse).abs() > 0.05);
    }
}

#[test]
fn equilibrium_maps_to_itself_across_the_refinement_interface_for_a_moving_fluid() {
    let (rho, u) = (1.3, [0.07, -0.04]);
    let feq = lbm_equilibrium(rho, u);
    let gri = GridRefinementInterface::new(4, 2, 0.8, 1);
    let down = gri.coarse_to_fine_rescale(&feq, rho, u);
    let up = gri.fine_to_coarse_rescale(&feq, rho, u);
    let tau_f = fine_relaxation_time(0.8, 2);
    let pop = coarse_to_fine_pop(&feq, rho, u, u, 0.8, tau_f, 2);
    for i in 0..9 {
        assert!(
            (down[i] - feq[i]).abs() < 1e-15,
            "coarse_to_fine_rescale {i}"
        );
        assert!((up[i] - feq[i]).abs() < 1e-15, "fine_to_coarse_rescale {i}");
        assert!((pop[i] - feq[i]).abs() < 1e-15, "coarse_to_fine_pop {i}");
    }
}

#[test]
fn coarse_to_fine_pop_uses_the_fine_velocity_and_never_clamps() {
    let rho = 1.0;
    let (u_coarse, u_fine) = ([0.3, 0.0], [-0.3, 0.0]);
    let (tau_c, ratio) = (0.8, 2usize);
    let tau_f = fine_relaxation_time(tau_c, ratio);
    let scale = coarse_to_fine_scale(tau_c, tau_f, ratio);
    // A fluid at rest is a strongly non-equilibrium state for u_coarse.
    let f_coarse: [f64; 9] = D2Q9_WEIGHTS.map(|w| w * rho);
    let feq_coarse = lbm_equilibrium(rho, u_coarse);
    let feq_fine = lbm_equilibrium(rho, u_fine);
    let result = coarse_to_fine_pop(&f_coarse, rho, u_coarse, u_fine, tau_c, tau_f, ratio);
    let mut negative_entries = 0;
    for i in 0..9 {
        let expected = feq_fine[i] + scale * (f_coarse[i] - feq_coarse[i]);
        assert!((result[i] - expected).abs() < 1e-15, "direction {i}");
        if result[i] < 0.0 {
            negative_entries += 1;
        }
    }
    assert!(
        negative_entries > 0,
        "the case must exercise negative results (no clamping to zero)"
    );
}

#[test]
fn quadtree_children_inherit_the_parent_velocity_and_nonequilibrium_part() {
    let (rho, u) = (1.2, [0.05, -0.03]);
    let tau_parent = 0.8;
    let mut mesh = AdaptiveMesh::new(1.0, 2, 0.01);
    mesh.nodes[0].rho = rho;
    mesh.nodes[0].u = u;
    mesh.nodes[0].f = lbm_equilibrium(rho, u);
    let eq_children = mesh.refine(0, tau_parent).expect("root is refinable");
    for &child in &eq_children {
        let node = &mesh.nodes[child];
        assert_eq!(node.u, u);
        assert_eq!(node.rho, rho);
        let feq = lbm_equilibrium(rho, u);
        for (value, expected) in node.f.iter().zip(feq.iter()) {
            assert!(
                (value - expected).abs() < 1e-15,
                "child {child}: an equilibrium parent must give equilibrium children"
            );
        }
    }

    let mut mesh = AdaptiveMesh::new(1.0, 2, 0.01);
    let parent_f = nonequilibrium_state(rho, u, 0.005);
    mesh.nodes[0].rho = rho;
    mesh.nodes[0].u = u;
    mesh.nodes[0].f = parent_f;
    let children = mesh.refine(0, tau_parent).expect("root is refinable");
    let tau_child = 2.0 * (tau_parent - 0.5) + 0.5;
    let scale = tau_child / (2.0 * tau_parent);
    let feq = lbm_equilibrium(rho, u);
    for &child in &children {
        let f_child = mesh.nodes[child].f;
        let (rho_c, j_c) = moments(&f_child);
        assert!((rho_c - rho).abs() < 1e-14, "child mass");
        assert!((j_c[0] - rho * u[0]).abs() < 1e-14 && (j_c[1] - rho * u[1]).abs() < 1e-14);
        for i in 0..9 {
            let expected = feq[i] + scale * (parent_f[i] - feq[i]);
            assert!((f_child[i] - expected).abs() < 1e-15);
        }
    }
}

// ---------------------------------------------------------------------------
// End-to-end check of the acoustic scaling with a decaying shear wave (F2632)
// ---------------------------------------------------------------------------

/// One periodic D2Q9 BGK step (collision followed by push streaming).
fn bgk_step(f: &mut [[f64; 9]], nx: usize, ny: usize, tau: f64) {
    let mut post = f.to_vec();
    for (cell, out) in f.iter().zip(post.iter_mut()) {
        let (rho, j) = moments(cell);
        let feq = lbm_equilibrium(rho, [j[0] / rho, j[1] / rho]);
        for i in 0..9 {
            out[i] = cell[i] - (cell[i] - feq[i]) / tau;
        }
    }
    for y in 0..ny {
        for x in 0..nx {
            for (i, c) in D2Q9_VELOCITIES.iter().enumerate() {
                let xn = (x as i64 + c[0] as i64).rem_euclid(nx as i64) as usize;
                let yn = (y as i64 + c[1] as i64).rem_euclid(ny as i64) as usize;
                f[yn * nx + xn][i] = post[y * nx + x][i];
            }
        }
    }
}

/// Relative amplitude of the `sin(2π (y + 1/2) / n)` mode of `u_x` after
/// `steps` BGK steps on an `n x n` periodic grid, starting from an
/// equilibrium shear wave of amplitude `u0`.
fn shear_wave_amplitude_ratio(n: usize, tau: f64, steps: usize) -> f64 {
    let u0 = 0.02;
    let phase = |y: usize| std::f64::consts::TAU * (y as f64 + 0.5) / n as f64;
    let mut f: Vec<[f64; 9]> = (0..n * n)
        .map(|k| lbm_equilibrium(1.0, [u0 * phase(k / n).sin(), 0.0]))
        .collect();
    for _ in 0..steps {
        bgk_step(&mut f, n, n, tau);
    }
    let mut projection = 0.0;
    for (k, cell) in f.iter().enumerate() {
        let (rho, j) = moments(cell);
        projection += j[0] / rho * phase(k / n).sin();
    }
    2.0 * projection / (n * n) as f64 / u0
}

#[test]
fn shear_wave_decays_at_the_same_physical_rate_on_both_refinement_levels() {
    let (n_coarse, ratio, coarse_steps) = (16usize, 2usize, 40usize);
    let tau_c = 0.8;
    let nu = CS2 * (tau_c - 0.5);
    let k = std::f64::consts::TAU / n_coarse as f64;
    let analytic = (-nu * k * k * coarse_steps as f64).exp();

    let coarse = shear_wave_amplitude_ratio(n_coarse, tau_c, coarse_steps);
    let tau_f = fine_relaxation_time(tau_c, ratio);
    let fine = shear_wave_amplitude_ratio(n_coarse * ratio, tau_f, coarse_steps * ratio);
    assert!(
        ((coarse - analytic) / analytic).abs() < 0.025,
        "coarse {coarse} vs analytic {analytic}"
    );
    assert!(
        ((fine - analytic) / analytic).abs() < 0.01,
        "fine {fine} vs analytic {analytic}"
    );
    assert!(
        ((fine - coarse) / coarse).abs() < 0.03,
        "fine {fine} vs coarse {coarse}"
    );

    // Negative control: the previous (inverted) formula τ_f = (τ_c - 1/2)/r + 1/2
    // makes the fine grid four times too inviscid in coarse units.
    let tau_inverted = (tau_c - 0.5) / ratio as f64 + 0.5;
    let wrong = shear_wave_amplitude_ratio(n_coarse * ratio, tau_inverted, coarse_steps * ratio);
    assert!(
        ((wrong - analytic) / analytic).abs() > 0.3,
        "the inverted scaling must be visibly wrong, got {wrong} vs {analytic}"
    );
}
