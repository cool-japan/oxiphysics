//! Regression tests for the D1Q3 lattice solvers of the multiscale module:
//! [`AdaptiveLbm`], [`MesoscaleLbm`] and [`HybridLbmNavier`].

use super::d1q3::{CS2, first_moment, hydrodynamic_equilibrium, zeroth_moment};
use super::{AdaptiveLbm, HybridLbmNavier, MesoscaleLbm};
use crate::error::Error;

/// Centroid and variance of a non-negative field sampled at `x_i = offset + i * dx`.
fn centroid_and_variance(field: &[f64], offset: f64, dx: f64) -> (f64, f64) {
    let mass: f64 = field.iter().sum();
    let position = |i: usize| offset + i as f64 * dx;
    let mean = field
        .iter()
        .enumerate()
        .map(|(i, &v)| v * position(i))
        .sum::<f64>()
        / mass;
    let variance = field
        .iter()
        .enumerate()
        .map(|(i, &v)| v * (position(i) - mean).powi(2))
        .sum::<f64>()
        / mass;
    (mean, variance)
}

/// Exact growth of the variance (in squared cells) after `steps` BGK steps of
/// the D1Q3 diffusion scheme (`u = 0`) started from equilibrium populations.
///
/// The second moment, the position-weighted flux and the second velocity moment
/// obey a closed linear recursion, which sums to
/// `cs² t + 2 cs² (1 - ω)/ω (t - (1 - (1 - ω)^t)/ω)`.
fn exact_variance_growth(omega: f64, steps: usize) -> f64 {
    let t = steps as f64;
    CS2 * t
        + 2.0 * CS2 * (1.0 - omega) / omega * (t - (1.0 - (1.0 - omega).powi(steps as i32)) / omega)
}

/// Gaussian profile `exp(-(i - centre)² / (2 σ²))` on `n` cells.
fn gaussian(n: usize, centre: f64, sigma: f64) -> Vec<f64> {
    (0..n)
        .map(|i| (-(i as f64 - centre).powi(2) / (2.0 * sigma * sigma)).exp())
        .collect()
}

/// Deterministic pseudo-random numbers in `[0, 1)` (64-bit LCG).
fn pseudo_random(seed: u64, count: usize) -> Vec<f64> {
    let mut state = seed;
    (0..count)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 11) as f64 / (1u64 << 53) as f64
        })
        .collect()
}

// ---------------------------------------------------------------------------
// AdaptiveLbm (F2639)
// ---------------------------------------------------------------------------

fn adaptive(n: usize, levels: &[usize], omega: f64) -> AdaptiveLbm {
    AdaptiveLbm::new(n, levels.to_vec(), omega).expect("valid adaptive LBM parameters")
}

/// Position of fine cell `i` of a level with ratio `r` in coarse-cell units.
fn fine_offset(r: usize) -> (f64, f64) {
    (0.5 / r as f64 - 0.5, 1.0 / r as f64)
}

#[test]
fn adaptive_new_validates_its_parameters() {
    assert!(AdaptiveLbm::new(16, vec![2, 4], 1.0).is_ok());
    assert!(AdaptiveLbm::new(16, vec![], 1.99).is_ok());
    assert_eq!(
        AdaptiveLbm::new(0, vec![2], 1.0).err(),
        Some(Error::InvalidGridDimensions { nx: 0, ny: 1 })
    );
    assert!(matches!(
        AdaptiveLbm::new(8, vec![2, 0], 1.0),
        Err(Error::IncompatibleParameters { .. })
    ));
    assert!(matches!(
        AdaptiveLbm::new(8, vec![2], 2.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    assert!(matches!(
        AdaptiveLbm::new(8, vec![2], 2.5),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    for bad in [0.0, -0.5, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            AdaptiveLbm::new(8, vec![2], bad),
            Err(Error::ParameterOutOfBounds { .. })
        ));
    }
}

#[test]
fn adaptive_levels_have_the_requested_sizes_and_start_empty() {
    let a = adaptive(16, &[2, 4], 1.0);
    assert_eq!(a.coarse_grid.len(), 16);
    assert_eq!(a.fine_grids[0].len(), 32);
    assert_eq!(a.fine_grids[1].len(), 64);
    assert!(a.coarse_grid.iter().all(|&v| v == 0.0));
    assert!(a.fine_grids.iter().flatten().all(|&v| v == 0.0));
}

#[test]
fn adaptive_relaxation_times_follow_acoustic_scaling() {
    let a = adaptive(16, &[2, 4], 1.25);
    assert!((a.relaxation_time() - 0.8).abs() < 1e-15);
    assert!((a.fine_relaxation_time(0).unwrap() - 1.1).abs() < 1e-14);
    assert!((a.fine_relaxation_time(1).unwrap() - 1.7).abs() < 1e-14);
    assert!(a.fine_relaxation_time(2).is_none());
    // The diffusivity in coarse units cs² (τ_f - 1/2) / r is the same on every level.
    let d = a.diffusivity();
    assert!((d - 0.1).abs() < 1e-15);
    for (level, &r) in a.refinement_levels.iter().enumerate() {
        let tau_f = a.fine_relaxation_time(level).unwrap();
        assert!((CS2 * (tau_f - 0.5) / r as f64 - d).abs() < 1e-15);
    }
}

#[test]
fn adaptive_spike_splits_two_thirds_one_sixth_whatever_omega() {
    // From equilibrium a D1Q3 BGK step returns the equilibrium, so one step
    // moves exactly 1/6 of the scalar to each neighbour for every ω.  A
    // relaxation towards zero (the previous implementation) cannot do this.
    for &omega in &[0.5, 1.0, 1.5, 1.9] {
        let mut a = adaptive(16, &[], omega);
        a.coarse_grid[5] = 1.0;
        a.step();
        for (i, &phi) in a.coarse_grid.iter().enumerate() {
            let expected = match i {
                4 | 6 => 1.0 / 6.0,
                5 => 2.0 / 3.0,
                _ => 0.0,
            };
            assert!((phi - expected).abs() < 1e-15, "ω={omega} cell {i}: {phi}");
        }
    }
}

#[test]
fn adaptive_spike_wraps_around_the_periodic_domain() {
    let mut a = adaptive(8, &[], 1.2);
    a.coarse_grid[0] = 3.0;
    a.step();
    assert!((a.coarse_grid[7] - 0.5).abs() < 1e-15);
    assert!((a.coarse_grid[0] - 2.0).abs() < 1e-15);
    assert!((a.coarse_grid[1] - 0.5).abs() < 1e-15);
}

#[test]
fn adaptive_step_advances_every_fine_level_by_exactly_r_substeps() {
    let n = 16;
    let spike = 8;
    let mut a = adaptive(n, &[2, 4], 1.25);
    a.coarse_grid[spike] = 1.0;
    a.refine_region(0, 0, n);
    a.refine_region(1, 0, n);
    for (level, &r) in a.refinement_levels.clone().iter().enumerate() {
        // The spike is prolonged onto the r fine cells of its coarse cell.
        let lo = spike * r;
        let hi = lo + r - 1;
        for (i, &phi) in a.fine_grids[level].iter().enumerate() {
            assert_eq!(phi != 0.0, (lo..=hi).contains(&i), "level {level} cell {i}");
        }
    }
    a.step();
    // The coarse level moved by one cell, a level with ratio r by r fine cells.
    let support = |field: &[f64]| -> (usize, usize) {
        let first = field.iter().position(|&v| v != 0.0).expect("non-empty");
        let last = field.iter().rposition(|&v| v != 0.0).expect("non-empty");
        (first, last)
    };
    assert_eq!(support(&a.coarse_grid), (spike - 1, spike + 1));
    for (level, &r) in a.refinement_levels.iter().enumerate() {
        assert_eq!(
            support(&a.fine_grids[level]),
            (spike * r - r, spike * r + 2 * r - 1),
            "level {level} (r = {r})"
        );
    }
}

#[test]
fn adaptive_conserves_scalar_mass_on_every_level() {
    let n = 48;
    let mut a = adaptive(n, &[2, 3], 1.4);
    a.set_advection_velocity(0.2)
        .expect("within the stable range");
    for (cell, value) in a.coarse_grid.iter_mut().zip(pseudo_random(5, n)) {
        *cell = 0.5 + value;
    }
    a.refine_region(0, 0, n);
    a.refine_region(1, 0, n);
    let coarse_mass = a.total_mass();
    for level in 0..2 {
        let mass = a.level_mass(level).expect("known level");
        assert!(
            (mass - coarse_mass).abs() < 1e-12 * coarse_mass,
            "prolongation must conserve mass: level {level} has {mass}, coarse {coarse_mass}"
        );
    }
    for _ in 0..200 {
        a.step();
    }
    assert!((a.total_mass() - coarse_mass).abs() < 1e-12 * coarse_mass);
    for level in 0..2 {
        let mass = a.level_mass(level).expect("known level");
        assert!(
            (mass - coarse_mass).abs() < 1e-12 * coarse_mass,
            "level {level}: {mass} vs {coarse_mass}"
        );
    }
    // Restriction conserves mass as well.
    a.coarsen_region(1, 0, n);
    assert!((a.total_mass() - coarse_mass).abs() < 1e-12 * coarse_mass);
}

#[test]
fn adaptive_gaussian_variance_grows_at_the_exact_lattice_diffusion_rate_on_every_level() {
    let n = 256;
    let omega = 1.25;
    let mut a = adaptive(n, &[2, 4], omega);
    a.coarse_grid.copy_from_slice(&gaussian(n, 128.0, 6.0));
    a.refine_region(0, 0, n);
    a.refine_region(1, 0, n);
    let steps = 120;

    let variance_of = |a: &AdaptiveLbm, level: Option<usize>| -> f64 {
        match level {
            None => centroid_and_variance(&a.coarse_grid, 0.0, 1.0).1,
            Some(l) => {
                let (offset, dx) = fine_offset(a.refinement_levels[l]);
                centroid_and_variance(&a.fine_grids[l], offset, dx).1
            }
        }
    };
    let initial = [
        variance_of(&a, None),
        variance_of(&a, Some(0)),
        variance_of(&a, Some(1)),
    ];
    for _ in 0..steps {
        a.step();
    }
    let growth = [
        variance_of(&a, None) - initial[0],
        variance_of(&a, Some(0)) - initial[1],
        variance_of(&a, Some(1)) - initial[2],
    ];

    // Coarse level: exact result of the scheme.
    let exact_coarse = exact_variance_growth(omega, steps);
    assert!(
        (growth[0] - exact_coarse).abs() < 1e-9,
        "coarse growth {} vs exact {exact_coarse}",
        growth[0]
    );
    // Refinement levels: r steps of size 1/r per coarse step with ω_f = 1/τ_f,
    // measured in squared coarse cells.
    for (level, &r) in a.refinement_levels.iter().enumerate() {
        let omega_f = 1.0 / a.fine_relaxation_time(level).expect("known level");
        let exact = exact_variance_growth(omega_f, steps * r) / (r * r) as f64;
        assert!(
            (growth[level + 1] - exact).abs() < 1e-9,
            "level {level}: growth {} vs exact {exact}",
            growth[level + 1]
        );
    }
    // All levels diffuse with the same physical diffusivity D: Var grows as 2 D t.
    let physical = 2.0 * a.diffusivity() * steps as f64;
    for (level, g) in growth.iter().enumerate() {
        assert!(
            ((g - physical) / physical).abs() < 0.02,
            "level {level}: variance growth {g} vs 2 D t = {physical}"
        );
    }
}

#[test]
fn adaptive_advection_translates_every_level_with_the_same_velocity() {
    let n = 256;
    let u = 0.1;
    let mut a = adaptive(n, &[2, 4], 1.5);
    a.set_advection_velocity(u)
        .expect("within the stable range");
    a.coarse_grid.copy_from_slice(&gaussian(n, 100.0, 5.0));
    a.refine_region(0, 0, n);
    a.refine_region(1, 0, n);
    let centroids = |a: &AdaptiveLbm| -> [f64; 3] {
        let coarse = centroid_and_variance(&a.coarse_grid, 0.0, 1.0).0;
        let mut out = [coarse, 0.0, 0.0];
        for (level, &r) in a.refinement_levels.iter().enumerate() {
            let (offset, dx) = fine_offset(r);
            out[level + 1] = centroid_and_variance(&a.fine_grids[level], offset, dx).0;
        }
        out
    };
    let start = centroids(&a);
    let steps = 150;
    for _ in 0..steps {
        a.step();
    }
    let end = centroids(&a);
    for level in 0..3 {
        let shift = end[level] - start[level];
        assert!(
            (shift - u * steps as f64).abs() < 1e-9,
            "level {level}: centroid moved by {shift}, expected {}",
            u * steps as f64
        );
    }
}

#[test]
fn adaptive_refine_then_coarsen_restores_the_coarse_populations() {
    let n = 64;
    let mut a = adaptive(n, &[2, 3], 1.25);
    a.set_advection_velocity(0.1)
        .expect("within the stable range");
    a.coarse_grid.copy_from_slice(&gaussian(n, 30.0, 4.0));
    // Develop a genuinely non-equilibrium state first.
    for _ in 0..25 {
        a.step();
    }
    let before = a.coarse_pops.clone();
    for level in 0..2 {
        let mut b = a.clone();
        b.refine_region(level, 0, n);
        b.coarsen_region(level, 0, n);
        for (cell, (restored, original)) in b.coarse_pops.iter().zip(before.iter()).enumerate() {
            for (x, y) in restored.iter().zip(original.iter()) {
                assert!(
                    (x - y).abs() < 1e-14,
                    "level {level} cell {cell}: {restored:?} vs {original:?}"
                );
            }
        }
    }
}

#[test]
fn adaptive_prolongation_scales_the_nonequilibrium_part_by_the_dupuis_chopard_factor() {
    let n = 64;
    let (omega, u) = (1.25, 0.1);
    let mut a = adaptive(n, &[2], omega);
    a.set_advection_velocity(u)
        .expect("within the stable range");
    a.coarse_grid.copy_from_slice(&gaussian(n, 30.0, 4.0));
    for _ in 0..25 {
        a.step();
    }
    let tau_c = a.relaxation_time();
    let tau_f = a.fine_relaxation_time(0).expect("known level");
    let scale = tau_f / (2.0 * tau_c);
    a.refine_region(0, 0, n);
    let equilibrium = |phi: f64| {
        [
            2.0 / 3.0 * phi,
            (1.0 + 3.0 * u) * phi / 6.0,
            (1.0 - 3.0 * u) * phi / 6.0,
        ]
    };
    let mut checked_non_trivial = 0;
    for ci in 0..n {
        let g = a.coarse_pops[ci];
        let neq_coarse: Vec<f64> = g
            .iter()
            .zip(equilibrium(zeroth_moment(&g)).iter())
            .map(|(x, e)| x - e)
            .collect();
        for k in 0..2 {
            let node = a.fine_pops[0][2 * ci + k];
            for (q, (x, e)) in node
                .iter()
                .zip(equilibrium(zeroth_moment(&node)).iter())
                .enumerate()
            {
                assert!(
                    (x - e - scale * neq_coarse[q]).abs() < 1e-14,
                    "cell {ci} sub-cell {k} direction {q}"
                );
            }
        }
        if neq_coarse.iter().any(|v| v.abs() > 1e-4) {
            checked_non_trivial += 1;
        }
    }
    assert!(
        checked_non_trivial > 5,
        "the test needs a non-equilibrium state"
    );
}

#[test]
fn adaptive_prolongation_keeps_the_physical_diffusive_flux_continuous() {
    // The diffusive flux seen by the scalar is (1 - 1/(2τ)) Σ c_i g_i^neq; the
    // Dupuis-Chopard factor together with τ_f = r (τ_c - 1/2) + 1/2 makes it the
    // same on the coarse and on every refined level.
    let n = 64;
    let (omega, u) = (1.25, 0.1);
    let mut a = adaptive(n, &[2, 3], omega);
    a.set_advection_velocity(u)
        .expect("within the stable range");
    a.coarse_grid.copy_from_slice(&gaussian(n, 30.0, 4.0));
    for _ in 0..25 {
        a.step();
    }
    let equilibrium = |phi: f64| {
        [
            2.0 / 3.0 * phi,
            (1.0 + 3.0 * u) * phi / 6.0,
            (1.0 - 3.0 * u) * phi / 6.0,
        ]
    };
    let flux_of = |g: &[f64; 3], tau: f64| {
        let eq = equilibrium(zeroth_moment(g));
        (1.0 - 0.5 / tau) * ((g[1] - eq[1]) - (g[2] - eq[2]))
    };
    let tau_c = a.relaxation_time();
    let coarse_pops = a.coarse_pops.clone();
    a.refine_region(0, 0, n);
    a.refine_region(1, 0, n);
    let mut largest = 0.0_f64;
    for (level, &r) in a.refinement_levels.iter().enumerate() {
        let tau_f = a.fine_relaxation_time(level).expect("known level");
        for (ci, g) in coarse_pops.iter().enumerate() {
            let flux_coarse = flux_of(g, tau_c);
            largest = largest.max(flux_coarse.abs());
            for k in 0..r {
                let flux_fine = flux_of(&a.fine_pops[level][r * ci + k], tau_f);
                assert!(
                    (flux_fine - flux_coarse).abs() < 1e-14,
                    "level {level} cell {ci}.{k}: {flux_fine} vs {flux_coarse}"
                );
            }
        }
    }
    assert!(largest > 1e-3, "the test needs a genuinely diffusing state");
}

#[test]
fn adaptive_prolongation_is_conservative_and_creates_no_new_extrema() {
    let mut a = adaptive(8, &[2, 4], 1.0);
    // Linear ramp: the limited slope is the exact slope in the interior.
    let ramp = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0];
    a.coarse_grid.copy_from_slice(&ramp);
    a.refine_region(0, 1, 7);
    for (ci, &value) in ramp.iter().enumerate().take(7).skip(1) {
        let fine = &a.fine_grids[0][2 * ci..2 * ci + 2];
        assert!((fine[0] - (value - 0.25)).abs() < 1e-14);
        assert!((fine[1] - (value + 0.25)).abs() < 1e-14);
    }
    // Local extremum and plateau: piecewise-constant, no overshoot.
    a.coarse_grid
        .copy_from_slice(&[0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 1.0, 1.0]);
    a.refine_region(1, 0, 8);
    let fine = &a.fine_grids[1];
    assert!(fine.iter().all(|&v| (0.0..=5.0).contains(&v)));
    for (ci, &expected) in [0.0, 0.0, 5.0, 0.0, 0.0, 1.0, 1.0, 1.0].iter().enumerate() {
        let cell_mean: f64 = fine[4 * ci..4 * ci + 4].iter().sum::<f64>() / 4.0;
        assert!((cell_mean - expected).abs() < 1e-14);
    }
}

#[test]
fn adaptive_levels_stay_consistent_after_many_steps() {
    // A smooth profile evolved on the coarse level and on a refined replica
    // agrees once the replica is restricted back to the coarse cells.
    let n = 128;
    let mut a = adaptive(n, &[2], 1.25);
    a.set_advection_velocity(0.05)
        .expect("within the stable range");
    a.coarse_grid.copy_from_slice(&gaussian(n, 64.0, 5.0));
    for _ in 0..40 {
        a.step();
    }
    a.refine_region(0, 0, n);
    for _ in 0..100 {
        a.step();
    }
    let peak = a.coarse_grid.iter().cloned().fold(0.0_f64, f64::max);
    let worst = a
        .coarse_grid
        .iter()
        .enumerate()
        .map(|(ci, &phi)| {
            let restricted = 0.5 * (a.fine_grids[0][2 * ci] + a.fine_grids[0][2 * ci + 1]);
            (phi - restricted).abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        worst < 3e-3 * peak,
        "coarse and restricted fine solutions differ by {worst} (peak {peak})"
    );
}

#[test]
fn adaptive_external_field_writes_are_imposed_on_the_populations() {
    let mut a = adaptive(4, &[2], 1.0);
    a.coarse_grid[1] = 5.0;
    a.refine_region(0, 1, 2);
    assert!((a.fine_grids[0][2] - 5.0).abs() < 1e-12);
    assert!((a.fine_grids[0][3] - 5.0).abs() < 1e-12);
    // Cells that were never refined stay empty.
    assert_eq!(a.fine_grids[0][0], 0.0);
    assert_eq!(a.fine_grids[0][4], 0.0);

    let mut b = adaptive(4, &[2], 1.0);
    b.fine_grids[0][0] = 2.0;
    b.fine_grids[0][1] = 4.0;
    b.coarsen_region(0, 0, 1);
    assert!((b.coarse_grid[0] - 3.0).abs() < 1e-12);
    assert_eq!(b.coarse_grid[1], 0.0);

    // A written value is not forgotten by the next step: it diffuses.
    let mut c = adaptive(8, &[], 1.0);
    c.coarse_grid[3] = 1.0;
    c.step();
    assert!((c.coarse_grid[3] - 2.0 / 3.0).abs() < 1e-15);
}

#[test]
fn adaptive_unknown_level_is_ignored_and_advection_velocity_is_validated() {
    let mut a = adaptive(8, &[2], 1.0);
    a.coarse_grid[2] = 1.0;
    a.refine_region(5, 0, 8);
    a.coarsen_region(5, 0, 8);
    assert!((a.coarse_grid[2] - 1.0).abs() < 1e-15);
    assert!(a.level_mass(5).is_none());
    assert!(a.set_advection_velocity(0.3).is_ok());
    assert!(a.set_advection_velocity(-1.0 / 3.0).is_ok());
    assert!(matches!(
        a.set_advection_velocity(0.4),
        Err(Error::ParameterOutOfBounds { .. })
    ));
    assert!((a.advection_velocity + 1.0 / 3.0).abs() < 1e-15);
}

#[test]
fn adaptive_time_and_spatial_interpolation_helpers() {
    let a = adaptive(4, &[3], 1.0);
    let old = [1.0, 2.0, 3.0, 4.0];
    let new = [5.0, 6.0, 7.0, 8.0];
    assert_eq!(a.time_interpolation(&old, &new, 0.0), old);
    assert_eq!(a.time_interpolation(&old, &new, 1.0), new);
    assert_eq!(a.time_interpolation(&old, &new, 2.0), new);
    assert_eq!(a.time_interpolation(&old, &new, 0.5), [3.0, 4.0, 5.0, 6.0]);
    let fine = a.spatial_interpolation(&old, 3);
    assert_eq!(fine.len(), 12);
    assert!((fine[0] - 1.0).abs() < 1e-15 && (fine[3] - 2.0).abs() < 1e-15);
    assert!((fine[1] - 4.0 / 3.0).abs() < 1e-15);
}

// ---------------------------------------------------------------------------
// MesoscaleLbm (F2640)
// ---------------------------------------------------------------------------

#[test]
fn mesoscale_new_cell_counts_and_defaults() {
    let m = MesoscaleLbm::new(10, 4);
    assert_eq!(m.meso_cells.len(), 10);
    assert_eq!(m.nano_cells.len(), 40);
    assert_eq!(m.n_meso(), 10);
    assert!(m.macro_density.iter().all(|&rho| rho == 1.0));
    assert!(m.macro_velocity.iter().all(|&u| u == 0.0));
    assert!((m.relaxation_time() - 2.0 / 3.0).abs() < 1e-15);
    assert!((m.kinematic_viscosity() - CS2 / 6.0).abs() < 1e-15);
}

#[test]
fn mesoscale_rest_state_is_invariant() {
    let mut m = MesoscaleLbm::new(8, 4);
    for _ in 0..20 {
        m.heterogeneous_step();
    }
    assert!(m.macro_density.iter().all(|&rho| rho == 1.0));
    assert!(m.macro_velocity.iter().all(|&u| u == 0.0));
    assert!(m.meso_cells.iter().all(|&v| v == 0.0));
}

#[test]
fn mesoscale_density_spike_splits_two_thirds_one_sixth_whatever_omega() {
    // A relaxation towards zero with ω = 1.5 would flip the sign of the spike.
    for &omega in &[0.6, 1.0, 1.5, 1.95] {
        let mut m = MesoscaleLbm::new(12, 2);
        m.set_relaxation_rate(omega)
            .expect("stable relaxation rate");
        m.macro_density[6] = 2.0;
        m.lbm_step();
        for (i, &rho) in m.macro_density.iter().enumerate() {
            let expected = match i {
                5 | 7 => 1.0 + 1.0 / 6.0,
                6 => 1.0 + 2.0 / 3.0,
                _ => 1.0,
            };
            assert!((rho - expected).abs() < 1e-14, "ω={omega} cell {i}: {rho}");
        }
        assert!((m.total_mass() - 13.0).abs() < 1e-13);
    }
}

#[test]
fn mesoscale_uniform_flow_is_preserved() {
    let mut m = MesoscaleLbm::new(16, 1);
    m.macro_density.iter_mut().for_each(|rho| *rho = 1.2);
    m.macro_velocity.iter_mut().for_each(|u| *u = 0.1);
    for _ in 0..100 {
        m.heterogeneous_step();
    }
    for (&rho, &u) in m.macro_density.iter().zip(m.macro_velocity.iter()) {
        assert!((rho - 1.2).abs() < 1e-13 && (u - 0.1).abs() < 1e-13);
    }
}

#[test]
fn mesoscale_conserves_mass_and_momentum() {
    let n = 64;
    let mut m = MesoscaleLbm::new(n, 1);
    m.set_relaxation_rate(1.2).expect("stable relaxation rate");
    for ((rho, u), (a, b)) in m
        .macro_density
        .iter_mut()
        .zip(m.macro_velocity.iter_mut())
        .zip(pseudo_random(1, n).into_iter().zip(pseudo_random(2, n)))
    {
        *rho = 0.9 + 0.2 * a;
        *u = 0.1 * (b - 0.5);
    }
    m.meso_to_macro();
    let mass = m.total_mass();
    let momentum = m.total_momentum();
    assert!(momentum.abs() > 1e-3, "the test needs a non-zero momentum");
    for _ in 0..300 {
        m.lbm_step();
    }
    assert!((m.total_mass() - mass).abs() < 1e-13, "mass drift");
    assert!(
        (m.total_momentum() - momentum).abs() < 1e-13,
        "momentum drift"
    );
}

#[test]
fn mesoscale_sound_pulses_travel_at_cs_and_broaden_with_the_lattice_viscosity() {
    let n = 512;
    let centre = 256usize;
    let omega = 1.25;
    let mut m = MesoscaleLbm::new(n, 1);
    m.set_relaxation_rate(omega)
        .expect("stable relaxation rate");
    let amplitude = 1e-4;
    for (rho, bump) in m
        .macro_density
        .iter_mut()
        .zip(gaussian(n, centre as f64, 8.0))
    {
        *rho = 1.0 + amplitude * bump;
    }
    m.meso_to_macro();
    let mass = m.total_mass();

    // Statistics of the right-going pulse: the density excess to the right of the source.
    let right_pulse = |m: &MesoscaleLbm| -> (f64, f64) {
        let excess: Vec<f64> = m.macro_density[centre + 1..]
            .iter()
            .map(|rho| rho - 1.0)
            .collect();
        centroid_and_variance(&excess, (centre + 1) as f64, 1.0)
    };
    let mut stats = Vec::new();
    for step in 1..=200 {
        m.lbm_step();
        if step == 100 || step == 200 {
            stats.push(right_pulse(&m));
        }
    }
    assert!((m.total_mass() - mass).abs() < 1e-11);

    let cs = CS2.sqrt();
    let travelled = stats[1].0 - stats[0].0;
    assert!(
        (travelled - 100.0 * cs).abs() < 0.05,
        "the pulse moved {travelled} cells in 100 steps, cs * 100 = {}",
        100.0 * cs
    );
    // Each sound pulse diffuses like a Gaussian with D = cs² (τ - 1/2).
    let rate = (stats[1].1 - stats[0].1) / 100.0;
    let expected = 2.0 * m.kinematic_viscosity();
    assert!(
        ((rate - expected) / expected).abs() < 0.01,
        "variance growth rate {rate} vs 2 ν = {expected}"
    );
    // The initial symmetric bump stays mirror symmetric.
    for k in 1..200 {
        assert!((m.macro_density[centre + k] - m.macro_density[centre - k]).abs() < 1e-13);
        assert!((m.macro_velocity[centre + k] + m.macro_velocity[centre - k]).abs() < 1e-13);
    }
}

#[test]
fn mesoscale_meso_to_macro_returns_the_moments_of_the_populations() {
    let mut m = MesoscaleLbm::new(5, 1);
    m.macro_density[2] = 1.3;
    m.macro_velocity[2] = 0.12;
    m.meso_to_macro();
    assert!((m.macro_density[2] - 1.3).abs() < 1e-14);
    assert!((m.macro_velocity[2] - 0.12).abs() < 1e-14);
    let f = m.meso_pops[2];
    let feq = hydrodynamic_equilibrium(1.3, 0.12);
    for (a, b) in f.iter().zip(feq.iter()) {
        assert!((a - b).abs() < 1e-14);
    }
    assert!((zeroth_moment(&f) - 1.3).abs() < 1e-14);
    assert!((first_moment(&f) - 1.3 * 0.12).abs() < 1e-14);
    // Second moment of the isothermal equilibrium: ρ (cs² + u²).
    let pressure_tensor = f[1] + f[2];
    assert!((pressure_tensor - 1.3 * (CS2 + 0.12 * 0.12)).abs() < 1e-14);
    // The momentum flux is not derived from a density gradient any more.
    assert!(m.macro_velocity[1].abs() < 1e-15 && m.macro_velocity[3].abs() < 1e-15);
}

#[test]
fn mesoscale_nano_values_are_averaged_and_transferred_as_mass() {
    let mut m = MesoscaleLbm::new(6, 4);
    for (k, v) in [1.0, 3.0, 5.0, 7.0].iter().enumerate() {
        m.nano_cells[2 * 4 + k] = *v;
    }
    m.nano_to_meso();
    assert!((m.meso_cells[2] - 4.0).abs() < 1e-15);
    assert!(
        m.meso_cells
            .iter()
            .enumerate()
            .all(|(i, &v)| i == 2 || v == 0.0)
    );
    assert!((m.macro_density[2] - 5.0).abs() < 1e-14);
    assert!((m.total_mass() - 10.0).abs() < 1e-13);
    // The transfer is incremental: nothing more is added for unchanged nano data.
    m.nano_to_meso();
    assert!((m.total_mass() - 10.0).abs() < 1e-13);
    // Removing the nano mass restores the reference density.
    m.nano_cells.iter_mut().for_each(|v| *v = 0.0);
    m.nano_to_meso();
    assert!((m.macro_density[2] - 1.0).abs() < 1e-13);
    assert!((m.total_mass() - 6.0).abs() < 1e-13);
}

#[test]
fn mesoscale_nano_mass_moves_with_the_local_fluid_velocity() {
    let mut m = MesoscaleLbm::new(6, 2);
    m.macro_velocity[3] = 0.2;
    m.meso_to_macro();
    m.nano_cells[6] = 0.5;
    m.nano_cells[7] = 0.5;
    m.nano_to_meso();
    assert!((m.macro_density[3] - 1.5).abs() < 1e-14);
    assert!((m.macro_velocity[3] - 0.2).abs() < 1e-14);
    assert!((m.total_momentum() - 1.5 * 0.2).abs() < 1e-14);
}

#[test]
fn mesoscale_nano_removal_never_drives_the_density_negative() {
    let mut m = MesoscaleLbm::new(4, 1);
    m.nano_cells[1] = -3.0;
    m.nano_to_meso();
    assert!(m.macro_density[1] > 0.0 && m.macro_density[1] < 1e-5);
    assert!(m.macro_velocity[1].is_finite());
    // The part that could not be removed stays pending and is given back.
    m.nano_cells[1] = 0.0;
    m.nano_to_meso();
    assert!((m.macro_density[1] - 1.0).abs() < 1e-12);
}

#[test]
fn mesoscale_heterogeneous_step_couples_the_three_scales() {
    // A nano density excess launches sound waves on the meso lattice; the macro
    // fields follow from the populations.
    let mut m = MesoscaleLbm::new(64, 2);
    m.nano_cells[64] = 2.0;
    m.nano_cells[65] = 2.0;
    for _ in 0..10 {
        m.heterogeneous_step();
    }
    assert!((m.total_mass() - (64.0 + 2.0)).abs() < 1e-12);
    assert!(m.total_momentum().abs() < 1e-12);
    // Waves travel away from the source cell 32 in both directions.
    assert!(m.macro_velocity[32 + 4] > 0.0 && m.macro_velocity[32 - 4] < 0.0);
    assert!(m.macro_density[32] < 3.0 && m.macro_density[32 + 5] > 1.0);
}

#[test]
fn mesoscale_relaxation_rate_is_validated() {
    let mut m = MesoscaleLbm::new(4, 1);
    assert!(m.set_relaxation_rate(1.25).is_ok());
    assert!((m.relaxation_time() - 0.8).abs() < 1e-15);
    assert!((m.kinematic_viscosity() - 0.1).abs() < 1e-15);
    assert!(matches!(
        m.set_relaxation_rate(2.0),
        Err(Error::UnstableRelaxationTime { .. })
    ));
    assert!(matches!(
        m.set_relaxation_rate(0.0),
        Err(Error::ParameterOutOfBounds { .. })
    ));
    assert!((m.omega - 1.25).abs() < 1e-15);
}

#[test]
fn mesoscale_nano_to_meso_average_and_macro_density_positive() {
    let mut m = MesoscaleLbm::new(2, 4);
    m.nano_cells[..4].copy_from_slice(&[1.0, 3.0, 5.0, 7.0]);
    m.nano_to_meso();
    assert!((m.meso_cells[0] - 4.0).abs() < 1e-12);
    m.meso_to_macro();
    assert!(m.macro_density.iter().all(|&d| d > 0.0));
    let mut m = MesoscaleLbm::new(5, 3);
    m.heterogeneous_step();
    assert_eq!(m.macro_velocity.len(), 5);
    assert_eq!(m.macro_density.len(), 5);
}

// ---------------------------------------------------------------------------
// HybridLbmNavier (F2640)
// ---------------------------------------------------------------------------

#[test]
fn hybrid_lbm_spike_splits_two_thirds_one_sixth_whatever_the_viscosity() {
    // A real D1Q3 step from equilibrium moves 1/6 of the momentum to each
    // neighbour independent of ν; a finite-difference Laplacian would keep
    // (1 - 2 ν dt) and send ν dt.
    for &nu in &[0.01, 0.1, 0.3] {
        let mut h = HybridLbmNavier::new(64, 4, nu);
        h.lbm_velocity[10] = 1.0;
        h.step(1.0);
        assert!((h.lbm_velocity[9] - 1.0 / 6.0).abs() < 1e-15, "ν={nu}");
        assert!((h.lbm_velocity[10] - 2.0 / 3.0).abs() < 1e-15, "ν={nu}");
        assert!((h.lbm_velocity[11] - 1.0 / 6.0).abs() < 1e-15, "ν={nu}");
        let others: f64 = h
            .lbm_velocity
            .iter()
            .enumerate()
            .filter(|(i, _)| !(9..=11).contains(i))
            .map(|(_, v)| v.abs())
            .sum();
        assert!(others < 1e-15, "ν={nu}: spurious velocity {others}");
    }
}

#[test]
fn hybrid_relaxation_time_follows_viscosity_time_step_and_cell_size() {
    let mut h = HybridLbmNavier::new(32, 2, 0.2);
    assert!((h.lbm_relaxation_time(1.0) - (0.5 + 0.6)).abs() < 1e-14);
    h.cell_size = 0.5;
    assert!((h.lbm_relaxation_time(0.05) - 0.62).abs() < 1e-14);
    h.lbm_viscosity = 0.0;
    assert!((h.lbm_relaxation_time(1.0) - 0.5).abs() < 1e-15);
}

#[test]
fn hybrid_lbm_gaussian_variance_matches_the_exact_lattice_diffusion_law() {
    let n = 300;
    let (nu, dx, dt) = (0.2, 0.5, 0.05);
    let mut h = HybridLbmNavier::new(n, 2, nu);
    h.cell_size = dx;
    let profile = gaussian(n, 70.0, 5.0);
    h.lbm_velocity.copy_from_slice(&profile);
    h.ns_velocity.copy_from_slice(&profile);
    let omega = 1.0 / h.lbm_relaxation_time(dt);
    let lbm_cells = h.lbm_region.iter().filter(|&&active| active).count();
    let (mean0, var0) = centroid_and_variance(&h.lbm_velocity[..lbm_cells], 0.0, 1.0);
    let mass0: f64 = h.lbm_velocity[..lbm_cells].iter().sum();
    let steps = 300;
    for _ in 0..steps {
        h.step(dt);
    }
    let (mean, var) = centroid_and_variance(&h.lbm_velocity[..lbm_cells], 0.0, 1.0);
    let mass: f64 = h.lbm_velocity[..lbm_cells].iter().sum();
    assert!(
        (mass - mass0).abs() < 1e-11,
        "momentum drift {}",
        mass - mass0
    );
    assert!(
        (mean - mean0).abs() < 1e-9,
        "centroid moved by {}",
        mean - mean0
    );
    let exact = exact_variance_growth(omega, steps);
    assert!(
        (var - var0 - exact).abs() < 1e-8,
        "variance growth {} vs exact {exact}",
        var - var0
    );
    // In physical units the variance grows as 2 ν t.
    let physical = 2.0 * nu * dt * steps as f64 / (dx * dx);
    assert!(((var - var0 - physical) / physical).abs() < 0.01);
}

#[test]
fn hybrid_decaying_sine_mode_matches_the_analytic_solution() {
    let n = 64;
    let nu = 0.1;
    let steps = 300;
    let length = (n - 1) as f64;
    let mode: Vec<f64> = (0..n)
        .map(|i| (std::f64::consts::PI * i as f64 / length).sin())
        .collect();
    let mut h = HybridLbmNavier::new(n, 4, nu);
    h.lbm_velocity.copy_from_slice(&mode);
    h.ns_velocity.copy_from_slice(&mode);
    for _ in 0..steps {
        h.step(1.0);
    }
    let decay = (-nu * (std::f64::consts::PI / length).powi(2) * steps as f64).exp();
    let mut worst = 0.0_f64;
    for (i, &shape) in mode.iter().enumerate() {
        let combined = if i < n / 2 {
            h.lbm_velocity[i]
        } else {
            h.ns_velocity[i]
        };
        worst = worst.max((combined - shape * decay).abs());
    }
    assert!(worst < 3e-3, "hybrid solution deviates by {worst}");
    assert_eq!(h.velocity_mismatch(), 0.0);
}

#[test]
fn hybrid_step_profile_tracks_the_finite_difference_reference_across_the_interface() {
    let n = 64;
    let (nu, steps) = (0.1, 400);
    let initial: Vec<f64> = (0..n).map(|i| if i < n / 2 { 0.0 } else { 1.0 }).collect();
    let mut h = HybridLbmNavier::new(n, 4, nu);
    h.lbm_velocity.copy_from_slice(&initial);
    h.ns_velocity.copy_from_slice(&initial);
    let mut reference = initial.clone();
    for _ in 0..steps {
        h.step(1.0);
        let previous = reference.clone();
        for i in 1..n - 1 {
            reference[i] += nu * (previous[i + 1] - 2.0 * previous[i] + previous[i - 1]);
        }
    }
    let mut worst = 0.0_f64;
    for (i, &expected) in reference.iter().enumerate() {
        let combined = if i < n / 2 {
            h.lbm_velocity[i]
        } else {
            h.ns_velocity[i]
        };
        worst = worst.max((combined - expected).abs());
    }
    assert!(
        worst < 8e-3,
        "hybrid deviates from the reference by {worst}"
    );
    // The information crossed the interface in both directions.
    assert!(h.lbm_velocity[n / 2 - 6] > 0.01);
    assert!(h.ns_velocity[n / 2 + 6] < 0.99);
}

#[test]
fn hybrid_written_velocities_and_schwarz_averages_are_imposed_on_the_lbm_populations() {
    // The LBM side must honour values written to `lbm_velocity` (initial
    // conditions, Schwarz averages): the populations are rebuilt from them.
    let mut h = HybridLbmNavier::new(40, 3, 0.1);
    for i in 0..40 {
        h.lbm_velocity[i] = 2.0;
        h.ns_velocity[i] = 0.0;
    }
    h.schwarz_coupling();
    assert_eq!(h.velocity_mismatch(), 0.0);
    let overlap_cell = (0..40)
        .find(|&i| h.lbm_region[i] && h.ns_region[i])
        .expect("the overlap is not empty");
    assert!((h.lbm_velocity[overlap_cell] - 1.0).abs() < 1e-15);
    for i in 5..8 {
        h.lbm_velocity[i] = 0.0;
    }
    h.step(1.0);
    // The first cell of the domain is a Dirichlet cell and keeps its value.
    assert!((h.lbm_velocity[0] - 2.0).abs() < 1e-14);
    // Cells 4 and 5 neighbour the zeroed block: 2/3 of their own value plus
    // 1/6 of each neighbour's value.
    assert!((h.lbm_velocity[4] - (2.0 / 3.0 * 2.0 + 1.0 / 6.0 * 2.0)).abs() < 1e-14);
    assert!((h.lbm_velocity[5] - 1.0 / 6.0 * 2.0).abs() < 1e-14);
    assert!(h.lbm_velocity[6].abs() < 1e-15);
}

#[test]
fn hybrid_old_interface_behaviour_is_preserved() {
    let h = HybridLbmNavier::new(100, 5, 0.1);
    assert_eq!(h.n_cells(), 100);
    let h = HybridLbmNavier::new(20, 2, 0.1);
    let both = (0..20)
        .filter(|&i| h.lbm_region[i] && h.ns_region[i])
        .count();
    assert!(both > 0, "Overlap must have cells active in both regions");
    let mut h = HybridLbmNavier::new(40, 3, 0.01);
    h.step(0.001);
    assert_eq!(h.lbm_velocity.len(), 40);
    assert_eq!(h.ns_velocity.len(), 40);
    assert!(h.flux_matching().is_finite());
    assert_eq!(h.velocity_mismatch(), 0.0);

    let mut h = HybridLbmNavier::new(20, 4, 0.1);
    for i in 0..20 {
        if h.lbm_region[i] {
            h.lbm_velocity[i] = 1.0;
        }
        if h.ns_region[i] {
            h.ns_velocity[i] = 0.0;
        }
    }
    let before = h.velocity_mismatch();
    h.schwarz_coupling();
    let after = h.velocity_mismatch();
    assert!(after < before, "Schwarz must reduce mismatch");
    assert_eq!(after, 0.0);
}

#[test]
fn hybrid_degenerate_domains_do_not_panic() {
    for n in 0..4 {
        let mut h = HybridLbmNavier::new(n, 1, 0.1);
        h.step(0.1);
        assert_eq!(h.n_cells(), n);
    }
}
