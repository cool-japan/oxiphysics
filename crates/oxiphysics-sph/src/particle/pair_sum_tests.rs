//! Regression tests for the cubic-spline gradient normalisation and the
//! neighbour-list based WCSPH pair sums (density, pressure acceleration and
//! Morris viscosity) of [`SphParticleSet`].

use super::functions::{cubic_spline_gradient, cubic_spline_kernel};
use super::types::SphParticleSet;
use std::f64::consts::PI;

/// Deterministic pseudo-random number in `[0, 1)` (64-bit LCG).
fn lcg_unit(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6_364_136_223_846_793_005)
        .wrapping_add(1_442_695_040_888_963_407);
    ((*state >> 11) as f64) / ((1_u64 << 53) as f64)
}

/// Relative comparison with an absolute floor for values near zero.
fn close(a: f64, b: f64, rel: f64, abs: f64) -> bool {
    (a - b).abs() <= abs + rel * b.abs()
}

/// Cubic lattice of `n³` particles with spacing `dx`, mass `rho0 dx³` and
/// uniform density `rho0` (set directly, no density summation).
fn lattice(n: usize, dx: f64, rho0: f64) -> SphParticleSet {
    let mut ps = SphParticleSet::new();
    let mass = rho0 * dx * dx * dx;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                ps.add_raw(
                    [i as f64 * dx, j as f64 * dx, k as f64 * dx],
                    [0.0; 3],
                    mass,
                    0.0,
                    false,
                );
            }
        }
    }
    for rho in ps.densities.iter_mut() {
        *rho = rho0;
    }
    ps
}

/// Random cloud of `n` fluid particles in `[0, extent]³`.
fn random_cloud(n: usize, extent: f64, seed: u64) -> SphParticleSet {
    let mut state = seed;
    let mut ps = SphParticleSet::new();
    for _ in 0..n {
        let pos = [
            extent * lcg_unit(&mut state),
            extent * lcg_unit(&mut state),
            extent * lcg_unit(&mut state),
        ];
        let mass = 0.5 + lcg_unit(&mut state);
        ps.add_raw(pos, [0.0; 3], mass, 0.0, false);
    }
    ps
}

/// Two fluid particles a distance `r` apart on the x axis.
fn pair(r: f64, mass: f64, rho: f64, pressure: f64) -> SphParticleSet {
    let mut ps = SphParticleSet::new();
    ps.add_raw([0.0; 3], [0.0; 3], mass, 0.0, false);
    ps.add_raw([r, 0.0, 0.0], [0.0; 3], mass, 0.0, false);
    for k in 0..2 {
        ps.densities[k] = rho;
        ps.pressures[k] = pressure;
    }
    ps
}

// ---------------------------------------------------------------------------
// Cubic-spline gradient versus kernel
// ---------------------------------------------------------------------------

#[test]
fn gradient_matches_finite_difference_of_kernel() {
    for &h in &[1.0_f64, 0.07, 3.3] {
        for &q in &[0.3_f64, 0.9, 1.2, 1.8] {
            let r = q * h;
            let eps = 1.0e-6 * h;
            let numeric =
                (cubic_spline_kernel(r + eps, h) - cubic_spline_kernel(r - eps, h)) / (2.0 * eps);
            let analytic = cubic_spline_gradient(r, h);
            assert!(
                close(analytic, numeric, 1.0e-6, 0.0),
                "h={h} q={q}: gradient {analytic} vs finite difference {numeric}"
            );
        }
    }
}

#[test]
fn gradient_has_the_closed_form_prefactor_one_over_pi_h4() {
    let h = 0.37_f64;
    // q = 1: (-3 + 2.25) / (π h⁴); q = 0.5: (-1.5 + 0.5625) / (π h⁴).
    let unit = 1.0 / (PI * h.powi(4));
    assert!(close(
        cubic_spline_gradient(h, h),
        -0.75 * unit,
        1.0e-13,
        0.0
    ));
    assert!(close(
        cubic_spline_gradient(0.5 * h, h),
        -0.9375 * unit,
        1.0e-13,
        0.0
    ));
    // q = 1.5: -0.75 (0.5)² / (π h⁴).
    assert!(close(
        cubic_spline_gradient(1.5 * h, h),
        -0.1875 * unit,
        1.0e-13,
        0.0
    ));
}

#[test]
fn gradient_satisfies_the_divergence_identity() {
    // ∫ r W'(r) 4π r² dr = −3 for a kernel normalised to one.
    let h = 0.8;
    let steps = 4000;
    let dr = 2.0 * h / steps as f64;
    let integral: f64 = (0..steps)
        .map(|i| {
            let r = (i as f64 + 0.5) * dr;
            r * cubic_spline_gradient(r, h) * 4.0 * PI * r * r * dr
        })
        .sum();
    assert!((integral + 3.0).abs() < 1.0e-3, "integral = {integral}");
    // The kernel itself is normalised to one.
    let mass: f64 = (0..steps)
        .map(|i| {
            let r = (i as f64 + 0.5) * dr;
            cubic_spline_kernel(r, h) * 4.0 * PI * r * r * dr
        })
        .sum();
    assert!((mass - 1.0).abs() < 1.0e-3, "∫W = {mass}");
}

#[test]
fn gradient_is_continuous_negative_inside_and_zero_outside_the_support() {
    let h = 0.5;
    let below = cubic_spline_gradient(h * (1.0 - 1.0e-12), h);
    let above = cubic_spline_gradient(h * (1.0 + 1.0e-12), h);
    assert!(close(below, above, 1.0e-9, 0.0), "{below} vs {above}");
    for i in 1..40 {
        let r = 2.0 * h * i as f64 / 40.0;
        assert!(cubic_spline_gradient(r, h) < 0.0, "r = {r}");
    }
    assert_eq!(cubic_spline_gradient(2.0 * h + 1.0e-9, h), 0.0);
    assert_eq!(cubic_spline_gradient(10.0 * h, h), 0.0);
}

// ---------------------------------------------------------------------------
// Neighbour lists
// ---------------------------------------------------------------------------

#[test]
fn support_neighbor_lists_match_brute_force_search() {
    let ps = random_cloud(300, 1.0, 7);
    let h = 0.09;
    let lists = ps.support_neighbor_lists(h);
    assert_eq!(lists.len(), ps.len());
    for (i, list) in lists.iter().enumerate() {
        let expected: Vec<usize> = (0..ps.len())
            .filter(|&j| {
                let d = [
                    ps.positions[i][0] - ps.positions[j][0],
                    ps.positions[i][1] - ps.positions[j][1],
                    ps.positions[i][2] - ps.positions[j][2],
                ];
                j != i && d[0] * d[0] + d[1] * d[1] + d[2] * d[2] <= 4.0 * h * h
            })
            .collect();
        assert_eq!(*list, expected, "neighbours of particle {i}");
    }
}

#[test]
fn support_neighbor_lists_of_a_lattice_have_the_expected_size() {
    // h = 1.3 dx → support radius 2.6 dx: 80 lattice neighbours of an interior
    // node (all (i, j, k) with i² + j² + k² ≤ 6, excluding the origin).
    let dx = 0.01;
    let ps = lattice(9, dx, 1000.0);
    let lists = ps.support_neighbor_lists(1.3 * dx);
    let interior = 4 * 81 + 4 * 9 + 4;
    assert_eq!(lists[interior].len(), 80);
    for (i, list) in lists.iter().enumerate() {
        assert!(
            !list.contains(&i),
            "particle {i} listed as its own neighbour"
        );
        assert!(list.windows(2).all(|w| w[0] < w[1]), "list {i} not sorted");
        for &j in list {
            assert!(lists[j].contains(&i), "neighbour relation not symmetric");
        }
    }
}

#[test]
fn support_neighbor_list_size_grows_linearly_with_the_particle_count() {
    // Doubling the lattice edge multiplies the particle count by 8 and the
    // number of neighbour pairs by roughly the same factor (an O(N²) search
    // would find 64 times as many candidate pairs).
    let dx = 0.01;
    let h = 1.3 * dx;
    let small = lattice(8, dx, 1000.0);
    let large = lattice(16, dx, 1000.0);
    let pairs =
        |ps: &SphParticleSet| -> usize { ps.support_neighbor_lists(h).iter().map(Vec::len).sum() };
    let (n_small, n_large) = (small.len() as f64, large.len() as f64);
    let per_particle_small = pairs(&small) as f64 / n_small;
    let per_particle_large = pairs(&large) as f64 / n_large;
    assert!(per_particle_large >= per_particle_small);
    assert!(
        per_particle_large <= 80.0 && per_particle_small > 20.0,
        "{per_particle_small} / {per_particle_large} pairs per particle"
    );
}

#[test]
fn support_neighbor_lists_are_empty_for_invalid_smoothing_lengths() {
    let ps = random_cloud(10, 1.0, 3);
    for &h in &[0.0, -0.1, f64::NAN, f64::INFINITY] {
        let lists = ps.support_neighbor_lists(h);
        assert_eq!(lists.len(), 10);
        assert!(lists.iter().all(Vec::is_empty), "h = {h}");
    }
}

// ---------------------------------------------------------------------------
// Density summation
// ---------------------------------------------------------------------------

fn brute_force_density(ps: &SphParticleSet, h: f64) -> Vec<f64> {
    (0..ps.len())
        .map(|i| {
            (0..ps.len())
                .map(|j| {
                    let d = [
                        ps.positions[i][0] - ps.positions[j][0],
                        ps.positions[i][1] - ps.positions[j][1],
                        ps.positions[i][2] - ps.positions[j][2],
                    ];
                    let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                    ps.masses[j] * cubic_spline_kernel(r, h)
                })
                .sum()
        })
        .collect()
}

#[test]
fn neighbor_list_density_matches_the_all_pairs_sum() {
    let mut ps = random_cloud(400, 1.0, 11);
    let h = 0.08;
    let expected = brute_force_density(&ps, h);
    ps.compute_sph_densities(h);
    for (i, (&rho, &reference)) in ps.densities.iter().zip(expected.iter()).enumerate() {
        assert!(
            close(rho, reference, 1.0e-12, 0.0),
            "i={i}: {rho} vs {reference}"
        );
    }
}

#[test]
fn density_of_an_isolated_particle_is_its_self_contribution() {
    let mut ps = SphParticleSet::new();
    ps.add_raw([3.0, -1.0, 2.0], [0.0; 3], 2.5, 0.0, false);
    ps.compute_sph_densities(0.2);
    assert!(close(
        ps.densities[0],
        2.5 / (PI * 0.2_f64.powi(3)),
        1.0e-14,
        0.0
    ));
}

#[test]
fn density_with_neighbors_ignores_self_entries_and_invalid_indices() {
    let mut ps = random_cloud(40, 0.5, 5);
    let h = 0.12;
    let expected = brute_force_density(&ps, h);
    let mut lists = ps.support_neighbor_lists(h);
    for (i, list) in lists.iter_mut().enumerate() {
        list.push(i);
        list.push(ps.len() + 17);
    }
    ps.compute_sph_densities_with_neighbors(h, &lists);
    for (&rho, &reference) in ps.densities.iter().zip(expected.iter()) {
        assert!(close(rho, reference, 1.0e-12, 0.0));
    }
    // Missing lists leave only the self contribution.
    ps.compute_sph_densities_with_neighbors(h, &[]);
    for (&rho, &m) in ps.densities.iter().zip(ps.masses.iter()) {
        assert!(close(rho, m / (PI * h.powi(3)), 1.0e-14, 0.0));
    }
}

#[test]
fn density_of_a_lattice_is_close_to_the_rest_density() {
    let dx = 0.01;
    let rho0 = 1000.0;
    let mut ps = lattice(11, dx, rho0);
    ps.compute_sph_densities(1.3 * dx);
    let center = 5 * 121 + 5 * 11 + 5;
    assert!(
        (ps.densities[center] - rho0).abs() < 0.03 * rho0,
        "interior density {}",
        ps.densities[center]
    );
}

// ---------------------------------------------------------------------------
// Pressure acceleration
// ---------------------------------------------------------------------------

#[test]
fn pressure_acceleration_of_a_fluid_pair_is_applied_exactly_once() {
    // a_i = −m (2P/ρ²) W'(h) r̂_ij with r̂_ij = −x̂ for the left particle.
    let h = 0.1;
    let (m, rho, p) = (1.0, 1000.0, 1.0e3);
    let mut ps = pair(h, m, rho, p);
    ps.accumulate_pressure_acceleration(h);
    let magnitude = m * (2.0 * p / (rho * rho)) * 0.75 / (PI * h.powi(4));
    assert!(
        close(ps.accelerations[0][0], -magnitude, 1.0e-12, 0.0),
        "{} vs {}",
        ps.accelerations[0][0],
        -magnitude
    );
    assert!(close(ps.accelerations[1][0], magnitude, 1.0e-12, 0.0));
    for k in 1..3 {
        assert_eq!(ps.accelerations[0][k], 0.0);
        assert_eq!(ps.accelerations[1][k], 0.0);
    }
    // Linear momentum: Σ m a = 0.
    let net = ps.masses[0] * ps.accelerations[0][0] + ps.masses[1] * ps.accelerations[1][0];
    assert!(net.abs() < 1.0e-12 * magnitude);
}

#[test]
fn fluid_fluid_and_fluid_boundary_pairs_give_the_same_acceleration() {
    // The pressure force on a fluid particle must not depend on whether its
    // neighbour is a fluid or a boundary particle (no double counting).
    let h = 0.1;
    let mut fluid_pair = pair(h, 1.0, 1000.0, 1.0e3);
    fluid_pair.accumulate_pressure_acceleration(h);

    let mut wall_pair = pair(h, 1.0, 1000.0, 1.0e3);
    wall_pair.is_boundary[1] = true;
    wall_pair.accumulate_pressure_acceleration(h);

    for k in 0..3 {
        assert!(close(
            wall_pair.accelerations[0][k],
            fluid_pair.accelerations[0][k],
            1.0e-14,
            0.0
        ));
    }
    assert_eq!(wall_pair.accelerations[1], [0.0; 3], "walls do not move");
}

#[test]
fn pressure_acceleration_matches_an_independent_all_pairs_sum() {
    let mut ps = random_cloud(120, 0.4, 21);
    let h = 0.1;
    let mut state = 99_u64;
    for i in 0..ps.len() {
        ps.densities[i] = 900.0 + 200.0 * lcg_unit(&mut state);
        ps.pressures[i] = 5.0e3 * lcg_unit(&mut state);
    }
    for i in (0..ps.len()).step_by(7) {
        ps.is_boundary[i] = true;
    }
    let n = ps.len();
    let mut expected = vec![[0.0_f64; 3]; n];
    for (i, expected_i) in expected.iter_mut().enumerate() {
        if ps.is_boundary[i] {
            continue;
        }
        for j in 0..n {
            if i == j {
                continue;
            }
            let d = [
                ps.positions[i][0] - ps.positions[j][0],
                ps.positions[i][1] - ps.positions[j][1],
                ps.positions[i][2] - ps.positions[j][2],
            ];
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            if r > 2.0 * h {
                continue;
            }
            let dw = cubic_spline_gradient(r, h);
            let sym = ps.pressures[i] / ps.densities[i].powi(2)
                + ps.pressures[j] / ps.densities[j].powi(2);
            for (a, dk) in expected_i.iter_mut().zip(d.iter()) {
                *a += -ps.masses[j] * sym * dw * dk / r;
            }
        }
    }
    ps.accumulate_pressure_acceleration(h);
    for (i, expected_i) in expected.iter().enumerate() {
        for (k, &reference) in expected_i.iter().enumerate() {
            assert!(
                close(ps.accelerations[i][k], reference, 1.0e-10, 1.0e-12),
                "a[{i}][{k}]: {} vs {reference}",
                ps.accelerations[i][k]
            );
        }
    }
}

#[test]
fn pressure_acceleration_conserves_momentum_for_a_fluid_only_cloud() {
    let mut ps = random_cloud(150, 0.5, 31);
    let h = 0.1;
    let mut state = 5_u64;
    for i in 0..ps.len() {
        ps.densities[i] = 1000.0;
        ps.pressures[i] = 1.0e4 * lcg_unit(&mut state);
    }
    ps.accumulate_pressure_acceleration(h);
    let mut net = [0.0_f64; 3];
    let mut scale = 0.0_f64;
    for (m, a) in ps.masses.iter().zip(ps.accelerations.iter()) {
        for k in 0..3 {
            net[k] += m * a[k];
            scale += (m * a[k]).abs();
        }
    }
    assert!(scale > 0.0);
    for component in net {
        assert!(component.abs() <= 1.0e-12 * scale, "net momentum {net:?}");
    }
}

#[test]
fn uniform_pressure_in_the_interior_of_a_lattice_gives_no_net_acceleration() {
    let dx = 0.01;
    let mut ps = lattice(11, dx, 1000.0);
    for p in ps.pressures.iter_mut() {
        *p = 2.0e3;
    }
    ps.accumulate_pressure_acceleration(1.3 * dx);
    let center = 5 * 121 + 5 * 11 + 5;
    for k in 0..3 {
        assert!(
            ps.accelerations[center][k].abs() < 1.0e-9,
            "a = {:?}",
            ps.accelerations[center]
        );
    }
}

// ---------------------------------------------------------------------------
// Morris viscosity
// ---------------------------------------------------------------------------

/// Two fluid particles at distance `r` along x; particle 0 moves with `u`.
fn shear_pair(r: f64, mass: f64, rho: f64, u: f64) -> SphParticleSet {
    let mut ps = pair(r, mass, rho, 0.0);
    ps.velocities[0] = [u, 0.0, 0.0];
    ps
}

#[test]
fn morris_viscosity_matches_the_analytic_two_particle_term() {
    // ρ = 1000 kg/m³, μ = 1e-3 Pa·s → ν = 1e-6 m²/s.
    let (h, r, m, rho, mu, u) = (0.1, 0.1, 2.0, 1000.0, 1.0e-3, 0.5);
    let mut ps = shear_pair(r, m, rho, u);
    ps.accumulate_morris_viscosity(h, mu);
    // a_i = m_j · 2μ/(ρ_i ρ_j) · (dW/dr)/r · (v_i − v_j), dW/dr(q = 1) = −0.75/(π h⁴).
    let dw = -0.75 / (PI * h.powi(4));
    let expected_i = m * 2.0 * mu / (rho * rho) * dw / r * u;
    assert!(expected_i < 0.0, "the faster particle is decelerated");
    assert!(
        close(ps.accelerations[0][0], expected_i, 1.0e-12, 0.0),
        "{} vs {expected_i}",
        ps.accelerations[0][0]
    );
    assert!(close(ps.accelerations[1][0], -expected_i, 1.0e-12, 0.0));
    assert_eq!(ps.accelerations[0][1], 0.0);
    assert_eq!(ps.accelerations[0][2], 0.0);
    // Equal masses: momentum is conserved.
    assert!((ps.accelerations[0][0] + ps.accelerations[1][0]).abs() < 1.0e-18);
}

#[test]
fn morris_viscosity_scales_with_dynamic_viscosity_over_density_squared() {
    let (h, r, m, u) = (0.1, 0.1, 2.0, 0.5);
    let acceleration = |rho: f64, mu: f64| -> f64 {
        let mut ps = shear_pair(r, m, rho, u);
        ps.accumulate_morris_viscosity(h, mu);
        ps.accelerations[0][0]
    };
    let base = acceleration(1000.0, 1.0e-3);
    // Linear in the dynamic viscosity.
    assert!(close(
        acceleration(1000.0, 3.0e-3),
        3.0 * base,
        1.0e-12,
        0.0
    ));
    // For fixed particle mass the acceleration falls with ρ⁻² (a kinematic
    // reading of μ would give ρ⁻¹).
    assert!(close(
        acceleration(2000.0, 1.0e-3),
        0.25 * base,
        1.0e-12,
        0.0
    ));
    // Scaling μ and ρ together at fixed ν = μ/ρ and m ∝ ρ leaves a unchanged.
    let scaled = {
        let mut ps = shear_pair(r, 10.0 * m, 10_000.0, u);
        ps.accumulate_morris_viscosity(h, 1.0e-2);
        ps.accelerations[0][0]
    };
    assert!(close(scaled, base, 1.0e-12, 0.0));
}

#[test]
fn morris_viscosity_is_a_no_op_for_zero_viscosity_and_skips_walls() {
    let mut ps = shear_pair(0.1, 1.0, 1000.0, 1.0);
    ps.accumulate_morris_viscosity(0.1, 0.0);
    assert_eq!(ps.accelerations, vec![[0.0; 3]; 2]);

    ps.is_boundary[1] = true;
    ps.accumulate_morris_viscosity(0.1, 1.0e-3);
    assert!(ps.accelerations[0][0] < 0.0);
    assert_eq!(
        ps.accelerations[1], [0.0; 3],
        "walls receive no viscous force"
    );
}

#[test]
fn morris_viscosity_reproduces_the_laplacian_of_a_parabolic_shear_flow() {
    // u_x = A y² has ν ∇²u = 2 A ν, i.e. a_x = 2 A μ / ρ₀ at every interior
    // particle.
    let dx = 0.01;
    let rho0 = 1000.0;
    let mu = 0.1;
    let a_coeff = 50.0;
    let n = 11;
    let mut ps = lattice(n, dx, rho0);
    for i in 0..ps.len() {
        let y = ps.positions[i][1] - 5.0 * dx;
        ps.velocities[i][0] = a_coeff * y * y;
    }
    ps.accumulate_morris_viscosity(1.3 * dx, mu);
    let expected = 2.0 * a_coeff * mu / rho0;
    let center = 5 * n * n + 5 * n + 5;
    let a = ps.accelerations[center];
    // The discrete lattice sum at h = 1.3 dx reproduces the continuum Laplacian
    // to better than 1 %; the tolerance leaves room for rounding only.
    assert!(
        (a[0] - expected).abs() < 0.02 * expected,
        "a_x = {} vs 2Aμ/ρ₀ = {expected}",
        a[0]
    );
    assert!(a[1].abs() < 1.0e-9 * expected && a[2].abs() < 1.0e-9 * expected);
    // The same holds one lattice site off-centre along x (translation
    // invariance of the interior).
    let shifted = (5 + 1) * n * n + 5 * n + 5;
    assert!(close(ps.accelerations[shifted][0], a[0], 1.0e-9, 0.0));
}

#[test]
fn morris_viscosity_damps_a_velocity_jump_towards_the_mean() {
    // The viscous acceleration must always oppose the relative velocity of a
    // neighbour pair, whatever the orientation.
    let mut state = 17_u64;
    for _ in 0..20 {
        let dir = [
            lcg_unit(&mut state) - 0.5,
            lcg_unit(&mut state) - 0.5,
            lcg_unit(&mut state) - 0.5,
        ];
        let norm = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        let mut ps = SphParticleSet::new();
        ps.add_raw([0.0; 3], [0.4, -0.2, 0.7], 1.0, 0.0, false);
        ps.add_raw(
            [
                0.07 * dir[0] / norm,
                0.07 * dir[1] / norm,
                0.07 * dir[2] / norm,
            ],
            [0.0; 3],
            1.0,
            0.0,
            false,
        );
        ps.densities = vec![1000.0; 2];
        ps.accumulate_morris_viscosity(0.05, 0.5);
        let dot: f64 = (0..3)
            .map(|k| ps.accelerations[0][k] * ps.velocities[0][k])
            .sum();
        assert!(dot < 0.0, "viscosity must decelerate the relative motion");
    }
}

// ---------------------------------------------------------------------------
// Wrappers versus the neighbour-list entry points
// ---------------------------------------------------------------------------

#[test]
fn plain_methods_equal_the_with_neighbors_variants() {
    let mut a = random_cloud(80, 0.4, 41);
    let mut b = a.clone();
    let h = 0.09;
    let mut state = 13_u64;
    for i in 0..a.len() {
        let v = [
            lcg_unit(&mut state) - 0.5,
            lcg_unit(&mut state) - 0.5,
            lcg_unit(&mut state) - 0.5,
        ];
        a.velocities[i] = v;
        b.velocities[i] = v;
    }
    a.compute_sph_densities(h);
    a.update_wcsph_pressure(1000.0, 40.0, 7.0);
    a.accumulate_pressure_acceleration(h);
    a.accumulate_morris_viscosity(h, 0.05);

    let lists = b.support_neighbor_lists(h);
    b.compute_sph_densities_with_neighbors(h, &lists);
    b.update_wcsph_pressure(1000.0, 40.0, 7.0);
    b.accumulate_pressure_acceleration_with_neighbors(h, &lists);
    b.accumulate_morris_viscosity_with_neighbors(h, 0.05, &lists);

    assert_eq!(a.densities, b.densities);
    assert_eq!(a.accelerations, b.accelerations);
}
