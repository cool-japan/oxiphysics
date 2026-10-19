//! Regression tests for the exact-kernel-gradient pair sums of the
//! pressure-solver types: the IISPH matrix, the continuity rate, the PCISPH
//! velocity correction and the divergence measures.

use super::support::{kernel_gradient, neighbor_row, separation, support_neighbors};
use super::types_advanced::IsphPressure;
use super::types_core::{
    DivergenceFreeProjector, IisphMatrixBuilder, PcisphCorrector, SphDensityIntegrator,
};
use crate::kernel::{CubicSplineKernel, SphKernel};

type Vec3 = [f64; 3];

const RHO0: f64 = 1000.0;

/// `dW/dr` of the cubic-spline kernel.
fn dw(r: f64, h: f64) -> f64 {
    CubicSplineKernel.grad_w(r, h)
}

/// `n_side³` lattice points with spacing `dx`, ordered `i`-major.
fn lattice(n_side: usize, dx: f64) -> Vec<Vec3> {
    let mut points = Vec::with_capacity(n_side * n_side * n_side);
    for i in 0..n_side {
        for j in 0..n_side {
            for k in 0..n_side {
                points.push([i as f64 * dx, j as f64 * dx, k as f64 * dx]);
            }
        }
    }
    points
}

/// Index of the centre site of [`lattice`] for odd `n_side`.
fn centre_index(n_side: usize) -> usize {
    let c = n_side / 2;
    (c * n_side + c) * n_side + c
}

/// Deterministic uniform pseudo-random numbers in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_f64()
    }

    fn vec3(&mut self, lo: f64, hi: f64) -> Vec3 {
        [self.range(lo, hi), self.range(lo, hi), self.range(lo, hi)]
    }
}

fn cloud(n: usize, seed: u64, extent: f64) -> Vec<Vec3> {
    let mut rng = Lcg(seed);
    (0..n).map(|_| rng.vec3(0.0, extent)).collect()
}

fn norm(v: Vec3) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn assert_close(actual: f64, expected: f64, rel: f64, what: &str) {
    let scale = expected.abs().max(f64::MIN_POSITIVE);
    assert!(
        (actual - expected).abs() <= rel * scale,
        "{what}: {actual} vs {expected}"
    );
}

// ---------------------------------------------------------------------------
// Kernel gradient helper
// ---------------------------------------------------------------------------

#[test]
fn kernel_gradient_vanishes_outside_the_support() {
    let h = 0.1;
    for &r in &[2.0 * h, 2.5 * h, 10.0 * h] {
        let g = kernel_gradient([r, 0.0, 0.0], r, h);
        assert_eq!(g, [0.0; 3], "r = {r}");
    }
}

#[test]
fn kernel_gradient_is_zero_for_degenerate_input() {
    assert_eq!(kernel_gradient([0.0; 3], 0.0, 0.1), [0.0; 3]);
    assert_eq!(kernel_gradient([0.05, 0.0, 0.0], f64::NAN, 0.1), [0.0; 3]);
    assert_eq!(kernel_gradient([0.05, 0.0, 0.0], 0.05, 0.0), [0.0; 3]);
    assert_eq!(kernel_gradient([0.05, 0.0, 0.0], 0.05, -0.1), [0.0; 3]);
    assert_eq!(kernel_gradient([0.05, 0.0, 0.0], 0.05, f64::NAN), [0.0; 3]);
}

#[test]
fn kernel_gradient_points_from_particle_i_towards_j() {
    let h = 0.1;
    // x_i = 0, x_j = +x: r_ij = x_i - x_j points to -x, W' < 0.
    let g = kernel_gradient([-0.05, 0.0, 0.0], 0.05, h);
    assert!(g[0] > 0.0, "gradient must point towards the neighbour");
    assert_close(g[0], -dw(0.05, h), 1e-14, "magnitude");
    assert_eq!(g[1], 0.0);
    assert_eq!(g[2], 0.0);
}

#[test]
fn kernel_gradient_matches_a_finite_difference_of_the_kernel() {
    let h = 0.12;
    let x_j = [0.31, -0.07, 0.19];
    let offsets: [Vec3; 5] = [
        [0.05, 0.02, -0.03],
        [-0.11, 0.07, 0.04],
        [0.02, -0.18, 0.09],
        [0.15, 0.1, 0.12],
        [-0.2, -0.05, 0.01],
    ];
    let eps = 1e-6 * h;
    for offset in offsets {
        let x_i = [x_j[0] + offset[0], x_j[1] + offset[1], x_j[2] + offset[2]];
        let (r_ij, r) = separation(x_i, x_j);
        assert!(
            r > 0.0 && r < 2.0 * h,
            "test point must lie inside the support"
        );
        let analytic = kernel_gradient(r_ij, r, h);
        let scale = dw(r, h).abs();
        for axis in 0..3 {
            let mut plus = x_i;
            let mut minus = x_i;
            plus[axis] += eps;
            minus[axis] -= eps;
            let (_, r_plus) = separation(plus, x_j);
            let (_, r_minus) = separation(minus, x_j);
            let numeric =
                (CubicSplineKernel.w(r_plus, h) - CubicSplineKernel.w(r_minus, h)) / (2.0 * eps);
            assert!(
                (analytic[axis] - numeric).abs() <= 1e-5 * scale,
                "axis {axis}, r = {r}: analytic {} vs finite difference {numeric}",
                analytic[axis]
            );
        }
    }
}

#[test]
fn kernel_gradient_does_not_blow_up_at_coincidence() {
    let h = 0.1;
    let tiny = 1e-6 * h;
    let near = norm(kernel_gradient([tiny, 0.0, 0.0], tiny, h));
    let mid = norm(kernel_gradient([0.5 * h, 0.0, 0.0], 0.5 * h, h));
    assert!(near.is_finite());
    assert!(near < 1e-5 * mid, "near = {near}, mid = {mid}");
}

#[test]
fn neighbour_rows_are_empty_outside_the_particle_range_and_for_degenerate_radii() {
    let positions = [[0.0; 3], [0.05, 0.0, 0.0], [5.0, 0.0, 0.0]];
    let lists = support_neighbors(&positions, 0.1);
    assert_eq!(neighbor_row(&lists, 0), &[1]);
    assert_eq!(neighbor_row(&lists, 1), &[0]);
    assert!(neighbor_row(&lists, 2).is_empty());
    assert!(neighbor_row(&lists, 3).is_empty());
    assert!(neighbor_row(&lists, usize::MAX).is_empty());
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let lists = support_neighbors(&positions, h);
        assert!(
            (0..3).all(|i| neighbor_row(&lists, i).is_empty()),
            "h = {h}"
        );
    }
}

// ---------------------------------------------------------------------------
// IisphMatrixBuilder::build_matrix
// ---------------------------------------------------------------------------

/// All-pairs reference implementation of the IISPH off-diagonal coefficients.
fn brute_force_matrix(
    builder: &IisphMatrixBuilder,
    positions: &[Vec3],
    masses: &[f64],
    densities: &[f64],
) -> Vec<Vec<(usize, f64)>> {
    let n = positions.len();
    let h = builder.h;
    let mut aij = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let (_, r) = separation(positions[i], positions[j]);
            if r >= 2.0 * h || r < 1e-14 {
                continue;
            }
            let g = dw(r, h);
            let coeff = -builder.dt
                * builder.dt
                * masses[j]
                * (masses[i] / (densities[i] * densities[i])
                    + masses[j] / (densities[j] * densities[j]))
                * g
                * g;
            aij[i].push((j, coeff));
        }
    }
    aij
}

fn pair_matrix(r: f64, h: f64) -> Vec<Vec<(usize, f64)>> {
    let builder = IisphMatrixBuilder::new(h, RHO0, 1.0e-3);
    builder.build_matrix(&[[0.0; 3], [r, 0.0, 0.0]], &[2.0, 3.0], &[1000.0, 900.0])
}

#[test]
fn iisph_matrix_uses_the_exact_kernel_gradient() {
    let (h, dt, r) = (0.1, 1.0e-3, 0.05);
    let builder = IisphMatrixBuilder::new(h, RHO0, dt);
    let positions = [[0.0; 3], [r, 0.0, 0.0]];
    let masses = [2.0, 3.0];
    let densities = [1000.0, 900.0];
    let aij = builder.build_matrix(&positions, &masses, &densities);
    let g2 = dw(r, h) * dw(r, h);
    let pair_weight = masses[0] / 1000.0_f64.powi(2) + masses[1] / 900.0_f64.powi(2);
    assert_eq!(aij[0].len(), 1);
    assert_eq!(aij[1].len(), 1);
    let (j0, c01) = aij[0][0];
    let (j1, c10) = aij[1][0];
    assert_eq!((j0, j1), (1, 0));
    assert_close(c01, -dt * dt * masses[1] * pair_weight * g2, 1e-12, "c01");
    assert_close(c10, -dt * dt * masses[0] * pair_weight * g2, 1e-12, "c10");
    assert!(c01 < 0.0 && c10 < 0.0, "off-diagonals are negative");
}

#[test]
fn iisph_matrix_has_no_entries_outside_the_support() {
    let h = 0.1;
    for &r in &[2.0 * h, 2.5 * h, 1.0] {
        let aij = pair_matrix(r, h);
        assert!(aij[0].is_empty() && aij[1].is_empty(), "r = {r}");
    }
    let inside = pair_matrix(1.99 * h, h);
    assert_eq!(inside[0].len(), 1);
}

#[test]
fn iisph_matrix_entries_vanish_as_particles_coincide() {
    let h = 0.1;
    let magnitude = |r: f64| pair_matrix(r, h)[0][0].1.abs();
    let reference = magnitude(0.5 * h);
    let mut previous = 0.0;
    for &r in &[1e-4 * h, 1e-3 * h, 1e-2 * h, 1e-1 * h, 0.5 * h] {
        let m = magnitude(r);
        assert!(
            m.is_finite() && m > previous,
            "r = {r}: {m} after {previous}"
        );
        previous = m;
    }
    assert!(magnitude(1e-4 * h) < 1e-6 * reference);
}

#[test]
fn iisph_matrix_matches_the_all_pairs_reference_on_a_cloud() {
    let h = 0.1;
    let builder = IisphMatrixBuilder::new(h, RHO0, 5.0e-4);
    let positions = cloud(250, 7, 0.5);
    let mut rng = Lcg(11);
    let masses: Vec<f64> = (0..positions.len()).map(|_| rng.range(0.5, 2.0)).collect();
    let densities: Vec<f64> = (0..positions.len())
        .map(|_| rng.range(900.0, 1100.0))
        .collect();
    let fast = builder.build_matrix(&positions, &masses, &densities);
    let slow = brute_force_matrix(&builder, &positions, &masses, &densities);
    assert_eq!(fast.len(), slow.len());
    let mut pairs = 0;
    for (i, (row_fast, row_slow)) in fast.iter().zip(slow.iter()).enumerate() {
        assert_eq!(row_fast.len(), row_slow.len(), "row {i}");
        for (&(j_f, c_f), &(j_s, c_s)) in row_fast.iter().zip(row_slow.iter()) {
            assert_eq!(j_f, j_s, "row {i}");
            assert_close(c_f, c_s, 1e-12, "coefficient");
            pairs += 1;
        }
    }
    assert!(pairs > 1000, "cloud too sparse: {pairs} pairs");
}

#[test]
fn iisph_matrix_is_symmetric_for_equal_particles_and_has_a_positive_diagonal() {
    let h = 0.1;
    let builder = IisphMatrixBuilder::new(h, RHO0, 1.0e-3);
    let positions = cloud(120, 3, 0.4);
    let masses = vec![1.5; positions.len()];
    let densities = vec![RHO0; positions.len()];
    let aij = builder.build_matrix(&positions, &masses, &densities);
    for (i, row) in aij.iter().enumerate() {
        for &(j, c_ij) in row {
            let c_ji = aij[j]
                .iter()
                .find(|&&(k, _)| k == i)
                .map(|&(_, c)| c)
                .unwrap_or(f64::NAN);
            assert_close(c_ij, c_ji, 1e-12, "symmetry");
        }
    }
    let diagonal = builder.diagonal(&aij);
    for (i, (&d, row)) in diagonal.iter().zip(aij.iter()).enumerate() {
        let row_sum: f64 = row.iter().map(|&(_, c)| c).sum();
        assert!((d + row_sum).abs() <= 1e-12 * d.abs().max(1e-30), "row {i}");
        assert!(d >= 0.0, "row {i}: diagonal {d}");
    }
}

// ---------------------------------------------------------------------------
// SphDensityIntegrator: continuity equation
// ---------------------------------------------------------------------------

fn approaching_pair(speed: f64) -> ([Vec3; 2], [Vec3; 2]) {
    (
        [[0.0; 3], [0.05, 0.0, 0.0]],
        [[speed, 0.0, 0.0], [-speed, 0.0, 0.0]],
    )
}

#[test]
fn continuity_rate_is_positive_for_approaching_particles() {
    let h = 0.1;
    let integrator = SphDensityIntegrator::new(RHO0, h);
    let (positions, velocities) = approaching_pair(1.0);
    let masses = [2.0, 2.0];
    // Dρ_0/Dt = m_1 (v_0 - v_1) . grad_0 W_01 with (v_0 - v_1) = (2, 0, 0) and
    // grad_0 W_01 = |W'| (+x) because W' < 0 and x_0 - x_1 points to -x ... the
    // neighbour approaches, so the density grows.
    let expected = masses[1] * 2.0 * dw(0.05, h).abs();
    for i in 0..2 {
        let rate = integrator.drho_dt(i, &positions, &velocities, &masses);
        assert!(rate > 0.0, "particle {i}: rate {rate}");
        assert_close(rate, expected, 1e-12, "closed form");
    }
}

#[test]
fn continuity_rate_is_negative_for_receding_particles() {
    let h = 0.1;
    let integrator = SphDensityIntegrator::new(RHO0, h);
    let (positions, velocities) = approaching_pair(-1.0);
    let masses = [2.0, 2.0];
    for i in 0..2 {
        let rate = integrator.drho_dt(i, &positions, &velocities, &masses);
        assert!(rate < 0.0, "particle {i}: rate {rate}");
    }
}

#[test]
fn continuity_rate_ignores_pairs_beyond_the_support() {
    let h = 0.1;
    let integrator = SphDensityIntegrator::new(RHO0, h);
    let positions = [[0.0; 3], [0.25, 0.0, 0.0]];
    let velocities = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
    let masses = [2.0, 2.0];
    assert_eq!(integrator.drho_dt(0, &positions, &velocities, &masses), 0.0);
    let on_the_boundary = [[0.0; 3], [2.0 * h, 0.0, 0.0]];
    assert_eq!(
        integrator.drho_dt(0, &on_the_boundary, &velocities, &masses),
        0.0
    );
}

#[test]
fn continuity_rate_of_a_uniform_expansion_is_minus_three_alpha_rho() {
    let (n_side, dx, alpha) = (7, 0.05, 0.8);
    let h = 1.3 * dx;
    let positions = lattice(n_side, dx);
    let c = positions[centre_index(n_side)];
    let velocities: Vec<Vec3> = positions
        .iter()
        .map(|p| {
            [
                alpha * (p[0] - c[0]),
                alpha * (p[1] - c[1]),
                alpha * (p[2] - c[2]),
            ]
        })
        .collect();
    let masses = vec![RHO0 * dx * dx * dx; positions.len()];
    let integrator = SphDensityIntegrator::new(RHO0, h);
    let rate = integrator.drho_dt(centre_index(n_side), &positions, &velocities, &masses);
    let expected = -3.0 * alpha * RHO0;
    assert!(
        (rate - expected).abs() < 0.05 * expected.abs(),
        "Dρ/Dt = {rate}, continuum value {expected}"
    );
}

#[test]
fn continuity_rates_of_all_particles_match_the_single_particle_queries() {
    let h = 0.1;
    let integrator = SphDensityIntegrator::new(RHO0, h);
    let positions = cloud(200, 5, 0.5);
    let mut rng = Lcg(9);
    let velocities: Vec<Vec3> = (0..positions.len()).map(|_| rng.vec3(-1.0, 1.0)).collect();
    let masses: Vec<f64> = (0..positions.len()).map(|_| rng.range(0.5, 2.0)).collect();
    let all = integrator.drho_dt_all(&positions, &velocities, &masses);
    assert_eq!(all.len(), positions.len());
    let scale = all.iter().fold(0.0_f64, |m, r| m.max(r.abs()));
    assert!(scale > 0.0);
    for (i, &rate) in all.iter().enumerate() {
        let single = integrator.drho_dt(i, &positions, &velocities, &masses);
        assert!(
            (rate - single).abs() <= 1e-12 * scale,
            "particle {i}: {rate} vs {single}"
        );
    }
}

#[test]
fn integrating_the_continuity_equation_compresses_approaching_particles() {
    let h = 0.1;
    let integrator = SphDensityIntegrator::new(RHO0, h);
    let (positions, velocities) = approaching_pair(1.0);
    let masses = [2.0, 2.0];
    let mut densities = vec![RHO0; 2];
    let rate = integrator.drho_dt(0, &positions, &velocities, &masses);
    integrator.integrate(&mut densities, &positions, &velocities, &masses, 1.0e-4);
    for &density in &densities {
        assert_close(density, RHO0 + rate * 1.0e-4, 1e-12, "Euler update");
        assert!(density > RHO0);
    }
}

// ---------------------------------------------------------------------------
// PcisphCorrector::correct_velocities
// ---------------------------------------------------------------------------

#[test]
fn pcisph_correction_pushes_a_pressurised_particle_away_from_its_neighbour() {
    let (h, dt, r) = (0.1, 1.0e-3, 0.05);
    let corrector = PcisphCorrector::new(dt, RHO0);
    let positions = [[0.0; 3], [r, 0.0, 0.0]];
    let pred_v = [[0.0; 3], [0.0; 3]];
    let pressures = [500.0, 0.0];
    let masses = [2.0, 3.0];
    let densities = [1000.0, 950.0];
    let v = corrector.correct_velocities(&positions, &pred_v, &pressures, &masses, &densities, h);
    // a_0 = -m_1 (p_0/rho_0^2) |W'| (+x): away from the neighbour at +x.
    let expected_0 = -dt * masses[1] * pressures[0] / (1000.0_f64.powi(2)) * dw(r, h).abs();
    assert_close(v[0][0], expected_0, 1e-12, "particle 0");
    assert!(v[0][0] < 0.0);
    // a_1 = +m_0 (p_0/rho_0^2) |W'|: pushed towards +x.
    let expected_1 = dt * masses[0] * pressures[0] / (1000.0_f64.powi(2)) * dw(r, h).abs();
    assert_close(v[1][0], expected_1, 1e-12, "particle 1");
    assert!(v[1][0] > 0.0);
    // The pair is aligned with x, so the other components stay untouched.
    assert_eq!(v[0][1..], [0.0, 0.0]);
    assert_eq!(v[1][1..], [0.0, 0.0]);
}

#[test]
fn pcisph_correction_is_an_acceleration_and_does_not_divide_by_the_own_mass() {
    let (h, dt) = (0.1, 1.0e-3);
    let corrector = PcisphCorrector::new(dt, RHO0);
    let positions = [[0.0; 3], [0.06, 0.01, 0.0]];
    let pred_v = [[0.0; 3], [0.0; 3]];
    let pressures = [400.0, 100.0];
    let densities = [1000.0, 1000.0];
    let light =
        corrector.correct_velocities(&positions, &pred_v, &pressures, &[1.0, 2.0], &densities, h);
    let heavy =
        corrector.correct_velocities(&positions, &pred_v, &pressures, &[7.0, 2.0], &densities, h);
    // Particle 0's acceleration depends on the neighbour mass only.
    for k in 0..3 {
        assert_close(heavy[0][k], light[0][k], 1e-14, "own mass must not enter");
    }
}

#[test]
fn pcisph_correction_conserves_momentum_for_unequal_masses() {
    let h = 0.1;
    let corrector = PcisphCorrector::new(1.0e-3, RHO0);
    let positions = cloud(200, 21, 0.5);
    let mut rng = Lcg(5);
    let n = positions.len();
    let masses: Vec<f64> = (0..n).map(|_| rng.range(0.3, 3.0)).collect();
    let densities: Vec<f64> = (0..n).map(|_| rng.range(900.0, 1100.0)).collect();
    let pressures: Vec<f64> = (0..n).map(|_| rng.range(0.0, 2000.0)).collect();
    let pred_v = vec![[0.0; 3]; n];
    let v = corrector.correct_velocities(&positions, &pred_v, &pressures, &masses, &densities, h);
    let mut total = [0.0; 3];
    let mut scale = 0.0;
    for (m, dv) in masses.iter().zip(v.iter()) {
        for k in 0..3 {
            total[k] += m * dv[k];
            scale += (m * dv[k]).abs();
        }
    }
    assert!(scale > 0.0);
    assert!(
        norm(total) <= 1e-10 * scale,
        "net momentum change {total:?} (scale {scale})"
    );
}

#[test]
fn pcisph_correction_matches_the_all_pairs_reference() {
    let h = 0.1;
    let dt = 2.0e-3;
    let corrector = PcisphCorrector::new(dt, RHO0);
    let positions = cloud(180, 33, 0.45);
    let mut rng = Lcg(77);
    let n = positions.len();
    let masses: Vec<f64> = (0..n).map(|_| rng.range(0.5, 2.0)).collect();
    let densities: Vec<f64> = (0..n).map(|_| rng.range(900.0, 1100.0)).collect();
    let pressures: Vec<f64> = (0..n).map(|_| rng.range(0.0, 1500.0)).collect();
    let pred_v: Vec<Vec3> = (0..n).map(|_| rng.vec3(-1.0, 1.0)).collect();
    let fast =
        corrector.correct_velocities(&positions, &pred_v, &pressures, &masses, &densities, h);
    let mut expected = pred_v.clone();
    let mut max_change = 0.0_f64;
    for i in 0..n {
        let mut accel = [0.0; 3];
        for j in 0..n {
            if i == j {
                continue;
            }
            let (r_ij, r) = separation(positions[i], positions[j]);
            if r >= 2.0 * h || r < 1e-14 {
                continue;
            }
            let s = dw(r, h) / r;
            let coeff = -masses[j]
                * (pressures[i] / (densities[i] * densities[i])
                    + pressures[j] / (densities[j] * densities[j]));
            for k in 0..3 {
                accel[k] += coeff * s * r_ij[k];
            }
        }
        for (k, a) in accel.iter().enumerate() {
            expected[i][k] += dt * a;
            max_change = max_change.max((dt * a).abs());
        }
    }
    assert!(
        max_change > 0.0,
        "the pressure field must change some velocity"
    );
    for i in 0..n {
        for k in 0..3 {
            assert!(
                (fast[i][k] - expected[i][k]).abs() <= 1e-12 * (1.0 + max_change),
                "particle {i} axis {k}: {} vs {}",
                fast[i][k],
                expected[i][k]
            );
        }
    }
}

#[test]
fn pcisph_correction_leaves_isolated_particles_alone() {
    let h = 0.1;
    let corrector = PcisphCorrector::new(1.0e-3, RHO0);
    let positions = [[0.0; 3], [0.5, 0.0, 0.0]];
    let pred_v = [[1.0, 2.0, 3.0], [-1.0, 0.0, 0.5]];
    let v = corrector.correct_velocities(
        &positions,
        &pred_v,
        &[1.0e5, 1.0e5],
        &[1.0, 1.0],
        &[1000.0, 1000.0],
        h,
    );
    assert_eq!(v[0], pred_v[0]);
    assert_eq!(v[1], pred_v[1]);
}

// ---------------------------------------------------------------------------
// Divergence measures
// ---------------------------------------------------------------------------

struct Expansion {
    positions: Vec<Vec3>,
    velocities: Vec<Vec3>,
    masses: Vec<f64>,
    densities: Vec<f64>,
    h: f64,
    centre: usize,
}

/// A lattice with the linear velocity field `v = alpha (x - x_centre)`.
fn expansion(alpha: f64) -> Expansion {
    let (n_side, dx) = (7, 0.05);
    let positions = lattice(n_side, dx);
    let centre = centre_index(n_side);
    let c = positions[centre];
    let velocities = positions
        .iter()
        .map(|p| {
            [
                alpha * (p[0] - c[0]),
                alpha * (p[1] - c[1]),
                alpha * (p[2] - c[2]),
            ]
        })
        .collect();
    let n = positions.len();
    Expansion {
        positions,
        velocities,
        masses: vec![RHO0 * dx * dx * dx; n],
        densities: vec![RHO0; n],
        h: 1.3 * dx,
        centre,
    }
}

#[test]
fn isph_divergence_of_an_expanding_lattice_is_three_alpha() {
    let alpha = 0.8;
    let e = expansion(alpha);
    let div = IsphPressure::divergence_velocity(
        &e.positions,
        &e.velocities,
        &e.masses,
        &e.densities,
        e.h,
        e.centre,
    );
    assert!(div > 0.0, "expansion must have positive divergence: {div}");
    assert!(
        (div - 3.0 * alpha).abs() < 0.05 * 3.0 * alpha,
        "div = {div}, continuum value {}",
        3.0 * alpha
    );
}

#[test]
fn isph_divergence_of_a_contracting_lattice_is_negative() {
    let e = expansion(-0.8);
    let div = IsphPressure::divergence_velocity(
        &e.positions,
        &e.velocities,
        &e.masses,
        &e.densities,
        e.h,
        e.centre,
    );
    assert!(
        div < 0.0,
        "contraction must have negative divergence: {div}"
    );
    assert!((div + 2.4).abs() < 0.05 * 2.4, "div = {div}");
}

#[test]
fn isph_divergence_counts_neighbours_out_to_twice_the_smoothing_length() {
    let h = 0.1;
    let positions = [[0.0; 3], [1.5 * h, 0.0, 0.0]];
    let velocities = [[0.0; 3], [1.0, 0.0, 0.0]];
    let div = IsphPressure::divergence_velocity(
        &positions,
        &velocities,
        &[1.0, 1.0],
        &[RHO0, RHO0],
        h,
        0,
    );
    // (m_1/rho_1) (v_1 - v_0) . grad_0 W_01 = (1/rho) * 1 * |W'(1.5 h)|.
    assert_close(div, dw(1.5 * h, h).abs() / RHO0, 1e-12, "pair divergence");
    assert!(div > 0.0);
}

#[test]
fn isph_divergence_of_a_rigid_rotation_vanishes() {
    let h = 0.1;
    let positions = cloud(150, 41, 0.4);
    let omega = [0.3, -0.7, 1.1];
    let velocities: Vec<Vec3> = positions
        .iter()
        .map(|p| {
            [
                omega[1] * p[2] - omega[2] * p[1],
                omega[2] * p[0] - omega[0] * p[2],
                omega[0] * p[1] - omega[1] * p[0],
            ]
        })
        .collect();
    let n = positions.len();
    let masses = vec![1.0; n];
    let densities = vec![RHO0; n];
    let all =
        IsphPressure::divergence_velocity_all(&positions, &velocities, &masses, &densities, h);
    let scale = dw(h * 2.0 / 3.0, h).abs() / RHO0 * norm(omega);
    for (i, &d) in all.iter().enumerate() {
        assert!(d.abs() <= 1e-9 * scale * n as f64, "particle {i}: {d}");
    }
}

#[test]
fn isph_divergence_of_uniform_flow_is_zero() {
    let e = expansion(0.0);
    let uniform = vec![[1.0, -2.0, 0.5]; e.positions.len()];
    let all =
        IsphPressure::divergence_velocity_all(&e.positions, &uniform, &e.masses, &e.densities, e.h);
    for (i, &d) in all.iter().enumerate() {
        assert!(d.abs() < 1e-9, "particle {i}: {d}");
    }
}

#[test]
fn isph_divergence_of_all_particles_matches_the_single_particle_queries() {
    let h = 0.1;
    let positions = cloud(160, 51, 0.45);
    let mut rng = Lcg(2);
    let n = positions.len();
    let velocities: Vec<Vec3> = (0..n).map(|_| rng.vec3(-1.0, 1.0)).collect();
    let masses: Vec<f64> = (0..n).map(|_| rng.range(0.5, 2.0)).collect();
    let densities: Vec<f64> = (0..n).map(|_| rng.range(900.0, 1100.0)).collect();
    let all =
        IsphPressure::divergence_velocity_all(&positions, &velocities, &masses, &densities, h);
    let scale = all.iter().fold(0.0_f64, |m, d| m.max(d.abs()));
    assert!(scale > 0.0);
    for (i, &d) in all.iter().enumerate() {
        let single =
            IsphPressure::divergence_velocity(&positions, &velocities, &masses, &densities, h, i);
        assert!((d - single).abs() <= 1e-12 * scale, "particle {i}");
    }
}

#[test]
fn projector_divergence_measure_is_the_negative_velocity_divergence() {
    let e = expansion(0.8);
    let projector = DivergenceFreeProjector::new(RHO0, 100, 1e-6, 0.7);
    let measure = projector.divergence_at(
        e.centre,
        &e.positions,
        &e.velocities,
        &e.masses,
        &e.densities,
        e.h,
    );
    // DFSPH convention (v_i - v_j) . grad W: negative for an expansion.
    assert!(measure < 0.0, "measure = {measure}");
    assert!((measure + 2.4).abs() < 0.05 * 2.4, "measure = {measure}");
    let compressing = expansion(-0.8);
    let converging = projector.divergence_at(
        compressing.centre,
        &compressing.positions,
        &compressing.velocities,
        &compressing.masses,
        &compressing.densities,
        compressing.h,
    );
    assert!(
        converging > 0.0,
        "converging flow must be positive: {converging}"
    );
    // It is exactly minus the ISPH divergence of the same field.
    for i in 0..e.positions.len() {
        let measure_i =
            projector.divergence_at(i, &e.positions, &e.velocities, &e.masses, &e.densities, e.h);
        let isph = IsphPressure::divergence_velocity(
            &e.positions,
            &e.velocities,
            &e.masses,
            &e.densities,
            e.h,
            i,
        );
        assert!(
            (measure_i + isph).abs() <= 1e-12 * (isph.abs() + 1.0),
            "particle {i}"
        );
    }
}

#[test]
fn projector_divergence_of_all_particles_matches_the_single_particle_queries() {
    let h = 0.1;
    let projector = DivergenceFreeProjector::new(RHO0, 100, 1e-6, 0.7);
    let positions = cloud(150, 61, 0.45);
    let mut rng = Lcg(14);
    let n = positions.len();
    let velocities: Vec<Vec3> = (0..n).map(|_| rng.vec3(-1.0, 1.0)).collect();
    let masses: Vec<f64> = (0..n).map(|_| rng.range(0.5, 2.0)).collect();
    let densities: Vec<f64> = (0..n).map(|_| rng.range(900.0, 1100.0)).collect();
    let all = projector.divergence_all(&positions, &velocities, &masses, &densities, h);
    let scale = all.iter().fold(0.0_f64, |m, d| m.max(d.abs()));
    assert!(scale > 0.0);
    for (i, &d) in all.iter().enumerate() {
        let single = projector.divergence_at(i, &positions, &velocities, &masses, &densities, h);
        assert!((d - single).abs() <= 1e-12 * scale, "particle {i}");
    }
    let max_div = projector.max_divergence(&positions, &velocities, &masses, &densities, h);
    assert!((max_div - scale).abs() <= 1e-12 * scale);
}

#[test]
fn projector_flags_an_expanding_field_and_accepts_uniform_flow() {
    let e = expansion(0.8);
    let strict = DivergenceFreeProjector::new(RHO0, 100, 1e-3, 0.7);
    assert!(!strict.is_divergence_free(&e.positions, &e.velocities, &e.masses, &e.densities, e.h));
    let uniform = vec![[2.0, 0.0, -1.0]; e.positions.len()];
    assert!(strict.is_divergence_free(&e.positions, &uniform, &e.masses, &e.densities, e.h));
    assert!(strict.max_divergence(&e.positions, &uniform, &e.masses, &e.densities, e.h) < 1e-9);
}
