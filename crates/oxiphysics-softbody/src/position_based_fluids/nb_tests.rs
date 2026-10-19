// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the spatial-hash neighbour search, the neighbour-list
//! (`*_nb`) variants of the PBF kernels and the Σ|∇C|² denominator of the
//! PBF multiplier.
//!
//! The expected values are either computed by an independent exhaustive scan
//! inside the test, derived in closed form by hand, or follow from a
//! conservation law of the scheme.

use super::*;
use std::f64::consts::PI;

const H: f64 = 0.1;
const RHO0: f64 = 1000.0;

// ------------------------------------------------------------------
// Test helpers
// ------------------------------------------------------------------

/// Deterministic linear congruential generator (reproducible clouds).
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Uniform sample in `[0, 1)`.
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// Uniform sample in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }
}

/// `n` points uniformly distributed in the cube `[centre − side/2, centre + side/2)³`.
fn random_cloud(n: usize, seed: u64, side: f64, centre: [f64; 3]) -> Vec<[f64; 3]> {
    let mut rng = Lcg::new(seed);
    (0..n)
        .map(|_| {
            [
                centre[0] + rng.range(-0.5 * side, 0.5 * side),
                centre[1] + rng.range(-0.5 * side, 0.5 * side),
                centre[2] + rng.range(-0.5 * side, 0.5 * side),
            ]
        })
        .collect()
}

/// Exhaustive O(n²) neighbour lists with the criterion `|r|² < h²`.
fn reference_lists(points: &[[f64; 3]], h: f64) -> Vec<Vec<usize>> {
    let h2 = h * h;
    points
        .iter()
        .enumerate()
        .map(|(i, pi)| {
            points
                .iter()
                .enumerate()
                .filter(|&(j, pj)| {
                    let d = [pi[0] - pj[0], pi[1] - pj[1], pi[2] - pj[2]];
                    j != i && d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < h2
                })
                .map(|(j, _)| j)
                .collect()
        })
        .collect()
}

/// Extended particles at `positions` (with `pred_pos` slightly displaced from
/// `pos` so that tests catch code that reads the wrong position field) and
/// masses uniformly drawn from `mass_range`.
fn ext_particles(positions: &[[f64; 3]], seed: u64, mass_range: (f64, f64)) -> Vec<PbfExtParticle> {
    let mut rng = Lcg::new(seed);
    positions
        .iter()
        .map(|&pos| {
            let mut p = PbfExtParticle::new(pos, rng.range(mass_range.0, mass_range.1));
            p.pred_pos = [
                pos[0] + rng.range(-0.01, 0.01),
                pos[1] + rng.range(-0.01, 0.01),
                pos[2] + rng.range(-0.01, 0.01),
            ];
            p
        })
        .collect()
}

fn pred_positions(particles: &[PbfExtParticle]) -> Vec<[f64; 3]> {
    particles.iter().map(|p| p.pred_pos).collect()
}

fn assert_close(a: f64, b: f64, rel: f64, what: &str) {
    let tol = rel * (1.0 + a.abs().max(b.abs()));
    assert!((a - b).abs() <= tol, "{what}: {a} vs {b} (tol {tol})");
}

fn assert_vec_close(a: [f64; 3], b: [f64; 3], rel: f64, what: &str) {
    for k in 0..3 {
        assert_close(a[k], b[k], rel, &format!("{what}[{k}]"));
    }
}

fn max_position_difference(a: &[PbfExtParticle], b: &[PbfExtParticle]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(p, q)| {
            (0..3)
                .map(|k| (p.pos[k] - q.pos[k]).abs())
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max)
}

/// Total linear momentum Σ m v of an extended particle set.
fn ext_momentum(particles: &[PbfExtParticle]) -> [f64; 3] {
    particles.iter().fold([0.0; 3], |acc, p| {
        [
            acc[0] + p.mass * p.vel[0],
            acc[1] + p.mass * p.vel[1],
            acc[2] + p.mass * p.vel[2],
        ]
    })
}

/// Cubic lattice with `n` particles per side centred at the origin.
fn lattice(n: usize, spacing: f64) -> Vec<[f64; 3]> {
    let offset = 0.5 * (n as f64 - 1.0) * spacing;
    let mut pts = Vec::with_capacity(n * n * n);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                pts.push([
                    ix as f64 * spacing - offset,
                    iy as f64 * spacing - offset,
                    iz as f64 * spacing - offset,
                ]);
            }
        }
    }
    pts
}

// ------------------------------------------------------------------
// Neighbour search (F4280)
// ------------------------------------------------------------------

#[test]
fn neighbor_lists_match_exhaustive_scan_on_random_clouds() {
    // (n, seed, side, centre, h): includes negative coordinates, a cloud far
    // from the origin and a sparse / dense pair of densities.
    let cases: [(usize, u64, f64, [f64; 3], f64); 5] = [
        (2, 1, 0.1, [0.0; 3], 0.1),
        (60, 7, 0.3, [0.0; 3], 0.1),
        (400, 11, 0.5, [-0.25, 0.1, -0.4], 0.12),
        (300, 23, 0.4, [-1000.3, 2000.7, -0.001], 0.15),
        (500, 5, 2.0, [0.0; 3], 0.05),
    ];
    for (n, seed, side, centre, h) in cases {
        let pts = random_cloud(n, seed, side, centre);
        let hashed = NeighborCsr::from_positions(&pts, h);
        let exhaustive = NeighborCsr::from_positions_all_pairs(&pts, h);
        assert_eq!(hashed, exhaustive, "n={n} seed={seed}");
        assert_eq!(hashed.to_lists(), reference_lists(&pts, h), "n={n}");
        assert_eq!(hashed.len(), n);
    }
    // The dense cases are not vacuous.
    let dense = NeighborCsr::from_positions(&random_cloud(400, 11, 0.5, [-0.25, 0.1, -0.4]), 0.12);
    assert!(dense.pair_count() > 400, "pairs: {}", dense.pair_count());
}

#[test]
fn neighbor_relation_is_symmetric_irreflexive_and_sorted() {
    let pts = random_cloud(350, 99, 0.35, [0.2, -0.1, 0.05]);
    let csr = NeighborCsr::from_positions(&pts, H);
    let mut directed_pairs = 0usize;
    for i in 0..csr.len() {
        let list = csr.neighbors(i);
        assert!(list.windows(2).all(|w| w[0] < w[1]), "list {i} not sorted");
        for &j in list {
            assert_ne!(i, j, "particle {i} lists itself");
            assert!(
                csr.neighbors(j).binary_search(&i).is_ok(),
                "{j} is a neighbour of {i} but not vice versa"
            );
            directed_pairs += 1;
        }
    }
    assert_eq!(directed_pairs, csr.pair_count());
    assert_eq!(directed_pairs % 2, 0);
    assert!(directed_pairs > 0);
}

#[test]
fn neighbor_criterion_is_strictly_less_than_h() {
    // Three collinear points exactly h apart (0.25 is exactly representable,
    // so |r|² == h² holds with no rounding): nobody is a neighbour.
    let pts = [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0], [0.5, 0.0, 0.0]];
    let touching = NeighborCsr::from_positions(&pts, 0.25);
    assert_eq!(touching.pair_count(), 0);
    assert_eq!(touching, NeighborCsr::from_positions_all_pairs(&pts, 0.25));
    // Marginally larger radius: the two adjacent pairs appear, the outer pair
    // (distance 0.5) stays out.
    let inside = NeighborCsr::from_positions(&pts, 0.250_000_1);
    assert_eq!(inside.to_lists(), vec![vec![1], vec![0, 2], vec![1]]);
}

#[test]
fn coincident_particles_are_neighbours() {
    let pts = [[0.3, -0.2, 0.1], [0.3, -0.2, 0.1], [5.0, 5.0, 5.0]];
    let csr = NeighborCsr::from_positions(&pts, H);
    assert_eq!(csr.to_lists(), vec![vec![1], vec![0], vec![]]);
}

#[test]
fn non_finite_particles_are_isolated() {
    let finite = random_cloud(70, 31, 0.3, [0.0; 3]);
    let mut all = finite.clone();
    let bad_nan = all.len();
    all.push([f64::NAN, 0.0, 0.0]);
    let bad_inf = all.len();
    all.push([0.0, f64::INFINITY, 0.0]);
    let bad_neg_inf = all.len();
    all.push([0.0, 0.0, f64::NEG_INFINITY]);

    let with_bad = NeighborCsr::from_positions(&all, H);
    assert_eq!(with_bad, NeighborCsr::from_positions_all_pairs(&all, H));
    assert_eq!(with_bad.len(), all.len());
    for bad in [bad_nan, bad_inf, bad_neg_inf] {
        assert!(with_bad.neighbors(bad).is_empty(), "bad particle {bad}");
    }
    // The finite particles see exactly the neighbours they see without the
    // non-finite ones present.
    let without_bad = NeighborCsr::from_positions(&finite, H);
    for i in 0..finite.len() {
        assert_eq!(
            with_bad.neighbors(i),
            without_bad.neighbors(i),
            "particle {i}"
        );
    }
    // A cloud made only of non-finite points is valid and has no pairs.
    let none = NeighborCsr::from_positions(&[[f64::NAN; 3], [f64::INFINITY; 3]], H);
    assert_eq!(none.pair_count(), 0);
    assert_eq!(none.len(), 2);
}

#[test]
fn invalid_radius_and_tiny_inputs_give_empty_lists() {
    let pts = random_cloud(40, 3, 0.2, [0.0; 3]);
    for h in [0.0, -0.1, f64::NAN] {
        let csr = NeighborCsr::from_positions(&pts, h);
        assert_eq!(csr.len(), pts.len());
        assert_eq!(csr.pair_count(), 0, "h={h}");
        assert_eq!(csr, NeighborCsr::from_positions_all_pairs(&pts, h));
    }
    let empty = NeighborCsr::from_positions(&[], H);
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert_eq!(empty.pair_count(), 0);
    assert!(empty.neighbors(0).is_empty());
    let single = NeighborCsr::from_positions(&[[0.0; 3]], H);
    assert_eq!(single.len(), 1);
    assert!(single.neighbors(0).is_empty());
    assert!(single.neighbors(17).is_empty(), "out-of-range index");
    assert!(NeighborCsr::default().is_empty());
}

#[test]
fn huge_extent_falls_back_to_the_exhaustive_scan() {
    // Extent / h = 1e15: far beyond what the hash path trusts. The close pairs
    // at both ends must still be found.
    let pts = [
        [0.0, 0.0, 0.0],
        [5.0e-4, 0.0, 0.0],
        [1.0e12, 0.0, 0.0],
        [1.0e12 + 5.0e-4, 0.0, 0.0],
    ];
    let csr = NeighborCsr::from_positions(&pts, 1.0e-3);
    assert_eq!(csr.to_lists(), vec![vec![1], vec![0], vec![3], vec![2]]);
    assert_eq!(csr, NeighborCsr::from_positions_all_pairs(&pts, 1.0e-3));
}

#[test]
fn find_neighbors_uses_predicted_positions() {
    // Positions are far apart, predicted positions are close (and vice versa).
    let mut a = PbfParticle::new([0.0, 0.0, 0.0], 1.0);
    let mut b = PbfParticle::new([1.0, 0.0, 0.0], 1.0);
    let mut c = PbfParticle::new([0.04, 0.0, 0.0], 1.0);
    a.predicted = [0.0, 0.0, 0.0];
    b.predicted = [0.05, 0.0, 0.0];
    c.predicted = [0.9, 0.0, 0.0];
    let lists = find_neighbors(&[a, b, c], H);
    assert_eq!(lists, vec![vec![1], vec![0], vec![]]);
}

#[test]
fn find_neighbors_matches_exhaustive_scan_on_random_predictions() {
    let pts = random_cloud(260, 41, 0.3, [0.0; 3]);
    let mut rng = Lcg::new(5);
    let particles: Vec<PbfParticle> = pts
        .iter()
        .map(|&p| {
            let mut q = PbfParticle::new(p, 1.0);
            q.predicted = [
                p[0] + rng.range(-0.05, 0.05),
                p[1] + rng.range(-0.05, 0.05),
                p[2] + rng.range(-0.05, 0.05),
            ];
            q
        })
        .collect();
    let predicted: Vec<[f64; 3]> = particles.iter().map(|p| p.predicted).collect();
    assert_eq!(
        find_neighbors(&particles, H),
        reference_lists(&predicted, H)
    );
}

#[test]
fn lattice_neighbour_counts_match_analytic_enumeration() {
    // 24³ = 13 824 particles on a lattice of spacing h/2.5. A particle sees
    // the lattice offsets (a, b, c) ≠ 0 with a² + b² + c² < 6.25, i.e. ≤ 6
    // (no integer sits near the cutoff, so rounding cannot matter). The
    // lattice is shifted off the cell grid on purpose.
    let n = 24usize;
    let spacing = H / 2.5;
    let shift = [-0.371, 0.113, 2.903];
    let pts: Vec<[f64; 3]> = lattice(n, spacing)
        .into_iter()
        .map(|p| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])
        .collect();

    let mut offsets: Vec<[i64; 3]> = Vec::new();
    for a in -2i64..=2 {
        for b in -2i64..=2 {
            for c in -2i64..=2 {
                if (a, b, c) != (0, 0, 0) && a * a + b * b + c * c <= 6 {
                    offsets.push([a, b, c]);
                }
            }
        }
    }
    assert_eq!(offsets.len(), 80, "full neighbour shell");

    let csr = NeighborCsr::from_positions(&pts, H);
    let index = |ix: i64, iy: i64, iz: i64| ((ix as usize) * n + iy as usize) * n + iz as usize;
    let range = 0..n as i64;
    let mut expected_pairs = 0usize;
    let mut interior = 0usize;
    for ix in 0..n as i64 {
        for iy in 0..n as i64 {
            for iz in 0..n as i64 {
                let mut expected: Vec<usize> = offsets
                    .iter()
                    .map(|o| (ix + o[0], iy + o[1], iz + o[2]))
                    .filter(|(x, y, z)| range.contains(x) && range.contains(y) && range.contains(z))
                    .map(|(x, y, z)| index(x, y, z))
                    .collect();
                expected.sort_unstable();
                expected_pairs += expected.len();
                let me = index(ix, iy, iz);
                assert_eq!(csr.neighbors(me), expected.as_slice(), "({ix},{iy},{iz})");
                if expected.len() == offsets.len() {
                    interior += 1;
                }
            }
        }
    }
    assert_eq!(csr.pair_count(), expected_pairs);
    assert_eq!(interior, (n - 4) * (n - 4) * (n - 4));
}

// ------------------------------------------------------------------
// Neighbour-list kernels equal the all-pairs kernels (F4283)
// ------------------------------------------------------------------

#[test]
fn nb_kernels_match_all_pairs_with_non_uniform_masses() {
    let pts = random_cloud(300, 17, 0.3, [0.1, -0.2, 0.3]);
    let particles = ext_particles(&pts, 18, (0.5, 2.0));
    let csr = NeighborCsr::from_positions(&pred_positions(&particles), H);
    assert!(csr.pair_count() > 3000, "pairs: {}", csr.pair_count());

    let mut rng = Lcg::new(77);
    let lambdas: Vec<f64> = (0..particles.len())
        .map(|_| rng.range(-0.02, 0.02))
        .collect();
    for i in 0..particles.len() {
        let nbrs = csr.neighbors(i);
        assert_close(
            pbf_density(&particles, i, H),
            pbf_density_nb(&particles, i, nbrs, H),
            1e-12,
            "density",
        );
        assert_close(
            pbf_grad_sum(&particles, i, H, RHO0),
            pbf_grad_sum_nb(&particles, i, nbrs, H, RHO0),
            1e-12,
            "grad_sum",
        );
        assert_vec_close(
            pbf_delta_pos(&particles, i, &lambdas, H, RHO0),
            pbf_delta_pos_nb(&particles, i, nbrs, &lambdas, H, RHO0),
            1e-12,
            "delta_pos",
        );
    }
}

#[test]
fn nb_kernels_ignore_self_and_distant_candidates() {
    // Feeding the *complete* index set as the candidate list must give the
    // all-pairs result: self entries are skipped and the r² < h² test rejects
    // the far particles.
    let pts = random_cloud(120, 8, 0.4, [0.0; 3]);
    let particles = ext_particles(&pts, 9, (1.0, 1.0));
    let everyone: Vec<usize> = (0..particles.len()).collect();
    let lambdas = vec![-0.004; particles.len()];
    for i in (0..particles.len()).step_by(7) {
        assert_close(
            pbf_density(&particles, i, H),
            pbf_density_nb(&particles, i, &everyone, H),
            1e-12,
            "density",
        );
        assert_close(
            pbf_grad_sum(&particles, i, H, RHO0),
            pbf_grad_sum_nb(&particles, i, &everyone, H, RHO0),
            1e-12,
            "grad_sum",
        );
        assert_vec_close(
            pbf_delta_pos(&particles, i, &lambdas, H, RHO0),
            pbf_delta_pos_nb(&particles, i, &everyone, &lambdas, H, RHO0),
            1e-12,
            "delta_pos",
        );
    }
}

#[test]
fn nb_kernels_are_total_on_degenerate_input() {
    let particles = ext_particles(&random_cloud(5, 1, 0.05, [0.0; 3]), 2, (1.0, 1.0));
    let near: Vec<usize> = (0..5).collect();
    // Unknown particle index, non-positive / NaN radius, missing lambdas.
    assert_eq!(pbf_density_nb(&particles, 99, &near, H), 0.0);
    assert_eq!(pbf_grad_sum_nb(&particles, 99, &near, H, RHO0), 0.0);
    for h in [0.0, -H, f64::NAN] {
        assert_eq!(pbf_density_nb(&particles, 0, &near, h), 0.0, "h={h}");
        assert_eq!(pbf_grad_sum_nb(&particles, 0, &near, h, RHO0), 0.0, "h={h}");
        assert_eq!(
            pbf_delta_pos_nb(&particles, 0, &near, &[0.1; 5], h, RHO0),
            [0.0; 3],
            "h={h}"
        );
    }
    assert_eq!(
        pbf_delta_pos_nb(&particles, 0, &near, &[], H, RHO0),
        [0.0; 3],
        "no lambdas"
    );
    // Neighbour indices outside the particle set are skipped, not a panic.
    assert_close(
        pbf_density_nb(&particles, 0, &[1, 2, 500], H),
        pbf_density_nb(&particles, 0, &[1, 2], H),
        1e-14,
        "out-of-range neighbour",
    );
    assert!(pbf_density_nb(&[], 0, &[], H) == 0.0);
}

// ------------------------------------------------------------------
// Σ_k |∇_{p_k} C_i|² and λ (F4284)
// ------------------------------------------------------------------

/// `∇W_spiky` (x component) at a separation `d` along +x seen from the
/// particle at the *origin*, i.e. for `r_vec = x_i − x_j = (−d, 0, 0)`:
/// `+45 (h − d)² / (π h⁶)`. Derived by hand from the Spiky gradient
/// `−45/(π h⁶) (h − r)² r_vec / r`.
fn spiky_gradient_x_of_origin_particle(d: f64, h: f64) -> f64 {
    45.0 * (h - d) * (h - d) / (PI * h.powi(6))
}

#[test]
fn grad_sum_of_a_pair_matches_the_closed_form_and_carries_the_neighbour_mass() {
    let d = 0.05;
    let (m_a, m_b) = (1.0, 2.0);
    let mut a = PbfExtParticle::new([0.0, 0.0, 0.0], m_a);
    let mut b = PbfExtParticle::new([d, 0.0, 0.0], m_b);
    a.pred_pos = a.pos;
    b.pred_pos = b.pos;
    let particles = [a, b];

    // For particle i the only neighbour is j: g = (m_j/ρ0) ∇W, the sum is
    // |g|² (k = i) + |g|² (k = j) = 2 |g|².
    let w = spiky_gradient_x_of_origin_particle(d, H);
    let g_a = m_b * w / RHO0; // ∇_{p_A} C_A, neighbour B has mass m_b
    let g_b = m_a * w / RHO0; // ∇_{p_B} C_B, neighbour A has mass m_a
    let expected_a = 2.0 * g_a * g_a;
    let expected_b = 2.0 * g_b * g_b;
    assert_close(pbf_grad_sum(&particles, 0, H, RHO0), expected_a, 1e-12, "A");
    assert_close(pbf_grad_sum(&particles, 1, H, RHO0), expected_b, 1e-12, "B");
    assert_close(
        pbf_grad_sum_nb(&particles, 0, &[1], H, RHO0),
        expected_a,
        1e-12,
        "A (nb)",
    );
    // The mass ratio enters squared.
    assert_close(expected_a / expected_b, (m_b / m_a).powi(2), 1e-12, "ratio");
    assert!(expected_a > 1.0, "non-trivial magnitude: {expected_a}");
}

#[test]
fn grad_sum_of_a_symmetric_triple_matches_the_closed_form() {
    // Particle 0 at the origin with neighbours at ±d along x (unit masses):
    // the k = i gradients cancel, leaving Σ_j |g_j|² = 2 g².
    let d = 0.04;
    let mut particles = vec![
        PbfExtParticle::new([0.0, 0.0, 0.0], 1.0),
        PbfExtParticle::new([d, 0.0, 0.0], 1.0),
        PbfExtParticle::new([-d, 0.0, 0.0], 1.0),
    ];
    for p in &mut particles {
        p.pred_pos = p.pos;
    }
    let g = spiky_gradient_x_of_origin_particle(d, H) / RHO0;
    assert_close(
        pbf_grad_sum(&particles, 0, H, RHO0),
        2.0 * g * g,
        1e-12,
        "centre",
    );
}

#[test]
fn grad_sum_vanishes_without_neighbours() {
    let lone = [PbfExtParticle::new([0.0; 3], 1.0)];
    assert_eq!(pbf_grad_sum(&lone, 0, H, RHO0), 0.0);
    let mut far = vec![
        PbfExtParticle::new([0.0; 3], 1.0),
        PbfExtParticle::new([2.0 * H, 0.0, 0.0], 1.0),
    ];
    for p in &mut far {
        p.pred_pos = p.pos;
    }
    assert_eq!(pbf_grad_sum(&far, 0, H, RHO0), 0.0);
}

#[test]
fn lambda_nb_equals_lambda_compute_for_unit_masses() {
    let pts = random_cloud(200, 55, 0.3, [0.0; 3]);
    let ext = ext_particles(&pts, 56, (1.0, 1.0));
    let classic: Vec<PbfParticle> = ext
        .iter()
        .map(|p| {
            let mut q = PbfParticle::new(p.pos, p.mass);
            q.predicted = p.pred_pos;
            q
        })
        .collect();
    let lists = find_neighbors(&classic, H);
    let csr = NeighborCsr::from_positions(&pred_positions(&ext), H);
    let eps = 1.0e-4;
    let mut with_neighbours = 0usize;
    for (i, list) in lists.iter().enumerate() {
        assert_eq!(list.as_slice(), csr.neighbors(i));
        let reference = lambda_compute(&classic, list, i, H, RHO0, eps);
        let candidate = pbf_lambda_nb(&ext, i, csr.neighbors(i), H, RHO0, eps);
        assert_close(reference, candidate, 1e-12, "lambda");
        if !list.is_empty() {
            with_neighbours += 1;
        }
    }
    assert!(with_neighbours > 150);
}

#[test]
fn lambda_denominator_is_the_gradient_sum_not_a_density_proxy() {
    // Regression: PbfSolver once used `ρ_i/ρ0² + 1e-6` as the denominator of
    // λ. Compare the stored multipliers of a real step with the definition.
    let pts = lattice(5, 0.04);
    let mut solver = PbfSolver::new(RHO0, H, 1);
    for (idx, p) in ext_particles(&pts, 3, (1.0, 1.4)).into_iter().enumerate() {
        let mut q = PbfExtParticle::new(pts[idx], p.mass);
        q.vel = [0.0; 3];
        solver.add_particle(q);
    }
    let before = solver.particles.clone();
    solver.step(1.0e-3, [0.0; 3]);

    // Reproduce the single projection iteration by hand.
    let mut predicted = before.clone();
    for p in &mut predicted {
        p.predict([0.0; 3], 1.0e-3);
    }
    let mut checked = 0usize;
    for i in 0..predicted.len() {
        let rho = pbf_density(&predicted, i, H);
        let grad_sum = pbf_grad_sum(&predicted, i, H, RHO0);
        let expected = -(rho / RHO0 - 1.0) / (grad_sum + solver.epsilon);
        assert_close(solver.particles[i].lambda, expected, 1e-12, "stored lambda");
        assert_close(solver.particles[i].rho, rho, 1e-12, "stored density");
        let fabricated = -(rho / RHO0 - 1.0) / (rho / (RHO0 * RHO0) + 1.0e-6 + solver.epsilon);
        assert!(
            (expected - fabricated).abs() > 0.5 * expected.abs(),
            "particle {i}: λ={expected} is indistinguishable from the old proxy {fabricated}"
        );
        checked += 1;
    }
    assert_eq!(checked, 125);
    // The gradient sum of a particle with neighbours is orders of magnitude
    // above the proxy, which is < 1e-2 here.
    assert!(pbf_grad_sum(&predicted, 62, H, RHO0) > 1.0e3);
}

#[test]
fn stored_multipliers_follow_the_definition_for_every_epsilon() {
    // λ_i = −C_i / (Σ_k |∇C_i|² + ε): a larger ε shrinks every multiplier of a
    // particle that has neighbours, and the stored values match the formula.
    let pts = random_cloud(90, 12, 0.16, [0.0; 3]);
    let solve_once = |epsilon: f64| {
        let mut solver = PbfSolver::new(RHO0, H, 1);
        solver.epsilon = epsilon;
        for p in ext_particles(&pts, 13, (0.8, 1.2)) {
            solver.add_particle(p);
        }
        let before = solver.particles.clone();
        solver.step(1.0e-3, [0.0; 3]);
        (before, solver)
    };
    let epsilons = [1.0e-4, 1.0, 1.0e3, 1.0e8];
    let mut previous: Option<Vec<f64>> = None;
    for epsilon in epsilons {
        let (before, solver) = solve_once(epsilon);
        let mut predicted = before.clone();
        for p in &mut predicted {
            p.predict([0.0; 3], 1.0e-3);
        }
        let lambdas: Vec<f64> = solver.particles.iter().map(|p| p.lambda).collect();
        let mut with_neighbours = 0usize;
        for (i, &lambda) in lambdas.iter().enumerate() {
            let rho = pbf_density(&predicted, i, H);
            let grad_sum = pbf_grad_sum(&predicted, i, H, RHO0);
            assert_close(
                lambda,
                -(rho / RHO0 - 1.0) / (grad_sum + epsilon),
                1e-12,
                "stored lambda",
            );
            if grad_sum > 0.0 {
                with_neighbours += 1;
            }
        }
        assert!(
            with_neighbours > 60,
            "particles with neighbours: {with_neighbours}"
        );
        if let Some(prev) = &previous {
            for (i, (&now, &then)) in lambdas.iter().zip(prev).enumerate() {
                assert!(
                    now.abs() <= then.abs() * (1.0 + 1e-12),
                    "particle {i}: {then} -> {now}"
                );
            }
        }
        previous = Some(lambdas);
    }
}

#[test]
fn default_epsilon_scales_with_the_inverse_square_of_the_smoothing_radius() {
    // κ / h² with κ = 6 gives Macklin & Müller's ε = 600 at h = 0.1.
    assert_close(default_cfm_epsilon(0.1), 600.0, 1e-12, "paper value");
    assert_close(
        default_cfm_epsilon(0.05),
        4.0 * default_cfm_epsilon(0.1),
        1e-12,
        "scaling",
    );
    for h in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert_eq!(default_cfm_epsilon(h), 0.0, "h={h}");
    }
    assert_eq!(PbfSolver::new(RHO0, H, 1).epsilon, default_cfm_epsilon(H));
    assert_eq!(
        PbfSolver::new(RHO0, 0.25, 3).epsilon,
        default_cfm_epsilon(0.25)
    );
}

#[test]
fn multiplier_is_zero_when_the_denominator_vanishes() {
    // No neighbours and no relaxation: −C/0 must not leak ±∞ / NaN.
    assert_eq!(pbf_lambda(1500.0, RHO0, 0.0, 0.0), 0.0);
    assert_eq!(pbf_lambda(500.0, RHO0, 0.0, 0.0), 0.0);
    assert_eq!(pbf_lambda(1500.0, RHO0, f64::NAN, 1.0), 0.0);
    let lone = [PbfParticle::new([0.0; 3], 1.0)];
    assert_eq!(lambda_compute(&lone, &[], 0, H, RHO0, 0.0), 0.0);
    // With relaxation the multiplier is the plain ratio.
    assert_close(
        pbf_lambda(1500.0, RHO0, 0.0, 2.0),
        -0.25,
        1e-14,
        "isolated, ε = 2",
    );
}

// ------------------------------------------------------------------
// Position correction with masses (F4281 overlap, closed form)
// ------------------------------------------------------------------

#[test]
fn delta_position_of_a_pair_matches_the_closed_form() {
    let d = 0.06;
    let (m_a, m_b) = (1.5, 0.5);
    let (lam_a, lam_b) = (-0.002, 0.003);
    let mut a = PbfParticle::new([0.0, 0.0, 0.0], m_a);
    let mut b = PbfParticle::new([d, 0.0, 0.0], m_b);
    a.lambda = lam_a;
    b.lambda = lam_b;
    let w = spiky_gradient_x_of_origin_particle(d, H);
    // Δp_A = (m_B λ_A + m_A λ_B) ∇W(x_A − x_B) / ρ0, Δp_B has the opposite sign.
    let expected_a = (m_b * lam_a + m_a * lam_b) * w / RHO0;
    let particles = [a, b];
    let dp_a = delta_position(&particles, &[1], 0, H, RHO0);
    let dp_b = delta_position(&particles, &[0], 1, H, RHO0);
    assert_vec_close(dp_a, [expected_a, 0.0, 0.0], 1e-12, "Δp_A");
    assert_vec_close(dp_b, [-expected_a, 0.0, 0.0], 1e-12, "Δp_B");

    let mut ext_a = PbfExtParticle::new([0.0; 3], m_a);
    let mut ext_b = PbfExtParticle::new([d, 0.0, 0.0], m_b);
    ext_a.pred_pos = ext_a.pos;
    ext_b.pred_pos = ext_b.pos;
    let ext = [ext_a, ext_b];
    let lambdas = [lam_a, lam_b];
    assert_vec_close(
        pbf_delta_pos(&ext, 0, &lambdas, H, RHO0),
        [expected_a, 0.0, 0.0],
        1e-12,
        "pbf_delta_pos",
    );
    assert_vec_close(
        pbf_delta_pos_nb(&ext, 1, &[0], &lambdas, H, RHO0),
        [-expected_a, 0.0, 0.0],
        1e-12,
        "pbf_delta_pos_nb",
    );
}

// ------------------------------------------------------------------
// Whole-step behaviour (F4280 / F4283 / F4284)
// ------------------------------------------------------------------

/// One `PbfSolver` step re-implemented with externally supplied neighbour
/// lists: `rebuild == false` builds them once after the prediction (Algorithm
/// 1 of the paper), `rebuild == true` rebuilds them every iteration.
fn reference_solver_step(
    solver: &PbfSolver,
    dt: f64,
    gravity: [f64; 3],
    rebuild: bool,
) -> Vec<PbfExtParticle> {
    let mut ps = solver.particles.clone();
    for p in &mut ps {
        p.predict(gravity, dt);
    }
    let mut lists = reference_lists(&pred_positions(&ps), solver.h);
    for _ in 0..solver.substeps {
        if rebuild {
            lists = reference_lists(&pred_positions(&ps), solver.h);
        }
        let n = ps.len();
        let densities: Vec<f64> = (0..n)
            .map(|i| pbf_density_nb(&ps, i, &lists[i], solver.h))
            .collect();
        let lambdas: Vec<f64> = (0..n)
            .map(|i| {
                let grad_sum = pbf_grad_sum_nb(&ps, i, &lists[i], solver.h, solver.rho0);
                pbf_lambda(densities[i], solver.rho0, grad_sum, solver.epsilon)
            })
            .collect();
        let deltas: Vec<[f64; 3]> = (0..n)
            .map(|i| pbf_delta_pos_nb(&ps, i, &lists[i], &lambdas, solver.h, solver.rho0))
            .collect();
        for (p, dp) in ps.iter_mut().zip(&deltas) {
            for (x, d) in p.pred_pos.iter_mut().zip(dp) {
                *x += d;
            }
        }
    }
    for p in &mut ps {
        p.update_velocity(dt);
    }
    ps
}

/// Dense, slightly overlapping cloud of unit-mass particles with a fixed seed.
fn dense_solver(substeps: usize, seed: u64) -> PbfSolver {
    let pts = random_cloud(220, seed, 0.2, [0.0; 3]);
    let mut solver = PbfSolver::new(RHO0, H, substeps);
    for p in ext_particles(&pts, seed + 1, (0.6, 1.4)) {
        let mut q = PbfExtParticle::new(p.pos, p.mass);
        q.vel = [0.0; 3];
        solver.add_particle(q);
    }
    solver
}

#[test]
fn single_iteration_step_matches_the_all_pairs_kernels() {
    // With one iteration the neighbour lists are valid for every evaluation, so
    // the solver must reproduce the all-pairs formulas exactly.
    let mut solver = dense_solver(1, 4);
    let (dt, g) = (2.0e-3, [0.0, -9.81, 0.0]);
    let mut expected = solver.particles.clone();
    for p in &mut expected {
        p.predict(g, dt);
    }
    let n = expected.len();
    let densities: Vec<f64> = (0..n).map(|i| pbf_density(&expected, i, H)).collect();
    let lambdas: Vec<f64> = (0..n)
        .map(|i| {
            let grad_sum = pbf_grad_sum(&expected, i, H, RHO0);
            pbf_lambda(densities[i], RHO0, grad_sum, solver.epsilon)
        })
        .collect();
    let deltas: Vec<[f64; 3]> = (0..n)
        .map(|i| pbf_delta_pos(&expected, i, &lambdas, H, RHO0))
        .collect();
    for (p, dp) in expected.iter_mut().zip(&deltas) {
        for (x, d) in p.pred_pos.iter_mut().zip(dp) {
            *x += d;
        }
    }
    for p in &mut expected {
        p.update_velocity(dt);
    }
    solver.step(dt, g);
    assert!(max_position_difference(&solver.particles, &expected) < 1e-12);
    for (p, q) in solver.particles.iter().zip(&expected) {
        assert_vec_close(p.vel, q.vel, 1e-9, "velocity");
    }
}

#[test]
fn multi_iteration_step_builds_the_neighbour_lists_once() {
    let solver = dense_solver(4, 21);
    let (dt, g) = (2.0e-3, [0.0, -9.81, 0.0]);
    let once = reference_solver_step(&solver, dt, g, false);
    let rebuilt = reference_solver_step(&solver, dt, g, true);
    // The scenario is sensitive to the policy ...
    let policy_gap = max_position_difference(&once, &rebuilt);
    assert!(policy_gap > 1.0e-10, "policy gap {policy_gap}");
    // ... and the solver follows Algorithm 1 (lists once, after the prediction).
    let mut stepped = solver.clone();
    stepped.step(dt, g);
    assert!(max_position_difference(&stepped.particles, &once) < 1e-12);
}

/// `pbf_step` re-implemented from the paper's Algorithm 1 with exhaustive
/// neighbour lists: `rebuild == false` builds them once after the prediction,
/// `rebuild == true` rebuilds them in every iteration.
fn reference_pbf_step(
    particles: &[PbfParticle],
    gravity: [f64; 3],
    dt: f64,
    h: f64,
    rho0: f64,
    epsilon: f64,
    iterations: usize,
    xsph_c: f64,
    rebuild: bool,
) -> Vec<PbfParticle> {
    let mut ps = particles.to_vec();
    for p in &mut ps {
        p.velocity = std::array::from_fn(|k| p.velocity[k] + gravity[k] * dt);
        p.predicted = std::array::from_fn(|k| p.position[k] + p.velocity[k] * dt);
    }
    let predicted = |ps: &[PbfParticle]| ps.iter().map(|p| p.predicted).collect::<Vec<_>>();
    let mut lists = reference_lists(&predicted(&ps), h);
    let once = lists.clone();
    for _ in 0..iterations {
        if rebuild {
            lists = reference_lists(&predicted(&ps), h);
        }
        let lambdas: Vec<f64> = (0..ps.len())
            .map(|i| lambda_compute(&ps, &lists[i], i, h, rho0, epsilon))
            .collect();
        for (p, lambda) in ps.iter_mut().zip(&lambdas) {
            p.lambda = *lambda;
        }
        let deltas: Vec<[f64; 3]> = (0..ps.len())
            .map(|i| delta_position(&ps, &lists[i], i, h, rho0))
            .collect();
        for (p, dp) in ps.iter_mut().zip(&deltas) {
            p.predicted = std::array::from_fn(|k| p.predicted[k] + dp[k]);
        }
    }
    let lists = if rebuild { once } else { lists };
    for p in &mut ps {
        p.velocity = std::array::from_fn(|k| (p.predicted[k] - p.position[k]) / dt);
    }
    let corrections: Vec<[f64; 3]> = (0..ps.len())
        .map(|i| xsph_viscosity(&ps, &lists[i], i, h, xsph_c))
        .collect();
    for (p, dv) in ps.iter_mut().zip(&corrections) {
        p.velocity = std::array::from_fn(|k| p.velocity[k] + dv[k]);
        p.position = p.predicted;
    }
    ps
}

#[test]
fn pbf_step_builds_the_neighbour_lists_once_per_step() {
    let pts = random_cloud(200, 71, 0.2, [0.0; 3]);
    let mut rng = Lcg::new(72);
    let start: Vec<PbfParticle> = pts
        .iter()
        .map(|&pos| {
            let mut p = PbfParticle::new(pos, rng.range(0.6, 1.4));
            p.velocity = [
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
                rng.range(-0.5, 0.5),
            ];
            p
        })
        .collect();
    let (g, dt, rho0, eps, iters, xsph) = ([0.0, -9.81, 0.0], 2.0e-3, 1000.0, 600.0, 4, 0.05);
    let once = reference_pbf_step(&start, g, dt, H, rho0, eps, iters, xsph, false);
    let rebuilt = reference_pbf_step(&start, g, dt, H, rho0, eps, iters, xsph, true);
    let spread = |a: &[PbfParticle], b: &[PbfParticle]| {
        a.iter()
            .zip(b)
            .map(|(p, q)| {
                (0..3)
                    .map(|k| (p.position[k] - q.position[k]).abs())
                    .fold(0.0, f64::max)
            })
            .fold(0.0, f64::max)
    };
    // The scenario distinguishes the two policies ...
    assert!(
        spread(&once, &rebuilt) > 1.0e-10,
        "policy gap {}",
        spread(&once, &rebuilt)
    );
    // ... and `pbf_step` follows Algorithm 1 (lists once) in positions and velocities.
    let mut stepped = start.clone();
    pbf_step(&mut stepped, g, dt, H, rho0, eps, iters, xsph);
    assert!(
        spread(&stepped, &once) < 1e-12,
        "gap {}",
        spread(&stepped, &once)
    );
    for (p, q) in stepped.iter().zip(&once) {
        assert_vec_close(p.velocity, q.velocity, 1e-9, "velocity");
    }
}

#[test]
fn solver_and_pbf_step_conserve_linear_momentum_for_equal_masses() {
    // Position corrections and XSPH are antisymmetric over particle pairs, so
    // with zero gravity the total momentum cannot change.
    let pts = random_cloud(180, 61, 0.2, [0.0; 3]);
    let mut rng = Lcg::new(62);
    let velocities: Vec<[f64; 3]> = pts
        .iter()
        .map(|_| {
            [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ]
        })
        .collect();

    let mut solver = PbfSolver::new(RHO0, H, 3);
    for (&pos, &vel) in pts.iter().zip(&velocities) {
        let mut p = PbfExtParticle::new(pos, 1.0);
        p.vel = vel;
        solver.add_particle(p);
    }
    let before = ext_momentum(&solver.particles);
    let before_com = solver.particles.iter().fold([0.0; 3], |acc, p| {
        [acc[0] + p.pos[0], acc[1] + p.pos[1], acc[2] + p.pos[2]]
    });
    for _ in 0..4 {
        solver.step(1.0e-3, [0.0; 3]);
    }
    let after = ext_momentum(&solver.particles);
    assert_vec_close(after, before, 1e-9, "solver momentum");
    // The velocity-integrated centre of mass moves by exactly t · p / M.
    let after_com = solver.particles.iter().fold([0.0; 3], |acc, p| {
        [acc[0] + p.pos[0], acc[1] + p.pos[1], acc[2] + p.pos[2]]
    });
    for k in 0..3 {
        assert_close(
            after_com[k] - before_com[k],
            4.0e-3 * before[k],
            1e-9,
            "centre-of-mass drift",
        );
    }

    let mut classic: Vec<PbfParticle> = pts
        .iter()
        .zip(&velocities)
        .map(|(&pos, &vel)| {
            let mut p = PbfParticle::new(pos, 1.0);
            p.velocity = vel;
            p
        })
        .collect();
    let momentum = |ps: &[PbfParticle]| {
        ps.iter().fold([0.0; 3], |acc, p| {
            [
                acc[0] + p.mass * p.velocity[0],
                acc[1] + p.mass * p.velocity[1],
                acc[2] + p.mass * p.velocity[2],
            ]
        })
    };
    let before_classic = momentum(&classic);
    for _ in 0..3 {
        pbf_step(&mut classic, [0.0; 3], 1.0e-3, H, RHO0, 1.0e-4, 3, 0.1);
    }
    assert_vec_close(
        momentum(&classic),
        before_classic,
        1e-9,
        "pbf_step momentum",
    );
    assert!(
        classic
            .iter()
            .all(|p| p.position.iter().all(|c| c.is_finite()))
    );
}

#[test]
fn xsph_acts_on_the_velocity_produced_by_the_constraint_solve() {
    // Two overlapping particles at rest are pushed apart by the density
    // constraint (A to −x, B to +x). XSPH, applied to the *corrected*
    // velocities, then shrinks their relative velocity by the factor
    // 1 − 2 c W(r_final); the pre-solve velocities (both zero) would give no
    // change at all.
    let step_with = |c: f64| {
        let mut ps = vec![
            PbfParticle::new([-0.01, 0.0, 0.0], 1.0),
            PbfParticle::new([0.01, 0.0, 0.0], 1.0),
        ];
        pbf_step(&mut ps, [0.0; 3], 1.0e-3, H, 100.0, 1.0e-4, 1, c);
        ps
    };
    let plain = step_with(0.0);
    let relative_plain = plain[1].velocity[0] - plain[0].velocity[0];
    assert!(
        relative_plain > 1.0,
        "the pair must be pushed apart: {relative_plain}"
    );
    assert!(plain[0].velocity[0] < 0.0 && plain[1].velocity[0] > 0.0);

    let c = 0.2;
    let smoothed = step_with(c);
    let relative_smoothed = smoothed[1].velocity[0] - smoothed[0].velocity[0];
    let r_final = plain[1].position[0] - plain[0].position[0];
    let w = poly6_kernel(r_final, H);
    assert!(w > 0.0);
    assert_close(
        relative_smoothed,
        relative_plain * (1.0 - 2.0 * c * w),
        1e-9,
        "XSPH relative velocity",
    );
    // Positions are not affected by the viscosity.
    assert_vec_close(smoothed[0].position, plain[0].position, 1e-14, "position A");
}

#[test]
fn compressed_block_relaxes_towards_the_rest_density() {
    // A 9³ lattice. ρ0 is calibrated to the density an *interior* particle
    // has at spacing s0 (a lone particle already has m·W(0) ≈ 1567 kg/m³ for
    // m = 1, so the module constants alone would leave everything overdense).
    let s0 = H / 2.5;
    let w_sum: f64 = {
        let mut sum = 0.0;
        for a in -3i64..=3 {
            for b in -3i64..=3 {
                for c in -3i64..=3 {
                    let r = s0 * ((a * a + b * b + c * c) as f64).sqrt();
                    sum += poly6_kernel(r, H);
                }
            }
        }
        sum
    };
    let mass = RHO0 / w_sum;

    let n = 9usize;
    let squeeze = 0.8;
    let mut solver = PbfSolver::new(RHO0, H, 4);
    for p in lattice(n, s0 * squeeze) {
        let mut q = PbfExtParticle::new(p, mass);
        q.vel = [0.0; 3];
        solver.add_particle(q);
    }
    let interior: Vec<usize> = (0..n * n * n)
        .filter(|&idx| {
            let (ix, iy, iz) = (idx / (n * n), (idx / n) % n, idx % n);
            [ix, iy, iz].iter().all(|&c| (3..=5).contains(&c))
        })
        .collect();
    assert_eq!(interior.len(), 27);
    let interior_error = |solver: &PbfSolver| {
        let csr = NeighborCsr::from_positions(&pred_positions(&solver.particles), solver.h);
        interior
            .iter()
            .map(|&i| {
                let rho = pbf_density_nb(&solver.particles, i, csr.neighbors(i), solver.h);
                (rho / solver.rho0 - 1.0).abs()
            })
            .sum::<f64>()
            / interior.len() as f64
    };

    let initial_interior = interior_error(&solver);
    let initial_total = solver.total_density_error();
    assert!(
        initial_interior > 0.8,
        "compression too weak: {initial_interior}"
    );
    let mut fastest: f64 = 0.0;
    for _ in 0..40 {
        solver.step(1.0e-3, [0.0; 3]);
        for p in &solver.particles {
            fastest = fastest
                .max((p.vel[0] * p.vel[0] + p.vel[1] * p.vel[1] + p.vel[2] * p.vel[2]).sqrt());
        }
    }
    let final_interior = interior_error(&solver);
    let final_total = solver.total_density_error();
    assert!(
        final_interior < 0.15 * initial_interior,
        "interior density error {initial_interior} -> {final_interior}"
    );
    assert!(
        final_total < 0.3 * initial_total,
        "total density error {initial_total} -> {final_total}"
    );
    // Stable relaxation: the over-pressure is released at O(10 m/s); an
    // unrelaxed iteration (ε ≪ Σ|∇C|²) ejects particles at > 10⁴ m/s.
    assert!(fastest < 500.0, "peak speed {fastest} m/s");
    assert!(
        solver
            .particles
            .iter()
            .all(|p| p.pos.iter().chain(p.vel.iter()).all(|c| c.is_finite()))
    );
    // No external forces and a symmetric start: the block cannot drift.
    let com = solver.particles.iter().fold([0.0; 3], |acc, p| {
        [acc[0] + p.pos[0], acc[1] + p.pos[1], acc[2] + p.pos[2]]
    });
    for c in com {
        assert!(c.abs() < 1e-9, "centre of mass drifted: {com:?}");
    }
}

#[test]
fn non_positive_time_steps_leave_the_particles_untouched() {
    let pts = random_cloud(20, 2, 0.1, [0.0; 3]);
    for dt in [0.0, -1.0e-3, f64::NAN] {
        let mut solver = PbfSolver::new(RHO0, H, 2);
        for p in ext_particles(&pts, 3, (1.0, 1.0)) {
            solver.add_particle(p);
        }
        let before = solver.particles.clone();
        solver.step(dt, [0.0, -9.81, 0.0]);
        assert_eq!(max_position_difference(&solver.particles, &before), 0.0);
        assert!(
            solver
                .particles
                .iter()
                .zip(&before)
                .all(|(p, q)| p.vel == q.vel)
        );

        let mut classic: Vec<PbfParticle> = pts.iter().map(|&p| PbfParticle::new(p, 1.0)).collect();
        pbf_step(&mut classic, [0.0, -9.81, 0.0], dt, H, RHO0, 1.0e-4, 2, 0.1);
        assert!(
            classic
                .iter()
                .zip(&pts)
                .all(|(p, q)| p.position == *q && p.velocity == [0.0; 3])
        );
    }
}
