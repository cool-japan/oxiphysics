//! Tests of the sparse topology filters against the former dense all-pairs loops.
//!
//! Every dense reference below is a verbatim copy of the loop that the sparse
//! filter replaced, so the comparisons pin the sparse results to the old output
//! (to rounding) while the scaling tests pin the new `O(n · k)` behaviour.

use super::*;
use crate::topology_opt::{
    ProjectionFilter, TopologyGrid, filter_sensitivities_gaussian, filter_sensitivities_mesh,
    heaviside_derivative,
};
use crate::{topology_optimization as classic, topology_optimization_ext as ext};

// ── deterministic data ───────────────────────────────────────────────────────

/// Small LCG so the tests need no external RNG and stay reproducible.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// Uniform sample in `[0, 1)`.
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn next_in(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }
}

fn random_cloud(n: usize, seed: u64, extent: f64) -> Vec<(f64, f64)> {
    let mut rng = Lcg::new(seed);
    (0..n)
        .map(|_| (rng.next_in(-extent, extent), rng.next_in(-extent, extent)))
        .collect()
}

fn random_field(n: usize, seed: u64, lo: f64, hi: f64) -> Vec<f64> {
    let mut rng = Lcg::new(seed);
    (0..n).map(|_| rng.next_in(lo, hi)).collect()
}

// ── dense references (the replaced all-pairs loops) ──────────────────────────

fn dense_hat(c: &[(f64, f64)], r: f64) -> Vec<Vec<f64>> {
    let n = c.len();
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let d = ((c[i].0 - c[j].0).powi(2) + (c[i].1 - c[j].1).powi(2)).sqrt();
                    (r - d).max(0.0)
                })
                .collect()
        })
        .collect()
}

fn dense_gaussian(c: &[(f64, f64)], r: f64) -> Vec<Vec<f64>> {
    let n = c.len();
    let sigma = r / 3.0;
    (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let d2 = (c[i].0 - c[j].0).powi(2) + (c[i].1 - c[j].1).powi(2);
                    if d2 < r * r {
                        (-d2 / (2.0 * sigma * sigma)).exp()
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect()
}

/// `Σⱼ wᵢⱼ vⱼ / Σⱼ wᵢⱼ` with the `den > 1e-30` fallback.
fn dense_average(w: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    (0..v.len())
        .map(|i| {
            let mut num = 0.0;
            let mut den = 0.0;
            for j in 0..v.len() {
                num += w[i][j] * v[j];
                den += w[i][j];
            }
            if den > 1e-30 { num / den } else { v[i] }
        })
        .collect()
}

/// `Σⱼ wᵢⱼ ρⱼ vⱼ / Σⱼ wᵢⱼ ρⱼ` with the `|den| > 1e-30` fallback.
fn dense_weighted_average(w: &[Vec<f64>], v: &[f64], rho: &[f64]) -> Vec<f64> {
    (0..v.len())
        .map(|i| {
            let mut num = 0.0;
            let mut den = 0.0;
            for j in 0..v.len() {
                num += w[i][j] * rho[j] * v[j];
                den += w[i][j] * rho[j];
            }
            if den.abs() > 1e-30 { num / den } else { v[i] }
        })
        .collect()
}

fn assert_close(got: &[f64], want: &[f64], tol: f64, what: &str) {
    assert_eq!(got.len(), want.len(), "{what}: length");
    for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
        assert!(
            (g - w).abs() <= tol * (1.0 + w.abs()),
            "{what}[{i}]: sparse {g} vs dense {w}"
        );
    }
}

// ── cell-list correctness ────────────────────────────────────────────────────

#[test]
fn sparse_weights_equal_dense_weights_on_random_point_clouds() {
    for (seed, n) in [(1_u64, 150_usize), (7, 220)] {
        let cloud = random_cloud(n, seed, 3.0);
        // radii below, around and far above the typical spacing
        for r in [0.05, 0.4, 1.0, 2.5, 20.0] {
            for kernel in [FilterKernel::Hat, FilterKernel::Gaussian] {
                let dense = match kernel {
                    FilterKernel::Hat => dense_hat(&cloud, r),
                    FilterKernel::Gaussian => dense_gaussian(&cloud, r),
                };
                let filter = SparseFilter::new(&cloud, r, kernel);
                let mut nonzero = 0;
                for i in 0..n {
                    let mut row_sum = 0.0;
                    for j in 0..n {
                        let got = filter.h.get(i, j);
                        assert!(
                            (got - dense[i][j]).abs() <= 1e-15 * (1.0 + dense[i][j]),
                            "{kernel:?} r={r} ({i},{j}): sparse {got} vs dense {}",
                            dense[i][j]
                        );
                        if dense[i][j] > 0.0 {
                            nonzero += 1;
                        }
                        row_sum += dense[i][j];
                    }
                    assert!(
                        (filter.row_sum[i] - row_sum).abs() <= 1e-13 * (1.0 + row_sum),
                        "{kernel:?} r={r}: row sum {i}"
                    );
                }
                assert_eq!(filter.nnz(), nonzero, "{kernel:?} r={r}: stored weights");
            }
        }
    }
}

#[test]
fn columns_are_ascending_within_each_row() {
    let cloud = random_cloud(300, 11, 2.0);
    let filter = SparseFilter::hat(&cloud, 0.7);
    for i in 0..cloud.len() {
        let cols = &filter.h.col_indices[filter.h.row_ptr[i]..filter.h.row_ptr[i + 1]];
        assert!(cols.windows(2).all(|w| w[0] < w[1]), "row {i}: {cols:?}");
        assert!(cols.contains(&i), "row {i} must contain its own weight");
    }
}

#[test]
fn kernels_are_truncated_at_the_radius() {
    let line: Vec<(f64, f64)> = (0..4).map(|i| (i as f64, 0.0)).collect();
    let hat = SparseFilter::hat(&line, 2.0);
    assert_eq!(hat.h.get(0, 0), 2.0);
    assert_eq!(hat.h.get(0, 1), 1.0);
    assert_eq!(hat.h.get(0, 2), 0.0); // d = r ⇒ zero weight
    assert_eq!(hat.row_sum[0], 3.0);
    assert_eq!(hat.row_sum[1], 4.0);
    let gauss = SparseFilter::new(&line, 2.0, FilterKernel::Gaussian);
    let sigma = 2.0 / 3.0;
    assert_eq!(gauss.h.get(0, 0), 1.0);
    assert!((gauss.h.get(0, 1) - (-1.0 / (2.0 * sigma * sigma)).exp()).abs() < 1e-15);
    assert_eq!(gauss.h.get(0, 2), 0.0); // d = r is excluded (strict cut-off)
}

// ── parity of every public filter with its dense predecessor ─────────────────

#[test]
fn filters_match_the_dense_reference_on_a_twenty_by_ten_grid() {
    let (nx, ny, r_min) = (20, 10, 2.5);
    let centroids = ext::regular_grid_centroids(nx, ny, 1.0);
    let n = centroids.len();
    let mut densities = random_field(n, 3, 0.05, 1.0);
    // a void pocket exercises the zero-denominator fallback of the weighted filter
    for j in 0..ny.min(6) {
        for i in 0..6 {
            densities[j * nx + i] = 0.0;
        }
    }
    let sens = random_field(n, 4, -5.0, -0.1);
    let w = dense_hat(&centroids, r_min);

    let want_avg_rho = dense_average(&w, &densities);
    let want_avg_sens = dense_average(&w, &sens);
    let want_weighted = dense_weighted_average(&w, &sens, &densities);

    assert_close(
        &ext::density_filter(&densities, &centroids, r_min),
        &want_avg_rho,
        1e-12,
        "ext::density_filter",
    );
    assert_close(
        &ext::filter_sensitivities_hat(&sens, &centroids, r_min),
        &want_avg_sens,
        1e-12,
        "ext::filter_sensitivities_hat",
    );
    assert_close(
        &ext::filter_sensitivities_density_weighted(&sens, &densities, &centroids, r_min),
        &want_weighted,
        1e-12,
        "ext::filter_sensitivities_density_weighted",
    );
    assert_close(
        &classic::filter_sensitivities_hat(&sens, &centroids, r_min),
        &want_avg_sens,
        1e-12,
        "classic::filter_sensitivities_hat",
    );
    assert_close(
        &classic::filter_sensitivities_density_weighted(&sens, &densities, &centroids, r_min),
        &want_weighted,
        1e-12,
        "classic::filter_sensitivities_density_weighted",
    );
    assert_close(
        &filter_sensitivities_mesh(&sens, &centroids, r_min),
        &want_avg_sens,
        1e-12,
        "filter_sensitivities_mesh",
    );

    // density filter + Heaviside projection
    let want_projected = classic::apply_heaviside_projection(&want_avg_rho, 4.0, 0.5);
    assert_close(
        &classic::density_filter_and_project(&densities, &centroids, r_min, 4.0, 0.5),
        &want_projected,
        1e-12,
        "density_filter_and_project",
    );

    // minimum length scale filter: clamped average, fallback unclamped
    let wild = random_field(n, 5, -0.5, 1.5);
    let want_min_scale: Vec<f64> = {
        let avg = dense_average(&w, &wild);
        avg.iter().map(|v| v.clamp(0.0, 1.0)).collect()
    };
    assert_close(
        &classic::minimum_length_scale_filter(&wild, &centroids, r_min),
        &want_min_scale,
        1e-12,
        "minimum_length_scale_filter",
    );
}

#[test]
fn density_filter_matrix_equals_the_normalised_dense_matrix() {
    let centroids = ext::regular_grid_centroids(20, 10, 1.0);
    let n = centroids.len();
    let r_min = 2.5;
    let h = ext::density_filter_matrix(&centroids, r_min);
    assert_eq!((h.nrows, h.ncols), (n, n));
    let w = dense_hat(&centroids, r_min);
    for i in 0..n {
        let row_sum: f64 = w[i].iter().sum();
        for j in 0..n {
            let want = w[i][j] / row_sum;
            assert!(
                (h.get(i, j) - want).abs() <= 1e-14,
                "H[{i}][{j}]: {} vs {want}",
                h.get(i, j)
            );
        }
        let stored: f64 = h.values[h.row_ptr[i]..h.row_ptr[i + 1]].iter().sum();
        assert!((stored - 1.0).abs() < 1e-12, "row {i} sums to {stored}");
    }
    // sparse: far fewer stored entries than the n² of the former dense matrix
    assert!(h.nnz() < n * 30, "nnz = {} for n = {n}", h.nnz());
    // and it applies like the averaging filter
    let field = random_field(n, 6, 0.0, 1.0);
    let via_matrix = h.mul_vec(&field);
    assert_close(
        &via_matrix,
        &ext::density_filter(&field, &centroids, r_min),
        1e-13,
        "normalised matrix · field",
    );
}

#[test]
fn grid_filters_match_the_dense_reference() {
    let mut grid = TopologyGrid::new_uniform(20, 10, 0.5, 1.0, 1.0, 1.0);
    let n = grid.elements.len();
    let rho = random_field(n, 8, 0.0, 1.0);
    let sens = random_field(n, 9, -4.0, 0.5);
    for (e, (&r, &s)) in grid.elements.iter_mut().zip(rho.iter().zip(sens.iter())) {
        e.rho = r;
        e.sensitivity = s;
    }
    let centroids = grid.centroids();
    for radius in [1.2, 2.2, 3.7] {
        let w = dense_hat(&centroids, radius);
        assert_close(
            &grid.filter_sensitivities(radius),
            &dense_weighted_average(&w, &sens, &rho),
            1e-12,
            "TopologyGrid::filter_sensitivities",
        );
        assert_close(
            &grid.filter_densities(radius),
            &dense_average(&w, &rho),
            1e-12,
            "TopologyGrid::filter_densities",
        );
        // projection filter chain rule
        let projection = ProjectionFilter::new(radius, 6.0, 0.5);
        let (rho_bar, _rho_phys) = projection.apply(&grid);
        let dc_phys = random_field(n, 10, -1.0, 1.0);
        let dc_times_dh: Vec<f64> = dc_phys
            .iter()
            .zip(rho_bar.iter())
            .map(|(&dc, &rb)| dc * heaviside_derivative(rb, 6.0, 0.5))
            .collect();
        assert_close(
            &projection.chain_sensitivity(&grid, &dc_phys, &rho_bar),
            &dense_average(&w, &dc_times_dh),
            1e-12,
            "ProjectionFilter::chain_sensitivity",
        );
        // Gaussian sensitivity filter
        let g = dense_gaussian(&centroids, radius);
        assert_close(
            &filter_sensitivities_gaussian(&grid, &sens, radius),
            &dense_weighted_average(&g, &sens, &rho),
            1e-12,
            "filter_sensitivities_gaussian",
        );
    }
}

#[test]
fn irregular_meshes_with_negative_coordinates_match_the_dense_reference() {
    let cloud = random_cloud(260, 21, 4.0);
    let rho = random_field(260, 22, 0.0, 1.0);
    let sens = random_field(260, 23, -3.0, 3.0);
    for r in [0.3, 0.9, 2.0, 12.0] {
        let w = dense_hat(&cloud, r);
        let filter = SparseFilter::hat(&cloud, r);
        assert_close(
            &filter.average(&sens),
            &dense_average(&w, &sens),
            1e-12,
            "average",
        );
        assert_close(
            &filter.density_weighted_average(&sens, &rho),
            &dense_weighted_average(&w, &sens, &rho),
            1e-12,
            "density_weighted_average",
        );
    }
}

// ── API behaviour ────────────────────────────────────────────────────────────

#[test]
fn normalised_matrix_applies_like_average_and_rows_sum_to_one() {
    let cloud = random_cloud(120, 31, 2.0);
    let field = random_field(120, 32, -1.0, 1.0);
    let filter = SparseFilter::hat(&cloud, 0.8);
    let m = filter.normalised();
    assert_close(&m.mul_vec(&field), &filter.average(&field), 1e-13, "H_n v");
    for i in 0..120 {
        let s: f64 = m.values[m.row_ptr[i]..m.row_ptr[i + 1]].iter().sum();
        assert!((s - 1.0).abs() < 1e-12, "row {i}: {s}");
    }
    // the raw matrix is left untouched
    assert!(filter.h.values.iter().any(|&v| v > 1.0 - 1e-12 && v <= 0.8));
}

#[test]
fn empty_and_degenerate_radii_give_identity_filters() {
    let empty = SparseFilter::hat(&[], 1.0);
    assert!(empty.is_empty());
    assert_eq!(empty.len(), 0);
    assert!(empty.apply(&[]).is_empty());
    assert_eq!(empty.normalised().nnz(), 0);

    let cloud = random_cloud(30, 41, 1.0);
    let v = random_field(30, 42, 0.0, 1.0);
    let rho = random_field(30, 43, 0.1, 1.0);
    for r in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let filter = SparseFilter::hat(&cloud, r);
        assert_eq!(filter.len(), 30);
        assert_eq!(filter.nnz(), 0, "radius {r}");
        assert!(filter.row_sum.iter().all(|&s| s == 0.0));
        assert_eq!(filter.average(&v), v);
        assert_eq!(filter.density_weighted_average(&v, &rho), v);
        assert!(filter.apply(&v).iter().all(|&y| y == 0.0));
    }
}

#[test]
fn non_finite_and_extreme_centroids_do_not_poison_the_filter() {
    let mut cloud = random_cloud(40, 51, 1.0);
    cloud[5] = (f64::NAN, 0.3);
    cloud[17] = (0.2, f64::INFINITY);
    cloud[23] = (1.0e300, -1.0e300);
    let v = random_field(40, 52, 0.0, 1.0);
    let filter = SparseFilter::hat(&cloud, 0.6);
    let w = dense_hat(&cloud, 0.6);
    assert_close(
        &filter.average(&v),
        &dense_average(&w, &v),
        1e-12,
        "average",
    );
    for bad in [5, 17] {
        assert_eq!(filter.row_sum[bad], 0.0);
        assert_eq!(filter.average(&v)[bad], v[bad]);
    }
}

#[test]
fn coincident_points_share_one_cell() {
    let cloud = vec![(0.25, -0.5); 9];
    let v: Vec<f64> = (0..9).map(|i| i as f64).collect();
    let filter = SparseFilter::hat(&cloud, 1.5);
    assert!(
        filter
            .row_sum
            .iter()
            .all(|&s| (s - 9.0 * 1.5).abs() < 1e-12)
    );
    let mean = 4.0;
    assert!(filter.average(&v).iter().all(|&a| (a - mean).abs() < 1e-12));
}

#[test]
#[should_panic(expected = "element count")]
fn applying_a_field_of_the_wrong_length_panics() {
    let cloud = random_cloud(10, 61, 1.0);
    let filter = SparseFilter::hat(&cloud, 0.5);
    let _ = filter.apply(&[1.0; 9]);
}

// ── scaling ──────────────────────────────────────────────────────────────────

#[test]
fn two_hundred_by_one_hundred_mesh_needs_only_linear_memory() {
    // The dense filter matrix of this mesh needs 20 000² · 8 B = 3.2 GB.
    let (nx, ny, r) = (200_usize, 100_usize, 3.0_f64);
    let centroids = ext::regular_grid_centroids(nx, ny, 1.0);
    let n = centroids.len();
    let filter = SparseFilter::hat(&centroids, r);
    assert_eq!(filter.len(), n);
    // at most one weight per lattice point inside the radius: π r² ≈ 28
    assert!(filter.nnz() <= 25 * n, "nnz = {}", filter.nnz());
    assert!(filter.nnz() >= 15 * n, "nnz = {}", filter.nnz());

    // an interior element sees the full stencil Σ (r − d) over the lattice offsets
    let mut stencil_sum = 0.0;
    for di in -3_i32..=3 {
        for dj in -3_i32..=3 {
            let d = f64::from(di * di + dj * dj).sqrt();
            stencil_sum += (r - d).max(0.0);
        }
    }
    let centre = 50 * nx + 100;
    assert!(
        (filter.row_sum[centre] - stencil_sum).abs() < 1e-9,
        "{} vs {stencil_sum}",
        filter.row_sum[centre]
    );

    // uniform fields are preserved, and linear ramps away from the boundary
    let uniform = vec![0.7; n];
    assert!(
        filter
            .average(&uniform)
            .iter()
            .all(|&v| (v - 0.7).abs() < 1e-12)
    );
    let ramp: Vec<f64> = centroids.iter().map(|c| c.0).collect();
    let smoothed = filter.average(&ramp);
    for j in 5..ny - 5 {
        for i in 5..nx - 5 {
            let idx = j * nx + i;
            assert!(
                (smoothed[idx] - ramp[idx]).abs() < 1e-9,
                "interior element ({i},{j}): {} vs {}",
                smoothed[idx],
                ramp[idx]
            );
        }
    }

    // the (row-parallel) SpMV equals a plain sequential CSR product bit for bit
    let field = random_field(n, 71, -1.0, 1.0);
    let parallel = filter.apply(&field);
    for (i, &p) in parallel.iter().enumerate() {
        let mut sum = 0.0;
        for idx in filter.h.row_ptr[i]..filter.h.row_ptr[i + 1] {
            sum += filter.h.values[idx] * field[filter.h.col_indices[idx]];
        }
        assert_eq!(p, sum, "row {i}");
    }
}

#[test]
fn mock_optimisation_driver_runs_with_the_reused_filter() {
    let cfg = ext::TopologyConfig {
        nx: 12,
        ny: 6,
        max_iter: 8,
        ..Default::default()
    };
    let result = ext::run_mock_topology_optimisation(&cfg);
    assert_eq!(result.densities.len(), 72);
    assert!(
        result
            .densities
            .iter()
            .all(|&r| (cfg.rho_min..=1.0).contains(&r))
    );
    assert!((result.volume_fraction - cfg.vf).abs() < 0.05);
}
