//! Regression tests for the Fast Marching Method.
//!
//! The update must solve the Eikonal equation with one, two or three upwind
//! neighbours, use accepted cells only and keep the sign of the seeds. An earlier
//! version always used the one-neighbour update (Manhattan distances) and let
//! negative seeds grow towards positive values.

use super::FastMarchingMethod;

const SQRT_2: f64 = std::f64::consts::SQRT_2;

fn euclid(a: [usize; 3], b: [usize; 3]) -> f64 {
    let d = |i: usize| a[i] as f64 - b[i] as f64;
    (d(0) * d(0) + d(1) * d(1) + d(2) * d(2)).sqrt()
}

/// FMM on an `n`-cube seeded with the single cell `centre` at distance zero.
fn point_source(n: usize, dx: f64) -> FastMarchingMethod {
    let c = n / 2;
    let mut fmm = FastMarchingMethod::new(n, n, n, dx);
    let seed = fmm.flat(c, c, c);
    fmm.set_known(&[(seed, 0.0)]);
    fmm.run();
    fmm
}

#[test]
fn axis_neighbours_get_exactly_the_one_term_update() {
    for dx in [1.0, 0.5, 0.125] {
        let fmm = point_source(9, dx);
        for k in 1..=4 {
            let expected = k as f64 * dx;
            for (ix, iy, iz) in [(4 + k, 4, 4), (4 - k, 4, 4), (4, 4 + k, 4), (4, 4, 4 - k)] {
                let d = fmm.distance_at(ix, iy, iz);
                assert!(
                    (d - expected).abs() < 1e-12,
                    "dx {dx}, cell ({ix},{iy},{iz}): {d} vs {expected}"
                );
            }
        }
    }
}

#[test]
fn face_diagonal_neighbour_uses_the_two_term_update() {
    // (t - 1)^2 + (t - 1)^2 = dx^2 with dx = 1 gives t = 1 + sqrt(2) / 2, not the
    // Manhattan value 2.
    for dx in [1.0, 0.3] {
        let fmm = point_source(9, dx);
        let expected = dx * (1.0 + SQRT_2 / 2.0);
        for (ix, iy, iz) in [(5, 5, 4), (3, 5, 4), (5, 4, 3), (4, 3, 3)] {
            let d = fmm.distance_at(ix, iy, iz);
            assert!(
                (d - expected).abs() < 1e-12,
                "dx {dx}, cell ({ix},{iy},{iz}): {d} vs {expected}"
            );
        }
    }
}

#[test]
fn space_diagonal_neighbour_uses_the_three_term_update() {
    // Three neighbours at 1 + sqrt(2) / 2 give t = 1 + sqrt(2) / 2 + sqrt(3) / 3
    // (not the Manhattan value 3).
    let fmm = point_source(9, 1.0);
    let expected = 1.0 + SQRT_2 / 2.0 + 3.0_f64.sqrt() / 3.0;
    for (ix, iy, iz) in [(5, 5, 5), (3, 3, 3), (5, 3, 5)] {
        let d = fmm.distance_at(ix, iy, iz);
        assert!(
            (d - expected).abs() < 1e-12,
            "cell ({ix},{iy},{iz}): {d} vs {expected}"
        );
    }
}

#[test]
fn point_source_distances_never_underestimate_and_are_not_manhattan() {
    let fmm = point_source(13, 1.0);
    let centre = [6, 6, 6];
    for iz in 0..13 {
        for iy in 0..13 {
            for ix in 0..13 {
                let reference = euclid([ix, iy, iz], centre);
                if reference == 0.0 {
                    continue;
                }
                let d = fmm.distance_at(ix, iy, iz);
                assert!(
                    d >= reference - 1e-9,
                    "({ix},{iy},{iz}): {d} underestimates {reference}"
                );
                // The one-sided first-order scheme overestimates by at most a third
                // next to the source and by about 10 % at the far corner; the
                // Manhattan distance would be 73 % too large on the space diagonal.
                assert!(
                    d <= 1.32 * reference,
                    "({ix},{iy},{iz}): {d} overestimates {reference}"
                );
            }
        }
    }
    let corner = fmm.distance_at(12, 12, 12);
    assert!(
        corner < 1.12 * 6.0 * 3.0_f64.sqrt(),
        "corner distance {corner}"
    );
}

#[test]
fn exact_ball_seeds_propagate_within_ten_percent() {
    let n = 13;
    let c = 6;
    let mut fmm = FastMarchingMethod::new(n, n, n, 1.0);
    let mut seeds = Vec::new();
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let r = euclid([ix, iy, iz], [c, c, c]);
                if r <= 2.5 {
                    seeds.push((fmm.flat(ix, iy, iz), r));
                }
            }
        }
    }
    fmm.set_known(&seeds);
    fmm.run();
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let r = euclid([ix, iy, iz], [c, c, c]);
                if r > 2.5 {
                    let d = fmm.distance_at(ix, iy, iz);
                    assert!(
                        d >= r - 1e-9 && d <= 1.10 * r,
                        "({ix},{iy},{iz}): {d} vs exact {r}"
                    );
                }
            }
        }
    }
    let k = 6;
    let diagonal = fmm.distance_at(c + k, c + k, c + k);
    let exact = k as f64 * 3.0_f64.sqrt();
    assert!(
        (diagonal - exact).abs() < 0.10 * exact,
        "diagonal {diagonal} vs {exact}"
    );
}

#[test]
fn planar_grid_uses_the_two_term_update() {
    let mut fmm = FastMarchingMethod::new(9, 9, 1, 1.0);
    let seed = fmm.flat(4, 4, 0);
    fmm.set_known(&[(seed, 0.0)]);
    fmm.run();
    assert!((fmm.distance_at(5, 5, 0) - (1.0 + SQRT_2 / 2.0)).abs() < 1e-12);
    assert!((fmm.distance_at(8, 4, 0) - 4.0).abs() < 1e-12);
    // Manhattan distance to (8, 8) would be 8; the Euclidean one is 4 sqrt(2).
    let far = fmm.distance_at(8, 8, 0);
    assert!(far < 1.2 * 4.0 * SQRT_2, "far corner {far}");
}

#[test]
fn every_cell_of_an_anisotropic_grid_is_reached() {
    let mut fmm = FastMarchingMethod::new(7, 4, 3, 0.5);
    let seed = fmm.flat(0, 0, 0);
    fmm.set_known(&[(seed, 0.0)]);
    fmm.run();
    for iz in 0..3 {
        for iy in 0..4 {
            for ix in 0..7 {
                let d = fmm.distance_at(ix, iy, iz);
                let reference = 0.5 * euclid([ix, iy, iz], [0, 0, 0]);
                assert!(d < f64::MAX, "({ix},{iy},{iz}) was not reached");
                assert!(d >= reference - 1e-9 && d <= 1.32 * reference + 1e-9);
            }
        }
    }
}

/// Sphere of radius 1.5 on a 21-cube with spacing 0.25, seeded with the exact
/// signed distance in a band around the surface.
fn sphere_band(half_width: f64) -> (FastMarchingMethod, Vec<f64>, Vec<bool>) {
    let n = 21;
    let h = 0.25;
    let radius = 1.5;
    let mut fmm = FastMarchingMethod::new(n, n, n, h);
    let mut exact = vec![0.0; n * n * n];
    let mut seeded = vec![false; n * n * n];
    let mut seeds = Vec::new();
    for iz in 0..n {
        for iy in 0..n {
            for ix in 0..n {
                let p = |i: usize| -2.5 + i as f64 * h;
                let signed = (p(ix).powi(2) + p(iy).powi(2) + p(iz).powi(2)).sqrt() - radius;
                let idx = fmm.flat(ix, iy, iz);
                exact[idx] = signed;
                if signed.abs() < half_width {
                    seeds.push((idx, signed));
                    seeded[idx] = true;
                }
            }
        }
    }
    fmm.set_known(&seeds);
    fmm.run();
    (fmm, exact, seeded)
}

#[test]
fn inside_seeds_yield_negative_interior_values() {
    let (fmm, exact, seeded) = sphere_band(0.4);
    let mut interior_cells = 0;
    for (idx, (&reference, &is_seed)) in exact.iter().zip(&seeded).enumerate() {
        if is_seed {
            continue;
        }
        let d = fmm.dist[idx];
        assert_eq!(
            d < 0.0,
            reference < 0.0,
            "cell {idx}: sign of {d} disagrees with the exact {reference}"
        );
        assert!(
            (d - reference).abs() < 0.25,
            "cell {idx}: {d} vs exact {reference}"
        );
        if reference < 0.0 {
            interior_cells += 1;
        }
    }
    assert!(interior_cells > 100, "only {interior_cells} interior cells");
    // The sphere centre is 1.5 below the surface.
    let centre = fmm.distance_at(10, 10, 10);
    assert!((centre + 1.5).abs() < 0.25, "centre value {centre}");
    // A far corner of the box is outside.
    let corner = fmm.distance_at(0, 0, 0);
    assert!(corner > 2.0, "corner value {corner}");
}

#[test]
fn a_negative_seed_keeps_its_sign_while_the_magnitude_grows() {
    let mut fmm = FastMarchingMethod::new(5, 5, 5, 1.0);
    let seed = fmm.flat(2, 2, 2);
    fmm.set_known(&[(seed, -0.5)]);
    fmm.run();
    assert!((fmm.distance_at(2, 2, 2) + 0.5).abs() < 1e-15);
    assert!((fmm.distance_at(3, 2, 2) + 1.5).abs() < 1e-12);
    assert!((fmm.distance_at(2, 2, 4) + 2.5).abs() < 1e-12);
    for iz in 0..5 {
        for iy in 0..5 {
            for ix in 0..5 {
                assert!(fmm.distance_at(ix, iy, iz) < 0.0);
            }
        }
    }
}

#[test]
fn opposite_signed_seeds_split_a_line_at_their_interface() {
    let mut fmm = FastMarchingMethod::new(9, 1, 1, 1.0);
    let inside = fmm.flat(3, 0, 0);
    let outside = fmm.flat(5, 0, 0);
    fmm.set_known(&[(inside, -0.5), (outside, 0.5)]);
    fmm.run();
    for (ix, expected) in [
        (2, -1.5),
        (1, -2.5),
        (0, -3.5),
        (6, 1.5),
        (7, 2.5),
        (8, 3.5),
    ] {
        let d = fmm.distance_at(ix, 0, 0);
        assert!((d - expected).abs() < 1e-12, "cell {ix}: {d} vs {expected}");
    }
}

#[test]
fn without_seeds_nothing_is_propagated() {
    let mut fmm = FastMarchingMethod::new(3, 3, 3, 1.0);
    fmm.run();
    assert!(fmm.dist.iter().all(|&d| d == f64::MAX));
}

#[test]
fn invalid_seeds_are_ignored() {
    let mut fmm = FastMarchingMethod::new(3, 3, 3, 1.0);
    let centre = fmm.flat(1, 1, 1);
    fmm.set_known(&[(1000, 0.0), (centre, f64::NAN), (centre, f64::INFINITY)]);
    fmm.run();
    assert!(fmm.dist.iter().all(|&d| d == f64::MAX));
}
