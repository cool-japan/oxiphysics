//! Regression tests for `FastMarchingMethod::reinitialize`.
//!
//! The quadratic level sets used here are deliberately not distance functions
//! (the gradient of `|x|^2 - R^2` has magnitude `2 |x|`): seeding the march with
//! the raw level-set values, or with an exact signed distance function, would
//! pass without ever exercising the `|phi| / |grad phi|` normalisation.

use super::FastMarchingMethod;
use crate::Error;

/// Cells per axis of the test grid.
const N: usize = 21;
/// Grid spacing.
const H: f64 = 0.25;
/// Coordinate of lattice index 0 on every axis (index 10 is the origin).
const ORIGIN: f64 = -2.5;
/// Sphere radius used by the quadratic level sets.
const RADIUS: f64 = 1.5;

fn coord(i: usize) -> f64 {
    ORIGIN + i as f64 * H
}

fn solver() -> FastMarchingMethod {
    FastMarchingMethod::new(N, N, N, H)
}

/// Evaluates `field` at every lattice point, in the flat grid layout (x fastest).
fn sample(field: impl Fn([f64; 3]) -> f64) -> Vec<f64> {
    let mut values = Vec::with_capacity(N * N * N);
    for iz in 0..N {
        for iy in 0..N {
            for ix in 0..N {
                values.push(field([coord(ix), coord(iy), coord(iz)]));
            }
        }
    }
    values
}

fn sphere_level_set(p: [f64; 3]) -> f64 {
    p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - RADIUS * RADIUS
}

/// Exact signed distance to the sphere of radius `RADIUS` around the origin.
fn sphere_distance(p: [f64; 3]) -> f64 {
    (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - RADIUS
}

fn position(index: usize) -> [f64; 3] {
    [
        coord(index % N),
        coord((index / N) % N),
        coord(index / (N * N)),
    ]
}

/// Whether a face neighbour of cell `index` has a strictly opposite sign.
fn has_sign_change(phi: &[f64], index: usize) -> bool {
    let (ix, iy, iz) = (index % N, (index / N) % N, index / (N * N));
    let mut neighbours = Vec::new();
    if ix > 0 {
        neighbours.push(index - 1);
    }
    if ix + 1 < N {
        neighbours.push(index + 1);
    }
    if iy > 0 {
        neighbours.push(index - N);
    }
    if iy + 1 < N {
        neighbours.push(index + N);
    }
    if iz > 0 {
        neighbours.push(index - N * N);
    }
    if iz + 1 < N {
        neighbours.push(index + N * N);
    }
    neighbours
        .iter()
        .any(|&other| (phi[other] < 0.0) != (phi[index] < 0.0) && phi[other] != 0.0)
}

#[test]
fn quadratic_level_set_gives_the_signed_distance_of_the_sphere() {
    let phi = sample(sphere_level_set);
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");

    // The centre is R below the surface. The first-order upwind update rounds off
    // the cone of the distance function there (the three-term update of six
    // neighbours at R - dx gives R - 0.1 or so), so the centre is a little shallower.
    let centre = fmm.distance_at(10, 10, 10);
    assert!(
        (centre + RADIUS).abs() < 0.25,
        "centre value {centre}, expected about {}",
        -RADIUS
    );

    // Outside the sphere the march overestimates by a few percent.
    let corner = fmm.distance_at(0, 0, 0);
    let exact_corner = sphere_distance(position(0));
    assert!(
        corner > exact_corner - 0.02 && corner < exact_corner * 1.08,
        "corner value {corner}, exact {exact_corner}"
    );

    // Manhattan-type distances would be off by more than a grid cell far from the
    // sphere; the Eikonal march stays within a quarter of the radius.
    for (index, &d) in fmm.dist.iter().enumerate() {
        let exact = sphere_distance(position(index));
        assert!(
            (d - exact).abs() < 0.25,
            "cell {index}: {d} vs exact {exact}"
        );
    }
}

#[test]
fn distance_signs_follow_the_level_set_everywhere() {
    let phi = sample(sphere_level_set);
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");
    let mut inside = 0;
    for (index, (&d, &p)) in fmm.dist.iter().zip(&phi).enumerate() {
        assert!(
            d.is_finite() && d != f64::MAX,
            "cell {index} was not reached"
        );
        assert_eq!(d < 0.0, p < 0.0, "cell {index}: distance {d}, phi {p}");
        if p < 0.0 {
            inside += 1;
        }
    }
    assert!(inside > 100, "only {inside} cells inside the sphere");
}

#[test]
fn positive_scaling_of_the_level_set_does_not_change_the_result() {
    let phi = sample(sphere_level_set);
    let mut reference = solver();
    reference.reinitialize(&phi).expect("grid-sized field");
    for factor in [10.0, 0.01, 1234.5] {
        let scaled: Vec<f64> = phi.iter().map(|p| p * factor).collect();
        let mut fmm = solver();
        fmm.reinitialize(&scaled).expect("grid-sized field");
        for (index, (&a, &b)) in fmm.dist.iter().zip(&reference.dist).enumerate() {
            assert!(
                (a - b).abs() < 1e-9,
                "factor {factor}, cell {index}: {a} vs {b}"
            );
        }
    }
}

#[test]
fn an_exact_distance_field_is_reproduced() {
    let phi = sample(sphere_distance);
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");
    for (index, &d) in fmm.dist.iter().enumerate() {
        let exact = sphere_distance(position(index));
        assert!(
            (d - exact).abs() < 0.25,
            "cell {index}: {d} vs exact {exact}"
        );
    }
}

#[test]
fn axis_aligned_plane_between_lattice_points_is_exact_everywhere() {
    // Plane x = 0.1 lies between the lattice planes x = 0 and x = 0.25; the
    // level set is not normalised (its gradient has magnitude 3).
    let phi = sample(|p| 3.0 * (p[0] - 0.1));
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");
    for (index, &d) in fmm.dist.iter().enumerate() {
        let exact = position(index)[0] - 0.1;
        assert!(
            (d - exact).abs() < 1e-12,
            "cell {index}: {d} vs exact {exact}"
        );
    }
}

#[test]
fn oblique_plane_seeds_carry_the_exact_perpendicular_distance() {
    // 2x + y = 0.3 is oblique to every grid axis. The axis-crossing distance
    // alone would overestimate; the gradient-based seed is exact for a plane.
    let phi = sample(|p| 2.0 * p[0] + p[1] - 0.3);
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");
    let norm = 5.0_f64.sqrt();
    let mut seeds = 0;
    for (index, (&d, &p)) in fmm.dist.iter().zip(&phi).enumerate() {
        let exact = p / norm;
        if has_sign_change(&phi, index) {
            seeds += 1;
            assert!(
                (d - exact).abs() < 1e-12,
                "seed {index}: {d} vs exact {exact}"
            );
        } else {
            // Marched cells never underestimate the distance to the plane and stay
            // within the first-order error of the scheme.
            assert!(
                d.abs() > exact.abs() - 1e-9 && d.abs() < 1.2 * exact.abs() + H,
                "cell {index}: {d} vs exact {exact}"
            );
            assert_eq!(d < 0.0, p < 0.0, "cell {index}");
        }
    }
    assert!(seeds > 40, "only {seeds} seed cells");
}

#[test]
fn cells_exactly_on_the_interface_are_seeded_with_zero() {
    // The plane x = 0 passes through the lattice plane ix = 10.
    let phi = sample(|p| p[0]);
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");
    for (index, &d) in fmm.dist.iter().enumerate() {
        let x = position(index)[0];
        assert!((d - x).abs() < 1e-12, "cell {index}: {d} vs {x}");
        if index % N == 10 {
            assert_eq!(d, 0.0, "cell {index} lies on the interface");
        }
    }
}

#[test]
fn reinitialize_discards_earlier_distances() {
    let mut used = solver();
    let seed = used.flat(3, 4, 5);
    used.set_known(&[(seed, 0.0)]);
    used.run();

    let phi = sample(sphere_level_set);
    used.reinitialize(&phi).expect("grid-sized field");
    let mut fresh = solver();
    fresh.reinitialize(&phi).expect("grid-sized field");
    assert_eq!(used.dist, fresh.dist);

    // The solver is reusable: a second call gives the same answer again.
    used.reinitialize(&phi).expect("grid-sized field");
    assert_eq!(used.dist, fresh.dist);
}

#[test]
fn wrong_field_length_is_rejected() {
    let mut fmm = solver();
    for length in [0, 5, N * N * N - 1, N * N * N + 1] {
        match fmm.reinitialize(&vec![1.0; length]) {
            Err(Error::LengthMismatch { expected, actual }) => {
                assert_eq!(expected, N * N * N);
                assert_eq!(actual, length);
            }
            other => panic!("expected a length mismatch, got {other:?}"),
        }
    }
    // A rejected call must not disturb the stored distances.
    assert!(fmm.dist.iter().all(|&d| d == f64::MAX));
}

#[test]
fn invalid_grid_spacing_is_rejected() {
    for dx in [0.0, -0.25, f64::NAN, f64::INFINITY] {
        let mut fmm = FastMarchingMethod::new(3, 3, 3, dx);
        assert!(
            fmm.reinitialize(&[1.0; 27]).is_err(),
            "spacing {dx} must be rejected"
        );
    }
}

#[test]
fn a_field_without_sign_change_has_no_interface() {
    let mut fmm = solver();
    fmm.reinitialize(&vec![2.5; N * N * N])
        .expect("grid-sized field");
    assert!(fmm.dist.iter().all(|&d| d == f64::MAX));
    fmm.reinitialize(&vec![-0.5; N * N * N])
        .expect("grid-sized field");
    assert!(fmm.dist.iter().all(|&d| d == f64::MAX));
}

#[test]
fn non_finite_values_are_never_seeded() {
    let phi = sample(sphere_level_set);
    let mut reference = solver();
    reference.reinitialize(&phi).expect("grid-sized field");

    let mut damaged = phi.clone();
    damaged[0] = f64::NAN;
    damaged[N * N * N - 1] = f64::INFINITY;
    let mut fmm = solver();
    fmm.reinitialize(&damaged).expect("grid-sized field");
    for (index, (&a, &b)) in fmm.dist.iter().zip(&reference.dist).enumerate() {
        assert!(
            a.is_finite() && a != f64::MAX,
            "cell {index} was not reached"
        );
        assert!((a - b).abs() < 1e-12, "cell {index}: {a} vs {b}");
    }
}

#[test]
fn a_single_sign_change_cell_is_enough_to_seed_the_march() {
    // One negative cell in a positive field: its six neighbours and itself form
    // the interface band, and every cell of the grid is reached.
    let mut phi = vec![1.0; N * N * N];
    let centre = 10 + N * (10 + N * 10);
    phi[centre] = -1.0;
    let mut fmm = solver();
    fmm.reinitialize(&phi).expect("grid-sized field");
    assert!(fmm.dist[centre] < 0.0);
    assert!(fmm.dist[centre - 1] > 0.0);
    assert!(fmm.dist.iter().all(|&d| d != f64::MAX));
    let far = fmm.distance_at(0, 0, 0);
    assert!(far > 4.0 && far < 5.0, "far corner {far}");
}

#[test]
fn a_point_seed_gives_euclidean_rather_than_manhattan_distances() {
    // A single zero-distance seed: the cell three steps along every axis is at
    // Euclidean distance sqrt(27) = 5.196 (the L1 distance would be 9). The
    // first-order scheme overestimates the diagonal by about 18 % three cells
    // from a point source and by about 8 % at ten cells.
    let mut fmm = solver();
    let seed = fmm.flat(10, 10, 10);
    fmm.set_known(&[(seed, 0.0)]);
    fmm.run();
    let d = fmm.distance_at(13, 13, 13);
    let exact = 27.0_f64.sqrt() * H;
    assert!(
        d >= exact - 1e-9 && d < 1.20 * exact,
        "diagonal distance {d}, Euclidean {exact}"
    );
    let far = fmm.distance_at(20, 20, 20);
    let far_exact = 300.0_f64.sqrt() * H;
    assert!(far < 1.10 * far_exact, "far diagonal {far} vs {far_exact}");
}

#[test]
fn grids_with_a_single_layer_or_row_are_supported() {
    // A 9 x 5 x 1 slab: the gradient only has x and y components and the
    // missing z neighbours must not disturb the seeds or the march.
    let (nx, ny, nz, h) = (9, 5, 1, 0.5);
    let mut fmm = FastMarchingMethod::new(nx, ny, nz, h);
    let phi: Vec<f64> = (0..nx * ny * nz)
        .map(|index| 2.0 * ((index % nx) as f64 * h - 0.3))
        .collect();
    fmm.reinitialize(&phi).expect("grid-sized field");
    for (index, &d) in fmm.dist.iter().enumerate() {
        let exact = (index % nx) as f64 * h - 0.3;
        assert!(
            (d - exact).abs() < 1e-12,
            "cell {index}: {d} vs exact {exact}"
        );
    }

    // A single row.
    let mut row = FastMarchingMethod::new(7, 1, 1, 1.0);
    let phi: Vec<f64> = (0..7).map(|i| i as f64 - 2.4).collect();
    row.reinitialize(&phi).expect("grid-sized field");
    for (i, &d) in row.dist.iter().enumerate() {
        assert!((d - (i as f64 - 2.4)).abs() < 1e-12, "cell {i}: {d}");
    }
}
