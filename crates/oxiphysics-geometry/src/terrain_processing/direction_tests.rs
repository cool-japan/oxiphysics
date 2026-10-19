//! Regression tests for the compass aspect ([`compute_aspect`]) and the
//! D-infinity flow direction ([`flow_direction_dinf`]): every cardinal and
//! diagonal orientation, every facet of the D-infinity table, the clamping
//! branches, the sentinel for cells without a downslope facet, and the solar
//! irradiance that depends on the aspect.

use super::tests::{cone_dem, flat_dem};
use super::*;

/// Smallest absolute difference between two angles on the circle.
fn angle_diff(a: f64, b: f64) -> f64 {
    let d = (a - b).rem_euclid(2.0 * PI);
    d.min(2.0 * PI - d)
}

/// Plane `z = dz_dcol * col + dz_drow * row` on a `rows x cols` grid.
fn plane_dem(rows: usize, cols: usize, cell_size: f64, dz_dcol: f64, dz_drow: f64) -> Dem {
    let heights: Vec<f64> = (0..rows)
        .flat_map(|r| (0..cols).map(move |c| dz_dcol * c as f64 + dz_drow * r as f64))
        .collect();
    Dem::new(heights, rows, cols, cell_size)
}

/// Plane that descends towards the compass-free math angle `theta`
/// (anticlockwise from east, rows increase southward) with the given
/// rise-over-run `gradient`.
fn plane_dipping_towards(rows: usize, cols: usize, theta: f64, gradient: f64) -> Dem {
    let cell_size = 10.0;
    plane_dem(
        rows,
        cols,
        cell_size,
        -gradient * cell_size * theta.cos(),
        gradient * cell_size * theta.sin(),
    )
}

#[test]
fn test_aspect_cardinal_directions() {
    // z = row rises towards the south, so the slope faces north (0).
    let cases = [
        (0.0, 1.0, 0.0, "rises south -> faces north"),
        (0.0, -1.0, PI, "rises north -> faces south"),
        (1.0, 0.0, 3.0 * PI / 2.0, "rises east -> faces west"),
        (-1.0, 0.0, PI / 2.0, "rises west -> faces east"),
    ];
    for &(dz_dcol, dz_drow, expected, label) in &cases {
        let dem = plane_dem(5, 5, 10.0, dz_dcol, dz_drow);
        let aspect = compute_aspect(&dem);
        for (i, &a) in aspect.iter().enumerate() {
            assert!(
                angle_diff(a, expected) < 1e-12 && (0.0..2.0 * PI).contains(&a),
                "{label}: cell {i} aspect {a}, expected {expected}"
            );
        }
    }
}

#[test]
fn test_aspect_diagonal_directions() {
    let cases = [
        (
            1.0,
            1.0,
            7.0 * PI / 4.0,
            "rises south-east -> faces north-west",
        ),
        (
            -1.0,
            -1.0,
            3.0 * PI / 4.0,
            "rises north-west -> faces south-east",
        ),
        (
            1.0,
            -1.0,
            5.0 * PI / 4.0,
            "rises north-east -> faces south-west",
        ),
        (-1.0, 1.0, PI / 4.0, "rises south-west -> faces north-east"),
    ];
    for &(dz_dcol, dz_drow, expected, label) in &cases {
        let dem = plane_dem(5, 5, 10.0, dz_dcol, dz_drow);
        let aspect = compute_aspect(&dem);
        let centre = aspect[2 * 5 + 2];
        assert!(
            angle_diff(centre, expected) < 1e-12,
            "{label}: aspect {centre}, expected {expected}"
        );
    }
}

#[test]
fn test_aspect_flat_is_minus_one() {
    let aspect = compute_aspect(&flat_dem(4, 4, 25.0));
    assert!(aspect.iter().all(|&a| a == -1.0));
}

#[test]
fn test_aspect_north_facing_slope_gets_more_morning_sun_from_the_north() {
    // Rising towards the south = facing north. A sun in the north (azimuth
    // 0) at 30 degrees elevation lights it far better than the mirrored
    // south-facing slope.
    let solar = SolarPosition::from_degrees(30.0, 0.0);
    let north_facing = plane_dem(5, 5, 10.0, 0.0, 5.0);
    let south_facing = plane_dem(5, 5, 10.0, 0.0, -5.0);
    let irradiance = |dem: &Dem| {
        let slope = compute_slope(dem);
        let aspect = compute_aspect(dem);
        terrain_solar_irradiance(dem, &slope, &aspect, &solar, 1000.0, 0.0)[2 * 5 + 2]
    };
    let beta = 0.5f64.atan();
    let zenith = solar.zenith;
    let expected_north = 1000.0 * (zenith - beta).cos();
    let expected_south = 1000.0 * (zenith + beta).cos();
    let north = irradiance(&north_facing);
    let south = irradiance(&south_facing);
    assert!(
        (north - expected_north).abs() < 1e-9,
        "north-facing {north}, expected {expected_north}"
    );
    assert!(
        (south - expected_south).abs() < 1e-9,
        "south-facing {south}, expected {expected_south}"
    );
    assert!(north > south + 500.0);

    // East/west mirror with the sun in the east.
    let east_sun = SolarPosition::from_degrees(30.0, 90.0);
    let faces_east = plane_dem(5, 5, 10.0, -5.0, 0.0);
    let faces_west = plane_dem(5, 5, 10.0, 5.0, 0.0);
    let irr_east_sun = |dem: &Dem| {
        let slope = compute_slope(dem);
        let aspect = compute_aspect(dem);
        terrain_solar_irradiance(dem, &slope, &aspect, &east_sun, 1000.0, 0.0)[2 * 5 + 2]
    };
    assert!(
        (irr_east_sun(&faces_east) - expected_north).abs() < 1e-9,
        "east-facing slope under an eastern sun"
    );
    assert!(
        (irr_east_sun(&faces_west) - expected_south).abs() < 1e-9,
        "west-facing slope under an eastern sun"
    );
}

#[test]
fn test_dinf_plane_dipping_towards_every_facet() {
    let mut thetas: Vec<f64> = (0..24).map(|k| (k as f64 * 15.0).to_radians()).collect();
    // Awkward angles (including the ones inside each odd/even facet).
    thetas.extend(
        [
            10.0f64, 30.0, 37.0, 120.0, 123.0, 210.0, 255.0, 300.0, 359.0,
        ]
        .iter()
        .map(|d| d.to_radians()),
    );
    for &theta in &thetas {
        let dem = plane_dipping_towards(7, 7, theta, 0.4);
        let fdir = flow_direction_dinf(&dem);
        for r in 1..6 {
            for c in 1..6 {
                let a = fdir[r * 7 + c];
                assert!(
                    (0.0..2.0 * PI).contains(&a),
                    "theta {theta}: cell ({r}, {c}) angle {a} out of range"
                );
                assert!(
                    angle_diff(a, theta) < 1e-9,
                    "theta {:.3} deg: cell ({r}, {c}) got {:.6} deg",
                    theta.to_degrees(),
                    a.to_degrees()
                );
            }
        }
    }
}

#[test]
fn test_dinf_plane_30_degrees_is_30_degrees() {
    let dem = plane_dipping_towards(5, 5, 30.0f64.to_radians(), 0.5);
    let fdir = flow_direction_dinf(&dem);
    assert!((fdir[2 * 5 + 2].to_degrees() - 30.0).abs() < 1e-9);
}

#[test]
fn test_dinf_clamps_to_cardinal_edge_when_angle_is_negative() {
    // Centre 0, east neighbour -1 but both east-side diagonals at 0: the
    // plane's steepest direction lies outside both facets, so the flow
    // follows the cardinal edge exactly. Everything else is uphill.
    let grid = [[5.0, 5.0, 0.0], [5.0, 0.0, -1.0], [5.0, 5.0, 0.0]];
    let dem = Dem::new(grid.concat(), 3, 3, 1.0);
    let fdir = flow_direction_dinf(&dem);
    assert_eq!(fdir[4], 0.0);
}

#[test]
fn test_dinf_clamps_to_diagonal_edge_when_angle_exceeds_facet() {
    // East neighbour barely lower, north-east neighbour much lower: the
    // steepest direction points beyond the diagonal, so flow follows the
    // diagonal exactly (north-east = π/4).
    let grid = [[5.0, 5.0, -1.0], [5.0, 0.0, -0.1], [5.0, 5.0, 5.0]];
    let dem = Dem::new(grid.concat(), 3, 3, 1.0);
    let fdir = flow_direction_dinf(&dem);
    assert!((fdir[4] - PI / 4.0).abs() < 1e-12, "got {}", fdir[4]);
}

#[test]
fn test_dinf_pits_and_flats_have_no_direction() {
    assert!(
        flow_direction_dinf(&flat_dem(4, 4, 10.0))
            .iter()
            .all(|&a| a == -1.0)
    );
    // A single-cell pit: every neighbour is higher.
    let mut heights = vec![10.0; 9];
    heights[4] = 1.0;
    let dem = Dem::new(heights, 3, 3, 1.0);
    assert_eq!(flow_direction_dinf(&dem)[4], -1.0);
    // A single column never has a full facet inside the grid.
    let column = Dem::new(vec![3.0, 2.0, 1.0], 3, 1, 1.0);
    assert!(flow_direction_dinf(&column).iter().all(|&a| a == -1.0));
}

#[test]
fn test_dinf_cone_flows_away_from_the_peak() {
    let dem = cone_dem(9, 9, 5.0);
    let fdir = flow_direction_dinf(&dem);
    // Two cells from the peak along each compass direction; the angle must
    // stay within one facet (π/4) of the radial direction.
    let cases = [
        (4usize, 6usize, 0.0, "east"),
        (2, 4, PI / 2.0, "north"),
        (4, 2, PI, "west"),
        (6, 4, 3.0 * PI / 2.0, "south"),
    ];
    for &(r, c, expected, label) in &cases {
        let a = fdir[r * 9 + c];
        assert!(
            angle_diff(a, expected) < PI / 4.0,
            "{label} of the peak: flow angle {:.1} deg",
            a.to_degrees()
        );
    }
}

#[test]
fn test_aspect_and_dinf_describe_the_same_downslope_direction() {
    // Compass aspect (clockwise from north) and D-infinity angle
    // (anticlockwise from east) must satisfy aspect + angle = π/2 (mod 2π).
    for k in 0..16 {
        let theta = (k as f64 * 22.5 + 7.0).to_radians();
        let dem = plane_dipping_towards(7, 7, theta, 0.3);
        let aspect = compute_aspect(&dem)[3 * 7 + 3];
        let angle = flow_direction_dinf(&dem)[3 * 7 + 3];
        assert!(
            angle_diff(aspect + angle, PI / 2.0) < 1e-9,
            "theta {:.1} deg: aspect {:.3} + dinf {:.3}",
            theta.to_degrees(),
            aspect,
            angle
        );
    }
}
