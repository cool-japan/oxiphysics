//! Regression tests for the ordered lattice packings.
//!
//! The HCP basis lives in the orthohexagonal cell, whose edges differ in length.
//! An earlier version repeated it on a cubic grid of edge `2*r*sqrt(2)`, which
//! made neighbouring spheres interpenetrate.

use super::functions::has_overlap;
use super::types::{LatticeType, OrderedPacking, PackingSphere};

const ALL_LATTICES: [LatticeType; 4] = [
    LatticeType::Fcc,
    LatticeType::Hcp,
    LatticeType::Bcc,
    LatticeType::SimpleCubic,
];

fn generated(lattice: LatticeType, radius: f64, n_cells: [usize; 3]) -> OrderedPacking {
    let mut packing = OrderedPacking::new(lattice, radius, n_cells);
    packing.generate();
    packing
}

/// Smallest centre-to-centre distance over all pairs of spheres.
fn minimum_pair_distance(spheres: &[PackingSphere]) -> f64 {
    let mut minimum = f64::INFINITY;
    for (i, a) in spheres.iter().enumerate() {
        for b in &spheres[i + 1..] {
            minimum = minimum.min(a.distance_to(b));
        }
    }
    minimum
}

#[test]
fn hcp_spheres_do_not_overlap() {
    for radius in [1.0, 0.37, 12.5] {
        let packing = generated(LatticeType::Hcp, radius, [3, 2, 3]);
        assert_eq!(packing.spheres.len(), 4 * 3 * 2 * 3);
        assert!(
            !has_overlap(&packing.spheres),
            "HCP spheres of radius {radius} interpenetrate"
        );
    }
}

#[test]
fn hcp_minimum_pair_distance_is_two_radii() {
    for radius in [1.0, 0.37, 12.5] {
        let packing = generated(LatticeType::Hcp, radius, [3, 2, 3]);
        let minimum = minimum_pair_distance(&packing.spheres);
        assert!(
            (minimum - 2.0 * radius).abs() < 1e-12 * radius,
            "minimum distance {minimum} for radius {radius}"
        );
    }
}

#[test]
fn hcp_packing_fraction_is_the_close_packed_limit() {
    let packing = generated(LatticeType::Hcp, 1.0, [2, 2, 2]);
    let measured = packing.actual_packing_fraction();
    assert!((measured - 0.7405).abs() < 1e-4, "fraction {measured}");
    assert!((measured - packing.theoretical_packing_fraction()).abs() < 1e-12);
}

#[test]
fn hcp_interior_sphere_has_twelve_nearest_neighbours() {
    let packing = generated(LatticeType::Hcp, 1.0, [4, 3, 4]);
    let [lx, ly, lz] = packing.cell_lengths();
    let centre = [2.0 * lx, 1.5 * ly, 2.0 * lz];
    let interior = packing
        .spheres
        .iter()
        .min_by(|a, b| {
            let da = (0..3)
                .map(|k| (a.center[k] - centre[k]).powi(2))
                .sum::<f64>();
            let db = (0..3)
                .map(|k| (b.center[k] - centre[k]).powi(2))
                .sum::<f64>();
            da.total_cmp(&db)
        })
        .expect("the packing contains spheres");
    let index = packing
        .spheres
        .iter()
        .position(|s| s.id == interior.id)
        .expect("the sphere belongs to the packing");
    let neighbours = packing.find_neighbors(index, 2.0 * packing.radius * (1.0 + 1e-9));
    assert_eq!(neighbours.len(), packing.coordination_number());
}

#[test]
fn hcp_cell_is_the_orthohexagonal_cell() {
    let radius = 1.5;
    let packing = OrderedPacking::new(LatticeType::Hcp, radius, [1, 1, 1]);
    assert!((packing.lattice_parameter - 2.0 * radius).abs() < 1e-15);
    let [lx, ly, lz] = packing.cell_lengths();
    assert!((lx - 2.0 * radius).abs() < 1e-12);
    assert!((ly - 2.0 * 3.0_f64.sqrt() * radius).abs() < 1e-12);
    assert!((lz - 2.0 * radius * (8.0_f64 / 3.0).sqrt()).abs() < 1e-12);
}

#[test]
fn hcp_layers_are_stacked_abab() {
    let packing = generated(LatticeType::Hcp, 1.0, [1, 1, 3]);
    let [_, _, lz] = packing.cell_lengths();
    let layer_of = |s: &PackingSphere| (s.center[2] / (0.5 * lz)).round() as i64;
    let in_layer = |layer: i64| -> Vec<[f64; 2]> {
        packing
            .spheres
            .iter()
            .filter(|s| layer_of(s) == layer)
            .map(|s| [s.center[0], s.center[1]])
            .collect()
    };
    let same_site =
        |a: [f64; 2], b: [f64; 2]| (a[0] - b[0]).abs() < 1e-9 && (a[1] - b[1]).abs() < 1e-9;
    // Layer 2 sits directly above layer 0 (A-B-A), layer 1 is shifted (B).
    for site in in_layer(0) {
        assert!(in_layer(2).iter().any(|other| same_site(site, *other)));
        assert!(!in_layer(1).iter().any(|other| same_site(site, *other)));
    }
}

#[test]
fn cubic_lattices_repeat_a_cube() {
    for lattice in [LatticeType::Fcc, LatticeType::Bcc, LatticeType::SimpleCubic] {
        let packing = OrderedPacking::new(lattice, 1.3, [1, 1, 1]);
        let a = packing.lattice_parameter;
        assert_eq!(packing.cell_lengths(), [a, a, a]);
    }
}

#[test]
fn every_lattice_reaches_its_theoretical_packing_fraction() {
    for lattice in ALL_LATTICES {
        let packing = generated(lattice, 1.0, [2, 2, 2]);
        let measured = packing.actual_packing_fraction();
        assert!(
            (measured - packing.theoretical_packing_fraction()).abs() < 1e-12,
            "{lattice:?}: measured {measured}, theory {}",
            packing.theoretical_packing_fraction()
        );
    }
}

#[test]
fn every_lattice_touches_at_exactly_two_radii() {
    for lattice in ALL_LATTICES {
        let packing = generated(lattice, 1.0, [2, 2, 2]);
        let minimum = minimum_pair_distance(&packing.spheres);
        assert!(
            (minimum - 2.0 * packing.radius).abs() < 1e-12,
            "{lattice:?}: minimum distance {minimum}"
        );
    }
}

#[test]
fn empty_cell_counts_generate_no_spheres() {
    let packing = generated(LatticeType::Hcp, 1.0, [0, 2, 2]);
    assert!(packing.spheres.is_empty());
    assert_eq!(packing.actual_packing_fraction(), 0.0);
}

#[test]
fn touching_spheres_do_not_overlap_but_penetration_does() {
    let a = PackingSphere::new([0.0; 3], 1.0, 0);
    let touching = PackingSphere::new([2.0, 0.0, 0.0], 1.0, 1);
    let rounded = PackingSphere::new([2.0 * (1.0 - 1e-14), 0.0, 0.0], 1.0, 2);
    let penetrating = PackingSphere::new([2.0 * (1.0 - 1e-6), 0.0, 0.0], 1.0, 3);
    assert!(!a.overlaps(&touching));
    assert!(!a.overlaps(&rounded));
    assert!(a.overlaps(&penetrating));
    assert!(!has_overlap(&[a, touching]));
}
