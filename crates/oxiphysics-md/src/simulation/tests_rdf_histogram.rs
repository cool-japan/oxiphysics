// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the radial distribution function accumulator
//! [`RdfHistogram`].
//!
//! The histogram counts *unique* pairs `i < j`, so the ideal-gas reference of one
//! frame is `N (N − 1) / 2 · v_shell / V` pairs per shell.  The former reference
//! `ρ · v_shell · n_frames · n_pairs` was a factor `N` too large: g(r) of an ideal
//! gas came out as `1 / N`, and the coordination number summed `ρ · v_shell ·
//! hist / norm` instead of counting the neighbours of an atom (and its
//! `n_atoms * (n_atoms − 1)` underflowed for an empty system).

use std::f64::consts::PI;

use super::stochastic::SimRng;
use super::*;

/// Accumulate `frames` frames of `n` uniformly random positions in a cube of edge
/// `edge` into an [`RdfHistogram`] with bins of width 1 Å up to 10 Å.
fn ideal_gas_histogram(n: usize, edge: f64, frames: usize, seed: u64) -> RdfHistogram {
    let mut rng = SimRng::new(seed);
    let mut rdf = RdfHistogram::new(10.0, 10, edge * edge * edge, n);
    for _ in 0..frames {
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|_| {
                [
                    rng.uniform() * edge,
                    rng.uniform() * edge,
                    rng.uniform() * edge,
                ]
            })
            .collect();
        rdf.accumulate(&positions, edge);
    }
    rdf
}

/// Volume of a ball of radius `r`.
fn ball_volume(r: f64) -> f64 {
    4.0 / 3.0 * PI * r.powi(3)
}

/// The `n³` sites of a simple cubic lattice with spacing `a`.
fn simple_cubic(n: usize, a: f64) -> Vec<[f64; 3]> {
    let mut sites = Vec::with_capacity(n * n * n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                sites.push([i as f64 * a, j as f64 * a, k as f64 * a]);
            }
        }
    }
    sites
}

#[test]
fn gofr_of_an_ideal_gas_is_one() {
    // Uncorrelated positions have g(r) = 1.  A reference that is a factor N too
    // large gives g = 1 / N = 0.001.
    let rdf = ideal_gas_histogram(1000, 30.0, 20, 2024);
    let g = rdf.gofr();
    let centres = rdf.bin_centres();
    assert_eq!(g.len(), 10);
    // Bin 0 (r < 1 Å) holds only ~1500 pairs; the others hold 10^4 and more.
    for (k, (&value, &r)) in g.iter().zip(centres.iter()).enumerate().skip(1) {
        assert!(
            (value - 1.0).abs() < 0.05,
            "bin {k}: g({r}) = {value}, expected 1 within 5 %"
        );
    }
}

#[test]
fn gofr_does_not_depend_on_the_number_of_atoms() {
    // Same number density, different N: g(r) of an ideal gas is 1 in both.
    let small = ideal_gas_histogram(200, 17.8, 400, 7);
    let large = ideal_gas_histogram(1000, 30.0, 20, 8);
    for rdf in [&small, &large] {
        let g = rdf.gofr();
        // 8 Å shells stay inside the 17.8 Å cube (r <= L / 2 = 8.9 Å).
        for (k, &value) in g.iter().enumerate().take(8).skip(1) {
            assert!(
                (value - 1.0).abs() < 0.06,
                "N = {}: bin {k}: g = {value}",
                rdf.n_atoms
            );
        }
    }
}

#[test]
fn coordination_number_of_an_ideal_gas_is_the_density_times_the_ball_volume() {
    // An atom has (N − 1) / V · 4π/3 r³ neighbours within r on average:
    // 9.92 within 4 Å for 1000 atoms in a 30 Å box.
    let (n, edge) = (1000_usize, 30.0_f64);
    let rdf = ideal_gas_histogram(n, edge, 20, 2024);
    let density = (n - 1) as f64 / edge.powi(3);
    for r_cut in [2.0, 3.0, 4.0, 5.0, 6.0] {
        let expected = density * ball_volume(r_cut);
        let coordination = rdf.coordination_number(r_cut);
        assert!(
            (coordination - expected).abs() < 0.03 * expected,
            "CN({r_cut}) = {coordination}, expected {expected}"
        );
    }
    let at_four = rdf.coordination_number(4.0);
    assert!((at_four - 9.92).abs() < 0.3, "CN(4) = {at_four}");
}

#[test]
fn coordination_number_counts_the_neighbours_of_a_simple_cubic_lattice() {
    // 6³ sites, spacing 2 Å, periodic box 12 Å.  The shells of the lattice hold
    // 6, 12, 8, 6, 24, 24, 12 neighbours at 2, 2.83, 3.46, 4, 4.47, 4.90 and 5.66 Å
    // (7 is not a sum of three squares, so there is no shell at 5.29 Å).  With bins
    // of 0.5 Å the bins [2, 2.5), [2.5, 3), [3, 3.5), [4, 4.5), [4.5, 5) and
    // [5.5, 6) hold 6, 12, 8, 30, 24 and 12 neighbours.
    let (sites, a) = (6, 2.0);
    let positions = simple_cubic(sites, a);
    let box_len = sites as f64 * a;
    let mut rdf = RdfHistogram::new(6.0, 12, box_len.powi(3), positions.len());
    rdf.accumulate(&positions, box_len);
    // Each unique pair is counted once, so the pair count of a shell is
    // n_atoms * neighbours / 2.
    let n_atoms = positions.len() as f64;
    let shell_pairs: u64 = rdf.histogram.iter().sum();
    assert_eq!(
        shell_pairs as f64,
        n_atoms * (6.0 + 12.0 + 8.0 + 6.0 + 24.0 + 24.0 + 12.0) / 2.0
    );
    let cumulative = [
        (2.0, 0.0),
        (2.5, 6.0),
        (3.0, 18.0),
        (3.5, 26.0),
        // The bin [4, 4.5) straddles a cutoff of 4 Å and is left out.
        (4.0, 26.0),
        (4.5, 56.0),
        (5.0, 80.0),
        (6.0, 92.0),
    ];
    for (r_cut, neighbours) in cumulative {
        let coordination = rdf.coordination_number(r_cut);
        assert!(
            (coordination - neighbours).abs() < 1e-12,
            "CN({r_cut}) = {coordination}, expected {neighbours}"
        );
    }
    // g(r) of the first shell: 648 pairs over N(N-1)/2 v_shell / V ideal ones.
    let g = rdf.gofr();
    let shell = ball_volume(2.5) - ball_volume(2.0);
    let ideal_pairs = n_atoms * (n_atoms - 1.0) / 2.0 * shell / box_len.powi(3);
    assert!(
        (g[4] - 648.0 / ideal_pairs).abs() < 1e-9,
        "g(2 Å) = {}, expected {}",
        g[4],
        648.0 / ideal_pairs
    );
    // Empty shells between the neighbour shells carry no pair.
    assert_eq!(g[0], 0.0);
    assert_eq!(g[3], 0.0);
}

#[test]
fn coordination_number_excludes_the_bin_that_straddles_the_cutoff() {
    let rdf = ideal_gas_histogram(300, 20.0, 5, 11);
    let on_edge = rdf.coordination_number(4.0);
    // 4.5 Å lies inside the bin [4, 5): only the bins below 4 Å are complete.
    assert_eq!(rdf.coordination_number(4.5), on_edge);
    assert_eq!(rdf.coordination_number(4.999), on_edge);
    // A cutoff on the next bin edge picks the bin up.
    assert!(rdf.coordination_number(5.0) > on_edge);
    // The coordination number grows monotonically with the cutoff.
    let mut previous = 0.0;
    for r_cut in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] {
        let value = rdf.coordination_number(r_cut);
        assert!(value >= previous, "CN({r_cut}) = {value} < {previous}");
        previous = value;
    }
}

#[test]
fn coordination_number_of_a_single_pair() {
    // One pair at 2.5 Å in every frame: both atoms have exactly one neighbour.
    let mut rdf = RdfHistogram::new(5.0, 5, 1000.0, 2);
    let frame = [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]];
    for _ in 0..7 {
        rdf.accumulate(&frame, 100.0);
    }
    assert_eq!(rdf.coordination_number(3.0), 1.0);
    assert_eq!(rdf.coordination_number(2.0), 0.0);
    // g(r) in the shell [2, 3) Å of volume 4π/3 (27 − 8) Å^3 is V / v_shell.
    let expected = 1000.0 / (4.0 / 3.0 * PI * (27.0 - 8.0));
    let g = rdf.gofr();
    assert!((g[2] - expected).abs() < 1e-12, "g = {}", g[2]);
    for (k, &value) in g.iter().enumerate() {
        if k != 2 {
            assert_eq!(value, 0.0, "bin {k}");
        }
    }
}

#[test]
fn degenerate_systems_have_zero_g_and_zero_coordination() {
    // No atom, one atom (no pair) and no frame: nothing to normalise.  The former
    // `n_atoms * (n_atoms - 1)` underflowed for n_atoms = 0.
    let mut empty = RdfHistogram::new(10.0, 10, 1000.0, 0);
    empty.accumulate(&[], 20.0);
    assert_eq!(empty.n_frames, 1);
    assert_eq!(empty.coordination_number(5.0), 0.0);
    assert!(empty.gofr().iter().all(|&g| g == 0.0));

    let mut single = RdfHistogram::new(10.0, 10, 1000.0, 1);
    single.accumulate(&[[1.0, 2.0, 3.0]], 20.0);
    assert_eq!(single.coordination_number(5.0), 0.0);
    assert!(single.gofr().iter().all(|&g| g == 0.0));

    let no_frames = RdfHistogram::new(10.0, 10, 1000.0, 50);
    assert_eq!(no_frames.coordination_number(5.0), 0.0);
    assert!(no_frames.gofr().iter().all(|&g| g == 0.0));

    // Without a volume the density is undefined: g(r) is zero.
    let mut no_volume = RdfHistogram::new(10.0, 10, 0.0, 2);
    no_volume.accumulate(&[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]], 20.0);
    assert!(no_volume.gofr().iter().all(|&g| g == 0.0));
}

#[test]
fn reset_restarts_the_accumulation() {
    let mut rdf = ideal_gas_histogram(50, 20.0, 3, 5);
    assert_eq!(rdf.n_frames, 3);
    assert!(rdf.coordination_number(10.0) > 0.0);
    rdf.reset();
    assert_eq!(rdf.n_frames, 0);
    assert_eq!(rdf.coordination_number(10.0), 0.0);
    assert!(rdf.gofr().iter().all(|&g| g == 0.0));
}
