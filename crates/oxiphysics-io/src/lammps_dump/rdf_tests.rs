// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Tests for the radial distribution function: normalisation, periodic boundaries
//! and agreement with a brute-force reference.

use super::{LammpsDumpFrame, LammpsDumpReader, radial_distribution_function};
use std::f64::consts::PI;
use std::io::Cursor;

/// Small deterministic xorshift generator (the crate has no `rand` dependency).
struct XorShift(u64);

impl XorShift {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform value in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn frame_in_box(
    positions: Vec<[f64; 3]>,
    types: Vec<u32>,
    bounds: [[f64; 2]; 3],
    tilt: Option<[f64; 3]>,
) -> LammpsDumpFrame {
    let n = positions.len();
    LammpsDumpFrame {
        timestep: 0,
        n_atoms: n,
        box_bounds: bounds,
        atom_ids: (1..=n as u32).collect(),
        types,
        positions,
        velocities: None,
        tilt_factors: tilt,
        custom_properties: Vec::new(),
    }
}

fn cubic_frame(l: f64, positions: Vec<[f64; 3]>, types: Vec<u32>) -> LammpsDumpFrame {
    frame_in_box(positions, types, [[0.0, l]; 3], None)
}

fn shell_volume(inner: f64, outer: f64) -> f64 {
    4.0 * PI / 3.0 * (outer.powi(3) - inner.powi(3))
}

fn assert_rel(got: f64, want: f64) {
    assert!(
        (got - want).abs() <= 1e-12 * want.abs().max(1.0),
        "got {got}, want {want}"
    );
}

#[test]
fn uniform_random_gas_has_unit_g() {
    let mut rng = XorShift(0x2545_F491_4F6C_DD1D);
    let positions: Vec<[f64; 3]> = (0..1000)
        .map(|_| [rng.unit() * 10.0, rng.unit() * 10.0, rng.unit() * 10.0])
        .collect();
    let frame = cubic_frame(10.0, positions, vec![1; 1000]);
    let (r, g) = radial_distribution_function(&frame, 0, 0, 5.0, 5);
    assert_eq!(r.len(), 5);
    for (ri, gi) in r.iter().zip(&g) {
        assert!((gi - 1.0).abs() < 0.1, "g({ri}) = {gi}");
    }
}

#[test]
fn simple_cubic_lattice_has_exact_neighbour_shells() {
    // 4 x 4 x 4 lattice, spacing 1, periodic: 6 neighbours at 1, 12 at sqrt(2),
    // 8 at sqrt(3).  Without periodic images the surface atoms would miss most.
    let n = 4;
    let mut positions = Vec::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                positions.push([i as f64, j as f64, k as f64]);
            }
        }
    }
    let count = positions.len() as f64;
    let frame = cubic_frame(4.0, positions, vec![1; 64]);
    let n_bins = 20;
    let (r, g) = radial_distribution_function(&frame, 0, 0, 2.0, n_bins);
    let dr = 0.1;
    let volume = 64.0;
    let expected = |shell_count: f64, bin: usize| -> f64 {
        let hist = count * shell_count;
        hist * volume
            / (count * (count - 1.0) * shell_volume(bin as f64 * dr, (bin + 1) as f64 * dr))
    };
    for (bin, gi) in g.iter().enumerate() {
        let want = match bin {
            10 => expected(6.0, bin),
            14 => expected(12.0, bin),
            17 => expected(8.0, bin),
            _ => 0.0,
        };
        assert_rel(*gi, want);
    }
    assert_rel(r[10], 1.05);
}

#[test]
fn pair_across_a_box_face_uses_the_minimum_image() {
    // Atoms at x = 0.10 and x = 9.85 in a box of 10 are only 0.25 apart.
    let frame = cubic_frame(10.0, vec![[0.10, 5.0, 5.0], [9.85, 5.0, 5.0]], vec![1, 1]);
    let (_, g) = radial_distribution_function(&frame, 0, 0, 5.0, 50);
    let dr = 0.1;
    // Two ordered pairs fall in bin 2 = [0.2, 0.3).
    let want = 2.0 * 1000.0 / (2.0 * shell_volume(2.0 * dr, 3.0 * dr));
    for (bin, gi) in g.iter().enumerate() {
        if bin == 2 {
            assert_rel(*gi, want);
        } else {
            assert_eq!(*gi, 0.0, "bin {bin}");
        }
    }
}

#[test]
fn pair_across_a_box_corner_uses_the_minimum_image() {
    let frame = cubic_frame(
        10.0,
        vec![[0.05, 0.05, 0.05], [9.95, 9.95, 9.95]],
        vec![1, 1],
    );
    // Separation sqrt(3) * 0.1 = 0.1732.
    let (_, g) = radial_distribution_function(&frame, 0, 0, 5.0, 50);
    assert!(g[1] > 0.0 && g.iter().enumerate().all(|(i, v)| i == 1 || *v == 0.0));
}

#[test]
fn two_types_are_normalised_by_the_product_of_their_counts() {
    // Two A atoms with two B neighbours each at distance 1; N_A = 2, N_B = 4.
    let positions = vec![
        [1.0, 1.0, 1.0],
        [5.0, 5.0, 5.0],
        [1.0, 1.0, 2.0],
        [1.0, 2.0, 1.0],
        [5.0, 5.0, 6.0],
        [5.0, 4.0, 5.0],
    ];
    let frame = cubic_frame(10.0, positions, vec![1, 1, 2, 2, 2, 2]);
    let (_, g_ab) = radial_distribution_function(&frame, 1, 2, 2.0, 20);
    let want = 4.0 * 1000.0 / (2.0 * 4.0 * shell_volume(1.0, 1.1));
    for (bin, gi) in g_ab.iter().enumerate() {
        if bin == 10 {
            assert_rel(*gi, want);
        } else {
            assert_eq!(*gi, 0.0, "bin {bin}");
        }
    }
    // g_AB(r) = g_BA(r): the same pairs, counted from the other side.
    let (_, g_ba) = radial_distribution_function(&frame, 2, 1, 2.0, 20);
    for (a, b) in g_ab.iter().zip(&g_ba) {
        assert_rel(*a, *b);
    }
}

#[test]
fn self_pairs_are_excluded_when_both_selections_overlap() {
    // Wildcards and an explicit type mean the same set of atoms.
    let mut rng = XorShift(99);
    let positions: Vec<[f64; 3]> = (0..60)
        .map(|_| [rng.unit() * 8.0, rng.unit() * 8.0, rng.unit() * 8.0])
        .collect();
    let frame = cubic_frame(8.0, positions, vec![3; 60]);
    let (_, g_wild) = radial_distribution_function(&frame, 0, 0, 3.0, 12);
    let (_, g_typed) = radial_distribution_function(&frame, 3, 3, 3.0, 12);
    let (_, g_half) = radial_distribution_function(&frame, 3, 0, 3.0, 12);
    assert_eq!(g_wild, g_typed);
    assert_eq!(g_wild, g_half);
}

#[test]
fn wildcard_centre_with_one_neighbour_type_removes_only_real_self_pairs() {
    // Centres: all 6 atoms.  Neighbours: the two type-2 atoms.  Ordered pairs:
    // 4 type-1 centres * 2 + 2 type-2 centres * 1 = 10 = 6 * 2 - 2.
    let positions = vec![
        [1.0, 1.0, 1.0],
        [3.0, 1.0, 1.0],
        [5.0, 1.0, 1.0],
        [7.0, 1.0, 1.0],
        [1.0, 3.0, 1.0],
        [1.0, 3.0, 2.0],
    ];
    let frame = cubic_frame(10.0, positions, vec![1, 1, 1, 1, 2, 2]);
    let (_, g) = radial_distribution_function(&frame, 0, 2, 4.0, 40);
    // Distances below 4 (bin width 0.1): atoms 4 and 5 are 1.0 apart in both
    // directions; atom 0 is 2.0 and 2.236 from them, atom 1 is 2.828 and 3.0 away.
    // Atoms 2 and 3 are 4.47 and 4.58 away, beyond the cutoff.
    let pair_density = 10.0 / 1000.0;
    let dr = 0.1;
    for (bin, gi) in g.iter().enumerate() {
        let pairs = match bin {
            10 => 2.0,
            20 | 22 | 28 | 30 => 1.0,
            _ => 0.0,
        };
        let want = pairs / (pair_density * shell_volume(bin as f64 * dr, (bin + 1) as f64 * dr));
        assert_rel(*gi, want);
    }
}

#[test]
fn cutoff_is_clamped_to_half_the_smallest_perpendicular_width() {
    let frame = cubic_frame(10.0, vec![[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], vec![1, 1]);
    let (r, g) = radial_distribution_function(&frame, 0, 0, 100.0, 10);
    assert_eq!(r.len(), 10);
    assert_eq!(g.len(), 10);
    // Bins span [0, 5]: centres 0.25, 0.75, ...
    assert_rel(r[0], 0.25);
    assert_rel(r[9], 4.75);
    // The one pair at distance 1 is in bin [1.0, 1.5).
    assert!(g[2] > 0.0);
    assert!(g.iter().all(|v| v.is_finite()));
}

#[test]
fn degenerate_inputs_give_zero_g_without_panicking() {
    let atoms = vec![[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]];
    let frame = cubic_frame(10.0, atoms.clone(), vec![1, 1]);
    let (r, g) = radial_distribution_function(&frame, 0, 0, 5.0, 0);
    assert!(r.is_empty() && g.is_empty());
    for r_max in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let (r, g) = radial_distribution_function(&frame, 0, 0, r_max, 8);
        assert_eq!(r.len(), 8);
        if r_max.is_infinite() {
            // An infinite cutoff is clamped to half the box.
            assert!(g.iter().any(|v| *v > 0.0));
        } else {
            assert!(g.iter().all(|v| *v == 0.0), "r_max {r_max}");
        }
    }
    // A box without volume.
    let flat = frame_in_box(atoms.clone(), vec![1, 1], [[0.0, 0.0]; 3], None);
    let (_, g) = radial_distribution_function(&flat, 0, 0, 5.0, 8);
    assert!(g.iter().all(|v| *v == 0.0));
    // A selection without atoms, a single atom, and an empty frame.
    let (_, g) = radial_distribution_function(&frame, 7, 0, 5.0, 8);
    assert!(g.iter().all(|v| *v == 0.0));
    let single = cubic_frame(10.0, vec![[1.0; 3]], vec![1]);
    let (_, g) = radial_distribution_function(&single, 0, 0, 5.0, 8);
    assert!(g.iter().all(|v| *v == 0.0));
    let empty = cubic_frame(10.0, Vec::new(), Vec::new());
    let (_, g) = radial_distribution_function(&empty, 0, 0, 5.0, 8);
    assert!(g.iter().all(|v| *v == 0.0));
}

#[test]
fn inconsistent_frames_and_non_finite_atoms_are_tolerated() {
    let mut frame = cubic_frame(
        10.0,
        vec![[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [3.0, 1.0, 1.0]],
        vec![1, 1, 1],
    );
    let (_, clean) = radial_distribution_function(&frame, 0, 0, 4.0, 8);
    // More types than positions: the extra entries are ignored.
    frame.types.push(1);
    frame.n_atoms = 99;
    let (_, with_extra_type) = radial_distribution_function(&frame, 0, 0, 4.0, 8);
    assert_eq!(clean, with_extra_type);
    // An atom with no position is dropped from both the counts and the pairs.
    frame.positions.push([f64::NAN, 0.0, 0.0]);
    let (_, with_nan) = radial_distribution_function(&frame, 0, 0, 4.0, 8);
    assert_eq!(clean, with_nan);
    // Fewer types than positions.
    frame.types.truncate(2);
    let (_, short) = radial_distribution_function(&frame, 0, 0, 4.0, 8);
    assert!(short.iter().all(|v| v.is_finite()));
}

#[test]
fn tiny_cutoff_does_not_exhaust_memory() {
    let mut rng = XorShift(7);
    let positions: Vec<[f64; 3]> = (0..200)
        .map(|_| [rng.unit() * 100.0, rng.unit() * 100.0, rng.unit() * 100.0])
        .collect();
    let frame = cubic_frame(100.0, positions, vec![1; 200]);
    let (_, g) = radial_distribution_function(&frame, 0, 0, 1e-6, 4);
    assert!(g.iter().all(|v| *v == 0.0));
}

/// Brute-force g(r): the shortest image over a block of lattice translations.
fn reference_g(
    frame: &LammpsDumpFrame,
    cell: [[f64; 3]; 3],
    type_a: u32,
    type_b: u32,
    r_max: f64,
    n_bins: usize,
) -> Vec<f64> {
    let volume = {
        let [a, b, c] = cell;
        a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
    };
    let in_a = |t: u32| type_a == 0 || t == type_a;
    let in_b = |t: u32| type_b == 0 || t == type_b;
    let n = frame.positions.len();
    let n_a = frame.types.iter().filter(|&&t| in_a(t)).count() as f64;
    let n_b = frame.types.iter().filter(|&&t| in_b(t)).count() as f64;
    let n_ab = frame.types.iter().filter(|&&t| in_a(t) && in_b(t)).count() as f64;
    let mut hist = vec![0u64; n_bins];
    for i in 0..n {
        for j in 0..n {
            if i == j || !in_a(frame.types[i]) || !in_b(frame.types[j]) {
                continue;
            }
            let d = [
                frame.positions[j][0] - frame.positions[i][0],
                frame.positions[j][1] - frame.positions[i][1],
                frame.positions[j][2] - frame.positions[i][2],
            ];
            let mut best = f64::INFINITY;
            for na in -3i32..=3 {
                for nb in -3i32..=3 {
                    for nc in -3i32..=3 {
                        let v = [
                            d[0] + f64::from(na) * cell[0][0]
                                + f64::from(nb) * cell[1][0]
                                + f64::from(nc) * cell[2][0],
                            d[1] + f64::from(na) * cell[0][1]
                                + f64::from(nb) * cell[1][1]
                                + f64::from(nc) * cell[2][1],
                            d[2] + f64::from(na) * cell[0][2]
                                + f64::from(nb) * cell[1][2]
                                + f64::from(nc) * cell[2][2],
                        ];
                        best = best.min(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
                    }
                }
            }
            let r = best.sqrt();
            if r < r_max {
                let bin = (r * n_bins as f64 / r_max) as usize;
                if bin < n_bins {
                    hist[bin] += 1;
                }
            }
        }
    }
    let dr = r_max / n_bins as f64;
    hist.iter()
        .enumerate()
        .map(|(i, &c)| {
            let shell = shell_volume(i as f64 * dr, (i + 1) as f64 * dr);
            c as f64 / ((n_a * n_b - n_ab) / volume * shell)
        })
        .collect()
}

/// Random atoms in a box given by its true cell; some atoms are displaced by
/// lattice vectors so they lie outside the primary cell.
fn random_atoms(
    rng: &mut XorShift,
    n: usize,
    cell: [[f64; 3]; 3],
    origin: [f64; 3],
) -> (Vec<[f64; 3]>, Vec<u32>) {
    let mut positions = Vec::new();
    let mut types = Vec::new();
    for idx in 0..n {
        let s = [rng.unit(), rng.unit(), rng.unit()];
        let shift = if idx.is_multiple_of(5) {
            [
                (rng.unit() * 3.0).floor() - 1.0,
                (rng.unit() * 3.0).floor() - 1.0,
                (rng.unit() * 3.0).floor() - 1.0,
            ]
        } else {
            [0.0; 3]
        };
        let frac = [s[0] + shift[0], s[1] + shift[1], s[2] + shift[2]];
        positions.push([
            origin[0] + frac[0] * cell[0][0] + frac[1] * cell[1][0] + frac[2] * cell[2][0],
            origin[1] + frac[0] * cell[0][1] + frac[1] * cell[1][1] + frac[2] * cell[2][1],
            origin[2] + frac[0] * cell[0][2] + frac[1] * cell[1][2] + frac[2] * cell[2][2],
        ]);
        types.push(1 + (idx % 3) as u32);
    }
    (positions, types)
}

#[test]
fn matches_the_brute_force_reference_in_orthogonal_boxes() {
    let mut rng = XorShift(0xDEAD_BEEF_1234_5678);
    let cell = [[10.0, 0.0, 0.0], [0.0, 12.0, 0.0], [0.0, 0.0, 9.0]];
    let origin = [-1.0, 2.0, 0.5];
    let (positions, types) = random_atoms(&mut rng, 300, cell, origin);
    let bounds = [[-1.0, 9.0], [2.0, 14.0], [0.5, 9.5]];
    let frame = frame_in_box(positions, types, bounds, None);
    // r_max 4.0 and the maximum 4.5 give 2 or 3 cells along z.
    for &(ta, tb, r_max) in &[(0u32, 0u32, 4.0), (1, 2, 4.5), (2, 2, 3.0), (0, 3, 4.5)] {
        let (_, g) = radial_distribution_function(&frame, ta, tb, r_max, 25);
        let want = reference_g(&frame, cell, ta, tb, r_max, 25);
        for (bin, (got, exp)) in g.iter().zip(&want).enumerate() {
            assert_rel(*got, *exp);
            assert!(got.is_finite(), "bin {bin}");
        }
        assert!(g.iter().any(|v| *v > 0.0));
    }
}

#[test]
fn matches_the_brute_force_reference_in_triclinic_boxes() {
    let mut rng = XorShift(0x0123_4567_89AB_CDEF);
    // True box: lx = 10, ly = 12, lz = 9, xy = 3, xz = -2, yz = 1.5.
    let (xy, xz, yz) = (3.0, -2.0, 1.5);
    let cell = [[10.0, 0.0, 0.0], [xy, 12.0, 0.0], [xz, yz, 9.0]];
    let origin = [1.0, -2.0, 0.5];
    let (positions, types) = random_atoms(&mut rng, 250, cell, origin);
    // LAMMPS bounding box of the cell.
    let min_x = 0.0_f64.min(xy).min(xz).min(xy + xz);
    let max_x = 0.0_f64.max(xy).max(xz).max(xy + xz);
    let bounds = [
        [origin[0] + min_x, origin[0] + 10.0 + max_x],
        [
            origin[1] + 0.0_f64.min(yz),
            origin[1] + 12.0 + 0.0_f64.max(yz),
        ],
        [origin[2], origin[2] + 9.0],
    ];
    let frame = frame_in_box(positions, types, bounds, Some([xy, xz, yz]));
    for &(ta, tb, r_max) in &[(0u32, 0u32, 3.5), (1, 1, 4.0), (2, 3, 1.5)] {
        let (_, g) = radial_distribution_function(&frame, ta, tb, r_max, 20);
        let want = reference_g(&frame, cell, ta, tb, r_max, 20);
        for (got, exp) in g.iter().zip(&want) {
            assert_rel(*got, *exp);
        }
        assert!(g.iter().any(|v| *v > 0.0), "types ({ta}, {tb})");
    }
}

#[test]
fn tilted_cell_pair_uses_the_triclinic_minimum_image() {
    // Box 10 x 10 x 10 with xy = 3: LAMMPS writes the bounding box x in [0, 13].
    let text = "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS xy xz yz pp pp pp\n0.0 13.0 3.0\n0.0 10.0 0.0\n0.0 10.0 0.0\nITEM: ATOMS id type x y z\n1 1 1.0 0.5 5.0\n2 1 8.5 9.5 5.0\n";
    let frame = LammpsDumpReader::read_frame(Cursor::new(text.as_bytes()))
        .expect("a frame")
        .unwrap();
    assert_eq!(frame.tilt_factors, Some([3.0, 0.0, 0.0]));
    // Direct difference (7.5, 9, 0); the nearest image is (4.5, -1, 0), 4.6098 away.
    // A box without tilt would see (-2.5, -1, 0) at 2.69.
    let (r, g) = radial_distribution_function(&frame, 0, 0, 4.7, 47);
    let dr = r[1] - r[0];
    let want = 2.0 * 1000.0 / (2.0 * shell_volume(46.0 * dr, 47.0 * dr));
    for (bin, gi) in g.iter().enumerate() {
        if bin == 46 {
            assert_rel(*gi, want);
        } else {
            assert_eq!(*gi, 0.0, "bin {bin}");
        }
    }
}
