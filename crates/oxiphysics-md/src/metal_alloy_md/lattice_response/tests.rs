// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the lattice response of cubic EAM crystals.
//!
//! The reference values do not use [`AlloySystem`] at all: the energy per atom of
//! a homogeneously deformed *infinite* lattice is summed directly over the
//! neighbours of the atom at the origin (no periodic box, no minimum image), and
//! the elastic constants, pressure and bulk modulus are formed from those
//! energies with the textbook formulas.  The old implementation returned
//! `C12 + C11 / 2` for `C12`, built an FCC crystal for BCC iron, and sheared
//! only the atoms of an orthorhombic periodic box, so its periodic images were
//! not sheared.

use super::*;
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

type Matrix = [[f64; 3]; 3];

const IDENTITY: Matrix = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

fn uniaxial(epsilon: f64) -> Matrix {
    [[1.0 + epsilon, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn biaxial(epsilon: f64) -> Matrix {
    [
        [1.0 + epsilon, 0.0, 0.0],
        [0.0, 1.0 + epsilon, 0.0],
        [0.0, 0.0, 1.0],
    ]
}

fn simple_shear(gamma: f64) -> Matrix {
    [[1.0, gamma, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
}

fn basis(structure: CrystalStructure) -> &'static [[f64; 3]] {
    match structure {
        CrystalStructure::FCC => &[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ],
        CrystalStructure::BCC => &[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        CrystalStructure::HCP => &[],
    }
}

fn atoms_in_cell(params: &EamParams) -> f64 {
    basis(params.structure).len() as f64
}

/// Energy per atom (eV) of the infinite crystal deformed by `deformation`
/// (`r -> F r`), summed over the neighbours of the atom at the origin with the
/// same cutoff conventions as [`AlloySystem`].
fn reference_energy_per_atom(params: &EamParams, lattice_a: f64, deformation: Matrix) -> f64 {
    let reach = (params.cutoff / (0.9 * lattice_a)).ceil() as i32 + 2;
    let (mut pair_sum, mut density) = (0.0, 0.0);
    for ix in -reach..=reach {
        for iy in -reach..=reach {
            for iz in -reach..=reach {
                for fraction in basis(params.structure) {
                    let p = [
                        (ix as f64 + fraction[0]) * lattice_a,
                        (iy as f64 + fraction[1]) * lattice_a,
                        (iz as f64 + fraction[2]) * lattice_a,
                    ];
                    let q = [
                        deformation[0][0] * p[0]
                            + deformation[0][1] * p[1]
                            + deformation[0][2] * p[2],
                        deformation[1][0] * p[0]
                            + deformation[1][1] * p[1]
                            + deformation[1][2] * p[2],
                        deformation[2][0] * p[0]
                            + deformation[2][1] * p[1]
                            + deformation[2][2] * p[2],
                    ];
                    let r = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
                    if r > 0.1 && r < params.cutoff {
                        pair_sum += params.pair_potential(r);
                        density += params.electron_density(r);
                    }
                }
            }
        }
    }
    0.5 * pair_sum + params.embedding_energy(density)
}

/// `(C11, C12, C44)` in GPa from the reference energies with strain amplitude
/// `epsilon`.
fn reference_constants(params: &EamParams, epsilon: f64) -> (f64, f64, f64) {
    let a = params.lattice_a;
    let volume = a.powi(3) / atoms_in_cell(params);
    let energy = |deformation: Matrix| reference_energy_per_atom(params, a, deformation);
    let e0 = energy(IDENTITY);
    let c11 = (energy(uniaxial(epsilon)) + energy(uniaxial(-epsilon)) - 2.0 * e0)
        / (volume * epsilon * epsilon);
    let c11_plus_c12 = (energy(biaxial(epsilon)) + energy(biaxial(-epsilon)) - 2.0 * e0)
        / (2.0 * volume * epsilon * epsilon);
    let gamma = 2.0 * epsilon;
    let c44 = (energy(simple_shear(gamma)) + energy(simple_shear(-gamma)) - 2.0 * e0)
        / (volume * gamma * gamma);
    (
        c11 * EV_PER_ANGSTROM3_TO_GPA,
        (c11_plus_c12 - c11) * EV_PER_ANGSTROM3_TO_GPA,
        c44 * EV_PER_ANGSTROM3_TO_GPA,
    )
}

/// `-dE/dV` in GPa at the lattice constant `lattice_a`.
fn reference_pressure(params: &EamParams, lattice_a: f64) -> f64 {
    let h = 1.0e-5 * lattice_a;
    let de_da = (reference_energy_per_atom(params, lattice_a + h, IDENTITY)
        - reference_energy_per_atom(params, lattice_a - h, IDENTITY))
        / (2.0 * h);
    let dv_da = 3.0 * lattice_a * lattice_a / atoms_in_cell(params);
    -de_da / dv_da * EV_PER_ANGSTROM3_TO_GPA
}

/// `V d²E/dV²` in GPa and `-dE/dV` in GPa of the energy per atom from a
/// three-point scan in the volume, `V = a^3 / atoms` at `a (1 +- h)`.
fn reference_volume_curvature(params: &EamParams, h: f64) -> (f64, f64) {
    let a = params.lattice_a;
    let v0 = a.powi(3) / atoms_in_cell(params);
    let sample = |scale: f64| {
        (
            v0 * scale.powi(3),
            reference_energy_per_atom(params, a * scale, IDENTITY),
        )
    };
    let (v_minus, e_minus) = sample(1.0 - h);
    let (_, e_zero) = sample(1.0);
    let (v_plus, e_plus) = sample(1.0 + h);
    let (h1, h2) = (v0 - v_minus, v_plus - v0);
    let first = -h2 / (h1 * (h1 + h2)) * e_minus
        + (h2 - h1) / (h1 * h2) * e_zero
        + h1 / (h2 * (h1 + h2)) * e_plus;
    let second =
        2.0 * (e_minus / (h1 * (h1 + h2)) - e_zero / (h1 * h2) + e_plus / (h2 * (h1 + h2)));
    (
        v0 * second * EV_PER_ANGSTROM3_TO_GPA,
        -first * EV_PER_ANGSTROM3_TO_GPA,
    )
}

/// The preset with its repulsion amplitude chosen so that `lattice_a` is the
/// zero-pressure lattice constant of the potential.
///
/// The tabulated presets are under a tensile pressure of 29 GPa (aluminium) to
/// 124 GPa (nickel) at their `lattice_a`, so their strain-energy curvatures are
/// negative.
fn zero_pressure(preset: EamParams) -> EamParams {
    preset
        .calibrated_to_zero_pressure()
        .expect("the presets are valid cubic crystals")
}

fn presets() -> [EamParams; 4] {
    [
        EamParams::copper(),
        EamParams::nickel(),
        EamParams::aluminum(),
        EamParams::iron(),
    ]
}

fn relative_difference(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(1.0e-12)
}

// ---------------------------------------------------------------------------
// SimBox with a tilt
// ---------------------------------------------------------------------------

fn lattice_vector(sim_box: &SimBox, n: [f64; 3]) -> [f64; 3] {
    [
        n[0] * sim_box.lengths[0] + n[1] * sim_box.xy,
        n[1] * sim_box.lengths[1],
        n[2] * sim_box.lengths[2],
    ]
}

fn norm(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// The periodic images `v + i a + j b + k c` for `i, j, k` in `-3..=3`.
fn nearest_image_length(sim_box: &SimBox, v: [f64; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for i in -3..=3 {
        for j in -3..=3 {
            for k in -3..=3 {
                let shift = lattice_vector(sim_box, [i as f64, j as f64, k as f64]);
                best = best.min(norm([v[0] + shift[0], v[1] + shift[1], v[2] + shift[2]]));
            }
        }
    }
    best
}

/// Whether `v - w` is an integer combination of the box vectors.
fn differs_by_a_lattice_vector(sim_box: &SimBox, v: [f64; 3], w: [f64; 3]) -> bool {
    let d = [v[0] - w[0], v[1] - w[1], v[2] - w[2]];
    let c = (d[2] / sim_box.lengths[2]).round();
    let b = (d[1] / sim_box.lengths[1]).round();
    let a = ((d[0] - b * sim_box.xy) / sim_box.lengths[0]).round();
    let shift = lattice_vector(sim_box, [a, b, c]);
    (0..3).all(|dim| (d[dim] - shift[dim]).abs() < 1.0e-9)
}

#[test]
fn tilted_box_min_image_returns_the_nearest_periodic_image() {
    let sim_box = SimBox::new(4.0, 5.0, 6.0).with_xy_tilt(1.2);
    let mut rng = StdRng::seed_from_u64(5);
    let mut checked = 0;
    for _ in 0..4000 {
        let dr = [
            8.0 * rng.random::<f64>() - 4.0,
            10.0 * rng.random::<f64>() - 5.0,
            12.0 * rng.random::<f64>() - 6.0,
        ];
        let image = sim_box.min_image(dr);
        assert!(differs_by_a_lattice_vector(&sim_box, dr, image));
        // Below half of the shortest edge the nearest image is unique and must
        // be the one that is returned.
        let best = nearest_image_length(&sim_box, dr);
        if best < 2.0 {
            checked += 1;
            assert!(
                (norm(image) - best).abs() < 1.0e-12,
                "{dr:?}: image {image:?} has length {} but the nearest image has {best}",
                norm(image)
            );
        }
    }
    assert!(checked > 300, "only {checked} displacements were checked");
}

#[test]
fn tilted_box_min_image_follows_the_tilt_across_the_y_boundary() {
    let sim_box = SimBox::new(4.0, 5.0, 6.0).with_xy_tilt(1.0);
    // Crossing the y boundary downwards drags x along by -xy: (3.5, 4.0)
    // becomes (3.5 - 1.0, 4.0 - 5.0) = (2.5, -1.0) and then x wraps to -1.5.
    let image = sim_box.min_image([3.5, 4.0, 0.0]);
    assert!((image[0] + 1.5).abs() < 1e-12, "{image:?}");
    assert!((image[1] + 1.0).abs() < 1e-12, "{image:?}");
    assert_eq!(image[2], 0.0);
    // Without the tilt the same displacement stays in the first x image.
    let straight = SimBox::new(4.0, 5.0, 6.0).min_image([3.5, 4.0, 0.0]);
    assert!((straight[0] + 0.5).abs() < 1e-12, "{straight:?}");
    assert!((straight[1] + 1.0).abs() < 1e-12, "{straight:?}");
}

#[test]
fn tilted_box_wrap_puts_positions_into_the_cell() {
    let sim_box = SimBox::new(4.0, 5.0, 6.0).with_xy_tilt(1.2);
    let mut rng = StdRng::seed_from_u64(6);
    for _ in 0..2000 {
        let position = [
            40.0 * rng.random::<f64>() - 20.0,
            40.0 * rng.random::<f64>() - 20.0,
            40.0 * rng.random::<f64>() - 20.0,
        ];
        let wrapped = sim_box.wrap(position);
        let sy = wrapped[1] / 5.0;
        let sx = (wrapped[0] - 1.2 * sy) / 4.0;
        let sz = wrapped[2] / 6.0;
        for (name, s) in [("x", sx), ("y", sy), ("z", sz)] {
            assert!(
                (-1.0e-12..1.0 + 1.0e-12).contains(&s),
                "fractional {name} = {s} of {position:?}"
            );
        }
        assert!(differs_by_a_lattice_vector(&sim_box, position, wrapped));
    }
}

#[test]
fn orthorhombic_box_keeps_the_original_image_conventions() {
    let sim_box = SimBox::new(3.0, 4.0, 5.0);
    assert_eq!(sim_box.xy, 0.0);
    assert_eq!(sim_box.volume(), 60.0);
    // Half a period is left alone, like the original incremental loops did.
    assert_eq!(sim_box.min_image([1.5, -2.0, 2.5]), [1.5, -2.0, 2.5]);
    assert_eq!(sim_box.min_image([1.6, 2.1, -2.6]), [-1.4, -1.9, 2.4]);
    assert_eq!(sim_box.min_image([7.5, 0.0, 0.0])[0], 1.5);
    assert_eq!(sim_box.wrap([3.0, -0.5, 5.0]), [0.0, 3.5, 0.0]);
    assert_eq!(sim_box.wrap([-6.5, 9.0, 11.0]), [2.5, 1.0, 1.0]);
    assert_eq!(SimBox::cubic(2.0).lengths, [2.0; 3]);
}

#[test]
fn wrap_never_returns_the_upper_edge_of_the_cell() {
    let sim_box = SimBox::cubic(3.6);
    // A coordinate a rounding error below zero is the lattice point itself.
    assert_eq!(sim_box.wrap([-1.0e-17, 3.6, -3.6]), [0.0, 0.0, 0.0]);
    let tilted = SimBox::cubic(3.6).with_xy_tilt(0.2);
    let wrapped = tilted.wrap([0.5, -1.0e-17, 1.0]);
    assert_eq!(wrapped[1], 0.0);
    assert!(
        (wrapped[0] - 0.5).abs() < 1e-15 && wrapped[2] == 1.0,
        "{wrapped:?}"
    );
}

#[test]
fn degenerate_boxes_and_coordinates_do_not_loop_forever() {
    // The original loops never terminated for a zero period or a non-finite value.
    let zero = SimBox::new(0.0, f64::NAN, -1.0);
    assert_eq!(zero.min_image([1.0, 2.0, 3.0]), [1.0, 2.0, 3.0]);
    assert_eq!(zero.wrap([1.0, 2.0, 3.0]), [1.0, 2.0, 3.0]);
    let sim_box = SimBox::cubic(2.0);
    let image = sim_box.min_image([f64::INFINITY, f64::NAN, 1.0]);
    assert!(image[0].is_infinite() && image[1].is_nan() && image[2] == 1.0);
    let wrapped = sim_box.wrap([f64::NEG_INFINITY, 3.0, f64::NAN]);
    assert!(wrapped[0].is_infinite() && wrapped[1] == 1.0 && wrapped[2].is_nan());
}

// ---------------------------------------------------------------------------
// Supercell sizing and parameter validation
// ---------------------------------------------------------------------------

#[test]
fn supercell_is_large_enough_for_the_minimum_image_convention() {
    for params in presets() {
        let n = minimum_supercell_cells(&params).expect("valid preset");
        let shrunk_edge = n as f64 * params.lattice_a * (1.0 - MAX_LATTICE_SHIFT);
        assert!(
            shrunk_edge > 2.0 * params.cutoff,
            "{}: {n} cells give {shrunk_edge} Å vs 2 rc = {}",
            params.label,
            2.0 * params.cutoff
        );
        // One cell fewer would violate the condition.
        let smaller = (n - 1) as f64 * params.lattice_a * (1.0 - MAX_LATTICE_SHIFT);
        assert!(smaller <= 2.0 * params.cutoff, "{}", params.label);
    }
    assert_eq!(minimum_supercell_cells(&EamParams::copper()).unwrap(), 4);
    assert_eq!(minimum_supercell_cells(&EamParams::iron()).unwrap(), 4);
}

#[test]
fn hcp_and_invalid_parameters_are_rejected() {
    let mut hcp = EamParams::copper();
    hcp.structure = CrystalStructure::HCP;
    for result in [
        elastic_constants_eam(&hcp).map(|_| ()),
        lattice_pressure_eam(&hcp).map(|_| ()),
        thermal_expansion_coefficient(&hcp, 300.0, 600.0, 4).map(|_| ()),
        minimum_supercell_cells(&hcp).map(|_| ()),
    ] {
        assert!(
            matches!(&result, Err(MdError::General(message)) if message.contains("HCP")),
            "{result:?}"
        );
    }
    let broken = |edit: fn(&mut EamParams)| {
        let mut params = EamParams::copper();
        edit(&mut params);
        params
    };
    for params in [
        broken(|p| p.lattice_a = f64::NAN),
        broken(|p| p.lattice_a = 0.0),
        broken(|p| p.lattice_a = -3.6),
        broken(|p| p.cutoff = 0.0),
        broken(|p| p.cutoff = f64::INFINITY),
        broken(|p| p.re = -1.0),
        broken(|p| p.cutoff = 500.0),
    ] {
        assert!(elastic_constants_eam(&params).is_err(), "{params:?}");
        assert!(lattice_pressure_eam(&params).is_err(), "{params:?}");
        assert!(minimum_supercell_cells(&params).is_err(), "{params:?}");
    }
}

// ---------------------------------------------------------------------------
// Sheared supercell and the independent lattice sum
// ---------------------------------------------------------------------------

fn deformation_matrix(deformation: Deformation) -> Matrix {
    [
        [deformation.xx, deformation.xy, 0.0],
        [0.0, deformation.yy, 0.0],
        [0.0, 0.0, deformation.zz],
    ]
}

#[test]
fn undeformed_supercell_energy_equals_the_lattice_sum_for_fcc_and_bcc() {
    for params in presets() {
        let n = minimum_supercell_cells(&params).unwrap();
        let system = crystal_system(&params, params.lattice_a, n, Deformation::IDENTITY);
        assert_eq!(
            system.num_atoms(),
            n.pow(3) * atoms_in_cell(&params) as usize,
            "{}",
            params.label
        );
        let reference = reference_energy_per_atom(&params, params.lattice_a, IDENTITY);
        assert!(
            relative_difference(energy_per_atom(&system), reference) < 1e-12,
            "{}: {} vs {reference}",
            params.label,
            energy_per_atom(&system)
        );
    }
}

#[test]
fn sheared_supercell_energy_equals_the_infinite_lattice_sum() {
    for params in [EamParams::copper(), EamParams::iron()] {
        let n = minimum_supercell_cells(&params).unwrap();
        for gamma in [0.004, 0.02, -0.03, 0.05] {
            let deformation = Deformation::simple_shear(gamma);
            let system = crystal_system(&params, params.lattice_a, n, deformation);
            let reference = reference_energy_per_atom(
                &params,
                params.lattice_a,
                deformation_matrix(deformation),
            );
            assert!(
                relative_difference(energy_per_atom(&system), reference) < 1e-11,
                "{} gamma {gamma}: {} vs {reference}",
                params.label,
                energy_per_atom(&system)
            );
        }
    }
}

#[test]
fn shearing_only_the_atoms_of_an_orthorhombic_box_breaks_the_periodic_images() {
    // What the original implementation did: the atoms are moved by F but the box
    // stays orthorhombic, so atoms near the y boundary see unsheared images.
    let params = zero_pressure(EamParams::copper());
    let a = params.lattice_a;
    let n = minimum_supercell_cells(&params).unwrap();
    let edge = n as f64 * a;
    let untilted_energy = |gamma: f64| {
        let mut system = AlloySystem::new(
            vec![params.clone()],
            SimBox::cubic(edge),
            MixingRule::Geometric,
            0.001,
        );
        for position in generate_crystal_positions(a, n, params.structure) {
            let sheared = [position[0] + gamma * position[1], position[1], position[2]];
            system.add_atom(system.sim_box.wrap(sheared), 0);
        }
        energy_per_atom(&system)
    };
    let tilted_energy = |gamma: f64| {
        energy_per_atom(&crystal_system(
            &params,
            a,
            n,
            Deformation::simple_shear(gamma),
        ))
    };
    let reference_energy = |gamma: f64| reference_energy_per_atom(&params, a, simple_shear(gamma));
    // Energy of a large shear: the tilted cell reproduces the lattice sum, the
    // orthorhombic one does not.
    let gamma = 0.05;
    assert!(relative_difference(tilted_energy(gamma), reference_energy(gamma)) < 1e-11);
    assert!(
        relative_difference(untilted_energy(gamma), reference_energy(gamma)) > 1e-5,
        "orthorhombic box {} eV vs lattice sum {} eV",
        untilted_energy(gamma),
        reference_energy(gamma)
    );
    // The effect on C44, which is the second difference at small shear.
    let volume = a.powi(3) / atoms_in_cell(&params);
    let small = 2.0 * ELASTIC_STRAIN;
    let c44_from = |energy: &dyn Fn(f64) -> f64| {
        (energy(small) + energy(-small) - 2.0 * energy(0.0)) / (volume * small * small)
            * EV_PER_ANGSTROM3_TO_GPA
    };
    let correct = c44_from(&reference_energy);
    assert!(relative_difference(c44_from(&tilted_energy), correct) < 1e-6);
    let wrong = c44_from(&untilted_energy);
    assert!(
        (wrong - correct).abs() > 0.05 * correct,
        "C44 from the orthorhombic box {wrong} GPa vs correct {correct} GPa"
    );
}

#[test]
fn strain_energy_is_the_same_along_x_y_and_z() {
    for params in [EamParams::copper(), EamParams::iron()] {
        let n = minimum_supercell_cells(&params).unwrap();
        let epsilon = 0.01;
        let energy = |deformation: Deformation| {
            crystal_system(&params, params.lattice_a, n, deformation).potential_energy()
        };
        let along_x = energy(Deformation {
            xx: 1.0 + epsilon,
            ..Deformation::IDENTITY
        });
        let along_y = energy(Deformation {
            yy: 1.0 + epsilon,
            ..Deformation::IDENTITY
        });
        let along_z = energy(Deformation {
            zz: 1.0 + epsilon,
            ..Deformation::IDENTITY
        });
        assert!(
            relative_difference(along_x, along_y) < 1e-10,
            "{}",
            params.label
        );
        assert!(
            relative_difference(along_x, along_z) < 1e-10,
            "{}",
            params.label
        );
        // Sign of the strain does matter: a genuine, non-degenerate comparison.
        let compressed = energy(Deformation {
            xx: 1.0 - epsilon,
            ..Deformation::IDENTITY
        });
        assert!(relative_difference(along_x, compressed) > 1e-4);
    }
}

// ---------------------------------------------------------------------------
// Elastic constants
// ---------------------------------------------------------------------------

#[test]
fn elastic_constants_equal_the_independent_lattice_sum_derivation() {
    for preset in presets() {
        let calibrated = zero_pressure(preset.clone());
        for params in [preset, calibrated] {
            let (c11, c12, c44) = elastic_constants_eam(&params).expect("cubic crystal");
            let (r11, r12, r44) = reference_constants(&params, ELASTIC_STRAIN);
            for (name, value, reference) in
                [("C11", c11, r11), ("C12", c12, r12), ("C44", c44, r44)]
            {
                assert!(
                    (value - reference).abs() <= 1.0e-6 * reference.abs().max(1.0),
                    "{} {name}: {value} GPa vs reference {reference} GPa",
                    params.label
                );
            }
        }
    }
}

#[test]
fn c12_comes_from_the_biaxial_energy_minus_c11() {
    let params = zero_pressure(EamParams::copper());
    let (c11, c12, _) = elastic_constants_eam(&params).unwrap();
    let a = params.lattice_a;
    let volume = a.powi(3) / atoms_in_cell(&params);
    let epsilon = ELASTIC_STRAIN;
    let energy = |deformation: Matrix| reference_energy_per_atom(&params, a, deformation);
    let biaxial_curvature = (energy(biaxial(epsilon)) + energy(biaxial(-epsilon))
        - 2.0 * energy(IDENTITY))
        / (volume * epsilon * epsilon);
    let c11_plus_c12 = biaxial_curvature / 2.0 * EV_PER_ANGSTROM3_TO_GPA;
    assert!((c11 + c12 - c11_plus_c12).abs() < 1e-5 * c11_plus_c12.abs());
    // The original code returned (C11 + C12) - C11 / 2: off by C11 / 2.
    let old_c12 = c11_plus_c12 - c11 / 2.0;
    assert!(
        (old_c12 - c12).abs() > 0.4 * c11,
        "C12 = {c12}, old formula {old_c12}, C11 = {c11}"
    );
}

#[test]
fn zero_pressure_potentials_are_mechanically_stable() {
    for preset in presets() {
        let params = zero_pressure(preset);
        assert!(
            lattice_pressure_eam(&params).unwrap().abs() < 1.0e-3,
            "{}",
            params.label
        );
        let (c11, c12, c44) = elastic_constants_eam(&params).unwrap();
        // Born stability of a cubic crystal; C11 > C12 > 0 and C44 > 0 for a
        // metal with a positive bulk modulus.
        assert!(
            c11 > c12 && c12 > 0.0 && c44 > 0.0,
            "{}: {c11} {c12} {c44}",
            params.label
        );
        assert!(c11 + 2.0 * c12 > 0.0);
    }
}

#[test]
fn bulk_modulus_of_the_strain_energy_matches_the_volume_scan() {
    for preset in presets() {
        // At zero pressure (C11 + 2 C12) / 3 = V d²E/dV².
        let params = zero_pressure(preset.clone());
        let (c11, c12, _) = elastic_constants_eam(&params).unwrap();
        let (bulk_from_volume, pressure) = reference_volume_curvature(&params, 5.0e-3);
        assert!(pressure.abs() < 0.2, "{} pressure {pressure}", params.label);
        let bulk_from_strain = (c11 + 2.0 * c12) / 3.0;
        assert!(
            (bulk_from_strain - bulk_from_volume).abs() < 0.02 * bulk_from_volume,
            "{}: strain {bulk_from_strain} GPa vs volume {bulk_from_volume} GPa",
            params.label
        );
        // Under a residual pressure P the identity gains -2 P / 3.
        let (c11, c12, _) = elastic_constants_eam(&preset).unwrap();
        let (curvature, pressure) = reference_volume_curvature(&preset, 5.0e-3);
        let bulk_from_strain = (c11 + 2.0 * c12) / 3.0;
        let expected = curvature - 2.0 * pressure / 3.0;
        // The two terms nearly cancel for some presets, so the tolerance refers to
        // the larger of them.
        let scale = curvature.abs().max((2.0 * pressure / 3.0).abs());
        assert!(
            (bulk_from_strain - expected).abs() < 0.02 * scale,
            "{} preset: strain {bulk_from_strain} GPa vs V E'' - 2P/3 = {expected} GPa",
            preset.label
        );
    }
}

#[test]
fn tabulated_presets_are_far_from_equilibrium_at_their_lattice_constant() {
    for preset in presets() {
        let pressure = lattice_pressure_eam(&preset).unwrap();
        let reference = reference_pressure(&preset, preset.lattice_a);
        assert!(
            (pressure - reference).abs() < 1e-5 * reference.abs(),
            "{}: {pressure} GPa vs {reference} GPa",
            preset.label
        );
        assert!(pressure < -10.0, "{}: {pressure} GPa", preset.label);
    }
}

#[test]
fn calibration_to_zero_pressure_is_confirmed_by_the_independent_lattice_sum() {
    for preset in presets() {
        let calibrated = preset.calibrated_to_zero_pressure().unwrap();
        // Only the repulsion amplitude changes, and it has to grow: the tabulated
        // presets are under tension.
        assert!(calibrated.a_coeff > preset.a_coeff, "{}", preset.label);
        assert_eq!(calibrated.lattice_a, preset.lattice_a);
        assert_eq!(calibrated.b_coeff, preset.b_coeff);
        assert_eq!(calibrated.f0, preset.f0);
        assert_eq!(calibrated.structure, preset.structure);
        let reference = reference_pressure(&calibrated, calibrated.lattice_a);
        assert!(
            reference.abs() < 1.0e-6,
            "{}: {reference} GPa",
            preset.label
        );
        let pressure = lattice_pressure_eam(&calibrated).unwrap();
        assert!(pressure.abs() < 1.0e-6, "{}: {pressure} GPa", preset.label);
    }
}

#[test]
fn calibration_rejects_unusable_parameters() {
    let mut hcp = EamParams::copper();
    hcp.structure = CrystalStructure::HCP;
    assert!(hcp.calibrated_to_zero_pressure().is_err());
    // A repulsive many-body term (negative `f0`) already pushes the lattice apart;
    // balancing it would need a negative pair-repulsion amplitude.
    let mut repulsive_embedding = EamParams::copper();
    repulsive_embedding.f0 = -50.0;
    let error = repulsive_embedding
        .calibrated_to_zero_pressure()
        .expect_err("a negative amplitude would be needed");
    assert!(
        matches!(&error, MdError::General(message) if message.contains("positive repulsion")),
        "{error:?}"
    );
}

#[test]
fn a_neighbour_shell_on_the_cutoff_is_reported_instead_of_returned() {
    // The fifth FCC shell lies at a sqrt(5/2); a cutoff just beyond it makes a
    // part of the shell leave the interaction range under +epsilon, so the energy
    // jumps and the second difference scales as 1 / epsilon².
    let mut params = zero_pressure(EamParams::copper());
    params.cutoff = params.lattice_a * 2.5_f64.sqrt() * (1.0 + 5.0e-5);
    let error = elastic_constants_eam(&params).expect_err("energy is discontinuous");
    assert!(
        matches!(&error, MdError::General(message) if message.contains("not smooth")),
        "{error:?}"
    );
}

// ---------------------------------------------------------------------------
// Quasi-harmonic Einstein model
// ---------------------------------------------------------------------------

#[test]
fn einstein_force_constant_matches_the_curvature_of_the_energy() {
    for preset in [EamParams::copper(), EamParams::iron()] {
        let params = zero_pressure(preset);
        let n = minimum_supercell_cells(&params).unwrap();
        for lattice_a in [params.lattice_a, 1.01 * params.lattice_a] {
            let k = einstein_force_constant(&params, lattice_a, n).unwrap();
            let displaced_energy = |shift: f64| {
                let mut system = crystal_system(&params, lattice_a, n, Deformation::IDENTITY);
                system.atoms[0].position[0] += shift;
                system.potential_energy()
            };
            let u = 2.0e-3;
            let curvature = (displaced_energy(u) + displaced_energy(-u)
                - 2.0 * displaced_energy(0.0))
                / (u * u);
            assert!(k > 1.0, "{}: k = {k} eV/Å²", params.label);
            assert!(
                (k - curvature).abs() < 1.0e-5 * k,
                "{} a = {lattice_a}: force {k} vs energy {curvature}",
                params.label
            );
        }
    }
}

#[test]
fn einstein_force_constant_is_isotropic_for_a_cubic_site() {
    let params = zero_pressure(EamParams::iron());
    let n = minimum_supercell_cells(&params).unwrap();
    let mut system = crystal_system(&params, params.lattice_a, n, Deformation::IDENTITY);
    let force_constant = |system: &mut AlloySystem, dim: usize| {
        let origin = system.atoms[0].position[dim];
        let mut force = |shift: f64| {
            system.atoms[0].position[dim] = origin + shift;
            system.compute_forces();
            system.atoms[0].force[dim]
        };
        let (plus, minus) = (force(1.0e-3), force(-1.0e-3));
        system.atoms[0].position[dim] = origin;
        -(plus - minus) / 2.0e-3
    };
    let kx = force_constant(&mut system, 0);
    let ky = force_constant(&mut system, 1);
    let kz = force_constant(&mut system, 2);
    assert!(relative_difference(kx, ky) < 1e-8 && relative_difference(kx, kz) < 1e-8);
}

#[test]
fn einstein_frequency_converts_units_and_rejects_unstable_lattices() {
    // k = 5 eV/Å², m = 63.546 amu: omega = sqrt(5 * 9648.533 / 63.546) rad/ps.
    let omega = einstein_frequency(5.0, 63.546).unwrap();
    assert!((omega - (5.0 * 9_648.533_216 / 63.546_f64).sqrt()).abs() < 1e-9);
    // Cu-like numbers give an Einstein temperature of a couple of hundred kelvin.
    let theta = HBAR_EV_PS * omega / KB_EV_PER_K;
    assert!((150.0..300.0).contains(&theta), "theta_E = {theta} K");
    assert!(einstein_frequency(0.0, 63.546).is_none());
    assert!(einstein_frequency(-1.0, 63.546).is_none());
    assert!(einstein_frequency(f64::NAN, 63.546).is_none());
}

#[test]
fn vibrational_free_energy_has_the_classical_and_zero_point_limits() {
    let kt = KB_EV_PER_K * 300.0;
    // hbar omega << kT: 3 kT ln(hbar omega / kT).
    let omega_low = 1.0e-3 * kt / HBAR_EV_PS;
    let classical = 3.0 * kt * (1.0e-3_f64).ln();
    assert!((vibrational_free_energy(omega_low, 300.0) - classical).abs() < 1e-5 * classical.abs());
    // hbar omega >> kT: the zero-point energy 3 hbar omega / 2.
    let omega_high = 60.0 * kt / HBAR_EV_PS;
    let zero_point = 1.5 * HBAR_EV_PS * omega_high;
    assert!((vibrational_free_energy(omega_high, 300.0) - zero_point).abs() < 1e-12 * zero_point);
    // dF/d omega = 3 hbar (1/2 + n_B) with the Bose occupation n_B.
    let omega = 30.0;
    let h = 1.0e-5;
    let slope = (vibrational_free_energy(omega + h, 400.0)
        - vibrational_free_energy(omega - h, 400.0))
        / (2.0 * h);
    let x = HBAR_EV_PS * omega / (KB_EV_PER_K * 400.0);
    let expected = 3.0 * HBAR_EV_PS * (0.5 + 1.0 / x.exp_m1());
    assert!(
        (slope - expected).abs() < 1e-8 * expected,
        "{slope} vs {expected}"
    );
}

#[test]
fn quasi_harmonic_lattice_constant_is_a_free_energy_minimum() {
    let params = zero_pressure(EamParams::nickel());
    let n = minimum_supercell_cells(&params).unwrap();
    for temperature in [100.0, 300.0, 800.0] {
        let a = quasi_harmonic_lattice_constant(&params, temperature, n).unwrap();
        let free_energy = |lattice_a: f64| free_energy_per_atom(&params, lattice_a, n, temperature);
        let delta = 2.0e-4 * a;
        let (below, here, above) = (
            free_energy(a - delta),
            free_energy(a),
            free_energy(a + delta),
        );
        assert!(
            here < below && here < above,
            "T = {temperature}: F({}) = {below}, F({a}) = {here}, F({}) = {above}",
            a - delta,
            a + delta
        );
        let slope = (above - below) / (2.0 * delta);
        assert!(
            slope.abs() < 2.0e-5,
            "T = {temperature}: dF/da = {slope} eV/Å"
        );
        // The static minimum is a0 itself, so the lattice only expands.
        assert!(a > params.lattice_a, "T = {temperature}: a = {a}");
    }
}

#[test]
fn thermal_expansion_matches_the_derivative_of_the_free_energy() {
    // da/dT = -(d²F/da dT) / (d²F/da²) from finite differences of the same free
    // energy, evaluated at the temperature in the middle of the interval.
    for preset in [EamParams::nickel(), EamParams::aluminum()] {
        let params = zero_pressure(preset);
        let n = minimum_supercell_cells(&params).unwrap();
        let (t_low, t_high) = (300.0, 310.0);
        let alpha = thermal_expansion_coefficient(&params, t_low, t_high, n).unwrap();
        let t_mid = 0.5 * (t_low + t_high);
        let a_mid = quasi_harmonic_lattice_constant(&params, t_mid, n).unwrap();
        let free_energy = |a: f64, t: f64| free_energy_per_atom(&params, a, n, t);
        let (h, dt) = (1.0e-3 * a_mid, 1.0);
        let f_aa = (free_energy(a_mid + h, t_mid) - 2.0 * free_energy(a_mid, t_mid)
            + free_energy(a_mid - h, t_mid))
            / (h * h);
        let f_at = (free_energy(a_mid + h, t_mid + dt)
            - free_energy(a_mid - h, t_mid + dt)
            - free_energy(a_mid + h, t_mid - dt)
            + free_energy(a_mid - h, t_mid - dt))
            / (4.0 * h * dt);
        let expected = -f_at / f_aa / a_mid;
        assert!(expected > 0.0, "{}", params.label);
        assert!(
            (alpha - expected).abs() < 0.01 * expected,
            "{}: secant {alpha} vs derivative {expected}",
            params.label
        );
    }
}

#[test]
fn lattice_expands_monotonically_with_temperature() {
    let params = zero_pressure(EamParams::copper());
    let n = minimum_supercell_cells(&params).unwrap();
    let mut previous = 0.0;
    for temperature in [10.0, 100.0, 300.0, 600.0, 900.0] {
        let a = quasi_harmonic_lattice_constant(&params, temperature, n).unwrap();
        assert!(a > previous, "a({temperature} K) = {a} after {previous}");
        previous = a;
    }
}

#[test]
fn thermal_expansion_of_the_zero_pressure_metals_is_physical() {
    for preset in presets() {
        let params = zero_pressure(preset);
        let n = minimum_supercell_cells(&params).unwrap();
        let alpha = thermal_expansion_coefficient(&params, 300.0, 600.0, n).unwrap();
        // Measured values of Cu, Ni, Al and Fe lie between 1e-5 and 2.4e-5 per kelvin.
        assert!(
            (1.0e-5..3.0e-5).contains(&alpha),
            "{}: alpha = {alpha} 1/K",
            params.label
        );
    }
}

#[test]
fn thermal_expansion_depends_on_the_temperatures() {
    // The original function ignored t_low and t_high and returned a constant
    // from a hard-coded Grüneisen parameter.
    let params = zero_pressure(EamParams::aluminum());
    let n = minimum_supercell_cells(&params).unwrap();
    let cold = thermal_expansion_coefficient(&params, 20.0, 60.0, n).unwrap();
    let warm = thermal_expansion_coefficient(&params, 300.0, 340.0, n).unwrap();
    let hot = thermal_expansion_coefficient(&params, 800.0, 840.0, n).unwrap();
    // The vibrational heat capacity freezes out below the Einstein temperature.
    assert!(cold < 0.5 * warm, "cold {cold}, warm {warm}");
    assert!(warm > 0.0 && hot >= warm * 0.99, "warm {warm}, hot {hot}");
    assert!(cold >= 0.0);
}

#[test]
fn invalid_temperatures_are_errors() {
    let params = zero_pressure(EamParams::copper());
    let n = minimum_supercell_cells(&params).unwrap();
    for (low, high) in [
        (300.0, 300.0),
        (600.0, 300.0),
        (0.0, 300.0),
        (-10.0, 300.0),
        (f64::NAN, 300.0),
        (300.0, f64::INFINITY),
    ] {
        let result = thermal_expansion_coefficient(&params, low, high, n);
        assert!(
            matches!(&result, Err(MdError::General(message)) if message.contains("t_low")),
            "({low}, {high}): {result:?}"
        );
    }
    for temperature in [0.0, -1.0, f64::NAN] {
        assert!(quasi_harmonic_lattice_constant(&params, temperature, n).is_err());
    }
    let mut massless = params.clone();
    massless.mass = 0.0;
    assert!(quasi_harmonic_lattice_constant(&massless, 300.0, n).is_err());
}

#[test]
fn too_small_or_too_large_a_supercell_is_an_error() {
    let params = zero_pressure(EamParams::copper());
    // Two cells (7.2 Å) cannot hold a cutoff sphere of radius 5.5 Å.
    let error = thermal_expansion_coefficient(&params, 300.0, 600.0, 2).unwrap_err();
    assert!(
        matches!(&error, MdError::General(message) if message.contains("too small")),
        "{error:?}"
    );
    let error = quasi_harmonic_lattice_constant(&params, 300.0, 40).unwrap_err();
    assert!(
        matches!(&error, MdError::General(message) if message.contains("maximum")),
        "{error:?}"
    );
}

#[test]
fn unstable_tabulated_potential_is_reported_with_its_pressure() {
    // The tabulated Cu parameters are under about -100 GPa at their lattice
    // constant, where the Einstein force constant is even negative.
    let error = thermal_expansion_coefficient(&EamParams::copper(), 300.0, 600.0, 4)
        .expect_err("no stable lattice at lattice_a");
    match &error {
        MdError::General(message) => {
            assert!(message.contains("dynamically unstable"), "{message}");
            assert!(message.contains("GPa"), "{message}");
        }
        other => panic!("unexpected error {other:?}"),
    }
}

#[test]
fn lattice_constant_far_from_equilibrium_has_no_free_energy_minimum_in_the_window() {
    // A lattice constant 5 % above the zero-pressure value is still dynamically
    // stable, but the minimum lies 4.8 % below it, outside the +-3 % window.
    let mut params = zero_pressure(EamParams::copper());
    params.lattice_a *= 1.05;
    let n = minimum_supercell_cells(&params).unwrap();
    let error = quasi_harmonic_lattice_constant(&params, 300.0, n).unwrap_err();
    match &error {
        MdError::General(message) => {
            assert!(message.contains("no minimum"), "{message}");
            assert!(message.contains("GPa"), "{message}");
        }
        other => panic!("unexpected error {other:?}"),
    }
}

#[test]
fn temperatures_beyond_the_search_window_are_an_error() {
    // At a million kelvin the vibrational pressure pushes the lattice far beyond +3 %.
    let params = zero_pressure(EamParams::copper());
    let n = minimum_supercell_cells(&params).unwrap();
    let error = quasi_harmonic_lattice_constant(&params, 1.0e6, n).unwrap_err();
    assert!(
        matches!(&error, MdError::General(message) if message.contains("no minimum")),
        "{error:?}"
    );
}

#[test]
fn nearest_interior_minimum_prefers_the_minimum_closest_to_the_centre() {
    // Two minima: index 1 (distance 4) and index 7 (distance 2).
    let values = [3.0, 1.0, 2.0, 3.0, 4.0, 3.0, 2.0, 0.5, 2.0, 5.0];
    assert_eq!(nearest_interior_minimum(&values, 5), Some(7));
    assert_eq!(nearest_interior_minimum(&values, 2), Some(1));
    // Monotonic data and minima at the ends are not interior minima.
    assert_eq!(nearest_interior_minimum(&[1.0, 2.0, 3.0, 4.0], 1), None);
    assert_eq!(nearest_interior_minimum(&[4.0, 3.0, 2.0, 1.0], 1), None);
    // Infinite neighbours are fine, an infinite candidate is not.
    assert_eq!(
        nearest_interior_minimum(&[f64::INFINITY, 1.0, 2.0], 1),
        Some(1)
    );
    assert_eq!(
        nearest_interior_minimum(&[1.0, f64::INFINITY, 2.0], 1),
        None
    );
}

#[test]
fn golden_section_finds_the_minimum_and_tolerates_infinite_regions() {
    let parabola = |x: f64| (x - 1.2345).powi(2) + 0.5;
    // A flat minimum above zero cannot be located better than sqrt(machine epsilon).
    let x = golden_section_minimum(&parabola, 0.0, 3.0, 1.0e-10);
    assert!((x - 1.2345).abs() < 1.0e-6, "{x}");
    let with_wall = |x: f64| {
        if x < 0.5 {
            f64::INFINITY
        } else {
            (x - 0.9).powi(2)
        }
    };
    let x = golden_section_minimum(&with_wall, 0.4, 2.0, 1.0e-10);
    assert!((x - 0.9).abs() < 1.0e-8, "{x}");
}

// ---------------------------------------------------------------------------
// Independent checks of the potential and of the end-to-end numbers
// ---------------------------------------------------------------------------

/// The reference tests above sum the energy with the methods of [`EamParams`],
/// so they cannot notice a slip in those methods.  These tests write the
/// documented formulas out once more and pin the final numbers to a separate
/// implementation.
#[test]
fn eam_functional_forms_match_the_documented_formulas() {
    for params in presets() {
        let EamParams {
            re,
            a_coeff,
            b_coeff,
            alpha,
            beta,
            fe,
            f0,
            cutoff,
            ..
        } = params.clone();
        let radii = [
            0.3,
            0.499,
            0.5,
            1.0,
            re,
            3.3,
            0.5 * (re + cutoff),
            cutoff - 1.0e-9,
            cutoff,
            cutoff + 1.0,
        ];
        for r in radii {
            let in_range = (0.5..cutoff).contains(&r);
            let x = r / re - 1.0;
            let repulsion = (-alpha * x).exp();
            let attraction = (-beta * x).exp();
            let (pair, pair_slope, density, density_slope) = if in_range {
                (
                    a_coeff * repulsion - b_coeff * attraction,
                    (-a_coeff * alpha * repulsion + b_coeff * beta * attraction) / re,
                    fe * attraction,
                    -fe * beta * attraction / re,
                )
            } else {
                (0.0, 0.0, 0.0, 0.0)
            };
            for (name, got, want) in [
                ("phi", params.pair_potential(r), pair),
                ("dphi/dr", params.pair_force(r), pair_slope),
                ("rho", params.electron_density(r), density),
                ("drho/dr", params.electron_density_deriv(r), density_slope),
            ] {
                assert!(
                    (got - want).abs() <= 1.0e-12 * want.abs().max(1.0e-3),
                    "{} {name}({r}) = {got}, documented {want}",
                    params.label
                );
            }
        }
        for rho in [-1.0_f64, 0.0, 1.0e-12, 0.5, 4.0, 18.6] {
            let (energy, slope) = if rho > 0.0 {
                (-f0 * rho.sqrt(), -f0 / (2.0 * rho.sqrt()))
            } else {
                (0.0, 0.0)
            };
            assert!((params.embedding_energy(rho) - energy).abs() <= 1.0e-13 * energy.abs());
            assert!(
                (params.embedding_deriv(rho) - slope).abs() <= 1.0e-12 * slope.abs(),
                "{} F'({rho})",
                params.label
            );
        }
    }
}

/// Reference numbers of a preset from a separate implementation of the
/// documented potential (a short script outside this repository): direct
/// lattice sums over an infinite crystal (no periodic box, no
/// [`AlloySystem`]), a two-point derivative of the energy for the pressure, the
/// analytic force on the displaced atom (written out from the formulas) for the
/// Einstein force constant, and a golden-section search of the free energy.
/// The numbers pin the whole chain from the potential to the thermal expansion;
/// they change only with the tabulated parameters or the model.
///
/// The Einstein force constant has to come from forces: the second difference of
/// the energy at the same displacement is limited by round-off to about 1e-9
/// relative, which makes the location of the free-energy minimum, and with it the
/// expansion coefficient, uncertain at the per-mille level.
struct IndependentValues {
    label: &'static str,
    /// Static pressure (GPa) of the tabulated preset at its lattice constant.
    stock_pressure: f64,
    /// Pair-repulsion amplitude (eV) that makes the lattice constant the
    /// zero-pressure one.
    a_coeff: f64,
    /// `C11`, `C12` and `C44` (GPa) at that amplitude.
    elastic: [f64; 3],
    /// Linear expansion coefficient (1/K) between 300 K and 600 K.
    expansion: f64,
}

const INDEPENDENT_VALUES: [IndependentValues; 4] = [
    IndependentValues {
        label: "Cu",
        stock_pressure: -97.863_158,
        a_coeff: 0.995_923_206_8,
        elastic: [148.988_756_55, 108.915_929_45, 51.088_533_32],
        expansion: 1.288_982_43e-5,
    },
    IndependentValues {
        label: "Ni",
        stock_pressure: -123.675_749,
        a_coeff: 1.135_697_458_6,
        elastic: [186.679_760_11, 137.816_681_76, 64.991_796_41],
        expansion: 1.124_142_61e-5,
    },
    IndependentValues {
        label: "Al",
        stock_pressure: -28.623_757,
        a_coeff: 0.572_763_364_8,
        elastic: [64.394_565_03, 39.583_989_04, 23.665_630_98],
        expansion: 1.816_058_87e-5,
    },
    IndependentValues {
        label: "Fe",
        stock_pressure: -78.902_692,
        a_coeff: 0.940_533_585_8,
        elastic: [134.309_675_32, 108.892_815_35, 54.096_282_31],
        expansion: 1.338_601_47e-5,
    },
];

#[test]
fn zero_pressure_presets_reproduce_an_independent_implementation() {
    for expected in &INDEPENDENT_VALUES {
        let preset = presets()
            .into_iter()
            .find(|candidate| candidate.label == expected.label)
            .expect("the table lists the presets");
        let pressure = lattice_pressure_eam(&preset).unwrap();
        assert!(
            relative_difference(pressure, expected.stock_pressure) < 1.0e-5,
            "{}: stock pressure {pressure} GPa, independent {} GPa",
            expected.label,
            expected.stock_pressure
        );
        let params = preset.calibrated_to_zero_pressure().unwrap();
        assert!(
            relative_difference(params.a_coeff, expected.a_coeff) < 1.0e-6,
            "{}: calibrated A = {} eV, independent {} eV",
            expected.label,
            params.a_coeff,
            expected.a_coeff
        );
        let (c11, c12, c44) = elastic_constants_eam(&params).unwrap();
        for (name, got, want) in [
            ("C11", c11, expected.elastic[0]),
            ("C12", c12, expected.elastic[1]),
            ("C44", c44, expected.elastic[2]),
        ] {
            assert!(
                relative_difference(got, want) < 1.0e-5,
                "{} {name} = {got} GPa, independent {want} GPa",
                expected.label
            );
        }
        let n = minimum_supercell_cells(&params).unwrap();
        let alpha = thermal_expansion_coefficient(&params, 300.0, 600.0, n).unwrap();
        assert!(
            relative_difference(alpha, expected.expansion) < 1.0e-4,
            "{}: alpha = {alpha} 1/K, independent {} 1/K",
            expected.label,
            expected.expansion
        );
    }
}

/// Positions (Å) of the neighbours of an atom of the perfect infinite crystal
/// that lie inside the cutoff.
fn lattice_neighbours(params: &EamParams, lattice_a: f64) -> Vec<[f64; 3]> {
    let reach = (params.cutoff / (0.9 * lattice_a)).ceil() as i32 + 2;
    let mut neighbours = Vec::new();
    for ix in -reach..=reach {
        for iy in -reach..=reach {
            for iz in -reach..=reach {
                for fraction in basis(params.structure) {
                    let q = [
                        (ix as f64 + fraction[0]) * lattice_a,
                        (iy as f64 + fraction[1]) * lattice_a,
                        (iz as f64 + fraction[2]) * lattice_a,
                    ];
                    let r = norm(q);
                    if r > 0.1 && r < params.cutoff {
                        neighbours.push(q);
                    }
                }
            }
        }
    }
    neighbours
}

/// Einstein force constant `d²E/du²` (eV/Å²) of the atom at the origin of the
/// infinite perfect crystal, from the analytic force on that atom displaced by
/// `±u` along x while every other atom stays on its site.
///
/// With the documented energy `E = Σ_j φ(r_j) + F(ρ_0) + Σ_j F(ρ_j)`, where
/// `r_j = |R_j - u x̂|`, `ρ_0 = Σ_j ρ(r_j)` and `ρ_j = ρ_rest - ρ(|R_j|) + ρ(r_j)`
/// (a neighbour only sees the displaced atom move), the force is
/// `F_x = Σ_j [φ'(r_j) + (F'(ρ_0) + F'(ρ_j)) ρ'(r_j)] (X_j - u) / r_j`.
fn reference_einstein_force_constant(params: &EamParams, lattice_a: f64) -> f64 {
    let neighbours = lattice_neighbours(params, lattice_a);
    let rest_density: f64 = neighbours
        .iter()
        .map(|q| params.electron_density(norm(*q)))
        .sum();
    let force_x = |shift: f64| -> f64 {
        let distances: Vec<f64> = neighbours
            .iter()
            .map(|q| norm([q[0] - shift, q[1], q[2]]))
            .collect();
        let displaced_density: f64 = distances.iter().map(|&r| params.electron_density(r)).sum();
        neighbours
            .iter()
            .zip(&distances)
            .map(|(q, &r)| {
                let neighbour_density =
                    rest_density - params.electron_density(norm(*q)) + params.electron_density(r);
                let slope = params.pair_force(r)
                    + (params.embedding_deriv(displaced_density)
                        + params.embedding_deriv(neighbour_density))
                        * params.electron_density_deriv(r);
                slope * (q[0] - shift) / r
            })
            .sum()
    };
    let u = 1.0e-3;
    -(force_x(u) - force_x(-u)) / (2.0 * u)
}

/// Helmholtz free energy per atom (eV) of the Einstein crystal from the lattice
/// sum and the reference force constant; `+inf` where the lattice is
/// dynamically unstable.
fn reference_free_energy(params: &EamParams, lattice_a: f64, temperature: f64) -> f64 {
    let force_constant = reference_einstein_force_constant(params, lattice_a);
    match einstein_frequency(force_constant, params.mass) {
        Some(omega) => {
            reference_energy_per_atom(params, lattice_a, IDENTITY)
                + vibrational_free_energy(omega, temperature)
        }
        None => f64::INFINITY,
    }
}

#[test]
fn quasi_harmonic_expansion_matches_the_displaced_atom_lattice_sum_model() {
    for preset in presets() {
        let params = zero_pressure(preset);
        let n = minimum_supercell_cells(&params).unwrap();
        let a0 = params.lattice_a;
        // The static minimum is a0 and the vibrations push the lattice outwards,
        // so the minimum of the free energy lies in [a0, 1.02 a0].
        let reference = |temperature: f64| {
            golden_section_minimum(
                &|a: f64| reference_free_energy(&params, a, temperature),
                a0,
                1.02 * a0,
                1.0e-9 * a0,
            )
        };
        let (reference_low, reference_high) = (reference(300.0), reference(600.0));
        let low = quasi_harmonic_lattice_constant(&params, 300.0, n).unwrap();
        let high = quasi_harmonic_lattice_constant(&params, 600.0, n).unwrap();
        for (temperature, got, want) in [(300.0, low, reference_low), (600.0, high, reference_high)]
        {
            assert!(
                relative_difference(got, want) < 1.0e-6,
                "{} at {temperature} K: a = {got} Å, reference {want} Å",
                params.label
            );
        }
        let alpha = thermal_expansion_coefficient(&params, 300.0, 600.0, n).unwrap();
        let expected =
            (reference_high - reference_low) / (0.5 * (reference_low + reference_high) * 300.0);
        assert!(
            relative_difference(alpha, expected) < 1.0e-4,
            "{}: alpha = {alpha} 1/K, reference {expected} 1/K",
            params.label
        );
    }
}
