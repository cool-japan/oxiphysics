// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests for the D3Q19 / D3Q27 Zou-He boundary conditions and the
//! D3Q19 moving-wall bounce-back.
//!
//! The oracles are explicit sums over the lattice velocity tables and
//! population-by-population reference formulas (Hecht & Harting 2010) rather
//! than the production code.  Most inputs are *non-equilibrium* states (an
//! equilibrium plus deterministic noise on every population): rest or
//! equilibrium inputs cannot tell mirror pairs from lattice opposites, and
//! they satisfy the mass and normal-momentum balance even with wrong
//! coefficients.  Equilibrium inputs are used for round trips, with the
//! unknown populations overwritten by junk first so that a function that does
//! nothing cannot pass.

use super::*;
use crate::lattice::{D3Q19_OPPOSITES, D3Q19_VELOCITIES, D3Q19_WEIGHTS, D3Q27_VELOCITIES};

/// Absolute tolerance for identities that hold up to rounding error.
const TOL: f64 = 1.0e-12;
/// Value written into populations that a boundary condition must overwrite.
const JUNK: f64 = 0.5;

/// Wall / inlet velocities, including components along every axis.
const WALL_VELOCITIES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [0.05, 0.0, 0.0],
    [0.0, -0.04, 0.03],
    [-0.03, 0.02, -0.06],
];

/// Tangential velocities for the pressure conditions.
const TANGENTIAL_VELOCITIES: [[f64; 2]; 3] = [[0.0, 0.0], [0.03, -0.02], [-0.05, 0.04]];

/// SplitMix64 generator, so the tests need no RNG crate and stay deterministic.
struct SplitMix64(u64);

impl SplitMix64 {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform sample in `[-1, 1)`.
    fn next_signed(&mut self) -> f64 {
        let unit = (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64;
        2.0 * unit - 1.0
    }
}

/// Add uniform noise of the given amplitude to every population.
fn with_noise<const Q: usize>(mut f: [f64; Q], amplitude: f64, seed: u64) -> [f64; Q] {
    let mut rng = SplitMix64(seed);
    for fi in &mut f {
        *fi += amplitude * rng.next_signed();
    }
    f
}

/// Velocity table and per-face unknown count of a lattice under test.
struct LatticeTables<const Q: usize> {
    velocities: &'static [[i32; 3]; Q],
    unknown_per_face: usize,
}

const D19: LatticeTables<19> = LatticeTables {
    velocities: &D3Q19_VELOCITIES,
    unknown_per_face: 5,
};
const D27: LatticeTables<27> = LatticeTables {
    velocities: &D3Q27_VELOCITIES,
    unknown_per_face: 9,
};

/// Inward unit normal of a face, written out independently of `ZouHeFace`.
fn inward_normal(face: ZouHeFace) -> [i32; 3] {
    match face {
        ZouHeFace::Left => [1, 0, 0],
        ZouHeFace::Right => [-1, 0, 0],
        ZouHeFace::Bottom => [0, 1, 0],
        ZouHeFace::Top => [0, -1, 0],
        ZouHeFace::Back => [0, 0, 1],
        ZouHeFace::Front => [0, 0, -1],
    }
}

/// Populations that stream into the domain through `face`, i.e. those whose
/// velocity has a positive component along the inward normal.
fn unknown_set<const Q: usize>(lattice: &LatticeTables<Q>, face: ZouHeFace) -> Vec<usize> {
    let n = inward_normal(face);
    (0..Q)
        .filter(|&i| {
            let c = lattice.velocities[i];
            c[0] * n[0] + c[1] * n[1] + c[2] * n[2] > 0
        })
        .collect()
}

/// Overwrite the listed populations with `value`.
fn overwrite<const Q: usize>(f: &mut [f64; Q], indices: &[usize], value: f64) {
    for &i in indices {
        f[i] = value;
    }
}

/// Density and momentum `(Σ f, Σ f c)` by explicit summation.
fn moments<const Q: usize>(f: &[f64; Q], velocities: &[[i32; 3]; Q]) -> (f64, [f64; 3]) {
    let rho: f64 = f.iter().sum();
    let mut momentum = [0.0; 3];
    for (fi, c) in f.iter().zip(velocities.iter()) {
        for (jd, &cd) in momentum.iter_mut().zip(c.iter()) {
            *jd += fi * f64::from(cd);
        }
    }
    (rho, momentum)
}

/// A velocity vector with `value` along `axis` and zero elsewhere.
fn along_axis(axis: usize, value: f64) -> [f64; 3] {
    let mut u = [0.0; 3];
    u[axis] = value;
    u
}

/// Assert the postconditions of a boundary condition on a node: the node
/// carries mass `rho` and momentum `rho u` in all three components, exactly
/// the unknown populations of `face` were rewritten, and the known ones are
/// bit-identical to the input.
fn assert_boundary_state<const Q: usize>(
    lattice: &LatticeTables<Q>,
    name: &str,
    face: ZouHeFace,
    before: &[f64; Q],
    after: &[f64; Q],
    target: (f64, [f64; 3]),
) {
    let (rho, u) = target;
    let (rho_sum, momentum) = moments(after, lattice.velocities);
    assert!(
        (rho_sum - rho).abs() < TOL,
        "{name} {face:?}: sum f = {rho_sum}, expected {rho}"
    );
    for (d, (&j, &ud)) in momentum.iter().zip(u.iter()).enumerate() {
        assert!(
            (j - rho * ud).abs() < TOL,
            "{name} {face:?}: momentum component {d} = {j}, expected {}",
            rho * ud
        );
    }
    let unknown = unknown_set(lattice, face);
    assert_eq!(
        unknown.len(),
        lattice.unknown_per_face,
        "{name} {face:?}: unknown population count"
    );
    for (i, (&b, &a)) in before.iter().zip(after.iter()).enumerate() {
        if !unknown.contains(&i) {
            assert_eq!(a, b, "{name} {face:?}: known population {i} was modified");
        }
    }
}

type D19Condition = fn(&mut [f64; 19], f64) -> f64;
type D19Entry = (&'static str, D19Condition, ZouHeFace);

/// Public D3Q19 velocity conditions `(name, function, face)`; the argument
/// is the velocity component normal to the face, the result the density.
const D19_VELOCITY_CONDITIONS: [D19Entry; 8] = [
    (
        "zou_he_inlet_left_d3q19",
        zou_he_inlet_left_d3q19,
        ZouHeFace::Left,
    ),
    (
        "apply_velocity_outlet_d3q19",
        apply_velocity_outlet_d3q19,
        ZouHeFace::Right,
    ),
    (
        "zou_he_top_face_d3q19",
        zou_he_top_face_d3q19,
        ZouHeFace::Top,
    ),
    (
        "zou_he_bottom_face_d3q19",
        zou_he_bottom_face_d3q19,
        ZouHeFace::Bottom,
    ),
    (
        "zou_he_front_face_d3q19",
        zou_he_front_face_d3q19,
        ZouHeFace::Front,
    ),
    (
        "zou_he_back_face_d3q19",
        zou_he_back_face_d3q19,
        ZouHeFace::Back,
    ),
    (
        "zou_he_velocity_front_d3q19",
        zou_he_velocity_front_d3q19,
        ZouHeFace::Front,
    ),
    (
        "zou_he_velocity_back_d3q19",
        zou_he_velocity_back_d3q19,
        ZouHeFace::Back,
    ),
];

/// Public D3Q19 pressure conditions `(name, function, face)`; the argument
/// is the density, the result the velocity component along the positive axis
/// normal to the face.
const D19_PRESSURE_CONDITIONS: [D19Entry; 9] = [
    (
        "apply_pressure_inlet_d3q19",
        apply_pressure_inlet_d3q19,
        ZouHeFace::Left,
    ),
    (
        "zou_he_outlet_right_d3q19",
        zou_he_outlet_right_d3q19,
        ZouHeFace::Right,
    ),
    (
        "zou_he_pressure_left_d3q19",
        zou_he_pressure_left_d3q19,
        ZouHeFace::Left,
    ),
    (
        "zou_he_pressure_right_d3q19",
        zou_he_pressure_right_d3q19,
        ZouHeFace::Right,
    ),
    (
        "zou_he_pressure_outlet_left_d3q19",
        zou_he_pressure_outlet_left_d3q19,
        ZouHeFace::Left,
    ),
    (
        "zou_he_pressure_top_d3q19",
        zou_he_pressure_top_d3q19,
        ZouHeFace::Top,
    ),
    (
        "zou_he_pressure_bottom_d3q19",
        zou_he_pressure_bottom_d3q19,
        ZouHeFace::Bottom,
    ),
    (
        "zou_he_pressure_front_d3q19",
        zou_he_pressure_front_d3q19,
        ZouHeFace::Front,
    ),
    (
        "zou_he_pressure_back_d3q19",
        zou_he_pressure_back_d3q19,
        ZouHeFace::Back,
    ),
];

// ---------------------------------------------------------------------------
// Face geometry
// ---------------------------------------------------------------------------

#[test]
fn face_geometry_matches_domain_layout() {
    let expected = [
        (ZouHeFace::Left, 0, 1, [1, 2]),
        (ZouHeFace::Right, 0, -1, [1, 2]),
        (ZouHeFace::Bottom, 1, 1, [0, 2]),
        (ZouHeFace::Top, 1, -1, [0, 2]),
        (ZouHeFace::Back, 2, 1, [0, 1]),
        (ZouHeFace::Front, 2, -1, [0, 1]),
    ];
    assert_eq!(ZouHeFace::ALL.len(), expected.len());
    for (face, axis, sign, tangential) in expected {
        assert_eq!(face.axis(), axis, "{face:?} axis");
        assert_eq!(face.inward_sign(), sign, "{face:?} inward sign");
        assert_eq!(
            face.tangential_axes(),
            tangential,
            "{face:?} tangential axes"
        );
        let n = inward_normal(face);
        assert_eq!(
            n[axis], sign,
            "{face:?} normal agrees with the explicit table"
        );
    }
    for (k, a) in ZouHeFace::ALL.iter().enumerate() {
        for b in &ZouHeFace::ALL[k + 1..] {
            assert_ne!(a, b, "ZouHeFace::ALL lists a face twice");
        }
    }
}

// ---------------------------------------------------------------------------
// D3Q19: generic face functions
// ---------------------------------------------------------------------------

#[test]
fn d3q19_velocity_face_imposes_full_velocity_on_every_face() {
    for (k, face) in ZouHeFace::ALL.into_iter().enumerate() {
        for (m, u) in WALL_VELOCITIES.into_iter().enumerate() {
            let seed = (10 * k + m) as u64;
            let before = with_noise(d3q19_equilibrium(1.03, 0.01, -0.02, 0.015), 0.004, seed);
            let mut after = before;
            let rho = zou_he_velocity_face_d3q19(&mut after, face, u);
            assert_boundary_state(&D19, "velocity", face, &before, &after, (rho, u));
        }
    }
}

#[test]
fn d3q19_pressure_face_imposes_density_and_tangential_velocity_on_every_face() {
    let rho_set = 1.02;
    for (k, face) in ZouHeFace::ALL.into_iter().enumerate() {
        for (m, ut) in TANGENTIAL_VELOCITIES.into_iter().enumerate() {
            let seed = (10 * k + m) as u64 + 500;
            let before = with_noise(d3q19_equilibrium(0.98, 0.0, 0.0, 0.0), 0.004, seed);
            let mut after = before;
            let un = zou_he_pressure_face_d3q19(&mut after, face, rho_set, ut);
            let [first, second] = face.tangential_axes();
            let mut u = along_axis(face.axis(), un);
            u[first] = ut[0];
            u[second] = ut[1];
            assert_boundary_state(&D19, "pressure", face, &before, &after, (rho_set, u));
        }
    }
}

#[test]
fn d3q19_face_functions_never_read_the_unknown_populations() {
    for face in ZouHeFace::ALL {
        let unknown = unknown_set(&D19, face);
        let base = with_noise(d3q19_equilibrium(1.01, 0.02, 0.01, -0.01), 0.004, 77);
        let mut first = base;
        let mut second = base;
        overwrite(&mut first, &unknown, JUNK);
        overwrite(&mut second, &unknown, -3.0);
        let rho_first = zou_he_velocity_face_d3q19(&mut first, face, [0.03, -0.02, 0.01]);
        let rho_second = zou_he_velocity_face_d3q19(&mut second, face, [0.03, -0.02, 0.01]);
        assert_eq!(
            rho_first, rho_second,
            "{face:?}: velocity density depends on unknowns"
        );
        assert_eq!(
            first, second,
            "{face:?}: velocity result depends on unknowns"
        );

        let mut first = base;
        let mut second = base;
        overwrite(&mut first, &unknown, JUNK);
        overwrite(&mut second, &unknown, -3.0);
        let un_first = zou_he_pressure_face_d3q19(&mut first, face, 1.02, [0.01, -0.01]);
        let un_second = zou_he_pressure_face_d3q19(&mut second, face, 1.02, [0.01, -0.01]);
        assert_eq!(
            un_first, un_second,
            "{face:?}: pressure velocity depends on unknowns"
        );
        assert_eq!(
            first, second,
            "{face:?}: pressure result depends on unknowns"
        );
    }
}

#[test]
fn d3q19_velocity_face_returns_equilibrium_for_consistent_input() {
    let rho0 = 1.02;
    for face in ZouHeFace::ALL {
        for u in WALL_VELOCITIES {
            let feq = d3q19_equilibrium(rho0, u[0], u[1], u[2]);
            let mut f = feq;
            overwrite(&mut f, &unknown_set(&D19, face), JUNK);
            let rho = zou_he_velocity_face_d3q19(&mut f, face, u);
            assert!((rho - rho0).abs() < TOL, "{face:?} {u:?}: density {rho}");
            for (i, (&got, &want)) in f.iter().zip(feq.iter()).enumerate() {
                assert!(
                    (got - want).abs() < TOL,
                    "{face:?} {u:?}: f[{i}] = {got}, equilibrium {want}"
                );
            }
        }
    }
}

#[test]
fn d3q19_pressure_face_returns_equilibrium_for_consistent_input() {
    let rho0 = 1.02;
    for face in ZouHeFace::ALL {
        for u in WALL_VELOCITIES {
            let feq = d3q19_equilibrium(rho0, u[0], u[1], u[2]);
            let mut f = feq;
            overwrite(&mut f, &unknown_set(&D19, face), JUNK);
            let [first, second] = face.tangential_axes();
            let un = zou_he_pressure_face_d3q19(&mut f, face, rho0, [u[first], u[second]]);
            assert!(
                (un - u[face.axis()]).abs() < TOL,
                "{face:?} {u:?}: normal velocity {un}"
            );
            for (i, (&got, &want)) in f.iter().zip(feq.iter()).enumerate() {
                assert!(
                    (got - want).abs() < TOL,
                    "{face:?} {u:?}: f[{i}] = {got}, equilibrium {want}"
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// D3Q19: the named per-face functions
// ---------------------------------------------------------------------------

#[test]
fn d3q19_velocity_conditions_impose_prescribed_velocity_without_transverse_momentum() {
    for (k, (name, condition, face)) in D19_VELOCITY_CONDITIONS.into_iter().enumerate() {
        for (m, speed) in [0.04, -0.03, 0.0].into_iter().enumerate() {
            let mut before = with_noise(
                d3q19_equilibrium(1.02, 0.0, 0.0, 0.0),
                0.004,
                (10 * k + m) as u64 + 900,
            );
            // Disturb the edge populations that carry transverse momentum.
            before[15] += 0.01;
            before[8] -= 0.007;
            let mut after = before;
            let rho = condition(&mut after, speed);
            let u = along_axis(face.axis(), speed);
            assert_boundary_state(&D19, name, face, &before, &after, (rho, u));
        }
    }
}

#[test]
fn d3q19_pressure_conditions_impose_density_without_transverse_momentum() {
    let rho_set = 1.02;
    for (k, (name, condition, face)) in D19_PRESSURE_CONDITIONS.into_iter().enumerate() {
        let mut before = with_noise(
            d3q19_equilibrium(0.99, 0.0, 0.0, 0.0),
            0.004,
            k as u64 + 1300,
        );
        before[15] += 0.01;
        before[8] -= 0.007;
        let mut after = before;
        let un = condition(&mut after, rho_set);
        let u = along_axis(face.axis(), un);
        assert_boundary_state(&D19, name, face, &before, &after, (rho_set, u));
    }
}

#[test]
fn d3q19_velocity_conditions_return_equilibrium_for_consistent_input() {
    let rho0 = 1.02;
    for (name, condition, face) in D19_VELOCITY_CONDITIONS {
        let speed = 0.045;
        let u = along_axis(face.axis(), speed);
        let feq = d3q19_equilibrium(rho0, u[0], u[1], u[2]);
        let mut f = feq;
        overwrite(&mut f, &unknown_set(&D19, face), JUNK);
        let rho = condition(&mut f, speed);
        assert!((rho - rho0).abs() < TOL, "{name}: density {rho}");
        for (i, (&got, &want)) in f.iter().zip(feq.iter()).enumerate() {
            assert!(
                (got - want).abs() < TOL,
                "{name}: f[{i}] = {got}, equilibrium {want}"
            );
        }
    }
}

#[test]
fn d3q19_pressure_conditions_return_equilibrium_for_consistent_input() {
    let rho0 = 1.02;
    for (name, condition, face) in D19_PRESSURE_CONDITIONS {
        let speed = 0.05;
        let u = along_axis(face.axis(), speed);
        let feq = d3q19_equilibrium(rho0, u[0], u[1], u[2]);
        let mut f = feq;
        overwrite(&mut f, &unknown_set(&D19, face), JUNK);
        let un = condition(&mut f, rho0);
        assert!((un - speed).abs() < TOL, "{name}: velocity {un}");
        for (i, (&got, &want)) in f.iter().zip(feq.iter()).enumerate() {
            assert!(
                (got - want).abs() < TOL,
                "{name}: f[{i}] = {got}, equilibrium {want}"
            );
        }
    }
}

#[test]
fn d3q19_named_functions_agree_with_the_generic_face_functions() {
    let before = with_noise(d3q19_equilibrium(1.01, 0.01, 0.02, -0.01), 0.004, 4242);
    for (name, condition, face) in D19_VELOCITY_CONDITIONS {
        let mut named = before;
        let mut generic = before;
        let rho_named = condition(&mut named, 0.037);
        let rho_generic =
            zou_he_velocity_face_d3q19(&mut generic, face, along_axis(face.axis(), 0.037));
        assert_eq!(rho_named, rho_generic, "{name}");
        assert_eq!(named, generic, "{name}");
    }
    for (name, condition, face) in D19_PRESSURE_CONDITIONS {
        let mut named = before;
        let mut generic = before;
        let un_named = condition(&mut named, 1.015);
        let un_generic = zou_he_pressure_face_d3q19(&mut generic, face, 1.015, [0.0, 0.0]);
        assert_eq!(un_named, un_generic, "{name}");
        assert_eq!(named, generic, "{name}");
    }
}

// ---------------------------------------------------------------------------
// D3Q19: Hecht-Harting reference formulas written out population by population
// ---------------------------------------------------------------------------

/// Hecht & Harting (2010) velocity inlet on the `x = 0` face for the
/// velocity `u`, written out population by population with their
/// `N_y`, `N_z` (which include the edge populations 15-18 and the
/// `rho u_t / 3` term).  Returns the new populations and the density.
fn hecht_harting_left_inlet(f: &[f64; 19], u: [f64; 3]) -> ([f64; 19], f64) {
    let known_mass = f[0]
        + f[3]
        + f[4]
        + f[5]
        + f[6]
        + f[15]
        + f[16]
        + f[17]
        + f[18]
        + 2.0 * (f[2] + f[8] + f[10] + f[12] + f[14]);
    let rho = known_mass / (1.0 - u[0]);
    let n_y = 0.5 * (f[3] + f[15] + f[17] - (f[4] + f[16] + f[18])) - rho * u[1] / 3.0;
    let n_z = 0.5 * (f[5] + f[15] + f[16] - (f[6] + f[17] + f[18])) - rho * u[2] / 3.0;
    let mut out = *f;
    out[1] = f[2] + rho * u[0] / 3.0;
    out[7] = f[10] + rho * (u[0] + u[1]) / 6.0 - n_y;
    out[9] = f[8] + rho * (u[0] - u[1]) / 6.0 + n_y;
    out[11] = f[14] + rho * (u[0] + u[2]) / 6.0 - n_z;
    out[13] = f[12] + rho * (u[0] - u[2]) / 6.0 + n_z;
    (out, rho)
}

/// The same scheme mirrored to the `x = nx - 1` face with a prescribed
/// density and tangential velocity `(uy, uz)`.  Returns the new populations
/// and the velocity `ux`.
fn hecht_harting_right_pressure(f: &[f64; 19], rho: f64, ut: [f64; 2]) -> ([f64; 19], f64) {
    let known_mass = f[0]
        + f[3]
        + f[4]
        + f[5]
        + f[6]
        + f[15]
        + f[16]
        + f[17]
        + f[18]
        + 2.0 * (f[1] + f[7] + f[9] + f[11] + f[13]);
    let ux = -1.0 + known_mass / rho;
    let n_y = 0.5 * (f[3] + f[15] + f[17] - (f[4] + f[16] + f[18])) - rho * ut[0] / 3.0;
    let n_z = 0.5 * (f[5] + f[15] + f[16] - (f[6] + f[17] + f[18])) - rho * ut[1] / 3.0;
    let mut out = *f;
    out[2] = f[1] - rho * ux / 3.0;
    out[8] = f[9] + rho * (-ux + ut[0]) / 6.0 - n_y;
    out[10] = f[7] + rho * (-ux - ut[0]) / 6.0 + n_y;
    out[12] = f[13] + rho * (-ux + ut[1]) / 6.0 - n_z;
    out[14] = f[11] + rho * (-ux - ut[1]) / 6.0 + n_z;
    (out, ux)
}

/// Assert that two populations arrays agree to the test tolerance.
fn assert_populations_close(context: &str, got: &[f64; 19], want: &[f64; 19]) {
    for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
        assert!(
            (g - w).abs() < TOL,
            "{context}: f[{i}] = {g}, reference {w}"
        );
    }
}

#[test]
fn d3q19_inlet_left_matches_hecht_harting_reference() {
    let velocities = [
        [0.05, 0.0, 0.0],
        [-0.03, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        [0.08, 0.0, 0.0],
        [0.04, 0.02, -0.03],
        [-0.02, -0.04, 0.05],
    ];
    for (seed, u) in (1u64..).zip(velocities) {
        let before = with_noise(d3q19_equilibrium(1.0, 0.0, 0.02, -0.01), 0.005, seed);
        let (expected, rho_expected) = hecht_harting_left_inlet(&before, u);

        let mut generic = before;
        let rho = zou_he_velocity_face_d3q19(&mut generic, ZouHeFace::Left, u);
        assert!(
            (rho - rho_expected).abs() < TOL,
            "seed {seed}: density {rho}, reference {rho_expected}"
        );
        assert_populations_close(&format!("seed {seed} {u:?}"), &generic, &expected);

        if u[1] == 0.0 && u[2] == 0.0 {
            let mut named = before;
            let rho_named = zou_he_inlet_left_d3q19(&mut named, u[0]);
            assert!(
                (rho_named - rho_expected).abs() < TOL,
                "seed {seed}: named density {rho_named}"
            );
            assert_populations_close(&format!("seed {seed} named"), &named, &expected);
        }
    }
}

#[test]
fn d3q19_pressure_right_matches_mirrored_hecht_harting_reference() {
    let cases = [
        (1.02, [0.0, 0.0]),
        (0.97, [0.0, 0.0]),
        (1.0, [0.0, 0.0]),
        (1.01, [0.03, -0.02]),
        (0.99, [-0.04, 0.05]),
    ];
    for (seed, (rho_set, ut)) in (11u64..).zip(cases) {
        let before = with_noise(d3q19_equilibrium(1.0, 0.03, 0.02, -0.01), 0.005, seed);
        let (expected, ux_expected) = hecht_harting_right_pressure(&before, rho_set, ut);

        let mut generic = before;
        let ux = zou_he_pressure_face_d3q19(&mut generic, ZouHeFace::Right, rho_set, ut);
        assert!((ux - ux_expected).abs() < TOL, "seed {seed}: ux {ux}");
        assert_populations_close(&format!("seed {seed} {ut:?}"), &generic, &expected);

        if ut == [0.0, 0.0] {
            for condition in [zou_he_outlet_right_d3q19, zou_he_pressure_right_d3q19] {
                let mut named = before;
                let ux_named = condition(&mut named, rho_set);
                assert!(
                    (ux_named - ux_expected).abs() < TOL,
                    "seed {seed}: named ux {ux_named}"
                );
                assert_populations_close(&format!("seed {seed} named"), &named, &expected);
            }
        }
    }
}

#[test]
fn d3q19_inlet_uses_d3q19_coefficients_and_true_opposites() {
    // At rest the transverse corrections vanish, so every unknown differs from
    // its lattice opposite by exactly rho * (c . u) / 6 for edge populations and
    // rho * ux / 3 for the axial one (the D2Q9 values 2/3 and 1/12 are wrong).
    let ux = 0.04;
    let mut f = d3q19_equilibrium(1.0, 0.0, 0.0, 0.0);
    let rho = zou_he_inlet_left_d3q19(&mut f, ux);
    assert!((rho - 1.0 / (1.0 - ux)).abs() < TOL, "density {rho}");
    assert!(
        (f[1] - f[2] - rho * ux / 3.0).abs() < TOL,
        "axial population"
    );
    for i in [7usize, 9, 11, 13] {
        let opposite = D3Q19_OPPOSITES[i];
        assert!(
            (f[i] - f[opposite] - rho * ux / 6.0).abs() < TOL,
            "edge population {i} against its opposite {opposite}"
        );
    }
}

#[test]
fn d3q19_inlet_removes_transverse_momentum_caused_by_edge_populations() {
    // Raising f15 alone puts momentum +1 in y and +1 in z on the node; the
    // boundary condition must compensate it with the unknown populations.
    let ux = 0.03;
    let mut f = d3q19_equilibrium(1.0, 0.0, 0.0, 0.0);
    f[15] += 0.01;
    let rho = zou_he_inlet_left_d3q19(&mut f, ux);
    let (rho_sum, momentum) = moments(&f, &D3Q19_VELOCITIES);
    assert!((rho_sum - rho).abs() < TOL, "mass {rho_sum} against {rho}");
    assert!((momentum[0] - rho * ux).abs() < TOL, "jx = {}", momentum[0]);
    assert!(momentum[1].abs() < TOL, "jy = {}", momentum[1]);
    assert!(momentum[2].abs() < TOL, "jz = {}", momentum[2]);
}

#[test]
fn d3q19_pressure_left_and_right_do_not_triple_count_known_populations() {
    let rho0 = 1.02;
    let ux0 = 0.05;
    let feq = d3q19_equilibrium(rho0, ux0, 0.0, 0.0);

    let mut left = feq;
    overwrite(&mut left, &[1, 7, 9, 11, 13], JUNK);
    let ux_left = zou_he_pressure_left_d3q19(&mut left, rho0);
    assert!((ux_left - ux0).abs() < TOL, "left ux = {ux_left}");
    for (i, (&got, &want)) in left.iter().zip(feq.iter()).enumerate() {
        assert!(
            (got - want).abs() < TOL,
            "left f[{i}] = {got}, equilibrium {want}"
        );
    }

    let mut right = feq;
    overwrite(&mut right, &[2, 8, 10, 12, 14], JUNK);
    let ux_right = zou_he_pressure_right_d3q19(&mut right, rho0);
    assert!((ux_right - ux0).abs() < TOL, "right ux = {ux_right}");
    for (i, (&got, &want)) in right.iter().zip(feq.iter()).enumerate() {
        assert!(
            (got - want).abs() < TOL,
            "right f[{i}] = {got}, equilibrium {want}"
        );
    }
}

#[test]
fn d3q19_pressure_outlet_left_only_rewrites_the_incoming_populations() {
    let before = with_noise(d3q19_equilibrium(1.0, -0.04, 0.0, 0.0), 0.004, 31);
    let rho_out = 0.99;
    let mut f = before;
    let ux = zou_he_pressure_outlet_left_d3q19(&mut f, rho_out);
    // The known -x populations 2, 8, 10, 12, 14 used to be overwritten.
    for i in [0usize, 2, 3, 4, 5, 6, 8, 10, 12, 14, 15, 16, 17, 18] {
        assert_eq!(f[i], before[i], "known population {i} was modified");
    }
    let mut written = 0;
    for i in [1usize, 7, 9, 11, 13] {
        if f[i] != before[i] {
            written += 1;
        }
    }
    assert_eq!(written, 5, "all five incoming populations must be rebuilt");
    let (rho_sum, momentum) = moments(&f, &D3Q19_VELOCITIES);
    assert!((rho_sum - rho_out).abs() < TOL, "density {rho_sum}");
    assert!(
        (momentum[0] - rho_out * ux).abs() < TOL,
        "jx = {}",
        momentum[0]
    );
    assert!(momentum[1].abs() < TOL && momentum[2].abs() < TOL);

    let mut inlet = before;
    let ux_inlet = apply_pressure_inlet_d3q19(&mut inlet, rho_out);
    assert_eq!(ux, ux_inlet);
    assert_eq!(f, inlet, "same condition as apply_pressure_inlet_d3q19");
}

// ---------------------------------------------------------------------------
// D3Q27
// ---------------------------------------------------------------------------

#[test]
fn d3q27_velocity_face_imposes_full_velocity_on_every_face() {
    for (k, face) in ZouHeFace::ALL.into_iter().enumerate() {
        for (m, u) in WALL_VELOCITIES.into_iter().enumerate() {
            let seed = (10 * k + m) as u64 + 2000;
            let before = with_noise(d3q27_equilibrium(1.03, 0.01, -0.02, 0.015), 0.001, seed);
            let mut after = before;
            let rho = zou_he_velocity_face_d3q27(&mut after, face, u);
            assert_boundary_state(&D27, "velocity", face, &before, &after, (rho, u));
        }
    }
}

#[test]
fn d3q27_pressure_face_imposes_density_and_tangential_velocity_on_every_face() {
    let rho_set = 1.02;
    for (k, face) in ZouHeFace::ALL.into_iter().enumerate() {
        for (m, ut) in TANGENTIAL_VELOCITIES.into_iter().enumerate() {
            let seed = (10 * k + m) as u64 + 2500;
            let before = with_noise(d3q27_equilibrium(0.98, 0.0, 0.0, 0.0), 0.001, seed);
            let mut after = before;
            let un = zou_he_pressure_face_d3q27(&mut after, face, rho_set, ut);
            let [first, second] = face.tangential_axes();
            let mut u = along_axis(face.axis(), un);
            u[first] = ut[0];
            u[second] = ut[1];
            assert_boundary_state(&D27, "pressure", face, &before, &after, (rho_set, u));
        }
    }
}

#[test]
fn d3q27_face_functions_return_equilibrium_for_consistent_input() {
    let rho0 = 1.02;
    for face in ZouHeFace::ALL {
        for u in WALL_VELOCITIES {
            let feq = d3q27_equilibrium(rho0, u[0], u[1], u[2]);

            let mut f = feq;
            overwrite(&mut f, &unknown_set(&D27, face), JUNK);
            let rho = zou_he_velocity_face_d3q27(&mut f, face, u);
            assert!((rho - rho0).abs() < TOL, "{face:?} {u:?}: density {rho}");
            for (i, (&got, &want)) in f.iter().zip(feq.iter()).enumerate() {
                assert!(
                    (got - want).abs() < TOL,
                    "velocity {face:?} {u:?}: f[{i}] = {got}, equilibrium {want}"
                );
            }

            let mut f = feq;
            overwrite(&mut f, &unknown_set(&D27, face), JUNK);
            let [first, second] = face.tangential_axes();
            let un = zou_he_pressure_face_d3q27(&mut f, face, rho0, [u[first], u[second]]);
            assert!(
                (un - u[face.axis()]).abs() < TOL,
                "{face:?} {u:?}: normal velocity {un}"
            );
            for (i, (&got, &want)) in f.iter().zip(feq.iter()).enumerate() {
                assert!(
                    (got - want).abs() < TOL,
                    "pressure {face:?} {u:?}: f[{i}] = {got}, equilibrium {want}"
                );
            }
        }
    }
}

#[test]
fn d3q27_velocity_inlet_imposes_density_and_all_velocity_components() {
    for (m, u) in WALL_VELOCITIES.into_iter().enumerate() {
        let before = with_noise(d3q27_equilibrium(1.03, 0.01, -0.02, 0.015), 0.001, m as u64);
        // Raise some of the populations that carry transverse momentum.
        let mut before = before;
        before[15] += 0.004;
        before[22] -= 0.002;
        let mut after = before;
        let rho = zou_he_velocity_inlet_d3q27(&mut after, u[0], u[1], u[2]);
        assert_boundary_state(
            &D27,
            "zou_he_velocity_inlet_d3q27",
            ZouHeFace::Left,
            &before,
            &after,
            (rho, u),
        );
    }
}

#[test]
fn d3q27_velocity_inlet_is_not_a_plain_copy_of_the_opposite_populations() {
    // At rest the unknowns differ from their opposites by 6 w rho (c . u).
    use crate::lattice::{D3Q27_OPPOSITES, D3Q27_WEIGHTS};
    let ux = 0.05;
    let mut f = d3q27_equilibrium(1.0, 0.0, 0.0, 0.0);
    let rho = zou_he_velocity_inlet_d3q27(&mut f, ux, 0.0, 0.0);
    assert!((rho - 1.0 / (1.0 - ux)).abs() < TOL, "density {rho}");
    for i in unknown_set(&D27, ZouHeFace::Left) {
        let opposite = D3Q27_OPPOSITES[i];
        let expected = 6.0 * D3Q27_WEIGHTS[i] * rho * f64::from(D3Q27_VELOCITIES[i][0]) * ux;
        assert!(
            (f[i] - f[opposite] - expected).abs() < TOL,
            "population {i} against opposite {opposite}: difference {}, expected {expected}",
            f[i] - f[opposite]
        );
        assert!(expected > 0.0);
    }
}

#[test]
fn d3q27_pressure_outlet_imposes_density_without_transverse_momentum() {
    let rho_out = 1.02;
    let mut before = with_noise(d3q27_equilibrium(0.99, 0.0, 0.0, 0.0), 0.001, 99);
    before[15] += 0.004;
    before[22] -= 0.002;
    let mut after = before;
    let ux = zou_he_pressure_outlet_d3q27(&mut after, rho_out);
    assert_boundary_state(
        &D27,
        "zou_he_pressure_outlet_d3q27",
        ZouHeFace::Right,
        &before,
        &after,
        (rho_out, [ux, 0.0, 0.0]),
    );
}

#[test]
fn d3q27_named_functions_return_equilibrium_for_consistent_input() {
    let rho0 = 1.02;
    let u = [0.04, -0.02, 0.03];
    let feq = d3q27_equilibrium(rho0, u[0], u[1], u[2]);
    let mut inlet = feq;
    overwrite(&mut inlet, &[1, 7, 9, 11, 13, 19, 21, 23, 25], JUNK);
    let rho = zou_he_velocity_inlet_d3q27(&mut inlet, u[0], u[1], u[2]);
    assert!((rho - rho0).abs() < TOL, "inlet density {rho}");
    for (i, (&got, &want)) in inlet.iter().zip(feq.iter()).enumerate() {
        assert!(
            (got - want).abs() < TOL,
            "inlet f[{i}] = {got}, equilibrium {want}"
        );
    }

    let feq = d3q27_equilibrium(rho0, u[0], 0.0, 0.0);
    let mut outlet = feq;
    overwrite(&mut outlet, &[2, 8, 10, 12, 14, 20, 22, 24, 26], JUNK);
    let ux = zou_he_pressure_outlet_d3q27(&mut outlet, rho0);
    assert!((ux - u[0]).abs() < TOL, "outlet ux = {ux}");
    for (i, (&got, &want)) in outlet.iter().zip(feq.iter()).enumerate() {
        assert!(
            (got - want).abs() < TOL,
            "outlet f[{i}] = {got}, equilibrium {want}"
        );
    }
}

// ---------------------------------------------------------------------------
// D3Q19 moving-wall bounce-back
// ---------------------------------------------------------------------------

#[test]
fn moving_top_wall_reflects_with_signed_linear_correction() {
    let rho = 1.03;
    let (ux_wall, uz_wall) = (0.05, -0.03);
    let before = with_noise(d3q19_equilibrium(rho, 0.0, 0.0, 0.0), 0.004, 7);
    let mut f = before;
    bounce_back_moving_top_d3q19(&mut f, rho, ux_wall, uz_wall);
    // (population towards the wall, reflected population, c_in . u_wall)
    let links = [
        (3usize, 4usize, 0.0),
        (7, 10, ux_wall),
        (8, 9, -ux_wall),
        (15, 18, uz_wall),
        (17, 16, -uz_wall),
    ];
    for (i_in, i_out, speed) in links {
        let expected = before[i_in] - 6.0 * D3Q19_WEIGHTS[i_in] * rho * speed;
        assert!(
            (f[i_out] - expected).abs() < 1.0e-15,
            "link {i_in} -> {i_out}: {}, expected {expected}",
            f[i_out]
        );
    }
    let written: Vec<usize> = links.iter().map(|link| link.1).collect();
    for i in 0..19 {
        if !written.contains(&i) {
            assert_eq!(f[i], before[i], "population {i} must not be touched");
        }
    }
    // f17 points towards the top wall and f16 away from it.
    assert_eq!(f[17], before[17]);
    assert_ne!(f[16], before[16]);
}

#[test]
fn moving_bottom_wall_reflects_with_signed_linear_correction() {
    let rho = 0.97;
    let (ux_wall, uz_wall) = (-0.04, 0.02);
    let before = with_noise(d3q19_equilibrium(rho, 0.0, 0.0, 0.0), 0.004, 8);
    let mut f = before;
    bounce_back_moving_bottom_d3q19(&mut f, rho, ux_wall, uz_wall);
    let links = [
        (4usize, 3usize, 0.0),
        (9, 8, ux_wall),
        (10, 7, -ux_wall),
        (16, 17, uz_wall),
        (18, 15, -uz_wall),
    ];
    for (i_in, i_out, speed) in links {
        let expected = before[i_in] - 6.0 * D3Q19_WEIGHTS[i_in] * rho * speed;
        assert!(
            (f[i_out] - expected).abs() < 1.0e-15,
            "link {i_in} -> {i_out}: {}, expected {expected}",
            f[i_out]
        );
    }
    let written: Vec<usize> = links.iter().map(|link| link.1).collect();
    for i in 0..19 {
        if !written.contains(&i) {
            assert_eq!(f[i], before[i], "population {i} must not be touched");
        }
    }
    assert_eq!(f[16], before[16]);
    assert_ne!(f[17], before[17]);
}

#[test]
fn moving_top_wall_from_rest_gives_the_documented_linear_offsets() {
    let rho = 1.0;
    let ux_wall = 0.05;
    let mut f = d3q19_equilibrium(rho, 0.0, 0.0, 0.0);
    let at_rest = f;
    bounce_back_moving_top_d3q19(&mut f, rho, ux_wall, 0.0);
    let offset = 6.0 * (1.0 / 36.0) * rho * ux_wall;
    assert!(
        (f[10] - f[7] + offset).abs() < 1.0e-15,
        "f10 - f7 = {}",
        f[10] - f[7]
    );
    assert!(
        (f[9] - f[8] - offset).abs() < 1.0e-15,
        "f9 - f8 = {}",
        f[9] - f[8]
    );
    // With no z motion the z links reduce to plain bounce-back.
    assert_eq!(f[18], at_rest[15]);
    assert_eq!(f[16], at_rest[17]);
}

#[test]
fn moving_wall_correction_is_linear_in_the_wall_velocity() {
    let rho = 1.0;
    let base = with_noise(d3q19_equilibrium(rho, 0.0, 0.0, 0.0), 0.004, 5);
    for bounce in [
        bounce_back_moving_top_d3q19,
        bounce_back_moving_bottom_d3q19,
    ] {
        let mut plain = base;
        bounce(&mut plain, rho, 0.0, 0.0);
        let deviation = |ux: f64, uz: f64| {
            let mut f = base;
            bounce(&mut f, rho, ux, uz);
            let mut d = [0.0; 19];
            for (di, (&fi, &pi)) in d.iter_mut().zip(f.iter().zip(plain.iter())) {
                *di = fi - pi;
            }
            d
        };
        let single = deviation(0.02, 0.01);
        let double = deviation(0.04, 0.02);
        let reversed = deviation(-0.02, -0.01);
        assert!(
            single.iter().any(|d| d.abs() > 1.0e-4),
            "wall velocity has no effect"
        );
        for i in 0..19 {
            assert!(
                (double[i] - 2.0 * single[i]).abs() < 1.0e-15,
                "doubling, population {i}"
            );
            assert!(
                (reversed[i] + single[i]).abs() < 1.0e-15,
                "reversal, population {i}"
            );
        }
    }
}

#[test]
fn moving_wall_injects_momentum_rho_u_over_three_parallel_to_the_wall() {
    // The reflected populations carry the extra momentum
    // sum 6 w rho (c . u) c = rho u / 3 along the wall and nothing normal to it.
    let rho = 1.07;
    let base = with_noise(d3q19_equilibrium(rho, 0.0, 0.0, 0.0), 0.004, 21);
    let (ux_wall, uz_wall) = (0.04, -0.025);
    for bounce in [
        bounce_back_moving_top_d3q19,
        bounce_back_moving_bottom_d3q19,
    ] {
        let mut plain = base;
        bounce(&mut plain, rho, 0.0, 0.0);
        let mut moving = base;
        bounce(&mut moving, rho, ux_wall, uz_wall);
        let mut difference = [0.0; 19];
        for (d, (&m, &p)) in difference.iter_mut().zip(moving.iter().zip(plain.iter())) {
            *d = m - p;
        }
        let (mass, momentum) = moments(&difference, &D3Q19_VELOCITIES);
        assert!(
            mass.abs() < 1.0e-15,
            "wall motion changed the mass by {mass}"
        );
        assert!(
            (momentum[0] - rho * ux_wall / 3.0).abs() < 1.0e-14,
            "jx = {}",
            momentum[0]
        );
        assert!(momentum[1].abs() < 1.0e-14, "jy = {}", momentum[1]);
        assert!(
            (momentum[2] - rho * uz_wall / 3.0).abs() < 1.0e-14,
            "jz = {}",
            momentum[2]
        );
    }
}

// ---------------------------------------------------------------------------
// Flow-level checks on one-dimensional D3Q19 columns
// ---------------------------------------------------------------------------

/// One BGK collision on a single D3Q19 node.
fn bgk_collide_d3q19(f: &mut [f64; 19], omega: f64) {
    let (rho, ux, uy, uz) = macroscopic_d3q19(f);
    let feq = d3q19_equilibrium(rho, ux, uy, uz);
    for (fi, &ei) in f.iter_mut().zip(feq.iter()) {
        *fi += omega * (ei - *fi);
    }
}

/// Stream a column of nodes along `axis`.  The two other directions are
/// periodic with period one, so populations that do not move along `axis`
/// stay on their node.  Populations entering from outside the column are
/// left at zero: they are the unknowns of the boundary conditions.
fn stream_column(f: &[[f64; 19]], axis: usize) -> Vec<[f64; 19]> {
    let n = f.len();
    let mut next = vec![[0.0; 19]; n];
    for (position, node) in f.iter().enumerate() {
        for (i, c) in D3Q19_VELOCITIES.iter().enumerate() {
            if let Some(target) = position.checked_add_signed(c[axis] as isize)
                && target < n
            {
                next[target][i] = node[i];
            }
        }
    }
    next
}

/// Plane Couette flow between a stationary bottom wall and a top wall moving
/// at `u_top`, with half-way bounce-back walls half a cell outside the first
/// and last fluid node.  Returns the populations after `steps` steps.
fn couette_with_moving_wall_bounce_back(
    ny: usize,
    tau: f64,
    u_top: f64,
    steps: usize,
) -> Vec<[f64; 19]> {
    let omega = 1.0 / tau;
    let mut f = vec![d3q19_equilibrium(1.0, 0.0, 0.0, 0.0); ny];
    for _ in 0..steps {
        for node in &mut f {
            bgk_collide_d3q19(node, omega);
        }
        let mut next = stream_column(&f, 1);
        let rho_top: f64 = f[ny - 1].iter().sum();
        let mut reflected_top = f[ny - 1];
        bounce_back_moving_top_d3q19(&mut reflected_top, rho_top, u_top, 0.0);
        let rho_bottom: f64 = f[0].iter().sum();
        let mut reflected_bottom = f[0];
        bounce_back_moving_bottom_d3q19(&mut reflected_bottom, rho_bottom, 0.0, 0.0);
        for (i, c) in D3Q19_VELOCITIES.iter().enumerate() {
            if c[1] == -1 {
                next[ny - 1][i] = reflected_top[i];
            }
            if c[1] == 1 {
                next[0][i] = reflected_bottom[i];
            }
        }
        f = next;
    }
    f
}

#[test]
fn couette_flow_with_moving_wall_bounce_back_is_linear() {
    let ny = 12;
    let u_top = 0.02;
    let f = couette_with_moving_wall_bounce_back(ny, 0.9, u_top, 6000);
    for (y, node) in f.iter().enumerate() {
        let (rho, ux, uy, uz) = macroscopic_d3q19(node);
        let expected = u_top * (y as f64 + 0.5) / ny as f64;
        assert!(
            (ux - expected).abs() < 1.0e-9 * u_top,
            "y = {y}: ux = {ux}, expected {expected}"
        );
        assert!(
            uy.abs() < 1.0e-12 && uz.abs() < 1.0e-12,
            "y = {y}: uy = {uy}, uz = {uz}"
        );
        assert!((rho - 1.0).abs() < 1.0e-9, "y = {y}: rho = {rho}");
    }
}

/// Plane Couette flow with on-node Zou-He walls: the bottom wall (`y = 0`) at
/// rest and the top wall (`y = ny - 1`) moving at `u_top` through the
/// tangential velocity of the generic face condition.
fn couette_with_zou_he_walls(ny: usize, tau: f64, u_top: f64, steps: usize) -> Vec<[f64; 19]> {
    let omega = 1.0 / tau;
    let mut f = vec![d3q19_equilibrium(1.0, 0.0, 0.0, 0.0); ny];
    for _ in 0..steps {
        for node in &mut f {
            bgk_collide_d3q19(node, omega);
        }
        f = stream_column(&f, 1);
        zou_he_bottom_face_d3q19(&mut f[0], 0.0);
        zou_he_velocity_face_d3q19(&mut f[ny - 1], ZouHeFace::Top, [u_top, 0.0, 0.0]);
    }
    f
}

#[test]
fn couette_flow_with_zou_he_walls_carries_the_tangential_wall_velocity() {
    let ny = 11;
    let u_top = 0.02;
    let f = couette_with_zou_he_walls(ny, 0.9, u_top, 6000);
    for (y, node) in f.iter().enumerate() {
        let (_, ux, uy, uz) = macroscopic_d3q19(node);
        let expected = u_top * y as f64 / (ny - 1) as f64;
        assert!(
            (ux - expected).abs() < 1.0e-9 * u_top,
            "y = {y}: ux = {ux}, expected {expected}"
        );
        assert!(
            uy.abs() < 1.0e-12 && uz.abs() < 1.0e-12,
            "y = {y}: uy = {uy}, uz = {uz}"
        );
    }
    let (_, ux_wall, _, _) = macroscopic_d3q19(&f[ny - 1]);
    assert!(
        (ux_wall - u_top).abs() < 1.0e-12,
        "wall node ux = {ux_wall}"
    );
}

/// Plug flow in a column along x with a Zou-He velocity inlet at `x = 0` and
/// a Zou-He pressure outlet at `x = nx - 1`.
fn channel_with_inlet_and_outlet(
    mut f: Vec<[f64; 19]>,
    tau: f64,
    u_in: f64,
    rho_out: f64,
    steps: usize,
) -> Vec<[f64; 19]> {
    let omega = 1.0 / tau;
    let last = f.len() - 1;
    for _ in 0..steps {
        for node in &mut f {
            bgk_collide_d3q19(node, omega);
        }
        f = stream_column(&f, 0);
        zou_he_inlet_left_d3q19(&mut f[0], u_in);
        zou_he_outlet_right_d3q19(&mut f[last], rho_out);
    }
    f
}

#[test]
fn uniform_plug_flow_is_a_fixed_point_of_inlet_and_outlet() {
    let (u_in, rho0) = (0.04, 1.0);
    let uniform = vec![d3q19_equilibrium(rho0, u_in, 0.0, 0.0); 9];
    let f = channel_with_inlet_and_outlet(uniform.clone(), 0.8, u_in, rho0, 200);
    for (x, (node, initial)) in f.iter().zip(uniform.iter()).enumerate() {
        for (i, (&got, &want)) in node.iter().zip(initial.iter()).enumerate() {
            assert!(
                (got - want).abs() < TOL,
                "x = {x}: f[{i}] = {got}, expected {want}"
            );
        }
    }
}

#[test]
fn inlet_and_outlet_drive_a_channel_to_uniform_plug_flow() {
    let (u_in, rho_out) = (0.03, 1.0);
    let rest = vec![d3q19_equilibrium(1.0, 0.0, 0.0, 0.0); 16];
    let f = channel_with_inlet_and_outlet(rest, 0.8, u_in, rho_out, 20000);
    for (x, node) in f.iter().enumerate() {
        let (rho, ux, uy, uz) = macroscopic_d3q19(node);
        assert!(
            (ux - u_in).abs() < 1.0e-9,
            "x = {x}: ux = {ux}, expected {u_in}"
        );
        assert!((rho - rho_out).abs() < 1.0e-9, "x = {x}: rho = {rho}");
        assert!(
            uy.abs() < 1.0e-12 && uz.abs() < 1.0e-12,
            "x = {x}: uy = {uy}, uz = {uz}"
        );
    }
}
