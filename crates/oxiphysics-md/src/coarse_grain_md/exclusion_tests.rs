//! Regression tests for the bonded exclusions and the Yukawa electrostatics of the
//! coarse-grained simulation: forces, potential energy and virial pressure share one
//! pair kernel and must agree with each other and with closed-form values.

use super::*;

/// Gas constant in kJ mol⁻¹ K⁻¹ (CODATA 2018), independent of the module constant.
const R_KJ_MOL_K: f64 = 8.314_462_618e-3;

fn bead(name: &str, type_: CgBeadType, pos: [f64; 3], charge: f64) -> CgBead {
    let mut b = CgBead::new(name, type_, 72.0, charge);
    b.pos = pos;
    b
}

fn single(name: &str, type_: CgBeadType, pos: [f64; 3], charge: f64) -> CgMolecule {
    let mut mol = CgMolecule::new();
    mol.beads.push(bead(name, type_, pos, charge));
    mol
}

/// Force field with all Lennard-Jones well depths zero: only the Yukawa term acts.
fn electrostatics_only_ff() -> CgForceField {
    let mut ff = CgForceField::default_martini();
    for row in &mut ff.lj_params {
        for entry in row.iter_mut() {
            entry.0 = 0.0;
        }
    }
    ff
}

/// Three-bead chain (bonds 0-1 and 1-2, angle 0-1-2) with the given charges, its
/// first bead at `origin`; the bonds are 0.50 and 0.43 nm long, the angle is 110°.
fn chain(origin: [f64; 3], charges: [f64; 3]) -> CgMolecule {
    let theta = 110.0_f64.to_radians();
    let middle = [origin[0] + 0.50, origin[1], origin[2]];
    let end = [
        middle[0] + 0.43 * (PI - theta).cos(),
        middle[1] + 0.43 * (PI - theta).sin(),
        middle[2],
    ];
    let mut mol = CgMolecule::new();
    mol.beads
        .push(bead("C0", CgBeadType::Backbone, origin, charges[0]));
    mol.beads
        .push(bead("C1", CgBeadType::Backbone, middle, charges[1]));
    mol.beads
        .push(bead("C2", CgBeadType::Backbone, end, charges[2]));
    mol.bonds.push((0, 1, 0.47, 1250.0));
    mol.bonds.push((1, 2, 0.47, 1250.0));
    mol.angles.push((0, 1, 2, 95.0_f64.to_radians(), 25.0));
    mol
}

/// A 4 nm box with a charged chain, charged single beads and a charged pair that
/// interacts through the periodic boundary in all three directions.  Every
/// non-excluded separation is above the 0.1 nm distance floor.
fn charged_system() -> CgSimulation {
    let molecules = vec![
        chain([1.0, 1.0, 1.0], [1.0, 0.0, -1.0]),
        single("P", CgBeadType::Polar, [1.2, 1.55, 1.1], 1.0),
        single("Q", CgBeadType::Charged, [0.7, 0.55, 1.35], -1.0),
        single("A", CgBeadType::Apolar, [1.8, 0.4, 0.8], 0.0),
        single("X", CgBeadType::Charged, [3.7, 3.8, 3.9], 0.5),
        single("Y", CgBeadType::Charged, [0.1, 0.2, 0.15], -0.5),
    ];
    CgSimulation::new(
        molecules,
        CgForceField::default_martini(),
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box")
}

#[test]
fn bonded_pairs_and_angle_ends_are_excluded_from_the_non_bonded_forces() {
    // A chain compressed well inside the LJ core (all separations ~0.30 nm < sigma =
    // 0.47 nm): without exclusions the 1-2 and 1-3 pairs would feel thousands of
    // kJ/mol/nm of LJ repulsion on top of the bonded terms.
    let mut mol = chain([1.0, 1.0, 1.0], [0.0; 3]);
    mol.beads[0].pos = [1.0, 1.0, 1.0];
    mol.beads[1].pos = [1.30, 1.0, 1.0];
    mol.beads[2].pos = [1.45, 1.26, 1.0];
    let bonded: Vec<[f64; 3]> = mol
        .bond_forces()
        .iter()
        .zip(mol.angle_forces().iter())
        .map(|(b, a)| add3(*b, *a))
        .collect();
    let sim = CgSimulation::new(
        vec![mol],
        CgForceField::default_martini(),
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box");
    let forces = sim.compute_all_forces();
    for (bead_index, (f, b)) in forces.iter().zip(bonded.iter()).enumerate() {
        for axis in 0..3 {
            assert_eq!(
                f[axis], b[axis],
                "bead {bead_index} axis {axis}: non-bonded force on an excluded pair"
            );
        }
    }
    // The non-bonded energy of the (fully excluded) molecule is zero.
    let bonded_energy = sim.molecules[0].bonded_energy();
    assert_eq!(sim.total_potential_energy(), bonded_energy);
    // The same beads in *different* molecules do interact.
    let split = CgSimulation::new(
        sim.molecules[0]
            .beads
            .iter()
            .map(|b| single(&b.name, b.type_, b.pos, b.charge))
            .collect(),
        CgForceField::default_martini(),
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box");
    let split_forces = split.compute_all_forces();
    assert!(
        split_forces[0][0] < -100.0,
        "unbonded beads at 0.30 nm repel strongly: {:?}",
        split_forces[0]
    );
}

#[test]
fn exclusions_use_global_bead_indices_across_molecules() {
    let sim = CgSimulation::new(
        vec![
            single("S", CgBeadType::Polar, [0.5; 3], 0.0),
            chain([1.0, 1.0, 1.0], [0.0; 3]),
            single("T", CgBeadType::Polar, [2.5; 3], 0.0),
            chain([2.0, 2.0, 2.0], [0.0; 3]),
        ],
        CgForceField::default_martini(),
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box");
    let exclusions = sim.exclusions();
    // Beads: 0 = S, 1-3 = first chain, 4 = T, 5-7 = second chain.
    let expected = [(1, 2), (2, 3), (1, 3), (5, 6), (6, 7), (5, 7)];
    for i in 0..8 {
        for j in (i + 1)..8 {
            assert_eq!(
                exclusions.contains(i, j),
                expected.contains(&(i, j)),
                "pair ({i}, {j})"
            );
        }
    }
    assert!(
        !exclusions.contains(3, 5),
        "chains do not exclude each other"
    );
    assert!(!exclusions.contains(100, 101));
}

#[test]
fn yukawa_force_and_energy_match_the_closed_form() {
    // Two beads of different molecules, only the Yukawa term active, r = 0.8 nm.
    let ff = electrostatics_only_ff();
    let lambda = ff.electrostatic_screening;
    let r = 0.8;
    let sim = CgSimulation::new(
        vec![
            single("P", CgBeadType::Charged, [1.0, 2.0, 2.0], 1.0),
            single("N", CgBeadType::Charged, [1.0 + r, 2.0, 2.0], -1.0),
        ],
        ff,
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box");
    let ke = 138.935;
    let screening = (-r / lambda).exp();
    let energy = -ke * screening / r;
    assert!((sim.total_potential_energy() - energy).abs() < 1e-12 * energy.abs());

    // -dU/dr = K q q e^{-r/lambda} (1/r^2 + 1/(lambda r)); opposite charges attract.
    let attraction = ke * screening * (1.0 / (r * r) + 1.0 / (lambda * r));
    let forces = sim.compute_all_forces();
    assert!(
        (forces[0][0] - attraction).abs() < 1e-10 * attraction,
        "force on the + bead {:?}, expected +{attraction} along x",
        forces[0]
    );
    assert!((forces[1][0] + attraction).abs() < 1e-10 * attraction);
    assert!(forces[0][1].abs() < 1e-12 && forces[0][2].abs() < 1e-12);

    // Like charges repel.
    let mut repelling = sim.clone();
    repelling.molecules[1].beads[0].charge = 1.0;
    let forces = repelling.compute_all_forces();
    assert!(forces[0][0] < 0.0 && forces[1][0] > 0.0);
    // A neutral partner carries no electrostatic force at all.
    repelling.molecules[1].beads[0].charge = 0.0;
    let forces = repelling.compute_all_forces();
    assert_eq!(forces[0], [0.0; 3]);
}

#[test]
fn electrostatics_use_the_minimum_image() {
    let ff = electrostatics_only_ff();
    let across = CgSimulation::new(
        vec![
            single("P", CgBeadType::Charged, [3.8, 2.0, 2.0], 1.0),
            single("N", CgBeadType::Charged, [0.2, 2.0, 2.0], -1.0),
        ],
        ff.clone(),
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box");
    let inside = CgSimulation::new(
        vec![
            single("P", CgBeadType::Charged, [1.8, 2.0, 2.0], 1.0),
            single("N", CgBeadType::Charged, [2.2, 2.0, 2.0], -1.0),
        ],
        ff,
        300.0,
        0.02,
        [4.0; 3],
    )
    .expect("valid box");
    let f_across = across.compute_all_forces();
    let f_inside = inside.compute_all_forces();
    for bead_index in 0..2 {
        for axis in 0..3 {
            assert!(
                (f_across[bead_index][axis] - f_inside[bead_index][axis]).abs() < 1e-9,
                "bead {bead_index} axis {axis}"
            );
        }
    }
    // The opposite charges attract through the boundary: the + bead at x = 3.8 is pulled
    // towards +x (towards the image of the - bead at x = 4.2).
    assert!(f_across[0][0] > 0.0);
    assert!((across.total_potential_energy() - inside.total_potential_energy()).abs() < 1e-9);
}

#[test]
fn forces_are_the_negative_gradient_of_the_total_potential_energy() {
    let sim = charged_system();
    let forces = sim.compute_all_forces();
    let step = 1.0e-6;
    let scale = forces
        .iter()
        .flat_map(|f| f.iter())
        .fold(1.0_f64, |m, v| m.max(v.abs()));
    assert!(scale > 10.0, "the test system must carry sizeable forces");

    let mut flat = 0;
    for mol_index in 0..sim.molecules.len() {
        for bead_index in 0..sim.molecules[mol_index].beads.len() {
            for (axis, &analytic) in forces[flat].iter().enumerate() {
                let mut plus = sim.clone();
                plus.molecules[mol_index].beads[bead_index].pos[axis] += step;
                let mut minus = sim.clone();
                minus.molecules[mol_index].beads[bead_index].pos[axis] -= step;
                let numeric = -(plus.total_potential_energy() - minus.total_potential_energy())
                    / (2.0 * step);
                assert!(
                    (analytic - numeric).abs() <= 1e-5 * scale,
                    "molecule {mol_index} bead {bead_index} axis {axis}: F = {analytic}, -dE/dx = {numeric}"
                );
            }
            flat += 1;
        }
    }
}

#[test]
fn total_force_vanishes_by_newtons_third_law() {
    let forces = charged_system().compute_all_forces();
    let mut total = [0.0_f64; 3];
    let mut scale = 0.0_f64;
    for f in &forces {
        for (sum, component) in total.iter_mut().zip(f) {
            *sum += component;
            scale = scale.max(component.abs());
        }
    }
    for (axis, component) in total.iter().enumerate() {
        assert!(
            component.abs() < 1e-9 * scale,
            "net force component {axis} = {component} (largest force {scale})"
        );
    }
}

/// The simulation with all positions and the box edges scaled by `s`.
fn scaled(sim: &CgSimulation, s: f64) -> CgSimulation {
    let mut out = sim.clone();
    for mol in &mut out.molecules {
        for b in &mut mol.beads {
            b.pos = scale3(b.pos, s);
        }
    }
    out.box_lengths = scale3(sim.box_lengths, s);
    out
}

#[test]
fn virial_pressure_includes_electrostatics_and_respects_the_exclusions() {
    // P = N k_B T / V - dU/dV for a uniform dilation of positions and box.
    let sim = charged_system();
    let delta = 1.0e-6;
    let u_plus = scaled(&sim, 1.0 + delta).total_potential_energy();
    let u_minus = scaled(&sim, 1.0 - delta).total_potential_energy();
    let volume = sim.volume();
    let d_u_d_volume = (u_plus - u_minus) / (2.0 * delta * 3.0 * volume);
    let n = sim.molecules.iter().map(|m| m.beads.len()).sum::<usize>() as f64;
    let ideal = n * R_KJ_MOL_K * 300.0 / volume;
    let expected = ideal - d_u_d_volume;
    let got = sim.pressure();
    assert!(
        ((got - expected) / expected.abs().max(1.0)).abs() < 1e-6,
        "virial pressure {got} vs ideal - dU/dV = {expected}"
    );
    assert!(
        (got - ideal).abs() > 0.05 * ideal,
        "the configurational part must matter in this test"
    );
}

#[test]
fn total_potential_energy_matches_a_hand_computed_pair_and_bond() {
    // A bonded pair (excluded from LJ) at 0.57 nm: only the harmonic bond acts.
    let mut mol = CgMolecule::new();
    mol.beads
        .push(bead("B0", CgBeadType::Apolar, [1.0, 1.0, 1.0], 0.0));
    mol.beads
        .push(bead("B1", CgBeadType::Apolar, [1.57, 1.0, 1.0], 0.0));
    mol.bonds.push((0, 1, 0.47, 3800.0));
    // A third bead of another molecule at 0.60 nm from B1 interacts through LJ only.
    let third = single("B2", CgBeadType::Apolar, [2.17, 1.0, 1.0], 0.0);
    let ff = CgForceField::default_martini();
    let (eps, sigma) = ff.lj_params[2][2];
    let sim = CgSimulation::new(vec![mol, third], ff, 300.0, 0.02, [4.0; 3]).expect("valid box");
    let bond = 0.5 * 3800.0 * 0.1_f64.powi(2);
    let lj = |r: f64| 4.0 * eps * ((sigma / r).powi(12) - (sigma / r).powi(6));
    let expected = bond + lj(0.60) + lj(1.17);
    let got = sim.total_potential_energy();
    assert!(
        (got - expected).abs() < 1e-9 * expected.abs(),
        "E = {got}, expected {expected}"
    );
}
