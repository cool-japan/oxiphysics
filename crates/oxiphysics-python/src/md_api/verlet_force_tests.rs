// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Regression tests of the Verlet-list Lennard-Jones forces, the Ewald real-space sum and
//! the velocity-Verlet dynamics of [`PyMdSimulation`](super::PyMdSimulation).

use super::*;

/// Deterministic generator of uniform numbers in `[0, 1)` (64-bit LCG).
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
    }
}

/// `count` atoms on a jittered simple-cubic lattice that fills the box, so
/// that no two atoms come closer than `spacing - 2 * jitter`.
fn jittered_lattice(count: usize, box_size: f64, jitter: f64, seed: u64) -> Vec<[f64; 3]> {
    let per_edge = (count as f64).cbrt().ceil() as usize;
    let spacing = box_size / per_edge as f64;
    let mut rng = Lcg(seed);
    let mut out = Vec::with_capacity(count);
    'fill: for ix in 0..per_edge {
        for iy in 0..per_edge {
            for iz in 0..per_edge {
                if out.len() == count {
                    break 'fill;
                }
                let mut p = [0.0; 3];
                for (k, index) in [ix, iy, iz].into_iter().enumerate() {
                    p[k] = (index as f64 + 0.5) * spacing + jitter * (2.0 * rng.next_unit() - 1.0);
                }
                out.push(p);
            }
        }
    }
    out
}

/// The original all-pairs implementation (single-shift minimum image),
/// kept as an independent reference for the neighbour-list forces.
fn all_pairs_reference(positions: &[[f64; 3]], cfg: &PyMdConfig) -> (Vec<[f64; 3]>, f64) {
    let n = positions.len();
    let mut forces = vec![[0.0f64; 3]; n];
    let eps = cfg.lj_epsilon;
    let sig = cfg.lj_sigma;
    let rc = cfg.cutoff;
    let rc2 = rc * rc;
    let box_size = cfg.box_size;
    let sig_rc2 = (sig / rc).powi(2);
    let sig_rc6 = sig_rc2 * sig_rc2 * sig_rc2;
    let u_shift = 4.0 * eps * sig_rc6 * (sig_rc6 - 1.0);
    let mut u_total = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let mut dr = [0.0f64; 3];
            for k in 0..3 {
                let mut d = positions[j][k] - positions[i][k];
                if d > 0.5 * box_size {
                    d -= box_size;
                } else if d < -0.5 * box_size {
                    d += box_size;
                }
                dr[k] = d;
            }
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if r2 >= rc2 || r2 < 1e-20 {
                continue;
            }
            let sig2_r2 = (sig * sig) / r2;
            let sig6_r6 = sig2_r2 * sig2_r2 * sig2_r2;
            let sig12_r12 = sig6_r6 * sig6_r6;
            let f_mag = 24.0 * eps / r2 * (2.0 * sig12_r12 - sig6_r6);
            for k in 0..3 {
                forces[i][k] -= f_mag * dr[k];
                forces[j][k] += f_mag * dr[k];
            }
            u_total += 4.0 * eps * sig6_r6 * (sig6_r6 - 1.0) - u_shift;
        }
    }
    (forces, u_total)
}

fn sim_with_atoms(cfg: &PyMdConfig, positions: &[[f64; 3]]) -> PyMdSimulation {
    let mut sim = PyMdSimulation::new(cfg.clone());
    sim.set_thermostat(false);
    for p in positions {
        sim.add_atom(p.to_vec(), 0);
    }
    sim
}

fn assert_matches_reference(sim: &PyMdSimulation, cfg: &PyMdConfig, context: &str) {
    let positions: Vec<[f64; 3]> = sim.atoms.iter().map(|a| a.position).collect();
    let (ref_forces, ref_energy) = all_pairs_reference(&positions, cfg);
    assert!(
        (sim.potential_energy() - ref_energy).abs() <= 1e-10 * ref_energy.abs().max(1.0),
        "{context}: potential energy {} vs reference {ref_energy}",
        sim.potential_energy()
    );
    for (i, (atom, expected)) in sim.atoms.iter().zip(&ref_forces).enumerate() {
        for k in 0..3 {
            assert!(
                (atom.force[k] - expected[k]).abs() <= 1e-10 * expected[k].abs().max(1.0),
                "{context}: force component {k} of atom {i}: {} vs reference {}",
                atom.force[k],
                expected[k]
            );
        }
    }
}

#[test]
fn test_neighbor_list_forces_match_all_pairs_for_500_atoms() {
    let cfg = PyMdConfig::default();
    let positions = jittered_lattice(500, cfg.box_size, 0.1, 7);
    let mut sim = sim_with_atoms(&cfg, &positions);
    sim.compute_forces();
    let state = sim
        .neighbor_list
        .as_ref()
        .expect("Verlet list must be used");
    assert!(
        state.list.pair_count() > 500,
        "the list must hold real pairs"
    );
    assert_matches_reference(&sim, &cfg, "fresh list");
}

#[test]
fn test_neighbor_list_forces_follow_moving_atoms() {
    let cfg = PyMdConfig::default();
    let positions = jittered_lattice(500, cfg.box_size, 0.1, 11);
    let mut sim = sim_with_atoms(&cfg, &positions);
    let mut rng = Lcg(5);
    for i in 0..sim.atom_count() {
        let v: Vec<f64> = (0..3).map(|_| 2.0 * rng.next_unit() - 1.0).collect();
        sim.set_velocity(i, v);
    }
    for step in 0..120 {
        sim.step(0.005);
        assert_matches_reference(&sim, &cfg, &format!("step {step}"));
    }
    let state = sim
        .neighbor_list
        .as_ref()
        .expect("Verlet list must be used");
    assert!(
        state.list.rebuild_count > 0,
        "atoms moved more than half the skin, so the list must have been rebuilt"
    );
    assert!(
        state.list.rebuild_count < 60,
        "the list must be reused between steps, rebuilt {} times in 120 steps",
        state.list.rebuild_count
    );
}

#[test]
fn test_small_box_falls_back_to_all_pairs() {
    let cfg = PyMdConfig {
        box_size: 6.0,
        ..PyMdConfig::default()
    };
    let positions = jittered_lattice(64, cfg.box_size, 0.05, 3);
    let mut sim = sim_with_atoms(&cfg, &positions);
    sim.compute_forces();
    assert!(
        sim.neighbor_list.is_none(),
        "6 / (2.5 + 0.3) < 3 cells per edge: no Verlet list"
    );
    assert_matches_reference(&sim, &cfg, "all-pairs fallback");
}

#[test]
fn test_add_atom_invalidates_the_neighbor_list() {
    let cfg = PyMdConfig::default();
    let positions = jittered_lattice(200, cfg.box_size, 0.1, 13);
    let mut sim = sim_with_atoms(&cfg, &positions);
    sim.compute_forces();
    assert!(sim.neighbor_list.is_some());
    // A new atom close to an existing one must feel (and exert) forces.
    let p = positions[17];
    sim.add_atom(vec![p[0] + 1.1, p[1], p[2]], 0);
    assert!(sim.neighbor_list.is_none());
    assert!(!sim.forces_valid);
    sim.ensure_forces();
    assert_matches_reference(&sim, &cfg, "after add_atom");
    let newest = sim.atoms.last().expect("atom was added");
    assert!(
        newest.force.iter().any(|f| f.abs() > 1e-6),
        "the added atom must interact with its neighbours"
    );
}

#[test]
fn test_first_step_uses_forces_at_the_initial_positions() {
    let mut sim = PyMdSimulation::new(PyMdConfig::default());
    sim.set_thermostat(false);
    sim.add_atom(vec![5.0, 5.0, 5.0], 0);
    sim.add_atom(vec![6.5, 5.0, 5.0], 0);
    let dt = 1.0e-3;
    sim.step(dt);
    // At r = 1.5 sigma the pair attracts with F = 24 (2 r^-13 - r^-7) = -1.158
    // (negative: attractive), so atom 0 gains v = F dt in +x over one step.
    let expected = 1.158_03 * dt;
    let v0 = sim.velocity(0).expect("atom 0")[0];
    let v1 = sim.velocity(1).expect("atom 1")[0];
    assert!(
        (v0 - expected).abs() < 0.02 * expected,
        "v0x = {v0}, expected {expected}"
    );
    assert!((v0 + v1).abs() < 1e-12, "equal and opposite momentum");
}

/// The original all-pairs Ewald real-space sum, kept as a reference.
fn ewald_all_pairs_reference(set: &AtomSet, alpha: f64, r_cut: f64) -> f64 {
    let n = set.len();
    let box_size = set.box_size;
    let mut energy = 0.0;
    for i in 0..n {
        let qi = set.charge(i);
        if qi == 0.0 {
            continue;
        }
        for j in (i + 1)..n {
            let qj = set.charge(j);
            if qj == 0.0 {
                continue;
            }
            let mut dr = [0.0f64; 3];
            for k in 0..3 {
                let mut d = set.positions[j][k] - set.positions[i][k];
                if d > 0.5 * box_size {
                    d -= box_size;
                } else if d < -0.5 * box_size {
                    d += box_size;
                }
                dr[k] = d;
            }
            let r2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
            if r2 >= r_cut * r_cut || r2 < 1e-20 {
                continue;
            }
            let r = r2.sqrt();
            energy += qi * qj * erfc_approx(alpha * r) / r;
        }
    }
    energy
}

/// A salt-like set: alternating Na+/Cl- on a jittered lattice, with some
/// uncharged argon atoms mixed in.
fn salt_set(count: usize, box_size: f64, seed: u64) -> AtomSet {
    let mut set = AtomSet::new(box_size);
    let na = set.register_type(AtomTypeDesc::sodium_ion());
    let cl = set.register_type(AtomTypeDesc::chloride_ion());
    let ar = set.register_type(AtomTypeDesc::argon());
    for (k, p) in jittered_lattice(count, box_size, 0.1, seed)
        .into_iter()
        .enumerate()
    {
        let species = match k % 3 {
            0 => na,
            1 => cl,
            _ => ar,
        };
        set.add_atom(p.to_vec(), species);
    }
    set
}

#[test]
fn test_ewald_real_space_cell_list_matches_all_pairs() {
    // 20 / 3 >= 3 cut-off lengths per edge: the cell list is used.
    let set = salt_set(343, 20.0, 31);
    let expected = ewald_all_pairs_reference(&set, 0.5, 3.0);
    let actual = ewald_real_space_energy(&set, 0.5, 3.0);
    assert!(expected.abs() > 1e-3, "the reference sum is not trivial");
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
        "{actual} vs reference {expected}"
    );
}

#[test]
fn test_ewald_real_space_small_box_falls_back_to_all_pairs() {
    // 8 / 3 < 3 cut-off lengths per edge: every charged pair is a candidate.
    let set = salt_set(64, 8.0, 41);
    let expected = ewald_all_pairs_reference(&set, 0.4, 3.0);
    let actual = ewald_real_space_energy(&set, 0.4, 3.0);
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
        "{actual} vs reference {expected}"
    );
}

#[test]
fn test_energy_is_conserved_by_microcanonical_dynamics() {
    let cfg = PyMdConfig::default();
    let positions = jittered_lattice(125, cfg.box_size, 0.1, 21);
    let mut sim = sim_with_atoms(&cfg, &positions);
    let mut rng = Lcg(9);
    for i in 0..sim.atom_count() {
        let v: Vec<f64> = (0..3)
            .map(|_| 0.8 * (2.0 * rng.next_unit() - 1.0))
            .collect();
        sim.set_velocity(i, v);
    }
    sim.ensure_forces();
    let e0 = sim.total_energy();
    for _ in 0..400 {
        sim.step(0.002);
    }
    let drift = (sim.total_energy() - e0).abs() / e0.abs().max(1.0);
    assert!(drift < 1e-3, "relative energy drift {drift}");
}
