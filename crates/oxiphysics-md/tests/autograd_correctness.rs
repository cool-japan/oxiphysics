// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the autograd bridge:
//! analytical gradient correctness verified against central finite differences.

use oxiphysics_core::math::Vec3;
use oxiphysics_md::autograd_bridge::{
    DifferentiableForceField, NnGrads, bp_descriptor_forces, nn_backward, nn_forward_with_tape,
};
use oxiphysics_md::ml_potential::{
    BehlerParrinelloDescriptor, BpNeuralNetworkPotential, NeuralNetworkLayer, identity_activation,
    identity_derivative, relu_activation, relu_derivative, tanh_activation, tanh_derivative,
};

// ---------------------------------------------------------------------------
// Test 1: activation derivative vs finite difference
// ---------------------------------------------------------------------------

#[test]
fn test_activation_derivatives_match_finite_diff() {
    let h = 1e-5_f64;
    let tol = 1e-7_f64;

    // Sample 16 points from -3.0 to 3.0 (step 0.4)
    let xs: Vec<f64> = (0..16).map(|i| -3.0 + 0.4 * i as f64).collect();

    // tanh
    for &x in &xs {
        let analytical = tanh_derivative(x);
        let numerical = (tanh_activation(x + h) - tanh_activation(x - h)) / (2.0 * h);
        assert!(
            (analytical - numerical).abs() < tol,
            "tanh_derivative at x={x}: analytical={analytical}, numerical={numerical}"
        );
    }

    // relu (skip x == 0, derivative undefined there)
    for &x in &xs {
        if x.abs() < 0.05 {
            continue; // skip near zero
        }
        let analytical = relu_derivative(x);
        let numerical = (relu_activation(x + h) - relu_activation(x - h)) / (2.0 * h);
        assert!(
            (analytical - numerical).abs() < tol,
            "relu_derivative at x={x}: analytical={analytical}, numerical={numerical}"
        );
    }

    // identity
    for &x in &xs {
        let analytical = identity_derivative(x);
        let numerical = (identity_activation(x + h) - identity_activation(x - h)) / (2.0 * h);
        assert!(
            (analytical - numerical).abs() < tol,
            "identity_derivative at x={x}: analytical={analytical}, numerical={numerical}"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 2: NeuralNetworkLayer::with_explicit_derivative matches new()
// ---------------------------------------------------------------------------

#[test]
fn test_nn_layer_with_explicit_derivative() {
    let weights = vec![vec![0.5, -0.3], vec![0.1, 0.7]];
    let biases = vec![0.2, -0.1];
    let input = vec![1.0, 0.5];

    let layer_auto = NeuralNetworkLayer::new(weights.clone(), biases.clone(), tanh_activation);
    let layer_explicit = NeuralNetworkLayer::with_explicit_derivative(
        weights.clone(),
        biases.clone(),
        tanh_activation,
        tanh_derivative,
    );

    let out_auto = layer_auto.forward(&input);
    let out_explicit = layer_explicit.forward(&input);

    assert_eq!(out_auto.len(), out_explicit.len());
    for (a, b) in out_auto.iter().zip(out_explicit.iter()) {
        assert!(
            (a - b).abs() < 1e-14,
            "forward mismatch: auto={a}, explicit={b}"
        );
    }

    // Also verify the dactivation field is the same function pointer
    let x = 0.7_f64;
    let d_auto = (layer_auto.dactivation)(x);
    let d_explicit = (layer_explicit.dactivation)(x);
    assert!(
        (d_auto - d_explicit).abs() < 1e-14,
        "dactivation mismatch: auto={d_auto}, explicit={d_explicit}"
    );
}

// ---------------------------------------------------------------------------
// Test 3: NN backward single-layer grad check (direct NN functions)
// ---------------------------------------------------------------------------

#[test]
fn test_nn_backward_single_layer_grad_check() {
    // 1 layer: 3 inputs → 2 outputs, tanh activation
    let weights = vec![vec![0.1, -0.2, 0.3], vec![0.4, 0.1, -0.1]];
    let biases = vec![0.0, 0.0];
    let layer = NeuralNetworkLayer::new(weights.clone(), biases.clone(), tanh_activation);
    let layers = vec![layer];
    let input = vec![1.0, 0.5, -0.3];

    let h = 1e-5_f64;
    let tol = 1e-5_f64;

    // Analytical gradients
    let (e0, tape) = nn_forward_with_tape(&layers, &input);
    let _ = e0;
    let grads: NnGrads = nn_backward(&layers, &tape);

    // Verify dW by finite difference
    for i in 0..2 {
        for j in 0..3 {
            let mut layers_plus = layers.clone();
            layers_plus[0].weights[i][j] += h;
            let (e_plus, _) = nn_forward_with_tape(&layers_plus, &input);

            let mut layers_minus = layers.clone();
            layers_minus[0].weights[i][j] -= h;
            let (e_minus, _) = nn_forward_with_tape(&layers_minus, &input);

            let fd = (e_plus - e_minus) / (2.0 * h);
            let analytical = grads.dw[0][i][j];
            assert!(
                (analytical - fd).abs() < tol,
                "dW[{i}][{j}]: analytical={analytical:.8}, fd={fd:.8}, diff={:.2e}",
                (analytical - fd).abs()
            );
        }
    }

    // Verify db by finite difference
    for i in 0..2 {
        let mut layers_plus = layers.clone();
        layers_plus[0].biases[i] += h;
        let (e_plus, _) = nn_forward_with_tape(&layers_plus, &input);

        let mut layers_minus = layers.clone();
        layers_minus[0].biases[i] -= h;
        let (e_minus, _) = nn_forward_with_tape(&layers_minus, &input);

        let fd = (e_plus - e_minus) / (2.0 * h);
        let analytical = grads.db[0][i];
        assert!(
            (analytical - fd).abs() < tol,
            "db[{i}]: analytical={analytical:.8}, fd={fd:.8}"
        );
    }
}

// ---------------------------------------------------------------------------
// Test 4: BP G2 gradient check (2-atom system)
// ---------------------------------------------------------------------------

#[test]
fn test_bp_g2_gradient_check() {
    // 2-atom system: atom 0 at origin, atom 1 at (1.5, 0, 0)
    let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
    let species = vec![1u32, 1];

    // Single G2 symmetry function
    let mut desc = BehlerParrinelloDescriptor::new(3.0);
    desc.add_g2(0.5, 0.0);

    // Simple 1→2→1 network
    let l1 = NeuralNetworkLayer::new(vec![vec![0.5], vec![-0.3]], vec![0.0, 0.0], tanh_activation);
    let l2 = NeuralNetworkLayer::new(vec![vec![0.4, 0.6]], vec![0.0], identity_activation);
    let pot = BpNeuralNetworkPotential::new(vec![l1, l2], desc);

    let (_, forces, _) = pot.energy_forces_param_grads(&positions, &species);

    // Verify forces via central finite difference on energy
    let h = 1e-5_f64;
    let tol = 1e-5_f64;
    let pos_arr: Vec<[f64; 3]> = positions.iter().map(|p| [p[0], p[1], p[2]]).collect();
    let neighbors: Vec<usize> = (0..2).collect();

    for atom in 0..2 {
        for alpha in 0..3 {
            let mut pos_plus = pos_arr.clone();
            pos_plus[atom][alpha] += h;
            let e_plus: f64 = (0..2)
                .map(|c| {
                    let d = pot.descriptor.compute(&pos_plus, c, &neighbors);
                    pot.energy(&d)
                })
                .sum();

            let mut pos_minus = pos_arr.clone();
            pos_minus[atom][alpha] -= h;
            let e_minus: f64 = (0..2)
                .map(|c| {
                    let d = pot.descriptor.compute(&pos_minus, c, &neighbors);
                    pot.energy(&d)
                })
                .sum();

            let fd_force = -(e_plus - e_minus) / (2.0 * h);
            let analytical_force = forces[atom][alpha];
            assert!(
                (analytical_force - fd_force).abs() < tol,
                "atom={atom} alpha={alpha}: analytical={analytical_force:.8}, fd={fd_force:.8}, diff={:.2e}",
                (analytical_force - fd_force).abs()
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Test 5: 4-atom cluster forces agree with finite diff
// ---------------------------------------------------------------------------

#[test]
fn test_bp_nn_forces_agree_with_fd() {
    // Tetrahedral cluster with bond length ~1.2 Å
    let s = (2.0_f64 / 3.0).sqrt() * 1.2;
    let positions = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.2, 0.0, 0.0),
        Vec3::new(0.6, s, 0.0),
        Vec3::new(
            0.6,
            s / 3.0,
            (1.2_f64.powi(2) - 0.6_f64.powi(2) - (s / 3.0).powi(2))
                .abs()
                .sqrt(),
        ),
    ];
    let species = vec![1u32; 4];

    // 2 G2 + 1 G4
    let mut desc = BehlerParrinelloDescriptor::new(3.0);
    desc.add_g2(0.5, 0.0);
    desc.add_g2(1.0, 0.0);
    desc.add_g4(0.1, 1.0, 1.0);

    // 3→4→1 with tanh
    let l1 = NeuralNetworkLayer::new(
        vec![
            vec![0.1, -0.2, 0.15],
            vec![0.3, 0.1, -0.1],
            vec![-0.05, 0.25, 0.2],
            vec![0.2, -0.15, 0.1],
        ],
        vec![0.0; 4],
        tanh_activation,
    );
    let l2 = NeuralNetworkLayer::new(
        vec![vec![0.1, 0.2, -0.1, 0.3]],
        vec![0.0],
        identity_activation,
    );
    let pot = BpNeuralNetworkPotential::new(vec![l1, l2], desc);

    let (_, forces, _) = pot.energy_forces_param_grads(&positions, &species);

    // Compare forces[0] with finite diff
    let h = 1e-5_f64;
    let tol = 1e-4_f64; // slightly looser for multi-body
    let pos_arr: Vec<[f64; 3]> = positions.iter().map(|p| [p[0], p[1], p[2]]).collect();
    let neighbors: Vec<usize> = (0..4).collect();

    for alpha in 0..3 {
        let mut pos_plus = pos_arr.clone();
        pos_plus[0][alpha] += h;
        let e_plus: f64 = (0..4)
            .map(|c| {
                let d = pot.descriptor.compute(&pos_plus, c, &neighbors);
                pot.energy(&d)
            })
            .sum();

        let mut pos_minus = pos_arr.clone();
        pos_minus[0][alpha] -= h;
        let e_minus: f64 = (0..4)
            .map(|c| {
                let d = pot.descriptor.compute(&pos_minus, c, &neighbors);
                pot.energy(&d)
            })
            .sum();

        let fd_force = -(e_plus - e_minus) / (2.0 * h);
        let analytical_force = forces[0][alpha];
        assert!(
            (analytical_force - fd_force).abs() < tol,
            "forces[0][{alpha}]: analytical={analytical_force:.8}, fd={fd_force:.8}, diff={:.2e}",
            (analytical_force - fd_force).abs()
        );
    }
}

// ---------------------------------------------------------------------------
// Test 6: apply_param_step changes energy
// ---------------------------------------------------------------------------

#[test]
fn test_apply_param_step_changes_energy() {
    let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
    let species = vec![1u32, 1];

    // Simple 1G2, [1→2→1] network
    let mut desc = BehlerParrinelloDescriptor::new(5.0);
    desc.add_g2(0.5, 0.0);
    let l1 = NeuralNetworkLayer::new(vec![vec![0.3], vec![-0.2]], vec![0.0, 0.0], tanh_activation);
    let l2 = NeuralNetworkLayer::new(vec![vec![0.5, 0.5]], vec![0.0], identity_activation);
    let mut pot = BpNeuralNetworkPotential::new(vec![l1, l2], desc);

    // Initial energy
    let (e0, _, _) = pot.energy_forces_param_grads(&positions, &species);

    // Apply a non-zero parameter step
    let n = pot.n_params();
    let step: Vec<f64> = vec![1e-3; n];
    pot.apply_param_step(&step);

    // New energy should differ
    let (e1, _, _) = pot.energy_forces_param_grads(&positions, &species);

    assert!(
        (e1 - e0).abs() > 1e-10,
        "Energy should change after param step: e0={e0}, e1={e1}"
    );
}

// ---------------------------------------------------------------------------
// Test 7: neighbour-list forces on a 30-atom cluster
// ---------------------------------------------------------------------------

/// Deterministic generator of uniform numbers in `[0, 1)`.
struct Lcg(u64);

impl Lcg {
    fn next_unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn next_in(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }
}

/// `n` points in the cube `[0, edge]^3` with pairwise separation at least
/// `min_dist` (rejection sampling).
fn random_cluster(n: usize, edge: f64, min_dist: f64, seed: u64) -> Vec<[f64; 3]> {
    let mut rng = Lcg(seed);
    let mut points: Vec<[f64; 3]> = Vec::with_capacity(n);
    while points.len() < n {
        let candidate = [
            rng.next_in(0.0, edge),
            rng.next_in(0.0, edge),
            rng.next_in(0.0, edge),
        ];
        let too_close = points.iter().any(|p| {
            let d2: f64 = (0..3).map(|a| (p[a] - candidate[a]).powi(2)).sum();
            d2 < min_dist * min_dist
        });
        if !too_close {
            points.push(candidate);
        }
    }
    points
}

/// BP potential with 3 G2 + 5 G4 symmetry functions (lambda = +-1, several
/// zeta) and a deterministic 8 -> 10 -> 1 network.
fn cluster_potential(r_cut: f64) -> BpNeuralNetworkPotential {
    let mut desc = BehlerParrinelloDescriptor::new(r_cut);
    desc.add_g2(0.3, 0.0);
    desc.add_g2(0.5, 1.5);
    desc.add_g2(1.0, 2.0);
    desc.add_g4(0.05, 1.0, 1.0);
    desc.add_g4(0.02, 2.0, 1.0);
    desc.add_g4(0.1, 1.0, -1.0);
    desc.add_g4(0.03, 4.0, -1.0);
    desc.add_g4(0.2, 2.0, 1.0);
    let n_in = desc.descriptor_len();
    let mut rng = Lcg(2024);
    let hidden = 10;
    let w1: Vec<Vec<f64>> = (0..hidden)
        .map(|_| (0..n_in).map(|_| rng.next_in(-0.4, 0.4)).collect())
        .collect();
    let b1: Vec<f64> = (0..hidden).map(|_| rng.next_in(-0.2, 0.2)).collect();
    let w2 = vec![(0..hidden).map(|_| rng.next_in(-0.5, 0.5)).collect()];
    let l1 = NeuralNetworkLayer::new(w1, b1, tanh_activation);
    let l2 = NeuralNetworkLayer::new(w2, vec![0.05], identity_activation);
    BpNeuralNetworkPotential::new(vec![l1, l2], desc)
}

/// Energy and forces with every atom treated as a neighbour of every centre
/// (the all-pairs evaluation that the cell-list path must reproduce).
fn all_atom_reference(
    pot: &BpNeuralNetworkPotential,
    positions: &[[f64; 3]],
) -> (f64, Vec<[f64; 3]>) {
    let n = positions.len();
    let everyone: Vec<usize> = (0..n).collect();
    let mut energy = 0.0;
    let mut forces = vec![[0.0_f64; 3]; n];
    for center in 0..n {
        let desc = pot.descriptor.compute(positions, center, &everyone);
        let (e_atom, tape) = nn_forward_with_tape(&pot.layers, &desc);
        energy += e_atom;
        let grads = nn_backward(&pot.layers, &tape);
        bp_descriptor_forces(
            &pot.descriptor,
            &grads.ddesc,
            positions,
            center,
            &everyone,
            &mut forces,
        );
    }
    (energy, forces)
}

fn to_vec3(positions: &[[f64; 3]]) -> Vec<Vec3> {
    positions
        .iter()
        .map(|p| Vec3::new(p[0], p[1], p[2]))
        .collect()
}

#[test]
fn test_bp_cell_list_forces_equal_all_atom_forces_on_30_atoms() {
    let positions = random_cluster(30, 7.0, 0.9, 11);
    let pot = cluster_potential(3.0);
    let species = vec![1u32; positions.len()];

    // The cluster is wider than the cutoff, so many pairs lie outside it.
    let far_pairs = (0..30)
        .flat_map(|i| (i + 1..30).map(move |j| (i, j)))
        .filter(|&(i, j)| {
            let d2: f64 = (0..3)
                .map(|a| (positions[i][a] - positions[j][a]).powi(2))
                .sum();
            d2 >= 9.0
        })
        .count();
    assert!(far_pairs > 100, "far pairs: {far_pairs}");

    let (energy, forces, _) = pot.energy_forces_param_grads(&to_vec3(&positions), &species);
    let (ref_energy, ref_forces) = all_atom_reference(&pot, &positions);

    assert!(
        (energy - ref_energy).abs() < 1e-12 * ref_energy.abs().max(1.0),
        "energy {energy} vs {ref_energy}"
    );
    let mut largest = 0.0_f64;
    for (atom, (f, r)) in forces.iter().zip(&ref_forces).enumerate() {
        for a in 0..3 {
            largest = largest.max(r[a].abs());
            assert!(
                (f[a] - r[a]).abs() < 1e-12,
                "atom {atom} component {a}: {} vs {}",
                f[a],
                r[a]
            );
        }
    }
    assert!(largest > 1e-4, "forces must be non-trivial, max {largest}");
}

#[test]
fn test_bp_forces_on_30_atoms_match_finite_differences_of_the_energy() {
    // The analytic forces (G2, G4 with several lambda/zeta/eta sets) against
    // central differences of the independent all-pairs energy.
    let positions = random_cluster(30, 7.0, 0.9, 11);
    let pot = cluster_potential(3.0);
    let species = vec![1u32; positions.len()];
    let (_, forces, _) = pot.energy_forces_param_grads(&to_vec3(&positions), &species);
    let everyone: Vec<usize> = (0..positions.len()).collect();
    let energy_at = |pos: &[[f64; 3]]| -> f64 {
        (0..pos.len())
            .map(|c| pot.energy(&pot.descriptor.compute(pos, c, &everyone)))
            .sum()
    };
    let h = 1e-5;
    for atom in [0usize, 7, 13, 21, 29] {
        for alpha in 0..3 {
            let mut plus = positions.clone();
            plus[atom][alpha] += h;
            let mut minus = positions.clone();
            minus[atom][alpha] -= h;
            let numeric = -(energy_at(&plus) - energy_at(&minus)) / (2.0 * h);
            assert!(
                (forces[atom][alpha] - numeric).abs() < 1e-7,
                "atom {atom} component {alpha}: analytic {} vs numeric {numeric}",
                forces[atom][alpha]
            );
        }
    }
}

#[test]
fn test_bp_forces_sum_to_zero_by_translation_invariance() {
    let positions = random_cluster(30, 7.0, 0.9, 5);
    let pot = cluster_potential(3.0);
    let species = vec![1u32; positions.len()];
    let (_, forces, _) = pot.energy_forces_param_grads(&to_vec3(&positions), &species);
    let mut total = [0.0_f64; 3];
    for f in &forces {
        for a in 0..3 {
            total[a] += f[a];
        }
    }
    for t in total {
        assert!(t.abs() < 1e-10, "net force {total:?}");
    }
}

#[test]
fn test_bp_cell_list_handles_a_long_chain_of_atoms() {
    // 400 atoms 1.3 A apart with r_cut = 3: every atom has 4 neighbours, so
    // the evaluation is linear in the number of atoms.  Spot-check forces
    // against finite differences of the (cell-list) energy, and the atoms
    // beyond the cutoff of an end atom must not influence it.
    let n = 400;
    let positions: Vec<[f64; 3]> = (0..n)
        .map(|i| {
            let s = i as f64;
            [1.3 * s, 0.25 * (0.7 * s).sin(), 0.2 * (0.45 * s).cos()]
        })
        .collect();
    let pot = cluster_potential(3.0);
    let species = vec![1u32; n];
    let (energy, forces, _) = pot.energy_forces_param_grads(&to_vec3(&positions), &species);
    assert!(energy.is_finite());
    assert!(forces.iter().all(|f| f.iter().all(|c| c.is_finite())));
    let net: Vec3 = forces.iter().fold(Vec3::zeros(), |acc, f| acc + f);
    assert!(net.norm() < 1e-9, "net force {net:?}");

    let energy_of =
        |pos: &[[f64; 3]]| -> f64 { pot.energy_forces_param_grads(&to_vec3(pos), &species).0 };
    let h = 1e-5;
    for atom in [0usize, 1, 150, 399] {
        for alpha in 0..3 {
            let mut plus = positions.clone();
            plus[atom][alpha] += h;
            let mut minus = positions.clone();
            minus[atom][alpha] -= h;
            let numeric = -(energy_of(&plus) - energy_of(&minus)) / (2.0 * h);
            assert!(
                (forces[atom][alpha] - numeric).abs() < 1e-7,
                "atom {atom} component {alpha}: analytic {} vs numeric {numeric}",
                forces[atom][alpha]
            );
        }
    }

    // Moving an atom that is farther than r_cut + r_cut from atom 0 leaves
    // the force on atom 0 bit-for-bit unchanged.
    let mut moved = positions.clone();
    moved[300][1] += 0.4;
    let (_, moved_forces, _) = pot.energy_forces_param_grads(&to_vec3(&moved), &species);
    assert_eq!(moved_forces[0], forces[0]);
}

#[test]
fn test_bp_descriptor_forces_accumulate_into_the_caller_buffer() {
    let positions = random_cluster(8, 3.5, 0.9, 3);
    let pot = cluster_potential(3.0);
    let everyone: Vec<usize> = (0..8).collect();
    let desc = pot.descriptor.compute(&positions, 2, &everyone);
    let (_, tape) = nn_forward_with_tape(&pot.layers, &desc);
    let grads = nn_backward(&pot.layers, &tape);

    let mut once = vec![[0.0_f64; 3]; 8];
    bp_descriptor_forces(
        &pot.descriptor,
        &grads.ddesc,
        &positions,
        2,
        &everyone,
        &mut once,
    );
    let mut twice = once.clone();
    bp_descriptor_forces(
        &pot.descriptor,
        &grads.ddesc,
        &positions,
        2,
        &everyone,
        &mut twice,
    );
    for (a, b) in once.iter().zip(&twice) {
        for c in 0..3 {
            assert!((2.0 * a[c] - b[c]).abs() < 1e-14, "{a:?} vs {b:?}");
        }
    }
    assert!(once.iter().any(|f| f.iter().any(|c| c.abs() > 1e-6)));

    // The centre itself, duplicated far atoms and out-of-range indices are
    // ignored: an unfiltered list with extras gives the same forces.
    let mut with_extras = vec![[0.0_f64; 3]; 8];
    let mut list = everyone.clone();
    list.push(99);
    bp_descriptor_forces(
        &pot.descriptor,
        &grads.ddesc,
        &positions,
        2,
        &list,
        &mut with_extras,
    );
    assert_eq!(with_extras, once);

    // Out-of-range centre or a short buffer: nothing happens.
    let mut untouched = vec![[0.0_f64; 3]; 8];
    bp_descriptor_forces(
        &pot.descriptor,
        &grads.ddesc,
        &positions,
        8,
        &everyone,
        &mut untouched,
    );
    let mut short = vec![[0.0_f64; 3]; 4];
    bp_descriptor_forces(
        &pot.descriptor,
        &grads.ddesc,
        &positions,
        2,
        &everyone,
        &mut short,
    );
    assert!(untouched.iter().chain(short.iter()).all(|f| *f == [0.0; 3]));
}
