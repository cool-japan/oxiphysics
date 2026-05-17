// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Integration tests for the autograd bridge:
//! analytical gradient correctness verified against central finite differences.

use oxiphysics_core::math::Vec3;
use oxiphysics_md::autograd_bridge::{
    DifferentiableForceField, NnGrads, nn_backward, nn_forward_with_tape,
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
