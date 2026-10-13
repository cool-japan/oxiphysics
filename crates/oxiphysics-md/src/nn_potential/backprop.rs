// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Reverse-mode differentiation of a stack of [`DenseLayer`]s.
//!
//! A forward pass keeps every layer's output ([`NetworkTrace`]); the backward
//! pass then turns the gradient with respect to the network output into the
//! gradient with respect to its input.  The derivative of each activation is
//! expressed through the stored layer output
//! ([`Activation::derivative_from_output`](super::types::Activation::derivative_from_output)),
//! so no pre-activations need to be cached.

use super::types::DenseLayer;

/// Layer outputs of one forward pass through a network.
///
/// `activations[0]` is the network input and `activations[l + 1]` the output
/// of layer `l`; the last entry is the network output.
#[derive(Debug, Clone)]
pub(super) struct NetworkTrace {
    activations: Vec<Vec<f64>>,
}

impl NetworkTrace {
    /// Run `input` through `layers`, recording every layer output.
    pub(super) fn forward(layers: &[DenseLayer], input: &[f64]) -> Self {
        let mut activations = Vec::with_capacity(layers.len() + 1);
        let mut current = input.to_vec();
        for layer in layers {
            let next = layer.forward(&current);
            activations.push(std::mem::replace(&mut current, next));
        }
        activations.push(current);
        Self { activations }
    }

    /// The network output (the input itself for an empty layer stack).
    pub(super) fn output(&self) -> &[f64] {
        self.activations.last().map(Vec::as_slice).unwrap_or(&[])
    }

    /// Gradient with respect to the network input, given the gradient
    /// `grad_output` with respect to the network output.
    pub(super) fn backward(&self, layers: &[DenseLayer], grad_output: &[f64]) -> Vec<f64> {
        let mut grad = grad_output.to_vec();
        for (l, layer) in layers.iter().enumerate().rev() {
            let output = &self.activations[l + 1];
            let input = &self.activations[l];
            let mut grad_in = vec![0.0; input.len()];
            for (i, (&y, &g)) in output.iter().zip(grad.iter()).enumerate() {
                let delta = g * layer.activation.derivative_from_output(y);
                for (gj, &w) in grad_in.iter_mut().zip(layer.weights[i].iter()) {
                    *gj += w * delta;
                }
            }
            grad = grad_in;
        }
        grad
    }
}

#[cfg(test)]
mod tests {
    use super::super::types::{Activation, DenseLayer};
    use super::*;

    fn test_network() -> Vec<DenseLayer> {
        vec![
            DenseLayer::new(
                vec![vec![0.3, -0.5, 0.2], vec![0.1, 0.4, -0.6]],
                vec![0.05, -0.1],
                Activation::Tanh,
            ),
            DenseLayer::new(
                vec![vec![0.7, -0.2], vec![-0.4, 0.9], vec![0.25, 0.35]],
                vec![0.0, 0.1, -0.2],
                Activation::Sigmoid,
            ),
            DenseLayer::new(vec![vec![1.1, -0.3, 0.6]], vec![0.02], Activation::Linear),
        ]
    }

    fn scalar_output(layers: &[DenseLayer], x: &[f64]) -> f64 {
        NetworkTrace::forward(layers, x).output().iter().sum()
    }

    #[test]
    fn backward_matches_central_differences_for_every_activation() {
        let layers = test_network();
        let x = [0.4, -0.7, 1.2];
        let trace = NetworkTrace::forward(&layers, &x);
        let grad = trace.backward(&layers, &[1.0]);
        assert_eq!(grad.len(), 3);
        let h = 1e-6;
        for k in 0..3 {
            let mut plus = x;
            let mut minus = x;
            plus[k] += h;
            minus[k] -= h;
            let numeric =
                (scalar_output(&layers, &plus) - scalar_output(&layers, &minus)) / (2.0 * h);
            assert!(
                (grad[k] - numeric).abs() < 1e-8,
                "input {k}: analytic {} vs numeric {numeric}",
                grad[k]
            );
        }
    }

    #[test]
    fn relu_uses_the_zero_subgradient_below_the_kink() {
        let layers = vec![DenseLayer::new(
            vec![vec![1.0, 1.0], vec![-1.0, -1.0]],
            vec![0.0, 0.0],
            Activation::ReLU,
        )];
        let trace = NetworkTrace::forward(&layers, &[0.5, 0.25]);
        assert_eq!(trace.output(), &[0.75, 0.0]);
        // Only the first (active) unit passes gradient back.
        let grad = trace.backward(&layers, &[1.0, 1.0]);
        assert_eq!(grad, vec![1.0, 1.0]);
    }

    #[test]
    fn vector_seed_weights_each_output_unit() {
        let layers = vec![DenseLayer::new(
            vec![vec![2.0, 0.0], vec![0.0, 3.0]],
            vec![0.0, 0.0],
            Activation::Linear,
        )];
        let trace = NetworkTrace::forward(&layers, &[1.0, 1.0]);
        assert_eq!(trace.backward(&layers, &[10.0, 100.0]), vec![20.0, 300.0]);
    }

    #[test]
    fn empty_stack_is_the_identity() {
        let trace = NetworkTrace::forward(&[], &[1.0, 2.0]);
        assert_eq!(trace.output(), &[1.0, 2.0]);
        assert_eq!(trace.backward(&[], &[5.0, 6.0]), vec![5.0, 6.0]);
    }
}
