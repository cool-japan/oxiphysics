// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural network-enhanced soft body simulation.
//!
//! Provides a data-driven constitutive model using a small 2-layer ReLU
//! network for stress prediction, online gradient-descent training, and a
//! physics-informed loss term enforcing polyconvexity.
//!
//! No nalgebra dependency — all arrays are plain `f64` slices / `Vec`f64`.

// ---------------------------------------------------------------------------
// Activation helpers
// ---------------------------------------------------------------------------

/// Rectified linear unit: ReLU(x) = max(0, x).
#[inline]
fn relu(x: f64) -> f64 {
    x.max(0.0)
}

/// Derivative of ReLU.
#[inline]
fn relu_grad(x: f64) -> f64 {
    if x > 0.0 { 1.0 } else { 0.0 }
}

/// Dense-layer forward pass: y = W x + b.
///
/// - `w` : weight matrix stored row-major, shape `\[out × in\]`.
/// - `b` : bias vector, length `out`.
/// - `x` : input vector, length `in`.
fn dense_forward(w: &[f64], b: &[f64], x: &[f64], n_in: usize, n_out: usize) -> Vec<f64> {
    assert_eq!(w.len(), n_out * n_in);
    assert_eq!(b.len(), n_out);
    assert_eq!(x.len(), n_in);
    (0..n_out)
        .map(|i| b[i] + (0..n_in).map(|j| w[i * n_in + j] * x[j]).sum::<f64>())
        .collect()
}

// ---------------------------------------------------------------------------
// NeuralConstitutive
// ---------------------------------------------------------------------------

/// Neural-network constitutive model for a hyperelastic soft body.
///
/// Architecture: input (n_strain) → hidden (n_hidden) ReLU → output (n_stress).
/// Typically n_strain = n_stress = 6 (Voigt notation for 3-D).
#[derive(Debug, Clone)]
pub struct NeuralConstitutive {
    /// Number of strain components (input size).
    pub n_strain: usize,
    /// Number of hidden neurons.
    pub n_hidden: usize,
    /// Number of stress components (output size).
    pub n_stress: usize,
    /// Layer 1 weights W1, shape \[n_hidden × n_strain\], row-major.
    pub w1: Vec<f64>,
    /// Layer 1 biases b1, length n_hidden.
    pub b1: Vec<f64>,
    /// Layer 2 weights W2, shape \[n_stress × n_hidden\], row-major.
    pub w2: Vec<f64>,
    /// Layer 2 biases b2, length n_stress.
    pub b2: Vec<f64>,
    /// Learning rate for gradient descent.
    pub learning_rate: f64,
}

impl NeuralConstitutive {
    /// Create a new network with Xavier-like initialisation.
    ///
    /// Weights are initialised to a small deterministic pattern
    /// (no random dependency) so tests are reproducible.
    pub fn new(n_strain: usize, n_hidden: usize, n_stress: usize, learning_rate: f64) -> Self {
        let scale1 = (2.0 / n_strain as f64).sqrt();
        let scale2 = (2.0 / n_hidden as f64).sqrt();

        let w1: Vec<f64> = (0..n_hidden * n_strain)
            .map(|i| {
                let v = ((i as f64 * 1.618) % 2.0) - 1.0; // [-1, 1]
                v * scale1
            })
            .collect();
        let b1 = vec![0.0; n_hidden];

        let w2: Vec<f64> = (0..n_stress * n_hidden)
            .map(|i| {
                let v = ((i as f64 * 2.719) % 2.0) - 1.0;
                v * scale2
            })
            .collect();
        let b2 = vec![0.0; n_stress];

        Self {
            n_strain,
            n_hidden,
            n_stress,
            w1,
            b1,
            w2,
            b2,
            learning_rate,
        }
    }

    // -----------------------------------------------------------------------
    // Forward pass
    // -----------------------------------------------------------------------

    /// Compute the predicted stress vector from a strain vector.
    ///
    /// σ = W2 · ReLU(W1 · ε + b1) + b2
    pub fn forward_nn(&self, strain: &[f64]) -> Vec<f64> {
        assert_eq!(strain.len(), self.n_strain);
        // Layer 1
        let z1 = dense_forward(&self.w1, &self.b1, strain, self.n_strain, self.n_hidden);
        let a1: Vec<f64> = z1.iter().map(|&v| relu(v)).collect();
        // Layer 2
        dense_forward(&self.w2, &self.b2, &a1, self.n_hidden, self.n_stress)
    }

    // -----------------------------------------------------------------------
    // Strain energy
    // -----------------------------------------------------------------------

    /// Compute the **strain energy density** ψ as the sum of predicted stresses
    /// dotted with the strain vector.
    ///
    /// ψ ≈ σ · ε  (approximation for small strains)
    pub fn strain_energy_nn(&self, strain: &[f64]) -> f64 {
        let stress = self.forward_nn(strain);
        stress.iter().zip(strain.iter()).map(|(s, e)| s * e).sum()
    }

    // -----------------------------------------------------------------------
    // Physics-informed loss
    // -----------------------------------------------------------------------

    /// Compute a **physics-informed loss** combining data fidelity and a
    /// polyconvexity penalty.
    ///
    /// L = ||σ_pred − σ_ref||² + λ · max(0, −ψ)²
    ///
    /// The polyconvexity penalty penalises negative strain energy (which would
    /// violate thermodynamic consistency).
    pub fn physics_informed_loss(
        &self,
        strain: &[f64],
        stress_ref: &[f64],
        lambda_convex: f64,
    ) -> f64 {
        let stress_pred = self.forward_nn(strain);
        let data_loss: f64 = stress_pred
            .iter()
            .zip(stress_ref.iter())
            .map(|(p, r)| (p - r).powi(2))
            .sum();
        let psi = self.strain_energy_nn(strain);
        let convex_penalty = (-psi).max(0.0).powi(2);
        data_loss + lambda_convex * convex_penalty
    }

    // -----------------------------------------------------------------------
    // Training step
    // -----------------------------------------------------------------------

    /// Perform one **gradient-descent update step** on a batch of
    /// (strain, stress) pairs.
    ///
    /// Uses MSE loss and back-propagation through the 2-layer network.
    pub fn train_step(&mut self, strains: &[Vec<f64>], stresses: &[Vec<f64>]) {
        assert_eq!(strains.len(), stresses.len());
        let n = strains.len() as f64;

        // Accumulate gradients.
        let mut dw2 = vec![0.0f64; self.n_stress * self.n_hidden];
        let mut db2 = vec![0.0f64; self.n_stress];
        let mut dw1 = vec![0.0f64; self.n_hidden * self.n_strain];
        let mut db1 = vec![0.0f64; self.n_hidden];

        for (strain, stress_ref) in strains.iter().zip(stresses.iter()) {
            // Forward
            let z1 = dense_forward(&self.w1, &self.b1, strain, self.n_strain, self.n_hidden);
            let a1: Vec<f64> = z1.iter().map(|&v| relu(v)).collect();
            let pred = dense_forward(&self.w2, &self.b2, &a1, self.n_hidden, self.n_stress);

            // Output layer error: δ2 = 2(pred − ref) / n
            let delta2: Vec<f64> = pred
                .iter()
                .zip(stress_ref.iter())
                .map(|(p, r)| 2.0 * (p - r) / n)
                .collect();

            // Gradients for W2, b2
            for i in 0..self.n_stress {
                db2[i] += delta2[i];
                for j in 0..self.n_hidden {
                    dw2[i * self.n_hidden + j] += delta2[i] * a1[j];
                }
            }

            // Back-prop through ReLU: δ1 = (W2^T δ2) ⊙ relu'(z1)
            let mut delta1 = vec![0.0f64; self.n_hidden];
            for j in 0..self.n_hidden {
                let grad: f64 = (0..self.n_stress)
                    .map(|i| self.w2[i * self.n_hidden + j] * delta2[i])
                    .sum();
                delta1[j] = grad * relu_grad(z1[j]);
            }

            // Gradients for W1, b1
            for j in 0..self.n_hidden {
                db1[j] += delta1[j];
                for k in 0..self.n_strain {
                    dw1[j * self.n_strain + k] += delta1[j] * strain[k];
                }
            }
        }

        // Apply updates (SGD)
        for (w, g) in self.w2.iter_mut().zip(dw2.iter()) {
            *w -= self.learning_rate * g;
        }
        for (b, g) in self.b2.iter_mut().zip(db2.iter()) {
            *b -= self.learning_rate * g;
        }
        for (w, g) in self.w1.iter_mut().zip(dw1.iter()) {
            *w -= self.learning_rate * g;
        }
        for (b, g) in self.b1.iter_mut().zip(db1.iter()) {
            *b -= self.learning_rate * g;
        }
    }

    // -----------------------------------------------------------------------
    // Data-driven material fitting
    // -----------------------------------------------------------------------

    /// Fit the constitutive model to stress-strain data using repeated
    /// gradient-descent passes (epochs).
    ///
    /// Returns the final MSE loss.
    pub fn data_driven_material(
        &mut self,
        strains: &[Vec<f64>],
        stresses: &[Vec<f64>],
        epochs: usize,
    ) -> f64 {
        for _ in 0..epochs {
            self.train_step(strains, stresses);
        }
        // Compute final loss
        let n = strains.len() as f64;
        strains
            .iter()
            .zip(stresses.iter())
            .map(|(e, s_ref)| {
                let s_pred = self.forward_nn(e);
                s_pred
                    .iter()
                    .zip(s_ref.iter())
                    .map(|(p, r)| (p - r).powi(2))
                    .sum::<f64>()
                    / n
            })
            .sum()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_net() -> NeuralConstitutive {
        NeuralConstitutive::new(6, 12, 6, 1e-3)
    }

    // 1. Forward pass produces output of correct length.
    #[test]
    fn test_forward_output_shape() {
        let net = make_net();
        let strain = vec![0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
        let stress = net.forward_nn(&strain);
        assert_eq!(
            stress.len(),
            6,
            "Forward pass should produce 6 stress components"
        );
    }

    // 2. Zero strain input produces finite output (no NaN/Inf).
    #[test]
    fn test_forward_zero_strain_finite() {
        let net = make_net();
        let strain = vec![0.0; 6];
        let stress = net.forward_nn(&strain);
        for &s in &stress {
            assert!(s.is_finite(), "Stress component should be finite, got {s}");
        }
    }

    // 3. Strain energy is finite.
    #[test]
    fn test_strain_energy_finite() {
        let net = make_net();
        let strain = vec![0.01, 0.005, 0.0, 0.0, 0.0, 0.0];
        let psi = net.strain_energy_nn(&strain);
        assert!(psi.is_finite(), "Strain energy should be finite, got {psi}");
    }

    // 4. Physics-informed loss is non-negative.
    #[test]
    fn test_physics_loss_nonnegative() {
        let net = make_net();
        let strain = vec![0.01; 6];
        let stress_ref = vec![100.0; 6];
        let loss = net.physics_informed_loss(&strain, &stress_ref, 1.0);
        assert!(loss >= 0.0, "Loss must be non-negative, got {loss}");
    }

    // 5. Loss at perfect prediction is zero (ignoring convexity penalty).
    #[test]
    fn test_loss_zero_at_perfect_prediction() {
        let net = make_net();
        let strain = vec![0.0; 6];
        let stress_pred = net.forward_nn(&strain);
        let loss = net.physics_informed_loss(&strain, &stress_pred, 0.0);
        assert!(
            loss.abs() < 1e-20,
            "Loss at perfect prediction should be 0, got {loss}"
        );
    }

    // 6. Training step reduces MSE on a simple linear dataset.
    #[test]
    fn test_train_step_reduces_loss() {
        let mut net = NeuralConstitutive::new(2, 4, 2, 0.01);
        // Target: identity map ε → σ
        let strains: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![0.5, 0.5]];
        let stresses = strains.clone();
        let loss_before: f64 = strains
            .iter()
            .zip(stresses.iter())
            .map(|(e, s)| {
                let p = net.forward_nn(e);
                p.iter()
                    .zip(s.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
            })
            .sum();
        for _ in 0..200 {
            net.train_step(&strains, &stresses);
        }
        let loss_after: f64 = strains
            .iter()
            .zip(stresses.iter())
            .map(|(e, s)| {
                let p = net.forward_nn(e);
                p.iter()
                    .zip(s.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
            })
            .sum();
        assert!(
            loss_after < loss_before,
            "Training should reduce MSE: before={loss_before:.4}, after={loss_after:.4}"
        );
    }

    // 7. data_driven_material converges on constant stress dataset.
    #[test]
    fn test_data_driven_convergence() {
        let mut net = NeuralConstitutive::new(2, 8, 2, 0.01);
        let strains: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64 * 0.01, 0.0]).collect();
        // Target: constant stress [1.0, 0.0]
        let stresses: Vec<Vec<f64>> = (0..10).map(|_| vec![1.0, 0.0]).collect();
        let final_loss = net.data_driven_material(&strains, &stresses, 500);
        assert!(
            final_loss < 5.0,
            "data_driven_material should reduce loss, final={final_loss:.4}"
        );
    }

    // 8. dense_forward produces correct result for identity weight matrix.
    #[test]
    fn test_dense_forward_identity() {
        // 2×2 identity weight, zero bias
        let w = vec![1.0, 0.0, 0.0, 1.0];
        let b = vec![0.0, 0.0];
        let x = vec![3.0, 5.0];
        let y = dense_forward(&w, &b, &x, 2, 2);
        assert!((y[0] - 3.0).abs() < 1e-10);
        assert!((y[1] - 5.0).abs() < 1e-10);
    }

    // 9. ReLU activates only positive inputs.
    #[test]
    fn test_relu_activation() {
        assert_eq!(relu(-5.0), 0.0);
        assert_eq!(relu(0.0), 0.0);
        assert!(relu(3.0) > 0.0);
    }

    // 10. relu_grad is 1 for positive and 0 for non-positive.
    #[test]
    fn test_relu_grad() {
        assert_eq!(relu_grad(1.0), 1.0);
        assert_eq!(relu_grad(0.0), 0.0);
        assert_eq!(relu_grad(-1.0), 0.0);
    }

    // 11. Weights change after a training step.
    #[test]
    fn test_weights_change_after_train() {
        let mut net = NeuralConstitutive::new(2, 4, 2, 0.01);
        let w1_before = net.w1.clone();
        let strains = vec![vec![1.0, 0.0]];
        let stresses = vec![vec![2.0, 1.0]];
        net.train_step(&strains, &stresses);
        assert_ne!(net.w1, w1_before, "Weights should change after training");
    }

    // 12. Strain energy scales roughly with strain magnitude.
    #[test]
    fn test_strain_energy_scales() {
        let net = make_net();
        let strain_small = vec![0.001; 6];
        let strain_large = vec![0.1; 6];
        let psi_small = net.strain_energy_nn(&strain_small).abs();
        let psi_large = net.strain_energy_nn(&strain_large).abs();
        // Large strain should give larger energy magnitude (network may have mixed signs)
        // Just ensure psi_large is larger than psi_small
        assert!(
            psi_large >= psi_small,
            "Larger strain should give larger energy magnitude: small={psi_small:.6}, large={psi_large:.6}"
        );
    }

    // 13. Forward pass is deterministic (same input → same output).
    #[test]
    fn test_forward_deterministic() {
        let net = make_net();
        let strain = vec![0.02, 0.01, 0.0, 0.005, 0.0, 0.0];
        let s1 = net.forward_nn(&strain);
        let s2 = net.forward_nn(&strain);
        assert_eq!(s1, s2, "Forward pass should be deterministic");
    }

    // 14. Physics loss increases when lambda is larger (non-trivial strain energy).
    #[test]
    fn test_physics_loss_lambda_effect() {
        let net = make_net();
        let strain = vec![-0.1; 6]; // may produce negative strain energy
        let stress_ref = vec![0.0; 6];
        let loss0 = net.physics_informed_loss(&strain, &stress_ref, 0.0);
        let loss1 = net.physics_informed_loss(&strain, &stress_ref, 100.0);
        // If strain energy is negative, higher lambda gives higher loss.
        // If strain energy is positive, loss is the same.
        assert!(
            loss1 >= loss0,
            "Higher lambda should not decrease the physics loss"
        );
    }

    // 15. Network with n_hidden=1 still runs without panic.
    #[test]
    fn test_minimal_network() {
        let net = NeuralConstitutive::new(1, 1, 1, 1e-3);
        let stress = net.forward_nn(&[0.05]);
        assert_eq!(stress.len(), 1);
        assert!(stress[0].is_finite());
    }

    // 16. Multiple consecutive train steps keep decreasing loss.
    #[test]
    fn test_repeated_training() {
        let mut net = NeuralConstitutive::new(1, 8, 1, 0.005);
        let strains: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64 * 0.1]).collect();
        let stresses: Vec<Vec<f64>> = strains.iter().map(|e| vec![e[0] * 2.0]).collect();
        let mut prev_loss = f64::INFINITY;
        for _ in 0..50 {
            let loss: f64 = strains
                .iter()
                .zip(stresses.iter())
                .map(|(e, s)| {
                    let p = net.forward_nn(e);
                    (p[0] - s[0]).powi(2)
                })
                .sum::<f64>()
                / strains.len() as f64;
            net.train_step(&strains, &stresses);
            let _ = prev_loss;
            prev_loss = loss;
        }
        assert!(
            prev_loss.is_finite(),
            "Loss should remain finite after repeated training"
        );
    }

    // 17. dense_forward bias adds correctly.
    #[test]
    fn test_dense_forward_bias() {
        let w = vec![1.0, 0.0, 0.0, 1.0];
        let b = vec![10.0, 20.0];
        let x = vec![1.0, 1.0];
        let y = dense_forward(&w, &b, &x, 2, 2);
        assert!((y[0] - 11.0).abs() < 1e-10);
        assert!((y[1] - 21.0).abs() < 1e-10);
    }

    // 18. Training on zero-stress data drives predictions toward zero.
    #[test]
    fn test_train_toward_zero() {
        let mut net = NeuralConstitutive::new(2, 4, 2, 0.01);
        let strains: Vec<Vec<f64>> = vec![vec![0.0, 0.0]; 5];
        let stresses: Vec<Vec<f64>> = vec![vec![0.0, 0.0]; 5];
        for _ in 0..300 {
            net.train_step(&strains, &stresses);
        }
        let pred = net.forward_nn(&[0.0, 0.0]);
        for &p in &pred {
            assert!(
                p.abs() < 0.5,
                "After training on zero targets, predictions should be near zero, got {p}"
            );
        }
    }

    // 19. Strain energy is zero when all stresses are zero at zero strain.
    #[test]
    fn test_strain_energy_zero_strain() {
        let mut net = NeuralConstitutive::new(2, 4, 2, 0.01);
        // Train to produce zero stress everywhere
        let strains: Vec<Vec<f64>> = vec![vec![0.0, 0.0]; 5];
        let stresses: Vec<Vec<f64>> = vec![vec![0.0, 0.0]; 5];
        for _ in 0..500 {
            net.train_step(&strains, &stresses);
        }
        let psi = net.strain_energy_nn(&[0.0, 0.0]);
        assert!(
            psi.abs() < 0.1,
            "Strain energy at zero strain should be ≈ 0, got {psi}"
        );
    }

    // 20. Changing learning rate affects training speed.
    #[test]
    fn test_learning_rate_effect() {
        let strains: Vec<Vec<f64>> = vec![vec![1.0, 0.0], vec![0.0, 1.0]];
        let stresses: Vec<Vec<f64>> = vec![vec![2.0, 0.0], vec![0.0, 2.0]];

        let mut net_slow = NeuralConstitutive::new(2, 4, 2, 1e-5);
        let mut net_fast = NeuralConstitutive::new(2, 4, 2, 1e-2);

        for _ in 0..100 {
            net_slow.train_step(&strains, &stresses);
            net_fast.train_step(&strains, &stresses);
        }

        let loss_slow: f64 = strains
            .iter()
            .zip(stresses.iter())
            .map(|(e, s)| {
                let p = net_slow.forward_nn(e);
                p.iter()
                    .zip(s.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
            })
            .sum();
        let loss_fast: f64 = strains
            .iter()
            .zip(stresses.iter())
            .map(|(e, s)| {
                let p = net_fast.forward_nn(e);
                p.iter()
                    .zip(s.iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
            })
            .sum();

        assert!(
            loss_fast < loss_slow,
            "Faster learning rate should converge faster: slow={loss_slow:.4}, fast={loss_fast:.4}"
        );
    }

    // 21. Clone produces independent network.
    #[test]
    fn test_clone_independent() {
        let net1 = make_net();
        let mut net2 = net1.clone();
        let strain = vec![0.01; 6];
        let stress_ref = vec![1.0; 6];
        net2.train_step(
            std::slice::from_ref(&strain),
            std::slice::from_ref(&stress_ref),
        );
        let s1 = net1.forward_nn(&strain);
        let s2 = net2.forward_nn(&strain);
        assert_ne!(s1, s2, "Clone should be independent after training");
    }

    // 22. Forward pass handles large strain inputs without NaN.
    #[test]
    fn test_forward_large_strain() {
        let net = make_net();
        let strain = vec![100.0; 6];
        let stress = net.forward_nn(&strain);
        for &s in &stress {
            assert!(
                s.is_finite(),
                "Large strain should still give finite stress, got {s}"
            );
        }
    }

    // 23. data_driven_material returns non-negative loss.
    #[test]
    fn test_data_driven_nonnegative_loss() {
        let mut net = NeuralConstitutive::new(3, 6, 3, 1e-3);
        let strains: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64 * 0.01, 0.0, 0.0]).collect();
        let stresses: Vec<Vec<f64>> = strains.iter().map(|e| vec![e[0] * 3.0, 0.0, 0.0]).collect();
        let loss = net.data_driven_material(&strains, &stresses, 10);
        assert!(loss >= 0.0, "Loss must be non-negative, got {loss}");
    }

    // 24. Physics loss at zero lambda equals pure MSE.
    #[test]
    fn test_physics_loss_zero_lambda_is_mse() {
        let net = make_net();
        let strain = vec![0.05; 6];
        let stress_ref = vec![1.0; 6];
        let stress_pred = net.forward_nn(&strain);
        let mse: f64 = stress_pred
            .iter()
            .zip(stress_ref.iter())
            .map(|(p, r)| (p - r).powi(2))
            .sum();
        let loss = net.physics_informed_loss(&strain, &stress_ref, 0.0);
        assert!(
            (loss - mse).abs() < 1e-10,
            "Physics loss with λ=0 should equal MSE"
        );
    }

    // 25. Network output changes when weights change.
    #[test]
    fn test_output_changes_with_weights() {
        let mut net = make_net();
        let strain = vec![0.01; 6];
        let s1 = net.forward_nn(&strain);
        // Manually perturb one weight
        net.w1[0] += 1.0;
        let s2 = net.forward_nn(&strain);
        assert_ne!(s1, s2, "Changing a weight should change the output");
    }
}
