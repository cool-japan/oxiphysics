//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::too_many_arguments)]
use crate::ml_potential::{MlPotential, cutoff_function};
use oxiphysics_core::math::Vec3;

use super::types::{AngularTriplet, EnsembleNnp, FeedForwardPotential, NnAtomisticSystem};

/// G5 angular symmetry function (different radial decay from G4).
///
/// G5 = 2^(1-zeta) * sum_{j,k≠j} (1 + lambda*cos(theta))^zeta
///        * exp(-eta*(r_ij^2 + r_ik^2)) * fc(r_ij) * fc(r_ik)
/// (No r_jk term, unlike G4.)
#[allow(dead_code)]
pub fn compute_g5(
    r_ij: f64,
    r_ik: f64,
    cos_theta: f64,
    eta: f64,
    zeta: f64,
    lambda: f64,
    rc: f64,
) -> f64 {
    use crate::ml_potential::cutoff_function;
    if r_ij >= rc || r_ik >= rc {
        return 0.0;
    }
    let prefactor = 2.0_f64.powf(1.0 - zeta);
    let angular = (1.0 + lambda * cos_theta).powf(zeta);
    let radial = (-eta * (r_ij * r_ij + r_ik * r_ik)).exp();
    let fc = cutoff_function(r_ij, rc) * cutoff_function(r_ik, rc);
    prefactor * angular * radial * fc
}
/// Radial power-type symmetry function (BP G3, less common):
/// G3 = sum_j cos(kappa * r_ij) * fc(r_ij)
#[allow(dead_code)]
pub fn compute_g3(r_ij: f64, kappa: f64, rc: f64) -> f64 {
    if r_ij >= rc {
        return 0.0;
    }
    (kappa * r_ij).cos() * cutoff_function(r_ij, rc)
}
/// Verify energy-force consistency via central finite differences.
///
/// Returns the maximum absolute difference between analytic and numerical
/// forces across all atoms and Cartesian components.
#[allow(dead_code)]
pub fn energy_force_consistency(
    potential: &FeedForwardPotential,
    positions: &[Vec3],
    species: &[u32],
    h: f64,
) -> f64 {
    let (_, forces_analytic) = potential.energy_and_forces(positions, species);
    let n = positions.len();
    let mut max_diff = 0.0_f64;
    let mut pos = positions.to_vec();
    for i in 0..n {
        for a in 0..3_usize {
            pos[i][a] += h;
            let (e_plus, _) = potential.energy_and_forces(&pos, species);
            pos[i][a] -= 2.0 * h;
            let (e_minus, _) = potential.energy_and_forces(&pos, species);
            pos[i][a] += h;
            let f_numeric = -(e_plus - e_minus) / (2.0 * h);
            let diff = (forces_analytic[i][a] - f_numeric).abs();
            if diff > max_diff {
                max_diff = diff;
            }
        }
    }
    max_diff
}
/// Compute Chebyshev polynomial descriptors for a pair distance `r`.
///
/// The distance is first mapped to \[-1, 1\] via `x = 2*(r/rc) - 1`, then
/// the first `n_basis` Chebyshev polynomials T_k(x) are evaluated.
///
/// Returns a vector of length `n_basis`.
#[allow(dead_code)]
pub fn chebyshev_descriptor(r: f64, rc: f64, n_basis: usize) -> Vec<f64> {
    if n_basis == 0 || r >= rc {
        return vec![0.0; n_basis];
    }
    let x = 2.0 * r / rc - 1.0;
    let mut t = vec![0.0; n_basis];
    t[0] = 1.0;
    if n_basis > 1 {
        t[1] = x;
    }
    for k in 2..n_basis {
        t[k] = 2.0 * x * t[k - 1] - t[k - 2];
    }
    let fc = cutoff_function(r, rc);
    t.iter().map(|&tk| tk * fc).collect()
}
/// Build a full pair descriptor: \[r, r/rc\] + Chebyshev expansion of length `n_cheb`.
///
/// Returns a vector of length `2 + n_cheb`.
#[allow(dead_code)]
pub fn pair_descriptor(r: f64, rc: f64, n_cheb: usize) -> Vec<f64> {
    let mut d = vec![r, r / rc];
    d.extend(chebyshev_descriptor(r, rc, n_cheb));
    d
}
/// G4 angular symmetry function (Behler-Parrinello).
///
/// ```text
/// G4 = 2^(1-zeta) * (1 + lambda*cos(theta))^zeta
///        * exp(-eta*(r_ij^2 + r_ik^2 + r_jk^2)) * fc(r_ij) * fc(r_ik) * fc(r_jk)
/// ```
#[allow(dead_code)]
pub fn compute_g4(
    r_ij: f64,
    r_ik: f64,
    r_jk: f64,
    cos_theta: f64,
    eta: f64,
    zeta: f64,
    lambda: f64,
    rc: f64,
) -> f64 {
    if r_ij >= rc || r_ik >= rc || r_jk >= rc {
        return 0.0;
    }
    let prefactor = 2.0_f64.powf(1.0 - zeta);
    let angular = (1.0 + lambda * cos_theta).powf(zeta);
    let radial = (-eta * (r_ij * r_ij + r_ik * r_ik + r_jk * r_jk)).exp();
    let fc = cutoff_function(r_ij, rc) * cutoff_function(r_ik, rc) * cutoff_function(r_jk, rc);
    prefactor * angular * radial * fc
}
/// Select the most uncertain candidate structure from a pool using query-by-committee.
///
/// For each candidate descriptor in `pool`, the disagreement is measured as the
/// standard deviation of energy predictions across ensemble members.
/// Returns the index of the candidate with the highest disagreement.
///
/// Returns `None` if the pool is empty.
#[allow(dead_code)]
pub fn query_by_committee(ensemble: &EnsembleNnp, pool: &[Vec<f64>]) -> Option<usize> {
    if pool.is_empty() {
        return None;
    }
    let mut max_std = -1.0_f64;
    let mut best_idx = 0_usize;
    for (idx, desc) in pool.iter().enumerate() {
        let (_mean, std) = ensemble.predict_energy(desc);
        if std > max_std {
            max_std = std;
            best_idx = idx;
        }
    }
    Some(best_idx)
}
/// Returns the disagreement (std) for every candidate in the pool.
///
/// Useful for ranking all candidates by uncertainty.
#[allow(dead_code)]
pub fn committee_disagreements(ensemble: &EnsembleNnp, pool: &[Vec<f64>]) -> Vec<f64> {
    pool.iter()
        .map(|desc| ensemble.predict_energy(desc).1)
        .collect()
}
/// Run a short NVE trajectory and measure energy drift.
///
/// Returns `(initial_energy, final_energy, max_drift)` where `max_drift` is the
/// maximum |E(t) - E(0)| normalised by |E(0)|.
#[allow(dead_code)]
pub fn energy_conservation_check(
    potential: &FeedForwardPotential,
    mut system: NnAtomisticSystem,
    n_steps: usize,
    dt: f64,
) -> (f64, f64, f64) {
    let e0 = system.total_energy(potential);
    let mut max_drift = 0.0_f64;
    for _ in 0..n_steps {
        system.step(potential, dt);
        let e = system.total_energy(potential);
        let drift = if e0.abs() > 1e-30 {
            (e - e0).abs() / e0.abs()
        } else {
            (e - e0).abs()
        };
        if drift > max_drift {
            max_drift = drift;
        }
    }
    let e_final = system.total_energy(potential);
    (e0, e_final, max_drift)
}
#[cfg(test)]
mod tests {
    use super::super::types::*;
    use super::*;
    use crate::Activation;
    use crate::nn_potential::BehlerParrinello;
    fn simple_potential() -> FeedForwardPotential {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff)
    }
    #[test]
    fn test_dense_layer_forward() {
        let layer = DenseLayer::new(
            vec![vec![1.0, 2.0], vec![3.0, 4.0]],
            vec![0.5, -0.5],
            Activation::Linear,
        );
        let out = layer.forward(&[1.0, 1.0]);
        assert!((out[0] - 3.5).abs() < 1e-12);
        assert!((out[1] - 6.5).abs() < 1e-12);
    }
    #[test]
    fn test_dense_layer_relu() {
        let layer = DenseLayer::new(
            vec![vec![1.0], vec![-1.0]],
            vec![0.0, 0.0],
            Activation::ReLU,
        );
        let out = layer.forward(&[2.0]);
        assert!((out[0] - 2.0).abs() < 1e-12);
        assert!(out[1].abs() < 1e-12);
    }
    #[test]
    fn test_potential_finite_energy() {
        let pot = simple_potential();
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(1.0, 1.5, 0.0),
        ];
        let species = vec![1, 1, 1];
        let (energy, _forces) = pot.energy_and_forces(&positions, &species);
        assert!(energy.is_finite(), "Energy should be finite, got {energy}");
    }
    #[test]
    fn test_forces_finite_and_correct_dimension() {
        let pot = simple_potential();
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1, 1];
        let (_energy, forces) = pot.energy_and_forces(&positions, &species);
        assert_eq!(forces.len(), 2);
        for f in &forces {
            assert!(f[0].is_finite());
            assert!(f[1].is_finite());
            assert!(f[2].is_finite());
        }
    }
    #[test]
    fn test_newton_third_law_approx() {
        let pot = simple_potential();
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(1.5, 0.0, 0.0)];
        let species = vec![1, 1];
        let (_energy, forces) = pot.energy_and_forces(&positions, &species);
        let total = forces[0] + forces[1];
        let mag = total.norm();
        assert!(
            mag < 1e-3,
            "Total force should be ~0 (Newton III), got magnitude {mag}"
        );
    }
    #[test]
    fn test_normalizer_fit_zero_mean() {
        let data = vec![
            vec![1.0, 2.0, 3.0],
            vec![3.0, 4.0, 5.0],
            vec![5.0, 6.0, 7.0],
        ];
        let norm = DescriptorNormalizer::fit(&data);
        assert!((norm.mean[0] - 3.0).abs() < 1e-12);
        assert!((norm.mean[1] - 4.0).abs() < 1e-12);
        assert!((norm.mean[2] - 5.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalizer_normalize_and_denormalize() {
        let norm = DescriptorNormalizer::new(vec![2.0, 4.0], vec![1.0, 2.0]);
        let input = vec![4.0, 8.0];
        let normalized = norm.normalize_vec(&input);
        assert!((normalized[0] - 2.0).abs() < 1e-12);
        assert!((normalized[1] - 2.0).abs() < 1e-12);
        let restored = norm.denormalize(&normalized);
        assert!((restored[0] - 4.0).abs() < 1e-12);
        assert!((restored[1] - 8.0).abs() < 1e-12);
    }
    #[test]
    fn test_normalizer_in_place() {
        let norm = DescriptorNormalizer::new(vec![0.0], vec![2.0]);
        let mut v = vec![6.0];
        norm.normalize(&mut v);
        assert!((v[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_ensemble_mean_energy_finite() {
        let pot1 = simple_potential();
        let pot2 = simple_potential();
        let ensemble = EnsembleNNP::new(vec![pot1, pot2]);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1, 1];
        let mean_e = ensemble.mean_energy(&positions, &species);
        assert!(
            mean_e.is_finite(),
            "mean energy should be finite, got {mean_e}"
        );
    }
    #[test]
    fn test_ensemble_uncertainty_zero_for_identical() {
        let pot1 = simple_potential();
        let pot2 = pot1.clone();
        let ensemble = EnsembleNNP::new(vec![pot1, pot2]);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1, 1];
        let unc = ensemble.uncertainty(&positions, &species);
        assert!(
            unc < 1e-10,
            "uncertainty should be ~0 for identical models, got {unc}"
        );
    }
    #[test]
    fn test_ensemble_n_models() {
        let pot = simple_potential();
        let ensemble = EnsembleNNP::new(vec![pot.clone(), pot.clone(), pot]);
        assert_eq!(ensemble.n_models(), 3);
    }
    #[test]
    fn test_transfer_learning_freeze() {
        let pot = simple_potential();
        let mut tl = TransferLearning::new(pot, vec![1]);
        assert!(tl.is_trainable(0), "unfrozen: all layers trainable");
        tl.freeze_base();
        assert!(
            !tl.is_trainable(0),
            "frozen: layer 0 should not be trainable"
        );
        assert!(tl.is_trainable(1), "frozen: layer 1 should be trainable");
    }
    #[test]
    fn test_transfer_learning_forward() {
        let pot = simple_potential();
        let tl = TransferLearning::new(pot, vec![]);
        let descriptor = vec![0.1, 0.2, 0.3, 0.4];
        let e = tl.forward(&descriptor);
        assert!(e.is_finite(), "forward output should be finite, got {e}");
    }
    #[test]
    fn test_transfer_learning_n_trainable() {
        let pot = simple_potential();
        let n_total = pot.layers.len();
        let mut tl = TransferLearning::new(pot, vec![1]);
        assert_eq!(tl.n_trainable(), n_total);
        tl.freeze_base();
        assert_eq!(tl.n_trainable(), 1);
    }
    #[test]
    fn test_transfer_learning_unfreeze() {
        let pot = simple_potential();
        let n_total = pot.layers.len();
        let mut tl = TransferLearning::new(pot, vec![0]);
        tl.freeze_base();
        assert_eq!(tl.n_trainable(), 1);
        tl.unfreeze();
        assert_eq!(tl.n_trainable(), n_total);
    }
    #[test]
    fn test_g5_zero_beyond_cutoff() {
        let v = compute_g5(6.0, 1.0, 0.5, 0.1, 1.0, 1.0, 5.0);
        assert_eq!(v, 0.0, "G5 should be 0 when r_ij >= rc");
    }
    #[test]
    fn test_g5_positive_inside_cutoff() {
        let v = compute_g5(1.0, 1.0, 0.5, 0.1, 1.0, 1.0, 5.0);
        assert!(v > 0.0, "G5 should be positive inside cutoff, got {v}");
    }
    #[test]
    fn test_g5_angular_sensitivity() {
        let v_pos = compute_g5(1.0, 1.0, 1.0, 0.1, 1.0, 1.0, 5.0);
        let v_neg = compute_g5(1.0, 1.0, -1.0, 0.1, 1.0, 1.0, 5.0);
        assert!(v_pos > v_neg, "G5 should be larger for cos=1 than cos=-1");
    }
    #[test]
    fn test_g3_zero_beyond_cutoff() {
        let v = compute_g3(6.0, 1.0, 5.0);
        assert_eq!(v, 0.0, "G3 should be 0 at/beyond cutoff");
    }
    #[test]
    fn test_g3_finite_inside() {
        let v = compute_g3(1.0, 1.0, 5.0);
        assert!(v.is_finite(), "G3 should be finite inside cutoff, got {v}");
    }
    #[test]
    fn test_nn_pair_potential_zero_beyond_cutoff() {
        let pot = NnPairPotential::new(4, 3.0);
        let e = pot.evaluate(5.0);
        assert_eq!(e, 0.0, "pair energy should be 0 beyond cutoff");
    }
    #[test]
    fn test_nn_pair_potential_finite_inside() {
        let pot = NnPairPotential::new(4, 3.0);
        let e = pot.evaluate(1.5);
        assert!(e.is_finite(), "pair energy should be finite, got {e}");
    }
    #[test]
    fn test_nn_pair_energy_and_force() {
        let pot = NnPairPotential::new(4, 3.0);
        let (e, f) = pot.energy_and_force(1.5);
        assert!(e.is_finite(), "energy = {e}");
        assert!(f.is_finite(), "force = {f}");
    }
    #[test]
    fn test_nn_pair_total_energy_two_atoms() {
        let pot = NnPairPotential::new(4, 5.0);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let e = pot.total_energy(&positions);
        assert!(e.is_finite(), "total energy = {e}");
    }
    #[test]
    fn test_energy_force_consistency_small_system() {
        let pot = simple_potential();
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1, 1];
        let max_diff = energy_force_consistency(&pot, &positions, &species, 1e-5);
        assert!(max_diff < 1e-3, "Force consistency error = {max_diff}");
    }
    #[test]
    fn test_descriptor_stats_mean() {
        let data = vec![vec![0.0f64, 2.0], vec![2.0, 4.0]];
        let stats = DescriptorStats::compute(&data);
        assert!(
            (stats.mean[0] - 1.0).abs() < 1e-12,
            "mean[0] = {}",
            stats.mean[0]
        );
        assert!(
            (stats.mean[1] - 3.0).abs() < 1e-12,
            "mean[1] = {}",
            stats.mean[1]
        );
    }
    #[test]
    fn test_descriptor_stats_min_max() {
        let data = vec![vec![1.0f64], vec![3.0], vec![2.0]];
        let stats = DescriptorStats::compute(&data);
        assert!((stats.min[0] - 1.0).abs() < 1e-12);
        assert!((stats.max[0] - 3.0).abs() < 1e-12);
    }
    #[test]
    fn test_descriptor_stats_range() {
        let data = vec![vec![0.0f64, 5.0], vec![4.0, 15.0]];
        let stats = DescriptorStats::compute(&data);
        let range = stats.range();
        assert!((range[0] - 4.0).abs() < 1e-12);
        assert!((range[1] - 10.0).abs() < 1e-12);
    }
    #[test]
    fn test_descriptor_stats_std_zero_constant() {
        let data = vec![vec![5.0f64]; 3];
        let stats = DescriptorStats::compute(&data);
        assert!(
            stats.std[0].abs() < 1e-10,
            "std of constant = {}",
            stats.std[0]
        );
    }
    #[test]
    fn test_nnp_normalizer_fit_mean() {
        let data = vec![vec![1.0, 10.0], vec![3.0, 20.0], vec![5.0, 30.0]];
        let norm = NnpNormalizer::fit(&data);
        assert!(
            (norm.mean[0] - 3.0).abs() < 1e-12,
            "mean[0] = {}",
            norm.mean[0]
        );
        assert!(
            (norm.mean[1] - 20.0).abs() < 1e-12,
            "mean[1] = {}",
            norm.mean[1]
        );
    }
    #[test]
    fn test_nnp_normalizer_transform_zero_mean() {
        let data = vec![vec![0.0], vec![2.0], vec![4.0]];
        let norm = NnpNormalizer::fit(&data);
        let transformed = norm.transform(&[2.0]);
        assert!(
            transformed[0].abs() < 1e-12,
            "transformed[0] = {}",
            transformed[0]
        );
    }
    #[test]
    fn test_nnp_normalizer_roundtrip() {
        let data = vec![vec![1.0, 5.0], vec![3.0, 9.0], vec![5.0, 13.0]];
        let norm = NnpNormalizer::fit(&data);
        let original = vec![2.5, 7.0];
        let transformed = norm.transform(&original);
        let restored = norm.inverse_transform(&transformed);
        for (i, (&o, &r)) in original.iter().zip(restored.iter()).enumerate() {
            assert!(
                (o - r).abs() < 1e-12,
                "roundtrip failed at index {i}: {o} != {r}"
            );
        }
    }
    #[test]
    fn test_nnp_normalizer_constant_feature_no_panic() {
        let data = vec![vec![3.0]; 5];
        let norm = NnpNormalizer::fit(&data);
        assert!(norm.std[0] >= 1e-12, "std floor must be >= 1e-12");
        let _ = norm.transform(&[3.0]);
    }
    #[test]
    fn test_ensemble_nnp_predict_energy_finite() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot1 = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym.clone(), cutoff);
        let pot2 = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot1, pot2]);
        let desc = vec![0.1, 0.2, 0.3, 0.4];
        let (mean, std) = ens.predict_energy(&desc);
        assert!(mean.is_finite(), "mean energy should be finite, got {mean}");
        assert!(
            std.is_finite() && std >= 0.0,
            "std should be non-negative finite, got {std}"
        );
    }
    #[test]
    fn test_ensemble_nnp_identical_models_zero_std() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot.clone(), pot]);
        let desc = vec![0.1, 0.2, 0.3, 0.4];
        let (_mean, std) = ens.predict_energy(&desc);
        assert!(std < 1e-12, "identical models should have std=0, got {std}");
    }
    #[test]
    fn test_ensemble_nnp_predict_forces_shape() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot]);
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let forces = ens.predict_forces(&positions);
        assert_eq!(forces.len(), 2, "force count should equal atom count");
        for f in &forces {
            assert!(
                f.iter().all(|&x| x.is_finite()),
                "all force components should be finite"
            );
        }
    }
    #[test]
    fn test_fine_tune_reduces_loss() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let mut tl = TransferLearning::new(pot, vec![1]);
        let desc = vec![0.5, 0.3, 0.2, 0.1];
        let target = 1.0;
        let e_before = tl.forward(&desc);
        tl.fine_tune(&[(desc.clone(), target)], 0.01, 5);
        let e_after = tl.forward(&desc);
        let loss_before = (e_before - target).abs();
        let loss_after = (e_after - target).abs();
        assert!(
            loss_after <= loss_before + 1e-6,
            "fine_tune should not increase loss: before={loss_before}, after={loss_after}"
        );
    }
    #[test]
    fn test_fine_tune_frozen_only_updates_fine_tune_layers() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let layer0_w_before = pot.layers[0].weights.clone();
        let mut tl = TransferLearning::new(pot, vec![1]);
        tl.freeze_base();
        let desc = vec![0.5, 0.3, 0.2, 0.1];
        tl.fine_tune(&[(desc, 1.0)], 0.01, 3);
        let layer0_w_after = &tl.base.layers[0].weights;
        assert_eq!(
            &layer0_w_before, layer0_w_after,
            "frozen layer 0 weights should not change during fine_tune"
        );
    }
    #[test]
    fn test_chebyshev_descriptor_zero_at_cutoff() {
        let d = chebyshev_descriptor(5.0, 5.0, 4);
        for (k, &v) in d.iter().enumerate() {
            assert!(v.abs() < 1e-12, "chebyshev[{k}] should be 0 at rc, got {v}");
        }
    }
    #[test]
    fn test_chebyshev_descriptor_first_term_at_origin() {
        let d = chebyshev_descriptor(0.0, 5.0, 3);
        assert!(
            (d[0] - 1.0).abs() < 1e-10,
            "T_0 * fc at r=0 should be 1.0, got {}",
            d[0]
        );
    }
    #[test]
    fn test_chebyshev_descriptor_length() {
        let n = 6;
        let d = chebyshev_descriptor(1.0, 5.0, n);
        assert_eq!(d.len(), n, "descriptor length should be n_basis");
    }
    #[test]
    fn test_pair_descriptor_length() {
        let d = pair_descriptor(1.5, 5.0, 4);
        assert_eq!(d.len(), 6, "pair descriptor: 2 raw + 4 Chebyshev");
    }
    #[test]
    fn test_pair_descriptor_first_element_is_r() {
        let r = 2.3_f64;
        let d = pair_descriptor(r, 5.0, 3);
        assert!((d[0] - r).abs() < 1e-14, "first element should be r");
    }
    #[test]
    fn test_g4_zero_beyond_cutoff() {
        let v = compute_g4(6.0, 1.0, 1.0, 0.5, 0.1, 1.0, 1.0, 5.0);
        assert_eq!(v, 0.0, "G4 should be 0 when any distance >= rc");
    }
    #[test]
    fn test_g4_positive_inside_cutoff() {
        let v = compute_g4(1.0, 1.0, 1.4, 0.5, 0.05, 1.0, 1.0, 5.0);
        assert!(v > 0.0, "G4 should be positive inside cutoff, got {v}");
    }
    #[test]
    fn test_g4_angular_discrimination() {
        let v_col = compute_g4(1.0, 1.0, 1.4, 1.0, 0.05, 2.0, 1.0, 5.0);
        let v_anti = compute_g4(1.0, 1.0, 1.4, -1.0, 0.05, 2.0, 1.0, 5.0);
        assert!(v_col > v_anti, "G4 should differentiate angle (lambda=+1)");
    }
    #[test]
    fn test_mpnn_total_energy_finite() {
        let mpnn = Mpnn::new(4, 2, 5.0);
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 1.5, 0.0]];
        let e = mpnn.total_energy(&positions);
        assert!(e.is_finite(), "MPNN total energy should be finite, got {e}");
    }
    #[test]
    fn test_mpnn_atom_energies_length() {
        let mpnn = Mpnn::new(4, 1, 5.0);
        let positions = vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let ae = mpnn.atom_energies(&positions);
        assert_eq!(ae.len(), 2, "atom_energies length should equal atom count");
    }
    #[test]
    fn test_mpnn_layer_forward() {
        let layer = MpnnLayer::random(2, Activation::Tanh, 42);
        let h_self = vec![1.0, 0.5];
        let neigh_sum = vec![0.3, 0.1];
        let out = layer.forward_node(&h_self, &neigh_sum);
        assert_eq!(out.len(), 2);
        for &v in &out {
            assert!(
                v.is_finite() && v.abs() <= 1.0 + 1e-10,
                "tanh output must be in [-1,1]"
            );
        }
    }
    #[test]
    fn test_query_by_committee_selects_index() {
        let cutoff = 5.0;
        let sym1 = FeedForwardPotential::default_symmetry_set(cutoff);
        let sym2 = FeedForwardPotential::default_symmetry_set(cutoff);
        let p1 = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym1, cutoff);
        let p2 = FeedForwardPotential::new(&[4, 8, 1], Activation::Sigmoid, sym2, cutoff);
        let ens = EnsembleNnp::new(vec![p1, p2]);
        let pool = vec![
            vec![0.1, 0.2, 0.3, 0.4],
            vec![1.0, 2.0, 3.0, 4.0],
            vec![0.0, 0.0, 0.0, 0.0],
        ];
        let idx = query_by_committee(&ens, &pool);
        assert!(idx.is_some(), "should return a valid index");
        assert!(
            idx.unwrap() < pool.len(),
            "index must be within pool bounds"
        );
    }
    #[test]
    fn test_query_by_committee_empty_pool() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot]);
        let idx = query_by_committee(&ens, &[]);
        assert!(idx.is_none(), "empty pool should return None");
    }
    #[test]
    fn test_committee_disagreements_length() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let p1 = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym.clone(), cutoff);
        let p2 = FeedForwardPotential::new(&[4, 4, 1], Activation::Sigmoid, sym, cutoff);
        let ens = EnsembleNnp::new(vec![p1, p2]);
        let pool = vec![vec![0.1, 0.2, 0.3, 0.4]; 5];
        let d = committee_disagreements(&ens, &pool);
        assert_eq!(d.len(), 5);
        for &v in &d {
            assert!(v >= 0.0, "disagreement must be non-negative");
        }
    }
    #[test]
    fn test_atomistic_system_kinetic_energy_at_rest() {
        let pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let sys = NnAtomisticSystem::new(pos, vec![1.0, 1.0], vec![1, 1]);
        assert!(sys.kinetic_energy().abs() < 1e-14, "KE at rest should be 0");
    }
    #[test]
    fn test_atomistic_system_step_finite() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let mut sys = NnAtomisticSystem::new(pos, vec![12.0, 12.0], vec![1, 1]);
        sys.step(&pot, 0.001);
        for p in &sys.positions {
            assert!(p[0].is_finite() && p[1].is_finite() && p[2].is_finite());
        }
    }
    #[test]
    fn test_energy_conservation_short_run() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let mut sys = NnAtomisticSystem::new(pos, vec![12.0, 12.0], vec![1, 1]);
        sys.velocities[0] = Vec3::new(0.01, 0.0, 0.0);
        sys.velocities[1] = Vec3::new(-0.01, 0.0, 0.0);
        let (_e0, _ef, max_drift) = energy_conservation_check(&pot, sys, 5, 0.0005);
        assert!(
            max_drift.is_finite(),
            "energy drift must be finite: {max_drift}"
        );
    }
    #[test]
    fn test_nnp_normalizer_fit_mean_correct() {
        let data = vec![vec![1.0, 4.0], vec![3.0, 2.0]];
        let norm = NnpNormalizer::fit(&data);
        assert!(
            (norm.mean[0] - 2.0).abs() < 1e-12,
            "mean[0] = {}",
            norm.mean[0]
        );
        assert!(
            (norm.mean[1] - 3.0).abs() < 1e-12,
            "mean[1] = {}",
            norm.mean[1]
        );
    }
    #[test]
    fn test_nnp_normalizer_fit_std_correct() {
        let data = vec![vec![0.0], vec![2.0]];
        let norm = NnpNormalizer::fit(&data);
        assert!(
            (norm.std[0] - 1.0).abs() < 1e-10,
            "std[0] = {}",
            norm.std[0]
        );
    }
    #[test]
    fn test_nnp_normalizer_fit_constant_std_clamped() {
        let data = vec![vec![5.0], vec![5.0], vec![5.0]];
        let norm = NnpNormalizer::fit(&data);
        assert!(norm.std[0] >= 1e-12, "std should be clamped to >= 1e-12");
    }
    #[test]
    fn test_nnp_normalizer_transform_zero_mean_v2() {
        let data = vec![vec![1.0, 2.0], vec![3.0, 4.0]];
        let norm = NnpNormalizer::fit(&data);
        let x = norm.transform(&[2.0, 3.0]);
        assert!(
            x[0].abs() < 1e-12,
            "transform(mean) should be 0, got {}",
            x[0]
        );
        assert!(
            x[1].abs() < 1e-12,
            "transform(mean) should be 0, got {}",
            x[1]
        );
    }
    #[test]
    fn test_nnp_normalizer_inverse_transform_roundtrip() {
        let data = vec![vec![0.0, 5.0], vec![4.0, 1.0], vec![2.0, 3.0]];
        let norm = NnpNormalizer::fit(&data);
        let original = vec![1.5, 2.7];
        let transformed = norm.transform(&original);
        let recovered = norm.inverse_transform(&transformed);
        for i in 0..2 {
            assert!(
                (recovered[i] - original[i]).abs() < 1e-10,
                "inverse roundtrip failed at index {i}: {} vs {}",
                recovered[i],
                original[i]
            );
        }
    }
    #[test]
    fn test_nnp_normalizer_transform_length_preserved() {
        let data = vec![vec![1.0, 2.0, 3.0]];
        let norm = NnpNormalizer::fit(&data);
        let x = norm.transform(&[1.0, 2.0, 3.0]);
        assert_eq!(x.len(), 3);
    }
    #[test]
    fn test_nnp_normalizer_inverse_transform_length_preserved() {
        let data = vec![vec![1.0, 2.0]];
        let norm = NnpNormalizer::fit(&data);
        let x = norm.inverse_transform(&[0.0, 1.0]);
        assert_eq!(x.len(), 2);
    }
    #[test]
    fn test_ensemble_nnp_predict_energy_mean_finite() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot]);
        let desc = vec![0.1, 0.2, 0.3, 0.4];
        let (mean, std) = ens.predict_energy(&desc);
        assert!(mean.is_finite(), "mean energy should be finite, got {mean}");
        assert!(std >= 0.0, "std should be >= 0, got {std}");
    }
    #[test]
    fn test_ensemble_nnp_std_zero_for_identical_models() {
        let cutoff = 5.0;
        let sym1 = FeedForwardPotential::default_symmetry_set(cutoff);
        let sym2 = FeedForwardPotential::default_symmetry_set(cutoff);
        let p1 = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym1, cutoff);
        let p2 = p1.clone();
        let p3 = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym2, cutoff);
        let _ = p3;
        let ens = EnsembleNnp::new(vec![p1.clone(), p1.clone()]);
        let desc = vec![0.5, 0.5, 0.5, 0.5];
        let (_mean, std) = ens.predict_energy(&desc);
        assert!(
            std < 1e-10,
            "identical models should give std ≈ 0, got {std}"
        );
        let _ = p2;
    }
    #[test]
    fn test_ensemble_nnp_n_models() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let p = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![p.clone(), p.clone(), p]);
        assert_eq!(ens.n_models(), 3);
    }
    #[test]
    fn test_ensemble_nnp_predict_forces_length() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot]);
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let forces = ens.predict_forces(&positions);
        assert_eq!(forces.len(), 2, "forces length should equal n_atoms");
    }
    #[test]
    fn test_ensemble_nnp_predict_forces_finite() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![pot]);
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let forces = ens.predict_forces(&positions);
        for f in &forces {
            for &c in f {
                assert!(c.is_finite(), "force component must be finite");
            }
        }
    }
    #[test]
    fn test_g4_beyond_cutoff_returns_zero() {
        let g = compute_g4(12.0, 2.0, 2.0, 0.5, 0.1, 1.0, 1.0, 10.0);
        assert_eq!(g, 0.0, "G4 with r_ij > rc should be 0");
    }
    #[test]
    fn test_g4_within_cutoff_finite() {
        let g = compute_g4(2.0, 2.0, 2.0, 0.5, 0.1, 1.0, 1.0, 10.0);
        assert!(g.is_finite(), "G4 within cutoff should be finite, got {g}");
    }
    #[test]
    fn test_g4_positive_for_lambda_positive() {
        let g = compute_g4(2.0, 2.0, 2.0, 0.5, 0.01, 1.0, 1.0, 10.0);
        assert!(g >= 0.0, "G4 with lambda=1, cos>0 should be >= 0, got {g}");
    }
    #[test]
    fn test_g4_symmetric_in_rij_rik() {
        let g1 = compute_g4(2.0, 3.0, 2.5, 0.5, 0.1, 1.0, 1.0, 10.0);
        let g2 = compute_g4(3.0, 2.0, 2.5, 0.5, 0.1, 1.0, 1.0, 10.0);
        assert!(
            (g1 - g2).abs() < 1e-12,
            "G4 should be symmetric in r_ij/r_ik: {g1} vs {g2}"
        );
    }
    #[test]
    fn test_g5_within_cutoff_finite() {
        let g = compute_g5(2.0, 2.0, 0.5, 0.1, 1.0, 1.0, 10.0);
        assert!(g.is_finite(), "G5 within cutoff should be finite, got {g}");
    }
    #[test]
    fn test_g5_beyond_cutoff_returns_zero() {
        let g = compute_g5(12.0, 2.0, 0.5, 0.1, 1.0, 1.0, 10.0);
        assert_eq!(g, 0.0, "G5 with r_ij > rc should be 0");
    }
    #[test]
    fn test_g3_within_cutoff_finite() {
        let g = compute_g3(3.0, 0.5, 10.0);
        assert!(g.is_finite(), "G3 within cutoff should be finite, got {g}");
    }
    #[test]
    fn test_g3_beyond_cutoff_zero() {
        let g = compute_g3(11.0, 0.5, 10.0);
        assert_eq!(g, 0.0, "G3 beyond cutoff should be 0");
    }
    #[test]
    fn test_g3_at_zero_r_is_one_times_fc() {
        let g = compute_g3(0.0, 0.5, 10.0);
        assert!(
            (g - 1.0).abs() < 1e-12,
            "G3(0, kappa, rc) should be 1.0, got {g}"
        );
    }
    #[test]
    fn test_chebyshev_descriptor_length_v2() {
        let rc = 5.0;
        let n_basis = 6;
        let d = chebyshev_descriptor(3.0, rc, n_basis);
        assert_eq!(d.len(), n_basis, "descriptor length should match n_basis");
    }
    #[test]
    fn test_chebyshev_descriptor_at_rc_is_zero() {
        let d = chebyshev_descriptor(5.0, 5.0, 4);
        for (i, &v) in d.iter().enumerate() {
            assert!(
                v.abs() < 1e-12,
                "Chebyshev component {i} at r=rc should be 0, got {v}"
            );
        }
    }
    #[test]
    fn test_chebyshev_descriptor_beyond_cutoff_zero() {
        let d = chebyshev_descriptor(6.0, 5.0, 4);
        for (i, &v) in d.iter().enumerate() {
            assert!(
                v.abs() < 1e-12,
                "Chebyshev component {i} beyond cutoff should be 0, got {v}"
            );
        }
    }
    #[test]
    fn test_chebyshev_descriptor_within_cutoff_nonzero() {
        let d = chebyshev_descriptor(2.5, 5.0, 4);
        let any_nonzero = d.iter().any(|&v| v.abs() > 1e-20);
        assert!(
            any_nonzero,
            "at least one component should be nonzero within cutoff"
        );
    }
    #[test]
    fn test_ensemble_nnp_capital_mean_energy_finite() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNNP::new(vec![pot]);
        let pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1u32, 1];
        let e = ens.mean_energy(&pos, &species);
        assert!(
            e.is_finite(),
            "EnsembleNNP mean_energy should be finite, got {e}"
        );
    }
    #[test]
    fn test_ensemble_nnp_capital_uncertainty_nonnegative() {
        let cutoff = 5.0;
        let sym1 = FeedForwardPotential::default_symmetry_set(cutoff);
        let sym2 = FeedForwardPotential::default_symmetry_set(cutoff);
        let p1 = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym1, cutoff);
        let p2 = FeedForwardPotential::new(&[4, 4, 1], Activation::Sigmoid, sym2, cutoff);
        let ens = EnsembleNNP::new(vec![p1, p2]);
        let pos = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1u32, 1];
        let u = ens.uncertainty(&pos, &species);
        assert!(u >= 0.0, "uncertainty must be non-negative, got {u}");
    }
    #[test]
    fn test_ensemble_nnp_capital_n_models() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let p = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNNP::new(vec![p.clone(), p]);
        assert_eq!(ens.n_models(), 2);
    }
    #[test]
    fn test_descriptor_stats_mean_correct() {
        let data = vec![vec![2.0, 4.0], vec![4.0, 6.0]];
        let stats = DescriptorStats::compute(&data);
        assert!(
            (stats.mean[0] - 3.0).abs() < 1e-12,
            "mean[0] = {}",
            stats.mean[0]
        );
        assert!(
            (stats.mean[1] - 5.0).abs() < 1e-12,
            "mean[1] = {}",
            stats.mean[1]
        );
    }
    #[test]
    fn test_descriptor_stats_std_nonnegative() {
        let data = vec![vec![1.0, 5.0], vec![3.0, 7.0], vec![5.0, 9.0]];
        let stats = DescriptorStats::compute(&data);
        for (i, &s) in stats.std.iter().enumerate() {
            assert!(s >= 0.0, "std[{i}] must be non-negative, got {s}");
        }
    }
    #[test]
    fn test_descriptor_stats_min_max_bounds() {
        let data = vec![vec![1.0, 10.0], vec![3.0, 6.0], vec![5.0, 8.0]];
        let stats = DescriptorStats::compute(&data);
        assert!((stats.min[0] - 1.0).abs() < 1e-12);
        assert!((stats.max[0] - 5.0).abs() < 1e-12);
        assert!((stats.min[1] - 6.0).abs() < 1e-12);
        assert!((stats.max[1] - 10.0).abs() < 1e-12);
    }
    /// Gradient output must have correct length: n_desc × N_atoms × 3.
    #[test]
    fn test_descriptor_gradient_output_length() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let positions = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(0.0, 2.0, 0.0),
        ];
        let species = vec![1u32; 3];
        let n_desc = pot.symmetry_functions.functions.len();
        let grad = pot.compute_descriptor_gradient(&positions, 0, &species, 1e-5);
        let expected_len = n_desc * positions.len() * 3;
        assert_eq!(
            grad.len(),
            expected_len,
            "gradient length: got {}, expected {expected_len}",
            grad.len()
        );
    }
    /// Gradient must be all finite values.
    #[test]
    fn test_descriptor_gradient_finite() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let species = vec![1u32; 2];
        let grad = pot.compute_descriptor_gradient(&positions, 0, &species, 1e-5);
        for (i, &g) in grad.iter().enumerate() {
            assert!(g.is_finite(), "gradient[{i}] is not finite: {g}");
        }
    }
    /// If all atoms are far beyond the cutoff, all gradient components must be zero.
    #[test]
    fn test_descriptor_gradient_beyond_cutoff_zero() {
        let cutoff = 3.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let pot = FeedForwardPotential::new(&[4, 8, 1], Activation::Tanh, sym, cutoff);
        let positions = vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(100.0, 0.0, 0.0)];
        let species = vec![1u32; 2];
        let grad = pot.compute_descriptor_gradient(&positions, 0, &species, 1e-5);
        for (i, &g) in grad.iter().enumerate() {
            assert!(
                g.abs() < 1e-10,
                "gradient[{i}] must be ~0 when neighbours are beyond cutoff: {g}"
            );
        }
    }
    /// Disagreement output length must equal descriptor dimension.
    #[test]
    fn test_ensemble_nnp_disagreement_length() {
        let cutoff = 5.0;
        let sym1 = FeedForwardPotential::default_symmetry_set(cutoff);
        let sym2 = FeedForwardPotential::default_symmetry_set(cutoff);
        let p1 = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym1, cutoff);
        let p2 = FeedForwardPotential::new(&[4, 4, 1], Activation::Sigmoid, sym2, cutoff);
        let ens = EnsembleNnp::new(vec![p1, p2]);
        let desc = vec![0.1, 0.5, 1.0, 2.0];
        let disag = ens.compute_disagreement(&desc, 1e-3);
        assert_eq!(disag.len(), desc.len(), "disagreement length mismatch");
    }
    /// Disagreement values must be non-negative.
    #[test]
    fn test_ensemble_nnp_disagreement_nonnegative() {
        let cutoff = 5.0;
        let sym1 = FeedForwardPotential::default_symmetry_set(cutoff);
        let sym2 = FeedForwardPotential::default_symmetry_set(cutoff);
        let p1 = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym1, cutoff);
        let p2 = FeedForwardPotential::new(&[4, 4, 1], Activation::ReLU, sym2, cutoff);
        let ens = EnsembleNnp::new(vec![p1, p2]);
        let desc = vec![0.2, 0.3, 0.4, 0.5];
        let disag = ens.compute_disagreement(&desc, 1e-3);
        for (k, &d) in disag.iter().enumerate() {
            assert!(d >= 0.0, "disagreement[{k}] must be non-negative, got {d}");
        }
    }
    /// A single-model ensemble must have zero disagreement.
    #[test]
    fn test_ensemble_nnp_single_model_zero_disagreement() {
        let cutoff = 5.0;
        let sym = FeedForwardPotential::default_symmetry_set(cutoff);
        let p = FeedForwardPotential::new(&[4, 4, 1], Activation::Tanh, sym, cutoff);
        let ens = EnsembleNnp::new(vec![p]);
        let desc = vec![0.1, 0.2, 0.3, 0.4];
        let disag = ens.compute_disagreement(&desc, 1e-3);
        for (k, &d) in disag.iter().enumerate() {
            assert!(
                d.abs() < 1e-12,
                "single model: disagreement[{k}] must be 0, got {d}"
            );
        }
    }
    /// Batch result must have N_atoms rows, each of length N_desc.
    #[test]
    fn test_bp_batch_output_shape() {
        let cutoff = 5.0;
        let bp = BehlerParrinello::default_radial(cutoff);
        let n_desc = bp.n_descriptors();
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let species = [1u32, 1, 1];
        let batch = bp.compute_symmetry_functions_batch(&positions, &species);
        assert_eq!(batch.len(), 3, "batch must have one row per atom");
        for (i, row) in batch.iter().enumerate() {
            assert_eq!(row.len(), n_desc, "batch[{i}] length mismatch");
        }
    }
    /// All descriptor values must be finite.
    #[test]
    fn test_bp_batch_finite() {
        let cutoff = 5.0;
        let bp = BehlerParrinello::default_radial(cutoff);
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0]];
        let species = [1u32, 1];
        let batch = bp.compute_symmetry_functions_batch(&positions, &species);
        for (i, row) in batch.iter().enumerate() {
            for (k, &v) in row.iter().enumerate() {
                assert!(v.is_finite(), "batch[{i}][{k}] not finite: {v}");
            }
        }
    }
    /// Isolated atom (all others beyond cutoff) must have descriptor = 0.
    #[test]
    fn test_bp_batch_isolated_atom_zero_descriptor() {
        let cutoff = 3.0;
        let bp = BehlerParrinello::default_radial(cutoff);
        let positions = [[0.0, 0.0, 0.0_f64], [100.0, 0.0, 0.0]];
        let species = [1u32, 1];
        let batch = bp.compute_symmetry_functions_batch(&positions, &species);
        for (k, &v) in batch[0].iter().enumerate() {
            assert!(
                v.abs() < 1e-15,
                "isolated atom descriptor[{k}] should be 0, got {v}"
            );
        }
    }
}
/// Collect all angular triplets (i, j, k) within cutoff.
#[allow(dead_code)]
pub fn collect_triplets(positions: &[[f64; 3]], cutoff: f64) -> Vec<AngularTriplet> {
    let n = positions.len();
    let mut triplets = Vec::new();
    for i in 0..n {
        let mut neighbours: Vec<usize> = Vec::new();
        for j in 0..n {
            if i == j {
                continue;
            }
            let dx = positions[j][0] - positions[i][0];
            let dy = positions[j][1] - positions[i][1];
            let dz = positions[j][2] - positions[i][2];
            if dx * dx + dy * dy + dz * dz < cutoff * cutoff {
                neighbours.push(j);
            }
        }
        for a in 0..neighbours.len() {
            for b in (a + 1)..neighbours.len() {
                triplets.push(AngularTriplet {
                    i,
                    j: neighbours[a],
                    k: neighbours[b],
                });
            }
        }
    }
    triplets
}
/// Compute cosine of the angle at atom i formed by atoms j and k.
#[allow(dead_code)]
pub fn triplet_cos_angle(positions: &[[f64; 3]], triplet: AngularTriplet) -> f64 {
    let ri = positions[triplet.i];
    let rj = positions[triplet.j];
    let rk = positions[triplet.k];
    let vj = [rj[0] - ri[0], rj[1] - ri[1], rj[2] - ri[2]];
    let vk = [rk[0] - ri[0], rk[1] - ri[1], rk[2] - ri[2]];
    let nj = (vj[0] * vj[0] + vj[1] * vj[1] + vj[2] * vj[2]).sqrt();
    let nk = (vk[0] * vk[0] + vk[1] * vk[1] + vk[2] * vk[2]).sqrt();
    if nj < 1e-15 || nk < 1e-15 {
        return 0.0;
    }
    let dot = vj[0] * vk[0] + vj[1] * vk[1] + vj[2] * vk[2];
    (dot / (nj * nk)).clamp(-1.0, 1.0)
}
/// Fourier-type angular basis: \[cos(m*theta)\] for m = 0..n_max.
#[allow(dead_code)]
pub fn angular_fourier_basis(cos_theta: f64, n_max: usize) -> Vec<f64> {
    let theta = cos_theta.clamp(-1.0, 1.0).acos();
    (0..n_max).map(|m| (m as f64 * theta).cos()).collect()
}
#[cfg(test)]
mod tests_deepmd_schnet {
    use super::super::types::*;
    use super::*;
    #[test]
    fn test_deepmd_descriptor_dim() {
        let d = DeepPotDescriptor::new(10, 6.0, 3.0);
        assert_eq!(d.descriptor_dim(), 40);
    }
    #[test]
    fn test_deepmd_descriptor_length() {
        let d = DeepPotDescriptor::new(5, 6.0, 3.0);
        let positions = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let desc = d.environment_matrix(&positions, 0);
        assert_eq!(desc.len(), d.descriptor_dim());
    }
    #[test]
    fn test_deepmd_descriptor_all_finite() {
        let d = DeepPotDescriptor::new(8, 6.0, 3.0);
        let positions = [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]];
        let desc = d.environment_matrix(&positions, 0);
        for (i, &v) in desc.iter().enumerate() {
            assert!(v.is_finite(), "desc[{i}] = {v}");
        }
    }
    #[test]
    fn test_deepmd_smooth_fn_zero_at_cutoff() {
        let d = DeepPotDescriptor::new(4, 5.0, 2.0);
        let s = d.smooth_fn(5.0);
        assert_eq!(s, 0.0, "smooth_fn at cutoff = 0");
    }
    #[test]
    fn test_deepmd_smooth_fn_positive_inside() {
        let d = DeepPotDescriptor::new(4, 5.0, 2.0);
        let s = d.smooth_fn(1.0);
        assert!(s > 0.0, "smooth_fn inside cutoff > 0, got {s}");
    }
    #[test]
    fn test_deepmd_potential_energy_finite() {
        let pot = DeepPotPotential::new(4, &[8, 8], &[16, 1], 5.0, 2.0);
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0]];
        let e = pot.total_energy(&positions);
        assert!(e.is_finite(), "energy = {e}");
    }
    #[test]
    fn test_deepmd_forces_length() {
        let pot = DeepPotPotential::new(4, &[8], &[8, 1], 5.0, 2.0);
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let forces = pot.forces(&positions);
        assert_eq!(forces.len(), 3);
    }
    #[test]
    fn test_deepmd_forces_all_finite() {
        let pot = DeepPotPotential::new(4, &[8], &[8, 1], 5.0, 2.0);
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0]];
        for f in pot.forces(&positions) {
            for &v in &f {
                assert!(v.is_finite(), "force = {v}");
            }
        }
    }
    #[test]
    fn test_gaussian_basis_length() {
        let basis = GaussianBasis::new(0.5, 5.5, 8, 0.5, 6.0);
        let vals = basis.evaluate(2.0);
        assert_eq!(vals.len(), 8);
    }
    #[test]
    fn test_gaussian_basis_zero_beyond_cutoff() {
        let basis = GaussianBasis::new(0.5, 5.5, 8, 0.5, 6.0);
        let vals = basis.evaluate(7.0);
        for &v in &vals {
            assert_eq!(v, 0.0);
        }
    }
    #[test]
    fn test_gaussian_basis_peaks_near_center() {
        let basis = GaussianBasis::new(1.0, 5.0, 5, 0.5, 6.0);
        let vals = basis.evaluate(1.0);
        assert!(vals[0] > 0.5, "Basis peak at first center: {}", vals[0]);
    }
    #[test]
    fn test_schnet_total_energy_finite() {
        let pot = SchNetPotential::new(4, 2, 4, 5.0);
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0]];
        let e = pot.total_energy(&positions);
        assert!(e.is_finite(), "SchNet energy = {e}");
    }
    #[test]
    fn test_schnet_interaction_output_size() {
        let basis = GaussianBasis::new(0.5, 5.0, 8, 0.5, 6.0);
        let layer = SchNetInteraction::new(4, basis);
        let features = vec![vec![0.1_f64; 4]; 3];
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]];
        let new_features = layer.forward(&features, &positions);
        assert_eq!(new_features.len(), 3);
        assert_eq!(new_features[0].len(), 4);
    }
    #[test]
    fn test_collect_triplets_known() {
        let positions = [[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triplets = collect_triplets(&positions, 2.0);
        assert!(!triplets.is_empty(), "Should find triplets");
    }
    #[test]
    fn test_triplet_cos_angle_right_angle() {
        let positions = [[1.0, 0.0, 0.0_f64], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let triplet = AngularTriplet { i: 1, j: 0, k: 2 };
        let cos_a = triplet_cos_angle(&positions, triplet);
        assert!(cos_a.abs() < 1e-10, "90-deg angle → cos=0, got {cos_a}");
    }
    #[test]
    fn test_angular_fourier_basis_length() {
        let basis = angular_fourier_basis(0.5, 6);
        assert_eq!(basis.len(), 6);
    }
    #[test]
    fn test_angular_fourier_basis_first_term_one() {
        let basis = angular_fourier_basis(0.3, 4);
        assert!(
            (basis[0] - 1.0).abs() < 1e-12,
            "cos(0) = 1, got {}",
            basis[0]
        );
    }
    #[test]
    fn test_dimenet_message_length() {
        let layer = DimeNetMessageLayer::new(6, 4, 5.0);
        let positions = [[0.0, 0.0, 0.0_f64], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]];
        let msg = layer.message(&positions, 0, 1);
        assert_eq!(msg.len(), 4, "message length should be n_out");
    }
    #[test]
    fn test_gat_total_energy_finite() {
        let pot = GatPotential::new(4, 2, 5.0);
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0]];
        let e = pot.total_energy(&positions);
        assert!(e.is_finite(), "GAT energy = {e}");
    }
    #[test]
    fn test_gat_layer_output_size() {
        let layer = GraphAttentionLayer::new(4, 4, 5.0);
        let features = vec![vec![0.1_f64; 4]; 2];
        let positions = [[0.0, 0.0, 0.0_f64], [2.0, 0.0, 0.0]];
        let new_feat = layer.forward(&features, &positions);
        assert_eq!(new_feat.len(), 2);
        assert_eq!(new_feat[0].len(), 4);
    }
    #[test]
    fn test_pes_stitching_blended_energy_finite() {
        let pes = PesStiching::new(vec![0.0, 5.0], vec![1.0, 1.0], vec![-10.0, -5.0]);
        let e = pes.blended_energy(2.5);
        assert!(e.is_finite(), "blended energy = {e}");
    }
    #[test]
    fn test_pes_stitching_weight_sum_to_one() {
        let pes = PesStiching::new(
            vec![0.0, 3.0, 6.0],
            vec![1.0, 1.0, 1.0],
            vec![0.0, 0.0, 0.0],
        );
        let total: f64 = (0..3).map(|k| pes.normalised_weight(k, 1.5)).sum();
        assert!((total - 1.0).abs() < 1e-10, "weights sum to 1: {total}");
    }
    #[test]
    fn test_pes_stitching_gradient_finite() {
        let pes = PesStiching::new(vec![0.0, 5.0], vec![2.0, 2.0], vec![-8.0, -3.0]);
        let g = pes.blended_energy_gradient(2.5);
        assert!(g.is_finite(), "gradient = {g}");
    }
    #[test]
    fn test_ecv_empty_ok() {
        let v = EnergyConservationValidator::new(0.01);
        assert!(v.validate().is_ok());
    }
    #[test]
    fn test_ecv_conserved_ok() {
        let mut v = EnergyConservationValidator::new(0.01);
        for _ in 0..10 {
            v.record(-100.0);
        }
        assert!(v.validate().is_ok(), "Constant energy should pass");
    }
    #[test]
    fn test_ecv_drift_detected() {
        let mut v = EnergyConservationValidator::new(0.001);
        v.record(-100.0);
        v.record(-95.0);
        assert!(v.validate().is_err(), "5% drift should fail 0.1% tolerance");
    }
    #[test]
    fn test_ecv_max_drift_nonneg() {
        let mut v = EnergyConservationValidator::new(0.01);
        v.record(-50.0);
        v.record(-51.0);
        assert!(v.max_drift() >= 0.0);
    }
    #[test]
    fn test_ecv_reset_clears() {
        let mut v = EnergyConservationValidator::new(0.01);
        v.record(-100.0);
        v.record(-200.0);
        v.reset();
        assert_eq!(v.n_steps(), 0);
    }
    #[test]
    fn test_ecv_n_steps() {
        let mut v = EnergyConservationValidator::new(0.01);
        for i in 0..7 {
            v.record(i as f64);
        }
        assert_eq!(v.n_steps(), 7);
    }
}
