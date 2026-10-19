// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Unit tests of the neural soft body module.

use super::*;

// -----------------------------------------------------------------------
// NeuralConstitutiveModel tests
// -----------------------------------------------------------------------

#[test]
fn test_constitutive_forward_zero_input() {
    let model = NeuralConstitutiveModel::new(16);
    let strain = [0.0f64; 6];
    let stress = model.forward(&strain);
    assert_eq!(stress.len(), 6);
    for s in &stress {
        assert!(
            s.abs() < 1e-12,
            "Zero strain must give zero stress (zero biases)"
        );
    }
}

#[test]
fn test_constitutive_forward_nonzero_weights() {
    let mut model = NeuralConstitutiveModel::linear_elastic(8, 1000.0);
    model.b2[0] = 0.5;
    let strain = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    let stress = model.forward(&strain);
    assert_eq!(stress.len(), 6);
    // With non-zero weights the stress should not all be zero
    let norm: f64 = stress.iter().map(|x| x * x).sum::<f64>().sqrt();
    assert!(
        norm > 0.0,
        "Non-zero weights should produce non-zero output"
    );
}

#[test]
fn test_constitutive_backward_gradient_shape() {
    let model = NeuralConstitutiveModel::new(8);
    let strain = [0.01, 0.02, -0.01, 0.0, 0.005, -0.003];
    let target = [100.0, 50.0, -20.0, 0.0, 10.0, -5.0];
    let (gw1, gb1, gw2, gb2) = model.backward(&strain, &target);
    assert_eq!(gw1.len(), model.n_hidden * 6);
    assert_eq!(gb1.len(), model.n_hidden);
    assert_eq!(gw2.len(), 6 * model.n_hidden);
    assert_eq!(gb2.len(), 6);
}

#[test]
fn test_constitutive_sgd_reduces_loss() {
    let mut model = NeuralConstitutiveModel::linear_elastic(16, 1.0);
    let strain = [0.01, 0.0, 0.0, 0.0, 0.0, 0.0];
    let target = [10.0, 5.0, 5.0, 0.0, 0.0, 0.0];

    let initial_loss: f64 = {
        let pred = model.forward(&strain);
        pred.iter()
            .zip(target.iter())
            .map(|(p, t)| (p - t).powi(2))
            .sum()
    };

    for _ in 0..200 {
        let (gw1, gb1, gw2, gb2) = model.backward(&strain, &target);
        model.sgd_step(&gw1, &gb1, &gw2, &gb2, 1e-3);
    }

    let final_loss: f64 = {
        let pred = model.forward(&strain);
        pred.iter()
            .zip(target.iter())
            .map(|(p, t)| (p - t).powi(2))
            .sum()
    };

    assert!(
        final_loss < initial_loss,
        "SGD should reduce loss: {initial_loss} -> {final_loss}"
    );
}

// -----------------------------------------------------------------------
// PhysicsInformedLoss tests
// -----------------------------------------------------------------------

#[test]
fn test_pinn_loss_zero_when_perfect() {
    let loss = PhysicsInformedLoss::new(0.0, 1e-5);
    let s = [1.0, 2.0, 3.0, 0.5, -0.5, 0.1];
    let l = loss.data_loss(&s, &s);
    assert!(l.abs() < 1e-12, "Loss should be zero when pred == target");
}

#[test]
fn test_pinn_loss_positive() {
    let loss = PhysicsInformedLoss::new(0.1, 1e-5);
    let pred = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let target = [0.0; 6];
    let l = loss.data_loss(&pred, &target);
    // MSE over six components: 1 / 6.
    assert!(
        (l - 1.0 / 6.0).abs() < 1e-15,
        "data MSE of a unit error in one of six components is 1/6, got {l}"
    );
}

#[test]
fn test_data_loss_is_independent_of_lambda_and_empty_input_is_zero() {
    let pred = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let target = [0.0; 6];
    let la = PhysicsInformedLoss::new(0.0, 1e-5).data_loss(&pred, &target);
    let lb = PhysicsInformedLoss::new(1000.0, 1e-5).data_loss(&pred, &target);
    assert_eq!(la, lb, "the data term must not contain the PDE weight");
    // Empty slices used to give 0/0 = NaN.
    let empty: [f64; 0] = [];
    let l = PhysicsInformedLoss::new(1.0, 1e-5).data_loss(&empty, &empty);
    assert_eq!(l, 0.0);
}

// -----------------------------------------------------------------------
// LearnedCollisionResponse tests
// -----------------------------------------------------------------------

#[test]
fn test_collision_response_output_shape() {
    let responder = LearnedCollisionResponse::new(8);
    let normal = [0.0, 1.0, 0.0];
    let impulse = responder.predict(&normal, -0.5, 0.01);
    assert_eq!(impulse.len(), 3);
}

#[test]
fn test_collision_response_zero_weights() {
    let responder = LearnedCollisionResponse::new(8);
    let normal = [0.0, 1.0, 0.0];
    let impulse = responder.predict(&normal, -0.5, 0.01);
    for x in &impulse {
        assert!(x.abs() < 1e-12, "Zero weights should give zero impulse");
    }
}

// -----------------------------------------------------------------------
// MaterialParameterIdentifier tests
// -----------------------------------------------------------------------

#[test]
fn test_material_id_lame_params() {
    let id = MaterialParameterIdentifier::new(1e6, 0.3, 1e-3, 0.0);
    let lam = id.lame_lambda();
    let mu = id.shear_modulus();
    assert!(lam > 0.0, "Lamé lambda should be positive");
    assert!(mu > 0.0, "Shear modulus should be positive");
}

#[test]
fn test_material_id_poisson_clamped() {
    // Poisson ratio should stay in [0, 0.499]
    let mut id = MaterialParameterIdentifier::new(1e6, 0.29, 1.0, 0.0);
    let sim_fn = |_e: f64, _nu: f64| -> Vec<[f64; 3]> { vec![[0.0; 3]; 1] };
    let obs = vec![[0.001, 0.0, 0.0]];
    id.step(&obs, &obs, &sim_fn);
    assert!(id.poisson_ratio >= 0.0 && id.poisson_ratio <= 0.499);
}

#[test]
fn test_material_id_young_positive() {
    let mut id = MaterialParameterIdentifier::new(1e6, 0.25, 1e-3, 0.0);
    let sim_fn = |_e: f64, _nu: f64| -> Vec<[f64; 3]> { vec![[0.0; 3]; 4] };
    let obs = vec![[0.001, 0.002, 0.0]; 4];
    for _ in 0..10 {
        id.step(&obs, &obs, &sim_fn);
    }
    assert!(
        id.young_modulus >= 1.0,
        "Young's modulus should stay positive"
    );
}

// -----------------------------------------------------------------------
// NeuralElasticityGraph tests
// -----------------------------------------------------------------------

#[test]
fn test_neural_graph_creation() {
    let nodes = vec![
        NeuralElasticNode::new([0.0, 0.0, 0.0], 4, 1.0),
        NeuralElasticNode::new([1.0, 0.0, 0.0], 4, 1.0),
    ];
    let edges = vec![NeuralElasticEdge::new(0, 1, 1.0)];
    let graph = NeuralElasticityGraph::new(nodes, edges, 4, 8);
    assert_eq!(graph.nodes.len(), 2);
    assert_eq!(graph.edges.len(), 1);
}

#[test]
fn test_neural_graph_zero_forces_with_zero_weights() {
    let nodes = vec![
        NeuralElasticNode::new([0.0, 0.0, 0.0], 4, 1.0),
        NeuralElasticNode::new([1.0, 0.0, 0.0], 4, 1.0),
        NeuralElasticNode::new([0.5, 1.0, 0.0], 4, 1.0),
    ];
    let edges = vec![
        NeuralElasticEdge::new(0, 1, 1.0),
        NeuralElasticEdge::new(1, 2, 1.0),
        NeuralElasticEdge::new(0, 2, 1.0),
    ];
    let graph = NeuralElasticityGraph::new(nodes, edges, 4, 8);
    let forces = graph.predict_forces();
    assert_eq!(forces.len(), 3);
    for f in &forces {
        for &c in f {
            assert!(c.abs() < 1e-12, "Zero weights should yield zero forces");
        }
    }
}

// -----------------------------------------------------------------------
// TrainingDataset tests
// -----------------------------------------------------------------------

#[test]
fn test_dataset_push_and_evict() {
    let mut ds = TrainingDataset::new(3);
    for i in 0..5u64 {
        ds.push(StrainStressSample {
            strain: [i as f64; 6],
            stress: [0.0; 6],
        });
    }
    assert_eq!(ds.len(), 3, "Dataset should not exceed max_size=3");
    assert_eq!(
        ds.samples[0].strain[0], 2.0,
        "Oldest samples should be evicted"
    );
}

#[test]
fn test_dataset_generate_linear_elastic() {
    let mut rng = rand::rng();
    let mut ds = TrainingDataset::new(100);
    ds.generate_linear_elastic(1e6, 0.3, 50, &mut rng);
    assert_eq!(ds.len(), 50);
    // Stress should not be all-zero for non-trivial strains
    let has_nonzero = ds
        .samples
        .iter()
        .any(|s| s.stress.iter().any(|&x| x != 0.0));
    assert!(has_nonzero, "Generated stress should be non-zero");
}

#[test]
fn test_dataset_is_empty() {
    let ds = TrainingDataset::new(10);
    assert!(ds.is_empty());
}

// -----------------------------------------------------------------------
// InferenceTimeCorrector tests
// -----------------------------------------------------------------------

#[test]
fn test_corrector_converges_for_identity_k() {
    let corrector = InferenceTimeCorrector::new(0.5, 100, 1e-8);
    let mut u = vec![2.0, -1.0, 0.5];
    let f = vec![0.0, 0.0, 0.0];
    // K = I => K*u = u, residual = u - f = u => correct to zero
    let iters = corrector.correct(&mut u, &f, &|x: &[f64]| x.to_vec());
    assert!(iters < 100, "Should converge before max iterations");
    for &ui in &u {
        assert!(ui.abs() < 1e-5, "u should converge to 0 = f when K=I");
    }
}

#[test]
fn test_corrector_already_converged() {
    let corrector = InferenceTimeCorrector::new(0.5, 100, 1e-8);
    let mut u = vec![0.0, 0.0, 0.0];
    let f = vec![0.0, 0.0, 0.0];
    let iters = corrector.correct(&mut u, &f, &|x: &[f64]| x.to_vec());
    assert_eq!(iters, 0, "Already converged should return 0 iterations");
}

// -----------------------------------------------------------------------
// NeuralPoseDeformer tests
// -----------------------------------------------------------------------

#[test]
fn test_pose_deformer_output_size() {
    let deformer = NeuralPoseDeformer::new(10, 4, 8);
    let weights = [0.25, 0.25, 0.25, 0.25];
    let result = deformer.deform(&weights);
    assert_eq!(result.len(), 30, "Output should be n_verts * 3");
}

#[test]
fn test_pose_deformer_zero_weights_zero_output() {
    let deformer = NeuralPoseDeformer::new(5, 3, 4);
    let weights = [0.0, 0.0, 0.0];
    let result = deformer.deform(&weights);
    for x in &result {
        assert!(
            x.abs() < 1e-12,
            "Zero pose weights should yield zero deformation"
        );
    }
}

// -----------------------------------------------------------------------
// LatentSpaceDeformer tests
// -----------------------------------------------------------------------

#[test]
fn test_latent_encode_output_size() {
    let lsd = LatentSpaceDeformer::new(5, 3);
    let pos = vec![0.0f64; 15];
    let z = lsd.encode(&pos);
    assert_eq!(z.len(), 3);
}

#[test]
fn test_latent_decode_output_size() {
    let lsd = LatentSpaceDeformer::new(5, 3);
    let z = vec![0.0f64; 3];
    let out = lsd.decode(&z);
    assert_eq!(out.len(), 15);
}

#[test]
fn test_latent_reconstruction_loss_zero_weights() {
    let lsd = LatentSpaceDeformer::new(4, 2);
    let pos = vec![
        1.0f64, 2.0, 3.0, -1.0, 0.5, 0.0, 0.1, -0.2, 2.0, 1.5, 0.0, 0.0,
    ];
    // With zero encoder weights encode gives tanh(0)=0, decode gives 0
    // Loss = MSE(pos, 0) > 0
    let loss = lsd.reconstruction_loss(&pos);
    assert!(
        loss > 0.0,
        "Reconstruction loss should be positive with random positions"
    );
}

// -----------------------------------------------------------------------
// RealTimeNeuralSoftBody tests
// -----------------------------------------------------------------------

#[test]
fn test_realtime_body_step_moves_dynamic_vertex() {
    let verts = vec![
        NeuralVertex::new_static([0.0, 0.0, 0.0]),
        NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
    ];
    let elems = vec![[0, 1, 2, 3]];
    let model = NeuralConstitutiveModel::new(8);
    let mut body = RealTimeNeuralSoftBody::new(verts, elems, model, 1000.0, 0.01);
    let gravity = [0.0, -9.81, 0.0];
    body.step(0.01, &gravity);
    // Dynamic vertex 1 should have moved in y
    assert!(
        body.vertices[1].position[1] < 0.0,
        "Dynamic vertex should fall under gravity"
    );
}

#[test]
fn test_realtime_body_static_vertex_does_not_move() {
    let verts = vec![
        NeuralVertex::new_static([0.5, 0.5, 0.5]),
        NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
    ];
    let elems = vec![[0, 1, 2, 3]];
    let model = NeuralConstitutiveModel::new(8);
    let mut body = RealTimeNeuralSoftBody::new(verts, elems, model, 1000.0, 0.01);
    let gravity = [0.0, -9.81, 0.0];
    for _ in 0..10 {
        body.step(0.01, &gravity);
    }
    assert_eq!(
        body.vertices[0].position,
        [0.5, 0.5, 0.5],
        "Static vertex should not move"
    );
}

// -----------------------------------------------------------------------
// NeuralPinnSolver tests
// -----------------------------------------------------------------------

#[test]
fn test_pinn_solver_train_step_returns_finite_loss() {
    let mut solver = NeuralPinnSolver::new(8, 0.01, 1e-4);
    let batch = vec![
        StrainStressSample {
            strain: [0.01, 0.0, 0.0, 0.0, 0.0, 0.0],
            stress: [1000.0, 300.0, 300.0, 0.0, 0.0, 0.0],
        },
        StrainStressSample {
            strain: [0.0, 0.02, 0.0, 0.0, 0.0, 0.0],
            stress: [300.0, 2000.0, 300.0, 0.0, 0.0, 0.0],
        },
    ];
    let loss = solver
        .train_step(&batch, [0.0; 3])
        .expect("independent samples never fail");
    assert!(loss.is_finite(), "Training loss should be finite");
    assert!(loss >= 0.0, "Training loss should be non-negative");
}

#[test]
fn test_pinn_solver_iteration_increments() {
    let mut solver = NeuralPinnSolver::new(4, 0.0, 1e-4);
    let batch = vec![StrainStressSample {
        strain: [0.0; 6],
        stress: [0.0; 6],
    }];
    for i in 1..=5 {
        solver
            .train_step(&batch, [0.0; 3])
            .expect("independent samples never fail");
        assert_eq!(solver.iteration, i);
    }
}

// -----------------------------------------------------------------------
// NeuralSoftBodySimulator tests
// -----------------------------------------------------------------------

#[test]
fn test_simulator_creation() {
    let verts = vec![
        NeuralVertex::new_static([0.0, 0.0, 0.0]),
        NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
    ];
    let elems = vec![[0, 1, 2, 3]];
    let sim = NeuralSoftBodySimulator::new(verts, elems, 8, 1000.0, 0.01, [0.0, -9.81, 0.0]);
    assert_eq!(sim.positions().len(), 4);
    assert_eq!(sim.rest_positions().len(), 4);
}

#[test]
fn test_simulator_step_no_panic() {
    let verts = vec![
        NeuralVertex::new_static([0.0, 0.0, 0.0]),
        NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
    ];
    let elems = vec![[0, 1, 2, 3]];
    let mut sim = NeuralSoftBodySimulator::new(verts, elems, 8, 1000.0, 0.01, [0.0, -9.81, 0.0]);
    for _ in 0..20 {
        let trained = sim.step(0.01, false).expect("stepping without training");
        assert!(trained.is_none());
    }
    let pos = sim.positions();
    assert_eq!(pos.len(), 4);
    // Positions should be finite
    for p in &pos {
        for &c in p {
            assert!(c.is_finite(), "Positions should remain finite");
        }
    }
}

// -----------------------------------------------------------------------
// Utility tests
// -----------------------------------------------------------------------

#[test]
fn test_vec3_norm() {
    let v = [3.0, 4.0, 0.0];
    assert!((vec3_norm(&v) - 5.0).abs() < 1e-12);
}

#[test]
fn test_vec3_normalise() {
    let v = [1.0, 0.0, 0.0];
    let n = vec3_normalise(&v);
    assert!((vec3_norm(&n) - 1.0).abs() < 1e-12);
}

#[test]
fn test_vec3_normalise_zero() {
    let v = [0.0, 0.0, 0.0];
    let n = vec3_normalise(&v);
    assert_eq!(n, [0.0, 0.0, 0.0]);
}

#[test]
fn test_vec3_cross() {
    let a = [1.0, 0.0, 0.0];
    let b = [0.0, 1.0, 0.0];
    let c = vec3_cross(&a, &b);
    assert!((c[0]).abs() < 1e-12);
    assert!((c[1]).abs() < 1e-12);
    assert!((c[2] - 1.0).abs() < 1e-12);
}

#[test]
fn test_vec3_dot() {
    let a = [1.0, 2.0, 3.0];
    let b = [4.0, 5.0, 6.0];
    assert!((vec3_dot(&a, &b) - 32.0).abs() < 1e-12);
}

#[test]
fn test_dense_forward_identity() {
    // 2×2 identity matrix
    let w = [1.0, 0.0, 0.0, 1.0];
    let b = [0.0, 0.0];
    let x = [3.0, 7.0];
    let y = dense_forward(&w, &b, &x, 2, 2);
    assert!((y[0] - 3.0).abs() < 1e-12);
    assert!((y[1] - 7.0).abs() < 1e-12);
}

#[test]
fn test_element_strain_rest_gives_zero() {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let rest = positions.clone();
    let elem = [0, 1, 2, 3];
    let strain = RealTimeNeuralSoftBody::element_strain(&positions, &rest, &elem);
    for &s in &strain {
        assert!(s.abs() < 1e-12, "No deformation => zero strain, got {s}");
    }
}

// ===========================================================================
// Regression tests: F4230 (equilibrium loss), F4231 (element forces) and
// F4232 (training reaches the body)
// ===========================================================================

/// Deterministic generator for test data.
struct Lcg(u64);

impl Lcg {
    fn unit(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn symmetric(&mut self) -> f64 {
        2.0 * self.unit() - 1.0
    }
}

fn linear_elastic_stress(young: f64, nu: f64, e: &[f64; 6]) -> [f64; 6] {
    let mu = young / (2.0 * (1.0 + nu));
    let lambda = young * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let trace = e[0] + e[1] + e[2];
    [
        lambda * trace + 2.0 * mu * e[0],
        lambda * trace + 2.0 * mu * e[1],
        lambda * trace + 2.0 * mu * e[2],
        2.0 * mu * e[3],
        2.0 * mu * e[4],
        2.0 * mu * e[5],
    ]
}

fn random_samples(
    count: usize,
    young: f64,
    nu: f64,
    amplitude: f64,
    seed: u64,
) -> Vec<StrainStressSample> {
    let mut rng = Lcg(seed);
    (0..count)
        .map(|_| {
            let mut strain = [0.0f64; 6];
            for s in strain.iter_mut() {
                *s = amplitude * rng.symmetric();
            }
            StrainStressSample {
                strain,
                stress: linear_elastic_stress(young, nu, &strain),
            }
        })
        .collect()
}

/// Sample a stress field on a regular grid (x fastest) with spacing `h`.
fn grid_field(dims: [usize; 3], h: f64, field: impl Fn([f64; 3]) -> [f64; 6]) -> Vec<[f64; 6]> {
    let mut out = Vec::with_capacity(dims[0] * dims[1] * dims[2]);
    for k in 0..dims[2] {
        for j in 0..dims[1] {
            for i in 0..dims[0] {
                out.push(field([i as f64 * h, j as f64 * h, k as f64 * h]));
            }
        }
    }
    out
}

fn random_stress_field(count: usize, seed: u64) -> Vec<[f64; 6]> {
    let mut rng = Lcg(seed);
    (0..count)
        .map(|_| {
            let mut s = [0.0f64; 6];
            for c in s.iter_mut() {
                *c = rng.symmetric();
            }
            s
        })
        .collect()
}

fn assert_vec3_close(actual: [f64; 3], expected: [f64; 3], tol: f64, what: &str) {
    for (a, e) in actual.iter().zip(expected.iter()) {
        assert!((a - e).abs() <= tol, "{what}: {actual:?} vs {expected:?}");
    }
}

// ---- F4230: equilibrium (PDE) loss ----------------------------------------

#[test]
fn pde_loss_of_a_constant_stress_field_vanishes() {
    let loss = PhysicsInformedLoss::new(1.0, 0.1);
    let dims = [5, 4, 3];
    let stress = vec![[1.0, 2.0, 3.0, 0.5, -0.5, 0.25]; 60];
    assert_eq!(loss.pde_loss(&stress, dims, [0.0; 3]).unwrap(), 0.0);
    let residuals = loss.equilibrium_residuals(&stress, dims, [0.0; 3]).unwrap();
    assert_eq!(residuals.len(), 3 * 2);
    assert!(residuals.iter().all(|r| *r == [0.0; 3]));
    // A constant stress cannot balance a body force: the residual is f itself.
    let f = [0.0, 0.0, -9.81];
    let residuals = loss.equilibrium_residuals(&stress, dims, f).unwrap();
    assert!(residuals.iter().all(|r| *r == f));
    let l = loss.pde_loss(&stress, dims, f).unwrap();
    assert!((l - 9.81 * 9.81).abs() < 1e-12);
}

#[test]
fn pde_loss_matches_the_analytic_divergence_of_a_linear_field() {
    // sigma_xx = 2x + 3y, sigma_yy = -2y + z, sigma_zz = z/2 + x,
    // sigma_xy = x - 4y + z, sigma_xz = 5z - x, sigma_yz = y + 3z
    // => div(sigma) = (2 - 4 + 5, 1 - 2 + 3, -1 + 1 + 1/2) = (3, 2, 1/2).
    let h = 0.25;
    let dims = [4, 5, 3];
    let stress = grid_field(dims, h, |p| {
        let [x, y, z] = p;
        [
            2.0 * x + 3.0 * y,
            -2.0 * y + z,
            0.5 * z + x,
            x - 4.0 * y + z,
            5.0 * z - x,
            y + 3.0 * z,
        ]
    });
    let loss = PhysicsInformedLoss::new(1.0, h);
    let residuals = loss.equilibrium_residuals(&stress, dims, [0.0; 3]).unwrap();
    assert_eq!(residuals.len(), 2 * 3);
    for r in &residuals {
        assert_vec3_close(*r, [3.0, 2.0, 0.5], 1e-12, "divergence");
    }
    let l = loss.pde_loss(&stress, dims, [0.0; 3]).unwrap();
    assert!((l - (9.0 + 4.0 + 0.25)).abs() < 1e-11, "L_pde = {l}");
    // The matching body force balances it exactly.
    let balanced = loss.pde_loss(&stress, dims, [-3.0, -2.0, -0.5]).unwrap();
    assert!(balanced < 1e-22, "balanced field residual {balanced}");
}

#[test]
fn pde_loss_vanishes_for_hydrostatic_equilibrium_and_handles_plane_problems() {
    // d(sigma_zz)/dz + f_z = 0 with sigma_zz = 5 z and f_z = -5.
    let h = 0.2;
    let dims = [3, 3, 6];
    let stress = grid_field(dims, h, |p| [0.0, 0.0, 5.0 * p[2], 0.0, 0.0, 0.0]);
    let loss = PhysicsInformedLoss::new(1.0, h);
    assert!(loss.pde_loss(&stress, dims, [0.0, 0.0, -5.0]).unwrap() < 1e-22);
    assert!((loss.pde_loss(&stress, dims, [0.0; 3]).unwrap() - 25.0).abs() < 1e-10);

    // Plane problem (a single node along y): only x and z derivatives count.
    // sigma_xx = 2x, sigma_xz = 3z, sigma_zz = 4z => div = (2 + 3, 0, 0 + 4).
    let dims = [6, 1, 5];
    let stress = grid_field(dims, h, |p| {
        [2.0 * p[0], 0.0, 4.0 * p[2], 0.0, 3.0 * p[2], 0.0]
    });
    let residuals = loss.equilibrium_residuals(&stress, dims, [0.0; 3]).unwrap();
    assert_eq!(residuals.len(), 4 * 3);
    for r in &residuals {
        assert_vec3_close(*r, [5.0, 0.0, 4.0], 1e-12, "plane divergence");
    }
}

#[test]
fn pde_gradient_matches_finite_differences() {
    let dims = [5, 4, 4];
    let h = 0.3;
    let loss = PhysicsInformedLoss::new(1.0, h);
    let body_force = [0.2, -0.1, 0.4];
    let stress = random_stress_field(80, 3);
    let gradient = loss.pde_gradient(&stress, dims, body_force).unwrap();
    assert_eq!(gradient.len(), 80);
    let eps = 1e-6;
    let mut max_abs = 0.0f64;
    for node in 0..80 {
        for component in 0..6 {
            let mut plus = stress.clone();
            let mut minus = stress.clone();
            plus[node][component] += eps;
            minus[node][component] -= eps;
            let fd = (loss.pde_loss(&plus, dims, body_force).unwrap()
                - loss.pde_loss(&minus, dims, body_force).unwrap())
                / (2.0 * eps);
            let analytic = gradient[node][component];
            max_abs = max_abs.max(analytic.abs());
            assert!(
                (analytic - fd).abs() <= 1e-7 * (1.0 + fd.abs()),
                "node {node} component {component}: analytic {analytic} vs fd {fd}"
            );
        }
    }
    assert!(max_abs > 1e-3, "the gradient check must not be vacuous");
    // Boundary nodes of a 3-D grid with a stencil of +-1 only see interior
    // neighbours: corner nodes are never part of a stencil.
    assert_eq!(gradient[0], [0.0; 6]);
}

#[test]
fn compute_field_combines_data_mse_and_weighted_pde_loss() {
    let dims = [4, 4, 3];
    let h = 0.5;
    let pred = random_stress_field(48, 9);
    let target = random_stress_field(48, 10);
    let f = [0.1, 0.2, -0.3];
    let loss = PhysicsInformedLoss::new(0.7, h);
    let flat_pred: Vec<f64> = pred.iter().flatten().copied().collect();
    let flat_target: Vec<f64> = target.iter().flatten().copied().collect();
    let expected =
        loss.data_loss(&flat_pred, &flat_target) + 0.7 * loss.pde_loss(&pred, dims, f).unwrap();
    let total = loss.compute_field(&pred, dims, &target, f).unwrap();
    assert!((total - expected).abs() < 1e-14 * (1.0 + expected.abs()));
    // A perfect, equilibrated prediction has zero loss.
    let constant = vec![[1.0, 2.0, 3.0, 0.1, 0.2, 0.3]; 48];
    assert_eq!(
        loss.compute_field(&constant, dims, &constant, [0.0; 3])
            .unwrap(),
        0.0
    );
    // lambda scales only the PDE part.
    let heavy = PhysicsInformedLoss::new(7.0, h);
    let heavy_total = heavy.compute_field(&pred, dims, &target, f).unwrap();
    let pde = loss.pde_loss(&pred, dims, f).unwrap();
    assert!((heavy_total - total - (7.0 - 0.7) * pde).abs() < 1e-12 * (1.0 + heavy_total));
}

#[test]
fn pde_operators_reject_invalid_grids() {
    let loss = PhysicsInformedLoss::new(1.0, 0.1);
    let stress = vec![[0.0f64; 6]; 24];
    // Control: 4 x 6 x 1 is a valid plane grid of 24 nodes.
    let valid = [4, 6, 1];
    assert!(loss.pde_loss(&stress, valid, [0.0; 3]).is_ok());
    // Length does not match the grid.
    assert!(loss.pde_loss(&stress, [3, 3, 3], [0.0; 3]).is_err());
    assert!(loss.pde_gradient(&stress, [5, 5, 5], [0.0; 3]).is_err());
    // An axis with exactly two nodes has no interior node.
    assert!(loss.pde_loss(&stress, [3, 2, 4], [0.0; 3]).is_err());
    // Empty grid.
    assert!(loss.pde_loss(&[], [0, 1, 1], [0.0; 3]).is_err());
    // Invalid spacing on an otherwise valid grid.
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let broken = PhysicsInformedLoss::new(1.0, bad);
        assert!(broken.pde_loss(&stress, valid, [0.0; 3]).is_err());
        assert!(broken.pde_gradient(&stress, valid, [0.0; 3]).is_err());
        assert!(
            broken
                .equilibrium_residuals(&stress, valid, [0.0; 3])
                .is_err()
        );
    }
    // Target of a different size.
    assert!(
        loss.compute_field(&stress, valid, &stress[..5], [0.0; 3])
            .is_err()
    );
    match loss.pde_loss(&stress, [3, 3, 3], [0.0; 3]) {
        Err(Error::InvalidParameter(message)) => assert!(message.contains("24")),
        other => panic!("expected InvalidParameter, got {other:?}"),
    }
}

fn total_loss(
    model: &NeuralConstitutiveModel,
    loss: &PhysicsInformedLoss,
    samples: &[StrainStressSample],
    dims: [usize; 3],
    body_force: [f64; 3],
) -> f64 {
    let predictions: Vec<[f64; 6]> = samples
        .iter()
        .map(|s| model.forward_voigt(&s.strain))
        .collect();
    let targets: Vec<[f64; 6]> = samples.iter().map(|s| s.stress).collect();
    loss.compute_field(&predictions, dims, &targets, body_force)
        .expect("valid grid")
}

fn flatten_parameters(model: &NeuralConstitutiveModel) -> Vec<f64> {
    let mut params = Vec::new();
    params.extend_from_slice(&model.w1);
    params.extend_from_slice(&model.b1);
    params.extend_from_slice(&model.w2);
    params.extend_from_slice(&model.b2);
    params
}

fn with_parameter_offset(
    model: &NeuralConstitutiveModel,
    index: usize,
    delta: f64,
) -> NeuralConstitutiveModel {
    let mut out = model.clone();
    let (n1, nb1, n2) = (out.w1.len(), out.b1.len(), out.w2.len());
    if index < n1 {
        out.w1[index] += delta;
    } else if index < n1 + nb1 {
        out.b1[index - n1] += delta;
    } else if index < n1 + nb1 + n2 {
        out.w2[index - n1 - nb1] += delta;
    } else {
        out.b2[index - n1 - nb1 - n2] += delta;
    }
    out
}

#[test]
fn train_step_gradient_includes_the_pde_term_and_matches_finite_differences() {
    let dims = [4, 4, 3];
    let h = 0.25;
    let body_force = [0.3, -0.2, 0.5];
    let lr = 0.01;
    let lambda = 0.5;
    let samples = random_samples(48, 5.0, 0.3, 0.5, 21);
    let mut solver = NeuralPinnSolver::new(6, lambda, lr);
    solver.pinn_loss = PhysicsInformedLoss::new(lambda, h).with_grid(dims);
    let initial = solver.model.clone();

    let reported = solver
        .train_step(&samples, body_force)
        .expect("valid grid batch");
    let expected = total_loss(&initial, &solver.pinn_loss, &samples, dims, body_force);
    assert!(
        (reported - expected).abs() < 1e-12 * (1.0 + expected.abs()),
        "reported loss {reported} vs recomputed {expected}"
    );
    // The PDE part is a significant share of the loss.
    let predictions: Vec<[f64; 6]> = samples
        .iter()
        .map(|s| initial.forward_voigt(&s.strain))
        .collect();
    let pde = solver
        .pinn_loss
        .pde_loss(&predictions, dims, body_force)
        .unwrap();
    assert!(lambda * pde > 1e-3, "PDE share {}", lambda * pde);

    // SGD applied w -= lr * grad, so the applied gradient is (before - after) / lr.
    let before = flatten_parameters(&initial);
    let after = flatten_parameters(&solver.model);
    let eps = 1e-6;
    let mut checked_nonzero = 0;
    for index in 0..before.len() {
        let analytic = (before[index] - after[index]) / lr;
        let plus = total_loss(
            &with_parameter_offset(&initial, index, eps),
            &solver.pinn_loss,
            &samples,
            dims,
            body_force,
        );
        let minus = total_loss(
            &with_parameter_offset(&initial, index, -eps),
            &solver.pinn_loss,
            &samples,
            dims,
            body_force,
        );
        let fd = (plus - minus) / (2.0 * eps);
        assert!(
            (analytic - fd).abs() <= 1e-5 * (1.0 + fd.abs()),
            "parameter {index}: analytic {analytic} vs finite difference {fd}"
        );
        if fd.abs() > 1e-6 {
            checked_nonzero += 1;
        }
    }
    assert!(checked_nonzero > before.len() / 2);
}

#[test]
fn lambda_pde_changes_the_trained_weights_on_a_grid() {
    let dims = [4, 4, 3];
    let samples = random_samples(48, 5.0, 0.3, 0.5, 22);
    let make = |lambda: f64, grid: bool| {
        let mut solver = NeuralPinnSolver::new(6, lambda, 0.01);
        let loss = PhysicsInformedLoss::new(lambda, 0.25);
        solver.pinn_loss = if grid { loss.with_grid(dims) } else { loss };
        solver
    };
    let mut data_only = make(0.0, true);
    let mut with_pde = make(1.0, true);
    let mut no_grid = make(1.0, false);
    let loss_a = data_only.train_step(&samples, [0.0; 3]).unwrap();
    let loss_b = with_pde.train_step(&samples, [0.0; 3]).unwrap();
    let loss_c = no_grid.train_step(&samples, [0.0; 3]).unwrap();
    assert!(loss_b > loss_a, "the PDE term adds a non-negative loss");
    assert_ne!(
        data_only.model, with_pde.model,
        "lambda_pde must reach the weights"
    );
    // Independent samples carry no spatial information: data term only.
    assert_eq!(loss_a, loss_c);
    assert_eq!(data_only.model, no_grid.model);
}

#[test]
fn pinn_term_reduces_the_equilibrium_residual_of_the_learned_field() {
    let dims = [6, 6, 1];
    let h = 0.5;
    let samples = random_samples(36, 1.0, 0.3, 0.5, 5);
    let train = |lambda: f64| {
        let mut solver = NeuralPinnSolver::new(8, lambda, 0.2);
        solver.pinn_loss = PhysicsInformedLoss::new(lambda, h).with_grid(dims);
        for _ in 0..200 {
            solver.train_step(&samples, [0.0; 3]).unwrap();
        }
        let predictions: Vec<[f64; 6]> = samples
            .iter()
            .map(|s| solver.model.forward_voigt(&s.strain))
            .collect();
        solver
            .pinn_loss
            .pde_loss(&predictions, dims, [0.0; 3])
            .unwrap()
    };
    let without = train(0.0);
    let with = train(1.0);
    assert!(
        without > 0.1,
        "a data-only fit leaves a large residual: {without}"
    );
    assert!(
        with < 0.1 * without,
        "lambda_pde = 1 must shrink the equilibrium residual: {with} vs {without}"
    );
}

#[test]
fn train_step_rejects_a_batch_that_does_not_match_the_grid() {
    let samples = random_samples(5, 1.0, 0.3, 0.1, 1);
    let mut solver = NeuralPinnSolver::new(4, 1.0, 0.1);
    solver.pinn_loss = PhysicsInformedLoss::new(1.0, 0.1).with_grid([3, 3, 3]);
    let before = solver.model.clone();
    assert!(solver.train_step(&samples, [0.0; 3]).is_err());
    assert_eq!(
        solver.model, before,
        "a rejected step must not touch the weights"
    );
    assert_eq!(solver.iteration, 0);
}

#[test]
fn train_step_on_an_empty_batch_is_a_no_op() {
    let mut solver = NeuralPinnSolver::new(4, 1.0, 0.1);
    let before = solver.model.clone();
    let loss = solver.train_step(&[], [0.0; 3]).unwrap();
    assert_eq!(loss, 0.0);
    assert_eq!(solver.model, before);
    assert_eq!(solver.iteration, 0);
}

// ---- F4232: initialisation and online training ---------------------------

#[test]
fn he_initialisation_is_deterministic_bounded_and_trainable() {
    let model = NeuralConstitutiveModel::new(8);
    assert_eq!(model, NeuralConstitutiveModel::new(8), "deterministic");
    assert_eq!(
        model,
        NeuralConstitutiveModel::with_seed(8, DEFAULT_INIT_SEED)
    );
    assert_ne!(model, NeuralConstitutiveModel::with_seed(8, 1));
    assert_ne!(
        NeuralConstitutiveModel::with_seed(8, 1),
        NeuralConstitutiveModel::with_seed(8, 2)
    );
    // He-uniform limits: sqrt(6 / fan_in).
    assert!(model.w1.iter().all(|w| w.abs() <= 1.0));
    let limit2 = (6.0f64 / 8.0).sqrt();
    assert!(model.w2.iter().all(|w| w.abs() <= limit2));
    assert!(model.w1.iter().any(|w| w.abs() > 0.1));
    assert!(model.w2.iter().any(|w| w.abs() > 0.1));
    assert!(model.b1.iter().chain(model.b2.iter()).all(|b| *b == 0.0));
    // Both layers use a symmetric distribution around zero.
    assert!(model.w1.iter().any(|w| *w > 0.0) && model.w1.iter().any(|w| *w < 0.0));
    // Zero strain still gives zero stress.
    assert_eq!(model.forward_voigt(&[0.0; 6]), [0.0; 6]);
}

#[test]
fn zero_initialised_model_has_no_hidden_layer_gradient_but_he_model_has() {
    let strain = [0.2, -0.1, 0.3, 0.05, -0.2, 0.1];
    let target = [1.0, -1.0, 0.5, 0.2, -0.3, 0.4];
    let zero = NeuralConstitutiveModel::zeros(8);
    let (gw1, gb1, gw2, gb2) = zero.backward(&strain, &target);
    assert!(
        gw1.iter()
            .chain(gb1.iter())
            .chain(gw2.iter())
            .all(|g| *g == 0.0)
    );
    assert!(gb2.iter().any(|g| *g != 0.0), "only the output bias learns");
    let he = NeuralConstitutiveModel::new(8);
    let (gw1, gb1, gw2, _gb2) = he.backward(&strain, &target);
    assert!(gw1.iter().any(|g| *g != 0.0));
    assert!(gb1.iter().any(|g| *g != 0.0));
    assert!(gw2.iter().any(|g| *g != 0.0));
}

#[test]
fn linear_elastic_model_is_exactly_linear_in_tension_and_compression() {
    let model = NeuralConstitutiveModel::linear_elastic(12, 7.5);
    for strain in [
        [0.01, -0.02, 0.03, -0.04, 0.05, -0.06],
        [-0.3, -0.1, 0.0, 0.2, -0.7, 0.4],
        [0.0; 6],
    ] {
        let stress = model.forward_voigt(&strain);
        for (s, e) in stress.iter().zip(strain.iter()) {
            assert!((s - 7.5 * e).abs() < 1e-14, "stress {s} vs {}", 7.5 * e);
        }
    }
    // Extra hidden units are inert; too few hidden units neither panic nor
    // produce stress for components without a unit pair.
    let wide = NeuralConstitutiveModel::linear_elastic(20, 2.0);
    assert_eq!(
        wide.forward_voigt(&[1.0, -1.0, 0.5, 0.0, 0.0, 0.0])[..3],
        [2.0, -2.0, 1.0]
    );
    let tiny = NeuralConstitutiveModel::linear_elastic(3, 2.0);
    assert_eq!(
        tiny.forward_voigt(&[1.0, -1.0, 0.5, 0.0, 0.0, 0.0]),
        [2.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    );
    let empty = NeuralConstitutiveModel::linear_elastic(0, 2.0);
    assert_eq!(empty.forward_voigt(&[1.0; 6]), [0.0; 6]);
    assert_eq!(
        NeuralConstitutiveModel::new(0).forward_voigt(&[1.0; 6]),
        [0.0; 6]
    );
}

fn tetrahedron_vertices() -> Vec<NeuralVertex> {
    vec![
        NeuralVertex::new_static([0.0, 0.0, 0.0]),
        NeuralVertex::new([1.0, 0.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 1.0, 0.0], 1.0),
        NeuralVertex::new([0.0, 0.0, 1.0], 1.0),
    ]
}

#[test]
fn online_training_reaches_the_body_and_reduces_the_loss() {
    let mut sim = NeuralSoftBodySimulator::new(
        tetrahedron_vertices(),
        vec![[0, 1, 2, 3]],
        8,
        1000.0,
        0.01,
        [0.0, -9.81, 0.0],
    );
    // The body starts with the solver's model.
    assert_eq!(sim.body.constitutive, sim.solver.model);
    for sample in random_samples(64, 1.0, 0.3, 0.5, 11) {
        sim.dataset.push(sample);
    }
    sim.solver.lr = 0.5;
    let initial_w1 = sim.body.constitutive.w1.clone();
    let mut losses = Vec::new();
    for _ in 0..100 {
        let loss = sim.step(0.01, true).expect("training step");
        losses.push(loss.expect("a non-empty dataset trains every step"));
    }
    assert!(
        losses[99] < 0.2 * losses[0],
        "loss must decrease: {} -> {}",
        losses[0],
        losses[99]
    );
    assert_ne!(
        sim.body.constitutive.w1, initial_w1,
        "the trained weights must reach the soft body"
    );
    assert_eq!(
        sim.body.constitutive, sim.solver.model,
        "model and body stay in sync"
    );
    assert_eq!(sim.solver.iteration, 100);
}

#[test]
fn simulator_without_training_leaves_the_model_alone_and_reports_grid_errors() {
    let mut sim = NeuralSoftBodySimulator::new(
        tetrahedron_vertices(),
        vec![[0, 1, 2, 3]],
        8,
        1000.0,
        0.01,
        [0.0, -9.81, 0.0],
    );
    for sample in random_samples(3, 1.0, 0.3, 0.1, 4) {
        sim.dataset.push(sample);
    }
    let model = sim.body.constitutive.clone();
    assert_eq!(sim.step(0.01, false).unwrap(), None);
    assert_eq!(sim.body.constitutive, model);
    // Empty dataset: nothing to train on.
    sim.dataset.samples.clear();
    assert_eq!(sim.step(0.01, true).unwrap(), None);
    // A grid layout that does not match the dataset is an error, not a panic.
    for sample in random_samples(3, 1.0, 0.3, 0.1, 4) {
        sim.dataset.push(sample);
    }
    sim.solver.pinn_loss = PhysicsInformedLoss::new(1.0, 0.1).with_grid([3, 3, 3]);
    assert!(sim.step(0.01, true).is_err());
}

#[test]
fn dataset_with_zero_capacity_ignores_samples() {
    let mut ds = TrainingDataset::new(0);
    ds.push(StrainStressSample {
        strain: [0.0; 6],
        stress: [0.0; 6],
    });
    assert!(ds.is_empty());
}

// ---- F4231: element strain and forces -------------------------------------

fn general_rest_tetrahedron() -> Vec<[f64; 3]> {
    vec![
        [0.1, 0.2, 0.3],
        [1.3, 0.1, -0.2],
        [0.4, 1.7, 0.5],
        [-0.3, 0.6, 2.1],
    ]
}

fn apply_linear_displacement(
    rest: &[[f64; 3]],
    gradient: [[f64; 3]; 3],
    translation: [f64; 3],
) -> Vec<[f64; 3]> {
    rest.iter()
        .map(|x| {
            let mut out = [0.0f64; 3];
            for (i, o) in out.iter_mut().enumerate() {
                *o = x[i]
                    + translation[i]
                    + gradient[i][0] * x[0]
                    + gradient[i][1] * x[1]
                    + gradient[i][2] * x[2];
            }
            out
        })
        .collect()
}

#[test]
fn element_strain_recovers_a_uniform_strain_on_any_tetrahedron() {
    let elem = [0, 1, 2, 3];
    // Symmetric part G and an arbitrary skew part W.
    let g = [
        [0.01, 0.004, -0.002],
        [0.004, -0.02, 0.003],
        [-0.002, 0.003, 0.015],
    ];
    let w = [[0.0, 0.05, -0.03], [-0.05, 0.0, 0.02], [0.03, -0.02, 0.0]];
    let mut full = [[0.0f64; 3]; 3];
    for ((row, g_row), w_row) in full.iter_mut().zip(g.iter()).zip(w.iter()) {
        for ((entry, gv), wv) in row.iter_mut().zip(g_row.iter()).zip(w_row.iter()) {
            *entry = gv + wv;
        }
    }
    let expected = [g[0][0], g[1][1], g[2][2], g[0][1], g[0][2], g[1][2]];
    for scale in [1.0, 2.0, 0.25] {
        let rest: Vec<[f64; 3]> = general_rest_tetrahedron()
            .iter()
            .map(|x| [x[0] * scale, x[1] * scale, x[2] * scale])
            .collect();
        let positions = apply_linear_displacement(&rest, full, [0.3, -0.2, 0.1]);
        let strain = RealTimeNeuralSoftBody::element_strain(&positions, &rest, &elem);
        for (s, e) in strain.iter().zip(expected.iter()) {
            assert!(
                (s - e).abs() < 1e-12,
                "scale {scale}: {strain:?} vs {expected:?}"
            );
        }
    }
}

#[test]
fn element_strain_of_a_rigid_translation_is_zero_and_degenerate_input_is_safe() {
    let rest = general_rest_tetrahedron();
    let elem = [0, 1, 2, 3];
    let moved: Vec<[f64; 3]> = rest
        .iter()
        .map(|x| [x[0] + 5.0, x[1] - 3.0, x[2] + 0.7])
        .collect();
    // Translating changes the coordinates by rounding noise only.
    let strain = RealTimeNeuralSoftBody::element_strain(&moved, &rest, &elem);
    assert!(strain.iter().all(|s| s.abs() < 1e-12), "strain {strain:?}");
    // Flat (coplanar) rest tetrahedron and out-of-range corners give zero strain.
    let flat = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
    ];
    assert_eq!(
        RealTimeNeuralSoftBody::element_strain(&flat, &flat, &elem),
        [0.0; 6]
    );
    assert_eq!(
        RealTimeNeuralSoftBody::element_strain(&rest, &rest, &[0, 1, 2, 9]),
        [0.0; 6]
    );
    assert!(TetRestShape::from_rest(&flat, &elem).is_none());
    let shape = TetRestShape::from_rest(&rest, &elem).expect("valid tetrahedron");
    assert!(shape.volume > 0.0);
}

fn body_from(
    rest: &[[f64; 3]],
    positions: &[[f64; 3]],
    model: NeuralConstitutiveModel,
) -> RealTimeNeuralSoftBody {
    let vertices: Vec<NeuralVertex> = rest.iter().map(|&x| NeuralVertex::new(x, 1.0)).collect();
    let mut body = RealTimeNeuralSoftBody::new(vertices, vec![[0, 1, 2, 3]], model, 1000.0, 0.0);
    for (v, p) in body.vertices.iter_mut().zip(positions.iter()) {
        v.position = *p;
    }
    body
}

#[test]
fn element_forces_sum_to_zero_for_a_random_deformation() {
    let rest = general_rest_tetrahedron();
    let mut rng = Lcg(77);
    for seed in 1..=5u64 {
        let positions: Vec<[f64; 3]> = rest
            .iter()
            .map(|x| {
                [
                    x[0] + 0.1 * rng.symmetric(),
                    x[1] + 0.1 * rng.symmetric(),
                    x[2] + 0.1 * rng.symmetric(),
                ]
            })
            .collect();
        let body = body_from(
            &rest,
            &positions,
            NeuralConstitutiveModel::with_seed(8, seed),
        );
        let forces = body.elastic_forces();
        let mut total = [0.0f64; 3];
        let mut largest = 0.0f64;
        for f in &forces {
            for (t, c) in total.iter_mut().zip(f.iter()) {
                *t += c;
                largest = largest.max(c.abs());
            }
        }
        assert!(largest > 1e-4, "seed {seed}: forces must be non-trivial");
        for t in total {
            assert!(
                t.abs() < 1e-12 * (1.0 + largest),
                "seed {seed}: net force {total:?}"
            );
        }
    }
}

#[test]
fn rigid_translation_gives_zero_force() {
    let rest = general_rest_tetrahedron();
    let moved: Vec<[f64; 3]> = rest
        .iter()
        .map(|x| [x[0] + 2.0, x[1] - 1.0, x[2] + 3.5])
        .collect();
    let body = body_from(&rest, &moved, NeuralConstitutiveModel::new(8));
    for f in body.elastic_forces() {
        assert_vec3_close(f, [0.0; 3], 1e-12, "a translated element is stress free");
    }
}

#[test]
fn compression_and_shear_of_a_unit_tetrahedron_give_the_analytic_nodal_forces() {
    let rest = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let stiffness = 100.0;
    let model = NeuralConstitutiveModel::linear_elastic(12, stiffness);
    let volume = 1.0 / 6.0;
    let delta = 0.01;

    // Vertex 3 pushed towards the base: sigma_zz = -k delta, which pushes the
    // apex back up and the base down.
    let mut squeezed = rest.clone();
    squeezed[3][2] -= delta;
    let forces = body_from(&rest, &squeezed, model.clone()).elastic_forces();
    assert_vec3_close(
        forces[3],
        [0.0, 0.0, volume * stiffness * delta],
        1e-12,
        "apex",
    );
    assert_vec3_close(
        forces[0],
        [0.0, 0.0, -volume * stiffness * delta],
        1e-12,
        "base",
    );
    assert_vec3_close(forces[1], [0.0; 3], 1e-12, "vertex 1");
    assert_vec3_close(forces[2], [0.0; 3], 1e-12, "vertex 2");

    // Vertex 3 sheared along x: eps_xz = delta / 2, which the old code dropped.
    let mut sheared = rest.clone();
    sheared[3][0] += delta;
    let forces = body_from(&rest, &sheared, model).elastic_forces();
    let shear_stress = stiffness * delta / 2.0;
    assert_vec3_close(
        forces[3],
        [-volume * shear_stress, 0.0, 0.0],
        1e-12,
        "sheared apex",
    );
    assert_vec3_close(
        forces[1],
        [0.0, 0.0, -volume * shear_stress],
        1e-12,
        "sheared vertex 1",
    );
    let mut net = [0.0f64; 3];
    for f in &forces {
        for (n, c) in net.iter_mut().zip(f.iter()) {
            *n += c;
        }
    }
    assert_vec3_close(net, [0.0; 3], 1e-12, "net shear force");
}

/// Elastic energy of a linear material `sigma = s * eps`:
/// `sum V0 / 2 * s * eps:eps` with tensor shear counted twice.
fn linear_energy(rest: &[[f64; 3]], positions: &[[f64; 3]], elems: &[[usize; 4]], s: f64) -> f64 {
    elems
        .iter()
        .map(|elem| {
            let strain = RealTimeNeuralSoftBody::element_strain(positions, rest, elem);
            let shape = TetRestShape::from_rest(rest, elem).expect("valid element");
            let norm2 = strain[0] * strain[0]
                + strain[1] * strain[1]
                + strain[2] * strain[2]
                + 2.0 * (strain[3] * strain[3] + strain[4] * strain[4] + strain[5] * strain[5]);
            0.5 * shape.volume * s * norm2
        })
        .sum()
}

#[test]
fn elastic_forces_are_the_negative_energy_gradient_of_a_linear_material() {
    // Two tetrahedra sharing the face (1, 2, 3).
    let rest = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.1, 0.0],
        [0.2, 1.1, 0.1],
        [0.1, 0.3, 0.9],
        [1.1, 1.2, 1.0],
    ];
    let elems = vec![[0, 1, 2, 3], [4, 1, 2, 3]];
    let s = 40.0;
    let vertices: Vec<NeuralVertex> = rest.iter().map(|&x| NeuralVertex::new(x, 1.0)).collect();
    let mut body = RealTimeNeuralSoftBody::new(
        vertices,
        elems.clone(),
        NeuralConstitutiveModel::linear_elastic(12, s),
        1000.0,
        0.0,
    );
    let mut rng = Lcg(31);
    let positions: Vec<[f64; 3]> = rest
        .iter()
        .map(|x| {
            [
                x[0] + 0.05 * rng.symmetric(),
                x[1] + 0.05 * rng.symmetric(),
                x[2] + 0.05 * rng.symmetric(),
            ]
        })
        .collect();
    for (v, p) in body.vertices.iter_mut().zip(positions.iter()) {
        v.position = *p;
    }
    let forces = body.elastic_forces();
    let eps = 1e-6;
    for (vertex, force) in forces.iter().enumerate() {
        for (axis, &component) in force.iter().enumerate() {
            let mut plus = positions.clone();
            let mut minus = positions.clone();
            plus[vertex][axis] += eps;
            minus[vertex][axis] -= eps;
            let gradient = (linear_energy(&rest, &plus, &elems, s)
                - linear_energy(&rest, &minus, &elems, s))
                / (2.0 * eps);
            assert!(
                (component + gradient).abs() < 1e-6 * (1.0 + gradient.abs()),
                "vertex {vertex} axis {axis}: force {component} vs -dE/dx {}",
                -gradient
            );
        }
    }
}

#[test]
fn a_free_body_conserves_momentum_and_a_pinned_body_relaxes() {
    // Free body, no gravity, no damping: the centre of mass cannot move.
    let rest = general_rest_tetrahedron();
    let mut rng = Lcg(5);
    let positions: Vec<[f64; 3]> = rest
        .iter()
        .map(|x| {
            [
                x[0] + 0.05 * rng.symmetric(),
                x[1] + 0.05 * rng.symmetric(),
                x[2] + 0.05 * rng.symmetric(),
            ]
        })
        .collect();
    let mut body = body_from(
        &rest,
        &positions,
        NeuralConstitutiveModel::linear_elastic(12, 50.0),
    );
    let com = |b: &RealTimeNeuralSoftBody| {
        let mut c = [0.0f64; 3];
        for v in &b.vertices {
            for (ck, pk) in c.iter_mut().zip(v.position.iter()) {
                *ck += pk / 4.0;
            }
        }
        c
    };
    let com0 = com(&body);
    for _ in 0..500 {
        body.step(1e-3, &[0.0; 3]);
    }
    let momentum: [f64; 3] = {
        let mut p = [0.0f64; 3];
        for v in &body.vertices {
            for (pk, vk) in p.iter_mut().zip(v.velocity.iter()) {
                *pk += v.mass * vk;
            }
        }
        p
    };
    assert_vec3_close(momentum, [0.0; 3], 1e-10, "momentum");
    assert_vec3_close(com(&body), com0, 1e-10, "centre of mass");
    let moved: f64 = body
        .vertices
        .iter()
        .zip(positions.iter())
        .map(|(v, p)| (v.position[0] - p[0]).abs() + (v.position[1] - p[1]).abs())
        .sum();
    assert!(moved > 1e-4, "the deformed body must actually move");

    // Pinned base, stretched apex, damped: the body resists the deformation and
    // returns to its rest shape.
    let mut vertices = tetrahedron_vertices();
    vertices[1] = NeuralVertex::new_static([1.0, 0.0, 0.0]);
    vertices[2] = NeuralVertex::new_static([0.0, 1.0, 0.0]);
    let mut pinned = RealTimeNeuralSoftBody::new(
        vertices,
        vec![[0, 1, 2, 3]],
        NeuralConstitutiveModel::linear_elastic(12, 100.0),
        1000.0,
        5.0,
    );
    pinned.vertices[0].position = [0.0, 0.0, 0.0];
    pinned.vertices[3].position = [0.0, 0.0, 1.05];
    for _ in 0..4000 {
        pinned.step(1e-3, &[0.0; 3]);
    }
    assert!(
        (pinned.vertices[3].position[2] - 1.0).abs() < 1e-3,
        "apex settled at z = {}",
        pinned.vertices[3].position[2]
    );
}

#[test]
fn rest_configuration_management_and_invalid_elements() {
    let vertices = tetrahedron_vertices();
    let mut body = RealTimeNeuralSoftBody::new(
        vertices,
        vec![[0, 1, 2, 3], [0, 1, 2, 7], [0, 1, 1, 3]],
        NeuralConstitutiveModel::new(4),
        1000.0,
        0.0,
    );
    // Out-of-range index and a repeated corner exert no force and never panic.
    assert_eq!(body.invalid_elements(), vec![1, 2]);
    body.step(0.01, &[0.0, -9.81, 0.0]);
    assert!(
        body.vertices
            .iter()
            .all(|v| v.position.iter().all(|c| c.is_finite()))
    );
    assert_eq!(body.rest_positions().len(), 4);
    assert!(body.set_rest_positions(vec![[0.0; 3]; 3]).is_err());
    // A new rest configuration is used for the strain.
    let scaled: Vec<[f64; 3]> = body
        .rest_positions()
        .iter()
        .map(|x| [2.0 * x[0], 2.0 * x[1], 2.0 * x[2]])
        .collect();
    body.set_rest_positions(scaled).unwrap();
    assert_eq!(body.rest_positions()[1], [2.0, 0.0, 0.0]);
    // Repairing the element list and refreshing recomputes the rest data.
    body.elements = vec![[0, 1, 2, 3]];
    body.refresh_rest_shapes();
    assert!(body.invalid_elements().is_empty());
}
