//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#[allow(unused_imports)]
use std::f64::consts::PI;

/// Dot product of two equal-length slices.
#[inline]
pub(super) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Euclidean norm of a slice.
#[inline]
pub(super) fn norm(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}
/// Element-wise addition.
#[inline]
pub(super) fn vadd(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
/// Element-wise subtraction.
#[inline]
pub(super) fn vsub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
/// Scale a vector.
#[inline]
pub(super) fn vscale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}
/// Rectified linear unit activation.
#[inline]
pub(super) fn relu(x: f64) -> f64 {
    x.max(0.0)
}
/// Hyperbolic tangent activation.
#[inline]
pub(super) fn tanh_act(x: f64) -> f64 {
    x.tanh()
}
/// Sigmoid activation.
#[inline]
pub(super) fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}
/// Softmax over a slice (returns normalised probabilities).
pub(super) fn softmax(v: &[f64]) -> Vec<f64> {
    let max = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = v.iter().map(|x| (x - max).exp()).collect();
    let sum: f64 = exps.iter().sum();
    exps.iter().map(|x| x / sum.max(1e-15)).collect()
}
/// Matrix–vector product: `W` is (out × in) row-major, returns length-`out` vec.
pub(super) fn matvec(w: &[Vec<f64>], v: &[f64]) -> Vec<f64> {
    w.iter().map(|row| dot(row, v)).collect()
}
/// Euclidean distance between two configurations.
pub(super) fn config_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y).powi(2))
        .sum::<f64>()
        .sqrt()
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::neural_control::AdaptiveControl;
    use crate::neural_control::Demonstration;
    use crate::neural_control::DenseLayer;
    use crate::neural_control::ImitationLearning;
    use crate::neural_control::ModelPredictiveControl;
    use crate::neural_control::MotionPlanning;
    use crate::neural_control::NeuralController;
    use crate::neural_control::NeuralDynamicsModel;
    use crate::neural_control::OUNoise;
    use crate::neural_control::PlannerNode;
    use crate::neural_control::QTable;
    use crate::neural_control::ReinforcementLearning;
    use crate::neural_control::ReplayBuffer;
    use crate::neural_control::RunningNormaliser;
    use crate::neural_control::Trajectory;
    use crate::neural_control::Transition;
    #[test]
    fn test_neural_controller_output_shape() {
        let nc = NeuralController::new(4, &[8, 8], 2, 1.0);
        let state = vec![0.1, -0.2, 0.3, 0.5];
        let action = nc.forward(&state);
        assert_eq!(action.len(), 2);
    }
    #[test]
    fn test_neural_controller_bounded_output() {
        let nc = NeuralController::new(3, &[16], 3, 5.0);
        let state = vec![1.0, -1.0, 0.5];
        let action = nc.forward(&state);
        for a in &action {
            assert!(
                *a > -5.0 - 1e-9 && *a < 5.0 + 1e-9,
                "action {a} out of bounds"
            );
        }
    }
    #[test]
    fn test_neural_controller_get_set_parameters() {
        let mut nc = NeuralController::new(2, &[4], 2, 1.0);
        let p = nc.get_parameters();
        assert!(!p.is_empty());
        let ones = vec![1.0_f64; p.len()];
        nc.set_parameters(&ones);
        let p2 = nc.get_parameters();
        for v in &p2 {
            assert!((*v - 1.0).abs() < 1e-10);
        }
    }
    #[test]
    fn test_neural_controller_sgd_changes_params() {
        let mut nc = NeuralController::new(2, &[4], 1, 1.0);
        let p_before = nc.get_parameters();
        let grad = vec![0.1_f64; p_before.len()];
        nc.sgd_step(&grad, 0.01);
        let p_after = nc.get_parameters();
        let changed = p_before
            .iter()
            .zip(p_after.iter())
            .any(|(a, b)| (a - b).abs() > 1e-12);
        assert!(changed, "SGD step must change parameters");
    }
    #[test]
    fn test_neural_controller_num_parameters() {
        let nc = NeuralController::new(4, &[8, 8], 2, 1.0);
        assert!(nc.num_parameters() > 0);
    }
    #[test]
    fn test_neural_controller_single_layer() {
        let nc = NeuralController::new(6, &[12], 3, 2.0);
        let state = vec![0.0; 6];
        let action = nc.forward(&state);
        assert_eq!(action.len(), 3);
    }
    #[test]
    fn test_neural_controller_three_hidden_layers() {
        let nc = NeuralController::new(8, &[16, 16, 8], 4, 1.0);
        let state = vec![0.5; 8];
        let action = nc.forward(&state);
        assert_eq!(action.len(), 4);
    }
    #[test]
    fn test_dense_layer_output_size() {
        let layer = DenseLayer::new(4, 8);
        let x = vec![1.0; 4];
        let y = layer.forward_relu(&x);
        assert_eq!(y.len(), 8);
    }
    #[test]
    fn test_dense_layer_relu_non_negative() {
        let layer = DenseLayer::new(4, 8);
        let x = vec![-1.0, -2.0, -3.0, -4.0];
        let y = layer.forward_relu(&x);
        for v in y {
            assert!(v >= 0.0);
        }
    }
    #[test]
    fn test_dense_layer_tanh_bounded() {
        let layer = DenseLayer::new(3, 5);
        let x = vec![1.0, -1.0, 0.0];
        let y = layer.forward_tanh(&x);
        for v in y {
            assert!(v > -1.0 - 1e-9 && v < 1.0 + 1e-9);
        }
    }
    #[test]
    fn test_replay_buffer_push_and_len() {
        let mut buf = ReplayBuffer::new(10);
        assert!(buf.is_empty());
        buf.push(Transition {
            state: vec![0.0],
            action: 0,
            reward: 1.0,
            next_state: vec![1.0],
            done: false,
        });
        assert_eq!(buf.len(), 1);
    }
    #[test]
    fn test_replay_buffer_wraps_around() {
        let mut buf = ReplayBuffer::new(3);
        for i in 0..6 {
            buf.push(Transition {
                state: vec![i as f64],
                action: i,
                reward: 0.0,
                next_state: vec![0.0],
                done: false,
            });
        }
        assert_eq!(buf.len(), 3);
    }
    #[test]
    fn test_replay_buffer_sample() {
        let mut buf = ReplayBuffer::new(100);
        for i in 0..50 {
            buf.push(Transition {
                state: vec![i as f64],
                action: 0,
                reward: 0.0,
                next_state: vec![0.0],
                done: false,
            });
        }
        let s = buf.sample(10);
        assert_eq!(s.len(), 10);
    }
    #[test]
    fn test_qtable_update_increases_value() {
        let mut qt = QTable::new(4, 2, 0.5, 0.9, 0.0);
        qt.update(0, 1, 10.0, 1, false);
        assert!(
            qt.q[0][1] > 0.0,
            "Q value should increase after positive reward"
        );
    }
    #[test]
    fn test_qtable_greedy_selects_best() {
        let mut qt = QTable::new(2, 3, 0.5, 0.9, 0.0);
        qt.q[0][2] = 100.0;
        assert_eq!(qt.select_action(0), 2);
    }
    #[test]
    fn test_qtable_epsilon_decay() {
        let mut qt = QTable::new(2, 2, 0.5, 0.9, 1.0);
        qt.decay_epsilon(0.9);
        assert!(qt.epsilon < 1.0);
        assert!(qt.epsilon >= 0.01);
    }
    #[test]
    fn test_rl_sample_action_valid_index() {
        let rl = ReinforcementLearning::new(4, &[8], 3, 0.99, 1e-3, 1e-3);
        let state = vec![0.0; 4];
        let (action, log_prob) = rl.sample_action(&state);
        assert!(action < 3);
        assert!(log_prob < 0.0, "log prob should be negative");
    }
    #[test]
    fn test_rl_compute_returns_discounted() {
        let mut rl = ReinforcementLearning::new(2, &[4], 2, 0.9, 1e-3, 1e-3);
        rl.record_step(vec![0.0; 2], 0, 0.0, 1.0);
        rl.record_step(vec![0.0; 2], 0, 0.0, 1.0);
        let ret = rl.compute_returns();
        assert_eq!(ret.len(), 2);
        assert!((ret[0] - 1.9).abs() < 1e-9, "ret[0]={}", ret[0]);
    }
    #[test]
    fn test_rl_reinforce_update_clears_trajectory() {
        let mut rl = ReinforcementLearning::new(2, &[4], 2, 0.99, 1e-4, 1e-4);
        for _ in 0..5 {
            rl.record_step(vec![0.0; 2], 0, -0.5, 1.0);
        }
        rl.reinforce_update();
        assert!(rl.trajectory.is_empty());
    }
    #[test]
    fn test_rl_a2c_update_runs() {
        let mut rl = ReinforcementLearning::new(4, &[8], 2, 0.99, 1e-4, 1e-4);
        let s = vec![0.1, -0.1, 0.2, -0.2];
        let ns = vec![0.2, -0.2, 0.3, -0.3];
        rl.a2c_update(&s, 0, 1.0, &ns, false);
    }
    #[test]
    fn test_dynamics_model_output_shape() {
        let m = NeuralDynamicsModel::new(4, 2, &[8]);
        let s = vec![0.0; 4];
        let a = vec![0.0; 2];
        let ns = m.predict(&s, &a);
        assert_eq!(ns.len(), 4);
    }
    #[test]
    fn test_dynamics_model_train_step_runs() {
        let mut m = NeuralDynamicsModel::new(3, 2, &[8]);
        let s = vec![1.0, 0.0, -1.0];
        let a = vec![0.5, -0.5];
        let ns = vec![0.9, 0.1, -0.9];
        m.train_step(&s, &a, &ns, 1e-4);
    }
    #[test]
    fn test_mpc_plan_returns_action() {
        let mpc = ModelPredictiveControl::new(2, 1, &[8], 5, 20, 0.9, 1.0);
        let state = vec![1.0, 0.0];
        let action = mpc.plan(&state, |s| s[0].powi(2) + s[1].powi(2));
        assert_eq!(action.len(), 1);
    }
    #[test]
    fn test_mpc_cem_plan_returns_action() {
        let mpc = ModelPredictiveControl::new(2, 1, &[8], 3, 20, 0.9, 1.0);
        let state = vec![0.5, -0.5];
        let action = mpc.cem_plan(&state, |s| s[0].powi(2), 3);
        assert_eq!(action.len(), 1);
    }
    #[test]
    fn test_adaptive_control_lyapunov_positive() {
        let ac = AdaptiveControl::new(2, 0.1);
        let state = vec![1.0, -1.0];
        let v = ac.lyapunov_value(&state);
        assert!(v >= 0.0, "Lyapunov value must be non-negative, got {v}");
    }
    #[test]
    fn test_adaptive_control_reference_model_step() {
        let mut ac = AdaptiveControl::new(2, 0.1);
        let r = vec![1.0, 1.0];
        ac.step_reference_model(&r, 0.01);
        let changed = ac.ref_state.iter().any(|&x| x.abs() > 1e-10);
        assert!(changed);
    }
    #[test]
    fn test_adaptive_control_parameter_update_runs() {
        let mut ac = AdaptiveControl::new(3, 0.1);
        let state = vec![1.0, -0.5, 0.2];
        let r = vec![0.0; 3];
        ac.update_parameters(&state, &r, 0.01);
    }
    #[test]
    fn test_adaptive_control_compute_action_shape() {
        let ac = AdaptiveControl::new(3, 0.1);
        let state = vec![1.0, -1.0, 0.5];
        let r = vec![0.1, 0.2, 0.3];
        let u = ac.compute_action(&state, &r);
        assert_eq!(u.len(), 3);
    }
    #[test]
    fn test_planner_random_config_in_bounds() {
        let planner = MotionPlanning::new(3, vec![-PI; 3], vec![PI; 3], 0.1, 100, 0.5, false);
        for _ in 0..10 {
            let q = planner.random_config();
            for (k, &qi) in q.iter().enumerate() {
                assert!(qi >= planner.lower[k] && qi <= planner.upper[k]);
            }
        }
    }
    #[test]
    fn test_planner_nearest_trivial() {
        let planner = MotionPlanning::new(2, vec![-1.0; 2], vec![1.0; 2], 0.1, 100, 0.1, false);
        let nodes = vec![
            PlannerNode::new(vec![0.0, 0.0], usize::MAX, 0.0),
            PlannerNode::new(vec![0.5, 0.5], usize::MAX, 0.0),
        ];
        let idx = planner.nearest(&nodes, &[0.6, 0.6]);
        assert_eq!(idx, 1);
    }
    #[test]
    fn test_planner_steer_no_overshoot() {
        let planner = MotionPlanning::new(2, vec![-1.0; 2], vec![1.0; 2], 0.2, 100, 0.1, false);
        let from = vec![0.0, 0.0];
        let to = vec![1.0, 0.0];
        let q_new = planner.steer(&from, &to);
        let d = config_dist(&q_new, &from);
        assert!(d <= planner.step_size + 1e-9);
    }
    #[test]
    fn test_planner_rrt_finds_trivial_path() {
        let planner = MotionPlanning::new(2, vec![0.0; 2], vec![1.0; 2], 0.2, 2000, 0.3, false);
        let start = vec![0.0, 0.0];
        let goal = vec![0.5, 0.5];
        let path = planner.rrt(&start, &goal);
        let _ = path;
    }
    #[test]
    fn test_planner_rrt_star_runs() {
        let planner = MotionPlanning::new(2, vec![0.0; 2], vec![1.0; 2], 0.15, 500, 0.2, true);
        let _ = planner.rrt(&[0.0, 0.0], &[0.8, 0.8]);
    }
    #[test]
    fn test_planner_prm_builds() {
        let planner = MotionPlanning::new(2, vec![0.0; 2], vec![1.0; 2], 0.2, 100, 0.1, false);
        let adj = planner.build_prm(20, 3);
        assert_eq!(adj.len(), 20);
    }
    #[test]
    fn test_planner_trajectory_optimise() {
        let planner = MotionPlanning::new(2, vec![0.0; 2], vec![10.0; 2], 0.5, 100, 0.2, false);
        let path = vec![vec![0.0, 0.0], vec![5.0, 5.0], vec![10.0, 10.0]];
        let opt = planner.optimise_trajectory(path.clone(), 10, 0.01);
        assert_eq!(opt.len(), path.len());
    }
    #[test]
    fn test_imitation_add_demonstration() {
        let mut il = ImitationLearning::new(4, 2, &[8], 1e-3);
        il.add_demonstration(Demonstration::new(vec![0.0; 4], vec![1.0; 2]));
        assert_eq!(il.dataset.len(), 1);
    }
    #[test]
    fn test_imitation_bc_step_runs() {
        let mut il = ImitationLearning::new(4, 2, &[8], 1e-3);
        for _ in 0..10 {
            il.add_demonstration(Demonstration::new(vec![0.0; 4], vec![0.0; 2]));
        }
        il.bc_step(4);
    }
    #[test]
    fn test_dagger_beta_decreases() {
        let mut il = ImitationLearning::new(2, 1, &[4], 1e-3);
        let b0 = il.dagger_beta();
        il.iter = 5;
        let b5 = il.dagger_beta();
        assert!(b5 < b0);
    }
    #[test]
    fn test_dagger_action_follows_expert_initially() {
        let mut il = ImitationLearning::new(2, 1, &[4], 1e-3);
        il.iter = 0;
        let state = vec![0.0; 2];
        let expert = vec![42.0];
        let action = il.dagger_action(&state, &expert);
        assert_eq!(action, expert);
    }
    #[test]
    fn test_irl_reward_dot_product() {
        let mut il = ImitationLearning::new(3, 1, &[4], 1e-3);
        il.reward_weights = vec![1.0, 0.0, 0.0];
        let r = il.irl_reward(&[5.0, 3.0, 1.0]);
        assert!((r - 5.0).abs() < 1e-9, "r={r}");
    }
    #[test]
    fn test_irl_update_normalises_weights() {
        let mut il = ImitationLearning::new(3, 1, &[4], 1e-3);
        let expert_fe = vec![1.0, 0.0, 0.0];
        let learner_fe = vec![0.0, 0.0, 0.0];
        il.irl_update(&learner_fe, &expert_fe, 0.1);
        let n = norm(&il.reward_weights);
        assert!(
            (n - 1.0).abs() < 1e-6,
            "weights should be unit norm, got {n}"
        );
    }
    #[test]
    fn test_ou_noise_dimension() {
        let mut ou = OUNoise::new(4, 0.15, 0.2);
        let s = ou.sample();
        assert_eq!(s.len(), 4);
    }
    #[test]
    fn test_ou_noise_reset() {
        let mut ou = OUNoise::new(3, 0.15, 0.2);
        for _ in 0..10 {
            ou.sample();
        }
        ou.reset();
        assert!(ou.state.iter().all(|&x| x == 0.0));
    }
    #[test]
    fn test_running_normaliser_mean() {
        let mut rn = RunningNormaliser::new(2);
        rn.update(&[1.0, 2.0]);
        rn.update(&[3.0, 4.0]);
        assert!((rn.mean[0] - 2.0).abs() < 1e-9);
        assert!((rn.mean[1] - 3.0).abs() < 1e-9);
    }
    #[test]
    fn test_running_normaliser_normalise() {
        let mut rn = RunningNormaliser::new(1);
        for i in 0..100 {
            rn.update(&[i as f64]);
        }
        let z = rn.normalise(&[rn.mean[0]]);
        assert!(
            z[0].abs() < 1e-6,
            "normalised mean should be ~0, got {}",
            z[0]
        );
    }
    #[test]
    fn test_trajectory_interpolate_endpoints() {
        let path = vec![vec![0.0, 0.0], vec![1.0, 1.0]];
        let traj = Trajectory::from_path(path, 1.0);
        let start = traj.interpolate(0.0);
        let end = traj.interpolate(1.0);
        assert!((start[0] - 0.0).abs() < 1e-9);
        assert!((end[0] - 1.0).abs() < 1e-9);
    }
    #[test]
    fn test_trajectory_interpolate_midpoint() {
        let path = vec![vec![0.0], vec![2.0]];
        let traj = Trajectory::from_path(path, 2.0);
        let mid = traj.interpolate(1.0);
        assert!((mid[0] - 1.0).abs() < 1e-6, "mid={}", mid[0]);
    }
    #[test]
    fn test_trajectory_duration() {
        let path = vec![vec![0.0], vec![1.0], vec![2.0]];
        let traj = Trajectory::from_path(path, 5.0);
        assert!((traj.duration() - 5.0).abs() < 1e-9);
    }
    #[test]
    fn test_softmax_sums_to_one() {
        let v = vec![1.0, 2.0, 3.0];
        let p = softmax(&v);
        let s: f64 = p.iter().sum();
        assert!((s - 1.0).abs() < 1e-9, "softmax sum={s}");
    }
    #[test]
    fn test_softmax_largest_wins() {
        let v = vec![0.0, 0.0, 10.0];
        let p = softmax(&v);
        assert!(p[2] > p[0] && p[2] > p[1]);
    }
}
