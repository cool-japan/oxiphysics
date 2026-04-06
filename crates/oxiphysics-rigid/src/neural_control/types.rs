//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop)]
#[allow(unused_imports)]
use super::functions::*;
use rand::RngExt;
/// Online running-mean / running-variance normaliser (Welford's algorithm).
#[derive(Debug, Clone)]
pub struct RunningNormaliser {
    /// Number of samples seen.
    pub count: usize,
    /// Running mean.
    pub mean: Vec<f64>,
    /// Running M2 (for variance).
    pub m2: Vec<f64>,
    /// Dimension.
    pub dim: usize,
}
impl RunningNormaliser {
    /// Create a normaliser for `dim`-dimensional inputs.
    pub fn new(dim: usize) -> Self {
        Self {
            count: 0,
            mean: vec![0.0; dim],
            m2: vec![0.0; dim],
            dim,
        }
    }
    /// Update statistics with a new sample.
    pub fn update(&mut self, x: &[f64]) {
        self.count += 1;
        for i in 0..self.dim {
            let delta = x[i] - self.mean[i];
            self.mean[i] += delta / self.count as f64;
            let delta2 = x[i] - self.mean[i];
            self.m2[i] += delta * delta2;
        }
    }
    /// Normalise a sample: (x - mean) / std.
    pub fn normalise(&self, x: &[f64]) -> Vec<f64> {
        x.iter()
            .enumerate()
            .map(|(i, &xi)| {
                let variance = if self.count > 1 {
                    self.m2[i] / (self.count - 1) as f64
                } else {
                    1.0
                };
                (xi - self.mean[i]) / variance.sqrt().max(1e-8)
            })
            .collect()
    }
}
/// Feedforward neural network controller that maps a state vector to an action vector.
///
/// Architecture: `state_dim → [hidden_1, hidden_2?] → action_dim`.
/// Output is passed through `tanh` and scaled by `action_scale` so that
/// actions stay bounded.
#[derive(Debug, Clone)]
pub struct NeuralController {
    /// Hidden layers (weights + biases).
    pub layers: Vec<DenseLayer>,
    /// Output layer.
    pub output_layer: DenseLayer,
    /// Hidden layer activation.
    pub activation: Activation,
    /// Multiplier applied to the tanh output.
    pub action_scale: f64,
    /// Dimensionality of the input state.
    pub state_dim: usize,
    /// Dimensionality of the action output.
    pub action_dim: usize,
}
impl NeuralController {
    /// Build a controller with one or two hidden layers.
    ///
    /// * `state_dim`    — length of the state vector
    /// * `hidden_sizes` — sizes of the hidden layers (1–2 entries)
    /// * `action_dim`   — length of the action vector
    /// * `action_scale` — output multiplier (e.g. max torque in N·m)
    pub fn new(
        state_dim: usize,
        hidden_sizes: &[usize],
        action_dim: usize,
        action_scale: f64,
    ) -> Self {
        let activation = Activation::Tanh;
        let mut layers = Vec::new();
        let mut prev = state_dim;
        for &h in hidden_sizes {
            layers.push(DenseLayer::new(prev, h));
            prev = h;
        }
        let output_layer = DenseLayer::new(prev, action_dim);
        Self {
            layers,
            output_layer,
            activation,
            action_scale,
            state_dim,
            action_dim,
        }
    }
    /// Run the network forward.  Returns bounded actions ∈ (-action_scale, +action_scale).
    pub fn forward(&self, state: &[f64]) -> Vec<f64> {
        let mut x = state.to_vec();
        for layer in &self.layers {
            x = match self.activation {
                Activation::Relu => layer.forward_relu(&x),
                Activation::Tanh => layer.forward_tanh(&x),
                Activation::Sigmoid => {
                    let pre = layer.forward_linear(&x);
                    pre.into_iter().map(sigmoid).collect()
                }
                Activation::Linear => layer.forward_linear(&x),
            };
        }
        let raw = self.output_layer.forward_linear(&x);
        raw.into_iter()
            .map(|v| tanh_act(v) * self.action_scale)
            .collect()
    }
    /// Collect all parameters (weights + biases) as a flat vector.
    pub fn get_parameters(&self) -> Vec<f64> {
        let mut params = Vec::new();
        for layer in &self.layers {
            for row in &layer.weights {
                params.extend_from_slice(row);
            }
            params.extend_from_slice(&layer.biases);
        }
        for row in &self.output_layer.weights {
            params.extend_from_slice(row);
        }
        params.extend_from_slice(&self.output_layer.biases);
        params
    }
    /// Set all parameters from a flat vector (same order as `get_parameters`).
    pub fn set_parameters(&mut self, params: &[f64]) {
        let mut idx = 0;
        for layer in &mut self.layers {
            for row in &mut layer.weights {
                for w in row.iter_mut() {
                    *w = params[idx];
                    idx += 1;
                }
            }
            for b in &mut layer.biases {
                *b = params[idx];
                idx += 1;
            }
        }
        for row in &mut self.output_layer.weights {
            for w in row.iter_mut() {
                *w = params[idx];
                idx += 1;
            }
        }
        for b in &mut self.output_layer.biases {
            *b = params[idx];
            idx += 1;
        }
    }
    /// Count total number of parameters.
    pub fn num_parameters(&self) -> usize {
        self.get_parameters().len()
    }
    /// Apply a gradient step (SGD) given a flat gradient vector and learning rate.
    pub fn sgd_step(&mut self, grad: &[f64], lr: f64) {
        let mut params = self.get_parameters();
        for (p, g) in params.iter_mut().zip(grad.iter()) {
            *p -= lr * g;
        }
        self.set_parameters(&params);
    }
}
/// Policy gradient (REINFORCE) and actor-critic reinforcement learning.
#[derive(Debug, Clone)]
pub struct ReinforcementLearning {
    /// Actor network (policy).
    pub actor: NeuralController,
    /// Critic network (value function).
    pub critic: NeuralController,
    /// Discount factor γ.
    pub gamma: f64,
    /// Actor learning rate.
    pub actor_lr: f64,
    /// Critic learning rate.
    pub critic_lr: f64,
    /// Episode trajectory storage: (state, action_index, log_prob, reward).
    pub trajectory: Vec<(Vec<f64>, usize, f64, f64)>,
}
impl ReinforcementLearning {
    /// Construct an RL agent with actor and critic networks.
    pub fn new(
        state_dim: usize,
        hidden: &[usize],
        action_dim: usize,
        gamma: f64,
        actor_lr: f64,
        critic_lr: f64,
    ) -> Self {
        let actor = NeuralController::new(state_dim, hidden, action_dim, 1.0);
        let critic = NeuralController::new(state_dim, hidden, 1, 1.0);
        Self {
            actor,
            critic,
            gamma,
            actor_lr,
            critic_lr,
            trajectory: Vec::new(),
        }
    }
    /// Sample an action from the actor's softmax policy.
    /// Returns (action index, log probability).
    pub fn sample_action(&self, state: &[f64]) -> (usize, f64) {
        use rand::RngExt;
        let logits = self.actor.forward(state);
        let probs = softmax(&logits);
        let mut rng = rand::rng();
        let sample: f64 = rng.random_range(0.0..1.0);
        let mut cumsum = 0.0;
        let mut action = probs.len() - 1;
        for (i, &p) in probs.iter().enumerate() {
            cumsum += p;
            if sample < cumsum {
                action = i;
                break;
            }
        }
        let log_prob = (probs[action] + 1e-8).ln();
        (action, log_prob)
    }
    /// Store a step in the current episode trajectory.
    pub fn record_step(&mut self, state: Vec<f64>, action: usize, log_prob: f64, reward: f64) {
        self.trajectory.push((state, action, log_prob, reward));
    }
    /// Compute discounted returns for the stored trajectory.
    pub fn compute_returns(&self) -> Vec<f64> {
        let n = self.trajectory.len();
        let mut returns = vec![0.0; n];
        let mut running = 0.0;
        for i in (0..n).rev() {
            running = self.trajectory[i].3 + self.gamma * running;
            returns[i] = running;
        }
        returns
    }
    /// REINFORCE policy gradient update step.
    /// Clears the trajectory after updating.
    pub fn reinforce_update(&mut self) {
        let returns = self.compute_returns();
        let mean_ret: f64 = returns.iter().sum::<f64>() / returns.len().max(1) as f64;
        let std_ret: f64 = {
            let var: f64 = returns.iter().map(|r| (r - mean_ret).powi(2)).sum::<f64>()
                / returns.len().max(1) as f64;
            var.sqrt().max(1e-8)
        };
        let n_params = self.actor.num_parameters();
        let mut grad = vec![0.0_f64; n_params];
        for (i, (state, action_idx, _log_prob, _reward)) in self.trajectory.iter().enumerate() {
            let logits = self.actor.forward(state);
            let probs = softmax(&logits);
            let advantage = (returns[i] - mean_ret) / std_ret;
            let eps = 1e-5;
            let params = self.actor.get_parameters();
            let mut perturbed = self.actor.clone();
            for j in 0..n_params {
                let mut p_plus = params.clone();
                p_plus[j] += eps;
                perturbed.set_parameters(&p_plus);
                let lp_plus = (softmax(&perturbed.forward(state))[*action_idx] + 1e-8).ln();
                let mut p_minus = params.clone();
                p_minus[j] -= eps;
                perturbed.set_parameters(&p_minus);
                let lp_minus = (softmax(&perturbed.forward(state))[*action_idx] + 1e-8).ln();
                perturbed.set_parameters(&params);
                let dlp_dw = (lp_plus - lp_minus) / (2.0 * eps);
                let _ = probs[*action_idx];
                grad[j] -= advantage * dlp_dw;
            }
        }
        let scale = 1.0 / self.trajectory.len().max(1) as f64;
        let grad_scaled: Vec<f64> = grad.iter().map(|g| g * scale).collect();
        self.actor.sgd_step(&grad_scaled, self.actor_lr);
        self.trajectory.clear();
    }
    /// Advantage actor-critic (A2C) update for a single transition.
    pub fn a2c_update(
        &mut self,
        state: &[f64],
        action: usize,
        reward: f64,
        next_state: &[f64],
        done: bool,
    ) {
        let value = self.critic.forward(state)[0];
        let next_value = if done {
            0.0
        } else {
            self.critic.forward(next_state)[0]
        };
        let td_target = reward + self.gamma * next_value;
        let advantage = td_target - value;
        let n_critic = self.critic.num_parameters();
        let mut c_grad = vec![0.0_f64; n_critic];
        let eps = 1e-5;
        let c_params = self.critic.get_parameters();
        let mut c_clone = self.critic.clone();
        for j in 0..n_critic {
            let mut p = c_params.clone();
            p[j] += eps;
            c_clone.set_parameters(&p);
            let v_plus = c_clone.forward(state)[0];
            p[j] -= 2.0 * eps;
            c_clone.set_parameters(&p);
            let v_minus = c_clone.forward(state)[0];
            c_clone.set_parameters(&c_params);
            let dv_dw = (v_plus - v_minus) / (2.0 * eps);
            c_grad[j] = -2.0 * advantage * dv_dw;
        }
        self.critic.sgd_step(&c_grad, self.critic_lr);
        let n_actor = self.actor.num_parameters();
        let mut a_grad = vec![0.0_f64; n_actor];
        let a_params = self.actor.get_parameters();
        let mut a_clone = self.actor.clone();
        for j in 0..n_actor {
            let mut p = a_params.clone();
            p[j] += eps;
            a_clone.set_parameters(&p);
            let lp_plus = (softmax(&a_clone.forward(state))[action] + 1e-8).ln();
            p[j] -= 2.0 * eps;
            a_clone.set_parameters(&p);
            let lp_minus = (softmax(&a_clone.forward(state))[action] + 1e-8).ln();
            a_clone.set_parameters(&a_params);
            let dlp = (lp_plus - lp_minus) / (2.0 * eps);
            a_grad[j] = -advantage * dlp;
        }
        self.actor.sgd_step(&a_grad, self.actor_lr);
    }
}
/// A single (s, a, r, s') transition.
#[derive(Debug, Clone)]
pub struct Transition {
    /// State before action.
    pub state: Vec<f64>,
    /// Action taken.
    pub action: usize,
    /// Reward received.
    pub reward: f64,
    /// Next state.
    pub next_state: Vec<f64>,
    /// Terminal flag.
    pub done: bool,
}
/// A single fully-connected layer with bias and activation.
#[derive(Debug, Clone)]
pub struct DenseLayer {
    /// Weight matrix (out_size × in_size).
    pub weights: Vec<Vec<f64>>,
    /// Bias vector (out_size).
    pub biases: Vec<f64>,
}
impl DenseLayer {
    /// Construct a layer with Xavier-uniform initialisation.
    pub fn new(in_size: usize, out_size: usize) -> Self {
        let limit = (6.0 / (in_size + out_size) as f64).sqrt();
        let mut rng = rand::rng();
        let weights = (0..out_size)
            .map(|_| {
                (0..in_size)
                    .map(|_| rng.random_range(-limit..limit))
                    .collect()
            })
            .collect();
        let biases = vec![0.0; out_size];
        Self { weights, biases }
    }
    /// Forward pass with ReLU activation.
    pub fn forward_relu(&self, input: &[f64]) -> Vec<f64> {
        let pre: Vec<f64> = matvec(&self.weights, input)
            .into_iter()
            .zip(self.biases.iter())
            .map(|(z, b)| relu(z + b))
            .collect();
        pre
    }
    /// Forward pass with tanh activation.
    pub fn forward_tanh(&self, input: &[f64]) -> Vec<f64> {
        matvec(&self.weights, input)
            .into_iter()
            .zip(self.biases.iter())
            .map(|(z, b)| tanh_act(z + b))
            .collect()
    }
    /// Forward pass with linear (identity) activation.
    pub fn forward_linear(&self, input: &[f64]) -> Vec<f64> {
        matvec(&self.weights, input)
            .into_iter()
            .zip(self.biases.iter())
            .map(|(z, b)| z + b)
            .collect()
    }
    /// Number of output neurons.
    pub fn out_size(&self) -> usize {
        self.biases.len()
    }
    /// Number of input neurons.
    pub fn in_size(&self) -> usize {
        if self.weights.is_empty() {
            0
        } else {
            self.weights[0].len()
        }
    }
}
/// Behavioural cloning, DAgger, and inverse reinforcement learning.
#[derive(Debug, Clone)]
pub struct ImitationLearning {
    /// Policy network being trained.
    pub policy: NeuralController,
    /// Dataset of demonstrations.
    pub dataset: Vec<Demonstration>,
    /// Learning rate for supervised updates.
    pub lr: f64,
    /// Iteration counter (used for DAgger mixing schedule).
    pub iter: usize,
    /// Reward weight vector for IRL (linear reward model: w · φ(s)).
    pub reward_weights: Vec<f64>,
    /// Feature dimension for IRL.
    pub feature_dim: usize,
}
impl ImitationLearning {
    /// Create an imitation learning agent.
    pub fn new(state_dim: usize, action_dim: usize, hidden: &[usize], lr: f64) -> Self {
        let policy = NeuralController::new(state_dim, hidden, action_dim, 1.0);
        let feature_dim = state_dim;
        let reward_weights = vec![1.0 / state_dim as f64; feature_dim];
        Self {
            policy,
            dataset: Vec::new(),
            lr,
            iter: 0,
            reward_weights,
            feature_dim,
        }
    }
    /// Add a demonstration to the dataset.
    pub fn add_demonstration(&mut self, demo: Demonstration) {
        self.dataset.push(demo);
    }
    /// Behavioural cloning: one mini-batch supervised step (MSE loss on actions).
    pub fn bc_step(&mut self, batch_size: usize) {
        if self.dataset.is_empty() {
            return;
        }
        let mut rng = rand::rng();
        let n = self.dataset.len();
        let n_params = self.policy.num_parameters();
        let mut grad = vec![0.0_f64; n_params];
        let eps = 1e-5;
        let params = self.policy.get_parameters();
        let mut clone = self.policy.clone();
        for _ in 0..batch_size {
            let demo = &self.dataset[rng.random_range(0..n)];
            let pred = self.policy.forward(&demo.state);
            let residual: Vec<f64> = vsub(&pred, &demo.action);
            for j in 0..n_params {
                let mut p = params.clone();
                p[j] += eps;
                clone.set_parameters(&p);
                let pp = clone.forward(&demo.state);
                p[j] -= 2.0 * eps;
                clone.set_parameters(&p);
                let pm = clone.forward(&demo.state);
                clone.set_parameters(&params);
                let dloss: f64 = residual
                    .iter()
                    .zip(pp.iter().zip(pm.iter()))
                    .map(|(r, (fp, fm))| r * (fp - fm) / (2.0 * eps))
                    .sum();
                grad[j] += dloss;
            }
        }
        let scale = 1.0 / batch_size.max(1) as f64;
        let scaled_grad: Vec<f64> = grad.iter().map(|g| g * scale).collect();
        self.policy.sgd_step(&scaled_grad, self.lr);
        self.iter += 1;
    }
    /// DAgger mixing probability: probability of following expert vs policy.
    /// Uses schedule β_i = p^i where p = 0.5.
    pub fn dagger_beta(&self) -> f64 {
        0.5_f64.powi(self.iter as i32)
    }
    /// Select action for DAgger: mix between expert and policy.
    pub fn dagger_action(&self, state: &[f64], expert_action: &[f64]) -> Vec<f64> {
        let mut rng = rand::rng();
        if rng.random_range(0.0..1.0) < self.dagger_beta() {
            expert_action.to_vec()
        } else {
            self.policy.forward(state)
        }
    }
    /// Aggregate a new demonstration from DAgger (expert labels policy trajectory).
    pub fn dagger_aggregate(&mut self, state: Vec<f64>, expert_action: Vec<f64>) {
        self.add_demonstration(Demonstration::new(state, expert_action));
    }
    /// IRL feature expectation under current policy (Monte Carlo estimate).
    pub fn feature_expectation(&self, trajectories: &[Vec<Vec<f64>>], gamma: f64) -> Vec<f64> {
        let mut fe = vec![0.0_f64; self.feature_dim];
        let n_traj = trajectories.len().max(1);
        for traj in trajectories {
            let mut disc = 1.0;
            for state in traj {
                for (k, &s) in state.iter().enumerate().take(self.feature_dim) {
                    fe[k] += disc * s;
                }
                disc *= gamma;
            }
        }
        fe.iter().map(|x| x / n_traj as f64).collect()
    }
    /// IRL expert feature expectation from the demonstration dataset.
    pub fn expert_feature_expectation(&self, gamma: f64) -> Vec<f64> {
        if self.dataset.is_empty() {
            return vec![0.0; self.feature_dim];
        }
        let traj: Vec<Vec<Vec<f64>>> = self.dataset.iter().map(|d| vec![d.state.clone()]).collect();
        self.feature_expectation(&traj, gamma)
    }
    /// IRL reward weight update: project-and-normalise gradient descent.
    pub fn irl_update(&mut self, learner_fe: &[f64], expert_fe: &[f64], irl_lr: f64) {
        let diff: Vec<f64> = vsub(expert_fe, learner_fe);
        for (w, d) in self.reward_weights.iter_mut().zip(diff.iter()) {
            *w += irl_lr * d;
        }
        let n = norm(&self.reward_weights).max(1e-12);
        for w in &mut self.reward_weights {
            *w /= n;
        }
    }
    /// Compute the IRL reward for a state using the current weight vector.
    pub fn irl_reward(&self, state: &[f64]) -> f64 {
        state
            .iter()
            .zip(self.reward_weights.iter())
            .map(|(s, w)| s * w)
            .sum()
    }
}
/// A time-parameterised trajectory waypoint.
#[derive(Debug, Clone)]
pub struct TrajectoryWaypoint {
    /// Time stamp.
    pub time: f64,
    /// Configuration at this waypoint.
    pub config: Vec<f64>,
    /// Optional velocity at waypoint.
    pub velocity: Option<Vec<f64>>,
}
impl TrajectoryWaypoint {
    /// Create a waypoint.
    pub fn new(time: f64, config: Vec<f64>) -> Self {
        Self {
            time,
            config,
            velocity: None,
        }
    }
}
/// A piecewise-linear trajectory for rigid body motion planning.
#[derive(Debug, Clone)]
pub struct Trajectory {
    /// Ordered waypoints.
    pub waypoints: Vec<TrajectoryWaypoint>,
}
impl Trajectory {
    /// Build from a path (equal time spacing).
    pub fn from_path(path: Vec<Vec<f64>>, total_time: f64) -> Self {
        let n = path.len();
        let dt = if n > 1 {
            total_time / (n - 1) as f64
        } else {
            0.0
        };
        let waypoints = path
            .into_iter()
            .enumerate()
            .map(|(i, c)| TrajectoryWaypoint::new(i as f64 * dt, c))
            .collect();
        Self { waypoints }
    }
    /// Linearly interpolate the configuration at time `t`.
    pub fn interpolate(&self, t: f64) -> Vec<f64> {
        if self.waypoints.is_empty() {
            return Vec::new();
        }
        if t <= self.waypoints[0].time {
            return self.waypoints[0].config.clone();
        }
        let last = self
            .waypoints
            .last()
            .expect("collection should not be empty");
        if t >= last.time {
            return last.config.clone();
        }
        for i in 0..(self.waypoints.len() - 1) {
            let t0 = self.waypoints[i].time;
            let t1 = self.waypoints[i + 1].time;
            if t >= t0 && t <= t1 {
                let alpha = (t - t0) / (t1 - t0).max(1e-15);
                return self.waypoints[i]
                    .config
                    .iter()
                    .zip(self.waypoints[i + 1].config.iter())
                    .map(|(a, b)| a + alpha * (b - a))
                    .collect();
            }
        }
        last.config.clone()
    }
    /// Total duration of the trajectory.
    pub fn duration(&self) -> f64 {
        self.waypoints.last().map(|w| w.time).unwrap_or(0.0)
    }
}
/// Neural dynamics model: predicts the next state given current state + action.
#[derive(Debug, Clone)]
pub struct NeuralDynamicsModel {
    /// Internal network.
    pub network: NeuralController,
    /// Dimension of state.
    pub state_dim: usize,
    /// Dimension of action.
    pub action_dim: usize,
}
impl NeuralDynamicsModel {
    /// Construct a neural dynamics model.
    pub fn new(state_dim: usize, action_dim: usize, hidden: &[usize]) -> Self {
        let network = NeuralController::new(state_dim + action_dim, hidden, state_dim, 1.0);
        Self {
            network,
            state_dim,
            action_dim,
        }
    }
    /// Predict next state: concatenate state and action, run forward pass.
    pub fn predict(&self, state: &[f64], action: &[f64]) -> Vec<f64> {
        let mut input = state.to_vec();
        input.extend_from_slice(action);
        let delta = self.network.forward(&input);
        vadd(state, &delta)
    }
    /// Train on a single (s, a, s') sample using MSE loss gradient.
    pub fn train_step(&mut self, state: &[f64], action: &[f64], next_state: &[f64], lr: f64) {
        let predicted = self.predict(state, action);
        let residual: Vec<f64> = vsub(&predicted, next_state);
        let n = self.network.num_parameters();
        let params = self.network.get_parameters();
        let mut grad = vec![0.0_f64; n];
        let eps = 1e-5;
        let mut clone = self.network.clone();
        let mut inp = state.to_vec();
        inp.extend_from_slice(action);
        for j in 0..n {
            let mut p = params.clone();
            p[j] += eps;
            clone.set_parameters(&p);
            let fwd_plus = {
                let raw = clone.forward(&inp);
                vadd(state, &raw)
            };
            p[j] -= 2.0 * eps;
            clone.set_parameters(&p);
            let fwd_minus = {
                let raw = clone.forward(&inp);
                vadd(state, &raw)
            };
            clone.set_parameters(&params);
            let dloss: f64 = residual
                .iter()
                .zip(fwd_plus.iter().zip(fwd_minus.iter()))
                .map(|(r, (fp, fm))| r * (fp - fm) / (2.0 * eps))
                .sum();
            grad[j] = dloss;
        }
        self.network.sgd_step(&grad, lr);
    }
}
/// Activation function used in hidden layers.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Activation {
    /// Rectified linear unit.
    Relu,
    /// Hyperbolic tangent.
    Tanh,
    /// Sigmoid (0–1).
    Sigmoid,
    /// Linear / identity.
    Linear,
}
/// MRAC (Model Reference Adaptive Control) for rigid body systems.
///
/// Tracks a reference model  ẋ_m = A_m x_m + B_m r  using an adaptive law
/// derived from Lyapunov stability analysis (MIT rule).
#[derive(Debug, Clone)]
pub struct AdaptiveControl {
    /// Reference model state matrix (diagonal approximation as vector).
    pub am_diag: Vec<f64>,
    /// Reference model input gain (diagonal).
    pub bm_diag: Vec<f64>,
    /// Adaptive parameter matrix (theta) — maps \[state; reference_input\] → action.
    pub theta: Vec<Vec<f64>>,
    /// Adaptation rate (γ).
    pub gamma_adapt: f64,
    /// State dimension.
    pub state_dim: usize,
    /// Reference state (x_m).
    pub ref_state: Vec<f64>,
    /// Estimated plant parameter vector for each output channel.
    pub param_estimate: Vec<f64>,
    /// Lyapunov matrix P (diagonal, symmetric positive definite).
    pub lyapunov_p: Vec<f64>,
}
impl AdaptiveControl {
    /// Initialise MRAC controller.
    pub fn new(state_dim: usize, gamma_adapt: f64) -> Self {
        let am_diag = vec![-1.0; state_dim];
        let bm_diag = vec![1.0; state_dim];
        let theta = vec![vec![0.0; state_dim + 1]; state_dim];
        let lyapunov_p = vec![1.0; state_dim];
        let param_estimate = vec![1.0; state_dim];
        Self {
            am_diag,
            bm_diag,
            theta,
            gamma_adapt,
            state_dim,
            ref_state: vec![0.0; state_dim],
            param_estimate,
            lyapunov_p,
        }
    }
    /// Integrate the reference model one step.
    pub fn step_reference_model(&mut self, r: &[f64], dt: f64) {
        for i in 0..self.state_dim {
            let x_dot = self.am_diag[i] * self.ref_state[i] + self.bm_diag[i] * r[i];
            self.ref_state[i] += x_dot * dt;
        }
    }
    /// Compute the adaptive control action for a given plant state.
    pub fn compute_action(&self, plant_state: &[f64], r: &[f64]) -> Vec<f64> {
        (0..self.state_dim)
            .map(|i| {
                let x_part: f64 = plant_state
                    .iter()
                    .enumerate()
                    .map(|(j, &x)| self.theta[i][j] * x)
                    .sum();
                let r_part = self.theta[i][self.state_dim] * r[i];
                x_part + r_part
            })
            .collect()
    }
    /// Update adaptive parameters using the MIT (gradient) rule:
    ///   dθ/dt = -γ · e · ∂y/∂θ
    pub fn update_parameters(&mut self, plant_state: &[f64], r: &[f64], dt: f64) {
        let e: Vec<f64> = plant_state
            .iter()
            .zip(self.ref_state.iter())
            .map(|(x, m)| x - m)
            .collect();
        for i in 0..self.state_dim {
            for j in 0..self.state_dim {
                let dtheta = -self.gamma_adapt * self.lyapunov_p[i] * e[i] * plant_state[j];
                self.theta[i][j] += dtheta * dt;
            }
            let dtheta_r = -self.gamma_adapt * self.lyapunov_p[i] * e[i] * r[i];
            self.theta[i][self.state_dim] += dtheta_r * dt;
        }
    }
    /// Lyapunov stability certificate: V = e^T P e (should be non-negative and decreasing).
    pub fn lyapunov_value(&self, plant_state: &[f64]) -> f64 {
        plant_state
            .iter()
            .zip(self.ref_state.iter())
            .zip(self.lyapunov_p.iter())
            .map(|((x, m), p)| p * (x - m).powi(2))
            .sum()
    }
    /// Estimate an unknown plant parameter via gradient descent on the tracking error.
    pub fn estimate_parameter(&mut self, error: f64, regressor: f64, dt: f64) {
        for p in &mut self.param_estimate {
            *p -= self.gamma_adapt * error * regressor * dt;
        }
    }
    /// Compute the Lyapunov derivative dV/dt ≈ (V(t) - V(t-dt)) / dt.
    pub fn lyapunov_derivative(&self, plant_state: &[f64], prev_v: f64, dt: f64) -> f64 {
        (self.lyapunov_value(plant_state) - prev_v) / dt
    }
}
/// A single demonstration: state → expert action (continuous or discrete).
#[derive(Debug, Clone)]
pub struct Demonstration {
    /// Observed state.
    pub state: Vec<f64>,
    /// Expert action (continuous).
    pub action: Vec<f64>,
}
impl Demonstration {
    /// Create a demonstration.
    pub fn new(state: Vec<f64>, action: Vec<f64>) -> Self {
        Self { state, action }
    }
}
/// A node in a motion planning tree / roadmap.
#[derive(Debug, Clone)]
pub struct PlannerNode {
    /// Configuration-space position (joint angles or Cartesian).
    pub config: Vec<f64>,
    /// Parent node index in the node list (usize::MAX = root).
    pub parent: usize,
    /// Cost-to-reach from root (used by RRT*).
    pub cost: f64,
}
impl PlannerNode {
    /// Create a new planner node.
    pub fn new(config: Vec<f64>, parent: usize, cost: f64) -> Self {
        Self {
            config,
            parent,
            cost,
        }
    }
}
/// Neural MPC controller with receding horizon optimisation.
#[derive(Debug, Clone)]
pub struct ModelPredictiveControl {
    /// Learned dynamics model.
    pub dynamics: NeuralDynamicsModel,
    /// Planning horizon (steps).
    pub horizon: usize,
    /// Number of random shooting candidates.
    pub num_candidates: usize,
    /// Discount factor over horizon.
    pub gamma: f64,
    /// State dimension.
    pub state_dim: usize,
    /// Action dimension.
    pub action_dim: usize,
    /// Action bound (±).
    pub action_bound: f64,
}
impl ModelPredictiveControl {
    /// Create an MPC controller.
    pub fn new(
        state_dim: usize,
        action_dim: usize,
        hidden: &[usize],
        horizon: usize,
        num_candidates: usize,
        gamma: f64,
        action_bound: f64,
    ) -> Self {
        let dynamics = NeuralDynamicsModel::new(state_dim, action_dim, hidden);
        Self {
            dynamics,
            horizon,
            num_candidates,
            gamma,
            state_dim,
            action_dim,
            action_bound,
        }
    }
    /// Random shooting: sample action sequences and return the first action of the best.
    ///
    /// The cost function is supplied as a closure `cost(state) -> f64` (lower is better).
    pub fn plan<F>(&self, state: &[f64], cost_fn: F) -> Vec<f64>
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut rng = rand::rng();
        let mut best_cost = f64::INFINITY;
        let mut best_first = vec![0.0_f64; self.action_dim];
        for _ in 0..self.num_candidates {
            let mut s = state.to_vec();
            let first_action: Vec<f64> = (0..self.action_dim)
                .map(|_| rng.random_range(-self.action_bound..self.action_bound))
                .collect();
            let mut total_cost = 0.0;
            let mut discount = 1.0;
            let mut action = first_action.clone();
            for _ in 0..self.horizon {
                s = self.dynamics.predict(&s, &action);
                total_cost += discount * cost_fn(&s);
                discount *= self.gamma;
                action = (0..self.action_dim)
                    .map(|_| rng.random_range(-self.action_bound..self.action_bound))
                    .collect();
            }
            if total_cost < best_cost {
                best_cost = total_cost;
                best_first = first_action;
            }
        }
        best_first
    }
    /// Cross-entropy method (CEM) planning: iteratively refine action distribution.
    pub fn cem_plan<F>(&self, state: &[f64], cost_fn: F, iterations: usize) -> Vec<f64>
    where
        F: Fn(&[f64]) -> f64,
    {
        let mut rng = rand::rng();
        let d = self.action_dim;
        let mut mean = vec![0.0_f64; d];
        let mut sigma = vec![self.action_bound; d];
        let elite_frac = 0.2_f64;
        for _ in 0..iterations {
            let mut samples: Vec<(Vec<f64>, f64)> = (0..self.num_candidates)
                .map(|_| {
                    let a: Vec<f64> = (0..d)
                        .map(|k| {
                            let z: f64 = rng.random_range(-1.0_f64..1.0_f64);
                            (mean[k] + sigma[k] * z).clamp(-self.action_bound, self.action_bound)
                        })
                        .collect();
                    let mut s = state.to_vec();
                    let mut cost = 0.0;
                    let mut disc = 1.0;
                    let mut act = a.clone();
                    for _ in 0..self.horizon {
                        s = self.dynamics.predict(&s, &act);
                        cost += disc * cost_fn(&s);
                        disc *= self.gamma;
                        act = (0..d)
                            .map(|_| rng.random_range(-self.action_bound..self.action_bound))
                            .collect();
                    }
                    (a, cost)
                })
                .collect();
            samples.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let n_elite = ((self.num_candidates as f64 * elite_frac) as usize).max(1);
            let elites: Vec<&Vec<f64>> = samples[..n_elite].iter().map(|(a, _)| a).collect();
            for k in 0..d {
                let m: f64 = elites.iter().map(|a| a[k]).sum::<f64>() / n_elite as f64;
                let v: f64 =
                    elites.iter().map(|a| (a[k] - m).powi(2)).sum::<f64>() / n_elite as f64;
                mean[k] = m;
                sigma[k] = v.sqrt().max(1e-4);
            }
        }
        mean
    }
}
/// Ornstein-Uhlenbeck noise for temporally correlated exploration.
#[derive(Debug, Clone)]
pub struct OUNoise {
    /// Current state.
    pub state: Vec<f64>,
    /// Mean reversion rate θ.
    pub theta: f64,
    /// Diffusion coefficient σ.
    pub sigma: f64,
    /// Long-run mean μ.
    pub mu: Vec<f64>,
}
impl OUNoise {
    /// Create a zero-mean OU process.
    pub fn new(dim: usize, theta: f64, sigma: f64) -> Self {
        Self {
            state: vec![0.0; dim],
            theta,
            sigma,
            mu: vec![0.0; dim],
        }
    }
    /// Step the process and return the current noise sample.
    pub fn sample(&mut self) -> Vec<f64> {
        let mut rng = rand::rng();
        for i in 0..self.state.len() {
            let dx = self.theta * (self.mu[i] - self.state[i])
                + self.sigma * rng.random_range(-1.0_f64..1.0_f64);
            self.state[i] += dx;
        }
        self.state.clone()
    }
    /// Reset to zero.
    pub fn reset(&mut self) {
        for s in &mut self.state {
            *s = 0.0;
        }
    }
}
/// Circular replay buffer for experience replay.
#[derive(Debug, Clone)]
pub struct ReplayBuffer {
    /// Stored transitions.
    pub buffer: Vec<Transition>,
    /// Maximum capacity.
    pub capacity: usize,
    /// Write position.
    pub write_pos: usize,
}
impl ReplayBuffer {
    /// Create a replay buffer with the given capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            buffer: Vec::with_capacity(capacity),
            capacity,
            write_pos: 0,
        }
    }
    /// Push a transition, overwriting the oldest if full.
    pub fn push(&mut self, t: Transition) {
        if self.buffer.len() < self.capacity {
            self.buffer.push(t);
        } else {
            self.buffer[self.write_pos] = t;
        }
        self.write_pos = (self.write_pos + 1) % self.capacity;
    }
    /// Current number of stored transitions.
    pub fn len(&self) -> usize {
        self.buffer.len()
    }
    /// True if empty.
    pub fn is_empty(&self) -> bool {
        self.buffer.is_empty()
    }
    /// Sample a mini-batch of size `n` (without replacement if possible).
    pub fn sample(&self, n: usize) -> Vec<&Transition> {
        let mut rng = rand::rng();
        let len = self.buffer.len();
        if len == 0 {
            return Vec::new();
        }
        (0..n.min(len))
            .map(|_| &self.buffer[rng.random_range(0..len)])
            .collect()
    }
}
/// RRT / RRT* planner for rigid body configuration spaces.
#[derive(Debug, Clone)]
pub struct MotionPlanning {
    /// Dimension of configuration space.
    pub config_dim: usize,
    /// Lower bounds of the configuration space.
    pub lower: Vec<f64>,
    /// Upper bounds of the configuration space.
    pub upper: Vec<f64>,
    /// Step size (extend length).
    pub step_size: f64,
    /// Maximum number of iterations.
    pub max_iters: usize,
    /// Goal tolerance.
    pub goal_tol: f64,
    /// Whether to use RRT* (optimal variant).
    pub use_rrt_star: bool,
    /// Rewiring radius for RRT*.
    pub rewire_radius: f64,
}
impl MotionPlanning {
    /// Construct a planner.
    pub fn new(
        config_dim: usize,
        lower: Vec<f64>,
        upper: Vec<f64>,
        step_size: f64,
        max_iters: usize,
        goal_tol: f64,
        use_rrt_star: bool,
    ) -> Self {
        let rewire_radius = step_size * 3.0;
        Self {
            config_dim,
            lower,
            upper,
            step_size,
            max_iters,
            goal_tol,
            use_rrt_star,
            rewire_radius,
        }
    }
    /// Sample a random configuration within bounds.
    pub fn random_config(&self) -> Vec<f64> {
        let mut rng = rand::rng();
        (0..self.config_dim)
            .map(|i| rng.random_range(self.lower[i]..self.upper[i]))
            .collect()
    }
    /// Find the nearest node index to `q` in `nodes`.
    pub fn nearest(&self, nodes: &[PlannerNode], q: &[f64]) -> usize {
        nodes
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                config_dist(&a.config, q)
                    .partial_cmp(&config_dist(&b.config, q))
                    .expect("operation should succeed")
            })
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Steer from `from` towards `to` by at most `step_size`.
    pub fn steer(&self, from: &[f64], to: &[f64]) -> Vec<f64> {
        let d = config_dist(from, to);
        if d <= self.step_size {
            to.to_vec()
        } else {
            let t = self.step_size / d;
            from.iter()
                .zip(to.iter())
                .map(|(a, b)| a + t * (b - a))
                .collect()
        }
    }
    /// Stub collision check (always free — replace with real geometry test).
    pub fn is_collision_free(&self, _config: &[f64]) -> bool {
        true
    }
    /// Run RRT and return the path from start to goal, or empty if not found.
    pub fn rrt(&self, start: &[f64], goal: &[f64]) -> Vec<Vec<f64>> {
        let mut rng = rand::rng();
        let mut nodes = vec![PlannerNode::new(start.to_vec(), usize::MAX, 0.0)];
        for _ in 0..self.max_iters {
            let q_rand = if rng.random_range(0.0..1.0) < 0.1 {
                goal.to_vec()
            } else {
                self.random_config()
            };
            let nearest_idx = self.nearest(&nodes, &q_rand);
            let q_new = self.steer(&nodes[nearest_idx].config.clone(), &q_rand);
            if !self.is_collision_free(&q_new) {
                continue;
            }
            let cost = nodes[nearest_idx].cost + config_dist(&nodes[nearest_idx].config, &q_new);
            let new_idx = nodes.len();
            nodes.push(PlannerNode::new(q_new.clone(), nearest_idx, cost));
            if config_dist(&q_new, goal) < self.goal_tol {
                let mut path = Vec::new();
                let mut idx = new_idx;
                while idx != usize::MAX {
                    path.push(nodes[idx].config.clone());
                    idx = nodes[idx].parent;
                    if idx == usize::MAX {
                        break;
                    }
                }
                path.reverse();
                path.push(goal.to_vec());
                return path;
            }
            if self.use_rrt_star {
                let near_indices: Vec<usize> = nodes
                    .iter()
                    .enumerate()
                    .filter(|(i, n)| {
                        *i < new_idx && config_dist(&n.config, &q_new) < self.rewire_radius
                    })
                    .map(|(i, _)| i)
                    .collect();
                for ni in near_indices {
                    let new_cost = nodes[ni].cost + config_dist(&nodes[ni].config, &q_new);
                    if new_cost < nodes[new_idx].cost {
                        nodes[new_idx].parent = ni;
                        nodes[new_idx].cost = new_cost;
                    }
                }
            }
        }
        Vec::new()
    }
    /// Probabilistic roadmap (PRM) construction: returns adjacency list.
    pub fn build_prm(&self, n_samples: usize, k_neighbours: usize) -> Vec<Vec<usize>> {
        let mut configs: Vec<Vec<f64>> = (0..n_samples).map(|_| self.random_config()).collect();
        configs[0] = configs[0].clone();
        let mut adj = vec![Vec::<usize>::new(); n_samples];
        for i in 0..n_samples {
            let mut dists: Vec<(usize, f64)> = (0..n_samples)
                .filter(|&j| j != i)
                .map(|j| (j, config_dist(&configs[i], &configs[j])))
                .collect();
            dists.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            for (j, _) in dists.iter().take(k_neighbours) {
                if self.is_collision_free(&configs[*j]) {
                    adj[i].push(*j);
                }
            }
        }
        adj
    }
    /// Trajectory optimisation: gradient descent on a path to minimise length.
    pub fn optimise_trajectory(
        &self,
        path: Vec<Vec<f64>>,
        iterations: usize,
        lr: f64,
    ) -> Vec<Vec<f64>> {
        let n = path.len();
        if n < 3 {
            return path;
        }
        let mut opt = path;
        for _ in 0..iterations {
            let mut new_path = opt.clone();
            for i in 1..(n - 1) {
                let grad: Vec<f64> = (0..self.config_dim)
                    .map(|k| 2.0 * opt[i][k] - opt[i - 1][k] - opt[i + 1][k])
                    .collect();
                for k in 0..self.config_dim {
                    new_path[i][k] -= lr * grad[k];
                    new_path[i][k] = new_path[i][k].clamp(self.lower[k], self.upper[k]);
                }
            }
            opt = new_path;
        }
        opt
    }
}
/// Q-table for discrete state–action spaces (tabular Q-learning).
#[derive(Debug, Clone)]
pub struct QTable {
    /// Q\[state\]\[action\].
    pub q: Vec<Vec<f64>>,
    /// Learning rate α.
    pub alpha: f64,
    /// Discount factor γ.
    pub gamma: f64,
    /// Exploration rate ε (ε-greedy).
    pub epsilon: f64,
}
impl QTable {
    /// Initialise a Q-table with zeros.
    pub fn new(
        num_states: usize,
        num_actions: usize,
        alpha: f64,
        gamma: f64,
        epsilon: f64,
    ) -> Self {
        Self {
            q: vec![vec![0.0; num_actions]; num_states],
            alpha,
            gamma,
            epsilon,
        }
    }
    /// ε-greedy action selection.
    pub fn select_action(&self, state: usize) -> usize {
        let mut rng = rand::rng();
        if rng.random_range(0.0..1.0) < self.epsilon {
            rng.random_range(0..self.q[state].len())
        } else {
            let row = &self.q[state];
            row.iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(i, _)| i)
                .unwrap_or(0)
        }
    }
    /// Tabular Q-learning update (SARSA-max / Q-learning).
    pub fn update(
        &mut self,
        state: usize,
        action: usize,
        reward: f64,
        next_state: usize,
        done: bool,
    ) {
        let max_next = if done {
            0.0
        } else {
            self.q[next_state]
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let td_target = reward + self.gamma * max_next;
        let td_error = td_target - self.q[state][action];
        self.q[state][action] += self.alpha * td_error;
    }
    /// Decay epsilon multiplicatively.
    pub fn decay_epsilon(&mut self, factor: f64) {
        self.epsilon = (self.epsilon * factor).max(0.01);
    }
}
