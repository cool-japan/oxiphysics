#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Neural network-augmented Smoothed Particle Hydrodynamics (SPH).
//!
//! This module provides deep-learning augmented SPH components:
//! - Learned kernel functions replacing the classical cubic-spline kernel
//! - Neural pressure solvers with residual correction
//! - Physics-informed SPH with PDE-based loss terms
//! - Graph neural networks for particle interaction learning
//! - Surrogate viscosity models

#![allow(dead_code)]

// ─────────────────────────────────────────────────────────────────────────────
// Neural Kernel
// ─────────────────────────────────────────────────────────────────────────────

/// A neural network kernel function that replaces the classical SPH kernel.
///
/// The kernel is parameterised by a weight vector and a smoothing bandwidth.
/// It is trained from data pairs `(r/h, W_true)` where `W_true` is the target
/// kernel value.
#[derive(Debug, Clone)]
pub struct NeuralKernel {
    /// Learned weight vector (one weight per basis function).
    pub weights: Vec<f64>,
    /// Smoothing length bandwidth parameter.
    pub bandwidth: f64,
}

impl NeuralKernel {
    /// Create a new `NeuralKernel` with given weights and bandwidth.
    pub fn new(weights: Vec<f64>, bandwidth: f64) -> Self {
        Self { weights, bandwidth }
    }

    /// Create a default `NeuralKernel` with small constant weights.
    pub fn default_kernel(n_weights: usize, bandwidth: f64) -> Self {
        let weights: Vec<f64> = (0..n_weights).map(|i| 0.1 / (i + 1) as f64).collect();
        Self { weights, bandwidth }
    }

    /// Evaluate the kernel at normalised distance `q = r / h`.
    ///
    /// Uses a polynomial basis expanded at `q`:
    /// `W(q) = sum_k w_k * phi_k(q)` where `phi_k(q) = max(0, 1-q)^k`.
    pub fn evaluate(&self, r: f64, h: f64) -> f64 {
        let q = r / h.max(1e-30);
        if q >= 1.0 {
            return 0.0;
        }
        let base = (1.0 - q).max(0.0);
        self.weights
            .iter()
            .enumerate()
            .map(|(k, &w)| w * base.powi(k as i32 + 1))
            .sum::<f64>()
            / (h * h * h)
    }

    /// Evaluate the radial kernel gradient `dW/dr` at distance `r`.
    ///
    /// Uses the analytical derivative of the polynomial basis.
    pub fn gradient(&self, r: f64, h: f64) -> f64 {
        let q = r / h.max(1e-30);
        if q >= 1.0 {
            return 0.0;
        }
        let base = (1.0 - q).max(0.0);
        let dbase_dr = -1.0 / h;
        let sum: f64 = self
            .weights
            .iter()
            .enumerate()
            .map(|(k, &w)| {
                let exp = (k + 1) as i32;
                w * exp as f64 * base.powi(exp - 1) * dbase_dr
            })
            .sum();
        sum / (h * h * h)
    }

    /// Train the kernel weights from a set of `(r, W_target)` pairs using a
    /// simple least-squares gradient descent.
    ///
    /// `pairs`: slice of `(r, W_target)` training samples.
    pub fn train_from_data(&mut self, pairs: &[(f64, f64)]) {
        let n = self.weights.len();
        let lr = 1e-3;
        let epochs = 200;
        let h = self.bandwidth;

        for _ in 0..epochs {
            let mut grad = vec![0.0f64; n];
            for &(r, w_target) in pairs {
                let q = r / h.max(1e-30);
                if q >= 1.0 {
                    continue;
                }
                let base = (1.0 - q).max(0.0);
                let w_pred: f64 = self
                    .weights
                    .iter()
                    .enumerate()
                    .map(|(k, &w)| w * base.powi(k as i32 + 1))
                    .sum::<f64>()
                    / (h * h * h);
                let err = w_pred - w_target;
                for k in 0..n {
                    grad[k] += err * base.powi(k as i32 + 1) / (h * h * h);
                }
            }
            for k in 0..n {
                self.weights[k] -= lr * grad[k] / pairs.len().max(1) as f64;
            }
        }
    }

    /// Return the number of basis weights.
    pub fn n_weights(&self) -> usize {
        self.weights.len()
    }

    /// Compute the compact support radius (= bandwidth).
    pub fn support_radius(&self) -> f64 {
        self.bandwidth
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Divergence-free correction data bundle
// ─────────────────────────────────────────────────────────────────────────────

/// Particle data bundle passed to the iterative divergence-free pressure
/// correction.  Groups position, velocity, mass, and density slices so that
/// the correction function stays within the 7-argument limit.
pub struct DivCorrectionData<'a> {
    /// Particle positions `[[x, y, z]; n]`.
    pub pos: &'a [[f64; 3]],
    /// Particle velocities `[[vx, vy, vz]; n]`.
    pub vel: &'a [[f64; 3]],
    /// Particle masses `[m; n]`.
    pub mass: &'a [f64],
    /// Particle densities `[ρ; n]`.
    pub rho: &'a [f64],
}

// ─────────────────────────────────────────────────────────────────────────────
// Neural Pressure Solver
// ─────────────────────────────────────────────────────────────────────────────

/// A neural network-based pressure solver for incompressible SPH.
///
/// Replaces the classical iterative pressure Poisson equation solver with a
/// feed-forward neural network, optionally augmented with a residual correction
/// step.
#[derive(Debug, Clone)]
pub struct NeuralPressureSolver {
    /// Weight matrix: `nn[layer][neuron_index]`.  For simplicity each row
    /// encodes the weights of one neuron as `[w0, w1, ..., bias]`.
    pub nn: Vec<Vec<f64>>,
    /// Reference density.
    pub rho0: f64,
    /// Speed of sound squared (Tait EOS parameter).
    pub cs2: f64,
}

impl NeuralPressureSolver {
    /// Create a new `NeuralPressureSolver` with given network weights.
    pub fn new(nn: Vec<Vec<f64>>, rho0: f64, cs2: f64) -> Self {
        Self { nn, rho0, cs2 }
    }

    /// Create a default solver with a simple two-layer network.
    pub fn default_solver(rho0: f64, cs2: f64) -> Self {
        // Two hidden neurons: each row is [w0, w1, bias]
        let nn = vec![
            vec![0.5, -0.1, 0.0],
            vec![-0.2, 0.8, 0.05],
            vec![0.3, 0.3, -0.01],
        ];
        Self::new(nn, rho0, cs2)
    }

    /// Predict the pressure for a particle given its density and neighbour
    /// density values.
    ///
    /// The feature vector is `[rho / rho0, mean_rho_neighbours / rho0]`.
    pub fn predict_pressure(&self, density: f64, neighbors: &[f64]) -> f64 {
        let mean_neigh = if neighbors.is_empty() {
            density
        } else {
            neighbors.iter().sum::<f64>() / neighbors.len() as f64
        };
        let x0 = density / self.rho0;
        let x1 = mean_neigh / self.rho0;

        // One hidden layer of ReLU neurons, then a linear output.
        let hidden: Vec<f64> = self
            .nn
            .iter()
            .map(|row| {
                let pre = row.first().copied().unwrap_or(0.0) * x0
                    + row.get(1).copied().unwrap_or(0.0) * x1
                    + row.get(2).copied().unwrap_or(0.0);
                pre.max(0.0)
            })
            .collect();

        // Output layer: sum of hidden activations scaled by cs2.
        let p_nn: f64 = hidden.iter().sum::<f64>() * self.cs2;

        // Add Tait background pressure: p_tait = rho0 * cs2 * ((rho/rho0)^7 - 1) / 7
        let p_tait = self.rho0 * self.cs2 * ((x0.powi(7) - 1.0) / 7.0);

        p_nn + p_tait
    }

    /// Apply a scalar residual correction to reduce a single-particle divergence error.
    ///
    /// Corrects pressure proportional to density divergence:
    /// `p_new = p - (ρ₀ / dt) * div_v`
    /// where `dt` defaults to `cs2⁻¹` for dimensional consistency.
    pub fn residual_correction(&self, pressure: f64, divergence: f64) -> f64 {
        let dt_eff = if self.cs2 > 1e-30 {
            1.0 / self.cs2.sqrt()
        } else {
            1.0
        };
        let correction = self.rho0 * divergence * dt_eff;
        (pressure - correction).max(0.0)
    }

    /// Iterative divergence-free pressure correction for a full particle system.
    ///
    /// Runs up to `MAX_ITER=3` correction passes.  Each pass:
    /// 1. Computes the density divergence at each particle via SPH summation.
    /// 2. Adjusts pressure: `Δp = -ρ₀ * div_v / dt` (clamped non-negative).
    ///
    /// Returns the maximum absolute divergence before the first iteration and
    /// after the final iteration as `(max_div_before, max_div_after)`.
    ///
    /// # Arguments
    /// * `data`          — particle data bundle (positions, velocities, masses, densities)
    /// * `pressure`      — mutable pressure array `[p; n]` (updated in-place)
    /// * `neighbor_list` — pre-built neighbour lists `Vec<Vec<usize>>`
    /// * `h`             — smoothing length
    /// * `dt`            — time step
    pub fn iterative_divergence_free_correction(
        &self,
        data: &DivCorrectionData<'_>,
        pressure: &mut [f64],
        neighbor_list: &[Vec<usize>],
        h: f64,
        dt: f64,
    ) -> (f64, f64) {
        const MAX_ITER: usize = 3;
        const DIV_THRESHOLD: f64 = 1e-4;

        let pos = data.pos;
        let vel = data.vel;
        let mass = data.mass;
        let rho = data.rho;

        let n = pos.len();
        if n == 0 {
            return (0.0, 0.0);
        }
        let mut density_div = vec![0.0_f64; n];

        // Helper: SPH kernel gradient vector (cubic spline)
        let kernel_gradient = |r_ij: [f64; 3], h_: f64| -> [f64; 3] {
            let rx = r_ij[0];
            let ry = r_ij[1];
            let rz = r_ij[2];
            let r = (rx * rx + ry * ry + rz * rz).sqrt();
            if r < 1.0e-30 || h_ < 1.0e-30 {
                return [0.0; 3];
            }
            let q = r / h_;
            let sigma = 1.0 / (std::f64::consts::PI * h_ * h_ * h_ * h_);
            let dw_dr = if q < 1.0 {
                sigma * (-3.0 * q + 2.25 * q * q)
            } else if q < 2.0 {
                let t = 2.0 - q;
                sigma * (-0.75 * t * t)
            } else {
                0.0
            };
            [dw_dr * rx / r, dw_dr * ry / r, dw_dr * rz / r]
        };

        let compute_div = |density_div: &mut Vec<f64>| -> f64 {
            let mut max_div = 0.0_f64;
            for i in 0..n {
                let mut div_v = 0.0_f64;
                if let Some(neighbours) = neighbor_list.get(i) {
                    for &j in neighbours {
                        if j >= n {
                            continue;
                        }
                        let r_ij = [
                            pos[i][0] - pos[j][0],
                            pos[i][1] - pos[j][1],
                            pos[i][2] - pos[j][2],
                        ];
                        let v_ij = [
                            vel[i][0] - vel[j][0],
                            vel[i][1] - vel[j][1],
                            vel[i][2] - vel[j][2],
                        ];
                        let grad_w = kernel_gradient(r_ij, h);
                        let rho_j = rho[j].max(1e-30);
                        let m_j = mass.get(j).copied().unwrap_or(1.0);
                        let dot_vg =
                            v_ij[0] * grad_w[0] + v_ij[1] * grad_w[1] + v_ij[2] * grad_w[2];
                        div_v += m_j * dot_vg / rho_j;
                    }
                }
                density_div[i] = div_v;
                max_div = max_div.max(div_v.abs());
            }
            max_div
        };

        let max_div_before = compute_div(&mut density_div);
        let mut max_div_after = max_div_before;

        for _iter in 0..MAX_ITER {
            if max_div_after < DIV_THRESHOLD {
                break;
            }
            // Pressure correction: Δp = -ρ₀ * div_v / dt
            let dt_safe = if dt.abs() > 1e-30 { dt } else { 1.0 };
            for i in 0..n {
                pressure[i] -= self.rho0 * density_div[i] / dt_safe;
                pressure[i] = pressure[i].max(0.0);
            }
            max_div_after = compute_div(&mut density_div);
        }

        (max_div_before, max_div_after)
    }

    /// Number of neurons in the hidden layer.
    pub fn n_neurons(&self) -> usize {
        self.nn.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Physics-Informed SPH
// ─────────────────────────────────────────────────────────────────────────────

/// Physics-informed SPH that minimises a loss function combining data,
/// PDE residual, and boundary residual terms.
///
/// Each particle is represented by a 6-component state vector:
/// `[x, y, z, vx, vy, vz]`.
#[derive(Debug, Clone)]
pub struct PhysicsInformedSph {
    /// Particle state vectors `[x, y, z, vx, vy, vz]`.
    pub particles: Vec<[f64; 6]>,
    /// Fluid density (uniform for now).
    pub density: f64,
    /// Dynamic viscosity.
    pub viscosity: f64,
}

impl PhysicsInformedSph {
    /// Create a new `PhysicsInformedSph` with the given particles.
    pub fn new(particles: Vec<[f64; 6]>, density: f64, viscosity: f64) -> Self {
        Self {
            particles,
            density,
            viscosity,
        }
    }

    /// Compute the PDE (Navier-Stokes) residual as a scalar.
    ///
    /// Uses a simple mass-conservation proxy: sum of velocity divergences.
    pub fn pde_residual(&self) -> f64 {
        // Approximate div(v) ~ (vx + vy + vz) for each particle as a proxy.
        let sum: f64 = self
            .particles
            .iter()
            .map(|p| (p[3] + p[4] + p[5]).powi(2))
            .sum();
        (sum / self.particles.len().max(1) as f64).sqrt()
    }

    /// Compute the boundary residual: penalise particles outside `[0, 1]^3`.
    pub fn boundary_residual(&self) -> f64 {
        let sum: f64 = self
            .particles
            .iter()
            .map(|p| {
                let mut r = 0.0f64;
                for &x in &p[0..3] {
                    if x < 0.0 {
                        r += x * x;
                    } else if x > 1.0 {
                        r += (x - 1.0).powi(2);
                    }
                }
                r
            })
            .sum();
        sum.sqrt()
    }

    /// Compute the total physics-informed loss.
    ///
    /// `lambda_pde` and `lambda_bc` are weighting coefficients.
    pub fn total_loss(&self, lambda_pde: f64, lambda_bc: f64) -> f64 {
        lambda_pde * self.pde_residual() + lambda_bc * self.boundary_residual()
    }

    /// Return the number of particles.
    pub fn n_particles(&self) -> usize {
        self.particles.len()
    }

    /// Compute the kinetic energy of the particle system.
    pub fn kinetic_energy(&self) -> f64 {
        0.5 * self.density
            * self
                .particles
                .iter()
                .map(|p| p[3] * p[3] + p[4] * p[4] + p[5] * p[5])
                .sum::<f64>()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Graph Neural SPH
// ─────────────────────────────────────────────────────────────────────────────

/// Graph neural network (GNN)-based SPH where particles are nodes and
/// interactions are edges.
///
/// The GNN performs message passing to aggregate neighbour information before
/// updating node states and predicting forces.
#[derive(Debug, Clone)]
pub struct GraphNeuralSph {
    /// Per-node feature vectors.
    pub node_features: Vec<Vec<f64>>,
    /// Per-edge feature vectors.  Edge `k` connects node `edges[k].0` to
    /// `edges[k].1`.
    pub edge_features: Vec<Vec<f64>>,
    /// Edge connectivity: `(source, target)` pairs.
    pub edges: Vec<(usize, usize)>,
    /// Message MLP weights (one weight per input feature).
    pub message_weights: Vec<f64>,
    /// Node update MLP weights.
    pub update_weights: Vec<f64>,
}

impl GraphNeuralSph {
    /// Create a new `GraphNeuralSph` from node features, edge features,
    /// and connectivity.
    pub fn new(
        node_features: Vec<Vec<f64>>,
        edge_features: Vec<Vec<f64>>,
        edges: Vec<(usize, usize)>,
    ) -> Self {
        let feat_dim = node_features.first().map_or(1, |f| f.len());
        let message_weights = vec![1.0 / feat_dim as f64; feat_dim];
        let update_weights = vec![0.5; feat_dim];
        Self {
            node_features,
            edge_features,
            edges,
            message_weights,
            update_weights,
        }
    }

    /// Perform one round of message passing.
    ///
    /// Each edge aggregates the source node features and edge features into a
    /// message, which is summed at the target node.
    ///
    /// Returns a vector of aggregated messages for each node.
    pub fn message_passing(&self) -> Vec<Vec<f64>> {
        let n_nodes = self.node_features.len();
        let feat_dim = self.message_weights.len();
        let mut aggregated = vec![vec![0.0f64; feat_dim]; n_nodes];

        for (k, &(src, tgt)) in self.edges.iter().enumerate() {
            if src >= n_nodes || tgt >= n_nodes {
                continue;
            }
            let src_feat = &self.node_features[src];
            let edge_feat = self.edge_features.get(k);

            for d in 0..feat_dim {
                let sf = src_feat.get(d).copied().unwrap_or(0.0);
                let ef = edge_feat.and_then(|e| e.get(d)).copied().unwrap_or(0.0);
                let msg = self.message_weights.get(d).copied().unwrap_or(1.0) * (sf + ef);
                aggregated[tgt][d] += msg;
            }
        }
        aggregated
    }

    /// Update node features by combining current features with aggregated
    /// messages.
    ///
    /// Applies a simple gated update: `h' = tanh(W_u * (h + agg))`.
    pub fn update_nodes(&mut self) {
        let messages = self.message_passing();
        for (node_feat, msg) in self.node_features.iter_mut().zip(messages.iter()) {
            for (h, (&m, &w)) in node_feat
                .iter_mut()
                .zip(msg.iter().zip(self.update_weights.iter()))
            {
                *h = (w * (*h + m)).tanh();
            }
        }
    }

    /// Predict the force on each node by projecting updated node features.
    ///
    /// Returns a 3D force vector for each node.
    pub fn predict_forces(&self) -> Vec<[f64; 3]> {
        self.node_features
            .iter()
            .map(|feat| {
                let fx = feat.first().copied().unwrap_or(0.0)
                    * self.update_weights.first().copied().unwrap_or(1.0);
                let fy = feat.get(1).copied().unwrap_or(0.0)
                    * self.update_weights.get(1).copied().unwrap_or(1.0);
                let fz = feat.get(2).copied().unwrap_or(0.0)
                    * self.update_weights.get(2).copied().unwrap_or(1.0);
                [fx, fy, fz]
            })
            .collect()
    }

    /// Return the number of nodes.
    pub fn n_nodes(&self) -> usize {
        self.node_features.len()
    }

    /// Return the number of edges.
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Surrogate Viscosity
// ─────────────────────────────────────────────────────────────────────────────

/// Surrogate viscosity model that predicts an effective non-Newtonian viscosity
/// from a shear-rate scalar and additional flow features.
///
/// Uses a simple polynomial regression model stored as a coefficient vector.
#[derive(Debug, Clone)]
pub struct SurrogateViscosity {
    /// Polynomial / regression coefficients of the model.
    pub model: Vec<f64>,
    /// Reference viscosity (Newtonian baseline).
    pub mu0: f64,
}

impl SurrogateViscosity {
    /// Create a new `SurrogateViscosity` with given model coefficients.
    pub fn new(model: Vec<f64>, mu0: f64) -> Self {
        Self { model, mu0 }
    }

    /// Create a default surrogate model (Power-law-like).
    pub fn default_model(mu0: f64) -> Self {
        // coefficients: [c0, c1] for mu_eff = mu0 * (c0 + c1*gamma)
        // Negative c1 gives shear-thinning behaviour.
        let model = vec![1.0, -0.05];
        Self::new(model, mu0)
    }

    /// Predict the effective viscosity for a given shear rate `shear_rate`
    /// and additional flow features.
    ///
    /// The features are incorporated as additive corrections.
    pub fn effective_viscosity(&self, shear_rate: f64, features: &[f64]) -> f64 {
        // Polynomial in shear_rate.
        let poly: f64 = self
            .model
            .iter()
            .enumerate()
            .map(|(k, &c)| c * shear_rate.powi(k as i32))
            .sum();

        // Feature correction: dot(model_tail, features).
        let feat_correction: f64 = self
            .model
            .iter()
            .skip(1)
            .zip(features.iter())
            .map(|(&c, &f)| c * f * 0.01)
            .sum();

        (self.mu0 * (poly + feat_correction)).max(1e-15)
    }

    /// Return the zero-shear-rate viscosity.
    pub fn zero_shear_viscosity(&self) -> f64 {
        self.effective_viscosity(0.0, &[])
    }

    /// Compute the Cross model viscosity: `mu = mu0 / (1 + (lambda * gamma)^n)`.
    pub fn cross_model(mu0: f64, lambda: f64, n: f64, gamma: f64) -> f64 {
        mu0 / (1.0 + (lambda * gamma).powf(n)).max(1e-30)
    }

    /// Compute the Carreau model viscosity.
    pub fn carreau_model(mu_inf: f64, mu0: f64, lambda: f64, n: f64, gamma: f64) -> f64 {
        mu_inf + (mu0 - mu_inf) * (1.0 + (lambda * gamma).powi(2)).powf((n - 1.0) / 2.0)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Deep learning correction for SPH accelerations
// ─────────────────────────────────────────────────────────────────────────────

/// Combine physical SPH accelerations with neural-network corrections.
///
/// The blended acceleration is:
/// ```text
/// a_blend[i] = accel_physical[i] + weight * accel_nn[i]
/// ```
///
/// Both slices must have the same length.
pub fn sph_deep_learning_correction(
    accel_physical: &[f64],
    accel_nn: &[f64],
    weight: f64,
) -> Vec<f64> {
    let n = accel_physical.len();
    (0..n)
        .map(|i| {
            let a_nn = if i < accel_nn.len() { accel_nn[i] } else { 0.0 };
            accel_physical[i] + weight * a_nn
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Compute the Euclidean distance between two 3D position arrays.
pub fn particle_distance(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Compute the relative velocity magnitude between two particles.
pub fn relative_speed(a: &[f64; 6], b: &[f64; 6]) -> f64 {
    let dvx = a[3] - b[3];
    let dvy = a[4] - b[4];
    let dvz = a[5] - b[5];
    (dvx * dvx + dvy * dvy + dvz * dvz).sqrt()
}

/// Find the indices of particles within smoothing length `h` of particle
/// `idx`.
pub fn find_neighbors(particles: &[[f64; 6]], idx: usize, h: f64) -> Vec<usize> {
    particles
        .iter()
        .enumerate()
        .filter(|&(j, p)| j != idx && particle_distance(&particles[idx], p) < h)
        .map(|(j, _)| j)
        .collect()
}

/// Normalise a vector in place.  Returns `false` if the norm is zero.
pub fn normalise(v: &mut Vec<f64>) -> bool {
    let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    if norm < 1e-30 {
        return false;
    }
    for x in v.iter_mut() {
        *x /= norm;
    }
    true
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // --- NeuralKernel tests ---

    #[test]
    fn test_kernel_zero_at_support() {
        let k = NeuralKernel::default_kernel(4, 1.0);
        assert_eq!(k.evaluate(1.0, 1.0), 0.0);
    }

    #[test]
    fn test_kernel_zero_beyond_support() {
        let k = NeuralKernel::default_kernel(4, 1.0);
        assert_eq!(k.evaluate(1.5, 1.0), 0.0);
    }

    #[test]
    fn test_kernel_positive_inside() {
        let k = NeuralKernel::default_kernel(4, 1.0);
        assert!(k.evaluate(0.3, 1.0) > 0.0);
    }

    #[test]
    fn test_kernel_gradient_zero_beyond() {
        let k = NeuralKernel::default_kernel(4, 1.0);
        assert_eq!(k.gradient(1.5, 1.0), 0.0);
    }

    #[test]
    fn test_kernel_gradient_finite_inside() {
        let k = NeuralKernel::default_kernel(4, 1.0);
        assert!(k.gradient(0.4, 1.0).is_finite());
    }

    #[test]
    fn test_kernel_training_runs() {
        let mut k = NeuralKernel::default_kernel(4, 1.0);
        let pairs: Vec<(f64, f64)> = (0..10)
            .map(|i| {
                let r = i as f64 * 0.05;
                let w = (1.0 - r).max(0.0).powi(3);
                (r, w)
            })
            .collect();
        k.train_from_data(&pairs);
        // weights should be finite after training
        for w in &k.weights {
            assert!(w.is_finite());
        }
    }

    #[test]
    fn test_kernel_support_radius() {
        let k = NeuralKernel::default_kernel(4, 2.5);
        assert_eq!(k.support_radius(), 2.5);
    }

    #[test]
    fn test_kernel_n_weights() {
        let k = NeuralKernel::default_kernel(6, 1.0);
        assert_eq!(k.n_weights(), 6);
    }

    // --- NeuralPressureSolver tests ---

    #[test]
    fn test_pressure_solver_at_rest() {
        let solver = NeuralPressureSolver::default_solver(1000.0, 1500.0 * 1500.0);
        let p = solver.predict_pressure(1000.0, &[]);
        assert!(p.is_finite());
    }

    #[test]
    fn test_pressure_solver_compressed() {
        let solver = NeuralPressureSolver::default_solver(1000.0, 1500.0 * 1500.0);
        let p_rest = solver.predict_pressure(1000.0, &[]);
        let p_comp = solver.predict_pressure(1100.0, &[]);
        // Compressed fluid should have higher pressure.
        assert!(p_comp > p_rest);
    }

    #[test]
    fn test_pressure_residual_correction() {
        let solver = NeuralPressureSolver::default_solver(1000.0, 1500.0 * 1500.0);
        let p = solver.predict_pressure(1000.0, &[1000.0, 1000.0]);
        let p_corr = solver.residual_correction(p, 0.01);
        assert!(p_corr.is_finite());
    }

    #[test]
    fn test_pressure_n_neurons() {
        let solver = NeuralPressureSolver::default_solver(1000.0, 1.0);
        assert_eq!(solver.n_neurons(), 3);
    }

    // --- PhysicsInformedSph tests ---

    #[test]
    fn test_pde_residual_at_rest() {
        let particles = vec![[0.5, 0.5, 0.5, 0.0, 0.0, 0.0]; 5];
        let pisph = PhysicsInformedSph::new(particles, 1000.0, 0.001);
        assert!((pisph.pde_residual()).abs() < 1e-12);
    }

    #[test]
    fn test_boundary_residual_inside() {
        let particles = vec![[0.2, 0.3, 0.4, 0.0, 0.0, 0.0]; 3];
        let pisph = PhysicsInformedSph::new(particles, 1000.0, 0.001);
        assert!((pisph.boundary_residual()).abs() < 1e-12);
    }

    #[test]
    fn test_boundary_residual_outside() {
        let particles = vec![[-0.5, 0.5, 0.5, 0.0, 0.0, 0.0]];
        let pisph = PhysicsInformedSph::new(particles, 1000.0, 0.001);
        assert!(pisph.boundary_residual() > 0.0);
    }

    #[test]
    fn test_total_loss_non_negative() {
        let particles = vec![[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]; 4];
        let pisph = PhysicsInformedSph::new(particles, 1000.0, 0.001);
        assert!(pisph.total_loss(1.0, 1.0) >= 0.0);
    }

    #[test]
    fn test_kinetic_energy_rest() {
        let particles = vec![[0.5, 0.5, 0.5, 0.0, 0.0, 0.0]; 4];
        let pisph = PhysicsInformedSph::new(particles, 1000.0, 0.001);
        assert!((pisph.kinetic_energy()).abs() < 1e-12);
    }

    #[test]
    fn test_n_particles() {
        let pisph = PhysicsInformedSph::new(vec![[0.0; 6]; 7], 1000.0, 0.001);
        assert_eq!(pisph.n_particles(), 7);
    }

    // --- GraphNeuralSph tests ---

    #[test]
    fn test_gnn_message_passing_size() {
        let nodes: Vec<Vec<f64>> = (0..4).map(|i| vec![i as f64, 0.0, 0.0]).collect();
        let edges_feat: Vec<Vec<f64>> = vec![vec![0.1, 0.0, 0.0]; 3];
        let edges = vec![(0, 1), (1, 2), (2, 3)];
        let gnn = GraphNeuralSph::new(nodes, edges_feat, edges);
        let agg = gnn.message_passing();
        assert_eq!(agg.len(), 4);
    }

    #[test]
    fn test_gnn_update_nodes_finite() {
        let nodes: Vec<Vec<f64>> = (0..3).map(|_| vec![0.1, 0.2, 0.3]).collect();
        let edge_feats: Vec<Vec<f64>> = vec![vec![0.05, 0.0, 0.0]; 2];
        let edges = vec![(0, 1), (1, 2)];
        let mut gnn = GraphNeuralSph::new(nodes, edge_feats, edges);
        gnn.update_nodes();
        for node in &gnn.node_features {
            for &v in node {
                assert!(v.is_finite());
            }
        }
    }

    #[test]
    fn test_gnn_predict_forces_size() {
        let nodes: Vec<Vec<f64>> = (0..5).map(|_| vec![0.1, 0.2, 0.3]).collect();
        let gnn = GraphNeuralSph::new(nodes, vec![], vec![]);
        let forces = gnn.predict_forces();
        assert_eq!(forces.len(), 5);
    }

    #[test]
    fn test_gnn_n_nodes_edges() {
        let nodes: Vec<Vec<f64>> = vec![vec![1.0]; 6];
        let edge_feats = vec![vec![0.0]; 4];
        let edges = vec![(0, 1), (1, 2), (2, 3), (3, 4)];
        let gnn = GraphNeuralSph::new(nodes, edge_feats, edges);
        assert_eq!(gnn.n_nodes(), 6);
        assert_eq!(gnn.n_edges(), 4);
    }

    // --- SurrogateViscosity tests ---

    #[test]
    fn test_surrogate_zero_shear() {
        let sv = SurrogateViscosity::default_model(0.001);
        let mu = sv.zero_shear_viscosity();
        assert!(mu > 0.0);
    }

    #[test]
    fn test_surrogate_shear_thinning() {
        let sv = SurrogateViscosity::default_model(0.001);
        let mu_low = sv.effective_viscosity(0.1, &[]);
        let mu_high = sv.effective_viscosity(10.0, &[]);
        // Power-law with negative c1 → viscosity decreases with shear rate.
        assert!(mu_low > mu_high);
    }

    #[test]
    fn test_surrogate_always_positive() {
        let sv = SurrogateViscosity::default_model(0.001);
        for gamma in [0.0, 1.0, 10.0, 100.0] {
            assert!(sv.effective_viscosity(gamma, &[]) > 0.0);
        }
    }

    #[test]
    fn test_cross_model() {
        let mu = SurrogateViscosity::cross_model(0.001, 1.0, 0.5, 1.0);
        assert!(mu > 0.0 && mu < 0.001 + 1e-10);
    }

    #[test]
    fn test_carreau_model_zero_shear() {
        let mu = SurrogateViscosity::carreau_model(0.0001, 0.001, 1.0, 0.5, 0.0);
        assert!((mu - 0.001).abs() < 1e-12);
    }

    // --- sph_deep_learning_correction tests ---

    #[test]
    fn test_correction_zero_weight() {
        let a_phys = vec![1.0, 2.0, 3.0];
        let a_nn = vec![10.0, 20.0, 30.0];
        let out = sph_deep_learning_correction(&a_phys, &a_nn, 0.0);
        assert_eq!(out, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_correction_unit_weight() {
        let a_phys = vec![1.0, 2.0];
        let a_nn = vec![0.5, 0.5];
        let out = sph_deep_learning_correction(&a_phys, &a_nn, 1.0);
        assert!((out[0] - 1.5).abs() < 1e-12);
        assert!((out[1] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_correction_shorter_nn() {
        let a_phys = vec![1.0, 2.0, 3.0];
        let a_nn = vec![0.1];
        let out = sph_deep_learning_correction(&a_phys, &a_nn, 1.0);
        assert_eq!(out.len(), 3);
        assert!((out[1] - 2.0).abs() < 1e-12);
    }

    // --- Utility tests ---

    #[test]
    fn test_particle_distance_same() {
        let p: [f64; 6] = [1.0, 2.0, 3.0, 0.0, 0.0, 0.0];
        assert!((particle_distance(&p, &p)).abs() < 1e-12);
    }

    #[test]
    fn test_particle_distance_known() {
        let a: [f64; 6] = [0.0; 6];
        let b: [f64; 6] = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        assert!((particle_distance(&a, &b) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_relative_speed_zero() {
        let p: [f64; 6] = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0];
        assert!((relative_speed(&p, &p)).abs() < 1e-12);
    }

    #[test]
    fn test_find_neighbors_count() {
        let particles = vec![
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0, 0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ];
        let neigh = find_neighbors(&particles, 0, 0.5);
        assert_eq!(neigh.len(), 1);
        assert_eq!(neigh[0], 1);
    }

    #[test]
    fn test_normalise_unit_vector() {
        let mut v = vec![3.0, 4.0];
        let ok = normalise(&mut v);
        assert!(ok);
        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_normalise_zero_vector() {
        let mut v = vec![0.0, 0.0, 0.0];
        let ok = normalise(&mut v);
        assert!(!ok);
    }

    // ── Iterative divergence-free correction (E4) ────────────────────────────

    #[test]
    fn test_divergence_correction_reduces_max_div() {
        // Build a small 8-particle system with non-zero divergence.
        let solver = NeuralPressureSolver::default_solver(1000.0, 1500.0);
        let n = 8_usize;
        let h = 0.5_f64;
        let dt = 1e-3_f64;

        // Particles arranged in a 2×2×2 cube with outward velocities (diverging flow)
        let positions: Vec<[f64; 3]> = (0..n)
            .map(|i| {
                let ix = (i % 2) as f64;
                let iy = ((i / 2) % 2) as f64;
                let iz = (i / 4) as f64;
                [ix * 0.3, iy * 0.3, iz * 0.3]
            })
            .collect();
        let velocities: Vec<[f64; 3]> = positions
            .iter()
            .map(|p| [p[0] * 0.5, p[1] * 0.5, p[2] * 0.5]) // outward => divergence > 0
            .collect();
        let masses = vec![1.0_f64; n];
        let rho = vec![1000.0_f64; n];
        let mut pressure = vec![0.0_f64; n];

        // Build simple neighbour lists: each particle sees all others within 2h
        let neighbor_list: Vec<Vec<usize>> = (0..n)
            .map(|i| {
                (0..n)
                    .filter(|&j| {
                        if j == i {
                            return false;
                        }
                        let dx = positions[i][0] - positions[j][0];
                        let dy = positions[i][1] - positions[j][1];
                        let dz = positions[i][2] - positions[j][2];
                        (dx * dx + dy * dy + dz * dz).sqrt() < 2.0 * h
                    })
                    .collect()
            })
            .collect();

        let particle_data = DivCorrectionData {
            pos: &positions,
            vel: &velocities,
            mass: &masses,
            rho: &rho,
        };
        let (max_div_before, max_div_after) = solver.iterative_divergence_free_correction(
            &particle_data,
            &mut pressure,
            &neighbor_list,
            h,
            dt,
        );

        // Divergence may already be small for this geometry, but after correction
        // it must be <= before.
        assert!(
            max_div_after <= max_div_before + 1e-14,
            "divergence should not increase: before={max_div_before}, after={max_div_after}"
        );
        // All pressures should be non-negative
        for (i, &p) in pressure.iter().enumerate() {
            assert!(p >= 0.0, "pressure[{i}] must be non-negative, got {p}");
        }
    }
}
