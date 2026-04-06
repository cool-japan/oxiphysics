//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

#![allow(clippy::needless_range_loop, clippy::ptr_arg)]
use super::functions::*;
use std::f64::consts::PI;

/// K-means clustering result.
#[derive(Clone, Debug)]
pub struct KMeansResult {
    /// Cluster centroids \[k × d\].
    pub centroids: Vec<Vec<f64>>,
    /// Cluster assignment for each data point.
    pub labels: Vec<usize>,
    /// Within-cluster sum of squares.
    pub inertia: f64,
}
/// A single node in a decision tree.
#[derive(Clone, Debug)]
pub enum TreeNode {
    /// Internal split node.
    Split {
        /// Feature index to split on.
        feature: usize,
        /// Threshold value.
        threshold: f64,
        /// Left child (feature <= threshold).
        left: Box<TreeNode>,
        /// Right child (feature > threshold).
        right: Box<TreeNode>,
    },
    /// Leaf node with a class label (most common label) or regression value.
    Leaf {
        /// Predicted value.
        value: f64,
    },
}
impl TreeNode {
    /// Predict for a single sample.
    pub fn predict(&self, x: &[f64]) -> f64 {
        match self {
            TreeNode::Leaf { value } => *value,
            TreeNode::Split {
                feature,
                threshold,
                left,
                right,
            } => {
                if x[*feature] <= *threshold {
                    left.predict(x)
                } else {
                    right.predict(x)
                }
            }
        }
    }
}
/// CART decision tree classifier.
#[derive(Clone, Debug)]
pub struct DecisionTree {
    pub(super) root: Option<TreeNode>,
    /// Maximum tree depth.
    pub max_depth: usize,
    /// Minimum samples to split.
    pub min_samples_split: usize,
}
impl DecisionTree {
    /// Create a new decision tree with specified hyperparameters.
    pub fn new(max_depth: usize, min_samples_split: usize) -> Self {
        Self {
            root: None,
            max_depth,
            min_samples_split,
        }
    }
    /// Fit the tree to data.
    pub fn fit(&mut self, x: &[Vec<f64>], y: &[f64]) {
        self.root = Some(build_tree(x, y, 0, self.max_depth, self.min_samples_split));
    }
    /// Predict a single sample.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        self.root.as_ref().map(|r| r.predict(x)).unwrap_or(0.0)
    }
    /// Predict a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
    /// Compute accuracy on a labelled dataset.
    pub fn accuracy(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let preds = self.predict(x);
        let correct = preds
            .iter()
            .zip(y.iter())
            .filter(|(p, t)| (*p - *t).abs() < 0.5)
            .count();
        correct as f64 / y.len() as f64
    }
}
/// Min-max scaler.
#[derive(Clone, Debug)]
pub struct MinMaxScaler {
    /// Minimum of each feature.
    pub min: Vec<f64>,
    /// Range (max - min) of each feature.
    pub range: Vec<f64>,
}
impl MinMaxScaler {
    /// Fit to data.
    pub fn fit(data: &[Vec<f64>]) -> Self {
        let d = data[0].len();
        let mut min = vec![f64::MAX; d];
        let mut max = vec![f64::MIN; d];
        for row in data {
            for j in 0..d {
                if row[j] < min[j] {
                    min[j] = row[j];
                }
                if row[j] > max[j] {
                    max[j] = row[j];
                }
            }
        }
        let range: Vec<f64> = min
            .iter()
            .zip(max.iter())
            .map(|(mn, mx)| (mx - mn).max(1e-300))
            .collect();
        Self { min, range }
    }
    /// Transform data.
    pub fn transform(&self, data: &[Vec<f64>]) -> Vec<Vec<f64>> {
        data.iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, &v)| (v - self.min[j]) / self.range[j])
                    .collect()
            })
            .collect()
    }
    /// Inverse transform.
    pub fn inverse_transform(&self, data: &[Vec<f64>]) -> Vec<Vec<f64>> {
        data.iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, &v)| v * self.range[j] + self.min[j])
                    .collect()
            })
            .collect()
    }
}
/// Standard scaler (z-score normalisation).
#[derive(Clone, Debug)]
pub struct StandardScaler {
    /// Feature means.
    pub mean: Vec<f64>,
    /// Feature standard deviations.
    pub std: Vec<f64>,
}
impl StandardScaler {
    /// Fit to data.
    pub fn fit(data: &[Vec<f64>]) -> Self {
        let mean = mean_vec(data);
        let var = variance_vec(data);
        let std: Vec<f64> = var.iter().map(|v| v.sqrt().max(1e-300)).collect();
        Self { mean, std }
    }
    /// Transform data.
    pub fn transform(&self, data: &[Vec<f64>]) -> Vec<Vec<f64>> {
        data.iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, &v)| (v - self.mean[j]) / self.std[j])
                    .collect()
            })
            .collect()
    }
    /// Inverse transform.
    pub fn inverse_transform(&self, data: &[Vec<f64>]) -> Vec<Vec<f64>> {
        data.iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .map(|(j, &v)| v * self.std[j] + self.mean[j])
                    .collect()
            })
            .collect()
    }
}
/// Distance metric for KNN.
#[derive(Clone, Debug, Copy, PartialEq)]
pub enum KnnMetric {
    /// Euclidean (L2) distance.
    Euclidean,
    /// Manhattan (L1) distance.
    Manhattan,
    /// Chebyshev (L∞) distance.
    Chebyshev,
}
/// Isolation Forest anomaly detector.
///
/// Anomaly score ≈ 1 means highly anomalous; ≈ 0.5 means normal.
#[derive(Clone, Debug)]
pub struct IsolationForest {
    /// Collection of isolation trees.
    pub trees: Vec<IsolationNode>,
    /// Subsample size used per tree.
    pub subsample_size: usize,
    /// Number of trees.
    pub n_trees: usize,
}
impl IsolationForest {
    /// Fit the isolation forest.
    pub fn fit(x: &[Vec<f64>], n_trees: usize, subsample_size: usize) -> Self {
        let n = x.len();
        let ss = subsample_size.min(n);
        let max_depth = ((ss as f64).log2().ceil() as usize).max(1);
        let mut trees = Vec::with_capacity(n_trees);
        let mut seed: u64 = 42_u64.wrapping_mul(6364136223846793005);
        for _ in 0..n_trees {
            let mut sub_indices: Vec<usize> = Vec::with_capacity(ss);
            for _ in 0..ss {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                sub_indices.push((seed >> 33) as usize % n);
            }
            let sub_x: Vec<Vec<f64>> = sub_indices.iter().map(|&i| x[i].clone()).collect();
            trees.push(build_isolation_tree(&sub_x, max_depth, 0, &mut seed));
        }
        Self {
            trees,
            subsample_size: ss,
            n_trees,
        }
    }
    /// Anomaly score for a single sample.
    ///
    /// Score in `(0, 1]`: higher values indicate anomalies.
    pub fn anomaly_score(&self, x: &[f64]) -> f64 {
        let avg_path: f64 = self
            .trees
            .iter()
            .map(|t| t.path_length(x, 0.0))
            .sum::<f64>()
            / self.n_trees as f64;
        let c = average_path_length(self.subsample_size);
        if c < 1e-10 {
            return 0.5;
        }
        0.5_f64.powf(avg_path / c)
    }
    /// Predict anomaly scores for a batch.
    pub fn score_samples(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.anomaly_score(row)).collect()
    }
    /// Classify as anomaly (true) if score > threshold.
    pub fn predict_anomaly(&self, x: &[Vec<f64>], threshold: f64) -> Vec<bool> {
        self.score_samples(x)
            .iter()
            .map(|&s| s > threshold)
            .collect()
    }
}
/// t-SNE dimensionality reduction result.
#[derive(Clone, Debug)]
pub struct TsneResult {
    /// 2D embedding coordinates (n × 2).
    pub embedding: Vec<[f64; 2]>,
    /// Final KL divergence cost.
    pub kl_divergence: f64,
}
/// A trainable dense (fully-connected) layer.
#[derive(Clone, Debug)]
pub struct DenseLayer {
    /// Weight matrix \[out × in\].
    pub weights: Matrix,
    /// Bias vector \[out\].
    pub bias: Vec<f64>,
    /// Cached input (for backprop).
    pub input_cache: Vec<f64>,
    /// Cached pre-activation (for backprop).
    pub z_cache: Vec<f64>,
}
impl DenseLayer {
    /// Construct with Xavier/Glorot uniform initialisation.
    pub fn new(in_dim: usize, out_dim: usize) -> Self {
        let limit = (6.0 / (in_dim + out_dim) as f64).sqrt();
        let mut data = vec![0.0f64; out_dim * in_dim];
        for (i, v) in data.iter_mut().enumerate() {
            let t = (i as f64 * 1.3).sin();
            *v = limit * t;
        }
        let weights = Matrix::from_vec(out_dim, in_dim, data);
        Self {
            weights,
            bias: vec![0.0; out_dim],
            input_cache: Vec::new(),
            z_cache: Vec::new(),
        }
    }
    /// Forward pass.  Returns the pre-activation output z = Wx + b.
    pub fn forward(&mut self, x: &[f64]) -> Vec<f64> {
        self.input_cache = x.to_vec();
        let z = self.weights.matvec(x);
        let z: Vec<f64> = z
            .iter()
            .zip(self.bias.iter())
            .map(|(zi, bi)| zi + bi)
            .collect();
        self.z_cache = z.clone();
        z
    }
    /// Backward pass.  Returns gradient w.r.t. input.
    pub fn backward(&mut self, grad_output: &[f64], lr: f64) -> Vec<f64> {
        let in_dim = self.weights.cols;
        let out_dim = self.weights.rows;
        for i in 0..out_dim {
            for j in 0..in_dim {
                let dw = grad_output[i] * self.input_cache[j];
                let w = self.weights.get(i, j) - lr * dw;
                self.weights.set(i, j, w);
            }
            self.bias[i] -= lr * grad_output[i];
        }
        let wt = self.weights.transpose();
        wt.matvec(grad_output)
    }
}
/// Out-of-bag error result for a random forest.
#[derive(Clone, Debug)]
pub struct OobResult {
    /// Out-of-bag predictions (one per training sample).
    pub oob_predictions: Vec<f64>,
    /// OOB accuracy (for classifiers: fraction of correct predictions).
    pub oob_accuracy: f64,
    /// Bootstrap index sets used per tree.
    pub bootstrap_sets: Vec<Vec<usize>>,
}
/// Random forest classifier.
#[derive(Clone, Debug)]
pub struct RandomForest {
    /// Individual decision trees.
    pub trees: Vec<DecisionTree>,
    /// Number of features to consider at each split (sqrt heuristic).
    pub n_features: usize,
    /// Number of trees.
    pub n_estimators: usize,
}
impl RandomForest {
    /// Create a new random forest.
    pub fn new(n_estimators: usize, max_depth: usize, min_samples_split: usize) -> Self {
        let trees = (0..n_estimators)
            .map(|_| DecisionTree::new(max_depth, min_samples_split))
            .collect();
        Self {
            trees,
            n_features: 0,
            n_estimators,
        }
    }
    /// Fit the random forest using bootstrap sampling.
    pub fn fit(&mut self, x: &[Vec<f64>], y: &[f64]) {
        let n = x.len();
        let d = x[0].len();
        self.n_features = ((d as f64).sqrt() as usize).max(1);
        for (t, tree) in self.trees.iter_mut().enumerate() {
            let mut indices = Vec::with_capacity(n);
            let mut seed = (t as u64 + 1).wrapping_mul(6364136223846793005);
            for _ in 0..n {
                seed = seed
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                indices.push((seed >> 33) as usize % n);
            }
            let bx: Vec<Vec<f64>> = indices.iter().map(|&i| x[i].clone()).collect();
            let by: Vec<f64> = indices.iter().map(|&i| y[i]).collect();
            tree.fit(&bx, &by);
        }
    }
    /// Predict a single sample by majority vote.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        let votes: Vec<f64> = self.trees.iter().map(|t| t.predict_one(x)).collect();
        majority_vote(&votes)
    }
    /// Predict a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
    /// Compute accuracy.
    pub fn accuracy(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let preds = self.predict(x);
        let correct = preds
            .iter()
            .zip(y.iter())
            .filter(|(p, t)| (*p - *t).abs() < 0.5)
            .count();
        correct as f64 / y.len() as f64
    }
}
/// Linear regression model.
#[derive(Clone, Debug)]
pub struct LinearRegression {
    /// Learned weights (coefficients).
    pub weights: Vec<f64>,
    /// Learned bias (intercept).
    pub bias: f64,
}
impl LinearRegression {
    /// Fit using normal equations: w = (X^T X)^{-1} X^T y.
    ///
    /// Adds a bias column automatically.
    pub fn fit_normal(x: &[Vec<f64>], y: &[f64]) -> Self {
        let n = x.len();
        let d = x[0].len() + 1;
        let mut xmat = Matrix::zeros(n, d);
        for (i, row) in x.iter().enumerate() {
            xmat.set(i, 0, 1.0);
            for (j, &v) in row.iter().enumerate() {
                xmat.set(i, j + 1, v);
            }
        }
        let xt = xmat.transpose();
        let xtx = xt.matmul(&xmat);
        let xty = xt.matvec(y);
        let w = solve_linear_system(&xtx, &xty).unwrap_or_else(|| vec![0.0; d]);
        Self {
            weights: w[1..].to_vec(),
            bias: w[0],
        }
    }
    /// Fit using gradient descent.
    pub fn fit_gradient_descent(x: &[Vec<f64>], y: &[f64], lr: f64, epochs: usize) -> Self {
        let n = x.len();
        let d = x[0].len();
        let mut w = vec![0.0f64; d];
        let mut b = 0.0f64;
        for _ in 0..epochs {
            let mut dw = vec![0.0f64; d];
            let mut db = 0.0f64;
            for (xi, &yi) in x.iter().zip(y.iter()) {
                let pred = dot(xi, &w) + b;
                let err = pred - yi;
                for j in 0..d {
                    dw[j] += err * xi[j];
                }
                db += err;
            }
            let inv_n = lr / n as f64;
            for j in 0..d {
                w[j] -= inv_n * dw[j];
            }
            b -= inv_n * db;
        }
        Self {
            weights: w,
            bias: b,
        }
    }
    /// Predict for a single sample.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        dot(x, &self.weights) + self.bias
    }
    /// Predict for a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
    /// Compute R² score.
    pub fn r2_score(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let preds = self.predict(x);
        let mean_y = y.iter().sum::<f64>() / y.len() as f64;
        let ss_res: f64 = preds
            .iter()
            .zip(y.iter())
            .map(|(p, t)| (p - t).powi(2))
            .sum();
        let ss_tot: f64 = y.iter().map(|t| (t - mean_y).powi(2)).sum();
        if ss_tot == 0.0 {
            1.0
        } else {
            1.0 - ss_res / ss_tot
        }
    }
    /// Mean squared error.
    pub fn mse(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let preds = self.predict(x);
        preds
            .iter()
            .zip(y.iter())
            .map(|(p, t)| (p - t).powi(2))
            .sum::<f64>()
            / y.len() as f64
    }
}
/// Gaussian Naive Bayes classifier.
#[derive(Clone, Debug)]
pub struct GaussianNaiveBayes {
    /// Class prior probabilities.
    pub priors: Vec<f64>,
    /// Per-class feature means \[class\]\[feature\].
    pub means: Vec<Vec<f64>>,
    /// Per-class feature variances \[class\]\[feature\].
    pub variances: Vec<Vec<f64>>,
    /// Number of classes.
    pub n_classes: usize,
}
impl GaussianNaiveBayes {
    /// Fit the model.
    pub fn fit(x: &[Vec<f64>], y: &[usize], n_classes: usize) -> Self {
        let n = x.len();
        let d = x[0].len();
        let mut means = vec![vec![0.0f64; d]; n_classes];
        let mut variances = vec![vec![0.0f64; d]; n_classes];
        let mut counts = vec![0usize; n_classes];
        for (xi, &yi) in x.iter().zip(y.iter()) {
            if yi < n_classes {
                counts[yi] += 1;
                for j in 0..d {
                    means[yi][j] += xi[j];
                }
            }
        }
        for c in 0..n_classes {
            let cnt = counts[c].max(1) as f64;
            for j in 0..d {
                means[c][j] /= cnt;
            }
        }
        for (xi, &yi) in x.iter().zip(y.iter()) {
            if yi < n_classes {
                for j in 0..d {
                    variances[yi][j] += (xi[j] - means[yi][j]).powi(2);
                }
            }
        }
        for c in 0..n_classes {
            let cnt = counts[c].max(1) as f64;
            for j in 0..d {
                variances[c][j] = variances[c][j] / cnt + 1e-9;
            }
        }
        let total = n as f64;
        let priors = counts.iter().map(|&c| c as f64 / total).collect();
        Self {
            priors,
            means,
            variances,
            n_classes,
        }
    }
    /// Predict log-probabilities for a single sample.
    pub fn predict_log_proba(&self, x: &[f64]) -> Vec<f64> {
        (0..self.n_classes)
            .map(|c| {
                let mut log_p = self.priors[c].ln();
                for j in 0..x.len() {
                    let mu = self.means[c][j];
                    let var = self.variances[c][j];
                    log_p += -0.5 * ((x[j] - mu).powi(2) / var + (2.0 * PI * var).ln());
                }
                log_p
            })
            .collect()
    }
    /// Predict class label for a single sample.
    pub fn predict_one(&self, x: &[f64]) -> usize {
        self.predict_log_proba(x)
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0)
    }
    /// Predict for a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<usize> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
}
/// PCA result.
#[derive(Clone, Debug)]
pub struct PcaResult {
    /// Principal components (eigenvectors) stored row-wise \[n_components × d\].
    pub components: Vec<Vec<f64>>,
    /// Eigenvalues (variance explained per component).
    pub eigenvalues: Vec<f64>,
    /// Mean of training data.
    pub mean: Vec<f64>,
}
impl PcaResult {
    /// Project data onto the principal components.
    pub fn transform(&self, data: &[Vec<f64>]) -> Vec<Vec<f64>> {
        data.iter()
            .map(|row| {
                let centered: Vec<f64> = row
                    .iter()
                    .zip(self.mean.iter())
                    .map(|(x, m)| x - m)
                    .collect();
                self.components
                    .iter()
                    .map(|pc| dot(&centered, pc))
                    .collect()
            })
            .collect()
    }
    /// Reconstruct data from projected coordinates.
    pub fn inverse_transform(&self, projected: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let d = self.mean.len();
        projected
            .iter()
            .map(|proj| {
                let mut rec = self.mean.clone();
                for (k, &coef) in proj.iter().enumerate() {
                    for j in 0..d {
                        rec[j] += coef * self.components[k][j];
                    }
                }
                rec
            })
            .collect()
    }
    /// Fraction of total variance explained by the retained components.
    pub fn explained_variance_ratio(&self) -> Vec<f64> {
        let total: f64 = self.eigenvalues.iter().sum();
        self.eigenvalues.iter().map(|e| e / total).collect()
    }
}
/// Gradient Boosted Decision Trees regressor (additive model with stumps).
///
/// Uses the squared-error loss function; negative gradient = residual.
#[derive(Clone, Debug)]
pub struct GradientBoosting {
    /// Ensemble of regression stumps.
    pub stumps: Vec<RegressionStump>,
    /// Initial prediction (mean of target).
    pub init_pred: f64,
    /// Shrinkage learning rate in `(0, 1]`.
    pub learning_rate: f64,
    /// Number of boosting rounds.
    pub n_estimators: usize,
}
impl GradientBoosting {
    /// Fit the GBDT model on `(x, y)` data.
    ///
    /// # Arguments
    /// * `x` — feature matrix (n × d)
    /// * `y` — target values
    /// * `n_estimators` — number of boosting rounds
    /// * `learning_rate` — shrinkage factor
    pub fn fit(x: &[Vec<f64>], y: &[f64], n_estimators: usize, learning_rate: f64) -> Self {
        let n = y.len();
        let init_pred = y.iter().sum::<f64>() / n as f64;
        let mut preds = vec![init_pred; n];
        let mut stumps = Vec::with_capacity(n_estimators);
        for _ in 0..n_estimators {
            let residuals: Vec<f64> = y.iter().zip(preds.iter()).map(|(t, p)| t - p).collect();
            let stump = RegressionStump::fit(x, &residuals);
            let step = stump.predict(x);
            for (p, s) in preds.iter_mut().zip(step.iter()) {
                *p += learning_rate * s;
            }
            stumps.push(stump);
        }
        Self {
            stumps,
            init_pred,
            learning_rate,
            n_estimators,
        }
    }
    /// Predict for a single sample.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        let mut pred = self.init_pred;
        for stump in &self.stumps {
            pred += self.learning_rate * stump.predict_one(x);
        }
        pred
    }
    /// Predict for a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
    /// Mean squared error on the training set.
    pub fn train_mse(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let preds = self.predict(x);
        mse(y, &preds)
    }
}
/// Binary logistic regression.
#[derive(Clone, Debug)]
pub struct LogisticRegression {
    /// Learned weights.
    pub weights: Vec<f64>,
    /// Learned bias.
    pub bias: f64,
}
impl LogisticRegression {
    /// Fit using gradient descent (binary cross-entropy loss).
    pub fn fit(x: &[Vec<f64>], y: &[f64], lr: f64, epochs: usize) -> Self {
        let n = x.len();
        let d = x[0].len();
        let mut w = vec![0.0f64; d];
        let mut b = 0.0f64;
        for _ in 0..epochs {
            let mut dw = vec![0.0f64; d];
            let mut db = 0.0f64;
            for (xi, &yi) in x.iter().zip(y.iter()) {
                let z = dot(xi, &w) + b;
                let pred = sigmoid(z);
                let err = pred - yi;
                for j in 0..d {
                    dw[j] += err * xi[j];
                }
                db += err;
            }
            let inv_n = lr / n as f64;
            for j in 0..d {
                w[j] -= inv_n * dw[j];
            }
            b -= inv_n * db;
        }
        Self {
            weights: w,
            bias: b,
        }
    }
    /// Predict probability for a single sample.
    pub fn predict_proba_one(&self, x: &[f64]) -> f64 {
        sigmoid(dot(x, &self.weights) + self.bias)
    }
    /// Predict class label (0 or 1) for a single sample.
    pub fn predict_one(&self, x: &[f64]) -> usize {
        if self.predict_proba_one(x) >= 0.5 {
            1
        } else {
            0
        }
    }
    /// Predict class labels for a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<usize> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
    /// Compute binary cross-entropy loss.
    pub fn log_loss(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let n = x.len() as f64;
        x.iter()
            .zip(y.iter())
            .map(|(xi, &yi)| {
                let p = self.predict_proba_one(xi).clamp(1e-15, 1.0 - 1e-15);
                -(yi * p.ln() + (1.0 - yi) * (1.0 - p).ln())
            })
            .sum::<f64>()
            / n
    }
    /// Classification accuracy.
    pub fn accuracy(&self, x: &[Vec<f64>], y: &[f64]) -> f64 {
        let preds = self.predict(x);
        let correct = preds
            .iter()
            .zip(y.iter())
            .filter(|(p, t)| **p as f64 == **t)
            .count();
        correct as f64 / y.len() as f64
    }
}
/// A simple dense matrix stored row-major.
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    /// Number of rows.
    pub rows: usize,
    /// Number of columns.
    pub cols: usize,
    /// Data stored row-major.
    pub data: Vec<f64>,
}
impl Matrix {
    /// Construct a zero matrix of shape `rows × cols`.
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {
            rows,
            cols,
            data: vec![0.0; rows * cols],
        }
    }
    /// Construct an identity matrix of size `n × n`.
    pub fn eye(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = 1.0;
        }
        m
    }
    /// Construct from row-major flat data.
    pub fn from_vec(rows: usize, cols: usize, data: Vec<f64>) -> Self {
        assert_eq!(data.len(), rows * cols);
        Self { rows, cols, data }
    }
    /// Return the element at `(row, col)`.
    #[inline]
    pub fn get(&self, row: usize, col: usize) -> f64 {
        self.data[row * self.cols + col]
    }
    /// Set the element at `(row, col)`.
    #[inline]
    pub fn set(&mut self, row: usize, col: usize, val: f64) {
        self.data[row * self.cols + col] = val;
    }
    /// Transpose.
    pub fn transpose(&self) -> Self {
        let mut out = Self::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                out.set(j, i, self.get(i, j));
            }
        }
        out
    }
    /// Matrix-matrix multiplication.
    pub fn matmul(&self, rhs: &Self) -> Self {
        assert_eq!(self.cols, rhs.rows, "Inner dimensions must match");
        let mut out = Self::zeros(self.rows, rhs.cols);
        for i in 0..self.rows {
            for k in 0..self.cols {
                let aik = self.get(i, k);
                for j in 0..rhs.cols {
                    let val = out.get(i, j) + aik * rhs.get(k, j);
                    out.set(i, j, val);
                }
            }
        }
        out
    }
    /// Matrix-vector multiplication.
    pub fn matvec(&self, v: &[f64]) -> Vec<f64> {
        assert_eq!(self.cols, v.len());
        (0..self.rows)
            .map(|i| (0..self.cols).map(|j| self.get(i, j) * v[j]).sum())
            .collect()
    }
    /// Add two matrices.
    pub fn add(&self, rhs: &Self) -> Self {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);
        let data = self
            .data
            .iter()
            .zip(rhs.data.iter())
            .map(|(a, b)| a + b)
            .collect();
        Self {
            rows: self.rows,
            cols: self.cols,
            data,
        }
    }
    /// Subtract two matrices.
    pub fn sub(&self, rhs: &Self) -> Self {
        assert_eq!(self.rows, rhs.rows);
        assert_eq!(self.cols, rhs.cols);
        let data = self
            .data
            .iter()
            .zip(rhs.data.iter())
            .map(|(a, b)| a - b)
            .collect();
        Self {
            rows: self.rows,
            cols: self.cols,
            data,
        }
    }
    /// Scale by a scalar.
    pub fn scale(&self, s: f64) -> Self {
        let data = self.data.iter().map(|x| x * s).collect();
        Self {
            rows: self.rows,
            cols: self.cols,
            data,
        }
    }
    /// Return a column as a `Vec`f64`.
    pub fn col(&self, j: usize) -> Vec<f64> {
        (0..self.rows).map(|i| self.get(i, j)).collect()
    }
    /// Return a row as a `Vec`f64`.
    pub fn row(&self, i: usize) -> Vec<f64> {
        self.data[i * self.cols..(i + 1) * self.cols].to_vec()
    }
    /// Frobenius norm.
    pub fn frob_norm(&self) -> f64 {
        self.data.iter().map(|x| x * x).sum::<f64>().sqrt()
    }
}
/// Gaussian Process Regression with RBF kernel.
///
/// Computes predictive mean and variance at test points, and the
/// log marginal likelihood of the observed data.
#[derive(Clone, Debug)]
pub struct GaussianProcessRegressor {
    /// Training inputs.
    pub x_train: Vec<Vec<f64>>,
    /// Training targets.
    pub y_train: Vec<f64>,
    /// Kernel amplitude (signal standard deviation).
    pub amplitude: f64,
    /// Kernel length scale.
    pub length_scale: f64,
    /// Noise standard deviation (observation noise).
    pub noise: f64,
    /// Cholesky factor L such that (K + σ²I) = L Lᵀ.
    pub l_factor: Matrix,
    /// α = L⁻ᵀ L⁻¹ y (precomputed).
    pub alpha: Vec<f64>,
}
impl GaussianProcessRegressor {
    /// Fit the GP to training data.
    ///
    /// Computes the Cholesky decomposition of K + σ²I.
    pub fn fit(
        x_train: Vec<Vec<f64>>,
        y_train: Vec<f64>,
        amplitude: f64,
        length_scale: f64,
        noise: f64,
    ) -> Self {
        let n = x_train.len();
        let mut km = kernel_matrix(&x_train, amplitude, length_scale);
        for i in 0..n {
            let v = km.get(i, i) + noise * noise;
            km.set(i, i, v);
        }
        let l = cholesky(&km);
        let ly = forward_substitution(&l, &y_train);
        let alpha = backward_substitution_t(&l, &ly);
        Self {
            x_train,
            y_train,
            amplitude,
            length_scale,
            noise,
            l_factor: l,
            alpha,
        }
    }
    /// Predict mean and variance at test point `x_star`.
    pub fn predict_one(&self, x_star: &[f64]) -> (f64, f64) {
        let k_star: Vec<f64> = self
            .x_train
            .iter()
            .map(|xi| rbf_kernel(xi, x_star, self.amplitude, self.length_scale))
            .collect();
        let mean = dot(&k_star, &self.alpha);
        let v = forward_substitution(&self.l_factor, &k_star);
        let k_ss = rbf_kernel(x_star, x_star, self.amplitude, self.length_scale);
        let var = (k_ss - dot(&v, &v)).max(0.0);
        (mean, var)
    }
    /// Predict means and variances for a batch.
    pub fn predict(&self, x_test: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>) {
        let (means, vars): (Vec<_>, Vec<_>) = x_test.iter().map(|x| self.predict_one(x)).unzip();
        (means, vars)
    }
    /// Log marginal likelihood: log p(y | X, θ).
    ///
    /// Used for hyperparameter optimisation.
    pub fn log_marginal_likelihood(&self) -> f64 {
        let n = self.y_train.len() as f64;
        let fit = -0.5 * dot(&self.y_train, &self.alpha);
        let log_det: f64 = (0..self.l_factor.rows)
            .map(|i| self.l_factor.get(i, i).max(1e-300).ln())
            .sum::<f64>()
            * 2.0;
        let complexity = -0.5 * log_det;
        let constant = -0.5 * n * (2.0 * PI).ln();
        fit + complexity + constant
    }
}
/// Batch normalisation layer.
#[derive(Clone, Debug)]
pub struct BatchNorm {
    /// Learnable scale parameter.
    pub gamma: Vec<f64>,
    /// Learnable shift parameter.
    pub beta: Vec<f64>,
    /// Running mean (used at inference).
    pub running_mean: Vec<f64>,
    /// Running variance (used at inference).
    pub running_var: Vec<f64>,
    /// Momentum for running statistics.
    pub momentum: f64,
    /// Small constant for numerical stability.
    pub eps: f64,
    /// Training mode flag.
    pub training: bool,
}
impl BatchNorm {
    /// Create with default parameters.
    pub fn new(n_features: usize) -> Self {
        Self {
            gamma: vec![1.0; n_features],
            beta: vec![0.0; n_features],
            running_mean: vec![0.0; n_features],
            running_var: vec![1.0; n_features],
            momentum: 0.1,
            eps: 1e-5,
            training: true,
        }
    }
    /// Forward pass for a batch of samples.
    pub fn forward_batch(&mut self, x: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = x.len() as f64;
        let d = x[0].len();
        let batch_mean: Vec<f64> = (0..d)
            .map(|j| x.iter().map(|r| r[j]).sum::<f64>() / n)
            .collect();
        let batch_var: Vec<f64> = (0..d)
            .map(|j| {
                x.iter()
                    .map(|r| (r[j] - batch_mean[j]).powi(2))
                    .sum::<f64>()
                    / n
            })
            .collect();
        if self.training {
            for j in 0..d {
                self.running_mean[j] =
                    (1.0 - self.momentum) * self.running_mean[j] + self.momentum * batch_mean[j];
                self.running_var[j] =
                    (1.0 - self.momentum) * self.running_var[j] + self.momentum * batch_var[j];
            }
        }
        x.iter()
            .map(|row| {
                (0..d)
                    .map(|j| {
                        let mean = if self.training {
                            batch_mean[j]
                        } else {
                            self.running_mean[j]
                        };
                        let var = if self.training {
                            batch_var[j]
                        } else {
                            self.running_var[j]
                        };
                        let norm = (row[j] - mean) / (var + self.eps).sqrt();
                        self.gamma[j] * norm + self.beta[j]
                    })
                    .collect()
            })
            .collect()
    }
}
/// RMSProp optimiser.
#[derive(Clone, Debug)]
pub struct RmsProp {
    /// Learning rate.
    pub lr: f64,
    /// Decay coefficient.
    pub alpha: f64,
    /// Epsilon.
    pub eps: f64,
    /// Running second moment.
    pub v: Vec<f64>,
}
impl RmsProp {
    /// Create RMSProp optimiser.
    pub fn new(n_params: usize, lr: f64, alpha: f64, eps: f64) -> Self {
        Self {
            lr,
            alpha,
            eps,
            v: vec![0.0; n_params],
        }
    }
    /// Perform one RMSProp step.
    pub fn step(&mut self, params: &mut Vec<f64>, grads: &[f64]) {
        for (i, (p, &g)) in params.iter_mut().zip(grads.iter()).enumerate() {
            self.v[i] = self.alpha * self.v[i] + (1.0 - self.alpha) * g * g;
            *p -= self.lr * g / (self.v[i].sqrt() + self.eps);
        }
    }
}
/// A versatile K-Nearest Neighbours model supporting both classification and regression.
///
/// Weighted voting uses inverse-distance weights: closer neighbours get higher weight.
#[derive(Clone, Debug)]
pub struct KnnModel {
    /// Training feature matrix.
    pub x_train: Vec<Vec<f64>>,
    /// Training targets (class labels or continuous values).
    pub y_train: Vec<f64>,
    /// Number of neighbours `k`.
    pub k: usize,
    /// Distance metric.
    pub metric: KnnMetric,
    /// Whether to use weighted (inverse-distance) voting.
    pub weighted: bool,
}
impl KnnModel {
    /// Create a new KNN model.
    pub fn new(k: usize, metric: KnnMetric, weighted: bool) -> Self {
        Self {
            x_train: Vec::new(),
            y_train: Vec::new(),
            k,
            metric,
            weighted,
        }
    }
    /// Store training data.
    pub fn fit(&mut self, x: Vec<Vec<f64>>, y: Vec<f64>) {
        self.x_train = x;
        self.y_train = y;
    }
    fn dist(&self, a: &[f64], b: &[f64]) -> f64 {
        match self.metric {
            KnnMetric::Euclidean => a
                .iter()
                .zip(b.iter())
                .map(|(u, v)| (u - v).powi(2))
                .sum::<f64>()
                .sqrt(),
            KnnMetric::Manhattan => a.iter().zip(b.iter()).map(|(u, v)| (u - v).abs()).sum(),
            KnnMetric::Chebyshev => a
                .iter()
                .zip(b.iter())
                .map(|(u, v)| (u - v).abs())
                .fold(0.0_f64, f64::max),
        }
    }
    /// Classify: predict class label using (weighted) majority vote.
    pub fn classify(&self, x: &[f64]) -> f64 {
        let mut dists: Vec<(f64, f64)> = self
            .x_train
            .iter()
            .zip(self.y_train.iter())
            .map(|(xi, &yi)| (self.dist(x, xi), yi))
            .collect();
        dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let neighbours = &dists[..self.k.min(dists.len())];
        if !self.weighted {
            majority_vote(&neighbours.iter().map(|(_, l)| *l).collect::<Vec<_>>())
        } else {
            let mut class_weights: std::collections::HashMap<i64, f64> =
                std::collections::HashMap::new();
            for &(d, label) in neighbours {
                let w = 1.0 / (d + 1e-10);
                *class_weights.entry(label as i64).or_insert(0.0) += w;
            }
            class_weights
                .into_iter()
                .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                .map(|(l, _)| l as f64)
                .unwrap_or(0.0)
        }
    }
    /// Regress: predict continuous value using (weighted) average.
    pub fn regress(&self, x: &[f64]) -> f64 {
        let mut dists: Vec<(f64, f64)> = self
            .x_train
            .iter()
            .zip(self.y_train.iter())
            .map(|(xi, &yi)| (self.dist(x, xi), yi))
            .collect();
        dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let neighbours = &dists[..self.k.min(dists.len())];
        if !self.weighted {
            neighbours.iter().map(|(_, v)| v).sum::<f64>() / neighbours.len() as f64
        } else {
            let mut sum_w = 0.0f64;
            let mut sum_wv = 0.0f64;
            for &(d, v) in neighbours {
                let w = 1.0 / (d + 1e-10);
                sum_w += w;
                sum_wv += w * v;
            }
            if sum_w > 0.0 { sum_wv / sum_w } else { 0.0 }
        }
    }
    /// Batch predict (classify).
    pub fn classify_batch(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.classify(row)).collect()
    }
    /// Batch predict (regress).
    pub fn regress_batch(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.regress(row)).collect()
    }
}
/// K-fold cross-validation result.
#[derive(Clone, Debug)]
pub struct CrossValResult {
    /// Score on each fold.
    pub fold_scores: Vec<f64>,
    /// Mean score.
    pub mean_score: f64,
    /// Standard deviation of scores.
    pub std_score: f64,
}
/// A single regression tree stump (depth-1 tree) used in boosting.
///
/// The stump splits on the single best feature/threshold pair to minimise
/// the mean-squared error of residuals.
#[derive(Clone, Debug)]
pub struct RegressionStump {
    /// Feature index used for the split.
    pub feature: usize,
    /// Threshold value.
    pub threshold: f64,
    /// Prediction for samples with feature <= threshold.
    pub left_value: f64,
    /// Prediction for samples with feature > threshold.
    pub right_value: f64,
}
impl RegressionStump {
    /// Fit a stump to residuals `r` given input `x`.
    pub fn fit(x: &[Vec<f64>], r: &[f64]) -> Self {
        let n = x.len();
        let n_features = if n == 0 { 0 } else { x[0].len() };
        let mut best_loss = f64::INFINITY;
        let mut best_feat = 0;
        let mut best_thresh = 0.0;
        let mut best_lv = 0.0;
        let mut best_rv = 0.0;
        for feat in 0..n_features {
            let mut vals: Vec<f64> = x.iter().map(|row| row[feat]).collect();
            vals.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            vals.dedup();
            for thresh in vals.windows(2).map(|w| (w[0] + w[1]) * 0.5) {
                let left_r: Vec<f64> = x
                    .iter()
                    .zip(r.iter())
                    .filter(|(row, _)| row[feat] <= thresh)
                    .map(|(_, &rv)| rv)
                    .collect();
                let right_r: Vec<f64> = x
                    .iter()
                    .zip(r.iter())
                    .filter(|(row, _)| row[feat] > thresh)
                    .map(|(_, &rv)| rv)
                    .collect();
                if left_r.is_empty() || right_r.is_empty() {
                    continue;
                }
                let lv = left_r.iter().sum::<f64>() / left_r.len() as f64;
                let rv_val = right_r.iter().sum::<f64>() / right_r.len() as f64;
                let loss: f64 = left_r.iter().map(|&v| (v - lv).powi(2)).sum::<f64>()
                    + right_r.iter().map(|&v| (v - rv_val).powi(2)).sum::<f64>();
                if loss < best_loss {
                    best_loss = loss;
                    best_feat = feat;
                    best_thresh = thresh;
                    best_lv = lv;
                    best_rv = rv_val;
                }
            }
        }
        Self {
            feature: best_feat,
            threshold: best_thresh,
            left_value: best_lv,
            right_value: best_rv,
        }
    }
    /// Predict the residual for one sample.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        if x[self.feature] <= self.threshold {
            self.left_value
        } else {
            self.right_value
        }
    }
    /// Predict for a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
}
/// K-Nearest Neighbours classifier.
#[derive(Clone, Debug)]
pub struct KNearestNeighbours {
    /// Training data.
    pub x_train: Vec<Vec<f64>>,
    /// Training labels.
    pub y_train: Vec<f64>,
    /// Number of neighbours.
    pub k: usize,
}
impl KNearestNeighbours {
    /// Create a new KNN classifier.
    pub fn new(k: usize) -> Self {
        Self {
            x_train: Vec::new(),
            y_train: Vec::new(),
            k,
        }
    }
    /// Store training data.
    pub fn fit(&mut self, x: Vec<Vec<f64>>, y: Vec<f64>) {
        self.x_train = x;
        self.y_train = y;
    }
    /// Euclidean distance between two points.
    fn distance(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(x, y)| (x - y).powi(2))
            .sum::<f64>()
            .sqrt()
    }
    /// Predict for a single sample.
    pub fn predict_one(&self, x: &[f64]) -> f64 {
        let mut dists: Vec<(f64, f64)> = self
            .x_train
            .iter()
            .zip(self.y_train.iter())
            .map(|(xi, &yi)| (Self::distance(x, xi), yi))
            .collect();
        dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let neighbours: Vec<f64> = dists[..self.k.min(dists.len())]
            .iter()
            .map(|(_, l)| *l)
            .collect();
        majority_vote(&neighbours)
    }
    /// Predict for a batch.
    pub fn predict(&self, x: &[Vec<f64>]) -> Vec<f64> {
        x.iter().map(|row| self.predict_one(row)).collect()
    }
}
/// Adam optimiser.
#[derive(Clone, Debug)]
pub struct Adam {
    /// Learning rate.
    pub lr: f64,
    /// Beta1 (first moment decay).
    pub beta1: f64,
    /// Beta2 (second moment decay).
    pub beta2: f64,
    /// Epsilon for numerical stability.
    pub eps: f64,
    /// First moment estimate.
    pub m: Vec<f64>,
    /// Second moment estimate.
    pub v: Vec<f64>,
    /// Time step.
    pub t: usize,
}
impl Adam {
    /// Create Adam optimiser.
    pub fn new(n_params: usize, lr: f64, beta1: f64, beta2: f64, eps: f64) -> Self {
        Self {
            lr,
            beta1,
            beta2,
            eps,
            m: vec![0.0; n_params],
            v: vec![0.0; n_params],
            t: 0,
        }
    }
    /// Perform one Adam update step.
    pub fn step(&mut self, params: &mut Vec<f64>, grads: &[f64]) {
        self.t += 1;
        let t = self.t as f64;
        let bias_corr1 = 1.0 - self.beta1.powf(t);
        let bias_corr2 = 1.0 - self.beta2.powf(t);
        for (i, (p, &g)) in params.iter_mut().zip(grads.iter()).enumerate() {
            self.m[i] = self.beta1 * self.m[i] + (1.0 - self.beta1) * g;
            self.v[i] = self.beta2 * self.v[i] + (1.0 - self.beta2) * g * g;
            let m_hat = self.m[i] / bias_corr1;
            let v_hat = self.v[i] / bias_corr2;
            *p -= self.lr * m_hat / (v_hat.sqrt() + self.eps);
        }
    }
}
/// A single isolation tree node.
#[derive(Clone, Debug)]
pub enum IsolationNode {
    /// External leaf: size of partition when isolated.
    Leaf {
        /// Size of the leaf partition.
        size: usize,
    },
    /// Internal split.
    Internal {
        /// Feature index split on.
        feature: usize,
        /// Split threshold.
        threshold: f64,
        /// Left subtree.
        left: Box<IsolationNode>,
        /// Right subtree.
        right: Box<IsolationNode>,
    },
}
impl IsolationNode {
    /// Path length for a sample (number of edges to reach an external node).
    pub fn path_length(&self, x: &[f64], depth: f64) -> f64 {
        match self {
            IsolationNode::Leaf { size } => depth + average_path_length(*size),
            IsolationNode::Internal {
                feature,
                threshold,
                left,
                right,
            } => {
                if x[*feature] <= *threshold {
                    left.path_length(x, depth + 1.0)
                } else {
                    right.path_length(x, depth + 1.0)
                }
            }
        }
    }
}
/// Activation functions.
pub struct Activation;
impl Activation {
    /// ReLU activation.
    pub fn relu(x: &[f64]) -> Vec<f64> {
        x.iter().map(|&v| v.max(0.0)).collect()
    }
    /// ReLU derivative.
    pub fn relu_deriv(x: &[f64]) -> Vec<f64> {
        x.iter().map(|&v| if v > 0.0 { 1.0 } else { 0.0 }).collect()
    }
    /// Sigmoid activation.
    pub fn sigmoid(x: &[f64]) -> Vec<f64> {
        x.iter().map(|&v| sigmoid(v)).collect()
    }
    /// Sigmoid derivative.
    pub fn sigmoid_deriv(x: &[f64]) -> Vec<f64> {
        x.iter()
            .map(|&v| {
                let s = sigmoid(v);
                s * (1.0 - s)
            })
            .collect()
    }
    /// Tanh activation.
    pub fn tanh(x: &[f64]) -> Vec<f64> {
        x.iter().map(|&v| v.tanh()).collect()
    }
    /// Tanh derivative.
    pub fn tanh_deriv(x: &[f64]) -> Vec<f64> {
        x.iter().map(|&v| 1.0 - v.tanh().powi(2)).collect()
    }
    /// Softmax activation.
    pub fn softmax(x: &[f64]) -> Vec<f64> {
        let max = x.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exps: Vec<f64> = x.iter().map(|v| (v - max).exp()).collect();
        let sum: f64 = exps.iter().sum();
        exps.iter().map(|e| e / sum).collect()
    }
    /// Leaky ReLU.
    pub fn leaky_relu(x: &[f64], alpha: f64) -> Vec<f64> {
        x.iter()
            .map(|&v| if v > 0.0 { v } else { alpha * v })
            .collect()
    }
    /// ELU (Exponential Linear Unit).
    pub fn elu(x: &[f64], alpha: f64) -> Vec<f64> {
        x.iter()
            .map(|&v| if v >= 0.0 { v } else { alpha * (v.exp() - 1.0) })
            .collect()
    }
}
/// Stochastic Gradient Descent with optional momentum.
#[derive(Clone, Debug)]
pub struct Sgd {
    /// Learning rate.
    pub lr: f64,
    /// Momentum coefficient.
    pub momentum: f64,
    /// Velocity for each parameter.
    pub velocity: Vec<f64>,
}
impl Sgd {
    /// Create SGD optimiser.
    pub fn new(n_params: usize, lr: f64, momentum: f64) -> Self {
        Self {
            lr,
            momentum,
            velocity: vec![0.0; n_params],
        }
    }
    /// Compute parameter update given gradients.  Modifies params in-place.
    pub fn step(&mut self, params: &mut Vec<f64>, grads: &[f64]) {
        for (i, (p, &g)) in params.iter_mut().zip(grads.iter()).enumerate() {
            self.velocity[i] = self.momentum * self.velocity[i] - self.lr * g;
            *p += self.velocity[i];
        }
    }
}
