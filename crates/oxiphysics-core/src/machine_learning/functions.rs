//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::types::{
    CrossValResult, DecisionTree, IsolationNode, KMeansResult, LinearRegression, Matrix, OobResult,
    PcaResult, TreeNode, TsneResult,
};

/// Compute the dot product of two equal-length slices.
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
/// Compute element-wise sum of two vectors.
pub fn vec_add(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x + y).collect()
}
/// Compute element-wise difference of two vectors.
pub fn vec_sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(x, y)| x - y).collect()
}
/// Scale a vector by a scalar.
pub fn vec_scale(a: &[f64], s: f64) -> Vec<f64> {
    a.iter().map(|x| x * s).collect()
}
/// Compute the Euclidean norm of a vector.
pub fn vec_norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}
/// Compute the component-wise mean of a collection of vectors.
pub fn mean_vec(data: &[Vec<f64>]) -> Vec<f64> {
    if data.is_empty() {
        return Vec::new();
    }
    let d = data[0].len();
    let n = data.len() as f64;
    let mut m = vec![0.0f64; d];
    for row in data {
        for (j, &v) in row.iter().enumerate() {
            m[j] += v;
        }
    }
    m.iter().map(|x| x / n).collect()
}
/// Compute the component-wise variance of a collection of vectors.
pub fn variance_vec(data: &[Vec<f64>]) -> Vec<f64> {
    let m = mean_vec(data);
    let n = data.len() as f64;
    let d = m.len();
    let mut var = vec![0.0f64; d];
    for row in data {
        for j in 0..d {
            var[j] += (row[j] - m[j]).powi(2);
        }
    }
    var.iter().map(|x| x / n).collect()
}
/// Standardize data to zero mean and unit variance per component.
///
/// Returns `(standardized, mean, std_dev)`.
pub fn standardize(data: &[Vec<f64>]) -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
    let mu = mean_vec(data);
    let var = variance_vec(data);
    let std: Vec<f64> = var.iter().map(|v| v.sqrt().max(1e-300)).collect();
    let normed = data
        .iter()
        .map(|row| {
            row.iter()
                .zip(mu.iter())
                .zip(std.iter())
                .map(|((x, m), s)| (x - m) / s)
                .collect()
        })
        .collect();
    (normed, mu, std)
}
/// Solve Ax = b via Gaussian elimination with partial pivoting.
pub fn solve_linear_system(a: &Matrix, b: &[f64]) -> Option<Vec<f64>> {
    let n = a.rows;
    assert_eq!(a.cols, n);
    assert_eq!(b.len(), n);
    let mut m = a.clone();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let max_row = (col..n).max_by(|&i, &j| {
            m.get(i, col)
                .abs()
                .partial_cmp(&m.get(j, col).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })?;
        for c in 0..n {
            let tmp = m.get(col, c);
            m.set(col, c, m.get(max_row, c));
            m.set(max_row, c, tmp);
        }
        rhs.swap(col, max_row);
        let pivot = m.get(col, col);
        if pivot.abs() < 1e-15 {
            return None;
        }
        for row in col + 1..n {
            let factor = m.get(row, col) / pivot;
            for c in col..n {
                let val = m.get(row, c) - factor * m.get(col, c);
                m.set(row, c, val);
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for (j, &xj) in x.iter().enumerate().skip(i + 1) {
            sum -= m.get(i, j) * xj;
        }
        x[i] = sum / m.get(i, i);
    }
    Some(x)
}
/// Compute the sigmoid (logistic) activation function.
pub fn sigmoid(z: f64) -> f64 {
    1.0 / (1.0 + (-z).exp())
}
/// Apply the sigmoid function element-wise to a slice.
pub fn sigmoid_vec(z: &[f64]) -> Vec<f64> {
    z.iter().map(|&v| sigmoid(v)).collect()
}
/// Run k-means clustering with k-means++ initialisation.
pub fn kmeans(data: &[Vec<f64>], k: usize, max_iter: usize, tol: f64) -> KMeansResult {
    let n = data.len();
    let d = data[0].len();
    assert!(k <= n, "k must be <= number of data points");
    let mut centroids: Vec<Vec<f64>> = Vec::with_capacity(k);
    centroids.push(data[0].clone());
    for _ in 1..k {
        let dists: Vec<f64> = data
            .iter()
            .map(|p| {
                centroids
                    .iter()
                    .map(|c| {
                        p.iter()
                            .zip(c.iter())
                            .map(|(a, b)| (a - b).powi(2))
                            .sum::<f64>()
                    })
                    .fold(f64::MAX, f64::min)
            })
            .collect();
        let total: f64 = dists.iter().sum();
        let max_idx = dists
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal))
            .map(|(i, _)| i)
            .unwrap_or(0);
        let _ = total;
        centroids.push(data[max_idx].clone());
    }
    let mut labels = vec![0usize; n];
    let mut prev_inertia = f64::MAX;
    for _ in 0..max_iter {
        let mut changed = false;
        for (i, point) in data.iter().enumerate() {
            let nearest = centroids
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| {
                    let da: f64 = point
                        .iter()
                        .zip(a.iter())
                        .map(|(x, y)| (x - y).powi(2))
                        .sum();
                    let db: f64 = point
                        .iter()
                        .zip(b.iter())
                        .map(|(x, y)| (x - y).powi(2))
                        .sum();
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|(i, _)| i)
                .unwrap_or(0);
            if labels[i] != nearest {
                changed = true;
            }
            labels[i] = nearest;
        }
        for c in centroids.iter_mut() {
            *c = vec![0.0; d];
        }
        let mut counts = vec![0usize; k];
        for (i, point) in data.iter().enumerate() {
            let c = labels[i];
            for j in 0..d {
                centroids[c][j] += point[j];
            }
            counts[c] += 1;
        }
        for (c, centroid) in centroids.iter_mut().enumerate() {
            let cnt = counts[c].max(1) as f64;
            for v in centroid.iter_mut() {
                *v /= cnt;
            }
        }
        let inertia: f64 = data
            .iter()
            .zip(labels.iter())
            .map(|(p, &l)| {
                p.iter()
                    .zip(centroids[l].iter())
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
            })
            .sum();
        if !changed || (prev_inertia - inertia).abs() < tol {
            break;
        }
        prev_inertia = inertia;
    }
    let inertia: f64 = data
        .iter()
        .zip(labels.iter())
        .map(|(p, &l)| {
            p.iter()
                .zip(centroids[l].iter())
                .map(|(a, b)| (a - b).powi(2))
                .sum::<f64>()
        })
        .sum();
    KMeansResult {
        centroids,
        labels,
        inertia,
    }
}
/// Fit PCA using the covariance matrix and power iteration.
///
/// `n_components`: number of principal components to retain.
pub fn pca_fit(data: &[Vec<f64>], n_components: usize) -> PcaResult {
    let n = data.len() as f64;
    let d = data[0].len();
    let mean = mean_vec(data);
    let centred: Vec<Vec<f64>> = data
        .iter()
        .map(|row| row.iter().zip(mean.iter()).map(|(x, m)| x - m).collect())
        .collect();
    let mut cov = Matrix::zeros(d, d);
    for row in &centred {
        for i in 0..d {
            for j in 0..d {
                let val = cov.get(i, j) + row[i] * row[j];
                cov.set(i, j, val);
            }
        }
    }
    for v in cov.data.iter_mut() {
        *v /= n;
    }
    let mut components = Vec::new();
    let mut eigenvalues = Vec::new();
    let mut deflated = cov.clone();
    for _ in 0..n_components.min(d) {
        let (eigval, eigvec) = power_iteration(&deflated, 100, 1e-10);
        eigenvalues.push(eigval);
        components.push(eigvec.clone());
        for i in 0..d {
            for j in 0..d {
                let val = deflated.get(i, j) - eigval * eigvec[i] * eigvec[j];
                deflated.set(i, j, val);
            }
        }
    }
    PcaResult {
        components,
        eigenvalues,
        mean,
    }
}
/// Power iteration to find the dominant eigenpair of a symmetric matrix.
pub fn power_iteration(m: &Matrix, max_iter: usize, tol: f64) -> (f64, Vec<f64>) {
    let n = m.rows;
    let mut v: Vec<f64> = vec![1.0; n];
    let norm = vec_norm(&v);
    for x in v.iter_mut() {
        *x /= norm;
    }
    let mut eigenval = 0.0f64;
    for _ in 0..max_iter {
        let mv = m.matvec(&v);
        let new_eigenval = dot(&v, &mv);
        let norm_mv = vec_norm(&mv).max(1e-300);
        let new_v: Vec<f64> = mv.iter().map(|x| x / norm_mv).collect();
        let diff = vec_sub(&new_v, &v);
        v = new_v;
        if (new_eigenval - eigenval).abs() < tol && vec_norm(&diff) < tol {
            eigenval = new_eigenval;
            break;
        }
        eigenval = new_eigenval;
    }
    (eigenval, v)
}
/// Compute Gini impurity of a label set.
pub fn gini_impurity(labels: &[f64]) -> f64 {
    if labels.is_empty() {
        return 0.0;
    }
    let n = labels.len() as f64;
    let mut counts = std::collections::HashMap::new();
    for &l in labels {
        *counts.entry(l as i64).or_insert(0usize) += 1;
    }
    1.0 - counts
        .values()
        .map(|&c| (c as f64 / n).powi(2))
        .sum::<f64>()
}
/// Find the best split for a set of samples.
pub fn best_split(x: &[Vec<f64>], y: &[f64], min_samples: usize) -> Option<(usize, f64, f64)> {
    let n_features = x[0].len();
    let n = x.len();
    let mut best_gain = -f64::INFINITY;
    let mut best_feat = 0;
    let mut best_thresh = 0.0;
    let parent_gini = gini_impurity(y);
    for feat in 0..n_features {
        let mut values: Vec<f64> = x.iter().map(|r| r[feat]).collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        values.dedup();
        for thresh in values.windows(2).map(|w| (w[0] + w[1]) / 2.0) {
            let left_y: Vec<f64> = x
                .iter()
                .zip(y.iter())
                .filter(|(r, _)| r[feat] <= thresh)
                .map(|(_, &l)| l)
                .collect();
            let right_y: Vec<f64> = x
                .iter()
                .zip(y.iter())
                .filter(|(r, _)| r[feat] > thresh)
                .map(|(_, &l)| l)
                .collect();
            if left_y.len() < min_samples || right_y.len() < min_samples {
                continue;
            }
            let gain = parent_gini
                - (left_y.len() as f64 / n as f64) * gini_impurity(&left_y)
                - (right_y.len() as f64 / n as f64) * gini_impurity(&right_y);
            if gain > best_gain {
                best_gain = gain;
                best_feat = feat;
                best_thresh = thresh;
            }
        }
    }
    if best_gain > 0.0 {
        Some((best_feat, best_thresh, best_gain))
    } else {
        None
    }
}
/// Return the most frequent class label from a slice of labels.
pub fn majority_vote(y: &[f64]) -> f64 {
    let mut counts = std::collections::HashMap::new();
    for &l in y {
        *counts.entry(l as i64).or_insert(0usize) += 1;
    }
    counts
        .into_iter()
        .max_by_key(|&(_, c)| c)
        .map(|(l, _)| l as f64)
        .unwrap_or(0.0)
}
/// Recursively build a decision tree node.
pub fn build_tree(
    x: &[Vec<f64>],
    y: &[f64],
    depth: usize,
    max_depth: usize,
    min_samples: usize,
) -> TreeNode {
    if depth >= max_depth || y.len() < 2 * min_samples || gini_impurity(y) == 0.0 {
        return TreeNode::Leaf {
            value: majority_vote(y),
        };
    }
    if let Some((feat, thresh, _)) = best_split(x, y, min_samples) {
        let (lx, ly): (Vec<_>, Vec<_>) = x
            .iter()
            .zip(y.iter())
            .filter(|(r, _)| r[feat] <= thresh)
            .map(|(r, &l)| (r.clone(), l))
            .unzip();
        let (rx, ry): (Vec<_>, Vec<_>) = x
            .iter()
            .zip(y.iter())
            .filter(|(r, _)| r[feat] > thresh)
            .map(|(r, &l)| (r.clone(), l))
            .unzip();
        if lx.is_empty() || rx.is_empty() {
            return TreeNode::Leaf {
                value: majority_vote(y),
            };
        }
        TreeNode::Split {
            feature: feat,
            threshold: thresh,
            left: Box::new(build_tree(&lx, &ly, depth + 1, max_depth, min_samples)),
            right: Box::new(build_tree(&rx, &ry, depth + 1, max_depth, min_samples)),
        }
    } else {
        TreeNode::Leaf {
            value: majority_vote(y),
        }
    }
}
/// Run k-fold cross-validation for a linear regression model (R² metric).
pub fn cross_val_linear_regression(x: &[Vec<f64>], y: &[f64], k: usize) -> CrossValResult {
    let n = x.len();
    let fold_size = n / k;
    let mut scores = Vec::new();
    for fold in 0..k {
        let val_start = fold * fold_size;
        let val_end = ((fold + 1) * fold_size).min(n);
        let train_x: Vec<Vec<f64>> = x[..val_start]
            .iter()
            .chain(x[val_end..].iter())
            .cloned()
            .collect();
        let train_y: Vec<f64> = y[..val_start]
            .iter()
            .chain(y[val_end..].iter())
            .cloned()
            .collect();
        let val_x = &x[val_start..val_end];
        let val_y = &y[val_start..val_end];
        if train_x.is_empty() || val_x.is_empty() {
            continue;
        }
        let model = LinearRegression::fit_normal(&train_x, &train_y);
        scores.push(model.r2_score(val_x, val_y));
    }
    let mean = scores.iter().sum::<f64>() / scores.len() as f64;
    let std = (scores.iter().map(|s| (s - mean).powi(2)).sum::<f64>() / scores.len() as f64).sqrt();
    CrossValResult {
        fold_scores: scores,
        mean_score: mean,
        std_score: std,
    }
}
/// Compute the confusion matrix for binary or multi-class classification.
///
/// Returns a 2-D vector `cm[true][pred]`.
pub fn confusion_matrix(y_true: &[usize], y_pred: &[usize], n_classes: usize) -> Vec<Vec<usize>> {
    let mut cm = vec![vec![0usize; n_classes]; n_classes];
    for (&t, &p) in y_true.iter().zip(y_pred.iter()) {
        if t < n_classes && p < n_classes {
            cm[t][p] += 1;
        }
    }
    cm
}
/// Compute precision, recall, and F1 score from a binary confusion matrix.
pub fn precision_recall_f1(cm: &[Vec<usize>]) -> (f64, f64, f64) {
    if cm.len() < 2 {
        return (0.0, 0.0, 0.0);
    }
    let tp = cm[1][1] as f64;
    let fp = cm[0][1] as f64;
    let fn_ = cm[1][0] as f64;
    let precision = if tp + fp > 0.0 { tp / (tp + fp) } else { 0.0 };
    let recall = if tp + fn_ > 0.0 { tp / (tp + fn_) } else { 0.0 };
    let f1 = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
    } else {
        0.0
    };
    (precision, recall, f1)
}
/// Compute overall accuracy from a confusion matrix.
pub fn accuracy_from_cm(cm: &[Vec<usize>]) -> f64 {
    let correct: usize = (0..cm.len()).map(|i| cm[i][i]).sum();
    let total: usize = cm.iter().flat_map(|r| r.iter()).sum();
    if total == 0 {
        0.0
    } else {
        correct as f64 / total as f64
    }
}
/// Compute the ROC AUC (Area Under the Curve) for binary classification.
///
/// `scores`: predicted probabilities. `labels`: true binary labels.
pub fn roc_auc(scores: &[f64], labels: &[usize]) -> f64 {
    let n = scores.len();
    let mut pairs: Vec<(f64, usize)> = scores.iter().cloned().zip(labels.iter().cloned()).collect();
    pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let n_pos = labels.iter().filter(|&&l| l == 1).count() as f64;
    let n_neg = (n - n_pos as usize) as f64;
    if n_pos == 0.0 || n_neg == 0.0 {
        return 0.5;
    }
    let mut auc = 0.0;
    let mut fp = 0.0;
    let mut tp_prev = 0.0;
    let mut fp_prev = 0.0;
    for &(_, label) in &pairs {
        if label == 1 {
            tp_prev += 1.0;
        } else {
            fp += 1.0;
        }
        auc += (fp / n_neg - fp_prev / n_neg) * tp_prev / n_pos;
        fp_prev = fp;
    }
    auc
}
/// Mean absolute error.
pub fn mae(y_true: &[f64], y_pred: &[f64]) -> f64 {
    y_true
        .iter()
        .zip(y_pred.iter())
        .map(|(t, p)| (t - p).abs())
        .sum::<f64>()
        / y_true.len() as f64
}
/// Mean squared error.
pub fn mse(y_true: &[f64], y_pred: &[f64]) -> f64 {
    y_true
        .iter()
        .zip(y_pred.iter())
        .map(|(t, p)| (t - p).powi(2))
        .sum::<f64>()
        / y_true.len() as f64
}
/// Root mean squared error.
pub fn rmse(y_true: &[f64], y_pred: &[f64]) -> f64 {
    mse(y_true, y_pred).sqrt()
}
/// R² (coefficient of determination).
pub fn r2(y_true: &[f64], y_pred: &[f64]) -> f64 {
    let mean_y = y_true.iter().sum::<f64>() / y_true.len() as f64;
    let ss_res: f64 = y_true
        .iter()
        .zip(y_pred.iter())
        .map(|(t, p)| (t - p).powi(2))
        .sum();
    let ss_tot: f64 = y_true.iter().map(|t| (t - mean_y).powi(2)).sum();
    if ss_tot == 0.0 {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::machine_learning::Activation;
    use crate::machine_learning::Adam;
    use crate::machine_learning::BatchNorm;
    use crate::machine_learning::DenseLayer;
    use crate::machine_learning::GaussianNaiveBayes;

    use crate::machine_learning::KNearestNeighbours;

    use crate::machine_learning::LogisticRegression;
    use crate::machine_learning::MinMaxScaler;
    use crate::machine_learning::RandomForest;

    use crate::machine_learning::RmsProp;
    use crate::machine_learning::StandardScaler;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    #[test]
    fn test_matrix_zeros() {
        let m = Matrix::zeros(2, 3);
        assert_eq!(m.rows, 2);
        assert_eq!(m.cols, 3);
        assert!(m.data.iter().all(|&v| v == 0.0));
    }
    #[test]
    fn test_matrix_eye() {
        let m = Matrix::eye(3);
        assert_eq!(m.get(0, 0), 1.0);
        assert_eq!(m.get(0, 1), 0.0);
        assert_eq!(m.get(1, 1), 1.0);
    }
    #[test]
    fn test_matrix_matmul() {
        let a = Matrix::from_vec(2, 2, vec![1.0, 2.0, 3.0, 4.0]);
        let b = Matrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, 1.0]);
        let c = a.matmul(&b);
        assert_eq!(c.data, a.data);
    }
    #[test]
    fn test_matrix_transpose() {
        let a = Matrix::from_vec(2, 3, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let at = a.transpose();
        assert_eq!(at.rows, 3);
        assert_eq!(at.cols, 2);
        assert_eq!(at.get(1, 0), 2.0);
    }
    #[test]
    fn test_matrix_matvec() {
        let a = Matrix::from_vec(2, 2, vec![1.0, 0.0, 0.0, 2.0]);
        let v = vec![3.0, 4.0];
        let out = a.matvec(&v);
        assert_eq!(out, vec![3.0, 8.0]);
    }
    #[test]
    fn test_linear_regression_perfect_fit() {
        let x: Vec<Vec<f64>> = (0..5).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = (0..5).map(|i| 2.0 * i as f64 + 1.0).collect();
        let model = LinearRegression::fit_normal(&x, &y);
        assert!(approx_eq(model.r2_score(&x, &y), 1.0, 1e-8));
    }
    #[test]
    fn test_linear_regression_gradient_descent() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64 / 10.0]).collect();
        let y: Vec<f64> = x.iter().map(|r| 3.0 * r[0] + 1.0).collect();
        let model = LinearRegression::fit_gradient_descent(&x, &y, 0.05, 2000);
        let r2 = model.r2_score(&x, &y);
        assert!(r2 > 0.95);
    }
    #[test]
    fn test_linear_regression_predict() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0]];
        let y = vec![0.0, 1.0, 2.0];
        let model = LinearRegression::fit_normal(&x, &y);
        let p = model.predict_one(&[3.0]);
        assert!(approx_eq(p, 3.0, 0.01));
    }
    #[test]
    fn test_logistic_regression_linearly_separable() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64 - 10.0]).collect();
        let y: Vec<f64> = (0..20).map(|i| if i < 10 { 0.0 } else { 1.0 }).collect();
        let model = LogisticRegression::fit(&x, &y, 0.5, 1000);
        let acc = model.accuracy(&x, &y);
        assert!(acc >= 0.9);
    }
    #[test]
    fn test_logistic_regression_sigmoid() {
        assert!(approx_eq(sigmoid(0.0), 0.5, 1e-10));
        assert!(sigmoid(10.0) > 0.99);
        assert!(sigmoid(-10.0) < 0.01);
    }
    #[test]
    fn test_kmeans_two_clusters() {
        let mut data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();
        data.extend((100..110).map(|i| vec![i as f64]));
        let result = kmeans(&data, 2, 100, 1e-6);
        assert_eq!(result.centroids.len(), 2);
        assert_eq!(result.labels.len(), 20);
        let c0 = result.labels[0];
        for &l in &result.labels[..10] {
            assert_eq!(l, c0);
        }
    }
    #[test]
    fn test_kmeans_inertia_positive() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, i as f64]).collect();
        let result = kmeans(&data, 3, 50, 1e-6);
        assert!(result.inertia >= 0.0);
    }
    #[test]
    fn test_pca_fit_shape() {
        let data: Vec<Vec<f64>> = (0..20)
            .map(|i| vec![i as f64, 2.0 * i as f64, 3.0 * i as f64])
            .collect();
        let pca = pca_fit(&data, 2);
        assert_eq!(pca.components.len(), 2);
        assert_eq!(pca.components[0].len(), 3);
    }
    #[test]
    fn test_pca_transform_shape() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, i as f64 * 2.0]).collect();
        let pca = pca_fit(&data, 1);
        let proj = pca.transform(&data);
        assert_eq!(proj.len(), 10);
        assert_eq!(proj[0].len(), 1);
    }
    #[test]
    fn test_relu() {
        let x = vec![-1.0, 0.0, 2.0, -3.0, 4.0];
        let y = Activation::relu(&x);
        assert_eq!(y, vec![0.0, 0.0, 2.0, 0.0, 4.0]);
    }
    #[test]
    fn test_softmax_sums_to_one() {
        let x = vec![1.0, 2.0, 3.0];
        let s = Activation::softmax(&x);
        let sum: f64 = s.iter().sum();
        assert!(approx_eq(sum, 1.0, 1e-10));
    }
    #[test]
    fn test_sigmoid_symmetry() {
        assert!(approx_eq(sigmoid(1.0) + sigmoid(-1.0), 1.0, 1e-10));
    }
    #[test]
    fn test_dense_layer_forward_shape() {
        let mut layer = DenseLayer::new(4, 3);
        let x = vec![1.0, 0.0, -1.0, 0.5];
        let out = layer.forward(&x);
        assert_eq!(out.len(), 3);
    }
    #[test]
    fn test_dense_layer_backward_shape() {
        let mut layer = DenseLayer::new(3, 2);
        let x = vec![1.0, 2.0, 3.0];
        let _ = layer.forward(&x);
        let grad_out = vec![0.1, -0.2];
        let grad_in = layer.backward(&grad_out, 0.01);
        assert_eq!(grad_in.len(), 3);
    }
    #[test]
    fn test_batch_norm_output_shape() {
        let mut bn = BatchNorm::new(3);
        let batch = vec![vec![1.0, 2.0, 3.0], vec![4.0, 5.0, 6.0]];
        let out = bn.forward_batch(&batch);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].len(), 3);
    }
    #[test]
    fn test_adam_step_decreases_params() {
        let mut params = vec![1.0, 1.0];
        let grads = vec![0.5, 0.5];
        let mut opt = Adam::new(2, 0.01, 0.9, 0.999, 1e-8);
        opt.step(&mut params, &grads);
        assert!(params[0] < 1.0);
    }
    #[test]
    fn test_rmsprop_step() {
        let mut params = vec![1.0];
        let grads = vec![1.0];
        let mut opt = RmsProp::new(1, 0.01, 0.9, 1e-8);
        opt.step(&mut params, &grads);
        assert!(params[0] < 1.0);
    }
    #[test]
    fn test_decision_tree_xor() {
        let x = vec![
            vec![0.0, 0.0],
            vec![0.5, 0.0],
            vec![0.0, 0.5],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![0.5, 1.0],
        ];
        let y = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
        let mut tree = DecisionTree::new(4, 1);
        tree.fit(&x, &y);
        let acc = tree.accuracy(&x, &y);
        assert!(acc >= 0.75);
    }
    #[test]
    fn test_decision_tree_linearly_separable() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = (0..20).map(|i| if i < 10 { 0.0 } else { 1.0 }).collect();
        let mut tree = DecisionTree::new(3, 1);
        tree.fit(&x, &y);
        let acc = tree.accuracy(&x, &y);
        assert!(acc >= 0.9);
    }
    #[test]
    fn test_random_forest_accuracy() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64, (i % 2) as f64]).collect();
        let y: Vec<f64> = (0..20).map(|i| if i < 10 { 0.0 } else { 1.0 }).collect();
        let mut rf = RandomForest::new(5, 3, 1);
        rf.fit(&x, &y);
        let acc = rf.accuracy(&x, &y);
        assert!(acc >= 0.8);
    }
    #[test]
    fn test_confusion_matrix() {
        let y_true = vec![0, 1, 1, 0, 1];
        let y_pred = vec![0, 1, 0, 0, 1];
        let cm = confusion_matrix(&y_true, &y_pred, 2);
        assert_eq!(cm[1][1], 2);
        assert_eq!(cm[1][0], 1);
    }
    #[test]
    fn test_precision_recall_f1() {
        let cm = vec![vec![2usize, 0], vec![0, 3]];
        let (p, r, f1) = precision_recall_f1(&cm);
        assert!(approx_eq(p, 1.0, 1e-10));
        assert!(approx_eq(r, 1.0, 1e-10));
        assert!(approx_eq(f1, 1.0, 1e-10));
    }
    #[test]
    fn test_accuracy_from_cm() {
        let cm = vec![vec![3usize, 1], vec![2, 4]];
        let acc = accuracy_from_cm(&cm);
        assert!(approx_eq(acc, 7.0 / 10.0, 1e-10));
    }
    #[test]
    fn test_mae() {
        let y = vec![1.0, 2.0, 3.0];
        let p = vec![1.0, 2.0, 4.0];
        assert!(approx_eq(mae(&y, &p), 1.0 / 3.0, 1e-10));
    }
    #[test]
    fn test_rmse() {
        let y = vec![0.0, 0.0];
        let p = vec![1.0, 1.0];
        assert!(approx_eq(rmse(&y, &p), 1.0, 1e-10));
    }
    #[test]
    fn test_r2_perfect() {
        let y = vec![1.0, 2.0, 3.0, 4.0];
        assert!(approx_eq(r2(&y, &y), 1.0, 1e-10));
    }
    #[test]
    fn test_min_max_scaler() {
        let data = vec![vec![0.0, 10.0], vec![2.0, 20.0], vec![4.0, 30.0]];
        let scaler = MinMaxScaler::fit(&data);
        let t = scaler.transform(&data);
        assert!(approx_eq(t[0][0], 0.0, 1e-10));
        assert!(approx_eq(t[2][0], 1.0, 1e-10));
        let inv = scaler.inverse_transform(&t);
        assert!(approx_eq(inv[1][1], 20.0, 1e-8));
    }
    #[test]
    fn test_standard_scaler() {
        let data = vec![vec![2.0], vec![4.0], vec![6.0]];
        let scaler = StandardScaler::fit(&data);
        let t = scaler.transform(&data);
        let mean: f64 = t.iter().map(|r| r[0]).sum::<f64>() / 3.0;
        assert!(mean.abs() < 1e-8);
    }
    #[test]
    fn test_gnb_fit_predict() {
        let x = vec![vec![1.0], vec![2.0], vec![10.0], vec![11.0]];
        let y = vec![0, 0, 1, 1];
        let gnb = GaussianNaiveBayes::fit(&x, &y, 2);
        assert_eq!(gnb.predict_one(&[1.5]), 0);
        assert_eq!(gnb.predict_one(&[10.5]), 1);
    }
    #[test]
    fn test_knn_predict() {
        let x_train = vec![vec![0.0], vec![1.0], vec![10.0], vec![11.0]];
        let y_train = vec![0.0, 0.0, 1.0, 1.0];
        let mut knn = KNearestNeighbours::new(2);
        knn.fit(x_train, y_train);
        assert!(approx_eq(knn.predict_one(&[0.5]), 0.0, 0.5));
        assert!(approx_eq(knn.predict_one(&[10.5]), 1.0, 0.5));
    }
    #[test]
    fn test_cross_val_linear() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = (0..20).map(|i| 2.0 * i as f64).collect();
        let result = cross_val_linear_regression(&x, &y, 4);
        assert!(result.mean_score > 0.9);
    }
    #[test]
    fn test_solve_linear_system() {
        let a = Matrix::from_vec(2, 2, vec![2.0, 1.0, 1.0, 3.0]);
        let b = vec![5.0, 10.0];
        let x = solve_linear_system(&a, &b).unwrap();
        assert!(approx_eq(x[0], 1.0, 1e-8));
        assert!(approx_eq(x[1], 3.0, 1e-8));
    }
    #[test]
    fn test_roc_auc_perfect() {
        let scores = vec![0.9, 0.8, 0.2, 0.1];
        let labels = vec![1, 1, 0, 0];
        let auc = roc_auc(&scores, &labels);
        assert!(approx_eq(auc, 1.0, 1e-8));
    }
    #[test]
    fn test_roc_auc_random() {
        let scores = vec![0.5, 0.5, 0.5, 0.5];
        let labels = vec![1, 0, 1, 0];
        let auc = roc_auc(&scores, &labels);
        assert!((0.0..=1.0).contains(&auc));
    }
}
/// Compute the out-of-bag error for a random forest.
///
/// Returns predictions for each sample using only trees that did NOT
/// see that sample during training (bootstrap).
pub fn random_forest_oob(x: &[Vec<f64>], y: &[f64], n_trees: usize, max_depth: usize) -> OobResult {
    let n = x.len();
    let mut all_votes: Vec<Vec<f64>> = vec![Vec::new(); n];
    let mut bootstrap_sets = Vec::with_capacity(n_trees);
    for t in 0..n_trees {
        let mut seed = (t as u64 + 17).wrapping_mul(6364136223846793005);
        let mut indices = Vec::with_capacity(n);
        for _ in 0..n {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            indices.push((seed >> 33) as usize % n);
        }
        let in_bag: std::collections::HashSet<usize> = indices.iter().cloned().collect();
        let bx: Vec<Vec<f64>> = indices.iter().map(|&i| x[i].clone()).collect();
        let by: Vec<f64> = indices.iter().map(|&i| y[i]).collect();
        let mut tree = DecisionTree::new(max_depth, 1);
        tree.fit(&bx, &by);
        for i in 0..n {
            if !in_bag.contains(&i) {
                all_votes[i].push(tree.predict_one(&x[i]));
            }
        }
        bootstrap_sets.push(indices);
    }
    let oob_predictions: Vec<f64> = all_votes
        .iter()
        .map(|votes| {
            if votes.is_empty() {
                0.0
            } else {
                majority_vote(votes)
            }
        })
        .collect();
    let correct = oob_predictions
        .iter()
        .zip(y.iter())
        .filter(|(p, t)| (*p - *t).abs() < 0.5)
        .count();
    let oob_accuracy = correct as f64 / n as f64;
    OobResult {
        oob_predictions,
        oob_accuracy,
        bootstrap_sets,
    }
}
/// Select a random subset of `k` feature indices from `0..d`.
///
/// Uses a deterministic LCG seeded by `seed`.
pub fn feature_bagging(d: usize, k: usize, seed: u64) -> Vec<usize> {
    let k = k.min(d);
    let mut rng_state = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    let mut indices: Vec<usize> = (0..d).collect();
    for i in 0..k {
        rng_state = rng_state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let j = i + (rng_state >> 33) as usize % (d - i);
        indices.swap(i, j);
    }
    indices[..k].to_vec()
}
/// RBF (squared-exponential) kernel.
///
/// k(x, x') = amplitude² · exp(-||x - x'||² / (2·length_scale²))
pub fn rbf_kernel(x1: &[f64], x2: &[f64], amplitude: f64, length_scale: f64) -> f64 {
    let sq_dist: f64 = x1.iter().zip(x2.iter()).map(|(a, b)| (a - b).powi(2)).sum();
    amplitude * amplitude * (-sq_dist / (2.0 * length_scale * length_scale)).exp()
}
/// Compute the kernel matrix K\[i, j\] = k(X\[i\], X\[j\]).
pub fn kernel_matrix(x: &[Vec<f64>], amplitude: f64, length_scale: f64) -> Matrix {
    let n = x.len();
    let mut km = Matrix::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            km.set(i, j, rbf_kernel(&x[i], &x[j], amplitude, length_scale));
        }
    }
    km
}
/// Cholesky decomposition (lower triangular factor).
///
/// Returns L such that A = L Lᵀ.  Assumes A is symmetric positive-definite.
pub fn cholesky(a: &Matrix) -> Matrix {
    let n = a.rows;
    let mut l = Matrix::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let mut sum = a.get(i, j);
            for k in 0..j {
                sum -= l.get(i, k) * l.get(j, k);
            }
            if i == j {
                l.set(i, j, (sum.max(1e-300)).sqrt());
            } else {
                let ljj = l.get(j, j).max(1e-300);
                l.set(i, j, sum / ljj);
            }
        }
    }
    l
}
/// Forward substitution: solve L x = b where L is lower-triangular.
pub fn forward_substitution(l: &Matrix, b: &[f64]) -> Vec<f64> {
    let n = l.rows;
    let mut x = vec![0.0f64; n];
    for i in 0..n {
        let mut s = b[i];
        for (j, &xj) in x.iter().enumerate().take(i) {
            s -= l.get(i, j) * xj;
        }
        x[i] = s / l.get(i, i).max(1e-300);
    }
    x
}
/// Backward substitution against Lᵀ: solve Lᵀ x = b.
pub fn backward_substitution_t(l: &Matrix, b: &[f64]) -> Vec<f64> {
    let n = l.rows;
    let mut x = vec![0.0f64; n];
    for i in (0..n).rev() {
        let mut s = b[i];
        for (j, &xj) in x.iter().enumerate().skip(i + 1) {
            s -= l.get(j, i) * xj;
        }
        x[i] = s / l.get(i, i).max(1e-300);
    }
    x
}
/// Expected path length of a random BST with `n` nodes.
///
/// Formula: c(n) = 2 H(n-1) - 2(n-1)/n  where H is the harmonic number.
pub fn average_path_length(n: usize) -> f64 {
    if n <= 1 {
        return 0.0;
    }
    let nf = n as f64;
    let h_n = (nf - 1.0).ln() + 0.5772156649;
    2.0 * h_n - 2.0 * (nf - 1.0) / nf
}
/// Build one isolation tree from a subsample.
pub fn build_isolation_tree(
    x: &[Vec<f64>],
    max_depth: usize,
    depth: usize,
    seed: &mut u64,
) -> IsolationNode {
    let n = x.len();
    if n <= 1 || depth >= max_depth {
        return IsolationNode::Leaf { size: n };
    }
    let d = x[0].len();
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let feat = (*seed >> 33) as usize % d;
    let mut min_v = x[0][feat];
    let mut max_v = min_v;
    for row in x.iter() {
        if row[feat] < min_v {
            min_v = row[feat];
        }
        if row[feat] > max_v {
            max_v = row[feat];
        }
    }
    if (max_v - min_v).abs() < 1e-15 {
        return IsolationNode::Leaf { size: n };
    }
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    let frac = (*seed >> 11) as f64 / (1u64 << 53) as f64;
    let threshold = min_v + frac * (max_v - min_v);
    let (left_x, right_x): (Vec<Vec<f64>>, Vec<Vec<f64>>) =
        x.iter().cloned().partition(|row| row[feat] <= threshold);
    IsolationNode::Internal {
        feature: feat,
        threshold,
        left: Box::new(build_isolation_tree(&left_x, max_depth, depth + 1, seed)),
        right: Box::new(build_isolation_tree(&right_x, max_depth, depth + 1, seed)),
    }
}
/// Compute the reachability distance: max(k-dist(o), dist(p, o)).
pub fn reach_dist(p: &[f64], o: &[f64], k_dist_o: f64) -> f64 {
    let d = p
        .iter()
        .zip(o.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt();
    d.max(k_dist_o)
}
/// Compute LOF scores for all points in `x`.
///
/// LOF ≈ 1 means typical point; LOF >> 1 means outlier.
pub fn lof_scores(x: &[Vec<f64>], k: usize) -> Vec<f64> {
    let n = x.len();
    if n == 0 {
        return Vec::new();
    }
    let dist = |a: &[f64], b: &[f64]| -> f64 {
        a.iter()
            .zip(b.iter())
            .map(|(u, v)| (u - v).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    let mut k_dists = vec![0.0f64; n];
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        let mut dists: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (dist(&x[i], &x[j]), j))
            .collect();
        dists.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let kk = k.min(dists.len());
        k_dists[i] = dists.get(kk.saturating_sub(1)).map(|d| d.0).unwrap_or(0.0);
        neighbours[i] = dists[..kk].iter().map(|d| d.1).collect();
    }
    let mut lrd = vec![0.0f64; n];
    for i in 0..n {
        let nb = &neighbours[i];
        if nb.is_empty() {
            lrd[i] = 1.0;
            continue;
        }
        let sum_reach: f64 = nb
            .iter()
            .map(|&j| reach_dist(&x[i], &x[j], k_dists[j]))
            .sum();
        lrd[i] = nb.len() as f64 / sum_reach.max(1e-300);
    }
    (0..n)
        .map(|i| {
            let nb = &neighbours[i];
            if nb.is_empty() {
                return 1.0;
            }
            let avg_lrd: f64 = nb.iter().map(|&j| lrd[j]).sum::<f64>() / nb.len() as f64;
            avg_lrd / lrd[i].max(1e-300)
        })
        .collect()
}
/// Compute pairwise affinities p_{ij} (symmetric, normalised) for t-SNE.
///
/// Uses a Gaussian kernel with a fixed perplexity converted to a bandwidth σ.
/// Bandwidth is estimated via binary search on the entropy condition.
pub fn compute_pairwise_affinities(x: &[Vec<f64>], perplexity: f64) -> Vec<Vec<f64>> {
    let n = x.len();
    let target_entropy = perplexity.ln();
    let mut p_matrix = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        let mut beta_lo = 0.0f64;
        let mut beta_hi = f64::INFINITY;
        let mut beta = 1.0f64;
        let dists: Vec<f64> = (0..n)
            .map(|j| {
                if j == i {
                    0.0
                } else {
                    x[i].iter()
                        .zip(x[j].iter())
                        .map(|(a, b)| (a - b).powi(2))
                        .sum()
                }
            })
            .collect();
        for _ in 0..50 {
            let exps: Vec<f64> = dists
                .iter()
                .enumerate()
                .map(|(j, &d)| if j == i { 0.0 } else { (-beta * d).exp() })
                .collect();
            let sum_exps: f64 = exps.iter().sum::<f64>().max(1e-300);
            let probs: Vec<f64> = exps.iter().map(|&e| e / sum_exps).collect();
            let h: f64 = probs
                .iter()
                .zip(dists.iter())
                .filter(|(p, _)| **p > 1e-300)
                .map(|(p, d)| -p * (p.ln()) + beta * p * d)
                .sum();
            if (h - target_entropy).abs() < 1e-5 {
                p_matrix[i][..n].copy_from_slice(&probs[..n]);
                break;
            }
            if h > target_entropy {
                beta_lo = beta;
                beta = if beta_hi.is_infinite() {
                    beta * 2.0
                } else {
                    (beta + beta_hi) / 2.0
                };
            } else {
                beta_hi = beta;
                beta = (beta + beta_lo) / 2.0;
            }
        }
        let exps: Vec<f64> = dists
            .iter()
            .enumerate()
            .map(|(j, &d)| if j == i { 0.0 } else { (-beta * d).exp() })
            .collect();
        let sum_exps: f64 = exps.iter().sum::<f64>().max(1e-300);
        for j in 0..n {
            p_matrix[i][j] = exps[j] / sum_exps;
        }
    }
    let total_n = (n * n) as f64;
    let mut p_sym = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            p_sym[i][j] = (p_matrix[i][j] + p_matrix[j][i]) / total_n;
        }
    }
    p_sym
}
/// Run a simplified t-SNE reduction to 2 dimensions.
///
/// Uses gradient descent with momentum.  Suitable for small datasets (n ≤ 500).
///
/// # Arguments
/// * `x` — input data (n × d)
/// * `perplexity` — typical range 5–50
/// * `n_iter` — number of gradient descent steps
/// * `lr` — learning rate (typical: 100–500)
pub fn tsne(x: &[Vec<f64>], perplexity: f64, n_iter: usize, lr: f64) -> TsneResult {
    let n = x.len();
    if n == 0 {
        return TsneResult {
            embedding: Vec::new(),
            kl_divergence: 0.0,
        };
    }
    let p = compute_pairwise_affinities(x, perplexity);
    let mut y: Vec<[f64; 2]> = (0..n)
        .map(|i| {
            let v = (i as f64 * 1.3).sin() * 0.01;
            let w = (i as f64 * 2.7).cos() * 0.01;
            [v, w]
        })
        .collect();
    let mut gains = vec![[1.0f64; 2]; n];
    let mut velocity = vec![[0.0f64; 2]; n];
    let momentum = 0.8;
    for _iter in 0..n_iter {
        let mut q_num = vec![vec![0.0f64; n]; n];
        let mut z = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let diff = [y[i][0] - y[j][0], y[i][1] - y[j][1]];
                let qval = 1.0 / (1.0 + diff[0] * diff[0] + diff[1] * diff[1]);
                q_num[i][j] = qval;
                q_num[j][i] = qval;
                z += 2.0 * qval;
            }
        }
        let z = z.max(1e-300);
        let mut grad = vec![[0.0f64; 2]; n];
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let q = q_num[i][j] / z;
                let factor = 4.0 * (p[i][j] - q) * q_num[i][j];
                let diff = [y[i][0] - y[j][0], y[i][1] - y[j][1]];
                grad[i][0] += factor * diff[0];
                grad[i][1] += factor * diff[1];
            }
        }
        for i in 0..n {
            for d in 0..2 {
                gains[i][d] = if (grad[i][d] > 0.0) == (velocity[i][d] > 0.0) {
                    (gains[i][d] * 0.8).max(0.01)
                } else {
                    gains[i][d] + 0.2
                };
                velocity[i][d] = momentum * velocity[i][d] - lr * gains[i][d] * grad[i][d];
                y[i][d] += velocity[i][d];
            }
        }
    }
    let mut q_num = vec![vec![0.0f64; n]; n];
    let mut z = 0.0f64;
    for i in 0..n {
        for j in (i + 1)..n {
            let diff = [y[i][0] - y[j][0], y[i][1] - y[j][1]];
            let qval = 1.0 / (1.0 + diff[0] * diff[0] + diff[1] * diff[1]);
            q_num[i][j] = qval;
            q_num[j][i] = qval;
            z += 2.0 * qval;
        }
    }
    let z = z.max(1e-300);
    let kl: f64 = (0..n)
        .flat_map(|i| (0..n).map(move |j| (i, j)))
        .filter(|&(i, j)| i != j && p[i][j] > 1e-300)
        .map(|(i, j)| {
            let q = (q_num[i][j] / z).max(1e-300);
            p[i][j] * (p[i][j] / q).ln()
        })
        .sum();
    TsneResult {
        embedding: y,
        kl_divergence: kl,
    }
}
/// Per-class precision, recall, and F1 for multi-class problems (one-vs-rest).
///
/// Returns a Vec of `(class, precision, recall, f1)` tuples.
pub fn multiclass_prf(cm: &[Vec<usize>]) -> Vec<(usize, f64, f64, f64)> {
    let n_classes = cm.len();
    (0..n_classes)
        .map(|c| {
            let tp = cm[c][c] as f64;
            let fp: f64 = (0..n_classes)
                .filter(|&r| r != c)
                .map(|r| cm[r][c] as f64)
                .sum();
            let fn_: f64 = (0..n_classes)
                .filter(|&col| col != c)
                .map(|col| cm[c][col] as f64)
                .sum();
            let prec = if tp + fp > 0.0 { tp / (tp + fp) } else { 0.0 };
            let rec = if tp + fn_ > 0.0 { tp / (tp + fn_) } else { 0.0 };
            let f1 = if prec + rec > 0.0 {
                2.0 * prec * rec / (prec + rec)
            } else {
                0.0
            };
            (c, prec, rec, f1)
        })
        .collect()
}
/// Macro-averaged F1 across all classes.
pub fn macro_f1(cm: &[Vec<usize>]) -> f64 {
    let prf = multiclass_prf(cm);
    if prf.is_empty() {
        return 0.0;
    }
    prf.iter().map(|(_, _, _, f1)| f1).sum::<f64>() / prf.len() as f64
}
/// Weighted-averaged F1 (weighted by class support).
pub fn weighted_f1(cm: &[Vec<usize>]) -> f64 {
    let prf = multiclass_prf(cm);
    let n_classes = cm.len();
    let supports: Vec<f64> = (0..n_classes)
        .map(|c| cm[c].iter().sum::<usize>() as f64)
        .collect();
    let total: f64 = supports.iter().sum();
    if total == 0.0 {
        return 0.0;
    }
    prf.iter()
        .zip(supports.iter())
        .map(|((_, _, _, f1), &s)| f1 * s)
        .sum::<f64>()
        / total
}
/// Compute the full ROC curve: returns `(fpr, tpr)` sorted by threshold.
pub fn roc_curve(scores: &[f64], labels: &[usize]) -> (Vec<f64>, Vec<f64>) {
    let n_pos = labels.iter().filter(|&&l| l == 1).count() as f64;
    let n_neg = (labels.len() as f64 - n_pos).max(1.0);
    let mut pairs: Vec<(f64, usize)> = scores.iter().cloned().zip(labels.iter().cloned()).collect();
    pairs.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut tpr = vec![0.0f64];
    let mut fpr = vec![0.0f64];
    let mut tp = 0.0f64;
    let mut fp = 0.0f64;
    for (_, label) in &pairs {
        if *label == 1 {
            tp += 1.0;
        } else {
            fp += 1.0;
        }
        tpr.push(tp / n_pos);
        fpr.push(fp / n_neg);
    }
    (fpr, tpr)
}
/// Compute AUC using the trapezoidal rule.
pub fn auc_trapezoid(x: &[f64], y: &[f64]) -> f64 {
    x.windows(2)
        .zip(y.windows(2))
        .map(|(xs, ys)| 0.5 * (xs[1] - xs[0]).abs() * (ys[0] + ys[1]))
        .sum()
}
#[cfg(test)]
mod expanded_ml_tests {

    use crate::machine_learning::GaussianProcessRegressor;
    use crate::machine_learning::GradientBoosting;
    use crate::machine_learning::IsolationForest;

    use crate::machine_learning::KnnMetric;
    use crate::machine_learning::KnnModel;

    use crate::machine_learning::RegressionStump;

    use crate::machine_learning::auc_trapezoid;
    use crate::machine_learning::average_path_length;
    use crate::machine_learning::feature_bagging;
    use crate::machine_learning::lof_scores;
    use crate::machine_learning::multiclass_prf;
    use crate::machine_learning::random_forest_oob;
    use crate::machine_learning::rbf_kernel;
    use crate::machine_learning::roc_curve;
    use crate::machine_learning::tsne;
    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }
    #[test]
    fn test_stump_constant_target() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]];
        let r = vec![1.0, 1.0, 1.0, 1.0];
        let stump = RegressionStump::fit(&x, &r);
        assert!(approx(stump.predict_one(&[0.5]), 1.0, 1e-6));
        assert!(approx(stump.predict_one(&[2.5]), 1.0, 1e-6));
    }
    #[test]
    fn test_stump_split() {
        let x = vec![vec![0.0], vec![1.0], vec![5.0], vec![6.0]];
        let r = vec![0.0, 0.0, 1.0, 1.0];
        let stump = RegressionStump::fit(&x, &r);
        assert!(stump.predict_one(&[0.5]) < 0.5);
        assert!(stump.predict_one(&[5.5]) > 0.5);
    }
    #[test]
    fn test_gbdt_mse_decreases() {
        let x: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = (0..20).map(|i| (i as f64) * 0.5).collect();
        let gb1 = GradientBoosting::fit(&x, &y, 5, 0.3);
        let gb2 = GradientBoosting::fit(&x, &y, 50, 0.3);
        assert!(gb2.train_mse(&x, &y) < gb1.train_mse(&x, &y) + 0.5);
    }
    #[test]
    fn test_gbdt_predict_len() {
        let x = vec![vec![1.0], vec![2.0], vec![3.0]];
        let y = vec![1.0, 2.0, 3.0];
        let gb = GradientBoosting::fit(&x, &y, 10, 0.1);
        let preds = gb.predict(&x);
        assert_eq!(preds.len(), 3);
    }
    #[test]
    fn test_feature_bagging_len() {
        let selected = feature_bagging(10, 3, 42);
        assert_eq!(selected.len(), 3);
    }
    #[test]
    fn test_feature_bagging_clamped() {
        let selected = feature_bagging(4, 10, 7);
        assert_eq!(selected.len(), 4);
    }
    #[test]
    fn test_oob_accuracy_separable() {
        let x: Vec<Vec<f64>> = (0..30).map(|i| vec![i as f64]).collect();
        let y: Vec<f64> = (0..30).map(|i| if i < 15 { 0.0 } else { 1.0 }).collect();
        let oob = random_forest_oob(&x, &y, 10, 4);
        assert!(oob.oob_accuracy >= 0.5, "oob acc={}", oob.oob_accuracy);
    }
    #[test]
    fn test_oob_len() {
        let x = vec![vec![0.0], vec![1.0], vec![2.0], vec![3.0]];
        let y = vec![0.0, 0.0, 1.0, 1.0];
        let oob = random_forest_oob(&x, &y, 4, 3);
        assert_eq!(oob.oob_predictions.len(), 4);
    }
    #[test]
    fn test_knn_model_classify_euclidean() {
        let mut knn = KnnModel::new(3, KnnMetric::Euclidean, false);
        knn.fit(
            vec![
                vec![0.0],
                vec![0.5],
                vec![1.0],
                vec![9.0],
                vec![9.5],
                vec![10.0],
            ],
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        );
        assert!(approx(knn.classify(&[0.2]), 0.0, 0.5));
        assert!(approx(knn.classify(&[9.8]), 1.0, 0.5));
    }
    #[test]
    fn test_knn_model_regress() {
        let mut knn = KnnModel::new(2, KnnMetric::Euclidean, false);
        knn.fit(vec![vec![0.0], vec![1.0], vec![2.0]], vec![0.0, 1.0, 2.0]);
        let pred = knn.regress(&[1.0]);
        assert!(approx(pred, 1.0, 0.6));
    }
    #[test]
    fn test_knn_model_weighted_classify() {
        let mut knn = KnnModel::new(3, KnnMetric::Euclidean, true);
        knn.fit(
            vec![vec![0.0], vec![1.0], vec![2.0], vec![10.0]],
            vec![0.0, 0.0, 0.0, 1.0],
        );
        let pred = knn.classify(&[0.1]);
        assert!(approx(pred, 0.0, 0.5));
    }
    #[test]
    fn test_knn_model_manhattan() {
        let mut knn = KnnModel::new(1, KnnMetric::Manhattan, false);
        knn.fit(vec![vec![0.0, 0.0], vec![3.0, 4.0]], vec![0.0, 1.0]);
        assert!(approx(knn.classify(&[0.1, 0.1]), 0.0, 0.5));
    }
    #[test]
    fn test_rbf_kernel_self() {
        let x = vec![1.0, 2.0, 3.0];
        let k = rbf_kernel(&x, &x, 1.0, 1.0);
        assert!(approx(k, 1.0, 1e-10));
    }
    #[test]
    fn test_rbf_kernel_decay() {
        let x1 = vec![0.0];
        let x2 = vec![1.0];
        let x3 = vec![5.0];
        let k12 = rbf_kernel(&x1, &x2, 1.0, 1.0);
        let k13 = rbf_kernel(&x1, &x3, 1.0, 1.0);
        assert!(k12 > k13);
    }
    #[test]
    fn test_gp_regression_interpolation() {
        let x_train: Vec<Vec<f64>> = vec![vec![0.0], vec![1.0], vec![2.0]];
        let y_train = vec![0.0, 1.0, 0.0];
        let gp = GaussianProcessRegressor::fit(x_train.clone(), y_train.clone(), 1.0, 1.0, 1e-3);
        let (means, _vars) = gp.predict(&x_train);
        for (m, t) in means.iter().zip(y_train.iter()) {
            assert!(approx(*m, *t, 0.3), "m={m}, t={t}");
        }
    }
    #[test]
    fn test_gp_variance_non_negative() {
        let x_train = vec![vec![0.0], vec![1.0]];
        let y_train = vec![0.0, 1.0];
        let gp = GaussianProcessRegressor::fit(x_train, y_train, 1.0, 1.0, 0.1);
        let x_test = vec![vec![0.5], vec![2.0]];
        let (_means, vars) = gp.predict(&x_test);
        for &v in &vars {
            assert!(v >= 0.0, "variance={v}");
        }
    }
    #[test]
    fn test_gp_log_marginal_likelihood() {
        let x_train = vec![vec![0.0], vec![1.0], vec![2.0]];
        let y_train = vec![0.0, 1.0, 4.0];
        let gp = GaussianProcessRegressor::fit(x_train, y_train, 1.0, 1.0, 0.1);
        let lml = gp.log_marginal_likelihood();
        assert!(lml.is_finite(), "lml={lml}");
    }
    #[test]
    fn test_isolation_forest_outlier_score() {
        let inliers: Vec<Vec<f64>> = (0..20).map(|i| vec![i as f64 * 0.1]).collect();
        let outlier = vec![vec![100.0_f64]];
        let mut all_data = inliers.clone();
        all_data.extend(outlier.clone());
        let iso = IsolationForest::fit(&all_data, 20, 10);
        let inlier_score = iso.anomaly_score(&inliers[0]);
        let outlier_score = iso.anomaly_score(&outlier[0]);
        assert!(
            outlier_score > inlier_score,
            "outlier={outlier_score}, inlier={inlier_score}"
        );
    }
    #[test]
    fn test_isolation_forest_score_range() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();
        let iso = IsolationForest::fit(&data, 5, 5);
        let scores = iso.score_samples(&data);
        for &s in &scores {
            assert!(s > 0.0 && s <= 1.0, "score={s}");
        }
    }
    #[test]
    fn test_avg_path_length_grows() {
        assert!(average_path_length(10) > average_path_length(2));
        assert!(average_path_length(100) > average_path_length(10));
    }
    #[test]
    fn test_lof_compact_cluster() {
        let data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64 * 0.01]).collect();
        let scores = lof_scores(&data, 3);
        for &s in &scores {
            assert!(s < 3.0, "lof={s}");
        }
    }
    #[test]
    fn test_lof_outlier() {
        let mut data: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64 * 0.1]).collect();
        data.push(vec![100.0]);
        let scores = lof_scores(&data, 3);
        let outlier_score = *scores.last().unwrap();
        let inlier_score = scores[0];
        assert!(
            outlier_score > inlier_score,
            "outlier={outlier_score}, inlier={inlier_score}"
        );
    }
    #[test]
    fn test_tsne_embedding_size() {
        let x: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64, (i % 3) as f64]).collect();
        let result = tsne(&x, 3.0, 20, 100.0);
        assert_eq!(result.embedding.len(), 10);
    }
    #[test]
    fn test_multiclass_prf_perfect() {
        let cm = vec![vec![5usize, 0, 0], vec![0, 4, 0], vec![0, 0, 3]];
        let prf = multiclass_prf(&cm);
        for (_, p, r, f1) in &prf {
            assert!(approx(*p, 1.0, 1e-10));
            assert!(approx(*r, 1.0, 1e-10));
            assert!(approx(*f1, 1.0, 1e-10));
        }
    }
    #[test]
    fn test_roc_auc_trapezoid() {
        let scores = vec![0.9, 0.8, 0.4, 0.2];
        let labels = vec![1usize, 1, 0, 0];
        let (fpr, tpr) = roc_curve(&scores, &labels);
        let auc = auc_trapezoid(&fpr, &tpr);
        assert!(approx(auc, 1.0, 1e-6), "auc={auc}");
    }
}
