// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Machine learning model I/O for the OxiPhysics engine.
//!
//! Covers:
//! - ML model serialization (weights/biases as binary)
//! - ONNX-like simplified format (op-graph with typed tensors)
//! - PyTorch-like state dict (string-keyed tensor store)
//! - Dataset I/O (training/validation split, shuffle)
//! - Feature normalization parameter storage
//! - Label encoding / decoding
//! - Confusion matrix export
//! - Training history (loss/accuracy per epoch)
//! - Hyperparameter configuration
//! - Model checkpoint with metadata

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Tensor — flat f64 storage with shape
// ---------------------------------------------------------------------------

/// A multi-dimensional tensor stored as a flat `Vec`f64`.
#[derive(Debug, Clone, PartialEq)]
pub struct Tensor {
    /// Shape of the tensor (row-major).
    pub shape: Vec<usize>,
    /// Flat data in row-major order.
    pub data: Vec<f64>,
}

impl Tensor {
    /// Construct a tensor with the given shape and flat data.
    ///
    /// # Panics
    /// Panics if the length of `data` does not match the product of `shape`.
    pub fn new(shape: Vec<usize>, data: Vec<f64>) -> Self {
        let expected: usize = shape.iter().product();
        assert_eq!(
            data.len(),
            expected,
            "data length {} does not match shape {:?} (product {})",
            data.len(),
            shape,
            expected
        );
        Tensor { shape, data }
    }

    /// Create a zero tensor with the given shape.
    pub fn zeros(shape: Vec<usize>) -> Self {
        let n: usize = shape.iter().product();
        Tensor {
            shape,
            data: vec![0.0; n],
        }
    }

    /// Total number of elements.
    pub fn numel(&self) -> usize {
        self.data.len()
    }

    /// Number of dimensions.
    pub fn ndim(&self) -> usize {
        self.shape.len()
    }

    /// Serialise to little-endian bytes: `\[ndim u64\]\[dim0 u64\]...\[elem f64\]...`
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::with_capacity(8 + 8 * self.shape.len() + 8 * self.data.len());
        buf.extend_from_slice(&(self.shape.len() as u64).to_le_bytes());
        for &d in &self.shape {
            buf.extend_from_slice(&(d as u64).to_le_bytes());
        }
        for &v in &self.data {
            buf.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        buf
    }

    /// Deserialise from bytes produced by [`Tensor::to_bytes`].
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() < 8 {
            return None;
        }
        let ndim = u64::from_le_bytes(bytes[0..8].try_into().ok()?) as usize;
        let header_len = 8 + 8 * ndim;
        if bytes.len() < header_len {
            return None;
        }
        let mut shape = Vec::with_capacity(ndim);
        for i in 0..ndim {
            let off = 8 + 8 * i;
            shape.push(u64::from_le_bytes(bytes[off..off + 8].try_into().ok()?) as usize);
        }
        let n: usize = shape.iter().product();
        if bytes.len() < header_len + 8 * n {
            return None;
        }
        let mut data = Vec::with_capacity(n);
        for i in 0..n {
            let off = header_len + 8 * i;
            let bits = u64::from_le_bytes(bytes[off..off + 8].try_into().ok()?);
            data.push(f64::from_bits(bits));
        }
        Some(Tensor { shape, data })
    }

    /// Element-wise add (must have identical shape).
    pub fn add(&self, other: &Tensor) -> Option<Tensor> {
        if self.shape != other.shape {
            return None;
        }
        let data = self
            .data
            .iter()
            .zip(&other.data)
            .map(|(a, b)| a + b)
            .collect();
        Some(Tensor {
            shape: self.shape.clone(),
            data,
        })
    }

    /// Scalar multiply.
    pub fn scale(&self, s: f64) -> Tensor {
        Tensor {
            shape: self.shape.clone(),
            data: self.data.iter().map(|v| v * s).collect(),
        }
    }

    /// Compute sum of all elements.
    pub fn sum(&self) -> f64 {
        self.data.iter().sum()
    }

    /// Compute mean of all elements.
    pub fn mean(&self) -> f64 {
        if self.data.is_empty() {
            return 0.0;
        }
        self.sum() / self.data.len() as f64
    }
}

// ---------------------------------------------------------------------------
// Layer — single dense layer
// ---------------------------------------------------------------------------

/// A single dense (fully connected) layer with weights and biases.
#[derive(Debug, Clone)]
pub struct DenseLayer {
    /// Layer name.
    pub name: String,
    /// Weight matrix: shape `\[out_features, in_features\]`.
    pub weights: Tensor,
    /// Bias vector: shape `[out_features]`.
    pub bias: Tensor,
    /// Activation function name (e.g. `"relu"`, `"sigmoid"`, `"tanh"`, `"linear"`).
    pub activation: String,
}

impl DenseLayer {
    /// Construct a zero-initialised dense layer.
    pub fn new(
        name: impl Into<String>,
        in_features: usize,
        out_features: usize,
        activation: impl Into<String>,
    ) -> Self {
        DenseLayer {
            name: name.into(),
            weights: Tensor::zeros(vec![out_features, in_features]),
            bias: Tensor::zeros(vec![out_features]),
            activation: activation.into(),
        }
    }

    /// Forward pass: `output\[i\] = sum_j(w\[i,j\] * input\[j\]) + bias\[i\]`, then activation.
    pub fn forward(&self, input: &[f64]) -> Vec<f64> {
        let in_feat = input.len();
        let out_feat = self.bias.data.len();
        let mut out = vec![0.0f64; out_feat];
        for (i, (o, &b)) in out.iter_mut().zip(self.bias.data.iter()).enumerate() {
            let mut acc = b;
            let cols = in_feat.min(self.weights.data.len() / out_feat);
            for (w, &x) in self.weights.data[i * in_feat..i * in_feat + cols]
                .iter()
                .zip(input[..cols].iter())
            {
                acc += w * x;
            }
            *o = apply_activation(acc, &self.activation);
        }
        out
    }

    /// Number of trainable parameters.
    pub fn param_count(&self) -> usize {
        self.weights.numel() + self.bias.numel()
    }

    /// Serialise to bytes (name length, name utf8, weights bytes, bias bytes).
    pub fn to_bytes(&self) -> Vec<u8> {
        let name_bytes = self.name.as_bytes();
        let act_bytes = self.activation.as_bytes();
        let mut buf = Vec::new();
        buf.extend_from_slice(&(name_bytes.len() as u64).to_le_bytes());
        buf.extend_from_slice(name_bytes);
        buf.extend_from_slice(&(act_bytes.len() as u64).to_le_bytes());
        buf.extend_from_slice(act_bytes);
        let wb = self.weights.to_bytes();
        buf.extend_from_slice(&(wb.len() as u64).to_le_bytes());
        buf.extend_from_slice(&wb);
        let bb = self.bias.to_bytes();
        buf.extend_from_slice(&(bb.len() as u64).to_le_bytes());
        buf.extend_from_slice(&bb);
        buf
    }
}

/// Apply a named activation function to a scalar.
pub fn apply_activation(x: f64, activation: &str) -> f64 {
    match activation {
        "relu" => x.max(0.0),
        "sigmoid" => 1.0 / (1.0 + (-x).exp()),
        "tanh" => x.tanh(),
        "softplus" => (1.0 + x.exp()).ln(),
        "elu" => {
            if x >= 0.0 {
                x
            } else {
                x.exp() - 1.0
            }
        }
        "leaky_relu" => {
            if x >= 0.0 {
                x
            } else {
                0.01 * x
            }
        }
        _ => x, // linear / identity
    }
}

// ---------------------------------------------------------------------------
// ModelWeights — named layers with binary I/O
// ---------------------------------------------------------------------------

/// A collection of named dense layers (binary-serialisable model weights).
#[derive(Debug, Clone, Default)]
pub struct ModelWeights {
    /// Layers in order.
    pub layers: Vec<DenseLayer>,
}

impl ModelWeights {
    /// Create an empty model.
    pub fn new() -> Self {
        ModelWeights { layers: Vec::new() }
    }

    /// Append a layer.
    pub fn add_layer(&mut self, layer: DenseLayer) {
        self.layers.push(layer);
    }

    /// Look up a layer by name.
    pub fn get_layer(&self, name: &str) -> Option<&DenseLayer> {
        self.layers.iter().find(|l| l.name == name)
    }

    /// Total trainable parameter count.
    pub fn total_params(&self) -> usize {
        self.layers.iter().map(|l| l.param_count()).sum()
    }

    /// Serialise all layers to a flat byte buffer.
    ///
    /// Format: `\[layer_count: u64\]\[layer_0_len: u64\][layer_0_bytes]...`
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&(self.layers.len() as u64).to_le_bytes());
        for layer in &self.layers {
            let lb = layer.to_bytes();
            buf.extend_from_slice(&(lb.len() as u64).to_le_bytes());
            buf.extend_from_slice(&lb);
        }
        buf
    }
}

// ---------------------------------------------------------------------------
// StateDict — PyTorch-like key-value tensor store
// ---------------------------------------------------------------------------

/// PyTorch-like state dict: a `HashMap<String, Tensor>`.
#[derive(Debug, Clone, Default)]
pub struct StateDict {
    /// The underlying key-value store.
    pub tensors: HashMap<String, Tensor>,
}

impl StateDict {
    /// Create an empty state dict.
    pub fn new() -> Self {
        StateDict {
            tensors: HashMap::new(),
        }
    }

    /// Insert a tensor under a key.
    pub fn insert(&mut self, key: impl Into<String>, tensor: Tensor) {
        self.tensors.insert(key.into(), tensor);
    }

    /// Retrieve a tensor by key.
    pub fn get(&self, key: &str) -> Option<&Tensor> {
        self.tensors.get(key)
    }

    /// Number of tensors.
    pub fn len(&self) -> usize {
        self.tensors.len()
    }

    /// Whether the dict is empty.
    pub fn is_empty(&self) -> bool {
        self.tensors.is_empty()
    }

    /// Total number of parameters (sum of all tensor element counts).
    pub fn total_params(&self) -> usize {
        self.tensors.values().map(|t| t.numel()).sum()
    }

    /// Serialise to bytes.
    ///
    /// Format: `\[entry_count: u64\](\[key_len: u64\][key_bytes]\[tensor_len: u64\][tensor_bytes])...`
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&(self.tensors.len() as u64).to_le_bytes());
        let mut keys: Vec<&String> = self.tensors.keys().collect();
        keys.sort(); // deterministic order
        for k in keys {
            let kb = k.as_bytes();
            buf.extend_from_slice(&(kb.len() as u64).to_le_bytes());
            buf.extend_from_slice(kb);
            let tb = self.tensors[k].to_bytes();
            buf.extend_from_slice(&(tb.len() as u64).to_le_bytes());
            buf.extend_from_slice(&tb);
        }
        buf
    }

    /// Deserialise from bytes produced by [`StateDict::to_bytes`].
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        let mut pos = 0usize;
        let n = read_u64(bytes, &mut pos)? as usize;
        let mut dict = StateDict::new();
        for _ in 0..n {
            let klen = read_u64(bytes, &mut pos)? as usize;
            if pos + klen > bytes.len() {
                return None;
            }
            let key = String::from_utf8(bytes[pos..pos + klen].to_vec()).ok()?;
            pos += klen;
            let tlen = read_u64(bytes, &mut pos)? as usize;
            if pos + tlen > bytes.len() {
                return None;
            }
            let tensor = Tensor::from_bytes(&bytes[pos..pos + tlen])?;
            pos += tlen;
            dict.insert(key, tensor);
        }
        Some(dict)
    }
}

/// Read a little-endian `u64` from `bytes` at `*pos`, advancing `*pos` by 8.
fn read_u64(bytes: &[u8], pos: &mut usize) -> Option<u64> {
    if *pos + 8 > bytes.len() {
        return None;
    }
    let v = u64::from_le_bytes(bytes[*pos..*pos + 8].try_into().ok()?);
    *pos += 8;
    Some(v)
}

// ---------------------------------------------------------------------------
// OnnxLikeGraph — simplified ONNX-style operator graph
// ---------------------------------------------------------------------------

/// A single operation node in an ONNX-like compute graph.
#[derive(Debug, Clone)]
pub struct OnnxNode {
    /// Unique node name.
    pub name: String,
    /// Operation type (e.g. `"MatMul"`, `"Relu"`, `"Add"`, `"Sigmoid"`).
    pub op_type: String,
    /// Names of input tensors.
    pub inputs: Vec<String>,
    /// Names of output tensors.
    pub outputs: Vec<String>,
    /// Optional scalar attributes (key → value).
    pub attributes: HashMap<String, f64>,
}

impl OnnxNode {
    /// Construct a new node with no attributes.
    pub fn new(
        name: impl Into<String>,
        op_type: impl Into<String>,
        inputs: Vec<String>,
        outputs: Vec<String>,
    ) -> Self {
        OnnxNode {
            name: name.into(),
            op_type: op_type.into(),
            inputs,
            outputs,
            attributes: HashMap::new(),
        }
    }

    /// Set an attribute.
    pub fn with_attr(mut self, key: impl Into<String>, value: f64) -> Self {
        self.attributes.insert(key.into(), value);
        self
    }
}

/// A simplified ONNX-like computation graph.
#[derive(Debug, Clone, Default)]
pub struct OnnxLikeGraph {
    /// Ordered list of operation nodes.
    pub nodes: Vec<OnnxNode>,
    /// Named initialiser tensors (model weights).
    pub initializers: StateDict,
    /// Input tensor names.
    pub inputs: Vec<String>,
    /// Output tensor names.
    pub outputs: Vec<String>,
    /// Graph name.
    pub name: String,
}

impl OnnxLikeGraph {
    /// Create an empty graph.
    pub fn new(name: impl Into<String>) -> Self {
        OnnxLikeGraph {
            name: name.into(),
            nodes: Vec::new(),
            initializers: StateDict::new(),
            inputs: Vec::new(),
            outputs: Vec::new(),
        }
    }

    /// Add an operation node.
    pub fn add_node(&mut self, node: OnnxNode) {
        self.nodes.push(node);
    }

    /// Add an initialiser.
    pub fn add_initializer(&mut self, name: impl Into<String>, tensor: Tensor) {
        self.initializers.insert(name, tensor);
    }

    /// Number of nodes in the graph.
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Count nodes by op type.
    pub fn count_op(&self, op: &str) -> usize {
        self.nodes.iter().filter(|n| n.op_type == op).count()
    }

    /// Topological order check: returns `true` if all node input names
    /// are either graph inputs or outputs of an earlier node.
    pub fn is_topologically_valid(&self) -> bool {
        let mut available: std::collections::HashSet<&str> =
            self.inputs.iter().map(|s| s.as_str()).collect();
        // Initialisers are always available
        for k in self.initializers.tensors.keys() {
            available.insert(k.as_str());
        }
        for node in &self.nodes {
            for inp in &node.inputs {
                if !available.contains(inp.as_str()) {
                    return false;
                }
            }
            for out in &node.outputs {
                available.insert(out.as_str());
            }
        }
        true
    }
}

// ---------------------------------------------------------------------------
// Dataset — rows of features with optional labels
// ---------------------------------------------------------------------------

/// A dataset row: a feature vector and an optional label index.
#[derive(Debug, Clone)]
pub struct DataRow {
    /// Feature values.
    pub features: Vec<f64>,
    /// Optional class label index.
    pub label: Option<usize>,
}

impl DataRow {
    /// Construct a labelled row.
    pub fn labelled(features: Vec<f64>, label: usize) -> Self {
        DataRow {
            features,
            label: Some(label),
        }
    }

    /// Construct an unlabelled row.
    pub fn unlabelled(features: Vec<f64>) -> Self {
        DataRow {
            features,
            label: None,
        }
    }
}

/// A dataset with optional train/validation split.
#[derive(Debug, Clone, Default)]
pub struct Dataset {
    /// All rows.
    pub rows: Vec<DataRow>,
    /// Feature names.
    pub feature_names: Vec<String>,
    /// Class names (index → name).
    pub class_names: Vec<String>,
}

impl Dataset {
    /// Create an empty dataset.
    pub fn new() -> Self {
        Dataset {
            rows: Vec::new(),
            feature_names: Vec::new(),
            class_names: Vec::new(),
        }
    }

    /// Add a row.
    pub fn push(&mut self, row: DataRow) {
        self.rows.push(row);
    }

    /// Number of rows.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    /// Whether the dataset is empty.
    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Number of features (from the first row, or 0).
    pub fn num_features(&self) -> usize {
        self.rows.first().map(|r| r.features.len()).unwrap_or(0)
    }

    /// Shuffle rows using Fisher-Yates with a simple LCG.
    pub fn shuffle(&mut self, seed: u64) {
        let n = self.rows.len();
        if n < 2 {
            return;
        }
        let mut rng = LcgRng::new(seed);
        for i in (1..n).rev() {
            let j = rng.next_usize_below(i + 1);
            self.rows.swap(i, j);
        }
    }

    /// Split into training and validation subsets.
    ///
    /// `val_fraction` is the fraction of rows reserved for validation.
    pub fn train_val_split(&self, val_fraction: f64) -> (Dataset, Dataset) {
        let val_count = ((self.rows.len() as f64) * val_fraction.clamp(0.0, 1.0)) as usize;
        let train_count = self.rows.len().saturating_sub(val_count);
        let mut train = Dataset {
            rows: self.rows[..train_count].to_vec(),
            feature_names: self.feature_names.clone(),
            class_names: self.class_names.clone(),
        };
        let mut val = Dataset {
            rows: self.rows[train_count..].to_vec(),
            feature_names: self.feature_names.clone(),
            class_names: self.class_names.clone(),
        };
        // suppress unused warnings
        let _ = &mut train;
        let _ = &mut val;
        (train, val)
    }

    /// Compute per-feature mean and standard deviation.
    ///
    /// Returns `(means, stds)` each of length `num_features`.
    pub fn feature_stats(&self) -> (Vec<f64>, Vec<f64>) {
        let nf = self.num_features();
        if nf == 0 || self.rows.is_empty() {
            return (vec![], vec![]);
        }
        let n = self.rows.len() as f64;
        let mut means = vec![0.0f64; nf];
        for row in &self.rows {
            for (k, &v) in row.features.iter().enumerate() {
                means[k] += v;
            }
        }
        for m in &mut means {
            *m /= n;
        }
        let mut stds = vec![0.0f64; nf];
        for row in &self.rows {
            for (k, &v) in row.features.iter().enumerate() {
                let d = v - means[k];
                stds[k] += d * d;
            }
        }
        for s in &mut stds {
            *s = (*s / n).sqrt();
        }
        (means, stds)
    }
}

// ---------------------------------------------------------------------------
// Minimal LCG RNG for shuffle (no rand dep in this module)
// ---------------------------------------------------------------------------

/// A minimal linear congruential generator used for dataset shuffling.
struct LcgRng {
    state: u64,
}

impl LcgRng {
    fn new(seed: u64) -> Self {
        LcgRng {
            state: seed ^ 0x1234_5678_9abc_def0,
        }
    }

    fn next_u64(&mut self) -> u64 {
        // Knuth's MMIX constants
        self.state = self
            .state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.state
    }

    fn next_usize_below(&mut self, n: usize) -> usize {
        if n == 0 {
            return 0;
        }
        (self.next_u64() % n as u64) as usize
    }
}

// ---------------------------------------------------------------------------
// NormalizationParams — feature normalization storage
// ---------------------------------------------------------------------------

/// Stored feature normalization parameters (mean and std for z-score normalization).
#[derive(Debug, Clone)]
pub struct NormalizationParams {
    /// Per-feature mean.
    pub means: Vec<f64>,
    /// Per-feature standard deviation.
    pub stds: Vec<f64>,
    /// Minimum values (for min-max normalization).
    pub mins: Vec<f64>,
    /// Maximum values (for min-max normalization).
    pub maxs: Vec<f64>,
}

impl NormalizationParams {
    /// Compute z-score normalization parameters from a dataset.
    pub fn from_dataset(dataset: &Dataset) -> Self {
        let (means, stds) = dataset.feature_stats();
        let nf = means.len();
        let mut mins = vec![f64::INFINITY; nf];
        let mut maxs = vec![f64::NEG_INFINITY; nf];
        for row in &dataset.rows {
            for (k, &v) in row.features.iter().enumerate() {
                if v < mins[k] {
                    mins[k] = v;
                }
                if v > maxs[k] {
                    maxs[k] = v;
                }
            }
        }
        NormalizationParams {
            means,
            stds,
            mins,
            maxs,
        }
    }

    /// Apply z-score normalization to a feature vector.
    pub fn normalize_zscore(&self, features: &[f64]) -> Vec<f64> {
        features
            .iter()
            .enumerate()
            .map(|(k, &v)| {
                let s = if k < self.stds.len() {
                    self.stds[k]
                } else {
                    1.0
                };
                let m = if k < self.means.len() {
                    self.means[k]
                } else {
                    0.0
                };
                if s.abs() < 1e-15 { 0.0 } else { (v - m) / s }
            })
            .collect()
    }

    /// Apply min-max normalization to a feature vector (maps to [0, 1]).
    pub fn normalize_minmax(&self, features: &[f64]) -> Vec<f64> {
        features
            .iter()
            .enumerate()
            .map(|(k, &v)| {
                let mn = if k < self.mins.len() {
                    self.mins[k]
                } else {
                    0.0
                };
                let mx = if k < self.maxs.len() {
                    self.maxs[k]
                } else {
                    1.0
                };
                let range = mx - mn;
                if range.abs() < 1e-15 {
                    0.0
                } else {
                    (v - mn) / range
                }
            })
            .collect()
    }

    /// Serialise to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        let write_vec = |buf: &mut Vec<u8>, v: &[f64]| {
            buf.extend_from_slice(&(v.len() as u64).to_le_bytes());
            for &x in v {
                buf.extend_from_slice(&x.to_bits().to_le_bytes());
            }
        };
        write_vec(&mut buf, &self.means);
        write_vec(&mut buf, &self.stds);
        write_vec(&mut buf, &self.mins);
        write_vec(&mut buf, &self.maxs);
        buf
    }
}

// ---------------------------------------------------------------------------
// LabelEncoder — integer ↔ class-name mapping
// ---------------------------------------------------------------------------

/// Encodes class labels as integers and decodes them back.
#[derive(Debug, Clone, Default)]
pub struct LabelEncoder {
    /// Class names in order (index 0 = first class).
    pub classes: Vec<String>,
    /// Reverse map: class name → index.
    index: HashMap<String, usize>,
}

impl LabelEncoder {
    /// Create an empty encoder.
    pub fn new() -> Self {
        LabelEncoder {
            classes: Vec::new(),
            index: HashMap::new(),
        }
    }

    /// Fit the encoder from a list of class names.
    pub fn fit(mut class_names: Vec<String>) -> Self {
        class_names.sort();
        class_names.dedup();
        let index = class_names
            .iter()
            .enumerate()
            .map(|(i, s)| (s.clone(), i))
            .collect();
        LabelEncoder {
            classes: class_names,
            index,
        }
    }

    /// Encode a class name to its index.
    pub fn encode(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    /// Decode an index to its class name.
    pub fn decode(&self, idx: usize) -> Option<&str> {
        self.classes.get(idx).map(|s| s.as_str())
    }

    /// Number of classes.
    pub fn num_classes(&self) -> usize {
        self.classes.len()
    }

    /// One-hot encode an index.
    pub fn one_hot(&self, idx: usize) -> Vec<f64> {
        let mut v = vec![0.0f64; self.num_classes()];
        if idx < v.len() {
            v[idx] = 1.0;
        }
        v
    }
}

// ---------------------------------------------------------------------------
// ConfusionMatrix
// ---------------------------------------------------------------------------

/// Confusion matrix for multi-class classification.
#[derive(Debug, Clone)]
pub struct ConfusionMatrix {
    /// Number of classes.
    pub num_classes: usize,
    /// Matrix data: `counts\[true_label * num_classes + predicted_label\]`.
    pub counts: Vec<u64>,
}

impl ConfusionMatrix {
    /// Create a zero-filled confusion matrix for `num_classes` classes.
    pub fn new(num_classes: usize) -> Self {
        ConfusionMatrix {
            num_classes,
            counts: vec![0; num_classes * num_classes],
        }
    }

    /// Record a single prediction.
    pub fn record(&mut self, true_label: usize, predicted: usize) {
        if true_label < self.num_classes && predicted < self.num_classes {
            self.counts[true_label * self.num_classes + predicted] += 1;
        }
    }

    /// Overall accuracy: fraction of correct predictions.
    pub fn accuracy(&self) -> f64 {
        let total: u64 = self.counts.iter().sum();
        if total == 0 {
            return 0.0;
        }
        let correct: u64 = (0..self.num_classes)
            .map(|i| self.counts[i * self.num_classes + i])
            .sum();
        correct as f64 / total as f64
    }

    /// Per-class precision: `TP / (TP + FP)`.
    pub fn precision(&self, class: usize) -> f64 {
        if class >= self.num_classes {
            return 0.0;
        }
        let tp = self.counts[class * self.num_classes + class] as f64;
        let fp: f64 = (0..self.num_classes)
            .filter(|&r| r != class)
            .map(|r| self.counts[r * self.num_classes + class] as f64)
            .sum();
        if tp + fp < 1e-15 { 0.0 } else { tp / (tp + fp) }
    }

    /// Per-class recall: `TP / (TP + FN)`.
    pub fn recall(&self, class: usize) -> f64 {
        if class >= self.num_classes {
            return 0.0;
        }
        let tp = self.counts[class * self.num_classes + class] as f64;
        let fn_: f64 = (0..self.num_classes)
            .filter(|&c| c != class)
            .map(|c| self.counts[class * self.num_classes + c] as f64)
            .sum();
        if tp + fn_ < 1e-15 {
            0.0
        } else {
            tp / (tp + fn_)
        }
    }

    /// Per-class F1 score.
    pub fn f1(&self, class: usize) -> f64 {
        let p = self.precision(class);
        let r = self.recall(class);
        if p + r < 1e-15 {
            0.0
        } else {
            2.0 * p * r / (p + r)
        }
    }

    /// Export as a CSV string.
    pub fn to_csv(&self) -> String {
        let mut s = String::new();
        // Header
        s.push_str("true\\pred");
        for j in 0..self.num_classes {
            s.push_str(&format!(",class_{j}"));
        }
        s.push('\n');
        for i in 0..self.num_classes {
            s.push_str(&format!("class_{i}"));
            for j in 0..self.num_classes {
                s.push_str(&format!(",{}", self.counts[i * self.num_classes + j]));
            }
            s.push('\n');
        }
        s
    }
}

// ---------------------------------------------------------------------------
// TrainingHistory
// ---------------------------------------------------------------------------

/// Per-epoch metrics.
#[derive(Debug, Clone)]
pub struct EpochRecord {
    /// Epoch number (0-indexed).
    pub epoch: usize,
    /// Training loss.
    pub train_loss: f64,
    /// Validation loss.
    pub val_loss: f64,
    /// Training accuracy.
    pub train_acc: f64,
    /// Validation accuracy.
    pub val_acc: f64,
    /// Learning rate used.
    pub learning_rate: f64,
}

/// Full training history for a model.
#[derive(Debug, Clone, Default)]
pub struct TrainingHistory {
    /// Records, one per epoch.
    pub records: Vec<EpochRecord>,
}

impl TrainingHistory {
    /// Create an empty history.
    pub fn new() -> Self {
        TrainingHistory {
            records: Vec::new(),
        }
    }

    /// Append an epoch record.
    pub fn push(&mut self, record: EpochRecord) {
        self.records.push(record);
    }

    /// Number of epochs recorded.
    pub fn num_epochs(&self) -> usize {
        self.records.len()
    }

    /// Best validation accuracy and the epoch at which it occurred.
    pub fn best_val_acc(&self) -> Option<(usize, f64)> {
        self.records
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| {
                a.val_acc
                    .partial_cmp(&b.val_acc)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, r)| (i, r.val_acc))
    }

    /// Best (lowest) validation loss and its epoch.
    pub fn best_val_loss(&self) -> Option<(usize, f64)> {
        self.records
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                a.val_loss
                    .partial_cmp(&b.val_loss)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, r)| (i, r.val_loss))
    }

    /// Export history as a CSV string.
    pub fn to_csv(&self) -> String {
        let mut s = String::from("epoch,train_loss,val_loss,train_acc,val_acc,lr\n");
        for r in &self.records {
            s.push_str(&format!(
                "{},{:.6},{:.6},{:.6},{:.6},{:.8}\n",
                r.epoch, r.train_loss, r.val_loss, r.train_acc, r.val_acc, r.learning_rate
            ));
        }
        s
    }
}

// ---------------------------------------------------------------------------
// HyperparamConfig — JSON-compatible hyperparameter store
// ---------------------------------------------------------------------------

/// Typed hyperparameter value.
#[derive(Debug, Clone, PartialEq)]
pub enum HpValue {
    /// Floating point (also covers int values stored as f64).
    Float(f64),
    /// Boolean flag.
    Bool(bool),
    /// String-valued parameter.
    Str(String),
}

impl HpValue {
    /// Return the f64 value if this is a `Float`, else `None`.
    pub fn as_float(&self) -> Option<f64> {
        if let HpValue::Float(v) = self {
            Some(*v)
        } else {
            None
        }
    }

    /// Return the bool value if this is a `Bool`, else `None`.
    pub fn as_bool(&self) -> Option<bool> {
        if let HpValue::Bool(v) = self {
            Some(*v)
        } else {
            None
        }
    }

    /// Return the string ref if this is a `Str`, else `None`.
    pub fn as_str(&self) -> Option<&str> {
        if let HpValue::Str(s) = self {
            Some(s.as_str())
        } else {
            None
        }
    }
}

/// Hyperparameter configuration container.
#[derive(Debug, Clone, Default)]
pub struct HyperparamConfig {
    /// Key-value map.
    pub params: HashMap<String, HpValue>,
}

impl HyperparamConfig {
    /// Create an empty config.
    pub fn new() -> Self {
        HyperparamConfig {
            params: HashMap::new(),
        }
    }

    /// Set a float hyperparameter.
    pub fn set_float(&mut self, key: impl Into<String>, value: f64) {
        self.params.insert(key.into(), HpValue::Float(value));
    }

    /// Set a boolean hyperparameter.
    pub fn set_bool(&mut self, key: impl Into<String>, value: bool) {
        self.params.insert(key.into(), HpValue::Bool(value));
    }

    /// Set a string hyperparameter.
    pub fn set_str(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.params.insert(key.into(), HpValue::Str(value.into()));
    }

    /// Get a float hyperparameter.
    pub fn get_float(&self, key: &str) -> Option<f64> {
        self.params.get(key)?.as_float()
    }

    /// Get a bool hyperparameter.
    pub fn get_bool(&self, key: &str) -> Option<bool> {
        self.params.get(key)?.as_bool()
    }

    /// Get a string hyperparameter.
    pub fn get_str(&self, key: &str) -> Option<&str> {
        self.params.get(key)?.as_str()
    }

    /// Serialise as a simple JSON string (no external dependencies).
    pub fn to_json(&self) -> String {
        let mut parts: Vec<String> = Vec::new();
        let mut keys: Vec<&String> = self.params.keys().collect();
        keys.sort();
        for k in keys {
            let v_str = match &self.params[k] {
                HpValue::Float(f) => format!("{f}"),
                HpValue::Bool(b) => format!("{b}"),
                HpValue::Str(s) => format!("\"{}\"", s.replace('"', "\\\"")),
            };
            parts.push(format!("\"{}\":{}", k.replace('"', "\\\""), v_str));
        }
        format!("{{{}}}", parts.join(","))
    }
}

// ---------------------------------------------------------------------------
// ModelCheckpoint
// ---------------------------------------------------------------------------

/// Metadata stored alongside a model checkpoint.
#[derive(Debug, Clone)]
pub struct CheckpointMeta {
    /// Epoch at which the checkpoint was saved.
    pub epoch: usize,
    /// Validation loss at checkpoint time.
    pub val_loss: f64,
    /// Validation accuracy at checkpoint time.
    pub val_acc: f64,
    /// Wall-clock training time in seconds (cumulative).
    pub train_time_secs: f64,
    /// Model architecture description.
    pub architecture: String,
    /// Framework version string.
    pub framework_version: String,
}

impl CheckpointMeta {
    /// Serialise as a simple text block.
    pub fn to_text(&self) -> String {
        format!(
            "epoch={}\nval_loss={:.8}\nval_acc={:.8}\ntrain_time_secs={:.3}\narchitecture={}\nframework_version={}\n",
            self.epoch,
            self.val_loss,
            self.val_acc,
            self.train_time_secs,
            self.architecture,
            self.framework_version
        )
    }
}

/// A model checkpoint: state dict + metadata + hyperparameters.
#[derive(Debug, Clone)]
pub struct ModelCheckpoint {
    /// Model weights.
    pub state: StateDict,
    /// Checkpoint metadata.
    pub meta: CheckpointMeta,
    /// Hyperparameters used when this checkpoint was saved.
    pub hparams: HyperparamConfig,
}

impl ModelCheckpoint {
    /// Create a new checkpoint.
    pub fn new(state: StateDict, meta: CheckpointMeta, hparams: HyperparamConfig) -> Self {
        ModelCheckpoint {
            state,
            meta,
            hparams,
        }
    }

    /// Serialise to a flat byte buffer.
    ///
    /// Format: `\[state_len: u64\][state_bytes]\[meta_text_len: u64\][meta_text_utf8]\[hp_json_len: u64\][hp_json_utf8]`
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        let sb = self.state.to_bytes();
        buf.extend_from_slice(&(sb.len() as u64).to_le_bytes());
        buf.extend_from_slice(&sb);
        let mt = self.meta.to_text();
        let mb = mt.as_bytes();
        buf.extend_from_slice(&(mb.len() as u64).to_le_bytes());
        buf.extend_from_slice(mb);
        let hp = self.hparams.to_json();
        let hb = hp.as_bytes();
        buf.extend_from_slice(&(hb.len() as u64).to_le_bytes());
        buf.extend_from_slice(hb);
        buf
    }

    /// Byte size of the serialised checkpoint.
    pub fn byte_size(&self) -> usize {
        self.to_bytes().len()
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Compute softmax of a slice.
pub fn softmax(logits: &[f64]) -> Vec<f64> {
    if logits.is_empty() {
        return vec![];
    }
    let max_v = logits.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = logits.iter().map(|&x| (x - max_v).exp()).collect();
    let sum: f64 = exps.iter().sum();
    if sum < 1e-15 {
        vec![1.0 / logits.len() as f64; logits.len()]
    } else {
        exps.iter().map(|e| e / sum).collect()
    }
}

/// Compute cross-entropy loss between `probs` and one-hot `targets`.
pub fn cross_entropy_loss(probs: &[f64], targets: &[f64]) -> f64 {
    probs
        .iter()
        .zip(targets)
        .map(|(&p, &t)| -t * (p.max(1e-15)).ln())
        .sum()
}

/// Argmax: index of the maximum value.
pub fn argmax(values: &[f64]) -> usize {
    values
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Compute mean squared error.
pub fn mse(predictions: &[f64], targets: &[f64]) -> f64 {
    if predictions.is_empty() {
        return 0.0;
    }
    let n = predictions.len().min(targets.len()) as f64;
    predictions
        .iter()
        .zip(targets)
        .map(|(&p, &t)| {
            let d = p - t;
            d * d
        })
        .sum::<f64>()
        / n
}

/// Compute mean absolute error.
pub fn mae(predictions: &[f64], targets: &[f64]) -> f64 {
    if predictions.is_empty() {
        return 0.0;
    }
    let n = predictions.len().min(targets.len()) as f64;
    predictions
        .iter()
        .zip(targets)
        .map(|(&p, &t)| (p - t).abs())
        .sum::<f64>()
        / n
}

// ---------------------------------------------------------------------------
// Unit tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Tensor ---

    #[test]
    fn test_tensor_new_shape_mismatch_panics() {
        let result = std::panic::catch_unwind(|| Tensor::new(vec![2, 3], vec![0.0; 5]));
        assert!(result.is_err());
    }

    #[test]
    fn test_tensor_zeros() {
        let t = Tensor::zeros(vec![3, 4]);
        assert_eq!(t.numel(), 12);
        assert!(t.data.iter().all(|&v| v == 0.0));
    }

    #[test]
    fn test_tensor_numel() {
        let t = Tensor::new(vec![2, 3], vec![1.0; 6]);
        assert_eq!(t.numel(), 6);
        assert_eq!(t.ndim(), 2);
    }

    #[test]
    fn test_tensor_sum_mean() {
        let t = Tensor::new(vec![4], vec![1.0, 2.0, 3.0, 4.0]);
        assert!((t.sum() - 10.0).abs() < 1e-12);
        assert!((t.mean() - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_tensor_scale() {
        let t = Tensor::new(vec![3], vec![1.0, 2.0, 3.0]);
        let t2 = t.scale(2.0);
        assert!((t2.data[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn test_tensor_add() {
        let a = Tensor::new(vec![3], vec![1.0, 2.0, 3.0]);
        let b = Tensor::new(vec![3], vec![4.0, 5.0, 6.0]);
        let c = a.add(&b).unwrap();
        assert!((c.data[2] - 9.0).abs() < 1e-12);
    }

    #[test]
    fn test_tensor_add_shape_mismatch() {
        let a = Tensor::new(vec![2], vec![1.0, 2.0]);
        let b = Tensor::new(vec![3], vec![1.0, 2.0, 3.0]);
        assert!(a.add(&b).is_none());
    }

    #[test]
    fn test_tensor_roundtrip_bytes() {
        let t = Tensor::new(vec![2, 3], vec![1.0, -2.5, 0.0, 3.125, 1e10, -1e-5]);
        let bytes = t.to_bytes();
        let t2 = Tensor::from_bytes(&bytes).unwrap();
        assert_eq!(t2.shape, t.shape);
        for (a, b) in t.data.iter().zip(&t2.data) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    #[test]
    fn test_tensor_from_bytes_empty_is_none() {
        assert!(Tensor::from_bytes(&[]).is_none());
    }

    // --- DenseLayer ---

    #[test]
    fn test_dense_layer_param_count() {
        let layer = DenseLayer::new("fc1", 4, 3, "relu");
        // weights: 3*4=12, bias: 3
        assert_eq!(layer.param_count(), 15);
    }

    #[test]
    fn test_dense_layer_forward_zero_weights() {
        let layer = DenseLayer::new("fc", 3, 2, "linear");
        let input = vec![1.0, 2.0, 3.0];
        let out = layer.forward(&input);
        assert_eq!(out.len(), 2);
        // All zero weights+bias → output is 0
        for v in &out {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn test_dense_layer_activation_relu() {
        assert!((apply_activation(-5.0, "relu")).abs() < 1e-12);
        assert!((apply_activation(3.0, "relu") - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_dense_layer_activation_sigmoid() {
        let v = apply_activation(0.0, "sigmoid");
        assert!((v - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_dense_layer_activation_tanh() {
        let v = apply_activation(0.0, "tanh");
        assert!(v.abs() < 1e-12);
    }

    // --- ModelWeights ---

    #[test]
    fn test_model_weights_add_and_get() {
        let mut model = ModelWeights::new();
        model.add_layer(DenseLayer::new("l1", 4, 8, "relu"));
        model.add_layer(DenseLayer::new("l2", 8, 2, "sigmoid"));
        assert_eq!(model.layers.len(), 2);
        assert!(model.get_layer("l1").is_some());
        assert!(model.get_layer("l3").is_none());
    }

    #[test]
    fn test_model_weights_total_params() {
        let mut model = ModelWeights::new();
        model.add_layer(DenseLayer::new("l1", 4, 3, "relu")); // 12+3=15
        model.add_layer(DenseLayer::new("l2", 3, 2, "linear")); // 6+2=8
        assert_eq!(model.total_params(), 23);
    }

    #[test]
    fn test_model_weights_to_bytes_nonempty() {
        let mut model = ModelWeights::new();
        model.add_layer(DenseLayer::new("l1", 2, 2, "relu"));
        let bytes = model.to_bytes();
        assert!(!bytes.is_empty());
    }

    // --- StateDict ---

    #[test]
    fn test_state_dict_insert_and_get() {
        let mut sd = StateDict::new();
        sd.insert("w1", Tensor::zeros(vec![4, 4]));
        assert_eq!(sd.len(), 1);
        assert_eq!(sd.get("w1").unwrap().numel(), 16);
    }

    #[test]
    fn test_state_dict_total_params() {
        let mut sd = StateDict::new();
        sd.insert("a", Tensor::zeros(vec![3, 3]));
        sd.insert("b", Tensor::zeros(vec![3]));
        assert_eq!(sd.total_params(), 12);
    }

    #[test]
    fn test_state_dict_roundtrip() {
        let mut sd = StateDict::new();
        sd.insert("w", Tensor::new(vec![2, 2], vec![1.0, 2.0, 3.0, 4.0]));
        sd.insert("b", Tensor::new(vec![2], vec![0.5, -0.5]));
        let bytes = sd.to_bytes();
        let sd2 = StateDict::from_bytes(&bytes).unwrap();
        assert_eq!(sd2.len(), 2);
        let w = sd2.get("w").unwrap();
        assert!((w.data[3] - 4.0).abs() < 1e-12);
    }

    // --- OnnxLikeGraph ---

    #[test]
    fn test_onnx_graph_node_count() {
        let mut g = OnnxLikeGraph::new("test_model");
        g.add_node(OnnxNode::new(
            "n0",
            "MatMul",
            vec!["x".into(), "w0".into()],
            vec!["h0".into()],
        ));
        g.add_node(OnnxNode::new(
            "n1",
            "Relu",
            vec!["h0".into()],
            vec!["h1".into()],
        ));
        assert_eq!(g.node_count(), 2);
        assert_eq!(g.count_op("Relu"), 1);
    }

    #[test]
    fn test_onnx_graph_topological_valid() {
        let mut g = OnnxLikeGraph::new("model");
        g.inputs.push("x".into());
        g.add_initializer("w0", Tensor::zeros(vec![4, 4]));
        g.add_node(OnnxNode::new(
            "mm",
            "MatMul",
            vec!["x".into(), "w0".into()],
            vec!["y".into()],
        ));
        g.add_node(OnnxNode::new(
            "act",
            "Relu",
            vec!["y".into()],
            vec!["z".into()],
        ));
        assert!(g.is_topologically_valid());
    }

    #[test]
    fn test_onnx_graph_topological_invalid() {
        let mut g = OnnxLikeGraph::new("model");
        g.inputs.push("x".into());
        // "undefined" is neither an input nor output of any prior node
        g.add_node(OnnxNode::new(
            "act",
            "Relu",
            vec!["undefined".into()],
            vec!["z".into()],
        ));
        assert!(!g.is_topologically_valid());
    }

    // --- Dataset ---

    #[test]
    fn test_dataset_len_and_features() {
        let mut ds = Dataset::new();
        ds.push(DataRow::labelled(vec![1.0, 2.0], 0));
        ds.push(DataRow::labelled(vec![3.0, 4.0], 1));
        assert_eq!(ds.len(), 2);
        assert_eq!(ds.num_features(), 2);
    }

    #[test]
    fn test_dataset_shuffle_changes_order() {
        let mut ds = Dataset::new();
        for i in 0..20 {
            ds.push(DataRow::labelled(vec![i as f64], 0));
        }
        let original: Vec<f64> = ds.rows.iter().map(|r| r.features[0]).collect();
        ds.shuffle(42);
        let shuffled: Vec<f64> = ds.rows.iter().map(|r| r.features[0]).collect();
        assert_ne!(original, shuffled);
    }

    #[test]
    fn test_dataset_train_val_split() {
        let mut ds = Dataset::new();
        for i in 0..100 {
            ds.push(DataRow::labelled(vec![i as f64], 0));
        }
        let (train, val) = ds.train_val_split(0.2);
        assert_eq!(train.len(), 80);
        assert_eq!(val.len(), 20);
    }

    #[test]
    fn test_dataset_feature_stats() {
        let mut ds = Dataset::new();
        ds.push(DataRow::labelled(vec![0.0, 10.0], 0));
        ds.push(DataRow::labelled(vec![2.0, 10.0], 1));
        let (means, _stds) = ds.feature_stats();
        assert!((means[0] - 1.0).abs() < 1e-12);
        assert!((means[1] - 10.0).abs() < 1e-12);
    }

    // --- NormalizationParams ---

    #[test]
    fn test_normalization_zscore() {
        let mut ds = Dataset::new();
        ds.push(DataRow::labelled(vec![0.0], 0));
        ds.push(DataRow::labelled(vec![2.0], 0));
        let norm = NormalizationParams::from_dataset(&ds);
        let z = norm.normalize_zscore(&[1.0]);
        // (1 - 1) / std = 0
        assert!(z[0].abs() < 1e-10);
    }

    #[test]
    fn test_normalization_minmax() {
        let mut ds = Dataset::new();
        ds.push(DataRow::labelled(vec![0.0], 0));
        ds.push(DataRow::labelled(vec![10.0], 0));
        let norm = NormalizationParams::from_dataset(&ds);
        let v = norm.normalize_minmax(&[5.0]);
        assert!((v[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_normalization_bytes_nonempty() {
        let mut ds = Dataset::new();
        ds.push(DataRow::labelled(vec![1.0, 2.0], 0));
        let norm = NormalizationParams::from_dataset(&ds);
        assert!(!norm.to_bytes().is_empty());
    }

    // --- LabelEncoder ---

    #[test]
    fn test_label_encoder_fit_and_encode() {
        let enc = LabelEncoder::fit(vec!["cat".into(), "dog".into(), "bird".into()]);
        assert_eq!(enc.num_classes(), 3);
        let i = enc.encode("dog").unwrap();
        assert_eq!(enc.decode(i), Some("dog"));
    }

    #[test]
    fn test_label_encoder_one_hot() {
        let enc = LabelEncoder::fit(vec!["a".into(), "b".into(), "c".into()]);
        let oh = enc.one_hot(enc.encode("b").unwrap());
        assert_eq!(oh.iter().filter(|&&v| v == 1.0).count(), 1);
        assert!((oh.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_label_encoder_unknown_returns_none() {
        let enc = LabelEncoder::fit(vec!["a".into()]);
        assert!(enc.encode("z").is_none());
    }

    // --- ConfusionMatrix ---

    #[test]
    fn test_confusion_matrix_accuracy() {
        let mut cm = ConfusionMatrix::new(2);
        cm.record(0, 0);
        cm.record(0, 0);
        cm.record(1, 1);
        cm.record(1, 0); // wrong
        assert!((cm.accuracy() - 0.75).abs() < 1e-12);
    }

    #[test]
    fn test_confusion_matrix_precision_recall() {
        let mut cm = ConfusionMatrix::new(2);
        cm.record(0, 0); // TP for class 0
        cm.record(0, 1); // FN for class 0
        cm.record(1, 0); // FP for class 0
        cm.record(1, 1); // TN for class 0
        let p = cm.precision(0);
        let r = cm.recall(0);
        assert!((p - 0.5).abs() < 1e-12);
        assert!((r - 0.5).abs() < 1e-12);
    }

    #[test]
    fn test_confusion_matrix_to_csv() {
        let mut cm = ConfusionMatrix::new(2);
        cm.record(0, 0);
        cm.record(1, 1);
        let csv = cm.to_csv();
        assert!(csv.contains("class_0"));
        assert!(csv.contains("class_1"));
    }

    // --- TrainingHistory ---

    #[test]
    fn test_training_history_best_val_acc() {
        let mut hist = TrainingHistory::new();
        for e in 0..5 {
            hist.push(EpochRecord {
                epoch: e,
                train_loss: 1.0 - e as f64 * 0.1,
                val_loss: 1.0 - e as f64 * 0.08,
                train_acc: e as f64 * 0.2,
                val_acc: e as f64 * 0.18,
                learning_rate: 0.001,
            });
        }
        let (best_epoch, best_acc) = hist.best_val_acc().unwrap();
        assert_eq!(best_epoch, 4);
        assert!((best_acc - 0.72).abs() < 1e-10);
    }

    #[test]
    fn test_training_history_to_csv() {
        let mut hist = TrainingHistory::new();
        hist.push(EpochRecord {
            epoch: 0,
            train_loss: 0.9,
            val_loss: 0.85,
            train_acc: 0.6,
            val_acc: 0.62,
            learning_rate: 0.01,
        });
        let csv = hist.to_csv();
        assert!(csv.starts_with("epoch,"));
        assert!(csv.contains("0,"));
    }

    // --- HyperparamConfig ---

    #[test]
    fn test_hyperparam_config_get_set() {
        let mut cfg = HyperparamConfig::new();
        cfg.set_float("lr", 0.001);
        cfg.set_bool("dropout", true);
        cfg.set_str("optimizer", "adam");
        assert!((cfg.get_float("lr").unwrap() - 0.001).abs() < 1e-15);
        assert!(cfg.get_bool("dropout").unwrap());
        assert_eq!(cfg.get_str("optimizer").unwrap(), "adam");
    }

    #[test]
    fn test_hyperparam_config_to_json() {
        let mut cfg = HyperparamConfig::new();
        cfg.set_float("lr", 0.01);
        let json = cfg.to_json();
        assert!(json.contains("lr"));
        assert!(json.starts_with('{'));
        assert!(json.ends_with('}'));
    }

    // --- ModelCheckpoint ---

    #[test]
    fn test_checkpoint_byte_size_nonzero() {
        let state = StateDict::new();
        let meta = CheckpointMeta {
            epoch: 10,
            val_loss: 0.1,
            val_acc: 0.95,
            train_time_secs: 3600.0,
            architecture: "MLP".into(),
            framework_version: "0.1.0".into(),
        };
        let hparams = HyperparamConfig::new();
        let ck = ModelCheckpoint::new(state, meta, hparams);
        assert!(ck.byte_size() > 0);
    }

    #[test]
    fn test_checkpoint_meta_to_text_contains_epoch() {
        let meta = CheckpointMeta {
            epoch: 42,
            val_loss: 0.05,
            val_acc: 0.98,
            train_time_secs: 100.0,
            architecture: "CNN".into(),
            framework_version: "0.1.0".into(),
        };
        let text = meta.to_text();
        assert!(text.contains("epoch=42"));
    }

    // --- utility functions ---

    #[test]
    fn test_softmax_sums_to_one() {
        let logits = vec![1.0, 2.0, 3.0];
        let probs = softmax(&logits);
        let total: f64 = probs.iter().sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_softmax_max_has_highest_prob() {
        let logits = vec![1.0, 5.0, 2.0];
        let probs = softmax(&logits);
        assert!(probs[1] > probs[0] && probs[1] > probs[2]);
    }

    #[test]
    fn test_cross_entropy_perfect_prediction() {
        let probs = vec![0.0, 1.0, 0.0];
        let targets = vec![0.0, 1.0, 0.0];
        let loss = cross_entropy_loss(&probs, &targets);
        assert!(loss < 1e-10);
    }

    #[test]
    fn test_argmax_basic() {
        let v = vec![0.1, 0.7, 0.2];
        assert_eq!(argmax(&v), 1);
    }

    #[test]
    fn test_mse_zero() {
        let p = vec![1.0, 2.0, 3.0];
        let t = vec![1.0, 2.0, 3.0];
        assert!(mse(&p, &t).abs() < 1e-12);
    }

    #[test]
    fn test_mse_known() {
        let p = vec![0.0, 0.0];
        let t = vec![1.0, 1.0];
        assert!((mse(&p, &t) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_mae_basic() {
        let p = vec![0.0, 1.0, 2.0];
        let t = vec![1.0, 1.0, 3.0];
        // |0-1| + |1-1| + |2-3| = 1+0+1 = 2, /3 = 0.666...
        let m = mae(&p, &t);
        assert!((m - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_apply_activation_leaky_relu() {
        assert!((apply_activation(-1.0, "leaky_relu") - (-0.01)).abs() < 1e-12);
        assert!((apply_activation(2.0, "leaky_relu") - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_apply_activation_elu() {
        let v = apply_activation(-1.0, "elu");
        // elu(-1) = e^(-1) - 1 ≈ -0.6321
        assert!(v < 0.0 && v > -1.0);
    }

    #[test]
    fn test_lcg_rng_produces_different_values() {
        let mut rng = LcgRng::new(1234);
        let a = rng.next_u64();
        let b = rng.next_u64();
        assert_ne!(a, b);
    }
}
