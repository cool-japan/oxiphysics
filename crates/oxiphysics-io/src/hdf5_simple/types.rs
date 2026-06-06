// Auto-generated module
//
// 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use super::functions::*;

/// XDMF topology type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum XdmfTopologyType {
    /// Triangle elements.
    Triangle,
    /// Tetrahedral elements.
    Tetrahedron,
    /// Hexahedral elements.
    Hexahedron,
    /// Unstructured elements (mixed).
    Mixed,
}
impl XdmfTopologyType {
    pub(super) fn as_str(self) -> &'static str {
        match self {
            Self::Triangle => "Triangle",
            Self::Tetrahedron => "Tetrahedron",
            Self::Hexahedron => "Hexahedron",
            Self::Mixed => "Mixed",
        }
    }
}
/// Metadata describing compression settings for a dataset.
#[derive(Debug, Clone)]
pub struct DeflateMetadata {
    /// Compression level.
    pub level: CompressionLevel,
    /// Whether shuffle filter is applied before compression.
    pub shuffle: bool,
    /// Chunk shape used during compression.
    pub chunk_shape: Vec<u64>,
    /// Compressed size in bytes (0 if not yet written).
    pub compressed_size: u64,
    /// Uncompressed size in bytes.
    pub uncompressed_size: u64,
}
impl DeflateMetadata {
    /// Create metadata for an uncompressed dataset.
    pub fn uncompressed(uncompressed_size: u64) -> Self {
        Self {
            level: CompressionLevel::None,
            shuffle: false,
            chunk_shape: Vec::new(),
            compressed_size: uncompressed_size,
            uncompressed_size,
        }
    }
    /// Estimated compression ratio (uncompressed / compressed).
    pub fn compression_ratio(&self) -> f64 {
        if self.compressed_size == 0 {
            return 1.0;
        }
        self.uncompressed_size as f64 / self.compressed_size as f64
    }
    /// Space savings fraction \[0, 1).
    pub fn space_savings(&self) -> f64 {
        if self.uncompressed_size == 0 {
            return 0.0;
        }
        1.0 - self.compressed_size as f64 / self.uncompressed_size as f64
    }
}
/// Compression level for SHDF datasets (analogous to HDF5 deflate filter).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompressionLevel {
    /// No compression.
    None,
    /// Fast compression (level 1).
    Fast,
    /// Balanced compression (level 5).
    Balanced,
    /// Maximum compression (level 9).
    Maximum,
}
impl CompressionLevel {
    /// Return the numeric compression level.
    pub fn level(self) -> u8 {
        match self {
            Self::None => 0,
            Self::Fast => 1,
            Self::Balanced => 5,
            Self::Maximum => 9,
        }
    }
    /// Whether compression is active.
    pub fn is_compressed(self) -> bool {
        !matches!(self, Self::None)
    }
}
/// A single named, shaped dataset inside a [`ShdfFile`].
#[derive(Debug, Clone)]
pub struct Dataset {
    /// Dataset name (must be unique within a file).
    pub name: String,
    /// Shape (dimension sizes), e.g. `[rows, cols]`.
    pub shape: Vec<usize>,
    /// Declared element type.
    pub dtype: DataType,
    /// Floating-point payload (used when `dtype` is Float64 or Float32).
    pub data_f64: Vec<f64>,
    /// Integer payload (used when `dtype` is Int32 or Int64).
    pub data_i32: Vec<i32>,
    /// Per-dataset key/value string attributes.
    pub attributes: Vec<(String, String)>,
}
/// Compression settings for a dataset.
#[derive(Debug, Clone)]
pub struct CompressionSettings {
    /// Which algorithm to use.
    pub algorithm: CompressionAlgorithm,
    /// Compression level (0-9, higher = more compression).
    pub level: u32,
}
impl CompressionSettings {
    /// No compression.
    pub fn none() -> Self {
        Self {
            algorithm: CompressionAlgorithm::None,
            level: 0,
        }
    }
    /// Delta encoding.
    pub fn delta() -> Self {
        Self {
            algorithm: CompressionAlgorithm::Delta,
            level: 1,
        }
    }
    /// Apply delta encoding to f64 data.
    pub fn delta_encode_f64(data: &[f64]) -> Vec<f64> {
        if data.is_empty() {
            return Vec::new();
        }
        let mut encoded = Vec::with_capacity(data.len());
        encoded.push(data[0]);
        for i in 1..data.len() {
            encoded.push(data[i] - data[i - 1]);
        }
        encoded
    }
    /// Decode delta-encoded f64 data.
    pub fn delta_decode_f64(encoded: &[f64]) -> Vec<f64> {
        if encoded.is_empty() {
            return Vec::new();
        }
        let mut decoded = Vec::with_capacity(encoded.len());
        decoded.push(encoded[0]);
        for i in 1..encoded.len() {
            decoded.push(decoded[i - 1] + encoded[i]);
        }
        decoded
    }
    /// Apply delta encoding to i32 data.
    pub fn delta_encode_i32(data: &[i32]) -> Vec<i32> {
        if data.is_empty() {
            return Vec::new();
        }
        let mut encoded = Vec::with_capacity(data.len());
        encoded.push(data[0]);
        for i in 1..data.len() {
            encoded.push(data[i] - data[i - 1]);
        }
        encoded
    }
    /// Decode delta-encoded i32 data.
    pub fn delta_decode_i32(encoded: &[i32]) -> Vec<i32> {
        if encoded.is_empty() {
            return Vec::new();
        }
        let mut decoded = Vec::with_capacity(encoded.len());
        decoded.push(encoded[0]);
        for i in 1..encoded.len() {
            decoded.push(decoded[i - 1] + encoded[i]);
        }
        decoded
    }
}
/// A field in a compound dataset (analogous to an HDF5 compound type member).
#[derive(Debug, Clone)]
pub struct CompoundField {
    /// Field name.
    pub name: String,
    /// Data type.
    pub dtype: DataType,
    /// Field values (all records, serialized as f64 regardless of dtype for simplicity).
    pub values: Vec<f64>,
}
impl CompoundField {
    /// Create a new compound field.
    pub fn new(name: impl Into<String>, dtype: DataType, values: Vec<f64>) -> Self {
        Self {
            name: name.into(),
            dtype,
            values,
        }
    }
}
/// Compression algorithm selection.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CompressionAlgorithm {
    /// No compression.
    None,
    /// Run-length encoding (simple).
    RunLength,
    /// Delta encoding (store differences).
    Delta,
}
/// A compound dataset (analogous to HDF5 compound dataset with multiple fields per record).
#[derive(Debug, Clone)]
pub struct CompoundDataset {
    /// Dataset name.
    pub name: String,
    /// Number of records.
    pub n_records: usize,
    /// Fields (each has `n_records` values).
    pub fields: Vec<CompoundField>,
    /// Attributes.
    pub attrs: Vec<(String, String)>,
}
impl CompoundDataset {
    /// Create a new compound dataset.
    pub fn new(name: impl Into<String>, n_records: usize) -> Self {
        Self {
            name: name.into(),
            n_records,
            fields: Vec::new(),
            attrs: Vec::new(),
        }
    }
    /// Add a field.
    pub fn add_field(&mut self, field: CompoundField) {
        assert_eq!(
            field.values.len(),
            self.n_records,
            "Field {} has {} values, expected {}",
            field.name,
            field.values.len(),
            self.n_records
        );
        self.fields.push(field);
    }
    /// Add an attribute.
    pub fn add_attr(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.attrs.push((key.into(), value.into()));
    }
    /// Get values for a named field.
    pub fn get_field(&self, name: &str) -> Option<&[f64]> {
        self.fields
            .iter()
            .find(|f| f.name == name)
            .map(|f| f.values.as_slice())
    }
    /// Number of fields.
    pub fn n_fields(&self) -> usize {
        self.fields.len()
    }
    /// Serialize the compound dataset to a flat CSV-like byte buffer (for debugging).
    pub fn to_csv_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let header: Vec<&str> = self.fields.iter().map(|f| f.name.as_str()).collect();
        out.extend_from_slice(header.join(",").as_bytes());
        out.push(b'\n');
        for rec in 0..self.n_records {
            let vals: Vec<String> = self
                .fields
                .iter()
                .map(|f| format!("{:.6}", f.values[rec]))
                .collect();
            out.extend_from_slice(vals.join(",").as_bytes());
            out.push(b'\n');
        }
        out
    }
}
/// Basic statistics computed over a [`Dataset`]'s raw f64 data.
#[derive(Debug, Clone)]
pub struct DatasetStats {
    /// Number of elements.
    pub count: usize,
    /// Minimum value.
    pub min: f64,
    /// Maximum value.
    pub max: f64,
    /// Arithmetic mean.
    pub mean: f64,
    /// Variance.
    pub variance: f64,
}
impl DatasetStats {
    /// Compute statistics from a slice of f64 values.
    pub fn from_slice(data: &[f64]) -> Option<Self> {
        if data.is_empty() {
            return None;
        }
        let count = data.len();
        let mut min = data[0];
        let mut max = data[0];
        let mut sum = 0.0_f64;
        for &v in data {
            if v < min {
                min = v;
            }
            if v > max {
                max = v;
            }
            sum += v;
        }
        let mean = sum / count as f64;
        let variance = data.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / count as f64;
        Some(Self {
            count,
            min,
            max,
            mean,
            variance,
        })
    }
    /// Standard deviation (sqrt of variance).
    pub fn std_dev(&self) -> f64 {
        self.variance.sqrt()
    }
    /// Range: max - min.
    pub fn range(&self) -> f64 {
        self.max - self.min
    }
}
/// Expected schema for an SHDF file.
#[derive(Debug, Clone)]
pub struct ShdfSchema {
    /// Expected dataset names and their types.
    pub expected_datasets: Vec<(String, DataType)>,
    /// Required global attributes.
    pub required_attributes: Vec<String>,
}
impl ShdfSchema {
    /// Create a new schema.
    pub fn new() -> Self {
        Self {
            expected_datasets: Vec::new(),
            required_attributes: Vec::new(),
        }
    }
    /// Add an expected dataset.
    pub fn expect_dataset(&mut self, name: &str, dtype: DataType) {
        self.expected_datasets.push((name.to_string(), dtype));
    }
    /// Add a required global attribute.
    pub fn require_attribute(&mut self, key: &str) {
        self.required_attributes.push(key.to_string());
    }
    /// Validate an SHDF file against this schema.
    ///
    /// Returns a list of validation errors (empty = valid).
    pub fn validate(&self, file: &ShdfFile) -> Vec<String> {
        let mut errors = Vec::new();
        for (name, dtype) in &self.expected_datasets {
            match file.datasets.iter().find(|d| &d.name == name) {
                None => errors.push(format!("Missing dataset: {name}")),
                Some(ds) => {
                    if ds.dtype != *dtype {
                        errors.push(format!(
                            "Dataset '{name}': expected {:?}, got {:?}",
                            dtype, ds.dtype
                        ));
                    }
                }
            }
        }
        for key in &self.required_attributes {
            if !file.global_attributes.iter().any(|(k, _)| k == key) {
                errors.push(format!("Missing global attribute: {key}"));
            }
        }
        errors
    }
}
impl Default for ShdfSchema {
    fn default() -> Self {
        Self::new()
    }
}
/// Chunk descriptor: defines chunk shape and offset within a dataset.
#[derive(Debug, Clone)]
pub struct ChunkDescriptor {
    /// Chunk dimensions (same rank as dataset).
    pub shape: Vec<u64>,
    /// Offset of this chunk in the dataset (per dimension).
    pub offset: Vec<u64>,
    /// Flat index of this chunk.
    pub index: u64,
}
/// Data type tag stored in each [`Dataset`].
#[derive(Debug, Clone, PartialEq)]
pub enum DataType {
    /// 64-bit IEEE 754 floating point.
    Float64,
    /// 32-bit IEEE 754 floating point.
    Float32,
    /// 32-bit signed integer.
    Int32,
    /// 64-bit signed integer.
    Int64,
}
/// A chunked dataset: stores data split into uniform-sized chunks.
#[derive(Debug, Clone)]
pub struct ChunkedDataset {
    /// Dataset name.
    pub name: String,
    /// Full dataset dimensions.
    pub dims: Vec<u64>,
    /// Chunk shape.
    pub chunk_shape: Vec<u64>,
    /// Data (f64), stored flat across all chunks.
    pub data: Vec<f64>,
    /// Attributes.
    pub attrs: Vec<(String, String)>,
}
impl ChunkedDataset {
    /// Create a new chunked dataset.
    pub fn new(name: impl Into<String>, dims: Vec<u64>, chunk_shape: Vec<u64>) -> Self {
        let total: u64 = dims.iter().product();
        Self {
            name: name.into(),
            dims,
            chunk_shape,
            data: vec![0.0; total as usize],
            attrs: Vec::new(),
        }
    }
    /// Total number of elements.
    pub fn n_elements(&self) -> usize {
        self.dims.iter().product::<u64>() as usize
    }
    /// Compute number of chunks per dimension.
    pub fn n_chunks_per_dim(&self) -> Vec<u64> {
        self.dims
            .iter()
            .zip(self.chunk_shape.iter())
            .map(|(&d, &c)| d.div_ceil(c))
            .collect()
    }
    /// Total number of chunks.
    pub fn total_chunks(&self) -> u64 {
        self.n_chunks_per_dim().iter().product()
    }
    /// Write data for a specific 1-D chunk (row-major slice of `data`).
    pub fn write_chunk_1d(&mut self, chunk_idx: u64, chunk_data: &[f64]) {
        let chunk_size = self.chunk_shape[0] as usize;
        let start = (chunk_idx as usize) * chunk_size;
        let end = (start + chunk_data.len()).min(self.data.len());
        let src_end = end - start;
        self.data[start..end].copy_from_slice(&chunk_data[..src_end]);
    }
    /// Add an attribute.
    pub fn add_attr(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.attrs.push((key.into(), value.into()));
    }
    /// Serialize the chunked dataset to bytes (SHDF-compatible format).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        let name_bytes = self.name.as_bytes();
        buf.extend_from_slice(&(name_bytes.len() as u32).to_le_bytes());
        buf.extend_from_slice(name_bytes);
        buf.push(0u8);
        buf.extend_from_slice(&(self.dims.len() as u32).to_le_bytes());
        for &d in &self.dims {
            buf.extend_from_slice(&d.to_le_bytes());
        }
        let n = self.n_elements() as u64;
        buf.extend_from_slice(&n.to_le_bytes());
        for &v in &self.data {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf.extend_from_slice(&(self.attrs.len() as u32).to_le_bytes());
        for (k, v) in &self.attrs {
            let kb = k.as_bytes();
            buf.extend_from_slice(&(kb.len() as u32).to_le_bytes());
            buf.extend_from_slice(kb);
            let vb = v.as_bytes();
            buf.extend_from_slice(&(vb.len() as u32).to_le_bytes());
            buf.extend_from_slice(vb);
        }
        buf
    }
}
/// Parameters for writing an XDMF file that references HDF5/SHDF data.
#[derive(Debug, Clone)]
pub struct XdmfParams {
    /// Path to the HDF5/SHDF file.
    pub hdf5_path: String,
    /// Dataset name for coordinates (within the HDF5 file).
    pub coords_dataset: String,
    /// Dataset name for connectivity (within the HDF5 file).
    pub connectivity_dataset: String,
    /// Number of nodes.
    pub n_nodes: usize,
    /// Number of elements.
    pub n_elements: usize,
    /// Nodes per element.
    pub nodes_per_element: usize,
    /// Topology type.
    pub topology: XdmfTopologyType,
    /// Named attribute datasets (attribute_name → dataset_path).
    pub attributes: Vec<(String, String)>,
}
/// Navigate an [`ShdfGroup`] hierarchy using HDF5-style slash-separated paths.
///
/// # Example
/// ```no_run
/// # use oxiphysics_io::hdf5_simple::*;
/// let mut root = ShdfGroup::new("root");
/// let mut sim  = ShdfGroup::new("simulation");
/// sim.add_dataset_f64("time", vec![3], vec![0.0, 0.5, 1.0]);
/// root.add_child(sim);
/// let nav = GroupNavigator::new(root);
/// assert!(nav.get_dataset("/simulation/time").is_some());
/// ```
pub struct GroupNavigator {
    /// Root group of the hierarchy.
    pub root: ShdfGroup,
}
impl GroupNavigator {
    /// Create a navigator wrapping a root group.
    pub fn new(root: ShdfGroup) -> Self {
        Self { root }
    }
    /// Resolve a slash-separated path such as `"/root/simulation/atoms/positions"`
    /// and return the terminal dataset if it exists.
    ///
    /// The first path component must match the root group name.
    pub fn get_dataset(&self, path: &str) -> Option<&Dataset> {
        let parts: Vec<&str> = path.trim_start_matches('/').splitn(64, '/').collect();
        if parts.is_empty() {
            return None;
        }
        let (ds_name, group_parts) = parts.split_last()?;
        let mut group = &self.root;
        let effective_parts = if group_parts.first().copied() == Some(self.root.name.as_str()) {
            &group_parts[1..]
        } else {
            group_parts
        };
        for &part in effective_parts {
            group = group.get_child(part)?;
        }
        group.get_dataset(ds_name)
    }
    /// Return all dataset paths reachable from the root, in DFS order.
    pub fn all_paths(&self) -> Vec<String> {
        let mut result = Vec::new();
        Self::collect_paths(&self.root, "", &mut result);
        result
    }
    fn collect_paths(group: &ShdfGroup, prefix: &str, out: &mut Vec<String>) {
        let base = if prefix.is_empty() {
            format!("/{}", group.name)
        } else {
            format!("{}/{}", prefix, group.name)
        };
        for ds in &group.datasets {
            out.push(format!("{}/{}", base, ds.name));
        }
        for child in &group.children {
            Self::collect_paths(child, &base, out);
        }
    }
    /// Count total datasets reachable from the root.
    pub fn total_datasets(&self) -> usize {
        self.root.total_datasets()
    }
}
/// Chunking configuration for a dataset.
#[derive(Debug, Clone)]
pub struct ChunkingConfig {
    /// Chunk dimensions. Must have the same number of dimensions as the dataset.
    pub chunk_dims: Vec<usize>,
}
impl ChunkingConfig {
    /// Create a chunking config with the given chunk dimensions.
    pub fn new(chunk_dims: Vec<usize>) -> Self {
        Self { chunk_dims }
    }
    /// Compute the number of chunks needed for a dataset with the given shape.
    pub fn n_chunks(&self, shape: &[usize]) -> usize {
        if shape.len() != self.chunk_dims.len() {
            return 0;
        }
        let mut total = 1_usize;
        for (s, c) in shape.iter().zip(self.chunk_dims.iter()) {
            if *c == 0 {
                return 0;
            }
            total *= (*s).div_ceil(*c);
        }
        total
    }
    /// Compute the linear index of the chunk containing the given element.
    pub fn chunk_index(&self, element_idx: &[usize], shape: &[usize]) -> usize {
        if shape.len() != self.chunk_dims.len() || element_idx.len() != shape.len() {
            return 0;
        }
        let mut idx = 0;
        let mut stride = 1;
        for d in (0..shape.len()).rev() {
            let chunk_pos = element_idx[d] / self.chunk_dims[d].max(1);
            let n_chunks_d = (shape[d] + self.chunk_dims[d] - 1) / self.chunk_dims[d].max(1);
            idx += chunk_pos * stride;
            stride *= n_chunks_d;
        }
        idx
    }
    /// Default chunking: chunk size of 64 in each dimension.
    pub fn default_for_shape(shape: &[usize]) -> Self {
        let chunk_dims: Vec<usize> = shape.iter().map(|&s| s.min(64)).collect();
        Self { chunk_dims }
    }
}
/// Extended attribute support with typed values.
#[derive(Debug, Clone, PartialEq)]
pub enum AttributeValue {
    /// String value.
    String(String),
    /// 64-bit float.
    Float64(f64),
    /// 32-bit integer.
    Int32(i32),
    /// Boolean.
    Bool(bool),
}
/// A group in the SHDF hierarchy, containing datasets and sub-groups.
#[derive(Debug, Clone)]
pub struct ShdfGroup {
    /// Group name.
    pub name: String,
    /// Datasets in this group.
    pub datasets: Vec<Dataset>,
    /// Sub-groups.
    pub children: Vec<ShdfGroup>,
    /// Group-level attributes.
    pub attributes: Vec<(String, String)>,
}
impl ShdfGroup {
    /// Create a new empty group.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            datasets: Vec::new(),
            children: Vec::new(),
            attributes: Vec::new(),
        }
    }
    /// Add a Float64 dataset to this group.
    pub fn add_dataset_f64(&mut self, name: &str, shape: Vec<usize>, data: Vec<f64>) {
        self.datasets.push(Dataset {
            name: name.to_string(),
            shape,
            dtype: DataType::Float64,
            data_f64: data,
            data_i32: Vec::new(),
            attributes: Vec::new(),
        });
    }
    /// Add an Int32 dataset to this group.
    pub fn add_dataset_i32(&mut self, name: &str, shape: Vec<usize>, data: Vec<i32>) {
        self.datasets.push(Dataset {
            name: name.to_string(),
            shape,
            dtype: DataType::Int32,
            data_f64: Vec::new(),
            data_i32: data,
            attributes: Vec::new(),
        });
    }
    /// Add a child group.
    pub fn add_child(&mut self, child: ShdfGroup) {
        self.children.push(child);
    }
    /// Add an attribute to this group.
    pub fn add_attribute(&mut self, key: &str, value: &str) {
        self.attributes.push((key.to_string(), value.to_string()));
    }
    /// Find a dataset by name.
    pub fn get_dataset(&self, name: &str) -> Option<&Dataset> {
        self.datasets.iter().find(|d| d.name == name)
    }
    /// Find a child group by name.
    pub fn get_child(&self, name: &str) -> Option<&ShdfGroup> {
        self.children.iter().find(|c| c.name == name)
    }
    /// Count total datasets (recursive).
    pub fn total_datasets(&self) -> usize {
        self.datasets.len()
            + self
                .children
                .iter()
                .map(|c| c.total_datasets())
                .sum::<usize>()
    }
    /// Generate a text summary of this group hierarchy.
    pub fn summary(&self, indent: usize) -> String {
        let prefix = " ".repeat(indent);
        let mut out = format!("{prefix}Group: {}\n", self.name);
        for (k, v) in &self.attributes {
            out.push_str(&format!("{prefix}  attr: {k} = {v}\n"));
        }
        for ds in &self.datasets {
            let shape_str: Vec<String> = ds.shape.iter().map(|s| s.to_string()).collect();
            out.push_str(&format!(
                "{prefix}  Dataset: {} shape=[{}] dtype={:?}\n",
                ds.name,
                shape_str.join("x"),
                ds.dtype,
            ));
        }
        for child in &self.children {
            out.push_str(&child.summary(indent + 2));
        }
        out
    }
}
/// Append new 1-D f64 frames to a growing time-series dataset.
///
/// Each call to [`TimeSeriesAppender::append`] concatenates data to the
/// internal buffer, tracking the number of appended frames.
pub struct TimeSeriesAppender {
    /// Name of the logical dataset.
    pub name: String,
    /// Accumulated sample data.
    pub data: Vec<f64>,
    /// Number of appended frames.
    pub n_frames: usize,
    /// Number of values per frame (fixed at construction time).
    pub frame_width: usize,
}
impl TimeSeriesAppender {
    /// Create a new appender.
    ///
    /// `frame_width` is the number of f64 values in each frame.
    pub fn new(name: impl Into<String>, frame_width: usize) -> Self {
        Self {
            name: name.into(),
            data: Vec::new(),
            n_frames: 0,
            frame_width,
        }
    }
    /// Append a single frame's worth of data.
    ///
    /// # Panics
    /// Panics if `frame.len() != self.frame_width`.
    pub fn append(&mut self, frame: &[f64]) {
        assert_eq!(
            frame.len(),
            self.frame_width,
            "frame length {} != frame_width {}",
            frame.len(),
            self.frame_width
        );
        self.data.extend_from_slice(frame);
        self.n_frames += 1;
    }
    /// Return the total number of f64 samples stored.
    pub fn total_samples(&self) -> usize {
        self.data.len()
    }
    /// Retrieve a specific frame by index (0-based).
    pub fn get_frame(&self, idx: usize) -> Option<&[f64]> {
        let start = idx * self.frame_width;
        let end = start + self.frame_width;
        self.data.get(start..end)
    }
    /// Export the accumulated data as a 2-D [`Dataset`] (shape `[n_frames, frame_width]`).
    pub fn to_dataset(&self) -> Dataset {
        let mut ds = Dataset {
            name: self.name.clone(),
            dtype: DataType::Float64,
            shape: vec![self.n_frames, self.frame_width],
            data_f64: self.data.clone(),
            data_i32: Vec::new(),
            attributes: Vec::new(),
        };
        ds.attributes
            .push(("n_frames".to_string(), self.n_frames.to_string()));
        ds.attributes
            .push(("frame_width".to_string(), self.frame_width.to_string()));
        ds
    }
}
/// Describes a virtual link from one dataset path to another — similar to
/// HDF5 virtual datasets or external links.
#[derive(Debug, Clone)]
pub struct VirtualLink {
    /// The logical path within this file (e.g. `"/virtual/positions"`).
    pub virtual_path: String,
    /// The source file (may be the same file or an external path).
    pub source_file: String,
    /// The dataset path inside the source file.
    pub source_path: String,
    /// Optional slice: `[start, stop, step]` for each dimension.
    pub slices: Vec<[usize; 3]>,
}
impl VirtualLink {
    /// Create a new virtual link with no slicing.
    pub fn new(
        virtual_path: impl Into<String>,
        source_file: impl Into<String>,
        source_path: impl Into<String>,
    ) -> Self {
        Self {
            virtual_path: virtual_path.into(),
            source_file: source_file.into(),
            source_path: source_path.into(),
            slices: Vec::new(),
        }
    }
    /// Add a per-dimension slice `[start, stop, step]`.
    pub fn with_slice(mut self, start: usize, stop: usize, step: usize) -> Self {
        self.slices.push([start, stop, step]);
        self
    }
    /// Return a CDL-style description of this virtual link.
    pub fn to_cdl(&self) -> String {
        let slice_str = if self.slices.is_empty() {
            "(:)".to_string()
        } else {
            let parts: Vec<String> = self
                .slices
                .iter()
                .map(|s| format!("{}:{}:{}", s[0], s[1], s[2]))
                .collect();
            format!("({})", parts.join(", "))
        };
        format!(
            "{} -> {}:{}{}",
            self.virtual_path, self.source_file, self.source_path, slice_str
        )
    }
}
/// An in-memory `.shdf` file.
#[derive(Debug, Clone)]
pub struct ShdfFile {
    /// Ordered list of datasets.
    pub datasets: Vec<Dataset>,
    /// File-level key/value string attributes.
    pub global_attributes: Vec<(String, String)>,
}
impl ShdfFile {
    /// Create an empty [`ShdfFile`].
    pub fn new() -> Self {
        ShdfFile {
            datasets: Vec::new(),
            global_attributes: Vec::new(),
        }
    }
    /// Append a Float64 dataset.
    pub fn add_dataset_f64(&mut self, name: &str, shape: Vec<usize>, data: Vec<f64>) {
        self.datasets.push(Dataset {
            name: name.to_string(),
            shape,
            dtype: DataType::Float64,
            data_f64: data,
            data_i32: Vec::new(),
            attributes: Vec::new(),
        });
    }
    /// Append an Int32 dataset.
    pub fn add_dataset_i32(&mut self, name: &str, shape: Vec<usize>, data: Vec<i32>) {
        self.datasets.push(Dataset {
            name: name.to_string(),
            shape,
            dtype: DataType::Int32,
            data_f64: Vec::new(),
            data_i32: data,
            attributes: Vec::new(),
        });
    }
    /// Add a global (file-level) key/value attribute.
    pub fn add_global_attr(&mut self, key: &str, value: &str) {
        self.global_attributes
            .push((key.to_string(), value.to_string()));
    }
    /// Look up a Float64 dataset by name and return its data slice.
    pub fn get_f64(&self, name: &str) -> Option<&[f64]> {
        self.datasets
            .iter()
            .find(|d| d.name == name)
            .map(|d| d.data_f64.as_slice())
    }
    /// Look up an Int32 dataset by name and return its data slice.
    pub fn get_i32(&self, name: &str) -> Option<&[i32]> {
        self.datasets
            .iter()
            .find(|d| d.name == name)
            .map(|d| d.data_i32.as_slice())
    }
    /// Serialize the file to a binary blob (little-endian).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out: Vec<u8> = Vec::new();
        out.extend_from_slice(MAGIC);
        out.extend_from_slice(&VERSION.to_le_bytes());
        out.extend_from_slice(&(self.global_attributes.len() as u32).to_le_bytes());
        for (k, v) in &self.global_attributes {
            out.extend_from_slice(&encode_string(k));
            out.extend_from_slice(&encode_string(v));
        }
        out.extend_from_slice(&(self.datasets.len() as u32).to_le_bytes());
        for ds in &self.datasets {
            out.extend_from_slice(&encode_string(&ds.name));
            let dtype_byte: u8 = match ds.dtype {
                DataType::Float64 => 0,
                DataType::Float32 => 1,
                DataType::Int32 => 2,
                DataType::Int64 => 3,
            };
            out.push(dtype_byte);
            out.extend_from_slice(&(ds.shape.len() as u32).to_le_bytes());
            for &dim in &ds.shape {
                out.extend_from_slice(&(dim as u64).to_le_bytes());
            }
            match ds.dtype {
                DataType::Float64 => {
                    out.extend_from_slice(&(ds.data_f64.len() as u64).to_le_bytes());
                    for &v in &ds.data_f64 {
                        out.extend_from_slice(&v.to_le_bytes());
                    }
                }
                DataType::Float32 => {
                    out.extend_from_slice(&(ds.data_f64.len() as u64).to_le_bytes());
                    for &v in &ds.data_f64 {
                        out.extend_from_slice(&(v as f32).to_le_bytes());
                    }
                }
                DataType::Int32 => {
                    out.extend_from_slice(&(ds.data_i32.len() as u64).to_le_bytes());
                    for &v in &ds.data_i32 {
                        out.extend_from_slice(&v.to_le_bytes());
                    }
                }
                DataType::Int64 => {
                    out.extend_from_slice(&(ds.data_i32.len() as u64).to_le_bytes());
                    for &v in &ds.data_i32 {
                        out.extend_from_slice(&(v as i64).to_le_bytes());
                    }
                }
            }
            out.extend_from_slice(&(ds.attributes.len() as u32).to_le_bytes());
            for (k, v) in &ds.attributes {
                out.extend_from_slice(&encode_string(k));
                out.extend_from_slice(&encode_string(v));
            }
        }
        out
    }
    /// Deserialize a binary blob written by [`ShdfFile::to_bytes`].
    ///
    /// Returns `Err(String)` on any format violation.
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        let mut pos: usize = 0;
        if data.len() < 4 {
            return Err("too short for magic".to_string());
        }
        if &data[pos..pos + 4] != MAGIC {
            return Err(format!("bad magic: {:?}", &data[pos..pos + 4]));
        }
        pos += 4;
        let version = read_u32(data, &mut pos)?;
        if version != VERSION {
            return Err(format!("unsupported version: {version}"));
        }
        let n_global = read_u32(data, &mut pos)? as usize;
        let mut global_attributes = Vec::with_capacity(n_global);
        for _ in 0..n_global {
            let k = decode_string(data, &mut pos)?;
            let v = decode_string(data, &mut pos)?;
            global_attributes.push((k, v));
        }
        let n_datasets = read_u32(data, &mut pos)? as usize;
        let mut datasets = Vec::with_capacity(n_datasets);
        for _ in 0..n_datasets {
            let name = decode_string(data, &mut pos)?;
            let dtype_byte = read_u8(data, &mut pos)?;
            let dtype = match dtype_byte {
                0 => DataType::Float64,
                1 => DataType::Float32,
                2 => DataType::Int32,
                3 => DataType::Int64,
                _ => return Err(format!("unknown dtype byte: {dtype_byte}")),
            };
            let n_dims = read_u32(data, &mut pos)? as usize;
            let mut shape = Vec::with_capacity(n_dims);
            for _ in 0..n_dims {
                shape.push(read_u64(data, &mut pos)? as usize);
            }
            let n_elems = read_u64(data, &mut pos)? as usize;
            let mut data_f64 = Vec::new();
            let mut data_i32 = Vec::new();
            match dtype {
                DataType::Float64 => {
                    data_f64.reserve(n_elems);
                    for _ in 0..n_elems {
                        data_f64.push(read_f64(data, &mut pos)?);
                    }
                }
                DataType::Float32 => {
                    data_f64.reserve(n_elems);
                    for _ in 0..n_elems {
                        data_f64.push(read_f32(data, &mut pos)? as f64);
                    }
                }
                DataType::Int32 => {
                    data_i32.reserve(n_elems);
                    for _ in 0..n_elems {
                        data_i32.push(read_i32(data, &mut pos)?);
                    }
                }
                DataType::Int64 => {
                    data_i32.reserve(n_elems);
                    for _ in 0..n_elems {
                        data_i32.push(read_i64(data, &mut pos)? as i32);
                    }
                }
            }
            let n_attrs = read_u32(data, &mut pos)? as usize;
            let mut attributes = Vec::with_capacity(n_attrs);
            for _ in 0..n_attrs {
                let k = decode_string(data, &mut pos)?;
                let v = decode_string(data, &mut pos)?;
                attributes.push((k, v));
            }
            datasets.push(Dataset {
                name,
                shape,
                dtype,
                data_f64,
                data_i32,
                attributes,
            });
        }
        Ok(ShdfFile {
            datasets,
            global_attributes,
        })
    }
    /// Return a human-readable summary of the file contents.
    pub fn write_to_text(&self) -> String {
        let mut out = String::new();
        out.push_str("=== SHDF File Summary ===\n");
        if !self.global_attributes.is_empty() {
            out.push_str("Global attributes:\n");
            for (k, v) in &self.global_attributes {
                out.push_str(&format!("  {k} = {v}\n"));
            }
        }
        out.push_str(&format!("Datasets: {}\n", self.datasets.len()));
        for ds in &self.datasets {
            let shape_str: Vec<String> = ds.shape.iter().map(|s| s.to_string()).collect();
            out.push_str(&format!(
                "  [{}] shape=[{}] dtype={:?}\n",
                ds.name,
                shape_str.join("×"),
                ds.dtype,
            ));
            let preview = match ds.dtype {
                DataType::Float64 | DataType::Float32 => {
                    let vals: Vec<String> = ds
                        .data_f64
                        .iter()
                        .take(5)
                        .map(|v| format!("{v:.6}"))
                        .collect();
                    vals.join(", ")
                }
                DataType::Int32 | DataType::Int64 => {
                    let vals: Vec<String> =
                        ds.data_i32.iter().take(5).map(|v| v.to_string()).collect();
                    vals.join(", ")
                }
            };
            if !preview.is_empty() {
                out.push_str(&format!("    first values: [{preview}]\n"));
            }
            if !ds.attributes.is_empty() {
                out.push_str("    attributes:\n");
                for (k, v) in &ds.attributes {
                    out.push_str(&format!("      {k} = {v}\n"));
                }
            }
        }
        out
    }
}
impl Default for ShdfFile {
    fn default() -> Self {
        Self::new()
    }
}
/// Helper for managing typed attributes.
pub struct AttributeHelper;
impl AttributeHelper {
    /// Serialize an attribute value to a string.
    pub fn to_string(val: &AttributeValue) -> String {
        match val {
            AttributeValue::String(s) => format!("s:{s}"),
            AttributeValue::Float64(f) => format!("f:{f}"),
            AttributeValue::Int32(i) => format!("i:{i}"),
            AttributeValue::Bool(b) => format!("b:{b}"),
        }
    }
    /// Deserialize an attribute value from a string.
    pub fn from_string(s: &str) -> AttributeValue {
        if let Some(rest) = s.strip_prefix("f:")
            && let Ok(f) = rest.parse::<f64>()
        {
            return AttributeValue::Float64(f);
        }
        if let Some(rest) = s.strip_prefix("i:")
            && let Ok(i) = rest.parse::<i32>()
        {
            return AttributeValue::Int32(i);
        }
        if let Some(rest) = s.strip_prefix("b:") {
            return AttributeValue::Bool(rest == "true");
        }
        if let Some(rest) = s.strip_prefix("s:") {
            return AttributeValue::String(rest.to_string());
        }
        AttributeValue::String(s.to_string())
    }
}
