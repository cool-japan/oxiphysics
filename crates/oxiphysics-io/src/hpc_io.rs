// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HPC and scientific data I/O.
//!
//! Provides mock implementations for HDF5, parallel netCDF, ADIOS2,
//! checkpoint management, distributed mesh I/O, performance logging,
//! and scientific JSON with base64-encoded float arrays.

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// Base64 encoding helper (no external crates)
// ---------------------------------------------------------------------------

const BASE64_CHARS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

/// Encode bytes to base64 string.
pub fn base64_encode(data: &[u8]) -> String {
    let mut result = String::new();
    let mut i = 0;
    while i + 2 < data.len() {
        let b0 = data[i] as usize;
        let b1 = data[i + 1] as usize;
        let b2 = data[i + 2] as usize;
        result.push(BASE64_CHARS[b0 >> 2] as char);
        result.push(BASE64_CHARS[((b0 & 3) << 4) | (b1 >> 4)] as char);
        result.push(BASE64_CHARS[((b1 & 0xf) << 2) | (b2 >> 6)] as char);
        result.push(BASE64_CHARS[b2 & 0x3f] as char);
        i += 3;
    }
    let remaining = data.len() - i;
    if remaining == 1 {
        let b0 = data[i] as usize;
        result.push(BASE64_CHARS[b0 >> 2] as char);
        result.push(BASE64_CHARS[(b0 & 3) << 4] as char);
        result.push('=');
        result.push('=');
    } else if remaining == 2 {
        let b0 = data[i] as usize;
        let b1 = data[i + 1] as usize;
        result.push(BASE64_CHARS[b0 >> 2] as char);
        result.push(BASE64_CHARS[((b0 & 3) << 4) | (b1 >> 4)] as char);
        result.push(BASE64_CHARS[(b1 & 0xf) << 2] as char);
        result.push('=');
    }
    result
}

/// Encode a slice of f64 values to base64.
pub fn f64_slice_to_base64(data: &[f64]) -> String {
    let bytes: Vec<u8> = data.iter().flat_map(|v| v.to_le_bytes().to_vec()).collect();
    base64_encode(&bytes)
}

// ---------------------------------------------------------------------------
// Hdf5Dataset
// ---------------------------------------------------------------------------

/// Supported HDF5-like data types.
#[derive(Clone, Debug, PartialEq)]
pub enum H5Dtype {
    /// 32-bit float.
    Float32,
    /// 64-bit float.
    Float64,
    /// 32-bit signed integer.
    Int32,
    /// 64-bit signed integer.
    Int64,
    /// Unsigned 8-bit integer.
    Uint8,
}

impl H5Dtype {
    /// Size in bytes.
    pub fn byte_size(&self) -> usize {
        match self {
            H5Dtype::Float32 => 4,
            H5Dtype::Float64 => 8,
            H5Dtype::Int32 => 4,
            H5Dtype::Int64 => 8,
            H5Dtype::Uint8 => 1,
        }
    }
}

/// Mock HDF5 dataset.
#[derive(Clone, Debug)]
pub struct Hdf5Dataset {
    /// Dataset name.
    pub name: String,
    /// Dimensions.
    pub dims: Vec<usize>,
    /// Data type.
    pub dtype: H5Dtype,
    /// Chunk size (for chunked storage).
    pub chunk_size: Option<Vec<usize>>,
    /// Compression level (0 = none, 9 = max).
    pub compression_level: u8,
    /// Dataset attributes.
    pub attributes: HashMap<String, String>,
    /// Raw data storage (f64 for simplicity).
    pub data: Vec<f64>,
}

impl Hdf5Dataset {
    /// Create a new dataset.
    pub fn new(name: &str, dims: Vec<usize>, dtype: H5Dtype) -> Self {
        let n: usize = dims.iter().product();
        Self {
            name: name.to_string(),
            dims,
            dtype,
            chunk_size: None,
            compression_level: 0,
            attributes: HashMap::new(),
            data: vec![0.0; n],
        }
    }

    /// Set chunk size.
    pub fn set_chunk_size(&mut self, chunk: Vec<usize>) {
        self.chunk_size = Some(chunk);
    }

    /// Set compression level.
    pub fn set_compression(&mut self, level: u8) {
        self.compression_level = level.min(9);
    }

    /// Set an attribute.
    pub fn set_attr(&mut self, key: &str, value: &str) {
        self.attributes.insert(key.to_string(), value.to_string());
    }

    /// Get an attribute.
    pub fn get_attr(&self, key: &str) -> Option<&str> {
        self.attributes.get(key).map(|s| s.as_str())
    }

    /// Total number of elements.
    pub fn n_elements(&self) -> usize {
        self.dims.iter().product()
    }

    /// Estimated memory footprint (bytes).
    pub fn memory_bytes(&self) -> usize {
        self.n_elements() * self.dtype.byte_size()
    }

    /// Write data slice.
    pub fn write_slice(&mut self, offset: usize, values: &[f64]) {
        let end = (offset + values.len()).min(self.data.len());
        for (i, &v) in values.iter().enumerate() {
            if offset + i < end {
                self.data[offset + i] = v;
            }
        }
    }

    /// Read data slice.
    pub fn read_slice(&self, offset: usize, length: usize) -> Vec<f64> {
        let end = (offset + length).min(self.data.len());
        self.data[offset..end].to_vec()
    }
}

// ---------------------------------------------------------------------------
// Hdf5Group
// ---------------------------------------------------------------------------

/// Mock HDF5 group.
#[derive(Clone, Debug)]
pub struct Hdf5Group {
    /// Group name.
    pub name: String,
    /// Child group names.
    pub children: Vec<String>,
    /// Group attributes.
    pub attributes: HashMap<String, String>,
    /// Datasets in this group.
    pub datasets: HashMap<String, Hdf5Dataset>,
}

impl Hdf5Group {
    /// Create a new group.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            children: Vec::new(),
            attributes: HashMap::new(),
            datasets: HashMap::new(),
        }
    }

    /// Add a child group name.
    pub fn add_child(&mut self, child_name: &str) {
        self.children.push(child_name.to_string());
    }

    /// Set a group attribute.
    pub fn set_attr(&mut self, key: &str, value: &str) {
        self.attributes.insert(key.to_string(), value.to_string());
    }

    /// Create a dataset in this group.
    pub fn create_dataset(&mut self, name: &str, dims: Vec<usize>, dtype: H5Dtype) {
        let ds = Hdf5Dataset::new(name, dims, dtype);
        self.datasets.insert(name.to_string(), ds);
    }

    /// Get a dataset.
    pub fn get_dataset(&self, name: &str) -> Option<&Hdf5Dataset> {
        self.datasets.get(name)
    }

    /// Get mutable dataset.
    pub fn get_dataset_mut(&mut self, name: &str) -> Option<&mut Hdf5Dataset> {
        self.datasets.get_mut(name)
    }

    /// List dataset names.
    pub fn dataset_names(&self) -> Vec<&str> {
        self.datasets.keys().map(|s| s.as_str()).collect()
    }
}

// ---------------------------------------------------------------------------
// Hdf5File
// ---------------------------------------------------------------------------

/// Mock HDF5 file.
#[derive(Clone, Debug)]
pub struct Hdf5File {
    /// File path.
    pub path: String,
    /// Root group.
    pub root: Hdf5Group,
    /// Whether file is open.
    pub is_open: bool,
    /// Read-only flag.
    pub read_only: bool,
}

impl Hdf5File {
    /// Open an existing file.
    pub fn open(path: &str) -> Self {
        Self {
            path: path.to_string(),
            root: Hdf5Group::new("/"),
            is_open: true,
            read_only: true,
        }
    }

    /// Create a new file.
    pub fn create(path: &str) -> Self {
        Self {
            path: path.to_string(),
            root: Hdf5Group::new("/"),
            is_open: true,
            read_only: false,
        }
    }

    /// Close the file.
    pub fn close(&mut self) {
        self.is_open = false;
    }

    /// Flush (no-op in mock).
    pub fn flush(&self) {}

    /// Create a group at the root.
    pub fn create_group(&mut self, name: &str) -> &mut Hdf5Group {
        self.root.add_child(name);
        // Insert if not present
        if !self.root.children.contains(&name.to_string()) {
            self.root.children.push(name.to_string());
        }
        &mut self.root
    }

    /// Get root group.
    pub fn root_group(&self) -> &Hdf5Group {
        &self.root
    }

    /// Get mutable root group.
    pub fn root_group_mut(&mut self) -> &mut Hdf5Group {
        &mut self.root
    }

    /// Create a dataset directly at root.
    pub fn create_dataset(&mut self, name: &str, dims: Vec<usize>, dtype: H5Dtype) {
        self.root.create_dataset(name, dims, dtype);
    }
}

// ---------------------------------------------------------------------------
// ParallelNetcdf
// ---------------------------------------------------------------------------

/// Parallel netCDF stub for distributed array I/O.
#[derive(Clone, Debug)]
pub struct ParallelNetcdf {
    /// File path.
    pub path: String,
    /// Global dimensions \[nx, ny, nz\].
    pub global_dims: Vec<usize>,
    /// Number of MPI ranks (simulated).
    pub n_ranks: usize,
    /// Current rank.
    pub rank: usize,
    /// Variables stored (name -> local data).
    pub variables: HashMap<String, Vec<f64>>,
    /// Variable attributes.
    pub var_attrs: HashMap<String, HashMap<String, String>>,
}

impl ParallelNetcdf {
    /// Create a new parallel netCDF file.
    pub fn new(path: &str, global_dims: Vec<usize>, n_ranks: usize, rank: usize) -> Self {
        Self {
            path: path.to_string(),
            global_dims,
            n_ranks,
            rank,
            variables: HashMap::new(),
            var_attrs: HashMap::new(),
        }
    }

    /// Define a variable.
    pub fn def_var(&mut self, name: &str, _dims: &[&str]) {
        self.variables.insert(name.to_string(), Vec::new());
        self.var_attrs.insert(name.to_string(), HashMap::new());
    }

    /// Write distributed array (local portion).
    pub fn put_var(&mut self, name: &str, data: Vec<f64>) {
        self.variables.insert(name.to_string(), data);
    }

    /// Read variable.
    pub fn get_var(&self, name: &str) -> Option<&[f64]> {
        self.variables.get(name).map(|v| v.as_slice())
    }

    /// Set variable attribute.
    pub fn set_var_attr(&mut self, var: &str, key: &str, value: &str) {
        if let Some(attrs) = self.var_attrs.get_mut(var) {
            attrs.insert(key.to_string(), value.to_string());
        }
    }

    /// Total global size.
    pub fn global_size(&self) -> usize {
        self.global_dims.iter().product()
    }

    /// Local size for this rank (simple block decomposition along first dim).
    pub fn local_size(&self) -> usize {
        if self.global_dims.is_empty() {
            return 0;
        }
        let n_global_first = self.global_dims[0];
        let local_first = n_global_first.div_ceil(self.n_ranks);
        let rest: usize = self.global_dims[1..].iter().product::<usize>().max(1);
        local_first.min(n_global_first) * rest
    }
}

// ---------------------------------------------------------------------------
// AdiosWriter
// ---------------------------------------------------------------------------

/// ADIOS2-style streaming I/O writer.
#[derive(Clone, Debug)]
pub struct AdiosWriter {
    /// Stream name.
    pub stream_name: String,
    /// Open flag.
    pub is_open: bool,
    /// Variables staged for write (name -> data).
    pub staged_vars: HashMap<String, Vec<f64>>,
    /// Variable metadata (shape, type).
    pub var_meta: HashMap<String, (Vec<usize>, String)>,
    /// Steps written.
    pub steps_written: usize,
}

impl AdiosWriter {
    /// Open an ADIOS2 output stream.
    pub fn open(stream_name: &str) -> Self {
        Self {
            stream_name: stream_name.to_string(),
            is_open: true,
            staged_vars: HashMap::new(),
            var_meta: HashMap::new(),
            steps_written: 0,
        }
    }

    /// Define a variable with shape.
    pub fn define_variable(&mut self, name: &str, shape: Vec<usize>, dtype: &str) {
        self.var_meta
            .insert(name.to_string(), (shape, dtype.to_string()));
    }

    /// Put (stage) a variable for writing.
    pub fn put_variable(&mut self, name: &str, data: Vec<f64>) {
        self.staged_vars.insert(name.to_string(), data);
    }

    /// Perform puts: flush staged variables as a step.
    pub fn perform_puts(&mut self) {
        // In real ADIOS2, this triggers I/O; here we just increment step
        self.steps_written += 1;
        self.staged_vars.clear();
    }

    /// Close the stream.
    pub fn close(&mut self) {
        self.is_open = false;
    }

    /// Check if variable is defined.
    pub fn has_variable(&self, name: &str) -> bool {
        self.var_meta.contains_key(name)
    }

    /// Get variable shape.
    pub fn variable_shape(&self, name: &str) -> Option<&[usize]> {
        self.var_meta.get(name).map(|(s, _)| s.as_slice())
    }
}

// ---------------------------------------------------------------------------
// CheckpointManager
// ---------------------------------------------------------------------------

/// A single checkpoint entry.
#[derive(Clone, Debug)]
pub struct CheckpointEntry {
    /// Checkpoint index.
    pub index: usize,
    /// Simulation step.
    pub step: u64,
    /// File path.
    pub path: String,
    /// Timestamp (simulated, seconds from epoch).
    pub timestamp: f64,
}

/// Versioned checkpoint manager.
#[derive(Clone, Debug)]
pub struct CheckpointManager {
    /// Base directory for checkpoints.
    pub base_dir: String,
    /// Maximum number of checkpoints to keep.
    pub keep_last_n: usize,
    /// List of checkpoint entries (oldest first).
    pub entries: Vec<CheckpointEntry>,
    /// Next checkpoint index.
    pub next_index: usize,
}

impl CheckpointManager {
    /// Create a new checkpoint manager.
    pub fn new(base_dir: &str, keep_last_n: usize) -> Self {
        Self {
            base_dir: base_dir.to_string(),
            keep_last_n,
            entries: Vec::new(),
            next_index: 0,
        }
    }

    /// Register a new checkpoint.
    pub fn register(&mut self, step: u64, timestamp: f64) -> String {
        let path = format!("{}/checkpoint_{:06}.bin", self.base_dir, self.next_index);
        let entry = CheckpointEntry {
            index: self.next_index,
            step,
            path: path.clone(),
            timestamp,
        };
        self.entries.push(entry);
        self.next_index += 1;

        // Remove old entries beyond keep_last_n
        while self.entries.len() > self.keep_last_n {
            self.entries.remove(0);
        }

        path
    }

    /// Get the latest checkpoint.
    pub fn latest(&self) -> Option<&CheckpointEntry> {
        self.entries.last()
    }

    /// Restore by index (returns path if found).
    pub fn restore_by_index(&self, index: usize) -> Option<&CheckpointEntry> {
        self.entries.iter().find(|e| e.index == index)
    }

    /// Number of stored checkpoints.
    pub fn n_checkpoints(&self) -> usize {
        self.entries.len()
    }
}

// ---------------------------------------------------------------------------
// RestartFile
// ---------------------------------------------------------------------------

/// Binary restart file for simulation state.
#[derive(Clone, Debug)]
pub struct RestartFile {
    /// File path.
    pub path: String,
    /// State vector.
    pub state: Vec<f64>,
    /// Metadata: simulation step.
    pub step: u64,
    /// Metadata: simulation time.
    pub time: f64,
    /// Metadata: arbitrary key-value pairs.
    pub metadata: HashMap<String, f64>,
}

impl RestartFile {
    /// Create a new restart file.
    pub fn new(path: &str) -> Self {
        Self {
            path: path.to_string(),
            state: Vec::new(),
            step: 0,
            time: 0.0,
            metadata: HashMap::new(),
        }
    }

    /// Write state to a byte buffer (mock serialization).
    pub fn write_to_buffer(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        // Magic number
        buf.extend_from_slice(b"OXIRS001");
        // Step (u64)
        buf.extend_from_slice(&self.step.to_le_bytes());
        // Time (f64)
        buf.extend_from_slice(&self.time.to_le_bytes());
        // State length (u64)
        let n = self.state.len() as u64;
        buf.extend_from_slice(&n.to_le_bytes());
        // State data
        for &v in &self.state {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    /// Read state from byte buffer (mock deserialization).
    pub fn read_from_buffer(path: &str, buf: &[u8]) -> Option<Self> {
        if buf.len() < 24 || &buf[..8] != b"OXIRS001" {
            return None;
        }
        let step = u64::from_le_bytes(buf[8..16].try_into().ok()?);
        let time = f64::from_le_bytes(buf[16..24].try_into().ok()?);
        let n = u64::from_le_bytes(buf[24..32].try_into().ok()?) as usize;
        let mut state = Vec::with_capacity(n);
        for i in 0..n {
            let off = 32 + i * 8;
            if off + 8 > buf.len() {
                break;
            }
            let v = f64::from_le_bytes(buf[off..off + 8].try_into().ok()?);
            state.push(v);
        }
        Some(Self {
            path: path.to_string(),
            state,
            step,
            time,
            metadata: HashMap::new(),
        })
    }

    /// Set a metadata value.
    pub fn set_meta(&mut self, key: &str, value: f64) {
        self.metadata.insert(key.to_string(), value);
    }
}

// ---------------------------------------------------------------------------
// DistributedMeshIO
// ---------------------------------------------------------------------------

/// I/O for a distributed mesh partition.
#[derive(Clone, Debug)]
pub struct DistributedMeshIO {
    /// Rank index.
    pub rank: usize,
    /// Total number of ranks.
    pub n_ranks: usize,
    /// Local node positions.
    pub local_nodes: Vec<[f64; 3]>,
    /// Ghost node positions (copies from neighbor ranks).
    pub ghost_nodes: Vec<[f64; 3]>,
    /// Global IDs of local nodes.
    pub local_global_ids: Vec<u64>,
    /// Global IDs of ghost nodes.
    pub ghost_global_ids: Vec<u64>,
    /// Connectivity (local element node indices).
    pub elements: Vec<Vec<usize>>,
}

impl DistributedMeshIO {
    /// Create a new distributed mesh partition.
    pub fn new(rank: usize, n_ranks: usize) -> Self {
        Self {
            rank,
            n_ranks,
            local_nodes: Vec::new(),
            ghost_nodes: Vec::new(),
            local_global_ids: Vec::new(),
            ghost_global_ids: Vec::new(),
            elements: Vec::new(),
        }
    }

    /// Add a local node.
    pub fn add_local_node(&mut self, pos: [f64; 3], global_id: u64) {
        self.local_nodes.push(pos);
        self.local_global_ids.push(global_id);
    }

    /// Add a ghost node.
    pub fn add_ghost_node(&mut self, pos: [f64; 3], global_id: u64) {
        self.ghost_nodes.push(pos);
        self.ghost_global_ids.push(global_id);
    }

    /// Add an element (connectivity).
    pub fn add_element(&mut self, nodes: Vec<usize>) {
        self.elements.push(nodes);
    }

    /// Total local + ghost nodes.
    pub fn total_nodes(&self) -> usize {
        self.local_nodes.len() + self.ghost_nodes.len()
    }

    /// Serialize local mesh to bytes (mock).
    pub fn serialize(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(b"OXIDMESH");
        let n_local = self.local_nodes.len() as u64;
        buf.extend_from_slice(&n_local.to_le_bytes());
        for node in &self.local_nodes {
            for &c in node {
                buf.extend_from_slice(&c.to_le_bytes());
            }
        }
        buf
    }

    /// Build a simple block-decomposed partition from a regular grid.
    pub fn partition_regular_grid(
        n_nodes_total: usize,
        positions: &[[f64; 3]],
        rank: usize,
        n_ranks: usize,
    ) -> Self {
        let block = n_nodes_total.div_ceil(n_ranks);
        let start = rank * block;
        let end = (start + block).min(n_nodes_total);
        let mut mesh = Self::new(rank, n_ranks);
        for (i, &pos) in positions[start..end].iter().enumerate() {
            mesh.add_local_node(pos, (start + i) as u64);
        }
        mesh
    }
}

// ---------------------------------------------------------------------------
// PerformanceLog
// ---------------------------------------------------------------------------

/// A single performance log entry.
#[derive(Clone, Debug)]
pub struct PerfEntry {
    /// Crate name.
    pub crate_name: String,
    /// Function name.
    pub function: String,
    /// Elapsed time in milliseconds.
    pub time_ms: f64,
    /// Memory usage in megabytes.
    pub memory_mb: f64,
    /// Number of threads.
    pub n_threads: usize,
    /// Additional notes.
    pub notes: String,
}

impl PerfEntry {
    /// Create a new performance entry.
    pub fn new(
        crate_name: &str,
        function: &str,
        time_ms: f64,
        memory_mb: f64,
        n_threads: usize,
    ) -> Self {
        Self {
            crate_name: crate_name.to_string(),
            function: function.to_string(),
            time_ms,
            memory_mb,
            n_threads,
            notes: String::new(),
        }
    }
}

/// Structured performance log.
#[derive(Clone, Debug)]
pub struct PerformanceLog {
    /// All recorded entries.
    pub entries: Vec<PerfEntry>,
    /// Run identifier.
    pub run_id: String,
}

impl PerformanceLog {
    /// Create a new performance log.
    pub fn new(run_id: &str) -> Self {
        Self {
            entries: Vec::new(),
            run_id: run_id.to_string(),
        }
    }

    /// Record an entry.
    pub fn record(&mut self, entry: PerfEntry) {
        self.entries.push(entry);
    }

    /// Total elapsed time across all entries.
    pub fn total_time_ms(&self) -> f64 {
        self.entries.iter().map(|e| e.time_ms).sum()
    }

    /// Mean time per entry.
    pub fn mean_time_ms(&self) -> f64 {
        if self.entries.is_empty() {
            return 0.0;
        }
        self.total_time_ms() / self.entries.len() as f64
    }

    /// Maximum memory usage.
    pub fn peak_memory_mb(&self) -> f64 {
        self.entries
            .iter()
            .map(|e| e.memory_mb)
            .fold(0.0f64, f64::max)
    }

    /// Filter entries by crate name.
    pub fn filter_by_crate(&self, crate_name: &str) -> Vec<&PerfEntry> {
        self.entries
            .iter()
            .filter(|e| e.crate_name == crate_name)
            .collect()
    }

    /// Serialize to CSV-like string.
    pub fn to_csv(&self) -> String {
        let mut s = "crate,function,time_ms,memory_mb,n_threads\n".to_string();
        for e in &self.entries {
            s.push_str(&format!(
                "{},{},{:.3},{:.3},{}\n",
                e.crate_name, e.function, e.time_ms, e.memory_mb, e.n_threads
            ));
        }
        s
    }

    /// Serialize to JSON string.
    pub fn to_json(&self) -> String {
        let entries_json: Vec<String> = self
            .entries
            .iter()
            .map(|e| {
                format!(
                    r#"{{"crate":"{}","function":"{}","time_ms":{:.3},"memory_mb":{:.3},"n_threads":{}}}"#,
                    e.crate_name, e.function, e.time_ms, e.memory_mb, e.n_threads
                )
            })
            .collect();
        format!(
            r#"{{"run_id":"{}","entries":[{}]}}"#,
            self.run_id,
            entries_json.join(",")
        )
    }
}

// ---------------------------------------------------------------------------
// ScientificJson
// ---------------------------------------------------------------------------

/// Large array JSON using base64-encoded float arrays to avoid huge string arrays.
#[derive(Clone, Debug)]
pub struct ScientificJson {
    /// Document name.
    pub name: String,
    /// Metadata fields (string).
    pub metadata: HashMap<String, String>,
    /// Named float arrays (stored as base64 in JSON).
    pub arrays: HashMap<String, Vec<f64>>,
    /// Named scalar values.
    pub scalars: HashMap<String, f64>,
}

impl ScientificJson {
    /// Create a new scientific JSON document.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            metadata: HashMap::new(),
            arrays: HashMap::new(),
            scalars: HashMap::new(),
        }
    }

    /// Set metadata.
    pub fn set_meta(&mut self, key: &str, value: &str) {
        self.metadata.insert(key.to_string(), value.to_string());
    }

    /// Set scalar.
    pub fn set_scalar(&mut self, key: &str, value: f64) {
        self.scalars.insert(key.to_string(), value);
    }

    /// Add a named float array.
    pub fn add_array(&mut self, name: &str, data: Vec<f64>) {
        self.arrays.insert(name.to_string(), data);
    }

    /// Get a named float array.
    pub fn get_array(&self, name: &str) -> Option<&[f64]> {
        self.arrays.get(name).map(|v| v.as_slice())
    }

    /// Serialize to JSON with base64-encoded arrays.
    pub fn to_json(&self) -> String {
        let mut parts: Vec<String> = Vec::new();

        // Metadata
        let meta_parts: Vec<String> = self
            .metadata
            .iter()
            .map(|(k, v)| format!(r#""{}":"{}""#, k, v))
            .collect();
        if !meta_parts.is_empty() {
            parts.push(format!(r#""metadata":{{{}}}"#, meta_parts.join(",")));
        }

        // Scalars
        let scalar_parts: Vec<String> = self
            .scalars
            .iter()
            .map(|(k, v)| format!(r#""{}":{}"#, k, v))
            .collect();
        if !scalar_parts.is_empty() {
            parts.push(format!(r#""scalars":{{{}}}"#, scalar_parts.join(",")));
        }

        // Arrays: base64-encoded
        let arr_parts: Vec<String> = self
            .arrays
            .iter()
            .map(|(k, v)| {
                let b64 = f64_slice_to_base64(v);
                format!(
                    r#""{}":{{"dtype":"f64","shape":[{}],"base64":"{}"}}"#,
                    k,
                    v.len(),
                    b64
                )
            })
            .collect();
        if !arr_parts.is_empty() {
            parts.push(format!(r#""arrays":{{{}}}"#, arr_parts.join(",")));
        }

        format!(r#"{{"name":"{}",{}}}"#, self.name, parts.join(","))
    }

    /// Parse a simple key-value from JSON string (basic mock parser).
    pub fn parse_scalar(json: &str, key: &str) -> Option<f64> {
        let search = format!("\"{}\":", key);
        let start = json.find(&search)? + search.len();
        let rest = &json[start..];
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        rest[..end].trim().parse().ok()
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base64_encode_empty() {
        assert_eq!(base64_encode(&[]), "");
    }

    #[test]
    fn test_base64_encode_hello() {
        let encoded = base64_encode(b"Hello");
        assert_eq!(encoded, "SGVsbG8=");
    }

    #[test]
    fn test_f64_slice_to_base64_roundtrip_length() {
        let data = vec![1.0, 2.0, 3.0];
        let b64 = f64_slice_to_base64(&data);
        // 3 f64 * 8 bytes = 24 bytes -> base64 len = ceil(24/3)*4 = 32
        assert_eq!(b64.len(), 32);
    }

    #[test]
    fn test_h5dtype_byte_size() {
        assert_eq!(H5Dtype::Float32.byte_size(), 4);
        assert_eq!(H5Dtype::Float64.byte_size(), 8);
        assert_eq!(H5Dtype::Uint8.byte_size(), 1);
    }

    #[test]
    fn test_hdf5_dataset_new() {
        let ds = Hdf5Dataset::new("data", vec![10, 10], H5Dtype::Float64);
        assert_eq!(ds.n_elements(), 100);
        assert_eq!(ds.data.len(), 100);
    }

    #[test]
    fn test_hdf5_dataset_memory_bytes() {
        let ds = Hdf5Dataset::new("data", vec![10], H5Dtype::Float64);
        assert_eq!(ds.memory_bytes(), 80);
    }

    #[test]
    fn test_hdf5_dataset_write_read_slice() {
        let mut ds = Hdf5Dataset::new("data", vec![10], H5Dtype::Float64);
        ds.write_slice(2, &[1.0, 2.0, 3.0]);
        let r = ds.read_slice(2, 3);
        assert_eq!(r, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_hdf5_dataset_set_attr() {
        let mut ds = Hdf5Dataset::new("data", vec![5], H5Dtype::Float32);
        ds.set_attr("units", "m/s");
        assert_eq!(ds.get_attr("units"), Some("m/s"));
    }

    #[test]
    fn test_hdf5_dataset_compression() {
        let mut ds = Hdf5Dataset::new("data", vec![5], H5Dtype::Float32);
        ds.set_compression(15);
        assert_eq!(ds.compression_level, 9); // capped
    }

    #[test]
    fn test_hdf5_group_new() {
        let g = Hdf5Group::new("sim");
        assert_eq!(g.name, "sim");
        assert!(g.datasets.is_empty());
    }

    #[test]
    fn test_hdf5_group_create_dataset() {
        let mut g = Hdf5Group::new("results");
        g.create_dataset("velocity", vec![100, 3], H5Dtype::Float64);
        assert!(g.get_dataset("velocity").is_some());
    }

    #[test]
    fn test_hdf5_group_dataset_names() {
        let mut g = Hdf5Group::new("results");
        g.create_dataset("pressure", vec![100], H5Dtype::Float64);
        g.create_dataset("density", vec![100], H5Dtype::Float64);
        let names = g.dataset_names();
        assert_eq!(names.len(), 2);
    }

    #[test]
    fn test_hdf5_file_create() {
        let f = Hdf5File::create("test.h5");
        assert!(f.is_open);
        assert!(!f.read_only);
    }

    #[test]
    fn test_hdf5_file_open() {
        let f = Hdf5File::open("test.h5");
        assert!(f.is_open);
        assert!(f.read_only);
    }

    #[test]
    fn test_hdf5_file_create_dataset() {
        let mut f = Hdf5File::create("test.h5");
        f.create_dataset("temperatures", vec![50], H5Dtype::Float32);
        assert!(f.root.get_dataset("temperatures").is_some());
    }

    #[test]
    fn test_hdf5_file_close() {
        let mut f = Hdf5File::create("test.h5");
        f.close();
        assert!(!f.is_open);
    }

    #[test]
    fn test_parallel_netcdf_new() {
        let nc = ParallelNetcdf::new("out.nc", vec![100, 100, 100], 4, 0);
        assert_eq!(nc.global_size(), 1_000_000);
    }

    #[test]
    fn test_parallel_netcdf_put_get_var() {
        let mut nc = ParallelNetcdf::new("out.nc", vec![10], 1, 0);
        nc.def_var("temp", &["x"]);
        nc.put_var("temp", vec![1.0, 2.0, 3.0]);
        let v = nc.get_var("temp").unwrap();
        assert_eq!(v, &[1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_parallel_netcdf_local_size() {
        let nc = ParallelNetcdf::new("out.nc", vec![100, 10], 4, 0);
        let ls = nc.local_size();
        assert!(ls > 0);
    }

    #[test]
    fn test_adios_writer_open() {
        let w = AdiosWriter::open("sim.bp");
        assert!(w.is_open);
        assert_eq!(w.steps_written, 0);
    }

    #[test]
    fn test_adios_writer_define_and_put() {
        let mut w = AdiosWriter::open("sim.bp");
        w.define_variable("velocity", vec![100, 3], "f64");
        assert!(w.has_variable("velocity"));
        w.put_variable("velocity", vec![1.0; 300]);
        assert!(w.staged_vars.contains_key("velocity"));
    }

    #[test]
    fn test_adios_writer_perform_puts() {
        let mut w = AdiosWriter::open("sim.bp");
        w.put_variable("p", vec![1.0]);
        w.perform_puts();
        assert_eq!(w.steps_written, 1);
        assert!(w.staged_vars.is_empty());
    }

    #[test]
    fn test_adios_writer_close() {
        let mut w = AdiosWriter::open("sim.bp");
        w.close();
        assert!(!w.is_open);
    }

    #[test]
    fn test_checkpoint_manager_register() {
        let base = std::env::temp_dir().join("hpc_ckpt_register");
        let mut cm = CheckpointManager::new(base.to_str().unwrap_or(""), 3);
        let path = cm.register(100, 1.0);
        assert!(path.contains("checkpoint_000000"));
        assert_eq!(cm.n_checkpoints(), 1);
    }

    #[test]
    fn test_checkpoint_manager_keep_last_n() {
        let base = std::env::temp_dir().join("hpc_ckpt_keepn");
        let mut cm = CheckpointManager::new(base.to_str().unwrap_or(""), 2);
        cm.register(0, 0.0);
        cm.register(1, 1.0);
        cm.register(2, 2.0);
        assert_eq!(cm.n_checkpoints(), 2);
    }

    #[test]
    fn test_checkpoint_manager_latest() {
        let base = std::env::temp_dir().join("hpc_ckpt_latest");
        let mut cm = CheckpointManager::new(base.to_str().unwrap_or(""), 5);
        cm.register(10, 1.0);
        cm.register(20, 2.0);
        assert_eq!(cm.latest().unwrap().step, 20);
    }

    #[test]
    fn test_checkpoint_manager_restore_by_index() {
        let base = std::env::temp_dir().join("hpc_ckpt_restore");
        let mut cm = CheckpointManager::new(base.to_str().unwrap_or(""), 5);
        cm.register(10, 1.0);
        let e = cm.restore_by_index(0);
        assert!(e.is_some());
        assert_eq!(e.unwrap().step, 10);
    }

    #[test]
    fn test_restart_file_write_read() {
        let mut rf = RestartFile::new("restart.bin");
        rf.state = vec![1.0, 2.0, 3.0];
        rf.step = 500;
        rf.time = 0.5;
        let buf = rf.write_to_buffer();
        let rf2 = RestartFile::read_from_buffer("restart.bin", &buf).unwrap();
        assert_eq!(rf2.step, 500);
        assert!((rf2.time - 0.5).abs() < 1e-10);
        assert_eq!(rf2.state, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_restart_file_invalid_magic() {
        let buf = vec![0u8; 64];
        let rf = RestartFile::read_from_buffer("x.bin", &buf);
        assert!(rf.is_none());
    }

    #[test]
    fn test_distributed_mesh_io_new() {
        let m = DistributedMeshIO::new(0, 4);
        assert_eq!(m.rank, 0);
        assert_eq!(m.total_nodes(), 0);
    }

    #[test]
    fn test_distributed_mesh_io_add_nodes() {
        let mut m = DistributedMeshIO::new(0, 4);
        m.add_local_node([0.0, 0.0, 0.0], 0);
        m.add_ghost_node([1.0, 0.0, 0.0], 100);
        assert_eq!(m.total_nodes(), 2);
        assert_eq!(m.local_nodes.len(), 1);
    }

    #[test]
    fn test_distributed_mesh_io_serialize() {
        let mut m = DistributedMeshIO::new(0, 1);
        m.add_local_node([1.0, 2.0, 3.0], 0);
        let buf = m.serialize();
        assert!(buf.starts_with(b"OXIDMESH"));
    }

    #[test]
    fn test_distributed_mesh_io_partition() {
        let positions: Vec<[f64; 3]> = (0..10).map(|i| [i as f64, 0.0, 0.0]).collect();
        let m = DistributedMeshIO::partition_regular_grid(10, &positions, 0, 2);
        assert!(m.local_nodes.len() >= 4);
    }

    #[test]
    fn test_performance_log_record() {
        let mut log = PerformanceLog::new("run001");
        log.record(PerfEntry::new("lbm", "stream", 10.0, 100.0, 4));
        assert_eq!(log.entries.len(), 1);
    }

    #[test]
    fn test_performance_log_total_time() {
        let mut log = PerformanceLog::new("run001");
        log.record(PerfEntry::new("lbm", "collide", 5.0, 50.0, 4));
        log.record(PerfEntry::new("lbm", "stream", 3.0, 50.0, 4));
        assert!((log.total_time_ms() - 8.0).abs() < 1e-10);
    }

    #[test]
    fn test_performance_log_peak_memory() {
        let mut log = PerformanceLog::new("run001");
        log.record(PerfEntry::new("md", "force", 20.0, 500.0, 8));
        log.record(PerfEntry::new("md", "integrate", 5.0, 100.0, 8));
        assert!((log.peak_memory_mb() - 500.0).abs() < 1e-10);
    }

    #[test]
    fn test_performance_log_filter_by_crate() {
        let mut log = PerformanceLog::new("run001");
        log.record(PerfEntry::new("lbm", "stream", 5.0, 50.0, 4));
        log.record(PerfEntry::new("md", "force", 10.0, 100.0, 8));
        let lbm_entries = log.filter_by_crate("lbm");
        assert_eq!(lbm_entries.len(), 1);
    }

    #[test]
    fn test_performance_log_to_csv() {
        let mut log = PerformanceLog::new("run001");
        log.record(PerfEntry::new("lbm", "stream", 5.0, 50.0, 4));
        let csv = log.to_csv();
        assert!(csv.contains("lbm"));
        assert!(csv.contains("stream"));
    }

    #[test]
    fn test_scientific_json_new() {
        let sj = ScientificJson::new("simulation");
        assert_eq!(sj.name, "simulation");
    }

    #[test]
    fn test_scientific_json_add_array() {
        let mut sj = ScientificJson::new("sim");
        sj.add_array("velocity", vec![1.0, 2.0, 3.0]);
        assert_eq!(sj.get_array("velocity"), Some([1.0, 2.0, 3.0].as_slice()));
    }

    #[test]
    fn test_scientific_json_to_json_contains_base64() {
        let mut sj = ScientificJson::new("sim");
        sj.add_array("pressure", vec![1.0, 2.0]);
        let json = sj.to_json();
        assert!(json.contains("base64"));
        assert!(json.contains("pressure"));
    }

    #[test]
    fn test_scientific_json_scalar() {
        let mut sj = ScientificJson::new("sim");
        sj.set_scalar("dt", 0.001);
        let json = sj.to_json();
        let val = ScientificJson::parse_scalar(&json, "dt");
        assert!(val.is_some());
        assert!((val.unwrap() - 0.001).abs() < 1e-10);
    }
}
