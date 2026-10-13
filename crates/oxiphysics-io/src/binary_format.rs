// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HDF5-inspired binary format for simulation checkpoints and time series data.
//!
//! The `OxiFile` format stores hierarchical groups of typed datasets with
//! attributes, serialized to a compact binary representation.

use std::fs;
use std::io::Write;

// ---------------------------------------------------------------------------
// DataType
// ---------------------------------------------------------------------------

/// Supported element types for datasets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DataType {
    Float32,
    Float64,
    Int32,
    Int64,
    UInt8,
    UInt32,
}

impl DataType {
    /// Returns the byte size of a single element of this type.
    pub fn size_bytes(&self) -> usize {
        match self {
            DataType::Float32 => 4,
            DataType::Float64 => 8,
            DataType::Int32 => 4,
            DataType::Int64 => 8,
            DataType::UInt8 => 1,
            DataType::UInt32 => 4,
        }
    }

    fn tag(&self) -> u8 {
        match self {
            DataType::Float32 => 0,
            DataType::Float64 => 1,
            DataType::Int32 => 2,
            DataType::Int64 => 3,
            DataType::UInt8 => 4,
            DataType::UInt32 => 5,
        }
    }

    fn from_tag(tag: u8) -> Result<Self, String> {
        match tag {
            0 => Ok(DataType::Float32),
            1 => Ok(DataType::Float64),
            2 => Ok(DataType::Int32),
            3 => Ok(DataType::Int64),
            4 => Ok(DataType::UInt8),
            5 => Ok(DataType::UInt32),
            _ => Err(format!("Unknown DataType tag: {}", tag)),
        }
    }
}

// ---------------------------------------------------------------------------
// DatasetShape
// ---------------------------------------------------------------------------

/// Describes the shape (dimensions) of a dataset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DatasetShape {
    pub dims: Vec<usize>,
}

impl DatasetShape {
    /// Creates a new shape from a list of dimensions.
    pub fn new(dims: Vec<usize>) -> Self {
        Self { dims }
    }

    /// Total number of elements (product of all dimensions).
    pub fn total_elements(&self) -> usize {
        if self.dims.is_empty() {
            return 1; // scalar
        }
        self.dims.iter().product()
    }

    /// Total number of elements, or `None` if the product of the dimensions
    /// overflows `usize`.
    ///
    /// Use this instead of [`DatasetShape::total_elements`] for shapes that
    /// come from untrusted input.
    pub fn checked_total_elements(&self) -> Option<usize> {
        self.dims
            .iter()
            .try_fold(1usize, |acc, &dim| acc.checked_mul(dim))
    }

    /// Returns `true` if the shape represents a scalar (no dimensions).
    pub fn is_scalar(&self) -> bool {
        self.dims.is_empty()
    }

    /// Number of dimensions.
    pub fn rank(&self) -> usize {
        self.dims.len()
    }
}

// ---------------------------------------------------------------------------
// Attribute / AttributeValue
// ---------------------------------------------------------------------------

/// A named metadata attribute attached to a dataset or group.
#[derive(Debug, Clone, PartialEq)]
pub struct Attribute {
    pub name: String,
    pub value: AttributeValue,
}

impl Attribute {
    /// Creates a new attribute with the given name and value.
    pub fn new(name: impl Into<String>, value: AttributeValue) -> Self {
        Self {
            name: name.into(),
            value,
        }
    }
}

/// Possible values for an attribute.
#[derive(Debug, Clone, PartialEq)]
pub enum AttributeValue {
    Int(i64),
    Float(f64),
    Text(String),
    IntArray(Vec<i64>),
    FloatArray(Vec<f64>),
}

impl AttributeValue {
    fn tag(&self) -> u8 {
        match self {
            AttributeValue::Int(_) => 0,
            AttributeValue::Float(_) => 1,
            AttributeValue::Text(_) => 2,
            AttributeValue::IntArray(_) => 3,
            AttributeValue::FloatArray(_) => 4,
        }
    }
}

// ---------------------------------------------------------------------------
// Dataset
// ---------------------------------------------------------------------------

/// A named, typed, shaped array of raw bytes with optional attributes.
#[derive(Debug, Clone)]
pub struct Dataset {
    pub name: String,
    pub dtype: DataType,
    pub shape: DatasetShape,
    /// Raw little-endian bytes of the dataset contents.
    pub data: Vec<u8>,
    pub attributes: Vec<Attribute>,
}

impl Dataset {
    /// Constructs a dataset from a slice of `f64` values.
    pub fn from_f64_slice(name: &str, data: &[f64], shape: DatasetShape) -> Self {
        let mut bytes = Vec::with_capacity(data.len() * 8);
        for &v in data {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        Self {
            name: name.to_string(),
            dtype: DataType::Float64,
            shape,
            data: bytes,
            attributes: Vec::new(),
        }
    }

    /// Constructs a dataset from a slice of `f32` values.
    pub fn from_f32_slice(name: &str, data: &[f32], shape: DatasetShape) -> Self {
        let mut bytes = Vec::with_capacity(data.len() * 4);
        for &v in data {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        Self {
            name: name.to_string(),
            dtype: DataType::Float32,
            shape,
            data: bytes,
            attributes: Vec::new(),
        }
    }

    /// Constructs a dataset from a slice of `i32` values.
    pub fn from_i32_slice(name: &str, data: &[i32], shape: DatasetShape) -> Self {
        let mut bytes = Vec::with_capacity(data.len() * 4);
        for &v in data {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        Self {
            name: name.to_string(),
            dtype: DataType::Int32,
            shape,
            data: bytes,
            attributes: Vec::new(),
        }
    }

    /// Decodes the raw bytes as consecutive `N`-byte little-endian elements.
    ///
    /// Fails if the shape does not describe exactly the stored bytes, which
    /// also covers shapes whose element count overflows `usize`.
    fn decode_elements<const N: usize, T>(
        &self,
        type_name: &str,
        decode: fn([u8; N]) -> T,
    ) -> Result<Vec<T>, String> {
        let mismatch = || {
            format!(
                "Data length mismatch: {} bytes for shape {:?} of {} elements",
                self.data.len(),
                self.shape.dims,
                type_name
            )
        };
        let expected = self
            .shape
            .checked_total_elements()
            .and_then(|n| n.checked_mul(N));
        if expected != Some(self.data.len()) {
            return Err(mismatch());
        }
        self.data
            .chunks_exact(N)
            .map(|chunk| <[u8; N]>::try_from(chunk).map(decode))
            .collect::<Result<Vec<T>, _>>()
            .map_err(|_| mismatch())
    }

    /// Decodes the raw bytes back into a `Vec`f64`.
    pub fn to_f64_vec(&self) -> Result<Vec<f64>, String> {
        if self.dtype != DataType::Float64 {
            return Err(format!("Expected Float64, got {:?}", self.dtype));
        }
        self.decode_elements::<8, f64>("f64", f64::from_le_bytes)
    }

    /// Decodes the raw bytes back into a `Vec`f32`.
    pub fn to_f32_vec(&self) -> Result<Vec<f32>, String> {
        if self.dtype != DataType::Float32 {
            return Err(format!("Expected Float32, got {:?}", self.dtype));
        }
        self.decode_elements::<4, f32>("f32", f32::from_le_bytes)
    }

    /// Decodes the raw bytes back into a `Vec`i32`.
    pub fn to_i32_vec(&self) -> Result<Vec<i32>, String> {
        if self.dtype != DataType::Int32 {
            return Err(format!("Expected Int32, got {:?}", self.dtype));
        }
        self.decode_elements::<4, i32>("i32", i32::from_le_bytes)
    }

    /// Appends an attribute to this dataset.
    pub fn add_attribute(&mut self, attr: Attribute) {
        self.attributes.push(attr);
    }

    /// Looks up an attribute by name.
    pub fn get_attribute(&self, name: &str) -> Option<&Attribute> {
        self.attributes.iter().find(|a| a.name == name)
    }
}

// ---------------------------------------------------------------------------
// Group
// ---------------------------------------------------------------------------

/// A named container that holds datasets, subgroups, and attributes.
#[derive(Debug, Clone)]
pub struct Group {
    pub name: String,
    pub datasets: Vec<Dataset>,
    pub subgroups: Vec<Group>,
    pub attributes: Vec<Attribute>,
}

impl Group {
    /// Creates an empty group with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            datasets: Vec::new(),
            subgroups: Vec::new(),
            attributes: Vec::new(),
        }
    }

    /// Adds a dataset to this group.
    pub fn add_dataset(&mut self, ds: Dataset) {
        self.datasets.push(ds);
    }

    /// Adds a subgroup to this group.
    pub fn add_subgroup(&mut self, g: Group) {
        self.subgroups.push(g);
    }

    /// Looks up a dataset by name.
    pub fn get_dataset(&self, name: &str) -> Option<&Dataset> {
        self.datasets.iter().find(|d| d.name == name)
    }

    /// Looks up a direct subgroup by name.
    pub fn get_subgroup(&self, name: &str) -> Option<&Group> {
        self.subgroups.iter().find(|g| g.name == name)
    }

    /// Appends an attribute to this group.
    pub fn add_attribute(&mut self, attr: Attribute) {
        self.attributes.push(attr);
    }

    /// Looks up an attribute by name.
    pub fn get_attribute(&self, name: &str) -> Option<&Attribute> {
        self.attributes.iter().find(|a| a.name == name)
    }
}

// ---------------------------------------------------------------------------
// OxiFile
// ---------------------------------------------------------------------------

const MAGIC: &[u8; 8] = b"OXIPHY01";

/// Root container for the OxiPhy binary file format.
///
/// Binary layout:
/// - 8 bytes magic: `OXIPHY01`
/// - 4 bytes version: `u32` little-endian
/// - Recursively serialized root group
#[derive(Debug, Clone)]
pub struct OxiFile {
    pub version: u32,
    pub root: Group,
}

impl OxiFile {
    /// Creates a new empty `OxiFile` at version 1.
    pub fn new() -> Self {
        Self {
            version: 1,
            root: Group::new("/"),
        }
    }

    /// Serializes the entire file to a byte vector.
    pub fn write_to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(MAGIC);
        write_u32(&mut buf, self.version);
        serialize_group(&mut buf, &self.root);
        buf
    }

    /// Deserializes an `OxiFile` from a byte slice.
    pub fn read_from_bytes(data: &[u8]) -> Result<Self, String> {
        if data.len() < 12 {
            return Err("Data too short to be a valid OxiFile".to_string());
        }
        if &data[0..8] != MAGIC {
            return Err("Invalid magic bytes: not an OxiFile".to_string());
        }
        let mut pos = 8usize;
        let version = read_u32(data, &mut pos)?;
        let root = deserialize_group(data, &mut pos, 0)?;
        Ok(Self { version, root })
    }

    /// Writes the file to disk at the given path.
    pub fn save(&self, path: &str) -> Result<(), String> {
        let bytes = self.write_to_bytes();
        let mut f =
            fs::File::create(path).map_err(|e| format!("Cannot create file '{}': {}", path, e))?;
        f.write_all(&bytes)
            .map_err(|e| format!("Write error: {}", e))?;
        Ok(())
    }

    /// Loads an `OxiFile` from disk.
    pub fn load(path: &str) -> Result<Self, String> {
        let bytes = fs::read(path).map_err(|e| format!("Cannot read file '{}': {}", path, e))?;
        Self::read_from_bytes(&bytes)
    }
}

impl Default for OxiFile {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Serialization helpers (public)
// ---------------------------------------------------------------------------

/// Appends a `u32` in little-endian format.
pub fn write_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

/// Appends a `u64` in little-endian format.
pub fn write_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

/// Appends a length-prefixed UTF-8 string (u32 length + bytes).
pub fn write_string(buf: &mut Vec<u8>, s: &str) {
    let bytes = s.as_bytes();
    write_u32(buf, bytes.len() as u32);
    buf.extend_from_slice(bytes);
}

/// Reads a `u32` in little-endian format, advancing `pos`.
pub fn read_u32(data: &[u8], pos: &mut usize) -> Result<u32, String> {
    Ok(u32::from_le_bytes(read_array::<4>(data, pos, "read_u32")?))
}

/// Reads a `u64` in little-endian format, advancing `pos`.
pub fn read_u64(data: &[u8], pos: &mut usize) -> Result<u64, String> {
    Ok(u64::from_le_bytes(read_array::<8>(data, pos, "read_u64")?))
}

/// Reads a length-prefixed UTF-8 string, advancing `pos`.
pub fn read_string(data: &[u8], pos: &mut usize) -> Result<String, String> {
    let len = read_u32(data, pos)? as usize;
    let start = *pos;
    let body = start
        .checked_add(len)
        .and_then(|end| data.get(start..end))
        .ok_or_else(|| format!("read_string: string body out of bounds at pos {start}"))?;
    let s = std::str::from_utf8(body)
        .map_err(|e| format!("Invalid UTF-8: {}", e))?
        .to_string();
    *pos = start + len;
    Ok(s)
}

// ---------------------------------------------------------------------------
// Internal serialization helpers
// ---------------------------------------------------------------------------

/// Error text shared by all bounds-checked readers.
fn eof_error(what: &str, pos: usize) -> String {
    format!("{what}: unexpected end of data at pos {pos}")
}

/// Reads exactly `N` bytes at `*pos` and advances `*pos`.
///
/// The end offset is computed with checked arithmetic, so a corrupt length
/// can never wrap around and defeat the bounds check.
fn read_array<const N: usize>(data: &[u8], pos: &mut usize, what: &str) -> Result<[u8; N], String> {
    let start = *pos;
    let bytes = start
        .checked_add(N)
        .and_then(|end| data.get(start..end))
        .ok_or_else(|| eof_error(what, start))?;
    let array = <[u8; N]>::try_from(bytes).map_err(|_| eof_error(what, start))?;
    *pos = start + N;
    Ok(array)
}

/// Largest element count that can still be present in `data[pos..]` when every
/// element occupies at least `min_item_bytes`.
///
/// Counts read from untrusted input must never be passed to
/// `Vec::with_capacity` directly; clamping them with this bound keeps a tiny
/// corrupt file from requesting gigabytes of memory.
fn bounded_capacity(count: usize, data: &[u8], pos: usize, min_item_bytes: usize) -> usize {
    count.min(data.len().saturating_sub(pos) / min_item_bytes.max(1))
}

fn write_i64(buf: &mut Vec<u8>, v: i64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn read_i64(data: &[u8], pos: &mut usize) -> Result<i64, String> {
    Ok(i64::from_le_bytes(read_array::<8>(data, pos, "read_i64")?))
}

fn write_f64(buf: &mut Vec<u8>, v: f64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn read_f64(data: &[u8], pos: &mut usize) -> Result<f64, String> {
    Ok(f64::from_le_bytes(read_array::<8>(data, pos, "read_f64")?))
}

// ---------------------------------------------------------------------------
// Attribute serialization
// ---------------------------------------------------------------------------

fn serialize_attribute(buf: &mut Vec<u8>, attr: &Attribute) {
    write_string(buf, &attr.name);
    buf.push(attr.value.tag());
    match &attr.value {
        AttributeValue::Int(v) => {
            write_i64(buf, *v);
        }
        AttributeValue::Float(v) => {
            write_f64(buf, *v);
        }
        AttributeValue::Text(s) => {
            write_string(buf, s);
        }
        AttributeValue::IntArray(arr) => {
            write_u64(buf, arr.len() as u64);
            for &v in arr {
                write_i64(buf, v);
            }
        }
        AttributeValue::FloatArray(arr) => {
            write_u64(buf, arr.len() as u64);
            for &v in arr {
                write_f64(buf, v);
            }
        }
    }
}

/// Reads an element count stored as `u64` and converts it to `usize`.
fn read_count_u64(data: &[u8], pos: &mut usize, what: &str) -> Result<usize, String> {
    usize::try_from(read_u64(data, pos)?)
        .map_err(|_| format!("{what}: count does not fit in usize"))
}

/// Smallest encoded attribute: name length prefix, type tag and an empty text.
const MIN_ATTRIBUTE_BYTES: usize = 4 + 1 + 4;

fn deserialize_attribute(data: &[u8], pos: &mut usize) -> Result<Attribute, String> {
    let name = read_string(data, pos)?;
    let tag = *data
        .get(*pos)
        .ok_or_else(|| "deserialize_attribute: missing type tag".to_string())?;
    *pos += 1;
    let value = match tag {
        0 => AttributeValue::Int(read_i64(data, pos)?),
        1 => AttributeValue::Float(read_f64(data, pos)?),
        2 => AttributeValue::Text(read_string(data, pos)?),
        3 => {
            let n = read_count_u64(data, pos, "deserialize_attribute")?;
            let mut arr = Vec::with_capacity(bounded_capacity(n, data, *pos, 8));
            for _ in 0..n {
                arr.push(read_i64(data, pos)?);
            }
            AttributeValue::IntArray(arr)
        }
        4 => {
            let n = read_count_u64(data, pos, "deserialize_attribute")?;
            let mut arr = Vec::with_capacity(bounded_capacity(n, data, *pos, 8));
            for _ in 0..n {
                arr.push(read_f64(data, pos)?);
            }
            AttributeValue::FloatArray(arr)
        }
        _ => return Err(format!("Unknown AttributeValue tag: {}", tag)),
    };
    Ok(Attribute { name, value })
}

// ---------------------------------------------------------------------------
// Dataset serialization
// ---------------------------------------------------------------------------

fn serialize_dataset(buf: &mut Vec<u8>, ds: &Dataset) {
    write_string(buf, &ds.name);
    buf.push(ds.dtype.tag());
    // shape
    write_u32(buf, ds.shape.dims.len() as u32);
    for &d in &ds.shape.dims {
        write_u64(buf, d as u64);
    }
    // raw data
    write_u64(buf, ds.data.len() as u64);
    buf.extend_from_slice(&ds.data);
    // attributes
    write_u32(buf, ds.attributes.len() as u32);
    for attr in &ds.attributes {
        serialize_attribute(buf, attr);
    }
}

/// Smallest encoded dataset: name length, dtype tag, dimension count, data
/// length and attribute count.
const MIN_DATASET_BYTES: usize = 4 + 1 + 4 + 8 + 4;

fn deserialize_dataset(data: &[u8], pos: &mut usize) -> Result<Dataset, String> {
    let name = read_string(data, pos)?;
    let tag = *data
        .get(*pos)
        .ok_or_else(|| "deserialize_dataset: missing dtype tag".to_string())?;
    let dtype = DataType::from_tag(tag)?;
    *pos += 1;
    let ndims = read_u32(data, pos)? as usize;
    let mut dims = Vec::with_capacity(bounded_capacity(ndims, data, *pos, 8));
    for _ in 0..ndims {
        dims.push(read_count_u64(data, pos, "deserialize_dataset")?);
    }
    let shape = DatasetShape { dims };
    let data_len = read_count_u64(data, pos, "deserialize_dataset")?;
    let body_start = *pos;
    let body = body_start
        .checked_add(data_len)
        .and_then(|end| data.get(body_start..end))
        .ok_or_else(|| {
            format!("deserialize_dataset: data body out of bounds at pos {body_start}")
        })?;
    let expected_len = shape
        .checked_total_elements()
        .and_then(|elements| elements.checked_mul(dtype.size_bytes()));
    if expected_len != Some(data_len) {
        return Err(format!(
            "deserialize_dataset: shape {:?} of {:?} elements does not match data length {}",
            shape.dims, dtype, data_len
        ));
    }
    let raw = body.to_vec();
    *pos = body_start + data_len;
    let n_attrs = read_u32(data, pos)? as usize;
    let mut attributes =
        Vec::with_capacity(bounded_capacity(n_attrs, data, *pos, MIN_ATTRIBUTE_BYTES));
    for _ in 0..n_attrs {
        attributes.push(deserialize_attribute(data, pos)?);
    }
    Ok(Dataset {
        name,
        dtype,
        shape,
        data: raw,
        attributes,
    })
}

// ---------------------------------------------------------------------------
// Group serialization
// ---------------------------------------------------------------------------

fn serialize_group(buf: &mut Vec<u8>, group: &Group) {
    write_string(buf, &group.name);
    // attributes
    write_u32(buf, group.attributes.len() as u32);
    for attr in &group.attributes {
        serialize_attribute(buf, attr);
    }
    // datasets
    write_u32(buf, group.datasets.len() as u32);
    for ds in &group.datasets {
        serialize_dataset(buf, ds);
    }
    // subgroups (recursive)
    write_u32(buf, group.subgroups.len() as u32);
    for sg in &group.subgroups {
        serialize_group(buf, sg);
    }
}

/// Maximum group nesting depth accepted by the decoder (the root group is at
/// depth 0).  An empty nested group costs only 16 bytes, so without a bound a
/// small corrupt file would overflow the stack through recursion.
const MAX_GROUP_DEPTH: u32 = 64;

/// Smallest encoded group: name length, attribute count, dataset count and
/// subgroup count.
const MIN_GROUP_BYTES: usize = 4 + 4 + 4 + 4;

fn deserialize_group(data: &[u8], pos: &mut usize, depth: u32) -> Result<Group, String> {
    if depth > MAX_GROUP_DEPTH {
        return Err(format!(
            "deserialize_group: group nesting too deep (limit {MAX_GROUP_DEPTH})"
        ));
    }
    let name = read_string(data, pos)?;
    let n_attrs = read_u32(data, pos)? as usize;
    let mut attributes =
        Vec::with_capacity(bounded_capacity(n_attrs, data, *pos, MIN_ATTRIBUTE_BYTES));
    for _ in 0..n_attrs {
        attributes.push(deserialize_attribute(data, pos)?);
    }
    let n_datasets = read_u32(data, pos)? as usize;
    let mut datasets =
        Vec::with_capacity(bounded_capacity(n_datasets, data, *pos, MIN_DATASET_BYTES));
    for _ in 0..n_datasets {
        datasets.push(deserialize_dataset(data, pos)?);
    }
    let n_subgroups = read_u32(data, pos)? as usize;
    let mut subgroups =
        Vec::with_capacity(bounded_capacity(n_subgroups, data, *pos, MIN_GROUP_BYTES));
    for _ in 0..n_subgroups {
        subgroups.push(deserialize_group(data, pos, depth + 1)?);
    }
    Ok(Group {
        name,
        datasets,
        subgroups,
        attributes,
    })
}

// ---------------------------------------------------------------------------
// SimulationCheckpoint
// ---------------------------------------------------------------------------

/// High-level helper for writing physics simulation checkpoint data.
pub struct SimulationCheckpoint;

impl SimulationCheckpoint {
    /// Creates a fresh `OxiFile` ready for checkpoint data.
    pub fn create() -> OxiFile {
        OxiFile::new()
    }

    /// Flattens a slice of 3-vectors and stores them as a `(N, 3)` dataset.
    pub fn add_positions(file: &mut OxiFile, group: &str, positions: &[[f64; 3]]) {
        let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        let shape = DatasetShape::new(vec![positions.len(), 3]);
        let ds = Dataset::from_f64_slice("positions", &flat, shape);
        Self::get_or_create_group(&mut file.root, group).add_dataset(ds);
    }

    /// Flattens a slice of 3-vectors and stores them as a `(N, 3)` dataset.
    pub fn add_velocities(file: &mut OxiFile, group: &str, velocities: &[[f64; 3]]) {
        let flat: Vec<f64> = velocities.iter().flat_map(|v| v.iter().copied()).collect();
        let shape = DatasetShape::new(vec![velocities.len(), 3]);
        let ds = Dataset::from_f64_slice("velocities", &flat, shape);
        Self::get_or_create_group(&mut file.root, group).add_dataset(ds);
    }

    /// Stores a 1-D scalar field as a `(N,)` dataset.
    pub fn add_scalar_field(file: &mut OxiFile, group: &str, name: &str, values: &[f64]) {
        let shape = DatasetShape::new(vec![values.len()]);
        let ds = Dataset::from_f64_slice(name, values, shape);
        Self::get_or_create_group(&mut file.root, group).add_dataset(ds);
    }

    /// Stores timestep metadata as attributes on the root group.
    pub fn add_timestep_metadata(file: &mut OxiFile, step: u64, time: f64, dt: f64) {
        file.root
            .add_attribute(Attribute::new("step", AttributeValue::Int(step as i64)));
        file.root
            .add_attribute(Attribute::new("time", AttributeValue::Float(time)));
        file.root
            .add_attribute(Attribute::new("dt", AttributeValue::Float(dt)));
    }

    // Returns a mutable reference to the named direct subgroup, creating it if absent.
    fn get_or_create_group<'a>(root: &'a mut Group, name: &str) -> &'a mut Group {
        let idx = match root.subgroups.iter().position(|g| g.name == name) {
            Some(idx) => idx,
            None => {
                root.subgroups.push(Group::new(name));
                root.subgroups.len() - 1
            }
        };
        &mut root.subgroups[idx]
    }
}

impl Default for SimulationCheckpoint {
    fn default() -> Self {
        Self
    }
}

// ---------------------------------------------------------------------------
// Endianness handling
// ---------------------------------------------------------------------------

/// Which byte order a binary stream uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Endianness {
    /// Least-significant byte first (x86/arm default).
    Little,
    /// Most-significant byte first (network order, big-endian).
    Big,
}

impl Endianness {
    /// Detect the native byte order at runtime.
    pub fn native() -> Self {
        if cfg!(target_endian = "little") {
            Endianness::Little
        } else {
            Endianness::Big
        }
    }

    /// Convert a `u32` to bytes in this byte order.
    pub fn u32_to_bytes(self, v: u32) -> [u8; 4] {
        match self {
            Endianness::Little => v.to_le_bytes(),
            Endianness::Big => v.to_be_bytes(),
        }
    }

    /// Convert a `u32` from bytes in this byte order.
    pub fn u32_from_bytes(self, b: [u8; 4]) -> u32 {
        match self {
            Endianness::Little => u32::from_le_bytes(b),
            Endianness::Big => u32::from_be_bytes(b),
        }
    }

    /// Convert a `f64` to bytes in this byte order.
    pub fn f64_to_bytes(self, v: f64) -> [u8; 8] {
        match self {
            Endianness::Little => v.to_le_bytes(),
            Endianness::Big => v.to_be_bytes(),
        }
    }

    /// Convert a `f64` from bytes in this byte order.
    pub fn f64_from_bytes(self, b: [u8; 8]) -> f64 {
        match self {
            Endianness::Little => f64::from_le_bytes(b),
            Endianness::Big => f64::from_be_bytes(b),
        }
    }

    /// Swap bytes of a `u32` from this endianness to native order.
    pub fn u32_to_native(self, v: u32) -> u32 {
        match self {
            Endianness::Little => u32::from_le_bytes(v.to_ne_bytes()),
            Endianness::Big => u32::from_be_bytes(v.to_ne_bytes()),
        }
    }
}

// ---------------------------------------------------------------------------
// Binary mesh format (simple triangle mesh)
// ---------------------------------------------------------------------------

/// A minimal binary triangle mesh representation.
///
/// Format:
/// - `u32` number of vertices
/// - `u32` number of triangles
/// - `3 * n_verts * f64` flat XYZ vertex positions
/// - `3 * n_tris * u32` triangle vertex indices
#[derive(Debug, Clone)]
pub struct BinaryMesh {
    /// Vertex positions stored as flat XYZ triples.
    pub vertices: Vec<[f64; 3]>,
    /// Triangle indices (0-based).
    pub triangles: Vec<[u32; 3]>,
}

impl BinaryMesh {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            triangles: Vec::new(),
        }
    }

    /// Serialize to bytes (little-endian).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        write_u32(&mut buf, self.vertices.len() as u32);
        write_u32(&mut buf, self.triangles.len() as u32);
        for v in &self.vertices {
            for component in v.iter() {
                buf.extend_from_slice(&component.to_le_bytes());
            }
        }
        for t in &self.triangles {
            for component in t.iter() {
                write_u32(&mut buf, *component);
            }
        }
        buf
    }

    /// Deserialize from bytes (little-endian).
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        let mut pos = 0usize;
        let n_verts = read_u32(data, &mut pos)? as usize;
        let n_tris = read_u32(data, &mut pos)? as usize;

        let mut vertices = Vec::with_capacity(bounded_capacity(n_verts, data, pos, 24));
        for _ in 0..n_verts {
            let mut xyz = [0.0_f64; 3];
            for component in xyz.iter_mut() {
                *component = read_f64(data, &mut pos)
                    .map_err(|_| "BinaryMesh: vertex data truncated".to_string())?;
            }
            vertices.push(xyz);
        }

        let mut triangles = Vec::with_capacity(bounded_capacity(n_tris, data, pos, 12));
        for _ in 0..n_tris {
            let i0 = read_u32(data, &mut pos)?;
            let i1 = read_u32(data, &mut pos)?;
            let i2 = read_u32(data, &mut pos)?;
            triangles.push([i0, i1, i2]);
        }

        Ok(Self {
            vertices,
            triangles,
        })
    }

    /// Number of vertices.
    pub fn n_vertices(&self) -> usize {
        self.vertices.len()
    }

    /// Number of triangles.
    pub fn n_triangles(&self) -> usize {
        self.triangles.len()
    }
}

impl Default for BinaryMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Binary particle data
// ---------------------------------------------------------------------------

/// Compact binary storage for particle data (positions + optional scalars).
///
/// Header: magic `b"OXIPART"` + u32 particle count + u32 field count.
/// Data: `n_particles * (3 + n_fields)` f64 values in row-major order.
#[derive(Debug, Clone)]
pub struct BinaryParticleData {
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Optional scalar fields (one Vec per field, each of length `n`).
    pub scalar_fields: Vec<Vec<f64>>,
    /// Field names corresponding to `scalar_fields`.
    pub field_names: Vec<String>,
}

const PARTICLE_MAGIC: &[u8; 7] = b"OXIPART";

impl BinaryParticleData {
    /// Create a new empty particle data container.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            scalar_fields: Vec::new(),
            field_names: Vec::new(),
        }
    }

    /// Add a scalar field by name.  Must have one value per particle.
    pub fn add_field(&mut self, name: &str, values: Vec<f64>) {
        assert_eq!(
            values.len(),
            self.positions.len(),
            "Field '{}' length {} != particle count {}",
            name,
            values.len(),
            self.positions.len()
        );
        self.field_names.push(name.to_string());
        self.scalar_fields.push(values);
    }

    /// Number of particles.
    pub fn n_particles(&self) -> usize {
        self.positions.len()
    }

    /// Number of scalar fields.
    pub fn n_fields(&self) -> usize {
        self.scalar_fields.len()
    }

    /// Serialize to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let n = self.positions.len();
        let nf = self.scalar_fields.len();
        let mut buf = Vec::new();
        buf.extend_from_slice(PARTICLE_MAGIC);
        write_u32(&mut buf, n as u32);
        write_u32(&mut buf, nf as u32);
        // Field names
        for name in &self.field_names {
            write_string(&mut buf, name);
        }
        // Positions
        for p in &self.positions {
            for component in p.iter() {
                buf.extend_from_slice(&component.to_le_bytes());
            }
        }
        // Scalar fields
        for field in &self.scalar_fields {
            for &v in field {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        buf
    }

    /// Deserialize from bytes.
    pub fn from_bytes(data: &[u8]) -> Result<Self, String> {
        if data.len() < 7 {
            return Err("BinaryParticleData: too short".to_string());
        }
        if &data[..7] != PARTICLE_MAGIC {
            return Err("BinaryParticleData: bad magic".to_string());
        }
        let mut pos = 7usize;
        let n = read_u32(data, &mut pos)? as usize;
        let nf = read_u32(data, &mut pos)? as usize;

        let mut field_names = Vec::with_capacity(bounded_capacity(nf, data, pos, 4));
        for _ in 0..nf {
            field_names.push(read_string(data, &mut pos)?);
        }

        let mut positions = Vec::with_capacity(bounded_capacity(n, data, pos, 24));
        for _ in 0..n {
            let mut xyz = [0.0_f64; 3];
            for component in xyz.iter_mut() {
                *component = read_f64(data, &mut pos)
                    .map_err(|_| "BinaryParticleData: positions truncated".to_string())?;
            }
            positions.push(xyz);
        }

        let mut scalar_fields = Vec::with_capacity(bounded_capacity(nf, data, pos, 8));
        for _ in 0..nf {
            let mut field = Vec::with_capacity(bounded_capacity(n, data, pos, 8));
            for _ in 0..n {
                let v = read_f64(data, &mut pos)
                    .map_err(|_| "BinaryParticleData: scalar field truncated".to_string())?;
                field.push(v);
            }
            scalar_fields.push(field);
        }

        Ok(Self {
            positions,
            scalar_fields,
            field_names,
        })
    }
}

impl Default for BinaryParticleData {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Compressed binary output (run-length encoding for f64 streams)
// ---------------------------------------------------------------------------

/// Very simple run-length compression for f64 arrays.
///
/// Consecutive equal values are stored as `(value, count)` pairs.
/// This is mainly useful for fields with large constant regions.
///
/// Format: `u32 n_runs || \[f64 value, u32 count\] × n_runs`
pub fn rle_compress_f64(values: &[f64]) -> Vec<u8> {
    if values.is_empty() {
        let mut buf = Vec::new();
        write_u32(&mut buf, 0);
        return buf;
    }

    let mut runs: Vec<(f64, u32)> = Vec::new();
    let mut cur = values[0];
    let mut cnt = 1u32;
    for &v in &values[1..] {
        if v.to_bits() == cur.to_bits() {
            cnt += 1;
        } else {
            runs.push((cur, cnt));
            cur = v;
            cnt = 1;
        }
    }
    runs.push((cur, cnt));

    let mut buf = Vec::new();
    write_u32(&mut buf, runs.len() as u32);
    for (val, count) in runs {
        buf.extend_from_slice(&val.to_le_bytes());
        write_u32(&mut buf, count);
    }
    buf
}

/// Decompress a run-length encoded f64 array.
pub fn rle_decompress_f64(data: &[u8]) -> Result<Vec<f64>, String> {
    let mut pos = 0usize;
    let n_runs = read_u32(data, &mut pos)? as usize;
    let mut result = Vec::new();
    for _ in 0..n_runs {
        let truncated = |_: String| "rle_decompress_f64: truncated".to_string();
        let val = read_f64(data, &mut pos).map_err(truncated)?;
        let count = read_u32(data, &mut pos).map_err(truncated)? as usize;
        for _ in 0..count {
            result.push(val);
        }
    }
    Ok(result)
}

// ---------------------------------------------------------------------------
// Binary format versioning
// ---------------------------------------------------------------------------

/// Supported OxiFile format versions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FormatVersion {
    /// Version 1: initial release.
    V1 = 1,
    /// Version 2: added compression flags.
    V2 = 2,
}

impl FormatVersion {
    /// Parse from a `u32` tag.
    pub fn from_u32(v: u32) -> Result<Self, String> {
        match v {
            1 => Ok(FormatVersion::V1),
            2 => Ok(FormatVersion::V2),
            _ => Err(format!("Unknown format version: {v}")),
        }
    }

    /// Convert to a `u32` tag.
    pub fn to_u32(self) -> u32 {
        self as u32
    }

    /// Whether this version supports compression metadata.
    pub fn supports_compression(self) -> bool {
        matches!(self, FormatVersion::V2)
    }
}

// ---------------------------------------------------------------------------
// Extended OxiFile helpers
// ---------------------------------------------------------------------------

impl OxiFile {
    /// Create an OxiFile at a specific format version.
    pub fn with_version(version: FormatVersion) -> Self {
        Self {
            version: version.to_u32(),
            root: Group::new("/"),
        }
    }

    /// Return the parsed `FormatVersion` if recognised, else an error.
    pub fn format_version(&self) -> Result<FormatVersion, String> {
        FormatVersion::from_u32(self.version)
    }

    /// Store a `BinaryMesh` under the given group name.
    pub fn add_binary_mesh(&mut self, group: &str, mesh: &BinaryMesh) {
        let bytes = mesh.to_bytes();
        let shape = DatasetShape::new(vec![bytes.len()]);
        let mut ds = Dataset {
            name: "mesh_binary".to_string(),
            dtype: DataType::UInt8,
            shape,
            data: bytes,
            attributes: Vec::new(),
        };
        ds.add_attribute(Attribute::new(
            "n_vertices",
            AttributeValue::Int(mesh.n_vertices() as i64),
        ));
        ds.add_attribute(Attribute::new(
            "n_triangles",
            AttributeValue::Int(mesh.n_triangles() as i64),
        ));
        SimulationCheckpoint::get_or_create_group(&mut self.root, group).add_dataset(ds);
    }

    /// Retrieve and decode a `BinaryMesh` from the given group.
    pub fn get_binary_mesh(&self, group: &str) -> Result<BinaryMesh, String> {
        let grp = self
            .root
            .get_subgroup(group)
            .ok_or_else(|| format!("Group '{}' not found", group))?;
        let ds = grp
            .get_dataset("mesh_binary")
            .ok_or_else(|| "Dataset 'mesh_binary' not found".to_string())?;
        BinaryMesh::from_bytes(&ds.data)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dataset_f64_round_trip() {
        let original = vec![1.0_f64, 2.5, -3.125, 0.0, 1e10];
        let shape = DatasetShape::new(vec![original.len()]);
        let ds = Dataset::from_f64_slice("test", &original, shape);
        let recovered = ds.to_f64_vec().expect("to_f64_vec failed");
        assert_eq!(original, recovered);
    }

    #[test]
    fn test_dataset_f32_round_trip() {
        let original = vec![1.0_f32, 2.5, -3.125, 0.0];
        let shape = DatasetShape::new(vec![original.len()]);
        let ds = Dataset::from_f32_slice("f32ds", &original, shape);
        let recovered = ds.to_f32_vec().expect("to_f32_vec failed");
        assert_eq!(original, recovered);
    }

    #[test]
    fn test_dataset_i32_round_trip() {
        let original = vec![0_i32, -1, 42, i32::MAX, i32::MIN];
        let shape = DatasetShape::new(vec![original.len()]);
        let ds = Dataset::from_i32_slice("i32ds", &original, shape);
        let recovered = ds.to_i32_vec().expect("to_i32_vec failed");
        assert_eq!(original, recovered);
    }

    #[test]
    fn test_group_add_get_dataset() {
        let mut g = Group::new("particles");
        let ds = Dataset::from_f64_slice("energy", &[1.0, 2.0, 3.0], DatasetShape::new(vec![3]));
        g.add_dataset(ds);
        let found = g.get_dataset("energy").expect("dataset not found");
        assert_eq!(found.name, "energy");
        assert!(g.get_dataset("missing").is_none());
    }

    #[test]
    fn test_oxifile_round_trip() {
        let mut file = OxiFile::new();
        let ds = Dataset::from_f64_slice("x", &[1.0, 2.0, 3.0], DatasetShape::new(vec![3]));
        file.root.add_dataset(ds);

        let bytes = file.write_to_bytes();
        let loaded = OxiFile::read_from_bytes(&bytes).expect("round-trip failed");
        assert_eq!(loaded.version, 1);
        let ds2 = loaded
            .root
            .get_dataset("x")
            .expect("dataset missing after round-trip");
        let vals = ds2.to_f64_vec().unwrap();
        assert_eq!(vals, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_read_from_bytes_invalid_magic() {
        let bad: Vec<u8> = b"BADMAGIC\x01\x00\x00\x00".to_vec();
        let result = OxiFile::read_from_bytes(&bad);
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("Invalid magic bytes"));
    }

    #[test]
    fn test_dataset_shape_total_elements() {
        let s = DatasetShape::new(vec![3, 4, 5]);
        assert_eq!(s.total_elements(), 60);
        assert_eq!(s.rank(), 3);
        assert!(!s.is_scalar());

        let scalar = DatasetShape::new(vec![]);
        assert_eq!(scalar.total_elements(), 1);
        assert!(scalar.is_scalar());
        assert_eq!(scalar.rank(), 0);
    }

    #[test]
    fn test_attribute_get_set() {
        let mut ds = Dataset::from_f64_slice("d", &[1.0], DatasetShape::new(vec![1]));
        ds.add_attribute(Attribute::new(
            "units",
            AttributeValue::Text("meters".to_string()),
        ));
        ds.add_attribute(Attribute::new("count", AttributeValue::Int(42)));

        let attr = ds.get_attribute("units").expect("units not found");
        assert_eq!(attr.value, AttributeValue::Text("meters".to_string()));
        assert!(ds.get_attribute("nope").is_none());
    }

    #[test]
    fn test_simulation_checkpoint_positions() {
        let mut file = SimulationCheckpoint::create();
        let positions = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        SimulationCheckpoint::add_positions(&mut file, "frame0", &positions);

        let grp = file.root.get_subgroup("frame0").expect("group missing");
        let ds = grp.get_dataset("positions").expect("positions missing");
        let vals = ds.to_f64_vec().expect("to_f64_vec failed");

        let expected: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        assert_eq!(vals, expected);
        assert_eq!(ds.shape.dims, vec![3, 3]);
    }

    #[test]
    fn test_simulation_checkpoint_round_trip() {
        let mut file = SimulationCheckpoint::create();
        let positions = [[0.1, 0.2, 0.3], [-1.0, 2.0, -3.0]];
        SimulationCheckpoint::add_positions(&mut file, "step1", &positions);
        SimulationCheckpoint::add_timestep_metadata(&mut file, 1, 0.01, 0.001);

        let bytes = file.write_to_bytes();
        let loaded = OxiFile::read_from_bytes(&bytes).expect("round-trip failed");

        let step_attr = loaded.root.get_attribute("step").expect("step missing");
        assert_eq!(step_attr.value, AttributeValue::Int(1));

        let grp = loaded.root.get_subgroup("step1").expect("subgroup missing");
        let ds = grp.get_dataset("positions").expect("dataset missing");
        let vals = ds.to_f64_vec().unwrap();
        let expected: Vec<f64> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        assert_eq!(vals, expected);
    }

    // -----------------------------------------------------------------------
    // Endianness tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_endianness_u32_round_trip() {
        let v: u32 = 0xDEAD_BEEF;
        for end in [Endianness::Little, Endianness::Big] {
            let bytes = end.u32_to_bytes(v);
            let back = end.u32_from_bytes(bytes);
            assert_eq!(back, v, "Endianness {:?} u32 round-trip failed", end);
        }
    }

    #[test]
    fn test_endianness_f64_round_trip() {
        let v = std::f64::consts::PI;
        for end in [Endianness::Little, Endianness::Big] {
            let bytes = end.f64_to_bytes(v);
            let back = end.f64_from_bytes(bytes);
            assert!(
                (back - v).abs() < 1e-15,
                "Endianness {:?} f64 round-trip failed",
                end
            );
        }
    }

    #[test]
    fn test_endianness_native() {
        let native = Endianness::native();
        assert!(native == Endianness::Little || native == Endianness::Big);
    }

    // -----------------------------------------------------------------------
    // BinaryMesh tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_binary_mesh_round_trip() {
        let mut mesh = BinaryMesh::new();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        mesh.triangles = vec![[0, 1, 2]];
        let bytes = mesh.to_bytes();
        let mesh2 = BinaryMesh::from_bytes(&bytes).expect("round-trip failed");
        assert_eq!(mesh2.n_vertices(), 3);
        assert_eq!(mesh2.n_triangles(), 1);
        assert!((mesh2.vertices[1][0] - 1.0).abs() < 1e-15);
        assert_eq!(mesh2.triangles[0], [0, 1, 2]);
    }

    #[test]
    fn test_binary_mesh_empty() {
        let mesh = BinaryMesh::new();
        let bytes = mesh.to_bytes();
        let mesh2 = BinaryMesh::from_bytes(&bytes).expect("empty mesh round-trip failed");
        assert_eq!(mesh2.n_vertices(), 0);
        assert_eq!(mesh2.n_triangles(), 0);
    }

    // -----------------------------------------------------------------------
    // BinaryParticleData tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_binary_particle_data_round_trip() {
        let mut pd = BinaryParticleData::new();
        pd.positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        pd.add_field("density", vec![1000.0, 1200.0]);
        pd.add_field("pressure", vec![101325.0, 202650.0]);

        let bytes = pd.to_bytes();
        let pd2 = BinaryParticleData::from_bytes(&bytes).expect("round-trip failed");
        assert_eq!(pd2.n_particles(), 2);
        assert_eq!(pd2.n_fields(), 2);
        assert_eq!(pd2.field_names[0], "density");
        assert!((pd2.scalar_fields[0][1] - 1200.0).abs() < 1e-12);
        assert!((pd2.positions[1][2] - 6.0).abs() < 1e-15);
    }

    #[test]
    fn test_binary_particle_data_bad_magic() {
        let bad: Vec<u8> = b"BADMAGIC".to_vec();
        assert!(BinaryParticleData::from_bytes(&bad).is_err());
    }

    // -----------------------------------------------------------------------
    // RLE compression tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_rle_compress_decompress_round_trip() {
        let original = vec![1.0, 1.0, 1.0, 2.5, 2.5, 3.0, 1.0];
        let compressed = rle_compress_f64(&original);
        let decompressed = rle_decompress_f64(&compressed).expect("decompression failed");
        assert_eq!(original.len(), decompressed.len());
        for (a, b) in original.iter().zip(decompressed.iter()) {
            assert!((a - b).abs() < 1e-15);
        }
    }

    #[test]
    fn test_rle_compress_empty() {
        let compressed = rle_compress_f64(&[]);
        let decompressed = rle_decompress_f64(&compressed).expect("empty decompression failed");
        assert!(decompressed.is_empty());
    }

    #[test]
    fn test_rle_compresses_constant_field() {
        let original = vec![3.125; 1000];
        let compressed = rle_compress_f64(&original);
        // Should be much smaller than the raw 8000 bytes.
        assert!(
            compressed.len() < 100,
            "RLE should compress constant field significantly"
        );
        let decompressed = rle_decompress_f64(&compressed).unwrap();
        assert_eq!(decompressed.len(), 1000);
    }

    // -----------------------------------------------------------------------
    // FormatVersion tests
    // -----------------------------------------------------------------------

    #[test]
    fn test_format_version_round_trip() {
        assert_eq!(FormatVersion::from_u32(1).unwrap(), FormatVersion::V1);
        assert_eq!(FormatVersion::from_u32(2).unwrap(), FormatVersion::V2);
        assert!(FormatVersion::from_u32(99).is_err());
    }

    #[test]
    fn test_format_version_supports_compression() {
        assert!(!FormatVersion::V1.supports_compression());
        assert!(FormatVersion::V2.supports_compression());
    }

    // -----------------------------------------------------------------------
    // OxiFile + BinaryMesh integration test
    // -----------------------------------------------------------------------

    #[test]
    fn test_oxifile_binary_mesh_store_retrieve() {
        let mut file = OxiFile::new();
        let mut mesh = BinaryMesh::new();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        mesh.triangles = vec![[0, 1, 2]];
        file.add_binary_mesh("geometry", &mesh);

        let mesh2 = file
            .get_binary_mesh("geometry")
            .expect("retrieve mesh failed");
        assert_eq!(mesh2.n_vertices(), 3);
        assert_eq!(mesh2.n_triangles(), 1);
    }

    #[test]
    fn test_oxifile_with_version() {
        let file = OxiFile::with_version(FormatVersion::V2);
        assert_eq!(file.version, 2);
        assert_eq!(file.format_version().unwrap(), FormatVersion::V2);
    }

    #[test]
    fn test_datatype_size_bytes() {
        assert_eq!(DataType::Float32.size_bytes(), 4);
        assert_eq!(DataType::Float64.size_bytes(), 8);
        assert_eq!(DataType::Int32.size_bytes(), 4);
        assert_eq!(DataType::Int64.size_bytes(), 8);
        assert_eq!(DataType::UInt8.size_bytes(), 1);
        assert_eq!(DataType::UInt32.size_bytes(), 4);
    }

    #[test]
    fn test_dataset_wrong_type_returns_err() {
        let ds = Dataset::from_f64_slice("d", &[1.0, 2.0], DatasetShape::new(vec![2]));
        assert!(ds.to_f32_vec().is_err());
        assert!(ds.to_i32_vec().is_err());
    }

    // -----------------------------------------------------------------------
    // Malformed-input hardening tests
    // -----------------------------------------------------------------------

    /// Magic bytes plus format version 1.
    fn oxifile_header() -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(MAGIC);
        write_u32(&mut buf, 1);
        buf
    }

    /// Header, root group name and the root attribute count.
    fn root_prefix(attribute_count: u32) -> Vec<u8> {
        let mut buf = oxifile_header();
        write_string(&mut buf, "");
        write_u32(&mut buf, attribute_count);
        buf
    }

    /// A file whose single leaf group sits `depth` levels below the root.
    ///
    /// Every ancestor is an empty group with one child and encodes in 16 bytes:
    /// name length, attribute count, dataset count and subgroup count.
    fn nested_oxifile_bytes(depth: usize) -> Vec<u8> {
        let mut buf = oxifile_header();
        for subgroups in std::iter::repeat_n(1u32, depth).chain(std::iter::once(0)) {
            write_string(&mut buf, "");
            write_u32(&mut buf, 0); // attributes
            write_u32(&mut buf, 0); // datasets
            write_u32(&mut buf, subgroups);
        }
        buf
    }

    /// Root group announcing one dataset, encoded up to and including the data
    /// length; the caller appends the data bytes and the trailing counts.
    fn dataset_prefix(dtype_tag: u8, dims: &[u64], data_len: u64) -> Vec<u8> {
        let mut buf = root_prefix(0);
        write_u32(&mut buf, 1); // datasets
        write_string(&mut buf, "d");
        buf.push(dtype_tag);
        write_u32(&mut buf, dims.len() as u32);
        for &dim in dims {
            write_u64(&mut buf, dim);
        }
        write_u64(&mut buf, data_len);
        buf
    }

    /// A file exercising every record type, used by the corruption tests.
    fn rich_oxifile() -> OxiFile {
        let mut file = OxiFile::new();
        file.root
            .add_attribute(Attribute::new("step", AttributeValue::Int(3)));
        file.root
            .add_attribute(Attribute::new("dt", AttributeValue::Float(0.5)));
        file.root.add_attribute(Attribute::new(
            "label",
            AttributeValue::Text("run".to_string()),
        ));
        file.root.add_attribute(Attribute::new(
            "ints",
            AttributeValue::IntArray(vec![1, -2, 3]),
        ));
        file.root.add_attribute(Attribute::new(
            "floats",
            AttributeValue::FloatArray(vec![0.25, 0.5]),
        ));
        let mut ds = Dataset::from_f64_slice("x", &[1.0, 2.0, 3.0], DatasetShape::new(vec![3]));
        ds.add_attribute(Attribute::new("units", AttributeValue::Text("m".into())));
        file.root.add_dataset(ds);
        let mut inner = Group::new("inner");
        inner.add_dataset(Dataset::from_i32_slice(
            "idx",
            &[1, 2, 3, 4],
            DatasetShape::new(vec![2, 2]),
        ));
        inner.add_subgroup(Group::new("leaf"));
        file.root.add_subgroup(inner);
        file
    }

    #[test]
    fn oxifile_deep_nesting_returns_err() {
        // 100 000 nested groups cost only 1.6 MB; unbounded recursion would
        // overflow the stack.
        let err = OxiFile::read_from_bytes(&nested_oxifile_bytes(100_000))
            .expect_err("deeply nested file must be rejected");
        assert!(err.contains("nesting too deep"), "unexpected error: {err}");
    }

    #[test]
    fn oxifile_nesting_limit_boundary() {
        let limit = MAX_GROUP_DEPTH as usize;
        let ok = OxiFile::read_from_bytes(&nested_oxifile_bytes(limit))
            .expect("nesting at the limit must decode");
        let mut depth = 0usize;
        let mut group = &ok.root;
        while let Some(child) = group.subgroups.first() {
            depth += 1;
            group = child;
        }
        assert_eq!(depth, limit);
        assert!(OxiFile::read_from_bytes(&nested_oxifile_bytes(limit + 1)).is_err());
    }

    #[test]
    fn oxifile_huge_data_len_returns_err() {
        // A data length close to usize::MAX used to wrap `pos + data_len`
        // past the bounds check and then panic on the slice.
        for data_len in [u64::MAX, u64::MAX - 3, u64::MAX / 2, 1u64 << 40] {
            let mut buf = dataset_prefix(1, &[1], data_len);
            buf.extend_from_slice(&[0u8; 16]);
            let err =
                OxiFile::read_from_bytes(&buf).expect_err("oversized data length must be rejected");
            assert!(
                err.contains("out of bounds") || err.contains("does not fit"),
                "data_len {data_len}: unexpected error {err}"
            );
        }
    }

    #[test]
    fn oxifile_shape_must_match_data_length() {
        // Four f64 elements need 32 bytes, but only 8 are stored.
        let mut buf = dataset_prefix(1, &[4], 8);
        buf.extend_from_slice(&[0u8; 8]);
        write_u32(&mut buf, 0); // dataset attributes
        write_u32(&mut buf, 0); // root subgroups
        let err = OxiFile::read_from_bytes(&buf).expect_err("shape/data mismatch");
        assert!(err.contains("does not match"), "unexpected error: {err}");

        // A dimension product that overflows usize is a mismatch, not a panic.
        let mut buf = dataset_prefix(1, &[u64::MAX / 2, 8], 16);
        buf.extend_from_slice(&[0u8; 16]);
        write_u32(&mut buf, 0);
        write_u32(&mut buf, 0);
        assert!(OxiFile::read_from_bytes(&buf).is_err());

        // The consistent encoding of a two-element dataset decodes.
        let mut buf = dataset_prefix(1, &[2], 16);
        buf.extend_from_slice(&[0u8; 16]);
        write_u32(&mut buf, 0);
        write_u32(&mut buf, 0);
        let file = OxiFile::read_from_bytes(&buf).expect("consistent dataset");
        assert_eq!(file.root.datasets[0].shape.dims, vec![2]);
        assert_eq!(file.root.datasets[0].to_f64_vec().unwrap(), vec![0.0, 0.0]);
    }

    #[test]
    fn oxifile_huge_counts_do_not_preallocate() {
        let huge = u32::MAX;
        // Root group announcing billions of attributes, datasets or subgroups.
        let mut attrs = root_prefix(huge);
        attrs.extend_from_slice(&[0u8; 8]);
        assert!(OxiFile::read_from_bytes(&attrs).is_err());

        let mut datasets = root_prefix(0);
        write_u32(&mut datasets, huge);
        datasets.extend_from_slice(&[0u8; 8]);
        assert!(OxiFile::read_from_bytes(&datasets).is_err());

        let mut subgroups = root_prefix(0);
        write_u32(&mut subgroups, 0); // datasets
        write_u32(&mut subgroups, huge);
        subgroups.extend_from_slice(&[0u8; 8]);
        assert!(OxiFile::read_from_bytes(&subgroups).is_err());

        // Dataset announcing billions of dimensions.
        let mut dims = root_prefix(0);
        write_u32(&mut dims, 1);
        write_string(&mut dims, "d");
        dims.push(1);
        write_u32(&mut dims, huge);
        assert!(OxiFile::read_from_bytes(&dims).is_err());

        // Attribute arrays announcing 2^64 - 1 elements.
        for tag in [3u8, 4u8] {
            let mut buf = root_prefix(1);
            write_string(&mut buf, "a");
            buf.push(tag);
            write_u64(&mut buf, u64::MAX);
            assert!(OxiFile::read_from_bytes(&buf).is_err());
        }
    }

    #[test]
    fn oxifile_truncated_at_every_offset_returns_err() {
        let bytes = rich_oxifile().write_to_bytes();
        let decoded = OxiFile::read_from_bytes(&bytes).expect("complete file decodes");
        assert_eq!(decoded.root.attributes.len(), 5);
        let inner = &decoded.root.subgroups[0];
        assert_eq!(inner.datasets[0].to_i32_vec().unwrap(), vec![1, 2, 3, 4]);
        assert_eq!(inner.subgroups[0].name, "leaf");
        for len in 0..bytes.len() {
            assert!(
                OxiFile::read_from_bytes(&bytes[..len]).is_err(),
                "prefix of {len} bytes must be rejected"
            );
        }
    }

    #[test]
    fn oxifile_corrupted_bytes_never_panic() {
        fn assert_shapes_consistent(group: &Group) {
            for ds in &group.datasets {
                let expected = ds
                    .shape
                    .checked_total_elements()
                    .and_then(|n| n.checked_mul(ds.dtype.size_bytes()));
                assert_eq!(expected, Some(ds.data.len()), "dataset '{}'", ds.name);
            }
            for sub in &group.subgroups {
                assert_shapes_consistent(sub);
            }
        }
        let bytes = rich_oxifile().write_to_bytes();
        let mut decoded_variants = 0usize;
        let mut rejected_variants = 0usize;
        for at in 0..bytes.len() {
            for fill in [0x00u8, 0x01, 0x7F, 0xFF] {
                let mut corrupt = bytes.clone();
                corrupt[at] = fill;
                match OxiFile::read_from_bytes(&corrupt) {
                    Ok(file) => {
                        decoded_variants += 1;
                        assert_shapes_consistent(&file.root);
                    }
                    Err(_) => rejected_variants += 1,
                }
            }
        }
        // Payload corruption (e.g. a float byte) still decodes, structural
        // corruption (counts, lengths, tags) is rejected.
        assert!(decoded_variants > 0);
        assert!(rejected_variants > 0);
        // Damage to the magic number is never accepted.
        for at in 0..MAGIC.len() {
            let mut corrupt = bytes.clone();
            corrupt[at] ^= 0xFF;
            assert!(
                OxiFile::read_from_bytes(&corrupt).is_err(),
                "magic byte {at} corrupted but the file was accepted"
            );
        }
    }

    #[test]
    fn read_helpers_reject_overflowing_offsets() {
        let data = [1u8, 2, 3, 4, 5, 6, 7, 8];
        for start in [usize::MAX, usize::MAX - 3, data.len() + 1] {
            let mut pos = start;
            assert!(read_u32(&data, &mut pos).is_err());
            assert!(read_u64(&data, &mut pos).is_err());
            assert!(read_string(&data, &mut pos).is_err());
            assert_eq!(pos, start, "failed reads must not advance the cursor");
        }
        // A string whose length prefix points far beyond the buffer.
        let mut buf = Vec::new();
        write_u32(&mut buf, u32::MAX);
        let mut pos = 0usize;
        let err = read_string(&buf, &mut pos).expect_err("oversized string");
        assert!(err.contains("out of bounds"), "unexpected error: {err}");
    }

    #[test]
    fn dataset_accessors_handle_overflowing_shapes() {
        let overflowing = Dataset {
            name: "bad".to_string(),
            dtype: DataType::Float64,
            shape: DatasetShape::new(vec![usize::MAX, 2]),
            data: vec![0u8; 16],
            attributes: Vec::new(),
        };
        assert_eq!(overflowing.shape.checked_total_elements(), None);
        assert!(overflowing.to_f64_vec().is_err());
        let overflowing_i32 = Dataset {
            dtype: DataType::Int32,
            shape: DatasetShape::new(vec![usize::MAX]),
            ..overflowing
        };
        assert!(overflowing_i32.to_i32_vec().is_err());
        let scalar = Dataset::from_f64_slice("s", &[2.5], DatasetShape::new(vec![]));
        assert_eq!(scalar.shape.checked_total_elements(), Some(1));
        assert_eq!(scalar.to_f64_vec().unwrap(), vec![2.5]);
    }

    #[test]
    fn binary_mesh_and_particle_decoders_survive_huge_counts() {
        let mut mesh = Vec::new();
        write_u32(&mut mesh, u32::MAX);
        write_u32(&mut mesh, u32::MAX);
        let err = BinaryMesh::from_bytes(&mesh).expect_err("mesh counts exceed data");
        assert!(err.contains("truncated"), "unexpected error: {err}");

        let mut particles = Vec::new();
        particles.extend_from_slice(PARTICLE_MAGIC);
        write_u32(&mut particles, u32::MAX);
        write_u32(&mut particles, 0);
        let err = BinaryParticleData::from_bytes(&particles).expect_err("particle count");
        assert!(err.contains("truncated"), "unexpected error: {err}");

        let mut fields = Vec::new();
        fields.extend_from_slice(PARTICLE_MAGIC);
        write_u32(&mut fields, 0);
        write_u32(&mut fields, u32::MAX);
        assert!(BinaryParticleData::from_bytes(&fields).is_err());
    }

    #[test]
    fn oxifile_crate_writers_produce_files_the_decoder_accepts() {
        // Every writer of this module must emit a shape that matches its data,
        // because the decoder rejects datasets whose shape and byte length
        // disagree.
        let mut file = SimulationCheckpoint::create();
        let positions = [[0.5, 1.5, 2.5], [3.5, 4.5, 5.5], [-6.5, 7.5, -8.5]];
        let velocities = [[-1.0, 0.0, 1.0], [2.0, 3.0, 4.0], [5.0, 6.0, 7.0]];
        SimulationCheckpoint::add_positions(&mut file, "frame0", &positions);
        SimulationCheckpoint::add_velocities(&mut file, "frame0", &velocities);
        SimulationCheckpoint::add_scalar_field(&mut file, "frame0", "rho", &[1.0, 2.0, 3.0]);
        SimulationCheckpoint::add_timestep_metadata(&mut file, 7, 0.07, 0.01);
        // Containers without entries describe a consistent zero-element shape.
        SimulationCheckpoint::add_positions(&mut file, "empty", &[]);
        SimulationCheckpoint::add_scalar_field(&mut file, "empty", "rho", &[]);
        // A scalar dataset has an empty shape and exactly one element.
        file.root.add_dataset(Dataset::from_f64_slice(
            "dt",
            &[0.01],
            DatasetShape::new(vec![]),
        ));
        file.root.add_dataset(Dataset::from_f32_slice(
            "f32_matrix",
            &[1.5, 2.5, 3.5, 4.5, 5.5, 6.5],
            DatasetShape::new(vec![2, 3]),
        ));
        file.root.add_dataset(Dataset::from_i32_slice(
            "i32_vector",
            &[-1, 2, -3, 4],
            DatasetShape::new(vec![4]),
        ));
        let mut mesh = BinaryMesh::new();
        mesh.vertices = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.5, 1.0, 0.0]];
        mesh.triangles = vec![[0, 1, 2]];
        file.add_binary_mesh("geometry", &mesh);

        let bytes = file.write_to_bytes();
        let loaded = OxiFile::read_from_bytes(&bytes)
            .expect("the decoder must accept what the crate's own writers produce");
        assert_eq!(loaded.write_to_bytes(), bytes, "re-encoding must be stable");

        let frame = loaded.root.get_subgroup("frame0").expect("frame0");
        let flat =
            |vectors: &[[f64; 3]]| -> Vec<f64> { vectors.iter().flatten().copied().collect() };
        let stored = |name: &str| -> Vec<f64> {
            frame
                .get_dataset(name)
                .expect("dataset")
                .to_f64_vec()
                .expect("f64 dataset")
        };
        assert_eq!(stored("positions"), flat(&positions));
        assert_eq!(stored("velocities"), flat(&velocities));
        assert_eq!(stored("rho"), vec![1.0, 2.0, 3.0]);
        let positions_ds = frame.get_dataset("positions").expect("positions");
        assert_eq!(positions_ds.shape.dims, vec![3, 3]);

        let empty = loaded.root.get_subgroup("empty").expect("empty group");
        for name in ["positions", "rho"] {
            let ds = empty.get_dataset(name).expect("empty dataset");
            assert!(ds.to_f64_vec().expect("f64 dataset").is_empty());
        }
        let dt = loaded.root.get_dataset("dt").expect("dt");
        assert!(dt.shape.is_scalar());
        assert_eq!(dt.to_f64_vec().expect("f64 dataset"), vec![0.01]);
        let matrix = loaded.root.get_dataset("f32_matrix").expect("matrix");
        assert_eq!(matrix.shape.dims, vec![2, 3]);
        assert_eq!(
            matrix.to_f32_vec().expect("f32 dataset"),
            vec![1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
        );
        let vector = loaded.root.get_dataset("i32_vector").expect("vector");
        assert_eq!(
            vector.to_i32_vec().expect("i32 dataset"),
            vec![-1, 2, -3, 4]
        );

        let restored = loaded.get_binary_mesh("geometry").expect("mesh");
        assert_eq!(restored.vertices, mesh.vertices);
        assert_eq!(restored.triangles, mesh.triangles);
        let blob = loaded
            .root
            .get_subgroup("geometry")
            .and_then(|group| group.get_dataset("mesh_binary"))
            .expect("mesh blob");
        assert_eq!(blob.dtype, DataType::UInt8);
        assert_eq!(
            blob.get_attribute("n_vertices").map(|a| a.value.clone()),
            Some(AttributeValue::Int(3))
        );
        assert_eq!(
            loaded.root.get_attribute("step").map(|a| a.value.clone()),
            Some(AttributeValue::Int(7))
        );
    }

    #[test]
    fn oxifile_save_and_load_round_trip_through_a_file() {
        let dir = std::env::temp_dir();
        let pid = std::process::id();
        let whole = dir.join(format!("oxiphysics_oxifile_{pid}_whole.oxi"));
        let cut = dir.join(format!("oxiphysics_oxifile_{pid}_cut.oxi"));
        let original = rich_oxifile();
        let bytes = original.write_to_bytes();

        let saved = original.save(whole.to_str().expect("UTF-8 temp path"));
        let loaded = OxiFile::load(whole.to_str().expect("UTF-8 temp path"));
        std::fs::write(&cut, &bytes[..bytes.len() / 2]).expect("write truncated copy");
        let truncated = OxiFile::load(cut.to_str().expect("UTF-8 temp path"));
        std::fs::remove_file(&whole).ok();
        std::fs::remove_file(&cut).ok();

        saved.expect("save");
        let loaded = loaded.expect("load");
        assert_eq!(loaded.write_to_bytes(), bytes);
        assert!(
            truncated.is_err(),
            "half of a file must not decode as a whole one"
        );
        assert!(OxiFile::load(whole.to_str().expect("UTF-8 temp path")).is_err());
    }
}
