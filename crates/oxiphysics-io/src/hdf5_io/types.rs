// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Core types for the HDF5 mock: error, dtype, storage, layout, attributes, links.

#![allow(dead_code)]

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// 1.  Error type
// ---------------------------------------------------------------------------

/// Errors that can occur when operating on the mock HDF5 store.
#[derive(Debug, Clone, PartialEq)]
pub enum Hdf5Error {
    /// A group, dataset or attribute was not found under the given path.
    NotFound(String),
    /// An item with the given name already exists.
    AlreadyExists(String),
    /// The requested hyperslab indices are out of range.
    HyperslabOutOfRange {
        /// The dimension index where the error occurred.
        dim: usize,
        /// The requested start offset.
        start: usize,
        /// The requested length.
        count: usize,
        /// The actual dimension size.
        size: usize,
    },
    /// The shapes or types are incompatible for the attempted operation.
    ShapeMismatch {
        /// Expected shape.
        expected: Vec<usize>,
        /// Provided shape.
        got: Vec<usize>,
    },
    /// The file is locked and cannot be written.
    FileLocked,
    /// A link target was not found.
    LinkTargetNotFound(String),
    /// Generic error carrying a human-readable message.
    Generic(String),
}

impl std::fmt::Display for Hdf5Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Hdf5Error::NotFound(p) => write!(f, "HDF5: not found: {p}"),
            Hdf5Error::AlreadyExists(p) => write!(f, "HDF5: already exists: {p}"),
            Hdf5Error::HyperslabOutOfRange {
                dim,
                start,
                count,
                size,
            } => write!(
                f,
                "HDF5: hyperslab out of range: dim={dim} start={start} count={count} size={size}"
            ),
            Hdf5Error::ShapeMismatch { expected, got } => {
                write!(
                    f,
                    "HDF5: shape mismatch: expected {expected:?}, got {got:?}"
                )
            }
            Hdf5Error::FileLocked => write!(f, "HDF5: file is locked"),
            Hdf5Error::LinkTargetNotFound(t) => {
                write!(f, "HDF5: link target not found: {t}")
            }
            Hdf5Error::Generic(msg) => write!(f, "HDF5: {msg}"),
        }
    }
}

/// Convenient alias for results from this module.
pub type Hdf5Result<T> = Result<T, Hdf5Error>;

// ---------------------------------------------------------------------------
// 2.  Datatypes
// ---------------------------------------------------------------------------

/// Scalar element type of an HDF5 dataset.
#[derive(Debug, Clone, PartialEq)]
pub enum Hdf5Dtype {
    /// 32-bit IEEE float.
    Float32,
    /// 64-bit IEEE float.
    Float64,
    /// 32-bit signed integer.
    Int32,
    /// 8-bit unsigned integer (byte).
    Uint8,
    /// Variable-length UTF-8 string.
    VlenString,
    /// Compound (struct-like) type described by a list of (name, dtype) fields.
    Compound(Vec<(String, Hdf5Dtype)>),
    /// A user-defined named type that aliases another dtype.
    Named {
        /// The user-supplied type name.
        name: String,
        /// The underlying base dtype.
        base: Box<Hdf5Dtype>,
    },
}

impl Hdf5Dtype {
    /// Return the byte-size of a single element (0 for variable-length types).
    pub fn element_size(&self) -> usize {
        match self {
            Hdf5Dtype::Float32 => 4,
            Hdf5Dtype::Float64 => 8,
            Hdf5Dtype::Int32 => 4,
            Hdf5Dtype::Uint8 => 1,
            Hdf5Dtype::VlenString => 0,
            Hdf5Dtype::Compound(fields) => fields.iter().map(|(_, dt)| dt.element_size()).sum(),
            Hdf5Dtype::Named { base, .. } => base.element_size(),
        }
    }
}

// ---------------------------------------------------------------------------
// 3.  Storage layout helpers
// ---------------------------------------------------------------------------

/// Chunked storage descriptor for a dataset.
#[derive(Debug, Clone)]
pub struct ChunkLayout {
    /// Chunk dimensions (same rank as dataset shape).
    pub chunk_shape: Vec<usize>,
    /// Gzip compression level (0 = no compression, 1-9).
    pub gzip_level: u8,
}

impl ChunkLayout {
    /// Create a new chunk descriptor.
    pub fn new(chunk_shape: Vec<usize>, gzip_level: u8) -> Self {
        Self {
            chunk_shape,
            gzip_level: gzip_level.min(9),
        }
    }

    /// Return the number of elements per chunk.
    pub fn chunk_volume(&self) -> usize {
        self.chunk_shape.iter().product()
    }
}

/// Hyperslab selection: a start offset + count (length) per dimension.
#[derive(Debug, Clone)]
pub struct Hyperslab {
    /// Start index per dimension.
    pub start: Vec<usize>,
    /// Number of elements to select per dimension.
    pub count: Vec<usize>,
}

impl Hyperslab {
    /// Create a new hyperslab.
    pub fn new(start: Vec<usize>, count: Vec<usize>) -> Self {
        assert_eq!(
            start.len(),
            count.len(),
            "Hyperslab: start and count must have equal rank"
        );
        Self { start, count }
    }

    /// Total number of elements in the selection.
    pub fn volume(&self) -> usize {
        self.count.iter().product()
    }

    /// Validate this selection against a dataset shape.
    pub fn validate(&self, shape: &[usize]) -> Hdf5Result<()> {
        if self.start.len() != shape.len() {
            return Err(Hdf5Error::Generic(format!(
                "Hyperslab rank {} != dataset rank {}",
                self.start.len(),
                shape.len()
            )));
        }
        for (dim, (&s, (&c, &sz))) in self
            .start
            .iter()
            .zip(self.count.iter().zip(shape.iter()))
            .enumerate()
        {
            if s + c > sz {
                return Err(Hdf5Error::HyperslabOutOfRange {
                    dim,
                    start: s,
                    count: c,
                    size: sz,
                });
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// 4.  Attribute value
// ---------------------------------------------------------------------------

/// The value stored in an HDF5 attribute.
#[derive(Debug, Clone, PartialEq)]
pub enum AttrValue {
    /// Scalar 64-bit float.
    Float64(f64),
    /// Scalar 32-bit float.
    Float32(f32),
    /// Scalar 32-bit integer.
    Int32(i32),
    /// UTF-8 string.
    String(String),
    /// Array of 64-bit floats.
    ArrayF64(Vec<f64>),
    /// Array of 32-bit floats.
    ArrayF32(Vec<f32>),
    /// Array of 32-bit integers.
    ArrayI32(Vec<i32>),
}

// ---------------------------------------------------------------------------
// 5.  Dataset storage
// ---------------------------------------------------------------------------

/// In-memory storage for a single dataset element type.
#[derive(Debug, Clone)]
pub enum DataStorage {
    /// Float32 array.
    Float32(Vec<f32>),
    /// Float64 array.
    Float64(Vec<f64>),
    /// Int32 array.
    Int32(Vec<i32>),
    /// Uint8 array.
    Uint8(Vec<u8>),
    /// Variable-length strings (one per element).
    VlenString(Vec<String>),
    /// Compound: each element is a map from field-name to float64 value.
    Compound(Vec<HashMap<String, f64>>),
}

impl DataStorage {
    /// Return the number of elements stored.
    pub fn len(&self) -> usize {
        match self {
            DataStorage::Float32(v) => v.len(),
            DataStorage::Float64(v) => v.len(),
            DataStorage::Int32(v) => v.len(),
            DataStorage::Uint8(v) => v.len(),
            DataStorage::VlenString(v) => v.len(),
            DataStorage::Compound(v) => v.len(),
        }
    }

    /// Return `true` if the storage is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// ---------------------------------------------------------------------------
// 6.  Link types
// ---------------------------------------------------------------------------

/// An HDF5 link (soft or hard) stored in a group.
#[derive(Debug, Clone)]
pub enum Hdf5Link {
    /// Soft (symbolic) link: stores the target path as a string.
    Soft(String),
    /// Hard link: stores the target path (both point to same object).
    Hard(String),
}

// ---------------------------------------------------------------------------
// 7.  External dataset reference
// ---------------------------------------------------------------------------

/// A reference to a dataset living in an external file.
#[derive(Debug, Clone)]
pub struct ExternalRef {
    /// Simulated external filename.
    pub filename: String,
    /// Path inside the external file.
    pub dataset_path: String,
    /// Byte offset (simulated, 64-bit for large-file support).
    pub byte_offset: u64,
}

// ---------------------------------------------------------------------------
// 8.  Dimension scale descriptor
// ---------------------------------------------------------------------------

/// A dimension-scale association following the HDF5 DimensionScales convention.
#[derive(Debug, Clone)]
pub struct DimScale {
    /// Path of the dataset that serves as a scale.
    pub scale_dataset: String,
    /// Axis index this scale is attached to.
    pub axis: usize,
    /// Human-readable label for this axis.
    pub label: String,
}

// ---------------------------------------------------------------------------
// 9.  Collective I/O metadata
// ---------------------------------------------------------------------------

/// Mock metadata for a collective-I/O operation (Parallel HDF5 style).
#[derive(Debug, Clone)]
pub struct CollectiveIoMeta {
    /// Number of MPI ranks participating.
    pub n_ranks: usize,
    /// Rank that initiated the write.
    pub root_rank: usize,
    /// Total bytes transferred across all ranks.
    pub total_bytes: u64,
    /// Simulated wall-clock time in seconds.
    pub wall_time_s: f64,
}

// ---------------------------------------------------------------------------
// File lock simulation
// ---------------------------------------------------------------------------

/// File-level lock state for exclusive write access simulation.
#[derive(Debug, Clone, PartialEq)]
pub enum LockState {
    /// No lock held.
    Unlocked,
    /// An exclusive write lock is held by `owner_id`.
    WriteLocked {
        /// Simulated process/thread identifier.
        owner_id: u64,
    },
    /// A shared read lock is held by `n_readers`.
    ReadLocked {
        /// Number of concurrent readers.
        n_readers: usize,
    },
}

// ---------------------------------------------------------------------------
// Parallel HDF5 metadata
// ---------------------------------------------------------------------------

/// Simulated parallel HDF5 (PHDF5) metadata for a multi-rank write.
#[derive(Debug, Clone)]
pub struct ParallelHdf5Meta {
    /// Total number of MPI ranks.
    pub n_ranks: usize,
    /// Per-rank byte counts.
    pub rank_byte_counts: Vec<u64>,
    /// Whether collective metadata writes are enabled.
    pub collective_metadata: bool,
    /// POSIX advisory lock held during the collective operation.
    pub lock_held: bool,
}

impl ParallelHdf5Meta {
    /// Create metadata for a `n_ranks`-rank job.
    pub fn new(n_ranks: usize) -> Self {
        Self {
            n_ranks,
            rank_byte_counts: vec![0; n_ranks],
            collective_metadata: true,
            lock_held: false,
        }
    }

    /// Record the number of bytes written by rank `rank`.
    pub fn record_rank_bytes(&mut self, rank: usize, bytes: u64) {
        if rank < self.n_ranks {
            self.rank_byte_counts[rank] = bytes;
        }
    }

    /// Total bytes across all ranks.
    pub fn total_bytes(&self) -> u64 {
        self.rank_byte_counts.iter().sum()
    }
}
