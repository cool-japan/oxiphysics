// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! High-level convenience functions, dimension scales, chunk iteration,
//! collective I/O, large-file support, checksums, string & group utilities,
//! matrix/tensor convenience, and annotations.

use super::dataset::Hdf5Dataset;
use super::file::Hdf5File;
use super::group::Hdf5Group;
use super::types::{
    AttrValue, CollectiveIoMeta, DataStorage, DimScale, Hdf5Dtype, Hdf5Error, Hdf5Result,
};

// ---------------------------------------------------------------------------
// High-level convenience functions
// ---------------------------------------------------------------------------

/// Write a 1-D slice of f64 data to the specified path in one call.
///
/// Creates intermediate groups and dataset as needed.
pub fn write_f64_dataset(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    data: &[f64],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    // Remove existing dataset if present (ignore error)
    let _ = file.create_dataset(group, name, vec![data.len()], Hdf5Dtype::Float64);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_f64(data)
}

/// Write a 1-D slice of f32 data.
pub fn write_f32_dataset(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    data: &[f32],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![data.len()], Hdf5Dtype::Float32);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_f32(data)
}

/// Write a 1-D slice of i32 data.
pub fn write_i32_dataset(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    data: &[i32],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![data.len()], Hdf5Dtype::Int32);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_i32(data)
}

/// Write a 1-D slice of u8 data.
pub fn write_u8_dataset(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    data: &[u8],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![data.len()], Hdf5Dtype::Uint8);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_u8(data)
}

/// Copy all datasets from a source group to a destination group in the same file.
pub fn copy_group_datasets(src_group: &Hdf5Group, dst_group: &mut Hdf5Group) -> Hdf5Result<()> {
    for (name, ds) in &src_group.datasets {
        if dst_group.datasets.contains_key(name) {
            return Err(Hdf5Error::AlreadyExists(name.clone()));
        }
        dst_group.datasets.insert(name.clone(), ds.clone());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Dimension-scale utilities
// ---------------------------------------------------------------------------

/// Create and register a 1-D dimension-scale dataset for axis coordinates.
pub fn create_dim_scale_1d(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    coords: &[f64],
    label: &str,
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![coords.len()], Hdf5Dtype::Float64);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_f64(coords)?;
    ds.make_dim_scale();
    ds.set_attr("CLASS", AttrValue::String("DIMENSION_SCALE".to_string()));
    ds.set_attr("NAME", AttrValue::String(name.to_string()));
    ds.set_attr("LABEL", AttrValue::String(label.to_string()));
    Ok(())
}

/// Attach a dimension scale to a dataset axis.
pub fn attach_dim_scale(
    file: &mut Hdf5File,
    group: &str,
    dataset: &str,
    scale_path: &str,
    axis: usize,
    label: &str,
) -> Hdf5Result<()> {
    let ds = file.open_dataset_mut(group, dataset)?;
    ds.attach_dim_scale(DimScale {
        scale_dataset: scale_path.to_string(),
        axis,
        label: label.to_string(),
    });
    Ok(())
}

// ---------------------------------------------------------------------------
// Chunk iterator
// ---------------------------------------------------------------------------

/// Iterator that yields the flat indices of all chunks in a dataset.
pub struct ChunkIterator {
    /// Dataset shape.
    shape: Vec<usize>,
    /// Chunk shape.
    chunk_shape: Vec<usize>,
    /// Current chunk multi-index.
    current: Vec<usize>,
    /// Whether the iteration is exhausted.
    done: bool,
}

impl ChunkIterator {
    /// Create an iterator over chunks of `shape` with step `chunk_shape`.
    pub fn new(shape: Vec<usize>, chunk_shape: Vec<usize>) -> Self {
        assert_eq!(shape.len(), chunk_shape.len());
        let ndim = shape.len();
        let done = shape.contains(&0);
        Self {
            shape,
            chunk_shape,
            current: vec![0; ndim],
            done,
        }
    }
}

impl Iterator for ChunkIterator {
    /// Each item is `(start_indices, actual_chunk_shape)`.
    type Item = (Vec<usize>, Vec<usize>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        let ndim = self.shape.len();
        let start = self.current.clone();
        // Compute actual chunk size (last chunk may be smaller).
        let actual: Vec<usize> = (0..ndim)
            .map(|d| (self.chunk_shape[d]).min(self.shape[d] - self.current[d]))
            .collect();

        // Advance to next chunk.
        let mut carry = true;
        for d in (0..ndim).rev() {
            if carry {
                self.current[d] += self.chunk_shape[d];
                if self.current[d] >= self.shape[d] {
                    self.current[d] = 0;
                } else {
                    carry = false;
                }
            }
        }
        if carry {
            self.done = true;
        }

        Some((start, actual))
    }
}

// ---------------------------------------------------------------------------
// Collective I/O simulation
// ---------------------------------------------------------------------------

/// Simulate a collective write operation across `n_ranks` ranks.
///
/// Each rank writes a contiguous slice of `data`; returns the metadata.
pub fn collective_write_f64(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    data: &[f64],
    n_ranks: usize,
) -> Hdf5Result<CollectiveIoMeta> {
    file.init_parallel(n_ranks);
    let chunk_size = data.len().div_ceil(n_ranks.max(1));
    for rank in 0..n_ranks {
        let bytes = (chunk_size.min(data.len().saturating_sub(rank * chunk_size)) * 8) as u64;
        file.record_rank_bytes(rank, bytes);
    }
    write_f64_dataset(file, group, name, data)?;
    let total_bytes = (data.len() * 8) as u64;
    Ok(CollectiveIoMeta {
        n_ranks,
        root_rank: 0,
        total_bytes,
        wall_time_s: total_bytes as f64 / (1024.0 * 1024.0 * 1024.0), // simulated at 1 GB/s
    })
}

// ---------------------------------------------------------------------------
// Large-file support helpers
// ---------------------------------------------------------------------------

/// Assign 64-bit byte offsets to all datasets in a file (simulation).
///
/// Traverses the group hierarchy and assigns sequential offsets,
/// simulating the behaviour of large-file-enabled HDF5.
pub fn assign_byte_offsets(file: &mut Hdf5File) {
    let mut offset: u64 = file.superblock.root_obj_header_offset + 512;
    assign_offsets_in_group(&mut file.root, &mut offset);
    file.update_eof(offset);
}

/// Recursive helper for `assign_byte_offsets`.
fn assign_offsets_in_group(group: &mut Hdf5Group, offset: &mut u64) {
    for ds in group.datasets.values_mut() {
        ds.byte_offset = *offset;
        let element_bytes = ds.dtype.element_size() as u64;
        let vol = ds.volume() as u64;
        *offset += element_bytes * vol + 64; // 64-byte header overhead
    }
    for child in group.groups.values_mut() {
        assign_offsets_in_group(child, offset);
    }
}

// ---------------------------------------------------------------------------
// MD5-like data integrity checksum (mock)
// ---------------------------------------------------------------------------

/// Compute a simple XOR-fold checksum over f64 data (not cryptographic).
///
/// Used to verify round-trip integrity in checkpoint tests.
pub fn data_checksum_f64(data: &[f64]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &v in data {
        let bits = v.to_bits();
        h ^= bits;
        h = h.wrapping_mul(0x00000100000001b3);
    }
    h
}

/// Verify that a dataset round-trips correctly through f64 write/read.
pub fn verify_roundtrip_f64(ds: &Hdf5Dataset, original: &[f64]) -> bool {
    if let Ok(v) = ds.read_f64() {
        v.len() == original.len() && v.iter().zip(original.iter()).all(|(a, b)| a == b)
    } else {
        false
    }
}

// ---------------------------------------------------------------------------
// String dataset utilities
// ---------------------------------------------------------------------------

/// Write a slice of variable-length strings to a dataset.
pub fn write_vlen_strings(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    strings: &[String],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![strings.len()], Hdf5Dtype::VlenString);
    let ds = file.open_dataset_mut(group, name)?;
    ds.data = DataStorage::VlenString(strings.to_vec());
    Ok(())
}

/// Read variable-length strings from a dataset.
pub fn read_vlen_strings(file: &Hdf5File, group: &str, name: &str) -> Hdf5Result<Vec<String>> {
    let ds = file.open_dataset(group, name)?;
    ds.read_vlen_strings()
}

// ---------------------------------------------------------------------------
// Additional group utilities
// ---------------------------------------------------------------------------

/// Count the total number of datasets in a group hierarchy (recursive).
pub fn count_datasets_recursive(group: &Hdf5Group) -> usize {
    let mut count = group.datasets.len();
    for child in group.groups.values() {
        count += count_datasets_recursive(child);
    }
    count
}

/// Collect all dataset paths in a group hierarchy (recursive).
pub fn list_datasets_recursive(group: &Hdf5Group, prefix: &str) -> Vec<String> {
    let mut paths = Vec::new();
    for name in group.datasets.keys() {
        paths.push(format!("{prefix}/{name}"));
    }
    for (name, child) in &group.groups {
        let child_prefix = format!("{prefix}/{name}");
        paths.extend(list_datasets_recursive(child, &child_prefix));
    }
    paths.sort();
    paths
}

// ---------------------------------------------------------------------------
// 2-D matrix convenience
// ---------------------------------------------------------------------------

/// Write a 2-D matrix (row-major) into the file.
pub fn write_matrix_f64(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    rows: usize,
    cols: usize,
    data: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(
        data.len(),
        rows * cols,
        "write_matrix_f64: data.len() != rows*cols"
    );
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![rows, cols], Hdf5Dtype::Float64);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_f64(data)
}

/// Read a 2-D matrix from the file, returning a `Vec<Vec`f64`>`.
pub fn read_matrix_f64(file: &Hdf5File, group: &str, name: &str) -> Hdf5Result<Vec<Vec<f64>>> {
    let ds = file.open_dataset(group, name)?;
    if ds.shape.len() != 2 {
        return Err(Hdf5Error::Generic(format!(
            "expected 2-D dataset, got {} dims",
            ds.shape.len()
        )));
    }
    let rows = ds.shape[0];
    let cols = ds.shape[1];
    let flat = ds.read_f64()?;
    let mut mat = Vec::with_capacity(rows);
    for r in 0..rows {
        mat.push(flat[r * cols..(r + 1) * cols].to_vec());
    }
    Ok(mat)
}

// ---------------------------------------------------------------------------
// 3-D tensor convenience
// ---------------------------------------------------------------------------

/// Write a 3-D tensor to the file.
pub fn write_tensor3_f64(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    d0: usize,
    d1: usize,
    d2: usize,
    data: &[f64],
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![d0, d1, d2], Hdf5Dtype::Float64);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_f64(data)
}

// ---------------------------------------------------------------------------
// Annotation helpers
// ---------------------------------------------------------------------------

/// Annotate a dataset with standard simulation metadata attributes.
pub fn annotate_dataset(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    units: &str,
    description: &str,
    n_atoms: usize,
    dt_ps: f64,
    creator: &str,
) -> Hdf5Result<()> {
    file.set_dataset_attr(group, name, "units", AttrValue::String(units.to_string()))?;
    file.set_dataset_attr(
        group,
        name,
        "description",
        AttrValue::String(description.to_string()),
    )?;
    file.set_dataset_attr(group, name, "n_atoms", AttrValue::Int32(n_atoms as i32))?;
    file.set_dataset_attr(group, name, "dt_ps", AttrValue::Float64(dt_ps))?;
    file.set_dataset_attr(
        group,
        name,
        "creator",
        AttrValue::String(creator.to_string()),
    )?;
    Ok(())
}
