// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HDF5 dataset: a named, typed, shaped array with optional metadata.

#![allow(dead_code)]

use std::collections::HashMap;

use super::types::{
    AttrValue, ChunkLayout, DataStorage, DimScale, ExternalRef, Hdf5Dtype, Hdf5Error, Hdf5Result,
    Hyperslab,
};

/// Compute row-major strides for a given shape.
pub(crate) fn compute_strides(shape: &[usize]) -> Vec<usize> {
    let ndim = shape.len();
    if ndim == 0 {
        return Vec::new();
    }
    let mut strides = vec![1_usize; ndim];
    for i in (0..ndim - 1).rev() {
        strides[i] = strides[i + 1] * shape[i + 1];
    }
    strides
}

/// An HDF5 dataset: a named, typed, shaped array with optional metadata.
#[derive(Debug, Clone)]
pub struct Hdf5Dataset {
    /// Dataset name (leaf).
    pub name: String,
    /// Shape (number of elements per dimension).
    pub shape: Vec<usize>,
    /// Element datatype.
    pub dtype: Hdf5Dtype,
    /// Flat element storage.
    pub data: DataStorage,
    /// Attributes attached to this dataset.
    pub attributes: HashMap<String, AttrValue>,
    /// Optional chunked storage descriptor.
    pub chunk_layout: Option<ChunkLayout>,
    /// Optional external dataset reference (virtual dataset).
    pub external_ref: Option<ExternalRef>,
    /// Dimension scale associations.
    pub dim_scales: Vec<DimScale>,
    /// Whether this dataset itself acts as a dimension scale.
    pub is_dim_scale: bool,
    /// 64-bit byte offset (simulated, large-file support).
    pub byte_offset: u64,
}

impl Hdf5Dataset {
    /// Return the total number of elements (`shape[0] * shape[1] * ...`).
    pub fn volume(&self) -> usize {
        if self.shape.is_empty() {
            0
        } else {
            self.shape.iter().product()
        }
    }

    /// Read the dataset as a flat `Vec`f64`.
    ///
    /// Returns an error if the storage variant does not contain float64 data.
    pub fn read_f64(&self) -> Hdf5Result<Vec<f64>> {
        match &self.data {
            DataStorage::Float64(v) => Ok(v.clone()),
            DataStorage::Float32(v) => Ok(v.iter().map(|&x| x as f64).collect()),
            _ => Err(Hdf5Error::Generic(format!(
                "dataset '{}' does not contain float64 data",
                self.name
            ))),
        }
    }

    /// Read the dataset as a flat `Vec`f32`.
    pub fn read_f32(&self) -> Hdf5Result<Vec<f32>> {
        match &self.data {
            DataStorage::Float32(v) => Ok(v.clone()),
            DataStorage::Float64(v) => Ok(v.iter().map(|&x| x as f32).collect()),
            _ => Err(Hdf5Error::Generic(format!(
                "dataset '{}' does not contain float32 data",
                self.name
            ))),
        }
    }

    /// Read the dataset as a flat `Vec`i32`.
    pub fn read_i32(&self) -> Hdf5Result<Vec<i32>> {
        match &self.data {
            DataStorage::Int32(v) => Ok(v.clone()),
            _ => Err(Hdf5Error::Generic(format!(
                "dataset '{}' does not contain int32 data",
                self.name
            ))),
        }
    }

    /// Read the dataset as a flat `Vec`u8`.
    pub fn read_u8(&self) -> Hdf5Result<Vec<u8>> {
        match &self.data {
            DataStorage::Uint8(v) => Ok(v.clone()),
            _ => Err(Hdf5Error::Generic(format!(
                "dataset '{}' does not contain uint8 data",
                self.name
            ))),
        }
    }

    /// Read variable-length strings from this dataset.
    pub fn read_vlen_strings(&self) -> Hdf5Result<Vec<String>> {
        match &self.data {
            DataStorage::VlenString(v) => Ok(v.clone()),
            _ => Err(Hdf5Error::Generic(format!(
                "dataset '{}' does not contain vlen-string data",
                self.name
            ))),
        }
    }

    /// Read a hyperslab from float64 data.
    ///
    /// Only 1-D and 2-D selections are fully supported; higher ranks are
    /// handled with a generic strided approach.
    pub fn read_hyperslab_f64(&self, slab: &Hyperslab) -> Hdf5Result<Vec<f64>> {
        slab.validate(&self.shape)?;
        let flat = self.read_f64()?;
        self.extract_hyperslab_f64(&flat, slab)
    }

    /// Internal helper: extract hyperslab from a flat buffer.
    fn extract_hyperslab_f64(&self, flat: &[f64], slab: &Hyperslab) -> Hdf5Result<Vec<f64>> {
        let ndim = self.shape.len();
        if ndim == 0 {
            return Ok(Vec::new());
        }
        let volume = slab.volume();
        let mut out = Vec::with_capacity(volume);

        // Compute strides for row-major layout.
        let strides = compute_strides(&self.shape);

        // Iterate over the selection using a flat counter.
        let mut indices = slab.start.clone();
        for _ in 0..volume {
            let flat_idx: usize = indices
                .iter()
                .zip(strides.iter())
                .map(|(&idx, &s)| idx * s)
                .sum();
            out.push(flat[flat_idx]);
            // Increment the multi-index from the last dimension.
            let mut carry = true;
            for d in (0..ndim).rev() {
                if carry {
                    indices[d] += 1;
                    if indices[d] >= slab.start[d] + slab.count[d] {
                        indices[d] = slab.start[d];
                    } else {
                        carry = false;
                    }
                }
            }
        }
        Ok(out)
    }

    /// Write float64 data replacing the entire dataset content.
    pub fn write_f64(&mut self, data: &[f64]) -> Hdf5Result<()> {
        let vol = self.volume();
        if data.len() != vol {
            return Err(Hdf5Error::ShapeMismatch {
                expected: self.shape.clone(),
                got: vec![data.len()],
            });
        }
        self.data = DataStorage::Float64(data.to_vec());
        Ok(())
    }

    /// Write float32 data replacing the entire dataset content.
    pub fn write_f32(&mut self, data: &[f32]) -> Hdf5Result<()> {
        let vol = self.volume();
        if data.len() != vol {
            return Err(Hdf5Error::ShapeMismatch {
                expected: self.shape.clone(),
                got: vec![data.len()],
            });
        }
        self.data = DataStorage::Float32(data.to_vec());
        Ok(())
    }

    /// Write i32 data replacing the entire dataset content.
    pub fn write_i32(&mut self, data: &[i32]) -> Hdf5Result<()> {
        let vol = self.volume();
        if data.len() != vol {
            return Err(Hdf5Error::ShapeMismatch {
                expected: self.shape.clone(),
                got: vec![data.len()],
            });
        }
        self.data = DataStorage::Int32(data.to_vec());
        Ok(())
    }

    /// Write u8 data replacing the entire dataset content.
    pub fn write_u8(&mut self, data: &[u8]) -> Hdf5Result<()> {
        let vol = self.volume();
        if data.len() != vol {
            return Err(Hdf5Error::ShapeMismatch {
                expected: self.shape.clone(),
                got: vec![data.len()],
            });
        }
        self.data = DataStorage::Uint8(data.to_vec());
        Ok(())
    }

    /// Set the value of a named attribute.
    pub fn set_attr(&mut self, name: &str, value: AttrValue) {
        self.attributes.insert(name.to_string(), value);
    }

    /// Get a reference to a named attribute.
    pub fn get_attr(&self, name: &str) -> Hdf5Result<&AttrValue> {
        self.attributes
            .get(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("attribute '{name}'")))
    }

    /// Attach a dimension scale to a specific axis.
    pub fn attach_dim_scale(&mut self, scale: DimScale) {
        self.dim_scales.push(scale);
    }

    /// Mark this dataset as a dimension scale with the given name.
    pub fn make_dim_scale(&mut self) {
        self.is_dim_scale = true;
    }

    /// Attach an external reference.
    pub fn set_external_ref(&mut self, ext: ExternalRef) {
        self.external_ref = Some(ext);
    }

    /// Set the chunk layout.
    pub fn set_chunk_layout(&mut self, layout: ChunkLayout) {
        self.chunk_layout = Some(layout);
    }

    /// List the names of all attributes on this dataset.
    pub fn attr_names(&self) -> Vec<String> {
        self.attributes.keys().cloned().collect()
    }
}
