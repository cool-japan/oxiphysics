// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HDF5 group: a container for sub-groups, datasets, links and attributes.

#![allow(dead_code)]

use std::collections::HashMap;

use super::dataset::Hdf5Dataset;
use super::types::{AttrValue, DataStorage, Hdf5Dtype, Hdf5Error, Hdf5Link, Hdf5Result};

/// An HDF5 group that contains sub-groups, datasets, links and attributes.
#[derive(Debug, Clone)]
pub struct Hdf5Group {
    /// Group name (leaf component).
    pub name: String,
    /// Child groups keyed by name.
    pub groups: HashMap<String, Hdf5Group>,
    /// Datasets keyed by name.
    pub datasets: HashMap<String, Hdf5Dataset>,
    /// Links (soft or hard) keyed by link name.
    pub links: HashMap<String, Hdf5Link>,
    /// Group-level attributes.
    pub attributes: HashMap<String, AttrValue>,
}

impl Hdf5Group {
    /// Create an empty group with the given name.
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            groups: HashMap::new(),
            datasets: HashMap::new(),
            links: HashMap::new(),
            attributes: HashMap::new(),
        }
    }

    /// Create a child group.  Returns an error if it already exists.
    pub fn create_group(&mut self, name: &str) -> Hdf5Result<()> {
        if self.groups.contains_key(name) {
            return Err(Hdf5Error::AlreadyExists(name.to_string()));
        }
        self.groups.insert(name.to_string(), Hdf5Group::new(name));
        Ok(())
    }

    /// Return a mutable reference to a child group.
    pub fn open_group_mut(&mut self, name: &str) -> Hdf5Result<&mut Hdf5Group> {
        self.groups
            .get_mut(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("group '{name}'")))
    }

    /// Return a shared reference to a child group.
    pub fn open_group(&self, name: &str) -> Hdf5Result<&Hdf5Group> {
        self.groups
            .get(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("group '{name}'")))
    }

    /// Create and register a dataset with the given shape and dtype.
    ///
    /// The data storage is initialised to zeroed values.
    pub fn create_dataset(
        &mut self,
        name: &str,
        shape: Vec<usize>,
        dtype: Hdf5Dtype,
    ) -> Hdf5Result<()> {
        if self.datasets.contains_key(name) {
            return Err(Hdf5Error::AlreadyExists(name.to_string()));
        }
        let volume: usize = if shape.is_empty() {
            0
        } else {
            shape.iter().product()
        };
        let data = match dtype {
            Hdf5Dtype::Float32 => DataStorage::Float32(vec![0.0_f32; volume]),
            Hdf5Dtype::Float64 => DataStorage::Float64(vec![0.0_f64; volume]),
            Hdf5Dtype::Int32 => DataStorage::Int32(vec![0_i32; volume]),
            Hdf5Dtype::Uint8 => DataStorage::Uint8(vec![0_u8; volume]),
            Hdf5Dtype::VlenString => DataStorage::VlenString(vec![String::new(); volume]),
            Hdf5Dtype::Compound(_) => DataStorage::Compound(vec![HashMap::new(); volume]),
            Hdf5Dtype::Named { ref base, .. } => match base.as_ref() {
                Hdf5Dtype::Float64 => DataStorage::Float64(vec![0.0_f64; volume]),
                Hdf5Dtype::Float32 => DataStorage::Float32(vec![0.0_f32; volume]),
                Hdf5Dtype::Int32 => DataStorage::Int32(vec![0_i32; volume]),
                Hdf5Dtype::Uint8 => DataStorage::Uint8(vec![0_u8; volume]),
                _ => DataStorage::Float64(vec![0.0_f64; volume]),
            },
        };
        let ds = Hdf5Dataset {
            name: name.to_string(),
            shape,
            dtype,
            data,
            attributes: HashMap::new(),
            chunk_layout: None,
            external_ref: None,
            dim_scales: Vec::new(),
            is_dim_scale: false,
            byte_offset: 0,
        };
        self.datasets.insert(name.to_string(), ds);
        Ok(())
    }

    /// Return a mutable reference to a dataset.
    pub fn open_dataset_mut(&mut self, name: &str) -> Hdf5Result<&mut Hdf5Dataset> {
        self.datasets
            .get_mut(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("dataset '{name}'")))
    }

    /// Return a shared reference to a dataset.
    pub fn open_dataset(&self, name: &str) -> Hdf5Result<&Hdf5Dataset> {
        self.datasets
            .get(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("dataset '{name}'")))
    }

    /// Create a soft link.
    pub fn create_soft_link(&mut self, link_name: &str, target: &str) -> Hdf5Result<()> {
        if self.links.contains_key(link_name) {
            return Err(Hdf5Error::AlreadyExists(link_name.to_string()));
        }
        self.links
            .insert(link_name.to_string(), Hdf5Link::Soft(target.to_string()));
        Ok(())
    }

    /// Create a hard link.
    pub fn create_hard_link(&mut self, link_name: &str, target: &str) -> Hdf5Result<()> {
        if self.links.contains_key(link_name) {
            return Err(Hdf5Error::AlreadyExists(link_name.to_string()));
        }
        self.links
            .insert(link_name.to_string(), Hdf5Link::Hard(target.to_string()));
        Ok(())
    }

    /// Set a group-level attribute.
    pub fn set_attr(&mut self, name: &str, value: AttrValue) {
        self.attributes.insert(name.to_string(), value);
    }

    /// Get a shared reference to a group-level attribute.
    pub fn get_attr(&self, name: &str) -> Hdf5Result<&AttrValue> {
        self.attributes
            .get(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("attribute '{name}'")))
    }

    /// List the names of all child groups, datasets and links.
    pub fn list_contents(&self) -> Vec<String> {
        let mut names: Vec<String> = self
            .groups
            .keys()
            .chain(self.datasets.keys())
            .chain(self.links.keys())
            .cloned()
            .collect();
        names.sort();
        names
    }

    /// Total number of datasets (excluding sub-groups).
    pub fn n_datasets(&self) -> usize {
        self.datasets.len()
    }

    /// Total number of child groups.
    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }
}
