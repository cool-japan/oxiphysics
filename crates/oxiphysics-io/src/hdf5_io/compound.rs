// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Compound datatype helpers.

#![allow(dead_code)]

use std::collections::HashMap;

use super::dataset::Hdf5Dataset;
use super::types::{DataStorage, Hdf5Error, Hdf5Result};

/// A single record in a compound dataset (struct-like element).
#[derive(Debug, Clone)]
pub struct CompoundRecord {
    /// Field values keyed by field name.
    pub fields: HashMap<String, f64>,
}

impl CompoundRecord {
    /// Create a new record with no fields.
    pub fn new() -> Self {
        Self {
            fields: HashMap::new(),
        }
    }

    /// Set a field value.
    pub fn set(&mut self, field: &str, value: f64) {
        self.fields.insert(field.to_string(), value);
    }

    /// Get a field value.
    pub fn get(&self, field: &str) -> Hdf5Result<f64> {
        self.fields
            .get(field)
            .copied()
            .ok_or_else(|| Hdf5Error::NotFound(format!("field '{field}'")))
    }
}

impl Default for CompoundRecord {
    fn default() -> Self {
        Self::new()
    }
}

/// Write a slice of `CompoundRecord`s into a dataset.
pub fn write_compound_records(
    ds: &mut Hdf5Dataset,
    records: Vec<CompoundRecord>,
) -> Hdf5Result<()> {
    let maps: Vec<HashMap<String, f64>> = records.into_iter().map(|r| r.fields).collect();
    ds.data = DataStorage::Compound(maps);
    Ok(())
}

/// Read compound records back from a dataset.
pub fn read_compound_records(ds: &Hdf5Dataset) -> Hdf5Result<Vec<CompoundRecord>> {
    match &ds.data {
        DataStorage::Compound(maps) => Ok(maps
            .iter()
            .map(|m| CompoundRecord { fields: m.clone() })
            .collect()),
        _ => Err(Hdf5Error::Generic(
            "dataset is not compound type".to_string(),
        )),
    }
}
