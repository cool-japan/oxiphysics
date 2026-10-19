// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Conversion between [`CdfDataset`] and the `f64`-valued [`NetCdfFile`].

use super::{CdfAttribute, CdfDataset, CdfVersion, NcValues, read_cdf_dataset, write_cdf_dataset};
use crate::error::{Error, Result};
use crate::netcdf::types::{NetCdfFile, NetCdfVariable, VariableAttribute};

impl CdfDataset {
    /// Build a typed dataset from a [`NetCdfFile`].
    ///
    /// Each variable is stored with its [`NetCdfVariable::nc_type`]; string
    /// attributes become `NC_CHAR` attributes.  The dimension named by
    /// [`NetCdfFile::unlimited_dim`] becomes the record dimension, its size
    /// being the number of records.
    ///
    /// # Errors
    ///
    /// Fails when the unlimited dimension is not declared, a value cannot be
    /// represented in its variable's type (for example a fractional value in
    /// an integer variable), or the result violates
    /// [`CdfDataset::validate`].
    pub fn from_netcdf_file(file: &NetCdfFile) -> Result<Self> {
        let mut ds = CdfDataset::new();
        if let Some(unlimited) = &file.unlimited_dim
            && !file.dimensions.iter().any(|(n, _)| n == unlimited)
        {
            return Err(Error::General(format!(
                "unlimited dimension '{unlimited}' is not declared"
            )));
        }
        for (name, size) in &file.dimensions {
            if file.unlimited_dim.as_deref() == Some(name.as_str()) {
                ds.add_unlimited_dimension(name, *size)?;
            } else {
                ds.add_dimension(name, *size)?;
            }
        }
        for (key, value) in &file.attributes {
            ds.set_attribute(CdfAttribute::text(key, value));
        }
        for var in &file.variables {
            let values = NcValues::from_f64(var.nc_type, &var.data)
                .map_err(|e| Error::General(format!("variable '{}': {e}", var.name)))?;
            let dims: Vec<&str> = var.dims.iter().map(String::as_str).collect();
            let index = ds.add_variable(&var.name, &dims, values)?;
            for attr in &var.attributes {
                ds.variables[index].set_attribute(CdfAttribute::text(&attr.key, &attr.value));
            }
        }
        ds.validate()?;
        Ok(ds)
    }

    /// Convert to a [`NetCdfFile`].
    ///
    /// Variable values are widened to `f64` (exact except for 64-bit integers
    /// above 2^53) and keep their type in [`NetCdfVariable::nc_type`];
    /// `NC_CHAR` variables hold their byte codes.  Text attributes keep their
    /// text (without trailing NUL padding); numeric attributes are rendered
    /// as `ncdump` does (`", "`-separated values), so their type is not
    /// preserved — use the [`CdfDataset`] itself when it matters.
    pub fn to_netcdf_file(&self) -> NetCdfFile {
        let mut file = NetCdfFile::new();
        for dim in &self.dimensions {
            if dim.unlimited {
                file.add_unlimited_dimension(&dim.name, dim.len);
            } else {
                file.add_dimension(&dim.name, dim.len);
            }
        }
        for attr in &self.attributes {
            file.add_attribute(&attr.name, &attr.values.to_display_string());
        }
        for var in &self.variables {
            let dims = var
                .dim_ids
                .iter()
                .filter_map(|&id| self.dimensions.get(id).map(|d| d.name.clone()))
                .collect();
            let mut out = NetCdfVariable::new(&var.name, dims, var.values.to_f64_vec())
                .with_nc_type(var.nc_type());
            out.attributes = var
                .attributes
                .iter()
                .map(|a| VariableAttribute {
                    key: a.name.clone(),
                    value: a.values.to_display_string(),
                })
                .collect();
            file.add_variable(out);
        }
        file
    }

    /// Serialise in the classic binary format `version`
    /// (see [`write_cdf_dataset`]).
    ///
    /// # Errors
    ///
    /// As [`write_cdf_dataset`].
    pub fn to_bytes(&self, version: CdfVersion) -> Result<Vec<u8>> {
        write_cdf_dataset(self, version)
    }

    /// Parse classic binary bytes (see [`read_cdf_dataset`]).
    ///
    /// # Errors
    ///
    /// As [`read_cdf_dataset`].
    pub fn from_bytes(bytes: &[u8]) -> Result<Self> {
        read_cdf_dataset(bytes)
    }
}

/// Serialise a [`NetCdfFile`] in the NetCDF classic binary format.
///
/// The result starts with the `CDF\x01` / `CDF\x02` / `CDF\x05` magic and is
/// readable by netCDF-C, `ncdump`, scipy, MDAnalysis, cpptraj and VMD.
/// Each variable is written with its [`NetCdfVariable::nc_type`].
///
/// # Errors
///
/// See [`CdfDataset::from_netcdf_file`] and [`write_cdf_dataset`].
pub fn write_cdf(file: &NetCdfFile, version: CdfVersion) -> Result<Vec<u8>> {
    write_cdf_dataset(&CdfDataset::from_netcdf_file(file)?, version)
}

/// Parse NetCDF classic binary bytes into a [`NetCdfFile`]
/// (see [`CdfDataset::to_netcdf_file`] for the conversion rules).
///
/// # Errors
///
/// See [`read_cdf_dataset`].
pub fn read_cdf(bytes: &[u8]) -> Result<NetCdfFile> {
    Ok(read_cdf_dataset(bytes)?.to_netcdf_file())
}

/// Write a [`NetCdfFile`] to `path` in the NetCDF classic binary format.
///
/// # Errors
///
/// Encoding errors as [`write_cdf`], plus I/O errors.
pub fn write_cdf_file<P: AsRef<std::path::Path>>(
    path: P,
    file: &NetCdfFile,
    version: CdfVersion,
) -> Result<()> {
    let bytes = write_cdf(file, version)?;
    std::fs::write(path, bytes)?;
    Ok(())
}

/// Read a NetCDF classic binary file from `path`.
///
/// # Errors
///
/// I/O errors, and parse errors as [`read_cdf`].
pub fn read_cdf_file<P: AsRef<std::path::Path>>(path: P) -> Result<NetCdfFile> {
    let bytes = std::fs::read(path)?;
    read_cdf(&bytes)
}
