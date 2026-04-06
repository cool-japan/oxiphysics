// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HDF5 file-image serialisation, region references, file statistics,
//! ring-buffer trajectories, merge helpers, snapshot utilities, and
//! extendable datasets.

#![allow(dead_code)]

use super::convenience::{
    list_datasets_recursive, read_vlen_strings, write_f64_dataset, write_matrix_f64,
    write_vlen_strings,
};
use super::file::Hdf5File;
use super::group::Hdf5Group;
use super::types::{AttrValue, Hdf5Dtype, Hdf5Error, Hdf5Result, Hyperslab};

// ===========================================================================
// HDF5 file-image serialisation (binary round-trip simulation)
// ===========================================================================

/// A byte-oriented image of an in-memory HDF5 file.
///
/// Encodes groups, datasets and attributes into a compact binary
/// representation for testing checkpoint / restart round-trip workflows.
#[derive(Debug, Clone)]
pub struct Hdf5FileImage {
    /// Version tag.
    pub version: u8,
    /// Serialised group segments: `(path, attr_json)`.
    pub group_segments: Vec<(String, String)>,
    /// Serialised dataset segments: `(group_path, name, shape, flat_f64)`.
    pub dataset_segments: Vec<(String, String, Vec<usize>, Vec<f64>)>,
    /// Serialised attribute segments: `(group, name, attr_name, value_str)`.
    pub attr_segments: Vec<(String, String, String, String)>,
}

impl Hdf5FileImage {
    /// Build an image from an `Hdf5File`.
    ///
    /// Only datasets containing f64 data are captured; other dtypes are skipped.
    pub fn from_file(file: &Hdf5File) -> Self {
        let mut img = Self {
            version: 1,
            group_segments: Vec::new(),
            dataset_segments: Vec::new(),
            attr_segments: Vec::new(),
        };
        Self::capture_group(file, &file.root, "", &mut img);
        img
    }

    fn capture_group(_file: &Hdf5File, group: &Hdf5Group, prefix: &str, img: &mut Self) {
        let path = if prefix.is_empty() {
            "/".to_string()
        } else {
            prefix.to_string()
        };
        let attrs_str = format!("{:?}", group.attributes.keys().collect::<Vec<_>>());
        img.group_segments.push((path.clone(), attrs_str));

        for (ds_name, ds) in &group.datasets {
            if let Ok(flat) = ds.read_f64() {
                img.dataset_segments
                    .push((path.clone(), ds_name.clone(), ds.shape.clone(), flat));
            }
            for (attr_name, attr_val) in &ds.attributes {
                let val_str = format!("{attr_val:?}");
                img.attr_segments
                    .push((path.clone(), ds_name.clone(), attr_name.clone(), val_str));
            }
        }
        for (child_name, child) in &group.groups {
            let child_prefix = if prefix.is_empty() {
                child_name.clone()
            } else {
                format!("{prefix}/{child_name}")
            };
            Self::capture_group(_file, child, &child_prefix, img);
        }
    }

    /// Restore datasets from this image into a (possibly empty) `Hdf5File`.
    pub fn restore_to_file(&self, file: &mut Hdf5File) -> Hdf5Result<()> {
        for (grp, name, shape, data) in &self.dataset_segments {
            let grp_path = if grp == "/" { "" } else { grp.as_str() };
            // Ensure the group exists.
            if !grp_path.is_empty() {
                file.create_group(grp_path)?;
            }
            let _ = file.create_dataset(grp_path, name, shape.clone(), Hdf5Dtype::Float64);
            let ds = file.open_dataset_mut(grp_path, name)?;
            ds.write_f64(data)?;
        }
        Ok(())
    }

    /// Total number of dataset segments in the image.
    pub fn n_datasets(&self) -> usize {
        self.dataset_segments.len()
    }
}

// ===========================================================================
// HDF5 region reference
// ===========================================================================

/// A region reference: points to a hyperslab inside a specific dataset.
#[derive(Debug, Clone)]
pub struct RegionReference {
    /// Group path containing the target dataset.
    pub group: String,
    /// Dataset name.
    pub dataset: String,
    /// Hyperslab selection.
    pub slab: Hyperslab,
}

impl RegionReference {
    /// Create a region reference.
    pub fn new(group: &str, dataset: &str, slab: Hyperslab) -> Self {
        Self {
            group: group.to_string(),
            dataset: dataset.to_string(),
            slab,
        }
    }

    /// Dereference: read the selected data from the file.
    pub fn dereference_f64(&self, file: &Hdf5File) -> Hdf5Result<Vec<f64>> {
        let ds = file.open_dataset(&self.group, &self.dataset)?;
        ds.read_hyperslab_f64(&self.slab)
    }
}

// ===========================================================================
// 3-D grid dataset helper
// ===========================================================================

/// Write a 3-D scalar field (e.g. electron density) to a dataset.
///
/// `data` has shape `[nx x ny x nz]` stored flat in row-major order.
pub fn write_3d_grid_f64(
    file: &mut Hdf5File,
    group: &str,
    name: &str,
    nx: usize,
    ny: usize,
    nz: usize,
    data: &[f64],
) -> Hdf5Result<()> {
    assert_eq!(data.len(), nx * ny * nz, "3D grid size mismatch");
    file.create_group(group)?;
    let _ = file.create_dataset(group, name, vec![nx, ny, nz], Hdf5Dtype::Float64);
    let ds = file.open_dataset_mut(group, name)?;
    ds.write_f64(data)?;
    ds.set_attr("grid_type", AttrValue::String("scalar_3d".to_string()));
    Ok(())
}

/// Read a 3-D scalar field, returning `(nx, ny, nz, data)`.
pub fn read_3d_grid_f64(
    file: &Hdf5File,
    group: &str,
    name: &str,
) -> Hdf5Result<(usize, usize, usize, Vec<f64>)> {
    let ds = file.open_dataset(group, name)?;
    if ds.shape.len() != 3 {
        return Err(Hdf5Error::Generic(format!(
            "expected 3D dataset, got {}D",
            ds.shape.len()
        )));
    }
    let data = ds.read_f64()?;
    Ok((ds.shape[0], ds.shape[1], ds.shape[2], data))
}

// ===========================================================================
// Force file (MD output)
// ===========================================================================

/// Write atomic forces `[n_atoms x n_frames x 3]` to a group.
pub fn write_forces(
    file: &mut Hdf5File,
    group: &str,
    forces: &[f64],
    n_frames: usize,
    n_atoms: usize,
) -> Hdf5Result<()> {
    file.create_group(group)?;
    let _ = file.create_dataset(
        group,
        "forces",
        vec![n_frames, n_atoms, 3],
        Hdf5Dtype::Float64,
    );
    let ds = file.open_dataset_mut(group, "forces")?;
    ds.write_f64(forces)?;
    ds.set_attr("units", AttrValue::String("kJ/mol/nm".to_string()));
    Ok(())
}

// ===========================================================================
// Energy dataset helper
// ===========================================================================

/// Write a per-frame energy array to a group.
pub fn write_energies(file: &mut Hdf5File, group: &str, energies: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "potential_energy", energies)?;
    file.set_dataset_attr(
        group,
        "potential_energy",
        "units",
        AttrValue::String("kJ/mol".to_string()),
    )?;
    Ok(())
}

// ===========================================================================
// Atom types string dataset
// ===========================================================================

/// Write atom-type strings (e.g. `["H", "H", "O"]`) to a dataset.
pub fn write_atom_types(file: &mut Hdf5File, group: &str, atom_types: &[String]) -> Hdf5Result<()> {
    write_vlen_strings(file, group, "atom_types", atom_types)
}

/// Read atom-type strings from a dataset.
pub fn read_atom_types(file: &Hdf5File, group: &str) -> Hdf5Result<Vec<String>> {
    read_vlen_strings(file, group, "atom_types")
}

// ===========================================================================
// Pairwise distance matrix helper
// ===========================================================================

/// Compute and write a pairwise distance matrix for `n_atoms` positions.
pub fn write_distance_matrix(
    file: &mut Hdf5File,
    group: &str,
    positions: &[[f64; 3]],
) -> Hdf5Result<()> {
    let n = positions.len();
    let mut mat = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let dx = positions[i][0] - positions[j][0];
            let dy = positions[i][1] - positions[j][1];
            let dz = positions[i][2] - positions[j][2];
            mat[i * n + j] = (dx * dx + dy * dy + dz * dz).sqrt();
        }
    }
    write_matrix_f64(file, group, "distance_matrix", n, n, &mat)
}

// ===========================================================================
// HDF5 iterator: depth-first walk
// ===========================================================================

/// Collect all dataset paths under a file root using depth-first traversal.
pub fn walk_datasets(file: &Hdf5File) -> Vec<String> {
    list_datasets_recursive(&file.root, "")
}

/// Collect all group paths using depth-first traversal.
pub fn walk_groups(group: &Hdf5Group, prefix: &str) -> Vec<String> {
    let mut paths = Vec::new();
    for (name, child) in &group.groups {
        let p = format!("{prefix}/{name}");
        paths.push(p.clone());
        paths.extend(walk_groups(child, &p));
    }
    paths.sort();
    paths
}

// ===========================================================================
// Metadata header for physics files
// ===========================================================================

/// Standard metadata written to the root group of a physics output file.
#[derive(Debug, Clone)]
pub struct PhysicsFileHeader {
    /// Name of the simulation code.
    pub code_name: String,
    /// Code version string.
    pub code_version: String,
    /// Title of the simulation.
    pub title: String,
    /// ISO-8601 creation timestamp (simulated).
    pub created: String,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Integration time step in picoseconds.
    pub dt_ps: f64,
    /// Total simulation time in picoseconds.
    pub total_time_ps: f64,
}

impl PhysicsFileHeader {
    /// Write all metadata as group-level attributes on the root group.
    pub fn write_to_file(&self, file: &mut Hdf5File) -> Hdf5Result<()> {
        // Write metadata dataset under "metadata"
        file.create_group("metadata")?;
        let grp = file.open_group_mut("metadata")?;
        grp.set_attr("code_name", AttrValue::String(self.code_name.clone()));
        grp.set_attr("code_version", AttrValue::String(self.code_version.clone()));
        grp.set_attr("title", AttrValue::String(self.title.clone()));
        grp.set_attr("created", AttrValue::String(self.created.clone()));
        grp.set_attr("n_atoms", AttrValue::Int32(self.n_atoms as i32));
        grp.set_attr("dt_ps", AttrValue::Float64(self.dt_ps));
        grp.set_attr("total_time_ps", AttrValue::Float64(self.total_time_ps));
        Ok(())
    }

    /// Read the metadata back from the file.
    pub fn read_from_file(file: &Hdf5File) -> Hdf5Result<Self> {
        let grp = file.open_group("metadata")?;
        let code_name = match grp.get_attr("code_name")? {
            AttrValue::String(s) => s.clone(),
            _ => String::new(),
        };
        let n_atoms = match grp.get_attr("n_atoms")? {
            AttrValue::Int32(v) => *v as usize,
            _ => 0,
        };
        let dt_ps = match grp.get_attr("dt_ps")? {
            AttrValue::Float64(v) => *v,
            _ => 0.0,
        };
        let total_time_ps = match grp.get_attr("total_time_ps")? {
            AttrValue::Float64(v) => *v,
            _ => 0.0,
        };
        Ok(Self {
            code_name,
            code_version: String::new(),
            title: String::new(),
            created: String::new(),
            n_atoms,
            dt_ps,
            total_time_ps,
        })
    }
}

// ===========================================================================
// Ring-buffer trajectory store
// ===========================================================================

/// A ring-buffer that stores only the last N frames of a trajectory.
#[derive(Debug, Clone)]
pub struct RingTrajectory {
    /// Maximum number of frames.
    pub capacity: usize,
    /// Flat position storage: `[capacity * n_atoms * 3]`.
    pub storage: Vec<f64>,
    /// Number of atoms.
    pub n_atoms: usize,
    /// Write head (index of next slot to write).
    pub head: usize,
    /// Total frames ever appended (for determining full/empty state).
    pub total_appended: usize,
}

impl RingTrajectory {
    /// Create a new ring-buffer trajectory.
    pub fn new(capacity: usize, n_atoms: usize) -> Self {
        Self {
            capacity,
            storage: vec![0.0_f64; capacity * n_atoms * 3],
            n_atoms,
            head: 0,
            total_appended: 0,
        }
    }

    /// Append a frame (positions `[n_atoms * 3]`).
    pub fn append(&mut self, positions: &[f64]) {
        assert_eq!(positions.len(), self.n_atoms * 3);
        let base = self.head * self.n_atoms * 3;
        self.storage[base..base + self.n_atoms * 3].copy_from_slice(positions);
        self.head = (self.head + 1) % self.capacity;
        self.total_appended += 1;
    }

    /// Number of frames currently stored (min of capacity and total appended).
    pub fn n_stored(&self) -> usize {
        self.total_appended.min(self.capacity)
    }

    /// Read the frame at logical index `i` (0 = oldest, n_stored-1 = newest).
    pub fn read_frame(&self, i: usize) -> Hdf5Result<Vec<[f64; 3]>> {
        let n = self.n_stored();
        if i >= n {
            return Err(Hdf5Error::NotFound(format!("ring frame {i}")));
        }
        let oldest_slot = if self.total_appended < self.capacity {
            0
        } else {
            self.head
        };
        let slot = (oldest_slot + i) % self.capacity;
        let base = slot * self.n_atoms * 3;
        let out: Vec<[f64; 3]> = (0..self.n_atoms)
            .map(|a| {
                let p = base + a * 3;
                [self.storage[p], self.storage[p + 1], self.storage[p + 2]]
            })
            .collect();
        Ok(out)
    }
}

// ===========================================================================
// HDF5 file merge helper
// ===========================================================================

/// Merge all datasets from `src` into `dst`, skipping duplicates.
pub fn merge_files(src: &Hdf5File, dst: &mut Hdf5File) -> Hdf5Result<usize> {
    let paths = list_datasets_recursive(&src.root, "");
    let mut merged = 0;
    for path in paths {
        // path = "/group/name" or "/name"
        let parts: Vec<&str> = path.trim_start_matches('/').rsplitn(2, '/').collect();
        let (name, group) = if parts.len() == 2 {
            (parts[0], parts[1])
        } else {
            (parts[0], "")
        };
        // Locate source dataset
        let src_grp = if group.is_empty() {
            &src.root
        } else {
            match src.open_group(group) {
                Ok(g) => g,
                Err(_) => continue,
            }
        };
        let src_ds = match src_grp.open_dataset(name) {
            Ok(ds) => ds.clone(),
            Err(_) => continue,
        };
        // Create in destination (skip if already exists)
        if !group.is_empty() {
            dst.create_group(group)?;
        }
        let dst_grp = if group.is_empty() {
            &mut dst.root
        } else {
            match dst.open_group_mut(group) {
                Ok(g) => g,
                Err(_) => continue,
            }
        };
        if dst_grp.datasets.contains_key(name) {
            continue;
        }
        dst_grp.datasets.insert(name.to_string(), src_ds);
        merged += 1;
    }
    Ok(merged)
}

// ===========================================================================
// Atom position round-trip helper
// ===========================================================================

/// Write atom positions and types together, then read them back.
pub fn write_snapshot(
    file: &mut Hdf5File,
    group: &str,
    positions: &[[f64; 3]],
    atom_types: &[String],
) -> Hdf5Result<()> {
    let n = positions.len();
    assert_eq!(n, atom_types.len());
    let flat: Vec<f64> = positions.iter().flat_map(|p| p.iter().cloned()).collect();
    file.create_group(group)?;
    file.create_dataset(group, "positions", vec![n, 3], Hdf5Dtype::Float64)?;
    file.open_dataset_mut(group, "positions")?
        .write_f64(&flat)?;
    write_vlen_strings(file, group, "atom_types", atom_types)?;
    Ok(())
}

/// Read snapshot positions and types back.
pub fn read_snapshot(file: &Hdf5File, group: &str) -> Hdf5Result<(Vec<[f64; 3]>, Vec<String>)> {
    let ds = file.open_dataset(group, "positions")?;
    let flat = ds.read_f64()?;
    let n = ds.shape[0];
    let positions: Vec<[f64; 3]> = (0..n)
        .map(|i| [flat[i * 3], flat[i * 3 + 1], flat[i * 3 + 2]])
        .collect();
    let types = read_vlen_strings(file, group, "atom_types")?;
    Ok((positions, types))
}

// ===========================================================================
// B-factor / temperature factor dataset
// ===========================================================================

/// Write per-atom B-factors (crystallographic temperature factors) to the file.
pub fn write_bfactors(file: &mut Hdf5File, group: &str, bfactors: &[f64]) -> Hdf5Result<()> {
    write_f64_dataset(file, group, "bfactor", bfactors)?;
    file.set_dataset_attr(
        group,
        "bfactor",
        "units",
        AttrValue::String("Angstrom^2".to_string()),
    )?;
    Ok(())
}

// ===========================================================================
// HDF5 file statistics
// ===========================================================================

/// Compute basic statistics over all f64 datasets in a file.
#[derive(Debug, Clone)]
pub struct FileStats {
    /// Total number of datasets.
    pub n_datasets: usize,
    /// Total number of f64 elements across all datasets.
    pub total_elements: usize,
    /// Global minimum value.
    pub global_min: f64,
    /// Global maximum value.
    pub global_max: f64,
    /// Global mean value.
    pub global_mean: f64,
}

impl FileStats {
    /// Compute statistics from a file.
    pub fn compute(file: &Hdf5File) -> Self {
        let paths = list_datasets_recursive(&file.root, "");
        let mut all_data: Vec<f64> = Vec::new();
        let mut n_ds = 0;

        for path in &paths {
            let parts: Vec<&str> = path.trim_start_matches('/').rsplitn(2, '/').collect();
            let (name, group) = if parts.len() == 2 {
                (parts[0], parts[1])
            } else {
                (parts[0], "")
            };
            let src_grp: &Hdf5Group = if group.is_empty() {
                &file.root
            } else {
                match file.open_group(group) {
                    Ok(g) => g,
                    Err(_) => continue,
                }
            };
            if let Ok(ds) = src_grp.open_dataset(name)
                && let Ok(data) = ds.read_f64()
            {
                all_data.extend_from_slice(&data);
                n_ds += 1;
            }
        }

        if all_data.is_empty() {
            return Self {
                n_datasets: n_ds,
                total_elements: 0,
                global_min: 0.0,
                global_max: 0.0,
                global_mean: 0.0,
            };
        }
        let min = all_data.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = all_data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mean = all_data.iter().sum::<f64>() / all_data.len() as f64;

        Self {
            n_datasets: n_ds,
            total_elements: all_data.len(),
            global_min: min,
            global_max: max,
            global_mean: mean,
        }
    }
}

// ===========================================================================
// Incremental dataset append (simulated extendable datasets)
// ===========================================================================

/// An incrementally-extendable 1-D dataset.
///
/// Simulates HDF5 unlimited-dimension datasets by re-creating the backing
/// vector on each extension.
#[derive(Debug, Clone)]
pub struct ExtendableDataset {
    /// Dataset name.
    pub name: String,
    /// Current data values.
    pub data: Vec<f64>,
    /// Chunk size for the underlying buffer.
    pub chunk_size: usize,
}

impl ExtendableDataset {
    /// Create a new extendable dataset.
    pub fn new(name: &str, chunk_size: usize) -> Self {
        Self {
            name: name.to_string(),
            data: Vec::new(),
            chunk_size: chunk_size.max(1),
        }
    }

    /// Append values to the dataset.
    pub fn extend(&mut self, values: &[f64]) {
        self.data.extend_from_slice(values);
    }

    /// Flush the current data to a group in an HDF5 file.
    pub fn flush(&self, file: &mut Hdf5File, group: &str) -> Hdf5Result<()> {
        write_f64_dataset(file, group, &self.name, &self.data)
    }

    /// Current length.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Return true if empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }
}
