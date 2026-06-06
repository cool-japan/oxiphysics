// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! HDF5 file: root container, superblock, named types, path utilities.

use std::collections::HashMap;

use super::dataset::Hdf5Dataset;
use super::group::Hdf5Group;
use super::types::{AttrValue, Hdf5Dtype, Hdf5Error, Hdf5Result, LockState, ParallelHdf5Meta};

// ---------------------------------------------------------------------------
// Superblock & object-header simulation
// ---------------------------------------------------------------------------

/// Simulated HDF5 superblock (format metadata).
///
/// In real HDF5 the superblock occupies bytes 0-511 (v2) of the file and
/// stores the root-group offset, free-space information, etc.  Here we store
/// the metadata in a plain struct.
#[derive(Debug, Clone)]
pub struct Hdf5Superblock {
    /// Format version (0, 1, 2 or 3).
    pub version: u8,
    /// Total simulated file size in bytes.
    pub file_size: u64,
    /// Byte offset of the root group object header.
    pub root_obj_header_offset: u64,
    /// Byte offset of the end-of-file.
    pub eof_address: u64,
    /// Size of lengths in bytes (4 or 8 for 32/64-bit offsets).
    pub size_of_lengths: u8,
    /// Size of offsets in bytes.
    pub size_of_offsets: u8,
}

impl Default for Hdf5Superblock {
    fn default() -> Self {
        Self {
            version: 2,
            file_size: 0,
            root_obj_header_offset: 512,
            eof_address: 0,
            size_of_lengths: 8,
            size_of_offsets: 8,
        }
    }
}

/// Simulated HDF5 object header (per-object metadata).
#[derive(Debug, Clone)]
pub struct Hdf5ObjectHeader {
    /// Object type: "group" or "dataset".
    pub object_type: String,
    /// Simulated byte address of this header.
    pub address: u64,
    /// Number of messages in the header.
    pub n_messages: u32,
    /// Total header size in simulated bytes.
    pub header_size: u32,
}

// ---------------------------------------------------------------------------
// Named datatype registry
// ---------------------------------------------------------------------------

/// Registry for named (committed) datatypes.
#[derive(Debug, Clone, Default)]
pub struct NamedDatatypeRegistry {
    /// Registered types: name -> dtype.
    types: HashMap<String, Hdf5Dtype>,
}

impl NamedDatatypeRegistry {
    /// Register a named type.  Returns an error if the name is already taken.
    pub fn register(&mut self, name: &str, dtype: Hdf5Dtype) -> Hdf5Result<()> {
        if self.types.contains_key(name) {
            return Err(Hdf5Error::AlreadyExists(name.to_string()));
        }
        self.types.insert(name.to_string(), dtype);
        Ok(())
    }

    /// Look up a named type by name.
    pub fn get(&self, name: &str) -> Hdf5Result<&Hdf5Dtype> {
        self.types
            .get(name)
            .ok_or_else(|| Hdf5Error::NotFound(format!("named type '{name}'")))
    }

    /// List all registered type names.
    pub fn names(&self) -> Vec<String> {
        let mut v: Vec<String> = self.types.keys().cloned().collect();
        v.sort();
        v
    }
}

// ---------------------------------------------------------------------------
// HDF5 File (root container)
// ---------------------------------------------------------------------------

/// Top-level HDF5 mock file containing a group hierarchy, named-type registry
/// and file-level metadata.
#[derive(Debug, Clone)]
pub struct Hdf5File {
    /// Virtual filename.
    pub filename: String,
    /// Root group ("/").
    pub root: Hdf5Group,
    /// Simulated superblock.
    pub superblock: Hdf5Superblock,
    /// Named (committed) datatype registry.
    pub named_types: NamedDatatypeRegistry,
    /// Current file lock state.
    pub lock_state: LockState,
    /// Optional parallel write metadata.
    pub parallel_meta: Option<ParallelHdf5Meta>,
    /// Simulation clock tick (incremented on each write).
    pub write_tick: u64,
}

impl Hdf5File {
    /// Create a new empty in-memory HDF5 file.
    pub fn new(filename: &str) -> Self {
        Self {
            filename: filename.to_string(),
            root: Hdf5Group::new("/"),
            superblock: Hdf5Superblock::default(),
            named_types: NamedDatatypeRegistry::default(),
            lock_state: LockState::Unlocked,
            parallel_meta: None,
            write_tick: 0,
        }
    }

    // -- lock / unlock --

    /// Acquire an exclusive write lock (simulated).
    ///
    /// Returns `Err(Hdf5Error::FileLocked)` if the file is already locked.
    pub fn lock_write(&mut self, owner_id: u64) -> Hdf5Result<()> {
        match self.lock_state {
            LockState::Unlocked => {
                self.lock_state = LockState::WriteLocked { owner_id };
                Ok(())
            }
            _ => Err(Hdf5Error::FileLocked),
        }
    }

    /// Release the write lock.
    pub fn unlock(&mut self) {
        self.lock_state = LockState::Unlocked;
    }

    /// Acquire a shared read lock (multiple readers allowed).
    pub fn lock_read(&mut self) -> Hdf5Result<()> {
        match self.lock_state {
            LockState::Unlocked => {
                self.lock_state = LockState::ReadLocked { n_readers: 1 };
                Ok(())
            }
            LockState::ReadLocked { n_readers } => {
                self.lock_state = LockState::ReadLocked {
                    n_readers: n_readers + 1,
                };
                Ok(())
            }
            LockState::WriteLocked { .. } => Err(Hdf5Error::FileLocked),
        }
    }

    /// Return `true` if the file is write-locked.
    pub fn is_locked(&self) -> bool {
        matches!(self.lock_state, LockState::WriteLocked { .. })
    }

    // -- superblock helpers --

    /// Simulate updating the superblock's EOF address after a write.
    pub fn update_eof(&mut self, new_eof: u64) {
        self.superblock.eof_address = new_eof;
        self.superblock.file_size = new_eof;
    }

    // -- group traversal --

    /// Create a group at the given slash-separated path (relative to root).
    ///
    /// Intermediate groups are created as needed (like `mkdir -p`).
    pub fn create_group(&mut self, path: &str) -> Hdf5Result<()> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        let parts = split_path(path);
        let mut current = &mut self.root;
        for part in parts {
            if !current.groups.contains_key(part) {
                current
                    .groups
                    .insert(part.to_string(), Hdf5Group::new(part));
            }
            // SAFETY: We just ensured the key exists above.
            current = current
                .groups
                .get_mut(part)
                .unwrap_or_else(|| unreachable!());
        }
        Ok(())
    }

    /// Return a shared reference to the group at `path`.
    pub fn open_group(&self, path: &str) -> Hdf5Result<&Hdf5Group> {
        let parts = split_path(path);
        let mut current = &self.root;
        for part in parts {
            current = current
                .groups
                .get(part)
                .ok_or_else(|| Hdf5Error::NotFound(format!("group '{path}'")))?;
        }
        Ok(current)
    }

    /// Return a mutable reference to the group at `path`.
    pub fn open_group_mut(&mut self, path: &str) -> Hdf5Result<&mut Hdf5Group> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        let parts = split_path(path);
        let mut current = &mut self.root;
        for part in parts {
            current = current
                .groups
                .get_mut(part)
                .ok_or_else(|| Hdf5Error::NotFound(format!("group '{path}'")))?;
        }
        Ok(current)
    }

    // -- dataset access --

    /// Create a dataset at `group_path/dataset_name`.
    pub fn create_dataset(
        &mut self,
        group_path: &str,
        name: &str,
        shape: Vec<usize>,
        dtype: Hdf5Dtype,
    ) -> Hdf5Result<()> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        self.write_tick += 1;
        let group = self.open_group_mut(group_path)?;
        group.create_dataset(name, shape, dtype)
    }

    /// Return a shared reference to a dataset at `group_path/dataset_name`.
    pub fn open_dataset(&self, group_path: &str, name: &str) -> Hdf5Result<&Hdf5Dataset> {
        let group = self.open_group(group_path)?;
        group.open_dataset(name)
    }

    /// Return a mutable reference to a dataset.
    pub fn open_dataset_mut(
        &mut self,
        group_path: &str,
        name: &str,
    ) -> Hdf5Result<&mut Hdf5Dataset> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        let group = self.open_group_mut(group_path)?;
        group.open_dataset_mut(name)
    }

    // -- attribute helpers --

    /// Set an attribute on a dataset.
    pub fn set_dataset_attr(
        &mut self,
        group_path: &str,
        dataset: &str,
        attr_name: &str,
        value: AttrValue,
    ) -> Hdf5Result<()> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        let ds = self.open_dataset_mut(group_path, dataset)?;
        ds.set_attr(attr_name, value);
        Ok(())
    }

    /// Get an attribute from a dataset.
    pub fn get_dataset_attr(
        &self,
        group_path: &str,
        dataset: &str,
        attr_name: &str,
    ) -> Hdf5Result<&AttrValue> {
        let ds = self.open_dataset(group_path, dataset)?;
        ds.get_attr(attr_name)
    }

    // -- named types --

    /// Register a named (committed) datatype.
    pub fn commit_datatype(&mut self, name: &str, dtype: Hdf5Dtype) -> Hdf5Result<()> {
        self.named_types.register(name, dtype)
    }

    /// Look up a named datatype.
    pub fn find_named_type(&self, name: &str) -> Hdf5Result<&Hdf5Dtype> {
        self.named_types.get(name)
    }

    // -- links --

    /// Create a soft link inside `group_path`.
    pub fn create_soft_link(
        &mut self,
        group_path: &str,
        link_name: &str,
        target: &str,
    ) -> Hdf5Result<()> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        let group = self.open_group_mut(group_path)?;
        group.create_soft_link(link_name, target)
    }

    /// Create a hard link inside `group_path`.
    pub fn create_hard_link(
        &mut self,
        group_path: &str,
        link_name: &str,
        target: &str,
    ) -> Hdf5Result<()> {
        if self.is_locked() {
            return Err(Hdf5Error::FileLocked);
        }
        let group = self.open_group_mut(group_path)?;
        group.create_hard_link(link_name, target)
    }

    // -- parallel metadata --

    /// Attach parallel HDF5 metadata for an N-rank job.
    pub fn init_parallel(&mut self, n_ranks: usize) {
        self.parallel_meta = Some(ParallelHdf5Meta::new(n_ranks));
    }

    /// Record bytes written by MPI rank `rank`.
    pub fn record_rank_bytes(&mut self, rank: usize, bytes: u64) {
        if let Some(ref mut meta) = self.parallel_meta {
            meta.record_rank_bytes(rank, bytes);
        }
    }
}

// ---------------------------------------------------------------------------
// Path utilities
// ---------------------------------------------------------------------------

/// Split a slash-separated HDF5 path into its components, skipping empty parts.
pub(crate) fn split_path(path: &str) -> Vec<&str> {
    path.split('/').filter(|s| !s.is_empty()).collect()
}
