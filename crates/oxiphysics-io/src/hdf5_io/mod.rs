// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Pure-Rust in-memory HDF5 mock for simulation I/O.
//!
//! This module simulates the HDF5 file-format API entirely in memory using
//! Rust data structures.  No C library (libhdf5) is required.  The API
//! surface mirrors the logical HDF5 model:
//!
//! * [`Hdf5File`] — the root container (analogous to an open HDF5 file).
//! * [`Hdf5Group`] — a named group that can hold sub-groups, datasets and
//!   attributes.
//! * [`Hdf5Dataset`] — a multi-dimensional typed array.
//! * `Hdf5Attribute` — scalar, string or array metadata.
//!
//! The implementation additionally covers:
//! - Chunked storage simulation
//! - Gzip-level compression placeholder
//! - Hyperslab (slice) selection
//! - Named datatypes
//! - Soft and hard links
//! - External dataset references
//! - Dimension scales (NetCDF-like convention)
//! - Multi-file virtual dataset simulation
//! - File lock/unlock simulation
//! - Collective I/O metadata
//! - Compound (struct-like) datatypes
//! - Variable-length strings
//! - Trajectory storage (atoms x time x 3)
//! - Checkpoint / restart round-trip
//! - Large-file offset simulation (64-bit)
//! - Parallel HDF5 metadata (mock)

#![allow(dead_code)]

// --- Sub-modules ---
pub mod compound;
pub mod convenience;
pub mod dataset;
pub mod encoding;
pub mod file;
pub mod file_image;
pub mod group;
pub mod science_io;
pub mod science_io_extra;
pub mod science_io_md;
pub mod science_io_micro;
pub mod trajectory;
pub mod types;

#[cfg(test)]
mod tests;

// --- Re-exports: preserve the original flat public API ---
pub use compound::*;
pub use convenience::*;
pub use dataset::Hdf5Dataset;
#[cfg(test)]
pub(crate) use dataset::compute_strides;
pub use encoding::*;
pub use file::*;
pub use file_image::*;
pub use group::Hdf5Group;
pub use science_io::*;
pub use science_io_extra::*;
pub use science_io_md::*;
pub use science_io_micro::*;
pub use trajectory::*;
pub use types::*;
