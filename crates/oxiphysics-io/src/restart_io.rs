// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation restart and checkpoint I/O.
//!
//! Provides writers and readers for simulation restart files in three
//! encodings that each round-trip every field of [`RestartData`]:
//!
//! * [`RestartFormat::Binary`] — compact little-endian binary,
//! * [`RestartFormat::Ascii`] — line-oriented, human-readable text,
//! * [`RestartFormat::Json`] — a JSON document.
//!
//! Floating-point values are stored losslessly: the binary format keeps the raw
//! IEEE-754 bits, and the text formats use the shortest decimal representation
//! that parses back to the identical `f64` (`-0.0`, subnormals and `1e300`
//! survive; NaN and infinities survive as values, but not their NaN payload
//! bits). Per-particle arrays other than `positions` may be absent, which every
//! encoding represents explicitly; see [`RestartData::validate`].
//!
//! [`RestartFormat::Hdf5Like`] and [`RestartFormat::MessagePack`] are the
//! binary encoding under a different 4-byte magic number; they contain no HDF5
//! or MessagePack structure.
//!
//! Also provided are a rolling [`CheckpointManager`] and a
//! [`RestartValidator`] for checksum verification.
//!
//! The implementation is split by concern: `binary`, `ascii` and `json` hold
//! the codecs, `writer` / `reader` the public entry points, `checkpoint` the
//! checkpoint bookkeeping and `convert` the one-way exports.

mod ascii;
mod binary;
mod checkpoint;
mod convert;
mod json;
mod reader;
mod writer;

#[cfg(test)]
mod tests;

pub use checkpoint::{CheckpointManager, IncrementalRestart, RestartValidator};
pub use convert::{restart_to_lammps_dump, restart_to_xyz};
pub use reader::RestartReader;
pub use writer::RestartWriter;

// ── Magic bytes ───────────────────────────────────────────────────────────────

/// Magic bytes written at the start of every binary restart file: "OXRS".
const BINARY_MAGIC: [u8; 4] = [0x4F, 0x58, 0x52, 0x53];

/// Magic bytes of the [`RestartFormat::Hdf5Like`] tag: "OXH5".
const HDF5LIKE_MAGIC: [u8; 4] = [0x4F, 0x58, 0x48, 0x35];

/// Magic bytes of the [`RestartFormat::MessagePack`] tag: "MPK1".
const MSGPACK_MAGIC: [u8; 4] = [0x4D, 0x50, 0x4B, 0x31];

/// Schema version written by every encoder.
///
/// Version 1 stored all five per-particle arrays unconditionally (binary
/// only). Version 2 adds the array-presence bitmask to the binary header so
/// absent arrays are skipped; the binary reader still accepts version 1.
const FORMAT_VERSION: u32 = 2;

// ── RestartFormat ─────────────────────────────────────────────────────────────

/// Serialization format for restart / checkpoint files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RestartFormat {
    /// Custom little-endian binary format (fastest, smallest).
    Binary,
    /// Human-readable ASCII key-value format.
    Ascii,
    /// The [`Binary`](Self::Binary) encoding tagged with the magic bytes
    /// `OXH5`. It is *not* an HDF5 file and has no group/dataset structure;
    /// the tag only lets callers label a file as belonging to an HDF5-style
    /// pipeline.
    Hdf5Like,
    /// JSON text format (portable, verbose).
    Json,
    /// The [`Binary`](Self::Binary) encoding tagged with the magic bytes
    /// `MPK1`. It is *not* MessagePack and cannot be read by MessagePack
    /// tools.
    MessagePack,
}

// ── RestartMetadata ───────────────────────────────────────────────────────────

/// Metadata stored alongside a restart snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct RestartMetadata {
    /// File-format version string (e.g. `"1.0"`).
    pub version: String,
    /// UNIX timestamp (seconds since epoch) when this file was written.
    pub timestamp: u64,
    /// Simulation step number.
    pub step: u64,
    /// Simulation time (in simulation units).
    pub time: f64,
    /// Name of the crate that produced this file.
    pub crate_name: String,
    /// Free-form human-readable description.
    pub description: String,
}

impl RestartMetadata {
    /// Construct a new `RestartMetadata`.
    pub fn new(
        version: impl Into<String>,
        timestamp: u64,
        step: u64,
        time: f64,
        crate_name: impl Into<String>,
        description: impl Into<String>,
    ) -> Self {
        Self {
            version: version.into(),
            timestamp,
            step,
            time,
            crate_name: crate_name.into(),
            description: description.into(),
        }
    }

    /// Return a default metadata for testing.
    pub fn default_test() -> Self {
        Self::new("1.0", 0, 0, 0.0, "oxiphysics", "test checkpoint")
    }
}

// ── RestartData ───────────────────────────────────────────────────────────────

/// Full simulation state stored in a restart file.
#[derive(Debug, Clone, PartialEq)]
pub struct RestartData {
    /// File metadata.
    pub metadata: RestartMetadata,
    /// Particle positions \[x, y, z\] (simulation units).
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities \[vx, vy, vz\].
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle forces \[fx, fy, fz\].
    pub forces: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// Per-particle type indices.
    pub types: Vec<u32>,
    /// Simulation box matrix (row-major: rows are box vectors a, b, c).
    pub box_matrix: [[f64; 3]; 3],
    /// Named extra scalar arrays `(name, values)`.
    pub extra_scalars: Vec<(String, Vec<f64>)>,
    /// Named extra vector arrays `(name, vectors)`.
    pub extra_vectors: Vec<(String, Vec<[f64; 3]>)>,
}

impl RestartData {
    /// Number of particles in this snapshot.
    pub fn n_particles(&self) -> usize {
        self.positions.len()
    }

    /// Construct a minimal empty `RestartData`.
    pub fn empty(metadata: RestartMetadata) -> Self {
        Self {
            metadata,
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            types: Vec::new(),
            box_matrix: [[0.0; 3]; 3],
            extra_scalars: Vec::new(),
            extra_vectors: Vec::new(),
        }
    }

    /// Build a simple single-particle `RestartData` for tests.
    pub fn single_particle_test() -> Self {
        let meta = RestartMetadata::default_test();
        let mut d = Self::empty(meta);
        d.positions = vec![[1.0, 2.0, 3.0]];
        d.velocities = vec![[0.1, 0.2, 0.3]];
        d.forces = vec![[0.0, -9.81, 0.0]];
        d.masses = vec![1.0];
        d.types = vec![0];
        d.box_matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        d
    }

    /// Check that the per-particle arrays are mutually consistent.
    ///
    /// `positions` defines the particle count `n`
    /// ([`n_particles`](Self::n_particles)). Each of `velocities`, `forces`,
    /// `masses` and `types` must either be absent (empty) or hold exactly `n`
    /// entries; any other length cannot be stored under a single particle
    /// count and is reported as an error. Every `RestartWriter` encoder calls
    /// this before producing any output.
    ///
    /// The extra scalar / vector arrays are free-form and are not checked.
    pub fn validate(&self) -> Result<(), String> {
        let n = self.positions.len();
        let arrays = [
            ("velocities", self.velocities.len()),
            ("forces", self.forces.len()),
            ("masses", self.masses.len()),
            ("types", self.types.len()),
        ];
        for (name, len) in arrays {
            if len != 0 && len != n {
                return Err(format!(
                    "inconsistent restart data: `{name}` has {len} entries but `positions` defines {n} particles (each per-particle array must be empty or have one entry per particle)"
                ));
            }
        }
        Ok(())
    }
}
