// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! [`RestartWriter`]: serialize [`RestartData`] to bytes or to a file.

use std::path::Path;

use super::{BINARY_MAGIC, HDF5LIKE_MAGIC, MSGPACK_MAGIC};
use super::{RestartData, RestartFormat, ascii, binary, json};

/// Writes simulation restart data to files or in-memory buffers.
///
/// Every encoder first calls [`RestartData::validate`], so inconsistent data
/// is rejected when it is written rather than discovered when the restart is
/// read.
#[derive(Debug, Clone)]
pub struct RestartWriter {
    /// Output file path.
    pub path: String,
    /// Serialization format to use.
    pub format: RestartFormat,
}

impl RestartWriter {
    /// Create a new `RestartWriter` targeting `path` with the given `format`.
    pub fn new(path: &str, format: RestartFormat) -> Self {
        Self {
            path: path.to_string(),
            format,
        }
    }

    /// Serialize `data` according to `self.format` into memory.
    ///
    /// Returns an error string if `data` is inconsistent
    /// ([`RestartData::validate`]) or cannot be encoded.
    pub fn encode(&self, data: &RestartData) -> Result<Vec<u8>, String> {
        data.validate()?;
        match &self.format {
            RestartFormat::Binary => Self::write_binary(data),
            RestartFormat::Ascii => Self::write_ascii(data).map(String::into_bytes),
            RestartFormat::Json => Self::write_json(data).map(String::into_bytes),
            RestartFormat::Hdf5Like => Self::write_hdf5like(data),
            RestartFormat::MessagePack => Self::write_msgpack(data),
        }
    }

    /// Serialize `data` according to `self.format` and write to `self.path`.
    ///
    /// The data is validated and fully encoded before the file system is
    /// touched. The bytes are written to a `<path>.tmp` sibling and renamed
    /// over `self.path`, so a failed or interrupted write never leaves a
    /// truncated file where a good one used to be.
    ///
    /// Returns an error string on validation, encoding or I/O failure.
    pub fn write(&self, data: &RestartData) -> Result<(), String> {
        let bytes = self.encode(data)?;
        write_atomic(&self.path, &bytes)
    }

    /// Encode `data` as a custom little-endian binary blob.
    ///
    /// Layout: `[MAGIC 4B][VERSION 4B][METADATA][N_PARTICLES 8B][PRESENCE 1B]
    ///          [POSITIONS…][VELOCITIES…][FORCES…][MASSES…][TYPES…]
    ///          [BOX_MATRIX 72B][N_EXTRA_SCALARS 4B][…][N_EXTRA_VECS 4B][…]`
    ///
    /// `PRESENCE` is a bitmask (bit 0 velocities, bit 1 forces, bit 2 masses,
    /// bit 3 types); arrays whose bit is clear are not stored.
    pub fn write_binary(data: &RestartData) -> Result<Vec<u8>, String> {
        binary::encode(data, BINARY_MAGIC)
    }

    /// Encode `data` as a human-readable ASCII text block.
    ///
    /// Floating-point values use the shortest decimal form that reads back
    /// bit-exactly.
    pub fn write_ascii(data: &RestartData) -> Result<String, String> {
        ascii::encode(data)
    }

    /// Encode `data` as a JSON document.
    ///
    /// Includes every field of [`RestartData`]; floating-point values use the
    /// shortest decimal form that reads back bit-exactly. Read it back with
    /// [`RestartReader::read_json`](super::RestartReader::read_json).
    pub fn write_json(data: &RestartData) -> Result<String, String> {
        json::encode(data)
    }

    /// Encode `data` as the binary format tagged with the `OXH5` magic.
    ///
    /// This is *not* an HDF5 file: the payload is identical to
    /// [`write_binary`](Self::write_binary), only the first four bytes differ.
    pub fn write_hdf5like(data: &RestartData) -> Result<Vec<u8>, String> {
        binary::encode(data, HDF5LIKE_MAGIC)
    }

    /// Encode `data` as the binary format tagged with the `MPK1` magic.
    ///
    /// This is *not* MessagePack: the payload is identical to
    /// [`write_binary`](Self::write_binary), only the first four bytes differ.
    pub fn write_msgpack(data: &RestartData) -> Result<Vec<u8>, String> {
        binary::encode(data, MSGPACK_MAGIC)
    }
}

/// Write `bytes` to `path` through a temporary sibling file and a rename.
fn write_atomic(path: &str, bytes: &[u8]) -> Result<(), String> {
    let target = Path::new(path);
    let file_name = target
        .file_name()
        .ok_or_else(|| format!("failed to write '{path}': path has no file name"))?;
    let mut tmp_name = file_name.to_os_string();
    tmp_name.push(".tmp");
    let tmp = target.with_file_name(tmp_name);

    if let Err(e) = std::fs::write(&tmp, bytes) {
        // Best-effort cleanup of a partially written temporary; the write
        // error below is what the caller needs to see.
        let _ = std::fs::remove_file(&tmp);
        return Err(format!("failed to write '{}': {e}", tmp.display()));
    }
    if let Err(e) = std::fs::rename(&tmp, target) {
        let _ = std::fs::remove_file(&tmp);
        return Err(format!("failed to write '{path}': {e}"));
    }
    Ok(())
}
