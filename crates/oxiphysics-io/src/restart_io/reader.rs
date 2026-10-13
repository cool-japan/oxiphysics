// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! [`RestartReader`]: parse [`RestartData`] from files or byte buffers.

use super::{BINARY_MAGIC, HDF5LIKE_MAGIC, MSGPACK_MAGIC};
use super::{RestartData, RestartFormat, ascii, binary, json};

/// Reads simulation restart data from files or in-memory buffers.
#[derive(Debug, Clone)]
pub struct RestartReader {
    /// Input file path.
    pub path: String,
}

impl RestartReader {
    /// Create a new `RestartReader` for the given file path.
    pub fn new(path: &str) -> Self {
        Self {
            path: path.to_string(),
        }
    }

    /// Read and parse a restart file from `self.path`.
    ///
    /// The format is detected with [`detect_format`](Self::detect_format) and
    /// routed to the matching decoder: JSON to
    /// [`read_json`](Self::read_json), ASCII to
    /// [`read_ascii`](Self::read_ascii), and the three binary tags to
    /// [`read_binary`](Self::read_binary).
    pub fn read(&self) -> Result<RestartData, String> {
        let bytes = std::fs::read(&self.path)
            .map_err(|e| format!("failed to read '{}': {e}", self.path))?;
        match Self::detect_format(&bytes) {
            RestartFormat::Ascii => Self::read_ascii(text_of(&bytes)?),
            RestartFormat::Json => Self::read_json(text_of(&bytes)?),
            RestartFormat::Binary | RestartFormat::Hdf5Like | RestartFormat::MessagePack => {
                Self::read_binary(&bytes)
            }
        }
    }

    /// Detect the serialization format by inspecting the leading bytes.
    ///
    /// The three binary magics (`OXRS`, `OXH5`, `MPK1`) are recognised first.
    /// Otherwise the content is text: a first non-whitespace character of `{`
    /// (after an optional UTF-8 byte-order mark) means JSON, anything else
    /// ASCII.
    pub fn detect_format(bytes: &[u8]) -> RestartFormat {
        if let Some(magic) = bytes.get(0..4) {
            if magic == BINARY_MAGIC {
                return RestartFormat::Binary;
            }
            if magic == HDF5LIKE_MAGIC {
                return RestartFormat::Hdf5Like;
            }
            if magic == MSGPACK_MAGIC {
                return RestartFormat::MessagePack;
            }
        }
        let text = bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes);
        match text.iter().find(|b| !b.is_ascii_whitespace()) {
            Some(b'{') => RestartFormat::Json,
            _ => RestartFormat::Ascii,
        }
    }

    /// Parse a `RestartData` from a raw binary blob (any of the binary formats).
    ///
    /// Accepts the `OXRS`, `OXH5` and `MPK1` magics and format versions 1
    /// (all five per-particle arrays always present) and 2 (array-presence
    /// bitmask). Every length read from the blob is checked against the bytes
    /// that remain before memory is reserved, so corrupted or malicious input
    /// returns an error instead of aborting.
    pub fn read_binary(bytes: &[u8]) -> Result<RestartData, String> {
        binary::decode(bytes)
    }

    /// Parse a `RestartData` from an ASCII restart text.
    ///
    /// Strict: malformed numbers, wrong column counts, unknown or unterminated
    /// sections and arrays whose length disagrees with `N_PARTICLES` are
    /// errors that carry the offending line number.
    pub fn read_ascii(text: &str) -> Result<RestartData, String> {
        ascii::decode(text)
    }

    /// Parse a `RestartData` from a JSON restart document (the output of
    /// [`RestartWriter::write_json`](super::RestartWriter::write_json)).
    pub fn read_json(text: &str) -> Result<RestartData, String> {
        json::decode(text)
    }
}

/// View `bytes` as UTF-8 text.
fn text_of(bytes: &[u8]) -> Result<&str, String> {
    std::str::from_utf8(bytes).map_err(|e| format!("UTF-8 error: {e}"))
}
