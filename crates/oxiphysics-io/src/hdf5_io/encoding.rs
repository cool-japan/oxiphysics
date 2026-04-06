// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Byte-layout simulation for object headers and superblock encoding.

#![allow(dead_code)]

use super::file::{Hdf5ObjectHeader, Hdf5Superblock};
use super::types::{Hdf5Error, Hdf5Result};

// ---------------------------------------------------------------------------
// Object header encoding
// ---------------------------------------------------------------------------

/// Simulate encoding an object header into a byte vector.
///
/// Real HDF5 uses a complex binary format; here we produce a human-readable
/// ASCII summary padded to a fixed size for unit-testing purposes.
pub fn encode_object_header(header: &Hdf5ObjectHeader) -> Vec<u8> {
    let raw = format!(
        "OBJ|type={}|addr={}|msgs={}|size={}",
        header.object_type, header.address, header.n_messages, header.header_size
    );
    let mut out = raw.into_bytes();
    out.resize(64, 0);
    out
}

/// Decode an object header from bytes (inverse of `encode_object_header`).
pub fn decode_object_header(bytes: &[u8]) -> Hdf5Result<Hdf5ObjectHeader> {
    let s = std::str::from_utf8(bytes)
        .map_err(|_| Hdf5Error::Generic("header bytes are not valid UTF-8".to_string()))?
        .trim_end_matches('\0');
    let parts: Vec<&str> = s.split('|').collect();
    if parts.len() < 5 || parts[0] != "OBJ" {
        return Err(Hdf5Error::Generic(format!("invalid header: {s}")));
    }
    let parse_kv =
        |kv: &str| -> Option<String> { kv.split_once('=').map(|x| x.1).map(|v| v.to_string()) };
    let object_type = parse_kv(parts[1]).unwrap_or_default();
    let address: u64 = parse_kv(parts[2]).and_then(|v| v.parse().ok()).unwrap_or(0);
    let n_messages: u32 = parse_kv(parts[3]).and_then(|v| v.parse().ok()).unwrap_or(0);
    let header_size: u32 = parse_kv(parts[4].trim_end_matches('\0'))
        .and_then(|v| v.parse().ok())
        .unwrap_or(0);
    Ok(Hdf5ObjectHeader {
        object_type,
        address,
        n_messages,
        header_size,
    })
}

// ---------------------------------------------------------------------------
// Superblock encoding
// ---------------------------------------------------------------------------

/// Encode a superblock as a fixed 64-byte byte vector (mock layout).
pub fn encode_superblock(sb: &Hdf5Superblock) -> Vec<u8> {
    let raw = format!(
        "SB|v={}|fs={}|roff={}|eof={}|sl={}|so={}",
        sb.version,
        sb.file_size,
        sb.root_obj_header_offset,
        sb.eof_address,
        sb.size_of_lengths,
        sb.size_of_offsets
    );
    let mut out = raw.into_bytes();
    out.resize(128, 0);
    out
}

/// Decode a superblock from bytes.
pub fn decode_superblock(bytes: &[u8]) -> Hdf5Result<Hdf5Superblock> {
    let s = std::str::from_utf8(bytes)
        .map_err(|_| Hdf5Error::Generic("superblock bytes are not valid UTF-8".to_string()))?
        .trim_end_matches('\0');
    let parts: Vec<&str> = s.split('|').collect();
    if parts.len() < 7 || parts[0] != "SB" {
        return Err(Hdf5Error::Generic(format!("invalid superblock: {s}")));
    }
    let parse_u8 = |kv: &str| -> u8 {
        kv.split_once('=')
            .map(|x| x.1)
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    };
    let parse_u64 = |kv: &str| -> u64 {
        kv.split_once('=')
            .map(|x| x.1)
            .and_then(|v| v.parse().ok())
            .unwrap_or(0)
    };
    Ok(Hdf5Superblock {
        version: parse_u8(parts[1]),
        file_size: parse_u64(parts[2]),
        root_obj_header_offset: parse_u64(parts[3]),
        eof_address: parse_u64(parts[4]),
        size_of_lengths: parse_u8(parts[5]),
        size_of_offsets: parse_u8(parts[6].trim_end_matches('\0')),
    })
}
