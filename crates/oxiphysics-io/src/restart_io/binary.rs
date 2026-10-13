// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Little-endian binary restart encoding.
//!
//! The same byte layout backs [`RestartFormat::Binary`](super::RestartFormat),
//! `Hdf5Like` and `MessagePack`; only the 4-byte magic differs.
//!
//! ```text
//! [MAGIC 4B][FORMAT_VERSION u32]
//! [METADATA: version str][timestamp u64][step u64][time f64][crate str][description str]
//! [N_PARTICLES u64][PRESENCE u8]                       (PRESENCE exists from version 2)
//! [POSITIONS n x 3 x f64]
//! [VELOCITIES n x 3 x f64]   if PRESENCE bit 0
//! [FORCES     n x 3 x f64]   if PRESENCE bit 1
//! [MASSES     n x f64]       if PRESENCE bit 2
//! [TYPES      n x u32]       if PRESENCE bit 3
//! [BOX_MATRIX 9 x f64]
//! [N_EXTRA_SCALARS u32] { [name str][count u64][count x f64] }*
//! [N_EXTRA_VECTORS u32] { [name str][count u64][count x 3 x f64] }*
//! ```
//!
//! A `str` is a `u32` byte length followed by UTF-8 bytes. Version 1 files
//! have no `PRESENCE` byte and always contain all five per-particle arrays.
//!
//! Every count read from the file is checked against the number of bytes that
//! remain *before* anything is allocated from it, so a corrupted count yields
//! an error instead of an allocation failure.

use super::{BINARY_MAGIC, FORMAT_VERSION, HDF5LIKE_MAGIC, MSGPACK_MAGIC};
use super::{RestartData, RestartMetadata};

/// `PRESENCE` bit: the `velocities` array is stored.
const PRESENT_VELOCITIES: u8 = 1 << 0;
/// `PRESENCE` bit: the `forces` array is stored.
const PRESENT_FORCES: u8 = 1 << 1;
/// `PRESENCE` bit: the `masses` array is stored.
const PRESENT_MASSES: u8 = 1 << 2;
/// `PRESENCE` bit: the `types` array is stored.
const PRESENT_TYPES: u8 = 1 << 3;
/// All `PRESENCE` bits this version knows about; also the implied value for
/// version 1 files.
const PRESENT_ALL: u8 = PRESENT_VELOCITIES | PRESENT_FORCES | PRESENT_MASSES | PRESENT_TYPES;

/// Encoded size of a 3-vector of `f64`.
const VEC3_BYTES: usize = 24;
/// Smallest encoded size of one extra-array record: name length + element count.
const EXTRA_HEADER_BYTES: usize = 4 + 8;

// ── Encoder ───────────────────────────────────────────────────────────────────

fn put_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn put_f64(buf: &mut Vec<u8>, v: f64) {
    buf.extend_from_slice(&v.to_bits().to_le_bytes());
}

fn put_vec3(buf: &mut Vec<u8>, v: &[f64; 3]) {
    for &c in v {
        put_f64(buf, c);
    }
}

fn put_str(buf: &mut Vec<u8>, s: &str) -> Result<(), String> {
    let len = u32::try_from(s.len()).map_err(|_| {
        format!(
            "string of {} bytes is too long for the restart format",
            s.len()
        )
    })?;
    put_u32(buf, len);
    buf.extend_from_slice(s.as_bytes());
    Ok(())
}

fn put_count(buf: &mut Vec<u8>, count: usize, what: &str) -> Result<(), String> {
    let count = u32::try_from(count)
        .map_err(|_| format!("too many {what} ({count}) for the restart format"))?;
    put_u32(buf, count);
    Ok(())
}

/// Encode `data` with the given 4-byte `magic`.
///
/// Fails if the data is inconsistent ([`RestartData::validate`]) or a string
/// or array count does not fit its `u32` length field.
pub(super) fn encode(data: &RestartData, magic: [u8; 4]) -> Result<Vec<u8>, String> {
    data.validate()?;
    let n = data.n_particles();
    let mut presence = 0u8;
    if !data.velocities.is_empty() {
        presence |= PRESENT_VELOCITIES;
    }
    if !data.forces.is_empty() {
        presence |= PRESENT_FORCES;
    }
    if !data.masses.is_empty() {
        presence |= PRESENT_MASSES;
    }
    if !data.types.is_empty() {
        presence |= PRESENT_TYPES;
    }

    let mut buf = Vec::with_capacity(128 + n * (VEC3_BYTES * 3 + 12));
    buf.extend_from_slice(&magic);
    put_u32(&mut buf, FORMAT_VERSION);
    // metadata
    put_str(&mut buf, &data.metadata.version)?;
    put_u64(&mut buf, data.metadata.timestamp);
    put_u64(&mut buf, data.metadata.step);
    put_f64(&mut buf, data.metadata.time);
    put_str(&mut buf, &data.metadata.crate_name)?;
    put_str(&mut buf, &data.metadata.description)?;
    // particle count + which optional arrays follow
    put_u64(&mut buf, n as u64);
    buf.push(presence);
    for p in &data.positions {
        put_vec3(&mut buf, p);
    }
    for v in &data.velocities {
        put_vec3(&mut buf, v);
    }
    for f in &data.forces {
        put_vec3(&mut buf, f);
    }
    for &m in &data.masses {
        put_f64(&mut buf, m);
    }
    for &t in &data.types {
        put_u32(&mut buf, t);
    }
    // box matrix (9 f64)
    for row in &data.box_matrix {
        for &c in row {
            put_f64(&mut buf, c);
        }
    }
    // extra scalars
    put_count(&mut buf, data.extra_scalars.len(), "extra scalar arrays")?;
    for (name, vals) in &data.extra_scalars {
        put_str(&mut buf, name)?;
        put_u64(&mut buf, vals.len() as u64);
        for &v in vals {
            put_f64(&mut buf, v);
        }
    }
    // extra vectors
    put_count(&mut buf, data.extra_vectors.len(), "extra vector arrays")?;
    for (name, vecs) in &data.extra_vectors {
        put_str(&mut buf, name)?;
        put_u64(&mut buf, vecs.len() as u64);
        for v in vecs {
            put_vec3(&mut buf, v);
        }
    }
    Ok(buf)
}

// ── Decoder ───────────────────────────────────────────────────────────────────

/// Bounds-checked forward reader over a byte slice.
///
/// Invariant: `off <= data.len()`.
struct Cursor<'a> {
    data: &'a [u8],
    off: usize,
}

impl<'a> Cursor<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self { data, off: 0 }
    }

    /// Bytes not yet consumed.
    fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.off)
    }

    fn take(&mut self, len: usize, what: &str) -> Result<&'a [u8], String> {
        let eof = || format!("unexpected EOF reading {what} at offset {}", self.off);
        let end = self.off.checked_add(len).ok_or_else(eof)?;
        let slice = self.data.get(self.off..end).ok_or_else(eof)?;
        self.off = end;
        Ok(slice)
    }

    fn array<const N: usize>(&mut self, what: &str) -> Result<[u8; N], String> {
        let slice = self.take(N, what)?;
        <[u8; N]>::try_from(slice).map_err(|_| format!("cannot read {what}"))
    }

    fn u8(&mut self, what: &str) -> Result<u8, String> {
        self.array::<1>(what).map(|b| b[0])
    }

    fn u32(&mut self, what: &str) -> Result<u32, String> {
        self.array::<4>(what).map(u32::from_le_bytes)
    }

    fn u64(&mut self, what: &str) -> Result<u64, String> {
        self.array::<8>(what).map(u64::from_le_bytes)
    }

    fn f64(&mut self, what: &str) -> Result<f64, String> {
        self.array::<8>(what)
            .map(|b| f64::from_bits(u64::from_le_bytes(b)))
    }

    fn vec3(&mut self, what: &str) -> Result<[f64; 3], String> {
        Ok([self.f64(what)?, self.f64(what)?, self.f64(what)?])
    }

    fn string(&mut self, what: &str) -> Result<String, String> {
        let len = self.u32(what)? as usize;
        let bytes = self.take(len, what)?;
        std::str::from_utf8(bytes)
            .map(str::to_owned)
            .map_err(|e| format!("UTF-8 error in {what}: {e}"))
    }

    /// Validate an element count read from the file.
    ///
    /// `raw` entries of `elem_bytes` bytes each must fit in the unread
    /// remainder; only then is the count returned as a `usize` that is safe to
    /// hand to `Vec::with_capacity`.
    fn checked_count(&self, raw: u64, elem_bytes: usize, what: &str) -> Result<usize, String> {
        let count = usize::try_from(raw)
            .map_err(|_| format!("{what}: count {raw} does not fit in memory"))?;
        let needed = count
            .checked_mul(elem_bytes)
            .ok_or_else(|| format!("{what}: count {raw} overflows the byte size"))?;
        if needed > self.remaining() {
            return Err(format!(
                "{what}: header declares {raw} entries ({needed} bytes) but only {} bytes remain",
                self.remaining()
            ));
        }
        Ok(count)
    }

    fn vec3_array(&mut self, count: usize, what: &str) -> Result<Vec<[f64; 3]>, String> {
        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            out.push(self.vec3(what)?);
        }
        Ok(out)
    }

    fn f64_array(&mut self, count: usize, what: &str) -> Result<Vec<f64>, String> {
        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            out.push(self.f64(what)?);
        }
        Ok(out)
    }

    fn u32_array(&mut self, count: usize, what: &str) -> Result<Vec<u32>, String> {
        let mut out = Vec::with_capacity(count);
        for _ in 0..count {
            out.push(self.u32(what)?);
        }
        Ok(out)
    }
}

/// Decode a binary restart blob (any of the three magics, versions 1 and 2).
pub(super) fn decode(bytes: &[u8]) -> Result<RestartData, String> {
    let mut cur = Cursor::new(bytes);
    let magic = cur.array::<4>("magic")?;
    if magic != BINARY_MAGIC && magic != HDF5LIKE_MAGIC && magic != MSGPACK_MAGIC {
        return Err(format!(
            "not a binary restart file: unknown magic {magic:02x?} (expected OXRS, OXH5 or MPK1)"
        ));
    }
    let file_version = cur.u32("format version")?;
    if file_version == 0 || file_version > FORMAT_VERSION {
        return Err(format!(
            "unsupported restart format version {file_version} (this reader supports 1..={FORMAT_VERSION})"
        ));
    }

    let version = cur.string("metadata version")?;
    let timestamp = cur.u64("timestamp")?;
    let step = cur.u64("step")?;
    let time = cur.f64("time")?;
    let crate_name = cur.string("crate name")?;
    let description = cur.string("description")?;
    let metadata = RestartMetadata {
        version,
        timestamp,
        step,
        time,
        crate_name,
        description,
    };

    let declared_n = cur.u64("particle count")?;
    let presence = if file_version >= 2 {
        cur.u8("array presence flags")?
    } else {
        PRESENT_ALL
    };
    if presence & !PRESENT_ALL != 0 {
        return Err(format!(
            "unknown array presence flags 0x{presence:02x} (a newer writer?)"
        ));
    }

    // Every particle needs at least its position plus each present array.
    let mut bytes_per_particle = VEC3_BYTES;
    if presence & PRESENT_VELOCITIES != 0 {
        bytes_per_particle += VEC3_BYTES;
    }
    if presence & PRESENT_FORCES != 0 {
        bytes_per_particle += VEC3_BYTES;
    }
    if presence & PRESENT_MASSES != 0 {
        bytes_per_particle += 8;
    }
    if presence & PRESENT_TYPES != 0 {
        bytes_per_particle += 4;
    }
    let n = cur.checked_count(declared_n, bytes_per_particle, "particle arrays")?;

    let positions = cur.vec3_array(n, "positions")?;
    let velocities = if presence & PRESENT_VELOCITIES != 0 {
        cur.vec3_array(n, "velocities")?
    } else {
        Vec::new()
    };
    let forces = if presence & PRESENT_FORCES != 0 {
        cur.vec3_array(n, "forces")?
    } else {
        Vec::new()
    };
    let masses = if presence & PRESENT_MASSES != 0 {
        cur.f64_array(n, "masses")?
    } else {
        Vec::new()
    };
    let types = if presence & PRESENT_TYPES != 0 {
        cur.u32_array(n, "types")?
    } else {
        Vec::new()
    };

    let mut box_matrix = [[0.0f64; 3]; 3];
    for row in &mut box_matrix {
        for c in row.iter_mut() {
            *c = cur.f64("box matrix")?;
        }
    }

    let n_scalar_arrays = cur.u32("extra scalar array count")?;
    let n_scalar_arrays = cur.checked_count(
        u64::from(n_scalar_arrays),
        EXTRA_HEADER_BYTES,
        "extra scalar arrays",
    )?;
    let mut extra_scalars = Vec::with_capacity(n_scalar_arrays);
    for _ in 0..n_scalar_arrays {
        let name = cur.string("extra scalar name")?;
        let raw = cur.u64("extra scalar length")?;
        let len = cur.checked_count(raw, 8, "extra scalar values")?;
        extra_scalars.push((name, cur.f64_array(len, "extra scalar values")?));
    }

    let n_vector_arrays = cur.u32("extra vector array count")?;
    let n_vector_arrays = cur.checked_count(
        u64::from(n_vector_arrays),
        EXTRA_HEADER_BYTES,
        "extra vector arrays",
    )?;
    let mut extra_vectors = Vec::with_capacity(n_vector_arrays);
    for _ in 0..n_vector_arrays {
        let name = cur.string("extra vector name")?;
        let raw = cur.u64("extra vector length")?;
        let len = cur.checked_count(raw, VEC3_BYTES, "extra vector values")?;
        extra_vectors.push((name, cur.vec3_array(len, "extra vector values")?));
    }

    Ok(RestartData {
        metadata,
        positions,
        velocities,
        forces,
        masses,
        types,
        box_matrix,
        extra_scalars,
        extra_vectors,
    })
}
