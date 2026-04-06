#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Binary format I/O for physics simulation data.
//!
//! Provides binary trajectory writers/readers, DCD format support,
//! minimal XTC-like compression, energy log writers, and checkpoint I/O.

#![allow(dead_code)]
#![allow(clippy::too_many_arguments)]

// ── BinaryHeader ─────────────────────────────────────────────────────────────

/// Magic bytes identifying the OxiPhysics binary trajectory format.
pub const OXIPHY_MAGIC: [u8; 4] = *b"OXIP";

/// Binary file header for OxiPhysics trajectory files.
///
/// Layout (little-endian):
/// - 4 bytes  : magic
/// - 4 bytes  : version (u32)
/// - 8 bytes  : n_particles (u64)
/// - 8 bytes  : n_frames (u64)
/// - 8 bytes  : dt (f64)
#[derive(Debug, Clone, PartialEq)]
pub struct BinaryHeader {
    /// Magic bytes — must equal `OXIPHY_MAGIC`.
    pub magic: [u8; 4],
    /// File format version.
    pub version: u32,
    /// Number of particles per frame.
    pub n_particles: u64,
    /// Number of frames stored.
    pub n_frames: u64,
    /// Integration time step (seconds).
    pub dt: f64,
}

impl BinaryHeader {
    /// Byte size of the serialised header.
    pub const SIZE: usize = 32;

    /// Create a new header with the given parameters.
    pub fn new(n_particles: u64, n_frames: u64, dt: f64) -> Self {
        Self {
            magic: OXIPHY_MAGIC,
            version: 1,
            n_particles,
            n_frames,
            dt,
        }
    }

    /// Validate that the magic bytes and version are recognised.
    ///
    /// Returns `true` when the header looks valid.
    pub fn validate(&self) -> bool {
        self.magic == OXIPHY_MAGIC && self.version >= 1
    }

    /// Serialise this header into `buf` (appends bytes).
    pub fn write(&self, buf: &mut Vec<u8>) {
        buf.extend_from_slice(&self.magic);
        buf.extend_from_slice(&self.version.to_le_bytes());
        buf.extend_from_slice(&self.n_particles.to_le_bytes());
        buf.extend_from_slice(&self.n_frames.to_le_bytes());
        buf.extend_from_slice(&self.dt.to_le_bytes());
    }

    /// Deserialise a header from a byte slice.
    ///
    /// Returns `None` if `buf` is shorter than [`BinaryHeader::SIZE`].
    pub fn read(buf: &[u8]) -> Option<Self> {
        if buf.len() < Self::SIZE {
            return None;
        }
        let magic = [buf[0], buf[1], buf[2], buf[3]];
        let version = u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]);
        let n_particles = u64::from_le_bytes(buf[8..16].try_into().ok()?);
        let n_frames = u64::from_le_bytes(buf[16..24].try_into().ok()?);
        let dt = f64::from_le_bytes(buf[24..32].try_into().ok()?);
        Some(Self {
            magic,
            version,
            n_particles,
            n_frames,
            dt,
        })
    }
}

// ── ParticleFrame ─────────────────────────────────────────────────────────────

/// A single snapshot of particle state at one simulation step.
#[derive(Debug, Clone, PartialEq)]
pub struct ParticleFrame {
    /// Simulation step index.
    pub step: u64,
    /// Simulation time (seconds).
    pub time: f64,
    /// Particle positions `[x, y, z]` in metres.
    pub positions: Vec<[f32; 3]>,
    /// Particle velocities `[vx, vy, vz]` in m/s.
    pub velocities: Vec<[f32; 3]>,
    /// Particle masses in kg.
    pub masses: Vec<f32>,
}

impl ParticleFrame {
    /// Create a new empty frame for `n` particles.
    pub fn new(step: u64, time: f64, n: usize) -> Self {
        Self {
            step,
            time,
            positions: vec![[0.0; 3]; n],
            velocities: vec![[0.0; 3]; n],
            masses: vec![1.0; n],
        }
    }

    /// Byte size of a serialised frame for `n` particles.
    pub fn serialized_size(n: usize) -> usize {
        // step(8) + time(8) + n_pos(n*12) + n_vel(n*12) + n_mass(n*4)
        8 + 8 + n * 12 + n * 12 + n * 4
    }

    /// Serialise this frame into `buf` (appends bytes).
    pub fn serialize(&self, buf: &mut Vec<u8>) {
        buf.extend_from_slice(&self.step.to_le_bytes());
        buf.extend_from_slice(&self.time.to_le_bytes());
        for p in &self.positions {
            buf.extend_from_slice(&p[0].to_le_bytes());
            buf.extend_from_slice(&p[1].to_le_bytes());
            buf.extend_from_slice(&p[2].to_le_bytes());
        }
        for v in &self.velocities {
            buf.extend_from_slice(&v[0].to_le_bytes());
            buf.extend_from_slice(&v[1].to_le_bytes());
            buf.extend_from_slice(&v[2].to_le_bytes());
        }
        for &m in &self.masses {
            buf.extend_from_slice(&m.to_le_bytes());
        }
    }

    /// Deserialise a frame from `buf` starting at byte 0, given `n` particles.
    ///
    /// Returns `None` if `buf` is too short.
    pub fn deserialize(buf: &[u8], n: usize) -> Option<Self> {
        let needed = Self::serialized_size(n);
        if buf.len() < needed {
            return None;
        }
        let step = u64::from_le_bytes(buf[0..8].try_into().ok()?);
        let time = f64::from_le_bytes(buf[8..16].try_into().ok()?);
        let mut offset = 16usize;
        let mut positions = Vec::with_capacity(n);
        for _ in 0..n {
            let x = f32::from_le_bytes(buf[offset..offset + 4].try_into().ok()?);
            let y = f32::from_le_bytes(buf[offset + 4..offset + 8].try_into().ok()?);
            let z = f32::from_le_bytes(buf[offset + 8..offset + 12].try_into().ok()?);
            positions.push([x, y, z]);
            offset += 12;
        }
        let mut velocities = Vec::with_capacity(n);
        for _ in 0..n {
            let x = f32::from_le_bytes(buf[offset..offset + 4].try_into().ok()?);
            let y = f32::from_le_bytes(buf[offset + 4..offset + 8].try_into().ok()?);
            let z = f32::from_le_bytes(buf[offset + 8..offset + 12].try_into().ok()?);
            velocities.push([x, y, z]);
            offset += 12;
        }
        let mut masses = Vec::with_capacity(n);
        for _ in 0..n {
            let m = f32::from_le_bytes(buf[offset..offset + 4].try_into().ok()?);
            masses.push(m);
            offset += 4;
        }
        Some(Self {
            step,
            time,
            positions,
            velocities,
            masses,
        })
    }
}

// ── BinaryTrajectoryWriter ────────────────────────────────────────────────────

/// Streaming writer for the OxiPhysics binary trajectory format.
///
/// Call [`write_frame`](BinaryTrajectoryWriter::write_frame) for each step,
/// then [`finalize`](BinaryTrajectoryWriter::finalize) to obtain the complete
/// byte buffer (header `n_frames` is patched at finalization).
#[derive(Debug)]
pub struct BinaryTrajectoryWriter {
    buf: Vec<u8>,
    /// Embedded file header.
    pub header: BinaryHeader,
    frames_written: u64,
}

impl BinaryTrajectoryWriter {
    /// Create a writer for the given number of particles and time step.
    pub fn new(n_particles: u64, dt: f64) -> Self {
        let header = BinaryHeader::new(n_particles, 0, dt);
        let mut buf = Vec::new();
        header.write(&mut buf); // placeholder; patched in finalize
        Self {
            buf,
            header,
            frames_written: 0,
        }
    }

    /// Append one frame to the internal buffer.
    pub fn write_frame(&mut self, frame: &ParticleFrame) {
        frame.serialize(&mut self.buf);
        self.frames_written += 1;
    }

    /// Finalise: patch `n_frames` in the header and return the buffer.
    pub fn finalize(mut self) -> Vec<u8> {
        // Patch n_frames at byte offset 16..24 in the header
        let nf = self.frames_written.to_le_bytes();
        self.buf[16..24].copy_from_slice(&nf);
        self.buf
    }

    /// Number of frames written so far.
    pub fn frame_count(&self) -> u64 {
        self.frames_written
    }
}

// ── BinaryTrajectoryReader ────────────────────────────────────────────────────

/// Reader for the OxiPhysics binary trajectory format.
#[derive(Debug)]
pub struct BinaryTrajectoryReader {
    buf: Vec<u8>,
    cursor: usize,
    /// Parsed file header.
    pub header: BinaryHeader,
}

impl BinaryTrajectoryReader {
    /// Create a reader from a byte buffer.
    ///
    /// Returns `None` if the header is invalid or magic bytes mismatch.
    pub fn new(buf: Vec<u8>) -> Option<Self> {
        let header = BinaryHeader::read(&buf)?;
        if !header.validate() {
            return None;
        }
        Some(Self {
            buf,
            cursor: BinaryHeader::SIZE,
            header,
        })
    }

    /// Read the next frame sequentially.
    ///
    /// Returns `None` when no more frames are available.
    pub fn read_frame(&mut self) -> Option<ParticleFrame> {
        let n = self.header.n_particles as usize;
        let needed = ParticleFrame::serialized_size(n);
        if self.cursor + needed > self.buf.len() {
            return None;
        }
        let frame = ParticleFrame::deserialize(&self.buf[self.cursor..], n)?;
        self.cursor += needed;
        Some(frame)
    }

    /// Seek to frame index `i` (0-based).
    ///
    /// Does nothing if `i` is out of range.
    pub fn seek_frame(&mut self, i: usize) {
        let n = self.header.n_particles as usize;
        let frame_size = ParticleFrame::serialized_size(n);
        let new_cursor = BinaryHeader::SIZE + i * frame_size;
        if new_cursor <= self.buf.len() {
            self.cursor = new_cursor;
        }
    }

    /// Return the number of frames encoded in the header.
    pub fn n_frames(&self) -> u64 {
        self.header.n_frames
    }
}

// ── DcdWriter ────────────────────────────────────────────────────────────────

/// Writer that produces CHARMM DCD trajectory files as a byte buffer.
///
/// The DCD format uses Fortran-style 4-byte record markers before and after
/// each block.
#[derive(Debug)]
pub struct DcdWriter {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of frames written so far.
    pub n_frames: usize,
    buf: Vec<u8>,
    /// Byte offset of the `NFILE` field so it can be patched at the end.
    nfile_offset: usize,
}

fn write_i32_le(buf: &mut Vec<u8>, v: i32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn write_f32_le(buf: &mut Vec<u8>, v: f32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn fortran_record(buf: &mut Vec<u8>, data: &[u8]) {
    let len = data.len() as i32;
    write_i32_le(buf, len);
    buf.extend_from_slice(data);
    write_i32_le(buf, len);
}

impl DcdWriter {
    /// Create a new DCD writer for `n_atoms` atoms and the given time step.
    pub fn new(n_atoms: usize, dt: f32) -> Self {
        let mut buf = Vec::new();
        // CORD block
        let mut hdr = Vec::new();
        hdr.extend_from_slice(b"CORD");
        // NFILE placeholder (patched later)
        let nfile_offset_in_hdr = hdr.len();
        write_i32_le(&mut hdr, 0); // NFILE (frames)
        write_i32_le(&mut hdr, 0); // ISTART
        write_i32_le(&mut hdr, 1); // NSAVC
        write_i32_le(&mut hdr, 0); // NSTEP
        write_i32_le(&mut hdr, 0); // 0
        write_i32_le(&mut hdr, 0); // 0
        write_i32_le(&mut hdr, 0); // 0
        write_i32_le(&mut hdr, 0); // 0
        write_i32_le(&mut hdr, 0); // 0
        write_f32_le(&mut hdr, dt); // DELTA
        for _ in 0..9 {
            write_i32_le(&mut hdr, 0); // padding
        }
        write_i32_le(&mut hdr, 24); // CHARMM version marker
        // The NFILE offset in the full buffer = 4 (Fortran record len) + 4 (CORD) + nfile_offset_in_hdr
        let nfile_offset = 4 + nfile_offset_in_hdr;
        fortran_record(&mut buf, &hdr);

        // TITLE block
        let mut title_block = Vec::new();
        write_i32_le(&mut title_block, 1); // NTITLE
        let title = b"OxiPhysics DCD                                                  ";
        title_block.extend_from_slice(&title[..80.min(title.len())]);
        if title.len() < 80 {
            title_block.extend(std::iter::repeat_n(b' ', 80 - title.len()));
        }
        fortran_record(&mut buf, &title_block);

        // NATOM block
        let mut natom_block = Vec::new();
        write_i32_le(&mut natom_block, n_atoms as i32);
        fortran_record(&mut buf, &natom_block);

        Self {
            n_atoms,
            n_frames: 0,
            buf,
            nfile_offset,
        }
    }

    /// Append one frame (x, y, z coordinate arrays) to the buffer.
    ///
    /// Each array must have length `n_atoms`.
    pub fn write_frame(&mut self, x: &[f32], y: &[f32], z: &[f32]) {
        let write_coord = |buf: &mut Vec<u8>, coords: &[f32]| {
            let mut data = Vec::with_capacity(coords.len() * 4);
            for &v in coords {
                write_f32_le(&mut data, v);
            }
            fortran_record(buf, &data);
        };
        write_coord(&mut self.buf, x);
        write_coord(&mut self.buf, y);
        write_coord(&mut self.buf, z);
        self.n_frames += 1;
    }

    /// Finalise and return the complete DCD byte buffer.
    ///
    /// Patches the `NFILE` field in the header with the actual frame count.
    pub fn finalize(mut self) -> Vec<u8> {
        let nf = (self.n_frames as i32).to_le_bytes();
        self.buf[self.nfile_offset..self.nfile_offset + 4].copy_from_slice(&nf);
        self.buf
    }
}

// ── DcdReader ────────────────────────────────────────────────────────────────

/// Reader for CHARMM DCD trajectory files (pure byte-slice, no I/O).
#[derive(Debug)]
pub struct DcdReader {
    buf: Vec<u8>,
    /// Number of atoms (from header).
    pub n_atoms: usize,
    /// Number of frames (from header).
    pub n_frames: usize,
    /// Byte offset where the first frame begins.
    frames_start: usize,
    /// Byte size of one frame.
    frame_size: usize,
}

impl DcdReader {
    /// Parse the DCD header from `buf` and return a reader, or `None` on error.
    pub fn parse_header(buf: Vec<u8>) -> Option<Self> {
        if buf.len() < 8 {
            return None;
        }
        // Fortran record: 4-byte length prefix
        let rec_len = i32::from_le_bytes(buf[0..4].try_into().ok()?) as usize;
        if buf.len() < 4 + rec_len + 4 {
            return None;
        }
        let hdr = &buf[4..4 + rec_len];
        if &hdr[0..4] != b"CORD" {
            return None;
        }
        let n_frames = i32::from_le_bytes(hdr[4..8].try_into().ok()?) as usize;
        // Skip to end of first Fortran record
        let mut offset = 4 + rec_len + 4;
        // Skip TITLE record
        if offset + 4 > buf.len() {
            return None;
        }
        let title_len = i32::from_le_bytes(buf[offset..offset + 4].try_into().ok()?) as usize;
        offset += 4 + title_len + 4;
        // NATOM record
        if offset + 4 > buf.len() {
            return None;
        }
        let natom_len = i32::from_le_bytes(buf[offset..offset + 4].try_into().ok()?) as usize;
        offset += 4;
        if natom_len < 4 || offset + 4 > buf.len() {
            return None;
        }
        let n_atoms = i32::from_le_bytes(buf[offset..offset + 4].try_into().ok()?) as usize;
        offset += natom_len + 4;
        let frames_start = offset;
        // Each frame: 3 coord arrays, each with Fortran record wrapper
        let coord_block = 4 + n_atoms * 4 + 4; // len_prefix + data + len_suffix
        let frame_size = 3 * coord_block;
        Some(Self {
            buf,
            n_atoms,
            n_frames,
            frames_start,
            frame_size,
        })
    }

    /// Read the coordinates for frame `frame_idx` (0-based).
    ///
    /// Returns `(x, y, z)` each as `Vec`f32`, or `None` if out of range.
    pub fn read_frame(&self, frame_idx: usize) -> Option<(Vec<f32>, Vec<f32>, Vec<f32>)> {
        if frame_idx >= self.n_frames {
            return None;
        }
        let start = self.frames_start + frame_idx * self.frame_size;
        let read_coord = |offset: usize| -> Option<Vec<f32>> {
            let _rec_len =
                i32::from_le_bytes(self.buf[offset..offset + 4].try_into().ok()?) as usize;
            let data_start = offset + 4;
            let mut coords = Vec::with_capacity(self.n_atoms);
            for i in 0..self.n_atoms {
                let o = data_start + i * 4;
                let v = f32::from_le_bytes(self.buf[o..o + 4].try_into().ok()?);
                coords.push(v);
            }
            Some(coords)
        };
        let coord_block = 4 + self.n_atoms * 4 + 4;
        let x = read_coord(start)?;
        let y = read_coord(start + coord_block)?;
        let z = read_coord(start + 2 * coord_block)?;
        Some((x, y, z))
    }
}

// ── XtcEncoder ───────────────────────────────────────────────────────────────

/// Minimal XTC-like lossy compression for particle positions.
///
/// Positions are quantised to integer multiples of `1e-3` nm (the default
/// XTC precision) and stored as big-endian i32 values, preceded by a 4-byte
/// particle count.
#[derive(Debug, Default)]
pub struct XtcEncoder;

impl XtcEncoder {
    /// Create a new encoder.
    pub fn new() -> Self {
        Self
    }

    /// Compress `positions` into a byte buffer using integer quantisation.
    ///
    /// The precision is fixed at 1000 (i.e. 3 decimal places in nm units).
    pub fn compress_frame(positions: &[[f32; 3]]) -> Vec<u8> {
        let precision: f32 = 1000.0;
        let mut buf = Vec::new();
        let n = positions.len() as u32;
        buf.extend_from_slice(&n.to_be_bytes());
        buf.extend_from_slice(&precision.to_bits().to_be_bytes());
        for p in positions {
            for &v in p.iter() {
                let q = (v * precision).round() as i32;
                buf.extend_from_slice(&q.to_be_bytes());
            }
        }
        buf
    }

    /// Decompress a byte buffer produced by [`compress_frame`](XtcEncoder::compress_frame).
    ///
    /// Returns `None` if `buf` is malformed.
    pub fn decompress_frame(buf: &[u8], n: usize) -> Option<Vec<[f32; 3]>> {
        if buf.len() < 8 {
            return None;
        }
        let _stored_n = u32::from_be_bytes(buf[0..4].try_into().ok()?);
        let precision = f32::from_bits(u32::from_be_bytes(buf[4..8].try_into().ok()?));
        if precision.abs() < 1e-9 {
            return None;
        }
        let needed = 8 + n * 12;
        if buf.len() < needed {
            return None;
        }
        let mut positions = Vec::with_capacity(n);
        let mut offset = 8usize;
        for _ in 0..n {
            let qx = i32::from_be_bytes(buf[offset..offset + 4].try_into().ok()?);
            let qy = i32::from_be_bytes(buf[offset + 4..offset + 8].try_into().ok()?);
            let qz = i32::from_be_bytes(buf[offset + 8..offset + 12].try_into().ok()?);
            offset += 12;
            positions.push([
                qx as f32 / precision,
                qy as f32 / precision,
                qz as f32 / precision,
            ]);
        }
        Some(positions)
    }
}

// ── EnergyLogWriter ──────────────────────────────────────────────────────────

/// A step entry in the energy log.
#[derive(Debug, Clone)]
pub struct EnergyEntry {
    /// Simulation step index.
    pub step: u64,
    /// Kinetic energy (J or reduced units).
    pub ke: f64,
    /// Potential energy (J or reduced units).
    pub pe: f64,
    /// Temperature (K).
    pub temp: f64,
    /// Pressure (Pa or reduced units).
    pub pressure: f64,
}

/// Writer for CSV-like energy time series.
///
/// Accumulates entries in memory and can produce a CSV string.
#[derive(Debug, Default)]
pub struct EnergyLogWriter {
    entries: Vec<EnergyEntry>,
}

impl EnergyLogWriter {
    /// Create a new empty energy log writer.
    pub fn new() -> Self {
        Self {
            entries: Vec::new(),
        }
    }

    /// Append one time step to the log.
    pub fn write_step(&mut self, step: u64, ke: f64, pe: f64, temp: f64, pressure: f64) {
        self.entries.push(EnergyEntry {
            step,
            ke,
            pe,
            temp,
            pressure,
        });
    }

    /// Render the log as a CSV string with a header row.
    pub fn to_csv_string(&self) -> String {
        let mut out = String::from("step,ke,pe,total_energy,temperature,pressure\n");
        for e in &self.entries {
            out.push_str(&format!(
                "{},{:.6},{:.6},{:.6},{:.6},{:.6}\n",
                e.step,
                e.ke,
                e.pe,
                e.ke + e.pe,
                e.temp,
                e.pressure
            ));
        }
        out
    }

    /// Number of entries recorded.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Returns `true` if no entries have been written.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

// ── CheckpointWriter ─────────────────────────────────────────────────────────

/// Serialises and deserialises simulation checkpoints as byte buffers.
///
/// A checkpoint stores positions, velocities, the current step, and the
/// current simulation time. All values are encoded as little-endian bytes.
#[derive(Debug, Default)]
pub struct CheckpointWriter;

impl CheckpointWriter {
    /// Create a new checkpoint writer/reader.
    pub fn new() -> Self {
        Self
    }

    /// Serialise a simulation state to a byte buffer.
    ///
    /// * `positions`  – `N` particle positions `\[x, y, z\]` (f64 metres)
    /// * `velocities` – `N` particle velocities `\[vx, vy, vz\]` (f64 m/s)
    /// * `step`       – current step index
    /// * `time`       – current simulation time (seconds)
    pub fn save_state(
        positions: &[[f64; 3]],
        velocities: &[[f64; 3]],
        step: u64,
        time: f64,
    ) -> Vec<u8> {
        let n = positions.len();
        // Header: magic(4) + n(8) + step(8) + time(8) = 28 bytes
        // Data: n * (3*8 + 3*8) = n * 48 bytes
        let mut buf = Vec::with_capacity(28 + n * 48);
        buf.extend_from_slice(b"OXCK"); // magic
        buf.extend_from_slice(&(n as u64).to_le_bytes());
        buf.extend_from_slice(&step.to_le_bytes());
        buf.extend_from_slice(&time.to_le_bytes());
        for p in positions {
            for &v in p.iter() {
                buf.extend_from_slice(&v.to_le_bytes());
            }
        }
        for v in velocities {
            for &c in v.iter() {
                buf.extend_from_slice(&c.to_le_bytes());
            }
        }
        buf
    }

    /// Deserialise a checkpoint from `buf`.
    ///
    /// Returns `(positions, velocities, step, time)` or panics if the buffer
    /// is too short or the magic bytes are wrong.
    pub fn load_state(buf: &[u8]) -> (Vec<[f64; 3]>, Vec<[f64; 3]>, u64, f64) {
        assert!(buf.len() >= 28, "checkpoint buffer too short");
        assert_eq!(&buf[0..4], b"OXCK", "bad checkpoint magic");
        let n =
            u64::from_le_bytes(buf[4..12].try_into().expect("slice length must match")) as usize;
        let step = u64::from_le_bytes(buf[12..20].try_into().expect("slice length must match"));
        let time = f64::from_le_bytes(buf[20..28].try_into().expect("slice length must match"));
        let mut offset = 28usize;
        let mut positions = Vec::with_capacity(n);
        for _ in 0..n {
            let x = f64::from_le_bytes(
                buf[offset..offset + 8]
                    .try_into()
                    .expect("slice length must match"),
            );
            let y = f64::from_le_bytes(
                buf[offset + 8..offset + 16]
                    .try_into()
                    .expect("slice length must match"),
            );
            let z = f64::from_le_bytes(
                buf[offset + 16..offset + 24]
                    .try_into()
                    .expect("slice length must match"),
            );
            positions.push([x, y, z]);
            offset += 24;
        }
        let mut velocities = Vec::with_capacity(n);
        for _ in 0..n {
            let vx = f64::from_le_bytes(
                buf[offset..offset + 8]
                    .try_into()
                    .expect("slice length must match"),
            );
            let vy = f64::from_le_bytes(
                buf[offset + 8..offset + 16]
                    .try_into()
                    .expect("slice length must match"),
            );
            let vz = f64::from_le_bytes(
                buf[offset + 16..offset + 24]
                    .try_into()
                    .expect("slice length must match"),
            );
            velocities.push([vx, vy, vz]);
            offset += 24;
        }
        (positions, velocities, step, time)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── BinaryHeader tests ────────────────────────────────────────────────

    #[test]
    fn test_header_new_fields() {
        let h = BinaryHeader::new(100, 50, 0.001);
        assert_eq!(h.magic, OXIPHY_MAGIC);
        assert_eq!(h.version, 1);
        assert_eq!(h.n_particles, 100);
        assert_eq!(h.n_frames, 50);
        assert!((h.dt - 0.001).abs() < 1e-15);
    }

    #[test]
    fn test_header_validate_ok() {
        let h = BinaryHeader::new(10, 5, 0.01);
        assert!(h.validate());
    }

    #[test]
    fn test_header_validate_bad_magic() {
        let mut h = BinaryHeader::new(10, 5, 0.01);
        h.magic = [0u8; 4];
        assert!(!h.validate());
    }

    #[test]
    fn test_header_validate_bad_version() {
        let mut h = BinaryHeader::new(10, 5, 0.01);
        h.version = 0;
        assert!(!h.validate());
    }

    #[test]
    fn test_header_write_read_roundtrip() {
        let h = BinaryHeader::new(42, 7, 0.002);
        let mut buf = Vec::new();
        h.write(&mut buf);
        assert_eq!(buf.len(), BinaryHeader::SIZE);
        let h2 = BinaryHeader::read(&buf).unwrap();
        assert_eq!(h, h2);
    }

    #[test]
    fn test_header_read_too_short() {
        let buf = vec![0u8; 10];
        assert!(BinaryHeader::read(&buf).is_none());
    }

    #[test]
    fn test_header_size_constant() {
        let mut buf = Vec::new();
        BinaryHeader::new(1, 1, 1.0).write(&mut buf);
        assert_eq!(buf.len(), BinaryHeader::SIZE);
    }

    // ── ParticleFrame tests ───────────────────────────────────────────────

    #[test]
    fn test_frame_new_sizes() {
        let f = ParticleFrame::new(0, 0.0, 5);
        assert_eq!(f.positions.len(), 5);
        assert_eq!(f.velocities.len(), 5);
        assert_eq!(f.masses.len(), 5);
    }

    #[test]
    fn test_frame_serialized_size() {
        // 8 + 8 + 3*12 + 3*12 + 3*4 = 16 + 36 + 36 + 12 = 100
        assert_eq!(ParticleFrame::serialized_size(3), 100);
    }

    #[test]
    fn test_frame_serialize_deserialize_roundtrip() {
        let mut f = ParticleFrame::new(7, 0.014, 3);
        f.positions[0] = [1.0, 2.0, 3.0];
        f.velocities[1] = [0.5, 0.6, 0.7];
        f.masses[2] = 4.0;
        let mut buf = Vec::new();
        f.serialize(&mut buf);
        let f2 = ParticleFrame::deserialize(&buf, 3).unwrap();
        assert_eq!(f2.step, 7);
        assert!((f2.time - 0.014).abs() < 1e-12);
        assert_eq!(f2.positions[0], [1.0f32, 2.0, 3.0]);
        assert!((f2.velocities[1][1] - 0.6).abs() < 1e-6);
        assert!((f2.masses[2] - 4.0).abs() < 1e-6);
    }

    #[test]
    fn test_frame_deserialize_too_short() {
        let buf = vec![0u8; 10];
        assert!(ParticleFrame::deserialize(&buf, 5).is_none());
    }

    #[test]
    fn test_frame_zero_particles() {
        let f = ParticleFrame::new(0, 0.0, 0);
        let mut buf = Vec::new();
        f.serialize(&mut buf);
        let f2 = ParticleFrame::deserialize(&buf, 0).unwrap();
        assert_eq!(f2.positions.len(), 0);
    }

    // ── BinaryTrajectoryWriter/Reader tests ───────────────────────────────

    #[test]
    fn test_trajectory_write_read_roundtrip() {
        let n_particles = 4u64;
        let dt = 0.001;
        let mut writer = BinaryTrajectoryWriter::new(n_particles, dt);
        for step in 0u64..3 {
            let mut frame = ParticleFrame::new(step, step as f64 * dt, n_particles as usize);
            for i in 0..n_particles as usize {
                frame.positions[i] = [i as f32, step as f32, 0.0];
            }
            writer.write_frame(&frame);
        }
        assert_eq!(writer.frame_count(), 3);
        let data = writer.finalize();

        let mut reader = BinaryTrajectoryReader::new(data).unwrap();
        assert_eq!(reader.n_frames(), 3);
        for step in 0u64..3 {
            let frame = reader.read_frame().unwrap();
            assert_eq!(frame.step, step);
        }
    }

    #[test]
    fn test_trajectory_reader_invalid_magic() {
        let mut buf = vec![0u8; BinaryHeader::SIZE + 100];
        buf[0] = b'X';
        assert!(BinaryTrajectoryReader::new(buf).is_none());
    }

    #[test]
    fn test_trajectory_seek_frame() {
        let n = 2u64;
        let mut writer = BinaryTrajectoryWriter::new(n, 0.01);
        for step in 0u64..5 {
            let frame = ParticleFrame::new(step, step as f64 * 0.01, n as usize);
            writer.write_frame(&frame);
        }
        let data = writer.finalize();
        let mut reader = BinaryTrajectoryReader::new(data).unwrap();
        reader.seek_frame(3);
        let frame = reader.read_frame().unwrap();
        assert_eq!(frame.step, 3);
    }

    #[test]
    fn test_trajectory_read_past_end() {
        let n = 2u64;
        let mut writer = BinaryTrajectoryWriter::new(n, 0.01);
        writer.write_frame(&ParticleFrame::new(0, 0.0, n as usize));
        let data = writer.finalize();
        let mut reader = BinaryTrajectoryReader::new(data).unwrap();
        assert!(reader.read_frame().is_some());
        assert!(reader.read_frame().is_none());
    }

    #[test]
    fn test_trajectory_frame_count() {
        let mut writer = BinaryTrajectoryWriter::new(1, 0.01);
        assert_eq!(writer.frame_count(), 0);
        writer.write_frame(&ParticleFrame::new(0, 0.0, 1));
        assert_eq!(writer.frame_count(), 1);
    }

    // ── DcdWriter/Reader tests ────────────────────────────────────────────

    #[test]
    fn test_dcd_write_read_roundtrip() {
        let n_atoms = 5;
        let mut writer = DcdWriter::new(n_atoms, 0.002);
        let x: Vec<f32> = (0..n_atoms).map(|i| i as f32).collect();
        let y: Vec<f32> = (0..n_atoms).map(|i| i as f32 * 2.0).collect();
        let z: Vec<f32> = (0..n_atoms).map(|i| i as f32 * 3.0).collect();
        writer.write_frame(&x, &y, &z);
        let buf = writer.finalize();

        let reader = DcdReader::parse_header(buf).unwrap();
        assert_eq!(reader.n_atoms, n_atoms);
        assert_eq!(reader.n_frames, 1);
        let (rx, ry, rz) = reader.read_frame(0).unwrap();
        assert!((rx[2] - 2.0).abs() < 1e-5);
        assert!((ry[2] - 4.0).abs() < 1e-5);
        assert!((rz[2] - 6.0).abs() < 1e-5);
    }

    #[test]
    fn test_dcd_multiple_frames() {
        let n_atoms = 3;
        let mut writer = DcdWriter::new(n_atoms, 0.001);
        for frame_i in 0..4 {
            let x: Vec<f32> = vec![frame_i as f32; n_atoms];
            let y: Vec<f32> = vec![0.0; n_atoms];
            let z: Vec<f32> = vec![0.0; n_atoms];
            writer.write_frame(&x, &y, &z);
        }
        let buf = writer.finalize();
        let reader = DcdReader::parse_header(buf).unwrap();
        assert_eq!(reader.n_frames, 4);
        let (rx, _, _) = reader.read_frame(3).unwrap();
        assert!((rx[0] - 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_dcd_out_of_range_frame() {
        let mut writer = DcdWriter::new(2, 0.001);
        writer.write_frame(&[1.0, 2.0], &[0.0, 0.0], &[0.0, 0.0]);
        let buf = writer.finalize();
        let reader = DcdReader::parse_header(buf).unwrap();
        assert!(reader.read_frame(99).is_none());
    }

    #[test]
    fn test_dcd_parse_header_too_short() {
        assert!(DcdReader::parse_header(vec![0u8; 3]).is_none());
    }

    #[test]
    fn test_dcd_zero_frames() {
        let writer = DcdWriter::new(4, 0.001);
        let buf = writer.finalize();
        let reader = DcdReader::parse_header(buf).unwrap();
        assert_eq!(reader.n_frames, 0);
        assert!(reader.read_frame(0).is_none());
    }

    // ── XtcEncoder tests ─────────────────────────────────────────────────

    #[test]
    fn test_xtc_compress_decompress_roundtrip() {
        let positions: Vec<[f32; 3]> =
            vec![[1.0, 2.0, 3.0], [4.5, -1.5, 0.0], [-3.0, 0.001, 100.0]];
        let buf = XtcEncoder::compress_frame(&positions);
        let decoded = XtcEncoder::decompress_frame(&buf, 3).unwrap();
        for (orig, dec) in positions.iter().zip(decoded.iter()) {
            for i in 0..3 {
                assert!(
                    (orig[i] - dec[i]).abs() < 0.002,
                    "mismatch at component {i}"
                );
            }
        }
    }

    #[test]
    fn test_xtc_empty_frame() {
        let buf = XtcEncoder::compress_frame(&[]);
        let decoded = XtcEncoder::decompress_frame(&buf, 0).unwrap();
        assert!(decoded.is_empty());
    }

    #[test]
    fn test_xtc_decompress_too_short() {
        assert!(XtcEncoder::decompress_frame(&[0u8; 5], 3).is_none());
    }

    #[test]
    fn test_xtc_single_particle() {
        let pos = vec![[0.123f32, -0.456, 7.89]];
        let buf = XtcEncoder::compress_frame(&pos);
        let dec = XtcEncoder::decompress_frame(&buf, 1).unwrap();
        assert!((dec[0][0] - 0.123).abs() < 0.002);
        assert!((dec[0][2] - 7.89).abs() < 0.002);
    }

    #[test]
    fn test_xtc_compressed_size() {
        let n = 10;
        let positions = vec![[0.0f32; 3]; n];
        let buf = XtcEncoder::compress_frame(&positions);
        // 4 (n) + 4 (precision) + n*12 = 8 + 120 = 128
        assert_eq!(buf.len(), 8 + n * 12);
    }

    // ── EnergyLogWriter tests ─────────────────────────────────────────────

    #[test]
    fn test_energy_log_empty() {
        let log = EnergyLogWriter::new();
        assert!(log.is_empty());
        let csv = log.to_csv_string();
        assert!(csv.starts_with("step,ke,pe,total_energy"));
    }

    #[test]
    fn test_energy_log_write_step() {
        let mut log = EnergyLogWriter::new();
        log.write_step(0, 10.0, -20.0, 300.0, 101325.0);
        assert_eq!(log.len(), 1);
    }

    #[test]
    fn test_energy_log_csv_total_energy() {
        let mut log = EnergyLogWriter::new();
        log.write_step(1, 5.0, -3.0, 200.0, 1.0);
        let csv = log.to_csv_string();
        assert!(csv.contains("2.000000"), "expected total energy 2.0 in csv");
    }

    #[test]
    fn test_energy_log_multiple_steps() {
        let mut log = EnergyLogWriter::new();
        for i in 0..10u64 {
            log.write_step(i, i as f64, -(i as f64), 300.0, 1.0);
        }
        assert_eq!(log.len(), 10);
        let csv = log.to_csv_string();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 11); // header + 10 data
    }

    #[test]
    fn test_energy_log_not_empty_after_write() {
        let mut log = EnergyLogWriter::new();
        log.write_step(0, 1.0, 2.0, 3.0, 4.0);
        assert!(!log.is_empty());
    }

    // ── CheckpointWriter tests ────────────────────────────────────────────

    #[test]
    fn test_checkpoint_roundtrip() {
        let positions = vec![[1.0f64, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let velocities = vec![[0.1f64, 0.2, 0.3], [0.4, 0.5, 0.6]];
        let step = 42u64;
        let time = 0.042;
        let buf = CheckpointWriter::save_state(&positions, &velocities, step, time);
        let (pos2, vel2, step2, time2) = CheckpointWriter::load_state(&buf);
        assert_eq!(step2, step);
        assert!((time2 - time).abs() < 1e-15);
        assert!((pos2[0][0] - 1.0).abs() < 1e-15);
        assert!((vel2[1][2] - 0.6).abs() < 1e-15);
    }

    #[test]
    fn test_checkpoint_zero_particles() {
        let buf = CheckpointWriter::save_state(&[], &[], 0, 0.0);
        let (pos, vel, step, time) = CheckpointWriter::load_state(&buf);
        assert!(pos.is_empty());
        assert!(vel.is_empty());
        assert_eq!(step, 0);
        assert!((time).abs() < 1e-15);
    }

    #[test]
    fn test_checkpoint_magic_bytes() {
        let buf = CheckpointWriter::save_state(&[], &[], 0, 0.0);
        assert_eq!(&buf[0..4], b"OXCK");
    }

    #[test]
    fn test_checkpoint_step_and_time() {
        let pos = vec![[0.0f64; 3]];
        let vel = vec![[0.0f64; 3]];
        let buf = CheckpointWriter::save_state(&pos, &vel, 9999, 99.99);
        let (_, _, s, t) = CheckpointWriter::load_state(&buf);
        assert_eq!(s, 9999);
        assert!((t - 99.99).abs() < 1e-10);
    }

    #[test]
    fn test_checkpoint_positions_preserved() {
        let pos: Vec<[f64; 3]> = (0..5).map(|i| [i as f64; 3]).collect();
        let vel = vec![[0.0f64; 3]; 5];
        let buf = CheckpointWriter::save_state(&pos, &vel, 0, 0.0);
        let (pos2, _, _, _) = CheckpointWriter::load_state(&buf);
        for i in 0..5 {
            assert!((pos2[i][0] - i as f64).abs() < 1e-15);
        }
    }

    #[test]
    fn test_checkpoint_velocities_preserved() {
        let pos = vec![[0.0f64; 3]; 3];
        let vel: Vec<[f64; 3]> = (0..3).map(|i| [i as f64 * 0.5; 3]).collect();
        let buf = CheckpointWriter::save_state(&pos, &vel, 0, 0.0);
        let (_, vel2, _, _) = CheckpointWriter::load_state(&buf);
        assert!((vel2[2][1] - 1.0).abs() < 1e-15);
    }

    // ── Integration / edge-case tests ─────────────────────────────────────

    #[test]
    fn test_full_pipeline_write_read() {
        // Write a trajectory with 2 particles and 5 frames, then read all back
        let n = 2u64;
        let dt = 0.005;
        let mut writer = BinaryTrajectoryWriter::new(n, dt);
        for step in 0u64..5 {
            let mut frame = ParticleFrame::new(step, step as f64 * dt, n as usize);
            frame.positions[0] = [step as f32, 0.0, 0.0];
            frame.positions[1] = [0.0, step as f32, 0.0];
            frame.masses = vec![1.0, 2.0];
            writer.write_frame(&frame);
        }
        let data = writer.finalize();
        let mut reader = BinaryTrajectoryReader::new(data).unwrap();
        assert_eq!(reader.n_frames(), 5);
        for step in 0u64..5 {
            let frame = reader.read_frame().unwrap();
            assert_eq!(frame.step, step);
            assert!((frame.positions[0][0] - step as f32).abs() < 1e-6);
            assert!((frame.masses[1] - 2.0).abs() < 1e-6);
        }
        assert!(reader.read_frame().is_none());
    }

    #[test]
    fn test_energy_log_csv_format() {
        let mut log = EnergyLogWriter::new();
        log.write_step(0, 1.0, -1.0, 300.0, 1.0);
        log.write_step(1, 2.0, -2.0, 310.0, 1.1);
        let csv = log.to_csv_string();
        let mut lines = csv.lines();
        let header = lines.next().unwrap();
        assert!(header.contains("step"));
        assert!(header.contains("ke"));
        assert!(header.contains("pe"));
        let first_data = lines.next().unwrap();
        assert!(first_data.starts_with("0,"));
    }

    #[test]
    fn test_dcd_n_atoms_in_reader() {
        let n_atoms = 7;
        let writer = DcdWriter::new(n_atoms, 0.002);
        let buf = writer.finalize();
        let reader = DcdReader::parse_header(buf).unwrap();
        assert_eq!(reader.n_atoms, n_atoms);
    }

    #[test]
    fn test_xtc_precision_loss_is_small() {
        let pos = vec![[3.15625f32, 2.71875, -1.40625]];
        let buf = XtcEncoder::compress_frame(&pos);
        let dec = XtcEncoder::decompress_frame(&buf, 1).unwrap();
        for i in 0..3 {
            assert!(
                (pos[0][i] - dec[0][i]).abs() < 0.002,
                "precision too large for component {i}"
            );
        }
    }

    #[test]
    fn test_binary_header_dt_preserved() {
        let h = BinaryHeader::new(10, 10, std::f64::consts::PI);
        let mut buf = Vec::new();
        h.write(&mut buf);
        let h2 = BinaryHeader::read(&buf).unwrap();
        assert!((h2.dt - std::f64::consts::PI).abs() < 1e-14);
    }

    #[test]
    fn test_particle_frame_step_preserved() {
        let f = ParticleFrame::new(12345, 1.23, 1);
        let mut buf = Vec::new();
        f.serialize(&mut buf);
        let f2 = ParticleFrame::deserialize(&buf, 1).unwrap();
        assert_eq!(f2.step, 12345);
    }
}
