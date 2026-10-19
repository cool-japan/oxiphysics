// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simplified binary trajectory format implementations for XTC and DCD.
//!
//! These are in-memory, simplified versions suitable for testing and
//! lightweight use. Real XTC uses lossy compression (XDR); real DCD is
//! a CHARMM/NAMD Fortran binary format. Both are faithfully mimicked
//! structurally but without external dependencies.
//!
//! Also provides:
//! - Frame seeking (by index)
//! - Trajectory metadata
//! - Coordinate compression (quantization)
//! - Box matrix utilities

// ─────────────────────────────────────────────
//  XTC  (simplified, lossless)
// ─────────────────────────────────────────────

/// Magic number written at the beginning of the simplified XTC byte stream.
const XTC_MAGIC: u32 = 0x5854_4300;

/// Smallest possible XTC frame: step (4) + time (4) + atom count (4) + box (9 × 4).
const XTC_MIN_FRAME_BYTES: usize = 48;

/// Bytes of one XTC atom position (3 × f32).
const XTC_ATOM_BYTES: usize = 12;

/// Nanometres to ångström (XTC stores nm, DCD stores Å).
const NM_TO_ANGSTROM: f32 = 10.0;

/// Femtoseconds per picosecond (XTC times are in ps, the DCD timestep is in fs).
const FS_PER_PS: f32 = 1000.0;

/// Bounds-checked little-endian reader over a byte slice.
///
/// A truncated stream is reported as `Err` instead of panicking, and
/// [`remaining`](Self::remaining) lets a reader compare a count taken from the
/// stream with the bytes that can actually follow it *before* it allocates
/// for that count.
struct LeCursor<'a> {
    data: &'a [u8],
    pos: usize,
    /// Format name used as the prefix of error messages (`"XTC"`, `"DCD"`).
    format: &'static str,
}

impl<'a> LeCursor<'a> {
    fn new(data: &'a [u8], format: &'static str) -> Self {
        Self {
            data,
            pos: 0,
            format,
        }
    }

    /// Number of bytes not yet consumed.
    fn remaining(&self) -> usize {
        self.data.len().saturating_sub(self.pos)
    }

    /// Consume the next `len` bytes.
    fn take(&mut self, len: usize) -> Result<&'a [u8], String> {
        let end = self
            .pos
            .checked_add(len)
            .filter(|&end| end <= self.data.len());
        match end.and_then(|end| self.data.get(self.pos..end)) {
            Some(bytes) => {
                self.pos += len;
                Ok(bytes)
            }
            None => Err(format!(
                "{}: unexpected EOF at offset {}",
                self.format, self.pos
            )),
        }
    }

    /// Consume the next `N` bytes as a fixed-size array.
    fn take_array<const N: usize>(&mut self) -> Result<[u8; N], String> {
        let offset = self.pos;
        let bytes = self.take(N)?;
        <[u8; N]>::try_from(bytes)
            .map_err(|_| format!("{}: unexpected EOF at offset {offset}", self.format))
    }

    fn read_u32(&mut self) -> Result<u32, String> {
        self.take_array::<4>().map(u32::from_le_bytes)
    }

    fn read_i32(&mut self) -> Result<i32, String> {
        self.take_array::<4>().map(i32::from_le_bytes)
    }

    fn read_f32(&mut self) -> Result<f32, String> {
        self.take_array::<4>().map(f32::from_le_bytes)
    }
}

/// A single frame of XTC trajectory data.
pub struct XtcFrame {
    /// Simulation step number.
    pub step: i32,
    /// Simulation time (ps).
    pub time: f32,
    /// 3×3 box matrix (nm).
    pub box_matrix: [[f32; 3]; 3],
    /// Atom positions `[x, y, z]` in nm.
    pub positions: Vec<[f32; 3]>,
}

impl XtcFrame {
    /// Number of atoms in this frame.
    pub fn n_atoms(&self) -> usize {
        self.positions.len()
    }

    /// Compute the center of mass of all positions (unweighted).
    pub fn center_of_geometry(&self) -> [f32; 3] {
        if self.positions.is_empty() {
            return [0.0; 3];
        }
        let n = self.positions.len() as f32;
        let mut cx = 0.0f32;
        let mut cy = 0.0f32;
        let mut cz = 0.0f32;
        for p in &self.positions {
            cx += p[0];
            cy += p[1];
            cz += p[2];
        }
        [cx / n, cy / n, cz / n]
    }

    /// Compute the bounding box of all positions: (min, max).
    pub fn bounding_box(&self) -> ([f32; 3], [f32; 3]) {
        if self.positions.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut min = self.positions[0];
        let mut max = self.positions[0];
        for p in &self.positions {
            for i in 0..3 {
                if p[i] < min[i] {
                    min[i] = p[i];
                }
                if p[i] > max[i] {
                    max[i] = p[i];
                }
            }
        }
        (min, max)
    }

    /// Box volume (product of diagonal elements for orthorhombic box).
    pub fn box_volume(&self) -> f32 {
        self.box_matrix[0][0] * self.box_matrix[1][1] * self.box_matrix[2][2]
    }
}

/// In-memory XTC-like writer.
pub struct SimpleXtcWriter {
    /// All frames accumulated so far.
    pub frames: Vec<XtcFrame>,
}

impl SimpleXtcWriter {
    /// Create an empty writer.
    pub fn new() -> Self {
        Self { frames: Vec::new() }
    }

    /// Append a frame.  The box matrix is set to an identity-like 3 nm cube.
    pub fn add_frame(&mut self, step: i32, time: f32, positions: Vec<[f32; 3]>) {
        let box_matrix = [[3.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 3.0]];
        self.frames.push(XtcFrame {
            step,
            time,
            box_matrix,
            positions,
        });
    }

    /// Append a frame with a custom box matrix.
    pub fn add_frame_with_box(
        &mut self,
        step: i32,
        time: f32,
        positions: Vec<[f32; 3]>,
        box_matrix: [[f32; 3]; 3],
    ) {
        self.frames.push(XtcFrame {
            step,
            time,
            box_matrix,
            positions,
        });
    }

    /// Number of frames currently stored.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Serialise to bytes.
    ///
    /// Binary layout:
    /// ```text
    /// [4]  magic  (u32 LE) = 0x58544300
    /// [4]  n_frames (u32 LE)
    /// per frame:
    ///   [4]  step     (i32 LE)
    ///   [4]  time     (f32 LE)
    ///   [4]  n_atoms  (u32 LE)
    ///   [36] box_matrix (9 × f32 LE, row-major)
    ///   [n_atoms*12] positions (n_atoms × 3 × f32 LE)
    /// ```
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();

        buf.extend_from_slice(&XTC_MAGIC.to_le_bytes());
        buf.extend_from_slice(&(self.frames.len() as u32).to_le_bytes());

        for frame in &self.frames {
            buf.extend_from_slice(&frame.step.to_le_bytes());
            buf.extend_from_slice(&frame.time.to_le_bytes());
            buf.extend_from_slice(&(frame.positions.len() as u32).to_le_bytes());

            for row in &frame.box_matrix {
                for &v in row {
                    buf.extend_from_slice(&v.to_le_bytes());
                }
            }

            for pos in &frame.positions {
                for &v in pos {
                    buf.extend_from_slice(&v.to_le_bytes());
                }
            }
        }

        buf
    }
}

impl Default for SimpleXtcWriter {
    fn default() -> Self {
        Self::new()
    }
}

/// In-memory XTC reader.
pub struct SimpleXtcReader;

impl SimpleXtcReader {
    /// Parse a byte slice produced by [`SimpleXtcWriter::to_bytes`].
    ///
    /// The frame count and every per-frame atom count are read from the
    /// stream, which may be corrupt or hostile, so each is checked against the
    /// bytes that remain *before* anything is allocated for it: a stream that
    /// declares more frames or atoms than it can possibly hold is rejected with
    /// an `Err` instead of aborting the process on a multi-gigabyte allocation.
    pub fn from_bytes(data: &[u8]) -> Result<Vec<XtcFrame>, String> {
        let mut cur = LeCursor::new(data, "XTC");

        let magic = cur.read_u32()?;
        if magic != XTC_MAGIC {
            return Err(format!("XTC: bad magic 0x{magic:08X}"));
        }

        let n_frames = cur.read_u32()? as usize;
        if n_frames > cur.remaining() / XTC_MIN_FRAME_BYTES {
            return Err(format!(
                "XTC: header declares {n_frames} frames of at least {XTC_MIN_FRAME_BYTES} bytes but only {} bytes follow",
                cur.remaining()
            ));
        }
        let mut frames = Vec::with_capacity(n_frames);

        for frame_idx in 0..n_frames {
            let step = cur.read_i32()?;
            let time = cur.read_f32()?;
            let n_atoms = cur.read_u32()? as usize;

            let mut box_matrix = [[0f32; 3]; 3];
            for row in &mut box_matrix {
                for v in row.iter_mut() {
                    *v = cur.read_f32()?;
                }
            }

            if n_atoms > cur.remaining() / XTC_ATOM_BYTES {
                return Err(format!(
                    "XTC: frame {frame_idx} declares {n_atoms} atoms of {XTC_ATOM_BYTES} bytes but only {} bytes follow",
                    cur.remaining()
                ));
            }
            let mut positions = Vec::with_capacity(n_atoms);
            for _ in 0..n_atoms {
                let x = cur.read_f32()?;
                let y = cur.read_f32()?;
                let z = cur.read_f32()?;
                positions.push([x, y, z]);
            }

            frames.push(XtcFrame {
                step,
                time,
                box_matrix,
                positions,
            });
        }

        Ok(frames)
    }

    /// Read only frame at a given index (seeks past earlier frames).
    pub fn read_frame_at(data: &[u8], frame_idx: usize) -> Result<XtcFrame, String> {
        let frames = Self::from_bytes(data)?;
        if frame_idx >= frames.len() {
            return Err(format!(
                "frame index {frame_idx} out of range (total {})",
                frames.len()
            ));
        }
        // We have to parse all frames in this simplified format,
        // but we return only the requested one.
        let mut frames = frames;
        Ok(frames.swap_remove(frame_idx.min(frames.len() - 1)))
    }
}

/// Trajectory metadata extracted from an XTC byte stream.
pub struct TrajectoryMetadata {
    /// Number of frames.
    pub n_frames: usize,
    /// Number of atoms (from first frame).
    pub n_atoms: usize,
    /// Time of first frame (ps).
    pub first_time: f32,
    /// Time of last frame (ps).
    pub last_time: f32,
    /// Step of first frame.
    pub first_step: i32,
    /// Step of last frame.
    pub last_step: i32,
}

impl TrajectoryMetadata {
    /// Extract metadata from parsed XTC frames.
    pub fn from_xtc_frames(frames: &[XtcFrame]) -> Option<Self> {
        if frames.is_empty() {
            return None;
        }
        Some(Self {
            n_frames: frames.len(),
            n_atoms: frames[0].positions.len(),
            first_time: frames[0].time,
            last_time: frames[frames.len() - 1].time,
            first_step: frames[0].step,
            last_step: frames[frames.len() - 1].step,
        })
    }

    /// Duration of the trajectory in ps.
    pub fn duration(&self) -> f32 {
        self.last_time - self.first_time
    }

    /// Average time step between frames.
    pub fn avg_dt(&self) -> f32 {
        if self.n_frames <= 1 {
            return 0.0;
        }
        self.duration() / (self.n_frames - 1) as f32
    }
}

// ─────────────────────────────────────────────
//  Coordinate compression (quantization)
// ─────────────────────────────────────────────

/// Simple lossy coordinate compression via fixed-point quantization.
///
/// Positions are quantized to a fixed number of decimal places,
/// reducing precision but allowing more compact storage.
pub struct CoordinateQuantizer {
    /// Number of decimal places to preserve.
    pub precision: u32,
    /// Scale factor = 10^precision.
    pub scale: f32,
}

impl CoordinateQuantizer {
    /// Create a quantizer with the given precision (decimal places).
    pub fn new(precision: u32) -> Self {
        let scale = 10.0_f32.powi(precision as i32);
        Self { precision, scale }
    }

    /// Quantize a position to fixed-point, then round back to f32.
    pub fn quantize(&self, pos: [f32; 3]) -> [f32; 3] {
        [
            (pos[0] * self.scale).round() / self.scale,
            (pos[1] * self.scale).round() / self.scale,
            (pos[2] * self.scale).round() / self.scale,
        ]
    }

    /// Quantize all positions in-place.
    pub fn quantize_frame(&self, positions: &mut [[f32; 3]]) {
        for p in positions.iter_mut() {
            *p = self.quantize(*p);
        }
    }

    /// Maximum quantization error per axis.
    pub fn max_error(&self) -> f32 {
        0.5 / self.scale
    }
}

// ─────────────────────────────────────────────
//  DCD  (simplified Fortran-record binary)
// ─────────────────────────────────────────────

/// Header metadata for a DCD file.
pub struct DcdHeader {
    /// Number of atoms.
    pub n_atoms: u32,
    /// Number of frames.
    pub n_frames: u32,
    /// Integration timestep (fs).
    pub timestep: f32,
    /// Title string (up to 80 characters in the real format).
    pub title: String,
}

impl DcdHeader {
    /// Duration of the trajectory in fs.
    pub fn duration(&self) -> f32 {
        self.timestep * self.n_frames as f32
    }
}

/// A single DCD coordinate frame.  CHARMM stores x, y, z as separate arrays.
pub struct DcdFrame {
    /// X coordinates (Å).
    pub x: Vec<f32>,
    /// Y coordinates (Å).
    pub y: Vec<f32>,
    /// Z coordinates (Å).
    pub z: Vec<f32>,
}

impl DcdFrame {
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.x.len()
    }

    /// Get position of atom i as \[x, y, z\].
    pub fn position(&self, i: usize) -> [f32; 3] {
        [self.x[i], self.y[i], self.z[i]]
    }

    /// Convert to interleaved positions.
    pub fn to_interleaved(&self) -> Vec<[f32; 3]> {
        (0..self.n_atoms()).map(|i| self.position(i)).collect()
    }

    /// Center of geometry (unweighted).
    pub fn center_of_geometry(&self) -> [f32; 3] {
        if self.x.is_empty() {
            return [0.0; 3];
        }
        let n = self.x.len() as f32;
        let cx: f32 = self.x.iter().sum::<f32>() / n;
        let cy: f32 = self.y.iter().sum::<f32>() / n;
        let cz: f32 = self.z.iter().sum::<f32>() / n;
        [cx, cy, cz]
    }
}

/// In-memory DCD writer.
pub struct SimpleDcdWriter {
    /// File header.
    pub header: DcdHeader,
    /// Accumulated frames.
    pub frames: Vec<DcdFrame>,
}

impl SimpleDcdWriter {
    /// Create a new writer.  `n_atoms` is fixed for the entire trajectory.
    pub fn new(n_atoms: u32, timestep: f32) -> Self {
        let header = DcdHeader {
            n_atoms,
            n_frames: 0,
            timestep,
            title: String::from("OxiPhysics simplified DCD"),
        };
        Self {
            header,
            frames: Vec::new(),
        }
    }

    /// Create a writer with a custom title.
    pub fn with_title(n_atoms: u32, timestep: f32, title: &str) -> Self {
        let header = DcdHeader {
            n_atoms,
            n_frames: 0,
            timestep,
            title: title.to_string(),
        };
        Self {
            header,
            frames: Vec::new(),
        }
    }

    /// Append a frame from interleaved `[x, y, z]` positions.
    pub fn add_frame(&mut self, positions: &[[f32; 3]]) {
        let x: Vec<f32> = positions.iter().map(|p| p[0]).collect();
        let y: Vec<f32> = positions.iter().map(|p| p[1]).collect();
        let z: Vec<f32> = positions.iter().map(|p| p[2]).collect();
        self.frames.push(DcdFrame { x, y, z });
        self.header.n_frames = self.frames.len() as u32;
    }

    /// Append a frame from separate x, y, z arrays.
    pub fn add_frame_xyz(&mut self, x: Vec<f32>, y: Vec<f32>, z: Vec<f32>) {
        self.frames.push(DcdFrame { x, y, z });
        self.header.n_frames = self.frames.len() as u32;
    }

    /// Number of frames.
    pub fn frame_count(&self) -> usize {
        self.frames.len()
    }

    /// Serialise to a simplified DCD byte stream.
    ///
    /// Each Fortran record is wrapped with a 4-byte integer giving the
    /// byte length of the payload (written before **and** after the block).
    ///
    /// Layout:
    /// ```text
    /// HEADER BLOCK:
    ///   [4] block_len
    ///   [4] "CORD" signature
    ///   [4] n_frames  (u32 LE)
    ///   [4] n_atoms   (u32 LE)
    ///   [4] timestep  (f32 LE)
    ///   [80] title (padded with spaces)
    ///   [4] block_len  (repeated)
    ///
    /// per frame:
    ///   X block: [4] len_x | x[0..n] as f32 LE | [4] len_x
    ///   Y block: [4] len_y | y[0..n] as f32 LE | [4] len_y
    ///   Z block: [4] len_z | z[0..n] as f32 LE | [4] len_z
    /// ```
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();

        // ── header block ──────────────────────────────────────────────────
        // payload: "CORD"(4) + n_frames(4) + n_atoms(4) + timestep(4) + title(80) = 96 bytes
        let hdr_payload_len: u32 = 4 + 4 + 4 + 4 + 80;
        buf.extend_from_slice(&hdr_payload_len.to_le_bytes());
        buf.extend_from_slice(b"CORD");
        buf.extend_from_slice(&self.header.n_frames.to_le_bytes());
        buf.extend_from_slice(&self.header.n_atoms.to_le_bytes());
        buf.extend_from_slice(&self.header.timestep.to_le_bytes());

        // title: exactly 80 bytes, space-padded / truncated
        let mut title_bytes = [b' '; 80];
        let src = self.header.title.as_bytes();
        let copy_len = src.len().min(80);
        title_bytes[..copy_len].copy_from_slice(&src[..copy_len]);
        buf.extend_from_slice(&title_bytes);

        buf.extend_from_slice(&hdr_payload_len.to_le_bytes()); // closing record length

        // ── coordinate frames ─────────────────────────────────────────────
        let coord_block_len = self.header.n_atoms * 4; // n_atoms × sizeof(f32)

        for frame in &self.frames {
            for coords in [&frame.x, &frame.y, &frame.z] {
                buf.extend_from_slice(&coord_block_len.to_le_bytes());
                for &v in coords.iter() {
                    buf.extend_from_slice(&v.to_le_bytes());
                }
                buf.extend_from_slice(&coord_block_len.to_le_bytes());
            }
        }

        buf
    }
}

/// In-memory DCD reader.
pub struct SimpleDcdReader;

impl SimpleDcdReader {
    /// Parse a byte slice produced by [`SimpleDcdWriter::to_bytes`].
    ///
    /// The frame and atom counts come from the header of a stream that may be
    /// corrupt or hostile. Every frame has a fixed size (three records of a
    /// 4-byte length, `4 * n_atoms` coordinate bytes and a 4-byte length), so
    /// the exact size the header declares is checked against the bytes that
    /// follow it *before* any allocation is sized from the counts.
    pub fn from_bytes(data: &[u8]) -> Result<(DcdHeader, Vec<DcdFrame>), String> {
        let mut cur = LeCursor::new(data, "DCD");

        // ── header block ──────────────────────────────────────────────────
        let hdr_len = cur.read_u32()? as usize;
        if hdr_len < 4 + 4 + 4 + 4 + 80 {
            return Err(format!("DCD: header block too short ({hdr_len})"));
        }

        if cur.take(4)? != b"CORD" {
            return Err("DCD: expected 'CORD' signature".to_string());
        }

        let n_frames = cur.read_u32()?;
        let n_atoms = cur.read_u32()?;
        let timestep = cur.read_f32()?;

        let title = String::from_utf8_lossy(cur.take(80)?)
            .trim_end()
            .to_string();

        // closing record length
        let hdr_len2 = cur.read_u32()?;
        if hdr_len2 as usize != hdr_len {
            return Err(format!(
                "DCD: header closing length mismatch ({hdr_len2} vs {hdr_len})"
            ));
        }

        let header = DcdHeader {
            n_atoms,
            n_frames,
            timestep,
            title,
        };
        let n = n_atoms as usize;

        // ── frames ────────────────────────────────────────────────────────
        // A frame is three records of (4-byte length, 4 * n bytes, 4-byte
        // length), so a stream that is shorter than the declared frames
        // cannot be valid. u128 arithmetic cannot overflow for u32 counts.
        let frame_bytes = 3 * (8 + 4 * u128::from(n_atoms));
        let declared_bytes = u128::from(n_frames) * frame_bytes;
        let available_bytes = cur.remaining() as u128;
        if declared_bytes > available_bytes {
            return Err(format!(
                "DCD: header declares {n_frames} frames of {n_atoms} atoms ({declared_bytes} bytes) but only {available_bytes} bytes follow"
            ));
        }
        let mut frames = Vec::with_capacity(n_frames as usize);

        for frame_idx in 0..n_frames {
            let mut xyz: [Vec<f32>; 3] = [
                Vec::with_capacity(n),
                Vec::with_capacity(n),
                Vec::with_capacity(n),
            ];

            for (dim, xyz_dim) in xyz.iter_mut().enumerate() {
                let block_len = cur.read_u32()? as usize;
                let expected = n * 4;
                if block_len != expected {
                    return Err(format!(
                        "DCD: frame {frame_idx} dim {dim}: block length {block_len} != {expected}"
                    ));
                }
                for _ in 0..n {
                    xyz_dim.push(cur.read_f32()?);
                }
                let block_len2 = cur.read_u32()? as usize;
                if block_len2 != block_len {
                    return Err(format!(
                        "DCD: frame {frame_idx} dim {dim}: closing block length mismatch"
                    ));
                }
            }

            let [x, y, z] = xyz;
            frames.push(DcdFrame { x, y, z });
        }

        Ok((header, frames))
    }

    /// Read metadata from a DCD byte stream without parsing all frames.
    pub fn read_metadata(data: &[u8]) -> Result<DcdHeader, String> {
        let (header, _) = Self::from_bytes(data)?;
        Ok(header)
    }
}

// ─────────────────────────────────────────────
//  DCD trajectory writing helpers
// ─────────────────────────────────────────────

/// Write a sequence of interleaved position frames to a DCD byte stream.
///
/// Each element in `frames` is a slice of `[f32; 3]` positions for one frame.
/// Returns the complete DCD byte stream.
pub fn write_dcd_trajectory(n_atoms: u32, timestep: f32, frames: &[Vec<[f32; 3]>]) -> Vec<u8> {
    let mut writer = SimpleDcdWriter::new(n_atoms, timestep);
    for frame in frames {
        writer.add_frame(frame);
    }
    writer.to_bytes()
}

/// Append a single frame to an existing DCD byte stream.
///
/// This reads the existing stream, appends the new frame, and returns the
/// updated byte stream.
pub fn append_dcd_frame(existing: &[u8], new_positions: &[[f32; 3]]) -> Result<Vec<u8>, String> {
    let (header, mut frames) = SimpleDcdReader::from_bytes(existing)?;
    let mut writer = SimpleDcdWriter::with_title(header.n_atoms, header.timestep, &header.title);
    for frame in &frames {
        writer.add_frame_xyz(frame.x.clone(), frame.y.clone(), frame.z.clone());
    }
    // Add new frame
    let new_frame: Vec<[f32; 3]> = new_positions.to_vec();
    writer.add_frame(&new_frame);
    // Suppress warning about unused mutable variable
    frames.clear();
    Ok(writer.to_bytes())
}

// ─────────────────────────────────────────────
//  XTC precision control
// ─────────────────────────────────────────────

/// Apply a specific precision level to all frames in a writer.
///
/// Precision is given as the number of decimal places to keep (1-6).
/// Higher precision retains more spatial detail but produces larger files.
pub fn apply_xtc_precision(writer: &mut SimpleXtcWriter, precision: u32) {
    let q = CoordinateQuantizer::new(precision);
    for frame in writer.frames.iter_mut() {
        q.quantize_frame(&mut frame.positions);
    }
}

/// Round-trip an XTC byte stream through a given precision level.
///
/// Returns the re-encoded bytes at the given coordinate precision.
pub fn xtc_recompress(data: &[u8], precision: u32) -> Result<Vec<u8>, String> {
    let frames = SimpleXtcReader::from_bytes(data)?;
    let q = CoordinateQuantizer::new(precision);
    let mut writer = SimpleXtcWriter::new();
    for mut frame in frames {
        q.quantize_frame(&mut frame.positions);
        writer.add_frame_with_box(frame.step, frame.time, frame.positions, frame.box_matrix);
    }
    Ok(writer.to_bytes())
}

/// Compute the root-mean-square deviation (RMSD) between two sets of positions.
pub fn compute_rmsd(a: &[[f32; 3]], b: &[[f32; 3]]) -> f32 {
    if a.is_empty() || a.len() != b.len() {
        return 0.0;
    }
    let n = a.len() as f32;
    let sum_sq: f32 = a
        .iter()
        .zip(b.iter())
        .map(|(ai, bi)| {
            let dx = ai[0] - bi[0];
            let dy = ai[1] - bi[1];
            let dz = ai[2] - bi[2];
            dx * dx + dy * dy + dz * dz
        })
        .sum();
    (sum_sq / n).sqrt()
}

// ─────────────────────────────────────────────
//  Trajectory format conversion (XTC ↔ DCD)
// ─────────────────────────────────────────────

/// Convert a simplified XTC byte stream to a DCD byte stream.
///
/// The two formats use different units, so the conversion rescales:
///
/// * positions: nanometres (XTC) → ångström (DCD), multiplied by 10;
/// * timestep: the average inter-frame interval in picoseconds (XTC) →
///   femtoseconds ([`DcdHeader::timestep`]), multiplied by 1000.
///
/// A DCD stream carries a single timestep, no start offset and no step
/// numbers, so frame times and step numbers survive an XTC → DCD → XTC round
/// trip only when the frames are uniformly spaced and start at time 0 (the
/// step number of frame `i` becomes `i`). For unevenly spaced frames the
/// average interval is stored.
///
/// The atom count is taken from the first XTC frame and a DCD stream has a
/// constant atom count, so an `Err` is returned if any frame has a different
/// number of atoms.
pub fn xtc_to_dcd(xtc_data: &[u8]) -> Result<Vec<u8>, String> {
    let frames = SimpleXtcReader::from_bytes(xtc_data)?;
    let Some(meta) = TrajectoryMetadata::from_xtc_frames(&frames) else {
        return Ok(SimpleDcdWriter::new(0, 0.0).to_bytes());
    };
    let n_atoms = meta.n_atoms;
    for (idx, frame) in frames.iter().enumerate() {
        if frame.n_atoms() != n_atoms {
            return Err(format!(
                "XTC frame {idx} has {} atoms but frame 0 has {n_atoms}; DCD requires a constant atom count",
                frame.n_atoms()
            ));
        }
    }
    let n_atoms_u32 = u32::try_from(n_atoms)
        .map_err(|_| format!("XTC atom count {n_atoms} does not fit a DCD header"))?;

    let timestep_fs = meta.avg_dt() * FS_PER_PS;
    let mut writer = SimpleDcdWriter::with_title(n_atoms_u32, timestep_fs, "Converted from XTC");
    for frame in &frames {
        let x: Vec<f32> = frame
            .positions
            .iter()
            .map(|p| p[0] * NM_TO_ANGSTROM)
            .collect();
        let y: Vec<f32> = frame
            .positions
            .iter()
            .map(|p| p[1] * NM_TO_ANGSTROM)
            .collect();
        let z: Vec<f32> = frame
            .positions
            .iter()
            .map(|p| p[2] * NM_TO_ANGSTROM)
            .collect();
        writer.add_frame_xyz(x, y, z);
    }
    Ok(writer.to_bytes())
}

/// Convert a simplified DCD byte stream to an XTC byte stream.
///
/// The inverse of [`xtc_to_dcd`]:
///
/// * positions: ångström (DCD) → nanometres (XTC), divided by 10;
/// * time: frame `i` gets the step number `i` and the time
///   `i * timestep`, with the DCD timestep converted from femtoseconds to
///   picoseconds (divided by 1000).
pub fn dcd_to_xtc(dcd_data: &[u8]) -> Result<Vec<u8>, String> {
    let (header, frames) = SimpleDcdReader::from_bytes(dcd_data)?;
    let timestep_ps = header.timestep / FS_PER_PS;
    let mut writer = SimpleXtcWriter::new();
    for (idx, frame) in frames.iter().enumerate() {
        let positions_nm: Vec<[f32; 3]> = (0..frame.n_atoms())
            .map(|i| {
                [
                    frame.x[i] / NM_TO_ANGSTROM,
                    frame.y[i] / NM_TO_ANGSTROM,
                    frame.z[i] / NM_TO_ANGSTROM,
                ]
            })
            .collect();
        let step = i32::try_from(idx)
            .map_err(|_| format!("DCD frame index {idx} does not fit an XTC step number"))?;
        let time_ps = timestep_ps * idx as f32;
        writer.add_frame(step, time_ps, positions_nm);
    }
    Ok(writer.to_bytes())
}

// ─────────────────────────────────────────────
//  Frame time extraction
// ─────────────────────────────────────────────

/// Extract the simulation time (in ps) for every frame in an XTC byte stream.
pub fn xtc_extract_frame_times(data: &[u8]) -> Result<Vec<f32>, String> {
    let frames = SimpleXtcReader::from_bytes(data)?;
    Ok(frames.iter().map(|f| f.time).collect())
}

/// Extract the simulation step number for every frame in an XTC byte stream.
pub fn xtc_extract_frame_steps(data: &[u8]) -> Result<Vec<i32>, String> {
    let frames = SimpleXtcReader::from_bytes(data)?;
    Ok(frames.iter().map(|f| f.step).collect())
}

/// Compute inter-frame time deltas (ps) from an XTC byte stream.
pub fn xtc_frame_time_deltas(data: &[u8]) -> Result<Vec<f32>, String> {
    let times = xtc_extract_frame_times(data)?;
    if times.len() < 2 {
        return Ok(Vec::new());
    }
    Ok(times.windows(2).map(|w| w[1] - w[0]).collect())
}

// ─────────────────────────────────────────────
//  Velocity data support
// ─────────────────────────────────────────────

/// An XTC-like frame augmented with per-atom velocity data.
pub struct XtcFrameWithVelocity {
    /// Underlying position frame.
    pub frame: XtcFrame,
    /// Per-atom velocities `[vx, vy, vz]` in nm/ps.
    pub velocities: Vec<[f32; 3]>,
}

impl XtcFrameWithVelocity {
    /// Create from a position-only frame and separately computed velocities.
    pub fn new(frame: XtcFrame, velocities: Vec<[f32; 3]>) -> Self {
        Self { frame, velocities }
    }

    /// Kinetic energy (sum of 0.5 * v^2 for each atom, unweighted).
    pub fn kinetic_energy_proxy(&self) -> f32 {
        self.velocities
            .iter()
            .map(|v| 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }

    /// RMS speed (root-mean-square over all atoms).
    pub fn rms_speed(&self) -> f32 {
        if self.velocities.is_empty() {
            return 0.0;
        }
        let n = self.velocities.len() as f32;
        let sum_sq: f32 = self
            .velocities
            .iter()
            .map(|v| v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
            .sum();
        (sum_sq / n).sqrt()
    }

    /// Maximum speed among all atoms.
    pub fn max_speed(&self) -> f32 {
        self.velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0f32, f32::max)
    }
}

/// A writer for XTC-like trajectories that includes velocity channels.
///
/// The binary layout appends a velocity block after each position block,
/// using the same length-prefix convention as DCD.
pub struct XtcVelocityWriter {
    /// Underlying XTC writer for positions.
    pub positions: SimpleXtcWriter,
    /// Velocity data per frame.
    pub velocities: Vec<Vec<[f32; 3]>>,
}

impl XtcVelocityWriter {
    /// Create an empty velocity writer.
    pub fn new() -> Self {
        Self {
            positions: SimpleXtcWriter::new(),
            velocities: Vec::new(),
        }
    }

    /// Add a frame with both positions and velocities.
    pub fn add_frame(
        &mut self,
        step: i32,
        time: f32,
        positions: Vec<[f32; 3]>,
        velocities: Vec<[f32; 3]>,
    ) {
        self.positions.add_frame(step, time, positions);
        self.velocities.push(velocities);
    }

    /// Number of frames stored.
    pub fn frame_count(&self) -> usize {
        self.positions.frame_count()
    }

    /// Serialize positions and velocities to bytes.
    ///
    /// Layout: XTC bytes + velocity magic (4) + velocity payload.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = self.positions.to_bytes();
        // Velocity block marker
        buf.extend_from_slice(b"OXVL");
        buf.extend_from_slice(&(self.velocities.len() as u32).to_le_bytes());
        for frame_vels in &self.velocities {
            buf.extend_from_slice(&(frame_vels.len() as u32).to_le_bytes());
            for v in frame_vels {
                for &c in v {
                    buf.extend_from_slice(&c.to_le_bytes());
                }
            }
        }
        buf
    }

    /// Compute kinetic energy proxy for a frame.
    pub fn frame_kinetic_energy(&self, frame_idx: usize) -> f32 {
        if frame_idx >= self.velocities.len() {
            return 0.0;
        }
        self.velocities[frame_idx]
            .iter()
            .map(|v| 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
            .sum()
    }
}

impl Default for XtcVelocityWriter {
    fn default() -> Self {
        Self::new()
    }
}

/// Compute finite-difference velocities between consecutive XTC frames.
///
/// Returns one fewer velocity frame than position frames.
/// `dt_ps` is the time step in picoseconds between consecutive frames.
pub fn compute_finite_difference_velocities(frames: &[XtcFrame], dt_ps: f32) -> Vec<Vec<[f32; 3]>> {
    if frames.len() < 2 || dt_ps < 1e-30 {
        return Vec::new();
    }
    frames
        .windows(2)
        .map(|w| {
            let n = w[0].positions.len().min(w[1].positions.len());
            (0..n)
                .map(|i| {
                    let p0 = w[0].positions[i];
                    let p1 = w[1].positions[i];
                    [
                        (p1[0] - p0[0]) / dt_ps,
                        (p1[1] - p0[1]) / dt_ps,
                        (p1[2] - p0[2]) / dt_ps,
                    ]
                })
                .collect()
        })
        .collect()
}

// ─────────────────────────────────────────────
//  Tests
// ─────────────────────────────────────────────

// Untrusted frame/atom counts in DCD and XTC headers (F2416).
#[cfg(test)]
mod dcd_header_tests;
#[cfg(test)]
mod xtc_header_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: f32, b: f32) -> bool {
        (a - b).abs() < 1e-5
    }

    // ── XTC tests ────────────────────────────────────────────────────────

    #[test]
    fn xtc_empty_writer() {
        let w = SimpleXtcWriter::new();
        assert_eq!(w.frame_count(), 0);
        let bytes = w.to_bytes();
        // magic (4) + n_frames (4) = 8 bytes
        assert_eq!(bytes.len(), 8);
    }

    #[test]
    fn xtc_round_trip_3frames_5atoms() {
        let mut writer = SimpleXtcWriter::new();

        let pos0: Vec<[f32; 3]> = (0..5)
            .map(|i| [i as f32 * 0.1, i as f32 * 0.2, i as f32 * 0.3])
            .collect();
        let pos1: Vec<[f32; 3]> = (0..5)
            .map(|i| [i as f32 * 0.4, i as f32 * 0.5, i as f32 * 0.6])
            .collect();
        let pos2: Vec<[f32; 3]> = (0..5)
            .map(|i| [i as f32 * 0.7, i as f32 * 0.8, i as f32 * 0.9])
            .collect();

        writer.add_frame(0, 0.0, pos0.clone());
        writer.add_frame(1, 0.5, pos1.clone());
        writer.add_frame(2, 1.0, pos2.clone());

        assert_eq!(writer.frame_count(), 3);

        let bytes = writer.to_bytes();
        let frames = SimpleXtcReader::from_bytes(&bytes).expect("XTC parse failed");

        assert_eq!(frames.len(), 3);

        for (f, original) in frames.iter().zip([&pos0, &pos1, &pos2]) {
            assert_eq!(f.positions.len(), 5);
            for (got, exp) in f.positions.iter().zip(original.iter()) {
                for k in 0..3 {
                    assert!(
                        approx_eq(got[k], exp[k]),
                        "mismatch: got {} exp {}",
                        got[k],
                        exp[k]
                    );
                }
            }
        }

        assert_eq!(frames[0].step, 0);
        assert!(approx_eq(frames[1].time, 0.5));
        assert_eq!(frames[2].step, 2);
    }

    #[test]
    fn xtc_frame_count_check() {
        let mut w = SimpleXtcWriter::new();
        for i in 0..10 {
            w.add_frame(i, i as f32 * 0.1, vec![[0.0; 3]]);
        }
        assert_eq!(w.frame_count(), 10);
    }

    // ── XTC frame properties ─────────────────────────────────────────────

    #[test]
    fn xtc_frame_center_of_geometry() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        };
        let cog = frame.center_of_geometry();
        assert!(approx_eq(cog[0], 2.0 / 3.0));
        assert!(approx_eq(cog[1], 2.0 / 3.0));
        assert!(approx_eq(cog[2], 0.0));
    }

    #[test]
    fn xtc_frame_bounding_box() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            positions: vec![[-1.0, 2.0, 0.0], [3.0, -1.0, 5.0]],
        };
        let (min, max) = frame.bounding_box();
        assert!(approx_eq(min[0], -1.0));
        assert!(approx_eq(max[0], 3.0));
        assert!(approx_eq(min[1], -1.0));
        assert!(approx_eq(max[1], 2.0));
    }

    #[test]
    fn xtc_frame_box_volume() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[3.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]],
            positions: vec![],
        };
        assert!(approx_eq(frame.box_volume(), 60.0));
    }

    // ── Frame seeking ────────────────────────────────────────────────────

    #[test]
    fn xtc_read_frame_at() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[1.0, 0.0, 0.0]]);
        writer.add_frame(1, 0.5, vec![[2.0, 0.0, 0.0]]);
        writer.add_frame(2, 1.0, vec![[3.0, 0.0, 0.0]]);

        let bytes = writer.to_bytes();
        let frame = SimpleXtcReader::read_frame_at(&bytes, 1).unwrap();
        assert_eq!(frame.step, 1);
    }

    #[test]
    fn xtc_read_frame_at_out_of_range() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[0.0; 3]]);
        let bytes = writer.to_bytes();
        assert!(SimpleXtcReader::read_frame_at(&bytes, 5).is_err());
    }

    // ── Trajectory metadata ──────────────────────────────────────────────

    #[test]
    fn trajectory_metadata_basic() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[0.0; 3]; 10]);
        writer.add_frame(100, 10.0, vec![[0.0; 3]; 10]);
        writer.add_frame(200, 20.0, vec![[0.0; 3]; 10]);

        let bytes = writer.to_bytes();
        let frames = SimpleXtcReader::from_bytes(&bytes).unwrap();
        let meta = TrajectoryMetadata::from_xtc_frames(&frames).unwrap();

        assert_eq!(meta.n_frames, 3);
        assert_eq!(meta.n_atoms, 10);
        assert!(approx_eq(meta.first_time, 0.0));
        assert!(approx_eq(meta.last_time, 20.0));
        assert!(approx_eq(meta.duration(), 20.0));
        assert!(approx_eq(meta.avg_dt(), 10.0));
    }

    #[test]
    fn trajectory_metadata_empty() {
        let frames: Vec<XtcFrame> = vec![];
        assert!(TrajectoryMetadata::from_xtc_frames(&frames).is_none());
    }

    // ── Coordinate quantizer ─────────────────────────────────────────────

    #[test]
    fn quantizer_precision_3() {
        let q = CoordinateQuantizer::new(3);
        let pos = [1.23456, 7.89012, -0.00050];
        let qp = q.quantize(pos);
        assert!(approx_eq(qp[0], 1.235));
        assert!(approx_eq(qp[1], 7.890));
        // -0.0005 rounds to -0.001 or 0.0 depending on float rounding
        assert!(qp[2].abs() < 0.001 + 1e-6);
    }

    #[test]
    fn quantizer_max_error() {
        let q = CoordinateQuantizer::new(3);
        assert!(approx_eq(q.max_error(), 0.0005));
    }

    #[test]
    fn quantizer_frame() {
        let q = CoordinateQuantizer::new(2);
        let mut positions = vec![[1.234, 5.678, 9.012]];
        q.quantize_frame(&mut positions);
        assert!(approx_eq(positions[0][0], 1.23));
        assert!(approx_eq(positions[0][1], 5.68));
    }

    // ── Custom box matrix ────────────────────────────────────────────────

    #[test]
    fn xtc_custom_box_matrix() {
        let mut writer = SimpleXtcWriter::new();
        let custom_box = [[5.0, 0.0, 0.0], [0.0, 6.0, 0.0], [0.0, 0.0, 7.0]];
        writer.add_frame_with_box(0, 0.0, vec![[0.0; 3]], custom_box);

        let bytes = writer.to_bytes();
        let frames = SimpleXtcReader::from_bytes(&bytes).unwrap();
        assert!(approx_eq(frames[0].box_matrix[0][0], 5.0));
        assert!(approx_eq(frames[0].box_matrix[1][1], 6.0));
        assert!(approx_eq(frames[0].box_matrix[2][2], 7.0));
    }

    // ── DCD tests ────────────────────────────────────────────────────────

    #[test]
    fn dcd_empty_writer() {
        let w = SimpleDcdWriter::new(4, 2.0);
        assert_eq!(w.frame_count(), 0);
        // Ensure to_bytes doesn't panic and produces valid header
        let bytes = w.to_bytes();
        assert!(!bytes.is_empty());
    }

    #[test]
    fn dcd_round_trip_2frames_4atoms() {
        let mut writer = SimpleDcdWriter::new(4, 2.0);

        let frame0: Vec<[f32; 3]> = vec![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ];
        let frame1: Vec<[f32; 3]> = vec![
            [0.1, 0.2, 0.3],
            [0.4, 0.5, 0.6],
            [0.7, 0.8, 0.9],
            [1.0, 1.1, 1.2],
        ];

        writer.add_frame(&frame0);
        writer.add_frame(&frame1);

        assert_eq!(writer.frame_count(), 2);
        assert_eq!(writer.header.n_frames, 2);

        let bytes = writer.to_bytes();
        let (header, frames) = SimpleDcdReader::from_bytes(&bytes).expect("DCD parse failed");

        assert_eq!(header.n_atoms, 4);
        assert_eq!(header.n_frames, 2);
        assert!(approx_eq(header.timestep, 2.0));
        assert_eq!(frames.len(), 2);

        // Verify X coordinates of frame 0
        let expected_x0: Vec<f32> = frame0.iter().map(|p| p[0]).collect();
        for (got, &exp) in frames[0].x.iter().zip(expected_x0.iter()) {
            assert!(approx_eq(*got, exp), "x mismatch: got {got} exp {exp}");
        }

        // Verify all coords of frame 1
        for (i, pos) in frame1.iter().enumerate() {
            assert!(approx_eq(frames[1].x[i], pos[0]));
            assert!(approx_eq(frames[1].y[i], pos[1]));
            assert!(approx_eq(frames[1].z[i], pos[2]));
        }
    }

    #[test]
    fn dcd_frame_count_check() {
        let mut w = SimpleDcdWriter::new(2, 1.0);
        for _ in 0..7 {
            w.add_frame(&[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        }
        assert_eq!(w.frame_count(), 7);
    }

    // ── DCD frame methods ────────────────────────────────────────────────

    #[test]
    fn dcd_frame_to_interleaved() {
        let frame = DcdFrame {
            x: vec![1.0, 4.0],
            y: vec![2.0, 5.0],
            z: vec![3.0, 6.0],
        };
        let interleaved = frame.to_interleaved();
        assert_eq!(interleaved.len(), 2);
        assert!(approx_eq(interleaved[0][0], 1.0));
        assert!(approx_eq(interleaved[1][2], 6.0));
    }

    #[test]
    fn dcd_frame_center_of_geometry() {
        let frame = DcdFrame {
            x: vec![0.0, 2.0],
            y: vec![0.0, 4.0],
            z: vec![0.0, 6.0],
        };
        let cog = frame.center_of_geometry();
        assert!(approx_eq(cog[0], 1.0));
        assert!(approx_eq(cog[1], 2.0));
        assert!(approx_eq(cog[2], 3.0));
    }

    // ── DCD custom title ─────────────────────────────────────────────────

    #[test]
    fn dcd_custom_title() {
        let mut writer = SimpleDcdWriter::with_title(2, 1.0, "Custom Title");
        writer.add_frame(&[[0.0; 3], [1.0; 3]]);

        let bytes = writer.to_bytes();
        let (header, _) = SimpleDcdReader::from_bytes(&bytes).unwrap();
        assert!(header.title.contains("Custom Title"));
    }

    // ── DCD add_frame_xyz ────────────────────────────────────────────────

    #[test]
    fn dcd_add_frame_xyz() {
        let mut writer = SimpleDcdWriter::new(3, 1.0);
        writer.add_frame_xyz(
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        );
        assert_eq!(writer.frame_count(), 1);

        let bytes = writer.to_bytes();
        let (_, frames) = SimpleDcdReader::from_bytes(&bytes).unwrap();
        assert_eq!(frames[0].x, vec![1.0, 2.0, 3.0]);
    }

    // ── DCD header duration ──────────────────────────────────────────────

    #[test]
    fn dcd_header_duration() {
        let header = DcdHeader {
            n_atoms: 10,
            n_frames: 100,
            timestep: 2.0,
            title: String::new(),
        };
        assert!(approx_eq(header.duration(), 200.0));
    }

    // ── DCD read_metadata ────────────────────────────────────────────────

    #[test]
    fn dcd_read_metadata() {
        let mut writer = SimpleDcdWriter::new(5, 3.0);
        writer.add_frame(&[[0.0; 3]; 5]);
        writer.add_frame(&[[1.0; 3]; 5]);

        let bytes = writer.to_bytes();
        let header = SimpleDcdReader::read_metadata(&bytes).unwrap();
        assert_eq!(header.n_atoms, 5);
        assert_eq!(header.n_frames, 2);
        assert!(approx_eq(header.timestep, 3.0));
    }

    // ── DCD trajectory writing helpers ───────────────────────────────────

    #[test]
    fn write_dcd_trajectory_basic() {
        let f0: Vec<[f32; 3]> = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let f1: Vec<[f32; 3]> = vec![[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]];
        let bytes = write_dcd_trajectory(2, 2.0, &[f0, f1]);
        let (header, frames) = SimpleDcdReader::from_bytes(&bytes).unwrap();
        assert_eq!(header.n_atoms, 2);
        assert_eq!(header.n_frames, 2);
        assert!(approx_eq(frames[0].x[0], 1.0));
        assert!(approx_eq(frames[1].z[1], 12.0));
    }

    #[test]
    fn write_dcd_trajectory_empty() {
        let bytes = write_dcd_trajectory(3, 1.0, &[]);
        let (header, frames) = SimpleDcdReader::from_bytes(&bytes).unwrap();
        assert_eq!(header.n_frames, 0);
        assert!(frames.is_empty());
    }

    #[test]
    fn append_dcd_frame_basic() {
        let mut writer = SimpleDcdWriter::new(2, 1.0);
        writer.add_frame(&[[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]);
        let original = writer.to_bytes();

        let new_frame: Vec<[f32; 3]> = vec![[2.0, 2.0, 2.0], [3.0, 3.0, 3.0]];
        let updated = append_dcd_frame(&original, &new_frame).unwrap();

        let (header, frames) = SimpleDcdReader::from_bytes(&updated).unwrap();
        assert_eq!(header.n_frames, 2);
        assert_eq!(frames.len(), 2);
        assert!(approx_eq(frames[1].x[0], 2.0));
    }

    // ── XTC precision control ─────────────────────────────────────────────

    #[test]
    fn apply_xtc_precision_basic() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[1.23456, 7.89012, -3.15625]]);
        apply_xtc_precision(&mut writer, 2);
        let p = &writer.frames[0].positions[0];
        assert!(approx_eq(p[0], 1.23) || (p[0] - 1.23).abs() < 0.005);
    }

    #[test]
    fn xtc_recompress_roundtrip() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[1.234, 5.678, 9.012]]);
        writer.add_frame(1, 1.0, vec![[0.111, 0.222, 0.333]]);
        let original_bytes = writer.to_bytes();
        let recompressed = xtc_recompress(&original_bytes, 3).unwrap();
        let frames = SimpleXtcReader::from_bytes(&recompressed).unwrap();
        assert_eq!(frames.len(), 2);
        // Values should be close but quantized
        assert!((frames[0].positions[0][0] - 1.234).abs() < 0.001);
    }

    #[test]
    fn compute_rmsd_identical() {
        let positions = vec![[1.0f32, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let rmsd = compute_rmsd(&positions, &positions);
        assert!(
            rmsd.abs() < 1e-6,
            "RMSD of identical sets should be 0, got {rmsd}"
        );
    }

    #[test]
    fn compute_rmsd_known() {
        let a = vec![[0.0f32, 0.0, 0.0]];
        let b = vec![[3.0f32, 4.0, 0.0]]; // distance = 5
        let rmsd = compute_rmsd(&a, &b);
        assert!(approx_eq(rmsd, 5.0), "RMSD should be 5.0, got {rmsd}");
    }

    #[test]
    fn compute_rmsd_different_lengths() {
        let a = vec![[0.0f32, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let b = vec![[0.0f32, 0.0, 0.0]];
        let rmsd = compute_rmsd(&a, &b);
        assert!(approx_eq(rmsd, 0.0)); // returns 0 on mismatch
    }

    // ── Format conversion (XTC ↔ DCD) ────────────────────────────────────

    #[test]
    fn xtc_to_dcd_basic() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[1.0, 2.0, 3.0]]);
        writer.add_frame(1, 1.0, vec![[4.0, 5.0, 6.0]]);
        let xtc_bytes = writer.to_bytes();

        let dcd_bytes = xtc_to_dcd(&xtc_bytes).unwrap();
        let (header, frames) = SimpleDcdReader::from_bytes(&dcd_bytes).unwrap();
        assert_eq!(header.n_atoms, 1);
        assert_eq!(header.n_frames, 2);
        assert_eq!(frames.len(), 2);
        // XTC positions are in nm; DCD stores ångström, i.e. nm * 10.
        assert!(approx_eq(frames[0].x[0], 10.0), "x={}", frames[0].x[0]);
        assert!(approx_eq(frames[0].y[0], 20.0), "y={}", frames[0].y[0]);
        assert!(approx_eq(frames[0].z[0], 30.0), "z={}", frames[0].z[0]);
        assert!(approx_eq(frames[1].x[0], 40.0), "x={}", frames[1].x[0]);
        assert!(approx_eq(frames[1].y[0], 50.0), "y={}", frames[1].y[0]);
        assert!(approx_eq(frames[1].z[0], 60.0), "z={}", frames[1].z[0]);
        // The 1 ps frame interval is stored as the DCD timestep in femtoseconds.
        assert!(
            approx_eq(header.timestep, 1000.0),
            "timestep={}",
            header.timestep
        );
    }

    #[test]
    fn xtc_to_dcd_timestep_is_in_femtoseconds() {
        let mut writer = SimpleXtcWriter::new();
        for i in 0..4 {
            writer.add_frame(i, 0.002 * i as f32, vec![[0.0; 3]]);
        }
        let dcd_bytes = xtc_to_dcd(&writer.to_bytes()).unwrap();
        let (header, _) = SimpleDcdReader::from_bytes(&dcd_bytes).unwrap();
        // 2 fs steps: the header duration is the trajectory length in fs.
        assert!(
            (header.timestep - 2.0).abs() < 1e-4,
            "dt={}",
            header.timestep
        );
        assert!((header.duration() - 8.0).abs() < 1e-3);
    }

    #[test]
    fn xtc_to_dcd_rejects_frames_with_different_atom_counts() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[0.0; 3]; 3]);
        writer.add_frame(1, 1.0, vec![[0.0; 3]; 2]);
        let err = xtc_to_dcd(&writer.to_bytes()).unwrap_err();
        assert!(err.contains("frame 1 has 2 atoms"), "{err}");
        assert!(err.contains("frame 0 has 3"), "{err}");
    }

    #[test]
    fn xtc_to_dcd_empty() {
        let writer = SimpleXtcWriter::new();
        let xtc_bytes = writer.to_bytes();
        let dcd_bytes = xtc_to_dcd(&xtc_bytes).unwrap();
        let (header, frames) = SimpleDcdReader::from_bytes(&dcd_bytes).unwrap();
        assert_eq!(header.n_frames, 0);
        assert!(frames.is_empty());
    }

    #[test]
    fn dcd_to_xtc_basic() {
        let mut writer = SimpleDcdWriter::new(2, 1.0);
        writer.add_frame(&[[10.0, 20.0, 30.0], [40.0, 50.0, 60.0]]);
        writer.add_frame(&[[11.0, 21.0, 31.0], [41.0, 51.0, 61.0]]);
        let dcd_bytes = writer.to_bytes();

        let xtc_bytes = dcd_to_xtc(&dcd_bytes).unwrap();
        let frames = SimpleXtcReader::from_bytes(&xtc_bytes).unwrap();
        assert_eq!(frames.len(), 2);
        // DCD is Å, XTC should be nm = /10
        assert!(approx_eq(frames[0].positions[0][0], 1.0));
        assert!(approx_eq(frames[0].positions[1][2], 6.0));
        assert!(approx_eq(frames[1].positions[0][0], 1.1));
    }

    #[test]
    fn dcd_to_xtc_timestep_is_converted_from_femtoseconds() {
        // 2 fs DCD timestep = 0.002 ps between XTC frames.
        let mut writer = SimpleDcdWriter::new(1, 2.0);
        for _ in 0..3 {
            writer.add_frame(&[[0.0; 3]]);
        }
        let frames = SimpleXtcReader::from_bytes(&dcd_to_xtc(&writer.to_bytes()).unwrap()).unwrap();
        assert_eq!(frames.len(), 3);
        for (idx, frame) in frames.iter().enumerate() {
            assert_eq!(frame.step, idx as i32);
            let expected = 0.002 * idx as f32;
            assert!(
                (frame.time - expected).abs() < 1e-6,
                "frame {idx}: time={} expected={expected}",
                frame.time
            );
        }
    }

    #[test]
    fn xtc_dcd_xtc_roundtrip() {
        // XTC (nm, ps) -> DCD (Å, fs) -> XTC must reproduce positions and times:
        // x10 then /10 on the positions, x1000 then /1000 on the timestep.
        let frame0 = vec![[10.0f32, 20.0, 30.0], [40.0, 50.0, 60.0]];
        let frame1 = vec![[10.5f32, 20.25, 29.0], [41.0, 50.125, 61.5]];
        let frame2 = vec![[11.0f32, 20.5, 28.0], [42.0, 50.25, 63.0]];
        let mut xtc_writer = SimpleXtcWriter::new();
        xtc_writer.add_frame(0, 0.0, frame0.clone());
        xtc_writer.add_frame(1, 0.5, frame1.clone());
        xtc_writer.add_frame(2, 1.0, frame2.clone());
        let xtc_bytes = xtc_writer.to_bytes();

        let dcd_bytes = xtc_to_dcd(&xtc_bytes).unwrap();
        let xtc2_bytes = dcd_to_xtc(&dcd_bytes).unwrap();
        let frames = SimpleXtcReader::from_bytes(&xtc2_bytes).unwrap();

        assert_eq!(frames.len(), 3);
        for (idx, original) in [frame0, frame1, frame2].iter().enumerate() {
            assert_eq!(frames[idx].step, idx as i32);
            assert!(
                (frames[idx].time - 0.5 * idx as f32).abs() < 1e-5,
                "frame {idx}: time {}",
                frames[idx].time
            );
            assert_eq!(frames[idx].positions.len(), original.len());
            for (got, exp) in frames[idx].positions.iter().zip(original.iter()) {
                for k in 0..3 {
                    assert!(
                        (got[k] - exp[k]).abs() < 1e-5,
                        "frame {idx}: got {} expected {}",
                        got[k],
                        exp[k]
                    );
                }
            }
        }
    }

    // ── Frame time extraction ─────────────────────────────────────────────

    #[test]
    fn xtc_extract_frame_times_basic() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[0.0; 3]]);
        writer.add_frame(1, 2.5, vec![[0.0; 3]]);
        writer.add_frame(2, 5.0, vec![[0.0; 3]]);
        let bytes = writer.to_bytes();
        let times = xtc_extract_frame_times(&bytes).unwrap();
        assert_eq!(times.len(), 3);
        assert!(approx_eq(times[0], 0.0));
        assert!(approx_eq(times[1], 2.5));
        assert!(approx_eq(times[2], 5.0));
    }

    #[test]
    fn xtc_extract_frame_steps_basic() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(10, 0.0, vec![[0.0; 3]]);
        writer.add_frame(20, 1.0, vec![[0.0; 3]]);
        let bytes = writer.to_bytes();
        let steps = xtc_extract_frame_steps(&bytes).unwrap();
        assert_eq!(steps, vec![10, 20]);
    }

    #[test]
    fn xtc_frame_time_deltas_basic() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[0.0; 3]]);
        writer.add_frame(1, 1.0, vec![[0.0; 3]]);
        writer.add_frame(2, 3.0, vec![[0.0; 3]]);
        let bytes = writer.to_bytes();
        let deltas = xtc_frame_time_deltas(&bytes).unwrap();
        assert_eq!(deltas.len(), 2);
        assert!(approx_eq(deltas[0], 1.0));
        assert!(approx_eq(deltas[1], 2.0));
    }

    #[test]
    fn xtc_frame_time_deltas_single_frame() {
        let mut writer = SimpleXtcWriter::new();
        writer.add_frame(0, 0.0, vec![[0.0; 3]]);
        let bytes = writer.to_bytes();
        let deltas = xtc_frame_time_deltas(&bytes).unwrap();
        assert!(deltas.is_empty());
    }

    // ── Velocity data support ─────────────────────────────────────────────

    #[test]
    fn xtc_frame_with_velocity_kinetic_energy() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[1.0; 3]; 3],
            positions: vec![[0.0; 3]],
        };
        let velocities = vec![[3.0f32, 4.0, 0.0]]; // speed = 5, KE = 0.5 * 25 = 12.5
        let wv = XtcFrameWithVelocity::new(frame, velocities);
        assert!((wv.kinetic_energy_proxy() - 12.5).abs() < 1e-5);
    }

    #[test]
    fn xtc_frame_with_velocity_rms_speed() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[1.0; 3]; 3],
            positions: vec![[0.0; 3]; 4],
        };
        // 4 atoms each moving at speed 2 in x
        let velocities = vec![[2.0f32, 0.0, 0.0]; 4];
        let wv = XtcFrameWithVelocity::new(frame, velocities);
        assert!(approx_eq(wv.rms_speed(), 2.0));
    }

    #[test]
    fn xtc_frame_with_velocity_max_speed() {
        let frame = XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[1.0; 3]; 3],
            positions: vec![[0.0; 3]; 3],
        };
        let velocities = vec![[1.0f32, 0.0, 0.0], [3.0, 4.0, 0.0], [0.1, 0.0, 0.0]];
        let wv = XtcFrameWithVelocity::new(frame, velocities);
        assert!(approx_eq(wv.max_speed(), 5.0)); // sqrt(9+16) = 5
    }

    #[test]
    fn xtc_velocity_writer_basic() {
        let mut writer = XtcVelocityWriter::new();
        let pos = vec![[1.0f32, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let vel = vec![[0.1f32, 0.0, 0.0], [0.0, 0.2, 0.0]];
        writer.add_frame(0, 0.0, pos, vel);
        assert_eq!(writer.frame_count(), 1);
        let bytes = writer.to_bytes();
        assert!(!bytes.is_empty());
        // Bytes should contain the velocity magic
        let magic_pos = bytes.windows(4).position(|w| w == b"OXVL");
        assert!(magic_pos.is_some(), "should contain OXVL magic");
    }

    #[test]
    fn xtc_velocity_writer_kinetic_energy() {
        let mut writer = XtcVelocityWriter::new();
        let pos = vec![[0.0f32; 3]];
        let vel = vec![[3.0f32, 4.0, 0.0]]; // KE = 0.5 * 25 = 12.5
        writer.add_frame(0, 0.0, pos, vel);
        let ke = writer.frame_kinetic_energy(0);
        assert!((ke - 12.5).abs() < 1e-5);
        // Out of range
        assert!(approx_eq(writer.frame_kinetic_energy(99), 0.0));
    }

    #[test]
    fn compute_finite_difference_velocities_basic() {
        let frames = vec![
            XtcFrame {
                step: 0,
                time: 0.0,
                box_matrix: [[1.0; 3]; 3],
                positions: vec![[0.0, 0.0, 0.0]],
            },
            XtcFrame {
                step: 1,
                time: 1.0,
                box_matrix: [[1.0; 3]; 3],
                positions: vec![[2.0, 0.0, 0.0]],
            },
        ];
        let vels = compute_finite_difference_velocities(&frames, 1.0);
        assert_eq!(vels.len(), 1);
        assert!(approx_eq(vels[0][0][0], 2.0)); // (2-0)/1 = 2 nm/ps
    }

    #[test]
    fn compute_finite_difference_velocities_three_frames() {
        let frames: Vec<XtcFrame> = (0..3)
            .map(|i| XtcFrame {
                step: i,
                time: i as f32,
                box_matrix: [[1.0; 3]; 3],
                positions: vec![[i as f32 * 3.0, 0.0, 0.0]],
            })
            .collect();
        let vels = compute_finite_difference_velocities(&frames, 1.0);
        // 3 frames → 2 velocity frames
        assert_eq!(vels.len(), 2);
        // Each frame: dx = 3 over dt = 1
        assert!(approx_eq(vels[0][0][0], 3.0));
        assert!(approx_eq(vels[1][0][0], 3.0));
    }

    #[test]
    fn compute_finite_difference_velocities_single_frame() {
        let frames = vec![XtcFrame {
            step: 0,
            time: 0.0,
            box_matrix: [[1.0; 3]; 3],
            positions: vec![[0.0; 3]],
        }];
        let vels = compute_finite_difference_velocities(&frames, 1.0);
        assert!(vels.is_empty());
    }
}
