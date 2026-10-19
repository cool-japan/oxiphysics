// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! CHARMM / NAMD / X-PLOR binary DCD trajectories.
//!
//! A DCD file is a sequence of Fortran unformatted records.  Every record is
//! `length | payload | length` with 32-bit lengths in bytes:
//!
//! | Record     | Payload                                                         |
//! |------------|-----------------------------------------------------------------|
//! | header     | 84 bytes: the magic `CORD` and the 20 32-bit control words `ICNTRL` |
//! | titles     | `4 + 80 * ntitle` bytes: `ntitle` and the 80-byte title lines   |
//! | atom count | 4 bytes: the number of atoms                                    |
//! | per frame  | an optional unit-cell record, then the `x`, `y` and `z` coordinates (and `w` in 4D files) |
//!
//! The control words (0-based, as in the VMD and NAMD sources) used here:
//!
//! | Word | Meaning                                                              |
//! |------|----------------------------------------------------------------------|
//! | 0    | `NSET`, number of frames                                             |
//! | 1    | `ISTART`, first step                                                 |
//! | 2    | `NSAVC`, steps between two frames                                    |
//! | 3    | `NSTEP`, last step                                                   |
//! | 8    | `NAMNF`, number of fixed atoms                                       |
//! | 9    | `DELTA`, the time step: an `f32` in CHARMM files, an `f64` spanning words 9 and 10 in X-PLOR files |
//! | 10   | unit-cell flag (CHARMM files)                                        |
//! | 11   | 4D flag (CHARMM files)                                               |
//! | 19   | CHARMM version, zero for X-PLOR files                                |
//!
//! The unit-cell record holds six `f64` in the order `A, cos(gamma), B,
//! cos(beta), cos(alpha), C` (older NAMD versions store the angles in degrees
//! instead of their cosines).  Each coordinate record holds `natom` `f32`.
//!
//! Both byte orders are read: a file whose first record marker reads `84`
//! only after a byte swap is big-endian.  Files that hold fixed atoms
//! (`NAMNF > 0`, the free-atom index table and reduced frames) and velocity
//! trajectories (magic `VELO`) are rejected with
//! [`io::ErrorKind::Unsupported`](std::io::ErrorKind::Unsupported).
//!
//! The reader derives the number of frames from the size of the file, since
//! every frame has the same size.  The `NSET` word of the header must not
//! exceed it, and the file must not end inside a frame; either is reported as
//! an error instead of silently returning fewer frames.  A header that was
//! never finalised (`NSET = 0`) or is out of date does not hide the frames
//! that follow.

use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use super::ParticleFrame;

/// Length in bytes of the first record: `CORD` and the 20 control words.
const HEADER_RECORD_BYTES: u32 = 84;
/// Number of 32-bit control words in the header record.
const ICNTRL_WORDS: usize = 20;
/// Length in bytes of one title line.
const TITLE_BYTES: u32 = 80;
/// Length in bytes of the unit-cell record: six `f64`.
const UNIT_CELL_RECORD_BYTES: u32 = 48;
/// Bytes of record markers around one record.
const MARKER_BYTES: u64 = 8;

// Positions of the control words in `ICNTRL`.
const IC_NSET: usize = 0;
const IC_ISTART: usize = 1;
const IC_NSAVC: usize = 2;
const IC_NSTEP: usize = 3;
const IC_NAMNF: usize = 8;
const IC_DELTA: usize = 9;
const IC_UNIT_CELL: usize = 10;
const IC_FOUR_D: usize = 11;
const IC_VERSION: usize = 19;

/// Version word written by [`DcdWriter`]; any non-zero value marks the CHARMM
/// layout (single-precision `DELTA`, unit-cell and 4D flags).
const CHARMM_VERSION: u32 = 24;
/// First step and step interval recorded by [`DcdWriter`].
const WRITER_FIRST_STEP: u32 = 0;
const WRITER_STEP_INTERVAL: u32 = 1;
/// File offsets of the words [`DcdWriter::finalize`] patches: the record
/// marker and `CORD` take eight bytes, then `ICNTRL[k]` lives at `8 + 4 k`.
const NSET_OFFSET: u64 = 8;
const NSTEP_OFFSET: u64 = 20;
/// Upper bound for the up-front reservation of the frame list.
const MAX_PREALLOCATED_FRAMES: usize = 4096;

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Header information parsed from a DCD file.
#[derive(Debug, Clone, Default)]
pub struct DcdHeader {
    /// Number of complete frames stored in the file.  It is derived from the
    /// size of the file; the `NSET` word of the header never exceeds it.
    pub n_frames: u32,
    /// First timestep number (`ISTART`).
    pub first_step: u32,
    /// Number of steps between saved frames (`NSAVC`).
    pub step_interval: u32,
    /// Number of atoms.
    pub n_atoms: u32,
    /// Title strings from the header, without their padding.
    pub titles: Vec<String>,
    /// Last timestep number (`NSTEP`).
    pub n_steps: u32,
    /// Integration time step as stored in the file (`DELTA`).  CHARMM and NAMD
    /// store AKMA time units; the value is not converted.
    pub delta: f64,
    /// `true` if every frame starts with a unit-cell record (CHARMM layout
    /// with the unit-cell flag set).
    pub has_unit_cell: bool,
    /// CHARMM version word; zero for X-PLOR files, which keep `DELTA` as a
    /// double and have no unit-cell flag.
    pub charmm_version: u32,
    /// `true` if the file was written in big-endian byte order.
    pub big_endian: bool,
}

/// The periodic cell stored in front of the coordinates of a frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DcdUnitCell {
    /// Edge lengths `[a, b, c]` (ångström).
    pub lengths: [f64; 3],
    /// Angles `[alpha, beta, gamma]` in degrees.
    pub angles: [f64; 3],
}

impl DcdUnitCell {
    /// Orthorhombic cell: all three angles are 90 degrees.
    pub fn orthorhombic(lengths: [f64; 3]) -> Self {
        Self {
            lengths,
            angles: [90.0; 3],
        }
    }

    /// Decode the six values of a unit-cell record, in file order
    /// `[A, gamma, B, beta, alpha, C]`.
    ///
    /// CHARMM and NAMD 2.5 or later store the cosines of the angles; older
    /// NAMD versions store the angles in degrees.  The record is read as
    /// cosines when all three angle entries lie in `[-1, 1]` (the rule VMD
    /// applies) and as degrees otherwise.
    pub fn from_record(record: [f64; 6]) -> Self {
        let [a, gamma, b, beta, alpha, c] = record;
        let cosines = [gamma, beta, alpha]
            .iter()
            .all(|value| (-1.0..=1.0).contains(value));
        let to_degrees = |value: f64| {
            if cosines {
                90.0 - value.asin().to_degrees()
            } else {
                value
            }
        };
        Self {
            lengths: [a, b, c],
            angles: [to_degrees(alpha), to_degrees(beta), to_degrees(gamma)],
        }
    }

    /// Encode the cell as a unit-cell record `[A, cos(gamma), B, cos(beta),
    /// cos(alpha), C]`, the inverse of [`DcdUnitCell::from_record`].
    pub fn to_record(&self) -> [f64; 6] {
        let cosine = |degrees: f64| {
            if degrees == 90.0 {
                0.0
            } else {
                degrees.to_radians().cos()
            }
        };
        let [a, b, c] = self.lengths;
        let [alpha, beta, gamma] = self.angles;
        [a, cosine(gamma), b, cosine(beta), cosine(alpha), c]
    }
}

/// A complete DCD file: header, frames and, for files with a unit-cell flag,
/// one cell per frame.
#[derive(Debug, Clone)]
pub struct DcdTrajectory {
    /// Parsed header.
    pub header: DcdHeader,
    /// Coordinate frames (`timestep` is the frame index).
    pub frames: Vec<ParticleFrame>,
    /// Periodic cell of each frame; empty unless `header.has_unit_cell`.
    pub unit_cells: Vec<DcdUnitCell>,
}

// ─────────────────────────────────────────────────────────────────────────────
// Reader
// ─────────────────────────────────────────────────────────────────────────────

/// DCD binary reader: parses the header and reads the coordinate frames.
///
/// Reads CHARMM, NAMD and X-PLOR trajectories in either byte order, including
/// the unit-cell record of periodic simulations.  See the
/// [module documentation](self) for the layout and the unsupported variants.
#[derive(Debug)]
pub struct DcdReader {
    /// Input path.
    pub path: std::path::PathBuf,
}

impl DcdReader {
    /// Create a new `DcdReader` for `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Parse only the DCD file header.
    ///
    /// # Errors
    ///
    /// Same as [`DcdReader::read_trajectory`]; the frame data is not decoded
    /// but its size is checked against the header.
    pub fn read_header(&self) -> io::Result<DcdHeader> {
        let file = File::open(&self.path)?;
        let len = file.metadata()?.len();
        let mut stream = DcdStream::new(BufReader::new(file), len);
        let layout = parse_dcd_header(&mut stream)?;
        let n_frames = layout.count_frames(len)?;
        layout.into_header(n_frames)
    }

    /// Read all coordinate frames.  Returns `(header, frames)`.
    ///
    /// Use [`DcdReader::read_trajectory`] to also get the unit cells.
    ///
    /// # Errors
    ///
    /// Same as [`DcdReader::read_trajectory`].
    pub fn read_all(&self) -> io::Result<(DcdHeader, Vec<ParticleFrame>)> {
        let trajectory = self.read_trajectory()?;
        Ok((trajectory.header, trajectory.frames))
    }

    /// Read the header, all coordinate frames and the unit cells.
    ///
    /// # Errors
    ///
    /// * `io::ErrorKind::InvalidData` if the first record is not the 84-byte
    ///   header record, the magic is not `CORD`, a record marker disagrees
    ///   with the size the layout requires (title block, atom count, unit
    ///   cell, or a coordinate record that is not `4 * natom` bytes long), or
    ///   the frame data does not fit the 32-bit frame count.
    /// * `io::ErrorKind::UnexpectedEof` if the file ends inside a record or a
    ///   frame, or holds fewer frames than the `NSET` word of the header
    ///   announces.
    /// * `io::ErrorKind::Unsupported` for velocity trajectories (`VELO`) and
    ///   for files with fixed atoms (`NAMNF > 0`).
    pub fn read_trajectory(&self) -> io::Result<DcdTrajectory> {
        let file = File::open(&self.path)?;
        let len = file.metadata()?.len();
        read_dcd_trajectory(BufReader::new(file), len)
    }
}

/// Byte order of a DCD file.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ByteOrder {
    Little,
    Big,
}

impl ByteOrder {
    fn u32(self, bytes: [u8; 4]) -> u32 {
        match self {
            Self::Little => u32::from_le_bytes(bytes),
            Self::Big => u32::from_be_bytes(bytes),
        }
    }

    fn f32(self, bytes: [u8; 4]) -> f32 {
        match self {
            Self::Little => f32::from_le_bytes(bytes),
            Self::Big => f32::from_be_bytes(bytes),
        }
    }

    fn f64(self, bytes: [u8; 8]) -> f64 {
        match self {
            Self::Little => f64::from_le_bytes(bytes),
            Self::Big => f64::from_be_bytes(bytes),
        }
    }
}

fn invalid_data(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

fn unsupported(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::Unsupported, message)
}

/// Append `context` to the message of `error`, keeping its kind.
fn with_context(error: io::Error, context: &str) -> io::Error {
    io::Error::new(error.kind(), format!("{error} ({context})"))
}

/// A byte stream that tracks how much of a file of known length is left, and
/// decodes the 32-bit record markers of a DCD file.
struct DcdStream<R: Read> {
    inner: R,
    order: ByteOrder,
    pos: u64,
    len: u64,
}

impl<R: Read> DcdStream<R> {
    fn new(inner: R, len: u64) -> Self {
        Self {
            inner,
            order: ByteOrder::Little,
            pos: 0,
            len,
        }
    }

    /// Bytes between the current position and the end of the file.
    fn remaining(&self) -> u64 {
        self.len.saturating_sub(self.pos)
    }

    /// Fill `buf`; a short file is reported as `UnexpectedEof` naming `what`.
    fn fill(&mut self, buf: &mut [u8], what: &str) -> io::Result<()> {
        match self.inner.read_exact(buf) {
            Ok(()) => {
                self.pos += buf.len() as u64;
                Ok(())
            }
            Err(error) if error.kind() == io::ErrorKind::UnexpectedEof => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "DCD: unexpected end of file while reading {what} at byte {}",
                    self.pos
                ),
            )),
            Err(error) => Err(error),
        }
    }

    fn read_u32(&mut self, what: &str) -> io::Result<u32> {
        let mut bytes = [0u8; 4];
        self.fill(&mut bytes, what)?;
        Ok(self.order.u32(bytes))
    }

    /// Read a record marker and require it to equal `expected`.
    fn expect_marker(&mut self, expected: u32, what: &str) -> io::Result<()> {
        let found = self.read_u32(what)?;
        if found == expected {
            Ok(())
        } else {
            Err(invalid_data(format!(
                "DCD: {what} is {found}, expected {expected}"
            )))
        }
    }
}

/// What the header says about the data that follows it.
struct DcdLayout {
    header: DcdHeader,
    /// `NSET` word of the header.
    declared_frames: u32,
    has_unit_cell: bool,
    four_d: bool,
    order: ByteOrder,
    /// Offset of the first frame.
    data_start: u64,
}

impl DcdLayout {
    /// Size in bytes of one coordinate record, markers included.
    fn coordinate_record_bytes(&self) -> u64 {
        MARKER_BYTES + 4 * u64::from(self.header.n_atoms)
    }

    /// Size in bytes of one frame.
    fn frame_bytes(&self) -> io::Result<u64> {
        let axes: u64 = if self.four_d { 4 } else { 3 };
        let cell = if self.has_unit_cell {
            MARKER_BYTES + u64::from(UNIT_CELL_RECORD_BYTES)
        } else {
            0
        };
        axes.checked_mul(self.coordinate_record_bytes())
            .and_then(|coordinates| coordinates.checked_add(cell))
            .ok_or_else(|| {
                invalid_data(format!(
                    "DCD: a frame of {} atoms does not fit a 64-bit size",
                    self.header.n_atoms
                ))
            })
    }

    /// Number of frames in a file of `file_len` bytes.
    ///
    /// The frames all have the same size, so the count follows from the file
    /// length.  A file that ends inside a frame, or that holds fewer frames
    /// than the header announces, is an error.
    fn count_frames(&self, file_len: u64) -> io::Result<u64> {
        let frame = self.frame_bytes()?;
        let data = file_len.saturating_sub(self.data_start);
        let complete = data / frame;
        let partial = data % frame;
        if partial != 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "DCD: the file ends inside frame {complete}: {partial} of the {frame} bytes \
                     of the frame are present"
                ),
            ));
        }
        if u64::from(self.declared_frames) > complete {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!(
                    "DCD: the header announces {} frames but the file holds only {complete}",
                    self.declared_frames
                ),
            ));
        }
        Ok(complete)
    }

    /// The public header, with the frame count `n_frames` found in the file.
    fn into_header(self, n_frames: u64) -> io::Result<DcdHeader> {
        let n_frames = u32::try_from(n_frames)
            .map_err(|_| invalid_data(format!("DCD: {n_frames} frames exceed the 32-bit count")))?;
        Ok(DcdHeader {
            n_frames,
            ..self.header
        })
    }
}

/// Read the header, title and atom-count records.
fn parse_dcd_header<R: Read>(stream: &mut DcdStream<R>) -> io::Result<DcdLayout> {
    // ── Record 1: 84 bytes, "CORD" and the 20 control words ──────────────────
    let mut marker = [0u8; 4];
    stream.fill(&mut marker, "the header record marker")?;
    stream.order = if u32::from_le_bytes(marker) == HEADER_RECORD_BYTES {
        ByteOrder::Little
    } else if u32::from_be_bytes(marker) == HEADER_RECORD_BYTES {
        ByteOrder::Big
    } else {
        return Err(invalid_data(format!(
            "DCD: the first record is {} bytes long, a DCD header record has \
             {HEADER_RECORD_BYTES}",
            u32::from_le_bytes(marker)
        )));
    };
    let order = stream.order;

    let mut magic = [0u8; 4];
    stream.fill(&mut magic, "the magic")?;
    match &magic {
        b"CORD" => {}
        b"VELO" => {
            return Err(unsupported(
                "DCD: velocity trajectories (magic VELO) are not supported".to_string(),
            ));
        }
        other => {
            return Err(invalid_data(format!(
                "DCD: bad magic {other:?}, expected b\"CORD\""
            )));
        }
    }
    let mut raw = [0u8; 4 * ICNTRL_WORDS];
    stream.fill(&mut raw, "the control words")?;
    stream.expect_marker(HEADER_RECORD_BYTES, "the header closing marker")?;

    let mut icntrl = [0u32; ICNTRL_WORDS];
    let (chunks, _) = raw.as_chunks::<4>();
    for (word, chunk) in icntrl.iter_mut().zip(chunks) {
        *word = order.u32(*chunk);
    }
    let charmm_version = icntrl[IC_VERSION];
    let charmm = charmm_version != 0;
    if icntrl[IC_NAMNF] != 0 {
        return Err(unsupported(format!(
            "DCD: files with fixed atoms (NAMNF = {}) are not supported",
            icntrl[IC_NAMNF]
        )));
    }
    // CHARMM files keep DELTA in one single-precision word and use the next
    // words as unit-cell and 4D flags.  X-PLOR files store DELTA as a double
    // that spans words 9 and 10, so word 10 is data, not a flag.
    let delta = if charmm {
        f64::from(f32::from_bits(icntrl[IC_DELTA]))
    } else {
        order.f64(std::array::from_fn(|i| raw[4 * IC_DELTA + i]))
    };
    let has_unit_cell = charmm && icntrl[IC_UNIT_CELL] != 0;
    let four_d = charmm && icntrl[IC_FOUR_D] == 1;

    // ── Record 2: titles ─────────────────────────────────────────────────────
    let title_bytes = stream.read_u32("the title record marker")?;
    let payload = title_bytes
        .checked_sub(4)
        .filter(|p| p.is_multiple_of(TITLE_BYTES));
    let Some(payload) = payload else {
        return Err(invalid_data(format!(
            "DCD: the title record is {title_bytes} bytes long, expected 4 + {TITLE_BYTES} per title"
        )));
    };
    let n_titles = payload / TITLE_BYTES;
    if u64::from(title_bytes) + 4 > stream.remaining() {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            format!(
                "DCD: the title record announces {title_bytes} bytes but only {} remain",
                stream.remaining()
            ),
        ));
    }
    let declared_titles = stream.read_u32("the number of titles")?;
    if declared_titles != n_titles {
        return Err(invalid_data(format!(
            "DCD: the title record holds {n_titles} lines but announces {declared_titles}"
        )));
    }
    let mut titles = Vec::new();
    for _ in 0..n_titles {
        let mut line = [0u8; TITLE_BYTES as usize];
        stream.fill(&mut line, "a title line")?;
        titles.push(
            String::from_utf8_lossy(&line)
                .trim_end_matches([' ', '\0'])
                .to_string(),
        );
    }
    stream.expect_marker(title_bytes, "the title closing marker")?;

    // ── Record 3: number of atoms ────────────────────────────────────────────
    stream.expect_marker(4, "the atom count record marker")?;
    let n_atoms = stream.read_u32("the number of atoms")?;
    stream.expect_marker(4, "the atom count closing marker")?;
    if n_atoms > u32::MAX / 4 {
        return Err(invalid_data(format!(
            "DCD: {n_atoms} atoms cannot be stored: a coordinate record would exceed the 32-bit \
             record length"
        )));
    }

    Ok(DcdLayout {
        header: DcdHeader {
            n_frames: 0,
            first_step: icntrl[IC_ISTART],
            step_interval: icntrl[IC_NSAVC],
            n_atoms,
            titles,
            n_steps: icntrl[IC_NSTEP],
            delta,
            has_unit_cell,
            charmm_version,
            big_endian: order == ByteOrder::Big,
        },
        declared_frames: icntrl[IC_NSET],
        has_unit_cell,
        four_d,
        order,
        data_start: stream.pos,
    })
}

/// Read one coordinate record of `n_atoms` values into `scratch`.
fn read_coordinate_record<R: Read>(
    stream: &mut DcdStream<R>,
    n_atoms: u32,
    scratch: &mut [u8],
    what: &str,
) -> io::Result<()> {
    // `4 * n_atoms` cannot overflow: the header parser rejects larger counts.
    let payload = 4 * n_atoms;
    stream.expect_marker(payload, what)?;
    stream.fill(scratch, what)?;
    stream.expect_marker(payload, what)
}

/// Parse one frame: the optional unit cell and the coordinate records.
fn parse_dcd_frame<R: Read>(
    stream: &mut DcdStream<R>,
    layout: &DcdLayout,
    index: u64,
    scratch: &mut Vec<u8>,
) -> io::Result<(ParticleFrame, Option<DcdUnitCell>)> {
    let order = layout.order;
    let n_atoms = layout.header.n_atoms;

    let unit_cell = if layout.has_unit_cell {
        stream.expect_marker(UNIT_CELL_RECORD_BYTES, "the unit cell record marker")?;
        let mut raw = [0u8; UNIT_CELL_RECORD_BYTES as usize];
        stream.fill(&mut raw, "the unit cell")?;
        stream.expect_marker(UNIT_CELL_RECORD_BYTES, "the unit cell closing marker")?;
        let mut values = [0.0f64; 6];
        let (chunks, _) = raw.as_chunks::<8>();
        for (value, chunk) in values.iter_mut().zip(chunks) {
            *value = order.f64(*chunk);
        }
        Some(DcdUnitCell::from_record(values))
    } else {
        None
    };

    // The frame list was sized from the file length, so `n_atoms` records of
    // this size are known to be present.
    scratch.resize((4 * n_atoms) as usize, 0);
    let mut positions = vec![[0.0f32; 3]; n_atoms as usize];
    for (axis, what) in ["x coordinates", "y coordinates", "z coordinates"]
        .into_iter()
        .enumerate()
    {
        read_coordinate_record(stream, n_atoms, scratch, what)?;
        let (chunks, _) = scratch.as_chunks::<4>();
        for (position, chunk) in positions.iter_mut().zip(chunks) {
            position[axis] = order.f32(*chunk);
        }
    }
    if layout.four_d {
        read_coordinate_record(stream, n_atoms, scratch, "fourth dimension coordinates")?;
    }

    let frame = ParticleFrame {
        timestep: index,
        n_particles: positions.len(),
        positions,
        velocities: None,
        types: vec![0u8; n_atoms as usize],
    };
    Ok((frame, unit_cell))
}

/// Read a whole DCD file of `len` bytes from `reader`.
pub(super) fn read_dcd_trajectory<R: Read>(reader: R, len: u64) -> io::Result<DcdTrajectory> {
    let mut stream = DcdStream::new(reader, len);
    let layout = parse_dcd_header(&mut stream)?;
    let n_frames = layout.count_frames(len)?;
    let capacity = usize::try_from(n_frames)
        .map_err(|_| invalid_data(format!("DCD: {n_frames} frames exceed the address space")))?
        .min(MAX_PREALLOCATED_FRAMES);

    let mut frames = Vec::with_capacity(capacity);
    let mut unit_cells = Vec::new();
    let mut scratch = Vec::new();
    for index in 0..n_frames {
        let (frame, cell) = parse_dcd_frame(&mut stream, &layout, index, &mut scratch)
            .map_err(|error| with_context(error, &format!("frame {index}")))?;
        frames.push(frame);
        unit_cells.extend(cell);
    }
    Ok(DcdTrajectory {
        header: layout.into_header(n_frames)?,
        frames,
        unit_cells,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Writer
// ─────────────────────────────────────────────────────────────────────────────

/// Write a slice of `f32` values as a Fortran-style binary record.
///
/// The record is bracketed by a leading and trailing 4-byte little-endian
/// integer whose value equals `data.len() * 4` (the byte count of the payload).
fn write_f32_fortran_block<W: Write>(w: &mut W, data: &[f32]) -> io::Result<()> {
    let byte_len = (data.len() * 4) as u32;
    w.write_all(&byte_len.to_le_bytes())?;
    for &v in data {
        w.write_all(&v.to_le_bytes())?;
    }
    w.write_all(&byte_len.to_le_bytes())?;
    Ok(())
}

/// DCD trajectory writer — produces files in the standard CHARMM/NAMD binary
/// DCD format that can be read by [`DcdReader`] and common molecular dynamics
/// visualisation tools (VMD, NAMD, etc.).
///
/// The header is the 84-byte CHARMM record: the control words hold the frame
/// count, first step (0), step interval (1) and last step, no fixed atoms, the
/// time step as a single-precision float, the unit-cell flag and the CHARMM
/// version 24.  Files are little-endian.
///
/// # Usage
///
/// ```rust,no_run
/// use oxiphysics_io::DcdWriter;
///
/// let path = std::env::temp_dir().join("traj.dcd");
/// let mut w = DcdWriter::new(&path, 10, 2.0e-3).unwrap();
/// // write frames …
/// w.finalize().unwrap();  // patches the frame count in the header
/// ```
#[derive(Debug)]
pub struct DcdWriter {
    /// Output path.
    pub path: std::path::PathBuf,
    /// Number of atoms per frame.
    pub n_atoms: u32,
    /// Integration time-step (ps).
    pub delta: f32,
    has_unit_cell: bool,
    frames_written: u32,
}

impl DcdWriter {
    /// Create a new `DcdWriter`, truncating / creating the output file and
    /// writing an initial DCD header.
    ///
    /// # Errors
    ///
    /// `io::ErrorKind::InvalidInput` if a coordinate record of `n_atoms`
    /// values would not fit the 32-bit record length, and any I/O error of
    /// creating the file.
    pub fn new(path: impl Into<std::path::PathBuf>, n_atoms: u32, delta: f32) -> io::Result<Self> {
        Self::create(path.into(), n_atoms, delta, false)
    }

    /// Like [`DcdWriter::new`], but every frame carries the periodic cell of
    /// the simulation box (the unit-cell flag of the header is set).  Frames
    /// are then written with [`DcdWriter::write_frame_with_cell`].
    ///
    /// # Errors
    ///
    /// Same as [`DcdWriter::new`].
    pub fn with_unit_cell(
        path: impl Into<std::path::PathBuf>,
        n_atoms: u32,
        delta: f32,
    ) -> io::Result<Self> {
        Self::create(path.into(), n_atoms, delta, true)
    }

    fn create(
        path: std::path::PathBuf,
        n_atoms: u32,
        delta: f32,
        has_unit_cell: bool,
    ) -> io::Result<Self> {
        if n_atoms.checked_mul(4).is_none() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("DCD: {n_atoms} atoms do not fit a 32-bit record length"),
            ));
        }
        let mut writer = Self {
            path,
            n_atoms,
            delta,
            has_unit_cell,
            frames_written: 0,
        };
        writer.write_header(0)?;
        Ok(writer)
    }

    /// Append one coordinate frame to the DCD file.
    ///
    /// The frame must contain at least `n_atoms` positions; extra positions are
    /// silently ignored, missing ones are padded with `0.0`.
    ///
    /// # Errors
    ///
    /// `io::ErrorKind::InvalidInput` if the writer was created with
    /// [`DcdWriter::with_unit_cell`], whose frames need a cell, and any I/O
    /// error of appending to the file.
    pub fn write_frame(&mut self, frame: &ParticleFrame) -> io::Result<()> {
        if self.has_unit_cell {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "DCD: this writer stores a unit cell per frame, use write_frame_with_cell",
            ));
        }
        self.append_frame(frame, None)
    }

    /// Append one coordinate frame together with its periodic cell.
    ///
    /// Positions are padded or truncated as in [`DcdWriter::write_frame`].
    ///
    /// # Errors
    ///
    /// `io::ErrorKind::InvalidInput` if the writer was not created with
    /// [`DcdWriter::with_unit_cell`], and any I/O error of appending to the
    /// file.
    pub fn write_frame_with_cell(
        &mut self,
        frame: &ParticleFrame,
        cell: &DcdUnitCell,
    ) -> io::Result<()> {
        if !self.has_unit_cell {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "DCD: this writer has no unit cell flag, create it with with_unit_cell",
            ));
        }
        self.append_frame(frame, Some(cell))
    }

    fn append_frame(
        &mut self,
        frame: &ParticleFrame,
        cell: Option<&DcdUnitCell>,
    ) -> io::Result<()> {
        let n = self.n_atoms as usize;

        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);
        let mut zs = Vec::with_capacity(n);
        for i in 0..n {
            if let Some(position) = frame.positions.get(i) {
                xs.push(position[0]);
                ys.push(position[1]);
                zs.push(position[2]);
            } else {
                xs.push(0.0_f32);
                ys.push(0.0_f32);
                zs.push(0.0_f32);
            }
        }

        let file = OpenOptions::new().append(true).open(&self.path)?;
        let mut w = BufWriter::new(file);
        if let Some(cell) = cell {
            w.write_all(&UNIT_CELL_RECORD_BYTES.to_le_bytes())?;
            for value in cell.to_record() {
                w.write_all(&value.to_le_bytes())?;
            }
            w.write_all(&UNIT_CELL_RECORD_BYTES.to_le_bytes())?;
        }
        write_f32_fortran_block(&mut w, &xs)?;
        write_f32_fortran_block(&mut w, &ys)?;
        write_f32_fortran_block(&mut w, &zs)?;
        w.flush()?;
        self.frames_written += 1;
        Ok(())
    }

    /// Patch the frame count (`NSET`, byte offset 8) and the last step
    /// (`NSTEP`, byte offset 20) of the DCD header with the actual number of
    /// frames written, then flush the file.
    ///
    /// Call this once after all frames have been written.
    pub fn finalize(&self) -> io::Result<()> {
        let mut file = OpenOptions::new().write(true).open(&self.path)?;
        file.seek(SeekFrom::Start(NSET_OFFSET))?;
        file.write_all(&self.frames_written.to_le_bytes())?;
        file.seek(SeekFrom::Start(NSTEP_OFFSET))?;
        file.write_all(&self.last_step(self.frames_written).to_le_bytes())?;
        file.flush()
    }

    /// Last step number after `n_frames` frames.
    fn last_step(&self, n_frames: u32) -> u32 {
        WRITER_FIRST_STEP.saturating_add(n_frames.saturating_mul(WRITER_STEP_INTERVAL))
    }

    /// Write the DCD file header.
    ///
    /// `n_frames` is the preliminary frame count; it is updated by
    /// [`finalize`](DcdWriter::finalize).
    fn write_header(&mut self, n_frames: u32) -> io::Result<()> {
        let mut f = BufWriter::new(File::create(&self.path)?);

        // ── Record 1: 84 bytes, "CORD" and the 20 control words ───────────────
        let mut icntrl = [0u32; ICNTRL_WORDS];
        icntrl[IC_NSET] = n_frames;
        icntrl[IC_ISTART] = WRITER_FIRST_STEP;
        icntrl[IC_NSAVC] = WRITER_STEP_INTERVAL;
        icntrl[IC_NSTEP] = self.last_step(n_frames);
        // Words 4..8 are unused and word 8 (NAMNF, fixed atoms) stays zero.
        icntrl[IC_DELTA] = self.delta.to_bits();
        icntrl[IC_UNIT_CELL] = u32::from(self.has_unit_cell);
        icntrl[IC_VERSION] = CHARMM_VERSION;
        f.write_all(&HEADER_RECORD_BYTES.to_le_bytes())?;
        f.write_all(b"CORD")?;
        for word in icntrl {
            f.write_all(&word.to_le_bytes())?;
        }
        f.write_all(&HEADER_RECORD_BYTES.to_le_bytes())?;

        // ── Record 2: one 80-byte title ────────────────────────────────────────
        let title_payload: u32 = 4 + TITLE_BYTES; // n_titles + 1 title
        f.write_all(&title_payload.to_le_bytes())?;
        f.write_all(&1u32.to_le_bytes())?; // n_titles = 1
        let mut title_buf = [b' '; TITLE_BYTES as usize];
        let label = b"Created by oxiphysics DcdWriter";
        let copy_len = label.len().min(title_buf.len());
        title_buf[..copy_len].copy_from_slice(&label[..copy_len]);
        f.write_all(&title_buf)?;
        f.write_all(&title_payload.to_le_bytes())?;

        // ── Record 3: n_atoms ─────────────────────────────────────────────────
        f.write_all(&4u32.to_le_bytes())?; // payload = 4 bytes
        f.write_all(&self.n_atoms.to_le_bytes())?;
        f.write_all(&4u32.to_le_bytes())?;

        f.flush()
    }
}
