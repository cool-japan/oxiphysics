// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle trajectory I/O formats.
//!
//! Provides readers and writers for common particle simulation formats:
//! - **XYZ** multi-frame text trajectory
//! - **LAMMPS dump** (`ITEM: ATOMS id type x y z [vx vy vz]`)
//! - **GROMACS GRO** (nm units, box vectors)
//! - **Binary frame** format (compact f32 arrays with magic header)
//! - **DCD** CHARMM/NAMD/X-PLOR binary trajectories (reader with unit cells,
//!   writer), see [`dcd`]
//!
//! [`ParticleTrajectory`] supports frame indexing, slicing, sub-sampling,
//! and analysis (MSD, RMSF, centre-of-mass drift).

use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};

use crate::lammps_dump::{LammpsDumpFrame, LammpsDumpReader as LammpsFrameReader};

pub mod dcd;
#[cfg(test)]
mod dcd_tests;

pub use dcd::{DcdHeader, DcdReader, DcdTrajectory, DcdUnitCell, DcdWriter};

// ─────────────────────────────────────────────────────────────────────────────
// ParticleFrame
// ─────────────────────────────────────────────────────────────────────────────

/// A single snapshot of a particle system at a given time-step.
#[derive(Debug, Clone, Default)]
pub struct ParticleFrame {
    /// Integer time-step index.
    pub timestep: u64,
    /// Number of particles (redundant with `positions.len()` but convenient).
    pub n_particles: usize,
    /// Positions of each particle `[x, y, z]` in simulation units.
    pub positions: Vec<[f32; 3]>,
    /// Velocities of each particle `[vx, vy, vz]`, if present.
    pub velocities: Option<Vec<[f32; 3]>>,
    /// Integer type/species identifier for each particle.
    pub types: Vec<u8>,
}

impl ParticleFrame {
    /// Create a new frame with positions and types.
    ///
    /// `velocities` is set to `None`; fill it manually if needed.
    pub fn new(timestep: u64, positions: Vec<[f32; 3]>, types: Vec<u8>) -> Self {
        let n = positions.len();
        Self {
            timestep,
            n_particles: n,
            positions,
            velocities: None,
            types,
        }
    }

    /// Create a frame where all particles have type `0` and no velocities.
    pub fn from_positions(timestep: u64, positions: Vec<[f32; 3]>) -> Self {
        let n = positions.len();
        let types = vec![0u8; n];
        Self::new(timestep, positions, types)
    }

    /// Centre of mass (unweighted, equal masses assumed).
    ///
    /// Returns `[0.0; 3]` for an empty frame.
    pub fn center_of_mass(&self) -> [f32; 3] {
        if self.positions.is_empty() {
            return [0.0; 3];
        }
        let n = self.positions.len() as f32;
        let mut com = [0.0f32; 3];
        for p in &self.positions {
            com[0] += p[0];
            com[1] += p[1];
            com[2] += p[2];
        }
        [com[0] / n, com[1] / n, com[2] / n]
    }

    /// Axis-aligned bounding box `(min, max)`.  Returns `([0;3],[0;3])` if empty.
    pub fn bounding_box(&self) -> ([f32; 3], [f32; 3]) {
        if self.positions.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut lo = self.positions[0];
        let mut hi = self.positions[0];
        for &p in &self.positions[1..] {
            for d in 0..3 {
                if p[d] < lo[d] {
                    lo[d] = p[d];
                }
                if p[d] > hi[d] {
                    hi[d] = p[d];
                }
            }
        }
        (lo, hi)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// XyzWriter
// ─────────────────────────────────────────────────────────────────────────────

/// Write XYZ trajectory files (multi-frame, one frame per call).
///
/// Each frame is written as:
/// ```text
/// `N`
/// Timestep=`ts` `comment`
/// `type` `x` `y` `z`
/// ...
/// ```
#[derive(Debug, Clone)]
pub struct XyzWriter {
    /// Output path.
    pub path: std::path::PathBuf,
}

impl XyzWriter {
    /// Create a new `XyzWriter` targeting `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Overwrite the file with a single frame.
    pub fn write_frame(&self, frame: &ParticleFrame, comment: &str) -> io::Result<()> {
        let f = File::create(&self.path)?;
        let mut w = BufWriter::new(f);
        write_xyz_frame_to(&mut w, frame, comment)
    }

    /// Append a frame to the file (creates if not present).
    pub fn append_frame(&self, frame: &ParticleFrame, comment: &str) -> io::Result<()> {
        let f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        let mut w = BufWriter::new(f);
        write_xyz_frame_to(&mut w, frame, comment)
    }

    /// Write a complete trajectory (overwrites).
    pub fn write_trajectory(&self, traj: &ParticleTrajectory) -> io::Result<()> {
        let f = File::create(&self.path)?;
        let mut w = BufWriter::new(f);
        for frame in &traj.frames {
            let comment = format!("Timestep={}", frame.timestep);
            write_xyz_frame_to(&mut w, frame, &comment)?;
        }
        Ok(())
    }
}

/// Write one XYZ frame to any `Write` implementor.
fn write_xyz_frame_to<W: Write>(w: &mut W, frame: &ParticleFrame, comment: &str) -> io::Result<()> {
    writeln!(w, "{}", frame.n_particles)?;
    writeln!(w, "Timestep={} {}", frame.timestep, comment)?;
    for (i, p) in frame.positions.iter().enumerate() {
        let t = frame.types.get(i).copied().unwrap_or(0);
        writeln!(w, "{} {:.8} {:.8} {:.8}", t, p[0], p[1], p[2])?;
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// XyzReader
// ─────────────────────────────────────────────────────────────────────────────

/// Read XYZ trajectory files, returning one or all frames.
#[derive(Debug, Clone)]
pub struct XyzReader {
    /// Input path.
    pub path: std::path::PathBuf,
}

impl XyzReader {
    /// Create a new `XyzReader` for `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Read the first frame from the file.
    pub fn read_frame(&self) -> io::Result<ParticleFrame> {
        let f = File::open(&self.path)?;
        let mut r = BufReader::new(f);
        parse_xyz_frame(&mut r)
    }

    /// Read all frames from the file.
    pub fn read_all(&self) -> io::Result<Vec<ParticleFrame>> {
        let f = File::open(&self.path)?;
        let mut r = BufReader::new(f);
        let mut frames = Vec::new();
        loop {
            match parse_xyz_frame(&mut r) {
                Ok(fr) => frames.push(fr),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
        }
        Ok(frames)
    }
}

/// Parse one XYZ frame from a `BufRead` source.
fn parse_xyz_frame<R: BufRead>(r: &mut R) -> io::Result<ParticleFrame> {
    let mut line = String::new();

    // Line 1: atom count.
    r.read_line(&mut line)?;
    if line.is_empty() {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "EOF"));
    }
    let n: usize = line
        .trim()
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "expected atom count"))?;
    line.clear();

    // Line 2: comment — extract timestep if present.
    r.read_line(&mut line)?;
    let timestep = extract_timestep_from_comment(&line);
    line.clear();

    // Lines 3..(3+n): atom records.
    let mut positions = Vec::with_capacity(n);
    let mut types = Vec::with_capacity(n);
    for _ in 0..n {
        r.read_line(&mut line)?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "short atom line",
            ));
        }
        let t: u8 = parts[0].parse().unwrap_or(0);
        let x: f32 = parts[1]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "bad x"))?;
        let y: f32 = parts[2]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "bad y"))?;
        let z: f32 = parts[3]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "bad z"))?;
        positions.push([x, y, z]);
        types.push(t);
        line.clear();
    }

    Ok(ParticleFrame {
        timestep,
        n_particles: n,
        positions,
        velocities: None,
        types,
    })
}

/// Extract `Timestep=`value` from an XYZ comment line.  Returns 0 on failure.
fn extract_timestep_from_comment(comment: &str) -> u64 {
    for token in comment.split_whitespace() {
        if let Some(rest) = token.strip_prefix("Timestep=")
            && let Ok(v) = rest.parse::<u64>()
        {
            return v;
        }
    }
    0
}

// ─────────────────────────────────────────────────────────────────────────────
// LammpsDumpWriter
// ─────────────────────────────────────────────────────────────────────────────

/// Write LAMMPS dump files in `ITEM: ATOMS id type x y z \[vx vy vz\]` format.
#[derive(Debug, Clone)]
pub struct LammpsDumpWriter {
    /// Output path.
    pub path: std::path::PathBuf,
    /// Whether to write velocities when available.
    pub write_velocities: bool,
}

impl LammpsDumpWriter {
    /// Create a new `LammpsDumpWriter`.  Velocities are included if present.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self {
            path: path.into(),
            write_velocities: true,
        }
    }

    /// Append one frame to the dump file.
    ///
    /// `box_lo` and `box_hi` are the simulation box boundaries `\[x, y, z\]`.
    pub fn write_frame(
        &self,
        frame: &ParticleFrame,
        box_lo: [f32; 3],
        box_hi: [f32; 3],
    ) -> io::Result<()> {
        let f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        let mut w = BufWriter::new(f);
        write_lammps_dump_frame(&mut w, frame, box_lo, box_hi, self.write_velocities)
    }
}

/// Write one LAMMPS dump frame to any `Write` implementor.
fn write_lammps_dump_frame<W: Write>(
    w: &mut W,
    frame: &ParticleFrame,
    box_lo: [f32; 3],
    box_hi: [f32; 3],
    with_vel: bool,
) -> io::Result<()> {
    writeln!(w, "ITEM: TIMESTEP")?;
    writeln!(w, "{}", frame.timestep)?;
    writeln!(w, "ITEM: NUMBER OF ATOMS")?;
    writeln!(w, "{}", frame.n_particles)?;
    writeln!(w, "ITEM: BOX BOUNDS pp pp pp")?;
    for d in 0..3 {
        writeln!(w, "{:.6} {:.6}", box_lo[d], box_hi[d])?;
    }
    let has_vel = with_vel && frame.velocities.is_some();
    if has_vel {
        writeln!(w, "ITEM: ATOMS id type x y z vx vy vz")?;
    } else {
        writeln!(w, "ITEM: ATOMS id type x y z")?;
    }
    let vels = frame.velocities.as_deref();
    for (i, p) in frame.positions.iter().enumerate() {
        let t = frame.types.get(i).copied().unwrap_or(0);
        if has_vel {
            let v = vels.map(|vs| vs[i]).unwrap_or([0.0; 3]);
            writeln!(
                w,
                "{} {} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8}",
                i + 1,
                t,
                p[0],
                p[1],
                p[2],
                v[0],
                v[1],
                v[2]
            )?;
        } else {
            writeln!(w, "{} {} {:.8} {:.8} {:.8}", i + 1, t, p[0], p[1], p[2])?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// LammpsDumpReader
// ─────────────────────────────────────────────────────────────────────────────

/// Read LAMMPS dump files, handling `ITEM: ATOMS` headers and multiple timesteps.
///
/// The columns are found by the names in the `ITEM: ATOMS` header, scaled
/// coordinates (`xs ys zs`, the default of `dump atom`) are converted to
/// Cartesian ones with the box of each frame and the atoms of every frame are
/// sorted by id.  See [`LammpsDumpReader::read_all`] for the error cases.
#[derive(Debug, Clone)]
pub struct LammpsDumpReader {
    /// Input path.
    pub path: std::path::PathBuf,
}

impl LammpsDumpReader {
    /// Create a new `LammpsDumpReader` for `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Read all frames from the dump file.
    ///
    /// # Errors
    ///
    /// Returns the I/O error of opening or reading the file, and
    /// `io::ErrorKind::InvalidData` for a malformed frame: no complete
    /// position columns, a missing `id` or `type` column, a row that is
    /// shorter than the header or holds a value that is not a number, fewer
    /// rows than `NUMBER OF ATOMS`, or an atom type above 255.
    pub fn read_all(&self) -> io::Result<Vec<ParticleFrame>> {
        let f = File::open(&self.path)?;
        let r = BufReader::new(f);
        parse_lammps_dump_all(r)
    }
}

/// Parse all LAMMPS dump frames from a `BufRead` source.
///
/// The frames are parsed by [`crate::lammps_dump::LammpsDumpReader`]: the
/// columns of the `ITEM: ATOMS` table are located by the names in its header
/// (any order), the position columns may be `x y z`, `xu yu zu`, `xs ys zs`
/// (the default of `dump atom`) or `xsu ysu zsu`, and scaled coordinates are
/// converted to Cartesian ones with the box of the frame, including the
/// tilts of a triclinic box.  The atoms of every frame are returned sorted by
/// their id: parallel LAMMPS runs dump the atoms in an arbitrary order that
/// changes between frames, and trajectory analysis (MSD, RMSF, ...) compares
/// the particle with the same index across frames.
///
/// # Errors
///
/// `io::ErrorKind::InvalidData` when a frame is malformed: no complete
/// position columns, a missing `id` or `type` column, a row that is shorter
/// than the header or holds a value that is not a number, fewer rows than
/// `NUMBER OF ATOMS`, or an atom type above 255 (the type is stored in a `u8`).
/// The frame parser reports a read error inside a frame as text, so that is
/// `InvalidData` as well; a read error between frames keeps its kind.
fn parse_lammps_dump_all<R: BufRead>(mut r: R) -> io::Result<Vec<ParticleFrame>> {
    const FRAME_START: &str = "ITEM: TIMESTEP";
    let mut frames = Vec::new();
    let mut line = String::new();
    loop {
        // Everything in front of a frame is skipped: blank lines and the
        // `ITEM: UNITS` / `ITEM: TIME` blocks of `dump_modify units/time yes`.
        let at_frame = loop {
            line.clear();
            if r.read_line(&mut line)? == 0 {
                break false;
            }
            if line.trim().starts_with(FRAME_START) {
                break true;
            }
        };
        if !at_frame {
            return Ok(frames);
        }
        // Hand the frame parser the `ITEM: TIMESTEP` line that was just
        // consumed in front of the rest of the stream.
        let frame_start = io::Cursor::new(b"ITEM: TIMESTEP\n".as_slice());
        match LammpsFrameReader::read_frame(frame_start.chain(&mut r)) {
            Some(Ok(dump_frame)) => frames.push(particle_frame_from_dump(dump_frame)?),
            Some(Err(message)) => {
                return Err(io::Error::new(io::ErrorKind::InvalidData, message));
            }
            None => return Ok(frames),
        }
    }
}

/// One atom of a parsed LAMMPS dump frame.
struct DumpAtom {
    id: u32,
    atom_type: u32,
    position: [f64; 3],
    velocity: Option<[f64; 3]>,
}

/// Narrow a vector to `f32`.
fn to_f32(v: [f64; 3]) -> [f32; 3] {
    [v[0] as f32, v[1] as f32, v[2] as f32]
}

/// Convert a parsed LAMMPS dump frame: narrow to `f32`, store the types as
/// `u8` and order the atoms by id (ties keep the file order).
fn particle_frame_from_dump(frame: LammpsDumpFrame) -> io::Result<ParticleFrame> {
    let LammpsDumpFrame {
        timestep,
        atom_ids,
        types,
        positions,
        velocities,
        ..
    } = frame;
    let has_velocities = velocities.is_some();
    let mut velocity_iter = velocities.map(Vec::into_iter);
    let mut atoms: Vec<DumpAtom> = Vec::with_capacity(positions.len());
    for ((id, atom_type), position) in atom_ids.into_iter().zip(types).zip(positions) {
        atoms.push(DumpAtom {
            id,
            atom_type,
            position,
            velocity: velocity_iter.as_mut().and_then(Iterator::next),
        });
    }
    atoms.sort_by_key(|atom| atom.id);

    let mut out_positions = Vec::with_capacity(atoms.len());
    let mut out_types = Vec::with_capacity(atoms.len());
    let mut out_velocities: Vec<[f32; 3]> = Vec::new();
    for atom in atoms {
        let atom_type = u8::try_from(atom.atom_type).map_err(|_| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "LAMMPS dump: atom {} has type {}, above the maximum of {} \
                     that a particle frame can store",
                    atom.id,
                    atom.atom_type,
                    u8::MAX
                ),
            )
        })?;
        out_positions.push(to_f32(atom.position));
        out_types.push(atom_type);
        if has_velocities && let Some(velocity) = atom.velocity {
            out_velocities.push(to_f32(velocity));
        }
    }
    Ok(ParticleFrame {
        timestep,
        n_particles: out_positions.len(),
        positions: out_positions,
        velocities: (!out_velocities.is_empty()).then_some(out_velocities),
        types: out_types,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// ParticleTrajectory
// ─────────────────────────────────────────────────────────────────────────────

/// An in-memory trajectory: an ordered sequence of [`ParticleFrame`]s.
#[derive(Debug, Clone, Default)]
pub struct ParticleTrajectory {
    /// Ordered frames (earliest first).
    pub frames: Vec<ParticleFrame>,
}

impl ParticleTrajectory {
    /// Create an empty trajectory.
    pub fn new() -> Self {
        Self { frames: Vec::new() }
    }

    /// Number of frames.
    pub fn n_frames(&self) -> usize {
        self.frames.len()
    }

    /// Number of particles (from the first frame).  Returns 0 if empty.
    pub fn n_particles(&self) -> usize {
        self.frames.first().map(|f| f.n_particles).unwrap_or(0)
    }

    /// Get frame `i` by index.
    pub fn get(&self, i: usize) -> Option<&ParticleFrame> {
        self.frames.get(i)
    }

    /// Return a sub-trajectory containing frames in `\[start, end)`.
    pub fn slice(&self, start: usize, end: usize) -> Self {
        let end = end.min(self.frames.len());
        Self {
            frames: self.frames[start..end].to_vec(),
        }
    }

    /// Return every `step`-th frame (stride sub-sampling).
    pub fn subsample(&self, step: usize) -> Self {
        if step == 0 {
            return Self::new();
        }
        Self {
            frames: self.frames.iter().step_by(step).cloned().collect(),
        }
    }

    /// Append a frame to the trajectory.
    pub fn push(&mut self, frame: ParticleFrame) {
        self.frames.push(frame);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// TrajectoryStats
// ─────────────────────────────────────────────────────────────────────────────

/// Compute analysis statistics over a [`ParticleTrajectory`].
#[derive(Debug, Clone, Default)]
pub struct TrajectoryStats {
    /// Mean square displacement per lag step (length = n_frames).
    pub msd: Vec<f64>,
    /// Root mean square fluctuation per particle (length = n_particles).
    pub rmsf: Vec<f64>,
    /// Centre-of-mass position per frame (length = n_frames).
    pub com_drift: Vec<[f64; 3]>,
}

impl TrajectoryStats {
    /// Compute all statistics from a trajectory.
    pub fn compute(traj: &ParticleTrajectory) -> Self {
        Self {
            msd: compute_msd(traj),
            rmsf: compute_rmsf(traj),
            com_drift: compute_com_drift(traj),
        }
    }
}

/// Compute MSD: `msd\[lag\] = <|r(t+lag) - r(t)|^2>` averaged over all
/// particles and all valid origins.
fn compute_msd(traj: &ParticleTrajectory) -> Vec<f64> {
    let nf = traj.n_frames();
    let np = traj.n_particles();
    if nf == 0 || np == 0 {
        return vec![];
    }
    let mut result = vec![0.0f64; nf];
    for (lag, res) in result.iter_mut().enumerate() {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for t in 0..(nf - lag) {
            let f0 = &traj.frames[t];
            let f1 = &traj.frames[t + lag];
            let n = f0.positions.len().min(f1.positions.len());
            for (r0, r1) in f0.positions[..n].iter().zip(f1.positions[..n].iter()) {
                let dr2 =
                    (r1[0] - r0[0]).powi(2) + (r1[1] - r0[1]).powi(2) + (r1[2] - r0[2]).powi(2);
                sum += dr2 as f64;
                count += 1;
            }
        }
        *res = if count > 0 { sum / count as f64 } else { 0.0 };
    }
    result
}

/// Compute RMSF: root-mean-square fluctuation of each particle about its
/// time-averaged position.
fn compute_rmsf(traj: &ParticleTrajectory) -> Vec<f64> {
    let nf = traj.n_frames();
    let np = traj.n_particles();
    if nf == 0 || np == 0 {
        return vec![];
    }
    // Compute mean position for each particle.
    let mut mean = vec![[0.0f64; 3]; np];
    for frame in &traj.frames {
        for (i, p) in frame.positions.iter().enumerate().take(np) {
            mean[i][0] += p[0] as f64;
            mean[i][1] += p[1] as f64;
            mean[i][2] += p[2] as f64;
        }
    }
    for m in &mut mean {
        m[0] /= nf as f64;
        m[1] /= nf as f64;
        m[2] /= nf as f64;
    }
    // Compute variance.
    let mut var = vec![0.0f64; np];
    for frame in &traj.frames {
        for (i, p) in frame.positions.iter().enumerate().take(np) {
            let dx = p[0] as f64 - mean[i][0];
            let dy = p[1] as f64 - mean[i][1];
            let dz = p[2] as f64 - mean[i][2];
            var[i] += dx * dx + dy * dy + dz * dz;
        }
    }
    var.iter().map(|&v| (v / nf as f64).sqrt()).collect()
}

/// Compute centre-of-mass position for each frame.
fn compute_com_drift(traj: &ParticleTrajectory) -> Vec<[f64; 3]> {
    traj.frames
        .iter()
        .map(|f| {
            if f.positions.is_empty() {
                [0.0; 3]
            } else {
                let n = f.positions.len() as f64;
                let mut com = [0.0f64; 3];
                for p in &f.positions {
                    com[0] += p[0] as f64;
                    com[1] += p[1] as f64;
                    com[2] += p[2] as f64;
                }
                [com[0] / n, com[1] / n, com[2] / n]
            }
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────────────────────
// BinaryFrameWriter / BinaryFrameReader
// ─────────────────────────────────────────────────────────────────────────────

/// Magic bytes used to identify the binary frame format.
const BINARY_MAGIC: u32 = 0x4F58_4950; // "OXIP"

/// Write compact binary frames: `MAGIC(4) | timestep(8) | n(4) | x0y0z0... | has_vel(1) | vx0vy0vz0...`.
#[derive(Debug, Clone)]
pub struct BinaryFrameWriter {
    /// Output path.
    pub path: std::path::PathBuf,
}

impl BinaryFrameWriter {
    /// Create a new `BinaryFrameWriter` targeting `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Append one frame in binary format.
    pub fn write_frame(&self, frame: &ParticleFrame) -> io::Result<()> {
        let f = std::fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&self.path)?;
        let mut w = BufWriter::new(f);
        write_binary_frame(&mut w, frame)
    }
}

/// Write one binary frame to any `Write` implementor.
fn write_binary_frame<W: Write>(w: &mut W, frame: &ParticleFrame) -> io::Result<()> {
    w.write_all(&BINARY_MAGIC.to_le_bytes())?;
    w.write_all(&frame.timestep.to_le_bytes())?;
    w.write_all(&(frame.n_particles as u32).to_le_bytes())?;
    for p in &frame.positions {
        w.write_all(&p[0].to_le_bytes())?;
        w.write_all(&p[1].to_le_bytes())?;
        w.write_all(&p[2].to_le_bytes())?;
    }
    // types
    for &t in &frame.types {
        w.write_all(&[t])?;
    }
    // velocities flag + data
    let has_vel = frame.velocities.is_some();
    w.write_all(&[has_vel as u8])?;
    if let Some(vels) = &frame.velocities {
        for v in vels {
            w.write_all(&v[0].to_le_bytes())?;
            w.write_all(&v[1].to_le_bytes())?;
            w.write_all(&v[2].to_le_bytes())?;
        }
    }
    Ok(())
}

/// Read compact binary frames written by [`BinaryFrameWriter`].
#[derive(Debug, Clone)]
pub struct BinaryFrameReader {
    /// Input path.
    pub path: std::path::PathBuf,
}

impl BinaryFrameReader {
    /// Create a new `BinaryFrameReader` for `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Read all frames from the binary file.
    pub fn read_all(&self) -> io::Result<Vec<ParticleFrame>> {
        let f = File::open(&self.path)?;
        let mut r = BufReader::new(f);
        let mut frames = Vec::new();
        loop {
            match read_binary_frame(&mut r) {
                Ok(fr) => frames.push(fr),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
        }
        Ok(frames)
    }
}

/// Read one binary frame from a `Read` implementor.
fn read_binary_frame<R: Read>(r: &mut R) -> io::Result<ParticleFrame> {
    let magic = read_u32_le(r)?;
    if magic != BINARY_MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "bad magic"));
    }
    let timestep = read_u64_le(r)?;
    let n = read_u32_le(r)? as usize;

    let mut positions = Vec::with_capacity(n);
    for _ in 0..n {
        let x = read_f32_le(r)?;
        let y = read_f32_le(r)?;
        let z = read_f32_le(r)?;
        positions.push([x, y, z]);
    }
    let mut types = vec![0u8; n];
    r.read_exact(&mut types)?;

    let mut has_vel_buf = [0u8; 1];
    r.read_exact(&mut has_vel_buf)?;
    let velocities = if has_vel_buf[0] != 0 {
        let mut vels = Vec::with_capacity(n);
        for _ in 0..n {
            let vx = read_f32_le(r)?;
            let vy = read_f32_le(r)?;
            let vz = read_f32_le(r)?;
            vels.push([vx, vy, vz]);
        }
        Some(vels)
    } else {
        None
    };

    Ok(ParticleFrame {
        timestep,
        n_particles: n,
        positions,
        velocities,
        types,
    })
}

// ── binary I/O helpers ────────────────────────────────────────────────────────

fn read_u32_le<R: Read>(r: &mut R) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn read_u64_le<R: Read>(r: &mut R) -> io::Result<u64> {
    let mut buf = [0u8; 8];
    r.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn read_f32_le<R: Read>(r: &mut R) -> io::Result<f32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(f32::from_le_bytes(buf))
}

// ─────────────────────────────────────────────────────────────────────────────
// GroWriter / GroReader
// ─────────────────────────────────────────────────────────────────────────────

/// Write GROMACS `.gro` format snapshots (coordinates in nm).
///
/// Positions are stored as-is; the caller is responsible for unit conversion.
#[derive(Debug, Clone)]
pub struct GroWriter {
    /// Output path.
    pub path: std::path::PathBuf,
}

impl GroWriter {
    /// Create a new `GroWriter` targeting `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Write a single frame in GRO format with optional box vectors.
    ///
    /// `box_vecs` should be `\[lx, ly, lz\]` in nm.
    pub fn write_frame(&self, frame: &ParticleFrame, box_vecs: [f32; 3]) -> io::Result<()> {
        let f = File::create(&self.path)?;
        let mut w = BufWriter::new(f);
        writeln!(w, "Generated by oxiphysics-io  t={}", frame.timestep)?;
        writeln!(w, "{}", frame.n_particles)?;
        for (i, p) in frame.positions.iter().enumerate() {
            let t = frame.types.get(i).copied().unwrap_or(0);
            // GRO: %5d%-5s%5s%5d%8.3f%8.3f%8.3f
            writeln!(
                w,
                "{:5}UNK  {:>4}{:5}{:8.3}{:8.3}{:8.3}",
                (i + 1) % 100_000,
                t,
                (i + 1) % 100_000,
                p[0],
                p[1],
                p[2]
            )?;
        }
        writeln!(
            w,
            "{:.5}  {:.5}  {:.5}",
            box_vecs[0], box_vecs[1], box_vecs[2]
        )?;
        Ok(())
    }
}

/// Read GROMACS `.gro` format snapshots.
#[derive(Debug, Clone)]
pub struct GroReader {
    /// Input path.
    pub path: std::path::PathBuf,
}

impl GroReader {
    /// Create a new `GroReader` for `path`.
    pub fn new(path: impl Into<std::path::PathBuf>) -> Self {
        Self { path: path.into() }
    }

    /// Read the first frame from the GRO file.
    pub fn read_frame(&self) -> io::Result<ParticleFrame> {
        let f = File::open(&self.path)?;
        let mut r = BufReader::new(f);
        parse_gro_frame(&mut r)
    }
}

/// Parse one GRO frame from a `BufRead` source.
fn parse_gro_frame<R: BufRead>(r: &mut R) -> io::Result<ParticleFrame> {
    let mut line = String::new();
    // Title line.
    r.read_line(&mut line)?;
    if line.is_empty() {
        return Err(io::Error::new(io::ErrorKind::UnexpectedEof, "EOF"));
    }
    line.clear();
    // Atom count.
    r.read_line(&mut line)?;
    let n: usize = line
        .trim()
        .parse()
        .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "expected atom count"))?;
    line.clear();

    let mut positions = Vec::with_capacity(n);
    let mut types = Vec::with_capacity(n);
    for _ in 0..n {
        r.read_line(&mut line)?;
        // GRO fixed-width: fields at specific columns. Use whitespace split as fallback.
        let parts: Vec<&str> = line.split_whitespace().collect();
        // fields: resid resname atomname atomid x y z [vx vy vz]
        if parts.len() < 5 {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "short GRO line"));
        }
        // atom type from atomname (try parse as u8, else 0)
        let t: u8 = parts[1].parse().unwrap_or(0);
        let x: f32 = parts[3]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "bad x"))?;
        let y: f32 = parts[4]
            .parse()
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "bad y"))?;
        let z: f32 = if parts.len() > 5 {
            parts[5].parse().unwrap_or(0.0)
        } else {
            0.0
        };
        positions.push([x, y, z]);
        types.push(t);
        line.clear();
    }
    // Box line (ignore for now).
    Ok(ParticleFrame {
        timestep: 0,
        n_particles: n,
        positions,
        velocities: None,
        types,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    // ── helpers ───────────────────────────────────────────────────────────────

    fn make_frame(n: usize) -> ParticleFrame {
        let positions: Vec<[f32; 3]> = (0..n).map(|i| [i as f32, 0.0, 0.0]).collect();
        let types = vec![1u8; n];
        ParticleFrame::new(0, positions, types)
    }

    fn make_traj(n_frames: usize, n_particles: usize) -> ParticleTrajectory {
        let mut traj = ParticleTrajectory::new();
        for t in 0..n_frames {
            let positions: Vec<[f32; 3]> = (0..n_particles).map(|_| [t as f32, 0.0, 0.0]).collect();
            traj.push(ParticleFrame::new(
                t as u64,
                positions,
                vec![0u8; n_particles],
            ));
        }
        traj
    }

    // ── ParticleFrame ─────────────────────────────────────────────────────────

    #[test]
    fn test_frame_new_sets_n_particles() {
        let f = make_frame(5);
        assert_eq!(f.n_particles, 5);
    }

    #[test]
    fn test_frame_from_positions() {
        let pos = vec![[1.0f32, 0.0, 0.0], [2.0, 0.0, 0.0]];
        let f = ParticleFrame::from_positions(10, pos);
        assert_eq!(f.n_particles, 2);
        assert_eq!(f.timestep, 10);
        assert_eq!(f.types, vec![0u8, 0]);
    }

    #[test]
    fn test_frame_center_of_mass_empty() {
        let f = make_frame(0);
        assert_eq!(f.center_of_mass(), [0.0; 3]);
    }

    #[test]
    fn test_frame_center_of_mass_single() {
        let f = ParticleFrame::from_positions(0, vec![[2.0f32, 4.0, 6.0]]);
        let com = f.center_of_mass();
        assert!((com[0] - 2.0).abs() < 1e-6);
        assert!((com[1] - 4.0).abs() < 1e-6);
        assert!((com[2] - 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_frame_center_of_mass_symmetric() {
        let pos = vec![[-1.0f32, 0.0, 0.0], [1.0, 0.0, 0.0]];
        let f = ParticleFrame::from_positions(0, pos);
        let com = f.center_of_mass();
        assert!(com[0].abs() < 1e-6);
    }

    #[test]
    fn test_frame_bounding_box_empty() {
        let f = make_frame(0);
        let (lo, hi) = f.bounding_box();
        assert_eq!(lo, [0.0; 3]);
        assert_eq!(hi, [0.0; 3]);
    }

    #[test]
    fn test_frame_bounding_box_basic() {
        let pos = vec![[1.0f32, -1.0, 0.0], [-1.0, 2.0, 3.0]];
        let f = ParticleFrame::from_positions(0, pos);
        let (lo, hi) = f.bounding_box();
        assert!((lo[0] - (-1.0)).abs() < 1e-6);
        assert!((hi[1] - 2.0).abs() < 1e-6);
        assert!((hi[2] - 3.0).abs() < 1e-6);
    }

    // ── XYZ round-trip ────────────────────────────────────────────────────────

    #[test]
    fn test_xyz_write_frame_internal() {
        let frame = make_frame(3);
        let mut buf = Vec::new();
        write_xyz_frame_to(&mut buf, &frame, "test").unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.starts_with("3\n"));
        assert!(s.contains("test"));
    }

    #[test]
    fn test_xyz_roundtrip_file() {
        let path = std::env::temp_dir().join("oxiphysics_pf_xyz_rt.xyz");
        let _ = std::fs::remove_file(&path);
        let pos = vec![[1.0f32, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let frame = ParticleFrame::new(42, pos.clone(), vec![1u8, 2]);
        XyzWriter::new(&path).write_frame(&frame, "test").unwrap();
        let loaded = XyzReader::new(&path).read_frame().unwrap();
        assert_eq!(loaded.n_particles, 2);
        assert_eq!(loaded.timestep, 42);
        assert!((loaded.positions[0][0] - 1.0).abs() < 1e-4);
        assert!((loaded.positions[1][2] - 6.0).abs() < 1e-4);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_xyz_append_and_read_all() {
        let path = std::env::temp_dir().join("oxiphysics_pf_xyz_all.xyz");
        let _ = std::fs::remove_file(&path);
        let w = XyzWriter::new(&path);
        for ts in 0..4u64 {
            let f = ParticleFrame::new(ts, vec![[ts as f32, 0.0, 0.0]], vec![0u8]);
            w.append_frame(&f, "").unwrap();
        }
        let frames = XyzReader::new(&path).read_all().unwrap();
        assert_eq!(frames.len(), 4);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_xyz_write_trajectory() {
        let path = std::env::temp_dir().join("oxiphysics_pf_xyz_traj.xyz");
        let _ = std::fs::remove_file(&path);
        let mut traj = ParticleTrajectory::new();
        for ts in 0..3u64 {
            traj.push(ParticleFrame::new(ts, vec![[0.0f32; 3]], vec![0u8]));
        }
        XyzWriter::new(&path).write_trajectory(&traj).unwrap();
        let frames = XyzReader::new(&path).read_all().unwrap();
        assert_eq!(frames.len(), 3);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_extract_timestep_from_comment() {
        assert_eq!(extract_timestep_from_comment("Timestep=100 other"), 100);
        assert_eq!(extract_timestep_from_comment("no timestep here"), 0);
    }

    // ── LAMMPS dump ───────────────────────────────────────────────────────────

    #[test]
    fn test_lammps_dump_write_internal() {
        let frame = make_frame(2);
        let mut buf = Vec::new();
        write_lammps_dump_frame(&mut buf, &frame, [0.0; 3], [10.0; 3], false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("ITEM: TIMESTEP"));
        assert!(s.contains("ITEM: ATOMS id type x y z"));
    }

    #[test]
    fn test_lammps_dump_write_with_velocities() {
        let mut frame = make_frame(2);
        frame.velocities = Some(vec![[1.0f32, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let mut buf = Vec::new();
        write_lammps_dump_frame(&mut buf, &frame, [0.0; 3], [10.0; 3], true).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("vx vy vz"));
    }

    #[test]
    fn test_lammps_dump_roundtrip_file() {
        let path = std::env::temp_dir().join("oxiphysics_pf_lmp.lammpstrj");
        let _ = std::fs::remove_file(&path);
        let w = LammpsDumpWriter::new(&path);
        for ts in 0..3u64 {
            let f = ParticleFrame::new(ts, vec![[ts as f32, 0.0, 0.0]], vec![1u8]);
            w.write_frame(&f, [0.0; 3], [10.0; 3]).unwrap();
        }
        let frames = LammpsDumpReader::new(&path).read_all().unwrap();
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[1].timestep, 1);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_lammps_dump_parse_inline() {
        let dump = "\
ITEM: TIMESTEP\n\
5\n\
ITEM: NUMBER OF ATOMS\n\
2\n\
ITEM: BOX BOUNDS pp pp pp\n\
0.0 10.0\n\
0.0 10.0\n\
0.0 10.0\n\
ITEM: ATOMS id type x y z\n\
1 1 1.0 2.0 3.0\n\
2 2 4.0 5.0 6.0\n";
        let r = BufReader::new(Cursor::new(dump));
        let frames = parse_lammps_dump_all(r).unwrap();
        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].timestep, 5);
        assert_eq!(frames[0].n_particles, 2);
    }

    // ── LAMMPS dump: column mapping, scaled coordinates, ids (F2253) ──────────

    fn parse_dump(dump: &str) -> io::Result<Vec<ParticleFrame>> {
        parse_lammps_dump_all(BufReader::new(Cursor::new(dump.to_owned())))
    }

    /// Assert that a position equals `expected` (every value used in these
    /// tests is exactly representable in `f32`).
    fn assert_position(actual: [f32; 3], expected: [f32; 3]) {
        assert_eq!(actual, expected, "position {actual:?} != {expected:?}");
    }

    #[test]
    fn lammps_dump_atom_style_scaled_columns_are_unscaled_and_sorted_by_id() {
        // `dump atom` writes scaled coordinates; ids arrive shuffled, as they do
        // from a multi-process run.  The box origin is not zero.
        let dump = "\
ITEM: TIMESTEP
100
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
-5.0 5.0
0.0 20.0
10.0 14.0
ITEM: ATOMS id type xs ys zs
3 2 0.25 0.5 1.0
1 1 0.0 0.0 0.0
2 1 0.5 0.25 0.5
";
        let frames = parse_dump(dump).unwrap();
        assert_eq!(frames.len(), 1);
        let f = &frames[0];
        assert_eq!(f.timestep, 100);
        assert_eq!(f.n_particles, 3);
        // id 1, 2, 3: lo + s * (hi - lo) per axis
        assert_position(f.positions[0], [-5.0, 0.0, 10.0]);
        assert_position(f.positions[1], [0.0, 5.0, 12.0]);
        assert_position(f.positions[2], [-2.5, 10.0, 14.0]);
        assert_eq!(f.types, vec![1, 1, 2]);
        assert!(f.velocities.is_none());
    }

    #[test]
    fn lammps_dump_unwrapped_and_scaled_unwrapped_columns_are_resolved() {
        let unwrapped = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type xu yu zu
1 1 12.5 -3.0 7.0
";
        let f = &parse_dump(unwrapped).unwrap()[0];
        assert_position(f.positions[0], [12.5, -3.0, 7.0]);

        // Scaled-unwrapped coordinates may leave [0, 1) (the image flag).
        let scaled_unwrapped = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
2.0 6.0
0.0 8.0
-1.0 1.0
ITEM: ATOMS id type xsu ysu zsu
1 1 1.25 -0.5 2.0
";
        let f = &parse_dump(scaled_unwrapped).unwrap()[0];
        assert_position(f.positions[0], [7.0, -4.0, 3.0]);
    }

    #[test]
    fn lammps_dump_columns_are_found_by_name_not_by_position() {
        let dump = "\
ITEM: TIMESTEP
7
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS z vz y vy id x vx type
3.0 0.3 2.0 0.2 2 1.0 0.1 5
6.0 0.6 5.0 0.5 1 4.0 0.4 4
";
        let f = &parse_dump(dump).unwrap()[0];
        assert_position(f.positions[0], [4.0, 5.0, 6.0]);
        assert_position(f.positions[1], [1.0, 2.0, 3.0]);
        assert_eq!(f.types, vec![4, 5]);
        let v = f.velocities.as_ref().expect("vx vy vz are present");
        assert_eq!(v[0], [0.4, 0.5, 0.6]);
        assert_eq!(v[1], [0.1, 0.2, 0.3]);
    }

    #[test]
    fn lammps_dump_triclinic_box_unscales_with_the_tilt_factors() {
        // Restricted triclinic box: xy = 2, xz = -1, yz = 0.5.  The bounds are
        // those of the bounding box, so
        //   xlo = -1 - min(0, 2, -1, 1) = 0,  xhi = 12 - max(0, 2, -1, 1) = 10,
        //   ylo =  0 - min(0, 0.5)      = 0,  yhi = 8.5 - max(0, 0.5)     = 8.
        let dump = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS xy xz yz pp pp pp
-1.0 12.0 2.0
0.0 8.5 -1.0
0.0 6.0 0.5
ITEM: ATOMS id type xs ys zs
2 1 1.0 1.0 1.0
1 1 0.5 0.25 0.5
";
        let f = &parse_dump(dump).unwrap()[0];
        // x = xlo + xs*lx + ys*xy + zs*xz, y = ylo + ys*ly + zs*yz, z = zlo + zs*lz
        assert_position(f.positions[0], [5.0, 2.25, 3.0]);
        assert_position(f.positions[1], [11.0, 8.5, 6.0]);
    }

    #[test]
    fn lammps_dump_atoms_are_matched_by_id_across_frames() {
        // The second frame lists the same two atoms in the opposite order.
        let dump = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 1.0 0.0 0.0
2 2 9.0 0.0 0.0
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
2 2 8.0 0.0 0.0
1 1 2.0 0.0 0.0
";
        let frames = parse_dump(dump).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[1].timestep, 10);
        for f in &frames {
            assert_eq!(f.types, vec![1, 2]);
        }
        assert_position(frames[0].positions[0], [1.0, 0.0, 0.0]);
        assert_position(frames[1].positions[0], [2.0, 0.0, 0.0]);
        assert_position(frames[0].positions[1], [9.0, 0.0, 0.0]);
        assert_position(frames[1].positions[1], [8.0, 0.0, 0.0]);

        let mut traj = ParticleTrajectory::new();
        for f in frames {
            traj.push(f);
        }
        // Atom 1 moves by +1, atom 2 by -1: mean squared displacement 1.
        let msd = TrajectoryStats::compute(&traj).msd;
        assert!((msd[1] - 1.0).abs() < 1e-9, "msd = {msd:?}");
    }

    #[test]
    fn lammps_dump_units_and_time_items_in_front_of_a_frame_are_skipped() {
        let dump = "\
ITEM: UNITS
metal
ITEM: TIME
0.5
ITEM: TIMESTEP
500
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0.0 1.0
0.0 1.0
0.0 1.0
ITEM: ATOMS id type x y z
1 1 0.5 0.5 0.5

ITEM: UNITS
metal
ITEM: TIME
1.0
ITEM: TIMESTEP
1000
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0.0 1.0
0.0 1.0
0.0 1.0
ITEM: ATOMS id type x y z
1 1 0.25 0.5 0.75
";
        let frames = parse_dump(dump).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].timestep, 500);
        assert_eq!(frames[1].timestep, 1000);
        assert_position(frames[1].positions[0], [0.25, 0.5, 0.75]);
    }

    #[test]
    fn lammps_dump_without_position_columns_is_invalid_data() {
        let dump = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
1
ITEM: BOX BOUNDS pp pp pp
0.0 1.0
0.0 1.0
0.0 1.0
ITEM: ATOMS id type vx vy vz
1 1 0.0 0.0 0.0
";
        let err = parse_dump(dump).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("position"), "{err}");

        // An incomplete family (x and y but no z) is not a position either.
        let partial = dump.replace("id type vx vy vz", "id type x y");
        let err = parse_dump(&partial).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn lammps_dump_unparsable_value_is_invalid_data() {
        let dump = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0.0 1.0
0.0 1.0
0.0 1.0
ITEM: ATOMS id type x y z
1 1 0.1 0.2 0.3
2 1 0.4 abc 0.6
";
        let err = parse_dump(dump).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("abc"), "{err}");

        // A row with fewer values than the header names columns.
        let short = dump.replace("2 1 0.4 abc 0.6", "2 1 0.4");
        let err = parse_dump(&short).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn lammps_dump_atom_type_above_255_is_invalid_data() {
        let build = |atom_type: u32| {
            format!(
                "ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS pp pp pp\n\
                 0.0 1.0\n0.0 1.0\n0.0 1.0\nITEM: ATOMS id type x y z\n1 {atom_type} 0.5 0.5 0.5\n"
            )
        };
        let f = &parse_dump(&build(255)).unwrap()[0];
        assert_eq!(f.types, vec![255]);
        let err = parse_dump(&build(256)).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("256"), "{err}");
    }

    #[test]
    fn lammps_dump_frame_with_fewer_rows_than_declared_is_an_error() {
        let dump = "\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 1.0
0.0 1.0
0.0 1.0
ITEM: ATOMS id type x y z
1 1 0.1 0.2 0.3
2 1 0.4 0.5 0.6
";
        let err = parse_dump(dump).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidData);
        assert!(err.to_string().contains("EOF"), "{err}");
    }

    // ── ParticleTrajectory ────────────────────────────────────────────────────

    #[test]
    fn test_trajectory_push_and_get() {
        let mut traj = ParticleTrajectory::new();
        traj.push(make_frame(3));
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(traj.n_particles(), 3);
        assert!(traj.get(0).is_some());
        assert!(traj.get(1).is_none());
    }

    #[test]
    fn test_trajectory_slice() {
        let traj = make_traj(10, 2);
        let sub = traj.slice(2, 5);
        assert_eq!(sub.n_frames(), 3);
    }

    #[test]
    fn test_trajectory_slice_beyond_end() {
        let traj = make_traj(5, 2);
        let sub = traj.slice(3, 100);
        assert_eq!(sub.n_frames(), 2);
    }

    #[test]
    fn test_trajectory_subsample() {
        let traj = make_traj(10, 1);
        let sub = traj.subsample(2);
        assert_eq!(sub.n_frames(), 5);
    }

    #[test]
    fn test_trajectory_subsample_zero_step() {
        let traj = make_traj(10, 1);
        let sub = traj.subsample(0);
        assert_eq!(sub.n_frames(), 0);
    }

    // ── TrajectoryStats ───────────────────────────────────────────────────────

    #[test]
    fn test_stats_msd_stationary() {
        let mut traj = ParticleTrajectory::new();
        for _ in 0..5 {
            traj.push(ParticleFrame::from_positions(0, vec![[0.0f32; 3]]));
        }
        let stats = TrajectoryStats::compute(&traj);
        for &m in &stats.msd {
            assert!(m.abs() < 1e-8, "msd = {m}");
        }
    }

    #[test]
    fn test_stats_msd_drift() {
        let mut traj = ParticleTrajectory::new();
        for i in 0..5usize {
            traj.push(ParticleFrame::from_positions(
                i as u64,
                vec![[i as f32, 0.0, 0.0]],
            ));
        }
        let stats = TrajectoryStats::compute(&traj);
        // msd[1] should be 1, msd[2] = 4, etc.
        assert!(
            (stats.msd[1] - 1.0).abs() < 1e-6,
            "msd[1] = {}",
            stats.msd[1]
        );
    }

    #[test]
    fn test_stats_rmsf_stationary() {
        let mut traj = ParticleTrajectory::new();
        for _ in 0..10 {
            traj.push(ParticleFrame::from_positions(0, vec![[1.0f32, 2.0, 3.0]]));
        }
        let stats = TrajectoryStats::compute(&traj);
        assert!(stats.rmsf[0].abs() < 1e-8);
    }

    #[test]
    fn test_stats_com_drift_single_frame() {
        let mut traj = ParticleTrajectory::new();
        traj.push(ParticleFrame::from_positions(0, vec![[2.0f32, 4.0, 6.0]]));
        let stats = TrajectoryStats::compute(&traj);
        let com = stats.com_drift[0];
        assert!((com[0] - 2.0).abs() < 1e-6);
    }

    #[test]
    fn test_stats_empty_trajectory() {
        let traj = ParticleTrajectory::new();
        let stats = TrajectoryStats::compute(&traj);
        assert!(stats.msd.is_empty());
        assert!(stats.rmsf.is_empty());
        assert!(stats.com_drift.is_empty());
    }

    // ── BinaryFrameWriter / BinaryFrameReader ─────────────────────────────────

    #[test]
    fn test_binary_frame_roundtrip_internal() {
        let frame = make_frame(4);
        let mut buf = Vec::new();
        write_binary_frame(&mut buf, &frame).unwrap();
        let mut cursor = Cursor::new(buf);
        let loaded = read_binary_frame(&mut cursor).unwrap();
        assert_eq!(loaded.n_particles, 4);
        for i in 0..4 {
            assert!((loaded.positions[i][0] - i as f32).abs() < 1e-6);
        }
    }

    #[test]
    fn test_binary_frame_with_velocities() {
        let mut frame = make_frame(2);
        frame.velocities = Some(vec![[1.0f32, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let mut buf = Vec::new();
        write_binary_frame(&mut buf, &frame).unwrap();
        let mut cursor = Cursor::new(buf);
        let loaded = read_binary_frame(&mut cursor).unwrap();
        assert!(loaded.velocities.is_some());
        let vels = loaded.velocities.unwrap();
        assert!((vels[0][0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_binary_frame_file_roundtrip() {
        let path = std::env::temp_dir().join("oxiphysics_pf_bin_rt.bin");
        let _ = std::fs::remove_file(&path);
        let w = BinaryFrameWriter::new(&path);
        for ts in 0..3u64 {
            let f = ParticleFrame::new(ts, vec![[ts as f32, 0.0, 0.0]], vec![0u8]);
            w.write_frame(&f).unwrap();
        }
        let frames = BinaryFrameReader::new(&path).read_all().unwrap();
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[2].timestep, 2);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_binary_frame_bad_magic() {
        let buf = vec![0xFFu8, 0xFF, 0xFF, 0xFF];
        let mut cursor = Cursor::new(buf);
        let result = read_binary_frame(&mut cursor);
        assert!(result.is_err());
    }

    // ── GroWriter / GroReader ─────────────────────────────────────────────────

    #[test]
    fn test_gro_write_creates_file() {
        let path = std::env::temp_dir().join("oxiphysics_pf_gro.gro");
        let frame = make_frame(3);
        GroWriter::new(&path)
            .write_frame(&frame, [10.0; 3])
            .unwrap();
        assert!(path.exists());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_gro_write_contains_box() {
        let path = std::env::temp_dir().join("oxiphysics_pf_gro_box.gro");
        let frame = make_frame(1);
        GroWriter::new(&path)
            .write_frame(&frame, [5.0f32, 5.0, 5.0])
            .unwrap();
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.contains("5.00000"));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_gro_roundtrip_positions() {
        let path = std::env::temp_dir().join("oxiphysics_pf_gro_rt.gro");
        let pos = vec![[1.5f32, 2.5, 3.5], [4.5, 5.5, 6.5]];
        let frame = ParticleFrame::new(0, pos.clone(), vec![0u8, 0]);
        GroWriter::new(&path)
            .write_frame(&frame, [10.0; 3])
            .unwrap();
        let loaded = GroReader::new(&path).read_frame().unwrap();
        assert_eq!(loaded.n_particles, 2);
        // GRO uses %.3f precision for coordinates
        assert!((loaded.positions[0][0] - 1.5).abs() < 0.01);
        let _ = std::fs::remove_file(&path);
    }
}
