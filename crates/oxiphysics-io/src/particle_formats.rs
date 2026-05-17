#![allow(clippy::needless_range_loop)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle trajectory I/O formats.
//!
//! Provides readers and writers for common particle simulation formats:
//! - **XYZ** multi-frame text trajectory
//! - **LAMMPS dump** (`ITEM: ATOMS id type x y z [vx vy vz]`)
//! - **GROMACS GRO** (nm units, box vectors)
//! - **Binary frame** format (compact f32 arrays with magic header)
//! - **DCD** stub reader (header parsing + coordinate frames)
//!
//! [`ParticleTrajectory`] supports frame indexing, slicing, sub-sampling,
//! and analysis (MSD, RMSF, centre-of-mass drift).

#![allow(dead_code)]

use std::fs::{File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

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
    pub fn read_all(&self) -> io::Result<Vec<ParticleFrame>> {
        let f = File::open(&self.path)?;
        let r = BufReader::new(f);
        parse_lammps_dump_all(r)
    }
}

/// Parse all LAMMPS dump frames from a `BufRead` source.
fn parse_lammps_dump_all<R: BufRead>(mut r: R) -> io::Result<Vec<ParticleFrame>> {
    let mut frames = Vec::new();
    let mut lines = Vec::<String>::new();
    // Collect all lines.
    {
        let mut line = String::new();
        while r.read_line(&mut line)? > 0 {
            lines.push(line.trim_end().to_string());
            line.clear();
        }
    }

    let mut i = 0;
    while i < lines.len() {
        if lines[i].starts_with("ITEM: TIMESTEP") {
            i += 1;
            let timestep: u64 = lines
                .get(i)
                .and_then(|l| l.trim().parse().ok())
                .unwrap_or(0);
            i += 1; // "ITEM: NUMBER OF ATOMS"
            i += 1;
            let n_atoms: usize = lines
                .get(i)
                .and_then(|l| l.trim().parse().ok())
                .unwrap_or(0);
            i += 1; // "ITEM: BOX BOUNDS ..."
            i += 4; // skip 3 bound lines + the item line itself

            // "ITEM: ATOMS ..."
            let header = lines.get(i).cloned().unwrap_or_default();
            let cols: Vec<&str> = header.split_whitespace().collect();
            // Find column indices for id, type, x, y, z, vx, vy, vz.
            let idx_id = col_idx(&cols, "id");
            let idx_type = col_idx(&cols, "type");
            let idx_x = col_idx(&cols, "x");
            let idx_y = col_idx(&cols, "y");
            let idx_z = col_idx(&cols, "z");
            let idx_vx = col_idx(&cols, "vx");
            let idx_vy = col_idx(&cols, "vy");
            let idx_vz = col_idx(&cols, "vz");
            let has_vel = idx_vx.is_some() && idx_vy.is_some() && idx_vz.is_some();
            // Skip ITEM: ATOMS line — offsets start at 2 after ITEM: ATOMS header.
            // Columns in header are offset by 2 ("ITEM:" "ATOMS" are first two).
            let data_offset = 2usize;
            i += 1;

            let mut positions = Vec::with_capacity(n_atoms);
            let mut types = Vec::with_capacity(n_atoms);
            let mut velocities: Vec<[f32; 3]> = Vec::with_capacity(n_atoms);

            for _ in 0..n_atoms {
                if i >= lines.len() {
                    break;
                }
                let parts: Vec<&str> = lines[i].split_whitespace().collect();

                let get = |opt_col: Option<usize>| -> f32 {
                    opt_col
                        .and_then(|c| {
                            let real_c = c.saturating_sub(data_offset);
                            parts.get(real_c).and_then(|s| s.parse().ok())
                        })
                        .unwrap_or(0.0)
                };
                let get_u8 = |opt_col: Option<usize>| -> u8 {
                    opt_col
                        .and_then(|c| {
                            let real_c = c.saturating_sub(data_offset);
                            parts.get(real_c).and_then(|s| s.parse::<u8>().ok())
                        })
                        .unwrap_or(0)
                };
                let _ = idx_id; // id unused but parsed
                let t = get_u8(idx_type);
                let x = get(idx_x);
                let y = get(idx_y);
                let z = get(idx_z);
                positions.push([x, y, z]);
                types.push(t);
                if has_vel {
                    let vx = get(idx_vx);
                    let vy = get(idx_vy);
                    let vz = get(idx_vz);
                    velocities.push([vx, vy, vz]);
                }
                i += 1;
            }

            let vel_opt = if has_vel && !velocities.is_empty() {
                Some(velocities)
            } else {
                None
            };

            frames.push(ParticleFrame {
                timestep,
                n_particles: positions.len(),
                positions,
                velocities: vel_opt,
                types,
            });
        } else {
            i += 1;
        }
    }
    Ok(frames)
}

/// Find the index of a column name in the LAMMPS ITEM:ATOMS header tokens.
/// The header starts with "ITEM:" "ATOMS", so data columns start at index 2.
fn col_idx(cols: &[&str], name: &str) -> Option<usize> {
    cols.iter().position(|&c| c == name)
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
    for lag in 0..nf {
        let mut sum = 0.0f64;
        let mut count = 0usize;
        for t in 0..(nf - lag) {
            let f0 = &traj.frames[t];
            let f1 = &traj.frames[t + lag];
            let n = f0.positions.len().min(f1.positions.len());
            for i in 0..n {
                let r0 = f0.positions[i];
                let r1 = f1.positions[i];
                let dr2 =
                    (r1[0] - r0[0]).powi(2) + (r1[1] - r0[1]).powi(2) + (r1[2] - r0[2]).powi(2);
                sum += dr2 as f64;
                count += 1;
            }
        }
        result[lag] = if count > 0 { sum / count as f64 } else { 0.0 };
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
// DcdReader / DcdWriter
// ─────────────────────────────────────────────────────────────────────────────

/// Header information parsed from a DCD file.
#[derive(Debug, Clone, Default)]
pub struct DcdHeader {
    /// Number of frames stored in the file.
    pub n_frames: u32,
    /// First timestep number.
    pub first_step: u32,
    /// Number of steps between saved frames.
    pub step_interval: u32,
    /// Number of atoms.
    pub n_atoms: u32,
    /// Title strings from the header.
    pub titles: Vec<String>,
}

/// Stub DCD binary reader: parses the header and reads coordinate frames.
///
/// Supports the standard CHARMM/NAMD little-endian DCD format.
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
    pub fn read_header(&self) -> io::Result<DcdHeader> {
        let f = File::open(&self.path)?;
        let mut r = BufReader::new(f);
        parse_dcd_header(&mut r)
    }

    /// Read all coordinate frames.  Returns `(header, frames)`.
    pub fn read_all(&self) -> io::Result<(DcdHeader, Vec<ParticleFrame>)> {
        let f = File::open(&self.path)?;
        let mut r = BufReader::new(f);
        let header = parse_dcd_header(&mut r)?;
        let n = header.n_atoms as usize;
        let mut frames = Vec::with_capacity(header.n_frames as usize);
        for step in 0..header.n_frames {
            match parse_dcd_frame(&mut r, n, step as u64) {
                Ok(fr) => frames.push(fr),
                Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => break,
                Err(e) => return Err(e),
            }
        }
        Ok((header, frames))
    }
}

/// Parse DCD header from a `Read + Seek` source.
fn parse_dcd_header<R: Read + Seek>(r: &mut R) -> io::Result<DcdHeader> {
    // Block 1: 4-byte Fortran record length, "CORD", 9 ints, delta, 10 ints.
    let rec_len = read_u32_le(r)?;
    if rec_len < 4 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "DCD block1 too short",
        ));
    }
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    // Expect "CORD" or "VELO"
    let n_frames = read_u32_le(r)?;
    let first_step = read_u32_le(r)?;
    let step_interval = read_u32_le(r)?;
    // Skip remaining ints + delta + padding to end of block.
    let remaining = rec_len as usize - 4 - 4 * 3;
    r.seek(SeekFrom::Current(remaining as i64))?;
    let _end_len = read_u32_le(r)?; // closing Fortran record

    // Block 2: title block.
    let title_block_len = read_u32_le(r)? as usize;
    let n_titles = read_u32_le(r)? as usize;
    let mut titles = Vec::new();
    for _ in 0..n_titles {
        let mut buf = vec![0u8; 80];
        r.read_exact(&mut buf)?;
        titles.push(String::from_utf8_lossy(&buf).trim_end().to_string());
    }
    let consumed = 4 + n_titles * 80;
    if consumed < title_block_len {
        r.seek(SeekFrom::Current((title_block_len - consumed) as i64))?;
    }
    let _end_title = read_u32_le(r)?;

    // Block 3: n_atoms.
    let _len3 = read_u32_le(r)?;
    let n_atoms = read_u32_le(r)?;
    let _end3 = read_u32_le(r)?;

    Ok(DcdHeader {
        n_frames,
        first_step,
        step_interval,
        n_atoms,
        titles,
    })
}

/// Parse one DCD coordinate frame.
fn parse_dcd_frame<R: Read>(r: &mut R, n_atoms: usize, step: u64) -> io::Result<ParticleFrame> {
    let read_f32_block = |r: &mut R| -> io::Result<Vec<f32>> {
        let len = read_u32_le(r)? as usize;
        let count = len / 4;
        let mut vals = Vec::with_capacity(count);
        for _ in 0..count {
            vals.push(read_f32_le(r)?);
        }
        let _end = read_u32_le(r)?;
        Ok(vals)
    };

    let xs = read_f32_block(r)?;
    let ys = read_f32_block(r)?;
    let zs = read_f32_block(r)?;

    let n = n_atoms.min(xs.len()).min(ys.len()).min(zs.len());
    let positions: Vec<[f32; 3]> = (0..n).map(|i| [xs[i], ys[i], zs[i]]).collect();
    let types = vec![0u8; n];

    Ok(ParticleFrame {
        timestep: step,
        n_particles: n,
        positions,
        velocities: None,
        types,
    })
}

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
/// # Usage
///
/// ```rust,no_run
/// use oxiphysics_io::DcdWriter;
///
/// let mut w = DcdWriter::new("/tmp/traj.dcd", 10, 2.0e-3).unwrap();
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
    frames_written: u32,
}

impl DcdWriter {
    /// Create a new `DcdWriter`, truncating / creating the output file and
    /// writing an initial DCD header.
    pub fn new(path: impl Into<std::path::PathBuf>, n_atoms: u32, delta: f32) -> io::Result<Self> {
        let path = path.into();
        let mut writer = Self {
            path,
            n_atoms,
            delta,
            frames_written: 0,
        };
        writer.write_header(0)?;
        Ok(writer)
    }

    /// Append one coordinate frame to the DCD file.
    ///
    /// The frame must contain at least `n_atoms` positions; extra positions are
    /// silently ignored, missing ones are padded with `0.0`.
    pub fn write_frame(&mut self, frame: &ParticleFrame) -> io::Result<()> {
        use std::io::Write as _;
        let n = self.n_atoms as usize;

        let mut xs = Vec::with_capacity(n);
        let mut ys = Vec::with_capacity(n);
        let mut zs = Vec::with_capacity(n);
        for i in 0..n {
            if i < frame.positions.len() {
                xs.push(frame.positions[i][0]);
                ys.push(frame.positions[i][1]);
                zs.push(frame.positions[i][2]);
            } else {
                xs.push(0.0_f32);
                ys.push(0.0_f32);
                zs.push(0.0_f32);
            }
        }

        let file = OpenOptions::new().append(true).open(&self.path)?;
        let mut w = BufWriter::new(file);
        write_f32_fortran_block(&mut w, &xs)?;
        write_f32_fortran_block(&mut w, &ys)?;
        write_f32_fortran_block(&mut w, &zs)?;
        w.flush()?;
        self.frames_written += 1;
        Ok(())
    }

    /// Patch the frame count stored in the DCD header (byte offset 8) with the
    /// actual number of frames written, then flush the file.
    ///
    /// Call this once after all frames have been written.
    pub fn finalize(&self) -> io::Result<()> {
        use std::io::{Seek, SeekFrom, Write as _};
        let mut file = OpenOptions::new().write(true).open(&self.path)?;
        // The frame count is stored at byte 8 (after the 4-byte Fortran record
        // length and the 4-byte "CORD" magic).
        file.seek(SeekFrom::Start(8))?;
        file.write_all(&self.frames_written.to_le_bytes())?;
        file.flush()
    }

    /// Write the DCD file header.
    ///
    /// `n_frames` is the preliminary frame count; it is updated by
    /// [`finalize`](DcdWriter::finalize).
    fn write_header(&mut self, n_frames: u32) -> io::Result<()> {
        use std::io::Write as _;
        let mut f = BufWriter::new(File::create(&self.path)?);

        // ── Block 1 (84 bytes payload) ────────────────────────────────────────
        // Fortran record: 84 bytes.
        let block1_payload: u32 = 84;
        f.write_all(&block1_payload.to_le_bytes())?; // leading sentinel
        f.write_all(b"CORD")?; // magic
        f.write_all(&n_frames.to_le_bytes())?; // nframes  (offset 8)
        f.write_all(&0u32.to_le_bytes())?; // first step
        f.write_all(&1u32.to_le_bytes())?; // step interval
        f.write_all(&0u32.to_le_bytes())?; // last step (unknown)
        f.write_all(&0u32.to_le_bytes())?; // reserved
        f.write_all(&0u32.to_le_bytes())?; // reserved
        f.write_all(&0u32.to_le_bytes())?; // reserved
        f.write_all(&0u32.to_le_bytes())?; // n_free_atoms
        f.write_all(&self.delta.to_le_bytes())?; // timestep (f32)
        // 10 padding integers (unit-cell flag + 9 zeros)
        for _ in 0..10u32 {
            f.write_all(&0u32.to_le_bytes())?;
        }
        // 4-byte CHARMM version word
        f.write_all(&0u32.to_le_bytes())?;
        f.write_all(&block1_payload.to_le_bytes())?; // trailing sentinel

        // ── Block 2: title block ──────────────────────────────────────────────
        // 1 title of 80 bytes.
        let title_payload: u32 = 4 + 80; // n_titles (u32) + 1 × 80-byte title
        f.write_all(&title_payload.to_le_bytes())?;
        f.write_all(&1u32.to_le_bytes())?; // n_titles = 1
        let mut title_buf = [b' '; 80];
        let label = b"Created by oxiphysics DcdWriter";
        let copy_len = label.len().min(80);
        title_buf[..copy_len].copy_from_slice(&label[..copy_len]);
        f.write_all(&title_buf)?;
        f.write_all(&title_payload.to_le_bytes())?;

        // ── Block 3: n_atoms ─────────────────────────────────────────────────
        f.write_all(&4u32.to_le_bytes())?; // payload = 4 bytes
        f.write_all(&self.n_atoms.to_le_bytes())?;
        f.write_all(&4u32.to_le_bytes())?;

        f.flush()
    }
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

    // ── DcdWriter / DcdReader round-trip tests ────────────────────────────────

    #[test]
    fn test_dcd_writer_creates_file() {
        let path = std::env::temp_dir().join("oxiphysics_pf_dcd_create.dcd");
        let _ = std::fs::remove_file(&path);
        let mut w = DcdWriter::new(&path, 3, 2.0e-3_f32).unwrap();
        let frame = make_frame(3);
        w.write_frame(&frame).unwrap();
        w.finalize().unwrap();
        assert!(path.exists());
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_dcd_roundtrip_single_frame() {
        let path = std::env::temp_dir().join("oxiphysics_pf_dcd_rt1.dcd");
        let _ = std::fs::remove_file(&path);
        let frame = make_frame(4);
        let mut w = DcdWriter::new(&path, 4, 1.0e-3_f32).unwrap();
        w.write_frame(&frame).unwrap();
        w.finalize().unwrap();

        let (header, frames) = DcdReader::new(&path).read_all().unwrap();
        assert_eq!(header.n_frames, 1);
        assert_eq!(header.n_atoms, 4);
        assert_eq!(frames.len(), 1);
        for i in 0..4 {
            assert!(
                (frames[0].positions[i][0] - i as f32).abs() < 1e-5,
                "position mismatch at atom {i}"
            );
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_dcd_roundtrip_multi_frame() {
        let path = std::env::temp_dir().join("oxiphysics_pf_dcd_rt3.dcd");
        let _ = std::fs::remove_file(&path);
        let n_atoms = 5u32;
        let mut w = DcdWriter::new(&path, n_atoms, 2.0e-3_f32).unwrap();
        for ts in 0..3u64 {
            let pos: Vec<[f32; 3]> = (0..n_atoms as usize)
                .map(|i| [ts as f32 + i as f32, 0.0, 0.0])
                .collect();
            let types = vec![0u8; n_atoms as usize];
            let frame = ParticleFrame::new(ts, pos, types);
            w.write_frame(&frame).unwrap();
        }
        w.finalize().unwrap();

        let (header, frames) = DcdReader::new(&path).read_all().unwrap();
        assert_eq!(header.n_frames, 3);
        assert_eq!(frames.len(), 3);
        // Verify frame 2, atom 0: x == 2.0
        assert!((frames[2].positions[0][0] - 2.0).abs() < 1e-5);
        let _ = std::fs::remove_file(&path);
    }
}
