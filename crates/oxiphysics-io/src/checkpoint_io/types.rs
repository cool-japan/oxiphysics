//! Auto-generated module
//!
//! 🤖 Generated with [SplitRS](https://github.com/cool-japan/splitrs)

use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

#[allow(unused_imports)]
use super::functions::*;
use super::functions::{
    FORMAT_VERSION, MAGIC, TAG_FOOTER, TAG_INTEGERS, TAG_POSITIONS, TAG_SCALARS, TAG_VELOCITIES,
};

/// A single simulation snapshot with version, timestamp, and raw state bytes.
///
/// `Checkpoint` is a self-contained snapshot that can be serialized to bytes,
/// compared with [`CheckpointDiff`], or stored in a [`CheckpointCatalog`].
#[derive(Debug, Clone, PartialEq)]
pub struct Checkpoint {
    /// Monotonically increasing version counter for this checkpoint sequence.
    pub version: u32,
    /// Unix timestamp (seconds since epoch) when this checkpoint was taken.
    pub timestamp: u64,
    /// Simulation step number.
    pub step: u64,
    /// Simulation time in seconds or simulation units.
    pub sim_time: f64,
    /// Arbitrary raw simulation state bytes.
    pub state: Vec<u8>,
    /// CRC-32 checksum of `state` (set by [`Checkpoint::compute_checksum`]).
    pub checksum: u32,
}
impl Checkpoint {
    /// Construct a new `Checkpoint`, computing the checksum automatically.
    pub fn new(version: u32, timestamp: u64, step: u64, sim_time: f64, state: Vec<u8>) -> Self {
        let checksum = compute_checksum(&state);
        Self {
            version,
            timestamp,
            step,
            sim_time,
            state,
            checksum,
        }
    }
    /// Recompute and store the checksum from the current `state`.
    pub fn compute_checksum(&mut self) {
        self.checksum = compute_checksum(&self.state);
    }
    /// Verify that the stored checksum matches the current `state`.
    pub fn verify(&self) -> bool {
        compute_checksum(&self.state) == self.checksum
    }
    /// Serialize this checkpoint to a flat byte buffer.
    ///
    /// Layout:
    /// - 4 bytes: version (LE u32)
    /// - 8 bytes: timestamp (LE u64)
    /// - 8 bytes: step (LE u64)
    /// - 8 bytes: sim_time (LE f64 bits)
    /// - 4 bytes: checksum (LE u32)
    /// - 8 bytes: state length (LE u64)
    /// - N bytes: state data
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&self.version.to_le_bytes());
        buf.extend_from_slice(&self.timestamp.to_le_bytes());
        buf.extend_from_slice(&self.step.to_le_bytes());
        buf.extend_from_slice(&self.sim_time.to_bits().to_le_bytes());
        buf.extend_from_slice(&self.checksum.to_le_bytes());
        buf.extend_from_slice(&(self.state.len() as u64).to_le_bytes());
        buf.extend_from_slice(&self.state);
        buf
    }
    /// Deserialize a `Checkpoint` from a byte slice.
    pub fn from_bytes(data: &[u8]) -> io::Result<Self> {
        let mut cursor = 0usize;
        let version = read_u32(data, &mut cursor)?;
        let timestamp = read_u64(data, &mut cursor)?;
        let step = read_u64(data, &mut cursor)?;
        let sim_time = read_f64(data, &mut cursor)?;
        let checksum = read_u32(data, &mut cursor)?;
        let state_len = read_u64(data, &mut cursor)? as usize;
        if cursor + state_len > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "state truncated",
            ));
        }
        let state = data[cursor..cursor + state_len].to_vec();
        Ok(Self {
            version,
            timestamp,
            step,
            sim_time,
            state,
            checksum,
        })
    }
}
/// Inspect checkpoint files in a directory without loading full particle data.
#[allow(dead_code)]
#[derive(Debug)]
pub struct CheckpointInspector {
    /// Directory to scan.
    pub base_dir: PathBuf,
}
impl CheckpointInspector {
    /// Create an inspector for `base_dir`.
    pub fn new(base_dir: impl Into<PathBuf>) -> Self {
        Self {
            base_dir: base_dir.into(),
        }
    }
    /// List all `.bin` checkpoint files sorted by name.
    pub fn list(&self) -> Vec<PathBuf> {
        let mut paths: Vec<PathBuf> = match fs::read_dir(&self.base_dir) {
            Ok(rd) => rd
                .filter_map(|e| e.ok())
                .map(|e| e.path())
                .filter(|p| p.extension().is_some_and(|ext| ext == "bin"))
                .collect(),
            Err(_) => Vec::new(),
        };
        paths.sort();
        paths
    }
    /// Read only the header from a checkpoint file (does not load state data).
    pub fn peek_header(&self, path: &Path) -> io::Result<Checkpoint> {
        let data = fs::read(path)?;
        Checkpoint::from_bytes(&data)
    }
    /// Return `(step, sim_time)` for each checkpoint in the directory.
    pub fn metadata_summary(&self) -> Vec<(u64, f64)> {
        self.list()
            .iter()
            .filter_map(|p| self.peek_header(p).ok())
            .map(|c| (c.step, c.sim_time))
            .collect()
    }
    /// Return the number of checkpoints in the directory.
    pub fn count(&self) -> usize {
        self.list().len()
    }
}
/// A complete in-memory simulation snapshot that can be saved to and loaded
/// from a checkpoint file.
#[derive(Debug, Clone)]
pub struct RestartFile {
    /// Checkpoint metadata.
    pub meta: CheckpointMetadata,
    /// Particle positions.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities.
    pub velocities: Vec<[f64; 3]>,
    /// Named scalar data arrays.
    pub scalars: HashMap<String, Vec<f64>>,
}
impl RestartFile {
    /// Construct a new `RestartFile`.
    pub fn new(
        meta: CheckpointMetadata,
        positions: Vec<[f64; 3]>,
        velocities: Vec<[f64; 3]>,
        scalars: HashMap<String, Vec<f64>>,
    ) -> Self {
        Self {
            meta,
            positions,
            velocities,
            scalars,
        }
    }
    /// Write the full restart snapshot to `path`.
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let writer = CheckpointWriter::new(path);
        writer.write_header(&self.meta)?;
        writer.write_positions(&self.positions)?;
        writer.write_velocities(&self.velocities)?;
        let mut keys: Vec<&String> = self.scalars.keys().collect();
        keys.sort();
        for k in keys {
            writer.write_scalars(k, &self.scalars[k])?;
        }
        writer.finalize()
    }
    /// Load a restart snapshot from `path`.
    pub fn load(path: &Path) -> io::Result<Self> {
        let reader = CheckpointReader::new(path);
        let meta = reader.read_metadata()?;
        let positions = reader.read_positions()?;
        let velocities = reader.read_velocities()?;
        let data = fs::read(path)?;
        let mut scalars: HashMap<String, Vec<f64>> = HashMap::new();
        let mut cursor = {
            let mut c = 0usize;
            let _magic = read_u32(&data, &mut c)?;
            let _version = read_u32(&data, &mut c)?;
            let meta_len = read_u32(&data, &mut c)? as usize;
            c += meta_len;
            c
        };
        while cursor < data.len() {
            let tag = data[cursor];
            cursor += 1;
            match tag {
                TAG_POSITIONS | TAG_VELOCITIES => {
                    let count = read_u64(&data, &mut cursor)? as usize;
                    cursor += count * 24;
                }
                TAG_SCALARS => {
                    let name = read_name(&data, &mut cursor)?;
                    let count = read_u64(&data, &mut cursor)? as usize;
                    let mut vals = Vec::with_capacity(count);
                    for _ in 0..count {
                        vals.push(read_f64(&data, &mut cursor)?);
                    }
                    scalars.insert(name, vals);
                }
                TAG_INTEGERS => {
                    let _name = read_name(&data, &mut cursor)?;
                    let count = read_u64(&data, &mut cursor)? as usize;
                    cursor += count * 4;
                }
                TAG_FOOTER => break,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("unknown tag 0x{tag:02X} while loading restart"),
                    ));
                }
            }
        }
        Ok(Self {
            meta,
            positions,
            velocities,
            scalars,
        })
    }
}
/// Supported checkpoint serialisation formats.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckpointFormat {
    /// Raw binary (little-endian f64 arrays).
    Binary,
    /// Human-readable JSON.
    Json,
    /// Run-length encoded binary (LZ4-style RLE).
    Compressed,
    /// HDF5-like hierarchical layout (header + named datasets).
    HDF5Like,
}
impl CheckpointFormat {
    /// File extension for this format.
    pub fn extension(self) -> &'static str {
        match self {
            Self::Binary => "bin",
            Self::Json => "json",
            Self::Compressed => "rle",
            Self::HDF5Like => "h5xt",
        }
    }
    /// Return `true` if the format stores human-readable text.
    pub fn is_text(self) -> bool {
        matches!(self, Self::Json)
    }
}
/// Header record written at the top of every checkpoint file.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct CheckpointHeader {
    /// File format version `[major, minor, patch]`.
    pub version: [u32; 3],
    /// Unix timestamp (seconds) of when the checkpoint was written.
    pub timestamp: u64,
    /// Simulation step number.
    pub step: u64,
    /// Producing crate name.
    pub crate_name: String,
    /// CRC32 checksum of the payload that follows.
    pub checksum: u32,
}
impl CheckpointHeader {
    /// Create a new `CheckpointHeader`.
    pub fn new(
        version: [u32; 3],
        timestamp: u64,
        step: u64,
        crate_name: impl Into<String>,
        checksum: u32,
    ) -> Self {
        Self {
            version,
            timestamp,
            step,
            crate_name: crate_name.into(),
            checksum,
        }
    }
    /// Serialise to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        for v in &self.version {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf.extend_from_slice(&self.timestamp.to_le_bytes());
        buf.extend_from_slice(&self.step.to_le_bytes());
        let name_bytes = self.crate_name.as_bytes();
        buf.extend_from_slice(&(name_bytes.len() as u32).to_le_bytes());
        buf.extend_from_slice(name_bytes);
        buf.extend_from_slice(&self.checksum.to_le_bytes());
        buf
    }
    /// Deserialise from bytes.
    pub fn from_bytes(data: &[u8]) -> io::Result<Self> {
        let mut c = 0usize;
        let v0 = read_u32(data, &mut c)?;
        let v1 = read_u32(data, &mut c)?;
        let v2 = read_u32(data, &mut c)?;
        let timestamp = read_u64(data, &mut c)?;
        let step = read_u64(data, &mut c)?;
        let name_len = read_u32(data, &mut c)? as usize;
        if c + name_len > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "header crate_name truncated",
            ));
        }
        let crate_name = std::str::from_utf8(&data[c..c + name_len])
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?
            .to_owned();
        c += name_len;
        let checksum = read_u32(data, &mut c)?;
        Ok(Self {
            version: [v0, v1, v2],
            timestamp,
            step,
            crate_name,
            checksum,
        })
    }
    /// Return `true` if this header's major version matches `expected`.
    pub fn version_compatible(&self, expected_major: u32) -> bool {
        self.version[0] == expected_major
    }
}
/// Writes a [`Checkpoint`] to a binary file and an accompanying JSON metadata
/// sidecar (`name`.json`).
#[derive(Debug, Clone)]
pub struct CheckpointFileWriter {
    /// Directory in which to write files.
    pub output_dir: PathBuf,
}
impl CheckpointFileWriter {
    /// Create a writer targeting `output_dir` (must exist).
    pub fn new(output_dir: impl Into<PathBuf>) -> Self {
        Self {
            output_dir: output_dir.into(),
        }
    }
    /// Write `checkpoint` to `output_dir`/checkpoint_NNNNNNNNNN.bin` and
    /// emit a JSON metadata sidecar `checkpoint_NNNNNNNNNN.json`.
    pub fn write(&self, checkpoint: &Checkpoint) -> io::Result<PathBuf> {
        let base_name = format!("checkpoint_{:010}", checkpoint.step);
        let bin_path = self.output_dir.join(format!("{base_name}.bin"));
        let json_path = self.output_dir.join(format!("{base_name}.json"));
        fs::write(&bin_path, checkpoint.to_bytes())?;
        let json = format!(
            r#"{{"version":{},"timestamp":{},"step":{},"sim_time":{},"state_len":{},"checksum":{}}}"#,
            checkpoint.version,
            checkpoint.timestamp,
            checkpoint.step,
            checkpoint.sim_time,
            checkpoint.state.len(),
            checkpoint.checksum
        );
        fs::write(&json_path, json.as_bytes())?;
        Ok(bin_path)
    }
}
/// Full simulation snapshot: positions, velocities, forces, and metadata.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct SimulationState {
    /// Particle positions `[x, y, z]`.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities `[vx, vy, vz]`.
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle forces `[fx, fy, fz]`.
    pub forces: Vec<[f64; 3]>,
    /// Arbitrary scalar metadata keyed by name.
    pub metadata: HashMap<String, f64>,
}
impl SimulationState {
    /// Create an empty `SimulationState`.
    pub fn new() -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            metadata: HashMap::new(),
        }
    }
    /// Number of particles in this state.
    pub fn len(&self) -> usize {
        self.positions.len()
    }
    /// Return `true` if there are no particles.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }
    /// Serialise to a flat byte buffer (positions + velocities + forces).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        let n = self.positions.len() as u64;
        buf.extend_from_slice(&n.to_le_bytes());
        for pos in &self.positions {
            for &c in pos {
                buf.extend_from_slice(&c.to_le_bytes());
            }
        }
        for vel in &self.velocities {
            for &c in vel {
                buf.extend_from_slice(&c.to_le_bytes());
            }
        }
        for frc in &self.forces {
            for &c in frc {
                buf.extend_from_slice(&c.to_le_bytes());
            }
        }
        buf
    }
    /// Deserialise from a flat byte buffer.
    pub fn from_bytes(data: &[u8]) -> io::Result<Self> {
        let mut cur = 0usize;
        let n = read_u64(data, &mut cur)? as usize;
        let mut positions = Vec::with_capacity(n);
        let mut velocities = Vec::with_capacity(n);
        let mut forces = Vec::with_capacity(n);
        for _ in 0..n {
            let x = read_f64(data, &mut cur)?;
            let y = read_f64(data, &mut cur)?;
            let z = read_f64(data, &mut cur)?;
            positions.push([x, y, z]);
        }
        for _ in 0..n {
            let x = read_f64(data, &mut cur)?;
            let y = read_f64(data, &mut cur)?;
            let z = read_f64(data, &mut cur)?;
            velocities.push([x, y, z]);
        }
        for _ in 0..n {
            let x = read_f64(data, &mut cur)?;
            let y = read_f64(data, &mut cur)?;
            let z = read_f64(data, &mut cur)?;
            forces.push([x, y, z]);
        }
        Ok(Self {
            positions,
            velocities,
            forces,
            metadata: HashMap::new(),
        })
    }
}
/// Strategy for selecting which checkpoint to restart from.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub enum RestartStrategy {
    /// Restart from the most recent checkpoint.
    FromLatest,
    /// Restart from the checkpoint at a specific step number.
    FromStep(u64),
    /// Restart from a specific file path.
    FromFile(String),
}
impl RestartStrategy {
    /// Resolve this strategy to a file path within `base_dir`.
    pub fn resolve(&self, manager: &CheckpointManager) -> Option<PathBuf> {
        match self {
            Self::FromLatest => manager.latest_checkpoint(),
            Self::FromStep(step) => {
                let p = manager.checkpoint_path(*step);
                if p.exists() { Some(p) } else { None }
            }
            Self::FromFile(path) => {
                let p = PathBuf::from(path);
                if p.exists() { Some(p) } else { None }
            }
        }
    }
    /// Return `true` if this strategy targets the latest checkpoint.
    pub fn is_latest(&self) -> bool {
        matches!(self, Self::FromLatest)
    }
}
/// Writes simulation snapshots to a binary checkpoint file.
///
/// Call methods in order:
/// 1. `write_header`(CheckpointWriter::write_header)
/// 2. [`write_positions`](CheckpointWriter::write_positions) (optional)
/// 3. [`write_velocities`](CheckpointWriter::write_velocities) (optional)
/// 4. Any number of [`write_scalars`](CheckpointWriter::write_scalars) /
///    [`write_integers`](CheckpointWriter::write_integers)
/// 5. [`finalize`](CheckpointWriter::finalize)
#[derive(Debug, Clone)]
pub struct CheckpointWriter {
    /// Destination file path.
    pub path: PathBuf,
    /// Whether to compress the output (currently a no-op placeholder).
    pub compress: bool,
}
impl CheckpointWriter {
    /// Create a new `CheckpointWriter` targeting `path`.
    pub fn new(path: impl Into<PathBuf>) -> Self {
        Self {
            path: path.into(),
            compress: false,
        }
    }
    /// Enable or disable compression (placeholder for future extension).
    pub fn with_compress(mut self, compress: bool) -> Self {
        self.compress = compress;
        self
    }
    /// Write the file header and metadata record.
    ///
    /// Creates (or truncates) the file and writes the magic number, format
    /// version, and serialized metadata.
    pub fn write_header(&self, meta: &CheckpointMetadata) -> io::Result<()> {
        let mut f = BufWriter::new(File::create(&self.path)?);
        f.write_all(&MAGIC.to_le_bytes())?;
        f.write_all(&FORMAT_VERSION.to_le_bytes())?;
        let meta_bytes = meta.to_bytes();
        f.write_all(&(meta_bytes.len() as u32).to_le_bytes())?;
        f.write_all(&meta_bytes)?;
        f.flush()
    }
    /// Append a positions block (`[[f64; 3\]]`) to the checkpoint file.
    ///
    /// Each position is written as three little-endian `f64` values.
    pub fn write_positions(&self, pos: &[[f64; 3]]) -> io::Result<()> {
        self.append_vec3_block(TAG_POSITIONS, pos)
    }
    /// Append a velocities block (`[[f64; 3\]]`) to the checkpoint file.
    pub fn write_velocities(&self, vel: &[[f64; 3]]) -> io::Result<()> {
        self.append_vec3_block(TAG_VELOCITIES, vel)
    }
    /// Append a named scalar (`f64`) array to the checkpoint file.
    pub fn write_scalars(&self, name: &str, data: &[f64]) -> io::Result<()> {
        let mut f = self.open_append()?;
        f.write_all(&[TAG_SCALARS])?;
        write_name(&mut f, name)?;
        f.write_all(&(data.len() as u64).to_le_bytes())?;
        for &v in data {
            f.write_all(&v.to_le_bytes())?;
        }
        f.flush()
    }
    /// Append a named integer (`i32`) array to the checkpoint file.
    pub fn write_integers(&self, name: &str, data: &[i32]) -> io::Result<()> {
        let mut f = self.open_append()?;
        f.write_all(&[TAG_INTEGERS])?;
        write_name(&mut f, name)?;
        f.write_all(&(data.len() as u64).to_le_bytes())?;
        for &v in data {
            f.write_all(&v.to_le_bytes())?;
        }
        f.flush()
    }
    /// Write the footer record containing the file checksum.
    ///
    /// Reads back the entire file written so far, computes the checksum, and
    /// appends a footer tag followed by the 4-byte checksum.
    pub fn finalize(&self) -> io::Result<()> {
        let existing = fs::read(&self.path)?;
        let csum = compute_checksum(&existing);
        let mut f = self.open_append()?;
        f.write_all(&[TAG_FOOTER])?;
        f.write_all(&csum.to_le_bytes())?;
        f.flush()
    }
    fn open_append(&self) -> io::Result<BufWriter<File>> {
        Ok(BufWriter::new(
            OpenOptions::new().append(true).open(&self.path)?,
        ))
    }
    fn append_vec3_block(&self, tag: u8, data: &[[f64; 3]]) -> io::Result<()> {
        let mut f = self.open_append()?;
        f.write_all(&[tag])?;
        f.write_all(&(data.len() as u64).to_le_bytes())?;
        for p in data {
            f.write_all(&p[0].to_le_bytes())?;
            f.write_all(&p[1].to_le_bytes())?;
            f.write_all(&p[2].to_le_bytes())?;
        }
        f.flush()
    }
}
/// Delta checkpoint that records only particles that changed since a base state.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub struct DeltaCheckpoint {
    /// Step number of the base state.
    pub base_step: u64,
    /// Step number of this delta.
    pub target_step: u64,
    /// Indices of changed particles.
    pub changed_indices: Vec<usize>,
    /// New positions for changed particles.
    pub positions: Vec<[f64; 3]>,
    /// New velocities for changed particles.
    pub velocities: Vec<[f64; 3]>,
}
impl DeltaCheckpoint {
    /// Compute a delta between `base` and `target` states.
    ///
    /// A particle is considered changed if any position or velocity component
    /// differs by more than `tol`.
    pub fn compute(
        base_step: u64,
        target_step: u64,
        base: &SimulationState,
        target: &SimulationState,
        tol: f64,
    ) -> Self {
        let n = base.positions.len().min(target.positions.len());
        let mut changed_indices = Vec::new();
        let mut positions = Vec::new();
        let mut velocities = Vec::new();
        for i in 0..n {
            let pos_changed = base.positions[i]
                .iter()
                .zip(target.positions[i].iter())
                .any(|(a, b)| (a - b).abs() > tol);
            let vel_changed = if i < base.velocities.len() && i < target.velocities.len() {
                base.velocities[i]
                    .iter()
                    .zip(target.velocities[i].iter())
                    .any(|(a, b)| (a - b).abs() > tol)
            } else {
                false
            };
            if pos_changed || vel_changed {
                changed_indices.push(i);
                positions.push(target.positions[i]);
                if i < target.velocities.len() {
                    velocities.push(target.velocities[i]);
                } else {
                    velocities.push([0.0; 3]);
                }
            }
        }
        Self {
            base_step,
            target_step,
            changed_indices,
            positions,
            velocities,
        }
    }
    /// Number of changed particles.
    pub fn num_changed(&self) -> usize {
        self.changed_indices.len()
    }
    /// Serialised size in bytes (approximate).
    pub fn byte_size(&self) -> usize {
        16 + self.changed_indices.len() * (8 + 3 * 8 + 3 * 8)
    }
    /// Apply this delta to `base`, returning a new reconstructed state.
    pub fn apply(&self, base: &SimulationState) -> SimulationState {
        let mut out = base.clone();
        for (k, &idx) in self.changed_indices.iter().enumerate() {
            if idx < out.positions.len() {
                out.positions[idx] = self.positions[k];
            }
            if idx < out.velocities.len() && k < self.velocities.len() {
                out.velocities[idx] = self.velocities[k];
            }
        }
        out
    }
}
/// Metadata stored at the head of a checkpoint file.
#[derive(Debug, Clone, PartialEq)]
pub struct CheckpointMetadata {
    /// Simulation step number at which this checkpoint was written.
    pub step: u64,
    /// Simulation time at this checkpoint (in seconds or simulation units).
    pub time: f64,
    /// Number of particles / atoms in the snapshot.
    pub n_particles: usize,
    /// Version of the producing crate `[major, minor, patch]`.
    pub crate_version: [u32; 3],
    /// Human-readable creation timestamp (ISO-8601 string).
    pub created_at: String,
}
impl CheckpointMetadata {
    /// Create a new `CheckpointMetadata`.
    pub fn new(
        step: u64,
        time: f64,
        n_particles: usize,
        crate_version: [u32; 3],
        created_at: impl Into<String>,
    ) -> Self {
        Self {
            step,
            time,
            n_particles,
            crate_version,
            created_at: created_at.into(),
        }
    }
    /// Serialize metadata to a `Vec`u8`.
    pub(crate) fn to_bytes(&self) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&self.step.to_le_bytes());
        buf.extend_from_slice(&self.time.to_le_bytes());
        buf.extend_from_slice(&(self.n_particles as u64).to_le_bytes());
        for v in &self.crate_version {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        let ts = self.created_at.as_bytes();
        buf.extend_from_slice(&(ts.len() as u32).to_le_bytes());
        buf.extend_from_slice(ts);
        buf
    }
    /// Deserialize metadata from a byte slice.
    pub(crate) fn from_bytes(data: &[u8]) -> io::Result<Self> {
        let mut cursor = 0usize;
        let step = read_u64(data, &mut cursor)?;
        let time = read_f64(data, &mut cursor)?;
        let n_particles = read_u64(data, &mut cursor)? as usize;
        let v0 = read_u32(data, &mut cursor)?;
        let v1 = read_u32(data, &mut cursor)?;
        let v2 = read_u32(data, &mut cursor)?;
        let ts_len = read_u32(data, &mut cursor)? as usize;
        if cursor + ts_len > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "created_at string truncated",
            ));
        }
        let created_at = String::from_utf8(data[cursor..cursor + ts_len].to_vec())
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("UTF-8 error: {e}")))?;
        Ok(Self {
            step,
            time,
            n_particles,
            crate_version: [v0, v1, v2],
            created_at,
        })
    }
}
/// Merge partial domain checkpoints from a parallel simulation.
///
/// Each domain produces a `SimulationState` covering a contiguous range of
/// particle indices.  `CheckpointMerger` concatenates the ranges in order.
#[allow(dead_code)]
#[derive(Debug, Default)]
pub struct CheckpointMerger {
    /// Partial states collected so far, each tagged with a domain index.
    pub(super) parts: Vec<(usize, SimulationState)>,
}
impl CheckpointMerger {
    /// Create an empty merger.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a partial state from `domain_id`.
    pub fn add_part(&mut self, domain_id: usize, state: SimulationState) {
        self.parts.push((domain_id, state));
    }
    /// Merge all added parts into one `SimulationState`.
    ///
    /// Parts are sorted by `domain_id` before merging.
    pub fn merge(&mut self) -> SimulationState {
        self.parts.sort_by_key(|(id, _)| *id);
        let mut merged = SimulationState::new();
        for (_, part) in &self.parts {
            merged.positions.extend_from_slice(&part.positions);
            merged.velocities.extend_from_slice(&part.velocities);
            merged.forces.extend_from_slice(&part.forces);
        }
        merged
    }
    /// Number of domain parts added.
    pub fn num_parts(&self) -> usize {
        self.parts.len()
    }
    /// Total particle count across all parts.
    pub fn total_particles(&self) -> usize {
        self.parts.iter().map(|(_, s)| s.len()).sum()
    }
}
/// Binary difference between two checkpoints.
///
/// Stores edit operations as a sequence of byte-level patches that transform
/// `base` into `target`.  Apply with [`CheckpointDiff::apply`].
#[derive(Debug, Clone)]
pub struct CheckpointDiff {
    /// Step number of the base checkpoint.
    pub base_step: u64,
    /// Step number of the target checkpoint.
    pub target_step: u64,
    /// Sequence of `(offset, old_byte, new_byte)` edits.
    pub edits: Vec<(usize, u8, u8)>,
}
impl CheckpointDiff {
    /// Compute the diff between `base` and `target` state byte slices.
    ///
    /// The diff records every position where `base` and `target` differ.
    pub fn compute(base_step: u64, base: &[u8], target_step: u64, target: &[u8]) -> Self {
        let len = base.len().max(target.len());
        let mut edits = Vec::new();
        for i in 0..len {
            let b = if i < base.len() { base[i] } else { 0 };
            let t = if i < target.len() { target[i] } else { 0 };
            if b != t {
                edits.push((i, b, t));
            }
        }
        Self {
            base_step,
            target_step,
            edits,
        }
    }
    /// Apply the diff to `base_state` to reconstruct `target_state`.
    ///
    /// Returns an error if any edit offset falls out of bounds.
    pub fn apply(&self, base_state: &[u8]) -> io::Result<Vec<u8>> {
        let mut out = base_state.to_vec();
        let max_off = self.edits.iter().map(|&(o, _, _)| o).max().unwrap_or(0);
        if max_off >= out.len() && !self.edits.is_empty() {
            out.resize(max_off + 1, 0);
        }
        for &(offset, _old, new) in &self.edits {
            if offset >= out.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("diff edit offset {offset} out of bounds"),
                ));
            }
            out[offset] = new;
        }
        Ok(out)
    }
    /// Number of bytes that differ between base and target.
    pub fn diff_size(&self) -> usize {
        self.edits.len()
    }
    /// Fraction of bytes that changed (relative to the base state length).
    pub fn change_ratio(&self, base_len: usize) -> f64 {
        if base_len == 0 {
            return 0.0_f64;
        }
        self.edits.len() as f64 / base_len as f64
    }
}
/// Manages a rolling set of checkpoint files in a directory.
///
/// Checkpoints are named `checkpoint_NNNNNNNNNN.bin` where `N` is the
/// zero-padded step number.
#[derive(Debug, Clone)]
pub struct CheckpointManager {
    /// Base directory for checkpoint files.
    pub base_dir: PathBuf,
    /// Maximum number of checkpoints to keep.
    pub max_checkpoints: usize,
    /// Write a checkpoint every `interval_steps` steps.
    pub interval_steps: u64,
}
impl CheckpointManager {
    /// Construct a new `CheckpointManager`.
    pub fn new(base_dir: impl Into<PathBuf>, max_checkpoints: usize, interval_steps: u64) -> Self {
        Self {
            base_dir: base_dir.into(),
            max_checkpoints,
            interval_steps,
        }
    }
    /// Return `true` when a checkpoint should be written at `step`.
    ///
    /// Always returns `true` for step 0 if `interval_steps > 0`.
    pub fn should_checkpoint(&self, step: u64) -> bool {
        if self.interval_steps == 0 {
            return false;
        }
        step.is_multiple_of(self.interval_steps)
    }
    /// Compute the canonical path for the checkpoint at `step`.
    pub fn checkpoint_path(&self, step: u64) -> PathBuf {
        self.base_dir.join(format!("checkpoint_{step:010}.bin"))
    }
    /// Return a sorted list of existing checkpoint file paths in `base_dir`.
    pub fn list_checkpoints(&self) -> Vec<PathBuf> {
        let Ok(entries) = fs::read_dir(&self.base_dir) else {
            return vec![];
        };
        let mut paths: Vec<PathBuf> = entries
            .flatten()
            .filter_map(|e| {
                let p = e.path();
                let name = p.file_name()?.to_string_lossy().into_owned();
                if name.starts_with("checkpoint_") && name.ends_with(".bin") {
                    Some(p)
                } else {
                    None
                }
            })
            .collect();
        paths.sort();
        paths
    }
    /// Return the path to the most recent (highest step) checkpoint, if any.
    pub fn latest_checkpoint(&self) -> Option<PathBuf> {
        self.list_checkpoints().into_iter().last()
    }
    /// Delete oldest checkpoints, keeping at most `max_checkpoints` files.
    pub fn prune_old_checkpoints(&self) -> io::Result<()> {
        let checkpoints = self.list_checkpoints();
        if checkpoints.len() <= self.max_checkpoints {
            return Ok(());
        }
        let to_delete = checkpoints.len() - self.max_checkpoints;
        for path in checkpoints.iter().take(to_delete) {
            fs::remove_file(path)?;
        }
        Ok(())
    }
}
/// LZ77-style run-length + match compressor for checkpoint state bytes.
///
/// The algorithm alternates between:
/// - Literal runs: a length byte followed by raw bytes.
/// - Back-references: a match entry `(offset, length)` pointing back into
///   the already-decoded output.
///
/// **Wire format** (simplified):
/// - `0x00 len data…` — literal run of `len` bytes.
/// - `0x01 offset_lo offset_hi len` — back-reference, copy `len` bytes from
///   `output\[pos - offset\]`.
#[derive(Debug, Clone, Default)]
pub struct CheckpointCompressor {
    /// Minimum match length to emit a back-reference (default: 3).
    pub min_match_len: usize,
    /// Maximum look-back distance for finding matches (default: 255).
    pub max_look_back: usize,
}
impl CheckpointCompressor {
    /// Create a new compressor with default settings.
    pub fn new() -> Self {
        Self {
            min_match_len: 3,
            max_look_back: 255,
        }
    }
    /// Compress `input` and return the compressed bytes.
    pub fn compress(&self, input: &[u8]) -> Vec<u8> {
        let min_match = self.min_match_len.max(1);
        let look_back = self.max_look_back.max(1);
        let mut out = Vec::new();
        let mut pos = 0usize;
        while pos < input.len() {
            let window_start = pos.saturating_sub(look_back);
            let mut best_off = 0usize;
            let mut best_len = 0usize;
            for start in window_start..pos {
                let mut len = 0usize;
                while pos + len < input.len() && input[start + len] == input[pos + len] && len < 255
                {
                    len += 1;
                    if start + len >= pos {
                        break;
                    }
                }
                if len > best_len && len >= min_match {
                    best_len = len;
                    best_off = pos - start;
                }
            }
            if best_len >= min_match {
                out.push(0x01);
                out.push((best_off & 0xFF) as u8);
                out.push(((best_off >> 8) & 0xFF) as u8);
                out.push(best_len as u8);
                pos += best_len;
            } else {
                let run_end = (pos + 255).min(input.len());
                let run_len = run_end - pos;
                out.push(0x00);
                out.push(run_len as u8);
                out.extend_from_slice(&input[pos..pos + run_len]);
                pos += run_len;
            }
        }
        out
    }
    /// Decompress bytes produced by [`compress`](CheckpointCompressor::compress).
    pub fn decompress(&self, input: &[u8]) -> io::Result<Vec<u8>> {
        let mut out: Vec<u8> = Vec::new();
        let mut i = 0usize;
        while i < input.len() {
            let tag = input[i];
            i += 1;
            match tag {
                0x00 => {
                    if i >= input.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::UnexpectedEof,
                            "literal run truncated",
                        ));
                    }
                    let run_len = input[i] as usize;
                    i += 1;
                    if i + run_len > input.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::UnexpectedEof,
                            "literal data truncated",
                        ));
                    }
                    out.extend_from_slice(&input[i..i + run_len]);
                    i += run_len;
                }
                0x01 => {
                    if i + 3 > input.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::UnexpectedEof,
                            "back-ref truncated",
                        ));
                    }
                    let off_lo = input[i] as usize;
                    let off_hi = input[i + 1] as usize;
                    let offset = off_lo | (off_hi << 8);
                    let length = input[i + 2] as usize;
                    i += 3;
                    if offset == 0 || offset > out.len() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("invalid back-ref offset {offset}"),
                        ));
                    }
                    let start = out.len() - offset;
                    for k in 0..length {
                        let byte = out[start + k];
                        out.push(byte);
                    }
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("unknown tag 0x{tag:02X}"),
                    ));
                }
            }
        }
        Ok(out)
    }
    /// Return the compression ratio (`compressed.len() / original.len()`).
    ///
    /// Returns `1.0_f64` for empty input (no compression).
    pub fn compression_ratio(original_len: usize, compressed_len: usize) -> f64 {
        if original_len == 0 {
            return 1.0_f64;
        }
        compressed_len as f64 / original_len as f64
    }
}
/// A catalog of available checkpoints with lazy loading from a base directory.
///
/// The catalog scans `base_dir` on construction and builds an index of
/// `(step, path)` pairs.  Actual checkpoint bytes are loaded on demand.
#[derive(Debug, Clone)]
pub struct CheckpointCatalog {
    /// Base directory containing checkpoint files.
    pub base_dir: PathBuf,
    /// Sorted list of `(step, path)` entries.
    pub entries: Vec<(u64, PathBuf)>,
}
impl CheckpointCatalog {
    /// Scan `base_dir` for checkpoint files and build the catalog.
    ///
    /// Files must be named `checkpoint_NNNNNNNNNN.bin` (10-digit zero-padded
    /// step number).
    pub fn scan(base_dir: impl Into<PathBuf>) -> Self {
        let base_dir: PathBuf = base_dir.into();
        let mut entries: Vec<(u64, PathBuf)> = Vec::new();
        if let Ok(dir_entries) = fs::read_dir(&base_dir) {
            for entry in dir_entries.flatten() {
                let path = entry.path();
                if let Some(name) = path.file_name().and_then(|n| n.to_str())
                    && name.starts_with("checkpoint_")
                    && name.ends_with(".bin")
                {
                    let step_str = &name[11..name.len() - 4];
                    if let Ok(step) = step_str.parse::<u64>() {
                        entries.push((step, path));
                    }
                }
            }
        }
        entries.sort_by_key(|(s, _)| *s);
        Self { base_dir, entries }
    }
    /// Return the number of checkpoints in the catalog.
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    /// Return `true` if the catalog is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
    /// Return the step numbers of all indexed checkpoints.
    pub fn steps(&self) -> Vec<u64> {
        self.entries.iter().map(|(s, _)| *s).collect()
    }
    /// Find the path for `step`, if it exists in the catalog.
    pub fn path_for_step(&self, step: u64) -> Option<&PathBuf> {
        self.entries
            .binary_search_by_key(&step, |(s, _)| *s)
            .ok()
            .map(|idx| &self.entries[idx].1)
    }
    /// Load and return the `Checkpoint` for the given `step`.
    ///
    /// Returns an error if the step is not catalogued or the file cannot be
    /// read / parsed.
    pub fn load_step(&self, step: u64) -> io::Result<Checkpoint> {
        let path = self.path_for_step(step).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::NotFound,
                format!("step {step} not in catalog"),
            )
        })?;
        let data = fs::read(path)?;
        Checkpoint::from_bytes(&data)
    }
    /// Return the latest (highest step) checkpoint path, if any.
    pub fn latest(&self) -> Option<&PathBuf> {
        self.entries.last().map(|(_, p)| p)
    }
    /// Return the earliest (lowest step) checkpoint path, if any.
    pub fn earliest(&self) -> Option<&PathBuf> {
        self.entries.first().map(|(_, p)| p)
    }
    /// Add a `Checkpoint` to the catalog by saving it to disk and indexing it.
    ///
    /// Uses `CheckpointManager`-compatible naming: `checkpoint_NNNNNNNNNN.bin`.
    pub fn add(&mut self, checkpoint: &Checkpoint) -> io::Result<()> {
        let path = self
            .base_dir
            .join(format!("checkpoint_{:010}.bin", checkpoint.step));
        let bytes = checkpoint.to_bytes();
        fs::write(&path, &bytes)?;
        let pos = self.entries.partition_point(|(s, _)| *s < checkpoint.step);
        self.entries.insert(pos, (checkpoint.step, path));
        Ok(())
    }
    /// Remove the catalog entry for `step` and delete the file from disk.
    pub fn remove_step(&mut self, step: u64) -> io::Result<()> {
        let pos = self
            .entries
            .binary_search_by_key(&step, |(s, _)| *s)
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::NotFound,
                    format!("step {step} not in catalog"),
                )
            })?;
        let (_, path) = self.entries.remove(pos);
        if path.exists() {
            fs::remove_file(&path)?;
        }
        Ok(())
    }
}
/// Reads simulation snapshots from a binary checkpoint file.
#[derive(Debug, Clone)]
pub struct CheckpointReader {
    /// Source file path.
    pub path: PathBuf,
}
impl CheckpointReader {
    /// Create a new `CheckpointReader` from a path.
    pub fn new(path: impl Into<PathBuf>) -> Self {
        Self { path: path.into() }
    }
    /// Read and return the checkpoint metadata from the file header.
    pub fn read_metadata(&self) -> io::Result<CheckpointMetadata> {
        let data = fs::read(&self.path)?;
        let mut cursor = 0usize;
        let magic = read_u32(&data, &mut cursor)?;
        if magic != MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "bad magic number",
            ));
        }
        let _version = read_u32(&data, &mut cursor)?;
        let meta_len = read_u32(&data, &mut cursor)? as usize;
        if cursor + meta_len > data.len() {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "metadata block truncated",
            ));
        }
        CheckpointMetadata::from_bytes(&data[cursor..cursor + meta_len])
    }
    /// Read and return all particle positions stored in the file.
    pub fn read_positions(&self) -> io::Result<Vec<[f64; 3]>> {
        self.read_vec3_block(TAG_POSITIONS)
    }
    /// Read and return all particle velocities stored in the file.
    pub fn read_velocities(&self) -> io::Result<Vec<[f64; 3]>> {
        self.read_vec3_block(TAG_VELOCITIES)
    }
    /// Read and return the named scalar array from the file.
    pub fn read_scalars(&self, name: &str) -> io::Result<Vec<f64>> {
        let data = fs::read(&self.path)?;
        let mut cursor = self.skip_header(&data)?;
        while cursor < data.len() {
            let tag = data[cursor];
            cursor += 1;
            match tag {
                TAG_SCALARS => {
                    let stored_name = read_name(&data, &mut cursor)?;
                    let count = read_u64(&data, &mut cursor)? as usize;
                    if stored_name == name {
                        let mut out = Vec::with_capacity(count);
                        for _ in 0..count {
                            out.push(read_f64(&data, &mut cursor)?);
                        }
                        return Ok(out);
                    } else {
                        cursor += count * 8;
                    }
                }
                TAG_POSITIONS | TAG_VELOCITIES => {
                    let count = read_u64(&data, &mut cursor)? as usize;
                    cursor += count * 24;
                }
                TAG_INTEGERS => {
                    let _n = read_name(&data, &mut cursor)?;
                    let count = read_u64(&data, &mut cursor)? as usize;
                    cursor += count * 4;
                }
                TAG_FOOTER => break,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("unknown tag 0x{tag:02X}"),
                    ));
                }
            }
        }
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            format!("scalar array '{name}' not found"),
        ))
    }
    fn skip_header(&self, data: &[u8]) -> io::Result<usize> {
        let mut cursor = 0usize;
        let _magic = read_u32(data, &mut cursor)?;
        let _version = read_u32(data, &mut cursor)?;
        let meta_len = read_u32(data, &mut cursor)? as usize;
        cursor += meta_len;
        Ok(cursor)
    }
    fn read_vec3_block(&self, target_tag: u8) -> io::Result<Vec<[f64; 3]>> {
        let data = fs::read(&self.path)?;
        let mut cursor = self.skip_header(&data)?;
        while cursor < data.len() {
            let tag = data[cursor];
            cursor += 1;
            match tag {
                t if t == target_tag => {
                    let count = read_u64(&data, &mut cursor)? as usize;
                    let mut out = Vec::with_capacity(count);
                    for _ in 0..count {
                        let x = read_f64(&data, &mut cursor)?;
                        let y = read_f64(&data, &mut cursor)?;
                        let z = read_f64(&data, &mut cursor)?;
                        out.push([x, y, z]);
                    }
                    return Ok(out);
                }
                TAG_POSITIONS | TAG_VELOCITIES => {
                    let count = read_u64(&data, &mut cursor)? as usize;
                    cursor += count * 24;
                }
                TAG_SCALARS | TAG_INTEGERS => {
                    let _n = read_name(&data, &mut cursor)?;
                    let count = read_u64(&data, &mut cursor)? as usize;
                    let elem_size = if tag == TAG_SCALARS { 8 } else { 4 };
                    cursor += count * elem_size;
                }
                TAG_FOOTER => break,
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("unknown tag 0x{tag:02X}"),
                    ));
                }
            }
        }
        Ok(vec![])
    }
}
/// Reads and validates a [`Checkpoint`] from a binary file.
#[derive(Debug, Clone)]
pub struct CheckpointFileReader {
    /// Path to the binary checkpoint file.
    pub path: PathBuf,
}
impl CheckpointFileReader {
    /// Create a reader for `path`.
    pub fn new(path: impl Into<PathBuf>) -> Self {
        Self { path: path.into() }
    }
    /// Read, parse, and validate the checkpoint.
    ///
    /// Returns an error if the file cannot be read, the byte format is
    /// invalid, or the embedded checksum does not match the state data.
    pub fn read_and_validate(&self) -> io::Result<Checkpoint> {
        let data = fs::read(&self.path)?;
        let ckpt = Checkpoint::from_bytes(&data)?;
        if !ckpt.verify() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "checkpoint checksum mismatch",
            ));
        }
        Ok(ckpt)
    }
}
