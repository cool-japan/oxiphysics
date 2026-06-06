// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation restart and checkpoint I/O.
//!
//! Provides writers and readers for simulation restart files in multiple
//! formats (binary, ASCII, JSON, simplified HDF5-like, MessagePack-like),
//! plus a rolling [`CheckpointManager`] and a [`RestartValidator`] for
//! checksum verification.

// ── Magic bytes ───────────────────────────────────────────────────────────────

/// Magic bytes written at the start of every binary restart file: "OXRS".
const BINARY_MAGIC: [u8; 4] = [0x4F, 0x58, 0x52, 0x53];

/// Format version embedded in binary headers.
const FORMAT_VERSION: u32 = 1;

// ── RestartFormat ─────────────────────────────────────────────────────────────

/// Serialization format for restart / checkpoint files.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RestartFormat {
    /// Custom little-endian binary format (fastest, smallest).
    Binary,
    /// Human-readable ASCII key-value format.
    Ascii,
    /// Simplified HDF5-like layout (group/dataset tags in binary).
    Hdf5Like,
    /// JSON text format (portable, verbose).
    Json,
    /// Simplified MessagePack-like binary encoding.
    MessagePack,
}

// ── RestartMetadata ───────────────────────────────────────────────────────────

/// Metadata stored alongside a restart snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct RestartMetadata {
    /// File-format version string (e.g. `"1.0"`).
    pub version: String,
    /// UNIX timestamp (seconds since epoch) when this file was written.
    pub timestamp: u64,
    /// Simulation step number.
    pub step: u64,
    /// Simulation time (in simulation units).
    pub time: f64,
    /// Name of the crate that produced this file.
    pub crate_name: String,
    /// Free-form human-readable description.
    pub description: String,
}

impl RestartMetadata {
    /// Construct a new `RestartMetadata`.
    pub fn new(
        version: impl Into<String>,
        timestamp: u64,
        step: u64,
        time: f64,
        crate_name: impl Into<String>,
        description: impl Into<String>,
    ) -> Self {
        Self {
            version: version.into(),
            timestamp,
            step,
            time,
            crate_name: crate_name.into(),
            description: description.into(),
        }
    }

    /// Return a default metadata for testing.
    pub fn default_test() -> Self {
        Self::new("1.0", 0, 0, 0.0, "oxiphysics", "test checkpoint")
    }
}

// ── RestartData ───────────────────────────────────────────────────────────────

/// Full simulation state stored in a restart file.
#[derive(Debug, Clone, PartialEq)]
pub struct RestartData {
    /// File metadata.
    pub metadata: RestartMetadata,
    /// Particle positions \[x, y, z\] (simulation units).
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities \[vx, vy, vz\].
    pub velocities: Vec<[f64; 3]>,
    /// Per-particle forces \[fx, fy, fz\].
    pub forces: Vec<[f64; 3]>,
    /// Per-particle masses.
    pub masses: Vec<f64>,
    /// Per-particle type indices.
    pub types: Vec<u32>,
    /// Simulation box matrix (row-major: rows are box vectors a, b, c).
    pub box_matrix: [[f64; 3]; 3],
    /// Named extra scalar arrays `(name, values)`.
    pub extra_scalars: Vec<(String, Vec<f64>)>,
    /// Named extra vector arrays `(name, vectors)`.
    pub extra_vectors: Vec<(String, Vec<[f64; 3]>)>,
}

impl RestartData {
    /// Number of particles in this snapshot.
    pub fn n_particles(&self) -> usize {
        self.positions.len()
    }

    /// Construct a minimal empty `RestartData`.
    pub fn empty(metadata: RestartMetadata) -> Self {
        Self {
            metadata,
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            types: Vec::new(),
            box_matrix: [[0.0; 3]; 3],
            extra_scalars: Vec::new(),
            extra_vectors: Vec::new(),
        }
    }

    /// Build a simple single-particle `RestartData` for tests.
    pub fn single_particle_test() -> Self {
        let meta = RestartMetadata::default_test();
        let mut d = Self::empty(meta);
        d.positions = vec![[1.0, 2.0, 3.0]];
        d.velocities = vec![[0.1, 0.2, 0.3]];
        d.forces = vec![[0.0, -9.81, 0.0]];
        d.masses = vec![1.0];
        d.types = vec![0];
        d.box_matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        d
    }
}

// ── Helper: encode / decode primitives (little-endian) ────────────────────────

fn encode_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn encode_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn encode_f64(buf: &mut Vec<u8>, v: f64) {
    buf.extend_from_slice(&v.to_bits().to_le_bytes());
}

fn encode_str(buf: &mut Vec<u8>, s: &str) {
    let bytes = s.as_bytes();
    encode_u32(buf, bytes.len() as u32);
    buf.extend_from_slice(bytes);
}

fn encode_vec3(buf: &mut Vec<u8>, v: &[f64; 3]) {
    encode_f64(buf, v[0]);
    encode_f64(buf, v[1]);
    encode_f64(buf, v[2]);
}

fn decode_u32(data: &[u8], offset: &mut usize) -> Result<u32, String> {
    if *offset + 4 > data.len() {
        return Err(format!("unexpected EOF at offset {}", *offset));
    }
    let v = u32::from_le_bytes(
        data[*offset..*offset + 4]
            .try_into()
            .expect("slice length must match"),
    );
    *offset += 4;
    Ok(v)
}

fn decode_u64(data: &[u8], offset: &mut usize) -> Result<u64, String> {
    if *offset + 8 > data.len() {
        return Err(format!("unexpected EOF at offset {}", *offset));
    }
    let v = u64::from_le_bytes(
        data[*offset..*offset + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *offset += 8;
    Ok(v)
}

fn decode_f64(data: &[u8], offset: &mut usize) -> Result<f64, String> {
    if *offset + 8 > data.len() {
        return Err(format!("unexpected EOF at offset {}", *offset));
    }
    let bits = u64::from_le_bytes(
        data[*offset..*offset + 8]
            .try_into()
            .expect("slice length must match"),
    );
    *offset += 8;
    Ok(f64::from_bits(bits))
}

fn decode_str(data: &[u8], offset: &mut usize) -> Result<String, String> {
    let len = decode_u32(data, offset)? as usize;
    if *offset + len > data.len() {
        return Err(format!("string extends past EOF at offset {}", *offset));
    }
    let s = std::str::from_utf8(&data[*offset..*offset + len])
        .map_err(|e| format!("UTF-8 error: {e}"))?
        .to_string();
    *offset += len;
    Ok(s)
}

fn decode_vec3(data: &[u8], offset: &mut usize) -> Result<[f64; 3], String> {
    let x = decode_f64(data, offset)?;
    let y = decode_f64(data, offset)?;
    let z = decode_f64(data, offset)?;
    Ok([x, y, z])
}

// ── RestartWriter ─────────────────────────────────────────────────────────────

/// Writes simulation restart data to files or in-memory buffers.
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

    /// Serialize `data` according to `self.format` and write to `self.path`.
    ///
    /// Returns an error string on I/O or encoding failure.
    pub fn write(&self, data: &RestartData) -> Result<(), String> {
        let bytes: Vec<u8> = match &self.format {
            RestartFormat::Binary => Self::write_binary(data),
            RestartFormat::Ascii => Self::write_ascii(data).into_bytes(),
            RestartFormat::Json => Self::write_json(data).into_bytes(),
            RestartFormat::Hdf5Like => Self::write_hdf5like(data),
            RestartFormat::MessagePack => Self::write_msgpack(data),
        };
        std::fs::write(&self.path, &bytes)
            .map_err(|e| format!("failed to write '{}': {e}", self.path))
    }

    /// Encode `data` as a custom little-endian binary blob.
    ///
    /// Layout: `[MAGIC 4B][VERSION 4B][METADATA][N_PARTICLES 8B]
    ///          [POSITIONS…][VELOCITIES…][FORCES…][MASSES…][TYPES…]
    ///          [BOX_MATRIX 72B][N_EXTRA_SCALARS 4B][…][N_EXTRA_VECS 4B][…]`
    pub fn write_binary(data: &RestartData) -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend_from_slice(&BINARY_MAGIC);
        encode_u32(&mut buf, FORMAT_VERSION);
        // metadata
        encode_str(&mut buf, &data.metadata.version);
        encode_u64(&mut buf, data.metadata.timestamp);
        encode_u64(&mut buf, data.metadata.step);
        encode_f64(&mut buf, data.metadata.time);
        encode_str(&mut buf, &data.metadata.crate_name);
        encode_str(&mut buf, &data.metadata.description);
        // particle count
        let n = data.n_particles() as u64;
        encode_u64(&mut buf, n);
        for p in &data.positions {
            encode_vec3(&mut buf, p);
        }
        for v in &data.velocities {
            encode_vec3(&mut buf, v);
        }
        for f in &data.forces {
            encode_vec3(&mut buf, f);
        }
        for m in &data.masses {
            encode_f64(&mut buf, *m);
        }
        for t in &data.types {
            encode_u32(&mut buf, *t);
        }
        // box matrix (9 f64)
        for row in &data.box_matrix {
            for &c in row {
                encode_f64(&mut buf, c);
            }
        }
        // extra scalars
        encode_u32(&mut buf, data.extra_scalars.len() as u32);
        for (name, vals) in &data.extra_scalars {
            encode_str(&mut buf, name);
            encode_u64(&mut buf, vals.len() as u64);
            for &v in vals {
                encode_f64(&mut buf, v);
            }
        }
        // extra vectors
        encode_u32(&mut buf, data.extra_vectors.len() as u32);
        for (name, vecs) in &data.extra_vectors {
            encode_str(&mut buf, name);
            encode_u64(&mut buf, vecs.len() as u64);
            for v in vecs {
                encode_vec3(&mut buf, v);
            }
        }
        buf
    }

    /// Encode `data` as a human-readable ASCII text block.
    pub fn write_ascii(data: &RestartData) -> String {
        let mut s = String::new();
        s.push_str("# OxiPhysics restart file\n");
        s.push_str(&format!("VERSION {}\n", data.metadata.version));
        s.push_str(&format!("TIMESTAMP {}\n", data.metadata.timestamp));
        s.push_str(&format!("STEP {}\n", data.metadata.step));
        s.push_str(&format!("TIME {:.6}\n", data.metadata.time));
        s.push_str(&format!("CRATE {}\n", data.metadata.crate_name));
        s.push_str(&format!("DESC {}\n", data.metadata.description));
        let n = data.n_particles();
        s.push_str(&format!("N_PARTICLES {n}\n"));
        s.push_str("BEGIN_POSITIONS\n");
        for p in &data.positions {
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", p[0], p[1], p[2]));
        }
        s.push_str("END_POSITIONS\n");
        s.push_str("BEGIN_VELOCITIES\n");
        for v in &data.velocities {
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", v[0], v[1], v[2]));
        }
        s.push_str("END_VELOCITIES\n");
        s.push_str("BEGIN_FORCES\n");
        for f in &data.forces {
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", f[0], f[1], f[2]));
        }
        s.push_str("END_FORCES\n");
        s.push_str("BEGIN_MASSES\n");
        for m in &data.masses {
            s.push_str(&format!("{:.6}\n", m));
        }
        s.push_str("END_MASSES\n");
        s.push_str("BEGIN_TYPES\n");
        for t in &data.types {
            s.push_str(&format!("{t}\n"));
        }
        s.push_str("END_TYPES\n");
        s.push_str("BEGIN_BOX\n");
        for row in &data.box_matrix {
            s.push_str(&format!("{:.6} {:.6} {:.6}\n", row[0], row[1], row[2]));
        }
        s.push_str("END_BOX\n");
        s
    }

    /// Encode `data` as a JSON string (hand-rolled, no external dependency).
    pub fn write_json(data: &RestartData) -> String {
        let mut s = String::new();
        s.push_str("{\n");
        s.push_str(&format!(
            "  \"version\": \"{}\",\n  \"timestamp\": {},\n  \"step\": {},\n  \"time\": {:.6},\n  \"crate\": \"{}\",\n  \"description\": \"{}\",\n",
            data.metadata.version, data.metadata.timestamp, data.metadata.step,
            data.metadata.time, data.metadata.crate_name, data.metadata.description
        ));
        s.push_str(&format!("  \"n_particles\": {},\n", data.n_particles()));
        // positions
        s.push_str("  \"positions\": [");
        for (i, p) in data.positions.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("[{:.6},{:.6},{:.6}]", p[0], p[1], p[2]));
        }
        s.push_str("],\n");
        // velocities
        s.push_str("  \"velocities\": [");
        for (i, v) in data.velocities.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("[{:.6},{:.6},{:.6}]", v[0], v[1], v[2]));
        }
        s.push_str("],\n");
        // masses
        s.push_str("  \"masses\": [");
        for (i, m) in data.masses.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("{:.6}", m));
        }
        s.push_str("],\n");
        // types
        s.push_str("  \"types\": [");
        for (i, t) in data.types.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&format!("{t}"));
        }
        s.push_str("]\n}\n");
        s
    }

    /// Encode `data` using a simplified HDF5-like tagged binary layout.
    pub fn write_hdf5like(data: &RestartData) -> Vec<u8> {
        // Reuse the binary writer with an alternate magic to mark it as hdf5-like.
        let mut buf = Self::write_binary(data);
        // Overwrite magic bytes to distinguish from plain Binary.
        buf[0] = 0x4F; // 'O'
        buf[1] = 0x58; // 'X'
        buf[2] = 0x48; // 'H'
        buf[3] = 0x35; // '5'
        buf
    }

    /// Encode `data` using a simplified MessagePack-like binary encoding.
    pub fn write_msgpack(data: &RestartData) -> Vec<u8> {
        // Reuse the binary writer with an alternate magic.
        let mut buf = Self::write_binary(data);
        buf[0] = 0x4D; // 'M'
        buf[1] = 0x50; // 'P'
        buf[2] = 0x4B; // 'K'
        buf[3] = 0x31; // '1'
        buf
    }
}

// ── RestartReader ─────────────────────────────────────────────────────────────

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
    pub fn read(&self) -> Result<RestartData, String> {
        let bytes = std::fs::read(&self.path)
            .map_err(|e| format!("failed to read '{}': {e}", self.path))?;
        let fmt = Self::detect_format(&bytes);
        match fmt {
            RestartFormat::Ascii => {
                let text = std::str::from_utf8(&bytes).map_err(|e| format!("UTF-8 error: {e}"))?;
                Self::read_ascii(text)
            }
            _ => Self::read_binary(&bytes),
        }
    }

    /// Detect the serialization format by inspecting magic bytes.
    pub fn detect_format(bytes: &[u8]) -> RestartFormat {
        if bytes.len() < 4 {
            return RestartFormat::Ascii;
        }
        match &bytes[0..4] {
            b"OXRS" => RestartFormat::Binary,
            b"OXHS" | [0x4F, 0x58, 0x48, 0x35] => RestartFormat::Hdf5Like,
            [0x4D, 0x50, 0x4B, 0x31] => RestartFormat::MessagePack,
            b"# Ox" | b"VERS" => RestartFormat::Ascii,
            _ if bytes.starts_with(b"{") => RestartFormat::Json,
            _ => RestartFormat::Ascii,
        }
    }

    /// Parse a `RestartData` from a raw binary blob (any of the binary formats).
    pub fn read_binary(bytes: &[u8]) -> Result<RestartData, String> {
        let mut off = 0usize;
        // magic (4) + version (4)
        if bytes.len() < 8 {
            return Err("binary too short".into());
        }
        off += 4; // skip magic
        let _file_version = decode_u32(bytes, &mut off)?;
        // metadata
        let version = decode_str(bytes, &mut off)?;
        let timestamp = decode_u64(bytes, &mut off)?;
        let step = decode_u64(bytes, &mut off)?;
        let time = decode_f64(bytes, &mut off)?;
        let crate_name = decode_str(bytes, &mut off)?;
        let description = decode_str(bytes, &mut off)?;
        let metadata = RestartMetadata {
            version,
            timestamp,
            step,
            time,
            crate_name,
            description,
        };
        // particle count
        let n = decode_u64(bytes, &mut off)? as usize;
        let mut positions = Vec::with_capacity(n);
        for _ in 0..n {
            positions.push(decode_vec3(bytes, &mut off)?);
        }
        let mut velocities = Vec::with_capacity(n);
        for _ in 0..n {
            velocities.push(decode_vec3(bytes, &mut off)?);
        }
        let mut forces = Vec::with_capacity(n);
        for _ in 0..n {
            forces.push(decode_vec3(bytes, &mut off)?);
        }
        let mut masses = Vec::with_capacity(n);
        for _ in 0..n {
            masses.push(decode_f64(bytes, &mut off)?);
        }
        let mut types = Vec::with_capacity(n);
        for _ in 0..n {
            types.push(decode_u32(bytes, &mut off)?);
        }
        // box matrix
        let mut box_matrix = [[0.0f64; 3]; 3];
        for row in &mut box_matrix {
            for c in row.iter_mut() {
                *c = decode_f64(bytes, &mut off)?;
            }
        }
        // extra scalars
        let n_es = decode_u32(bytes, &mut off)? as usize;
        let mut extra_scalars = Vec::with_capacity(n_es);
        for _ in 0..n_es {
            let name = decode_str(bytes, &mut off)?;
            let count = decode_u64(bytes, &mut off)? as usize;
            let mut vals = Vec::with_capacity(count);
            for _ in 0..count {
                vals.push(decode_f64(bytes, &mut off)?);
            }
            extra_scalars.push((name, vals));
        }
        // extra vectors
        let n_ev = decode_u32(bytes, &mut off)? as usize;
        let mut extra_vectors = Vec::with_capacity(n_ev);
        for _ in 0..n_ev {
            let name = decode_str(bytes, &mut off)?;
            let count = decode_u64(bytes, &mut off)? as usize;
            let mut vecs = Vec::with_capacity(count);
            for _ in 0..count {
                vecs.push(decode_vec3(bytes, &mut off)?);
            }
            extra_vectors.push((name, vecs));
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

    /// Parse a `RestartData` from an ASCII restart text.
    pub fn read_ascii(text: &str) -> Result<RestartData, String> {
        let mut version = String::new();
        let mut timestamp = 0u64;
        let mut step = 0u64;
        let mut time = 0.0f64;
        let mut crate_name = String::new();
        let mut description = String::new();
        let mut positions: Vec<[f64; 3]> = Vec::new();
        let mut velocities: Vec<[f64; 3]> = Vec::new();
        let mut forces: Vec<[f64; 3]> = Vec::new();
        let mut masses: Vec<f64> = Vec::new();
        let mut types: Vec<u32> = Vec::new();
        let mut box_matrix = [[0.0f64; 3]; 3];

        #[derive(PartialEq)]
        enum Section {
            None,
            Positions,
            Velocities,
            Forces,
            Masses,
            Types,
            Box,
        }
        let mut section = Section::None;
        let mut box_row = 0usize;

        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            match line {
                "BEGIN_POSITIONS" => {
                    section = Section::Positions;
                    continue;
                }
                "END_POSITIONS" => {
                    section = Section::None;
                    continue;
                }
                "BEGIN_VELOCITIES" => {
                    section = Section::Velocities;
                    continue;
                }
                "END_VELOCITIES" => {
                    section = Section::None;
                    continue;
                }
                "BEGIN_FORCES" => {
                    section = Section::Forces;
                    continue;
                }
                "END_FORCES" => {
                    section = Section::None;
                    continue;
                }
                "BEGIN_MASSES" => {
                    section = Section::Masses;
                    continue;
                }
                "END_MASSES" => {
                    section = Section::None;
                    continue;
                }
                "BEGIN_TYPES" => {
                    section = Section::Types;
                    continue;
                }
                "END_TYPES" => {
                    section = Section::None;
                    continue;
                }
                "BEGIN_BOX" => {
                    section = Section::Box;
                    box_row = 0;
                    continue;
                }
                "END_BOX" => {
                    section = Section::None;
                    continue;
                }
                _ => {}
            }
            match section {
                Section::Positions | Section::Velocities | Section::Forces => {
                    let nums: Vec<f64> = line
                        .split_whitespace()
                        .map(|s| s.parse::<f64>().unwrap_or(0.0))
                        .collect();
                    if nums.len() >= 3 {
                        let arr = [nums[0], nums[1], nums[2]];
                        match section {
                            Section::Positions => positions.push(arr),
                            Section::Velocities => velocities.push(arr),
                            Section::Forces => forces.push(arr),
                            _ => {}
                        }
                    }
                }
                Section::Masses => {
                    if let Ok(m) = line.parse::<f64>() {
                        masses.push(m);
                    }
                }
                Section::Types => {
                    if let Ok(t) = line.parse::<u32>() {
                        types.push(t);
                    }
                }
                Section::Box => {
                    if box_row < 3 {
                        let nums: Vec<f64> = line
                            .split_whitespace()
                            .map(|s| s.parse::<f64>().unwrap_or(0.0))
                            .collect();
                        if nums.len() >= 3 {
                            box_matrix[box_row] = [nums[0], nums[1], nums[2]];
                            box_row += 1;
                        }
                    }
                }
                Section::None => {
                    // key-value header lines
                    if let Some(rest) = line.strip_prefix("VERSION ") {
                        version = rest.trim().to_string();
                    } else if let Some(rest) = line.strip_prefix("TIMESTAMP ") {
                        timestamp = rest.trim().parse().unwrap_or(0);
                    } else if let Some(rest) = line.strip_prefix("STEP ") {
                        step = rest.trim().parse().unwrap_or(0);
                    } else if let Some(rest) = line.strip_prefix("TIME ") {
                        time = rest.trim().parse().unwrap_or(0.0);
                    } else if let Some(rest) = line.strip_prefix("CRATE ") {
                        crate_name = rest.trim().to_string();
                    } else if let Some(rest) = line.strip_prefix("DESC ") {
                        description = rest.trim().to_string();
                    }
                }
            }
        }
        let metadata = RestartMetadata {
            version,
            timestamp,
            step,
            time,
            crate_name,
            description,
        };
        Ok(RestartData {
            metadata,
            positions,
            velocities,
            forces,
            masses,
            types,
            box_matrix,
            extra_scalars: Vec::new(),
            extra_vectors: Vec::new(),
        })
    }
}

// ── CheckpointManager ─────────────────────────────────────────────────────────

/// Manages a rolling set of checkpoint files in a directory.
#[derive(Debug, Clone)]
pub struct CheckpointManager {
    /// Base directory where checkpoint files are stored.
    pub base_dir: String,
    /// Maximum number of checkpoints to retain.
    pub max_checkpoints: usize,
    /// In-memory list of `(step, file_path)` known checkpoints.
    checkpoints: Vec<(u64, String)>,
}

impl CheckpointManager {
    /// Create a new manager rooted at `base_dir`, keeping at most `max_checkpoints` files.
    pub fn new(base_dir: &str, max_checkpoints: usize) -> Self {
        Self {
            base_dir: base_dir.to_string(),
            max_checkpoints: max_checkpoints.max(1),
            checkpoints: Vec::new(),
        }
    }

    /// Serialize `data` and write it as `checkpoint_`step`.bin` under `base_dir`.
    ///
    /// Returns the full path of the newly written file.
    pub fn save_checkpoint(&mut self, data: &RestartData, step: u64) -> String {
        let filename = format!("{}/checkpoint_{step:010}.bin", self.base_dir);
        let writer = RestartWriter::new(&filename, RestartFormat::Binary);
        let _ = std::fs::create_dir_all(&self.base_dir);
        let _ = writer.write(data);
        self.checkpoints.push((step, filename.clone()));
        self.prune_old_checkpoints();
        filename
    }

    /// Load and return the most recent checkpoint, or `None` if none exist.
    pub fn load_latest(&self) -> Option<RestartData> {
        let (_, path) = self.checkpoints.last()?;
        let reader = RestartReader::new(path);
        reader.read().ok()
    }

    /// Return all known checkpoints as `(step, path)` pairs, oldest first.
    pub fn list_checkpoints(&self) -> Vec<(u64, String)> {
        self.checkpoints.clone()
    }

    /// Remove the oldest checkpoints so that at most `max_checkpoints` remain.
    pub fn prune_old_checkpoints(&mut self) {
        while self.checkpoints.len() > self.max_checkpoints {
            let (_, path) = self.checkpoints.remove(0);
            let _ = std::fs::remove_file(&path);
        }
    }
}

// ── IncrementalRestart ────────────────────────────────────────────────────────

/// Tracks which particle indices have changed since the last checkpoint.
///
/// Only changed particles are serialized, reducing I/O cost for large systems
/// with sparse updates.
#[derive(Debug, Clone, Default)]
pub struct IncrementalRestart {
    /// Indices of particles that have been marked as modified.
    pub changed: Vec<usize>,
}

impl IncrementalRestart {
    /// Create an empty tracker.
    pub fn new() -> Self {
        Self::default()
    }

    /// Mark particle `idx` as changed.
    pub fn mark_changed(&mut self, idx: usize) {
        if !self.changed.contains(&idx) {
            self.changed.push(idx);
        }
    }

    /// Clear the changed set (call after a checkpoint is written).
    pub fn reset(&mut self) {
        self.changed.clear();
    }

    /// Extract only the changed particles from `full` into a new `RestartData`.
    ///
    /// Metadata is copied from `full`; all other arrays are filtered to `changed` indices.
    pub fn extract_delta(&self, full: &RestartData) -> RestartData {
        let indices = &self.changed;
        let positions = indices
            .iter()
            .filter_map(|&i| full.positions.get(i).copied())
            .collect();
        let velocities = indices
            .iter()
            .filter_map(|&i| full.velocities.get(i).copied())
            .collect();
        let forces = indices
            .iter()
            .filter_map(|&i| full.forces.get(i).copied())
            .collect();
        let masses = indices
            .iter()
            .filter_map(|&i| full.masses.get(i).copied())
            .collect();
        let types = indices
            .iter()
            .filter_map(|&i| full.types.get(i).copied())
            .collect();
        RestartData {
            metadata: full.metadata.clone(),
            positions,
            velocities,
            forces,
            masses,
            types,
            box_matrix: full.box_matrix,
            extra_scalars: Vec::new(),
            extra_vectors: Vec::new(),
        }
    }
}

// ── RestartValidator ──────────────────────────────────────────────────────────

/// Computes and verifies simple checksums over restart binary blobs.
#[derive(Debug, Clone, Default)]
pub struct RestartValidator;

impl RestartValidator {
    /// Create a new validator.
    pub fn new() -> Self {
        Self
    }

    /// Compute a simple additive sum checksum over `bytes`.
    pub fn checksum_sum(bytes: &[u8]) -> u64 {
        bytes.iter().map(|&b| b as u64).sum()
    }

    /// Compute a byte-wise XOR checksum over `bytes`.
    pub fn checksum_xor(bytes: &[u8]) -> u8 {
        bytes.iter().fold(0u8, |acc, &b| acc ^ b)
    }

    /// Verify that `bytes` have not changed since `expected_sum` was computed.
    pub fn verify_sum(bytes: &[u8], expected_sum: u64) -> bool {
        Self::checksum_sum(bytes) == expected_sum
    }

    /// Verify that `bytes` have not changed since `expected_xor` was computed.
    pub fn verify_xor(bytes: &[u8], expected_xor: u8) -> bool {
        Self::checksum_xor(bytes) == expected_xor
    }
}

// ── Conversion utilities ──────────────────────────────────────────────────────

/// Convert a `RestartData` snapshot to an XYZ-format string.
///
/// The XYZ file contains a single frame with element symbol "XX" for all particles.
pub fn restart_to_xyz(data: &RestartData) -> String {
    let n = data.n_particles();
    let mut s = String::new();
    s.push_str(&format!("{n}\n"));
    s.push_str(&format!(
        "Restart step={} time={:.6}\n",
        data.metadata.step, data.metadata.time
    ));
    for (i, p) in data.positions.iter().enumerate() {
        let sym = if let Some(&t) = data.types.get(i) {
            match t {
                0 => "H",
                1 => "C",
                2 => "N",
                3 => "O",
                _ => "X",
            }
        } else {
            "X"
        };
        s.push_str(&format!("{sym} {:.6} {:.6} {:.6}\n", p[0], p[1], p[2]));
    }
    s
}

/// Convert a `RestartData` snapshot to a LAMMPS dump-format string.
///
/// Fields: `id type x y z vx vy vz fx fy fz mass`
pub fn restart_to_lammps_dump(data: &RestartData) -> String {
    let n = data.n_particles();
    let step = data.metadata.step;
    let mut s = String::new();
    s.push_str(&format!("ITEM: TIMESTEP\n{step}\n"));
    s.push_str(&format!("ITEM: NUMBER OF ATOMS\n{n}\n"));
    let bm = &data.box_matrix;
    s.push_str("ITEM: BOX BOUNDS pp pp pp\n");
    s.push_str(&format!(
        "0.0 {:.6}\n0.0 {:.6}\n0.0 {:.6}\n",
        bm[0][0], bm[1][1], bm[2][2]
    ));
    s.push_str("ITEM: ATOMS id type x y z vx vy vz fx fy fz mass\n");
    for i in 0..n {
        let p = data.positions.get(i).copied().unwrap_or([0.0; 3]);
        let v = data.velocities.get(i).copied().unwrap_or([0.0; 3]);
        let f = data.forces.get(i).copied().unwrap_or([0.0; 3]);
        let m = data.masses.get(i).copied().unwrap_or(1.0);
        let tp = data.types.get(i).copied().unwrap_or(0);
        s.push_str(&format!(
            "{} {} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6} {:.6}\n",
            i + 1,
            tp,
            p[0],
            p[1],
            p[2],
            v[0],
            v[1],
            v[2],
            f[0],
            f[1],
            f[2],
            m
        ));
    }
    s
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_data(n: usize) -> RestartData {
        let meta = RestartMetadata::new("1.0", 1234567890, 42, 4.2, "oxiphysics-md", "unit test");
        let positions = (0..n).map(|i| [i as f64, i as f64 * 0.5, 0.0]).collect();
        let velocities = (0..n).map(|i| [0.1 * i as f64, 0.0, 0.0]).collect();
        let forces = (0..n).map(|_| [0.0, -9.81, 0.0]).collect();
        let masses = (0..n).map(|_| 1.0).collect();
        let types = (0..n).map(|i| (i % 3) as u32).collect();
        let box_matrix = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]];
        RestartData {
            metadata: meta,
            positions,
            velocities,
            forces,
            masses,
            types,
            box_matrix,
            extra_scalars: Vec::new(),
            extra_vectors: Vec::new(),
        }
    }

    // -- RestartMetadata --

    #[test]
    fn test_metadata_new_fields() {
        let m = RestartMetadata::new("2.0", 100, 5, 0.5, "crate_x", "desc");
        assert_eq!(m.version, "2.0");
        assert_eq!(m.timestamp, 100);
        assert_eq!(m.step, 5);
        assert!((m.time - 0.5).abs() < 1e-12);
        assert_eq!(m.crate_name, "crate_x");
        assert_eq!(m.description, "desc");
    }

    #[test]
    fn test_metadata_default_test() {
        let m = RestartMetadata::default_test();
        assert_eq!(m.version, "1.0");
        assert_eq!(m.step, 0);
    }

    // -- RestartData --

    #[test]
    fn test_restart_data_n_particles() {
        let d = make_data(7);
        assert_eq!(d.n_particles(), 7);
    }

    #[test]
    fn test_restart_data_empty() {
        let d = RestartData::empty(RestartMetadata::default_test());
        assert_eq!(d.n_particles(), 0);
        assert!(d.positions.is_empty());
    }

    #[test]
    fn test_restart_data_single_particle_test() {
        let d = RestartData::single_particle_test();
        assert_eq!(d.n_particles(), 1);
        assert!((d.positions[0][0] - 1.0).abs() < 1e-12);
    }

    // -- Binary round-trip --

    #[test]
    fn test_binary_roundtrip_zero_particles() {
        let d = make_data(0);
        let bytes = RestartWriter::write_binary(&d);
        let recovered = RestartReader::read_binary(&bytes).unwrap();
        assert_eq!(recovered.n_particles(), 0);
        assert_eq!(recovered.metadata.step, 42);
    }

    #[test]
    fn test_binary_roundtrip_positions() {
        let d = make_data(5);
        let bytes = RestartWriter::write_binary(&d);
        let r = RestartReader::read_binary(&bytes).unwrap();
        for i in 0..5 {
            assert!((r.positions[i][0] - d.positions[i][0]).abs() < 1e-12);
            assert!((r.positions[i][1] - d.positions[i][1]).abs() < 1e-12);
        }
    }

    #[test]
    fn test_binary_roundtrip_metadata() {
        let d = make_data(3);
        let bytes = RestartWriter::write_binary(&d);
        let r = RestartReader::read_binary(&bytes).unwrap();
        assert_eq!(r.metadata.version, "1.0");
        assert_eq!(r.metadata.timestamp, 1234567890);
        assert_eq!(r.metadata.step, 42);
        assert!((r.metadata.time - 4.2).abs() < 1e-10);
        assert_eq!(r.metadata.crate_name, "oxiphysics-md");
    }

    #[test]
    fn test_binary_roundtrip_box_matrix() {
        let d = make_data(2);
        let bytes = RestartWriter::write_binary(&d);
        let r = RestartReader::read_binary(&bytes).unwrap();
        assert!((r.box_matrix[0][0] - 10.0).abs() < 1e-12);
        assert!((r.box_matrix[1][1] - 10.0).abs() < 1e-12);
        assert!((r.box_matrix[2][2] - 10.0).abs() < 1e-12);
        assert!(r.box_matrix[0][1].abs() < 1e-12);
    }

    #[test]
    fn test_binary_roundtrip_types() {
        let d = make_data(6);
        let bytes = RestartWriter::write_binary(&d);
        let r = RestartReader::read_binary(&bytes).unwrap();
        for i in 0..6 {
            assert_eq!(r.types[i], (i % 3) as u32);
        }
    }

    #[test]
    fn test_binary_roundtrip_extra_scalars() {
        let mut d = make_data(3);
        d.extra_scalars
            .push(("charge".into(), vec![0.1, -0.2, 0.3]));
        let bytes = RestartWriter::write_binary(&d);
        let r = RestartReader::read_binary(&bytes).unwrap();
        assert_eq!(r.extra_scalars.len(), 1);
        assert_eq!(r.extra_scalars[0].0, "charge");
        assert!((r.extra_scalars[0].1[1] - (-0.2)).abs() < 1e-12);
    }

    #[test]
    fn test_binary_roundtrip_extra_vectors() {
        let mut d = make_data(2);
        d.extra_vectors
            .push(("spin".into(), vec![[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]]));
        let bytes = RestartWriter::write_binary(&d);
        let r = RestartReader::read_binary(&bytes).unwrap();
        assert_eq!(r.extra_vectors.len(), 1);
        assert_eq!(r.extra_vectors[0].0, "spin");
        assert!((r.extra_vectors[0].1[0][2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_binary_magic_bytes() {
        let d = make_data(1);
        let bytes = RestartWriter::write_binary(&d);
        assert_eq!(&bytes[0..4], b"OXRS");
    }

    // -- ASCII round-trip --

    #[test]
    fn test_ascii_roundtrip_basic() {
        let d = make_data(4);
        let text = RestartWriter::write_ascii(&d);
        let r = RestartReader::read_ascii(&text).unwrap();
        assert_eq!(r.n_particles(), 4);
        assert_eq!(r.metadata.step, 42);
    }

    #[test]
    fn test_ascii_roundtrip_positions() {
        let d = make_data(3);
        let text = RestartWriter::write_ascii(&d);
        let r = RestartReader::read_ascii(&text).unwrap();
        for i in 0..3 {
            assert!((r.positions[i][0] - d.positions[i][0]).abs() < 1e-4);
        }
    }

    #[test]
    fn test_ascii_roundtrip_box() {
        let d = make_data(1);
        let text = RestartWriter::write_ascii(&d);
        let r = RestartReader::read_ascii(&text).unwrap();
        assert!((r.box_matrix[0][0] - 10.0).abs() < 1e-4);
    }

    #[test]
    fn test_ascii_contains_keywords() {
        let d = make_data(2);
        let text = RestartWriter::write_ascii(&d);
        assert!(text.contains("BEGIN_POSITIONS"));
        assert!(text.contains("END_POSITIONS"));
        assert!(text.contains("VERSION"));
        assert!(text.contains("STEP"));
    }

    // -- JSON --

    #[test]
    fn test_json_contains_particles_key() {
        let d = make_data(3);
        let j = RestartWriter::write_json(&d);
        assert!(j.contains("n_particles"));
        assert!(j.contains("positions"));
        assert!(j.contains("velocities"));
    }

    #[test]
    fn test_json_step_present() {
        let d = make_data(1);
        let j = RestartWriter::write_json(&d);
        assert!(j.contains("\"step\": 42"));
    }

    // -- detect_format --

    #[test]
    fn test_detect_format_binary() {
        let d = make_data(1);
        let bytes = RestartWriter::write_binary(&d);
        assert_eq!(RestartReader::detect_format(&bytes), RestartFormat::Binary);
    }

    #[test]
    fn test_detect_format_ascii() {
        let text = "# OxiPhysics restart\nVERSION 1.0\n";
        assert_eq!(
            RestartReader::detect_format(text.as_bytes()),
            RestartFormat::Ascii
        );
    }

    #[test]
    fn test_detect_format_json() {
        let j = "{\"n_particles\": 0}";
        assert_eq!(
            RestartReader::detect_format(j.as_bytes()),
            RestartFormat::Json
        );
    }

    // -- Validator --

    #[test]
    fn test_validator_sum_consistent() {
        let bytes = b"hello world";
        let sum = RestartValidator::checksum_sum(bytes);
        assert!(RestartValidator::verify_sum(bytes, sum));
    }

    #[test]
    fn test_validator_xor_consistent() {
        let bytes = b"test data";
        let xor = RestartValidator::checksum_xor(bytes);
        assert!(RestartValidator::verify_xor(bytes, xor));
    }

    #[test]
    fn test_validator_sum_detects_corruption() {
        let bytes = b"original";
        let sum = RestartValidator::checksum_sum(bytes);
        let corrupt = b"0riginal";
        assert!(!RestartValidator::verify_sum(corrupt, sum));
    }

    #[test]
    fn test_validator_xor_empty() {
        let xor = RestartValidator::checksum_xor(b"");
        assert_eq!(xor, 0);
    }

    // -- IncrementalRestart --

    #[test]
    fn test_incremental_mark_and_extract() {
        let full = make_data(5);
        let mut inc = IncrementalRestart::new();
        inc.mark_changed(1);
        inc.mark_changed(3);
        let delta = inc.extract_delta(&full);
        assert_eq!(delta.n_particles(), 2);
        assert!((delta.positions[0][0] - full.positions[1][0]).abs() < 1e-12);
    }

    #[test]
    fn test_incremental_reset_clears_changes() {
        let mut inc = IncrementalRestart::new();
        inc.mark_changed(0);
        inc.mark_changed(2);
        inc.reset();
        assert!(inc.changed.is_empty());
    }

    #[test]
    fn test_incremental_no_duplicates() {
        let mut inc = IncrementalRestart::new();
        inc.mark_changed(0);
        inc.mark_changed(0);
        inc.mark_changed(0);
        assert_eq!(inc.changed.len(), 1);
    }

    // -- Conversion utilities --

    #[test]
    fn test_restart_to_xyz_header_line_count() {
        let d = make_data(4);
        let xyz = restart_to_xyz(&d);
        let lines: Vec<&str> = xyz.lines().collect();
        assert!(lines.len() >= 6); // count + comment + 4 atoms
        assert_eq!(lines[0], "4");
    }

    #[test]
    fn test_restart_to_lammps_dump_has_timestep() {
        let d = make_data(2);
        let dump = restart_to_lammps_dump(&d);
        assert!(dump.contains("ITEM: TIMESTEP"));
        assert!(dump.contains("ITEM: NUMBER OF ATOMS"));
        assert!(dump.contains("ITEM: ATOMS"));
    }

    #[test]
    fn test_restart_to_lammps_dump_atom_count() {
        let d = make_data(3);
        let dump = restart_to_lammps_dump(&d);
        let atom_lines: Vec<&str> = dump
            .lines()
            .skip_while(|l| !l.starts_with("ITEM: ATOMS"))
            .skip(1)
            .collect();
        assert_eq!(atom_lines.len(), 3);
    }

    // -- CheckpointManager (in-memory only, uses /tmp) --

    #[test]
    fn test_checkpoint_manager_list_empty() {
        let dir = std::env::temp_dir().join("oxi_test_ckpt_empty");
        let mgr = CheckpointManager::new(dir.to_str().unwrap_or(""), 3);
        assert!(mgr.list_checkpoints().is_empty());
    }

    #[test]
    fn test_checkpoint_manager_prune_keeps_max() {
        let dir = std::env::temp_dir().join("oxi_test_ckpt_prune");
        let _ = std::fs::remove_dir_all(&dir);
        let mut mgr = CheckpointManager::new(dir.to_str().unwrap_or(""), 2);
        let d = make_data(1);
        mgr.save_checkpoint(&d, 1);
        mgr.save_checkpoint(&d, 2);
        mgr.save_checkpoint(&d, 3);
        assert_eq!(mgr.list_checkpoints().len(), 2);
        // Latest should be step 3.
        let latest = mgr.list_checkpoints().last().unwrap().0;
        assert_eq!(latest, 3);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
