// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Parallel I/O utilities for large-scale distributed physics simulations.
//!
//! This module provides:
//! - **Domain decomposition metadata** – rank, ghost cells, partition boundaries
//! - **Chunked file I/O** – fixed-size chunks with parallel writer state tracking
//! - **Parallel VTK output** – PVTU / PVD index files for distributed meshes
//! - **MPI-like rank / domain storage** – load-balanced partition export
//! - **Ghost cell metadata** – halo layers for stencil communication
//! - **Distributed array serialization** – flat binary layout for parallel restart
//! - **Checkpoint / restart** – step-indexed checkpoint files
//! - **Merge of distributed outputs** – reassembling per-rank files
//! - **Parallel HDF5-like layout** – dataset headers + contiguous data blocks

use std::fmt;
use std::fs;
use std::io::{self, BufRead, Write};

// ─────────────────────────────────────────────────────────────────────────────
// Domain decomposition metadata
// ─────────────────────────────────────────────────────────────────────────────

/// Axis-aligned bounding box in 3-D used to describe a subdomain.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainBox {
    /// Minimum corner `[x_min, y_min, z_min]`.
    pub min: [f64; 3],
    /// Maximum corner `[x_max, y_max, z_max]`.
    pub max: [f64; 3],
}

impl DomainBox {
    /// Construct a new domain box.
    pub fn new(min: [f64; 3], max: [f64; 3]) -> Self {
        Self { min, max }
    }

    /// Volume of the bounding box.
    pub fn volume(&self) -> f64 {
        (self.max[0] - self.min[0]).max(0.0)
            * (self.max[1] - self.min[1]).max(0.0)
            * (self.max[2] - self.min[2]).max(0.0)
    }

    /// Return the centre of the box.
    pub fn centre(&self) -> [f64; 3] {
        [
            0.5 * (self.min[0] + self.max[0]),
            0.5 * (self.min[1] + self.max[1]),
            0.5 * (self.min[2] + self.max[2]),
        ]
    }

    /// Test whether the point `p` lies strictly inside the box.
    pub fn contains(&self, p: [f64; 3]) -> bool {
        p[0] > self.min[0]
            && p[0] < self.max[0]
            && p[1] > self.min[1]
            && p[1] < self.max[1]
            && p[2] > self.min[2]
            && p[2] < self.max[2]
    }

    /// Return `true` if this box overlaps (strictly) with `other`.
    pub fn overlaps(&self, other: &DomainBox) -> bool {
        self.min[0] < other.max[0]
            && self.max[0] > other.min[0]
            && self.min[1] < other.max[1]
            && self.max[1] > other.min[1]
            && self.min[2] < other.max[2]
            && self.max[2] > other.min[2]
    }
}

/// Per-rank subdomain metadata for a domain-decomposed simulation.
#[derive(Debug, Clone)]
pub struct RankDomain {
    /// MPI rank index (0-based).
    pub rank: usize,
    /// Total number of ranks in the communicator.
    pub n_ranks: usize,
    /// Spatial bounding box of the owned (non-ghost) region.
    pub owned_box: DomainBox,
    /// Number of ghost cell layers on each face `[x-, x+, y-, y+, z-, z+]`.
    pub ghost_layers: [usize; 6],
    /// Number of owned cells in each direction `[nx, ny, nz]`.
    pub owned_cells: [usize; 3],
    /// Global offset of this rank's first cell `[ox, oy, oz]`.
    pub global_offset: [usize; 3],
}

impl RankDomain {
    /// Construct domain metadata for a single rank.
    pub fn new(
        rank: usize,
        n_ranks: usize,
        owned_box: DomainBox,
        ghost_layers: [usize; 6],
        owned_cells: [usize; 3],
        global_offset: [usize; 3],
    ) -> Self {
        Self {
            rank,
            n_ranks,
            owned_box,
            ghost_layers,
            owned_cells,
            global_offset,
        }
    }

    /// Total number of locally owned cells.
    pub fn n_owned_cells(&self) -> usize {
        self.owned_cells[0] * self.owned_cells[1] * self.owned_cells[2]
    }

    /// Total cells including ghost layers (all six faces).
    pub fn n_total_cells(&self) -> usize {
        let ext_x = self.owned_cells[0] + self.ghost_layers[0] + self.ghost_layers[1];
        let ext_y = self.owned_cells[1] + self.ghost_layers[2] + self.ghost_layers[3];
        let ext_z = self.owned_cells[2] + self.ghost_layers[4] + self.ghost_layers[5];
        ext_x * ext_y * ext_z
    }

    /// Fraction of the global workload assigned to this rank.
    ///
    /// Assumes uniform distribution; override if load varies.
    pub fn load_fraction(&self) -> f64 {
        if self.n_ranks == 0 {
            return 1.0;
        }
        self.n_owned_cells() as f64
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Ghost cell metadata
// ─────────────────────────────────────────────────────────────────────────────

/// Cardinal face directions for a 3-D structured mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Face {
    /// Negative-X face.
    XMinus,
    /// Positive-X face.
    XPlus,
    /// Negative-Y face.
    YMinus,
    /// Positive-Y face.
    YPlus,
    /// Negative-Z face.
    ZMinus,
    /// Positive-Z face.
    ZPlus,
}

impl fmt::Display for Face {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            Face::XMinus => "XMinus",
            Face::XPlus => "XPlus",
            Face::YMinus => "YMinus",
            Face::YPlus => "YPlus",
            Face::ZMinus => "ZMinus",
            Face::ZPlus => "ZPlus",
        };
        write!(f, "{s}")
    }
}

/// Describes the ghost cell layer on a single face.
#[derive(Debug, Clone)]
pub struct GhostLayer {
    /// Which face this ghost layer belongs to.
    pub face: Face,
    /// Owning rank on the other side of this interface.
    pub neighbour_rank: usize,
    /// Number of ghost cell layers.
    pub depth: usize,
    /// Number of cells in the two tangential directions `[t1, t2]`.
    pub tangential_cells: [usize; 2],
}

impl GhostLayer {
    /// Total number of ghost cells in this layer.
    pub fn n_cells(&self) -> usize {
        self.depth * self.tangential_cells[0] * self.tangential_cells[1]
    }
}

/// Full set of ghost-cell metadata for one rank.
#[derive(Debug, Clone, Default)]
pub struct GhostCellMeta {
    /// One entry per active face interface.
    pub layers: Vec<GhostLayer>,
}

impl GhostCellMeta {
    /// Create an empty metadata container.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a ghost layer descriptor.
    pub fn add_layer(&mut self, layer: GhostLayer) {
        self.layers.push(layer);
    }

    /// Total ghost cells summed across all faces.
    pub fn total_ghost_cells(&self) -> usize {
        self.layers.iter().map(|l| l.n_cells()).sum()
    }

    /// List all neighbour ranks (deduplicated).
    pub fn neighbour_ranks(&self) -> Vec<usize> {
        let mut ranks: Vec<usize> = self.layers.iter().map(|l| l.neighbour_rank).collect();
        ranks.sort_unstable();
        ranks.dedup();
        ranks
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Chunk splitting
// ─────────────────────────────────────────────────────────────────────────────

/// Specification for a single data chunk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChunkSpec {
    /// Inclusive start index of this chunk.
    pub start: usize,
    /// Exclusive end index of this chunk.
    pub end: usize,
    /// Zero-based chunk index.
    pub chunk_idx: usize,
}

impl ChunkSpec {
    /// Number of elements in this chunk.
    pub fn size(&self) -> usize {
        self.end.saturating_sub(self.start)
    }
}

/// Divide `total` elements into at most `n_chunks` contiguous chunks.
///
/// Returns a `Vec`ChunkSpec` covering `\[0, total)`. The last chunk absorbs
/// any remainder so that all elements are covered.
pub fn split_into_chunks(total: usize, n_chunks: usize) -> Vec<ChunkSpec> {
    if total == 0 || n_chunks == 0 {
        return Vec::new();
    }
    let actual_chunks = n_chunks.min(total);
    let base = total / actual_chunks;
    let remainder = total % actual_chunks;
    let mut specs = Vec::with_capacity(actual_chunks);
    let mut start = 0;
    for idx in 0..actual_chunks {
        let extra = if idx < remainder { 1 } else { 0 };
        let end = start + base + extra;
        specs.push(ChunkSpec {
            start,
            end,
            chunk_idx: idx,
        });
        start = end;
    }
    specs
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel writer state
// ─────────────────────────────────────────────────────────────────────────────

/// Tracks parallel write progress across multiple chunk files.
#[derive(Debug, Clone)]
pub struct ParallelWriter {
    /// Base filename (without chunk suffix).
    pub filename: String,
    /// Number of chunks that have been written so far.
    pub chunks_written: usize,
    /// Total number of chunks expected.
    pub total_chunks: usize,
}

impl ParallelWriter {
    /// Create a new `ParallelWriter` for `total_chunks` chunks.
    pub fn new(filename: impl Into<String>, total_chunks: usize) -> Self {
        Self {
            filename: filename.into(),
            chunks_written: 0,
            total_chunks,
        }
    }

    /// Record that one more chunk has been written.
    pub fn mark_chunk_written(&mut self) {
        if self.chunks_written < self.total_chunks {
            self.chunks_written += 1;
        }
    }

    /// Returns `true` when all chunks have been written.
    pub fn is_complete(&self) -> bool {
        self.chunks_written >= self.total_chunks
    }

    /// Fraction of chunks completed, in `\[0.0, 1.0\]`.
    pub fn progress(&self) -> f64 {
        if self.total_chunks == 0 {
            return 1.0;
        }
        self.chunks_written as f64 / self.total_chunks as f64
    }

    /// Reset the writer to zero written chunks.
    pub fn reset(&mut self) {
        self.chunks_written = 0;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel VTK output — PVTU / PVD
// ─────────────────────────────────────────────────────────────────────────────

/// Metadata for a single piece in a parallel VTU (PVTU) file.
#[derive(Debug, Clone)]
pub struct PvtuPiece {
    /// Relative path to the per-rank `.vtu` file.
    pub file: String,
    /// MPI rank that owns this piece.
    pub rank: usize,
}

/// Generates the XML content of a `.pvtu` (Parallel VTK Unstructured) index file.
///
/// * `pieces` – per-rank VTU file references.
/// * `point_arrays` – names of point data arrays present in each piece.
/// * `cell_arrays` – names of cell data arrays present in each piece.
pub fn write_pvtu_xml(pieces: &[PvtuPiece], point_arrays: &[&str], cell_arrays: &[&str]) -> String {
    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\"?>\n");
    s.push_str(
        "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n",
    );
    s.push_str("  <PUnstructuredGrid GhostLevel=\"1\">\n");
    s.push_str("    <PPoints>\n");
    s.push_str("      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n");
    s.push_str("    </PPoints>\n");
    if !point_arrays.is_empty() {
        s.push_str("    <PPointData>\n");
        for name in point_arrays {
            s.push_str(&format!(
                "      <PDataArray type=\"Float64\" Name=\"{name}\"/>\n"
            ));
        }
        s.push_str("    </PPointData>\n");
    }
    if !cell_arrays.is_empty() {
        s.push_str("    <PCellData>\n");
        for name in cell_arrays {
            s.push_str(&format!(
                "      <PDataArray type=\"Float64\" Name=\"{name}\"/>\n"
            ));
        }
        s.push_str("    </PCellData>\n");
    }
    for piece in pieces {
        s.push_str(&format!("    <Piece Source=\"{}\"/>\n", piece.file));
    }
    s.push_str("  </PUnstructuredGrid>\n");
    s.push_str("</VTKFile>\n");
    s
}

/// Generates the XML content of a `.pvd` (ParaView Data) time-series file.
///
/// * `entries` – list of `(time_value, pvtu_file)` pairs.
pub fn write_pvd_xml(entries: &[(f64, &str)]) -> String {
    let mut s = String::new();
    s.push_str("<?xml version=\"1.0\"?>\n");
    s.push_str("<VTKFile type=\"Collection\" version=\"0.1\">\n");
    s.push_str("  <Collection>\n");
    for (t, file) in entries {
        s.push_str(&format!(
            "    <DataSet timestep=\"{t:.6e}\" group=\"\" part=\"0\" file=\"{file}\"/>\n"
        ));
    }
    s.push_str("  </Collection>\n");
    s.push_str("</VTKFile>\n");
    s
}

// ─────────────────────────────────────────────────────────────────────────────
// Load-balanced partition export
// ─────────────────────────────────────────────────────────────────────────────

/// Partition assignment for load-balanced domain decomposition.
#[derive(Debug, Clone)]
pub struct PartitionMap {
    /// For each global cell index, the owning rank.
    pub cell_to_rank: Vec<usize>,
    /// Total number of ranks.
    pub n_ranks: usize,
}

impl PartitionMap {
    /// Create a `PartitionMap` from a flat assignment vector.
    pub fn new(cell_to_rank: Vec<usize>, n_ranks: usize) -> Self {
        Self {
            cell_to_rank,
            n_ranks,
        }
    }

    /// Number of global cells.
    pub fn n_cells(&self) -> usize {
        self.cell_to_rank.len()
    }

    /// Number of cells assigned to the given rank.
    pub fn cells_on_rank(&self, rank: usize) -> usize {
        self.cell_to_rank.iter().filter(|&&r| r == rank).count()
    }

    /// Load imbalance factor: `max_cells / ideal_cells`.
    ///
    /// Returns `1.0` for a perfectly balanced partition; higher values indicate
    /// greater imbalance.  Returns `1.0` if there are no cells or ranks.
    pub fn imbalance(&self) -> f64 {
        if self.n_cells() == 0 || self.n_ranks == 0 {
            return 1.0;
        }
        let ideal = self.n_cells() as f64 / self.n_ranks as f64;
        let max_load = (0..self.n_ranks)
            .map(|r| self.cells_on_rank(r))
            .max()
            .unwrap_or(0) as f64;
        if ideal < 1e-12 { 1.0 } else { max_load / ideal }
    }

    /// Build a simple round-robin partition for `n_cells` cells over `n_ranks` ranks.
    pub fn round_robin(n_cells: usize, n_ranks: usize) -> Self {
        if n_ranks == 0 {
            return Self::new(vec![], 0);
        }
        let map = (0..n_cells).map(|i| i % n_ranks).collect();
        Self::new(map, n_ranks)
    }

    /// Serialize to a simple text format (one rank per line).
    pub fn to_text(&self) -> String {
        let mut s = format!("n_ranks={}\n", self.n_ranks);
        for (i, &r) in self.cell_to_rank.iter().enumerate() {
            s.push_str(&format!("{i} {r}\n"));
        }
        s
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Distributed array serialization
// ─────────────────────────────────────────────────────────────────────────────

/// A distributed scalar array owned by one rank.
#[derive(Debug, Clone)]
pub struct DistributedArray {
    /// MPI rank that owns this shard.
    pub rank: usize,
    /// Global length of the full array.
    pub global_len: usize,
    /// Global start index of this shard.
    pub offset: usize,
    /// Local data values.
    pub data: Vec<f64>,
}

impl DistributedArray {
    /// Create a new shard.
    pub fn new(rank: usize, global_len: usize, offset: usize, data: Vec<f64>) -> Self {
        Self {
            rank,
            global_len,
            offset,
            data,
        }
    }

    /// Local length of this shard.
    pub fn local_len(&self) -> usize {
        self.data.len()
    }

    /// Serialise to a compact binary-hex text format.
    ///
    /// Format: `rank offset global_len n_values <hex-encoded f64s>\n`
    pub fn serialize_text(&self) -> String {
        let hex: String = self
            .data
            .iter()
            .flat_map(|&v| v.to_le_bytes())
            .map(|b| format!("{b:02x}"))
            .collect();
        format!(
            "{} {} {} {} {}\n",
            self.rank,
            self.offset,
            self.global_len,
            self.data.len(),
            hex
        )
    }

    /// Deserialise from the text format produced by `serialize_text`.
    pub fn deserialize_text(s: &str) -> Option<Self> {
        let mut it = s.split_whitespace();
        let rank: usize = it.next()?.parse().ok()?;
        let offset: usize = it.next()?.parse().ok()?;
        let global_len: usize = it.next()?.parse().ok()?;
        let n: usize = it.next()?.parse().ok()?;
        let hex = it.next().unwrap_or("");
        if hex.len() != n * 16 {
            return None;
        }
        let data: Vec<f64> = (0..n)
            .map(|i| {
                let chunk = &hex[i * 16..(i + 1) * 16];
                let mut bytes = [0u8; 8];
                for (j, byte) in bytes.iter_mut().enumerate() {
                    *byte = u8::from_str_radix(&chunk[j * 2..j * 2 + 2], 16).unwrap_or(0);
                }
                f64::from_le_bytes(bytes)
            })
            .collect();
        Some(Self {
            rank,
            global_len,
            offset,
            data,
        })
    }

    /// Reconstruct the global array by assembling shards in offset order.
    ///
    /// All shards must cover the entire global range without gaps or overlaps.
    pub fn assemble(shards: &[DistributedArray]) -> Option<Vec<f64>> {
        if shards.is_empty() {
            return Some(vec![]);
        }
        let global_len = shards[0].global_len;
        let mut result = vec![0.0_f64; global_len];
        for shard in shards {
            let end = shard.offset + shard.local_len();
            if end > global_len {
                return None;
            }
            result[shard.offset..end].copy_from_slice(&shard.data);
        }
        Some(result)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Checkpoint / restart
// ─────────────────────────────────────────────────────────────────────────────

/// Header record written at the start of every checkpoint file.
#[derive(Debug, Clone)]
pub struct CheckpointHeader {
    /// Simulation step index.
    pub step: u64,
    /// Simulation time.
    pub time: f64,
    /// Number of MPI ranks used when the checkpoint was written.
    pub n_ranks: usize,
    /// Number of particles/cells in this checkpoint.
    pub n_items: usize,
    /// A user-defined tag (simulation name, version, etc.).
    pub tag: String,
}

impl CheckpointHeader {
    /// Create a new checkpoint header.
    pub fn new(step: u64, time: f64, n_ranks: usize, n_items: usize, tag: &str) -> Self {
        Self {
            step,
            time,
            n_ranks,
            n_items,
            tag: tag.to_owned(),
        }
    }

    /// Serialize to a single-line text format.
    pub fn to_line(&self) -> String {
        format!(
            "CHECKPOINT step={} time={:.10e} n_ranks={} n_items={} tag={}\n",
            self.step, self.time, self.n_ranks, self.n_items, self.tag
        )
    }

    /// Parse a checkpoint header from a text line.
    pub fn from_line(line: &str) -> Option<Self> {
        if !line.trim_start().starts_with("CHECKPOINT") {
            return None;
        }
        let mut step = 0u64;
        let mut time = 0.0f64;
        let mut n_ranks = 1usize;
        let mut n_items = 0usize;
        let mut tag = String::new();
        for token in line.split_whitespace().skip(1) {
            if let Some(rest) = token.strip_prefix("step=") {
                step = rest.parse().unwrap_or(0);
            } else if let Some(rest) = token.strip_prefix("time=") {
                time = rest.parse().unwrap_or(0.0);
            } else if let Some(rest) = token.strip_prefix("n_ranks=") {
                n_ranks = rest.parse().unwrap_or(1);
            } else if let Some(rest) = token.strip_prefix("n_items=") {
                n_items = rest.parse().unwrap_or(0);
            } else if let Some(rest) = token.strip_prefix("tag=") {
                tag = rest.to_owned();
            }
        }
        Some(Self {
            step,
            time,
            n_ranks,
            n_items,
            tag,
        })
    }
}

/// Write a checkpoint file containing a header and particle positions.
///
/// The file is created at `path`. Each position is written as three
/// space-separated f64 values, one atom per line.
pub fn write_checkpoint(
    path: &str,
    header: &CheckpointHeader,
    positions: &[[f64; 3]],
) -> Result<(), io::Error> {
    let mut file = fs::File::create(path)?;
    write!(file, "{}", header.to_line())?;
    for pos in positions {
        writeln!(file, "{:.10e} {:.10e} {:.10e}", pos[0], pos[1], pos[2])?;
    }
    Ok(())
}

/// Read a checkpoint file written by [`write_checkpoint`].
///
/// Returns the header and a vector of positions.
pub fn read_checkpoint(path: &str) -> Result<(CheckpointHeader, Vec<[f64; 3]>), io::Error> {
    let file = fs::File::open(path)?;
    let reader = io::BufReader::new(file);
    let mut lines_iter = reader.lines();
    let header_line = lines_iter
        .next()
        .ok_or_else(|| io::Error::new(io::ErrorKind::UnexpectedEof, "empty checkpoint"))??;
    let header = CheckpointHeader::from_line(&header_line)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "bad checkpoint header"))?;
    let mut positions = Vec::new();
    for line in lines_iter {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() >= 3 {
            let x: f64 = parts[0].parse().unwrap_or(0.0);
            let y: f64 = parts[1].parse().unwrap_or(0.0);
            let z: f64 = parts[2].parse().unwrap_or(0.0);
            positions.push([x, y, z]);
        }
    }
    Ok((header, positions))
}

// ─────────────────────────────────────────────────────────────────────────────
// File merging
// ─────────────────────────────────────────────────────────────────────────────

/// Concatenate `n_chunks` chunk files into a single output file.
///
/// Chunk files are expected at `{base_path}.chunk{i}` for `i` in `0..n_chunks`.
/// Each chunk file is appended in order to `output_path`.
pub fn merge_chunk_files(
    base_path: &str,
    n_chunks: usize,
    output_path: &str,
) -> Result<(), io::Error> {
    let mut out = fs::File::create(output_path)?;
    for i in 0..n_chunks {
        let chunk_path = format!("{base_path}.chunk{i}");
        let data = fs::read(&chunk_path)?;
        out.write_all(&data)?;
    }
    Ok(())
}

/// Merge per-rank checkpoint position files into a single file.
///
/// Per-rank files are expected at `{base_path}.rank{r}` for `r` in `0..n_ranks`.
pub fn merge_rank_files(
    base_path: &str,
    n_ranks: usize,
    output_path: &str,
) -> Result<(), io::Error> {
    let mut out = fs::File::create(output_path)?;
    for r in 0..n_ranks {
        let rank_path = format!("{base_path}.rank{r}");
        if let Ok(data) = fs::read(&rank_path) {
            out.write_all(&data)?;
        }
    }
    Ok(())
}

// ─────────────────────────────────────────────────────────────────────────────
// Bandwidth estimation
// ─────────────────────────────────────────────────────────────────────────────

/// Estimate I/O bandwidth in MB/s.
///
/// * `bytes` — number of bytes transferred.
/// * `elapsed_secs` — elapsed time in seconds.
pub fn estimate_io_bandwidth(bytes: usize, elapsed_secs: f64) -> f64 {
    if elapsed_secs <= 0.0 {
        return 0.0;
    }
    (bytes as f64) / elapsed_secs / 1e6
}

// ─────────────────────────────────────────────────────────────────────────────
// Chunked trajectory
// ─────────────────────────────────────────────────────────────────────────────

/// A trajectory divided into fixed-size chunks of frames.
#[derive(Debug, Clone)]
pub struct ChunkedTrajectory {
    /// Total number of frames in the full trajectory.
    pub n_frames: usize,
    /// Number of atoms per frame.
    pub n_atoms: usize,
    /// Number of frames per chunk.
    pub chunk_size: usize,
    /// All frames stored as `Vec<\[f64; 3\]>` (one entry per atom).
    pub frames: Vec<Vec<[f64; 3]>>,
}

impl ChunkedTrajectory {
    /// Create an empty trajectory with the given parameters.
    pub fn new(n_atoms: usize, chunk_size: usize) -> Self {
        Self {
            n_frames: 0,
            n_atoms,
            chunk_size: chunk_size.max(1),
            frames: Vec::new(),
        }
    }

    /// Append a frame to the trajectory.
    pub fn push_frame(&mut self, frame: Vec<[f64; 3]>) {
        self.frames.push(frame);
        self.n_frames = self.frames.len();
    }

    /// Number of chunks needed to cover all frames.
    pub fn n_chunks(&self) -> usize {
        if self.n_frames == 0 {
            return 0;
        }
        self.n_frames.div_ceil(self.chunk_size)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// XYZ chunk I/O
// ─────────────────────────────────────────────────────────────────────────────

/// Write a `ChunkedTrajectory` to XYZ files, one file per chunk.
///
/// Output files are named `{base_path}.chunk{i}.xyz`.
pub fn write_chunked_xyz(traj: &ChunkedTrajectory, base_path: &str) -> Result<(), io::Error> {
    let n_chunks = traj.n_chunks();
    for chunk_idx in 0..n_chunks {
        let start_frame = chunk_idx * traj.chunk_size;
        let end_frame = (start_frame + traj.chunk_size).min(traj.n_frames);
        let chunk_path = format!("{base_path}.chunk{chunk_idx}.xyz");
        let mut file = fs::File::create(&chunk_path)?;
        for frame_idx in start_frame..end_frame {
            let frame = &traj.frames[frame_idx];
            writeln!(file, "{}", frame.len())?;
            writeln!(file, "Frame {frame_idx}")?;
            for pos in frame {
                writeln!(file, "X {:.6} {:.6} {:.6}", pos[0], pos[1], pos[2])?;
            }
        }
    }
    Ok(())
}

/// Read chunked XYZ files written by [`write_chunked_xyz`].
///
/// Reads `{base_path}.chunk{i}.xyz` for `i` in `0..n_chunks`.
pub fn read_chunked_xyz(base_path: &str, n_chunks: usize) -> Result<ChunkedTrajectory, io::Error> {
    let mut frames: Vec<Vec<[f64; 3]>> = Vec::new();
    for chunk_idx in 0..n_chunks {
        let chunk_path = format!("{base_path}.chunk{chunk_idx}.xyz");
        let file = fs::File::open(&chunk_path)?;
        let reader = io::BufReader::new(file);
        let lines: Vec<String> = reader.lines().collect::<Result<Vec<_>, _>>()?;
        let mut idx = 0;
        while idx < lines.len() {
            let n_atoms: usize = lines[idx].trim().parse().unwrap_or(0);
            idx += 1;
            if idx >= lines.len() {
                break;
            }
            idx += 1; // skip comment line
            let mut frame = Vec::with_capacity(n_atoms);
            for _ in 0..n_atoms {
                if idx >= lines.len() {
                    break;
                }
                let parts: Vec<&str> = lines[idx].split_whitespace().collect();
                if parts.len() >= 4 {
                    let x: f64 = parts[1].parse().unwrap_or(0.0);
                    let y: f64 = parts[2].parse().unwrap_or(0.0);
                    let z: f64 = parts[3].parse().unwrap_or(0.0);
                    frame.push([x, y, z]);
                }
                idx += 1;
            }
            if !frame.is_empty() {
                frames.push(frame);
            }
        }
    }
    let n_frames = frames.len();
    let n_atoms = frames.first().map(|f| f.len()).unwrap_or(0);
    Ok(ChunkedTrajectory {
        n_frames,
        n_atoms,
        chunk_size: n_frames.max(1),
        frames,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Parallel HDF5-like layout
// ─────────────────────────────────────────────────────────────────────────────

/// Data-type descriptor for a dataset column.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DsetType {
    /// 64-bit float.
    Float64,
    /// 32-bit float.
    Float32,
    /// 64-bit unsigned integer.
    Uint64,
    /// 32-bit unsigned integer.
    Uint32,
}

impl fmt::Display for DsetType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            DsetType::Float64 => "float64",
            DsetType::Float32 => "float32",
            DsetType::Uint64 => "uint64",
            DsetType::Uint32 => "uint32",
        };
        write!(f, "{s}")
    }
}

impl DsetType {
    /// Size in bytes of a single element.
    pub fn byte_size(&self) -> usize {
        match self {
            DsetType::Float64 | DsetType::Uint64 => 8,
            DsetType::Float32 | DsetType::Uint32 => 4,
        }
    }
}

/// Header for a single named dataset within a parallel file.
#[derive(Debug, Clone)]
pub struct DatasetHeader {
    /// Human-readable dataset name.
    pub name: String,
    /// Element type.
    pub dtype: DsetType,
    /// Number of components per element (e.g. 3 for a 3-D vector).
    pub n_components: usize,
    /// Total number of elements in the global dataset.
    pub global_n: usize,
    /// Byte offset of this dataset's data block within the file body.
    pub byte_offset: usize,
}

impl DatasetHeader {
    /// Construct a new dataset header.
    pub fn new(
        name: &str,
        dtype: DsetType,
        n_components: usize,
        global_n: usize,
        byte_offset: usize,
    ) -> Self {
        Self {
            name: name.to_owned(),
            dtype,
            n_components,
            global_n,
            byte_offset,
        }
    }

    /// Total byte size of this dataset's data block.
    pub fn data_size_bytes(&self) -> usize {
        self.global_n * self.n_components * self.dtype.byte_size()
    }

    /// Serialize to a text header line.
    pub fn to_header_line(&self) -> String {
        format!(
            "DATASET name={} dtype={} n_comp={} global_n={} offset={}\n",
            self.name, self.dtype, self.n_components, self.global_n, self.byte_offset
        )
    }
}

/// A lightweight parallel-HDF5-like file layout description.
#[derive(Debug, Clone, Default)]
pub struct ParallelFileLayout {
    /// All dataset headers in declaration order.
    pub datasets: Vec<DatasetHeader>,
}

impl ParallelFileLayout {
    /// Create an empty layout.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register a dataset and compute its byte offset automatically.
    pub fn add_dataset(
        &mut self,
        name: &str,
        dtype: DsetType,
        n_components: usize,
        global_n: usize,
    ) {
        let offset = self
            .datasets
            .last()
            .map(|d| d.byte_offset + d.data_size_bytes())
            .unwrap_or(0);
        self.datasets.push(DatasetHeader::new(
            name,
            dtype,
            n_components,
            global_n,
            offset,
        ));
    }

    /// Total byte size of all dataset data blocks.
    pub fn total_data_bytes(&self) -> usize {
        self.datasets
            .last()
            .map(|d| d.byte_offset + d.data_size_bytes())
            .unwrap_or(0)
    }

    /// Serialize all headers to a string.
    pub fn header_block(&self) -> String {
        self.datasets.iter().map(|d| d.to_header_line()).collect()
    }

    /// Look up a dataset by name.
    pub fn find(&self, name: &str) -> Option<&DatasetHeader> {
        self.datasets.iter().find(|d| d.name == name)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Delta encoding
// ─────────────────────────────────────────────────────────────────────────────

/// Compress trajectory frames using delta encoding.
///
/// The first frame is stored as-is; subsequent frames store the difference
/// from the previous frame.
pub fn compress_trajectory_delta(frames: &[Vec<[f64; 3]>]) -> Vec<Vec<[f64; 3]>> {
    if frames.is_empty() {
        return Vec::new();
    }
    let mut result = Vec::with_capacity(frames.len());
    result.push(frames[0].clone());
    for i in 1..frames.len() {
        let prev = &frames[i - 1];
        let curr = &frames[i];
        let n = curr.len().min(prev.len());
        let delta: Vec<[f64; 3]> = curr
            .iter()
            .enumerate()
            .map(|(j, &c)| {
                if j < n {
                    [c[0] - prev[j][0], c[1] - prev[j][1], c[2] - prev[j][2]]
                } else {
                    c
                }
            })
            .collect();
        result.push(delta);
    }
    result
}

/// Reconstruct trajectory frames from delta-encoded data.
pub fn decompress_trajectory_delta(deltas: &[Vec<[f64; 3]>]) -> Vec<Vec<[f64; 3]>> {
    if deltas.is_empty() {
        return Vec::new();
    }
    let mut result = Vec::with_capacity(deltas.len());
    result.push(deltas[0].clone());
    for i in 1..deltas.len() {
        let prev = &result[i - 1];
        let delta = &deltas[i];
        let frame: Vec<[f64; 3]> = delta
            .iter()
            .enumerate()
            .map(|(j, &d)| {
                if j < prev.len() {
                    [prev[j][0] + d[0], prev[j][1] + d[1], prev[j][2] + d[2]]
                } else {
                    d
                }
            })
            .collect();
        result.push(frame);
    }
    result
}

// ─────────────────────────────────────────────────────────────────────────────
// I/O statistics
// ─────────────────────────────────────────────────────────────────────────────

/// Return a formatted I/O statistics summary string.
///
/// * `bytes_read` — total bytes read.
/// * `bytes_written` — total bytes written.
/// * `frames` — number of trajectory frames processed.
pub fn io_stats(bytes_read: usize, bytes_written: usize, frames: usize) -> String {
    format!("IO Stats: read={bytes_read}B written={bytes_written}B frames={frames}")
}

/// Aggregated I/O performance counters.
#[derive(Debug, Clone, Default)]
pub struct IoCounters {
    /// Total bytes read across all operations.
    pub bytes_read: u64,
    /// Total bytes written across all operations.
    pub bytes_written: u64,
    /// Number of read operations.
    pub read_ops: u64,
    /// Number of write operations.
    pub write_ops: u64,
    /// Cumulative wall-clock time spent in read calls (seconds).
    pub read_time_s: f64,
    /// Cumulative wall-clock time spent in write calls (seconds).
    pub write_time_s: f64,
}

impl IoCounters {
    /// Create zeroed counters.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a completed read operation.
    pub fn record_read(&mut self, bytes: u64, elapsed_s: f64) {
        self.bytes_read += bytes;
        self.read_ops += 1;
        self.read_time_s += elapsed_s;
    }

    /// Record a completed write operation.
    pub fn record_write(&mut self, bytes: u64, elapsed_s: f64) {
        self.bytes_written += bytes;
        self.write_ops += 1;
        self.write_time_s += elapsed_s;
    }

    /// Average read bandwidth in MB/s.
    pub fn avg_read_bw_mb(&self) -> f64 {
        if self.read_time_s < 1e-15 {
            0.0
        } else {
            self.bytes_read as f64 / self.read_time_s / 1e6
        }
    }

    /// Average write bandwidth in MB/s.
    pub fn avg_write_bw_mb(&self) -> f64 {
        if self.write_time_s < 1e-15 {
            0.0
        } else {
            self.bytes_written as f64 / self.write_time_s / 1e6
        }
    }

    /// Reset all counters to zero.
    pub fn reset(&mut self) {
        *self = Self::default();
    }

    /// Merge another counter set into this one.
    pub fn merge(&mut self, other: &IoCounters) {
        self.bytes_read += other.bytes_read;
        self.bytes_written += other.bytes_written;
        self.read_ops += other.read_ops;
        self.write_ops += other.write_ops;
        self.read_time_s += other.read_time_s;
        self.write_time_s += other.write_time_s;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── DomainBox ──────────────────────────────────────────────────────────────

    #[test]
    fn test_domain_box_volume() {
        let b = DomainBox::new([0.0, 0.0, 0.0], [2.0, 3.0, 4.0]);
        assert!((b.volume() - 24.0).abs() < 1e-10);
    }

    #[test]
    fn test_domain_box_volume_zero() {
        let b = DomainBox::new([1.0, 1.0, 1.0], [1.0, 2.0, 2.0]);
        assert_eq!(b.volume(), 0.0);
    }

    #[test]
    fn test_domain_box_centre() {
        let b = DomainBox::new([0.0, 0.0, 0.0], [2.0, 4.0, 6.0]);
        let c = b.centre();
        assert!((c[0] - 1.0).abs() < 1e-10);
        assert!((c[1] - 2.0).abs() < 1e-10);
        assert!((c[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_domain_box_contains() {
        let b = DomainBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        assert!(b.contains([0.5, 0.5, 0.5]));
        assert!(!b.contains([1.5, 0.5, 0.5]));
        assert!(!b.contains([0.0, 0.5, 0.5])); // boundary is excluded
    }

    #[test]
    fn test_domain_box_overlaps() {
        let a = DomainBox::new([0.0, 0.0, 0.0], [1.0, 1.0, 1.0]);
        let b = DomainBox::new([0.5, 0.5, 0.5], [1.5, 1.5, 1.5]);
        assert!(a.overlaps(&b));
        let c = DomainBox::new([2.0, 2.0, 2.0], [3.0, 3.0, 3.0]);
        assert!(!a.overlaps(&c));
    }

    // ── RankDomain ────────────────────────────────────────────────────────────

    #[test]
    fn test_rank_domain_n_owned_cells() {
        let rd = RankDomain::new(
            0,
            4,
            DomainBox::new([0.0; 3], [1.0; 3]),
            [1; 6],
            [4, 4, 4],
            [0; 3],
        );
        assert_eq!(rd.n_owned_cells(), 64);
    }

    #[test]
    fn test_rank_domain_n_total_cells() {
        let rd = RankDomain::new(
            0,
            1,
            DomainBox::new([0.0; 3], [1.0; 3]),
            [1, 1, 1, 1, 1, 1],
            [4, 4, 4],
            [0; 3],
        );
        // ext: (4+1+1)^3 = 216
        assert_eq!(rd.n_total_cells(), 216);
    }

    // ── GhostLayer / GhostCellMeta ────────────────────────────────────────────

    #[test]
    fn test_ghost_layer_n_cells() {
        let gl = GhostLayer {
            face: Face::XMinus,
            neighbour_rank: 1,
            depth: 2,
            tangential_cells: [8, 8],
        };
        assert_eq!(gl.n_cells(), 128);
    }

    #[test]
    fn test_ghost_cell_meta_total() {
        let mut meta = GhostCellMeta::new();
        meta.add_layer(GhostLayer {
            face: Face::XPlus,
            neighbour_rank: 1,
            depth: 1,
            tangential_cells: [4, 4],
        });
        meta.add_layer(GhostLayer {
            face: Face::XMinus,
            neighbour_rank: 2,
            depth: 1,
            tangential_cells: [4, 4],
        });
        assert_eq!(meta.total_ghost_cells(), 32);
    }

    #[test]
    fn test_ghost_cell_meta_neighbour_ranks() {
        let mut meta = GhostCellMeta::new();
        meta.add_layer(GhostLayer {
            face: Face::XPlus,
            neighbour_rank: 3,
            depth: 1,
            tangential_cells: [2, 2],
        });
        meta.add_layer(GhostLayer {
            face: Face::YMinus,
            neighbour_rank: 1,
            depth: 1,
            tangential_cells: [2, 2],
        });
        meta.add_layer(GhostLayer {
            face: Face::ZPlus,
            neighbour_rank: 3,
            depth: 1,
            tangential_cells: [2, 2],
        });
        let ranks = meta.neighbour_ranks();
        assert_eq!(ranks, vec![1, 3]);
    }

    #[test]
    fn test_face_display() {
        assert_eq!(format!("{}", Face::XMinus), "XMinus");
        assert_eq!(format!("{}", Face::ZPlus), "ZPlus");
    }

    // ── ChunkSpec ─────────────────────────────────────────────────────────────

    #[test]
    fn test_chunk_spec_size() {
        let c = ChunkSpec {
            start: 5,
            end: 10,
            chunk_idx: 0,
        };
        assert_eq!(c.size(), 5);
    }

    #[test]
    fn test_chunk_spec_size_zero() {
        let c = ChunkSpec {
            start: 10,
            end: 5,
            chunk_idx: 0,
        };
        assert_eq!(c.size(), 0);
    }

    #[test]
    fn test_split_into_chunks_exact() {
        let chunks = split_into_chunks(12, 3);
        assert_eq!(chunks.len(), 3);
        assert_eq!(chunks[0].size(), 4);
        assert_eq!(chunks[1].size(), 4);
        assert_eq!(chunks[2].size(), 4);
    }

    #[test]
    fn test_split_into_chunks_remainder() {
        let chunks = split_into_chunks(10, 3);
        assert_eq!(chunks.len(), 3);
        let total: usize = chunks.iter().map(|c| c.size()).sum();
        assert_eq!(total, 10);
    }

    #[test]
    fn test_split_into_chunks_zero_total() {
        assert!(split_into_chunks(0, 4).is_empty());
    }

    #[test]
    fn test_split_into_chunks_zero_n() {
        assert!(split_into_chunks(10, 0).is_empty());
    }

    #[test]
    fn test_split_cover_all_elements() {
        for total in [1, 5, 100, 1000] {
            for n_chunks in [1, 3, 7, 10] {
                let chunks = split_into_chunks(total, n_chunks);
                let covered: usize = chunks.iter().map(|c| c.size()).sum();
                assert_eq!(covered, total, "total={total} n_chunks={n_chunks}");
            }
        }
    }

    // ── ParallelWriter ────────────────────────────────────────────────────────

    #[test]
    fn test_parallel_writer_progress() {
        let mut pw = ParallelWriter::new("test.bin", 4);
        pw.mark_chunk_written();
        pw.mark_chunk_written();
        assert!((pw.progress() - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_parallel_writer_complete() {
        let mut pw = ParallelWriter::new("test.bin", 2);
        pw.mark_chunk_written();
        pw.mark_chunk_written();
        assert!(pw.is_complete());
    }

    #[test]
    fn test_parallel_writer_no_overflow() {
        let mut pw = ParallelWriter::new("test.bin", 1);
        pw.mark_chunk_written();
        pw.mark_chunk_written();
        assert_eq!(pw.chunks_written, 1);
    }

    #[test]
    fn test_parallel_writer_zero_chunks() {
        let pw = ParallelWriter::new("test.bin", 0);
        assert!((pw.progress() - 1.0).abs() < 1e-10);
        assert!(pw.is_complete());
    }

    #[test]
    fn test_parallel_writer_reset() {
        let mut pw = ParallelWriter::new("f", 3);
        pw.mark_chunk_written();
        pw.mark_chunk_written();
        pw.reset();
        assert_eq!(pw.chunks_written, 0);
        assert!(!pw.is_complete());
    }

    // ── PVTU / PVD XML ────────────────────────────────────────────────────────

    #[test]
    fn test_pvtu_xml_contains_pieces() {
        let pieces = vec![
            PvtuPiece {
                file: "rank0.vtu".into(),
                rank: 0,
            },
            PvtuPiece {
                file: "rank1.vtu".into(),
                rank: 1,
            },
        ];
        let xml = write_pvtu_xml(&pieces, &["velocity"], &["pressure"]);
        assert!(xml.contains("rank0.vtu"));
        assert!(xml.contains("rank1.vtu"));
        assert!(xml.contains("velocity"));
        assert!(xml.contains("pressure"));
        assert!(xml.contains("PUnstructuredGrid"));
    }

    #[test]
    fn test_pvtu_xml_no_arrays() {
        let pieces = vec![PvtuPiece {
            file: "r0.vtu".into(),
            rank: 0,
        }];
        let xml = write_pvtu_xml(&pieces, &[], &[]);
        assert!(xml.contains("PUnstructuredGrid"));
        assert!(!xml.contains("PPointData"));
        assert!(!xml.contains("PCellData"));
    }

    #[test]
    fn test_pvd_xml_entries() {
        let entries = vec![(0.0, "t0.pvtu"), (0.1, "t1.pvtu"), (0.2, "t2.pvtu")];
        let xml = write_pvd_xml(&entries);
        assert!(xml.contains("t0.pvtu"));
        assert!(xml.contains("t2.pvtu"));
        assert!(xml.contains("Collection"));
    }

    // ── PartitionMap ──────────────────────────────────────────────────────────

    #[test]
    fn test_partition_map_round_robin() {
        let pm = PartitionMap::round_robin(9, 3);
        assert_eq!(pm.cells_on_rank(0), 3);
        assert_eq!(pm.cells_on_rank(1), 3);
        assert_eq!(pm.cells_on_rank(2), 3);
    }

    #[test]
    fn test_partition_map_imbalance_perfect() {
        let pm = PartitionMap::round_robin(6, 3);
        assert!((pm.imbalance() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_partition_map_imbalance_skewed() {
        // All cells on rank 0 → high imbalance
        let pm = PartitionMap::new(vec![0; 9], 3);
        assert!(pm.imbalance() > 1.5);
    }

    #[test]
    fn test_partition_map_to_text_roundtrip_count() {
        let pm = PartitionMap::round_robin(4, 2);
        let text = pm.to_text();
        // Should have n_ranks line + 4 cell lines
        assert_eq!(text.lines().count(), 5);
    }

    // ── DistributedArray ──────────────────────────────────────────────────────

    #[test]
    fn test_distributed_array_serialize_deserialize() {
        let arr = DistributedArray::new(0, 6, 0, vec![1.0, 2.0, 3.0]);
        let text = arr.serialize_text();
        let recovered = DistributedArray::deserialize_text(&text).unwrap();
        assert_eq!(recovered.rank, 0);
        assert_eq!(recovered.offset, 0);
        assert_eq!(recovered.global_len, 6);
        assert!((recovered.data[0] - 1.0).abs() < 1e-10);
        assert!((recovered.data[2] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_distributed_array_assemble() {
        let shard0 = DistributedArray::new(0, 6, 0, vec![1.0, 2.0, 3.0]);
        let shard1 = DistributedArray::new(1, 6, 3, vec![4.0, 5.0, 6.0]);
        let global = DistributedArray::assemble(&[shard0, shard1]).unwrap();
        assert_eq!(global.len(), 6);
        for (i, &v) in global.iter().enumerate() {
            assert!((v - (i + 1) as f64).abs() < 1e-10);
        }
    }

    #[test]
    fn test_distributed_array_assemble_empty() {
        let result = DistributedArray::assemble(&[]).unwrap();
        assert!(result.is_empty());
    }

    // ── CheckpointHeader ──────────────────────────────────────────────────────

    #[test]
    fn test_checkpoint_header_roundtrip() {
        let h = CheckpointHeader::new(42, 3.125, 4, 1000, "sim_v1");
        let line = h.to_line();
        let h2 = CheckpointHeader::from_line(&line).unwrap();
        assert_eq!(h2.step, 42);
        assert!((h2.time - 3.125).abs() < 1e-6);
        assert_eq!(h2.n_ranks, 4);
        assert_eq!(h2.n_items, 1000);
        assert_eq!(h2.tag, "sim_v1");
    }

    #[test]
    fn test_checkpoint_header_bad_line() {
        assert!(CheckpointHeader::from_line("garbage").is_none());
    }

    // ── write_checkpoint / read_checkpoint ────────────────────────────────────

    #[test]
    fn test_checkpoint_write_read_roundtrip() {
        let path = std::env::temp_dir().join("test_oxiphysics_ckpt.txt");
        let header = CheckpointHeader::new(10, 1.5, 2, 3, "test");
        let positions = vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        write_checkpoint(path.to_str().unwrap_or(""), &header, &positions).unwrap();
        let (h2, pos2) = read_checkpoint(path.to_str().unwrap_or("")).unwrap();
        assert_eq!(h2.step, 10);
        assert_eq!(pos2.len(), 3);
        assert!((pos2[0][0] - 1.0).abs() < 1e-6);
        assert!((pos2[2][2] - 9.0).abs() < 1e-6);
        let _ = fs::remove_file(&path);
    }

    // ── merge_chunk_files / merge_rank_files ──────────────────────────────────

    #[test]
    fn test_merge_chunk_files() {
        let tmpdir = std::env::temp_dir();
        let base = tmpdir
            .join("test_merge_par_chunks")
            .to_str()
            .unwrap_or("")
            .to_string();
        let output = tmpdir.join("test_merge_par_output.bin");
        for i in 0..3usize {
            fs::write(format!("{base}.chunk{i}"), format!("chunk{i}\n")).unwrap();
        }
        merge_chunk_files(&base, 3, output.to_str().unwrap_or("")).unwrap();
        let merged = fs::read_to_string(&output).unwrap();
        assert!(merged.contains("chunk0"));
        assert!(merged.contains("chunk2"));
        for i in 0..3usize {
            let _ = fs::remove_file(format!("{base}.chunk{i}"));
        }
        let _ = fs::remove_file(&output);
    }

    #[test]
    fn test_merge_rank_files() {
        let tmpdir = std::env::temp_dir();
        let base = tmpdir
            .join("test_merge_par_ranks")
            .to_str()
            .unwrap_or("")
            .to_string();
        let output = tmpdir.join("test_merge_par_rank_output.bin");
        for r in 0..2usize {
            fs::write(format!("{base}.rank{r}"), format!("rank{r}\n")).unwrap();
        }
        merge_rank_files(&base, 2, output.to_str().unwrap_or("")).unwrap();
        let merged = fs::read_to_string(&output).unwrap();
        assert!(merged.contains("rank0"));
        assert!(merged.contains("rank1"));
        for r in 0..2usize {
            let _ = fs::remove_file(format!("{base}.rank{r}"));
        }
        let _ = fs::remove_file(&output);
    }

    // ── estimate_io_bandwidth ─────────────────────────────────────────────────

    #[test]
    fn test_estimate_io_bandwidth_normal() {
        let bw = estimate_io_bandwidth(1_000_000, 1.0);
        assert!((bw - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_estimate_io_bandwidth_zero_time() {
        assert_eq!(estimate_io_bandwidth(1000, 0.0), 0.0);
    }

    #[test]
    fn test_bandwidth_large() {
        let bw = estimate_io_bandwidth(1_000_000_000, 0.5);
        assert!((bw - 2000.0).abs() < 1.0);
    }

    // ── ChunkedTrajectory ─────────────────────────────────────────────────────

    #[test]
    fn test_chunked_trajectory_n_chunks() {
        let mut traj = ChunkedTrajectory::new(2, 3);
        for _ in 0..7 {
            traj.push_frame(vec![[0.0; 3]; 2]);
        }
        assert_eq!(traj.n_chunks(), 3);
    }

    #[test]
    fn test_write_read_chunked_xyz_roundtrip() {
        let base = std::env::temp_dir()
            .join("test_chunked_traj_par2")
            .to_str()
            .unwrap_or("")
            .to_string();
        let mut traj = ChunkedTrajectory::new(2, 2);
        traj.push_frame(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        traj.push_frame(vec![[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]);
        traj.push_frame(vec![[9.0, 8.0, 7.0], [6.0, 5.0, 4.0]]);
        write_chunked_xyz(&traj, &base).unwrap();
        let loaded = read_chunked_xyz(&base, traj.n_chunks()).unwrap();
        assert_eq!(loaded.n_frames, 3);
        assert!((loaded.frames[0][0][0] - 1.0).abs() < 1e-4);
        for i in 0..traj.n_chunks() {
            let _ = fs::remove_file(format!("{base}.chunk{i}.xyz"));
        }
    }

    // ── Delta encoding ────────────────────────────────────────────────────────

    #[test]
    fn test_compress_decompress_delta_roundtrip() {
        let frames = vec![
            vec![[1.0_f64, 2.0, 3.0], [4.0, 5.0, 6.0]],
            vec![[1.1, 2.1, 3.1], [4.2, 5.2, 6.2]],
            vec![[1.5, 2.5, 3.5], [5.0, 5.5, 6.5]],
        ];
        let deltas = compress_trajectory_delta(&frames);
        let recovered = decompress_trajectory_delta(&deltas);
        assert_eq!(recovered.len(), frames.len());
        for (f, r) in frames.iter().zip(recovered.iter()) {
            for (fa, ra) in f.iter().zip(r.iter()) {
                assert!((fa[0] - ra[0]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_compress_delta_single_frame() {
        let frames = vec![vec![[1.0_f64, 2.0, 3.0]]];
        let d = compress_trajectory_delta(&frames);
        assert_eq!(d.len(), 1);
        assert!((d[0][0][0] - 1.0).abs() < 1e-10);
    }

    // ── ParallelFileLayout ────────────────────────────────────────────────────

    #[test]
    fn test_parallel_file_layout_offsets() {
        let mut layout = ParallelFileLayout::new();
        layout.add_dataset("positions", DsetType::Float64, 3, 100);
        layout.add_dataset("velocity", DsetType::Float64, 3, 100);
        let pos = layout.find("positions").unwrap();
        assert_eq!(pos.byte_offset, 0);
        let vel = layout.find("velocity").unwrap();
        // 100 * 3 * 8 = 2400 bytes
        assert_eq!(vel.byte_offset, 2400);
    }

    #[test]
    fn test_parallel_file_layout_total_bytes() {
        let mut layout = ParallelFileLayout::new();
        layout.add_dataset("x", DsetType::Float32, 1, 50);
        layout.add_dataset("y", DsetType::Uint32, 1, 50);
        // 50*4 + 50*4 = 400
        assert_eq!(layout.total_data_bytes(), 400);
    }

    #[test]
    fn test_dataset_type_byte_sizes() {
        assert_eq!(DsetType::Float64.byte_size(), 8);
        assert_eq!(DsetType::Float32.byte_size(), 4);
        assert_eq!(DsetType::Uint64.byte_size(), 8);
        assert_eq!(DsetType::Uint32.byte_size(), 4);
    }

    #[test]
    fn test_dataset_header_serialise() {
        let h = DatasetHeader::new("temp", DsetType::Float32, 1, 256, 0);
        let line = h.to_header_line();
        assert!(line.contains("DATASET"));
        assert!(line.contains("temp"));
        assert!(line.contains("float32"));
    }

    // ── IoCounters ────────────────────────────────────────────────────────────

    #[test]
    fn test_io_counters_read_bw() {
        let mut c = IoCounters::new();
        c.record_read(10_000_000, 0.01);
        let bw = c.avg_read_bw_mb();
        assert!((bw - 1000.0).abs() < 1.0);
    }

    #[test]
    fn test_io_counters_merge() {
        let mut a = IoCounters::new();
        a.record_write(1000, 0.001);
        let mut b = IoCounters::new();
        b.record_write(2000, 0.002);
        a.merge(&b);
        assert_eq!(a.bytes_written, 3000);
        assert_eq!(a.write_ops, 2);
    }

    #[test]
    fn test_io_counters_reset() {
        let mut c = IoCounters::new();
        c.record_read(999, 1.0);
        c.reset();
        assert_eq!(c.bytes_read, 0);
        assert_eq!(c.read_ops, 0);
    }

    #[test]
    fn test_io_stats_format() {
        let s = io_stats(1024, 2048, 10);
        assert!(s.contains("1024"));
        assert!(s.contains("2048"));
        assert!(s.contains("10"));
    }
}
