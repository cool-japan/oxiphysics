// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Rolling checkpoint management, incremental (delta) extraction and simple
//! checksums for restart blobs.

use std::path::Path;

use super::{RestartData, RestartFormat, RestartReader, RestartWriter};

// ── CheckpointManager ─────────────────────────────────────────────────────────

/// Manages a rolling set of checkpoint files in a directory.
#[derive(Debug, Clone)]
pub struct CheckpointManager {
    /// Base directory where checkpoint files are stored.
    pub base_dir: String,
    /// Maximum number of checkpoints to retain.
    pub max_checkpoints: usize,
    /// In-memory list of `(step, file_path)` known checkpoints, oldest first.
    checkpoints: Vec<(u64, String)>,
}

impl CheckpointManager {
    /// Create a new manager rooted at `base_dir`, keeping at most `max_checkpoints` files.
    ///
    /// `max_checkpoints` is raised to at least 1 so the checkpoint that was
    /// just written is never pruned away. Nothing is created on disk until the
    /// first [`save_checkpoint`](Self::save_checkpoint).
    pub fn new(base_dir: &str, max_checkpoints: usize) -> Self {
        Self {
            base_dir: base_dir.to_string(),
            max_checkpoints: max_checkpoints.max(1),
            checkpoints: Vec::new(),
        }
    }

    /// Serialize `data` and write it as `checkpoint_<step>.bin` under `base_dir`,
    /// with `<step>` zero-padded to ten digits.
    ///
    /// Returns the full path of the newly written file. The order of
    /// operations protects the checkpoints that already exist: the directory is
    /// created and the new file written (and verified to have been written)
    /// first; only then is it recorded and the oldest checkpoints beyond
    /// `max_checkpoints` deleted. If anything before that point fails, nothing
    /// is recorded and nothing is deleted. Saving a step that is already
    /// recorded replaces its file and keeps a single entry.
    ///
    /// # Errors
    ///
    /// Returns an error string if the directory cannot be created, the data is
    /// inconsistent ([`RestartData::validate`]) or the file cannot be written,
    /// or — after the new checkpoint is safely on disk and recorded — if an
    /// old checkpoint cannot be deleted (see
    /// [`prune_old_checkpoints`](Self::prune_old_checkpoints)).
    pub fn save_checkpoint(&mut self, data: &RestartData, step: u64) -> Result<String, String> {
        let path = Path::new(&self.base_dir).join(format!("checkpoint_{step:010}.bin"));
        let filename = path.to_string_lossy().into_owned();
        std::fs::create_dir_all(&self.base_dir).map_err(|e| {
            format!(
                "failed to create checkpoint directory '{}': {e}",
                self.base_dir
            )
        })?;
        RestartWriter::new(&filename, RestartFormat::Binary).write(data)?;
        self.checkpoints.retain(|(_, known)| *known != filename);
        self.checkpoints.push((step, filename.clone()));
        self.prune_old_checkpoints().map_err(|e| {
            format!("checkpoint '{filename}' was written but pruning old checkpoints failed: {e}")
        })?;
        Ok(filename)
    }

    /// Load and return the most recent checkpoint, or `Ok(None)` if none exist.
    ///
    /// # Errors
    ///
    /// Returns an error string if the most recent checkpoint file cannot be
    /// read or parsed.
    pub fn load_latest(&self) -> Result<Option<RestartData>, String> {
        let Some((_, path)) = self.checkpoints.last() else {
            return Ok(None);
        };
        RestartReader::new(path).read().map(Some)
    }

    /// Return all known checkpoints as `(step, path)` pairs, oldest first.
    pub fn list_checkpoints(&self) -> Vec<(u64, String)> {
        self.checkpoints.clone()
    }

    /// Remove the oldest checkpoints so that at most `max_checkpoints` remain.
    ///
    /// The newest checkpoint is always kept, even if the public
    /// `max_checkpoints` field was set to 0 after construction.
    ///
    /// A checkpoint is dropped from the list only after its file has been
    /// deleted (a file that is already gone counts as deleted). If a file
    /// cannot be deleted the error is returned and that checkpoint stays
    /// listed, so a later call retries it.
    pub fn prune_old_checkpoints(&mut self) -> Result<(), String> {
        let keep = self.max_checkpoints.max(1);
        while self.checkpoints.len() > keep {
            let Some((_, path)) = self.checkpoints.first() else {
                break;
            };
            match std::fs::remove_file(path) {
                Ok(()) => {}
                Err(e) if e.kind() == std::io::ErrorKind::NotFound => {}
                Err(e) => {
                    return Err(format!("failed to remove checkpoint '{path}': {e}"));
                }
            }
            self.checkpoints.remove(0);
        }
        Ok(())
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
    /// Metadata and the box are copied from `full`; the per-particle arrays are
    /// filtered to the `changed` indices (indices beyond the end of an array
    /// are skipped, and an array that is absent in `full` stays absent in the
    /// delta). The extra scalar / vector arrays are not carried over.
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
