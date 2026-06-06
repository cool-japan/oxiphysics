// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Trajectory storage, checkpoints, virtual datasets.

use std::collections::HashMap;

use super::file::Hdf5File;
use super::types::{AttrValue, ExternalRef, Hdf5Dtype, Hdf5Error, Hdf5Result, Hyperslab};

// ---------------------------------------------------------------------------
// Trajectory storage
// ---------------------------------------------------------------------------

/// Multi-frame trajectory container: `(n_atoms x n_frames x 3)`.
///
/// Follows the HDF5MD convention used by MDAnalysis / h5md.
#[derive(Debug, Clone)]
pub struct TrajectoryStore {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Number of time frames stored.
    pub n_frames: usize,
    /// Flat storage: positions\[frame * n_atoms * 3 + atom * 3 + xyz\].
    pub positions: Vec<f64>,
    /// Flat storage: velocities (same layout as positions, optional).
    pub velocities: Vec<f64>,
    /// Flat storage: forces (same layout, optional).
    pub forces: Vec<f64>,
    /// Time value for each frame (in picoseconds).
    pub times: Vec<f64>,
    /// Per-frame box vectors \[frame * 9\]: row-major 3x3.
    pub box_vectors: Vec<f64>,
    /// Whether velocities are stored.
    pub has_velocities: bool,
    /// Whether forces are stored.
    pub has_forces: bool,
}

impl TrajectoryStore {
    /// Create an empty trajectory for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            n_atoms,
            n_frames: 0,
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            times: Vec::new(),
            box_vectors: Vec::new(),
            has_velocities: false,
            has_forces: false,
        }
    }

    /// Append a single frame from a slice of positions `[n_atoms * 3]`.
    pub fn append_frame(&mut self, time: f64, positions: &[f64], box_vecs: &[f64; 9]) {
        assert_eq!(
            positions.len(),
            self.n_atoms * 3,
            "append_frame: positions must have n_atoms*3 elements"
        );
        self.times.push(time);
        self.positions.extend_from_slice(positions);
        self.box_vectors.extend_from_slice(box_vecs);
        self.n_frames += 1;
    }

    /// Append velocities for the last frame.
    ///
    /// Must be called immediately after `append_frame` and before the next one.
    pub fn append_velocities(&mut self, velocities: &[f64]) {
        assert_eq!(velocities.len(), self.n_atoms * 3);
        self.velocities.extend_from_slice(velocities);
        self.has_velocities = true;
    }

    /// Append forces for the last frame.
    pub fn append_forces(&mut self, forces: &[f64]) {
        assert_eq!(forces.len(), self.n_atoms * 3);
        self.forces.extend_from_slice(forces);
        self.has_forces = true;
    }

    /// Read back the positions for frame `frame_id`.
    pub fn read_positions(&self, frame_id: usize) -> Hdf5Result<Vec<[f64; 3]>> {
        if frame_id >= self.n_frames {
            return Err(Hdf5Error::NotFound(format!("frame {frame_id}")));
        }
        let base = frame_id * self.n_atoms * 3;
        let out: Vec<[f64; 3]> = (0..self.n_atoms)
            .map(|a| {
                let i = base + a * 3;
                [
                    self.positions[i],
                    self.positions[i + 1],
                    self.positions[i + 2],
                ]
            })
            .collect();
        Ok(out)
    }

    /// Read back the velocities for frame `frame_id`.
    pub fn read_velocities(&self, frame_id: usize) -> Hdf5Result<Vec<[f64; 3]>> {
        if !self.has_velocities {
            return Err(Hdf5Error::Generic(
                "trajectory has no velocities".to_string(),
            ));
        }
        if frame_id >= self.n_frames {
            return Err(Hdf5Error::NotFound(format!("frame {frame_id}")));
        }
        let base = frame_id * self.n_atoms * 3;
        let out: Vec<[f64; 3]> = (0..self.n_atoms)
            .map(|a| {
                let i = base + a * 3;
                [
                    self.velocities[i],
                    self.velocities[i + 1],
                    self.velocities[i + 2],
                ]
            })
            .collect();
        Ok(out)
    }

    /// Return the total trajectory duration in picoseconds.
    pub fn total_time(&self) -> f64 {
        self.times.last().copied().unwrap_or(0.0) - self.times.first().copied().unwrap_or(0.0)
    }

    /// Flush the trajectory into an HDF5 file under `/trajectory`.
    ///
    /// Creates two datasets: `positions` of shape `[n_frames, n_atoms, 3]`
    /// and `time` of shape `[n_frames]`.
    pub fn flush_to_file(&self, file: &mut Hdf5File) -> Hdf5Result<()> {
        file.create_group("trajectory")?;
        file.create_dataset(
            "trajectory",
            "positions",
            vec![self.n_frames, self.n_atoms, 3],
            Hdf5Dtype::Float64,
        )?;
        let ds = file.open_dataset_mut("trajectory", "positions")?;
        ds.write_f64(&self.positions)?;

        // time dataset
        let group = file.open_group_mut("trajectory")?;
        group.create_dataset("time", vec![self.n_frames], Hdf5Dtype::Float64)?;
        let tds = group.open_dataset_mut("time")?;
        tds.write_f64(&self.times)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Checkpoint / restart
// ---------------------------------------------------------------------------

/// Simulation checkpoint: serialises/deserialises the key state vectors.
#[derive(Debug, Clone)]
pub struct Checkpoint {
    /// Checkpoint label (e.g. `"step_1000"`).
    pub label: String,
    /// Atom positions `[n_atoms * 3]`.
    pub positions: Vec<f64>,
    /// Atom velocities `[n_atoms * 3]`.
    pub velocities: Vec<f64>,
    /// Simulation step number.
    pub step: u64,
    /// Simulation time (ps).
    pub time: f64,
    /// Box vectors (9 elements, row-major 3x3).
    pub box_vectors: [f64; 9],
    /// Optional additional scalar fields.
    pub scalars: HashMap<String, f64>,
}

impl Checkpoint {
    /// Create a new checkpoint.
    pub fn new(label: &str, step: u64, time: f64) -> Self {
        Self {
            label: label.to_string(),
            positions: Vec::new(),
            velocities: Vec::new(),
            step,
            time,
            box_vectors: [0.0; 9],
            scalars: HashMap::new(),
        }
    }

    /// Write this checkpoint into an HDF5 file under `/checkpoints/`label`.
    pub fn write_to_file(&self, file: &mut Hdf5File) -> Hdf5Result<()> {
        let base = format!("checkpoints/{}", self.label);
        file.create_group("checkpoints")?;
        file.create_group(&base)?;

        let n_atoms = self.positions.len() / 3;

        // positions
        file.create_dataset(&base, "positions", vec![n_atoms, 3], Hdf5Dtype::Float64)?;
        file.open_dataset_mut(&base, "positions")?
            .write_f64(&self.positions)?;

        // velocities
        file.create_dataset(&base, "velocities", vec![n_atoms, 3], Hdf5Dtype::Float64)?;
        file.open_dataset_mut(&base, "velocities")?
            .write_f64(&self.velocities)?;

        // step / time attributes
        file.set_dataset_attr(
            &base,
            "positions",
            "step",
            AttrValue::Int32(self.step as i32),
        )?;
        file.set_dataset_attr(&base, "positions", "time", AttrValue::Float64(self.time))?;

        // box vectors
        file.create_dataset(&base, "box", vec![3, 3], Hdf5Dtype::Float64)?;
        file.open_dataset_mut(&base, "box")?
            .write_f64(&self.box_vectors)?;

        Ok(())
    }

    /// Restore a checkpoint from an HDF5 file.
    pub fn read_from_file(file: &Hdf5File, label: &str) -> Hdf5Result<Self> {
        let base = format!("checkpoints/{label}");
        let group = file.open_group(&base)?;

        let pos_ds = group.open_dataset("positions")?;
        let positions = pos_ds.read_f64()?;
        let velocities = group.open_dataset("velocities")?.read_f64()?;
        let box_flat = group.open_dataset("box")?.read_f64()?;

        // Read step attribute
        let step = match pos_ds.get_attr("step")? {
            AttrValue::Int32(v) => *v as u64,
            _ => 0,
        };
        let time = match pos_ds.get_attr("time")? {
            AttrValue::Float64(v) => *v,
            _ => 0.0,
        };

        let mut box_vectors = [0.0_f64; 9];
        for (i, &v) in box_flat.iter().take(9).enumerate() {
            box_vectors[i] = v;
        }

        Ok(Self {
            label: label.to_string(),
            positions,
            velocities,
            step,
            time,
            box_vectors,
            scalars: HashMap::new(),
        })
    }
}

// ---------------------------------------------------------------------------
// Virtual dataset (multi-file)
// ---------------------------------------------------------------------------

/// A virtual dataset that concatenates slices from multiple source files.
#[derive(Debug, Clone)]
pub struct VirtualDataset {
    /// Final shape of the combined virtual dataset.
    pub shape: Vec<usize>,
    /// Sources in order: `(external_ref, hyperslab_into_vds)`.
    pub sources: Vec<(ExternalRef, Hyperslab)>,
    /// Datatype of the virtual dataset.
    pub dtype: Hdf5Dtype,
}

impl VirtualDataset {
    /// Create a new virtual dataset definition.
    pub fn new(shape: Vec<usize>, dtype: Hdf5Dtype) -> Self {
        Self {
            shape,
            sources: Vec::new(),
            dtype,
        }
    }

    /// Add a source contributing to a slice of the virtual dataset.
    pub fn add_source(&mut self, ext_ref: ExternalRef, slab: Hyperslab) {
        self.sources.push((ext_ref, slab));
    }

    /// Simulate resolving the virtual dataset (returns a description string).
    pub fn resolve_description(&self) -> String {
        let src_strs: Vec<String> = self
            .sources
            .iter()
            .map(|(r, s)| {
                format!(
                    "{}:{}@offset={}->slab_vol={}",
                    r.filename,
                    r.dataset_path,
                    r.byte_offset,
                    s.volume()
                )
            })
            .collect();
        format!(
            "VDS shape={:?} sources=[{}]",
            self.shape,
            src_strs.join(", ")
        )
    }
}

// ---------------------------------------------------------------------------
// Checkpoint manager
// ---------------------------------------------------------------------------

/// Manages a rolling window of checkpoints in a single HDF5 file.
#[derive(Debug, Clone)]
pub struct CheckpointManager {
    /// Maximum number of checkpoints to keep.
    pub max_checkpoints: usize,
    /// Labels of currently stored checkpoints (oldest first).
    pub checkpoint_labels: Vec<String>,
}

impl CheckpointManager {
    /// Create a new manager with the given capacity.
    pub fn new(max_checkpoints: usize) -> Self {
        Self {
            max_checkpoints,
            checkpoint_labels: Vec::new(),
        }
    }

    /// Write a checkpoint, evicting the oldest if capacity is exceeded.
    pub fn write(&mut self, file: &mut Hdf5File, ckpt: &Checkpoint) -> Hdf5Result<()> {
        ckpt.write_to_file(file)?;
        self.checkpoint_labels.push(ckpt.label.clone());
        if self.checkpoint_labels.len() > self.max_checkpoints {
            let _old = self.checkpoint_labels.remove(0);
            // In a real implementation we would delete the group from the file.
        }
        Ok(())
    }

    /// Read the most recent checkpoint.
    pub fn read_latest(&self, file: &Hdf5File) -> Hdf5Result<Checkpoint> {
        let label = self
            .checkpoint_labels
            .last()
            .ok_or_else(|| Hdf5Error::Generic("no checkpoints stored".to_string()))?;
        Checkpoint::read_from_file(file, label)
    }

    /// Number of stored checkpoints.
    pub fn count(&self) -> usize {
        self.checkpoint_labels.len()
    }
}
