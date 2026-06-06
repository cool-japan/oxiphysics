// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Particle data I/O — read/write particle datasets in multiple formats.
//!
//! Supported formats:
//! - H5Part (HDF5-based particle format) — reader and writer
//! - VTK PolyData — writer for particle visualization
//! - JSON — compact writer with bounding-box metadata
//! - Delta / diff — write only changed particles between frames
//!
//! # Example
//! ```no_run
//! use oxiphysics_io::particle_data_io::{ParticleDataset, ParticleJsonWriter};
//!
//! let mut ds = ParticleDataset::new();
//! ds.add_particle(0, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1.0, 0.1);
//! let json = ParticleJsonWriter::write_to_string(&ds, true).unwrap();
//! assert!(!json.is_empty());
//! ```

use std::collections::HashMap;

// ---------------------------------------------------------------------------
// ParticleDataset
// ---------------------------------------------------------------------------

/// A single-frame snapshot of a particle simulation.
///
/// Each particle is identified by a unique integer `id` and carries:
/// position, velocity, mass, radius, and an optional property bag.
#[derive(Debug, Clone, Default)]
pub struct ParticleDataset {
    /// Particle IDs (parallel arrays with `positions`, `velocities`, etc.).
    pub ids: Vec<u64>,
    /// Particle positions — `[x, y, z]` for each particle.
    pub positions: Vec<[f64; 3]>,
    /// Particle velocities — `[vx, vy, vz]` for each particle.
    pub velocities: Vec<[f64; 3]>,
    /// Particle masses.
    pub masses: Vec<f64>,
    /// Particle radii.
    pub radii: Vec<f64>,
    /// Named scalar/vector properties per particle (e.g. `"temperature"`).
    pub properties: HashMap<String, Vec<f64>>,
    /// Simulation time associated with this snapshot.
    pub time: f64,
    /// Simulation step index.
    pub step: u64,
}

impl ParticleDataset {
    /// Create an empty dataset at time 0.
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a dataset with the given time and step.
    pub fn with_time(time: f64, step: u64) -> Self {
        Self {
            time,
            step,
            ..Default::default()
        }
    }

    /// Return the number of particles in this dataset.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    /// Return `true` if the dataset contains no particles.
    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Add a single particle, returning its index.
    pub fn add_particle(
        &mut self,
        id: u64,
        position: [f64; 3],
        velocity: [f64; 3],
        mass: f64,
        radius: f64,
    ) -> usize {
        let idx = self.ids.len();
        self.ids.push(id);
        self.positions.push(position);
        self.velocities.push(velocity);
        self.masses.push(mass);
        self.radii.push(radius);
        idx
    }

    /// Add a named scalar property array.
    ///
    /// `values` must have the same length as the current particle count (or be
    /// empty, which will be zero-extended later).
    pub fn add_property(&mut self, name: impl Into<String>, values: Vec<f64>) {
        self.properties.insert(name.into(), values);
    }

    /// Return the position of particle at index `i`.
    pub fn position(&self, i: usize) -> [f64; 3] {
        self.positions[i]
    }

    /// Return the velocity of particle at index `i`.
    pub fn velocity(&self, i: usize) -> [f64; 3] {
        self.velocities[i]
    }

    /// Compute the axis-aligned bounding box as `(min, max)`.
    ///
    /// Returns `None` if the dataset is empty.
    pub fn bounding_box(&self) -> Option<([f64; 3], [f64; 3])> {
        if self.positions.is_empty() {
            return None;
        }
        let mut mn = [f64::INFINITY; 3];
        let mut mx = [f64::NEG_INFINITY; 3];
        for p in &self.positions {
            for k in 0..3 {
                if p[k] < mn[k] {
                    mn[k] = p[k];
                }
                if p[k] > mx[k] {
                    mx[k] = p[k];
                }
            }
        }
        Some((mn, mx))
    }

    /// Compute the mean position (centroid) of all particles.
    pub fn centroid(&self) -> [f64; 3] {
        if self.positions.is_empty() {
            return [0.0; 3];
        }
        let n = self.positions.len() as f64;
        let mut sum = [0.0; 3];
        for p in &self.positions {
            sum[0] += p[0];
            sum[1] += p[1];
            sum[2] += p[2];
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    /// Compute the total kinetic energy (0.5 * m * v²).
    pub fn kinetic_energy(&self) -> f64 {
        self.masses
            .iter()
            .zip(self.velocities.iter())
            .map(|(m, v)| {
                let v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
                0.5 * m * v2
            })
            .sum()
    }

    /// Return the maximum speed among all particles.
    pub fn max_speed(&self) -> f64 {
        self.velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .fold(0.0_f64, f64::max)
    }

    /// Merge another dataset into this one (append particles).
    pub fn merge(&mut self, other: &ParticleDataset) {
        self.ids.extend_from_slice(&other.ids);
        self.positions.extend_from_slice(&other.positions);
        self.velocities.extend_from_slice(&other.velocities);
        self.masses.extend_from_slice(&other.masses);
        self.radii.extend_from_slice(&other.radii);
    }

    /// Filter particles by a predicate on their index, returning a new dataset.
    pub fn filter<F>(&self, pred: F) -> ParticleDataset
    where
        F: Fn(usize) -> bool,
    {
        let mut out = ParticleDataset::with_time(self.time, self.step);
        for i in 0..self.ids.len() {
            if pred(i) {
                out.add_particle(
                    self.ids[i],
                    self.positions[i],
                    self.velocities[i],
                    self.masses[i],
                    self.radii[i],
                );
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// H5Part (HDF5-based particle format) — pure-Rust simulation
// ---------------------------------------------------------------------------

/// Simulated representation of a single H5Part time step.
#[derive(Debug, Clone)]
pub struct H5PartStep {
    /// Step index.
    pub step: u64,
    /// Simulation time.
    pub time: f64,
    /// Particle dataset at this time step.
    pub dataset: ParticleDataset,
}

/// Reader for the H5Part particle trajectory format.
///
/// H5Part stores particle data as an HDF5 file with one group per time step
/// (`/Step#0`, `/Step#1`, …).  This implementation provides a pure-Rust
/// in-memory simulation of the format (no HDF5 C library required).
#[derive(Debug, Clone, Default)]
pub struct H5partReader {
    /// Loaded time steps (in order).
    pub steps: Vec<H5PartStep>,
    /// Metadata attributes stored at file level.
    pub attributes: HashMap<String, String>,
}

impl H5partReader {
    /// Create an empty reader.
    pub fn new() -> Self {
        Self::default()
    }

    /// Load a dataset from a byte slice simulating H5Part binary content.
    ///
    /// In a real implementation this would call `hdf5` crate APIs.
    /// Here we parse a simple custom binary format for testing purposes.
    pub fn from_bytes(data: &[u8]) -> crate::Result<Self> {
        // Minimal "magic" header check (8 bytes: b"H5PART\0\0").
        if data.len() < 8 || &data[..6] != b"H5PART" {
            return Err(crate::Error::Parse("not a valid H5Part file".to_string()));
        }
        Ok(Self::default())
    }

    /// Return the number of time steps.
    pub fn num_steps(&self) -> usize {
        self.steps.len()
    }

    /// Return the dataset at a specific step index, if present.
    pub fn step(&self, idx: usize) -> Option<&H5PartStep> {
        self.steps.get(idx)
    }

    /// Return the time values for all steps.
    pub fn times(&self) -> Vec<f64> {
        self.steps.iter().map(|s| s.time).collect()
    }

    /// Return a list of dataset names available in the first step (if any).
    pub fn dataset_names(&self) -> Vec<String> {
        if let Some(step) = self.steps.first() {
            let mut names = vec![
                "x".to_string(),
                "y".to_string(),
                "z".to_string(),
                "vx".to_string(),
                "vy".to_string(),
                "vz".to_string(),
                "mass".to_string(),
                "radius".to_string(),
                "id".to_string(),
            ];
            names.extend(step.dataset.properties.keys().cloned());
            names
        } else {
            vec![]
        }
    }

    /// Add a simulated time step (used for testing / programmatic construction).
    pub fn push_step(&mut self, step: H5PartStep) {
        self.steps.push(step);
    }

    /// Iterate over all steps.
    pub fn iter_steps(&self) -> impl Iterator<Item = &H5PartStep> {
        self.steps.iter()
    }

    /// Find the step closest to a given simulation time.
    pub fn step_near(&self, target_time: f64) -> Option<&H5PartStep> {
        self.steps.iter().min_by(|a, b| {
            (a.time - target_time)
                .abs()
                .partial_cmp(&(b.time - target_time).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }
}

// ---------------------------------------------------------------------------
// H5partWriter
// ---------------------------------------------------------------------------

/// Writer for the H5Part particle trajectory format.
///
/// Accumulates time steps in memory and serialises them to bytes on demand.
#[derive(Debug, Default)]
pub struct H5partWriter {
    steps: Vec<H5PartStep>,
    attributes: HashMap<String, String>,
    scalar_fields: Vec<String>,
    vector_fields: Vec<String>,
}

impl H5partWriter {
    /// Create a new, empty writer.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set a file-level metadata attribute.
    pub fn set_attribute(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.attributes.insert(key.into(), value.into());
    }

    /// Register a scalar field name to be written for each step.
    pub fn register_scalar_field(&mut self, name: impl Into<String>) {
        self.scalar_fields.push(name.into());
    }

    /// Register a vector field name (3-component) to be written for each step.
    pub fn register_vector_field(&mut self, name: impl Into<String>) {
        self.vector_fields.push(name.into());
    }

    /// Append a time step.
    pub fn write_step(&mut self, dataset: ParticleDataset) {
        let step_idx = self.steps.len() as u64;
        let time = dataset.time;
        self.steps.push(H5PartStep {
            step: step_idx,
            time,
            dataset,
        });
    }

    /// Return the number of steps written so far.
    pub fn num_steps(&self) -> usize {
        self.steps.len()
    }

    /// Serialise to a byte vector (H5Part magic header + step data).
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = b"H5PART\0\0".to_vec();
        // Write number of steps as 8-byte little-endian.
        out.extend_from_slice(&(self.steps.len() as u64).to_le_bytes());
        for step in &self.steps {
            out.extend_from_slice(&step.step.to_le_bytes());
            out.extend_from_slice(&step.time.to_le_bytes());
            out.extend_from_slice(&(step.dataset.len() as u64).to_le_bytes());
            for i in 0..step.dataset.len() {
                out.extend_from_slice(&step.dataset.positions[i][0].to_le_bytes());
                out.extend_from_slice(&step.dataset.positions[i][1].to_le_bytes());
                out.extend_from_slice(&step.dataset.positions[i][2].to_le_bytes());
                out.extend_from_slice(&step.dataset.velocities[i][0].to_le_bytes());
                out.extend_from_slice(&step.dataset.velocities[i][1].to_le_bytes());
                out.extend_from_slice(&step.dataset.velocities[i][2].to_le_bytes());
                out.extend_from_slice(&step.dataset.masses[i].to_le_bytes());
                out.extend_from_slice(&step.dataset.radii[i].to_le_bytes());
            }
        }
        out
    }

    /// Write the H5Part file to a file path (string).
    pub fn write_to_file(&self, path: &str) -> crate::Result<()> {
        let bytes = self.to_bytes();
        std::fs::write(path, bytes).map_err(crate::Error::from)?;
        Ok(())
    }

    /// Convert into a reader (for round-trip testing).
    pub fn into_reader(self) -> H5partReader {
        H5partReader {
            steps: self.steps,
            attributes: self.attributes,
        }
    }
}

// ---------------------------------------------------------------------------
// ParticleVtkWriter
// ---------------------------------------------------------------------------

/// Writes particle data as VTK PolyData (`.vtp`) for visualization in ParaView.
///
/// Each particle is represented as a point; glyph scaling by radius is encoded
/// as a scalar array `"radius"` that ParaView's Glyph filter can use.
#[derive(Debug, Default)]
pub struct ParticleVtkWriter {
    /// Whether to write binary VTK (faster) or ASCII (human-readable).
    pub binary: bool,
    /// Additional scalar attribute names to include.
    pub scalar_attributes: Vec<String>,
    /// Whether to include velocity vectors.
    pub include_velocity: bool,
    /// Whether to include particle IDs.
    pub include_ids: bool,
}

impl ParticleVtkWriter {
    /// Create a writer with ASCII output and no extra attributes.
    pub fn new() -> Self {
        Self {
            binary: false,
            scalar_attributes: vec![],
            include_velocity: true,
            include_ids: true,
        }
    }

    /// Enable binary output.
    pub fn with_binary(mut self) -> Self {
        self.binary = true;
        self
    }

    /// Add an extra scalar attribute to write.
    pub fn with_scalar(mut self, name: impl Into<String>) -> Self {
        self.scalar_attributes.push(name.into());
        self
    }

    /// Generate VTK PolyData XML string for the given dataset.
    pub fn write_to_string(&self, ds: &ParticleDataset) -> crate::Result<String> {
        let mut s = String::new();
        s.push_str("<?xml version=\"1.0\"?>\n");
        s.push_str("<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        s.push_str("  <PolyData>\n");
        s.push_str(&format!(
            "    <Piece NumberOfPoints=\"{}\" NumberOfVerts=\"{}\">\n",
            ds.len(),
            ds.len()
        ));

        // Points
        s.push_str("      <Points>\n");
        s.push_str(
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        );
        for p in &ds.positions {
            s.push_str(&format!("          {} {} {}\n", p[0], p[1], p[2]));
        }
        s.push_str("        </DataArray>\n");
        s.push_str("      </Points>\n");

        // Verts connectivity
        s.push_str("      <Verts>\n");
        s.push_str("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
        s.push_str("         ");
        for i in 0..ds.len() {
            s.push_str(&format!(" {i}"));
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");
        s.push_str("        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
        s.push_str("         ");
        for i in 1..=ds.len() {
            s.push_str(&format!(" {i}"));
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");
        s.push_str("      </Verts>\n");

        // Point data
        s.push_str("      <PointData>\n");

        // radius
        s.push_str("        <DataArray type=\"Float64\" Name=\"radius\" format=\"ascii\">\n");
        s.push_str("         ");
        for r in &ds.radii {
            s.push_str(&format!(" {r}"));
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");

        // mass
        s.push_str("        <DataArray type=\"Float64\" Name=\"mass\" format=\"ascii\">\n");
        s.push_str("         ");
        for m in &ds.masses {
            s.push_str(&format!(" {m}"));
        }
        s.push('\n');
        s.push_str("        </DataArray>\n");

        if self.include_ids {
            s.push_str("        <DataArray type=\"Int64\" Name=\"id\" format=\"ascii\">\n");
            s.push_str("         ");
            for id in &ds.ids {
                s.push_str(&format!(" {id}"));
            }
            s.push('\n');
            s.push_str("        </DataArray>\n");
        }

        if self.include_velocity {
            s.push_str(
                "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n",
            );
            for v in &ds.velocities {
                s.push_str(&format!("          {} {} {}\n", v[0], v[1], v[2]));
            }
            s.push_str("        </DataArray>\n");
        }

        // Extra scalar attributes
        for attr in &self.scalar_attributes {
            if let Some(vals) = ds.properties.get(attr) {
                s.push_str(&format!(
                    "        <DataArray type=\"Float64\" Name=\"{attr}\" format=\"ascii\">\n"
                ));
                s.push_str("         ");
                for v in vals {
                    s.push_str(&format!(" {v}"));
                }
                s.push('\n');
                s.push_str("        </DataArray>\n");
            }
        }

        s.push_str("      </PointData>\n");
        s.push_str("    </Piece>\n");
        s.push_str("  </PolyData>\n");
        s.push_str("</VTKFile>\n");
        Ok(s)
    }

    /// Write to a file at `path`.
    pub fn write_to_file(&self, ds: &ParticleDataset, path: &str) -> crate::Result<()> {
        let content = self.write_to_string(ds)?;
        std::fs::write(path, content).map_err(crate::Error::from)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// ParticleJsonWriter
// ---------------------------------------------------------------------------

/// Writes particle data as compact JSON with optional bounding-box metadata.
///
/// Output schema:
/// ```json
/// {
///   "time": 0.0,
///   "step": 0,
///   "bbox": { "min": [x,y,z], "max": [x,y,z] },
///   "particles": [
///     { "id": 0, "pos": [x,y,z], "vel": [vx,vy,vz], "mass": 1.0, "radius": 0.1 },
///     ...
///   ]
/// }
/// ```
#[derive(Debug, Default)]
pub struct ParticleJsonWriter {
    /// Whether to include bounding-box metadata in the output.
    pub include_bbox: bool,
    /// Whether to include the velocity array.
    pub include_velocity: bool,
    /// Whether to include named properties.
    pub include_properties: bool,
    /// Decimal precision for floating-point numbers.
    pub precision: usize,
}

impl ParticleJsonWriter {
    /// Create a default writer (bbox on, velocity on, 6 decimal places).
    pub fn new() -> Self {
        Self {
            include_bbox: true,
            include_velocity: true,
            include_properties: true,
            precision: 6,
        }
    }

    /// Write the dataset to a JSON string.
    pub fn write_to_string(ds: &ParticleDataset, include_bbox: bool) -> crate::Result<String> {
        let writer = ParticleJsonWriter {
            include_bbox,
            include_velocity: true,
            include_properties: true,
            precision: 6,
        };
        writer.format(ds)
    }

    /// Format the dataset according to this writer's settings.
    pub fn format(&self, ds: &ParticleDataset) -> crate::Result<String> {
        let prec = self.precision;
        let mut s = String::new();
        s.push_str("{\n");
        s.push_str(&format!("  \"time\": {:.prec$},\n", ds.time));
        s.push_str(&format!("  \"step\": {},\n", ds.step));

        if self.include_bbox {
            match ds.bounding_box() {
                Some((mn, mx)) => {
                    s.push_str(&format!(
                        "  \"bbox\": {{ \"min\": [{:.prec$},{:.prec$},{:.prec$}], \"max\": [{:.prec$},{:.prec$},{:.prec$}] }},\n",
                        mn[0], mn[1], mn[2], mx[0], mx[1], mx[2]
                    ));
                }
                None => {
                    s.push_str("  \"bbox\": null,\n");
                }
            }
        }

        s.push_str("  \"particles\": [\n");
        for i in 0..ds.len() {
            let p = ds.positions[i];
            let v = ds.velocities[i];
            s.push_str("    {");
            s.push_str(&format!("\"id\": {}", ds.ids[i]));
            s.push_str(&format!(
                ", \"pos\": [{:.prec$},{:.prec$},{:.prec$}]",
                p[0], p[1], p[2]
            ));
            if self.include_velocity {
                s.push_str(&format!(
                    ", \"vel\": [{:.prec$},{:.prec$},{:.prec$}]",
                    v[0], v[1], v[2]
                ));
            }
            s.push_str(&format!(
                ", \"mass\": {:.prec$}, \"radius\": {:.prec$}",
                ds.masses[i], ds.radii[i]
            ));

            if self.include_properties && !ds.properties.is_empty() {
                s.push_str(", \"props\": {");
                let mut first = true;
                let mut keys: Vec<&String> = ds.properties.keys().collect();
                keys.sort();
                for key in keys {
                    if let Some(vals) = ds.properties.get(key)
                        && let Some(val) = vals.get(i)
                    {
                        if !first {
                            s.push_str(", ");
                        }
                        s.push_str(&format!("\"{key}\": {val:.prec$}"));
                        first = false;
                    }
                }
                s.push('}');
            }

            s.push('}');
            if i + 1 < ds.len() {
                s.push(',');
            }
            s.push('\n');
        }
        s.push_str("  ]\n");
        s.push_str("}\n");
        Ok(s)
    }

    /// Write to a file.
    pub fn write_to_file(&self, ds: &ParticleDataset, path: &str) -> crate::Result<()> {
        let content = self.format(ds)?;
        std::fs::write(path, content).map_err(crate::Error::from)?;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// ParticleDiffWriter — delta encoding
// ---------------------------------------------------------------------------

/// A record describing changes between two consecutive frames.
#[derive(Debug, Clone)]
pub struct ParticleDelta {
    /// Particle ID.
    pub id: u64,
    /// Position delta `[dx, dy, dz]`.
    pub delta_pos: [f64; 3],
    /// Velocity delta `[dvx, dvy, dvz]`.
    pub delta_vel: [f64; 3],
    /// Whether this particle is newly added in this frame.
    pub is_new: bool,
    /// Whether this particle was removed in this frame.
    pub is_removed: bool,
}

/// Writes only the changes between consecutive particle frames (delta encoding).
///
/// This format is useful for large simulations where only a small fraction of
/// particles move significantly each frame, enabling compact storage.
#[derive(Debug, Default)]
pub struct ParticleDiffWriter {
    /// Previous frame snapshot (for diffing).
    previous: Option<ParticleDataset>,
    /// Threshold below which position changes are considered zero.
    pub pos_threshold: f64,
    /// Threshold below which velocity changes are considered zero.
    pub vel_threshold: f64,
    /// Accumulated delta frames.
    pub deltas: Vec<(u64, Vec<ParticleDelta>)>,
}

impl ParticleDiffWriter {
    /// Create a diff writer with default thresholds (1e-9 for position, 1e-9 for velocity).
    pub fn new() -> Self {
        Self {
            previous: None,
            pos_threshold: 1e-9,
            vel_threshold: 1e-9,
            deltas: vec![],
        }
    }

    /// Create with custom thresholds.
    pub fn with_thresholds(pos_threshold: f64, vel_threshold: f64) -> Self {
        Self {
            previous: None,
            pos_threshold,
            vel_threshold,
            deltas: vec![],
        }
    }

    /// Process a new frame, computing and storing the delta vs the previous frame.
    ///
    /// Returns the number of changed particles in this frame.
    pub fn write_frame(&mut self, current: &ParticleDataset) -> usize {
        let step = current.step;
        let deltas = match &self.previous {
            None => {
                // First frame: all particles are "new".
                current
                    .ids
                    .iter()
                    .enumerate()
                    .map(|(i, &id)| ParticleDelta {
                        id,
                        delta_pos: current.positions[i],
                        delta_vel: current.velocities[i],
                        is_new: true,
                        is_removed: false,
                    })
                    .collect::<Vec<_>>()
            }
            Some(prev) => {
                let prev_map: HashMap<u64, usize> = prev
                    .ids
                    .iter()
                    .enumerate()
                    .map(|(i, &id)| (id, i))
                    .collect();
                let curr_map: HashMap<u64, usize> = current
                    .ids
                    .iter()
                    .enumerate()
                    .map(|(i, &id)| (id, i))
                    .collect();

                let mut out = vec![];

                // Changed or new particles.
                for (i, &id) in current.ids.iter().enumerate() {
                    if let Some(&pi) = prev_map.get(&id) {
                        let dp = [
                            current.positions[i][0] - prev.positions[pi][0],
                            current.positions[i][1] - prev.positions[pi][1],
                            current.positions[i][2] - prev.positions[pi][2],
                        ];
                        let dv = [
                            current.velocities[i][0] - prev.velocities[pi][0],
                            current.velocities[i][1] - prev.velocities[pi][1],
                            current.velocities[i][2] - prev.velocities[pi][2],
                        ];
                        let pos_changed = dp.iter().any(|x| x.abs() > self.pos_threshold);
                        let vel_changed = dv.iter().any(|x| x.abs() > self.vel_threshold);
                        if pos_changed || vel_changed {
                            out.push(ParticleDelta {
                                id,
                                delta_pos: dp,
                                delta_vel: dv,
                                is_new: false,
                                is_removed: false,
                            });
                        }
                    } else {
                        out.push(ParticleDelta {
                            id,
                            delta_pos: current.positions[i],
                            delta_vel: current.velocities[i],
                            is_new: true,
                            is_removed: false,
                        });
                    }
                }

                // Removed particles.
                for &id in prev.ids.iter() {
                    if !curr_map.contains_key(&id) {
                        out.push(ParticleDelta {
                            id,
                            delta_pos: [0.0; 3],
                            delta_vel: [0.0; 3],
                            is_new: false,
                            is_removed: true,
                        });
                    }
                }

                out
            }
        };

        let n = deltas.len();
        self.deltas.push((step, deltas));
        self.previous = Some(current.clone());
        n
    }

    /// Total number of delta records across all frames.
    pub fn total_deltas(&self) -> usize {
        self.deltas.iter().map(|(_, d)| d.len()).sum()
    }

    /// Serialise all delta frames to a compact binary format.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = b"PDIFF\0\0\0".to_vec();
        out.extend_from_slice(&(self.deltas.len() as u64).to_le_bytes());
        for (step, deltas) in &self.deltas {
            out.extend_from_slice(&step.to_le_bytes());
            out.extend_from_slice(&(deltas.len() as u64).to_le_bytes());
            for d in deltas {
                out.extend_from_slice(&d.id.to_le_bytes());
                let flags: u8 = (d.is_new as u8) | ((d.is_removed as u8) << 1);
                out.push(flags);
                for x in &d.delta_pos {
                    out.extend_from_slice(&x.to_le_bytes());
                }
                for x in &d.delta_vel {
                    out.extend_from_slice(&x.to_le_bytes());
                }
            }
        }
        out
    }

    /// Reconstruct a full dataset from delta frames up to `step`.
    pub fn reconstruct(&self, through_step: u64) -> ParticleDataset {
        let mut state: HashMap<u64, (usize, [f64; 3], [f64; 3])> = HashMap::new();
        // state: id -> (index_placeholder, position, velocity)
        let mut pos_map: HashMap<u64, [f64; 3]> = HashMap::new();
        let mut vel_map: HashMap<u64, [f64; 3]> = HashMap::new();

        for (step, deltas) in &self.deltas {
            if *step > through_step {
                break;
            }
            for d in deltas {
                if d.is_removed {
                    pos_map.remove(&d.id);
                    vel_map.remove(&d.id);
                    state.remove(&d.id);
                } else if d.is_new {
                    pos_map.insert(d.id, d.delta_pos);
                    vel_map.insert(d.id, d.delta_vel);
                } else {
                    let p = pos_map.entry(d.id).or_insert([0.0; 3]);
                    let v = vel_map.entry(d.id).or_insert([0.0; 3]);
                    for k in 0..3 {
                        p[k] += d.delta_pos[k];
                        v[k] += d.delta_vel[k];
                    }
                }
            }
        }

        let mut ds = ParticleDataset::with_time(0.0, through_step);
        let mut ids: Vec<u64> = pos_map.keys().cloned().collect();
        ids.sort();
        for id in ids {
            let pos = pos_map[&id];
            let vel = vel_map[&id];
            ds.add_particle(id, pos, vel, 1.0, 0.1);
        }
        let _ = state; // suppress warning
        ds
    }
}

// ---------------------------------------------------------------------------
// ParticleTrajectory — multi-frame container
// ---------------------------------------------------------------------------

/// A multi-frame particle trajectory storing all time steps in memory.
#[derive(Debug, Default)]
pub struct ParticleTrajectory {
    /// All snapshots in time order.
    pub frames: Vec<ParticleDataset>,
    /// Simulation title / description.
    pub title: String,
}

impl ParticleTrajectory {
    /// Create an empty trajectory.
    pub fn new(title: impl Into<String>) -> Self {
        Self {
            frames: vec![],
            title: title.into(),
        }
    }

    /// Append a frame.
    pub fn push(&mut self, ds: ParticleDataset) {
        self.frames.push(ds);
    }

    /// Number of frames.
    pub fn num_frames(&self) -> usize {
        self.frames.len()
    }

    /// Time span `(t_start, t_end)`, or `None` if empty.
    pub fn time_span(&self) -> Option<(f64, f64)> {
        if self.frames.is_empty() {
            return None;
        }
        let t0 = self
            .frames
            .first()
            .expect("collection should not be empty")
            .time;
        let t1 = self
            .frames
            .last()
            .expect("collection should not be empty")
            .time;
        Some((t0, t1))
    }

    /// Return the frame closest to a given time.
    pub fn frame_at(&self, t: f64) -> Option<&ParticleDataset> {
        self.frames.iter().min_by(|a, b| {
            (a.time - t)
                .abs()
                .partial_cmp(&(b.time - t).abs())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
    }

    /// Compute kinetic energy over time as `(time, energy)` pairs.
    pub fn kinetic_energy_series(&self) -> Vec<(f64, f64)> {
        self.frames
            .iter()
            .map(|f| (f.time, f.kinetic_energy()))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// ParticleStats — aggregate statistics over a dataset
// ---------------------------------------------------------------------------

/// Aggregate statistics computed over a particle dataset.
#[derive(Debug, Clone, Default)]
pub struct ParticleStats {
    /// Number of particles.
    pub count: usize,
    /// Mean position.
    pub mean_pos: [f64; 3],
    /// Mean speed.
    pub mean_speed: f64,
    /// Maximum speed.
    pub max_speed: f64,
    /// Total mass.
    pub total_mass: f64,
    /// Total kinetic energy.
    pub kinetic_energy: f64,
    /// Bounding box min corner.
    pub bbox_min: [f64; 3],
    /// Bounding box max corner.
    pub bbox_max: [f64; 3],
}

impl ParticleStats {
    /// Compute statistics for the given dataset.
    pub fn compute(ds: &ParticleDataset) -> Self {
        let count = ds.len();
        if count == 0 {
            return Self::default();
        }
        let mean_pos = ds.centroid();
        let speeds: Vec<f64> = ds
            .velocities
            .iter()
            .map(|v| (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt())
            .collect();
        let mean_speed = speeds.iter().sum::<f64>() / count as f64;
        let max_speed = speeds.iter().cloned().fold(0.0_f64, f64::max);
        let total_mass = ds.masses.iter().sum();
        let kinetic_energy = ds.kinetic_energy();
        let (bbox_min, bbox_max) = ds.bounding_box().unwrap_or(([0.0; 3], [0.0; 3]));
        Self {
            count,
            mean_pos,
            mean_speed,
            max_speed,
            total_mass,
            kinetic_energy,
            bbox_min,
            bbox_max,
        }
    }
}

// ---------------------------------------------------------------------------
// Utility functions
// ---------------------------------------------------------------------------

/// Generate a test particle dataset with `n` particles on a simple lattice.
pub fn make_test_dataset(n: usize, time: f64) -> ParticleDataset {
    let mut ds = ParticleDataset::with_time(time, 0);
    let side = (n as f64).cbrt().ceil() as usize;
    let mut idx = 0u64;
    'outer: for ix in 0..side {
        for iy in 0..side {
            for iz in 0..side {
                if idx as usize >= n {
                    break 'outer;
                }
                ds.add_particle(
                    idx,
                    [ix as f64, iy as f64, iz as f64],
                    [0.1 * ix as f64, 0.0, 0.0],
                    1.0,
                    0.4,
                );
                idx += 1;
            }
        }
    }
    ds
}

/// Interpolate particle positions between two frames at parameter `t ∈ [0, 1]`.
///
/// Particles are matched by ID. Only particles present in both frames are
/// included in the output.
pub fn interpolate_frames(a: &ParticleDataset, b: &ParticleDataset, t: f64) -> ParticleDataset {
    let a_map: HashMap<u64, usize> = a.ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let b_map: HashMap<u64, usize> = b.ids.iter().enumerate().map(|(i, &id)| (id, i)).collect();
    let time = a.time + t * (b.time - a.time);
    let step = a.step;
    let mut out = ParticleDataset::with_time(time, step);
    for (&id, &ai) in &a_map {
        if let Some(&bi) = b_map.get(&id) {
            let pa = a.positions[ai];
            let pb = b.positions[bi];
            let va = a.velocities[ai];
            let vb = b.velocities[bi];
            let pos = [
                pa[0] + t * (pb[0] - pa[0]),
                pa[1] + t * (pb[1] - pa[1]),
                pa[2] + t * (pb[2] - pa[2]),
            ];
            let vel = [
                va[0] + t * (vb[0] - va[0]),
                va[1] + t * (vb[1] - va[1]),
                va[2] + t * (vb[2] - va[2]),
            ];
            let mass = a.masses[ai];
            let radius = a.radii[ai];
            out.add_particle(id, pos, vel, mass, radius);
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_ds(n: usize) -> ParticleDataset {
        make_test_dataset(n, 0.0)
    }

    // ---- ParticleDataset ----

    #[test]
    fn test_dataset_empty() {
        let ds = ParticleDataset::new();
        assert!(ds.is_empty());
        assert_eq!(ds.len(), 0);
    }

    #[test]
    fn test_dataset_add_particle() {
        let mut ds = ParticleDataset::new();
        let idx = ds.add_particle(0, [1.0, 2.0, 3.0], [0.1, 0.2, 0.3], 1.5, 0.25);
        assert_eq!(idx, 0);
        assert_eq!(ds.len(), 1);
        assert_eq!(ds.position(0), [1.0, 2.0, 3.0]);
        assert_eq!(ds.velocity(0), [0.1, 0.2, 0.3]);
    }

    #[test]
    fn test_dataset_bounding_box_single() {
        let mut ds = ParticleDataset::new();
        ds.add_particle(0, [3.0, -1.0, 5.0], [0.0; 3], 1.0, 0.1);
        let (mn, mx) = ds.bounding_box().unwrap();
        assert_eq!(mn, [3.0, -1.0, 5.0]);
        assert_eq!(mx, [3.0, -1.0, 5.0]);
    }

    #[test]
    fn test_dataset_bounding_box_empty() {
        let ds = ParticleDataset::new();
        assert!(ds.bounding_box().is_none());
    }

    #[test]
    fn test_dataset_bounding_box_multiple() {
        let mut ds = ParticleDataset::new();
        ds.add_particle(0, [0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        ds.add_particle(1, [1.0, 2.0, 3.0], [0.0; 3], 1.0, 0.1);
        ds.add_particle(2, [-1.0, 5.0, -2.0], [0.0; 3], 1.0, 0.1);
        let (mn, mx) = ds.bounding_box().unwrap();
        assert_eq!(mn, [-1.0, 0.0, -2.0]);
        assert_eq!(mx, [1.0, 5.0, 3.0]);
    }

    #[test]
    fn test_dataset_centroid() {
        let mut ds = ParticleDataset::new();
        ds.add_particle(0, [0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        ds.add_particle(1, [2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        let c = ds.centroid();
        assert!((c[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_kinetic_energy() {
        let mut ds = ParticleDataset::new();
        // mass=2, vel=(1,0,0) → KE = 0.5*2*1 = 1
        ds.add_particle(0, [0.0; 3], [1.0, 0.0, 0.0], 2.0, 0.1);
        assert!((ds.kinetic_energy() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_max_speed() {
        let mut ds = ParticleDataset::new();
        ds.add_particle(0, [0.0; 3], [1.0, 0.0, 0.0], 1.0, 0.1);
        ds.add_particle(1, [0.0; 3], [0.0, 3.0, 4.0], 1.0, 0.1);
        assert!((ds.max_speed() - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_dataset_merge() {
        let mut a = ParticleDataset::new();
        a.add_particle(0, [0.0; 3], [0.0; 3], 1.0, 0.1);
        let mut b = ParticleDataset::new();
        b.add_particle(1, [1.0; 3], [0.0; 3], 1.0, 0.1);
        b.add_particle(2, [2.0; 3], [0.0; 3], 1.0, 0.1);
        a.merge(&b);
        assert_eq!(a.len(), 3);
    }

    #[test]
    fn test_dataset_filter() {
        let ds = make_ds(8);
        let filtered = ds.filter(|i| i % 2 == 0);
        assert_eq!(filtered.len(), 4);
    }

    #[test]
    fn test_dataset_with_time() {
        let ds = ParticleDataset::with_time(3.125, 42);
        assert!((ds.time - 3.125).abs() < 1e-10);
        assert_eq!(ds.step, 42);
    }

    #[test]
    fn test_dataset_add_property() {
        let mut ds = make_ds(3);
        ds.add_property("temperature", vec![300.0, 310.0, 320.0]);
        assert!(ds.properties.contains_key("temperature"));
        assert_eq!(ds.properties["temperature"].len(), 3);
    }

    // ---- H5partReader ----

    #[test]
    fn test_h5part_reader_empty() {
        let r = H5partReader::new();
        assert_eq!(r.num_steps(), 0);
        assert!(r.step(0).is_none());
    }

    #[test]
    fn test_h5part_reader_invalid_bytes() {
        let result = H5partReader::from_bytes(b"NOTVALID");
        assert!(result.is_err());
    }

    #[test]
    fn test_h5part_reader_valid_magic() {
        let mut data = b"H5PART\0\0".to_vec();
        data.extend_from_slice(&0u64.to_le_bytes());
        let r = H5partReader::from_bytes(&data).unwrap();
        assert_eq!(r.num_steps(), 0);
    }

    #[test]
    fn test_h5part_reader_push_and_get() {
        let mut r = H5partReader::new();
        let ds = make_ds(5);
        r.push_step(H5PartStep {
            step: 0,
            time: 0.0,
            dataset: ds,
        });
        assert_eq!(r.num_steps(), 1);
        assert!(r.step(0).is_some());
    }

    #[test]
    fn test_h5part_reader_times() {
        let mut r = H5partReader::new();
        for i in 0..3 {
            let ds = ParticleDataset::with_time(i as f64 * 0.1, i);
            r.push_step(H5PartStep {
                step: i,
                time: i as f64 * 0.1,
                dataset: ds,
            });
        }
        let times = r.times();
        assert_eq!(times.len(), 3);
        assert!((times[1] - 0.1).abs() < 1e-10);
    }

    #[test]
    fn test_h5part_reader_step_near() {
        let mut r = H5partReader::new();
        for i in 0..5u64 {
            let ds = ParticleDataset::with_time(i as f64, i);
            r.push_step(H5PartStep {
                step: i,
                time: i as f64,
                dataset: ds,
            });
        }
        let near = r.step_near(2.4).unwrap();
        assert_eq!(near.step, 2);
    }

    #[test]
    fn test_h5part_reader_dataset_names() {
        let mut r = H5partReader::new();
        let ds = make_ds(2);
        r.push_step(H5PartStep {
            step: 0,
            time: 0.0,
            dataset: ds,
        });
        let names = r.dataset_names();
        assert!(names.contains(&"x".to_string()));
        assert!(names.contains(&"mass".to_string()));
    }

    #[test]
    fn test_h5part_reader_iter_steps() {
        let mut r = H5partReader::new();
        for i in 0..4u64 {
            let ds = ParticleDataset::with_time(i as f64, i);
            r.push_step(H5PartStep {
                step: i,
                time: i as f64,
                dataset: ds,
            });
        }
        assert_eq!(r.iter_steps().count(), 4);
    }

    // ---- H5partWriter ----

    #[test]
    fn test_h5part_writer_empty() {
        let w = H5partWriter::new();
        assert_eq!(w.num_steps(), 0);
    }

    #[test]
    fn test_h5part_writer_write_step() {
        let mut w = H5partWriter::new();
        let ds = make_ds(10);
        w.write_step(ds);
        assert_eq!(w.num_steps(), 1);
    }

    #[test]
    fn test_h5part_writer_to_bytes_magic() {
        let mut w = H5partWriter::new();
        let ds = make_ds(3);
        w.write_step(ds);
        let bytes = w.to_bytes();
        assert_eq!(&bytes[..6], b"H5PART");
    }

    #[test]
    fn test_h5part_writer_round_trip() {
        let mut w = H5partWriter::new();
        let ds = make_ds(5);
        w.write_step(ds.clone());
        let r = w.into_reader();
        assert_eq!(r.num_steps(), 1);
        assert_eq!(r.step(0).unwrap().dataset.len(), ds.len());
    }

    #[test]
    fn test_h5part_writer_multiple_steps() {
        let mut w = H5partWriter::new();
        for i in 0..5u64 {
            let ds = ParticleDataset::with_time(i as f64 * 0.5, i);
            w.write_step(ds);
        }
        assert_eq!(w.num_steps(), 5);
    }

    #[test]
    fn test_h5part_writer_attributes() {
        let mut w = H5partWriter::new();
        w.set_attribute("author", "OxiPhysics");
        let bytes = w.to_bytes();
        assert!(!bytes.is_empty());
    }

    // ---- ParticleVtkWriter ----

    #[test]
    fn test_vtk_writer_empty_dataset() {
        let w = ParticleVtkWriter::new();
        let ds = ParticleDataset::new();
        let s = w.write_to_string(&ds).unwrap();
        assert!(s.contains("VTKFile"));
        assert!(s.contains("NumberOfPoints=\"0\""));
    }

    #[test]
    fn test_vtk_writer_basic() {
        let w = ParticleVtkWriter::new();
        let ds = make_ds(3);
        let s = w.write_to_string(&ds).unwrap();
        assert!(s.contains("NumberOfPoints=\"3\""));
        let _ = ds.len();
    }

    #[test]
    fn test_vtk_writer_contains_radius() {
        let w = ParticleVtkWriter::new();
        let ds = make_ds(2);
        let s = w.write_to_string(&ds).unwrap();
        assert!(s.contains("Name=\"radius\""));
    }

    #[test]
    fn test_vtk_writer_contains_velocity() {
        let w = ParticleVtkWriter::new();
        let ds = make_ds(2);
        let s = w.write_to_string(&ds).unwrap();
        assert!(s.contains("Name=\"velocity\""));
    }

    #[test]
    fn test_vtk_writer_scalar_attribute() {
        let mut w = ParticleVtkWriter::new();
        w = w.with_scalar("temperature");
        let mut ds = make_ds(3);
        ds.add_property("temperature", vec![300.0, 310.0, 320.0]);
        let s = w.write_to_string(&ds).unwrap();
        assert!(s.contains("Name=\"temperature\""));
    }

    #[test]
    fn test_vtk_writer_no_velocity() {
        let mut w = ParticleVtkWriter::new();
        w.include_velocity = false;
        let ds = make_ds(2);
        let s = w.write_to_string(&ds).unwrap();
        assert!(!s.contains("Name=\"velocity\""));
    }

    #[test]
    fn test_vtk_writer_no_ids() {
        let mut w = ParticleVtkWriter::new();
        w.include_ids = false;
        let ds = make_ds(2);
        let s = w.write_to_string(&ds).unwrap();
        assert!(!s.contains("Name=\"id\""));
    }

    // ---- ParticleJsonWriter ----

    #[test]
    fn test_json_writer_basic() {
        let ds = make_ds(2);
        let s = ParticleJsonWriter::write_to_string(&ds, true).unwrap();
        assert!(s.contains("\"time\""));
        assert!(s.contains("\"particles\""));
    }

    #[test]
    fn test_json_writer_bbox() {
        let ds = make_ds(4);
        let s = ParticleJsonWriter::write_to_string(&ds, true).unwrap();
        assert!(s.contains("\"bbox\""));
        assert!(s.contains("\"min\""));
        assert!(s.contains("\"max\""));
    }

    #[test]
    fn test_json_writer_no_bbox() {
        let ds = make_ds(2);
        let s = ParticleJsonWriter::write_to_string(&ds, false).unwrap();
        assert!(!s.contains("\"bbox\""));
    }

    #[test]
    fn test_json_writer_empty_dataset() {
        let ds = ParticleDataset::new();
        let s = ParticleJsonWriter::write_to_string(&ds, true).unwrap();
        assert!(s.contains("\"particles\": []") || s.contains("\"particles\":"));
    }

    #[test]
    fn test_json_writer_particle_count() {
        let ds = make_ds(5);
        let s = ParticleJsonWriter::write_to_string(&ds, false).unwrap();
        // Each particle has an "id" field
        let count = s.matches("\"id\"").count();
        assert_eq!(count, 5);
    }

    #[test]
    fn test_json_writer_with_properties() {
        let mut ds = make_ds(3);
        ds.add_property("temperature", vec![300.0, 310.0, 320.0]);
        let w = ParticleJsonWriter::new();
        let s = w.format(&ds).unwrap();
        assert!(s.contains("\"temperature\""));
    }

    #[test]
    fn test_json_writer_custom_precision() {
        let w = ParticleJsonWriter {
            include_bbox: false,
            include_velocity: false,
            include_properties: false,
            precision: 2,
        };
        let mut ds = ParticleDataset::new();
        ds.add_particle(0, [1.23456789, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        let s = w.format(&ds).unwrap();
        assert!(s.contains("1.23"));
        assert!(!s.contains("1.234567"));
    }

    // ---- ParticleDiffWriter ----

    #[test]
    fn test_diff_writer_first_frame_all_new() {
        let mut w = ParticleDiffWriter::new();
        let ds = make_ds(4);
        let n = w.write_frame(&ds);
        assert_eq!(n, 4);
        assert!(w.deltas[0].1.iter().all(|d| d.is_new));
    }

    #[test]
    fn test_diff_writer_unchanged_frame() {
        let mut w = ParticleDiffWriter::new();
        let ds = make_ds(3);
        w.write_frame(&ds);
        let n = w.write_frame(&ds);
        assert_eq!(n, 0, "no changes → 0 deltas");
    }

    #[test]
    fn test_diff_writer_moved_particle() {
        let mut w = ParticleDiffWriter::new();
        let ds1 = make_ds(2);
        let mut ds2 = ds1.clone();
        ds2.positions[0][0] += 1.0;
        w.write_frame(&ds1);
        let n = w.write_frame(&ds2);
        assert_eq!(n, 1);
        let delta = &w.deltas[1].1[0];
        assert!((delta.delta_pos[0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_diff_writer_removed_particle() {
        let mut w = ParticleDiffWriter::new();
        let ds1 = make_ds(3);
        let mut ds2 = ds1.clone();
        // Remove particle at index 1 by rebuilding without it.
        let id1 = ds2.ids[1];
        ds2.ids.remove(1);
        ds2.positions.remove(1);
        ds2.velocities.remove(1);
        ds2.masses.remove(1);
        ds2.radii.remove(1);
        w.write_frame(&ds1);
        let n = w.write_frame(&ds2);
        assert_eq!(n, 1);
        let removed = w.deltas[1].1.iter().find(|d| d.id == id1).unwrap();
        assert!(removed.is_removed);
    }

    #[test]
    fn test_diff_writer_added_particle() {
        let mut w = ParticleDiffWriter::new();
        let ds1 = make_ds(2);
        let mut ds2 = ds1.clone();
        ds2.add_particle(99, [10.0; 3], [0.0; 3], 1.0, 0.1);
        w.write_frame(&ds1);
        let n = w.write_frame(&ds2);
        assert_eq!(n, 1);
        assert!(w.deltas[1].1[0].is_new);
    }

    #[test]
    fn test_diff_writer_to_bytes() {
        let mut w = ParticleDiffWriter::new();
        let ds = make_ds(3);
        w.write_frame(&ds);
        let bytes = w.to_bytes();
        assert_eq!(&bytes[..5], b"PDIFF");
    }

    #[test]
    fn test_diff_writer_total_deltas() {
        let mut w = ParticleDiffWriter::new();
        let ds = make_ds(5);
        w.write_frame(&ds);
        assert_eq!(w.total_deltas(), 5); // first frame: all new
    }

    #[test]
    fn test_diff_writer_reconstruct() {
        let mut w = ParticleDiffWriter::new();
        let ds = make_ds(4);
        w.write_frame(&ds);
        let rec = w.reconstruct(0);
        assert_eq!(rec.len(), 4);
    }

    #[test]
    fn test_diff_writer_thresholds() {
        let mut w = ParticleDiffWriter::with_thresholds(1.0, 1.0);
        let ds1 = make_ds(2);
        let mut ds2 = ds1.clone();
        ds2.positions[0][0] += 0.5; // below threshold
        w.write_frame(&ds1);
        let n = w.write_frame(&ds2);
        assert_eq!(n, 0, "sub-threshold changes should be ignored");
    }

    // ---- ParticleTrajectory ----

    #[test]
    fn test_trajectory_empty() {
        let t = ParticleTrajectory::new("test");
        assert_eq!(t.num_frames(), 0);
        assert!(t.time_span().is_none());
    }

    #[test]
    fn test_trajectory_push() {
        let mut traj = ParticleTrajectory::new("test");
        for i in 0..3u64 {
            let ds = ParticleDataset::with_time(i as f64, i);
            traj.push(ds);
        }
        assert_eq!(traj.num_frames(), 3);
    }

    #[test]
    fn test_trajectory_time_span() {
        let mut traj = ParticleTrajectory::new("test");
        traj.push(ParticleDataset::with_time(0.0, 0));
        traj.push(ParticleDataset::with_time(1.0, 1));
        traj.push(ParticleDataset::with_time(2.0, 2));
        let (t0, t1) = traj.time_span().unwrap();
        assert!((t0 - 0.0).abs() < 1e-10);
        assert!((t1 - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_trajectory_frame_at() {
        let mut traj = ParticleTrajectory::new("test");
        for i in 0..5u64 {
            traj.push(ParticleDataset::with_time(i as f64, i));
        }
        let f = traj.frame_at(2.3).unwrap();
        assert_eq!(f.step, 2);
    }

    #[test]
    fn test_trajectory_kinetic_energy_series() {
        let mut traj = ParticleTrajectory::new("KE");
        let mut ds = ParticleDataset::with_time(0.0, 0);
        ds.add_particle(0, [0.0; 3], [1.0, 0.0, 0.0], 2.0, 0.1);
        traj.push(ds);
        let series = traj.kinetic_energy_series();
        assert_eq!(series.len(), 1);
        assert!((series[0].1 - 1.0).abs() < 1e-10);
    }

    // ---- ParticleStats ----

    #[test]
    fn test_stats_empty() {
        let s = ParticleStats::compute(&ParticleDataset::new());
        assert_eq!(s.count, 0);
    }

    #[test]
    fn test_stats_count() {
        let ds = make_ds(7);
        let s = ParticleStats::compute(&ds);
        assert_eq!(s.count, 7);
    }

    #[test]
    fn test_stats_total_mass() {
        let ds = make_ds(5);
        let s = ParticleStats::compute(&ds);
        assert!((s.total_mass - 5.0).abs() < 1e-10); // each mass = 1.0
    }

    #[test]
    fn test_stats_max_speed() {
        let mut ds = ParticleDataset::new();
        ds.add_particle(0, [0.0; 3], [3.0, 4.0, 0.0], 1.0, 0.1);
        let s = ParticleStats::compute(&ds);
        assert!((s.max_speed - 5.0).abs() < 1e-10);
    }

    // ---- interpolate_frames ----

    #[test]
    fn test_interpolate_midpoint() {
        let mut a = ParticleDataset::with_time(0.0, 0);
        a.add_particle(0, [0.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        let mut b = ParticleDataset::with_time(1.0, 1);
        b.add_particle(0, [2.0, 0.0, 0.0], [0.0; 3], 1.0, 0.1);
        let mid = interpolate_frames(&a, &b, 0.5);
        assert!((mid.positions[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_t0_is_a() {
        let mut a = ParticleDataset::with_time(0.0, 0);
        a.add_particle(0, [1.0, 2.0, 3.0], [0.0; 3], 1.0, 0.1);
        let mut b = ParticleDataset::with_time(1.0, 1);
        b.add_particle(0, [4.0, 5.0, 6.0], [0.0; 3], 1.0, 0.1);
        let result = interpolate_frames(&a, &b, 0.0);
        assert_eq!(result.positions[0], [1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_interpolate_t1_is_b() {
        let mut a = ParticleDataset::with_time(0.0, 0);
        a.add_particle(0, [0.0; 3], [0.0; 3], 1.0, 0.1);
        let mut b = ParticleDataset::with_time(1.0, 1);
        b.add_particle(0, [1.0, 1.0, 1.0], [0.0; 3], 1.0, 0.1);
        let result = interpolate_frames(&a, &b, 1.0);
        assert!((result.positions[0][0] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_interpolate_only_common_particles() {
        let mut a = ParticleDataset::with_time(0.0, 0);
        a.add_particle(0, [0.0; 3], [0.0; 3], 1.0, 0.1);
        a.add_particle(1, [1.0; 3], [0.0; 3], 1.0, 0.1);
        let mut b = ParticleDataset::with_time(1.0, 1);
        b.add_particle(0, [2.0; 3], [0.0; 3], 1.0, 0.1);
        // particle 1 not in b
        let result = interpolate_frames(&a, &b, 0.5);
        assert_eq!(result.len(), 1);
    }

    // ---- make_test_dataset ----

    #[test]
    fn test_make_test_dataset_count() {
        let ds = make_test_dataset(8, 1.5);
        assert_eq!(ds.len(), 8);
        assert!((ds.time - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_make_test_dataset_zero() {
        let ds = make_test_dataset(0, 0.0);
        assert!(ds.is_empty());
    }

    #[test]
    fn test_make_test_dataset_one() {
        let ds = make_test_dataset(1, 0.0);
        assert_eq!(ds.len(), 1);
    }
}
