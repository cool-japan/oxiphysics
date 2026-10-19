// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Delta encoding of particle frames.
//!
//! [`ParticleDiffWriter`] stores, for every frame, only the particles whose
//! state differs from what a reader reconstructs from the frames written so
//! far.  The encoder tracks the *reconstructed* state, not the raw previous
//! frame, so motion below the configured thresholds is never lost: it keeps
//! accumulating against the last state that was actually written until it
//! exceeds the threshold and is emitted.  The reconstruction error of every
//! position and velocity component therefore never exceeds its threshold.

use super::ParticleDataset;
use std::collections::{BTreeMap, HashSet};

/// Magic bytes that open the binary delta stream.  The last byte is the
/// format version (2: every record carries mass and radius, every frame its
/// simulation time).
const PDIFF_MAGIC: [u8; 8] = *b"PDIFF\0\0\x02";

/// One record of a delta frame.
///
/// For a particle that is [new](Self::is_new) the position and velocity fields
/// hold absolute values; for every other record they hold the change relative
/// to the state a reader has reconstructed from the earlier records.  `mass`
/// and `radius` are always absolute.
#[derive(Debug, Clone)]
pub struct ParticleDelta {
    /// Particle ID.
    pub id: u64,
    /// Position delta `[dx, dy, dz]` (absolute position when `is_new`).
    pub delta_pos: [f64; 3],
    /// Velocity delta `[dvx, dvy, dvz]` (absolute velocity when `is_new`).
    pub delta_vel: [f64; 3],
    /// Absolute particle mass in this frame (the last known mass for a removed
    /// particle).
    pub mass: f64,
    /// Absolute particle radius in this frame (the last known radius for a
    /// removed particle).
    pub radius: f64,
    /// Whether this particle is newly added in this frame.
    pub is_new: bool,
    /// Whether this particle was removed in this frame.
    pub is_removed: bool,
}

/// Particle state as a reader reconstructs it from the records seen so far.
///
/// The encoder and [`ParticleDiffWriter::reconstruct`] both advance this state
/// with [`ReconstructedState::apply`], so the encoder's bookkeeping is
/// bit-identical to what a reader computes.
#[derive(Debug, Clone, Copy)]
struct ReconstructedState {
    position: [f64; 3],
    velocity: [f64; 3],
    mass: f64,
    radius: f64,
}

impl ReconstructedState {
    /// State introduced by a record that is flagged as new.
    fn from_new_record(record: &ParticleDelta) -> Self {
        Self {
            position: record.delta_pos,
            velocity: record.delta_vel,
            mass: record.mass,
            radius: record.radius,
        }
    }

    /// Advance the state by a record of an already known particle.
    fn apply(&mut self, record: &ParticleDelta) {
        for (p, d) in self.position.iter_mut().zip(&record.delta_pos) {
            *p += d;
        }
        for (v, d) in self.velocity.iter_mut().zip(&record.delta_vel) {
            *v += d;
        }
        self.mass = record.mass;
        self.radius = record.radius;
    }
}

/// Component-wise `a - b`.
fn sub3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// Writes only the changes between consecutive particle frames (delta encoding).
///
/// This format is useful for large simulations where only a small fraction of
/// particles move significantly each frame, enabling compact storage.
///
/// A particle is written again when its position or velocity has moved by more
/// than the respective threshold (in any component) since the state a reader
/// reconstructs, or when its mass or radius changed.  Positions and velocities
/// reconstructed with [`reconstruct`](Self::reconstruct) differ from the
/// recorded datasets by at most the thresholds per component, however many
/// frames are written; mass, radius and the frame time are reproduced exactly.
/// Named per-particle properties of [`ParticleDataset`] are not recorded.
#[derive(Debug, Default)]
pub struct ParticleDiffWriter {
    /// State of every live particle exactly as a reader reconstructs it.
    state: BTreeMap<u64, ReconstructedState>,
    /// Threshold below which position changes are considered zero.
    pub pos_threshold: f64,
    /// Threshold below which velocity changes are considered zero.
    pub vel_threshold: f64,
    /// Accumulated delta frames as `(step, simulation time, records)`.
    pub deltas: Vec<(u64, f64, Vec<ParticleDelta>)>,
}

impl ParticleDiffWriter {
    /// Create a diff writer with default thresholds (1e-9 for position, 1e-9 for velocity).
    pub fn new() -> Self {
        Self::with_thresholds(1e-9, 1e-9)
    }

    /// Create with custom thresholds.
    pub fn with_thresholds(pos_threshold: f64, vel_threshold: f64) -> Self {
        Self {
            state: BTreeMap::new(),
            pos_threshold,
            vel_threshold,
            deltas: Vec::new(),
        }
    }

    /// Process a new frame, computing and storing the delta against the state
    /// a reader has reconstructed from the earlier frames.
    ///
    /// The particles of `current` are read from the parallel arrays of the
    /// dataset; if the arrays differ in length the shortest one defines the
    /// particle count.  Returns the number of records written for this frame
    /// (changed, new and removed particles).
    pub fn write_frame(&mut self, current: &ParticleDataset) -> usize {
        let pos_threshold = self.pos_threshold;
        let vel_threshold = self.vel_threshold;
        let mut records: Vec<ParticleDelta> = Vec::new();
        let mut present: HashSet<u64> = HashSet::with_capacity(current.ids.len());

        let particles = current
            .ids
            .iter()
            .zip(&current.positions)
            .zip(&current.velocities)
            .zip(&current.masses)
            .zip(&current.radii);
        for ((((&id, &position), &velocity), &mass), &radius) in particles {
            present.insert(id);
            match self.state.get_mut(&id) {
                None => {
                    let record = ParticleDelta {
                        id,
                        delta_pos: position,
                        delta_vel: velocity,
                        mass,
                        radius,
                        is_new: true,
                        is_removed: false,
                    };
                    self.state
                        .insert(id, ReconstructedState::from_new_record(&record));
                    records.push(record);
                }
                Some(known) => {
                    let delta_pos = sub3(position, known.position);
                    let delta_vel = sub3(velocity, known.velocity);
                    let moved = delta_pos.iter().any(|x| x.abs() > pos_threshold);
                    let accelerated = delta_vel.iter().any(|x| x.abs() > vel_threshold);
                    let reshaped = mass.to_bits() != known.mass.to_bits()
                        || radius.to_bits() != known.radius.to_bits();
                    if moved || accelerated || reshaped {
                        let record = ParticleDelta {
                            id,
                            delta_pos,
                            delta_vel,
                            mass,
                            radius,
                            is_new: false,
                            is_removed: false,
                        };
                        known.apply(&record);
                        records.push(record);
                    }
                }
            }
        }

        // Particles that vanished since the previous frame, in ascending ID order.
        let removed: Vec<u64> = self
            .state
            .keys()
            .filter(|id| !present.contains(*id))
            .copied()
            .collect();
        for id in removed {
            if let Some(last) = self.state.remove(&id) {
                records.push(ParticleDelta {
                    id,
                    delta_pos: [0.0; 3],
                    delta_vel: [0.0; 3],
                    mass: last.mass,
                    radius: last.radius,
                    is_new: false,
                    is_removed: true,
                });
            }
        }

        let written = records.len();
        self.deltas.push((current.step, current.time, records));
        written
    }

    /// Total number of delta records across all frames.
    pub fn total_deltas(&self) -> usize {
        self.deltas.iter().map(|(_, _, d)| d.len()).sum()
    }

    /// Serialise all delta frames to a compact binary format.
    ///
    /// All numbers are little-endian:
    ///
    /// ```text
    /// magic            : b"PDIFF\0\0\x02"       (8 bytes, last byte = version)
    /// num_frames       : u64
    /// per frame:
    ///     step         : u64
    ///     time         : f64
    ///     num_records  : u64
    ///     per record:
    ///         id           : u64
    ///         flags        : u8   (bit 0 = new, bit 1 = removed)
    ///         delta_pos    : 3 x f64
    ///         delta_vel    : 3 x f64
    ///         mass         : f64
    ///         radius       : f64
    /// ```
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = PDIFF_MAGIC.to_vec();
        out.extend_from_slice(&(self.deltas.len() as u64).to_le_bytes());
        for (step, time, records) in &self.deltas {
            out.extend_from_slice(&step.to_le_bytes());
            out.extend_from_slice(&time.to_le_bytes());
            out.extend_from_slice(&(records.len() as u64).to_le_bytes());
            for d in records {
                out.extend_from_slice(&d.id.to_le_bytes());
                let flags: u8 = u8::from(d.is_new) | (u8::from(d.is_removed) << 1);
                out.push(flags);
                for x in d.delta_pos.iter().chain(&d.delta_vel) {
                    out.extend_from_slice(&x.to_le_bytes());
                }
                out.extend_from_slice(&d.mass.to_le_bytes());
                out.extend_from_slice(&d.radius.to_le_bytes());
            }
        }
        out
    }

    /// Reconstruct the full dataset after the last recorded frame whose step is
    /// not greater than `through_step`.
    ///
    /// Frames are applied in the order they were written and application stops
    /// at the first frame beyond `through_step`.  Positions and velocities are
    /// the accumulated deltas; mass and radius are taken from the latest record
    /// of each particle; the dataset's `time` and `step` are those of the last
    /// frame applied.  Particles are sorted by ID.  Records that refer to a
    /// particle that was never introduced are ignored.  When no frame lies at
    /// or before `through_step` an empty dataset at time 0 and step
    /// `through_step` is returned.
    pub fn reconstruct(&self, through_step: u64) -> ParticleDataset {
        let mut state: BTreeMap<u64, ReconstructedState> = BTreeMap::new();
        let mut applied: Option<(u64, f64)> = None;

        for (step, time, records) in &self.deltas {
            if *step > through_step {
                break;
            }
            for record in records {
                if record.is_removed {
                    state.remove(&record.id);
                } else if record.is_new {
                    state.insert(record.id, ReconstructedState::from_new_record(record));
                } else if let Some(known) = state.get_mut(&record.id) {
                    known.apply(record);
                }
            }
            applied = Some((*step, *time));
        }

        let (step, time) = applied.unwrap_or((through_step, 0.0));
        let mut dataset = ParticleDataset::with_time(time, step);
        for (id, particle) in state {
            dataset.add_particle(
                id,
                particle.position,
                particle.velocity,
                particle.mass,
                particle.radius,
            );
        }
        dataset
    }
}
