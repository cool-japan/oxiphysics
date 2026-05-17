// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! Simulation state export/import types and utilities.

use serde::{Deserialize, Serialize};
use wasm_bindgen::prelude::*;

use crate::wasm_helpers::err_to_jsvalue;

// ---------------------------------------------------------------------------
// ExportFormat
// ---------------------------------------------------------------------------

/// Export format selection.
///
/// Exposed to JavaScript as `SimExportFormat` to avoid colliding with
/// [`crate::io_bridge::ExportFormat`].
#[wasm_bindgen(js_name = "SimExportFormat")]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ExportFormat {
    /// JSON text format.
    Json,
    /// Binary (CBOR-like) format.
    Binary,
    /// VTK legacy ASCII format for ParaView.
    Vtk,
    /// CSV comma-separated values.
    Csv,
}

// ---------------------------------------------------------------------------
// SimulationSnapshot
// ---------------------------------------------------------------------------

/// Serialised snapshot of one simulation frame for export.
///
/// Exposed to JavaScript as `SimExportSnapshot` to avoid colliding with the
/// snapshot type defined in [`crate::body_query`].
#[wasm_bindgen(js_name = "SimExportSnapshot")]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationSnapshot {
    /// Simulation time (s).
    pub time: f64,
    /// Flat body positions \[x,y,z, ...\].
    #[wasm_bindgen(skip)]
    pub positions: Vec<f64>,
    /// Flat body velocities \[vx,vy,vz, ...\].
    #[wasm_bindgen(skip)]
    pub velocities: Vec<f64>,
    /// Flat body quaternions \[qx,qy,qz,qw, ...\].
    #[wasm_bindgen(skip)]
    pub orientations: Vec<f64>,
    /// Flat angular velocities \[wx,wy,wz, ...\].
    #[wasm_bindgen(skip)]
    pub angular_velocities: Vec<f64>,
    /// Body masses.
    #[wasm_bindgen(skip)]
    pub masses: Vec<f64>,
}

impl SimulationSnapshot {
    /// Create an empty snapshot at the given time.
    pub fn new(time: f64) -> Self {
        SimulationSnapshot {
            time,
            positions: Vec::new(),
            velocities: Vec::new(),
            orientations: Vec::new(),
            angular_velocities: Vec::new(),
            masses: Vec::new(),
        }
    }

    /// Number of bodies in the snapshot.
    pub fn body_count(&self) -> usize {
        self.masses.len()
    }
}

#[wasm_bindgen(js_class = "SimExportSnapshot")]
impl SimulationSnapshot {
    /// Construct an empty snapshot at simulation time `time`.
    #[wasm_bindgen(constructor)]
    pub fn new_js(time: f64) -> SimulationSnapshot {
        SimulationSnapshot::new(time)
    }

    /// Number of bodies (derived from the mass array length).
    #[wasm_bindgen(js_name = "body_count")]
    pub fn body_count_js(&self) -> u32 {
        self.body_count() as u32
    }

    /// Return a clone of the position array.
    #[wasm_bindgen(js_name = "positions")]
    pub fn positions_js(&self) -> Vec<f64> {
        self.positions.clone()
    }
    /// Return a clone of the velocity array.
    #[wasm_bindgen(js_name = "velocities")]
    pub fn velocities_js(&self) -> Vec<f64> {
        self.velocities.clone()
    }
    /// Return a clone of the orientation array.
    #[wasm_bindgen(js_name = "orientations")]
    pub fn orientations_js(&self) -> Vec<f64> {
        self.orientations.clone()
    }
    /// Return a clone of the angular velocity array.
    #[wasm_bindgen(js_name = "angular_velocities")]
    pub fn angular_velocities_js(&self) -> Vec<f64> {
        self.angular_velocities.clone()
    }
    /// Return a clone of the mass array.
    #[wasm_bindgen(js_name = "masses")]
    pub fn masses_js(&self) -> Vec<f64> {
        self.masses.clone()
    }
    /// Replace the position array with the supplied flat data.
    #[wasm_bindgen(js_name = "set_positions")]
    pub fn set_positions_js(&mut self, positions: Vec<f64>) {
        self.positions = positions;
    }
    /// Replace the velocity array with the supplied flat data.
    #[wasm_bindgen(js_name = "set_velocities")]
    pub fn set_velocities_js(&mut self, velocities: Vec<f64>) {
        self.velocities = velocities;
    }
    /// Replace the orientation array with the supplied flat data.
    #[wasm_bindgen(js_name = "set_orientations")]
    pub fn set_orientations_js(&mut self, orientations: Vec<f64>) {
        self.orientations = orientations;
    }
    /// Replace the angular-velocity array with the supplied flat data.
    #[wasm_bindgen(js_name = "set_angular_velocities")]
    pub fn set_angular_velocities_js(&mut self, angular_velocities: Vec<f64>) {
        self.angular_velocities = angular_velocities;
    }
    /// Replace the mass array with the supplied data.
    #[wasm_bindgen(js_name = "set_masses")]
    pub fn set_masses_js(&mut self, masses: Vec<f64>) {
        self.masses = masses;
    }
}

// ---------------------------------------------------------------------------
// SimulationExporter
// ---------------------------------------------------------------------------

/// Utility for exporting and importing simulation state.
#[wasm_bindgen]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SimulationExporter {
    /// Cached snapshots (e.g., for VTK series).
    #[wasm_bindgen(skip)]
    pub snapshots: Vec<SimulationSnapshot>,
}

impl SimulationExporter {
    /// Create a new empty exporter.
    pub fn new() -> Self {
        Self::default()
    }

    /// Export a single snapshot to JSON.
    pub fn export_json(&self, snapshot: &SimulationSnapshot) -> String {
        serde_json::to_string_pretty(snapshot).unwrap_or_default()
    }

    /// Export a single snapshot to a simple binary format (little-endian f64 stream).
    pub fn export_binary(&self, snapshot: &SimulationSnapshot) -> Vec<u8> {
        let mut buf = Vec::new();
        // Header: time as 8 bytes
        buf.extend_from_slice(&snapshot.time.to_le_bytes());
        let n = snapshot.positions.len() as u64;
        buf.extend_from_slice(&n.to_le_bytes());
        for &v in &snapshot.positions {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        for &v in &snapshot.velocities {
            buf.extend_from_slice(&v.to_le_bytes());
        }
        buf
    }

    /// Export a snapshot in VTK legacy ASCII format.
    pub fn export_vtk(&self, snapshot: &SimulationSnapshot, frame_idx: usize) -> String {
        let n = snapshot.positions.len() / 3;
        let mut out = format!(
            "# vtk DataFile Version 3.0\nFrame {frame_idx} t={:.6}\nASCII\nDATASET UNSTRUCTURED_GRID\n",
            snapshot.time
        );
        out.push_str(&format!("POINTS {n} double\n"));
        for i in 0..n {
            let x = snapshot.positions.get(i * 3).copied().unwrap_or(0.0);
            let y = snapshot.positions.get(i * 3 + 1).copied().unwrap_or(0.0);
            let z = snapshot.positions.get(i * 3 + 2).copied().unwrap_or(0.0);
            out.push_str(&format!("{x:.6} {y:.6} {z:.6}\n"));
        }
        out
    }

    /// Export all cached snapshots as a VTK series (one string per frame).
    pub fn export_vtk_series(&self) -> Vec<String> {
        self.snapshots
            .iter()
            .enumerate()
            .map(|(i, s)| self.export_vtk(s, i))
            .collect()
    }

    /// Import a snapshot from JSON.
    pub fn import_state(&self, json: &str) -> Result<SimulationSnapshot, String> {
        serde_json::from_str(json).map_err(|e| e.to_string())
    }

    /// Add a snapshot to the internal series.
    pub fn push_snapshot(&mut self, snapshot: SimulationSnapshot) {
        self.snapshots.push(snapshot);
    }

    /// Clear cached snapshots.
    pub fn clear(&mut self) {
        self.snapshots.clear();
    }
}

#[wasm_bindgen]
impl SimulationExporter {
    /// Construct a new empty exporter.
    #[wasm_bindgen(constructor)]
    pub fn new_js() -> SimulationExporter {
        SimulationExporter::default()
    }

    /// Number of cached snapshots.
    #[wasm_bindgen(js_name = "snapshot_count")]
    pub fn snapshot_count_js(&self) -> u32 {
        self.snapshots.len() as u32
    }

    /// Export a single snapshot to a (pretty) JSON string.
    #[wasm_bindgen(js_name = "export_json")]
    pub fn export_json_js(&self, snapshot: &SimulationSnapshot) -> String {
        self.export_json(snapshot)
    }

    /// Export a single snapshot to a packed little-endian binary buffer.
    #[wasm_bindgen(js_name = "export_binary")]
    pub fn export_binary_js(&self, snapshot: &SimulationSnapshot) -> Vec<u8> {
        self.export_binary(snapshot)
    }

    /// Export a single snapshot to a VTK legacy ASCII string.
    #[wasm_bindgen(js_name = "export_vtk")]
    pub fn export_vtk_js(&self, snapshot: &SimulationSnapshot, frame_idx: u32) -> String {
        self.export_vtk(snapshot, frame_idx as usize)
    }

    /// Export all cached snapshots as a VTK series; returns a JSON-encoded
    /// array of strings (one per frame) for transfer to JavaScript.
    #[wasm_bindgen(js_name = "export_vtk_series_json")]
    pub fn export_vtk_series_json(&self) -> Result<String, JsValue> {
        let series = self.export_vtk_series();
        serde_json::to_string(&series).map_err(err_to_jsvalue)
    }

    /// Import a snapshot from a JSON string.
    #[wasm_bindgen(js_name = "import_state")]
    pub fn import_state_js(&self, json: &str) -> Result<SimulationSnapshot, JsValue> {
        self.import_state(json).map_err(err_to_jsvalue)
    }

    /// Add a snapshot to the internal series (consumes a clone of `snapshot`).
    #[wasm_bindgen(js_name = "push_snapshot")]
    pub fn push_snapshot_js(&mut self, snapshot: &SimulationSnapshot) {
        self.push_snapshot(snapshot.clone());
    }

    /// Clear all cached snapshots.
    #[wasm_bindgen(js_name = "clear")]
    pub fn clear_js(&mut self) {
        self.clear();
    }
}
