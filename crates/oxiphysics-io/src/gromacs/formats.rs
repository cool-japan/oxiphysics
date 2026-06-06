// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GROMACS auxiliary file formats: XTC trajectory, XVG data, energy records,
//! ITP atom types, residue library, and analysis helpers.

// ============================================================================
// XTC Trajectory
// ============================================================================

/// A minimal XTC-like compressed trajectory frame (positions only, no C-library).
#[derive(Debug, Clone, Default)]
pub struct XtcTrajectory {
    /// Number of atoms (must be constant across all frames).
    pub n_atoms: usize,
    /// Frames: each frame is a flat list of (x, y, z) in nm stored as f32.
    pub frames: Vec<Vec<f32>>,
    /// Simulation time for each frame (ps).
    pub times: Vec<f32>,
}

impl XtcTrajectory {
    /// Create a new trajectory for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            n_atoms,
            frames: Vec::new(),
            times: Vec::new(),
        }
    }

    /// Add a frame of positions (flat \[x,y,z,...\] per atom).
    pub fn add_frame(&mut self, positions: &[[f32; 3]], time_ps: f32) -> Result<(), String> {
        if positions.len() != self.n_atoms {
            return Err(format!(
                "expected {} atoms, got {}",
                self.n_atoms,
                positions.len()
            ));
        }
        let flat: Vec<f32> = positions.iter().flat_map(|p| p.iter().copied()).collect();
        self.frames.push(flat);
        self.times.push(time_ps);
        Ok(())
    }

    /// Get position of atom `atom_idx` in frame `frame_idx`.
    pub fn get_position(&self, frame_idx: usize, atom_idx: usize) -> Option<[f32; 3]> {
        let frame = self.frames.get(frame_idx)?;
        let base = atom_idx * 3;
        if base + 2 >= frame.len() {
            return None;
        }
        Some([frame[base], frame[base + 1], frame[base + 2]])
    }

    /// Number of frames.
    pub fn n_frames(&self) -> usize {
        self.frames.len()
    }
}

// ============================================================================
// XTC Writer
// ============================================================================

/// An extended XTC trajectory writer that produces compact binary frames.
pub struct XtcWriter {
    /// Number of atoms (constant across frames).
    pub n_atoms: usize,
    /// Serialised frame data.
    pub data: Vec<u8>,
    /// Number of frames written.
    pub n_frames: u32,
}

impl XtcWriter {
    /// Create a new XTC writer for `n_atoms` atoms.
    pub fn new(n_atoms: usize) -> Self {
        Self {
            n_atoms,
            data: Vec::new(),
            n_frames: 0,
        }
    }

    /// Append a compressed frame to the internal buffer.
    pub fn write_compressed_frame(
        &mut self,
        positions: &[[f32; 3]],
        time_ps: f32,
    ) -> std::result::Result<(), String> {
        if positions.len() != self.n_atoms {
            return Err(format!(
                "expected {} atoms, got {}",
                self.n_atoms,
                positions.len()
            ));
        }
        self.data.extend_from_slice(&self.n_frames.to_le_bytes());
        self.data.extend_from_slice(&time_ps.to_le_bytes());
        self.data
            .extend_from_slice(&(self.n_atoms as u32).to_le_bytes());
        for p in positions {
            for &v in p.iter() {
                self.data.extend_from_slice(&v.to_le_bytes());
            }
        }
        self.n_frames += 1;
        Ok(())
    }

    /// Read back a frame from the internal buffer by frame index.
    pub fn read_frame(
        &self,
        frame_idx: usize,
    ) -> std::result::Result<(f32, Vec<[f32; 3]>), String> {
        let frame_bytes = 4 + 4 + 4 + self.n_atoms * 3 * 4;
        let offset = frame_idx * frame_bytes;
        if offset + frame_bytes > self.data.len() {
            return Err(format!("frame {} out of range", frame_idx));
        }
        let d = &self.data[offset..];
        let time_ps = f32::from_le_bytes([d[4], d[5], d[6], d[7]]);
        let mut positions = Vec::with_capacity(self.n_atoms);
        let coord_start = 12;
        for i in 0..self.n_atoms {
            let base = coord_start + i * 12;
            let x = f32::from_le_bytes([d[base], d[base + 1], d[base + 2], d[base + 3]]);
            let y = f32::from_le_bytes([d[base + 4], d[base + 5], d[base + 6], d[base + 7]]);
            let z = f32::from_le_bytes([d[base + 8], d[base + 9], d[base + 10], d[base + 11]]);
            positions.push([x, y, z]);
        }
        Ok((time_ps, positions))
    }

    /// Total bytes written to the internal buffer.
    pub fn byte_len(&self) -> usize {
        self.data.len()
    }
}

// ============================================================================
// Energy Records
// ============================================================================

/// A single energy term record.
#[derive(Debug, Clone)]
pub struct EnergyRecord {
    /// Name of the energy term.
    pub name: String,
    /// Value (kJ/mol).
    pub value: f64,
    /// Standard deviation (kJ/mol), if available.
    pub std_dev: Option<f64>,
}

/// A simple GROMACS energy file entry (one time step).
#[derive(Debug, Clone, Default)]
pub struct EnergyEntry {
    /// Simulation time (ps).
    pub time: f64,
    /// Energy terms.
    pub terms: Vec<EnergyRecord>,
}

/// Parse a GROMACS-like energy summary in plain-text format.
pub fn parse_energy_summary(s: &str) -> Vec<EnergyEntry> {
    let mut entries: Vec<EnergyEntry> = Vec::new();
    let mut current: Option<EnergyEntry> = None;

    for line in s.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with(';') {
            continue;
        }
        if let Some(t) = trimmed.strip_prefix("Time:") {
            if let Some(entry) = current.take() {
                entries.push(entry);
            }
            let time: f64 = t.trim().parse().unwrap_or(0.0);
            current = Some(EnergyEntry {
                time,
                terms: Vec::new(),
            });
        } else if let Some(ref mut entry) = current {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.is_empty() {
                continue;
            }
            let name = parts[0].to_string();
            let value: f64 = parts.get(1).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            let std_dev: Option<f64> = parts.get(2).and_then(|s| s.parse().ok());
            entry.terms.push(EnergyRecord {
                name,
                value,
                std_dev,
            });
        }
    }
    if let Some(entry) = current {
        entries.push(entry);
    }
    entries
}

// ============================================================================
// XVG File
// ============================================================================

/// A parsed GROMACS XVG (Grace/xmgrace) data file.
#[derive(Debug, Clone, Default)]
pub struct XvgFile {
    /// Title from `@ title` directive.
    pub title: String,
    /// X-axis label.
    pub xaxis_label: String,
    /// Y-axis label.
    pub yaxis_label: String,
    /// Legend entries.
    pub legend: Vec<String>,
    /// Data columns: each row is `[x, y0, y1, ...]`.
    pub data: Vec<Vec<f64>>,
}

impl XvgFile {
    /// Parse an XVG file from a string.
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut xvg = XvgFile::default();

        for line in s.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if let Some(after_at) = trimmed.strip_prefix('@') {
                let directive = &after_at.trim_start();
                if let Some(rest) = directive.strip_prefix("title") {
                    xvg.title = rest.trim().trim_matches('"').to_string();
                } else if directive.contains("xaxis") && directive.contains("label") {
                    if let Some(pos) = directive.find('"') {
                        let end = directive[pos + 1..].find('"').map(|e| pos + 1 + e);
                        if let Some(end_pos) = end {
                            xvg.xaxis_label = directive[pos + 1..end_pos].to_string();
                        }
                    }
                } else if directive.contains("yaxis") && directive.contains("label") {
                    if let Some(pos) = directive.find('"') {
                        let end = directive[pos + 1..].find('"').map(|e| pos + 1 + e);
                        if let Some(end_pos) = end {
                            xvg.yaxis_label = directive[pos + 1..end_pos].to_string();
                        }
                    }
                } else if let Some(rest) = directive.strip_prefix("legend") {
                    let entry = rest.trim().trim_matches('"').to_string();
                    if !entry.is_empty() {
                        xvg.legend.push(entry);
                    }
                }
                continue;
            }
            let values: Vec<f64> = trimmed
                .split_whitespace()
                .filter_map(|tok| tok.parse::<f64>().ok())
                .collect();
            if !values.is_empty() {
                xvg.data.push(values);
            }
        }
        Ok(xvg)
    }

    /// Number of data rows.
    pub fn n_rows(&self) -> usize {
        self.data.len()
    }
    /// Number of columns.
    pub fn n_cols(&self) -> usize {
        self.data.first().map(|r| r.len()).unwrap_or(0)
    }
    /// Extract the X column (column 0).
    pub fn x_values(&self) -> Vec<f64> {
        self.data
            .iter()
            .filter_map(|row| row.first().copied())
            .collect()
    }
    /// Extract column `col`.
    pub fn column(&self, col: usize) -> Vec<f64> {
        self.data
            .iter()
            .filter_map(|row| row.get(col).copied())
            .collect()
    }
    /// Compute mean of column `col`.
    pub fn column_mean(&self, col: usize) -> Option<f64> {
        let vals = self.column(col);
        if vals.is_empty() {
            return None;
        }
        Some(vals.iter().sum::<f64>() / vals.len() as f64)
    }
    /// Find the minimum value in column `col`.
    pub fn column_min(&self, col: usize) -> Option<f64> {
        self.column(col).into_iter().reduce(f64::min)
    }
    /// Find the maximum value in column `col`.
    pub fn column_max(&self, col: usize) -> Option<f64> {
        self.column(col).into_iter().reduce(f64::max)
    }

    /// Write the XVG data back to a string in Grace format.
    pub fn to_xvg_string(&self) -> String {
        let mut s = String::new();
        if !self.title.is_empty() {
            s.push_str(&format!("@ title \"{}\"\n", self.title));
        }
        if !self.xaxis_label.is_empty() {
            s.push_str(&format!("@ xaxis label \"{}\"\n", self.xaxis_label));
        }
        if !self.yaxis_label.is_empty() {
            s.push_str(&format!("@ yaxis label \"{}\"\n", self.yaxis_label));
        }
        for legend in &self.legend {
            s.push_str(&format!("@ legend \"{}\"\n", legend));
        }
        for row in &self.data {
            let vals: Vec<String> = row.iter().map(|v| format!("{:.6}", v)).collect();
            s.push_str(&vals.join("  "));
            s.push('\n');
        }
        s
    }
}

impl std::str::FromStr for XvgFile {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ============================================================================
// ITP Atom Types
// ============================================================================

/// An atom type definition for GROMACS.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomTypeDef {
    /// Atom type name.
    pub name: String,
    /// Atomic mass (amu).
    pub mass: f64,
    /// Charge (e).
    pub charge: f64,
    /// Lennard-Jones sigma (nm).
    pub sigma: f64,
    /// Lennard-Jones epsilon (kJ/mol).
    pub epsilon: f64,
}

/// A GROMACS force-field include file (`.itp`) atom type section entry.
#[derive(Debug, Clone, Default)]
pub struct ItpAtomTypes {
    /// Parsed atom type definitions.
    pub atom_types: Vec<AtomTypeDef>,
}

impl ItpAtomTypes {
    /// Parse atom type definitions from a GROMACS `[ atomtypes ]` section string.
    pub fn parse(s: &str) -> Result<Self, String> {
        let mut atom_types = Vec::new();
        let mut in_section = false;

        for line in s.lines() {
            let trimmed = line.trim();
            if trimmed.is_empty() || trimmed.starts_with(';') {
                continue;
            }
            if let Some(inner) = trimmed.strip_prefix('[').and_then(|s| s.strip_suffix(']')) {
                let sec_name = inner.trim().to_lowercase();
                in_section = sec_name == "atomtypes";
                continue;
            }
            if !in_section {
                continue;
            }
            let line_data = if let Some(sc) = trimmed.find(';') {
                trimmed[..sc].trim()
            } else {
                trimmed
            };
            let parts: Vec<&str> = line_data.split_whitespace().collect();
            if parts.len() < 6 {
                continue;
            }
            let name = parts[0].to_string();
            let mass: f64 = parts[1].parse().unwrap_or(0.0);
            let charge: f64 = parts[2].parse().unwrap_or(0.0);
            let sigma: f64 = parts[4].parse().unwrap_or(0.0);
            let epsilon: f64 = parts[5].parse().unwrap_or(0.0);
            atom_types.push(AtomTypeDef {
                name,
                mass,
                charge,
                sigma,
                epsilon,
            });
        }
        Ok(ItpAtomTypes { atom_types })
    }

    /// Find an atom type by name.
    pub fn get(&self, name: &str) -> Option<&AtomTypeDef> {
        self.atom_types.iter().find(|a| a.name == name)
    }
    /// Number of atom types.
    pub fn len(&self) -> usize {
        self.atom_types.len()
    }
    /// Whether the table is empty.
    pub fn is_empty(&self) -> bool {
        self.atom_types.is_empty()
    }
}

impl std::str::FromStr for ItpAtomTypes {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ============================================================================
// Residue Library
// ============================================================================

/// A minimal residue template entry.
#[derive(Debug, Clone, PartialEq)]
pub struct ResidueAtomTemplate {
    /// Atom name.
    pub atom_name: String,
    /// Atom type string.
    pub atom_type: String,
    /// Partial charge (e).
    pub charge: f64,
    /// Mass (amu).
    pub mass: f64,
}

/// A residue template.
#[derive(Debug, Clone)]
pub struct ResidueTemplate {
    /// Three-letter residue code.
    pub name: String,
    /// Atom definitions for this residue.
    pub atoms: Vec<ResidueAtomTemplate>,
}

impl ResidueTemplate {
    /// Create a new residue template.
    pub fn new(name: &str) -> Self {
        ResidueTemplate {
            name: name.to_string(),
            atoms: Vec::new(),
        }
    }
    /// Add an atom to this residue.
    pub fn add_atom(&mut self, atom_name: &str, atom_type: &str, charge: f64, mass: f64) {
        self.atoms.push(ResidueAtomTemplate {
            atom_name: atom_name.to_string(),
            atom_type: atom_type.to_string(),
            charge,
            mass,
        });
    }
    /// Look up an atom by name.
    pub fn get_atom(&self, name: &str) -> Option<&ResidueAtomTemplate> {
        self.atoms.iter().find(|a| a.atom_name == name)
    }
    /// Total charge.
    pub fn total_charge(&self) -> f64 {
        self.atoms.iter().map(|a| a.charge).sum()
    }
    /// Total mass.
    pub fn total_mass(&self) -> f64 {
        self.atoms.iter().map(|a| a.mass).sum()
    }
    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.atoms.len()
    }
}

/// A library of residue templates.
#[derive(Debug, Clone, Default)]
pub struct ResidueLibrary {
    /// Residue templates indexed by name.
    pub templates: Vec<ResidueTemplate>,
}

impl ResidueLibrary {
    /// Create an empty residue library.
    pub fn new() -> Self {
        Self::default()
    }
    /// Add a residue template.
    pub fn add(&mut self, template: ResidueTemplate) {
        self.templates.push(template);
    }
    /// Look up a residue by name.
    pub fn get(&self, name: &str) -> Option<&ResidueTemplate> {
        self.templates.iter().find(|t| t.name == name)
    }
    /// Look up a residue by name (case-insensitive).
    pub fn get_ci(&self, name: &str) -> Option<&ResidueTemplate> {
        let upper = name.to_uppercase();
        self.templates
            .iter()
            .find(|t| t.name.to_uppercase() == upper)
    }
    /// Number of residue templates.
    pub fn len(&self) -> usize {
        self.templates.len()
    }
    /// Whether the library is empty.
    pub fn is_empty(&self) -> bool {
        self.templates.is_empty()
    }

    /// Build a small GROMACS/AMBER-compatible residue library with backbone atoms.
    pub fn amber_backbone() -> Self {
        let mut lib = ResidueLibrary::new();
        let mut ala = ResidueTemplate::new("ALA");
        ala.add_atom("N", "N", -0.4157, 14.007);
        ala.add_atom("H", "H", 0.2719, 1.008);
        ala.add_atom("CA", "CT", 0.0337, 12.011);
        ala.add_atom("HA", "H1", 0.0823, 1.008);
        ala.add_atom("CB", "CT", -0.1825, 12.011);
        ala.add_atom("C", "C", 0.5973, 12.011);
        ala.add_atom("O", "O", -0.5679, 15.999);
        lib.add(ala);
        let mut gly = ResidueTemplate::new("GLY");
        gly.add_atom("N", "N", -0.4157, 14.007);
        gly.add_atom("H", "H", 0.2719, 1.008);
        gly.add_atom("CA", "CT", 0.0, 12.011);
        gly.add_atom("HA2", "H1", 0.0698, 1.008);
        gly.add_atom("C", "C", 0.5973, 12.011);
        gly.add_atom("O", "O", -0.5679, 15.999);
        lib.add(gly);
        let mut wat = ResidueTemplate::new("HOH");
        wat.add_atom("O", "OW", -0.8476, 15.999);
        wat.add_atom("H1", "HW", 0.4238, 1.008);
        wat.add_atom("H2", "HW", 0.4238, 1.008);
        lib.add(wat);
        lib
    }
}

// ============================================================================
// Analysis Helpers
// ============================================================================

/// Compute the running average of a data column.
pub fn running_average(data: &[f64]) -> Vec<f64> {
    let mut result = Vec::with_capacity(data.len());
    let mut sum = 0.0;
    for (i, &v) in data.iter().enumerate() {
        sum += v;
        result.push(sum / (i + 1) as f64);
    }
    result
}

/// Compute the block average of `data` using blocks of size `block_size`.
pub fn block_average(data: &[f64], block_size: usize) -> Vec<f64> {
    if block_size == 0 || data.is_empty() {
        return Vec::new();
    }
    data.chunks(block_size)
        .map(|chunk| chunk.iter().sum::<f64>() / chunk.len() as f64)
        .collect()
}

/// Compute the autocorrelation function of `data` at lag `lag`.
pub fn autocorrelation(data: &[f64], lag: usize) -> f64 {
    let n = data.len();
    if n <= lag {
        return 0.0;
    }
    let mean = data.iter().sum::<f64>() / n as f64;
    let var: f64 = data.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    if var < 1e-30 {
        return 0.0;
    }
    let cov: f64 = (0..n - lag)
        .map(|i| (data[i] - mean) * (data[i + lag] - mean))
        .sum::<f64>()
        / (n - lag) as f64;
    cov / var
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests_formats {
    use super::*;

    #[test]
    fn test_xvg_parse_basic() {
        let xvg_str = "# comment\n@ title \"Total Energy\"\n@ xaxis label \"Time (ps)\"\n@ yaxis label \"E (kJ/mol)\"\n0.000  -1234.56  567.89\n0.002  -1230.11  568.22\n";
        let xvg = XvgFile::parse(xvg_str).expect("parse");
        assert_eq!(xvg.n_rows(), 2);
        assert_eq!(xvg.title, "Total Energy");
    }

    #[test]
    fn test_xvg_column_mean() {
        let xvg_str = "0.0  10.0\n0.1  20.0\n0.2  30.0\n";
        let xvg = XvgFile::parse(xvg_str).expect("parse");
        let mean = xvg.column_mean(1).expect("mean");
        assert!((mean - 20.0).abs() < 1e-9);
    }

    #[test]
    fn test_energy_summary_parse() {
        let data = "Time: 0.0\nLJ-14  100.0  2.5\nCoulomb  -200.0\n";
        let entries = parse_energy_summary(data);
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].terms.len(), 2);
    }

    #[test]
    fn test_xtc_trajectory_add_frame() {
        let mut traj = XtcTrajectory::new(2);
        let positions = [[0.0f32, 0.0, 0.0], [1.0, 2.0, 3.0]];
        traj.add_frame(&positions, 0.0).expect("add_frame");
        assert_eq!(traj.n_frames(), 1);
        let p = traj.get_position(0, 1).expect("get_position");
        assert!((p[0] - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_xtc_writer_roundtrip() {
        let mut xtc = XtcWriter::new(3);
        let positions = [[1.0f32, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        xtc.write_compressed_frame(&positions, 0.1).expect("write");
        let (t, pos) = xtc.read_frame(0).expect("read");
        assert!((t - 0.1f32).abs() < 1e-5);
        assert!((pos[0][0] - 1.0f32).abs() < 1e-5);
    }

    #[test]
    fn test_itp_atom_types_parse() {
        let itp = "[ atomtypes ]\nCT   12.011  0.000  A  3.39967e-01  4.57730e-01\n";
        let types = ItpAtomTypes::parse(itp).expect("parse");
        assert_eq!(types.len(), 1);
    }

    #[test]
    fn test_residue_library_amber() {
        let lib = ResidueLibrary::amber_backbone();
        assert!(lib.get("ALA").is_some());
        assert!(lib.get("GLY").is_some());
    }

    #[test]
    fn test_running_average() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ra = running_average(&data);
        assert!((ra[0] - 1.0).abs() < 1e-10);
        assert!((ra[4] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_block_average() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let ba = block_average(&data, 2);
        assert_eq!(ba.len(), 3);
        assert!((ba[0] - 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_autocorrelation_zero_lag() {
        let data = [1.0, 2.0, 3.0, 4.0, 5.0];
        let ac = autocorrelation(&data, 0);
        assert!((ac - 1.0).abs() < 1e-10);
    }
}
