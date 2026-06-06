// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! GROMACS GRO format reader, writer, trajectory, and utility functions.

use std::io::{BufRead, Write};

/// A single atom record in a GRO file.
#[derive(Debug, Clone, PartialEq)]
pub struct GroAtom {
    /// Residue (molecule) number.
    pub residue_id: i32,
    /// Residue name (max 5 chars).
    pub residue_name: String,
    /// Atom name (max 5 chars).
    pub atom_name: String,
    /// Atom number.
    pub atom_id: i32,
    /// Position in nm.
    pub position: [f64; 3],
    /// Velocity in nm/ps (optional).
    pub velocity: Option<[f64; 3]>,
}

/// Parsed representation of a GROMACS GRO file.
#[derive(Debug, Clone)]
pub struct GroFile {
    /// Title / comment line.
    pub title: String,
    /// Atom records.
    pub atoms: Vec<GroAtom>,
    /// Orthorhombic box vectors in nm.
    pub box_vectors: [f64; 3],
}

impl GroFile {
    /// Parse a GRO file from any `BufRead` source.
    pub fn read(reader: impl BufRead) -> Result<Self, String> {
        let mut lines = reader.lines();

        // Line 1: title
        let title = lines
            .next()
            .ok_or_else(|| "GRO: missing title line".to_string())?
            .map_err(|e| e.to_string())?;

        // Line 2: number of atoms
        let n_atoms: usize = lines
            .next()
            .ok_or_else(|| "GRO: missing atom count line".to_string())?
            .map_err(|e| e.to_string())?
            .trim()
            .parse()
            .map_err(|e: std::num::ParseIntError| format!("GRO: bad atom count: {e}"))?;

        let mut atoms = Vec::with_capacity(n_atoms);

        for i in 0..n_atoms {
            let raw = lines
                .next()
                .ok_or_else(|| format!("GRO: missing atom line {i}"))?
                .map_err(|e| e.to_string())?;

            // Fixed-width columns per GRO spec:
            //  0..5   residue id
            //  5..10  residue name
            // 10..15  atom name
            // 15..20  atom id
            // 20..28  x
            // 28..36  y
            // 36..44  z
            // 44..52  vx (optional)
            // 52..60  vy (optional)
            // 60..68  vz (optional)
            if raw.len() < 44 {
                return Err(format!(
                    "GRO: atom line {i} too short ({} chars)",
                    raw.len()
                ));
            }
            let residue_id: i32 =
                raw[0..5]
                    .trim()
                    .parse()
                    .map_err(|e: std::num::ParseIntError| {
                        format!("GRO: bad residue_id on line {i}: {e}")
                    })?;
            let residue_name = raw[5..10].trim().to_string();
            let atom_name = raw[10..15].trim().to_string();
            let atom_id: i32 =
                raw[15..20]
                    .trim()
                    .parse()
                    .map_err(|e: std::num::ParseIntError| {
                        format!("GRO: bad atom_id on line {i}: {e}")
                    })?;
            let x: f64 = raw[20..28]
                .trim()
                .parse()
                .map_err(|e: std::num::ParseFloatError| format!("GRO: bad x on line {i}: {e}"))?;
            let y: f64 = raw[28..36]
                .trim()
                .parse()
                .map_err(|e: std::num::ParseFloatError| format!("GRO: bad y on line {i}: {e}"))?;
            let z: f64 = raw[36..44]
                .trim()
                .parse()
                .map_err(|e: std::num::ParseFloatError| format!("GRO: bad z on line {i}: {e}"))?;

            let velocity = if raw.len() >= 68 {
                let vx: f64 =
                    raw[44..52]
                        .trim()
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| {
                            format!("GRO: bad vx on line {i}: {e}")
                        })?;
                let vy: f64 =
                    raw[52..60]
                        .trim()
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| {
                            format!("GRO: bad vy on line {i}: {e}")
                        })?;
                let vz: f64 =
                    raw[60..68]
                        .trim()
                        .parse()
                        .map_err(|e: std::num::ParseFloatError| {
                            format!("GRO: bad vz on line {i}: {e}")
                        })?;
                Some([vx, vy, vz])
            } else {
                None
            };

            atoms.push(GroAtom {
                residue_id,
                residue_name,
                atom_name,
                atom_id,
                position: [x, y, z],
                velocity,
            });
        }

        // Last line: box vectors
        let box_line = lines
            .next()
            .ok_or_else(|| "GRO: missing box vector line".to_string())?
            .map_err(|e| e.to_string())?;
        let parts: Vec<&str> = box_line.split_whitespace().collect();
        if parts.len() < 3 {
            return Err(format!("GRO: box line has only {} tokens", parts.len()));
        }
        let bx: f64 = parts[0]
            .parse()
            .map_err(|e: std::num::ParseFloatError| format!("GRO: bad box_x: {e}"))?;
        let by: f64 = parts[1]
            .parse()
            .map_err(|e: std::num::ParseFloatError| format!("GRO: bad box_y: {e}"))?;
        let bz: f64 = parts[2]
            .parse()
            .map_err(|e: std::num::ParseFloatError| format!("GRO: bad box_z: {e}"))?;

        Ok(GroFile {
            title,
            atoms,
            box_vectors: [bx, by, bz],
        })
    }

    /// Serialize this `GroFile` to any `Write` sink.
    pub fn write(&self, mut writer: impl Write) -> Result<(), String> {
        writeln!(writer, "{}", self.title).map_err(|e| e.to_string())?;
        writeln!(writer, "{:5}", self.atoms.len()).map_err(|e| e.to_string())?;

        for atom in &self.atoms {
            let res_name = truncate(&atom.residue_name, 5);
            let at_name = truncate(&atom.atom_name, 5);
            let [x, y, z] = atom.position;

            if let Some([vx, vy, vz]) = atom.velocity {
                writeln!(
                    writer,
                    "{:5}{:<5}{:<5}{:5}{:8.3}{:8.3}{:8.3}{:8.4}{:8.4}{:8.4}",
                    atom.residue_id, res_name, at_name, atom.atom_id, x, y, z, vx, vy, vz,
                )
                .map_err(|e| e.to_string())?;
            } else {
                writeln!(
                    writer,
                    "{:5}{:<5}{:<5}{:5}{:8.3}{:8.3}{:8.3}",
                    atom.residue_id, res_name, at_name, atom.atom_id, x, y, z,
                )
                .map_err(|e| e.to_string())?;
            }
        }

        writeln!(
            writer,
            "{:10.5} {:10.5} {:10.5}",
            self.box_vectors[0], self.box_vectors[1], self.box_vectors[2],
        )
        .map_err(|e| e.to_string())?;

        Ok(())
    }

    /// Create a minimal `GroFile` from a slice of positions (nm).
    ///
    /// All atoms get `residue_name = "MOL"`, `atom_name = "X"`, `residue_id = 1`.
    pub fn from_xyz(positions: &[[f64; 3]], box_size: [f64; 3]) -> Self {
        let atoms = positions
            .iter()
            .enumerate()
            .map(|(i, &pos)| GroAtom {
                residue_id: 1,
                residue_name: "MOL".to_string(),
                atom_name: "X".to_string(),
                atom_id: (i + 1) as i32,
                position: pos,
                velocity: None,
            })
            .collect();

        GroFile {
            title: "Generated by oxiphysics-io".to_string(),
            atoms,
            box_vectors: box_size,
        }
    }

    /// Return all atom positions converted from nm to Angstrom (x10).
    pub fn positions_angstrom(&self) -> Vec<[f64; 3]> {
        self.atoms
            .iter()
            .map(|a| {
                [
                    a.position[0] * 10.0,
                    a.position[1] * 10.0,
                    a.position[2] * 10.0,
                ]
            })
            .collect()
    }
}

/// Truncate a string slice to at most `n` chars without panicking.
fn truncate(s: &str, n: usize) -> &str {
    let mut idx = s.len();
    while !s.is_char_boundary(idx) || idx > n {
        idx -= 1;
    }
    &s[..idx]
}

// ─── GroFile from_str / to_string ────────────────────────────────────────────

impl GroFile {
    /// Parse a GRO file from a string.
    pub fn parse(data: &str) -> Result<Self, String> {
        let cursor = std::io::Cursor::new(data.as_bytes());
        Self::read(cursor)
    }

    /// Serialize to a String.
    pub fn to_gro_string(&self) -> Result<String, String> {
        let mut buf: Vec<u8> = Vec::new();
        self.write(&mut buf)?;
        String::from_utf8(buf).map_err(|e| e.to_string())
    }
}

impl std::str::FromStr for GroFile {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ============================================================================
// GroFrame type alias
// ============================================================================

/// A single frame of a GROMACS trajectory (title + atoms + box).
///
/// This is a re-export-compatible type alias for [`GroFile`], providing
/// the name used in the task specification.
pub type GroFrame = GroFile;

// ============================================================================
// GroTrajectory
// ============================================================================

/// A multi-frame GRO trajectory (sequence of GRO snapshots).
#[derive(Debug, Clone, Default)]
pub struct GroTrajectory {
    /// Individual frames.
    pub frames: Vec<GroFile>,
}

impl GroTrajectory {
    /// Create an empty trajectory.
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse multiple GRO frames concatenated in a single string.
    pub fn parse(data: &str) -> Result<Self, String> {
        let mut frames = Vec::new();
        let mut remaining = data;

        while !remaining.trim().is_empty() {
            let lines: Vec<&str> = remaining.lines().collect();
            if lines.len() < 3 {
                break;
            }
            // Line 0: title, Line 1: N atoms
            let n_atoms: usize = match lines[1].trim().parse() {
                Ok(n) => n,
                Err(_) => break,
            };
            // Total lines for one frame: 1 (title) + 1 (count) + n_atoms + 1 (box)
            let frame_lines = 3 + n_atoms;
            if lines.len() < frame_lines {
                break;
            }
            let frame_text: String = lines[..frame_lines].join("\n");
            let frame = GroFile::parse(&frame_text)?;
            frames.push(frame);

            // Advance
            let consumed: usize = lines[..frame_lines]
                .iter()
                .map(|l| l.len() + 1) // +1 for newline
                .sum();
            if consumed >= remaining.len() {
                break;
            }
            remaining = &remaining[consumed..];
        }

        Ok(GroTrajectory { frames })
    }

    /// Number of frames.
    pub fn n_frames(&self) -> usize {
        self.frames.len()
    }

    /// Extract the trajectory of a single atom across all frames.
    pub fn atom_trajectory(&self, atom_index: usize) -> Vec<[f64; 3]> {
        self.frames
            .iter()
            .filter_map(|f| f.atoms.get(atom_index).map(|a| a.position))
            .collect()
    }

    /// Mean box vectors across all frames.
    pub fn mean_box(&self) -> [f64; 3] {
        if self.frames.is_empty() {
            return [0.0; 3];
        }
        let n = self.frames.len() as f64;
        let mut sum = [0.0; 3];
        for f in &self.frames {
            sum[0] += f.box_vectors[0];
            sum[1] += f.box_vectors[1];
            sum[2] += f.box_vectors[2];
        }
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
}

impl std::str::FromStr for GroTrajectory {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::parse(s)
    }
}

// ============================================================================
// GroWriter
// ============================================================================

/// A builder for writing GROMACS GRO files with velocities.
///
/// Wraps a `GroFile` and forces all atoms to have their velocity field set
/// before serialisation, padding missing velocities with zero.
pub struct GroWriter {
    inner: GroFile,
}

impl GroWriter {
    /// Create a new `GroWriter` from an existing `GroFile`.
    pub fn new(gro: GroFile) -> Self {
        Self { inner: gro }
    }

    /// Attach velocities to all atoms.
    ///
    /// `velocities` must have the same length as the number of atoms; returns
    /// an error otherwise.
    pub fn write_velocities(&mut self, velocities: &[[f64; 3]]) -> std::result::Result<(), String> {
        if velocities.len() != self.inner.atoms.len() {
            return Err(format!(
                "velocity count {} != atom count {}",
                velocities.len(),
                self.inner.atoms.len()
            ));
        }
        for (atom, &vel) in self.inner.atoms.iter_mut().zip(velocities.iter()) {
            atom.velocity = Some(vel);
        }
        Ok(())
    }

    /// Serialise the GRO file (with velocities) to a string.
    pub fn to_gro_string(&self) -> std::result::Result<String, String> {
        self.inner.to_gro_string()
    }

    /// Access the inner `GroFile`.
    pub fn gro(&self) -> &GroFile {
        &self.inner
    }

    /// Number of atoms.
    pub fn n_atoms(&self) -> usize {
        self.inner.atoms.len()
    }

    /// Returns `true` if every atom has a velocity set.
    pub fn all_have_velocities(&self) -> bool {
        self.inner.atoms.iter().all(|a| a.velocity.is_some())
    }
}

// ============================================================================
// GRO Utility Functions
// ============================================================================

/// Merge two `GroFile` structures into one, renumbering atoms and residues.
///
/// Atoms from `other` are appended after `base`, with residue IDs offset by
/// the maximum residue ID in `base` + 1.
pub fn merge_gro_files(base: &GroFile, other: &GroFile) -> GroFile {
    let max_res_id = base.atoms.iter().map(|a| a.residue_id).max().unwrap_or(0);
    let max_atom_id = base.atoms.iter().map(|a| a.atom_id).max().unwrap_or(0);

    let mut atoms = base.atoms.clone();
    for (i, a) in other.atoms.iter().enumerate() {
        atoms.push(GroAtom {
            residue_id: a.residue_id + max_res_id,
            residue_name: a.residue_name.clone(),
            atom_name: a.atom_name.clone(),
            atom_id: max_atom_id + i as i32 + 1,
            position: a.position,
            velocity: a.velocity,
        });
    }

    let box_vectors = [
        base.box_vectors[0].max(other.box_vectors[0]),
        base.box_vectors[1].max(other.box_vectors[1]),
        base.box_vectors[2].max(other.box_vectors[2]),
    ];

    GroFile {
        title: format!("{} + {}", base.title, other.title),
        atoms,
        box_vectors,
    }
}

/// Translate all atom positions by `delta`.
pub fn translate_gro(gro: &mut GroFile, delta: [f64; 3]) {
    for atom in &mut gro.atoms {
        atom.position[0] += delta[0];
        atom.position[1] += delta[1];
        atom.position[2] += delta[2];
    }
}

/// Compute the centre of geometry of all atoms.
pub fn gro_centre_of_geometry(gro: &GroFile) -> [f64; 3] {
    if gro.atoms.is_empty() {
        return [0.0; 3];
    }
    let n = gro.atoms.len() as f64;
    let mut sum = [0.0; 3];
    for a in &gro.atoms {
        sum[0] += a.position[0];
        sum[1] += a.position[1];
        sum[2] += a.position[2];
    }
    [sum[0] / n, sum[1] / n, sum[2] / n]
}

/// Translate the system so that the centre of geometry is at the origin.
pub fn centre_gro(gro: &mut GroFile) {
    let cog = gro_centre_of_geometry(gro);
    translate_gro(gro, [-cog[0], -cog[1], -cog[2]]);
}

/// Wrap all atom positions into the periodic box `[0, box_x) x [0, box_y) x [0, box_z)`.
pub fn wrap_into_box(gro: &mut GroFile) {
    for atom in &mut gro.atoms {
        for d in 0..3 {
            let l = gro.box_vectors[d];
            if l > 0.0 {
                atom.position[d] = ((atom.position[d] % l) + l) % l;
            }
        }
    }
}

/// Filter atoms by residue name, returning a new `GroFile`.
pub fn filter_by_residue(gro: &GroFile, residue_name: &str) -> GroFile {
    let atoms: Vec<GroAtom> = gro
        .atoms
        .iter()
        .filter(|a| a.residue_name == residue_name)
        .cloned()
        .collect();
    GroFile {
        title: format!("{} [{}]", gro.title, residue_name),
        atoms,
        box_vectors: gro.box_vectors,
    }
}

/// Compute all pairwise distances between atoms.
pub fn gro_pairwise_distances(gro: &GroFile) -> Vec<f64> {
    let n = gro.atoms.len();
    let mut dists = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let a = &gro.atoms[i];
            let b = &gro.atoms[j];
            let dx = a.position[0] - b.position[0];
            let dy = a.position[1] - b.position[1];
            let dz = a.position[2] - b.position[2];
            dists.push((dx * dx + dy * dy + dz * dz).sqrt());
        }
    }
    dists
}

// ============================================================================
// Tests
// ============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn write_to_string(gro: &GroFile) -> String {
        let mut buf: Vec<u8> = Vec::new();
        gro.write(&mut buf).expect("write failed");
        String::from_utf8(buf).expect("utf8")
    }

    #[test]
    fn test_gro_write_read_roundtrip() {
        let original = GroFile {
            title: "Test system".to_string(),
            atoms: vec![
                GroAtom {
                    residue_id: 1,
                    residue_name: "SOL".to_string(),
                    atom_name: "OW".to_string(),
                    atom_id: 1,
                    position: [0.100, 0.200, 0.300],
                    velocity: Some([0.001, -0.002, 0.003]),
                },
                GroAtom {
                    residue_id: 1,
                    residue_name: "SOL".to_string(),
                    atom_name: "HW1".to_string(),
                    atom_id: 2,
                    position: [0.110, 0.210, 0.310],
                    velocity: None,
                },
            ],
            box_vectors: [3.0, 3.0, 3.0],
        };
        let text = write_to_string(&original);
        let parsed = GroFile::read(Cursor::new(text.as_bytes())).expect("parse");
        assert_eq!(parsed.atoms.len(), 2);
        assert_eq!(parsed.title.trim(), "Test system");
        assert!((parsed.box_vectors[0] - 3.0).abs() < 1e-3);
    }

    #[test]
    fn test_gro_from_xyz() {
        let gro = GroFile::from_xyz(&[[1.0, 2.0, 3.0]], [5.0, 5.0, 5.0]);
        assert_eq!(gro.atoms.len(), 1);
        assert_eq!(gro.atoms[0].residue_name, "MOL");
    }

    #[test]
    fn test_positions_angstrom() {
        let gro = GroFile::from_xyz(&[[1.0, 0.0, 0.0]], [5.0, 5.0, 5.0]);
        let ang = gro.positions_angstrom();
        assert!((ang[0][0] - 10.0).abs() < 1e-9);
    }
}

#[cfg(test)]
mod tests_gro_trajectory {
    use super::*;

    fn make_gro_text(n_atoms: usize) -> String {
        let mut s = format!("Test frame\n{n_atoms}\n");
        for i in 0..n_atoms {
            s.push_str(&format!(
                "{:5}{:<5}{:<5}{:5}{:8.3}{:8.3}{:8.3}\n",
                1,
                "SOL",
                "OW",
                i + 1,
                i as f64 * 0.1,
                0.0,
                0.0
            ));
        }
        s.push_str("3.0 3.0 3.0\n");
        s
    }

    #[test]
    fn test_gro_trajectory_single_frame() {
        let text = make_gro_text(3);
        let traj = GroTrajectory::parse(&text).expect("parse");
        assert_eq!(traj.n_frames(), 1);
        assert_eq!(traj.frames[0].atoms.len(), 3);
    }

    #[test]
    fn test_gro_trajectory_two_frames() {
        let mut text = make_gro_text(2);
        text.push_str(&make_gro_text(2));
        let traj = GroTrajectory::parse(&text).expect("parse");
        assert_eq!(traj.n_frames(), 2);
    }

    #[test]
    fn test_gro_trajectory_atom_path() {
        let mut text = make_gro_text(2);
        text.push_str(&make_gro_text(2));
        let traj = GroTrajectory::parse(&text).expect("parse");
        let path = traj.atom_trajectory(0);
        assert_eq!(path.len(), 2);
    }

    #[test]
    fn test_gro_trajectory_mean_box() {
        let text = make_gro_text(1);
        let traj = GroTrajectory::parse(&text).expect("parse");
        let mb = traj.mean_box();
        assert!((mb[0] - 3.0).abs() < 1e-5);
    }

    #[test]
    fn test_gro_trajectory_empty() {
        let traj = GroTrajectory::parse("").expect("parse");
        assert_eq!(traj.n_frames(), 0);
        assert_eq!(traj.mean_box(), [0.0; 3]);
    }
}

#[cfg(test)]
mod tests_gro_utilities {
    use super::*;

    fn simple_gro(n: usize) -> GroFile {
        let atoms: Vec<GroAtom> = (0..n)
            .map(|i| GroAtom {
                residue_id: i as i32 + 1,
                residue_name: "SOL".to_string(),
                atom_name: "OW".to_string(),
                atom_id: i as i32 + 1,
                position: [i as f64 * 0.3, 0.0, 0.0],
                velocity: None,
            })
            .collect();
        GroFile {
            title: "test".to_string(),
            atoms,
            box_vectors: [3.0, 3.0, 3.0],
        }
    }

    #[test]
    fn merge_gro_atom_count() {
        let a = simple_gro(3);
        let b = simple_gro(2);
        let merged = merge_gro_files(&a, &b);
        assert_eq!(merged.atoms.len(), 5);
    }

    #[test]
    fn merge_gro_residue_ids_no_collision() {
        let a = simple_gro(2);
        let b = simple_gro(2);
        let merged = merge_gro_files(&a, &b);
        let ids: Vec<i32> = merged.atoms.iter().map(|x| x.residue_id).collect();
        let unique: std::collections::HashSet<_> = ids.iter().collect();
        assert_eq!(unique.len(), ids.len());
    }

    #[test]
    fn merge_gro_title_combines() {
        let mut a = simple_gro(1);
        a.title = "water".to_string();
        let mut b = simple_gro(1);
        b.title = "ions".to_string();
        let m = merge_gro_files(&a, &b);
        assert!(m.title.contains("water"));
        assert!(m.title.contains("ions"));
    }

    #[test]
    fn translate_gro_shifts_all_atoms() {
        let mut g = simple_gro(3);
        let delta = [1.0, 2.0, 3.0];
        translate_gro(&mut g, delta);
        for (i, a) in g.atoms.iter().enumerate() {
            assert!((a.position[0] - (i as f64 * 0.3 + 1.0)).abs() < 1e-12);
            assert!((a.position[1] - 2.0).abs() < 1e-12);
            assert!((a.position[2] - 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn cog_symmetric_atoms() {
        let mut g = simple_gro(0);
        g.atoms.push(GroAtom {
            residue_id: 1,
            residue_name: "X".to_string(),
            atom_name: "A".to_string(),
            atom_id: 1,
            position: [-1.0, 0.0, 0.0],
            velocity: None,
        });
        g.atoms.push(GroAtom {
            residue_id: 2,
            residue_name: "X".to_string(),
            atom_name: "A".to_string(),
            atom_id: 2,
            position: [1.0, 0.0, 0.0],
            velocity: None,
        });
        let cog = gro_centre_of_geometry(&g);
        assert!(cog[0].abs() < 1e-12);
    }

    #[test]
    fn cog_empty_returns_zero() {
        let g = GroFile {
            title: "".to_string(),
            atoms: vec![],
            box_vectors: [0.0; 3],
        };
        assert_eq!(gro_centre_of_geometry(&g), [0.0; 3]);
    }

    #[test]
    fn centre_gro_places_cog_at_origin() {
        let mut g = simple_gro(4);
        centre_gro(&mut g);
        let cog = gro_centre_of_geometry(&g);
        for c in cog {
            assert!(c.abs() < 1e-12, "cog component {c} not zero");
        }
    }

    #[test]
    fn wrap_into_box_positions_in_range() {
        let mut g = simple_gro(5);
        g.atoms[4].position = [4.0, 4.0, 4.0]; // box is 3x3x3
        wrap_into_box(&mut g);
        for a in &g.atoms {
            assert!(a.position[0] >= 0.0 && a.position[0] < 3.0);
            assert!(a.position[1] >= 0.0 && a.position[1] < 3.0);
            assert!(a.position[2] >= 0.0 && a.position[2] < 3.0);
        }
    }

    #[test]
    fn filter_by_residue_keeps_matching() {
        let mut g = simple_gro(3);
        g.atoms[0].residue_name = "NA".to_string();
        let filtered = filter_by_residue(&g, "SOL");
        assert_eq!(filtered.atoms.len(), 2);
    }

    #[test]
    fn filter_by_residue_none_match_empty() {
        let g = simple_gro(3);
        let filtered = filter_by_residue(&g, "CL");
        assert!(filtered.atoms.is_empty());
    }

    #[test]
    fn pairwise_distances_count() {
        let g = simple_gro(4);
        let d = gro_pairwise_distances(&g);
        assert_eq!(d.len(), 6);
    }

    #[test]
    fn pairwise_distances_positive() {
        let g = simple_gro(3);
        for d in gro_pairwise_distances(&g) {
            assert!(d >= 0.0);
        }
    }

    #[test]
    fn test_gro_writer_write_velocities_ok() {
        let gro = GroFile::from_xyz(&[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], [5.0, 5.0, 5.0]);
        let mut writer = GroWriter::new(gro);
        let vels = [[0.1, 0.0, 0.0], [0.2, 0.0, 0.0]];
        assert!(writer.write_velocities(&vels).is_ok());
        assert!(writer.all_have_velocities());
    }

    #[test]
    fn test_gro_writer_write_velocities_mismatch_error() {
        let gro = GroFile::from_xyz(&[[0.0; 3]], [1.0, 1.0, 1.0]);
        let mut writer = GroWriter::new(gro);
        let result = writer.write_velocities(&[[0.0; 3], [0.0; 3]]);
        assert!(result.is_err());
    }

    #[test]
    fn test_gro_writer_to_string_contains_velocity() {
        let gro = GroFile::from_xyz(&[[0.5, 0.5, 0.5]], [2.0, 2.0, 2.0]);
        let mut writer = GroWriter::new(gro);
        writer
            .write_velocities(&[[1.0, 2.0, 3.0]])
            .expect("write velocities");
        let text = writer.to_gro_string().expect("to_gro_string");
        assert!(!text.is_empty());
        assert_eq!(writer.n_atoms(), 1);
    }

    #[test]
    fn test_gro_frame_type_alias() {
        let frame: GroFrame = GroFile::from_xyz(&[[0.0, 0.0, 0.0]], [3.0, 3.0, 3.0]);
        assert_eq!(frame.atoms.len(), 1);
        let s = frame.to_gro_string().expect("to_gro_string");
        let parsed: GroFrame = GroFile::parse(&s).expect("parse");
        assert_eq!(parsed.atoms.len(), 1);
    }
}

#[cfg(test)]
mod tests_gro_extra {
    use super::*;
    use std::io::BufReader;

    fn minimal_gro(title: &str, n: usize) -> String {
        let mut s = format!("{}\n{}\n", title, n);
        for i in 0..n {
            s.push_str(&format!(
                "{:>5}{:<5}{:<5}{:>5}{:>8.3}{:>8.3}{:>8.3}\n",
                1,
                "SOL",
                "OW",
                i + 1,
                0.0_f64,
                0.0_f64,
                0.0_f64
            ));
        }
        s.push_str("  5.0   5.0   5.0\n");
        s
    }

    #[test]
    fn gro_file_from_str_basic() {
        let data = minimal_gro("Water box", 2);
        let gro = GroFile::parse(&data).expect("parse");
        assert_eq!(gro.atoms.len(), 2);
        assert_eq!(gro.title.trim(), "Water box");
    }

    #[test]
    fn gro_file_box_vectors_parsed() {
        let data = minimal_gro("test", 1);
        let gro = GroFile::parse(&data).expect("parse");
        assert!((gro.box_vectors[0] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn gro_file_roundtrip_via_string() {
        let data = minimal_gro("Round trip", 3);
        let gro = GroFile::parse(&data).expect("parse");
        let s = gro.to_gro_string().expect("to_gro_string");
        let gro2 = GroFile::parse(&s).expect("reparse");
        assert_eq!(gro2.atoms.len(), 3);
    }

    #[test]
    fn gro_file_from_xyz_creates_atoms() {
        let positions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let gro = GroFile::from_xyz(&positions, [10.0, 10.0, 10.0]);
        assert_eq!(gro.atoms.len(), 3);
        assert!((gro.atoms[0].position[0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn gro_file_positions_angstrom_conversion() {
        let gro = GroFile::from_xyz(&[[1.0, 0.0, 0.0]], [5.0, 5.0, 5.0]);
        let ang = gro.positions_angstrom();
        assert!((ang[0][0] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn gro_file_read_from_bufreader() {
        let data = minimal_gro("BufRead test", 2);
        let reader = BufReader::new(data.as_bytes());
        let gro = GroFile::read(reader).expect("parse");
        assert_eq!(gro.atoms.len(), 2);
    }

    #[test]
    fn gro_file_atoms_have_correct_residue_name() {
        let data = minimal_gro("check", 2);
        let gro = GroFile::parse(&data).expect("parse");
        assert_eq!(gro.atoms[0].residue_name, "SOL");
    }

    #[test]
    fn gro_file_atom_id_sequential() {
        let data = minimal_gro("seq", 4);
        let gro = GroFile::parse(&data).expect("parse");
        for (i, atom) in gro.atoms.iter().enumerate() {
            assert_eq!(atom.atom_id, (i + 1) as i32);
        }
    }
}
