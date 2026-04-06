#![allow(clippy::should_implement_trait)]
// Copyright 2026 COOLJAPAN OU (Team KitaSan)
// SPDX-License-Identifier: Apache-2.0

//! LAMMPS dump format reader and writer (frame-oriented API).
//!
//! Supports both positions-only and positions+velocities frames,
//! multi-frame reading, atom property extraction, triclinic box bounds,
//! dump file writing with custom properties, and timestep extraction.
//!
//! Format:
//! ```text
//! ITEM: TIMESTEP
//! 0
//! ITEM: NUMBER OF ATOMS
//! N
//! ITEM: BOX BOUNDS pp pp pp
//! xlo xhi
//! ylo yhi
//! zlo zhi
//! ITEM: ATOMS id type x y z [vx vy vz]
//! 1 1 x y z ...
//! ```

use std::io::{BufRead, Write};

/// A single simulation frame from a LAMMPS dump file.
#[derive(Debug, Clone)]
pub struct LammpsDumpFrame {
    /// Simulation timestep.
    pub timestep: u64,
    /// Number of atoms in this frame.
    pub n_atoms: usize,
    /// Box bounds: `[[xlo, xhi\], [ylo, yhi], [zlo, zhi]]`.
    pub box_bounds: [[f64; 2]; 3],
    /// Atom IDs.
    pub atom_ids: Vec<u32>,
    /// Atom type IDs.
    pub types: Vec<u32>,
    /// Atom positions.
    pub positions: Vec<[f64; 3]>,
    /// Atom velocities (present when the dump includes `vx vy vz` columns).
    pub velocities: Option<Vec<[f64; 3]>>,
    /// Triclinic tilt factors \[xy, xz, yz\], if present.
    pub tilt_factors: Option<[f64; 3]>,
    /// Custom per-atom properties (name -> values).
    pub custom_properties: Vec<(String, Vec<f64>)>,
}

impl LammpsDumpFrame {
    /// Extract positions for atoms of a given type.
    #[allow(dead_code)]
    pub fn positions_by_type(&self, atom_type: u32) -> Vec<[f64; 3]> {
        self.types
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| {
                if t == atom_type {
                    Some(self.positions[i])
                } else {
                    None
                }
            })
            .collect()
    }

    /// Count atoms of a given type.
    #[allow(dead_code)]
    pub fn count_by_type(&self, atom_type: u32) -> usize {
        self.types.iter().filter(|&&t| t == atom_type).count()
    }

    /// Get the box dimensions \[lx, ly, lz\].
    #[allow(dead_code)]
    pub fn box_dimensions(&self) -> [f64; 3] {
        [
            self.box_bounds[0][1] - self.box_bounds[0][0],
            self.box_bounds[1][1] - self.box_bounds[1][0],
            self.box_bounds[2][1] - self.box_bounds[2][0],
        ]
    }

    /// Get the box center \[cx, cy, cz\].
    #[allow(dead_code)]
    pub fn box_center(&self) -> [f64; 3] {
        [
            (self.box_bounds[0][0] + self.box_bounds[0][1]) * 0.5,
            (self.box_bounds[1][0] + self.box_bounds[1][1]) * 0.5,
            (self.box_bounds[2][0] + self.box_bounds[2][1]) * 0.5,
        ]
    }

    /// Get box volume (orthogonal case).
    #[allow(dead_code)]
    pub fn box_volume(&self) -> f64 {
        let dims = self.box_dimensions();
        dims[0] * dims[1] * dims[2]
    }

    /// Compute the center of mass (equal-mass atoms).
    #[allow(dead_code)]
    pub fn center_of_mass(&self) -> [f64; 3] {
        if self.positions.is_empty() {
            return [0.0; 3];
        }
        let n = self.positions.len() as f64;
        let mut com = [0.0; 3];
        for p in &self.positions {
            com[0] += p[0];
            com[1] += p[1];
            com[2] += p[2];
        }
        [com[0] / n, com[1] / n, com[2] / n]
    }

    /// Get a custom property by name.
    #[allow(dead_code)]
    pub fn get_custom_property(&self, name: &str) -> Option<&[f64]> {
        self.custom_properties
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_slice())
    }

    /// Get all distinct atom types in this frame.
    #[allow(dead_code)]
    pub fn atom_types(&self) -> Vec<u32> {
        let mut types: Vec<u32> = self.types.clone();
        types.sort_unstable();
        types.dedup();
        types
    }
}

/// Reader for LAMMPS dump files; reads one frame at a time.
pub struct LammpsDumpReader;

impl LammpsDumpReader {
    /// Attempt to read one frame from the given `BufRead` source.
    ///
    /// Returns `None` when there is no more data, or
    /// `Some(Err(...))` on a parse / format error.
    pub fn read_frame(reader: impl BufRead) -> Option<Result<LammpsDumpFrame, String>> {
        let mut lines = reader.lines().peekable();

        // Drain any leading blank lines / EOF
        loop {
            match lines.peek() {
                None => return None,
                Some(Ok(l)) if l.trim().is_empty() => {
                    lines.next();
                }
                Some(Err(_)) => {
                    let e = lines
                        .next()
                        .expect("iterator should have elements")
                        .expect_err("matched Err(_) arm so value must be Err");
                    return Some(Err(e.to_string()));
                }
                _ => break,
            }
        }

        // ---- helpers ---------------------------------------------------
        macro_rules! next_line {
            ($label:expr) => {{
                match lines.next() {
                    None => return Some(Err(format!("LAMMPS dump: unexpected EOF at {}", $label))),
                    Some(Err(e)) => return Some(Err(e.to_string())),
                    Some(Ok(l)) => l,
                }
            }};
        }

        macro_rules! expect_keyword {
            ($kw:expr) => {{
                let l = next_line!($kw);
                if !l.trim().starts_with($kw) {
                    return Some(Err(format!(
                        "LAMMPS dump: expected '{}', got '{}'",
                        $kw,
                        l.trim()
                    )));
                }
                l
            }};
        }
        // ----------------------------------------------------------------

        expect_keyword!("ITEM: TIMESTEP");
        let ts_str = next_line!("timestep value");
        let timestep: u64 = match ts_str.trim().parse() {
            Ok(v) => v,
            Err(e) => return Some(Err(format!("LAMMPS dump: bad timestep: {e}"))),
        };

        expect_keyword!("ITEM: NUMBER OF ATOMS");
        let na_str = next_line!("n_atoms value");
        let n_atoms: usize = match na_str.trim().parse() {
            Ok(v) => v,
            Err(e) => return Some(Err(format!("LAMMPS dump: bad n_atoms: {e}"))),
        };

        let box_header = expect_keyword!("ITEM: BOX BOUNDS");
        let is_triclinic = box_header.contains("xy xz yz");

        let mut box_bounds = [[0.0f64; 2]; 3];
        let mut tilt_factors: Option<[f64; 3]> = None;

        if is_triclinic {
            let mut tilts = [0.0f64; 3];
            for (i, row) in box_bounds.iter_mut().enumerate() {
                let bl = next_line!(format!("box bound {i}"));
                let parts: Vec<&str> = bl.split_whitespace().collect();
                if parts.len() < 3 {
                    return Some(Err(format!(
                        "LAMMPS dump: triclinic box line {i} too short"
                    )));
                }
                let lo: f64 = match parts[0].parse() {
                    Ok(v) => v,
                    Err(e) => return Some(Err(format!("LAMMPS dump: box lo[{i}]: {e}"))),
                };
                let hi: f64 = match parts[1].parse() {
                    Ok(v) => v,
                    Err(e) => return Some(Err(format!("LAMMPS dump: box hi[{i}]: {e}"))),
                };
                let tilt: f64 = match parts[2].parse() {
                    Ok(v) => v,
                    Err(e) => return Some(Err(format!("LAMMPS dump: tilt[{i}]: {e}"))),
                };
                *row = [lo, hi];
                tilts[i] = tilt;
            }
            tilt_factors = Some(tilts);
        } else {
            for (i, row) in box_bounds.iter_mut().enumerate() {
                let bl = next_line!(format!("box bound {i}"));
                let parts: Vec<&str> = bl.split_whitespace().collect();
                if parts.len() < 2 {
                    return Some(Err(format!("LAMMPS dump: box line {i} too short")));
                }
                let lo: f64 = match parts[0].parse() {
                    Ok(v) => v,
                    Err(e) => return Some(Err(format!("LAMMPS dump: box lo[{i}]: {e}"))),
                };
                let hi: f64 = match parts[1].parse() {
                    Ok(v) => v,
                    Err(e) => return Some(Err(format!("LAMMPS dump: box hi[{i}]: {e}"))),
                };
                *row = [lo, hi];
            }
        }

        // ITEM: ATOMS id type x y z [vx vy vz] [custom...]
        let header_line = expect_keyword!("ITEM: ATOMS");
        let has_velocities = header_line.contains("vx");

        // Parse column names from the header
        let col_names: Vec<&str> = header_line
            .trim()
            .strip_prefix("ITEM: ATOMS")
            .unwrap_or("")
            .split_whitespace()
            .collect();

        // Identify custom property columns (beyond id, type, x, y, z, vx, vy, vz)
        let standard_cols = ["id", "type", "x", "y", "z", "vx", "vy", "vz"];
        let custom_col_indices: Vec<(usize, String)> = col_names
            .iter()
            .enumerate()
            .filter(|(_, name)| !standard_cols.contains(name))
            .map(|(i, name)| (i, name.to_string()))
            .collect();

        let mut atom_ids = Vec::with_capacity(n_atoms);
        let mut types = Vec::with_capacity(n_atoms);
        let mut positions = Vec::with_capacity(n_atoms);
        let mut velocities: Vec<[f64; 3]> = if has_velocities {
            Vec::with_capacity(n_atoms)
        } else {
            Vec::new()
        };
        let mut custom_data: Vec<Vec<f64>> = custom_col_indices
            .iter()
            .map(|_| Vec::with_capacity(n_atoms))
            .collect();

        for i in 0..n_atoms {
            let al = next_line!(format!("atom {i}"));
            let parts: Vec<&str> = al.split_whitespace().collect();
            let min_cols = if has_velocities { 8 } else { 5 };
            if parts.len() < min_cols {
                return Some(Err(format!(
                    "LAMMPS dump: atom line {i} has {} cols, need {}",
                    parts.len(),
                    min_cols
                )));
            }
            macro_rules! parse_col {
                ($idx:expr, $ty:ty, $label:expr) => {
                    match parts[$idx].parse::<$ty>() {
                        Ok(v) => v,
                        Err(e) => {
                            return Some(Err(format!(
                                "LAMMPS dump: atom {i} col {} ({}): {e}",
                                $idx, $label
                            )))
                        }
                    }
                };
            }
            let id = parse_col!(0, u32, "id");
            let ty = parse_col!(1, u32, "type");
            let x = parse_col!(2, f64, "x");
            let y = parse_col!(3, f64, "y");
            let z = parse_col!(4, f64, "z");
            atom_ids.push(id);
            types.push(ty);
            positions.push([x, y, z]);
            if has_velocities {
                let vx = parse_col!(5, f64, "vx");
                let vy = parse_col!(6, f64, "vy");
                let vz = parse_col!(7, f64, "vz");
                velocities.push([vx, vy, vz]);
            }

            // Parse custom properties
            for (ci, (col_idx, _)) in custom_col_indices.iter().enumerate() {
                if let Some(val_str) = parts.get(*col_idx) {
                    let val = val_str.parse::<f64>().unwrap_or(0.0);
                    custom_data[ci].push(val);
                }
            }
        }

        let custom_properties: Vec<(String, Vec<f64>)> = custom_col_indices
            .into_iter()
            .zip(custom_data)
            .map(|((_, name), data)| (name, data))
            .collect();

        Some(Ok(LammpsDumpFrame {
            timestep,
            n_atoms,
            box_bounds,
            atom_ids,
            types,
            positions,
            velocities: if has_velocities {
                Some(velocities)
            } else {
                None
            },
            tilt_factors,
            custom_properties,
        }))
    }

    /// Read all frames from a `BufRead` source.
    #[allow(dead_code)]
    pub fn read_all_frames(content: &str) -> Result<Vec<LammpsDumpFrame>, String> {
        let mut frames = Vec::new();
        // Split by "ITEM: TIMESTEP" to get individual frame chunks
        let chunks: Vec<&str> = content
            .split("ITEM: TIMESTEP\n")
            .filter(|s| !s.trim().is_empty())
            .collect();

        for chunk in chunks {
            let full = format!("ITEM: TIMESTEP\n{chunk}");
            let cursor = std::io::Cursor::new(full.as_bytes());
            match Self::read_frame(cursor) {
                Some(Ok(frame)) => frames.push(frame),
                Some(Err(e)) => return Err(e),
                None => {}
            }
        }
        Ok(frames)
    }

    /// Extract all timesteps from a dump content without fully parsing frames.
    #[allow(dead_code)]
    pub fn extract_timesteps(content: &str) -> Vec<u64> {
        let mut timesteps = Vec::new();
        let mut lines = content.lines();
        while let Some(line) = lines.next() {
            if line.trim() == "ITEM: TIMESTEP"
                && let Some(ts_line) = lines.next()
                && let Ok(ts) = ts_line.trim().parse::<u64>()
            {
                timesteps.push(ts);
            }
        }
        timesteps
    }
}

/// Writer for LAMMPS dump files.
pub struct LammpsDumpWriter;

impl LammpsDumpWriter {
    /// Write one frame to any `Write` sink.
    pub fn write_frame(mut writer: impl Write, frame: &LammpsDumpFrame) -> Result<(), String> {
        writeln!(writer, "ITEM: TIMESTEP").map_err(|e| e.to_string())?;
        writeln!(writer, "{}", frame.timestep).map_err(|e| e.to_string())?;
        writeln!(writer, "ITEM: NUMBER OF ATOMS").map_err(|e| e.to_string())?;
        writeln!(writer, "{}", frame.n_atoms).map_err(|e| e.to_string())?;

        if let Some(tilts) = frame.tilt_factors {
            writeln!(writer, "ITEM: BOX BOUNDS xy xz yz pp pp pp").map_err(|e| e.to_string())?;
            for (i, b) in frame.box_bounds.iter().enumerate() {
                writeln!(writer, "{} {} {}", b[0], b[1], tilts[i]).map_err(|e| e.to_string())?;
            }
        } else {
            writeln!(writer, "ITEM: BOX BOUNDS pp pp pp").map_err(|e| e.to_string())?;
            for b in &frame.box_bounds {
                writeln!(writer, "{} {}", b[0], b[1]).map_err(|e| e.to_string())?;
            }
        }

        let has_vel = frame.velocities.is_some();
        let custom_names: Vec<&str> = frame
            .custom_properties
            .iter()
            .map(|(n, _)| n.as_str())
            .collect();

        // Build header
        let mut header = String::from("ITEM: ATOMS id type x y z");
        if has_vel {
            header.push_str(" vx vy vz");
        }
        for name in &custom_names {
            header.push(' ');
            header.push_str(name);
        }
        writeln!(writer, "{}", header).map_err(|e| e.to_string())?;

        let vels = frame.velocities.as_deref();
        for i in 0..frame.n_atoms {
            let [x, y, z] = frame.positions[i];
            let mut line = if has_vel {
                let [vx, vy, vz] = vels.expect("value should be present")[i];
                format!(
                    "{} {} {} {} {} {} {} {}",
                    frame.atom_ids[i], frame.types[i], x, y, z, vx, vy, vz,
                )
            } else {
                format!("{} {} {} {} {}", frame.atom_ids[i], frame.types[i], x, y, z,)
            };

            for (_, values) in &frame.custom_properties {
                if i < values.len() {
                    line.push_str(&format!(" {}", values[i]));
                }
            }

            writeln!(writer, "{}", line).map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    /// Write multiple frames to a single `Write` sink.
    #[allow(dead_code)]
    pub fn write_frames(writer: impl Write, frames: &[LammpsDumpFrame]) -> Result<(), String> {
        let mut writer = writer;
        for frame in frames {
            Self::write_frame(&mut writer, frame)?;
        }
        Ok(())
    }
}

/// Compute the number density (atoms / volume) for a frame.
#[allow(dead_code)]
pub fn number_density(frame: &LammpsDumpFrame) -> f64 {
    let vol = frame.box_volume();
    if vol.abs() < 1e-30 {
        return 0.0;
    }
    frame.n_atoms as f64 / vol
}

/// Compute mean kinetic energy per atom (assuming unit mass).
/// KE = 0.5 * v^2 per atom.
#[allow(dead_code)]
pub fn mean_kinetic_energy(frame: &LammpsDumpFrame) -> Option<f64> {
    let vels = frame.velocities.as_ref()?;
    if vels.is_empty() {
        return Some(0.0);
    }
    let total: f64 = vels
        .iter()
        .map(|v| 0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]))
        .sum();
    Some(total / vels.len() as f64)
}

// ─── Per-atom stress extraction ───────────────────────────────────────────────

/// Extract per-atom stress tensor components from custom properties.
///
/// LAMMPS may output the 6 independent Voigt components:
/// `c_stress[1..6]` or `sxx syy szz sxy sxz syz`.
#[allow(dead_code)]
pub fn extract_stress_tensor(frame: &LammpsDumpFrame) -> Option<Vec<[f64; 6]>> {
    // Try Voigt-named properties
    let names = ["sxx", "syy", "szz", "sxy", "sxz", "syz"];
    let props: Vec<Option<&[f64]>> = names.iter().map(|n| frame.get_custom_property(n)).collect();

    if props.iter().all(|p| p.is_some()) {
        let n = frame.n_atoms;
        let result: Vec<[f64; 6]> = (0..n)
            .map(|i| {
                [
                    props[0].expect("value should be present")[i],
                    props[1].expect("value should be present")[i],
                    props[2].expect("value should be present")[i],
                    props[3].expect("value should be present")[i],
                    props[4].expect("value should be present")[i],
                    props[5].expect("value should be present")[i],
                ]
            })
            .collect();
        return Some(result);
    }

    None
}

/// Compute the von Mises stress from a Voigt stress tensor \[s11,s22,s33,s12,s13,s23\].
#[allow(dead_code)]
pub fn von_mises_stress(s: [f64; 6]) -> f64 {
    let [s11, s22, s33, s12, s13, s23] = s;
    (0.5 * ((s11 - s22).powi(2) + (s22 - s33).powi(2) + (s33 - s11).powi(2))
        + 3.0 * (s12 * s12 + s13 * s13 + s23 * s23))
        .sqrt()
}

/// Hydrostatic pressure from a Voigt stress tensor: p = -(s11 + s22 + s33) / 3.
#[allow(dead_code)]
pub fn hydrostatic_pressure(s: [f64; 6]) -> f64 {
    -(s[0] + s[1] + s[2]) / 3.0
}

// ─── Dump frame time series ────────────────────────────────────────────────────

/// Analysis utilities for a time series of LAMMPS dump frames.
pub struct DumpTimeSeries {
    /// All frames in the series.
    pub frames: Vec<LammpsDumpFrame>,
}

impl DumpTimeSeries {
    /// Construct from a vector of frames.
    #[allow(dead_code)]
    pub fn new(frames: Vec<LammpsDumpFrame>) -> Self {
        Self { frames }
    }

    /// Parse from a multi-frame dump string.
    #[allow(dead_code)]
    pub fn from_str(content: &str) -> Result<Self, String> {
        let frames = LammpsDumpReader::read_all_frames(content)?;
        Ok(Self { frames })
    }

    /// Number of frames.
    #[allow(dead_code)]
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether there are no frames.
    #[allow(dead_code)]
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Extract all timesteps in order.
    #[allow(dead_code)]
    pub fn timesteps(&self) -> Vec<u64> {
        self.frames.iter().map(|f| f.timestep).collect()
    }

    /// Extract mean kinetic energy per frame.
    #[allow(dead_code)]
    pub fn mean_ke_series(&self) -> Vec<Option<f64>> {
        self.frames.iter().map(mean_kinetic_energy).collect()
    }

    /// Extract number density per frame.
    #[allow(dead_code)]
    pub fn number_density_series(&self) -> Vec<f64> {
        self.frames.iter().map(number_density).collect()
    }

    /// Get the frame at a given timestep, if present.
    #[allow(dead_code)]
    pub fn frame_at_timestep(&self, ts: u64) -> Option<&LammpsDumpFrame> {
        self.frames.iter().find(|f| f.timestep == ts)
    }

    /// Mean position of atom type `atom_type` averaged over all frames.
    #[allow(dead_code)]
    pub fn mean_position_by_type(&self, atom_type: u32) -> Option<[f64; 3]> {
        let mut sum = [0.0f64; 3];
        let mut count = 0usize;
        for frame in &self.frames {
            for pos in frame.positions_by_type(atom_type) {
                sum[0] += pos[0];
                sum[1] += pos[1];
                sum[2] += pos[2];
                count += 1;
            }
        }
        if count == 0 {
            return None;
        }
        Some([
            sum[0] / count as f64,
            sum[1] / count as f64,
            sum[2] / count as f64,
        ])
    }
}

// ─── LAMMPS atom styles ────────────────────────────────────────────────────────

/// Supported LAMMPS atom styles for dump parsing.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)]
pub enum LammpsAtomStyle {
    /// Positions only (atomic style).
    Atomic,
    /// Positions + molecule ID + charge.
    Full,
    /// Positions + molecule ID.
    Molecular,
    /// Custom (user-defined columns).
    Custom,
}

/// Parse the atom style from the ITEM: ATOMS header line.
#[allow(dead_code)]
pub fn detect_atom_style(header: &str) -> LammpsAtomStyle {
    let lower = header.to_lowercase();
    let has_q = lower.contains(" q ");
    let has_mol = lower.contains(" mol ") || lower.contains("mol\n");
    match (has_mol, has_q) {
        (true, true) => LammpsAtomStyle::Full,
        (true, false) => LammpsAtomStyle::Molecular,
        _ => LammpsAtomStyle::Atomic,
    }
}

// ─── Radial distribution function (g(r)) ─────────────────────────────────────

/// Compute the radial distribution function g(r) from a single frame.
///
/// Returns `(r_bins, g_of_r)` where `r_bins[i]` is the centre of bin `i`
/// and `g_of_r[i]` is g(r) at that bin.
///
/// Only considers atoms of `type_a` (centre) and `type_b` (neighbours).
/// Use `type_a = type_b = 0` to include all atoms.
#[allow(dead_code)]
pub fn radial_distribution_function(
    frame: &LammpsDumpFrame,
    type_a: u32,
    type_b: u32,
    r_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<f64>) {
    let dr = r_max / n_bins as f64;
    let mut hist = vec![0u64; n_bins];

    let n = frame.n_atoms;
    let mut n_a = 0usize;
    let mut n_b = 0usize;

    for i in 0..n {
        let ti = frame.types[i];
        if type_a != 0 && ti != type_a {
            continue;
        }
        n_a += 1;
        for j in 0..n {
            if i == j {
                continue;
            }
            let tj = frame.types[j];
            if type_b != 0 && tj != type_b {
                continue;
            }
            n_b += 1;
            let pi = frame.positions[i];
            let pj = frame.positions[j];
            let dx = pi[0] - pj[0];
            let dy = pi[1] - pj[1];
            let dz = pi[2] - pj[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            if r < r_max {
                let bin = (r / dr) as usize;
                if bin < n_bins {
                    hist[bin] += 1;
                }
            }
        }
    }

    let vol = frame.box_volume();
    let rho = if n_b > 0 && vol > 0.0 {
        n_b as f64 / vol
    } else {
        1.0
    };

    let r_bins: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * dr).collect();
    let g_of_r: Vec<f64> = r_bins
        .iter()
        .zip(hist.iter())
        .map(|(&r, &count)| {
            let shell_vol = 4.0 * std::f64::consts::PI * r * r * dr;
            let n_ideal = rho * shell_vol;
            if n_ideal < 1e-30 || n_a == 0 {
                return 0.0;
            }
            count as f64 / (n_a as f64 * n_ideal)
        })
        .collect();

    (r_bins, g_of_r)
}

// ─── Per-atom force extraction ────────────────────────────────────────────────

/// Extract per-atom force vectors from custom properties `fx`, `fy`, `fz`.
#[allow(dead_code)]
pub fn extract_forces(frame: &LammpsDumpFrame) -> Option<Vec<[f64; 3]>> {
    let fx = frame.get_custom_property("fx")?;
    let fy = frame.get_custom_property("fy")?;
    let fz = frame.get_custom_property("fz")?;
    let n = frame.n_atoms;
    let forces: Vec<[f64; 3]> = (0..n).map(|i| [fx[i], fy[i], fz[i]]).collect();
    Some(forces)
}

/// Mean force magnitude per atom.
#[allow(dead_code)]
pub fn mean_force_magnitude(frame: &LammpsDumpFrame) -> Option<f64> {
    let forces = extract_forces(frame)?;
    if forces.is_empty() {
        return Some(0.0);
    }
    let total: f64 = forces
        .iter()
        .map(|f| (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt())
        .sum();
    Some(total / forces.len() as f64)
}

// ─── Triclinic box utilities ──────────────────────────────────────────────────

/// Convert triclinic box bounds + tilt factors to orthogonal lattice vectors.
///
/// Returns `(a, b, c)` where each is a `[f64; 3]` Cartesian vector.
///
/// LAMMPS triclinic convention:
/// a = \[lx,  0,  0\]
/// b = \[xy, ly,  0\]
/// c = \[xz, yz, lz\]
#[allow(dead_code)]
pub fn triclinic_lattice_vectors(
    box_bounds: [[f64; 2]; 3],
    tilts: [f64; 3],
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let lx = box_bounds[0][1] - box_bounds[0][0];
    let ly = box_bounds[1][1] - box_bounds[1][0];
    let lz = box_bounds[2][1] - box_bounds[2][0];
    let [xy, xz, yz] = tilts;
    let a = [lx, 0.0, 0.0];
    let b = [xy, ly, 0.0];
    let c = [xz, yz, lz];
    (a, b, c)
}

/// Triclinic box volume from lattice vectors.
#[allow(dead_code)]
pub fn triclinic_volume(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    // V = a · (b × c)
    let bc = [
        b[1] * c[2] - b[2] * c[1],
        b[2] * c[0] - b[0] * c[2],
        b[0] * c[1] - b[1] * c[0],
    ];
    a[0] * bc[0] + a[1] * bc[1] + a[2] * bc[2]
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn make_frame(timestep: u64, with_velocities: bool) -> LammpsDumpFrame {
        let positions = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 1.0],
        ];
        let velocities = if with_velocities {
            Some(vec![
                [0.1, 0.2, 0.3],
                [-0.1, 0.0, 0.1],
                [0.0, -0.2, 0.0],
                [0.5, -0.5, 0.5],
            ])
        } else {
            None
        };
        LammpsDumpFrame {
            timestep,
            n_atoms: 4,
            box_bounds: [[0.0, 5.0], [0.0, 5.0], [0.0, 5.0]],
            atom_ids: vec![1, 2, 3, 4],
            types: vec![1, 1, 2, 2],
            positions,
            velocities,
            tilt_factors: None,
            custom_properties: Vec::new(),
        }
    }

    #[test]
    fn test_lammps_write_read_roundtrip() {
        let frame = make_frame(100, false);
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, &frame).expect("write failed");

        let parsed = LammpsDumpReader::read_frame(Cursor::new(&buf))
            .expect("no frame found")
            .expect("parse error");

        assert_eq!(parsed.timestep, 100);
        assert_eq!(parsed.n_atoms, 4);
        assert_eq!(parsed.atom_ids, vec![1, 2, 3, 4]);
        assert_eq!(parsed.types, vec![1, 1, 2, 2]);
        assert!((parsed.positions[0][0]).abs() < 1e-10);
        assert!((parsed.positions[1][0] - 1.0).abs() < 1e-10);
        assert!((parsed.positions[3][2] - 1.0).abs() < 1e-10);
        assert!(parsed.velocities.is_none());
    }

    #[test]
    fn test_lammps_multi_frame() {
        let mut buf: Vec<u8> = Vec::new();
        for ts in [0u64, 100, 200] {
            let frame = make_frame(ts, false);
            LammpsDumpWriter::write_frame(&mut buf, &frame).expect("write failed");
        }

        // Read all three frames from the combined buffer
        let text = String::from_utf8(buf).expect("utf8");
        let mut timesteps_found = Vec::new();

        // Split by "ITEM: TIMESTEP" to get individual frame chunks
        let frames_text: Vec<&str> = text
            .split("ITEM: TIMESTEP\n")
            .filter(|s| !s.trim().is_empty())
            .collect();
        assert_eq!(frames_text.len(), 3, "expected 3 frames");

        for chunk in &frames_text {
            let full = format!("ITEM: TIMESTEP\n{chunk}");
            let result = LammpsDumpReader::read_frame(Cursor::new(full.as_bytes()))
                .expect("no frame")
                .expect("parse error");
            timesteps_found.push(result.timestep);
        }

        assert_eq!(timesteps_found, vec![0, 100, 200]);
    }

    #[test]
    fn test_lammps_box_bounds() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 1,
            box_bounds: [[-2.5, 2.5], [0.0, 10.0], [1.0, 9.0]],
            atom_ids: vec![1],
            types: vec![1],
            positions: vec![[0.0, 5.0, 5.0]],
            velocities: None,
            tilt_factors: None,
            custom_properties: Vec::new(),
        };
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, &frame).expect("write");

        let parsed = LammpsDumpReader::read_frame(Cursor::new(&buf))
            .expect("no frame")
            .expect("parse");

        assert!((parsed.box_bounds[0][0] - (-2.5)).abs() < 1e-10);
        assert!((parsed.box_bounds[0][1] - 2.5).abs() < 1e-10);
        assert!((parsed.box_bounds[1][0]).abs() < 1e-10);
        assert!((parsed.box_bounds[1][1] - 10.0).abs() < 1e-10);
        assert!((parsed.box_bounds[2][0] - 1.0).abs() < 1e-10);
        assert!((parsed.box_bounds[2][1] - 9.0).abs() < 1e-10);
    }

    #[test]
    fn test_lammps_with_velocities() {
        let frame = make_frame(50, true);
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, &frame).expect("write");

        let parsed = LammpsDumpReader::read_frame(Cursor::new(&buf))
            .expect("no frame")
            .expect("parse");

        assert_eq!(parsed.timestep, 50);
        let vels = parsed.velocities.expect("velocities missing");
        assert_eq!(vels.len(), 4);
        assert!((vels[0][0] - 0.1).abs() < 1e-10);
        assert!((vels[0][1] - 0.2).abs() < 1e-10);
        assert!((vels[0][2] - 0.3).abs() < 1e-10);
        assert!((vels[1][0] - (-0.1)).abs() < 1e-10);
        assert!((vels[3][0] - 0.5).abs() < 1e-10);
        assert!((vels[3][2] - 0.5).abs() < 1e-10);
    }

    // ─── New tests ────────────────────────────────────────────────────

    #[test]
    fn test_positions_by_type() {
        let frame = make_frame(0, false);
        let type1 = frame.positions_by_type(1);
        assert_eq!(type1.len(), 2);
        let type2 = frame.positions_by_type(2);
        assert_eq!(type2.len(), 2);
        let type3 = frame.positions_by_type(3);
        assert!(type3.is_empty());
    }

    #[test]
    fn test_count_by_type() {
        let frame = make_frame(0, false);
        assert_eq!(frame.count_by_type(1), 2);
        assert_eq!(frame.count_by_type(2), 2);
        assert_eq!(frame.count_by_type(99), 0);
    }

    #[test]
    fn test_box_dimensions() {
        let frame = make_frame(0, false);
        let dims = frame.box_dimensions();
        assert!((dims[0] - 5.0).abs() < 1e-12);
        assert!((dims[1] - 5.0).abs() < 1e-12);
        assert!((dims[2] - 5.0).abs() < 1e-12);
    }

    #[test]
    fn test_box_center() {
        let frame = make_frame(0, false);
        let center = frame.box_center();
        assert!((center[0] - 2.5).abs() < 1e-12);
        assert!((center[1] - 2.5).abs() < 1e-12);
        assert!((center[2] - 2.5).abs() < 1e-12);
    }

    #[test]
    fn test_box_volume() {
        let frame = make_frame(0, false);
        assert!((frame.box_volume() - 125.0).abs() < 1e-10);
    }

    #[test]
    fn test_center_of_mass() {
        let frame = make_frame(0, false);
        let com = frame.center_of_mass();
        // Positions: [0,0,0], [1,0,0], [0,1,0], [1,1,1]
        // Mean: [0.5, 0.5, 0.25]
        assert!((com[0] - 0.5).abs() < 1e-12);
        assert!((com[1] - 0.5).abs() < 1e-12);
        assert!((com[2] - 0.25).abs() < 1e-12);
    }

    #[test]
    fn test_center_of_mass_empty() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 0,
            box_bounds: [[0.0, 1.0]; 3],
            atom_ids: vec![],
            types: vec![],
            positions: vec![],
            velocities: None,
            tilt_factors: None,
            custom_properties: Vec::new(),
        };
        let com = frame.center_of_mass();
        assert_eq!(com, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_atom_types() {
        let frame = make_frame(0, false);
        let types = frame.atom_types();
        assert_eq!(types, vec![1, 2]);
    }

    #[test]
    fn test_read_all_frames() {
        let mut buf: Vec<u8> = Vec::new();
        for ts in [0u64, 100, 200] {
            LammpsDumpWriter::write_frame(&mut buf, &make_frame(ts, false)).unwrap();
        }
        let text = String::from_utf8(buf).unwrap();
        let frames = LammpsDumpReader::read_all_frames(&text).unwrap();
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[0].timestep, 0);
        assert_eq!(frames[1].timestep, 100);
        assert_eq!(frames[2].timestep, 200);
    }

    #[test]
    fn test_extract_timesteps() {
        let mut buf: Vec<u8> = Vec::new();
        for ts in [10u64, 20, 30] {
            LammpsDumpWriter::write_frame(&mut buf, &make_frame(ts, false)).unwrap();
        }
        let text = String::from_utf8(buf).unwrap();
        let timesteps = LammpsDumpReader::extract_timesteps(&text);
        assert_eq!(timesteps, vec![10, 20, 30]);
    }

    #[test]
    fn test_triclinic_write_read() {
        let frame = LammpsDumpFrame {
            timestep: 42,
            n_atoms: 2,
            box_bounds: [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0]],
            atom_ids: vec![1, 2],
            types: vec![1, 1],
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: None,
            tilt_factors: Some([0.5, 0.0, -0.3]),
            custom_properties: Vec::new(),
        };
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, &frame).unwrap();

        let parsed = LammpsDumpReader::read_frame(Cursor::new(&buf))
            .unwrap()
            .unwrap();

        assert_eq!(parsed.timestep, 42);
        let tilts = parsed.tilt_factors.expect("tilt factors missing");
        assert!((tilts[0] - 0.5).abs() < 1e-10);
        assert!((tilts[1] - 0.0).abs() < 1e-10);
        assert!((tilts[2] - (-0.3)).abs() < 1e-10);
    }

    #[test]
    fn test_number_density() {
        let frame = make_frame(0, false);
        let rho = number_density(&frame);
        // 4 atoms / 125.0 volume = 0.032
        assert!((rho - 0.032).abs() < 1e-10);
    }

    #[test]
    fn test_mean_kinetic_energy() {
        let frame = make_frame(0, true);
        let ke = mean_kinetic_energy(&frame).unwrap();
        // vels: [0.1,0.2,0.3], [-0.1,0.0,0.1], [0.0,-0.2,0.0], [0.5,-0.5,0.5]
        // KE per atom: 0.5*(0.01+0.04+0.09)=0.07, 0.5*(0.01+0+0.01)=0.01,
        //              0.5*(0+0.04+0)=0.02, 0.5*(0.25+0.25+0.25)=0.375
        // total = 0.475, mean = 0.475/4 = 0.11875
        assert!((ke - 0.11875).abs() < 1e-10);
    }

    #[test]
    fn test_mean_kinetic_energy_no_velocities() {
        let frame = make_frame(0, false);
        assert!(mean_kinetic_energy(&frame).is_none());
    }

    #[test]
    fn test_write_frames_multiple() {
        let frames: Vec<LammpsDumpFrame> = (0..3).map(|i| make_frame(i * 50, false)).collect();
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frames(&mut buf, &frames).unwrap();
        let text = String::from_utf8(buf).unwrap();
        let timesteps = LammpsDumpReader::extract_timesteps(&text);
        assert_eq!(timesteps, vec![0, 50, 100]);
    }

    #[test]
    fn test_custom_properties_roundtrip() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 2,
            box_bounds: [[0.0, 5.0]; 3],
            atom_ids: vec![1, 2],
            types: vec![1, 1],
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: None,
            tilt_factors: None,
            custom_properties: vec![("energy".to_string(), vec![10.5, 20.3])],
        };
        let mut buf: Vec<u8> = Vec::new();
        LammpsDumpWriter::write_frame(&mut buf, &frame).unwrap();
        let text = String::from_utf8(buf).unwrap();
        // Verify the header includes the custom property
        assert!(text.contains("energy"));
    }

    #[test]
    fn test_get_custom_property() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 2,
            box_bounds: [[0.0, 1.0]; 3],
            atom_ids: vec![1, 2],
            types: vec![1, 1],
            positions: vec![[0.0; 3]; 2],
            velocities: None,
            tilt_factors: None,
            custom_properties: vec![("charge".to_string(), vec![1.0, -1.0])],
        };
        let charge = frame.get_custom_property("charge").unwrap();
        assert_eq!(charge.len(), 2);
        assert!((charge[0] - 1.0).abs() < 1e-12);
        assert!((charge[1] - (-1.0)).abs() < 1e-12);
        assert!(frame.get_custom_property("missing").is_none());
    }

    #[test]
    fn test_box_dimensions_asymmetric() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 0,
            box_bounds: [[-5.0, 5.0], [0.0, 20.0], [10.0, 30.0]],
            atom_ids: vec![],
            types: vec![],
            positions: vec![],
            velocities: None,
            tilt_factors: None,
            custom_properties: Vec::new(),
        };
        let dims = frame.box_dimensions();
        assert!((dims[0] - 10.0).abs() < 1e-12);
        assert!((dims[1] - 20.0).abs() < 1e-12);
        assert!((dims[2] - 20.0).abs() < 1e-12);
    }

    // ─── Tests for new functions ───────────────────────────────────────

    fn make_frame_with_stress(timestep: u64) -> LammpsDumpFrame {
        LammpsDumpFrame {
            timestep,
            n_atoms: 2,
            box_bounds: [[0.0, 10.0]; 3],
            atom_ids: vec![1, 2],
            types: vec![1, 1],
            positions: vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            velocities: None,
            tilt_factors: None,
            custom_properties: vec![
                ("sxx".to_string(), vec![1.0, 2.0]),
                ("syy".to_string(), vec![3.0, 4.0]),
                ("szz".to_string(), vec![5.0, 6.0]),
                ("sxy".to_string(), vec![0.5, 0.5]),
                ("sxz".to_string(), vec![0.2, 0.2]),
                ("syz".to_string(), vec![0.1, 0.1]),
            ],
        }
    }

    #[test]
    fn test_extract_stress_tensor_present() {
        let frame = make_frame_with_stress(0);
        let tensors = extract_stress_tensor(&frame).expect("stress tensors should be present");
        assert_eq!(tensors.len(), 2);
        assert!((tensors[0][0] - 1.0).abs() < 1e-12); // sxx for atom 0
        assert!((tensors[1][0] - 2.0).abs() < 1e-12); // sxx for atom 1
        assert!((tensors[0][2] - 5.0).abs() < 1e-12); // szz for atom 0
    }

    #[test]
    fn test_extract_stress_tensor_missing() {
        let frame = make_frame(0, false);
        assert!(extract_stress_tensor(&frame).is_none());
    }

    #[test]
    fn test_von_mises_stress_isotropic() {
        // Pure hydrostatic: s11=s22=s33=p, shear=0 → von Mises = 0
        let s = [5.0, 5.0, 5.0, 0.0, 0.0, 0.0];
        let vm = von_mises_stress(s);
        assert!(vm.abs() < 1e-10);
    }

    #[test]
    fn test_von_mises_stress_uniaxial() {
        // s11=1, rest=0 → VM = sqrt(0.5*((1-0)^2+(0-0)^2+(0-1)^2)) = sqrt(1) = 1
        let s = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let vm = von_mises_stress(s);
        assert!((vm - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hydrostatic_pressure() {
        let s = [3.0, 3.0, 3.0, 0.0, 0.0, 0.0];
        let p = hydrostatic_pressure(s);
        assert!((p - (-3.0)).abs() < 1e-12);
    }

    #[test]
    fn test_hydrostatic_pressure_compression() {
        // Negative diagonal = compression → positive pressure
        let s = [-10.0, -10.0, -10.0, 0.0, 0.0, 0.0];
        let p = hydrostatic_pressure(s);
        assert!((p - 10.0).abs() < 1e-12);
    }

    #[test]
    fn test_dump_time_series_new_len() {
        let frames: Vec<LammpsDumpFrame> = (0..5).map(|i| make_frame(i * 100, false)).collect();
        let ts = DumpTimeSeries::new(frames);
        assert_eq!(ts.len(), 5);
        assert!(!ts.is_empty());
    }

    #[test]
    fn test_dump_time_series_empty() {
        let ts = DumpTimeSeries::new(vec![]);
        assert_eq!(ts.len(), 0);
        assert!(ts.is_empty());
    }

    #[test]
    fn test_dump_time_series_timesteps() {
        let frames: Vec<LammpsDumpFrame> = (0..3).map(|i| make_frame(i * 50, false)).collect();
        let ts = DumpTimeSeries::new(frames);
        assert_eq!(ts.timesteps(), vec![0, 50, 100]);
    }

    #[test]
    fn test_dump_time_series_from_str() {
        let mut buf: Vec<u8> = Vec::new();
        for i in 0..3u64 {
            LammpsDumpWriter::write_frame(&mut buf, &make_frame(i * 10, false)).unwrap();
        }
        let content = String::from_utf8(buf).unwrap();
        let ts = DumpTimeSeries::from_str(&content).unwrap();
        assert_eq!(ts.len(), 3);
        assert_eq!(ts.timesteps(), vec![0, 10, 20]);
    }

    #[test]
    fn test_dump_time_series_mean_ke_series() {
        let frames: Vec<LammpsDumpFrame> = vec![
            make_frame(0, true),
            make_frame(1, false),
            make_frame(2, true),
        ];
        let ts = DumpTimeSeries::new(frames);
        let ke_series = ts.mean_ke_series();
        assert_eq!(ke_series.len(), 3);
        assert!(ke_series[0].is_some());
        assert!(ke_series[1].is_none());
        assert!(ke_series[2].is_some());
    }

    #[test]
    fn test_dump_time_series_number_density_series() {
        let frames: Vec<LammpsDumpFrame> = (0..3).map(|i| make_frame(i, false)).collect();
        let ts = DumpTimeSeries::new(frames);
        let rho_series = ts.number_density_series();
        // Each frame has 4 atoms in 5^3 = 125 volume
        assert_eq!(rho_series.len(), 3);
        for rho in &rho_series {
            assert!((rho - 0.032).abs() < 1e-10);
        }
    }

    #[test]
    fn test_dump_time_series_frame_at_timestep() {
        let frames: Vec<LammpsDumpFrame> = (0..5).map(|i| make_frame(i * 100, false)).collect();
        let ts = DumpTimeSeries::new(frames);
        assert!(ts.frame_at_timestep(200).is_some());
        assert_eq!(ts.frame_at_timestep(200).unwrap().timestep, 200);
        assert!(ts.frame_at_timestep(999).is_none());
    }

    #[test]
    fn test_dump_time_series_mean_position_by_type() {
        let frames: Vec<LammpsDumpFrame> = (0..2).map(|i| make_frame(i, false)).collect();
        let ts = DumpTimeSeries::new(frames);
        // type 1 positions across all frames: [0,0,0] and [1,0,0] repeated twice
        // mean = [0.5, 0, 0]
        let mean = ts
            .mean_position_by_type(1)
            .expect("type 1 should have positions");
        assert!((mean[0] - 0.5).abs() < 1e-12);
        assert!((mean[1] - 0.0).abs() < 1e-12);
        assert!((mean[2] - 0.0).abs() < 1e-12);
        assert!(ts.mean_position_by_type(99).is_none());
    }

    #[test]
    fn test_detect_atom_style_atomic() {
        let style = detect_atom_style("ITEM: ATOMS id type x y z");
        assert_eq!(style, LammpsAtomStyle::Atomic);
    }

    #[test]
    fn test_detect_atom_style_full() {
        let style = detect_atom_style("ITEM: ATOMS id mol type q x y z");
        assert_eq!(style, LammpsAtomStyle::Full);
    }

    #[test]
    fn test_detect_atom_style_molecular() {
        let style = detect_atom_style("ITEM: ATOMS id mol type x y z");
        assert_eq!(style, LammpsAtomStyle::Molecular);
    }

    #[test]
    fn test_rdf_length() {
        let frame = make_frame(0, false);
        let n_bins = 20;
        let (r_bins, g) = radial_distribution_function(&frame, 0, 0, 2.5, n_bins);
        assert_eq!(r_bins.len(), n_bins);
        assert_eq!(g.len(), n_bins);
    }

    #[test]
    fn test_rdf_r_bins_positive() {
        let frame = make_frame(0, false);
        let (r_bins, _) = radial_distribution_function(&frame, 0, 0, 2.0, 10);
        for &r in &r_bins {
            assert!(r > 0.0);
        }
    }

    #[test]
    fn test_rdf_nonnegative() {
        let frame = make_frame(0, false);
        let (_, g) = radial_distribution_function(&frame, 0, 0, 2.5, 20);
        for &v in &g {
            assert!(v >= 0.0);
        }
    }

    #[test]
    fn test_extract_forces_present() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 2,
            box_bounds: [[0.0, 5.0]; 3],
            atom_ids: vec![1, 2],
            types: vec![1, 1],
            positions: vec![[0.0; 3]; 2],
            velocities: None,
            tilt_factors: None,
            custom_properties: vec![
                ("fx".to_string(), vec![1.0, -1.0]),
                ("fy".to_string(), vec![0.5, -0.5]),
                ("fz".to_string(), vec![0.0, 0.0]),
            ],
        };
        let forces = extract_forces(&frame).expect("forces should be present");
        assert_eq!(forces.len(), 2);
        assert!((forces[0][0] - 1.0).abs() < 1e-12);
        assert!((forces[1][0] - (-1.0)).abs() < 1e-12);
    }

    #[test]
    fn test_extract_forces_missing() {
        let frame = make_frame(0, false);
        assert!(extract_forces(&frame).is_none());
    }

    #[test]
    fn test_mean_force_magnitude() {
        let frame = LammpsDumpFrame {
            timestep: 0,
            n_atoms: 1,
            box_bounds: [[0.0, 5.0]; 3],
            atom_ids: vec![1],
            types: vec![1],
            positions: vec![[0.0; 3]],
            velocities: None,
            tilt_factors: None,
            custom_properties: vec![
                ("fx".to_string(), vec![3.0]),
                ("fy".to_string(), vec![4.0]),
                ("fz".to_string(), vec![0.0]),
            ],
        };
        let mean_f = mean_force_magnitude(&frame).expect("should have forces");
        assert!((mean_f - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_triclinic_lattice_vectors_orthogonal() {
        let bounds = [[0.0, 10.0], [0.0, 8.0], [0.0, 6.0]];
        let tilts = [0.0, 0.0, 0.0];
        let (a, b, c) = triclinic_lattice_vectors(bounds, tilts);
        assert!((a[0] - 10.0).abs() < 1e-12);
        assert!((a[1]).abs() < 1e-12);
        assert!((b[0]).abs() < 1e-12);
        assert!((b[1] - 8.0).abs() < 1e-12);
        assert!((c[2] - 6.0).abs() < 1e-12);
    }

    #[test]
    fn test_triclinic_lattice_vectors_tilted() {
        let bounds = [[0.0, 10.0], [0.0, 10.0], [0.0, 10.0]];
        let tilts = [2.0, 1.0, 0.5];
        let (a, b, c) = triclinic_lattice_vectors(bounds, tilts);
        assert!((a[0] - 10.0).abs() < 1e-12);
        assert!((b[0] - 2.0).abs() < 1e-12); // xy
        assert!((c[0] - 1.0).abs() < 1e-12); // xz
        assert!((c[1] - 0.5).abs() < 1e-12); // yz
    }

    #[test]
    fn test_triclinic_volume_orthogonal() {
        let a = [10.0, 0.0, 0.0];
        let b = [0.0, 8.0, 0.0];
        let c = [0.0, 0.0, 6.0];
        let vol = triclinic_volume(a, b, c);
        assert!((vol - 480.0).abs() < 1e-10);
    }

    #[test]
    fn test_triclinic_volume_tilted() {
        // Cube with tilt: volume should still be lx*ly*lz for small tilts
        let bounds = [[0.0, 10.0]; 3];
        let tilts = [1.0, 0.0, 0.0];
        let (a, b, c) = triclinic_lattice_vectors(bounds, tilts);
        let vol = triclinic_volume(a, b, c);
        // V = |a·(b×c)| = lx*ly*lz = 1000 regardless of xy tilt
        assert!((vol - 1000.0).abs() < 1e-8);
    }
}
