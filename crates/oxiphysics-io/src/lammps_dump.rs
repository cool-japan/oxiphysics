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
//!
//! The columns of the `ITEM: ATOMS` table are located by name, so any column
//! order is read, including the scaled coordinates (`xs ys zs`) of `dump atom`.

use std::io::{BufRead, Write};

#[cfg(test)]
mod frame_stream_tests;
mod layout;
#[cfg(test)]
mod layout_tests;
mod rdf;
#[cfg(test)]
mod rdf_tests;
#[cfg(test)]
mod writer_tests;

use layout::{AtomTable, DumpBox, DumpLayout};
pub use rdf::radial_distribution_function;

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
    ///
    /// A type without a matching entry in [`LammpsDumpFrame::positions`] (the
    /// arrays disagree) is skipped.
    pub fn positions_by_type(&self, atom_type: u32) -> Vec<[f64; 3]> {
        self.types
            .iter()
            .zip(&self.positions)
            .filter(|&(&t, _)| t == atom_type)
            .map(|(_, &position)| position)
            .collect()
    }

    /// Count atoms of a given type.
    pub fn count_by_type(&self, atom_type: u32) -> usize {
        self.types.iter().filter(|&&t| t == atom_type).count()
    }

    /// Get the box dimensions \[lx, ly, lz\].
    pub fn box_dimensions(&self) -> [f64; 3] {
        [
            self.box_bounds[0][1] - self.box_bounds[0][0],
            self.box_bounds[1][1] - self.box_bounds[1][0],
            self.box_bounds[2][1] - self.box_bounds[2][0],
        ]
    }

    /// Get the box center \[cx, cy, cz\].
    pub fn box_center(&self) -> [f64; 3] {
        [
            (self.box_bounds[0][0] + self.box_bounds[0][1]) * 0.5,
            (self.box_bounds[1][0] + self.box_bounds[1][1]) * 0.5,
            (self.box_bounds[2][0] + self.box_bounds[2][1]) * 0.5,
        ]
    }

    /// Get box volume (orthogonal case).
    pub fn box_volume(&self) -> f64 {
        let dims = self.box_dimensions();
        dims[0] * dims[1] * dims[2]
    }

    /// Compute the center of mass (equal-mass atoms).
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
    pub fn get_custom_property(&self, name: &str) -> Option<&[f64]> {
        self.custom_properties
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, v)| v.as_slice())
    }

    /// Get all distinct atom types in this frame.
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
    ///
    /// The columns of the `ITEM: ATOMS` table are located by the names in its
    /// header, in any order.  `id`, `type` and a complete position triple
    /// (`x y z`, `xu yu zu`, `xs ys zs` or `xsu ysu zsu`) are required; scaled
    /// coordinates are converted to Cartesian ones with the (triclinic-aware)
    /// box.  `vx vy vz` become [`LammpsDumpFrame::velocities`] when all three are
    /// present.  Every other numeric column, such as `mol`, `q` or `fx`, is
    /// returned in [`LammpsDumpFrame::custom_properties`]; an `element` column
    /// holds text and is skipped.  A row with fewer values than the header names
    /// columns, or a value that is not a number, is an error.
    ///
    /// Blank lines and the `ITEM: UNITS` / `ITEM: TIME` blocks that
    /// `dump_modify units yes time yes` writes in front of a frame are skipped.
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
                    return match lines.next() {
                        Some(Err(e)) => Some(Err(e.to_string())),
                        _ => None,
                    };
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

        // An `ITEM: UNITS` or `ITEM: TIME` block is the item line and one value line.
        loop {
            let item = next_line!("ITEM: TIMESTEP");
            match item.trim() {
                "ITEM: UNITS" | "ITEM: TIME" => {
                    next_line!(item.trim());
                }
                other if other.starts_with("ITEM: TIMESTEP") => break,
                other => {
                    return Some(Err(format!(
                        "LAMMPS dump: expected 'ITEM: TIMESTEP', got '{other}'"
                    )));
                }
            }
        }
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

        // ITEM: ATOMS <column names>: the columns are found by name, in any order.
        let header_line = expect_keyword!("ITEM: ATOMS");
        let col_names: Vec<&str> = header_line
            .trim()
            .strip_prefix("ITEM: ATOMS")
            .unwrap_or("")
            .split_whitespace()
            .collect();
        let layout = match DumpLayout::from_header(&col_names) {
            Ok(layout) => layout,
            Err(e) => return Some(Err(e)),
        };
        let dump_box = DumpBox::from_bounds(box_bounds, tilt_factors);
        let mut table = AtomTable::new(&layout, n_atoms);
        for i in 0..n_atoms {
            let al = next_line!(format!("atom {i}"));
            if let Err(e) = table.push_row(&layout, &dump_box, &al, i) {
                return Some(Err(e));
            }
        }
        let custom_properties: Vec<(String, Vec<f64>)> = layout
            .custom_names()
            .map(str::to_string)
            .zip(table.custom)
            .collect();

        Some(Ok(LammpsDumpFrame {
            timestep,
            n_atoms,
            box_bounds,
            atom_ids: table.ids,
            types: table.types,
            positions: table.positions,
            velocities: table.velocities,
            tilt_factors,
            custom_properties,
        }))
    }

    /// Read every frame of the dump text `content`, in file order.
    ///
    /// The frames are read one after another with
    /// [`LammpsDumpReader::read_frame`], so line endings may be `\n` or `\r\n`,
    /// and blank lines and `ITEM: UNITS` / `ITEM: TIME` blocks between frames
    /// are skipped.  Anything else between two frames, or a frame that is cut
    /// short, is an error.
    pub fn read_all_frames(content: &str) -> Result<Vec<LammpsDumpFrame>, String> {
        let mut stream = std::io::Cursor::new(content.as_bytes());
        let mut frames = Vec::new();
        while let Some(frame) = Self::read_frame(&mut stream) {
            frames.push(frame?);
        }
        Ok(frames)
    }

    /// Extract all timesteps from a dump content without fully parsing frames.
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

/// The first `n` entries of a per-atom array, or an error naming the array when
/// it is shorter.
fn first_entries<'a, T>(values: &'a [T], n: usize, what: &str) -> Result<&'a [T], String> {
    values.get(..n).ok_or_else(|| {
        format!(
            "LAMMPS dump: frame declares {n} atoms but {what} holds only {} entries",
            values.len()
        )
    })
}

impl LammpsDumpWriter {
    /// Write one frame to any `Write` sink.
    ///
    /// `atom_ids`, `types`, `positions`, `velocities` (when present) and every
    /// custom property must hold at least `n_atoms` values.  A shorter array is
    /// reported as an error before anything is written; entries beyond
    /// `n_atoms` are not written.
    pub fn write_frame(mut writer: impl Write, frame: &LammpsDumpFrame) -> Result<(), String> {
        let n_atoms = frame.n_atoms;
        let atom_ids = first_entries(&frame.atom_ids, n_atoms, "atom_ids")?;
        let types = first_entries(&frame.types, n_atoms, "types")?;
        let positions = first_entries(&frame.positions, n_atoms, "positions")?;
        let velocities = match &frame.velocities {
            Some(values) => Some(first_entries(values, n_atoms, "velocities")?),
            None => None,
        };
        let mut custom_values: Vec<&[f64]> = Vec::with_capacity(frame.custom_properties.len());
        for (name, values) in &frame.custom_properties {
            custom_values.push(first_entries(
                values,
                n_atoms,
                &format!("custom property '{name}'"),
            )?);
        }
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

        let custom_names: Vec<&str> = frame
            .custom_properties
            .iter()
            .map(|(n, _)| n.as_str())
            .collect();

        // Build header
        let mut header = String::from("ITEM: ATOMS id type x y z");
        if velocities.is_some() {
            header.push_str(" vx vy vz");
        }
        for name in &custom_names {
            header.push(' ');
            header.push_str(name);
        }
        writeln!(writer, "{}", header).map_err(|e| e.to_string())?;

        // Every slice holds exactly `n_atoms` entries, so `i` is in range for all of them.
        let rows = atom_ids.iter().zip(types).zip(positions).enumerate();
        for (i, ((id, atom_type), [x, y, z])) in rows {
            let mut line = format!("{id} {atom_type} {x} {y} {z}");
            if let Some(velocities) = velocities {
                let [vx, vy, vz] = velocities[i];
                line.push_str(&format!(" {vx} {vy} {vz}"));
            }
            for values in &custom_values {
                line.push_str(&format!(" {}", values[i]));
            }
            writeln!(writer, "{line}").map_err(|e| e.to_string())?;
        }
        Ok(())
    }

    /// Write multiple frames to a single `Write` sink.
    pub fn write_frames(writer: impl Write, frames: &[LammpsDumpFrame]) -> Result<(), String> {
        let mut writer = writer;
        for frame in frames {
            Self::write_frame(&mut writer, frame)?;
        }
        Ok(())
    }
}

/// Compute the number density (atoms / volume) for a frame.
pub fn number_density(frame: &LammpsDumpFrame) -> f64 {
    let vol = frame.box_volume();
    if vol.abs() < 1e-30 {
        return 0.0;
    }
    frame.n_atoms as f64 / vol
}

/// Compute mean kinetic energy per atom (assuming unit mass).
/// KE = 0.5 * v^2 per atom.
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
///
/// Returns `None` unless all six named components (`sxx syy szz sxy sxz syz`)
/// are present and hold a value for each of the `n_atoms` atoms.
pub fn extract_stress_tensor(frame: &LammpsDumpFrame) -> Option<Vec<[f64; 6]>> {
    let sxx = frame.get_custom_property("sxx")?;
    let syy = frame.get_custom_property("syy")?;
    let szz = frame.get_custom_property("szz")?;
    let sxy = frame.get_custom_property("sxy")?;
    let sxz = frame.get_custom_property("sxz")?;
    let syz = frame.get_custom_property("syz")?;
    (0..frame.n_atoms)
        .map(|i| {
            Some([
                *sxx.get(i)?,
                *syy.get(i)?,
                *szz.get(i)?,
                *sxy.get(i)?,
                *sxz.get(i)?,
                *syz.get(i)?,
            ])
        })
        .collect()
}

/// Compute the von Mises stress from a Voigt stress tensor \[s11,s22,s33,s12,s13,s23\].
pub fn von_mises_stress(s: [f64; 6]) -> f64 {
    let [s11, s22, s33, s12, s13, s23] = s;
    (0.5 * ((s11 - s22).powi(2) + (s22 - s33).powi(2) + (s33 - s11).powi(2))
        + 3.0 * (s12 * s12 + s13 * s13 + s23 * s23))
        .sqrt()
}

/// Hydrostatic pressure from a Voigt stress tensor: p = -(s11 + s22 + s33) / 3.
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
    pub fn new(frames: Vec<LammpsDumpFrame>) -> Self {
        Self { frames }
    }

    /// Number of frames.
    pub fn len(&self) -> usize {
        self.frames.len()
    }

    /// Whether there are no frames.
    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// Extract all timesteps in order.
    pub fn timesteps(&self) -> Vec<u64> {
        self.frames.iter().map(|f| f.timestep).collect()
    }

    /// Extract mean kinetic energy per frame.
    pub fn mean_ke_series(&self) -> Vec<Option<f64>> {
        self.frames.iter().map(mean_kinetic_energy).collect()
    }

    /// Extract number density per frame.
    pub fn number_density_series(&self) -> Vec<f64> {
        self.frames.iter().map(number_density).collect()
    }

    /// Get the frame at a given timestep, if present.
    pub fn frame_at_timestep(&self, ts: u64) -> Option<&LammpsDumpFrame> {
        self.frames.iter().find(|f| f.timestep == ts)
    }

    /// Mean position of atom type `atom_type` averaged over all frames.
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

impl std::str::FromStr for DumpTimeSeries {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let frames = LammpsDumpReader::read_all_frames(s)?;
        Ok(Self { frames })
    }
}

// ─── LAMMPS atom styles ────────────────────────────────────────────────────────

/// Supported LAMMPS atom styles for dump parsing.
#[derive(Debug, Clone, Copy, PartialEq)]
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

// ─── Per-atom force extraction ────────────────────────────────────────────────

/// Extract per-atom force vectors from custom properties `fx`, `fy`, `fz`.
///
/// Returns `None` unless all three components are present and hold a value for
/// each of the `n_atoms` atoms.
pub fn extract_forces(frame: &LammpsDumpFrame) -> Option<Vec<[f64; 3]>> {
    let fx = frame.get_custom_property("fx")?;
    let fy = frame.get_custom_property("fy")?;
    let fz = frame.get_custom_property("fz")?;
    (0..frame.n_atoms)
        .map(|i| Some([*fx.get(i)?, *fy.get(i)?, *fz.get(i)?]))
        .collect()
}

/// Mean force magnitude per atom.
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
        let ts = content.parse::<DumpTimeSeries>().unwrap();
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
